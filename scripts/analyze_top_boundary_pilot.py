"""Frozen analysis and decision pipeline for the top-boundary causal pilot."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from compare_lateral_domain_sequences import read_sequence
from storm_dynamics.diagnostic_capture import curl


PRIMARY = {
    "zeta_max_s-1": 1,
    "circulation_4200_m2_s": 1,
    "vtheta_max_m_s": 1,
    "zeta_halfmax_width_m": -1,
    "convergence_max_s-1": 1,
    "stretching_conditional_mean_s-2": 1,
}
PROGNOSTIC = ("u", "v", "w", "p", "theta", "qv", "ql", "qi", "qr", "qs", "qg", "qh")


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def geometric_gate(path, damping_faces=4):
    frames = []
    all_finite = True
    with h5py.File(path, "r") as f:
        z = np.asarray(f["grid/zc"]); zf = np.asarray(f["grid/zf"])
        damp_bottom = float(zf[-damping_faces])
        for name, g in f["snapshots"].items():
            t = float(g.attrs["time_s"])
            if not 2790 <= t <= 3301:
                continue
            all_finite &= all(np.isfinite(g[field][:]).all() for field in PROGNOSTIC)
            qcond = sum(np.asarray(g[n]) for n in ("ql", "qi", "qr", "qs", "qg", "qh"))
            w = 0.5 * (np.asarray(g["w"][:, :, :-1]) + np.asarray(g["w"][:, :, 1:]))

            def top_and_margin(mask):
                occupied = mask.any(axis=(0, 1))
                top = float(z[occupied].max()) if occupied.any() else None
                return top, (None if top is None else damp_bottom - top)

            qtop, qmargin = top_and_margin(qcond > 1e-5)
            w5top, w5margin = top_and_margin(w > 5.0)
            w1top, w1margin = top_and_margin(np.abs(w) > 1.0)
            frames.append({
                "snapshot": name, "time_s": t,
                "qcond_1e-5_top_m": qtop, "qcond_1e-5_margin_to_damping_m": qmargin,
                "qcond_1e-5_top_cell_area_km2": float((qcond[:, :, -1] > 1e-5).sum()
                                                        * (f["grid/xc"][1] - f["grid/xc"][0])
                                                        * (f["grid/yc"][1] - f["grid/yc"][0]) / 1e6),
                "w_gt5_top_m": w5top, "w_gt5_margin_to_damping_m": w5margin,
                "abs_w_gt1_top_m": w1top, "abs_w_gt1_margin_to_damping_m": w1margin,
                "w_gt5_in_damping": bool((w[:, :, z >= damp_bottom] > 5.0).any()),
                "abs_w_gt1_in_damping": bool((np.abs(w[:, :, z >= damp_bottom]) > 1.0).any()),
            })
    margins = [r["abs_w_gt1_margin_to_damping_m"] for r in frames
               if r["abs_w_gt1_margin_to_damping_m"] is not None]
    gate = {
        "condensate_not_in_top_cell": all(r["qcond_1e-5_top_cell_area_km2"] == 0 for r in frames),
        "updraft_gt5_outside_damping": not any(r["w_gt5_in_damping"] for r in frames),
        "abs_w_gt1_positive_margin": not any(r["abs_w_gt1_in_damping"] for r in frames),
        "all_prognostic_fields_finite": bool(all_finite),
        "minimum_abs_w_gt1_margin_m": float(min(margins)) if margins else None,
        "damping_bottom_m": damp_bottom,
    }
    gate["status"] = "PASS" if all(v for k, v in gate.items()
                                      if k not in ("minimum_abs_w_gt1_margin_m", "damping_bottom_m", "status")) else "FAIL"
    return frames, gate


def budget_closure(path):
    rows = []
    with h5py.File(path, "r") as f:
        nk = int(f.attrs["budget_nz_with_halo"])
        coords = tuple(np.asarray(f[f"grid/{n}"]) for n in ("xc", "yc", "zc"))
        coords = (*coords[:2], coords[2][:nk])
        snapshots = list(f["snapshots"].values())

        def velocity(group):
            return tuple(np.asarray(group[n][:, :, :nk + (n == "w")]) for n in ("u", "v", "w"))

        for a, b in zip(snapshots, snapshots[1:]):
            total = sum((curl(velocity(stage), coords) for stage in b["increments"].values()))
            change = curl(velocity(b), coords) - curl(velocity(a), coords)
            residual = change - total
            scale = np.sqrt(np.mean(change * change))
            rows.append({
                "time_s": float(b.attrs["time_s"]),
                "absolute_rms_s-1": float(np.sqrt(np.mean(residual * residual))),
                "relative_rms": float(np.sqrt(np.mean(residual * residual)) / max(scale, 1e-30)),
                "max_abs_s-1": float(np.max(np.abs(residual))),
            })
    maximum = max((r["relative_rms"] for r in rows), default=float("inf"))
    return rows, {"maximum_relative_rms": maximum, "threshold": 1e-10,
                  "status": "PASS" if maximum <= 1e-10 else "FAIL"}


def tracking_gate(rows_a, rows_b, domain_m=72000.0):
    def periodic_distance(a, b):
        dx = abs(a["center_x_m"] - b["center_x_m"]); dx = min(dx, domain_m - dx)
        dy = abs(a["center_y_m"] - b["center_y_m"]); dy = min(dy, domain_m - dy)
        return float(np.hypot(dx, dy))

    jumps_a = [periodic_distance(a, b) for a, b in zip(rows_a, rows_a[1:])]
    jumps_b = [periodic_distance(a, b) for a, b in zip(rows_b, rows_b[1:])]
    cross = [periodic_distance(a, b) for a, b in zip(rows_a, rows_b)]
    status = (len(rows_a) == len(rows_b) == 18 and max(jumps_a, default=0) < 6000.0
              and max(jumps_b, default=0) < 6000.0 and max(cross, default=0) < 12000.0)
    return {"frames_each": [len(rows_a), len(rows_b)],
            "max_control_jump_m": max(jumps_a, default=None),
            "max_high_top_jump_m": max(jumps_b, default=None),
            "max_cross_branch_separation_m": max(cross, default=None),
            "status": "PASS" if status else "FAIL"}


def aggregate_top_observer(path):
    with h5py.File(path, "r") as f:
        if f.attrs["status"] != "complete":
            raise RuntimeError(f"incomplete top observer: {path}")
        step = np.asarray(f["step"], dtype=int)
        time = np.asarray(f["time_s"], dtype=float)
        ordinal = np.asarray(f["ordinal"], dtype=int)
        expected = np.tile(np.arange(4), len(step) // 4)
        ordered = len(step) % 4 == 0 and np.array_equal(ordinal, expected)
        energy = np.asarray(f["energy_removed_J"]).sum(axis=1)
        flux_delta = (np.asarray(f["pressure_flux_after_W"])
                      - np.asarray(f["pressure_flux_before_W"])).sum(axis=1)
        curl_driver = (np.asarray(f["delta_xi_rms"]).mean(axis=1)
                       + np.asarray(f["delta_eta_rms"]).mean(axis=1))
        z = 0.5 * (np.asarray(f["zf_m"][:-1]) + np.asarray(f["zf_m"][1:]))
        low = z <= 2000.0
        conv = np.asarray(f["convergence_rms"])[:, low].mean(axis=1)
        wrms = np.asarray(f["w_rms"])[:, low].mean(axis=1)
        rows = []
        for s in np.unique(step):
            take = step == s
            post = take & (ordinal == 3)
            rows.append({
                "step": int(s), "time_s": float(time[take][0]),
                "energy_removed_J": float(energy[take].sum()),
                "pressure_flux_change_W": float(flux_delta[take].sum()),
                "curl_driver_s-1": float(curl_driver[take].sum()),
                "low_level_convergence_rms_s-1": float(conv[post][0]),
                "low_level_w_rms_m_s": float(wrms[post][0]),
            })
        integrity = {
            "calls": len(step), "steps": len(rows), "exactly_four_calls_per_step": len(step) == 4 * len(rows),
            "ordered_contexts": bool(ordered), "all_finite": bool(all(np.isfinite(np.asarray(f[n])).all()
                for n in ("energy_removed_J", "delta_xi_rms", "delta_eta_rms",
                          "pressure_flux_before_W", "pressure_flux_after_W", "w_rms", "convergence_rms"))),
            "full_field_calls": len(f["full_fields/step"]) if "full_fields" in f else 0,
        }
        integrity["status"] = "PASS" if (integrity["exactly_four_calls_per_step"]
            and integrity["ordered_contexts"] and integrity["all_finite"]
            and integrity["full_field_calls"] > 0) else "FAIL"
    return rows, integrity


def lag_correlation(rows, driver, response, max_lag_s=600.0, spacing_s=1.0):
    time = np.asarray([r["time_s"] for r in rows], dtype=float)
    keep = (time >= 2790.0) & (time <= 3300.1)
    time = time[keep]
    x = np.asarray([r[driver] for r in rows], dtype=float)[keep]
    y = np.asarray([r[response] for r in rows], dtype=float)[keep]
    grid = np.arange(np.ceil(time[0]), np.floor(time[-1]) + spacing_s / 2, spacing_s)
    x = np.interp(grid, time, x); y = np.interp(grid, time, y)
    q = np.arange(len(grid), dtype=float)
    x = x - np.polyval(np.polyfit(q, x, 1), q)
    y = y - np.polyval(np.polyfit(q, y, 1), q)
    rows_out = []
    for lag in np.arange(0.0, max_lag_s + spacing_s / 2, spacing_s):
        n = int(round(lag / spacing_s))
        xa, ya = (x, y) if n == 0 else (x[:-n], y[n:])
        corr = float(np.corrcoef(xa, ya)[0, 1]) if len(xa) >= 20 and xa.std() and ya.std() else float("nan")
        rows_out.append({"driver": driver, "response": response, "lag_s": float(lag), "correlation": corr})
    finite = [r for r in rows_out if np.isfinite(r["correlation"])]
    peak = max(finite, key=lambda r: abs(r["correlation"])) if finite else None
    return rows_out, peak


def metric_comparison(control, high):
    ta = np.asarray([r["time_s"] for r in control]); tb = np.asarray([r["time_s"] for r in high])
    if len(ta) != 18 or not np.array_equal(ta, tb):
        raise RuntimeError("the two arms must contain the same 18 native output times")
    metrics = sorted((set(control[0]) & set(high[0])) - {"time_s", "center_x_m", "center_y_m"})
    stats, paired = [], []
    for metric in metrics:
        a = np.asarray([r[metric] for r in control], dtype=float)
        b = np.asarray([r[metric] for r in high], dtype=float)
        valid = np.isfinite(a) & np.isfinite(b)
        delta = b - a; scale = np.nanmedian(np.abs(a[valid]))
        frac = float(np.nanmedian(delta[valid]) / scale) if scale else float("nan")
        row = {"metric": metric, "median_control15": float(np.nanmedian(a)),
               "median_hightop20": float(np.nanmedian(b)), "median_delta": float(np.nanmedian(delta)),
               "median_fractional_delta": frac, "high_greater_pairs": int(np.sum(b[valid] > a[valid])),
               "high_less_pairs": int(np.sum(b[valid] < a[valid])), "valid_pairs": int(valid.sum())}
        stats.append(row)
        paired.extend({"time_s": float(t), "metric": metric, "control15": float(x),
                       "hightop20": float(y), "delta": float(y - x)} for t, x, y in zip(ta, a, b))
    primary = []
    for metric, direction in PRIMARY.items():
        row = next(r for r in stats if r["metric"] == metric)
        favorable = row["high_greater_pairs"] if direction > 0 else row["high_less_pairs"]
        primary.append({**row, "favorable_direction": direction, "favorable_pairs": favorable,
                        "coherent_12_of_18": favorable >= 12,
                        "material_20_percent": direction * row["median_fractional_delta"] >= 0.20})
    material = sum(r["coherent_12_of_18"] and r["material_20_percent"] for r in primary)
    classification = ("material favorable top response" if material >= 3 else
                      "no dominant favorable top effect" if material <= 1 else "mixed top response")
    return stats, paired, primary, material, classification


def figures(out, paired, primary, observer):
    import matplotlib.pyplot as plt

    labels = [r["metric"].replace("_", "\n") for r in primary]
    values = [100 * r["median_fractional_delta"] * r["favorable_direction"] for r in primary]
    fig, ax = plt.subplots(figsize=(10, 5)); ax.bar(labels, values)
    ax.axhline(20, color="k", ls="--", label="material threshold")
    ax.set_ylabel("favorable signed median change (%)"); ax.legend(); fig.tight_layout()
    fig.savefig(out / "primary_effects.png", dpi=160); plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for metric, ax in zip(("zeta_max_s-1", "circulation_4200_m2_s"), axes):
        q = [r for r in paired if r["metric"] == metric]
        ax.plot([r["time_s"] for r in q], [r["control15"] for r in q], label="CONTROL-15")
        ax.plot([r["time_s"] for r in q], [r["hightop20"] for r in q], label="HIGH-TOP-20")
        ax.set_ylabel(metric); ax.legend()
    axes[-1].set_xlabel("time (s)"); fig.tight_layout()
    fig.savefig(out / "vortex_timeseries.png", dpi=160); plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    for arm, rows in observer.items():
        take = [r for r in rows if 2790 <= r["time_s"] <= 3301]
        ax.plot([r["time_s"] for r in take], [r["energy_removed_J"] for r in take], label=arm)
    ax.set_xlabel("time (s)"); ax.set_ylabel("damping energy removed per step (J)"); ax.legend()
    fig.tight_layout(); fig.savefig(out / "damping_energy.png", dpi=160); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--control", type=Path, required=True)
    ap.add_argument("--high", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(); args.out.mkdir(parents=True, exist_ok=False)
    control, gc, _ = read_sequence(args.control / "sequence.h5")
    high, gh, _ = read_sequence(args.high / "sequence.h5")
    if not (gc["Lz_m"] == 15000.0 and gc["nz"] == 48 and gh["Lz_m"] == 20000.0 and gh["nz"] == 54):
        raise RuntimeError(f"unexpected geometries: {gc}, {gh}")
    stats, paired, primary, material, classification = metric_comparison(control, high)
    tracking = tracking_gate(control, high)
    geom_rows = {}; geom = {}; closures = {}; closure_rows = {}; obs = {}; obs_integrity = {}
    for arm, directory in (("control15", args.control), ("hightop20", args.high)):
        geom_rows[arm], geom[arm] = geometric_gate(directory / "sequence.h5")
        closure_rows[arm], closures[arm] = budget_closure(directory / "sequence.h5")
        obs[arm], obs_integrity[arm] = aggregate_top_observer(directory / "top_boundary_calls.h5")
        write_csv(args.out / f"{arm}_geometry.csv", geom_rows[arm])
        write_csv(args.out / f"{arm}_budget_closure.csv", closure_rows[arm])
        write_csv(args.out / f"{arm}_top_step_series.csv", obs[arm])
    lag_rows = []; lag_peaks = []
    for arm in ("control15", "hightop20"):
        for driver in ("energy_removed_J", "pressure_flux_change_W", "curl_driver_s-1"):
            for response in ("low_level_convergence_rms_s-1", "low_level_w_rms_m_s"):
                rows, peak = lag_correlation(obs[arm], driver, response)
                lag_rows.extend({"arm": arm, **r} for r in rows)
                if peak: lag_peaks.append({"arm": arm, **peak})
    write_csv(args.out / "native_diagnostics.csv",
              [{"arm": "control15", **r} for r in control] + [{"arm": "hightop20", **r} for r in high])
    write_csv(args.out / "metric_comparison.csv", stats)
    write_csv(args.out / "paired_timeseries.csv", paired)
    write_csv(args.out / "lag_correlations.csv", lag_rows)
    write_csv(args.out / "lag_peaks.csv", lag_peaks)
    gates = {
        "hightop_geometry": geom["hightop20"], "budget_closure": closures,
        "observer_integrity": obs_integrity, "tracking": tracking,
    }
    gates_pass = (geom["hightop20"]["status"] == "PASS"
                  and all(v["status"] == "PASS" for v in closures.values())
                  and all(v["status"] == "PASS" for v in obs_integrity.values())
                  and tracking["status"] == "PASS")
    final_classification = classification if gates_pass else "UNRESOLVED: one or more gates failed"
    result = {
        "status": "complete", "geometry": {"control15": gc, "hightop20": gh},
        "common_times_s": [r["time_s"] for r in control], "primary_metrics": primary,
        "all_metrics": stats, "materially_favorable_primary_metrics": material,
        "screening_classification_before_gates": classification,
        "classification": final_classification, "all_gates_pass": gates_pass,
        "gates": gates, "lag_peaks": lag_peaks,
        "limitations": [
            "Single deterministic pair; ensemble uncertainty is not estimated.",
            "The low-memory projector imposes w=0 at the top; this is a rigid lid plus internal sponge.",
            "The retained reversed damping ramp and 5.51% layer-thickness difference limit external validity.",
            "Lag correlations show association/propagation proxies, not a modal proof of reflection.",
        ],
    }
    (args.out / "summary.json").write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
    figures(args.out, paired, primary, obs)
    artifacts = {}
    for path in sorted(args.out.iterdir()):
        if path.is_file() and path.name != "manifest.json":
            artifacts[path.name] = {"bytes": path.stat().st_size,
                                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    (args.out / "manifest.json").write_text(json.dumps({"artifacts": artifacts}, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    main()
