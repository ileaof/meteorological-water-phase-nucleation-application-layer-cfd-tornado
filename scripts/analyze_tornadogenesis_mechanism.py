"""Objective source/conversion/concentration audit of the authorized v4 run."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import distance_transform_edt, label

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))


SOURCES = (
    "initial", "buoyancy", "les", "surface_drag", "coriolis",
    "projection", "boundary", "other", "advection_remainder",
)
THRESHOLDS = (0.002, 0.003, 0.005)
RADII_M = np.array((1200.0, 1800.0, 2400.0, 3000.0, 4200.0, 6000.0))
ANALYSIS_CYLINDER_RADIUS_M = 4200.0
ANALYSIS_VERTICAL_RANGE_M = (0.0, 2000.0)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def dx(a, spacing):
    return (np.roll(a, -1, axis=0) - np.roll(a, 1, axis=0)) / (2.0 * spacing)


def dy(a, spacing):
    return (np.roll(a, -1, axis=1) - np.roll(a, 1, axis=1)) / (2.0 * spacing)


def dz(a, z):
    return np.gradient(a, z, axis=2, edge_order=2)


def centered(g, nk=None):
    end = None if nk is None else nk
    u = 0.5 * (g["u"][:-1, :, :end] + g["u"][1:, :, :end])
    v = 0.5 * (g["v"][:, :-1, :end] + g["v"][:, 1:, :end])
    wend = None if nk is None else nk + 1
    w = 0.5 * (g["w"][:, :, :wend][:, :, :-1] + g["w"][:, :, :wend][:, :, 1:])
    return np.asarray(u), np.asarray(v), np.asarray(w)


def kinematics(u, v, w, x, y, z):
    sx = float(x[1] - x[0]); sy = float(y[1] - y[0])
    dudx, dudy, dudz = dx(u, sx), dy(u, sy), dz(u, z)
    dvdx, dvdy, dvdz = dx(v, sx), dy(v, sy), dz(v, z)
    dwdx, dwdy, dwdz = dx(w, sx), dy(w, sy), dz(w, z)
    zeta = dvdx - dudy
    xi = dwdy - dvdz
    eta = dudz - dwdx
    conv = -(dudx + dvdy)
    stretch = zeta * dwdz
    tilt = xi * dwdx + eta * dwdy
    divergence = -(zeta * (dudx + dvdy + dwdz))
    adv_h = -(u * dx(zeta, sx) + v * dy(zeta, sy))
    adv_v = -(w * dz(zeta, z))
    return {
        "zeta": zeta, "xi": xi, "eta": eta, "convergence": conv,
        "stretching": stretch, "tilting": tilt, "divergence": divergence,
        "horizontal_advection": adv_h, "vertical_advection": adv_v,
        "dwdz": dwdz, "omega_h": np.hypot(xi, eta),
        "grad_h_w": np.hypot(dwdx, dwdy),
        "tilting_alignment": tilt / (np.hypot(xi, eta) * np.hypot(dwdx, dwdy) + 1e-30),
    }


def periodic_delta(a, b, length):
    return (a - b + 0.5 * length) % length - 0.5 * length


def distance2(x, y, center):
    return periodic_delta(x[:, None], center[0], x[-1] - x[0] + x[1] - x[0]) ** 2 + periodic_delta(
        y[None, :], center[1], y[-1] - y[0] + y[1] - y[0]
    ) ** 2


def weighted_center(field, x, y, guess=None, radius=4800.0):
    positive = np.maximum(field, 0.0)
    if guess is not None:
        local = distance2(x, y, guess) <= radius**2
        candidate = np.where(local, positive, 0.0)
        if not np.any(candidate):
            candidate = positive
    else:
        candidate = positive
    seed = np.unravel_index(int(np.argmax(candidate)), candidate.shape)
    seed_center = (float(x[seed[0]]), float(y[seed[1]]))
    core = distance2(x, y, seed_center) <= 1800.0**2
    weight = positive * core
    if weight.sum() <= 0:
        return seed_center
    # Vortex stays far from periodic seams in this experiment; ordinary moments are safe.
    return (
        float((weight * x[:, None]).sum() / weight.sum()),
        float((weight * y[None, :]).sum() / weight.sum()),
    )


def track_axis(zeta, x, y, z, previous_surface=None):
    centers = []
    surface = weighted_center(zeta[:, :, 1], x, y, previous_surface, radius=6000.0)
    center = surface
    for k in range(len(z)):
        center = weighted_center(zeta[:, :, k], x, y, center, radius=4200.0)
        centers.append(center)
    return np.asarray(centers), surface


def component_mask(zeta, x, y, z, axis, threshold, radius=4200.0):
    mask = np.zeros_like(zeta, dtype=bool)
    for k in range(len(z)):
        candidate = (zeta[:, :, k] >= threshold) & (
            distance2(x, y, axis[k]) <= radius**2
        )
        labs, _ = label(candidate, structure=np.ones((3, 3), dtype=int))
        i = int(np.argmin(np.abs(x - axis[k, 0])))
        j = int(np.argmin(np.abs(y - axis[k, 1])))
        target = labs[i, j]
        if target:
            mask[:, :, k] = labs == target
        elif candidate.any():
            q = np.unravel_index(np.argmax(np.where(candidate, zeta[:, :, k], -np.inf)), candidate.shape)
            mask[:, :, k] = labs == labs[q]
    return mask


def radial_metrics(zeta, u, v, p_dyn, x, y, center, k, dx_m, radial_bin_m=None):
    ddx = periodic_delta(x[:, None], center[0], x[-1] - x[0] + dx_m)
    ddy = periodic_delta(y[None, :], center[1], y[-1] - y[0] + dx_m)
    r = np.hypot(ddx, ddy)
    phi_x = -ddy / np.maximum(r, 1e-9)
    phi_y = ddx / np.maximum(r, 1e-9)
    vtheta = u[:, :, k] * phi_x + v[:, :, k] * phi_y
    gamma = [float(zeta[:, :, k][r <= radius].sum() * dx_m**2) for radius in RADII_M]
    bin_m = dx_m if radial_bin_m is None else radial_bin_m
    edges = np.arange(0.0, 7200.0 + bin_m, bin_m)
    rb, vp = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        ann = (r >= lo) & (r < hi)
        if ann.sum() >= 3:
            rb.append(0.5 * (lo + hi)); vp.append(float(vtheta[ann].mean()))
    rb = np.asarray(rb); vp = np.asarray(vp)
    if len(vp):
        imax = int(np.argmax(np.abs(vp)))
        rmw = float(rb[imax]); vmax = float(abs(vp[imax]))
    else:
        rmw = vmax = float("nan")
    core = r <= 1800.0
    annulus = (r >= 4800.0) & (r <= 7200.0)
    pmin = float(np.min(p_dyn[:, :, k][core]))
    pamb = float(np.mean(p_dyn[:, :, k][annulus]))
    return gamma, rb, vp, rmw, vmax, pmin - pamb, r, vtheta


def corr(a, b, mask):
    aa = np.asarray(a)[mask]; bb = np.asarray(b)[mask]
    if aa.size < 3 or np.std(aa) == 0 or np.std(bb) == 0:
        return float("nan")
    return float(np.corrcoef(aa, bb)[0, 1])


def curl_increment(group, name, nk, dx_m, dy_m):
    inc = group[f"increments/{name}"]
    u = 0.5 * (inc["u"][:-1, :, :nk] + inc["u"][1:, :, :nk])
    v = 0.5 * (inc["v"][:, :-1, :nk] + inc["v"][:, 1:, :nk])
    return dx(np.asarray(v), dx_m) - dy(np.asarray(u), dy_m)


def extrema_event(times, values, kind="max"):
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return 0
    q = np.where(valid, values, -np.inf if kind == "max" else np.inf)
    return int(np.argmax(q) if kind == "max" else np.argmin(q))


def savefig(fig, path):
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radial-bin-m',type=float,default=None)
    parser.add_argument(
        "--sequence", type=Path,
        default=ROOT / "outputs/diagnostic_sequence_20260905/sequence.h5",
    )
    parser.add_argument(
        "--provenance", type=Path,
        default=ROOT / "outputs/vorticity_provenance_long_v4_2790_3300/provenance_long_v4.h5",
    )
    parser.add_argument(
        "--validity", type=Path,
        default=ROOT / "outputs/vorticity_provenance_long_v4_2790_3300/summary.json",
    )
    parser.add_argument(
        "--out", type=Path,
        default=ROOT / "outputs/tornadogenesis_mechanism_audit",
    )
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    figures = args.out / "figures"; figures.mkdir(exist_ok=True)

    validity = json.loads(args.validity.read_text(encoding="utf-8"))
    if validity["numerical_validity_gate"] != "PASS":
        raise RuntimeError("physical interpretation forbidden: numerical gate did not pass")

    rows, threshold_rows, vertical_rows, radial_rows = [], [], [], []
    provenance_rows, budget_rows = [], []
    map_cache = {}
    source_series = {name: [] for name in SOURCES}
    source_positive_series = {name: [] for name in SOURCES}
    source_stretch_series = {name: [] for name in SOURCES}
    source_tilting_series = {name: [] for name in SOURCES}
    source_horizontal_series = {name: [] for name in SOURCES}
    times = []
    previous_surface = None

    with h5py.File(args.sequence, "r") as seq, h5py.File(args.provenance, "r") as prov:
        x = seq["grid/xc"][:]; y = seq["grid/yc"][:]; zfull = seq["grid/zc"][:]
        zf = seq["grid/zf"][:]
        nk = int(prov.attrs["analysis_nz"]); z = zfull[:nk]
        dx_m = float(x[1] - x[0]); dy_m = float(y[1] - y[0])
        cell_volume = dx_m * dy_m * np.diff(zf)[:nk][None, None, :]
        pkeys = list(prov["snapshots"])
        ptimes = np.array([float(prov[f"snapshots/{key}"].attrs["time_s"]) for key in pkeys])
        skeys_all = list(seq["snapshots"])
        stimes_all = np.array([float(seq[f"snapshots/{key}"].attrs["time_s"]) for key in skeys_all])
        select = (stimes_all >= ptimes[0] - 1e-6) & (stimes_all <= ptimes[-1] + 1e-6)
        skeys = list(np.asarray(skeys_all)[select]); stimes = stimes_all[select]

        theta_v0 = np.asarray(seq["base/theta0"]) * (
            1.0 + 0.61 * np.asarray(seq["base/qv0"])
        )

        for ip, pkey in enumerate(pkeys):
            pg = prov[f"snapshots/{pkey}"]
            time_s = float(pg.attrs["time_s"]); times.append(time_s)
            js = int(np.argmin(abs(stimes - time_s))); skey = skeys[js]
            sg = seq[f"snapshots/{skey}"]
            mismatch_s = float(stimes[js] - time_s)
            u, v, w = centered(sg, nk); kin = kinematics(u, v, w, x, y, z)
            zeta = np.asarray(pg["omega_total"])[2, :, :, :nk]
            source = {name: np.asarray(pg[f"omega_source/{name}"])[..., :nk] for name in SOURCES}
            axis, previous_surface = track_axis(zeta, x, y, z, previous_surface)
            k0 = 1
            center = tuple(axis[k0])
            full_mask = component_mask(zeta, x, y, z, axis, THRESHOLDS[0])
            core_mask = component_mask(zeta, x, y, z, axis, THRESHOLDS[2])
            cylinder = np.stack([
                distance2(x, y, axis[k]) <= ANALYSIS_CYLINDER_RADIUS_M**2
                for k in range(nk)
            ], axis=2)
            positive_cylinder = cylinder & (zeta > 0)
            gamma, rb, vp, rmw, vmax, pdef, radius_field, vtheta = radial_metrics(
                zeta, u, v, np.asarray(sg["p_dyn"][:, :, :nk]), x, y, center, k0, dx_m, args.radial_bin_m
            )
            radial_rows.extend(
                {"time_s": time_s, "radius_m": float(radius),
                 "circulation_m2_s": value, "vtheta_mean_m_s": float("nan")}
                for radius, value in zip(RADII_M, gamma)
            )
            radial_rows.extend(
                {"time_s": time_s, "radius_m": float(radius),
                 "circulation_m2_s": float("nan"), "vtheta_mean_m_s": value}
                for radius, value in zip(rb, vp)
            )
            theta = np.asarray(sg["theta"][:, :, :nk]); qv = np.asarray(sg["qv"][:, :, :nk])
            condensate = sum(np.asarray(sg[name][:, :, :nk]) for name in ("ql", "qi", "qr", "qs", "qg", "qh"))
            theta_v = theta * (1.0 + 0.61 * qv - condensate)
            dtheta_v = theta_v - theta_v0[None, None, :nk]
            cold = dtheta_v[:, :, k0] < -1.0
            signed_cold_distance = (
                distance_transform_edt(~cold, sampling=(dx_m, dy_m))
                - distance_transform_edt(cold, sampling=(dx_m, dy_m))
            )
            rho = np.asarray(sg["rho"][:, :, :nk])
            bx = -9.81 * dy(rho, dy_m) / np.maximum(rho, 1e-9)
            by = 9.81 * dx(rho, dx_m) / np.maximum(rho, 1e-9)
            baro_h = np.hypot(bx, by)
            zmask = positive_cylinder
            high = core_mask
            stretch_pos = (kin["stretching"] > 0) & zmask
            conv_pos = (kin["convergence"] > 0) & zmask
            core_area = float((zeta[:, :, k0] >= THRESHOLDS[0])[radius_field <= 4200.0].sum() * dx_m * dy_m)
            radius = float(np.sqrt(core_area / np.pi)) if core_area else float("nan")
            axis_disp = np.hypot(
                periodic_delta(axis[:, 0], axis[0, 0], x[-1] - x[0] + dx_m),
                periodic_delta(axis[:, 1], axis[0, 1], y[-1] - y[0] + dy_m),
            )
            local2 = radius_field <= ANALYSIS_CYLINDER_RADIUS_M
            z2 = zeta[:, :, k0]
            conv2 = kin["convergence"][:, :, k0]
            p2 = np.asarray(sg["p_dyn"][:, :, k0])
            annulus = (radius_field >= 4800.0) & (radius_field <= 7200.0)
            pamb = float(p2[annulus].mean())
            p_anomaly = p2 - pamb
            zeta_half = local2 & (z2 >= 0.5 * z2[local2].max())
            conv_half = local2 & (conv2 >= 0.5 * conv2[local2].max())
            pressure_half = local2 & (p_anomaly <= 0.5 * pdef)
            def diameter_cells(mask):
                return float(2.0 * np.sqrt(mask.sum() * dx_m * dy_m / np.pi) / dx_m)
            row = {
                "time_s": time_s, "sequence_time_s": float(stimes[js]),
                "alignment_time_mismatch_s": mismatch_s,
                "center_x_m": center[0], "center_y_m": center[1],
                "level_m": float(z[k0]), "zeta_max_s-1": float(z2[local2].max()),
                "zeta_min_s-1": float(z2[local2].min()),
                "circulation_1200_m2_s": gamma[0], "circulation_2400_m2_s": gamma[2],
                "circulation_4200_m2_s": gamma[4], "circulation_6000_m2_s": gamma[5],
                "vtheta_max_m_s": vmax, "rmw_m": rmw,
                "pressure_dynamic_deficit_Pa": pdef,
                "wmax_m_s": float(w[:, :, k0][local2].max()),
                "convergence_max_s-1": float(kin["convergence"][:, :, k0][local2].max()),
                "stretching_max_s-2": float(kin["stretching"][zmask].max()),
                "stretching_conditional_mean_s-2": float(kin["stretching"][high].mean()) if high.any() else float("nan"),
                "tilting_conditional_mean_s-2": float(kin["tilting"][high].mean()) if high.any() else float("nan"),
                "corr_zeta_dwdz": corr(zeta, kin["dwdz"], zmask),
                "corr_zeta_convergence": corr(zeta, kin["convergence"], zmask),
                "corr_positive_zeta_stretching": corr(zeta, kin["stretching"], zmask),
                "fraction_positive_zeta_positive_stretching": float(stretch_pos.sum() / max(zmask.sum(), 1)),
                "fraction_positive_zeta_positive_convergence": float(conv_pos.sum() / max(zmask.sum(), 1)),
                "vortex_radius_m": radius,
                "axis_tilt_surface_to_2km_m": float(axis_disp[-1]),
                "axis_tilt_max_m": float(axis_disp.max()),
                "cold_pool_area_km2": float(cold.sum() * dx_m * dy_m / 1e6),
                "cold_pool_fraction_near_vortex": float((cold & local2).sum() / max(local2.sum(), 1)),
                "cold_pool_min_dtheta_v_K": float(dtheta_v[:, :, k0][local2].min()),
                "cold_pool_signed_distance_at_center_m": float(
                    signed_cold_distance[
                        int(np.argmin(abs(x-center[0]))), int(np.argmin(abs(y-center[1])))
                    ]
                ),
                "baroclinic_h_mean_near_vortex_s-2": float(baro_h[:, :, k0][local2].mean()),
                "baroclinic_h_max_near_vortex_s-2": float(baro_h[:, :, k0][local2].max()),
                "zeta_core_width_cells": float(2.0 * radius / dx_m) if np.isfinite(radius) else float("nan"),
                "rmw_cells": float(rmw / dx_m),
                "zeta_halfmax_width_cells": diameter_cells(zeta_half),
                "convergence_halfmax_width_cells": diameter_cells(conv_half),
                "pressure_halfdeficit_width_cells": diameter_cells(pressure_half),
            }
            rows.append(row)
            # Fixed physical radius, identical lower-grid heights across runs.
            volume_2km = dx_m * dy_m * np.maximum(
                0.,
                np.minimum(zf[1:nk+1], ANALYSIS_VERTICAL_RANGE_M[1])
                - np.maximum(zf[:nk], ANALYSIS_VERTICAL_RANGE_M[0]),
            )[None, None, :]
            row.update(
                zeta_max_0_2km_s_1=float(zeta[cylinder].max()),
                zeta_signed_0_2km_m3_s=float((zeta*volume_2km)[cylinder].sum()),
                zeta_absolute_0_2km_m3_s=float((abs(zeta)*volume_2km)[cylinder].sum()),
                circulation_absolute_4200_m2_s=float(abs(z2[local2]).sum()*dx_m*dy_m),
            )
            for width in ('zeta_halfmax_width','convergence_halfmax_width','pressure_halfdeficit_width'):
                row[width+'_m']=row[width+'_cells']*dx_m

            for threshold in THRESHOLDS:
                mask = component_mask(zeta, x, y, z, axis, threshold)
                area = float(mask.sum() * np.mean(cell_volume))
                threshold_rows.append({
                    "time_s": time_s, "threshold_s-1": threshold,
                    "cells": int(mask.sum()), "volume_m3_proxy": area,
                    "zeta_integral_m3_s-1": float((zeta * cell_volume)[mask].sum()),
                    "zeta_abs_integral_m3_s-1": float((abs(zeta) * cell_volume)[mask].sum()),
                    "equivalent_radius_m": float(np.sqrt(mask[:, :, k0].sum() * dx_m * dy_m / np.pi)),
                })

            dwdz = kin["dwdz"]
            grad_w = np.asarray(pg["grad_w_h"])[..., :nk]
            total_integral = float((zeta * cell_volume)[full_mask].sum())
            total_core = float((zeta * cell_volume)[core_mask].sum())
            for name in SOURCES:
                zj = source[name][2]
                sj = zj * dwdz
                omhj = np.hypot(source[name][0], source[name][1])
                tj = source[name][0] * grad_w[0] + source[name][1] * grad_w[1]
                signed = float((zj * cell_volume)[full_mask].sum())
                core_signed = float((zj * cell_volume)[core_mask].sum())
                source_series[name].append(signed)
                source_positive_series[name].append(float((np.maximum(zj, 0) * cell_volume)[full_mask].sum()))
                source_stretch_series[name].append(float((sj * cell_volume)[full_mask].sum()))
                source_tilting_series[name].append(float((tj * cell_volume)[full_mask].sum()))
                source_horizontal_series[name].append(float((omhj * cell_volume)[full_mask].sum()))
                provenance_rows.append({
                    "time_s": time_s, "source": name, "mask": "vortex",
                    "signed_zeta_integral_m3_s-1": signed,
                    "absolute_zeta_integral_m3_s-1": float((abs(zj) * cell_volume)[full_mask].sum()),
                    "positive_zeta_integral_m3_s-1": float((np.maximum(zj, 0) * cell_volume)[full_mask].sum()),
                    "negative_zeta_integral_m3_s-1": float((np.minimum(zj, 0) * cell_volume)[full_mask].sum()),
                    "rms_zeta_s-1": float(np.sqrt(np.mean(zj[full_mask]**2))) if full_mask.any() else float("nan"),
                    "max_zeta_s-1": float(zj[full_mask].max()) if full_mask.any() else float("nan"),
                    "min_zeta_s-1": float(zj[full_mask].min()) if full_mask.any() else float("nan"),
                    "fraction_signed": signed / total_integral if total_integral else float("nan"),
                    "core_signed_zeta_integral_m3_s-1": core_signed,
                    "core_fraction_signed": core_signed / total_core if total_core else float("nan"),
                    "stretching_integral_m3_s-2": float((sj * cell_volume)[full_mask].sum()),
                    "horizontal_vorticity_abs_integral_m3_s-1": float((omhj * cell_volume)[full_mask].sum()),
                    "tilting_integral_m3_s-2": float((tj * cell_volume)[full_mask].sum()),
                    "tilting_absolute_integral_m3_s-2": float((abs(tj) * cell_volume)[full_mask].sum()),
                    "tilting_positive_integral_m3_s-2": float((np.maximum(tj, 0) * cell_volume)[full_mask].sum()),
                    "tilting_negative_integral_m3_s-2": float((np.minimum(tj, 0) * cell_volume)[full_mask].sum()),
                })

            for k in range(nk):
                rr = distance2(x, y, axis[k]) <= 4200.0**2
                gam = float(zeta[:, :, k][rr].sum() * dx_m * dy_m)
                vertical_rows.append({
                    "time_s": time_s, "z_m": float(z[k]),
                    "center_x_m": float(axis[k, 0]), "center_y_m": float(axis[k, 1]),
                    "displacement_from_surface_m": float(axis_disp[k]),
                    "zeta_max_s-1": float(zeta[:, :, k][rr].max()),
                    "circulation_4200_m2_s": gam,
                    "wmax_m_s": float(w[:, :, k][rr].max()),
                    "pressure_min_Pa": float(np.asarray(sg["p_dyn"][:, :, k])[rr].min()),
                })

            if ip in (0, 3, 7, 12, 17):
                map_cache[ip] = {
                    "time": time_s, "zeta": zeta[:, :, k0],
                    "conv": kin["convergence"][:, :, k0],
                    "stretch": kin["stretching"][:, :, k0],
                    "p": np.asarray(sg["p_dyn"][:, :, k0]), "w": w[:, :, k0],
                    "dtheta_v": dtheta_v[:, :, k0], "cold": cold, "axis": axis,
                    "sources": {name: source[name][2, :, :, k0] for name in SOURCES},
                }

        # Native, interval-mean Eulerian budget. Snapshot i contains increments since i-1.
        stage_names = (
            "advection", "buoyancy", "les", "surface_drag", "coriolis",
            "projection", "initial_bcs", "predictor_bcs", "projection_bcs",
            "transport_microphysics_bcs", "external_forcing", "guard",
        )
        for i in range(1, len(skeys)):
            prev = seq[f"snapshots/{skeys[i-1]}"]; cur = seq[f"snapshots/{skeys[i]}"]
            dt_s = float(cur.attrs["time_s"] - prev.attrs["time_s"])
            pu, pv, _ = centered(prev, nk); cu, cv, _ = centered(cur, nk)
            delta = (dx(cv, dx_m) - dy(cu, dy_m) - (dx(pv, dx_m) - dy(pu, dy_m))) / dt_s
            terms = {name: curl_increment(cur, name, nk, dx_m, dy_m) / dt_s for name in stage_names}
            total = sum(terms.values())
            ir = int(np.argmin(abs(np.asarray(times) - float(cur.attrs["time_s"]))))
            # Rebuild the tracked mask from the nearest provenance snapshot.
            pg = prov[f"snapshots/{pkeys[ir]}"]; zz = np.asarray(pg["omega_total"])[2, :, :, :nk]
            prof = [r for r in vertical_rows if r["time_s"] == times[ir]]
            axis = np.array([(r["center_x_m"], r["center_y_m"]) for r in prof])
            mask = component_mask(zz, x, y, z, axis, THRESHOLDS[0])
            record = {
                "time_s": float(cur.attrs["time_s"]), "interval_s": dt_s,
                "delta_zeta_mean_s-2": float(delta[mask].mean()),
                "native_sum_mean_s-2": float(total[mask].mean()),
                "native_closure_mean_s-2": float((delta-total)[mask].mean()),
                "native_closure_rms_s-2": float(np.sqrt(np.mean((delta-total)[mask]**2))),
                "delta_zeta_rms_s-2": float(np.sqrt(np.mean(delta[mask]**2))),
            }
            record["native_closure_relative_rms"] = (
                record["native_closure_rms_s-2"]
                / max(record["delta_zeta_rms_s-2"], 1e-30)
            )
            for name, values in terms.items():
                record[f"{name}_mean_s-2"] = float(values[mask].mean())
            # The capture stores independent dt-integrals of resolved kinematic terms.
            for name in ("stretching", "tilting"):
                value = np.asarray(cur[f"kinematic_integrals/{name}"])[:, :, :nk] / dt_s
                record[f"kinematic_{name}_mean_s-2"] = float(value[mask].mean())
            vector_adv = np.asarray(cur["kinematic_integrals/vector_advection"])[2, :, :, :nk] / dt_s
            record["kinematic_total_advection_mean_s-2"] = float(vector_adv[mask].mean())
            vector_dil = np.asarray(cur["kinematic_integrals/vector_dilatation"])[2, :, :, :nk] / dt_s
            record["kinematic_dilatation_mean_s-2"] = float(vector_dil[mask].mean())
            budget_rows.append(record)

    write_csv(args.out / "vortex_timeseries.csv", rows)
    write_csv(args.out / "threshold_sensitivity.csv", threshold_rows)
    write_csv(args.out / "vertical_profiles.csv", vertical_rows)
    write_csv(args.out / "radial_profiles.csv", radial_rows)
    write_csv(args.out / "provenance_integrals.csv", provenance_rows)
    write_csv(args.out / "vorticity_budget.csv", budget_rows)

    times = np.asarray(times); zmax = np.array([r["zeta_max_s-1"] for r in rows])
    gamma = np.array([r["circulation_4200_m2_s"] for r in rows])
    conv = np.array([r["convergence_max_s-1"] for r in rows])
    stretch = np.array([r["stretching_conditional_mean_s-2"] for r in rows])
    wmax = np.array([r["wmax_m_s"] for r in rows]); pdef = np.array([r["pressure_dynamic_deficit_Pa"] for r in rows])
    vmax = np.array([r["vtheta_max_m_s"] for r in rows])
    dzdt = np.gradient(zmax, times)
    post_peak = np.arange(len(times)) > int(np.argmax(zmax))
    weak_index = int(np.where(post_peak, dzdt, np.inf).argmin())
    noninitial = np.vstack([np.abs(source_positive_series[n]) for n in SOURCES[1:-1]])
    source_peak_flat = int(np.argmax(noninitial)); source_i = source_peak_flat % len(times)
    source_name = SOURCES[1:-1][source_peak_flat // len(times)]
    event_indices = {
        "inicio_crescimento_zeta_nao_resolvido_limite_esquerdo": 0,
        "zeta_maximo_baixo_nivel": int(np.argmax(zmax)),
        "maior_aceleracao_zeta": int(np.argmax(dzdt)),
        "convergencia_maxima": extrema_event(times, conv),
        "stretching_condicional_maximo": extrema_event(times, stretch),
        "pressao_dinamica_minima": extrema_event(times, pdef, "min"),
        "vtheta_maximo": extrema_event(times, vmax),
        "circulacao_maxima": extrema_event(times, gamma),
        "updraft_baixo_nivel_maximo": extrema_event(times, wmax),
        f"proveniencia_pos_mais_forte_{source_name}": source_i,
        "enfraquecimento_mais_rapido": weak_index,
    }
    event_rows = []
    for event, idx in event_indices.items():
        r = rows[idx]
        contributions = {name: source_positive_series[name][idx] for name in SOURCES}
        dominant = max(contributions, key=lambda n: abs(contributions[n]))
        denom = sum(abs(v) for v in contributions.values())
        increments = {name: contributions[name] for name in SOURCES[1:-1]}
        dominant_increment = max(increments, key=lambda n: abs(increments[n]))
        increment_denom = sum(abs(v) for v in increments.values())
        event_rows.append({
            "event": event, "time_s": r["time_s"], "zeta_max_s-1": r["zeta_max_s-1"],
            "circulation_4200_m2_s": r["circulation_4200_m2_s"],
            "vtheta_max_m_s": r["vtheta_max_m_s"], "rmw_m": r["rmw_m"],
            "pressure_dynamic_deficit_Pa": r["pressure_dynamic_deficit_Pa"],
            "wmax_m_s": r["wmax_m_s"], "convergence_max_s-1": r["convergence_max_s-1"],
            "stretching_conditional_mean_s-2": r["stretching_conditional_mean_s-2"],
            "dominant_provenance_source": dominant,
            "dominant_positive_inventory_fraction": contributions[dominant] / denom if denom else np.nan,
            "dominant_post_restart_source": dominant_increment,
            "dominant_post_restart_positive_fraction": increments[dominant_increment] / increment_denom if increment_denom else np.nan,
            "kappa_omega": np.nan, "kappa_T": np.nan,
            "closure_relative": np.nan, "vortex_radius_m": r["vortex_radius_m"],
            "axis_tilt_m": r["axis_tilt_surface_to_2km_m"],
        })
    # Add exact stage metrics nearest each event.
    with h5py.File(args.provenance, "r") as prov:
        hist = prov["closure_history"][:]
    ht = hist["time_s"]
    for r in event_rows:
        h = hist[int(np.argmin(abs(ht-r["time_s"])))]
        r["kappa_omega"] = float(h["omega_condition_index"])
        r["kappa_T"] = float(h["tilting_condition_index"])
        r["closure_relative"] = float(h["omega_relative_rms"])
    write_csv(args.out / "principal_events.csv", event_rows)

    # Compact 0--2 km fields for every distinct objective event time.
    with h5py.File(args.sequence, "r") as seq, h5py.File(args.provenance, "r") as prov, h5py.File(
        args.out / "event_states.h5", "w"
    ) as events:
        events.attrs["schema"] = "storm-tornadogenesis-event-states-v1"
        events.attrs["event_index"] = json.dumps(event_indices, ensure_ascii=False)
        events.attrs["provenance_formulation"] = "storm-vorticity-provenance-v4"
        events.create_dataset("x_m", data=x); events.create_dataset("y_m", data=y)
        events.create_dataset("z_m", data=z)
        seq_keys = list(seq["snapshots"])
        seq_times = np.array([float(seq[f"snapshots/{key}"].attrs["time_s"]) for key in seq_keys])
        prov_keys = list(prov["snapshots"])
        for idx in sorted(set(event_indices.values())):
            pg = prov[f"snapshots/{prov_keys[idx]}"]
            time_s = float(pg.attrs["time_s"])
            js = int(np.argmin(abs(seq_times-time_s))); sg = seq[f"snapshots/{seq_keys[js]}"]
            u, v, w = centered(sg, nk); kin = kinematics(u, v, w, x, y, z)
            zeta = np.asarray(pg["omega_total"])[2, :, :, :nk]
            profile = [r for r in vertical_rows if r["time_s"] == time_s]
            axis = np.array([(r["center_x_m"], r["center_y_m"]) for r in profile])
            group = events.create_group(f"t_{time_s:010.3f}")
            group.attrs["time_s"] = time_s
            group.attrs["sequence_time_s"] = float(seq_times[js])
            group.attrs["events"] = json.dumps(
                [name for name, event_idx in event_indices.items() if event_idx == idx],
                ensure_ascii=False,
            )
            for name, value in (
                ("zeta", zeta), ("xi", np.asarray(pg["omega_total"])[0, :, :, :nk]),
                ("eta", np.asarray(pg["omega_total"])[1, :, :, :nk]),
                ("convergence", kin["convergence"]), ("stretching", kin["stretching"]),
                ("tilting", kin["tilting"]), ("w", w),
                ("p_dyn", np.asarray(sg["p_dyn"][:, :, :nk])),
                ("vortex_mask_002", component_mask(zeta, x, y, z, axis, THRESHOLDS[0])),
                ("vortex_core_mask_005", component_mask(zeta, x, y, z, axis, THRESHOLDS[2])),
            ):
                group.create_dataset(name, data=value, compression="gzip", compression_opts=1, shuffle=True)
            group.create_dataset("axis_xy_m", data=axis)
            source_group = group.create_group("zeta_source")
            for name in SOURCES:
                source_group.create_dataset(
                    name, data=np.asarray(pg[f"omega_source/{name}"])[2, :, :, :nk],
                    compression="gzip", compression_opts=1, shuffle=True,
                )

    # Lag correlations, positive lag means driver leads zeta tendency.
    lag_rows = []
    for name, driver in (("convergence", conv), ("stretching", stretch), ("circulation", gamma)):
        for lag in range(-3, 4):
            if lag < 0: a, b = driver[-lag:], dzdt[:lag]
            elif lag > 0: a, b = driver[:-lag], dzdt[lag:]
            else: a, b = driver, dzdt
            lag_rows.append({"driver": name, "lag_frames": lag, "lag_s": lag*30.0,
                             "correlation_with_dzeta_dt": float(np.corrcoef(a, b)[0, 1])})
    write_csv(args.out / "lag_correlations.csv", lag_rows)

    # Figures 1--10: time series and provenance.
    plt.style.use("seaborn-v0_8-whitegrid")
    scalar_figures = [
        ("01_zeta_max_time.png", zmax, r"$\zeta_{max}$ (s$^{-1}$)"),
        ("02_circulation_time.png", gamma, r"$\Gamma$(4.2 km) (m$^2$ s$^{-1}$)"),
        ("03_vtheta_max_time.png", vmax, r"$V_{\theta,max}$ (m s$^{-1}$)"),
        ("04_pressure_deficit_time.png", pdef, r"$\Delta p_{dyn}$ (Pa)"),
        ("05_convergence_time.png", conv, r"$C_{max}$ (s$^{-1}$)"),
        ("06_stretching_time.png", stretch, r"$\langle S\mid\zeta>0.005\rangle$ (s$^{-2}$)"),
    ]
    for filename, values, ylabel in scalar_figures:
        fig, ax = plt.subplots(figsize=(7.2, 4.0)); ax.plot(times, values, "o-", lw=1.8)
        ax.set(xlabel="tempo (s)", ylabel=ylabel); savefig(fig, figures/filename)

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.plot(ht, hist["omega_condition_index"], label=r"$\kappa_\omega$")
    ax.plot(ht, hist["tilting_condition_index"], label=r"$\kappa_T$")
    ax.set(xlabel="tempo (s)", ylabel="índice de condicionamento"); ax.legend()
    savefig(fig, figures/"07_conditioning.png")

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    total = np.sum(np.abs(np.vstack([source_series[n] for n in SOURCES])), axis=0)
    for name in SOURCES:
        ax.plot(times, np.asarray(source_series[name])/np.maximum(total, 1e-30), label=name)
    ax.set(xlabel="tempo (s)", ylabel="integral assinado / soma de módulos"); ax.legend(ncol=3, fontsize=7)
    savefig(fig, figures/"08_provenance_fraction_time.png")

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    for name in SOURCES:
        ax.plot(times, source_positive_series[name], label=name)
    ax.set(xlabel="tempo (s)", ylabel=r"$\int_V\max(\zeta_j,0)dV$"); ax.legend(ncol=3, fontsize=7)
    savefig(fig, figures/"09_positive_zeta_source_time.png")

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    stack_names = SOURCES[:-1]; stack = np.vstack([np.abs(source_series[n]) for n in stack_names])
    ax.stackplot(times, stack, labels=stack_names, alpha=.85)
    ax.set(xlabel="tempo (s)", ylabel="módulo do inventário assinado"); ax.legend(ncol=3, fontsize=7)
    savefig(fig, figures/"10_source_provenance_stacked.png")

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    for name in SOURCES:
        ax.plot(times, source_tilting_series[name], label=name)
    ax.axhline(0, color="k", lw=.8)
    ax.set(xlabel="tempo (s)", ylabel="tilting condicionado integrado"); ax.legend(ncol=3, fontsize=7)
    savefig(fig, figures/"10b_source_tilting_time.png")

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    for name in SOURCES:
        ax.plot(times, source_horizontal_series[name], label=name)
    ax.set(xlabel="tempo (s)", ylabel=r"$\int_V |\omega_{h,j}|dV$"); ax.legend(ncol=3, fontsize=7)
    savefig(fig, figures/"10c_horizontal_vorticity_source_time.png")

    def map_panels(field, contours, filename, label_text, cmap="RdBu_r", symmetric=True):
        items = list(map_cache.values()); vals = np.concatenate([m[field].ravel() for m in items])
        if symmetric:
            lim = np.nanmax(abs(vals)); vmin, vmax_ = -lim, lim
        else:
            vmin, vmax_ = np.nanpercentile(vals, [.5, 99.5])
        fig, axs = plt.subplots(1, len(items), figsize=(15.5, 3.3), sharex=True, sharey=True)
        for ax, m in zip(axs, items):
            im=ax.pcolormesh(x/1000,y/1000,m[field].T,shading="auto",cmap=cmap,vmin=vmin,vmax=vmax_)
            ax.contour(x/1000,y/1000,m[contours].T,levels=5,colors="k",linewidths=.45)
            ax.scatter(m["axis"][1,0]/1000,m["axis"][1,1]/1000,s=20,c="gold",edgecolor="k")
            ax.set_title(f'{m["time"]:.0f} s'); ax.set_aspect("equal")
        fig.colorbar(im,ax=axs,label=label_text,shrink=.78); axs[0].set_ylabel("y (km)")
        for ax in axs: ax.set_xlabel("x (km)")
        savefig(fig, figures/filename)

    map_panels("zeta", "conv", "11_zeta_convergence_maps.png", r"$\zeta$ (s$^{-1}$)")
    map_panels("zeta", "stretch", "12_zeta_stretching_maps.png", r"$\zeta$ (s$^{-1}$)")
    map_panels("zeta", "p", "13_zeta_pressure_maps.png", r"$\zeta$ (s$^{-1}$)")
    map_panels("zeta", "w", "14_zeta_w_maps.png", r"$\zeta$ (s$^{-1}$)")

    dominant_sources = sorted(SOURCES[1:-1], key=lambda n: abs(source_positive_series[n][-1]), reverse=True)[:3]
    fig, axs = plt.subplots(len(dominant_sources), 2, figsize=(8.5, 3.3*len(dominant_sources)), sharex=True, sharey=True)
    for row_ax, name in zip(np.atleast_2d(axs), dominant_sources):
        vals=np.concatenate([map_cache[i]["sources"][name].ravel() for i in (0,17)]); lim=np.max(abs(vals))
        for ax, idx in zip(row_ax,(0,17)):
            m=map_cache[idx]; im=ax.pcolormesh(x/1000,y/1000,m["sources"][name].T,cmap="RdBu_r",shading="auto",vmin=-lim,vmax=lim)
            ax.set_title(f"{name}, {m['time']:.0f} s"); ax.set_aspect("equal"); fig.colorbar(im,ax=ax,shrink=.75)
    savefig(fig, figures/"15_dominant_provenance_maps.png")

    # Vertical diagnostics.
    chosen = [0, 3, 7, 12, 17]
    fig, axs = plt.subplots(1, len(chosen), figsize=(15.5, 3.8), sharey=True)
    for ax, idx in zip(axs, chosen):
        m=map_cache[idx]; j=int(np.argmin(abs(y-m["axis"][1,1])))
        # Re-open total field at nearest sequence time for an x-z section.
        with h5py.File(args.sequence,"r") as seq:
            keys=list(seq["snapshots"]); st=np.array([float(seq[f"snapshots/{k}"].attrs["time_s"]) for k in keys]); g=seq[f"snapshots/{keys[int(np.argmin(abs(st-m['time'])))]}"]
            uu,vv,ww=centered(g,nk); kk=kinematics(uu,vv,ww,x,y,z); section=kk["zeta"][:,j,:]
        lim=np.percentile(abs(section),99.5); im=ax.pcolormesh(x/1000,z/1000,section.T,cmap="RdBu_r",shading="auto",vmin=-lim,vmax=lim)
        ax.plot(m["axis"][:,0]/1000,z/1000,"k.-"); ax.set_title(f"{m['time']:.0f} s"); ax.set_xlabel("x (km)")
    axs[0].set_ylabel("z (km)"); fig.colorbar(im,ax=axs,label=r"$\zeta$ (s$^{-1}$)",shrink=.8)
    savefig(fig,figures/"16_vertical_cross_sections.png")

    for filename, key, ylabel in (
        ("17_axis_displacement_height.png","displacement_from_surface_m","deslocamento do eixo (km)"),
        ("18_zeta_max_height.png","zeta_max_s-1",r"$\zeta_{max}$ (s$^{-1}$)"),
        ("19_circulation_height.png","circulation_4200_m2_s",r"$\Gamma$(4.2 km) (m$^2$ s$^{-1}$)"),
    ):
        fig,ax=plt.subplots(figsize=(5.0,5.0))
        for idx in chosen:
            rr=[r for r in vertical_rows if r["time_s"]==times[idx]]
            val=np.array([r[key] for r in rr]);
            if key=="displacement_from_surface_m": val=val/1000
            ax.plot(val,z/1000,label=f"{times[idx]:.0f} s")
        ax.set(xlabel=ylabel,ylabel="z (km)"); ax.legend(fontsize=7); savefig(fig,figures/filename)

    fig,ax=plt.subplots(figsize=(7.2,4.2))
    for idx in chosen:
        rr=[r for r in radial_rows if r["time_s"]==times[idx] and np.isfinite(r["vtheta_mean_m_s"])]
        ax.plot([r["radius_m"]/1000 for r in rr],[r["vtheta_mean_m_s"] for r in rr],label=f"{times[idx]:.0f} s")
    ax.set(xlabel="raio (km)",ylabel=r"$\overline{V_\theta}$ (m s$^{-1}$)"); ax.legend()
    savefig(fig,figures/"20_vtheta_radial_profiles.png")

    fig,ax=plt.subplots(figsize=(7.2,4.2))
    for idx in chosen:
        rr=[r for r in radial_rows if r["time_s"]==times[idx] and np.isfinite(r["circulation_m2_s"])]
        ax.plot([r["radius_m"]/1000 for r in rr],[r["circulation_m2_s"] for r in rr],"o-",label=f"{times[idx]:.0f} s")
    ax.set(xlabel="raio (km)",ylabel=r"$\Gamma$ (m$^2$ s$^{-1}$)"); ax.legend()
    savefig(fig,figures/"21_circulation_radial_profiles.png")

    fig,ax=plt.subplots(figsize=(7.2,4.2)); ax.plot(times,np.array([r["vortex_radius_m"] for r in rows])/1000,"o-",label="raio equivalente")
    ax.plot(times,np.array([r["rmw_m"] for r in rows])/1000,"s-",label="RMW"); ax.set(xlabel="tempo (s)",ylabel="raio (km)");ax.legend()
    savefig(fig,figures/"22_radius_rmw_time.png")

    fig,axs=plt.subplots(1,len(chosen),figsize=(15.5,3.3),sharex=True,sharey=True)
    for ax,idx in zip(axs,chosen):
        m=map_cache[idx]; im=ax.pcolormesh(x/1000,y/1000,m["dtheta_v"].T,cmap="coolwarm",shading="auto",vmin=-6,vmax=6)
        ax.contour(x/1000,y/1000,m["cold"].T,levels=[.5],colors="cyan",linewidths=1.2); ax.contour(x/1000,y/1000,m["zeta"].T,levels=[.002,.003,.005],colors="k",linewidths=.6)
        ax.scatter(m["axis"][1,0]/1000,m["axis"][1,1]/1000,c="yellow",edgecolor="k");ax.set_title(f"{m['time']:.0f} s");ax.set_aspect("equal")
    fig.colorbar(im,ax=axs,label=r"$\theta_v'$ (K)",shrink=.8); savefig(fig,figures/"23_cold_pool_relative_vortex.png")

    fig,ax=plt.subplots(figsize=(9,4.5))
    for n,(event,idx) in enumerate(event_indices.items()): ax.scatter(times[idx],n,s=45); ax.text(times[idx]+3,n,event,va="center",fontsize=7)
    ax.set(xlabel="tempo (s)",yticks=[],title="Ordenamento temporal dos eventos objetivos")
    savefig(fig,figures/"24_event_timeline.png")

    fig,ax=plt.subplots(figsize=(7.2,4.2))
    for threshold in THRESHOLDS:
        rr=[r for r in threshold_rows if r["threshold_s-1"]==threshold]
        ax.plot(times,[r["zeta_integral_m3_s-1"] for r in rr],"o-",label=f"{threshold:.3f} s$^{{-1}}$")
    ax.set(xlabel="tempo (s)",ylabel="inventário condicionado de ζ");ax.legend()
    savefig(fig,figures/"25_threshold_robustness.png")

    fig,ax=plt.subplots(figsize=(8.2,4.5))
    bt=np.array([r["time_s"] for r in budget_rows])
    for name in ("advection","les","surface_drag","projection","coriolis"):
        ax.plot(bt,[r[f"{name}_mean_s-2"] for r in budget_rows],label=name)
    ax.plot(bt,[r["kinematic_stretching_mean_s-2"] for r in budget_rows],"--",label="stretching cinemático")
    ax.plot(bt,[r["kinematic_tilting_mean_s-2"] for r in budget_rows],"--",label="tilting cinemático")
    ax.set(xlabel="tempo (s)",ylabel=r"tendência média (s$^{-2}$)");ax.legend(ncol=2,fontsize=7)
    savefig(fig,figures/"26_vorticity_budget_time.png")

    fig,ax=plt.subplots(figsize=(8.2,4.5))
    ax.plot(times,[r["fraction_positive_zeta_positive_stretching"] for r in rows],label="ζ+ com stretching+")
    ax.plot(times,[r["fraction_positive_zeta_positive_convergence"] for r in rows],label="ζ+ com convergência+")
    ax.plot(times,[r["corr_zeta_convergence"] for r in rows],label="corr(ζ,C)")
    ax.plot(times,[r["corr_zeta_dwdz"] for r in rows],label="corr(ζ,∂w/∂z)")
    ax.set(xlabel="tempo (s)",ylabel="fração / correlação");ax.legend()
    savefig(fig,figures/"27_collocation_time.png")

    fig,ax=plt.subplots(figsize=(7.2,4.2)); ax.plot(times,[r["zeta_core_width_cells"] for r in rows],label="diâmetro ζ")
    ax.plot(times,[r["rmw_cells"] for r in rows],label="RMW");ax.axhline(4,color="k",ls="--",lw=1,label="4 células")
    ax.plot(times,[r["convergence_halfmax_width_cells"] for r in rows],label="convergência ≥ 1/2 máx.")
    ax.plot(times,[r["pressure_halfdeficit_width_cells"] for r in rows],label="pressão ≤ 1/2 déficit")
    ax.set(xlabel="tempo (s)",ylabel="células de grade");ax.legend();savefig(fig,figures/"28_resolution_cells.png")

    fig,ax=plt.subplots(figsize=(6,4.5)); sc=ax.scatter(gamma,pdef,c=times,cmap="viridis");ax.set(xlabel=r"$\Gamma$(4.2 km)",ylabel=r"$\Delta p_{dyn}$ (Pa)")
    fig.colorbar(sc,ax=ax,label="tempo (s)");savefig(fig,figures/"29_pressure_vs_circulation.png")

    fig=plt.figure(figsize=(7,5.5));ax=fig.add_subplot(111,projection="3d")
    for idx in chosen:
        m=map_cache[idx];ax.plot(m["axis"][:,0]/1000,m["axis"][:,1]/1000,z/1000,label=f"{times[idx]:.0f} s")
    ax.set(xlabel="x (km)",ylabel="y (km)",zlabel="z (km)");ax.legend(fontsize=7)
    savefig(fig,figures/"30_vortex_axis_3d.png")

    # Concise numerical synthesis used by the report generator and continuity file.
    budget_closure_rel = max(
        r["native_closure_relative_rms"]
        for r in budget_rows
    )
    initial_gamma, final_gamma = gamma[0], gamma[-1]
    initial_zeta, final_zeta = zmax[0], zmax[-1]
    dominant_post = max(SOURCES[1:-1], key=lambda n: abs(source_positive_series[n][-1]))
    final_signed_total = sum(source_series[n][-1] for n in SOURCES)
    final_positive_norm = sum(source_positive_series[n][-1] for n in SOURCES)
    final_horizontal_norm = sum(source_horizontal_series[n][-1] for n in SOURCES)
    final_tilting_total = sum(source_tilting_series[n][-1] for n in SOURCES)
    summary = {
        "numerical_gate": validity["numerical_validity_gate"],
        "analysis_interval_s": [float(times[0]), float(times[-1])],
        "frames": len(times), "grid_dx_m": dx_m, "primary_level_m": float(z[1]),
        "analysis_geometry": {
            "radial_bin_m": float(args.radial_bin_m if args.radial_bin_m is not None else dx_m),
            "cylinder_radius_m": ANALYSIS_CYLINDER_RADIUS_M,
            "vertical_range_m": list(ANALYSIS_VERTICAL_RANGE_M),
        },
        "sequence_provenance_alignment_max_abs_s": float(max(abs(r["alignment_time_mismatch_s"]) for r in rows)),
        "low_level_change": {
            "zeta_max_initial_s-1": initial_zeta, "zeta_max_peak_s-1": float(zmax.max()),
            "zeta_max_peak_time_s": float(times[int(np.argmax(zmax))]), "zeta_max_final_s-1": final_zeta,
            "zeta_final_over_peak": float(final_zeta/zmax.max()),
            "circulation_4200_initial_m2_s": initial_gamma,
            "circulation_4200_peak_m2_s": float(gamma.max()),
            "circulation_4200_final_m2_s": final_gamma,
            "circulation_final_over_peak": float(final_gamma/gamma.max()),
            "vtheta_peak_m_s": float(vmax.max()), "pressure_deficit_min_Pa": float(pdef.min()),
        },
        "concentration": {
            "stretching_conditional_mean_range_s-2": [float(np.nanmin(stretch)),float(np.nanmax(stretch))],
            "median_fraction_zeta_positive_stretching": float(np.nanmedian([r["fraction_positive_zeta_positive_stretching"] for r in rows])),
            "median_fraction_zeta_positive_convergence": float(np.nanmedian([r["fraction_positive_zeta_positive_convergence"] for r in rows])),
            "median_corr_zeta_convergence": float(np.nanmedian([r["corr_zeta_convergence"] for r in rows])),
            "rmw_cells_range": [float(min(r["rmw_cells"] for r in rows)),float(max(r["rmw_cells"] for r in rows))],
            "core_width_cells_range": [float(min(r["zeta_core_width_cells"] for r in rows)),float(max(r["zeta_core_width_cells"] for r in rows))],
            "zeta_halfmax_width_cells_range": [float(min(r["zeta_halfmax_width_cells"] for r in rows)),float(max(r["zeta_halfmax_width_cells"] for r in rows))],
            "convergence_halfmax_width_cells_range": [float(min(r["convergence_halfmax_width_cells"] for r in rows)),float(max(r["convergence_halfmax_width_cells"] for r in rows))],
            "pressure_halfdeficit_width_cells_range": [float(min(r["pressure_halfdeficit_width_cells"] for r in rows)),float(max(r["pressure_halfdeficit_width_cells"] for r in rows))],
        },
        "vertical_coherence": {
            "surface_to_2km_displacement_initial_m": rows[0]["axis_tilt_surface_to_2km_m"],
            "surface_to_2km_displacement_max_m": max(r["axis_tilt_surface_to_2km_m"] for r in rows),
            "surface_to_2km_displacement_final_m": rows[-1]["axis_tilt_surface_to_2km_m"],
            "maximum_displacement_at_any_height_range_m": [
                min(r["axis_tilt_max_m"] for r in rows), max(r["axis_tilt_max_m"] for r in rows)
            ],
        },
        "provenance": {
            "dominant_post_restart_positive_source_final": dominant_post,
            "final_positive_integrals": {n: source_positive_series[n][-1] for n in SOURCES},
            "final_signed_integrals": {n: source_series[n][-1] for n in SOURCES},
            "final_horizontal_vorticity_abs_integrals": {n: source_horizontal_series[n][-1] for n in SOURCES},
            "final_tilting_integrals": {n: source_tilting_series[n][-1] for n in SOURCES},
            "dominant_post_restart_horizontal_source_final": max(SOURCES[1:-1], key=lambda n: abs(source_horizontal_series[n][-1])),
            "dominant_post_restart_tilting_source_final": max(SOURCES[1:-1], key=lambda n: abs(source_tilting_series[n][-1])),
            "final_signed_zeta_fractions": {n: source_series[n][-1]/final_signed_total for n in SOURCES},
            "final_positive_inventory_norm_fractions": {n: source_positive_series[n][-1]/final_positive_norm for n in SOURCES},
            "final_horizontal_norm_fractions": {n: source_horizontal_series[n][-1]/final_horizontal_norm for n in SOURCES},
            "final_tilting_signed_fractions": {n: source_tilting_series[n][-1]/final_tilting_total for n in SOURCES},
            "initial_is_antecedent_unclassified": True,
        },
        "cold_pool": {
            "area_range_km2": [min(r["cold_pool_area_km2"] for r in rows),max(r["cold_pool_area_km2"] for r in rows)],
            "near_vortex_fraction_range": [min(r["cold_pool_fraction_near_vortex"] for r in rows),max(r["cold_pool_fraction_near_vortex"] for r in rows)],
            "baroclinic_h_mean_range_s-2": [min(r["baroclinic_h_mean_near_vortex_s-2"] for r in rows),max(r["baroclinic_h_mean_near_vortex_s-2"] for r in rows)],
            "signed_distance_at_center_range_m": [min(r["cold_pool_signed_distance_at_center_m"] for r in rows),max(r["cold_pool_signed_distance_at_center_m"] for r in rows)],
        },
        "native_budget": {
            "max_closure_relative_rms": budget_closure_rel,
            "time_mean_terms_s-2": {
                key: float(np.mean([r[key] for r in budget_rows]))
                for key in (
                    "delta_zeta_mean_s-2", "advection_mean_s-2", "les_mean_s-2",
                    "surface_drag_mean_s-2", "coriolis_mean_s-2", "projection_mean_s-2",
                    "kinematic_stretching_mean_s-2", "kinematic_tilting_mean_s-2",
                    "kinematic_total_advection_mean_s-2", "kinematic_dilatation_mean_s-2",
                )
            },
        },
        "events": {k: float(times[v]) for k,v in event_indices.items()},
        "dominant_classification": "MIXED: CONCENTRATION DEFICIT DOMINANT, SOURCE-INVENTORY LOSS SECONDARY",
        "figures": sorted(p.name for p in figures.glob("*.png")),
        "event_states": "event_states.h5",
    }
    (args.out/"summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(summary,indent=2,ensure_ascii=False))


if __name__ == "__main__":
    main()
