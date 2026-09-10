#!/usr/bin/env python3
"""Diagnose the stopped vorticity-provenance experiment without physical attribution.

This is an offline, read-only analysis of the shared and partial CONTROL archives.
It quantifies closure, cancellation/conditioning and the exact automatic-stop event.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


LIMITS = {
    "omega_relative_rms": 1.0e-10,
    "omega_max_abs": 1.0e-10,
    "tilting_relative_rms": 1.0e-10,
    "tilting_max_abs": 1.0e-12,
}


def rms(value: np.ndarray) -> float:
    value = np.asarray(value, dtype=np.float64)
    return float(np.sqrt(np.mean(value * value)))


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    av = np.asarray(a, dtype=np.float64).ravel()
    bv = np.asarray(b, dtype=np.float64).ravel()
    scale = np.linalg.norm(av) * np.linalg.norm(bv)
    return float(np.dot(av, bv) / scale) if scale else math.nan


def snapshot_rows(path: Path) -> tuple[list[dict], list[str]]:
    rows: list[dict] = []
    with h5py.File(path, "r") as handle:
        names = tuple(json.loads(handle.attrs["source_names"]))
        nk = int(handle.attrs["analysis_nz"])
        for key in sorted(handle["snapshots"], key=int):
            group = handle[f"snapshots/{key}"]
            total = np.asarray(group["omega_total_h"][:, :, :, :nk], dtype=np.float64)
            grad = np.asarray(group["grad_w_h"][:, :, :, :nk], dtype=np.float64)
            source = {
                name: np.asarray(
                    group[f"omega_source_h/{name}"][:, :, :, :nk], dtype=np.float64
                )
                for name in names
            }
            total_t = np.sum(total * grad, axis=0)
            source_t = {name: np.sum(value * grad, axis=0) for name, value in source.items()}
            sum_source = np.sum(np.stack(list(source.values())), axis=0)
            sum_t = np.sum(np.stack(list(source_t.values())), axis=0)
            omega_total_rms = rms(total)
            t_total_rms = rms(total_t)
            omega_norm_sum = sum(rms(value) for value in source.values())
            t_norm_sum = sum(rms(value) for value in source_t.values())
            remainder = source["advection_remainder"]
            other = sum_source - remainder
            remainder_t = source_t["advection_remainder"]
            other_t = sum_t - remainder_t
            row = {
                "time_s": float(group.attrs["time_s"]),
                "step": int(group.attrs["step"]),
                "omega_total_rms": omega_total_rms,
                "omega_source_norm_sum": omega_norm_sum,
                "omega_condition_index": omega_norm_sum / max(omega_total_rms, 1e-300),
                "tilting_total_rms": t_total_rms,
                "tilting_source_norm_sum": t_norm_sum,
                "tilting_condition_index": t_norm_sum / max(t_total_rms, 1e-300),
                "omega_remainder_other_cosine": cosine(remainder, other),
                "tilting_remainder_other_cosine": cosine(remainder_t, other_t),
                "omega_remainder_other_norm_ratio": rms(remainder) / max(rms(other), 1e-300),
                "tilting_remainder_other_norm_ratio": rms(remainder_t) / max(rms(other_t), 1e-300),
                "omega_reconstruction_rms": rms(total - sum_source),
                "tilting_reconstruction_rms": rms(total_t - sum_t),
            }
            for name in names:
                row[f"omega_{name}_rms"] = rms(source[name])
                row[f"tilting_{name}_rms"] = rms(source_t[name])
            rows.append(row)
    return rows, list(names)


def closure_rows(path: Path) -> np.ndarray:
    with h5py.File(path, "r") as handle:
        rows = handle["closure_history"][...]
    return rows[rows["stage"] == b"transport_microphysics_bcs"]


def is_bad(row) -> bool:
    return any(float(row[key]) > limit for key, limit in LIMITS.items())


def first_persistent_failure(rows: np.ndarray, length: int = 3) -> tuple[int, int] | None:
    run = 0
    for index, row in enumerate(rows):
        run = run + 1 if is_bad(row) else 0
        if run >= length:
            return index - length + 1, index
    return None


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_closure(rows: np.ndarray, out: Path) -> None:
    t = rows["time_s"]
    fig, axes = plt.subplots(2, 1, figsize=(9.0, 7.2), sharex=True)
    specs = (
        ("omega_relative_rms", "RMS relativo de vorticidade", LIMITS["omega_relative_rms"]),
        ("tilting_relative_rms", "RMS relativo de tilting", LIMITS["tilting_relative_rms"]),
    )
    for ax, (field, label, threshold) in zip(axes, specs):
        ax.semilogy(t, rows[field], color="#1f4e79", lw=1.2, label=label)
        ax.axhline(threshold, color="#b22222", ls="--", lw=1.2, label="limite")
        ax.grid(True, which="both", alpha=0.25)
        ax.set_ylabel(label)
        ax.legend(loc="upper left")
    axes[-1].set_xlabel("tempo (s)")
    fig.suptitle("Fechamento online no ramo CONTROL")
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_conditioning(rows: list[dict], out: Path) -> None:
    rows = [row for row in rows if row["time_s"] >= 2370.0]
    t = np.array([r["time_s"] for r in rows])
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    ax.semilogy(t, [r["omega_condition_index"] for r in rows], "o-", label=r"$\omega_h$")
    ax.semilogy(t, [r["tilting_condition_index"] for r in rows], "s-", label=r"$T_z$")
    ax.set_xlabel("tempo (s)")
    ax.set_ylabel("soma dos RMS das parcelas / RMS do total")
    ax.set_title("Índice de cancelamento da decomposição (instantâneos 2370–2981 s)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_source_growth(rows: list[dict], names: list[str], out: Path) -> None:
    rows = [row for row in rows if row["time_s"] >= 2370.0]
    t = np.array([r["time_s"] for r in rows])
    fig, axes = plt.subplots(2, 1, figsize=(10.0, 8.0), sharex=True)
    axes[0].semilogy(t, [r["omega_total_rms"] for r in rows], color="black", lw=2.2, label="total")
    axes[1].semilogy(t, [r["tilting_total_rms"] for r in rows], color="black", lw=2.2, label="total")
    for name in names:
        omega_values = [r[f"omega_{name}_rms"] for r in rows]
        tilting_values = [r[f"tilting_{name}_rms"] for r in rows]
        if max(omega_values) > 0.0:
            axes[0].semilogy(t, omega_values, label=name)
        if max(tilting_values) > 0.0:
            axes[1].semilogy(t, tilting_values, label=name)
    axes[0].set_ylabel(r"RMS de $\omega_h$ (s$^{-1}$)")
    axes[1].set_ylabel(r"RMS de $T_z$ (s$^{-2}$)")
    axes[1].set_xlabel("tempo (s)")
    axes[0].set_title("Crescimento das parcelas e permanência do campo total (2370–2981 s)")
    for ax in axes:
        ax.grid(True, which="both", alpha=0.2)
    axes[0].legend(ncol=5, fontsize=7)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_centered_euler(out: Path) -> None:
    phase = np.linspace(0.0, np.pi, 500)
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for courant in (0.01, 0.10, 0.30):
        amplification = np.sqrt(1.0 + courant**2 * np.sin(phase) ** 2)
        ax.plot(phase / np.pi, amplification - 1.0, label=f"C={courant:.2f}")
    ax.axhline(0.0, color="black", lw=0.8)
    ax.set_xlabel(r"número de onda normalizado $k\Delta x/\pi$")
    ax.set_ylabel(r"$|G|-1$")
    ax.set_title("Euler explícito + diferença centrada: amplificação por passo")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shared", type=Path, default=Path("outputs/vorticity_provenance_shared/shared_provenance_sparse.h5"))
    parser.add_argument("--control", type=Path, default=Path("outputs/vorticity_provenance_control_v2/provenance.h5"))
    parser.add_argument("--out", type=Path, default=Path("outputs/vorticity_provenance_closure_failure"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    shared, names = snapshot_rows(args.shared)
    control, control_names = snapshot_rows(args.control)
    if names != control_names:
        raise ValueError("source labels differ between archives")
    combined = shared + [row for row in control if row["time_s"] > shared[-1]["time_s"] + 1e-9]
    combined.sort(key=lambda row: row["time_s"])
    closure = closure_rows(args.control)
    failure = first_persistent_failure(closure)
    if failure is None:
        raise RuntimeError("no persistent three-step closure failure found")
    first, last = failure

    write_csv(args.out / "conditioning_timeseries.csv", combined)
    plot_closure(closure, args.out / "closure_thresholds.png")
    plot_conditioning(combined, args.out / "conditioning_growth.png")
    plot_source_growth(combined, names, args.out / "source_rms_growth.png")
    plot_centered_euler(args.out / "centered_euler_amplification.png")

    start = min(combined, key=lambda row: abs(row["time_s"] - 2370.2115521238206))
    peak = min(combined, key=lambda row: abs(row["time_s"] - 2790.2534902770794))
    final = combined[-1]
    summary = {
        "classification": "numerical-instrument failure diagnosis; physical attribution prohibited",
        "archives": {"shared": str(args.shared.resolve()), "control": str(args.control.resolve())},
        "limits": LIMITS,
        "persistent_failure": {
            "required_consecutive_steps": 3,
            "first_index": int(first),
            "last_index": int(last),
            "rows": [
                {name: (value.decode() if isinstance(value, bytes) else value.item()) for name, value in zip(closure.dtype.names, closure[i])}
                for i in range(first, last + 1)
            ],
        },
        "conditioning": {
            "at_2370_s": {key: start[key] for key in start if "condition_index" in key},
            "at_2790_s": {key: peak[key] for key in peak if "condition_index" in key},
            "at_stop": {key: final[key] for key in final if "condition_index" in key},
            "stop_cancellation": {
                key: final[key]
                for key in final
                if "remainder_other" in key or key.endswith("reconstruction_rms")
            },
        },
        "growth_2370_to_stop": {
            key: final[key] / max(start[key], 1e-300)
            for key in final
            if (key.startswith("omega_") or key.startswith("tilting_")) and key.endswith("_rms")
        },
        "mathematical_diagnosis": {
            "tracer_transport": "forward Euler in time with centered spatial derivatives",
            "one_dimensional_uniform_advection_amplification": "|G|=sqrt(1+C^2 sin^2(k dx)) > 1 for C != 0 and nontrivial modes",
            "solver_transport": "native flux-form MUSCL update",
            "conclusion": "the labelled homogeneous propagator is not a stable, discretely consistent decomposition of the native momentum update",
        },
    }
    (args.out / "failure_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
