"""Short ON/OFF validation of the passive vorticity-provenance observer.

This is deliberately a tiny solver run for instrumentation validation.  It is
not a storm experiment and does not alter any physical parameter in the model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics.vorticity_provenance import VorticityProvenanceTracer


STATE_NAMES = (
    "u", "v", "w", "p", "theta", "qv", "ql", "qi", "qr", "qs", "qg",
    "qh", "T", "rho", "RH_w", "P_total",
)


def production_cost(nx=120, ny=120, nz=48, output_nz=20, sources=9,
                    itemsize=8, frames=32):
    cells = nx * ny * nz
    face_velocity = ((nx + 1) * ny * nz + nx * (ny + 1) * nz
                     + nx * ny * (nz + 1))
    persistent = (sources * 3 * cells + 3 * cells + face_velocity) * itemsize
    saved_scalars = 2 * sources + 2 + 2  # source omega_h, total omega_h, grad_h w
    raw_snapshot = nx * ny * output_nz * saved_scalars * itemsize
    return {
        "persistent_bytes": int(persistent),
        "raw_snapshot_bytes": int(raw_snapshot),
        "raw_32_frame_file_bytes": int(raw_snapshot * frames),
        "raw_three_case_bytes": int(raw_snapshot * frames * 3),
        "note": "working temporaries are additional and reused one source at a time; gzip size is data dependent",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "outputs" / "vorticity_provenance_validation")
    parser.add_argument("--steps", type=int, default=6)
    parser.add_argument("--dt", type=float, default=0.1)
    args = parser.parse_args()
    args.out = args.out.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    h5_path = args.out / "provenance_short.h5"
    checkpoint_path = args.out / "provenance_checkpoint_short.h5"
    for path in (h5_path, checkpoint_path):
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite {path}")

    cfg = build_storm_config(
        nx=10, ny=10, nz=8, Lx=8000.0, Ly=8000.0, Lz=4000.0,
        z_stretch=1.05, duration=args.steps * args.dt, device="cpu",
    )
    plain = StormSimulation(cfg)
    instrumented = StormSimulation(cfg)
    tracer = VorticityProvenanceTracer(
        instrumented, h5_path,
        metadata={"purpose": "short passive-instrumentation validation"},
        interval=0.2, top=2000.0, output_halo=2, storage_dtype="f8",
        record_stage_metrics=True,
    )
    instrumented.diagnostic_observer = tracer

    elapsed_plain = 0.0
    elapsed_instrumented = 0.0
    bitwise = {name: True for name in STATE_NAMES}
    max_abs = {name: 0.0 for name in STATE_NAMES}
    for _ in range(args.steps):
        started = time.perf_counter()
        plain._step(args.dt)
        elapsed_plain += time.perf_counter() - started
        plain.step += 1
        plain.t = float(plain.state.t)

        started = time.perf_counter()
        instrumented._step(args.dt)
        elapsed_instrumented += time.perf_counter() - started
        instrumented.step += 1
        instrumented.t = float(instrumented.state.t)

        for name in STATE_NAMES:
            a = plain.backend.to_cpu(getattr(plain.state, name))
            b = instrumented.backend.to_cpu(getattr(instrumented.state, name))
            bitwise[name] &= bool(np.array_equal(a, b))
            max_abs[name] = max(max_abs[name], float(np.max(np.abs(a - b))))

    final_closure = tracer._closure_metrics()
    source_max = {
        name: float(np.max(np.abs(instrumented.backend.to_cpu(
            tracer.omega_sources[tracer.source_index[name], :, :, :, :tracer.analysis_nk]
        ))))
        for name in tracer.source_names
    }
    stage_closure_max = {
        key: max(row[key] for row in tracer.history)
        for key in (
            "omega_max_abs", "omega_rms", "omega_relative_rms",
            "tilting_max_abs", "tilting_rms", "tilting_relative_rms",
        )
    }
    local_memory = tracer.memory_estimate()
    tracer.write_full_checkpoint(checkpoint_path, instrumented)
    restored = VorticityProvenanceTracer(instrumented)
    restored_closure = restored.restore_full_checkpoint(checkpoint_path, instrumented)
    checkpoint_bitwise = bool(np.array_equal(
        instrumented.backend.to_cpu(restored.omega_sources),
        instrumented.backend.to_cpu(tracer.omega_sources),
    ))
    tracer.close(instrumented)

    figure_path = args.out / "closure_residuals.png"
    event = np.arange(len(tracer.history))
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    axes[0].semilogy(event, np.maximum(
        [row["omega_max_abs"] for row in tracer.history], 1e-30
    ), label="máximo")
    axes[0].semilogy(event, np.maximum(
        [row["omega_rms"] for row in tracer.history], 1e-30
    ), label="RMS")
    axes[0].set_ylabel("resíduo de ω (s⁻¹)")
    axes[0].legend()
    axes[1].semilogy(event, np.maximum(
        [row["tilting_max_abs"] for row in tracer.history], 1e-30
    ), label="máximo")
    axes[1].semilogy(event, np.maximum(
        [row["tilting_rms"] for row in tracer.history], 1e-30
    ), label="RMS")
    axes[1].set_ylabel("resíduo de Tz (s⁻²)")
    axes[1].set_xlabel("evento de operador na validação curta")
    axes[1].legend()
    fig.suptitle("Fechamento dos traçadores passivos por estágio")
    fig.tight_layout()
    fig.savefig(figure_path, dpi=160)
    plt.close(fig)

    summary = {
        "classification": "instrumentation validation only; not a storm experiment",
        "steps": args.steps,
        "dt_s": args.dt,
        "final_time_s": instrumented.t,
        "grid": {"nx": 10, "ny": 10, "nz": 8},
        "bitwise_identity_all_prognostics": bool(all(bitwise.values())),
        "bitwise_identity": bitwise,
        "max_abs_difference": max_abs,
        "final_closure": final_closure,
        "max_closure_over_all_stages": stage_closure_max,
        "checkpoint_roundtrip": {
            "bitwise_sources": checkpoint_bitwise,
            "closure": restored_closure,
            "path": str(checkpoint_path.relative_to(ROOT)),
            "bytes": checkpoint_path.stat().st_size,
            "sha256": hashlib.sha256(checkpoint_path.read_bytes()).hexdigest(),
        },
        "source_max_abs_s-1": source_max,
        "stages_seen": sorted({row["stage"] for row in tracer.history}),
        "short_run_memory": local_memory,
        "production_120x120x48_cost_f64": production_cost(),
        "timing_s": {
            "plain": elapsed_plain,
            "instrumented": elapsed_instrumented,
            "ratio": elapsed_instrumented / max(elapsed_plain, 1e-30),
            "warning": "tiny CPU timing is dominated by overhead and is not a production benchmark",
        },
        "archive": {
            "path": str(h5_path.relative_to(ROOT)),
            "bytes": h5_path.stat().st_size,
            "sha256": hashlib.sha256(h5_path.read_bytes()).hexdigest(),
        },
        "closure_figure": {
            "path": str(figure_path.relative_to(ROOT)),
            "bytes": figure_path.stat().st_size,
            "sha256": hashlib.sha256(figure_path.read_bytes()).hexdigest(),
        },
        "source_sha256": hashlib.sha256(
            (ROOT / "src" / "storm_dynamics" / "vorticity_provenance.py").read_bytes()
        ).hexdigest(),
        "full_3900_s_run_executed": False,
        "causal_matrix_executed": False,
    }
    (args.out / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
