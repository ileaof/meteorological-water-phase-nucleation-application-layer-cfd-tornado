"""Run the validated v4 provenance instrument over the mature storm interval.

The instrumented state is advanced beside an uninstrumented control using the
archived native time-step schedule.  Every prognostic field is compared after
every step; any loss of bitwise neutrality aborts the run.  The restart-time
field is labelled ``initial`` (antecedent), so source attribution applies only
to increments accumulated after the selected restart.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from run_causal_evaporation_branch import PROGNOSTIC, restore
from run_diagnostic_sequence import build
from storm_dynamics.vorticity_provenance import VorticityProvenanceTracer


def _sha256(path: Path, chunk_size: int = 16 * 1024**2) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def _source_hashes() -> dict[str, str]:
    paths = (
        ROOT / "src/storm_dynamics/vorticity_provenance.py",
        ROOT / "src/storm_dynamics/core.py",
        ROOT / "src/storm_dynamics/momentum.py",
        ROOT / "scripts/run_vorticity_provenance_long.py",
        ROOT / "scripts/run_diagnostic_sequence.py",
        ROOT / "scripts/run_causal_evaporation_branch.py",
    )
    return {str(path.relative_to(ROOT)): _sha256(path) for path in paths}


def _scalar(xp, value) -> float:
    return float(value.item() if hasattr(value, "item") else value)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sequence", type=Path,
        default=ROOT / "outputs/diagnostic_sequence_20260905/sequence.h5",
    )
    parser.add_argument("--snapshot", type=int, default=93)
    parser.add_argument("--target-snapshot", type=int, default=110)
    parser.add_argument("--interval", type=float, default=30.0)
    parser.add_argument("--device", choices=("cpu", "gpu"), default="gpu")
    parser.add_argument(
        "--out", type=Path,
        default=ROOT / "outputs/vorticity_provenance_long_v4_2790_3300",
    )
    args = parser.parse_args()
    args.sequence = args.sequence.resolve()
    args.out = args.out.resolve()
    if args.out.exists() and any(args.out.iterdir()):
        raise FileExistsError(f"Refusing nonempty output directory: {args.out}")
    args.out.mkdir(parents=True, exist_ok=True)
    archive = args.out / "provenance_long_v4.h5"

    with h5py.File(args.sequence, "r") as handle:
        start_group = handle[f"snapshots/{args.snapshot:05d}"]
        target_group = handle[f"snapshots/{args.target_snapshot:05d}"]
        start_time = float(start_group.attrs["time_s"])
        start_step = int(start_group.attrs["step"])
        target_time = float(target_group.attrs["time_s"])
        target_step = int(target_group.attrs["step"])
        all_schedule = handle["steps"][:]
        nx = len(handle["grid/xc"])
        ny = len(handle["grid/yc"])
        nz = len(handle["grid/zc"])
        if nx != ny:
            raise ValueError("build requires a square horizontal grid")
    schedule = all_schedule[
        (all_schedule[:, 2] >= start_step) & (all_schedule[:, 2] < target_step)
    ]
    if len(schedule) != target_step - start_step:
        raise ValueError("archived native time-step schedule is incomplete")
    if abs(start_time + float(schedule[:, 1].sum()) - target_time) > 2e-8:
        raise ValueError("archived schedule does not land on target snapshot")

    plain, storm_motion = build(args.device, nx, nz, target_time, 1.0)
    observed, _ = build(args.device, nx, nz, target_time, 1.0)
    restore(plain, args.sequence, args.snapshot)
    restore(observed, args.sequence, args.snapshot)
    xp = observed.grid.xp

    metadata = {
        "purpose": "authorized long v4 tornadogenesis provenance run",
        "source_sequence": str(args.sequence),
        "restart_snapshot": args.snapshot,
        "target_snapshot": args.target_snapshot,
        "start_time_s": start_time,
        "target_time_s": target_time,
        "initial_label_semantics": (
            "all velocity/vorticity antecedent to restart; post-restart source "
            "labels do not retrospectively classify it"
        ),
        "provenance_formulation": "storm-vorticity-provenance-v4-staggered-muscl",
        "storm_motion_ground_ms": storm_motion,
        "rain_evaporation_factor": 1.0,
        "git_head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "git_status": subprocess.check_output(
            ["git", "status", "--short"], cwd=ROOT, text=True
        ),
        "source_sha256": _source_hashes(),
        "command": sys.argv,
    }
    tracer = VorticityProvenanceTracer(
        observed, archive, metadata=metadata, interval=args.interval,
        top=2000.0, output_halo=3, storage_dtype="f8",
        record_stage_metrics=True, write_full_vorticity=True,
    )
    observed.diagnostic_observer = tracer

    neutrality = {
        name: {
            "bitwise": True, "max_abs": 0.0, "max_rms": 0.0,
            "max_relative_rms": 0.0,
        }
        for name in PROGNOSTIC
    }
    progress = []
    started = time.perf_counter()
    completed = False
    try:
        for number, (t0, dt, archived_step) in enumerate(schedule, start=1):
            history_start = len(tracer.history)
            for sim in (plain, observed):
                if abs(float(sim.t) - float(t0)) > 2e-8:
                    raise RuntimeError(
                        f"archived schedule mismatch {sim.t} versus {t0}"
                    )
                sim._step(float(dt))
                sim.step += 1
                sim.t = float(sim.state.t)

            mismatch = False
            for name in PROGNOSTIC:
                a = getattr(plain.state, name)
                b = getattr(observed.state, name)
                equal = bool(xp.array_equal(a, b))
                neutrality[name]["bitwise"] &= equal
                if not equal:
                    delta = a - b
                    max_abs = _scalar(xp, xp.max(xp.abs(delta)))
                    rms = _scalar(xp, xp.sqrt(xp.mean(delta * delta)))
                    scale = _scalar(xp, xp.sqrt(xp.mean(a * a)))
                    row = neutrality[name]
                    row["max_abs"] = max(row["max_abs"], max_abs)
                    row["max_rms"] = max(row["max_rms"], rms)
                    row["max_relative_rms"] = max(
                        row["max_relative_rms"], rms / max(scale, 1e-30)
                    )
                    mismatch = True
            if mismatch:
                raise RuntimeError(
                    f"prognostic neutrality failed at t={observed.t:.9f} s"
                )

            # Original audit tolerances, applied to every recorded stage.
            for metrics in tracer.history[history_start:]:
                limits = {"omega_relative_rms": 1e-10, "omega_max_abs": 1e-10,
                          "tilting_relative_rms": 1e-10, "tilting_max_abs": 1e-12}
                if not all(np.isfinite(value) for key, value in metrics.items()
                           if key not in ("stage",)):
                    raise RuntimeError("Nonfinite provenance metric")
                for key, limit in limits.items():
                    if metrics[key] > limit:
                        raise RuntimeError(f"Provenance gate failed: {key}={metrics[key]}")

            if number == 1 or number % 50 == 0 or number == len(schedule):
                closure = tracer._closure_metrics()
                row = {
                    "time_s": float(observed.t), "step": int(observed.step),
                    "schedule_index": number, "neutral_bitwise": True,
                    **closure,
                }
                progress.append(row)
                message = (
                    f"t={observed.t:.3f}s step={observed.step} "
                    f"({number}/{len(schedule)}) kappa_omega="
                    f"{closure['omega_condition_index']:.6g} kappa_T="
                    f"{closure['tilting_condition_index']:.6g} closure="
                    f"{closure['omega_relative_rms']:.3e} wall="
                    f"{time.perf_counter()-started:.1f}s"
                )
                print(message, flush=True)
                with (args.out / "progress.log").open("a", encoding="utf-8") as f:
                    f.write(message + "\n")
        completed = True
    finally:
        final_closure = tracer._closure_metrics()
        history = list(tracer.history)
        tracer.close(observed, completed=completed)

        closure_keys = (
            "omega_max_abs", "omega_rms", "omega_relative_rms",
            "tilting_max_abs", "tilting_rms", "tilting_relative_rms",
            "omega_condition_index", "tilting_condition_index",
            "omega_advection_remainder_fraction",
            "omega_advection_remainder_source_fraction",
            "tilting_advection_remainder_fraction",
            "tilting_advection_remainder_source_fraction",
        )
        extrema = {
            key: {
                "min": min(float(row[key]) for row in history),
                "max": max(float(row[key]) for row in history),
                "final": float(final_closure[key]),
            }
            for key in closure_keys
        }
        summary = {
            "status": "complete" if completed else "interrupted",
            "numerical_validity_gate": "PASS" if completed else "FAIL",
            "start_time_s": start_time,
            "end_time_s": float(observed.t),
            "target_time_s": target_time,
            "native_steps": len(schedule),
            "restart_snapshot": args.snapshot,
            "target_snapshot": args.target_snapshot,
            "device": args.device,
            "grid": {"nx": nx, "ny": ny, "nz": nz},
            "neutrality": neutrality,
            "closure_conditioning": extrema,
            "monitor_samples": progress,
            "wall_clock_s": time.perf_counter() - started,
            "archive": {
                "path": str(archive.relative_to(ROOT)),
                "bytes": archive.stat().st_size if archive.exists() else None,
                "sha256": _sha256(archive) if archive.exists() else None,
            },
            "metadata": metadata,
        }
        (args.out / "summary.json").write_text(
            json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
