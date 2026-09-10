"""Finite-window provenance validation from an archived mature-storm state."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sequence", type=Path,
        default=ROOT / "outputs" / "diagnostic_sequence_20260905" / "sequence.h5",
    )
    parser.add_argument("--snapshot", type=int, default=93)
    parser.add_argument("--steps", type=int, default=2)
    parser.add_argument("--device", choices=("cpu", "gpu"), default="gpu")
    parser.add_argument(
        "--out", type=Path,
        default=ROOT / "outputs" / "vorticity_provenance_validation_mature",
    )
    args = parser.parse_args()
    args.sequence = args.sequence.resolve()
    args.out = args.out.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    provenance_path = args.out / "provenance_mature_short.h5"
    if provenance_path.exists():
        raise FileExistsError(f"Refusing to overwrite {provenance_path}")

    with h5py.File(args.sequence, "r") as handle:
        group = handle[f"snapshots/{args.snapshot:05d}"]
        start_time = float(group.attrs["time_s"])
        start_step = int(group.attrs["step"])
        schedule = handle["steps"][:]
    schedule = schedule[schedule[:, 2] >= start_step][:args.steps]
    if len(schedule) != args.steps:
        raise ValueError("archived dt schedule is too short")
    end_time = float(start_time + schedule[:, 1].sum())

    plain, _ = build(args.device, 120, 48, end_time, 1.0)
    observed, _ = build(args.device, 120, 48, end_time, 1.0)
    restore(plain, args.sequence, args.snapshot)
    restore(observed, args.sequence, args.snapshot)
    cp = None
    if args.device == "gpu":
        import cupy as cp

    def gpu_memory():
        if cp is None:
            return None
        free, total = cp.cuda.runtime.memGetInfo()
        return {
            "free_bytes": int(free),
            "total_bytes": int(total),
            "device_used_bytes": int(total - free),
            "pool_used_bytes": int(cp.get_default_memory_pool().used_bytes()),
            "pool_reserved_bytes": int(cp.get_default_memory_pool().total_bytes()),
        }

    memory_before_tracer = gpu_memory()
    tracer = VorticityProvenanceTracer(
        observed, provenance_path,
        metadata={
            "purpose": "finite-window mature-state instrumentation validation",
            "source_sequence": str(args.sequence),
            "snapshot": args.snapshot,
            "initial_label_semantics": "antecedent vorticity at validation restart",
        },
        interval=30.0, top=2000.0, output_halo=3, storage_dtype="f8",
        record_stage_metrics=True,
    )
    observed.diagnostic_observer = tracer
    memory_after_tracer = gpu_memory()
    peak_device_used = 0 if memory_after_tracer is None else memory_after_tracer["device_used_bytes"]
    peak_pool_used = 0 if memory_after_tracer is None else memory_after_tracer["pool_used_bytes"]
    peak_pool_reserved = 0 if memory_after_tracer is None else memory_after_tracer["pool_reserved_bytes"]

    timings = {"plain": 0.0, "instrumented": 0.0}
    bitwise = {name: True for name in PROGNOSTIC}
    max_abs = {name: 0.0 for name in PROGNOSTIC}
    for t0, dt, _ in schedule:
        for name, sim in (("plain", plain), ("instrumented", observed)):
            if abs(float(sim.t) - float(t0)) > 2e-8:
                raise RuntimeError(f"archived schedule mismatch: {sim.t} versus {t0}")
            started = time.perf_counter()
            sim._step(float(dt))
            timings[name] += time.perf_counter() - started
            sim.step += 1
            sim.t = float(sim.state.t)
            if name == "instrumented" and cp is not None:
                cp.cuda.Stream.null.synchronize()
                current_memory = gpu_memory()
                peak_device_used = max(peak_device_used, current_memory["device_used_bytes"])
                peak_pool_used = max(peak_pool_used, current_memory["pool_used_bytes"])
                peak_pool_reserved = max(peak_pool_reserved, current_memory["pool_reserved_bytes"])
        for name in PROGNOSTIC:
            a = plain.backend.to_cpu(getattr(plain.state, name))
            b = observed.backend.to_cpu(getattr(observed.state, name))
            bitwise[name] &= bool(np.array_equal(a, b))
            max_abs[name] = max(max_abs[name], float(np.max(np.abs(a - b))))

    final_closure = tracer._closure_metrics()
    stage_closure = {
        key: max(row[key] for row in tracer.history)
        for key in (
            "omega_max_abs", "omega_rms", "omega_relative_rms",
            "tilting_max_abs", "tilting_rms", "tilting_relative_rms",
        )
    }
    memory = tracer.memory_estimate()
    tracer.close(observed)

    summary = {
        "classification": "finite-window mature-state instrumentation validation only",
        "source_sequence": str(args.sequence),
        "snapshot": args.snapshot,
        "start_time_s": start_time,
        "end_time_s": observed.t,
        "native_steps": args.steps,
        "native_dt_s": schedule[:, 1].tolist(),
        "device": args.device,
        "grid": {"nx": 120, "ny": 120, "nz": 48},
        "bitwise_identity_all_prognostics": bool(all(bitwise.values())),
        "bitwise_identity": bitwise,
        "max_abs_difference": max_abs,
        "final_closure": final_closure,
        "max_closure_over_all_stages": stage_closure,
        "memory": memory,
        "gpu_memory": {
            "before_tracer_two_simulations": memory_before_tracer,
            "after_tracer_two_simulations": memory_after_tracer,
            "peak_device_used_bytes": peak_device_used,
            "peak_pool_used_bytes": peak_pool_used,
            "peak_pool_reserved_bytes": peak_pool_reserved,
            "measurement_note": "two prognostic simulations coexist for bitwise OFF/ON comparison; production runs one branch at a time",
        },
        "timing_s": {
            **timings,
            "ratio": timings["instrumented"] / max(timings["plain"], 1e-30),
            "scope": "requested native-step window; includes device synchronisation from closure metrics",
        },
        "archive": {
            "path": str(provenance_path.relative_to(ROOT)),
            "bytes": provenance_path.stat().st_size,
            "sha256": hashlib.sha256(provenance_path.read_bytes()).hexdigest(),
        },
        "physical_interpretation_allowed": False,
        "full_run_executed": False,
    }
    (args.out / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
