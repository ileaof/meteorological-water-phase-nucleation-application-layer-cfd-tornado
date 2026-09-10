"""Run the approved shared-prefix/fork vorticity-provenance experiment.

Modes are separate so every long product is restartable and reviewable:

* ``shared``: exact archived dt schedule from 0 to snapshot 79 (~2370 s),
  passive provenance active from t=0, then full prognostic/provenance checkpoints.
* ``branch``: restore both checkpoints and integrate one approved evaporation
  branch to snapshot 110 (~3300 s), with conventional and provenance capture.
"""
from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from run_diagnostic_sequence import build
from storm_dynamics.diagnostic_capture import DiagnosticCapture
from storm_dynamics.vorticity_provenance import (
    DiagnosticObserverChain,
    VorticityProvenanceTracer,
)


PROGNOSTIC = ("u", "v", "w", "p", "theta", "qv", "ql", "qi", "qr", "qs", "qg", "qh")
FACTORS = {"WEAK-EVAP": 0.95, "CONTROL": 1.0, "STRONG-EVAP": 1.05}
MIN_FREE_BYTES = 20 * 1024 ** 3


def sha256_file(path, chunk=16 * 1024 ** 2):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while data := stream.read(chunk):
            digest.update(data)
    return digest.hexdigest()


def source_hashes():
    paths = [
        ROOT / "src" / "storm_dynamics" / "vorticity_provenance.py",
        ROOT / "src" / "storm_dynamics" / "core.py",
        ROOT / "src" / "storm_dynamics" / "momentum.py",
        ROOT / "src" / "storm_dynamics" / "turbulence.py",
        ROOT / "src" / "storm_dynamics" / "surface_drag.py",
        ROOT / "src" / "storm_dynamics" / "coriolis.py",
        ROOT / "src" / "meteorological_flow" / "buoyancy.py",
        ROOT / "src" / "meteorological_flow" / "pressure_solver.py",
        ROOT / "scripts" / "run_vorticity_provenance_experiment.py",
    ]
    return {str(path.relative_to(ROOT)): sha256_file(path) for path in paths}


def base_metadata(sim, command, source_sequence):
    config = dataclasses.asdict(sim.scfg)
    return {
        "command": command,
        "python": sys.version,
        "numpy": np.__version__,
        "platform": platform.platform(),
        "backend": sim.backend.device_info(),
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "git_status": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True),
        "config": config,
        "config_sha256": hashlib.sha256(
            json.dumps(config, sort_keys=True, separators=(",", ":"), default=str).encode()
        ).hexdigest(),
        "source_sha256": source_hashes(),
        "source_sequence": str(Path(source_sequence).resolve()),
        "provenance_formulation": "storm-vorticity-provenance-v4-staggered-muscl",
        "closure_limits": {
            "omega_relative_rms": 1e-10,
            "omega_max_abs_s-1": 1e-10,
            "tilting_relative_rms": 1e-10,
            "tilting_max_abs_s-2": 1e-12,
        },
        "purpose": "diagnostic provenance; not tornado tuning",
    }


def save_state_checkpoint(path, sim, next_dt):
    path = Path(path)
    with h5py.File(path, "x") as handle:
        handle.attrs["schema"] = "storm-prognostic-checkpoint-v1"
        handle.attrs["time_s"] = float(sim.t)
        handle.attrs["state_time_s"] = float(sim.state.t)
        handle.attrs["step"] = int(sim.step)
        handle.attrs["next_archived_dt_s"] = float(next_dt)
        handle.attrs["random_seed"] = int(sim.cfg.random_seed)
        state = handle.create_group("state")
        for name in PROGNOSTIC:
            state.create_dataset(
                name, data=sim.grid.backend.to_cpu(getattr(sim.state, name)),
                compression="gzip", compression_opts=1, shuffle=True,
            )
        if hasattr(sim.state, "p_dyn"):
            state.create_dataset(
                "p_dyn", data=sim.grid.backend.to_cpu(sim.state.p_dyn),
                compression="gzip", compression_opts=1, shuffle=True,
            )
        precip = state.create_group("surface_precip")
        for name, value in sim.state.surface_precip.items():
            precip.create_dataset(name, data=sim.grid.backend.to_cpu(value), compression="gzip")
        integrator = handle.create_group("integrator")
        integrator.attrs["last_residual"] = float(sim._last_res)
        integrator.attrs["last_iterations"] = int(sim._last_iters)
        if getattr(sim, "_Km", None) is not None:
            integrator.create_dataset("Km", data=sim.grid.backend.to_cpu(sim._Km), compression="gzip")
        if getattr(sim, "_tke", None) is not None:
            integrator.create_dataset("tke", data=sim.grid.backend.to_cpu(sim._tke), compression="gzip")
        grid = handle.create_group("grid")
        for name in ("xc", "yc", "zc", "xf", "yf", "zf"):
            grid.create_dataset(name, data=sim.grid.backend.to_cpu(getattr(sim.grid, name)))
        base = handle.create_group("base")
        for name in ("zc", "theta0", "qv0", "p0", "T0", "rho0", "u0", "v0"):
            base.create_dataset(name, data=np.asarray(getattr(sim.base, name)))


def restore_state_checkpoint(path, sim):
    with h5py.File(path, "r") as handle:
        if handle.attrs.get("schema") != "storm-prognostic-checkpoint-v1":
            raise ValueError("unsupported prognostic checkpoint")
        for name in ("xc", "yc", "zc", "xf", "yf", "zf"):
            if not np.array_equal(sim.grid.backend.to_cpu(getattr(sim.grid, name)), handle[f"grid/{name}"][:]):
                raise ValueError(f"checkpoint grid mismatch: {name}")
        for name in ("zc", "theta0", "qv0", "p0", "T0", "rho0", "u0", "v0"):
            if not np.array_equal(np.asarray(getattr(sim.base, name)), handle[f"base/{name}"][:]):
                raise ValueError(f"checkpoint base mismatch: {name}")
        xp = sim.grid.xp
        for name in PROGNOSTIC:
            setattr(sim.state, name, xp.asarray(handle[f"state/{name}"][:]))
        if "p_dyn" in handle["state"]:
            sim.state.p_dyn = xp.asarray(handle["state/p_dyn"][:])
        sim.state.surface_precip = {
            name: xp.asarray(value[:]) for name, value in handle["state/surface_precip"].items()
        }
        sim.t = float(handle.attrs["time_s"])
        sim.state.t = float(handle.attrs["state_time_s"])
        sim.step = int(handle.attrs["step"])
        sim._last_res = float(handle["integrator"].attrs["last_residual"])
        sim._last_iters = int(handle["integrator"].attrs["last_iterations"])
        sim._Km = xp.asarray(handle["integrator/Km"][:]) if "Km" in handle["integrator"] else None
        if "tke" in handle["integrator"]:
            sim._tke = xp.asarray(handle["integrator/tke"][:])
        sim.state.diagnose(sim.cfg)
        return float(handle.attrs["next_archived_dt_s"])


def closure_failure(metrics):
    return (
        metrics["omega_relative_rms"] > 1e-10
        or metrics["omega_max_abs"] > 1e-10
        or metrics["tilting_relative_rms"] > 1e-10
        or metrics["tilting_max_abs"] > 1e-12
    )


def schedule_and_target(sequence, start_snapshot, end_snapshot):
    with h5py.File(sequence, "r") as handle:
        start = handle[f"snapshots/{start_snapshot:05d}"]
        end = handle[f"snapshots/{end_snapshot:05d}"]
        start_time = float(start.attrs["time_s"])
        start_step = int(start.attrs["step"])
        end_time = float(end.attrs["time_s"])
        end_step = int(end.attrs["step"])
        steps = handle["steps"][:]
    schedule = steps[(steps[:, 2] >= start_step) & (steps[:, 2] < end_step)]
    return start_time, start_step, end_time, end_step, schedule


def run_schedule(sim, tracer, schedule, out, capture=None):
    last_report = float(sim.t) - 31.0
    failures = 0
    started = time.time()
    for t0, dt, step_index in schedule:
        if abs(float(sim.t) - float(t0)) > 2e-8 or int(sim.step) != int(step_index):
            raise RuntimeError(f"archived schedule mismatch at {sim.t}/{sim.step}: {t0}/{step_index}")
        sim._step(float(dt))
        sim.step += 1
        sim.t = float(sim.state.t)
        metrics = tracer._closure_metrics()
        failures = failures + 1 if closure_failure(metrics) else 0
        if failures >= 3:
            raise RuntimeError(f"persistent provenance closure failure: {metrics}")
        if sim.t - last_report >= 30.0 or sim.t >= float(schedule[-1, 0] + schedule[-1, 1]) - 1e-8:
            if not bool(sim.grid.xp.isfinite(sim.state.w).all()):
                raise RuntimeError("nonfinite velocity")
            free = shutil.disk_usage(out).free
            if free < MIN_FREE_BYTES:
                raise RuntimeError(f"free disk below 20 GiB: {free}")
            line = (
                f"t={sim.t:.6f} step={sim.step} wall={time.time()-started:.1f}s "
                f"Romega={metrics['omega_max_abs']:.3e} RT={metrics['tilting_max_abs']:.3e} "
                f"Komega={metrics['omega_condition_index']:.3e} "
                f"KT={metrics['tilting_condition_index']:.3e} "
                f"FremT={metrics['tilting_advection_remainder_fraction']:.3e} "
                f"prov_frames={tracer.count} capture_frames={capture.count if capture else 0} "
                f"free_GiB={free/1024**3:.2f}"
            )
            print(line, flush=True)
            with (Path(out) / "progress.log").open("a", encoding="utf-8") as stream:
                stream.write(line + "\n")
            last_report = sim.t
    return time.time() - started


def validate_against_snapshot(sim, sequence, snapshot):
    result = {"bitwise": True, "max_abs": 0.0, "fields": {}}
    with h5py.File(sequence, "r") as handle:
        group = handle[f"snapshots/{snapshot:05d}"]
        for name in PROGNOSTIC:
            current = sim.grid.backend.to_cpu(getattr(sim.state, name))
            reference = group[name][:]
            same = bool(np.array_equal(current, reference))
            error = float(np.max(np.abs(current - reference)))
            result["fields"][name] = {"bitwise": same, "max_abs": error}
            result["bitwise"] &= same
            result["max_abs"] = max(result["max_abs"], error)
    return result


def run_shared(args):
    args.out.mkdir(parents=True, exist_ok=False)
    with h5py.File(args.sequence, "r") as handle:
        target = handle[f"snapshots/{args.shared_snapshot:05d}"]
        target_time = float(target.attrs["time_s"])
        target_step = int(target.attrs["step"])
        steps = handle["steps"][:]
    schedule = steps[steps[:, 2] < target_step]
    if len(schedule) != target_step:
        raise RuntimeError("shared archived schedule is incomplete")
    sim, motion = build(args.device, 120, 48, target_time, 1.0)
    metadata = base_metadata(sim, sys.argv, args.sequence)
    metadata.update({
        "mode": "shared", "start_time_s": 0.0, "target_time_s": target_time,
        "target_snapshot": args.shared_snapshot, "storm_motion_ground_ms": motion,
        "rain_evaporation_factor": 1.0,
    })
    tracer = VorticityProvenanceTracer(
        sim, args.out / "shared_provenance_sparse.h5", metadata=metadata,
        interval=target_time + 1.0, top=2000.0, output_halo=3,
        storage_dtype="f8", record_stage_metrics=True,
    )
    sim.diagnostic_observer = tracer
    wall = run_schedule(sim, tracer, schedule, args.out)
    validation = validate_against_snapshot(sim, args.sequence, args.shared_snapshot)
    if not validation["bitwise"]:
        tracer.close(sim, completed=False)
        raise RuntimeError(f"shared state does not reproduce archived checkpoint: {validation}")
    next_rows = steps[steps[:, 2] >= target_step]
    next_dt = float(next_rows[0, 1])
    state_path = args.out / "state_checkpoint.h5"
    provenance_path = args.out / "provenance_checkpoint.h5"
    save_state_checkpoint(state_path, sim, next_dt)
    tracer.write_full_checkpoint(provenance_path, sim)
    final_closure = tracer._closure_metrics()
    tracer.close(sim)
    metadata.update({
        "status": "complete", "wall_clock_s": wall, "final_time_s": sim.t,
        "final_step": sim.step, "next_archived_dt_s": next_dt,
        "historical_validation": validation, "final_closure": final_closure,
        "state_checkpoint": {"path": str(state_path), "bytes": state_path.stat().st_size, "sha256": sha256_file(state_path)},
        "provenance_checkpoint": {"path": str(provenance_path), "bytes": provenance_path.stat().st_size, "sha256": sha256_file(provenance_path)},
        "historical_sequence": {"bytes": args.sequence.stat().st_size, "sha256": sha256_file(args.sequence)},
    })
    (args.out / "metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    print(json.dumps({k: metadata[k] for k in ("status", "wall_clock_s", "final_time_s", "final_step", "historical_validation", "final_closure")}, indent=2), flush=True)


def run_branch(args):
    if args.case not in FACTORS:
        raise ValueError(f"case must be one of {tuple(FACTORS)}")
    args.out.mkdir(parents=True, exist_ok=False)
    factor = FACTORS[args.case]
    start_time, start_step, end_time, end_step, schedule = schedule_and_target(
        args.sequence, args.shared_snapshot, args.end_snapshot
    )
    sim, motion = build(args.device, 120, 48, end_time, factor)
    next_dt = restore_state_checkpoint(args.shared / "state_checkpoint.h5", sim)
    if abs(sim.t - start_time) > 1e-8 or sim.step != start_step:
        raise RuntimeError("shared checkpoint does not match branch schedule")
    if abs(next_dt - float(schedule[0, 1])) > 1e-15:
        raise RuntimeError("stored next dt does not match archived schedule")
    metadata = base_metadata(sim, sys.argv, args.sequence)
    metadata.update({
        "mode": "branch", "case": args.case, "rain_evaporation_factor": factor,
        "start_time_s": start_time, "start_step": start_step,
        "end_time_s": end_time, "end_step": end_step,
        "shared_checkpoint": str(args.shared.resolve()),
        "timestep_semantics": "exact archived CONTROL dt schedule imposed on all branches",
        "storm_motion_ground_ms": motion,
        "only_deliberate_branch_difference": "rain_evaporation_factor",
    })
    capture = DiagnosticCapture(
        sim, args.out / "sequence.h5", metadata, interval=30.0,
        top=2000.0, full_column=False,
    )
    provenance = VorticityProvenanceTracer(
        sim, args.out / "provenance.h5", metadata=metadata, interval=30.0,
        top=2000.0, output_halo=3, storage_dtype="f8",
        record_stage_metrics=True,
        source_checkpoint=args.shared / "provenance_checkpoint.h5",
    )
    absolute_next = float(np.ceil(start_time / 30.0) * 30.0)
    capture.next_output = absolute_next
    provenance.next_output = absolute_next
    sim.diagnostic_observer = DiagnosticObserverChain(capture, provenance)
    completed = False
    try:
        wall = run_schedule(sim, provenance, schedule, args.out, capture=capture)
        completed = abs(sim.t - end_time) < 1e-7 and sim.step == end_step
    finally:
        sim.diagnostic_observer.close(sim, completed=completed)
    if not completed:
        raise RuntimeError("branch did not reach approved endpoint")
    final_closure = provenance._closure_metrics()
    metadata.update({
        "status": "complete", "wall_clock_s": wall, "final_time_s": sim.t,
        "final_step": sim.step, "frames": capture.count,
        "provenance_frames": provenance.count, "final_closure": final_closure,
        "sequence": {"bytes": (args.out / "sequence.h5").stat().st_size, "sha256": sha256_file(args.out / "sequence.h5")},
        "provenance": {"bytes": (args.out / "provenance.h5").stat().st_size, "sha256": sha256_file(args.out / "provenance.h5")},
    })
    (args.out / "metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    print(json.dumps({k: metadata[k] for k in ("case", "status", "wall_clock_s", "final_time_s", "final_step", "frames", "provenance_frames", "final_closure")}, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("shared", "branch"))
    parser.add_argument("--sequence", type=Path, default=ROOT / "outputs" / "diagnostic_sequence_20260905" / "sequence.h5")
    parser.add_argument("--shared", type=Path, default=ROOT / "outputs" / "vorticity_provenance_shared")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--case", choices=tuple(FACTORS))
    parser.add_argument("--shared-snapshot", type=int, default=79)
    parser.add_argument("--end-snapshot", type=int, default=110)
    parser.add_argument("--device", choices=("cpu", "gpu"), default="gpu")
    args = parser.parse_args()
    args.sequence = args.sequence.resolve()
    args.shared = args.shared.resolve()
    args.out = args.out.resolve()
    if shutil.disk_usage(args.out.parent).free < MIN_FREE_BYTES:
        raise RuntimeError("less than 20 GiB free before run")
    if args.mode == "shared":
        run_shared(args)
    else:
        run_branch(args)


if __name__ == "__main__":
    main()
