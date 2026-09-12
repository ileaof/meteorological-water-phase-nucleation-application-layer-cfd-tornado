"""Run one frozen arm of the CONTROL-15/HIGH-TOP-20 causal pilot."""
from __future__ import annotations

import argparse
import csv
import dataclasses
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from scripts.run_diagnostic_sequence import build
from storm_dynamics.diagnostic_capture import DiagnosticCapture
from storm_dynamics.top_boundary import TopBoundaryObserver


SCHEDULE_PATH = ROOT / "outputs/diagnostic_sequence_20260905/sequence.h5"
SCHEDULE_SHA256 = "6ebc5b14ac31b0f11819340e56331913cb233bd41453d617e4296501e742d5a3"
END_S = 3300.023494218
CAPTURE_START_S = 2790.0
REFERENCE_T_K = 266.84789304648075
REFERENCE_QV = 0.006483976752963724
FIELDS = ("u", "v", "w", "p", "theta", "qv", "ql", "qi", "qr", "qs", "qg", "qh")


def load_schedule():
    with h5py.File(SCHEDULE_PATH, "r") as f:
        steps = np.asarray(f["steps"])
    steps = steps[steps[:, 0] < END_S - 1e-12].astype("<f8", copy=True)
    digest = hashlib.sha256(steps.tobytes()).hexdigest()
    if digest != SCHEDULE_SHA256:
        raise RuntimeError(f"registered schedule hash mismatch: {digest}")
    if abs(float(steps[-1, 0] + steps[-1, 1]) - END_S) > 1e-9:
        raise RuntimeError("registered schedule does not reach the frozen end time")
    return steps


def estimate_output_storage(sim, steps):
    """Return a conservative, dimensioned upper bound for one arm's output."""
    g = sim.grid
    nx, ny, nz = int(g.nx), int(g.ny), int(g.nz)
    itemsize = max(np.dtype(getattr(sim.state.w, "dtype", np.float64)).itemsize, 8)
    z = np.asarray(sim.backend.to_cpu(g.zc))
    nk = min(nz, int(np.searchsorted(z, 2000.0, side="right")) + 3)

    # Reproduce DiagnosticCapture's time triggering conservatively, including
    # its initial frame and a possible final partial-interval frame.
    first_capture = int(np.flatnonzero(steps[:, 0] + steps[:, 1] >= CAPTURE_START_S)[0])
    capture_t0 = float(steps[first_capture, 0] + steps[first_capture, 1])
    diagnostic_frames = int(np.ceil((END_S - capture_t0) / 30.0)) + 1
    diagnostic_intervals = max(diagnostic_frames - 1, 0)

    center = nx * ny * nz
    native_full = (nx + 1) * ny * nz + nx * (ny + 1) * nz + nx * ny * (nz + 1)
    center_low = nx * ny * nk
    native_low = ((nx + 1) * ny * nk + nx * (ny + 1) * nk
                  + nx * ny * (nk + 1))
    # 13 centered state arrays plus u/v/w; 12 native momentum stages; and
    # three vector plus two scalar kinematic integrals (11 centered arrays).
    diagnostic_state = diagnostic_frames * itemsize * (13 * center + native_full)
    diagnostic_increments = diagnostic_intervals * itemsize * (12 * native_low + 11 * center_low)

    full_steps = int(np.count_nonzero(steps[:, 0] >= CAPTURE_START_S))
    full_calls = 4 * full_steps
    top_spatial = full_calls * itemsize * (2 * nx * ny * 4)
    # Eight full-column profiles plus damping/curl/flux vectors and scalars.
    top_profiles = len(steps) * 4 * 8 * (10 * nz + 80)
    audit_tables = len(steps) * 6 * 8

    raw = diagnostic_state + diagnostic_increments + top_spatial + top_profiles + audit_tables
    estimated = int(np.ceil(raw * 1.25 + 256 * 2**20))
    return {
        "method": "uncompressed arrays plus 25% HDF5 overhead and 256 MiB fixed allowance",
        "itemsize_bytes": itemsize,
        "diagnostic_nz": nk,
        "diagnostic_frames_upper_bound": diagnostic_frames,
        "top_full_field_steps": full_steps,
        "top_full_field_calls": full_calls,
        "components_bytes": {
            "diagnostic_states": diagnostic_state,
            "diagnostic_increments_and_kinematics": diagnostic_increments,
            "top_spatial_fields": top_spatial,
            "top_profiles": top_profiles,
            "audit_tables": audit_tables,
        },
        "raw_bytes": raw,
        "estimated_bytes": estimated,
        "estimated_gib": estimated / 2**30,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", choices=("control15", "hightop20"), required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--device", choices=("cpu", "gpu"), default="gpu")
    parser.add_argument("--reserve-free-gib", type=float, default=3.0)
    args = parser.parse_args()
    steps = load_schedule()
    common = dict(
        device=args.device, duration=END_S, damping_faces=4,
        reference_T=REFERENCE_T_K, reference_qv=REFERENCE_QV,
    )
    sim, motion = build(**common, high_top_m=(20000.0 if args.arm == "hightop20" else None))
    storage = estimate_output_storage(sim, steps)
    free_before = shutil.disk_usage(args.out.parent).free
    required = storage["estimated_bytes"] + int(args.reserve_free_gib * 2**30)
    if free_before < required:
        raise RuntimeError(
            f"pilot requires {required / 2**30:.2f} GiB free "
            f"({storage['estimated_gib']:.2f} GiB estimate + "
            f"{args.reserve_free_gib:g} GiB reserve), found {free_before / 2**30:.2f} GiB"
        )
    args.out.mkdir(parents=True, exist_ok=False)
    hashes = {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for directory in ("src", "scripts") for p in (ROOT / directory).rglob("*.py")
    }
    metadata = dict(
        status="running", arm=args.arm, registered_experiment="docs/TOP_BOUNDARY_CAUSAL_EXPERIMENT.md",
        config=dataclasses.asdict(sim.scfg), source_sha256=hashes,
        git_head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        git_diff=subprocess.check_output(["git", "diff"], cwd=ROOT, text=True),
        command=sys.argv, backend=sim.backend.device_info(), storm_motion_ground_ms=motion,
        shared_reference_T_K=REFERENCE_T_K, shared_reference_qv=REFERENCE_QV,
        schedule_source=str(SCHEDULE_PATH), schedule_sha256=SCHEDULE_SHA256,
        schedule_rows=len(steps), capture_start_target_s=CAPTURE_START_S,
        final_time_target_s=END_S, cfl_tolerance_fraction=1e-12,
        top_observer_full_field_start_s=CAPTURE_START_S,
        storage_estimate=storage, reserve_free_gib=args.reserve_free_gib,
        disk_free_gib_before=free_before / 2**30,
    )
    (args.out / "metadata.json").write_text(json.dumps(metadata, indent=2, default=str), encoding="utf-8")
    observer = TopBoundaryObserver(
        sim, args.out / "top_boundary_calls.h5", metadata=metadata,
        flush_every=64, full_field_start_s=CAPTURE_START_S,
    )
    sim.top_boundary_observer = observer
    capture = None
    completed = False
    cfl_rows = []
    start = time.time()
    last_report = -1e30
    failure = None
    try:
        for row_index, (t_start, imposed_dt, archived_step) in enumerate(steps):
            if abs(sim.t - float(t_start)) > 2e-9:
                raise RuntimeError(f"schedule time mismatch at row {row_index}: {sim.t} vs {t_start}")
            own_limit = float(sim._dt())
            ratio = float(imposed_dt) / own_limit
            cfl_rows.append((row_index, sim.step, sim.t, imposed_dt, own_limit, ratio))
            if ratio > 1.0 + 1e-12:
                raise RuntimeError(
                    f"CFL_SCHEDULE_VIOLATION row={row_index} t={sim.t:.15g} "
                    f"imposed={imposed_dt:.15g} own={own_limit:.15g} ratio={ratio:.15g}"
                )
            sim._step(float(imposed_dt))
            sim.step += 1
            sim.t = float(sim.state.t)
            if capture is None and sim.t >= CAPTURE_START_S:
                capture = DiagnosticCapture(sim, args.out / "sequence.h5", metadata, interval=30.0)
                sim.diagnostic_observer = capture
            if sim.t - last_report >= 30.0 or row_index == len(steps) - 1:
                finite = {name: bool(sim.grid.xp.isfinite(getattr(sim.state, name)).all())
                          for name in FIELDS}
                if not all(finite.values()):
                    raise RuntimeError(f"nonfinite prognostic field: {finite}")
                free_gib = shutil.disk_usage(args.out).free / 2**30
                if free_gib < 3.0:
                    raise RuntimeError("capture stopped: less than 3 GiB free")
                message = (
                    f"arm={args.arm} t={sim.t:.6f}s step={sim.step}/{len(steps)} "
                    f"frames={capture.count if capture else 0} ratio={ratio:.9f} "
                    f"free={free_gib:.2f}GiB wall={time.time()-start:.1f}s"
                )
                print(message, flush=True)
                with (args.out / "progress.log").open("a", encoding="utf-8") as stream:
                    stream.write(message + "\n")
                last_report = sim.t
        completed = True
    except Exception as exc:
        failure = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        if capture is not None:
            capture.close(sim, completed=completed)
        observer.close(completed=completed)
        with (args.out / "cfl_audit.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(("row", "step", "time_s", "imposed_dt_s", "own_limit_s", "ratio"))
            writer.writerows(cfl_rows)
        metadata.update(
            status="complete" if completed else "interrupted", failure=failure,
            final_time_s=sim.t, steps=sim.step, frames=capture.count if capture else 0,
            cfl_max_ratio=max((r[-1] for r in cfl_rows), default=None),
            cfl_violations=sum(r[-1] > 1.0 + 1e-12 for r in cfl_rows),
            wall_clock_s=time.time() - start, disk_free_gib=shutil.disk_usage(args.out).free / 2**30,
        )
        (args.out / "metadata.json").write_text(json.dumps(metadata, indent=2, default=str), encoding="utf-8")


if __name__ == "__main__":
    main()
