"""Preflight the frozen CONTROL-15/HIGH-TOP-20 causal pilot."""
from __future__ import annotations

import dataclasses
import gc
import hashlib
import json
from pathlib import Path
import shutil
import sys

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from scripts.run_diagnostic_sequence import build
from scripts.run_top_boundary_pilot import estimate_output_storage
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics.top_boundary import TopBoundaryObserver
from storm_dynamics.top_boundary import configure_high_top
from storm_dynamics.vorticity_provenance import VorticityProvenanceTracer


FIELDS = ("u", "v", "w", "p", "theta", "qv", "ql", "qi", "qr", "qs", "qg", "qh")
SCHEDULE = ROOT / "outputs/diagnostic_sequence_20260905/sequence.h5"
END = 3300.023494218
EXPECTED_SCHEDULE_SHA = "6ebc5b14ac31b0f11819340e56331913cb233bd41453d617e4296501e742d5a3"
REFERENCE_T_K = 266.84789304648075
REFERENCE_QV = 0.006483976752963724


def sha_array(a):
    a = np.ascontiguousarray(a)
    return hashlib.sha256(a.view(np.uint8)).hexdigest()


def schedule():
    with h5py.File(SCHEDULE, "r") as f:
        steps = np.asarray(f["steps"])
    steps = steps[steps[:, 0] < END - 1e-12].astype("<f8", copy=True)
    digest = hashlib.sha256(steps.tobytes()).hexdigest()
    if digest != EXPECTED_SCHEDULE_SHA:
        raise RuntimeError(f"schedule hash mismatch: {digest}")
    if abs(float(steps[-1, 0] + steps[-1, 1]) - END) > 1e-9:
        raise RuntimeError("schedule does not terminate at the registered end time")
    return steps, digest


def release(sim):
    del sim
    gc.collect()
    try:
        import cupy as cp
        cp.get_default_memory_pool().free_all_blocks()
    except Exception:
        pass


def state_hashes(sim):
    return {name: sha_array(sim.backend.to_cpu(getattr(sim.state, name))) for name in FIELDS}


def lower_state(sim, nz=48):
    result = {}
    for name in FIELDS:
        value = sim.backend.to_cpu(getattr(sim.state, name))
        result[name] = np.array(value[:, :, :nz + (name == "w")], copy=True)
    return result


def main():
    out = ROOT / "outputs/top_boundary_preflight_v3_20260911"
    out.mkdir(parents=True, exist_ok=False)
    steps, schedule_sha = schedule()
    report = {
        "status": "running",
        "schedule": {
            "source": str(SCHEDULE), "rows": len(steps), "sha256": schedule_sha,
            "first": steps[0].tolist(), "last": steps[-1].tolist(),
            "dt_min_s": float(steps[:, 1].min()), "dt_max_s": float(steps[:, 1].max()),
            "end_s": float(steps[-1, 0] + steps[-1, 1]),
        },
    }
    (out / "summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    common = dict(device="gpu", duration=END, damping_faces=4,
                  reference_T=REFERENCE_T_K, reference_qv=REFERENCE_QV)
    control, motion_control = build(**common)
    to = control.backend.to_cpu
    control_zf = np.asarray(to(control.grid.zf))
    control_zc = np.asarray(to(control.grid.zc))
    control_dzc = np.asarray(to(control.grid.dz_c))
    control_lower = lower_state(control)
    control_dt0 = float(control._dt())
    control_storage = estimate_output_storage(control, steps)
    control_cfg = dataclasses.asdict(control.scfg)
    device = control.backend.device_info()
    release(control)

    high, motion_high = build(**common, high_top_m=20000.0)
    high_zf = np.asarray(to(high.grid.zf))
    lower_equal = {name: bool(np.array_equal(control_lower[name],
                                             np.asarray(to(getattr(high.state, name)))[:, :, :48 + (name == "w")]))
                   for name in FIELDS}
    grid_checks = {
        "zf_prefix_bitwise": bool(np.array_equal(control_zf, high_zf[:49])),
        "zc_prefix_bitwise": bool(np.array_equal(control_zc, np.asarray(to(high.grid.zc))[:48])),
        "dz_c_prefix_bitwise": bool(np.array_equal(control_dzc, np.asarray(to(high.grid.dz_c))[:48])),
        "strictly_increasing": bool(np.all(np.diff(high_zf) > 0)),
        "positive_volumes": bool(np.all(np.asarray(to(high.grid.cell_vol_c)) > 0)),
        "top_m": float(high_zf[-1]), "nz": high.grid.nz,
        "control_reference_dz_m": 312.5, "high_reference_dz_m": float(high.grid.dz),
        "new_faces_m": high_zf[49:].tolist(),
        "damping_faces_m": high_zf[-4:].tolist(),
    }
    high_cfg = dataclasses.asdict(high.scfg)
    own_dt0 = float(high._dt())
    high_storage = estimate_output_storage(high, steps)
    imposed_dt0 = float(steps[0, 1])
    if imposed_dt0 > own_dt0 * (1.0 + 1e-12):
        raise RuntimeError("first registered time step violates HIGH-TOP CFL")
    high._step(imposed_dt0)
    finite_plain = {name: bool(np.isfinite(high.backend.to_cpu(getattr(high.state, name))).all())
                    for name in FIELDS}
    hashes_plain = state_hashes(high)
    release(high)

    observed, motion_observed = build(**common, high_top_m=20000.0)
    observer = TopBoundaryObserver(
        observed, out / "observer_one_step.h5", metadata={"purpose": "actual-grid passivity"},
        flush_every=4, full_field_start_s=0.0,
    )
    observed.top_boundary_observer = observer
    observed._step(imposed_dt0)
    hashes_observed = state_hashes(observed)
    observer.close(completed=True)
    passivity = {name: hashes_plain[name] == hashes_observed[name] for name in FIELDS}
    with h5py.File(out / "observer_one_step.h5", "r") as f:
        observer_checks = {
            "calls": int(f.attrs["calls"]),
            "contexts": list(f["context"].asstr()[:]),
            "energy_finite": bool(np.isfinite(f["energy_removed_J"][:]).all()),
            "full_fields_calls": len(f["full_fields/step"]),
        }
    release(observed)

    small_cfg = build_storm_config(
        preset="storm", nx=4, ny=4, nz=48, Lx=4000.0, Ly=4000.0, Lz=15000.0,
        duration=1.0, z_stretch=1.05, device="cpu",
    )
    small_cfg.sim.physics.bubble_dtheta = 5.0
    small = StormSimulation(configure_high_top(small_cfg))
    tracer = VorticityProvenanceTracer(small)
    small.diagnostic_observer = tracer
    small._step(0.01)
    v4_closure = tracer._closure_metrics()

    allowed = {
        "sim.domain.Lz", "sim.grid.nz", "sim.grid.z_faces_m",
        "sim.grid.vertical_reference_dz_m",
    }

    def flatten(value, prefix=""):
        if isinstance(value, dict):
            result = {}
            for key, item in value.items():
                result.update(flatten(item, f"{prefix}.{key}" if prefix else key))
            return result
        return {prefix: value}

    fc, fh = flatten(control_cfg), flatten(high_cfg)
    differences = sorted(k for k in set(fc) | set(fh) if fc.get(k) != fh.get(k))
    report.update({
        "status": "PASS",
        "device": device,
        "storm_motion_control_ms": motion_control,
        "storm_motion_high_top_ms": motion_high,
        "storm_motion_observed_ms": motion_observed,
        "storm_motion_bitwise": sha_array(np.asarray(motion_control)) == sha_array(np.asarray(motion_high)),
        "shared_reference_values": {"T_ref_K": REFERENCE_T_K, "qv_ref": REFERENCE_QV},
        "storage_estimates": {
            "control15": control_storage,
            "hightop20": high_storage,
            "reserve_free_gib_per_arm": 3.0,
            "required_free_before_control_gib": control_storage["estimated_gib"] + 3.0,
            "required_free_before_hightop_gib": high_storage["estimated_gib"] + 3.0,
        },
        "config_differences": differences,
        "config_differences_allowed": sorted(allowed),
        "only_registered_config_differences": set(differences) == allowed,
        "grid_checks": grid_checks,
        "lower_initial_fields_bitwise": lower_equal,
        "first_step_cfl": {
            "imposed_dt_s": imposed_dt0, "control_limit_s": control_dt0,
            "high_top_limit_s": own_dt0,
            "ratio_imposed_to_high_limit": imposed_dt0 / own_dt0,
        },
        "one_step_fields_finite": finite_plain,
        "actual_grid_observer_passivity_bitwise": passivity,
        "observer_checks": observer_checks,
        "v4_54_level_closure": v4_closure,
        "disk_free_gib_after": shutil.disk_usage(out).free / 2**30,
    })
    required = [
        report["only_registered_config_differences"], report["storm_motion_bitwise"],
        all(grid_checks[k] for k in ("zf_prefix_bitwise", "zc_prefix_bitwise",
                                     "dz_c_prefix_bitwise", "strictly_increasing", "positive_volumes")),
        all(lower_equal.values()), all(finite_plain.values()), all(passivity.values()),
        observer_checks["calls"] == 4,
        observer_checks["contexts"] == list(TopBoundaryObserver.contexts),
        v4_closure["omega_relative_rms"] <= 1e-10,
        v4_closure["tilting_relative_rms"] <= 1e-10,
    ]
    if not all(required):
        report["status"] = "FAIL"
    (out / "summary.json").write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    print(json.dumps(report, indent=2, default=str))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
