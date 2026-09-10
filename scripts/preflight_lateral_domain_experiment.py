"""Preflight the isolated 72-to-120 km lateral-domain experiment."""
from __future__ import annotations

import dataclasses
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from meteorological_flow.config import estimate_memory_gb, geometry
from run_diagnostic_sequence import build


FIELDS = ("u", "v", "w", "p", "theta", "qv", "ql", "qi", "qr", "qs", "qg", "qh")
ALLOWED_CONFIG_DIFFERENCES = {
    "sim.domain.Lx", "sim.domain.Ly", "sim.grid.nx", "sim.grid.ny",
}


def flatten(value, prefix=""):
    if isinstance(value, dict):
        out = {}
        for key, child in value.items():
            out.update(flatten(child, f"{prefix}.{key}" if prefix else key))
        return out
    return {prefix: value}


def crop(array, field, offset, n):
    if field == "u":
        return array[offset:offset+n+1, offset:offset+n, :]
    if field == "v":
        return array[offset:offset+n, offset:offset+n+1, :]
    return array[offset:offset+n, offset:offset+n, :]


def gpu_memory():
    import cupy as cp
    free, total = cp.cuda.runtime.memGetInfo()
    pool = cp.get_default_memory_pool()
    return {
        "device_used_bytes": int(total-free),
        "device_free_bytes": int(free),
        "pool_used_bytes": int(pool.used_bytes()),
        "pool_reserved_bytes": int(pool.total_bytes()),
    }


def main():
    out = ROOT / "outputs/domain_extent_120km_preflight_20260910"
    out.mkdir(parents=True, exist_ok=True)
    common = dict(device="cpu", nz=48, duration=3300.023494218,
                  rain_evaporation_factor=1.0, Lz_m=15000.)
    control, motion_c = build(nx=120, ny=120, Lx_m=72000., Ly_m=72000., **common)
    lateral, motion_l = build(nx=200, ny=200, Lx_m=120000., Ly_m=120000., **common)

    cfg_c = flatten(dataclasses.asdict(control.scfg))
    cfg_l = flatten(dataclasses.asdict(lateral.scfg))
    config_differences = sorted(key for key in cfg_c | cfg_l if cfg_c.get(key) != cfg_l.get(key))
    unexpected = sorted(set(config_differences) - ALLOWED_CONFIG_DIFFERENCES)
    assert not unexpected, unexpected
    assert set(config_differences) == ALLOWED_CONFIG_DIFFERENCES, config_differences
    assert motion_c == motion_l
    assert control.grid.dx == lateral.grid.dx == 600.
    assert control.grid.dy == lateral.grid.dy == 600.
    zc_c = control.grid.backend.to_cpu(control.grid.zc)
    zc_l = lateral.grid.backend.to_cpu(lateral.grid.zc)
    zf_c = control.grid.backend.to_cpu(control.grid.zf)
    zf_l = lateral.grid.backend.to_cpu(lateral.grid.zf)
    assert np.array_equal(zc_c, zc_l) and np.array_equal(zf_c, zf_l)

    offset = (lateral.grid.nx-control.grid.nx)//2
    initial_identity = {}
    for field in FIELDS:
        a = np.asarray(control.grid.backend.to_cpu(getattr(control.state, field)))
        b = np.asarray(lateral.grid.backend.to_cpu(getattr(lateral.state, field)))
        bcrop = crop(b, field, offset, control.grid.nx)
        initial_identity[field] = {
            "bitwise": bool(np.array_equal(a, bcrop)),
            "max_abs": float(np.max(np.abs(a-bcrop))),
        }
    assert all(v["bitwise"] for v in initial_identity.values()), initial_identity

    control_geometry = geometry(control.cfg)
    lateral_geometry = geometry(lateral.cfg)
    estimates = {
        "control_field_memory_gb": estimate_memory_gb(control.cfg),
        "lateral_field_memory_gb": estimate_memory_gb(lateral.cfg),
        "selected_sequence_gb": 6.994295972399414 * (200/240)**2,
        "selected_provenance_gb": 1.592700738 * (200/120)**2,
        "disk_free_gb": shutil.disk_usage(ROOT).free/1e9,
    }
    del control, lateral

    before = gpu_memory()
    start = time.time()
    sim, _ = build(device="gpu", nx=200, ny=200, nz=48,
                   Lx_m=120000., Ly_m=120000., Lz_m=15000., duration=1.0)
    after_initialization = gpu_memory()
    dt = min(float(sim._dt()), 1.0-sim.t)
    sim._step(dt)
    after_one_step = gpu_memory()
    finite = all(bool(sim.grid.xp.isfinite(getattr(sim.state, f)).all()) for f in FIELDS)
    assert finite
    one_step = {
        "dt_s": dt,
        "wall_clock_s": time.time()-start,
        "finite": finite,
        "w_absmax_m_s": float(sim.grid.xp.abs(sim.state.w).max()),
    }

    summary = {
        "status": "PASS",
        "experiment": "isolated lateral extension",
        "control_geometry": control_geometry,
        "lateral_geometry": lateral_geometry,
        "allowed_config_differences": config_differences,
        "vertical_grid_sha256": hashlib.sha256(zc_c.tobytes()+zf_c.tobytes()).hexdigest(),
        "storm_motion_ground_ms": list(motion_c),
        "centered_initial_state_identity": initial_identity,
        "estimates": estimates,
        "gpu_memory": {"before": before, "after_initialization": after_initialization,
                       "after_one_step": after_one_step},
        "one_step": one_step,
        "authorized_production_command": [
            sys.executable, "scripts/run_diagnostic_sequence.py",
            "--out", "outputs/domain_extent_120km_600m_20260910",
            "--device", "gpu", "--nx", "200", "--ny", "200",
            "--Lx-m", "120000", "--Ly-m", "120000", "--Lz-m", "15000",
            "--nz", "48", "--duration", "3300.023494218",
            "--interval", "30", "--capture-start", "2790",
        ],
        "interpretation": (
            "Only lateral extent changes. The centered control subdomain is bitwise-identical "
            "at initialization; subsequent differences include the intended global-domain response."
        ),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
