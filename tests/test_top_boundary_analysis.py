"""Synthetic checks for the frozen top-boundary decision pipeline."""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np

from scripts.analyze_top_boundary_pilot import (
    PRIMARY,
    budget_closure,
    geometric_gate,
    lag_correlation,
    metric_comparison,
    tracking_gate,
)
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics.top_boundary import configure_high_top
from storm_dynamics.vorticity_provenance import VorticityProvenanceTracer


def _metric_rows(scale=1.0):
    rows = []
    for i in range(18):
        row = {"time_s": 2820.0 + 30 * i, "center_x_m": 1000.0 + 100 * i,
               "center_y_m": 2000.0 + 50 * i}
        for metric, direction in PRIMARY.items():
            row[metric] = scale * (1.0 + 0.01 * i) * (1.0 if direction > 0 else 1000.0)
        row["secondary"] = scale * (2.0 + i)
        rows.append(row)
    return rows


def test_metric_thresholds_and_tracking_are_frozen():
    control = _metric_rows()
    high = _metric_rows(1.25)
    # For lower-is-favourable width, make HIGH 25% smaller instead.
    for a, b in zip(control, high):
        b["zeta_halfmax_width_m"] = 0.75 * a["zeta_halfmax_width_m"]
    _, _, primary, material, classification = metric_comparison(control, high)
    assert material == 6
    assert classification == "material favorable top response"
    assert all(r["coherent_12_of_18"] and r["material_20_percent"] for r in primary)
    assert tracking_gate(control, high)["status"] == "PASS"
    high[5]["center_x_m"] += 20000.0
    assert tracking_gate(control, high)["status"] == "FAIL"


def test_lag_correlation_recovers_planted_positive_delay():
    rng = np.random.default_rng(91)
    driver = rng.standard_normal(700)
    response = np.r_[np.zeros(37), driver[:-37]]
    rows = [{"time_s": 2790.0 + i, "driver": driver[i], "response": response[i]}
            for i in range(511)]
    _, peak = lag_correlation(rows, "driver", "response", max_lag_s=100.0)
    assert peak is not None
    assert peak["lag_s"] == 37.0
    assert peak["correlation"] > 0.99


def _synthetic_sequence(path: Path):
    nx = ny = 4; nz = 6
    zf = np.array([0., 100., 300., 600., 1000., 1500., 2000.])
    zc = 0.5 * (zf[:-1] + zf[1:])
    with h5py.File(path, "w") as f:
        f.attrs["budget_nz_with_halo"] = nz
        grid = f.create_group("grid")
        grid["xc"] = np.arange(nx) * 600. + 300.; grid["yc"] = np.arange(ny) * 600. + 300.
        grid["zc"] = zc; grid["zf"] = zf
        snaps = f.create_group("snapshots")
        previous = None
        for j, t in enumerate((2790., 2820.)):
            g = snaps.create_group(f"{j:05d}"); g.attrs["time_s"] = t
            u = np.zeros((nx + 1, ny, nz)); v = np.zeros((nx, ny + 1, nz))
            w = np.zeros((nx, ny, nz + 1)); w[:, :, 1:3] = j * 0.1
            fields = {"u": u, "v": v, "w": w, "p": np.zeros((nx, ny, nz)),
                      "theta": np.full((nx, ny, nz), 300.), "qv": np.zeros((nx, ny, nz)),
                      "ql": np.zeros((nx, ny, nz)), "qi": np.zeros((nx, ny, nz)),
                      "qr": np.zeros((nx, ny, nz)), "qs": np.zeros((nx, ny, nz)),
                      "qg": np.zeros((nx, ny, nz)), "qh": np.zeros((nx, ny, nz))}
            fields["ql"][:, :, 1] = 2e-5
            for name, value in fields.items(): g[name] = value
            if previous is not None:
                inc = g.create_group("increments").create_group("synthetic")
                for name in ("u", "v", "w"): inc[name] = fields[name] - previous[name]
            previous = fields


def test_geometry_and_budget_gates_on_synthetic_sequence(tmp_path):
    path = tmp_path / "sequence.h5"; _synthetic_sequence(path)
    rows, gate = geometric_gate(path, damping_faces=4)
    assert len(rows) == 2
    assert gate["status"] == "PASS"
    closure, closure_gate = budget_closure(path)
    assert len(closure) == 1
    assert closure_gate["status"] == "PASS"
    assert closure_gate["maximum_relative_rms"] == 0.0


def test_v4_closure_on_54_level_high_top_grid():
    cfg = build_storm_config(
        preset="storm", nx=4, ny=4, nz=48, Lx=4000., Ly=4000., Lz=15000.,
        duration=1., z_stretch=1.05, device="cpu",
    )
    cfg.sim.physics.bubble_dtheta = 5.0
    high = configure_high_top(cfg)
    sim = StormSimulation(high)
    tracer = VorticityProvenanceTracer(sim)
    sim.diagnostic_observer = tracer
    sim._step(0.01)
    metrics = tracer._closure_metrics()
    assert sim.grid.nz == 54
    assert metrics["omega_relative_rms"] <= 1e-10
    assert metrics["tilting_relative_rms"] <= 1e-10

