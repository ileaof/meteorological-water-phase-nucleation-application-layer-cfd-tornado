"""Mathematical closure and passivity tests for provenance diagnostics."""
import numpy as np
import pytest

from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics.momentum import momentum_advection_tendency
from storm_dynamics.vorticity_provenance import (
    SOURCE_NAMES,
    VorticityProvenanceTracer,
    discrete_homogeneous_vorticity_tendency,
    frozen_muscl_advection_tendency,
    frozen_momentum_advection_tendency,
    homogeneous_vorticity_tendency,
    velocity_kinematics,
)


def _small_sim(device="cpu"):
    cfg = build_storm_config(
        nx=10, ny=10, nz=8, Lx=8000.0, Ly=8000.0, Lz=4000.0,
        z_stretch=1.05, duration=1.0, device=device,
    )
    return StormSimulation(cfg)


def test_homogeneous_provenance_operator_is_linear():
    sim = _small_sim()
    tracer = VorticityProvenanceTracer(sim)
    centered, gradient, divergence, _ = velocity_kinematics(
        tracer.previous_velocity, sim.grid
    )
    rng = np.random.default_rng(4)
    a = rng.normal(size=tracer.total_omega.shape)
    b = rng.normal(size=tracer.total_omega.shape)
    lhs = homogeneous_vorticity_tendency(a + b, centered, gradient, divergence, sim.grid)
    rhs = (
        homogeneous_vorticity_tendency(a, centered, gradient, divergence, sim.grid)
        + homogeneous_vorticity_tendency(b, centered, gradient, divergence, sim.grid)
    )
    assert np.allclose(lhs, rhs, rtol=2e-15, atol=2e-15)


def test_frozen_muscl_operator_is_linear_and_reconstructs_total():
    sim = _small_sim()
    grid = sim.grid
    rng = np.random.default_rng(41)
    a = rng.normal(size=(grid.nx, grid.ny, grid.nz))
    b = rng.normal(size=a.shape)
    total = a + b
    velocity = (sim.state.u, sim.state.v, sim.state.w)
    whole = frozen_muscl_advection_tendency(total, total, velocity, grid)
    pieces = (
        frozen_muscl_advection_tendency(a, total, velocity, grid)
        + frozen_muscl_advection_tendency(b, total, velocity, grid)
    )
    assert np.allclose(whole, pieces, rtol=3e-14, atol=3e-14)


def test_frozen_muscl_uniform_advection_is_bounded():
    sim = _small_sim()
    grid = sim.grid
    x = np.arange(grid.nx, dtype=float)[:, None, None]
    q = np.sin(2.0 * np.pi * x / grid.nx) * np.ones(
        (1, grid.ny, grid.nz), dtype=float
    )
    u = np.full(grid.u_shape, 10.0)
    v = np.zeros(grid.v_shape)
    w = np.zeros(grid.w_shape)
    dt = 0.45 * grid.dx / 10.0
    initial_bound = float(np.max(np.abs(q)))
    for _ in range(250):
        q += dt * frozen_muscl_advection_tendency(q, q, (u, v, w), grid)
    assert np.all(np.isfinite(q))
    assert float(np.max(np.abs(q))) <= initial_bound + 1e-12


def test_discrete_vorticity_operator_reconstructs_split_field():
    sim = _small_sim()
    tracer = VorticityProvenanceTracer(sim)
    _, gradient, _, _ = velocity_kinematics(tracer.previous_velocity, sim.grid)
    rng = np.random.default_rng(42)
    a = rng.normal(size=tracer.total_omega.shape)
    b = tracer.total_omega - a
    whole = discrete_homogeneous_vorticity_tendency(
        tracer.total_omega,
        tracer.total_omega,
        tracer.previous_velocity,
        gradient,
        sim.grid,
    )
    pieces = (
        discrete_homogeneous_vorticity_tendency(
            a, tracer.total_omega, tracer.previous_velocity, gradient, sim.grid
        )
        + discrete_homogeneous_vorticity_tendency(
            b, tracer.total_omega, tracer.previous_velocity, gradient, sim.grid
        )
    )
    assert np.allclose(whole, pieces, rtol=5e-13, atol=5e-13)


def test_frozen_staggered_momentum_matches_native_and_split_sum():
    sim = _small_sim()
    grid = sim.grid
    total = (sim.state.u.copy(), sim.state.v.copy(), sim.state.w.copy())
    native = momentum_advection_tendency(sim.state, grid, order=2)
    frozen_total = frozen_momentum_advection_tendency(total, total, grid, order=2)
    for expected, actual in zip(native, frozen_total):
        assert np.array_equal(expected, actual)

    rng = np.random.default_rng(43)
    a = tuple(rng.normal(size=value.shape) for value in total)
    b = tuple(value - part for value, part in zip(total, a))
    ta = frozen_momentum_advection_tendency(a, total, grid, order=2)
    tb = frozen_momentum_advection_tendency(b, total, grid, order=2)
    for expected, left, right in zip(native, ta, tb):
        assert np.allclose(expected, left + right, rtol=2e-13, atol=2e-13)


def test_buoyancy_injects_horizontal_but_no_vertical_vorticity():
    sim = _small_sim()
    tracer = VorticityProvenanceTracer(sim)
    sim.diagnostic_observer = tracer
    sim._predictor(0.1)
    buoy = tracer.omega_sources[tracer.source_index["buoyancy"]]
    assert np.max(np.abs(buoy[:2])) > 0.0
    assert np.array_equal(buoy[2], np.zeros_like(buoy[2]))


def test_optional_full_vorticity_output_preserves_all_components(tmp_path):
    h5py = pytest.importorskip("h5py")
    sim = _small_sim()
    path = tmp_path / "full_vorticity.h5"
    tracer = VorticityProvenanceTracer(
        sim, path, interval=1.0, top=1000.0, output_halo=2,
        write_full_vorticity=True,
    )
    tracer.close(sim)
    with h5py.File(path, "r") as handle:
        assert bool(handle.attrs["write_full_vorticity"])
        sample = handle["snapshots/00000"]
        assert sample["omega_total"].shape[0] == 3
        assert set(sample["omega_source"]) == set(SOURCE_NAMES)
        reconstructed = sum(sample["omega_source"][name][:] for name in SOURCE_NAMES)
        assert np.max(np.abs(reconstructed - sample["omega_total"][:])) < 2e-15
        assert np.array_equal(sample["omega_total"][:2], sample["omega_total_h"][:])


def test_close_time_snapshot_uses_completed_driver_step(tmp_path):
    h5py = pytest.importorskip("h5py")
    sim = _small_sim()
    path = tmp_path / "close_step.h5"
    tracer = VorticityProvenanceTracer(sim, path, interval=10.0)
    sim.diagnostic_observer = tracer
    sim._step(0.1)
    sim.step += 1
    sim.t = float(sim.state.t)
    tracer.close(sim)
    with h5py.File(path, "r") as handle:
        assert int(handle["snapshots/00001"].attrs["step"]) == sim.step


@pytest.mark.parametrize("device", ["cpu", "gpu"])
def test_provenance_is_bitwise_passive_and_closes(tmp_path, device):
    h5py = pytest.importorskip("h5py")
    if device == "gpu":
        cp = pytest.importorskip("cupy")
        try:
            cp.zeros(1)
        except Exception as exc:
            pytest.skip(str(exc))
    plain = _small_sim(device)
    observed = _small_sim(device)
    path = tmp_path / "provenance.h5"
    tracer = VorticityProvenanceTracer(
        observed, path, metadata={"purpose": "unit test"}, interval=0.2,
        top=1000.0, output_halo=2,
    )
    observed.diagnostic_observer = tracer

    for _ in range(4):
        for sim in (plain, observed):
            sim._step(0.1)
            sim.step += 1
            sim.t = float(sim.state.t)
        for name in (
            "u", "v", "w", "p", "theta", "qv", "ql", "qi",
            "qr", "qs", "qg", "qh", "T", "rho", "P_total",
        ):
            a = plain.backend.to_cpu(getattr(plain.state, name))
            b = observed.backend.to_cpu(getattr(observed.state, name))
            assert np.array_equal(a, b), name
        if hasattr(plain.state, "p_dyn") or hasattr(observed.state, "p_dyn"):
            assert np.array_equal(
                getattr(plain.state, "p_dyn"), getattr(observed.state, "p_dyn")
            )

    final_metrics = tracer._closure_metrics()
    assert final_metrics["omega_max_abs"] < 2e-15
    assert final_metrics["tilting_max_abs"] < 2e-17
    assert any(row["stage"] == "advection" for row in tracer.history)
    checkpoint = tmp_path / "provenance_checkpoint.h5"
    tracer.write_full_checkpoint(checkpoint, observed)
    restored = VorticityProvenanceTracer(observed)
    restored_metrics = restored.restore_full_checkpoint(checkpoint, observed)
    assert restored_metrics["omega_max_abs"] < 2e-15
    assert np.array_equal(
        observed.backend.to_cpu(restored.omega_sources),
        observed.backend.to_cpu(tracer.omega_sources),
    )
    forked_path = tmp_path / "forked_provenance.h5"
    forked = VorticityProvenanceTracer(
        observed, forked_path, interval=1.0, source_checkpoint=checkpoint
    )
    assert np.array_equal(
        observed.backend.to_cpu(forked.omega_sources),
        observed.backend.to_cpu(tracer.omega_sources),
    )
    forked.close(observed)
    tracer.close(observed)

    with h5py.File(path, "r") as handle:
        assert handle.attrs["schema"] == "storm-vorticity-provenance-v4"
        assert bool(handle.attrs["passive"])
        assert bool(handle.attrs["internal_full_column"])
        assert handle.attrs["status"] == "complete"
        assert len(handle["snapshots"]) == 3
        sample = handle["snapshots/00001"]
        assert set(sample["omega_source_h"]) == set(SOURCE_NAMES)
        reconstructed = sum(sample["omega_source_h"][name][:] for name in SOURCE_NAMES)
        assert np.max(np.abs(reconstructed - sample["omega_total_h"][:])) < 2e-15
        assert "closure_history" in handle
