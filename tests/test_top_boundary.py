"""Explicit high-top geometry and passive top-damping instrumentation."""
from __future__ import annotations

import numpy as np
import pytest

from meteorological_flow.boundary_conditions import apply_velocity_bcs
from meteorological_flow.grid import Grid
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics.top_boundary import (
    TopBoundaryObserver,
    configure_high_top,
    continued_top_faces,
)


def _control_config(device="cpu", nx=8, ny=8):
    return build_storm_config(
        preset="storm", nx=nx, ny=ny, nz=48,
        Lx=8000.0, Ly=8000.0, Lz=15000.0,
        duration=1.0, z_stretch=1.05, device=device,
    )


def test_legacy_grid_path_is_bitwise_unchanged():
    expected_w = 1.05 ** np.arange(48)
    expected_zf = 15000.0 * np.concatenate(([0.0], np.cumsum(expected_w))) / expected_w.sum()
    legacy = Grid(8, 8, 48, 8000.0, 8000.0, 15000.0, z_stretch=1.05)
    explicit_none = Grid(
        8, 8, 48, 8000.0, 8000.0, 15000.0,
        z_stretch=1.05, z_faces_m=None, vertical_reference_dz_m=None,
    )
    assert np.array_equal(legacy.zf, expected_zf)
    assert np.array_equal(legacy.zf, explicit_none.zf)
    assert np.array_equal(legacy.zc, explicit_none.zc)
    assert np.array_equal(legacy.dz_c, explicit_none.dz_c)
    assert legacy.dz == explicit_none.dz == 312.5


@pytest.mark.parametrize("faces", [
    [0.0, 1.0],
    [0.0, 1.0, 0.5, 3.0],
    [0.0, 1.0, np.nan, 3.0],
    [0.0, 1.0, 2.0, 2.5],
])
def test_explicit_faces_reject_invalid_geometry(faces):
    with pytest.raises(ValueError):
        Grid(3, 3, 3, 3.0, 3.0, 3.0, z_faces_m=faces)


def test_high_top_preserves_all_control_faces_and_lower_cells_bitwise():
    cfg = _control_config()
    high = configure_high_top(cfg)
    control_grid = Grid(8, 8, 48, 8000.0, 8000.0, 15000.0, z_stretch=1.05)
    high_grid = Grid(
        8, 8, high.sim.grid.nz, 8000.0, 8000.0, high.sim.domain.Lz,
        z_stretch=1.05, z_faces_m=high.sim.grid.z_faces_m,
        vertical_reference_dz_m=high.sim.grid.vertical_reference_dz_m,
    )
    assert high.sim.grid.nz == 54
    assert high.sim.boundaries.damping_faces == 4
    assert np.array_equal(high_grid.zf[:49], control_grid.zf)
    assert np.array_equal(high_grid.zc[:48], control_grid.zc)
    assert np.array_equal(high_grid.dz_c[:48], control_grid.dz_c)
    assert high_grid.zf[-1] == 20000.0
    assert high_grid.dz == control_grid.dz == 312.5
    assert np.all(np.diff(high_grid.zf) > 0)
    assert np.all(high_grid.cell_vol_c > 0)
    # Four retained sponge faces are all strictly above the old top except the
    # no-op top coefficient; no face at or below 15 km is directly multiplied.
    affected = high_grid.zf[-high.sim.boundaries.damping_faces:]
    assert np.all(affected > 15000.0)


def test_continued_faces_reach_exact_top_with_positive_final_cell():
    control = Grid(4, 4, 48, 4000.0, 4000.0, 15000.0, z_stretch=1.05)
    extended = continued_top_faces(control.zf, 20000.0)
    assert np.array_equal(extended[:49], control.zf)
    assert extended[-1] == 20000.0
    assert np.all(np.diff(extended) > 0)


def test_damping_observer_is_bitwise_passive_and_records_four_calls(tmp_path):
    h5py = pytest.importorskip("h5py")
    cfg = _control_config()
    plain, observed = StormSimulation(cfg), StormSimulation(cfg)
    observer = TopBoundaryObserver(
        observed, tmp_path / "top.h5", flush_every=2, full_field_start_s=0.0
    )
    observed.top_boundary_observer = observer
    plain._step(0.1)
    observed._step(0.1)
    names = ("u", "v", "w", "p", "theta", "qv", "ql", "qi", "qr", "qs", "qg", "qh")
    for name in names:
        assert np.array_equal(getattr(plain.state, name), getattr(observed.state, name)), name
    observer.close()
    with h5py.File(tmp_path / "top.h5") as f:
        assert f.attrs["status"] == "complete"
        assert f.attrs["passive"]
        assert int(f.attrs["calls"]) == 4
        assert list(f["context"].asstr()[:]) == list(TopBoundaryObserver.contexts)
        assert np.array_equal(f["ordinal"][:], np.arange(4))
        assert np.all(f["energy_removed_J"][:] >= -1e-12)
        assert f["delta_w_rms"].shape == (4, 4)
        assert f["w_rms"].shape == (4, 48)
        assert f["full_fields/w_before"].shape == (4, 8, 8, 4)
        reconstructed = f["full_fields/w_before"][0] + f["full_fields/delta_w"][0]
        assert np.all(np.isfinite(reconstructed))


def test_energy_and_curl_are_exact_for_one_analytic_damping_call(tmp_path):
    h5py = pytest.importorskip("h5py")
    cfg = _control_config(nx=4, ny=5)
    sim = StormSimulation(cfg)
    sim.cfg.boundaries.damping_faces = 4
    x = np.asarray(sim.grid.xc)[:, None, None]
    y = np.asarray(sim.grid.yc)[None, :, None]
    sim.state.w[:] = 0.0
    for k in range(sim.grid.nz - 3, sim.grid.nz + 1):
        sim.state.w[:, :, k] = 2.0 + 1e-4 * x[:, :, 0] + 2e-4 * y[:, :, 0]
    before = sim.state.w.copy()
    observer = TopBoundaryObserver(sim, tmp_path / "analytic.h5", flush_every=1)
    apply_velocity_bcs(
        sim.state, sim.grid, sim.cfg, damping_observer=observer,
        context="pre_predictor", step=0, time_s=0.0, dt=0.2,
    )
    observer.close()
    with h5py.File(tmp_path / "analytic.h5") as f:
        idx = f["face_indices"][0].astype(int)
        multipliers = f["multipliers"][0]
        rho = np.asarray(sim.rho0_wface)[idx]
        dz = np.asarray(sim.grid.dz_c)
        dual = np.r_[0.5 * dz[0], 0.5 * (dz[:-1] + dz[1:]), 0.5 * dz[-1]]
        expected = []
        for j, k in enumerate(idx):
            field = before[:, :, k]
            expected.append(0.5 * rho[j] * sim.grid.dx * sim.grid.dy * dual[k]
                            * np.sum(field * field * (1.0 - multipliers[j] ** 2)))
        assert np.allclose(f["energy_removed_J"][0], expected, rtol=2e-14, atol=1e-9)
        xi = np.asarray(f["delta_xi_rms"][0])
        eta = np.asarray(f["delta_eta_rms"][0])
        assert np.any(xi > 0.0) and np.any(eta > 0.0)
        assert np.all(np.isfinite(xi)) and np.all(np.isfinite(eta))


@pytest.mark.parametrize("device", ["cpu", "gpu"])
def test_high_top_one_step_is_finite_and_lower_initial_state_matches(tmp_path, device):
    if device == "gpu":
        cp = pytest.importorskip("cupy")
        try:
            cp.zeros(1)
        except Exception as exc:
            pytest.skip(str(exc))
    control_cfg = _control_config(device=device)
    high_cfg = configure_high_top(control_cfg)
    control, high = StormSimulation(control_cfg), StormSimulation(high_cfg)
    to = high.backend.to_cpu
    assert np.array_equal(to(control.grid.zf), to(high.grid.zf[:49]))
    # Horizontal faces differ only in the extra vertical suffix; centered scalar
    # and velocity values below 15 km are otherwise exactly initialized.
    for name in ("theta", "qv", "ql", "qi", "qr", "qs", "qg", "qh"):
        assert np.array_equal(to(getattr(control.state, name)), to(getattr(high.state, name)[:, :, :48])), name
    high._step(0.05)
    for name in ("u", "v", "w", "p", "theta", "qv", "ql", "qi", "qr", "qs", "qg", "qh"):
        assert np.all(np.isfinite(to(getattr(high.state, name)))), name
