"""Analytic checks for the tornado-precursor scoring."""
import numpy as np
import pytest

from meteorological_flow.grid import Grid
from storm_dynamics.soundings import build_sounding
from storm_dynamics.config import HodographConfig
from storm_dynamics.vorticity_budget import horizontal_vorticity, streamwise_crosswise
from storm_dynamics.tornado_precursor import (GATES, environment_precursor, state_precursor)


def _grid(nx=16, ny=16, nz=8, Lz=2000.0):
    return Grid(nx=nx, ny=ny, nz=nz, Lx=8000.0, Ly=8000.0, Lz=Lz, periodic=True)


def _fields(grid, u=None, v=None, w=None):
    shape = (grid.nx, grid.ny, grid.nz)
    zero = np.zeros(shape)
    return (zero.copy() if u is None else u, zero.copy() if v is None else v,
            zero.copy() if w is None else w)


def test_streamwise_crosswise_is_a_rotation():
    """The split is a rotation of (xi, eta), so the magnitude must be conserved exactly."""
    grid = _grid()
    rng = np.random.default_rng(20260913)
    shape = (grid.nx, grid.ny, grid.nz)
    uc, vc, wc = (rng.normal(size=shape) for _ in range(3))
    xi, eta = horizontal_vorticity(uc, vc, wc, grid)
    sw, cw = streamwise_crosswise(uc, vc, xi, eta, grid, storm_motion=(3.0, -2.0))
    assert np.allclose(sw ** 2 + cw ** 2, xi ** 2 + eta ** 2, rtol=1e-12, atol=1e-14)


def test_unidirectional_shear_is_crosswise_not_streamwise():
    """u(z) only, storm drifting along x: omega_h points along +y, the flow along x.

    Perpendicular vectors -> the streamwise component is zero.  This is the textbook
    statement that straight-line shear gives a splitting pair, not a rotating updraft.
    """
    grid = _grid()
    z = np.asarray(grid.zc)
    uc = np.broadcast_to((5.0 + 0.01 * z)[None, None, :], (grid.nx, grid.ny, grid.nz)).copy()
    vc = np.zeros_like(uc)
    wc = np.zeros_like(uc)
    xi, eta = horizontal_vorticity(uc, vc, wc, grid)
    sw, cw = streamwise_crosswise(uc, vc, xi, eta, grid, storm_motion=(0.0, 0.0))
    assert np.allclose(sw, 0.0, atol=1e-12)
    assert np.abs(cw).max() > 1e-3


def test_alignment_is_one_when_updraft_gradient_follows_the_vortex_lines():
    """omega_h along +y by construction; make grad_h w also along +y -> cos(theta) = 1."""
    grid = _grid()
    z = np.asarray(grid.zc)
    y = np.asarray(grid.yc)
    uc = np.broadcast_to((0.01 * z)[None, None, :], (grid.nx, grid.ny, grid.nz)).copy()
    vc = np.zeros_like(uc)
    wc = np.broadcast_to((2e-3 * y)[None, :, None], uc.shape).copy()
    score = state_precursor(uc, vc, wc, grid, z_lo=float(z[1]), z_hi=float(z[-2]),
                            w_updraft=-1e9)
    assert score.values['alignment_mean'] == pytest.approx(1.0, abs=2e-2)

    wc_perp = np.broadcast_to((2e-3 * np.asarray(grid.xc))[:, None, None], uc.shape).copy()
    perp = state_precursor(uc, vc, wc_perp, grid, z_lo=float(z[1]), z_hi=float(z[-2]),
                           w_updraft=-1e9)
    assert abs(perp.values['alignment_mean']) < 2e-2


def test_anticyclonic_alignment_is_reported_negative_and_fails_the_gate():
    grid = _grid()
    z = np.asarray(grid.zc)
    y = np.asarray(grid.yc)
    uc = np.broadcast_to((0.01 * z)[None, None, :], (grid.nx, grid.ny, grid.nz)).copy()
    vc = np.zeros_like(uc)
    wc = np.broadcast_to((-2e-3 * y)[None, :, None], uc.shape).copy()   # gradient reversed
    score = state_precursor(uc, vc, wc, grid, z_lo=float(z[1]), z_hi=float(z[-2]),
                            w_updraft=-1e9)
    assert score.values['alignment_mean'] < 0
    assert score.gates['alignment'] is False
    assert 'NEGATIVE' in score.note


def test_environment_scores_a_quarter_circle_hodograph_above_a_straight_one():
    """A turning hodograph is the verbal precursor; it must score higher than straight shear."""
    grid = _grid(nz=48, Lz=15000.0)
    turning = build_sounding(grid, HodographConfig(kind='quarter_circle', U_max=30.0))
    straight = build_sounding(grid, HodographConfig(kind='unidirectional', U_max=30.0))
    a = environment_precursor(turning)
    b = environment_precursor(straight)
    assert a.values['srh_0_1km_m2_s2'] > b.values['srh_0_1km_m2_s2']
    assert a.values['streamwise_fraction_0_1km'] > b.values['streamwise_fraction_0_1km']


def test_environment_report_states_that_it_is_not_sufficient():
    """The note is load-bearing: attempt E raised SRH 254->648 for ~0% change."""
    grid = _grid(nz=48, Lz=15000.0)
    base = build_sounding(grid, HodographConfig(kind='quarter_circle', U_max=30.0))
    score = environment_precursor(base)
    assert 'not' in score.note and 'sufficient' in score.note
    assert score.kind == 'environment'
    assert 'alignment' not in score.gates          # cannot be evaluated without an updraft


def test_no_updraft_makes_the_conditional_scores_undefined_and_says_so():
    grid = _grid()
    uc, vc, wc = _fields(grid)
    uc = uc + np.asarray(grid.zc)[None, None, :] * 0.01
    score = state_precursor(uc, vc, wc, grid, w_updraft=1.0)
    assert not np.isfinite(score.values['alignment_updraft'])
    assert score.gates['alignment'] is False
    assert 'undefined' in score.note


def test_gates_are_documented_as_model_calibrated():
    assert set(GATES) >= {'omega_h_0_1km_s', 'streamwise_fraction', 'alignment_updraft'}
    import storm_dynamics.tornado_precursor as mod
    assert 'not' in mod.__doc__ and 'climatolog' in mod.__doc__
