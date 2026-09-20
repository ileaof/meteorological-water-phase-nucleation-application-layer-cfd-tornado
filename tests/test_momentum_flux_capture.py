"""Optional flux outputs must preserve tendencies and reconstruct their sum."""
from types import SimpleNamespace

import numpy as np
import pytest

from meteorological_flow.grid import Grid
from storm_dynamics.momentum import _u_tendency, _v_tendency
from storm_dynamics.vorticity_provenance import _frozen_u_tendency, _frozen_v_tendency


def reconstruct(flux, grid, component, periodic):
    shape = grid.u_shape if component == 'u' else grid.v_shape
    result = np.zeros(shape)
    normal_axis = 0 if component == 'u' else 1
    for axis, direction, spacing in ((0, 'x', grid.dx), (1, 'y', grid.dy)):
        f = np.moveaxis(flux[f'F{direction}_{component}'], axis, 0)
        target = np.moveaxis(result, axis, 0)
        if axis == normal_axis:
            target[1:-1] -= (f[1:] - f[:-1]) / spacing
            if periodic:
                wrap = -(f[0] - f[-1]) / spacing
                target[0] += wrap
                target[-1] += wrap
        elif periodic:
            target -= (np.roll(f, -1, axis=0) - f) / spacing
        else:
            target[1:-1] -= (f[2:-1] - f[1:-2]) / spacing
    f = flux[f'Fz_{component}']
    dz = grid.dz_c[None, None, :] if grid.stretched else grid.dz
    result -= (f[:, :, 1:] - f[:, :, :-1]) / dz
    return result


@pytest.mark.parametrize('periodic', [False, True])
@pytest.mark.parametrize('order', [1, 2])
@pytest.mark.parametrize('component', ['u', 'v'])
def test_optional_fluxes_are_passive_and_complete(periodic, order, component):
    grid = Grid(8, 7, 6, 800., 700., 600., periodic=periodic, z_stretch=1.05)
    rng = np.random.default_rng(31)
    total = tuple(rng.normal(size=shape) for shape in (grid.u_shape, grid.v_shape, grid.w_shape))
    label = tuple(rng.normal(size=shape) for shape in (grid.u_shape, grid.v_shape, grid.w_shape))
    originals = [a.copy() for a in (*total, *label)]
    state = SimpleNamespace(**dict(zip(('u', 'v', 'w'), total)))
    native = _u_tendency if component == 'u' else _v_tendency
    frozen = _frozen_u_tendency if component == 'u' else _frozen_v_tendency
    for function, args in ((native, (state, grid, order, periodic)),
                           (frozen, (label, total, grid, order, periodic))):
        expected = function(*args)
        fluxes = {}
        captured = function(*args, fluxes=fluxes)
        assert np.array_equal(expected, captured)
        reconstructed = reconstruct(fluxes, grid, component, periodic)
        np.testing.assert_allclose(reconstructed, captured, rtol=1e-13, atol=1e-14)
        assert np.all(fluxes[f'Fz_{component}'][:, :, (0, -1)] == 0)
    for before, after in zip(originals, (*total, *label)):
        assert np.array_equal(before, after)
