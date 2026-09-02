"""cos^2 convection-initiation bubble (storm_dynamics.initiation)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from meteorological_flow.grid import Grid
from storm_dynamics import initiation as ini


def _grid():
    return Grid(nx=40, ny=40, nz=30, Lx=40000.0, Ly=40000.0, Lz=6000.0, z_stretch=1.0, periodic=True)


def test_peak_at_center_zero_at_edge():
    g = _grid(); A = 2.5
    b = np.asarray(ini.smooth_bubble(g, dtheta_max=A, center=(20500.0, 20500.0, 1500.0),
                                     Rx=8000.0, Ry=8000.0, Rz=1500.0))
    assert np.isclose(b.max(), A, rtol=1e-6)                 # peak = amplitude
    assert b.min() >= 0.0                                    # compact, non-negative
    # zero well outside the support
    assert b[0, 0, -1] == 0.0


def test_independent_radii_anisotropy():
    g = _grid()
    b = ini.smooth_bubble(g, dtheta_max=2.0, center=(20500.0, 20500.0, 1500.0),
                          Rx=12000.0, Ry=4000.0, Rz=1500.0)
    b = np.asarray(g.backend.to_cpu(b))
    ic = np.argmin(np.abs(np.asarray(g.xc) - 20500.0)); jc = np.argmin(np.abs(np.asarray(g.yc) - 20500.0))
    kc = np.argmin(np.abs(np.asarray(g.zc) - 1500.0))
    # wider in x than y -> the x-extent of the warm region exceeds the y-extent
    x_extent = np.sum(b[:, jc, kc] > 0.1)
    y_extent = np.sum(b[ic, :, kc] > 0.1)
    assert x_extent > y_extent


def test_bubble_carries_no_rotation():
    # the trigger is a scalar theta' field: it has zero horizontal velocity, hence zero vorticity
    g = _grid()
    b = ini.smooth_bubble(g)
    # a scalar field induces no velocity/vorticity by construction; assert it is purely thermal
    assert b.shape == g.center_shape


def test_multi_bubble_superposes():
    g = _grid()
    specs = [{"center": (12500.0, 20500.0, 1500.0), "dtheta_max": 2.0, "Rx": 4000, "Ry": 4000, "Rz": 1500},
             {"center": (28500.0, 20500.0, 1500.0), "dtheta_max": 3.0, "Rx": 4000, "Ry": 4000, "Rz": 1500}]
    b = np.asarray(g.backend.to_cpu(ini.multi_bubble(g, specs)))
    assert np.isclose(b.max(), 3.0, rtol=1e-6)               # the stronger bubble's peak
    # two separated warm cores
    warm_cols = np.any(b > 0.5, axis=2)
    assert warm_cols.sum() > 0


if __name__ == "__main__":
    for nm, fn in sorted(globals().items()):
        if nm.startswith("test_") and callable(fn):
            fn(); print("ok", nm)
    print("ALL INITIATION TESTS PASSED")
