"""Cold-pool diagnostics (storm_dynamics.coldpool) on a controlled cold slab."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from meteorological_flow.grid import Grid
from meteorological_flow.state import FlowState
from storm_dynamics import coldpool as cp

_G = 9.81


def _grid(backend=None):
    return Grid(nx=32, ny=32, nz=30, Lx=6400.0, Ly=6400.0, Lz=3000.0, z_stretch=1.0,
                periodic=True, backend=backend)


def _state_with_cold_patch(g, dT=3.0, qv=0.01, h=1000.0):
    st = FlowState.zeros(g)
    xp = g.xp
    st.theta = 300.0 + g.zeros_c()
    st.qv = qv + g.zeros_c()
    z = np.asarray(g.backend.to_cpu(g.zc))
    kmask = xp.asarray((z < h), dtype=float)[None, None, :]
    # cold patch over a small central footprint (so the horizontal mean stays ~ambient)
    patch = g.zeros_c()
    patch[12:20, 12:20, :] = 1.0
    st.theta = st.theta - dT * patch * kmask
    st.diagnose = lambda *a, **k: None       # not needed; coldpool uses theta/qv directly
    return st


def test_thetav_pert_and_intensity():
    g = _grid(); dT, qv, h = 3.0, 0.01, 1000.0
    st = _state_with_cold_patch(g, dT=dT, qv=qv, h=h)
    _, tvp = cp.coldpool_buoyancy(st, g)
    tvmin = float(np.asarray(g.backend.to_cpu(tvp)).min())
    assert tvmin < -2.5           # ~ -dT*(1+0.61 qv), reduced slightly by the horizontal mean
    C = cp.coldpool_intensity(st, g)
    Cmax = float(np.asarray(g.backend.to_cpu(C)).max())
    thetav0 = 300.0 * (1.0 + 0.61 * qv)
    C_expected = np.sqrt(2.0 * _G * (dT * (1 + 0.61 * qv) / thetav0) * h)   # ~14 m/s
    assert np.isclose(Cmax, C_expected, rtol=0.15)


def test_report_keys_and_cold_area():
    g = _grid()
    st = _state_with_cold_patch(g)
    rep = cp.coldpool_report(st, g)
    for k in ("coldpool_min_thetav_pert_K", "coldpool_intensity_C_m_s", "coldpool_area_fraction",
              "gust_front_thetav_grad_K_m", "gust_front_convergence_s", "downdraft_w_min_m_s"):
        assert k in rep
    assert 0.0 < rep["coldpool_area_fraction"] < 0.3          # small central patch
    assert rep["coldpool_min_thetav_pert_K"] < -2.5


def test_no_cold_pool_when_uniform():
    g = _grid()
    st = FlowState.zeros(g); st.theta = 300.0 + g.zeros_c(); st.qv = 0.01 + g.zeros_c()
    rep = cp.coldpool_report(st, g)
    assert abs(rep["coldpool_min_thetav_pert_K"]) < 1e-6
    assert rep["coldpool_intensity_C_m_s"] < 1e-6


if __name__ == "__main__":
    for nm, fn in sorted(globals().items()):
        if nm.startswith("test_") and callable(fn):
            fn(); print("ok", nm)
    print("ALL COLD-POOL TESTS PASSED")
