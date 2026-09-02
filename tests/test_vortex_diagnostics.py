"""Analytic verification of the vortex diagnostics (storm_dynamics.vortex_diagnostics) using a
Lamb-Oseen vortex (smooth, so finite differences are clean)."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from meteorological_flow.grid import Grid
from storm_dynamics import vortex_diagnostics as vd


def _grid(backend=None):
    return Grid(nx=80, ny=80, nz=4, Lx=4000.0, Ly=4000.0, Lz=400.0, z_stretch=1.0,
                periodic=True, backend=backend)


def _lamb_oseen(g, Gamma=2.0e4, rc=200.0):
    xp = g.xp
    xc = yc = 2000.0
    X = xp.asarray(g.xc)[:, None]; Y = xp.asarray(g.yc)[None, :]
    dx = X - xc; dy = Y - yc
    r = xp.sqrt(dx * dx + dy * dy) + 1e-9
    vth = (Gamma / (2 * xp.pi * r)) * (1.0 - xp.exp(-(r * r) / (rc * rc)))
    phi = xp.arctan2(dy, dx)
    uc2 = -vth * xp.sin(phi); vc2 = vth * xp.cos(phi)
    zeta2 = g._central_x(vc2[:, :, None])[:, :, 0] - g._central_y(uc2[:, :, None])[:, :, 0]
    return uc2, vc2, zeta2, xc, yc


def test_center_at_vortex():
    g = _grid()
    _, _, zeta2, xc, yc = _lamb_oseen(g)
    ic, jc, cx, cy = vd.find_vortex_center(zeta2, g)
    assert abs(cx - xc) <= 2 * g.dx and abs(cy - yc) <= 2 * g.dy


def test_total_circulation():
    g = _grid(); Gamma = 2.0e4
    _, _, zeta2, xc, yc = _lamb_oseen(g, Gamma=Gamma)
    # a disk large enough to enclose essentially all the vorticity -> Gamma_inf
    gam = vd.circulation(zeta2, g, (0, 0, xc, yc), radius_m=1600.0)
    assert np.isclose(gam, Gamma, rtol=0.08)


def test_tangential_peak_and_core():
    g = _grid(); Gamma, rc = 2.0e4, 200.0
    uc2, vc2, zeta2, xc, yc = _lamb_oseen(g, Gamma=Gamma, rc=rc)
    vth, vrad, r = vd.tangential_radial(uc2, vc2, g, (0, 0, xc, yc))
    _, _, vth_max, core = vd.tangential_profile(vth, r, g, r_max_m=1200.0, nbins=40)
    vth_expected = 0.638 * Gamma / (2 * np.pi * rc)      # Lamb-Oseen peak ~10.15 m/s
    assert np.isclose(vth_max, vth_expected, rtol=0.12)
    assert np.isclose(core, 1.12 * rc, rtol=0.35)         # peak at r ~ 1.12 rc


def test_pressure_deficit_negative_core():
    g = _grid(); xp = g.xp
    X = xp.asarray(g.xc)[:, None]; Y = xp.asarray(g.yc)[None, :]
    r2 = (X - 2000.0) ** 2 + (Y - 2000.0) ** 2
    p2 = -300.0 * xp.exp(-r2 / (300.0 ** 2))              # a 300 Pa low at the centre
    ddef = vd.pressure_deficit(p2, g, center=(40, 40))
    assert ddef < -250.0


def _gpu_available():
    try:
        from meteorological_flow.backend import get_backend
        return get_backend("gpu") is not None
    except Exception:
        return False


@pytest.mark.skipif(not _gpu_available(), reason="no GPU backend")
def test_cpu_gpu_parity():
    from meteorological_flow.backend import get_backend
    res = []
    for backend in (None, get_backend("gpu")):
        g = _grid(backend=backend)
        uc2, vc2, zeta2, xc, yc = _lamb_oseen(g)
        ic, jc, cx, cy = vd.find_vortex_center(zeta2, g)      # exercises the GPU argmax path
        gam = vd.circulation(zeta2, g, (0, 0, xc, yc), radius_m=1600.0)
        vth, _, r = vd.tangential_radial(uc2, vc2, g, (0, 0, xc, yc))
        _, _, vmax, core = vd.tangential_profile(vth, r, g, r_max_m=1200.0, nbins=40)
        res.append((cx, cy, gam, vmax, core))
    assert np.allclose(res[0], res[1], rtol=1e-4)


if __name__ == "__main__":
    for nm, fn in sorted(globals().items()):
        if nm.startswith("test_") and callable(fn) and "gpu" not in nm:
            fn(); print("ok", nm)
    print("ALL VORTEX-DIAGNOSTIC TESTS PASSED")
