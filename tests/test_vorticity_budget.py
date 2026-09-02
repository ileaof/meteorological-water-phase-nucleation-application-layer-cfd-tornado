"""Analytic verification of the vertical-vorticity budget (storm_dynamics.vorticity_budget).

Each production term is checked against a hand-computed value on a controlled centred-velocity
field, then CPU/GPU parity.  Central differences are exact for these low-order polynomial fields at
interior points, so tolerances are tight."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from meteorological_flow.grid import Grid
from storm_dynamics import vorticity_budget as vb


def _grid(nx=32, backend=None):
    return Grid(nx=nx, ny=nx, nz=24, Lx=3200.0, Ly=3200.0, Lz=2400.0, z_stretch=1.0,
                periodic=True, backend=backend)


def _mesh(g):
    xp = g.xp
    xc = xp.asarray(g.xc); yc = xp.asarray(g.yc); zc = xp.asarray(g.zc)
    X = xc[:, None, None] + xp.zeros((g.nx, g.ny, g.nz))
    Y = yc[None, :, None] + xp.zeros((g.nx, g.ny, g.nz))
    Z = zc[None, None, :] + xp.zeros((g.nx, g.ny, g.nz))
    return X, Y, Z


def _interior(a):
    return np.asarray(a)[6:-6, 6:-6, 6:-6]


def test_zeta_solid_body_rotation():
    g = _grid()
    X, Y, Z = _mesh(g); xp = g.xp
    xc0, yc0 = 1600.0, 1600.0
    Om = 3e-3
    uc = -Om * (Y - yc0); vc = Om * (X - xc0); wc = xp.zeros_like(X)
    t = vb.budget_from_velocity(uc, vc, wc, g)
    # zeta = dv/dx - du/dy = Om - (-Om) = 2 Om, uniform
    assert np.allclose(_interior(t["zeta"]), 2 * Om, rtol=1e-6, atol=1e-9)


def test_stretching_term():
    g = _grid()
    X, Y, Z = _mesh(g); xp = g.xp
    Om, a = 3e-3, 2e-3
    uc = -Om * (Y - 1600.0); vc = Om * (X - 1600.0); wc = a * Z          # dw/dz = a
    t = vb.budget_from_velocity(uc, vc, wc, g)
    # stretching = zeta * dw/dz = 2 Om * a
    assert np.allclose(_interior(t["stretching"]), 2 * Om * a, rtol=1e-6, atol=1e-12)


def test_tilting_term():
    g = _grid()
    X, Y, Z = _mesh(g); xp = g.xp
    V0, A = 4e-3, 5e-3
    uc = xp.zeros_like(X); vc = V0 * Z; wc = A * X                        # dv/dz=V0, dw/dx=A
    t = vb.budget_from_velocity(uc, vc, wc, g)
    # xi = dw/dy - dv/dz = -V0 ; eta = du/dz - dw/dx = -A
    # tilting = xi dw/dx + eta dw/dy = (-V0)(A) + (-A)(0) = -V0 A
    assert np.allclose(_interior(t["tilting"]), -V0 * A, rtol=1e-6, atol=1e-12)


def test_baroclinic_term():
    g = _grid()
    X, Y, Z = _mesh(g); xp = g.xp
    rho0, drho, dp, L = 1.0, 0.02, 50.0, 3200.0
    uc = xp.zeros_like(X); vc = xp.zeros_like(X); wc = xp.zeros_like(X)
    rho = rho0 + drho * (X - 1600.0) / L                                 # drho/dx = drho/L
    p = 1.0e5 + dp * (Y - 1600.0) / L                                    # dp/dy = dp/L
    t = vb.budget_from_velocity(uc, vc, wc, g, rho=rho, p=p)
    # B = (1/rho^2)(drho/dx dp/dy - drho/dy dp/dx) = (drho*dp/L^2)/rho^2 (rho ~ rho0 near centre)
    expected = (drho * dp / L**2)
    got = _interior(t["baroclinic"]) * (_interior(rho)**2)               # multiply back rho^2
    assert np.allclose(got, expected, rtol=2e-2)                          # rho varies mildly across cell


def test_advection_term():
    g = _grid()
    X, Y, Z = _mesh(g); xp = g.xp
    # zeta field from vc = 0.5*k*x^2 -> dv/dx = k x -> zeta = k x ; uniform u advects it
    k, U = 1e-6, 7.0
    uc = U + xp.zeros_like(X); vc = 0.5 * k * X * X; wc = xp.zeros_like(X)
    t = vb.budget_from_velocity(uc, vc, wc, g)
    # zeta = k x ; advection = -(u dzeta/dx) = -U*k
    assert np.allclose(_interior(t["advection"]), -U * k, rtol=1e-5, atol=1e-12)


def test_streamwise_crosswise_alignment():
    g = _grid()
    X, Y, Z = _mesh(g); xp = g.xp
    # uniform flow in +x (vc=0 everywhere, no cross-flow), horizontal vorticity along +x via
    # dw/dy: xi = dw/dy - dv/dz = W, eta = du/dz - dw/dx = 0 -> purely streamwise.
    W = 3e-3
    uc = 10.0 + xp.zeros_like(X); vc = xp.zeros_like(X); wc = W * Y
    xi, eta = vb.horizontal_vorticity(uc, vc, wc, g)                      # full form
    sw, cw = vb.streamwise_crosswise(uc, vc, xi, eta, g, storm_motion=(0.0, 0.0))
    assert np.allclose(_interior(sw), W, rtol=1e-6, atol=1e-9)
    assert np.allclose(_interior(cw), 0.0, atol=1e-9)


def test_tendency_is_sum_of_terms():
    g = _grid()
    X, Y, Z = _mesh(g); xp = g.xp
    uc = -3e-3 * (Y - 1600.0) + 2.0; vc = 3e-3 * (X - 1600.0); wc = 1e-3 * Z
    t = vb.budget_from_velocity(uc, vc, wc, g, Km=30.0)
    s = t["advection"] + t["stretching"] + t["tilting"] + t["baroclinic"] + t["divergence"] \
        + t["diffusion"] + t["friction"]
    assert np.allclose(np.asarray(t["tendency"]), np.asarray(s), rtol=1e-10, atol=1e-14)


def test_tilting_efficiency_factorisation():
    g = _grid()
    X, Y, Z = _mesh(g); xp = g.xp
    # omega_h along +x (xi=W via dw/dy) and grad_h w along +x -> perfectly ALIGNED (cos=1)
    W, A = 3e-3, 4e-3
    uc = xp.zeros_like(X); vc = xp.zeros_like(X); wc = W * Y + A * X
    e = vb.tilting_efficiency(uc, vc, wc, g)
    # xi = dw/dy = W, eta = -dw/dx = -A ; grad_h w = (A, W)
    # tilting = xi*A + eta*W = W*A - A*W = 0 -> orthogonal, alignment ~ 0
    assert abs(float(np.asarray(_interior(e["alignment"])).mean())) < 1e-6
    # genuine alignment: shear vorticity xi=W along +x, with a WEAK dw/dx (A2 << W) so the
    # w-gradient's own contribution to eta (=-dw/dx) does not rotate omega_h away from +x.
    # Then alignment = W/sqrt(W^2+A2^2) -> ~1.
    A2 = 1e-4
    wc2 = A2 * X
    vc2 = -W * Z                     # dv/dz = -W -> xi = dw/dy - dv/dz = 0 + W = W
    e2 = vb.tilting_efficiency(uc, vc2, wc2, g)
    expected = W / np.sqrt(W ** 2 + A2 ** 2)
    got = float(np.asarray(_interior(e2["alignment"])).mean())
    assert got > 0.99 and abs(got - expected) < 1e-3


def test_baroclinic_horizontal_generation():
    g = _grid()
    X, Y, Z = _mesh(g); xp = g.xp
    rho0, drho, L = 1.1, 0.05, 3200.0
    rho = rho0 + drho * (X - 1600.0) / L                 # drho/dx = drho/L, drho/dy = 0
    Gx, Gy, Gmag = vb.baroclinic_horizontal_generation(rho, g)
    # Gx = -g/rho drho/dy = 0 ; Gy = g/rho drho/dx = 9.81/rho * drho/L
    assert np.allclose(_interior(Gx), 0.0, atol=1e-9)
    expected = 9.81 / _interior(rho) * (drho / L)
    assert np.allclose(_interior(Gy), expected, rtol=1e-3)


def test_dominant_mechanism_picks_baroclinic():
    g = _grid()
    X, Y, Z = _mesh(g); xp = g.xp
    uc = xp.zeros_like(X); vc = xp.zeros_like(X); wc = xp.zeros_like(X)
    rho = 1.0 + 0.05 * (X - 1600.0) / 3200.0
    p = 1.0e5 + 80.0 * (Y - 1600.0) / 3200.0
    t = vb.budget_from_velocity(uc, vc, wc, g, rho=rho, p=p)
    name, _ = vb.dominant_mechanism(t, g, 0.0, 2400.0)
    assert name == "baroclinic"


def _gpu_available():
    try:
        from meteorological_flow.backend import get_backend
        return get_backend("gpu") is not None
    except Exception:
        return False


@pytest.mark.skipif(not _gpu_available(), reason="no GPU backend")
def test_cpu_gpu_parity():
    from meteorological_flow.backend import get_backend
    gc = _grid(); gg = _grid(backend=get_backend("gpu"))
    for g in (gc, gg):
        X, Y, Z = _mesh(g); xp = g.xp
        uc = -3e-3 * (Y - 1600.0) + 2.0; vc = 3e-3 * (X - 1600.0); wc = 1e-3 * Z
        rho = 1.0 + 0.02 * (X - 1600.0) / 3200.0; p = 1e5 + 40.0 * (Y - 1600.0) / 3200.0
        g._budget = vb.budget_from_velocity(uc, vc, wc, g, rho=rho, p=p, Km=30.0)
    for k in ("advection", "stretching", "tilting", "baroclinic", "divergence", "diffusion", "tendency"):
        a = np.asarray(gc._budget[k]); b = np.asarray(gg.backend.to_cpu(gg._budget[k]))
        assert np.allclose(a, b, rtol=1e-5, atol=1e-10), k


if __name__ == "__main__":
    for nm, fn in sorted(globals().items()):
        if nm.startswith("test_") and callable(fn) and "gpu" not in nm:
            fn(); print("ok", nm)
    print("ALL VORTICITY-BUDGET TESTS PASSED")
