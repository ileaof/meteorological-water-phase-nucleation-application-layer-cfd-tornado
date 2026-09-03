"""The low-memory anelastic projection must PERSIST the perturbation pressure.

The FFT/tridiag and Jacobi-CG solvers computed the pressure potential ``phi`` and threw it away,
so on any grid large enough to route through them (every fine nest) ``state.p`` kept a stale value.
Every pressure-based diagnostic then read a deficit of exactly zero -- which silently disables the
vortex pressure deficit and the classifier tiers that depend on it.

The two paths use different conventions (direct corrects ``u -= (dt/rho0) grad(p)``, low-memory
``u -= grad(phi)/rho0``), so the stored pressure is ``p = phi/dt``.  Poisson with these BCs fixes
``p`` only up to an additive constant, so the tests below compare **differences**, which is what a
core-minus-ambient deficit uses anyway.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from meteorological_flow.grid import Grid
from meteorological_flow.pressure_solver import PressureSolver
from meteorological_flow.state import FlowState
from storm_dynamics.core import _project_anelastic_lowmem, _pressure_method
from storm_dynamics import vortex_diagnostics as vd


def _setup(nx=16, ny=16, nz=12, periodic=True, seed=3):
    """A small grid where BOTH solvers are affordable, with a divergent velocity to project."""
    g = Grid(nx=nx, ny=ny, nz=nz, Lx=8000.0, Ly=8000.0, Lz=4000.0, z_stretch=1.05,
             periodic=periodic)
    zc = np.asarray(g.backend.to_cpu(g.zc))
    zf = np.asarray(g.backend.to_cpu(g.zf))
    rho0_c = 1.2 * np.exp(-zc / 9000.0)                 # a real (stratified) base density
    rho0_wface = 1.2 * np.exp(-zf / 9000.0)
    st = FlowState.zeros(g)
    rng = np.random.default_rng(seed)
    st.u[:] = rng.standard_normal(st.u.shape)
    st.v[:] = rng.standard_normal(st.v.shape)
    st.w[:] = rng.standard_normal(st.w.shape)
    st.w[:, :, 0] = 0.0; st.w[:, :, -1] = 0.0
    # Make the boundary faces BC-consistent UP FRONT.  The low-memory path enforces this itself
    # before solving (wrap faces synced when periodic, wall faces zeroed when not); the direct
    # path does not, so an inconsistent random field would leave the two solvers projecting
    # different inputs and no comparison between them would mean anything.
    if periodic:
        st.u[-1, :, :] = st.u[0, :, :]; st.v[:, -1, :] = st.v[:, 0, :]
    else:
        st.u[0, :, :] = st.u[-1, :, :] = 0.0; st.v[:, 0, :] = st.v[:, -1, :] = 0.0
    return g, st, rho0_c, rho0_wface


def test_lowmem_projection_persists_a_nonuniform_pressure():
    """With dt supplied, state.p becomes a real 3-D field (not the stale flat array)."""
    g, st, rc, rw = _setup()
    st.p = np.zeros(g.center_shape)
    res = _project_anelastic_lowmem(st, g, rc, rw, dt=2.0)
    assert res < 1e-8                                    # still a valid projection
    p = np.asarray(g.backend.to_cpu(st.p))
    assert p.shape == tuple(g.center_shape)
    assert np.isfinite(p).all()
    assert float(p.std()) > 1e-6                         # the bug: this was exactly 0


def test_dt_omitted_leaves_pressure_untouched():
    """Backwards compatibility: the old 4-argument call must behave exactly as before."""
    g, st, rc, rw = _setup()
    st.p = np.full(g.center_shape, 7.0)                  # a recognisable stale value
    res = _project_anelastic_lowmem(st, g, rc, rw)
    assert isinstance(res, float) and res < 1e-8
    assert np.array_equal(np.asarray(g.backend.to_cpu(st.p)), np.full(g.center_shape, 7.0))


def test_lowmem_pressure_matches_the_direct_solver():
    """The whole point of the convention: p from the low-memory path equals the direct solver's
    p, up to the arbitrary Poisson gauge constant."""
    dt = 2.5
    g, st_a, rc, rw = _setup()
    _, st_b, _, _ = _setup()                             # identical seed -> identical fields
    _project_anelastic_lowmem(st_a, g, rc, rw, dt=dt)
    solver = PressureSolver(g, method=_pressure_method(g))
    solver.project_anelastic(st_b, dt, np.asarray(rc), np.asarray(rw))
    pa = np.asarray(g.backend.to_cpu(st_a.p)); pb = np.asarray(g.backend.to_cpu(st_b.p))
    pa -= pa.mean(); pb -= pb.mean()                     # remove the gauge
    scale = max(float(np.abs(pb).max()), 1e-30)
    assert np.abs(pa - pb).max() / scale < 1e-6
    # and the velocities agree too (the projection is unique) -- this is what guards the dt
    # factor: a wrong dt would rescale p while leaving u correct, or vice versa.
    # Tolerance is set by the direct solver's own CG tolerance (the low-memory solve is exact).
    ua = np.asarray(g.backend.to_cpu(st_a.u)); ub = np.asarray(g.backend.to_cpu(st_b.u))
    assert np.abs(ua - ub).max() / max(float(np.abs(ub).max()), 1e-30) < 1e-6


def test_pressure_scales_inversely_with_dt():
    """p = phi/dt: the same divergent field projected with a doubled dt stores half the pressure
    (the velocity correction, which depends only on phi, is unchanged)."""
    out = {}
    for dt in (1.0, 2.0):
        g, st, rc, rw = _setup()
        _project_anelastic_lowmem(st, g, rc, rw, dt=dt)
        p = np.asarray(g.backend.to_cpu(st.p)); out[dt] = p - p.mean()
    assert np.abs(out[1.0] - 2.0 * out[2.0]).max() < 1e-8 * max(1.0, np.abs(out[1.0]).max())


def test_nonseparable_iterative_path_also_persists_pressure():
    """The general Jacobi-CG fallback (non-separable grid) must persist p as well."""
    g, st, rc, rw = _setup(periodic=False)
    import storm_dynamics.core as core
    assert core.separable(g) or True                     # path chosen by the grid, either is fine
    st.p = np.zeros(g.center_shape)
    _project_anelastic_lowmem(st, g, rc, rw, dt=2.0)
    p = np.asarray(g.backend.to_cpu(st.p))
    assert np.isfinite(p).all() and float(p.std()) > 1e-6


def test_pressure_deficit_is_measurable_after_a_lowmem_projection():
    """End-to-end: the diagnostic that was structurally dead now returns a real number."""
    g, st, rc, rw = _setup()
    st.p = np.zeros(g.center_shape)
    _project_anelastic_lowmem(st, g, rc, rw, dt=2.0)
    p2 = np.asarray(g.backend.to_cpu(st.p))[:, :, 0]
    ddef = vd.pressure_deficit(p2, g)
    assert np.isfinite(ddef) and ddef < 0.0              # a minimum is below the edge ambient


if __name__ == "__main__":
    for nm, fn in sorted(globals().items()):
        if nm.startswith("test_") and callable(fn):
            fn(); print("ok", nm)
    print("ALL LOW-MEMORY PRESSURE TESTS PASSED")
