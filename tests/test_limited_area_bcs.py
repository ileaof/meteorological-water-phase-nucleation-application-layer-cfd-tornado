"""Limited-area lateral boundary conditions actually reaching the run.

A real case is a limited-area problem: the lateral boundaries must be driven toward the
analysed environment by the Davies zone.  Three separate defects each defeated that on their
own, and all three were silent:

1. ``driver._make_grid`` built ``periodic=True``, so advection wrapped and the storm ingested
   its own cold pool instead of environmental air.
2. Fixing (1) was cosmetic -- ``StormSimulation`` builds its OWN grid from ``build_storm_config``
   (the preprocess grid only supplies dimensions), so the run stayed periodic anyway.
3. ``run_multilevel_real_case`` never applied the lateral relaxation at all: ``run_case`` calls
   it every step, but the cascade path delegates to ``run_multilevel_nest``, whose loop had no
   hook.  So every real-case CASCADE ran with unconstrained laterals.

And the engine had no open ``y`` boundary -- only ``free_slip``/``wall`` (both pin v=0) or
``periodic`` -- so a Davies zone nudging v toward an environmental profile with v != 0 was
fighting a wall, and meridional inflow was impossible in principle.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from meteorological_flow.boundary_conditions import apply_velocity_bcs
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import Grid, StormSimulation
from storm_dynamics import nesting as nst


def _sim(**kw):
    cfg = build_storm_config(preset="storm", nx=16, ny=16, nz=12, Lx=16000.0, Ly=16000.0,
                             Lz=8000.0, duration=1.0, dt_max=3.0, device="cpu", **kw)
    return StormSimulation(cfg)


# ------------------------------------------------------------------ open y boundary
def test_y_outflow_is_zero_gradient_not_a_wall():
    """y='outflow' copies the inner slab; free_slip/wall pin v=0 (which fights the Davies zone)."""
    sim = _sim(periodic=False)
    g = sim.grid
    xp = g.xp
    sim.state.v[...] = xp.asarray(np.linspace(1.0, 5.0, g.v_shape[1])[None, :, None]
                                  * np.ones(g.v_shape))

    sim.cfg.boundaries.y = "free_slip"
    apply_velocity_bcs(sim.state, g, sim.cfg)
    v_wall = np.asarray(g.backend.to_cpu(sim.state.v))
    assert v_wall[:, 0, :].max() == 0.0 and v_wall[:, -1, :].max() == 0.0   # pinned to zero

    sim.state.v[...] = xp.asarray(np.linspace(1.0, 5.0, g.v_shape[1])[None, :, None]
                                  * np.ones(g.v_shape))
    sim.cfg.boundaries.y = "outflow"
    apply_velocity_bcs(sim.state, g, sim.cfg)
    v_open = np.asarray(g.backend.to_cpu(sim.state.v))
    assert v_open[:, 0, :].max() > 0.0                                      # NOT pinned
    assert np.allclose(v_open[:, 0, :], v_open[:, 1, :])                    # zero normal gradient
    assert np.allclose(v_open[:, -1, :], v_open[:, -2, :])


def test_y_default_is_unchanged():
    """The new option must not move the default: y stays free_slip."""
    from meteorological_flow.config import BoundaryConfig
    assert BoundaryConfig().y == "free_slip"


# ------------------------------------------------------------------ the cascade hook
def test_run_multilevel_nest_calls_the_parent_hook_every_parent_step():
    """Without this hook a real-case cascade steps the parent with NO lateral BC at all."""
    sim = _sim()
    for _ in range(3):
        sim._step(sim._dt())

    calls = []

    def hook(parent, dt):
        calls.append((float(dt), int(parent.grid.nx)))

    spec = lambda g: nst.NestSpec.aligned(g, i0=2, j0=2, ncx=8, ncy=8, refine=3, nz=g.nz)
    window = 4.0 * float(sim._dt())
    nst.run_multilevel_nest(sim, [spec], window=window, parent_hook=hook)

    assert calls, "parent_hook was never called"
    assert all(dt > 0 for dt, _ in calls)
    assert all(nx == sim.grid.nx for _, nx in calls)     # called with the PARENT, not a nest


def test_parent_hook_defaults_to_none_and_changes_nothing():
    """Opt-in: omitting the hook must reproduce the previous trajectory exactly."""
    a = _sim(); b = _sim()
    for s in (a, b):
        for _ in range(2):
            s._step(s._dt())
    spec = lambda g: nst.NestSpec.aligned(g, i0=2, j0=2, ncx=8, ncy=8, refine=3, nz=g.nz)
    w = 3.0 * float(a._dt())
    sims_a, _ = nst.run_multilevel_nest(a, [spec], window=w)
    sims_b, _ = nst.run_multilevel_nest(b, [spec], window=w, parent_hook=None)
    ua = np.asarray(sims_a[-1].grid.backend.to_cpu(sims_a[-1].state.u))
    ub = np.asarray(sims_b[-1].grid.backend.to_cpu(sims_b[-1].state.u))
    assert np.array_equal(ua, ub)


# ------------------------------------------------------------------ the real-case wiring
def test_real_case_simulation_is_limited_area_not_periodic():
    """The grid the simulation RUNS on must be non-periodic with open lateral faces."""
    cfg = build_storm_config(preset="storm", nx=16, ny=16, nz=12, Lx=16000.0, Ly=16000.0,
                             Lz=8000.0, duration=1.0, dt_max=3.0, device="cpu", periodic=False)
    cfg.sim.boundaries.x_west = cfg.sim.boundaries.x_east = "outflow"
    cfg.sim.boundaries.y = "outflow"
    sim = StormSimulation(cfg)
    assert getattr(sim.grid, "periodic", True) is False
    assert sim.cfg.boundaries.x_west == "outflow" and sim.cfg.boundaries.y == "outflow"


# ------------------------------------------------------------------ open-boundary projection
def _poisson_case(seed=0, nx=32, ny=32, nz=24):
    rng = np.random.default_rng(seed)
    dx = dy = 600.0
    zf = np.linspace(0.0, 15000.0, nz + 1)
    dzc = np.diff(zf); zc = 0.5 * (zf[1:] + zf[:-1]); dzf = np.diff(zc)
    rc = 1.2 * np.exp(-zc / 8500.0); rw = np.interp(zf, zc, rc)
    u = rng.normal(0, 5, (nx + 1, ny, nz))
    v = rng.normal(0, 5, (nx, ny + 1, nz))
    w = rng.normal(0, 2, (nx, ny, nz + 1))
    return u, v, w, rc, rw, dx, dy, dzc, dzf


def _scaled_divergence(u, v, w, rc, rw, dx, dy, dzc):
    from storm_dynamics.pressure_fft import anelastic_divergence
    d = np.abs(anelastic_divergence(u, v, w, rc, rw, dx, dy, dzc))
    return d / (rc.max() * max(np.abs(u).max(), 1e-9) / dx)


def test_open_lateral_projection_is_divergence_free_including_the_boundary_face():
    """The lowmem solver assumed SOLID WALLS on every non-periodic grid.

    That is right for a nest (walls + sponge) but wrong for a limited-area parent, where the
    boundary normal velocity carries the environmental inflow: zeroing it each step and letting
    the BC restore it left that face divergent (measured 3.5e-01 against 5.9e-06 periodic).
    """
    from storm_dynamics.pressure_fft import project_anelastic_fft
    u, v, w, rc, rw, dx, dy, dzc, dzf = _poisson_case()
    project_anelastic_fft(u, v, w, rc, rw, dx, dy, dzc, dzf, periodic_h=False, lateral="open")
    d = _scaled_divergence(u, v, w, rc, rw, dx, dy, dzc)
    assert d.max() < 1e-12, "open-boundary projection must be divergence-free everywhere"
    assert d[0, :, :].max() < 1e-12 and d[-1, :, :].max() < 1e-12       # the boundary faces too


def test_open_lateral_default_is_wall_and_byte_identical():
    """Opt-in: the default must reproduce the previous (solid-wall) behaviour exactly."""
    from storm_dynamics.pressure_fft import project_anelastic_fft
    a = _poisson_case(); b = _poisson_case()
    project_anelastic_fft(*a[:5], *a[5:], periodic_h=False)
    project_anelastic_fft(*b[:5], *b[5:], periodic_h=False, lateral="wall")
    for x, y in zip(a[:3], b[:3]):
        assert np.array_equal(x, y)


def test_open_lateral_preserves_the_inflow_profile():
    """A solid-wall projection zeroes the boundary normal velocity; 'open' must not -- it may
    only remove a UNIFORM offset (the Neumann solvability condition), keeping the structure."""
    from storm_dynamics.pressure_fft import project_anelastic_fft
    u, v, w, rc, rw, dx, dy, dzc, dzf = _poisson_case()
    inflow_before = u[0, :, :].copy()
    project_anelastic_fft(u, v, w, rc, rw, dx, dy, dzc, dzf, periodic_h=False, lateral="open")
    assert np.abs(u[0, :, :]).max() > 0.0                       # NOT zeroed
    shifted = u[0, :, :] - inflow_before
    assert np.ptp(shifted) < 1e-9                                # differs only by a constant


def test_limited_area_grid_is_flagged_for_the_open_projection():
    """The routing: an outflow, non-periodic config must reach the projection as 'open'."""
    cfg = build_storm_config(preset="storm", nx=16, ny=16, nz=12, Lx=16000.0, Ly=16000.0,
                             Lz=8000.0, duration=1.0, dt_max=3.0, device="cpu", periodic=False)
    cfg.sim.boundaries.x_west = cfg.sim.boundaries.x_east = "outflow"
    cfg.sim.boundaries.y = "outflow"
    assert StormSimulation(cfg).grid._lateral_open is True
    # a nest-style / periodic domain must NOT be flagged
    cfg2 = build_storm_config(preset="storm", nx=16, ny=16, nz=12, Lx=16000.0, Ly=16000.0,
                              Lz=8000.0, duration=1.0, dt_max=3.0, device="cpu", periodic=True)
    assert StormSimulation(cfg2).grid._lateral_open is False


# ------------------------------------- the Davies zone in METRES (audit row 14, last of the class)
def test_davies_zone_cell_width_rescales_but_physical_width_does_not():
    """`width` is a CELL COUNT, so the band's physical width changes with the mesh -- the same
    defect class as NestSpec.relax_width.  It had not corrupted a result only because parent dx
    was constant within every run; a parent-RESOLUTION study is exactly what it breaks."""
    from storm_dynamics.limited_area import lateral_relaxation_weight
    reach_cells, reach_phys = {}, {}
    for nx in (40, 80):                                   # dx = 600 m and 300 m over 24 km
        g = Grid(nx=nx, ny=nx, nz=8, Lx=24000.0, Ly=24000.0, Lz=2000.0, periodic=False)
        mid = nx // 2                                     # slice mid-domain: the corner is all-band
        a = np.asarray(g.backend.to_cpu(lateral_relaxation_weight(g, width=8)))[:, mid, 0]
        b = np.asarray(g.backend.to_cpu(
            lateral_relaxation_weight(g, width_m=4800.0)))[:, mid, 0]
        reach_cells[g.dx] = float((a[:mid] > 0).sum() * g.dx)
        reach_phys[g.dx] = float((b[:mid] > 0).sum() * g.dx)
    # the DEFECT: a cell-count band halves in physical width when the mesh refines
    assert reach_cells[600.0] > 1.8 * reach_cells[300.0], reach_cells
    # the FIX: a metre-specified band is mesh-independent
    assert abs(reach_phys[600.0] - reach_phys[300.0]) <= 600.0 + 1e-9, reach_phys
    assert all(abs(v - 4800.0) <= 600.0 for v in reach_phys.values()), reach_phys


def test_davies_zone_default_is_unchanged():
    """Opt-in: width_m=None must reproduce the cell-count behaviour exactly."""
    from storm_dynamics.limited_area import lateral_relaxation_weight
    g = Grid(nx=40, ny=40, nz=8, Lx=24000.0, Ly=24000.0, Lz=2000.0, periodic=False)
    a = np.asarray(g.backend.to_cpu(lateral_relaxation_weight(g, width=8)))
    b = np.asarray(g.backend.to_cpu(lateral_relaxation_weight(g, width=8, width_m=None)))
    assert np.array_equal(a, b)
