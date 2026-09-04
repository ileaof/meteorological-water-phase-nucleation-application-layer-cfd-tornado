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
from storm_dynamics.core import StormSimulation
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
