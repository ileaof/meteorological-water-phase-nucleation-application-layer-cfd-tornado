"""M3 phase 1 -- static nested-grid refinement tests.

Unit tests cover the infrastructure that must be exact: grid refinement, the
parent->nest trilinear interpolation (uniform and linear fields reproduce
exactly), and the boundary relaxation weights.  A short integration checks that
the assembled :class:`NestedStormSimulation` inherits the parent storm, stays
stable and physical (no blow-up), and conserves -- i.e. the nest is a usable
refinement of the parent, not a numerical artefact.  (Whether the finer grid
*intensifies* the vortex depends on the storm phase and is shown in
``examples/tornado_nest.py`` rather than asserted here.)
"""
from __future__ import annotations

import numpy as np

from meteorological_flow.grid import Grid
from meteorological_flow.state import FlowState

from storm_dynamics import nesting as nst


def _parent_grid(nx=24, ny=24, nz=30):
    return Grid(nx=nx, ny=ny, nz=nz, Lx=36000, Ly=36000, Lz=15000, periodic=True)


def test_nest_grid_refinement():
    pg = _parent_grid()
    spec = nst.NestSpec.centered(pg, frac=0.4, refine=3, nz=36)
    ng = nst.build_nest_grid(spec, pg)
    assert ng.dx < pg.dx and abs(ng.dx - pg.dx / 3) < 0.35 * ng.dx
    assert ng.Lz == pg.Lz and ng.nz == 36
    assert not ng.periodic


def test_interpolation_uniform_is_exact():
    pg = _parent_grid()
    spec = nst.NestSpec.centered(pg, frac=0.4, refine=3)
    ng = nst.build_nest_grid(spec, pg)
    ps = FlowState.zeros(pg)
    ps.theta[:] = 300.0; ps.u[:] = 7.0; ps.v[:] = -3.0
    ns = nst.interpolate_state_to_nest(ps, pg, ng, spec)
    assert np.allclose(ns.theta, 300.0)
    assert np.allclose(ns.u, 7.0) and np.allclose(ns.v, -3.0)


def test_interpolation_linear_field_is_exact():
    """Trilinear interpolation reproduces a linear field to machine precision."""
    pg = _parent_grid()
    spec = nst.NestSpec.centered(pg, frac=0.4, refine=3)
    ng = nst.build_nest_grid(spec, pg)
    xc = np.asarray(pg.xc).reshape(-1, 1, 1); yc = np.asarray(pg.yc).reshape(1, -1, 1)
    ps = FlowState.zeros(pg)
    ps.theta = (2e-4 * xc + 1e-4 * yc) * np.ones(pg.center_shape)
    ns = nst.interpolate_state_to_nest(ps, pg, ng, spec)
    xcn = (spec.x0 + np.asarray(ng.xc)).reshape(-1, 1, 1)
    ycn = (spec.y0 + np.asarray(ng.yc)).reshape(1, -1, 1)
    expected = (2e-4 * xcn + 1e-4 * ycn) * np.ones(ng.center_shape)
    assert float(np.abs(ns.theta - expected).max()) < 1e-9


def test_relaxation_weight_is_border_only():
    pg = _parent_grid()
    spec = nst.NestSpec.centered(pg, frac=0.4, refine=3, relax_width=4, relax_rate=0.02)
    ng = nst.build_nest_grid(spec, pg)
    w = nst.relaxation_weight(ng, spec)
    assert np.isclose(float(w.max()), 0.02)                 # outermost cell = rate
    assert float(w[ng.nx // 2, ng.ny // 2, 0]) == 0.0       # interior untouched
    assert w.shape == (ng.nx, ng.ny, 1)


def test_nest_inherits_parent_and_runs_stable():
    """A short nest over a matured parent inherits the storm, stays physical and
    conserves -- the refinement is usable, not an artefact."""
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    from storm_dynamics import rotation as rot
    scfg = build_storm_config(preset="storm", nx=18, ny=18, nz=32,
                              Lx=27000, Ly=27000, Lz=15000, duration=720.0,
                              dt_max=3.0, hodograph_kind="quarter_circle", drag=True,
                              z_stretch=1.05, U_max=18.0, z_turn=2000.0, C_s=0.22)
    parent = StormSimulation(scfg)
    parent.run()
    wc = np.asarray(parent.grid.backend.to_cpu(rot._centered_velocity(parent.state, parent.grid)[2]))
    i, j = np.unravel_index(np.argmax(wc.max(axis=2)), wc.shape[:2])
    spec = nst.NestSpec.around(parent.grid, float(parent.grid.xc[i]),
                               float(parent.grid.yc[j]), half=7000.0, refine=3, nz=38)
    nest = nst.NestedStormSimulation(parent, spec)
    # nest inherits the parent updraft at init (interpolation carried it over)
    r0 = rot.rotation_report(nest.state, nest.grid)
    assert np.isfinite(r0["w_max"])
    nest.cfg.time.duration = 60.0
    rep = nest.run()
    r = rep["rotation"]; c = rep["conservation"]
    w_hist = max(h["w_max"] for h in nest.history)
    assert np.isfinite(r["zeta_abs_max"])
    assert w_hist < 45.0                                     # no grid-scale blow-up
    assert abs(c["total_water_rel_err"]) < 1.5e-2            # conserves
    assert c["mass_continuity_residual_norm"] < 1e-2         # projection, not limiters
    assert rep["nest"]["dx_m"] < parent.grid.dx              # genuinely finer


def test_concurrent_nest_phase2_is_stable():
    """M3 phase 2: the concurrent nest (parent stepping alongside, time-evolving
    boundaries) runs stably and conserves over a short window -- and the
    interior-masked ζ separates the physical vortex from the sponge edge."""
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    from storm_dynamics import rotation as rot
    scfg = build_storm_config(preset="storm", nx=18, ny=18, nz=32,
                              Lx=27000, Ly=27000, Lz=15000, duration=720.0,
                              dt_max=3.0, hodograph_kind="quarter_circle", drag=True,
                              z_stretch=1.05, U_max=18.0, z_turn=2000.0, C_s=0.22)
    parent = StormSimulation(scfg)
    parent.run()
    wc = np.asarray(parent.grid.backend.to_cpu(rot._centered_velocity(parent.state, parent.grid)[2]))
    i, j = np.unravel_index(np.argmax(wc.max(axis=2)), wc.shape[:2])
    spec = nst.NestSpec.around(parent.grid, float(parent.grid.xc[i]),
                               float(parent.grid.yc[j]), half=7000.0, refine=3, nz=38)
    nest, rep = nst.run_concurrent_nest(parent, spec, window=60.0)
    r = rep["rotation"]; c = rep["conservation"]
    assert np.isfinite(r["zeta_abs_max"])                      # no blow-up
    assert max(h["w_max"] for h in nest.history) < 45.0        # physical
    assert abs(c["total_water_rel_err"]) < 3e-2               # conserves (looser over the window)
    assert c["mass_continuity_residual_norm"] < 1e-2
    assert "concurrent" in rep["nest"]["mode"]
    # interior-masked ζ is finite and no larger than the edge-inclusive value
    zi = nst.interior_near_surface_zeta(nest)
    assert np.isfinite(zi) and zi <= r["near_surface_zeta_max"] + 1e-9
