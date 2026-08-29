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
import pytest

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


def test_storm_following_nest_phase2b_sustains():
    """M3 phase 2b: the storm-following (storm-relative) nest keeps the cell centred
    so the updraft is sustained -- not decayed toward zero as a fixed nest would --
    over a window, running stably and conserving."""
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
    nest, rep = nst.run_concurrent_nest(parent, spec, window=180.0, follow=True)
    r = rep["rotation"]; c = rep["conservation"]
    w_hist = [h["w_max"] for h in nest.history]
    assert np.isfinite(r["zeta_abs_max"]) and max(w_hist) < 45.0   # stable, physical
    assert w_hist[-1] > 2.0                                        # updraft sustained (not collapsed)
    assert abs(c["total_water_rel_err"]) < 1.5e-2                  # conserves well when centred
    assert "storm-following" in rep["nest"]["mode"]
    assert len(rep["nest"]["storm_motion"]) == 2                   # C recorded


def test_conservative_restriction_preserves_overlap_integral():
    """The first rigorous AMR conservation piece: average-down (conservative
    restriction) of a cell-aligned, matched-z nest preserves the scalar integral
    over the overlap EXACTLY (to machine precision) -- unlike phase-3a injection."""
    from types import SimpleNamespace
    pg = Grid(nx=24, ny=24, nz=20, Lx=24000, Ly=24000, Lz=12000, periodic=True, z_stretch=1.05)
    spec = nst.NestSpec.aligned(pg, i0=6, j0=6, ncx=8, ncy=8, refine=3)
    ng = nst.build_nest_grid(spec, pg)
    assert ng.nz == pg.nz and np.allclose(np.asarray(ng.zc), np.asarray(pg.zc))  # matched z
    assert ng.nx == 8 * 3 and ng.ny == 8 * 3                                     # aligned
    parent = SimpleNamespace(grid=pg, state=FlowState.zeros(pg))
    nest = SimpleNamespace(grid=ng, state=FlowState.zeros(ng))
    xc = np.asarray(ng.xc).reshape(-1, 1, 1); yc = np.asarray(ng.yc).reshape(1, -1, 1)
    zc = np.asarray(ng.zc).reshape(1, 1, -1); ones = np.ones(ng.center_shape)
    nest.state.theta = (300 + 5 * np.sin(2 * np.pi * xc / ng.Lx) *
                        np.cos(2 * np.pi * yc / ng.Ly) * np.exp(-zc / 6000)) * ones
    nest.state.qv = 1e-3 * (1 + 0.5 * np.cos(3 * np.pi * xc / ng.Lx)) * ones
    res = nst.conservative_restrict(nest, parent, spec)
    for nm, v in res.items():
        assert abs(v["coarse_after_minus_fine"]) < 1e-6, (nm, v)   # exact conservation


def test_amr_refluxing_conserves_across_interface():
    """AMR Milestone 1: Berger-Colella refluxing restores exact conservation across
    a static coarse-fine interface -- WITHOUT it the total mass drifts (interface
    leak); WITH it the drift is machine precision."""
    from storm_dynamics import amr
    d = amr.demo(nsteps=40)
    assert d["no_reflux"] > 1e-6, d          # a real leak exists to fix
    assert d["reflux"] < 1e-12, d            # refluxing conserves to machine precision
    assert d["reflux"] < d["no_reflux"] / 1e6
    assert d["free_stream"] < 1e-12, d       # a uniform field stays uniform (correctness)


def test_amr_port_scaffold_conserves_on_multifab():
    """Port scaffold (M3): the field on an AMReX MultiFab, halo-exchanged by AMReX,
    stepped by our flux-form physics, conserves mass to machine precision.  Skipped
    where pyAMReX is unavailable (it lives in the WSL amr312 env only)."""
    from storm_dynamics import amr_port
    if not amr_port.have_amrex():
        pytest.skip("pyAMReX not installed (build via scripts/build_pyamrex_wsl.sh)")
    d = amr_port.demo(n=24, nsteps=30)
    assert d["mass_rel_drift"] < 1e-12, d


def test_poisson_multigrid_converges_and_is_second_order():
    """Our geometric-multigrid Poisson (the AMR projection kernel pyAMReX doesn't
    provide) converges h-independently and is 2nd-order accurate on a manufactured
    solution phi = sin(2pi x) sin(2pi y)."""
    from storm_dynamics import poisson_mg as mg
    errs, ncyc, last = {}, {}, None
    for n in (64, 128):
        h = 1.0 / n
        xs = (np.arange(n) + 0.5) * h
        X, Y = np.meshgrid(xs, xs, indexing="ij")
        exact = np.sin(2 * np.pi * X) * np.sin(2 * np.pi * Y)
        phi, hist = mg.solve(-8 * np.pi ** 2 * exact, h)
        errs[n] = float(np.abs((phi - phi.mean()) - (exact - exact.mean())).max())
        ncyc[n] = len(hist); last = hist[-1]
    assert last < 1e-8                              # converged
    assert ncyc[64] <= 15 and ncyc[128] <= 15       # h-independent V-cycle count
    ratio = errs[64] / errs[128]
    assert 3.5 < ratio < 4.5, ratio                 # error ∝ h^2 (2nd order)


def test_composite_poisson_1d_second_order_across_interface():
    """The composite 2-level Poisson (coarse periodic grid + fine patch) with the
    2nd-order coarse-fine interface stencil is 2nd-order accurate -- the AMR-
    projection crux, verified in 1-D on a manufactured solution."""
    from storm_dynamics import composite_poisson as cp
    e1, _ = cp.manufactured_error(96)
    e2, _ = cp.manufactured_error(192)
    assert 3.5 < e1 / e2 < 4.5, (e1, e2)            # error ∝ h^2 through the interface
    assert e2 < 2e-4, e2                            # accurate


def test_composite_poisson_2d_second_order_with_corners():
    """The 2-D composite Poisson (fine rectangular patch in a periodic coarse grid,
    with the tangential coarse interpolation and the four corners) is 2nd-order
    accurate through the coarse-fine interface."""
    from storm_dynamics import composite_poisson as cp
    e1, _ = cp.manufactured_error_2d(48)
    e2, _ = cp.manufactured_error_2d(96)
    assert 3.5 < e1 / e2 < 4.5, (e1, e2)            # 2nd order incl. corners
    assert e2 < 1e-3, e2


def test_composite_poisson_3d_second_order_with_edges_and_corners():
    """The 3-D composite Poisson (fine box patch in a periodic coarse grid, with the
    bilinear tangential coarse interpolation and the box edges/corners) is 2nd-order
    accurate through the coarse-fine interface."""
    from storm_dynamics import composite_poisson as cp
    e1, _ = cp.manufactured_error_3d(10)           # doubling 10->20 (direct solve)
    e2, _ = cp.manufactured_error_3d(20)
    assert 3.3 < e1 / e2 < 5.2, (e1, e2)            # ~h^2 across the interface
    assert e2 < 1.5e-2, e2


def test_composite_projection_2d_divergence_free_across_interface():
    """The two-level MAC projection built on the composite Poisson makes a random
    face-flux velocity discretely divergence-free to the solve tolerance, INCLUDING
    at the coarse-fine interface (divergence, gradient and the Laplacian share the
    single-valued interface flux, so div(grad p) = L p exactly)."""
    from storm_dynamics import composite_poisson as cp
    for nc in (12, 24):
        dc, df, di = cp.project_divergence_2d(nc, 2, seed=1)
        assert dc < 1e-9 and df < 1e-9 and di < 1e-9, (nc, dc, df, di)


def test_amr_conservative_prolong_is_inverse_of_restrict():
    """The regridding operator (coarse->fine) is conservative and inverts
    average-down on a constant block: restrict(prolong(x)) == x."""
    import numpy as _np
    from storm_dynamics import amr
    x = _np.random.default_rng(0).random((5, 7))
    f = amr.conservative_prolong(x, 3)
    assert f.shape == (15, 21)
    back = f.reshape(5, 3, 7, 3).mean(axis=(1, 3))
    assert _np.allclose(back, x)             # exact inverse -> conservative


def test_two_way_feedback_phase3a_influences_parent():
    """M3 phase 3a: approximate two-way feedback runs stably and the nest's finer
    solution measurably changes the parent overlap (a closed parent<->nest loop)."""
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    from storm_dynamics import rotation as rot

    def mature():
        scfg = build_storm_config(preset="storm", nx=18, ny=18, nz=32,
                                  Lx=27000, Ly=27000, Lz=15000, duration=720.0,
                                  dt_max=3.0, hodograph_kind="quarter_circle", drag=True,
                                  z_stretch=1.05, U_max=18.0, z_turn=2000.0, C_s=0.22)
        p = StormSimulation(scfg); p.run(); return p

    def _spec(p):
        wc = np.asarray(p.grid.backend.to_cpu(rot._centered_velocity(p.state, p.grid)[2]))
        i, j = np.unravel_index(np.argmax(wc.max(axis=2)), wc.shape[:2])
        return nst.NestSpec.around(p.grid, float(p.grid.xc[i]), float(p.grid.yc[j]),
                                   half=7000.0, refine=3, nz=38)

    p_tw = mature()
    _, rep = nst.run_concurrent_nest(p_tw, _spec(p_tw), window=150.0, follow=True, two_way=True)
    w_tw = float(p_tw.grid.xp.abs(p_tw.state.w).max())
    assert np.isfinite(rep["rotation"]["zeta_abs_max"])            # stable loop, no blow-up
    assert abs(rep["conservation"]["total_water_rel_err"]) < 1.5e-2
    assert rep["nest"]["two_way"] and "two-way" in rep["nest"]["mode"]

    p_ctl = mature()
    nst.run_concurrent_nest(p_ctl, _spec(p_ctl), window=150.0, follow=True, two_way=False)
    w_ctl = float(p_ctl.grid.xp.abs(p_ctl.state.w).max())
    # the feedback measurably changed the parent (same maturation, seed) -> two-way active
    assert abs(w_tw - w_ctl) > 0.3, (w_tw, w_ctl)
