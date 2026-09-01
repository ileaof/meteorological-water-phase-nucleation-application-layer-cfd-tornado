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


def test_composite_project_two_level_call_site_divergence_free():
    """The composite-projection CALL SITE: one solve over the parent(coarse)+nest(fine)
    mass fluxes makes div(rho0 u)=0 across the coarse-fine interface, operating on real
    staggered FlowState arrays with a stretched anelastic density profile."""
    from types import SimpleNamespace
    pg = Grid(nx=18, ny=18, nz=12, Lx=18000, Ly=18000, Lz=12000, periodic=True, z_stretch=1.05)
    spec = nst.NestSpec.aligned(pg, i0=6, j0=6, ncx=6, ncy=6, refine=2)
    ng = nst.build_nest_grid(spec, pg)
    zc = np.asarray(pg.zc); zf = np.asarray(pg.zf)
    rho0_c = np.exp(-zc / 8000.0)
    rho0_wface = np.interp(zf, zc, rho0_c)
    rng = np.random.default_rng(0)
    ps = FlowState.zeros(pg); ns = FlowState.zeros(ng)
    ps.u = rng.standard_normal(pg.u_shape); ps.v = rng.standard_normal(pg.v_shape)
    ps.w = rng.standard_normal(pg.w_shape)
    ns.u = rng.standard_normal(ng.u_shape); ns.v = rng.standard_normal(ng.v_shape)
    ns.w = rng.standard_normal(ng.w_shape)
    parent = SimpleNamespace(grid=pg, state=ps, dynamics="anelastic",
                             rho0_c=rho0_c, rho0_wface=rho0_wface)
    nest = SimpleNamespace(grid=ng, state=ns)
    res = nst.composite_project_two_level(parent, nest, spec)
    assert res["div_coarse"] < 1e-9, res
    assert res["div_fine"] < 1e-9, res
    assert res["div_interface"] < 1e-9, res


def test_composite_projection_in_time_loop():
    """§1 (docs/ROADMAP.md): with ``composite_projection=True`` the stepping driver
    replaces the two per-level projections with ONE composite solve, so
    ``div(rho0 u) = 0`` holds across the interface EVERY sub-step (not per level
    in isolation).  Asserts the worst-case interface divergence over the window
    stays ~solve tolerance, the run stays stable/physical and conserves."""
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    from storm_dynamics import rotation as rot
    scfg = build_storm_config(preset="storm", nx=14, ny=14, nz=16,
                              Lx=24000, Ly=24000, Lz=12000, duration=180.0,
                              dt_max=3.0, hodograph_kind="quarter_circle", drag=True,
                              z_stretch=1.05, U_max=14.0, z_turn=2000.0, C_s=0.22)
    parent = StormSimulation(scfg)
    parent.run()
    wc = np.asarray(parent.grid.backend.to_cpu(rot._centered_velocity(parent.state, parent.grid)[2]))
    i, j = np.unravel_index(np.argmax(wc.max(axis=2)), wc.shape[:2])
    spec = nst.NestSpec.aligned(parent.grid, i0=max(0, i - 3), j0=max(0, j - 3),
                                ncx=6, ncy=6, refine=2)
    nest, rep = nst.run_concurrent_nest(parent, spec, window=18.0,
                                        composite_projection=True)
    r = rep["rotation"]; c = rep["conservation"]
    # definition of done: |div(rho0 u)| at the interface ~ solve tolerance each step
    assert rep["nest"]["composite_div_interface"] < 1e-9, rep["nest"]
    assert rep["n_steps"] > 0 and nest._composite
    assert np.isfinite(r["zeta_abs_max"])                      # no blow-up
    assert max(h["w_max"] for h in nest.history) < 45.0        # physical
    assert abs(c["total_water_rel_err"]) < 3e-2
    assert c["mass_continuity_residual_norm"] < 1e-2
    assert "composite two-level projection" in rep["nest"]["mode"]
    # and it intensifies at least as much as the per-level sponge path
    nest_sponge, _ = nst.run_concurrent_nest(parent, spec, window=18.0)
    assert max(x["w_max"] for x in nest.history) >= max(
        x["w_max"] for x in nest_sponge.history) - 1.0  # at least as strong


def test_adaptive_regridding_primitives_track_the_vortex():
    """ROADMAP §2a: tag_cells + cluster_to_box + regrid_spec detect the rotating region
    and return an aligned nest footprint that contains the vortex (data-driven follow)."""
    from types import SimpleNamespace
    nx = ny = 24; nz = 8
    pg = Grid(nx=nx, ny=ny, nz=nz, Lx=24000.0, Ly=24000.0, Lz=8000.0, periodic=True)
    st = FlowState.zeros(pg)
    ic, jc = 16, 7                                          # vortex centre (coarse column)
    x0 = float(np.asarray(pg.xc)[ic]); y0 = float(np.asarray(pg.yc)[jc]); sig = 2500.0
    Xu, Yu = np.meshgrid(np.asarray(pg.xf), np.asarray(pg.yc), indexing="ij")   # u faces
    Xv, Yv = np.meshgrid(np.asarray(pg.xc), np.asarray(pg.yf), indexing="ij")   # v faces
    Om = 0.02
    wu = np.exp(-((Xu - x0) ** 2 + (Yu - y0) ** 2) / (2 * sig ** 2))
    wv = np.exp(-((Xv - x0) ** 2 + (Yv - y0) ** 2) / (2 * sig ** 2))
    st.u = np.broadcast_to((-Om * (Yu - y0) * wu)[:, :, None], pg.u_shape).copy()
    st.v = np.broadcast_to((+Om * (Xv - x0) * wv)[:, :, None], pg.v_shape).copy()
    parent = SimpleNamespace(state=st, grid=pg)

    tags = nst.tag_cells(st, pg, field="zeta", frac=0.5)
    assert tags[ic, jc], "the vortex column should be tagged"
    i0, j0, ncx, ncy = nst.cluster_to_box(tags, margin=2)
    assert i0 <= ic < i0 + ncx and j0 <= jc < j0 + ncy, "box must contain the vortex"
    assert nst.cluster_to_box(np.zeros((nx, ny), bool)) is None      # nothing tagged -> None

    spec = nst.regrid_spec(parent, refine=3, field="zeta", frac=0.5, margin=2)
    fi0 = int(round(spec.x0 / pg.dx)); fj0 = int(round(spec.y0 / pg.dy))
    fncx = int(round(spec.Lx / pg.dx)); fncy = int(round(spec.Ly / pg.dy))
    assert fi0 <= ic < fi0 + fncx and fj0 <= jc < fj0 + fncy, "aligned footprint must contain it"
    assert spec.refine == 3 and spec.nz == nz               # aligned => matched-z


def test_regrid_nest_preserves_fine_structure_in_overlap():
    """ROADMAP §2a increment 2: regrid_nest re-creates the nest at a shifted (aligned)
    footprint, preserving the old nest's fine field EXACTLY in the overlap (integer
    fine-cell shift) and filling the newly-exposed strip from the parent."""
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    scfg = build_storm_config(preset="storm", nx=18, ny=18, nz=16, Lx=18000.0, Ly=18000.0,
                              Lz=10000.0, duration=1.0, device="cpu")
    parent = StormSimulation(scfg)
    spec = nst.NestSpec.aligned(parent.grid, i0=6, j0=6, ncx=6, ncy=6, refine=3)
    nest = nst.NestedStormSimulation(parent, spec)
    nfx, nfy, nz = nest.grid.nx, nest.grid.ny, nest.grid.nz
    # stamp a distinctive fine pattern (unique per fine column) into theta
    a = np.arange(nfx)[:, None, None]; b = np.arange(nfy)[None, :, None]
    nest.state.theta = (300.0 + 1.0 * a + 0.1 * b) * np.ones((nfx, nfy, nz))
    old_theta = nest.state.theta.copy()

    new_spec = nst.NestSpec.aligned(parent.grid, i0=7, j0=6, ncx=6, ncy=6, refine=3)  # +1 coarse cell in x
    new = nst.regrid_nest(nest, parent, new_spec)
    d0 = (7 - 6) * 3                                        # exact fine-cell shift
    # overlap: new.theta[:nfx-d0] must equal the old field shifted by d0 (structure preserved)
    assert np.allclose(new.state.theta[:nfx - d0], old_theta[d0:], atol=1e-9), "overlap not preserved"
    # the newly-exposed strip came from the parent (horizontally uniform base state,
    # NOT the stamped x-ramp) -> its columns should be ~constant in x, unlike old_theta
    strip = np.asarray(new.state.theta[nfx - d0:])
    assert np.ptp(strip[:, 0, nz // 2]) < 0.5, "exposed strip should be parent-filled (uniform), not the ramp"
    assert int(round(new.spec.x0 / parent.grid.dx)) == 7 and new.grid.nx == nfx


def test_regrid_interval_recentres_nest_on_the_vortex():
    """ROADMAP §2a increment 2 (loop wiring): with regrid_interval, a nest that starts
    off the vortex hops (data-driven, ground frame) to re-centre on the tagged rotation."""
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    scfg = build_storm_config(preset="storm", nx=18, ny=18, nz=12, Lx=18000.0, Ly=18000.0,
                              Lz=9000.0, duration=1.0, dt_max=2.0, device="cpu")
    parent = StormSimulation(scfg); pg = parent.grid
    ic, jc = 14, 14                                        # vortex far from the nest's start
    x0 = float(np.asarray(pg.xc)[ic]); y0 = float(np.asarray(pg.yc)[jc]); sig = 2500.0
    Xu, Yu = np.meshgrid(np.asarray(pg.xf), np.asarray(pg.yc), indexing="ij")
    Xv, Yv = np.meshgrid(np.asarray(pg.xc), np.asarray(pg.yf), indexing="ij")
    Om = 0.03
    wu = np.exp(-((Xu - x0) ** 2 + (Yu - y0) ** 2) / (2 * sig ** 2))
    wv = np.exp(-((Xv - x0) ** 2 + (Yv - y0) ** 2) / (2 * sig ** 2))
    parent.state.u = parent.state.u + np.broadcast_to((-Om * (Yu - y0) * wu)[:, :, None], pg.u_shape)
    parent.state.v = parent.state.v + np.broadcast_to((+Om * (Xv - x0) * wv)[:, :, None], pg.v_shape)
    spec = nst.NestSpec.aligned(pg, i0=1, j0=1, ncx=6, ncy=6, refine=3)     # start in a corner
    nest, rep = nst.run_concurrent_nest(parent, spec, window=2.5, regrid_interval=1,
                                        regrid_field="zeta", regrid_frac=0.5)
    assert rep["nest"]["regrids"] >= 1, rep["nest"]
    fi0, fj0 = rep["nest"]["final_footprint"]
    assert fi0 >= 6 and fj0 >= 6, (fi0, fj0)              # hopped toward the vortex at (14,14)
    assert np.isfinite(float(np.max(np.abs(nest.state.w))))


def test_fft_tridiag_anelastic_projector_divergence_free():
    """ROADMAP §3f — the low-memory fine-nest pressure solver. The FFT/DCT-in-(x,y) +
    tridiagonal-in-z projector makes div(rho0 u)=0 to round-off on a STRETCHED grid, for
    both the parent (periodic x,y) and the nest (wall x,y), with no stored LU factorisation
    (the direct splu that OOMs at ~48³ is replaced by O(N) transforms + Thomas)."""
    from storm_dynamics.pressure_fft import project_anelastic_fft, anelastic_divergence
    nx = 24; nz = 40; dx = dy = 500.0
    dz = 1.05 ** np.arange(nz); dz *= 15000.0 / dz.sum()
    zf = np.concatenate([[0.0], np.cumsum(dz)]); zc = 0.5 * (zf[:-1] + zf[1:])
    dzc = zf[1:] - zf[:-1]; dzf = np.diff(zc)
    rc = np.exp(-zc / 8000.0); rw = np.interp(zf, zc, rc)
    for periodic_h in (True, False):
        rng = np.random.default_rng(0)
        u = rng.standard_normal((nx + 1, nx, nz)); v = rng.standard_normal((nx, nx + 1, nz))
        w = rng.standard_normal((nx, nx, nz + 1))
        before = np.abs(anelastic_divergence(u, v, w, rc, rw, dx, dy, dzc)).max()
        res = project_anelastic_fft(u, v, w, rc, rw, dx, dy, dzc, dzf, periodic_h=periodic_h)
        assert before > 1e-3 and res < 1e-10, (periodic_h, before, res)   # div driven to round-off


def test_iterative_cg_anelastic_projector_matches_fft_and_is_divergence_free():
    """ROADMAP §3f increment 3 — the GENERAL (non-separable) low-memory solver.  The
    Jacobi-preconditioned CG on the SPD (negated, volume-weighted) anelastic operator makes
    div(rho0 u)=0 to round-off on a STRETCHED grid, periodic AND wall, with NO factorisation
    (splu OOMs, ILU breaks CG's SPD requirement, the wrong-sign operator makes CG stall) — and
    it AGREES with the independent FFT+tridiag projector, so the two cross-validate."""
    from storm_dynamics.pressure_iterative import project_anelastic_iterative
    from storm_dynamics.pressure_fft import project_anelastic_fft, anelastic_divergence
    nx = 24; nz = 40; dx = dy = 500.0
    dz = 1.05 ** np.arange(nz); dz *= 15000.0 / dz.sum()
    zf = np.concatenate([[0.0], np.cumsum(dz)]); zc = 0.5 * (zf[:-1] + zf[1:])
    dzc = zf[1:] - zf[:-1]; dzf = np.diff(zc)
    rc = np.exp(-zc / 8000.0); rw = np.interp(zf, zc, rc)
    for periodic_h in (True, False):
        rng = np.random.default_rng(0)
        u = rng.standard_normal((nx + 1, nx, nz)); v = rng.standard_normal((nx, nx + 1, nz))
        w = rng.standard_normal((nx, nx, nz + 1))
        uf, vf, wf = u.copy(), v.copy(), w.copy()
        before = np.abs(anelastic_divergence(u, v, w, rc, rw, dx, dy, dzc)).max()
        res = project_anelastic_iterative(u, v, w, rc, rw, dx, dy, dzc, dzf, periodic_h=periodic_h)
        assert before > 1e-3 and res < 1e-8, (periodic_h, before, res)    # div driven to the CG tol
        project_anelastic_fft(uf, vf, wf, rc, rw, dx, dy, dzc, dzf, periodic_h=periodic_h)
        du = max(np.abs(u - uf).max(), np.abs(v - vf).max(), np.abs(w - wf).max())
        assert du < 1e-6, (periodic_h, du)                                # two solvers agree


def test_lowmem_pressure_router_matches_direct_and_is_divergence_free():
    """ROADMAP §3f increment 2 — the router `_project_anelastic_lowmem` picks the low-memory
    solver by grid structure (`separable` -> FFT+tridiag, else Jacobi-CG) and is physically
    transparent: on a properly periodic / zero-wall-normal state (what the BCs maintain) it
    drives div(rho0 u) to round-off AND lands on the same interior field the exact direct
    solve does, up to the discrete null-space gauge.  So switching large nests onto it does
    not change the physics."""
    from meteorological_flow.pressure_solver import PressureSolver
    from storm_dynamics.core import _project_anelastic_lowmem, separable
    for periodic in (True, False):
        g = Grid(nx=16, ny=16, nz=20, Lx=8000.0, Ly=8000.0, Lz=12000.0,
                 z_stretch=1.05, periodic=periodic)
        assert separable(g)                                    # storm grids are separable -> FFT
        rc = np.exp(-np.asarray(g.zc) / 8000.0); rw = np.interp(np.asarray(g.zf), np.asarray(g.zc), rc)
        rng = np.random.default_rng(0); st = FlowState.zeros(g)
        st.u = rng.standard_normal(g.u_shape); st.v = rng.standard_normal(g.v_shape)
        st.w = rng.standard_normal(g.w_shape); st.w[:, :, 0] = 0.0; st.w[:, :, -1] = 0.0
        if periodic:
            st.u[-1] = st.u[0]; st.v[:, -1] = st.v[:, 0]       # proper periodic wrap faces
        else:
            st.u[0] = st.u[-1] = 0.0; st.v[:, 0] = st.v[:, -1] = 0.0   # solid wall normals
        sd, sl = st.copy(), st.copy()
        PressureSolver(g, method="direct").project_anelastic(sd, 1.0, rc, rw)
        res = _project_anelastic_lowmem(sl, g, rc, rw)
        assert res < 1e-10, (periodic, res)                    # div driven to round-off
        iu = slice(0, -1) if periodic else slice(1, -1)        # exclude redundant/wall faces
        du = np.abs(np.asarray(sd.u)[iu] - np.asarray(sl.u)[iu]).max()
        assert du < 1e-4, (periodic, du)                       # matches direct in the interior


def test_lowmem_pressure_path_steps_stably_in_the_time_loop():
    """ROADMAP §3f increment 2 — the wired `_project` low-memory branch (forced here on a small
    grid so it is cheap) steps the full storm stably: the projection keeps the mass-continuity
    residual near round-off and the rotation diagnostics finite, exactly as the direct solve."""
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    scfg = build_storm_config(preset="storm", nx=16, ny=16, nz=20, Lx=8000.0, Ly=8000.0,
                              Lz=12000.0, duration=1.0, device="cpu")
    sim = StormSimulation(scfg)
    assert sim.dynamics == "anelastic"
    sim._lowmem_pressure = True                                # force the §3f route (small grid)
    sim.cfg.time.duration = 3 * float(sim._dt())
    rep = sim.run()
    assert np.isfinite(rep["rotation"]["zeta_abs_max"])
    # residual is measured after the full predictor->project->transport step (not right after
    # the projection), so it is small but not round-off -- same order the direct solve gives.
    assert rep["conservation"]["mass_continuity_residual_norm"] < 1e-3, rep["conservation"]


def test_pressure_cg_does_not_converge_on_stretched_operator():
    """ROADMAP §2b/§3f — documents *why* the fine-nest memory fix is not just "switch to CG".
    The direct (splu) solve is exact on the stretched anelastic Poisson; the existing
    Jacobi-preconditioned **CG does NOT converge** on it, so `_pressure_method` must keep
    stretched grids on `direct` (and the real low-memory fix is a proper solver — FFT/DST +
    tridiagonal-in-z, or ILU/multigrid-preconditioned CG).  If CG is ever fixed, this test
    fails and should be updated to the matches-direct assertion."""
    import warnings
    from meteorological_flow.pressure_solver import PressureSolver
    from storm_dynamics.core import _pressure_method
    g = Grid(nx=12, ny=12, nz=16, Lx=12000.0, Ly=12000.0, Lz=10000.0, periodic=False, z_stretch=1.06)
    assert _pressure_method(g) == "direct"                        # stretched -> direct (never CG)
    rng = np.random.default_rng(0)
    st = FlowState.zeros(g)
    st.u = rng.standard_normal(g.u_shape); st.v = rng.standard_normal(g.v_shape)
    st.w = rng.standard_normal(g.w_shape); st.w[:, :, 0] = 0.0; st.w[:, :, -1] = 0.0
    sd, sc = st.copy(), st.copy()
    rho0_c = np.exp(-np.asarray(g.zc) / 8000.0)
    rho0_wf = np.interp(np.asarray(g.zf), np.asarray(g.zc), rho0_c)
    resd, _ = PressureSolver(g, method="direct").project_anelastic(sd, 1.0, rho0_c, rho0_wf)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rescg, _ = PressureSolver(g, method="cg", tol=1e-10, maxiter=5000).project_anelastic(sc, 1.0, rho0_c, rho0_wf)
    assert resd < 1e-6, resd                                      # direct is exact
    du = np.abs(np.asarray(g.backend.to_cpu(sd.u)) - np.asarray(g.backend.to_cpu(sc.u))).max()
    assert not (du < 1e-3), "CG unexpectedly matches direct -- update _pressure_method + this test"


def test_recursive_nesting_second_level_refines_further():
    """ROADMAP §2b (increment 1 — the funnel-resolution path): a nest OF a nest composes
    (`NestedStormSimulation(nest1, spec2)`), giving Δx = parent.dx / r² and running stably.
    Two refinements take 1.3 km-class parents to O(100 m), the tornado-funnel scale."""
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    scfg = build_storm_config(preset="storm", nx=18, ny=18, nz=16, Lx=18000.0, Ly=18000.0,
                              Lz=10000.0, duration=1.0, device="cpu")
    parent = StormSimulation(scfg)
    spec1 = nst.NestSpec.aligned(parent.grid, i0=6, j0=6, ncx=6, ncy=6, refine=3)
    nest1 = nst.NestedStormSimulation(parent, spec1)
    assert abs(nest1.grid.dx - parent.grid.dx / 3) < 1e-6 * parent.grid.dx   # aligned => exact /3
    spec2 = nst.NestSpec.around(nest1.grid, xc=nest1.grid.Lx / 2, yc=nest1.grid.Ly / 2,
                                half=nest1.grid.Lx * 0.3, refine=3, nz=nest1.grid.nz,
                                z_stretch=parent.grid.z_stretch)
    nest2 = nst.NestedStormSimulation(nest1, spec2)                 # nest of a nest
    # ~parent.dx/9 (integer-nx rounding on the small sub-nest gives a few % slack)
    assert abs(nest2.grid.dx - parent.grid.dx / 9) < 0.15 * (parent.grid.dx / 9)
    assert nest2.grid.dx < nest1.grid.dx < parent.grid.dx           # each level strictly finer
    assert nest2.grid.nz == parent.grid.nz                          # matched z through the stack
    nest2.cfg.time.duration = 3 * float(nest2._dt())
    rep = nest2.run()
    assert np.isfinite(rep["rotation"]["zeta_abs_max"])             # stable / finite


def test_restrict_velocity_face_average_down():
    """ROADMAP §2b: restrict_velocity conservatively averages the fine staggered velocity
    onto the coincident coarse faces of the overlap (the mass-flux-preserving momentum
    feedback up), coarse x-face = block-mean of the r fine faces it tiles."""
    from types import SimpleNamespace
    r = 3
    pg = Grid(nx=12, ny=12, nz=6, Lx=12000.0, Ly=12000.0, Lz=6000.0, periodic=True)
    spec = nst.NestSpec.aligned(pg, i0=4, j0=4, ncx=4, ncy=4, refine=r)
    ng = nst.build_nest_grid(spec, pg)
    parent = SimpleNamespace(grid=pg, state=FlowState.zeros(pg))
    nest = SimpleNamespace(grid=ng, state=FlowState.zeros(ng))
    # stamp u = (fine y-index): a coarse x-face is the block-mean of the r fine y-faces
    yidx = np.arange(ng.ny)[None, :, None].astype(float)
    nest.state.u = np.broadcast_to(yidx, ng.u_shape).copy()
    chg = nst.restrict_velocity(nest, parent, spec)
    cu = np.asarray(pg.backend.to_cpu(parent.state.u))
    for B in range(4):
        assert np.allclose(cu[4, 4 + B, :], B * r + (r - 1) / 2.0), (B, cu[4, 4 + B, 0])
    assert chg > 0.0


def test_multilevel_concurrent_driver_sustains_finest_level():
    """ROADMAP §2b increment 2: run_multilevel_nest drives a 2-level stack — each level
    sub-cycles under the one above (boundaries fed down, conservative scalar restriction
    up) — so the finest level (Δx = parent.dx / r²) is integrated over the window, stable."""
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    scfg = build_storm_config(preset="storm", nx=18, ny=18, nz=16, Lx=18000.0, Ly=18000.0,
                              Lz=10000.0, duration=1.0, dt_max=2.0, device="cpu")
    parent = StormSimulation(scfg)
    spec1 = nst.NestSpec.aligned(parent.grid, i0=6, j0=6, ncx=6, ncy=6, refine=3)
    spec2 = lambda g: nst.NestSpec.around(g, g.Lx / 2, g.Ly / 2, half=g.Lx * 0.3, refine=3,
                                          nz=g.nz, z_stretch=parent.grid.z_stretch)
    sims, rep = nst.run_multilevel_nest(parent, [spec1, spec2], window=4.0)
    assert len(sims) == 3                                       # parent + 2 nest levels
    assert rep["nest"]["levels"] == 2 and rep["nest"]["total_refine"] == 9
    assert abs(sims[-1].grid.dx - parent.grid.dx / 9) < 0.15 * (parent.grid.dx / 9)
    assert sims[-1].step > 0                                    # the finest level was sub-cycled
    assert np.isfinite(rep["rotation"]["zeta_abs_max"])         # stable / finite


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


def test_amr_momentum_reflux_conserves_nonlinear_flux():
    """ROADMAP §2b — the flux-register reflux extended to **momentum** (a NONLINEAR flux
    u²/2, inviscid Burgers): the coarse interface flux differs from the average of the fine
    fluxes across a resolved gradient, so WITHOUT reflux the total momentum ∑u drifts; WITH
    reflux the interface mismatch is applied back to the coarse cell outside the patch and
    momentum is conserved to machine precision.  This is the momentum-conservation algorithm
    the storm's staggered interface reflux ports onto a FluxRegister (`TwoLevelBurgersReflux`),
    nailed down in pure NumPy first -- and the rigorous complement to `restrict_velocity`,
    which only averages the overlap, not the interface flux."""
    from storm_dynamics.amr import TwoLevelBurgersReflux
    no = TwoLevelBurgersReflux().run(nsteps=40, reflux=False)
    yes = TwoLevelBurgersReflux().run(nsteps=40, reflux=True)
    fs = TwoLevelBurgersReflux().free_stream_error(nsteps=30, reflux=True)
    assert no > 1e-6, no                     # a real nonlinear interface leak exists
    assert yes < 1e-12, yes                  # refluxing conserves momentum to machine precision
    assert yes < no / 1e6, (yes, no)
    assert fs < 1e-12, fs                    # uniform momentum stays uniform


def test_amr_momentum_reflux_ports_storm_flux_and_conserves():
    """ROADMAP §2b — the momentum reflux **ported onto the storm's staggered flux**.  Viewed on
    the u-point grid, the storm's x-momentum self-advection is a conservation law with the flux
    `Uc·u` (`Uc=½(uᵢ+uᵢ₊₁)`); `amr._momentum_upwind_x` uses that exact flux (tied here to
    `momentum._u_tendency`'s `Fx_u` on a v=w=0 slab, bit-for-bit), and the flux-register reflux
    conserves the domain x-momentum ∑u across a coarse-fine interface to machine precision where
    it otherwise leaks -- the correction `restrict_velocity` (overlap average-down) omits."""
    from storm_dynamics import amr, momentum as mom
    # (1) conservation of the ported reflux
    no = amr.TwoLevelMomentumReflux().run(nsteps=40, reflux=False)
    yes = amr.TwoLevelMomentumReflux().run(nsteps=40, reflux=True)
    assert no > 1e-6, no                     # a real interface momentum leak exists
    assert yes < 1e-12, yes                  # refluxing conserves ∑u to machine precision
    assert yes < no / 1e6, (yes, no)
    # (2) tie the reference flux to the REAL storm operator (identical, v=w=0 slab)
    g = Grid(nx=24, ny=4, nz=1, Lx=24.0, Ly=4.0, Lz=1.0, periodic=True)
    prof = 1.0 + 0.3 * np.sin(2 * np.pi * np.arange(24) / 24)
    u = np.zeros((25, 4, 1)); u[:24, :, 0] = prof[:, None]; u[24] = u[0]
    st = FlowState.zeros(g); st.u = u; st.v = np.zeros(g.v_shape); st.w = np.zeros(g.w_shape)
    _, _, _, fx = mom.momentum_advection_tendency(st, g, order=1, periodic=True, return_fluxes=True)
    _, Fx_mine = amr._momentum_upwind_x(prof[:, None] * np.ones((1, 4)), 0.0, g.dx, periodic=True)
    assert np.abs(Fx_mine[1:] - fx["Fx_u"][:, :, 0]).max() < 1e-14   # bit-for-bit storm flux


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


def test_composite_solid_wall_bc_second_order_and_projection():
    """The solid-wall (Neumann) composite solve is 2nd-order, and the two-level
    projection is divergence-free across the interface with walls too -- the storm's
    boundary condition (plumbing step-2 item a)."""
    from storm_dynamics import composite_poisson as cp
    e1, _ = cp.manufactured_error_2d_wall(48)
    e2, _ = cp.manufactured_error_2d_wall(96)
    assert 3.5 < e1 / e2 < 4.5, (e1, e2)            # 2nd order with walls
    for nc in (12, 24):
        dc, df, di = cp.project_divergence_2d(nc, 2, seed=2, periodic=False)
        assert dc < 1e-9 and df < 1e-9 and di < 1e-9, (nc, dc, df, di)


def test_composite_massflux_bridge_storm_arrays_divergence_free():
    """Step-2 item (b): the face-array bridge reads the storm's staggered C-grid
    mass-flux arrays (u:(n+1,n), v:(n,n+1)) for coarse parent + fine nest, projects,
    and writes back so that div(m) recomputed straight from the written-back arrays is
    ~0 across the interface -- with the solid-wall BC the storm uses (and periodic)."""
    import numpy as np
    from storm_dynamics import composite_poisson as cp
    for periodic in (True, False):
        for nc in (12, 24):
            r = 2; ci0, ci1, cj0, cj1 = nc // 3, 2 * nc // 3, nc // 3, 2 * nc // 3
            nfx, nfy = r * (ci1 - ci0), r * (cj1 - cj0)
            rng = np.random.default_rng(3)
            mu_c = rng.standard_normal((nc + 1, nc)); mv_c = rng.standard_normal((nc, nc + 1))
            mu_f = rng.standard_normal((nfx + 1, nfy)); mv_f = rng.standard_normal((nfx, nfy + 1))
            dc, df, di = cp.composite_project_massflux_2d(
                mu_c, mv_c, mu_f, mv_f, nc, r, ci0, ci1, cj0, cj1, periodic=periodic)
            assert dc < 1e-9 and df < 1e-9 and di < 1e-9, (periodic, nc, dc, df, di)


def test_composite_stretched_vertical_metric_second_order():
    """Step-2 item (c): the stretched vertical metric (variable-dz finite volume, walls)
    composes with the horizontal composite interface. Uniform z is clean 2nd order;
    moderate stretching is supraconvergent (~1.8-2), the standard non-uniform-FV order."""
    from storm_dynamics import composite_poisson as cp
    e1 = cp.manufactured_error_metric_z(24, 24, s=1.0)      # uniform z -> clean 2nd order
    e2 = cp.manufactured_error_metric_z(48, 48, s=1.0)
    assert 3.6 < e1 / e2 < 4.4, (e1, e2)
    s1 = cp.manufactured_error_metric_z(24, 24, s=1.05)     # stretched -> supraconvergent
    s2 = cp.manufactured_error_metric_z(48, 48, s=1.05)
    assert 3.0 < s1 / s2 < 4.4 and s2 < 1e-2, (s1, s2)


def test_composite_hz_unified_operator_second_order():
    """FINAL ASSEMBLY: the unified 3-D operator (horizontal composite interface at each
    z-level + variable-dz vertical, walls) is 2nd-order on a manufactured solution."""
    from storm_dynamics import composite_poisson as cp
    e1 = cp.manufactured_error_hz(12, 12, s=1.0)
    e2 = cp.manufactured_error_hz(24, 24, s=1.0)
    assert 3.6 < e1 / e2 < 4.4 and e2 < 5e-2, (e1, e2)


def test_composite_projection_hz_full_storm_divergence_free():
    """FINAL ASSEMBLY: the full 3-D storm mass-flux projection (horizontal composite
    interface + variable-dz vertical + walls) makes div(m)=0 across the coarse-fine
    interface, recomputed independently from the written-back staggered arrays, for both
    the nest (walls) and the parent (periodic horizontal)."""
    import numpy as np
    from storm_dynamics import composite_poisson as cp
    for periodic_h in (False, True):
        nc, nz, r = 12, 8, 2
        ci0, ci1, cj0, cj1 = nc // 3, 2 * nc // 3, nc // 3, 2 * nc // 3
        nfx, nfy = r * (ci1 - ci0), r * (cj1 - cj0)
        dz = 1.05 ** np.arange(nz); dz *= 2.0 / dz.sum()
        zf = np.concatenate([[0.0], np.cumsum(dz)]); zc = 0.5 * (zf[:-1] + zf[1:])
        dzc = zf[1:] - zf[:-1]; dzf = np.diff(zc)
        rng = np.random.default_rng(4)
        mu_c = rng.standard_normal((nc + 1, nc, nz)); mv_c = rng.standard_normal((nc, nc + 1, nz))
        mw_c = rng.standard_normal((nc, nc, nz + 1))
        mu_f = rng.standard_normal((nfx + 1, nfy, nz)); mv_f = rng.standard_normal((nfx, nfy + 1, nz))
        mw_f = rng.standard_normal((nfx, nfy, nz + 1))
        dc, df, di = cp.composite_project_massflux_hz(
            mu_c, mv_c, mw_c, mu_f, mv_f, mw_f, nc, nz, r, ci0, ci1, cj0, cj1,
            dzc, dzf, periodic_h=periodic_h)
        assert dc < 1e-9 and df < 1e-9 and di < 1e-9, (periodic_h, dc, df, di)


def test_composite_hz_anisotropic_second_order_and_projection():
    """Anisotropic (dx != dy, ncx != ncy): the unified operator is still 2nd-order and the
    full 3-D mass-flux projection is divergence-free across the interface -- rectangular /
    non-square parents (removes the square-grid restriction)."""
    import numpy as np
    from storm_dynamics import composite_poisson as cp
    e1 = cp.manufactured_error_hz(12, 12, ncy=18)          # hx=1/12 != hy=1/18
    e2 = cp.manufactured_error_hz(24, 24, ncy=36)
    assert 3.5 < e1 / e2 < 4.5 and e2 < 5e-2, (e1, e2)     # 2nd order, anisotropic
    ncx, ncy, nz, r = 16, 10, 8, 2
    ci0, ci1, cj0, cj1 = 5, 11, 3, 7
    nfx, nfy = r * (ci1 - ci0), r * (cj1 - cj0)
    dz = 1.05 ** np.arange(nz); dz *= 12000.0 / dz.sum()
    zf = np.concatenate([[0.0], np.cumsum(dz)]); zc = 0.5 * (zf[:-1] + zf[1:])
    dzc = zf[1:] - zf[:-1]; dzf = np.diff(zc)
    rng = np.random.default_rng(5)
    mu_c = rng.standard_normal((ncx + 1, ncy, nz)); mv_c = rng.standard_normal((ncx, ncy + 1, nz))
    mw_c = rng.standard_normal((ncx, ncy, nz + 1))
    mu_f = rng.standard_normal((nfx + 1, nfy, nz)); mv_f = rng.standard_normal((nfx, nfy + 1, nz))
    mw_f = rng.standard_normal((nfx, nfy, nz + 1))
    for periodic_h in (False, True):
        dc, df, di = cp.composite_project_massflux_hz(
            mu_c.copy(), mv_c.copy(), mw_c.copy(), mu_f.copy(), mv_f.copy(), mw_f.copy(),
            ncx, nz, r, ci0, ci1, cj0, cj1, dzc, dzf, periodic_h=periodic_h,
            hx=400.0, hy=700.0, ncy=ncy)
        assert dc < 1e-9 and df < 1e-9 and di < 1e-9, (periodic_h, dc, df, di)


def test_composite_projection_3d_divergence_free_across_interface():
    """The 3-D two-level MAC projection (built on solve_3d) makes a random face-flux
    velocity discretely divergence-free including at the coarse-fine interface."""
    from storm_dynamics import composite_poisson as cp
    for nc in (8, 12):
        dc, df, di = cp.project_divergence_3d(nc, 2, seed=1)
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
