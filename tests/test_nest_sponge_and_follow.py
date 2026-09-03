"""Nest boundary sponge width and moving-nest tracker robustness.

Two defects found on the 4.0 km matched-domain resolution pair (67 m vs 22 m):

1. ``relax_width`` is a CELL COUNT, so the sponge's physical width shrinks by ``refine``
   at every cascade level (4 cells = 267 m at dx=67 m but only 89 m at dx=22 m).  The
   damping crossing the band goes like ``relax_rate * width / U``, so the fine nest damped
   boundary error ~3x less: measured edge/interior |zeta| was 0.9 at 67 m but 5.0 at 22 m.
2. ``tag_cells`` thresholds at ``frac * domain max``, so when the tracked grid is itself a
   nest its sponge-generated edge vorticity both inflates the normaliser and drags the
   tagged cluster into the wall -- the moving nest chases its own artifact.

Both fixes are OPT-IN; the defaults must stay byte-identical.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from storm_dynamics.config import build_storm_config
from storm_dynamics.core import Grid, StormSimulation
from storm_dynamics import nesting as nst


def _grid(nx, L=4000.0):
    return Grid(nx=nx, ny=nx, nz=8, Lx=L, Ly=L, Lz=2000.0, periodic=False)


# --------------------------------------------------------------------------- sponge width
def test_effective_relax_width_defaults_to_the_cell_count():
    """relax_width_m=None -> exactly the historical cell-count behaviour."""
    for nx in (60, 180):
        g = _grid(nx)
        spec = nst.NestSpec(x0=0.0, y0=0.0, Lx=4000.0, Ly=4000.0)
        assert nst.effective_relax_width(spec, g) == spec.relax_width == 4


def test_cell_count_sponge_shrinks_physically_but_metres_do_not():
    """The defect and the fix, side by side, at the two meshes actually compared."""
    coarse, fine = _grid(60), _grid(180)          # dx = 66.7 m and 22.2 m over the SAME 4.0 km
    old = nst.NestSpec(x0=0.0, y0=0.0, Lx=4000.0, Ly=4000.0)
    new = nst.NestSpec(x0=0.0, y0=0.0, Lx=4000.0, Ly=4000.0, relax_width_m=267.0)

    # DEFECT: same cell count -> the physical band collapses by the refinement factor
    w_old = [nst.effective_relax_width(old, g) * g.dx for g in (coarse, fine)]
    assert w_old[0] > 2.5 * w_old[1]

    # FIX: the physical band is mesh-independent (within one cell)
    w_new = [nst.effective_relax_width(new, g) * g.dx for g in (coarse, fine)]
    assert abs(w_new[0] - w_new[1]) < max(coarse.dx, fine.dx)
    # ...and at the coarse mesh the fix is a no-op, so the control run is a clean comparison
    assert nst.effective_relax_width(new, coarse) == nst.effective_relax_width(old, coarse)


def test_relaxation_weight_default_is_byte_identical():
    g = _grid(60)
    spec = nst.NestSpec(x0=0.0, y0=0.0, Lx=4000.0, Ly=4000.0)
    got = np.asarray(g.backend.to_cpu(nst.relaxation_weight(g, spec)))[:, :, 0]
    i = np.arange(60)
    dist = np.minimum(np.minimum(i[:, None], 59 - i[:, None]),
                      np.minimum(i[None, :], 59 - i[None, :]))
    ref = 0.02 * np.clip(1.0 - dist / 4.0, 0.0, 1.0) ** 2
    assert np.array_equal(got, ref)


def test_relaxation_weight_physical_band_matches_across_meshes():
    """The nudged band covers the same PHYSICAL distance at 67 m and at 22 m."""
    reach = []
    for nx in (60, 180):
        g = _grid(nx)
        spec = nst.NestSpec(x0=0.0, y0=0.0, Lx=4000.0, Ly=4000.0, relax_width_m=267.0)
        w = np.asarray(g.backend.to_cpu(nst.relaxation_weight(g, spec)))[:, 0, 0]
        reach.append(int(np.count_nonzero(w[: nx // 2])) * g.dx)
    assert abs(reach[0] - reach[1]) < 80.0        # < ~4 fine cells apart


# --------------------------------------------------------------------------- tracker
def _two_vortex_sim(edge_strength=3.0, interior_strength=1.0):
    """A sim whose |zeta| MAXIMUM sits on the boundary (the sponge artifact) with a weaker,
    genuine vortex in the interior -- the situation that steered the 22 m nest into its wall."""
    scfg = build_storm_config(preset="storm", nx=40, ny=40, nz=16, Lx=40000.0, Ly=40000.0,
                              Lz=12000.0, duration=1.0, dt_max=3.0, device="cpu")
    sim = StormSimulation(scfg)
    g = sim.grid
    xc = np.asarray(g.backend.to_cpu(g.xc)); yc = np.asarray(g.backend.to_cpu(g.yc))
    u = np.zeros_like(np.asarray(g.backend.to_cpu(sim.state.u)))
    v = np.zeros_like(np.asarray(g.backend.to_cpu(sim.state.v)))
    s = 2.0 * g.dx

    def add(cx, cy, amp, U, V):
        X = xc[:, None] - cx; Y = yc[None, :] - cy
        e = amp * np.exp(-(X ** 2 + Y ** 2) / (2 * s ** 2))
        U[: len(xc), : len(yc), :] += (-e * Y / s)[:, :, None]
        V[: len(xc), : len(yc), :] += (e * X / s)[:, :, None]

    add(xc[2], yc[2], edge_strength, u, v)                 # boundary artifact
    add(xc[20], yc[20], interior_strength, u, v)           # the real vortex
    sim.state.u[...] = g.backend.asarray(u[: sim.state.u.shape[0], : sim.state.u.shape[1]])
    sim.state.v[...] = g.backend.asarray(v[: sim.state.v.shape[0], : sim.state.v.shape[1]])
    return sim


def test_tag_cells_border_excludes_the_sponge_artifact():
    sim = _two_vortex_sim()
    g = sim.grid
    box0 = nst.cluster_to_box(nst.tag_cells(sim.state, g, field="zeta", frac=0.5), margin=0)
    box6 = nst.cluster_to_box(nst.tag_cells(sim.state, g, field="zeta", frac=0.5, border=6),
                              margin=0)
    assert box0 is not None and box6 is not None
    c0 = box0[0] + box0[2] // 2
    c6 = box6[0] + box6[2] // 2
    assert c0 < 8, "default tracker should latch onto the boundary artifact (the defect)"
    assert c6 > 12, "border= must move the tracked centre to the interior vortex (the fix)"


def test_follow_spec_border_steers_away_from_the_wall():
    sim = _two_vortex_sim()
    # start OFF-centre from both vortices, so each tracker has somewhere to move
    old = nst.NestSpec.aligned(sim.grid, i0=10, j0=10, ncx=12, ncy=12, refine=3)
    near = nst.follow_spec(old, sim, field="zeta", frac=0.5, alpha=1.0, border=0)
    away = nst.follow_spec(old, sim, field="zeta", frac=0.5, alpha=1.0, border=6)
    assert away is not None
    # tracking the artifact drags the box toward the SW corner; excluding it does not
    if near is not None:
        assert away.x0 > near.x0 and away.y0 > near.y0
    assert away.Lx == old.Lx and away.Ly == old.Ly        # size preserved (exact integer shift)


def test_tag_cells_border_is_ignored_when_it_would_blank_the_grid():
    sim = _two_vortex_sim()
    tags = nst.tag_cells(sim.state, sim.grid, field="zeta", frac=0.5, border=50)
    assert tags.shape == (sim.grid.nx, sim.grid.ny)       # no crash, falls back to no exclusion
    assert tags.any()


def test_follow_spec_carries_the_sponge_across_a_move():
    """A moving nest must KEEP its sponge configuration.

    Rebuilding the spec with only (refine, nz, z_stretch) silently reverted
    relax_width/relax_rate/relax_width_m to their defaults on the FIRST re-centring, so a
    moving nest's 267 m physical band dropped back to 4 cells (89 m at dx=22 m) after one
    move -- i.e. the sponge fix was inactive for all but the first seconds of a run.
    """
    sim = _two_vortex_sim()
    old = nst.NestSpec.aligned(sim.grid, i0=10, j0=10, ncx=12, ncy=12, refine=3,
                               relax_width_m=267.0, relax_rate=0.05, relax_width=6)
    moved = nst.follow_spec(old, sim, field="zeta", frac=0.5, alpha=1.0, border=6)
    assert moved is not None and (moved.x0, moved.y0) != (old.x0, old.y0)   # it really moved
    assert moved.relax_width_m == old.relax_width_m == 267.0
    assert moved.relax_rate == old.relax_rate == 0.05
    assert moved.relax_width == old.relax_width == 6
    # and the band the moved nest actually builds is still the PHYSICAL one
    g_old = nst.build_nest_grid(old, sim.grid)
    g_new = nst.build_nest_grid(moved, sim.grid)
    assert nst.effective_relax_width(moved, g_new) == nst.effective_relax_width(old, g_old)
