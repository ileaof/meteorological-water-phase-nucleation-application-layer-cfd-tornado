"""storm_dynamics milestone regression tests (heavier integrations).

M1 -- rotating supercell: unidirectional shear produces storm splitting (cyclonic
AND anticyclonic vertical vorticity) and a mid-level mesocyclone (Klemp-Wilhelmson
/ Weisman-Klemp).  M2 -- low-level rotation: a curved hodograph + surface drag +
evaporative cold pool spins up near-surface vertical vorticity (tornadogenesis
proxy).  These run a coarse demonstration grid for a few (wall-clock) tens of
seconds each; they assert the qualitative rotation signatures + conservation, not
quantitative magnitudes (the vortex is under-resolved -- see
docs/storm_dynamics_guide.md).
"""
from __future__ import annotations

import numpy as np

from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation


def _run(hodograph_kind, drag, nx=24, ny=24, nz=36, duration=1200.0,
         Lx=36000, Ly=36000, Lz=15000, **kw):
    scfg = build_storm_config(preset="storm", nx=nx, ny=ny, nz=nz,
                              Lx=Lx, Ly=Ly, Lz=Lz, duration=duration,
                              dt_max=3.0, hodograph_kind=hodograph_kind, drag=drag, **kw)
    sim = StormSimulation(scfg)
    report = sim.run()
    return sim, report


def test_M1_rotating_supercell_splits_and_forms_mesocyclone():
    sim, report = _run("unidirectional", drag=False)
    r = report["rotation"]
    pk = report["rotation_peak"]
    c = report["conservation"]
    # a deep updraft develops
    assert r["w_max"] > 8.0, r
    # storm splitting: BOTH cyclonic and anticyclonic vertical vorticity
    assert r["zeta_max"] > 3e-3, r
    assert r["zeta_min"] < -3e-3, r
    # mid-level mesocyclone above threshold (~O(1e-2) at demo resolution)
    assert pk["peak_midlevel_mesocyclone"] > 4e-3, pk
    # rotating updraft (updraft helicity clearly supercellular)
    assert pk["peak_updraft_helicity"] > 80.0, pk
    # conservation (repo standard) -- projection enforces continuity, not limiters
    assert abs(c["total_water_rel_err"]) < 1e-2, c
    assert c["mass_continuity_residual_norm"] < 1e-2, c


def test_M1_environment_is_supercellular():
    """The M1 environment carries supercell-strength deep-layer shear."""
    from storm_dynamics import soundings as snd
    sim, _ = _run("unidirectional", drag=False, duration=1.0)   # just build the env
    assert snd.bulk_shear(sim.base, 0.0, 6000.0) > 18.0
    assert snd.storm_relative_helicity(sim.base, 3000.0) > 0.0  # right-mover


def test_M2_low_level_rotation_under_curved_hodograph():
    """M2: curved hodograph + surface drag + evaporative cold pool spins up and
    SUSTAINS near-surface vertical vorticity (the tornadogenesis proxy), with a
    physical (non-blown-up) updraft.  Near-surface levels are clustered
    (z_stretch) so the drag / corner-flow layer is resolved."""
    sim, report = _run("quarter_circle", drag=True, nx=24, ny=24, nz=38,
                       Lx=32000, Ly=32000, Lz=15000, duration=1800.0,
                       z_stretch=1.05, U_max=18.0, z_turn=2000.0, C_s=0.22)
    r = report["rotation"]
    pk = report["rotation_peak"]
    c = report["conservation"]
    w_hist = max(h["w_max"] for h in sim.history)
    # the updraft stays PHYSICAL -- no grid-scale blow-up (parcel limit ~ sqrt(2 CAPE))
    assert 6.0 < w_hist < 45.0, w_hist
    # near-surface vertical vorticity develops (present) ...
    assert pk["peak_near_surface_zeta"] > 1.8e-3, pk
    # ... and is SUSTAINED, not a transient spike that decays away
    assert r["near_surface_zeta_max"] > 1.5e-3, r
    # conservation holds (projection, not limiters)
    assert abs(c["total_water_rel_err"]) < 1.5e-2, c
    assert c["mass_continuity_residual_norm"] < 1e-2, c


def test_M2_curved_hodograph_has_low_level_srh():
    """The M2 environment carries the low-level storm-relative helicity that the
    straight-hodograph M1 environment lacks."""
    from storm_dynamics import soundings as snd
    sim, _ = _run("quarter_circle", drag=True, duration=1.0,
                  U_max=18.0, z_turn=2000.0, z_stretch=1.05)
    straight, _ = _run("unidirectional", drag=False, duration=1.0)
    assert snd.storm_relative_helicity(sim.base, 1000.0) > \
        snd.storm_relative_helicity(straight.base, 1000.0) + 50.0
