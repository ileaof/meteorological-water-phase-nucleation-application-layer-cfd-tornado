"""Physical scales must survive a change of mesh -- the study's most damaging failure mode.

Three independent defects in this project shared one shape: a quantity that is physically a LENGTH
was stored or converted as a CELL COUNT, so it silently rescaled with dx and corrupted every
cross-resolution comparison.  The worst was

    R = max(3, int(radius_m / grid.dx))          # surface_connection_report

which turned a requested 400 m sampling radius into 1800 m at dx=600 m and 900 m at dx=300 m --
a 2x mismatch in exactly the quantity a resolution study measures, biased toward the coarse mesh.

These tests pin the contract: a physical request is honoured or explicitly refused, never
silently substituted, and never enlarged.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from meteorological_flow.state import FlowState
from storm_dynamics import scales
from storm_dynamics import vortex_diagnostics as vd
from storm_dynamics.core import Grid

# the four meshes this study actually used
STUDY_DX = (600.0, 300.0, 66.7, 22.2)


# --------------------------------------------------------------------- the conversion contract
def test_never_enlarges_a_physical_region():
    """floor, not round: a represented region may be smaller than requested, never larger.
    Rounding up is how the old code inflated a coarse-grid measurement."""
    for dx in STUDY_DX:
        for req in (400.0, 1000.0, 3000.0, 137.0):
            sc = scales.cells_for_length(req, dx, min_cells=1)
            assert sc.represented_m <= req + 1e-9, sc.describe()
            assert sc.relative_error <= 1e-12


def test_under_resolved_is_reported_not_substituted():
    """The defect: 400 m at dx=600 silently became 3 cells = 1800 m.  It must now be refused."""
    sc = scales.cells_for_length(400.0, 600.0, min_cells=3, name="V_rot radius")
    assert sc.resolved is False and sc.status == "under_resolved"
    assert sc.cells == 0
    assert sc.represented_m == 0.0          # NOT 1800.0 -- no silent substitution
    assert bool(sc) is False                # truthiness follows resolution


def test_strict_mode_raises_with_an_informative_message():
    with pytest.raises(scales.UnderResolvedError) as e:
        scales.cells_for_length(400.0, 600.0, min_cells=3, name="V_rot radius", strict=True)
    msg = str(e.value)
    assert "UNDER-RESOLVED" in msg and "400" in msg and "1800" in msg   # says what WOULD work


def test_same_request_is_physically_comparable_across_meshes():
    """REQUIRED TEST 1: the same physical region on several grids.  Where the request is
    resolvable the represented length must agree across meshes to within one coarse cell."""
    req = 3000.0
    got = {}
    for dx in STUDY_DX:
        sc = scales.cells_for_length(req, dx, min_cells=3)
        assert sc.resolved, sc.describe()
        got[dx] = sc.represented_m
    spread = max(got.values()) - min(got.values())
    assert spread <= max(STUDY_DX) + 1e-9, got
    # and every mesh is within a cell of the request
    for dx, rep in got.items():
        assert req - rep < dx + 1e-9


def test_common_comparison_length_is_bound_by_the_coarsest_mesh():
    c_ok = scales.common_comparison_length(3000.0, [600.0, 300.0], min_cells=3)
    assert c_ok.resolved and c_ok.dx_m == 600.0
    c_bad = scales.common_comparison_length(400.0, [600.0, 300.0], min_cells=3)
    assert c_bad.resolved is False            # the 600/300 study CANNOT measure 400 m fairly


def test_smallest_resolvable_length():
    assert scales.smallest_resolvable_length(600.0, 3) == 1800.0
    assert scales.smallest_resolvable_length(22.2, 3) == pytest.approx(66.6)


def test_rejects_nonsense_input():
    with pytest.raises(ValueError):
        scales.cells_for_length(100.0, 0.0)
    with pytest.raises(ValueError):
        scales.cells_for_length(-1.0, 10.0)


# ------------------------------------------------------- the diagnostic that carried the defect
def _vortex_state(g, amp=3.0, rc=300.0, cx=None, cy=None):
    xp = g.xp
    cx = 0.5 * g.Lx if cx is None else cx
    cy = 0.5 * g.Ly if cy is None else cy
    X = xp.asarray(g.xc)[:, None, None] + g.zeros_c()
    Y = xp.asarray(g.yc)[None, :, None] + g.zeros_c()
    r = xp.sqrt((X - cx) ** 2 + (Y - cy) ** 2) + 1e-9
    phi = xp.arctan2(Y - cy, X - cx)
    vth = amp * (r / rc) * xp.exp(-(r * r) / (2 * rc * rc))
    st = FlowState.zeros(g)
    uc = -vth * xp.sin(phi); vc = vth * xp.cos(phi)
    st.u[1:-1, :, :] = 0.5 * (uc[:-1] + uc[1:])
    st.v[:, 1:-1, :] = 0.5 * (vc[:, :-1] + vc[:, 1:])
    return st


def test_surface_connection_report_refuses_an_under_resolved_radius():
    """REQUIRED TEST 1 (diagnostic side): a 400 m radius must be honoured or explicitly rejected,
    never silently widened to 3 cells."""
    g = Grid(nx=24, ny=24, nz=20, Lx=14400.0, Ly=14400.0, Lz=3000.0, z_stretch=1.05, periodic=True)
    assert g.dx == pytest.approx(600.0)
    rep = vd.surface_connection_report(_vortex_state(g, rc=2000.0), g, radius_m=400.0)
    assert rep["valid"] is False
    assert rep["resolution"]["status"] == "under_resolved"
    assert rep["resolution"]["represented_m"] == 0.0        # not 1800
    assert all(np.isnan(p["v_rot_m_s"]) for p in rep["profile"])
    assert np.isnan(rep["surface_aloft_ratio"])
    assert rep["surface_connected"] is False                # never "connected" on invalid data
    # zeta is a point maximum and does not depend on the radius -- still reported
    assert all(np.isfinite(p["zeta_max_s"]) for p in rep["profile"])


def test_surface_connection_report_strict_raises():
    g = Grid(nx=24, ny=24, nz=20, Lx=14400.0, Ly=14400.0, Lz=3000.0, z_stretch=1.05, periodic=True)
    with pytest.raises(scales.UnderResolvedError):
        vd.surface_connection_report(_vortex_state(g, rc=2000.0), g, radius_m=400.0, strict=True)


def test_surface_connection_report_radius_is_comparable_across_two_meshes():
    """The same physical domain and vortex at two resolutions must report the SAME represented
    sampling radius -- the property the old max(3, int(...)) floor destroyed."""
    reps = {}
    for nx in (40, 80):                                   # dx = 300 m and 150 m over 12 km
        g = Grid(nx=nx, ny=nx, nz=20, Lx=12000.0, Ly=12000.0, Lz=3000.0, z_stretch=1.05,
                 periodic=True)
        rep = vd.surface_connection_report(_vortex_state(g, rc=1500.0), g, radius_m=3000.0)
        assert rep["valid"] is True, rep["resolution"]
        reps[g.dx] = rep["resolution"]["represented_m"]
    assert abs(reps[300.0] - reps[150.0]) <= 300.0 + 1e-9, reps
    # the OLD formula would have given max(3,int(3000/300))*300 = 3000 vs
    # max(3,int(3000/150))*150 = 3000 here; the failure appears when radius < 3*dx, asserted above


def test_physical_interior_margin_is_mesh_independent():
    """border_m excludes the same PHYSICAL width on every mesh; border_frac does not."""
    got = {}
    for nx in (40, 80):
        g = Grid(nx=nx, ny=nx, nz=20, Lx=12000.0, Ly=12000.0, Lz=3000.0, z_stretch=1.05,
                 periodic=True)
        rep = vd.surface_connection_report(_vortex_state(g, rc=1500.0), g, radius_m=3000.0,
                                           border_m=1200.0)
        got[g.dx] = rep["interior_margin"]["represented_m"]
    assert abs(got[300.0] - got[150.0]) <= 300.0 + 1e-9, got
    assert all(abs(v - 1200.0) <= 300.0 for v in got.values()), got


def test_report_carries_discretisation_provenance():
    """A cross-resolution claim must be checkable from the result alone."""
    g = Grid(nx=40, ny=40, nz=20, Lx=12000.0, Ly=12000.0, Lz=3000.0, z_stretch=1.05, periodic=True)
    rep = vd.surface_connection_report(_vortex_state(g, rc=1500.0), g, radius_m=3000.0)
    res = rep["resolution"]
    for key in ("requested_m", "represented_m", "cells", "dx_m", "min_cells", "resolved", "status"):
        assert key in res, key
    assert res["requested_m"] == 3000.0
    assert res["cells"] == 10 and res["dx_m"] == pytest.approx(300.0)


# ------------------------------------------------- A8: nest pressure must not drive classification
def test_pressure_branch_cannot_promote_on_a_nest():
    """A8: nest `p_dyn` carries a boundary-driven component ~1/dt (measured -1418, -3250, -5161,
    +7412 Pa against a cyclostrophic scale ~-120 Pa).  A -200 Pa threshold is meaningless there,
    so the pressure branch must not be able to reach a tornado-like tier on a nest."""
    from storm_dynamics import classification as cl
    diag = {"w_max": 40.0, "midlevel_mesocyclone": 1.0, "updraft_helicity_2_5km": 1e4,
            "near_surface_zeta_max": 1.0, "v_theta_max_m_s": 5.0,      # BELOW the v_theta gate
            "circulation_m2_s": 1.0e5, "pressure_deficit_Pa": -3000.0,  # contaminated nest value
            "level_m": 50.0, "gust_front_convergence_s": 1.0}
    with_p = cl.classify(diag, allow_pressure=True)
    without_p = cl.classify(diag, allow_pressure=False)
    assert with_p["criteria"]["tornado_like_via_pressure"] is True
    assert without_p["criteria"]["tornado_like_via_pressure"] is False
    assert without_p["criteria"]["pressure_branch_available"] is False
    assert without_p["rank"] < with_p["rank"], (without_p["category"], with_p["category"])


def test_v_theta_branch_is_unaffected_by_pressure_gating():
    """Gating pressure must not weaken the honest branch."""
    from storm_dynamics import classification as cl
    diag = {"w_max": 40.0, "midlevel_mesocyclone": 1.0, "updraft_helicity_2_5km": 1e4,
            "near_surface_zeta_max": 1.0, "v_theta_max_m_s": 25.0,     # ABOVE the v_theta gate
            "circulation_m2_s": 0.0, "pressure_deficit_Pa": 0.0,
            "level_m": 50.0, "gust_front_convergence_s": 1.0}
    a = cl.classify(diag, allow_pressure=True)
    b = cl.classify(diag, allow_pressure=False)
    assert a["category"] == b["category"] and a["criteria"]["tornado_like"] is True


def test_missing_pressure_is_not_treated_as_zero_deficit():
    """A missing deficit must not silently satisfy `<= -200 Pa` via a 0.0 default."""
    from storm_dynamics import classification as cl
    diag = {"w_max": 40.0, "midlevel_mesocyclone": 1.0, "updraft_helicity_2_5km": 1e4,
            "near_surface_zeta_max": 1.0, "v_theta_max_m_s": 5.0,
            "circulation_m2_s": 1.0e5, "level_m": 50.0, "gust_front_convergence_s": 1.0}
    out = cl.classify(diag, allow_pressure=True)
    assert out["criteria"]["pressure_branch_available"] is False
    assert out["criteria"]["tornado_like_via_pressure"] is False


# ------------------------------- the tangential-profile / classification defect (external audit)
def _vth_r(g, st, z_m=100.0):
    """v_theta and r about the domain centre at a level -- matches tangential_radial's API."""
    from storm_dynamics.rotation import _centered_velocity
    import numpy as _np
    zc = _np.asarray(g.backend.to_cpu(g.zc))
    k = int(_np.argmin(_np.abs(zc - z_m)))
    uc, vc, _ = _centered_velocity(st, g)
    vth, _vr, r = vd.tangential_radial(uc[:, :, k], vc[:, :, k], g, (0.5 * g.Lx, 0.5 * g.Ly))
    return vth, r


def test_tangential_profile_bin_width_tracks_the_MESH_not_a_fixed_count():
    """The defect: nbins was fixed at 24, so bin width = r_max/24 at EVERY resolution.  With the
    classifier's radius_m=1500 that is 62.5 m on any mesh -- and `core_radius_m` came back as
    exactly 62.5 m at dx=600 AND dx=300, i.e. it reported the binning, not the vortex.

    Correct behaviour: the number of bins spanning a FIXED physical radius must GROW as the mesh
    refines, and every bin must span at least `min_cells_per_bin` cells.  On a mesh too coarse to
    support even one such bin the profile collapses to a single bin -- which is the honest answer,
    not a 24-bin fiction.
    """
    got = {}
    for nx in (20, 40, 80):                                  # dx = 600, 300, 150 m over 12 km
        g = Grid(nx=nx, ny=nx, nz=8, Lx=12000.0, Ly=12000.0, Lz=2000.0, periodic=True)
        vth, r = _vth_r(g, _vortex_state(g, rc=1500.0))
        rc, prof, vmax, core = vd.tangential_profile(vth, r, g, r_max_m=1500.0)
        w = float(rc[1] - rc[0]) if len(rc) > 1 else float("nan")
        got[g.dx] = (len(rc), w)
        if len(rc) > 1:
            assert w >= 3.0 * g.dx - 1e-6, (g.dx, w)          # >= min_cells_per_bin cells
    # bin COUNT over a fixed physical radius must increase as the mesh refines
    counts = [got[dx][0] for dx in (600.0, 300.0, 150.0)]
    assert counts[-1] > counts[0], got
    # and the old constant 62.5 m core must be gone: core radius now varies with the mesh
    cores = set()
    for nx in (20, 40, 80):
        g = Grid(nx=nx, ny=nx, nz=8, Lx=12000.0, Ly=12000.0, Lz=2000.0, periodic=True)
        vth, r = _vth_r(g, _vortex_state(g, rc=1500.0))
        cores.add(round(vd.tangential_profile(vth, r, g, r_max_m=1500.0)[3], 1))
    assert len(cores) > 1, cores


def test_tangential_profile_excludes_r_zero_and_reports_nan_for_empty_bins():
    """At r=0, phi=arctan2(0,0)=0 so v_theta collapses to the raw meridional wind of ONE cell.
    That single number crossed the tornado gate.  It must be excluded, and unsampled radii must
    be NaN rather than a 0.0 that looks like a measurement."""
    g = Grid(nx=40, ny=40, nz=8, Lx=12000.0, Ly=12000.0, Lz=2000.0, periodic=True)
    st = _vortex_state(g, rc=1500.0)
    vth, r = _vth_r(g, st)
    rc, prof, vmax, core = vd.tangential_profile(vth, r, g, r_max_m=1500.0)
    assert rc[0] >= g.dx, "the centre cell (r=0) must be excluded"
    assert not np.any(prof == 0.0), "empty bins must be NaN, never 0.0"
    assert np.isfinite(vmax) and np.isfinite(core)


def test_classifier_refuses_a_tornado_tier_when_the_core_is_unresolved():
    """The live false claim: both uniform runs were classified SURFACE_CONNECTED_TORNADO_LIKE_
    VORTEX on a 600 m mesh with a 62.5 m 'core' -- 0.1 cells."""
    from storm_dynamics import classification as cl
    base = {"w_max": 40.0, "midlevel_mesocyclone": 1.0, "updraft_helicity_2_5km": 1e4,
            "near_surface_zeta_max": 1.0, "v_theta_max_m_s": 18.19, "circulation_m2_s": 0.0,
            "level_m": 50.0, "gust_front_convergence_s": 1.0}
    unresolved = dict(base, core_radius_m=62.5, dx_m=600.0)      # 0.10 cells
    out = cl.classify(unresolved)
    assert out["criteria"]["core_resolved"] is False
    assert out["category"] == "LOW_LEVEL_MESOCYCLONE"

    # the gate must NOT suppress a legitimately resolved core
    resolved = dict(base, core_radius_m=2400.0, dx_m=600.0)      # 4.0 cells
    assert cl.classify(resolved)["category"] == "SURFACE_CONNECTED_TORNADO_LIKE_VORTEX"
    same_core_fine_mesh = dict(base, core_radius_m=62.5, dx_m=20.0)   # 3.1 cells
    assert cl.classify(same_core_fine_mesh)["category"] == "SURFACE_CONNECTED_TORNADO_LIKE_VORTEX"


def test_classifier_gate_is_inert_without_dx_information():
    """Back-compatible: a diag dict with no dx_m cannot be gated, so it must classify as before."""
    from storm_dynamics import classification as cl
    d = {"w_max": 40.0, "midlevel_mesocyclone": 1.0, "updraft_helicity_2_5km": 1e4,
         "near_surface_zeta_max": 1.0, "v_theta_max_m_s": 18.19, "level_m": 50.0,
         "gust_front_convergence_s": 1.0}
    assert cl.classify(d)["criteria"]["core_resolved"] is True
