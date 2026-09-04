"""Height-consistent (neutral log-law / MOST) surface drag closure.

A fixed bulk C_d is calibrated for one first-cell height; the surface-sensitivity study showed it
over-damps the lowest level once the corner-flow layer is resolved.  The log-law closure ties C_d to
the roughness and the ACTUAL first cell-centre height."""
import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from meteorological_flow.grid import Grid
from meteorological_flow.state import FlowState
from storm_dynamics.config import SurfaceDragConfig
from storm_dynamics import surface_drag as sd


def test_log_law_formula():
    cd = sd.log_law_drag_coefficient(10.0, 0.1)
    assert np.isclose(cd, (0.4 / math.log(100.0)) ** 2, rtol=1e-12)


def test_cd_increases_as_mesh_refines():
    """Closer to the ground the log-law C_d is LARGER (the log shrinks) -- the resolution
    dependence a fixed bulk constant cannot express."""
    cd_50 = sd.log_law_drag_coefficient(50.0, 0.1)
    cd_6 = sd.log_law_drag_coefficient(6.0, 0.1)
    assert cd_6 > cd_50 > 0.0


def test_rougher_surface_gives_larger_cd():
    assert sd.log_law_drag_coefficient(10.0, 0.5) > sd.log_law_drag_coefficient(10.0, 0.01)


def test_effective_coefficient_switches():
    g = Grid(nx=8, ny=8, nz=20, Lx=800.0, Ly=800.0, Lz=2000.0, z_stretch=1.05, periodic=True)
    bulk = SurfaceDragConfig(enabled=True, C_d=0.012)
    assert sd.effective_drag_coefficient(g, bulk) == 0.012            # unchanged default
    loglaw = SurfaceDragConfig(enabled=True, C_d=0.012, use_log_law=True, roughness_length_m=0.1)
    z1 = float(np.asarray(g.backend.to_cpu(g.zc))[0])
    assert np.isclose(sd.effective_drag_coefficient(g, loglaw),
                      sd.log_law_drag_coefficient(z1, 0.1), rtol=1e-12)


def test_default_bulk_path_is_unchanged():
    """use_log_law=False must retard the lowest level exactly as before."""
    g = Grid(nx=8, ny=8, nz=20, Lx=800.0, Ly=800.0, Lz=2000.0, z_stretch=1.05, periodic=True)
    out = []
    for cfg in (SurfaceDragConfig(enabled=True, C_d=0.012),
                SurfaceDragConfig(enabled=True, C_d=0.012)):
        st = FlowState.zeros(g); st.u[:] = 10.0; st.v[:] = 4.0
        sd.apply_surface_drag(st, g, 2.0, cfg)
        out.append(np.asarray(st.u[:, :, 0]).copy())
    assert np.array_equal(out[0], out[1])
    # and the drag actually retards
    st = FlowState.zeros(g); st.u[:] = 10.0
    sd.apply_surface_drag(st, g, 2.0, SurfaceDragConfig(enabled=True, C_d=0.012))
    assert float(np.asarray(st.u[:, :, 0]).max()) < 10.0


def test_log_law_retards_more_on_a_refined_mesh():
    """Same roughness, finer first cell -> stronger relative retardation of the lowest level."""
    frac = {}
    for zs in (1.05, 1.20):
        g = Grid(nx=8, ny=8, nz=24, Lx=800.0, Ly=800.0, Lz=2000.0, z_stretch=zs, periodic=True)
        st = FlowState.zeros(g); st.u[:] = 10.0
        cfg = SurfaceDragConfig(enabled=True, use_log_law=True, roughness_length_m=0.1)
        sd.apply_surface_drag(st, g, 2.0, cfg)
        frac[zs] = float(np.asarray(st.u[:, :, 0]).max()) / 10.0
    assert frac[1.20] < frac[1.05]           # the refined mesh feels the (larger) log-law C_d


def test_stress_divergence_is_mesh_independent():
    """The defect the surface study found: the bulk sink rate C_d|V|/dz1 blows up as the mesh is
    refined (stripping the tangential wind).  The stress-divergence form spreads the stress over a
    PHYSICAL depth h, so the lowest cell's retention is independent of dz1."""
    ret_bulk, ret_div = {}, {}
    for zs in (1.05, 1.20):
        g = Grid(nx=8, ny=8, nz=26, Lx=800.0, Ly=800.0, Lz=2000.0, z_stretch=zs, periodic=True)
        for name, cfg, store in (
            ("bulk", SurfaceDragConfig(enabled=True, C_d=0.012), ret_bulk),
            ("div", SurfaceDragConfig(enabled=True, C_d=0.012, stress_divergence=True,
                                      surface_layer_depth_m=150.0), ret_div)):
            st = FlowState.zeros(g); st.u[:] = 10.0
            sd.apply_surface_drag(st, g, 2.0, cfg)
            store[zs] = float(np.asarray(st.u[:, :, 0]).max()) / 10.0
    # bulk: refining the mesh damps the lowest cell much harder
    assert ret_bulk[1.20] < ret_bulk[1.05] - 0.02
    # stress-divergence: essentially unchanged by the mesh
    assert abs(ret_div[1.20] - ret_div[1.05]) < 1e-3


def test_stress_divergence_spreads_through_the_layer():
    """The bulk law touches only k=0; the divergence form retards every level inside h."""
    g = Grid(nx=8, ny=8, nz=30, Lx=800.0, Ly=800.0, Lz=3000.0, z_stretch=1.05, periodic=True)
    zc = np.asarray(g.backend.to_cpu(g.zc))
    h = 150.0
    st_b = FlowState.zeros(g); st_b.u[:] = 10.0
    sd.apply_surface_drag(st_b, g, 5.0, SurfaceDragConfig(enabled=True, C_d=0.012))
    st_d = FlowState.zeros(g); st_d.u[:] = 10.0
    sd.apply_surface_drag(st_d, g, 5.0, SurfaceDragConfig(enabled=True, C_d=0.012,
                                                          stress_divergence=True,
                                                          surface_layer_depth_m=h))
    k_in = int(np.where(zc < h)[0][-1])            # a level inside the layer, above k=0
    assert k_in >= 1
    assert np.isclose(float(np.asarray(st_b.u[:, :, k_in]).max()), 10.0)      # bulk: untouched
    assert float(np.asarray(st_d.u[:, :, k_in]).max()) < 10.0                 # divergence: retarded
    k_out = int(np.where(zc >= h)[0][0])
    assert np.isclose(float(np.asarray(st_d.u[:, :, k_out]).max()), 10.0)     # nothing above h


if __name__ == "__main__":
    for nm, fn in sorted(globals().items()):
        if nm.startswith("test_") and callable(fn):
            fn(); print("ok", nm)
    print("ALL LOG-LAW DRAG TESTS PASSED")


def test_log_law_reference_height_removes_the_dz1_confound():
    """C_d = (kappa/ln(z1/z0))^2 genuinely depends on the reference height -- correct physics, and
    the whole point of the log-law closure.  But it means a dz1 sweep varies the DRAG at the same
    time as the RESOLUTION: measured 0.0045 at z1=39.9 m vs 0.0146 at 2.73 m, a 3.2x spread.  A
    fine-dz1 run then came back with a much lower peak surface wind than its coarse counterpart,
    and the tripled drag is at least as plausible an explanation as the mesh.

    Pinning the evaluation height makes every member of the sweep apply an IDENTICAL C_d.
    """
    import numpy as np
    from storm_dynamics import surface_drag as sd
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation

    default_cd, pinned_cd = [], []
    for nz, zs in ((48, 1.05), (64, 1.06), (64, 1.09)):
        cfg = build_storm_config(preset="storm", nx=16, ny=16, nz=nz, Lx=72000.0, Ly=72000.0,
                                 Lz=15000.0, duration=1.0, dt_max=3.0, drag=True,
                                 z_stretch=zs, device="cpu")
        cfg.dyn.drag.use_log_law = True
        cfg.dyn.drag.roughness_length_m = 0.1
        g = StormSimulation(cfg).grid
        default_cd.append(sd.effective_drag_coefficient(g, cfg.dyn.drag))
        cfg.dyn.drag.log_law_reference_height_m = 10.0
        pinned_cd.append(sd.effective_drag_coefficient(g, cfg.dyn.drag))

    # the CONFOUND: mesh-following C_d spreads widely across a dz1 sweep
    assert max(default_cd) / min(default_cd) > 3.0, default_cd
    # the FIX: pinned C_d is identical on every mesh
    assert max(pinned_cd) - min(pinned_cd) < 1e-12, pinned_cd
    # and it equals the log law evaluated at the reference height
    assert pinned_cd[0] == pytest.approx(sd.log_law_drag_coefficient(10.0, 0.1, 0.4))


def test_log_law_reference_height_default_is_mesh_following():
    """Opt-in: None must reproduce the previous behaviour (evaluate at the first cell centre)."""
    import numpy as np
    from storm_dynamics import surface_drag as sd
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    cfg = build_storm_config(preset="storm", nx=16, ny=16, nz=64, Lx=72000.0, Ly=72000.0,
                             Lz=15000.0, duration=1.0, dt_max=3.0, drag=True, z_stretch=1.09,
                             device="cpu")
    cfg.dyn.drag.use_log_law = True
    cfg.dyn.drag.roughness_length_m = 0.1
    assert cfg.dyn.drag.log_law_reference_height_m is None      # default unchanged
    g = StormSimulation(cfg).grid
    z1 = float(np.asarray(g.backend.to_cpu(g.zc))[0])
    assert sd.effective_drag_coefficient(g, cfg.dyn.drag) == pytest.approx(
        sd.log_law_drag_coefficient(z1, 0.1, 0.4))
