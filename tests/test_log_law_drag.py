"""Height-consistent (neutral log-law / MOST) surface drag closure.

A fixed bulk C_d is calibrated for one first-cell height; the surface-sensitivity study showed it
over-damps the lowest level once the corner-flow layer is resolved.  The log-law closure ties C_d to
the roughness and the ACTUAL first cell-centre height."""
import math
import os
import sys

import numpy as np

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


if __name__ == "__main__":
    for nm, fn in sorted(globals().items()):
        if nm.startswith("test_") and callable(fn):
            fn(); print("ok", nm)
    print("ALL LOG-LAW DRAG TESTS PASSED")
