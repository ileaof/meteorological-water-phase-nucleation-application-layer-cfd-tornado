"""Surface / boundary-layer sensitivity for the SURFACE-CONNECTION gap.

Attempt J showed the resolved vortex is ELEVATED (V_rot grows with height; it does not reach the
ground).  The remaining physics is the near-surface corner-flow: surface drag drives the convergent
inflow that concentrates and stretches vorticity at the lowest levels, and the first cell-centre
height sets what the mesh can even represent.  This runs the freely-evolving supercell across a
surface matrix and measures, with `vortex_diagnostics.surface_connection_report`, whether the vortex
descends (surface/aloft V_rot ratio) -- the spec's roughness/drag sensitivity analysis.

    python examples/surface_sensitivity.py            # fast (small grid, short runs)
    python examples/surface_sensitivity.py --full     # larger grid, longer runs, GPU
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np

from meteorological_flow.base_state import BaseState
from storm_dynamics.config import build_storm_config, SurfaceFluxConfig
from storm_dynamics.core import StormSimulation
from storm_dynamics import soundings as snd, vortex_diagnostics as vd


# (name, C_d, z_stretch, surface fluxes on?) -- drag & near-surface resolution sensitivity
CASES = [
    ("baseline",        dict(C_d=0.012, z_stretch=1.05)),
    ("no_drag",         dict(C_d=0.0,   z_stretch=1.05)),   # control: drag removed
    ("rough_surface",   dict(C_d=0.024, z_stretch=1.05)),   # z0 ~ rougher
    ("smooth_surface",  dict(C_d=0.004, z_stretch=1.05)),
    ("coarse+loglaw",   dict(C_d=0.012, z_stretch=1.05, log_law=True)),
    ("fine_near_sfc",   dict(C_d=0.012, z_stretch=1.12)),   # much finer first cells, bulk C_d
    ("fine+loglaw",     dict(C_d=0.012, z_stretch=1.12, log_law=True)),   # the physical combination
    ("fine+no_drag",    dict(C_d=0.0,   z_stretch=1.12)),
    ("drag+fluxes",     dict(C_d=0.012, z_stretch=1.05, fluxes=True)),
]


def _run(name, ov, nx, nz, steps, device):
    ov = dict(ov); fluxes = ov.pop("fluxes", False); log_law = ov.pop("log_law", False)
    scfg = build_storm_config(preset="storm", nx=nx, ny=nx, nz=nz, Lx=nx * 600.0, Ly=nx * 600.0,
                              Lz=15000.0, duration=1.0, dt_max=3.0, drag=(ov["C_d"] > 0),
                              z_stretch=ov["z_stretch"], C_s=0.20,
                              hodograph_kind="quarter_circle", U_max=30.0, device=device)
    scfg.sim.physics.bubble_dtheta = 5.0
    if ov["C_d"] > 0:
        scfg.dyn.drag.C_d = ov["C_d"]
    if log_law:
        scfg.dyn.drag.use_log_law = True
        scfg.dyn.drag.roughness_length_m = 0.1
    if fluxes:
        scfg.dyn.fluxes = SurfaceFluxConfig(enabled=True, C_h=1.2e-3, C_q=1.2e-3,
                                            dtheta_sfc_K=2.0, saturate_surface=True)
    s0 = StormSimulation(scfg); b = s0.base
    cx, cy = snd.bunkers_storm_motion(b)
    base = BaseState(zc=b.zc, theta0=b.theta0, qv0=b.qv0, p0=b.p0, T0=b.T0, rho0=b.rho0,
                     u0=b.u0 - cx, v0=b.v0 - cy)                    # storm-relative
    sim = StormSimulation(scfg, base=base)
    for _ in range(steps):
        sim._step(sim._dt()); sim.step += 1; sim.t = float(sim.state.t)
        if not np.isfinite(np.asarray(sim.grid.backend.to_cpu(sim.state.w))).all():
            break
    sim.state.diagnose(sim.cfg)
    from storm_dynamics.surface_drag import effective_drag_coefficient
    cd_eff = effective_drag_coefficient(sim.grid, scfg.dyn.drag) if scfg.dyn.drag.enabled else 0.0
    r = vd.surface_connection_report(sim.state, sim.grid)
    p = r["profile"]
    return {"name": name, "C_d": cd_eff, "dz1_m": r["first_cell_height_m"],
            "v_sfc": p[0]["v_rot_m_s"], "v_aloft": max(q["v_rot_m_s"] for q in p),
            "ratio": r["surface_aloft_ratio"], "connected": r["surface_connected"],
            "conv": r["near_surface_convergence_s"]}


def run_matrix(nx=40, nz=44, steps=250, device="cpu"):
    rows = [_run(n, o, nx, nz, steps, device) for n, o in CASES]
    hdr = (f"{'case':<16}{'C_d_eff':>9}{'dz1[m]':>8}{'V_sfc':>8}{'V_aloft':>9}"
           f"{'sfc/aloft':>11}{'conv':>10}  connected")
    print(hdr); print("-" * len(hdr))
    for r in rows:
        print(f"{r['name']:<16}{r['C_d']:>9.4f}{r['dz1_m']:>8.1f}{r['v_sfc']:>8.2f}"
              f"{r['v_aloft']:>9.2f}{r['ratio']:>11.2f}{r['conv']:>10.2e}  {r['connected']}")
    print("\nsurface/aloft ratio > 0.8 => surface-connected; < 0.8 => ELEVATED vortex.")
    print("NOTE: ground contact is never claimed below the first cell-centre height dz1.")
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--full", action="store_true")
    p.add_argument("--device", default="cpu")
    a = p.parse_args(argv)
    if a.full:
        run_matrix(nx=96, nz=48, steps=1200, device=a.device)
    else:
        run_matrix(nx=40, nz=44, steps=250, device=a.device)
    return 0


if __name__ == "__main__":
    sys.exit(main())
