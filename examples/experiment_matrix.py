"""Parametric + control experiment matrix for tornadogenesis mechanisms.

Runs short idealised simulations that vary one ingredient at a time and classifies each with the
objective diagnostics, to show *which* mechanisms are necessary for a supercell / mesocyclone /
low-level rotation.  Controls include: no shear, no surface drag, weaker CAPE, coarse resolution.
Nothing is imposed -- rotation, where it appears, is grown from the equations.

    python examples/experiment_matrix.py            # fast demo (small grid, short runs)
    python examples/experiment_matrix.py --full     # larger grid, longer runs

These are deliberately short illustrative runs (the framework + control logic), not converged
supercells; see docs/tornadogenesis_modeling.md for full runs.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np

from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics import rotation as rot, classification as cl


# (name, overrides) -- overrides go to build_storm_config; U_max=0 removes the shear (control)
EXPERIMENTS = [
    ("baseline",        dict()),
    ("no_shear",        dict(U_max=0.0)),
    ("no_surface_drag", dict(drag=False)),
    ("weak_bubble",     dict(_bubble=1.0)),
    ("coarse",          dict(_coarse=True)),
    ("tke15_closure",   dict(les_model="tke15")),
]


def _run_one(name, ov, nx, nz, steps, device):
    ov = dict(ov)
    bubble = ov.pop("_bubble", None)
    coarse = ov.pop("_coarse", False)
    n = max(16, nx // 2) if coarse else nx
    L = 40000.0
    scfg = build_storm_config(preset="storm", nx=n, ny=n, nz=nz, Lx=L, Ly=L, Lz=15000.0,
                              duration=1.0, dt_max=4.0, z_stretch=1.05, device=device, **ov)
    if bubble is not None:
        scfg.sim.physics.bubble_dtheta = bubble
    sim = StormSimulation(scfg)
    for _ in range(steps):
        sim._step(sim._dt()); sim.step += 1; sim.t = float(sim.state.t)
    sim.state.diagnose(sim.cfg)
    rep = rot.rotation_report(sim.state, sim.grid, base=sim.base)
    cat = cl.classify_simulation(sim, z_surface_m=200.0)["category"]
    return {"name": name, "dx_m": round(sim.grid.dx), "w_max": rep["w_max"],
            "midlevel_meso": rep["midlevel_mesocyclone"],
            "near_surface_zeta": rep["near_surface_zeta_max"],
            "UH_2_5km": rep["updraft_helicity_max"], "category": cat}


def run_matrix(nx=32, nz=40, steps=80, device="cpu"):
    rows = [_run_one(nm, ov, nx, nz, steps, device) for nm, ov in EXPERIMENTS]
    hdr = f"{'experiment':<16}{'dx':>6}{'w_max':>8}{'meso':>10}{'sfc_zeta':>10}{'UH':>8}  category"
    print(hdr); print("-" * len(hdr))
    for r in rows:
        print(f"{r['name']:<16}{r['dx_m']:>6}{r['w_max']:>8.1f}{r['midlevel_meso']:>10.2e}"
              f"{r['near_surface_zeta']:>10.2e}{r['UH_2_5km']:>8.1f}  {r['category']}")
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--full", action="store_true")
    p.add_argument("--device", default="cpu")
    a = p.parse_args(argv)
    if a.full:
        run_matrix(nx=64, nz=44, steps=300, device=a.device)
    else:
        run_matrix(nx=32, nz=40, steps=80, device=a.device)
    return 0


if __name__ == "__main__":
    sys.exit(main())
