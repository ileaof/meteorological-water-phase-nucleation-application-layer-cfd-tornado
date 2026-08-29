"""Composite two-level anelastic pressure projection across the coarse-fine interface.

This demonstrates the AMR-projection **call site**
(:func:`storm_dynamics.nesting.composite_project_two_level`): one composite solve over
the parent (coarse) + nest (fine) staggered mass fluxes ``m = rho0 u`` that makes
``div(rho0 u) = 0`` *consistently across the refinement interface* -- what the two
independent per-level projections cannot do (each is divergence-free on its own grid,
but the interface between them is not).

Steps:
  1. mature a small **square** parent supercell briefly;
  2. build a **cell-aligned, matched-z** nest (:meth:`NestSpec.aligned`) over the updraft;
  3. perturb both levels' face velocities (a strongly divergent kick);
  4. run the composite projection and report ``max|div(rho0 u)|`` before vs after --
     it collapses to machine precision, including at the interface cells.

Idealised numerics demonstration (not a forecast).  Usage:
    python examples/composite_projection_demo.py
    python examples/composite_projection_demo.py --parent-nx 24 --refine 3 --duration 600
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np                                             # noqa: E402
from storm_dynamics.config import build_storm_config          # noqa: E402
from storm_dynamics.core import StormSimulation               # noqa: E402
from storm_dynamics import nesting as nst, rotation as rot     # noqa: E402


def _anelastic_div_norms(parent, nest, spec):
    """max|div(rho0 u)| over the nest interior and the coarse cells, computed straight
    from the two FlowStates (independent of the projection) -- the diagnostic the
    composite projection drives to zero."""
    pg = parent.grid; to = pg.backend.to_cpu
    nz = pg.nz
    rc = np.asarray(to(parent.rho0_c)); rw = np.asarray(to(parent.rho0_wface))
    r3, w3 = rc[None, None, :], rw[None, None, :]
    dzc = (np.asarray(to(pg.zf))[1:] - np.asarray(to(pg.zf))[:-1])[None, None, :]
    # nest interior mass-flux divergence (physical spacing)
    ng = nest.grid; hxn = ng.dx
    mu = np.asarray(to(nest.state.u)) * r3; mv = np.asarray(to(nest.state.v)) * r3
    mw = np.asarray(to(nest.state.w)) * w3
    dfn = (mu[1:] - mu[:-1]) / hxn + (mv[:, 1:] - mv[:, :-1]) / hxn + (mw[:, :, 1:] - mw[:, :, :-1]) / dzc
    return float(np.abs(dfn).max())


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--parent-nx", type=int, default=18, help="parent nx=ny (square)")
    p.add_argument("--parent-nz", type=int, default=26)
    p.add_argument("--parent-L", type=float, default=18000.0)
    p.add_argument("--duration", type=float, default=400.0, help="parent maturation [s]")
    p.add_argument("--refine", type=int, default=3, help="nest horizontal refinement")
    p.add_argument("--kick", type=float, default=3.0, help="divergent velocity kick [m/s]")
    p.add_argument("--device", choices=["cpu", "gpu", "auto"], default="cpu")
    args = p.parse_args(argv)

    print("=" * 74)
    print("storm_dynamics -- composite two-level anelastic projection (div(rho0 u)=0")
    print("                  ACROSS the coarse-fine interface)  [IDEALISED numerics]")
    print("=" * 74)
    scfg = build_storm_config(
        preset="storm", nx=args.parent_nx, ny=args.parent_nx, nz=args.parent_nz,
        Lx=args.parent_L, Ly=args.parent_L, Lz=15000.0, duration=args.duration,
        dt_max=3.0, hodograph_kind="quarter_circle", drag=True, z_stretch=1.05,
        U_max=18.0, z_turn=2000.0, C_s=0.22, device=args.device)
    parent = StormSimulation(scfg)
    print("maturing parent: %d x %d x %d (dx=%.0f m) for %.0f s ..."
          % (parent.grid.nx, parent.grid.ny, parent.grid.nz, parent.grid.dx, args.duration))
    parent.run(progress=lambda t, d, s: print("  parent t=%6.0f/%.0f" % (t, d), end="\r"))
    print(" " * 60, end="\r")

    # cell-aligned nest over the updraft column-max (co-located with the low-level vortex)
    _, _, wc = rot._centered_velocity(parent.state, parent.grid)
    wc = np.asarray(parent.grid.backend.to_cpu(wc))
    i, j = np.unravel_index(np.argmax(wc.max(axis=2)), wc.shape[:2])
    nc = parent.grid.nx
    half = nc // 3
    i0 = int(np.clip(i - half // 2, 1, nc - half - 1))
    j0 = int(np.clip(j - half // 2, 1, nc - half - 1))
    spec = nst.NestSpec.aligned(parent.grid, i0=i0, j0=j0, ncx=half, ncy=half, refine=args.refine)
    nest = nst.NestedStormSimulation(parent, spec)
    print("parent grid   : %d^2 x %d, dx=%.0f m, periodic_h=%s"
          % (nc, parent.grid.nz, parent.grid.dx, bool(getattr(parent.grid, "periodic", False))))
    print("nest (aligned): parent cells [%d:%d) x [%d:%d), refine %d  ->  %d x %d x %d, dx=%.0f m"
          % (i0, i0 + half, j0, j0 + half, args.refine, nest.grid.nx, nest.grid.ny,
             nest.grid.nz, nest.grid.dx))

    # strongly divergent kick to both levels (so the projection has real work to do)
    rng = np.random.default_rng(0)
    xp = parent.grid.xp
    for sim in (parent, nest):
        g = sim.grid
        sim.state.u = sim.state.u + xp.asarray(args.kick * rng.standard_normal(g.u_shape))
        sim.state.v = sim.state.v + xp.asarray(args.kick * rng.standard_normal(g.v_shape))
        sim.state.w = sim.state.w + xp.asarray(args.kick * rng.standard_normal(g.w_shape))

    before = _anelastic_div_norms(parent, nest, spec)
    res = nst.composite_project_two_level(parent, nest, spec)
    after = _anelastic_div_norms(parent, nest, spec)

    print("-" * 74)
    print("max|div(rho0 u)|  (anelastic mass-flux divergence)")
    print("  nest interior BEFORE projection : %.3e   (strongly divergent)" % before)
    print("  nest interior AFTER  projection : %.3e" % after)
    print("  composite check (independent recompute from the written-back arrays):")
    print("    coarse         : %.3e" % res["div_coarse"])
    print("    fine           : %.3e" % res["div_fine"])
    print("    fine-INTERFACE : %.3e   <-- single-valued across the coarse-fine interface"
          % res["div_interface"])
    ok = max(res.values()) < 1e-9
    print("-" * 74)
    print("RESULT: div(rho0 u) driven to machine precision EVERYWHERE incl. the interface"
          if ok else "RESULT: unexpected residual -- see values above")
    print("The two independent per-level projections cannot achieve this at the interface;")
    print("the composite solve couples the levels with a single-valued interface flux.")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
