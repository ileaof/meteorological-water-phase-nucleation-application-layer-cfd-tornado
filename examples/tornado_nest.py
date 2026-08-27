"""M3 (phase 1) -- static nested-grid refinement of the low-level vortex.

Mature the storm on the coarse **parent** domain, then integrate a finer **nest**
over the updraft / low-level-rotation region for a short window, with the nest's
border relaxed toward the (frozen) parent state.  The finer grid *intensifies*
the near-surface vertical vorticity -- the vortex sharpening under refinement,
which is what M3 is about.

Honest scope (see docs/storm_dynamics_guide.md): this is a **one-way, static,
short-window** nest (the parent is frozen at the nest border), a *demonstration
of the method* -- NOT full AMR, not two-way, and at O(100 m--1 km) only
approaching a resolved vortex, never a forecast.

Usage:
    python examples/tornado_nest.py --device gpu --plots --animate
    python examples/tornado_nest.py --parent-duration 1500 --refine 4 --window 120
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np                                              # noqa: E402
from storm_dynamics.config import build_storm_config           # noqa: E402
from storm_dynamics.core import StormSimulation                # noqa: E402
from storm_dynamics import nesting as nst, rotation as rot      # noqa: E402


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--parent-nx", type=int, default=24)
    p.add_argument("--parent-nz", type=int, default=38)
    p.add_argument("--parent-L", type=float, default=32000.0)
    p.add_argument("--parent-duration", type=float, default=1200.0)
    p.add_argument("--refine", type=int, default=3, help="nest horizontal refinement factor")
    p.add_argument("--half", type=float, default=8000.0, help="nest half-width [m]")
    p.add_argument("--nest-nz", type=int, default=44)
    p.add_argument("--window", type=float, default=120.0, help="nest integration window [s]")
    p.add_argument("--device", choices=["cpu", "gpu", "auto"], default="cpu")
    p.add_argument("--plots", action="store_true")
    p.add_argument("--animate", action="store_true")
    p.add_argument("--outdir", default="outputs/tornado_nest")
    args = p.parse_args(argv)

    print("=" * 72)
    print("storm_dynamics -- M3 phase 1: static nested-grid refinement (IDEALISED)")
    print("=" * 72)
    scfg = build_storm_config(
        preset="storm", nx=args.parent_nx, ny=args.parent_nx, nz=args.parent_nz,
        Lx=args.parent_L, Ly=args.parent_L, Lz=15000.0,
        duration=args.parent_duration, dt_max=3.0, hodograph_kind="quarter_circle",
        drag=True, z_stretch=1.05, U_max=18.0, z_turn=2000.0, C_s=0.22, device=args.device)
    parent = StormSimulation(scfg)
    print("maturing parent: %d x %d x %d (dx=%.0f m) for %.0f s ..."
          % (parent.grid.nx, parent.grid.ny, parent.grid.nz, parent.grid.dx,
             args.parent_duration))
    parent.run(progress=lambda t, d, s: print("  parent t=%6.0f/%.0f" % (t, d), end="\r"))
    print(" " * 60, end="\r")

    # place the nest on the updraft column-max (co-located with the low-level vortex)
    _, _, wc = rot._centered_velocity(parent.state, parent.grid)
    wc = np.asarray(parent.grid.backend.to_cpu(wc))
    i, j = np.unravel_index(np.argmax(wc.max(axis=2)), wc.shape[:2])
    xc = float(parent.grid.xc[i]); yc = float(parent.grid.yc[j])
    r_par = rot.rotation_report(parent.state, parent.grid)

    spec = nst.NestSpec.around(parent.grid, xc, yc, half=args.half,
                               refine=args.refine, nz=args.nest_nz, z_stretch=1.06)
    nest = nst.NestedStormSimulation(parent, spec)
    nest.cfg.time.duration = args.window
    r_nest0 = rot.rotation_report(nest.state, nest.grid)
    print("nest region   : centred (%.1f, %.1f) km, half %.0f km" % (xc / 1000, yc / 1000, args.half / 1000))
    print("nest grid     : %d x %d x %d  dx=%.0f m (parent %.0f m, %dx finer), dz0=%.0f m"
          % (nest.grid.nx, nest.grid.ny, nest.grid.nz, nest.grid.dx, parent.grid.dx,
             args.refine, float(nest.grid.dz_c[0])))
    print("-" * 72)
    rep = nest.run(progress=lambda t, d, s: print("  nest   t=%6.0f/%.0f" % (t, d), end="\r"),
                   capture_frames=args.animate)
    print(" " * 60, end="\r")

    r = rep["rotation"]; pk = rep["rotation_peak"]; c = rep["conservation"]
    print("near-surface zeta [s^-1]:")
    print("  parent (coarse, dx=%.0f m)        : %.2e" % (parent.grid.dx, r_par["near_surface_zeta_max"]))
    print("  nest initial (interpolated)        : %.2e" % r_nest0["near_surface_zeta_max"])
    print("  nest after %.0f s (dx=%.0f m)   : %.2e (peak %.2e)"
          % (args.window, nest.grid.dx, r["near_surface_zeta_max"], pk["peak_near_surface_zeta"]))
    intens = pk["peak_near_surface_zeta"] / max(r_nest0["near_surface_zeta_max"], 1e-12)
    print("  -> intensification over the window : %.2fx" % intens)
    print("w_max: parent %.1f -> nest init %.1f -> nest final %.1f m/s"
          % (r_par["w_max"], r_nest0["w_max"], r["w_max"]))
    print("conservation (nest): water %.2e | mass-continuity |div| %.2e"
          % (c["total_water_rel_err"], c["mass_continuity_residual_norm"]))
    print("=" * 72)
    for lim in rep["limitations"]:
        print("NOTE:", lim)
    print("NOTE: one-way STATIC nest, frozen parent border, short window -- a "
          "demonstration of the refinement method, not full AMR (M3 phase 1).")

    if args.plots or args.animate:
        from storm_dynamics import plotting as splt
        os.makedirs(args.outdir, exist_ok=True)
        print("-" * 72)
        if args.plots:
            print("figure slices -> %s" % splt.plot_rotation_slices(nest, args.outdir, tag="nest"))
        if args.animate:
            print("animation     -> %s" % splt.animate_rotation(nest, args.outdir, tag="nest"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
