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
    p.add_argument("--regrid-interval", type=int, default=None,
                   help="M3/§2a adaptive regridding: every N parent steps, re-centre the "
                        "fixed-size nest on the tagged vortex (data-driven follow, ground "
                        "frame; not combinable with --follow)")
    p.add_argument("--concurrent", action="store_true",
                   help="M3 phase 2: step the parent alongside the nest so the nest "
                        "border sees time-evolving parent boundaries (sustains the "
                        "storm beyond the frozen-boundary short window)")
    p.add_argument("--follow", action="store_true",
                   help="M3 phase 2b: storm-following nest (storm-relative frame) so "
                        "the cell stays centred and the vortex is sustained + "
                        "intensified over a long window (implies --concurrent)")
    p.add_argument("--two-way", action="store_true",
                   help="M3 phase 3a: approximate two-way feedback -- blend the nest's "
                        "finer solution back onto the parent overlap (implies --concurrent)")
    p.add_argument("--composite", action="store_true",
                   help="docs/ROADMAP.md section 1: replace the two per-level anelastic "
                        "projections with ONE composite parent+nest mass-flux solve so "
                        "div(rho0 u)=0 holds across the coarse-fine interface every "
                        "sub-step (implies --concurrent; the nest footprint is snapped "
                        "to parent cells / matched-z)")
    p.add_argument("--device", choices=["cpu", "gpu", "auto"], default="cpu")
    p.add_argument("--u-max", type=float, default=18.0,
                   help="hodograph magnitude [m/s] (larger = stronger shear -> stronger storm)")
    p.add_argument("--les-boost", type=float, default=1.25,
                   help="nest SGS dissipation factor (raise for stability at high --refine)")
    p.add_argument("--cfl", type=float, default=0.25, help="nest CFL cap (lower = more stable)")
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
        drag=True, z_stretch=1.05, U_max=args.u_max, z_turn=2000.0, C_s=0.22, device=args.device)
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
    if args.composite or args.regrid_interval:
        # composite projection AND adaptive regridding need a cell-aligned, matched-z
        # nest: snap the footprint to parent cells + use the parent's vertical grid
        i0 = int(round(spec.x0 / parent.grid.dx))
        j0 = int(round(spec.y0 / parent.grid.dy))
        ncx = int(round(spec.Lx / parent.grid.dx))
        ncy = int(round(spec.Ly / parent.grid.dy))
        if args.nest_nz not in (None, parent.grid.nz):
            print("note        : --composite snaps the nest to the parent's vertical "
                  "grid (matched z): nest nz %d -> %d" % (args.nest_nz, parent.grid.nz))
        spec = nst.NestSpec.aligned(parent.grid, i0=i0, j0=j0, ncx=ncx, ncy=ncy,
                                    refine=args.refine)
    concurrent = args.concurrent or args.follow or args.two_way or args.composite or args.regrid_interval
    mode = ("phase 3a (two-way feedback" + (", storm-following)" if args.follow else ")") if args.two_way else
            "phase 2b (storm-following, concurrent)" if args.follow else
            "phase 2 (concurrent, time-evolving parent boundary)" if concurrent else
            "phase 1 (frozen parent boundary)")
    if args.composite:
        mode += " + composite two-level projection"
    print("nest region   : centred (%.1f, %.1f) km, half %.0f km" % (xc / 1000, yc / 1000, args.half / 1000))
    print("mode          : %s" % mode)
    if concurrent:
        nest = nst.NestedStormSimulation(parent, spec)
        r_nest0 = rot.rotation_report(nest.state, nest.grid)
        z0_int = nst.interior_near_surface_zeta(nest)
        print("nest grid     : %d x %d x %d  dx=%.0f m (parent %.0f m, %dx finer), dz0=%.0f m"
              % (nest.grid.nx, nest.grid.ny, nest.grid.nz, nest.grid.dx, parent.grid.dx,
                 args.refine, float(nest.grid.dz_c[0])))
        print("-" * 72)
        nest, rep = nst.run_concurrent_nest(
            parent, spec, window=args.window, capture_frames=args.animate,
            follow=args.follow, two_way=args.two_way, composite_projection=args.composite,
            les_boost=args.les_boost, cfl=args.cfl, regrid_interval=args.regrid_interval,
            progress=lambda t, d, s: print("  nest   t=%6.0f/%.0f" % (t, d), end="\r"))
    else:
        nest = nst.NestedStormSimulation(parent, spec)
        nest.cfg.time.duration = args.window
        r_nest0 = rot.rotation_report(nest.state, nest.grid)
        z0_int = nst.interior_near_surface_zeta(nest)
        print("nest grid     : %d x %d x %d  dx=%.0f m (parent %.0f m, %dx finer), dz0=%.0f m"
              % (nest.grid.nx, nest.grid.ny, nest.grid.nz, nest.grid.dx, parent.grid.dx,
                 args.refine, float(nest.grid.dz_c[0])))
        print("-" * 72)
        rep = nest.run(progress=lambda t, d, s: print("  nest   t=%6.0f/%.0f" % (t, d), end="\r"),
                       capture_frames=args.animate)
    print(" " * 60, end="\r")

    r = rep["rotation"]; pk = rep["rotation_peak"]; c = rep["conservation"]
    z_int = nst.interior_near_surface_zeta(nest)
    print("near-surface zeta [s^-1]  (interior = physical vortex, excludes the sponge band):")
    print("  parent (coarse, dx=%.0f m)        : %.2e" % (parent.grid.dx, r_par["near_surface_zeta_max"]))
    print("  nest initial (interior)            : %.2e" % z0_int)
    print("  nest final   (interior)            : %.2e" % z_int)
    print("  nest final   (incl. sponge edge)   : %.2e" % r["near_surface_zeta_max"])
    intens = z_int / max(z0_int, 1e-12)
    print("  -> interior intensification        : %.2fx" % intens)
    print("w_max: parent %.1f -> nest init %.1f -> nest final %.1f m/s"
          % (r_par["w_max"], r_nest0["w_max"], r["w_max"]))
    print("conservation (nest): water %.2e | mass-continuity |div| %.2e"
          % (c["total_water_rel_err"], c["mass_continuity_residual_norm"]))
    print("=" * 72)
    for lim in rep["limitations"]:
        print("NOTE:", lim)
    if args.two_way:
        print("NOTE: approximate TWO-WAY feedback (M3 phase 3a): the nest's finer "
              "solution is blended back onto the parent overlap each parent step, so "
              "the parent is improved by the nest (a closed pai<->nest loop). This is "
              "*injection* feedback, NOT rigorous flux-conservative refluxing -- strict "
              "interface conservation (Berger-Colella refluxing + multilevel Poisson) is "
              "the full-AMR project. Adaptive (dynamic) refinement also remains future work.")
    elif args.follow:
        print("NOTE: storm-FOLLOWING nest (M3 phase 2b): the nest runs in the "
              "storm-relative frame (motion C=%s m/s) so the cell stays centred; the "
              "vortex is sustained and INTENSIFIED over the window (vs the fixed nest, "
              "which decays). Still one-way + fixed refinement; higher refinement "
              "(O(10-100 m)) and two-way/adaptive (AMR) nesting remain future work."
              % ([round(v, 1) for v in rep["nest"]["storm_motion"]]))
    elif concurrent:
        print("NOTE: one-way CONCURRENT nest (M3 phase 2): the parent steps alongside "
              "the nest and feeds time-evolving boundaries, so the storm is sustained "
              "beyond the frozen-boundary window. Still one-way + fixed refinement; a "
              "storm-following MOVING nest, higher refinement and two-way/adaptive "
              "(AMR) nesting remain future work.")
    else:
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
