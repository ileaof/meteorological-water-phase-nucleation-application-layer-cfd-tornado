"""Idealised rotating supercell / tornadogenesis demonstration (storm_dynamics).

Runs the rotating deep-convection core on a demonstration grid and prints the
rotation diagnostics -- vertical vorticity (zeta), updraft helicity, the
mid-level mesocyclone and near-surface zeta trackers, and the environmental
storm-relative helicity (SRH) / bulk shear of the hodograph.

Two scenarios (choose with --scenario):

* ``supercell`` (M1): unidirectional shear -> a warm bubble splits into left- and
  right-moving cells with a mid-level mesocyclone (Klemp-Wilhelmson /
  Weisman-Klemp storm-splitting reference).
* ``tornadogenesis`` (M2): curved (quarter-circle) hodograph + surface drag +
  evaporative cold pool -> near-surface vertical vorticity (the tornadogenesis
  proxy) on the forward-flank / cold-pool interface.

IDEALISED, NOT a forecast: no data assimilation, no real event, no observational
verification (see docs/storm_dynamics_guide.md).

Usage:
    PYTHONPATH=src python examples/supercell_tornadogenesis.py --scenario supercell
    PYTHONPATH=src python examples/supercell_tornadogenesis.py --scenario tornadogenesis \
        --nx 32 --ny 32 --nz 40 --duration 2400
"""
from __future__ import annotations

import argparse
import os
import sys

# run directly (`python examples/supercell_tornadogenesis.py`) without needing
# PYTHONPATH=src: put the repo's src/ on the path, like the other examples do.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from storm_dynamics.config import build_storm_config   # noqa: E402
from storm_dynamics.core import StormSimulation         # noqa: E402
from meteorological_flow.backend import get_backend      # noqa: E402


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scenario", choices=["supercell", "tornadogenesis"],
                   default="supercell")
    p.add_argument("--nx", type=int, default=32)
    p.add_argument("--ny", type=int, default=32)
    p.add_argument("--nz", type=int, default=40)
    p.add_argument("--Lx", type=float, default=40000.0)
    p.add_argument("--Ly", type=float, default=40000.0)
    p.add_argument("--Lz", type=float, default=16000.0)
    p.add_argument("--duration", type=float, default=2400.0)
    p.add_argument("--dt-max", type=float, default=3.0)
    p.add_argument("--kernel-nucleation", action="store_true",
                   help="couple the validated nucleation kernel as the microphysics "
                        "embryo source (builds a lookup table; slower)")
    p.add_argument("--device", choices=["auto", "cpu", "gpu"], default="cpu",
                   help="compute backend (default cpu, matching the prior "
                        "hardcoded behaviour); auto uses GPU when available "
                        "else CPU, gpu fails loudly if unavailable")
    p.add_argument("--plots", action="store_true",
                   help="write rotation figures (hodograph, slices, time series) as PNGs")
    p.add_argument("--outdir", default=None,
                   help="output directory for figures (default outputs/storm_<scenario>)")
    args = p.parse_args(argv)

    tornado = args.scenario == "tornadogenesis"
    hodo = "quarter_circle" if tornado else "unidirectional"
    # the tornadogenesis scenario clusters vertical levels near the surface and
    # moderates the hodograph so the drag / corner-flow layer is resolved and the
    # updraft stays physical (a uniform coarse dz blows up under strong low-level
    # shear + surface drag); the supercell (M1) scenario uses a uniform grid.
    scfg = build_storm_config(
        preset="storm", nx=args.nx, ny=args.ny, nz=args.nz,
        Lx=args.Lx, Ly=args.Ly, Lz=args.Lz,
        duration=args.duration, dt_max=args.dt_max,
        hodograph_kind=hodo, drag=tornado,
        z_stretch=1.05 if tornado else None,
        U_max=18.0 if tornado else None,
        z_turn=2000.0 if tornado else None,
        C_s=0.22 if tornado else None,
        couple_nucleation=args.kernel_nucleation)

    if args.kernel_nucleation:
        cache = os.path.join(scfg.sim.output.outdir, "nucleation_lookup.npz")
        if not os.path.exists(cache):
            print("[kernel] building the nucleation lookup table (one-time, slow -- "
                  "several minutes); it is cached at %s for later runs..." % cache,
                  flush=True)
    backend = get_backend(args.device)
    sim = StormSimulation(scfg, backend=backend)
    from storm_dynamics import soundings as snd
    print("=" * 72)
    print("storm_dynamics -- %s scenario (IDEALISED, not a forecast)" % args.scenario)
    print("=" * 72)
    print("backend         : %s%s" % (backend.name,
          "  (%s)" % backend.fallback_reason if backend.fallback_reason else ""))
    print("grid            : %d x %d x %d  (dx=%.0f m, dz=%.0f m)"
          % (sim.grid.nx, sim.grid.ny, sim.grid.nz, sim.grid.dx, sim.grid.dz))
    print("hodograph       : %s" % hodo)
    print("f (Coriolis)    : %.3e 1/s  (lat %.1f N)" % (sim.f, scfg.dyn.latitude_deg))
    print("kernel coupling : %s" % ("ON (validated nucleation kernel -> embryo source)"
                                    if sim.couple_nucleation else "off"))
    print("SRH 0-1 / 0-3 km: %.0f / %.0f m^2/s^2"
          % (snd.storm_relative_helicity(sim.base, 1000.0),
             snd.storm_relative_helicity(sim.base, 3000.0)))
    print("shear 0-1/0-6 km: %.1f / %.1f m/s"
          % (snd.bulk_shear(sim.base, 0, 1000), snd.bulk_shear(sim.base, 0, 6000)))
    print("-" * 72)

    def progress(t, dur, step):
        print("  t=%6.0f s / %.0f  (step %d)" % (t, dur, step), end="\r")

    report = sim.run(progress=progress)
    print(" " * 72, end="\r")

    rp = report["rotation"]
    pk = report["rotation_peak"]
    cons = report["conservation"]
    print("final time      : %.0f s  (%d steps, %.1f s wall)"
          % (report["final_time"], report["n_steps"], report["wall_clock_s"]))
    print("w_max / w_min   : %+.1f / %+.1f m/s" % (rp["w_max"], rp["w_min"]))
    print("zeta max / min  : %+.2e / %+.2e 1/s  (storm splitting = both signs)"
          % (rp["zeta_max"], rp["zeta_min"]))
    print("mid-level meso  : %.2e 1/s   (peak %.2e)"
          % (rp["midlevel_mesocyclone"], pk["peak_midlevel_mesocyclone"]))
    print("near-surface zeta: %.2e 1/s  (peak %.2e)  <- tornadogenesis proxy"
          % (rp["near_surface_zeta_max"], pk["peak_near_surface_zeta"]))
    print("updraft helicity: %.1f m^2/s^2  (peak %.1f)"
          % (rp["updraft_helicity_max"], pk["peak_updraft_helicity"]))
    print("-" * 72)
    print("conservation    : water rel-err %.2e | energy rel-err %.2e | "
          "mass-continuity |div| %.2e"
          % (cons["total_water_rel_err"], cons["total_energy_rel_err"],
             cons["mass_continuity_residual_abs"]))
    print("=" * 72)
    for lim in report["limitations"]:
        print("NOTE:", lim)

    if args.plots:
        from storm_dynamics import plotting as splt
        outdir = args.outdir or ("outputs/storm_%s" % args.scenario)
        paths = splt.plot_all(sim, outdir, tag=args.scenario)
        print("-" * 72)
        for name, path in paths.items():
            print("figure %-11s -> %s" % (name, path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
