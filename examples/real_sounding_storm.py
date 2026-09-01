"""ROADMAP §3a -- initialise the storm from a REAL observed sounding.

Ingest a radiosonde column (here the bundled illustrative supercell profile, but any observed
sounding works) into the base state via ``soundings.from_observed_sounding``, print the
environment diagnostics (CAPE/CIN/LCL/shear/SRH), then run a short storm on it and report the
rotation -- the first concrete step from the analytic Weisman-Klemp column toward a
real-atmosphere environment.

    python examples/real_sounding_storm.py            # bundled example sounding
    python examples/real_sounding_storm.py --device gpu

To use your OWN sounding, replace ``obs`` below with your radiosonde columns (see
``from_observed_sounding`` for the accepted units / keys).
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np

from meteorological_flow.grid import Grid
from meteorological_flow.base_state import sounding_diagnostics
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics import soundings as snd


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--device", choices=["cpu", "gpu", "auto"], default="cpu")
    p.add_argument("--nx", type=int, default=24)
    p.add_argument("--nz", type=int, default=40)
    p.add_argument("--duration", type=float, default=600.0)
    args = p.parse_args(argv)

    Lx = Ly = 32000.0; Lz = 15000.0; zst = 1.05
    scfg = build_storm_config(preset="storm", nx=args.nx, ny=args.nx, nz=args.nz,
                              Lx=Lx, Ly=Ly, Lz=Lz, duration=args.duration, dt_max=3.0,
                              drag=True, z_stretch=zst, device=args.device)
    g = Grid(nx=args.nx, ny=args.nx, nz=args.nz, Lx=Lx, Ly=Ly, Lz=Lz,
             z_stretch=zst, periodic=True)

    obs = snd.example_observed_sounding()
    base = snd.from_observed_sounding(g, **obs)
    d = sounding_diagnostics(base)
    print("=" * 66)
    print("REAL observed sounding -> storm base state (ROADMAP §3a)")
    print("=" * 66)
    print("  levels ingested : %d  (surface %.0f hPa -> top %.0f hPa)"
          % (len(obs["pressure_hPa"]), obs["pressure_hPa"][0], obs["pressure_hPa"][-1]))
    print("  CAPE            : %6.0f J/kg" % d["CAPE_J_kg"])
    print("  CIN             : %6.0f J/kg" % d["CIN_J_kg"])
    print("  LCL / LFC / EL  : %.0f / %.0f / %.0f m" % (d["LCL_m"], d["LFC_m"], d["EL_m"]))
    print("  0-6 km shear    : %6.1f m/s" % d["shear_0_6km_m_s"])
    print("  0-3 km SRH      : %6.0f m^2/s^2" % snd.storm_relative_helicity(base))
    cx, cy = snd.bunkers_storm_motion(base)
    print("  Bunkers motion  : (%.1f, %.1f) m/s" % (cx, cy))
    print("  parcel w_max    : %6.1f m/s  (sqrt(2 CAPE), thermodynamic ceiling)"
          % d["w_max_parcel_m_s"])

    sim = StormSimulation(scfg, base=base)
    print("\nrunning %d^2 x %d storm for %.0f s (device=%s) ..."
          % (args.nx, args.nz, args.duration, sim.grid.backend.name))
    rep = sim.run(progress=lambda t, dur, s: print("  t=%6.0f/%.0f" % (t, dur), end="\r"))
    print(" " * 40, end="\r")
    rot = rep["rotation"]
    print("  zeta_abs_max    : %.3e s^-1" % rot["zeta_abs_max"])
    print("  updraft w_max   : %.1f m/s" % rot.get("w_max", float("nan")))
    print("  mass residual   : %.2e (normalised)"
          % rep["conservation"]["mass_continuity_residual_norm"])
    print("\nIDEALISED, NOT a forecast: a real sounding sets the ENVIRONMENT, but this is still "
          "an idealised\nwarm-bubble simulation without data assimilation or observational "
          "verification (ROADMAP §3a-e).")


if __name__ == "__main__":
    main()
