"""ROADMAP §3e -- the Straka et al. (1993) density-current benchmark.

A cold bubble (Delta_theta = -15 K) collapses in a neutral, dry, non-rotating atmosphere at a
fixed viscosity nu = 75 m^2/s.  This is a standard code-verification problem: the leading edge
of the surface cold outflow reaches ~15.5 km from the centre at 900 s (100 m reference), with
three Kelvin-Helmholtz rotors and peak velocities ~15-20 m/s.

    python examples/straka_benchmark.py                 # 200 m, ~30 s CPU
    python examples/straka_benchmark.py --ref --device gpu   # 100 m reference
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np

from storm_dynamics.benchmarks import straka_simulation, straka_front_position


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref", action="store_true", help="100 m reference (nz=64) instead of 200 m")
    p.add_argument("--device", choices=["cpu", "gpu", "auto"], default="cpu")
    p.add_argument("--duration", type=float, default=900.0)
    args = p.parse_args(argv)
    nz = 64 if args.ref else 32
    dx = 51200.0 / (8 * nz)

    print("=" * 68)
    print("Straka (1993) density current -- %.0f m resolution, nu=75 m^2/s, %.0f s"
          % (dx, args.duration))
    print("=" * 68)
    sim = straka_simulation(nz=nz, ny=4, duration=args.duration, device=args.device)
    sim.run(progress=lambda t, d, s: print("  t=%4.0f/%.0f s" % (t, d), end="\r"))
    print(" " * 30, end="\r")

    to = sim.grid.backend.to_cpu
    dth = np.asarray(to(sim.state.theta - sim.theta0_field))
    w = np.asarray(to(sim.state.w))
    x = np.asarray(to(sim.grid.xc)); xc = 0.5 * float(sim.grid.Lx)
    front = straka_front_position(sim)
    cold = np.where(dth[:, 0, 0] < -1.0)[0]
    left = (xc - float(x[cold].min())) / 1000.0 if cold.size else 0.0
    right = (float(x[cold].max()) - xc) / 1000.0 if cold.size else 0.0
    print("  front (from centre) : %.2f km   (reference ~15.5 km at 100 m, 900 s)" % (front / 1000.0))
    print("  symmetry  L / R     : %.2f / %.2f km" % (left, right))
    print("  cold-pool min theta': %.1f K   (started at -15 K, diffused by nu)" % dth.min())
    print("  vertical velocity   : w in [%.1f, %.1f] m/s (KH rotors ~15-20)" % (w.min(), w.max()))
    print("\nStandard idealised BENCHMARK (code verification), not a forecast.")


if __name__ == "__main__":
    main()
