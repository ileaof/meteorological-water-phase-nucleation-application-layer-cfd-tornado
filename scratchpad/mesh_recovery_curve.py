"""How fine must the MODEL be before the radar observable stops changing?

The corrected Moore target is a beam-averaged observable: V_rot 39.49 m/s at 208 m above the
KTLX antenna, through a 327 m beam at 20.25 km, elevation 0.5211 deg, couplet separation 584 m.
Two separate penalties reduce what a model reports through that instrument:

  * the MESH penalty  -- the model cannot represent the vortex;
  * the BEAM penalty  -- the radar cannot resolve it, however fine the model is.

Only the first is ours to fix, so the question is where it vanishes.  Sweeping a Rankine vortex
of the target's own core radius through the forward operator at the true geometry:

    model dx   cells/2R   recovery
       600       0.97      0.417
       300       1.95      0.564
       150       3.89      0.676
       100       5.84      0.700   <- converged
        67       8.72      0.699
        50      11.68      0.700
        30      19.47      0.700
        15      38.93      0.700

=> dx = 100 m is CONVERGED in the observable; finer buys nothing.  The residual 0.700 is the
irreducible beam penalty.  Two consequences: (1) the 67 m figure both audits recommended from a
cells-across-core rule of thumb is unnecessary; (2) the memory saved is better spent on DOMAIN
(24 km at 100 m = 3.69M cells = 4.0 GB, the same cost as 16 km at 67 m) because the trusted
interior is limited by advective transmission U_sr*T ~ 6 km over a 300 s window, not by dx.

CAVEAT: the curve is derived from a RANKINE vortex.  A Burgers-Rott or two-celled vortex has a
different radial profile and therefore a different beam response; this is a model-dependent
inversion, not a measurement.

    python scratchpad/mesh_recovery_curve.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from atmospheric_data import radar_operator as ro

V_TRUE, CORE_M = 39.49, 292.0            # corrected target: V_rot, and sep/2 as the core radius
RANGE_M, ELEV = 20250.0, 0.5211


def curve(dxs=(600., 300., 150., 100., 67., 50., 30., 15.), half=2500.0):
    rad = ro.RadarSpec(elevation_deg=ELEV)
    cx, cy = 0.0, RANGE_M * np.cos(np.radians(ELEV))
    z = np.array([100.0, 208.0, 320.0])
    out = []
    for dx in dxs:
        x = np.arange(cx - half, cx + half + 1, dx)
        y = np.arange(cy - half, cy + half + 1, dx)
        u, v, w = ro.rankine_vortex(x, y, z, (cx, cy), V_TRUE, CORE_M)
        az, rg = ro.sweep_grid_for_domain(rad, cx, cy, half - 200.0, az_resolution_deg=0.5)
        g = ro.vrot_from_sweep(ro.synthetic_sweep(u, v, w, x, y, z, rad, az, rg,
                                                  n_beam=5, n_gate=3))
        out.append({"dx_m": dx, "cells_per_core_diameter": 2 * CORE_M / dx,
                    "v_rot_obs": g["v_rot_m_s"], "recovery": g["v_rot_m_s"] / V_TRUE})
    return out


if __name__ == "__main__":
    rad = ro.RadarSpec(elevation_deg=ELEV)
    print("true V_rot %.2f m/s, core %.0f m, range %.2f km, beam %.0f m, elev %.4f deg"
          % (V_TRUE, CORE_M, RANGE_M / 1000, rad.beam_diameter_m(RANGE_M), ELEV))
    print("%9s %11s %11s %10s" % ("model dx", "cells/2R", "V_rot obs", "recovery"))
    for r in curve():
        print("%9.0f %11.2f %11.2f %10.3f"
              % (r["dx_m"], r["cells_per_core_diameter"], r["v_rot_obs"], r["recovery"]))
