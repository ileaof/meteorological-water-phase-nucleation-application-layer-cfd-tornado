"""ROUND-TRIP FOLD TEST -- put an error bar on the dealiased Moore target (run on WSL2).

The corrected target V_rot = 39.49 m/s rests on Py-ART's region-based dealiasing of a field whose
Nyquist velocity is 26.12 m/s.  The audit verified that unfold is internally consistent (integer
folds only, |k| <= 1, adjacent-gate discontinuities 1051 -> 10, an independent BFS agreeing on
99.745% of gates) but flagged six gates at the vortex core where two algorithms disagree, which
is precisely where the couplet is defined.

Internal consistency is not an error bar.  This script makes one: build a synthetic sweep from an
analytic vortex whose V_rot is KNOWN, fold it at the observation's own Nyquist, run the IDENTICAL
`pyart.correct.dealias_region_based` call, and measure what is recovered.

    truth  ->  fold at +-26.12  ->  Py-ART dealias  ->  measured V_rot
                                                        vs known V_rot

The difference is the dealiasing's contribution to the target's uncertainty.  Anything the
algorithm cannot recover here it also could not recover from the real sweep.

    python3 deploy/wsl2_fold_roundtrip.py
"""
import os
import sys

import numpy as np

OUTDIR = ("/mnt/c/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado/"
          "outputs/nexrad_moore")
REPO = "/mnt/c/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"

NYQUIST = 26.12          # the real sweep's own Nyquist velocity
ELEV = 0.5211
RANGE_M = 20250.0        # range to the Moore couplet
GATE_M = 250.0
AZ_RES = 0.5


def fold(v, nyq):
    """Wrap velocities into [-nyq, +nyq) exactly as a pulse-pair estimator would."""
    return ((np.asarray(v) + nyq) % (2.0 * nyq)) - nyq


def build(v_true_rot, core_m, ambient_u=-2.0, ambient_v=10.0, ngates=180, nrays=360):
    """A PPI sweep containing a Rankine vortex of known V_rot plus a uniform ambient wind."""
    import pyart
    radar = pyart.testing.make_empty_ppi_radar(ngates, nrays, 1)
    radar.range["data"] = np.arange(ngates, dtype="float64") * GATE_M + 2000.0
    radar.azimuth["data"] = np.arange(nrays, dtype="float64") * (360.0 / nrays)
    radar.elevation["data"] = np.full(nrays, ELEV)
    radar.fixed_angle["data"] = np.array([ELEV])
    radar.latitude["data"] = np.array([35.3331]); radar.longitude["data"] = np.array([-97.2775])
    radar.altitude["data"] = np.array([370.0])
    radar.init_gate_x_y_z()

    gx = np.asarray(radar.gate_x["data"]); gy = np.asarray(radar.gate_y["data"])
    gz = np.asarray(radar.gate_z["data"])
    # vortex centre due north at the observation range
    cx, cy = 0.0, RANGE_M * np.cos(np.radians(ELEV))
    dx = gx - cx; dy = gy - cy
    r = np.sqrt(dx * dx + dy * dy) + 1e-9
    R = float(core_m)
    vth = np.where(r <= R, v_true_rot * r / R, v_true_rot * R / r)
    u = ambient_u - vth * dy / r
    v = ambient_v + vth * dx / r
    rr = np.sqrt(gx * gx + gy * gy + gz * gz) + 1e-9
    v_r = (u * gx + v * gy) / rr                       # w = 0; horizontal wind only
    return radar, v_r, (cx, cy)


def couplet(gx, gy, vr, centre, radius_m=3000.0):
    d = np.hypot(gx - centre[0], gy - centre[1])
    sel = d <= radius_m
    if sel.sum() < 10:
        return None
    xs, ys, vs = gx[sel], gy[sel], vr[sel]
    i_hi = int(np.nanargmax(vs)); i_lo = int(np.nanargmin(vs))
    return {"v_rot": 0.5 * float(vs[i_hi] - vs[i_lo]),
            "sep": float(np.hypot(xs[i_hi] - xs[i_lo], ys[i_hi] - ys[i_lo])),
            "hi": float(vs[i_hi]), "lo": float(vs[i_lo])}


def main():
    import pyart
    print("ROUND-TRIP FOLD TEST  (Nyquist %.2f m/s, elev %.4f deg, range %.2f km)"
          % (NYQUIST, ELEV, RANGE_M / 1000))
    print("%10s %9s | %10s %10s %10s | %9s %9s"
          % ("true Vrot", "core", "folded", "dealiased", "recovered", "err", "sep err"))
    rows = []
    for v_true, core in ((39.49, 292.0), (39.49, 584.0), (30.0, 292.0),
                         (45.0, 292.0), (60.0, 292.0)):
        radar, v_r, centre = build(v_true, core)
        gx = np.asarray(radar.gate_x["data"]).ravel()
        gy = np.asarray(radar.gate_y["data"]).ravel()
        v_folded = fold(v_r, NYQUIST)
        n_folded = int(np.sum(np.abs(v_r) > NYQUIST))

        radar.add_field("velocity", {"data": np.ma.array(v_folded),
                                     "units": "m/s", "long_name": "folded"},
                        replace_existing=True)
        radar.instrument_parameters = {
            "nyquist_velocity": {"data": np.full(radar.nrays, NYQUIST),
                                 "units": "m/s", "long_name": "nyquist"}}
        gf = pyart.filters.GateFilter(radar)
        gf.exclude_masked("velocity")
        dz = pyart.correct.dealias_region_based(radar, vel_field="velocity", gatefilter=gf,
                                                nyquist_vel=NYQUIST, keep_original=False)
        v_dea = np.ma.filled(dz["data"], np.nan).ravel()

        c_true = couplet(gx, gy, v_r.ravel(), centre)
        c_fold = couplet(gx, gy, v_folded.ravel(), centre)
        c_dea = couplet(gx, gy, v_dea, centre)
        rows.append((v_true, core, c_true, c_fold, c_dea, n_folded))
        print("%10.2f %9.0f | %10.2f %10.2f %10.2f | %+8.2f %+8.0f   (%d gates folded)"
              % (v_true, core, c_fold["v_rot"], c_dea["v_rot"], c_dea["v_rot"],
                 c_dea["v_rot"] - c_true["v_rot"], c_dea["sep"] - c_true["sep"], n_folded))

    print()
    print("truth here is the OPERATOR-SAMPLED field before folding, so this isolates the")
    print("dealiasing alone -- beam and mesh penalties are NOT part of this error bar.")
    errs = [abs(r[4]["v_rot"] - r[2]["v_rot"]) for r in rows]
    print("dealiasing error on V_rot: max %.2f m/s, mean %.2f m/s over %d cases"
          % (max(errs), float(np.mean(errs)), len(errs)))
    np.savez_compressed(os.path.join(OUTDIR, "fold_roundtrip.npz"),
                        v_true=[r[0] for r in rows], core_m=[r[1] for r in rows],
                        vrot_true=[r[2]["v_rot"] for r in rows],
                        vrot_folded=[r[3]["v_rot"] for r in rows],
                        vrot_dealiased=[r[4]["v_rot"] for r in rows],
                        nyquist=NYQUIST)
    print("wrote %s/fold_roundtrip.npz" % OUTDIR)
    return 0


if __name__ == "__main__":
    sys.exit(main())
