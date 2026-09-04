"""DEALIAS the KTLX Moore 2013 velocity and re-derive the validation target (run on WSL2).

Why: the cached extraction (outputs/nexrad_moore/ktlx_velocity.npz) saturates at exactly
+-26.000 m/s in EVERY sub-region, and 9 adjacent-gate pairs differ by >40 m/s (7 in the 45-52
band = 2*Nyquist).  So the recorded "V_rot = 26 m/s" is the NYQUIST CEILING, not a measurement,
and every ratio the study quoted against it was measured against an instrument limit.
Separately, the saved couplet was extracted 4.9 km from the radar -- 21 km away from Moore --
because the old script took the strongest inbound/outbound pair over the whole cropped field.

This script fixes both:
  1. dealias the Doppler sweep (Py-ART region-based) using the radar's own Nyquist velocity;
  2. constrain the couplet search to the MOORE mesocyclone, not the whole field;
  3. report the target WITH its scan geometry (range, beam width, sampling height), so the
     forward operator can be pointed at the right place.

    python3 deploy/wsl2_dealias_moore.py [YYYY MM DD HH MM] [STATION]
"""
import os
import sys
import datetime as dt

import numpy as np

MOORE_LAT, MOORE_LON = 35.34, -97.49          # tornado track near Moore, OK
OUTDIR = ("/mnt/c/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado/"
          "outputs/nexrad_moore")


def couplet(x, y, v, sel, label, radar_lat=None):
    """Strongest inbound/outbound pair within `sel`, with geometry."""
    if sel.sum() < 10:
        print("   %-22s: too few gates (%d)" % (label, sel.sum()))
        return None
    xs, ys, vs = x[sel], y[sel], v[sel]
    i_hi = int(np.nanargmax(vs)); i_lo = int(np.nanargmin(vs))
    sep = float(np.hypot(xs[i_hi] - xs[i_lo], ys[i_hi] - ys[i_lo]))
    rng = 0.5 * (float(np.hypot(xs[i_hi], ys[i_hi])) + float(np.hypot(xs[i_lo], ys[i_lo])))
    dv = float(vs[i_hi] - vs[i_lo])
    print("   %-22s: dV=%6.1f  V_rot=%6.2f  sep=%6.0f m  range=%5.2f km  "
          "(in %+.1f / out %+.1f)"
          % (label, dv, 0.5 * dv, sep, rng / 1000.0, float(vs[i_lo]), float(vs[i_hi])))
    return {"dv": dv, "v_rot": 0.5 * dv, "separation_m": sep, "range_m": rng,
            "inbound": float(vs[i_lo]), "outbound": float(vs[i_hi]),
            "x_hi": float(xs[i_hi]), "y_hi": float(ys[i_hi]),
            "x_lo": float(xs[i_lo]), "y_lo": float(ys[i_lo])}


def main():
    argv = sys.argv[1:]
    y, mo, d, hh, mi = (int(argv[0]), int(argv[1]), int(argv[2]), int(argv[3]), int(argv[4])) \
        if len(argv) >= 5 else (2013, 5, 20, 20, 20)
    station = argv[5] if len(argv) >= 6 else "KTLX"
    os.makedirs(OUTDIR, exist_ok=True)
    tmp = os.path.expanduser("~/nexrad_dl"); os.makedirs(tmp, exist_ok=True)

    import nexradaws
    import pyart
    conn = nexradaws.NexradAwsInterface()
    scans = [s for s in conn.get_avail_scans(y, mo, d, station) if s.scan_time]
    target = dt.datetime(y, mo, d, hh, mi)
    pick = min(scans, key=lambda s: abs((s.scan_time.replace(tzinfo=None) - target).total_seconds()))
    local = os.path.join(tmp, pick.filename)
    if not os.path.exists(local):
        conn.download(pick, tmp)
    fp = local if os.path.exists(local) else conn.download(pick, tmp).success[0].filepath
    radar = pyart.io.read_nexrad_archive(fp)
    print("scan %s  (%s)  sweeps=%d" % (pick.filename, pick.scan_time, radar.nsweeps))

    # radar position and Moore in the radar-relative frame
    rlat = float(radar.latitude["data"][0]); rlon = float(radar.longitude["data"][0])
    import math
    mx = (MOORE_LON - rlon) * 111320.0 * math.cos(math.radians(rlat))
    my = (MOORE_LAT - rlat) * 110540.0
    print("radar at (%.4f, %.4f);  Moore at (%+.0f, %+.0f) m -> range %.1f km"
          % (rlat, rlon, mx, my, math.hypot(mx, my) / 1000.0))

    # lowest sweep that actually carries velocity
    vfield = "velocity" if "velocity" in radar.fields else "VEL"
    sw = None
    for s_ in range(radar.nsweeps):
        v = np.ma.filled(radar.fields[vfield]["data"][radar.get_slice(s_)], np.nan)
        if np.isfinite(v).any():
            sw = s_; break
    sl = radar.get_slice(sw)
    elev = float(np.mean(radar.elevation["data"][sl]))
    nyq = None
    try:
        nyq = float(np.mean(radar.instrument_parameters["nyquist_velocity"]["data"][sl]))
    except Exception:
        pass
    print("velocity sweep %d, elevation %.4f deg, NYQUIST = %s m/s"
          % (sw, elev, ("%.2f" % nyq) if nyq else "unknown"))

    # ---- dealias -------------------------------------------------------------------
    gf = pyart.filters.GateFilter(radar)
    gf.exclude_masked(vfield)
    gf.exclude_invalid(vfield)
    try:
        dealiased = pyart.correct.dealias_region_based(
            radar, vel_field=vfield, gatefilter=gf,
            nyquist_vel=nyq, keep_original=False)
        radar.add_field("velocity_dealiased", dealiased, replace_existing=True)
        ok = True
    except Exception as e:
        print("DEALIASING FAILED: %r" % (e,))
        ok = False

    x = np.asarray(radar.gate_x["data"][sl]).ravel()
    yg = np.asarray(radar.gate_y["data"][sl]).ravel()
    z = np.asarray(radar.gate_altitude["data"][sl]).ravel() - float(radar.altitude["data"][0])
    v_raw = np.ma.filled(radar.fields[vfield]["data"][sl], np.nan).ravel()
    v_dea = (np.ma.filled(radar.fields["velocity_dealiased"]["data"][sl], np.nan).ravel()
             if ok else v_raw)
    good = np.isfinite(v_raw)
    x, yg, z, v_raw, v_dea = x[good], yg[good], z[good], v_raw[good], v_dea[good]
    print("gates with velocity: %d" % good.sum())
    print("  RAW       range: %+.2f .. %+.2f m/s" % (np.nanmin(v_raw), np.nanmax(v_raw)))
    print("  DEALIASED range: %+.2f .. %+.2f m/s" % (np.nanmin(v_dea), np.nanmax(v_dea)))
    if ok and np.nanmax(np.abs(v_dea)) > np.nanmax(np.abs(v_raw)) + 1.0:
        print("  => dealiasing RECOVERED velocities beyond the Nyquist ceiling")

    d_moore = np.hypot(x - mx, yg - my)
    print()
    print("COUPLET, constrained to the Moore mesocyclone (was: whole field, 21 km away):")
    out = {}
    for R in (3000.0, 5000.0, 8000.0):
        print("  within %.0f km of Moore:" % (R / 1000.0))
        out["raw_%dkm" % (R / 1000)] = couplet(x, yg, v_raw, d_moore <= R, "RAW (aliased)")
        out["dea_%dkm" % (R / 1000)] = couplet(x, yg, v_dea, d_moore <= R, "DEALIASED")

    # geometry at the dealiased 5 km couplet
    c = out.get("dea_5km")
    if c:
        rng = c["range_m"]
        beam = rng * math.radians(0.925)
        hh_ = np.interp(rng, np.sort(np.hypot(x, yg)), z[np.argsort(np.hypot(x, yg))])
        print()
        print("CORRECTED TARGET GEOMETRY (dealiased, Moore-constrained):")
        print("   V_rot            = %.2f m/s" % c["v_rot"])
        print("   delta-v          = %.1f m/s" % c["dv"])
        print("   couplet sep      = %.0f m" % c["separation_m"])
        print("   range from KTLX  = %.2f km" % (rng / 1000.0))
        print("   beam diameter    = %.0f m   (0.925 deg)" % beam)
        print("   sampling height  = %.0f m AGL" % hh_)
        print("   elevation        = %.4f deg" % elev)

    np.savez_compressed(os.path.join(OUTDIR, "ktlx_velocity_dealiased.npz"),
                        x_m=x, y_m=yg, z_m=z, vr_raw=v_raw, vr_dealiased=v_dea,
                        elevation_deg=elev, nyquist_m_s=(nyq if nyq else np.nan),
                        scan_time=str(pick.scan_time), radar_lat=rlat, radar_lon=rlon,
                        moore_xy=np.array([mx, my]), dealiased_ok=ok)
    print()
    print("wrote %s/ktlx_velocity_dealiased.npz" % OUTDIR)
    return 0


if __name__ == "__main__":
    sys.exit(main())
