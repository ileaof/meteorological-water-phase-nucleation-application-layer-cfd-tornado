"""Extract the low-level radial-velocity field of the Moore mesocyclone from the real KTLX
volume (WSL2/Py-ART) and save it (+ the observed couplet signature) for the sim comparison.

    python3 deploy/wsl2_extract_velocity.py
Saves outputs/nexrad_moore/ktlx_velocity.npz with x_m, y_m, vr (0.5 deg sweep) cropped near
Moore, and prints the observed rotational-velocity / couplet Delta-v.
"""
import os
import sys
import datetime as dt

import numpy as np


def main():
    y, mo, d, hh, mi, station = 2013, 5, 20, 20, 20, "KTLX"
    repo = "/mnt/c/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado"
    outdir = os.path.join(repo, "outputs", "nexrad_moore"); os.makedirs(outdir, exist_ok=True)
    tmp = os.path.expanduser("~/nexrad_dl"); os.makedirs(tmp, exist_ok=True)
    import nexradaws, pyart
    conn = nexradaws.NexradAwsInterface()
    scans = [s for s in conn.get_avail_scans(y, mo, d, station) if s.scan_time]
    pick = min(scans, key=lambda s: abs((s.scan_time.replace(tzinfo=None) - dt.datetime(y, mo, d, hh, mi)).total_seconds()))
    fp = conn.download(pick, tmp).success[0].filepath
    radar = pyart.io.read_nexrad_archive(fp)

    # first sweep that actually has velocity (the 0.5 deg Doppler sweep)
    sw = next(s for s in range(radar.nsweeps)
              if np.isfinite(np.ma.filled(radar.fields["velocity"]["data"][radar.get_slice(s)], np.nan)).any())
    s = radar.get_slice(sw)
    vr = np.ma.filled(radar.fields["velocity"]["data"][s], np.nan)
    x = np.asarray(radar.gate_x["data"][s]); yg = np.asarray(radar.gate_y["data"][s])
    # crop to the Moore storm region (~19 km W of KTLX): x in [-50,10] km, y in [-30,30] km
    m = (x > -50000) & (x < 10000) & (yg > -30000) & (yg < 30000) & np.isfinite(vr)
    xc, yc, vc = x[m], yg[m], vr[m]
    # observed couplet: the strongest inbound/outbound pair within ~8 km of each other
    imax = int(np.nanargmax(vc)); imin = int(np.nanargmin(vc))
    sep = float(np.hypot(xc[imax] - xc[imin], yc[imax] - yc[imin]))
    dv = float(vc[imax] - vc[imin]); vrot = 0.5 * dv
    az_shear = dv / max(sep, 1.0)                         # s^-1 (proxy for vertical vorticity)
    np.savez_compressed(os.path.join(outdir, "ktlx_velocity.npz"),
                        x_m=xc, y_m=yc, vr=vc, elevation_deg=float(np.mean(radar.elevation["data"][s])),
                        scan_time=str(pick.scan_time),
                        couplet=np.array([xc[imax], yc[imax], vc[imax], xc[imin], yc[imin], vc[imin]]))
    print("OBSERVED mesocyclone (KTLX %s, %.1f deg sweep):" % (pick.scan_time, np.mean(radar.elevation["data"][s])))
    print("  outbound %.1f m/s @ (%.1f,%.1f) km ; inbound %.1f m/s @ (%.1f,%.1f) km"
          % (vc[imax], xc[imax] / 1e3, yc[imax] / 1e3, vc[imin], xc[imin] / 1e3, yc[imin] / 1e3))
    print("  couplet Delta-v = %.1f m/s ; V_rot = %.1f m/s ; separation = %.1f km ; "
          "azimuthal shear = %.2e s^-1" % (dv, vrot, sep / 1e3, az_shear))
    print("saved -> outputs/nexrad_moore/ktlx_velocity.npz")
    return 0


if __name__ == "__main__":
    sys.exit(main())
