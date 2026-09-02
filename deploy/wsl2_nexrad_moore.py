"""Read the REAL NEXRAD KTLX Level II volume of the Moore 2013 tornado (run on WSL2).

Needs Py-ART + nexradaws (pip install --user arm_pyart nexradaws).  Downloads the KTLX scan
nearest the requested time from AWS (no credentials), reads reflectivity + radial velocity, and
reports the observed signature (the velocity couplet = the mesocyclone/TVS) + saves a figure.

    python3 deploy/wsl2_nexrad_moore.py [YYYY MM DD HH MM] [STATION]
"""
import os
import sys
import datetime as dt

import numpy as np


def main():
    argv = sys.argv[1:]
    y, mo, d, hh, mi = (int(argv[0]), int(argv[1]), int(argv[2]), int(argv[3]), int(argv[4])) \
        if len(argv) >= 5 else (2013, 5, 20, 20, 20)
    station = argv[5] if len(argv) >= 6 else "KTLX"
    outdir = "/mnt/c/Users/ileao/OneDrive/Documentos/met_h2o_nucleation_cfd_tornado/outputs/nexrad_moore"
    os.makedirs(outdir, exist_ok=True)
    tmp = os.path.expanduser("~/nexrad_dl"); os.makedirs(tmp, exist_ok=True)

    import nexradaws
    import pyart
    conn = nexradaws.NexradAwsInterface()
    scans = [s for s in conn.get_avail_scans(y, mo, d, station) if s.scan_time]
    if not scans:
        print("no scans for %s on %04d-%02d-%02d" % (station, y, mo, d)); return 1
    target = dt.datetime(y, mo, d, hh, mi)
    pick = min(scans, key=lambda s: abs((s.scan_time.replace(tzinfo=None) - target).total_seconds()))
    print("nearest scan to %sZ: %s (%s)" % (target, pick.filename, pick.scan_time))
    res = conn.download(pick, tmp)
    fp = res.success[0].filepath
    radar = pyart.io.read_nexrad_archive(fp)
    print("REAL KTLX volume read: %d sweeps, fields=%s" % (radar.nsweeps, list(radar.fields)))

    s = radar.get_slice(0)
    refl = np.ma.filled(radar.fields["reflectivity"]["data"][s], np.nan)
    vel = None
    for nm in ("velocity", "VEL"):
        if nm in radar.fields:
            # velocity is usually on a separate (Doppler) sweep; find the first one that has it
            for sw in range(radar.nsweeps):
                ss = radar.get_slice(sw)
                v = np.ma.filled(radar.fields[nm]["data"][ss], np.nan)
                if np.isfinite(v).any():
                    vel = v; break
            break
    print("OBSERVED reflectivity: max=%.1f dBZ" % np.nanmax(refl))
    if vel is not None:
        vmax, vmin = np.nanmax(vel), np.nanmin(vel)
        print("OBSERVED radial velocity couplet: %.1f (outbound) / %.1f (inbound) m/s -> Δ=%.1f "
              "m/s (the mesocyclone/TVS signature)" % (vmax, vmin, vmax - vmin))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        disp = pyart.graph.RadarDisplay(radar)
        fig = plt.figure(figsize=(13, 6))
        ax1 = fig.add_subplot(121)
        disp.plot("reflectivity", 0, ax=ax1, vmin=-8, vmax=64, cmap="NWSRef",
                  title="KTLX reflectivity  %s" % pick.scan_time)
        disp.set_limits((-150, 150), (-150, 150), ax=ax1)
        if vel is not None:
            ax2 = fig.add_subplot(122)
            vsweep = next((sw for sw in range(radar.nsweeps)
                           if "velocity" in radar.fields and
                           np.isfinite(np.ma.filled(radar.fields["velocity"]["data"][radar.get_slice(sw)], np.nan)).any()), 1)
            disp.plot("velocity", vsweep, ax=ax2, vmin=-40, vmax=40, cmap="NWSVel",
                      title="KTLX radial velocity (couplet = rotation)")
            disp.set_limits((-150, 150), (-150, 150), ax=ax2)
        out = os.path.join(outdir, "ktlx_moore_%04d%02d%02d_%02d%02d.png" % (y, mo, d, hh, mi))
        fig.tight_layout(); fig.savefig(out, dpi=120)
        print("figure ->", out)
    except Exception as e:
        print("figure skipped:", e)
    print("\nREAL radar of the Moore tornado ingested. Compare synthetic V_r from the CFD with "
          "atmospheric_data.radial (same beam geometry) for validation.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
