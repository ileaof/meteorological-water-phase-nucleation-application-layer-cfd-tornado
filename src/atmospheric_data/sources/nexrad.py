"""NEXRAD Level II ingestion (ROADMAP §3a, source C): radar volumes for VALIDATION.

Reads Level II via ``Py-ART`` (or ``xradar``) into the same volume dict the synthetic radar
produces (gate lat/lon/alt, azimuth/elevation, reflectivity, radial velocity, ...).  Download
is from the **AWS Open Data** bucket ``noaa-nexrad-level2`` (no credentials).  Requires
``arm_pyart`` (optional) -- absent -> ``SourceUnavailable``.

CRUCIAL (task): radial velocity is ``V_r = V . r_hat``, NOT the 3-D vector.  This module only
*reads* observations; the synthetic-radial operator that projects the CFD (u,v,w) onto the beam
is in :mod:`atmospheric_data.radial`.
"""
from __future__ import annotations

import os

import numpy as np

from .base import require, SourceUnavailable, optional_import


def download(cfg, cache, when=None):
    """Fetch the nearest Level II volume for the station/time from the AWS archive.

    The ``noaa-nexrad-level2`` bucket does NOT allow anonymous listing, so we use **nexradaws**
    (``pip install nexradaws``) to query the archive and pick the scan nearest the case time --
    the standard tool.  ``requests`` is the fallback only when the exact object key is known.
    No credentials needed either way."""
    st = cfg.data.radar_station
    date = cfg.case.date.replace("-", "")
    hh = cfg.case.start_time_utc.split(":")[0].zfill(2)
    key = "%s_%s_%s" % (st, date, hh)
    path = cache.path("nexrad", key, ".ar2v")
    cache.require_offline_ok("nexrad", key, ".ar2v")
    if cache.has("nexrad", key, ".ar2v"):
        return path
    nexradaws = optional_import("nexradaws")
    if nexradaws is None:
        raise SourceUnavailable(
            "listing the NEXRAD archive needs 'nexradaws' (pip install nexradaws); the bucket "
            "does not allow anonymous listing. Install it (works on WSL2) or place the .ar2v "
            "file in the cache and run --offline.")
    import datetime as dt
    conn = nexradaws.NexradAwsInterface()
    y, mo, d = int(date[:4]), int(date[4:6]), int(date[6:8])
    start = dt.datetime(y, mo, d, int(hh), 0)
    scans = conn.get_avail_scans(y, mo, d, st)
    if not scans:
        raise SourceUnavailable("no NEXRAD scans for %s on %s" % (st, cfg.case.date))
    pick = min(scans, key=lambda s: abs((s.scan_time.replace(tzinfo=None) - start).total_seconds())
               if s.scan_time else 1e18)
    res = conn.download(pick, os.path.dirname(path))
    got = res.success[0].filepath if getattr(res, "success", None) else pick.filename
    if os.path.abspath(got) != os.path.abspath(path):
        import shutil; shutil.copy(got, path)
    cache.record("nexrad", key, path, {"scan": pick.filename})
    return path


def load(cfg, cache, sweep=0):
    """Read a Level II sweep -> radar volume dict (same schema as ``synthetic_radar``)."""
    pyart = require("pyart", "NEXRAD reader", "(pip install arm_pyart) or xradar")
    path = download(cfg, cache)
    radar = pyart.io.read_nexrad_archive(path)
    s = radar.get_slice(sweep)
    az = np.asarray(radar.azimuth["data"][s]); rng = np.asarray(radar.range["data"])
    el = float(np.mean(radar.elevation["data"][s]))
    lat = np.asarray(radar.gate_latitude["data"][s]); lon = np.asarray(radar.gate_longitude["data"][s])
    alt = np.asarray(radar.gate_altitude["data"][s])
    x = np.asarray(radar.gate_x["data"][s]); y = np.asarray(radar.gate_y["data"][s])
    def fld(*names):
        for n in names:
            if n in radar.fields:
                return np.ma.filled(radar.fields[n]["data"][s], np.nan)
        return None
    return {"azimuth_deg": az, "range_m": rng, "elevation_deg": el,
            "lat": lat, "lon": lon, "alt_m": alt, "x_m": x, "y_m": y,
            "reflectivity": fld("reflectivity", "REF"),
            "radial_velocity": fld("velocity", "VEL"),
            "spectrum_width": fld("spectrum_width", "SW"),
            "differential_reflectivity": fld("differential_reflectivity", "ZDR"),
            "cross_correlation_ratio": fld("cross_correlation_ratio", "RHO"),
            "differential_phase": fld("differential_phase", "PHI"),
            "radar_lat": float(radar.latitude["data"][0]),
            "radar_lon": float(radar.longitude["data"][0]),
            "radar_alt_m": float(radar.altitude["data"][0]),
            "station": cfg.data.radar_station, "source": "nexrad"}


def available():
    from .base import available as _av
    return _av(["pyart"])
