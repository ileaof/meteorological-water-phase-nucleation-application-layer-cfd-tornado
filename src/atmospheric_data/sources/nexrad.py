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

import numpy as np

from .base import require, SourceUnavailable


def download(cfg, cache, when=None):
    """Fetch the nearest Level II volume for the station/time from AWS (no creds)."""
    import requests
    st = cfg.data.radar_station
    date = cfg.case.date.replace("-", "")
    hh = cfg.case.start_time_utc.split(":")[0].zfill(2)
    key = "%s_%s_%s" % (st, date, hh)
    path = cache.path("nexrad", key, ".ar2v")
    cache.require_offline_ok("nexrad", key, ".ar2v")
    if cache.has("nexrad", key, ".ar2v"):
        return path
    prefix = "%s/%s/%s/%s/%s" % (date[:4], date[4:6], date[6:8], st, st)  # list bucket
    idx = ("https://noaa-nexrad-level2.s3.amazonaws.com/?list-type=2&prefix=%s"
           % ("/".join([date[:4], date[4:6], date[6:8], st, "%s%s_%s" % (st, date, hh)])))
    r = requests.get(idx, timeout=60); r.raise_for_status()
    import re
    keys = re.findall(r"<Key>([^<]+)</Key>", r.text)
    if not keys:
        raise SourceUnavailable("no NEXRAD volume found for %s at %sZ" % (st, hh))
    url = "https://noaa-nexrad-level2.s3.amazonaws.com/" + keys[0]
    rr = requests.get(url, stream=True, timeout=180); rr.raise_for_status()
    with open(path, "wb") as f:
        for chunk in rr.iter_content(1 << 20):
            f.write(chunk)
    cache.record("nexrad", key, path, {"url": url})
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
