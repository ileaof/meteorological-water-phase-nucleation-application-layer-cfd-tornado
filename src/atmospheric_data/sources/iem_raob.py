"""Real radiosonde download (ROADMAP §3a, source D) — Iowa Environmental Mesonet (IEM) RAOB.

Downloads a **real observed sounding** (no credentials, no heavy dependencies) from the IEM
archive JSON and returns an SI profile dict for :func:`sources.sounding.profile_to_basestate`.
This is genuine observational data — e.g. KOUN (Norman, OK) near the Moore 2013 tornado.

    https://mesonet.agron.iastate.edu/json/raob.py?ts=<ISO8601Z>&station=<ID>

Missing levels/winds are handled explicitly (masked; winds re-interpolated onto the
thermodynamic levels) — never invented.
"""
from __future__ import annotations

import numpy as np

from .. import units


def download_sounding(station="KOUN", when="2013-05-21T00:00:00Z", cache=None, offline=False,
                      timeout=45):
    """Fetch a real RAOB and return an SI profile dict (height AGL, pressure, T, Td, u, v).

    ``cache`` (optional) stores the raw JSON; ``offline`` uses only the cached copy."""
    import json
    key = "%s_%s" % (station, when.replace(":", "").replace("-", ""))
    raw = None
    if cache is not None:
        path = cache.path("iem_raob", key, ".json")
        cache.require_offline_ok("iem_raob", key, ".json")
        if cache.has("iem_raob", key, ".json"):
            with open(path, "r", encoding="utf-8") as f:
                raw = json.load(f)
    if raw is None:
        if offline:
            raise FileNotFoundError("offline: RAOB %s %s not cached" % (station, when))
        import requests
        url = "https://mesonet.agron.iastate.edu/json/raob.py?ts=%s&station=%s" % (when, station)
        r = requests.get(url, timeout=timeout, headers={"User-Agent": "met_h2o research (educational)"})
        r.raise_for_status()
        raw = r.json()
        if cache is not None:
            with open(cache.path("iem_raob", key, ".json"), "w", encoding="utf-8") as f:
                json.dump(raw, f)
            cache.record("iem_raob", key, cache.path("iem_raob", key, ".json"), {"url": url})
    prof = raw["profiles"][0] if "profiles" in raw else raw
    return _to_si(prof)


def _to_si(prof):
    rows = prof["profile"]
    col = lambda k: np.array([row.get(k, None) for row in rows], dtype=object)
    to_f = lambda a: np.array([np.nan if v is None else float(v) for v in a], float)
    p = to_f(col("pres")); z = to_f(col("hght")); T = to_f(col("tmpc")); Td = to_f(col("dwpc"))
    drct = to_f(col("drct")); sknt = to_f(col("sknt"))
    good = np.isfinite(p) & np.isfinite(z) & np.isfinite(T) & np.isfinite(Td)
    p, z, T, Td = p[good], z[good], T[good], Td[good]
    dr, sk = drct[good], sknt[good]
    # winds: interpolate over the levels that HAVE a report (fill gaps on the thermo heights)
    wgood = np.isfinite(dr) & np.isfinite(sk)
    if wgood.sum() >= 2:
        u_w, v_w = units.wind_dir_speed_to_uv(dr[wgood], units.knots_to_ms(sk[wgood]))
        u = np.interp(z, z[wgood], u_w); v = np.interp(z, z[wgood], v_w)
    else:
        u = np.zeros_like(z); v = np.zeros_like(z)
    order = np.argsort(z)
    return {"height_m": (z - z[0])[order], "pressure_Pa": units.hpa_to_pa(p)[order],
            "temperature_K": units.celsius_to_kelvin(T)[order],
            "dewpoint_K": units.celsius_to_kelvin(Td)[order],
            "u_ms": u[order], "v_ms": v[order],
            "station": prof.get("station", ""), "valid": prof.get("valid", "")}
