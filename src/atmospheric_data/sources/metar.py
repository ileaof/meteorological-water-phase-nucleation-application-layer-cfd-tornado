"""METAR/ASOS surface observations (ROADMAP §3a, source E) -- for SURFACE validation only.

Reads a CSV of station reports; converts to SI; used to check the model's near-surface fields
(T, Td, p, wind), never to build the flow.  Columns auto-detected: ``station lat lon
elevation_m valid|time temperature_C|_K dewpoint_C|_K pressure_hPa|_Pa wind_dir_deg
wind_speed_ms|_kt wind_gust visibility_m precip_mm``.

:func:`download_metar` adds a live fetch (no credentials) from the NOAA Aviation Weather
API, which -- unlike the IEM RAOB archive used by :mod:`.iem_raob` -- has **global**
coverage, including Brazilian stations (SBCG, SBBR, SBCY, SBGO, ...).  Surface reports
cannot give a wind *profile*, so shear and SRH stay out of reach; what they do give is
boundary-layer moisture, hence the LCL (:func:`lcl_estimate_m`), which is one of the
discriminators the tornadogenesis study isolated.
"""
from __future__ import annotations

import numpy as np

from .. import units

#: Live METAR endpoint (no credentials, global coverage).
METAR_API = "https://aviationweather.gov/api/data/metar"


def read_metar_csv(path):
    """Read METAR/ASOS CSV -> pandas DataFrame in SI (adds ``*_si`` columns where converted)."""
    import pandas as pd
    df = pd.read_csv(path, comment="#")
    df.columns = [c.strip().lower() for c in df.columns]
    if "temperature_c" in df:
        df["temperature_k"] = units.celsius_to_kelvin(df["temperature_c"].to_numpy(float))
    if "dewpoint_c" in df:
        df["dewpoint_k"] = units.celsius_to_kelvin(df["dewpoint_c"].to_numpy(float))
    if "pressure_hpa" in df:
        df["pressure_pa"] = units.hpa_to_pa(df["pressure_hpa"].to_numpy(float))
    if "wind_speed_kt" in df:
        df["wind_speed_ms"] = units.knots_to_ms(df["wind_speed_kt"].to_numpy(float))
    if "wind_dir_deg" in df and "wind_speed_ms" in df:
        u, v = units.wind_dir_speed_to_uv(df["wind_dir_deg"].to_numpy(float),
                                          df["wind_speed_ms"].to_numpy(float))
        df["u_ms"] = u; df["v_ms"] = v
    return df


# ---------------------------------------------------------------------------
# live download (global coverage, no credentials)
# ---------------------------------------------------------------------------
def _num(value):
    """Coerce an API field to float, or ``None``.

    The API mixes numbers with sentinel strings -- ``wdir`` is ``"VRB"`` for a variable
    wind, ``visib`` is ``"6+"`` -- so a plain ``float()`` would raise on real reports.
    """
    if value is None or isinstance(value, bool):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def download_metar(stations=("KOUN",), cache=None, offline=False, timeout=45, hours=None):
    """Fetch current METARs and return a list of SI report dicts (newest first per API).

    ``stations`` is an iterable of ICAO ids (``"SBCG"``, ``"KOUN"``, ...) or a single
    comma-separated string.  ``hours`` requests that many hours of history instead of
    only the latest cycle.  ``cache`` stores the raw JSON; ``offline`` uses only the
    cached copy.

    Each returned dict carries ``station, valid_utc, lat_deg, lon_deg, elevation_m,
    temperature_K, dewpoint_K, pressure_Pa, wind_dir_deg, wind_speed_ms, u_ms, v_ms,
    lcl_estimate_m, raw`` -- with ``None`` wherever the report omits the field or
    reports it non-numerically (a variable wind leaves ``wind_dir_deg``/``u_ms``/
    ``v_ms`` as ``None`` rather than inventing a direction).
    """
    import json
    ids = stations if isinstance(stations, str) else ",".join(stations)
    key = "%s_h%s" % (ids.replace(",", "-"), hours if hours is not None else "latest")
    raw = None
    if cache is not None:
        path = cache.path("metar", key, ".json")
        cache.require_offline_ok("metar", key, ".json")
        if cache.has("metar", key, ".json"):
            with open(path, "r", encoding="utf-8") as f:
                raw = json.load(f)
    if raw is None:
        if offline:
            raise FileNotFoundError("offline: METAR %s not cached" % ids)
        import requests

        from . import _tls
        url = "%s?ids=%s&format=json" % (METAR_API, ids)
        if hours is not None:
            url += "&hours=%d" % int(hours)
        r = requests.get(url, timeout=timeout,
                         headers={"User-Agent": "met_h2o research (educational)"},
                         verify=_tls.ca_bundle())
        r.raise_for_status()
        raw = r.json()
        if cache is not None:
            with open(cache.path("metar", key, ".json"), "w", encoding="utf-8") as f:
                json.dump(raw, f)
            cache.record("metar", key, cache.path("metar", key, ".json"), {"url": url})
    if not isinstance(raw, list):
        raise ValueError("unexpected METAR payload: %r" % type(raw))
    return [_report_to_si(rec) for rec in raw]


def _report_to_si(rec):
    T = _num(rec.get("temp"))
    Td = _num(rec.get("dewp"))
    drct = _num(rec.get("wdir"))          # None when the report says "VRB"
    sknt = _num(rec.get("wspd"))
    speed = units.knots_to_ms(sknt) if sknt is not None else None
    if drct is not None and speed is not None:
        u, v = units.wind_dir_speed_to_uv(np.array([drct]), np.array([speed]))
        u, v = float(u[0]), float(v[0])
    else:
        u = v = None
    altim = _num(rec.get("altim"))
    T_K = units.celsius_to_kelvin(T) if T is not None else None
    Td_K = units.celsius_to_kelvin(Td) if Td is not None else None
    return {
        "station": rec.get("icaoId", ""),
        "name": rec.get("name", ""),
        "valid_utc": rec.get("reportTime", ""),
        "lat_deg": _num(rec.get("lat")), "lon_deg": _num(rec.get("lon")),
        "elevation_m": _num(rec.get("elev")),
        "temperature_K": T_K, "dewpoint_K": Td_K,
        "pressure_Pa": units.hpa_to_pa(altim) if altim is not None else None,
        "wind_dir_deg": drct, "wind_speed_ms": speed, "u_ms": u, "v_ms": v,
        "gust_ms": units.knots_to_ms(_num(rec.get("wgst"))) if _num(rec.get("wgst")) is not None else None,
        "lcl_estimate_m": lcl_estimate_m(T_K, Td_K),
        "raw": rec.get("rawOb", ""),
    }


def lcl_estimate_m(temperature_K, dewpoint_K):
    """LCL height [m] from the surface dewpoint depression, or ``None``.

    Lawrence (2005) approximation ``z_LCL ~ 125 * (T - Td)`` with T, Td in K or C.
    Good to roughly +/-100 m for a well-mixed boundary layer, which is all a *surface*
    report can support -- a real LCL needs the column, i.e.
    :func:`sources.sounding.profile_diagnostics`.  Use it as a moisture discriminator,
    not as a precise cloud base.
    """
    if temperature_K is None or dewpoint_K is None:
        return None
    return float(125.0 * max(0.0, float(temperature_K) - float(dewpoint_K)))
