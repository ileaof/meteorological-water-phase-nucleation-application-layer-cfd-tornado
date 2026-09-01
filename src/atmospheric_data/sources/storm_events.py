"""NOAA Storm Events / SWDI (ROADMAP §3a, source F) -- case selection + track VALIDATION.

Reads CSV/JSON of tornado events (time, start/end lat-lon, EF rating, track length/width).
Used ONLY to pick a case and to validate the simulated track's location/timing -- NEVER to
build the flow fields (task: "nunca para construir artificialmente os campos de escoamento").
Shapefile support is optional (needs geopandas); CSV/JSON always work.
"""
from __future__ import annotations

import json

import numpy as np


def read_storm_events(path):
    """Read tornado events -> list of dicts with SI-ish fields.  Supports CSV, JSON, and
    (optional) shapefile.  Each event: ``begin_time end_time begin_lat begin_lon end_lat
    end_lon ef_rating track_length_km track_width_m``."""
    p = str(path).lower()
    if p.endswith(".json"):
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        rows = data if isinstance(data, list) else data.get("events", [])
        return [_norm(r) for r in rows]
    if p.endswith((".shp",)):
        from .base import require
        gpd = require("geopandas", "Storm Events shapefile", "or use CSV/JSON")
        gdf = gpd.read_file(path)
        return [_norm(r) for r in gdf.to_dict("records")]
    import pandas as pd
    df = pd.read_csv(path, comment="#")
    df.columns = [c.strip().lower() for c in df.columns]
    return [_norm(r) for r in df.to_dict("records")]


def _norm(r):
    r = {str(k).strip().lower(): v for k, v in r.items()}
    g = lambda *ks, d=None: next((r[k] for k in ks if k in r and r[k] == r[k]), d)  # skip NaN
    return {"begin_time": g("begin_time", "begin_date_time", "start_time_utc"),
            "end_time": g("end_time", "end_date_time", "end_time_utc"),
            "begin_lat": _f(g("begin_lat", "begin_latitude", "start_lat")),
            "begin_lon": _f(g("begin_lon", "begin_longitude", "start_lon")),
            "end_lat": _f(g("end_lat", "end_latitude")),
            "end_lon": _f(g("end_lon", "end_longitude")),
            "ef_rating": g("ef_rating", "tor_f_scale", "magnitude"),
            "track_length_km": _f(g("track_length_km", "length_km", "tor_length")),
            "track_width_m": _f(g("track_width_m", "width_m", "tor_width"))}


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None
