"""METAR/ASOS surface observations (ROADMAP §3a, source E) -- for SURFACE validation only.

Reads a CSV of station reports; converts to SI; used to check the model's near-surface fields
(T, Td, p, wind), never to build the flow.  Columns auto-detected: ``station lat lon
elevation_m valid|time temperature_C|_K dewpoint_C|_K pressure_hPa|_Pa wind_dir_deg
wind_speed_ms|_kt wind_gust visibility_m precip_mm``.
"""
from __future__ import annotations

import numpy as np

from .. import units


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
