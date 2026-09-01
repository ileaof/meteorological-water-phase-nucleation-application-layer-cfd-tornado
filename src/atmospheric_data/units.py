"""Unit conversions to the model's SI core (ROADMAP §3a).

Real datasets arrive in mixed units (hPa, degC, knots, geopotential, wind dir/speed).  Every
conversion here is explicit and logged by the caller so the metadata records exactly what was
applied (task requirement 7: "registre ... unidades e transformações aplicadas").
"""
from __future__ import annotations

import numpy as np

KT_TO_MS = 0.514444
DEG2RAD = np.pi / 180.0


def hpa_to_pa(x):
    return np.asarray(x, float) * 100.0


def pa_to_hpa(x):
    return np.asarray(x, float) / 100.0


def celsius_to_kelvin(x):
    return np.asarray(x, float) + 273.15


def kelvin_to_celsius(x):
    return np.asarray(x, float) - 273.15


def knots_to_ms(x):
    return np.asarray(x, float) * KT_TO_MS


def wind_dir_speed_to_uv(direction_deg, speed_ms):
    """Meteorological wind (direction FROM, degrees; speed m/s) -> (u_east, v_north) [m/s].

    ``u = -spd sin(dir)``, ``v = -spd cos(dir)`` (the standard convention: a 270 deg wind at
    10 m/s is a westerly, u=+10, v=0)."""
    d = np.asarray(direction_deg, float) * DEG2RAD
    s = np.asarray(speed_ms, float)
    return -s * np.sin(d), -s * np.cos(d)


def uv_to_dir_speed(u, v):
    """(u_east, v_north) [m/s] -> (direction_from_deg, speed_ms)."""
    u = np.asarray(u, float); v = np.asarray(v, float)
    spd = np.hypot(u, v)
    d = (np.degrees(np.arctan2(-u, -v))) % 360.0
    return d, spd


# canonical SI units for the internal standard variables (used to tag metadata)
SI_UNITS = {
    "T": "K", "p": "Pa", "rho": "kg m-3", "theta": "K", "theta_v": "K",
    "qv": "kg kg-1", "qc": "kg kg-1", "qr": "kg kg-1", "qi": "kg kg-1",
    "qs": "kg kg-1", "qg": "kg kg-1",
    "u": "m s-1", "v": "m s-1", "w": "m s-1",
    "terrain": "m", "reflectivity": "dBZ",
    "z": "m", "y": "m", "x": "m",
}
