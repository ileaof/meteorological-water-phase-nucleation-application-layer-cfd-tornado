"""Radiosonde ingestion (ROADMAP §3a, source D): CSV / whitespace text profiles.

Accepted columns (any subset; SI or common units auto-detected by name):
``height_m pressure_Pa|pressure_hPa temperature_K|temperature_C dewpoint_K|dewpoint_C
relative_humidity specific_humidity mixing_ratio u_ms v_ms wind_dir_deg wind_speed_ms``.
Returns an SI profile dict; :func:`profile_to_basestate` builds the model base state (reusing
``storm_dynamics.soundings.from_observed_sounding``), and :func:`profile_diagnostics` returns
CAPE/CIN/LCL/LFC/EL/shear/SRH.
"""
from __future__ import annotations

import numpy as np

from .. import thermo, units


def read_sounding(path):
    """Read a sounding CSV/text into an SI profile dict (ascending height)."""
    import pandas as pd
    sep = None if str(path).lower().endswith(".csv") else r"\s+"
    df = pd.read_csv(path, sep=sep, engine="python", comment="#")
    df.columns = [c.strip().lower() for c in df.columns]
    g = lambda *names: next((df[n].to_numpy(float) for n in names if n in df.columns), None)

    z = g("height_m", "height", "z", "geopotential_height_m")
    p = g("pressure_pa")
    if p is None and g("pressure_hpa") is not None:
        p = units.hpa_to_pa(g("pressure_hpa"))
    T = g("temperature_k")
    if T is None and g("temperature_c") is not None:
        T = units.celsius_to_kelvin(g("temperature_c"))
    Td = g("dewpoint_k")
    if Td is None and g("dewpoint_c") is not None:
        Td = units.celsius_to_kelvin(g("dewpoint_c"))
    rh = g("relative_humidity", "rh")
    if rh is not None and np.nanmax(rh) > 1.5:
        rh = rh / 100.0
    q = g("specific_humidity", "qv")
    r = g("mixing_ratio")
    u = g("u_ms", "u"); v = g("v_ms", "v")
    wdir = g("wind_dir_deg", "wind_direction", "drct")
    wspd = g("wind_speed_ms", "wind_speed", "sknt")
    if u is None and wdir is not None and wspd is not None:
        if g("sknt") is not None and wspd is g("sknt"):
            wspd = units.knots_to_ms(wspd)
        u, v = units.wind_dir_speed_to_uv(wdir, wspd)
    if q is None and r is not None:
        q = thermo.specific_humidity_from_mixing_ratio(r)
    if q is None and Td is not None and p is not None:
        q = thermo.specific_humidity_from_dewpoint(Td, p)
    if q is None and rh is not None and T is not None and p is not None:
        q = thermo.specific_humidity_from_rh(rh, T, p)

    prof = {"height_m": z, "pressure_Pa": p, "temperature_K": T, "dewpoint_K": Td,
            "relative_humidity": rh, "specific_humidity": q, "u_ms": u, "v_ms": v}
    prof = {k: val for k, val in prof.items() if val is not None}
    if "height_m" in prof:                                   # sort by height, ascending
        order = np.argsort(prof["height_m"])
        prof = {k: np.asarray(val)[order] for k, val in prof.items()}
    return prof


def profile_to_basestate(grid, prof):
    """Build a model ``BaseState`` from an SI profile dict (delegates to the verified
    ``storm_dynamics.soundings.from_observed_sounding``)."""
    from storm_dynamics.soundings import from_observed_sounding
    kw = dict(pressure_hPa=units.pa_to_hpa(prof["pressure_Pa"]),
              height_m=prof["height_m"],
              temperature_C=units.kelvin_to_celsius(prof["temperature_K"]))
    if "specific_humidity" in prof:
        kw["qv_kgkg"] = prof["specific_humidity"]
    elif "dewpoint_K" in prof:
        kw["dewpoint_C"] = units.kelvin_to_celsius(prof["dewpoint_K"])
    if "u_ms" in prof and "v_ms" in prof:
        kw["u_ms"] = prof["u_ms"]; kw["v_ms"] = prof["v_ms"]
    else:
        kw["wind_dir_deg"] = prof.get("wind_dir_deg", np.zeros_like(prof["height_m"]))
        kw["wind_speed_ms"] = prof.get("wind_speed_ms", np.zeros_like(prof["height_m"]))
    return from_observed_sounding(grid, **kw)


def profile_diagnostics(grid, prof):
    """CAPE/CIN/LCL/LFC/EL/shear (+SRH) for the profile (reuses the engine diagnostics)."""
    from meteorological_flow.base_state import sounding_diagnostics
    from storm_dynamics import soundings as snd
    base = profile_to_basestate(grid, prof)
    d = sounding_diagnostics(base)
    d["SRH_0_3km_m2_s2"] = snd.storm_relative_helicity(base)
    return {k: (float(v) if np.isscalar(v) or (hasattr(v, "ndim") and v.ndim == 0) else v)
            for k, v in d.items() if k in ("CAPE_J_kg", "CIN_J_kg", "LCL_m", "LFC_m", "EL_m",
                                           "shear_0_6km_m_s", "SRH_0_3km_m2_s2", "freezing_level_m")}
