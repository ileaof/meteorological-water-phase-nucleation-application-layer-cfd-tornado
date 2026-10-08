"""Real-time global model (ROADMAP §3a, source G) -- NOAA GFS via THREDDS NCSS.

The gap this closes: every other gridded source here is **US-only** (HRRR, NEXRAD) or runs
days behind (ERA5), and the IEM radiosonde archive carries no stations outside the US.  So
for a non-US case there was no wind *profile* and no real-time IC/BC at all -- surface
METARs give moisture and nothing above it.

Why NCSS and not GRIB or OPeNDAP:

* **GRIB2** (NOMADS filter, AWS ``noaa-gfs-bdp-pds``) needs ``cfgrib``/``eccodes``, which are
  optional and absent on plain Windows installs.
* **NOMADS OPeNDAP is retired** (NWS SCN25-81) -- ``/dods/`` now serves an HTML notice.
* **THREDDS NCSS** subsets server-side and returns CSV or NetCDF over plain HTTPS, so it
  needs only ``requests`` (profiles) or ``netCDF4`` (grids), both already required.

Profiles come back in the dict shape :func:`sources.sounding.profile_to_basestate` and
:func:`sources.sounding.profile_diagnostics` already consume, so CAPE/CIN/LCL and the
low-level geometry diagnostics work unchanged.

**Southern-hemisphere caution.** GFS reports negative ``Storm_relative_helicity`` south of
the equator because cyclonic rotation there is *clockwise* and the favoured supercell is the
**left**-mover.  ``storm_dynamics.soundings.bunkers_storm_motion`` deviates to the right by
default (northern-hemisphere convention); pass ``hemisphere="south"`` -- or let
:func:`storm_dynamics.soundings.bunkers_storm_motion` infer it from ``latitude_deg`` -- or
every storm-relative quantity is computed for the wrong storm.
"""
from __future__ import annotations

import csv
import io

import numpy as np

from .. import thermo, units

#: THREDDS NetCDF Subset Service for the GFS 0.25-degree "Best" (analysis + forecast) series.
GFS_NCSS = "https://thredds.ucar.edu/thredds/ncss/grid/grib/NCEP/GFS/Global_0p25deg/Best"

#: Isobaric variables that make up a profile.
PROFILE_VARS = (
    "Temperature_isobaric",
    "Relative_humidity_isobaric",
    "u-component_of_wind_isobaric",
    "v-component_of_wind_isobaric",
    "Geopotential_height_isobaric",
)

#: Single-level fields fetched alongside: terrain and surface pressure are needed to put the
#: profile on height-above-ground and to drop levels that lie *below* ground (at Brasilia,
#: 1051 m, the 1000 hPa level is underground and GFS reports an extrapolated value there).
#: The rest are GFS's own convective diagnostics, kept purely as an independent cross-check
#: against the ones this project computes.
SURFACE_VARS = (
    "Pressure_surface",
    "Geopotential_height_surface",
    "Temperature_height_above_ground",
    "Dewpoint_temperature_height_above_ground",
    "u-component_of_wind_height_above_ground",
    "v-component_of_wind_height_above_ground",
    "Convective_available_potential_energy_surface",
    "Convective_inhibition_surface",
    "Storm_relative_helicity_height_above_ground_layer",
)

#: Height [m AGL] assigned to the screen-level (2 m) surface observation when it is
#: prepended to the profile.
SURFACE_LEVEL_HEIGHT_M = 2.0


def _get(url, timeout, cache=None, source="gfs", key=None):
    """HTTPS GET returning text, with optional caching."""
    if cache is not None and key is not None and cache.has(source, key, ".csv"):
        with open(cache.path(source, key, ".csv"), "r", encoding="utf-8") as f:
            return f.read()
    import requests

    from . import _tls
    r = requests.get(url, timeout=timeout,
                     headers={"User-Agent": "met_h2o research (educational)"},
                     verify=_tls.ca_bundle())
    r.raise_for_status()
    text = r.text
    if cache is not None and key is not None:
        with open(cache.path(source, key, ".csv"), "w", encoding="utf-8") as f:
            f.write(text)
        cache.record(source, key, cache.path(source, key, ".csv"), {"url": url})
    return text


def _parse_ncss_csv(text):
    """NCSS CSV -> (list of row dicts, column order).

    Headers carry their units inline (``Temperature_isobaric[unit="K"]``); the unit is
    stripped from the key and the bare name kept.
    """
    rows = list(csv.reader(io.StringIO(text.strip())))
    if len(rows) < 2:
        raise ValueError("NCSS returned no data rows (got %d line(s))" % len(rows))
    header = [h.split("[")[0].strip() for h in rows[0]]
    out = []
    for raw in rows[1:]:
        if not raw or len(raw) != len(header):
            continue
        out.append(dict(zip(header, raw)))
    if not out:
        raise ValueError("NCSS returned a header but no usable rows")
    return out, header


def _to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def download_surface(lat, lon, when="present", cache=None, offline=False, timeout=90,
                     variables=SURFACE_VARS):
    """Single-level GFS fields at a point -> ``{variable: value}`` (SI as GFS reports them).

    Each variable is requested separately: NCSS rejects some combinations of single-level
    variables in one request, and one bad name would otherwise lose the whole batch.
    """
    out = {}
    for var in variables:
        key = "sfc_%s_%.3f_%.3f_%s" % (var, lat, lon, str(when).replace(":", ""))
        if offline and not (cache is not None and cache.has("gfs", key, ".csv")):
            raise FileNotFoundError("offline: GFS surface %s not cached" % var)
        url = "%s?var=%s&latitude=%s&longitude=%s&time=%s&accept=csv" % (
            GFS_NCSS, var, lat, lon, when)
        try:
            rows, header = _parse_ncss_csv(_get(url, timeout, cache, key=key))
        except Exception:
            out[var] = float("nan")       # a missing diagnostic must not kill the profile
            continue
        out[var] = _to_float(rows[0][header[-1]])
        out.setdefault("valid_utc", rows[0].get("time", ""))
    return out


def _prepend_surface_level(z, p, T, rh, u, v, sfc):
    """Put the real screen-level observation at the bottom of the profile.

    Dropping below-ground isobaric levels is correct but costs the *surface parcel*: over
    terrain the lowest surviving level can sit 100-200 m up, and a parcel lifted from there
    systematically under-reads CAPE compared with a true surface-based parcel.  GFS carries
    the 2 m temperature/dewpoint, 10 m wind and surface pressure, so the bottom of the
    column is an observation-equivalent, not an extrapolation.

    Returns the arrays unchanged when the surface fields are missing.
    """
    p_sfc = sfc.get("Pressure_surface", float("nan"))
    T_sfc = sfc.get("Temperature_height_above_ground", float("nan"))
    Td_sfc = sfc.get("Dewpoint_temperature_height_above_ground", float("nan"))
    if not (np.isfinite(p_sfc) and np.isfinite(T_sfc) and np.isfinite(Td_sfc)):
        return z, p, T, rh, u, v
    if len(p) and p_sfc <= p.max():
        return z, p, T, rh, u, v          # nothing below the lowest kept level
    # RH from the dewpoint via the project's own saturation formula (no second convention)
    e = thermo.saturation_vapor_pressure(np.array([Td_sfc], float))
    es = thermo.saturation_vapor_pressure(np.array([T_sfc], float))
    rh_sfc = float(np.clip(e[0] / es[0], 0.0, 1.0)) if es[0] > 0 else float("nan")
    if not np.isfinite(rh_sfc):
        return z, p, T, rh, u, v
    u_sfc = sfc.get("u-component_of_wind_height_above_ground", float("nan"))
    v_sfc = sfc.get("v-component_of_wind_height_above_ground", float("nan"))
    if not (np.isfinite(u_sfc) and np.isfinite(v_sfc)):
        u_sfc, v_sfc = u[0], v[0]         # reuse the lowest level rather than invent calm
    return (np.concatenate([[SURFACE_LEVEL_HEIGHT_M], z]),
            np.concatenate([[p_sfc], p]),
            np.concatenate([[T_sfc], T]),
            np.concatenate([[rh_sfc], rh]),
            np.concatenate([[u_sfc], u]),
            np.concatenate([[v_sfc], v]))


def download_profile(lat, lon, when="present", cache=None, offline=False, timeout=90,
                     drop_below_ground=True, include_surface=True):
    """GFS vertical profile at a point -> SI profile dict for ``sounding.profile_to_basestate``.

    Returns ``height_m`` (**above ground**), ``pressure_Pa``, ``temperature_K``,
    ``relative_humidity`` (0-1), ``specific_humidity``, ``u_ms``, ``v_ms``, ascending in
    height, plus ``station``/``valid``/``latitude_deg``/``longitude_deg``/``terrain_m``/
    ``surface_pressure_Pa`` and a ``gfs_diagnostics`` sub-dict (GFS's own CAPE/CIN/SRH,
    for cross-checking -- never as a substitute for computing them from the column).

    With ``drop_below_ground`` (default) isobaric levels whose pressure exceeds the surface
    pressure are discarded rather than trusted: GFS fills them by extrapolation, so over
    terrain they are not atmosphere.  With ``include_surface`` (default) the real screen-level
    observation is then prepended, so the column still starts at the ground -- without it a
    parcel would be lifted from the lowest surviving isobaric level (134 m AGL at Campo
    Grande) and CAPE would read systematically low.

    **Vertical-resolution warning.** CAPE from this profile is sensitive to the model grid it
    is interpolated onto: measured at Dourados, ``dz = 200 m`` gave 100 J/kg against 729 J/kg
    at ``dz <= 100 m`` -- a 7x under-read. Use ``dz <= 100 m`` for CAPE, and treat a value
    from a coarse grid as meaningless rather than conservative.
    """
    key = "prof_%.3f_%.3f_%s" % (lat, lon, str(when).replace(":", ""))
    if offline and not (cache is not None and cache.has("gfs", key, ".csv")):
        raise FileNotFoundError("offline: GFS profile %.3f,%.3f not cached" % (lat, lon))
    url = "%s?%s&latitude=%s&longitude=%s&time=%s&accept=csv" % (
        GFS_NCSS, "&".join("var=%s" % v for v in PROFILE_VARS), lat, lon, when)
    rows, _header = _parse_ncss_csv(_get(url, timeout, cache, key=key))

    p = np.array([_to_float(r.get("alt")) for r in rows], float)
    T = np.array([_to_float(r.get("Temperature_isobaric")) for r in rows], float)
    rh = np.array([_to_float(r.get("Relative_humidity_isobaric")) for r in rows], float) / 100.0
    u = np.array([_to_float(r.get("u-component_of_wind_isobaric")) for r in rows], float)
    v = np.array([_to_float(r.get("v-component_of_wind_isobaric")) for r in rows], float)
    z_msl = np.array([_to_float(r.get("Geopotential_height_isobaric")) for r in rows], float)

    sfc = download_surface(lat, lon, when, cache=cache, offline=offline, timeout=timeout)
    terrain = sfc.get("Geopotential_height_surface", float("nan"))
    p_sfc = sfc.get("Pressure_surface", float("nan"))

    good = np.isfinite(p) & np.isfinite(T) & np.isfinite(z_msl) & np.isfinite(u) & np.isfinite(v)
    if drop_below_ground and np.isfinite(p_sfc):
        good &= p <= p_sfc
    if good.sum() < 3:
        raise ValueError("GFS profile at %.3f,%.3f has only %d usable level(s) "
                         "(surface pressure %.0f Pa)" % (lat, lon, int(good.sum()), p_sfc))
    p, T, rh, u, v, z_msl = p[good], T[good], rh[good], u[good], v[good], z_msl[good]

    z_agl = z_msl - (terrain if np.isfinite(terrain) else z_msl.min())
    rh = np.clip(rh, 0.0, 1.0)
    if include_surface:
        order0 = np.argsort(z_agl)
        z_agl, p, T, rh, u, v = _prepend_surface_level(
            z_agl[order0], p[order0], T[order0], rh[order0], u[order0], v[order0], sfc)
        rh = np.clip(rh, 0.0, 1.0)
    qv = thermo.specific_humidity_from_rh(rh, T, p)

    order = np.argsort(z_agl)
    return {
        "height_m": z_agl[order], "pressure_Pa": p[order], "temperature_K": T[order],
        "relative_humidity": rh[order], "specific_humidity": qv[order],
        "u_ms": u[order], "v_ms": v[order],
        "station": "GFS_%.3f_%.3f" % (lat, lon), "valid": sfc.get("valid_utc", ""),
        "latitude_deg": float(lat), "longitude_deg": float(lon),
        "terrain_m": float(terrain), "surface_pressure_Pa": float(p_sfc),
        "gfs_diagnostics": {
            "CAPE_J_kg": sfc.get("Convective_available_potential_energy_surface", float("nan")),
            "CIN_J_kg": sfc.get("Convective_inhibition_surface", float("nan")),
            "SRH_m2_s2": sfc.get("Storm_relative_helicity_height_above_ground_layer", float("nan")),
        },
    }


def download_grid(out_path, north, south, west, east, when="present", variables=PROFILE_VARS,
                  timeout=300, cache=None, offline=False, accept="netcdf4"):
    """Download a lat/lon box of GFS as NetCDF (for IC/BC) and return the written path.

    Server-side subsetting keeps this small; reading the result needs only ``netCDF4``.
    The file is left in the model's own hands -- ``atmospheric_data.ic_bc`` consumes gridded
    fields, and this function does not pretend to do that coupling itself.
    """
    import os
    key = "grid_%.2f_%.2f_%.2f_%.2f_%s" % (north, south, west, east, str(when).replace(":", ""))
    if offline:
        if os.path.exists(out_path):
            return out_path
        raise FileNotFoundError("offline: GFS grid %s not present at %s" % (key, out_path))
    url = ("%s?%s&north=%s&south=%s&west=%s&east=%s&horizStride=1&time=%s&accept=%s"
           % (GFS_NCSS, "&".join("var=%s" % v for v in variables),
              north, south, west, east, when, accept))
    import requests

    from . import _tls
    r = requests.get(url, timeout=timeout, stream=True,
                     headers={"User-Agent": "met_h2o research (educational)"},
                     verify=_tls.ca_bundle())
    r.raise_for_status()
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    with open(out_path, "wb") as f:
        for chunk in r.iter_content(chunk_size=1 << 20):
            if chunk:
                f.write(chunk)
    if cache is not None:
        cache.record("gfs", key, out_path, {"url": url})
    return out_path


__all__ = ["GFS_NCSS", "PROFILE_VARS", "SURFACE_VARS",
           "download_profile", "download_surface", "download_grid"]
