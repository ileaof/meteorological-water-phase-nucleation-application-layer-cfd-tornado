"""NOAA HRRR ingestion (ROADMAP §3a, source A): GRIB2 -> unified internal format.

Reads HRRR pressure-level GRIB2 (via ``cfgrib``/``xarray``) and maps its variables to the
model's standard names in SI.  Download is from the **AWS Open Data** bucket (no credentials):
``s3://noaa-hrrr-bdp-pds``.  Requires ``cfgrib`` + ``eccodes`` (optional): absent -> a clear
``SourceUnavailable`` that never affects the idealized mode.

HRRR variables span several GRIB level types (isobaric, surface, 2 m, 10 m); this reads the
isobaric group for the 3-D fields and the surface group for terrain/fluxes.  Missing variables
are SKIPPED (never invented; task requirement 8).
"""
from __future__ import annotations

import numpy as np

from .base import require, SourceUnavailable
from .. import thermo, units
from ..internal import AtmosphericState
from ..project import Projection

# HRRR isobaric GRIB shortName -> internal standard name (+ needed conversion)
_ISOBARIC = {"t": "T", "gh": "_gh", "q": "qv", "u": "u", "v": "v", "w": "w", "pres": "p"}


def download(cfg, cache):
    """Fetch the HRRR pressure-level GRIB2 for the case hour from AWS Open Data (no creds)."""
    import requests
    date = cfg.case.date.replace("-", "")
    hour = cfg.case.start_time_utc.split(":")[0].zfill(2)
    key = "%s_t%sz_wrfprs" % (date, hour)
    path = cache.path("hrrr", key, ".grib2")
    cache.require_offline_ok("hrrr", key, ".grib2")
    if cache.has("hrrr", key, ".grib2"):
        return path
    url = ("https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.%s/conus/"
           "hrrr.t%sz.wrfprsf00.grib2" % (date, hour))
    r = requests.get(url, stream=True, timeout=120)
    r.raise_for_status()
    with open(path, "wb") as f:
        for chunk in r.iter_content(1 << 20):
            f.write(chunk)
    cache.record("hrrr", key, path, {"url": url})
    return path


def load(cfg, cache):
    """Read HRRR GRIB2 -> :class:`AtmosphericState` cropped to the case domain (SI)."""
    require("cfgrib", "HRRR reader", "and the ecCodes C library (conda install -c conda-forge eccodes cfgrib)")
    import xarray as xr
    path = download(cfg, cache)
    try:
        iso = xr.open_dataset(path, engine="cfgrib",
                              backend_kwargs={"filter_by_keys": {"typeOfLevel": "isobaricInhPa"}})
    except Exception as e:
        raise SourceUnavailable("could not read HRRR GRIB2 %s: %s" % (path, e))
    dom = cfg.domain
    proj = Projection(dom.center_lat, dom.center_lon, dom.projection)
    lat = np.asarray(iso["latitude"].values); lon = np.asarray(iso["longitude"].values)
    lon = np.where(lon > 180, lon - 360, lon)
    xx, yy = proj.to_xy(lat, lon)
    half = 0.5 * dom.width_km * 1000.0
    mask = (np.abs(xx) <= half) & (np.abs(yy) <= half)
    if not mask.any():
        raise SourceUnavailable("HRRR grid does not cover the requested domain")
    # (this crop/regrid to the model mesh is finished in interpolate.regrid_to_model; here we
    #  hand back the standardised native-grid state for that step)
    st = _to_state(iso, cfg, proj)
    iso.close()
    from ._common import to_height_levels
    return to_height_levels(st)                              # pressure proxy -> geometric height


def _to_state(iso, cfg, proj):
    plev = np.asarray(iso["isobaricInhPa"].values, float)                # hPa, descending
    order = np.argsort(plev)                                             # ascending pressure
    z_proxy = -np.log(plev[order])                                       # monotone vertical coord
    lat = np.asarray(iso["latitude"].values); lon = np.asarray(iso["longitude"].values)
    lon = np.where(lon > 180, lon - 360, lon)
    xx, yy = proj.to_xy(lat, lon)
    ny, nx = lat.shape if lat.ndim == 2 else (lat.size, 1)
    times = np.atleast_1d(np.asarray(iso.get("valid_time", iso.get("time")).values))
    st = AtmosphericState.new(times, z_proxy, yy[:, 0] if yy.ndim == 2 else yy,
                              xx[0] if xx.ndim == 2 else xx,
                              projection=cfg.domain.projection, source="hrrr")
    p_full = np.broadcast_to((plev[order] * 100.0)[None, :, None, None],
                             (times.size, order.size, ny, nx))
    st.add("p", p_full, source="hrrr", original_name="isobaricInhPa", units="Pa")
    for short, name in _ISOBARIC.items():
        if short not in iso or name in ("p",):
            continue
        arr = np.asarray(iso[short].values, float)
        arr = arr[..., order, :, :] if arr.ndim >= 3 else arr
        arr = np.atleast_1d(arr)[None] if arr.ndim == 3 else arr
        if name == "_gh":                                                # geopotential height -> skip (vertical is pressure)
            continue
        st.add(name, arr, source="hrrr", original_name=short, valid_time=str(times[0]))
    if st.has("T") and st.has("p"):
        st.add("theta", thermo.potential_temperature(st.var("T"), st.var("p")),
               source="hrrr", original_name="derived", interpolation_method="theta(T,p)")
    return st


def available():
    from .base import available as _av
    return _av(["cfgrib"])
