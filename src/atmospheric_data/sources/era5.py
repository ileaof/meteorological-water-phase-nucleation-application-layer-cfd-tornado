"""ERA5 / Copernicus ingestion (ROADMAP §3a, source B): the global synoptic environment.

Reading a **NetCDF** ERA5 file needs only ``xarray`` (always available); **downloading** needs
``cdsapi`` + the user's CDS credentials (``~/.cdsapirc`` -- never stored in the repo).  ERA5
gives the synoptic setting, NOT the tornado vortex (documented limitation).  Vertical velocity
in ERA5 is pressure-velocity omega [Pa/s]; converted to w [m/s] via ``w = -omega/(rho g)``.
"""
from __future__ import annotations

import numpy as np

from .base import require, SourceUnavailable, optional_import
from .. import thermo, units
from ..internal import AtmosphericState
from ..project import Projection

# ERA5 short name -> internal standard name
_MAP = {"t": "T", "q": "qv", "u": "u", "v": "v"}


def download(cfg, cache):
    """Retrieve ERA5 pressure levels via the CDS API (needs cdsapi + ~/.cdsapirc)."""
    cdsapi = require("cdsapi", "ERA5 download",
                     "and put your key in ~/.cdsapirc (https://cds.climate.copernicus.eu)")
    key = "%s_%s" % (cfg.case.date, cfg.case.start_time_utc.replace(":", ""))
    path = cache.path("era5", key, ".nc")
    cache.require_offline_ok("era5", key, ".nc")
    if cache.has("era5", key, ".nc"):
        return path
    dom = cfg.domain
    dlat = dom.height_km / 111.0; dlon = dom.width_km / (111.0 * np.cos(np.radians(dom.center_lat)))
    c = cdsapi.Client()
    c.retrieve("reanalysis-era5-pressure-levels", {
        "product_type": "reanalysis", "format": "netcdf",
        "variable": ["temperature", "specific_humidity", "u_component_of_wind",
                     "v_component_of_wind", "vertical_velocity", "geopotential"],
        "pressure_level": ["1000", "925", "850", "700", "500", "400", "300", "250",
                           "200", "150", "100"],
        "year": cfg.case.date[:4], "month": cfg.case.date[5:7], "day": cfg.case.date[8:10],
        "time": cfg.case.start_time_utc,
        "area": [dom.center_lat + dlat, dom.center_lon - dlon,
                 dom.center_lat - dlat, dom.center_lon + dlon]}, path)
    cache.record("era5", key, path)
    return path


def load(cfg, cache):
    """Read an ERA5 NetCDF -> :class:`AtmosphericState` (SI)."""
    import xarray as xr
    key = "%s_%s" % (cfg.case.date, cfg.case.start_time_utc.replace(":", ""))
    path = cache.path("era5", key, ".nc")
    if not cache.has("era5", key, ".nc"):
        path = download(cfg, cache)
    ds = xr.open_dataset(path)
    dom = cfg.domain
    proj = Projection(dom.center_lat, dom.center_lon, dom.projection)
    lat = np.asarray(ds["latitude"].values); lon = np.asarray(ds["longitude"].values)
    lon = np.where(lon > 180, lon - 360, lon)
    LAT, LON = np.meshgrid(lat, lon, indexing="ij")
    xx, yy = proj.to_xy(LAT, LON)
    plev_name = "level" if "level" in ds.dims else ("pressure_level" if "pressure_level" in ds.dims else "isobaricInhPa")
    plev = np.asarray(ds[plev_name].values, float)
    order = np.argsort(plev)
    times = np.atleast_1d(np.asarray(ds[("valid_time" if "valid_time" in ds else "time")].values))
    st = AtmosphericState.new(times, -np.log(plev[order]), yy[:, 0], xx[0],
                              projection=dom.projection, source="era5")
    def field(nm):
        a = np.asarray(ds[nm].values, float)
        while a.ndim < 4:
            a = a[None]
        return a[:, order]
    p_full = np.broadcast_to((plev[order] * 100.0)[None, :, None, None], field("t").shape)
    st.add("p", p_full, source="era5", original_name=plev_name, units="Pa")
    for short, name in _MAP.items():
        if short in ds:
            st.add(name, field(short), source="era5", original_name=short, valid_time=str(times[0]))
    if "w" in ds:                                            # omega [Pa/s] -> w [m/s]
        omega = field("w")
        Tv = thermo.virtual_temperature(st.var("T"), st.var("qv")) if st.has("qv") else st.var("T")
        rho = thermo.density(p_full, Tv)
        st.add("w", -omega / (rho * thermo.g0), source="era5", original_name="w (omega)",
               interpolation_method="w=-omega/(rho g)")
    if st.has("T"):
        st.add("theta", thermo.potential_temperature(st.var("T"), p_full), source="era5",
               original_name="derived", interpolation_method="theta(T,p)")
    ds.close()
    from ._common import to_height_levels
    return to_height_levels(st)                              # pressure proxy -> geometric height


def available():
    return True                                              # NetCDF read needs only xarray
