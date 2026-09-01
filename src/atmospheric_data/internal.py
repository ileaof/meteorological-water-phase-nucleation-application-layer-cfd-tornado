"""The unified internal format (ROADMAP §3a): a CF-conventions NetCDF4 dataset.

Every source (HRRR, ERA5, NEXRAD, soundings, ...) is converted into ONE standard container --
an :class:`xarray.Dataset` with dimensions ``(time, z, y, x)``, standardised variable names,
SI units, and provenance attributes on every variable.  Downstream (interpolation, base state,
IC/BC, QC) reads only this container, so adding a source never touches the physics.

Standard variable names (task spec): ``T p rho theta theta_v qv qc qr qi qs qg u v w terrain
reflectivity``.  Per-variable attributes: ``units long_name source original_name
interpolation_method valid_time projection missing_value``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .units import SI_UNITS

# standard internal variables -> (long_name)
STANDARD_VARS = {
    "T": "air temperature", "p": "air pressure", "rho": "air density",
    "theta": "potential temperature", "theta_v": "virtual potential temperature",
    "qv": "water vapour specific humidity", "qc": "cloud liquid water mixing ratio",
    "qr": "rain water mixing ratio", "qi": "cloud ice mixing ratio",
    "qs": "snow mixing ratio", "qg": "graupel mixing ratio",
    "u": "eastward wind", "v": "northward wind", "w": "upward air velocity",
    "terrain": "surface altitude", "reflectivity": "radar reflectivity",
}
DIMS = ("time", "z", "y", "x")


@dataclass
class AtmosphericState:
    """A CF-NetCDF container of standardised atmospheric fields on ``(time, z, y, x)``.

    Build incrementally with :meth:`add`, then :meth:`to_netcdf`.  Requires an xarray backend;
    the class is a thin, well-documented wrapper so the rest of the pipeline is xarray-native."""
    ds: "object" = None                 # xarray.Dataset
    provenance: dict = field(default_factory=dict)

    @classmethod
    def new(cls, time, z, y, x, *, projection="none", source="unknown"):
        """Create an empty state on the given coordinate arrays (SI: z,y,x in m; time as
        datetime64 or seconds).  ``projection`` names the horizontal projection used for x,y."""
        import xarray as xr
        ds = xr.Dataset(coords={"time": np.atleast_1d(time), "z": np.asarray(z, float),
                                "y": np.asarray(y, float), "x": np.asarray(x, float)})
        for c, u in (("z", "m"), ("y", "m"), ("x", "m")):
            ds[c].attrs.update(units=u, long_name={"z": "height", "y": "y (projected)",
                                                   "x": "x (projected)"}[c], axis=c.upper())
        ds.attrs.update(Conventions="CF-1.8", source=source, projection=projection,
                        title="unified atmospheric state (met_h2o real-case ingestion)")
        return cls(ds=ds, provenance={"source": source, "projection": projection, "vars": {}})

    def add(self, name, data, *, source, original_name="", interpolation_method="none",
            valid_time="", missing_value=np.nan, units=None, long_name=None):
        """Add a standard variable (name in :data:`STANDARD_VARS`).  ``data`` shape must be
        broadcastable to ``(time, z, y, x)`` (surface fields may omit z -> stored as (time,y,x))."""
        import xarray as xr
        if name not in STANDARD_VARS:
            raise KeyError("unknown standard variable %r; allowed: %s"
                           % (name, sorted(STANDARD_VARS)))
        data = np.asarray(data, float)
        if data.ndim == 4:
            dims = DIMS
        elif data.ndim == 3:
            dims = ("time", "y", "x")            # a time-varying surface field
        elif data.ndim == 2:
            dims = ("y", "x")                    # a time-invariant surface field (e.g. terrain)
        elif data.ndim == 1 and data.shape[0] == self.ds.sizes["z"]:
            dims = ("z",)                         # a 1-D profile (soundings)
        else:
            raise ValueError("%s: shape %s not (t,z,y,x)/(t,y,x)/(z,)" % (name, data.shape))
        attrs = dict(units=units or SI_UNITS.get(name, ""), long_name=long_name or STANDARD_VARS[name],
                     source=str(source), original_name=str(original_name),
                     interpolation_method=str(interpolation_method), valid_time=str(valid_time),
                     projection=self.ds.attrs.get("projection", "none"),
                     missing_value=float(missing_value))
        self.ds[name] = xr.DataArray(data, dims=dims, attrs=attrs)
        self.provenance.setdefault("vars", {})[name] = {"source": str(source),
                                                        "original_name": str(original_name),
                                                        "interpolation_method": str(interpolation_method)}
        return self

    def has(self, name):
        return self.ds is not None and name in self.ds

    def var(self, name):
        return np.asarray(self.ds[name].values, float)

    def to_netcdf(self, path):
        """Write CF-NetCDF4.  Falls back to the scipy backend if netCDF4/h5netcdf absent."""
        import os
        os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
        try:
            self.ds.to_netcdf(path, format="NETCDF4")
        except Exception:
            self.ds.to_netcdf(path)             # let xarray pick any available backend
        return path

    @classmethod
    def from_netcdf(cls, path):
        import xarray as xr
        ds = xr.open_dataset(path)
        return cls(ds=ds, provenance={"source": ds.attrs.get("source", "unknown"),
                                       "projection": ds.attrs.get("projection", "none")})

    def summary(self):
        """A compact dict of what's present -- variables, shapes, units, sources."""
        if self.ds is None:
            return {}
        return {"dims": dict(self.ds.sizes),
                "projection": self.ds.attrs.get("projection", "none"),
                "variables": {n: {"shape": tuple(self.ds[n].shape),
                                  "units": self.ds[n].attrs.get("units", ""),
                                  "source": self.ds[n].attrs.get("source", "")}
                              for n in self.ds.data_vars}}
