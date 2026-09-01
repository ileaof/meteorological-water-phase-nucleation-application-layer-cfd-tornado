"""Initial + lateral/top boundary conditions from real data (ROADMAP §3a).

Writes the CF-NetCDF products the real_case driver consumes:
``initial_conditions.nc`` (full model-grid state at t0), ``boundary_{west,east,south,north,
top}.nc`` (edge time series for the Davies relaxation zone -- see
:mod:`storm_dynamics.limited_area`), and ``surface_forcing.nc`` (terrain + surface fluxes).
Lateral boundaries carry a TIME axis so the nudging can interpolate between analysis times.
"""
from __future__ import annotations

import os

import numpy as np

_STATE_VARS = ("u", "v", "w", "theta", "qv", "p")


def _ds(coords, data_vars, attrs):
    import xarray as xr
    return xr.Dataset(data_vars, coords=coords, attrs=attrs)


def write_initial_conditions(fields, x, y, z, path, it=0, meta=None):
    """Write the initial model-grid state (cell-centred u,v,w,theta,qv,p) to NetCDF."""
    import xarray as xr
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    dv = {}
    for nm in _STATE_VARS:
        if nm in fields:
            dv[nm] = (("z", "y", "x"), np.asarray(fields[nm][it], float))
    ds = xr.Dataset(dv, coords={"z": z, "y": y, "x": x},
                    attrs={"Conventions": "CF-1.8", "title": "real-case initial conditions",
                           **(meta or {})})
    ds.to_netcdf(path)
    return path


def write_boundaries(fields, x, y, z, times, outdir, meta=None):
    """Write the five boundary files (time series of edge slices) for the relaxation zone."""
    import xarray as xr
    os.makedirs(outdir, exist_ok=True)
    edges = {"west": (slice(None), slice(None), 0), "east": (slice(None), slice(None), -1),
             "south": (slice(None), 0, slice(None)), "north": (slice(None), -1, slice(None)),
             "top": (-1, slice(None), slice(None))}
    paths = {}
    for name, sl in edges.items():
        dv = {}
        for nm in _STATE_VARS:
            if nm not in fields:
                continue
            arr = np.asarray(fields[nm])[(slice(None),) + sl]              # (time, ...) edge slice
            if name in ("west", "east"):
                dims = ("time", "z", "y")
            elif name in ("south", "north"):
                dims = ("time", "z", "x")
            else:
                dims = ("time", "y", "x")
            dv[nm] = (dims, arr)
        coords = {"time": times, "z": z, "y": y, "x": x}
        ds = xr.Dataset(dv, coords={k: coords[k] for k in set(sum([list(v[0]) for v in dv.values()], []))},
                        attrs={"Conventions": "CF-1.8", "boundary": name, **(meta or {})})
        p = os.path.join(outdir, "boundary_%s.nc" % name)
        ds.to_netcdf(p); paths[name] = p
    return paths


def write_surface_forcing(terrain, x, y, path, fluxes=None, meta=None):
    """Write terrain (+ optional surface fluxes) on the MODEL grid to surface_forcing.nc.

    ``terrain`` is a ``(ny,nx)`` array already regridded to the model axes ``x,y`` (use
    :func:`interpolate.regrid_surface`)."""
    import xarray as xr
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    dv = {}
    if terrain is not None:
        dv["terrain"] = (("y", "x"), np.asarray(terrain, float))
    for nm, arr in (fluxes or {}).items():
        dv[nm] = (("y", "x"), np.asarray(arr, float))
    ds = xr.Dataset(dv, coords={"y": y, "x": x},
                    attrs={"Conventions": "CF-1.8", "title": "surface forcing", **(meta or {})})
    ds.to_netcdf(path)
    return path
