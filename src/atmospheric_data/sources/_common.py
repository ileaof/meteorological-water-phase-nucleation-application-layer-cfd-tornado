"""Shared reader post-processing (ROADMAP §3a, task-2 refinement).

``to_height_levels`` converts a pressure-level :class:`AtmosphericState` (as produced by the
HRRR/ERA5 readers, whose vertical is a pressure proxy) onto the model's **geometric height**
coordinate, using the hypsometric equation on the domain-mean virtual-temperature profile.
This makes the vertical map height-correct (surface-first, ascending) before regridding.
Per-column height (vs the domain-mean used here) is a further refinement.
"""
from __future__ import annotations

import numpy as np

from .. import thermo


def to_height_levels(state, z_surface=0.0):
    """Return a copy of ``state`` with the vertical coordinate replaced by geometric height
    [m] (ascending, surface first).  Requires ``p`` and ``T`` (and uses ``qv`` if present).
    Reorders every 4-D field to the new (surface-first) level order."""
    ds = state.ds
    if "p" not in ds or "T" not in ds:
        return state                                    # nothing to convert (e.g. synthetic on height)
    p = np.asarray(ds["p"].values, float)               # (t,z,y,x)
    T = np.asarray(ds["T"].values, float)
    qv = np.asarray(ds["qv"].values, float) if "qv" in ds else np.zeros_like(T)
    p_col = np.nanmean(p, axis=(0, 2, 3))               # domain-mean pressure per level
    Tv_col = np.nanmean(thermo.virtual_temperature(T, qv), axis=(0, 2, 3))
    order = np.argsort(-p_col)                          # descending pressure = surface first
    p_sorted = p_col[order]; Tv_sorted = Tv_col[order]
    height = thermo.hypsometric_height(p_sorted, Tv_sorted, z0=z_surface)   # ascending
    # rebuild the state on height levels, reordering all z-indexed variables
    import xarray as xr
    new = ds.isel(z=order).assign_coords(z=("z", height))
    new["z"].attrs.update(units="m", long_name="geometric height", axis="Z")
    state.ds = new
    state.provenance.setdefault("vars", {})
    state.provenance["vertical_conversion"] = "pressure->height (hypsometric, domain-mean Tv)"
    return state
