"""Physical length <-> cell-count conversion, with explicit under-resolution reporting.

**Why this module exists.**  Three independent defects in this project all had the same shape: a
quantity that is physically a LENGTH was stored or converted as a CELL COUNT, so it silently
rescaled with the mesh and corrupted every comparison made across two resolutions.

1. ``NestSpec.relax_width = 4`` cells -- a 267 m Davies sponge at dx=67 m became 89 m at
   dx=22 m, so the finer nest damped incoming boundary error ~3x less.
2. ``R = max(3, int(radius_m / grid.dx))`` in ``surface_connection_report`` -- a requested 400 m
   diagnostic radius silently became 1800 m at dx=600 m and 900 m at dx=300 m (the ``max(3, ...)``
   floor), a 2x mismatch in exactly the quantity a resolution study was measuring.
3. Per-level ``radius_m`` / interior-margin conventions that differed between cascade levels,
   making V_rot incomparable across the nest ladder.

None was visible in code review; each appeared only when a physical quantity was compared across
two meshes.  The rule this module enforces:

    A physical scale is requested in METRES.  The discretisation is reported, never hidden.
    A scale that cannot be represented is an EXPLICIT failure, never a silent substitution.

**Design decisions.**

* ``floor``, not ``round``: a represented region must never be LARGER than requested.  Rounding
  up would silently enlarge a measurement disk by up to dx/2, which is precisely how defect (2)
  inflated the coarse-grid V_rot.
* A minimum cell count is a *numerical support* requirement (you cannot estimate a gradient or a
  rotational velocity from one cell), and it is reported separately from the physical request.
  When the request cannot meet it the result is ``resolved=False`` -- the caller decides whether
  to return NaN, warn, or raise.  This module never substitutes a different physical scale.
* ``strict=True`` raises, for call sites where silently continuing is never acceptable.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import warnings


class UnderResolvedError(ValueError):
    """A requested physical scale cannot be represented on this mesh."""


@dataclass(frozen=True)
class PhysicalScale:
    """The result of asking for ``requested_m`` on a mesh of spacing ``dx_m``.

    ``cells`` is the number of whole cells spanning the request (floor -- never enlarged).
    ``represented_m = cells * dx_m`` is what the mesh actually delivers.
    ``resolved`` is False when ``cells < min_cells``, i.e. there is not enough numerical support;
    the caller must then decide what to do -- this object never silently substitutes a value.
    """
    name: str
    requested_m: float
    dx_m: float
    cells: int
    represented_m: float
    min_cells: int
    resolved: bool
    relative_error: float

    @property
    def status(self) -> str:
        return "ok" if self.resolved else "under_resolved"

    def __bool__(self) -> bool:
        return self.resolved

    def describe(self) -> str:
        return ("%s: requested %.4g m on dx=%.4g m -> %d cells = %.4g m "
                "(rel. err %+.1f%%, min_cells=%d, %s)"
                % (self.name, self.requested_m, self.dx_m, self.cells, self.represented_m,
                   100.0 * self.relative_error, self.min_cells, self.status))

    def as_dict(self) -> dict:
        return {"name": self.name, "requested_m": self.requested_m, "dx_m": self.dx_m,
                "cells": self.cells, "represented_m": self.represented_m,
                "min_cells": self.min_cells, "resolved": self.resolved,
                "relative_error": self.relative_error, "status": self.status}


def cells_for_length(requested_m, dx_m, min_cells: int = 3, name: str = "length",
                     strict: bool = False, warn: bool = False) -> PhysicalScale:
    """Convert a physical length to a whole number of cells, reporting what was achieved.

    Parameters
    ----------
    requested_m : float   the PHYSICAL scale being asked for [m]
    dx_m        : float   mesh spacing [m]
    min_cells   : int     numerical support required for the quantity being computed (e.g. a
                          rotational velocity needs a few cells across the sampling window).
                          This is a support requirement, NOT a physical scale: falling below it
                          sets ``resolved=False`` rather than substituting a bigger radius.
    strict      : bool    raise :class:`UnderResolvedError` instead of returning unresolved.
    warn        : bool    emit a :class:`RuntimeWarning` when unresolved.

    Notes
    -----
    Uses ``floor`` so the represented region is never larger than requested.  A request smaller
    than one cell yields ``cells=0`` and ``resolved=False`` -- it is not rounded up to 1.
    """
    dx_m = float(dx_m)
    requested_m = float(requested_m)
    if not (dx_m > 0.0):
        raise ValueError("dx_m must be positive, got %r" % (dx_m,))
    if requested_m < 0.0:
        raise ValueError("requested_m must be non-negative, got %r" % (requested_m,))
    cells = int(math.floor(requested_m / dx_m))
    represented = cells * dx_m
    resolved = cells >= int(min_cells)
    rel = (represented - requested_m) / requested_m if requested_m > 0.0 else 0.0
    sc = PhysicalScale(name=str(name), requested_m=requested_m, dx_m=dx_m, cells=cells,
                       represented_m=represented, min_cells=int(min_cells), resolved=resolved,
                       relative_error=rel)
    if not resolved:
        msg = ("%s is UNDER-RESOLVED: %.4g m needs >= %d cells at dx=%.4g m but spans only %d. "
               "The measurement is not physically comparable at this resolution; it has NOT been "
               "silently widened. Use a coarser-mesh-compatible scale (>= %.4g m) or a finer mesh."
               % (sc.name, requested_m, min_cells, dx_m, cells, min_cells * dx_m))
        if strict:
            raise UnderResolvedError(msg)
        if warn:
            warnings.warn(msg, RuntimeWarning, stacklevel=2)
    return sc


def smallest_resolvable_length(dx_m, min_cells: int = 3) -> float:
    """The smallest physical scale this mesh can support at ``min_cells`` -- i.e. the floor a
    cross-resolution comparison must respect.  For a study spanning several dx, the comparison
    scale must be >= ``max(smallest_resolvable_length(dx) for dx in all_meshes)``."""
    return float(min_cells) * float(dx_m)


def common_comparison_length(requested_m, dx_list, min_cells: int = 3,
                             name: str = "comparison length") -> PhysicalScale:
    """The scale to use when the SAME physical region must be measured on several meshes.

    Returns the request evaluated on the COARSEST mesh, which is the binding constraint.  If that
    is under-resolved the whole comparison is invalid -- refusing here is the point, because this
    is exactly the situation that produced "1800 m vs 900 m" while both call sites said 400 m.
    """
    dxs = [float(d) for d in dx_list]
    if not dxs:
        raise ValueError("dx_list is empty")
    return cells_for_length(requested_m, max(dxs), min_cells=min_cells, name=name)


__all__ = ["PhysicalScale", "UnderResolvedError", "cells_for_length",
           "smallest_resolvable_length", "common_comparison_length"]
