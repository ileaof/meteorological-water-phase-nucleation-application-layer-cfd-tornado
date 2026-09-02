"""Convection initiation -- how the storm is STARTED (never how it rotates).

Three modes are supported across the package: (1) a **warm bubble** (this module + the engine's
Gaussian ``warm_bubble``), (2) **sustained low-level ascent / convergence** forcing
(:mod:`storm_dynamics.forcing`), and (3) **real atmospheric fields** (:mod:`atmospheric_data`).

The idealised thermal perturbation must only *trigger* convection -- it carries NO rotation (it is a
scalar theta' field; any rotation must emerge later from tilting/stretching of the environmental and
baroclinic vorticity).  This module provides the spec's smooth cos^2 bubble with INDEPENDENT
horizontal/vertical radii and multi-bubble support:

    theta'(x,y,z) = dtheta_max * cos^2(pi r / 2) for r < 1,   r^2 = sum((x-xc)/Rx)^2 (+y,+z),

which (unlike a Gaussian) has compact support and a zero-gradient edge.  `grid.xp`-generic.
"""
from __future__ import annotations

import numpy as np


def smooth_bubble(grid, dtheta_max=2.0, center=None, Rx=10000.0, Ry=10000.0, Rz=1500.0):
    """A single cos^2 thermal bubble theta' [K] field with independent radii Rx, Ry, Rz [m]."""
    xp = grid.xp
    cx, cy, cz = center if center is not None else (0.5 * grid.Lx, 0.5 * grid.Ly, 0.0)
    X = xp.asarray(grid.xc)[:, None, None]; Y = xp.asarray(grid.yc)[None, :, None]
    Z = xp.asarray(grid.zc)[None, None, :]
    r = xp.sqrt(((X - cx) / Rx) ** 2 + ((Y - cy) / Ry) ** 2 + ((Z - cz) / Rz) ** 2)
    return xp.where(r < 1.0, dtheta_max * xp.cos(0.5 * xp.pi * r) ** 2, 0.0)


def multi_bubble(grid, specs):
    """Sum several cos^2 bubbles.  ``specs`` is a list of dicts with keys among
    ``dtheta_max, center, Rx, Ry, Rz`` (defaults as :func:`smooth_bubble`)."""
    xp = grid.xp
    out = grid.zeros_c()
    for s in specs:
        out = out + smooth_bubble(grid, **s)
    return out


def apply_bubble_to_state(state, grid, **kw):
    """Add a cos^2 bubble to ``state.theta`` in place (idealised trigger; no rotation)."""
    state.theta = state.theta + smooth_bubble(grid, **kw)
    return state


__all__ = ["smooth_bubble", "multi_bubble", "apply_bubble_to_state"]
