"""Limited-area lateral boundary conditions -- a Davies (1976) relaxation zone (ROADMAP §3a).

A limited-area run is not periodic: its lateral boundaries must be tied to an **external
driving state** (a coarser model, or a reanalysis field), or interior disturbances reflect off
the walls.  The standard treatment is a *relaxation / Davies zone*: over a band of ``width``
cells at every lateral edge, the prognostic fields are nudged toward the external target with a
weight that is strongest at the outermost cell and ramps to zero in the interior --

    phi <- phi + w(i,j) * dt * (phi_target - phi) ,   w = rate * [max(0, 1 - d/width)]^2 ,

``d`` the distance (in cells) to the nearest lateral edge.  The **AMR nest machinery already
uses this pattern** to tie a nest to its parent (:func:`nesting.relaxation_weight`,
``_relax_to_parent``); here it is exposed for the *outermost* domain driven by a prescribed
environment.  The target may be **time-dependent** (swap it each step, e.g. interpolating
between reanalysis times) -- exactly how a real limited-area forecast ingests its lateral
boundary conditions.  Reanalysis/gridded *ingestion itself* (reading ERA5/HRRR files, building
the time-varying target) is the remaining data-plumbing; this module is the BC operator it
feeds.
"""
from __future__ import annotations

import numpy as np

from meteorological_flow.base_state import BaseState
from meteorological_flow.grid import Grid
from meteorological_flow.state import FlowState


def lateral_relaxation_weight(grid: Grid, width: int = 8, rate: float = 1.0 / 300.0):
    """The Davies-zone nudging weight ``(nx,ny,1)``: ``rate`` [1/s] at the outermost cell,
    ramping quadratically to 0 at ``width`` cells in from any lateral edge; 0 in the interior."""
    xp = grid.xp
    nx, ny = grid.nx, grid.ny
    ix = xp.arange(nx); iy = xp.arange(ny)
    d = xp.minimum(xp.minimum(ix, nx - 1 - ix)[:, None],
                   xp.minimum(iy, ny - 1 - iy)[None, :])       # distance to nearest lateral edge
    w = xp.clip(1.0 - d / max(int(width), 1), 0.0, 1.0) ** 2
    return (float(rate) * w)[:, :, None]


def environment_target(grid: Grid, base: BaseState) -> dict:
    """Build a full-field lateral-BC target from an environmental :class:`BaseState`
    (the driving profile): horizontally-uniform ``u,v,theta,qv`` from the sounding, ``w=0``.
    For a time-dependent driver, rebuild (or interpolate) this each step."""
    xp = grid.xp
    u0 = xp.asarray(np.asarray(base.u0, float)); v0 = xp.asarray(np.asarray(base.v0, float))
    return {
        "u": xp.broadcast_to(u0.reshape(1, 1, -1), grid.u_shape).copy(),
        "v": xp.broadcast_to(v0.reshape(1, 1, -1), grid.v_shape).copy(),
        "w": xp.zeros(grid.w_shape),
        "theta": base.field(base.theta0, grid.center_shape, xp=xp),
        "qv": base.field(base.qv0, grid.center_shape, xp=xp),
    }


def apply_lateral_relaxation(state: FlowState, grid: Grid, target: dict, dt: float,
                             weight=None, width: int = 8, rate: float = 1.0 / 300.0):
    """Nudge the lateral boundary band of ``state`` toward ``target`` (a dict with ``u,v,w,
    theta`` and optional ``qv``) by one step ``dt`` -- the limited-area lateral BC.  Modifies
    ``state`` in place; returns the centre-cell weight (cache it and pass as ``weight`` to skip
    rebuilding).  The centre weight is averaged onto the staggered faces for ``u,v,w``."""
    xp = grid.xp
    wu = weight if weight is not None else lateral_relaxation_weight(grid, width, rate)
    wru = xp.zeros(grid.u_shape); wru[1:-1] = 0.5 * (wu[:-1] + wu[1:]); wru[0] = wu[0]; wru[-1] = wu[-1]
    wrv = xp.zeros(grid.v_shape); wrv[:, 1:-1] = 0.5 * (wu[:, :-1] + wu[:, 1:]); wrv[:, 0] = wu[:, 0]; wrv[:, -1] = wu[:, -1]
    state.u += wru * dt * (target["u"] - state.u)
    state.v += wrv * dt * (target["v"] - state.v)
    state.w += wu * dt * (target["w"] - state.w)
    state.theta += wu * dt * (target["theta"] - state.theta)
    if "qv" in target:
        state.qv = xp.maximum(state.qv + wu * dt * (target["qv"] - state.qv), 0.0)
    return wu


__all__ = ["lateral_relaxation_weight", "environment_target", "apply_lateral_relaxation"]
