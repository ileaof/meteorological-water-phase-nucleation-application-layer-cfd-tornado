"""Vortex-scale diagnostics -- circulation, tangential velocity, the vortex centre, the
perturbation-pressure deficit, core radius and a swirl proxy.

These turn a rotating column into the quantities a tornado is actually judged by, WITHOUT ever
imposing a vortex: everything is measured from the resolved fields.  All `grid.xp`-generic.

Conventions: 2-D horizontal fields at a chosen model level; the azimuthal (tangential) velocity
about a centre ``(xc, yc)`` is ``v_theta = -u sin(phi) + v cos(phi)`` with ``phi=atan2(y-yc,x-xc)``
(the spec's form); the radial velocity is ``v_r = u cos(phi) + v sin(phi)``.  Velocities are made
storm-relative (mean or supplied storm motion removed) so rotation is not masked by translation.
"""
from __future__ import annotations

import numpy as np

from . import scales as _scales
from .rotation import _centered_velocity, vertical_vorticity


def _level_index(grid, z_m):
    z = np.asarray(grid.backend.to_cpu(grid.zc))
    return int(np.argmin(np.abs(z - z_m)))


def find_vortex_center(zeta2d, grid, p2d=None, border_frac=1.0 / 6.0, sign=1):
    """Locate the vortex core on a 2-D level: the interior extremum of ``sign*zeta`` (cyclonic
    ``sign=+1``), optionally snapped to a nearby perturbation-pressure minimum.  Returns
    ``(ic, jc, xc, yc)``."""
    xp = grid.xp
    nx, ny = zeta2d.shape
    b = max(1, int(border_frac * nx))
    sub = (sign * zeta2d)[b:-b, b:-b]
    ii, jj = np.unravel_index(int(grid.backend.to_cpu(xp.argmax(sub))), sub.shape)
    ic, jc = ii + b, jj + b
    if p2d is not None:                      # refine to the local p'-min within a small window
        w = max(2, nx // 12)
        i0, j0 = max(0, ic - w), max(0, jc - w)
        win = p2d[i0:ic + w + 1, j0:jc + w + 1]
        di, dj = np.unravel_index(int(grid.backend.to_cpu(xp.argmin(win))), win.shape)
        ic, jc = i0 + di, j0 + dj
    xc = float(np.asarray(grid.backend.to_cpu(grid.xc))[ic])
    yc = float(np.asarray(grid.backend.to_cpu(grid.yc))[jc])
    return ic, jc, xc, yc


def _phi_r(grid, xc, yc):
    xp = grid.xp
    X = xp.asarray(grid.xc)[:, None]; Y = xp.asarray(grid.yc)[None, :]
    dx = X - xc; dy = Y - yc
    r = xp.sqrt(dx * dx + dy * dy)
    phi = xp.arctan2(dy, dx)
    return r, phi


def circulation(zeta2d, grid, center, radius_m):
    """Circulation Gamma = area-integral of zeta over the disk of ``radius_m`` about ``center``
    (equal to the line integral of v_theta by Stokes) [m^2/s]."""
    xp = grid.xp
    _, _, xc, yc = center if len(center) == 4 else (0, 0, center[0], center[1])
    r, _ = _phi_r(grid, xc, yc)
    dA = grid.dx * grid.dy
    return float(xp.sum(xp.where(r <= radius_m, zeta2d, 0.0)) * dA)


def tangential_radial(uc2d, vc2d, grid, center, storm_motion=(0.0, 0.0)):
    """Storm-relative tangential ``v_theta`` and radial ``v_r`` 2-D fields about ``center``."""
    xp = grid.xp
    _, _, xc, yc = center if len(center) == 4 else (0, 0, center[0], center[1])
    r, phi = _phi_r(grid, xc, yc)
    ur = uc2d - storm_motion[0]; vr = vc2d - storm_motion[1]
    vth = -ur * xp.sin(phi) + vr * xp.cos(phi)
    vrad = ur * xp.cos(phi) + vr * xp.sin(phi)
    return vth, vrad, r


def tangential_profile(vth2d, r2d, grid, r_max_m=None, nbins=None, min_cells_per_bin=3):
    """Azimuthally-averaged v_theta(r); returns (r_centres, vth_mean, vth_max, core_radius_m).

    RESOLUTION SEMANTICS (CHANGED -- see docs/REVIEW_REQUEST.md, the audit of 2026-09-04).
    ``nbins`` was fixed at 24, so the bin width was ``r_max/24`` INDEPENDENT of the mesh: with
    the classifier's hardwired ``radius_m=1500`` that is 62.5 m at every resolution.  Two
    consequences, both observed:

    * the innermost bin contained a single cell -- the vortex centre, where ``r = 0`` and
      ``phi = arctan2(0,0) = 0``, so ``v_theta = -u sin0 + v cos0 = v``, i.e. the raw meridional
      wind at ONE cell.  That number crossed the ``v_theta >= 15 m/s`` tornado gate;
    * ``core_radius_m`` came back as exactly 62.5 m at dx=600 m AND at dx=300 m -- a constant
      that reports the binning, not the vortex.

    Now the bin width is tied to the mesh (>= ``min_cells_per_bin`` cells), radii below one cell
    are excluded (no v_theta is defined at r=0), and EMPTY bins return NaN rather than 0.0 so an
    unsampled radius cannot masquerade as a measured zero.  ``nbins=None`` derives the count.
    """
    r = np.asarray(grid.backend.to_cpu(r2d)).ravel()
    v = np.asarray(grid.backend.to_cpu(vth2d)).ravel()
    dx = float(grid.dx)
    r_max = r_max_m if r_max_m is not None else float(np.percentile(r, 40))
    # a bin must span enough cells to average over; and r=0 carries no tangential direction
    bin_w = max(float(min_cells_per_bin) * dx, r_max / 64.0)
    if nbins is None:
        nbins = max(1, int(np.floor((r_max - dx) / bin_w)))
    edges = np.linspace(dx, dx + nbins * bin_w, nbins + 1)
    keep = r >= dx                                     # exclude the centre cell itself
    r = r[keep]; v = v[keep]
    idx = np.digitize(r, edges) - 1
    prof = np.array([v[idx == k].mean() if np.any(idx == k) else np.nan for k in range(nbins)])
    rc = 0.5 * (edges[:-1] + edges[1:])
    if not np.isfinite(prof).any():
        return rc, prof, float("nan"), float("nan")
    kmax = int(np.nanargmax(prof))
    return rc, prof, float(prof[kmax]), float(rc[kmax])


def pressure_deficit(p2d, grid, center=None, window=2, ambient_frac=None):
    """Perturbation-pressure deficit = p'_core_min - p'_ambient [Pa] (negative for a low-pressure
    vortex core).

    ``ambient_frac`` chooses the reference ring: ``None`` (default) uses the domain-edge mean.  On a
    NEST that edge IS the boundary-relaxation zone, where the imposed inflow leaves a large
    projection pressure -- referencing against it corrupts the deficit.  Passing a fraction (0.2,
    matching the search margin) instead averages an INTERIOR ring at that inset."""
    xp = grid.xp
    if ambient_frac:
        nb = max(1, int(ambient_frac * p2d.shape[0]))
        q = p2d[nb:-nb, nb:-nb]
        ambient = 0.25 * (float(xp.mean(q[0, :])) + float(xp.mean(q[-1, :]))
                          + float(xp.mean(q[:, 0])) + float(xp.mean(q[:, -1])))
    else:
        ambient = 0.25 * (float(xp.mean(p2d[0, :])) + float(xp.mean(p2d[-1, :]))
                          + float(xp.mean(p2d[:, 0])) + float(xp.mean(p2d[:, -1])))
    if center is not None:
        ic, jc = center[0], center[1]
        w = window
        core = p2d[max(0, ic - w):ic + w + 1, max(0, jc - w):jc + w + 1]
        pmin = float(xp.min(core))
    else:
        pmin = float(xp.min(p2d))
    return pmin - ambient


def vortex_report(state, grid, z_m=100.0, storm_motion=(0.0, 0.0), radius_m=1500.0,
                  border_frac=1.0 / 6.0):
    """Full vortex diagnosis at level ``z_m`` from a live :class:`FlowState` (never imposed).

    Returns centre (x,y), circulation, peak tangential velocity, core radius, pressure deficit,
    peak vertical vorticity, radial-inflow magnitude and an (approximate) swirl ratio.

    ``border_frac`` is the interior margin of the centre search: on a NEST it must be wide enough
    to exclude the boundary-relaxation zone (0.2 works), or the centre latches onto a
    boundary artifact and every quantity derived from it is meaningless -- see the same note on
    :func:`surface_connection_report`."""
    xp = grid.xp
    k = _level_index(grid, z_m)
    uc, vc, wc = _centered_velocity(state, grid)
    uc2 = uc[:, :, k]; vc2 = vc[:, :, k]
    zeta2 = vertical_vorticity(state, grid)[:, :, k]
    # Prefer the projection pressure stored by the low-memory solver (``p_dyn``, present on
    # nests); fall back to ``state.p``, which the direct solver sets on a parent.  ``p_dyn`` is
    # kept out of ``state.p`` on purpose -- see _project_anelastic_lowmem: ``p`` feeds P_total
    # and hence the thermodynamics, so writing a nest's boundary-influenced projection pressure
    # there changes the run instead of merely reporting on it.
    p2 = getattr(state, "p_dyn", None)
    if p2 is None:
        p2 = getattr(state, "p", None)
    p2 = p2[:, :, k] if p2 is not None else None
    # A CONSTANT pressure field carries no information -- it means the solver never wrote one
    # (a nest's stale zeros, or a level freshly rebuilt by regrid_nest that has not stepped yet).
    # Report that as "unmeasured" (NaN), never as a deficit of 0.0, which would read as a real
    # measurement and quietly pass through into the classification.
    if p2 is not None and float(xp.max(p2) - xp.min(p2)) == 0.0:
        p2 = None
    center = find_vortex_center(zeta2, grid, p2d=p2, border_frac=border_frac)
    gamma = circulation(zeta2, grid, center, radius_m)
    vth, vrad, r = tangential_radial(uc2, vc2, grid, center, storm_motion)
    rc, prof, vth_max, core = tangential_profile(vth, r, grid, r_max_m=radius_m)
    inflow = float(xp.max(-vrad * xp.asarray((np.asarray(grid.backend.to_cpu(r)) <= max(core, grid.dx)))))
    # reference the deficit against the same interior region the centre search trusts
    ddef = (pressure_deficit(p2, grid, center, ambient_frac=border_frac)
            if p2 is not None else float("nan"))
    swirl = (vth_max / inflow) if inflow > 1e-6 else float("inf")
    return {
        "level_m": float(np.asarray(grid.backend.to_cpu(grid.zc))[k]),
        "center_x_m": center[2], "center_y_m": center[3],
        "circulation_m2_s": gamma,
        "v_theta_max_m_s": vth_max, "core_radius_m": core,
        "radial_inflow_m_s": inflow, "swirl_ratio": swirl,
        "pressure_deficit_Pa": ddef,
        "zeta_core_s": float(zeta2[center[0], center[1]]),
    }


def surface_connection_report(state, grid, storm_motion=(0.0, 0.0),
                              z_levels=(50.0, 200.0, 500.0, 1000.0, 1500.0), radius_m=1000.0,
                              border_frac=0.2, border_m=None, min_cells=3, strict=False):
    """Is the vortex SURFACE-CONNECTED or ELEVATED?  Returns V_rot and peak |zeta| on a height
    ladder (interior only), the **surface-to-aloft ratio** (>1 ⇒ surface-intensified/descending,
    <1 ⇒ elevated), the near-surface convergence at the vortex, and — as the spec requires — the
    **first cell-centre height**, so ground contact is never claimed at a height the mesh cannot
    resolve.

    ``border_frac`` sets the interior margin excluded from the search.  It must be generous enough to
    stay outside a nest's boundary-relaxation zone: a too-narrow margin picks up boundary-contaminated
    cells and reports a spuriously huge near-surface V_rot (observed: 37.9 vs the true ~9.0).
    ``border_m`` (PREFERRED) sets that margin as a physical width in metres instead, so the
    excluded region is the same on every mesh; it overrides ``border_frac`` when given.

    RESOLUTION SEMANTICS (CHANGED -- see docs/REVIEW_REQUEST.md A7).  ``radius_m`` is a PHYSICAL
    sampling radius.  The previous implementation did ``R = max(3, int(radius_m / grid.dx))``,
    which silently substituted a 3-cell window whenever the request was smaller than 3*dx: a
    nominal 400 m radius became 1800 m at dx=600 m and 900 m at dx=300 m, a 2x mismatch in exactly
    the quantity a resolution study measures, biased toward the coarse mesh.  It now REFUSES
    instead: when the radius cannot be spanned by ``min_cells`` cells the V_rot entries are NaN,
    ``valid`` is False and ``surface_connected`` is False (``strict=True`` raises
    :class:`scales.UnderResolvedError`).  The ``resolution`` block reports the requested radius,
    the represented radius, the cell count and the status, so a cross-resolution comparison can be
    CHECKED rather than assumed.  ``zeta_max_s`` is still reported when under-resolved: it is a
    point maximum and does not depend on the sampling radius."""
    xp = grid.xp
    zc = np.asarray(grid.backend.to_cpu(grid.zc))
    uc, vc, wc = _centered_velocity(state, grid)
    zeta3 = vertical_vorticity(state, grid)
    # interior margin: PHYSICAL when border_m is given, else the legacy domain fraction
    if border_m is not None:
        _mb = _scales.cells_for_length(border_m, grid.dx, min_cells=1, name="interior margin",
                                       strict=strict)
        nb = max(2, _mb.cells)
        margin_info = _mb.as_dict()
    else:
        nb = max(2, int(border_frac * grid.nx))
        margin_info = {"name": "interior margin (border_frac)", "requested_m": float(nb * grid.dx),
                       "dx_m": float(grid.dx), "cells": int(nb),
                       "represented_m": float(nb * grid.dx), "min_cells": 2, "resolved": True,
                       "relative_error": 0.0, "status": "ok (domain fraction, NOT a fixed length)"}

    # PHYSICAL sampling radius -- refuses rather than silently substituting 3 cells (A7)
    rad = _scales.cells_for_length(radius_m, grid.dx, min_cells=min_cells,
                                   name="V_rot sampling radius", strict=strict, warn=True)
    R = rad.cells
    prof = []
    for zt in z_levels:
        k = int(np.argmin(np.abs(zc - zt)))
        Zi = np.abs(np.asarray(grid.backend.to_cpu(zeta3))[nb:-nb, nb:-nb, k])
        ii, jj = np.unravel_index(int(np.argmax(Zi)), Zi.shape); ii += nb; jj += nb
        u2 = np.asarray(grid.backend.to_cpu(uc))[:, :, k]; v2 = np.asarray(grid.backend.to_cpu(vc))[:, :, k]
        if rad.resolved:
            us = u2[max(0, ii - R):ii + R, max(0, jj - R):jj + R]
            vs = v2[max(0, ii - R):ii + R, max(0, jj - R):jj + R]
            vr = float(np.sqrt((us - us.mean()) ** 2 + (vs - vs.mean()) ** 2).max())
        else:
            vr = float("nan")      # under-resolved: not comparable, and NOT silently widened
        prof.append({"z_m": float(zc[k]), "v_rot_m_s": vr, "zeta_max_s": float(Zi.max()),
                     "edge_cells": int(min(ii, jj, grid.nx - 1 - ii, grid.ny - 1 - jj))})
    v_sfc = prof[0]["v_rot_m_s"]; v_aloft = max(p["v_rot_m_s"] for p in prof)
    # near-surface convergence at the lowest sampled level
    k0 = int(np.argmin(np.abs(zc - z_levels[0])))
    conv = -(grid._central_x(uc[:, :, k0:k0 + 1])[:, :, 0] + grid._central_y(vc[:, :, k0:k0 + 1])[:, :, 0])
    _ok = bool(rad.resolved and v_sfc == v_sfc and v_aloft == v_aloft)
    return {
        "first_cell_height_m": float(zc[0]),
        "profile": prof,
        "surface_aloft_ratio": (float(v_sfc / v_aloft) if (_ok and v_aloft > 1e-9)
                                else (0.0 if _ok else float("nan"))),
        "surface_connected": bool(_ok and v_sfc >= 0.8 * v_aloft),
        "near_surface_convergence_s": float(xp.max(conv)),
        # discretisation provenance, so a cross-resolution comparison can be CHECKED not assumed
        "resolution": rad.as_dict(),
        "interior_margin": margin_info,
        "valid": _ok,
    }


__all__ = [
    "find_vortex_center", "circulation", "tangential_radial", "tangential_profile",
    "pressure_deficit", "vortex_report", "surface_connection_report",
]
