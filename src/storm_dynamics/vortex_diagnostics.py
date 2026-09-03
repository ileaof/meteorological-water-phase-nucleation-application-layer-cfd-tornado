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


def tangential_profile(vth2d, r2d, grid, r_max_m=None, nbins=24):
    """Azimuthally-averaged v_theta(r); returns (r_centres, vth_mean, vth_max, core_radius_m)."""
    xp = grid.xp
    r = np.asarray(grid.backend.to_cpu(r2d)).ravel()
    v = np.asarray(grid.backend.to_cpu(vth2d)).ravel()
    r_max = r_max_m if r_max_m is not None else float(np.percentile(r, 40))
    edges = np.linspace(0.0, r_max, nbins + 1)
    idx = np.clip(np.digitize(r, edges) - 1, 0, nbins - 1)
    prof = np.array([v[idx == k].mean() if np.any(idx == k) else 0.0 for k in range(nbins)])
    rc = 0.5 * (edges[:-1] + edges[1:])
    kmax = int(np.argmax(prof))
    return rc, prof, float(prof[kmax]), float(rc[kmax])


def pressure_deficit(p2d, grid, center=None, window=2):
    """Perturbation-pressure deficit = p'_core_min - p'_ambient (domain-edge mean) [Pa] (negative
    for a low-pressure vortex core)."""
    xp = grid.xp
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


def vortex_report(state, grid, z_m=100.0, storm_motion=(0.0, 0.0), radius_m=1500.0):
    """Full vortex diagnosis at level ``z_m`` from a live :class:`FlowState` (never imposed).

    Returns centre (x,y), circulation, peak tangential velocity, core radius, pressure deficit,
    peak vertical vorticity, radial-inflow magnitude and an (approximate) swirl ratio."""
    xp = grid.xp
    k = _level_index(grid, z_m)
    uc, vc, wc = _centered_velocity(state, grid)
    uc2 = uc[:, :, k]; vc2 = vc[:, :, k]
    zeta2 = vertical_vorticity(state, grid)[:, :, k]
    p2 = getattr(state, "p", None)
    p2 = p2[:, :, k] if p2 is not None else None
    center = find_vortex_center(zeta2, grid, p2d=p2)
    gamma = circulation(zeta2, grid, center, radius_m)
    vth, vrad, r = tangential_radial(uc2, vc2, grid, center, storm_motion)
    rc, prof, vth_max, core = tangential_profile(vth, r, grid, r_max_m=radius_m)
    inflow = float(xp.max(-vrad * xp.asarray((np.asarray(grid.backend.to_cpu(r)) <= max(core, grid.dx)))))
    ddef = pressure_deficit(p2, grid, center) if p2 is not None else float("nan")
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
                              border_frac=0.2):
    """Is the vortex SURFACE-CONNECTED or ELEVATED?  Returns V_rot and peak |zeta| on a height
    ladder (interior only), the **surface-to-aloft ratio** (>1 ⇒ surface-intensified/descending,
    <1 ⇒ elevated), the near-surface convergence at the vortex, and — as the spec requires — the
    **first cell-centre height**, so ground contact is never claimed at a height the mesh cannot
    resolve.

    ``border_frac`` sets the interior margin excluded from the search.  It must be generous enough to
    stay outside a nest's boundary-relaxation zone: a too-narrow margin picks up boundary-contaminated
    cells and reports a spuriously huge near-surface V_rot (observed: 37.9 vs the true ~9.0)."""
    xp = grid.xp
    zc = np.asarray(grid.backend.to_cpu(grid.zc))
    uc, vc, wc = _centered_velocity(state, grid)
    zeta3 = vertical_vorticity(state, grid)
    nb = max(2, int(border_frac * grid.nx))
    prof = []
    for zt in z_levels:
        k = int(np.argmin(np.abs(zc - zt)))
        Zi = np.abs(np.asarray(grid.backend.to_cpu(zeta3))[nb:-nb, nb:-nb, k])
        ii, jj = np.unravel_index(int(np.argmax(Zi)), Zi.shape); ii += nb; jj += nb
        u2 = np.asarray(grid.backend.to_cpu(uc))[:, :, k]; v2 = np.asarray(grid.backend.to_cpu(vc))[:, :, k]
        R = max(3, int(radius_m / grid.dx))
        us = u2[max(0, ii - R):ii + R, max(0, jj - R):jj + R]
        vs = v2[max(0, ii - R):ii + R, max(0, jj - R):jj + R]
        vr = float(np.sqrt((us - us.mean()) ** 2 + (vs - vs.mean()) ** 2).max())
        prof.append({"z_m": float(zc[k]), "v_rot_m_s": vr, "zeta_max_s": float(Zi.max()),
                     "edge_cells": int(min(ii, jj, grid.nx - 1 - ii, grid.ny - 1 - jj))})
    v_sfc = prof[0]["v_rot_m_s"]; v_aloft = max(p["v_rot_m_s"] for p in prof)
    # near-surface convergence at the lowest sampled level
    k0 = int(np.argmin(np.abs(zc - z_levels[0])))
    conv = -(grid._central_x(uc[:, :, k0:k0 + 1])[:, :, 0] + grid._central_y(vc[:, :, k0:k0 + 1])[:, :, 0])
    return {
        "first_cell_height_m": float(zc[0]),
        "profile": prof,
        "surface_aloft_ratio": float(v_sfc / v_aloft) if v_aloft > 1e-9 else 0.0,
        "surface_connected": bool(v_sfc >= 0.8 * v_aloft),
        "near_surface_convergence_s": float(xp.max(conv)),
    }


__all__ = [
    "find_vortex_center", "circulation", "tangential_radial", "tangential_profile",
    "pressure_deficit", "vortex_report", "surface_connection_report",
]
