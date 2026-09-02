"""Cold-pool & downdraft diagnostics -- the low-level BAROCLINIC engine of tornadogenesis.

The evaporative cold pool exists in the model as a *process* (rain evaporation cools the sub-cloud
air); this module measures it: the virtual-potential-temperature deficit, the buoyancy, the
density-current intensity ``C``, the cold-pool footprint, the gust-front strength (the leading
theta_v gradient / convergence), and the downdraft field.  These set the horizontal density
gradient that generates low-level horizontal vorticity baroclinically.  All `grid.xp`-generic.

The rear-flank downdraft is *not* forced cold here -- its thermal deficit is whatever the resolved
thermodynamics + microphysics produce, and is simply reported.
"""
from __future__ import annotations

import numpy as np

from .rotation import _centered_velocity

_G = 9.81


def virtual_potential_temperature(state, grid):
    """theta_v = theta (1 + 0.61 q_v - q_condensate) at cell centres [K]."""
    xp = grid.xp
    qc = state.qv * 0.0
    for nm in ("ql", "qi", "qr", "qs", "qg", "qh"):
        a = getattr(state, nm, None)
        if a is not None:
            qc = qc + a
    return state.theta * (1.0 + 0.61 * state.qv - qc)


def _hmean(a, xp):
    return xp.mean(a, axis=(0, 1), keepdims=True)


def coldpool_buoyancy(state, grid, thetav=None):
    """Buoyancy b = g * theta_v' / <theta_v> relative to the horizontal mean at each level [m/s^2]
    (negative in the cold pool).  Returns (b, thetav_pert)."""
    xp = grid.xp
    tv = virtual_potential_temperature(state, grid) if thetav is None else thetav
    tv0 = _hmean(tv, xp)
    tvp = tv - tv0
    return _G * tvp / tv0, tvp


def coldpool_intensity(state, grid, z_top_m=3000.0):
    """Cold-pool intensity C = sqrt(2 * integral_0^h max(-b,0) dz) [m/s] as a 2-D field -- the
    theoretical density-current propagation speed; take ``.max()`` for the storm's cold pool."""
    xp = grid.xp
    b, _ = coldpool_buoyancy(state, grid)
    z = np.asarray(grid.backend.to_cpu(grid.zc))
    mask = xp.asarray((z < z_top_m), dtype=float)[None, None, :]
    dz = grid.dz if not getattr(grid, "stretched", False) else grid.dz_c[None, None, :]
    integ = xp.sum(xp.maximum(-b, 0.0) * mask * dz, axis=2)
    return xp.sqrt(2.0 * integ)


def coldpool_mask(state, grid, thresh_K=-1.0, z_m=100.0):
    """Boolean 2-D footprint where the near-surface theta_v' is below ``thresh_K`` [K]."""
    xp = grid.xp
    _, tvp = coldpool_buoyancy(state, grid)
    k = int(np.argmin(np.abs(np.asarray(grid.backend.to_cpu(grid.zc)) - z_m)))
    return tvp[:, :, k] < thresh_K


def gust_front_strength(state, grid, z_m=100.0):
    """Peak horizontal |grad theta_v'| near the surface [K/m] -- the leading edge / gust front,
    and the peak low-level convergence -(du/dx+dv/dy) [1/s] there."""
    xp = grid.xp
    _, tvp = coldpool_buoyancy(state, grid)
    k = int(np.argmin(np.abs(np.asarray(grid.backend.to_cpu(grid.zc)) - z_m)))
    tv2 = tvp[:, :, k]
    gx = grid._central_x(tv2[:, :, None])[:, :, 0]; gy = grid._central_y(tv2[:, :, None])[:, :, 0]
    grad = xp.sqrt(gx * gx + gy * gy)
    uc, vc, _ = _centered_velocity(state, grid)
    conv = -(grid._central_x(uc[:, :, k:k + 1])[:, :, 0] + grid._central_y(vc[:, :, k:k + 1])[:, :, 0])
    return float(xp.max(grad)), float(xp.max(conv))


def coldpool_report(state, grid, z_m=100.0):
    """Compact cold-pool / downdraft summary for the JSON report."""
    xp = grid.xp
    _, tvp = coldpool_buoyancy(state, grid)
    C = coldpool_intensity(state, grid)
    _, _, wc = _centered_velocity(state, grid)
    k = int(np.argmin(np.abs(np.asarray(grid.backend.to_cpu(grid.zc)) - z_m)))
    tv_low = tvp[:, :, k]
    cold = tv_low < -1.0
    area_frac = float(xp.mean(xp.asarray(cold, dtype=float)))
    gf_grad, gf_conv = gust_front_strength(state, grid, z_m)
    return {
        "coldpool_min_thetav_pert_K": float(xp.min(tv_low)),
        "coldpool_intensity_C_m_s": float(xp.max(C)),
        "coldpool_area_fraction": area_frac,
        "gust_front_thetav_grad_K_m": gf_grad,
        "gust_front_convergence_s": gf_conv,
        "downdraft_w_min_m_s": float(xp.min(wc)),
        "downdraft_area_fraction": float(xp.mean(xp.asarray(wc[:, :, k] < -1.0, dtype=float))),
    }


__all__ = [
    "virtual_potential_temperature", "coldpool_buoyancy", "coldpool_intensity",
    "coldpool_mask", "gust_front_strength", "coldpool_report",
]
