"""Rotation diagnostics (item 7) -- the numbers that prove the storm rotated.

Given a :class:`~meteorological_flow.state.FlowState` this computes the 3-D
vorticity, the **vertical vorticity** zeta = dv/dx - du/dy (the mesocyclone /
tornado indicator), **updraft helicity** UH = integral of w*zeta over a mid-level
layer (the standard supercell rotation-track measure), and column trackers for
the **mid-level mesocyclone** (max zeta at ~3-6 km) and the **near-surface**
zeta (max at ~0-1 km -- the tornadogenesis proxy).  Environmental storm-relative
helicity (SRH) and bulk shear come from the base-state hodograph
(:mod:`storm_dynamics.soundings`).

Vorticity is evaluated at cell centres from the centred velocity via the grid's
central-difference stencils (periodic-aware in x, y; one-sided at the z walls) --
adequate for diagnostics; the prognostic core keeps the fields on the C-grid.
"""
from __future__ import annotations

import numpy as np

from meteorological_flow.grid import Grid
from meteorological_flow.state import FlowState


def _centered_velocity(state, grid):
    uc = 0.5 * (state.u[:-1, :, :] + state.u[1:, :, :])
    vc = 0.5 * (state.v[:, :-1, :] + state.v[:, 1:, :])
    wc = 0.5 * (state.w[:, :, :-1] + state.w[:, :, 1:])
    return uc, vc, wc


def vorticity_3d(state: FlowState, grid: Grid):
    """Return (xi, eta, zeta) vorticity components at cell centres [1/s]."""
    uc, vc, wc = _centered_velocity(state, grid)
    dwdy = grid._central_y(wc); dvdz = grid._central_z(vc)
    dudz = grid._central_z(uc); dwdx = grid._central_x(wc)
    dvdx = grid._central_x(vc); dudy = grid._central_y(uc)
    xi = dwdy - dvdz           # x-vorticity
    eta = dudz - dwdx          # y-vorticity
    zeta = dvdx - dudy         # vertical vorticity
    return xi, eta, zeta


def vertical_vorticity(state: FlowState, grid: Grid):
    """Vertical vorticity zeta = dv/dx - du/dy at cell centres [1/s]."""
    uc, vc, _ = _centered_velocity(state, grid)
    return grid._central_x(vc) - grid._central_y(uc)


def _layer_mask(grid, z_lo, z_hi):
    z = grid.backend.to_cpu(grid.zc)
    return (z >= z_lo) & (z < z_hi)


def updraft_helicity(state: FlowState, grid: Grid, z_lo=2000.0, z_hi=5000.0):
    """Updraft helicity UH = integral_{z_lo}^{z_hi} w*zeta dz (w>0, zeta>0) [m^2/s^2].

    Returns the 2-D field over (x, y); take ``.max()`` for the storm's rotation
    track intensity.
    """
    xp = grid.xp
    _, _, wc = _centered_velocity(state, grid)
    zeta = vertical_vorticity(state, grid)
    dz = grid.dz if not getattr(grid, "stretched", False) else grid.dz_c[None, None, :]
    integrand = xp.where((wc > 0.0) & (zeta > 0.0), wc * zeta, 0.0)
    mask = xp.asarray(_layer_mask(grid, z_lo, z_hi), dtype=float)[None, None, :]
    return xp.sum(integrand * mask * dz, axis=2)


def _max_zeta_in_layer(zeta, grid, z_lo, z_hi):
    xp = grid.xp
    mask = _layer_mask(grid, z_lo, z_hi)
    if not mask.any():
        return 0.0
    sub = zeta[:, :, xp.asarray(mask)] if hasattr(mask, "shape") else zeta
    zmask = xp.asarray(mask)
    layer = zeta[:, :, zmask]
    return float(xp.max(xp.abs(layer))) if layer.size else 0.0


def rotation_report(state: FlowState, grid: Grid, base=None) -> dict:
    """Compact rotation summary for the JSON report / regression tests."""
    xp = grid.xp
    xi, eta, zeta = vorticity_3d(state, grid)
    _, _, wc = _centered_velocity(state, grid)
    uh = updraft_helicity(state, grid)
    # signed mid-level extrema (storm splitting -> a cyclonic AND an anticyclonic core)
    ml = _layer_mask(grid, 3000.0, 6000.0)
    zmid = zeta[:, :, xp.asarray(ml)] if ml.any() else zeta
    ll = _layer_mask(grid, 0.0, 1000.0)
    zlow = zeta[:, :, xp.asarray(ll)] if ll.any() else zeta
    out = {
        "zeta_max": float(xp.max(zeta)),
        "zeta_min": float(xp.min(zeta)),
        "zeta_abs_max": float(xp.max(xp.abs(zeta))),
        "midlevel_zeta_max": float(xp.max(zmid)) if zmid.size else 0.0,
        "midlevel_zeta_min": float(xp.min(zmid)) if zmid.size else 0.0,
        "midlevel_mesocyclone": float(xp.max(xp.abs(zmid))) if zmid.size else 0.0,
        "near_surface_zeta_max": float(xp.max(zlow)) if zlow.size else 0.0,
        "near_surface_zeta_abs_max": float(xp.max(xp.abs(zlow))) if zlow.size else 0.0,
        "updraft_helicity_max": float(xp.max(uh)),          # 2-5 km (legacy key)
        "updraft_helicity_2_5km": float(xp.max(uh)),
        "updraft_helicity_0_1km": float(xp.max(updraft_helicity(state, grid, 0.0, 1000.0))),
        "updraft_helicity_0_3km": float(xp.max(updraft_helicity(state, grid, 0.0, 3000.0))),
        "w_max": float(xp.max(wc)),
        "w_min": float(xp.min(wc)),
        "horiz_vort_max": float(xp.max(xp.sqrt(xi**2 + eta**2))),
    }
    if base is not None:
        from . import soundings as snd
        out["env_SRH_0_1km"] = snd.storm_relative_helicity(base, z_top=1000.0)
        out["env_SRH_0_3km"] = snd.storm_relative_helicity(base, z_top=3000.0)
        out["env_shear_0_1km"] = snd.bulk_shear(base, 0.0, 1000.0)
        out["env_shear_0_6km"] = snd.bulk_shear(base, 0.0, 6000.0)
        out["env_BRN"] = snd.bulk_richardson_number(base)
    return out


class TornadogenesisTracker:
    """Track the near-surface zeta maximum (the tornadogenesis proxy) over time,
    plus mid-level mesocyclone history, for the M1/M2 diagnostics."""

    def __init__(self):
        self.t = []
        self.near_surface_zeta_max = []
        self.midlevel_meso = []
        self.updraft_helicity_max = []
        self.w_max = []

    def update(self, t, state, grid):
        rep = rotation_report(state, grid)
        self.t.append(float(t))
        self.near_surface_zeta_max.append(rep["near_surface_zeta_max"])
        self.midlevel_meso.append(rep["midlevel_mesocyclone"])
        self.updraft_helicity_max.append(rep["updraft_helicity_max"])
        self.w_max.append(rep["w_max"])
        return rep

    def peak(self) -> dict:
        f = lambda a: (max(a) if a else 0.0)
        return {
            "peak_near_surface_zeta": f(self.near_surface_zeta_max),
            "peak_midlevel_mesocyclone": f(self.midlevel_meso),
            "peak_updraft_helicity": f(self.updraft_helicity_max),
            "peak_w": f(self.w_max),
        }


__all__ = [
    "vorticity_3d", "vertical_vorticity", "updraft_helicity",
    "rotation_report", "TornadogenesisTracker",
]
