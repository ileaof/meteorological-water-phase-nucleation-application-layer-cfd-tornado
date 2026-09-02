"""Sustained low-level mesoscale-ascent forcing -- a dryline/convergence proxy.

A one-shot warm bubble in a *really* capped environment (a CIN inherited from real
ERA5 / radiosonde data) pulses once and decays: nothing renews the lifting through the
cap, so the incipient updraft detrains and dies (exactly what the real-data Moore runs
showed -- w_max 9 -> 1 over 28 min).  The real supercell broke the cap by **sustained
mesoscale ascent** -- dryline convergence plus afternoon boundary-layer heating -- not a
single thermal.

This module represents that with a smooth, low-level heating (+ optional moistening)
cylinder held for a fixed ``duration_s`` so parcels are *continuously* lifted through
the cap and a supercell can establish; afterwards the forcing is removed and the storm
must sustain on its own dynamics.  It is **additive and opt-in**
(``StormDynamicsConfig.forcing.enabled`` -- off by default, so idealised runs keep the
one-shot bubble and every existing test is untouched)."""
from __future__ import annotations

import numpy as np


def meso_forcing_mask(grid, cfg_forcing, center=None, xp=np):
    """Smooth [0, 1] cylinder: ~1 in the low-level core, cosine-tapering to 0 at
    ``radius_m`` horizontally and at ``z_top_m`` in the vertical (peak at the surface).
    """
    fc = cfg_forcing
    cx, cy = center if center is not None else (0.5 * grid.Lx, 0.5 * grid.Ly)
    xc = xp.asarray(grid.xc); yc = xp.asarray(grid.yc); zc = xp.asarray(grid.zc)
    X = xc[:, None, None]; Y = yc[None, :, None]; Z = zc[None, None, :]
    r = xp.sqrt((X - cx) ** 2 + (Y - cy) ** 2)
    fr = 0.5 * (1.0 + xp.cos(xp.pi * xp.clip(r / fc.radius_m, 0.0, 1.0)))
    fz = xp.where(Z < fc.z_top_m,
                  0.5 * (1.0 + xp.cos(xp.pi * xp.clip(Z / fc.z_top_m, 0.0, 1.0))), 0.0)
    return fr * fz


def apply_meso_forcing(state, grid, cfg_forcing, t, dt, center=None, xp=np):
    """Add the sustained heating(+moistening) tendency while ``t < duration_s``.

    Returns ``True`` when the forcing was active this step (so callers can log it),
    ``False`` once it has switched off or is disabled."""
    fc = cfg_forcing
    if not getattr(fc, "enabled", False) or t >= fc.duration_s:
        return False
    if not (fc.heat_rate_K_s or fc.moist_rate_kgkg_s):
        return False
    m = meso_forcing_mask(grid, fc, center=center, xp=xp)
    if fc.heat_rate_K_s:
        state.theta = state.theta + dt * fc.heat_rate_K_s * m
    if fc.moist_rate_kgkg_s:
        state.qv = xp.maximum(state.qv + dt * fc.moist_rate_kgkg_s * m, 0.0)
    return True
