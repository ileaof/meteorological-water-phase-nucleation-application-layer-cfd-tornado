"""Environmental soundings for the rotating storm -- curved-hodograph shear.

The thermodynamic column (theta0, qv0, p0, T0, rho0) is reused *unchanged* from
:func:`meteorological_flow.base_state.weisman_klemp` (the standard CAPE ~2000-2500
J/kg idealised sounding).  Only the winds are replaced: the demonstration
:func:`weisman_klemp` imposes a purely *unidirectional* ``u0(z)`` (``v0 = 0``),
which carries no streamwise (storm-relative) vorticity, so it can split a storm
but cannot preferentially spin up low-level rotation.

Here we populate BOTH ``u0(z)`` and ``v0(z)`` to give a **curved** hodograph:

* ``"unidirectional"`` -- straight hodograph (``v0 = 0``); the Klemp-Wilhelmson /
  Weisman-Klemp storm-splitting reference (M1).
* ``"quarter_circle"`` -- the wind vector turns through a quarter circle over the
  lowest ``z_turn`` (veering with height), then is constant aloft.  The curvature
  is the source of storm-relative helicity that selects the right-moving
  supercell and feeds low-level rotation (Rotunno & Klemp 1982, 1985) -- the M2
  target.

References: Weisman & Klemp (1982); Rotunno & Klemp (1982, 1985); Davies-Jones
(1984, storm-relative helicity); Bunkers et al. (2000, storm-motion estimate).
"""
from __future__ import annotations

import numpy as np

from meteorological_flow.base_state import BaseState, weisman_klemp

from .config import HodographConfig


def _quarter_circle_winds(z: np.ndarray, hodo: HodographConfig):
    """Quarter-circle hodograph winds (u0, v0) [m/s] on heights ``z``.

    The wind traces a quarter circle of radius ``U_max`` from the origin at the
    surface to ``(U_max, -U_max)`` at ``z_turn``.  The wind vector veers (turns
    clockwise) with height, giving strong **positive** storm-relative helicity --
    the Northern-Hemisphere right-moving supercell configuration -- then holds
    constant above ``z_turn``.
    """
    R = hodo.U_max
    zt = max(hodo.z_turn, 1.0)
    beta = 0.5 * np.pi * np.clip(z / zt, 0.0, 1.0)
    u0 = R * np.sin(beta)
    v0 = -R * (1.0 - np.cos(beta))
    # above the turning layer the winds are constant (u,v held at their z_turn value)
    aloft = z > zt
    u0[aloft] = R * np.sin(0.5 * np.pi)
    v0[aloft] = -R * (1.0 - np.cos(0.5 * np.pi))
    return u0, v0


def _unidirectional_winds(z: np.ndarray, hodo: HodographConfig):
    """Straight (unidirectional) hodograph: u0 = U_max tanh(z/u_half), v0 = 0."""
    u0 = hodo.U_max * np.tanh(z / max(hodo.u_half, 1.0))
    v0 = np.zeros_like(z)
    return u0, v0


def build_sounding(grid, hodo: HodographConfig | None = None, **wk_kwargs) -> BaseState:
    """Return a :class:`BaseState` with the WK thermodynamics + a curved hodograph.

    ``wk_kwargs`` are forwarded to :func:`weisman_klemp` (e.g. ``qv_sfc``,
    ``theta_sfc``) so CAPE can be tuned; the wind arguments of ``weisman_klemp``
    are overridden here.
    """
    hodo = hodo or HodographConfig()
    base = weisman_klemp(grid, u_shear=0.0, **wk_kwargs)   # thermodynamics only
    z = np.asarray(base.zc, dtype=float)
    if hodo.kind == "unidirectional":
        u0, v0 = _unidirectional_winds(z, hodo)
    elif hodo.kind == "quarter_circle":
        u0, v0 = _quarter_circle_winds(z, hodo)
    else:
        raise ValueError("unknown hodograph kind %r" % hodo.kind)
    base.u0 = u0
    base.v0 = v0
    return base


# ---------------------------------------------------------------------------
# hodograph diagnostics (bulk shear + storm-relative helicity)
# ---------------------------------------------------------------------------
def bulk_shear(base: BaseState, z1: float = 0.0, z2: float = 6000.0) -> float:
    """|V(z2) - V(z1)| bulk vertical wind shear [m/s]."""
    z = np.asarray(base.zc)
    u = np.interp([z1, z2], z, np.asarray(base.u0))
    v = np.interp([z1, z2], z, np.asarray(base.v0))
    return float(np.hypot(u[1] - u[0], v[1] - v[0]))


def bunkers_storm_motion(base: BaseState, deviation: float = 7.5):
    """Right-moving supercell motion estimate (Bunkers et al. 2000, simplified).

    Mean 0-6 km wind plus a ``deviation`` [m/s] offset to the RIGHT of the 0-6 km
    shear vector.  Returns ``(cx, cy)`` [m/s].
    """
    z = np.asarray(base.zc)
    sel = z <= 6000.0
    if sel.sum() < 2:
        sel = slice(None)
    um = float(np.mean(np.asarray(base.u0)[sel]))
    vm = float(np.mean(np.asarray(base.v0)[sel]))
    # 0-6 km shear vector
    u = np.interp([0.0, 6000.0], z, np.asarray(base.u0))
    v = np.interp([0.0, 6000.0], z, np.asarray(base.v0))
    shx, shy = u[1] - u[0], v[1] - v[0]
    smag = np.hypot(shx, shy) + 1e-9
    # unit vector 90 deg to the RIGHT of the shear (clockwise): (shy, -shx)/|sh|
    cx = um + deviation * (shy / smag)
    cy = vm + deviation * (-shx / smag)
    return cx, cy


def storm_relative_helicity(base: BaseState, z_top: float = 3000.0,
                            storm_motion=None) -> float:
    """Storm-relative helicity SRH [m^2/s^2] over 0..``z_top`` (Davies-Jones 1984).

        SRH = - integral_0^h k . [ (V - C) x dV/dz ]  dz
            =   integral_0^h [ (v-cy) du/dz - (u-cx) dv/dz ] dz

    ``storm_motion`` ``(cx, cy)``; None -> :func:`bunkers_storm_motion` (right-mover).
    Positive SRH (right-moving supercell in the NH) for a clockwise-curving
    (veering) hodograph.
    """
    z = np.asarray(base.zc, dtype=float)
    u = np.asarray(base.u0, dtype=float)
    v = np.asarray(base.v0, dtype=float)
    cx, cy = storm_motion if storm_motion is not None else bunkers_storm_motion(base)
    sel = z <= z_top
    if sel.sum() < 2:
        return 0.0
    zc = z[sel]; uc = u[sel]; vc = v[sel]
    dudz = np.gradient(uc, zc)
    dvdz = np.gradient(vc, zc)
    integrand = (vc - cy) * dudz - (uc - cx) * dvdz
    _trapz = getattr(np, "trapezoid", np.trapz)
    return float(_trapz(integrand, zc))


__all__ = [
    "build_sounding", "bulk_shear", "bunkers_storm_motion",
    "storm_relative_helicity",
]
