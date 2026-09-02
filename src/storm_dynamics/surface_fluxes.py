"""Surface sensible + latent heat fluxes -- the boundary-layer heat/moisture source.

Complements the momentum drag (:mod:`surface_drag`): a bulk-aerodynamic sensible-heat flux relaxes
the lowest model level's potential temperature toward a (warm) surface value, and a latent-heat
flux relaxes its humidity toward the surface saturation (moist ground), each at rate
``C |V| / dz0``.  This deepens the mixed layer, sustains the inflow's theta-e, and lets the cold
pool recover -- all of which modulate low-level vorticity.  Opt-in (``StormDynamicsConfig.fluxes``);
off by default so existing runs are byte-unchanged.  `grid.xp`-generic; implicit (unconditionally
stable) like the drag.
"""
from __future__ import annotations

from meteorological_flow import thermodynamics as th

CP_D = 1004.0
LV = 2.5e6


def _dz0(grid):
    return grid.dz if not getattr(grid, "stretched", False) else float(grid.dz_c[0])


def surface_targets(state, grid, sfc, base=None):
    """Resolve the surface potential-temperature and humidity targets (2-D or scalar)."""
    xp = grid.xp
    if sfc.theta_sfc_K is not None:
        th_sfc = float(sfc.theta_sfc_K)
    elif base is not None:
        th_sfc = float(xp.asarray(base.theta0)[0]) + sfc.dtheta_sfc_K
    else:
        th_sfc = float(xp.mean(state.theta[:, :, 0])) + sfc.dtheta_sfc_K
    if sfc.saturate_surface and state.T is not None and state.P_total is not None:
        qv_sfc = th.q_v_from_p_v(th.psat_water(state.T[:, :, 0], xp=xp), state.P_total[:, :, 0], xp=xp)
    elif sfc.qv_sfc_kgkg is not None:
        qv_sfc = float(sfc.qv_sfc_kgkg)
    else:
        qv_sfc = None
    return th_sfc, qv_sfc


def apply_surface_fluxes(state, grid, dt, sfc, base=None) -> None:
    """Implicitly relax the lowest level toward the surface targets, in place."""
    if not getattr(sfc, "enabled", False):
        return
    xp = grid.xp
    dz0 = _dz0(grid)
    uc0 = 0.5 * (state.u[:-1, :, 0] + state.u[1:, :, 0])
    vc0 = 0.5 * (state.v[:, :-1, 0] + state.v[:, 1:, 0])
    speed = xp.maximum(xp.sqrt(uc0 * uc0 + vc0 * vc0), sfc.U_min)
    th_sfc, qv_sfc = surface_targets(state, grid, sfc, base=base)
    if sfc.C_h > 0.0:
        rh = sfc.C_h * speed / dz0
        state.theta[:, :, 0] = (state.theta[:, :, 0] + rh * dt * th_sfc) / (1.0 + rh * dt)
    if qv_sfc is not None and sfc.C_q > 0.0:
        rq = sfc.C_q * speed / dz0
        state.qv[:, :, 0] = xp.maximum((state.qv[:, :, 0] + rq * dt * qv_sfc) / (1.0 + rq * dt), 0.0)


def surface_flux_report(state, grid, sfc, base=None) -> dict:
    """Diagnose the bulk surface fluxes H (sensible) and LE (latent) [W/m^2]."""
    xp = grid.xp
    if not getattr(sfc, "enabled", False):
        return {"sensible_heat_flux_W_m2": 0.0, "latent_heat_flux_W_m2": 0.0}
    uc0 = 0.5 * (state.u[:-1, :, 0] + state.u[1:, :, 0])
    vc0 = 0.5 * (state.v[:, :-1, 0] + state.v[:, 1:, 0])
    speed = xp.maximum(xp.sqrt(uc0 * uc0 + vc0 * vc0), sfc.U_min)
    th_sfc, qv_sfc = surface_targets(state, grid, sfc, base=base)
    rho0 = state.rho[:, :, 0] if state.rho is not None else 1.1
    H = rho0 * CP_D * sfc.C_h * speed * (th_sfc - state.theta[:, :, 0])
    out = {"sensible_heat_flux_W_m2": float(xp.mean(H))}
    if qv_sfc is not None:
        LE = rho0 * LV * sfc.C_q * speed * (qv_sfc - state.qv[:, :, 0])
        out["latent_heat_flux_W_m2"] = float(xp.mean(LE))
    else:
        out["latent_heat_flux_W_m2"] = 0.0
    return out


def neutral_drag_coefficient(z1_m, z0_m, kappa=0.4):
    """Neutral bulk transfer coefficient from the log-law: C = (kappa/ln(z1/z0))^2."""
    import math
    return (kappa / math.log(max(z1_m, 2 * z0_m) / z0_m)) ** 2


__all__ = ["surface_targets", "apply_surface_fluxes", "surface_flux_report",
           "neutral_drag_coefficient"]
