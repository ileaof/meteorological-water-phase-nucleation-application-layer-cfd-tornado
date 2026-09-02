"""LES subgrid-scale closure (item 3) -- replaces the demonstration Rayleigh drag.

The ``meteorological_flow`` storm bounds its otherwise-unbounded Boussinesq plume
with a linear Rayleigh drag (``gamma_damp``) plus a hard 120 m/s velocity clip.
Both are unphysical for a vortex: the Rayleigh drag damps the very rotation we
want, and the clip caps peak winds artificially.  Here they are replaced by a
**Smagorinsky** eddy-viscosity closure, so the subgrid dissipation is set by the
resolved strain rate:

    K_m = (C_s Delta)^2 |S| * f_stab(Ri) ,   |S| = sqrt(2 S_ij S_ij) ,
    Delta = (dx dy dz)^(1/3) ,

with a stability correction ``f_stab = sqrt(max(0, 1 - Ri/Ri_c))`` that damps the
mixing in statically stable layers (Lilly 1962).  Momentum feels
``d u_i/dt = d/dx_j ( K_m (du_i/dx_j) )`` (gradient-diffusion form) and scalars
feel ``K_h = K_m / Pr_t``.  A small background viscosity ``nu_background`` keeps
the operator well behaved where the resolved strain vanishes.

The only remaining velocity bound is an *extreme numerical guard*
(``StormDynamicsConfig.v_guard``, applied in the core), documented as a safety
rail that should never bite in a resolved simulation -- not a physical cap.
"""
from __future__ import annotations

import numpy as np

from meteorological_flow import thermodynamics as th
from meteorological_flow.grid import Grid
from meteorological_flow.state import FlowState

from .config import LESConfig


def _centered_velocity(state, grid):
    u, v, w = state.u, state.v, state.w
    uc = 0.5 * (u[:-1, :, :] + u[1:, :, :])
    vc = 0.5 * (v[:, :-1, :] + v[:, 1:, :])
    wc = 0.5 * (w[:, :, :-1] + w[:, :, 1:])
    return uc, vc, wc


def strain_and_viscosity(state: FlowState, grid: Grid, les: LESConfig,
                         theta0=None):
    """Return the eddy viscosity ``K_m`` (cell centres, (nx,ny,nz)).

    ``les.model``: ``"smagorinsky"`` (implemented), ``"none"`` (constant
    background viscosity only), or ``"tke15"`` (prognostic Deardorff TKE-1.5 --
    implemented; called statelessly here it returns the equilibrium ``K_m``, while the
    time loop evolves the TKE field via :func:`deardorff_tke_step`).
    """
    xp = grid.xp
    if les.model == "none":
        return xp.full(grid.center_shape, les.nu_background, dtype=float)
    if les.model == "tke15":
        # prognostic Deardorff closure; called statelessly here -> return the equilibrium K_m
        # (the time loop uses deardorff_tke_step to evolve the TKE field, see core._predictor).
        return deardorff_tke_step(state, grid, les, None, 0.0, theta0)[0]
    if les.model != "smagorinsky":
        raise ValueError("unknown LES model %r" % les.model)
    uc, vc, wc = _centered_velocity(state, grid)
    # all nine gradients via the grid's central-difference operators (periodic-aware
    # in x,y through the Laplacian/central stencils; one-sided at z walls)
    dudx = grid._central_x(uc); dudy = grid._central_y(uc); dudz = grid._central_z(uc)
    dvdx = grid._central_x(vc); dvdy = grid._central_y(vc); dvdz = grid._central_z(vc)
    dwdx = grid._central_x(wc); dwdy = grid._central_y(wc); dwdz = grid._central_z(wc)
    Sxx, Syy, Szz = dudx, dvdy, dwdz
    Sxy = 0.5 * (dudy + dvdx)
    Sxz = 0.5 * (dudz + dwdx)
    Syz = 0.5 * (dvdz + dwdy)
    modS = xp.sqrt(2.0 * (Sxx**2 + Syy**2 + Szz**2)
                   + 4.0 * (Sxy**2 + Sxz**2 + Syz**2)) + 1e-12
    dz = grid.dz if not getattr(grid, "stretched", False) else float(grid.dz_c.mean())
    Delta = (grid.dx * grid.dy * dz) ** (1.0 / 3.0)
    Km = (les.C_s * Delta) ** 2 * modS
    # stability (Richardson) correction: damp mixing in stable stratification
    if theta0 is not None and les.Ri_c > 0:
        thref = theta0
        dthdz = grid._central_z(state.theta - thref) + grid._central_z(thref)
        N2 = th.g0 / xp.clip(thref, 1e-3, None) * dthdz
        Ri = N2 / (modS ** 2)
        fstab = xp.sqrt(xp.clip(1.0 - xp.clip(Ri, 0.0, None) / les.Ri_c, 0.0, 1.0))
        Km = Km * fstab
    return Km + les.nu_background


def _div_k_grad(f, K, grid, periodic):
    """div( K grad f ) for a cell-centred field with cell-centred K.

    Face diffusivities are the arithmetic mean of the neighbouring cells; periodic
    wrap in x,y (storm) and zero-flux (Neumann) at the z walls.
    """
    xp = grid.xp
    dz = grid.dz if not getattr(grid, "stretched", False) else grid.dz_c[None, None, :]
    out = xp.zeros_like(f)
    # x
    if periodic:
        Kx = 0.5 * (K + xp.roll(K, -1, axis=0))            # face i+1/2
        Fx = Kx * (xp.roll(f, -1, axis=0) - f) / grid.dx
        out += (Fx - xp.roll(Fx, 1, axis=0)) / grid.dx
    else:
        Kx = 0.5 * (K[:-1, :, :] + K[1:, :, :])
        Fx = Kx * (f[1:, :, :] - f[:-1, :, :]) / grid.dx
        out[1:-1, :, :] += (Fx[1:, :, :] - Fx[:-1, :, :]) / grid.dx
    # y
    if periodic:
        Ky = 0.5 * (K + xp.roll(K, -1, axis=1))
        Fy = Ky * (xp.roll(f, -1, axis=1) - f) / grid.dy
        out += (Fy - xp.roll(Fy, 1, axis=1)) / grid.dy
    else:
        Ky = 0.5 * (K[:, :-1, :] + K[:, 1:, :])
        Fy = Ky * (f[:, 1:, :] - f[:, :-1, :]) / grid.dy
        out[:, 1:-1, :] += (Fy[:, 1:, :] - Fy[:, :-1, :]) / grid.dy
    # z (walls: zero flux)
    Kz = 0.5 * (K[:, :, :-1] + K[:, :, 1:])
    dzf = dz if not getattr(grid, "stretched", False) else grid.dzc_f[None, None, 1:-1]
    Fz = Kz * (f[:, :, 1:] - f[:, :, :-1]) / (dzf if getattr(grid, "stretched", False) else grid.dz)
    out[:, :, 1:-1] += (Fz[:, :, 1:] - Fz[:, :, :-1]) / (dz[:, :, 1:-1] if getattr(grid, "stretched", False) else grid.dz)
    return out


def _center_to_faces(tend_c, grid, axis):
    xp = grid.xp
    if axis == 0:
        out = xp.zeros(grid.u_shape)
        out[1:-1, :, :] = 0.5 * (tend_c[:-1, :, :] + tend_c[1:, :, :])
        if getattr(grid, "periodic", False):
            wrap = 0.5 * (tend_c[-1, :, :] + tend_c[0, :, :])
            out[0, :, :] = wrap; out[-1, :, :] = wrap
        else:
            out[0, :, :] = tend_c[0, :, :]; out[-1, :, :] = tend_c[-1, :, :]
    elif axis == 1:
        out = xp.zeros(grid.v_shape)
        out[:, 1:-1, :] = 0.5 * (tend_c[:, :-1, :] + tend_c[:, 1:, :])
        if getattr(grid, "periodic", False):
            wrap = 0.5 * (tend_c[:, -1, :] + tend_c[:, 0, :])
            out[:, 0, :] = wrap; out[:, -1, :] = wrap
        else:
            out[:, 0, :] = tend_c[:, 0, :]; out[:, -1, :] = tend_c[:, -1, :]
    else:
        out = xp.zeros(grid.w_shape)
        out[:, :, 1:-1] = 0.5 * (tend_c[:, :, :-1] + tend_c[:, :, 1:])
        out[:, :, 0] = 0.0; out[:, :, -1] = tend_c[:, :, -1]
    return out


def apply_les_momentum(state: FlowState, grid: Grid, Km, dt: float) -> None:
    """Apply the SGS momentum diffusion ``d u_i/dt = div(K_m grad u_i)`` in place."""
    xp = grid.xp
    periodic = getattr(grid, "periodic", False)
    uc, vc, wc = _centered_velocity(state, grid)
    du = _div_k_grad(uc, Km, grid, periodic)
    dv = _div_k_grad(vc, Km, grid, periodic)
    dw = _div_k_grad(wc, Km, grid, periodic)
    state.u += dt * _center_to_faces(du, grid, 0)
    state.v += dt * _center_to_faces(dv, grid, 1)
    state.w += dt * _center_to_faces(dw, grid, 2)


# ---------------------------------------------------------------------------
# TKE-1.5 (Deardorff 1980) prognostic subgrid closure -- item 3c.
# A prognostic subgrid TKE ``e`` sets the eddy viscosity, instead of the diagnostic
# Smagorinsky ``K_m = (C_s Delta)^2 |S|``:
#     K_m = C_k l sqrt(e),   K_h = (1 + 2 l/Delta) K_m,
#     de/dt = P_shear + P_buoy - eps + div(2 K_m grad e),
#     P_shear = K_m |S|^2,  P_buoy = -K_h N^2,  eps = C_eps e^{3/2}/l,
#     l = Delta (reduced to 0.76 sqrt(e)/N in stable stratification),
#     C_eps = 0.19 + 0.51 l/Delta.
# Reference: Deardorff (1980, BLM 18, 495); Moeng & Wyngaard (1988).
# ---------------------------------------------------------------------------
_CK = 0.10
TKE_MIN = 1e-6


def _strain_and_stability(state, grid, theta0):
    """Return ``(modS2, N2, Delta)`` at cell centres: squared strain, Brunt-Vaisala N^2, and
    the filter width -- the inputs both the Smagorinsky and Deardorff closures share."""
    xp = grid.xp
    uc, vc, wc = _centered_velocity(state, grid)
    dudx = grid._central_x(uc); dudy = grid._central_y(uc); dudz = grid._central_z(uc)
    dvdx = grid._central_x(vc); dvdy = grid._central_y(vc); dvdz = grid._central_z(vc)
    dwdx = grid._central_x(wc); dwdy = grid._central_y(wc); dwdz = grid._central_z(wc)
    Sxy = 0.5 * (dudy + dvdx); Sxz = 0.5 * (dudz + dwdx); Syz = 0.5 * (dvdz + dwdy)
    modS2 = (2.0 * (dudx ** 2 + dvdy ** 2 + dwdz ** 2) + 4.0 * (Sxy ** 2 + Sxz ** 2 + Syz ** 2))
    dz = grid.dz if not getattr(grid, "stretched", False) else float(grid.dz_c.mean())
    Delta = (grid.dx * grid.dy * dz) ** (1.0 / 3.0)
    dthdz = grid._central_z(state.theta - theta0) + grid._central_z(theta0) if theta0 is not None \
        else grid._central_z(state.theta)
    thref = theta0 if theta0 is not None else state.theta
    N2 = th.g0 / xp.clip(thref, 1e-3, None) * dthdz
    return modS2, N2, Delta


def deardorff_viscosity(tke, N2, Delta, les, xp=np):
    """Deardorff eddy viscosity/diffusivity from subgrid TKE: returns ``(K_m, K_h, l, eps)``."""
    e = xp.clip(tke, TKE_MIN, None)
    Ns = xp.sqrt(xp.clip(N2, 1e-12, None))
    l = xp.where(N2 > 1e-12, xp.minimum(Delta, 0.76 * xp.sqrt(e) / Ns), Delta)
    Km = _CK * l * xp.sqrt(e)
    Kh = (1.0 + 2.0 * l / Delta) * Km
    Ceps = 0.19 + 0.51 * l / Delta
    eps = Ceps * e ** 1.5 / xp.clip(l, 1e-3, None)
    return Km, Kh, l, eps


def deardorff_tke_step(state, grid, les, tke, dt, theta0):
    """Advance the prognostic subgrid TKE one step and return ``(K_m, tke_new)``.

    ``tke=None`` initialises ``e`` from the local Smagorinsky-equilibrium (so the first call is
    well posed); ``dt=0`` returns the equilibrium viscosity without stepping (diagnostic use)."""
    xp = grid.xp
    modS2, N2, Delta = _strain_and_stability(state, grid, theta0)
    if tke is None:
        Km_s = (les.C_s * Delta) ** 2 * xp.sqrt(modS2 + 1e-12)
        tke = xp.clip((Km_s / (_CK * Delta)) ** 2, TKE_MIN, None)
    Km, Kh, l, eps = deardorff_viscosity(tke, N2, Delta, les, xp=xp)
    if dt > 0.0:
        periodic = getattr(grid, "periodic", False)
        Ps = Km * modS2                                    # shear production (>= 0)
        Pb = -Kh * N2                                      # buoyancy: <0 stable, >0 unstable
        diff = _div_k_grad(tke, 2.0 * Km, grid, periodic)  # TKE self-diffusion
        tke = xp.clip(tke + dt * (Ps + Pb - eps + diff), TKE_MIN, None)
    return Km + les.nu_background, tke


def les_scalar_diffusion(field, Km, grid: Grid, les: LESConfig, dt: float,
                         base=None):
    """Return ``field + dt * div(K_h grad field')`` with ``K_h = K_m / Pr_t``.

    ``base`` (e.g. theta0(z)/qv0(z)) diffuses only the perturbation, exactly as
    :func:`meteorological_flow.diffusion.diffuse_center` does, so the stratified
    reference stays a discrete equilibrium.
    """
    Kh = Km / max(les.Pr_t, 1e-6)
    f = field if base is None else (field - base)
    return field + dt * _div_k_grad(f, Kh, grid, getattr(grid, "periodic", False))


__all__ = [
    "strain_and_viscosity", "apply_les_momentum", "les_scalar_diffusion",
    "deardorff_viscosity", "deardorff_tke_step",
]
