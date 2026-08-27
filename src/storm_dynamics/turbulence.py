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
    documented as future work, raises :class:`NotImplementedError`).
    """
    xp = grid.xp
    if les.model == "none":
        return xp.full(grid.center_shape, les.nu_background, dtype=float)
    if les.model == "tke15":
        raise NotImplementedError(
            "TKE-1.5 (Deardorff) closure is a documented future option; use "
            "les_model='smagorinsky' (see docs/storm_dynamics_guide.md).")
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
]
