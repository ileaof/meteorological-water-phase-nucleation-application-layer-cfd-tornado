"""Conservative flux-form momentum advection on the staggered Arakawa C-grid.

This is the *enabling* piece of the rotating core (item 1): the demonstration
``meteorological_flow`` solver leaves momentum advection off (a documented v1
simplification), so ``(u.grad)u`` -- and with it the **tilting** and **stretching**
of vorticity -- is absent and the flow cannot spin up a vortex.  Here we compute

    du_i/dt = - d/dx_j ( u_j u_i )        (flux / divergence form)

for each staggered component on its own control volume, using an upwind-biased
2nd-order MUSCL (minmod) reconstruction in the spirit of
:func:`meteorological_flow.advection.advect_center_massflux`.  Because the
transporting velocity is the projection's discretely divergence-free field, the
flux form is equivalent to the advective form ``u_j d u_i/dx_j`` while
telescoping to the boundary flux -- so the domain-integrated momentum of each
component is conserved (to the boundary flux) at any intensity, with no velocity
clip.

Grid / boundaries.  The storm runs with **periodic** lateral (x, y) boundaries
(the mean-wind shear is ingested through them), so horizontal neighbours wrap
(``numpy.roll``); the vertical is bounded (``w = 0`` at the ground, outflow/damping
at the top), so the z-fluxes use the physical wall value (no wrap).  A
non-periodic lateral configuration falls back to a zero-gradient (clamped) edge,
which carries zero advective flux through a closed wall.
"""
from __future__ import annotations

import numpy as np

from meteorological_flow.advection import _minmod
from meteorological_flow.grid import Grid
from meteorological_flow.state import FlowState


def _slope(q, axis, xp, periodic):
    """minmod-limited centred slope of ``q`` along ``axis`` (0 at hard edges)."""
    if periodic:
        fwd = xp.roll(q, -1, axis=axis) - q
        bwd = q - xp.roll(q, 1, axis=axis)
        return _minmod(bwd, fwd, xp=xp)
    s = xp.zeros_like(q)
    sl = [slice(None)] * q.ndim
    slm = [slice(None)] * q.ndim
    slp = [slice(None)] * q.ndim
    sl[axis] = slice(1, -1); slm[axis] = slice(0, -2); slp[axis] = slice(2, None)
    fwd = q[tuple(slp)] - q[tuple(sl)]
    bwd = q[tuple(sl)] - q[tuple(slm)]
    s[tuple(sl)] = _minmod(bwd, fwd, xp=xp)
    return s


def _face_upwind(qm, qp, vel, sm, sp, xp, order):
    """Reconstruct the advected value at a flux face between minus/plus states."""
    if order >= 2:
        left = qm + 0.5 * sm
        right = qp - 0.5 * sp
        return xp.where(vel > 0.0, left, right)
    return xp.where(vel > 0.0, qm, qp)


def _u_tendency(st, g, order, periodic, fluxes=None):
    """-div(u u) on the u-faces (nx+1, ny, nz).  If ``fluxes`` is a dict, the cell-centred
    x-flux ``Fx_u`` (the normal u-momentum flux ``Uc·u``, at (nx,ny,nz)) is recorded into it
    for the interface reflux — the same array the divergence below uses (single source)."""
    xp = g.xp
    u, v, w = st.u, st.v, st.w
    dz = g.dz if not getattr(g, "stretched", False) else g.dz_c[None, None, :]
    tend = xp.zeros_like(u)

    # --- x-flux at cell centres (nx,ny,nz): F = Uc * u_recon, Uc = 0.5(u[i]+u[i+1])
    Uc = 0.5 * (u[:-1, :, :] + u[1:, :, :])                 # (nx,ny,nz)
    su = _slope(u, 0, xp, periodic)                         # slope of u along x
    Fx = Uc * _face_upwind(u[:-1, :, :], u[1:, :, :], Uc,
                           su[:-1, :, :], su[1:, :, :], xp, order)   # (nx,ny,nz)
    if fluxes is not None:
        fluxes["Fx_u"] = Fx
    # tendency at face I = -(Fx[I]-Fx[I-1])/dx.  Interior faces 1..nx-1:
    tend[1:-1, :, :] += -(Fx[1:, :, :] - Fx[:-1, :, :]) / g.dx
    if periodic:
        wrap = -(Fx[0, :, :] - Fx[-1, :, :]) / g.dx
        tend[0, :, :] += wrap
        tend[-1, :, :] += wrap

    # --- y-flux at (x-face, y-face) corners: transporting v averaged onto x-faces
    # Vx[I,J] = 0.5(v[I-1,J]+v[I,J]); advected u from u[:,J-1],u[:,J] upwind by Vx.
    if periodic:
        v_xm = xp.roll(v, 1, axis=0)                        # v[i-1] in x (wrap)
    else:
        v_xm = xp.zeros_like(v); v_xm[1:, :, :] = v[:-1, :, :]
    Vx = 0.5 * (v_xm + v)                                   # (nx,ny+1,nz) at x-faces i=0..nx-1
    # bring Vx to the full (nx+1) x-face set: face nx == face 0 (periodic) else 0
    Vxf = xp.zeros((g.nx + 1, g.ny + 1, g.nz))
    Vxf[:-1, :, :] = Vx
    Vxf[-1, :, :] = Vx[0, :, :] if periodic else 0.0
    suy = _slope(u, 1, xp, periodic)                        # slope of u along y
    if periodic:
        u_ym = xp.roll(u, 1, axis=1); s_ym = xp.roll(suy, 1, axis=1)
        # corner value at J=0..ny-1 (between cells J-1,J); Vxf interior J=1..ny-1 + wrap
        Vc = Vxf[:, :-1, :]                                 # (nx+1,ny,nz) corners J=0..ny-1
        u_corner = _face_upwind(u_ym, u, Vc, s_ym, suy, xp, order)   # (nx+1,ny,nz)
        Fy = Vc * u_corner                                 # corner flux at J=0..ny-1
        # tendency at u-cell j: -(Fy[j+1]-Fy[j])/dy with periodic wrap in J
        tend += -(xp.roll(Fy, -1, axis=1) - Fy) / g.dy
    else:
        # clamped: corners at interior J only; edge cells get zero-gradient (no flux)
        u_ym = xp.zeros_like(u); u_ym[:, 1:, :] = u[:, :-1, :]
        s_ym = xp.zeros_like(suy); s_ym[:, 1:, :] = suy[:, :-1, :]
        Vc = Vxf[:, :-1, :]
        u_corner = _face_upwind(u_ym, u, Vc, s_ym, suy, xp, order)
        Fy_full = xp.zeros((g.nx + 1, g.ny + 1, g.nz))
        Fy_full[:, :-1, :] = Vc * u_corner
        tend[:, 1:-1, :] += -(Fy_full[:, 2:-1, :] - Fy_full[:, 1:-2, :]) / g.dy

    # --- z-flux at (x-face, z-face) corners: transporting w averaged onto x-faces
    if periodic:
        w_xm = xp.roll(w, 1, axis=0)
    else:
        w_xm = xp.zeros_like(w); w_xm[1:, :, :] = w[:-1, :, :]
    Wx = 0.5 * (w_xm + w)                                   # (nx,ny,nz+1) at x-faces i=0..nx-1
    Wxf = xp.zeros((g.nx + 1, g.ny, g.nz + 1))
    Wxf[:-1, :, :] = Wx
    Wxf[-1, :, :] = Wx[0, :, :] if periodic else 0.0
    suz = _slope(u, 2, xp, periodic=False)                 # z never periodic
    u_zm = xp.zeros_like(u); u_zm[:, :, 1:] = u[:, :, :-1]
    s_zm = xp.zeros_like(suz); s_zm[:, :, 1:] = suz[:, :, :-1]
    Wc = Wxf[:, :, 1:-1]                                    # interior z-faces K=1..nz-1
    u_zcorner = _face_upwind(u_zm[:, :, 1:], u[:, :, 1:], Wc,
                             s_zm[:, :, 1:], suz[:, :, 1:], xp, order)
    Fz = xp.zeros((g.nx + 1, g.ny, g.nz + 1))
    Fz[:, :, 1:-1] = Wc * u_zcorner                        # zero flux through z-walls
    tend += -(Fz[:, :, 1:] - Fz[:, :, :-1]) / dz
    return tend


def _v_tendency(st, g, order, periodic, fluxes=None):
    """-div(u v) on the v-faces (nx, ny+1, nz).  Mirror of :func:`_u_tendency`; records the
    cell-centred normal v-momentum flux ``Fy_v`` when ``fluxes`` is given."""
    xp = g.xp
    u, v, w = st.u, st.v, st.w
    dz = g.dz if not getattr(g, "stretched", False) else g.dz_c[None, None, :]
    tend = xp.zeros_like(v)

    # y-flux at cell centres: Vc = 0.5(v[j]+v[j+1]); advect v along y
    Vc = 0.5 * (v[:, :-1, :] + v[:, 1:, :])                 # (nx,ny,nz)
    sv = _slope(v, 1, xp, periodic)
    Fy = Vc * _face_upwind(v[:, :-1, :], v[:, 1:, :], Vc,
                           sv[:, :-1, :], sv[:, 1:, :], xp, order)
    if fluxes is not None:
        fluxes["Fy_v"] = Fy
    tend[:, 1:-1, :] += -(Fy[:, 1:, :] - Fy[:, :-1, :]) / g.dy
    if periodic:
        wrap = -(Fy[:, 0, :] - Fy[:, -1, :]) / g.dy
        tend[:, 0, :] += wrap
        tend[:, -1, :] += wrap

    # x-flux at corners: transporting u averaged onto y-faces
    if periodic:
        u_ym = xp.roll(u, 1, axis=1)
    else:
        u_ym = xp.zeros_like(u); u_ym[:, 1:, :] = u[:, :-1, :]
    Uy = 0.5 * (u_ym + u)                                   # (nx+1,ny,nz) at y-faces j=0..ny-1
    Uyf = xp.zeros((g.nx + 1, g.ny + 1, g.nz))
    Uyf[:, :-1, :] = Uy
    Uyf[:, -1, :] = Uy[:, 0, :] if periodic else 0.0
    svx = _slope(v, 0, xp, periodic)
    if periodic:
        v_xm = xp.roll(v, 1, axis=0); s_xm = xp.roll(svx, 1, axis=0)
        Uc2 = Uyf[:-1, :, :]                               # corners I=0..nx-1
        v_corner = _face_upwind(v_xm, v, Uc2, s_xm, svx, xp, order)
        Fx = Uc2 * v_corner
        tend += -(xp.roll(Fx, -1, axis=0) - Fx) / g.dx
    else:
        v_xm = xp.zeros_like(v); v_xm[1:, :, :] = v[:-1, :, :]
        s_xm = xp.zeros_like(svx); s_xm[1:, :, :] = svx[:-1, :, :]
        Uc2 = Uyf[:-1, :, :]
        v_corner = _face_upwind(v_xm, v, Uc2, s_xm, svx, xp, order)
        Fx_full = xp.zeros((g.nx + 1, g.ny + 1, g.nz))
        Fx_full[:-1, :, :] = Uc2 * v_corner
        tend[1:-1, :, :] += -(Fx_full[2:-1, :, :] - Fx_full[1:-2, :, :]) / g.dx

    # z-flux at corners: transporting w averaged onto y-faces
    if periodic:
        w_ym = xp.roll(w, 1, axis=1)
    else:
        w_ym = xp.zeros_like(w); w_ym[:, 1:, :] = w[:, :-1, :]
    Wy = 0.5 * (w_ym + w)                                   # (nx,ny,nz+1) at y-faces j=0..ny-1
    Wyf = xp.zeros((g.nx, g.ny + 1, g.nz + 1))
    Wyf[:, :-1, :] = Wy
    Wyf[:, -1, :] = Wy[:, 0, :] if periodic else 0.0
    svz = _slope(v, 2, xp, periodic=False)
    v_zm = xp.zeros_like(v); v_zm[:, :, 1:] = v[:, :, :-1]
    s_zm = xp.zeros_like(svz); s_zm[:, :, 1:] = svz[:, :, :-1]
    Wc = Wyf[:, :, 1:-1]
    v_zcorner = _face_upwind(v_zm[:, :, 1:], v[:, :, 1:], Wc,
                             s_zm[:, :, 1:], svz[:, :, 1:], xp, order)
    Fz = xp.zeros((g.nx, g.ny + 1, g.nz + 1))
    Fz[:, :, 1:-1] = Wc * v_zcorner
    tend += -(Fz[:, :, 1:] - Fz[:, :, :-1]) / dz
    return tend


def _w_tendency(st, g, order, periodic):
    """-div(u w) on the w-faces (nx, ny, nz+1)."""
    xp = g.xp
    u, v, w = st.u, st.v, st.w
    dz = g.dz if not getattr(g, "stretched", False) else g.dz_c[None, None, :]
    dzc_f = g.dzc_f[None, None, :] if getattr(g, "stretched", False) else g.dz
    tend = xp.zeros_like(w)

    # z-flux at cell centres: Wc = 0.5(w[k]+w[k+1]); advect w along z
    Wc = 0.5 * (w[:, :, :-1] + w[:, :, 1:])                 # (nx,ny,nz)
    sw = _slope(w, 2, xp, periodic=False)
    Fz = Wc * _face_upwind(w[:, :, :-1], w[:, :, 1:], Wc,
                           sw[:, :, :-1], sw[:, :, 1:], xp, order)
    # interior w-faces K=1..nz-1 (K=0 ground, K=nz top set by BC).  The w-face
    # control volume spans two half-cells, so divide by the centre-to-centre
    # spacing dzc_f (== dz for a uniform grid).
    _zspacing = dzc_f[:, :, 1:-1] if getattr(g, "stretched", False) else g.dz
    tend[:, :, 1:-1] += -(Fz[:, :, 1:] - Fz[:, :, :-1]) / _zspacing

    # x-flux at corners: transporting u interpolated onto w's z-faces
    # (average the two z-cells around each interior z-face; clamp at the walls)
    Uzf = xp.zeros((g.nx + 1, g.ny, g.nz + 1))
    Uzf[:, :, 1:-1] = 0.5 * (u[:, :, :-1] + u[:, :, 1:])
    Uzf[:, :, 0] = u[:, :, 0]; Uzf[:, :, -1] = u[:, :, -1]
    swx = _slope(w, 0, xp, periodic)
    if periodic:
        w_xm = xp.roll(w, 1, axis=0); s_xm = xp.roll(swx, 1, axis=0)
        Uc2 = Uzf[:-1, :, :]                               # corners I=0..nx-1
        w_corner = _face_upwind(w_xm, w, Uc2, s_xm, swx, xp, order)
        Fx = Uc2 * w_corner
        tend += -(xp.roll(Fx, -1, axis=0) - Fx) / g.dx
    else:
        w_xm = xp.zeros_like(w); w_xm[1:, :, :] = w[:-1, :, :]
        s_xm = xp.zeros_like(swx); s_xm[1:, :, :] = swx[:-1, :, :]
        Uc2 = Uzf[:-1, :, :]
        w_corner = _face_upwind(w_xm, w, Uc2, s_xm, swx, xp, order)
        Fx_full = xp.zeros((g.nx + 1, g.ny, g.nz + 1))
        Fx_full[:-1, :, :] = Uc2 * w_corner
        tend[1:-1, :, :] += -(Fx_full[2:-1, :, :] - Fx_full[1:-2, :, :]) / g.dx

    # y-flux at corners: transporting v averaged onto z-faces
    Vzf = xp.zeros((g.nx, g.ny + 1, g.nz + 1))
    Vzf[:, :, 1:-1] = 0.5 * (v[:, :, :-1] + v[:, :, 1:])
    Vzf[:, :, 0] = v[:, :, 0]; Vzf[:, :, -1] = v[:, :, -1]
    swy = _slope(w, 1, xp, periodic)
    if periodic:
        w_ym = xp.roll(w, 1, axis=1); s_ym = xp.roll(swy, 1, axis=1)
        Vc2 = Vzf[:, :-1, :]                               # corners J=0..ny-1
        w_corner = _face_upwind(w_ym, w, Vc2, s_ym, swy, xp, order)
        Fy = Vc2 * w_corner
        tend += -(xp.roll(Fy, -1, axis=1) - Fy) / g.dy
    else:
        w_ym = xp.zeros_like(w); w_ym[:, 1:, :] = w[:, :-1, :]
        s_ym = xp.zeros_like(swy); s_ym[:, 1:, :] = swy[:, :-1, :]
        Vc2 = Vzf[:, :-1, :]
        w_corner = _face_upwind(w_ym, w, Vc2, s_ym, swy, xp, order)
        Fy_full = xp.zeros((g.nx, g.ny + 1, g.nz + 1))
        Fy_full[:, :-1, :] = Vc2 * w_corner
        tend[:, 1:-1, :] += -(Fy_full[:, 2:-1, :] - Fy_full[:, 1:-2, :]) / g.dy

    tend[:, :, 0] = 0.0        # ground: w pinned by BC
    return tend


def momentum_advection_tendency(state: FlowState, grid: Grid, order: int = 2,
                                periodic: bool | None = None, return_fluxes: bool = False):
    """Return ``(du_dt, dv_dt, dw_dt)`` = -div(u u_i) on the staggered faces.

    ``periodic`` defaults to the grid's own lateral periodicity.  ``order`` 1
    (upwind) or 2 (MUSCL/minmod).  With ``return_fluxes=True`` also returns a dict of the
    cell-centred **normal momentum fluxes** (``Fx_u`` = ``Uc·u`` at (nx,ny,nz), ``Fy_v`` =
    ``Vc·v``) — the divergence-form fluxes the interface reflux (``momentum_reflux``) needs;
    they are the very arrays the tendency divergences use, so recording them is exact."""
    if periodic is None:
        periodic = getattr(grid, "periodic", False)
    fluxes = {} if return_fluxes else None
    du = _u_tendency(state, grid, order, periodic, fluxes)
    dv = _v_tendency(state, grid, order, periodic, fluxes)
    dw = _w_tendency(state, grid, order, periodic)
    if return_fluxes:
        return du, dv, dw, fluxes
    return du, dv, dw


def add_momentum_advection(state: FlowState, grid: Grid, dt: float,
                           order: int = 2, periodic: bool | None = None) -> None:
    """Advance u, v, w in place by the explicit momentum-advection tendency."""
    du, dv, dw = momentum_advection_tendency(state, grid, order, periodic)
    state.u += dt * du
    state.v += dt * dv
    state.w += dt * dw


__all__ = ["momentum_advection_tendency", "add_momentum_advection"]
