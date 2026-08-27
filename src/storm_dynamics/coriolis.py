"""f-plane Coriolis force (item 2).

Idealised **f-plane**: a constant Coriolis parameter ``f = 2 Omega sin(lat)`` for a
Tornado-Alley latitude (~36 N), NOT a beta-plane or spherical geometry.  The
horizontal momentum feels

    du/dt = + f v' ,     dv/dt = - f u'

The force is applied to the **perturbation** wind ``(u', v') = (u - u0, v - v0)``,
treating the environmental hodograph ``u0(z), v0(z)`` as a steady (geostrophically
balanced) base state.  Applying f to the full wind instead would make the
unbalanced base-state wind inertially oscillate, since this idealised setup
carries no imposed synoptic pressure-gradient force to balance it (the standard
choice in idealised supercell models, e.g. CM1/Bryan).  Coriolis is a weak effect
at storm scale over the ~1-2 h integration but is included for completeness and
because it sets the sign/º preference consistent with the hemisphere.
"""
from __future__ import annotations

from meteorological_flow.grid import Grid
from meteorological_flow.state import FlowState


def add_coriolis(state: FlowState, grid: Grid, dt: float, f: float,
                 u0_face=None, v0_face=None) -> None:
    """Add the f-plane Coriolis tendency to u, v in place (perturbation wind)."""
    if f == 0.0:
        return
    xp = grid.xp
    up = state.u if u0_face is None else state.u - u0_face
    vp = state.v if v0_face is None else state.v - v0_face
    # interpolate the perturbation cross-component onto each momentum's own faces.
    # v' (on v-faces, nx,ny+1,nz) -> u-faces (nx+1,ny,nz)
    vc = 0.5 * (vp[:, :-1, :] + vp[:, 1:, :])              # (nx,ny,nz) cell centres
    v_at_u = xp.zeros(grid.u_shape)
    v_at_u[1:-1, :, :] = 0.5 * (vc[:-1, :, :] + vc[1:, :, :])
    if getattr(grid, "periodic", False):
        wrap = 0.5 * (vc[-1, :, :] + vc[0, :, :])
        v_at_u[0, :, :] = wrap; v_at_u[-1, :, :] = wrap
    else:
        v_at_u[0, :, :] = vc[0, :, :]; v_at_u[-1, :, :] = vc[-1, :, :]
    # u' (on u-faces) -> v-faces
    uc = 0.5 * (up[:-1, :, :] + up[1:, :, :])             # (nx,ny,nz)
    u_at_v = xp.zeros(grid.v_shape)
    u_at_v[:, 1:-1, :] = 0.5 * (uc[:, :-1, :] + uc[:, 1:, :])
    if getattr(grid, "periodic", False):
        wrap = 0.5 * (uc[:, -1, :] + uc[:, 0, :])
        u_at_v[:, 0, :] = wrap; u_at_v[:, -1, :] = wrap
    else:
        u_at_v[:, 0, :] = uc[:, 0, :]; u_at_v[:, -1, :] = uc[:, -1, :]
    state.u += dt * f * v_at_u
    state.v -= dt * f * u_at_v


__all__ = ["add_coriolis"]
