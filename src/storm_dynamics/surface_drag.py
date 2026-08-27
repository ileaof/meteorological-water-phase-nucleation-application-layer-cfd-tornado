"""Surface layer -- bulk aerodynamic drag law (item 4).

The demonstration core uses free-slip / no-slip at the ground.  Neither produces
the frictional inflow that is essential to low-level rotation: surface drag
retards the near-ground wind, tilting the horizontal (frictionally generated)
vorticity and driving the convergent **corner flow** that intensifies a
near-surface vortex (Rotunno & Klemp 1985; the tornado corner-flow region).

A bulk drag law applies a stress to the lowest model level proportional to the
square of the near-surface wind:

    tau = rho C_d |V| V ,   du/dt|_{k=0} = - C_d |V| u / dz_0 ,

with ``|V| = max(sqrt(u^2+v^2), U_min)`` so a weak near-surface flow still feels
some friction.  ``C_d ~ 0.01-0.02`` over land.  Only the lowest cell is retarded;
the interior LES closure carries the momentum flux upward.  This is treated
implicitly (a per-step exponential decay factor) so a strong low-level wind at a
coarse ``dz_0`` cannot overshoot into instability.
"""
from __future__ import annotations

from meteorological_flow.grid import Grid
from meteorological_flow.state import FlowState

from .config import SurfaceDragConfig


def apply_surface_drag(state: FlowState, grid: Grid, dt: float,
                       drag: SurfaceDragConfig) -> None:
    """Retard the lowest model level by the bulk drag law (implicit, in place)."""
    if not drag.enabled or drag.C_d <= 0.0:
        return
    xp = grid.xp
    dz0 = grid.dz if not getattr(grid, "stretched", False) else float(grid.dz_c[0])
    # lowest-level speed at cell centres
    uc0 = 0.5 * (state.u[:-1, :, 0] + state.u[1:, :, 0])          # (nx,ny)
    vc0 = 0.5 * (state.v[:, :-1, 0] + state.v[:, 1:, 0])          # (nx,ny)
    speed_c = xp.sqrt(uc0 ** 2 + vc0 ** 2)
    speed_c = xp.maximum(speed_c, drag.U_min)
    # drag rate r = C_d |V| / dz0  -> implicit decay exp(-r dt) ~ 1/(1+r dt)
    # interpolate the centre speed onto the u- and v-faces of the lowest level
    su = xp.zeros_like(state.u[:, :, 0])
    su[1:-1, :] = 0.5 * (speed_c[:-1, :] + speed_c[1:, :])
    if getattr(grid, "periodic", False):
        w = 0.5 * (speed_c[-1, :] + speed_c[0, :]); su[0, :] = w; su[-1, :] = w
    else:
        su[0, :] = speed_c[0, :]; su[-1, :] = speed_c[-1, :]
    sv = xp.zeros_like(state.v[:, :, 0])
    sv[:, 1:-1] = 0.5 * (speed_c[:, :-1] + speed_c[:, 1:])
    if getattr(grid, "periodic", False):
        w = 0.5 * (speed_c[:, -1] + speed_c[:, 0]); sv[:, 0] = w; sv[:, -1] = w
    else:
        sv[:, 0] = speed_c[:, 0]; sv[:, -1] = speed_c[:, -1]
    ru = drag.C_d * su / dz0
    rv = drag.C_d * sv / dz0
    state.u[:, :, 0] /= (1.0 + ru * dt)
    state.v[:, :, 0] /= (1.0 + rv * dt)


__all__ = ["apply_surface_drag"]
