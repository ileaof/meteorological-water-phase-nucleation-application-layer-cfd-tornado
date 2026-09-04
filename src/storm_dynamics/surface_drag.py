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


def log_law_drag_coefficient(z1_m: float, z0_m: float, kappa: float = 0.4) -> float:
    """Neutral log-law drag coefficient ``C_d = (kappa / ln(z1/z0))^2`` at the first cell-centre
    height ``z1``.  Unlike a fixed bulk ``C_d`` (calibrated for one z1), this stays consistent as
    the near-surface mesh is refined toward the corner-flow layer -- the surface-sensitivity study
    showed a fixed bulk coefficient over-damps the lowest level once z1 is resolved."""
    import math
    z1 = max(float(z1_m), 2.0 * float(z0_m))          # keep the log positive/well-posed
    return (kappa / math.log(z1 / float(z0_m))) ** 2


def effective_drag_coefficient(grid: Grid, drag: SurfaceDragConfig) -> float:
    """The C_d actually applied: the log-law value at this mesh's first cell-centre height when
    ``drag.use_log_law``, else the configured bulk constant.

    **Height dependence is correct physics, and a confound for a resolution experiment.**
    ``C_d = (kappa/ln(z1/z0))^2`` genuinely depends on the reference height, which is the whole
    point of the log-law closure -- it is what keeps the surface stress consistent as the mesh is
    refined.  But it means a dz1 sweep varies the drag coefficient at the same time as the
    resolution: measured at z0 = 0.1 m, C_d is 0.0045 at z1 = 39.9 m (nz=48, zs=1.05), 0.0072 at
    11.1 m (nz=64, zs=1.06) and 0.0146 at 2.73 m (nz=64, zs=1.09) -- a 3.2x spread.  A run at
    dz1 = 5.5 m therefore came back with a much lower peak surface wind than its dz1 = 79.8 m
    counterpart, and the tripled drag is at least as plausible an explanation as the resolution.

    ``drag.log_law_reference_height_m`` pins the evaluation height, so every member of a dz1
    sweep applies an IDENTICAL C_d and the only variable is the mesh.  ``None`` (default) keeps
    the mesh-following behaviour, which remains the right choice for a production run."""
    if not getattr(drag, "use_log_law", False):
        return float(drag.C_d)
    z_ref = getattr(drag, "log_law_reference_height_m", None)
    if z_ref is None:
        z_ref = (float(grid.zc[0]) if not hasattr(grid.zc, "get")
                 else float(grid.backend.to_cpu(grid.zc)[0]))
    return log_law_drag_coefficient(float(z_ref), drag.roughness_length_m,
                                    getattr(drag, "kappa", 0.4))


def _face_speed(state: FlowState, grid: Grid, k: int, U_min: float):
    """Lowest-common helper: |V| at level ``k`` interpolated onto the u- and v-faces."""
    xp = grid.xp
    uc = 0.5 * (state.u[:-1, :, k] + state.u[1:, :, k])
    vc = 0.5 * (state.v[:, :-1, k] + state.v[:, 1:, k])
    sc = xp.maximum(xp.sqrt(uc ** 2 + vc ** 2), U_min)
    su = xp.zeros_like(state.u[:, :, k]); sv = xp.zeros_like(state.v[:, :, k])
    su[1:-1, :] = 0.5 * (sc[:-1, :] + sc[1:, :])
    sv[:, 1:-1] = 0.5 * (sc[:, :-1] + sc[:, 1:])
    if getattr(grid, "periodic", False):
        w = 0.5 * (sc[-1, :] + sc[0, :]); su[0, :] = w; su[-1, :] = w
        w = 0.5 * (sc[:, -1] + sc[:, 0]); sv[:, 0] = w; sv[:, -1] = w
    else:
        su[0, :] = sc[0, :]; su[-1, :] = sc[-1, :]
        sv[:, 0] = sc[:, 0]; sv[:, -1] = sc[:, -1]
    return su, sv


def apply_surface_stress_divergence(state: FlowState, grid: Grid, dt: float,
                                    drag: SurfaceDragConfig, C_d: float) -> None:
    """Surface stress as a **flux divergence through a resolved surface layer** (in place).

    The bulk form puts the whole stress into the lowest cell, so the sink rate ``C_d|V|/dz1``
    diverges as the near-surface mesh is refined and strips the *tangential* wind — measured as the
    blocker of surface connection.  Here the stress decays linearly over a **physical** depth ``h``,
    ``tau(z) = tau_s (1 - z/h)``, so the tendency ``-(1/rho) dtau/dz = -tau_s/h`` is uniform through
    the layer and **independent of the mesh**.  Applied implicitly per level for stability."""
    xp = grid.xp
    h = max(float(drag.surface_layer_depth_m), 1.0)
    zc = grid.zc if not hasattr(grid.zc, "get") else grid.backend.to_cpu(grid.zc)
    import numpy as _np
    zc = _np.asarray(zc)
    for k in range(grid.nz):
        z = float(zc[k])
        if z >= h:
            break
        su, sv = _face_speed(state, grid, k, drag.U_min)
        r = C_d / h                                   # tau_s/h with tau_s = C_d |V| V
        state.u[:, :, k] /= (1.0 + r * su * dt)
        state.v[:, :, k] /= (1.0 + r * sv * dt)


def apply_surface_drag(state: FlowState, grid: Grid, dt: float,
                       drag: SurfaceDragConfig) -> None:
    """Retard the lowest model level by the bulk drag law (implicit, in place).

    With ``drag.use_log_law`` the coefficient is the height-consistent neutral log-law value at the
    actual first cell-centre height (see :func:`effective_drag_coefficient`)."""
    if not drag.enabled:
        return
    C_d = effective_drag_coefficient(grid, drag)
    if C_d <= 0.0:
        return
    if getattr(drag, "stress_divergence", False):
        apply_surface_stress_divergence(state, grid, dt, drag, C_d)
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
    ru = C_d * su / dz0
    rv = C_d * sv / dz0
    state.u[:, :, 0] /= (1.0 + ru * dt)
    state.v[:, :, 0] /= (1.0 + rv * dt)


__all__ = ["apply_surface_drag", "log_law_drag_coefficient", "effective_drag_coefficient",
           "apply_surface_stress_divergence"]
