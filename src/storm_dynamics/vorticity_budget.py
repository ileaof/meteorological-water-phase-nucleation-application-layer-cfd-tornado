"""Vertical-vorticity BUDGET -- the instrument that says *why* a column rotates.

The prognostic core carries tilting/stretching/baroclinicity only *implicitly* inside the
`(u.grad)u` advection; this module makes each production term of the vertical-vorticity equation an
explicit, diagnosable field so tornadogenesis can be attributed to a mechanism rather than asserted.

For vertical vorticity ``zeta = dv/dx - du/dy`` the material budget is

    D zeta / Dt = ADVECTION + STRETCHING + TILTING + BAROCLINIC + DIVERGENCE + DIFFUSION (+ FRICTION)

with, at cell centres (all `grid.xp`-generic; periodic-aware in x,y, one-sided at the z walls):

* **advection**   A  = -( u dzeta/dx + v dzeta/dy + w dzeta/dz )
* **stretching**  S  =  (zeta + f) dw/dz                     (the spec's S_zeta, planetary f optional)
* **tilting**     T  =  xi dw/dx + eta dw/dy                 (tilts horizontal vorticity into vertical)
* **baroclinic**  B  =  (1/rho^2) ( drho/dx dp/dy - drho/dy dp/dx )     (vertical component of
                                                             (1/rho^2) grad rho x grad p)
* **divergence**  Dv = -(zeta + f)( du/dx + dv/dy + dw/dz )  (compressibility; ~0 incompressible,
                                                             small anelastic residual)
* **diffusion**   K  =  div( Km grad zeta )                  (SGS-stress proxy of curl(div tau / rho))
* **friction**    Fr =  curl_z( friction_force )            (surface drag, when a force field is given)

Note on the tilting sign: the spec writes ``T_zeta = dw/dx dv/dz - dw/dy du/dz``.  With
``zeta = dv/dx - du/dy`` the budget-consistent tilting term is ``xi dw/dx + eta dw/dy`` which equals
``-T_zeta(spec)``; :func:`tilting_Tzeta` returns the spec's literal expression for reference, while
the budget uses the sign that makes the terms sum to D zeta/Dt (verified by the closure test).
"""
from __future__ import annotations

import numpy as np

from .rotation import _centered_velocity, vorticity_3d, vertical_vorticity

_TERMS = ("advection", "stretching", "tilting", "baroclinic", "divergence", "diffusion", "friction")


def _grads(field, grid):
    return grid._central_x(field), grid._central_y(field), grid._central_z(field)


def horizontal_vorticity(uc, vc, wc, grid, shear_only=False):
    """Horizontal vorticity (xi, eta) at cell centres.  ``shear_only`` gives the reduced
    environmental form ``omega_h = (-dv/dz, du/dz)`` (drops the dw terms)."""
    dvdz = grid._central_z(vc); dudz = grid._central_z(uc)
    if shear_only:
        return -dvdz, dudz
    dwdy = grid._central_y(wc); dwdx = grid._central_x(wc)
    return dwdy - dvdz, dudz - dwdx


def streamwise_crosswise(uc, vc, xi, eta, grid, storm_motion=(0.0, 0.0)):
    """Decompose the horizontal vorticity (xi, eta) into STREAMWISE (along the storm-relative
    horizontal wind) and CROSSWISE (perpendicular) components [1/s].  Returns (streamwise,
    crosswise), each a 3-D field; streamwise vorticity is the ingredient supercells tilt into a
    *rotating* (rather than merely tilting) updraft."""
    xp = grid.xp
    cu, cv = storm_motion
    ur = uc - cu; vr = vc - cv
    speed = xp.sqrt(ur * ur + vr * vr) + 1e-12
    sx = ur / speed; sy = vr / speed               # unit vector along storm-relative flow
    streamwise = xi * sx + eta * sy                # projection onto the flow
    crosswise = -xi * sy + eta * sx                # projection onto the perpendicular
    return streamwise, crosswise


def tilting_Tzeta(uc, vc, wc, grid):
    """The spec's literal tilting expression T_zeta = dw/dx dv/dz - dw/dy du/dz [1/s^2]
    (see the module note on sign; the budget's tilting term is its negative)."""
    dwdx = grid._central_x(wc); dwdy = grid._central_y(wc)
    dvdz = grid._central_z(vc); dudz = grid._central_z(uc)
    return dwdx * dvdz - dwdy * dudz


def tilting_efficiency(uc, vc, wc, grid):
    """Factorise the tilting term to explain *why* horizontal vorticity does (or does not) become
    vertical rotation.  Tilting = omega_h . grad_h w = |omega_h| |grad_h w| cos(theta), so a weak
    tilting rate is one of three distinct failures:

    * ``omega_h`` small  -> no horizontal vorticity available (weak shear / weak cold pool),
    * ``grad_h w`` small -> no low-level updraft gradient to do the tilting,
    * ``alignment`` ~ 0  -> both present but geometrically MISALIGNED (the updraft gradient is not
      oriented to lift the vortex lines).

    Returns a dict of fields: ``omega_h``, ``grad_h_w``, ``tilting``, ``alignment`` (cos theta)."""
    xp = grid.xp
    xi, eta = horizontal_vorticity(uc, vc, wc, grid)
    dwdx = grid._central_x(wc); dwdy = grid._central_y(wc)
    omh = xp.sqrt(xi * xi + eta * eta)
    gw = xp.sqrt(dwdx * dwdx + dwdy * dwdy)
    tilt = xi * dwdx + eta * dwdy
    return {"omega_h": omh, "grad_h_w": gw, "tilting": tilt,
            "alignment": tilt / (omh * gw + 1e-20)}


def baroclinic_horizontal_generation(rho, grid, g0=9.81):
    """Hydrostatic baroclinic generation of **horizontal** vorticity from the density (cold-pool)
    field: with the solenoidal term (1/rho^2) grad rho x grad p and hydrostatic
    ``grad p ~ (0,0,-rho g)`` this is ``(-g/rho drho/dy, g/rho drho/dx)`` [1/s^2].

    This is the cold pool's vorticity SOURCE -- the horizontal vorticity generated along the
    buoyancy (density) gradient, which the updraft then TILTS into the vertical.  It is the
    physically-meaningful baroclinic diagnostic for tornadogenesis: the *direct* vertical-vorticity
    baroclinic term (:func:`budget_from_velocity`'s ``baroclinic``) vanishes under hydrostatic
    balance, and needs the dynamic pressure which the low-memory anelastic solver does not persist
    in ``state.p`` on large (nest) grids.  Computable from ``rho`` alone.  Returns (Gx, Gy, |G|)."""
    xp = grid.xp
    rx = grid._central_x(rho); ry = grid._central_y(rho)
    Gx = -g0 * ry / rho; Gy = g0 * rx / rho
    return Gx, Gy, xp.sqrt(Gx * Gx + Gy * Gy)


def budget_from_velocity(uc, vc, wc, grid, *, rho=None, p=None, Km=None, f=0.0,
                         friction_force=None):
    """Vertical-vorticity budget term FIELDS from centred velocity (+ optional rho, p, Km).

    ``rho``/``p`` (cell-centred total density & pressure) enable the BAROCLINIC term; without them
    it is returned as zeros (so the routine still works on saved velocity-only fields).  ``Km`` (an
    eddy-viscosity field or scalar) enables the DIFFUSION proxy.  ``f`` adds planetary vorticity to
    the stretching/divergence terms.  ``friction_force`` (fx, fy) [m/s^2] enables the FRICTION term.
    Returns a dict of 3-D arrays keyed by the term name plus ``zeta`` and ``tendency`` (the RHS sum).
    """
    xp = grid.xp
    zeta = grid._central_x(vc) - grid._central_y(uc)
    zx, zy, zz = _grads(zeta, grid)
    dudx = grid._central_x(uc); dvdy = grid._central_y(vc); dwdz = grid._central_z(wc)
    dwdx = grid._central_x(wc); dwdy = grid._central_y(wc)
    dvdz = grid._central_z(vc); dudz = grid._central_z(uc)
    xi = dwdy - dvdz; eta = dudz - dwdx
    zabs = zeta + f

    terms = {}
    terms["advection"] = -(uc * zx + vc * zy + wc * zz)
    terms["stretching"] = zabs * dwdz
    terms["tilting"] = xi * dwdx + eta * dwdy
    terms["divergence"] = -zabs * (dudx + dvdy + dwdz)
    if rho is not None and p is not None:
        rx = grid._central_x(rho); ry = grid._central_y(rho)
        px = grid._central_x(p); py = grid._central_y(p)
        terms["baroclinic"] = (rx * py - ry * px) / (rho * rho)
    else:
        terms["baroclinic"] = xp.zeros_like(zeta)
    if Km is not None:
        # div(Km grad zeta): allow Km scalar or field
        Kx = Km if xp.isscalar(Km) or getattr(Km, "ndim", 0) == 0 else Km
        gx, gy, gz = _grads(zeta, grid)
        terms["diffusion"] = (grid._central_x(Kx * gx) + grid._central_y(Kx * gy)
                              + grid._central_z(Kx * gz))
    else:
        terms["diffusion"] = xp.zeros_like(zeta)
    if friction_force is not None:
        fx, fy = friction_force
        terms["friction"] = grid._central_x(fy) - grid._central_y(fx)
    else:
        terms["friction"] = xp.zeros_like(zeta)

    terms["zeta"] = zeta
    terms["tendency"] = sum(terms[k] for k in _TERMS)
    return terms


def zeta_budget(state, grid, *, Km=None, f=0.0, friction_force=None):
    """Convenience: the vertical-vorticity budget from a live :class:`FlowState`.

    Uses the diagnosed ``state.rho`` and ``state.P_total`` for the baroclinic term (call
    ``state.diagnose(cfg)`` first).  ``Km`` may be the simulation's ``self._Km`` eddy viscosity."""
    uc, vc, wc = _centered_velocity(state, grid)
    rho = getattr(state, "rho", None)
    p = getattr(state, "P_total", None)
    return budget_from_velocity(uc, vc, wc, grid, rho=rho, p=p, Km=Km, f=f,
                                friction_force=friction_force)


def _layer_mask(grid, z_lo, z_hi):
    z = np.asarray(grid.backend.to_cpu(grid.zc))
    return (z >= z_lo) & (z < z_hi)


def budget_layer_summary(terms, grid, z_lo=0.0, z_hi=1000.0):
    """Per-term signed extremum + magnitude mean within a height layer -- the compact numbers for
    the JSON report: which mechanism dominates the vertical-vorticity production there."""
    xp = grid.xp
    mask = xp.asarray(_layer_mask(grid, z_lo, z_hi))
    out = {"z_lo_m": float(z_lo), "z_hi_m": float(z_hi)}
    for k in _TERMS + ("tendency",):
        layer = terms[k][:, :, mask]
        if layer.size:
            out[k + "_max"] = float(xp.max(layer))
            out[k + "_min"] = float(xp.min(layer))
            out[k + "_absmean"] = float(xp.mean(xp.abs(layer)))
        else:
            out[k + "_max"] = out[k + "_min"] = out[k + "_absmean"] = 0.0
    return out


def dominant_mechanism(terms, grid, z_lo=0.0, z_hi=1000.0):
    """Name the production term with the largest layer-mean magnitude (baroclinic vs tilting vs
    stretching ...) -- the one-word answer to 'what is making this column spin here?'."""
    s = budget_layer_summary(terms, grid, z_lo, z_hi)
    prod = {k: s[k + "_absmean"] for k in ("stretching", "tilting", "baroclinic", "friction")}
    name = max(prod, key=prod.get)
    return name, prod


__all__ = [
    "horizontal_vorticity", "streamwise_crosswise", "tilting_Tzeta",
    "baroclinic_horizontal_generation", "tilting_efficiency",
    "budget_from_velocity", "zeta_budget", "budget_layer_summary", "dominant_mechanism",
]
