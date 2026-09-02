"""Two-scale temperature gradient: the resolved (macro) CFD gradient vs the sub-grid (micro)
gradient that the phase-change / nucleation physics actually feels.

The macroscopic gradient ``|grad T|_macro`` is the cell-scale finite difference (fronts, cold
pools, cloud edges).  The nucleation kernel, however, responds to the temperature gradient in the
immediate neighbourhood of droplets / ice interfaces, which lives far below the grid scale.  We do
NOT paste the micro gradient into the resolved energy equation; instead we provide a documented
**sub-grid closure** that estimates it from resolved quantities, to drive the shifted-equilibrium
nucleation term only:

    |grad T|_micro = F(|grad T|_macro, eps, q_v, q_c, q_i, r, T, p).

**Closure (Obukhov-Corrsin inertial-convective scaling with a Batchelor cutoff).**  In the
inertial-convective subrange a passive scalar's increments scale as ``dT(r) ~ r^{1/3}`` so its
gradient scales as ``r^{-2/3}``.  Extrapolating the resolved gradient at the grid scale ``Delta``
down to a microphysical scale ``r_eff`` gives

    |grad T|_micro = |grad T|_macro * (Delta / r_eff)^{2/3},

with ``r_eff = max(r_particle, eta_B)`` where ``eta_B = (nu^3/eps)^{1/4} * sqrt(kappa_T/nu)`` is the
Batchelor scale (below which molecular diffusion erases scalar gradients).  ``eps`` (turbulent
dissipation) enters through ``eta_B``; the hydrometeor field sets ``r_particle``; ``T, p`` set the
molecular ``nu, kappa_T``.  The enhancement is floored at 1 and capped (documented) to stay finite.
`grid.xp`-generic.  A pure diagnostic -- it changes no prognostic field here.
"""
from __future__ import annotations

import numpy as np

_NU = 1.5e-5          # kinematic viscosity of air [m^2/s]
_KAPPA_T = 2.0e-5     # thermal diffusivity of air [m^2/s]
_MAX_ENHANCE = 1.0e4  # cap on the micro/macro gradient ratio (documented guard)


def macro_temperature_gradient(state, grid):
    """Resolved cell-scale |grad T| [K/m] (the CFD macroscopic gradient)."""
    if getattr(state, "gradT_mag", None) is not None:
        return state.gradT_mag
    return grid.grad_magnitude(state.T)


def _characteristic_radius(state, grid):
    """A characteristic hydrometeor/interface radius [m] from the local condensate: larger drops
    where rain/graupel dominate, ~cloud-droplet scale in cloud, a molecular floor in clear air."""
    xp = grid.xp
    r = xp.full(grid.center_shape, 1.0e-5)                 # ~10 micron cloud droplet default
    qc = getattr(state, "ql", None)
    if qc is not None:
        r = xp.where(qc > 1e-6, 1.5e-5, r)                 # cloud
    for nm, rad in (("qr", 5e-4), ("qs", 3e-4), ("qg", 2e-3), ("qh", 5e-3)):
        a = getattr(state, nm, None)
        if a is not None:
            r = xp.where(a > 1e-6, rad, r)                 # precipitation -> larger particles
    return r


def batchelor_scale(eps, nu=_NU, kappa_T=_KAPPA_T):
    """Batchelor scale eta_B = (nu^3/eps)^{1/4} * sqrt(kappa_T/nu) [m] (scalar-gradient cutoff)."""
    import numpy as _np
    eps = _np.maximum(eps, 1e-9)
    eta_k = (nu ** 3 / eps) ** 0.25
    return eta_k * (kappa_T / nu) ** 0.5


def micro_temperature_gradient(state, grid, eps=1.0e-3, r_particle=None):
    """Sub-grid micro temperature gradient [K/m] from the closure above.  ``eps`` [m^2/s^3] is the
    turbulent dissipation (scalar or field; the LES can supply a better estimate); ``r_particle``
    overrides the condensate-derived radius.  Returns (grad_micro, enhancement)."""
    xp = grid.xp
    macro = macro_temperature_gradient(state, grid)
    r_part = _characteristic_radius(state, grid) if r_particle is None else r_particle
    eta_b = batchelor_scale(eps)
    r_eff = xp.maximum(r_part, xp.asarray(eta_b))
    Delta = (grid.dx * grid.dy * (grid.dz if not getattr(grid, "stretched", False)
                                  else float(grid.dz_c[0]))) ** (1.0 / 3.0)
    enhancement = xp.clip((Delta / r_eff) ** (2.0 / 3.0), 1.0, _MAX_ENHANCE)
    return macro * enhancement, enhancement


def nucleation_diagnostics(state, grid, eps=1.0e-3) -> dict:
    """Named two-scale + nucleation diagnostics for the report / w-vs-gradient analyses:
    macro & micro temperature gradient, local supersaturation, a nucleation-rate proxy, and the
    latent-heat-release proxy.  The nucleation-rate PROXY is monotone in supersaturation and the
    micro gradient; the *validated* kernel is in meteorological_flow.nucleation_adapter."""
    xp = grid.xp
    macro = macro_temperature_gradient(state, grid)
    micro, enh = micro_temperature_gradient(state, grid, eps=eps)
    S_w = getattr(state, "S_w", None)
    ss = (S_w - 1.0) if S_w is not None else xp.zeros(grid.center_shape)
    ss_pos = xp.maximum(ss, 0.0)
    # proxy: nucleation is steep in supersaturation and enhanced by sharper local gradients
    nuc_proxy = ss_pos ** 3 * (1.0 + micro / (macro + 1e-12) * 0.0 + xp.log1p(enh))
    # latent-heat proxy: condensation rate ~ supersaturation removal (report the available signal)
    lh_proxy = ss_pos * (getattr(state, "ql", 0.0) + getattr(state, "qi", 0.0)
                         if getattr(state, "ql", None) is not None else ss_pos)
    return {
        "macro_temperature_gradient_max_K_m": float(xp.max(macro)),
        "micro_temperature_gradient_max_K_m": float(xp.max(micro)),
        "micro_macro_enhancement_max": float(xp.max(enh)),
        "local_supersaturation_max": float(xp.max(ss)),
        "nucleation_rate_proxy_max": float(xp.max(nuc_proxy)),
        "latent_heat_release_proxy_max": float(xp.max(lh_proxy)),
    }


__all__ = [
    "macro_temperature_gradient", "micro_temperature_gradient", "batchelor_scale",
    "nucleation_diagnostics",
]
