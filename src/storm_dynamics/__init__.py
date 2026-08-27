"""storm_dynamics -- idealised rotating deep-convection (supercell / tornadogenesis) core.

A *fork* of the ``meteorological_flow`` dynamical core that adds the physics a
storm needs in order to **rotate**: fully conservative flux-form momentum
advection (so vorticity can tilt and stretch), f-plane Coriolis, an LES subgrid
closure (in place of the demonstration Rayleigh drag + velocity clip), a
bulk-drag surface layer, and curved-hodograph environmental shear.

Everything *else* is imported unchanged from ``meteorological_flow`` and
``precip_microphysics`` and used as a physics library: the staggered
:class:`~meteorological_flow.grid.Grid`, the Chorin/anelastic
:class:`~meteorological_flow.pressure_solver.PressureSolver`, conservative scalar
transport, moist buoyancy, the bulk microphysics (evaporative cold pool =
baroclinic vorticity source), and the validated nucleation kernel.  The immutable
``met_water_nucleation._engine`` is never touched, so ``--validate`` stays green.

Scope (be honest): this is **idealised** simulation, not operational forecasting.
No data assimilation, no real-event initial/boundary conditions, no observational
verification.  The deliverable is a rotational dynamical core validated against
the *classical* supercell results (Klemp & Wilhelmson 1978; Weisman & Klemp 1982;
Rotunno & Klemp 1985), reusing the repo's already-validated physics.
"""
from __future__ import annotations

__version__ = "0.1.0"

__all__ = ["__version__"]
