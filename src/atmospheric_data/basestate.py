"""Build the model base state from real data (ROADMAP §3a).

The anelastic core evolves perturbations about a hydrostatic, horizontally-uniform base state
``(theta0, qv0, p0, rho0, u0, v0)(z)``.  From the regridded real fields we take the
**horizontal mean** profile at the initial time and re-integrate ``p0`` hydrostatically
(``dp0/dz = -rho0 g``) so the base is discretely consistent with the solver -- then the full
real field decomposes as ``phi = phi0 + phi'`` (task: perturbation decomposition, anelastic
formulation, no mixing of formulations).
"""
from __future__ import annotations

import numpy as np

from . import thermo


def base_state_from_fields(fields, z_model, it=0):
    """Return a ``meteorological_flow.base_state.BaseState`` from regridded model-grid fields
    (dict of ``(time,nz,ny,nx)`` arrays incl. ``theta,qv,u,v``) at time index ``it``."""
    from meteorological_flow.base_state import BaseState
    z = np.asarray(z_model, float)
    theta0 = fields["theta"][it].mean(axis=(1, 2))
    qv0 = (fields["qv"][it].mean(axis=(1, 2)) if "qv" in fields else np.zeros_like(z))
    u0 = (fields["u"][it].mean(axis=(1, 2)) if "u" in fields else np.zeros_like(z))
    v0 = (fields["v"][it].mean(axis=(1, 2)) if "v" in fields else np.zeros_like(z))
    p_sfc = (float(fields["p"][it, 0].mean()) if "p" in fields else 1.0e5)
    p0, T0, rho0 = thermo.hydrostatic_base_pressure(z, theta0, qv0, p_sfc)
    return BaseState(zc=z, theta0=theta0, qv0=qv0, p0=p0, T0=T0, rho0=rho0, u0=u0, v0=v0)


def decompose_perturbations(fields, base, it=0):
    """Split the real fields into base + perturbation about ``base`` at time ``it``.

    Returns ``{theta_prime, qv_prime, u_prime, v_prime, w}`` (w has no base, w0=0), each
    ``(nz,ny,nx)`` -- the initial state the anelastic core steps."""
    z0 = lambda a: np.asarray(a, float)[:, None, None]
    return {"theta_prime": fields["theta"][it] - z0(base.theta0),
            "qv_prime": (fields["qv"][it] - z0(base.qv0)) if "qv" in fields else np.zeros_like(fields["theta"][it]),
            "u_prime": (fields["u"][it] - z0(base.u0)) if "u" in fields else np.zeros_like(fields["theta"][it]),
            "v_prime": (fields["v"][it] - z0(base.v0)) if "v" in fields else np.zeros_like(fields["theta"][it]),
            "w": fields["w"][it] if "w" in fields else np.zeros_like(fields["theta"][it])}
