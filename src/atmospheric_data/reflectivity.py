"""Diagnostic radar reflectivity from model hydrometeors (ROADMAP §3a validation).

The idealised microphysics does not emit a calibrated reflectivity, so we compute an
**approximate Rayleigh (10 cm) equivalent reflectivity** `Z_e` from the hydrometeor mixing
ratios and compare it to NEXRAD in dBZ.  This is a *diagnostic*, with the assumptions stated
below — not a validated forward operator.

Assumptions (documented, task requirement):
* Marshall-Palmer exponential size distributions, fixed intercepts, Rayleigh scattering, so
  `Z_e = Gamma(7)/(pi^1.75 rho_x^1.75 N0_x^0.75) (rho q_x)^1.75 x 1e18` [mm^6 m^-3] with
  `rho q_x` the water content in **kg m^-3**.  This gives the species coefficients below and,
  e.g., ~41 dBZ for 1 g m^-3 of rain (the expected value).
* Rain N0=8e6 m^-4, rho_w=1000; ice species use the ice dielectric factor (~0.21) with denser,
  fewer particles (dry snow / graupel) -> the tabulated `a_x`.
* Dry snow/graupel (no wet-growth brightband), no attenuation, no melting layer.
Reference: Smith (1984); Marshall & Palmer (1948); WRF `REFL_10CM` family (simplified).
"""
from __future__ import annotations

import numpy as np

# a_x for Z_e[mm^6 m^-3] = a_x (rho q_x [kg m^-3])^1.75  (rain, dry snow, graupel)
_A_RAIN, _A_SNOW, _A_GRAUPEL = 2.49e9, 1.0e9, 5.0e9


def reflectivity_dbz(qr=None, qs=None, qg=None, rho=1.0, floor_dbz=-30.0):
    """Equivalent reflectivity [dBZ] from rain/snow/graupel mixing ratios [kg/kg] and air
    density ``rho`` [kg/m^3].  Missing species are treated as zero (never invented)."""
    def _z(q, a):
        if q is None:
            return 0.0
        wc = np.asarray(rho, float) * np.clip(np.asarray(q, float), 0.0, None)   # kg/m^3
        return a * wc ** 1.75
    Ze = _z(qr, _A_RAIN) + _z(qs, _A_SNOW) + _z(qg, _A_GRAUPEL)     # mm^6 / m^3
    dbz = 10.0 * np.log10(np.clip(Ze, 1e-3, None))
    return np.maximum(dbz, floor_dbz)


def cfd_reflectivity_field(state, grid, rho0_c):
    """3-D simulated reflectivity [dBZ] on the model (nx,ny,nz) grid from the state's
    hydrometeors (``qr,qs,qg`` if present) and the anelastic reference density ``rho0_c(z)``."""
    to = grid.backend.to_cpu
    rho = np.asarray(to(rho0_c), float)[None, None, :]
    get = lambda nm: (np.asarray(to(getattr(state, nm))) if getattr(state, nm, None) is not None else None)
    return reflectivity_dbz(get("qr"), get("qs"), get("qg"), rho=rho)


def reflectivity_at_gates(dbz_xyz, x_model, y_model, z_model, radar):
    """Interpolate a model (x,y,z) reflectivity field to the radar gate positions (same shape
    as ``radar['reflectivity']``) for a gate-for-gate comparison."""
    from scipy.interpolate import RegularGridInterpolator
    gx = np.asarray(radar["x_m"], float) + radar.get("radar_x_m", 0.0)
    gy = np.asarray(radar["y_m"], float) + radar.get("radar_y_m", 0.0)
    gz = np.asarray(radar["alt_m"], float)
    pts = np.stack([np.clip(gz, z_model.min(), z_model.max()).ravel(),
                    np.clip(gy, y_model.min(), y_model.max()).ravel(),
                    np.clip(gx, x_model.min(), x_model.max()).ravel()], -1)
    f = RegularGridInterpolator((z_model, y_model, x_model), np.transpose(dbz_xyz, (2, 1, 0)),
                                bounds_error=False, fill_value=None)
    return f(pts).reshape(gx.shape)
