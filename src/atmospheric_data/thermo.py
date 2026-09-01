"""SI thermodynamic transforms for real-data ingestion (ROADMAP §3a — real cases).

All functions are pure NumPy in **SI units** (K, Pa, kg/kg, m, m/s) and reuse the model's own
constants from :mod:`meteorological_flow.thermodynamics`, so the ingested base state is exactly
consistent with the solver.  Nothing here touches the engine; it converts real fields
(HRRR/ERA5/soundings) into the model's prognostic variables.

Key relations (documented in the task spec):
    rho   = p / (R_d T_v)
    T_v   = T (1 + 0.61 q_v - q_l - q_i)
    theta = T (p0/p)^(R_d/c_p)
"""
from __future__ import annotations

import numpy as np

from meteorological_flow import thermodynamics as th

Rd = float(th.R_d); cp = float(th.cp_d); g0 = float(th.g0)
P0 = float(th.P0_REF); Lv = float(th.Lv); EPS = float(th.EPS)
KAPPA = Rd / cp


def virtual_temperature(T, qv, ql=0.0, qi=0.0):
    """Virtual temperature [K]: ``T_v = T (1 + 0.61 q_v - q_l - q_i)`` (SI mixing ratios)."""
    return np.asarray(T, float) * (1.0 + 0.61 * np.asarray(qv, float)
                                   - np.asarray(ql, float) - np.asarray(qi, float))


def density(p, Tv):
    """Moist air density [kg/m^3] from pressure [Pa] and virtual temperature [K]."""
    return np.asarray(p, float) / (Rd * np.asarray(Tv, float))


def potential_temperature(T, p):
    """Potential temperature [K]: ``theta = T (P0/p)^kappa``."""
    return np.asarray(T, float) * (P0 / np.asarray(p, float)) ** KAPPA


def temperature_from_theta(theta, p):
    """Absolute temperature [K] from potential temperature and pressure [Pa]."""
    return np.asarray(theta, float) * (np.asarray(p, float) / P0) ** KAPPA


def theta_v(theta, qv, ql=0.0, qi=0.0):
    """Virtual potential temperature [K]."""
    return np.asarray(theta, float) * (1.0 + 0.61 * np.asarray(qv, float)
                                       - np.asarray(ql, float) - np.asarray(qi, float))


def saturation_vapor_pressure(T):
    """Saturation vapour pressure over water [Pa] (engine's formula, vectorised)."""
    return th.psat_water(np.asarray(T, float))


def specific_humidity_from_vapor_pressure(e, p):
    """Specific humidity [kg/kg] from vapour pressure ``e`` and pressure ``p`` [Pa].

    ``q_v = eps e / (p - (1-eps) e)`` (exact specific humidity, not the mixing-ratio form)."""
    e = np.minimum(np.asarray(e, float), 0.99 * np.asarray(p, float))
    return EPS * e / (np.asarray(p, float) - (1.0 - EPS) * e)


def specific_humidity_from_dewpoint(Td, p):
    """Specific humidity [kg/kg] from dew-point temperature [K] and pressure [Pa]."""
    return specific_humidity_from_vapor_pressure(saturation_vapor_pressure(Td), p)


def specific_humidity_from_rh(rh, T, p):
    """Specific humidity [kg/kg] from relative humidity (fraction 0..1) at T [K], p [Pa]."""
    return specific_humidity_from_vapor_pressure(np.clip(np.asarray(rh, float), 0.0, 1.0)
                                                 * saturation_vapor_pressure(T), p)


def mixing_ratio_from_specific_humidity(q):
    """Water-vapour mixing ratio [kg/kg] from specific humidity: ``r = q/(1-q)``."""
    q = np.asarray(q, float)
    return q / np.clip(1.0 - q, 1e-12, None)


def specific_humidity_from_mixing_ratio(r):
    """Specific humidity [kg/kg] from mixing ratio: ``q = r/(1+r)``."""
    r = np.asarray(r, float)
    return r / (1.0 + r)


def relative_humidity(qv, T, p):
    """Relative humidity (fraction) from specific humidity, T [K], p [Pa]."""
    e = np.asarray(qv, float) * np.asarray(p, float) / (EPS + (1.0 - EPS) * np.asarray(qv, float))
    return e / saturation_vapor_pressure(T)


def geopotential_to_height(geopotential):
    """Geometric height [m] from geopotential [m^2/s^2] (divide by g0).  Accepts geopotential
    HEIGHT [m] already (returned unchanged) only if the caller passes it as such -- callers
    must know which they have; this does the ``phi/g`` conversion."""
    return np.asarray(geopotential, float) / g0


def brunt_vaisala_squared(theta, z):
    """``N^2 = (g/theta) d(theta)/dz`` [1/s^2] on heights ``z`` [m] (1-D profile)."""
    theta = np.asarray(theta, float); z = np.asarray(z, float)
    return g0 / np.clip(theta, 1e-3, None) * np.gradient(theta, z)


def hypsometric_height(p_desc, Tv, z0=0.0):
    """Geometric height [m] of pressure levels via the hypsometric equation.

    ``p_desc`` pressure [Pa] ordered **descending** (surface first, p high -> top, p low);
    ``Tv`` virtual temperature [K] per level; ``z0`` surface height.  Ascending height:
    ``z_k = z_{k-1} + (R_d T_v_mean / g) ln(p_{k-1}/p_k)`` -- the physically correct
    pressure->height map used to put pressure-level analyses on the model's height coordinate."""
    p = np.asarray(p_desc, float); Tv = np.asarray(Tv, float)
    z = np.empty(p.size); z[0] = float(z0)
    for k in range(1, p.size):
        Tv_mean = 0.5 * (Tv[k] + Tv[k - 1])
        z[k] = z[k - 1] + Rd * Tv_mean / g0 * np.log(p[k - 1] / p[k])
    return z


def hydrostatic_base_pressure(z, theta0, qv0, p_sfc):
    """Integrate the hydrostatic base pressure upward for a given theta0(z), qv0(z) column.

    ``dp0/dz = -rho0 g``, ``rho0 = p0/(R_d T_v0)``, ``T0 = theta0 (p0/P0)^kappa`` -- the same
    fixed-point integration the model's :func:`meteorological_flow.base_state.build_base_state`
    uses, so a base state built here is discretely consistent.  Returns ``(p0, T0, rho0)``."""
    z = np.asarray(z, float); theta0 = np.asarray(theta0, float); qv0 = np.asarray(qv0, float)
    nz = z.size
    p0 = np.empty(nz); T0 = np.empty(nz)
    p_prev, z_prev = float(p_sfc), 0.0
    for k in range(nz):
        dz = z[k] - z_prev; p_new = p_prev
        for _ in range(4):
            Tk = float(theta0[k]) * (p_new / P0) ** KAPPA
            Tv = Tk * (1.0 + 0.61 * float(qv0[k]))
            p_new = p_prev * float(np.exp(-g0 * dz / (Rd * Tv)))
        p0[k] = p_new
        T0[k] = float(theta0[k]) * (p0[k] / P0) ** KAPPA
        p_prev, z_prev = p_new, z[k]
    rho0 = p0 / (Rd * T0 * (1.0 + 0.61 * qv0))
    return p0, T0, rho0
