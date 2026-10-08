"""Environmental soundings for the rotating storm -- curved-hodograph shear.

The thermodynamic column (theta0, qv0, p0, T0, rho0) is reused *unchanged* from
:func:`meteorological_flow.base_state.weisman_klemp` (the standard CAPE ~2000-2500
J/kg idealised sounding).  Only the winds are replaced: the demonstration
:func:`weisman_klemp` imposes a purely *unidirectional* ``u0(z)`` (``v0 = 0``),
which carries no streamwise (storm-relative) vorticity, so it can split a storm
but cannot preferentially spin up low-level rotation.

Here we populate BOTH ``u0(z)`` and ``v0(z)`` to give a **curved** hodograph:

* ``"unidirectional"`` -- straight hodograph (``v0 = 0``); the Klemp-Wilhelmson /
  Weisman-Klemp storm-splitting reference (M1).
* ``"quarter_circle"`` -- the wind vector turns through a quarter circle over the
  lowest ``z_turn`` (veering with height), then is constant aloft.  The curvature
  is the source of storm-relative helicity that selects the right-moving
  supercell and feeds low-level rotation (Rotunno & Klemp 1982, 1985) -- the M2
  target.

References: Weisman & Klemp (1982); Rotunno & Klemp (1982, 1985); Davies-Jones
(1984, storm-relative helicity); Bunkers et al. (2000, storm-motion estimate).
"""
from __future__ import annotations

import numpy as np

from meteorological_flow import thermodynamics as th
from meteorological_flow.base_state import BaseState, weisman_klemp

from .config import HodographConfig


def _quarter_circle_winds(z: np.ndarray, hodo: HodographConfig):
    """Quarter-circle hodograph winds (u0, v0) [m/s] on heights ``z``.

    The wind traces a quarter circle of radius ``U_max`` from the origin at the
    surface to ``(U_max, -U_max)`` at ``z_turn``.  The wind vector veers (turns
    clockwise) with height, giving strong **positive** storm-relative helicity --
    the Northern-Hemisphere right-moving supercell configuration -- then holds
    constant above ``z_turn``.
    """
    R = hodo.U_max
    zt = max(hodo.z_turn, 1.0)
    beta = 0.5 * np.pi * np.clip(z / zt, 0.0, 1.0)
    u0 = R * np.sin(beta)
    v0 = -R * (1.0 - np.cos(beta))
    # above the turning layer the winds are constant (u,v held at their z_turn value)
    aloft = z > zt
    u0[aloft] = R * np.sin(0.5 * np.pi)
    v0[aloft] = -R * (1.0 - np.cos(0.5 * np.pi))
    return u0, v0


def _unidirectional_winds(z: np.ndarray, hodo: HodographConfig):
    """Straight (unidirectional) hodograph: u0 = U_max tanh(z/u_half), v0 = 0."""
    u0 = hodo.U_max * np.tanh(z / max(hodo.u_half, 1.0))
    v0 = np.zeros_like(z)
    return u0, v0


def build_sounding(grid, hodo: HodographConfig | None = None, **wk_kwargs) -> BaseState:
    """Return a :class:`BaseState` with the WK thermodynamics + a curved hodograph.

    ``wk_kwargs`` are forwarded to :func:`weisman_klemp` (e.g. ``qv_sfc``,
    ``theta_sfc``) so CAPE can be tuned; the wind arguments of ``weisman_klemp``
    are overridden here.
    """
    hodo = hodo or HodographConfig()
    base = weisman_klemp(grid, u_shear=0.0, **wk_kwargs)   # thermodynamics only
    z = np.asarray(base.zc, dtype=float)
    if hodo.kind == "unidirectional":
        u0, v0 = _unidirectional_winds(z, hodo)
    elif hodo.kind == "quarter_circle":
        u0, v0 = _quarter_circle_winds(z, hodo)
    else:
        raise ValueError("unknown hodograph kind %r" % hodo.kind)
    base.u0 = u0
    base.v0 = v0
    return base


# ---------------------------------------------------------------------------
# ROADMAP §3a -- REAL initial conditions: ingest an observed sounding
# ---------------------------------------------------------------------------
def from_observed_sounding(grid, *, pressure_hPa, height_m, temperature_C,
                           dewpoint_C=None, qv_kgkg=None,
                           wind_dir_deg=None, wind_speed_ms=None, u_ms=None, v_ms=None,
                           p_sfc_Pa=None) -> BaseState:
    """Build a hydrostatic :class:`BaseState` on ``grid`` from a **real observed sounding**
    (ROADMAP §3a -- the first step toward a real-atmosphere environment, replacing the
    analytic Weisman-Klemp column).

    Radiosonde columns (each a 1-D array over the reported levels):

    * ``pressure_hPa`` [hPa], ``height_m`` [m AGL], ``temperature_C`` [degC];
    * moisture as EITHER ``dewpoint_C`` [degC] OR ``qv_kgkg`` [kg/kg];
    * winds as EITHER (``wind_dir_deg`` [deg, met. convention -- direction FROM],
      ``wind_speed_ms`` [m/s]) OR (``u_ms``, ``v_ms``) [m/s].

    The observed column is converted to ``(theta, qv, u, v)`` and interpolated onto the model
    heights; the pressure is then **re-integrated hydrostatically on the model grid** (exactly
    as :func:`meteorological_flow.base_state.build_base_state` / ``weisman_klemp`` do), so the
    base state is *discretely balanced on the stretched mesh* rather than a raw interpolation
    of the reported pressures.  Returns a ``BaseState`` ready for ``StormSimulation(scfg,
    base=...)``.  Winds above the top reported level are held constant (``np.interp`` clamp)."""
    p_obs = np.asarray(pressure_hPa, float) * 100.0                 # hPa -> Pa
    z_obs = np.asarray(height_m, float)
    T_obs = np.asarray(temperature_C, float) + 273.15              # degC -> K
    order = np.argsort(z_obs)                                      # ascending height
    p_obs, z_obs, T_obs = p_obs[order], z_obs[order], T_obs[order]
    n = z_obs.size
    theta_obs = np.array([th.theta_from_T(float(T_obs[k]), float(p_obs[k]), th.P0_REF)
                          for k in range(n)])
    if qv_kgkg is not None:
        qv_obs = np.asarray(qv_kgkg, float)[order]
    elif dewpoint_C is not None:
        Td = np.asarray(dewpoint_C, float)[order] + 273.15
        qv_obs = np.array([th.q_v_from_p_v(min(float(th.psat_water(float(Td[k]))),
                                               0.99 * float(p_obs[k])), float(p_obs[k]))
                           for k in range(n)])
    else:
        raise ValueError("supply moisture as dewpoint_C or qv_kgkg")
    if u_ms is not None and v_ms is not None:
        u_obs = np.asarray(u_ms, float)[order]; v_obs = np.asarray(v_ms, float)[order]
    elif wind_dir_deg is not None and wind_speed_ms is not None:
        d = np.deg2rad(np.asarray(wind_dir_deg, float)[order])
        s = np.asarray(wind_speed_ms, float)[order]
        u_obs = -s * np.sin(d); v_obs = -s * np.cos(d)            # met. convention (FROM)
    else:
        raise ValueError("supply winds as (wind_dir_deg, wind_speed_ms) or (u_ms, v_ms)")

    zc = np.asarray(grid.backend.to_cpu(grid.zc), float)
    theta0 = np.interp(zc, z_obs, theta_obs)
    qv0 = np.maximum(np.interp(zc, z_obs, qv_obs), 0.0)
    u0 = np.interp(zc, z_obs, u_obs); v0 = np.interp(zc, z_obs, v_obs)

    p_sfc = float(p_sfc_Pa) if p_sfc_Pa is not None else float(np.interp(0.0, z_obs, p_obs))
    nz = zc.size; p0 = np.empty(nz); T0 = np.empty(nz)
    p_prev, z_prev = p_sfc, 0.0
    for k in range(nz):                                           # hydrostatic re-integration
        dz = zc[k] - z_prev; p_new = p_prev
        for _ in range(3):                                        # fixed point: T_v(p)
            Tk = float(th.T_from_theta(theta0[k], p_new, th.P0_REF))
            Tv = Tk * (1.0 + 0.61 * qv0[k])
            p_new = p_prev * float(np.exp(-th.g0 * dz / (th.R_d * Tv)))
        p0[k] = p_new; T0[k] = float(th.T_from_theta(theta0[k], p0[k], th.P0_REF))
        p_prev, z_prev = p_new, zc[k]
    rho0 = p0 / (th.R_d * T0 * (1.0 + 0.61 * qv0))
    return BaseState(zc=zc, theta0=theta0, qv0=qv0, p0=p0, T0=T0, rho0=rho0, u0=u0, v0=v0)


def example_observed_sounding() -> dict:
    """An illustrative severe-weather proximity sounding in observed (radiosonde) format --
    a moist, unstable boundary layer under drier mid-levels with a veering (clockwise-turning)
    hodograph, the classic supercell environment.  Columns are keyword args for
    :func:`from_observed_sounding`.  ILLUSTRATIVE, not a specific archived observation."""
    return {  # CAPE ~3300 J/kg, 0-6 km shear ~32 m/s, SRH(0-3 km) ~230 m2/s2 -- strong supercell
        "pressure_hPa":   [1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100],
        "height_m":       [100, 780, 1520, 3160, 4420, 5880, 7620, 9640, 10900, 12400, 14200, 16600],
        "temperature_C":  [28.0, 22.0, 16.5, 7.5, 0.0, -9.0, -21.0, -37.0, -47.0, -53.0, -58.0, -63.0],
        "dewpoint_C":     [21.0, 17.5, 12.0, -6.0, -16.0, -27.0, -40.0, -55.0, -63.0, -70.0, -76.0, -80.0],
        "wind_dir_deg":   [160, 180, 200, 230, 245, 255, 260, 265, 270, 270, 270, 270],
        "wind_speed_ms":  [8, 12, 16, 22, 26, 30, 36, 44, 48, 46, 40, 34],
    }


# ---------------------------------------------------------------------------
# hodograph diagnostics (bulk shear + storm-relative helicity)
# ---------------------------------------------------------------------------
def bulk_shear(base: BaseState, z1: float = 0.0, z2: float = 6000.0) -> float:
    """|V(z2) - V(z1)| bulk vertical wind shear [m/s]."""
    z = np.asarray(base.zc)
    u = np.interp([z1, z2], z, np.asarray(base.u0))
    v = np.interp([z1, z2], z, np.asarray(base.v0))
    return float(np.hypot(u[1] - u[0], v[1] - v[0]))


def bulk_richardson_number(base: BaseState, cape: float = None) -> float:
    """Bulk Richardson number BRN = CAPE / (0.5 * U^2), with the BRN shear
    U = |mean(0-6 km) wind - mean(0-500 m) wind| [m/s] (Weisman & Klemp 1982).
    BRN ~ 10-45 favours supercells; large BRN -> multicell, small -> shear-dominated.
    ``cape`` may be supplied to avoid recomputing it."""
    z = np.asarray(base.zc); u = np.asarray(base.u0); v = np.asarray(base.v0)
    lo = z <= 500.0; hi = z <= 6000.0
    if lo.sum() < 1 or hi.sum() < 2:
        return float("nan")
    du = float(np.mean(u[hi]) - np.mean(u[lo])); dv = float(np.mean(v[hi]) - np.mean(v[lo]))
    shear2 = 0.5 * (du * du + dv * dv)
    if cape is None:
        from meteorological_flow.base_state import sounding_diagnostics
        cape = sounding_diagnostics(base)["CAPE_J_kg"]
    return float(cape / shear2) if shear2 > 1e-6 else float("inf")


def _hemisphere_sign(hemisphere=None, latitude_deg=None) -> float:
    """``+1`` for a right-deviating (northern) mover, ``-1`` for a left-deviating (southern).

    Default is ``+1`` -- the northern-hemisphere convention this module was written with --
    so omitting both arguments reproduces the original behaviour exactly.
    """
    if hemisphere is not None:
        h = str(hemisphere).strip().lower()
        if h.startswith("s"):
            return -1.0
        if h.startswith("n"):
            return 1.0
        raise ValueError("hemisphere must be 'north' or 'south', got %r" % (hemisphere,))
    if latitude_deg is not None and float(latitude_deg) < 0.0:
        return -1.0
    return 1.0


def bunkers_storm_motion(base: BaseState, deviation: float = 7.5, hemisphere=None,
                         latitude_deg=None):
    """Deviant supercell motion estimate (Bunkers et al. 2000, simplified).

    Mean 0-6 km wind plus a ``deviation`` [m/s] offset perpendicular to the 0-6 km shear
    vector: to the **RIGHT** in the northern hemisphere (the default, unchanged), to the
    **LEFT** in the southern, where cyclonic rotation is clockwise and the favoured
    supercell is the left-mover.

    Give ``hemisphere="south"`` explicitly, or ``latitude_deg`` to have it inferred
    (negative latitude -> south).  Getting this wrong does not just flip a sign: it
    computes every storm-relative quantity -- SRH, critical angle, streamwise fraction --
    for the *opposite* storm, i.e. for the member that is not the one being studied.
    """
    sign = _hemisphere_sign(hemisphere, latitude_deg)
    z = np.asarray(base.zc)
    sel = z <= 6000.0
    if sel.sum() < 2:
        sel = slice(None)
    um = float(np.mean(np.asarray(base.u0)[sel]))
    vm = float(np.mean(np.asarray(base.v0)[sel]))
    # 0-6 km shear vector
    u = np.interp([0.0, 6000.0], z, np.asarray(base.u0))
    v = np.interp([0.0, 6000.0], z, np.asarray(base.v0))
    shx, shy = u[1] - u[0], v[1] - v[0]
    smag = np.hypot(shx, shy) + 1e-9
    # unit vector 90 deg to the RIGHT of the shear (clockwise): (shy, -shx)/|sh|
    # sign = -1 mirrors it to the LEFT for the southern-hemisphere mover.
    cx = um + sign * deviation * (shy / smag)
    cy = vm + sign * deviation * (-shx / smag)
    return cx, cy


def storm_relative_helicity(base: BaseState, z_top: float = 3000.0,
                            storm_motion=None, z_bot: float = 0.0) -> float:
    """Storm-relative helicity SRH [m^2/s^2] over ``z_bot``..``z_top`` (Davies-Jones 1984).

        SRH = - integral_0^h k . [ (V - C) x dV/dz ]  dz
            =   integral_0^h [ (v-cy) du/dz - (u-cx) dv/dz ] dz

    ``storm_motion`` ``(cx, cy)``; None -> :func:`bunkers_storm_motion` (right-mover).
    Positive SRH (right-moving supercell in the NH) for a clockwise-curving
    (veering) hodograph.

    ``z_bot`` (default 0, i.e. unchanged) opens up *shallow* layers.  That matters here:
    attempt E measured that raising 0-3 km SRH by +67% moved low-level ``V_rot`` by ~0%,
    while the decisive misalignment diagnosed in the vorticity budget sat in the lowest
    few hundred metres -- so ``z_top=500`` is the layer to interrogate, not 3000.
    """
    z = np.asarray(base.zc, dtype=float)
    u = np.asarray(base.u0, dtype=float)
    v = np.asarray(base.v0, dtype=float)
    cx, cy = storm_motion if storm_motion is not None else bunkers_storm_motion(base)
    sel = (z <= z_top) & (z >= z_bot)
    if sel.sum() < 2:
        return 0.0
    zc = z[sel]; uc = u[sel]; vc = v[sel]
    dudz = np.gradient(uc, zc)
    dvdz = np.gradient(vc, zc)
    integrand = (vc - cy) * dudz - (uc - cx) * dvdz
    _trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return float(_trapz(integrand, zc))


# ---------------------------------------------------------------------------
# low-level hodograph GEOMETRY (the quantity the A-L study isolated)
# ---------------------------------------------------------------------------
# The investigation eliminated the *amount* of environmental helicity as the lever
# (attempt E: +67% SRH -> ~0% change in low-level rotation) and localised the
# bottleneck to GEOMETRY: tilting = |omega_h| |grad_h w| cos(theta), with cos(theta)
# measured at only +0.04..+0.09 and NEGATIVE (-0.146) right at the surface, so tilting
# there was generating ANTIcyclonic vorticity.  The free-evolution test (attempt I)
# then showed the streamwise fraction climbing 0.40 -> 0.64 as the storm occluded.
#
# cos(theta) itself needs the 3-D w field (model or dual-Doppler) and has no sounding
# proxy.  The *streamwise fraction* does: for shear-generated horizontal vorticity,
# omega_h is the layer shear rotated 90 deg, so the angle between the shear vector and
# the storm-relative wind -- Esterheld & Giuliano's (2008) "critical angle" -- fixes it
# exactly.  90 deg means the low-level vorticity is purely streamwise.
def _layer_shear_and_sr_wind(base: BaseState, z_top: float, storm_motion=None):
    """``(shear_vector, storm_relative_surface_wind)`` for the 0..``z_top`` layer."""
    z = np.asarray(base.zc, dtype=float)
    u = np.asarray(base.u0, dtype=float)
    v = np.asarray(base.v0, dtype=float)
    cx, cy = storm_motion if storm_motion is not None else bunkers_storm_motion(base)
    u_lo, u_hi = np.interp([0.0, z_top], z, u)
    v_lo, v_hi = np.interp([0.0, z_top], z, v)
    shear = np.array([u_hi - u_lo, v_hi - v_lo])
    sr_wind = np.array([u_lo - cx, v_lo - cy])
    return shear, sr_wind


#: Below this layer shear [m/s] the shear *direction* is numerically degenerate, so the
#: critical angle carries no information at all and NaN is returned.
MIN_SHEAR_MS = 0.5

#: Below this layer shear [m/s] the angle is computable but ILL-CONDITIONED: a 1-2 m/s
#: shear vector's direction swings wildly on tiny perturbations, so a confident-looking
#: angle would be noise.  Esterheld & Giuliano's critical angle was established on
#: supercell environments with substantial low-level shear, not on near-calm hodographs.
WELL_POSED_SHEAR_MS = 2.0

#: Same idea for the storm-relative inflow: a near-zero one makes "streamwise" undefined.
WELL_POSED_SR_WIND_MS = 2.0


def critical_angle(base: BaseState, z_top: float = 500.0, storm_motion=None,
                   min_shear_ms: float = MIN_SHEAR_MS) -> float:
    """Critical angle [deg] -- angle between the 0..``z_top`` shear vector and the
    storm-relative *surface* wind (Esterheld & Giuliano 2008).

    ~90 deg means the storm-relative inflow is parallel to the low-level horizontal
    vorticity, i.e. that vorticity is purely **streamwise** and is tilted straight into
    the updraft.  Departures from 90 deg leave a crosswise component that tilting cannot
    convert into cyclonic vertical vorticity.

    Returns NaN when the hodograph is degenerate -- layer shear below ``min_shear_ms`` or a
    vanishing storm-relative wind.  **A finite return is not a promise the number is
    meaningful**: between ``min_shear_ms`` and :data:`WELL_POSED_SHEAR_MS` the angle is
    ill-conditioned.  :func:`low_level_geometry_report` flags that case explicitly; prefer
    it over calling this directly.
    """
    shear, sr_wind = _layer_shear_and_sr_wind(base, z_top, storm_motion)
    ns, nw = np.hypot(*shear), np.hypot(*sr_wind)
    if ns < min_shear_ms or nw < 1e-9:
        return float("nan")
    cos = float(np.dot(shear, sr_wind) / (ns * nw))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def streamwise_vorticity_fraction(base: BaseState, z_top: float = 500.0,
                                  storm_motion=None,
                                  min_shear_ms: float = MIN_SHEAR_MS) -> float:
    """Fraction of the 0..``z_top`` horizontal vorticity that is **streamwise**, in [0, 1].

    Shear-generated horizontal vorticity is ``omega_h = (-dv/dz, du/dz)`` -- the shear
    vector rotated 90 deg -- so this is ``|cos(angle(omega_h, V_storm_relative))|``, which
    reduces exactly to ``|sin(critical_angle)|``.  It is the sounding-side counterpart of
    the streamwise fraction that :func:`vorticity_budget.tilting_efficiency` measures on
    model fields (0.40 -> 0.64 as the attempt-I supercell occluded); 1.0 is purely
    streamwise, 0.0 purely crosswise.

    Same conditioning caveat as :func:`critical_angle` -- NaN when degenerate, and a finite
    value between ``min_shear_ms`` and :data:`WELL_POSED_SHEAR_MS` is noise.
    """
    shear, sr_wind = _layer_shear_and_sr_wind(base, z_top, storm_motion)
    ns, nw = np.hypot(*shear), np.hypot(*sr_wind)
    if ns < min_shear_ms or nw < 1e-9:
        return float("nan")
    omega_h = np.array([-shear[1], shear[0]])          # shear rotated +90 deg
    return float(abs(np.dot(omega_h, sr_wind)) / (ns * nw))


def low_level_geometry_report(base: BaseState, z_top: float = 500.0, storm_motion=None,
                              hemisphere=None, latitude_deg=None):
    """Bundle the low-level geometry discriminators for one sounding.

    ``hemisphere``/``latitude_deg`` select the deviant mover (see
    :func:`bunkers_storm_motion`); both omitted keeps the northern-hemisphere default.

    SRH is reported twice on purpose.  ``SRH_*_m2_s2`` is the raw Davies-Jones integral,
    which for the favoured southern-hemisphere (left-moving) supercell comes out
    **negative** -- GFS reports it that way too.  ``SRH_*_cyclonic_m2_s2`` flips that sign
    in the southern hemisphere so the number is directly comparable with the
    northern-hemisphere values quoted in the literature and elsewhere in this project.
    Compare magnitudes across hemispheres, never raw signs.
    """
    import warnings

    sign = _hemisphere_sign(hemisphere, latitude_deg)
    cx, cy = (storm_motion if storm_motion is not None
              else bunkers_storm_motion(base, hemisphere=hemisphere, latitude_deg=latitude_deg))
    srh_layer = storm_relative_helicity(base, z_top=z_top, storm_motion=(cx, cy))
    srh_3km = storm_relative_helicity(base, z_top=3000.0, storm_motion=(cx, cy))

    shear_vec, sr_wind = _layer_shear_and_sr_wind(base, z_top, (cx, cy))
    shear_mag = float(np.hypot(*shear_vec))
    sr_mag = float(np.hypot(*sr_wind))
    reasons = []
    if shear_mag < WELL_POSED_SHEAR_MS:
        reasons.append("0-%.0f m shear is only %.2f m/s (< %.1f)"
                       % (z_top, shear_mag, WELL_POSED_SHEAR_MS))
    if sr_mag < WELL_POSED_SR_WIND_MS:
        reasons.append("storm-relative inflow is only %.2f m/s (< %.1f)"
                       % (sr_mag, WELL_POSED_SR_WIND_MS))
    if reasons:
        warnings.warn(
            "low-level geometry is ILL-CONDITIONED: %s. The critical angle and streamwise "
            "fraction are reported but carry no information at this shear -- the shear "
            "direction swings on noise. They have NOT been silently presented as valid; "
            "check 'geometry_is_ill_conditioned' before quoting them."
            % "; ".join(reasons), RuntimeWarning, stacklevel=2)

    return {
        "layer_top_m": float(z_top),
        "hemisphere": "south" if sign < 0 else "north",
        "storm_motion_ms": (float(cx), float(cy)),
        "critical_angle_deg": critical_angle(base, z_top, (cx, cy)),
        "streamwise_fraction": streamwise_vorticity_fraction(base, z_top, (cx, cy)),
        "SRH_layer_m2_s2": srh_layer,
        "SRH_0_3km_m2_s2": srh_3km,
        "SRH_layer_cyclonic_m2_s2": sign * srh_layer,
        "SRH_0_3km_cyclonic_m2_s2": sign * srh_3km,
        "shear_layer_m_s": bulk_shear(base, 0.0, z_top),
        "shear_0_6km_m_s": bulk_shear(base, 0.0, 6000.0),
        "storm_relative_inflow_m_s": sr_mag,
        "geometry_is_ill_conditioned": bool(reasons),
        "ill_conditioned_reason": "; ".join(reasons) if reasons else "",
    }


__all__ = [
    "build_sounding", "from_observed_sounding", "example_observed_sounding",
    "bulk_shear", "bunkers_storm_motion", "storm_relative_helicity",
    "bulk_richardson_number", "critical_angle", "streamwise_vorticity_fraction",
    "low_level_geometry_report",
]
