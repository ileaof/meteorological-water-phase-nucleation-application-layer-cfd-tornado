"""Objective, resolution-aware storm classification.

Maps a bundle of diagnostics (rotation report + vortex report + cold-pool report) onto a discrete
category ladder -- WITHOUT ever imposing a vortex; the class is *read* from the resolved fields:

    NO_DEEP_CONVECTION
    ORDINARY_CONVECTION
    SUPERCELL
    LOW_LEVEL_MESOCYCLONE
    TORNADO_LIKE_VORTEX
    SURFACE_CONNECTED_TORNADO_LIKE_VORTEX

The decisive thresholds are on *resolution-robust* quantities -- peak tangential velocity,
circulation, pressure deficit, updraft helicity -- rather than raw zeta (which grows as dx shrinks),
so a coarse and a fine run of the same storm classify consistently.  "TORNADO_LIKE_VORTEX" is
deliberate: at these grid spacings we resolve a tornado-*scale* circulation, not the true core.
Every threshold is overridable.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

CATEGORIES = [
    "NO_DEEP_CONVECTION", "ORDINARY_CONVECTION", "SUPERCELL",
    "LOW_LEVEL_MESOCYCLONE", "TORNADO_LIKE_VORTEX", "SURFACE_CONNECTED_TORNADO_LIKE_VORTEX",
]


@dataclass
class ClassThresholds:
    w_deep_m_s: float = 10.0                 # deep-convection updraft
    meso_zeta_s: float = 5.0e-3              # mid-level mesocyclone vertical vorticity
    meso_uh_m2_s2: float = 50.0             # mid-level updraft helicity
    lowlevel_zeta_s: float = 3.0e-3         # near-surface vertical vorticity
    lowlevel_vtheta_m_s: float = 8.0        # near-surface tangential velocity
    tlv_vtheta_m_s: float = 15.0            # tornado-like tangential velocity
    tlv_circulation_m2_s: float = 4.0e4     # tornado-like circulation
    tlv_pdeficit_Pa: float = -200.0         # tornado-like pressure deficit (<=)
    surface_level_m: float = 250.0          # "surface-connected" if the vortex reaches this low
    surface_convergence_s: float = 5.0e-3   # low-level convergence at the vortex
    min_persistence_s: float = 120.0        # required rotation lifetime for the tornado-like tiers


def classify(diag: dict, thresholds: ClassThresholds = None, persistence_s: float = None,
             allow_pressure: bool = True) -> dict:
    """Classify a storm from a merged diagnostics dict. Recognised keys (all optional, default 0):
    ``w_max``, ``midlevel_mesocyclone``, ``updraft_helicity_2_5km``, ``near_surface_zeta_max``,
    ``v_theta_max_m_s``, ``circulation_m2_s``, ``pressure_deficit_Pa``, vortex ``level_m``,
    ``gust_front_convergence_s``.  ``persistence_s`` (from a tracker) gates the tornado-like tiers.
    Returns ``{category, criteria, rank}``."""
    th = thresholds or ClassThresholds()
    g = lambda k, d=0.0: float(diag.get(k, d) if diag.get(k, d) is not None else d)

    w = g("w_max")
    meso = (g("midlevel_mesocyclone") >= th.meso_zeta_s) and (g("updraft_helicity_2_5km") >= th.meso_uh_m2_s2)
    lowlevel = (g("near_surface_zeta_max") >= th.lowlevel_zeta_s) or (g("v_theta_max_m_s") >= th.lowlevel_vtheta_m_s)
    # PRESSURE BRANCH GATING (REVIEW_REQUEST.md A8).  On a NEST the projection absorbs the
    # imposed-inflow imbalance, so `phi` (hence p_dyn and any pressure deficit) carries a large
    # boundary-driven component ~1/dt: measured local values on one 67 m field ran -1418, -3250,
    # -5161 and +7412 Pa against a cyclostrophic scale of ~-120 Pa.  A -200 Pa threshold is
    # meaningless against that noise, so the pressure branch must NOT be able to promote a nest
    # to a tornado-like tier.  `allow_pressure=False` disables that branch; the v_theta branch is
    # unaffected.  Historically the dP branch was dead anyway (A-K failed the v_theta branch
    # honestly at ~7.8 m/s), so gating it changes no past verdict -- it prevents a FUTURE
    # false positive.
    _p_ok = bool(allow_pressure) and diag.get("pressure_deficit_Pa") is not None
    tlv_pressure = (_p_ok and g("circulation_m2_s") >= th.tlv_circulation_m2_s
                    and g("pressure_deficit_Pa", 0.0) <= th.tlv_pdeficit_Pa)
    tlv = (g("v_theta_max_m_s") >= th.tlv_vtheta_m_s) or tlv_pressure
    persist_ok = (persistence_s is None) or (persistence_s >= th.min_persistence_s)
    surface = tlv and (g("level_m", 1e9) <= th.surface_level_m) \
        and (g("gust_front_convergence_s") >= th.surface_convergence_s) and persist_ok

    criteria = {"deep_convection": w >= th.w_deep_m_s, "mesocyclone": meso,
                "low_level_rotation": lowlevel, "tornado_like": tlv and persist_ok,
                "surface_connected": surface,
                "pressure_branch_available": _p_ok, "tornado_like_via_pressure": tlv_pressure}

    if w < th.w_deep_m_s:
        cat = "NO_DEEP_CONVECTION"
    elif not meso:
        cat = "ORDINARY_CONVECTION"
    elif not lowlevel:
        cat = "SUPERCELL"
    elif not (tlv and persist_ok):
        cat = "LOW_LEVEL_MESOCYCLONE"
    elif not surface:
        cat = "TORNADO_LIKE_VORTEX"
    else:
        cat = "SURFACE_CONNECTED_TORNADO_LIKE_VORTEX"
    return {"category": cat, "rank": CATEGORIES.index(cat), "criteria": criteria,
            "thresholds": asdict(th)}


def classify_simulation(sim, z_surface_m=150.0, storm_motion=(0.0, 0.0),
                        thresholds: ClassThresholds = None, persistence_s: float = None,
                        allow_pressure: bool = None) -> dict:
    """Assemble the diagnostics from a live simulation and classify it. Combines
    ``rotation.rotation_report``, ``vortex_diagnostics.vortex_report`` (near-surface) and
    ``coldpool.coldpool_report``.

    ``allow_pressure=None`` (default) AUTO-DISABLES the pressure branch of the tornado-like test
    on a NEST, where the projection absorbs the imposed-inflow imbalance and the pressure deficit
    is not trustworthy (REVIEW_REQUEST.md A8).  Pass True to force it on (not recommended on a
    nest) or False to force it off."""
    from . import rotation as rot
    from . import vortex_diagnostics as vd
    from . import coldpool as cp
    is_nest = hasattr(sim, "spec")                    # NestedStormSimulation carries its NestSpec
    if allow_pressure is None:
        allow_pressure = not is_nest
    diag = dict(rot.rotation_report(sim.state, sim.grid, base=getattr(sim, "base", None)))
    vr = vd.vortex_report(sim.state, sim.grid, z_m=z_surface_m, storm_motion=storm_motion)
    diag.update(vr)                                   # adds v_theta_max_m_s, circulation, level_m, ...
    diag.update(cp.coldpool_report(sim.state, sim.grid, z_m=z_surface_m))
    out = classify(diag, thresholds=thresholds, persistence_s=persistence_s,
                   allow_pressure=allow_pressure)
    out["diagnostics"] = diag
    out["pressure_branch"] = {
        "allowed": bool(allow_pressure), "is_nest": bool(is_nest),
        "reason": ("nest pressure is boundary-contaminated (A8): projection absorbs the imposed "
                   "inflow imbalance, deficit ~1/dt" if is_nest and not allow_pressure
                   else "parent domain: pressure branch enabled" if allow_pressure
                   else "explicitly disabled by caller")}
    return out


__all__ = ["CATEGORIES", "ClassThresholds", "classify", "classify_simulation"]
