"""Unified `tornadogenesis:` configuration overlay (additive).

The repo has two working config schemas (the idealised ``storm_dynamics`` flat YAML and the
real-case ``atmospheric_data`` structured YAML).  This overlay adds the spec's block-structured
schema WITHOUT replacing either: it parses ``tornadogenesis / environment / convection / surface /
microphysics / diagnostics / nesting / domain`` and maps them onto the existing
:class:`StormConfig` dataclasses (via :func:`build_storm_config`) plus the diagnostic/nesting knobs.

The invariant ``tornadogenesis.impose_vortex`` must be ``false`` -- the model never inserts a
Rankine/Lamb-Oseen vortex; rotation emerges from the equations.  A ``true`` value raises.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

from .config import build_storm_config, StormConfig, SurfaceFluxConfig, MesoForcingConfig


@dataclass
class TornadoRunConfig:
    storm: StormConfig
    environment: dict = field(default_factory=dict)
    convection: dict = field(default_factory=dict)
    diagnostics: dict = field(default_factory=dict)
    nesting: dict = field(default_factory=dict)
    impose_vortex: bool = False       # enforced False


_DEFAULT_DIAGNOSTICS = {
    "vorticity_budget": True, "updraft_helicity": True, "circulation": True,
    "pressure_deficit": True, "cold_pool": True, "radar_operator": True,
    "classification": True, "macro_micro_gradient": True,
}


def load_tornadogenesis_config(source) -> TornadoRunConfig:
    """Build a :class:`TornadoRunConfig` from a spec-shaped dict or YAML path."""
    if isinstance(source, str):
        import yaml
        with open(source) as f:
            d = yaml.safe_load(f) or {}
    else:
        d = dict(source or {})

    tg = d.get("tornadogenesis", {}) or {}
    if tg.get("impose_vortex", False):
        raise ValueError("tornadogenesis.impose_vortex must be false -- the tornado is never "
                         "imposed; rotation must emerge from the governing equations.")
    env = d.get("environment", {}) or {}
    conv = d.get("convection", {}) or {}
    surf = d.get("surface", {}) or {}
    micro = d.get("microphysics", {}) or {}
    diag = dict(_DEFAULT_DIAGNOSTICS); diag.update(d.get("diagnostics", {}) or {})
    nest = d.get("nesting", {}) or {}
    dom = d.get("domain", {}) or {}

    # domain: explicit nx/Lx, or derive Lx from parent_dx_m * nx
    nx = int(dom.get("nx", 64)); ny = int(dom.get("ny", nx)); nz = int(dom.get("nz", 44))
    parent_dx = float(nest.get("parent_dx_m", dom.get("dx_m", 1000.0)))
    Lx = float(dom.get("Lx_m", parent_dx * nx)); Ly = float(dom.get("Ly_m", parent_dx * ny))
    Lz = float(dom.get("Lz_m", 15000.0))

    scfg = build_storm_config(
        preset="storm", nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz,
        duration=float(d.get("duration_s", 1.0)), dt_max=d.get("dt_max_s", None),
        z_stretch=dom.get("z_stretch", 1.05), device=d.get("device", "cpu"),
        drag=bool(surf.get("drag_enabled", True)),
        les_model=str(micro.get("les_model", d.get("les_model", "smagorinsky"))),
        coriolis=bool(env.get("coriolis", True)),
        couple_nucleation=bool(micro.get("shifted_equilibrium_nucleation", False)),
    )
    # convection trigger (idealised warm bubble; NO rotation)
    if str(conv.get("initiation", "warm_bubble")) == "warm_bubble":
        scfg.sim.physics.bubble_dtheta = float(conv.get("bubble_theta_amplitude_K", 2.0))
    # sustained convergence forcing (optional)
    if str(conv.get("initiation", "")) == "surface_convergence" or conv.get("sustained_forcing"):
        scfg.dyn.forcing = MesoForcingConfig(
            enabled=True, heat_rate_K_s=float(conv.get("forcing_heat_K_s", 0.006)),
            moist_rate_kgkg_s=float(conv.get("forcing_moist_kgkg_s", 3e-6)),
            radius_m=float(conv.get("forcing_radius_m", 7000.0)),
            duration_s=float(conv.get("forcing_duration_s", 1500.0)))
    # surface heat/moisture fluxes (optional)
    if surf.get("sensible_heat_flux") or surf.get("latent_heat_flux"):
        scfg.dyn.fluxes = SurfaceFluxConfig(
            enabled=True,
            C_h=float(surf.get("C_h", 1.2e-3)) if surf.get("sensible_heat_flux") else 0.0,
            C_q=float(surf.get("C_q", 1.2e-3)) if surf.get("latent_heat_flux") else 0.0,
            saturate_surface=bool(surf.get("latent_heat_flux", False)),
            dtheta_sfc_K=float(surf.get("surface_theta_excess_K", 0.0)),
            roughness_length_m=float(surf.get("roughness_length_m", 0.1)))
    if "roughness_length_m" in surf and scfg.dyn.drag.enabled:
        pass  # z0 is carried on the flux config; drag keeps its bulk C_d

    return TornadoRunConfig(storm=scfg, environment=env, convection=conv, diagnostics=diag,
                            nesting=nest, impose_vortex=False)


def run_diagnostics(sim, cfg: TornadoRunConfig, z_surface_m=150.0, storm_motion=(0.0, 0.0)) -> dict:
    """Run the diagnostics enabled in ``cfg.diagnostics`` on a live simulation and return one dict
    (the spec's diagnostic bundle: rotation, vorticity budget, vortex, cold pool, classification,
    macro/micro gradient)."""
    from . import rotation as rot, vorticity_budget as vb, vortex_diagnostics as vd
    from . import coldpool as cp, classification as cl, micro_gradient as mg
    d = cfg.diagnostics
    out = {}
    sim.state.diagnose(sim.cfg)
    if d.get("updraft_helicity", True) or True:
        out["rotation"] = rot.rotation_report(sim.state, sim.grid, base=getattr(sim, "base", None))
    if d.get("vorticity_budget", True):
        terms = vb.zeta_budget(sim.state, sim.grid, Km=getattr(sim, "_Km", None))
        out["vorticity_budget_low"] = vb.budget_layer_summary(terms, sim.grid, 0.0, 1000.0)
        out["dominant_low_mechanism"] = vb.dominant_mechanism(terms, sim.grid, 0.0, 1000.0)[0]
    if d.get("circulation", True) or d.get("pressure_deficit", True):
        out["vortex"] = vd.vortex_report(sim.state, sim.grid, z_m=z_surface_m, storm_motion=storm_motion)
    if d.get("cold_pool", True):
        out["cold_pool"] = cp.coldpool_report(sim.state, sim.grid, z_m=z_surface_m)
    if d.get("macro_micro_gradient", True):
        out["micro_gradient"] = mg.nucleation_diagnostics(sim.state, sim.grid)
    if d.get("classification", True):
        out["classification"] = cl.classify_simulation(sim, z_surface_m=z_surface_m,
                                                       storm_motion=storm_motion)["category"]
    return out


__all__ = ["TornadoRunConfig", "load_tornadogenesis_config", "run_diagnostics"]
