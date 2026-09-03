"""Configuration for the rotating storm_dynamics core.

A :class:`StormConfig` carries the *reused* :class:`meteorological_flow.config.
SimulationConfig` (``.sim`` -- passed straight to the reused grid / projection /
buoyancy / transport / microphysics helpers, which already expect that object)
plus a :class:`StormDynamicsConfig` (``.dyn``) holding the knobs that are new to
this fork: the Coriolis parameter, the LES closure, the surface bulk-drag law,
and the curved-hodograph environmental shear.

The demonstration core's Rayleigh drag (``flow.gamma_damp``) and the hard 120 m/s
velocity clip are *not* used here -- the LES closure provides the physically
motivated dissipation instead.  ``flow.gamma_damp`` is forced to 0 by
:func:`build_storm_config` so a stray value cannot re-damp the vortex; only a
documented extreme numerical guard (``dyn.v_guard``) remains.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from meteorological_flow.config import SimulationConfig, apply_overrides


# Earth rotation rate [rad/s] (sidereal); f = 2 Omega sin(lat).
_OMEGA_EARTH = 7.2921159e-5


def coriolis_f(latitude_deg: float) -> float:
    """f-plane Coriolis parameter [1/s] at a latitude.  Tornado Alley ~ 35-40N."""
    return 2.0 * _OMEGA_EARTH * math.sin(math.radians(latitude_deg))


@dataclass
class HodographConfig:
    """Environmental wind profile (curved hodograph -> storm-relative helicity).

    ``kind``:
      * ``"unidirectional"`` -- straight hodograph (Weisman-Klemp ``u_shear``);
        the M1 storm-splitting reference (no preferred left/right member).
      * ``"quarter_circle"`` -- a Klemp/Rotunno quarter-circle turning of the wind
        with height through the ``z_turn`` layer, then straight above; the curved
        low-level hodograph that gives streamwise (storm-relative) vorticity and
        favours the right-moving supercell + low-level rotation (M2).
    """
    kind: str = "quarter_circle"
    U_max: float = 25.0          # asymptotic hodograph speed [m/s]
    z_turn: float = 3000.0       # depth of the curved (turning) layer [m]
    u_half: float = 3000.0       # tanh ramp half-depth for the straight part [m]
    storm_motion: tuple = None   # (cx, cy) [m/s]; None -> estimated from the profile


@dataclass
class LESConfig:
    """Subgrid-scale closure (replaces the demonstration Rayleigh drag + clip)."""
    model: str = "smagorinsky"   # smagorinsky (implemented) | none (background only)
                                 # | tke15 (Deardorff, documented future work)
    C_s: float = 0.2             # Smagorinsky constant (0.18-0.23 typical)
    Pr_t: float = 1.0 / 3.0      # turbulent Prandtl number (K_h = K_m / Pr_t)
    Ri_c: float = 0.25           # critical Richardson number (stability damping)
    nu_background: float = 1.0   # small background viscosity floor [m^2/s]


@dataclass
class SurfaceDragConfig:
    """Bottom bulk aerodynamic drag law (corner-flow friction; enables low-level
    rotation).  tau = rho C_d |V| V applied to the lowest model level."""
    enabled: bool = True
    C_d: float = 0.012           # bulk drag coefficient (~0.01-0.02 over land)
    U_min: float = 1.0           # floor on |V| so weak flow still feels some drag [m/s]
    # Height-consistent (neutral log-law / MOST) closure: C_d = (kappa/ln(z1/z0))^2 evaluated at the
    # ACTUAL first cell-centre height z1.  A fixed bulk C_d is only calibrated for one z1 -- at a
    # refined near-surface mesh (z1 ~ 5-10 m, the corner-flow layer) it over-damps the lowest level.
    # Off by default (fixed C_d, byte-identical); on -> C_d follows the mesh and the roughness.
    use_log_law: bool = False
    roughness_length_m: float = 0.1   # aerodynamic roughness z0 [m] (~0.1 over open land)
    kappa: float = 0.4                # von Karman constant
    # Surface-layer STRESS-DIVERGENCE form: instead of damping only the lowest cell (sink rate
    # C_d|V|/dz1, which BLOWS UP as the mesh is refined and strips the tangential wind), apply the
    # stress as a flux divergence spread over a PHYSICAL surface-layer depth: tau(z)=tau_s(1-z/h)
    # => du/dt = -tau_s/h, uniform through the layer and resolution-independent.  This retains the
    # near-surface tangential wind while still producing the drag-induced imbalance that drives the
    # convergent corner flow.  Off by default (lowest-cell damping, byte-identical).
    stress_divergence: bool = False
    surface_layer_depth_m: float = 150.0


@dataclass
class SurfaceFluxConfig:
    """Bulk surface sensible + latent heat fluxes (opt-in; complements the momentum drag).
    Relaxes the lowest model level toward a surface potential temperature and humidity via bulk
    transfer coefficients ``H = rho cp C_h |V| (theta_sfc - theta)``, ``E = rho Lv C_q |V|
    (q_sfc - q)`` -- the boundary-layer heat/moisture source that deepens the mixed layer and can
    modulate cold-pool recovery.  Off by default (no behaviour change)."""
    enabled: bool = False
    C_h: float = 1.2e-3           # sensible-heat bulk transfer coefficient
    C_q: float = 1.2e-3           # latent-heat (moisture) bulk transfer coefficient
    U_min: float = 1.0            # gustiness floor on |V| [m/s]
    theta_sfc_K: float = None     # surface potential temperature; None -> base surface theta0 + dtheta
    dtheta_sfc_K: float = 0.0     # warm-ground excess over the environmental surface theta [K]
    saturate_surface: bool = True # drive q_v toward q_sat at the surface (moist ground)
    qv_sfc_kgkg: float = None     # explicit surface q_v if not saturating
    roughness_length_m: float = 0.1


@dataclass
class MesoForcingConfig:
    """Sustained low-level mesoscale-ascent forcing (a dryline/convergence proxy):
    a heating(+moistening) cylinder that keeps lifting parcels through the CIN cap for
    ``duration_s``, so a supercell can establish from a *real* (capped) environment
    instead of a single warm bubble that pulses and decays.  Off by default -- idealised
    runs keep the one-shot bubble (see :mod:`storm_dynamics.forcing`)."""
    enabled: bool = False
    heat_rate_K_s: float = 0.0        # potential-temperature source in the core [K/s]
    moist_rate_kgkg_s: float = 0.0    # water-vapour source in the core [kg/kg/s]
    radius_m: float = 6000.0          # horizontal radius of the forced cylinder
    z_top_m: float = 2500.0           # forcing confined below this height [m]
    duration_s: float = 1200.0        # active for the first duration_s, then removed
    center: tuple = None              # (x, y) [m]; None -> domain centre


@dataclass
class StormDynamicsConfig:
    latitude_deg: float = 36.0                       # Tornado Alley (f-plane)
    coriolis: bool = True
    f: float | None = None                            # explicit f [1/s]; None -> from latitude
    momentum_advection: bool = True                   # the enabling item (1)
    momentum_order: int = 2                            # 1 upwind | 2 MUSCL/minmod
    les: LESConfig = field(default_factory=LESConfig)
    drag: SurfaceDragConfig = field(default_factory=SurfaceDragConfig)
    fluxes: SurfaceFluxConfig = field(default_factory=SurfaceFluxConfig)
    forcing: MesoForcingConfig = field(default_factory=MesoForcingConfig)
    hodograph: HodographConfig = field(default_factory=HodographConfig)
    v_guard: float = 150.0     # extreme numerical guard only (NOT a physical cap);
                               # documented, should never bite in a resolved vortex.

    def f_value(self) -> float:
        if not self.coriolis:
            return 0.0
        return self.f if self.f is not None else coriolis_f(self.latitude_deg)


@dataclass
class StormConfig:
    sim: SimulationConfig = field(default_factory=SimulationConfig)
    dyn: StormDynamicsConfig = field(default_factory=StormDynamicsConfig)


def build_storm_config(*, preset: str = "storm", nx=None, ny=None, nz=None,
                       Lx=None, Ly=None, Lz=None, duration=None,
                       dynamics: str = "anelastic", periodic: bool = True,
                       hodograph_kind: str = "quarter_circle",
                       latitude_deg: float = 36.0,
                       coriolis: bool = True,
                       drag: bool = True,
                       les_model: str = "smagorinsky",
                       z_stretch=None, U_max=None, z_turn=None, C_s=None,
                       couple_nucleation: bool = False, nucleation_method="lookup",
                       device: str = "cpu",
                       cfl=None, dt_max=None, outdir=None) -> StormConfig:
    """Assemble a :class:`StormConfig` for the rotating storm.

    Starts from the repo's deep-convection storm preset (stratified base state,
    periodic lateral BCs for the mean-wind shear, anelastic core, two-way
    microphysics) then swaps in the storm_dynamics physics: no Rayleigh drag /
    velocity clip (LES instead), Coriolis on, surface drag on.
    """
    sim = apply_overrides(SimulationConfig(), preset=preset, periodic=periodic,
                          dynamics=dynamics)
    # the fork provides its own dissipation (LES) and momentum transport, so the
    # demonstration Rayleigh drag must be OFF -- otherwise it relaxes the vortex.
    sim.flow.gamma_damp = 0.0
    # 2nd-order MUSCL scalar transport keeps the sounding / cold-pool gradients sharp.
    sim.flow.advection_order = 2
    # optional: feed the validated nucleation kernel rate J as the microphysics
    # embryo source (eq39 pathway), exactly as meteorological_flow does.  Off by
    # default (building the lookup table is slow; the CCN/IN activation already
    # drives the cold pool).
    sim.nucleation.couple_kernel = bool(couple_nucleation)
    sim.nucleation.method = str(nucleation_method)
    # compute backend: "cpu" (default, preserves prior behaviour), "gpu" (fails
    # loudly if CuPy/CUDA absent), or "auto" (GPU when available, else CPU).
    sim.performance.device = str(device)
    if nx is not None:
        sim.grid.nx = int(nx)
    if ny is not None:
        sim.grid.ny = int(ny)
    if nz is not None:
        sim.grid.nz = int(nz)
    if Lx is not None:
        sim.domain.Lx = float(Lx)
    if Ly is not None:
        sim.domain.Ly = float(Ly)
    if Lz is not None:
        sim.domain.Lz = float(Lz)
    if duration is not None:
        sim.time.duration = float(duration)
    if z_stretch is not None:
        sim.grid.z_stretch = float(z_stretch)
    if cfl is not None:
        sim.time.cfl = float(cfl)
    if dt_max is not None:
        sim.time.dt_max = float(dt_max)
    if outdir is not None:
        sim.output.outdir = outdir

    les = LESConfig(model=les_model)
    if C_s is not None:
        les.C_s = float(C_s)
    hodo = HodographConfig(kind=hodograph_kind)
    if U_max is not None:
        hodo.U_max = float(U_max)
    if z_turn is not None:
        hodo.z_turn = float(z_turn)
    dyn = StormDynamicsConfig(
        latitude_deg=latitude_deg,
        coriolis=coriolis,
        les=les,
        drag=SurfaceDragConfig(enabled=drag),
        hodograph=hodo,
    )
    return StormConfig(sim=sim, dyn=dyn)


def storm_config_from_dict(d: dict) -> StormConfig:
    """Build a :class:`StormConfig` from a plain (YAML) dict.

    Recognised keys (all optional; defaults mirror :func:`build_storm_config`)::

        preset, nx, ny, nz, Lx, Ly, Lz, duration, dt_max, cfl, dynamics,
        periodic, latitude_deg, coriolis, drag, les_model, hodograph_kind, outdir
    """
    d = d or {}
    return build_storm_config(
        preset=str(d.get("preset", "storm")),
        nx=d.get("nx"), ny=d.get("ny"), nz=d.get("nz"),
        Lx=d.get("Lx"), Ly=d.get("Ly"), Lz=d.get("Lz"),
        duration=d.get("duration"), dt_max=d.get("dt_max"), cfl=d.get("cfl"),
        dynamics=str(d.get("dynamics", "anelastic")),
        periodic=bool(d.get("periodic", True)),
        hodograph_kind=str(d.get("hodograph_kind", "quarter_circle")),
        latitude_deg=float(d.get("latitude_deg", 36.0)),
        coriolis=bool(d.get("coriolis", True)),
        drag=bool(d.get("drag", True)),
        les_model=str(d.get("les_model", "smagorinsky")),
        z_stretch=d.get("z_stretch"), U_max=d.get("U_max"),
        z_turn=d.get("z_turn"), C_s=d.get("C_s"),
        device=str(d.get("device", "cpu")),
        outdir=d.get("outdir"),
    )


def storm_config_from_yaml(path: str) -> StormConfig:
    import yaml
    with open(path, "r", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    return storm_config_from_dict(raw if isinstance(raw, dict) else {})


__all__ = [
    "HodographConfig", "LESConfig", "SurfaceDragConfig", "SurfaceFluxConfig", "MesoForcingConfig",
    "StormDynamicsConfig", "StormConfig",
    "build_storm_config", "storm_config_from_dict", "storm_config_from_yaml",
    "coriolis_f",
]
