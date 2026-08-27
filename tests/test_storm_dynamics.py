"""storm_dynamics -- unit + regression tests for the rotating core.

Fast unit tests cover each new dynamical piece in isolation (momentum
conservation, Coriolis, LES viscosity, surface drag, curved-hodograph SRH); the
short-integration tests check that the assembled :class:`StormSimulation` keeps
water / mass conservation and (M1) spins up rotation.  The heavier full-milestone
runs (storm splitting; low-level rotation) live in
``test_storm_milestones.py`` so the default suite stays quick.
"""
from __future__ import annotations

import numpy as np
import pytest

from meteorological_flow.grid import Grid
from meteorological_flow.state import FlowState

from storm_dynamics import soundings as snd
from storm_dynamics.config import (
    HodographConfig, LESConfig, StormDynamicsConfig, SurfaceDragConfig,
    build_storm_config, coriolis_f,
)
from storm_dynamics.coriolis import add_coriolis
from storm_dynamics.momentum import momentum_advection_tendency
from storm_dynamics.surface_drag import apply_surface_drag
from storm_dynamics.turbulence import strain_and_viscosity


def _periodic_grid(nx=12, ny=10, nz=8):
    return Grid(nx=nx, ny=ny, nz=nz, Lx=1200, Ly=1000, Lz=800, periodic=True)


# --------------------------------------------------------------------------
# item 1: conservative flux-form momentum advection
# --------------------------------------------------------------------------
def test_momentum_uniform_flow_zero_tendency():
    """A spatially uniform horizontal flow has no advective tendency."""
    g = _periodic_grid()
    st = FlowState.zeros(g)
    st.u[:] = 7.3; st.v[:] = -2.1
    du, dv, dw = momentum_advection_tendency(st, g, order=2, periodic=True)
    assert float(np.abs(du).max()) < 1e-12
    assert float(np.abs(dv).max()) < 1e-12
    assert float(np.abs(dw).max()) < 1e-12


def test_momentum_globally_conservative():
    """Flux form telescopes: the domain-integrated momentum tendency is ~0
    (periodic laterals, solid z-walls) -- momentum is conserved, no clip needed."""
    g = _periodic_grid()
    st = FlowState.zeros(g)
    xc = g.xc.reshape(-1, 1, 1); yc = g.yc.reshape(1, -1, 1)
    st.u[:] = np.sin(2 * np.pi * np.linspace(0, 1, g.nx + 1)).reshape(-1, 1, 1) \
        * np.cos(2 * np.pi * yc / g.Ly)
    st.v[:] = np.cos(2 * np.pi * np.linspace(0, 1, g.ny + 1)).reshape(1, -1, 1) \
        * np.sin(2 * np.pi * xc / g.Lx)
    du, dv, dw = momentum_advection_tendency(st, g, order=1, periodic=True)
    # independent faces only (drop the duplicated periodic face)
    assert abs(float(du[:-1].sum())) < 1e-10
    assert abs(float(dv[:, :-1].sum())) < 1e-10
    assert abs(float(dw.sum())) < 1e-10


def test_momentum_produces_finite_tendency():
    g = _periodic_grid()
    st = FlowState.zeros(g)
    st.u[:] = 5.0; st.u[5, :, :] += 3.0; st.w[:, :, 3] = 1.0
    du, dv, dw = momentum_advection_tendency(st, g, order=2, periodic=True)
    for a in (du, dv, dw):
        assert np.isfinite(a).all()


# --------------------------------------------------------------------------
# item 2: f-plane Coriolis
# --------------------------------------------------------------------------
def test_coriolis_f_value():
    assert coriolis_f(0.0) == 0.0
    f36 = coriolis_f(36.0)
    assert 8e-5 < f36 < 9e-5           # ~8.6e-5 at 36 N
    assert coriolis_f(-36.0) == -f36   # sign flips in the SH


def test_coriolis_rotates_perturbation_only():
    """With u0=v0=perturbation, f turns u into v; a pure base-state wind (no
    perturbation) feels no force."""
    g = _periodic_grid()
    st = FlowState.zeros(g)
    u0 = np.zeros(g.u_shape); v0 = np.zeros(g.v_shape)
    st.u[:] = 10.0                                     # perturbation (u0=0)
    add_coriolis(st, g, dt=1.0, f=1e-3, u0_face=u0, v0_face=v0)
    assert float(st.v.mean()) < 0.0                    # dv/dt = -f u  < 0
    # base-state wind only -> no force
    st2 = FlowState.zeros(g); u0b = np.full(g.u_shape, 8.0)
    st2.u[:] = 8.0
    v_before = st2.v.copy()
    add_coriolis(st2, g, dt=1.0, f=1e-3, u0_face=u0b, v0_face=np.zeros(g.v_shape))
    assert np.allclose(st2.v, v_before)


# --------------------------------------------------------------------------
# item 3: LES closure
# --------------------------------------------------------------------------
def test_les_viscosity_positive_and_strain_dependent():
    g = _periodic_grid()
    st = FlowState.zeros(g)
    les = LESConfig(model="smagorinsky", nu_background=0.5)
    Km_rest = strain_and_viscosity(st, g, les)
    assert float(Km_rest.min()) >= 0.5 - 1e-9         # background floor at rest
    # add shear -> more eddy viscosity
    st.u[:] = np.linspace(0, 20, g.nz).reshape(1, 1, -1)
    Km_shear = strain_and_viscosity(st, g, les)
    assert float(Km_shear.max()) > float(Km_rest.max())


def test_les_model_option_is_honest():
    """The advertised LES models behave as documented: 'none' -> background only,
    'tke15' -> NotImplemented (documented future work), unknown -> error."""
    import pytest
    g = _periodic_grid()
    st = FlowState.zeros(g)
    st.u[:] = np.linspace(0, 20, g.nz).reshape(1, 1, -1)
    none = strain_and_viscosity(st, g, LESConfig(model="none", nu_background=0.7))
    assert float(none.max()) == float(none.min()) == 0.7   # constant background
    with pytest.raises(NotImplementedError):
        strain_and_viscosity(st, g, LESConfig(model="tke15"))
    with pytest.raises(ValueError):
        strain_and_viscosity(st, g, LESConfig(model="bogus"))


# --------------------------------------------------------------------------
# item 4: surface drag
# --------------------------------------------------------------------------
def test_surface_drag_retards_lowest_level_only():
    g = _periodic_grid()
    st = FlowState.zeros(g)
    st.u[:] = 15.0
    drag = SurfaceDragConfig(enabled=True, C_d=0.02)
    u0_before = float(np.abs(st.u[:, :, 0]).mean())
    u_top_before = float(np.abs(st.u[:, :, -1]).mean())
    apply_surface_drag(st, g, dt=10.0, drag=drag)
    assert float(np.abs(st.u[:, :, 0]).mean()) < u0_before     # retarded at ground
    assert np.isclose(float(np.abs(st.u[:, :, -1]).mean()), u_top_before)  # aloft untouched


def test_surface_drag_disabled_is_noop():
    g = _periodic_grid()
    st = FlowState.zeros(g); st.u[:] = 12.0
    before = st.u.copy()
    apply_surface_drag(st, g, 10.0, SurfaceDragConfig(enabled=False))
    assert np.array_equal(st.u, before)


# --------------------------------------------------------------------------
# item 5: curved-hodograph environment (SRH / shear)
# --------------------------------------------------------------------------
def _deep_grid(nz=60):
    return Grid(nx=4, ny=4, nz=nz, Lx=16000, Ly=16000, Lz=18000)


def test_curved_hodograph_has_more_srh_than_straight():
    g = _deep_grid()
    straight = snd.build_sounding(g, HodographConfig(kind="unidirectional"))
    curved = snd.build_sounding(g, HodographConfig(kind="quarter_circle"))
    srh_s = snd.storm_relative_helicity(straight, 3000.0)
    srh_c = snd.storm_relative_helicity(curved, 3000.0)
    assert srh_s > 0.0 and srh_c > 0.0            # right-mover: positive SRH
    assert srh_c > srh_s + 100.0                  # curvature adds SRH


def test_curved_hodograph_supercell_ranges():
    g = _deep_grid()
    b = snd.build_sounding(g, HodographConfig(kind="quarter_circle"))
    assert 100.0 < snd.storm_relative_helicity(b, 1000.0) < 800.0
    assert 150.0 < snd.storm_relative_helicity(b, 3000.0) < 900.0
    assert snd.bulk_shear(b, 0, 6000) > 20.0      # supercell 0-6 km shear
    # v0 populated (curved) vs unidirectional (v0 == 0)
    assert np.any(np.abs(b.v0) > 1.0)


def test_unidirectional_v0_is_zero():
    g = _deep_grid()
    b = snd.build_sounding(g, HodographConfig(kind="unidirectional"))
    assert np.allclose(b.v0, 0.0)


# --------------------------------------------------------------------------
# assembled core: short integration keeps water / mass conservation
# --------------------------------------------------------------------------
def test_storm_short_run_conserves_and_is_finite():
    scfg = build_storm_config(preset="storm-quick", nx=14, ny=14, nz=28, Lz=16000,
                              duration=45.0, dt_max=3.0)
    from storm_dynamics.core import StormSimulation
    sim = StormSimulation(scfg)
    report = sim.run()
    cons = report["conservation"]
    assert abs(cons["total_water_rel_err"]) < 1e-2          # water closes (repo std)
    assert cons["mass_continuity_residual_norm"] < 1e-2     # projection, not limiters
    assert np.isfinite(report["rotation"]["zeta_abs_max"])


def test_storm_kernel_coupling_reuses_validated_nucleation():
    """Optional: feed the validated nucleation kernel rate J as the microphysics
    embryo source (item 6, as meteorological_flow does).  A tiny lookup keeps the
    build fast; the coupled run must stay finite and conserve water."""
    from storm_dynamics.core import StormSimulation
    scfg = build_storm_config(preset="storm-quick", nx=12, ny=12, nz=24, Lz=15000,
                              duration=20.0, dt_max=3.0, couple_nucleation=True,
                              nucleation_method="lookup", outdir="outputs/_test_kernel")
    lk = scfg.sim.nucleation.lookup
    lk.n_T = 6; lk.n_pv = 5; lk.n_grad = 3; lk.scan_resolution = 12; lk.rebuild = True
    sim = StormSimulation(scfg)
    assert sim.couple_nucleation and sim.adapter is not None and sim.lookup is not None
    report = sim.run()
    assert report["kernel_coupled"] is True
    assert hasattr(sim, "last_nf")                         # kernel evaluated each step
    assert np.isfinite(report["rotation"]["zeta_abs_max"])
    assert abs(report["conservation"]["total_water_rel_err"]) < 1e-2


def test_storm_default_has_no_kernel_coupling():
    scfg = build_storm_config(preset="storm-quick", nx=12, ny=12, nz=24, duration=5.0)
    from storm_dynamics.core import StormSimulation
    sim = StormSimulation(scfg)
    assert sim.couple_nucleation is False and sim.adapter is None


def test_storm_config_disables_rayleigh_and_clip_path():
    """The fork must not carry the demonstration Rayleigh drag (it would damp the
    vortex); the LES closure is the dissipation instead."""
    scfg = build_storm_config(preset="storm")
    assert scfg.sim.flow.gamma_damp == 0.0
    assert scfg.dyn.les.model == "smagorinsky"
    assert scfg.dyn.momentum_advection is True


# --------------------------------------------------------------------------
# compute backend (CPU default; optional GPU equivalence)
# --------------------------------------------------------------------------
def test_storm_defaults_to_cpu_backend():
    scfg = build_storm_config(preset="storm-quick", nx=12, ny=12, nz=24, duration=5.0)
    assert scfg.sim.performance.device == "cpu"
    from storm_dynamics.core import StormSimulation
    sim = StormSimulation(scfg)
    assert sim.backend.name == "cpu"
    assert sim.run()["backend"]["name"] == "cpu"


def _gpu_available() -> bool:
    try:
        import cupy
        _ = cupy.cuda.Device(0).compute_capability
        return True
    except Exception:
        return False


@pytest.mark.skipif(not _gpu_available(),
                    reason="no working CUDA/CuPy GPU in this environment")
def test_storm_gpu_matches_cpu():
    """The rotating core is backend-agnostic (grid.xp): a short run on the GPU
    reproduces the CPU rotation diagnostics and conservation."""
    from storm_dynamics.core import StormSimulation

    def run(device):
        scfg = build_storm_config(preset="storm", nx=20, ny=20, nz=32,
                                  Lx=32000, Ly=32000, Lz=15000, duration=120.0,
                                  dt_max=3.0, hodograph_kind="unidirectional",
                                  drag=False, device=device)
        rep = StormSimulation(scfg).run()
        return rep

    cpu = run("cpu")
    gpu = run("gpu")
    assert gpu["backend"]["name"] == "gpu"
    for k in ("w_max", "zeta_abs_max", "midlevel_mesocyclone", "updraft_helicity_max"):
        assert np.isclose(cpu["rotation"][k], gpu["rotation"][k], rtol=1e-6, atol=1e-9), k
    assert np.isclose(cpu["conservation"]["total_water_rel_err"],
                      gpu["conservation"]["total_water_rel_err"], rtol=1e-6, atol=1e-9)
