"""Surface sensible + latent heat fluxes (storm_dynamics.surface_fluxes)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from storm_dynamics.config import build_storm_config, SurfaceFluxConfig
from storm_dynamics.core import StormSimulation
from storm_dynamics import surface_fluxes as sfl


def _sim(**flux):
    scfg = build_storm_config(preset="storm", nx=20, ny=20, nz=28, Lx=20000.0, Ly=20000.0,
                              Lz=12000.0, duration=1.0, dt_max=3.0, device="cpu")
    if flux:
        scfg.dyn.fluxes = SurfaceFluxConfig(**flux)
    sim = StormSimulation(scfg)
    sim.state.diagnose(sim.cfg)
    return sim


def test_disabled_is_noop():
    sim = _sim()
    assert sim.dyn.fluxes.enabled is False
    th0 = np.asarray(sim.state.theta[:, :, 0]).copy()
    sfl.apply_surface_fluxes(sim.state, sim.grid, 3.0, sim.dyn.fluxes, base=sim.base)
    assert np.array_equal(np.asarray(sim.state.theta[:, :, 0]), th0)


def test_warm_ground_heats_lowest_level():
    sim = _sim(enabled=True, C_h=2e-2, dtheta_sfc_K=6.0, saturate_surface=False)
    th_before = float(np.mean(np.asarray(sim.state.theta[:, :, 0])))
    prev = th_before
    for _ in range(15):
        sfl.apply_surface_fluxes(sim.state, sim.grid, 3.0, sim.dyn.fluxes, base=sim.base)
        cur = float(np.mean(np.asarray(sim.state.theta[:, :, 0])))
        assert cur >= prev - 1e-12               # monotone relaxation toward the warmer surface
        prev = cur
    assert prev > th_before + 1e-3               # net heating


def test_moist_ground_moistens_lowest_level():
    sim = _sim(enabled=True, C_q=5e-3, saturate_surface=True)
    qv_before = float(np.mean(np.asarray(sim.state.qv[:, :, 0])))
    for _ in range(10):
        sim.state.diagnose(sim.cfg)             # refresh T/P for the saturation target
        sfl.apply_surface_fluxes(sim.state, sim.grid, 3.0, sim.dyn.fluxes, base=sim.base)
    qv_after = float(np.mean(np.asarray(sim.state.qv[:, :, 0])))
    assert qv_after >= qv_before                 # driven toward surface saturation (or unchanged)


def test_flux_report_positive_for_warm_ground():
    sim = _sim(enabled=True, C_h=5e-3, dtheta_sfc_K=6.0, saturate_surface=False)
    rep = sfl.surface_flux_report(sim.state, sim.grid, sim.dyn.fluxes, base=sim.base)
    assert rep["sensible_heat_flux_W_m2"] > 0.0


def test_neutral_drag_coefficient_monotone():
    c_smooth = sfl.neutral_drag_coefficient(10.0, 0.01)
    c_rough = sfl.neutral_drag_coefficient(10.0, 0.5)
    assert c_rough > c_smooth > 0.0              # rougher surface -> larger transfer coefficient


if __name__ == "__main__":
    for nm, fn in sorted(globals().items()):
        if nm.startswith("test_") and callable(fn):
            fn(); print("ok", nm)
    print("ALL SURFACE-FLUX TESTS PASSED")
