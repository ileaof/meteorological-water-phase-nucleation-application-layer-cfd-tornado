"""Unified tornadogenesis config overlay (storm_dynamics.tornado_config)."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from storm_dynamics import tornado_config as tc

CFG = os.path.join(os.path.dirname(__file__), "..", "config", "tornadogenesis.yaml")


def test_impose_vortex_true_is_rejected():
    with pytest.raises(ValueError):
        tc.load_tornadogenesis_config({"tornadogenesis": {"impose_vortex": True}})


def test_load_example_yaml_builds_storm_config():
    cfg = tc.load_tornadogenesis_config(CFG)
    assert cfg.impose_vortex is False
    assert cfg.storm.sim.grid.nx == 64 and cfg.storm.sim.grid.nz == 44
    assert cfg.diagnostics["vorticity_budget"] is True
    assert cfg.nesting.get("two_way") is True


def test_surface_fluxes_wired_when_enabled():
    d = {"surface": {"sensible_heat_flux": True, "latent_heat_flux": True, "C_h": 2e-3},
         "domain": {"nx": 16, "ny": 16, "nz": 20}}
    cfg = tc.load_tornadogenesis_config(d)
    assert cfg.storm.dyn.fluxes.enabled is True
    assert cfg.storm.dyn.fluxes.C_h == 2e-3
    assert cfg.storm.dyn.fluxes.saturate_surface is True


def test_defaults_no_fluxes():
    cfg = tc.load_tornadogenesis_config({"domain": {"nx": 16, "ny": 16, "nz": 20}})
    assert cfg.storm.dyn.fluxes.enabled is False       # off unless requested


def test_run_diagnostics_bundle():
    from storm_dynamics.core import StormSimulation
    cfg = tc.load_tornadogenesis_config({"domain": {"nx": 20, "ny": 20, "nz": 28, "Lz_m": 12000.0},
                                         "device": "cpu"})
    sim = StormSimulation(cfg.storm)
    for _ in range(3):
        sim._step(sim._dt()); sim.step += 1; sim.t = float(sim.state.t)
    out = tc.run_diagnostics(sim, cfg)
    for k in ("rotation", "vorticity_budget_low", "vortex", "cold_pool", "micro_gradient",
              "classification", "dominant_low_mechanism"):
        assert k in out
    assert out["classification"] in __import__("storm_dynamics.classification",
                                               fromlist=["CATEGORIES"]).CATEGORIES


if __name__ == "__main__":
    for nm, fn in sorted(globals().items()):
        if nm.startswith("test_") and callable(fn):
            fn(); print("ok", nm)
    print("ALL TORNADO-CONFIG TESTS PASSED")
