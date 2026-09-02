"""Sustained mesoscale-ascent forcing (storm_dynamics.forcing) -- the dryline/convergence
proxy that lets a supercell establish from a real (capped) environment.  Additive + opt-in:
default OFF must leave the step bit-for-bit unchanged."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from storm_dynamics import forcing as frc
from storm_dynamics.config import MesoForcingConfig, build_storm_config
from storm_dynamics.core import StormSimulation


class _G:
    Lx = Ly = 24000.0
    xc = np.linspace(0, 24000.0, 24)
    yc = np.linspace(0, 24000.0, 24)
    zc = np.linspace(0, 12000.0, 24)
    xp = np


def _sim(nx=24):
    scfg = build_storm_config(preset="storm", nx=nx, ny=nx, nz=nx, Lx=24000.0, Ly=24000.0,
                              Lz=12000.0, duration=1.0, dt_max=3.0, device="cpu")
    return scfg


def test_default_off_is_a_noop():
    """forcing.enabled defaults False; a step must not touch theta via the forcing path."""
    scfg = _sim()
    assert scfg.dyn.forcing.enabled is False
    sim = StormSimulation(scfg)
    th = np.asarray(sim.state.theta).copy()
    # apply_meso_forcing with a disabled config returns False and changes nothing
    changed = frc.apply_meso_forcing(sim.state, sim.grid, scfg.dyn.forcing, 0.0, 3.0, xp=np)
    assert changed is False
    assert np.array_equal(np.asarray(sim.state.theta), th)


def test_mask_shape_peak_and_taper():
    fc = MesoForcingConfig(enabled=True, heat_rate_K_s=0.006, radius_m=6000.0, z_top_m=2500.0)
    m = frc.meso_forcing_mask(_G, fc, xp=np)
    assert m.shape == (24, 24, 24)
    assert 0.0 <= m.min() and m.max() <= 1.0
    # peak near the centre-surface, ~0 at the top and the horizontal edges
    assert m[12, 12, 0] > 0.8
    assert m[12, 12, -1] == 0.0
    assert m[0, 0, 0] < 1e-6


def test_enabled_heating_and_moistening_raise_core_scalars():
    fc = MesoForcingConfig(enabled=True, heat_rate_K_s=0.01, moist_rate_kgkg_s=5e-6,
                           radius_m=6000.0, z_top_m=2500.0, duration_s=1200.0)
    scfg = _sim()
    sim = StormSimulation(scfg)
    th0 = np.asarray(sim.state.theta)[12, 12, 1]
    qv0 = np.asarray(sim.state.qv)[12, 12, 1]
    dt = 5.0
    for _ in range(3):
        changed = frc.apply_meso_forcing(sim.state, sim.grid, fc, 0.0, dt, xp=np)
        assert changed is True
    assert np.asarray(sim.state.theta)[12, 12, 1] > th0 + 0.1        # ~0.01*15 = 0.15 K
    assert np.asarray(sim.state.qv)[12, 12, 1] > qv0


def test_forcing_switches_off_after_duration():
    fc = MesoForcingConfig(enabled=True, heat_rate_K_s=0.01, duration_s=600.0)
    scfg = _sim()
    sim = StormSimulation(scfg)
    assert frc.apply_meso_forcing(sim.state, sim.grid, fc, 599.0, 3.0, xp=np) is True
    assert frc.apply_meso_forcing(sim.state, sim.grid, fc, 600.0, 3.0, xp=np) is False
    assert frc.apply_meso_forcing(sim.state, sim.grid, fc, 900.0, 3.0, xp=np) is False


if __name__ == "__main__":
    for nm, fn in sorted(globals().items()):
        if nm.startswith("test_") and callable(fn):
            fn(); print("ok", nm)
    print("ALL FORCING TESTS PASSED")
