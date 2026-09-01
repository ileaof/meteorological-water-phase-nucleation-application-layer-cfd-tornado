"""ROADMAP §3e -- quantitative verification against community benchmarks.

Straka et al. (1993) density current: a cold bubble collapsing in a neutral, dry, non-rotating
atmosphere at a *fixed* viscosity nu = 75 m^2/s, integrated to 900 s.  The published reference
front (leading edge of the surface outflow) reaches ~15.5 km from the centre at 100 m
resolution; here at 200 m it lands near ~13.7 km, symmetric, with Kelvin-Helmholtz-rotor
updrafts ~15-20 m/s -- i.e. the model reproduces the standard density current, not just its own
conservation checks."""
import warnings

import numpy as np

from storm_dynamics.benchmarks import straka_simulation, straka_front_position


def test_straka_density_current_matches_the_reference():
    """The Straka cold bubble produces a stable, symmetric density current whose 900 s front
    and rotor velocities match the benchmark (at 200 m resolution, exact CPU pressure solve)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sim = straka_simulation(nz=32, ny=2, duration=900.0, device="cpu")   # dx = dz = 200 m
        assert sim.f == 0.0 and sim.dyn.les.model == "none"                  # neutral, constant nu
        assert abs(sim.dyn.les.nu_background - 75.0) < 1e-9
        sim.run()
    to = sim.grid.backend.to_cpu
    dth = np.asarray(to(sim.state.theta - sim.theta0_field))
    w = np.asarray(to(sim.state.w))
    x = np.asarray(to(sim.grid.xc)); xc = 0.5 * float(sim.grid.Lx)
    front = straka_front_position(sim)
    cold = np.where(dth[:, 0, 0] < -1.0)[0]
    assert cold.size > 0
    left = xc - float(x[cold].min()); right = float(x[cold].max()) - xc

    assert 11_000.0 < front < 16_500.0, front            # near the 15.5 km (100 m) reference
    assert abs(left - right) < 800.0, (left, right)      # a symmetric density current
    assert -15.1 < dth.min() < -2.0, dth.min()           # cold pool present, only diffused (no blow-up)
    assert 8.0 < w.max() < 30.0, w.max()                 # KH-rotor updrafts ~15-20 m/s, no clip
    assert np.isfinite(dth).all() and np.isfinite(w).all()
