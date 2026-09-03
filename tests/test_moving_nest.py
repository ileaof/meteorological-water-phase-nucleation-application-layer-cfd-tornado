"""Moving nest in the live multi-level cascade (criterion 8): the refined domain tracks the storm
along a FILTERED trajectory, keeping its size so the fine field transfers by exact integer shift."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics import nesting as nst


def _matured(steps=30):
    scfg = build_storm_config(preset="storm", nx=28, ny=28, nz=30, Lx=28000.0, Ly=28000.0,
                              Lz=12000.0, duration=1.0, dt_max=3.0, device="cpu")
    scfg.sim.physics.bubble_dtheta = 5.0
    sim = StormSimulation(scfg)
    for _ in range(steps):
        sim._step(sim._dt()); sim.step += 1; sim.t = float(sim.state.t)
    return sim


def test_follow_spec_keeps_size_and_filters():
    sim = _matured()
    old = nst.NestSpec.aligned(sim.grid, i0=2, j0=2, ncx=10, ncy=10, refine=3)
    ns = nst.follow_spec(old, sim, field="zeta", frac=0.3, alpha=1.0)
    if ns is not None:                                    # something was tagged
        assert np.isclose(ns.Lx, old.Lx) and np.isclose(ns.Ly, old.Ly)   # SAME SIZE
        assert ns.refine == old.refine
        # a filtered (alpha=0.5) move must land between the old origin and the alpha=1.0 target
        half = nst.follow_spec(old, sim, field="zeta", frac=0.3, alpha=0.5)
        if half is not None:
            lo, hi = sorted((old.x0, ns.x0))
            assert lo - 1e-6 <= half.x0 <= hi + 1e-6


def test_follow_spec_returns_none_when_nothing_tagged():
    sim = _matured(steps=1)
    sim.state.u[:] = 0.0; sim.state.v[:] = 0.0; sim.state.w[:] = 0.0   # no rotation at all
    old = nst.NestSpec.aligned(sim.grid, i0=2, j0=2, ncx=10, ncy=10, refine=3)
    assert nst.follow_spec(old, sim, field="zeta", frac=0.5) is None


def test_multilevel_moving_nest_runs_and_reports():
    sim = _matured()
    spec = nst.NestSpec.aligned(sim.grid, i0=8, j0=8, ncx=10, ncy=10, refine=3)
    sims, rep = nst.run_multilevel_nest(sim, [spec], window=8.0, cfl=0.2,
                                        follow_interval=1, follow_field="zeta",
                                        follow_frac=0.3, follow_filter=0.5)
    assert rep["nest"]["moving_nest"] is True
    assert rep["nest"]["nest_moves"] >= 0
    assert np.isfinite(np.asarray(sims[-1].state.w)).all()     # stays stable while moving


def test_default_is_fixed_frame_and_deterministic():
    a = _matured(); spec_a = nst.NestSpec.aligned(a.grid, i0=8, j0=8, ncx=10, ncy=10, refine=3)
    _, rep_a = nst.run_multilevel_nest(a, [spec_a], window=8.0, cfl=0.2)
    b = _matured(); spec_b = nst.NestSpec.aligned(b.grid, i0=8, j0=8, ncx=10, ncy=10, refine=3)
    nst.run_multilevel_nest(b, [spec_b], window=8.0, cfl=0.2)
    assert rep_a["nest"]["moving_nest"] is False and rep_a["nest"]["nest_moves"] == 0
    assert np.allclose(np.asarray(a.state.theta), np.asarray(b.state.theta), rtol=1e-12, atol=1e-12)


if __name__ == "__main__":
    for nm, fn in sorted(globals().items()):
        if nm.startswith("test_") and callable(fn):
            fn(); print("ok", nm)
    print("ALL MOVING-NEST TESTS PASSED")
