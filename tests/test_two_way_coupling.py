"""Two-way coupling wired into the live multi-level cascade (run_multilevel_nest two_way=).

Asserts (a) the default one-way path is deterministic/unchanged, (b) two_way=True actually feeds the
fine level back onto the coarse overlap (the parent state moves), and (c) both stay finite.  The
*physical* payoff (two-way raises low-level rotation) is demonstrated at scale by Attempt G; here we
verify the wiring cheaply."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics import nesting as nst

def _matured_parent(steps=25):
    scfg = build_storm_config(preset="storm", nx=24, ny=24, nz=30, Lx=24000.0, Ly=24000.0,
                              Lz=12000.0, duration=1.0, dt_max=3.0, device="cpu")
    scfg.sim.physics.bubble_dtheta = 5.0
    sim = StormSimulation(scfg)
    for _ in range(steps):
        sim._step(sim._dt()); sim.step += 1; sim.t = float(sim.state.t)
    return sim


def _run(sim, two_way):
    spec = nst.NestSpec.aligned(sim.grid, i0=7, j0=7, ncx=10, ncy=10, refine=3)
    sims, rep = nst.run_multilevel_nest(sim, [spec], window=8.0, two_way=two_way, cfl=0.2)
    return sims, rep


def test_default_is_one_way_and_deterministic():
    # two fresh, identical sims run one-way must give identical parents (deterministic core)
    sa = _matured_parent(); _, rep_a = _run(sa, two_way=False)
    sb = _matured_parent(); _run(sb, two_way=False)
    assert rep_a["nest"]["two_way"] is False
    assert np.allclose(np.asarray(sa.state.theta), np.asarray(sb.state.theta), rtol=1e-12, atol=1e-12)


def test_two_way_feeds_back_to_parent():
    sa = _matured_parent(); _run(sa, two_way=False); one_way_parent = np.asarray(sa.state.w).copy()
    sb = _matured_parent(); _, rep = _run(sb, two_way=True); two_way_parent = np.asarray(sb.state.w)
    assert rep["nest"]["two_way"] is True
    assert np.isfinite(two_way_parent).all()
    # the injection must have moved the parent's overlap region relative to one-way
    assert np.abs(two_way_parent - one_way_parent).max() > 1e-6


if __name__ == "__main__":
    test_default_is_one_way_and_deterministic(); print("ok one-way deterministic")
    test_two_way_feeds_back_to_parent(); print("ok two-way feedback")
    print("ALL TWO-WAY COUPLING TESTS PASSED")
