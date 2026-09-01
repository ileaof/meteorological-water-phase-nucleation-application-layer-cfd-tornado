"""ROADMAP §3a -- limited-area lateral boundary conditions (Davies relaxation zone)."""
import numpy as np

from meteorological_flow.grid import Grid
from meteorological_flow.state import FlowState
from meteorological_flow.base_state import build_base_state
from storm_dynamics.limited_area import (
    lateral_relaxation_weight, environment_target, apply_lateral_relaxation)


def _grid():
    return Grid(nx=40, ny=40, nz=16, Lx=40000.0, Ly=40000.0, Lz=12000.0, periodic=False)


def test_lateral_relaxation_is_consistent_and_absorbs_at_the_boundary():
    """The Davies zone (a) leaves a state that already equals the target unchanged
    (consistency), and (b) drives a boundary perturbation to the target while leaving the
    interior untouched -- the sponge that stops reflections in a limited-area run."""
    g = _grid()
    base = build_base_state(g)
    target = environment_target(g, base)

    # (a) consistency: state == target -> zero tendency
    st = FlowState.zeros(g)
    st.u = target["u"].copy(); st.v = target["v"].copy(); st.w = target["w"].copy()
    st.theta = target["theta"].copy(); st.qv = target["qv"].copy()
    apply_lateral_relaxation(st, g, target, dt=1.0, width=8, rate=1.0 / 50.0)
    assert np.abs(np.asarray(st.theta) - np.asarray(target["theta"])).max() < 1e-12
    assert np.abs(np.asarray(st.u) - np.asarray(target["u"])).max() < 1e-12

    # (b) sponge: start target + a uniform perturbation; relax many steps
    delta = 5.0
    st.theta = target["theta"] + delta
    w = lateral_relaxation_weight(g, width=8, rate=1.0 / 50.0)
    for _ in range(400):
        apply_lateral_relaxation(st, g, target, dt=1.0, weight=w)
    dev = np.asarray(st.theta - target["theta"])
    edge = abs(float(dev[0, g.ny // 2, 0]))           # outermost cell (full weight)
    interior = abs(float(dev[g.nx // 2, g.ny // 2, 0]))   # centre (zero weight)
    assert edge < 0.05 * delta, edge                   # boundary nudged to the target
    assert interior > 0.95 * delta, interior           # interior untouched by the lateral BC
    assert np.all(np.isfinite(dev))


def test_lateral_relaxation_weight_profile():
    """The weight is `rate` at the edge, decreasing monotonically to 0 by `width` cells in."""
    g = _grid()
    w = np.asarray(lateral_relaxation_weight(g, width=6, rate=0.02))[:, g.ny // 2, 0]
    assert abs(w[0] - 0.02) < 1e-12                     # rate at the outermost cell
    assert w[6] == 0.0 and w[g.nx // 2] == 0.0          # zero by `width` cells in / interior
    assert np.all(np.diff(w[:7]) <= 1e-15)             # monotonically ramping down inward
