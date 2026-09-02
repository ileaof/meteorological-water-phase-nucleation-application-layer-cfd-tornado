"""The diagnostic figure suite renders a PNG from a live sim (examples/tornado_diagnostic_figures)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "examples"))


def test_make_diagnostic_figure(tmp_path):
    import tornado_diagnostic_figures as tdf
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    scfg = build_storm_config(preset="storm", nx=16, ny=16, nz=24, Lx=16000.0, Ly=16000.0,
                              Lz=12000.0, duration=1.0, dt_max=4.0, z_stretch=1.05, device="cpu")
    scfg.sim.physics.bubble_dtheta = 5.0
    sim = StormSimulation(scfg)
    for _ in range(12):
        sim._step(sim._dt()); sim.step += 1; sim.t = float(sim.state.t)
    out = str(tmp_path / "diag.png")
    tdf.make_diagnostic_figure(sim, out)
    assert os.path.exists(out) and os.path.getsize(out) > 5000


if __name__ == "__main__":
    import tempfile

    class _P:
        def __truediv__(self, n): return os.path.join(tempfile.mkdtemp(), n)
    test_make_diagnostic_figure(_P()); print("ALL DIAGNOSTIC-FIGURE TESTS PASSED")
