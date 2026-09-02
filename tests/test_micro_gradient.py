"""Two-scale temperature gradient + closure (storm_dynamics.micro_gradient)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from meteorological_flow.grid import Grid
from meteorological_flow.state import FlowState
from storm_dynamics import micro_gradient as mg


def _grid(backend=None):
    return Grid(nx=24, ny=24, nz=20, Lx=2400.0, Ly=2400.0, Lz=2000.0, z_stretch=1.0,
                periodic=True, backend=backend)


def _state_linear_T(g, slope=0.01):
    st = FlowState.zeros(g)
    xp = g.xp
    X = xp.asarray(g.xc)[:, None, None] + g.zeros_c()
    st.T = 290.0 + slope * X                 # dT/dx = slope -> |grad T|_macro = slope
    st.gradT_mag = None                      # force recomputation from T
    st.S_w = 1.0 + 0.02 + g.zeros_c()        # 2% supersaturation
    return st


def test_macro_gradient_of_linear_field():
    g = _grid(); slope = 0.01
    st = _state_linear_T(g, slope)
    macro = np.asarray(g.backend.to_cpu(mg.macro_temperature_gradient(st, g)))
    assert np.allclose(macro[4:-4, 4:-4, 4:-4], slope, rtol=1e-6)


def test_micro_exceeds_macro_and_scales():
    g = _grid(); slope = 0.01
    st = _state_linear_T(g, slope)
    micro, enh = mg.micro_temperature_gradient(st, g, eps=1e-3)
    micro = np.asarray(g.backend.to_cpu(micro)); enh = np.asarray(g.backend.to_cpu(enh))
    assert (enh > 1.0).all()                                   # sub-grid sharpening
    macro = np.asarray(g.backend.to_cpu(mg.macro_temperature_gradient(st, g)))
    assert np.allclose(micro, macro * enh, rtol=1e-10)          # micro = macro * enhancement


def test_stronger_turbulence_sharpens_gradient():
    g = _grid()
    st = _state_linear_T(g)
    _, enh_lo = mg.micro_temperature_gradient(st, g, eps=1e-5)
    _, enh_hi = mg.micro_temperature_gradient(st, g, eps=1e-1)
    # larger eps -> smaller Batchelor scale -> larger enhancement
    assert float(np.asarray(g.backend.to_cpu(enh_hi)).max()) > float(np.asarray(g.backend.to_cpu(enh_lo)).max())


def test_batchelor_scale_monotone():
    assert mg.batchelor_scale(1e-5) > mg.batchelor_scale(1e-1)  # weaker turbulence -> larger scale


def test_named_diagnostics():
    g = _grid()
    st = _state_linear_T(g)
    d = mg.nucleation_diagnostics(st, g)
    for k in ("macro_temperature_gradient_max_K_m", "micro_temperature_gradient_max_K_m",
              "micro_macro_enhancement_max", "local_supersaturation_max",
              "nucleation_rate_proxy_max", "latent_heat_release_proxy_max"):
        assert k in d
    assert d["micro_temperature_gradient_max_K_m"] > d["macro_temperature_gradient_max_K_m"]
    assert abs(d["local_supersaturation_max"] - 0.02) < 1e-6


if __name__ == "__main__":
    for nm, fn in sorted(globals().items()):
        if nm.startswith("test_") and callable(fn):
            fn(); print("ok", nm)
    print("ALL MICRO-GRADIENT TESTS PASSED")
