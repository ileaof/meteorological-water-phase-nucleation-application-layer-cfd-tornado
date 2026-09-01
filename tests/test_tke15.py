"""TKE-1.5 (Deardorff 1980) prognostic subgrid closure (ROADMAP §3c).

`storm_dynamics/turbulence.py` used to raise NotImplementedError for `tke15`; it now evolves a
prognostic subgrid TKE that sets the eddy viscosity.  The idealised default (`smagorinsky`) is
unchanged.
"""
import warnings

import numpy as np
import pytest

from storm_dynamics import turbulence as tb
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation


def _visc(e, n2, Delta):
    return tb.deardorff_viscosity(np.asarray(e), N2=np.full_like(np.asarray(e), n2),
                                  Delta=Delta, les=None, xp=np)


def test_deardorff_viscosity_relations():
    """K_m = C_k l sqrt(e); K_h >= K_m; dissipation eps > 0; length reduced in stable N^2."""
    e = np.array([[[0.25]]]); Delta = 100.0
    Km, Kh, l, eps = _visc(e, 0.0, Delta)
    assert Km[0, 0, 0] == pytest.approx(0.10 * Delta * 0.5, rel=1e-6)   # C_k=0.1, l=Delta, sqrt(e)=0.5
    assert Kh[0, 0, 0] >= Km[0, 0, 0] and eps[0, 0, 0] > 0.0
    _, _, l_stable, _ = _visc(e, 4e-3, Delta)              # stable N^2 shortens the mixing length
    assert l_stable[0, 0, 0] < l[0, 0, 0]


def test_tke15_runs_stably_and_evolves_a_positive_tke_field():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scfg = build_storm_config(preset="storm", nx=16, ny=16, nz=20, Lx=16000.0, Ly=16000.0,
                                  Lz=10000.0, duration=1.0, dt_max=3.0, drag=True,
                                  z_stretch=1.05, les_model="tke15", device="cpu")
        sim = StormSimulation(scfg)
        assert sim.dyn.les.model == "tke15"
        sim.cfg.time.duration = 6 * float(sim._dt())
        rep = sim.run()
    to = sim.grid.backend.to_cpu
    assert np.isfinite(rep["rotation"]["zeta_abs_max"])
    tke = np.asarray(to(sim._tke))
    assert tke.min() >= 0.0 and tke.max() > 1e-3            # a real, positive subgrid TKE field
    Km = np.asarray(to(sim._Km))
    assert Km.max() > Km.min()                             # strain-dependent eddy viscosity
    assert rep["conservation"]["mass_continuity_residual_norm"] < 1e-2


def test_smagorinsky_default_unchanged():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scfg = build_storm_config(preset="storm", nx=14, ny=14, nz=16, Lx=14000.0, Ly=14000.0,
                                  Lz=9000.0, duration=1.0, dt_max=3.0, les_model="smagorinsky",
                                  device="cpu")
        sim = StormSimulation(scfg)
        sim.cfg.time.duration = 4 * float(sim._dt())
        rep = sim.run()
    assert np.isfinite(rep["rotation"]["zeta_abs_max"]) and getattr(sim, "_tke", None) is None
