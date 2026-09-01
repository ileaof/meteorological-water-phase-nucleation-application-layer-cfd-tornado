"""ROADMAP §3a -- real initial conditions: ingest an observed sounding into the storm's
base state (`storm_dynamics.soundings.from_observed_sounding`).  The first concrete step from
the analytic Weisman-Klemp column toward a real-atmosphere environment."""
import warnings

import numpy as np

from meteorological_flow.grid import Grid
from meteorological_flow.base_state import weisman_klemp, sounding_diagnostics
from storm_dynamics import soundings as snd


def _grid():
    return Grid(nx=8, ny=8, nz=40, Lx=32000.0, Ly=32000.0, Lz=15000.0,
                z_stretch=1.05, periodic=True)


def test_observed_sounding_ingestor_round_trips_a_known_atmosphere():
    """Feeding a known analytic atmosphere (Weisman-Klemp) back through the observed-sounding
    ingestor as radiosonde columns recovers it: theta/qv/winds to round-off and -- with the
    true surface pressure -- the hydrostatic p0 and CAPE to <1%.  Verifies the conversion and
    the hydrostatic re-integration are correct, with no external data."""
    g = _grid()
    wk = weisman_klemp(g)
    rec = snd.from_observed_sounding(
        g, pressure_hPa=np.asarray(wk.p0) / 100.0, height_m=np.asarray(wk.zc),
        temperature_C=np.asarray(wk.T0) - 273.15, qv_kgkg=np.asarray(wk.qv0),
        u_ms=np.asarray(wk.u0), v_ms=np.asarray(wk.v0), p_sfc_Pa=1.0e5)
    assert np.abs(rec.theta0 - wk.theta0).max() < 1e-9        # theta recovered to round-off
    assert np.abs(rec.qv0 - wk.qv0).max() < 1e-12             # qv exact (passed through)
    assert np.abs(rec.u0 - wk.u0).max() < 1e-12 and np.abs(rec.v0 - wk.v0).max() < 1e-12
    assert np.abs(rec.p0 - wk.p0).max() < 1.0                 # hydrostatic p0 to ~1 Pa
    c_wk = sounding_diagnostics(wk)["CAPE_J_kg"]; c_rec = sounding_diagnostics(rec)["CAPE_J_kg"]
    assert abs(c_rec - c_wk) < 0.01 * c_wk + 1.0              # CAPE recovered to <1%


def test_example_observed_sounding_is_a_supercell_environment():
    """The bundled illustrative radiosonde sounding builds a physically sane, supercell-
    supporting base state: sizeable CAPE, strong deep-layer shear and 0-3 km SRH, a
    monotonically decreasing density and non-negative moisture."""
    g = _grid()
    base = snd.from_observed_sounding(g, **snd.example_observed_sounding())
    d = sounding_diagnostics(base)
    assert 1500.0 < d["CAPE_J_kg"] < 5000.0, d["CAPE_J_kg"]   # strong but not unphysical
    assert d["shear_0_6km_m_s"] > 20.0, d["shear_0_6km_m_s"]  # deep-layer shear (supercell)
    assert snd.storm_relative_helicity(base) > 100.0          # veering hodograph -> +SRH
    assert d["LCL_m"] is not None and d["LFC_m"] is not None
    assert np.all(np.diff(base.rho0) < 0.0)                   # density decreases with height
    assert np.all(base.qv0 >= 0.0) and np.all(np.isfinite(base.p0))


def test_storm_runs_stably_on_the_ingested_observed_base():
    """A StormSimulation initialised from the ingested observed sounding (passed as
    `base=`) steps stably -- finite rotation, small mass-continuity residual -- so real
    initial conditions drop into the existing solver unchanged."""
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    scfg = build_storm_config(preset="storm", nx=16, ny=16, nz=40, Lx=32000.0, Ly=32000.0,
                              Lz=15000.0, duration=1.0, dt_max=3.0, drag=True,
                              z_stretch=1.05, device="cpu")
    g = Grid(nx=16, ny=16, nz=40, Lx=32000.0, Ly=32000.0, Lz=15000.0,
             z_stretch=1.05, periodic=True)
    base = snd.from_observed_sounding(g, **snd.example_observed_sounding())
    assert np.allclose(np.asarray(base.zc), np.asarray(g.zc))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sim = StormSimulation(scfg, base=base)
        sim.cfg.time.duration = 4 * float(sim._dt())
        rep = sim.run()
    assert np.isfinite(rep["rotation"]["zeta_abs_max"])
    assert rep["conservation"]["mass_continuity_residual_norm"] < 1e-3, rep["conservation"]
