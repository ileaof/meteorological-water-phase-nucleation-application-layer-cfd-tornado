"""Tests for the real-data ingestion module (ROADMAP §3a).

Covers the task's checklist without large downloads: synthetic samples + tiny temp files
exercise the whole pipeline; the GRIB/NEXRAD real readers are guarded (skipped when the
optional library is absent) but their unavailability path is asserted.
"""
import json
import os
import tempfile
import warnings

import numpy as np
import pytest

import atmospheric_data as ad
from atmospheric_data import (thermo, units, interpolate, basestate, ic_bc, radial,
                              validation, qc)
from atmospheric_data.sources import synthetic, sounding, metar, storm_events
from atmospheric_data.sources.base import require, SourceUnavailable, available


# ---- 5,6 units + thermodynamics -------------------------------------------------
def test_thermo_transforms_roundtrip_and_hydrostatic():
    T = np.array([300.0, 275.0]); p = np.array([1.0e5, 7e4]); qv = np.array([0.015, 0.006])
    assert np.abs(thermo.temperature_from_theta(thermo.potential_temperature(T, p), p) - T).max() < 1e-9
    rho = thermo.density(p, thermo.virtual_temperature(T, qv))
    assert 1.0 < rho[0] < 1.3
    z = np.linspace(0, 15000, 40); th0 = 300 + 0.004 * z; q0 = 0.012 * np.exp(-z / 3000)
    p0, T0, r0 = thermo.hydrostatic_base_pressure(z, th0, q0, 1.0e5)
    assert p0[0] == pytest.approx(1.0e5, rel=1e-6) and np.all(np.diff(r0) < 0)
    assert thermo.specific_humidity_from_rh(0.5, 300.0, 1e5) > 0


def test_unit_conversions():
    u, v = units.wind_dir_speed_to_uv(270.0, 10.0)      # westerly
    assert u == pytest.approx(10.0, abs=1e-6) and abs(v) < 1e-6
    assert units.hpa_to_pa(1000.0) == 1.0e5
    assert units.celsius_to_kelvin(0.0) == pytest.approx(273.15)


# ---- 2 internal NetCDF format ---------------------------------------------------
def test_internal_format_roundtrip_and_standard_names():
    z = np.linspace(0, 12000, 20)
    st = ad.AtmosphericState.new(np.datetime64("2013-05-20T18:00"), z,
                                 np.linspace(0, 1e5, 8), np.linspace(0, 1e5, 9),
                                 projection="lambert_conformal", source="synthetic")
    st.add("theta", np.full((1, 20, 8, 9), 305.0), source="synthetic", original_name="pt")
    with pytest.raises(KeyError):
        st.add("not_a_standard_name", np.zeros((1, 20, 8, 9)), source="x")
    p = os.path.join(tempfile.mkdtemp(), "s.nc"); st.to_netcdf(p)
    st2 = ad.AtmosphericState.from_netcdf(p)
    assert "theta" in st2.ds and st2.ds["theta"].attrs["units"] == "K"
    assert st2.ds["theta"].attrs["source"] == "synthetic"


# ---- config + cache + offline ---------------------------------------------------
def test_config_validation_and_cache_offline():
    cfg = ad.CaseConfig(); cfg.validate()
    cfg.model.input_mode = "bogus"
    with pytest.raises(ValueError):
        cfg.validate()
    cache = ad.Cache(tempfile.mkdtemp(), offline=True)
    with pytest.raises(FileNotFoundError):
        cache.require_offline_ok("hrrr", "missing_key", ".grib2")


# ---- 3,4 source readers (synthetic + CSV) --------------------------------------
def test_synthetic_sources_and_csv_readers():
    cfg = ad.CaseConfig()
    st = synthetic.synthetic_atmosphere(cfg)
    assert {"T", "theta", "qv", "u", "v", "w", "terrain"} <= set(st.ds.data_vars)
    rad = synthetic.synthetic_radar(cfg)
    assert rad["radial_velocity"].shape == rad["reflectivity"].shape
    tmp = tempfile.mkdtemp()
    import pandas as pd
    sp = os.path.join(tmp, "s.csv")
    pd.DataFrame({"height_m": [100, 1500, 6000], "pressure_hPa": [1000, 850, 500],
                  "temperature_C": [28, 16, -9], "dewpoint_C": [21, 12, -27],
                  "wind_dir_deg": [160, 200, 255], "wind_speed_ms": [8, 16, 30]}).to_csv(sp, index=False)
    prof = sounding.read_sounding(sp)
    assert "specific_humidity" in prof and prof["u_ms"].size == 3
    ep = os.path.join(tmp, "e.json")
    json.dump([{"begin_lat": 35.2, "begin_lon": -97.5, "ef_rating": "EF5"}], open(ep, "w"))
    assert storm_events.read_storm_events(ep)[0]["ef_rating"] == "EF5"


# ---- 7,8 interpolation (+ no silent extrapolation) ------------------------------
def test_interpolation_horizontal_vertical_and_clamp():
    x = np.linspace(-1, 1, 5); y = np.linspace(-1, 1, 5)
    f = np.add.outer(y, x)                                # (y,x), value = y + x
    XT, YT = np.meshgrid([0.0, 2.0], [0.0], indexing="xy")   # x=2 is outside -> clamped to 1
    out = interpolate.horizontal_interp(f, x, y, XT, YT)
    assert out[0, 0] == pytest.approx(0.0) and out[0, 1] == pytest.approx(1.0)  # clamped, not extrapolated
    col = np.array([1.0, 2.0, 3.0]); zs = np.array([0.0, 1.0, 2.0])
    assert interpolate.vertical_remap(col, zs, np.array([0.5]))[0] == pytest.approx(1.5)


# ---- pressure -> geometric height (task-2 refinement) ---------------------------
def test_pressure_to_height_conversion():
    p = np.array([1.0e5, 5.0e4]); Tv = np.array([288.0, 265.0])
    z = thermo.hypsometric_height(p, Tv)
    assert 5000 < z[1] < 6000                              # 1000->500 hPa ~ 5.5 km
    from atmospheric_data.sources._common import to_height_levels
    st = synthetic.synthetic_pressure_level_state(ad.CaseConfig())
    st2 = to_height_levels(st)
    zz = np.asarray(st2.ds["z"].values)
    assert np.all(np.diff(zz) > 0) and zz[0] == pytest.approx(0.0) and zz[-1] > 10000
    assert "hypsometric" in st2.provenance["vertical_conversion"]


# ---- 9 base state + decomposition -----------------------------------------------
def test_basestate_and_perturbation_decomposition():
    cfg = ad.CaseConfig()
    st = synthetic.synthetic_atmosphere(cfg)
    xm = np.linspace(-40000, 40000, 12); zm = np.linspace(0, 15000, 20)
    fields = interpolate.regrid_to_model(st, xm, xm, zm, conservative=True)
    base = basestate.base_state_from_fields(fields, zm)
    assert np.all(np.diff(base.rho0) < 0) and base.theta0.shape == (20,)
    pert = basestate.decompose_perturbations(fields, base)
    assert pert["theta_prime"].shape == (20, 12, 12) and np.isfinite(pert["w"]).all()


# ---- 11,12 IC/BC generation -----------------------------------------------------
def test_ic_bc_generation_writes_cf_netcdf():
    cfg = ad.CaseConfig()
    st = synthetic.synthetic_atmosphere(cfg)
    xm = np.linspace(-40000, 40000, 12); zm = np.linspace(0, 15000, 20)
    fields = interpolate.regrid_to_model(st, xm, xm, zm)
    tmp = tempfile.mkdtemp()
    ic_bc.write_initial_conditions(fields, xm, xm, zm, os.path.join(tmp, "ic.nc"))
    bcs = ic_bc.write_boundaries(fields, xm, xm, zm, st.ds["time"].values, tmp)
    assert set(bcs) == {"west", "east", "south", "north", "top"}
    import xarray as xr
    assert "theta" in xr.open_dataset(os.path.join(tmp, "ic.nc"))


# ---- 10 radial-velocity observation operator ------------------------------------
def test_radial_velocity_operator_geometry():
    # a purely eastward wind: Vr = +u east of radar, -u west of radar, 0 north/south
    vr = radial.project_to_radial(u=10.0, v=0.0, w=0.0, gate_x=1000.0, gate_y=0.0, gate_z=0.0,
                                  radar_xyz=(0.0, 0.0, 0.0))
    assert vr == pytest.approx(10.0)
    vr_w = radial.project_to_radial(10.0, 0.0, 0.0, -1000.0, 0.0, 0.0, (0.0, 0.0, 0.0))
    assert vr_w == pytest.approx(-10.0)
    vr_n = radial.project_to_radial(10.0, 0.0, 0.0, 0.0, 1000.0, 0.0, (0.0, 0.0, 0.0))
    assert abs(vr_n) < 1e-9


# ---- 6 validation metrics -------------------------------------------------------
def test_validation_metrics():
    o = np.array([1.0, 2.0, 3.0, 4.0]); s = o + 1.0
    assert validation.bias(o, s) == pytest.approx(1.0)
    assert validation.rmse(o, s) == pytest.approx(1.0)
    assert validation.correlation(o, s) == pytest.approx(1.0)
    assert validation.critical_success_index(o, s, 2.5) >= 0.0
    field = np.zeros((10, 10)); field[5, 5] = 10.0
    assert validation.displacement_error(field, np.roll(field, 2, 0), 5.0, spacing=1.0) == pytest.approx(2.0)


# ---- quality control ------------------------------------------------------------
def test_quality_control_passes_clean_and_flags_bad():
    cfg = ad.CaseConfig()
    st = synthetic.synthetic_atmosphere(cfg)
    rep = qc.quality_control(st)
    assert rep["summary"]["ok"] and rep["summary"]["failed"] == 0
    st.ds["T"].values[:] = 999.0                          # unphysical temperature
    rep2 = qc.quality_control(st)
    assert not rep2["summary"]["ok"]
    tmp = tempfile.mkdtemp()
    jp, mp = qc.write_reports(rep, os.path.join(tmp, "q.json"), os.path.join(tmp, "q.md"))
    assert os.path.exists(jp) and os.path.exists(mp)


# ---- 13,15 driver on CPU, offline ----------------------------------------------
def test_driver_preprocess_build_run_compare_cpu_offline():
    from atmospheric_data import driver
    cfg = ad.CaseConfig(); cfg.offline = True
    cfg.domain.width_km = 120.0; cfg.domain.height_km = 15.0
    cfg.model.parent_dx_m = 6000.0; cfg.model.execution_backend = "cpu"
    cache = ad.Cache(tempfile.mkdtemp(), offline=True); out = tempfile.mkdtemp()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pre = driver.preprocess(cfg, cache, out, logger=lambda *a: None, max_n=18)
        assert pre["qc"]["summary"]["ok"]
        assert os.path.exists(pre["initial_conditions"]) and len(pre["boundaries"]) == 5
        sim = driver.run_case(cfg, pre, steps=3, logger=lambda *a: None)
        assert sim.grid.backend.name == "cpu"
        to = sim.grid.backend.to_cpu
        assert np.isfinite(np.asarray(to(sim.state.theta))).all()
        res = driver.compare_radar(cfg, pre, cache, sim=sim, logger=lambda *a: None)
        assert "rmse" in res["metrics"]["radial_velocity"]


# ---- real IC -> multi-level AMR nest cascade ------------------------------------
def test_real_case_multilevel_cascade_offline():
    from atmospheric_data import driver
    cfg = ad.CaseConfig(); cfg.offline = True
    cfg.domain.width_km = 90.0; cfg.domain.height_km = 15.0
    cfg.model.parent_dx_m = 6000.0; cfg.model.nest_dx_m = 2000.0; cfg.model.fine_dx_m = 1000.0
    cfg.model.execution_backend = "cpu"
    cache = ad.Cache(tempfile.mkdtemp(), offline=True); out = tempfile.mkdtemp()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pre = driver.preprocess(cfg, cache, out, logger=lambda *a: None, max_n=16)
        sims, rep = driver.run_multilevel_real_case(cfg, pre, logger=lambda *a: None)
    dxs = [s.grid.dx for s in sims]
    assert dxs[0] > dxs[1] > dxs[2]                        # each level strictly finer
    assert sims[-1].grid.dx < 1500.0                       # reached the fine level from real IC
    assert np.isfinite(rep["rotation"]["zeta_abs_max"])


# ---- 14 backend auto falls back to CPU when no GPU -----------------------------
def test_execution_backend_auto_falls_back_to_cpu():
    from atmospheric_data import driver
    cfg = ad.CaseConfig(); cfg.domain.width_km = 60.0; cfg.model.parent_dx_m = 6000.0
    cfg.model.execution_backend = "auto"
    g = driver._make_grid(cfg, max_n=12)
    assert g.backend.name in ("cpu", "gpu")               # never crashes; CPU if no GPU


# ---- 17 corrupted/incomplete file handling -------------------------------------
def test_corrupted_input_raises_cleanly():
    tmp = tempfile.mkdtemp()
    bad = os.path.join(tmp, "bad.csv")
    with open(bad, "w") as f:
        f.write("not,a,valid\nsounding,file,\x00\x01")
    with pytest.raises(Exception):
        sounding.read_sounding(bad)


# ---- optional-dependency guards (1,3 real readers) -----------------------------
def test_optional_dependency_guard_raises_sourceunavailable():
    with pytest.raises(SourceUnavailable):
        require("a_module_that_does_not_exist_123", "test source", "hint")


@pytest.mark.skipif(not available(["cfgrib"]), reason="cfgrib (HRRR GRIB2) not installed")
def test_hrrr_reader_available():                          # pragma: no cover (dep-gated)
    from atmospheric_data.sources import hrrr
    assert hrrr.available()


@pytest.mark.skipif(not available(["pyart"]), reason="Py-ART (NEXRAD) not installed")
def test_nexrad_reader_available():                        # pragma: no cover (dep-gated)
    from atmospheric_data.sources import nexrad
    assert nexrad.available()
