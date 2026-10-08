"""Low-level hodograph geometry + the data plumbing that feeds it.

Covers the two additions made for the real-case (Brazil) path:

* ``storm_dynamics.soundings``: configurable-layer SRH, the Esterheld & Giuliano
  critical angle, and the streamwise vorticity fraction -- the sounding-side proxy for
  the geometry the A-L study isolated as the tornadogenesis bottleneck.
* ``atmospheric_data.sources``: TLS trust resolution (a TLS-inspecting middlebox broke
  every ``requests`` download while ``curl`` worked) and the live global METAR fetch.

The geometry tests are analytic -- hodographs whose answer is known by construction --
so they pin the sign and normalisation conventions, not just "it runs".
"""
from __future__ import annotations

import os
from types import SimpleNamespace

import numpy as np
import pytest

from storm_dynamics import soundings as snd


def _hodo(u0, v0, z=None):
    """Minimal base-state stand-in: the geometry functions read only zc/u0/v0."""
    z = np.linspace(0.0, 6000.0, 61) if z is None else np.asarray(z, float)
    return SimpleNamespace(zc=z, u0=np.asarray(u0, float), v0=np.asarray(v0, float))


def _linear_u(z, rate=0.02):
    """Unidirectional shear along +x: shear vector is (+, 0), so omega_h is (0, +)."""
    return _hodo(rate * z, np.zeros_like(z), z)


# ---------------------------------------------------------------------------
# critical angle / streamwise fraction -- analytic
# ---------------------------------------------------------------------------
def test_perpendicular_inflow_is_purely_streamwise():
    """Storm-relative inflow along +y against +x shear: omega_h is along +y too."""
    z = np.linspace(0.0, 6000.0, 61)
    base = _linear_u(z)
    # surface wind is (0,0), so storm_motion (0,-10) gives sr_wind (0,+10) || omega_h
    assert snd.critical_angle(base, 500.0, (0.0, -10.0)) == pytest.approx(90.0, abs=1e-6)
    assert snd.streamwise_vorticity_fraction(base, 500.0, (0.0, -10.0)) == pytest.approx(1.0, abs=1e-9)


def test_parallel_inflow_is_purely_crosswise():
    """Inflow along the shear vector leaves vorticity entirely crosswise."""
    z = np.linspace(0.0, 6000.0, 61)
    base = _linear_u(z)
    assert snd.critical_angle(base, 500.0, (-10.0, 0.0)) == pytest.approx(0.0, abs=1e-6)
    assert snd.streamwise_vorticity_fraction(base, 500.0, (-10.0, 0.0)) == pytest.approx(0.0, abs=1e-9)


def test_antiparallel_inflow_is_also_crosswise():
    """180 deg is as crosswise as 0 deg -- the fraction is a magnitude."""
    z = np.linspace(0.0, 6000.0, 61)
    base = _linear_u(z)
    assert snd.critical_angle(base, 500.0, (10.0, 0.0)) == pytest.approx(180.0, abs=1e-6)
    assert snd.streamwise_vorticity_fraction(base, 500.0, (10.0, 0.0)) == pytest.approx(0.0, abs=1e-9)


@pytest.mark.parametrize("theta_deg", [10.0, 35.0, 60.0, 90.0, 125.0, 160.0])
def test_fraction_equals_sine_of_critical_angle(theta_deg):
    """The two diagnostics are one quantity: fraction == |sin(critical angle)|.

    omega_h is the shear vector rotated 90 deg, so the angle it makes with the inflow is
    the critical angle minus 90 deg.  Holding that identity is what makes the fraction
    comparable with vorticity_budget.tilting_efficiency's streamwise share.
    """
    z = np.linspace(0.0, 6000.0, 61)
    base = _linear_u(z)                     # shear along +x
    th = np.radians(theta_deg)
    # storm motion chosen so the storm-relative surface wind sits at theta from +x
    sr = np.array([np.cos(th), np.sin(th)]) * 12.0
    motion = (-sr[0], -sr[1])               # surface wind is (0,0)
    ca = snd.critical_angle(base, 500.0, motion)
    fr = snd.streamwise_vorticity_fraction(base, 500.0, motion)
    assert ca == pytest.approx(theta_deg, abs=1e-6)
    assert fr == pytest.approx(abs(np.sin(th)), abs=1e-9)


def test_fraction_is_bounded():
    z = np.linspace(0.0, 6000.0, 61)
    base = _linear_u(z)
    for motion in [(3.0, -7.0), (-11.0, 2.0), (0.5, 0.5), (-4.0, -9.0)]:
        fr = snd.streamwise_vorticity_fraction(base, 500.0, motion)
        assert 0.0 - 1e-12 <= fr <= 1.0 + 1e-12


def test_degenerate_hodographs_give_nan_not_a_number():
    """Zero shear or zero storm-relative wind must not fake a geometry."""
    z = np.linspace(0.0, 6000.0, 61)
    calm = _hodo(np.zeros_like(z), np.zeros_like(z), z)
    assert np.isnan(snd.critical_angle(calm, 500.0, (5.0, 5.0)))
    assert np.isnan(snd.streamwise_vorticity_fraction(calm, 500.0, (5.0, 5.0)))
    base = _linear_u(z)
    # storm motion equal to the surface wind -> zero storm-relative inflow
    assert np.isnan(snd.critical_angle(base, 500.0, (0.0, 0.0)))
    assert np.isnan(snd.streamwise_vorticity_fraction(base, 500.0, (0.0, 0.0)))


def test_veering_hodograph_is_more_streamwise_than_straight_one():
    """The physics check: curvature is what buys streamwise vorticity."""
    z = np.linspace(0.0, 6000.0, 121)
    straight = _hodo(0.004 * z, np.zeros_like(z), z)
    turn = np.clip(z / 1000.0, 0.0, 1.0) * (np.pi / 2.0)
    veering = _hodo(20.0 * np.sin(turn), 20.0 * (np.cos(turn) - 1.0), z)
    fr_straight = snd.streamwise_vorticity_fraction(straight, 500.0)
    fr_veering = snd.streamwise_vorticity_fraction(veering, 500.0)
    assert np.isfinite(fr_straight) and np.isfinite(fr_veering)
    assert fr_veering > fr_straight


# ---------------------------------------------------------------------------
# SRH: the new z_bot must not disturb the established default
# ---------------------------------------------------------------------------
def test_srh_default_is_byte_identical_to_the_old_selection():
    """Adding z_bot kept the default path unchanged (heights are AGL, so z>=0 always)."""
    from meteorological_flow.grid import Grid
    base = snd.build_sounding(Grid(nx=4, ny=4, nz=64, Lx=16000, Ly=16000, Lz=16000))
    motion = snd.bunkers_storm_motion(base)

    z = np.asarray(base.zc, float)
    u = np.asarray(base.u0, float); v = np.asarray(base.v0, float)
    sel = z <= 3000.0                                   # the pre-change selection
    zc, uc, vc = z[sel], u[sel], v[sel]
    integrand = (vc - motion[1]) * np.gradient(uc, zc) - (uc - motion[0]) * np.gradient(vc, zc)
    trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    expected = float(trapz(integrand, zc))

    assert snd.storm_relative_helicity(base, storm_motion=motion) == pytest.approx(expected, rel=1e-12)
    assert (snd.storm_relative_helicity(base, storm_motion=motion)
            == snd.storm_relative_helicity(base, storm_motion=motion, z_bot=0.0))


def test_shallow_srh_layer_is_a_subset_of_the_deep_one():
    from meteorological_flow.grid import Grid
    base = snd.build_sounding(Grid(nx=4, ny=4, nz=64, Lx=16000, Ly=16000, Lz=16000))
    motion = snd.bunkers_storm_motion(base)
    srh_500 = snd.storm_relative_helicity(base, z_top=500.0, storm_motion=motion)
    srh_3k = snd.storm_relative_helicity(base, z_top=3000.0, storm_motion=motion)
    assert abs(srh_500) < abs(srh_3k)       # a 500 m slice cannot exceed the 3 km one
    assert snd.storm_relative_helicity(base, z_top=3000.0, z_bot=2900.0,
                                       storm_motion=motion) != srh_3k


def test_geometry_report_keys_and_consistency():
    from meteorological_flow.grid import Grid
    base = snd.build_sounding(Grid(nx=4, ny=4, nz=64, Lx=16000, Ly=16000, Lz=16000))
    rep = snd.low_level_geometry_report(base, z_top=500.0)
    for k in ("layer_top_m", "storm_motion_ms", "critical_angle_deg", "streamwise_fraction",
              "SRH_layer_m2_s2", "SRH_0_3km_m2_s2", "shear_layer_m_s", "shear_0_6km_m_s"):
        assert k in rep
    assert rep["layer_top_m"] == 500.0
    assert rep["streamwise_fraction"] == pytest.approx(
        abs(np.sin(np.radians(rep["critical_angle_deg"]))), abs=1e-9)
    assert rep["shear_layer_m_s"] <= rep["shear_0_6km_m_s"] + 1e-9


# ---------------------------------------------------------------------------
# TLS trust resolution (no network)
# ---------------------------------------------------------------------------
def test_ca_bundle_never_disables_verification():
    """It must return a usable bundle path or True -- never False/None, which would
    silently turn off certificate checking."""
    from atmospheric_data.sources import _tls
    bundle = _tls.ca_bundle(refresh=True)
    assert bundle is not False and bundle is not None
    assert bundle is True or isinstance(bundle, str)
    if isinstance(bundle, str):
        assert os.path.exists(bundle)
        with open(bundle, "r", encoding="utf-8") as f:
            assert "BEGIN CERTIFICATE" in f.read(200000)


def test_explicit_env_bundle_wins(tmp_path, monkeypatch):
    from atmospheric_data.sources import _tls
    pem = tmp_path / "custom.pem"
    pem.write_text("-----BEGIN CERTIFICATE-----\nnot-a-real-cert\n-----END CERTIFICATE-----\n")
    monkeypatch.setenv("MET_H2O_CA_BUNDLE", str(pem))
    assert _tls.ca_bundle() == str(pem)


def test_extra_ca_dir_is_merged(tmp_path, monkeypatch):
    from atmospheric_data.sources import _tls
    extra = tmp_path / "extra"; extra.mkdir()
    marker = "-----BEGIN CERTIFICATE-----\nMARKER-EXTRA-CA\n-----END CERTIFICATE-----"
    (extra / "site.pem").write_text(marker + "\n")
    monkeypatch.delenv("MET_H2O_CA_BUNDLE", raising=False)
    monkeypatch.setenv("MET_H2O_EXTRA_CA_DIR", str(extra))
    bundle = _tls.ca_bundle(refresh=True, cache_dir=str(tmp_path / "cache"))
    assert isinstance(bundle, str)
    with open(bundle, "r", encoding="utf-8") as f:
        assert "MARKER-EXTRA-CA" in f.read()


def test_describe_reports_the_resolution():
    from atmospheric_data.sources import _tls
    d = _tls.describe()
    for k in ("bundle", "certifi", "os_roots_found", "enum_certificates"):
        assert k in d
    assert isinstance(d["os_roots_found"], int) and d["os_roots_found"] >= 0


# ---------------------------------------------------------------------------
# live METAR parsing (no network -- a real captured report)
# ---------------------------------------------------------------------------
#: A genuine report, kept verbatim: SBBR carries wdir "VRB" and visib "6+", the two
#: non-numeric sentinels that a naive float() would crash on.
_SBBR = {
    "icaoId": "SBBR", "reportTime": "2026-10-08T10:00:00.000Z",
    "temp": 20, "dewp": 14, "wdir": "VRB", "wspd": 2, "visib": "6+", "altim": 1020,
    "lat": -15.867, "lon": -47.933, "elev": 1051, "name": "Brasilia Intl, DF, BR",
    "rawOb": "METAR SBBR 081000Z VRB02KT CAVOK 20/14 Q1020",
}
_SBCG = {
    "icaoId": "SBCG", "reportTime": "2026-10-08T10:00:00.000Z",
    "temp": 25, "dewp": 19, "wdir": 50, "wspd": 8, "altim": 1016,
    "lat": -20.47, "lon": -54.67, "elev": 556, "rawOb": "METAR SBCG 081000Z 05008KT 25/19 Q1016",
}


def test_num_tolerates_the_api_sentinels():
    from atmospheric_data.sources.metar import _num
    assert _num("VRB") is None
    assert _num("6+") is None
    assert _num(None) is None
    assert _num(True) is None            # bool must not slip through as 1.0
    assert _num(50) == 50.0
    assert _num("12.5") == 12.5


def test_variable_wind_yields_no_invented_direction():
    from atmospheric_data.sources.metar import _report_to_si
    r = _report_to_si(_SBBR)
    assert r["wind_dir_deg"] is None and r["u_ms"] is None and r["v_ms"] is None
    assert r["wind_speed_ms"] == pytest.approx(2 * 0.514444, rel=1e-4)


def test_report_to_si_conversions():
    from atmospheric_data.sources.metar import _report_to_si
    r = _report_to_si(_SBCG)
    assert r["station"] == "SBCG"
    assert r["temperature_K"] == pytest.approx(298.15, abs=1e-6)
    assert r["dewpoint_K"] == pytest.approx(292.15, abs=1e-6)
    assert r["pressure_Pa"] == pytest.approx(101600.0, abs=1e-6)
    assert r["wind_speed_ms"] == pytest.approx(8 * 0.514444, rel=1e-4)
    assert r["u_ms"] is not None and r["v_ms"] is not None
    assert np.hypot(r["u_ms"], r["v_ms"]) == pytest.approx(r["wind_speed_ms"], rel=1e-6)
    assert r["elevation_m"] == 556.0


def test_lcl_estimate_tracks_dewpoint_depression():
    from atmospheric_data.sources.metar import _report_to_si, lcl_estimate_m
    assert lcl_estimate_m(None, 290.0) is None
    assert lcl_estimate_m(300.0, 300.0) == 0.0
    assert lcl_estimate_m(300.0, 294.0) == pytest.approx(750.0)
    assert lcl_estimate_m(300.0, 310.0) == 0.0      # clamped, never negative
    # SBCG: 25/19 -> 6 K depression -> ~750 m; SBBR: 20/14 -> also 6 K
    assert _report_to_si(_SBCG)["lcl_estimate_m"] == pytest.approx(750.0)
    assert _report_to_si(_SBBR)["lcl_estimate_m"] == pytest.approx(750.0)


def test_metar_offline_without_cache_refuses_rather_than_downloading():
    from atmospheric_data.sources.metar import download_metar
    with pytest.raises(FileNotFoundError):
        download_metar(("SBCG",), offline=True)


def test_metar_reads_from_cache_without_network(tmp_path):
    """A cached payload must satisfy an offline request end to end."""
    import json

    from atmospheric_data.cache import Cache
    from atmospheric_data.sources.metar import download_metar
    cache = Cache(directory=str(tmp_path), offline=True)
    key = "SBCG-SBBR_hlatest"
    with open(cache.path("metar", key, ".json"), "w", encoding="utf-8") as f:
        json.dump([_SBCG, _SBBR], f)
    reports = download_metar(("SBCG", "SBBR"), cache=cache, offline=True)
    assert [r["station"] for r in reports] == ["SBCG", "SBBR"]
    assert reports[0]["temperature_K"] == pytest.approx(298.15, abs=1e-6)
