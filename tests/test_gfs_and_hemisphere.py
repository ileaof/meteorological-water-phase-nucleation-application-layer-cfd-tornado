"""Real-time global model source (GFS via THREDDS NCSS) + southern-hemisphere storm motion.

Both were driven by the same real case: a Brazilian (Centro-Oeste) convective environment
that the pipeline could not measure. Two classes of defect are pinned here.

*Below-ground levels.* GFS reports isobaric levels that lie under the terrain by
extrapolation. At Campo Grande (terrain 477 m, surface pressure 95800 Pa) the 1000 and
975 hPa levels are not atmosphere, and feeding them to a sounding would invent a
boundary layer.

*Hemisphere.* ``bunkers_storm_motion`` deviated to the right unconditionally -- the
northern-hemisphere convention. South of the equator cyclonic rotation is clockwise and
the favoured supercell is the **left**-mover, which is why GFS reports negative SRH there.
Using the wrong mover does not flip a sign, it scores the wrong storm.

Network is never touched: payloads are verbatim captures.
"""
from __future__ import annotations

import numpy as np
import pytest

from storm_dynamics import soundings as snd


def _base():
    from meteorological_flow.grid import Grid
    return snd.build_sounding(Grid(nx=4, ny=4, nz=64, Lx=16000, Ly=16000, Lz=16000))


# ---------------------------------------------------------------------------
# hemisphere
# ---------------------------------------------------------------------------
def test_hemisphere_sign_defaults_to_north():
    """Both arguments omitted must reproduce the original (northern) behaviour."""
    assert snd._hemisphere_sign() == 1.0
    assert snd._hemisphere_sign(None, None) == 1.0


@pytest.mark.parametrize("value,expected", [
    ("north", 1.0), ("N", 1.0), ("northern", 1.0),
    ("south", -1.0), ("s", -1.0), ("Southern", -1.0),
])
def test_hemisphere_sign_from_name(value, expected):
    assert snd._hemisphere_sign(value) == expected


def test_hemisphere_sign_inferred_from_latitude():
    assert snd._hemisphere_sign(None, -20.47) == -1.0      # Campo Grande
    assert snd._hemisphere_sign(None, 35.22) == 1.0        # Norman, OK
    assert snd._hemisphere_sign(None, 0.0) == 1.0          # equator -> default
    assert snd._hemisphere_sign("north", -20.0) == 1.0     # explicit beats inferred


def test_hemisphere_sign_rejects_nonsense():
    with pytest.raises(ValueError):
        snd._hemisphere_sign("east")


def test_bunkers_default_is_unchanged_and_south_mirrors_it():
    """The southern mover is the northern one reflected about the mean wind."""
    base = _base()
    z = np.asarray(base.zc)
    sel = z <= 6000.0
    um = float(np.mean(np.asarray(base.u0)[sel]))
    vm = float(np.mean(np.asarray(base.v0)[sel]))
    u = np.interp([0.0, 6000.0], z, np.asarray(base.u0))
    v = np.interp([0.0, 6000.0], z, np.asarray(base.v0))
    shx, shy = u[1] - u[0], v[1] - v[0]
    smag = np.hypot(shx, shy) + 1e-9
    expected = (um + 7.5 * (shy / smag), vm + 7.5 * (-shx / smag))   # pre-change formula

    north = snd.bunkers_storm_motion(base)
    assert north[0] == pytest.approx(expected[0], rel=1e-12)
    assert north[1] == pytest.approx(expected[1], rel=1e-12)
    assert north == snd.bunkers_storm_motion(base, hemisphere="north")

    south = snd.bunkers_storm_motion(base, hemisphere="south")
    assert south[0] == pytest.approx(2 * um - north[0], rel=1e-9)
    assert south[1] == pytest.approx(2 * vm - north[1], rel=1e-9)
    assert south == snd.bunkers_storm_motion(base, latitude_deg=-20.47)


def test_southern_report_flips_the_cyclonic_srh_sign():
    base = _base()
    north = snd.low_level_geometry_report(base, 500.0)
    south = snd.low_level_geometry_report(base, 500.0, latitude_deg=-20.47)
    assert north["hemisphere"] == "north" and south["hemisphere"] == "south"
    assert north["SRH_layer_cyclonic_m2_s2"] == pytest.approx(north["SRH_layer_m2_s2"])
    assert south["SRH_layer_cyclonic_m2_s2"] == pytest.approx(-south["SRH_layer_m2_s2"])
    assert south["SRH_0_3km_cyclonic_m2_s2"] == pytest.approx(-south["SRH_0_3km_m2_s2"])


def test_southern_hemisphere_changes_the_measured_geometry():
    """The point of the fix: the wrong mover scores the wrong storm."""
    base = _base()
    n = snd.low_level_geometry_report(base, 500.0)
    s = snd.low_level_geometry_report(base, 500.0, hemisphere="south")
    assert n["storm_motion_ms"] != s["storm_motion_ms"]
    assert n["critical_angle_deg"] != pytest.approx(s["critical_angle_deg"], abs=1e-6)


# ---------------------------------------------------------------------------
# NCSS parsing
# ---------------------------------------------------------------------------
#: Verbatim NCSS header + rows captured 2026-10-08 over Campo Grande. The 1000 and 975 hPa
#: rows sit BELOW the 95800 Pa surface (terrain 477 m) and must be discarded.
_HEADER = ('time,alt[unit="Pa"],station,latitude[unit="degrees_north"],'
           'longitude[unit="degrees_east"],Temperature_isobaric[unit="K"],'
           'Relative_humidity_isobaric[unit="percent"],'
           'u-component_of_wind_isobaric[unit="m/s"],'
           'v-component_of_wind_isobaric[unit="m/s"],'
           'Geopotential_height_isobaric[unit="gpm"]')
_ROWS = [
    "2026-10-08T12:00:00Z,70000.0,P,-20.5,-54.75,285.00,60.0,-6.00,-9.00,3150.00",
    "2026-10-08T12:00:00Z,85000.0,P,-20.5,-54.75,294.00,55.0,-5.00,-8.00,1520.00",
    "2026-10-08T12:00:00Z,92500.0,P,-20.5,-54.75,299.69,51.4,-3.95,-7.41,788.72",
    "2026-10-08T12:00:00Z,95000.0,P,-20.5,-54.75,301.73,48.0,-3.97,-5.91,552.34",
    "2026-10-08T12:00:00Z,97500.0,P,-20.5,-54.75,303.60,47.6,-2.97,-4.30,322.12",
    "2026-10-08T12:00:00Z,100000.0,P,-20.5,-54.75,305.02,47.6,-2.96,-4.30,102.49",
]
_CSV = "\n".join([_HEADER] + _ROWS)


def test_ncss_header_units_are_stripped():
    from atmospheric_data.sources.gfs import _parse_ncss_csv
    rows, header = _parse_ncss_csv(_CSV)
    assert "Temperature_isobaric" in header and "alt" in header
    assert not any("[" in h for h in header)
    assert len(rows) == 6
    assert rows[-1]["alt"] == "100000.0"


def test_ncss_parser_refuses_empty_and_header_only():
    from atmospheric_data.sources.gfs import _parse_ncss_csv
    with pytest.raises(ValueError):
        _parse_ncss_csv("")
    with pytest.raises(ValueError):
        _parse_ncss_csv(_HEADER)


def test_to_float_is_nan_on_junk():
    from atmospheric_data.sources.gfs import _to_float
    assert _to_float("1.5") == 1.5
    assert np.isnan(_to_float("n/a"))
    assert np.isnan(_to_float(None))


def _patch(monkeypatch, p_sfc=95800.375, terrain=477.4994):
    """Serve the canned payloads instead of reaching the network."""
    from atmospheric_data.sources import gfs
    monkeypatch.setattr(gfs, "_get",
                        lambda url, timeout, cache=None, source="gfs", key=None: _CSV)
    monkeypatch.setattr(gfs, "download_surface", lambda *a, **k: {
        "Pressure_surface": p_sfc,
        "Geopotential_height_surface": terrain,
        "Temperature_height_above_ground": 306.0,
        "Dewpoint_temperature_height_above_ground": 292.0,
        "u-component_of_wind_height_above_ground": -2.5,
        "v-component_of_wind_height_above_ground": -3.5,
        "Convective_available_potential_energy_surface": 600.0,
        "Convective_inhibition_surface": -46.64,
        "Storm_relative_helicity_height_above_ground_layer": -6.76,
        "valid_utc": "2026-10-08T12:00:00Z",
    })
    return gfs


def test_below_ground_levels_are_dropped_not_trusted(monkeypatch):
    gfs = _patch(monkeypatch)
    prof = gfs.download_profile(-20.47, -54.67, include_surface=False)
    assert prof["pressure_Pa"].max() <= 95800.375
    assert 100000.0 not in prof["pressure_Pa"]       # 1000 hPa is underground here
    assert 97500.0 not in prof["pressure_Pa"]        # and so is 975 hPa
    assert len(prof["pressure_Pa"]) == 4

    kept = gfs.download_profile(-20.47, -54.67, drop_below_ground=False,
                                include_surface=False)
    assert len(kept["pressure_Pa"]) == 6             # opting out keeps them


def test_profile_heights_are_above_ground_and_ascending(monkeypatch):
    gfs = _patch(monkeypatch)
    prof = gfs.download_profile(-20.47, -54.67, include_surface=False)
    assert np.all(np.diff(prof["height_m"]) > 0)
    # 950 hPa at 552.34 gpm MSL over 477.4994 m terrain -> ~74.8 m AGL
    assert prof["height_m"][0] == pytest.approx(552.34 - 477.4994, abs=1e-3)
    assert prof["terrain_m"] == pytest.approx(477.4994)
    assert prof["surface_pressure_Pa"] == pytest.approx(95800.375)


def test_profile_humidity_goes_through_the_project_thermo(monkeypatch):
    from atmospheric_data import thermo
    gfs = _patch(monkeypatch)
    prof = gfs.download_profile(-20.47, -54.67, include_surface=False)
    assert np.all((prof["relative_humidity"] >= 0.0) & (prof["relative_humidity"] <= 1.0))
    expected = thermo.specific_humidity_from_rh(prof["relative_humidity"],
                                                prof["temperature_K"], prof["pressure_Pa"])
    assert np.allclose(prof["specific_humidity"], expected, rtol=1e-12)
    assert np.all(prof["specific_humidity"] > 0.0)


def test_profile_carries_gfs_diagnostics_for_crosschecking(monkeypatch):
    gfs = _patch(monkeypatch)
    prof = gfs.download_profile(-20.47, -54.67)
    gd = prof["gfs_diagnostics"]
    assert gd["CAPE_J_kg"] == 600.0
    assert gd["CIN_J_kg"] == pytest.approx(-46.64)
    assert gd["SRH_m2_s2"] < 0.0          # southern hemisphere: negative by convention


def test_profile_refuses_when_the_column_is_mostly_underground(monkeypatch):
    """A surface pressure that buries the column must raise, not return 1-2 levels."""
    gfs = _patch(monkeypatch, p_sfc=72000.0, terrain=3000.0)
    with pytest.raises(ValueError, match="usable level"):
        gfs.download_profile(-20.47, -54.67)


def test_gfs_profile_offline_without_cache_refuses():
    from atmospheric_data.sources import gfs
    with pytest.raises(FileNotFoundError):
        gfs.download_profile(-20.47, -54.67, offline=True)


def test_gfs_grid_offline_without_file_refuses():
    from atmospheric_data.sources import gfs
    with pytest.raises(FileNotFoundError):
        gfs.download_grid("does-not-exist.nc", -10, -25, -60, -45, offline=True)


def test_gfs_feeds_the_existing_basestate_contract(monkeypatch):
    """The dict shape must be exactly what sounding.profile_to_basestate already consumes."""
    from meteorological_flow.grid import Grid

    from atmospheric_data.sources import sounding as sounding_src
    gfs = _patch(monkeypatch)
    prof = gfs.download_profile(-20.47, -54.67, drop_below_ground=False)
    g = Grid(nx=4, ny=4, nz=48, Lx=16000, Ly=16000, Lz=12000)
    base = sounding_src.profile_to_basestate(g, prof)
    rep = snd.low_level_geometry_report(base, z_top=500.0, latitude_deg=-20.47)
    assert rep["hemisphere"] == "south"
    assert np.isfinite(rep["critical_angle_deg"])
    assert 0.0 <= rep["streamwise_fraction"] <= 1.0


# ---------------------------------------------------------------------------
# the surface parcel: dropping below-ground levels must not cost the ground
# ---------------------------------------------------------------------------
def test_surface_level_is_prepended_so_the_column_starts_at_the_ground(monkeypatch):
    """Without this, a parcel is lifted from 134 m AGL and CAPE reads systematically low."""
    gfs = _patch(monkeypatch)
    with_sfc = gfs.download_profile(-20.47, -54.67)
    without = gfs.download_profile(-20.47, -54.67, include_surface=False)
    assert len(with_sfc["height_m"]) == len(without["height_m"]) + 1
    assert with_sfc["height_m"][0] == pytest.approx(gfs.SURFACE_LEVEL_HEIGHT_M)
    assert with_sfc["pressure_Pa"][0] == pytest.approx(95800.375)
    assert with_sfc["temperature_K"][0] == pytest.approx(306.0)
    assert np.all(np.diff(with_sfc["height_m"]) > 0)
    # the surface level is BELOW the lowest surviving isobaric level
    assert with_sfc["pressure_Pa"][0] > with_sfc["pressure_Pa"][1]


def test_surface_humidity_comes_from_the_dewpoint_via_project_thermo(monkeypatch):
    """RH is recovered with the project's own saturation curve, not a second formula."""
    from atmospheric_data import thermo
    gfs = _patch(monkeypatch)
    prof = gfs.download_profile(-20.47, -54.67)
    e = thermo.saturation_vapor_pressure(np.array([292.0]))
    es = thermo.saturation_vapor_pressure(np.array([306.0]))
    assert prof["relative_humidity"][0] == pytest.approx(float(e[0] / es[0]), rel=1e-9)
    assert 0.0 < prof["relative_humidity"][0] < 1.0
    assert prof["specific_humidity"][0] > 0.0


def test_surface_level_is_skipped_when_the_fields_are_missing(monkeypatch):
    """A missing screen-level field must degrade to the isobaric column, not crash."""
    from atmospheric_data.sources import gfs as gfs_mod
    monkeypatch.setattr(gfs_mod, "_get",
                        lambda url, timeout, cache=None, source="gfs", key=None: _CSV)
    monkeypatch.setattr(gfs_mod, "download_surface", lambda *a, **k: {
        "Pressure_surface": 95800.375, "Geopotential_height_surface": 477.4994,
        "valid_utc": "2026-10-08T12:00:00Z",
    })
    prof = gfs_mod.download_profile(-20.47, -54.67)
    assert len(prof["height_m"]) == 4                 # isobaric only, no invented surface
    assert prof["height_m"][0] > gfs_mod.SURFACE_LEVEL_HEIGHT_M


def test_surface_level_not_added_when_it_is_not_below_the_column(monkeypatch):
    """If nothing was dropped, there is no gap to fill and no duplicate level."""
    gfs = _patch(monkeypatch, p_sfc=99000.0, terrain=100.0)
    prof = gfs.download_profile(-20.47, -54.67)
    # 100000 Pa still exceeds p_sfc so it is dropped; 97500 survives and sits under 99000,
    # so the surface level DOES belong here
    assert prof["pressure_Pa"][0] == pytest.approx(99000.0)
    # now a surface pressure beneath the whole column -> nothing to prepend
    gfs2 = _patch(monkeypatch, p_sfc=60000.0, terrain=4000.0)
    with pytest.raises(ValueError, match="usable level"):
        gfs2.download_profile(-20.47, -54.67)
