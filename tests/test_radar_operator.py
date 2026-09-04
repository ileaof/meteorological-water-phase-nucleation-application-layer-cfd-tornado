"""Radar forward operator -- validated against analytic vortices with KNOWN answers.

The Moore 2013 target, V_rot = 26 m/s, is a beam-averaged Doppler observable, not a wind speed.
Comparing a model grid-point wind against it is a category error.  These tests pin the operator
that makes the comparison well posed, and pin the physical behaviour it must reproduce: a radar
UNDER-READS a vortex smaller than its beam, monotonically, converging to the truth only when the
core is much larger than the beam.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from atmospheric_data import radar_operator as ro


# ------------------------------------------------------------------------- beam geometry
def test_beam_height_matches_wsr88d_reference_values():
    """0.5 deg under 4/3-earth refraction: ~315 m AGL at 30 km is the standard textbook value."""
    assert ro.beam_height_m(30000.0, 0.5) == pytest.approx(315.0, abs=15.0)
    assert ro.beam_height_m(10000.0, 0.5) == pytest.approx(93.0, abs=10.0)
    # monotone in range and in elevation
    assert ro.beam_height_m(40000.0, 0.5) > ro.beam_height_m(30000.0, 0.5)
    assert ro.beam_height_m(30000.0, 1.5) > ro.beam_height_m(30000.0, 0.5)


def test_range_for_beam_height_inverts_beam_height():
    for h in (200.0, 460.0, 1000.0):
        r = ro.range_for_beam_height(h, 0.5)
        assert ro.beam_height_m(r, 0.5) == pytest.approx(h, rel=1e-6)


def test_beam_diameter_scales_with_range():
    rad = ro.RadarSpec(beamwidth_deg=0.925)
    assert rad.beam_diameter_m(30000.0) == pytest.approx(484.0, abs=10.0)
    assert rad.beam_diameter_m(60000.0) == pytest.approx(2 * rad.beam_diameter_m(30000.0), rel=1e-9)
    arr = rad.beam_diameter_m(np.array([10000.0, 20000.0]))
    assert arr.shape == (2,)                       # array in, array out (the sweep needs this)


# ------------------------------------------------------------------------- analytic vortex
def _sweep_vrot(core_m, v_max=26.0, range_m=None, n_beam=5, n_gate=3, dx=20.0):
    rad = ro.RadarSpec(elevation_deg=0.5)
    rng = range_m if range_m is not None else ro.range_for_beam_height(460.0, 0.5)
    cx, cy = 0.0, rng * np.cos(np.radians(0.5))
    half = max(2500.0, 3.0 * core_m)
    x = np.arange(cx - half, cx + half + 1, dx)
    y = np.arange(cy - half, cy + half + 1, dx)
    z = np.array([300.0, 460.0, 620.0])
    u, v, w = ro.rankine_vortex(x, y, z, (cx, cy), v_max, core_m)
    az, rg = ro.sweep_grid_for_domain(rad, cx, cy, half - 200.0, az_resolution_deg=0.25)
    sw = ro.synthetic_sweep(u, v, w, x, y, z, rad, az, rg, n_beam=n_beam, n_gate=n_gate)
    return ro.vrot_from_sweep(sw, max_separation_m=4.0 * half), rad.beam_diameter_m(rng)


def test_rankine_vortex_has_the_expected_analytic_structure():
    x = np.arange(-1000.0, 1001.0, 25.0); y = x.copy(); z = np.array([0.0])
    u, v, _ = ro.rankine_vortex(x, y, z, (0.0, 0.0), 26.0, 250.0)
    speed = np.sqrt(u ** 2 + v ** 2)[:, :, 0]
    assert speed.max() == pytest.approx(26.0, rel=0.05)      # peak equals V_max
    c = len(x) // 2
    assert speed[c, c] < 3.0                                  # ~zero at the centre


def test_operator_converges_to_truth_when_the_core_is_much_larger_than_the_beam():
    """The essential correctness property: a well-resolved vortex must be read accurately."""
    got, D = _sweep_vrot(core_m=2000.0, v_max=26.0)
    assert 2 * 2000.0 / D > 4.0                               # core really is >> beam
    assert got["v_rot_m_s"] == pytest.approx(26.0, rel=0.15)


def test_operator_under_reads_a_sub_beam_core_and_does_so_monotonically():
    """A vortex smaller than the beam MUST be under-reported, more so as it shrinks.  This is
    the physical reason a 26 m/s observation does not mean a 26 m/s wind."""
    fracs = {}
    for R in (125.0, 250.0, 500.0, 1000.0, 2000.0):
        got, _ = _sweep_vrot(core_m=R, v_max=26.0)
        fracs[R] = got["v_rot_m_s"] / 26.0
    seq = [fracs[R] for R in (125.0, 250.0, 500.0, 1000.0, 2000.0)]
    assert all(b > a for a, b in zip(seq, seq[1:])), fracs     # strictly monotone
    assert fracs[125.0] < 0.5                                  # a 125 m core is badly under-read
    assert fracs[2000.0] > 0.8                                 # a large core is read well


def test_beam_averaging_reads_lower_than_point_sampling():
    """Volume averaging is what causes the under-reading; point sampling must read higher."""
    pt, _ = _sweep_vrot(core_m=125.0, n_beam=1, n_gate=1)
    bm, _ = _sweep_vrot(core_m=125.0, n_beam=5, n_gate=3)
    assert bm["v_rot_m_s"] < pt["v_rot_m_s"]


def test_under_reading_worsens_with_range_because_the_beam_widens():
    near, _ = _sweep_vrot(core_m=125.0, range_m=10000.0)
    far, _ = _sweep_vrot(core_m=125.0, range_m=40000.0)
    assert far["v_rot_m_s"] < near["v_rot_m_s"]


def test_vrot_rejects_extrema_that_are_not_a_couplet():
    """max and min far apart are two unrelated features, not a couplet -- must be refused."""
    got, _ = _sweep_vrot(core_m=125.0)
    tight = ro.vrot_from_sweep({"v_r": np.array([[-10.0, 10.0]]),
                                "x_m": np.array([[0.0, 50000.0]]),
                                "y_m": np.array([[0.0, 0.0]]),
                                "beam_diameter_m": np.array([[500.0, 500.0]]),
                                "height_m": np.array([[460.0, 460.0]])},
                               max_separation_m=4000.0)
    assert tight["valid"] is False and np.isnan(tight["v_rot_m_s"])
    assert "not a couplet" in tight["reason"]


def test_sweep_reports_its_own_geometry_provenance():
    got, _ = _sweep_vrot(core_m=500.0)
    for key in ("couplet_separation_m", "beam_diameter_m", "sample_height_m", "valid"):
        assert key in got, key
    assert got["sample_height_m"] > 0.0 and got["beam_diameter_m"] > 0.0
