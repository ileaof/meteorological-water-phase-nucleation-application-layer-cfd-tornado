"""Objective storm classification (storm_dynamics.classification): synthetic storms must land in
each category, and the ladder must be monotone."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from storm_dynamics import classification as cl


def _base():
    return {"w_max": 0.0, "midlevel_mesocyclone": 0.0, "updraft_helicity_2_5km": 0.0,
            "near_surface_zeta_max": 0.0, "v_theta_max_m_s": 0.0, "circulation_m2_s": 0.0,
            "pressure_deficit_Pa": 0.0, "level_m": 5000.0, "gust_front_convergence_s": 0.0}


def test_no_deep_convection():
    d = _base(); d["w_max"] = 3.0
    assert cl.classify(d)["category"] == "NO_DEEP_CONVECTION"


def test_ordinary_convection():
    d = _base(); d.update(w_max=22.0, midlevel_mesocyclone=1e-3)
    assert cl.classify(d)["category"] == "ORDINARY_CONVECTION"


def test_supercell():
    d = _base(); d.update(w_max=28.0, midlevel_mesocyclone=9e-3, updraft_helicity_2_5km=250.0)
    assert cl.classify(d)["category"] == "SUPERCELL"


def test_low_level_mesocyclone():
    d = _base(); d.update(w_max=28.0, midlevel_mesocyclone=9e-3, updraft_helicity_2_5km=250.0,
                          near_surface_zeta_max=6e-3, v_theta_max_m_s=10.0)
    assert cl.classify(d)["category"] == "LOW_LEVEL_MESOCYCLONE"


def test_tornado_like_vortex():
    d = _base(); d.update(w_max=32.0, midlevel_mesocyclone=1.2e-2, updraft_helicity_2_5km=400.0,
                          near_surface_zeta_max=1e-2, v_theta_max_m_s=22.0, level_m=800.0)
    # tornado-like intensity but aloft / no surface convergence -> not surface-connected
    assert cl.classify(d, persistence_s=200.0)["category"] == "TORNADO_LIKE_VORTEX"


def test_surface_connected():
    d = _base(); d.update(w_max=35.0, midlevel_mesocyclone=1.5e-2, updraft_helicity_2_5km=600.0,
                          near_surface_zeta_max=2e-2, v_theta_max_m_s=28.0, circulation_m2_s=6e4,
                          pressure_deficit_Pa=-400.0, level_m=120.0, gust_front_convergence_s=1e-2)
    out = cl.classify(d, persistence_s=200.0)
    assert out["category"] == "SURFACE_CONNECTED_TORNADO_LIKE_VORTEX"
    assert out["rank"] == len(cl.CATEGORIES) - 1


def test_persistence_gates_tornado_tier():
    d = _base(); d.update(w_max=35.0, midlevel_mesocyclone=1.5e-2, updraft_helicity_2_5km=600.0,
                          near_surface_zeta_max=2e-2, v_theta_max_m_s=28.0, level_m=120.0,
                          gust_front_convergence_s=1e-2)
    # a transient (short-lived) intense vortex must not reach the tornado-like tiers
    assert cl.classify(d, persistence_s=10.0)["category"] == "LOW_LEVEL_MESOCYCLONE"


if __name__ == "__main__":
    for nm, fn in sorted(globals().items()):
        if nm.startswith("test_") and callable(fn):
            fn(); print("ok", nm)
    print("ALL CLASSIFICATION TESTS PASSED")
