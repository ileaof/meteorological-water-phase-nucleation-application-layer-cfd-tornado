"""The experiment matrix runs and its shear control behaves physically (shear -> rotation)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "examples"))

import experiment_matrix as em


def test_run_one_returns_expected_keys():
    r = em._run_one("baseline", {}, nx=20, nz=28, steps=30, device="cpu")
    for k in ("name", "dx_m", "w_max", "midlevel_meso", "near_surface_zeta", "UH_2_5km", "category"):
        assert k in r


def test_shear_control_produces_more_rotation():
    base = em._run_one("baseline", {}, nx=20, nz=28, steps=60, device="cpu")
    nosh = em._run_one("no_shear", {"U_max": 0.0}, nx=20, nz=28, steps=60, device="cpu")
    # removing the vertical shear must strongly reduce the mid-level rotation it can tilt/stretch
    assert base["midlevel_meso"] > 3.0 * nosh["midlevel_meso"]


if __name__ == "__main__":
    test_run_one_returns_expected_keys(); print("ok keys")
    test_shear_control_produces_more_rotation(); print("ok shear control")
    print("ALL EXPERIMENT-MATRIX TESTS PASSED")
