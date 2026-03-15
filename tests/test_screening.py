"""Tests for the Yukawa screening model."""

import math
import pytest
from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.screening import (
    compute_screening_parameters,
    compute_G_eff,
    compute_screening_profile,
    classify_measurements,
)


@pytest.fixture
def constants():
    return get_constants("CODATA 2018")


@pytest.fixture
def screening_params(constants):
    return compute_screening_parameters(constants)


class TestComputeScreeningParameters:
    """Tests for compute_screening_parameters."""

    def test_g_vacuum_matches_predict_g(self, constants, screening_params):
        from alpha_ladder_core.predict_g import predict_G, get_bridge_candidates
        bridges = get_bridge_candidates(constants)
        expected = float(predict_G(bridges["\u03c6\u00b2/2"], constants))
        assert screening_params["G_vacuum"] == pytest.approx(expected, rel=1e-10)

    def test_g_vacuum_order_of_magnitude(self, screening_params):
        assert 6e-11 < screening_params["G_vacuum"] < 7e-11

    def test_g_lab_is_codata(self, screening_params):
        assert screening_params["G_lab"] == pytest.approx(6.67430e-11, rel=1e-4)

    def test_alpha_screening_range(self, screening_params):
        assert 1e-4 < screening_params["alpha_screening"] < 3e-4

    def test_ppm_excess_around_160(self, screening_params):
        # alpha_screening ~ 1.6e-4, so ppm = 1.6e-4 * 1e6 = ~160
        assert 100 < screening_params["ppm_excess"] < 250

    def test_lambda_dilaton_positive(self, screening_params):
        assert screening_params["lambda_dilaton_m"] > 0

    def test_lambda_dilaton_near_au(self, screening_params):
        AU = 1.496e11
        # Should be on the order of AU (within a few orders of magnitude)
        ratio = screening_params["lambda_dilaton_m"] / AU
        assert 0.01 < ratio < 100

    def test_dilaton_mass_positive(self, screening_params):
        assert screening_params["dilaton_mass_kg"] > 0
        assert screening_params["dilaton_mass_eV"] > 0


class TestComputeGEff:
    """Tests for compute_G_eff."""

    def test_near_zero_gives_max(self, screening_params):
        G_near = compute_G_eff(0.01, screening_params)
        G_vac = screening_params["G_vacuum"]
        alpha_s = screening_params["alpha_screening"]
        expected = G_vac * (1 + alpha_s)
        assert G_near == pytest.approx(expected, rel=1e-6)

    def test_large_distance_converges_to_vacuum(self, screening_params):
        G_far = compute_G_eff(1e15, screening_params)
        G_vac = screening_params["G_vacuum"]
        assert G_far == pytest.approx(G_vac, rel=1e-10)

    def test_monotonically_decreasing(self, screening_params):
        distances = [0.01, 1.0, 1e6, 1e11, 1e15]
        G_values = [compute_G_eff(r, screening_params) for r in distances]
        for i in range(len(G_values) - 1):
            assert G_values[i] >= G_values[i + 1]

    def test_always_above_vacuum(self, screening_params):
        G_vac = screening_params["G_vacuum"]
        for r in [0.01, 1.0, 1e6, 1e11, 1e15]:
            assert compute_G_eff(r, screening_params) >= G_vac


class TestComputeScreeningProfile:
    """Tests for compute_screening_profile."""

    def test_correct_n_points(self, screening_params):
        profile = compute_screening_profile(screening_params, n_points=100)
        assert len(profile["r_meters"]) == 100
        assert len(profile["G_eff"]) == 100

    def test_default_n_points(self, screening_params):
        profile = compute_screening_profile(screening_params)
        assert len(profile["r_meters"]) == 500

    def test_r_meters_sorted(self, screening_params):
        profile = compute_screening_profile(screening_params, n_points=50)
        for i in range(len(profile["r_meters"]) - 1):
            assert profile["r_meters"][i] < profile["r_meters"][i + 1]

    def test_all_g_eff_positive(self, screening_params):
        profile = compute_screening_profile(screening_params, n_points=50)
        assert all(g > 0 for g in profile["G_eff"])

    def test_landmarks_present(self, screening_params):
        profile = compute_screening_profile(screening_params)
        assert "Lab bench (1 m)" in profile["landmarks"]
        assert "1 AU (1.496e11 m)" in profile["landmarks"]
        assert len(profile["landmarks"]) == 4


class TestClassifyMeasurements:
    """Tests for classify_measurements."""

    def test_total_count(self, screening_params, constants):
        result = classify_measurements(screening_params, constants)
        total = len(result["high_cluster"]) + len(result["low_cluster"])
        assert total == 7

    def test_cluster_split(self, screening_params, constants):
        result = classify_measurements(screening_params, constants)
        assert len(result["high_cluster"]) == 6
        assert len(result["low_cluster"]) == 1

    def test_rosi_in_low_cluster(self, screening_params, constants):
        result = classify_measurements(screening_params, constants)
        low_names = [m["experiment"] for m in result["low_cluster"]]
        assert any("Rosi" in n for n in low_names)

    def test_high_mean_greater_than_low(self, screening_params, constants):
        result = classify_measurements(screening_params, constants)
        assert result["high_mean"] > result["low_mean"]

    def test_all_have_cluster_tag(self, screening_params, constants):
        result = classify_measurements(screening_params, constants)
        for m in result["measurements"]:
            assert m["cluster"] in ("high", "low")
