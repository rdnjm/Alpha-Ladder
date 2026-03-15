"""Tests for alpha_ladder_core.predict_g module."""

import pytest
from decimal import Decimal, getcontext

getcontext().prec = 50

from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.predict_g import (
    get_bridge_candidates,
    predict_G,
    predict_G_hierarchy,
    get_G_measurements,
    compare_prediction,
    summarize_predictions,
)


@pytest.fixture
def constants():
    return get_constants("CODATA 2018")


class TestGetBridgeCandidates:
    def test_returns_three_candidates(self, constants):
        candidates = get_bridge_candidates(constants)
        assert len(candidates) == 3

    def test_candidate_names(self, constants):
        candidates = get_bridge_candidates(constants)
        assert "φ²/2" in candidates
        assert "(5/12) · π" in candidates
        assert "√e / ³√2" in candidates

    def test_candidate_values_are_decimal(self, constants):
        candidates = get_bridge_candidates(constants)
        for name, val in candidates.items():
            assert isinstance(val, Decimal), f"{name} should be Decimal"

    def test_phi_squared_over_2_value(self, constants):
        candidates = get_bridge_candidates(constants)
        val = float(candidates["φ²/2"])
        # phi^2/2 = (phi+1)/2 ~ 1.309
        assert abs(val - 1.309) < 0.01


class TestPredictG:
    def test_phi_squared_over_2_gives_correct_G(self, constants):
        candidates = get_bridge_candidates(constants)
        G_pred = predict_G(candidates["φ²/2"], constants)
        G_float = float(G_pred)
        # Should be approximately 6.67323e-11 (close to measured G)
        assert abs(G_float - 6.67e-11) / 6.67e-11 < 0.01  # within 1%

    def test_returns_decimal(self, constants):
        candidates = get_bridge_candidates(constants)
        G_pred = predict_G(candidates["φ²/2"], constants)
        assert isinstance(G_pred, Decimal)

    def test_different_coefficients_give_different_G(self, constants):
        candidates = get_bridge_candidates(constants)
        values = [predict_G(c, constants) for c in candidates.values()]
        assert len(set(values)) == 3


class TestGetGMeasurements:
    def test_returns_seven_measurements(self):
        measurements = get_G_measurements()
        assert len(measurements) == 7

    def test_all_measurements_are_decimal_tuples(self):
        measurements = get_G_measurements()
        for name, (val, unc) in measurements.items():
            assert isinstance(val, Decimal), f"{name} value should be Decimal"
            assert isinstance(unc, Decimal), f"{name} uncertainty should be Decimal"

    def test_codata_2018_present(self):
        measurements = get_G_measurements()
        assert "CODATA 2018 recommended" in measurements

    def test_all_expected_experiments_present(self):
        measurements = get_G_measurements()
        expected = [
            "CODATA 2018 recommended",
            "Quinn et al. 2013 (BIPM)",
            "Rosi et al. 2014 (atom interf.)",
            "Newman et al. 2014",
            "Li et al. 2018 (HUST-A)",
            "Li et al. 2018 (HUST-B)",
            "CODATA 2014 recommended",
        ]
        for name in expected:
            assert name in measurements, f"Missing measurement: {name}"

    def test_G_values_in_expected_range(self):
        measurements = get_G_measurements()
        for name, (val, unc) in measurements.items():
            G_float = float(val)
            assert 6.67e-11 < G_float < 6.68e-11, f"{name} G value out of range"

    def test_hust_uncertainties_match_published(self):
        measurements = get_G_measurements()
        expected_unc = Decimal("0.000078e-11")
        _, unc_a = measurements["Li et al. 2018 (HUST-A)"]
        _, unc_b = measurements["Li et al. 2018 (HUST-B)"]
        assert unc_a == expected_unc, f"HUST-A uncertainty {unc_a} != {expected_unc}"
        assert unc_b == expected_unc, f"HUST-B uncertainty {unc_b} != {expected_unc}"


class TestComparePrediction:
    def test_returns_correct_structure(self, constants):
        candidates = get_bridge_candidates(constants)
        G_pred = predict_G(candidates["φ²/2"], constants)
        measurements = get_G_measurements()
        results = compare_prediction(G_pred, measurements)

        assert len(results) == 7
        for r in results:
            assert "experiment" in r
            assert "G_exp" in r
            assert "G_unc" in r
            assert "sigma" in r
            assert "direction" in r

    def test_sigma_values_are_positive(self, constants):
        candidates = get_bridge_candidates(constants)
        G_pred = predict_G(candidates["φ²/2"], constants)
        measurements = get_G_measurements()
        results = compare_prediction(G_pred, measurements)

        for r in results:
            assert r["sigma"] >= 0

    def test_direction_is_plus_or_minus(self, constants):
        candidates = get_bridge_candidates(constants)
        G_pred = predict_G(candidates["φ²/2"], constants)
        measurements = get_G_measurements()
        results = compare_prediction(G_pred, measurements)

        for r in results:
            assert r["direction"] in ("+", "-")


class TestSummarizePredictions:
    def test_returns_all_three_bridges(self, constants):
        summary = summarize_predictions(constants)
        assert len(summary) == 3
        assert "φ²/2" in summary
        assert "(5/12) · π" in summary
        assert "√e / ³√2" in summary

    def test_summary_structure(self, constants):
        summary = summarize_predictions(constants)
        for name, data in summary.items():
            assert "G_pred" in data
            assert "avg_sigma" in data
            assert "comparisons" in data
            assert isinstance(data["comparisons"], list)
            assert len(data["comparisons"]) == 7


class TestPredictGHierarchy:
    def test_returns_dict_with_all_expected_keys(self, constants):
        result = predict_G_hierarchy(constants)
        expected_keys = {"alpha_g", "G_predicted", "mu", "exponent", "residual_ppm"}
        assert set(result.keys()) == expected_keys

    def test_exponent_is_24(self, constants):
        result = predict_G_hierarchy(constants)
        assert result["exponent"] == 24

    def test_mu_approximately_1836(self, constants):
        result = predict_G_hierarchy(constants)
        mu_float = float(result["mu"])
        assert abs(mu_float - 1836.15) < 0.01, f"mu = {mu_float}, expected ~1836.15"

    def test_G_predicted_in_valid_range(self, constants):
        result = predict_G_hierarchy(constants)
        G_float = float(result["G_predicted"])
        assert 6.67e-11 <= G_float <= 6.68e-11, f"G_predicted = {G_float} out of range"

    def test_G_predicted_approximately_hierarchy_value(self, constants):
        result = predict_G_hierarchy(constants)
        G_float = float(result["G_predicted"])
        expected = 6.6789e-11
        relative_error = abs(G_float - expected) / expected
        assert relative_error < 0.001, (
            f"G_predicted = {G_float}, expected ~{expected} (within 0.1%)"
        )

    def test_residual_ppm_approximately_688(self, constants):
        result = predict_G_hierarchy(constants)
        assert abs(result["residual_ppm"] - 688) < 50, (
            f"residual_ppm = {result['residual_ppm']}, expected ~688 (within 50)"
        )

    def test_alpha_g_is_decimal(self, constants):
        result = predict_G_hierarchy(constants)
        assert isinstance(result["alpha_g"], Decimal)

    def test_G_predicted_is_decimal(self, constants):
        result = predict_G_hierarchy(constants)
        assert isinstance(result["G_predicted"], Decimal)
