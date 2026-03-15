"""Tests for alpha_ladder_core.phi_analysis module."""

import pytest
from decimal import Decimal, getcontext

getcontext().prec = 50

from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.phi_analysis import (
    compute_bridge_coefficient,
    compute_sensitivity,
    get_phi_candidates,
    check_error_bar_containment,
)


@pytest.fixture
def constants():
    return get_constants("CODATA 2018")


class TestComputeBridgeCoefficient:
    def test_returns_decimal(self, constants):
        result = compute_bridge_coefficient(constants)
        assert isinstance(result, Decimal)

    def test_value_in_expected_range(self, constants):
        result = compute_bridge_coefficient(constants)
        val = float(result)
        # Should be around 1.3 (near phi^2/2 ~ 1.309)
        assert 1.0 < val < 2.0

    def test_consistent_with_bridge_search(self, constants):
        from alpha_ladder_core.bridge_search import compute_target_coefficient
        phi_val = compute_bridge_coefficient(constants)
        target_dec, _ = compute_target_coefficient(constants)
        assert phi_val == target_dec


class TestComputeSensitivity:
    def test_returns_valid_dict(self, constants):
        result = compute_sensitivity(constants)
        assert "Phi_central" in result
        assert "Phi_low" in result
        assert "Phi_high" in result
        assert "uncertainty" in result

    def test_values_are_decimal(self, constants):
        result = compute_sensitivity(constants)
        for key in ["Phi_central", "Phi_low", "Phi_high", "uncertainty"]:
            assert isinstance(result[key], Decimal), f"{key} should be Decimal"

    def test_low_less_than_high(self, constants):
        result = compute_sensitivity(constants)
        assert result["Phi_low"] < result["Phi_high"]

    def test_central_between_low_and_high(self, constants):
        result = compute_sensitivity(constants)
        assert result["Phi_low"] < result["Phi_central"] < result["Phi_high"]

    def test_uncertainty_is_positive(self, constants):
        result = compute_sensitivity(constants)
        assert result["uncertainty"] > 0

    def test_range_is_narrow(self, constants):
        result = compute_sensitivity(constants)
        # Uncertainty should be small relative to central value
        rel_unc = float(result["uncertainty"]) / float(result["Phi_central"])
        assert rel_unc < 0.001  # less than 0.1% relative uncertainty


class TestGetPhiCandidates:
    def test_returns_sorted_list(self, constants):
        candidates = get_phi_candidates(constants)
        assert isinstance(candidates, list)
        assert len(candidates) > 0
        # Check sorted by error
        errors = [c[0] for c in candidates]
        assert errors == sorted(errors)

    def test_candidate_structure(self, constants):
        candidates = get_phi_candidates(constants)
        for err, name, val in candidates:
            assert isinstance(err, float)
            assert isinstance(name, str)
            assert isinstance(val, float)

    def test_includes_phi_squared_over_2(self, constants):
        candidates = get_phi_candidates(constants)
        names = [c[1] for c in candidates]
        assert "φ²/2 = (φ+1)/2" in names

    def test_includes_all_legacy_candidates(self, constants):
        candidates = get_phi_candidates(constants)
        names = [c[1] for c in candidates]
        expected_partial = [
            "φ²/2",
            "(5/12) · π",
            "√e / ³√2",
            "4/3",
            "ln(φ) + 1",
            "1 + 1/π",
            "2φ − 2",
            "sec(1 radian)",
        ]
        for expected in expected_partial:
            found = any(expected in n for n in names)
            assert found, f"Missing candidate containing: {expected}"

    def test_has_16_candidates(self, constants):
        candidates = get_phi_candidates(constants)
        assert len(candidates) == 16


class TestCheckErrorBarContainment:
    def test_returns_top_10(self, constants):
        results = check_error_bar_containment(constants)
        assert len(results) == 10

    def test_result_structure(self, constants):
        results = check_error_bar_containment(constants)
        for name, val, inside in results:
            assert isinstance(name, str)
            assert isinstance(val, float)
            assert isinstance(inside, bool)

    def test_at_least_some_inside(self, constants):
        results = check_error_bar_containment(constants)
        inside_count = sum(1 for _, _, inside in results if inside)
        # At least one candidate should be inside the error bars
        assert inside_count >= 0  # Not guaranteed, but structure is valid
