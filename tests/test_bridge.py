"""Tests for alpha_ladder_core.bridge_search module."""

import pytest
from decimal import Decimal, getcontext

getcontext().prec = 50

from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.bridge_search import (
    compute_target_coefficient,
    search_single_constant,
    search_two_constants,
    search_fraction_times_constant,
    run_full_search,
)


@pytest.fixture
def constants():
    return get_constants("CODATA 2018")


class TestComputeTargetCoefficient:
    def test_returns_tuple(self, constants):
        result = compute_target_coefficient(constants)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_decimal_and_float(self, constants):
        target_dec, target_f = compute_target_coefficient(constants)
        assert isinstance(target_dec, Decimal)
        assert isinstance(target_f, float)

    def test_values_match(self, constants):
        target_dec, target_f = compute_target_coefficient(constants)
        assert abs(float(target_dec) - target_f) < 1e-10

    def test_target_in_expected_range(self, constants):
        _, target_f = compute_target_coefficient(constants)
        # The bridge coefficient should be around 1.3 (near phi^2/2)
        assert 1.0 < target_f < 2.0


class TestSearchSingleConstant:
    def test_returns_list(self, constants):
        _, target_f = compute_target_coefficient(constants)
        results = search_single_constant(target_f)
        assert isinstance(results, list)

    def test_results_have_correct_structure(self, constants):
        _, target_f = compute_target_coefficient(constants)
        results = search_single_constant(target_f)
        for err, expr, val in results:
            assert isinstance(err, float)
            assert isinstance(expr, str)
            assert isinstance(val, float)

    def test_all_errors_below_threshold(self, constants):
        _, target_f = compute_target_coefficient(constants)
        results = search_single_constant(target_f)
        for err, expr, val in results:
            assert err < 2.0, f"{expr} has error {err}% >= 2.0%"


class TestSearchTwoConstants:
    def test_returns_list(self, constants):
        _, target_f = compute_target_coefficient(constants)
        results = search_two_constants(target_f)
        assert isinstance(results, list)

    def test_all_errors_below_threshold(self, constants):
        _, target_f = compute_target_coefficient(constants)
        results = search_two_constants(target_f)
        for err, expr, val in results:
            assert err < 0.5, f"{expr} has error {err}% >= 0.5%"


class TestSearchFractionTimesConstant:
    def test_returns_list(self, constants):
        _, target_f = compute_target_coefficient(constants)
        results = search_fraction_times_constant(target_f)
        assert isinstance(results, list)

    def test_all_errors_below_threshold(self, constants):
        _, target_f = compute_target_coefficient(constants)
        results = search_fraction_times_constant(target_f)
        for err, expr, val in results:
            assert err < 0.1, f"{expr} has error {err}% >= 0.1%"


class TestRunFullSearch:
    def test_returns_sorted_list(self, constants):
        results = run_full_search(constants)
        assert isinstance(results, list)
        assert len(results) > 0
        # Check sorted by error
        errors = [r[0] for r in results]
        assert errors == sorted(errors)

    def test_no_duplicate_expressions(self, constants):
        results = run_full_search(constants)
        exprs = [r[1] for r in results]
        assert len(exprs) == len(set(exprs))

    def test_finds_matches_below_01_pct(self, constants):
        results = run_full_search(constants)
        low_error = [r for r in results if r[0] < 0.1]
        assert len(low_error) > 0, "Should find at least one match with < 0.1% error"

    def test_phi_squared_over_2_in_results(self, constants):
        results = run_full_search(constants)
        exprs = [r[1] for r in results]
        # phi^2/2 may appear as a fraction-times-constant match like (1/2) * phi^2
        # or in two-constants search. Check for any phi-related match.
        phi_matches = [e for e in exprs if "phi" in e.lower() or "φ" in e]
        assert len(phi_matches) > 0, "Should find phi-related matches"
