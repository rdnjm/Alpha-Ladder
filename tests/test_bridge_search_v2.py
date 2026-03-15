"""Tests for alpha_ladder_core.bridge_search_v2 module."""

import pytest
import math
from decimal import Decimal, getcontext

getcontext().prec = 50

from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.bridge_search_v2 import (
    compute_exact_bridge,
    search_simple_expressions,
    search_alpha_mu_bridges,
    search_hybrid_bridges,
    summarize_bridge_search,
)


PHI = (1 + math.sqrt(5)) / 2
PHI_SQ_OVER_2 = PHI ** 2 / 2  # ~1.30902


@pytest.fixture
def constants():
    return get_constants("CODATA 2018")


# ---------------------------------------------------------------------------
# TestComputeExactBridge
# ---------------------------------------------------------------------------
class TestComputeExactBridge:
    def test_returns_dict(self, constants):
        result = compute_exact_bridge(constants)
        assert isinstance(result, dict)

    def test_c_exact_is_decimal(self, constants):
        result = compute_exact_bridge(constants)
        assert isinstance(result["C_exact"], Decimal)

    def test_c_float_is_float(self, constants):
        result = compute_exact_bridge(constants)
        assert isinstance(result["C_float"], float)

    def test_c_float_in_expected_range(self, constants):
        result = compute_exact_bridge(constants)
        assert 1.309 <= result["C_float"] <= 1.310

    def test_c_float_above_phi_sq_over_2(self, constants):
        result = compute_exact_bridge(constants)
        assert result["C_float"] > PHI_SQ_OVER_2

    def test_c_exact_high_precision(self, constants):
        result = compute_exact_bridge(constants)
        # At least 10 significant digits: string repr should have 10+ digits
        digits = str(result["C_exact"]).replace(".", "").replace("-", "").lstrip("0")
        assert len(digits) >= 10


# ---------------------------------------------------------------------------
# TestSearchSimpleExpressions
# ---------------------------------------------------------------------------
class TestSearchSimpleExpressions:
    def test_returns_list(self, constants):
        results = search_simple_expressions(constants)
        assert isinstance(results, list)

    def test_entries_have_required_keys(self, constants):
        results = search_simple_expressions(constants)
        required = {"name", "value", "residual_ppm", "has_fitted_params", "category"}
        for entry in results:
            assert required.issubset(entry.keys()), f"Missing keys in {entry}"

    def test_sorted_by_abs_residual_ppm(self, constants):
        results = search_simple_expressions(constants)
        residuals = [abs(r["residual_ppm"]) for r in results]
        assert residuals == sorted(residuals)

    def test_all_within_max_ppm(self, constants):
        """Default max_ppm is 1000; all entries should be within that."""
        results = search_simple_expressions(constants)
        for r in results:
            assert abs(r["residual_ppm"]) < 1000, (
                f"{r['name']} has |residual_ppm| = {abs(r['residual_ppm'])} >= 1000"
            )

    def test_custom_max_ppm_filters(self, constants):
        all_results = search_simple_expressions(constants)
        filtered = search_simple_expressions(constants, max_ppm=200)
        assert len(filtered) <= len(all_results)
        for r in filtered:
            assert abs(r["residual_ppm"]) < 200

    def test_phi_sq_over_2_present(self, constants):
        results = search_simple_expressions(constants)
        names = [r["name"] for r in results]
        phi_matches = [n for n in names if "phi" in n.lower()]
        assert len(phi_matches) > 0, "phi^2/2 should appear in simple expression results"

    def test_5_12_pi_present(self, constants):
        results = search_simple_expressions(constants)
        names = [r["name"] for r in results]
        pi_matches = [n for n in names if "pi" in n.lower()]
        assert len(pi_matches) > 0, "(5/12)*pi should appear in simple expression results"

    def test_valid_categories(self, constants):
        results = search_simple_expressions(constants)
        valid = {"pure_math", "physical", "mixed"}
        for r in results:
            assert r["category"] in valid, f"Invalid category '{r['category']}' in {r['name']}"


# ---------------------------------------------------------------------------
# TestSearchAlphaMuBridges
# ---------------------------------------------------------------------------
class TestSearchAlphaMuBridges:
    def test_returns_list(self, constants):
        results = search_alpha_mu_bridges(constants)
        assert isinstance(results, list)

    def test_entries_have_required_keys(self, constants):
        results = search_alpha_mu_bridges(constants)
        required = {
            "a", "b", "p", "q",
            "total_exponent", "mu_exponent",
            "residual_ppm", "formula_str",
        }
        for entry in results:
            assert required.issubset(entry.keys()), f"Missing keys in {entry}"

    def test_sorted_by_abs_residual_ppm(self, constants):
        results = search_alpha_mu_bridges(constants)
        residuals = [abs(r["residual_ppm"]) for r in results]
        assert residuals == sorted(residuals)

    def test_all_within_max_ppm(self, constants):
        results = search_alpha_mu_bridges(constants)
        for r in results:
            assert abs(r["residual_ppm"]) < 1000, (
                f"Entry has |residual_ppm| = {abs(r['residual_ppm'])} >= 1000"
            )

    def test_alpha3_mu2_present(self, constants):
        """The case alpha^3 * mu^2 (a=3, b=2) should appear with ~688 ppm."""
        results = search_alpha_mu_bridges(constants)
        match = [r for r in results if r["a"] == 3 and r["b"] == 2]
        assert len(match) > 0, "alpha^3 * mu^2 bridge should be in results"
        assert abs(match[0]["residual_ppm"]) < 1000

    def test_total_exponent_formula(self, constants):
        """total_exponent should equal 21 + a for each entry."""
        results = search_alpha_mu_bridges(constants)
        for r in results:
            assert r["total_exponent"] == 21 + r["a"], (
                f"total_exponent={r['total_exponent']} != 21 + a={21 + r['a']}"
            )

    def test_mu_exponent_equals_b(self, constants):
        """mu_exponent should equal b for each entry."""
        results = search_alpha_mu_bridges(constants)
        for r in results:
            assert r["mu_exponent"] == r["b"], (
                f"mu_exponent={r['mu_exponent']} != b={r['b']}"
            )


# ---------------------------------------------------------------------------
# TestSearchHybridBridges
# ---------------------------------------------------------------------------
class TestSearchHybridBridges:
    def test_returns_list(self, constants):
        results = search_hybrid_bridges(constants)
        assert isinstance(results, list)

    def test_entries_have_required_keys_when_nonempty(self, constants):
        results = search_hybrid_bridges(constants, max_ppm=1000)
        for entry in results:
            assert "residual_ppm" in entry
            assert "formula_str" in entry

    def test_sorted_by_abs_residual_ppm(self, constants):
        results = search_hybrid_bridges(constants, max_ppm=1000)
        residuals = [abs(r["residual_ppm"]) for r in results]
        assert residuals == sorted(residuals)

    def test_all_within_max_ppm(self, constants):
        results = search_hybrid_bridges(constants, max_ppm=1000)
        for r in results:
            assert abs(r["residual_ppm"]) < 1000, (
                f"Entry has |residual_ppm| = {abs(r['residual_ppm'])} >= 1000"
            )

    def test_custom_max_ppm_filters(self, constants):
        wide = search_hybrid_bridges(constants, max_ppm=1000)
        narrow = search_hybrid_bridges(constants, max_ppm=100)
        assert len(narrow) <= len(wide)

    def test_default_max_ppm_is_strict(self, constants):
        """Default max_ppm=50 is strict; results may be empty."""
        results = search_hybrid_bridges(constants)
        for r in results:
            assert abs(r["residual_ppm"]) < 50


# ---------------------------------------------------------------------------
# TestSummarizeBridgeSearch
# ---------------------------------------------------------------------------
class TestSummarizeBridgeSearch:
    @pytest.fixture(autouse=True)
    def _summary(self, constants):
        self.summary = summarize_bridge_search(constants)

    def test_returns_dict(self):
        assert isinstance(self.summary, dict)

    def test_has_required_keys(self):
        required = {
            "C_exact", "C_float",
            "best_pure_math", "best_physical", "best_hybrid",
            "top_20", "key_insight",
        }
        assert required.issubset(self.summary.keys())

    def test_top_20_is_list(self):
        top = self.summary["top_20"]
        assert isinstance(top, list)
        assert len(top) <= 20

    def test_top_20_sorted(self):
        top = self.summary["top_20"]
        residuals = [abs(r["residual_ppm"]) for r in top]
        assert residuals == sorted(residuals)

    def test_key_insight_non_empty(self):
        assert isinstance(self.summary["key_insight"], str)
        assert len(self.summary["key_insight"]) > 0

    def test_best_pure_math_has_fields(self):
        bpm = self.summary["best_pure_math"]
        assert bpm is not None, "Expected at least one pure_math result"
        assert "name" in bpm
        assert "residual_ppm" in bpm

    def test_best_physical_has_fields(self):
        bp = self.summary["best_physical"]
        assert bp is not None, "Expected at least one physical result"
        assert "name" in bp
        assert "residual_ppm" in bp

    def test_best_hybrid_is_none_or_dict(self):
        """best_hybrid may be None if no hybrid results are within threshold."""
        bh = self.summary["best_hybrid"]
        if bh is not None:
            assert "residual_ppm" in bh
