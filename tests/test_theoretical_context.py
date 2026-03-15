"""
Tests for alpha_ladder_core/theoretical_context.py
"""

import math
import pytest
from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.theoretical_context import (
    compute_screening_discrepancy,
    position_in_literature,
    analyze_anomaly_status,
    summarize_theoretical_status,
)


@pytest.fixture
def constants():
    return get_constants()


# -----------------------------------------------------------------------
# TestScreeningDiscrepancy
# -----------------------------------------------------------------------
class TestScreeningDiscrepancy:
    def test_ratio_order_of_magnitude(self, constants):
        """Ratio should be ~10^3.5 (around 3000-4000)."""
        result = compute_screening_discrepancy(constants)
        assert 1000 < abs(result["ratio"]) < 10000

    def test_alpha_tree_order(self, constants):
        """alpha_tree ~ 0.6."""
        result = compute_screening_discrepancy(constants)
        assert 0.1 < result["alpha_tree"] < 1.0

    def test_alpha_empirical_order(self, constants):
        """alpha_empirical ~ 1.6e-4."""
        result = compute_screening_discrepancy(constants)
        assert 1e-5 < abs(result["alpha_empirical"]) < 1e-3

    def test_log10_ratio(self, constants):
        result = compute_screening_discrepancy(constants)
        expected = math.log10(abs(result["ratio"]))
        assert pytest.approx(result["log10_ratio"], rel=1e-5) == expected

    def test_resolutions_nonempty(self, constants):
        result = compute_screening_discrepancy(constants)
        assert len(result["resolutions"]) >= 3

    def test_resolutions_have_required_keys(self, constants):
        result = compute_screening_discrepancy(constants)
        for res in result["resolutions"]:
            assert "name" in res
            assert "mechanism" in res
            assert "plausibility" in res
            assert "description" in res

    def test_resolutions_no_f2_high_plausibility(self, constants):
        """No resolution involving F^2 or trace anomaly should have high plausibility."""
        result = compute_screening_discrepancy(constants)
        for res in result["resolutions"]:
            if res["plausibility"] == "high":
                name_lower = res["name"].lower()
                mech_lower = res["mechanism"].lower()
                assert "f^2" not in name_lower and "f^2" not in mech_lower, (
                    f"Resolution '{res['name']}' with high plausibility should not claim F^2 suppression"
                )
                assert "trace anomaly" not in name_lower and "trace anomaly" not in mech_lower, (
                    f"Resolution '{res['name']}' with high plausibility should not claim trace anomaly suppression"
                )

    def test_heavy_dilaton_resolution_present(self, constants):
        """A 'Heavy dilaton' resolution must exist with high plausibility."""
        result = compute_screening_discrepancy(constants)
        heavy = [r for r in result["resolutions"] if "heavy dilaton" in r["name"].lower()]
        assert len(heavy) == 1, "Expected exactly one 'Heavy dilaton' resolution"
        assert heavy[0]["plausibility"] == "high"


# -----------------------------------------------------------------------
# TestLiteraturePosition
# -----------------------------------------------------------------------
class TestLiteraturePosition:
    def test_four_frameworks(self):
        result = position_in_literature()
        assert len(result["comparisons"]) == 4

    def test_each_has_similarities_and_differences(self):
        result = position_in_literature()
        for comp in result["comparisons"]:
            assert len(comp["similarities"]) >= 1
            assert len(comp["differences"]) >= 1

    def test_salam_sezgin_present(self):
        result = position_in_literature()
        names = [c["framework"] for c in result["comparisons"]]
        assert any("Salam" in n for n in names)

    def test_each_has_reference(self):
        result = position_in_literature()
        for comp in result["comparisons"]:
            assert "key_reference" in comp
            assert len(comp["key_reference"]) > 0

    def test_kklt_present(self):
        result = position_in_literature()
        names = [c["framework"] for c in result["comparisons"]]
        assert any("KKLT" in n for n in names)


# -----------------------------------------------------------------------
# TestAnomalyStatus
# -----------------------------------------------------------------------
class TestAnomalyStatus:
    def test_returns_dict(self):
        result = analyze_anomaly_status()
        assert isinstance(result, dict)

    def test_has_anomaly_types(self):
        result = analyze_anomaly_status()
        assert "anomaly_types" in result
        assert len(result["anomaly_types"]) >= 2

    def test_has_interpretation(self):
        result = analyze_anomaly_status()
        assert "interpretation" in result
        assert isinstance(result["interpretation"], str)
        assert len(result["interpretation"]) > 0

    def test_green_schwarz_field(self):
        result = analyze_anomaly_status()
        assert "green_schwarz_required" in result
        assert isinstance(result["green_schwarz_required"], bool)


# -----------------------------------------------------------------------
# TestSummarize
# -----------------------------------------------------------------------
class TestSummarize:
    def test_strengths_nonempty(self, constants):
        result = summarize_theoretical_status(constants)
        assert len(result["strengths"]) >= 3

    def test_open_problems_nonempty(self, constants):
        result = summarize_theoretical_status(constants)
        assert len(result["open_problems"]) >= 3

    def test_contains_discrepancy(self, constants):
        result = summarize_theoretical_status(constants)
        assert "discrepancy" in result

    def test_contains_literature(self, constants):
        result = summarize_theoretical_status(constants)
        assert "literature" in result

    def test_contains_anomalies(self, constants):
        result = summarize_theoretical_status(constants)
        assert "anomalies" in result
