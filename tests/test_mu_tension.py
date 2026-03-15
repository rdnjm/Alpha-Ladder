"""Tests for alpha_ladder_core.mu_tension module."""

import pytest
from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.mu_tension import (
    compute_leading_order_identity,
    compute_formula_comparison,
    solve_mu_from_bridge,
    find_exact_c2,
    find_exact_c3,
    verify_unification,
    analyze_correction_hierarchy,
    predict_mu_from_alpha_phi,
    summarize_mu_tension,
)


@pytest.fixture
def constants():
    return get_constants("CODATA 2018")


# ---------------------------------------------------------------------------
# TestLeadingOrderIdentity
# ---------------------------------------------------------------------------
class TestLeadingOrderIdentity:
    def test_returns_dict(self, constants):
        result = compute_leading_order_identity(constants)
        assert isinstance(result, dict)

    def test_required_keys(self, constants):
        result = compute_leading_order_identity(constants)
        expected = {
            "mu_sq_alpha_cubed", "phi_sq_over_2", "ratio",
            "residual_ppm", "mu_predicted_leading", "assessment",
        }
        assert expected.issubset(set(result.keys()))

    def test_ratio_close_to_one(self, constants):
        result = compute_leading_order_identity(constants)
        assert abs(result["ratio"] - 1) < 0.001

    def test_residual_ppm_between_800_and_900(self, constants):
        result = compute_leading_order_identity(constants)
        assert 800 < abs(result["residual_ppm"]) < 900

    def test_mu_predicted_leading_in_range(self, constants):
        result = compute_leading_order_identity(constants)
        assert 1835 < result["mu_predicted_leading"] < 1837

    def test_assessment_nonempty(self, constants):
        result = compute_leading_order_identity(constants)
        assert len(result["assessment"]) > 0


# ---------------------------------------------------------------------------
# TestFormulaComparison
# ---------------------------------------------------------------------------
class TestFormulaComparison:
    def test_returns_dict(self, constants):
        result = compute_formula_comparison(constants)
        assert isinstance(result, dict)

    def test_required_keys(self, constants):
        result = compute_formula_comparison(constants)
        expected = {
            "C_mu", "C_bridge_lo", "C_bridge_c2", "C_bridge_nlo",
            "residuals", "brackets_exact", "assessment",
        }
        assert expected.issubset(set(result.keys()))

    def test_brackets_exact(self, constants):
        result = compute_formula_comparison(constants)
        assert result["brackets_exact"] is True

    def test_C_mu_positive(self, constants):
        result = compute_formula_comparison(constants)
        assert result["C_mu"] > 0

    def test_c2_residual_positive_and_small(self, constants):
        result = compute_formula_comparison(constants)
        r = result["residuals"]["c2_only"]
        assert 0 < r < 1

    def test_nlo_residual_negative_and_small(self, constants):
        result = compute_formula_comparison(constants)
        r = result["residuals"]["nlo"]
        assert -1 < r < 0

    def test_assessment_nonempty(self, constants):
        result = compute_formula_comparison(constants)
        assert len(result["assessment"]) > 0


# ---------------------------------------------------------------------------
# TestSolveMuFromBridge
# ---------------------------------------------------------------------------
class TestSolveMuFromBridge:
    def test_returns_dict(self, constants):
        result = solve_mu_from_bridge(constants=constants)
        assert isinstance(result, dict)

    def test_required_keys(self, constants):
        result = solve_mu_from_bridge(constants=constants)
        expected = {
            "c2", "c3", "F", "mu_predicted", "mu_measured",
            "difference", "residual_ppm", "tension_sigma", "assessment",
        }
        assert expected.issubset(set(result.keys()))

    def test_mu_predicted_in_range(self, constants):
        result = solve_mu_from_bridge(constants=constants)
        assert 1835 < result["mu_predicted"] < 1837

    def test_c2_3_c3_0_residual_less_than_1_ppm(self, constants):
        result = solve_mu_from_bridge(c2=3.0, c3=0.0, constants=constants)
        assert abs(result["residual_ppm"]) < 1

    def test_c2_3_c3_1_6_residual_less_than_1_ppm(self, constants):
        result = solve_mu_from_bridge(c2=3.0, c3=1.6, constants=constants)
        assert abs(result["residual_ppm"]) < 1

    def test_no_correction_residual_large(self, constants):
        result = solve_mu_from_bridge(c2=0.0, c3=0.0, constants=constants)
        assert abs(result["residual_ppm"]) > 50

    def test_tension_sigma_large_for_c2_3_c3_0(self, constants):
        result = solve_mu_from_bridge(c2=3.0, c3=0.0, constants=constants)
        assert result["tension_sigma"] > 100

    def test_assessment_nonempty(self, constants):
        result = solve_mu_from_bridge(constants=constants)
        assert len(result["assessment"]) > 0


# ---------------------------------------------------------------------------
# TestFindExactC2
# ---------------------------------------------------------------------------
class TestFindExactC2:
    def test_returns_dict(self, constants):
        result = find_exact_c2(constants)
        assert isinstance(result, dict)

    def test_required_keys(self, constants):
        result = find_exact_c2(constants)
        expected = {
            "c2_exact", "c2_integer_part", "c2_excess",
            "c2_excess_over_alpha", "assessment",
        }
        assert expected.issubset(set(result.keys()))

    def test_c2_exact_between_3_and_3_01(self, constants):
        result = find_exact_c2(constants)
        assert 3.0 < result["c2_exact"] < 3.01

    def test_c2_integer_part_is_3(self, constants):
        result = find_exact_c2(constants)
        assert result["c2_integer_part"] == 3

    def test_c2_excess_in_range(self, constants):
        result = find_exact_c2(constants)
        assert 0.005 < result["c2_excess"] < 0.007

    def test_c2_excess_over_alpha_in_range(self, constants):
        result = find_exact_c2(constants)
        assert 0.7 < result["c2_excess_over_alpha"] < 0.9


# ---------------------------------------------------------------------------
# TestFindExactC3
# ---------------------------------------------------------------------------
class TestFindExactC3:
    def test_returns_dict(self, constants):
        result = find_exact_c3(constants=constants)
        assert isinstance(result, dict)

    def test_required_keys(self, constants):
        result = find_exact_c3(constants=constants)
        expected = {
            "c3_exact", "c3_bridge", "difference",
            "c3_exact_over_phi", "c3_clean_candidates", "assessment",
        }
        assert expected.issubset(set(result.keys()))

    def test_c3_exact_between_0_7_and_0_9(self, constants):
        result = find_exact_c3(constants=constants)
        assert 0.7 < result["c3_exact"] < 0.9

    def test_c3_bridge_is_1_6(self, constants):
        result = find_exact_c3(constants=constants)
        assert result["c3_bridge"] == 1.6

    def test_difference_is_negative(self, constants):
        result = find_exact_c3(constants=constants)
        assert result["difference"] < 0

    def test_c3_clean_candidates_has_entries(self, constants):
        result = find_exact_c3(constants=constants)
        assert len(result["c3_clean_candidates"]) >= 3

    def test_assessment_nonempty(self, constants):
        result = find_exact_c3(constants=constants)
        assert len(result["assessment"]) > 0


# ---------------------------------------------------------------------------
# TestVerifyUnification
# ---------------------------------------------------------------------------
class TestVerifyUnification:
    def test_returns_dict(self, constants):
        result = verify_unification(constants)
        assert isinstance(result, dict)

    def test_required_keys(self, constants):
        result = verify_unification(constants)
        expected = {
            "G_bridge", "G_mu_structure", "G_measured",
            "bridge_residual_ppm", "mu_structure_residual_ppm",
            "difference_ppm", "unified", "assessment",
        }
        assert expected.issubset(set(result.keys()))

    def test_unified(self, constants):
        result = verify_unification(constants)
        assert result["unified"] is True

    def test_G_bridge_in_range(self, constants):
        result = verify_unification(constants)
        assert result["G_bridge"] > 6.6e-11

    def test_residuals_approximately_equal(self, constants):
        result = verify_unification(constants)
        diff = abs(result["bridge_residual_ppm"] - result["mu_structure_residual_ppm"])
        assert diff < 0.01

    def test_assessment_nonempty(self, constants):
        result = verify_unification(constants)
        assert len(result["assessment"]) > 0


# ---------------------------------------------------------------------------
# TestCorrectionHierarchy
# ---------------------------------------------------------------------------
class TestCorrectionHierarchy:
    def test_returns_dict(self, constants):
        result = analyze_correction_hierarchy(constants)
        assert isinstance(result, dict)

    def test_required_keys(self, constants):
        result = analyze_correction_hierarchy(constants)
        expected = {"terms", "converges", "convergence_ratio", "assessment"}
        assert expected.issubset(set(result.keys()))

    def test_terms_has_three_entries(self, constants):
        result = analyze_correction_hierarchy(constants)
        assert len(result["terms"]) >= 3

    def test_converges(self, constants):
        result = analyze_correction_hierarchy(constants)
        assert result["converges"] is True

    def test_convergence_ratio_small(self, constants):
        result = analyze_correction_hierarchy(constants)
        assert result["convergence_ratio"] < 0.01

    def test_first_term_largest_contribution(self, constants):
        result = analyze_correction_hierarchy(constants)
        terms = result["terms"]
        first_abs = abs(terms[0]["contribution_ppm"])
        for t in terms[1:]:
            assert first_abs > abs(t["contribution_ppm"])


# ---------------------------------------------------------------------------
# TestPredictMuFromAlphaPhi
# ---------------------------------------------------------------------------
class TestPredictMuFromAlphaPhi:
    def test_returns_dict(self, constants):
        result = predict_mu_from_alpha_phi(constants)
        assert isinstance(result, dict)

    def test_required_keys(self, constants):
        result = predict_mu_from_alpha_phi(constants)
        expected = {
            "mu_predicted", "mu_measured", "residual_ppm",
            "is_prediction", "n_fitted_params", "assessment",
        }
        assert expected.issubset(set(result.keys()))

    def test_mu_predicted_in_range(self, constants):
        result = predict_mu_from_alpha_phi(constants)
        assert 1835 < result["mu_predicted"] < 1837

    def test_residual_ppm_less_than_1(self, constants):
        result = predict_mu_from_alpha_phi(constants)
        assert abs(result["residual_ppm"]) < 1

    def test_is_prediction(self, constants):
        result = predict_mu_from_alpha_phi(constants)
        assert result["is_prediction"] is True

    def test_n_fitted_params_zero(self, constants):
        result = predict_mu_from_alpha_phi(constants)
        assert result["n_fitted_params"] == 0

    def test_assessment_nonempty(self, constants):
        result = predict_mu_from_alpha_phi(constants)
        assert len(result["assessment"]) > 0


# ---------------------------------------------------------------------------
# TestSummarizeMuTension
# ---------------------------------------------------------------------------
class TestSummarizeMuTension:
    def test_returns_dict(self, constants):
        result = summarize_mu_tension(constants)
        assert isinstance(result, dict)

    def test_all_sub_keys_present(self, constants):
        result = summarize_mu_tension(constants)
        expected = {
            "leading_order", "formula_comparison", "mu_from_bridge",
            "exact_c2", "exact_c3", "unification", "correction_hierarchy",
            "mu_prediction", "key_finding", "honest_assessment",
        }
        assert expected.issubset(set(result.keys()))

    def test_key_finding_nonempty(self, constants):
        result = summarize_mu_tension(constants)
        assert len(result["key_finding"]) > 0

    def test_honest_assessment_nonempty(self, constants):
        result = summarize_mu_tension(constants)
        assert len(result["honest_assessment"]) > 0

    def test_sub_results_are_dicts(self, constants):
        result = summarize_mu_tension(constants)
        for key in ["leading_order", "formula_comparison", "mu_from_bridge",
                     "exact_c2", "exact_c3", "unification",
                     "correction_hierarchy", "mu_prediction"]:
            assert isinstance(result[key], dict), f"{key} should be a dict"

    def test_mu_prediction_in_result(self, constants):
        result = summarize_mu_tension(constants)
        assert "mu_prediction" in result
        assert "mu_predicted" in result["mu_prediction"]
