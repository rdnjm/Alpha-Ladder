"""Tests for alpha_ladder_core.dimension_uniqueness module.

Validates the three-gap resolution: dimension uniqueness, c3 = phi/2,
and mu gap closure.
"""

import math
import unittest

from alpha_ladder_core.dimension_uniqueness import (
    scan_exponent_constraint,
    scan_volume_cancellation,
    scan_vacuum_polynomial,
    prove_uniqueness,
    derive_c3_phi_half,
    predict_mu_unified,
    predict_G_unified,
    compute_complete_formula,
    summarize_dimension_uniqueness,
)


class TestScanExponentConstraint(unittest.TestCase):
    """Tests for scan_exponent_constraint."""

    def test_returns_dict(self):
        result = scan_exponent_constraint()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = scan_exponent_constraint()
        for key in ["scan_results", "sub_1000_pairs", "n_sub_1000", "assessment"]:
            self.assertIn(key, result)

    def test_scan_results_has_25_entries(self):
        result = scan_exponent_constraint()
        self.assertEqual(len(result["scan_results"]), 25)

    def test_sub_1000_pairs_has_2_entries(self):
        result = scan_exponent_constraint()
        self.assertEqual(result["n_sub_1000"], 2)

    def test_both_sub_1000_have_exponent_24(self):
        result = scan_exponent_constraint()
        for d, n in result["sub_1000_pairs"]:
            D = d + n
            self.assertEqual(d * D, 24,
                             f"(d={d},n={n}) has exponent {d*D}, expected 24")

    def test_assessment_nonempty(self):
        result = scan_exponent_constraint()
        self.assertIsInstance(result["assessment"], str)
        self.assertGreater(len(result["assessment"]), 0)


class TestScanVolumeCancellation(unittest.TestCase):
    """Tests for scan_volume_cancellation."""

    def test_returns_dict(self):
        result = scan_volume_cancellation()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = scan_volume_cancellation()
        for key in ["scan_results", "unique_n", "assessment"]:
            self.assertIn(key, result)

    def test_unique_n_is_2(self):
        result = scan_volume_cancellation()
        self.assertEqual(result["unique_n"], 2)

    def test_scan_results_has_6_entries(self):
        result = scan_volume_cancellation()
        self.assertEqual(len(result["scan_results"]), 6)

    def test_only_one_gives_product_1(self):
        result = scan_volume_cancellation()
        count = sum(1 for entry in result["scan_results"]
                    if abs(entry["product"] - 1.0) < 1e-12)
        self.assertEqual(count, 1)


class TestScanVacuumPolynomial(unittest.TestCase):
    """Tests for scan_vacuum_polynomial."""

    def test_returns_dict(self):
        result = scan_vacuum_polynomial()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = scan_vacuum_polynomial()
        for key in ["scan_results", "phi_pairs", "assessment"]:
            self.assertIn(key, result)

    def test_phi_pairs_includes_4_2(self):
        result = scan_vacuum_polynomial()
        self.assertIn((4, 2), result["phi_pairs"])

    def test_phi_pairs_excludes_3_5(self):
        result = scan_vacuum_polynomial()
        self.assertNotIn((3, 5), result["phi_pairs"])

    def test_d4_D6_involves_phi(self):
        result = scan_vacuum_polynomial()
        entry = next(e for e in result["scan_results"]
                     if e["d"] == 4 and e["n"] == 2)
        self.assertTrue(entry["involves_phi"])

    def test_d3_D8_does_not_involve_phi(self):
        result = scan_vacuum_polynomial()
        entry = next(e for e in result["scan_results"]
                     if e["d"] == 3 and e["n"] == 5)
        self.assertFalse(entry["involves_phi"])


class TestProveUniqueness(unittest.TestCase):
    """Tests for prove_uniqueness."""

    def test_returns_dict(self):
        result = prove_uniqueness()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = prove_uniqueness()
        for key in ["exponent_constraint", "volume_constraint",
                     "polynomial_constraint", "intersection",
                     "unique", "unique_pair", "assessment"]:
            self.assertIn(key, result)

    def test_unique_is_true(self):
        result = prove_uniqueness()
        self.assertTrue(result["unique"])

    def test_unique_pair_is_d4_n2_D6(self):
        result = prove_uniqueness()
        self.assertEqual(result["unique_pair"]["d"], 4)
        self.assertEqual(result["unique_pair"]["n"], 2)
        self.assertEqual(result["unique_pair"]["D"], 6)

    def test_intersection_has_exactly_1_element(self):
        result = prove_uniqueness()
        self.assertEqual(len(result["intersection"]), 1)

    def test_assessment_nonempty(self):
        result = prove_uniqueness()
        self.assertIsInstance(result["assessment"], str)
        self.assertGreater(len(result["assessment"]), 0)


class TestDeriveC3PhiHalf(unittest.TestCase):
    """Tests for derive_c3_phi_half."""

    def test_returns_dict(self):
        result = derive_c3_phi_half()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = derive_c3_phi_half()
        for key in ["c2", "c2_origin", "c3_phi_half", "c3_exact",
                     "c3_residual_percent", "F_with_phi_half", "assessment"]:
            self.assertIn(key, result)

    def test_c3_phi_half_close_to_0_809(self):
        result = derive_c3_phi_half()
        self.assertAlmostEqual(result["c3_phi_half"], 0.809, delta=0.001)

    def test_c3_residual_percent_under_1(self):
        result = derive_c3_phi_half()
        self.assertLess(result["c3_residual_percent"], 1.0)

    def test_c2_is_3(self):
        result = derive_c3_phi_half()
        self.assertEqual(result["c2"], 3)

    def test_c2_origin_mentions_spatial_or_d_minus_1(self):
        result = derive_c3_phi_half()
        origin = result["c2_origin"].lower()
        self.assertTrue("spatial" in origin or "d-1" in origin or "d - 1" in origin,
                        f"c2_origin should mention spatial or d-1, got: {result['c2_origin']}")


class TestPredictMuUnified(unittest.TestCase):
    """Tests for predict_mu_unified."""

    def test_returns_dict(self):
        result = predict_mu_unified()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = predict_mu_unified()
        for key in ["F", "mu_predicted", "mu_measured",
                     "residual_ppm", "sigma_tension", "assessment"]:
            self.assertIn(key, result)

    def test_mu_predicted_in_range(self):
        result = predict_mu_unified()
        self.assertGreater(result["mu_predicted"], 1835)
        self.assertLess(result["mu_predicted"], 1837)

    def test_residual_ppm_magnitude_under_0_01(self):
        result = predict_mu_unified()
        self.assertLess(abs(result["residual_ppm"]), 0.01)

    def test_sigma_tension_under_5(self):
        result = predict_mu_unified()
        self.assertLess(result["sigma_tension"], 5)

    def test_assessment_nonempty(self):
        result = predict_mu_unified()
        self.assertIsInstance(result["assessment"], str)
        self.assertGreater(len(result["assessment"]), 0)


class TestPredictGUnified(unittest.TestCase):
    """Tests for predict_G_unified."""

    def test_returns_dict(self):
        result = predict_G_unified()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = predict_G_unified()
        for key in ["G_unified", "G_mu_structure", "G_measured",
                     "residual_unified_ppm", "residual_mu_ppm",
                     "difference_ppm", "assessment"]:
            self.assertIn(key, result)

    def test_G_unified_in_range(self):
        result = predict_G_unified()
        self.assertGreater(result["G_unified"], 6.6e-11)
        self.assertLess(result["G_unified"], 6.8e-11)

    def test_residual_unified_magnitude_under_1_ppm(self):
        result = predict_G_unified()
        self.assertLess(abs(result["residual_unified_ppm"]), 1.0)

    def test_difference_ppm_under_0_1(self):
        result = predict_G_unified()
        self.assertLess(result["difference_ppm"], 0.1)

    def test_assessment_nonempty(self):
        result = predict_G_unified()
        self.assertIsInstance(result["assessment"], str)
        self.assertGreater(len(result["assessment"]), 0)


class TestComputeCompleteFormula(unittest.TestCase):
    """Tests for compute_complete_formula."""

    def test_returns_dict(self):
        result = compute_complete_formula()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_complete_formula()
        for key in ["formula_components", "n_fitted_params",
                     "n_derived", "G_predicted", "residual_ppm", "assessment"]:
            self.assertIn(key, result)

    def test_n_fitted_params_is_0(self):
        result = compute_complete_formula()
        self.assertEqual(result["n_fitted_params"], 0)

    def test_formula_components_has_5_plus_entries(self):
        result = compute_complete_formula()
        self.assertGreaterEqual(len(result["formula_components"]), 5)

    def test_G_predicted_in_range(self):
        result = compute_complete_formula()
        self.assertGreater(result["G_predicted"], 6.6e-11)
        self.assertLess(result["G_predicted"], 6.8e-11)

    def test_residual_ppm_magnitude_under_1(self):
        result = compute_complete_formula()
        self.assertLess(abs(result["residual_ppm"]), 1.0)


class TestSummarize(unittest.TestCase):
    """Tests for summarize_dimension_uniqueness."""

    def test_returns_dict(self):
        result = summarize_dimension_uniqueness()
        self.assertIsInstance(result, dict)

    def test_all_sub_keys_present(self):
        result = summarize_dimension_uniqueness()
        for key in ["exponent_scan", "volume_scan", "polynomial_scan",
                     "uniqueness_proof", "c3_derivation", "mu_prediction",
                     "G_prediction", "complete_formula",
                     "key_finding", "honest_assessment"]:
            self.assertIn(key, result)

    def test_key_finding_nonempty(self):
        result = summarize_dimension_uniqueness()
        self.assertIsInstance(result["key_finding"], str)
        self.assertGreater(len(result["key_finding"]), 0)

    def test_honest_assessment_nonempty(self):
        result = summarize_dimension_uniqueness()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 0)

    def test_sub_results_are_dicts(self):
        result = summarize_dimension_uniqueness()
        for key in ["exponent_scan", "volume_scan", "polynomial_scan",
                     "uniqueness_proof", "c3_derivation", "mu_prediction",
                     "G_prediction", "complete_formula"]:
            self.assertIsInstance(result[key], dict)

    def test_uniqueness_confirmed(self):
        result = summarize_dimension_uniqueness()
        self.assertTrue(result["uniqueness_proof"]["unique"])


if __name__ == "__main__":
    unittest.main()
