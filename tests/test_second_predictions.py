"""Tests for alpha_ladder_core.second_predictions module.

Covers all 7 public functions with ~45 focused unit tests verifying
return structure, physical constraints, and numerical predictions.
"""

import math
import unittest

from alpha_ladder_core.second_predictions import (
    compute_time_variation_coefficients,
    compute_current_bounds,
    compute_m6_from_a0,
    scan_a0_m6_landscape,
    analyze_correlated_variations,
    analyze_mu_tension,
    summarize_second_predictions,
)


# ---------------------------------------------------------------------------
# 1. Time variation coefficients
# ---------------------------------------------------------------------------

class TestTimeVariationCoefficients(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_time_variation_coefficients()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_time_variation_coefficients()
        for key in ("A_coefficient", "B_coefficient", "formula_description",
                     "comparison_with_other_frameworks", "assessment"):
            self.assertIn(key, result)

    def test_A_coefficient_close_to_24(self):
        result = compute_time_variation_coefficients()
        self.assertAlmostEqual(result["A_coefficient"], 24.0, delta=0.01)

    def test_B_coefficient_close_to_2(self):
        result = compute_time_variation_coefficients()
        self.assertAlmostEqual(result["B_coefficient"], 2.0, delta=0.01)

    def test_comparison_has_five_frameworks(self):
        result = compute_time_variation_coefficients()
        comp = result["comparison_with_other_frameworks"]
        self.assertGreaterEqual(len(comp), 5)

    def test_assessment_nonempty(self):
        result = compute_time_variation_coefficients()
        self.assertTrue(len(result["assessment"]) > 0)


# ---------------------------------------------------------------------------
# 2. Current bounds
# ---------------------------------------------------------------------------

class TestCurrentBounds(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_current_bounds()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_current_bounds()
        for key in ("alpha_bound", "mu_bound", "G_bound",
                     "predicted_dG_from_alpha_bound",
                     "ratio_to_current_sensitivity", "testable"):
            self.assertIn(key, result)

    def test_predicted_dG_less_than_G_bound(self):
        result = compute_current_bounds()
        self.assertLess(result["predicted_dG_from_alpha_bound"],
                        result["G_bound"])

    def test_not_currently_testable(self):
        result = compute_current_bounds()
        self.assertFalse(result["testable"])

    def test_ratio_below_one(self):
        result = compute_current_bounds()
        self.assertLess(result["ratio_to_current_sensitivity"], 1.0)


# ---------------------------------------------------------------------------
# 3. M_6 from a_0
# ---------------------------------------------------------------------------

class TestComputeM6FromA0(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_m6_from_a0(50e-6)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_m6_from_a0(50e-6)
        for key in ("a0_meters", "M_6_kg", "M_6_eV", "M_6_TeV",
                     "ratio_to_electron_mass", "ratio_to_proton_mass",
                     "nearest_known_scale", "alpha_power"):
            self.assertIn(key, result)

    def test_M6_positive(self):
        result = compute_m6_from_a0(50e-6)
        self.assertGreater(result["M_6_TeV"], 0.0)

    def test_M6_at_30um_between_3_and_6_TeV(self):
        result = compute_m6_from_a0(30e-6)
        self.assertGreater(result["M_6_TeV"], 3.0)
        self.assertLess(result["M_6_TeV"], 6.0)

    def test_M6_at_planck_length_near_planck_mass(self):
        """At a_0 = l_Pl, M_6 should be close to M_Pl energy."""
        from alpha_ladder_core.constants import get_constants
        c = get_constants("CODATA 2018")
        hbar = float(c.hbar)
        cc = float(c.c)
        G = float(c.G)
        l_Pl = math.sqrt(hbar * G / cc ** 3)
        result = compute_m6_from_a0(l_Pl)
        # M_Pl ~ 1.22e19 GeV = 1.22e16 TeV
        # M_6 at l_Pl should be within an order of magnitude of M_Pl
        M_Pl_TeV = 1.22e16
        log_ratio = abs(math.log10(result["M_6_TeV"] / M_Pl_TeV))
        self.assertLess(log_ratio, 1.5)

    def test_alpha_power_is_finite(self):
        result = compute_m6_from_a0(50e-6)
        self.assertTrue(math.isfinite(result["alpha_power"]))


# ---------------------------------------------------------------------------
# 4. Scan a_0 - M_6 landscape
# ---------------------------------------------------------------------------

class TestScanA0M6Landscape(unittest.TestCase):

    def test_returns_dict(self):
        result = scan_a0_m6_landscape()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = scan_a0_m6_landscape()
        for key in ("scan_results", "survival_window", "assessment"):
            self.assertIn(key, result)

    def test_scan_has_enough_entries(self):
        result = scan_a0_m6_landscape()
        self.assertGreaterEqual(len(result["scan_results"]), 15)

    def test_survival_window_has_required_keys(self):
        result = scan_a0_m6_landscape()
        sw = result["survival_window"]
        for key in ("a0_min", "a0_max", "M6_min_TeV", "M6_max_TeV"):
            self.assertIn(key, sw)

    def test_survival_M6_min_positive(self):
        result = scan_a0_m6_landscape()
        self.assertGreater(result["survival_window"]["M6_min_TeV"], 0.0)


# ---------------------------------------------------------------------------
# 5. Correlated variations
# ---------------------------------------------------------------------------

class TestAnalyzeCorrelatedVariations(unittest.TestCase):

    def test_returns_dict(self):
        result = analyze_correlated_variations()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = analyze_correlated_variations()
        for key in ("dalpha_alpha", "predicted_dG_G",
                     "predicted_dG_G_with_qcd", "R_qcd",
                     "violates_llr_bound", "assessment"):
            self.assertIn(key, result)

    def test_predicted_dG_G_positive_for_positive_dalpha(self):
        result = analyze_correlated_variations(dalpha_alpha=1e-17)
        self.assertGreater(result["predicted_dG_G"], 0.0)

    def test_qcd_enhanced_larger_than_alpha_only(self):
        result = analyze_correlated_variations(dalpha_alpha=1e-17)
        self.assertGreater(abs(result["predicted_dG_G_with_qcd"]),
                           abs(result["predicted_dG_G"]))

    def test_does_not_violate_llr_at_1e17(self):
        result = analyze_correlated_variations(dalpha_alpha=1e-17)
        self.assertFalse(result["violates_llr_bound"])


# ---------------------------------------------------------------------------
# 6. Mu tension
# ---------------------------------------------------------------------------

class TestAnalyzeMuTension(unittest.TestCase):

    def test_returns_dict(self):
        result = analyze_mu_tension()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = analyze_mu_tension()
        for key in ("mu_measured", "mu_predicted_from_bridge",
                     "tension_ppm", "tension_sigma",
                     "interpretation", "possible_resolutions"):
            self.assertIn(key, result)

    def test_mu_measured_close_to_1836(self):
        result = analyze_mu_tension()
        self.assertAlmostEqual(result["mu_measured"], 1836.15, delta=0.1)

    def test_tension_sigma_above_10(self):
        result = analyze_mu_tension()
        self.assertGreater(result["tension_sigma"], 10.0)

    def test_mu_predicted_differs_from_measured(self):
        result = analyze_mu_tension()
        # They differ at the ~ppm level (4th-5th decimal place)
        self.assertNotAlmostEqual(result["mu_predicted_from_bridge"],
                                   result["mu_measured"], places=5)

    def test_possible_resolutions_has_items(self):
        result = analyze_mu_tension()
        self.assertIsInstance(result["possible_resolutions"], list)
        self.assertGreaterEqual(len(result["possible_resolutions"]), 2)


# ---------------------------------------------------------------------------
# 7. Summary
# ---------------------------------------------------------------------------

class TestSummarizeSecondPredictions(unittest.TestCase):

    def test_returns_dict(self):
        result = summarize_second_predictions()
        self.assertIsInstance(result, dict)

    def test_all_sub_keys_present(self):
        result = summarize_second_predictions()
        for key in ("time_variation", "current_bounds", "m6_landscape",
                     "correlated_variations", "mu_tension",
                     "strongest_prediction", "honest_assessment"):
            self.assertIn(key, result)

    def test_strongest_prediction_nonempty(self):
        result = summarize_second_predictions()
        self.assertTrue(len(result["strongest_prediction"]) > 0)

    def test_honest_assessment_nonempty(self):
        result = summarize_second_predictions()
        self.assertTrue(len(result["honest_assessment"]) > 0)

    def test_sub_results_are_dicts(self):
        result = summarize_second_predictions()
        for key in ("time_variation", "current_bounds", "m6_landscape",
                     "correlated_variations", "mu_tension"):
            self.assertIsInstance(result[key], dict)

    def test_time_variation_in_result(self):
        result = summarize_second_predictions()
        self.assertIn("A_coefficient", result["time_variation"])

    def test_mu_tension_in_result(self):
        result = summarize_second_predictions()
        self.assertIn("tension_sigma", result["mu_tension"])


if __name__ == "__main__":
    unittest.main()
