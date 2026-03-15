"""Tests for alpha_ladder_core.testable_predictions module.

Covers all 6 public functions with ~45 focused unit tests verifying
return structure, physical constraints, and numerical predictions.
"""

import math
import unittest

from alpha_ladder_core.testable_predictions import (
    predict_G_precision,
    predict_mu_from_consistency,
    predict_cosmological_constant,
    predict_fifth_force_signal,
    predict_G_multiple_experiments,
    summarize_testable_predictions,
)


# ---------------------------------------------------------------------------
# 1. G precision prediction from mu-structure formula
# ---------------------------------------------------------------------------

class TestPredictGPrecision(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.result = predict_G_precision()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {
            "G_predicted", "G_measured", "residual_ppm",
            "within_current_uncertainty", "predicted_uncertainty_ppb",
            "current_uncertainty_ppb", "improvement_factor",
            "G_prediction_string", "assessment",
        }
        for key in required:
            self.assertIn(key, self.result, f"Missing key: {key}")

    def test_G_predicted_reasonable(self):
        """G_predicted should lie between 6.67e-11 and 6.68e-11."""
        G = self.result["G_predicted"]
        self.assertGreater(G, 6.67e-11)
        self.assertLess(G, 6.68e-11)

    def test_residual_within_uncertainty(self):
        """The prediction should be within current experimental uncertainty."""
        self.assertTrue(self.result["within_current_uncertainty"])

    def test_residual_magnitude(self):
        """|residual_ppm| should be between 1 and 20."""
        mag = abs(self.result["residual_ppm"])
        self.assertGreater(mag, 1)
        self.assertLess(mag, 20)

    def test_predicted_uncertainty_sub_10ppb(self):
        """Predicted uncertainty should be below 10 ppb."""
        self.assertLess(self.result["predicted_uncertainty_ppb"], 10)

    def test_improvement_factor_large(self):
        """Improvement factor over current measurements should exceed 1000."""
        self.assertGreater(self.result["improvement_factor"], 1000)

    def test_prediction_string_nonempty(self):
        s = self.result["G_prediction_string"]
        self.assertIsInstance(s, str)
        self.assertTrue(len(s) > 0)

    def test_assessment_nonempty(self):
        s = self.result["assessment"]
        self.assertIsInstance(s, str)
        self.assertTrue(len(s) > 0)


# ---------------------------------------------------------------------------
# 2. mu consistency prediction
# ---------------------------------------------------------------------------

class TestPredictMuConsistency(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.result = predict_mu_from_consistency()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {
            "mu_predicted", "mu_measured", "residual_ppm",
            "sigma_tension", "falsified", "interpretation", "assessment",
        }
        for key in required:
            self.assertIn(key, self.result, f"Missing key: {key}")

    def test_mu_predicted_close(self):
        """mu_predicted should be between 1836.0 and 1836.5."""
        mu = self.result["mu_predicted"]
        self.assertGreater(mu, 1836.0)
        self.assertLess(mu, 1836.5)

    def test_residual_positive(self):
        """residual_ppm should be positive (predicted > measured)."""
        self.assertGreater(self.result["residual_ppm"], 0)

    def test_residual_small(self):
        """|residual_ppm| should be less than 10."""
        self.assertLess(abs(self.result["residual_ppm"]), 10)

    def test_sigma_tension_large(self):
        """sigma_tension should exceed 1000 (extremely precise measurement)."""
        self.assertGreater(self.result["sigma_tension"], 1000)

    def test_falsified(self):
        """Consistency check should be falsified (mu too precisely known)."""
        self.assertTrue(self.result["falsified"])

    def test_assessment_nonempty(self):
        s = self.result["assessment"]
        self.assertIsInstance(s, str)
        self.assertTrue(len(s) > 0)


# ---------------------------------------------------------------------------
# 3. Cosmological constant prediction
# ---------------------------------------------------------------------------

class TestPredictCosmologicalConstant(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.result = predict_cosmological_constant()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {
            "Lambda_predicted", "Lambda_observed", "ratio",
            "log10_ratio", "naive_discrepancy_orders",
            "our_discrepancy_orders", "orders_improvement",
            "missing_factor", "assessment",
        }
        for key in required:
            self.assertIn(key, self.result, f"Missing key: {key}")

    def test_Lambda_predicted_positive(self):
        self.assertGreater(self.result["Lambda_predicted"], 0)

    def test_Lambda_observed_positive(self):
        self.assertGreater(self.result["Lambda_observed"], 0)

    def test_ratio_order_one(self):
        """The ratio Lambda_predicted/Lambda_observed should be O(1)."""
        ratio = self.result["ratio"]
        self.assertGreater(ratio, 0.1)
        self.assertLess(ratio, 10)

    def test_log10_ratio_small(self):
        """|log10_ratio| should be less than 1 (within an order of magnitude)."""
        self.assertLess(abs(self.result["log10_ratio"]), 1)

    def test_orders_improvement_huge(self):
        """Improvement over naive QFT prediction should exceed 100 orders."""
        self.assertGreater(self.result["orders_improvement"], 100)

    def test_assessment_nonempty(self):
        s = self.result["assessment"]
        self.assertIsInstance(s, str)
        self.assertTrue(len(s) > 0)


# ---------------------------------------------------------------------------
# 4. Fifth force signal prediction
# ---------------------------------------------------------------------------

class TestPredictFifthForce(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.result = predict_fifth_force_signal()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {
            "dilaton_mass_eV", "force_range_m", "signal_strength",
            "observable", "testable_window_exists",
            "eot_wash_optimal_range_um", "eot_wash_optimal_mass_meV",
            "assessment",
        }
        for key in required:
            self.assertIn(key, self.result, f"Missing key: {key}")

    def test_dilaton_mass_positive(self):
        self.assertGreater(self.result["dilaton_mass_eV"], 0)

    def test_observable_false(self):
        """Planck-mass dilaton should not produce an observable signal."""
        self.assertFalse(self.result["observable"])

    def test_testable_window_exists(self):
        """A testable window should exist in principle."""
        self.assertTrue(self.result["testable_window_exists"])

    def test_eot_wash_range(self):
        """Eot-Wash optimal range should be between 10 and 100 um."""
        r = self.result["eot_wash_optimal_range_um"]
        self.assertGreater(r, 10)
        self.assertLess(r, 100)

    def test_assessment_nonempty(self):
        s = self.result["assessment"]
        self.assertIsInstance(s, str)
        self.assertTrue(len(s) > 0)


# ---------------------------------------------------------------------------
# 5. G comparison against multiple experiments
# ---------------------------------------------------------------------------

class TestPredictGExperiments(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.result = predict_G_multiple_experiments()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {
            "G_predicted", "comparisons", "best_agreement",
            "worst_agreement", "n_within_1sigma", "n_within_2sigma",
            "assessment",
        }
        for key in required:
            self.assertIn(key, self.result, f"Missing key: {key}")

    def test_seven_comparisons(self):
        """There should be exactly 7 experimental comparisons."""
        self.assertEqual(len(self.result["comparisons"]), 7)

    def test_each_comparison_keys(self):
        """Every comparison entry must have the required fields."""
        required_comp_keys = {"experiment", "G_exp", "G_unc", "sigma", "direction"}
        for comp in self.result["comparisons"]:
            for key in required_comp_keys:
                self.assertIn(key, comp, f"Missing key in comparison: {key}")

    def test_sigma_positive(self):
        """All sigma values must be positive."""
        for comp in self.result["comparisons"]:
            self.assertGreater(comp["sigma"], 0,
                               f"Non-positive sigma for {comp.get('experiment', '?')}")

    def test_some_within_2sigma(self):
        """At least one experiment should agree within 2 sigma."""
        self.assertGreaterEqual(self.result["n_within_2sigma"], 1)

    def test_best_sigma_less_than_worst(self):
        """Best-agreement sigma should be smaller than worst-agreement sigma."""
        best = self.result["best_agreement"]["sigma"]
        worst = self.result["worst_agreement"]["sigma"]
        self.assertLess(best, worst)

    def test_assessment_nonempty(self):
        s = self.result["assessment"]
        self.assertIsInstance(s, str)
        self.assertTrue(len(s) > 0)


# ---------------------------------------------------------------------------
# 6. Summarize all testable predictions
# ---------------------------------------------------------------------------

class TestSummarizePredictions(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.result = summarize_testable_predictions()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {
            "g_precision", "mu_consistency", "cosmological_constant",
            "fifth_force", "g_experiments", "predictions_summary",
            "key_finding", "honest_assessment",
        }
        for key in required:
            self.assertIn(key, self.result, f"Missing key: {key}")

    def test_predictions_summary_list(self):
        """predictions_summary should be a list with at least 4 entries."""
        summary = self.result["predictions_summary"]
        self.assertIsInstance(summary, list)
        self.assertGreaterEqual(len(summary), 4)

    def test_key_finding_nonempty(self):
        s = self.result["key_finding"]
        self.assertIsInstance(s, str)
        self.assertTrue(len(s) > 0)

    def test_honest_assessment_nonempty(self):
        s = self.result["honest_assessment"]
        self.assertIsInstance(s, str)
        self.assertTrue(len(s) > 0)


if __name__ == "__main__":
    unittest.main()
