"""Tests for alpha_ladder_core.fifth_force_predictions module.

Covers all 8 public functions with ~50 focused unit tests verifying
return structure, physical constraints, and numerical predictions.
"""

import math
import unittest

from alpha_ladder_core.fifth_force_predictions import (
    compute_yukawa_signal,
    compute_signal_vs_distance,
    compute_exclusion_map,
    compute_eot_wash_prediction,
    compute_casimir_overlap,
    compute_discovery_reach,
    compute_alpha_ladder_prediction_line,
    summarize_fifth_force_predictions,
)


# ---------------------------------------------------------------------------
# 1. Yukawa signal at a single distance
# ---------------------------------------------------------------------------

class TestYukawaSignal(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_yukawa_signal(1e-4, 1e-4)
        self.assertIsInstance(result, dict)
        for key in ("r_meters", "lambda_m", "alpha",
                     "fractional_deviation", "exponential_factor", "detectable"):
            self.assertIn(key, result)

    def test_at_r_equals_lambda(self):
        """At r = lambda the signal should be alpha * 2 * exp(-1) ~ 0.455."""
        lam = 1e-4
        result = compute_yukawa_signal(lam, lam)
        expected = 0.618 * 2.0 * math.exp(-1.0)
        self.assertAlmostEqual(result["fractional_deviation"], expected, places=3)

    def test_at_small_r(self):
        """For r << lambda the signal approaches alpha = 0.618."""
        result = compute_yukawa_signal(1e-9, 1e-3)
        self.assertAlmostEqual(result["fractional_deviation"], 0.618, places=2)

    def test_at_large_r(self):
        """For r >> lambda the signal is exponentially suppressed."""
        result = compute_yukawa_signal(1.0, 1e-6)
        self.assertLess(result["fractional_deviation"], 1e-10)

    def test_alpha_is_0618(self):
        result = compute_yukawa_signal(1e-4, 1e-4)
        self.assertEqual(result["alpha"], 0.618)

    def test_detectable_true_at_small_r(self):
        result = compute_yukawa_signal(1e-9, 1e-3)
        self.assertTrue(result["detectable"])

    def test_fractional_deviation_positive(self):
        result = compute_yukawa_signal(5e-5, 1e-4)
        self.assertGreater(result["fractional_deviation"], 0.0)


# ---------------------------------------------------------------------------
# 2. Signal vs distance profile
# ---------------------------------------------------------------------------

class TestSignalVsDistance(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_signal_vs_distance(50e-6)
        self.assertIsInstance(result, dict)
        for key in ("r_meters", "fractional_deviation", "lambda_m",
                     "alpha", "peak_distance", "peak_signal"):
            self.assertIn(key, result)

    def test_default_100_points(self):
        result = compute_signal_vs_distance(50e-6)
        self.assertEqual(len(result["r_meters"]), 100)
        self.assertEqual(len(result["fractional_deviation"]), 100)

    def test_peak_signal_near_alpha(self):
        """Peak signal should be close to alpha (0.618), reached at r << lambda."""
        result = compute_signal_vs_distance(50e-6, r_min=1e-9, r_max=1e-3)
        self.assertGreater(result["peak_signal"], 0.5)
        self.assertLessEqual(result["peak_signal"], 0.618 + 0.001)

    def test_deviations_all_non_negative(self):
        result = compute_signal_vs_distance(50e-6)
        for d in result["fractional_deviation"]:
            self.assertGreaterEqual(d, 0.0)

    def test_deviations_decrease_at_large_r(self):
        result = compute_signal_vs_distance(50e-6, r_min=1e-6, r_max=10.0)
        devs = result["fractional_deviation"]
        # Last 5 entries should be smaller than the first 5 entries (on average)
        self.assertLess(sum(devs[-5:]), sum(devs[:5]))

    def test_peak_distance_near_lambda(self):
        """For default range starting at 1e-6, peak should be at r_min (smallest r)."""
        lam = 50e-6
        result = compute_signal_vs_distance(lam, r_min=1e-6, r_max=1.0)
        # Peak is at smallest r since signal = alpha*(1+r/lam)*exp(-r/lam)
        # is monotonically decreasing for r > 0
        self.assertLess(result["peak_distance"], lam)


# ---------------------------------------------------------------------------
# 3. Exclusion map
# ---------------------------------------------------------------------------

class TestExclusionMap(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_exclusion_map()
        self.assertIsInstance(result, dict)

    def test_seven_experiments(self):
        result = compute_exclusion_map()
        self.assertEqual(len(result["experiments"]), 7)

    def test_alpha_ladder_line_0618(self):
        result = compute_exclusion_map()
        self.assertEqual(result["alpha_ladder_line"], 0.618)

    def test_survival_window_exists(self):
        result = compute_exclusion_map()
        window = result["survival_window_m"]
        self.assertIsNotNone(window)
        self.assertEqual(len(window), 2)
        self.assertLess(window[0], window[1])

    def test_survival_lower_bound_near_30um(self):
        result = compute_exclusion_map()
        lower_um = result["survival_window_um"][0]
        self.assertGreaterEqual(lower_um, 20.0)
        self.assertLessEqual(lower_um, 60.0)

    def test_survival_upper_bound_near_71um(self):
        result = compute_exclusion_map()
        upper_um = result["survival_window_um"][1]
        self.assertGreaterEqual(upper_um, 50.0)
        self.assertLessEqual(upper_um, 100.0)

    def test_some_experiments_exclude(self):
        result = compute_exclusion_map()
        self.assertGreater(result["n_experiments_excluding"], 0)

    def test_some_experiments_dont_exclude(self):
        result = compute_exclusion_map()
        self.assertGreater(result["n_experiments_not_excluding"], 0)

    def test_has_honest_assessment(self):
        result = compute_exclusion_map()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 50)


# ---------------------------------------------------------------------------
# 4. Eot-Wash torsion balance prediction
# ---------------------------------------------------------------------------

class TestEotWashPrediction(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_eot_wash_prediction()
        self.assertIsInstance(result, dict)

    def test_default_a0_30um(self):
        result = compute_eot_wash_prediction()
        self.assertEqual(result["a_0_m"], 30e-6)

    def test_gap_distance_52um(self):
        result = compute_eot_wash_prediction()
        self.assertAlmostEqual(result["gap_distance_m"], 52e-6, places=10)

    def test_signal_between_01_and_05(self):
        result = compute_eot_wash_prediction(a_0_m=30e-6)
        self.assertGreater(result["predicted_signal"], 0.1)
        self.assertLess(result["predicted_signal"], 0.5)

    def test_detectable_true(self):
        result = compute_eot_wash_prediction(a_0_m=30e-6)
        self.assertTrue(result["detectable"])

    def test_signal_to_sensitivity_greater_than_1(self):
        result = compute_eot_wash_prediction(a_0_m=30e-6)
        self.assertGreater(result["signal_to_sensitivity_ratio"], 1.0)

    def test_has_signal_description(self):
        result = compute_eot_wash_prediction()
        self.assertIsInstance(result["signal_description"], str)
        self.assertGreater(len(result["signal_description"]), 20)

    def test_has_experiment_params(self):
        result = compute_eot_wash_prediction()
        params = result["experiment_params"]
        self.assertIsInstance(params, dict)
        self.assertIn("experiment", params)
        self.assertIn("reference", params)


# ---------------------------------------------------------------------------
# 5. Casimir overlap analysis
# ---------------------------------------------------------------------------

class TestCasimirOverlap(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_casimir_overlap()
        self.assertIsInstance(result, dict)

    def test_casimir_pressure_positive(self):
        result = compute_casimir_overlap()
        self.assertGreater(result["casimir_pressure_Pa"], 0.0)

    def test_yukawa_pressure_positive(self):
        result = compute_casimir_overlap()
        self.assertGreater(result["yukawa_pressure_Pa"], 0.0)

    def test_casimir_dominates_true(self):
        result = compute_casimir_overlap()
        self.assertTrue(result["casimir_dominates"])

    def test_ratio_much_less_than_1(self):
        result = compute_casimir_overlap()
        self.assertLess(result["ratio_yukawa_to_casimir"], 1.0)

    def test_has_honest_assessment(self):
        result = compute_casimir_overlap()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 50)


# ---------------------------------------------------------------------------
# 6. Discovery reach across experiments
# ---------------------------------------------------------------------------

class TestDiscoveryReach(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_discovery_reach()
        self.assertIsInstance(result, dict)

    def test_default_six_experiments(self):
        result = compute_discovery_reach()
        self.assertEqual(len(result["experiments"]), 6)

    def test_best_experiment_not_none(self):
        result = compute_discovery_reach()
        self.assertIsNotNone(result["best_experiment"])

    def test_all_have_signal_at_optimal(self):
        """At optimal lambda, signal = alpha * 2 * exp(-1) ~ 0.455."""
        result = compute_discovery_reach()
        expected = 0.618 * 2.0 * math.exp(-1.0)
        for exp in result["experiments"]:
            self.assertAlmostEqual(exp["signal_at_optimal"], expected, places=3)

    def test_eot_wash_detectable_at_optimal(self):
        result = compute_discovery_reach()
        eot_wash = [e for e in result["experiments"]
                    if "Eot-Wash 2006" in e["name"]][0]
        self.assertTrue(eot_wash["detectable_at_optimal"])

    def test_casimir_gap_0_2um(self):
        result = compute_discovery_reach()
        casimir = [e for e in result["experiments"]
                   if "Casimir" in e["name"]][0]
        self.assertAlmostEqual(casimir["gap_m"], 0.2e-6, places=10)

    def test_overall_range_exists(self):
        result = compute_discovery_reach()
        self.assertIsNotNone(result["overall_discovery_range_m"])
        rng = result["overall_discovery_range_m"]
        self.assertLess(rng[0], rng[1])

    def test_has_honest_assessment(self):
        result = compute_discovery_reach()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 50)


# ---------------------------------------------------------------------------
# 7. Prediction line for exclusion plots
# ---------------------------------------------------------------------------

class TestPredictionLine(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_alpha_ladder_prediction_line()
        self.assertIsInstance(result, dict)

    def test_default_200_points(self):
        result = compute_alpha_ladder_prediction_line()
        self.assertEqual(len(result["a_0_values"]), 200)

    def test_all_alpha_0618(self):
        result = compute_alpha_ladder_prediction_line()
        for a in result["alpha_values"]:
            self.assertEqual(a, 0.618)

    def test_lambda_equals_a0(self):
        result = compute_alpha_ladder_prediction_line()
        for a0, lam in zip(result["a_0_values"], result["lambda_values"]):
            self.assertEqual(a0, lam)

    def test_m_phi_positive(self):
        result = compute_alpha_ladder_prediction_line()
        for m in result["m_phi_eV_values"]:
            self.assertGreater(m, 0.0)

    def test_has_survival_window_flags(self):
        result = compute_alpha_ladder_prediction_line()
        flags = result["in_survival_window"]
        self.assertIn(True, flags)
        self.assertIn(False, flags)

    def test_eot_wash_signal_values_length(self):
        result = compute_alpha_ladder_prediction_line()
        self.assertEqual(len(result["eot_wash_signal_values"]),
                         len(result["a_0_values"]))


# ---------------------------------------------------------------------------
# 8. Summary / dashboard entry point
# ---------------------------------------------------------------------------

class TestSummarize(unittest.TestCase):

    def test_returns_dict(self):
        result = summarize_fifth_force_predictions()
        self.assertIsInstance(result, dict)

    def test_has_exclusion_map(self):
        result = summarize_fifth_force_predictions()
        self.assertIn("exclusion_map", result)
        self.assertIsInstance(result["exclusion_map"], dict)

    def test_has_eot_wash_prediction(self):
        result = summarize_fifth_force_predictions()
        self.assertIn("eot_wash_prediction", result)
        self.assertIsInstance(result["eot_wash_prediction"], dict)

    def test_has_casimir_overlap(self):
        result = summarize_fifth_force_predictions()
        self.assertIn("casimir_overlap", result)
        self.assertIsInstance(result["casimir_overlap"], dict)

    def test_has_discovery_reach(self):
        result = summarize_fifth_force_predictions()
        self.assertIn("discovery_reach", result)
        self.assertIsInstance(result["discovery_reach"], dict)

    def test_has_key_prediction(self):
        result = summarize_fifth_force_predictions()
        self.assertIsInstance(result["key_prediction"], str)
        self.assertGreater(len(result["key_prediction"]), 20)

    def test_has_survival_window_um(self):
        result = summarize_fifth_force_predictions()
        window = result["survival_window_um"]
        self.assertIsNotNone(window)
        self.assertEqual(len(window), 2)

    def test_has_optimal_experiment(self):
        result = summarize_fifth_force_predictions()
        self.assertIsNotNone(result["optimal_experiment"])
        self.assertIsInstance(result["optimal_experiment"], str)


if __name__ == "__main__":
    unittest.main()
