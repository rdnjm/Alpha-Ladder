"""Tests for bridge_significance module."""

import unittest
from alpha_ladder_core.constants import get_constants
from alpha_ladder_core import bridge_significance as bs


class TestEnumerateExpressions(unittest.TestCase):
    """Tests for enumerate_all_expression_values."""

    @classmethod
    def setUpClass(cls):
        cls.data = bs.enumerate_all_expression_values()

    def test_total_count_reasonable(self):
        """Total should be in range 100K-200K."""
        self.assertGreater(self.data["total_count"], 100_000)
        self.assertLess(self.data["total_count"], 200_000)

    def test_unique_less_than_total(self):
        self.assertLess(len(self.data["unique_values"]), self.data["total_count"])

    def test_phase_counts_sum(self):
        pc = self.data["phase_counts"]
        self.assertEqual(
            pc["phase1"] + pc["phase2"] + pc["phase3"],
            self.data["total_count"],
        )

    def test_all_values_finite(self):
        import math
        for v in self.data["values"][:1000]:  # spot check first 1000
            self.assertTrue(math.isfinite(v))

    def test_all_values_positive(self):
        for v in self.data["values"][:1000]:
            self.assertGreater(v, 0)

    def test_phi_expressions_nonempty(self):
        self.assertGreater(self.data["phi_count"], 0)

    def test_labels_match_values(self):
        self.assertEqual(len(self.data["values"]), len(self.data["labels"]))


class TestComputeCoverage(unittest.TestCase):
    """Tests for compute_coverage."""

    @classmethod
    def setUpClass(cls):
        cls.data = bs.enumerate_all_expression_values()
        cls.unique = cls.data["unique_values"]

    def test_coverage_between_0_and_1(self):
        result = bs.compute_coverage(self.unique)
        self.assertGreaterEqual(result["coverage_fraction"], 0.0)
        self.assertLessEqual(result["coverage_fraction"], 1.0)

    def test_coverage_increases_with_epsilon(self):
        r1 = bs.compute_coverage(self.unique, epsilon=0.0001)
        r2 = bs.compute_coverage(self.unique, epsilon=0.001)
        self.assertLessEqual(r1["coverage_fraction"], r2["coverage_fraction"])

    def test_zero_epsilon_gives_zero(self):
        result = bs.compute_coverage(self.unique, epsilon=0.0)
        self.assertEqual(result["coverage_fraction"], 0.0)

    def test_single_value_known_result(self):
        """Single value at 1.0 in [0.5, 2.0] with eps=0.1 should cover 0.2/1.5."""
        result = bs.compute_coverage([1.0], interval_min=0.5, interval_max=2.0, epsilon=0.1)
        # Band: [0.9, 1.1], length 0.2, interval length 1.5
        expected = 0.2 / 1.5
        self.assertAlmostEqual(result["coverage_fraction"], expected, places=5)

    def test_covered_length_positive(self):
        result = bs.compute_coverage(self.unique)
        self.assertGreater(result["covered_length"], 0)

    def test_values_in_interval_count(self):
        result = bs.compute_coverage(self.unique)
        self.assertGreater(result["values_in_interval"], 0)

    def test_interval_length_correct(self):
        result = bs.compute_coverage(self.unique, interval_min=0.5, interval_max=2.0)
        self.assertAlmostEqual(result["interval_length"], 1.5, places=10)


class TestBonferroniPvalue(unittest.TestCase):
    """Tests for compute_bonferroni_pvalue."""

    def test_p_capped_at_1(self):
        result = bs.compute_bonferroni_pvalue(n_trials=1_000_000)
        self.assertLessEqual(result["p_corrected"], 1.0)

    def test_small_n_correct(self):
        result = bs.compute_bonferroni_pvalue(n_trials=1, epsilon=0.01)
        self.assertAlmostEqual(result["p_single"], 0.02)
        self.assertAlmostEqual(result["p_corrected"], 0.02)

    def test_large_n_saturates(self):
        result = bs.compute_bonferroni_pvalue(n_trials=135000, epsilon=0.00016)
        self.assertEqual(result["p_corrected"], 1.0)

    def test_significant_flag(self):
        result = bs.compute_bonferroni_pvalue(n_trials=10, epsilon=0.001)
        # p = 10 * 0.002 = 0.02 < 0.05
        self.assertTrue(result["significant_at_005"])


class TestEmpiricalPvalue(unittest.TestCase):
    """Tests for compute_empirical_pvalue."""

    @classmethod
    def setUpClass(cls):
        cls.data = bs.enumerate_all_expression_values()
        cls.unique = cls.data["unique_values"]

    def test_p_empirical_less_than_bonferroni(self):
        result = bs.compute_empirical_pvalue(self.unique)
        self.assertLessEqual(result["p_empirical"], result["p_bonferroni"])

    def test_interpretation_nonempty(self):
        result = bs.compute_empirical_pvalue(self.unique)
        self.assertGreater(len(result["interpretation"]), 0)

    def test_p_empirical_is_coverage(self):
        result = bs.compute_empirical_pvalue(self.unique)
        cov = bs.compute_coverage(self.unique)
        self.assertAlmostEqual(result["p_empirical"], cov["coverage_fraction"], places=10)

    def test_returns_all_keys(self):
        result = bs.compute_empirical_pvalue(self.unique)
        for key in ["p_empirical", "p_bonferroni", "interpretation"]:
            self.assertIn(key, result)


class TestMonteCarloValidation(unittest.TestCase):
    """Tests for run_monte_carlo_validation."""

    @classmethod
    def setUpClass(cls):
        cls.data = bs.enumerate_all_expression_values()
        cls.unique = cls.data["unique_values"]

    def test_reproducible_with_seed(self):
        r1 = bs.run_monte_carlo_validation(self.unique, n_samples=1000, seed=42)
        r2 = bs.run_monte_carlo_validation(self.unique, n_samples=1000, seed=42)
        self.assertEqual(r1["p_monte_carlo"], r2["p_monte_carlo"])

    def test_different_seeds_may_differ(self):
        r1 = bs.run_monte_carlo_validation(self.unique, n_samples=1000, seed=42)
        r2 = bs.run_monte_carlo_validation(self.unique, n_samples=1000, seed=99)
        # They CAN be equal but usually won't be - just check they both run
        self.assertIn("p_monte_carlo", r1)
        self.assertIn("p_monte_carlo", r2)

    def test_p_mc_between_0_and_1(self):
        result = bs.run_monte_carlo_validation(self.unique, n_samples=1000, seed=42)
        self.assertGreaterEqual(result["p_monte_carlo"], 0.0)
        self.assertLessEqual(result["p_monte_carlo"], 1.0)

    def test_consistent_with_coverage(self):
        result = bs.run_monte_carlo_validation(self.unique, n_samples=5000, seed=42)
        self.assertTrue(result["consistent"])

    def test_uncertainty_positive(self):
        result = bs.run_monte_carlo_validation(self.unique, n_samples=1000, seed=42)
        self.assertGreater(result["p_monte_carlo_uncertainty"], 0)

    def test_n_hits_reasonable(self):
        result = bs.run_monte_carlo_validation(self.unique, n_samples=1000, seed=42)
        self.assertGreaterEqual(result["n_hits"], 0)
        self.assertLessEqual(result["n_hits"], 1000)


class TestAnalyzePhiSubset(unittest.TestCase):
    """Tests for analyze_phi_subset."""

    @classmethod
    def setUpClass(cls):
        cls.data = bs.enumerate_all_expression_values()

    def test_phi_count_reasonable(self):
        result = bs.analyze_phi_subset(self.data)
        self.assertGreater(result["phi_count"], 100)

    def test_best_phi_match_is_phi_sq_over_2(self):
        result = bs.analyze_phi_subset(self.data, target_value=1.309)
        label = result["best_phi_match"]["label"].lower()
        self.assertTrue("phi" in label)

    def test_algebraic_note_mentions_sqrt5(self):
        result = bs.analyze_phi_subset(self.data)
        self.assertIn("sqrt(5)", result["algebraic_field_note"])

    def test_phi_coverage_between_0_and_1(self):
        result = bs.analyze_phi_subset(self.data)
        self.assertGreaterEqual(result["phi_coverage_fraction"], 0.0)
        self.assertLessEqual(result["phi_coverage_fraction"], 1.0)

    def test_residual_ppm_positive(self):
        result = bs.analyze_phi_subset(self.data, target_value=1.309)
        self.assertGreater(result["best_phi_match"]["residual_ppm"], 0)


class TestSummarize(unittest.TestCase):
    """Tests for summarize_bridge_significance."""

    @classmethod
    def setUpClass(cls):
        cls.constants = get_constants("CODATA 2018")
        cls.result = bs.summarize_bridge_significance(cls.constants)

    def test_all_keys_present(self):
        expected = [
            "search_space", "target_value", "residual_ppm",
            "coverage", "bonferroni", "empirical_pvalue",
            "monte_carlo", "phi_analysis", "histogram",
            "honest_assessment",
        ]
        for key in expected:
            self.assertIn(key, self.result)

    def test_target_around_1_3(self):
        self.assertGreater(self.result["target_value"], 1.0)
        self.assertLess(self.result["target_value"], 2.0)

    def test_residual_less_than_200_ppm(self):
        self.assertLess(self.result["residual_ppm"], 200)

    def test_honest_assessment_nonempty(self):
        self.assertGreater(len(self.result["honest_assessment"]), 50)

    def test_search_space_has_counts(self):
        ss = self.result["search_space"]
        self.assertIn("total_count", ss)
        self.assertIn("unique_count", ss)


if __name__ == "__main__":
    unittest.main()
