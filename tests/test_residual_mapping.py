"""
Tests for alpha_ladder_core/residual_mapping.py

Maps the residual delta = sqrt(phi) - k_exact against Standard Model constants.
"""

import math
import unittest

from alpha_ladder_core.residual_mapping import (
    compute_delta,
    scan_anomalous_magnetic_moment,
    scan_zeta_functions,
    scan_composite_expressions,
    scan_k_closed_forms,
    analyze_delta_structure,
    summarize_residual_mapping,
)


class TestComputeDelta(unittest.TestCase):
    """Tests for compute_delta."""

    def test_returns_dict(self):
        result = compute_delta()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_delta()
        for key in ("k_exact", "sqrt_phi", "delta", "delta_over_alpha", "mu", "alpha_g_measured"):
            self.assertIn(key, result)

    def test_k_exact_approx(self):
        """k_exact should be approximately 1.262."""
        result = compute_delta()
        self.assertAlmostEqual(result["k_exact"], 1.262, delta=0.01)

    def test_sqrt_phi_approx(self):
        """sqrt(phi) should be approximately 1.27202."""
        result = compute_delta()
        self.assertAlmostEqual(result["sqrt_phi"], 1.27202, delta=0.001)

    def test_delta_positive(self):
        """delta = sqrt(phi) - k_exact should be positive (sqrt(phi) > k_exact)."""
        result = compute_delta()
        self.assertGreater(result["delta"], 0)

    def test_delta_order_of_magnitude(self):
        """delta should be O(0.01), i.e., between 0.001 and 0.1."""
        result = compute_delta()
        self.assertGreater(result["delta"], 0.001)
        self.assertLess(result["delta"], 0.1)

    def test_mu_approx(self):
        """mu should be approximately 1836.15."""
        result = compute_delta()
        self.assertAlmostEqual(result["mu"], 1836.15, delta=0.1)

    def test_alpha_g_measured_order(self):
        """alpha_g_measured should be O(10^-45)."""
        result = compute_delta()
        self.assertGreater(result["alpha_g_measured"], 1e-46)
        self.assertLess(result["alpha_g_measured"], 1e-44)


class TestScanAnomalousMagneticMoment(unittest.TestCase):
    """Tests for scan_anomalous_magnetic_moment."""

    def test_returns_dict(self):
        result = scan_anomalous_magnetic_moment()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = scan_anomalous_magnetic_moment()
        for key in ("candidates", "best_match", "assessment"):
            self.assertIn(key, result)

    def test_candidates_nonempty(self):
        result = scan_anomalous_magnetic_moment()
        self.assertGreater(len(result["candidates"]), 0)

    def test_candidate_keys(self):
        result = scan_anomalous_magnetic_moment()
        for c in result["candidates"]:
            for key in ("name", "expression", "value", "ratio"):
                self.assertIn(key, c)

    def test_schwinger_term_present(self):
        """Should include the Schwinger alpha/(2*pi) term."""
        result = scan_anomalous_magnetic_moment()
        names = [c["name"] for c in result["candidates"]]
        self.assertTrue(any("Schwinger" in n for n in names))

    def test_assessment_nonempty(self):
        result = scan_anomalous_magnetic_moment()
        self.assertIsInstance(result["assessment"], str)
        self.assertGreater(len(result["assessment"]), 0)


class TestScanZetaFunctions(unittest.TestCase):
    """Tests for scan_zeta_functions."""

    def test_returns_dict(self):
        result = scan_zeta_functions()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = scan_zeta_functions()
        for key in ("candidates", "best_match", "assessment"):
            self.assertIn(key, result)

    def test_candidates_nonempty(self):
        result = scan_zeta_functions()
        self.assertGreater(len(result["candidates"]), 0)

    def test_zeta3_present(self):
        """Should include zeta(3) / Apery's constant."""
        result = scan_zeta_functions()
        names = [c["name"] for c in result["candidates"]]
        self.assertTrue(any("zeta(3)" in n for n in names))

    def test_euler_gamma_present(self):
        """Should include Euler-Mascheroni gamma."""
        result = scan_zeta_functions()
        names = [c["name"] for c in result["candidates"]]
        self.assertTrue(any("gamma" in n.lower() for n in names))

    def test_assessment_nonempty(self):
        result = scan_zeta_functions()
        self.assertIsInstance(result["assessment"], str)
        self.assertGreater(len(result["assessment"]), 0)


class TestScanCompositeExpressions(unittest.TestCase):
    """Tests for scan_composite_expressions."""

    def test_returns_dict(self):
        result = scan_composite_expressions()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = scan_composite_expressions()
        for key in ("candidates", "best_matches", "n_within_100ppm", "n_within_1000ppm", "assessment"):
            self.assertIn(key, result)

    def test_best_matches_sorted(self):
        """best_matches should be sorted by residual_ppm ascending."""
        result = scan_composite_expressions()
        ppms = [c["residual_ppm"] for c in result["best_matches"]]
        self.assertEqual(ppms, sorted(ppms))

    def test_best_matches_max_10(self):
        result = scan_composite_expressions()
        self.assertLessEqual(len(result["best_matches"]), 10)

    def test_n_within_counters_consistent(self):
        """n_within_100ppm <= n_within_1000ppm."""
        result = scan_composite_expressions()
        self.assertLessEqual(result["n_within_100ppm"], result["n_within_1000ppm"])

    def test_candidate_keys(self):
        result = scan_composite_expressions()
        if result["best_matches"]:
            c = result["best_matches"][0]
            for key in ("name", "value", "residual_ppm"):
                self.assertIn(key, c)

    def test_assessment_nonempty(self):
        result = scan_composite_expressions()
        self.assertIsInstance(result["assessment"], str)
        self.assertGreater(len(result["assessment"]), 0)


class TestScanKClosedForms(unittest.TestCase):
    """Tests for scan_k_closed_forms."""

    def test_returns_dict(self):
        result = scan_k_closed_forms()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = scan_k_closed_forms()
        for key in ("candidates", "best_matches", "assessment"):
            self.assertIn(key, result)

    def test_sqrt_phi_in_candidates(self):
        """sqrt(phi) should be one of the candidates."""
        result = scan_k_closed_forms()
        names = [c["name"] for c in result["candidates"]]
        self.assertIn("sqrt(phi)", names)

    def test_sqrt_phi_residual_matches_mu_structure(self):
        """sqrt(phi) residual should be ~7800 ppm (since k_exact != sqrt(phi))."""
        result = scan_k_closed_forms()
        sqrt_phi_entry = [c for c in result["candidates"] if c["name"] == "sqrt(phi)"][0]
        self.assertGreater(sqrt_phi_entry["residual_ppm"], 5000)
        self.assertLess(sqrt_phi_entry["residual_ppm"], 12000)

    def test_best_matches_sorted(self):
        result = scan_k_closed_forms()
        ppms = [c["residual_ppm"] for c in result["best_matches"]]
        self.assertEqual(ppms, sorted(ppms))

    def test_best_matches_max_10(self):
        result = scan_k_closed_forms()
        self.assertLessEqual(len(result["best_matches"]), 10)

    def test_assessment_nonempty(self):
        result = scan_k_closed_forms()
        self.assertIsInstance(result["assessment"], str)
        self.assertGreater(len(result["assessment"]), 0)


class TestAnalyzeDeltaStructure(unittest.TestCase):
    """Tests for analyze_delta_structure."""

    def test_returns_dict(self):
        result = analyze_delta_structure()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = analyze_delta_structure()
        for key in ("delta", "delta_over_alpha", "delta_over_alpha_sq",
                     "delta_times_mu", "delta_over_schwinger",
                     "delta_as_alpha_power", "interpretation", "assessment"):
            self.assertIn(key, result)

    def test_delta_as_alpha_power_near_1(self):
        """delta should be approximately alpha^0.9-1.1, i.e., O(alpha)."""
        result = analyze_delta_structure()
        self.assertGreater(result["delta_as_alpha_power"], 0.8)
        self.assertLess(result["delta_as_alpha_power"], 1.2)

    def test_delta_times_mu_approx(self):
        """delta * mu should be approximately 18."""
        result = analyze_delta_structure()
        self.assertAlmostEqual(result["delta_times_mu"], 18.0, delta=1.0)

    def test_interpretation_nonempty(self):
        result = analyze_delta_structure()
        self.assertIsInstance(result["interpretation"], str)
        self.assertGreater(len(result["interpretation"]), 0)

    def test_assessment_nonempty(self):
        result = analyze_delta_structure()
        self.assertIsInstance(result["assessment"], str)
        self.assertGreater(len(result["assessment"]), 0)


class TestSummarizeResidualMapping(unittest.TestCase):
    """Tests for summarize_residual_mapping."""

    def test_returns_dict(self):
        result = summarize_residual_mapping()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = summarize_residual_mapping()
        for key in ("delta_info", "g2_scan", "zeta_scan", "composite_scan",
                     "k_closed_forms", "delta_structure", "key_finding",
                     "honest_assessment"):
            self.assertIn(key, result)

    def test_key_finding_nonempty(self):
        result = summarize_residual_mapping()
        self.assertIsInstance(result["key_finding"], str)
        self.assertGreater(len(result["key_finding"]), 0)

    def test_honest_assessment_nonempty(self):
        result = summarize_residual_mapping()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 0)

    def test_sub_results_are_dicts(self):
        result = summarize_residual_mapping()
        for key in ("delta_info", "g2_scan", "zeta_scan", "composite_scan",
                     "k_closed_forms", "delta_structure"):
            self.assertIsInstance(result[key], dict)


if __name__ == "__main__":
    unittest.main()
