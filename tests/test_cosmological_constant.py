"""Tests for the cosmological constant module."""

import unittest
import math


class TestVacuumEnergy(unittest.TestCase):
    """Tests for extract_vacuum_energy."""

    def test_returns_dict(self):
        from alpha_ladder_core.cosmological_constant import extract_vacuum_energy
        result = extract_vacuum_energy()
        self.assertIsInstance(result, dict)

    def test_v_min_finite(self):
        from alpha_ladder_core.cosmological_constant import extract_vacuum_energy
        result = extract_vacuum_energy(N=1)
        self.assertIsNotNone(result["V_min_planck"])
        self.assertTrue(math.isfinite(result["V_min_planck"]))

    def test_v_min_nonzero(self):
        from alpha_ladder_core.cosmological_constant import extract_vacuum_energy
        result = extract_vacuum_energy(N=1)
        self.assertNotEqual(result["V_min_planck"], 0.0)

    def test_minimum_exists(self):
        from alpha_ladder_core.cosmological_constant import extract_vacuum_energy
        result = extract_vacuum_energy(N=1)
        self.assertTrue(result["minimum_exists"])

    def test_v_min_eV4_computed(self):
        from alpha_ladder_core.cosmological_constant import extract_vacuum_energy
        result = extract_vacuum_energy(N=1)
        self.assertIsNotNone(result["V_min_eV4"])

    def test_sign_is_string(self):
        from alpha_ladder_core.cosmological_constant import extract_vacuum_energy
        result = extract_vacuum_energy(N=1)
        self.assertIn(result["sign"], ["positive", "negative", "zero"])

    def test_log10_magnitude_computed(self):
        from alpha_ladder_core.cosmological_constant import extract_vacuum_energy
        result = extract_vacuum_energy(N=1)
        self.assertIsNotNone(result["log10_magnitude"])

    def test_consistent_with_flux_stabilization(self):
        from alpha_ladder_core.cosmological_constant import extract_vacuum_energy
        from alpha_ladder_core.flux_stabilization import find_flux_minimum
        vacuum = extract_vacuum_energy(N=1, a_0=1.0)
        flux = find_flux_minimum(N=1, a_0=1.0)
        self.assertAlmostEqual(vacuum["V_min_planck"], flux["V_at_minimum"], places=10)


class TestCompareObservation(unittest.TestCase):
    """Tests for compare_with_observation."""

    def test_returns_dict(self):
        from alpha_ladder_core.cosmological_constant import compare_with_observation
        result = compare_with_observation()
        self.assertIsInstance(result, dict)

    def test_lambda_obs_correct_order(self):
        from alpha_ladder_core.cosmological_constant import compare_with_observation
        result = compare_with_observation()
        log_lambda = math.log10(result["Lambda_obs_planck"])
        self.assertAlmostEqual(log_lambda, -122, delta=1)

    def test_ratio_huge(self):
        from alpha_ladder_core.cosmological_constant import compare_with_observation
        result = compare_with_observation()
        self.assertIsNotNone(result["ratio"])
        self.assertGreater(result["ratio"], 1e100)

    def test_log10_ratio_about_122(self):
        from alpha_ladder_core.cosmological_constant import compare_with_observation
        result = compare_with_observation()
        self.assertIsNotNone(result["log10_ratio"])
        self.assertGreater(result["log10_ratio"], 100)

    def test_discrepancy_orders(self):
        from alpha_ladder_core.cosmological_constant import compare_with_observation
        result = compare_with_observation()
        self.assertIsNotNone(result["discrepancy_orders"])
        self.assertGreater(result["discrepancy_orders"], 100)

    def test_is_fine_tuning(self):
        from alpha_ladder_core.cosmological_constant import compare_with_observation
        result = compare_with_observation()
        self.assertTrue(result["is_fine_tuning"])

    def test_honest_assessment_nonempty(self):
        from alpha_ladder_core.cosmological_constant import compare_with_observation
        result = compare_with_observation()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 50)

    def test_lambda_obs_eV4_tiny(self):
        from alpha_ladder_core.cosmological_constant import compare_with_observation
        result = compare_with_observation()
        # Lambda_obs in eV^4 should be extremely small
        self.assertGreater(result["Lambda_obs_eV4"], 0)


class TestCCMechanisms(unittest.TestCase):
    """Tests for analyze_cc_mechanisms."""

    def test_returns_dict(self):
        from alpha_ladder_core.cosmological_constant import analyze_cc_mechanisms
        result = analyze_cc_mechanisms()
        self.assertIsInstance(result, dict)

    def test_at_least_three_mechanisms(self):
        from alpha_ladder_core.cosmological_constant import analyze_cc_mechanisms
        result = analyze_cc_mechanisms()
        self.assertGreaterEqual(result["n_mechanisms"], 3)

    def test_no_resolution(self):
        from alpha_ladder_core.cosmological_constant import analyze_cc_mechanisms
        result = analyze_cc_mechanisms()
        self.assertFalse(result["any_resolution"])

    def test_universal_problem(self):
        from alpha_ladder_core.cosmological_constant import analyze_cc_mechanisms
        result = analyze_cc_mechanisms()
        self.assertTrue(result["universal_problem"])

    def test_mechanisms_have_keys(self):
        from alpha_ladder_core.cosmological_constant import analyze_cc_mechanisms
        result = analyze_cc_mechanisms()
        for m in result["mechanisms"]:
            self.assertIn("name", m)
            self.assertIn("applicable_here", m)
            self.assertIn("resolves_problem", m)

    def test_none_resolve_problem(self):
        from alpha_ladder_core.cosmological_constant import analyze_cc_mechanisms
        result = analyze_cc_mechanisms()
        for m in result["mechanisms"]:
            self.assertFalse(m["resolves_problem"])

    def test_none_applicable(self):
        from alpha_ladder_core.cosmological_constant import analyze_cc_mechanisms
        result = analyze_cc_mechanisms()
        for m in result["mechanisms"]:
            self.assertFalse(m["applicable_here"])

    def test_honest_assessment_nonempty(self):
        from alpha_ladder_core.cosmological_constant import analyze_cc_mechanisms
        result = analyze_cc_mechanisms()
        self.assertGreater(len(result["honest_assessment"]), 20)


class TestCCScan(unittest.TestCase):
    """Tests for compute_cc_scan."""

    def test_returns_dict(self):
        from alpha_ladder_core.cosmological_constant import compute_cc_scan
        result = compute_cc_scan(N_max=3)
        self.assertIsInstance(result, dict)

    def test_scan_results_length(self):
        from alpha_ladder_core.cosmological_constant import compute_cc_scan
        result = compute_cc_scan(N_max=5)
        self.assertEqual(len(result["scan_results"]), 5)

    def test_all_O1_planck(self):
        from alpha_ladder_core.cosmological_constant import compute_cc_scan
        result = compute_cc_scan(N_max=5)
        self.assertTrue(result["all_O1_planck"])

    def test_none_close_to_observation(self):
        from alpha_ladder_core.cosmological_constant import compute_cc_scan
        result = compute_cc_scan(N_max=5)
        self.assertFalse(result["any_close_to_observation"])

    def test_log10_ratios_all_large(self):
        from alpha_ladder_core.cosmological_constant import compute_cc_scan
        result = compute_cc_scan(N_max=5)
        for entry in result["scan_results"]:
            if entry["log10_ratio"] is not None:
                self.assertGreater(entry["log10_ratio"], 100)

    def test_min_log10_ratio_large(self):
        from alpha_ladder_core.cosmological_constant import compute_cc_scan
        result = compute_cc_scan(N_max=5)
        self.assertGreater(result["min_log10_ratio"], 100)

    def test_max_log10_ratio_large(self):
        from alpha_ladder_core.cosmological_constant import compute_cc_scan
        result = compute_cc_scan(N_max=5)
        self.assertGreater(result["max_log10_ratio"], 100)

    def test_honest_assessment_nonempty(self):
        from alpha_ladder_core.cosmological_constant import compute_cc_scan
        result = compute_cc_scan(N_max=3)
        self.assertGreater(len(result["honest_assessment"]), 20)


class TestNoGo(unittest.TestCase):
    """Tests for assess_no_go_theorem."""

    def test_returns_dict(self):
        from alpha_ladder_core.cosmological_constant import assess_no_go_theorem
        result = assess_no_go_theorem()
        self.assertIsInstance(result, dict)

    def test_applies_to_framework(self):
        from alpha_ladder_core.cosmological_constant import assess_no_go_theorem
        result = assess_no_go_theorem()
        self.assertTrue(result["applies_to_framework"])

    def test_all_conditions_met(self):
        from alpha_ladder_core.cosmological_constant import assess_no_go_theorem
        result = assess_no_go_theorem()
        self.assertTrue(result["all_conditions_met"])

    def test_smooth_potential(self):
        from alpha_ladder_core.cosmological_constant import assess_no_go_theorem
        result = assess_no_go_theorem()
        self.assertTrue(result["smooth_potential"])

    def test_lorentz_invariant(self):
        from alpha_ladder_core.cosmological_constant import assess_no_go_theorem
        result = assess_no_go_theorem()
        self.assertTrue(result["lorentz_invariant"])

    def test_standard_kinetic(self):
        from alpha_ladder_core.cosmological_constant import assess_no_go_theorem
        result = assess_no_go_theorem()
        self.assertTrue(result["standard_kinetic"])

    def test_conditions_list(self):
        from alpha_ladder_core.cosmological_constant import assess_no_go_theorem
        result = assess_no_go_theorem()
        self.assertIsInstance(result["conditions"], list)
        self.assertEqual(len(result["conditions"]), 3)

    def test_reference_weinberg(self):
        from alpha_ladder_core.cosmological_constant import assess_no_go_theorem
        result = assess_no_go_theorem()
        self.assertIn("Weinberg", result["reference"])


if __name__ == "__main__":
    unittest.main()
