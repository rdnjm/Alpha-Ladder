"""Tests for the hierarchy derivation module."""

import unittest
import math

from alpha_ladder_core.hierarchy_derivation import (
    compute_formula_basics,
    analyze_kk_volume_suppression,
    analyze_metric_component_counting,
    analyze_power_law_running,
    analyze_volume_mass_relation,
    analyze_induced_gravity,
    scan_dimension_pairs,
    analyze_residual,
    summarize_hierarchy_derivation,
    analyze_swampland_distance,
    analyze_emergence_tower,
    analyze_one_loop_matching,
)


class TestFormulaBasics(unittest.TestCase):
    """Tests for compute_formula_basics."""

    def test_returns_dict(self):
        result = compute_formula_basics()
        self.assertIsInstance(result, dict)

    def test_exponent_is_24(self):
        result = compute_formula_basics()
        self.assertEqual(result["exponent"], 24)

    def test_residual_about_688_ppm(self):
        result = compute_formula_basics()
        self.assertGreater(abs(result["residual_ppm"]), 600)
        self.assertLess(abs(result["residual_ppm"]), 800)

    def test_ratio_near_one(self):
        result = compute_formula_basics()
        self.assertLess(abs(result["ratio"] - 1), 0.001)

    def test_G_predicted_in_range(self):
        result = compute_formula_basics()
        self.assertGreater(result["G_predicted"], 6.6e-11)
        self.assertLess(result["G_predicted"], 6.8e-11)

    def test_G_measured_in_range(self):
        result = compute_formula_basics()
        self.assertGreater(result["G_measured"], 6.6e-11)
        self.assertLess(result["G_measured"], 6.8e-11)

    def test_mu_about_1836(self):
        result = compute_formula_basics()
        self.assertLess(abs(result["mu"] - 1836.15), 1.0)


class TestKKVolumeSuppression(unittest.TestCase):
    """Tests for analyze_kk_volume_suppression."""

    def test_returns_dict(self):
        result = analyze_kk_volume_suppression()
        self.assertIsInstance(result, dict)

    def test_works_is_false(self):
        result = analyze_kk_volume_suppression()
        self.assertFalse(result["works"])

    def test_M_D_finite(self):
        result = analyze_kk_volume_suppression()
        self.assertIsNotNone(result["M_D"])
        if isinstance(result["M_D"], (int, float)):
            self.assertTrue(math.isfinite(result["M_D"]))

    def test_gap_description_nonempty(self):
        result = analyze_kk_volume_suppression()
        self.assertGreater(len(result["gap_description"]), 0)

    def test_n_powers_less_than_dD(self):
        result = analyze_kk_volume_suppression()
        self.assertLess(result["n_powers_of_radius"], result["required_powers"])


class TestMetricComponentCounting(unittest.TestCase):
    """Tests for analyze_metric_component_counting."""

    def test_returns_dict(self):
        result = analyze_metric_component_counting()
        self.assertIsInstance(result, dict)

    def test_symmetric_components_21(self):
        result = analyze_metric_component_counting()
        self.assertEqual(result["symmetric_components"], 21)

    def test_dD_product_24(self):
        result = analyze_metric_component_counting()
        self.assertEqual(result["dD_product"], 24)

    def test_kk_decomposition_sums_to_21(self):
        result = analyze_metric_component_counting()
        kk = result["kk_decomposition"]
        total = kk["graviton"] + kk["vectors"] + kk["scalars"]
        self.assertEqual(total, 21)

    def test_dD_not_equal_symmetric(self):
        result = analyze_metric_component_counting()
        self.assertNotEqual(result["dD_product"], result["symmetric_components"])

    def test_mismatch_is_3(self):
        result = analyze_metric_component_counting()
        self.assertEqual(result["mismatch"], 3)


class TestPowerLawRunning(unittest.TestCase):
    """Tests for analyze_power_law_running."""

    def test_returns_dict(self):
        result = analyze_power_law_running()
        self.assertIsInstance(result, dict)

    def test_effective_exponent_finite(self):
        result = analyze_power_law_running()
        self.assertTrue(math.isfinite(result["effective_exponent"]))

    def test_effective_exponent_positive(self):
        result = analyze_power_law_running()
        self.assertGreater(result["effective_exponent"], 0)

    def test_description_nonempty(self):
        result = analyze_power_law_running()
        has_content = (
            len(result.get("description", "")) > 0
            or len(result.get("honest", "")) > 0
        )
        self.assertTrue(has_content)


class TestVolumeMassRelation(unittest.TestCase):
    """Tests for analyze_volume_mass_relation."""

    def test_returns_dict(self):
        result = analyze_volume_mass_relation()
        self.assertIsInstance(result, dict)

    def test_mu_alpha_exponent_about_1_53(self):
        result = analyze_volume_mass_relation()
        self.assertLess(abs(result["mu_alpha_exponent"] - 1.53), 0.1)

    def test_chain_works_false(self):
        result = analyze_volume_mass_relation()
        self.assertFalse(result["chain_works"])

    def test_chain_description_nonempty(self):
        result = analyze_volume_mass_relation()
        self.assertGreater(len(result["chain_description"]), 0)


class TestInducedGravity(unittest.TestCase):
    """Tests for analyze_induced_gravity."""

    def test_returns_dict(self):
        result = analyze_induced_gravity()
        self.assertIsInstance(result, dict)

    def test_N_species_24(self):
        result = analyze_induced_gravity()
        self.assertEqual(result["N_species"], 24)

    def test_Lambda_UV_finite(self):
        result = analyze_induced_gravity()
        self.assertTrue(math.isfinite(result["Lambda_UV"]))

    def test_Lambda_UV_positive(self):
        result = analyze_induced_gravity()
        self.assertGreater(result["Lambda_UV"], 0)

    def test_works_is_false(self):
        result = analyze_induced_gravity()
        self.assertFalse(result["works"])


class TestDimensionScan(unittest.TestCase):
    """Tests for scan_dimension_pairs."""

    def test_returns_dict(self):
        result = scan_dimension_pairs()
        self.assertIsInstance(result, dict)

    def test_scan_has_25_entries(self):
        result = scan_dimension_pairs()
        self.assertEqual(len(result["scan_results"]), 25)

    def test_best_pair_is_4_2(self):
        result = scan_dimension_pairs()
        self.assertEqual(result["best_pair"], (4, 2))

    def test_best_ppm_about_688(self):
        result = scan_dimension_pairs()
        self.assertLess(abs(result["best_ppm"] - 688), 100)

    def test_unique_sub_1000(self):
        result = scan_dimension_pairs()
        self.assertTrue(result["unique_sub_1000"])

    def test_runner_up_much_larger(self):
        result = scan_dimension_pairs()
        self.assertGreater(result["runner_up_ppm"], 1000)

    def test_all_entries_have_ppm(self):
        result = scan_dimension_pairs()
        self.assertTrue(all("ppm" in e for e in result["scan_results"]))


class TestResidualAnalysis(unittest.TestCase):
    """Tests for analyze_residual."""

    def test_returns_dict(self):
        result = analyze_residual()
        self.assertIsInstance(result, dict)

    def test_raw_ppm_about_neg_688(self):
        result = analyze_residual()
        self.assertGreater(abs(result["raw_ppm"]), 600)
        self.assertLess(abs(result["raw_ppm"]), 800)

    def test_corrections_nonempty(self):
        result = analyze_residual()
        self.assertGreater(len(result["corrections"]), 0)

    def test_best_correction_improves(self):
        result = analyze_residual()
        self.assertIn("best_correction", result)
        self.assertIn("corrected_ppm", result["best_correction"])
        self.assertLess(
            abs(result["best_correction"]["corrected_ppm"]),
            abs(result["raw_ppm"]),
        )

    def test_qed_candidate_present(self):
        result = analyze_residual()
        self.assertTrue(any("QED" in c["name"] for c in result["corrections"]))


class TestSummarize(unittest.TestCase):
    """Tests for summarize."""

    def test_returns_dict(self):
        result = summarize_hierarchy_derivation()
        self.assertIsInstance(result, dict)

    def test_has_formula_basics(self):
        result = summarize_hierarchy_derivation()
        self.assertIn("formula_basics", result)

    def test_honest_assessment_is_string(self):
        result = summarize_hierarchy_derivation()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 0)

    def test_angles_summary_has_10(self):
        result = summarize_hierarchy_derivation()
        self.assertEqual(len(result["angles_summary"]), 10)

    def test_all_angle_keys_present(self):
        result = summarize_hierarchy_derivation()
        for i in range(1, 11):
            self.assertIn(f"angle_{i}", result)


class TestSwamplandDistance(unittest.TestCase):
    """Tests for analyze_swampland_distance."""

    def test_returns_dict(self):
        result = analyze_swampland_distance()
        self.assertIsInstance(result, dict)

    def test_lambda_EH_value(self):
        result = analyze_swampland_distance()
        # lambda_EH = 1/(2*sqrt(2)) for K_EH = -4, d=4, n=2
        # lambda_EH = 1/sqrt(2*|K_EH|) = 1/sqrt(8) = 1/(2*sqrt(2))
        self.assertAlmostEqual(result["lambda_EH"], 1 / (2 * math.sqrt(2)), places=10)

    def test_lambda_bound(self):
        result = analyze_swampland_distance()
        # lambda_bound = 1/sqrt(d-2) = 1/sqrt(2)
        self.assertAlmostEqual(result["lambda_bound"], 1 / math.sqrt(2), places=10)

    def test_sdc_violated_EH(self):
        result = analyze_swampland_distance()
        self.assertFalse(result["sdc_satisfied_EH"])

    def test_sdc_satisfied_GB(self):
        result = analyze_swampland_distance()
        self.assertTrue(result["sdc_satisfied_GB"])

    def test_lambda_GB_exceeds_bound(self):
        result = analyze_swampland_distance()
        self.assertGreater(result["lambda_GB"], result["lambda_bound"])

    def test_K_EH_is_minus_4(self):
        result = analyze_swampland_distance()
        self.assertAlmostEqual(result["K_EH"], -4.0, places=10)

    def test_K_GB_magnitude(self):
        result = analyze_swampland_distance()
        # abs_K_GB = omega_target * 2 * n^2 * beta^2
        # For d=4, n=2: beta = -1, so beta^2 = 1
        # abs_K_GB = omega_target * 2 * 4 * 1 = 8 * omega_target
        # omega_target ~ 0.118, so abs_K_GB ~ 0.944
        self.assertAlmostEqual(result["abs_K_GB"], 0.944, delta=0.01)

    def test_omega_saturation(self):
        result = analyze_swampland_distance()
        # omega_sat = abs_K_sat / (2 * n^2 * beta^2)
        # abs_K_sat = (d-2)/2 = 1
        # omega_sat = 1 / (2 * 4 * 1) = 0.125
        self.assertAlmostEqual(result["omega_saturation"], 0.125, places=10)

    def test_near_saturation(self):
        result = analyze_swampland_distance()
        # omega_target ~ 0.118, omega_sat = 0.125, diff ~ 0.007 < 0.01
        self.assertTrue(result["near_saturation"])

    def test_margin_positive(self):
        result = analyze_swampland_distance()
        self.assertGreater(result["margin"], 0)

    def test_assessment_nonempty(self):
        result = analyze_swampland_distance()
        self.assertGreater(len(result["assessment"]), 10)


class TestEmergenceTower(unittest.TestCase):
    """Tests for analyze_emergence_tower."""

    def test_returns_dict(self):
        result = analyze_emergence_tower()
        self.assertIsInstance(result, dict)

    def test_has_canonical_R_results(self):
        result = analyze_emergence_tower()
        self.assertIn("canonical_R_results", result)
        self.assertGreaterEqual(len(result["canonical_R_results"]), 3)

    def test_works_is_false(self):
        result = analyze_emergence_tower()
        self.assertFalse(result["works"])

    def test_honest_assessment_nonempty(self):
        result = analyze_emergence_tower()
        self.assertGreater(len(result["honest_assessment"]), 10)

    def test_scaling_present(self):
        result = analyze_emergence_tower()
        self.assertIn("scaling_large_L", result)

    def test_canonical_results_have_keys(self):
        result = analyze_emergence_tower()
        for r in result["canonical_R_results"]:
            self.assertIn("name", r)
            self.assertIn("R_m", r)

    def test_N_pol_default(self):
        result = analyze_emergence_tower()
        self.assertEqual(result["N_pol_default"], 21)

    def test_R_needed_finite(self):
        result = analyze_emergence_tower()
        if result.get("R_needed_for_match_m") is not None:
            self.assertGreater(result["R_needed_for_match_m"], 0)


class TestOneLoopMatching(unittest.TestCase):
    """Tests for analyze_one_loop_matching."""

    def test_returns_dict(self):
        result = analyze_one_loop_matching()
        self.assertIsInstance(result, dict)

    def test_uv_divergence_degree(self):
        result = analyze_one_loop_matching()
        # D-2 = 6-2 = 4
        self.assertEqual(result["uv_divergence_degree"], 4)

    def test_specific_cases_present(self):
        result = analyze_one_loop_matching()
        self.assertIn("specific_cases", result)
        self.assertGreaterEqual(len(result["specific_cases"]), 2)

    def test_works_is_false(self):
        result = analyze_one_loop_matching()
        self.assertFalse(result["works"])

    def test_mu_exponent_problem(self):
        result = analyze_one_loop_matching()
        self.assertIn("mu_exponent_problem", result)
        self.assertIn("mu", result["mu_exponent_problem"].lower())

    def test_honest_assessment_nonempty(self):
        result = analyze_one_loop_matching()
        self.assertGreater(len(result["honest_assessment"]), 10)

    def test_analytic_scaling(self):
        result = analyze_one_loop_matching()
        self.assertIn("analytic_scaling", result)

    def test_alpha_exponent_scan_present(self):
        result = analyze_one_loop_matching()
        self.assertIn("alpha_exponent_scan", result)


if __name__ == "__main__":
    unittest.main()
