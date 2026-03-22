"""Tests for phi_vev_golden_ratio.py -- dilaton vev vs golden ratio relationship check."""
import unittest
import math



class TestModuleConstants(unittest.TestCase):
    """Test module-level constants."""

    def test_phi_vev_value(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import PHI_VEV
        self.assertAlmostEqual(PHI_VEV, -0.5973, delta=0.01)

    def test_phi_golden_value(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import PHI_GOLDEN
        self.assertAlmostEqual(PHI_GOLDEN, (1 + math.sqrt(5)) / 2, places=10)

    def test_c0_value(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import C_0, PHI_GOLDEN
        self.assertAlmostEqual(C_0, PHI_GOLDEN ** 2 / 2, places=10)

    def test_phi_golden_satisfies_quadratic(self):
        """phi_golden is algebraic: x^2 - x - 1 = 0."""
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import PHI_GOLDEN
        self.assertAlmostEqual(PHI_GOLDEN ** 2 - PHI_GOLDEN - 1.0, 0.0, places=10)


class TestExponentialPowerRelations(unittest.TestCase):
    """Test check_exponential_power_relations function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_exponential_power_relations
        result = check_exponential_power_relations()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_exponential_power_relations
        result = check_exponential_power_relations()
        for key in ["matches", "closest", "any_clean_match", "n_matches_within_1pct"]:
            self.assertIn(key, result)

    def test_no_clean_match(self):
        """No exp(k*phi_vev) = phi_golden^n match within 100 ppm."""
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_exponential_power_relations
        result = check_exponential_power_relations()
        self.assertFalse(result["any_clean_match"])

    def test_matches_sorted_by_ppm(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_exponential_power_relations
        result = check_exponential_power_relations()
        for i in range(len(result["matches"]) - 1):
            self.assertLessEqual(result["matches"][i]["ppm_off"],
                                 result["matches"][i + 1]["ppm_off"])

    def test_closest_has_k_and_n(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_exponential_power_relations
        result = check_exponential_power_relations()
        if result["closest"] is not None:
            self.assertIn("k", result["closest"])
            self.assertIn("n", result["closest"])
            self.assertIn("ppm_off", result["closest"])

    def test_ppm_values_positive(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_exponential_power_relations
        result = check_exponential_power_relations()
        for m in result["matches"]:
            self.assertGreaterEqual(m["ppm_off"], 0)

    def test_has_description(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_exponential_power_relations
        result = check_exponential_power_relations()
        self.assertIn("description", result)

    def test_exp_val_computation(self):
        """Verify that exp_val = exp(k * phi_vev)."""
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_exponential_power_relations, PHI_VEV
        result = check_exponential_power_relations()
        if result["matches"]:
            m = result["matches"][0]
            expected = math.exp(m["k"] * PHI_VEV)
            self.assertAlmostEqual(m["exp_val"], expected, places=6)


class TestAlgebraicCombinations(unittest.TestCase):
    """Test check_algebraic_combinations function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_algebraic_combinations
        result = check_algebraic_combinations()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_algebraic_combinations
        result = check_algebraic_combinations()
        for key in ["checks", "any_match", "closest"]:
            self.assertIn(key, result)

    def test_checks_list_nonempty(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_algebraic_combinations
        result = check_algebraic_combinations()
        self.assertGreater(len(result["checks"]), 5)

    def test_each_check_has_formula(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_algebraic_combinations
        result = check_algebraic_combinations()
        for c in result["checks"]:
            self.assertIn("formula", c)
            self.assertIn("value", c)
            self.assertIn("target", c)
            self.assertIn("ppm_off", c)

    def test_e_phi_plus_1_not_phi_golden(self):
        """e^(phi_vev) + 1 is NOT the golden ratio."""
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_algebraic_combinations, PHI_VEV, PHI_GOLDEN
        val = math.exp(PHI_VEV) + 1
        ppm = abs(val / PHI_GOLDEN - 1) * 1e6
        self.assertGreater(ppm, 1000)  # Not close

    def test_phi_vev_squared_not_1_over_phi_golden(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import PHI_VEV, PHI_GOLDEN
        val = PHI_VEV ** 2
        target = 1 / PHI_GOLDEN
        ppm = abs(val / target - 1) * 1e6
        self.assertGreater(ppm, 1000)


class TestBridgeVsGaugeMatching(unittest.TestCase):
    """Test check_bridge_vs_gauge_matching function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_bridge_vs_gauge_matching
        result = check_bridge_vs_gauge_matching()
        self.assertIsInstance(result, dict)

    def test_c0_value(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_bridge_vs_gauge_matching
        result = check_bridge_vs_gauge_matching()
        self.assertAlmostEqual(result["C_0"], 1.309, delta=0.01)

    def test_exp_4phi_value(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_bridge_vs_gauge_matching, ALPHA
        result = check_bridge_vs_gauge_matching()
        self.assertAlmostEqual(result["exp_4phi"], 4 * math.pi * ALPHA, delta=0.01)

    def test_ratio_not_exactly_14(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_bridge_vs_gauge_matching
        result = check_bridge_vs_gauge_matching()
        self.assertGreater(result["ratio_close_to_14"], 0.5)  # More than 0.5% off

    def test_ratio_not_exactly_4pi(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_bridge_vs_gauge_matching
        result = check_bridge_vs_gauge_matching()
        self.assertGreater(result["ratio_close_to_4pi"], 0.5)

    def test_has_description(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import check_bridge_vs_gauge_matching
        result = check_bridge_vs_gauge_matching()
        self.assertIn("description", result)


class TestMathematicalIncompatibility(unittest.TestCase):
    """Test mathematical_incompatibility function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import mathematical_incompatibility
        result = mathematical_incompatibility()
        self.assertIsInstance(result, dict)

    def test_phi_vev_transcendental(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import mathematical_incompatibility
        result = mathematical_incompatibility()
        self.assertEqual(result["phi_vev_type"], "transcendental")

    def test_phi_golden_algebraic(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import mathematical_incompatibility
        result = mathematical_incompatibility()
        self.assertEqual(result["phi_golden_type"], "algebraic")

    def test_phi_golden_satisfies_quadratic(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import mathematical_incompatibility
        result = mathematical_incompatibility()
        self.assertEqual(result["phi_golden_satisfies"], "x^2 - x - 1 = 0")

    def test_lindemann_weierstrass_present(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import mathematical_incompatibility
        result = mathematical_incompatibility()
        self.assertIn("lindemann_weierstrass", result)
        self.assertIn("Lindemann-Weierstrass", result["lindemann_weierstrass"])

    def test_phi_vev_involves_ln_pi_alpha(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import mathematical_incompatibility
        result = mathematical_incompatibility()
        involves = result["phi_vev_involves"]
        self.assertIn("ln", involves)
        self.assertIn("pi", involves)

    def test_connection_would_require_present(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import mathematical_incompatibility
        result = mathematical_incompatibility()
        self.assertIn("connection_would_require", result)
        self.assertIsInstance(result["connection_would_require"], str)


class TestSummarizePhiVevGoldenRatio(unittest.TestCase):
    """Test summarize_phi_vev_golden_ratio function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import summarize_phi_vev_golden_ratio
        result = summarize_phi_vev_golden_ratio()
        self.assertIsInstance(result, dict)

    def test_no_relationship_found(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import summarize_phi_vev_golden_ratio
        result = summarize_phi_vev_golden_ratio()
        self.assertFalse(result["relationship_found"])

    def test_has_all_sub_results(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import summarize_phi_vev_golden_ratio
        result = summarize_phi_vev_golden_ratio()
        for key in ["exponential_power_search", "algebraic_combinations",
                     "bridge_vs_gauge", "mathematical_argument"]:
            self.assertIn(key, result)

    def test_has_key_messages(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import summarize_phi_vev_golden_ratio
        result = summarize_phi_vev_golden_ratio()
        self.assertIn("key_messages", result)
        self.assertGreaterEqual(len(result["key_messages"]), 5)

    def test_has_honest_assessment(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import summarize_phi_vev_golden_ratio
        result = summarize_phi_vev_golden_ratio()
        self.assertIn("honest_assessment", result)
        self.assertIn("No clean relationship", result["honest_assessment"])

    def test_constants_in_summary(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import summarize_phi_vev_golden_ratio, PHI_VEV, PHI_GOLDEN
        result = summarize_phi_vev_golden_ratio()
        self.assertAlmostEqual(result["phi_vev"], PHI_VEV, places=10)
        self.assertAlmostEqual(result["phi_golden"], PHI_GOLDEN, places=10)

    def test_key_messages_mention_transcendental(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import summarize_phi_vev_golden_ratio
        result = summarize_phi_vev_golden_ratio()
        messages_text = " ".join(result["key_messages"])
        self.assertIn("transcendental", messages_text)

    def test_key_messages_mention_algebraic(self):
        from alpha_ladder_core.gauge_phi_vev_golden_ratio import summarize_phi_vev_golden_ratio
        result = summarize_phi_vev_golden_ratio()
        messages_text = " ".join(result["key_messages"])
        self.assertIn("algebraic", messages_text)


if __name__ == "__main__":
    unittest.main()
