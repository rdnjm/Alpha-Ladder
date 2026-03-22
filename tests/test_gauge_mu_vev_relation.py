"""Tests for mu_vev_relation.py -- proton-to-electron mass ratio vs dilaton vev."""
import unittest
import math



class TestModuleConstants(unittest.TestCase):
    """Test module-level constants."""

    def test_phi_vev_value(self):
        from alpha_ladder_core.gauge_mu_vev_relation import PHI_VEV
        self.assertAlmostEqual(PHI_VEV, -0.5973, delta=0.01)

    def test_mu_value(self):
        from alpha_ladder_core.gauge_mu_vev_relation import MU
        self.assertAlmostEqual(MU, 1836.15267344, places=5)

    def test_phi_golden_value(self):
        from alpha_ladder_core.gauge_mu_vev_relation import PHI_GOLDEN
        self.assertAlmostEqual(PHI_GOLDEN, (1 + math.sqrt(5)) / 2, places=10)

    def test_alpha_value(self):
        from alpha_ladder_core.gauge_mu_vev_relation import ALPHA
        self.assertAlmostEqual(ALPHA, 0.0072973525693, places=10)


class TestExponentialRelations(unittest.TestCase):
    """Test check_exponential_relations function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_exponential_relations
        result = check_exponential_relations()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_exponential_relations
        result = check_exponential_relations()
        for key in ["relations", "best_k", "best_ratio", "exact_k", "clean_match"]:
            self.assertIn(key, result)

    def test_no_clean_match(self):
        """phi_vev and mu are numerically independent -- no clean exp(k*phi_vev) = mu."""
        from alpha_ladder_core.gauge_mu_vev_relation import check_exponential_relations
        result = check_exponential_relations()
        self.assertFalse(result["clean_match"])

    def test_exact_k_not_simple(self):
        """The exact k for exp(k*phi_vev) = mu is not a simple number."""
        from alpha_ladder_core.gauge_mu_vev_relation import check_exponential_relations
        result = check_exponential_relations()
        k = result["exact_k"]
        # k should not be a simple integer or half-integer
        self.assertFalse(abs(k - round(k)) < 0.01)

    def test_relations_list_nonempty(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_exponential_relations
        result = check_exponential_relations()
        self.assertGreater(len(result["relations"]), 10)

    def test_each_relation_has_keys(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_exponential_relations
        result = check_exponential_relations()
        for r in result["relations"]:
            self.assertIn("k", r)
            self.assertIn("exp_val", r)
            self.assertIn("ratio_to_mu", r)

    def test_exp_val_computation(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_exponential_relations, PHI_VEV
        result = check_exponential_relations()
        for r in result["relations"][:5]:
            expected = math.exp(r["k"] * PHI_VEV)
            self.assertAlmostEqual(r["exp_val"], expected, places=8)

    def test_has_description(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_exponential_relations
        result = check_exponential_relations()
        self.assertIn("description", result)
        self.assertIsInstance(result["description"], str)


class TestCombinationRelations(unittest.TestCase):
    """Test check_combination_relations function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_combination_relations
        result = check_combination_relations()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_combination_relations
        result = check_combination_relations()
        for key in ["combinations", "any_match", "closest_name", "closest_ratio"]:
            self.assertIn(key, result)

    def test_combinations_list(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_combination_relations
        result = check_combination_relations()
        self.assertGreater(len(result["combinations"]), 10)

    def test_each_combination_has_name(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_combination_relations
        result = check_combination_relations()
        for c in result["combinations"]:
            self.assertIn("name", c)
            self.assertIn("value", c)
            self.assertIn("ratio_to_mu", c)

    def test_mu_times_alpha_identity(self):
        """mu * alpha should give a well-known value."""
        from alpha_ladder_core.gauge_mu_vev_relation import check_combination_relations, MU, ALPHA
        result = check_combination_relations()
        mu_alpha = [c for c in result["combinations"] if c["name"] == "mu * alpha"]
        self.assertEqual(len(mu_alpha), 1)
        self.assertAlmostEqual(mu_alpha[0]["value"], MU * ALPHA, places=5)


class TestBridgeConsistency(unittest.TestCase):
    """Test check_bridge_consistency function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_bridge_consistency
        result = check_bridge_consistency()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_bridge_consistency
        result = check_bridge_consistency()
        for key in ["bridge_lhs", "bridge_rhs", "match_ppm", "consistent"]:
            self.assertIn(key, result)

    def test_match_within_0_01_ppm(self):
        """Bridge consistency: 0.00 ppm match."""
        from alpha_ladder_core.gauge_mu_vev_relation import check_bridge_consistency
        result = check_bridge_consistency()
        self.assertLess(result["match_ppm"], 0.01)

    def test_consistent_flag(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_bridge_consistency
        result = check_bridge_consistency()
        self.assertTrue(result["consistent"])

    def test_lhs_positive(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_bridge_consistency
        result = check_bridge_consistency()
        self.assertGreater(result["bridge_lhs"], 0)

    def test_rhs_positive(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_bridge_consistency
        result = check_bridge_consistency()
        self.assertGreater(result["bridge_rhs"], 0)

    def test_lhs_rhs_close(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_bridge_consistency
        result = check_bridge_consistency()
        ratio = result["bridge_lhs"] / result["bridge_rhs"]
        self.assertAlmostEqual(ratio, 1.0, places=5)


class TestSpecificRelations(unittest.TestCase):
    """Test check_specific_relations function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_specific_relations
        result = check_specific_relations()
        self.assertIsInstance(result, dict)

    def test_mu_times_e_phi_not_clean(self):
        """mu * e^(phi_vev) is not a clean number."""
        from alpha_ladder_core.gauge_mu_vev_relation import check_specific_relations
        result = check_specific_relations()
        val = result["mu_times_e_phi"]
        # Should be ~1010, not a recognizable number
        self.assertGreater(val, 900)
        self.assertLess(val, 1200)

    def test_four_pi_over_alpha_close_to_mu(self):
        """4*pi/alpha = 1722 vs mu = 1836, ~6% off."""
        from alpha_ladder_core.gauge_mu_vev_relation import check_specific_relations
        result = check_specific_relations()
        ratio = result["ratio_4pi_alpha_to_mu"]
        self.assertAlmostEqual(ratio, 0.94, delta=0.1)

    def test_has_description(self):
        from alpha_ladder_core.gauge_mu_vev_relation import check_specific_relations
        result = check_specific_relations()
        self.assertIn("description", result)


class TestSummarizeMuVevAnalysis(unittest.TestCase):
    """Test summarize_mu_vev_analysis function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_mu_vev_relation import summarize_mu_vev_analysis
        result = summarize_mu_vev_analysis()
        self.assertIsInstance(result, dict)

    def test_no_relation_found(self):
        from alpha_ladder_core.gauge_mu_vev_relation import summarize_mu_vev_analysis
        result = summarize_mu_vev_analysis()
        self.assertFalse(result["any_relation_found"])

    def test_has_sub_results(self):
        from alpha_ladder_core.gauge_mu_vev_relation import summarize_mu_vev_analysis
        result = summarize_mu_vev_analysis()
        for key in ["exponential_relations", "combination_relations",
                     "bridge_consistency", "specific_relations"]:
            self.assertIn(key, result)

    def test_has_key_messages(self):
        from alpha_ladder_core.gauge_mu_vev_relation import summarize_mu_vev_analysis
        result = summarize_mu_vev_analysis()
        self.assertIn("key_messages", result)
        self.assertGreaterEqual(len(result["key_messages"]), 5)

    def test_has_honest_assessment(self):
        from alpha_ladder_core.gauge_mu_vev_relation import summarize_mu_vev_analysis
        result = summarize_mu_vev_analysis()
        self.assertIn("honest_assessment", result)
        self.assertIn("independent", result["honest_assessment"])

    def test_constants_in_summary(self):
        from alpha_ladder_core.gauge_mu_vev_relation import summarize_mu_vev_analysis, PHI_VEV, MU
        result = summarize_mu_vev_analysis()
        self.assertAlmostEqual(result["phi_vev"], PHI_VEV, places=10)
        self.assertAlmostEqual(result["mu"], MU, places=5)

    def test_bridge_consistency_in_summary(self):
        from alpha_ladder_core.gauge_mu_vev_relation import summarize_mu_vev_analysis
        result = summarize_mu_vev_analysis()
        self.assertTrue(result["bridge_consistency"]["consistent"])


if __name__ == "__main__":
    unittest.main()
