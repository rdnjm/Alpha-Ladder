"""Tests for the anomaly cancellation module."""

import unittest


class TestGravitationalAnomaly(unittest.TestCase):
    """Tests for compute_gravitational_anomaly_polynomial."""

    def test_returns_dict(self):
        from alpha_ladder_core.anomaly_cancellation import compute_gravitational_anomaly_polynomial
        result = compute_gravitational_anomaly_polynomial()
        self.assertIsInstance(result, dict)

    def test_pure_gravity_anomaly_free(self):
        from alpha_ladder_core.anomaly_cancellation import compute_gravitational_anomaly_polynomial
        result = compute_gravitational_anomaly_polynomial()
        self.assertTrue(result["pure_gravity_anomaly_free"])

    def test_gravitino_coefficient(self):
        from alpha_ladder_core.anomaly_cancellation import compute_gravitational_anomaly_polynomial
        result = compute_gravitational_anomaly_polynomial()
        self.assertAlmostEqual(result["gravitino_coefficient"], 1.0 / 5760.0, places=10)

    def test_graviton_zero_coefficient(self):
        from alpha_ladder_core.anomaly_cancellation import compute_gravitational_anomaly_polynomial
        result = compute_gravitational_anomaly_polynomial()
        self.assertEqual(result["fields"]["graviton"]["I8_coefficient"], 0.0)

    def test_graviton_not_chiral(self):
        from alpha_ladder_core.anomaly_cancellation import compute_gravitational_anomaly_polynomial
        result = compute_gravitational_anomaly_polynomial()
        self.assertFalse(result["fields"]["graviton"]["chiral"])

    def test_weyl_fermion_is_chiral(self):
        from alpha_ladder_core.anomaly_cancellation import compute_gravitational_anomaly_polynomial
        result = compute_gravitational_anomaly_polynomial()
        self.assertTrue(result["fields"]["weyl_fermion"]["chiral"])

    def test_matter_anomalous(self):
        from alpha_ladder_core.anomaly_cancellation import compute_gravitational_anomaly_polynomial
        result = compute_gravitational_anomaly_polynomial()
        self.assertTrue(result["matter_anomalous"])

    def test_anomaly_condition_present(self):
        from alpha_ladder_core.anomaly_cancellation import compute_gravitational_anomaly_polynomial
        result = compute_gravitational_anomaly_polynomial()
        self.assertIn("273", result["anomaly_condition"])


class TestGreenSchwarz(unittest.TestCase):
    """Tests for check_green_schwarz_factorization."""

    def test_e8_factorizes(self):
        from alpha_ladder_core.anomaly_cancellation import check_green_schwarz_factorization
        result = check_green_schwarz_factorization("E8_x_E8")
        self.assertTrue(result["factorizes"])

    def test_so32_factorizes(self):
        from alpha_ladder_core.anomaly_cancellation import check_green_schwarz_factorization
        result = check_green_schwarz_factorization("SO_32")
        self.assertTrue(result["factorizes"])

    def test_sm_does_not_factorize(self):
        from alpha_ladder_core.anomaly_cancellation import check_green_schwarz_factorization
        result = check_green_schwarz_factorization("SU3_x_SU2_x_U1")
        self.assertFalse(result["factorizes"])

    def test_e8_gs_applies(self):
        from alpha_ladder_core.anomaly_cancellation import check_green_schwarz_factorization
        result = check_green_schwarz_factorization("E8_x_E8")
        self.assertTrue(result["gs_mechanism_applies"])

    def test_sm_gs_does_not_apply(self):
        from alpha_ladder_core.anomaly_cancellation import check_green_schwarz_factorization
        result = check_green_schwarz_factorization("SU3_x_SU2_x_U1")
        self.assertFalse(result["gs_mechanism_applies"])

    def test_e8_anomaly_condition(self):
        from alpha_ladder_core.anomaly_cancellation import check_green_schwarz_factorization
        result = check_green_schwarz_factorization("E8_x_E8")
        self.assertTrue(result["anomaly_condition_met"])

    def test_so32_anomaly_condition(self):
        from alpha_ladder_core.anomaly_cancellation import check_green_schwarz_factorization
        result = check_green_schwarz_factorization("SO_32")
        self.assertTrue(result["anomaly_condition_met"])

    def test_unknown_group(self):
        from alpha_ladder_core.anomaly_cancellation import check_green_schwarz_factorization
        result = check_green_schwarz_factorization("SU_999")
        self.assertFalse(result["factorizes"])

    def test_e7_factorizes(self):
        from alpha_ladder_core.anomaly_cancellation import check_green_schwarz_factorization
        result = check_green_schwarz_factorization("E7_x_E7")
        self.assertTrue(result["factorizes"])

    def test_description_nonempty(self):
        from alpha_ladder_core.anomaly_cancellation import check_green_schwarz_factorization
        result = check_green_schwarz_factorization("E8_x_E8")
        self.assertIsInstance(result["description"], str)
        self.assertGreater(len(result["description"]), 10)


class TestScanGroups(unittest.TestCase):
    """Tests for scan_anomaly_free_groups."""

    def test_returns_dict(self):
        from alpha_ladder_core.anomaly_cancellation import scan_anomaly_free_groups
        result = scan_anomaly_free_groups()
        self.assertIsInstance(result, dict)

    def test_at_least_three_anomaly_free(self):
        from alpha_ladder_core.anomaly_cancellation import scan_anomaly_free_groups
        result = scan_anomaly_free_groups()
        self.assertGreaterEqual(result["n_anomaly_free"], 3)

    def test_any_contain_sm(self):
        from alpha_ladder_core.anomaly_cancellation import scan_anomaly_free_groups
        result = scan_anomaly_free_groups()
        self.assertTrue(result["any_contain_sm"])

    def test_minimal_group_exists(self):
        from alpha_ladder_core.anomaly_cancellation import scan_anomaly_free_groups
        result = scan_anomaly_free_groups()
        self.assertIsNotNone(result["minimal_group"])

    def test_groups_have_required_keys(self):
        from alpha_ladder_core.anomaly_cancellation import scan_anomaly_free_groups
        result = scan_anomaly_free_groups()
        for g in result["groups"]:
            self.assertIn("name", g)
            self.assertIn("gs_factorizes", g)
            self.assertIn("contains_sm", g)
            self.assertIn("rank", g)

    def test_sm_not_anomaly_free_alone(self):
        from alpha_ladder_core.anomaly_cancellation import scan_anomaly_free_groups
        result = scan_anomaly_free_groups()
        sm_entry = None
        for g in result["groups"]:
            if "SU(3)" in g["name"]:
                sm_entry = g
                break
        self.assertIsNotNone(sm_entry)
        self.assertFalse(sm_entry["gs_factorizes"])

    def test_n_contain_sm_positive(self):
        from alpha_ladder_core.anomaly_cancellation import scan_anomaly_free_groups
        result = scan_anomaly_free_groups()
        self.assertGreater(result["n_contain_sm"], 0)

    def test_honest_assessment_nonempty(self):
        from alpha_ladder_core.anomaly_cancellation import scan_anomaly_free_groups
        result = scan_anomaly_free_groups()
        self.assertGreater(len(result["honest_assessment"]), 20)


class TestSMEmbedding(unittest.TestCase):
    """Tests for check_sm_embedding."""

    def test_e8_contains_sm(self):
        from alpha_ladder_core.anomaly_cancellation import check_sm_embedding
        result = check_sm_embedding("E8_x_E8")
        self.assertTrue(result["contains_sm"])

    def test_so32_contains_sm(self):
        from alpha_ladder_core.anomaly_cancellation import check_sm_embedding
        result = check_sm_embedding("SO_32")
        self.assertTrue(result["contains_sm"])

    def test_e8_embedding_chain_nonempty(self):
        from alpha_ladder_core.anomaly_cancellation import check_sm_embedding
        result = check_sm_embedding("E8_x_E8")
        self.assertGreater(len(result["embedding_chain"]), 0)

    def test_e8_multiple_steps(self):
        from alpha_ladder_core.anomaly_cancellation import check_sm_embedding
        result = check_sm_embedding("E8_x_E8")
        self.assertGreater(result["n_steps"], 1)

    def test_branching_rules_known(self):
        from alpha_ladder_core.anomaly_cancellation import check_sm_embedding
        result = check_sm_embedding("E8_x_E8")
        self.assertTrue(result["branching_rules_known"])

    def test_unknown_group_no_embedding(self):
        from alpha_ladder_core.anomaly_cancellation import check_sm_embedding
        result = check_sm_embedding("SU_999")
        self.assertFalse(result["contains_sm"])

    def test_description_nonempty(self):
        from alpha_ladder_core.anomaly_cancellation import check_sm_embedding
        result = check_sm_embedding("E8_x_E8")
        self.assertIsInstance(result["description"], str)


class TestAlphaLadderConstraints(unittest.TestCase):
    """Tests for compute_anomaly_constraints_on_alpha_ladder."""

    def test_returns_dict(self):
        from alpha_ladder_core.anomaly_cancellation import compute_anomaly_constraints_on_alpha_ladder
        result = compute_anomaly_constraints_on_alpha_ladder()
        self.assertIsInstance(result, dict)

    def test_g_prediction_unaffected(self):
        from alpha_ladder_core.anomaly_cancellation import compute_anomaly_constraints_on_alpha_ladder
        result = compute_anomaly_constraints_on_alpha_ladder()
        self.assertTrue(result["g_prediction_unaffected"])

    def test_gravity_sector_safe(self):
        from alpha_ladder_core.anomaly_cancellation import compute_anomaly_constraints_on_alpha_ladder
        result = compute_anomaly_constraints_on_alpha_ladder()
        self.assertTrue(result["gravity_sector_anomaly_free"])

    def test_matter_constrained(self):
        from alpha_ladder_core.anomaly_cancellation import compute_anomaly_constraints_on_alpha_ladder
        result = compute_anomaly_constraints_on_alpha_ladder()
        self.assertTrue(result["matter_sector_constrained"])

    def test_no_constraints_on_g(self):
        from alpha_ladder_core.anomaly_cancellation import compute_anomaly_constraints_on_alpha_ladder
        result = compute_anomaly_constraints_on_alpha_ladder()
        self.assertFalse(result["constraints_on_g_prediction"])

    def test_detailed_reasoning_nonempty(self):
        from alpha_ladder_core.anomaly_cancellation import compute_anomaly_constraints_on_alpha_ladder
        result = compute_anomaly_constraints_on_alpha_ladder()
        self.assertIsInstance(result["detailed_reasoning"], list)
        self.assertGreater(len(result["detailed_reasoning"]), 3)

    def test_honest_assessment_nonempty(self):
        from alpha_ladder_core.anomaly_cancellation import compute_anomaly_constraints_on_alpha_ladder
        result = compute_anomaly_constraints_on_alpha_ladder()
        self.assertGreater(len(result["honest_assessment"]), 20)


if __name__ == "__main__":
    unittest.main()
