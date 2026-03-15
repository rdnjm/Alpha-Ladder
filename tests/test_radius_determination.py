"""Tests for the radius determination module."""

import unittest


class TestScalingSymmetry(unittest.TestCase):
    """Tests for prove_scaling_symmetry."""

    def test_default_returns_dict(self):
        from alpha_ladder_core.radius_determination import prove_scaling_symmetry
        result = prove_scaling_symmetry()
        self.assertIsInstance(result, dict)

    def test_total_dimensions(self):
        from alpha_ladder_core.radius_determination import prove_scaling_symmetry
        result = prove_scaling_symmetry(d=4, n=2)
        self.assertEqual(result["D"], 6)

    def test_eh_exponent_zero_for_n2(self):
        from alpha_ladder_core.radius_determination import prove_scaling_symmetry
        result = prove_scaling_symmetry(d=4, n=2)
        self.assertEqual(result["eh_scaling_exponent"], 0)

    def test_eh_scale_invariant_for_n2(self):
        from alpha_ladder_core.radius_determination import prove_scaling_symmetry
        result = prove_scaling_symmetry(d=4, n=2)
        self.assertTrue(result["eh_is_scale_invariant"])

    def test_gb_topological_for_n2(self):
        from alpha_ladder_core.radius_determination import prove_scaling_symmetry
        result = prove_scaling_symmetry(d=4, n=2)
        self.assertTrue(result["gb_is_topological"])

    def test_is_scale_invariant_for_n2(self):
        from alpha_ladder_core.radius_determination import prove_scaling_symmetry
        result = prove_scaling_symmetry(d=4, n=2)
        self.assertTrue(result["is_scale_invariant"])

    def test_not_scale_invariant_for_n3(self):
        from alpha_ladder_core.radius_determination import prove_scaling_symmetry
        result = prove_scaling_symmetry(d=4, n=3)
        self.assertFalse(result["is_scale_invariant"])

    def test_eh_exponent_nonzero_for_n1(self):
        from alpha_ladder_core.radius_determination import prove_scaling_symmetry
        result = prove_scaling_symmetry(d=4, n=1)
        self.assertEqual(result["eh_scaling_exponent"], -1)
        self.assertFalse(result["eh_is_scale_invariant"])

    def test_proof_steps_nonempty(self):
        from alpha_ladder_core.radius_determination import prove_scaling_symmetry
        result = prove_scaling_symmetry()
        self.assertIsInstance(result["proof_steps"], list)
        self.assertGreater(len(result["proof_steps"]), 3)

    def test_gb_value_for_n2(self):
        from alpha_ladder_core.radius_determination import prove_scaling_symmetry
        result = prove_scaling_symmetry(d=4, n=2)
        self.assertEqual(result["gb_value"], "4*pi*chi(M_2)")


class TestMechanismCatalog(unittest.TestCase):
    """Tests for catalog_radius_mechanisms."""

    def test_returns_dict(self):
        from alpha_ladder_core.radius_determination import catalog_radius_mechanisms
        result = catalog_radius_mechanisms()
        self.assertIsInstance(result, dict)

    def test_mechanisms_is_list(self):
        from alpha_ladder_core.radius_determination import catalog_radius_mechanisms
        result = catalog_radius_mechanisms()
        self.assertIsInstance(result["mechanisms"], list)

    def test_at_least_four_mechanisms(self):
        from alpha_ladder_core.radius_determination import catalog_radius_mechanisms
        result = catalog_radius_mechanisms()
        self.assertGreaterEqual(result["n_mechanisms"], 4)

    def test_mechanism_has_required_keys(self):
        from alpha_ladder_core.radius_determination import catalog_radius_mechanisms
        result = catalog_radius_mechanisms()
        for m in result["mechanisms"]:
            self.assertIn("name", m)
            self.assertIn("description", m)
            self.assertIn("status", m)
            self.assertIn("computed", m)

    def test_flux_casimir_mechanism_present(self):
        from alpha_ladder_core.radius_determination import catalog_radius_mechanisms
        result = catalog_radius_mechanisms()
        names = [m["name"] for m in result["mechanisms"]]
        self.assertTrue(any("flux" in n.lower() or "casimir" in n.lower() for n in names))

    def test_at_least_one_computed(self):
        from alpha_ladder_core.radius_determination import catalog_radius_mechanisms
        result = catalog_radius_mechanisms()
        self.assertTrue(result["any_computed"])

    def test_none_successful(self):
        from alpha_ladder_core.radius_determination import catalog_radius_mechanisms
        result = catalog_radius_mechanisms()
        self.assertFalse(result["any_successful"])


class TestFluxCasimirBalance(unittest.TestCase):
    """Tests for compute_flux_casimir_balance."""

    def test_returns_dict(self):
        from alpha_ladder_core.radius_determination import compute_flux_casimir_balance
        result = compute_flux_casimir_balance(N=1)
        self.assertIsInstance(result, dict)

    def test_no_solution_for_n1(self):
        from alpha_ladder_core.radius_determination import compute_flux_casimir_balance
        result = compute_flux_casimir_balance(N=1)
        self.assertIsNone(result["a_0_solution"])

    def test_v_min_values_nonempty(self):
        from alpha_ladder_core.radius_determination import compute_flux_casimir_balance
        result = compute_flux_casimir_balance(N=1)
        self.assertIsInstance(result["V_min_values"], list)
        self.assertGreater(len(result["V_min_values"]), 0)

    def test_v_min_values_have_a0(self):
        from alpha_ladder_core.radius_determination import compute_flux_casimir_balance
        result = compute_flux_casimir_balance(N=1)
        for entry in result["V_min_values"]:
            self.assertIn("a_0", entry)
            self.assertIn("V_min", entry)

    def test_is_monotonic(self):
        from alpha_ladder_core.radius_determination import compute_flux_casimir_balance
        result = compute_flux_casimir_balance(N=1)
        self.assertIsInstance(result["is_monotonic"], bool)

    def test_honest_result_nonempty(self):
        from alpha_ladder_core.radius_determination import compute_flux_casimir_balance
        result = compute_flux_casimir_balance(N=1)
        self.assertIsInstance(result["honest_result"], str)
        self.assertGreater(len(result["honest_result"]), 10)

    def test_n_preserved(self):
        from alpha_ladder_core.radius_determination import compute_flux_casimir_balance
        result = compute_flux_casimir_balance(N=3)
        self.assertEqual(result["N"], 3)

    def test_method_is_flux_casimir(self):
        from alpha_ladder_core.radius_determination import compute_flux_casimir_balance
        result = compute_flux_casimir_balance(N=1)
        self.assertEqual(result["method"], "flux_casimir_balance")

    def test_no_solution_for_n5(self):
        from alpha_ladder_core.radius_determination import compute_flux_casimir_balance
        result = compute_flux_casimir_balance(N=5)
        self.assertIsNone(result["a_0_solution"])

    def test_all_minima_exist(self):
        from alpha_ladder_core.radius_determination import compute_flux_casimir_balance
        result = compute_flux_casimir_balance(N=1)
        for entry in result["V_min_values"]:
            self.assertTrue(entry["minimum_exists"])


class TestRadiusLandscape(unittest.TestCase):
    """Tests for compute_radius_landscape."""

    def test_returns_dict(self):
        from alpha_ladder_core.radius_determination import compute_radius_landscape
        result = compute_radius_landscape()
        self.assertIsInstance(result, dict)

    def test_landscape_is_list(self):
        from alpha_ladder_core.radius_determination import compute_radius_landscape
        result = compute_radius_landscape()
        self.assertIsInstance(result["landscape"], list)

    def test_default_ten_entries(self):
        from alpha_ladder_core.radius_determination import compute_radius_landscape
        result = compute_radius_landscape()
        self.assertEqual(len(result["landscape"]), 10)

    def test_custom_n_values(self):
        from alpha_ladder_core.radius_determination import compute_radius_landscape
        result = compute_radius_landscape(N_values=[1, 2, 3])
        self.assertEqual(len(result["landscape"]), 3)

    def test_landscape_entries_have_keys(self):
        from alpha_ladder_core.radius_determination import compute_radius_landscape
        result = compute_radius_landscape(N_values=[1])
        entry = result["landscape"][0]
        self.assertIn("N", entry)
        self.assertIn("sigma_0", entry)
        self.assertIn("m_phi_eV", entry)
        self.assertIn("classification", entry)
        self.assertIn("minimum_exists", entry)

    def test_any_testable_is_bool(self):
        from alpha_ladder_core.radius_determination import compute_radius_landscape
        result = compute_radius_landscape()
        self.assertIsInstance(result["any_testable"], bool)

    def test_all_invisible_at_planck(self):
        from alpha_ladder_core.radius_determination import compute_radius_landscape
        result = compute_radius_landscape()
        for entry in result["landscape"]:
            self.assertEqual(entry["classification"], "invisible_planck")

    def test_n_values_preserved(self):
        from alpha_ladder_core.radius_determination import compute_radius_landscape
        result = compute_radius_landscape(N_values=[2, 4, 6])
        self.assertEqual(result["N_values"], [2, 4, 6])

    def test_honest_assessment_nonempty(self):
        from alpha_ladder_core.radius_determination import compute_radius_landscape
        result = compute_radius_landscape()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 10)

    def test_all_minima_exist(self):
        from alpha_ladder_core.radius_determination import compute_radius_landscape
        result = compute_radius_landscape()
        for entry in result["landscape"]:
            self.assertTrue(entry["minimum_exists"])


class TestSummarize(unittest.TestCase):
    """Tests for summarize_radius_determination."""

    def test_returns_dict(self):
        from alpha_ladder_core.radius_determination import summarize_radius_determination
        result = summarize_radius_determination()
        self.assertIsInstance(result, dict)

    def test_has_scaling_symmetry(self):
        from alpha_ladder_core.radius_determination import summarize_radius_determination
        result = summarize_radius_determination()
        self.assertIn("scaling_symmetry", result)

    def test_has_mechanisms(self):
        from alpha_ladder_core.radius_determination import summarize_radius_determination
        result = summarize_radius_determination()
        self.assertIn("mechanisms", result)

    def test_has_flux_casimir_balance(self):
        from alpha_ladder_core.radius_determination import summarize_radius_determination
        result = summarize_radius_determination()
        self.assertIn("flux_casimir_balance", result)

    def test_has_landscape(self):
        from alpha_ladder_core.radius_determination import summarize_radius_determination
        result = summarize_radius_determination()
        self.assertIn("landscape", result)

    def test_a0_not_determined(self):
        from alpha_ladder_core.radius_determination import summarize_radius_determination
        result = summarize_radius_determination()
        self.assertFalse(result["a_0_determined"])

    def test_scaling_symmetry_proven(self):
        from alpha_ladder_core.radius_determination import summarize_radius_determination
        result = summarize_radius_determination()
        self.assertTrue(result["scaling_symmetry"]["is_scale_invariant"])

    def test_overall_assessment_nonempty(self):
        from alpha_ladder_core.radius_determination import summarize_radius_determination
        result = summarize_radius_determination()
        self.assertIsInstance(result["overall_assessment"], str)
        self.assertGreater(len(result["overall_assessment"]), 50)

    def test_honest_assessment_nonempty(self):
        from alpha_ladder_core.radius_determination import summarize_radius_determination
        result = summarize_radius_determination()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 50)

    def test_no_false_determination_claim(self):
        from alpha_ladder_core.radius_determination import summarize_radius_determination
        result = summarize_radius_determination()
        assessment = result["overall_assessment"].lower()
        self.assertNotIn("determined from first principles", assessment)


if __name__ == "__main__":
    unittest.main()
