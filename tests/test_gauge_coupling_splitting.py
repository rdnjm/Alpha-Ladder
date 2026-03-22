"""Tests for gauge_coupling_splitting.py -- SO(3)/SO(2) Gauge Fields on S^2."""
import unittest
import math



ALPHA_EM = 1.0 / 137.035999084
PHI_VEV = 0.25 * math.log(4.0 * math.pi * ALPHA_EM)


class TestGaugeKineticMatrix(unittest.TestCase):
    """Tests for gauge_kinetic_matrix()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_coupling_splitting import gauge_kinetic_matrix
        result = gauge_kinetic_matrix()
        self.assertIsInstance(result, dict)

    def test_eigenvalues_0_1_1(self):
        from alpha_ladder_core.gauge_coupling_splitting import gauge_kinetic_matrix
        result = gauge_kinetic_matrix()
        self.assertEqual(result["eigenvalues"], [0.0, 1.0, 1.0])

    def test_rank_2(self):
        from alpha_ladder_core.gauge_coupling_splitting import gauge_kinetic_matrix
        result = gauge_kinetic_matrix()
        self.assertEqual(result["rank"], 2)

    def test_matrix_is_3x3(self):
        from alpha_ladder_core.gauge_coupling_splitting import gauge_kinetic_matrix
        result = gauge_kinetic_matrix()
        self.assertEqual(len(result["matrix"]), 3)
        for row in result["matrix"]:
            self.assertEqual(len(row), 3)

    def test_matrix_diagonal(self):
        from alpha_ladder_core.gauge_coupling_splitting import gauge_kinetic_matrix
        result = gauge_kinetic_matrix()
        K = result["matrix"]
        self.assertAlmostEqual(K[0][0], 0.0, places=10)
        self.assertAlmostEqual(K[1][1], 1.0, places=10)
        self.assertAlmostEqual(K[2][2], 1.0, places=10)

    def test_matrix_off_diagonal_zero(self):
        from alpha_ladder_core.gauge_coupling_splitting import gauge_kinetic_matrix
        result = gauge_kinetic_matrix()
        K = result["matrix"]
        for i in range(3):
            for j in range(3):
                if i != j:
                    self.assertAlmostEqual(K[i][j], 0.0, places=10)

    def test_so2_direction_is_0(self):
        from alpha_ladder_core.gauge_coupling_splitting import gauge_kinetic_matrix
        result = gauge_kinetic_matrix()
        self.assertEqual(result["so2_direction"], 0)

    def test_coset_directions(self):
        from alpha_ladder_core.gauge_coupling_splitting import gauge_kinetic_matrix
        result = gauge_kinetic_matrix()
        self.assertEqual(result["coset_directions"], [1, 2])

    def test_eigenvectors_present(self):
        from alpha_ladder_core.gauge_coupling_splitting import gauge_kinetic_matrix
        result = gauge_kinetic_matrix()
        self.assertEqual(len(result["eigenvectors"]), 3)

    def test_description_present(self):
        from alpha_ladder_core.gauge_coupling_splitting import gauge_kinetic_matrix
        result = gauge_kinetic_matrix()
        self.assertIsInstance(result["description"], str)
        self.assertGreater(len(result["description"]), 50)


class TestTreeLevelCouplings(unittest.TestCase):
    """Tests for tree_level_couplings()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_coupling_splitting import tree_level_couplings
        result = tree_level_couplings()
        self.assertIsInstance(result, dict)

    def test_alpha_1_equals_alpha_2(self):
        """No splitting at tree level: alpha_1 = alpha_2."""
        from alpha_ladder_core.gauge_coupling_splitting import tree_level_couplings
        result = tree_level_couplings()
        self.assertAlmostEqual(result["alpha_1"], result["alpha_2"], places=15)

    def test_are_equal_flag(self):
        from alpha_ladder_core.gauge_coupling_splitting import tree_level_couplings
        result = tree_level_couplings()
        self.assertTrue(result["are_equal"])

    def test_alpha_matches_alpha_em(self):
        from alpha_ladder_core.gauge_coupling_splitting import tree_level_couplings
        result = tree_level_couplings()
        self.assertAlmostEqual(result["alpha_1"], ALPHA_EM, places=12)

    def test_g_values_positive(self):
        from alpha_ladder_core.gauge_coupling_splitting import tree_level_couplings
        result = tree_level_couplings()
        self.assertGreater(result["g_1"], 0.0)
        self.assertGreater(result["g_2"], 0.0)

    def test_g_equals(self):
        from alpha_ladder_core.gauge_coupling_splitting import tree_level_couplings
        result = tree_level_couplings()
        self.assertAlmostEqual(result["g_1"], result["g_2"], places=15)

    def test_custom_phi_vev(self):
        from alpha_ladder_core.gauge_coupling_splitting import tree_level_couplings
        result = tree_level_couplings(phi_vev=0.0)
        expected = 1.0 / (4.0 * math.pi)
        self.assertAlmostEqual(result["alpha_1"], expected, places=12)

    def test_reason_is_string(self):
        from alpha_ladder_core.gauge_coupling_splitting import tree_level_couplings
        result = tree_level_couplings()
        self.assertIsInstance(result["reason"], str)

    def test_g_alpha_relation(self):
        """g^2 / (4*pi) = alpha."""
        from alpha_ladder_core.gauge_coupling_splitting import tree_level_couplings
        result = tree_level_couplings()
        self.assertAlmostEqual(
            result["g_1"] ** 2 / (4.0 * math.pi),
            result["alpha_1"],
            places=12,
        )


class TestOneLoopBetaFunctions(unittest.TestCase):
    """Tests for one_loop_beta_functions()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_coupling_splitting import one_loop_beta_functions
        result = one_loop_beta_functions()
        self.assertIsInstance(result, dict)

    def test_b2_zero(self):
        from alpha_ladder_core.gauge_coupling_splitting import one_loop_beta_functions
        result = one_loop_beta_functions()
        self.assertEqual(result["b_2"], 0.0)

    def test_b3_zero(self):
        from alpha_ladder_core.gauge_coupling_splitting import one_loop_beta_functions
        result = one_loop_beta_functions()
        self.assertEqual(result["b_3"], 0.0)

    def test_are_equal(self):
        from alpha_ladder_core.gauge_coupling_splitting import one_loop_beta_functions
        result = one_loop_beta_functions()
        self.assertTrue(result["are_equal"])

    def test_contributions_present(self):
        from alpha_ladder_core.gauge_coupling_splitting import one_loop_beta_functions
        result = one_loop_beta_functions()
        contribs = result["contributions"]
        for key in ["graviton_loops", "dilaton_loops", "contact_interaction",
                     "kk_tower_below_mkk", "kk_tower_above_mkk"]:
            self.assertIn(key, contribs)

    def test_graviton_contribution_zero(self):
        from alpha_ladder_core.gauge_coupling_splitting import one_loop_beta_functions
        result = one_loop_beta_functions()
        self.assertEqual(result["contributions"]["graviton_loops"]["b"], 0.0)

    def test_dilaton_contribution_zero(self):
        from alpha_ladder_core.gauge_coupling_splitting import one_loop_beta_functions
        result = one_loop_beta_functions()
        self.assertEqual(result["contributions"]["dilaton_loops"]["b"], 0.0)

    def test_m_kk_positive(self):
        from alpha_ladder_core.gauge_coupling_splitting import one_loop_beta_functions
        result = one_loop_beta_functions()
        self.assertGreater(result["m_kk_eV"], 0.0)

    def test_alpha_at_tree(self):
        from alpha_ladder_core.gauge_coupling_splitting import one_loop_beta_functions
        result = one_loop_beta_functions()
        self.assertAlmostEqual(result["alpha_at_tree"], ALPHA_EM, places=12)


class TestCouplingEvolution(unittest.TestCase):
    """Tests for coupling_evolution()."""

    def test_returns_list(self):
        from alpha_ladder_core.gauge_coupling_splitting import coupling_evolution
        result = coupling_evolution()
        self.assertIsInstance(result, list)

    def test_default_scales(self):
        from alpha_ladder_core.gauge_coupling_splitting import coupling_evolution
        result = coupling_evolution()
        self.assertEqual(len(result), 9)

    def test_all_equal(self):
        from alpha_ladder_core.gauge_coupling_splitting import coupling_evolution
        result = coupling_evolution()
        for entry in result:
            self.assertTrue(entry["are_equal"])
            self.assertAlmostEqual(entry["alpha_2"], entry["alpha_3"], places=15)

    def test_no_running_below_mkk(self):
        """Below m_KK, alpha should equal tree-level value."""
        from alpha_ladder_core.gauge_coupling_splitting import coupling_evolution
        result = coupling_evolution()
        alpha_tree = result[0]["alpha_2"]
        for entry in result:
            self.assertAlmostEqual(entry["alpha_2"], alpha_tree, places=15)

    def test_custom_scales(self):
        from alpha_ladder_core.gauge_coupling_splitting import coupling_evolution
        scales = [("low", 1.0), ("high", 1e12)]
        result = coupling_evolution(mu_values=scales)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0]["label"], "low")

    def test_regime_labels(self):
        from alpha_ladder_core.gauge_coupling_splitting import coupling_evolution
        result = coupling_evolution()
        for entry in result:
            self.assertIn("regime", entry)
            self.assertIsInstance(entry["regime"], str)

    def test_n_active_kk_increases(self):
        """Higher energies should have more active KK modes."""
        from alpha_ladder_core.gauge_coupling_splitting import coupling_evolution
        result = coupling_evolution()
        low_kk = result[0]["n_active_KK"]
        high_kk = result[-1]["n_active_KK"]
        self.assertGreaterEqual(high_kk, low_kk)

    def test_l_max_present(self):
        from alpha_ladder_core.gauge_coupling_splitting import coupling_evolution
        result = coupling_evolution()
        for entry in result:
            self.assertIn("l_max", entry)


class TestUnificationAnalysis(unittest.TestCase):
    """Tests for unification_analysis()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_coupling_splitting import unification_analysis
        result = unification_analysis()
        self.assertIsInstance(result, dict)

    def test_unification_scale_positive(self):
        from alpha_ladder_core.gauge_coupling_splitting import unification_analysis
        result = unification_analysis()
        self.assertGreater(result["unification_scale_eV"], 0.0)

    def test_coupling_at_unification_is_alpha_em(self):
        from alpha_ladder_core.gauge_coupling_splitting import unification_analysis
        result = unification_analysis()
        self.assertAlmostEqual(
            result["coupling_at_unification"], ALPHA_EM, places=12
        )

    def test_much_below_gut_scale(self):
        from alpha_ladder_core.gauge_coupling_splitting import unification_analysis
        result = unification_analysis()
        self.assertGreater(result["ratio_to_GUT"], 1e10)

    def test_orders_of_magnitude_apart(self):
        from alpha_ladder_core.gauge_coupling_splitting import unification_analysis
        result = unification_analysis()
        orders = result["comparison_with_GUT"]["orders_of_magnitude_apart"]
        self.assertGreater(orders, 10)

    def test_description_present(self):
        from alpha_ladder_core.gauge_coupling_splitting import unification_analysis
        result = unification_analysis()
        self.assertIsInstance(result["description"], str)
        self.assertIn("geometric", result["description"].lower())


class TestTwoBraneSplitting(unittest.TestCase):
    """Tests for two_brane_splitting()."""

    def test_returns_list(self):
        from alpha_ladder_core.gauge_coupling_splitting import two_brane_splitting
        result = two_brane_splitting()
        self.assertIsInstance(result, list)

    def test_default_8_theta_values(self):
        from alpha_ladder_core.gauge_coupling_splitting import two_brane_splitting
        result = two_brane_splitting()
        self.assertEqual(len(result), 8)

    def test_no_splitting_at_pole(self):
        """At theta_2 = 0 (same pole as brane 1), no splitting."""
        from alpha_ladder_core.gauge_coupling_splitting import two_brane_splitting
        result = two_brane_splitting()
        pole_entry = result[0]
        self.assertAlmostEqual(pole_entry["theta_2"], 0.0, places=10)
        self.assertAlmostEqual(pole_entry["splitting"], 0.0, places=10)
        self.assertTrue(pole_entry["are_equal_brane2"])

    def test_no_splitting_at_south_pole(self):
        """At theta_2 = pi (south pole, antipodal), no splitting."""
        from alpha_ladder_core.gauge_coupling_splitting import two_brane_splitting
        result = two_brane_splitting()
        south = result[-1]
        self.assertAlmostEqual(south["theta_2"], math.pi, places=10)
        self.assertAlmostEqual(south["splitting"], 0.0, places=10)
        self.assertTrue(south["are_equal_brane2"])

    def test_max_splitting_at_equator(self):
        """Maximum splitting at theta_2 = pi/2 (equator)."""
        from alpha_ladder_core.gauge_coupling_splitting import two_brane_splitting
        result = two_brane_splitting()
        splittings = [r["splitting"] for r in result]
        max_idx = splittings.index(max(splittings))
        self.assertAlmostEqual(
            result[max_idx]["theta_2"], math.pi / 2.0, places=10
        )

    def test_equator_splitting_value(self):
        """At equator: sin^2(pi/2)/(1+cos^2(pi/2)) = 1/1 = 1.0."""
        from alpha_ladder_core.gauge_coupling_splitting import two_brane_splitting
        result = two_brane_splitting()
        equator = [r for r in result if abs(r["theta_2"] - math.pi / 2.0) < 0.01][0]
        self.assertAlmostEqual(equator["splitting"], 1.0, places=10)

    def test_brane1_always_equal(self):
        """On Brane 1 (north pole), couplings always equal."""
        from alpha_ladder_core.gauge_coupling_splitting import two_brane_splitting
        result = two_brane_splitting()
        for r in result:
            self.assertTrue(r["are_equal_brane1"])
            self.assertAlmostEqual(
                r["alpha_2_brane1"], r["alpha_3_brane1"], places=15
            )

    def test_K_squared_at_equator(self):
        """At theta=pi/2: K2^2 = cos^2(pi/2) = 0, K3^2 = 1."""
        from alpha_ladder_core.gauge_coupling_splitting import two_brane_splitting
        result = two_brane_splitting()
        equator = [r for r in result if abs(r["theta_2"] - math.pi / 2.0) < 0.01][0]
        self.assertAlmostEqual(equator["K2_squared"], 0.0, places=10)
        self.assertAlmostEqual(equator["K3_squared"], 1.0, places=10)

    def test_custom_theta_values(self):
        from alpha_ladder_core.gauge_coupling_splitting import two_brane_splitting
        thetas = [0.0, math.pi / 4.0, math.pi]
        result = two_brane_splitting(theta_2_values=thetas)
        self.assertEqual(len(result), 3)

    def test_splitting_symmetric_about_equator(self):
        """Splitting at pi/3 and 2*pi/3 should be the same."""
        from alpha_ladder_core.gauge_coupling_splitting import two_brane_splitting
        thetas = [math.pi / 3.0, 2.0 * math.pi / 3.0]
        result = two_brane_splitting(theta_2_values=thetas)
        # Both have |cos^2(theta)| but cos(pi/3) = 0.5, cos(2pi/3) = -0.5
        # K2^2 = cos^2(theta), so both give K2^2 = 0.25
        self.assertAlmostEqual(
            result[0]["splitting"], result[1]["splitting"], places=10
        )

    def test_degrees_conversion(self):
        from alpha_ladder_core.gauge_coupling_splitting import two_brane_splitting
        result = two_brane_splitting()
        for r in result:
            self.assertAlmostEqual(
                r["theta_2_deg"], math.degrees(r["theta_2"]), places=6
            )


class TestElectroweakAnalogy(unittest.TestCase):
    """Tests for electroweak_analogy()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_coupling_splitting import electroweak_analogy
        result = electroweak_analogy()
        self.assertIsInstance(result, dict)

    def test_structural_analogy_true(self):
        from alpha_ladder_core.gauge_coupling_splitting import electroweak_analogy
        result = electroweak_analogy()
        self.assertTrue(result["structural_analogy"])

    def test_physical_analogy_false(self):
        from alpha_ladder_core.gauge_coupling_splitting import electroweak_analogy
        result = electroweak_analogy()
        self.assertFalse(result["physical_analogy"])

    def test_comparison_has_entries(self):
        from alpha_ladder_core.gauge_coupling_splitting import electroweak_analogy
        result = electroweak_analogy()
        self.assertEqual(len(result["comparison"]), 8)

    def test_key_difference_present(self):
        from alpha_ladder_core.gauge_coupling_splitting import electroweak_analogy
        result = electroweak_analogy()
        self.assertIsInstance(result["key_difference"], str)
        self.assertGreater(len(result["key_difference"]), 50)

    def test_so3_generators(self):
        from alpha_ladder_core.gauge_coupling_splitting import electroweak_analogy
        result = electroweak_analogy()
        self.assertEqual(result["so3_group_theory"]["generators"], 3)

    def test_su2_generators(self):
        from alpha_ladder_core.gauge_coupling_splitting import electroweak_analogy
        result = electroweak_analogy()
        self.assertEqual(result["su2_group_theory"]["generators"], 4)

    def test_both_unbroken_u1(self):
        """Both breaking patterns leave a U(1)."""
        from alpha_ladder_core.gauge_coupling_splitting import electroweak_analogy
        result = electroweak_analogy()
        self.assertEqual(result["so3_group_theory"]["unbroken_generators"], 1)
        self.assertEqual(result["su2_group_theory"]["unbroken_generators"], 1)


class TestSummarizeGaugeSplitting(unittest.TestCase):
    """Tests for summarize_gauge_splitting()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_coupling_splitting import summarize_gauge_splitting
        result = summarize_gauge_splitting()
        self.assertIsInstance(result, dict)

    def test_all_sections_present(self):
        from alpha_ladder_core.gauge_coupling_splitting import summarize_gauge_splitting
        result = summarize_gauge_splitting()
        for key in ["kinetic_matrix", "tree_level", "beta_functions",
                     "evolution", "unification", "two_brane", "electroweak",
                     "summary_text"]:
            self.assertIn(key, result)

    def test_no_splitting_single_brane(self):
        from alpha_ladder_core.gauge_coupling_splitting import summarize_gauge_splitting
        result = summarize_gauge_splitting()
        self.assertFalse(result["any_splitting_single_brane"])

    def test_max_two_brane_splitting(self):
        from alpha_ladder_core.gauge_coupling_splitting import summarize_gauge_splitting
        result = summarize_gauge_splitting()
        self.assertAlmostEqual(result["max_two_brane_splitting"], 1.0, places=6)

    def test_summary_text_is_string(self):
        from alpha_ladder_core.gauge_coupling_splitting import summarize_gauge_splitting
        result = summarize_gauge_splitting()
        self.assertIsInstance(result["summary_text"], str)
        self.assertGreater(len(result["summary_text"]), 500)

    def test_summary_mentions_no_splitting(self):
        from alpha_ladder_core.gauge_coupling_splitting import summarize_gauge_splitting
        result = summarize_gauge_splitting()
        self.assertIn("NO splitting", result["summary_text"])


class TestHelpers(unittest.TestCase):
    """Tests for helper functions."""

    def test_eigenvalues_identity(self):
        from alpha_ladder_core.gauge_coupling_splitting import _eigenvalues_3x3_symmetric
        I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        evals = _eigenvalues_3x3_symmetric(I)
        for ev in evals:
            self.assertAlmostEqual(ev, 1.0, places=8)

    def test_eigenvalues_zero_matrix(self):
        from alpha_ladder_core.gauge_coupling_splitting import _eigenvalues_3x3_symmetric
        Z = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        evals = _eigenvalues_3x3_symmetric(Z)
        for ev in evals:
            self.assertAlmostEqual(ev, 0.0, places=8)

    def test_eigenvalues_diagonal(self):
        from alpha_ladder_core.gauge_coupling_splitting import _eigenvalues_3x3_symmetric
        D = [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
        evals = _eigenvalues_3x3_symmetric(D)
        self.assertAlmostEqual(evals[0], 1.0, places=8)
        self.assertAlmostEqual(evals[1], 2.0, places=8)
        self.assertAlmostEqual(evals[2], 3.0, places=8)

    def test_mat_multiply(self):
        from alpha_ladder_core.gauge_coupling_splitting import _mat_multiply
        A = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        B = [[2, 0, 0], [0, 3, 0], [0, 0, 4]]
        C = _mat_multiply(A, B)
        self.assertAlmostEqual(C[0][0], 2.0)
        self.assertAlmostEqual(C[1][1], 3.0)
        self.assertAlmostEqual(C[2][2], 4.0)

    def test_mat_transpose(self):
        from alpha_ladder_core.gauge_coupling_splitting import _mat_transpose
        A = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        AT = _mat_transpose(A)
        self.assertEqual(AT[0][1], 4)
        self.assertEqual(AT[1][0], 2)

    def test_m_kk_positive(self):
        from alpha_ladder_core.gauge_coupling_splitting import _m_kk
        m = _m_kk(28e-6)
        self.assertGreater(m, 0.0)

    def test_R_phys_positive(self):
        from alpha_ladder_core.gauge_coupling_splitting import _R_phys
        R = _R_phys(28e-6)
        self.assertGreater(R, 0.0)


if __name__ == "__main__":
    unittest.main()
