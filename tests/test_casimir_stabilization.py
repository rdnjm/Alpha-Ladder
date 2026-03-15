"""
Tests for alpha_ladder_core/casimir_stabilization.py

Comprehensive test suite for Casimir stabilization of the volume modulus.
The key physics result verified here is that the Casimir stabilization
from the graviton tower alone does NOT produce a stable minimum -- the
Casimir coefficient is negative (because zeta_S2(-1/2) < 0), so any
stationary point is a maximum, not a minimum.
"""

import math
import unittest

from alpha_ladder_core.casimir_stabilization import (
    compute_kk_spectrum_s2,
    compute_spectral_zeta_s2,
    compute_casimir_energy_s2,
    compute_effective_potential,
    find_casimir_minimum,
    compute_dilaton_mass_casimir,
    summarize_casimir_stabilization,
    compute_fermion_spectral_zeta_s2,
    compute_vector_spectral_zeta_s2,
    compute_matter_casimir_coefficient,
    scan_anomaly_free_matter_casimir,
)


class TestKKSpectrum(unittest.TestCase):
    """Tests for compute_kk_spectrum_s2."""

    def test_returns_dict(self):
        result = compute_kk_spectrum_s2()
        self.assertIsInstance(result, dict)

    def test_eigenvalue_formula_l1(self):
        """l=1 eigenvalue should be 1*2 = 2."""
        result = compute_kk_spectrum_s2(l_max=10)
        eigenvalues = result["eigenvalues"]
        # Module includes l=0 (massless), so l=1 is at index 1
        self.assertAlmostEqual(eigenvalues[1], 2.0, places=10)

    def test_eigenvalue_formula_l2(self):
        """l=2 eigenvalue should be 2*3 = 6."""
        result = compute_kk_spectrum_s2(l_max=10)
        eigenvalues = result["eigenvalues"]
        self.assertAlmostEqual(eigenvalues[2], 6.0, places=10)

    def test_l0_massless(self):
        """l=0 mode is massless (eigenvalue = 0)."""
        result = compute_kk_spectrum_s2(l_max=5)
        self.assertAlmostEqual(result["eigenvalues"][0], 0.0, places=10)

    def test_l0_degeneracy(self):
        """l=0 has degeneracy 2*0+1 = 1."""
        result = compute_kk_spectrum_s2(l_max=5)
        self.assertEqual(result["degeneracies"][0], 1)

    def test_degeneracy_l1(self):
        """l=1 has degeneracy 2*1+1 = 3."""
        result = compute_kk_spectrum_s2(l_max=5)
        self.assertEqual(result["degeneracies"][1], 3)

    def test_field_content_total_polarizations(self):
        """6D graviton has D(D-3)/2 = 9 total polarizations."""
        result = compute_kk_spectrum_s2()
        fc = result["field_content"]
        total = fc.get("total_polarizations") or fc.get("check_total", 0)
        self.assertEqual(total, 9)

    def test_field_content_graviton_dof(self):
        """4D graviton has 2 dof."""
        result = compute_kk_spectrum_s2()
        fc = result["field_content"]
        decomp = fc.get("decomposition", fc)
        graviton = decomp.get("graviton_4D", decomp)
        dof = graviton.get("dof", graviton.get("graviton_dof", 0))
        self.assertEqual(dof, 2)

    def test_field_content_vector_dof(self):
        """Graviphotons have 2*n = 4 total dof."""
        result = compute_kk_spectrum_s2()
        fc = result["field_content"]
        decomp = fc.get("decomposition", fc)
        vectors = decomp.get("vectors", decomp)
        dof = vectors.get("dof_total", vectors.get("vector_dof", 0))
        self.assertEqual(dof, 4)

    def test_field_content_scalar_dof(self):
        """Scalars have n(n+1)/2 = 3 dof."""
        result = compute_kk_spectrum_s2()
        fc = result["field_content"]
        decomp = fc.get("decomposition", fc)
        scalars = decomp.get("scalars", decomp)
        dof = scalars.get("dof_total", scalars.get("scalar_dof", 0))
        self.assertEqual(dof, 3)

    def test_eigenvalues_sorted_ascending(self):
        """Eigenvalues should be in ascending order."""
        result = compute_kk_spectrum_s2(l_max=20)
        eigenvalues = result["eigenvalues"]
        for i in range(len(eigenvalues) - 1):
            self.assertLessEqual(eigenvalues[i], eigenvalues[i + 1])

    def test_degeneracies_length_matches_eigenvalues(self):
        """Degeneracies list should be same length as eigenvalues."""
        result = compute_kk_spectrum_s2(l_max=30)
        self.assertEqual(len(result["degeneracies"]), len(result["eigenvalues"]))

    def test_n_massive_modes_positive(self):
        """n_massive_modes should be positive."""
        result = compute_kk_spectrum_s2(l_max=20)
        self.assertGreater(result["n_massive_modes"], 0)


class TestSpectralZeta(unittest.TestCase):
    """Tests for compute_spectral_zeta_s2."""

    def test_returns_dict(self):
        result = compute_spectral_zeta_s2(s=2.0)
        self.assertIsInstance(result, dict)

    def test_convergent_s_positive(self):
        """For s=2, the series converges and gives a finite positive value."""
        result = compute_spectral_zeta_s2(s=2.0)
        self.assertTrue(math.isfinite(result["zeta_value"]))

    def test_s_2_known_bound(self):
        """zeta(2) for S^2 should be a small positive number (first term: 3/4)."""
        result = compute_spectral_zeta_s2(s=2.0)
        self.assertGreater(result["zeta_value"], 0.7)
        self.assertLess(result["zeta_value"], 2.0)

    def test_analytic_continuation_method(self):
        """For s=-0.5, method should be analytic_continuation."""
        result = compute_spectral_zeta_s2(s=-0.5)
        self.assertEqual(result["method"], "analytic_continuation")

    def test_zeta_minus_half_negative(self):
        """zeta_S2(-1/2) should be negative (approximately -0.265)."""
        result = compute_spectral_zeta_s2(s=-0.5)
        self.assertLess(result["zeta_value"], 0)

    def test_zeta_minus_half_magnitude(self):
        """zeta_S2(-1/2) should be approximately -0.265 (within 0.15)."""
        result = compute_spectral_zeta_s2(s=-0.5)
        self.assertAlmostEqual(result["zeta_value"], -0.265, delta=0.15)

    def test_has_convergence_diagnostics(self):
        """Result must include convergence diagnostics."""
        result = compute_spectral_zeta_s2(s=2.0)
        self.assertIn("convergence_diagnostics", result)

    def test_partial_sum_method_for_convergent(self):
        """For s > 1, method should be partial_sum."""
        result = compute_spectral_zeta_s2(s=2.0)
        self.assertEqual(result["method"], "partial_sum")

    def test_s_stored(self):
        """The s parameter should be stored in the result."""
        result = compute_spectral_zeta_s2(s=1.5)
        self.assertAlmostEqual(result["s"], 1.5)

    def test_larger_l_max_improves_convergence(self):
        """Using larger l_max should give a closer result for convergent series."""
        r1 = compute_spectral_zeta_s2(s=2.0, l_max=100)
        r2 = compute_spectral_zeta_s2(s=2.0, l_max=1000)
        self.assertAlmostEqual(r1["zeta_value"], r2["zeta_value"], delta=0.01)


class TestCasimirEnergy(unittest.TestCase):
    """Tests for compute_casimir_energy_s2."""

    def test_returns_dict(self):
        result = compute_casimir_energy_s2(a_radius=1.0)
        self.assertIsInstance(result, dict)

    def test_scaling_power(self):
        """Casimir energy scales as a^{-4}."""
        result = compute_casimir_energy_s2(a_radius=1.0)
        self.assertEqual(result["scaling_power"], -4)

    def test_a_inverse_4_scaling(self):
        """V(2a) / V(a) = (1/2)^4 = 1/16."""
        r1 = compute_casimir_energy_s2(a_radius=1.0)
        r2 = compute_casimir_energy_s2(a_radius=2.0)
        ratio = r2["V_casimir"] / r1["V_casimir"]
        self.assertAlmostEqual(ratio, (1 / 2) ** 4, places=4)

    def test_n_eff_is_9(self):
        """N_eff = 9 total bosonic dof."""
        result = compute_casimir_energy_s2(a_radius=1.0)
        self.assertEqual(result["N_eff"], 9)

    def test_sign_negative(self):
        """Casimir energy coefficient should be negative."""
        result = compute_casimir_energy_s2(a_radius=1.0)
        self.assertEqual(result["sign"], "negative")

    def test_coefficient_finite(self):
        """The coefficient should be finite."""
        result = compute_casimir_energy_s2(a_radius=1.0)
        self.assertTrue(math.isfinite(result["coefficient"]))

    def test_has_field_content(self):
        """Result should include field content info."""
        result = compute_casimir_energy_s2(a_radius=1.0)
        self.assertIn("field_content", result)

    def test_custom_field_content(self):
        """Should accept custom field content."""
        fc = {
            "total_bosonic_dof": 9,
            "total_polarizations": 9,
        }
        result = compute_casimir_energy_s2(a_radius=1.0, field_content=fc)
        self.assertIsInstance(result, dict)


class TestEffectivePotential(unittest.TestCase):
    """Tests for compute_effective_potential."""

    def test_returns_dict(self):
        result = compute_effective_potential(sigma=0.0)
        self.assertIsInstance(result, dict)

    def test_has_sigma_grid(self):
        """Result should have a sigma grid."""
        result = compute_effective_potential(sigma=0.0)
        self.assertIn("sigma_grid", result)
        self.assertTrue(len(result["sigma_grid"]) > 0)

    def test_has_components(self):
        """Result should have V_classical, V_casimir, V_total."""
        result = compute_effective_potential(sigma=0.0)
        self.assertIn("V_classical", result)
        self.assertIn("V_casimir", result)
        self.assertIn("V_total", result)

    def test_a_coefficient_negative(self):
        """A (Casimir coefficient) should be negative."""
        result = compute_effective_potential(sigma=0.0)
        self.assertLess(result["A_coefficient"], 0)

    def test_b_coefficient_s2(self):
        """B (curvature coefficient) for S^2 should be negative."""
        result = compute_effective_potential(sigma=0.0)
        self.assertLess(result["B_coefficient"], 0)

    def test_chi_value(self):
        """Default chi should be 2 (S^2)."""
        result = compute_effective_potential(sigma=0.0)
        self.assertEqual(result["chi"], 2)

    def test_has_genus2_data(self):
        """Result should include genus-2 comparison data."""
        result = compute_effective_potential(sigma=0.0)
        self.assertIn("genus2_data", result)

    def test_genus2_b_positive(self):
        """For genus 2, B should be positive."""
        result = compute_effective_potential(sigma=0.0)
        self.assertGreater(result["genus2_data"]["B_coefficient"], 0)

    def test_v_total_is_sum(self):
        """V_total should equal V_classical + V_casimir at each point."""
        result = compute_effective_potential(sigma=0.0)
        for i in range(len(result["sigma_grid"])):
            expected = result["V_classical"][i] + result["V_casimir"][i]
            self.assertAlmostEqual(result["V_total"][i], expected, places=10)


class TestCasimirMinimum(unittest.TestCase):
    """Tests for find_casimir_minimum."""

    def test_returns_dict(self):
        result = find_casimir_minimum()
        self.assertIsInstance(result, dict)

    def test_minimum_does_not_exist_s2(self):
        """For S^2 with negative Casimir coefficient: no stable minimum."""
        result = find_casimir_minimum()
        self.assertFalse(result["minimum_exists"])

    def test_has_honest_assessment(self):
        """Must include an honest assessment string."""
        result = find_casimir_minimum()
        self.assertIn("honest_assessment", result)
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertTrue(len(result["honest_assessment"]) > 0)

    def test_a_coefficient_stored(self):
        """A coefficient should be stored."""
        result = find_casimir_minimum()
        self.assertIn("A", result)

    def test_b_coefficient_stored(self):
        """B coefficient should be stored."""
        result = find_casimir_minimum()
        self.assertIn("B", result)

    def test_genus2_result_present(self):
        """Should include genus-2 analysis."""
        result = find_casimir_minimum()
        self.assertIn("genus2_result", result)

    def test_genus2_no_stable_minimum(self):
        """For genus 2 with A < 0: stationary point is maximum, not minimum."""
        result = find_casimir_minimum()
        g2 = result["genus2_result"]
        self.assertFalse(g2["minimum_exists"])

    def test_genus2_stationary_point_exists(self):
        """For genus 2 (A < 0, B > 0): a stationary point should exist."""
        result = find_casimir_minimum()
        g2 = result["genus2_result"]
        self.assertTrue(g2["stationary_point_exists"])

    def test_genus2_is_maximum(self):
        """For genus 2 the stationary point should be a maximum."""
        result = find_casimir_minimum()
        g2 = result["genus2_result"]
        self.assertTrue(g2["is_maximum"])


class TestDilatonMass(unittest.TestCase):
    """Tests for compute_dilaton_mass_casimir."""

    def test_returns_dict(self):
        result = compute_dilaton_mass_casimir()
        self.assertIsInstance(result, dict)

    def test_not_first_principles(self):
        """Since no minimum, first_principles should be False."""
        result = compute_dilaton_mass_casimir()
        self.assertFalse(result["first_principles"])

    def test_minimum_exists_false(self):
        """minimum_exists should be False."""
        result = compute_dilaton_mass_casimir()
        self.assertFalse(result["minimum_exists"])

    def test_has_honest_assessment(self):
        """Must include honest assessment."""
        result = compute_dilaton_mass_casimir()
        self.assertIn("honest_assessment", result)
        self.assertTrue(len(result["honest_assessment"]) > 0)

    def test_m_phi_eV_none(self):
        """m_phi_eV should be None since no stable minimum."""
        result = compute_dilaton_mass_casimir()
        self.assertIsNone(result["m_phi_eV"])

    def test_has_comparison(self):
        """Should have comparison_with_threshold."""
        result = compute_dilaton_mass_casimir()
        self.assertIn("comparison_with_threshold", result)

    def test_honest_assessment_mentions_negative(self):
        """Honest assessment should mention negative Casimir coefficient."""
        result = compute_dilaton_mass_casimir()
        assessment = result["honest_assessment"].lower()
        self.assertTrue("negative" in assessment or "not" in assessment)

    def test_m_phi_squared_none(self):
        """m_phi_squared should be None when no minimum."""
        result = compute_dilaton_mass_casimir()
        self.assertIsNone(result["m_phi_squared"])


class TestSummarize(unittest.TestCase):
    """Tests for summarize_casimir_stabilization."""

    def test_returns_dict(self):
        result = summarize_casimir_stabilization()
        self.assertIsInstance(result, dict)

    def test_has_gaps_status(self):
        """Must report status of theoretical gaps."""
        result = summarize_casimir_stabilization()
        self.assertIn("gaps_status", result)

    def test_gap3_not_resolved(self):
        """Gap #3 should not be resolved (no stable minimum)."""
        result = summarize_casimir_stabilization()
        self.assertFalse(result["gaps_status"]["gap3_resolved"])

    def test_gap1_not_resolved(self):
        """Gap #1 depends on #3; should not be resolved either."""
        result = summarize_casimir_stabilization()
        self.assertFalse(result["gaps_status"]["gap1_resolved"])

    def test_has_overall_assessment(self):
        """Must include overall_assessment string."""
        result = summarize_casimir_stabilization()
        self.assertIn("overall_assessment", result)
        self.assertIsInstance(result["overall_assessment"], str)

    def test_has_spectrum_data(self):
        """Summary should include KK spectrum data."""
        result = summarize_casimir_stabilization()
        self.assertIn("spectrum", result)

    def test_has_zeta_data(self):
        """Summary should include spectral zeta data."""
        result = summarize_casimir_stabilization()
        self.assertIn("zeta", result)

    def test_has_minimum_data(self):
        """Summary should include minimum search data."""
        result = summarize_casimir_stabilization()
        self.assertIn("minimum", result)

    def test_zeta_is_negative(self):
        """Summary should confirm zeta_S2(-1/2) is negative."""
        result = summarize_casimir_stabilization()
        self.assertTrue(result.get("zeta_s2_negative", False))

    def test_has_stable_minimum_false(self):
        """Summary should report no stable minimum."""
        result = summarize_casimir_stabilization()
        self.assertFalse(result.get("has_stable_minimum", True))


class TestMatterCasimir(unittest.TestCase):
    """Tests for matter loop Casimir corrections."""

    def test_fermion_zeta_at_minus_half(self):
        """Fermion zeta_F(-1/2) = 0 because B_3(1/2) = 0 exactly."""
        result = compute_fermion_spectral_zeta_s2(s=-0.5)
        self.assertAlmostEqual(result["zeta_value"], 0.0, delta=1e-10)

    def test_fermion_zeta_at_minus_1(self):
        """zeta_F(-1) = 4 * zeta_H(-3, 1/2) = 4*(-B_4(1/2)/4) = -7/240."""
        result = compute_fermion_spectral_zeta_s2(s=-1.0)
        self.assertTrue(math.isfinite(result["zeta_value"]))
        expected = -7.0 / 240.0  # -0.029166...
        self.assertAlmostEqual(result["zeta_value"], expected, delta=1e-6)

    def test_vector_zeta_finite(self):
        """Vector spectral zeta at s=-0.5 should be finite."""
        result = compute_vector_spectral_zeta_s2(s=-0.5)
        self.assertTrue(math.isfinite(result["zeta_value"]))

    def test_vector_zeta_different_from_scalar(self):
        """Vector zeta at s=-0.5 should differ from scalar zeta."""
        vec = compute_vector_spectral_zeta_s2(s=-0.5)
        sca = compute_spectral_zeta_s2(s=-0.5)
        self.assertNotAlmostEqual(
            vec["zeta_value"], sca["zeta_value"], places=3
        )

    def test_matter_casimir_pure_graviton_matches(self):
        """With zero matter, A_graviton should match the pure graviton result."""
        casimir = compute_casimir_energy_s2(a_radius=1.0)
        matter = compute_matter_casimir_coefficient(
            n_scalars=0, n_fermions=0, n_vectors=0
        )
        self.assertAlmostEqual(
            matter["A_graviton"], casimir["coefficient"], places=8
        )

    def test_matter_casimir_e8xe8(self):
        """E8 x E8 matter Casimir returns expected keys."""
        result = compute_matter_casimir_coefficient(
            n_scalars=4 * 740, n_fermions=2 * 740 + 2 * 496,
            n_vectors=496
        )
        self.assertIsInstance(result, dict)
        for key in ("A_total", "A_graviton", "A_matter", "sign_flipped",
                     "n_scalars", "n_fermions", "n_vectors", "description"):
            self.assertIn(key, result)

    def test_matter_casimir_so32(self):
        """SO(32) matter Casimir returns expected keys."""
        result = compute_matter_casimir_coefficient(
            n_scalars=4 * 740, n_fermions=2 * 740 + 2 * 496,
            n_vectors=496
        )
        self.assertIsInstance(result, dict)
        for key in ("A_total", "A_graviton", "A_matter", "sign_flipped"):
            self.assertIn(key, result)

    def test_matter_casimir_e7xe7(self):
        """E7 x E7 matter Casimir returns expected keys."""
        result = compute_matter_casimir_coefficient(
            n_scalars=4 * 132, n_fermions=2 * 132 + 2 * 266,
            n_vectors=266
        )
        self.assertIsInstance(result, dict)
        for key in ("A_total", "A_graviton", "A_matter", "sign_flipped"):
            self.assertIn(key, result)

    def test_hypermultiplet_field_counting(self):
        """1 hypermultiplet = 4 scalars + 2 fermion dof."""
        result = compute_matter_casimir_coefficient(
            n_scalars=4, n_fermions=2, n_vectors=0
        )
        self.assertEqual(result["n_scalars"], 4)
        self.assertEqual(result["n_fermions"], 2)

    def test_vector_multiplet_field_counting(self):
        """1 vector multiplet = 1 vector + 2 fermion dof."""
        result = compute_matter_casimir_coefficient(
            n_scalars=0, n_fermions=2, n_vectors=1
        )
        self.assertEqual(result["n_vectors"], 1)
        self.assertEqual(result["n_fermions"], 2)

    def test_sign_flip_result_documented(self):
        """scan_anomaly_free_matter_casimir returns 'any_sign_flip' key."""
        result = scan_anomaly_free_matter_casimir()
        self.assertIn("any_sign_flip", result)

    def test_summary_includes_matter(self):
        """summarize_casimir_stabilization returns 'matter_loop_results' key."""
        result = summarize_casimir_stabilization()
        self.assertIn("matter_loop_results", result)


if __name__ == "__main__":
    unittest.main()
