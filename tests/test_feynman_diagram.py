"""
Tests for alpha_ladder_core/feynman_diagram.py

Explicit Feynman diagram calculation for the (1-alpha) correction:
sigma -> KK photon loop -> sigma on S^2.
"""

import math
import unittest

from alpha_ladder_core.feynman_diagram import (
    compute_gauge_kinetic_coupling,
    compute_kk_photon_spectrum,
    compute_passarino_veltman_b0,
    compute_one_loop_self_energy,
    compute_mass_correction,
    verify_diagram_consistency,
    analyze_scheme_dependence,
    summarize_feynman_diagram,
)


class TestComputeGaugeKineticCoupling(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_gauge_kinetic_coupling()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_gauge_kinetic_coupling()
        for key in ("n", "vertex_factor", "gauge_kinetic_function",
                     "coupling_description"):
            self.assertIn(key, result)

    def test_vertex_factor_2_for_s2(self):
        """For S^2 (n=2), vertex factor should be 2."""
        result = compute_gauge_kinetic_coupling(n=2)
        self.assertEqual(result["vertex_factor"], 2)

    def test_vertex_factor_1_for_s1(self):
        """For S^1 (n=1), vertex factor should be 1."""
        result = compute_gauge_kinetic_coupling(n=1)
        self.assertEqual(result["vertex_factor"], 1)

    def test_n_preserved(self):
        result = compute_gauge_kinetic_coupling(n=3)
        self.assertEqual(result["n"], 3)


class TestComputeKKPhotonSpectrum(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_kk_photon_spectrum()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_kk_photon_spectrum()
        for key in ("n", "modes", "total_modes_counted"):
            self.assertIn(key, result)

    def test_first_mode_l1_mass_sq_is_2(self):
        """For S^2, l=1 mode has mass^2/R^2 = 1*(1+1) = 2."""
        result = compute_kk_photon_spectrum(n=2)
        self.assertEqual(result["modes"][0]["l"], 1)
        self.assertEqual(result["modes"][0]["mass_sq_over_R2"], 2)

    def test_l1_degeneracy_is_3(self):
        """For S^2, l=1 mode has degeneracy 2*1+1 = 3."""
        result = compute_kk_photon_spectrum(n=2)
        self.assertEqual(result["modes"][0]["degeneracy"], 3)

    def test_correct_number_of_modes(self):
        """With l_max=5, should have 5 modes (l=1..5)."""
        result = compute_kk_photon_spectrum(n=2, l_max=5)
        self.assertEqual(len(result["modes"]), 5)


class TestPassarinoVeltmanB0(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_passarino_veltman_b0(1.0)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_passarino_veltman_b0(1.0)
        for key in ("m_sq", "mu_R_sq", "b0_finite", "prefactor", "log_term"):
            self.assertIn(key, result)

    def test_prefactor_is_one_over_16pi2(self):
        """Prefactor should be 1/(16*pi^2)."""
        result = compute_passarino_veltman_b0(1.0)
        expected = 1.0 / (16 * math.pi ** 2)
        self.assertAlmostEqual(result["prefactor"], expected, places=15)

    def test_b0_negative_for_m_greater_than_mu(self):
        """B0 should be negative when m^2 > mu_R^2 (log is negative)."""
        result = compute_passarino_veltman_b0(10.0, mu_R_sq=1.0)
        self.assertLess(result["b0_finite"], 0)

    def test_log_term_negative_for_m_greater_than_mu(self):
        """log_term = -ln(m^2/mu_R^2) should be negative when m^2 > mu_R^2."""
        result = compute_passarino_veltman_b0(10.0, mu_R_sq=1.0)
        self.assertLess(result["log_term"], 0)


class TestComputeOneLoopSelfEnergy(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_one_loop_self_energy()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_one_loop_self_energy()
        for key in ("coupling_factor", "loop_factor", "volume_factor",
                     "product", "effective_correction", "mode_sum_info",
                     "assessment"):
            self.assertIn(key, result)

    def test_product_is_one_for_s2(self):
        """Loop factor * volume factor = 1 for S^2."""
        result = compute_one_loop_self_energy(n=2)
        self.assertAlmostEqual(result["product"], 1.0, places=10)

    def test_effective_correction_is_alpha_for_s2(self):
        """Effective correction should be alpha for n=2."""
        from alpha_ladder_core.constants import get_constants
        constants = get_constants("CODATA 2018")
        alpha = float(constants.alpha)
        result = compute_one_loop_self_energy(n=2, constants=constants)
        self.assertAlmostEqual(result["effective_correction"], alpha, places=10)

    def test_coupling_factor_is_n(self):
        """Coupling factor should be n."""
        result = compute_one_loop_self_energy(n=2)
        self.assertEqual(result["coupling_factor"], 2)

    def test_volume_cancellation_occurs(self):
        """1/(4*pi) * 4*pi = 1 for S^2."""
        result = compute_one_loop_self_energy(n=2)
        loop = result["loop_factor"]
        vol = result["volume_factor"]
        self.assertAlmostEqual(loop * vol, 1.0, places=10)

    def test_assessment_nonempty(self):
        result = compute_one_loop_self_energy()
        self.assertGreater(len(result["assessment"]), 0)


class TestComputeMassCorrection(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_mass_correction()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_mass_correction()
        for key in ("correction_to_mass", "correction_to_ratio", "effective_k",
                     "k_bare", "k_corrected"):
            self.assertIn(key, result)

    def test_correction_is_negative(self):
        """Mass correction should be negative (self-energy reduces mass)."""
        result = compute_mass_correction()
        self.assertLess(result["correction_to_mass"], 0)

    def test_k_corrected_less_than_k_bare(self):
        """k_corrected = sqrt(phi)*(1-alpha) < sqrt(phi) = k_bare."""
        result = compute_mass_correction()
        self.assertLess(result["k_corrected"], result["k_bare"])

    def test_correction_magnitude_is_alpha(self):
        """The magnitude of the correction should be alpha."""
        from alpha_ladder_core.constants import get_constants
        constants = get_constants("CODATA 2018")
        alpha = float(constants.alpha)
        result = compute_mass_correction(constants=constants)
        self.assertAlmostEqual(abs(result["correction_to_mass"]), alpha, places=10)

    def test_effective_k_equals_sqrt_phi_times_one_minus_alpha(self):
        """effective_k should be approximately sqrt(phi)*(1-alpha)."""
        from alpha_ladder_core.constants import get_constants
        constants = get_constants("CODATA 2018")
        alpha = float(constants.alpha)
        phi = constants.phi
        sqrt_phi = float(phi.sqrt())
        expected = sqrt_phi * (1 - alpha)
        result = compute_mass_correction(constants=constants)
        self.assertAlmostEqual(result["effective_k"], expected, places=10)


class TestVerifyDiagramConsistency(unittest.TestCase):

    def test_returns_dict(self):
        result = verify_diagram_consistency()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = verify_diagram_consistency()
        for key in ("diagram_correction", "volume_cancellation_correction",
                     "consistent", "G_predicted", "residual_ppm"):
            self.assertIn(key, result)

    def test_consistent_true(self):
        """Diagram correction should match volume cancellation correction."""
        result = verify_diagram_consistency()
        self.assertTrue(result["consistent"])

    def test_residual_sub_ppm(self):
        """Residual should be less than 1 ppm."""
        result = verify_diagram_consistency()
        self.assertLess(abs(result["residual_ppm"]), 1.0)

    def test_G_in_reasonable_range(self):
        """G_predicted should be in [6.6e-11, 6.8e-11]."""
        result = verify_diagram_consistency()
        G = float(result["G_predicted"])
        self.assertGreater(G, 6.6e-11)
        self.assertLess(G, 6.8e-11)


class TestAnalyzeSchemeDependence(unittest.TestCase):

    def test_returns_dict(self):
        result = analyze_scheme_dependence()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = analyze_scheme_dependence()
        for key in ("scheme_independent_part", "scheme_dependent_part",
                     "conclusion"):
            self.assertIn(key, result)

    def test_scheme_independent_describes_coefficient(self):
        """scheme_independent_part should mention the coefficient."""
        result = analyze_scheme_dependence()
        self.assertIn("coefficient", result["scheme_independent_part"].lower())

    def test_conclusion_nonempty(self):
        result = analyze_scheme_dependence()
        self.assertGreater(len(result["conclusion"]), 0)


class TestSummarizeFeynmanDiagram(unittest.TestCase):

    def test_returns_dict(self):
        result = summarize_feynman_diagram()
        self.assertIsInstance(result, dict)

    def test_all_sub_keys_present(self):
        result = summarize_feynman_diagram()
        for key in ("gauge_coupling", "kk_spectrum", "self_energy",
                     "mass_correction", "consistency", "scheme_analysis",
                     "diagram_description", "honest_assessment"):
            self.assertIn(key, result)

    def test_diagram_description_nonempty(self):
        result = summarize_feynman_diagram()
        self.assertIsInstance(result["diagram_description"], str)
        self.assertGreater(len(result["diagram_description"]), 0)

    def test_honest_assessment_nonempty(self):
        result = summarize_feynman_diagram()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 0)

    def test_sub_results_are_dicts(self):
        result = summarize_feynman_diagram()
        for key in ("gauge_coupling", "kk_spectrum", "self_energy",
                     "mass_correction", "consistency", "scheme_analysis"):
            self.assertIsInstance(result[key], dict)

    def test_diagram_mentions_sigma(self):
        """The diagram description should mention the sigma field."""
        result = summarize_feynman_diagram()
        self.assertIn("sigma", result["diagram_description"])


if __name__ == "__main__":
    unittest.main()
