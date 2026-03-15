"""
Tests for alpha_ladder_core/one_alpha_derivation.py

Derivation of the (1-alpha) correction from S^2 volume cancellation.
"""

import math
import unittest

from alpha_ladder_core.one_alpha_derivation import (
    compute_loop_factor,
    compute_s2_volume_factor,
    analyze_volume_cancellation,
    scan_candidate_mechanisms,
    analyze_uniqueness_of_s2,
    derive_sign,
    derive_degeneracy,
    compute_mode_sum,
    compute_corrected_prediction,
    summarize_one_alpha_derivation,
)


class TestComputeLoopFactor(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_loop_factor()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_loop_factor()
        for key in ("d", "loop_factor_full", "mass_correction_factor", "formula_description"):
            self.assertIn(key, result)

    def test_d4_mass_correction_factor(self):
        """In d=4, mass correction factor should be 1/(4*pi)."""
        result = compute_loop_factor(d=4)
        expected = 1.0 / (4 * math.pi)
        self.assertAlmostEqual(result["mass_correction_factor"], expected, places=10)

    def test_d_preserved(self):
        result = compute_loop_factor(d=6)
        self.assertEqual(result["d"], 6)

    def test_description_nonempty(self):
        result = compute_loop_factor()
        self.assertGreater(len(result["formula_description"]), 0)


class TestComputeS2VolumeFactor(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_s2_volume_factor()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_s2_volume_factor()
        for key in ("n", "volume_over_R_n", "volume_formula", "cancellation_factor"):
            self.assertIn(key, result)

    def test_s2_volume(self):
        """Vol(S^2)/R^2 = 4*pi."""
        result = compute_s2_volume_factor(n=2)
        self.assertAlmostEqual(result["volume_over_R_n"], 4 * math.pi, places=10)

    def test_s1_volume(self):
        """Vol(S^1)/R = 2*pi."""
        result = compute_s2_volume_factor(n=1)
        self.assertAlmostEqual(result["volume_over_R_n"], 2 * math.pi, places=10)

    def test_s3_volume(self):
        """Vol(S^3)/R^3 = 2*pi^2."""
        result = compute_s2_volume_factor(n=3)
        self.assertAlmostEqual(result["volume_over_R_n"], 2 * math.pi**2, places=10)

    def test_n_preserved(self):
        result = compute_s2_volume_factor(n=4)
        self.assertEqual(result["n"], 4)


class TestAnalyzeVolumeCancellation(unittest.TestCase):

    def test_returns_dict(self):
        result = analyze_volume_cancellation()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = analyze_volume_cancellation()
        for key in ("d", "n", "loop_factor", "volume_factor", "product",
                     "effective_correction", "ratio_to_alpha", "matches_alpha",
                     "mechanism", "assessment"):
            self.assertIn(key, result)

    def test_product_is_one_for_s2(self):
        """The key result: loop_factor * volume_factor = 1 for S^2 in d=4."""
        result = analyze_volume_cancellation(d=4, n=2)
        self.assertAlmostEqual(result["product"], 1.0, places=10)

    def test_matches_alpha_true_for_s2(self):
        result = analyze_volume_cancellation(d=4, n=2)
        self.assertTrue(result["matches_alpha"])

    def test_ratio_to_alpha_is_one(self):
        result = analyze_volume_cancellation(d=4, n=2)
        self.assertAlmostEqual(result["ratio_to_alpha"], 1.0, places=10)

    def test_no_cancellation_for_s3(self):
        result = analyze_volume_cancellation(d=4, n=3)
        self.assertFalse(result["matches_alpha"])

    def test_assessment_nonempty(self):
        result = analyze_volume_cancellation()
        self.assertGreater(len(result["assessment"]), 0)


class TestScanCandidateMechanisms(unittest.TestCase):

    def test_returns_dict(self):
        result = scan_candidate_mechanisms()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = scan_candidate_mechanisms()
        for key in ("mechanisms", "best_match", "n_matches", "assessment"):
            self.assertIn(key, result)

    def test_multiple_mechanisms(self):
        result = scan_candidate_mechanisms()
        self.assertGreater(len(result["mechanisms"]), 5)

    def test_exactly_one_match(self):
        """Only one mechanism should produce bare alpha."""
        result = scan_candidate_mechanisms()
        self.assertEqual(result["n_matches"], 1)

    def test_best_match_is_kk_graviton(self):
        result = scan_candidate_mechanisms()
        self.assertIsNotNone(result["best_match"])
        self.assertIn("S^2", result["best_match"]["name"])

    def test_qed_vertex_does_not_match(self):
        result = scan_candidate_mechanisms()
        qed = [m for m in result["mechanisms"] if "QED vertex" in m["name"]]
        self.assertTrue(len(qed) > 0)
        self.assertFalse(qed[0]["matches_alpha"])

    def test_assessment_nonempty(self):
        result = scan_candidate_mechanisms()
        self.assertGreater(len(result["assessment"]), 0)


class TestAnalyzeUniquenessOfS2(unittest.TestCase):

    def test_returns_dict(self):
        result = analyze_uniqueness_of_s2()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = analyze_uniqueness_of_s2()
        for key in ("scan_results", "unique_n", "s2_is_unique", "assessment"):
            self.assertIn(key, result)

    def test_six_spheres_scanned(self):
        result = analyze_uniqueness_of_s2()
        self.assertEqual(len(result["scan_results"]), 6)

    def test_unique_n_is_2(self):
        """S^2 should be the unique sphere with volume cancellation."""
        result = analyze_uniqueness_of_s2()
        self.assertEqual(result["unique_n"], 2)

    def test_s2_is_unique_true(self):
        result = analyze_uniqueness_of_s2()
        self.assertTrue(result["s2_is_unique"])

    def test_only_s2_gives_bare_alpha(self):
        result = analyze_uniqueness_of_s2()
        matches = [r for r in result["scan_results"] if r["gives_bare_alpha"]]
        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0]["n"], 2)


class TestDeriveSign(unittest.TestCase):
    """Tests for derive_sign."""

    def test_returns_dict(self):
        result = derive_sign()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = derive_sign()
        for key in ("sign", "correction_form", "physical_argument",
                     "numerical_cross_check", "assessment"):
            self.assertIn(key, result)

    def test_sign_is_negative(self):
        """The sign should be -1 (mass reduction)."""
        result = derive_sign()
        self.assertEqual(result["sign"], -1)

    def test_correction_form(self):
        result = derive_sign()
        self.assertEqual(result["correction_form"], "(1 - alpha)")

    def test_minus_moves_toward_data(self):
        """(1-alpha) should move G toward measured value."""
        result = derive_sign()
        self.assertTrue(result["numerical_cross_check"]["minus_moves_toward_data"])

    def test_assessment_nonempty(self):
        result = derive_sign()
        self.assertGreater(len(result["assessment"]), 0)


class TestDeriveDegeneracy(unittest.TestCase):
    """Tests for derive_degeneracy."""

    def test_returns_dict(self):
        result = derive_degeneracy()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = derive_degeneracy()
        for key in ("kk_spectrum", "l0_degeneracy", "l1_degeneracy",
                     "physical_argument", "why_not_3", "assessment"):
            self.assertIn(key, result)

    def test_l0_degeneracy_is_1(self):
        """l=0 sector should have degeneracy 1."""
        result = derive_degeneracy()
        self.assertEqual(result["l0_degeneracy"], 1)

    def test_l1_degeneracy_is_3(self):
        """l=1 sector should have degeneracy 3."""
        result = derive_degeneracy()
        self.assertEqual(result["l1_degeneracy"], 3)

    def test_spectrum_has_entries(self):
        result = derive_degeneracy()
        self.assertGreater(len(result["kk_spectrum"]), 3)

    def test_why_not_3_nonempty(self):
        result = derive_degeneracy()
        self.assertGreater(len(result["why_not_3"]), 0)


class TestComputeModeSum(unittest.TestCase):
    """Tests for compute_mode_sum."""

    def test_returns_dict(self):
        result = compute_mode_sum(l_max=100)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_mode_sum(l_max=100)
        for key in ("spectral_zeta_minus1", "polynomial_identity_verified",
                     "gravitational_vs_em", "assessment"):
            self.assertIn(key, result)

    def test_spectral_zeta_minus1_is_zero(self):
        """Zeta-regularized sum should be exactly zero."""
        result = compute_mode_sum(l_max=100)
        self.assertEqual(result["spectral_zeta_minus1"], 0)

    def test_polynomial_identity_verified(self):
        result = compute_mode_sum(l_max=100)
        self.assertTrue(result["polynomial_identity_verified"])

    def test_gravitational_too_small(self):
        """Gravitational corrections should be ~43 orders below EM."""
        result = compute_mode_sum(l_max=100)
        gap = result["gravitational_vs_em"]["orders_of_magnitude_gap"]
        self.assertGreater(gap, 40)

    def test_assessment_nonempty(self):
        result = compute_mode_sum(l_max=100)
        self.assertGreater(len(result["assessment"]), 0)


class TestComputeCorrectedPrediction(unittest.TestCase):

    def test_returns_dict(self):
        result = compute_corrected_prediction()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_corrected_prediction()
        for key in ("G_predicted", "G_measured", "residual_ppm",
                     "k_offset", "sqrt_phi", "derivation_chain", "assessment"):
            self.assertIn(key, result)

    def test_residual_sub_ppm(self):
        result = compute_corrected_prediction()
        self.assertLess(abs(result["residual_ppm"]), 1.0)

    def test_derivation_chain_has_steps(self):
        result = compute_corrected_prediction()
        self.assertGreater(len(result["derivation_chain"]), 5)

    def test_G_predicted_reasonable(self):
        result = compute_corrected_prediction()
        self.assertGreater(float(result["G_predicted"]), 6.6e-11)
        self.assertLess(float(result["G_predicted"]), 6.8e-11)


class TestSummarizeOneAlphaDerivation(unittest.TestCase):

    def test_returns_dict(self):
        result = summarize_one_alpha_derivation()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = summarize_one_alpha_derivation()
        for key in ("volume_cancellation", "mechanisms", "uniqueness",
                     "sign_analysis", "degeneracy_analysis", "mode_sum",
                     "prediction", "key_finding", "honest_assessment"):
            self.assertIn(key, result)

    def test_key_finding_nonempty(self):
        result = summarize_one_alpha_derivation()
        self.assertIsInstance(result["key_finding"], str)
        self.assertGreater(len(result["key_finding"]), 0)

    def test_honest_assessment_nonempty(self):
        result = summarize_one_alpha_derivation()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 0)

    def test_sub_results_are_dicts(self):
        result = summarize_one_alpha_derivation()
        for key in ("volume_cancellation", "mechanisms", "uniqueness",
                     "sign_analysis", "degeneracy_analysis", "mode_sum",
                     "prediction"):
            self.assertIsInstance(result[key], dict)


if __name__ == "__main__":
    unittest.main()
