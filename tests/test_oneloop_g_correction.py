"""
Tests for alpha_ladder_core/oneloop_g_correction.py

Comprehensive test suite for one-loop quantum corrections to Newton's constant G
from Kaluza-Klein towers on S^2. The key physics result verified here is that
the spectral zeta function at s=-1 is zero (or very close), because the
weighted sum sum_{l=1}^{L} (2l+1)*l*(l+1) is an exact polynomial in L with no
constant term. This means naive zeta regularization does NOT produce the
desired 3*alpha^2 correction to G.
"""

import math
import unittest

from alpha_ladder_core.oneloop_g_correction import (
    kk_spectrum_s2,
    compute_oneloop_correction,
    scan_spin_coefficients,
    analyze_so3_contribution,
    compute_spectral_zeta_s2,
    summarize_oneloop_calculation,
)


# Fine-structure constant for reference
ALPHA = 0.0072973525693


class TestKKSpectrumS2(unittest.TestCase):
    """Tests for kk_spectrum_s2."""

    def test_returns_dict(self):
        result = kk_spectrum_s2(l_max=5)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = kk_spectrum_s2(l_max=5)
        for key in ("levels", "sectors", "total_dof_per_level"):
            self.assertIn(key, result, f"Missing key: {key}")

    def test_levels_length(self):
        """Number of levels should equal l_max."""
        for l_max in (5, 10, 20):
            result = kk_spectrum_s2(l_max=l_max)
            self.assertEqual(len(result["levels"]), l_max)

    def test_l1_mass_and_degeneracy(self):
        """l=1: mass_sq/R^2 = l(l+1) = 2, degeneracy = 2l+1 = 3."""
        result = kk_spectrum_s2(l_max=5)
        level_1 = result["levels"][0]
        self.assertAlmostEqual(level_1["mass_sq_over_R2"], 2.0, places=10)
        self.assertEqual(level_1["degeneracy"], 3)
        self.assertEqual(level_1["l"], 1)

    def test_l2_mass_and_degeneracy(self):
        """l=2: mass_sq/R^2 = l(l+1) = 6, degeneracy = 2l+1 = 5."""
        result = kk_spectrum_s2(l_max=5)
        level_2 = result["levels"][1]
        self.assertAlmostEqual(level_2["mass_sq_over_R2"], 6.0, places=10)
        self.assertEqual(level_2["degeneracy"], 5)
        self.assertEqual(level_2["l"], 2)

    def test_sectors_keys(self):
        """sectors should have scalars, vectors, tensors sub-dicts."""
        result = kk_spectrum_s2(l_max=5)
        sectors = result["sectors"]
        for key in ("scalars", "vectors", "tensors"):
            self.assertIn(key, sectors, f"Missing sector: {key}")

    def test_sector_field_counts(self):
        """N_scalar = 3, N_vector = 3, N_tensor = 1."""
        result = kk_spectrum_s2(l_max=5)
        sectors = result["sectors"]
        self.assertEqual(sectors["scalars"]["n_fields"], 3)
        self.assertEqual(sectors["vectors"]["n_fields"], 3)
        self.assertEqual(sectors["tensors"]["n_fields"], 1)

    def test_total_dof_positive(self):
        result = kk_spectrum_s2(l_max=5)
        self.assertGreater(result["total_dof_per_level"], 0)


class TestComputeSpectralZetaS2(unittest.TestCase):
    """Tests for compute_spectral_zeta_s2."""

    def test_returns_dict(self):
        result = compute_spectral_zeta_s2(2.0, l_max=100)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_spectral_zeta_s2(2.0, l_max=100)
        for key in ("s", "value", "method", "l_max_used", "converged"):
            self.assertIn(key, result, f"Missing key: {key}")

    def test_s2_convergent_positive(self):
        """zeta(2) on S^2 should be a finite positive number."""
        result = compute_spectral_zeta_s2(2.0, l_max=1000)
        self.assertGreater(result["value"], 0)
        self.assertTrue(math.isfinite(result["value"]))

    def test_s3_smaller_than_s2(self):
        """zeta(3) < zeta(2) since higher s damps more."""
        r2 = compute_spectral_zeta_s2(2.0, l_max=1000)
        r3 = compute_spectral_zeta_s2(3.0, l_max=1000)
        self.assertGreater(r2["value"], r3["value"])

    def test_converged_for_large_s(self):
        """Should have converged key for Re(s) > 1."""
        for s in (2.0, 3.0, 5.0):
            result = compute_spectral_zeta_s2(s, l_max=1000)
            self.assertIn("converged", result)

    def test_method_direct_for_convergent(self):
        """For s > 1 should use direct summation."""
        result = compute_spectral_zeta_s2(2.0, l_max=100)
        self.assertEqual(result["method"], "direct")

    def test_s_minus_1_close_to_zero(self):
        """
        The key physics result: zeta(-1) on S^2 should be 0 or very close.

        The weighted sum sum_{l=1}^{L} (2l+1)*[l(l+1)]^s at s=-1 is
        sum (2l+1)/(l(l+1)) which telescopes... BUT the full sum
        sum (2l+1)*l*(l+1) = L(L+1)^2(L+2)/2 is an exact polynomial,
        so its analytic continuation (zeta regularization) at s=-1
        should give approximately zero.
        """
        result = compute_spectral_zeta_s2(-1.0, l_max=10000)
        self.assertAlmostEqual(result["value"], 0.0, delta=1e-6)

    def test_l_max_used_matches_input(self):
        result = compute_spectral_zeta_s2(2.0, l_max=500)
        self.assertEqual(result["l_max_used"], 500)


class TestComputeOneloopCorrection(unittest.TestCase):
    """Tests for compute_oneloop_correction."""

    def test_returns_dict(self):
        result = compute_oneloop_correction(1.0, l_max=50)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_oneloop_correction(1.0, l_max=50)
        for key in ("R", "l_max", "raw_sum", "subtracted_divergences",
                     "finite_part", "sector_contributions",
                     "total_delta_G_over_G"):
            self.assertIn(key, result, f"Missing key: {key}")

    def test_R_preserved(self):
        result = compute_oneloop_correction(2.5, l_max=50)
        self.assertAlmostEqual(result["R"], 2.5)

    def test_finite_part_not_none(self):
        result = compute_oneloop_correction(1.0, l_max=50)
        self.assertIsNotNone(result["finite_part"])

    def test_sector_contributions_keys(self):
        """sector_contributions should have scalar, vector, tensor."""
        result = compute_oneloop_correction(1.0, l_max=50)
        sc = result["sector_contributions"]
        self.assertIsInstance(sc, dict)
        # At minimum we expect keys for the three sectors
        self.assertTrue(len(sc) >= 3,
                        "Expected at least 3 sector contributions")

    def test_different_R_values(self):
        """Should work for different R without error."""
        for R in (0.1, 1.0, 10.0):
            result = compute_oneloop_correction(R, l_max=20)
            self.assertIsInstance(result["total_delta_G_over_G"], (int, float))


class TestScanSpinCoefficients(unittest.TestCase):
    """Tests for scan_spin_coefficients."""

    def test_returns_dict(self):
        result = scan_spin_coefficients()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = scan_spin_coefficients()
        for key in ("target_correction", "N_scalar", "N_vector", "N_tensor",
                     "best_combinations", "pure_vector_match", "assessment"):
            self.assertIn(key, result, f"Missing key: {key}")

    def test_field_counts(self):
        """N_scalar = 3, N_vector = 3, N_tensor = 1."""
        result = scan_spin_coefficients()
        self.assertEqual(result["N_scalar"], 3)
        self.assertEqual(result["N_vector"], 3)
        self.assertEqual(result["N_tensor"], 1)

    def test_target_correction_approx_3_alpha_sq(self):
        """target should be approximately 3 * alpha^2 ~ 1.6e-4."""
        result = scan_spin_coefficients()
        expected = 3.0 * ALPHA ** 2
        self.assertAlmostEqual(result["target_correction"], expected,
                               delta=1e-6)

    def test_best_combinations_is_list(self):
        result = scan_spin_coefficients()
        self.assertIsInstance(result["best_combinations"], list)

    def test_assessment_nonempty(self):
        result = scan_spin_coefficients()
        self.assertIsInstance(result["assessment"], str)
        self.assertGreater(len(result["assessment"]), 0)


class TestAnalyzeSO3Contribution(unittest.TestCase):
    """Tests for analyze_so3_contribution."""

    def test_returns_dict(self):
        result = analyze_so3_contribution()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = analyze_so3_contribution()
        for key in ("n_vector_fields", "vector_coefficient",
                     "spectral_zeta_minus_1", "delta_G_over_G_from_vectors",
                     "matches_3_alpha_sq", "assessment"):
            self.assertIn(key, result, f"Missing key: {key}")

    def test_n_vector_fields_is_3(self):
        """Should have 3 vector fields matching the 3 in 3*alpha^2."""
        result = analyze_so3_contribution()
        self.assertEqual(result["n_vector_fields"], 3)

    def test_spectral_zeta_minus_1_close_to_zero(self):
        """Spectral zeta at s=-1 should be close to 0."""
        result = analyze_so3_contribution()
        self.assertAlmostEqual(result["spectral_zeta_minus_1"], 0.0,
                               delta=1e-6)

    def test_vector_coefficient_present(self):
        result = analyze_so3_contribution()
        self.assertIn("vector_coefficient", result)
        self.assertIsNotNone(result["vector_coefficient"])

    def test_assessment_nonempty(self):
        result = analyze_so3_contribution()
        self.assertIsInstance(result["assessment"], str)
        self.assertGreater(len(result["assessment"]), 0)


class TestSummarizeOneloop(unittest.TestCase):
    """Tests for summarize_oneloop_calculation."""

    def test_returns_dict(self):
        result = summarize_oneloop_calculation()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = summarize_oneloop_calculation()
        for key in ("spectrum", "spectral_zeta", "so3_contribution",
                     "spin_scan", "total_correction", "matches_3_alpha_sq",
                     "honest_assessment"):
            self.assertIn(key, result, f"Missing key: {key}")

    def test_honest_assessment_nonempty(self):
        """Should contain a non-empty honest assessment string."""
        result = summarize_oneloop_calculation()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 0)

    def test_spectrum_is_dict(self):
        result = summarize_oneloop_calculation()
        self.assertIsInstance(result["spectrum"], dict)

    def test_spectral_zeta_is_dict(self):
        result = summarize_oneloop_calculation()
        self.assertIsInstance(result["spectral_zeta"], dict)

    def test_so3_n_vector_fields(self):
        """so3_contribution should report n_vector_fields = 3."""
        result = summarize_oneloop_calculation()
        so3 = result["so3_contribution"]
        self.assertEqual(so3["n_vector_fields"], 3)

    def test_total_correction_is_number(self):
        result = summarize_oneloop_calculation()
        self.assertIsInstance(result["total_correction"], (int, float))


class TestPolynomialSumFormula(unittest.TestCase):
    """
    Verify the exact polynomial identity that drives the key physics result.

    sum_{l=1}^{L} (2l+1)*l*(l+1) = L*(L+1)^2*(L+2)/2

    This is a degree-4 polynomial in L with no constant term, which means
    the zeta-regularized value (analytic continuation to L -> infinity)
    should give zero, making naive one-loop zeta regularization unable
    to produce a nonzero 3*alpha^2 correction.
    """

    def test_polynomial_identity_small(self):
        """Verify the polynomial identity for small L values."""
        for L in range(1, 20):
            direct = sum((2 * l + 1) * l * (l + 1) for l in range(1, L + 1))
            formula = L * (L + 1) ** 2 * (L + 2) // 2
            self.assertEqual(direct, formula,
                             f"Polynomial identity fails at L={L}")

    def test_polynomial_no_constant_term(self):
        """L=0 gives 0, confirming no constant term."""
        L = 0
        formula = L * (L + 1) ** 2 * (L + 2) // 2
        self.assertEqual(formula, 0)

    def test_polynomial_leading_coefficient(self):
        """Leading term is L^4/2 for large L."""
        L = 10000
        formula = L * (L + 1) ** 2 * (L + 2) / 2
        leading = L ** 4 / 2
        ratio = formula / leading
        self.assertAlmostEqual(ratio, 1.0, delta=1e-3)


if __name__ == "__main__":
    unittest.main()
