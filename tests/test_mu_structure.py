"""Tests for alpha_ladder_core.mu_structure module.

Comprehensive tests for the mu-structure analysis: the formula
alpha^24 * mu * (mu - sqrt(phi)) and its comparison to the hierarchy
formula alpha^24 * mu^2.
"""

import unittest
from decimal import Decimal

from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.mu_structure import (
    compute_mu_structure,
    compute_mu_structure_refined,
    scan_mu_offsets,
    analyze_sqrt_phi_origin,
    compare_all_bridges,
    analyze_codata_stability_mu,
    summarize_mu_structure,
)

ALPHA = 0.0072973525693


class TestComputeMuStructure(unittest.TestCase):
    """Tests for compute_mu_structure."""

    @classmethod
    def setUpClass(cls):
        cls.constants = get_constants("CODATA 2018")
        cls.result = compute_mu_structure(cls.constants)

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        expected = {
            "mu", "sqrt_phi", "mu_minus_sqrt_phi",
            "alpha_g_mu_sq", "alpha_g_mu_structure", "alpha_g_measured",
            "G_predicted", "G_measured",
            "residual_ppm", "hierarchy_residual_ppm", "improvement_factor",
            "within_codata_uncertainty", "sqrt_phi_m_e_MeV",
        }
        self.assertTrue(expected.issubset(self.result.keys()),
                        f"Missing keys: {expected - self.result.keys()}")

    def test_mu_value(self):
        mu = float(self.result["mu"])
        self.assertAlmostEqual(mu, 1836.15, delta=0.01)

    def test_sqrt_phi_value(self):
        sqrt_phi = float(self.result["sqrt_phi"])
        self.assertAlmostEqual(sqrt_phi, 1.2720, delta=0.001)

    def test_residual_ppm_within_codata(self):
        self.assertLess(abs(self.result["residual_ppm"]), 22.0)

    def test_residual_ppm_sign(self):
        # Predicted < measured, so residual should be negative
        self.assertLess(self.result["residual_ppm"], 0)

    def test_residual_magnitude(self):
        abs_res = abs(self.result["residual_ppm"])
        self.assertGreater(abs_res, 1.0)
        self.assertLess(abs_res, 10.0)

    def test_hierarchy_residual_larger(self):
        self.assertGreater(
            abs(self.result["hierarchy_residual_ppm"]),
            abs(self.result["residual_ppm"]),
        )

    def test_improvement_factor(self):
        self.assertGreater(self.result["improvement_factor"], 50)

    def test_sqrt_phi_m_e_MeV(self):
        val = self.result["sqrt_phi_m_e_MeV"]
        self.assertGreater(val, 0.6)
        self.assertLess(val, 0.7)


class TestComputeMuStructureRefined(unittest.TestCase):
    """Tests for compute_mu_structure_refined -- the (1-alpha) corrected formula."""

    def test_returns_dict(self):
        result = compute_mu_structure_refined()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        result = compute_mu_structure_refined()
        for key in ("alpha_g", "G_predicted", "mu", "sqrt_phi", "k_offset",
                     "residual_ppm", "bare_residual_ppm", "improvement_over_bare",
                     "assessment"):
            self.assertIn(key, result)

    def test_residual_sub_ppm(self):
        """Refined formula should achieve sub-ppm residual."""
        result = compute_mu_structure_refined()
        self.assertLess(abs(result["residual_ppm"]), 1.0)

    def test_residual_approx_minus_0_31(self):
        """Residual should be approximately -0.31 ppm."""
        result = compute_mu_structure_refined()
        self.assertAlmostEqual(result["residual_ppm"], -0.31, delta=0.1)

    def test_improvement_over_bare(self):
        """Should improve over bare sqrt(phi) by at least 10x."""
        result = compute_mu_structure_refined()
        self.assertGreater(result["improvement_over_bare"], 10)

    def test_k_offset_less_than_sqrt_phi(self):
        """k_offset = sqrt(phi)*(1-alpha) < sqrt(phi)."""
        result = compute_mu_structure_refined()
        self.assertLess(result["k_offset"], result["sqrt_phi"])

    def test_k_offset_close_to_sqrt_phi(self):
        """k_offset should be within 1% of sqrt(phi)."""
        result = compute_mu_structure_refined()
        ratio = result["k_offset"] / result["sqrt_phi"]
        self.assertAlmostEqual(ratio, 1.0, delta=0.01)

    def test_G_predicted_reasonable(self):
        """G should be in the right ballpark."""
        result = compute_mu_structure_refined()
        self.assertGreater(float(result["G_predicted"]), 6.6e-11)
        self.assertLess(float(result["G_predicted"]), 6.8e-11)

    def test_bare_residual_larger(self):
        """Bare residual should be larger in magnitude than refined."""
        result = compute_mu_structure_refined()
        self.assertGreater(abs(result["bare_residual_ppm"]), abs(result["residual_ppm"]))

    def test_assessment_nonempty(self):
        result = compute_mu_structure_refined()
        self.assertIsInstance(result["assessment"], str)
        self.assertGreater(len(result["assessment"]), 0)

    def test_zero_fitted_params(self):
        """The assessment should mention zero fitted parameters."""
        result = compute_mu_structure_refined()
        self.assertIn("Zero fitted", result["assessment"])


class TestScanMuOffsets(unittest.TestCase):
    """Tests for scan_mu_offsets."""

    @classmethod
    def setUpClass(cls):
        cls.constants = get_constants("CODATA 2018")
        cls.result = scan_mu_offsets(cls.constants)

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        expected = {"offsets", "best_offset", "assessment"}
        self.assertTrue(expected.issubset(self.result.keys()))

    def test_offsets_is_list(self):
        offsets = self.result["offsets"]
        self.assertIsInstance(offsets, list)
        self.assertGreaterEqual(len(offsets), 8)

    def test_k_zero_is_mu_sq(self):
        """k=0 gives the original mu^2 formula with ~688 ppm."""
        offsets = self.result["offsets"]
        k0 = [o for o in offsets if float(o["k_value"]) == 0.0]
        self.assertTrue(len(k0) > 0, "k=0 offset not found")
        residual = abs(k0[0]["residual_ppm"])
        self.assertGreater(residual, 500)
        self.assertLess(residual, 900)

    def test_k_one_is_mu_minus_one(self):
        """k=1 should give a residual around 100-200 ppm."""
        offsets = self.result["offsets"]
        k1 = [o for o in offsets if abs(float(o["k_value"]) - 1.0) < 0.01]
        self.assertTrue(len(k1) > 0, "k=1 offset not found")
        residual = abs(k1[0]["residual_ppm"])
        self.assertGreater(residual, 50)
        self.assertLess(residual, 300)

    def test_sqrt_phi_best(self):
        """best_offset should be sqrt(phi) or very close."""
        best = self.result["best_offset"]
        k_val = float(best["k_value"])
        sqrt_phi = (1 + 5**0.5) / 2
        sqrt_phi = sqrt_phi**0.5
        self.assertAlmostEqual(k_val, sqrt_phi, delta=0.01)

    def test_best_within_codata(self):
        best = self.result["best_offset"]
        self.assertLess(abs(best["residual_ppm"]), 22.0)

    def test_assessment_nonempty(self):
        self.assertIsInstance(self.result["assessment"], str)
        self.assertGreater(len(self.result["assessment"]), 20)


class TestAnalyzeSqrtPhiOrigin(unittest.TestCase):
    """Tests for analyze_sqrt_phi_origin."""

    @classmethod
    def setUpClass(cls):
        cls.constants = get_constants("CODATA 2018")
        cls.result = analyze_sqrt_phi_origin(cls.constants)

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        expected = {
            "vacuum_polynomial_roots", "phi_from_roots", "sqrt_phi_value",
            "sqrt_phi_m_e_MeV", "subtracted_mass_MeV", "proton_mass_MeV",
            "fraction_subtracted", "physical_interpretations", "assessment",
        }
        self.assertTrue(expected.issubset(self.result.keys()),
                        f"Missing keys: {expected - self.result.keys()}")

    def test_vacuum_polynomial_roots(self):
        """Two roots of x^2 + 6x + 4 = 0: -3 +/- sqrt(5)."""
        roots = self.result["vacuum_polynomial_roots"]
        self.assertEqual(len(roots), 2)
        # -3 + sqrt(5) ~ -0.764, -3 - sqrt(5) ~ -5.236
        roots_sorted = sorted(roots)
        self.assertAlmostEqual(roots_sorted[0], -5.236, delta=0.01)
        self.assertAlmostEqual(roots_sorted[1], -0.764, delta=0.01)

    def test_phi_from_roots(self):
        text = self.result["phi_from_roots"]
        self.assertIsInstance(text, str)
        self.assertGreater(len(text), 10)

    def test_sqrt_phi_value(self):
        self.assertAlmostEqual(self.result["sqrt_phi_value"], 1.2720, delta=0.001)

    def test_subtracted_mass_positive(self):
        self.assertGreater(self.result["subtracted_mass_MeV"], 0)

    def test_physical_interpretations_list(self):
        interps = self.result["physical_interpretations"]
        self.assertIsInstance(interps, list)
        self.assertGreaterEqual(len(interps), 3)

    def test_assessment_nonempty(self):
        self.assertIsInstance(self.result["assessment"], str)
        self.assertGreater(len(self.result["assessment"]), 20)


class TestCompareAllBridges(unittest.TestCase):
    """Tests for compare_all_bridges."""

    @classmethod
    def setUpClass(cls):
        cls.constants = get_constants("CODATA 2018")
        cls.result = compare_all_bridges(cls.constants)

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        expected = {"formulas", "best_zero_param", "best_overall", "assessment"}
        self.assertTrue(expected.issubset(self.result.keys()))

    def test_formulas_list_length(self):
        formulas = self.result["formulas"]
        self.assertIsInstance(formulas, list)
        self.assertGreaterEqual(len(formulas), 7)

    def test_sorted_by_abs_residual(self):
        formulas = self.result["formulas"]
        abs_residuals = [f["abs_residual_ppm"] for f in formulas]
        self.assertEqual(abs_residuals, sorted(abs_residuals))

    def test_hierarchy_present(self):
        labels = [f["formula_label"] for f in self.result["formulas"]]
        hierarchy_found = any("mu^2" in label and "hierarchy" in label
                              for label in labels)
        self.assertTrue(hierarchy_found,
                        f"mu^2 hierarchy formula not found in: {labels}")

    def test_sqrt_phi_present(self):
        labels = [f["formula_label"] for f in self.result["formulas"]]
        sqrt_phi_found = any("sqrt(phi)" in label for label in labels)
        self.assertTrue(sqrt_phi_found,
                        f"sqrt(phi) formula not found in: {labels}")

    def test_refined_mu_structure_best_overall(self):
        """Refined mu-structure (1-alpha) should have the smallest abs residual overall."""
        best = self.result["best_overall"]
        self.assertIn("1-alpha", best["formula_label"])

    def test_sqrt_phi_best_zero_param(self):
        """sqrt(phi) should be the best zero-parameter formula."""
        best_zp = self.result["best_zero_param"]
        self.assertIn("sqrt(phi)", best_zp["formula_label"])

    def test_each_has_required_keys(self):
        required = {
            "formula_label", "alpha_g", "G_predicted", "residual_ppm",
            "abs_residual_ppm", "n_fitted_params", "within_codata_uncertainty",
        }
        for f in self.result["formulas"]:
            self.assertTrue(required.issubset(f.keys()),
                            f"Missing keys in formula: {required - f.keys()}")

    def test_refined_formula_present(self):
        """Refined mu-structure formula should be in the comparison."""
        result = compare_all_bridges()
        labels = [f["formula_label"] for f in result["formulas"]]
        self.assertTrue(any("1-alpha" in l for l in labels))

    def test_refined_formula_sub_ppm(self):
        """Refined formula should be sub-ppm in the comparison."""
        result = compare_all_bridges()
        refined = [f for f in result["formulas"] if "1-alpha" in f["formula_label"]]
        self.assertTrue(len(refined) > 0)
        self.assertLess(refined[0]["abs_residual_ppm"], 1.0)

    def test_assessment_nonempty(self):
        self.assertIsInstance(self.result["assessment"], str)
        self.assertGreater(len(self.result["assessment"]), 20)


class TestCodataStability(unittest.TestCase):
    """Tests for analyze_codata_stability_mu."""

    @classmethod
    def setUpClass(cls):
        cls.result = analyze_codata_stability_mu()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_has_both_editions(self):
        self.assertIn("CODATA 2014", self.result)
        self.assertIn("CODATA 2018", self.result)

    def test_each_edition_keys(self):
        expected = {"G_predicted", "G_measured", "residual_ppm", "within_uncertainty"}
        for edition in ("CODATA 2014", "CODATA 2018"):
            self.assertTrue(
                expected.issubset(self.result[edition].keys()),
                f"Missing keys in {edition}: {expected - self.result[edition].keys()}",
            )

    def test_residuals_opposite_sign(self):
        """CODATA 2014 and 2018 bracket the prediction (opposite signs).

        2014 gives G=6.67408e-11 (lower), 2018 gives G=6.67430e-11 (higher).
        The mu-structure prediction sits between them, so the signed residual
        flips sign across editions.
        """
        r14 = self.result["CODATA 2014"]["residual_ppm"]
        r18 = self.result["CODATA 2018"]["residual_ppm"]
        self.assertNotEqual(
            r14 < 0, r18 < 0,
            f"Expected opposite signs: 2014={r14:.2f}, 2018={r18:.2f}",
        )

    def test_residuals_similar_magnitude(self):
        r14 = self.result["CODATA 2014"]["residual_ppm"]
        r18 = self.result["CODATA 2018"]["residual_ppm"]
        self.assertLess(abs(r14 - r18), 50)

    def test_within_uncertainty_2018(self):
        self.assertTrue(self.result["CODATA 2018"]["within_uncertainty"])


class TestSummarizeMuStructure(unittest.TestCase):
    """Tests for summarize_mu_structure."""

    @classmethod
    def setUpClass(cls):
        cls.constants = get_constants("CODATA 2018")
        cls.result = summarize_mu_structure(cls.constants)

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        expected = {
            "mu_structure", "offset_scan", "sqrt_phi_origin",
            "formula_comparison", "codata_stability",
            "key_finding", "honest_assessment",
        }
        self.assertTrue(expected.issubset(self.result.keys()),
                        f"Missing keys: {expected - self.result.keys()}")

    def test_key_finding_nonempty(self):
        self.assertIsInstance(self.result["key_finding"], str)
        self.assertGreater(len(self.result["key_finding"]), 20)

    def test_honest_assessment_nonempty(self):
        self.assertIsInstance(self.result["honest_assessment"], str)
        self.assertGreater(len(self.result["honest_assessment"]), 20)

    def test_mu_structure_is_dict(self):
        self.assertIsInstance(self.result["mu_structure"], dict)

    def test_offset_scan_is_dict(self):
        self.assertIsInstance(self.result["offset_scan"], dict)

    def test_formula_comparison_is_dict(self):
        self.assertIsInstance(self.result["formula_comparison"], dict)

    def test_refined_key_present(self):
        """Summary should include the refined formula result."""
        result = summarize_mu_structure()
        self.assertIn("refined", result)
        self.assertIsInstance(result["refined"], dict)

    def test_honest_assessment_mentions_sqrt_phi(self):
        text = self.result["honest_assessment"].lower()
        self.assertTrue(
            "sqrt" in text or "phi" in text or "golden" in text,
            f"honest_assessment should mention sqrt, phi, or golden ratio: {text[:200]}",
        )


if __name__ == "__main__":
    unittest.main()
