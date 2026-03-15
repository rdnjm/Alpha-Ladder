"""Tests for unified_formula module."""

import unittest
from alpha_ladder_core.constants import get_constants
from alpha_ladder_core import unified_formula as uf


class TestExactBridgeCoefficient(unittest.TestCase):
    """Tests for compute_exact_bridge_coefficient."""

    @classmethod
    def setUpClass(cls):
        cls.result = uf.compute_exact_bridge_coefficient()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {
            "C_exact", "C_tree", "fractional_correction",
            "c2_exact", "c3_from_exact", "k_exact",
        }
        self.assertTrue(required.issubset(self.result.keys()))

    def test_C_exact_between_tree_and_hierarchy(self):
        """C_tree (phi^2/2) < C_exact < alpha^3*mu^2."""
        C_exact = float(self.result["C_exact"])
        C_tree = float(self.result["C_tree"])
        self.assertGreater(C_exact, C_tree)
        # C_exact ~ 1.3092, should be near 1.309
        self.assertAlmostEqual(C_exact, 1.3092, delta=0.001)

    def test_c2_exact_near_3(self):
        """c2_exact should be between 3.0 and 3.1."""
        c2 = self.result["c2_exact"]
        self.assertGreater(c2, 3.0)
        self.assertLess(c2, 3.1)

    def test_c3_near_8_5(self):
        """c3_from_exact should be between 1.5 and 1.7 (near 8/5 = 1.6)."""
        c3 = self.result["c3_from_exact"]
        self.assertGreater(c3, 1.5)
        self.assertLess(c3, 1.7)

    def test_k_exact_less_than_sqrt_phi(self):
        """k_exact should be less than sqrt(phi) = 1.2720."""
        k = self.result["k_exact"]
        self.assertLess(k, 1.2720)
        # Also should be positive and close to sqrt(phi)
        self.assertGreater(k, 1.25)

    def test_fractional_correction_small(self):
        """Fractional correction should be small (< 0.001)."""
        self.assertLess(abs(self.result["fractional_correction"]), 0.001)


class TestScanResummedBridges(unittest.TestCase):
    """Tests for scan_resummed_bridges."""

    @classmethod
    def setUpClass(cls):
        cls.result = uf.scan_resummed_bridges()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {"candidates", "best_candidate", "assessment"}
        self.assertTrue(required.issubset(self.result.keys()))

    def test_candidates_list(self):
        """Should have at least 5 candidates."""
        self.assertGreaterEqual(len(self.result["candidates"]), 5)

    def test_best_sub_ppm(self):
        """Best candidate should have residual < 1 ppm."""
        best = self.result["best_candidate"]
        self.assertLess(best["abs_residual_ppm"], 1.0)

    def test_polynomial_nlo_near_top(self):
        """The NLO polynomial (1+3a^2+8/5*a^3) should be among the top 3."""
        candidates = self.result["candidates"]
        nlo_rank = None
        for i, c in enumerate(candidates):
            if "8/5" in c["label"] and "exp" not in c["label"] and "(1+" not in c["label"].split("8/5")[0][-3:]:
                # The plain polynomial NLO
                if "1 + 3*a^2 + 8/5*a^3" in c["label"]:
                    nlo_rank = i + 1
                    break
        self.assertIsNotNone(nlo_rank, "NLO polynomial not found in candidates")
        self.assertLessEqual(nlo_rank, 3)

    def test_each_candidate_keys(self):
        """Each candidate should have required keys."""
        required = {"label", "C_value", "residual_ppm", "implied_k",
                    "n_empirical_coefficients"}
        for c in self.result["candidates"]:
            self.assertTrue(
                required.issubset(c.keys()),
                f"Candidate '{c.get('label', '?')}' missing keys: "
                f"{required - set(c.keys())}",
            )

    def test_assessment_nonempty(self):
        self.assertGreater(len(self.result["assessment"]), 0)


class TestScanMuStructureCorrections(unittest.TestCase):
    """Tests for scan_mu_structure_corrections."""

    @classmethod
    def setUpClass(cls):
        cls.result = uf.scan_mu_structure_corrections()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {"candidates", "best_candidate", "assessment"}
        self.assertTrue(required.issubset(self.result.keys()))

    def test_candidates_list(self):
        """Should have at least 8 candidates."""
        self.assertGreaterEqual(len(self.result["candidates"]), 8)

    def test_sqrt_phi_near_best(self):
        """sqrt(phi) baseline should be among the top 5."""
        candidates = self.result["candidates"]
        sqrt_phi_rank = None
        for i, c in enumerate(candidates):
            if "baseline" in c["label"]:
                sqrt_phi_rank = i + 1
                break
        self.assertIsNotNone(sqrt_phi_rank, "sqrt(phi) baseline not found")
        self.assertLessEqual(sqrt_phi_rank, 5)

    def test_each_candidate_keys(self):
        """Each candidate should have required keys."""
        required = {"label", "k_value", "C_value", "residual_ppm"}
        for c in self.result["candidates"]:
            self.assertTrue(
                required.issubset(c.keys()),
                f"Candidate '{c.get('label', '?')}' missing keys: "
                f"{required - set(c.keys())}",
            )

    def test_assessment_nonempty(self):
        self.assertGreater(len(self.result["assessment"]), 0)


class TestAnalyzeGapStructure(unittest.TestCase):
    """Tests for analyze_gap_structure."""

    @classmethod
    def setUpClass(cls):
        cls.result = uf.analyze_gap_structure()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {
            "gap_ppm", "c2_exact", "c3_exact", "k_exact", "k_shortfall",
            "path_a_viable", "path_b_viable", "path_c_candidates",
            "fundamental_tension", "assessment",
        }
        self.assertTrue(required.issubset(self.result.keys()))

    def test_gap_ppm_positive(self):
        """Gap should be between 3 and 8 ppm."""
        gap = self.result["gap_ppm"]
        self.assertGreater(gap, 3.0)
        self.assertLess(gap, 8.0)

    def test_c2_near_3(self):
        """c2_exact should be between 3.0 and 3.02."""
        c2 = self.result["c2_exact"]
        self.assertGreater(c2, 3.0)
        self.assertLess(c2, 3.02)

    def test_c3_near_1_6(self):
        """c3_exact should be between 1.5 and 1.7."""
        c3 = self.result["c3_exact"]
        self.assertGreater(c3, 1.5)
        self.assertLess(c3, 1.7)

    def test_k_shortfall_positive(self):
        """sqrt(phi) > k_exact, so k_shortfall > 0."""
        self.assertGreater(self.result["k_shortfall"], 0)

    def test_path_a_viable(self):
        """Path A (3 = d-1) should be viable."""
        self.assertTrue(self.result["path_a_viable"])

    def test_assessment_nonempty(self):
        self.assertGreater(len(self.result["assessment"]), 0)


class TestPredictMuFromUnified(unittest.TestCase):
    """Tests for predict_mu_from_unified."""

    @classmethod
    def setUpClass(cls):
        cls.result = uf.predict_mu_from_unified()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {"predictions", "best_prediction", "original_tension",
                    "assessment"}
        self.assertTrue(required.issubset(self.result.keys()))

    def test_predictions_list(self):
        """Should have at least 3 predictions."""
        self.assertGreaterEqual(len(self.result["predictions"]), 3)

    def test_each_prediction_keys(self):
        """Each prediction should have required keys."""
        required = {"formula_label", "mu_predicted", "residual_ppm",
                    "sigma_tension"}
        for p in self.result["predictions"]:
            self.assertTrue(
                required.issubset(p.keys()),
                f"Prediction '{p.get('formula_label', '?')}' missing keys: "
                f"{required - set(p.keys())}",
            )

    def test_mu_predicted_range(self):
        """All mu predictions should be between 1836.0 and 1836.5."""
        for p in self.result["predictions"]:
            self.assertGreater(p["mu_predicted"], 1836.0)
            self.assertLess(p["mu_predicted"], 1836.5)

    def test_original_tension_value(self):
        """Original tension should be between 2 and 3 ppm."""
        tension = self.result["original_tension"]
        self.assertIsInstance(tension, dict)
        ppm = tension["residual_ppm"]
        self.assertGreater(ppm, 2.0)
        self.assertLess(ppm, 3.0)

    def test_assessment_nonempty(self):
        self.assertGreater(len(self.result["assessment"]), 0)


class TestSummarizeUnifiedFormula(unittest.TestCase):
    """Tests for summarize_unified_formula."""

    @classmethod
    def setUpClass(cls):
        cls.result = uf.summarize_unified_formula()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {
            "exact_coefficient", "resummed_bridges", "mu_corrections",
            "gap_analysis", "mu_predictions", "key_finding",
            "honest_assessment",
        }
        self.assertEqual(required, set(self.result.keys()))

    def test_key_finding_nonempty(self):
        self.assertGreater(len(self.result["key_finding"]), 0)

    def test_honest_assessment_nonempty(self):
        self.assertGreater(len(self.result["honest_assessment"]), 0)

    def test_all_sub_dicts(self):
        """Each sub-result should be a dict."""
        for key in ("exact_coefficient", "resummed_bridges",
                    "mu_corrections", "gap_analysis", "mu_predictions"):
            self.assertIsInstance(
                self.result[key], dict,
                f"{key} should be a dict",
            )


if __name__ == "__main__":
    unittest.main()
