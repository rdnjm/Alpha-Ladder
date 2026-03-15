"""Tests for literature_comparison module."""

import math
import unittest
from alpha_ladder_core.constants import get_constants
from alpha_ladder_core import literature_comparison as lc


class TestBeckCosmologicalConstant(unittest.TestCase):
    """Tests for compare_beck_cosmological_constant."""

    @classmethod
    def setUpClass(cls):
        cls.result = lc.compare_beck_cosmological_constant()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {
            "paper", "arxiv_id", "formula_label", "Lambda_beck",
            "Lambda_observed", "ratio", "log10_ratio", "beck_agreement",
            "alpha_ladder_connection", "assessment",
        }
        self.assertTrue(required.issubset(self.result.keys()),
                        f"Missing keys: {required - self.result.keys()}")

    def test_paper_name(self):
        self.assertEqual(self.result["paper"], "Beck (2008)")

    def test_arxiv_id(self):
        self.assertEqual(self.result["arxiv_id"], "0810.0752")

    def test_lambda_beck_positive(self):
        self.assertGreater(float(self.result["Lambda_beck"]), 0)

    def test_lambda_observed_positive(self):
        self.assertGreater(float(self.result["Lambda_observed"]), 0)

    def test_ratio_finite(self):
        ratio = float(self.result["ratio"])
        self.assertTrue(math.isfinite(ratio),
                        f"ratio should be finite, got {ratio}")

    def test_assessment_nonempty(self):
        self.assertIsInstance(self.result["assessment"], str)
        self.assertGreater(len(self.result["assessment"]), 0)


class TestAlexanderHierarchy(unittest.TestCase):
    """Tests for compare_alexander_hierarchy."""

    @classmethod
    def setUpClass(cls):
        cls.result = lc.compare_alexander_hierarchy()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {
            "paper", "arxiv_id", "bound", "alpha_g_over_alpha",
            "log10_ratio", "satisfies_bound", "hierarchy_prediction",
            "log10_hierarchy", "margin_orders", "connection_to_ladder",
            "assessment",
        }
        self.assertTrue(required.issubset(self.result.keys()),
                        f"Missing keys: {required - self.result.keys()}")

    def test_satisfies_bound(self):
        self.assertTrue(self.result["satisfies_bound"])

    def test_log10_ratio_range(self):
        log10_ratio = float(self.result["log10_ratio"])
        self.assertGreater(log10_ratio, -50)
        self.assertLess(log10_ratio, -35)

    def test_margin_orders_positive(self):
        self.assertGreater(float(self.result["margin_orders"]), 0)

    def test_hierarchy_prediction_positive(self):
        self.assertGreater(float(self.result["hierarchy_prediction"]), 0)

    def test_bound_value(self):
        self.assertAlmostEqual(float(self.result["bound"]), 1e-34)

    def test_assessment_nonempty(self):
        self.assertIsInstance(self.result["assessment"], str)
        self.assertGreater(len(self.result["assessment"]), 0)


class TestEavesLogarithmic(unittest.TestCase):
    """Tests for compare_eaves_logarithmic."""

    @classmethod
    def setUpClass(cls):
        cls.result = lc.compare_eaves_logarithmic()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {
            "paper", "arxiv_id", "formula_label", "r_e", "l_P",
            "V_e", "V_P", "ln_ratio", "alpha_inverse",
            "fractional_difference", "connection_to_ladder", "assessment",
        }
        self.assertTrue(required.issubset(self.result.keys()),
                        f"Missing keys: {required - self.result.keys()}")

    def test_r_e_positive(self):
        self.assertGreater(float(self.result["r_e"]), 0)

    def test_l_P_positive(self):
        self.assertGreater(float(self.result["l_P"]), 0)

    def test_r_e_greater_than_l_P(self):
        self.assertGreater(float(self.result["r_e"]),
                           float(self.result["l_P"]))

    def test_alpha_inverse_near_137(self):
        alpha_inv = float(self.result["alpha_inverse"])
        self.assertGreater(alpha_inv, 137.03)
        self.assertLess(alpha_inv, 137.04)

    def test_ln_ratio_order_of_magnitude(self):
        ln_ratio = float(self.result["ln_ratio"])
        self.assertGreater(ln_ratio, 100)

    def test_assessment_nonempty(self):
        self.assertIsInstance(self.result["assessment"], str)
        self.assertGreater(len(self.result["assessment"]), 0)


class TestBlauVisserWipf(unittest.TestCase):
    """Tests for compare_blau_visser_wipf."""

    @classmethod
    def setUpClass(cls):
        cls.result = lc.compare_blau_visser_wipf()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {
            "paper", "arxiv_id", "method", "our_result",
            "polynomial_identity", "polynomial_verified", "implication",
            "methodological_agreement", "assessment",
        }
        self.assertTrue(required.issubset(self.result.keys()),
                        f"Missing keys: {required - self.result.keys()}")

    def test_polynomial_verified(self):
        self.assertTrue(self.result["polynomial_verified"])

    def test_method_mentions_zeta(self):
        self.assertIn("zeta", self.result["method"].lower())

    def test_our_result_mentions_zero(self):
        self.assertIn("0", str(self.result["our_result"]))

    def test_assessment_nonempty(self):
        self.assertIsInstance(self.result["assessment"], str)
        self.assertGreater(len(self.result["assessment"]), 0)


class TestEddingtonDirac(unittest.TestCase):
    """Tests for compare_eddington_dirac."""

    @classmethod
    def setUpClass(cls):
        cls.result = lc.compare_eddington_dirac()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {
            "hypothesis", "N_ED", "log10_N_ED", "N_ED_via_ladder",
            "log10_ladder", "agreement_with_measured",
            "ladder_provides_explanation", "assessment",
        }
        self.assertTrue(required.issubset(self.result.keys()),
                        f"Missing keys: {required - self.result.keys()}")

    def test_N_ED_large(self):
        self.assertGreater(float(self.result["N_ED"]), 1e35)

    def test_log10_N_ED_range(self):
        log10_val = float(self.result["log10_N_ED"])
        self.assertGreater(log10_val, 40)
        self.assertLess(log10_val, 45)

    def test_log10_ladder_close(self):
        diff = abs(float(self.result["log10_ladder"]) -
                   float(self.result["log10_N_ED"]))
        self.assertLess(diff, 1,
                        "Ladder log10 should be within 1 order of N_ED")

    def test_agreement_small(self):
        self.assertLess(float(self.result["agreement_with_measured"]), 1)

    def test_assessment_nonempty(self):
        self.assertIsInstance(self.result["assessment"], str)
        self.assertGreater(len(self.result["assessment"]), 0)


class TestSummarizeLiterature(unittest.TestCase):
    """Tests for summarize_literature_comparison."""

    @classmethod
    def setUpClass(cls):
        cls.result = lc.summarize_literature_comparison()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_required_keys(self):
        required = {
            "beck", "alexander", "eaves", "blau_visser_wipf",
            "eddington_dirac", "papers_analyzed", "key_finding",
            "honest_assessment",
        }
        self.assertTrue(required.issubset(self.result.keys()),
                        f"Missing keys: {required - self.result.keys()}")

    def test_papers_analyzed_length(self):
        self.assertGreaterEqual(len(self.result["papers_analyzed"]), 5)

    def test_key_finding_nonempty(self):
        self.assertIsInstance(self.result["key_finding"], str)
        self.assertGreater(len(self.result["key_finding"]), 0)

    def test_honest_assessment_nonempty(self):
        self.assertIsInstance(self.result["honest_assessment"], str)
        self.assertGreater(len(self.result["honest_assessment"]), 0)


if __name__ == "__main__":
    unittest.main()
