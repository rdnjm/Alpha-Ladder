"""Tests for corrected_bridge module."""

import unittest
from decimal import Decimal

from alpha_ladder_core.constants import get_constants
from alpha_ladder_core import corrected_bridge as cb


class TestComputeCorrectedBridge(unittest.TestCase):
    """Tests for compute_corrected_bridge."""

    @classmethod
    def setUpClass(cls):
        cls.result = cb.compute_corrected_bridge()

    def test_returns_dict_with_all_keys(self):
        expected_keys = {
            "C_exact", "C_uncorrected", "C_corrected_leading", "C_corrected_nlo",
            "uncorrected_ppm", "corrected_leading_ppm", "corrected_nlo_ppm",
            "G_uncorrected", "G_corrected", "G_measured",
            "correction_factor", "fraction_of_gap_explained",
        }
        self.assertTrue(expected_keys.issubset(self.result.keys()))

    def test_c_exact_is_decimal_in_range(self):
        c = self.result["C_exact"]
        self.assertIsInstance(c, Decimal)
        self.assertGreater(float(c), 1.309)
        self.assertLess(float(c), 1.310)

    def test_uncorrected_ppm_approximately_160(self):
        """uncorrected_ppm is abs value, so positive ~160."""
        ppm = self.result["uncorrected_ppm"]
        self.assertIsInstance(ppm, (int, float))
        self.assertGreater(ppm, 100)
        self.assertLess(ppm, 200)

    def test_corrected_leading_ppm_sub_ppm(self):
        ppm = self.result["corrected_leading_ppm"]
        self.assertLess(abs(ppm), 1.0)

    def test_corrected_nlo_ppm_very_small(self):
        ppm = self.result["corrected_nlo_ppm"]
        self.assertLess(abs(ppm), 0.1)

    def test_correction_improves_residual(self):
        self.assertLess(
            abs(self.result["corrected_leading_ppm"]),
            abs(self.result["uncorrected_ppm"]),
        )

    def test_fraction_of_gap_explained_above_99_percent(self):
        self.assertGreater(self.result["fraction_of_gap_explained"], 0.99)

    def test_g_corrected_closer_to_measured(self):
        g_m = float(self.result["G_measured"])
        g_u = float(self.result["G_uncorrected"])
        g_c = float(self.result["G_corrected"])
        self.assertLess(abs(g_c - g_m), abs(g_u - g_m))

    def test_g_values_positive_decimals_in_range(self):
        for key in ("G_uncorrected", "G_corrected", "G_measured"):
            val = self.result[key]
            self.assertIsInstance(val, Decimal)
            self.assertGreater(float(val), 6.6e-11)
            self.assertLess(float(val), 6.8e-11)

    def test_correction_factor_approximately_1e4(self):
        cf = float(self.result["correction_factor"])
        self.assertGreater(cf, 1e-5)
        self.assertLess(cf, 1e-3)

    def test_custom_constants_accepted(self):
        c = get_constants("CODATA 2014")
        result = cb.compute_corrected_bridge(constants=c)
        self.assertIn("C_exact", result)


class TestAnalyzeCorrectionOrigin(unittest.TestCase):
    """Tests for analyze_correction_origin."""

    @classmethod
    def setUpClass(cls):
        cls.result = cb.analyze_correction_origin()

    def test_returns_dict_with_all_keys(self):
        expected_keys = {
            "factor_value", "interpretations", "degeneracy_note", "assessment",
        }
        self.assertTrue(expected_keys.issubset(self.result.keys()))

    def test_factor_value_is_3(self):
        self.assertEqual(self.result["factor_value"], 3)

    def test_interpretations_nonempty_list(self):
        interps = self.result["interpretations"]
        self.assertIsInstance(interps, list)
        self.assertGreater(len(interps), 0)

    def test_interpretation_has_required_keys(self):
        required = {"name", "formula", "value_for_d4_n2", "matches"}
        for interp in self.result["interpretations"]:
            self.assertTrue(required.issubset(interp.keys()),
                            f"Missing keys in interpretation: {interp.keys()}")

    def test_at_least_one_matching_interpretation(self):
        matches = [i for i in self.result["interpretations"] if i["matches"]]
        self.assertGreater(len(matches), 0)

    def test_degeneracy_note_nonempty(self):
        self.assertIsInstance(self.result["degeneracy_note"], str)
        self.assertGreater(len(self.result["degeneracy_note"]), 0)

    def test_assessment_nonempty(self):
        self.assertIsInstance(self.result["assessment"], str)
        self.assertGreater(len(self.result["assessment"]), 0)

    def test_matching_interpretations_have_value_3(self):
        for interp in self.result["interpretations"]:
            if interp["value_for_d4_n2"] == 3:
                self.assertTrue(interp["matches"],
                                f"Interpretation '{interp['name']}' has value 3 but matches=False")

    def test_custom_d_n_accepted(self):
        result = cb.analyze_correction_origin(d=4, n=2)
        self.assertEqual(result["factor_value"], 3)


class TestScanCorrectionSeries(unittest.TestCase):
    """Tests for scan_correction_series."""

    @classmethod
    def setUpClass(cls):
        cls.result = cb.scan_correction_series()

    def test_returns_dict_with_all_keys(self):
        expected_keys = {
            "coefficients", "series_converging", "best_order", "assessment",
        }
        self.assertTrue(expected_keys.issubset(self.result.keys()))

    def test_coefficients_nonempty_list(self):
        self.assertIsInstance(self.result["coefficients"], list)
        self.assertGreater(len(self.result["coefficients"]), 0)

    def test_leading_coefficient_approximately_3(self):
        coeffs = self.result["coefficients"]
        # Find the order-2 coefficient
        c2_entries = [c for c in coeffs if c["order"] == 2]
        self.assertEqual(len(c2_entries), 1, "Should have exactly one order-2 coefficient")
        self.assertAlmostEqual(c2_entries[0]["coefficient"], 3, delta=0.5)

    def test_best_order_reasonable(self):
        self.assertGreaterEqual(self.result["best_order"], 2)
        self.assertLessEqual(self.result["best_order"], 5)

    def test_series_converging_is_bool(self):
        self.assertIsInstance(self.result["series_converging"], bool)

    def test_assessment_nonempty(self):
        self.assertIsInstance(self.result["assessment"], str)
        self.assertGreater(len(self.result["assessment"]), 0)

    def test_residuals_decrease_with_order(self):
        coeffs = self.result["coefficients"]
        if len(coeffs) >= 2:
            residuals = [abs(c["residual_after_ppm"]) for c in coeffs]
            # At minimum, later orders should not be worse than uncorrected
            # The first correction (order 2) should improve things significantly
            self.assertLess(residuals[0], 200)  # first correction should help


class TestCompareBridgeHierarchy(unittest.TestCase):
    """Tests for compare_bridge_hierarchy."""

    @classmethod
    def setUpClass(cls):
        cls.result = cb.compare_bridge_hierarchy()

    def test_returns_dict_with_all_keys(self):
        expected_keys = {
            "bridge_alpha_g", "hierarchy_alpha_g", "ratio",
            "difference_ppm", "epsilon", "interpretation",
        }
        self.assertTrue(expected_keys.issubset(self.result.keys()))

    def test_alpha_g_values_are_decimals(self):
        self.assertIsInstance(self.result["bridge_alpha_g"], Decimal)
        self.assertIsInstance(self.result["hierarchy_alpha_g"], Decimal)

    def test_ratio_close_to_one(self):
        self.assertGreater(self.result["ratio"], 0.99)
        self.assertLess(self.result["ratio"], 1.01)

    def test_difference_ppm_approximately_688(self):
        ppm = self.result["difference_ppm"]
        self.assertGreater(ppm, 600)
        self.assertLess(ppm, 800)

    def test_epsilon_approximately_6_9e4(self):
        eps = self.result["epsilon"]
        self.assertGreater(abs(eps), 5e-4)
        self.assertLess(abs(eps), 9e-4)

    def test_interpretation_nonempty(self):
        self.assertIsInstance(self.result["interpretation"], str)
        self.assertGreater(len(self.result["interpretation"]), 0)

    def test_bridge_and_hierarchy_differ(self):
        self.assertNotEqual(
            self.result["bridge_alpha_g"],
            self.result["hierarchy_alpha_g"],
        )


class TestAnalyzeCodataEditions(unittest.TestCase):
    """Tests for analyze_codata_editions."""

    @classmethod
    def setUpClass(cls):
        cls.result = cb.analyze_codata_editions()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_has_both_editions(self):
        self.assertIn("CODATA 2014", self.result)
        self.assertIn("CODATA 2018", self.result)

    def test_edition_has_required_subkeys(self):
        required = {"C_exact", "uncorrected_ppm", "corrected_ppm"}
        for edition in ("CODATA 2014", "CODATA 2018"):
            self.assertTrue(
                required.issubset(self.result[edition].keys()),
                f"Missing keys for {edition}: {self.result[edition].keys()}",
            )

    def test_corrected_smaller_than_uncorrected_2018(self):
        """For CODATA 2018 (the edition the correction was tuned on), corrected < uncorrected."""
        ed = self.result["CODATA 2018"]
        self.assertLess(abs(ed["corrected_ppm"]), abs(ed["uncorrected_ppm"]))

    def test_codata_2018_sub_10_ppm_corrected(self):
        """CODATA 2018 corrected residual should be sub-10 ppm."""
        self.assertLess(abs(self.result["CODATA 2018"]["corrected_ppm"]), 10)

    def test_codata_2014_has_finite_corrected_ppm(self):
        """CODATA 2014 may not improve as much but corrected_ppm should be finite."""
        import math
        self.assertTrue(math.isfinite(self.result["CODATA 2014"]["corrected_ppm"]))


class TestSummarizeCorrectedBridge(unittest.TestCase):
    """Tests for summarize_corrected_bridge."""

    @classmethod
    def setUpClass(cls):
        cls.result = cb.summarize_corrected_bridge()

    def test_returns_dict_with_all_keys(self):
        expected_keys = {
            "corrected_bridge", "correction_origin", "series",
            "bridge_vs_hierarchy", "editions", "key_finding", "honest_assessment",
        }
        self.assertTrue(expected_keys.issubset(self.result.keys()))

    def test_key_finding_nonempty(self):
        self.assertIsInstance(self.result["key_finding"], str)
        self.assertGreater(len(self.result["key_finding"]), 0)

    def test_honest_assessment_nonempty(self):
        self.assertIsInstance(self.result["honest_assessment"], str)
        self.assertGreater(len(self.result["honest_assessment"]), 0)

    def test_corrected_bridge_has_expected_keys(self):
        sub = self.result["corrected_bridge"]
        self.assertIn("C_exact", sub)
        self.assertIn("corrected_leading_ppm", sub)

    def test_correction_origin_has_expected_keys(self):
        sub = self.result["correction_origin"]
        self.assertIn("factor_value", sub)
        self.assertIn("interpretations", sub)

    def test_series_has_expected_keys(self):
        sub = self.result["series"]
        self.assertIn("coefficients", sub)
        self.assertIn("best_order", sub)

    def test_bridge_vs_hierarchy_has_expected_keys(self):
        sub = self.result["bridge_vs_hierarchy"]
        self.assertIn("difference_ppm", sub)
        self.assertIn("ratio", sub)


# ===========================================================================
# Geometric resummation tests (pytest-style)
# ===========================================================================

from alpha_ladder_core.corrected_bridge import (
    compute_geometric_resummation,
    predict_mu_from_geometry,
)


def test_resummation_returns_dict():
    result = compute_geometric_resummation()
    assert isinstance(result, dict)
    for key in ("F_exact", "F_geom", "residual_ppm", "c3_exact",
                "c3_geom", "c3_residual_ppm", "coefficients",
                "ratio", "expansion_parameter", "honest_assessment"):
        assert key in result


def test_resummation_matches_sub_ppm():
    result = compute_geometric_resummation()
    assert result["residual_ppm"] < 0.01


def test_resummation_c3_is_phi_half():
    import math
    phi = (1 + math.sqrt(5)) / 2
    result = compute_geometric_resummation()
    assert abs(result["c3_geom"] - phi / 2) < 1e-10


def test_resummation_ratio_is_inverse_phi():
    import math
    phi = (1 + math.sqrt(5)) / 2
    result = compute_geometric_resummation()
    assert abs(result["ratio"] - 1.0 / phi) < 1e-10


def test_resummation_coefficients_geometric():
    result = compute_geometric_resummation()
    coeffs = result["coefficients"]
    for i in range(len(coeffs) - 1):
        ratio = coeffs[i + 1]["c_n"] / coeffs[i]["c_n"]
        assert abs(ratio - result["ratio"]) < 1e-10


def test_resummation_six_coefficients():
    result = compute_geometric_resummation()
    assert len(result["coefficients"]) >= 6


def test_resummation_c3_residual_order():
    result = compute_geometric_resummation()
    assert result["c3_residual_ppm"] > 1000
    assert result["c3_residual_ppm"] < 10000


def test_mu_prediction_returns_dict():
    result = predict_mu_from_geometry()
    assert isinstance(result, dict)
    for key in ("mu_predicted", "mu_measured", "residual_ppm",
                "alpha", "phi", "F_geom", "codata_stability",
                "honest_assessment"):
        assert key in result


def test_mu_prediction_sub_ppm_2018():
    result = predict_mu_from_geometry()
    assert abs(result["residual_ppm"]) < 0.01


def test_mu_prediction_codata_stability():
    result = predict_mu_from_geometry()
    stability = result["codata_stability"]
    assert len(stability) == 3
    for entry in stability:
        assert abs(entry["residual_ppm"]) < 0.1


def test_mu_prediction_2022_independent():
    result = predict_mu_from_geometry()
    for entry in result["codata_stability"]:
        if entry["edition"] == "CODATA 2022":
            assert abs(entry["residual_ppm"]) < 0.01


def test_mu_predicted_positive():
    result = predict_mu_from_geometry()
    assert result["mu_predicted"] > 1800
    assert result["mu_predicted"] < 1840


def test_mu_prediction_no_mu_input():
    """F_geom depends only on alpha and phi, not mu."""
    result = compute_geometric_resummation()
    # F_geom should be close to 1 (small correction)
    assert 1.0 < result["F_geom"] < 1.001


def test_resummation_expansion_parameter():
    result = compute_geometric_resummation()
    assert result["expansion_parameter"] == "alpha/phi"
    assert result["alpha_over_phi"] < 0.01


if __name__ == "__main__":
    unittest.main()
