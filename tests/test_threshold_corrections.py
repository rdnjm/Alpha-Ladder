"""Tests for alpha_ladder_core.threshold_corrections module.

Validates threshold correction physics: KK log-sum structure,
matching to 3*alpha^2 correction, scale scanning, and overall
threshold analysis summary.
"""

import math
import unittest

from alpha_ladder_core.threshold_corrections import (
    compute_threshold_sum,
    compute_threshold_correction_to_G,
    analyze_log_sum_structure,
    scan_matching_scales,
    summarize_threshold_analysis,
)


# ---------------------------------------------------------------------------
# TestComputeThresholdSum
# ---------------------------------------------------------------------------
class TestComputeThresholdSum(unittest.TestCase):
    """Tests for compute_threshold_sum."""

    def test_returns_dict_with_required_keys(self):
        result = compute_threshold_sum(L_max=10)
        self.assertIsInstance(result, dict)
        for key in ("L_max", "S_log", "asymptotic_terms", "finite_part",
                     "convergence_data", "converged"):
            self.assertIn(key, result, f"Missing key: {key}")

    def test_s_log_is_positive(self):
        """S_log is a sum of positive terms since l*(l+1) > 1 for l >= 1."""
        result = compute_threshold_sum(L_max=10)
        self.assertGreater(result["S_log"], 0)

    def test_s_log_monotonically_increasing(self):
        """S_log(100) > S_log(10) because we are adding positive terms."""
        r10 = compute_threshold_sum(L_max=10)
        r100 = compute_threshold_sum(L_max=100)
        self.assertGreater(r100["S_log"], r10["S_log"])

    def test_finite_part_is_finite(self):
        result = compute_threshold_sum(L_max=100)
        self.assertTrue(math.isfinite(result["finite_part"]),
                        "finite_part should not be inf or nan")

    def test_convergence_data_nonempty(self):
        result = compute_threshold_sum(L_max=50)
        self.assertIsInstance(result["convergence_data"], list)
        self.assertGreater(len(result["convergence_data"]), 0)

    def test_l_max_matches_input(self):
        result = compute_threshold_sum(L_max=42)
        self.assertEqual(result["L_max"], 42)

    def test_s_log_10_reasonable_magnitude(self):
        """For L_max=10, S_log should be a moderate positive number.

        Manual partial check: l=1 contributes 3*2*ln(2) ~ 4.16,
        l=2 contributes 5*6*ln(6) ~ 53.75, so total > 50.
        """
        result = compute_threshold_sum(L_max=10)
        self.assertGreater(result["S_log"], 50)

    def test_asymptotic_terms_is_dict(self):
        result = compute_threshold_sum(L_max=100)
        self.assertIsInstance(result["asymptotic_terms"], dict)

    def test_converged_is_bool(self):
        result = compute_threshold_sum(L_max=100)
        self.assertIsInstance(result["converged"], bool)


# ---------------------------------------------------------------------------
# TestComputeThresholdCorrectionToG
# ---------------------------------------------------------------------------
class TestComputeThresholdCorrectionToG(unittest.TestCase):
    """Tests for compute_threshold_correction_to_G."""

    def test_returns_dict_with_required_keys(self):
        result = compute_threshold_correction_to_G()
        self.assertIsInstance(result, dict)
        for key in ("finite_part_S_log", "R_choices", "sector_results",
                     "target", "best_match", "residual_ppm", "assessment"):
            self.assertIn(key, result, f"Missing key: {key}")

    def test_target_approximately_3_alpha_sq(self):
        """target should be approximately 3*alpha^2 ~ 1.6e-4."""
        result = compute_threshold_correction_to_G()
        target = result["target"]
        # Within a factor of 10 of 1.6e-4
        self.assertGreater(target, 1e-5)
        self.assertLess(target, 2e-3)

    def test_sector_results_nonempty(self):
        result = compute_threshold_correction_to_G()
        self.assertIsInstance(result["sector_results"], dict)
        self.assertGreater(len(result["sector_results"]), 0)

    def test_assessment_nonempty_string(self):
        result = compute_threshold_correction_to_G()
        self.assertIsInstance(result["assessment"], str)
        self.assertGreater(len(result["assessment"]), 0)

    def test_finite_part_s_log_is_finite(self):
        result = compute_threshold_correction_to_G()
        self.assertTrue(math.isfinite(result["finite_part_S_log"]))

    def test_r_choices_has_entries(self):
        result = compute_threshold_correction_to_G()
        self.assertIsInstance(result["R_choices"], dict)
        self.assertGreaterEqual(len(result["R_choices"]), 1)

    def test_residual_ppm_is_finite(self):
        result = compute_threshold_correction_to_G()
        self.assertTrue(math.isfinite(result["residual_ppm"]))


# ---------------------------------------------------------------------------
# TestAnalyzeLogSumStructure
# ---------------------------------------------------------------------------
class TestAnalyzeLogSumStructure(unittest.TestCase):
    """Tests for analyze_log_sum_structure."""

    def test_returns_dict_with_required_keys(self):
        result = analyze_log_sum_structure(L_max=100)
        self.assertIsInstance(result, dict)
        for key in ("values", "effective_log", "finite_part_vs_L",
                     "assessment"):
            self.assertIn(key, result, f"Missing key: {key}")

    def test_values_nonempty_list(self):
        result = analyze_log_sum_structure(L_max=100)
        self.assertIsInstance(result["values"], list)
        self.assertGreater(len(result["values"]), 0)

    def test_effective_log_is_finite(self):
        result = analyze_log_sum_structure(L_max=100)
        self.assertTrue(math.isfinite(result["effective_log"]))

    def test_assessment_nonempty_string(self):
        result = analyze_log_sum_structure(L_max=100)
        self.assertIsInstance(result["assessment"], str)
        self.assertGreater(len(result["assessment"]), 0)

    def test_finite_part_vs_l_is_list(self):
        result = analyze_log_sum_structure(L_max=100)
        self.assertIsInstance(result["finite_part_vs_L"], list)


# ---------------------------------------------------------------------------
# TestScanMatchingScales
# ---------------------------------------------------------------------------
class TestScanMatchingScales(unittest.TestCase):
    """Tests for scan_matching_scales."""

    def test_returns_dict_with_required_keys(self):
        result = scan_matching_scales()
        self.assertIsInstance(result, dict)
        for key in ("matching_scales", "best_scale", "assessment"):
            self.assertIn(key, result, f"Missing key: {key}")

    def test_matching_scales_has_entries(self):
        result = scan_matching_scales()
        self.assertIsInstance(result["matching_scales"], dict)
        self.assertGreaterEqual(len(result["matching_scales"]), 2)

    def test_best_scale_present(self):
        result = scan_matching_scales()
        self.assertIsNotNone(result["best_scale"])

    def test_assessment_nonempty_string(self):
        result = scan_matching_scales()
        self.assertIsInstance(result["assessment"], str)
        self.assertGreater(len(result["assessment"]), 0)

    def test_matching_scales_values_are_numeric_or_dict(self):
        """Each matching scale entry should contain numeric or dict data."""
        result = scan_matching_scales()
        for key, val in result["matching_scales"].items():
            self.assertTrue(
                isinstance(val, (int, float, dict)),
                f"matching_scales['{key}'] should be numeric or dict, "
                f"got {type(val)}"
            )


# ---------------------------------------------------------------------------
# TestSummarizeThresholdAnalysis
# ---------------------------------------------------------------------------
class TestSummarizeThresholdAnalysis(unittest.TestCase):
    """Tests for summarize_threshold_analysis."""

    def test_returns_dict_with_required_keys(self):
        result = summarize_threshold_analysis()
        self.assertIsInstance(result, dict)
        for key in ("threshold_sum", "correction_to_G", "log_structure",
                     "matching_scales", "matches_3_alpha_sq", "key_finding",
                     "honest_assessment"):
            self.assertIn(key, result, f"Missing key: {key}")

    def test_key_finding_nonempty_string(self):
        result = summarize_threshold_analysis()
        self.assertIsInstance(result["key_finding"], str)
        self.assertGreater(len(result["key_finding"]), 0)

    def test_honest_assessment_nonempty_string(self):
        result = summarize_threshold_analysis()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 0)

    def test_matches_3_alpha_sq_is_bool(self):
        result = summarize_threshold_analysis()
        self.assertIsInstance(result["matches_3_alpha_sq"], bool)

    def test_threshold_sum_sub_dict_has_s_log(self):
        result = summarize_threshold_analysis()
        self.assertIn("S_log", result["threshold_sum"])

    def test_correction_to_g_sub_dict_has_target(self):
        result = summarize_threshold_analysis()
        self.assertIn("target", result["correction_to_G"])

    def test_all_sub_dicts_are_dicts(self):
        result = summarize_threshold_analysis()
        for key in ("threshold_sum", "correction_to_G", "log_structure",
                     "matching_scales"):
            self.assertIsInstance(result[key], dict,
                                 f"{key} should be a dict")


# ---------------------------------------------------------------------------
# TestLogSumVsPowerSum
# ---------------------------------------------------------------------------
class TestLogSumVsPowerSum(unittest.TestCase):
    """Compare S_log vs the closed-form power sum S_power = L(L+1)^2(L+2)/2.

    The log sum includes an extra ln[l(l+1)] factor, so it grows faster
    than the power sum for large L.
    """

    @staticmethod
    def _power_sum(L):
        """Closed-form: sum_{l=1}^{L} (2l+1)*l*(l+1) = L*(L+1)^2*(L+2)/2."""
        return L * (L + 1) ** 2 * (L + 2) / 2

    @staticmethod
    def _manual_s_log(L):
        """Direct computation of S_log for validation."""
        return sum(
            (2 * l + 1) * l * (l + 1) * math.log(l * (l + 1))
            for l in range(1, L + 1)
        )

    def test_s_log_greater_than_power_sum_for_large_l(self):
        """S_log(L) > S_power(L) for L > 1, due to ln factor > 1."""
        for L in [5, 10, 50]:
            result = compute_threshold_sum(L_max=L)
            s_log = result["S_log"]
            s_power = self._power_sum(L)
            self.assertGreater(
                s_log, s_power,
                f"S_log({L}) = {s_log} should be > S_power({L}) = {s_power}"
            )

    def test_ratio_s_log_over_power_sum_increases(self):
        """The ratio S_log/S_power should increase with L because ln grows."""
        ratios = []
        for L in [10, 50, 100]:
            result = compute_threshold_sum(L_max=L)
            s_log = result["S_log"]
            s_power = self._power_sum(L)
            ratios.append(s_log / s_power)
        for i in range(len(ratios) - 1):
            self.assertGreater(ratios[i + 1], ratios[i],
                               "S_log/S_power ratio should increase with L")

    def test_manual_s_log_matches_module(self):
        """Cross-check module S_log against manual computation for L=10."""
        result = compute_threshold_sum(L_max=10)
        manual = self._manual_s_log(10)
        self.assertAlmostEqual(result["S_log"], manual, places=4,
                               msg="Module S_log should match manual sum")


if __name__ == "__main__":
    unittest.main()
