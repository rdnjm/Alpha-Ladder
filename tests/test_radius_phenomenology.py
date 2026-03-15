"""Tests for alpha_ladder_core.radius_phenomenology module.

Validates dilaton mass vs radius mapping, screening amplitude,
experimental bounds, and the testable window for extra dimensions.
"""

import math
import unittest

from alpha_ladder_core.radius_phenomenology import (
    compute_dilaton_mass_vs_radius,
    compute_screening_amplitude_vs_radius,
    compute_experimental_bounds,
    compute_testable_window,
    compute_phenomenology_at_radius,
    summarize_radius_phenomenology,
)


class TestDilatonMassVsRadius(unittest.TestCase):
    """Tests for compute_dilaton_mass_vs_radius."""

    def test_returns_dict(self):
        result = compute_dilaton_mass_vs_radius()
        self.assertIsInstance(result, dict)

    def test_a0_values_is_list(self):
        result = compute_dilaton_mass_vs_radius()
        self.assertIn("a_0_values", result)
        self.assertIsInstance(result["a_0_values"], list)
        self.assertGreater(len(result["a_0_values"]), 1)

    def test_m_phi_same_length_as_a0(self):
        result = compute_dilaton_mass_vs_radius()
        self.assertEqual(len(result["m_phi_eV"]),
                         len(result["a_0_values"]))

    def test_mass_decreases_with_radius(self):
        result = compute_dilaton_mass_vs_radius()
        a0 = result["a_0_values"]
        mass = result["m_phi_eV"]
        # Find pairs where a_0 increases and verify mass decreases
        for i in range(len(a0) - 1):
            if a0[i] < a0[i + 1]:
                self.assertGreater(mass[i], mass[i + 1],
                                   "Mass should decrease as radius increases")

    def test_planck_radius_gives_planck_mass(self):
        result = compute_dilaton_mass_vs_radius()
        l_pl = 1.616e-35
        m_pl_eV = 1.22e28
        # Find entry closest to Planck length
        a0 = result["a_0_values"]
        mass = result["m_phi_eV"]
        idx = min(range(len(a0)), key=lambda i: abs(a0[i] - l_pl))
        ratio = mass[idx] / m_pl_eV
        self.assertGreater(ratio, 0.1, "At Planck radius, mass ~ M_Pl")
        self.assertLess(ratio, 10, "At Planck radius, mass ~ M_Pl")

    def test_sub_mm_radius_gives_mev_mass(self):
        result = compute_dilaton_mass_vs_radius()
        a0_target = 1e-4  # 0.1 mm
        expected_mass = 2e-3  # ~2 meV
        a0 = result["a_0_values"]
        mass = result["m_phi_eV"]
        idx = min(range(len(a0)), key=lambda i: abs(a0[i] - a0_target))
        ratio = mass[idx] / expected_mass
        self.assertGreater(ratio, 0.1,
                           "At 0.1 mm, mass ~ 2 meV within order of mag")
        self.assertLess(ratio, 10,
                        "At 0.1 mm, mass ~ 2 meV within order of mag")

    def test_lambda_compton_same_length(self):
        result = compute_dilaton_mass_vs_radius()
        self.assertIn("lambda_compton_m", result)
        self.assertEqual(len(result["lambda_compton_m"]),
                         len(result["a_0_values"]))


class TestScreeningAmplitude(unittest.TestCase):
    """Tests for compute_screening_amplitude_vs_radius."""

    def test_returns_dict(self):
        result = compute_screening_amplitude_vs_radius()
        self.assertIsInstance(result, dict)

    def test_has_alpha_eff(self):
        result = compute_screening_amplitude_vs_radius()
        self.assertIn("alpha_eff", result)
        self.assertIsInstance(result["alpha_eff"], list)

    def test_alpha_fifth_value(self):
        result = compute_screening_amplitude_vs_radius()
        self.assertAlmostEqual(result["alpha_fifth"], 0.618, places=2)

    def test_planck_radius_fully_screened(self):
        result = compute_screening_amplitude_vs_radius()
        a0 = result["a_0_values"]
        alpha = result["alpha_eff"]
        # At Planck scale, dilaton is super-massive -> fully screened
        idx = min(range(len(a0)), key=lambda i: abs(a0[i] - 1.6e-35))
        self.assertLess(alpha[idx], 1e-10,
                        "At Planck radius, alpha_eff ~ 0 (fully screened)")

    def test_large_radius_less_screened(self):
        result = compute_screening_amplitude_vs_radius()
        a0 = result["a_0_values"]
        alpha = result["alpha_eff"]
        # At large radius (1 mm), screening is weaker than at Planck scale
        # But at default r_test=0.1m, even 1mm Compton still gives some suppression
        idx_large = min(range(len(a0)), key=lambda i: abs(a0[i] - 1e-3))
        idx_small = min(range(len(a0)), key=lambda i: abs(a0[i] - 1e-30))
        self.assertGreaterEqual(alpha[idx_large], alpha[idx_small])

    def test_has_r_test(self):
        result = compute_screening_amplitude_vs_radius()
        self.assertIn("r_test", result)


class TestExperimentalBounds(unittest.TestCase):
    """Tests for compute_experimental_bounds."""

    def test_returns_dict(self):
        result = compute_experimental_bounds()
        self.assertIsInstance(result, dict)

    def test_has_bounds(self):
        result = compute_experimental_bounds()
        self.assertIn("bounds", result)
        self.assertIsInstance(result["bounds"], dict)

    def test_at_least_three_experiments(self):
        result = compute_experimental_bounds()
        self.assertGreaterEqual(len(result["bounds"]), 3)

    def test_eot_wash_entry(self):
        result = compute_experimental_bounds()
        # Key may be eot_wash or eot_wash_2006
        has_eot = any("eot_wash" in k for k in result["bounds"])
        self.assertTrue(has_eot, "Should have an Eot-Wash entry")
        eot_key = [k for k in result["bounds"] if "eot_wash" in k][0]
        self.assertGreater(result["bounds"][eot_key]["lambda_min_m"], 0)

    def test_cassini_entry(self):
        result = compute_experimental_bounds()
        self.assertIn("cassini", result["bounds"])
        cassini = result["bounds"]["cassini"]
        self.assertAlmostEqual(cassini["alpha_max"], 2.3e-5, places=8)


class TestTestableWindow(unittest.TestCase):
    """Tests for compute_testable_window."""

    def test_returns_dict(self):
        result = compute_testable_window()
        self.assertIsInstance(result, dict)

    def test_window_exists(self):
        result = compute_testable_window()
        self.assertTrue(result["window_exists"])

    def test_a0_min_near_planck(self):
        result = compute_testable_window()
        l_pl = 1.616e-35
        ratio = result["a_0_min"] / l_pl
        self.assertGreater(ratio, 0.1)
        self.assertLess(ratio, 10)

    def test_a0_max_near_sub_mm(self):
        result = compute_testable_window()
        target = 1e-4  # 0.1 mm
        ratio = result["a_0_max_testable"] / target
        self.assertGreater(ratio, 0.5, "a_0_max within factor of 2 of 0.1 mm")
        self.assertLess(ratio, 2.0, "a_0_max within factor of 2 of 0.1 mm")

    def test_max_greater_than_min(self):
        result = compute_testable_window()
        self.assertGreater(result["a_0_max_testable"], result["a_0_min"])

    def test_has_experiments_in_window(self):
        result = compute_testable_window()
        self.assertIn("experiments_in_window", result)
        self.assertIsInstance(result["experiments_in_window"], list)

    def test_window_description_nonempty(self):
        result = compute_testable_window()
        self.assertIn("window_description", result)
        self.assertIsInstance(result["window_description"], str)
        self.assertGreater(len(result["window_description"]), 0)


class TestPhenomenologyAtRadius(unittest.TestCase):
    """Tests for compute_phenomenology_at_radius."""

    def test_planck_radius_invisible(self):
        result = compute_phenomenology_at_radius(a_0=1e-35)
        self.assertIsInstance(result, dict)
        self.assertEqual(result["classification"], "invisible")
        self.assertTrue(result["passes_cassini"])
        self.assertTrue(result["passes_eot_wash"])

    def test_50um_radius_sub_mm(self):
        result = compute_phenomenology_at_radius(a_0=50e-6)
        self.assertIsInstance(result, dict)
        self.assertIn(result["classification"], ["sub_mm", "sub-mm", "marginal"])

    def test_01mm_radius_marginal(self):
        result = compute_phenomenology_at_radius(a_0=1e-4)
        self.assertIsInstance(result, dict)
        self.assertEqual(result["classification"], "marginal")

    def test_1cm_radius_excluded(self):
        result = compute_phenomenology_at_radius(a_0=1e-2)
        self.assertIsInstance(result, dict)
        self.assertEqual(result["classification"], "excluded")

    def test_has_mass_key(self):
        result = compute_phenomenology_at_radius(a_0=1e-35)
        self.assertIn("m_phi_eV", result)

    def test_has_passes_cassini(self):
        result = compute_phenomenology_at_radius(a_0=1e-35)
        self.assertIn("passes_cassini", result)

    def test_has_passes_eot_wash(self):
        result = compute_phenomenology_at_radius(a_0=1e-35)
        self.assertIn("passes_eot_wash", result)

    def test_has_screening_at_1au(self):
        result = compute_phenomenology_at_radius(a_0=1e-35)
        # Key may be screening_at_1AU or alpha_at_1AU
        has_screening = "screening_at_1AU" in result or "alpha_at_1AU" in result
        self.assertTrue(has_screening)


class TestSummarize(unittest.TestCase):
    """Tests for summarize_radius_phenomenology."""

    def test_returns_dict(self):
        result = summarize_radius_phenomenology()
        self.assertIsInstance(result, dict)

    def test_has_sub_dicts(self):
        result = summarize_radius_phenomenology()
        for key in ["mass_curve", "testable_window",
                     "key_radii", "experimental_bounds"]:
            self.assertIn(key, result)

    def test_has_overall_assessment(self):
        result = summarize_radius_phenomenology()
        self.assertIn("overall_assessment", result)
        self.assertIsInstance(result["overall_assessment"], str)

    def test_has_framework_prediction(self):
        result = summarize_radius_phenomenology()
        self.assertIn("framework_prediction", result)
        self.assertIsInstance(result["framework_prediction"], str)

    def test_key_radii_entries(self):
        result = summarize_radius_phenomenology()
        kr = result["key_radii"]
        self.assertIn("planck", kr)
        has_eot = any("eot_wash" in k for k in kr)
        self.assertTrue(has_eot, "Should have an eot_wash-related key in key_radii")


if __name__ == "__main__":
    unittest.main()
