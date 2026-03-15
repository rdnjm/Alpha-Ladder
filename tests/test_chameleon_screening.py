"""
Tests for alpha_ladder_core.chameleon_screening.

Covers all 7 public functions with ~45 focused tests across 7 test classes.
"""

import unittest

from alpha_ladder_core.chameleon_screening import (
    compute_effective_potential_with_matter,
    find_density_dependent_minimum,
    compute_chameleon_profile,
    check_kk_truncation_validity,
    assess_chameleon_self_consistency,
    compute_meV_dark_sector_potential,
    summarize_chameleon_screening,
)


class TestEffectivePotentialWithMatter(unittest.TestCase):
    """Tests for compute_effective_potential_with_matter."""

    def test_returns_dict(self):
        result = compute_effective_potential_with_matter()
        self.assertIsInstance(result, dict)

    def test_default_grid_length(self):
        result = compute_effective_potential_with_matter()
        self.assertEqual(len(result["sigma_grid"]), 81)
        self.assertEqual(len(result["V_flux"]), 81)
        self.assertEqual(len(result["V_matter"]), 81)
        self.assertEqual(len(result["V_total"]), 81)

    def test_zero_density_matches_flux(self):
        result = compute_effective_potential_with_matter(rho=0.0)
        for v_m in result["V_matter"]:
            self.assertEqual(v_m, 0.0)
        for v_f, v_t in zip(result["V_flux"], result["V_total"]):
            self.assertAlmostEqual(v_f, v_t, places=15)

    def test_matter_term_positive(self):
        result = compute_effective_potential_with_matter(rho=1e-10)
        # For sigma > -350 (where exp(2*sigma) is nonzero), V_matter >= 0
        for sig, v_m in zip(result["sigma_grid"], result["V_matter"]):
            if sig > -350:
                self.assertGreaterEqual(v_m, 0.0)

    def test_has_all_keys(self):
        result = compute_effective_potential_with_matter()
        expected_keys = {
            "sigma_grid", "V_flux", "V_matter", "V_total",
            "A", "B", "C", "rho", "beta_matter",
        }
        self.assertEqual(set(result.keys()), expected_keys)

    def test_custom_grid(self):
        grid = [-5.0, -2.0, 0.0, 1.0, 3.0]
        result = compute_effective_potential_with_matter(sigma_grid=grid)
        self.assertEqual(len(result["sigma_grid"]), 5)
        self.assertEqual(result["sigma_grid"], grid)


class TestDensityDependentMinimum(unittest.TestCase):
    """Tests for find_density_dependent_minimum."""

    def test_returns_dict(self):
        result = find_density_dependent_minimum(rho=0.0)
        self.assertIsInstance(result, dict)

    def test_minimum_exists_at_zero_density(self):
        result = find_density_dependent_minimum(rho=0.0)
        self.assertTrue(result["minimum_exists"])

    def test_b_eff_equals_b_plus_rho(self):
        rho = 1e-50
        beta = 2.5
        result = find_density_dependent_minimum(rho=rho, beta_matter=beta)
        expected_b_eff = result["B"] + beta * rho
        self.assertAlmostEqual(result["B_eff"], expected_b_eff, places=15)

    def test_positive_density_still_has_minimum(self):
        result = find_density_dependent_minimum(rho=1e-80)
        self.assertTrue(result["minimum_exists"])

    def test_has_sigma_0(self):
        result = find_density_dependent_minimum(rho=0.0)
        self.assertIsNotNone(result["sigma_0"])
        self.assertIsInstance(result["sigma_0"], float)

    def test_has_a_effective(self):
        result = find_density_dependent_minimum(rho=0.0)
        self.assertIsNotNone(result["a_effective"])
        self.assertGreater(result["a_effective"], 0.0)

    def test_discriminant_positive(self):
        result = find_density_dependent_minimum(rho=0.0)
        self.assertGreater(result["discriminant"], 0.0)


class TestChameleonProfile(unittest.TestCase):
    """Tests for compute_chameleon_profile."""

    def test_returns_dict(self):
        result = compute_chameleon_profile()
        self.assertIsInstance(result, dict)

    def test_default_48_points(self):
        result = compute_chameleon_profile()
        self.assertEqual(len(result["rho_values_si"]), 48)

    def test_rho_values_si_matches_length(self):
        result = compute_chameleon_profile()
        self.assertEqual(
            len(result["rho_values_si"]),
            len(result["rho_values_planck"]),
        )

    def test_rho_values_planck_matches_length(self):
        result = compute_chameleon_profile()
        self.assertEqual(
            len(result["rho_values_planck"]),
            len(result["sigma_0_values"]),
        )

    def test_all_lists_same_length(self):
        result = compute_chameleon_profile()
        n = len(result["rho_values_si"])
        self.assertEqual(len(result["rho_values_planck"]), n)
        self.assertEqual(len(result["sigma_0_values"]), n)
        self.assertEqual(len(result["a_eff_values"]), n)
        self.assertEqual(len(result["m_phi_eff_values"]), n)
        self.assertEqual(len(result["m_phi_eff_eV_values"]), n)

    def test_sigma_0_values_are_numeric(self):
        result = compute_chameleon_profile()
        numeric_count = sum(
            1 for v in result["sigma_0_values"] if v is not None
        )
        self.assertGreater(numeric_count, 0)

    def test_custom_rho_values(self):
        rho_vals = [1e-20, 1e-10, 1.0, 1e10]
        result = compute_chameleon_profile(rho_values=rho_vals)
        self.assertEqual(len(result["rho_values_si"]), 4)
        self.assertEqual(result["rho_values_si"], rho_vals)


class TestKKTruncationValidity(unittest.TestCase):
    """Tests for check_kk_truncation_validity."""

    def test_returns_dict(self):
        result = check_kk_truncation_validity(1e-10)
        self.assertIsInstance(result, dict)

    def test_planck_length_valid(self):
        result = check_kk_truncation_validity(1.6e-35)
        self.assertEqual(result["exclusion_status"], "valid")
        self.assertTrue(result["kk_truncation_valid"])

    def test_30um_valid(self):
        result = check_kk_truncation_validity(30e-6)
        self.assertEqual(result["exclusion_status"], "valid")
        self.assertTrue(result["kk_truncation_valid"])

    def test_100um_marginal_or_excluded(self):
        result = check_kk_truncation_validity(100e-6)
        self.assertIn(result["exclusion_status"], ("marginal", "excluded"))
        self.assertFalse(result["kk_truncation_valid"])

    def test_1m_excluded(self):
        result = check_kk_truncation_validity(1.0)
        self.assertEqual(result["exclusion_status"], "excluded")
        self.assertFalse(result["kk_truncation_valid"])

    def test_10000km_excluded(self):
        result = check_kk_truncation_validity(1e7)
        self.assertEqual(result["exclusion_status"], "excluded")
        self.assertFalse(result["kk_truncation_valid"])
        self.assertIn("catastrophically", result["exclusion_reason"])

    def test_m_kk_first_positive(self):
        result = check_kk_truncation_validity(30e-6)
        self.assertGreater(result["m_kk_first_eV"], 0.0)

    def test_gravity_6d_threshold_equals_a_eff(self):
        a_eff = 42e-6
        result = check_kk_truncation_validity(a_eff)
        self.assertEqual(result["gravity_becomes_6d_below_m"], a_eff)


class TestSelfConsistency(unittest.TestCase):
    """Tests for assess_chameleon_self_consistency."""

    def test_returns_dict(self):
        result = assess_chameleon_self_consistency()
        self.assertIsInstance(result, dict)

    def test_has_four_environments(self):
        result = assess_chameleon_self_consistency()
        self.assertEqual(len(result["environment_results"]), 4)

    def test_environment_names(self):
        result = assess_chameleon_self_consistency()
        names = [e["name"] for e in result["environment_results"]]
        self.assertIn("cosmic_void", names)
        self.assertIn("interplanetary", names)
        self.assertIn("laboratory", names)
        self.assertIn("earth_surface", names)

    def test_fuzzy_dm_not_viable(self):
        result = assess_chameleon_self_consistency()
        self.assertFalse(result["fuzzy_dm_viable"])

    def test_has_overall_verdict(self):
        result = assess_chameleon_self_consistency()
        self.assertIsInstance(result["overall_verdict"], str)
        self.assertGreater(len(result["overall_verdict"]), 0)

    def test_has_no_go_reason(self):
        result = assess_chameleon_self_consistency()
        self.assertIsInstance(result["no_go_reason"], str)
        self.assertGreater(len(result["no_go_reason"]), 0)

    def test_mev_regime_viable_is_bool(self):
        result = assess_chameleon_self_consistency()
        self.assertIsInstance(result["mev_regime_viable"], bool)


class TestMeVDarkSector(unittest.TestCase):
    """Tests for compute_meV_dark_sector_potential."""

    def test_returns_dict(self):
        result = compute_meV_dark_sector_potential()
        self.assertIsInstance(result, dict)

    def test_default_a0_30um(self):
        result = compute_meV_dark_sector_potential()
        self.assertAlmostEqual(result["a_0_m"], 30e-6, places=10)

    def test_m_phi_order_of_magnitude(self):
        result = compute_meV_dark_sector_potential(a_0_m=30e-6)
        self.assertGreater(result["m_phi_eV"], 1e-3)
        self.assertLess(result["m_phi_eV"], 1e-1)

    def test_not_fuzzy_dm(self):
        result = compute_meV_dark_sector_potential()
        self.assertFalse(result["is_fuzzy_dm_candidate"])

    def test_not_quintessence(self):
        result = compute_meV_dark_sector_potential()
        self.assertFalse(result["is_quintessence_candidate"])

    def test_is_fifth_force(self):
        result = compute_meV_dark_sector_potential()
        self.assertTrue(result["is_testable_fifth_force"])

    def test_has_honest_assessment(self):
        result = compute_meV_dark_sector_potential()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 0)


class TestSummarize(unittest.TestCase):
    """Tests for summarize_chameleon_screening."""

    def test_returns_dict(self):
        result = summarize_chameleon_screening()
        self.assertIsInstance(result, dict)

    def test_has_chameleon_profile(self):
        result = summarize_chameleon_screening()
        self.assertIn("chameleon_profile", result)
        self.assertIsInstance(result["chameleon_profile"], dict)

    def test_has_self_consistency(self):
        result = summarize_chameleon_screening()
        self.assertIn("self_consistency", result)
        self.assertIsInstance(result["self_consistency"], dict)

    def test_has_mev_analysis(self):
        result = summarize_chameleon_screening()
        self.assertIn("mev_analysis", result)
        self.assertIsInstance(result["mev_analysis"], dict)

    def test_fuzzy_dm_excluded(self):
        result = summarize_chameleon_screening()
        self.assertTrue(result["fuzzy_dm_excluded"])

    def test_has_key_finding(self):
        result = summarize_chameleon_screening()
        self.assertIsInstance(result["key_finding"], str)
        self.assertGreater(len(result["key_finding"]), 0)

    def test_has_mev_regime_status(self):
        result = summarize_chameleon_screening()
        self.assertIn("mev_regime_status", result)
        self.assertIsInstance(result["mev_regime_status"], str)
        self.assertGreater(len(result["mev_regime_status"]), 0)


if __name__ == "__main__":
    unittest.main()
