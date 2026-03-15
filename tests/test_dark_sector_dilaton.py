"""
Tests for the NEW dilaton phenomenology functions in alpha_ladder_core.dark_sector.

These tests cover the six new functions added for dilaton dark sector analysis:
  - compute_relic_abundance
  - compute_self_interaction
  - compute_equation_of_state
  - compute_fuzzy_dm_constraints
  - compute_dark_sector_landscape
  - summarize_dark_sector

The original backward-compatible functions (compute_dark_sector, compute_wave_profile,
compute_alps_simulation, get_experimental_bounds) are tested in test_dark_sector.py.
"""

import math
import unittest

from alpha_ladder_core.dark_sector import (
    compute_relic_abundance,
    compute_self_interaction,
    compute_equation_of_state,
    compute_fuzzy_dm_constraints,
    compute_dark_sector_landscape,
    summarize_dark_sector,
    M_Pl_eV,
)


# Default meV dilaton mass used throughout tests
M_MEV = 7e-3  # 7 meV in eV


class TestRelicAbundance(unittest.TestCase):
    """Tests for compute_relic_abundance (misalignment mechanism)."""

    def test_returns_dict(self):
        result = compute_relic_abundance(M_MEV)
        self.assertIsInstance(result, dict)
        expected_keys = {
            "m_phi_eV", "f_initial_eV", "Omega_phi_h2", "overproduced",
            "ratio_to_observed", "f_required_for_dm_eV",
            "f_required_over_M_Pl", "fine_tuning_required",
            "honest_assessment",
        }
        self.assertEqual(set(result.keys()), expected_keys)

    def test_mev_overproduces(self):
        """meV dilaton with f=M_Pl massively overproduces dark matter."""
        result = compute_relic_abundance(M_MEV)
        self.assertTrue(result["overproduced"])

    def test_ratio_much_greater_than_one(self):
        """Omega_phi / Omega_DM >> 1 for meV dilaton with f=M_Pl."""
        result = compute_relic_abundance(M_MEV)
        self.assertGreater(result["ratio_to_observed"], 1e6)

    def test_f_required_less_than_mpl(self):
        """f_required / M_Pl << 1 for meV dilaton -- needs fine tuning."""
        result = compute_relic_abundance(M_MEV)
        self.assertLess(result["f_required_over_M_Pl"], 1.0)

    def test_fine_tuning_required(self):
        """meV dilaton requires severe fine tuning (f/M_Pl < 0.01)."""
        result = compute_relic_abundance(M_MEV)
        self.assertTrue(result["fine_tuning_required"])

    def test_has_honest_assessment(self):
        result = compute_relic_abundance(M_MEV)
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 50)

    def test_custom_f_initial(self):
        """Custom f_initial should be used instead of M_Pl default."""
        f_custom = 1e20  # much less than M_Pl
        result = compute_relic_abundance(M_MEV, f_initial_eV=f_custom)
        self.assertEqual(result["f_initial_eV"], f_custom)
        # With a smaller f, Omega should be smaller than with f=M_Pl
        result_default = compute_relic_abundance(M_MEV)
        self.assertLess(result["Omega_phi_h2"], result_default["Omega_phi_h2"])


class TestSelfInteraction(unittest.TestCase):
    """Tests for compute_self_interaction (bullet cluster and MW bounds)."""

    def test_returns_dict(self):
        result = compute_self_interaction(M_MEV)
        self.assertIsInstance(result, dict)
        expected_keys = {
            "m_phi_eV", "alpha_coupling", "sigma_cm2",
            "sigma_over_m_cm2_g", "bullet_cluster_bound",
            "passes_bullet_cluster", "milky_way_bound",
            "passes_milky_way", "honest_assessment",
        }
        self.assertEqual(set(result.keys()), expected_keys)

    def test_sigma_positive(self):
        result = compute_self_interaction(M_MEV)
        self.assertGreater(result["sigma_cm2"], 0)

    def test_sigma_over_m_positive(self):
        result = compute_self_interaction(M_MEV)
        self.assertGreater(result["sigma_over_m_cm2_g"], 0)

    def test_mev_sigma_over_m_finite(self):
        """meV dilaton has a finite sigma/m value."""
        result = compute_self_interaction(M_MEV)
        self.assertNotEqual(result["sigma_over_m_cm2_g"], float('inf'))
        self.assertIsInstance(result["passes_bullet_cluster"], bool)

    def test_mev_passes_consistency(self):
        """passes_milky_way implies passes_bullet_cluster (tighter bound)."""
        result = compute_self_interaction(M_MEV)
        if result["passes_milky_way"]:
            self.assertTrue(result["passes_bullet_cluster"])

    def test_bounds_are_correct(self):
        """Verify the numerical values of the observational bounds."""
        result = compute_self_interaction(M_MEV)
        self.assertEqual(result["bullet_cluster_bound"], 1.0)
        self.assertEqual(result["milky_way_bound"], 0.1)

    def test_has_honest_assessment(self):
        result = compute_self_interaction(M_MEV)
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 50)


class TestEquationOfState(unittest.TestCase):
    """Tests for compute_equation_of_state (w(z) for dilaton field)."""

    def test_returns_dict(self):
        result = compute_equation_of_state(M_MEV)
        self.assertIsInstance(result, dict)
        expected_keys = {
            "m_phi_eV", "z_values", "w_values", "z_transition",
            "H_at_transition_eV", "w_today", "behaves_as_matter",
            "behaves_as_dark_energy", "honest_assessment",
        }
        self.assertEqual(set(result.keys()), expected_keys)

    def test_mev_w_today_zero(self):
        """meV dilaton has w_today ~ 0 (oscillating matter, not dark energy)."""
        result = compute_equation_of_state(M_MEV)
        self.assertAlmostEqual(result["w_today"], 0.0, places=5)

    def test_mev_behaves_as_matter(self):
        """meV dilaton began oscillating at very high z -- behaves as matter."""
        result = compute_equation_of_state(M_MEV)
        self.assertTrue(result["behaves_as_matter"])

    def test_mev_not_dark_energy(self):
        """meV dilaton is NOT dark energy."""
        result = compute_equation_of_state(M_MEV)
        self.assertFalse(result["behaves_as_dark_energy"])

    def test_z_transition_large(self):
        """meV dilaton transition redshift is very large (> 10^6)."""
        result = compute_equation_of_state(M_MEV)
        self.assertGreater(result["z_transition"], 1e6)

    def test_default_z_values_length(self):
        """Default z_values list has 13 entries."""
        result = compute_equation_of_state(M_MEV)
        self.assertEqual(len(result["z_values"]), 13)

    def test_w_values_same_length_as_z(self):
        """w_values and z_values must have the same length."""
        result = compute_equation_of_state(M_MEV)
        self.assertEqual(len(result["w_values"]), len(result["z_values"]))

    def test_has_honest_assessment(self):
        result = compute_equation_of_state(M_MEV)
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 50)


class TestFuzzyDMConstraints(unittest.TestCase):
    """Tests for compute_fuzzy_dm_constraints (ultralight DM bounds)."""

    def test_returns_dict(self):
        result = compute_fuzzy_dm_constraints(M_MEV)
        self.assertIsInstance(result, dict)
        expected_keys = {
            "m_phi_eV", "constraints", "passes_all_bounds",
            "de_broglie_wavelength_m", "de_broglie_wavelength_kpc",
            "is_fuzzy_candidate", "kk_truncation_valid",
            "honest_assessment",
        }
        self.assertEqual(set(result.keys()), expected_keys)

    def test_mev_passes_all_bounds(self):
        """meV dilaton trivially passes all fuzzy DM mass bounds (it is heavier)."""
        result = compute_fuzzy_dm_constraints(M_MEV)
        self.assertTrue(result["passes_all_bounds"])

    def test_mev_not_fuzzy_candidate(self):
        """meV dilaton is NOT a fuzzy DM candidate (de Broglie wavelength too small)."""
        result = compute_fuzzy_dm_constraints(M_MEV)
        self.assertFalse(result["is_fuzzy_candidate"])

    def test_mev_kk_valid(self):
        """meV dilaton has kk_truncation_valid = True (a_0 < 56 um)."""
        result = compute_fuzzy_dm_constraints(M_MEV)
        self.assertTrue(result["kk_truncation_valid"])

    def test_has_4_constraints(self):
        """There are exactly 4 observational constraints."""
        result = compute_fuzzy_dm_constraints(M_MEV)
        self.assertEqual(len(result["constraints"]), 4)

    def test_de_broglie_sub_kpc(self):
        """meV dilaton de Broglie wavelength is far below kpc scale (not fuzzy)."""
        result = compute_fuzzy_dm_constraints(M_MEV)
        self.assertLess(result["de_broglie_wavelength_kpc"], 0.1)

    def test_has_honest_assessment(self):
        result = compute_fuzzy_dm_constraints(M_MEV)
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 50)


class TestDarkSectorLandscape(unittest.TestCase):
    """Tests for compute_dark_sector_landscape (mass-scale survey)."""

    def test_returns_dict(self):
        result = compute_dark_sector_landscape()
        self.assertIsInstance(result, dict)
        expected_keys = {
            "mass_points", "landscape", "viable_window",
            "honest_assessment",
        }
        self.assertEqual(set(result.keys()), expected_keys)

    def test_default_11_mass_points(self):
        """Default mass_points list has 11 entries."""
        result = compute_dark_sector_landscape()
        self.assertEqual(len(result["mass_points"]), 11)

    def test_landscape_same_length_as_mass_points(self):
        """Landscape list matches mass_points in length."""
        result = compute_dark_sector_landscape()
        self.assertEqual(len(result["landscape"]), len(result["mass_points"]))

    def test_planck_mass_invisible(self):
        """Planck-mass entry (1e28 eV) is classified as invisible_planck."""
        result = compute_dark_sector_landscape()
        planck_entry = [e for e in result["landscape"] if e["m_phi_eV"] == 1e28]
        self.assertEqual(len(planck_entry), 1)
        self.assertEqual(planck_entry[0]["classification"], "invisible_planck")

    def test_has_viable_window(self):
        """Result contains a viable_window tuple."""
        result = compute_dark_sector_landscape()
        self.assertIsNotNone(result["viable_window"])
        self.assertEqual(len(result["viable_window"]), 2)

    def test_each_entry_has_classification(self):
        """Every landscape entry must have a classification string."""
        result = compute_dark_sector_landscape()
        for entry in result["landscape"]:
            self.assertIn("classification", entry)
            self.assertIsInstance(entry["classification"], str)

    def test_custom_mass_points(self):
        """Custom mass_points list is used correctly."""
        custom = [1e-3, 1e0, 1e10]
        result = compute_dark_sector_landscape(mass_points=custom)
        self.assertEqual(len(result["landscape"]), 3)
        self.assertEqual(result["mass_points"], custom)

    def test_has_honest_assessment(self):
        result = compute_dark_sector_landscape()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 50)


class TestSummarizeDarkSector(unittest.TestCase):
    """Tests for summarize_dark_sector (dashboard entry point)."""

    @classmethod
    def setUpClass(cls):
        cls.result = summarize_dark_sector()

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)
        expected_keys = {
            "mev_analysis", "relic_abundance", "self_interaction",
            "equation_of_state", "fuzzy_constraints", "landscape",
            "key_finding", "overall_verdict", "framework_position",
        }
        self.assertEqual(set(self.result.keys()), expected_keys)

    def test_has_mev_analysis(self):
        mev = self.result["mev_analysis"]
        self.assertIsInstance(mev, dict)
        self.assertIn("m_phi_eV", mev)
        self.assertIn("a_0_m", mev)

    def test_mev_not_dark_matter(self):
        self.assertFalse(self.result["mev_analysis"]["is_dark_matter"])

    def test_mev_not_dark_energy(self):
        self.assertFalse(self.result["mev_analysis"]["is_dark_energy"])

    def test_mev_is_fifth_force(self):
        self.assertTrue(self.result["mev_analysis"]["is_fifth_force"])

    def test_has_key_finding(self):
        self.assertIsInstance(self.result["key_finding"], str)
        self.assertGreater(len(self.result["key_finding"]), 50)

    def test_has_overall_verdict(self):
        self.assertIsInstance(self.result["overall_verdict"], str)
        self.assertGreater(len(self.result["overall_verdict"]), 50)

    def test_has_framework_position(self):
        self.assertIsInstance(self.result["framework_position"], str)
        self.assertGreater(len(self.result["framework_position"]), 50)


if __name__ == "__main__":
    unittest.main()
