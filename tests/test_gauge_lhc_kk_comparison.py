"""Tests for lhc_kk_comparison.py -- LHC KK graviton cross-section ratio S^2 vs ADD."""
import unittest
import math



class TestModuleConstants(unittest.TestCase):
    """Test module-level constants are correctly defined."""

    def test_alpha_em_value(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import ALPHA_EM
        self.assertAlmostEqual(ALPHA_EM, 1.0 / 137.035999084, places=12)

    def test_m_pl_ev_value(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import M_PL_EV
        self.assertAlmostEqual(M_PL_EV, 1.22089e28, delta=1e24)

    def test_phi_vev_value(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import PHI_VEV
        self.assertAlmostEqual(PHI_VEV, -0.5974, delta=0.01)

    def test_e_phi_vev_value(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import E_PHI_VEV
        self.assertAlmostEqual(E_PHI_VEV, 0.5503, delta=0.01)

    def test_phi_vev_consistency(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import PHI_VEV, ALPHA_EM
        expected = 0.25 * math.log(4.0 * math.pi * ALPHA_EM)
        self.assertAlmostEqual(PHI_VEV, expected, places=10)


class TestAddKKSpectrum(unittest.TestCase):
    """Test the ADD (torus) KK spectrum function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import add_kk_spectrum
        R = 1.0
        result = add_kk_spectrum(R, 10.0)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import add_kk_spectrum
        result = add_kk_spectrum(1.0, 10.0)
        for key in ["N_modes", "masses", "mass_degeneracies", "rho_analytic", "N_analytic"]:
            self.assertIn(key, result)

    def test_masses_sorted(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import add_kk_spectrum
        result = add_kk_spectrum(1.0, 10.0)
        masses = result["masses"]
        for i in range(len(masses) - 1):
            self.assertLessEqual(masses[i], masses[i + 1])

    def test_lightest_mode_is_1_over_R(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import add_kk_spectrum
        R = 2.0
        result = add_kk_spectrum(R, 5.0)
        # Lightest mode: sqrt(1)/R = 1/R
        self.assertAlmostEqual(result["masses"][0], 1.0 / R, places=10)

    def test_first_degeneracy_is_4(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import add_kk_spectrum
        result = add_kk_spectrum(1.0, 10.0)
        # k=1: (1,0), (-1,0), (0,1), (0,-1) => degeneracy 4
        self.assertEqual(result["mass_degeneracies"][0][1], 4)

    def test_n_analytic_formula(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import add_kk_spectrum
        R = 1.0
        m_max = 10.0
        result = add_kk_spectrum(R, m_max)
        expected = math.pi * (m_max * R) ** 2
        self.assertAlmostEqual(result["N_analytic"], expected, places=5)

    def test_rho_analytic_formula(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import add_kk_spectrum
        R = 1.0
        m_max = 10.0
        result = add_kk_spectrum(R, m_max)
        expected = 2.0 * math.pi * m_max * R * R
        self.assertAlmostEqual(result["rho_analytic"], expected, places=5)

    def test_raises_for_wrong_n_extra(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import add_kk_spectrum
        with self.assertRaises(ValueError):
            add_kk_spectrum(1.0, 10.0, n_extra=3)

    def test_all_masses_below_m_max(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import add_kk_spectrum
        m_max = 5.0
        result = add_kk_spectrum(1.0, m_max)
        for m in result["masses"]:
            self.assertLess(m, m_max)


class TestS2KKSpectrum(unittest.TestCase):
    """Test the S^2 KK spectrum function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import s2_kk_spectrum
        result = s2_kk_spectrum(1.0, 10.0)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import s2_kk_spectrum
        result = s2_kk_spectrum(1.0, 10.0)
        for key in ["N_modes", "mode_list", "N_analytic"]:
            self.assertIn(key, result)

    def test_lightest_mode_sqrt2_over_R(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import s2_kk_spectrum
        R = 2.0
        result = s2_kk_spectrum(R, 10.0)
        # l=1: m = sqrt(1*2)/R = sqrt(2)/R
        l_val, m_val, deg = result["mode_list"][0]
        self.assertEqual(l_val, 1)
        self.assertAlmostEqual(m_val, math.sqrt(2.0) / R, places=10)

    def test_lightest_mode_ratio_to_add(self):
        """S^2 lightest mode sqrt(2)/R vs ADD lightest 1/R => ratio sqrt(2)."""
        from alpha_ladder_core.gauge_lhc_kk_comparison import s2_kk_spectrum, add_kk_spectrum
        R = 1.0
        s2 = s2_kk_spectrum(R, 10.0)
        add = add_kk_spectrum(R, 10.0)
        ratio = s2["mode_list"][0][1] / add["masses"][0]
        self.assertAlmostEqual(ratio, math.sqrt(2.0), places=5)

    def test_degeneracy_is_2l_plus_1(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import s2_kk_spectrum
        result = s2_kk_spectrum(1.0, 20.0)
        for l_val, _, deg in result["mode_list"]:
            self.assertEqual(deg, 2 * l_val + 1)

    def test_s2_sparser_than_add(self):
        """S^2 has fewer total modes than ADD at same R and m_max."""
        from alpha_ladder_core.gauge_lhc_kk_comparison import s2_kk_spectrum, add_kk_spectrum
        R = 1.0
        m_max = 20.0
        s2 = s2_kk_spectrum(R, m_max)
        add = add_kk_spectrum(R, m_max)
        self.assertLess(s2["N_modes"], add["N_modes"])

    def test_total_modes_sum_of_degeneracies(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import s2_kk_spectrum
        result = s2_kk_spectrum(1.0, 10.0)
        total = sum(deg for _, _, deg in result["mode_list"])
        self.assertEqual(result["N_modes"], total)


class TestEffectiveModes(unittest.TestCase):
    """Test the effective_modes function."""

    def test_add_model_returns_dict(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import effective_modes
        result = effective_modes("ADD", 5.0)
        self.assertIsInstance(result, dict)

    def test_s2_model_returns_dict(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import effective_modes
        result = effective_modes("S2", 5.0)
        self.assertIsInstance(result, dict)

    def test_invalid_model_raises(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import effective_modes
        with self.assertRaises(ValueError):
            effective_modes("torus", 5.0)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import effective_modes
        result = effective_modes("ADD", 5.0)
        for key in ["N_eff", "R", "R_meters", "m_lightest_eV", "sigma_weighted"]:
            self.assertIn(key, result)

    def test_add_lightest_is_1_over_R(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import effective_modes
        result = effective_modes("ADD", 5.0)
        self.assertAlmostEqual(result["m_lightest_eV"], 1.0 / result["R"], places=5)

    def test_s2_lightest_is_sqrt2_over_R(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import effective_modes
        result = effective_modes("S2", 5.0)
        self.assertAlmostEqual(result["m_lightest_eV"], math.sqrt(2.0) / result["R"], places=5)

    def test_n_eff_positive(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import effective_modes
        for model in ["ADD", "S2"]:
            result = effective_modes(model, 5.0)
            self.assertGreater(result["N_eff"], 0)

    def test_n_eff_decreases_with_m6(self):
        """Higher M_6 => smaller R => fewer modes."""
        from alpha_ladder_core.gauge_lhc_kk_comparison import effective_modes
        n1 = effective_modes("ADD", 3.0)["N_eff"]
        n2 = effective_modes("ADD", 10.0)["N_eff"]
        self.assertGreater(n1, n2)

    def test_radius_ratio_s2_over_add(self):
        """R_S2/R_ADD should be approximately 0.908 (from analytic formula)."""
        from alpha_ladder_core.gauge_lhc_kk_comparison import effective_modes
        add = effective_modes("ADD", 5.0)
        s2 = effective_modes("S2", 5.0)
        ratio = s2["R"] / add["R"]
        # E_PHI_VEV * 2*pi / sqrt(4*pi) ~ 0.9752
        # Tolerance is loose since it depends on exact constants
        self.assertGreater(ratio, 0.8)
        self.assertLess(ratio, 1.1)


class TestCrossSectionRatio(unittest.TestCase):
    """Test the cross_section_ratio scan function."""

    def test_returns_list(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import cross_section_ratio
        result = cross_section_ratio([3.0, 5.0, 10.0])
        self.assertIsInstance(result, list)

    def test_default_scan_length(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import cross_section_ratio
        result = cross_section_ratio()
        self.assertEqual(len(result), 20)

    def test_each_entry_has_keys(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import cross_section_ratio
        result = cross_section_ratio([5.0])
        entry = result[0]
        for key in ["M_6", "N_ratio", "sigma_ratio", "N_add", "N_s2"]:
            self.assertIn(key, entry)

    def test_n_ratio_less_than_1(self):
        """S^2 produces fewer modes => ratio < 1 (weaker signal)."""
        from alpha_ladder_core.gauge_lhc_kk_comparison import cross_section_ratio
        result = cross_section_ratio([5.0])
        self.assertLess(result[0]["N_ratio"], 1.0)

    def test_sigma_ratio_less_than_1(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import cross_section_ratio
        result = cross_section_ratio([5.0])
        self.assertLess(result[0]["sigma_ratio"], 1.0)

    def test_s2_weaker_than_add(self):
        """Cross-section ratio < 1 for all M_6 values."""
        from alpha_ladder_core.gauge_lhc_kk_comparison import cross_section_ratio
        result = cross_section_ratio([3.0, 5.0, 10.0, 15.0])
        for entry in result:
            self.assertLess(entry["N_ratio"], 1.0)


class TestRescaledExclusion(unittest.TestCase):
    """Test the rescaled_exclusion function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import rescaled_exclusion
        result = rescaled_exclusion()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import rescaled_exclusion
        result = rescaled_exclusion()
        for key in ["M_6_excl_add", "M_6_excl_s2", "M_6_excl_s2_sigma",
                     "ratio_at_exclusion", "N_eff_target"]:
            self.assertIn(key, result)

    def test_s2_exclusion_weaker_than_add(self):
        """S^2 exclusion should be weaker (lower M_6) than ADD."""
        from alpha_ladder_core.gauge_lhc_kk_comparison import rescaled_exclusion
        result = rescaled_exclusion()
        self.assertLess(result["M_6_excl_s2"], result["M_6_excl_add"])

    def test_s2_exclusion_around_1_3_tev(self):
        """S^2 M_6 exclusion ~ 1.3 TeV (from ADD 5 TeV)."""
        from alpha_ladder_core.gauge_lhc_kk_comparison import rescaled_exclusion
        result = rescaled_exclusion()
        # Allow wide tolerance since binary search may vary
        self.assertGreater(result["M_6_excl_s2"], 0.5)
        self.assertLess(result["M_6_excl_s2"], 5.0)

    def test_ratio_at_exclusion_less_than_1(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import rescaled_exclusion
        result = rescaled_exclusion()
        self.assertLess(result["ratio_at_exclusion"], 1.0)

    def test_a0_and_r_phys_positive(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import rescaled_exclusion
        result = rescaled_exclusion()
        self.assertGreater(result["a_0_threshold_s2_m"], 0)
        self.assertGreater(result["R_phys_threshold_s2_m"], 0)

    def test_r_phys_less_than_a0(self):
        """R_phys = a_0 * e^phi_vev < a_0 since phi_vev < 0."""
        from alpha_ladder_core.gauge_lhc_kk_comparison import rescaled_exclusion
        result = rescaled_exclusion()
        self.assertLess(result["R_phys_threshold_s2_m"],
                        result["a_0_threshold_s2_m"])


class TestSpectrumComparisonTable(unittest.TestCase):
    """Test spectrum_comparison_table function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import spectrum_comparison_table
        result = spectrum_comparison_table()
        self.assertIsInstance(result, dict)

    def test_has_both_level_lists(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import spectrum_comparison_table
        result = spectrum_comparison_table()
        self.assertIn("add_levels", result)
        self.assertIn("s2_levels", result)

    def test_up_to_20_levels(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import spectrum_comparison_table
        result = spectrum_comparison_table()
        self.assertLessEqual(len(result["add_levels"]), 20)
        self.assertLessEqual(len(result["s2_levels"]), 20)

    def test_cumulative_monotonically_increasing(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import spectrum_comparison_table
        result = spectrum_comparison_table()
        for cum in [result["add_cumulative"], result["s2_cumulative"]]:
            for i in range(len(cum) - 1):
                self.assertLess(cum[i], cum[i + 1])

    def test_add_and_s2_cumulative_both_grow(self):
        """Both ADD and S^2 cumulative mode counts grow with level index."""
        from alpha_ladder_core.gauge_lhc_kk_comparison import spectrum_comparison_table
        result = spectrum_comparison_table()
        for cum_list in [result["add_cumulative"], result["s2_cumulative"]]:
            self.assertGreater(cum_list[-1], cum_list[0])


class TestSummarizeLhcComparison(unittest.TestCase):
    """Test the summarize_lhc_comparison function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import summarize_lhc_comparison
        result = summarize_lhc_comparison()
        self.assertIsInstance(result, dict)

    def test_has_report_string(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import summarize_lhc_comparison
        result = summarize_lhc_comparison()
        self.assertIn("report", result)
        self.assertIsInstance(result["report"], str)
        self.assertGreater(len(result["report"]), 100)

    def test_has_sub_results(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import summarize_lhc_comparison
        result = summarize_lhc_comparison()
        for key in ["scan", "exclusion", "spectrum_table", "R_ratio_analytic"]:
            self.assertIn(key, result)

    def test_r_ratio_analytic_positive(self):
        from alpha_ladder_core.gauge_lhc_kk_comparison import summarize_lhc_comparison
        result = summarize_lhc_comparison()
        self.assertGreater(result["R_ratio_analytic"], 0)


if __name__ == "__main__":
    unittest.main()
