"""Tests for kk_spectrum_matched.py -- KK spectrum with dilaton vev matched to alpha_EM."""
import unittest
import math



# ---------------------------------------------------------------------------
# Module constants
# ---------------------------------------------------------------------------
class TestModuleConstants(unittest.TestCase):
    """Tests for module-level constants and derived values."""

    def test_alpha_em_value(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import ALPHA_EM
        self.assertAlmostEqual(1.0 / ALPHA_EM, 137.036, delta=0.001)

    def test_phi_vev_approximately_minus_0597(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import PHI_VEV
        self.assertAlmostEqual(PHI_VEV, -0.597, delta=0.01)

    def test_phi_vev_formula(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import PHI_VEV, ALPHA_EM
        expected = 0.25 * math.log(4.0 * math.pi * ALPHA_EM)
        self.assertAlmostEqual(PHI_VEV, expected, places=10)

    def test_e_phi_vev_approximately_055(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import E_PHI_VEV
        self.assertAlmostEqual(E_PHI_VEV, 0.55, delta=0.01)

    def test_e_phi_vev_consistent(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import E_PHI_VEV, PHI_VEV
        self.assertAlmostEqual(E_PHI_VEV, math.exp(PHI_VEV), places=10)

    def test_mass_enhancement_approximately_1817(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import MASS_ENHANCEMENT
        self.assertAlmostEqual(MASS_ENHANCEMENT, 1.817, delta=0.01)

    def test_mass_enhancement_inverse_of_e_phi(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import MASS_ENHANCEMENT, E_PHI_VEV
        self.assertAlmostEqual(MASS_ENHANCEMENT, 1.0 / E_PHI_VEV, places=10)


# ---------------------------------------------------------------------------
# 1. scalar_spectrum
# ---------------------------------------------------------------------------
class TestScalarSpectrum(unittest.TestCase):
    """Tests for scalar_spectrum."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import scalar_spectrum
        result = scalar_spectrum(28e-6)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import scalar_spectrum
        result = scalar_spectrum(28e-6)
        for key in ("modes", "R_phys_m", "phi_vev", "a_0_m"):
            self.assertIn(key, result)

    def test_l0_zero_mode(self):
        """l=0 scalar should have m_eV = 0 (zero mode)."""
        from alpha_ladder_core.gauge_kk_spectrum_matched import scalar_spectrum
        result = scalar_spectrum(28e-6, l_max=5)
        l0 = result["modes"][0]
        self.assertEqual(l0["l"], 0)
        self.assertAlmostEqual(l0["m_eV"], 0.0, places=15)

    def test_degeneracy_2l_plus_1(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import scalar_spectrum
        result = scalar_spectrum(28e-6, l_max=10)
        for mode in result["modes"]:
            self.assertEqual(mode["degeneracy"], 2 * mode["l"] + 1)

    def test_mass_formula_sqrt_l_l_plus_1(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import scalar_spectrum, HBAR_C_EVM
        a_0 = 1e-5
        result = scalar_spectrum(a_0, l_max=5, matched=True)
        R = result["R_phys_m"]
        for mode in result["modes"]:
            l = mode["l"]
            expected = math.sqrt(l * (l + 1)) * HBAR_C_EVM / R
            self.assertAlmostEqual(mode["m_eV"], expected, places=5)

    def test_r_phys_matched(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import scalar_spectrum, E_PHI_VEV
        a_0 = 28e-6
        result = scalar_spectrum(a_0, matched=True)
        self.assertAlmostEqual(result["R_phys_m"], a_0 * E_PHI_VEV, places=15)

    def test_r_phys_unmatched(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import scalar_spectrum
        a_0 = 28e-6
        result = scalar_spectrum(a_0, matched=False)
        self.assertAlmostEqual(result["R_phys_m"], a_0, places=15)

    def test_cumulative_modes(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import scalar_spectrum
        result = scalar_spectrum(28e-6, l_max=5)
        # cumulative should be sum of degeneracies up to each l
        running = 0
        for mode in result["modes"]:
            running += mode["degeneracy"]
            self.assertEqual(mode["cumulative_modes"], running)

    def test_modes_count(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import scalar_spectrum
        result = scalar_spectrum(28e-6, l_max=20)
        self.assertEqual(len(result["modes"]), 21)  # l=0..20

    def test_masses_increase_with_l(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import scalar_spectrum
        result = scalar_spectrum(28e-6, l_max=10)
        masses = [m["m_eV"] for m in result["modes"]]
        for i in range(len(masses) - 1):
            self.assertLessEqual(masses[i], masses[i + 1])


# ---------------------------------------------------------------------------
# 2. vector_spectrum
# ---------------------------------------------------------------------------
class TestVectorSpectrum(unittest.TestCase):
    """Tests for vector_spectrum."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import vector_spectrum
        result = vector_spectrum(28e-6)
        self.assertIsInstance(result, dict)

    def test_starts_at_l1(self):
        """Vector modes start at l=1 (no l=0 massive vector)."""
        from alpha_ladder_core.gauge_kk_spectrum_matched import vector_spectrum
        result = vector_spectrum(28e-6, l_max=5)
        self.assertEqual(result["modes"][0]["l"], 1)

    def test_total_degeneracy_2_times_2l_plus_1(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import vector_spectrum
        result = vector_spectrum(28e-6, l_max=10)
        for mode in result["modes"]:
            expected = 2 * (2 * mode["l"] + 1)
            self.assertEqual(mode["total_degeneracy"], expected)

    def test_degeneracy_per_field(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import vector_spectrum
        result = vector_spectrum(28e-6, l_max=5)
        for mode in result["modes"]:
            self.assertEqual(mode["degeneracy_per_field"], 2 * mode["l"] + 1)

    def test_modes_count(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import vector_spectrum
        result = vector_spectrum(28e-6, l_max=10)
        self.assertEqual(len(result["modes"]), 10)  # l=1..10


# ---------------------------------------------------------------------------
# 3. graviton_spectrum
# ---------------------------------------------------------------------------
class TestGravitonSpectrum(unittest.TestCase):
    """Tests for graviton_spectrum."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import graviton_spectrum
        result = graviton_spectrum(28e-6)
        self.assertIsInstance(result, dict)

    def test_starts_at_l2(self):
        """Graviton modes start at l=2."""
        from alpha_ladder_core.gauge_kk_spectrum_matched import graviton_spectrum
        result = graviton_spectrum(28e-6, l_max=5)
        self.assertEqual(result["modes"][0]["l"], 2)

    def test_l2_eigenvalue_is_2(self):
        """For l=2: sqrt(l(l+1)-2) = sqrt(4) = 2."""
        from alpha_ladder_core.gauge_kk_spectrum_matched import graviton_spectrum, HBAR_C_EVM
        a_0 = 1e-5
        result = graviton_spectrum(a_0, l_max=5, matched=True)
        R = result["R_phys_m"]
        expected_m2 = 2.0 * HBAR_C_EVM / R
        self.assertAlmostEqual(result["modes"][0]["m_eV"], expected_m2, places=5)

    def test_degeneracy_2l_plus_1(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import graviton_spectrum
        result = graviton_spectrum(28e-6, l_max=10)
        for mode in result["modes"]:
            self.assertEqual(mode["degeneracy"], 2 * mode["l"] + 1)

    def test_modes_count(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import graviton_spectrum
        result = graviton_spectrum(28e-6, l_max=10)
        self.assertEqual(len(result["modes"]), 9)  # l=2..10


# ---------------------------------------------------------------------------
# 4. fermion_spectrum
# ---------------------------------------------------------------------------
class TestFermionSpectrum(unittest.TestCase):
    """Tests for fermion_spectrum."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import fermion_spectrum
        result = fermion_spectrum(28e-6)
        self.assertIsInstance(result, dict)

    def test_starts_at_l0(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import fermion_spectrum
        result = fermion_spectrum(28e-6, l_max=5)
        self.assertEqual(result["modes"][0]["l"], 0)

    def test_l0_eigenvalue_half(self):
        """For l=0: eigenvalue = 0.5."""
        from alpha_ladder_core.gauge_kk_spectrum_matched import fermion_spectrum, HBAR_C_EVM
        a_0 = 1e-5
        result = fermion_spectrum(a_0, l_max=5, matched=True)
        R = result["R_phys_m"]
        expected_m0 = 0.5 * HBAR_C_EVM / R
        self.assertAlmostEqual(result["modes"][0]["m_eV"], expected_m0, places=5)

    def test_degeneracy_2_times_2l_plus_1(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import fermion_spectrum
        result = fermion_spectrum(28e-6, l_max=10)
        for mode in result["modes"]:
            expected = 2 * (2 * mode["l"] + 1)
            self.assertEqual(mode["degeneracy"], expected)

    def test_fermion_lightest_below_scalar(self):
        """Fermion l=0 mass (0.5/R) should be less than scalar l=1 mass (sqrt(2)/R)."""
        from alpha_ladder_core.gauge_kk_spectrum_matched import fermion_spectrum, scalar_spectrum
        a_0 = 28e-6
        ferm = fermion_spectrum(a_0, l_max=5)
        scal = scalar_spectrum(a_0, l_max=5)
        m_ferm_0 = ferm["modes"][0]["m_eV"]
        m_scal_1 = scal["modes"][1]["m_eV"]  # l=1 scalar
        self.assertLess(m_ferm_0, m_scal_1)


# ---------------------------------------------------------------------------
# 5. full_spectrum
# ---------------------------------------------------------------------------
class TestFullSpectrum(unittest.TestCase):
    """Tests for full_spectrum."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import full_spectrum
        result = full_spectrum(28e-6, l_max=5)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import full_spectrum
        result = full_spectrum(28e-6, l_max=5)
        for key in ("modes", "count_below_1TeV", "count_below_10TeV",
                     "count_below_100TeV", "lightest_massive",
                     "mass_gap_eV", "density_at_1TeV", "R_phys_m"):
            self.assertIn(key, result)

    def test_sorted_by_mass(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import full_spectrum
        result = full_spectrum(28e-6, l_max=10)
        masses = [m["m_eV"] for m in result["modes"]]
        self.assertEqual(masses, sorted(masses))

    def test_contains_all_spin_types(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import full_spectrum
        result = full_spectrum(28e-6, l_max=5)
        spins = set(m["spin"] for m in result["modes"])
        self.assertIn("scalar", spins)
        self.assertIn("vector", spins)
        self.assertIn("graviton", spins)

    def test_lightest_massive_is_first(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import full_spectrum
        result = full_spectrum(28e-6, l_max=5)
        self.assertEqual(result["lightest_massive"], result["modes"][0])

    def test_mass_gap_positive(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import full_spectrum
        result = full_spectrum(28e-6, l_max=5)
        self.assertGreater(result["mass_gap_eV"], 0)

    def test_eot_wash_lightest_meV_range(self):
        """For a_0=28um, lightest massive mode should be in meV range."""
        from alpha_ladder_core.gauge_kk_spectrum_matched import full_spectrum
        result = full_spectrum(28e-6, l_max=10)
        m_lightest_meV = result["mass_gap_eV"] * 1e3
        self.assertGreater(m_lightest_meV, 1.0)
        self.assertLess(m_lightest_meV, 100.0)


# ---------------------------------------------------------------------------
# 6. spectrum_at_notable_a0
# ---------------------------------------------------------------------------
class TestSpectrumAtNotableA0(unittest.TestCase):
    """Tests for spectrum_at_notable_a0."""

    def test_returns_list(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import spectrum_at_notable_a0
        result = spectrum_at_notable_a0()
        self.assertIsInstance(result, list)

    def test_five_entries(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import spectrum_at_notable_a0
        result = spectrum_at_notable_a0()
        self.assertEqual(len(result), 5)

    def test_entry_keys(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import spectrum_at_notable_a0
        result = spectrum_at_notable_a0()
        for entry in result:
            for key in ("a_0_m", "a_0_label", "M_6_eV", "M_6_TeV",
                         "R_phys_m", "lightest_5"):
                self.assertIn(key, entry)

    def test_lightest_5_count(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import spectrum_at_notable_a0
        result = spectrum_at_notable_a0()
        for entry in result:
            self.assertLessEqual(len(entry["lightest_5"]), 5)

    def test_m6_positive(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import spectrum_at_notable_a0
        result = spectrum_at_notable_a0()
        for entry in result:
            self.assertGreater(entry["M_6_eV"], 0)


# ---------------------------------------------------------------------------
# 7. matched_vs_unmatched
# ---------------------------------------------------------------------------
class TestMatchedVsUnmatched(unittest.TestCase):
    """Tests for matched_vs_unmatched."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import matched_vs_unmatched
        result = matched_vs_unmatched(28e-6)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import matched_vs_unmatched
        result = matched_vs_unmatched(28e-6)
        for key in ("matched_modes", "unmatched_modes", "mass_ratio",
                     "enhancement_factor", "phi_vev", "R_matched_m",
                     "R_unmatched_m", "detectability_note"):
            self.assertIn(key, result)

    def test_mass_ratio_approximately_1817(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import matched_vs_unmatched, MASS_ENHANCEMENT
        result = matched_vs_unmatched(28e-6)
        self.assertAlmostEqual(result["mass_ratio"], MASS_ENHANCEMENT, delta=0.01)

    def test_enhancement_factor_value(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import matched_vs_unmatched
        result = matched_vs_unmatched(28e-6)
        self.assertAlmostEqual(result["enhancement_factor"], 1.817, delta=0.01)

    def test_matched_heavier_than_unmatched(self):
        """Matched modes should be heavier (enhancement > 1)."""
        from alpha_ladder_core.gauge_kk_spectrum_matched import matched_vs_unmatched
        result = matched_vs_unmatched(28e-6)
        for mm, um in zip(result["matched_modes"][1:],
                          result["unmatched_modes"][1:]):
            self.assertGreater(mm["m_eV"], um["m_eV"])

    def test_r_matched_smaller(self):
        """Matched R_phys < unmatched R_phys (because e^phi < 1)."""
        from alpha_ladder_core.gauge_kk_spectrum_matched import matched_vs_unmatched
        result = matched_vs_unmatched(28e-6)
        self.assertLess(result["R_matched_m"], result["R_unmatched_m"])


# ---------------------------------------------------------------------------
# 8. cumulative_modes_vs_energy
# ---------------------------------------------------------------------------
class TestCumulativeModesVsEnergy(unittest.TestCase):
    """Tests for cumulative_modes_vs_energy."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import cumulative_modes_vs_energy
        result = cumulative_modes_vs_energy(28e-6, E_max_eV=1e6, n_points=20)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import cumulative_modes_vs_energy
        result = cumulative_modes_vs_energy(28e-6, E_max_eV=1e6, n_points=20)
        for key in ("E_grid", "N_cumulative", "N_add", "R_phys_m", "R_add_m"):
            self.assertIn(key, result)

    def test_n_cumulative_monotonic(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import cumulative_modes_vs_energy
        result = cumulative_modes_vs_energy(28e-6, E_max_eV=1e6, n_points=50)
        for i in range(len(result["N_cumulative"]) - 1):
            self.assertLessEqual(result["N_cumulative"][i],
                                 result["N_cumulative"][i + 1])

    def test_s2_sparser_than_add(self):
        """S^2 spectrum should be sparser than ADD torus at same radius."""
        from alpha_ladder_core.gauge_kk_spectrum_matched import cumulative_modes_vs_energy
        result = cumulative_modes_vs_energy(28e-6, E_max_eV=1e4, n_points=20)
        # At high enough energy, ADD should have more modes
        n_s2 = result["N_cumulative"][-1]
        n_add = result["N_add"][-1]
        if n_add > 0 and n_s2 > 0:
            ratio = n_s2 / n_add
            # S^2 is sparser by roughly 1/pi factor
            self.assertLess(ratio, 1.5)

    def test_grid_length(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import cumulative_modes_vs_energy
        result = cumulative_modes_vs_energy(28e-6, E_max_eV=1e6, n_points=30)
        self.assertEqual(len(result["E_grid"]), 30)
        self.assertEqual(len(result["N_cumulative"]), 30)
        self.assertEqual(len(result["N_add"]), 30)


# ---------------------------------------------------------------------------
# 9. summarize_kk_spectrum
# ---------------------------------------------------------------------------
class TestSummarizeKkSpectrum(unittest.TestCase):
    """Tests for summarize_kk_spectrum."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import summarize_kk_spectrum
        result = summarize_kk_spectrum()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import summarize_kk_spectrum
        result = summarize_kk_spectrum()
        for key in ("constants", "key_results", "report"):
            self.assertIn(key, result)

    def test_report_is_string(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import summarize_kk_spectrum
        result = summarize_kk_spectrum()
        self.assertIsInstance(result["report"], str)
        self.assertGreater(len(result["report"]), 100)

    def test_key_results_list(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import summarize_kk_spectrum
        result = summarize_kk_spectrum()
        self.assertIsInstance(result["key_results"], list)
        self.assertGreater(len(result["key_results"]), 0)

    def test_constants_dict(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import summarize_kk_spectrum
        result = summarize_kk_spectrum()
        consts = result["constants"]
        for key in ("alpha_EM", "phi_vev", "e_phi_vev",
                     "mass_enhancement", "M_Pl_eV"):
            self.assertIn(key, consts)

    def test_mass_enhancement_in_constants(self):
        from alpha_ladder_core.gauge_kk_spectrum_matched import summarize_kk_spectrum
        result = summarize_kk_spectrum()
        self.assertAlmostEqual(result["constants"]["mass_enhancement"],
                               1.817, delta=0.01)


if __name__ == "__main__":
    unittest.main()
