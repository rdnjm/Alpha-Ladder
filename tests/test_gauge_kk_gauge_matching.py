"""Tests for kk_gauge_matching.py -- KK gauge matching and Coleman-Weinberg potential."""
import unittest
import math



# ---------------------------------------------------------------------------
# 1. tree_level_matching
# ---------------------------------------------------------------------------
class TestTreeLevelMatching(unittest.TestCase):
    """Tests for tree_level_matching."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import tree_level_matching
        result = tree_level_matching()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import tree_level_matching
        result = tree_level_matching()
        for key in ("phi_vev", "g_KK", "alpha_KK", "alpha_EM",
                     "e_to_phi", "e_to_2phi", "e_to_4phi", "verification",
                     "description"):
            self.assertIn(key, result)

    def test_phi_vev_value(self):
        """phi_vev = (1/4)*ln(4*pi*alpha_EM) should be approximately -0.597."""
        from alpha_ladder_core.gauge_kk_gauge_matching import tree_level_matching
        result = tree_level_matching()
        self.assertAlmostEqual(result["phi_vev"], -0.597, delta=0.01)

    def test_phi_vev_negative(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import tree_level_matching
        result = tree_level_matching()
        self.assertLess(result["phi_vev"], 0)

    def test_alpha_kk_equals_alpha_em(self):
        """alpha_KK should equal alpha_EM to high precision."""
        from alpha_ladder_core.gauge_kk_gauge_matching import tree_level_matching, ALPHA_EM
        result = tree_level_matching()
        self.assertAlmostEqual(result["alpha_KK"], ALPHA_EM, places=12)

    def test_verification_true(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import tree_level_matching
        result = tree_level_matching()
        self.assertTrue(result["verification"])

    def test_custom_alpha(self):
        """Providing a custom alpha should give a different phi_vev."""
        from alpha_ladder_core.gauge_kk_gauge_matching import tree_level_matching
        result = tree_level_matching(alpha_em=0.01)
        expected = 0.25 * math.log(4.0 * math.pi * 0.01)
        self.assertAlmostEqual(result["phi_vev"], expected, places=10)

    def test_custom_alpha_verification(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import tree_level_matching
        result = tree_level_matching(alpha_em=0.01)
        self.assertTrue(result["verification"])

    def test_e_to_phi_consistent(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import tree_level_matching
        result = tree_level_matching()
        self.assertAlmostEqual(result["e_to_phi"],
                               math.exp(result["phi_vev"]), places=10)

    def test_e_to_4phi_consistent(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import tree_level_matching
        result = tree_level_matching()
        self.assertAlmostEqual(result["e_to_4phi"],
                               math.exp(4 * result["phi_vev"]), places=10)


# ---------------------------------------------------------------------------
# 2. kk_mass_spectrum
# ---------------------------------------------------------------------------
class TestKkMassSpectrum(unittest.TestCase):
    """Tests for kk_mass_spectrum."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import kk_mass_spectrum
        result = kk_mass_spectrum(-0.597, a_0=1.0, l_max=10)
        self.assertIsInstance(result, dict)

    def test_modes_count(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import kk_mass_spectrum
        result = kk_mass_spectrum(-0.597, a_0=1.0, l_max=10)
        self.assertEqual(len(result["modes"]), 10)

    def test_modes_start_at_l1(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import kk_mass_spectrum
        result = kk_mass_spectrum(-0.597, a_0=1.0, l_max=5)
        self.assertEqual(result["modes"][0]["l"], 1)

    def test_mass_formula(self):
        """m_l = sqrt(l(l+1)) / R_phys for each mode."""
        from alpha_ladder_core.gauge_kk_gauge_matching import kk_mass_spectrum
        phi = -0.5
        a_0 = 2.0
        result = kk_mass_spectrum(phi, a_0=a_0, l_max=5)
        R_phys = a_0 * math.exp(phi)
        for mode in result["modes"]:
            l = mode["l"]
            expected_m = math.sqrt(l * (l + 1)) / R_phys
            self.assertAlmostEqual(mode["m_l"], expected_m, places=10)

    def test_degeneracy_2l_plus_1(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import kk_mass_spectrum
        result = kk_mass_spectrum(-0.597, a_0=1.0, l_max=10)
        for mode in result["modes"]:
            self.assertEqual(mode["degeneracy"], 2 * mode["l"] + 1)

    def test_m_l_over_m1_ratio(self):
        """The ratio m_l/m_1 should be sqrt(l(l+1)/2)."""
        from alpha_ladder_core.gauge_kk_gauge_matching import kk_mass_spectrum
        result = kk_mass_spectrum(-0.5, a_0=1.0, l_max=5)
        for mode in result["modes"]:
            l = mode["l"]
            expected = math.sqrt(l * (l + 1) / 2.0)
            self.assertAlmostEqual(mode["m_l_over_m1"], expected, places=8)

    def test_r_phys_stored(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import kk_mass_spectrum
        phi = -0.5
        a_0 = 3.0
        result = kk_mass_spectrum(phi, a_0=a_0)
        self.assertAlmostEqual(result["R_phys"], a_0 * math.exp(phi), places=10)

    def test_masses_increase_with_l(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import kk_mass_spectrum
        result = kk_mass_spectrum(-0.5, a_0=1.0, l_max=10)
        masses = [m["m_l"] for m in result["modes"]]
        for i in range(len(masses) - 1):
            self.assertLess(masses[i], masses[i + 1])


# ---------------------------------------------------------------------------
# 3. coleman_weinberg_potential
# ---------------------------------------------------------------------------
class TestColemanWeinbergPotential(unittest.TestCase):
    """Tests for coleman_weinberg_potential."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import coleman_weinberg_potential
        result = coleman_weinberg_potential(-0.5, a_0=1.0)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import coleman_weinberg_potential
        result = coleman_weinberg_potential(-0.5)
        for key in ("V_CW", "dV_dphi", "phi"):
            self.assertIn(key, result)

    def test_phi_stored(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import coleman_weinberg_potential
        result = coleman_weinberg_potential(-0.3)
        self.assertAlmostEqual(result["phi"], -0.3, places=10)

    def test_v_finite(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import coleman_weinberg_potential
        result = coleman_weinberg_potential(-1.0, a_0=1.0, l_max=20)
        self.assertTrue(math.isfinite(result["V_CW"]))
        self.assertTrue(math.isfinite(result["dV_dphi"]))

    def test_v_changes_with_phi(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import coleman_weinberg_potential
        v1 = coleman_weinberg_potential(-2.0)["V_CW"]
        v2 = coleman_weinberg_potential(0.0)["V_CW"]
        self.assertNotAlmostEqual(v1, v2, places=5)


# ---------------------------------------------------------------------------
# 3b. coleman_weinberg_scan
# ---------------------------------------------------------------------------
class TestColemanWeinberyScan(unittest.TestCase):
    """Tests for coleman_weinberg_scan."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import coleman_weinberg_scan
        result = coleman_weinberg_scan(n_points=20)
        self.assertIsInstance(result, dict)

    def test_grid_length(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import coleman_weinberg_scan
        result = coleman_weinberg_scan(n_points=30)
        self.assertEqual(len(result["phi_grid"]), 30)
        self.assertEqual(len(result["V_CW"]), 30)
        self.assertEqual(len(result["dV_dphi"]), 30)

    def test_has_minimum_key(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import coleman_weinberg_scan
        result = coleman_weinberg_scan(n_points=50)
        self.assertIn("has_minimum", result)
        self.assertIsInstance(result["has_minimum"], bool)


# ---------------------------------------------------------------------------
# 4. find_cw_minimum
# ---------------------------------------------------------------------------
class TestFindCwMinimum(unittest.TestCase):
    """Tests for find_cw_minimum."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import find_cw_minimum
        result = find_cw_minimum(a_0=1.0, n_charged=1, l_max=20)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import find_cw_minimum
        result = find_cw_minimum(a_0=1.0, n_charged=1, l_max=20)
        for key in ("phi_min", "phi_vev", "ratio", "m_phi_squared",
                     "m_phi", "coincidence_with_alpha_matching",
                     "description", "honest_assessment"):
            self.assertIn(key, result)

    def test_phi_vev_matches_tree_level(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import find_cw_minimum, tree_level_matching
        tree = tree_level_matching()
        result = find_cw_minimum(a_0=1.0, n_charged=1, l_max=20)
        self.assertAlmostEqual(result["phi_vev"], tree["phi_vev"], places=10)

    def test_honest_assessment_present(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import find_cw_minimum
        result = find_cw_minimum(a_0=1.0, n_charged=1, l_max=20)
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 20)


# ---------------------------------------------------------------------------
# 5. scan_n_charged
# ---------------------------------------------------------------------------
class TestScanNCharged(unittest.TestCase):
    """Tests for scan_n_charged."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import scan_n_charged
        result = scan_n_charged(n_values=[1, 2], l_max=20)
        self.assertIsInstance(result, dict)

    def test_results_length(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import scan_n_charged
        result = scan_n_charged(n_values=[1, 5, 10], l_max=20)
        self.assertEqual(len(result["results"]), 3)

    def test_phi_vev_present(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import scan_n_charged
        result = scan_n_charged(n_values=[1], l_max=20)
        self.assertIn("phi_vev", result)

    def test_result_row_keys(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import scan_n_charged
        result = scan_n_charged(n_values=[1], l_max=20)
        row = result["results"][0]
        for key in ("n_charged", "phi_min", "ratio_to_vev",
                     "m_phi_squared", "deviation_percent"):
            self.assertIn(key, row)


# ---------------------------------------------------------------------------
# 6. alpha_running (gauge module version)
# ---------------------------------------------------------------------------
class TestAlphaRunningGauge(unittest.TestCase):
    """Tests for alpha_running in kk_gauge_matching module."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import alpha_running
        result = alpha_running(-0.597, a_0=1.0, n_steps=10)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import alpha_running
        result = alpha_running(-0.597, a_0=1.0, n_steps=10)
        for key in ("running", "alpha_at_lab", "alpha_at_compactification",
                     "description"):
            self.assertIn(key, result)

    def test_running_length(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import alpha_running
        result = alpha_running(-0.597, a_0=1.0, n_steps=20)
        self.assertEqual(len(result["running"]), 20)

    def test_running_point_keys(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import alpha_running
        result = alpha_running(-0.597, a_0=1.0, n_steps=10)
        pt = result["running"][0]
        for key in ("mu", "alpha_at_mu", "n_active_modes"):
            self.assertIn(key, pt)

    def test_alpha_at_compactification(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import alpha_running, ALPHA_EM
        result = alpha_running(-0.597, a_0=1.0, n_steps=10)
        self.assertAlmostEqual(result["alpha_at_compactification"],
                               ALPHA_EM, places=10)


# ---------------------------------------------------------------------------
# 7. consistency_with_G
# ---------------------------------------------------------------------------
class TestConsistencyWithG(unittest.TestCase):
    """Tests for consistency_with_G."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import consistency_with_G
        result = consistency_with_G(-0.597)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import consistency_with_G
        result = consistency_with_G(-0.597)
        for key in ("phi_vev", "phi_golden", "relationships_checked",
                     "any_match", "honest_assessment"):
            self.assertIn(key, result)

    def test_phi_golden_value(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import consistency_with_G, PHI_GOLDEN
        result = consistency_with_G(-0.597)
        self.assertAlmostEqual(result["phi_golden"], PHI_GOLDEN, places=10)

    def test_relationships_checked_nonempty(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import consistency_with_G
        result = consistency_with_G(-0.597)
        self.assertGreater(len(result["relationships_checked"]), 0)

    def test_honest_assessment_string(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import consistency_with_G
        result = consistency_with_G(-0.597)
        self.assertIsInstance(result["honest_assessment"], str)


# ---------------------------------------------------------------------------
# 8. summarize_gauge_matching
# ---------------------------------------------------------------------------
class TestSummarizeGaugeMatching(unittest.TestCase):
    """Tests for summarize_gauge_matching."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import summarize_gauge_matching
        result = summarize_gauge_matching()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import summarize_gauge_matching
        result = summarize_gauge_matching()
        for key in ("tree_level", "kk_spectrum", "cw_scan", "cw_minimum",
                     "n_charged_scan", "alpha_running", "G_consistency",
                     "summary"):
            self.assertIn(key, result)

    def test_summary_keys(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import summarize_gauge_matching
        result = summarize_gauge_matching()
        summary = result["summary"]
        for key in ("key_result_1", "key_result_2", "key_result_3",
                     "gap_1_status", "gap_2_status", "honest_caveat"):
            self.assertIn(key, summary)

    def test_tree_level_verification(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import summarize_gauge_matching
        result = summarize_gauge_matching()
        self.assertTrue(result["tree_level"]["verification"])


# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------
class TestGaugeModuleConstants(unittest.TestCase):
    """Tests for module-level constants."""

    def test_alpha_em(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import ALPHA_EM
        self.assertAlmostEqual(1.0 / ALPHA_EM, 137.036, delta=0.001)

    def test_phi_golden(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import PHI_GOLDEN
        self.assertAlmostEqual(PHI_GOLDEN, (1 + math.sqrt(5)) / 2, places=10)

    def test_planck_mass(self):
        from alpha_ladder_core.gauge_kk_gauge_matching import M_PL_EV
        self.assertAlmostEqual(M_PL_EV, 1.22089e28, delta=1e24)


if __name__ == "__main__":
    unittest.main()
