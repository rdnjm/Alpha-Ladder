"""Tests for alpha_running.py -- Running of alpha_EM with KK thresholds."""
import unittest
import math



# ---------------------------------------------------------------------------
# 1. sm_alpha_running
# ---------------------------------------------------------------------------
class TestSmAlphaRunning(unittest.TestCase):
    """Tests for sm_alpha_running."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_alpha_running import sm_alpha_running
        result = sm_alpha_running(1e9)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_alpha_running import sm_alpha_running
        result = sm_alpha_running(1e9)
        for key in ("alpha", "inv_alpha", "active", "mu_eV", "description"):
            self.assertIn(key, result)

    def test_alpha_at_me_is_thomson(self):
        """alpha(m_e) should be approximately 1/137.036."""
        from alpha_ladder_core.gauge_alpha_running import sm_alpha_running, M_E
        result = sm_alpha_running(M_E)
        self.assertAlmostEqual(result["inv_alpha"], 137.036, delta=0.01)

    def test_alpha_at_mz_approximately_128(self):
        """alpha(m_Z) should be around 1/128 (one-loop QED approximation)."""
        from alpha_ladder_core.gauge_alpha_running import sm_alpha_running, M_Z
        result = sm_alpha_running(M_Z)
        self.assertGreater(result["inv_alpha"], 125.0)
        self.assertLess(result["inv_alpha"], 135.0)

    def test_alpha_below_me_is_frozen(self):
        """Below electron mass, alpha is frozen at Thomson limit."""
        from alpha_ladder_core.gauge_alpha_running import sm_alpha_running, ALPHA_EM_0
        result = sm_alpha_running(1.0)
        self.assertAlmostEqual(result["alpha"], ALPHA_EM_0, places=12)

    def test_active_fermions_empty_below_me(self):
        from alpha_ladder_core.gauge_alpha_running import sm_alpha_running
        result = sm_alpha_running(1.0)
        self.assertEqual(result["active"], [])

    def test_active_fermions_at_1gev(self):
        """At 1 GeV, electron, up, down, muon, strange should be active."""
        from alpha_ladder_core.gauge_alpha_running import sm_alpha_running
        result = sm_alpha_running(1e9)
        self.assertIn("electron", result["active"])
        self.assertIn("up", result["active"])
        self.assertIn("down", result["active"])

    def test_all_fermions_active_above_top(self):
        """Above top mass, all 9 SM fermions should be active."""
        from alpha_ladder_core.gauge_alpha_running import sm_alpha_running
        result = sm_alpha_running(1e12)
        self.assertEqual(len(result["active"]), 9)

    def test_mu_ev_stored(self):
        from alpha_ladder_core.gauge_alpha_running import sm_alpha_running
        result = sm_alpha_running(42.0)
        self.assertEqual(result["mu_eV"], 42.0)

    def test_alpha_increases_with_energy(self):
        """alpha should increase (1/alpha decreases) at higher energy."""
        from alpha_ladder_core.gauge_alpha_running import sm_alpha_running
        low = sm_alpha_running(1e6)
        high = sm_alpha_running(1e11)
        self.assertGreater(low["inv_alpha"], high["inv_alpha"])

    def test_zero_scale(self):
        """Scale <= 0 should return Thomson limit."""
        from alpha_ladder_core.gauge_alpha_running import sm_alpha_running, ALPHA_EM_0
        result = sm_alpha_running(0.0)
        self.assertAlmostEqual(result["alpha"], ALPHA_EM_0, places=12)

    def test_negative_scale(self):
        from alpha_ladder_core.gauge_alpha_running import sm_alpha_running, ALPHA_EM_0
        result = sm_alpha_running(-100.0)
        self.assertAlmostEqual(result["alpha"], ALPHA_EM_0, places=12)

    def test_description_is_string(self):
        from alpha_ladder_core.gauge_alpha_running import sm_alpha_running
        result = sm_alpha_running(1e9)
        self.assertIsInstance(result["description"], str)


# ---------------------------------------------------------------------------
# 2. kk_alpha_running
# ---------------------------------------------------------------------------
class TestKkAlphaRunning(unittest.TestCase):
    """Tests for kk_alpha_running."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_alpha_running import kk_alpha_running
        result = kk_alpha_running(1e9, 28e-6)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_alpha_running import kk_alpha_running
        result = kk_alpha_running(1e9, 28e-6)
        for key in ("alpha", "inv_alpha", "n_kk_active", "m_kk_lightest_eV",
                     "active_sm", "mu_eV", "description"):
            self.assertIn(key, result)

    def test_minimal_framework_no_kk_correction(self):
        """n_charged=0 should yield exactly the SM result."""
        from alpha_ladder_core.gauge_alpha_running import kk_alpha_running, sm_alpha_running
        mu = 1e9
        sm = sm_alpha_running(mu)
        kk = kk_alpha_running(mu, 28e-6, n_charged=0)
        self.assertAlmostEqual(kk["alpha"], sm["alpha"], places=12)
        self.assertEqual(kk["n_kk_active"], 0)

    def test_lightest_kk_mass_meV_range_28um(self):
        """For a_0=28 um, lightest KK mass should be in the meV range."""
        from alpha_ladder_core.gauge_alpha_running import kk_alpha_running
        result = kk_alpha_running(1e9, 28e-6, n_charged=1)
        m_kk = result["m_kk_lightest_eV"]
        self.assertGreater(m_kk, 1e-4)
        self.assertLess(m_kk, 1.0)

    def test_kk_shift_huge_for_28um(self):
        """With charged scalars at a_0=28um, the KK shift should be enormous."""
        from alpha_ladder_core.gauge_alpha_running import kk_alpha_running, sm_alpha_running
        mu = 1e9
        sm = sm_alpha_running(mu)
        kk = kk_alpha_running(mu, 28e-6, n_charged=1)
        shift_ppm = abs(kk["inv_alpha"] - sm["inv_alpha"]) / sm["inv_alpha"] * 1e6
        self.assertGreater(shift_ppm, 1e10)

    def test_many_kk_modes_active_at_high_energy(self):
        """At high energy with large a_0, many KK modes should be active."""
        from alpha_ladder_core.gauge_alpha_running import kk_alpha_running
        result = kk_alpha_running(1e11, 28e-6, n_charged=1)
        self.assertGreater(result["n_kk_active"], 1000)

    def test_below_kk_threshold_no_modes(self):
        """Below the KK scale, no KK modes should be active."""
        from alpha_ladder_core.gauge_alpha_running import kk_alpha_running
        result = kk_alpha_running(1e-6, 1e-35, n_charged=1)
        self.assertEqual(result["n_kk_active"], 0)

    def test_unmatched_radius(self):
        """With matched=False, should use bare radius."""
        from alpha_ladder_core.gauge_alpha_running import kk_alpha_running
        matched = kk_alpha_running(1e9, 28e-6, n_charged=1, matched=True)
        unmatched = kk_alpha_running(1e9, 28e-6, n_charged=1, matched=False)
        self.assertNotAlmostEqual(matched["m_kk_lightest_eV"],
                                  unmatched["m_kk_lightest_eV"], places=2)

    def test_l_max_override(self):
        """Explicit l_max should cap the number of modes."""
        from alpha_ladder_core.gauge_alpha_running import kk_alpha_running
        result = kk_alpha_running(1e15, 28e-6, n_charged=1, l_max=5)
        self.assertLessEqual(result["n_kk_active"], 5 * (5 + 2))

    def test_higher_charge_increases_shift(self):
        """Higher Q_charge should produce a larger shift."""
        from alpha_ladder_core.gauge_alpha_running import kk_alpha_running
        q1 = kk_alpha_running(1e9, 28e-6, n_charged=1, Q_charge=1.0)
        q2 = kk_alpha_running(1e9, 28e-6, n_charged=1, Q_charge=2.0)
        self.assertGreater(abs(q2["inv_alpha"] - q1["inv_alpha"]), 0)


# ---------------------------------------------------------------------------
# 3. compare_sm_vs_kk
# ---------------------------------------------------------------------------
class TestCompareSmVsKk(unittest.TestCase):
    """Tests for compare_sm_vs_kk."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_alpha_running import compare_sm_vs_kk
        result = compare_sm_vs_kk(28e-6)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_alpha_running import compare_sm_vs_kk
        result = compare_sm_vs_kk(28e-6)
        for key in ("table", "a_0_m", "m_kk_eV", "M_6_eV", "description"):
            self.assertIn(key, result)

    def test_table_has_entries(self):
        from alpha_ladder_core.gauge_alpha_running import compare_sm_vs_kk
        result = compare_sm_vs_kk(28e-6)
        self.assertGreater(len(result["table"]), 0)

    def test_table_row_keys(self):
        from alpha_ladder_core.gauge_alpha_running import compare_sm_vs_kk
        result = compare_sm_vs_kk(28e-6)
        row = result["table"][0]
        for key in ("label", "mu_eV", "alpha_SM", "alpha_KK",
                     "delta_alpha_ppm", "n_kk_active", "inv_alpha_SM",
                     "inv_alpha_KK"):
            self.assertIn(key, row)

    def test_table_sorted_by_energy(self):
        from alpha_ladder_core.gauge_alpha_running import compare_sm_vs_kk
        result = compare_sm_vs_kk(28e-6)
        energies = [row["mu_eV"] for row in result["table"]]
        self.assertEqual(energies, sorted(energies))

    def test_delta_ppm_zero_for_minimal(self):
        """With n_charged=0, delta_alpha_ppm should be zero at all scales."""
        from alpha_ladder_core.gauge_alpha_running import compare_sm_vs_kk
        result = compare_sm_vs_kk(28e-6, n_charged=0)
        for row in result["table"]:
            self.assertAlmostEqual(row["delta_alpha_ppm"], 0.0, places=6)

    def test_m6_positive(self):
        from alpha_ladder_core.gauge_alpha_running import compare_sm_vs_kk
        result = compare_sm_vs_kk(28e-6)
        self.assertGreater(result["M_6_eV"], 0)


# ---------------------------------------------------------------------------
# 4. does_it_match_observation
# ---------------------------------------------------------------------------
class TestDoesItMatchObservation(unittest.TestCase):
    """Tests for does_it_match_observation."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_alpha_running import does_it_match_observation
        result = does_it_match_observation(28e-6)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_alpha_running import does_it_match_observation
        result = does_it_match_observation(28e-6)
        for key in ("alpha_predicted_at_me", "alpha_predicted_at_mZ",
                     "alpha_observed_at_me", "alpha_observed_at_mZ",
                     "shift_at_me_ppm", "shift_at_mZ_ppm", "consistent",
                     "m_kk_eV", "description"):
            self.assertIn(key, result)

    def test_minimal_is_consistent(self):
        """n_charged=0 must always be consistent."""
        from alpha_ladder_core.gauge_alpha_running import does_it_match_observation
        result = does_it_match_observation(28e-6, n_charged=0)
        self.assertTrue(result["consistent"])

    def test_n1_28um_inconsistent(self):
        """n_charged=1 at a_0=28um should be inconsistent (huge shift)."""
        from alpha_ladder_core.gauge_alpha_running import does_it_match_observation
        result = does_it_match_observation(28e-6, n_charged=1)
        self.assertFalse(result["consistent"])

    def test_shift_at_me_huge_for_28um(self):
        """The shift at m_e should be enormous for charged KK at 28um."""
        from alpha_ladder_core.gauge_alpha_running import does_it_match_observation
        result = does_it_match_observation(28e-6, n_charged=1)
        self.assertGreater(abs(result["shift_at_me_ppm"]), 1e10)

    def test_observed_alpha_values_correct(self):
        from alpha_ladder_core.gauge_alpha_running import does_it_match_observation, ALPHA_EM_0, ALPHA_EM_MZ
        result = does_it_match_observation(28e-6)
        self.assertAlmostEqual(result["alpha_observed_at_me"], ALPHA_EM_0, places=12)
        self.assertAlmostEqual(result["alpha_observed_at_mZ"], ALPHA_EM_MZ, places=12)

    def test_description_mentions_inconsistent(self):
        from alpha_ladder_core.gauge_alpha_running import does_it_match_observation
        result = does_it_match_observation(28e-6, n_charged=1)
        self.assertIn("INCONSISTENT", result["description"])

    def test_minimal_description_mentions_consistent(self):
        from alpha_ladder_core.gauge_alpha_running import does_it_match_observation
        result = does_it_match_observation(28e-6, n_charged=0)
        self.assertIn("consistent", result["description"].lower())


# ---------------------------------------------------------------------------
# 5. running_profile
# ---------------------------------------------------------------------------
class TestRunningProfile(unittest.TestCase):
    """Tests for running_profile."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_alpha_running import running_profile
        result = running_profile(28e-6, n_charged=0, n_points=10)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_alpha_running import running_profile
        result = running_profile(28e-6, n_charged=0, n_points=10)
        for key in ("mu_eV", "inv_alpha_SM", "inv_alpha_KK", "m_kk_eV",
                     "n_points", "description"):
            self.assertIn(key, result)

    def test_grid_length(self):
        from alpha_ladder_core.gauge_alpha_running import running_profile
        result = running_profile(28e-6, n_charged=0, n_points=50)
        self.assertEqual(len(result["mu_eV"]), 50)
        self.assertEqual(len(result["inv_alpha_SM"]), 50)
        self.assertEqual(len(result["inv_alpha_KK"]), 50)

    def test_mu_monotonically_increasing(self):
        from alpha_ladder_core.gauge_alpha_running import running_profile
        result = running_profile(28e-6, n_charged=0, n_points=20)
        for i in range(len(result["mu_eV"]) - 1):
            self.assertLess(result["mu_eV"][i], result["mu_eV"][i + 1])

    def test_sm_equals_kk_for_minimal(self):
        """With n_charged=0, SM and KK profiles should be identical."""
        from alpha_ladder_core.gauge_alpha_running import running_profile
        result = running_profile(28e-6, n_charged=0, n_points=20)
        for sm, kk in zip(result["inv_alpha_SM"], result["inv_alpha_KK"]):
            self.assertAlmostEqual(sm, kk, places=10)

    def test_m_kk_positive(self):
        from alpha_ladder_core.gauge_alpha_running import running_profile
        result = running_profile(28e-6, n_charged=1, n_points=10)
        self.assertGreater(result["m_kk_eV"], 0)


# ---------------------------------------------------------------------------
# 6. summarize_alpha_running
# ---------------------------------------------------------------------------
class TestSummarizeAlphaRunning(unittest.TestCase):
    """Tests for summarize_alpha_running."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_alpha_running import summarize_alpha_running
        result = summarize_alpha_running()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_alpha_running import summarize_alpha_running
        result = summarize_alpha_running()
        for key in ("sm_check", "minimal_check", "charged_check",
                     "comparison", "key_messages", "report"):
            self.assertIn(key, result)

    def test_sm_check_inv_alpha_me(self):
        from alpha_ladder_core.gauge_alpha_running import summarize_alpha_running
        result = summarize_alpha_running()
        self.assertAlmostEqual(result["sm_check"]["inv_alpha_at_me"],
                               137.036, delta=0.01)

    def test_minimal_check_is_consistent(self):
        from alpha_ladder_core.gauge_alpha_running import summarize_alpha_running
        result = summarize_alpha_running()
        self.assertTrue(result["minimal_check"]["consistent"])

    def test_charged_check_is_inconsistent(self):
        from alpha_ladder_core.gauge_alpha_running import summarize_alpha_running
        result = summarize_alpha_running()
        self.assertFalse(result["charged_check"]["consistent"])

    def test_key_messages_list(self):
        from alpha_ladder_core.gauge_alpha_running import summarize_alpha_running
        result = summarize_alpha_running()
        self.assertIsInstance(result["key_messages"], list)
        self.assertGreater(len(result["key_messages"]), 0)

    def test_report_is_string(self):
        from alpha_ladder_core.gauge_alpha_running import summarize_alpha_running
        result = summarize_alpha_running()
        self.assertIsInstance(result["report"], str)
        self.assertGreater(len(result["report"]), 100)

    def test_comparison_has_table(self):
        from alpha_ladder_core.gauge_alpha_running import summarize_alpha_running
        result = summarize_alpha_running()
        self.assertIn("table", result["comparison"])
        self.assertGreater(len(result["comparison"]["table"]), 0)


# ---------------------------------------------------------------------------
# Module-level constants and helpers
# ---------------------------------------------------------------------------
class TestModuleConstants(unittest.TestCase):
    """Tests for module-level constants."""

    def test_alpha_em_0_value(self):
        from alpha_ladder_core.gauge_alpha_running import ALPHA_EM_0
        self.assertAlmostEqual(1.0 / ALPHA_EM_0, 137.036, delta=0.001)

    def test_alpha_em_mz_value(self):
        from alpha_ladder_core.gauge_alpha_running import ALPHA_EM_MZ
        self.assertAlmostEqual(1.0 / ALPHA_EM_MZ, 127.944, delta=0.001)

    def test_sm_fermions_sorted_by_mass(self):
        from alpha_ladder_core.gauge_alpha_running import SM_FERMIONS
        masses = [f[1] for f in SM_FERMIONS]
        self.assertEqual(masses, sorted(masses))

    def test_sm_fermions_count(self):
        from alpha_ladder_core.gauge_alpha_running import SM_FERMIONS
        self.assertEqual(len(SM_FERMIONS), 9)

    def test_golden_ratio(self):
        from alpha_ladder_core.gauge_alpha_running import PHI_GOLDEN
        self.assertAlmostEqual(PHI_GOLDEN, (1 + math.sqrt(5)) / 2, places=10)

    def test_b_coefficient(self):
        from alpha_ladder_core.gauge_alpha_running import _b_coefficient
        # electron: Q=1, Nc=1 -> 1
        self.assertAlmostEqual(_b_coefficient(1.0, 1), 1.0)
        # up quark: Q=2/3, Nc=3 -> 3*(4/9) = 4/3
        self.assertAlmostEqual(_b_coefficient(2.0 / 3, 3), 4.0 / 3, places=10)


if __name__ == "__main__":
    unittest.main()
