"""Tests for m6_parameter_space.py -- (a_0, M_6) parameter space with dilaton vev constraint."""
import unittest
import math



class TestModuleConstants(unittest.TestCase):
    """Test module-level constants."""

    def test_phi_vev_value(self):
        from alpha_ladder_core.gauge_m6_parameter_space import PHI_VEV
        self.assertAlmostEqual(PHI_VEV, -0.5974, delta=0.01)

    def test_e_phi_vev_value(self):
        from alpha_ladder_core.gauge_m6_parameter_space import E_PHI_VEV
        self.assertAlmostEqual(E_PHI_VEV, 0.5503, delta=0.01)

    def test_phi_vev_from_alpha(self):
        from alpha_ladder_core.gauge_m6_parameter_space import PHI_VEV, ALPHA_EM
        expected = 0.25 * math.log(4.0 * math.pi * ALPHA_EM)
        self.assertAlmostEqual(PHI_VEV, expected, places=10)

    def test_e_phi_vev_consistency(self):
        from alpha_ladder_core.gauge_m6_parameter_space import PHI_VEV, E_PHI_VEV
        self.assertAlmostEqual(E_PHI_VEV, math.exp(PHI_VEV), places=10)

    def test_planck_length_value(self):
        from alpha_ladder_core.gauge_m6_parameter_space import L_PL
        self.assertAlmostEqual(L_PL, 1.61625e-35, delta=1e-38)


class TestComputeParameterPoint(unittest.TestCase):
    """Test compute_parameter_point function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_m6_parameter_space import compute_parameter_point
        result = compute_parameter_point(28e-6)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_m6_parameter_space import compute_parameter_point
        result = compute_parameter_point(28e-6)
        for key in ["a_0_m", "M_6_eV", "M_6_TeV", "R_phys_m", "R_bare_m",
                     "m_KK_eV", "m_KK_l1_eV", "regime"]:
            self.assertIn(key, result)

    def test_r_phys_equals_a0_times_e_phi(self):
        from alpha_ladder_core.gauge_m6_parameter_space import compute_parameter_point, E_PHI_VEV
        a_0 = 28e-6
        result = compute_parameter_point(a_0)
        self.assertAlmostEqual(result["R_phys_m"], a_0 * E_PHI_VEV, places=15)

    def test_r_phys_less_than_a0(self):
        """phi_vev < 0 => R_phys < a_0."""
        from alpha_ladder_core.gauge_m6_parameter_space import compute_parameter_point
        result = compute_parameter_point(28e-6)
        self.assertLess(result["R_phys_m"], result["R_bare_m"])

    def test_m6_from_planck_mass_relation(self):
        """M_6^4 = M_Pl^2 / (4*pi*a_0^2)."""
        from alpha_ladder_core.gauge_m6_parameter_space import compute_parameter_point, M_PL_EV, HBAR_C_EVM
        a_0_m = 28e-6
        result = compute_parameter_point(a_0_m)
        a_0_nat = a_0_m / HBAR_C_EVM
        expected_m6_4 = M_PL_EV ** 2 / (4.0 * math.pi * a_0_nat ** 2)
        expected_m6 = expected_m6_4 ** 0.25
        self.assertAlmostEqual(result["M_6_eV"] / expected_m6, 1.0, places=5)

    def test_m_kk_l1_is_sqrt2_times_m_kk(self):
        from alpha_ladder_core.gauge_m6_parameter_space import compute_parameter_point
        result = compute_parameter_point(28e-6)
        self.assertAlmostEqual(result["m_KK_l1_eV"],
                               math.sqrt(2.0) * result["m_KK_eV"], places=5)

    def test_r_phys_over_r_bare_is_e_phi(self):
        from alpha_ladder_core.gauge_m6_parameter_space import compute_parameter_point, E_PHI_VEV
        result = compute_parameter_point(28e-6)
        self.assertAlmostEqual(result["R_phys_over_R_bare"], E_PHI_VEV, places=10)

    def test_m_kk_enhancement_factor(self):
        from alpha_ladder_core.gauge_m6_parameter_space import compute_parameter_point, E_PHI_VEV
        result = compute_parameter_point(28e-6)
        self.assertAlmostEqual(result["m_KK_over_m_KK_bare"],
                               1.0 / E_PHI_VEV, places=10)

    def test_planck_scale_regime(self):
        from alpha_ladder_core.gauge_m6_parameter_space import compute_parameter_point, L_PL
        result = compute_parameter_point(L_PL)
        # At Planck length, M_6 is very high
        self.assertIn(result["regime"], ["Planck", "above_LHC"])

    def test_tev_regime(self):
        from alpha_ladder_core.gauge_m6_parameter_space import compute_parameter_point
        result = compute_parameter_point(28e-6)
        # At ~28 um, M_6 ~ few TeV
        self.assertIn(result["regime"], ["TeV", "above_LHC"])

    def test_a0_over_l_pl_positive(self):
        from alpha_ladder_core.gauge_m6_parameter_space import compute_parameter_point
        result = compute_parameter_point(28e-6)
        self.assertGreater(result["a_0_over_l_Pl"], 1.0)


class TestParameterSpaceScan(unittest.TestCase):
    """Test parameter_space_scan function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_m6_parameter_space import parameter_space_scan
        result = parameter_space_scan()
        self.assertIsInstance(result, dict)

    def test_has_points_and_summary(self):
        from alpha_ladder_core.gauge_m6_parameter_space import parameter_space_scan
        result = parameter_space_scan()
        self.assertIn("points", result)
        self.assertIn("summary", result)

    def test_default_30_points(self):
        from alpha_ladder_core.gauge_m6_parameter_space import parameter_space_scan
        result = parameter_space_scan()
        self.assertEqual(len(result["points"]), 30)

    def test_custom_values(self):
        from alpha_ladder_core.gauge_m6_parameter_space import parameter_space_scan
        vals = [1e-6, 10e-6, 100e-6]
        result = parameter_space_scan(vals)
        self.assertEqual(len(result["points"]), 3)

    def test_m6_monotonically_decreases_with_a0(self):
        """Larger a_0 => smaller M_6."""
        from alpha_ladder_core.gauge_m6_parameter_space import parameter_space_scan
        result = parameter_space_scan()
        pts = result["points"]
        for i in range(len(pts) - 1):
            if pts[i]["a_0_m"] < pts[i + 1]["a_0_m"]:
                self.assertGreater(pts[i]["M_6_eV"], pts[i + 1]["M_6_eV"])

    def test_summary_has_ranges(self):
        from alpha_ladder_core.gauge_m6_parameter_space import parameter_space_scan
        result = parameter_space_scan()
        s = result["summary"]
        for key in ["a_0_min_m", "a_0_max_m", "M_6_min_eV", "M_6_max_eV",
                     "m_KK_min_eV", "m_KK_max_eV"]:
            self.assertIn(key, s)


class TestNotablePoints(unittest.TestCase):
    """Test notable_points function."""

    def test_returns_list(self):
        from alpha_ladder_core.gauge_m6_parameter_space import notable_points
        result = notable_points()
        self.assertIsInstance(result, list)

    def test_has_11_entries(self):
        from alpha_ladder_core.gauge_m6_parameter_space import notable_points
        result = notable_points()
        self.assertEqual(len(result), 11)

    def test_each_has_description(self):
        from alpha_ladder_core.gauge_m6_parameter_space import notable_points
        result = notable_points()
        for pt in result:
            self.assertIn("description", pt)
            self.assertIsInstance(pt["description"], str)

    def test_planck_scale_first(self):
        from alpha_ladder_core.gauge_m6_parameter_space import notable_points, L_PL
        result = notable_points()
        self.assertAlmostEqual(result[0]["a_0_m"], L_PL, places=40)

    def test_eotwash_point_present(self):
        from alpha_ladder_core.gauge_m6_parameter_space import notable_points
        result = notable_points()
        eotwash = [p for p in result if "28.2" in p.get("description", "")]
        self.assertEqual(len(eotwash), 1)

    def test_1mm_point_present(self):
        from alpha_ladder_core.gauge_m6_parameter_space import notable_points
        result = notable_points()
        mm = [p for p in result if p.get("description", "").startswith("1 mm")]
        self.assertEqual(len(mm), 1)


class TestExperimentalReach(unittest.TestCase):
    """Test experimental_reach function."""

    def test_returns_list(self):
        from alpha_ladder_core.gauge_m6_parameter_space import experimental_reach
        result = experimental_reach()
        self.assertIsInstance(result, list)

    def test_has_at_least_5_experiments(self):
        from alpha_ladder_core.gauge_m6_parameter_space import experimental_reach
        result = experimental_reach()
        self.assertGreaterEqual(len(result), 5)

    def test_each_has_experiment_name(self):
        from alpha_ladder_core.gauge_m6_parameter_space import experimental_reach
        result = experimental_reach()
        for ex in result:
            self.assertIn("experiment", ex)
            self.assertIsInstance(ex["experiment"], str)

    def test_each_has_status(self):
        from alpha_ladder_core.gauge_m6_parameter_space import experimental_reach
        result = experimental_reach()
        for ex in result:
            self.assertIn("status", ex)
            self.assertIn(ex["status"], ["current", "planned"])

    def test_lhc_entry_exists(self):
        from alpha_ladder_core.gauge_m6_parameter_space import experimental_reach
        result = experimental_reach()
        lhc = [e for e in result if "LHC" in e["experiment"]]
        self.assertGreaterEqual(len(lhc), 1)

    def test_eotwash_entry_exists(self):
        from alpha_ladder_core.gauge_m6_parameter_space import experimental_reach
        result = experimental_reach()
        ew = [e for e in result if "Eot-Wash" in e["experiment"]]
        self.assertGreaterEqual(len(ew), 1)

    def test_cassini_has_no_threshold(self):
        from alpha_ladder_core.gauge_m6_parameter_space import experimental_reach
        result = experimental_reach()
        cassini = [e for e in result if "Cassini" in e["experiment"]]
        if cassini:
            self.assertIsNone(cassini[0]["a_0_threshold_m"])


class TestDilatonVevEffects(unittest.TestCase):
    """Test dilaton_vev_effects function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_m6_parameter_space import dilaton_vev_effects
        result = dilaton_vev_effects()
        self.assertIsInstance(result, dict)

    def test_phi_vev_is_negative(self):
        from alpha_ladder_core.gauge_m6_parameter_space import dilaton_vev_effects
        result = dilaton_vev_effects()
        self.assertLess(result["phi_vev"], 0)

    def test_radius_reduction_positive(self):
        from alpha_ladder_core.gauge_m6_parameter_space import dilaton_vev_effects
        result = dilaton_vev_effects()
        self.assertGreater(result["radius_reduction_percent"], 0)

    def test_m_kk_enhancement_greater_than_1(self):
        from alpha_ladder_core.gauge_m6_parameter_space import dilaton_vev_effects
        result = dilaton_vev_effects()
        self.assertGreater(result["m_KK_enhancement_factor"], 1.0)

    def test_m6_unchanged_by_vev(self):
        from alpha_ladder_core.gauge_m6_parameter_space import dilaton_vev_effects
        result = dilaton_vev_effects()
        self.assertTrue(result["M_6_unchanged"])

    def test_alpha_matched_is_alpha_em(self):
        from alpha_ladder_core.gauge_m6_parameter_space import dilaton_vev_effects, ALPHA_EM
        result = dilaton_vev_effects()
        self.assertAlmostEqual(result["alpha_matched"], ALPHA_EM, places=12)

    def test_alpha_unmatched_is_1_over_4pi(self):
        from alpha_ladder_core.gauge_m6_parameter_space import dilaton_vev_effects
        result = dilaton_vev_effects()
        self.assertAlmostEqual(result["alpha_unmatched"],
                               1.0 / (4.0 * math.pi), places=10)


class TestCanAnythingFixA0(unittest.TestCase):
    """Test can_anything_fix_a0 function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_m6_parameter_space import can_anything_fix_a0
        result = can_anything_fix_a0()
        self.assertIsInstance(result, dict)

    def test_8_mechanisms(self):
        from alpha_ladder_core.gauge_m6_parameter_space import can_anything_fix_a0
        result = can_anything_fix_a0()
        self.assertEqual(result["n_tried"], 8)

    def test_none_succeed(self):
        from alpha_ladder_core.gauge_m6_parameter_space import can_anything_fix_a0
        result = can_anything_fix_a0()
        self.assertEqual(result["n_succeed"], 0)

    def test_all_fail(self):
        from alpha_ladder_core.gauge_m6_parameter_space import can_anything_fix_a0
        result = can_anything_fix_a0()
        for mech in result["mechanisms"]:
            self.assertFalse(mech["can_fix_a0"])

    def test_each_has_verdict(self):
        from alpha_ladder_core.gauge_m6_parameter_space import can_anything_fix_a0
        result = can_anything_fix_a0()
        for mech in result["mechanisms"]:
            self.assertIn("verdict", mech)
            self.assertIsInstance(mech["verdict"], str)

    def test_overall_assessment_present(self):
        from alpha_ladder_core.gauge_m6_parameter_space import can_anything_fix_a0
        result = can_anything_fix_a0()
        self.assertIn("overall_assessment", result)
        self.assertIn("0 succeed", result["overall_assessment"])


class TestSummarizeParameterSpace(unittest.TestCase):
    """Test summarize_parameter_space function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_m6_parameter_space import summarize_parameter_space
        result = summarize_parameter_space()
        self.assertIsInstance(result, dict)

    def test_has_key_messages(self):
        from alpha_ladder_core.gauge_m6_parameter_space import summarize_parameter_space
        result = summarize_parameter_space()
        self.assertIn("key_messages", result)
        self.assertGreaterEqual(len(result["key_messages"]), 3)

    def test_has_sub_results(self):
        from alpha_ladder_core.gauge_m6_parameter_space import summarize_parameter_space
        result = summarize_parameter_space()
        for key in ["dilaton_vev_effects", "notable_points",
                     "experimental_reach", "mechanisms_to_fix_a0"]:
            self.assertIn(key, result)

    def test_notable_points_match(self):
        from alpha_ladder_core.gauge_m6_parameter_space import summarize_parameter_space, notable_points
        result = summarize_parameter_space()
        direct = notable_points()
        self.assertEqual(len(result["notable_points"]), len(direct))


class TestHelperFunctions(unittest.TestCase):
    """Test helper and formatting functions."""

    def test_classify_regime_planck(self):
        from alpha_ladder_core.gauge_m6_parameter_space import _classify_regime, EV_PER_TEV
        self.assertEqual(_classify_regime(2e16 * EV_PER_TEV), "Planck")

    def test_classify_regime_tev(self):
        from alpha_ladder_core.gauge_m6_parameter_space import _classify_regime, EV_PER_TEV
        self.assertEqual(_classify_regime(5 * EV_PER_TEV), "TeV")

    def test_classify_regime_gev(self):
        from alpha_ladder_core.gauge_m6_parameter_space import _classify_regime, EV_PER_GEV
        self.assertEqual(_classify_regime(5 * EV_PER_GEV), "GeV")

    def test_format_ev_tev(self):
        from alpha_ladder_core.gauge_m6_parameter_space import _format_eV, EV_PER_TEV
        result = _format_eV(5 * EV_PER_TEV)
        self.assertIn("TeV", result)

    def test_format_m_mm(self):
        from alpha_ladder_core.gauge_m6_parameter_space import _format_m
        result = _format_m(1e-3)
        self.assertIn("mm", result)


if __name__ == "__main__":
    unittest.main()
