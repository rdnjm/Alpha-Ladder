"""Tests for two_brane_physics.py -- two-brane physics on S^2."""
import unittest
import math



class TestModuleConstants(unittest.TestCase):
    """Test module-level constants."""

    def test_phi_vev_value(self):
        from alpha_ladder_core.gauge_two_brane_physics import PHI_VEV
        self.assertAlmostEqual(PHI_VEV, -0.5974, delta=0.01)

    def test_e_phi_vev_value(self):
        from alpha_ladder_core.gauge_two_brane_physics import E_PHI_VEV
        self.assertAlmostEqual(E_PHI_VEV, 0.5503, delta=0.01)

    def test_zeta_s2_minus_half(self):
        from alpha_ladder_core.gauge_two_brane_physics import ZETA_S2_MINUS_HALF
        self.assertAlmostEqual(ZETA_S2_MINUS_HALF, -0.25, places=10)


class TestLegendrePolynomials(unittest.TestCase):
    """Test _legendre_P helper."""

    def test_p0_is_1(self):
        from alpha_ladder_core.gauge_two_brane_physics import _legendre_P
        P = _legendre_P(0, 0.5)
        self.assertAlmostEqual(P[0], 1.0, places=10)

    def test_p1_is_x(self):
        from alpha_ladder_core.gauge_two_brane_physics import _legendre_P
        x = 0.7
        P = _legendre_P(1, x)
        self.assertAlmostEqual(P[1], x, places=10)

    def test_p2_formula(self):
        """P_2(x) = (3x^2 - 1) / 2."""
        from alpha_ladder_core.gauge_two_brane_physics import _legendre_P
        x = 0.6
        P = _legendre_P(2, x)
        expected = (3 * x ** 2 - 1) / 2
        self.assertAlmostEqual(P[2], expected, places=10)

    def test_p_l_at_1(self):
        """P_l(1) = 1 for all l."""
        from alpha_ladder_core.gauge_two_brane_physics import _legendre_P
        P = _legendre_P(10, 1.0)
        for l in range(11):
            self.assertAlmostEqual(P[l], 1.0, places=8)

    def test_p_l_at_minus1(self):
        """P_l(-1) = (-1)^l."""
        from alpha_ladder_core.gauge_two_brane_physics import _legendre_P
        P = _legendre_P(10, -1.0)
        for l in range(11):
            self.assertAlmostEqual(P[l], (-1.0) ** l, places=8)

    def test_empty_for_negative_lmax(self):
        from alpha_ladder_core.gauge_two_brane_physics import _legendre_P
        P = _legendre_P(-1, 0.5)
        self.assertEqual(len(P), 0)


class TestTwoBraneGeometry(unittest.TestCase):
    """Test two_brane_geometry function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_two_brane_physics import two_brane_geometry
        result = two_brane_geometry(1e48, 1e48, 3e12)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_two_brane_physics import two_brane_geometry
        result = two_brane_geometry(1e48, 1e48, 3e12)
        for key in ["delta_1_rad", "delta_2_rad", "Omega", "is_physical",
                     "M_Pl_eff_eV", "M_Pl_ratio"]:
            self.assertIn(key, result)

    def test_deficit_formula(self):
        """delta_i = 2*pi*T_i / M_6^4."""
        from alpha_ladder_core.gauge_two_brane_physics import two_brane_geometry
        T = 1e48
        M6 = 3e12
        result = two_brane_geometry(T, T, M6)
        expected_delta = 2 * math.pi * T / (M6 ** 4)
        self.assertAlmostEqual(result["delta_1_rad"], expected_delta, places=5)

    def test_omega_formula(self):
        """Omega = 4*pi - delta_1 - delta_2."""
        from alpha_ladder_core.gauge_two_brane_physics import two_brane_geometry
        T = 1e48
        M6 = 3e12
        result = two_brane_geometry(T, T, M6)
        expected = 4 * math.pi - result["delta_1_rad"] - result["delta_2_rad"]
        self.assertAlmostEqual(result["Omega"], expected, places=5)

    def test_small_tension_nearly_full_sphere(self):
        """For T << M_6^4, Omega ~ 4*pi."""
        from alpha_ladder_core.gauge_two_brane_physics import two_brane_geometry
        T = 1.0  # tiny tension
        M6 = 1e12
        result = two_brane_geometry(T, T, M6)
        self.assertAlmostEqual(result["fraction_of_sphere"], 1.0, places=5)

    def test_large_tension_unphysical(self):
        """If delta_1 + delta_2 > 4*pi, geometry is unphysical."""
        from alpha_ladder_core.gauge_two_brane_physics import two_brane_geometry
        M6 = 1e3
        T = M6 ** 4  # T = M_6^4
        result = two_brane_geometry(T, T, M6)
        # delta = 2*pi*T/M_6^4 = 2*pi per brane, total = 4*pi => Omega = 0
        self.assertFalse(result["is_physical"])

    def test_zero_m6_returns_error(self):
        from alpha_ladder_core.gauge_two_brane_physics import two_brane_geometry
        result = two_brane_geometry(1.0, 1.0, 0.0)
        self.assertIn("error", result)

    def test_fraction_near_1_for_small_tension(self):
        """Small tension means almost no deficit -- nearly full sphere."""
        from alpha_ladder_core.gauge_two_brane_physics import two_brane_geometry
        result = two_brane_geometry(1.0, 1.0, 1e12)
        self.assertAlmostEqual(result["fraction_of_sphere"], 1.0, places=5)


class TestSplitCouplings(unittest.TestCase):
    """Test split_couplings function."""

    def test_returns_list(self):
        from alpha_ladder_core.gauge_two_brane_physics import split_couplings
        result = split_couplings()
        self.assertIsInstance(result, list)

    def test_default_8_angles(self):
        from alpha_ladder_core.gauge_two_brane_physics import split_couplings
        result = split_couplings()
        self.assertEqual(len(result), 8)

    def test_single_angle(self):
        from alpha_ladder_core.gauge_two_brane_physics import split_couplings
        result = split_couplings(theta_2=math.pi / 4)
        self.assertEqual(len(result), 1)

    def test_no_splitting_at_north_pole(self):
        """At theta_2 = 0 (same pole), no splitting."""
        from alpha_ladder_core.gauge_two_brane_physics import split_couplings
        result = split_couplings(theta_2=0.0)
        self.assertAlmostEqual(result[0]["splitting_param"], 0.0, places=10)
        self.assertTrue(result[0]["couplings_equal_b2"])

    def test_no_splitting_at_south_pole(self):
        """At theta_2 = pi (rugby ball), no splitting."""
        from alpha_ladder_core.gauge_two_brane_physics import split_couplings
        result = split_couplings(theta_2=math.pi)
        self.assertAlmostEqual(result[0]["splitting_param"], 0.0, places=10)
        self.assertTrue(result[0]["couplings_equal_b2"])

    def test_max_splitting_at_equator(self):
        """At theta_2 = pi/2 (equator), cos^2(pi/2) = 0 => maximal splitting."""
        from alpha_ladder_core.gauge_two_brane_physics import split_couplings
        result = split_couplings(theta_2=math.pi / 2)
        self.assertAlmostEqual(result[0]["splitting_param"], 1.0, places=10)

    def test_brane1_always_equal(self):
        """Brane 1 at north pole always has equal couplings."""
        from alpha_ladder_core.gauge_two_brane_physics import split_couplings
        result = split_couplings()
        for s in result:
            self.assertTrue(s["couplings_equal_b1"])
            self.assertAlmostEqual(s["alpha_2_b1"], s["alpha_3_b1"], places=10)

    def test_k2_squared_is_cos2_theta(self):
        from alpha_ladder_core.gauge_two_brane_physics import split_couplings
        theta = math.pi / 3
        result = split_couplings(theta_2=theta)
        self.assertAlmostEqual(result[0]["K2_squared"],
                               math.cos(theta) ** 2, places=10)

    def test_k3_squared_is_1(self):
        """K3^2 = 1 for all non-polar angles."""
        from alpha_ladder_core.gauge_two_brane_physics import split_couplings
        result = split_couplings(theta_2=math.pi / 4)
        self.assertAlmostEqual(result[0]["K3_squared"], 1.0, places=10)

    def test_splitting_zero_at_poles_nonzero_elsewhere(self):
        from alpha_ladder_core.gauge_two_brane_physics import split_couplings
        result = split_couplings()
        for s in result:
            if abs(s["theta_2_rad"]) < 1e-10 or abs(s["theta_2_rad"] - math.pi) < 1e-10:
                self.assertAlmostEqual(s["splitting_param"], 0.0, places=10)
            else:
                self.assertGreater(s["splitting_param"], 0.0)


class TestInterBraneForce(unittest.TestCase):
    """Test inter_brane_force function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_two_brane_physics import inter_brane_force
        result = inter_brane_force(math.pi)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_two_brane_physics import inter_brane_force
        result = inter_brane_force(math.pi)
        for key in ["theta_2_rad", "R_phys_m", "d_arc_m", "V_results",
                     "m_lightest_KK_eV", "G_S2_angular"]:
            self.assertIn(key, result)

    def test_arc_distance_formula(self):
        """Arc distance = theta * R_phys."""
        from alpha_ladder_core.gauge_two_brane_physics import inter_brane_force
        theta = math.pi / 2
        result = inter_brane_force(theta)
        self.assertAlmostEqual(result["d_arc_m"],
                               theta * result["R_phys_m"], places=15)

    def test_exponential_suppression_at_macroscopic_distances(self):
        """At r = 1 cm, inter-brane force should be extremely suppressed."""
        from alpha_ladder_core.gauge_two_brane_physics import inter_brane_force
        result = inter_brane_force(math.pi)
        # Find 1 cm entry
        for vr in result["V_results"]:
            if vr["label"] == "1 cm":
                # V_4D should be extremely small or zero
                self.assertAlmostEqual(vr["V_4D_eV"], 0.0, places=30)

    def test_v_results_has_4_distances(self):
        from alpha_ladder_core.gauge_two_brane_physics import inter_brane_force
        result = inter_brane_force(math.pi)
        self.assertEqual(len(result["V_results"]), 4)

    def test_lightest_kk_mass_sqrt2_over_r(self):
        from alpha_ladder_core.gauge_two_brane_physics import inter_brane_force, HBAR_C_EVM
        result = inter_brane_force(math.pi)
        R_phys = result["R_phys_m"]
        R_eV_inv = R_phys / HBAR_C_EVM
        expected = math.sqrt(2.0) / R_eV_inv
        self.assertAlmostEqual(result["m_lightest_KK_eV"], expected, places=5)


class TestBraneAntibraneStabilization(unittest.TestCase):
    """Test brane_antibrane_stabilization function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_two_brane_physics import brane_antibrane_stabilization
        result = brane_antibrane_stabilization()
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_two_brane_physics import brane_antibrane_stabilization
        result = brane_antibrane_stabilization()
        for key in ["n_ratios_scanned", "n_with_minimum", "verdict",
                     "C_cas", "C_flux", "results"]:
            self.assertIn(key, result)

    def test_default_21_ratios(self):
        from alpha_ladder_core.gauge_two_brane_physics import brane_antibrane_stabilization
        result = brane_antibrane_stabilization()
        self.assertEqual(result["n_ratios_scanned"], 21)

    def test_custom_ratios(self):
        from alpha_ladder_core.gauge_two_brane_physics import brane_antibrane_stabilization
        result = brane_antibrane_stabilization(T_ratio_values=[-1.5, -1.0, -0.5])
        self.assertEqual(result["n_ratios_scanned"], 3)

    def test_c_cas_negative(self):
        """Casimir coefficient is negative (zeta_S2(-1/2) = -0.25)."""
        from alpha_ladder_core.gauge_two_brane_physics import brane_antibrane_stabilization
        result = brane_antibrane_stabilization()
        self.assertLess(result["C_cas"], 0)

    def test_c_flux_positive(self):
        from alpha_ladder_core.gauge_two_brane_physics import brane_antibrane_stabilization
        result = brane_antibrane_stabilization()
        self.assertGreater(result["C_flux"], 0)

    def test_each_result_has_t_ratio(self):
        from alpha_ladder_core.gauge_two_brane_physics import brane_antibrane_stabilization
        result = brane_antibrane_stabilization()
        for r in result["results"]:
            self.assertIn("T_ratio", r)
            self.assertIn("has_minimum", r)

    def test_has_note(self):
        from alpha_ladder_core.gauge_two_brane_physics import brane_antibrane_stabilization
        result = brane_antibrane_stabilization()
        self.assertIn("note", result)
        self.assertIn("speculative", result["note"])


class TestBrane2Physics(unittest.TestCase):
    """Test brane_2_physics function."""

    def test_returns_list(self):
        from alpha_ladder_core.gauge_two_brane_physics import brane_2_physics
        result = brane_2_physics()
        self.assertIsInstance(result, list)

    def test_default_8_angles(self):
        from alpha_ladder_core.gauge_two_brane_physics import brane_2_physics
        result = brane_2_physics()
        self.assertEqual(len(result), 8)

    def test_single_angle(self):
        from alpha_ladder_core.gauge_two_brane_physics import brane_2_physics
        result = brane_2_physics(theta_2=math.pi / 4)
        self.assertEqual(len(result), 1)

    def test_brane2_at_pole_same_physics(self):
        """At theta_2 = 0, Brane 2 has same physics as Brane 1."""
        from alpha_ladder_core.gauge_two_brane_physics import brane_2_physics
        result = brane_2_physics(theta_2=0.0)
        self.assertAlmostEqual(result[0]["alpha_ratio_b2_over_b1"], 1.0, places=10)
        self.assertAlmostEqual(result[0]["E_hydrogen_ratio"], 1.0, places=10)
        self.assertAlmostEqual(result[0]["bohr_radius_ratio"], 1.0, places=10)

    def test_brane2_at_south_pole_same_physics(self):
        """At theta_2 = pi, also same physics (rugby ball)."""
        from alpha_ladder_core.gauge_two_brane_physics import brane_2_physics
        result = brane_2_physics(theta_2=math.pi)
        self.assertAlmostEqual(result[0]["alpha_ratio_b2_over_b1"], 1.0, places=10)

    def test_equator_alpha_eff_zero(self):
        """At theta_2 = pi/2, alpha_2 = 0 => alpha_eff = 0 (pathological)."""
        from alpha_ladder_core.gauge_two_brane_physics import brane_2_physics
        result = brane_2_physics(theta_2=math.pi / 2)
        self.assertAlmostEqual(result[0]["alpha_eff_brane2"], 0.0, places=10)

    def test_equator_no_hydrogen(self):
        """At equator, no bound atoms possible."""
        from alpha_ladder_core.gauge_two_brane_physics import brane_2_physics
        result = brane_2_physics(theta_2=math.pi / 2)
        self.assertAlmostEqual(result[0]["E_hydrogen_eV_b2"], 0.0, places=10)

    def test_theta_60_weaker_binding(self):
        """At theta_2 = 60 deg, alpha_eff < alpha_EM => weaker binding."""
        from alpha_ladder_core.gauge_two_brane_physics import brane_2_physics
        result = brane_2_physics(theta_2=math.radians(60))
        self.assertLess(result[0]["alpha_ratio_b2_over_b1"], 1.0)
        self.assertLess(result[0]["E_hydrogen_ratio"], 1.0)
        self.assertGreater(result[0]["bohr_radius_ratio"], 1.0)

    def test_hydrogen_ratio_is_alpha_ratio_squared(self):
        """E_H ~ alpha^2, so E_H_ratio = alpha_ratio^2."""
        from alpha_ladder_core.gauge_two_brane_physics import brane_2_physics
        result = brane_2_physics(theta_2=math.radians(45))
        r = result[0]
        if r["alpha_ratio_b2_over_b1"] > 0:
            expected = r["alpha_ratio_b2_over_b1"] ** 2
            self.assertAlmostEqual(r["E_hydrogen_ratio"], expected, places=8)

    def test_bohr_ratio_is_inverse_alpha_ratio(self):
        """a_B ~ 1/alpha, so a_B_ratio = 1/alpha_ratio."""
        from alpha_ladder_core.gauge_two_brane_physics import brane_2_physics
        result = brane_2_physics(theta_2=math.radians(30))
        r = result[0]
        if r["alpha_ratio_b2_over_b1"] > 0:
            expected = 1.0 / r["alpha_ratio_b2_over_b1"]
            self.assertAlmostEqual(r["bohr_radius_ratio"], expected, places=8)

    def test_thomson_ratio_equals_alpha_ratio_squared(self):
        from alpha_ladder_core.gauge_two_brane_physics import brane_2_physics
        result = brane_2_physics(theta_2=math.radians(45))
        r = result[0]
        if r["alpha_ratio_b2_over_b1"] > 0:
            expected = r["alpha_ratio_b2_over_b1"] ** 2
            self.assertAlmostEqual(r["thomson_ratio"], expected, places=8)


class TestRugbyBallSpecialCase(unittest.TestCase):
    """Test rugby_ball_special_case function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_two_brane_physics import rugby_ball_special_case
        result = rugby_ball_special_case(1e48)
        self.assertIsInstance(result, dict)

    def test_required_keys(self):
        from alpha_ladder_core.gauge_two_brane_physics import rugby_ball_special_case
        result = rugby_ball_special_case(1e48)
        for key in ["geometry", "spectrum_shift", "M_Pl_ratio",
                     "V_dilaton_shift", "coupling_check", "is_physical"]:
            self.assertIn(key, result)

    def test_no_coupling_splitting(self):
        """Rugby ball: no splitting (both branes at poles)."""
        from alpha_ladder_core.gauge_two_brane_physics import rugby_ball_special_case
        result = rugby_ball_special_case(1e48)
        self.assertAlmostEqual(result["coupling_check"]["splitting_brane1"], 0.0)
        self.assertAlmostEqual(result["coupling_check"]["splitting_brane2"], 0.0)

    def test_total_deficit_is_2_delta(self):
        from alpha_ladder_core.gauge_two_brane_physics import rugby_ball_special_case
        result = rugby_ball_special_case(1e48)
        geom = result["geometry"]
        self.assertAlmostEqual(geom["total_deficit_rad"],
                               2.0 * geom["delta_rad"], places=10)

    def test_omega_formula(self):
        from alpha_ladder_core.gauge_two_brane_physics import rugby_ball_special_case
        result = rugby_ball_special_case(1e48)
        geom = result["geometry"]
        expected = 4 * math.pi - 2 * geom["delta_rad"]
        self.assertAlmostEqual(geom["Omega"], expected, places=5)

    def test_dilaton_potential_doubles(self):
        """Two equal branes => V_rugby = 2 * V_single."""
        from alpha_ladder_core.gauge_two_brane_physics import rugby_ball_special_case
        result = rugby_ball_special_case(1e48)
        self.assertAlmostEqual(result["V_dilaton_shift"]["ratio"], 2.0, places=10)

    def test_mass_shift_factor_greater_than_1(self):
        """Deficit makes KK masses heavier."""
        from alpha_ladder_core.gauge_two_brane_physics import rugby_ball_special_case
        T = (1e12) ** 4
        result = rugby_ball_special_case(T)
        if result["is_physical"]:
            self.assertGreaterEqual(
                result["spectrum_shift"]["mass_shift_factor"], 1.0)

    def test_spectrum_has_10_modes(self):
        from alpha_ladder_core.gauge_two_brane_physics import rugby_ball_special_case
        result = rugby_ball_special_case(1e48)
        self.assertEqual(len(result["spectrum_shift"]["modes"]), 10)

    def test_small_tension_physical(self):
        from alpha_ladder_core.gauge_two_brane_physics import rugby_ball_special_case
        result = rugby_ball_special_case(1.0)  # tiny tension
        self.assertTrue(result["is_physical"])

    def test_has_description(self):
        from alpha_ladder_core.gauge_two_brane_physics import rugby_ball_special_case
        result = rugby_ball_special_case(1e48)
        self.assertIn("description", result)
        self.assertIn("Rugby ball", result["description"])


class TestSummarizeTwoBrane(unittest.TestCase):
    """Test summarize_two_brane function."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_two_brane_physics import summarize_two_brane
        result = summarize_two_brane()
        self.assertIsInstance(result, dict)

    def test_has_report(self):
        from alpha_ladder_core.gauge_two_brane_physics import summarize_two_brane
        result = summarize_two_brane()
        self.assertIn("report", result)
        self.assertIsInstance(result["report"], str)
        self.assertGreater(len(result["report"]), 500)

    def test_has_all_sub_analyses(self):
        from alpha_ladder_core.gauge_two_brane_physics import summarize_two_brane
        result = summarize_two_brane()
        for key in ["geometry", "split_couplings", "inter_brane_force",
                     "brane_antibrane", "brane_2_physics", "rugby_ball"]:
            self.assertIn(key, result)

    def test_geometry_is_physical(self):
        from alpha_ladder_core.gauge_two_brane_physics import summarize_two_brane
        result = summarize_two_brane()
        self.assertTrue(result["geometry"]["is_physical"])

    def test_report_mentions_key_conclusions(self):
        from alpha_ladder_core.gauge_two_brane_physics import summarize_two_brane
        result = summarize_two_brane()
        self.assertIn("KEY CONCLUSIONS", result["report"])
        self.assertIn("COUPLING SPLITTING", result["report"])
        self.assertIn("RUGBY BALL", result["report"])


if __name__ == "__main__":
    unittest.main()
