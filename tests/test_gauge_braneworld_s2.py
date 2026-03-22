"""Tests for braneworld_s2.py -- Braneworld Physics on S^2."""
import unittest
import math



# ---------------------------------------------------------------------------
# Constants used across tests (duplicated from module for validation)
# ---------------------------------------------------------------------------
ALPHA_EM = 1.0 / 137.035999084
PHI_VEV = 0.25 * math.log(4.0 * math.pi * ALPHA_EM)
E_PHI_VEV = math.exp(PHI_VEV)
ZETA_S2_MINUS_HALF = -0.25


class TestBraneTensionPotential(unittest.TestCase):
    """Tests for brane_tension_potential()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_tension_potential
        result = brane_tension_potential(T=1.0, phi=0.0)
        self.assertIsInstance(result, dict)

    def test_V_at_phi_zero(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_tension_potential
        result = brane_tension_potential(T=1.0, phi=0.0, n_branes=1)
        self.assertAlmostEqual(result["V"], 1.0, places=10)

    def test_V_functional_form(self):
        """V_brane = n*T*exp(4*phi)."""
        from alpha_ladder_core.gauge_braneworld_s2 import brane_tension_potential
        T, phi, n = 2.5, 0.3, 2
        result = brane_tension_potential(T=T, phi=phi, n_branes=n)
        expected = n * T * math.exp(4.0 * phi)
        self.assertAlmostEqual(result["V"], expected, places=10)

    def test_derivative(self):
        """dV/dphi = 4*n*T*exp(4*phi)."""
        from alpha_ladder_core.gauge_braneworld_s2 import brane_tension_potential
        T, phi, n = 3.0, -1.0, 1
        result = brane_tension_potential(T=T, phi=phi, n_branes=n)
        expected_dV = 4.0 * n * T * math.exp(4.0 * phi)
        self.assertAlmostEqual(result["dV_dphi"], expected_dV, places=10)

    def test_derivative_is_4_times_V(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_tension_potential
        result = brane_tension_potential(T=1.0, phi=0.5)
        self.assertAlmostEqual(result["dV_dphi"], 4.0 * result["V"], places=10)

    def test_two_branes(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_tension_potential
        r1 = brane_tension_potential(T=1.0, phi=0.0, n_branes=1)
        r2 = brane_tension_potential(T=1.0, phi=0.0, n_branes=2)
        self.assertAlmostEqual(r2["V"], 2.0 * r1["V"], places=10)

    def test_negative_tension(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_tension_potential
        result = brane_tension_potential(T=-1.0, phi=0.0)
        self.assertAlmostEqual(result["V"], -1.0, places=10)

    def test_keys_present(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_tension_potential
        result = brane_tension_potential(T=1.0, phi=0.0)
        for key in ["V", "dV_dphi", "n_branes", "T", "phi", "note"]:
            self.assertIn(key, result)

    def test_V_positive_for_positive_T(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_tension_potential
        for phi in [-2.0, -1.0, 0.0, 1.0, 2.0]:
            result = brane_tension_potential(T=1.0, phi=phi)
            self.assertGreater(result["V"], 0.0)


class TestCasimirPotential(unittest.TestCase):
    """Tests for casimir_potential()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_braneworld_s2 import casimir_potential
        result = casimir_potential(phi=0.0)
        self.assertIsInstance(result, dict)

    def test_zeta_value(self):
        from alpha_ladder_core.gauge_braneworld_s2 import casimir_potential
        result = casimir_potential(phi=0.0)
        self.assertEqual(result["zeta_value"], -0.25)

    def test_C_cas_negative(self):
        from alpha_ladder_core.gauge_braneworld_s2 import casimir_potential
        result = casimir_potential(phi=0.0, a_0=1.0)
        self.assertLess(result["C_cas"], 0.0)

    def test_C_cas_formula(self):
        """C_cas = zeta(-1/2) / (4*pi*a_0^4) = -1/(16*pi*a_0^4)."""
        from alpha_ladder_core.gauge_braneworld_s2 import casimir_potential
        a_0 = 2.0
        result = casimir_potential(phi=0.0, a_0=a_0)
        expected = -1.0 / (16.0 * math.pi * a_0 ** 4)
        self.assertAlmostEqual(result["C_cas"], expected, places=12)

    def test_V_cas_negative(self):
        """V_cas < 0 for all phi (since C_cas < 0 and exp(-4*phi) > 0)."""
        from alpha_ladder_core.gauge_braneworld_s2 import casimir_potential
        for phi in [-3.0, 0.0, 3.0]:
            result = casimir_potential(phi=phi)
            self.assertLess(result["V_cas"], 0.0)

    def test_V_cas_formula(self):
        from alpha_ladder_core.gauge_braneworld_s2 import casimir_potential
        phi = 1.5
        result = casimir_potential(phi=phi, a_0=1.0)
        expected = result["C_cas"] * math.exp(-4.0 * phi)
        self.assertAlmostEqual(result["V_cas"], expected, places=12)

    def test_derivative_sign(self):
        """dV_cas/dphi = -4*C_cas*exp(-4*phi). Since C_cas < 0, dV > 0."""
        from alpha_ladder_core.gauge_braneworld_s2 import casimir_potential
        result = casimir_potential(phi=0.0)
        self.assertGreater(result["dV_cas_dphi"], 0.0)

    def test_keys_present(self):
        from alpha_ladder_core.gauge_braneworld_s2 import casimir_potential
        result = casimir_potential(phi=0.0)
        for key in ["V_cas", "dV_cas_dphi", "C_cas", "zeta_value", "phi", "a_0"]:
            self.assertIn(key, result)


class TestCombinedPotential(unittest.TestCase):
    """Tests for combined_potential()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_braneworld_s2 import combined_potential
        result = combined_potential(phi=0.0, T=1.0)
        self.assertIsInstance(result, dict)

    def test_V_is_sum_of_components(self):
        from alpha_ladder_core.gauge_braneworld_s2 import combined_potential
        result = combined_potential(phi=0.5, T=1.0, N_flux=1)
        total = result["V_brane"] + result["V_casimir"] + result["V_flux"]
        self.assertAlmostEqual(result["V"], total, places=10)

    def test_no_flux(self):
        from alpha_ladder_core.gauge_braneworld_s2 import combined_potential
        result = combined_potential(phi=0.0, T=1.0, N_flux=0)
        self.assertAlmostEqual(result["V_flux"], 0.0, places=12)

    def test_brane_component_matches(self):
        from alpha_ladder_core.gauge_braneworld_s2 import combined_potential
        phi, T = 0.3, 2.0
        result = combined_potential(phi=phi, T=T, N_flux=0)
        expected_brane = T * math.exp(4.0 * phi)
        self.assertAlmostEqual(result["V_brane"], expected_brane, places=10)

    def test_no_minimum_positive_T_no_flux(self):
        """For T > 0, C_cas < 0, no flux: no analytic minimum (ratio < 0)."""
        from alpha_ladder_core.gauge_braneworld_s2 import combined_potential
        result = combined_potential(phi=0.0, T=1.0, N_flux=0)
        self.assertFalse(result["has_minimum_analytic_no_flux"])

    def test_minimum_negative_T_no_flux(self):
        """For T < 0, C_cas < 0: ratio = C_cas/T > 0, minimum exists."""
        from alpha_ladder_core.gauge_braneworld_s2 import combined_potential
        result = combined_potential(phi=0.0, T=-1.0, N_flux=0)
        self.assertTrue(result["has_minimum_analytic_no_flux"])
        self.assertIsNotNone(result["phi_min_analytic_no_flux"])

    def test_analytic_minimum_formula(self):
        """phi_min = ln(C_cas/T) / 8 for no-flux case with T < 0."""
        from alpha_ladder_core.gauge_braneworld_s2 import combined_potential
        T = -0.5
        result = combined_potential(phi=0.0, T=T, N_flux=0, a_0=1.0)
        C_cas = result["C_cas"]
        expected_phi = math.log(C_cas / T) / 8.0
        self.assertAlmostEqual(result["phi_min_analytic_no_flux"], expected_phi, places=10)

    def test_flux_coefficient(self):
        from alpha_ladder_core.gauge_braneworld_s2 import combined_potential
        result = combined_potential(phi=0.0, T=1.0, N_flux=2, a_0=1.0)
        expected_C_fl = 4.0 / (32.0 * math.pi ** 2)
        self.assertAlmostEqual(result["C_flux"], expected_C_fl, places=12)

    def test_keys_present(self):
        from alpha_ladder_core.gauge_braneworld_s2 import combined_potential
        result = combined_potential(phi=0.0, T=1.0)
        for key in ["V", "dV_dphi", "V_brane", "V_casimir", "V_flux",
                     "C_cas", "C_flux", "T", "phi", "N_flux"]:
            self.assertIn(key, result)


class TestFindBraneStabilization(unittest.TestCase):
    """Tests for find_brane_stabilization()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_braneworld_s2 import find_brane_stabilization
        result = find_brane_stabilization(a_0=1.0, N_flux=1)
        self.assertIsInstance(result, dict)

    def test_scans_many_T_values(self):
        from alpha_ladder_core.gauge_braneworld_s2 import find_brane_stabilization
        result = find_brane_stabilization()
        self.assertGreater(result["n_T_scanned"], 100)

    def test_finds_minima_with_flux(self):
        from alpha_ladder_core.gauge_braneworld_s2 import find_brane_stabilization
        result = find_brane_stabilization(a_0=1.0, N_flux=1)
        self.assertGreater(result["n_minima_found"], 0)

    def test_custom_T_values(self):
        from alpha_ladder_core.gauge_braneworld_s2 import find_brane_stabilization
        T_vals = [-1.0, -0.1, -0.01, 0.01, 0.1, 1.0]
        result = find_brane_stabilization(a_0=1.0, T_values=T_vals, N_flux=1)
        self.assertEqual(result["n_T_scanned"], len(T_vals))

    def test_phi_vev_target(self):
        from alpha_ladder_core.gauge_braneworld_s2 import find_brane_stabilization, PHI_VEV
        result = find_brane_stabilization()
        self.assertAlmostEqual(result["phi_vev_target"], PHI_VEV, places=6)

    def test_verdict_is_string(self):
        from alpha_ladder_core.gauge_braneworld_s2 import find_brane_stabilization
        result = find_brane_stabilization()
        self.assertIsInstance(result["verdict"], str)
        self.assertGreater(len(result["verdict"]), 10)

    def test_keys_present(self):
        from alpha_ladder_core.gauge_braneworld_s2 import find_brane_stabilization
        result = find_brane_stabilization()
        for key in ["n_T_scanned", "n_minima_found", "best_match",
                     "verdict", "phi_vev_target", "N_flux", "a_0"]:
            self.assertIn(key, result)


class TestDeficitAngle(unittest.TestCase):
    """Tests for deficit_angle()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_braneworld_s2 import deficit_angle
        result = deficit_angle(T=1.0, M_6_eV=1.0)
        self.assertIsInstance(result, dict)

    def test_delta_formula(self):
        """delta = 2*pi*T/M_6^4."""
        from alpha_ladder_core.gauge_braneworld_s2 import deficit_angle
        T, M_6 = 1e48, 1e12
        result = deficit_angle(T, M_6)
        expected = 2.0 * math.pi * T / M_6 ** 4
        self.assertAlmostEqual(result["delta_rad"], expected, places=6)

    def test_delta_degrees(self):
        from alpha_ladder_core.gauge_braneworld_s2 import deficit_angle
        result = deficit_angle(T=1e48, M_6_eV=1e12)
        self.assertAlmostEqual(
            result["delta_deg"], math.degrees(result["delta_rad"]), places=6
        )

    def test_zero_M6_returns_error(self):
        from alpha_ladder_core.gauge_braneworld_s2 import deficit_angle
        result = deficit_angle(T=1.0, M_6_eV=0.0)
        self.assertIn("error", result)

    def test_physical_range_small_T(self):
        """Small T/M_6^4 gives physical deficit (0 < delta < 2*pi)."""
        from alpha_ladder_core.gauge_braneworld_s2 import deficit_angle
        M_6 = 1e12
        T = 0.01 * M_6 ** 4
        result = deficit_angle(T, M_6)
        self.assertTrue(result["physical_single"])

    def test_effective_solid_angle_single_brane(self):
        from alpha_ladder_core.gauge_braneworld_s2 import deficit_angle
        M_6 = 1e12
        T = 0.1 * M_6 ** 4
        result = deficit_angle(T, M_6)
        expected = 4.0 * math.pi - result["delta_rad"]
        self.assertAlmostEqual(result["effective_solid_angle_1brane"], expected, places=6)

    def test_effective_solid_angle_rugby(self):
        from alpha_ladder_core.gauge_braneworld_s2 import deficit_angle
        M_6 = 1e12
        T = 0.1 * M_6 ** 4
        result = deficit_angle(T, M_6)
        expected = 4.0 * math.pi - 2.0 * result["delta_rad"]
        self.assertAlmostEqual(result["effective_solid_angle_2brane"], expected, places=6)

    def test_fraction_bounds(self):
        from alpha_ladder_core.gauge_braneworld_s2 import deficit_angle
        M_6 = 1e12
        T = 0.05 * M_6 ** 4
        result = deficit_angle(T, M_6)
        self.assertGreater(result["fraction_of_sphere_1brane"], 0.0)
        self.assertLess(result["fraction_of_sphere_1brane"], 1.0)


class TestBraneLocalizedCoupling(unittest.TestCase):
    """Tests for brane_localized_coupling()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_localized_coupling
        result = brane_localized_coupling()
        self.assertIsInstance(result, dict)

    def test_alpha_brane_matches_alpha_EM(self):
        """alpha_brane = exp(4*phi_vev)/(4*pi) = alpha_EM by construction."""
        from alpha_ladder_core.gauge_braneworld_s2 import brane_localized_coupling
        result = brane_localized_coupling()
        self.assertAlmostEqual(
            result["alpha_brane"], result["alpha_EM"], places=12
        )

    def test_alpha_match_ppm_near_zero(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_localized_coupling
        result = brane_localized_coupling()
        self.assertAlmostEqual(result["alpha_match_ppm"], 0.0, places=3)

    def test_custom_phi_vev(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_localized_coupling
        result = brane_localized_coupling(phi_vev=0.0)
        expected = math.exp(0.0) / (4.0 * math.pi)
        self.assertAlmostEqual(result["alpha_brane"], expected, places=12)

    def test_G_4D_positive(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_localized_coupling
        result = brane_localized_coupling()
        self.assertGreater(result["G_4D"], 0.0)

    def test_keys_present(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_localized_coupling
        result = brane_localized_coupling()
        for key in ["phi_vev", "alpha_brane", "alpha_EM",
                     "alpha_match_ppm", "G_4D"]:
            self.assertIn(key, result)


class TestBraneMatterSpectrum(unittest.TestCase):
    """Tests for brane_matter_spectrum()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_matter_spectrum
        result = brane_matter_spectrum()
        self.assertIsInstance(result, dict)

    def test_10_graviton_modes(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_matter_spectrum
        result = brane_matter_spectrum()
        self.assertEqual(result["n_graviton_modes"], 10)

    def test_10_gauge_modes(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_matter_spectrum
        result = brane_matter_spectrum()
        self.assertEqual(result["n_gauge_modes"], 10)

    def test_graviton_starts_at_l2(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_matter_spectrum
        result = brane_matter_spectrum()
        self.assertEqual(result["graviton_modes"][0]["l"], 2)

    def test_gauge_starts_at_l1(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_matter_spectrum
        result = brane_matter_spectrum()
        self.assertEqual(result["gauge_modes"][0]["l"], 1)

    def test_graviton_mass_formula(self):
        """m_l = sqrt(l(l+1)-2) * hbar_c / R_phys."""
        from alpha_ladder_core.gauge_braneworld_s2 import brane_matter_spectrum, HBAR_C_EVM
        a_0_m = 28e-6
        result = brane_matter_spectrum(a_0_m=a_0_m)
        R = result["R_phys_m"]
        l = 2
        expected = math.sqrt(l * (l + 1) - 2) * HBAR_C_EVM / R
        self.assertAlmostEqual(
            result["graviton_modes"][0]["m_eV"], expected, places=5
        )

    def test_gauge_mass_formula(self):
        """m_l = sqrt(l(l+1)) * hbar_c / R_phys."""
        from alpha_ladder_core.gauge_braneworld_s2 import brane_matter_spectrum, HBAR_C_EVM
        a_0_m = 28e-6
        result = brane_matter_spectrum(a_0_m=a_0_m)
        R = result["R_phys_m"]
        l = 1
        expected = math.sqrt(l * (l + 1)) * HBAR_C_EVM / R
        self.assertAlmostEqual(
            result["gauge_modes"][0]["m_eV"], expected, places=5
        )

    def test_lightest_gauge_lighter_than_graviton(self):
        """Gauge l=1 (sqrt(2)) < graviton l=2 (sqrt(4)=2)."""
        from alpha_ladder_core.gauge_braneworld_s2 import brane_matter_spectrum
        result = brane_matter_spectrum()
        self.assertLess(
            result["gauge_modes"][0]["m_eV"],
            result["graviton_modes"][0]["m_eV"],
        )

    def test_lightest_bulk_mode_is_gauge(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_matter_spectrum
        result = brane_matter_spectrum()
        self.assertEqual(result["lightest_bulk_mode"]["l"], 1)

    def test_degeneracy(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_matter_spectrum
        result = brane_matter_spectrum()
        for mode in result["graviton_modes"]:
            self.assertEqual(mode["degeneracy"], 2 * mode["l"] + 1)

    def test_R_phys_positive(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_matter_spectrum
        result = brane_matter_spectrum()
        self.assertGreater(result["R_phys_m"], 0.0)

    def test_M_6_positive(self):
        from alpha_ladder_core.gauge_braneworld_s2 import brane_matter_spectrum
        result = brane_matter_spectrum()
        self.assertGreater(result["M_6_eV"], 0.0)


class TestDoesBraneCloseGap1(unittest.TestCase):
    """Tests for does_brane_close_gap1()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_braneworld_s2 import does_brane_close_gap1
        result = does_brane_close_gap1()
        self.assertIsInstance(result, dict)

    def test_gap1_not_closed(self):
        """The brane mechanism does NOT close Gap 1."""
        from alpha_ladder_core.gauge_braneworld_s2 import does_brane_close_gap1
        result = does_brane_close_gap1()
        self.assertFalse(result["polynomial_emerges"])

    def test_gap1_status_open(self):
        from alpha_ladder_core.gauge_braneworld_s2 import does_brane_close_gap1
        result = does_brane_close_gap1()
        self.assertEqual(result["gap1_status"], "OPEN")

    def test_target_polynomial(self):
        from alpha_ladder_core.gauge_braneworld_s2 import does_brane_close_gap1
        result = does_brane_close_gap1()
        self.assertEqual(result["target_polynomial"], "x^2 + 6x + 4 = 0")

    def test_target_roots(self):
        from alpha_ladder_core.gauge_braneworld_s2 import does_brane_close_gap1
        result = does_brane_close_gap1()
        x_plus = -3.0 + math.sqrt(5.0)
        x_minus = -3.0 - math.sqrt(5.0)
        self.assertAlmostEqual(result["target_roots"]["x_plus"], x_plus, places=6)
        self.assertAlmostEqual(result["target_roots"]["x_minus"], x_minus, places=6)

    def test_x_plus_is_negative(self):
        """x = -3+sqrt(5) < 0, cannot be exp(2*phi)."""
        from alpha_ladder_core.gauge_braneworld_s2 import does_brane_close_gap1
        result = does_brane_close_gap1()
        self.assertLess(result["target_roots"]["x_plus"], 0.0)

    def test_reason_is_string(self):
        from alpha_ladder_core.gauge_braneworld_s2 import does_brane_close_gap1
        result = does_brane_close_gap1()
        self.assertIsInstance(result["reason"], str)
        self.assertGreater(len(result["reason"]), 50)

    def test_C_cas_present(self):
        from alpha_ladder_core.gauge_braneworld_s2 import does_brane_close_gap1
        result = does_brane_close_gap1()
        self.assertIn("C_cas", result)
        self.assertLess(result["C_cas"], 0.0)


class TestSummarizeBraneworld(unittest.TestCase):
    """Tests for summarize_braneworld()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_braneworld_s2 import summarize_braneworld
        result = summarize_braneworld()
        self.assertIsInstance(result, dict)

    def test_all_sections_present(self):
        from alpha_ladder_core.gauge_braneworld_s2 import summarize_braneworld
        result = summarize_braneworld()
        for key in ["brane_potential", "casimir_potential", "combined_potential",
                     "stabilization_scan", "deficit_angle", "brane_coupling",
                     "matter_spectrum", "gap1_analysis", "report"]:
            self.assertIn(key, result)

    def test_report_is_string(self):
        from alpha_ladder_core.gauge_braneworld_s2 import summarize_braneworld
        result = summarize_braneworld()
        self.assertIsInstance(result["report"], str)
        self.assertGreater(len(result["report"]), 500)

    def test_report_mentions_gap1(self):
        from alpha_ladder_core.gauge_braneworld_s2 import summarize_braneworld
        result = summarize_braneworld()
        self.assertIn("GAP 1", result["report"])

    def test_gap1_analysis_not_closed(self):
        from alpha_ladder_core.gauge_braneworld_s2 import summarize_braneworld
        result = summarize_braneworld()
        self.assertFalse(result["gap1_analysis"]["polynomial_emerges"])


class TestHelpers(unittest.TestCase):
    """Tests for helper functions and module constants."""

    def test_phi_vev_value(self):
        from alpha_ladder_core.gauge_braneworld_s2 import PHI_VEV
        expected = 0.25 * math.log(4.0 * math.pi / 137.035999084)
        self.assertAlmostEqual(PHI_VEV, expected, places=10)

    def test_phi_vev_negative(self):
        from alpha_ladder_core.gauge_braneworld_s2 import PHI_VEV
        self.assertLess(PHI_VEV, 0.0)

    def test_zeta_value(self):
        from alpha_ladder_core.gauge_braneworld_s2 import ZETA_S2_MINUS_HALF
        self.assertEqual(ZETA_S2_MINUS_HALF, -0.25)

    def test_R_phys_helper(self):
        from alpha_ladder_core.gauge_braneworld_s2 import _R_phys, E_PHI_VEV
        a_0 = 28e-6
        self.assertAlmostEqual(_R_phys(a_0), a_0 * E_PHI_VEV, places=15)

    def test_C_casimir_helper(self):
        from alpha_ladder_core.gauge_braneworld_s2 import _C_casimir
        a_0 = 1.0
        expected = -0.25 / (4.0 * math.pi)
        self.assertAlmostEqual(_C_casimir(a_0), expected, places=12)

    def test_C_flux_zero_N(self):
        from alpha_ladder_core.gauge_braneworld_s2 import _C_flux
        self.assertEqual(_C_flux(0, 1.0), 0.0)

    def test_C_flux_positive(self):
        from alpha_ladder_core.gauge_braneworld_s2 import _C_flux
        self.assertGreater(_C_flux(1, 1.0), 0.0)


if __name__ == "__main__":
    unittest.main()
