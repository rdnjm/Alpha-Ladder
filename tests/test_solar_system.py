"""
Tests for alpha_ladder_core/solar_system.py
"""

import math
import pytest
from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.solar_system import (
    compute_ppn_parameters,
    compute_ppn_profile,
    compute_fifth_force_point,
    generate_exclusion_bounds,
    check_dilaton_exclusion,
    compute_minimum_dilaton_mass,
    summarize_solar_system,
)


@pytest.fixture
def constants():
    return get_constants()


@pytest.fixture
def omega():
    """The Brans-Dicke parameter from phi^2/2 bridge (~0.118)."""
    from alpha_ladder_core.dilaton import compute_bd_parameter
    bd = compute_bd_parameter(get_constants())
    return bd["omega"]


@pytest.fixture
def min_mass(omega):
    """Minimum dilaton mass from Cassini bound."""
    return compute_minimum_dilaton_mass(omega)


@pytest.fixture
def m_phi_eV(omega):
    """Dilaton mass from dilaton module (corrected)."""
    from alpha_ladder_core.dilaton import compute_bd_parameter
    bd = compute_bd_parameter(get_constants())
    return bd["dilaton_mass_min_eV"]


# -----------------------------------------------------------------------
# TestPPNParameters
# -----------------------------------------------------------------------
class TestPPNParameters:
    def test_gamma_massless_formula(self, omega, m_phi_eV):
        """gamma_massless = (1 + omega) / (2 + omega)."""
        result = compute_ppn_parameters(omega, m_phi_eV)
        expected = (1.0 + omega) / (2.0 + omega)
        assert pytest.approx(result["gamma_PPN_massless"], rel=1e-10) == expected

    def test_gamma_massless_approaches_1_for_large_omega(self):
        """For omega -> infinity (GR limit), gamma_massless -> 1."""
        result = compute_ppn_parameters(1e8, 1e-10)
        assert abs(result["gamma_PPN_massless"] - 1.0) < 1e-6

    def test_gamma_at_cassini_near_1(self, omega, m_phi_eV):
        """For massive dilaton, gamma at 1 AU should be very close to 1."""
        result = compute_ppn_parameters(omega, m_phi_eV)
        assert abs(result["gamma_PPN_at_cassini"] - 1.0) < 0.01

    def test_passes_cassini(self, omega, m_phi_eV):
        """The corrected dilaton mass should pass Cassini."""
        result = compute_ppn_parameters(omega, m_phi_eV)
        assert result["passes_cassini"] is True

    def test_gamma_deviation_below_cassini_bound(self, omega, m_phi_eV):
        """gamma deviation must be below 2.3e-5."""
        result = compute_ppn_parameters(omega, m_phi_eV)
        assert result["gamma_deviation_at_cassini"] < result["cassini_bound"]

    def test_returns_required_keys(self, omega, m_phi_eV):
        result = compute_ppn_parameters(omega, m_phi_eV)
        for key in ["omega_BD", "m_phi_eV", "beta_coupling", "alpha_fifth",
                     "lambda_compton_m", "gamma_PPN_massless", "gamma_PPN_at_cassini",
                     "gamma_deviation_at_cassini", "cassini_bound", "passes_cassini",
                     "required_min_mass_eV"]:
            assert key in result

    def test_alpha_fifth_positive(self, omega, m_phi_eV):
        result = compute_ppn_parameters(omega, m_phi_eV)
        assert result["alpha_fifth"] > 0

    def test_beta_coupling_positive(self, omega, m_phi_eV):
        result = compute_ppn_parameters(omega, m_phi_eV)
        assert result["beta_coupling"] > 0

    def test_alpha_fifth_equals_2_over_2w3(self, omega, m_phi_eV):
        result = compute_ppn_parameters(omega, m_phi_eV)
        expected = 2.0 / (2.0 * omega + 3.0)
        assert pytest.approx(result["alpha_fifth"], rel=1e-10) == expected

    def test_required_min_mass_positive(self, omega, m_phi_eV):
        result = compute_ppn_parameters(omega, m_phi_eV)
        assert result["required_min_mass_eV"] > 0


# -----------------------------------------------------------------------
# TestPPNProfile
# -----------------------------------------------------------------------
class TestPPNProfile:
    def test_correct_n_points(self, omega, m_phi_eV):
        result = compute_ppn_profile(omega, m_phi_eV, n_points=100)
        assert len(result["r_meters"]) == 100
        assert len(result["gamma_PPN"]) == 100
        assert len(result["gamma_deviation"]) == 100

    def test_default_n_points(self, omega, m_phi_eV):
        result = compute_ppn_profile(omega, m_phi_eV)
        assert len(result["r_meters"]) == 500

    def test_gamma_bounded(self, omega, m_phi_eV):
        """gamma should be between (1 - alpha_fifth) and 1."""
        result = compute_ppn_profile(omega, m_phi_eV)
        ppn = compute_ppn_parameters(omega, m_phi_eV)
        gamma_min = 1.0 - ppn["alpha_fifth"]  # limit as r -> 0
        for g in result["gamma_PPN"]:
            assert gamma_min - 1e-10 <= g <= 1.0 + 1e-10

    def test_gamma_monotonically_approaches_1(self, omega, m_phi_eV):
        """gamma deviation should decrease with distance."""
        result = compute_ppn_profile(omega, m_phi_eV)
        devs = result["gamma_deviation"]
        for i in range(1, len(devs)):
            assert devs[i] <= devs[i - 1] + 1e-15

    def test_transition_radius_exists(self, omega, m_phi_eV):
        """Transition radius should exist and be positive."""
        result = compute_ppn_profile(omega, m_phi_eV)
        tr = result["transition_radius_m"]
        assert tr is not None
        assert tr > 0

    def test_landmarks_present(self, omega, m_phi_eV):
        result = compute_ppn_profile(omega, m_phi_eV)
        landmarks = result["landmarks"]
        assert "lab" in landmarks
        assert "1_AU" in landmarks
        assert len(landmarks) >= 4

    def test_landmarks_have_gamma(self, omega, m_phi_eV):
        result = compute_ppn_profile(omega, m_phi_eV)
        for key, lm in result["landmarks"].items():
            assert "gamma_PPN" in lm
            assert "gamma_deviation" in lm
            assert "r_meters" in lm


# -----------------------------------------------------------------------
# TestFifthForce
# -----------------------------------------------------------------------
class TestFifthForce:
    def test_alpha_fifth_positive(self, omega, m_phi_eV):
        result = compute_fifth_force_point(omega, m_phi_eV)
        assert result["alpha_fifth"] > 0

    def test_lambda_fifth_positive(self, omega, m_phi_eV):
        result = compute_fifth_force_point(omega, m_phi_eV)
        assert result["lambda_fifth_m"] > 0

    def test_log10_values(self, omega, m_phi_eV):
        result = compute_fifth_force_point(omega, m_phi_eV)
        assert result["log10_alpha"] == pytest.approx(math.log10(result["alpha_fifth"]), rel=1e-10)
        assert result["log10_lambda_m"] == pytest.approx(math.log10(result["lambda_fifth_m"]), rel=1e-10)

    def test_bounds_have_names(self):
        bounds = generate_exclusion_bounds()
        for b in bounds:
            assert "name" in b
            assert isinstance(b["name"], str)
            assert len(b["name"]) > 0

    def test_at_least_4_bounds(self):
        bounds = generate_exclusion_bounds()
        assert len(bounds) >= 4

    def test_bounds_have_boundary_data(self):
        bounds = generate_exclusion_bounds()
        for b in bounds:
            assert "boundary" in b
            assert len(b["boundary"]) >= 3

    def test_dilaton_in_allowed_region(self, omega, m_phi_eV):
        result = check_dilaton_exclusion(omega, m_phi_eV)
        assert result["allowed"] is True

    def test_exclusion_returns_required_keys(self, omega, m_phi_eV):
        result = check_dilaton_exclusion(omega, m_phi_eV)
        assert "dilaton_point" in result
        assert "bounds" in result
        assert "excluded_by" in result
        assert "allowed" in result


# -----------------------------------------------------------------------
# TestMinimumDilatonMass
# -----------------------------------------------------------------------
class TestMinimumDilatonMass:
    def test_min_mass_greater_than_old_estimate(self, omega):
        """Corrected mass should be higher than the naive hbar/(c*AU)."""
        result = compute_minimum_dilaton_mass(omega)
        hbar = 1.054571817e-34
        c = 2.99792458e8
        AU = 1.496e11
        eV_to_J = 1.602176634e-19
        old_mass_eV = hbar / (c * AU) * c**2 / eV_to_J
        assert result["m_phi_min_eV"] > old_mass_eV

    def test_lambda_max_less_than_AU(self, omega):
        """Maximum Compton wavelength should be less than 1 AU."""
        result = compute_minimum_dilaton_mass(omega)
        AU = 1.496e11
        assert result["lambda_max_m"] < AU

    def test_min_mass_passes_cassini(self, omega):
        """At the minimum mass, gamma deviation should be at or below Cassini bound."""
        result = compute_minimum_dilaton_mass(omega)
        ppn = compute_ppn_parameters(omega, result["m_phi_min_eV"])
        # At the boundary, deviation should be approximately equal to Cassini bound
        assert ppn["gamma_deviation_at_cassini"] <= ppn["cassini_bound"] * 1.01

    def test_suppression_factor_greater_than_10(self, omega):
        result = compute_minimum_dilaton_mass(omega)
        assert result["suppression_factor"] > 10

    def test_returns_required_keys(self, omega):
        result = compute_minimum_dilaton_mass(omega)
        for key in ["m_phi_min_eV", "m_phi_min_kg", "lambda_max_m",
                     "lambda_max_au", "suppression_factor"]:
            assert key in result


# -----------------------------------------------------------------------
# TestSummarize
# -----------------------------------------------------------------------
class TestSummarize:
    def test_returns_dict(self, constants):
        result = summarize_solar_system(constants)
        assert isinstance(result, dict)

    def test_passes_all_bool(self, constants):
        result = summarize_solar_system(constants)
        assert isinstance(result["passes_all"], bool)

    def test_summary_string(self, constants):
        result = summarize_solar_system(constants)
        assert isinstance(result["summary"], str)
        assert len(result["summary"]) > 0

    def test_contains_ppn(self, constants):
        result = summarize_solar_system(constants)
        assert "ppn" in result
        assert isinstance(result["ppn"], dict)

    def test_contains_fifth_force(self, constants):
        result = summarize_solar_system(constants)
        assert "fifth_force" in result
        assert isinstance(result["fifth_force"], dict)
