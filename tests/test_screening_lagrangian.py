"""Tests for the screening Lagrangian derivation functions."""

import math
import pytest
from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.screening import (
    derive_dilaton_lagrangian,
    derive_field_equation,
    derive_yukawa_profile,
    verify_screening_consistency,
    compute_screening_parameters,
)


@pytest.fixture
def constants():
    return get_constants("CODATA 2018")


@pytest.fixture
def screening_params(constants):
    return compute_screening_parameters(constants)


class TestDeriveDilatonLagrangian:
    """Tests for derive_dilaton_lagrangian."""

    def test_returns_dict(self):
        result = derive_dilaton_lagrangian(omega_BD=1000)
        assert isinstance(result, dict)

    def test_has_action_terms(self):
        result = derive_dilaton_lagrangian(omega_BD=1000)
        assert "action_terms" in result
        assert isinstance(result["action_terms"], list)

    def test_four_action_terms(self):
        result = derive_dilaton_lagrangian(omega_BD=1000)
        assert len(result["action_terms"]) == 4

    def test_beta_coupling_positive_omega(self):
        """For large positive omega=1000, beta = 1/sqrt(2003)."""
        result = derive_dilaton_lagrangian(omega_BD=1000)
        expected_beta = 1.0 / math.sqrt(2 * 1000 + 3)
        assert result["beta_coupling"] == pytest.approx(expected_beta, rel=1e-10)
        assert result["beta_coupling"] == pytest.approx(0.02235, rel=1e-3)

    def test_beta_coupling_formula(self):
        """For omega=0, beta = 1/sqrt(3)."""
        result = derive_dilaton_lagrangian(omega_BD=0)
        expected_beta = 1.0 / math.sqrt(3)
        assert result["beta_coupling"] == pytest.approx(expected_beta, rel=1e-10)
        assert result["beta_coupling"] == pytest.approx(0.5774, rel=1e-3)

    def test_omega_bd_preserved(self):
        omega_val = 42.5
        result = derive_dilaton_lagrangian(omega_BD=omega_val)
        assert result["omega_BD"] == omega_val

    def test_has_conformal_function(self):
        result = derive_dilaton_lagrangian(omega_BD=100)
        assert "conformal_function" in result

    def test_has_canonical_normalization(self):
        result = derive_dilaton_lagrangian(omega_BD=100)
        assert "canonical_normalization" in result

    def test_mass_none_symbolic(self):
        """When m_phi_eV is None, potential term should be symbolic (no numeric mass)."""
        result = derive_dilaton_lagrangian(omega_BD=100, m_phi_eV=None)
        # Find the potential term among the action terms
        potential_terms = [t for t in result["action_terms"] if "potential" in t.get("name", "").lower()
                          or "potential" in t.get("type", "").lower()
                          or "V(" in str(t.get("expression", ""))]
        assert len(potential_terms) > 0
        # The potential term should indicate symbolic / no numeric mass
        pot = potential_terms[0]
        pot_str = str(pot)
        # Should not contain a specific numeric mass value
        assert "1e-05" not in pot_str or "symbolic" in pot_str.lower() or "m_phi" in pot_str

    def test_mass_provided(self):
        """When m_phi_eV is given, potential should reference the numeric value."""
        result = derive_dilaton_lagrangian(omega_BD=100, m_phi_eV=1e-5)
        potential_terms = [t for t in result["action_terms"] if "potential" in t.get("name", "").lower()
                          or "potential" in t.get("type", "").lower()
                          or "V(" in str(t.get("expression", ""))]
        assert len(potential_terms) > 0
        pot = potential_terms[0]
        pot_str = str(pot)
        # Should contain numeric mass info
        assert "1e-05" in pot_str or "1.0e-05" in pot_str or pot.get("m_phi_eV") == pytest.approx(1e-5)


class TestDeriveFieldEquation:
    """Tests for derive_field_equation."""

    def test_returns_dict(self):
        result = derive_field_equation(omega_BD=100)
        assert isinstance(result, dict)

    def test_has_equation(self):
        result = derive_field_equation(omega_BD=100)
        assert "equation" in result
        assert isinstance(result["equation"], str)

    def test_has_green_function(self):
        result = derive_field_equation(omega_BD=100)
        assert "green_function" in result
        assert isinstance(result["green_function"], dict)

    def test_beta_coupling_matches_lagrangian(self):
        """Beta from field equation should match beta from Lagrangian for the same omega."""
        omega = 500
        lagrangian_result = derive_dilaton_lagrangian(omega_BD=omega)
        field_eq_result = derive_field_equation(omega_BD=omega)
        assert field_eq_result["beta_coupling"] == pytest.approx(
            lagrangian_result["beta_coupling"], rel=1e-10
        )

    def test_has_linearization_condition(self):
        result = derive_field_equation(omega_BD=100)
        assert "linearization_condition" in result
        assert isinstance(result["linearization_condition"], str)


class TestDeriveYukawaProfile:
    """Tests for derive_yukawa_profile."""

    def test_returns_dict(self):
        result = derive_yukawa_profile(omega_BD=100)
        assert isinstance(result, dict)

    def test_has_derivation_steps(self):
        result = derive_yukawa_profile(omega_BD=100)
        assert "derivation_steps" in result
        assert isinstance(result["derivation_steps"], list)

    def test_five_derivation_steps(self):
        result = derive_yukawa_profile(omega_BD=100)
        assert len(result["derivation_steps"]) == 5

    def test_alpha_screening_formula(self):
        """alpha_screening_from_lagrangian = 2 * beta^2."""
        omega = 100
        result = derive_yukawa_profile(omega_BD=omega)
        beta = 1.0 / math.sqrt(2 * omega + 3)
        expected_alpha = 2.0 * beta ** 2
        assert result["alpha_screening_from_lagrangian"] == pytest.approx(expected_alpha, rel=1e-10)

    def test_alpha_screening_positive(self):
        """Alpha screening should be positive for any real omega where 2*omega+3 != 0."""
        for omega in [-1.0, 0.0, 0.118, 1.0, 100.0, 40000.0]:
            result = derive_yukawa_profile(omega_BD=omega)
            assert result["alpha_screening_from_lagrangian"] > 0, (
                f"alpha_screening should be positive for omega={omega}"
            )

    def test_has_yukawa_formula(self):
        result = derive_yukawa_profile(omega_BD=100)
        assert "yukawa_formula" in result
        assert isinstance(result["yukawa_formula"], str)

    def test_large_omega_small_screening(self):
        """For Cassini-bound omega=40000, alpha_screening should be very small."""
        result = derive_yukawa_profile(omega_BD=40000)
        # beta = 1/sqrt(80003) ~ 0.003536, alpha = 2*beta^2 ~ 2.5e-5
        assert result["alpha_screening_from_lagrangian"] == pytest.approx(2.5e-5, rel=0.01)

    def test_omega_zero_screening(self):
        """For omega=0, beta=1/sqrt(3), alpha=2/3."""
        result = derive_yukawa_profile(omega_BD=0)
        assert result["alpha_screening_from_lagrangian"] == pytest.approx(2.0 / 3.0, rel=1e-10)


class TestVerifyScreeningConsistency:
    """Tests for verify_screening_consistency."""

    def test_returns_dict(self, constants):
        result = verify_screening_consistency(constants)
        assert isinstance(result, dict)

    def test_has_all_keys(self, constants):
        result = verify_screening_consistency(constants)
        expected_keys = [
            "omega_BD",
            "beta_coupling",
            "alpha_screening_lagrangian",
            "alpha_screening_empirical",
            "ratio",
            "consistent",
            "interpretation",
        ]
        for key in expected_keys:
            assert key in result, f"Missing key: {key}"

    def test_alpha_screening_lagrangian_positive(self, constants):
        result = verify_screening_consistency(constants)
        assert result["alpha_screening_lagrangian"] > 0

    def test_alpha_screening_empirical_positive(self, constants):
        result = verify_screening_consistency(constants)
        assert result["alpha_screening_empirical"] > 0

    def test_interpretation_string(self, constants):
        result = verify_screening_consistency(constants)
        assert isinstance(result["interpretation"], str)
        assert len(result["interpretation"]) > 0

    def test_omega_bd_matches_dilaton(self, constants):
        """omega_BD should match the value from the dilaton module."""
        from alpha_ladder_core.dilaton import compute_bd_parameter

        bd = compute_bd_parameter(constants)
        result = verify_screening_consistency(constants)
        assert result["omega_BD"] == pytest.approx(bd["omega"], rel=1e-10)
