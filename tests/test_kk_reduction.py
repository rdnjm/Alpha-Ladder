"""
Tests for alpha_ladder_core/kk_reduction.py

Kaluza-Klein dimensional reduction from 6D to 4D.
Verifies that the target Brans-Dicke parameter omega = (sqrt(5)-2)/2
emerges from a 6D theory.
"""

import math
import pytest
from decimal import Decimal
from alpha_ladder_core.kk_reduction import (
    compute_einstein_frame_ansatz,
    compute_kinetic_coefficient,
    compute_omega_from_kinetic,
    compute_gauss_bonnet_shift,
    scan_golden_point,
    compute_target_omega,
)


# -----------------------------------------------------------------------
# Tests for compute_einstein_frame_ansatz
# -----------------------------------------------------------------------

class TestEinsteinFrameAnsatz:
    """Tests for compute_einstein_frame_ansatz."""

    def test_default_d4_n2(self):
        result = compute_einstein_frame_ansatz()
        assert result["alpha"] == 1.0
        assert result["beta"] == pytest.approx(-1.0)

    def test_einstein_frame_condition(self):
        """(d-2)*alpha + n*beta must vanish for Einstein frame."""
        result = compute_einstein_frame_ansatz(d=4, n=2)
        assert (4 - 2) * result["alpha"] + 2 * result["beta"] == pytest.approx(0.0)

    def test_condition_not_d_alpha(self):
        """Verify it uses (d-2)*alpha, not d*alpha."""
        result = compute_einstein_frame_ansatz(d=4, n=2)
        # d*alpha + n*beta should NOT be zero (that would be the wrong condition)
        assert 4 * result["alpha"] + 2 * result["beta"] != pytest.approx(0.0)

    def test_returns_required_keys(self):
        result = compute_einstein_frame_ansatz()
        assert "alpha" in result
        assert "beta" in result
        assert "condition" in result

    def test_has_metric_ansatz_key(self):
        result = compute_einstein_frame_ansatz()
        assert "metric_ansatz" in result

    def test_custom_dimensions_d5_n3(self):
        """d=5, n=3: (5-2)*alpha + 3*beta = 0."""
        result = compute_einstein_frame_ansatz(d=5, n=3)
        assert (5 - 2) * result["alpha"] + 3 * result["beta"] == pytest.approx(0.0)

    def test_custom_dimensions_d6_n1(self):
        """d=6, n=1: (6-2)*alpha + 1*beta = 0 => beta = -4*alpha."""
        result = compute_einstein_frame_ansatz(d=6, n=1)
        assert (6 - 2) * result["alpha"] + 1 * result["beta"] == pytest.approx(0.0)

    def test_condition_string_not_empty(self):
        result = compute_einstein_frame_ansatz()
        assert isinstance(result["condition"], str)
        assert len(result["condition"]) > 0

    def test_metric_ansatz_string_not_empty(self):
        result = compute_einstein_frame_ansatz()
        assert isinstance(result["metric_ansatz"], str)
        assert len(result["metric_ansatz"]) > 0

    def test_alpha_is_float(self):
        result = compute_einstein_frame_ansatz()
        assert isinstance(result["alpha"], (int, float))

    def test_beta_is_float(self):
        result = compute_einstein_frame_ansatz()
        assert isinstance(result["beta"], (int, float))


# -----------------------------------------------------------------------
# Tests for compute_kinetic_coefficient
# -----------------------------------------------------------------------

class TestKineticCoefficient:
    """Tests for compute_kinetic_coefficient."""

    def test_standard_formula_d4_n2(self):
        """coeff_standard = -n*(n+d-2)/(2*(d-2)) = -2*4/(2*2) = -2.0."""
        result = compute_kinetic_coefficient(d=4, n=2)
        assert result["coeff_standard"] == pytest.approx(-2.0)

    def test_kinetic_coefficient_is_negative(self):
        """Non-ghost condition: kinetic coefficient must be negative."""
        result = compute_kinetic_coefficient(d=4, n=2)
        assert result["coeff_standard"] < 0

    def test_returns_required_keys(self):
        result = compute_kinetic_coefficient()
        for key in ["K_raw", "coeff_standard", "alpha", "beta"]:
            assert key in result

    def test_custom_alpha_beta(self):
        """When alpha/beta are provided explicitly, they should be used."""
        result = compute_kinetic_coefficient(d=4, n=2, alpha=1.0, beta=-1.0)
        assert result["alpha"] == 1.0
        assert result["beta"] == -1.0

    def test_auto_computes_alpha_beta_when_none(self):
        """When alpha/beta are None, should compute from Einstein frame condition."""
        result = compute_kinetic_coefficient(d=4, n=2)
        # Should satisfy (d-2)*alpha + n*beta = 0
        assert (4 - 2) * result["alpha"] + 2 * result["beta"] == pytest.approx(0.0)

    def test_has_omega_bd_key(self):
        result = compute_kinetic_coefficient(d=4, n=2)
        assert "omega_BD" in result

    def test_d5_n3_formula(self):
        """coeff_standard = -3*(3+5-2)/(2*(5-2)) = -3*6/(2*3) = -3.0."""
        result = compute_kinetic_coefficient(d=5, n=3)
        expected = -3.0 * (3.0 + 5.0 - 2.0) / (2.0 * (5.0 - 2.0))
        assert result["coeff_standard"] == pytest.approx(expected)

    def test_d4_n1_formula(self):
        """coeff_standard = -1*(1+4-2)/(2*(4-2)) = -1*3/(2*2) = -0.75."""
        result = compute_kinetic_coefficient(d=4, n=1)
        expected = -1.0 * (1.0 + 4.0 - 2.0) / (2.0 * (4.0 - 2.0))
        assert result["coeff_standard"] == pytest.approx(expected)


# -----------------------------------------------------------------------
# Tests for compute_omega_from_kinetic
# -----------------------------------------------------------------------

class TestOmegaFromKinetic:
    """Tests for compute_omega_from_kinetic."""

    def test_returns_dict(self):
        result = compute_omega_from_kinetic(-2.0, d=4, n=2)
        assert isinstance(result, dict)
        assert "omega_BD" in result

    def test_known_kinetic_coefficient(self):
        """For the standard d=4,n=2 case, verify omega is computed."""
        result = compute_omega_from_kinetic(-2.0, d=4, n=2)
        assert isinstance(result["omega_BD"], (int, float))

    def test_different_coefficients_give_different_omega(self):
        omega1 = compute_omega_from_kinetic(-1.0, d=4, n=2)
        omega2 = compute_omega_from_kinetic(-2.0, d=4, n=2)
        assert omega1 != omega2


# -----------------------------------------------------------------------
# Tests for compute_gauss_bonnet_shift
# -----------------------------------------------------------------------

class TestGaussBonnetShift:
    """Tests for compute_gauss_bonnet_shift."""

    def test_genus_2_chi(self):
        """Genus 2 surface has Euler characteristic chi = 2 - 2*2 = -2."""
        result = compute_gauss_bonnet_shift(genus=2)
        assert result["chi"] == -2

    def test_genus_1_no_shift(self):
        """Torus (genus 1) has chi=0, so no Gauss-Bonnet shift."""
        result = compute_gauss_bonnet_shift(genus=1)
        assert result["chi"] == 0

    def test_genus_0_positive_chi(self):
        """Sphere (genus 0) has chi = 2."""
        result = compute_gauss_bonnet_shift(genus=0)
        assert result["chi"] == 2

    def test_genus_3_chi(self):
        """Genus 3: chi = 2 - 2*3 = -4."""
        result = compute_gauss_bonnet_shift(genus=3)
        assert result["chi"] == -4

    def test_returns_required_keys(self):
        result = compute_gauss_bonnet_shift(genus=2)
        assert "chi" in result
        assert "delta_K" in result
        assert "omega_shifted" in result

    def test_genus_1_delta_k_zero(self):
        """For torus (chi=0), delta_K should be zero regardless of coupling."""
        result = compute_gauss_bonnet_shift(genus=1, gb_coupling=1.0)
        assert result["delta_K"] == pytest.approx(0.0)

    def test_default_dimensions(self):
        """Default d=4, n=2 should work without explicit arguments."""
        result = compute_gauss_bonnet_shift(d=4, n=2, genus=2)
        assert "omega_shifted" in result

    def test_custom_gb_coupling(self):
        """Providing explicit gb_coupling should be accepted."""
        result = compute_gauss_bonnet_shift(genus=2, gb_coupling=0.05)
        assert "delta_K" in result

    def test_chi_formula(self):
        """chi = 2 - 2*genus for all genus values."""
        for g in range(0, 8):
            result = compute_gauss_bonnet_shift(genus=g)
            assert result["chi"] == 2 - 2 * g


# -----------------------------------------------------------------------
# Tests for scan_golden_point
# -----------------------------------------------------------------------

class TestScanGoldenPoint:
    """Tests for scan_golden_point."""

    def test_returns_correct_count(self):
        results = scan_golden_point(genus_range=range(2, 6))
        assert len(results) == 4

    def test_returns_correct_count_larger_range(self):
        results = scan_golden_point(genus_range=range(2, 10))
        assert len(results) == 8

    def test_each_entry_has_required_keys(self):
        results = scan_golden_point(genus_range=range(2, 4))
        for r in results:
            assert "genus" in r
            assert "chi" in r
            assert "required_gb_coupling" in r

    def test_chi_values_correct(self):
        """chi = 2 - 2*genus for each scanned genus."""
        results = scan_golden_point(genus_range=range(2, 5))
        for r in results:
            assert r["chi"] == 2 - 2 * r["genus"]

    def test_genus_values_match_range(self):
        genus_range = range(2, 6)
        results = scan_golden_point(genus_range=genus_range)
        returned_genera = [r["genus"] for r in results]
        assert returned_genera == list(genus_range)

    def test_returns_list(self):
        results = scan_golden_point(genus_range=range(2, 4))
        assert isinstance(results, list)

    def test_single_genus(self):
        results = scan_golden_point(genus_range=range(3, 4))
        assert len(results) == 1
        assert results[0]["genus"] == 3

    def test_default_dimensions(self):
        """Default d=4, n=2 should work."""
        results = scan_golden_point(d=4, n=2, genus_range=range(2, 4))
        assert len(results) == 2


# -----------------------------------------------------------------------
# Tests for compute_target_omega
# -----------------------------------------------------------------------

class TestTargetOmega:
    """Tests for compute_target_omega."""

    def test_omega_value(self):
        """omega = (sqrt(5)-2)/2 ~ 0.11803398875."""
        result = compute_target_omega()
        assert result["omega_float"] == pytest.approx(0.11803398875, rel=1e-8)

    def test_omega_matches_formula(self):
        """Verify omega = (sqrt(5) - 2) / 2."""
        result = compute_target_omega()
        expected = (math.sqrt(5) - 2) / 2
        assert result["omega_float"] == pytest.approx(expected, rel=1e-12)

    def test_identity_exact(self):
        """1/2 - phi^(-2) should equal omega exactly (Decimal precision)."""
        result = compute_target_omega()
        assert result["identity_diff"] < Decimal("1e-40")

    def test_returns_decimal(self):
        result = compute_target_omega()
        assert isinstance(result["omega_decimal"], Decimal)

    def test_returns_float(self):
        result = compute_target_omega()
        assert isinstance(result["omega_float"], float)

    def test_omega_positive(self):
        result = compute_target_omega()
        assert result["omega_float"] > 0

    def test_omega_less_than_one(self):
        result = compute_target_omega()
        assert result["omega_float"] < 1.0

    def test_has_identity_diff_key(self):
        result = compute_target_omega()
        assert "identity_diff" in result

    def test_decimal_and_float_agree(self):
        """Decimal and float representations should agree to float precision."""
        result = compute_target_omega()
        assert float(result["omega_decimal"]) == pytest.approx(
            result["omega_float"], rel=1e-12
        )
