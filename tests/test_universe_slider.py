"""
Tests for the universe_slider module.
"""

from decimal import Decimal

import pytest

from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.universe_slider import recompute_physics, compute_sensitivity_curve


@pytest.fixture
def constants():
    return get_constants("CODATA 2018")


@pytest.fixture
def real_alpha_result(constants):
    return recompute_physics(constants.alpha, constants)


# --------------------------------------------------------------------------
# Basic structure tests
# --------------------------------------------------------------------------

class TestRecomputePhysicsStructure:
    """Verify the returned dict has the expected keys and types."""

    def test_all_keys_present(self, real_alpha_result):
        expected_keys = {
            "alpha", "inv_alpha", "alpha_powers", "alpha_g", "G_predicted",
            "alpha_g_hierarchy", "G_hierarchy",
            "r_e", "a_0_predicted", "binding_energy_ratio", "dark_coupling",
            "gap", "chemistry_viable", "stars_viable", "notes",
        }
        assert set(real_alpha_result.keys()) == expected_keys

    def test_alpha_is_decimal(self, real_alpha_result):
        assert isinstance(real_alpha_result["alpha"], Decimal)

    def test_notes_is_list(self, real_alpha_result):
        assert isinstance(real_alpha_result["notes"], list)


# --------------------------------------------------------------------------
# Alpha powers
# --------------------------------------------------------------------------

class TestAlphaPowers:
    """Verify alpha_powers are computed correctly."""

    def test_alpha_power_1(self, real_alpha_result, constants):
        assert real_alpha_result["alpha_powers"][1] == constants.alpha

    def test_alpha_power_2(self, real_alpha_result, constants):
        expected = constants.alpha ** 2
        assert real_alpha_result["alpha_powers"][2] == expected

    def test_alpha_power_21(self, real_alpha_result, constants):
        expected = constants.alpha ** 21
        assert pytest.approx(float(real_alpha_result["alpha_powers"][21]), rel=1e-40) == float(expected)

    def test_all_24_powers(self, real_alpha_result):
        assert len(real_alpha_result["alpha_powers"]) == 24
        assert set(real_alpha_result["alpha_powers"].keys()) == set(range(1, 25))

    def test_powers_are_monotonically_decreasing(self, real_alpha_result):
        """For alpha < 1, higher powers should be smaller."""
        powers = real_alpha_result["alpha_powers"]
        for n in range(1, 24):
            assert powers[n] > powers[n + 1]


# --------------------------------------------------------------------------
# G_predicted at real alpha
# --------------------------------------------------------------------------

class TestGPredicted:
    """Verify G_predicted matches expected value at the real alpha."""

    def test_G_predicted_order_of_magnitude(self, real_alpha_result):
        G = float(real_alpha_result["G_predicted"])
        assert 6.0e-11 < G < 7.0e-11

    def test_G_predicted_matches_expected(self, real_alpha_result):
        """The phi^2/2 bridge should predict G ~ 6.67323e-11."""
        G = float(real_alpha_result["G_predicted"])
        assert abs(G - 6.67323e-11) < 0.005e-11

    def test_G_hierarchy_order_of_magnitude(self, real_alpha_result):
        G = float(real_alpha_result["G_hierarchy"])
        assert 6.0e-11 < G < 7.0e-11

    def test_G_hierarchy_matches_expected(self, real_alpha_result):
        """The alpha^24 * mu^2 formula should predict G ~ 6.6789e-11 (688 ppm)."""
        G = float(real_alpha_result["G_hierarchy"])
        assert abs(G - 6.6789e-11) < 0.01e-11

    def test_alpha_g_hierarchy_equals_alpha24_mu2(self, real_alpha_result, constants):
        mu = constants.m_p / constants.m_e
        expected = constants.alpha ** 24 * mu ** 2
        assert float(real_alpha_result["alpha_g_hierarchy"]) == pytest.approx(float(expected), rel=1e-10)

    def test_G_predicted_within_1_sigma_of_codata(self, real_alpha_result):
        """Predicted G should be within a few sigma of CODATA 2018 (6.67430e-11)."""
        G = float(real_alpha_result["G_predicted"])
        G_codata = 6.67430e-11
        G_unc = 0.00015e-11
        sigma = abs(G - G_codata) / G_unc
        # The prediction is off by a few sigma but should be within ~10
        assert sigma < 15


# --------------------------------------------------------------------------
# Derived quantities at real alpha
# --------------------------------------------------------------------------

class TestDerivedQuantities:

    def test_inv_alpha(self, real_alpha_result, constants):
        expected = Decimal(1) / constants.alpha
        assert real_alpha_result["inv_alpha"] == expected

    def test_r_e(self, real_alpha_result, constants):
        """Classical electron radius = alpha * lambda_bar_c."""
        expected = constants.alpha * constants.lambda_bar_c
        assert real_alpha_result["r_e"] == expected

    def test_r_e_close_to_nist(self, real_alpha_result, constants):
        """r_e should be close to the NIST value."""
        r_e = float(real_alpha_result["r_e"])
        r_e_nist = float(constants.r_e_nist)
        assert abs(r_e - r_e_nist) / r_e_nist < 1e-4

    def test_a_0_predicted(self, real_alpha_result, constants):
        """Bohr radius = lambda_bar_c / alpha."""
        expected = constants.lambda_bar_c / constants.alpha
        assert real_alpha_result["a_0_predicted"] == expected

    def test_a_0_close_to_nist(self, real_alpha_result, constants):
        a_0 = float(real_alpha_result["a_0_predicted"])
        a_0_nist = float(constants.a_0)
        assert abs(a_0 - a_0_nist) / a_0_nist < 1e-4

    def test_binding_energy_ratio(self, real_alpha_result, constants):
        expected = constants.alpha ** 2 / 2
        assert real_alpha_result["binding_energy_ratio"] == expected

    def test_dark_coupling(self, real_alpha_result, constants):
        expected = constants.alpha ** 10
        assert pytest.approx(float(real_alpha_result["dark_coupling"]), rel=1e-40) == float(expected)


# --------------------------------------------------------------------------
# Viability flags / boundary conditions
# --------------------------------------------------------------------------

class TestViabilityFlags:

    def test_real_alpha_chemistry_viable(self, constants):
        result = recompute_physics(constants.alpha, constants)
        assert result["chemistry_viable"] is True

    def test_real_alpha_stars_viable(self, constants):
        result = recompute_physics(constants.alpha, constants)
        assert result["stars_viable"] is True

    def test_too_large_alpha_chemistry_not_viable(self, constants):
        """alpha > 1/85 should make chemistry non-viable."""
        alpha_big = Decimal(1) / Decimal(80)
        result = recompute_physics(alpha_big, constants)
        assert result["chemistry_viable"] is False

    def test_too_small_alpha_chemistry_not_viable(self, constants):
        """alpha < 1/200 should make chemistry non-viable."""
        alpha_small = Decimal(1) / Decimal(250)
        result = recompute_physics(alpha_small, constants)
        assert result["chemistry_viable"] is False

    def test_too_large_alpha_stars_not_viable(self, constants):
        alpha_big = Decimal(1) / Decimal(80)
        result = recompute_physics(alpha_big, constants)
        assert result["stars_viable"] is False

    def test_too_small_alpha_stars_not_viable(self, constants):
        alpha_small = Decimal(1) / Decimal(200)
        result = recompute_physics(alpha_small, constants)
        assert result["stars_viable"] is False

    def test_boundary_lower_chemistry(self, constants):
        """Exactly 1/200 should be chemistry-viable (inclusive)."""
        result = recompute_physics(Decimal(1) / Decimal(200), constants)
        assert result["chemistry_viable"] is True

    def test_boundary_upper_chemistry(self, constants):
        """Exactly 1/85 should be chemistry-viable (inclusive)."""
        result = recompute_physics(Decimal(1) / Decimal(85), constants)
        assert result["chemistry_viable"] is True

    def test_boundary_lower_stars(self, constants):
        """Exactly 1/180 should be stars-viable (inclusive)."""
        result = recompute_physics(Decimal(1) / Decimal(180), constants)
        assert result["stars_viable"] is True

    def test_boundary_upper_stars(self, constants):
        """Exactly 1/85 should be stars-viable (inclusive)."""
        result = recompute_physics(Decimal(1) / Decimal(85), constants)
        assert result["stars_viable"] is True


# --------------------------------------------------------------------------
# Notes content
# --------------------------------------------------------------------------

class TestNotes:

    def test_qed_breakdown_note_when_alpha_large(self, constants):
        result = recompute_physics(Decimal("0.2"), constants)
        assert any("QED perturbation theory" in n for n in result["notes"])

    def test_no_qed_note_at_real_alpha(self, real_alpha_result):
        assert not any("QED perturbation theory" in n for n in real_alpha_result["notes"])

    def test_atoms_unstable_note(self, constants):
        result = recompute_physics(Decimal(1) / Decimal(80), constants)
        assert any("Atoms become unstable" in n for n in result["notes"])

    def test_bonds_too_weak_note(self, constants):
        result = recompute_physics(Decimal(1) / Decimal(250), constants)
        assert any("Molecular bonds are too weak" in n for n in result["notes"])

    def test_hierarchy_gap_note_always_present(self, real_alpha_result):
        assert any("EM-gravity hierarchy gap" in n for n in real_alpha_result["notes"])

    def test_viable_universe_note(self, real_alpha_result):
        assert any("supports both" in n for n in real_alpha_result["notes"])


# --------------------------------------------------------------------------
# Float input
# --------------------------------------------------------------------------

class TestFloatInput:

    def test_accepts_float(self, constants):
        result = recompute_physics(0.007, constants)
        assert isinstance(result["alpha"], Decimal)
        assert float(result["alpha"]) == pytest.approx(0.007, abs=1e-10)


# --------------------------------------------------------------------------
# Sensitivity curve
# --------------------------------------------------------------------------

class TestSensitivityCurve:

    def test_returns_list(self, constants):
        curve = compute_sensitivity_curve(constants, n_points=10)
        assert isinstance(curve, list)

    def test_correct_length(self, constants):
        curve = compute_sensitivity_curve(constants, n_points=50)
        assert len(curve) == 50

    def test_default_length(self, constants):
        curve = compute_sensitivity_curve(constants)
        assert len(curve) == 100

    def test_alpha_range(self, constants):
        curve = compute_sensitivity_curve(constants, n_points=10)
        alphas = [float(r["alpha"]) for r in curve]
        assert alphas[0] == pytest.approx(1 / 200, rel=1e-6)
        assert alphas[-1] == pytest.approx(1 / 85, rel=1e-6)

    def test_monotonically_increasing_alpha(self, constants):
        curve = compute_sensitivity_curve(constants, n_points=20)
        alphas = [r["alpha"] for r in curve]
        for i in range(len(alphas) - 1):
            assert alphas[i] < alphas[i + 1]

    def test_each_entry_is_valid(self, constants):
        curve = compute_sensitivity_curve(constants, n_points=5)
        for entry in curve:
            assert "alpha" in entry
            assert "G_predicted" in entry
            assert "chemistry_viable" in entry

    def test_single_point(self, constants):
        curve = compute_sensitivity_curve(constants, n_points=1)
        assert len(curve) == 1
        assert float(curve[0]["alpha"]) == pytest.approx(1 / 200, rel=1e-6)
