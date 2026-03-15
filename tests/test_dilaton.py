"""
Tests for alpha_ladder_core/dilaton.py
"""

import pytest
from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.dilaton import decompose_21, compute_bd_parameter, reconcile_6d


@pytest.fixture
def constants():
    return get_constants()


# -----------------------------------------------------------------------
# Tests for decompose_21
# -----------------------------------------------------------------------

class TestDecompose21:
    def test_returns_required_keys(self, constants):
        result = decompose_21(constants)
        assert "alpha_20" in result
        assert "alpha_1" in result
        assert "bridge" in result
        assert "alpha_21" in result
        assert "alternative_splits" in result

    def test_alpha_20_times_alpha_1_equals_alpha_21(self, constants):
        result = decompose_21(constants)
        product = result["alpha_20"] * result["alpha_1"]
        assert pytest.approx(product, rel=1e-10) == result["alpha_21"]

    def test_alpha_1_equals_alpha(self, constants):
        result = decompose_21(constants)
        assert pytest.approx(result["alpha_1"], rel=1e-12) == float(constants.alpha)

    def test_bridge_is_phi_squared_over_2(self, constants):
        import math
        phi = (1 + math.sqrt(5)) / 2
        expected = phi ** 2 / 2
        result = decompose_21(constants)
        assert pytest.approx(result["bridge"], rel=1e-10) == expected

    def test_alternative_splits_all_sum_to_21(self, constants):
        result = decompose_21(constants)
        for split in result["alternative_splits"]:
            assert split["a"] + split["b"] == 21

    def test_alternative_splits_count(self, constants):
        result = decompose_21(constants)
        # Legacy has 8 alternative splits
        assert len(result["alternative_splits"]) == 8

    def test_20_1_split_has_riemann_meaning(self, constants):
        result = decompose_21(constants)
        split_20_1 = [s for s in result["alternative_splits"] if s["a"] == 20 and s["b"] == 1]
        assert len(split_20_1) == 1
        assert "Riemann" in split_20_1[0]["meaning"]

    def test_alpha_values_are_positive(self, constants):
        result = decompose_21(constants)
        assert result["alpha_20"] > 0
        assert result["alpha_1"] > 0
        assert result["alpha_21"] > 0
        for split in result["alternative_splits"]:
            assert split["alpha_a"] > 0
            assert split["alpha_b"] > 0


# -----------------------------------------------------------------------
# Tests for compute_bd_parameter
# -----------------------------------------------------------------------

class TestComputeBDParameter:
    def test_returns_required_keys(self, constants):
        result = compute_bd_parameter(constants)
        assert "omega" in result
        assert "bare_ratio" in result
        assert "phi_squared_over_2" in result
        assert "dilaton_mass_min_kg" in result
        assert "dilaton_mass_min_eV" in result
        assert "dark_scale_eV" in result

    def test_omega_is_small_positive(self, constants):
        """omega ~ 0.118 for phi^2/2 bridge, positive but far below Cassini bound."""
        result = compute_bd_parameter(constants)
        assert result["omega"] is not None
        assert result["omega"] > 0
        assert result["omega"] < 1

    def test_omega_excluded_for_massless(self, constants):
        """Legacy shows |omega| << 40000 (Cassini bound)."""
        result = compute_bd_parameter(constants)
        assert result["omega_excluded_massless"] is True
        assert result["cassini_bound"] == 40000

    def test_bare_ratio_matches_phi_squared_over_2(self, constants):
        import math
        phi = (1 + math.sqrt(5)) / 2
        expected = phi ** 2 / 2
        result = compute_bd_parameter(constants)
        assert pytest.approx(result["bare_ratio"], rel=1e-10) == expected
        assert pytest.approx(result["phi_squared_over_2"], rel=1e-10) == expected

    def test_omega_formula(self, constants):
        """Verify omega = (3*phi^2/2 - 4) / (2 - phi^2) from legacy."""
        import math
        phi = (1 + math.sqrt(5)) / 2
        phi2 = phi ** 2
        expected_omega = (3 * phi2 / 2 - 4) / (2 - phi2)
        result = compute_bd_parameter(constants)
        assert pytest.approx(result["omega"], rel=1e-10) == expected_omega

    def test_dilaton_mass_min_is_positive(self, constants):
        result = compute_bd_parameter(constants)
        assert result["dilaton_mass_min_kg"] > 0
        assert result["dilaton_mass_min_eV"] > 0

    def test_dark_scale_positive(self, constants):
        result = compute_bd_parameter(constants)
        assert result["dark_scale_eV"] > 0


# -----------------------------------------------------------------------
# Tests for reconcile_6d
# -----------------------------------------------------------------------

class TestReconcile6D:
    def test_returns_required_keys(self, constants):
        result = reconcile_6d(constants)
        assert "g_uv_4d" in result
        assert "g_ab_2d" in result
        assert "g_ua_mixed" in result
        assert "total_6d_metric" in result
        assert "riemann_4d" in result
        assert "scalar_dilaton" in result
        assert "riemann_plus_scalar" in result
        assert "views" in result

    def test_6d_metric_sums_to_21(self, constants):
        result = reconcile_6d(constants)
        total = result["g_uv_4d"] + result["g_ab_2d"] + result["g_ua_mixed"]
        assert total == 21
        assert result["total_6d_metric"] == 21

    def test_riemann_plus_scalar_is_21(self, constants):
        result = reconcile_6d(constants)
        assert result["riemann_4d"] + result["scalar_dilaton"] == 21
        assert result["riemann_plus_scalar"] == 21

    def test_4d_metric_components(self, constants):
        """4D symmetric metric has 4*5/2 = 10 components."""
        result = reconcile_6d(constants)
        assert result["g_uv_4d"] == 10

    def test_2d_metric_components(self, constants):
        """2D symmetric metric has 2*3/2 = 3 components."""
        result = reconcile_6d(constants)
        assert result["g_ab_2d"] == 3

    def test_mixed_components(self, constants):
        """Mixed components: 4 * 2 = 8."""
        result = reconcile_6d(constants)
        assert result["g_ua_mixed"] == 8

    def test_riemann_4d_is_20(self, constants):
        result = reconcile_6d(constants)
        assert result["riemann_4d"] == 20

    def test_views_present(self, constants):
        result = reconcile_6d(constants)
        assert "bottom_up" in result["views"]
        assert "top_down" in result["views"]
