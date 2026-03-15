"""Tests for alpha_ladder_core.experimental module."""

import math
import sys
import os
import pytest
from decimal import Decimal

# Ensure the project root is on the path so imports work without install.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.experimental import (
    strategy_mass_ratios,
    strategy_multiple_paths,
    strategy_dark_sector,
    strategy_muon_g2,
    strategy_experimental_approaches,
    get_all_strategies,
)


@pytest.fixture
def constants():
    return get_constants("CODATA 2018")


# ------------------------------------------------------------------
# strategy_mass_ratios
# ------------------------------------------------------------------

class TestStrategyMassRatios:
    def test_returns_list_of_three(self, constants):
        result = strategy_mass_ratios(constants)
        assert isinstance(result, list)
        assert len(result) == 3

    def test_each_entry_has_required_keys(self, constants):
        result = strategy_mass_ratios(constants)
        required = {"name", "ratio", "n_exact", "n_round", "coeff_at_round", "nearby_rungs"}
        for entry in result:
            assert required.issubset(entry.keys()), f"Missing keys in {entry['name']}"

    def test_names_correct(self, constants):
        result = strategy_mass_ratios(constants)
        names = [r["name"] for r in result]
        assert names == ["m_p / m_e", "m_mu / m_e", "m_p / m_mu"]

    def test_mp_me_ratio_value(self, constants):
        result = strategy_mass_ratios(constants)
        mp_me = result[0]
        assert mp_me["name"] == "m_p / m_e"
        assert pytest.approx(mp_me["ratio"], rel=1e-3) == 1836.15

    def test_mmu_me_ratio_value(self, constants):
        result = strategy_mass_ratios(constants)
        mmu_me = result[1]
        assert mmu_me["name"] == "m_mu / m_e"
        assert pytest.approx(mmu_me["ratio"], rel=1e-3) == 206.77

    def test_mp_mmu_ratio_value(self, constants):
        result = strategy_mass_ratios(constants)
        mp_mmu = result[2]
        assert mp_mmu["name"] == "m_p / m_mu"
        assert pytest.approx(mp_mmu["ratio"], rel=1e-2) == 8.88

    def test_all_ratios_positive(self, constants):
        result = strategy_mass_ratios(constants)
        for entry in result:
            assert entry["ratio"] > 0, f"Ratio for {entry['name']} should be positive"

    def test_all_n_exact_positive(self, constants):
        result = strategy_mass_ratios(constants)
        for entry in result:
            assert entry["n_exact"] > 0, f"n_exact for {entry['name']} should be positive"

    def test_mp_me_n_exact_approximate(self, constants):
        result = strategy_mass_ratios(constants)
        mp_me = result[0]
        assert pytest.approx(mp_me["n_exact"], abs=0.1) == 1.53

    def test_mmu_me_n_exact_approximate(self, constants):
        result = strategy_mass_ratios(constants)
        mmu_me = result[1]
        assert pytest.approx(mmu_me["n_exact"], abs=0.1) == 1.08

    def test_mp_mmu_n_exact_approximate(self, constants):
        result = strategy_mass_ratios(constants)
        mp_mmu = result[2]
        assert pytest.approx(mp_mmu["n_exact"], abs=0.1) == 0.44

    def test_n_round_is_integer(self, constants):
        result = strategy_mass_ratios(constants)
        for entry in result:
            assert isinstance(entry["n_round"], int), f"n_round for {entry['name']} should be int"

    def test_n_round_is_nearest_integer_to_n_exact(self, constants):
        result = strategy_mass_ratios(constants)
        for entry in result:
            assert entry["n_round"] == round(entry["n_exact"])

    def test_coeff_at_round_is_float(self, constants):
        result = strategy_mass_ratios(constants)
        for entry in result:
            assert isinstance(entry["coeff_at_round"], float)

    def test_nearby_rungs_has_three_entries(self, constants):
        result = strategy_mass_ratios(constants)
        for entry in result:
            assert len(entry["nearby_rungs"]) == 3

    def test_nearby_rungs_have_n_and_coeff(self, constants):
        result = strategy_mass_ratios(constants)
        for entry in result:
            for rung in entry["nearby_rungs"]:
                assert "n" in rung
                assert "coeff" in rung

    def test_nearby_rungs_n_values_consecutive(self, constants):
        """The three nearby rungs should cover n_round-1, n_round, n_round+1 (negated)."""
        result = strategy_mass_ratios(constants)
        for entry in result:
            n_vals = [r["n"] for r in entry["nearby_rungs"]]
            n_round = entry["n_round"]
            expected = [-(n_round - 1), -n_round, -(n_round + 1)]
            assert n_vals == expected

    def test_nearby_rungs_coeff_positive(self, constants):
        result = strategy_mass_ratios(constants)
        for entry in result:
            for rung in entry["nearby_rungs"]:
                assert rung["coeff"] > 0

    def test_ratio_consistency_mp_me_times_mmu_me(self, constants):
        """m_p/m_e should equal (m_p/m_mu) * (m_mu/m_e)."""
        result = strategy_mass_ratios(constants)
        mp_me = result[0]["ratio"]
        mmu_me = result[1]["ratio"]
        mp_mmu = result[2]["ratio"]
        assert pytest.approx(mp_me, rel=1e-10) == mp_mmu * mmu_me


# ------------------------------------------------------------------
# strategy_multiple_paths
# ------------------------------------------------------------------

class TestStrategyMultiplePaths:
    def test_returns_dict(self, constants):
        result = strategy_multiple_paths(constants)
        assert isinstance(result, dict)

    def test_has_all_paths(self, constants):
        result = strategy_multiple_paths(constants)
        assert {"path_A", "path_B", "path_C", "path_D"}.issubset(result.keys())

    # -- Path A --
    def test_path_a_has_description(self, constants):
        pa = strategy_multiple_paths(constants)["path_A"]
        assert "description" in pa

    def test_path_a_has_alpha_g(self, constants):
        pa = strategy_multiple_paths(constants)["path_A"]
        assert "alpha_G" in pa
        assert isinstance(pa["alpha_G"], float)
        assert pa["alpha_G"] > 0

    def test_path_a_g_predicted(self, constants):
        pa = strategy_multiple_paths(constants)["path_A"]
        assert "G_predicted" in pa
        assert pytest.approx(pa["G_predicted"], rel=1e-3) == 6.673e-11

    def test_path_a_g_predicted_order_of_magnitude(self, constants):
        pa = strategy_multiple_paths(constants)["path_A"]
        assert 1e-12 < pa["G_predicted"] < 1e-9

    # -- Path B --
    def test_path_b_has_required_keys(self, constants):
        pb = strategy_multiple_paths(constants)["path_B"]
        for key in ["alpha_G_proton", "n_proton", "nearest_rung", "coeff_at_rung"]:
            assert key in pb, f"Missing key {key} in path_B"

    def test_path_b_alpha_g_proton_positive(self, constants):
        pb = strategy_multiple_paths(constants)["path_B"]
        assert pb["alpha_G_proton"] > 0

    def test_path_b_n_proton_finite(self, constants):
        pb = strategy_multiple_paths(constants)["path_B"]
        assert math.isfinite(pb["n_proton"])

    def test_path_b_nearest_rung_is_integer(self, constants):
        pb = strategy_multiple_paths(constants)["path_B"]
        assert isinstance(pb["nearest_rung"], int)

    def test_path_b_coeff_at_rung_positive(self, constants):
        pb = strategy_multiple_paths(constants)["path_B"]
        assert pb["coeff_at_rung"] > 0

    # -- Path C --
    def test_path_c_has_required_keys(self, constants):
        pc = strategy_multiple_paths(constants)["path_C"]
        for key in ["m_e_over_m_Pl", "alpha_G_from_planck", "n_planck", "n_alpha_G"]:
            assert key in pc, f"Missing key {key} in path_C"

    def test_path_c_m_e_over_m_pl_very_small(self, constants):
        pc = strategy_multiple_paths(constants)["path_C"]
        assert pc["m_e_over_m_Pl"] < 1e-20

    def test_path_c_n_alpha_g_is_twice_n_planck(self, constants):
        pc = strategy_multiple_paths(constants)["path_C"]
        assert pytest.approx(pc["n_alpha_G"], rel=1e-12) == 2 * pc["n_planck"]

    def test_path_c_alpha_g_from_planck_positive(self, constants):
        pc = strategy_multiple_paths(constants)["path_C"]
        assert pc["alpha_G_from_planck"] > 0

    def test_path_c_n_planck_finite(self, constants):
        pc = strategy_multiple_paths(constants)["path_C"]
        assert math.isfinite(pc["n_planck"])

    # -- Path D --
    def test_path_d_has_required_keys(self, constants):
        pd = strategy_multiple_paths(constants)["path_D"]
        for key in ["alpha_G", "G_predicted", "residual_ppm"]:
            assert key in pd, f"Missing key {key} in path_D"

    def test_path_d_g_predicted(self, constants):
        pd = strategy_multiple_paths(constants)["path_D"]
        assert pytest.approx(pd["G_predicted"], rel=1e-2) == 6.679e-11

    def test_path_d_g_predicted_order_of_magnitude(self, constants):
        pd = strategy_multiple_paths(constants)["path_D"]
        assert 1e-12 < pd["G_predicted"] < 1e-9

    def test_path_d_residual_ppm_approximate(self, constants):
        pd = strategy_multiple_paths(constants)["path_D"]
        assert pytest.approx(pd["residual_ppm"], rel=0.1) == 688

    def test_path_d_residual_ppm_positive(self, constants):
        pd = strategy_multiple_paths(constants)["path_D"]
        assert pd["residual_ppm"] > 0

    def test_path_d_description_mentions_zero_free(self, constants):
        pd = strategy_multiple_paths(constants)["path_D"]
        assert "zero free" in pd["description"].lower()

    # -- Cross-path checks --
    def test_both_g_predictions_same_order_of_magnitude(self, constants):
        result = strategy_multiple_paths(constants)
        g_a = result["path_A"]["G_predicted"]
        g_d = result["path_D"]["G_predicted"]
        assert abs(math.log10(g_a) - math.log10(g_d)) < 1


# ------------------------------------------------------------------
# strategy_dark_sector
# ------------------------------------------------------------------

class TestStrategyDarkSector:
    def test_returns_dict(self, constants):
        result = strategy_dark_sector(constants)
        assert isinstance(result, dict)

    def test_has_required_keys(self, constants):
        result = strategy_dark_sector(constants)
        required = {
            "alpha_10", "epsilon_predicted", "experimental_bounds",
            "below_current_bounds", "orders_below_sensitivity",
        }
        assert required.issubset(result.keys())

    def test_alpha_10_value(self, constants):
        result = strategy_dark_sector(constants)
        assert pytest.approx(result["alpha_10"], rel=0.1) == 4.4e-22

    def test_alpha_10_is_alpha_to_the_tenth(self, constants):
        alpha = float(constants.alpha)
        result = strategy_dark_sector(constants)
        assert pytest.approx(result["alpha_10"], rel=1e-8) == alpha ** 10

    def test_epsilon_predicted_value(self, constants):
        result = strategy_dark_sector(constants)
        assert pytest.approx(result["epsilon_predicted"], rel=0.1) == 2.1e-11

    def test_epsilon_is_sqrt_of_alpha_10(self, constants):
        result = strategy_dark_sector(constants)
        assert pytest.approx(result["epsilon_predicted"], rel=1e-8) == math.sqrt(result["alpha_10"])

    def test_epsilon_equals_alpha_to_the_fifth(self, constants):
        alpha = float(constants.alpha)
        result = strategy_dark_sector(constants)
        assert pytest.approx(result["epsilon_predicted"], rel=1e-8) == alpha ** 5

    def test_below_current_bounds_true(self, constants):
        result = strategy_dark_sector(constants)
        assert result["below_current_bounds"] is True

    def test_orders_below_sensitivity_positive(self, constants):
        result = strategy_dark_sector(constants)
        assert result["orders_below_sensitivity"] > 0

    def test_orders_below_sensitivity_greater_than_7(self, constants):
        result = strategy_dark_sector(constants)
        assert result["orders_below_sensitivity"] > 7

    def test_orders_below_sensitivity_value(self, constants):
        """epsilon ~ 2e-11 vs bound 1e-3 => about 7.7 orders."""
        result = strategy_dark_sector(constants)
        expected = math.log10(1e-3 / result["epsilon_predicted"])
        assert pytest.approx(result["orders_below_sensitivity"], rel=1e-6) == expected

    def test_experimental_bounds_is_dict(self, constants):
        result = strategy_dark_sector(constants)
        assert isinstance(result["experimental_bounds"], dict)

    def test_experimental_bounds_has_entries(self, constants):
        result = strategy_dark_sector(constants)
        assert len(result["experimental_bounds"]) >= 3

    def test_experimental_bounds_entries_have_epsilon_max(self, constants):
        result = strategy_dark_sector(constants)
        for name, bound in result["experimental_bounds"].items():
            assert "epsilon_max" in bound, f"Missing epsilon_max in {name}"

    def test_epsilon_below_all_bounds(self, constants):
        result = strategy_dark_sector(constants)
        epsilon = result["epsilon_predicted"]
        for name, bound in result["experimental_bounds"].items():
            assert epsilon < bound["epsilon_max"], f"epsilon not below {name}"

    def test_alpha_10_positive(self, constants):
        result = strategy_dark_sector(constants)
        assert result["alpha_10"] > 0

    def test_epsilon_predicted_positive(self, constants):
        result = strategy_dark_sector(constants)
        assert result["epsilon_predicted"] > 0


# ------------------------------------------------------------------
# strategy_muon_g2
# ------------------------------------------------------------------

class TestStrategyMuonG2:
    def test_returns_dict(self, constants):
        result = strategy_muon_g2(constants)
        assert isinstance(result, dict)

    def test_has_required_keys(self, constants):
        result = strategy_muon_g2(constants)
        required = {
            "delta_a_mu", "a_mu_exp", "fractional_anomaly",
            "n_anomaly", "nearest_rung_anomaly", "anomaly_coeff",
            "n_fractional", "nearest_rung_fractional",
        }
        assert required.issubset(result.keys())

    def test_delta_a_mu_value(self, constants):
        result = strategy_muon_g2(constants)
        assert pytest.approx(result["delta_a_mu"], rel=1e-6) == 2.51e-9

    def test_a_mu_exp_value(self, constants):
        result = strategy_muon_g2(constants)
        assert pytest.approx(result["a_mu_exp"], rel=1e-3) == 1.166e-3

    def test_a_mu_exp_approximately_one_over_860(self, constants):
        result = strategy_muon_g2(constants)
        assert pytest.approx(result["a_mu_exp"], rel=0.01) == 1.0 / 860

    def test_fractional_anomaly_positive(self, constants):
        result = strategy_muon_g2(constants)
        assert result["fractional_anomaly"] > 0

    def test_fractional_anomaly_is_ratio(self, constants):
        result = strategy_muon_g2(constants)
        expected = result["delta_a_mu"] / result["a_mu_exp"]
        assert pytest.approx(result["fractional_anomaly"], rel=1e-6) == expected

    def test_n_anomaly_finite(self, constants):
        result = strategy_muon_g2(constants)
        assert math.isfinite(result["n_anomaly"])

    def test_nearest_rung_anomaly_is_integer(self, constants):
        result = strategy_muon_g2(constants)
        assert isinstance(result["nearest_rung_anomaly"], int)

    def test_nearest_rung_anomaly_close_to_n_anomaly(self, constants):
        result = strategy_muon_g2(constants)
        assert abs(result["nearest_rung_anomaly"] - result["n_anomaly"]) <= 0.5

    def test_anomaly_coeff_is_float(self, constants):
        result = strategy_muon_g2(constants)
        assert isinstance(result["anomaly_coeff"], float)

    def test_n_fractional_finite(self, constants):
        result = strategy_muon_g2(constants)
        assert math.isfinite(result["n_fractional"])

    def test_nearest_rung_fractional_is_integer(self, constants):
        result = strategy_muon_g2(constants)
        assert isinstance(result["nearest_rung_fractional"], int)

    def test_nearest_rung_fractional_close_to_n_fractional(self, constants):
        result = strategy_muon_g2(constants)
        assert abs(result["nearest_rung_fractional"] - result["n_fractional"]) <= 0.5

    def test_delta_a_mu_much_smaller_than_a_mu_exp(self, constants):
        result = strategy_muon_g2(constants)
        assert result["delta_a_mu"] < result["a_mu_exp"] * 0.01


# ------------------------------------------------------------------
# strategy_experimental_approaches
# ------------------------------------------------------------------

class TestStrategyExperimentalApproaches:
    def test_returns_dict(self, constants):
        result = strategy_experimental_approaches(constants)
        assert isinstance(result, dict)

    def test_has_required_keys(self, constants):
        result = strategy_experimental_approaches(constants)
        required = {
            "ladder_prediction", "hierarchy_prediction", "codata_2018",
            "difference", "difference_ppm", "hierarchy_difference_ppm",
            "required_precision_ppm", "approaches",
        }
        assert required.issubset(result.keys())

    def test_ladder_prediction_value(self, constants):
        result = strategy_experimental_approaches(constants)
        assert pytest.approx(result["ladder_prediction"], rel=1e-3) == 6.673e-11

    def test_ladder_prediction_order_of_magnitude(self, constants):
        result = strategy_experimental_approaches(constants)
        assert 1e-12 < result["ladder_prediction"] < 1e-9

    def test_codata_2018_is_dict_with_value(self, constants):
        result = strategy_experimental_approaches(constants)
        codata = result["codata_2018"]
        assert isinstance(codata, dict)
        assert "value" in codata
        assert "uncertainty" in codata

    def test_codata_2018_value(self, constants):
        result = strategy_experimental_approaches(constants)
        assert pytest.approx(result["codata_2018"]["value"], rel=1e-3) == 6.674e-11

    def test_codata_2018_uncertainty_positive(self, constants):
        result = strategy_experimental_approaches(constants)
        assert result["codata_2018"]["uncertainty"] > 0

    def test_hierarchy_prediction_value(self, constants):
        result = strategy_experimental_approaches(constants)
        assert pytest.approx(result["hierarchy_prediction"], rel=1e-2) == 6.679e-11

    def test_hierarchy_prediction_order_of_magnitude(self, constants):
        result = strategy_experimental_approaches(constants)
        assert 1e-12 < result["hierarchy_prediction"] < 1e-9

    def test_hierarchy_difference_ppm_approximate(self, constants):
        result = strategy_experimental_approaches(constants)
        assert pytest.approx(result["hierarchy_difference_ppm"], rel=0.1) == 688

    def test_difference_ppm_value(self, constants):
        result = strategy_experimental_approaches(constants)
        assert result["difference_ppm"] == 160

    def test_required_precision_ppm_value(self, constants):
        result = strategy_experimental_approaches(constants)
        assert result["required_precision_ppm"] == 5

    def test_required_precision_less_than_difference(self, constants):
        result = strategy_experimental_approaches(constants)
        assert result["required_precision_ppm"] < result["difference_ppm"]

    def test_approaches_has_three_entries(self, constants):
        result = strategy_experimental_approaches(constants)
        assert len(result["approaches"]) == 3

    def test_approaches_is_dict(self, constants):
        result = strategy_experimental_approaches(constants)
        assert isinstance(result["approaches"], dict)

    def test_approaches_all_have_description(self, constants):
        result = strategy_experimental_approaches(constants)
        for name, approach in result["approaches"].items():
            assert "description" in approach, f"Missing description in {name}"

    def test_approaches_all_have_technique(self, constants):
        result = strategy_experimental_approaches(constants)
        for name, approach in result["approaches"].items():
            assert "technique" in approach or "note" in approach, \
                f"Missing technique/note in {name}"

    def test_atom_interferometry_present(self, constants):
        result = strategy_experimental_approaches(constants)
        assert "A_atom_interferometry" in result["approaches"]

    def test_atom_interferometry_has_groups(self, constants):
        result = strategy_experimental_approaches(constants)
        ai = result["approaches"]["A_atom_interferometry"]
        assert "groups" in ai
        assert len(ai["groups"]) >= 2

    def test_difference_consistent_with_values(self, constants):
        result = strategy_experimental_approaches(constants)
        diff = abs(result["ladder_prediction"] - result["codata_2018"]["value"])
        assert pytest.approx(result["difference"], rel=0.1) == diff


# ------------------------------------------------------------------
# get_all_strategies
# ------------------------------------------------------------------

class TestGetAllStrategies:
    def test_returns_list(self, constants):
        result = get_all_strategies(constants)
        assert isinstance(result, list)

    def test_returns_five_strategies(self, constants):
        result = get_all_strategies(constants)
        assert len(result) == 5

    def test_each_entry_has_strategy_and_data(self, constants):
        result = get_all_strategies(constants)
        for entry in result:
            assert "strategy" in entry
            assert "data" in entry

    def test_strategy_names(self, constants):
        result = get_all_strategies(constants)
        names = [s["strategy"] for s in result]
        assert names == [
            "mass_ratios",
            "multiple_paths",
            "dark_sector",
            "muon_g2",
            "experimental_approaches",
        ]

    def test_mass_ratios_data_matches_direct_call(self, constants):
        all_strats = get_all_strategies(constants)
        direct = strategy_mass_ratios(constants)
        assert all_strats[0]["data"] == direct

    def test_multiple_paths_data_matches_direct_call(self, constants):
        all_strats = get_all_strategies(constants)
        direct = strategy_multiple_paths(constants)
        assert all_strats[1]["data"] == direct

    def test_dark_sector_data_matches_direct_call(self, constants):
        all_strats = get_all_strategies(constants)
        direct = strategy_dark_sector(constants)
        assert all_strats[2]["data"] == direct

    def test_muon_g2_data_matches_direct_call(self, constants):
        all_strats = get_all_strategies(constants)
        direct = strategy_muon_g2(constants)
        assert all_strats[3]["data"] == direct

    def test_experimental_approaches_data_matches_direct_call(self, constants):
        all_strats = get_all_strategies(constants)
        direct = strategy_experimental_approaches(constants)
        assert all_strats[4]["data"] == direct

    def test_all_data_non_empty(self, constants):
        result = get_all_strategies(constants)
        for entry in result:
            assert entry["data"], f"Empty data for {entry['strategy']}"
