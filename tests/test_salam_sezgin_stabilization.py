"""
Tests for alpha_ladder_core/salam_sezgin_stabilization.py

Comprehensive test suite for the Salam-Sezgin stabilization of the S^2 radius.
The mechanism adds Lambda_6 > 0 to the 6D action, creating a minimum in the
radial potential V(R) = Lambda_6*R^2 - 2/R^2 + n^2/(2*R^4).

Key physics:
    - Gauge-matched radius R_phys = 0.550293 l_Pl (from alpha_em)
    - Required Lambda_6 ~ 14.2 (O(1), natural, no fine-tuning)
    - Dilaton mass ~ Planck scale (completely decoupled)
    - 4D CC ~ O(1) M_Pl^4 (122 orders from observed -- standard CC problem)
    - G prediction at -0.31 ppm is independent of all this
"""

import math
import pytest

from alpha_ladder_core.salam_sezgin_stabilization import (
    compute_salam_sezgin_potential,
    find_ss_minimum,
    compute_ss_dilaton_mass,
    compute_4d_cosmological_constant,
    scan_lambda6,
    find_lambda6_for_target_radius,
    gauge_matching_constraint,
    summarize_salam_sezgin,
)


# ===========================================================================
# 1. compute_salam_sezgin_potential (~8 tests)
# ===========================================================================


def test_potential_returns_dict_with_expected_keys():
    result = compute_salam_sezgin_potential([0.5, 1.0, 2.0], Lambda_6=10.0)
    assert isinstance(result, dict)
    for key in ("R_grid", "V_cc", "V_curv", "V_flux", "V_total",
                "Lambda_6", "n", "formula"):
        assert key in result, f"Missing key: {key}"


def test_potential_grid_lengths_match():
    R_grid = [0.3, 0.5, 0.7, 1.0, 1.5]
    result = compute_salam_sezgin_potential(R_grid, Lambda_6=10.0)
    for key in ("V_cc", "V_curv", "V_flux", "V_total"):
        assert len(result[key]) == len(R_grid)


def test_potential_components_at_R1():
    result = compute_salam_sezgin_potential([1.0], Lambda_6=10.0, n=1)
    assert abs(result["V_cc"][0] - 10.0) < 1e-12
    assert abs(result["V_curv"][0] - (-2.0)) < 1e-12
    assert abs(result["V_flux"][0] - 0.5) < 1e-12
    assert abs(result["V_total"][0] - 8.5) < 1e-12


def test_potential_flux_scales_with_n_squared():
    r1 = compute_salam_sezgin_potential([1.0], Lambda_6=10.0, n=1)
    r2 = compute_salam_sezgin_potential([1.0], Lambda_6=10.0, n=2)
    assert abs(r2["V_flux"][0] / r1["V_flux"][0] - 4.0) < 1e-12


def test_potential_cc_term_grows_with_R():
    result = compute_salam_sezgin_potential([0.5, 1.0, 2.0], Lambda_6=10.0)
    assert result["V_cc"][0] < result["V_cc"][1] < result["V_cc"][2]


def test_potential_curvature_term_negative():
    result = compute_salam_sezgin_potential([0.3, 1.0, 5.0], Lambda_6=10.0)
    for v in result["V_curv"]:
        assert v < 0


def test_potential_handles_zero_R():
    result = compute_salam_sezgin_potential([0.0, 1.0], Lambda_6=10.0)
    assert not math.isfinite(result["V_cc"][0]) or math.isnan(result["V_total"][0])
    assert math.isfinite(result["V_total"][1])


def test_potential_has_minimum_in_interior():
    R_grid = [0.1 + 0.05 * i for i in range(40)]
    result = compute_salam_sezgin_potential(R_grid, Lambda_6=14.0, n=1)
    V = result["V_total"]
    has_interior_min = any(
        V[i] < V[i - 1] and V[i] < V[i + 1]
        for i in range(1, len(V) - 1)
    )
    assert has_interior_min


# ===========================================================================
# 2. find_ss_minimum (~10 tests)
# ===========================================================================


def test_minimum_returns_dict_with_expected_keys():
    result = find_ss_minimum(Lambda_6=10.0)
    assert isinstance(result, dict)
    for key in ("minimum_exists", "R_0", "R_0_squared", "V_at_minimum",
                "V_double_prime", "is_minimum", "Lambda_6", "n",
                "cubic_equation", "status"):
        assert key in result, f"Missing key: {key}"


def test_minimum_exists_for_positive_lambda6():
    result = find_ss_minimum(Lambda_6=10.0)
    assert result["minimum_exists"] is True
    assert result["status"] == "DERIVED"


def test_minimum_fails_for_negative_lambda6():
    result = find_ss_minimum(Lambda_6=-1.0)
    assert result["minimum_exists"] is False
    assert "FAILED" in result["status"]


def test_minimum_fails_for_zero_lambda6():
    result = find_ss_minimum(Lambda_6=0.0)
    assert result["minimum_exists"] is False


def test_minimum_R0_positive():
    result = find_ss_minimum(Lambda_6=14.0)
    assert result["R_0"] > 0


def test_minimum_stationarity_condition():
    Lambda_6 = 14.0
    result = find_ss_minimum(Lambda_6, n=1)
    R0 = result["R_0"]
    R2 = R0 * R0
    R6 = R2 * R2 * R2
    residual = Lambda_6 * R6 + 2.0 * R2 - 1.0
    assert abs(residual) < 1e-10


def test_minimum_second_derivative_positive():
    result = find_ss_minimum(Lambda_6=14.0)
    assert result["V_double_prime"] > 0


def test_minimum_second_derivative_value():
    gm = gauge_matching_constraint()
    result = find_ss_minimum(gm["Lambda_6"])
    assert abs(result["V_double_prime"] - 257.7) / 257.7 < 0.05


def test_minimum_larger_lambda6_gives_smaller_R0():
    r1 = find_ss_minimum(Lambda_6=5.0)
    r2 = find_ss_minimum(Lambda_6=50.0)
    assert r2["R_0"] < r1["R_0"]


def test_minimum_larger_n_gives_larger_R0():
    r1 = find_ss_minimum(Lambda_6=14.0, n=1)
    r2 = find_ss_minimum(Lambda_6=14.0, n=2)
    assert r2["R_0"] > r1["R_0"]


# ===========================================================================
# 3. compute_ss_dilaton_mass (~7 tests)
# ===========================================================================


def test_dilaton_mass_returns_dict():
    result = compute_ss_dilaton_mass(Lambda_6=14.0)
    assert isinstance(result, dict)
    for key in ("m_phi_squared", "m_phi_planck", "m_phi_eV",
                "lambda_compton_m", "R_0", "mass_scale", "status"):
        assert key in result


def test_dilaton_mass_positive():
    result = compute_ss_dilaton_mass(Lambda_6=14.0)
    assert result["m_phi_eV"] > 0
    assert result["m_phi_planck"] > 0


def test_dilaton_mass_planck_scale():
    gm = gauge_matching_constraint()
    result = compute_ss_dilaton_mass(gm["Lambda_6"])
    assert result["mass_scale"] == "planck_scale"
    assert result["m_phi_eV"] > 1e29
    assert result["m_phi_eV"] < 1e30


def test_dilaton_mass_compton_wavelength_tiny():
    gm = gauge_matching_constraint()
    result = compute_ss_dilaton_mass(gm["Lambda_6"])
    assert result["lambda_compton_m"] < 1e-30


def test_dilaton_mass_fails_for_negative_lambda6():
    result = compute_ss_dilaton_mass(Lambda_6=-1.0)
    assert result["m_phi_eV"] is None
    assert "FAILED" in result["status"]


def test_dilaton_mass_increases_with_lambda6():
    r1 = compute_ss_dilaton_mass(Lambda_6=5.0)
    r2 = compute_ss_dilaton_mass(Lambda_6=50.0)
    assert r2["m_phi_eV"] > r1["m_phi_eV"]


def test_dilaton_mass_squared_equals_V_double_prime():
    Lambda_6 = 14.0
    mass = compute_ss_dilaton_mass(Lambda_6)
    minimum = find_ss_minimum(Lambda_6)
    assert abs(mass["m_phi_squared"] - minimum["V_double_prime"]) < 1e-12


# ===========================================================================
# 4. compute_4d_cosmological_constant (~7 tests)
# ===========================================================================


def test_cc_returns_dict():
    result = compute_4d_cosmological_constant(Lambda_6=14.0)
    assert isinstance(result, dict)
    for key in ("Lambda_4", "Lambda_4_obs", "ratio", "log10_ratio",
                "discrepancy_orders", "sign", "honest_assessment", "status"):
        assert key in result


def test_cc_positive_for_gauge_matched():
    gm = gauge_matching_constraint()
    result = compute_4d_cosmological_constant(gm["Lambda_6"])
    assert result["Lambda_4"] > 0
    assert result["sign"] == "positive"


def test_cc_order_of_magnitude():
    gm = gauge_matching_constraint()
    result = compute_4d_cosmological_constant(gm["Lambda_6"])
    assert 1.0 < result["Lambda_4"] < 10.0


def test_cc_discrepancy_122_orders():
    gm = gauge_matching_constraint()
    result = compute_4d_cosmological_constant(gm["Lambda_6"])
    assert result["discrepancy_orders"] is not None
    assert abs(result["discrepancy_orders"] - 122) < 2


def test_cc_fails_for_negative_lambda6():
    result = compute_4d_cosmological_constant(Lambda_6=-1.0)
    assert result["Lambda_4"] is None
    assert result["status"] == "FAILED"


def test_cc_observed_value_stored():
    result = compute_4d_cosmological_constant(Lambda_6=14.0)
    assert abs(result["Lambda_4_obs"] - 2.888e-122) / 2.888e-122 < 0.01


def test_cc_honest_assessment_mentions_standard_problem():
    gm = gauge_matching_constraint()
    result = compute_4d_cosmological_constant(gm["Lambda_6"])
    honest = result["honest_assessment"].lower()
    assert "122" in result["honest_assessment"] or "cosmological" in honest


# ===========================================================================
# 5. scan_lambda6 (~5 tests)
# ===========================================================================


def test_scan_returns_dict():
    result = scan_lambda6([1.0, 10.0, 100.0])
    assert isinstance(result, dict)
    for key in ("results", "n", "n_scanned", "all_stable",
                "R_0_range", "mass_range_eV"):
        assert key in result


def test_scan_correct_count():
    values = [1.0, 5.0, 10.0, 50.0, 100.0]
    result = scan_lambda6(values)
    assert result["n_scanned"] == len(values)
    assert len(result["results"]) == len(values)


def test_scan_all_stable_for_positive_lambda6():
    result = scan_lambda6([1.0, 10.0, 100.0])
    assert result["all_stable"] is True


def test_scan_R0_range_monotonic():
    values = [1.0, 5.0, 10.0, 50.0, 100.0]
    result = scan_lambda6(values)
    radii = [r["R_0"] for r in result["results"]]
    for i in range(len(radii) - 1):
        assert radii[i] > radii[i + 1]


def test_scan_mass_range_exists():
    result = scan_lambda6([1.0, 10.0, 100.0])
    assert result["mass_range_eV"] is not None
    assert result["mass_range_eV"][0] > 0


# ===========================================================================
# 6. find_lambda6_for_target_radius (~8 tests)
# ===========================================================================


def test_find_lambda6_returns_dict():
    result = find_lambda6_for_target_radius(R_target=0.55)
    assert isinstance(result, dict)
    for key in ("Lambda_6", "R_target", "n", "R_max", "is_natural",
                "naturalness_ratio", "V_at_minimum", "m_phi_eV",
                "formula", "status"):
        assert key in result


def test_find_lambda6_gauge_matched_value():
    result = find_lambda6_for_target_radius(R_target=0.550293)
    assert result["Lambda_6"] is not None
    assert abs(result["Lambda_6"] - 14.2) / 14.2 < 0.02


def test_find_lambda6_is_natural():
    result = find_lambda6_for_target_radius(R_target=0.550293)
    assert result["is_natural"] is True


def test_find_lambda6_inverts_minimum():
    R_target = 0.5
    result = find_lambda6_for_target_radius(R_target)
    minimum = find_ss_minimum(result["Lambda_6"])
    assert abs(minimum["R_0"] - R_target) / R_target < 1e-8


def test_find_lambda6_R_max():
    result = find_lambda6_for_target_radius(R_target=0.5)
    assert abs(result["R_max"] - 1.0 / math.sqrt(2.0)) < 1e-10


def test_find_lambda6_fails_for_R_above_Rmax():
    R_max = 1.0 / math.sqrt(2.0)
    result = find_lambda6_for_target_radius(R_target=R_max + 0.1)
    assert result["Lambda_6"] <= 0
    assert "FAILED" in result["status"]


def test_find_lambda6_fails_for_R_zero():
    result = find_lambda6_for_target_radius(R_target=0.0)
    assert "FAILED" in result["status"]


def test_find_lambda6_algebraic_formula():
    R = 0.4
    result = find_lambda6_for_target_radius(R, n=1)
    expected = (1.0 - 2.0 * R * R) / (R ** 6)
    assert abs(result["Lambda_6"] - expected) < 1e-12


# ===========================================================================
# 7. gauge_matching_constraint (~7 tests)
# ===========================================================================


def test_gauge_matching_returns_dict():
    result = gauge_matching_constraint()
    assert isinstance(result, dict)
    for key in ("R_phys", "alpha_em", "phi_vev", "g_KK", "Lambda_6",
                "is_natural", "m_phi_eV", "Lambda_4",
                "honest_assessment", "status"):
        assert key in result


def test_gauge_matching_default_alpha():
    result = gauge_matching_constraint()
    assert abs(result["alpha_em"] - 1.0 / 137.036) < 1e-10


def test_gauge_matching_R_phys():
    result = gauge_matching_constraint()
    assert abs(result["R_phys"] - 0.550293) < 0.001


def test_gauge_matching_phi_vev():
    result = gauge_matching_constraint()
    expected = 0.25 * math.log(4.0 * math.pi / 137.036)
    assert abs(result["phi_vev"] - expected) < 1e-10


def test_gauge_matching_lambda6_natural():
    result = gauge_matching_constraint()
    assert result["is_natural"] is True
    assert 1.0 < result["Lambda_6"] < 100.0


def test_gauge_matching_custom_alpha():
    r1 = gauge_matching_constraint(alpha_em=1.0 / 137.036)
    r2 = gauge_matching_constraint(alpha_em=1.0 / 100.0)
    assert r2["R_phys"] > r1["R_phys"]


def test_gauge_matching_R_below_Rmax():
    result = gauge_matching_constraint()
    R_max = 1.0 / math.sqrt(2.0)
    assert result["R_phys"] < R_max


# ===========================================================================
# 8. summarize_salam_sezgin (~8 tests)
# ===========================================================================


def test_summary_returns_dict():
    result = summarize_salam_sezgin()
    assert isinstance(result, dict)
    for key in ("gauge_match", "scan", "potential", "dilaton_mass",
                "cc_analysis", "gap1_status", "radius_fixed",
                "overall_assessment", "honest_assessment", "first_principles"):
        assert key in result


def test_summary_radius_fixed():
    result = summarize_salam_sezgin()
    assert result["radius_fixed"] is True


def test_summary_gap1_resolved():
    result = summarize_salam_sezgin()
    assert "RESOLVED" in result["gap1_status"]


def test_summary_potential_computed():
    result = summarize_salam_sezgin()
    pot = result["potential"]
    assert pot is not None
    assert len(pot["R_grid"]) > 0


def test_summary_dilaton_mass_computed():
    result = summarize_salam_sezgin()
    dm = result["dilaton_mass"]
    assert dm is not None
    assert dm["m_phi_eV"] > 0


def test_summary_cc_analysis_computed():
    result = summarize_salam_sezgin()
    cc = result["cc_analysis"]
    assert cc is not None
    assert cc["Lambda_4"] is not None


def test_summary_honest_assessment_mentions_trade():
    result = summarize_salam_sezgin()
    honest = result["honest_assessment"].lower()
    assert "lambda" in honest or "parameter" in honest


def test_summary_scan_has_50_points():
    result = summarize_salam_sezgin()
    assert result["scan"]["n_scanned"] == 50


# ===========================================================================
# Cross-function consistency (~5 tests)
# ===========================================================================


def test_gauge_match_lambda6_reproduces_R_phys():
    gm = gauge_matching_constraint()
    minimum = find_ss_minimum(gm["Lambda_6"])
    assert abs(minimum["R_0"] - gm["R_phys"]) / gm["R_phys"] < 1e-6


def test_potential_minimum_matches_find_ss_minimum():
    Lambda_6 = 14.0
    minimum = find_ss_minimum(Lambda_6)
    pot = compute_salam_sezgin_potential([minimum["R_0"]], Lambda_6)
    assert abs(pot["V_total"][0] - minimum["V_at_minimum"]) < 1e-10


def test_dilaton_mass_consistent_with_minimum():
    Lambda_6 = 14.0
    mass = compute_ss_dilaton_mass(Lambda_6)
    minimum = find_ss_minimum(Lambda_6)
    assert abs(mass["R_0"] - minimum["R_0"]) < 1e-12


def test_V_prime_zero_at_minimum_direct():
    Lambda_6 = 20.0
    result = find_ss_minimum(Lambda_6, n=1)
    R = result["R_0"]
    V_prime = 2.0 * Lambda_6 * R + 4.0 / (R ** 3) - 2.0 / (R ** 5)
    assert abs(V_prime) < 1e-8


def test_roundtrip_radius_to_lambda6_and_back():
    for R_target in [0.3, 0.4, 0.5, 0.6]:
        result = find_lambda6_for_target_radius(R_target)
        if result["Lambda_6"] is not None and result["Lambda_6"] > 0:
            minimum = find_ss_minimum(result["Lambda_6"])
            assert abs(minimum["R_0"] - R_target) / R_target < 1e-6, (
                f"Roundtrip failed for R={R_target}: got R_0={minimum['R_0']}"
            )
