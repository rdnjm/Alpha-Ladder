"""
Tests for alpha_ladder_core/charged_matter_loops.py

One-loop corrections from charged matter fields on S^2 with Dirac monopole
background.  The key physics result: the monopole flips the sign of the
spectral zeta function from -17/480 (negative, pure gravity) to +1/10
(positive, charged q=1, n=1).  This is the first calculation in the
framework producing a positive one-loop correction to G.

Key values:
    - zeta_neutral(-1) = -17/480 = -0.035416... (pure gravity)
    - zeta_charged(-1, q=1, n=1) = 1/10 = 0.1 exactly (sign flip!)
    - j_min = |qn|/2 = 1/2 for q=1, n=1 (monopole harmonics)
    - E8xE8: 3456 scalars + 2472 fermions + 496 vectors (overshoots 3*alpha^2)
"""

import math
import pytest

from alpha_ladder_core.charged_matter_loops import (
    compute_monopole_kk_spectrum,
    compute_monopole_spectral_zeta,
    compute_charged_cw_potential,
    compute_charged_oneloop_g_correction,
    scan_anomaly_free_matter_loops,
    summarize_charged_matter_loops,
)


# ===========================================================================
# 1. compute_monopole_kk_spectrum (~6 tests)
# ===========================================================================


def test_spectrum_returns_dict():
    result = compute_monopole_kk_spectrum(charge_q=1, n_monopole=1)
    assert isinstance(result, dict)
    for key in ("charge_q", "n_monopole", "j_min", "levels",
                "n_levels", "spectrum_shift"):
        assert key in result


def test_spectrum_j_min_neutral():
    result = compute_monopole_kk_spectrum(charge_q=0, n_monopole=1)
    assert result["j_min"] == 0.0


def test_spectrum_j_min_charged():
    result = compute_monopole_kk_spectrum(charge_q=1, n_monopole=1)
    assert result["j_min"] == 0.5


def test_spectrum_j_min_charge2():
    result = compute_monopole_kk_spectrum(charge_q=2, n_monopole=1)
    assert result["j_min"] == 1.0


def test_spectrum_eigenvalues_positive():
    result = compute_monopole_kk_spectrum(charge_q=1, n_monopole=1, l_max=20)
    for lev in result["levels"]:
        assert lev["eigenvalue"] >= 0


def test_spectrum_first_eigenvalue_q1_n1():
    result = compute_monopole_kk_spectrum(charge_q=1, n_monopole=1, l_max=5)
    first = result["levels"][0]
    assert first["j"] == 0.5
    # eigenvalue = j(j+1) - (qn/2)^2 = 0.5*1.5 - 0.25 = 0.5
    assert abs(first["eigenvalue"] - 0.5) < 1e-12
    assert first["degeneracy"] == 2


# ===========================================================================
# 2. compute_monopole_spectral_zeta (~10 tests)
# ===========================================================================


def test_zeta_neutral_matches_known():
    """Neutral zeta at s=-1 should match -17/480."""
    result = compute_monopole_spectral_zeta(charge_q=0, n_monopole=1, s=-1)
    expected = -17.0 / 480.0
    assert abs(result["value"] - expected) < 1e-10


def test_zeta_charged_q1_n1_exact():
    """Charged zeta(q=1, n=1, s=-1) = 1/10 exactly."""
    result = compute_monopole_spectral_zeta(charge_q=1, n_monopole=1, s=-1)
    assert abs(result["value"] - 0.1) < 1e-10


def test_zeta_sign_flip():
    """The monopole flips the spectral zeta from negative to positive."""
    z0 = compute_monopole_spectral_zeta(charge_q=0, n_monopole=1, s=-1)
    z1 = compute_monopole_spectral_zeta(charge_q=1, n_monopole=1, s=-1)
    assert z0["value"] < 0
    assert z1["value"] > 0


def test_zeta_charged_q1_n2():
    """Charged zeta(q=1, n=2) should match q=2, n=1 (same |qn|)."""
    z12 = compute_monopole_spectral_zeta(charge_q=1, n_monopole=2, s=-1)
    z21 = compute_monopole_spectral_zeta(charge_q=2, n_monopole=1, s=-1)
    assert abs(z12["value"] - z21["value"]) < 1e-10


def test_zeta_method_hurwitz():
    """s=-1 should use Hurwitz analytic continuation."""
    result = compute_monopole_spectral_zeta(charge_q=1, n_monopole=1, s=-1)
    assert result["method"] == "hurwitz_analytic_continuation"
    assert result["converged"] is True


def test_zeta_converges_for_s_gt_1():
    """For s > 1, direct summation should converge."""
    result = compute_monopole_spectral_zeta(charge_q=1, n_monopole=1, s=2)
    assert result["method"] == "direct"
    assert result["converged"] is True
    assert result["value"] > 0


def test_zeta_comparison_neutral_provided():
    """The comparison_neutral field should be populated for s=-1."""
    result = compute_monopole_spectral_zeta(charge_q=1, n_monopole=1, s=-1)
    assert result["comparison_neutral"] is not None
    assert abs(result["comparison_neutral"] - (-17.0 / 480.0)) < 1e-10


def test_zeta_increases_with_qn():
    """Larger |qn| should give larger zeta (more positive)."""
    z1 = compute_monopole_spectral_zeta(charge_q=1, n_monopole=1, s=-1)
    z2 = compute_monopole_spectral_zeta(charge_q=1, n_monopole=2, s=-1)
    assert z2["value"] > z1["value"]


def test_zeta_q1_n2_value():
    """zeta(q=1, n=2, s=-1) = 141/160 (from Hurwitz)."""
    result = compute_monopole_spectral_zeta(charge_q=1, n_monopole=2, s=-1)
    assert abs(result["value"] - 141.0 / 160.0) < 1e-8


def test_zeta_neutral_comparison_returned():
    """comparison_neutral should be -17/480 for charged fields."""
    result = compute_monopole_spectral_zeta(charge_q=1, n_monopole=1, s=-1)
    assert abs(result["comparison_neutral"] - (-17.0 / 480.0)) < 1e-10


# ===========================================================================
# 3. compute_charged_cw_potential (~5 tests)
# ===========================================================================


def test_cw_returns_dict():
    result = compute_charged_cw_potential(n_scalars=4, n_fermions=2, n_vectors=1)
    assert isinstance(result, dict)
    for key in ("V_scalars", "V_fermions", "V_vectors", "V_total",
                "n_scalars", "n_fermions", "n_vectors",
                "charge_q", "n_monopole", "R", "j_min",
                "effective_lambda6", "assessment"):
        assert key in result


def test_cw_scalars_positive():
    """Scalar CW potential should be positive (bosonic)."""
    result = compute_charged_cw_potential(n_scalars=4, n_fermions=0, n_vectors=0)
    assert result["V_scalars"] > 0


def test_cw_fermions_negative():
    """Fermion CW potential should be negative (fermionic sign)."""
    result = compute_charged_cw_potential(n_scalars=0, n_fermions=4, n_vectors=0)
    assert result["V_fermions"] < 0


def test_cw_default_R_gauge_matched():
    """Default R should be gauge-matched value."""
    result = compute_charged_cw_potential(n_scalars=1, n_fermions=0, n_vectors=0)
    assert abs(result["R"] - 0.550293) < 0.001


def test_cw_effective_lambda6_computed():
    """effective_lambda6 should be V_total / R^2."""
    result = compute_charged_cw_potential(n_scalars=4, n_fermions=2, n_vectors=1)
    expected = result["V_total"] / (result["R"] ** 2)
    assert abs(result["effective_lambda6"] - expected) < abs(expected) * 1e-10


# ===========================================================================
# 4. compute_charged_oneloop_g_correction (~7 tests)
# ===========================================================================


def test_g_correction_returns_dict():
    result = compute_charged_oneloop_g_correction(
        n_scalars=4, n_fermions=2, n_vectors=1
    )
    assert isinstance(result, dict)
    for key in ("zeta_neutral", "zeta_charged", "delta_zeta",
                "delta_G_over_G_gravity", "delta_G_over_G_matter",
                "delta_G_over_G_total", "target_3_alpha_sq",
                "ratio_to_target", "assessment"):
        assert key in result


def test_g_correction_gravity_negative():
    """Pure gravity contribution should be negative."""
    result = compute_charged_oneloop_g_correction(
        n_scalars=0, n_fermions=0, n_vectors=0
    )
    assert result["delta_G_over_G_gravity"] < 0


def test_g_correction_zeta_neutral_value():
    """zeta_neutral should be -17/480."""
    result = compute_charged_oneloop_g_correction(
        n_scalars=1, n_fermions=0, n_vectors=0
    )
    assert abs(result["zeta_neutral"] - (-17.0 / 480.0)) < 1e-10


def test_g_correction_zeta_charged_positive():
    """Charged zeta should be positive for q=1, n=1."""
    result = compute_charged_oneloop_g_correction(
        n_scalars=1, n_fermions=0, n_vectors=0,
        charge_q=1, n_monopole=1
    )
    assert result["zeta_charged"] > 0


def test_g_correction_matter_contribution_sign():
    """With positive zeta and positive c_matter, matter should give positive delta_G."""
    result = compute_charged_oneloop_g_correction(
        n_scalars=4, n_fermions=0, n_vectors=1,
        charge_q=1, n_monopole=1
    )
    assert result["delta_G_over_G_matter"] > 0


def test_g_correction_e8_overshoots():
    """E8xE8 matter content should massively overshoot 3*alpha^2."""
    result = compute_charged_oneloop_g_correction(
        n_scalars=4 * 740 + 496,
        n_fermions=2 * 740 + 2 * 496,
        n_vectors=496,
        charge_q=1, n_monopole=1
    )
    assert result["ratio_to_target"] > 1000


def test_g_correction_right_order_of_magnitude():
    """With small matter content, correction should be O(10^{-4}) to O(10^{-3})."""
    result = compute_charged_oneloop_g_correction(
        n_scalars=4, n_fermions=2, n_vectors=1,
        charge_q=1, n_monopole=1
    )
    total = abs(result["delta_G_over_G_total"])
    assert 1e-5 < total < 1e-2


# ===========================================================================
# 5. scan_anomaly_free_matter_loops (~4 tests)
# ===========================================================================


def test_scan_returns_dict():
    result = scan_anomaly_free_matter_loops()
    assert isinstance(result, dict)
    for key in ("results", "n_groups", "any_match_3alpha2",
                "any_significant_lambda6", "assessment"):
        assert key in result


def test_scan_has_groups():
    result = scan_anomaly_free_matter_loops()
    assert result["n_groups"] >= 3


def test_scan_no_match_3alpha2():
    """No anomaly-free group should match 3*alpha^2 (all overshoot)."""
    result = scan_anomaly_free_matter_loops()
    assert result["any_match_3alpha2"] is False


def test_scan_all_groups_overshoot():
    """All anomaly-free groups have too many fields, ratio >> 1."""
    result = scan_anomaly_free_matter_loops()
    for r in result["results"]:
        assert abs(r["g_correction"]["ratio_to_target"]) > 100


# ===========================================================================
# 6. summarize_charged_matter_loops (~5 tests)
# ===========================================================================


def test_summary_returns_dict():
    result = summarize_charged_matter_loops()
    assert isinstance(result, dict)
    for key in ("monopole_spectrum", "spectral_zeta", "e8_g_correction",
                "e8_cw_potential", "group_scan", "pure_gravity_zeta",
                "monopole_zeta", "zeta_shift", "gap2_impact",
                "honest_assessment"):
        assert key in result


def test_summary_zeta_values():
    result = summarize_charged_matter_loops()
    assert abs(result["pure_gravity_zeta"] - (-17.0 / 480.0)) < 1e-10
    assert abs(result["monopole_zeta"] - 0.1) < 1e-10


def test_summary_zeta_shift_positive():
    result = summarize_charged_matter_loops()
    assert result["zeta_shift"] > 0


def test_summary_e8_computed():
    result = summarize_charged_matter_loops()
    assert result["e8_g_correction"]["delta_G_over_G_total"] != 0


def test_summary_honest_assessment_present():
    result = summarize_charged_matter_loops()
    assert len(result["honest_assessment"]) > 50


# ===========================================================================
# Cross-checks (~3 tests)
# ===========================================================================


def test_zeta_shift_is_difference():
    z0 = compute_monopole_spectral_zeta(charge_q=0, n_monopole=1, s=-1)
    z1 = compute_monopole_spectral_zeta(charge_q=1, n_monopole=1, s=-1)
    result = summarize_charged_matter_loops()
    assert abs(result["zeta_shift"] - (z1["value"] - z0["value"])) < 1e-10


def test_gravity_correction_independent_of_matter():
    """Pure gravity part should be the same regardless of matter content."""
    g1 = compute_charged_oneloop_g_correction(n_scalars=1, n_fermions=0, n_vectors=0)
    g2 = compute_charged_oneloop_g_correction(n_scalars=100, n_fermions=50, n_vectors=10)
    assert abs(g1["delta_G_over_G_gravity"] - g2["delta_G_over_G_gravity"]) < 1e-15


def test_bernoulli_polynomial_consistency():
    """Verify B_4(3/2) used in zeta_H(-3, 3/2) is correct."""
    from alpha_ladder_core.charged_matter_loops import _bernoulli_polynomial
    B4_1_5 = _bernoulli_polynomial(4, 1.5)
    expected = (1.5) ** 4 - 2 * (1.5) ** 3 + (1.5) ** 2 - 1.0 / 30.0
    assert abs(B4_1_5 - expected) < 1e-12
