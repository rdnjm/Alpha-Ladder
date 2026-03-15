"""
Comprehensive tests for alpha_ladder_core.rung_spacing module.

Tests cover:
- compute_rungs
- score_spacing
- search_rational_spacings
- search_irrational_spacings
- search_continuous_optimum
- get_best_fit_details
"""

import math
import pytest
from types import SimpleNamespace

from alpha_ladder_core.constants import get_constants, get_particle_masses
from alpha_ladder_core.rung_spacing import (
    compute_rungs,
    score_spacing,
    search_rational_spacings,
    search_irrational_spacings,
    search_continuous_optimum,
    get_best_fit_details,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def constants():
    """CODATA 2018 constants as SimpleNamespace."""
    return get_constants("CODATA 2018")


@pytest.fixture
def rung_dict(constants):
    """Rung dictionary computed from CODATA 2018 constants."""
    return compute_rungs(constants)


@pytest.fixture
def rung_values(rung_dict):
    """List of rung values (floats) for continuous optimum search."""
    return list(rung_dict.values())


# ---------------------------------------------------------------------------
# Tests: compute_rungs
# ---------------------------------------------------------------------------

class TestComputeRungs:
    def test_returns_dict(self, constants):
        result = compute_rungs(constants)
        assert isinstance(result, dict)

    def test_particle_count(self, constants):
        result = compute_rungs(constants)
        assert len(result) == 12

    def test_contains_all_particles(self, constants):
        result = compute_rungs(constants)
        masses = get_particle_masses()
        expected = {n for n in masses if n != "Electron"}
        assert set(result.keys()) == expected

    def test_all_positive_rungs(self, rung_dict):
        """All particles are heavier than the electron, so all rungs positive."""
        for name, rung in rung_dict.items():
            assert rung > 0, f"{name} rung should be positive"

    def test_proton_rung_approximate(self, rung_dict):
        assert abs(rung_dict["Proton"] - 1.53) < 0.05

    def test_top_quark_rung_approximate(self, rung_dict):
        assert abs(rung_dict["Top quark"] - 2.59) < 0.05

    def test_muon_rung_less_than_proton(self, rung_dict):
        assert rung_dict["Muon"] < rung_dict["Proton"]

    def test_heavier_particle_has_higher_rung(self, rung_dict):
        assert rung_dict["W boson"] < rung_dict["Z boson"]
        assert rung_dict["Z boson"] < rung_dict["Higgs boson"]
        assert rung_dict["Higgs boson"] < rung_dict["Top quark"]

    def test_up_quark_smallest_rung(self, rung_dict):
        """Up quark is lightest, should have smallest rung."""
        up_rung = rung_dict["Up quark"]
        for name, rung in rung_dict.items():
            if name != "Up quark":
                assert rung > up_rung, f"{name} should have higher rung than Up quark"

    def test_accepts_dict_constants(self):
        """compute_rungs should also accept dict-style constants."""
        constants_dict = {"alpha": 7.2973525693e-3}
        result = compute_rungs(constants_dict)
        assert len(result) == 12
        assert all(v > 0 for v in result.values())

    def test_accepts_simplenamespace_constants(self, constants):
        """SimpleNamespace access path works."""
        assert hasattr(constants, "alpha")
        result = compute_rungs(constants)
        assert len(result) == 12

    def test_rung_values_are_float(self, rung_dict):
        for v in rung_dict.values():
            assert isinstance(v, float)

    def test_rung_values_finite(self, rung_dict):
        for name, v in rung_dict.items():
            assert math.isfinite(v), f"{name} rung should be finite"


# ---------------------------------------------------------------------------
# Tests: score_spacing
# ---------------------------------------------------------------------------

class TestScoreSpacing:
    def test_returns_tuple_of_three(self, rung_dict):
        result = score_spacing(0.5, rung_dict)
        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_avg_closeness_range(self, rung_dict):
        avg, _, _ = score_spacing(0.5, rung_dict)
        assert 0 <= avg <= 0.5

    def test_n_matches_nonnegative(self, rung_dict):
        _, n_matches, _ = score_spacing(0.5, rung_dict)
        assert n_matches >= 0

    def test_n_matches_at_most_12(self, rung_dict):
        _, n_matches, _ = score_spacing(0.5, rung_dict)
        assert n_matches <= 12

    def test_details_length(self, rung_dict):
        _, _, details = score_spacing(0.5, rung_dict)
        assert len(details) == 12

    def test_detail_tuple_structure(self, rung_dict):
        _, _, details = score_spacing(0.5, rung_dict)
        for d in details:
            assert len(d) == 7
            name, n, nearest_rung, nearest_k, delta, frac, match = d
            assert isinstance(name, str)
            assert isinstance(n, float)
            assert isinstance(nearest_rung, float)
            assert isinstance(nearest_k, int)
            assert isinstance(delta, float)
            assert isinstance(frac, float)
            assert isinstance(match, bool)

    def test_frac_between_0_and_half(self, rung_dict):
        _, _, details = score_spacing(0.5, rung_dict)
        for d in details:
            frac = d[5]
            assert 0 <= frac <= 0.5 + 1e-10

    def test_delta_nonnegative(self, rung_dict):
        _, _, details = score_spacing(0.5, rung_dict)
        for d in details:
            assert d[4] >= 0

    def test_match_reflects_tolerance(self, rung_dict):
        _, _, details = score_spacing(0.5, rung_dict, tolerance_frac=0.15)
        for d in details:
            frac, match = d[5], d[6]
            assert match == (frac < 0.15)

    def test_custom_tolerance(self, rung_dict):
        _, n_loose, _ = score_spacing(0.5, rung_dict, tolerance_frac=0.50)
        _, n_tight, _ = score_spacing(0.5, rung_dict, tolerance_frac=0.01)
        assert n_loose >= n_tight

    def test_perfect_spacing_for_single_particle(self):
        """If spacing equals the rung value, frac should be ~0."""
        rung = {"Test": 1.5}
        avg, n_matches, details = score_spacing(0.75, rung)
        # 1.5 / 0.75 = 2.0 exactly
        assert details[0][5] < 1e-10  # frac ~0
        assert details[0][6] is True

    def test_n_matches_equals_count_of_matched(self, rung_dict):
        _, n_matches, details = score_spacing(0.5, rung_dict)
        counted = sum(1 for d in details if d[6])
        assert n_matches == counted

    def test_avg_equals_mean_frac(self, rung_dict):
        avg, _, details = score_spacing(0.5, rung_dict)
        mean = sum(d[5] for d in details) / len(details)
        assert abs(avg - mean) < 1e-12


# ---------------------------------------------------------------------------
# Tests: search_rational_spacings
# ---------------------------------------------------------------------------

class TestSearchRationalSpacings:
    def test_returns_list(self, rung_dict):
        result = search_rational_spacings(rung_dict)
        assert isinstance(result, list)

    def test_default_24_entries(self, rung_dict):
        result = search_rational_spacings(rung_dict)
        assert len(result) == 24

    def test_custom_k_max(self, rung_dict):
        result = search_rational_spacings(rung_dict, k_max=10)
        assert len(result) == 10

    def test_entry_structure(self, rung_dict):
        result = search_rational_spacings(rung_dict)
        for entry in result:
            assert len(entry) == 5
            k, spacing, avg, matches, details = entry
            assert isinstance(k, int)
            assert isinstance(spacing, float)
            assert isinstance(avg, float)
            assert isinstance(matches, int)
            assert isinstance(details, list)

    def test_sorted_by_matches_descending(self, rung_dict):
        result = search_rational_spacings(rung_dict)
        matches_list = [entry[3] for entry in result]
        for i in range(len(matches_list) - 1):
            assert matches_list[i] >= matches_list[i + 1]

    def test_spacings_are_reciprocals(self, rung_dict):
        result = search_rational_spacings(rung_dict)
        k_values = sorted([entry[0] for entry in result])
        assert k_values == list(range(1, 25))
        for entry in result:
            k, spacing = entry[0], entry[1]
            assert abs(spacing - 1.0 / k) < 1e-15

    def test_all_matches_within_range(self, rung_dict):
        result = search_rational_spacings(rung_dict)
        for entry in result:
            assert 0 <= entry[3] <= 12

    def test_details_have_12_particles(self, rung_dict):
        result = search_rational_spacings(rung_dict)
        for entry in result:
            assert len(entry[4]) == 12


# ---------------------------------------------------------------------------
# Tests: search_irrational_spacings
# ---------------------------------------------------------------------------

class TestSearchIrrationalSpacings:
    def test_returns_list(self, rung_dict):
        result = search_irrational_spacings(rung_dict)
        assert isinstance(result, list)

    def test_15_entries(self, rung_dict):
        result = search_irrational_spacings(rung_dict)
        assert len(result) == 15

    def test_entry_structure(self, rung_dict):
        result = search_irrational_spacings(rung_dict)
        for entry in result:
            assert len(entry) == 5
            name, spacing, avg, matches, details = entry
            assert isinstance(name, str)
            assert isinstance(spacing, float)
            assert isinstance(avg, float)
            assert isinstance(matches, int)
            assert isinstance(details, list)

    def test_sorted_by_matches_descending(self, rung_dict):
        result = search_irrational_spacings(rung_dict)
        matches_list = [entry[3] for entry in result]
        for i in range(len(matches_list) - 1):
            assert matches_list[i] >= matches_list[i + 1]

    def test_all_spacings_positive(self, rung_dict):
        result = search_irrational_spacings(rung_dict)
        for entry in result:
            assert entry[1] > 0

    def test_all_spacings_finite(self, rung_dict):
        result = search_irrational_spacings(rung_dict)
        for entry in result:
            assert math.isfinite(entry[1])

    def test_details_have_12_particles(self, rung_dict):
        result = search_irrational_spacings(rung_dict)
        for entry in result:
            assert len(entry[4]) == 12

    def test_known_spacing_names(self, rung_dict):
        result = search_irrational_spacings(rung_dict)
        names = {entry[0] for entry in result}
        # Check a few expected names
        assert "1/\u03c6" in names or "ln2" in names  # at least some known names present
        assert len(names) == 15  # all unique


# ---------------------------------------------------------------------------
# Tests: search_continuous_optimum
# ---------------------------------------------------------------------------

class TestSearchContinuousOptimum:
    def test_returns_tuple_of_three(self, rung_values):
        result = search_continuous_optimum(rung_values)
        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_best_spacing_in_range(self, rung_values):
        best_spacing, _, _ = search_continuous_optimum(rung_values)
        assert 0.01 <= best_spacing <= 2.0

    def test_best_score_positive_finite(self, rung_values):
        _, best_score, _ = search_continuous_optimum(rung_values)
        assert best_score > 0
        assert math.isfinite(best_score)

    def test_local_mins_is_list(self, rung_values):
        _, _, local_mins = search_continuous_optimum(rung_values)
        assert isinstance(local_mins, list)

    def test_local_mins_sorted_by_rms_ascending(self, rung_values):
        _, _, local_mins = search_continuous_optimum(rung_values)
        rms_values = [lm[1] for lm in local_mins]
        for i in range(len(rms_values) - 1):
            assert rms_values[i] <= rms_values[i + 1]

    def test_local_mins_entries_are_tuples(self, rung_values):
        _, _, local_mins = search_continuous_optimum(rung_values)
        for lm in local_mins:
            assert len(lm) == 2
            spacing, rms = lm
            assert isinstance(spacing, float)
            assert isinstance(rms, float)

    def test_best_score_leq_all_local_mins(self, rung_values):
        _, best_score, local_mins = search_continuous_optimum(rung_values)
        for lm in local_mins:
            assert best_score <= lm[1] + 1e-12

    def test_local_mins_spacings_in_range(self, rung_values):
        _, _, local_mins = search_continuous_optimum(rung_values)
        for lm in local_mins:
            assert 0.01 <= lm[0] <= 2.0

    def test_best_spacing_is_among_scanned(self, rung_values):
        best_spacing, _, _ = search_continuous_optimum(rung_values)
        # Should be a multiple of 0.001
        k = round(best_spacing * 1000)
        assert abs(best_spacing - k / 1000.0) < 1e-10

    def test_at_least_one_local_min(self, rung_values):
        _, _, local_mins = search_continuous_optimum(rung_values)
        assert len(local_mins) >= 1


# ---------------------------------------------------------------------------
# Tests: get_best_fit_details
# ---------------------------------------------------------------------------

class TestGetBestFitDetails:
    def test_returns_list(self, rung_dict):
        result = get_best_fit_details(0.5, rung_dict)
        assert isinstance(result, list)

    def test_12_entries(self, rung_dict):
        result = get_best_fit_details(0.5, rung_dict)
        assert len(result) == 12

    def test_sorted_by_rung_ascending(self, rung_dict):
        result = get_best_fit_details(0.5, rung_dict)
        rung_vals = [entry[1] for entry in result]
        for i in range(len(rung_vals) - 1):
            assert rung_vals[i] <= rung_vals[i + 1]

    def test_detail_tuple_structure(self, rung_dict):
        result = get_best_fit_details(0.5, rung_dict)
        for d in result:
            assert len(d) == 7
            name, n, nearest_rung, k, delta, frac, match = d
            assert isinstance(name, str)
            assert isinstance(frac, float)
            assert isinstance(match, bool)

    def test_first_entry_is_lightest_particle(self, rung_dict):
        """Sorted by rung, first should be lightest (Up quark)."""
        result = get_best_fit_details(0.5, rung_dict)
        assert result[0][0] == "Up quark"

    def test_last_entry_is_heaviest_particle(self, rung_dict):
        """Sorted by rung, last should be heaviest (Top quark)."""
        result = get_best_fit_details(0.5, rung_dict)
        assert result[-1][0] == "Top quark"

    def test_consistent_with_score_spacing(self, rung_dict):
        """Details should contain same data as score_spacing, just sorted."""
        details_sorted = get_best_fit_details(0.5, rung_dict)
        _, _, details_unsorted = score_spacing(0.5, rung_dict)
        # Same elements, possibly different order
        sorted_names = {d[0] for d in details_sorted}
        unsorted_names = {d[0] for d in details_unsorted}
        assert sorted_names == unsorted_names

    def test_different_spacings_give_different_results(self, rung_dict):
        r1 = get_best_fit_details(0.3, rung_dict)
        r2 = get_best_fit_details(0.7, rung_dict)
        # nearest_k values should differ for at least one particle
        k_vals_1 = [d[3] for d in r1]
        k_vals_2 = [d[3] for d in r2]
        assert k_vals_1 != k_vals_2


# ---------------------------------------------------------------------------
# Tests: edge cases and cross-function consistency
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_score_spacing_single_particle(self):
        rung = {"Electron": 0.0}
        # spacing 1.0, rung 0.0 -> nearest_k=0, delta=0, frac=0
        avg, n_matches, details = score_spacing(1.0, rung)
        assert avg == 0.0
        assert n_matches == 1
        assert details[0][5] == 0.0

    def test_score_spacing_very_small_spacing(self, rung_dict):
        avg, n_matches, _ = score_spacing(0.001, rung_dict, tolerance_frac=0.49)
        # Very small spacing means frac can still be up to 0.5, use loose tolerance
        assert n_matches >= 6

    def test_rational_and_irrational_cover_different_spacings(self, rung_dict):
        rational = search_rational_spacings(rung_dict)
        irrational = search_irrational_spacings(rung_dict)
        rational_spacings = {round(r[1], 10) for r in rational}
        irrational_spacings = {round(r[1], 10) for r in irrational}
        # No overlap expected (rational are 1/k, irrational are transcendental)
        overlap = rational_spacings & irrational_spacings
        assert len(overlap) == 0

    def test_continuous_optimum_with_trivial_input(self):
        """Single rung value at 1.0."""
        best_spacing, best_score, local_mins = search_continuous_optimum([1.0])
        assert best_score < 0.01  # should find near-perfect fit
        assert 0.01 <= best_spacing <= 2.0

    def test_compute_rungs_deterministic(self, constants):
        r1 = compute_rungs(constants)
        r2 = compute_rungs(constants)
        for name in r1:
            assert r1[name] == r2[name]
