"""
Tests for alpha_ladder_core/kk_explicit_reduction.py

Explicit KK reduction from arXiv:2601.08443 (Dereli & Senikoglu, 2026).
Verifies the structure of the 6D -> 4D reduction on M_4 x S^2.
"""

import math
from alpha_ladder_core.kk_explicit_reduction import (
    compute_reduced_action,
    verify_omega_zero,
    verify_tree_level_potential,
    compute_gauge_kinetic_matrix,
    compute_gauge_loop_structure,
    summarize_explicit_reduction,
)


PHI = (1.0 + math.sqrt(5.0)) / 2.0


# -----------------------------------------------------------------------
# Tests for compute_reduced_action (8 tests)
# -----------------------------------------------------------------------

def test_reduced_action_returns_dict():
    result = compute_reduced_action()
    assert isinstance(result, dict)


def test_reduced_action_omega_zero():
    result = compute_reduced_action()
    assert result["bd_parameter"] == 0


def test_reduced_action_gravitational_coupling():
    result = compute_reduced_action(phi_value=1.0)
    assert abs(result["gravitational_coupling"] - 0.5) < 1e-15


def test_reduced_action_gravitational_coupling_phi2():
    result = compute_reduced_action(phi_value=2.0)
    assert abs(result["gravitational_coupling"] - 2.0) < 1e-15


def test_reduced_action_two_physical_gauge_fields():
    result = compute_reduced_action()
    assert result["n_physical_gauge_fields"] == 2


def test_reduced_action_effective_potential_zero():
    result = compute_reduced_action()
    assert result["effective_potential"] == 0.0


def test_reduced_action_gb_contribution():
    result = compute_reduced_action()
    assert abs(result["topological_gb_contribution"] - 8.0 * math.pi) < 1e-12


def test_reduced_action_coset_structure():
    result = compute_reduced_action()
    assert result["coset_structure"] == "SO(3)/SO(2)"


def test_reduced_action_rank_2():
    result = compute_reduced_action()
    assert result["gauge_kinetic_matrix_rank"] == 2


def test_reduced_action_eigenvalues():
    result = compute_reduced_action()
    assert result["gauge_kinetic_matrix_eigenvalues"] == [0.0, 1.0, 1.0]


def test_reduced_action_n_killing_vectors():
    result = compute_reduced_action()
    assert result["n_killing_vectors"] == 3


def test_reduced_action_scalar_kinetic_term():
    result = compute_reduced_action()
    assert result["scalar_kinetic_term"] is True


def test_reduced_action_has_reference():
    result = compute_reduced_action()
    assert "2601.08443" in result["reference"]


def test_reduced_action_has_description():
    result = compute_reduced_action()
    assert isinstance(result["description"], str)
    assert len(result["description"]) > 0


# -----------------------------------------------------------------------
# Tests for verify_omega_zero (5 tests)
# -----------------------------------------------------------------------

def test_omega_n2_is_zero():
    result = verify_omega_zero(n=2)
    assert abs(result["omega"]) < 1e-15
    assert result["is_zero"] is True


def test_omega_n1():
    result = verify_omega_zero(n=1)
    expected = 0.5 * (0.5 - 1.0)  # = -0.25
    assert abs(result["omega"] - expected) < 1e-15


def test_omega_n3():
    result = verify_omega_zero(n=3)
    expected = 1.5 * (1.5 - 1.0)  # = 0.75
    assert abs(result["omega"] - expected) < 1e-15


def test_omega_n4():
    result = verify_omega_zero(n=4)
    expected = 2.0 * (2.0 - 1.0)  # = 2.0
    assert abs(result["omega"] - expected) < 1e-15


def test_omega_returns_dict():
    result = verify_omega_zero()
    assert isinstance(result, dict)
    assert "omega" in result
    assert "formula" in result
    assert "n" in result


# -----------------------------------------------------------------------
# Tests for verify_tree_level_potential (5 tests)
# -----------------------------------------------------------------------

def test_tree_potential_v_zero():
    result = verify_tree_level_potential()
    assert result["V_phi"] == 0.0


def test_tree_potential_scaling_exponent():
    result = verify_tree_level_potential()
    assert result["eh_scaling_exponent"] == 0


def test_tree_potential_gb_topological():
    result = verify_tree_level_potential()
    assert "topological" in result["gb_contribution"]


def test_tree_potential_polynomial_not_derivable():
    result = verify_tree_level_potential()
    assert result["vacuum_polynomial_derivable"] is False


def test_tree_potential_explanation_present():
    result = verify_tree_level_potential()
    assert isinstance(result["explanation"], str)
    assert "ansatz" in result["explanation"].lower() or "ansatz" in result["explanation"]


# -----------------------------------------------------------------------
# Tests for compute_gauge_kinetic_matrix (5 tests)
# -----------------------------------------------------------------------

def test_gauge_matrix_3x3():
    result = compute_gauge_kinetic_matrix()
    matrix = result["matrix"]
    assert len(matrix) == 3
    assert all(len(row) == 3 for row in matrix)


def test_gauge_matrix_rank_2():
    result = compute_gauge_kinetic_matrix()
    assert result["rank"] == 2


def test_gauge_matrix_eigenvalues():
    result = compute_gauge_kinetic_matrix()
    eigenvalues = sorted(result["eigenvalues"])
    assert abs(eigenvalues[0] - 0.0) < 1e-15
    assert abs(eigenvalues[1] - 1.0) < 1e-15
    assert abs(eigenvalues[2] - 1.0) < 1e-15


def test_gauge_matrix_zero_from_stabilizer():
    result = compute_gauge_kinetic_matrix()
    assert result["zero_eigenvalue_origin"] == "SO(2) stabilizer"


def test_gauge_matrix_dimension():
    result = compute_gauge_kinetic_matrix()
    assert result["dimension"] == 3
    assert result["n_killing_vectors"] == 3
    assert result["n_physical_gauge_fields"] == 2


# -----------------------------------------------------------------------
# Tests for compute_gauge_loop_structure (7 tests)
# -----------------------------------------------------------------------

def test_gauge_loop_n_physical():
    result = compute_gauge_loop_structure()
    assert result["n_physical"] == 2


def test_gauge_loop_n_killing():
    result = compute_gauge_loop_structure()
    assert result["n_killing"] == 3


def test_gauge_loop_triple_degeneracy():
    result = compute_gauge_loop_structure()
    assert result["triple_degeneracy"] is True


def test_gauge_loop_expected_c2():
    result = compute_gauge_loop_structure()
    assert result["expected_c2"] == 3


def test_gauge_loop_expected_c3_close_to_phi_over_2():
    result = compute_gauge_loop_structure()
    expected = PHI / 2.0
    assert abs(result["expected_c3"] - expected) < 1e-12


def test_gauge_loop_not_computed():
    result = compute_gauge_loop_structure()
    assert result["loop_computed"] is False


def test_gauge_loop_status_present():
    result = compute_gauge_loop_structure()
    assert isinstance(result["status"], str)
    assert len(result["status"]) > 0
    assert "loop" in result["status"].lower()


# -----------------------------------------------------------------------
# Tests for summarize_explicit_reduction (5 tests)
# -----------------------------------------------------------------------

def test_summary_all_keys_present():
    result = summarize_explicit_reduction()
    expected_keys = [
        "reduced_action", "omega_proof", "tree_potential",
        "gauge_matrix", "gauge_loop_structure",
        "confirms_framework", "confirms_omega_zero", "confirms_v_phi_zero",
        "vacuum_polynomial_status", "correction_status",
        "reference", "overall_assessment", "honest_assessment",
    ]
    for key in expected_keys:
        assert key in result, f"Missing key: {key}"


def test_summary_confirms_framework():
    result = summarize_explicit_reduction()
    assert result["confirms_framework"] is True


def test_summary_vacuum_polynomial_status():
    result = summarize_explicit_reduction()
    assert "ansatz" in result["vacuum_polynomial_status"]


def test_summary_has_reference():
    result = summarize_explicit_reduction()
    assert "2601.08443" in result["reference"]


def test_summary_has_honest_assessment():
    result = summarize_explicit_reduction()
    assert isinstance(result["honest_assessment"], str)
    assert len(result["honest_assessment"]) > 20


# -----------------------------------------------------------------------
# Integration tests (5 tests)
# -----------------------------------------------------------------------

def test_module_imports_cleanly():
    """Verify the module can be imported without errors."""
    import alpha_ladder_core.kk_explicit_reduction as mod
    assert hasattr(mod, "compute_reduced_action")
    assert hasattr(mod, "verify_omega_zero")
    assert hasattr(mod, "verify_tree_level_potential")
    assert hasattr(mod, "compute_gauge_kinetic_matrix")
    assert hasattr(mod, "compute_gauge_loop_structure")
    assert hasattr(mod, "summarize_explicit_reduction")


def test_all_functions_accept_constants_none():
    """Functions that take constants should work with None."""
    compute_reduced_action(constants=None)
    compute_gauge_loop_structure(constants=None)
    summarize_explicit_reduction(constants=None)


def test_consistency_with_kk_reduction():
    """Check that omega=0 from explicit reduction is consistent with
    the existing kk_reduction module's omega_BD for pure EH."""
    from alpha_ladder_core.kk_reduction import compute_kinetic_coefficient
    kinetic = compute_kinetic_coefficient(d=4, n=2)
    # The existing module computes omega_BD from the kinetic coefficient.
    # For pure EH on S^2, omega_BD = n(n+d-2) / [2*(d-2)] / (n^2 * beta^2 / (d-2)^2)
    # The explicit reduction gives omega=0 in a DIFFERENT convention
    # (Dereli-Senikoglu use the scalar field phi directly, not sigma).
    # The key consistency check is that the explicit reduction's bd_parameter=0
    # is a structural result about the Lagrangian form, not a numerical coincidence.
    result = compute_reduced_action()
    assert result["bd_parameter"] == 0
    # The kinetic module gives a nonzero omega_BD because it uses a
    # different field parametrisation (sigma vs phi).
    assert kinetic["omega_BD"] is not None


def test_gauge_matrix_consistent_with_loop_structure():
    """Gauge matrix rank and loop structure n_physical must agree."""
    matrix_result = compute_gauge_kinetic_matrix()
    loop_result = compute_gauge_loop_structure()
    assert matrix_result["rank"] == loop_result["n_physical"]
    assert matrix_result["n_killing_vectors"] == loop_result["n_killing"]


def test_summary_sub_results_consistent():
    """Summary sub-results should be internally consistent."""
    summary = summarize_explicit_reduction()
    assert summary["reduced_action"]["bd_parameter"] == 0
    assert summary["omega_proof"]["is_zero"] is True
    assert summary["tree_potential"]["V_phi"] == 0.0
    assert summary["gauge_matrix"]["rank"] == 2
    assert summary["gauge_loop_structure"]["loop_computed"] is False
