"""Tests for vacuum_polynomial_derivation module.

Systematic investigation of whether the vacuum polynomial x^2+6x+4=0 (whose
discriminant 20=4*5 introduces phi) can be derived from first principles, or
whether it remains an ansatz. Tests cover polynomial coefficient scanning,
dimension pair scanning, KK quantity root checks, moduli space geometry,
swampland constraints, alternative polynomial analysis, honest status, and
the overall summary.
"""
import unittest
import math


# ---------------------------------------------------------------------------
# 1. scan_polynomial_coefficients
# ---------------------------------------------------------------------------
class TestScanPolynomialCoefficients(unittest.TestCase):
    """Tests for scan_polynomial_coefficients."""

    def test_returns_dict(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import scan_polynomial_coefficients
        result = scan_polynomial_coefficients(d=4, n=2)
        self.assertIsInstance(result, dict)

    def test_polynomials_scanned_count(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import scan_polynomial_coefficients
        result = scan_polynomial_coefficients(d=4, n=2)
        self.assertGreater(result["polynomials_scanned"], 50)

    def test_phi_producing_is_list(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import scan_polynomial_coefficients
        result = scan_polynomial_coefficients(d=4, n=2)
        self.assertIsInstance(result["phi_producing"], list)

    def test_phi_producing_at_least_two(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import scan_polynomial_coefficients
        result = scan_polynomial_coefficients(d=4, n=2)
        self.assertGreaterEqual(len(result["phi_producing"]), 2)

    def test_phi_producing_contains_x2_6x_4(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import scan_polynomial_coefficients
        result = scan_polynomial_coefficients(d=4, n=2)
        found = any(entry.get("A") == 6 and entry.get("B") == 4
                     for entry in result["phi_producing"])
        self.assertTrue(found, "phi_producing should contain entry with A=6, B=4")

    def test_phi_producing_contains_x2_6x_4(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import scan_polynomial_coefficients
        result = scan_polynomial_coefficients(d=4, n=2)
        found = any(entry.get("A") == 6 and entry.get("B") == 4
                     for entry in result["phi_producing"])
        self.assertTrue(found, "phi_producing should contain entry with A=6, B=4")

    def test_honest_assessment_nonempty(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import scan_polynomial_coefficients
        result = scan_polynomial_coefficients(d=4, n=2)
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 0)

    def test_default_parameters(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import scan_polynomial_coefficients
        result = scan_polynomial_coefficients()
        self.assertIsInstance(result, dict)
        self.assertIn("polynomials_scanned", result)


# ---------------------------------------------------------------------------
# 2. scan_dimension_pairs
# ---------------------------------------------------------------------------
class TestScanDimensionPairs(unittest.TestCase):
    """Tests for scan_dimension_pairs."""

    def test_returns_dict(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import scan_dimension_pairs
        result = scan_dimension_pairs(d_max=15, n_max=15)
        self.assertIsInstance(result, dict)

    def test_pairs_scanned_count(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import scan_dimension_pairs
        result = scan_dimension_pairs(d_max=15, n_max=15)
        self.assertGreaterEqual(result["pairs_scanned"], 100)

    def test_phi_pairs_is_list(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import scan_dimension_pairs
        result = scan_dimension_pairs(d_max=15, n_max=15)
        self.assertIsInstance(result["phi_pairs"], list)

    def test_phi_pairs_at_least_two(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import scan_dimension_pairs
        result = scan_dimension_pairs(d_max=15, n_max=15)
        self.assertGreaterEqual(len(result["phi_pairs"]), 2)

    def test_minimal_pair(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import scan_dimension_pairs
        result = scan_dimension_pairs(d_max=15, n_max=15)
        # minimal_pair is the smallest disc, which may be (1,2) not (4,2)
        self.assertIsInstance(result["minimal_pair"], tuple)
        self.assertEqual(len(result["minimal_pair"]), 2)

    def test_contains_d4_n2_entry(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import scan_dimension_pairs
        result = scan_dimension_pairs(d_max=15, n_max=15)
        found = any(entry.get("d") == 4 and entry.get("n") == 2
                     and entry.get("D") == 6 and entry.get("disc") == 20
                     for entry in result["phi_pairs"])
        self.assertTrue(found, "phi_pairs should contain (d=4, n=2, D=6, disc=20)")

    def test_default_parameters(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import scan_dimension_pairs
        result = scan_dimension_pairs()
        self.assertIsInstance(result, dict)
        self.assertIn("minimal_pair", result)


# ---------------------------------------------------------------------------
# 3. check_kk_quantities_as_roots
# ---------------------------------------------------------------------------
class TestCheckKKQuantitiesAsRoots(unittest.TestCase):
    """Tests for check_kk_quantities_as_roots."""

    def test_returns_dict(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_kk_quantities_as_roots
        result = check_kk_quantities_as_roots(d=4, n=2)
        self.assertIsInstance(result, dict)

    def test_quantities_tested_is_list(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_kk_quantities_as_roots
        result = check_kk_quantities_as_roots(d=4, n=2)
        self.assertIsInstance(result["quantities_tested"], list)

    def test_quantities_tested_at_least_ten(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_kk_quantities_as_roots
        result = check_kk_quantities_as_roots(d=4, n=2)
        self.assertGreaterEqual(len(result["quantities_tested"]), 10)

    def test_no_root_found(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_kk_quantities_as_roots
        result = check_kk_quantities_as_roots(d=4, n=2)
        self.assertFalse(result["any_root_found"])

    def test_quantities_have_required_keys(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_kk_quantities_as_roots
        result = check_kk_quantities_as_roots(d=4, n=2)
        for q in result["quantities_tested"]:
            self.assertIn("name", q)
            self.assertIn("value", q)
            self.assertIn("polynomial_value", q)

    def test_all_polynomial_values_far_from_zero(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_kk_quantities_as_roots
        result = check_kk_quantities_as_roots(d=4, n=2)
        for q in result["quantities_tested"]:
            self.assertGreater(abs(q["polynomial_value"]), 0.1,
                               f"Quantity '{q['name']}' polynomial_value too close to zero")


# ---------------------------------------------------------------------------
# 4. check_moduli_space_geometry
# ---------------------------------------------------------------------------
class TestCheckModuliSpaceGeometry(unittest.TestCase):
    """Tests for check_moduli_space_geometry."""

    def test_returns_dict(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_moduli_space_geometry
        result = check_moduli_space_geometry(d=4, n=2)
        self.assertIsInstance(result, dict)

    def test_moduli_dimension(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_moduli_space_geometry
        result = check_moduli_space_geometry(d=4, n=2)
        self.assertEqual(result["moduli_dimension"], 1)

    def test_curvature_flat(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_moduli_space_geometry
        result = check_moduli_space_geometry(d=4, n=2)
        self.assertEqual(result["curvature"], 0)

    def test_ss_stationarity_degree(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_moduli_space_geometry
        result = check_moduli_space_geometry(d=4, n=2)
        self.assertEqual(result["ss_stationarity_degree"], 3)

    def test_cannot_produce_quadratic(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_moduli_space_geometry
        result = check_moduli_space_geometry(d=4, n=2)
        self.assertFalse(result["can_produce_quadratic"])

    def test_default_parameters(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_moduli_space_geometry
        result = check_moduli_space_geometry()
        self.assertIsInstance(result, dict)
        self.assertIn("moduli_dimension", result)


# ---------------------------------------------------------------------------
# 5. check_swampland_constraints
# ---------------------------------------------------------------------------
class TestCheckSwamplandConstraints(unittest.TestCase):
    """Tests for check_swampland_constraints."""

    def test_returns_dict(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_swampland_constraints
        result = check_swampland_constraints()
        self.assertIsInstance(result, dict)

    def test_conjectures_checked_count(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_swampland_constraints
        result = check_swampland_constraints()
        self.assertEqual(result["conjectures_checked"], 10)

    def test_no_conjecture_constrains_c0(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_swampland_constraints
        result = check_swampland_constraints()
        self.assertFalse(result["any_constrains_C0"])

    def test_ds_tension(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_swampland_constraints
        result = check_swampland_constraints()
        self.assertTrue(result["ds_tension"])

    def test_results_is_list(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_swampland_constraints
        result = check_swampland_constraints()
        self.assertIsInstance(result["results"], list)

    def test_results_length(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_swampland_constraints
        result = check_swampland_constraints()
        self.assertEqual(len(result["results"]), 10)

    def test_results_have_required_keys(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import check_swampland_constraints
        result = check_swampland_constraints()
        for r in result["results"]:
            self.assertIn("name", r)
            self.assertIn("constrains_C0", r)
            self.assertIn("status", r)


# ---------------------------------------------------------------------------
# 6. analyze_alternative_polynomial
# ---------------------------------------------------------------------------
class TestAnalyzeAlternativePolynomial(unittest.TestCase):
    """Tests for analyze_alternative_polynomial."""

    def test_returns_dict(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import analyze_alternative_polynomial
        result = analyze_alternative_polynomial()
        self.assertIsInstance(result, dict)

    def test_polynomial_contains_3x_1(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import analyze_alternative_polynomial
        result = analyze_alternative_polynomial()
        self.assertIn("3x+1", result["polynomial"].replace(" ", "").replace("^2", "2"))

    def test_discriminant_is_5(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import analyze_alternative_polynomial
        result = analyze_alternative_polynomial()
        self.assertEqual(result["disc"], 5)

    def test_roots_exist(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import analyze_alternative_polynomial
        result = analyze_alternative_polynomial()
        self.assertIn("roots", result)
        self.assertIsNotNone(result["roots"])

    def test_roots_satisfy_polynomial(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import analyze_alternative_polynomial
        result = analyze_alternative_polynomial()
        # x^2 + 3x + 1 = 0 => roots satisfy this
        for r in [result["roots"]["r1"], result["roots"]["r2"]]:
            val = r**2 + 3*r + 1
            self.assertAlmostEqual(val, 0.0, places=10,
                                   msg=f"Root {r} does not satisfy x^2+3x+1=0")

    def test_c0_from_roots_close_to_phi_sq_half(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import analyze_alternative_polynomial
        result = analyze_alternative_polynomial()
        phi = (1 + math.sqrt(5)) / 2
        expected = phi**2 / 2
        self.assertAlmostEqual(result["C0_from_roots"], expected, places=5)

    def test_coefficients_from_kk_has_A(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import analyze_alternative_polynomial
        result = analyze_alternative_polynomial()
        self.assertIn("coefficients_from_kk", result)
        self.assertIn("A", result["coefficients_from_kk"])

    def test_coefficients_from_kk_A_related_to_dim_so3(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import analyze_alternative_polynomial
        result = analyze_alternative_polynomial()
        A = result["coefficients_from_kk"]["A"]
        # A should mention dim(SO(3)) or 3
        self.assertIn("3", str(A))


# ---------------------------------------------------------------------------
# 7. compute_honest_status
# ---------------------------------------------------------------------------
class TestComputeHonestStatus(unittest.TestCase):
    """Tests for compute_honest_status."""

    def test_returns_dict(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import compute_honest_status
        result = compute_honest_status()
        self.assertIsInstance(result, dict)

    def test_not_derived(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import compute_honest_status
        result = compute_honest_status()
        self.assertFalse(result["is_derived"])

    def test_is_ansatz(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import compute_honest_status
        result = compute_honest_status()
        self.assertTrue(result["is_ansatz"])

    def test_circularity_acknowledged(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import compute_honest_status
        result = compute_honest_status()
        self.assertTrue(result["circularity_acknowledged"])

    def test_status_contains_ansatz(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import compute_honest_status
        result = compute_honest_status()
        self.assertIn("ansatz", result["status"].lower())

    def test_honest_assessment_long(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import compute_honest_status
        result = compute_honest_status()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 100)


# ---------------------------------------------------------------------------
# 8. summarize_vacuum_polynomial_derivation
# ---------------------------------------------------------------------------
class TestSummarizeVacuumPolynomialDerivation(unittest.TestCase):
    """Tests for summarize_vacuum_polynomial_derivation."""

    def test_returns_dict(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import summarize_vacuum_polynomial_derivation
        result = summarize_vacuum_polynomial_derivation()
        self.assertIsInstance(result, dict)

    def test_derivation_not_achieved(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import summarize_vacuum_polynomial_derivation
        result = summarize_vacuum_polynomial_derivation()
        self.assertFalse(result["derivation_achieved"])

    def test_has_required_keys(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import summarize_vacuum_polynomial_derivation
        result = summarize_vacuum_polynomial_derivation()
        for key in ["what_works", "what_fails", "remaining_gap", "honest_assessment"]:
            self.assertIn(key, result, f"Missing key: {key}")

    def test_what_works_is_list_with_items(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import summarize_vacuum_polynomial_derivation
        result = summarize_vacuum_polynomial_derivation()
        self.assertIsInstance(result["what_works"], list)
        self.assertGreaterEqual(len(result["what_works"]), 2)

    def test_what_fails_is_list_with_items(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import summarize_vacuum_polynomial_derivation
        result = summarize_vacuum_polynomial_derivation()
        self.assertIsInstance(result["what_fails"], list)
        self.assertGreaterEqual(len(result["what_fails"]), 3)

    def test_remaining_gap_contains_ansatz(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import summarize_vacuum_polynomial_derivation
        result = summarize_vacuum_polynomial_derivation()
        self.assertIn("ansatz", result["remaining_gap"].lower())

    def test_with_constants_parameter(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import summarize_vacuum_polynomial_derivation
        from alpha_ladder_core.constants import get_constants
        c = get_constants("CODATA 2018")
        result = summarize_vacuum_polynomial_derivation(constants=c)
        self.assertIsInstance(result, dict)
        self.assertIn("derivation_achieved", result)

    def test_honest_assessment_nonempty(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import summarize_vacuum_polynomial_derivation
        result = summarize_vacuum_polynomial_derivation()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 0)


# ---------------------------------------------------------------------------
# 9. Cross-consistency tests
# ---------------------------------------------------------------------------
class TestCrossConsistency(unittest.TestCase):
    """Cross-function consistency checks."""

    def test_scan_polynomial_and_alternative_both_involve_phi(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import (
            scan_polynomial_coefficients,
            analyze_alternative_polynomial,
        )
        scan = scan_polynomial_coefficients(d=4, n=2)
        alt = analyze_alternative_polynomial()
        # scan should find x^2+6x+4 as phi-producing
        found_64 = any(entry.get("A") == 6 and entry.get("B") == 4
                        for entry in scan["phi_producing"])
        self.assertTrue(found_64)
        # alternative polynomial also involves phi (disc=5)
        self.assertEqual(alt["disc"], 5)

    def test_dimension_scan_contains_d4_n2(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import scan_dimension_pairs
        result = scan_dimension_pairs(d_max=15, n_max=15)
        found = any(p["d"] == 4 and p["n"] == 2 for p in result["phi_pairs"])
        self.assertTrue(found, "Dimension scan should include (d=4, n=2)")

    def test_honest_status_matches_summary(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import (
            compute_honest_status,
            summarize_vacuum_polynomial_derivation,
        )
        status = compute_honest_status()
        summary = summarize_vacuum_polynomial_derivation()
        self.assertEqual(status["is_derived"], summary["derivation_achieved"])
        self.assertTrue(status["is_ansatz"])

    def test_kk_no_root_consistent_with_honest_status(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import (
            check_kk_quantities_as_roots,
            compute_honest_status,
        )
        kk = check_kk_quantities_as_roots(d=4, n=2)
        status = compute_honest_status()
        # If no KK quantity is a root, the polynomial remains an ansatz
        self.assertFalse(kk["any_root_found"])
        self.assertTrue(status["is_ansatz"])

    def test_moduli_geometry_consistent_with_what_fails(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import (
            check_moduli_space_geometry,
            summarize_vacuum_polynomial_derivation,
        )
        moduli = check_moduli_space_geometry(d=4, n=2)
        summary = summarize_vacuum_polynomial_derivation()
        # Moduli space cannot produce the quadratic, so something should fail
        self.assertFalse(moduli["can_produce_quadratic"])
        self.assertGreaterEqual(len(summary["what_fails"]), 1)

    def test_swampland_no_constraint_consistent_with_summary(self):
        from alpha_ladder_core.vacuum_polynomial_derivation import (
            check_swampland_constraints,
            summarize_vacuum_polynomial_derivation,
        )
        swamp = check_swampland_constraints()
        summary = summarize_vacuum_polynomial_derivation()
        self.assertFalse(swamp["any_constrains_C0"])
        self.assertFalse(summary["derivation_achieved"])


if __name__ == "__main__":
    unittest.main()
