"""
Tests for alpha_ladder_core/vacuum_polynomial.py

Vacuum polynomial analysis for Kaluza-Klein compactification.
Verifies that the characteristic polynomial x^2 + Dx + d for (d=4, n=2)
has discriminant 20, placing its roots in the golden field Q(sqrt(5)),
and that lambda_GB = sqrt(5) - 3 is the unique minimal solution.
"""

import math
import unittest

from alpha_ladder_core.vacuum_polynomial import (
    compute_dimensional_polynomial,
    verify_lambda_is_root,
    scan_discriminant_field,
    prove_golden_uniqueness,
    derive_algebraic_closure,
    summarize_vacuum_polynomial,
)

SQRT5 = math.sqrt(5)
LAMBDA_GB = SQRT5 - 3  # approx -0.7639320225


# -----------------------------------------------------------------------
# Tests for compute_dimensional_polynomial
# -----------------------------------------------------------------------

class TestDimensionalPolynomial(unittest.TestCase):
    """Tests for compute_dimensional_polynomial."""

    def test_d4_n2_discriminant(self):
        """Discriminant = D^2 - 4*d = 36 - 16 = 20."""
        result = compute_dimensional_polynomial(d=4, n=2)
        self.assertEqual(result["discriminant"], 20)

    def test_d4_n2_D(self):
        """D = d + n = 4 + 2 = 6."""
        result = compute_dimensional_polynomial(d=4, n=2)
        self.assertEqual(result["D"], 6)

    def test_d4_n2_roots(self):
        """Roots are (-D +/- sqrt(discriminant)) / 2 = -3 +/- sqrt(5)."""
        result = compute_dimensional_polynomial(d=4, n=2)
        expected_plus = -3 + SQRT5
        expected_minus = -3 - SQRT5
        self.assertAlmostEqual(result["roots"][0], expected_plus, places=10)
        self.assertAlmostEqual(result["roots"][1], expected_minus, places=10)

    def test_d4_n2_root_plus(self):
        """root_plus = sqrt(5) - 3 ~ -0.7639320225."""
        result = compute_dimensional_polynomial(d=4, n=2)
        self.assertAlmostEqual(result["root_plus"], SQRT5 - 3, places=10)

    def test_d4_n2_golden_field(self):
        """Squarefree part of discriminant 20 is 5, so is_golden_field = True."""
        result = compute_dimensional_polynomial(d=4, n=2)
        self.assertTrue(result["is_golden_field"])

    def test_d4_n3_not_golden(self):
        """For n=3: D=7, discriminant = 49 - 16 = 33, squarefree(33) = 33 != 5."""
        result = compute_dimensional_polynomial(d=4, n=3)
        self.assertEqual(result["discriminant"], 33)
        self.assertFalse(result["is_golden_field"])

    def test_d4_n4_not_golden(self):
        """For n=4: D=8, discriminant = 64 - 16 = 48, squarefree(48) = 3 != 5."""
        result = compute_dimensional_polynomial(d=4, n=4)
        self.assertEqual(result["discriminant"], 48)
        self.assertFalse(result["is_golden_field"])

    def test_coefficients(self):
        """Polynomial x^2 + 6x + 4 has coefficients [1, 6, 4]."""
        result = compute_dimensional_polynomial(d=4, n=2)
        self.assertEqual(result["coefficients"], [1, 6, 4])

    def test_returns_dict(self):
        result = compute_dimensional_polynomial(d=4, n=2)
        self.assertIsInstance(result, dict)

    def test_squarefree_20_is_5(self):
        """The squarefree part of 20 = 4*5 should be 5."""
        result = compute_dimensional_polynomial(d=4, n=2)
        # discriminant is 20, and is_golden_field is True iff squarefree(20) == 5
        self.assertEqual(result["discriminant"], 20)
        self.assertTrue(result["is_golden_field"])


# -----------------------------------------------------------------------
# Tests for verify_lambda_is_root
# -----------------------------------------------------------------------

class TestVerifyLambdaIsRoot(unittest.TestCase):
    """Tests for verify_lambda_is_root."""

    def test_is_root(self):
        """lambda_GB = sqrt(5) - 3 is a root of x^2 + 6x + 4 = 0."""
        result = verify_lambda_is_root(d=4, n=2)
        self.assertTrue(result["is_root"])

    def test_polynomial_value_near_zero(self):
        """Substituting lambda_GB into x^2 + 6x + 4 should give ~0."""
        result = verify_lambda_is_root(d=4, n=2)
        self.assertAlmostEqual(result["polynomial_value"], 0.0, places=12)

    def test_vieta_sum(self):
        """Sum of roots = -D = -6."""
        result = verify_lambda_is_root(d=4, n=2)
        root_sum = result["root_sum"] if "root_sum" in result else None
        if root_sum is not None:
            self.assertAlmostEqual(root_sum, -6.0, places=10)
        # Also check the boolean
        self.assertTrue(result["vieta_sum_check"])

    def test_vieta_product(self):
        """Product of roots = d = 4."""
        result = verify_lambda_is_root(d=4, n=2)
        root_product = result.get("root_product")
        if root_product is not None:
            self.assertAlmostEqual(root_product, 4.0, places=10)
        # Also check the boolean
        self.assertTrue(result["vieta_product_check"])

    def test_vieta_sum_check(self):
        result = verify_lambda_is_root(d=4, n=2)
        self.assertTrue(result["vieta_sum_check"])

    def test_vieta_product_check(self):
        result = verify_lambda_is_root(d=4, n=2)
        self.assertTrue(result["vieta_product_check"])

    def test_custom_lambda(self):
        """With lambda_GB = 0.5, it should NOT be a root of x^2 + 6x + 4."""
        result = verify_lambda_is_root(d=4, n=2, lambda_GB=0.5)
        self.assertFalse(result["is_root"])

    def test_default_lambda(self):
        """Default (lambda_GB=None) should use sqrt(5) - 3."""
        result = verify_lambda_is_root(d=4, n=2)
        self.assertTrue(result["is_root"])


# -----------------------------------------------------------------------
# Tests for scan_discriminant_field
# -----------------------------------------------------------------------

class TestScanDiscriminantField(unittest.TestCase):
    """Tests for scan_discriminant_field."""

    def test_n2_is_golden(self):
        """n=2 gives discriminant 20 with squarefree part 5 => Q(sqrt(5))."""
        result = scan_discriminant_field(d=4, n_max=20)
        self.assertIn(2, result["golden_field_solutions"])

    def test_n3_not_golden(self):
        """n=3 gives discriminant 33, squarefree(33) = 33 != 5."""
        result = scan_discriminant_field(d=4, n_max=20)
        self.assertNotIn(3, result["golden_field_solutions"])

    def test_minimal_n_is_2(self):
        """The smallest n producing Q(sqrt(5)) is n=2."""
        result = scan_discriminant_field(d=4, n_max=20)
        self.assertEqual(result["minimal_n"], 2)

    def test_n10_is_golden(self):
        """n=10: D=14, discriminant = 196 - 16 = 180, squarefree(180) = 5."""
        result = scan_discriminant_field(d=4, n_max=20)
        self.assertIn(10, result["golden_field_solutions"])

    def test_unique_minimal(self):
        """n=2 is the unique minimal solution."""
        result = scan_discriminant_field(d=4, n_max=20)
        self.assertTrue(result["is_n2_unique_minimal"])

    def test_results_length(self):
        """Results should have n_max - 1 entries (n=2 through n_max)."""
        n_max = 20
        result = scan_discriminant_field(d=4, n_max=n_max)
        self.assertEqual(len(result["results"]), n_max - 1)

    def test_all_viable_for_n_ge_2(self):
        """All entries with n >= 2 should have gb_viable = True."""
        result = scan_discriminant_field(d=4, n_max=10)
        for entry in result["results"]:
            if entry["n"] >= 2:
                self.assertTrue(entry["gb_viable"])


# -----------------------------------------------------------------------
# Tests for prove_golden_uniqueness
# -----------------------------------------------------------------------

class TestProveGoldenUniqueness(unittest.TestCase):
    """Tests for prove_golden_uniqueness via Diophantine equation n(n+8) = 5*m^2."""

    def test_minimal_solution_n2_m2(self):
        """Minimal solution is (n, m) = (2, 2): 2*10 = 20 = 5*4."""
        result = prove_golden_uniqueness(d=4, n_max=50)
        self.assertEqual(result["minimal_solution"], (2, 2))

    def test_n2_satisfies_diophantine(self):
        """Verify 2*(2+8) = 20 = 5*2^2."""
        n, m = 2, 2
        self.assertEqual(n * (n + 8), 5 * m**2)

    def test_next_solution_n10(self):
        """Next solution is (n, m) = (10, 6): 10*18 = 180 = 5*36."""
        result = prove_golden_uniqueness(d=4, n_max=50)
        self.assertEqual(result["next_solution"], (10, 6))

    def test_n2_is_minimal(self):
        result = prove_golden_uniqueness(d=4, n_max=50)
        self.assertTrue(result["n2_is_minimal"])

    def test_has_physical_arguments(self):
        """Physical arguments should be a non-empty list."""
        result = prove_golden_uniqueness(d=4, n_max=50)
        self.assertIsInstance(result["physical_arguments"], list)
        self.assertGreater(len(result["physical_arguments"]), 0)


# -----------------------------------------------------------------------
# Tests for derive_algebraic_closure
# -----------------------------------------------------------------------

class TestAlgebraicClosure(unittest.TestCase):
    """Tests for derive_algebraic_closure."""

    def test_all_in_same_field(self):
        result = derive_algebraic_closure(d=4, n=2)
        self.assertTrue(result["all_in_same_field"])

    def test_phi_in_constants(self):
        """The golden ratio phi should appear in the constants list."""
        result = derive_algebraic_closure(d=4, n=2)
        names = [c["name"].lower() for c in result["constants"]]
        self.assertTrue(
            any("golden" in name or "phi" in name for name in names),
            f"Expected 'golden ratio' or 'phi' in constants, got {names}",
        )

    def test_lambda_gb_in_constants(self):
        """lambda_GB should appear in the constants list."""
        result = derive_algebraic_closure(d=4, n=2)
        names = [c["name"].lower() for c in result["constants"]]
        self.assertTrue(
            any("lambda" in name for name in names),
            f"Expected 'lambda' in constants, got {names}",
        )

    def test_omega_target_in_constants(self):
        """omega_target should appear in the constants list."""
        result = derive_algebraic_closure(d=4, n=2)
        names = [c["name"].lower() for c in result["constants"]]
        self.assertTrue(
            any("omega" in name for name in names),
            f"Expected 'omega' in constants, got {names}",
        )

    def test_bridge_in_constants(self):
        """phi^2/2 (the bridge constant) should appear."""
        result = derive_algebraic_closure(d=4, n=2)
        names = [c["name"].lower() for c in result["constants"]]
        self.assertTrue(
            any("bridge" in name or "phi^2/2" in name or "phi²/2" in name for name in names),
            f"Expected bridge constant in constants, got {names}",
        )

    def test_all_marked_in_field(self):
        """Every constant should have in_field = True."""
        result = derive_algebraic_closure(d=4, n=2)
        for c in result["constants"]:
            self.assertTrue(
                c["in_field"],
                f"Expected {c['name']} to be in Q(sqrt(5))",
            )

    def test_at_least_5_constants(self):
        result = derive_algebraic_closure(d=4, n=2)
        self.assertGreaterEqual(len(result["constants"]), 5)

    def test_phi_minimal_poly(self):
        """phi's minimal polynomial is x^2 - x - 1."""
        result = derive_algebraic_closure(d=4, n=2)
        phi_entries = [
            c for c in result["constants"]
            if "phi" in c["name"].lower() or "golden" in c["name"].lower()
        ]
        self.assertGreater(len(phi_entries), 0)
        phi_entry = phi_entries[0]
        min_poly = phi_entry.get("minimal_polynomial", "")
        self.assertTrue(
            "x^2 - x - 1" in min_poly or "x² - x - 1" in min_poly,
            f"Expected 'x^2 - x - 1' in minimal polynomial, got '{min_poly}'",
        )

    def test_lambda_minimal_poly_coefficients(self):
        """lambda_GB's minimal polynomial should involve D and d."""
        result = derive_algebraic_closure(d=4, n=2)
        lambda_entries = [
            c for c in result["constants"]
            if "lambda" in c["name"].lower()
        ]
        self.assertGreater(len(lambda_entries), 0)
        min_poly = lambda_entries[0].get("minimal_polynomial", "")
        self.assertGreater(len(min_poly), 0)


# -----------------------------------------------------------------------
# Tests for summarize_vacuum_polynomial
# -----------------------------------------------------------------------

class TestSummarizeVacuumPolynomial(unittest.TestCase):
    """Tests for summarize_vacuum_polynomial."""

    def test_lambda_is_root(self):
        result = summarize_vacuum_polynomial(d=4, n=2)
        self.assertTrue(result["lambda_GB_is_root"])

    def test_golden_field(self):
        result = summarize_vacuum_polynomial(d=4, n=2)
        self.assertTrue(result["golden_field"])

    def test_n2_unique_minimal(self):
        result = summarize_vacuum_polynomial(d=4, n=2)
        self.assertTrue(result["n2_unique_minimal"])

    def test_has_summary(self):
        result = summarize_vacuum_polynomial(d=4, n=2)
        self.assertIsInstance(result["summary"], str)
        self.assertGreater(len(result["summary"]), 0)

    def test_has_key_result(self):
        result = summarize_vacuum_polynomial(d=4, n=2)
        self.assertIsInstance(result["key_result"], str)
        self.assertGreater(len(result["key_result"]), 0)

    def test_returns_dict(self):
        result = summarize_vacuum_polynomial(d=4, n=2)
        self.assertIsInstance(result, dict)


if __name__ == "__main__":
    unittest.main()
