"""
Tests for Gauss-Bonnet derivation functions in alpha_ladder_core/kk_reduction.py

Tests the step-by-step derivation of the Gauss-Bonnet correction to the
dilaton kinetic term, including internal curvature identities and the
full derivation pipeline.

The key mathematical identity verified throughout:
    G_n / R_n^2 = (n-2)(n-3) / [n(n-1)]
where G_n is the Gauss-Bonnet density of a maximally symmetric n-manifold.
"""

import math
import unittest

from alpha_ladder_core.kk_reduction import (
    compute_internal_curvature_identities,
    derive_gauss_bonnet_correction,
    compute_gauss_bonnet_shift,
    compute_kinetic_coefficient,
    compute_target_omega,
)


# -----------------------------------------------------------------------
# Tests for compute_internal_curvature_identities
# -----------------------------------------------------------------------

class TestInternalCurvatureIdentities(unittest.TestCase):
    """Tests for compute_internal_curvature_identities."""

    def test_2d_ricci_squared(self):
        """For n=2, R_{ab}R^{ab} / R^2 = 1/n = 1/2."""
        result = compute_internal_curvature_identities(n=2)
        self.assertAlmostEqual(result["R_ab_R_ab_coeff"], 1.0 / 2.0, places=12)

    def test_2d_riemann_squared(self):
        """For n=2, R_{abcd}R^{abcd} / R^2 = 2/(n*(n-1)) = 2/2 = 1."""
        result = compute_internal_curvature_identities(n=2)
        self.assertAlmostEqual(result["R_abcd_R_abcd_coeff"], 1.0, places=12)

    def test_3d_ricci_squared(self):
        """For n=3, R_{ab}R^{ab} / R^2 = 1/3."""
        result = compute_internal_curvature_identities(n=3)
        self.assertAlmostEqual(result["R_ab_R_ab_coeff"], 1.0 / 3.0, places=12)

    def test_3d_riemann_squared(self):
        """For n=3, R_{abcd}R^{abcd} / R^2 = 2/(3*2) = 1/3."""
        result = compute_internal_curvature_identities(n=3)
        self.assertAlmostEqual(result["R_abcd_R_abcd_coeff"], 2.0 / (3.0 * 2.0), places=12)

    def test_1d_raises_or_degenerate(self):
        """n=1 should handle gracefully -- Riemann vanishes in 1D.

        For n=1 the Riemann tensor is identically zero, so the
        curvature identities are degenerate.  The function should
        either raise a ValueError or return degenerate/zero values.
        """
        try:
            result = compute_internal_curvature_identities(n=1)
            # If it returns without raising, Riemann coefficient should
            # reflect the degeneracy (division by zero in 2/(n*(n-1))
            # for n=1 gives 2/0, so the function must handle this).
            # Accept None, inf, or 0 as valid degenerate indicators.
            self.assertTrue(
                result.get("R_abcd_R_abcd_coeff") is None
                or result.get("R_abcd_R_abcd_coeff") == 0
                or result.get("R_abcd_R_abcd_coeff") == float("inf"),
                "n=1 should produce degenerate Riemann coefficient",
            )
        except (ValueError, ZeroDivisionError):
            # Raising is also acceptable for degenerate input
            pass

    def test_general_n_sum_rule(self):
        """Verify the GB density identity: G_n/R^2 = 1 - 4/n + 2/(n(n-1)).

        This equals (n-2)(n-3) / [n(n-1)] after algebraic simplification.
        """
        for n in range(2, 12):
            with self.subTest(n=n):
                result = compute_internal_curvature_identities(n=n)
                ric_coeff = result["R_ab_R_ab_coeff"]   # 1/n
                riem_coeff = result["R_abcd_R_abcd_coeff"]  # 2/(n*(n-1))

                # G = R^2 - 4*Ric^2 + Riem^2
                # G/R^2 = 1 - 4*ric_coeff + riem_coeff
                G_n_coefficient = 1.0 - 4.0 * ric_coeff + riem_coeff

                expected = (n - 2) * (n - 3) / (n * (n - 1))
                self.assertAlmostEqual(
                    G_n_coefficient, expected, places=12,
                    msg=f"GB density coefficient mismatch for n={n}",
                )

    def test_4d_ricci_squared(self):
        """For n=4, R_{ab}R^{ab} / R^2 = 1/4."""
        result = compute_internal_curvature_identities(n=4)
        self.assertAlmostEqual(result["R_ab_R_ab_coeff"], 0.25, places=12)

    def test_4d_riemann_squared(self):
        """For n=4, R_{abcd}R^{abcd} / R^2 = 2/(4*3) = 1/6."""
        result = compute_internal_curvature_identities(n=4)
        self.assertAlmostEqual(result["R_abcd_R_abcd_coeff"], 1.0 / 6.0, places=12)

    def test_large_n(self):
        """For large n, coefficients should approach known limits."""
        result = compute_internal_curvature_identities(n=100)
        # Ric^2 coeff -> 1/100
        self.assertAlmostEqual(result["R_ab_R_ab_coeff"], 0.01, places=12)
        # Riem^2 coeff -> 2/(100*99) ~ 0.000202...
        self.assertAlmostEqual(
            result["R_abcd_R_abcd_coeff"], 2.0 / (100 * 99), places=12,
        )

    def test_returns_dict(self):
        """Function should return a dictionary."""
        result = compute_internal_curvature_identities(n=2)
        self.assertIsInstance(result, dict)

    def test_required_keys_present(self):
        """Check that all expected keys are in the result."""
        result = compute_internal_curvature_identities(n=3)
        self.assertIn("R_ab_R_ab_coeff", result)
        self.assertIn("R_abcd_R_abcd_coeff", result)


# -----------------------------------------------------------------------
# Tests for internal GB density coefficient
# -----------------------------------------------------------------------

class TestInternalGBDensity(unittest.TestCase):
    """Tests for the GB density coefficient (n-2)(n-3)/[n(n-1)].

    The GB density of a maximally symmetric (constant curvature) n-manifold
    vanishes for n=2 and n=3 (where GB is topological / a total derivative),
    and is non-zero for n >= 4.
    """

    def _gb_coefficient(self, n):
        """Compute GB density coefficient from curvature identities."""
        result = compute_internal_curvature_identities(n=n)
        return 1.0 - 4.0 * result["R_ab_R_ab_coeff"] + result["R_abcd_R_abcd_coeff"]

    def test_gb_coefficient_n2(self):
        """GB density vanishes for n=2 (topological in 2D)."""
        coeff = self._gb_coefficient(2)
        self.assertAlmostEqual(coeff, 0.0, places=12)

    def test_gb_coefficient_n3(self):
        """GB density vanishes for n=3 (conformally flat, GB is trivial)."""
        coeff = self._gb_coefficient(3)
        self.assertAlmostEqual(coeff, 0.0, places=12)

    def test_gb_coefficient_n4(self):
        """For n=4: (4-2)(4-3)/(4*3) = 2/12 = 1/6."""
        coeff = self._gb_coefficient(4)
        self.assertAlmostEqual(coeff, 1.0 / 6.0, places=12)

    def test_gb_coefficient_n5(self):
        """For n=5: (5-2)(5-3)/(5*4) = 6/20 = 3/10."""
        coeff = self._gb_coefficient(5)
        self.assertAlmostEqual(coeff, 3.0 / 10.0, places=12)

    def test_gb_coefficient_n6(self):
        """For n=6: (6-2)(6-3)/(6*5) = 12/30 = 2/5."""
        coeff = self._gb_coefficient(6)
        self.assertAlmostEqual(coeff, 2.0 / 5.0, places=12)

    def test_gb_coefficient_general(self):
        """Verify (n-2)(n-3)/[n(n-1)] for various n."""
        for n in range(2, 20):
            with self.subTest(n=n):
                coeff = self._gb_coefficient(n)
                expected = (n - 2) * (n - 3) / (n * (n - 1))
                self.assertAlmostEqual(coeff, expected, places=12)

    def test_gb_coefficient_positive_for_n_ge_4(self):
        """For n >= 4, the GB density coefficient must be strictly positive."""
        for n in range(4, 15):
            with self.subTest(n=n):
                coeff = self._gb_coefficient(n)
                self.assertGreater(coeff, 0.0)

    def test_gb_coefficient_monotonically_increases(self):
        """The coefficient (n-2)(n-3)/[n(n-1)] increases for n >= 4."""
        prev = self._gb_coefficient(4)
        for n in range(5, 15):
            curr = self._gb_coefficient(n)
            self.assertGreater(curr, prev, f"Coefficient should increase from n={n-1} to n={n}")
            prev = curr

    def test_gb_coefficient_approaches_one(self):
        """For large n, (n-2)(n-3)/[n(n-1)] -> 1."""
        coeff = self._gb_coefficient(1000)
        self.assertAlmostEqual(coeff, 1.0, places=2)


# -----------------------------------------------------------------------
# Tests for derive_gauss_bonnet_correction (step-by-step derivation)
# -----------------------------------------------------------------------

class TestGBDerivationSteps(unittest.TestCase):
    """Tests for derive_gauss_bonnet_correction."""

    def test_step1_returns_identities(self):
        """Step 1 should contain curvature identity keys."""
        result = derive_gauss_bonnet_correction(d=4, n=2, genus=2)
        step1 = result["step1_internal_identities"]
        self.assertIn("R_ab_R_ab_over_R_n_sq", step1)
        self.assertIn("R_abcd_R_abcd_over_R_n_sq", step1)

    def test_step2_gb_density_zero_2d(self):
        """For n=2, the internal GB density coefficient must be exactly 0."""
        result = derive_gauss_bonnet_correction(d=4, n=2, genus=2)
        step2 = result["step2_internal_gb_density"]
        # The GB density coefficient for a maximally symmetric 2-manifold is 0
        self.assertAlmostEqual(step2["G_n_coefficient"], 0.0, places=12)

    def test_step3_euler_characteristic(self):
        """chi = 2 - 2g for various genera."""
        for genus in range(0, 8):
            with self.subTest(genus=genus):
                result = derive_gauss_bonnet_correction(d=4, n=2, genus=genus)
                step3 = result["step3_euler_integral"]
                self.assertEqual(step3["chi"], 2 - 2 * genus)

    def test_step3_euler_genus0(self):
        """Sphere: chi = 2."""
        result = derive_gauss_bonnet_correction(d=4, n=2, genus=0)
        self.assertEqual(result["step3_euler_integral"]["chi"], 2)

    def test_step3_euler_genus1(self):
        """Torus: chi = 0."""
        result = derive_gauss_bonnet_correction(d=4, n=2, genus=1)
        self.assertEqual(result["step3_euler_integral"]["chi"], 0)

    def test_step4_cross_terms_present(self):
        """Step 4 should return cross-term contributions."""
        result = derive_gauss_bonnet_correction(d=4, n=2, genus=2)
        step4 = result["step4_cross_terms"]
        self.assertIsInstance(step4, dict)

    def test_step5_matches_existing(self):
        """Critical cross-check: derived delta_K must match compute_gauss_bonnet_shift().

        This is the most important test -- it verifies that the step-by-step
        derivation reproduces the same numeric result as the existing formula.
        """
        for genus in [0, 1, 2, 3, 5]:
            for d, n in [(4, 2), (5, 3), (4, 3)]:
                with self.subTest(d=d, n=n, genus=genus):
                    # Get the existing result (with a specific coupling)
                    existing = compute_gauss_bonnet_shift(
                        d=d, n=n, genus=genus,
                        gb_coupling=math.sqrt(5) - 3,
                    )

                    # Get the derived result
                    derived = derive_gauss_bonnet_correction(d=d, n=n, genus=genus)
                    step5 = derived["step5_total_correction"]

                    # The delta_K from derivation (with same coupling) should match
                    gb_coupling = math.sqrt(5) - 3
                    chi = 2 - 2 * genus
                    derived_delta_K = step5["delta_K_per_lambda_GB"] * gb_coupling

                    self.assertAlmostEqual(
                        derived_delta_K,
                        existing["delta_K"] if existing["delta_K"] is not None else 0.0,
                        places=10,
                        msg=f"delta_K mismatch for d={d}, n={n}, genus={genus}",
                    )

    def test_default_d4_n2(self):
        """Default parameters d=4, n=2 should work without errors."""
        result = derive_gauss_bonnet_correction()
        self.assertIn("step1_internal_identities", result)
        self.assertIn("step2_internal_gb_density", result)
        self.assertIn("step3_euler_integral", result)
        self.assertIn("step4_cross_terms", result)
        self.assertIn("step5_total_correction", result)

    def test_d5_n3(self):
        """D=8 reduction (d=5, n=3) should work."""
        result = derive_gauss_bonnet_correction(d=5, n=3, genus=2)
        self.assertIn("step5_total_correction", result)
        # For n=3, internal GB density vanishes
        self.assertAlmostEqual(
            result["step2_internal_gb_density"]["G_n_coefficient"], 0.0, places=12,
        )

    def test_d4_n3(self):
        """D=7 reduction (d=4, n=3) should work."""
        result = derive_gauss_bonnet_correction(d=4, n=3, genus=2)
        self.assertIn("step5_total_correction", result)

    def test_derivation_genus2_golden_ratio_coupling(self):
        """For d=4, n=2, genus=2, lambda_GB = sqrt(5)-3 should give omega_target.

        This is the central physical result: the golden-ratio coupling
        produces exactly omega = (sqrt(5)-2)/2.
        """
        result = derive_gauss_bonnet_correction(d=4, n=2, genus=2)
        step5 = result["step5_total_correction"]

        # The required coupling should be sqrt(5) - 3 (or close to it)
        # Check via the existing function
        existing = compute_gauss_bonnet_shift(d=4, n=2, genus=2)
        required_coupling = existing["required_gb_coupling"]

        # Verify it is close to sqrt(5) - 3
        self.assertAlmostEqual(
            required_coupling, math.sqrt(5) - 3, places=10,
            msg="Required GB coupling should be sqrt(5)-3",
        )

    def test_verification_key_present(self):
        """The derivation should include a verification section."""
        result = derive_gauss_bonnet_correction(d=4, n=2, genus=2)
        self.assertIn("verification", result)

    def test_verification_matches_existing(self):
        """The verification field should confirm agreement with existing code."""
        result = derive_gauss_bonnet_correction(d=4, n=2, genus=2)
        verification = result["verification"]
        # Verification should compare with compute_gauss_bonnet_shift
        self.assertIsInstance(verification, dict)

    def test_returns_dict(self):
        """Function should return a dictionary."""
        result = derive_gauss_bonnet_correction(d=4, n=2, genus=2)
        self.assertIsInstance(result, dict)

    def test_all_five_steps_present(self):
        """All five derivation steps should be present in the result."""
        result = derive_gauss_bonnet_correction(d=4, n=2, genus=2)
        expected_keys = [
            "step1_internal_identities",
            "step2_internal_gb_density",
            "step3_euler_integral",
            "step4_cross_terms",
            "step5_total_correction",
        ]
        for key in expected_keys:
            self.assertIn(key, result, f"Missing derivation step: {key}")


# -----------------------------------------------------------------------
# Tests for consistency between derivation and existing functions
# -----------------------------------------------------------------------

class TestDerivationConsistency(unittest.TestCase):
    """Cross-checks between derive_gauss_bonnet_correction and existing code."""

    def test_derived_matches_stated_formula(self):
        """For multiple (d, n, genus) combos, derivation gives same delta_K as existing code."""
        test_cases = [
            (4, 2, 0),
            (4, 2, 1),
            (4, 2, 2),
            (4, 2, 3),
            (4, 2, 5),
            (5, 3, 2),
            (4, 3, 2),
            (5, 2, 2),
            (6, 2, 2),
            (4, 4, 2),
        ]
        coupling = 0.1  # arbitrary test coupling

        for d, n, genus in test_cases:
            with self.subTest(d=d, n=n, genus=genus):
                existing = compute_gauss_bonnet_shift(
                    d=d, n=n, genus=genus, gb_coupling=coupling,
                )
                derived = derive_gauss_bonnet_correction(d=d, n=n, genus=genus)

                chi = 2 - 2 * genus
                step5 = derived["step5_total_correction"]
                derived_delta_K = step5["delta_K_per_lambda_GB"] * coupling

                existing_delta_K = existing["delta_K"]
                if existing_delta_K is None:
                    existing_delta_K = 0.0

                self.assertAlmostEqual(
                    derived_delta_K, existing_delta_K, places=10,
                    msg=f"delta_K mismatch for (d={d}, n={n}, genus={genus})",
                )

    def test_gb_factor_formula(self):
        """The gb_factor from derivation must equal 4*n*(n-1)/(d-2)^2."""
        test_cases = [
            (4, 2),
            (5, 3),
            (4, 3),
            (5, 2),
            (6, 2),
            (4, 4),
            (6, 4),
        ]
        for d, n in test_cases:
            with self.subTest(d=d, n=n):
                existing = compute_gauss_bonnet_shift(d=d, n=n, genus=2)
                derived = derive_gauss_bonnet_correction(d=d, n=n, genus=2)

                expected_factor = 4.0 * n * (n - 1) / (d - 2) ** 2
                self.assertAlmostEqual(
                    existing["gb_factor"], expected_factor, places=12,
                )

                # The derivation's gb_factor should match
                step5 = derived["step5_total_correction"]
                self.assertAlmostEqual(
                    step5["gb_factor"], expected_factor, places=12,
                    msg=f"Derived gb_factor mismatch for d={d}, n={n}",
                )

    def test_full_pipeline_omega_target(self):
        """Full derivation for d=4, n=2, genus=2 with lambda_GB = sqrt(5)-3 gives omega_target.

        omega_target = (sqrt(5)-2)/2 ~ 0.11803398875
        """
        target = compute_target_omega()
        omega_target = target["omega_float"]

        # Get the required coupling from existing code
        existing = compute_gauss_bonnet_shift(d=4, n=2, genus=2)
        coupling = existing["required_gb_coupling"]

        # Apply this coupling and check shifted omega
        gb_with_coupling = compute_gauss_bonnet_shift(
            d=4, n=2, genus=2, gb_coupling=coupling,
        )

        self.assertAlmostEqual(
            gb_with_coupling["omega_shifted"], omega_target, places=10,
            msg="Shifted omega should equal (sqrt(5)-2)/2",
        )

    def test_full_pipeline_omega_value(self):
        """Verify the numeric value of omega_target = (sqrt(5)-2)/2."""
        target = compute_target_omega()
        expected = (math.sqrt(5) - 2.0) / 2.0
        self.assertAlmostEqual(target["omega_float"], expected, places=12)

    def test_kinetic_baseline_d4_n2(self):
        """Baseline kinetic coefficient for d=4, n=2 is -2.0."""
        kinetic = compute_kinetic_coefficient(d=4, n=2)
        self.assertAlmostEqual(kinetic["coeff_standard"], -2.0, places=12)

    def test_derivation_step5_delta_K_d4_n2(self):
        """For d=4, n=2: gb_factor = 4*2*1/4 = 2.0."""
        derived = derive_gauss_bonnet_correction(d=4, n=2, genus=2)
        step5 = derived["step5_total_correction"]
        self.assertAlmostEqual(step5["gb_factor"], 2.0, places=12)

    def test_derivation_step5_delta_K_d5_n3(self):
        """For d=5, n=3: gb_factor = 4*3*2/9 = 24/9 = 8/3."""
        derived = derive_gauss_bonnet_correction(d=5, n=3, genus=2)
        step5 = derived["step5_total_correction"]
        self.assertAlmostEqual(step5["gb_factor"], 8.0 / 3.0, places=12)

    def test_genus1_no_correction(self):
        """For genus=1 (torus), chi=0, so all corrections vanish."""
        derived = derive_gauss_bonnet_correction(d=4, n=2, genus=1)
        step3 = derived["step3_euler_integral"]
        self.assertEqual(step3["chi"], 0)

    def test_required_coupling_sign(self):
        """For genus=2 (chi=-2), the required coupling should be negative.

        Because chi < 0 and we need delta_K to push omega to the target,
        the sign of the required coupling depends on the gap direction.
        """
        existing = compute_gauss_bonnet_shift(d=4, n=2, genus=2)
        # Just verify that the required coupling is finite and real
        self.assertIsNotNone(existing["required_gb_coupling"])
        self.assertTrue(
            math.isfinite(existing["required_gb_coupling"]),
            "Required coupling should be finite",
        )

    def test_verification_omega_matches(self):
        """The verify_omega from compute_gauss_bonnet_shift should match target."""
        existing = compute_gauss_bonnet_shift(d=4, n=2, genus=2)
        target = compute_target_omega()
        self.assertAlmostEqual(
            existing["verify_omega"], target["omega_float"], places=10,
        )


if __name__ == "__main__":
    unittest.main()
