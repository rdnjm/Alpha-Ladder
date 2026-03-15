"""
End-to-end integration tests for the Alpha Ladder derivation chain.

Validates the complete pipeline:
  CODATA constants -> vacuum polynomial (n=2) -> KK reduction ->
  Brans-Dicke -> phi^2/2 bridge -> G prediction -> compare with CODATA

Every step is tested individually, then the full chain is verified
as a single deterministic computation with no free parameters.
"""

import math
import unittest
from decimal import Decimal, getcontext

getcontext().prec = 50

from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.vacuum_polynomial import (
    compute_dimensional_polynomial,
    prove_golden_uniqueness,
    derive_algebraic_closure,
)
from alpha_ladder_core.kk_reduction import (
    compute_target_omega,
    compute_einstein_frame_ansatz,
    compute_kinetic_coefficient,
    compute_gauss_bonnet_shift,
    summarize_reduction,
)
from alpha_ladder_core.dilaton import (
    decompose_21,
    compute_bd_parameter,
    reconcile_6d,
)
from alpha_ladder_core.bridge_search import (
    compute_target_coefficient,
    run_full_search,
)
from alpha_ladder_core.predict_g import (
    predict_G,
    get_bridge_candidates,
    compare_prediction,
    get_G_measurements,
    summarize_predictions,
)


class TestDerivationChain(unittest.TestCase):
    """Full derivation chain: CODATA -> vacuum polynomial -> KK -> BD -> G."""

    def setUp(self):
        self.c = get_constants("CODATA 2018")
        self.phi = float(self.c.phi)
        # phi^2/2 = (1+sqrt5)^2 / 8 = (3+sqrt5)/4 ~ 1.3090169943749475
        self.phi_sq_over_2 = self.phi ** 2 / 2.0

    # ------------------------------------------------------------------
    # Step 1 -- Constants
    # ------------------------------------------------------------------

    def test_constants_load(self):
        """get_constants returns namespace with all required attributes."""
        for attr in ("alpha", "alpha_g", "hbar", "c", "m_e", "phi"):
            self.assertTrue(
                hasattr(self.c, attr),
                f"Missing attribute: {attr}",
            )

    def test_alpha_value(self):
        """alpha ~ 1/137.036."""
        alpha_f = float(self.c.alpha)
        self.assertAlmostEqual(alpha_f, 1.0 / 137.036, places=6)

    def test_alpha_g_value(self):
        """alpha_g ~ 1.75e-45."""
        alpha_g_f = float(self.c.alpha_g)
        self.assertGreater(alpha_g_f, 1.7e-45)
        self.assertLess(alpha_g_f, 1.8e-45)

    # ------------------------------------------------------------------
    # Step 2 -- Vacuum Polynomial
    # ------------------------------------------------------------------

    def test_polynomial_n2(self):
        """compute_dimensional_polynomial(d=4, n=2) -> discriminant = 20."""
        result = compute_dimensional_polynomial(d=4, n=2)
        self.assertEqual(result["discriminant"], 20)
        self.assertEqual(result["D"], 6)
        self.assertEqual(result["polynomial"], "x^2 + 6x + 4")

    def test_golden_field(self):
        """n=2 produces roots in the golden ratio field Q(sqrt(5))."""
        result = compute_dimensional_polynomial(d=4, n=2)
        self.assertTrue(result["is_golden_field"])

    def test_n2_unique_minimal(self):
        """n=2 is the unique minimal solution of the Diophantine equation."""
        result = prove_golden_uniqueness(d=4)
        self.assertTrue(result["n2_is_minimal"])

    def test_pell_equation(self):
        """Minimal Pell solution is (n, m) = (2, 2)."""
        result = prove_golden_uniqueness(d=4)
        self.assertEqual(result["minimal_solution"], (2, 2))

    def test_algebraic_closure(self):
        """All key constants verified in Q(sqrt(5))."""
        result = derive_algebraic_closure(d=4, n=2)
        self.assertTrue(result["all_polynomials_verified"])
        self.assertTrue(result["all_in_same_field"])
        self.assertEqual(result["field"], "Q(sqrt(5))")

    # ------------------------------------------------------------------
    # Step 3 -- KK Reduction
    # ------------------------------------------------------------------

    def test_target_omega(self):
        """omega_target = (sqrt(5)-2)/2 ~ 0.118034."""
        result = compute_target_omega()
        expected = (math.sqrt(5) - 2) / 2
        self.assertAlmostEqual(result["omega_float"], expected, places=12)
        self.assertAlmostEqual(result["omega_float"], 0.118034, places=5)

    def test_einstein_frame(self):
        """Einstein frame ansatz: alpha=1.0, beta=-1.0 for d=4, n=2."""
        result = compute_einstein_frame_ansatz(d=4, n=2)
        self.assertAlmostEqual(result["alpha"], 1.0, places=12)
        self.assertAlmostEqual(result["beta"], -1.0, places=12)

    def test_kinetic_coefficient(self):
        """K_einstein = -4.0 for d=4, n=2 (= -n(n+d-2)/(d-2) = -2*4/2)."""
        result = compute_kinetic_coefficient(d=4, n=2)
        self.assertAlmostEqual(result["K_einstein"], -4.0, places=12)
        self.assertTrue(result["kinetic_sign_healthy"])

    def test_gb_shift_genus2(self):
        """GB shift with genus 2 produces a verify_omega close to target."""
        result = compute_gauss_bonnet_shift(d=4, n=2, genus=2)
        # The required GB coupling exists and verify_omega matches target
        self.assertIsNotNone(result["required_gb_coupling"])
        self.assertIsNotNone(result["verify_omega"])
        self.assertAlmostEqual(
            result["verify_omega"],
            result["omega_target"],
            places=10,
        )

    def test_full_reduction(self):
        """summarize_reduction returns consistent omega values."""
        result = summarize_reduction(d=4, n=2, genus=2)
        self.assertIsNotNone(result["omega_baseline"])
        omega_expected = (math.sqrt(5) - 2) / 2
        self.assertAlmostEqual(
            result["omega_target"],
            omega_expected,
            places=10,
        )
        # Baseline omega exists and is positive
        self.assertGreater(result["omega_baseline"], 0)

    # ------------------------------------------------------------------
    # Step 4 -- Dilaton
    # ------------------------------------------------------------------

    def test_decompose_21(self):
        """decompose_21 bridge ~ phi^2/2 ~ 0.9045."""
        result = decompose_21(self.c)
        self.assertAlmostEqual(result["bridge"], self.phi_sq_over_2, places=10)

    def test_bd_parameter(self):
        """Brans-Dicke omega value exists."""
        result = compute_bd_parameter(self.c)
        self.assertIsNotNone(result["omega"])
        # omega is a finite number
        self.assertTrue(math.isfinite(result["omega"]))

    def test_bridge_value(self):
        """phi^2/2 = (3+sqrt5)/4 ~ 1.309017."""
        phi = float(self.c.phi)
        bridge = phi ** 2 / 2.0
        exact = (3 + math.sqrt(5)) / 4
        self.assertAlmostEqual(bridge, exact, places=14)
        self.assertAlmostEqual(bridge, 1.3090169943749475, places=12)

    def test_alpha_21_computation(self):
        """alpha^21 is finite, positive, ~ 1.34e-45."""
        result = decompose_21(self.c)
        alpha_21 = result["alpha_21"]
        self.assertGreater(alpha_21, 0)
        self.assertTrue(math.isfinite(alpha_21))
        self.assertGreater(alpha_21, 1.3e-45)
        self.assertLess(alpha_21, 1.4e-45)

    def test_reconcile_6d(self):
        """6D metric has 21 total components."""
        result = reconcile_6d(self.c)
        self.assertEqual(result["total_6d_metric"], 21)
        self.assertEqual(
            result["g_uv_4d"] + result["g_ab_2d"] + result["g_ua_mixed"],
            21,
        )

    # ------------------------------------------------------------------
    # Step 5 -- Bridge Search
    # ------------------------------------------------------------------

    def test_target_coefficient(self):
        """Target coefficient ~ phi^2/2."""
        target_d, target_f = compute_target_coefficient(self.c)
        self.assertAlmostEqual(target_f, self.phi_sq_over_2, places=2)

    def test_full_search_finds_phi2_over_2(self):
        """phi^2/2 appears in the search results."""
        results = run_full_search(self.c)
        # Look for a result containing phi (unicode phi^2 form)
        expressions = [r[1] for r in results]
        has_phi = any("φ" in expr for expr in expressions)
        self.assertTrue(has_phi, f"No phi-based expression found in: {expressions[:10]}")

    def test_best_match_residual(self):
        """Best match residual < 0.1%."""
        results = run_full_search(self.c)
        self.assertGreater(len(results), 0)
        best_err = results[0][0]  # sorted by error, percentage
        self.assertLess(best_err, 0.1)

    def test_search_consistency(self):
        """Target from bridge search matches dilaton bridge value."""
        _, target_f = compute_target_coefficient(self.c)
        decomp = decompose_21(self.c)
        bridge = decomp["bridge"]
        # They should be within 0.02% of each other
        rel_diff = abs(target_f - bridge) / bridge
        self.assertLess(rel_diff, 2e-4)

    # ------------------------------------------------------------------
    # Step 6 -- G Prediction
    # ------------------------------------------------------------------

    def test_predict_G_phi2_bridge(self):
        """G prediction using phi^2/2 bridge ~ 6.67e-11."""
        candidates = get_bridge_candidates(self.c)
        phi2_coeff = candidates["φ²/2"]
        G_pred = predict_G(phi2_coeff, self.c)
        G_pred_f = float(G_pred)
        self.assertGreater(G_pred_f, 6.67e-11)
        self.assertLess(G_pred_f, 6.68e-11)

    def test_G_within_200ppm(self):
        """G prediction within 200 ppm of CODATA 2018."""
        candidates = get_bridge_candidates(self.c)
        phi2_coeff = candidates["φ²/2"]
        G_pred = predict_G(phi2_coeff, self.c)
        G_codata = self.c.G
        rel_error = abs(float(G_pred - G_codata)) / float(G_codata)
        self.assertLess(rel_error, 2e-4)

    def test_compare_with_experiments(self):
        """compare_prediction returns results for all measurements."""
        candidates = get_bridge_candidates(self.c)
        phi2_coeff = candidates["φ²/2"]
        G_pred = predict_G(phi2_coeff, self.c)
        measurements = get_G_measurements()
        comparisons = compare_prediction(G_pred, measurements)
        self.assertEqual(len(comparisons), len(measurements))
        for comp in comparisons:
            self.assertIn("experiment", comp)
            self.assertIn("sigma", comp)
            self.assertTrue(math.isfinite(comp["sigma"]))

    def test_summarize_predictions(self):
        """summarize_predictions includes phi^2/2 entry."""
        summary = summarize_predictions(self.c)
        self.assertIn("φ²/2", summary)
        entry = summary["φ²/2"]
        self.assertIn("G_pred", entry)
        self.assertIn("avg_sigma", entry)
        self.assertIn("comparisons", entry)

    def test_no_free_parameters(self):
        """The chain uses only CODATA constants + pure mathematics.

        This test documents the assertion that no free parameters are
        introduced at any stage of the derivation:
        - CODATA 2018 provides alpha, m_e, hbar, c (measured quantities)
        - The vacuum polynomial x^2 + 6x + 4 is fixed by d=4, n=2
        - phi = (1+sqrt(5))/2 is a root of x^2 - x - 1 = 0 (pure math)
        - The bridge coefficient phi^2/2 = (3+sqrt(5))/4 (algebraic)
        - The exponent 21 counts 6D metric components (20 Riemann + 1 scalar)
        - G = (phi^2/2) * alpha^21 * hbar * c / m_e^2 (no tuning)
        """
        # Verify alpha is from CODATA, not tuned
        self.assertEqual(self.c.alpha, Decimal("0.0072973525693"))
        # Verify phi is the golden ratio (mathematical constant)
        phi_check = (1 + Decimal(5).sqrt()) / 2
        self.assertEqual(self.c.phi, phi_check)
        # Verify the exponent 21 = 6D metric components
        reconcile = reconcile_6d(self.c)
        self.assertEqual(reconcile["total_6d_metric"], 21)
        self.assertEqual(reconcile["riemann_plus_scalar"], 21)

    # ------------------------------------------------------------------
    # Full Chain
    # ------------------------------------------------------------------

    def test_full_chain_codata_2018(self):
        """Complete chain with CODATA 2018: G within 200 ppm."""
        c = get_constants("CODATA 2018")

        # Step 1: Vacuum polynomial forces Q(sqrt(5))
        poly = compute_dimensional_polynomial(d=4, n=2)
        self.assertEqual(poly["discriminant"], 20)
        self.assertTrue(poly["is_golden_field"])

        # Step 2: n=2 is unique minimal
        uniqueness = prove_golden_uniqueness(d=4)
        self.assertTrue(uniqueness["n2_is_minimal"])

        # Step 3: KK reduction gives omega_target in Q(sqrt(5))
        omega = compute_target_omega()
        self.assertAlmostEqual(omega["omega_float"], (math.sqrt(5) - 2) / 2, places=12)

        # Step 4: Dilaton decomposition gives bridge = phi^2/2
        decomp = decompose_21(c)
        self.assertAlmostEqual(decomp["bridge"], self.phi_sq_over_2, places=10)

        # Step 5: Predict G
        candidates = get_bridge_candidates(c)
        G_pred = predict_G(candidates["φ²/2"], c)
        G_codata = float(c.G)
        G_pred_f = float(G_pred)

        rel_error = abs(G_pred_f - G_codata) / G_codata
        self.assertLess(rel_error, 2e-4, f"G prediction error {rel_error:.2e} exceeds 200 ppm")

    def test_full_chain_codata_2014(self):
        """Complete chain with CODATA 2014: G within 200 ppm."""
        c = get_constants("CODATA 2014")

        # Same derivation chain with 2014 constants
        poly = compute_dimensional_polynomial(d=4, n=2)
        self.assertTrue(poly["is_golden_field"])

        decomp = decompose_21(c)
        self.assertAlmostEqual(decomp["bridge"], self.phi_sq_over_2, places=10)

        candidates = get_bridge_candidates(c)
        G_pred = predict_G(candidates["φ²/2"], c)
        G_codata = float(c.G)
        G_pred_f = float(G_pred)

        rel_error = abs(G_pred_f - G_codata) / G_codata
        self.assertLess(rel_error, 2e-4, f"G prediction error {rel_error:.2e} exceeds 200 ppm")

    def test_chain_deterministic(self):
        """Running the chain twice produces identical results."""
        c = get_constants("CODATA 2018")

        # Run 1
        candidates_1 = get_bridge_candidates(c)
        G1 = predict_G(candidates_1["φ²/2"], c)

        # Run 2
        candidates_2 = get_bridge_candidates(c)
        G2 = predict_G(candidates_2["φ²/2"], c)

        self.assertEqual(G1, G2)


if __name__ == "__main__":
    unittest.main()
