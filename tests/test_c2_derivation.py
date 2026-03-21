"""Tests for c2_derivation module.

Systematic derivation attempts for the coefficient c_2 = 3 appearing in the
corrected bridge formula phi^2/2 * (1 + 3*alpha^2). Tests cover heat kernel,
Euler characteristic, conformal anomaly, Seeley-DeWitt coefficients, spectral
zeta values, Killing vector analysis, sphere scan, and loop attempts.
"""
import unittest
import math


# ---------------------------------------------------------------------------
# 1. compute_heat_kernel_a1
# ---------------------------------------------------------------------------
class TestHeatKernelA1(unittest.TestCase):
    """Tests for compute_heat_kernel_a1."""

    def test_returns_dict(self):
        from alpha_ladder_core.c2_derivation import compute_heat_kernel_a1
        result = compute_heat_kernel_a1(n=2)
        self.assertIsInstance(result, dict)

    def test_s2_a1_density(self):
        from alpha_ladder_core.c2_derivation import compute_heat_kernel_a1
        result = compute_heat_kernel_a1(n=2)
        self.assertAlmostEqual(result["a1_density"], 1.0 / 3.0, places=10)

    def test_s2_inverse_a1(self):
        from alpha_ladder_core.c2_derivation import compute_heat_kernel_a1
        result = compute_heat_kernel_a1(n=2)
        self.assertEqual(result["inverse_a1"], 3)

    def test_s2_matches_c2(self):
        from alpha_ladder_core.c2_derivation import compute_heat_kernel_a1
        result = compute_heat_kernel_a1(n=2)
        self.assertTrue(result["matches_c2"])

    def test_s3_a1_density(self):
        from alpha_ladder_core.c2_derivation import compute_heat_kernel_a1
        result = compute_heat_kernel_a1(n=3)
        self.assertAlmostEqual(result["a1_density"], 1.0, places=10)

    def test_s3_inverse_a1(self):
        from alpha_ladder_core.c2_derivation import compute_heat_kernel_a1
        result = compute_heat_kernel_a1(n=3)
        self.assertAlmostEqual(result["inverse_a1"], 1.0, places=10)

    def test_s3_does_not_match_c2(self):
        from alpha_ladder_core.c2_derivation import compute_heat_kernel_a1
        result = compute_heat_kernel_a1(n=3)
        self.assertFalse(result["matches_c2"])

    def test_s1_handles_gracefully(self):
        from alpha_ladder_core.c2_derivation import compute_heat_kernel_a1
        result = compute_heat_kernel_a1(n=1)
        self.assertIsInstance(result, dict)
        # S^1 has R=0, so a1=0, inverse undefined
        self.assertAlmostEqual(result["a1_density"], 0.0, places=10)

    def test_s2_ricci_scalar(self):
        from alpha_ladder_core.c2_derivation import compute_heat_kernel_a1
        result = compute_heat_kernel_a1(n=2)
        self.assertAlmostEqual(result["R_Sn"], 2.0, places=10)

    def test_s2_volume(self):
        from alpha_ladder_core.c2_derivation import compute_heat_kernel_a1
        result = compute_heat_kernel_a1(n=2)
        self.assertAlmostEqual(result["Vol_Sn"], 4.0 * math.pi, delta=1e-10)


# ---------------------------------------------------------------------------
# 2. compute_chi_tangent_bundle
# ---------------------------------------------------------------------------
class TestChiTangentBundle(unittest.TestCase):
    """Tests for compute_chi_tangent_bundle."""

    def test_returns_dict(self):
        from alpha_ladder_core.c2_derivation import compute_chi_tangent_bundle
        result = compute_chi_tangent_bundle(n=2)
        self.assertIsInstance(result, dict)

    def test_s2_chi_t(self):
        from alpha_ladder_core.c2_derivation import compute_chi_tangent_bundle
        result = compute_chi_tangent_bundle(n=2)
        self.assertEqual(result["chi_T"], 3)

    def test_s2_h0(self):
        from alpha_ladder_core.c2_derivation import compute_chi_tangent_bundle
        result = compute_chi_tangent_bundle(n=2)
        self.assertEqual(result["h0"], 3)

    def test_s2_h1(self):
        from alpha_ladder_core.c2_derivation import compute_chi_tangent_bundle
        result = compute_chi_tangent_bundle(n=2)
        self.assertEqual(result["h1"], 0)

    def test_s2_is_complex(self):
        from alpha_ladder_core.c2_derivation import compute_chi_tangent_bundle
        result = compute_chi_tangent_bundle(n=2)
        self.assertTrue(result["is_complex_manifold"])

    def test_s3_not_complex(self):
        from alpha_ladder_core.c2_derivation import compute_chi_tangent_bundle
        result = compute_chi_tangent_bundle(n=3)
        self.assertFalse(result["is_complex_manifold"])

    def test_s2_hrr_decomposition(self):
        from alpha_ladder_core.c2_derivation import compute_chi_tangent_bundle
        result = compute_chi_tangent_bundle(n=2)
        self.assertIn("2w", result["hrr_decomposition"]["euler_part"])
        self.assertIn("Td", result["hrr_decomposition"]["todd_part"])

    def test_s2_killing_vectors(self):
        from alpha_ladder_core.c2_derivation import compute_chi_tangent_bundle
        result = compute_chi_tangent_bundle(n=2)
        self.assertEqual(result["killing_vectors"], 3)


# ---------------------------------------------------------------------------
# 3. compute_gauge_matching_identity
# ---------------------------------------------------------------------------
class TestGaugeMatchingIdentity(unittest.TestCase):
    """Tests for compute_gauge_matching_identity."""

    def test_returns_dict(self):
        from alpha_ladder_core.c2_derivation import compute_gauge_matching_identity
        result = compute_gauge_matching_identity()
        self.assertIsInstance(result, dict)

    def test_g_kk_squared_value(self):
        from alpha_ladder_core.c2_derivation import compute_gauge_matching_identity
        result = compute_gauge_matching_identity()
        alpha = 7.2973525693e-3
        expected = 4.0 * math.pi * alpha
        self.assertAlmostEqual(float(result["g_KK_squared"]), expected, places=6)

    def test_ratio_matches_alpha_squared(self):
        from alpha_ladder_core.c2_derivation import compute_gauge_matching_identity
        result = compute_gauge_matching_identity()
        self.assertLess(abs(result["residual_ppm"]), 0.01)

    def test_pi_cancellation(self):
        from alpha_ladder_core.c2_derivation import compute_gauge_matching_identity
        result = compute_gauge_matching_identity()
        self.assertTrue(result["pi_cancellation"])

    def test_residual_small(self):
        from alpha_ladder_core.c2_derivation import compute_gauge_matching_identity
        result = compute_gauge_matching_identity()
        self.assertLess(result["residual_ppm"], 0.01)

    def test_with_constants_parameter(self):
        from alpha_ladder_core.c2_derivation import compute_gauge_matching_identity
        from alpha_ladder_core.constants import get_constants
        c = get_constants("CODATA 2018")
        result = compute_gauge_matching_identity(constants=c)
        self.assertIsInstance(result, dict)
        self.assertIn("g_KK_squared", result)


# ---------------------------------------------------------------------------
# 4. compute_conformal_anomaly
# ---------------------------------------------------------------------------
class TestConformalAnomaly(unittest.TestCase):
    """Tests for compute_conformal_anomaly."""

    def test_returns_dict(self):
        from alpha_ladder_core.c2_derivation import compute_conformal_anomaly
        result = compute_conformal_anomaly(n=2)
        self.assertIsInstance(result, dict)

    def test_central_charge(self):
        from alpha_ladder_core.c2_derivation import compute_conformal_anomaly
        result = compute_conformal_anomaly(n=2)
        self.assertEqual(result["central_charge"], 1)

    def test_chi_s2(self):
        from alpha_ladder_core.c2_derivation import compute_conformal_anomaly
        result = compute_conformal_anomaly(n=2)
        self.assertEqual(result["chi"], 2)

    def test_integrated_anomaly(self):
        from alpha_ladder_core.c2_derivation import compute_conformal_anomaly
        result = compute_conformal_anomaly(n=2)
        self.assertAlmostEqual(result["integrated_anomaly"], 1.0 / 3.0, places=10)

    def test_inverse_anomaly(self):
        from alpha_ladder_core.c2_derivation import compute_conformal_anomaly
        result = compute_conformal_anomaly(n=2)
        self.assertEqual(result["inverse_anomaly"], 3)

    def test_matches_c2(self):
        from alpha_ladder_core.c2_derivation import compute_conformal_anomaly
        result = compute_conformal_anomaly(n=2)
        self.assertTrue(result["matches_c2"])


# ---------------------------------------------------------------------------
# 5. compute_seeley_dewitt_coefficients
# ---------------------------------------------------------------------------
class TestSeeleyDeWittCoefficients(unittest.TestCase):
    """Tests for compute_seeley_dewitt_coefficients."""

    def test_returns_dict(self):
        from alpha_ladder_core.c2_derivation import compute_seeley_dewitt_coefficients
        result = compute_seeley_dewitt_coefficients(n=2)
        self.assertIsInstance(result, dict)

    def test_scalar_a0(self):
        from alpha_ladder_core.c2_derivation import compute_seeley_dewitt_coefficients
        result = compute_seeley_dewitt_coefficients(n=2, field_type="scalar")
        self.assertAlmostEqual(result["a0"], 1.0, places=10)

    def test_scalar_a1_density(self):
        from alpha_ladder_core.c2_derivation import compute_seeley_dewitt_coefficients
        result = compute_seeley_dewitt_coefficients(n=2, field_type="scalar")
        self.assertAlmostEqual(result["a1_density"], 1.0 / 3.0, places=10)

    def test_kretschner_s2(self):
        from alpha_ladder_core.c2_derivation import compute_seeley_dewitt_coefficients
        result = compute_seeley_dewitt_coefficients(n=2, field_type="scalar")
        self.assertAlmostEqual(result["Kretschner"], 4.0, places=10)

    def test_ric_squared_s2(self):
        from alpha_ladder_core.c2_derivation import compute_seeley_dewitt_coefficients
        result = compute_seeley_dewitt_coefficients(n=2, field_type="scalar")
        self.assertAlmostEqual(result["Ric_squared"], 2.0, places=10)

    def test_scalar_a2_density(self):
        from alpha_ladder_core.c2_derivation import compute_seeley_dewitt_coefficients
        result = compute_seeley_dewitt_coefficients(n=2, field_type="scalar")
        # (1/180)(2*4 - 2*2 + 5*4) = (1/180)(8 - 4 + 20) = 24/180 = 2/15
        self.assertAlmostEqual(result["a2_density"], 2.0 / 15.0, places=10)

    def test_vector_components(self):
        from alpha_ladder_core.c2_derivation import compute_seeley_dewitt_coefficients
        result = compute_seeley_dewitt_coefficients(n=2, field_type="vector")
        self.assertEqual(result["n_components"], 2)

    def test_graviton_components(self):
        from alpha_ladder_core.c2_derivation import compute_seeley_dewitt_coefficients
        result = compute_seeley_dewitt_coefficients(n=2, field_type="graviton")
        self.assertEqual(result["n_components"], 3)


# ---------------------------------------------------------------------------
# 6. compute_spectral_zeta_values
# ---------------------------------------------------------------------------
class TestSpectralZetaValues(unittest.TestCase):
    """Tests for compute_spectral_zeta_values."""

    def test_returns_dict(self):
        from alpha_ladder_core.c2_derivation import compute_spectral_zeta_values
        result = compute_spectral_zeta_values(n=2)
        self.assertIsInstance(result, dict)

    def test_zeta_0_value(self):
        from alpha_ladder_core.c2_derivation import compute_spectral_zeta_values
        result = compute_spectral_zeta_values(n=2)
        # zeta(0) = -11/12 for scalar Laplacian on S^2 (excluding zero mode)
        self.assertAlmostEqual(result["zeta_0"], -11.0 / 12.0, places=6)

    def test_zeta_neg1_finite(self):
        from alpha_ladder_core.c2_derivation import compute_spectral_zeta_values
        result = compute_spectral_zeta_values(n=2)
        self.assertTrue(math.isfinite(result["zeta_neg1"]))

    def test_zeta_2_value(self):
        from alpha_ladder_core.c2_derivation import compute_spectral_zeta_values
        result = compute_spectral_zeta_values(n=2)
        # Telescoping identity: sum (2l+1)/[l(l+1)]^2 = 1
        self.assertAlmostEqual(result["zeta_2"], 1.0, places=6)

    def test_has_required_keys(self):
        from alpha_ladder_core.c2_derivation import compute_spectral_zeta_values
        result = compute_spectral_zeta_values(n=2)
        self.assertIn("zeta_0", result)
        self.assertIn("zeta_neg1", result)
        self.assertIn("zeta_2", result)


# ---------------------------------------------------------------------------
# 7. compute_killing_vector_analysis
# ---------------------------------------------------------------------------
class TestKillingVectorAnalysis(unittest.TestCase):
    """Tests for compute_killing_vector_analysis."""

    def test_returns_dict(self):
        from alpha_ladder_core.c2_derivation import compute_killing_vector_analysis
        result = compute_killing_vector_analysis(n=2)
        self.assertIsInstance(result, dict)

    def test_s2_killing_count(self):
        from alpha_ladder_core.c2_derivation import compute_killing_vector_analysis
        result = compute_killing_vector_analysis(n=2)
        self.assertEqual(result["n_killing"], 3)

    def test_s3_killing_count(self):
        from alpha_ladder_core.c2_derivation import compute_killing_vector_analysis
        result = compute_killing_vector_analysis(n=3)
        self.assertEqual(result["n_killing"], 6)

    def test_s4_killing_count(self):
        from alpha_ladder_core.c2_derivation import compute_killing_vector_analysis
        result = compute_killing_vector_analysis(n=4)
        self.assertEqual(result["n_killing"], 10)

    def test_s2_isometry_group(self):
        from alpha_ladder_core.c2_derivation import compute_killing_vector_analysis
        result = compute_killing_vector_analysis(n=2)
        self.assertIn("SO(3)", result["isometry_group"])

    def test_not_zero_mode(self):
        from alpha_ladder_core.c2_derivation import compute_killing_vector_analysis
        result = compute_killing_vector_analysis(n=2)
        self.assertFalse(result["is_zero_mode"])

    def test_not_harmonic(self):
        from alpha_ladder_core.c2_derivation import compute_killing_vector_analysis
        result = compute_killing_vector_analysis(n=2)
        self.assertFalse(result["is_harmonic"])

    def test_eigenvalue_positive_finite(self):
        from alpha_ladder_core.c2_derivation import compute_killing_vector_analysis
        result = compute_killing_vector_analysis(n=2)
        self.assertGreater(result["eigenvalue"], 0)
        self.assertTrue(math.isfinite(result["eigenvalue"]))


# ---------------------------------------------------------------------------
# 8. scan_spheres
# ---------------------------------------------------------------------------
class TestScanSpheres(unittest.TestCase):
    """Tests for scan_spheres."""

    def test_returns_dict(self):
        from alpha_ladder_core.c2_derivation import scan_spheres
        result = scan_spheres(n_max=8)
        self.assertIsInstance(result, dict)

    def test_scan_results_length(self):
        from alpha_ladder_core.c2_derivation import scan_spheres
        result = scan_spheres(n_max=8)
        self.assertGreaterEqual(len(result["scan_results"]), 7)

    def test_n2_all_candidates_equal_3(self):
        from alpha_ladder_core.c2_derivation import scan_spheres
        result = scan_spheres(n_max=8)
        n2_entry = None
        for entry in result["scan_results"]:
            if entry.get("n") == 2:
                n2_entry = entry
                break
        self.assertIsNotNone(n2_entry)
        # All candidate c_2 values should equal 3 for n=2
        for val in n2_entry.get("candidates", {}).values():
            self.assertAlmostEqual(val, 3.0, places=6)

    def test_n3_degeneracy_breaks(self):
        from alpha_ladder_core.c2_derivation import scan_spheres
        result = scan_spheres(n_max=8)
        n3_entry = None
        for entry in result["scan_results"]:
            if entry.get("n") == 3:
                n3_entry = entry
                break
        self.assertIsNotNone(n3_entry)
        # At least two candidates should disagree for n=3
        vals = list(n3_entry.get("candidates", {}).values())
        self.assertGreaterEqual(len(vals), 2)
        self.assertFalse(all(abs(v - vals[0]) < 1e-6 for v in vals))

    def test_unique_to_n2(self):
        from alpha_ladder_core.c2_derivation import scan_spheres
        result = scan_spheres(n_max=8)
        self.assertTrue(result["unique_to_n2"])

    def test_degeneracy_count(self):
        from alpha_ladder_core.c2_derivation import scan_spheres
        result = scan_spheres(n_max=8)
        self.assertGreaterEqual(result["degeneracy_count"], 5)


# ---------------------------------------------------------------------------
# 9. compute_loop_attempts
# ---------------------------------------------------------------------------
class TestComputeLoopAttempts(unittest.TestCase):
    """Tests for compute_loop_attempts."""

    def test_returns_dict(self):
        from alpha_ladder_core.c2_derivation import compute_loop_attempts
        result = compute_loop_attempts()
        self.assertIsInstance(result, dict)

    def test_at_least_five_attempts(self):
        from alpha_ladder_core.c2_derivation import compute_loop_attempts
        result = compute_loop_attempts()
        self.assertGreaterEqual(len(result["attempts"]), 5)

    def test_attempts_have_required_keys(self):
        from alpha_ladder_core.c2_derivation import compute_loop_attempts
        result = compute_loop_attempts()
        for attempt in result["attempts"]:
            self.assertIn("name", attempt)
            self.assertIn("result", attempt)
            self.assertIn("target", attempt)
            self.assertIn("status", attempt)
            self.assertIn("gap", attempt)

    def test_at_least_one_matches(self):
        from alpha_ladder_core.c2_derivation import compute_loop_attempts
        result = compute_loop_attempts()
        statuses = [a["status"].lower() for a in result["attempts"]]
        has_match = any("correct" in s or "matches" in s or "match" in s for s in statuses)
        self.assertTrue(has_match, f"No matching attempt found. Statuses: {statuses}")

    def test_at_least_one_fails(self):
        from alpha_ladder_core.c2_derivation import compute_loop_attempts
        result = compute_loop_attempts()
        statuses = [a["status"].lower() for a in result["attempts"]]
        has_fail = any("wrong" in s or "fail" in s or "miss" in s for s in statuses)
        self.assertTrue(has_fail, f"No failing attempt found. Statuses: {statuses}")

    def test_monopole_zeta_attempt(self):
        from alpha_ladder_core.c2_derivation import compute_loop_attempts
        result = compute_loop_attempts()
        monopole_found = False
        for attempt in result["attempts"]:
            if "monopole" in attempt["name"].lower() and "charged" in attempt["name"].lower():
                self.assertAlmostEqual(attempt["result"], 0.1, delta=0.05)
                monopole_found = True
                break
        self.assertTrue(monopole_found, "No monopole/zeta attempt found")

    def test_with_constants_parameter(self):
        from alpha_ladder_core.c2_derivation import compute_loop_attempts
        from alpha_ladder_core.constants import get_constants
        c = get_constants("CODATA 2018")
        result = compute_loop_attempts(constants=c)
        self.assertIsInstance(result, dict)


# ---------------------------------------------------------------------------
# 10. summarize_c2_derivation
# ---------------------------------------------------------------------------
class TestSummarizeC2Derivation(unittest.TestCase):
    """Tests for summarize_c2_derivation."""

    def test_returns_dict(self):
        from alpha_ladder_core.c2_derivation import summarize_c2_derivation
        result = summarize_c2_derivation()
        self.assertIsInstance(result, dict)

    def test_has_required_keys(self):
        from alpha_ladder_core.c2_derivation import summarize_c2_derivation
        result = summarize_c2_derivation()
        for key in ["best_candidate", "derivation_achieved", "gap_status",
                     "honest_assessment", "what_works", "what_fails",
                     "remaining_gap"]:
            self.assertIn(key, result, f"Missing key: {key}")

    def test_derivation_not_achieved(self):
        from alpha_ladder_core.c2_derivation import summarize_c2_derivation
        result = summarize_c2_derivation()
        self.assertFalse(result["derivation_achieved"])

    def test_best_candidate_references_chi_or_hrr_or_tangent(self):
        from alpha_ladder_core.c2_derivation import summarize_c2_derivation
        result = summarize_c2_derivation()
        bc = result["best_candidate"].lower()
        self.assertTrue(
            "chi" in bc or "hrr" in bc or "tangent" in bc or "hirzebruch" in bc,
            f"best_candidate '{result['best_candidate']}' does not reference expected terms"
        )

    def test_what_works_list(self):
        from alpha_ladder_core.c2_derivation import summarize_c2_derivation
        result = summarize_c2_derivation()
        self.assertIsInstance(result["what_works"], list)
        self.assertGreaterEqual(len(result["what_works"]), 3)

    def test_what_fails_list(self):
        from alpha_ladder_core.c2_derivation import summarize_c2_derivation
        result = summarize_c2_derivation()
        self.assertIsInstance(result["what_fails"], list)
        self.assertGreaterEqual(len(result["what_fails"]), 2)

    def test_honest_assessment_nonempty(self):
        from alpha_ladder_core.c2_derivation import summarize_c2_derivation
        result = summarize_c2_derivation()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 20)

    def test_with_constants_parameter(self):
        from alpha_ladder_core.c2_derivation import summarize_c2_derivation
        from alpha_ladder_core.constants import get_constants
        c = get_constants("CODATA 2018")
        result = summarize_c2_derivation(constants=c)
        self.assertIsInstance(result, dict)


# ---------------------------------------------------------------------------
# 11. Cross-consistency tests
# ---------------------------------------------------------------------------
class TestCrossConsistency(unittest.TestCase):
    """Cross-module consistency: various quantities must agree at n=2."""

    def test_inverse_a1_equals_chi_T(self):
        from alpha_ladder_core.c2_derivation import (
            compute_heat_kernel_a1,
            compute_chi_tangent_bundle,
        )
        hk = compute_heat_kernel_a1(n=2)
        chi = compute_chi_tangent_bundle(n=2)
        self.assertAlmostEqual(hk["inverse_a1"], chi["chi_T"], places=10)

    def test_inverse_a1_not_dim_so_for_s3(self):
        from alpha_ladder_core.c2_derivation import compute_heat_kernel_a1
        hk = compute_heat_kernel_a1(n=3)
        # dim SO(4) = 6, but 1/a1 = 1 for S^3
        self.assertNotAlmostEqual(hk["inverse_a1"], 6.0, places=5)

    def test_conformal_inverse_equals_inverse_a1(self):
        from alpha_ladder_core.c2_derivation import (
            compute_heat_kernel_a1,
            compute_conformal_anomaly,
        )
        hk = compute_heat_kernel_a1(n=2)
        ca = compute_conformal_anomaly(n=2)
        self.assertAlmostEqual(ca["inverse_anomaly"], hk["inverse_a1"], places=10)

    def test_gauge_matching_residual_small(self):
        from alpha_ladder_core.c2_derivation import compute_gauge_matching_identity
        result = compute_gauge_matching_identity()
        self.assertLess(abs(result["residual_ppm"]), 0.01)

    def test_seeley_dewitt_a1_matches_heat_kernel(self):
        from alpha_ladder_core.c2_derivation import (
            compute_heat_kernel_a1,
            compute_seeley_dewitt_coefficients,
        )
        hk = compute_heat_kernel_a1(n=2)
        sd = compute_seeley_dewitt_coefficients(n=2, field_type="scalar")
        self.assertAlmostEqual(hk["a1_density"], sd["a1_density"], places=10)


# ---------------------------------------------------------------------------
# 12. Generalization tests
# ---------------------------------------------------------------------------
class TestGeneralization(unittest.TestCase):
    """Test 1/a_1 = 6/(n(n-1)) sequence and degeneracy breaking."""

    def test_inverse_a1_sequence(self):
        from alpha_ladder_core.c2_derivation import compute_heat_kernel_a1
        expected = {2: 3.0, 3: 1.0, 4: 0.5, 5: 3.0 / 10.0}
        for n, val in expected.items():
            result = compute_heat_kernel_a1(n=n)
            self.assertAlmostEqual(
                result["inverse_a1"], val, places=6,
                msg=f"1/a_1 for S^{n}: expected {val}, got {result['inverse_a1']}"
            )

    def test_s2_triple_degeneracy(self):
        """For n=2: n(n+1)/2 = 1/a_1 = n+1 = d-1 (all = 3)."""
        from alpha_ladder_core.c2_derivation import compute_heat_kernel_a1
        n = 2
        result = compute_heat_kernel_a1(n=n)
        self.assertEqual(n * (n + 1) // 2, 3)
        self.assertAlmostEqual(result["inverse_a1"], 3.0, places=10)
        self.assertEqual(n + 1, 3)
        # d = n + 2 = 4, so d - 1 = 3
        self.assertEqual(n + 2 - 1, 3)

    def test_s3_degeneracy_broken(self):
        """For n=3: n(n+1)/2 = 6, 1/a_1 = 1, n+1 = 4 -- all different."""
        from alpha_ladder_core.c2_derivation import compute_heat_kernel_a1
        n = 3
        result = compute_heat_kernel_a1(n=n)
        vals = [n * (n + 1) // 2, result["inverse_a1"], n + 1]
        # All three should be different
        self.assertEqual(len(set([round(v, 6) for v in vals])), 3,
                         f"Expected all different but got {vals}")

    def test_killing_formula_n2_to_n5(self):
        """n_killing = (n+1)(n+2)/2 - 1 = dim SO(n+1) for S^n."""
        from alpha_ladder_core.c2_derivation import compute_killing_vector_analysis
        expected = {2: 3, 3: 6, 4: 10, 5: 15}
        for n, val in expected.items():
            result = compute_killing_vector_analysis(n=n)
            self.assertEqual(
                result["n_killing"], val,
                msg=f"Killing vectors for S^{n}: expected {val}, got {result['n_killing']}"
            )


# ---------------------------------------------------------------------------
# 13. compute_c3_vacuum_polynomial
# ---------------------------------------------------------------------------
class TestC3VacuumPolynomial(unittest.TestCase):
    """Tests for compute_c3_vacuum_polynomial."""

    def test_returns_dict(self):
        from alpha_ladder_core.c2_derivation import compute_c3_vacuum_polynomial
        result = compute_c3_vacuum_polynomial(d=4, D=6)
        self.assertIsInstance(result, dict)

    def test_c3_equals_phi_half(self):
        from alpha_ladder_core.c2_derivation import compute_c3_vacuum_polynomial
        result = compute_c3_vacuum_polynomial(d=4, D=6)
        phi_half = (1 + math.sqrt(5)) / 4
        self.assertAlmostEqual(result["c3_value"], phi_half, places=10)

    def test_r_plus(self):
        from alpha_ladder_core.c2_derivation import compute_c3_vacuum_polynomial
        result = compute_c3_vacuum_polynomial(d=4, D=6)
        self.assertAlmostEqual(result["r_plus"], -3 + math.sqrt(5), places=10)

    def test_r_minus(self):
        from alpha_ladder_core.c2_derivation import compute_c3_vacuum_polynomial
        result = compute_c3_vacuum_polynomial(d=4, D=6)
        self.assertAlmostEqual(result["r_minus"], -3 - math.sqrt(5), places=10)

    def test_C_0(self):
        from alpha_ladder_core.c2_derivation import compute_c3_vacuum_polynomial
        result = compute_c3_vacuum_polynomial(d=4, D=6)
        phi = (1 + math.sqrt(5)) / 2
        self.assertAlmostEqual(result["C_0"], phi**2 / 2, places=10)

    def test_c3_equals_phi_half_flag(self):
        from alpha_ladder_core.c2_derivation import compute_c3_vacuum_polynomial
        result = compute_c3_vacuum_polynomial(d=4, D=6)
        self.assertTrue(result["c3_equals_phi_half"])

    def test_ratio_c3_to_C0(self):
        from alpha_ladder_core.c2_derivation import compute_c3_vacuum_polynomial
        result = compute_c3_vacuum_polynomial(d=4, D=6)
        phi = (1 + math.sqrt(5)) / 2
        self.assertAlmostEqual(result["ratio_c3_to_C0"], 1.0 / phi, places=10)

    def test_formula_contains_expected_terms(self):
        from alpha_ladder_core.c2_derivation import compute_c3_vacuum_polynomial
        result = compute_c3_vacuum_polynomial(d=4, D=6)
        self.assertIn("r_+", result["formula"])
        self.assertIn("d", result["formula"])

    def test_status_contains_identified(self):
        from alpha_ladder_core.c2_derivation import compute_c3_vacuum_polynomial
        result = compute_c3_vacuum_polynomial(d=4, D=6)
        self.assertIn("IDENTIFIED", result["status"])

    def test_d3_D5_c3_equals_phi_half_check(self):
        from alpha_ladder_core.c2_derivation import compute_c3_vacuum_polynomial
        # D=5, d=3: x^2+5x+3=0, discriminant=13, NOT 5, so phi should not appear
        # But the formula (r_+ + d)/d may still give phi/2 if discriminant happens to match
        result = compute_c3_vacuum_polynomial(d=3, D=5)
        # Check that c3_equals_phi_half flag correctly reports whether it matches
        self.assertIsInstance(result["c3_equals_phi_half"], bool)


# ---------------------------------------------------------------------------
# 14. compute_pentagonal_polynomial
# ---------------------------------------------------------------------------
class TestPentagonalPolynomial(unittest.TestCase):
    """Tests for compute_pentagonal_polynomial."""

    def test_returns_dict(self):
        from alpha_ladder_core.c2_derivation import compute_pentagonal_polynomial
        result = compute_pentagonal_polynomial(d=4, D=6)
        self.assertIsInstance(result, dict)

    def test_transformed_poly_contains_expected_terms(self):
        from alpha_ladder_core.c2_derivation import compute_pentagonal_polynomial
        result = compute_pentagonal_polynomial(d=4, D=6)
        poly = result["transformed_poly"]
        self.assertIn("4f", poly)
        self.assertIn("2f", poly)
        self.assertIn("1", poly)

    def test_f_plus_is_phi_half(self):
        from alpha_ladder_core.c2_derivation import compute_pentagonal_polynomial
        result = compute_pentagonal_polynomial(d=4, D=6)
        self.assertTrue(result["f_plus_is_phi_half"])

    def test_is_cos_pi_5(self):
        from alpha_ladder_core.c2_derivation import compute_pentagonal_polynomial
        result = compute_pentagonal_polynomial(d=4, D=6)
        self.assertTrue(result["is_cos_pi_5"])

    def test_verification_residual(self):
        from alpha_ladder_core.c2_derivation import compute_pentagonal_polynomial
        result = compute_pentagonal_polynomial(d=4, D=6)
        self.assertLess(result["verification_residual"], 1e-15)

    def test_f_plus_value(self):
        from alpha_ladder_core.c2_derivation import compute_pentagonal_polynomial
        result = compute_pentagonal_polynomial(d=4, D=6)
        phi_half = (1 + math.sqrt(5)) / 4
        self.assertAlmostEqual(result["roots"]["f_plus"], phi_half, places=10)

    def test_f_minus_value(self):
        from alpha_ladder_core.c2_derivation import compute_pentagonal_polynomial
        result = compute_pentagonal_polynomial(d=4, D=6)
        f_minus = (1 - math.sqrt(5)) / 4
        self.assertAlmostEqual(result["roots"]["f_minus"], f_minus, places=10)

    def test_scan_is_list(self):
        from alpha_ladder_core.c2_derivation import compute_pentagonal_polynomial
        result = compute_pentagonal_polynomial(d=4, D=6)
        self.assertIsInstance(result["scan"], list)

    def test_scan_includes_d4_D6_phi(self):
        from alpha_ladder_core.c2_derivation import compute_pentagonal_polynomial
        result = compute_pentagonal_polynomial(d=4, D=6)
        found = False
        for entry in result["scan"]:
            if entry.get("D") == 6 and entry.get("d") == 4:
                found = True
                break
        self.assertTrue(found, "Scan should include (D=6, d=4) entry")

    def test_original_polynomial(self):
        from alpha_ladder_core.c2_derivation import compute_pentagonal_polynomial
        result = compute_pentagonal_polynomial(d=4, D=6)
        self.assertEqual(result["original_poly"], "x^2 + 6x + 4 = 0")


# ---------------------------------------------------------------------------
# 15. compute_rationality_nogo
# ---------------------------------------------------------------------------
class TestRationalityNogo(unittest.TestCase):
    """Tests for compute_rationality_nogo."""

    def test_returns_dict(self):
        from alpha_ladder_core.c2_derivation import compute_rationality_nogo
        result = compute_rationality_nogo()
        self.assertIsInstance(result, dict)

    def test_phi_half_irrational(self):
        from alpha_ladder_core.c2_derivation import compute_rationality_nogo
        result = compute_rationality_nogo()
        self.assertTrue(result["phi_half_irrational"])

    def test_conclusion_contains_cannot(self):
        from alpha_ladder_core.c2_derivation import compute_rationality_nogo
        result = compute_rationality_nogo()
        self.assertIn("cannot", result["conclusion"].lower())

    def test_source_contains_vacuum_and_polynomial(self):
        from alpha_ladder_core.c2_derivation import compute_rationality_nogo
        result = compute_rationality_nogo()
        src = result["source"].lower()
        self.assertIn("vacuum", src)
        self.assertIn("polynomial", src)

    def test_status_contains_proven(self):
        from alpha_ladder_core.c2_derivation import compute_rationality_nogo
        result = compute_rationality_nogo()
        self.assertIn("PROVEN", result["status"])

    def test_s2_spectral_data_is_dict(self):
        from alpha_ladder_core.c2_derivation import compute_rationality_nogo
        result = compute_rationality_nogo()
        self.assertIsInstance(result["s2_spectral_data"], dict)

    def test_s2_spectral_a_0(self):
        from alpha_ladder_core.c2_derivation import compute_rationality_nogo
        result = compute_rationality_nogo()
        self.assertEqual(result["s2_spectral_data"]["a_0"], 1)

    def test_s2_spectral_a_1(self):
        from alpha_ladder_core.c2_derivation import compute_rationality_nogo
        result = compute_rationality_nogo()
        self.assertAlmostEqual(result["s2_spectral_data"]["a_1"], 1.0 / 3.0, places=10)

    def test_s2_spectral_chi(self):
        from alpha_ladder_core.c2_derivation import compute_rationality_nogo
        result = compute_rationality_nogo()
        self.assertEqual(result["s2_spectral_data"]["chi"], 2)


# ---------------------------------------------------------------------------
# 16. compute_c3_coefficient_structure
# ---------------------------------------------------------------------------
class TestC3CoefficientStructure(unittest.TestCase):
    """Tests for compute_c3_coefficient_structure."""

    def test_returns_dict(self):
        from alpha_ladder_core.c2_derivation import compute_c3_coefficient_structure
        result = compute_c3_coefficient_structure()
        self.assertIsInstance(result, dict)

    def test_c2_value(self):
        from alpha_ladder_core.c2_derivation import compute_c3_coefficient_structure
        result = compute_c3_coefficient_structure()
        self.assertEqual(result["c2"], 3)

    def test_c3_close_to_phi_half(self):
        from alpha_ladder_core.c2_derivation import compute_c3_coefficient_structure
        result = compute_c3_coefficient_structure()
        self.assertAlmostEqual(result["c3"], 0.80902, places=5)

    def test_C_0_close_to_phi_sq_half(self):
        from alpha_ladder_core.c2_derivation import compute_c3_coefficient_structure
        result = compute_c3_coefficient_structure()
        self.assertAlmostEqual(result["C_0"], 1.30902, places=5)

    def test_ratio_c3_C0_close_to_inv_phi(self):
        from alpha_ladder_core.c2_derivation import compute_c3_coefficient_structure
        result = compute_c3_coefficient_structure()
        self.assertAlmostEqual(result["ratio_c3_C0"], 0.61803, places=5)

    def test_origin_c2(self):
        from alpha_ladder_core.c2_derivation import compute_c3_coefficient_structure
        result = compute_c3_coefficient_structure()
        origin = result["origin_c2"].lower()
        self.assertTrue("topology" in origin or "s^2" in origin,
                        f"origin_c2 '{result['origin_c2']}' missing expected terms")

    def test_origin_c3(self):
        from alpha_ladder_core.c2_derivation import compute_c3_coefficient_structure
        result = compute_c3_coefficient_structure()
        origin = result["origin_c3"].lower()
        self.assertTrue("vacuum" in origin or "polynomial" in origin,
                        f"origin_c3 '{result['origin_c3']}' missing expected terms")

    def test_suppression_ratio_small(self):
        from alpha_ladder_core.c2_derivation import compute_c3_coefficient_structure
        result = compute_c3_coefficient_structure()
        self.assertLess(result["suppression_ratio"], 0.01)

    def test_F_value_range(self):
        from alpha_ladder_core.c2_derivation import compute_c3_coefficient_structure
        result = compute_c3_coefficient_structure()
        self.assertGreater(result["F_value"], 1.0)
        self.assertLess(result["F_value"], 1.001)

    def test_alpha3_contribution_ppm(self):
        from alpha_ladder_core.c2_derivation import compute_c3_coefficient_structure
        result = compute_c3_coefficient_structure()
        self.assertGreater(result["alpha3_contribution_ppm"], 0)
        self.assertLess(result["alpha3_contribution_ppm"], 1.0)

    def test_honest_assessment_nonempty(self):
        from alpha_ladder_core.c2_derivation import compute_c3_coefficient_structure
        result = compute_c3_coefficient_structure()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 0)

    def test_with_constants_parameter(self):
        from alpha_ladder_core.c2_derivation import compute_c3_coefficient_structure
        from alpha_ladder_core.constants import get_constants
        c = get_constants("CODATA 2018")
        result = compute_c3_coefficient_structure(constants=c)
        self.assertIsInstance(result, dict)
        self.assertIn("c3", result)


# ---------------------------------------------------------------------------
# 17. c3 cross-consistency tests
# ---------------------------------------------------------------------------
class TestC3CrossConsistency(unittest.TestCase):
    """Cross-checks between c_3 and c_2 results."""

    def test_c3_vacuum_matches_c3_structure(self):
        from alpha_ladder_core.c2_derivation import (
            compute_c3_vacuum_polynomial,
            compute_c3_coefficient_structure,
        )
        vac = compute_c3_vacuum_polynomial(d=4, D=6)
        struct = compute_c3_coefficient_structure()
        self.assertAlmostEqual(vac["c3_value"], struct["c3"], places=10)

    def test_pentagonal_residual_tiny(self):
        from alpha_ladder_core.c2_derivation import compute_pentagonal_polynomial
        result = compute_pentagonal_polynomial(d=4, D=6)
        self.assertLess(result["verification_residual"], 1e-14)

    def test_inverse_a1_not_c3(self):
        from alpha_ladder_core.c2_derivation import (
            compute_heat_kernel_a1,
            compute_c3_vacuum_polynomial,
        )
        hk = compute_heat_kernel_a1(n=2)
        vac = compute_c3_vacuum_polynomial(d=4, D=6)
        # 1/a1 = 3, c3 = phi/2 ~ 0.809 -- they should NOT be equal
        self.assertNotAlmostEqual(1.0 / hk["a1_density"], vac["c3_value"], places=5)

    def test_summary_has_c3_key(self):
        from alpha_ladder_core.c2_derivation import summarize_c2_derivation
        result = summarize_c2_derivation()
        self.assertIn("c3_vacuum_polynomial", result)

    def test_summary_what_works_expanded(self):
        from alpha_ladder_core.c2_derivation import summarize_c2_derivation
        result = summarize_c2_derivation()
        self.assertGreaterEqual(len(result["what_works"]), 10)


if __name__ == "__main__":
    unittest.main()
