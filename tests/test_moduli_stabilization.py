"""
Tests for alpha_ladder_core/moduli_stabilization.py

Moduli stabilization for compact internal manifolds in the Alpha Ladder
framework.  Verifies that all moduli (volume and shape) acquire positive
mass-squared through Gauss-Bonnet and curvature mechanisms, ensuring a
stable compactification.
"""

import math
import unittest

from alpha_ladder_core.moduli_stabilization import (
    describe_moduli_space,
    compute_volume_stabilization,
    compute_shape_stabilization,
    compute_moduli_mass_spectrum,
    summarize_stabilization,
)


# -----------------------------------------------------------------------
# Tests for describe_moduli_space
# -----------------------------------------------------------------------

class TestDescribeModuliSpace(unittest.TestCase):
    """Tests for describe_moduli_space."""

    def test_genus2_real_moduli(self):
        """Genus-2 surface has 6*2 - 6 = 6 real moduli."""
        result = describe_moduli_space(genus=2)
        self.assertEqual(result["real_moduli"], 6)

    def test_genus2_complex_moduli(self):
        """Genus-2 surface has 3*2 - 3 = 3 complex structure moduli."""
        result = describe_moduli_space(genus=2)
        self.assertEqual(result["complex_structure_moduli"], 3)

    def test_genus2_total(self):
        """Total moduli for genus-2: 6 real + 1 volume = 7."""
        result = describe_moduli_space(genus=2)
        self.assertEqual(result["total_moduli"], 7)

    def test_genus2_hyperelliptic(self):
        """Every genus-2 Riemann surface is hyperelliptic."""
        result = describe_moduli_space(genus=2)
        self.assertTrue(result["hyperelliptic"])

    def test_genus3_real_moduli(self):
        """Genus-3 surface has 6*3 - 6 = 12 real moduli."""
        result = describe_moduli_space(genus=3)
        self.assertEqual(result["real_moduli"], 12)

    def test_genus3_not_necessarily_hyperelliptic(self):
        """For genus >= 3, the surface is not necessarily hyperelliptic."""
        result = describe_moduli_space(genus=3)
        # hyperelliptic may be False for g >= 3
        self.assertFalse(result["hyperelliptic"])

    def test_genus2_teichmuller_dim(self):
        """Teichmuller dimension for genus-2 is 6*2 - 6 = 6."""
        result = describe_moduli_space(genus=2)
        self.assertEqual(result["teichmuller_dim"], 6)

    def test_genus2_volume_modulus(self):
        """There is always exactly 1 volume modulus."""
        result = describe_moduli_space(genus=2)
        self.assertEqual(result["volume_modulus"], 1)

    def test_returns_dict(self):
        """describe_moduli_space must return a dict."""
        result = describe_moduli_space(genus=2)
        self.assertIsInstance(result, dict)


# -----------------------------------------------------------------------
# Tests for compute_volume_stabilization
# -----------------------------------------------------------------------

class TestVolumeStabilization(unittest.TestCase):
    """Tests for compute_volume_stabilization."""

    def test_is_stable(self):
        """Volume modulus must be stabilized."""
        result = compute_volume_stabilization(genus=2)
        self.assertTrue(result["stable"])

    def test_sigma_vev_zero(self):
        """The dilaton VEV sigma should be stabilized at zero."""
        result = compute_volume_stabilization(genus=2)
        self.assertAlmostEqual(result["sigma_vev"], 0.0, places=10)

    def test_has_mechanism(self):
        """Result must include a stabilization mechanism description."""
        result = compute_volume_stabilization(genus=2)
        self.assertIn("mechanism", result)
        self.assertIsInstance(result["mechanism"], str)

    def test_with_mass(self):
        """When m_phi_eV is provided, mass info should be numeric."""
        result = compute_volume_stabilization(genus=2, m_phi_eV=1e-5)
        # The result should contain some mass-related numeric value
        has_numeric = False
        for key, val in result.items():
            if isinstance(val, (int, float)) and val != 0:
                has_numeric = True
                break
        self.assertTrue(has_numeric)

    def test_without_mass(self):
        """When m_phi_eV is None, still returns a valid result."""
        result = compute_volume_stabilization(genus=2, m_phi_eV=None)
        self.assertIsInstance(result, dict)
        self.assertIn("stable", result)

    def test_returns_dict(self):
        """compute_volume_stabilization must return a dict."""
        result = compute_volume_stabilization(genus=2)
        self.assertIsInstance(result, dict)

    def test_casimir_result_parameter_accepted(self):
        """compute_volume_stabilization should accept casimir_result parameter."""
        result = compute_volume_stabilization(genus=2, casimir_result={"minimum_exists": False, "honest_assessment": "test"})
        self.assertIsInstance(result, dict)

    def test_casimir_result_none(self):
        """casimir_result=None should work (backward compatible)."""
        result = compute_volume_stabilization(genus=2, casimir_result=None)
        self.assertEqual(result["mechanism"], "dilaton_mass")

    def test_casimir_result_no_minimum(self):
        """When casimir_result has minimum_exists=False, mechanism stays dilaton_mass."""
        result = compute_volume_stabilization(
            genus=2, casimir_result={"minimum_exists": False, "honest_assessment": "No stable minimum found"}
        )
        self.assertEqual(result["mechanism"], "dilaton_mass")
        self.assertIn("casimir_attempted", result)


# -----------------------------------------------------------------------
# Tests for compute_shape_stabilization
# -----------------------------------------------------------------------

class TestShapeStabilization(unittest.TestCase):
    """Tests for compute_shape_stabilization."""

    def test_genus2_stable(self):
        """Shape moduli for genus-2 must be stable."""
        result = compute_shape_stabilization(genus=2)
        self.assertTrue(result["stable"])

    def test_all_positive(self):
        """All shape moduli mass-squared values must be positive."""
        result = compute_shape_stabilization(genus=2)
        self.assertTrue(result["all_positive"])

    def test_constant_curvature_minimum(self):
        """Constant curvature metric must be a minimum of the action."""
        result = compute_shape_stabilization(genus=2)
        self.assertTrue(result["constant_curvature_is_minimum"])

    def test_lichnerowicz_bound(self):
        """Lichnerowicz eigenvalue bound must be >= 2.0."""
        result = compute_shape_stabilization(genus=2)
        self.assertGreaterEqual(result["lichnerowicz_bound"], 2.0)

    def test_six_moduli_masses(self):
        """Genus-2 has 6*2 - 6 = 6 shape moduli masses."""
        result = compute_shape_stabilization(genus=2)
        self.assertEqual(len(result["moduli_masses"]), 6)

    def test_all_masses_positive(self):
        """Every individual shape modulus mass-squared must be > 0."""
        result = compute_shape_stabilization(genus=2)
        for entry in result["moduli_masses"]:
            self.assertGreater(entry["mass_squared"], 0)

    def test_hessian_positive_definite(self):
        """The Hessian of the effective potential must be positive definite."""
        result = compute_shape_stabilization(genus=2)
        self.assertTrue(result["hessian_positive_definite"])

    def test_genus3_has_12_moduli(self):
        """Genus-3 has 6*3 - 6 = 12 shape moduli masses."""
        result = compute_shape_stabilization(genus=3)
        self.assertEqual(len(result["moduli_masses"]), 12)

    def test_default_lambda_gb(self):
        """Works correctly with lambda_GB=None (default)."""
        result = compute_shape_stabilization(genus=2, lambda_GB=None)
        self.assertIsInstance(result, dict)
        self.assertIn("stable", result)

    def test_specific_lambda_gb(self):
        """Works correctly with lambda_GB = sqrt(5) - 3."""
        lambda_val = math.sqrt(5) - 3
        result = compute_shape_stabilization(genus=2, lambda_GB=lambda_val)
        self.assertIsInstance(result, dict)
        self.assertIn("stable", result)


# -----------------------------------------------------------------------
# Tests for compute_moduli_mass_spectrum
# -----------------------------------------------------------------------

class TestModuliMassSpectrum(unittest.TestCase):
    """Tests for compute_moduli_mass_spectrum."""

    def test_genus2_total_moduli(self):
        """Genus-2 has 6 shape + 1 volume = 7 total moduli."""
        result = compute_moduli_mass_spectrum(genus=2)
        self.assertEqual(result["total_moduli"], 7)

    def test_all_stable(self):
        """All moduli must be stable."""
        result = compute_moduli_mass_spectrum(genus=2)
        self.assertTrue(result["all_stable"])

    def test_has_volume_modulus(self):
        """Result must contain a volume_modulus entry."""
        result = compute_moduli_mass_spectrum(genus=2)
        self.assertIn("volume_modulus", result)

    def test_has_shape_moduli(self):
        """Result must contain shape_moduli with 6 entries for genus-2."""
        result = compute_moduli_mass_spectrum(genus=2)
        self.assertIn("shape_moduli", result)
        self.assertEqual(len(result["shape_moduli"]), 6)

    def test_lightest_modulus_positive(self):
        """The lightest modulus must have positive mass-squared."""
        result = compute_moduli_mass_spectrum(genus=2)
        all_masses = []
        vol_msq = result.get("volume_modulus", {}).get("mass_squared")
        if vol_msq is not None:
            all_masses.append(vol_msq)
        for entry in result.get("shape_moduli", []):
            msq = entry.get("mass_squared")
            if msq is not None:
                all_masses.append(msq)
        self.assertTrue(len(all_masses) > 0)
        self.assertGreater(min(all_masses), 0)

    def test_mass_hierarchy_finite(self):
        """Mass hierarchy (heaviest/lightest) must be a finite positive number."""
        result = compute_moduli_mass_spectrum(genus=2)
        hierarchy = result.get("mass_hierarchy")
        self.assertIsNotNone(hierarchy)
        self.assertTrue(math.isfinite(hierarchy))
        self.assertGreater(hierarchy, 0)

    def test_default_parameters(self):
        """Works with all default parameters."""
        result = compute_moduli_mass_spectrum()
        self.assertIsInstance(result, dict)
        self.assertIn("total_moduli", result)


# -----------------------------------------------------------------------
# Tests for summarize_stabilization
# -----------------------------------------------------------------------

class TestSummarizeStabilization(unittest.TestCase):
    """Tests for summarize_stabilization."""

    def test_all_stable(self):
        """Overall stabilization summary must report all_stable True."""
        result = summarize_stabilization(genus=2)
        self.assertTrue(result["all_stable"])

    def test_volume_stabilized(self):
        """Volume modulus must be reported as stabilized."""
        result = summarize_stabilization(genus=2)
        self.assertTrue(result["volume_stabilized"])

    def test_shape_stabilized(self):
        """Shape moduli must be reported as stabilized."""
        result = summarize_stabilization(genus=2)
        self.assertTrue(result["shape_stabilized"])

    def test_has_summary(self):
        """Result must include a non-empty summary string."""
        result = summarize_stabilization(genus=2)
        self.assertIn("summary", result)
        self.assertIsInstance(result["summary"], str)
        self.assertTrue(len(result["summary"]) > 0)

    def test_genus2(self):
        """Works correctly for genus=2."""
        result = summarize_stabilization(genus=2)
        self.assertIsInstance(result, dict)

    def test_returns_dict(self):
        """summarize_stabilization must return a dict."""
        result = summarize_stabilization(genus=2)
        self.assertIsInstance(result, dict)


if __name__ == "__main__":
    unittest.main()
