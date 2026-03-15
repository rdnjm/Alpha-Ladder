"""Tests for alpha_ladder_core.flux_stabilization module.

Validates flux stabilization physics: minimum exists for all N >= 1,
produces Planck-scale masses, and resolves theoretical gaps.
"""

import math
import unittest

from alpha_ladder_core.flux_stabilization import (
    compute_flux_coefficient,
    compute_flux_potential,
    find_flux_minimum,
    compute_flux_dilaton_mass,
    scan_flux_quanta,
    compute_flux_gap_closure,
    summarize_flux_stabilization,
)


class TestFluxCoefficient(unittest.TestCase):
    """Tests for compute_flux_coefficient."""

    def test_returns_dict(self):
        result = compute_flux_coefficient(N=1)
        self.assertIsInstance(result, dict)

    def test_c_positive_for_n_ge_1(self):
        for n in [1, 2, 3, 5, 10, 100]:
            result = compute_flux_coefficient(N=n)
            self.assertGreater(result["C"], 0, f"C should be > 0 for N={n}")

    def test_c_scales_as_n_squared(self):
        r1 = compute_flux_coefficient(N=1)
        r2 = compute_flux_coefficient(N=2)
        ratio = r2["C"] / r1["C"]
        self.assertAlmostEqual(ratio, 4.0, places=5,
                               msg="C should scale as N^2")

    def test_c_zero_for_n_zero(self):
        result = compute_flux_coefficient(N=0)
        self.assertEqual(result["C"], 0)

    def test_has_n_and_a0_keys(self):
        result = compute_flux_coefficient(N=3)
        self.assertIn("N", result)
        self.assertIn("a_0", result)


class TestFluxPotential(unittest.TestCase):
    """Tests for compute_flux_potential."""

    def test_returns_dict(self):
        result = compute_flux_potential(N=1)
        self.assertIsInstance(result, dict)

    def test_has_sigma_grid(self):
        result = compute_flux_potential(N=1)
        self.assertIn("sigma_grid", result)
        self.assertGreater(len(result["sigma_grid"]), 1)

    def test_has_potential_components(self):
        result = compute_flux_potential(N=1)
        for key in ["V_casimir", "V_curvature", "V_flux", "V_total"]:
            self.assertIn(key, result)
            self.assertIsInstance(result[key], list)

    def test_v_total_is_sum_of_components(self):
        result = compute_flux_potential(N=1)
        for i in range(len(result["V_total"])):
            expected = (result["V_casimir"][i]
                        + result["V_curvature"][i]
                        + result["V_flux"][i])
            self.assertAlmostEqual(result["V_total"][i], expected, places=10,
                                   msg=f"V_total[{i}] != sum of components")

    def test_a_coefficient_negative(self):
        result = compute_flux_potential(N=1)
        self.assertLess(result["A"], 0, "Casimir coefficient A should be < 0")

    def test_b_coefficient_negative(self):
        result = compute_flux_potential(N=1)
        self.assertLess(result["B"], 0,
                        "Curvature coefficient B for S^2 should be < 0")

    def test_c_coefficient_positive(self):
        result = compute_flux_potential(N=1)
        self.assertGreater(result["C"], 0,
                           "Flux coefficient C should be > 0")


class TestFluxMinimum(unittest.TestCase):
    """Tests for find_flux_minimum."""

    def test_returns_dict(self):
        result = find_flux_minimum(N=1)
        self.assertIsInstance(result, dict)

    def test_minimum_exists_n1(self):
        result = find_flux_minimum(N=1)
        self.assertTrue(result["minimum_exists"])

    def test_minimum_exists_n5(self):
        result = find_flux_minimum(N=5)
        self.assertTrue(result["minimum_exists"])

    def test_sigma_0_is_finite(self):
        result = find_flux_minimum(N=1)
        self.assertTrue(math.isfinite(result["sigma_0"]))

    def test_v_double_prime_positive(self):
        result = find_flux_minimum(N=1)
        self.assertGreater(result["V_double_prime"], 0,
                           "V'' > 0 confirms a minimum, not a maximum")

    def test_a_stabilized_positive(self):
        result = find_flux_minimum(N=1)
        self.assertGreater(result["a_stabilized"], 0)

    def test_coefficients_stored(self):
        result = find_flux_minimum(N=1)
        for key in ["A", "B", "C"]:
            self.assertIn(key, result)

    def test_discriminant_positive(self):
        result = find_flux_minimum(N=1)
        self.assertGreater(result["discriminant"], 0)


class TestFluxDilatonMass(unittest.TestCase):
    """Tests for compute_flux_dilaton_mass."""

    def test_returns_dict(self):
        result = compute_flux_dilaton_mass(N=1)
        self.assertIsInstance(result, dict)

    def test_first_principles_true(self):
        result = compute_flux_dilaton_mass(N=1)
        self.assertTrue(result["first_principles"],
                        "Minimum exists so first_principles should be True")

    def test_mass_positive(self):
        result = compute_flux_dilaton_mass(N=1)
        self.assertGreater(result["m_phi_eV"], 0)

    def test_mass_planck_scale(self):
        result = compute_flux_dilaton_mass(N=1)
        self.assertGreater(result["m_phi_eV"], 1e25,
                           "Flux-stabilized mass should be Planck scale")

    def test_compton_wavelength_tiny(self):
        result = compute_flux_dilaton_mass(N=1)
        self.assertLess(result["lambda_compton_m"], 1e-30,
                        "Compton wavelength should be sub-Planckian")

    def test_exceeds_threshold_true(self):
        result = compute_flux_dilaton_mass(N=1)
        self.assertTrue(result["exceeds_threshold"],
                        "Planck-scale mass trivially exceeds any threshold")

    def test_mass_scale_planck(self):
        result = compute_flux_dilaton_mass(N=1)
        self.assertEqual(result["mass_scale"], "planck_scale")


class TestScanFluxQuanta(unittest.TestCase):
    """Tests for scan_flux_quanta."""

    def test_returns_dict(self):
        result = scan_flux_quanta(N_max=5)
        self.assertIsInstance(result, dict)

    def test_results_length(self):
        n_max = 5
        result = scan_flux_quanta(N_max=n_max)
        self.assertGreaterEqual(len(result["results"]), n_max)

    def test_all_minima_exist(self):
        result = scan_flux_quanta(N_max=5)
        for entry in result["results"]:
            self.assertTrue(entry["minimum_exists"],
                            f"Minimum should exist for N={entry.get('N')}")

    def test_all_masses_positive(self):
        result = scan_flux_quanta(N_max=5)
        for entry in result["results"]:
            self.assertGreater(entry["m_phi_eV"], 0)

    def test_mass_decreases_with_n(self):
        result = scan_flux_quanta(N_max=5)
        masses = [entry["m_phi_eV"] for entry in result["results"]]
        for i in range(len(masses) - 1):
            self.assertGreater(masses[i], masses[i + 1],
                               "Mass should decrease with increasing N")

    def test_all_stable(self):
        result = scan_flux_quanta(N_max=5)
        self.assertTrue(result["all_stable"])


class TestFluxGapClosure(unittest.TestCase):
    """Tests for compute_flux_gap_closure."""

    def test_returns_dict(self):
        result = compute_flux_gap_closure()
        self.assertIsInstance(result, dict)

    def test_gap3_resolved(self):
        result = compute_flux_gap_closure()
        self.assertTrue(result["gap3_resolved"])

    def test_gap1_resolved(self):
        result = compute_flux_gap_closure()
        self.assertTrue(result["gap1_resolved"])

    def test_dilaton_decoupled(self):
        result = compute_flux_gap_closure()
        self.assertTrue(result["dilaton_effectively_decoupled"])

    def test_has_mechanism_strings(self):
        result = compute_flux_gap_closure()
        self.assertIsInstance(result["gap1_mechanism"], str)
        self.assertIsInstance(result["gap3_mechanism"], str)
        self.assertGreater(len(result["gap1_mechanism"]), 0)
        self.assertGreater(len(result["gap3_mechanism"]), 0)

    def test_has_honest_assessment(self):
        result = compute_flux_gap_closure()
        self.assertIsInstance(result["honest_assessment"], str)
        self.assertGreater(len(result["honest_assessment"]), 0)


class TestSummarize(unittest.TestCase):
    """Tests for summarize_flux_stabilization."""

    def test_returns_dict(self):
        result = summarize_flux_stabilization()
        self.assertIsInstance(result, dict)

    def test_has_overall_assessment(self):
        result = summarize_flux_stabilization()
        self.assertIn("overall_assessment", result)
        self.assertIsInstance(result["overall_assessment"], str)

    def test_first_principles_true(self):
        result = summarize_flux_stabilization()
        self.assertTrue(result["first_principles"])

    def test_has_sub_dicts(self):
        result = summarize_flux_stabilization()
        for key in ["potential", "minimum", "mass", "scan", "gap_closure"]:
            self.assertIn(key, result)
            self.assertIsInstance(result[key], dict)

    def test_gap_closure_both_resolved(self):
        result = summarize_flux_stabilization()
        gc = result["gap_closure"]
        self.assertTrue(gc["gap1_resolved"])
        self.assertTrue(gc["gap3_resolved"])

    def test_has_flux_quantum_n(self):
        result = summarize_flux_stabilization()
        self.assertIn("flux_quantum_N", result)


if __name__ == "__main__":
    unittest.main()
