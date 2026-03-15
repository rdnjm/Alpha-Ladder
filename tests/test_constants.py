"""Tests for alpha_ladder_core.constants module."""

import unittest
from decimal import Decimal, getcontext

getcontext().prec = 50

from alpha_ladder_core.constants import (
    CODATA_EDITIONS,
    DEFAULT_EDITION,
    available_editions,
    get_constants,
    get_particle_masses,
)


class TestAvailableEditions(unittest.TestCase):
    """Verify edition discovery."""

    def test_returns_list(self):
        result = available_editions()
        self.assertIsInstance(result, list)

    def test_contains_both_editions(self):
        editions = available_editions()
        self.assertIn("CODATA 2014", editions)
        self.assertIn("CODATA 2018", editions)

    def test_default_edition_is_2014(self):
        self.assertEqual(DEFAULT_EDITION, "CODATA 2014")


class TestGetConstantsLoads(unittest.TestCase):
    """Both editions should load without error and expose all attributes."""

    _REQUIRED_PHYSICAL = [
        "alpha", "alpha_g", "m_e", "hbar", "c", "G",
        "m_p", "m_mu", "k_B", "e_charge", "epsilon_0",
        "lambda_bar_c", "a_0", "r_e_nist",
    ]

    _REQUIRED_MATH = ["pi", "phi", "e", "ln2", "sqrt2", "sqrt3", "sqrt5"]

    def _check_edition(self, edition):
        ns = get_constants(edition)
        for attr in self._REQUIRED_PHYSICAL + self._REQUIRED_MATH:
            self.assertTrue(
                hasattr(ns, attr),
                f"{edition} missing attribute {attr!r}",
            )
            val = getattr(ns, attr)
            self.assertIsInstance(val, Decimal, f"{attr} should be Decimal")

    def test_codata_2014_loads(self):
        self._check_edition("CODATA 2014")

    def test_codata_2018_loads(self):
        self._check_edition("CODATA 2018")

    def test_unknown_edition_raises(self):
        with self.assertRaises(ValueError):
            get_constants("CODATA 1999")

    def test_default_edition_loads(self):
        ns = get_constants()
        self.assertIsInstance(ns.alpha, Decimal)


class TestPhysicalValues(unittest.TestCase):
    """Spot-check selected constant values."""

    def test_alpha_2014(self):
        ns = get_constants("CODATA 2014")
        self.assertEqual(ns.alpha, Decimal('0.0072973525664'))

    def test_alpha_2018(self):
        ns = get_constants("CODATA 2018")
        self.assertEqual(ns.alpha, Decimal('0.0072973525693'))

    def test_c_is_exact(self):
        for edition in available_editions():
            ns = get_constants(edition)
            self.assertEqual(ns.c, Decimal('299792458'))

    def test_alpha_reciprocal_near_137(self):
        ns = get_constants()
        recip = float(1 / ns.alpha)
        self.assertAlmostEqual(recip, 137.036, places=2)


class TestMathConstants(unittest.TestCase):
    """Mathematical constants must be consistent across editions
    and satisfy known identities."""

    def test_phi_squared_equals_phi_plus_one(self):
        ns = get_constants()
        self.assertAlmostEqual(
            float(ns.phi ** 2), float(ns.phi + 1), places=40
        )

    def test_phi_value(self):
        ns = get_constants()
        self.assertAlmostEqual(float(ns.phi), 1.6180339887, places=9)

    def test_pi_value(self):
        ns = get_constants()
        self.assertAlmostEqual(float(ns.pi), 3.14159265358979, places=14)

    def test_euler_e_value(self):
        ns = get_constants()
        self.assertAlmostEqual(float(ns.e), 2.71828182845904, places=14)

    def test_ln2_value(self):
        ns = get_constants()
        self.assertAlmostEqual(float(ns.ln2), 0.693147180559945, places=15)

    def test_sqrt2_squared(self):
        ns = get_constants()
        self.assertAlmostEqual(float(ns.sqrt2 ** 2), 2.0, places=40)

    def test_sqrt3_squared(self):
        ns = get_constants()
        self.assertAlmostEqual(float(ns.sqrt3 ** 2), 3.0, places=40)

    def test_sqrt5_squared(self):
        ns = get_constants()
        self.assertAlmostEqual(float(ns.sqrt5 ** 2), 5.0, places=40)

    def test_math_constants_same_across_editions(self):
        ns14 = get_constants("CODATA 2014")
        ns18 = get_constants("CODATA 2018")
        for attr in ("pi", "phi", "e", "ln2", "sqrt2", "sqrt3", "sqrt5"):
            self.assertEqual(
                getattr(ns14, attr),
                getattr(ns18, attr),
                f"{attr} differs across editions",
            )


class TestAlphaGComputation(unittest.TestCase):
    """For CODATA 2018, alpha_g must be computed from G, m_e, hbar, c."""

    def test_alpha_g_computed_from_G(self):
        ns = get_constants("CODATA 2018")
        expected = ns.G * ns.m_e ** 2 / (ns.hbar * ns.c)
        self.assertEqual(ns.alpha_g, expected)

    def test_alpha_g_order_of_magnitude(self):
        ns = get_constants("CODATA 2018")
        val = float(ns.alpha_g)
        self.assertGreater(val, 1e-46)
        self.assertLess(val, 1e-44)

    def test_alpha_g_2014_is_literal(self):
        ns = get_constants("CODATA 2014")
        self.assertEqual(ns.alpha_g, Decimal('1.7512e-45'))


class TestGetParticleMasses(unittest.TestCase):
    """Tests for the shared particle mass dictionary."""

    def setUp(self):
        self.masses = get_particle_masses()

    def test_returns_dict(self):
        self.assertIsInstance(self.masses, dict)

    def test_has_13_entries(self):
        self.assertEqual(len(self.masses), 13)

    def test_all_values_positive(self):
        for name, mass in self.masses.items():
            with self.subTest(particle=name):
                self.assertGreater(mass, 0.0)

    def test_electron_mass(self):
        self.assertAlmostEqual(self.masses["Electron"], 0.511, places=3)

    def test_proton_mass(self):
        self.assertAlmostEqual(self.masses["Proton"], 938.27, places=1)

    def test_returns_fresh_copy(self):
        m1 = get_particle_masses()
        m2 = get_particle_masses()
        self.assertIsNot(m1, m2)

    def test_all_values_are_floats(self):
        for name, mass in self.masses.items():
            with self.subTest(particle=name):
                self.assertIsInstance(mass, float)

    def test_expected_particle_names(self):
        expected = {
            "Electron", "Muon", "Tau",
            "Up quark", "Down quark", "Strange quark",
            "Charm quark", "Bottom quark", "Top quark",
            "Proton", "W boson", "Z boson", "Higgs boson",
        }
        self.assertEqual(set(self.masses.keys()), expected)


if __name__ == "__main__":
    unittest.main()
