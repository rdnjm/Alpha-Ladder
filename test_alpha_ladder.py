"""Tests for the Alpha Ladder geometric relationships."""

import unittest
from decimal import Decimal, getcontext

getcontext().prec = 50

# Import everything from the main module
from alpha_ladder import (
    alpha, alpha_g, pi, phi, e,
    lambda_bar_c, a_0,
    r_e_from_compton, r_e_from_bohr, r_e_nist,
    calculate_geometric_rungs,
)


class TestFoundationalConstants(unittest.TestCase):
    """Verify constants are correctly defined."""

    def test_alpha_value(self):
        self.assertAlmostEqual(float(alpha), 0.0072973525664, places=13)

    def test_alpha_reciprocal(self):
        recip = float(1 / alpha)
        self.assertAlmostEqual(recip, 137.036, places=2)

    def test_phi_identity(self):
        # phi^2 = phi + 1 (defining property)
        self.assertAlmostEqual(float(phi ** 2), float(phi + 1), places=40)

    def test_phi_value(self):
        self.assertAlmostEqual(float(phi), 1.6180339887, places=9)

    def test_pi_value(self):
        self.assertAlmostEqual(float(pi), 3.14159265358979, places=14)


class TestElectronGeometry(unittest.TestCase):
    """r_e = alpha * lambda_bar_c = alpha^2 * a_0"""

    def test_compton_path(self):
        # alpha * lambda_bar_c should give classical electron radius
        self.assertAlmostEqual(float(r_e_from_compton), float(r_e_nist), places=19)

    def test_bohr_path(self):
        # alpha^2 * a_0 should give classical electron radius
        self.assertAlmostEqual(float(r_e_from_bohr), float(r_e_nist), places=19)

    def test_both_paths_agree(self):
        self.assertAlmostEqual(float(r_e_from_compton), float(r_e_from_bohr), places=19)


class TestThe42OrderGap(unittest.TestCase):
    """alpha / alpha_G ~ 4.17e42"""

    def test_gap_order_of_magnitude(self):
        gap = float(alpha / alpha_g)
        self.assertGreater(gap, 1e42)
        self.assertLess(gap, 1e43)

    def test_gap_value(self):
        gap = float(alpha / alpha_g)
        self.assertAlmostEqual(gap / 1e42, 4.17, places=1)


class TestGeometricLadder(unittest.TestCase):
    """Verify the ladder structure and rung labels."""

    def setUp(self):
        self.rungs, self.bridges, self.alpha_dm, self.pi_bridge = calculate_geometric_rungs()

    def test_rung_count(self):
        self.assertEqual(len(self.rungs), 24)

    def test_rungs_decrease_monotonically(self):
        for i in range(1, len(self.rungs)):
            self.assertLess(self.rungs[i]['value'], self.rungs[i - 1]['value'])

    def test_dark_matter_label_at_10(self):
        rung_10 = self.rungs[9]  # 0-indexed
        self.assertEqual(rung_10['power'], 10)
        self.assertIn("Dark Matter", rung_10['label'])

    def test_gravity_label_at_21(self):
        rung_21 = self.rungs[20]
        self.assertEqual(rung_21['power'], 21)
        self.assertIn("Gravity", rung_21['label'])

    def test_alpha_power_correctness(self):
        # Spot-check: alpha^1 = alpha
        self.assertEqual(self.rungs[0]['value'], alpha)
        # alpha^2
        self.assertEqual(self.rungs[1]['value'], alpha ** 2)

    def test_ratio_to_gravity(self):
        for r in self.rungs:
            expected_ratio = r['value'] / alpha_g
            self.assertEqual(r['ratio_to_gravity'], expected_ratio)


class TestDarkSector(unittest.TestCase):
    """alpha_DM = alpha^10"""

    def setUp(self):
        _, _, self.alpha_dm, _ = calculate_geometric_rungs()

    def test_dark_matter_value(self):
        self.assertEqual(self.alpha_dm, alpha ** 10)

    def test_dark_matter_order(self):
        val = float(self.alpha_dm)
        self.assertGreater(val, 1e-23)
        self.assertLess(val, 1e-21)


class TestGeometricBridge(unittest.TestCase):
    """alpha_G = (phi^2 / 2) * alpha^21  -- the corrected bridge."""

    def setUp(self):
        _, self.bridges, _, _ = calculate_geometric_rungs()

    def test_original_phi_bridge_residual_is_large(self):
        pred = self.bridges["phi * alpha^21  (original)"]
        residual = abs(float((alpha_g - pred) / alpha_g)) * 100
        self.assertGreater(residual, 20.0, "Original bridge should have ~23% residual")

    def test_phi_squared_half_bridge(self):
        pred = self.bridges["(phi^2 / 2) * alpha^21"]
        residual = abs(float((alpha_g - pred) / alpha_g)) * 100
        self.assertLess(residual, 0.02, "phi^2/2 bridge should be < 0.02% residual")

    def test_five_twelfths_pi_bridge(self):
        pred = self.bridges["(5/12) * pi * alpha^21"]
        residual = abs(float((alpha_g - pred) / alpha_g)) * 100
        self.assertLess(residual, 0.02)

    def test_sqrt_e_cbrt2_bridge(self):
        pred = self.bridges["sqrt(e) / cbrt(2) * alpha^21"]
        residual = abs(float((alpha_g - pred) / alpha_g)) * 100
        self.assertLess(residual, 0.02)

    def test_phi_squared_half_is_best_phi_form(self):
        # phi^2/2 should be much better than bare phi
        orig = self.bridges["phi * alpha^21  (original)"]
        corrected = self.bridges["(phi^2 / 2) * alpha^21"]
        err_orig = abs(float(alpha_g - orig))
        err_corr = abs(float(alpha_g - corrected))
        self.assertLess(err_corr, err_orig / 1000)


class TestPiBridge(unittest.TestCase):
    """Pi-Bridge: alpha^21 * pi^2"""

    def setUp(self):
        _, _, _, self.pi_bridge = calculate_geometric_rungs()

    def test_pi_bridge_value(self):
        expected = alpha ** 21 * pi ** 2
        self.assertEqual(self.pi_bridge, expected)

    def test_pi_bridge_order(self):
        val = float(self.pi_bridge)
        self.assertGreater(val, 1e-45)
        self.assertLess(val, 1e-43)


class TestInternalConsistency(unittest.TestCase):
    """Cross-checks between different parts of the ladder."""

    def test_alpha_dm_times_alpha_11_near_alpha_g(self):
        # alpha^10 * alpha^11 = alpha^21, so alpha_DM * alpha^11 ~ alpha^21
        _, _, alpha_dm, _ = calculate_geometric_rungs()
        product = alpha_dm * alpha ** 11
        self.assertEqual(product, alpha ** 21)

    def test_gap_spans_21_rungs(self):
        # alpha / alpha_G = alpha^1 / (C * alpha^21) = 1/(C * alpha^20)
        # So gap * alpha^20 = 1/C where C ~ 1.309
        gap = float(alpha / alpha_g)
        alpha_20 = float(alpha ** 20)
        inv_c = gap * alpha_20
        # inv_c should be 1/1.309 ~ 0.764
        self.assertAlmostEqual(inv_c, 1.0 / 1.3088, places=3)


if __name__ == "__main__":
    unittest.main()
