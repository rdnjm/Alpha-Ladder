"""
Comprehensive test suite for alpha_ladder_core.dark_sector module.

Tests cover:
- compute_dark_sector: dark matter coupling constants and derived quantities
- compute_wave_profile: wavefunction probability density profiles
- compute_alps_simulation: ALPS-style light-shining-through-wall parameters
- Cross-module consistency with ladder and experimental modules
"""

import math
import unittest
from decimal import Decimal

from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.dark_sector import (
    compute_alps_simulation,
    compute_dark_sector,
    compute_wave_profile,
    get_experimental_bounds,
)
from alpha_ladder_core.experimental import strategy_dark_sector
from alpha_ladder_core.ladder import calculate_geometric_rungs


class TestComputeDarkSector(unittest.TestCase):
    """Tests for compute_dark_sector return values and physical consistency."""

    @classmethod
    def setUpClass(cls):
        cls.constants = get_constants("CODATA 2018")
        cls.result = compute_dark_sector(cls.constants)

    def test_alpha_10_value(self):
        """alpha_10 (alpha^10) should be approximately 4.2821e-22."""
        alpha_10 = float(self.result["alpha_10"])
        self.assertTrue(
            math.isclose(alpha_10, 4.2821e-22, rel_tol=1e-3),
            f"alpha_10 = {alpha_10}, expected ~4.2821e-22",
        )

    def test_epsilon_value(self):
        """Kinetic mixing epsilon (alpha^5) should be approximately 2.0693e-11."""
        epsilon = float(self.result["epsilon"])
        self.assertTrue(
            math.isclose(epsilon, 2.0693e-11, rel_tol=1e-3),
            f"epsilon = {epsilon}, expected ~2.0693e-11",
        )

    def test_m_axion_eV(self):
        """Axion mass in eV should be approximately 2.188e-16 eV."""
        m_axion = float(self.result["m_axion_eV"])
        self.assertTrue(
            math.isclose(m_axion, 2.188e-16, rel_tol=1e-2),
            f"m_axion_eV = {m_axion}, expected ~2.188e-16",
        )

    def test_lambda_bar_km(self):
        """Reduced Compton wavelength in km should be approximately 901806 km."""
        lambda_bar = float(self.result["lambda_bar_km"])
        self.assertAlmostEqual(lambda_bar, 901806, delta=10,
                               msg=f"lambda_bar_km = {lambda_bar}, expected ~901806 km")

    def test_lambda_ratio(self):
        """lambda_ratio should be approximately 2.346."""
        ratio = float(self.result["lambda_ratio"])
        self.assertTrue(
            math.isclose(ratio, 2.346, rel_tol=1e-2),
            f"lambda_ratio = {ratio}, expected ~2.346",
        )

    def test_ratio_to_em(self):
        """ratio_to_em should be approximately 5.868e-20."""
        ratio = float(self.result["ratio_to_em"])
        self.assertTrue(
            math.isclose(ratio, 5.868e-20, rel_tol=1e-2),
            f"ratio_to_em = {ratio}, expected ~5.868e-20",
        )

    def test_ratio_to_gravity(self):
        """ratio_to_gravity should be greater than 1e23."""
        ratio = float(self.result["ratio_to_gravity"])
        self.assertGreater(ratio, 1e23,
                           f"ratio_to_gravity = {ratio}, expected > 1e23")

    def test_log10_values(self):
        """Logarithmic coupling values should match expected ranges."""
        log10_em = float(self.result["log10_em"])
        log10_dark = float(self.result["log10_dark"])
        log10_grav = float(self.result["log10_grav"])

        self.assertAlmostEqual(log10_em, -2.137, delta=0.1,
                               msg=f"log10_em = {log10_em}, expected ~-2.137")
        self.assertAlmostEqual(log10_dark, -21.37, delta=0.1,
                               msg=f"log10_dark = {log10_dark}, expected ~-21.37")
        self.assertAlmostEqual(log10_grav, -44.8, delta=0.1,
                               msg=f"log10_grav = {log10_grav}, expected ~-44.8")

    def test_return_keys(self):
        """All expected keys should be present in the returned dict."""
        expected_keys = {
            "alpha_10",
            "epsilon",
            "m_axion_eV",
            "lambda_bar_km",
            "lambda_ratio",
            "ratio_to_em",
            "ratio_to_gravity",
            "log10_em",
            "log10_dark",
            "log10_grav",
        }
        self.assertTrue(
            expected_keys.issubset(set(self.result.keys())),
            f"Missing keys: {expected_keys - set(self.result.keys())}",
        )


class TestComputeWaveProfile(unittest.TestCase):
    """Tests for compute_wave_profile wavefunction output."""

    @classmethod
    def setUpClass(cls):
        cls.constants = get_constants("CODATA 2018")
        cls.result = compute_wave_profile(cls.constants)

    def test_array_lengths(self):
        """Default n_points=500 should produce arrays of length 500."""
        self.assertEqual(len(self.result["x_km"]), 500)
        self.assertEqual(len(self.result["psi_squared"]), 500)

    def test_custom_n_points(self):
        """Custom n_points=100 should produce arrays of length 100."""
        result = compute_wave_profile(self.constants, n_points=100)
        self.assertEqual(len(result["x_km"]), 100)
        self.assertEqual(len(result["psi_squared"]), 100)

    def test_psi_squared_non_negative(self):
        """Probability density psi_squared must be non-negative everywhere."""
        for i, val in enumerate(self.result["psi_squared"]):
            self.assertGreaterEqual(
                float(val), 0.0,
                f"psi_squared[{i}] = {val} is negative",
            )

    def test_psi_squared_max_is_one(self):
        """Peak of normalized psi_squared should be approximately 1.0."""
        peak = max(float(v) for v in self.result["psi_squared"])
        self.assertAlmostEqual(peak, 1.0, delta=0.01,
                               msg=f"max(psi_squared) = {peak}, expected ~1.0")

    def test_wavelength_matches_compute(self):
        """lambda_bar_km from wave profile should match compute_dark_sector."""
        ds_result = compute_dark_sector(self.constants)
        wp_lambda = float(self.result["lambda_bar_km"])
        ds_lambda = float(ds_result["lambda_bar_km"])
        self.assertTrue(
            math.isclose(wp_lambda, ds_lambda, rel_tol=1e-10),
            f"Wave profile lambda={wp_lambda} != dark sector lambda={ds_lambda}",
        )

    def test_alpha_override(self):
        """Overriding alpha should produce a different wavelength."""
        result_override = compute_wave_profile(
            self.constants, alpha_override=Decimal("0.008")
        )
        default_lambda = float(self.result["lambda_bar_km"])
        override_lambda = float(result_override["lambda_bar_km"])
        self.assertNotAlmostEqual(
            default_lambda, override_lambda, places=3,
            msg="alpha_override should change the wavelength",
        )


class TestComputeAlpsSimulation(unittest.TestCase):
    """Tests for compute_alps_simulation ALPS experiment parameters."""

    @classmethod
    def setUpClass(cls):
        cls.constants = get_constants("CODATA 2018")
        cls.result = compute_alps_simulation(cls.constants)
        cls.alpha_float = float(cls.constants.alpha)

    def test_epsilon_is_alpha_5(self):
        """ALPS epsilon should equal alpha^5."""
        epsilon = float(self.result["epsilon"])
        expected = self.alpha_float ** 5
        self.assertTrue(
            math.isclose(epsilon, expected, rel_tol=1e-10),
            f"epsilon = {epsilon}, expected alpha^5 = {expected}",
        )

    def test_P_tunneling_is_alpha_20(self):
        """Tunneling probability should equal alpha^20."""
        p_tunnel = float(self.result["P_tunneling"])
        expected = self.alpha_float ** 20
        self.assertTrue(
            math.isclose(p_tunnel, expected, rel_tol=1e-6),
            f"P_tunneling = {p_tunnel}, expected alpha^20 = {expected}",
        )

    def test_photons_needed_positive(self):
        """photons_needed should be a very large positive number."""
        photons = float(self.result["photons_needed"])
        self.assertGreater(photons, 1e10,
                           f"photons_needed = {photons}, expected very large positive")

    def test_alps_parameters(self):
        """ALPS magnet parameters should be B=5.3 T and L=106 m."""
        self.assertAlmostEqual(float(self.result["B_tesla"]), 5.3, places=1,
                               msg="B_tesla should be 5.3")
        self.assertAlmostEqual(float(self.result["L_meters"]), 106, places=0,
                               msg="L_meters should be 106")

    def test_return_keys(self):
        """All expected keys should be present in the returned dict."""
        expected_keys = {
            "epsilon",
            "P_tunneling",
            "photons_needed",
            "B_tesla",
            "L_meters",
            "laser_power_W",
            "photons_per_sec",
            "wait_years",
            "log10_wait_years",
            "alps_ii_sensitivity",
            "orders_below_alps",
        }
        self.assertTrue(
            expected_keys.issubset(set(self.result.keys())),
            f"Missing keys: {expected_keys - set(self.result.keys())}",
        )

    def test_wait_years_astronomically_large(self):
        """Wait time should be vastly longer than the age of the universe."""
        log10_wait = self.result["log10_wait_years"]
        # Age of universe is ~10^10 years; wait should be >> that
        self.assertGreater(log10_wait, 15,
                           f"log10(wait_years) = {log10_wait}, expected >> 10")

    def test_orders_below_alps_positive(self):
        """Rung 10 prediction should be many orders below ALPS II sensitivity."""
        orders = self.result["orders_below_alps"]
        self.assertGreater(orders, 5,
                           f"orders_below_alps = {orders}, expected > 5")


class TestCrossModuleConsistency(unittest.TestCase):
    """Tests verifying consistency between dark_sector, ladder, and experimental modules."""

    @classmethod
    def setUpClass(cls):
        cls.constants = get_constants("CODATA 2018")

    def test_alpha_dm_matches_ladder(self):
        """alpha_dm from ladder module should match alpha_10 from dark_sector."""
        rungs, bridges, alpha_dm, pi_bridge = calculate_geometric_rungs(self.constants)
        ds_result = compute_dark_sector(self.constants)

        ladder_val = float(alpha_dm)
        dark_val = float(ds_result["alpha_10"])

        self.assertTrue(
            math.isclose(ladder_val, dark_val, rel_tol=1e-10),
            f"ladder alpha_dm={ladder_val} != dark_sector alpha_10={dark_val}",
        )

    def test_epsilon_matches_experimental(self):
        """epsilon from dark_sector should match experimental module's epsilon_predicted."""
        ds_result = compute_dark_sector(self.constants)
        exp_result = strategy_dark_sector(self.constants)

        ds_epsilon = float(ds_result["epsilon"])
        exp_epsilon = float(exp_result["epsilon_predicted"])

        self.assertTrue(
            math.isclose(ds_epsilon, exp_epsilon, rel_tol=1e-6),
            f"dark_sector epsilon={ds_epsilon} != experimental epsilon={exp_epsilon}",
        )


if __name__ == "__main__":
    unittest.main()
