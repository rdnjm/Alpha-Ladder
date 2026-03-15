"""Tests for alpha_ladder_core.ladder module.

Migrated from legacy/test_alpha_ladder.py, adapted to use
get_constants() for the constants namespace.
"""

import unittest
from decimal import Decimal, getcontext

getcontext().prec = 50

from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.ladder import (
    calculate_geometric_rungs,
    compute_electron_geometry,
    compute_gap,
    compute_physical_rung_map,
)


class TestElectronGeometry(unittest.TestCase):
    """r_e = alpha * lambda_bar_c = alpha^2 * a_0"""

    def setUp(self):
        self.c = get_constants("CODATA 2014")
        self.geom = compute_electron_geometry(self.c)

    def test_compton_path(self):
        self.assertAlmostEqual(
            float(self.geom["r_e_from_compton"]),
            float(self.geom["r_e_nist"]),
            places=19,
        )

    def test_bohr_path(self):
        self.assertAlmostEqual(
            float(self.geom["r_e_from_bohr"]),
            float(self.geom["r_e_nist"]),
            places=19,
        )

    def test_both_paths_agree(self):
        self.assertAlmostEqual(
            float(self.geom["r_e_from_compton"]),
            float(self.geom["r_e_from_bohr"]),
            places=19,
        )

    def test_returns_all_keys(self):
        for key in ("r_e_from_compton", "r_e_from_bohr", "r_e_nist"):
            self.assertIn(key, self.geom)
            self.assertIsInstance(self.geom[key], Decimal)


class TestThe42OrderGap(unittest.TestCase):
    """alpha / alpha_G ~ 4.17e42"""

    def setUp(self):
        self.c = get_constants("CODATA 2014")

    def test_gap_order_of_magnitude(self):
        gap = float(compute_gap(self.c))
        self.assertGreater(gap, 1e42)
        self.assertLess(gap, 1e43)

    def test_gap_value(self):
        gap = float(compute_gap(self.c))
        self.assertAlmostEqual(gap / 1e42, 4.17, places=1)


class TestGeometricLadder(unittest.TestCase):
    """Verify the ladder structure and rung labels."""

    def setUp(self):
        self.c = get_constants("CODATA 2014")
        self.rungs, self.bridges, self.alpha_dm, self.pi_bridge = (
            calculate_geometric_rungs(self.c)
        )

    def test_rung_count(self):
        self.assertEqual(len(self.rungs), 24)

    def test_rungs_decrease_monotonically(self):
        for i in range(1, len(self.rungs)):
            self.assertLess(
                self.rungs[i]["value"], self.rungs[i - 1]["value"]
            )

    def test_dark_matter_label_at_10(self):
        rung_10 = self.rungs[9]
        self.assertEqual(rung_10["power"], 10)
        self.assertIn("Dark Matter", rung_10["label"])

    def test_gravity_label_at_21(self):
        rung_21 = self.rungs[20]
        self.assertEqual(rung_21["power"], 21)
        self.assertIn("Gravity", rung_21["label"])

    def test_alpha_power_correctness(self):
        self.assertEqual(self.rungs[0]["value"], self.c.alpha)
        self.assertEqual(self.rungs[1]["value"], self.c.alpha ** 2)

    def test_ratio_to_gravity(self):
        for r in self.rungs:
            expected = r["value"] / self.c.alpha_g
            self.assertEqual(r["ratio_to_gravity"], expected)


class TestDarkSector(unittest.TestCase):
    """alpha_DM = alpha^10"""

    def setUp(self):
        self.c = get_constants("CODATA 2014")
        _, _, self.alpha_dm, _ = calculate_geometric_rungs(self.c)

    def test_dark_matter_value(self):
        self.assertEqual(self.alpha_dm, self.c.alpha ** 10)

    def test_dark_matter_order(self):
        val = float(self.alpha_dm)
        self.assertGreater(val, 1e-23)
        self.assertLess(val, 1e-21)


class TestGeometricBridge(unittest.TestCase):
    """Bridge residual tests (using CODATA 2014 values to match legacy)."""

    def setUp(self):
        self.c = get_constants("CODATA 2014")
        _, self.bridges, _, _ = calculate_geometric_rungs(self.c)

    def test_original_phi_bridge_residual_is_large(self):
        pred = self.bridges["φ · α²¹  (original)"]
        residual = abs(float((self.c.alpha_g - pred) / self.c.alpha_g)) * 100
        self.assertGreater(residual, 20.0)

    def test_phi_squared_half_bridge(self):
        pred = self.bridges["(φ²/2) · α²¹"]
        residual = abs(float((self.c.alpha_g - pred) / self.c.alpha_g)) * 100
        self.assertLess(residual, 0.02)

    def test_five_twelfths_pi_bridge(self):
        pred = self.bridges["(5/12) · π · α²¹"]
        residual = abs(float((self.c.alpha_g - pred) / self.c.alpha_g)) * 100
        self.assertLess(residual, 0.02)

    def test_sqrt_e_cbrt2_bridge(self):
        pred = self.bridges["√e / ³√2 · α²¹"]
        residual = abs(float((self.c.alpha_g - pred) / self.c.alpha_g)) * 100
        self.assertLess(residual, 0.02)

    def test_phi_squared_half_is_best_phi_form(self):
        orig = self.bridges["φ · α²¹  (original)"]
        corrected = self.bridges["(φ²/2) · α²¹"]
        err_orig = abs(float(self.c.alpha_g - orig))
        err_corr = abs(float(self.c.alpha_g - corrected))
        self.assertLess(err_corr, err_orig / 1000)


class TestPiBridge(unittest.TestCase):
    """Pi-Bridge: alpha^21 * pi^2"""

    def setUp(self):
        self.c = get_constants("CODATA 2014")
        _, _, _, self.pi_bridge = calculate_geometric_rungs(self.c)

    def test_pi_bridge_value(self):
        expected = self.c.alpha ** 21 * self.c.pi ** 2
        self.assertEqual(self.pi_bridge, expected)

    def test_pi_bridge_order(self):
        val = float(self.pi_bridge)
        self.assertGreater(val, 1e-45)
        self.assertLess(val, 1e-43)


class TestInternalConsistency(unittest.TestCase):
    """Cross-checks between different parts of the ladder."""

    def setUp(self):
        self.c = get_constants("CODATA 2014")

    def test_alpha_dm_times_alpha_11_near_alpha_21(self):
        _, _, alpha_dm, _ = calculate_geometric_rungs(self.c)
        product = alpha_dm * self.c.alpha ** 11
        self.assertEqual(product, self.c.alpha ** 21)

    def test_gap_spans_21_rungs(self):
        gap = float(self.c.alpha / self.c.alpha_g)
        alpha_20 = float(self.c.alpha ** 20)
        inv_c = gap * alpha_20
        self.assertAlmostEqual(inv_c, 1.0 / 1.3088, places=3)


class TestWithCODATA2018(unittest.TestCase):
    """Ensure the ladder works with CODATA 2018 constants too."""

    def setUp(self):
        self.c = get_constants("CODATA 2018")

    def test_rungs_build(self):
        rungs, bridges, alpha_dm, pi_bridge = calculate_geometric_rungs(self.c)
        self.assertEqual(len(rungs), 24)
        self.assertEqual(len(bridges), 4)

    def test_electron_geometry(self):
        geom = compute_electron_geometry(self.c)
        self.assertAlmostEqual(
            float(geom["r_e_from_compton"]),
            float(geom["r_e_nist"]),
            places=19,
        )

    def test_gap_order(self):
        gap = float(compute_gap(self.c))
        self.assertGreater(gap, 1e42)
        self.assertLess(gap, 1e43)


class TestPhysicalRungMap(unittest.TestCase):
    """Tests for the physical rung map."""

    def setUp(self):
        self.c = get_constants("CODATA 2018")
        self.result = compute_physical_rung_map(self.c)

    def test_returns_dict(self):
        self.assertIsInstance(self.result, dict)

    def test_has_31_rungs(self):
        # k=0..30 = 31 entries
        self.assertEqual(len(self.result["rungs"]), 31)

    def test_rung_0_is_bohr_radius(self):
        rung_0 = self.result["rungs"][0]
        self.assertEqual(rung_0["k"], 0)
        # length should be ~5.29e-11 m (Bohr radius)
        self.assertAlmostEqual(rung_0["length_m"], 5.29e-11, delta=0.01e-11)

    def test_rung_1_is_compton(self):
        rung_1 = self.result["rungs"][1]
        # length should be ~3.86e-13 m (reduced Compton wavelength)
        self.assertAlmostEqual(rung_1["length_m"], 3.86e-13, delta=0.01e-13)

    def test_rung_2_is_classical_radius(self):
        rung_2 = self.result["rungs"][2]
        # length should be ~2.82e-15 m (classical electron radius)
        self.assertAlmostEqual(rung_2["length_m"], 2.82e-15, delta=0.01e-15)

    def test_lengths_decrease_monotonically(self):
        rungs = self.result["rungs"]
        for i in range(1, len(rungs)):
            self.assertLess(rungs[i]["length_m"], rungs[i-1]["length_m"])

    def test_energies_increase_monotonically(self):
        rungs = self.result["rungs"]
        for i in range(1, len(rungs)):
            self.assertGreater(rungs[i]["energy_eV"], rungs[i-1]["energy_eV"])

    def test_particle_rungs_has_13_entries(self):
        self.assertEqual(len(self.result["particle_rungs"]), 13)

    def test_electron_at_rung_0(self):
        electron = next(p for p in self.result["particle_rungs"] if p["name"] == "Electron")
        self.assertAlmostEqual(electron["rung_position"], 0.0, delta=0.01)

    def test_proton_at_rung_1_53(self):
        proton = next(p for p in self.result["particle_rungs"] if p["name"] == "Proton")
        self.assertAlmostEqual(proton["rung_position"], 1.53, delta=0.05)

    def test_gut_scale_at_rung_9(self):
        gut = next(s for s in self.result["scale_rungs"] if "GUT" in s["name"])
        self.assertAlmostEqual(gut["rung_position"], 9.0, delta=0.2)

    def test_planck_length_at_rung_11_5(self):
        # k_planck is for Planck LENGTH: l_Pl lands at rung ~11.47
        # (Planck MASS in energy units is at ~10.47, but k_planck uses length)
        self.assertAlmostEqual(self.result["k_planck"], 11.5, delta=0.2)

    def test_gravity_at_rung_21(self):
        self.assertAlmostEqual(self.result["k_gravity"], 20.95, delta=0.1)

    def test_key_matches_nonempty(self):
        self.assertGreater(len(self.result["key_matches"]), 0)

    def test_base_identity_holds(self):
        bi = self.result["base_identity"]
        self.assertAlmostEqual(bi["r_e_from_compton"], bi["r_e_from_bohr"], delta=1e-25)

    def test_electroweak_cluster_at_2_5(self):
        # W, Z, Higgs, Top should all be between rung 2.3 and 2.7
        ew_particles = ["W boson", "Z boson", "Higgs boson", "Top quark"]
        for name in ew_particles:
            p = next(pr for pr in self.result["particle_rungs"] if pr["name"] == name)
            self.assertGreater(p["rung_position"], 2.3)
            self.assertLess(p["rung_position"], 2.7)


if __name__ == "__main__":
    unittest.main()
