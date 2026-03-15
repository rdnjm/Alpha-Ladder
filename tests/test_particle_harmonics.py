"""Tests for alpha_ladder_core.particle_harmonics module.

Validates particle mass dictionary and alpha-power harmonic rung
computations against known physical values.
"""

import math
import unittest

from alpha_ladder_core.constants import get_constants, get_particle_masses
from alpha_ladder_core.particle_harmonics import compute_harmonics


class TestGetParticleMasses(unittest.TestCase):
    """Tests for the get_particle_masses helper."""

    def setUp(self):
        self.masses = get_particle_masses()

    def test_returns_dict(self):
        self.assertIsInstance(self.masses, dict)

    def test_has_13_particles(self):
        self.assertEqual(len(self.masses), 13)

    def test_electron_mass(self):
        self.assertAlmostEqual(self.masses["Electron"], 0.511, places=3)

    def test_proton_mass(self):
        self.assertAlmostEqual(self.masses["Proton"], 938.272, places=2)

    def test_top_quark_mass(self):
        self.assertEqual(self.masses["Top quark"], 172760.0)

    def test_all_masses_positive(self):
        for name, m in self.masses.items():
            with self.subTest(particle=name):
                self.assertGreater(m, 0.0)

    def test_expected_particle_names(self):
        expected = {
            "Electron", "Muon", "Tau",
            "Up quark", "Down quark", "Strange quark",
            "Charm quark", "Bottom quark", "Top quark",
            "W boson", "Z boson", "Higgs boson", "Proton",
        }
        self.assertEqual(set(self.masses.keys()), expected)

    def test_returns_fresh_copy(self):
        m1 = get_particle_masses()
        m2 = get_particle_masses()
        self.assertIsNot(m1, m2)

    def test_all_values_are_floats(self):
        for name, m in self.masses.items():
            with self.subTest(particle=name):
                self.assertIsInstance(m, float)


class TestComputeHarmonicsWithNamespace(unittest.TestCase):
    """Tests for compute_harmonics using SimpleNamespace constants."""

    def setUp(self):
        self.constants = get_constants("CODATA 2018")
        self.results = compute_harmonics(self.constants)

    def test_returns_list(self):
        self.assertIsInstance(self.results, list)

    def test_returns_13_entries(self):
        self.assertEqual(len(self.results), 13)

    def test_each_entry_has_required_keys(self):
        required = {"name", "mass", "ratio", "rung_n", "nearest_half", "closeness", "is_match"}
        for entry in self.results:
            with self.subTest(particle=entry["name"]):
                self.assertEqual(set(entry.keys()), required)

    def test_electron_ratio_is_one(self):
        electron = self._find("Electron")
        self.assertAlmostEqual(electron["ratio"], 1.0, places=10)

    def test_electron_rung_n_is_zero(self):
        electron = self._find("Electron")
        self.assertAlmostEqual(electron["rung_n"], 0.0, places=10)

    def test_electron_nearest_half_is_zero(self):
        electron = self._find("Electron")
        self.assertAlmostEqual(electron["nearest_half"], 0.0, places=10)

    def test_electron_closeness_is_zero(self):
        electron = self._find("Electron")
        self.assertAlmostEqual(electron["closeness"], 0.0, places=10)

    def test_electron_is_match(self):
        electron = self._find("Electron")
        self.assertTrue(electron["is_match"])

    def test_proton_ratio(self):
        proton = self._find("Proton")
        masses = get_particle_masses()
        self.assertAlmostEqual(proton["ratio"], masses["Proton"] / masses["Electron"], places=1)

    def test_proton_rung_n(self):
        proton = self._find("Proton")
        masses = get_particle_masses()
        expected_n = math.log(masses["Proton"] / masses["Electron"]) / math.log(1.0 / float(self.constants.alpha))
        self.assertAlmostEqual(proton["rung_n"], expected_n, places=5)
        self.assertAlmostEqual(proton["rung_n"], 1.53, delta=0.01)

    def test_all_ratios_ge_one(self):
        for entry in self.results:
            with self.subTest(particle=entry["name"]):
                self.assertGreaterEqual(entry["ratio"], 1.0)

    def test_all_rung_n_ge_zero(self):
        for entry in self.results:
            with self.subTest(particle=entry["name"]):
                self.assertGreaterEqual(entry["rung_n"], 0.0)

    def test_nearest_half_is_half_integer(self):
        for entry in self.results:
            with self.subTest(particle=entry["name"]):
                doubled = entry["nearest_half"] * 2
                self.assertAlmostEqual(doubled, round(doubled), places=10)

    def test_closeness_is_nonnegative(self):
        for entry in self.results:
            with self.subTest(particle=entry["name"]):
                self.assertGreaterEqual(entry["closeness"], 0.0)

    def test_closeness_formula(self):
        for entry in self.results:
            with self.subTest(particle=entry["name"]):
                expected = abs(entry["rung_n"] - entry["nearest_half"])
                self.assertAlmostEqual(entry["closeness"], expected, places=10)

    def test_is_match_threshold(self):
        for entry in self.results:
            with self.subTest(particle=entry["name"]):
                expected_match = entry["closeness"] < 0.05
                self.assertEqual(entry["is_match"], expected_match)

    def test_rung_n_formula(self):
        alpha = float(self.constants.alpha)
        inv_alpha = 1.0 / alpha
        masses = get_particle_masses()
        m_e = masses["Electron"]
        for entry in self.results:
            with self.subTest(particle=entry["name"]):
                expected_n = math.log(entry["mass"] / m_e) / math.log(inv_alpha)
                self.assertAlmostEqual(entry["rung_n"], expected_n, places=10)

    def test_top_quark_mass_in_results(self):
        top = self._find("Top quark")
        self.assertEqual(top["mass"], 172760.0)

    def test_names_match_mass_dict(self):
        masses = get_particle_masses()
        result_names = {e["name"] for e in self.results}
        self.assertEqual(result_names, set(masses.keys()))

    def test_masses_match_mass_dict(self):
        masses = get_particle_masses()
        for entry in self.results:
            with self.subTest(particle=entry["name"]):
                self.assertEqual(entry["mass"], masses[entry["name"]])

    def _find(self, name):
        for entry in self.results:
            if entry["name"] == name:
                return entry
        self.fail(f"Particle '{name}' not found in results")


class TestComputeHarmonicsWithDict(unittest.TestCase):
    """Tests for compute_harmonics using a plain dict for constants."""

    def test_accepts_dict_with_alpha_key(self):
        constants = {"alpha": 0.00729735}
        results = compute_harmonics(constants)
        self.assertEqual(len(results), 13)

    def test_dict_produces_same_results_as_namespace(self):
        ns = get_constants("CODATA 2018")
        dict_c = {"alpha": float(ns.alpha)}
        results_ns = compute_harmonics(ns)
        results_dict = compute_harmonics(dict_c)
        for r_ns, r_dict in zip(results_ns, results_dict):
            self.assertAlmostEqual(r_ns["rung_n"], r_dict["rung_n"], places=10)
            self.assertEqual(r_ns["is_match"], r_dict["is_match"])


class TestComputeHarmonicsOrdering(unittest.TestCase):
    """Tests for result ordering and structure."""

    def setUp(self):
        self.constants = get_constants("CODATA 2018")
        self.results = compute_harmonics(self.constants)

    def test_order_matches_mass_dict_insertion_order(self):
        masses = get_particle_masses()
        expected_order = list(masses.keys())
        actual_order = [e["name"] for e in self.results]
        self.assertEqual(actual_order, expected_order)

    def test_electron_is_first(self):
        self.assertEqual(self.results[0]["name"], "Electron")

    def test_last_matches_dict_order(self):
        masses = get_particle_masses()
        expected_last = list(masses.keys())[-1]
        self.assertEqual(self.results[-1]["name"], expected_last)


class TestComputeHarmonicsEdgeCases(unittest.TestCase):
    """Edge-case and boundary tests."""

    def test_is_match_values_are_booleans(self):
        results = compute_harmonics(get_constants("CODATA 2018"))
        for entry in results:
            with self.subTest(particle=entry["name"]):
                self.assertIsInstance(entry["is_match"], bool)

    def test_at_least_one_match_exists(self):
        results = compute_harmonics(get_constants("CODATA 2018"))
        matches = [e for e in results if e["is_match"]]
        self.assertGreater(len(matches), 0, "Expected at least one alpha-harmonic match")

    def test_closeness_bounded_by_quarter(self):
        """Closeness to nearest half-integer should never exceed 0.25."""
        results = compute_harmonics(get_constants("CODATA 2018"))
        for entry in results:
            with self.subTest(particle=entry["name"]):
                self.assertLessEqual(entry["closeness"], 0.25 + 1e-10)

    def test_codata_2014_also_works(self):
        results = compute_harmonics(get_constants("CODATA 2014"))
        self.assertEqual(len(results), 13)
        electron = next(e for e in results if e["name"] == "Electron")
        self.assertAlmostEqual(electron["rung_n"], 0.0, places=10)


if __name__ == "__main__":
    unittest.main()
