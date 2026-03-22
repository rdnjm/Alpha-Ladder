"""Tests for escape_routes.py -- Escape Routes for KK Gauge Matching vs Alpha Running."""
import unittest
import math



ALPHA_EM = 1.0 / 137.035999084
PHI_VEV = 0.25 * math.log(4.0 * math.pi * ALPHA_EM)


class TestRoutePlanckScale(unittest.TestCase):
    """Tests for route_planck_scale()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_escape_routes import route_planck_scale
        result = route_planck_scale()
        self.assertIsInstance(result, dict)

    def test_running_safe(self):
        from alpha_ladder_core.gauge_escape_routes import route_planck_scale
        result = route_planck_scale()
        self.assertTrue(result["running_safe"])

    def test_eotwash_not_viable(self):
        from alpha_ladder_core.gauge_escape_routes import route_planck_scale
        result = route_planck_scale()
        self.assertFalse(result["eotwash_viable"])

    def test_lhc_not_viable(self):
        from alpha_ladder_core.gauge_escape_routes import route_planck_scale
        result = route_planck_scale()
        self.assertFalse(result["lhc_viable"])

    def test_a_0_is_planck_length(self):
        from alpha_ladder_core.gauge_escape_routes import route_planck_scale, L_PL
        result = route_planck_scale()
        self.assertEqual(result["a_0"], L_PL)

    def test_m_kk_very_large(self):
        """KK mass at Planck length should be near Planck scale."""
        from alpha_ladder_core.gauge_escape_routes import route_planck_scale, M_PL_EV
        result = route_planck_scale()
        self.assertGreater(result["m_KK_1"], 1e25)

    def test_R_phys_smaller_than_a_0(self):
        """R = a_0 * exp(phi_vev) < a_0 since phi_vev < 0."""
        from alpha_ladder_core.gauge_escape_routes import route_planck_scale
        result = route_planck_scale()
        self.assertLess(result["R_phys"], result["a_0"])

    def test_a_0_max_safe_positive(self):
        from alpha_ladder_core.gauge_escape_routes import route_planck_scale
        result = route_planck_scale()
        self.assertGreater(result["a_0_max_safe"], 0.0)

    def test_verdict_mentions_boring(self):
        from alpha_ladder_core.gauge_escape_routes import route_planck_scale
        result = route_planck_scale()
        self.assertIn("boring", result["verdict"].lower())

    def test_description_is_string(self):
        from alpha_ladder_core.gauge_escape_routes import route_planck_scale
        result = route_planck_scale()
        self.assertIsInstance(result["description"], str)
        self.assertGreater(len(result["description"]), 100)

    def test_keys_present(self):
        from alpha_ladder_core.gauge_escape_routes import route_planck_scale
        result = route_planck_scale()
        for key in ["a_0", "R_phys", "m_KK_1", "M_6", "running_safe",
                     "eotwash_viable", "a_0_max_safe", "verdict"]:
            self.assertIn(key, result)


class TestRouteBraneMatter(unittest.TestCase):
    """Tests for route_brane_matter()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_escape_routes import route_brane_matter
        result = route_brane_matter()
        self.assertIsInstance(result, dict)

    def test_running_safe(self):
        from alpha_ladder_core.gauge_escape_routes import route_brane_matter
        result = route_brane_matter()
        self.assertTrue(result["running_safe"])

    def test_matching_works(self):
        from alpha_ladder_core.gauge_escape_routes import route_brane_matter
        result = route_brane_matter()
        self.assertTrue(result["matching_works"])

    def test_eotwash_viable(self):
        from alpha_ladder_core.gauge_escape_routes import route_brane_matter
        result = route_brane_matter()
        self.assertTrue(result["eotwash_viable"])

    def test_model_type_braneworld(self):
        from alpha_ladder_core.gauge_escape_routes import route_brane_matter
        result = route_brane_matter()
        self.assertEqual(result["model_type"], "braneworld")

    def test_alpha_kk_matches_alpha_em(self):
        from alpha_ladder_core.gauge_escape_routes import route_brane_matter, ALPHA_EM
        result = route_brane_matter()
        self.assertAlmostEqual(
            result["alpha_kk_verified"], ALPHA_EM, places=12
        )

    def test_verdict_viable(self):
        from alpha_ladder_core.gauge_escape_routes import route_brane_matter
        result = route_brane_matter()
        self.assertIn("VIABLE", result["verdict"])

    def test_bulk_fields_present(self):
        from alpha_ladder_core.gauge_escape_routes import route_brane_matter
        result = route_brane_matter()
        self.assertIsInstance(result["bulk_fields"], list)
        self.assertGreater(len(result["bulk_fields"]), 0)

    def test_brane_fields_present(self):
        from alpha_ladder_core.gauge_escape_routes import route_brane_matter
        result = route_brane_matter()
        self.assertIsInstance(result["brane_fields"], list)
        self.assertGreater(len(result["brane_fields"]), 0)

    def test_inv_alpha_me(self):
        from alpha_ladder_core.gauge_escape_routes import route_brane_matter
        result = route_brane_matter()
        self.assertAlmostEqual(result["inv_alpha_me"], 137.036, places=2)

    def test_inv_alpha_mz_observed(self):
        from alpha_ladder_core.gauge_escape_routes import route_brane_matter
        result = route_brane_matter()
        self.assertAlmostEqual(result["inv_alpha_mz_observed"], 127.944, places=2)


class TestRouteTinyCharge(unittest.TestCase):
    """Tests for route_tiny_charge()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_escape_routes import route_tiny_charge
        result = route_tiny_charge()
        self.assertIsInstance(result, dict)

    def test_matching_fails(self):
        from alpha_ladder_core.gauge_escape_routes import route_tiny_charge
        result = route_tiny_charge()
        self.assertFalse(result["matching_works"])

    def test_running_safe_at_Q_max(self):
        from alpha_ladder_core.gauge_escape_routes import route_tiny_charge
        result = route_tiny_charge()
        self.assertTrue(result["running_safe_at_Q_max"])

    def test_shift_ppm_huge(self):
        """Shift at Q=1 should be enormously large (>10^10 ppm)."""
        from alpha_ladder_core.gauge_escape_routes import route_tiny_charge
        result = route_tiny_charge()
        self.assertGreater(abs(result["shift_ppm_Q1"]), 1e10)

    def test_Q_max_tiny(self):
        """Q_max for safe running is extremely small."""
        from alpha_ladder_core.gauge_escape_routes import route_tiny_charge
        result = route_tiny_charge()
        self.assertLess(result["Q_max_safe"], 1e-5)

    def test_alpha_at_Q_max_much_smaller_than_alpha_em(self):
        from alpha_ladder_core.gauge_escape_routes import route_tiny_charge
        result = route_tiny_charge()
        self.assertLess(result["alpha_at_Q_max"], result["alpha_EM"] * 1e-10)

    def test_log10_ratio_negative(self):
        """Matching fails by many orders -- ratio << 1."""
        from alpha_ladder_core.gauge_escape_routes import route_tiny_charge
        result = route_tiny_charge()
        self.assertLess(result["log10_ratio"], -10)

    def test_fails_by_about_16_orders(self):
        from alpha_ladder_core.gauge_escape_routes import route_tiny_charge
        result = route_tiny_charge()
        self.assertGreater(abs(result["log10_ratio"]), 14)

    def test_verdict_dead(self):
        from alpha_ladder_core.gauge_escape_routes import route_tiny_charge
        result = route_tiny_charge()
        self.assertIn("DEAD", result["verdict"])

    def test_L_max_modes_large(self):
        from alpha_ladder_core.gauge_escape_routes import route_tiny_charge
        result = route_tiny_charge()
        self.assertGreater(result["L_max_modes"], 1000)

    def test_n_modes_large(self):
        from alpha_ladder_core.gauge_escape_routes import route_tiny_charge
        result = route_tiny_charge()
        self.assertGreater(result["n_modes_total"], 1e6)


class TestRouteHeavyBulkMatter(unittest.TestCase):
    """Tests for route_heavy_bulk_matter()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_escape_routes import route_heavy_bulk_matter
        result = route_heavy_bulk_matter()
        self.assertIsInstance(result, dict)

    def test_table_has_5_entries(self):
        from alpha_ladder_core.gauge_escape_routes import route_heavy_bulk_matter
        result = route_heavy_bulk_matter()
        self.assertEqual(len(result["table"]), 5)

    def test_kk_spacing_meV_order(self):
        """KK spacing at a_0=28um should be O(10) meV."""
        from alpha_ladder_core.gauge_escape_routes import route_heavy_bulk_matter
        result = route_heavy_bulk_matter()
        self.assertGreater(result["kk_spacing_meV"], 1.0)
        self.assertLess(result["kk_spacing_meV"], 100.0)

    def test_m_e_not_safe(self):
        """With M_bulk = m_e, running is NOT safe."""
        from alpha_ladder_core.gauge_escape_routes import route_heavy_bulk_matter
        result = route_heavy_bulk_matter()
        me_row = result["table"][0]
        self.assertEqual(me_row["label"], "m_e")
        self.assertFalse(me_row["running_safe"])

    def test_M6_safe(self):
        """With M_bulk = M_6, running IS safe (no modes)."""
        from alpha_ladder_core.gauge_escape_routes import route_heavy_bulk_matter
        result = route_heavy_bulk_matter()
        m6_row = result["table"][4]
        self.assertEqual(m6_row["label"], "M_6")
        self.assertTrue(m6_row["running_safe"])

    def test_M6_not_useful(self):
        """With M_bulk = M_6, the matter is not useful."""
        from alpha_ladder_core.gauge_escape_routes import route_heavy_bulk_matter
        result = route_heavy_bulk_matter()
        m6_row = result["table"][4]
        self.assertFalse(m6_row["useful"])

    def test_verdict_not_real_escape(self):
        from alpha_ladder_core.gauge_escape_routes import route_heavy_bulk_matter
        result = route_heavy_bulk_matter()
        self.assertIn("NOT a real escape", result["verdict"])

    def test_M_6_positive(self):
        from alpha_ladder_core.gauge_escape_routes import route_heavy_bulk_matter
        result = route_heavy_bulk_matter()
        self.assertGreater(result["M_6_eV"], 0.0)

    def test_custom_a_0(self):
        from alpha_ladder_core.gauge_escape_routes import route_heavy_bulk_matter
        result = route_heavy_bulk_matter(a_0_m=1e-6)
        self.assertAlmostEqual(result["a_0_m"], 1e-6, places=12)


class TestRouteComparisonTable(unittest.TestCase):
    """Tests for route_comparison_table()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_escape_routes import route_comparison_table
        result = route_comparison_table()
        self.assertIsInstance(result, dict)

    def test_table_has_4_routes(self):
        from alpha_ladder_core.gauge_escape_routes import route_comparison_table
        result = route_comparison_table()
        self.assertEqual(len(result["table"]), 4)

    def test_brane_is_best(self):
        from alpha_ladder_core.gauge_escape_routes import route_comparison_table
        result = route_comparison_table()
        brane_row = result["table"][1]
        self.assertEqual(brane_row["viable"], "BEST")

    def test_brane_all_true(self):
        from alpha_ladder_core.gauge_escape_routes import route_comparison_table
        result = route_comparison_table()
        brane = result["table"][1]
        self.assertTrue(brane["running_safe"])
        self.assertTrue(brane["matching_works"])
        self.assertTrue(brane["eotwash_open"])

    def test_tiny_q_matching_false(self):
        from alpha_ladder_core.gauge_escape_routes import route_comparison_table
        result = route_comparison_table()
        tiny_q = result["table"][2]
        self.assertFalse(tiny_q["matching_works"])

    def test_description_string(self):
        from alpha_ladder_core.gauge_escape_routes import route_comparison_table
        result = route_comparison_table()
        self.assertIsInstance(result["description"], str)


class TestBraneworldImplications(unittest.TestCase):
    """Tests for braneworld_implications()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_escape_routes import braneworld_implications
        result = braneworld_implications()
        self.assertIsInstance(result, dict)

    def test_not_position_dependent(self):
        from alpha_ladder_core.gauge_escape_routes import braneworld_implications
        result = braneworld_implications()
        self.assertFalse(result["position_dependent"])

    def test_consistency_verified(self):
        from alpha_ladder_core.gauge_escape_routes import braneworld_implications
        result = braneworld_implications()
        self.assertTrue(result["consistency"])

    def test_all_positions_agree(self):
        from alpha_ladder_core.gauge_escape_routes import braneworld_implications
        result = braneworld_implications()
        self.assertTrue(result["all_positions_agree"])

    def test_many_positions_checked(self):
        from alpha_ladder_core.gauge_escape_routes import braneworld_implications
        result = braneworld_implications()
        self.assertEqual(result["brane_positions_checked"], 28)

    def test_alpha_at_brane_matches_alpha_em(self):
        from alpha_ladder_core.gauge_escape_routes import braneworld_implications
        result = braneworld_implications()
        self.assertAlmostEqual(
            result["alpha_at_brane"], result["alpha_EM"], places=12
        )

    def test_bulk_content_keys(self):
        from alpha_ladder_core.gauge_escape_routes import braneworld_implications
        result = braneworld_implications()
        self.assertIn("metric_graviton", result["bulk_content"])
        self.assertIn("dilaton", result["bulk_content"])
        self.assertIn("gauge_fields", result["bulk_content"])

    def test_brane_content_keys(self):
        from alpha_ladder_core.gauge_escape_routes import braneworld_implications
        result = braneworld_implications()
        self.assertIn("fermions", result["brane_content"])
        self.assertIn("higgs", result["brane_content"])


class TestSummarizeEscapeRoutes(unittest.TestCase):
    """Tests for summarize_escape_routes()."""

    def test_returns_dict(self):
        from alpha_ladder_core.gauge_escape_routes import summarize_escape_routes
        result = summarize_escape_routes()
        self.assertIsInstance(result, dict)

    def test_all_routes_present(self):
        from alpha_ladder_core.gauge_escape_routes import summarize_escape_routes
        result = summarize_escape_routes()
        for key in ["route_1_planck", "route_2_brane", "route_3_tiny_q",
                     "route_4_heavy_bulk", "comparison", "braneworld"]:
            self.assertIn(key, result)

    def test_key_messages_list(self):
        from alpha_ladder_core.gauge_escape_routes import summarize_escape_routes
        result = summarize_escape_routes()
        self.assertIsInstance(result["key_messages"], list)
        self.assertEqual(len(result["key_messages"]), 6)

    def test_report_is_long_string(self):
        from alpha_ladder_core.gauge_escape_routes import summarize_escape_routes
        result = summarize_escape_routes()
        self.assertIsInstance(result["report"], str)
        self.assertGreater(len(result["report"]), 1000)

    def test_report_mentions_braneworld(self):
        from alpha_ladder_core.gauge_escape_routes import summarize_escape_routes
        result = summarize_escape_routes()
        self.assertIn("BRANEWORLD", result["report"])

    def test_route2_is_viable(self):
        from alpha_ladder_core.gauge_escape_routes import summarize_escape_routes
        result = summarize_escape_routes()
        self.assertIn("VIABLE", result["route_2_brane"]["verdict"])

    def test_route3_is_dead(self):
        from alpha_ladder_core.gauge_escape_routes import summarize_escape_routes
        result = summarize_escape_routes()
        self.assertIn("DEAD", result["route_3_tiny_q"]["verdict"])


if __name__ == "__main__":
    unittest.main()
