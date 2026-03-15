"""Tests for alpha_ladder_core.alpha_units module."""

import math
import sys
import os
import pytest
from decimal import Decimal

# Ensure the project root is on the path so imports work without install.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.alpha_units import (
    get_alpha_units,
    si_to_alpha,
    alpha_to_si,
    express_in_alpha_powers,
)


@pytest.fixture
def constants():
    return get_constants()


@pytest.fixture
def units(constants):
    return get_alpha_units(constants)


# ------------------------------------------------------------------
# get_alpha_units
# ------------------------------------------------------------------

class TestGetAlphaUnits:
    def test_returns_all_base_keys(self, units):
        expected = {"length", "time", "mass", "energy", "charge", "force", "temperature"}
        assert expected == set(units.keys())

    def test_length_is_r_e(self, constants, units):
        alpha = float(constants.alpha)
        a_0 = float(constants.a_0)
        assert units["length"] == pytest.approx(alpha ** 2 * a_0, rel=1e-10)

    def test_time_is_r_e_over_c(self, constants, units):
        c = float(constants.c)
        assert units["time"] == pytest.approx(units["length"] / c, rel=1e-12)

    def test_mass_is_m_e(self, constants, units):
        assert units["mass"] == pytest.approx(float(constants.m_e), rel=1e-12)

    def test_energy_is_m_e_c2(self, constants, units):
        m_e = float(constants.m_e)
        c = float(constants.c)
        assert units["energy"] == pytest.approx(m_e * c ** 2, rel=1e-10)

    def test_charge_is_e(self, constants, units):
        assert units["charge"] == pytest.approx(float(constants.e_charge), rel=1e-12)

    def test_force_is_energy_over_length(self, units):
        assert units["force"] == pytest.approx(units["energy"] / units["length"], rel=1e-10)

    def test_temperature_is_energy_over_kB(self, constants, units):
        k_B = float(constants.k_B)
        assert units["temperature"] == pytest.approx(units["energy"] / k_B, rel=1e-10)

    def test_all_values_positive(self, units):
        for key, val in units.items():
            assert val > 0, f"{key} should be positive, got {val}"


# ------------------------------------------------------------------
# Round-trip conversions: SI -> alpha -> SI
# ------------------------------------------------------------------

class TestRoundTrip:
    @pytest.mark.parametrize("unit_type", [
        "length", "time", "mass", "energy", "charge", "force", "temperature",
    ])
    def test_round_trip_float(self, constants, unit_type):
        original = 1.23456e-10
        alpha_val = si_to_alpha(original, unit_type, constants)
        recovered = alpha_to_si(alpha_val, unit_type, constants)
        assert recovered == pytest.approx(original, rel=1e-12)

    @pytest.mark.parametrize("unit_type", [
        "length", "time", "mass", "energy", "charge", "force", "temperature",
    ])
    def test_round_trip_decimal(self, constants, unit_type):
        original = Decimal("5.6789e-20")
        alpha_val = si_to_alpha(original, unit_type, constants)
        recovered = alpha_to_si(alpha_val, unit_type, constants)
        assert recovered == pytest.approx(float(original), rel=1e-12)


# ------------------------------------------------------------------
# Identity tests: natural quantities equal 1.0 in alpha-units
# ------------------------------------------------------------------

class TestIdentityValues:
    def test_electron_mass_is_one(self, constants):
        m_e = float(constants.m_e)
        assert si_to_alpha(m_e, "mass", constants) == pytest.approx(1.0, rel=1e-12)

    def test_r_e_is_one(self, constants):
        alpha = float(constants.alpha)
        a_0 = float(constants.a_0)
        r_e = alpha ** 2 * a_0
        assert si_to_alpha(r_e, "length", constants) == pytest.approx(1.0, rel=1e-12)

    def test_e_charge_is_one(self, constants):
        e = float(constants.e_charge)
        assert si_to_alpha(e, "charge", constants) == pytest.approx(1.0, rel=1e-12)

    def test_rest_energy_is_one(self, constants):
        m_e = float(constants.m_e)
        c = float(constants.c)
        E = m_e * c ** 2
        assert si_to_alpha(E, "energy", constants) == pytest.approx(1.0, rel=1e-12)

    def test_time_unit_is_one(self, constants):
        alpha = float(constants.alpha)
        a_0 = float(constants.a_0)
        c = float(constants.c)
        r_e = alpha ** 2 * a_0
        t = r_e / c
        assert si_to_alpha(t, "time", constants) == pytest.approx(1.0, rel=1e-12)


# ------------------------------------------------------------------
# Error handling
# ------------------------------------------------------------------

class TestErrors:
    def test_si_to_alpha_bad_unit(self, constants):
        with pytest.raises(ValueError, match="Unknown unit_type"):
            si_to_alpha(1.0, "luminosity", constants)

    def test_alpha_to_si_bad_unit(self, constants):
        with pytest.raises(ValueError, match="Unknown unit_type"):
            alpha_to_si(1.0, "luminosity", constants)

    def test_express_non_positive_raises(self, constants):
        with pytest.raises(ValueError, match="non-positive"):
            express_in_alpha_powers(-1.0, "mass", constants)


# ------------------------------------------------------------------
# express_in_alpha_powers
# ------------------------------------------------------------------

class TestExpressInAlphaPowers:
    def test_alpha_unit_value_matches_si_to_alpha(self, constants):
        si_val = float(constants.m_e)
        result = express_in_alpha_powers(si_val, "mass", constants)
        expected = si_to_alpha(si_val, "mass", constants)
        assert result["alpha_unit_value"] == pytest.approx(expected, rel=1e-12)

    def test_electron_mass_power_is_zero(self, constants):
        """m_e in alpha mass units is 1.0, and alpha^0 = 1, so n = 0."""
        m_e = float(constants.m_e)
        result = express_in_alpha_powers(m_e, "mass", constants)
        assert result["nearest_power"] == 0
        assert result["exact_power"] == pytest.approx(0.0, abs=1e-10)
        assert result["residual_pct"] < 1e-6

    def test_known_power_for_alpha_squared_length(self, constants):
        """If we feed r_e * alpha^3 as an SI length, the alpha-unit value
        should be alpha^3, giving exact_power ~ 3."""
        alpha = float(constants.alpha)
        a_0 = float(constants.a_0)
        r_e = alpha ** 2 * a_0
        si_val = r_e * alpha ** 3  # should be alpha^3 in alpha length units
        result = express_in_alpha_powers(si_val, "length", constants)
        assert result["nearest_power"] == 3
        assert result["exact_power"] == pytest.approx(3.0, abs=1e-8)
        assert result["residual_pct"] < 1e-6

    def test_exact_power_is_consistent(self, constants):
        """alpha ^ exact_power should reproduce the alpha-unit value."""
        alpha = float(constants.alpha)
        si_val = 1e-20
        result = express_in_alpha_powers(si_val, "length", constants)
        reconstructed = alpha ** result["exact_power"]
        assert reconstructed == pytest.approx(
            result["alpha_unit_value"], rel=1e-8
        )

    def test_residual_pct_is_non_negative(self, constants):
        result = express_in_alpha_powers(1e-15, "length", constants)
        assert result["residual_pct"] >= 0.0
