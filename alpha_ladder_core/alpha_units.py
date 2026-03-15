"""
Conversions between SI units and alpha-natural units.

In the alpha-natural unit system the fine-structure constant is the
fundamental scaling parameter.  Base units are defined in terms of
alpha and electron properties:

    length      classical electron radius  r_e = alpha^2 * a_0
    time        r_e / c
    mass        m_e  (electron mass)
    energy      m_e * c^2  (electron rest energy)
    charge      e  (elementary charge)
    force       m_e * c^2 / r_e
    temperature m_e * c^2 / k_B  (if k_B available)
"""

import math
from decimal import Decimal


# ------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------

def _to_float(v):
    """Coerce Decimal or float to float."""
    return float(v) if isinstance(v, Decimal) else v


# ------------------------------------------------------------------
# public API
# ------------------------------------------------------------------

def get_alpha_units(constants):
    """Return the base alpha-unit system.

    Parameters
    ----------
    constants : SimpleNamespace
        Object returned by ``get_constants()``.

    Returns
    -------
    dict
        Mapping of unit type names to their SI equivalents:

        * ``length``  -- r_e in metres
        * ``time``    -- r_e / c in seconds
        * ``mass``    -- m_e in kg
        * ``energy``  -- m_e * c**2 in Joules
        * ``charge``  -- e in Coulombs
        * ``force``   -- m_e * c**2 / r_e in Newtons
        * ``temperature`` -- m_e * c**2 / k_B in Kelvin (only when
          ``k_B`` is present in *constants*)
    """
    alpha = _to_float(constants.alpha)
    a_0 = _to_float(constants.a_0)
    c = _to_float(constants.c)
    m_e = _to_float(constants.m_e)
    e_charge = _to_float(constants.e_charge)

    r_e = alpha ** 2 * a_0          # classical electron radius
    t_unit = r_e / c                # time unit
    energy = m_e * c ** 2           # electron rest energy
    force = energy / r_e            # force unit

    units = {
        "length": r_e,
        "time": t_unit,
        "mass": m_e,
        "energy": energy,
        "charge": e_charge,
        "force": force,
    }

    if hasattr(constants, "k_B") and constants.k_B is not None:
        k_B = _to_float(constants.k_B)
        units["temperature"] = energy / k_B

    return units


def si_to_alpha(value, unit_type, constants):
    """Convert an SI value to alpha-units.

    Parameters
    ----------
    value : float or Decimal
        The value expressed in SI units.
    unit_type : str
        One of ``'length'``, ``'time'``, ``'mass'``, ``'energy'``,
        ``'charge'``, ``'force'``, ``'temperature'``.
    constants : SimpleNamespace
        Object returned by ``get_constants()``.

    Returns
    -------
    float
        The equivalent value in alpha-units.

    Raises
    ------
    ValueError
        If *unit_type* is not recognised.
    """
    units = get_alpha_units(constants)
    if unit_type not in units:
        raise ValueError(
            f"Unknown unit_type {unit_type!r}. "
            f"Available: {sorted(units)}"
        )
    return _to_float(value) / units[unit_type]


def alpha_to_si(value, unit_type, constants):
    """Convert an alpha-unit value back to SI.

    Parameters
    ----------
    value : float or Decimal
        The value expressed in alpha-units.
    unit_type : str
        One of ``'length'``, ``'time'``, ``'mass'``, ``'energy'``,
        ``'charge'``, ``'force'``, ``'temperature'``.
    constants : SimpleNamespace
        Object returned by ``get_constants()``.

    Returns
    -------
    float
        The equivalent value in SI units.

    Raises
    ------
    ValueError
        If *unit_type* is not recognised.
    """
    units = get_alpha_units(constants)
    if unit_type not in units:
        raise ValueError(
            f"Unknown unit_type {unit_type!r}. "
            f"Available: {sorted(units)}"
        )
    return _to_float(value) * units[unit_type]


def express_in_alpha_powers(si_value, unit_type, constants):
    """Express an SI quantity as alpha**n in alpha-units.

    First the *si_value* is converted to alpha-units, then the result
    is decomposed into powers of the fine-structure constant.

    Parameters
    ----------
    si_value : float or Decimal
        The SI quantity.
    unit_type : str
        Unit type string (see ``si_to_alpha``).
    constants : SimpleNamespace
        Object returned by ``get_constants()``.

    Returns
    -------
    dict
        * ``alpha_unit_value`` -- the value in alpha units
        * ``nearest_power``    -- closest integer *n* where
          ``alpha**n`` approximates the value
        * ``exact_power``      -- float *n* where ``alpha**n == value``
          exactly
        * ``residual_pct``     -- percentage deviation between
          ``alpha**nearest_power`` and the alpha-unit value
    """
    alpha = _to_float(constants.alpha)
    alpha_val = si_to_alpha(si_value, unit_type, constants)

    if alpha_val <= 0:
        raise ValueError(
            "Cannot express a non-positive value as a power of alpha."
        )

    exact_power = math.log(alpha_val) / math.log(alpha)
    nearest_power = round(exact_power)
    approx = alpha ** nearest_power
    residual_pct = abs(approx - alpha_val) / abs(alpha_val) * 100.0

    return {
        "alpha_unit_value": alpha_val,
        "nearest_power": nearest_power,
        "exact_power": exact_power,
        "residual_pct": residual_pct,
    }
