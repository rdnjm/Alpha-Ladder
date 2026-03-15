"""
Single source of truth for all physical and mathematical constants
used across the Alpha Ladder project.

Each CODATA edition stores Decimal values for fundamental constants.
Mathematical constants are shared across all editions.
"""

from decimal import Decimal, getcontext
from types import SimpleNamespace

getcontext().prec = 50

# ---------------------------------------------------------------------------
# Mathematical constants (edition-independent)
# ---------------------------------------------------------------------------
_MATH_CONSTANTS = {
    "pi": Decimal('3.14159265358979323846264338327950288419716939937510'),
    "phi": (1 + Decimal(5).sqrt()) / 2,
    "e": Decimal('2.71828182845904523536028747135266249775724709369995'),
    "ln2": Decimal('0.69314718055994530941723212145817656807550013436026'),
    "sqrt2": Decimal(2).sqrt(),
    "sqrt3": Decimal(3).sqrt(),
    "sqrt5": Decimal(5).sqrt(),
}

# ---------------------------------------------------------------------------
# CODATA editions
# ---------------------------------------------------------------------------
_CODATA_2014_PHYS = {
    "alpha": Decimal('0.0072973525664'),
    "alpha_g": Decimal('1.7512e-45'),
    "m_e": Decimal('9.10938356e-31'),
    "hbar": Decimal('1.0545718e-34'),
    "c": Decimal('299792458'),
    "m_p": Decimal('1.67262192369e-27'),
    "m_mu": Decimal('1.883531627e-28'),
    "k_B": Decimal('1.380649e-23'),
    "e_charge": Decimal('1.602176634e-19'),
    "epsilon_0": Decimal('8.8541878128e-12'),
    "G": Decimal('6.67408e-11'),
    "lambda_bar_c": Decimal('3.8615926764e-13'),
    "a_0": Decimal('5.2917721067e-11'),
    "r_e_nist": Decimal('2.8179403227e-15'),
}

_CODATA_2018_PHYS = {
    "alpha": Decimal('0.0072973525693'),
    "m_e": Decimal('9.1093837015e-31'),
    "hbar": Decimal('1.054571817e-34'),
    "c": Decimal('299792458'),
    "G": Decimal('6.67430e-11'),
    "m_p": Decimal('1.67262192369e-27'),
    "m_mu": Decimal('1.883531627e-28'),
    "k_B": Decimal('1.380649e-23'),
    "e_charge": Decimal('1.602176634e-19'),
    "epsilon_0": Decimal('8.8541878128e-12'),
    "lambda_bar_c": Decimal('3.8615926764e-13'),
    "a_0": Decimal('5.2917721067e-11'),
    "r_e_nist": Decimal('2.8179403227e-15'),
}

# Compute alpha_g from G for CODATA 2018:  alpha_g = G * m_e^2 / (hbar * c)
_CODATA_2018_PHYS["alpha_g"] = (
    _CODATA_2018_PHYS["G"]
    * _CODATA_2018_PHYS["m_e"] ** 2
    / (_CODATA_2018_PHYS["hbar"] * _CODATA_2018_PHYS["c"])
)

_CODATA_2022_PHYS = {
    "alpha": Decimal('0.0072973525643'),
    "m_e": Decimal('9.1093837139e-31'),
    "hbar": Decimal('1.054571817e-34'),       # exact (2019 SI)
    "c": Decimal('299792458'),                 # exact
    "G": Decimal('6.67430e-11'),
    "m_p": Decimal('1.67262192595e-27'),
    "m_mu": Decimal('1.883531627e-28'),
    "k_B": Decimal('1.380649e-23'),            # exact (2019 SI)
    "e_charge": Decimal('1.602176634e-19'),    # exact (2019 SI)
    "epsilon_0": Decimal('8.8541878128e-12'),
    "lambda_bar_c": Decimal('3.8615926764e-13'),
    "a_0": Decimal('5.2917721067e-11'),
    "r_e_nist": Decimal('2.8179403227e-15'),
}

# Compute alpha_g from G for CODATA 2022:  alpha_g = G * m_e^2 / (hbar * c)
_CODATA_2022_PHYS["alpha_g"] = (
    _CODATA_2022_PHYS["G"]
    * _CODATA_2022_PHYS["m_e"] ** 2
    / (_CODATA_2022_PHYS["hbar"] * _CODATA_2022_PHYS["c"])
)

CODATA_EDITIONS = {
    "CODATA 2014": _CODATA_2014_PHYS,
    "CODATA 2018": _CODATA_2018_PHYS,
    "CODATA 2022": _CODATA_2022_PHYS,
}

DEFAULT_EDITION = "CODATA 2022"


def available_editions():
    """Return a list of available CODATA edition names."""
    return list(CODATA_EDITIONS.keys())


def get_constants(edition=DEFAULT_EDITION):
    """Return a SimpleNamespace containing all physical and mathematical
    constants for the requested CODATA edition.

    Parameters
    ----------
    edition : str
        One of the keys in ``CODATA_EDITIONS``.

    Returns
    -------
    SimpleNamespace
        Namespace with attributes for every physical constant in the
        chosen edition plus all mathematical constants.

    Raises
    ------
    ValueError
        If *edition* is not a recognised edition name.
    """
    if edition not in CODATA_EDITIONS:
        raise ValueError(
            f"Unknown edition {edition!r}. "
            f"Available: {available_editions()}"
        )
    combined = {}
    combined.update(CODATA_EDITIONS[edition])
    combined.update(_MATH_CONSTANTS)
    return SimpleNamespace(**combined)


def get_particle_masses():
    """Return dict of particle name to mass in MeV.

    Single source of truth for all particle masses used across
    the Alpha Ladder project. Values from PDG 2024.

    Returns
    -------
    dict
        Mapping particle name -> mass in MeV.
    """
    return {
        "Electron": 0.51099895000,
        "Muon": 105.6583755,
        "Tau": 1776.86,
        "Up quark": 2.16,
        "Down quark": 4.67,
        "Strange quark": 93.4,
        "Charm quark": 1270.0,
        "Bottom quark": 4180.0,
        "Top quark": 172760.0,
        "Proton": 938.272088,
        "W boson": 80377.0,
        "Z boson": 91187.6,
        "Higgs boson": 125100.0,
    }
