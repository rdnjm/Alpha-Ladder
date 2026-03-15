"""
Investigate the bridge coefficient Phi = alpha_G / alpha^21.

Refactored from legacy/phi_capital.py.
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50


def compute_bridge_coefficient(constants):
    """Compute Phi = alpha_g / alpha^21 as a Decimal."""
    alpha = constants.alpha
    alpha_g = constants.alpha_g
    alpha_21 = alpha ** 21
    return alpha_g / alpha_21


def compute_sensitivity(constants):
    """Compute sensitivity of Phi to G uncertainty.

    Uses G +/- 1 sigma to compute Phi range.
    Returns dict with keys: Phi_central, Phi_low, Phi_high, uncertainty.
    """
    alpha = constants.alpha
    m_e = constants.m_e
    hbar = constants.hbar
    c = constants.c
    alpha_21 = alpha ** 21

    G_central = Decimal("6.67430e-11")
    G_low = Decimal("6.67415e-11")       # -1 sigma
    G_high = Decimal("6.67445e-11")      # +1 sigma

    def alpha_g_from_G(G):
        return G * m_e ** 2 / (hbar * c)

    Phi_central = alpha_g_from_G(G_central) / alpha_21
    Phi_low = alpha_g_from_G(G_low) / alpha_21
    Phi_high = alpha_g_from_G(G_high) / alpha_21

    return {
        "Phi_central": Phi_central,
        "Phi_low": Phi_low,
        "Phi_high": Phi_high,
        "uncertainty": (Phi_high - Phi_low) / 2,
    }


def get_phi_candidates(constants):
    """Return list of (err_pct, name, value) tuples sorted by error.

    Includes all candidates from the legacy phi_capital.py script.
    """
    alpha = constants.alpha
    alpha_g = constants.alpha_g

    pi = Decimal("3.14159265358979323846264338327950288419716939937510")
    phi = (1 + Decimal(5).sqrt()) / 2
    e = Decimal("2.71828182845904523536028747135266249775724709369995")

    alpha_21 = alpha ** 21
    Phi = alpha_g / alpha_21
    Phi_f = float(Phi)

    candidates = {
        "φ²/2 = (φ+1)/2": float(phi ** 2 / 2),
        "(5/12) · π": 5.0 / 12.0 * float(pi),
        "√e / ³√2": float(e.sqrt() / Decimal(2) ** (Decimal(1) / Decimal(3))),
        "4/3  (exact)": 4.0 / 3.0,
        "ln(φ) + 1": float(phi.ln() + 1),
        "1 + 1/π": 1.0 + 1.0 / float(pi),
        "√φ · φ⁻¹ᐟ⁴": float(phi) ** 0.5 * float(phi) ** (-0.25),
        "2φ − 2": 2.0 * float(phi) - 2.0,
        "φ / √(φ+1)": float(phi / (phi + 1).sqrt()),
        "e / (e − φ/π)": float(e) / (float(e) - float(phi) / float(pi)),
        "7/4 − φ/4": 7.0 / 4.0 - float(phi) / 4.0,
        "1 + π/e²": 1.0 + float(pi) / float(e) ** 2,
        "sec(1 radian)": 1.0 / math.cos(1.0),
        "2·cos(π/5)": 2.0 * math.cos(math.pi / 5),
        "√φ + 1/φ²": math.sqrt(float(phi)) + 1.0 / float(phi) ** 2,
        "1 + Catalan/π": 1.0 + 0.9159655941 / math.pi,
    }

    ranked = []
    for name, val in candidates.items():
        err = abs(val - Phi_f) / Phi_f * 100
        ranked.append((err, name, val))

    ranked.sort()
    return ranked


def check_error_bar_containment(constants):
    """Check which of the top 10 candidates fall within G's uncertainty window.

    Returns list of (name, value, inside_bool) for top 10 candidates.
    """
    sensitivity = compute_sensitivity(constants)
    Phi_low = float(sensitivity["Phi_low"])
    Phi_high = float(sensitivity["Phi_high"])

    ranked = get_phi_candidates(constants)

    results = []
    for err, name, val in ranked[:10]:
        inside = Phi_low <= val <= Phi_high
        results.append((name, val, inside))

    return results
