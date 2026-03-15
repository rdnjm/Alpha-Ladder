"""
Dilaton interpretation of the gravity rung (alpha^21).
Refactored from legacy/dilaton_check.py.
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50


def decompose_21(constants: dict) -> dict:
    """
    Decompose alpha^21 = alpha^20 * alpha^1 and compute alternative splits.
    Returns dict with alpha_20, alpha_1, bridge values and all alternative splits.
    """
    alpha_val = float(constants.alpha) if hasattr(constants, 'alpha') else constants["alpha"]
    alpha = Decimal(str(alpha_val))
    phi = (1 + Decimal(5).sqrt()) / Decimal(2)

    alpha_20 = alpha ** 20
    alpha_1 = alpha
    bridge = phi ** 2 / Decimal(2)

    # Alternative splits from legacy (identical)
    splits = {
        (20, 1): "Riemann(4D) + scalar",
        (19, 2): "???",
        (18, 3): "???",
        (15, 6): "Riemann(3D)=6 + 15 ???",
        (14, 7): "G2 Lie algebra + 7 ???",
        (12, 9): "SM gauge bosons + 9 ???",
        (10, 11): "Ricci(4D) or Weyl(4D) + 11",
        (6, 15):  "metric(3D) + Riemann(3D)=6 ???",
    }

    alternative_splits = []
    for (a, b), meaning in splits.items():
        alternative_splits.append({
            "a": a,
            "b": b,
            "alpha_a": float(alpha ** a),
            "alpha_b": float(alpha ** b),
            "meaning": meaning,
        })

    return {
        "alpha_20": float(alpha_20),
        "alpha_1": float(alpha_1),
        "bridge": float(bridge),
        "alpha_21": float(alpha ** 21),
        "alternative_splits": alternative_splits,
    }


def compute_bd_parameter(constants: dict) -> dict:
    """
    Compute the Brans-Dicke parameter implied by the phi^2/2 bridge coefficient.
    Returns dict with omega, bare_ratio, dilaton_mass info.
    IDENTICAL formulas to legacy.
    """
    alpha_val = float(constants.alpha) if hasattr(constants, 'alpha') else constants["alpha"]
    alpha = Decimal(str(alpha_val))
    phi_d = (1 + Decimal(5).sqrt()) / Decimal(2)
    phi_f = float(phi_d)

    bridge = phi_f ** 2 / 2.0
    phi2 = phi_f ** 2

    # BD parameter: (2w+4)/(2w+3) = phi^2/2
    # => w = (3*phi^2/2 - 4) / (2 - phi^2)
    numerator = 3 * phi2 / 2 - 4
    denominator = 2 - phi2

    omega = None
    if denominator != 0:
        omega = numerator / denominator

    # Minimum dilaton mass to evade Cassini bound
    hbar_f = 1.054571817e-34
    c_f = 2.99792458e8
    AU = 1.496e11
    cassini_limit = 2.3e-5

    # Fifth-force coupling: alpha_fifth = 2 / (2*omega + 3)
    denom_fifth = 2.0 * omega + 3.0 if omega is not None else 3.0
    alpha_fifth = 2.0 / denom_fifth if denom_fifth != 0 else 0.0

    # Required: alpha_fifth * exp(-AU / lambda_max) < cassini_limit
    # => lambda_max = AU / ln(alpha_fifth / cassini_limit)
    if alpha_fifth > cassini_limit:
        suppression = math.log(alpha_fifth / cassini_limit)
        lambda_max = AU / suppression
    else:
        suppression = 0.0
        lambda_max = AU  # fallback: already passes

    m_d_min = hbar_f / (c_f * lambda_max)
    m_d_min_eV = m_d_min * c_f**2 / 1.602176634e-19

    # Dark/dilaton mass scale from the ladder
    m_Pl_eV = 1.22089e28  # eV
    m_dark_scale = float(alpha) ** 5 * m_Pl_eV

    return {
        "omega": omega,
        "bare_ratio": bridge,
        "phi_squared_over_2": bridge,
        "omega_negative": omega is not None and omega < 0,
        "cassini_bound": 40000,
        "omega_excluded_massless": omega is not None and abs(omega) < 40000,
        "dilaton_mass_min_kg": m_d_min,
        "dilaton_mass_min_eV": m_d_min_eV,
        "cassini_lambda_max_m": lambda_max,
        "cassini_suppression_factor": suppression,
        "dark_scale_eV": m_dark_scale,
    }


def compute_casimir_dilaton_mass(constants=None):
    """
    Compute the dilaton mass from Casimir stabilization on S^2.

    This wraps the casimir_stabilization module to provide the dilaton
    mass prediction alongside the Cassini lower bound from compute_bd_parameter().

    Parameters
    ----------
    constants : SimpleNamespace or None
        Physical constants (needed for compute_bd_parameter comparison).

    Returns
    -------
    dict with:
        casimir_mass_eV : float or None (from Casimir stabilization)
        cassini_lower_bound_eV : float (from Cassini constraint)
        casimir_available : bool
        casimir_minimum_exists : bool
        first_principles : bool
        comparison : str (text comparing the two)
        honest_assessment : str
    """
    # Cassini bound
    cassini_mass = None
    if constants is not None:
        bd = compute_bd_parameter(constants)
        cassini_mass = bd.get("dilaton_mass_min_eV")

    # Casimir prediction
    casimir_available = False
    casimir_mass = None
    casimir_min_exists = False
    honest = "Casimir stabilization module not available."
    try:
        from alpha_ladder_core.casimir_stabilization import compute_dilaton_mass_casimir
        cas = compute_dilaton_mass_casimir()
        casimir_available = True
        casimir_min_exists = cas.get("minimum_exists", False)
        casimir_mass = cas.get("m_phi_eV")
        honest = cas.get("honest_assessment", "")
    except ImportError:
        pass

    comparison = ""
    if casimir_mass is not None and cassini_mass is not None:
        if casimir_mass > cassini_mass:
            comparison = f"Casimir mass ({casimir_mass:.2e} eV) exceeds Cassini lower bound ({cassini_mass:.2e} eV)."
        else:
            comparison = f"Casimir mass ({casimir_mass:.2e} eV) below Cassini lower bound ({cassini_mass:.2e} eV)."
    elif casimir_mass is None:
        comparison = "Casimir mass not computable (no stable minimum). Cassini lower bound still applies."

    return {
        "casimir_mass_eV": casimir_mass,
        "cassini_lower_bound_eV": cassini_mass,
        "casimir_available": casimir_available,
        "casimir_minimum_exists": casimir_min_exists,
        "first_principles": casimir_min_exists and casimir_mass is not None,
        "comparison": comparison,
        "honest_assessment": honest,
    }


def reconcile_6d(constants: dict) -> dict:
    """
    Reconcile the 20+1 decomposition with 6D metric.
    Returns dict with the 6D metric decomposition info.
    """
    # 6D symmetric rank-2 tensor: 6*7/2 = 21 components
    # Decomposes as: g_uv(4D,10) + g_ab(2D,3) + g_ua(mixed,8) = 21
    g_uv_4d = 10   # 4D metric components (symmetric 4x4)
    g_ab_2d = 3    # 2D internal metric (symmetric 2x2)
    g_ua_mixed = 8  # mixed components (4 x 2)
    total = g_uv_4d + g_ab_2d + g_ua_mixed

    # Riemann tensor in 4D has 20 independent components
    riemann_4d = 20
    scalar_dilaton = 1
    riemann_plus_scalar = riemann_4d + scalar_dilaton

    return {
        "g_uv_4d": g_uv_4d,
        "g_ab_2d": g_ab_2d,
        "g_ua_mixed": g_ua_mixed,
        "total_6d_metric": total,
        "riemann_4d": riemann_4d,
        "scalar_dilaton": scalar_dilaton,
        "riemann_plus_scalar": riemann_plus_scalar,
        "views": {
            "bottom_up": "4D gravity needs 20+1 to couple to the SM",
            "top_down": "6D gravity has 21 metric components naturally",
        },
    }
