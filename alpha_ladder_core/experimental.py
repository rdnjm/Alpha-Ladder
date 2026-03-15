"""
Experimental strategy for testing the Alpha Ladder.
Refactored from legacy/experimental_strategy.py.
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50


def _get_alpha(constants):
    """Extract alpha as a float from either a SimpleNamespace or dict."""
    return float(constants.alpha) if hasattr(constants, 'alpha') else constants["alpha"]


def strategy_mass_ratios(constants) -> list:
    """
    Strategy 1: Do particle mass ratios fall on alpha-power rungs?
    Returns list of dicts with ratio analysis for each mass pair.
    """
    alpha_f = _get_alpha(constants)
    alpha_d = Decimal(str(alpha_f))

    m_e = constants.m_e
    m_p = constants.m_p
    m_mu = constants.m_mu

    ratios = {
        "m_p / m_e":  float(m_p / m_e),
        "m_mu / m_e": float(m_mu / m_e),
        "m_p / m_mu": float(m_p / m_mu),
    }

    results = []
    for name, ratio in ratios.items():
        n_exact = math.log(ratio) / math.log(1 / alpha_f)
        n_round = round(n_exact)
        coeff = ratio * alpha_f ** n_round

        nearby = []
        for n in [n_round - 1, n_round, n_round + 1]:
            coeff_n = ratio * alpha_f ** n
            nearby.append({"n": -n, "coeff": coeff_n})

        results.append({
            "name": name,
            "ratio": ratio,
            "n_exact": n_exact,
            "n_round": n_round,
            "coeff_at_round": coeff,
            "nearby_rungs": nearby,
        })

    return results


def strategy_multiple_paths(constants) -> dict:
    """
    Strategy 2: Multiple independent paths to alpha_G.
    Returns dict with Path A, B, C analyses.
    """
    alpha_f = _get_alpha(constants)
    alpha_d = Decimal(str(alpha_f))

    phi_d = constants.phi
    hbar = constants.hbar
    c = constants.c
    m_e = constants.m_e
    m_p = constants.m_p
    G = constants.G

    alpha_21 = alpha_d ** 21

    # Path A: phi^2/2 * alpha^21
    alpha_g_A = (phi_d ** 2 / Decimal(2)) * alpha_21
    G_A = alpha_g_A * hbar * c / m_e ** 2

    # Path B: proton-based alpha_G
    alpha_g_proton = G * m_p ** 2 / (hbar * c)
    n_proton = math.log(float(alpha_g_proton)) / math.log(float(alpha_d))
    proton_rung = round(n_proton)
    proton_coeff = float(alpha_g_proton) / float(alpha_d) ** proton_rung

    # Path C: Planck mass
    hbar_c_over_G = hbar * c / G
    m_Pl = hbar_c_over_G.sqrt()
    ratio_e_Pl = m_e / m_Pl
    n_planck = math.log(float(ratio_e_Pl)) / math.log(float(alpha_d))

    # Path D: alpha^24 * mu^2 (zero free parameters)
    mu = m_p / m_e
    alpha_g_D = alpha_d ** 24 * mu ** 2
    G_D = alpha_g_D * hbar * c / m_e ** 2

    return {
        "path_A": {
            "description": "phi^2/2 * alpha^21",
            "alpha_G": float(alpha_g_A),
            "G_predicted": float(G_A),
        },
        "path_B": {
            "description": "proton-based alpha_G",
            "alpha_G_proton": float(alpha_g_proton),
            "n_proton": n_proton,
            "nearest_rung": proton_rung,
            "coeff_at_rung": proton_coeff,
        },
        "path_C": {
            "description": "Planck mass relationship",
            "m_e_over_m_Pl": float(ratio_e_Pl),
            "alpha_G_from_planck": float(ratio_e_Pl ** 2),
            "n_planck": n_planck,
            "n_alpha_G": 2 * n_planck,
        },
        "path_D": {
            "description": "alpha^24 * mu^2 (zero free parameters)",
            "alpha_G": float(alpha_g_D),
            "G_predicted": float(G_D),
            "residual_ppm": float((G_D - G) / G * Decimal('1000000')),
        },
    }


def strategy_dark_sector(constants) -> dict:
    """
    Strategy 3: The dark matter coupling (alpha^10).
    Returns dict with dark photon coupling predictions and experimental bounds.
    """
    alpha_f = _get_alpha(constants)
    alpha_10 = alpha_f ** 10
    epsilon = math.sqrt(alpha_10)

    bounds = {
        "BaBar_2017": {"epsilon_max": 1e-3, "mass_range": "1-10 GeV"},
        "NA64_2019": {"epsilon_max": 1e-4, "mass_range": "1-100 MeV"},
        "LDMX_projected": {"epsilon_max": 1e-6, "mass_range": "1-100 MeV"},
    }

    below_current = epsilon < 1e-3
    if below_current:
        orders_below = math.log10(1e-3 / epsilon)
    else:
        orders_below = None

    return {
        "alpha_10": alpha_10,
        "epsilon_predicted": epsilon,
        "experimental_bounds": bounds,
        "below_current_bounds": below_current,
        "orders_below_sensitivity": orders_below,
    }


def strategy_muon_g2(constants) -> dict:
    """
    Strategy 4: Muon g-2 anomaly analysis.
    Returns dict with anomaly analysis.
    """
    alpha_f = _get_alpha(constants)
    alpha_d = Decimal(str(alpha_f))

    delta_a_mu = Decimal('2.51e-9')
    a_mu_exp = Decimal('1.16592061e-3')

    n_anomaly = math.log(float(delta_a_mu)) / math.log(float(alpha_d))
    nearest_rung_anomaly = round(n_anomaly)
    anomaly_coeff = float(delta_a_mu) / float(alpha_d) ** nearest_rung_anomaly

    frac_anomaly = float(delta_a_mu / a_mu_exp)
    n_frac = math.log(frac_anomaly) / math.log(float(alpha_d))
    nearest_rung_frac = round(n_frac)

    return {
        "delta_a_mu": float(delta_a_mu),
        "a_mu_exp": float(a_mu_exp),
        "fractional_anomaly": frac_anomaly,
        "n_anomaly": n_anomaly,
        "nearest_rung_anomaly": nearest_rung_anomaly,
        "anomaly_coeff": anomaly_coeff,
        "n_fractional": n_frac,
        "nearest_rung_fractional": nearest_rung_frac,
    }


def strategy_experimental_approaches(constants) -> dict:
    """
    Strategy 5: Accessible experimental approaches to G.
    Returns dict with experimental approach descriptions and data.
    """
    alpha_d = Decimal(str(_get_alpha(constants)))
    phi = constants.phi
    hbar = constants.hbar
    c = constants.c
    m_e = constants.m_e
    m_p = constants.m_p
    G_measured = constants.G

    # Bridge prediction
    alpha_g_bridge = (phi ** 2 / 2) * alpha_d ** 21
    G_bridge = float(alpha_g_bridge * hbar * c / m_e ** 2)

    # Hierarchy prediction
    mu = m_p / m_e
    alpha_g_hier = alpha_d ** 24 * mu ** 2
    G_hierarchy = float(alpha_g_hier * hbar * c / m_e ** 2)

    G_meas = float(G_measured)
    diff_bridge_ppm = abs(G_bridge - G_meas) / G_meas * 1e6
    diff_hier_ppm = abs(G_hierarchy - G_meas) / G_meas * 1e6

    return {
        "ladder_prediction": G_bridge,
        "hierarchy_prediction": G_hierarchy,
        "codata_2018": {"value": G_meas, "uncertainty": 0.00015e-11},
        "difference": abs(G_bridge - G_meas),
        "difference_ppm": round(diff_bridge_ppm),
        "hierarchy_difference_ppm": round(diff_hier_ppm),
        "required_precision_ppm": 5,
        "approaches": {
            "A_atom_interferometry": {
                "description": "Atom interferometry (most promising)",
                "rosi_2014_value": 6.67191e-11,
                "technique": "Cold atom clouds as test masses",
                "projected_precision_ppm": "5-10",
                "groups": ["Stanford (Kasevich)", "Florence (Tino)", "Wuhan (Zhu)"],
            },
            "B_MEMS_oscillators": {
                "description": "MEMS oscillators",
                "technique": "Microscale torsion oscillators on silicon chips",
                "current_precision_ppm": 100,
            },
            "C_targeted_reanalysis": {
                "description": "Targeted re-analysis of Rosi/Tino data",
                "rosi_uncertainty_ppm": 99,
                "prediction_sigma_from_rosi": 1.3,
                "note": "Only existing measurement consistent with the ladder",
            },
        },
    }


def get_all_strategies(constants) -> list:
    """
    Run all 5 strategies and return results as a list.
    """
    return [
        {"strategy": "mass_ratios", "data": strategy_mass_ratios(constants)},
        {"strategy": "multiple_paths", "data": strategy_multiple_paths(constants)},
        {"strategy": "dark_sector", "data": strategy_dark_sector(constants)},
        {"strategy": "muon_g2", "data": strategy_muon_g2(constants)},
        {"strategy": "experimental_approaches", "data": strategy_experimental_approaches(constants)},
    ]
