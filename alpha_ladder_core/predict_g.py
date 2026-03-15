"""
Predict Newton's gravitational constant G from the Alpha Ladder.

Bridge formula
--------------
If alpha_G = C * alpha^21, and alpha_G = G * m_e^2 / (hbar * c), then:
    G = C * alpha^21 * hbar * c / m_e^2

Hierarchy formula (no fitted parameters)
-----------------------------------------
alpha_G = alpha^24 * mu^2,  where mu = m_p / m_e, giving:
    G = alpha^24 * mu^2 * hbar * c / m_e^2

Refactored from legacy/predict_G.py.
"""

from decimal import Decimal, getcontext

getcontext().prec = 50


def get_bridge_candidates(constants):
    """Return dict mapping candidate name to Decimal bridge coefficient value.

    The three candidates from the legacy script:
      - phi^2 / 2
      - (5/12) * pi
      - sqrt(e) / cbrt(2)
    """
    phi = constants.phi
    pi = constants.pi
    e = constants.e

    return {
        "φ²/2": phi ** 2 / 2,
        "(5/12) · π": Decimal(5) / Decimal(12) * pi,
        "√e / ³√2": e.sqrt() / Decimal(2) ** (Decimal(1) / Decimal(3)),
    }


def predict_G(bridge_coeff, constants):
    """Predict G given a bridge coefficient and physical constants.

    Formula: G = bridge_coeff * alpha^21 * hbar * c / m_e^2
    """
    alpha = constants.alpha
    hbar = constants.hbar
    c = constants.c
    m_e = constants.m_e

    alpha_21 = alpha ** 21
    alpha_g = bridge_coeff * alpha_21
    return alpha_g * hbar * c / m_e ** 2


def get_G_measurements():
    """Return dict mapping experiment name to (G_value, G_uncertainty) Decimal tuples.

    All 7 measurements from the legacy script.
    """
    return {
        "CODATA 2018 recommended": (Decimal("6.67430e-11"), Decimal("0.00015e-11")),
        "Quinn et al. 2013 (BIPM)": (Decimal("6.67545e-11"), Decimal("0.00018e-11")),
        "Rosi et al. 2014 (atom interf.)": (Decimal("6.67191e-11"), Decimal("0.00099e-11")),
        "Newman et al. 2014": (Decimal("6.67435e-11"), Decimal("0.00013e-11")),
        "Li et al. 2018 (HUST-A)": (Decimal("6.674184e-11"), Decimal("0.000078e-11")),
        "Li et al. 2018 (HUST-B)": (Decimal("6.674484e-11"), Decimal("0.000078e-11")),
        "CODATA 2014 recommended": (Decimal("6.67408e-11"), Decimal("0.00031e-11")),
    }


def compare_prediction(G_pred, G_measurements):
    """Compare a predicted G value against experimental measurements.

    Returns list of dicts with keys: experiment, G_exp, G_unc, sigma, direction.
    """
    results = []
    for exp_name, (G_exp, G_unc) in G_measurements.items():
        diff = G_pred - G_exp
        sigma = float(abs(diff) / G_unc) if G_unc != 0 else float("inf")
        direction = "+" if diff > 0 else "-"
        results.append({
            "experiment": exp_name,
            "G_exp": G_exp,
            "G_unc": G_unc,
            "sigma": sigma,
            "direction": direction,
        })
    return results


def predict_G_hierarchy(constants):
    """Predict G from the hierarchy formula with no fitted parameters.

    Formula: alpha_G = alpha^24 * mu^2,  G = alpha_G * hbar * c / m_e^2
    where mu = m_p / m_e.

    Returns dict with keys: alpha_g, G_predicted, mu, exponent, residual_ppm.
    """
    alpha = constants.alpha
    hbar = constants.hbar
    c = constants.c
    m_e = constants.m_e
    m_p = constants.m_p
    G_measured = constants.G

    mu = m_p / m_e
    exponent = 24
    alpha_g = alpha ** exponent * mu ** 2
    G_predicted = alpha_g * hbar * c / m_e ** 2

    residual_ppm = float(abs(G_predicted - G_measured) / G_measured * Decimal("1e6"))

    return {
        "alpha_g": alpha_g,
        "G_predicted": G_predicted,
        "mu": mu,
        "exponent": exponent,
        "residual_ppm": residual_ppm,
    }


def predict_G_mu_structure(constants):
    """Predict G from the mu-structure formula with no fitted parameters.

    Formula: alpha_G = alpha^24 * mu * (mu - sqrt(phi))
    where mu = m_p / m_e and phi = golden ratio.

    Returns dict with keys: alpha_g, G_predicted, mu, sqrt_phi, residual_ppm.
    """
    alpha = constants.alpha
    hbar = constants.hbar
    c = constants.c
    m_e = constants.m_e
    m_p = constants.m_p
    phi = constants.phi
    G_measured = constants.G

    mu = m_p / m_e
    sqrt_phi = phi.sqrt()
    alpha_g = alpha ** 24 * mu * (mu - sqrt_phi)
    G_predicted = alpha_g * hbar * c / m_e ** 2

    residual_ppm = float((G_predicted - G_measured) / G_measured * Decimal("1e6"))

    return {
        "alpha_g": alpha_g,
        "G_predicted": G_predicted,
        "mu": mu,
        "sqrt_phi": sqrt_phi,
        "residual_ppm": residual_ppm,
    }


def predict_G_mu_structure_refined(constants):
    """Predict G from the refined mu-structure formula with (1-alpha) correction.

    Formula: alpha_G = alpha^24 * mu * (mu - sqrt(phi) * (1 - alpha))
    where mu = m_p / m_e and phi = golden ratio.

    Returns dict with keys: alpha_g, G_predicted, mu, sqrt_phi, k_offset, residual_ppm.
    """
    alpha = constants.alpha
    hbar = constants.hbar
    c = constants.c
    m_e = constants.m_e
    m_p = constants.m_p
    phi = constants.phi
    G_measured = constants.G

    mu = m_p / m_e
    sqrt_phi = phi.sqrt()
    k_offset = sqrt_phi * (1 - alpha)
    alpha_g = alpha ** 24 * mu * (mu - k_offset)
    G_predicted = alpha_g * hbar * c / m_e ** 2

    residual_ppm = float((G_predicted - G_measured) / G_measured * Decimal("1e6"))

    return {
        "alpha_g": alpha_g,
        "G_predicted": G_predicted,
        "mu": float(mu),
        "sqrt_phi": float(sqrt_phi),
        "k_offset": float(k_offset),
        "residual_ppm": residual_ppm,
    }


def summarize_predictions(constants):
    """Summarize predictions for all bridge candidates.

    Returns dict mapping bridge name to {G_pred, avg_sigma, comparisons}.
    """
    candidates = get_bridge_candidates(constants)
    measurements = get_G_measurements()
    summary = {}

    for name, coeff in candidates.items():
        G_pred = predict_G(coeff, constants)
        comparisons = compare_prediction(G_pred, measurements)
        total_sigma = sum(c["sigma"] for c in comparisons)
        avg_sigma = total_sigma / len(comparisons)
        summary[name] = {
            "G_pred": G_pred,
            "avg_sigma": avg_sigma,
            "comparisons": comparisons,
        }

    return summary
