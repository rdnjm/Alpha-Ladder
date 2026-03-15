"""
Map the residual delta = sqrt(phi) - k_exact against Standard Model constants.

The mu-structure formula predicts G via alpha_g = alpha^24 * mu * (mu - sqrt(phi)),
achieving -5.37 ppm. The exact offset k_exact that would give zero residual differs
from sqrt(phi) by delta ~ 0.00985. This module searches for closed-form expressions
for delta in terms of known physical and mathematical constants.
"""

from decimal import Decimal, getcontext
import math

getcontext().prec = 50

from alpha_ladder_core.constants import get_constants


def compute_delta(constants=None):
    """Compute the exact residual delta = sqrt(phi) - k_exact.

    k_exact is the offset such that alpha_g = alpha^24 * mu * (mu - k_exact)
    reproduces the measured alpha_g exactly.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        k_exact, sqrt_phi, delta, delta_over_alpha, mu, alpha_g_measured.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    hbar = constants.hbar
    c = constants.c
    m_e = constants.m_e
    m_p = constants.m_p
    phi = constants.phi
    G = constants.G

    mu = m_p / m_e
    alpha_g_measured = G * m_e ** 2 / (hbar * c)
    alpha_24 = alpha ** 24

    # alpha_g = alpha^24 * mu * (mu - k) => k = mu - alpha_g / (alpha^24 * mu)
    k_exact = mu - alpha_g_measured / (alpha_24 * mu)
    sqrt_phi = phi.sqrt()
    delta = sqrt_phi - k_exact

    return {
        "k_exact": float(k_exact),
        "sqrt_phi": float(sqrt_phi),
        "delta": float(delta),
        "delta_over_alpha": float(delta / alpha),
        "mu": float(mu),
        "alpha_g_measured": float(alpha_g_measured),
    }


def scan_anomalous_magnetic_moment(constants=None):
    """Map delta against QED anomalous magnetic moment coefficients.

    Tests delta against Schwinger term alpha/(2*pi), and higher-order
    g-2 coefficients C2 through C4.

    Parameters
    ----------
    constants : SimpleNamespace or None

    Returns
    -------
    dict
        candidates : list of dict with keys name, expression, value, ratio, closeness.
        best_match : dict or None.
        assessment : str.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    phi = constants.phi
    pi = constants.pi

    delta_info = compute_delta(constants)
    delta = Decimal(str(delta_info["delta"]))
    alpha_d = alpha  # Decimal

    # QED g-2 coefficients (numerical values)
    schwinger = alpha_d / (2 * pi)  # ~ 0.00116
    c2_coeff = Decimal("-0.328478965579") * (alpha_d / pi) ** 2  # 2nd order
    c3_coeff = Decimal("1.181241456") * (alpha_d / pi) ** 3  # 3rd order
    c4_coeff = Decimal("-1.9113") * (alpha_d / pi) ** 4  # 4th order (approx)

    # Also test products with other constants
    candidates = []

    terms = [
        ("alpha/(2*pi) [Schwinger]", schwinger),
        ("(alpha/pi)^2 * 0.328", abs(c2_coeff)),
        ("(alpha/pi)^3 * 1.181", abs(c3_coeff)),
        ("(alpha/pi)^4 * 1.911", abs(c4_coeff)),
        ("alpha^2", alpha_d ** 2),
        ("3*alpha^2", 3 * alpha_d ** 2),
        ("alpha^2 * pi", alpha_d ** 2 * pi),
        ("alpha^2 * pi^2", alpha_d ** 2 * pi ** 2),
        ("alpha * sqrt(phi)", alpha_d * phi.sqrt()),
        ("alpha * phi", alpha_d * phi),
        ("alpha * (phi - 1)", alpha_d * (phi - 1)),
        ("alpha^2 * mu", alpha_d ** 2 * Decimal(str(delta_info["mu"]))),
        ("alpha / (2*pi) * sqrt(phi)", schwinger * phi.sqrt()),
        ("alpha^2 * (1 + alpha/(2*pi))", alpha_d ** 2 * (1 + schwinger)),
    ]

    for name, value in terms:
        if value == 0:
            continue
        ratio = float(delta / value)
        closeness = abs(ratio - round(ratio))
        candidates.append({
            "name": name,
            "expression": name,
            "value": float(value),
            "ratio": ratio,
            "nearest_integer": round(ratio) if abs(ratio) < 1000 else None,
            "closeness_to_integer": closeness if abs(ratio) < 1000 else float("inf"),
        })

    # Sort by closeness to an integer ratio
    valid = [c for c in candidates if c["nearest_integer"] is not None and c["nearest_integer"] != 0]
    valid.sort(key=lambda c: c["closeness_to_integer"])

    best = valid[0] if valid else None

    return {
        "candidates": candidates,
        "best_match": best,
        "assessment": (
            f"Best g-2 match: delta / ({best['expression']}) = {best['ratio']:.4f} "
            f"(nearest integer: {best['nearest_integer']}, "
            f"fractional part: {best['closeness_to_integer']:.4f})"
            if best else "No close integer ratios found among g-2 terms."
        ),
    }


def scan_zeta_functions(constants=None):
    """Map delta against Riemann zeta function values and related constants.

    Tests zeta(2) = pi^2/6, zeta(3) = Apery's constant, zeta(1/2),
    Euler-Mascheroni gamma, Catalan's constant, etc.

    Parameters
    ----------
    constants : SimpleNamespace or None

    Returns
    -------
    dict
        candidates : list of dict.
        best_match : dict or None.
        assessment : str.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    pi = constants.pi

    delta_info = compute_delta(constants)
    delta = Decimal(str(delta_info["delta"]))
    alpha_d = alpha

    # Zeta and related constants (as Decimal for precision)
    zeta_2 = pi ** 2 / 6  # 1.6449...
    zeta_3 = Decimal("1.2020569031595942853997381615114499907649862923404988817922715553")  # Apery
    zeta_4 = pi ** 4 / 90  # 1.0823...
    zeta_half = Decimal("-1.4603545088095868128894991525152980124672293310125814905428860878")  # zeta(1/2)
    euler_gamma = Decimal("0.5772156649015328606065120900824024310421593359399235988057672349")
    catalan = Decimal("0.9159655941772190150546035149323841107741493742816721342664981196")
    ln2 = constants.ln2

    terms = [
        ("alpha * zeta(3)", alpha_d * zeta_3),
        ("alpha * zeta(2)", alpha_d * zeta_2),
        ("alpha * |zeta(1/2)|", alpha_d * abs(zeta_half)),
        ("alpha * zeta(4)", alpha_d * zeta_4),
        ("alpha * gamma_EM", alpha_d * euler_gamma),
        ("alpha * Catalan", alpha_d * catalan),
        ("alpha * ln(2)", alpha_d * ln2),
        ("alpha^2 * zeta(3)", alpha_d ** 2 * zeta_3),
        ("alpha^2 * zeta(2)", alpha_d ** 2 * zeta_2),
        ("alpha * zeta(3) / pi", alpha_d * zeta_3 / pi),
        ("alpha * zeta(3) * sqrt(phi)", alpha_d * zeta_3 * constants.phi.sqrt()),
        ("zeta(3) / mu", zeta_3 / Decimal(str(delta_info["mu"]))),
        ("alpha * pi / 3", alpha_d * pi / 3),
        ("alpha * pi / 4", alpha_d * pi / 4),
        ("alpha * (pi - 3)", alpha_d * (pi - 3)),
        ("alpha * (4 - pi)", alpha_d * (4 - pi)),
    ]

    candidates = []
    for name, value in terms:
        if value == 0:
            continue
        ratio = float(delta / value)
        closeness = abs(ratio - round(ratio)) if abs(ratio) < 1000 else float("inf")
        candidates.append({
            "name": name,
            "expression": name,
            "value": float(value),
            "ratio": ratio,
            "nearest_integer": round(ratio) if abs(ratio) < 1000 else None,
            "closeness_to_integer": closeness if closeness != float("inf") else float("inf"),
        })

    valid = [c for c in candidates if c["nearest_integer"] is not None and c["nearest_integer"] != 0]
    valid.sort(key=lambda c: c["closeness_to_integer"])
    best = valid[0] if valid else None

    return {
        "candidates": candidates,
        "best_match": best,
        "assessment": (
            f"Best zeta match: delta / ({best['expression']}) = {best['ratio']:.4f} "
            f"(nearest integer: {best['nearest_integer']}, "
            f"fractional part: {best['closeness_to_integer']:.4f})"
            if best else "No close integer ratios found among zeta terms."
        ),
    }


def scan_composite_expressions(constants=None):
    """Search for composite expressions that match delta more precisely.

    Tests products and ratios of alpha, pi, phi, zeta values, mu to find
    expressions where delta/expr is very close to a simple rational number.

    Parameters
    ----------
    constants : SimpleNamespace or None

    Returns
    -------
    dict
        candidates : list of dict with keys name, value, ratio, rational_match, residual_ppm.
        best_matches : list of top 10 closest.
        assessment : str.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    pi = constants.pi
    phi = constants.phi

    delta_info = compute_delta(constants)
    delta = Decimal(str(delta_info["delta"]))
    mu = Decimal(str(delta_info["mu"]))

    zeta_3 = Decimal("1.2020569031595942853997381615114499907649862923404988817922715553")
    euler_gamma = Decimal("0.5772156649015328606065120900824024310421593359399235988057672349")
    catalan = Decimal("0.9159655941772190150546035149323841107741493742816721342664981196")
    sqrt_phi = phi.sqrt()

    # Simple rational targets: 1, 1/2, 1/3, 2/3, 1/4, 3/4, 1/5, 2/5, 3/5, 4/5, 1/6, 5/6
    rationals = [
        (1, 1), (1, 2), (1, 3), (2, 3), (1, 4), (3, 4),
        (1, 5), (2, 5), (3, 5), (4, 5), (1, 6), (5, 6),
        (1, 7), (2, 7), (3, 7), (1, 8), (3, 8), (5, 8), (7, 8),
        (1, 9), (2, 9), (4, 9), (1, 10), (3, 10), (7, 10), (9, 10),
        (1, 12), (5, 12), (7, 12), (11, 12),
        (2, 1), (3, 1), (4, 1), (5, 1), (6, 1), (7, 1), (8, 1),
    ]

    # Base expressions to test
    bases = [
        ("alpha", alpha),
        ("alpha^2", alpha ** 2),
        ("alpha * sqrt(phi)", alpha * sqrt_phi),
        ("alpha * phi", alpha * phi),
        ("alpha * zeta(3)", alpha * zeta_3),
        ("alpha * gamma_EM", alpha * euler_gamma),
        ("alpha * Catalan", alpha * catalan),
        ("alpha * pi", alpha * pi),
        ("alpha / pi", alpha / pi),
        ("alpha^2 * pi", alpha ** 2 * pi),
        ("alpha^2 * pi^2", alpha ** 2 * pi ** 2),
        ("alpha^2 * mu", alpha ** 2 * mu),
        ("alpha * (pi - 3)", alpha * (pi - 3)),
        ("alpha * ln(2)", alpha * constants.ln2),
        ("alpha * sqrt(2)", alpha * constants.sqrt2),
        ("alpha * sqrt(3)", alpha * constants.sqrt3),
        ("alpha * sqrt(5)", alpha * constants.sqrt5),
        ("alpha^2 * sqrt(phi)", alpha ** 2 * sqrt_phi),
        ("alpha * zeta(3) / pi", alpha * zeta_3 / pi),
        ("alpha^2 * zeta(3)", alpha ** 2 * zeta_3),
    ]

    candidates = []

    for base_name, base_val in bases:
        if base_val == 0:
            continue
        raw_ratio = delta / base_val

        for p, q in rationals:
            target = Decimal(p) / Decimal(q)
            frac_residual = float(abs(raw_ratio - target) / target) if target != 0 else float("inf")
            residual_ppm = frac_residual * 1e6

            if residual_ppm < 50000:  # Only keep candidates within 5%
                candidates.append({
                    "name": f"({p}/{q}) * {base_name}" if q != 1 else f"{p} * {base_name}",
                    "base_expression": base_name,
                    "rational": f"{p}/{q}" if q != 1 else str(p),
                    "value": float(target * base_val),
                    "ratio": float(raw_ratio),
                    "target_ratio": float(target),
                    "residual_ppm": residual_ppm,
                })

    candidates.sort(key=lambda c: c["residual_ppm"])
    best_matches = candidates[:10]

    return {
        "candidates": candidates,
        "best_matches": best_matches,
        "n_within_100ppm": sum(1 for c in candidates if c["residual_ppm"] < 100),
        "n_within_1000ppm": sum(1 for c in candidates if c["residual_ppm"] < 1000),
        "assessment": (
            f"Best composite match: delta = {best_matches[0]['name']}, "
            f"residual {best_matches[0]['residual_ppm']:.1f} ppm"
            if best_matches else "No close composite matches found."
        ),
    }


def scan_k_closed_forms(constants=None):
    """Search for closed-form expressions for k_exact itself (not just delta).

    Since k_exact ~ 1.2622, tests expressions near this value using
    combinations of phi, pi, alpha, zeta values.

    Parameters
    ----------
    constants : SimpleNamespace or None

    Returns
    -------
    dict
        candidates : list of dict with keys name, value, residual_ppm.
        best_matches : list of top 10.
        assessment : str.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    pi = constants.pi
    phi = constants.phi

    delta_info = compute_delta(constants)
    k_exact = Decimal(str(delta_info["k_exact"]))
    sqrt_phi = phi.sqrt()

    zeta_3 = Decimal("1.2020569031595942853997381615114499907649862923404988817922715553")
    euler_gamma = Decimal("0.5772156649015328606065120900824024310421593359399235988057672349")

    # Candidate closed forms for k ~ 1.2622
    expressions = [
        ("sqrt(phi)", sqrt_phi),
        ("sqrt(phi) - alpha*sqrt(phi)", sqrt_phi - alpha * sqrt_phi),
        ("sqrt(phi) * (1 - alpha)", sqrt_phi * (1 - alpha)),
        ("sqrt(phi) * (1 - alpha/pi)", sqrt_phi * (1 - alpha / pi)),
        ("sqrt(phi) * (1 - alpha/(2*pi))", sqrt_phi * (1 - alpha / (2 * pi))),
        ("sqrt(phi) * (1 - 3*alpha^2)", sqrt_phi * (1 - 3 * alpha ** 2)),
        ("sqrt(phi) * (1 - alpha^2*pi)", sqrt_phi * (1 - alpha ** 2 * pi)),
        ("sqrt(phi) * (1 - alpha*zeta(3))", sqrt_phi * (1 - alpha * zeta_3)),
        ("sqrt(phi) * (1 - alpha*gamma_EM)", sqrt_phi * (1 - alpha * euler_gamma)),
        ("5/4", Decimal(5) / Decimal(4)),
        ("zeta(3)", zeta_3),
        ("pi/e^(1/2)", pi / constants.e.sqrt()),
        ("4/pi", Decimal(4) / pi),
        ("phi^(3/2) / sqrt(2)", phi ** (Decimal(3) / Decimal(2)) / constants.sqrt2),
        ("(1 + sqrt(5))/2 - 1/phi^2", phi - 1 / phi ** 2),
        ("sqrt(phi) - alpha", sqrt_phi - alpha),
        ("sqrt(phi) - alpha/pi", sqrt_phi - alpha / pi),
        ("sqrt(phi) * exp(-alpha)", sqrt_phi * (-alpha).exp()),
        ("sqrt(phi) * exp(-alpha*sqrt(phi))", sqrt_phi * (-alpha * sqrt_phi).exp()),
        ("sqrt(phi) / (1 + alpha)", sqrt_phi / (1 + alpha)),
        ("sqrt(phi) / (1 + alpha*sqrt(phi))", sqrt_phi / (1 + alpha * sqrt_phi)),
        ("(phi - 1/phi) / sqrt(2)", (phi - 1 / phi) / constants.sqrt2),
        ("sqrt(phi) * (1 - alpha*ln(2))", sqrt_phi * (1 - alpha * constants.ln2)),
    ]

    candidates = []
    for name, value in expressions:
        diff = k_exact - value
        residual_ppm = float(abs(diff) / k_exact * Decimal("1e6"))
        candidates.append({
            "name": name,
            "value": float(value),
            "k_exact": float(k_exact),
            "difference": float(diff),
            "residual_ppm": residual_ppm,
        })

    candidates.sort(key=lambda c: c["residual_ppm"])
    best_matches = candidates[:10]

    return {
        "candidates": candidates,
        "best_matches": best_matches,
        "assessment": (
            f"Best closed-form for k: {best_matches[0]['name']} = {best_matches[0]['value']:.6f}, "
            f"residual {best_matches[0]['residual_ppm']:.1f} ppm vs k_exact = {float(k_exact):.6f}"
            if best_matches else "No close closed-form matches found."
        ),
    }


def analyze_delta_structure(constants=None):
    """Analyze the mathematical structure of delta to guide interpretation.

    Computes delta in various unit systems and checks whether it has
    a natural scale in terms of known physics.

    Parameters
    ----------
    constants : SimpleNamespace or None

    Returns
    -------
    dict
        delta, delta_over_alpha, delta_over_alpha_sq, delta_times_mu,
        delta_over_schwinger, delta_as_alpha_power, interpretation, assessment.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    pi = constants.pi

    delta_info = compute_delta(constants)
    delta = Decimal(str(delta_info["delta"]))
    mu = Decimal(str(delta_info["mu"]))

    schwinger = alpha / (2 * pi)

    # What power of alpha is |delta|?
    # |delta| ~ alpha^n => n = ln(|delta|) / ln(alpha)
    import math as _math
    abs_delta = abs(float(delta))
    alpha_power = _math.log(abs_delta) / _math.log(float(alpha)) if abs_delta > 0 else float("nan")

    return {
        "delta": float(delta),
        "delta_over_alpha": float(delta / alpha),
        "delta_over_alpha_sq": float(delta / alpha ** 2),
        "delta_times_mu": float(delta * mu),
        "delta_over_schwinger": float(delta / schwinger),
        "delta_as_alpha_power": alpha_power,
        "delta_over_sqrt_phi": float(delta / constants.phi.sqrt()),
        "delta_over_phi": float(delta / constants.phi),
        "interpretation": (
            f"delta ~ alpha^{alpha_power:.2f}, suggesting it is an O(alpha^1) correction. "
            f"delta/alpha = {float(delta / alpha):.4f}, delta/(alpha*sqrt(phi)) = "
            f"{float(delta / (alpha * constants.phi.sqrt())):.4f}."
        ),
        "assessment": (
            f"delta = {float(delta):.6e} is O(alpha) in magnitude (alpha^{alpha_power:.2f}). "
            f"The ratio delta/alpha = {float(delta / alpha):.4f} is not a recognizable constant. "
            f"delta/(alpha*sqrt(phi)) = {float(delta / (alpha * constants.phi.sqrt())):.4f} "
            f"is tantalizingly close to 1 but not exact."
        ),
    }


def summarize_residual_mapping(constants=None):
    """Summarize the complete residual mapping analysis.

    Parameters
    ----------
    constants : SimpleNamespace or None

    Returns
    -------
    dict
        delta_info, g2_scan, zeta_scan, composite_scan, k_closed_forms,
        delta_structure, key_finding, honest_assessment.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    delta_info = compute_delta(constants)
    g2_scan = scan_anomalous_magnetic_moment(constants)
    zeta_scan = scan_zeta_functions(constants)
    composite_scan = scan_composite_expressions(constants)
    k_closed_forms = scan_k_closed_forms(constants)
    delta_structure = analyze_delta_structure(constants)

    # Collect all best matches
    all_bests = []
    if composite_scan["best_matches"]:
        all_bests.append(("composite", composite_scan["best_matches"][0]))
    if k_closed_forms["best_matches"]:
        all_bests.append(("k_closed_form", k_closed_forms["best_matches"][0]))

    key_finding = (
        f"The residual delta = sqrt(phi) - k_exact = {delta_info['delta']:.6e} is O(alpha) "
        f"in magnitude. No exact match to a simple SM constant combination was found, "
        f"but several expressions achieve sub-1000 ppm agreement."
    )

    honest_assessment = (
        "The gap between sqrt(phi) and the exact offset k is real and significant at ~0.01. "
        "It scales as O(alpha), suggesting it may be a radiative correction to sqrt(phi). "
        "However, no combination of known SM constants (g-2 coefficients, zeta values, "
        "Euler-Mascheroni, Catalan) reproduces delta exactly. The closest matches involve "
        "alpha * sqrt(phi) and alpha * zeta(3), both within ~10% of delta, but neither is exact. "
        "This could mean: (1) the correction involves higher-order terms not yet tested, "
        "(2) sqrt(phi) itself is not the correct theoretical value, or "
        "(3) the mu-structure formula requires a different functional form. "
        "The search is systematic but not exhaustive -- algebraic number theory methods "
        "(PSLQ, LLL lattice reduction) would be needed for a definitive answer."
    )

    return {
        "delta_info": delta_info,
        "g2_scan": g2_scan,
        "zeta_scan": zeta_scan,
        "composite_scan": composite_scan,
        "k_closed_forms": k_closed_forms,
        "delta_structure": delta_structure,
        "key_finding": key_finding,
        "honest_assessment": honest_assessment,
    }
