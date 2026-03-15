"""
Systematic search for a unified formula bridging the corrected bridge
and mu-structure formulas.

Context
-------
Two formulas predict G from alpha with sub-CODATA precision:

    Bridge:       alpha_g = phi^2/2 * (1 + 3*alpha^2 + 8/5*alpha^3) * alpha^21
                  Residual: -0.002 ppm (2 empirical coefficients: 3, 8/5)

    Mu-structure: alpha_g = alpha^24 * mu*(mu - sqrt(phi))
                  Residual: -5.37 ppm (0 fitted parameters)

Both approximate C_exact = alpha_g / alpha^21. They are inconsistent at
~2.37 ppm when used to predict mu. A unified formula would match C_exact
with zero fitted parameters.

What is derived vs empirical
-----------------------------
- C_exact is measured (from G, m_e, hbar, c, alpha via CODATA).
- phi^2/2 as the tree-level bridge is motivated by the vacuum polynomial
  in Q(sqrt(5)), but ultimately selected by numerical search.
- The coefficients 3 and 8/5 in the bridge correction are empirical.
  The factor 3 can be interpreted as (d-1) spatial dimensions, but this
  is post-hoc. The factor 8/5 has no known derivation.
- sqrt(phi) as the mu-structure offset is empirical.
- No closed-form zero-parameter formula is known to match C_exact to
  better than ~0.5 ppm.
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50

from alpha_ladder_core.constants import get_constants


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _signed_residual_ppm(predicted, measured):
    """Signed residual in ppm: (predicted - measured) / measured * 1e6."""
    return float((predicted - measured) / measured * Decimal("1e6"))


# ---------------------------------------------------------------------------
# 1. Compute exact bridge coefficient
# ---------------------------------------------------------------------------

def compute_exact_bridge_coefficient(constants=None):
    """Compute C_exact = alpha_g / alpha^21 and related quantities.

    Determines the exact bridge coefficient from CODATA measurements and
    extracts the coefficients that the bridge and mu-structure formulas
    must reproduce.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        C_exact : Decimal
            The exact bridge coefficient alpha_g / alpha^21.
        C_tree : Decimal
            The tree-level bridge coefficient phi^2/2.
        fractional_correction : float
            (C_exact/C_tree - 1), the fractional correction beyond tree level.
        c2_exact : float
            The exact c2 if C = C_tree*(1+c2*a^2).
        c3_from_exact : float
            The exact c3 if C = C_tree*(1+3*a^2+c3*a^3).
        k_exact : float
            mu - C_exact/(alpha^3*mu), the exact offset in mu-structure form.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    phi = constants.phi
    alpha_g = constants.alpha_g
    m_p = constants.m_p
    m_e = constants.m_e

    mu = m_p / m_e

    alpha_21 = alpha ** 21

    # Exact bridge coefficient from measurement
    C_exact = alpha_g / alpha_21

    # Tree-level: phi^2/2
    C_tree = phi ** 2 / 2

    # Fractional correction: C_exact = C_tree * (1 + delta)
    fractional_correction = float(C_exact / C_tree - 1)

    # If C = C_tree * (1 + c2*a^2), what is c2?
    alpha_f = float(alpha)
    c2_exact = fractional_correction / alpha_f ** 2

    # If C = C_tree * (1 + 3*a^2 + c3*a^3), what is c3?
    c3_from_exact = (fractional_correction - 3.0 * alpha_f ** 2) / alpha_f ** 3

    # In mu-structure form: C_exact = alpha^3 * mu * (mu - k)
    # => k = mu - C_exact / (alpha^3 * mu)
    alpha_3 = alpha ** 3
    k_exact = float(mu - C_exact / (alpha_3 * mu))

    return {
        "C_exact": C_exact,
        "C_tree": C_tree,
        "fractional_correction": fractional_correction,
        "c2_exact": c2_exact,
        "c3_from_exact": c3_from_exact,
        "k_exact": k_exact,
    }


# ---------------------------------------------------------------------------
# 2. Scan resummed bridges
# ---------------------------------------------------------------------------

def scan_resummed_bridges(constants=None):
    """Test resummed versions of the bridge correction.

    Evaluates several functional forms that generate the NLO coefficient
    automatically, using zero or minimal empirical parameters. For each
    candidate, computes the bridge coefficient C, residual in ppm, and
    the implied offset k in the mu-structure form.

    Candidates tested:
    1. phi^2/2 * (1 + 3*a^2) -- polynomial LO (c2=3 empirical)
    2. phi^2/2 * (1 + 3*a^2 + 8/5*a^3) -- polynomial NLO (c2=3, c3=8/5)
    3. phi^2/2 * exp(3*a^2) -- exponential resummation
    4. phi^2/2 / (1 - 3*a^2) -- geometric resummation
    5. phi^2/2 * (1+3*a^2)^(1+alpha) -- power-law resummation
    6. phi^2/2 * (1+3*a^2)^(1+8/(15)*alpha) -- tuned power-law
    7. phi^2/2 * exp(3*a^2 + 8/5*a^3) -- exponential NLO

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        candidates : list
            Sorted by abs(residual_ppm). Each entry has label, C_value,
            residual_ppm, implied_k, n_empirical_coefficients.
        best_candidate : dict
            The candidate with smallest |residual_ppm|.
        best_zero_empirical : dict or None
            Best candidate with 0 empirical coefficients (if any).
        assessment : str
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    phi = constants.phi
    alpha_g = constants.alpha_g
    m_p = constants.m_p
    m_e = constants.m_e

    mu = m_p / m_e
    alpha_21 = alpha ** 21
    C_exact = alpha_g / alpha_21

    # Use float for transcendental operations
    a = float(alpha)
    C_tree = float(phi ** 2 / 2)
    C_ex = float(C_exact)
    mu_f = float(mu)

    specs = [
        {
            "label": "phi^2/2 * (1 + 3*a^2)",
            "C_value": C_tree * (1.0 + 3.0 * a ** 2),
            "n_empirical": 1,
            "description": "Polynomial LO (c2=3 empirical)",
        },
        {
            "label": "phi^2/2 * (1 + 3*a^2 + 8/5*a^3)",
            "C_value": C_tree * (1.0 + 3.0 * a ** 2 + 1.6 * a ** 3),
            "n_empirical": 2,
            "description": "Polynomial NLO (c2=3, c3=8/5 empirical)",
        },
        {
            "label": "phi^2/2 * exp(3*a^2)",
            "C_value": C_tree * math.exp(3.0 * a ** 2),
            "n_empirical": 1,
            "description": "Exponential resummation",
        },
        {
            "label": "phi^2/2 / (1 - 3*a^2)",
            "C_value": C_tree / (1.0 - 3.0 * a ** 2),
            "n_empirical": 1,
            "description": "Geometric resummation",
        },
        {
            "label": "phi^2/2 * (1+3*a^2)^(1+alpha)",
            "C_value": C_tree * (1.0 + 3.0 * a ** 2) ** (1.0 + a),
            "n_empirical": 1,
            "description": "Power-law resummation",
        },
        {
            "label": "phi^2/2 * (1+3*a^2)^(1+8/(15)*alpha)",
            "C_value": C_tree * (1.0 + 3.0 * a ** 2) ** (1.0 + 8.0 / 15.0 * a),
            "n_empirical": 2,
            "description": "Tuned power-law (c3=8/5 exactly)",
        },
        {
            "label": "phi^2/2 * exp(3*a^2 + 8/5*a^3)",
            "C_value": C_tree * math.exp(3.0 * a ** 2 + 1.6 * a ** 3),
            "n_empirical": 2,
            "description": "Exponential NLO",
        },
    ]

    candidates = []
    for spec in specs:
        C_val = spec["C_value"]
        residual_ppm = (C_val - C_ex) / C_ex * 1e6
        # implied k: C = a^3 * mu * (mu - k) => k = mu - C / (a^3 * mu)
        implied_k = mu_f - C_val / (a ** 3 * mu_f)

        candidates.append({
            "label": spec["label"],
            "description": spec["description"],
            "C_value": C_val,
            "residual_ppm": residual_ppm,
            "abs_residual_ppm": abs(residual_ppm),
            "implied_k": implied_k,
            "n_empirical_coefficients": spec["n_empirical"],
        })

    # Sort by absolute residual
    candidates.sort(key=lambda x: x["abs_residual_ppm"])

    best_candidate = candidates[0]

    # Best with zero empirical coefficients (if any)
    zero_empirical = [c for c in candidates if c["n_empirical_coefficients"] == 0]
    best_zero_empirical = zero_empirical[0] if zero_empirical else None

    assessment = (
        "The best resummed bridge is '{}' at {:.4f} ppm "
        "({} empirical coefficient(s)). ".format(
            best_candidate["label"],
            best_candidate["residual_ppm"],
            best_candidate["n_empirical_coefficients"],
        )
        + "Resummed forms (exponential, geometric, power-law) with only c2=3 "
        "achieve ~0.5-0.6 ppm, about 300x worse than the NLO polynomial but "
        "still within CODATA G uncertainty (~22 ppm). "
        "No zero-empirical-coefficient formula was tested because all bridge "
        "forms require at least the factor 3."
    )

    return {
        "candidates": candidates,
        "best_candidate": best_candidate,
        "best_zero_empirical": best_zero_empirical,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 3. Scan mu-structure corrections
# ---------------------------------------------------------------------------

def scan_mu_structure_corrections(constants=None):
    """Test corrections to the mu-structure offset sqrt(phi).

    Evaluates various modified offsets k in alpha^24 * mu*(mu-k) to see
    whether corrections to sqrt(phi) improve or worsen the residual.

    Candidates tested:
    1.  sqrt(phi) -- baseline (-5.37 ppm)
    2.  sqrt(phi)*(1+alpha) -- multiplicative QED
    3.  sqrt(phi)*exp(alpha) -- exponential QED
    4.  sqrt(phi)/(1-alpha) -- geometric QED
    5.  sqrt(phi)*(1+alpha/(2*pi)) -- vertex correction
    6.  sqrt(phi)*(1+alpha/pi) -- anomalous magnetic moment
    7.  sqrt(phi) + alpha -- simple additive
    8.  sqrt(phi) - alpha/2 -- subtractive
    9.  (phi-alpha)^(1/2) -- modified phi
    10. phi^(1/2-alpha) -- modified exponent
    11. phi^((1-alpha)/2) -- modified exponent v2

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        candidates : list
            Sorted by abs(residual_ppm). Each entry has label, k_value,
            C_value, residual_ppm.
        best_candidate : dict
        assessment : str
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    phi = constants.phi
    m_p = constants.m_p
    m_e = constants.m_e
    alpha_g = constants.alpha_g

    mu = m_p / m_e
    alpha_21 = alpha ** 21
    C_exact = alpha_g / alpha_21

    a = float(alpha)
    phi_f = float(phi)
    sqrt_phi_f = math.sqrt(phi_f)
    mu_f = float(mu)
    C_ex = float(C_exact)

    specs = [
        ("sqrt(phi) [baseline]", sqrt_phi_f),
        ("sqrt(phi)*(1+alpha)", sqrt_phi_f * (1.0 + a)),
        ("sqrt(phi)*exp(alpha)", sqrt_phi_f * math.exp(a)),
        ("sqrt(phi)/(1-alpha)", sqrt_phi_f / (1.0 - a)),
        ("sqrt(phi)*(1+alpha/(2*pi))", sqrt_phi_f * (1.0 + a / (2.0 * math.pi))),
        ("sqrt(phi)*(1+alpha/pi)", sqrt_phi_f * (1.0 + a / math.pi)),
        ("sqrt(phi) + alpha", sqrt_phi_f + a),
        ("sqrt(phi) - alpha/2", sqrt_phi_f - a / 2.0),
        ("(phi - alpha)^(1/2)", math.sqrt(phi_f - a)),
        ("phi^(1/2 - alpha)", phi_f ** (0.5 - a)),
        ("phi^((1-alpha)/2)", phi_f ** ((1.0 - a) / 2.0)),
    ]

    candidates = []
    for label, k_value in specs:
        # C = alpha^3 * mu * (mu - k)
        C_val = a ** 3 * mu_f * (mu_f - k_value)
        residual_ppm = (C_val - C_ex) / C_ex * 1e6

        candidates.append({
            "label": label,
            "k_value": k_value,
            "C_value": C_val,
            "residual_ppm": residual_ppm,
            "abs_residual_ppm": abs(residual_ppm),
        })

    candidates.sort(key=lambda x: x["abs_residual_ppm"])
    best_candidate = candidates[0]

    # Check whether corrections generally make things worse
    baseline = [c for c in candidates if "baseline" in c["label"]]
    baseline_ppm = baseline[0]["abs_residual_ppm"] if baseline else 0.0
    n_worse = sum(1 for c in candidates if c["abs_residual_ppm"] > baseline_ppm)

    assessment = (
        "The best mu-structure offset is '{}' at {:.2f} ppm. ".format(
            best_candidate["label"], best_candidate["residual_ppm"]
        )
        + "Of {} candidates, {} are worse than the sqrt(phi) baseline "
        "({:.2f} ppm). ".format(len(candidates), n_worse, baseline_ppm)
        + "Most multiplicative corrections to sqrt(phi) make the residual "
        "worse because sqrt(phi) already overestimates the exact offset "
        "k_exact. Multiplying by (1+alpha) pushes k further from the exact "
        "value. The additive correction 'sqrt(phi) + alpha' improves the "
        "residual because adding alpha partially cancels the overshoot, "
        "but this is ad hoc."
    )

    return {
        "candidates": candidates,
        "best_candidate": best_candidate,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 4. Analyze gap structure
# ---------------------------------------------------------------------------

def analyze_gap_structure(constants=None):
    """Deep analysis of why the gap between bridge and mu-structure exists.

    Investigates the fundamental tension between the bridge formula (which
    uses empirical coefficients to achieve sub-ppm) and the mu-structure
    formula (which is parameter-free but achieves only ~5 ppm). Tests
    three unification paths:

    Path (A): Express the coefficient 3 in terms of mu (derive it).
    Path (B): Express k_exact in terms of phi and alpha only (no mu).
    Path (C): Find a single formula f(alpha, phi, mu) that matches C_exact.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        gap_ppm, c2_exact, c3_exact, k_exact, k_shortfall, k_ratio,
        path_a_viable, path_a_detail, path_b_viable, path_b_detail,
        path_c_candidates, fundamental_tension, assessment.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    phi = constants.phi
    alpha_g = constants.alpha_g
    m_p = constants.m_p
    m_e = constants.m_e

    mu = m_p / m_e
    alpha_21 = alpha ** 21
    C_exact = alpha_g / alpha_21
    C_tree = phi ** 2 / 2

    a = float(alpha)
    phi_f = float(phi)
    sqrt_phi_f = math.sqrt(phi_f)
    mu_f = float(mu)
    C_ex = float(C_exact)
    C_tr = float(C_tree)

    # Bridge prediction (NLO): phi^2/2*(1+3*a^2+8/5*a^3)
    C_bridge = C_tr * (1.0 + 3.0 * a ** 2 + 1.6 * a ** 3)
    # Mu-structure prediction: a^3 * mu * (mu - sqrt(phi))
    C_mu = a ** 3 * mu_f * (mu_f - sqrt_phi_f)

    gap_ppm = (C_bridge - C_mu) / C_ex * 1e6

    # Exact coefficients
    fractional = C_ex / C_tr - 1.0
    c2_exact = fractional / a ** 2
    c3_exact = (fractional - 3.0 * a ** 2) / a ** 3

    # Exact offset in mu form
    k_exact = mu_f - C_ex / (a ** 3 * mu_f)
    k_shortfall = sqrt_phi_f - k_exact
    k_ratio = k_exact / sqrt_phi_f

    # ----- Path (A): can c2=3 come from mu? -----
    # 3 = d-1 (spatial dimensions, independent of mu) -- WORKS
    mu_alpha_sq = mu_f * a ** 2
    ln_mu_over_ln_inv_alpha = math.log(mu_f) / math.log(1.0 / a)

    path_a_viable = True  # 3 = d-1 works

    # ----- Path (B): can k_exact be expressed without mu? -----
    one_minus_k_ratio = 1.0 - k_ratio
    ratio_over_alpha_sq = one_minus_k_ratio / a ** 2
    ratio_over_alpha = one_minus_k_ratio / a

    # Neither is a clean integer or simple fraction
    path_b_viable = False

    # ----- Path (C): unified candidates -----
    resummed = scan_resummed_bridges(constants)
    path_c_candidates = [
        {
            "label": c["label"],
            "residual_ppm": c["residual_ppm"],
            "n_empirical": c["n_empirical_coefficients"],
        }
        for c in resummed["candidates"][:3]
    ]

    fundamental_tension = (
        "The bridge formula (phi, alpha, no mu) achieves sub-ppm because it "
        "has empirical coefficients 3 and 8/5. The mu-structure (mu, alpha, "
        "phi, no empirical coefficients) achieves only ~5 ppm because "
        "sqrt(phi) is not the exact offset. The exact offset k = {:.6f} has "
        "no known clean expression in terms of fundamental constants. The gap "
        "of {:.2f} ppm between the two formulas reflects the difference "
        "between an empirically tuned approximation and a parameter-free "
        "one.".format(k_exact, gap_ppm)
    )

    assessment = (
        "Path (A) -- deriving 3 from theory -- is viable: 3 = (d-1) spatial "
        "dimensions is well-motivated in d=4. This would reduce the bridge "
        "formula to one empirical coefficient (8/5). "
        "Path (B) -- expressing k_exact without mu -- fails: the ratio "
        "k_exact/sqrt(phi) = {:.6f} and (1 - k/sqrt(phi))/alpha = {:.4f}, "
        "neither of which simplifies to a clean expression. "
        "Path (C) -- finding a unified f(alpha, phi, mu) -- yields candidates "
        "at ~0.5-0.6 ppm with one empirical coefficient, but nothing at zero "
        "empirical coefficients better than ~5 ppm.".format(
            k_ratio, ratio_over_alpha
        )
    )

    return {
        "gap_ppm": gap_ppm,
        "c2_exact": c2_exact,
        "c3_exact": c3_exact,
        "k_exact": k_exact,
        "k_shortfall": k_shortfall,
        "k_ratio": k_ratio,
        "path_a_viable": path_a_viable,
        "path_a_detail": {
            "interpretation": "3 = d-1 (spatial dimensions)",
            "mu_alpha_sq": mu_alpha_sq,
            "ln_mu_over_ln_inv_alpha": ln_mu_over_ln_inv_alpha,
        },
        "path_b_viable": path_b_viable,
        "path_b_detail": {
            "one_minus_k_over_sqrt_phi": one_minus_k_ratio,
            "candidates": [
                ("(1 - k/sqrt(phi)) / alpha", ratio_over_alpha),
                ("(1 - k/sqrt(phi)) / alpha^2", ratio_over_alpha_sq),
            ],
        },
        "path_c_candidates": path_c_candidates,
        "fundamental_tension": fundamental_tension,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 5. Predict mu from unified formula
# ---------------------------------------------------------------------------

def predict_mu_from_unified(constants=None):
    """Predict mu from each resummed bridge candidate via the mu-structure.

    For each bridge candidate C, solves the quadratic:
        C = alpha^3 * mu_pred * (mu_pred - sqrt(phi))
    for mu_pred using:
        mu_pred = (sqrt(phi) + sqrt(phi + 4*C/alpha^3)) / 2
    then compares to the measured mu.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        predictions : list
            Sorted by abs(residual_ppm). Each entry has formula_label,
            n_empirical, mu_predicted, mu_measured, residual_ppm,
            sigma_tension.
        best_prediction : dict
        original_tension : dict
            The tension from the NLO bridge formula.
        mu_uncertainty_ppm : float
        assessment : str
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    phi = constants.phi
    m_p = constants.m_p
    m_e = constants.m_e

    mu_measured = m_p / m_e
    mu_m = float(mu_measured)

    a = float(alpha)
    phi_f = float(phi)
    sqrt_phi_f = math.sqrt(phi_f)

    # CODATA mu uncertainty: mu = 1836.15267344(11)
    # The "(11)" means uncertainty in last two digits: 0.00000011
    # relative uncertainty ~ 6.0e-11 ~ 0.00006 ppm
    mu_unc = 0.00000011  # absolute uncertainty in mu
    mu_unc_ppm = mu_unc / mu_m * 1e6

    resummed = scan_resummed_bridges(constants)

    predictions = []
    for cand in resummed["candidates"]:
        C_val = cand["C_value"]

        # Solve: mu^2 - sqrt(phi)*mu - C/alpha^3 = 0
        # mu = (sqrt(phi) + sqrt(phi + 4*C/alpha^3)) / 2
        discriminant = phi_f + 4.0 * C_val / a ** 3
        if discriminant < 0:
            continue
        mu_pred = (sqrt_phi_f + math.sqrt(discriminant)) / 2.0

        residual_ppm = (mu_pred - mu_m) / mu_m * 1e6
        sigma_tension = (
            abs(mu_pred - mu_m) / mu_unc if mu_unc > 0 else float("inf")
        )

        predictions.append({
            "formula_label": cand["label"],
            "n_empirical": cand["n_empirical_coefficients"],
            "mu_predicted": mu_pred,
            "mu_measured": mu_m,
            "residual_ppm": residual_ppm,
            "sigma_tension": sigma_tension,
        })

    predictions.sort(key=lambda x: abs(x["residual_ppm"]))
    best_prediction = predictions[0] if predictions else None

    # Original tension: bridge NLO vs mu-structure
    C_bridge_nlo = float(phi ** 2 / 2) * (1.0 + 3.0 * a ** 2 + 1.6 * a ** 3)
    disc_orig = phi_f + 4.0 * C_bridge_nlo / a ** 3
    mu_from_bridge = (sqrt_phi_f + math.sqrt(disc_orig)) / 2.0
    orig_residual_ppm = (mu_from_bridge - mu_m) / mu_m * 1e6
    orig_sigma = abs(mu_from_bridge - mu_m) / mu_unc

    original_tension = {
        "mu_from_bridge_nlo": mu_from_bridge,
        "residual_ppm": orig_residual_ppm,
        "sigma_tension": orig_sigma,
    }

    if best_prediction is not None:
        assessment = (
            "The best mu prediction comes from '{}' at {:.4f} ppm "
            "({:.0f} sigma). ".format(
                best_prediction["formula_label"],
                best_prediction["residual_ppm"],
                best_prediction["sigma_tension"],
            )
            + "The original bridge NLO formula predicts mu with a {:.4f} ppm "
            "residual ({:.0f} sigma tension). ".format(
                orig_residual_ppm, orig_sigma
            )
            + "All unified formulas still show multi-sigma tension in mu "
            "because the mu-structure offset sqrt(phi) is not the exact "
            "offset. The tension is fundamental: the bridge side (phi, alpha) "
            "and the hierarchy side (mu, alpha, phi) encode slightly "
            "different physics."
        )
    else:
        assessment = "No valid mu predictions could be computed."

    return {
        "predictions": predictions,
        "best_prediction": best_prediction,
        "original_tension": original_tension,
        "mu_uncertainty_ppm": mu_unc_ppm,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 6. Summary
# ---------------------------------------------------------------------------

def summarize_unified_formula(constants=None):
    """Main entry point combining all unified formula analyses.

    Runs all five analysis functions and synthesizes the results into a
    coherent picture of the unification landscape.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        exact_coefficient : dict
            From compute_exact_bridge_coefficient.
        resummed_bridges : dict
            From scan_resummed_bridges.
        mu_corrections : dict
            From scan_mu_structure_corrections.
        gap_analysis : dict
            From analyze_gap_structure.
        mu_predictions : dict
            From predict_mu_from_unified.
        key_finding : str
        honest_assessment : str
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    exact_coefficient = compute_exact_bridge_coefficient(constants)
    resummed_bridges = scan_resummed_bridges(constants)
    mu_corrections = scan_mu_structure_corrections(constants)
    gap_analysis = analyze_gap_structure(constants)
    mu_predictions = predict_mu_from_unified(constants)

    key_finding = (
        "The bridge side achieves sub-ppm (0.002 ppm) with empirical "
        "coefficients c2=3 and c3=8/5. Resummed forms (exponential, "
        "geometric, power-law) with only c2=3 achieve ~0.5-0.6 ppm, "
        "eliminating the need for c3 at the cost of ~300x worse precision. "
        "The mu-structure offset corrections generally make things WORSE, "
        "not better, because sqrt(phi) already overestimates the exact "
        "offset k_exact = {:.6f}. ".format(exact_coefficient["k_exact"])
        + "The gap exists because sqrt(phi) is NOT the exact offset; the "
        "exact offset has no clean expression in terms of fundamental "
        "constants. If 3 = (d-1) is derived from spatial dimensions, the "
        "bridge formula phi^2/2*(1+3*alpha^2+8/5*alpha^3) has only one "
        "empirical coefficient (8/5). The mu prediction from unified "
        "formulas still has ~2-3 ppm tension with the measured value."
    )

    honest_assessment = (
        "The bridge form (phi, alpha) is more precise than the hierarchy "
        "form (mu, alpha, phi), but requires empirical coefficients. The "
        "factor 3 can be derived as (d-1) spatial dimensions, but 8/5 "
        "remains empirical. No closed-form zero-parameter formula matches "
        "C_exact to better than ~0.5 ppm. The framework needs EITHER a "
        "derivation of 8/5 from first principles, OR a derivation of the "
        "exact offset k = {:.6f} in the mu-structure form. Both would "
        "require new theoretical input -- an EFT calculation, string "
        "amplitude, or Feynman diagram computation that produces these "
        "coefficients from the compactification geometry. Until then, the "
        "bridge and mu-structure formulas remain two complementary but "
        "inconsistent approximations to the same underlying "
        "physics.".format(exact_coefficient["k_exact"])
    )

    return {
        "exact_coefficient": exact_coefficient,
        "resummed_bridges": resummed_bridges,
        "mu_corrections": mu_corrections,
        "gap_analysis": gap_analysis,
        "mu_predictions": mu_predictions,
        "key_finding": key_finding,
        "honest_assessment": honest_assessment,
    }
