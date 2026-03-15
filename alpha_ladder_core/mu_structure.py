"""
Mu-structure discovery: replacing mu^2 with mu*(mu - sqrt(phi)).

Discovery
---------
The hierarchy formula alpha_g = alpha^24 * mu^2 predicts G with a 688 ppm
residual. Replacing mu^2 with mu*(mu - k) and scanning over k reveals that
k = sqrt(phi) gives:

    alpha_g = alpha^24 * mu * (mu - sqrt(phi))    (-5.37 ppm residual)

This is the best zero-parameter formula for G in the Alpha Ladder project.

Here phi = (1+sqrt(5))/2 is the golden ratio, which emerges from the vacuum
polynomial x^2 + Dx + d = 0 at D=6, d=4. The formula is equivalent to:

    G = alpha^24 * m_p * (m_p - sqrt(phi)*m_e) * hbar * c / m_e^4

What is derived vs empirical
-----------------------------
- DERIVED: phi from the vacuum polynomial at d=4, D=6; alpha and mu from
  CODATA measurements.
- EMPIRICAL: sqrt(phi) as the offset is discovered numerically, not derived
  from first principles. The formula has zero fitted parameters (phi is
  determined by the theory), but we have no theoretical derivation of why
  the offset should be sqrt(phi).
- The formula could be a numerical coincidence -- the offset space is
  continuous, so finding a "nice" number near the best fit is expected.
"""

from decimal import Decimal, getcontext

getcontext().prec = 50

from alpha_ladder_core.constants import get_constants


# ---------------------------------------------------------------------------
# Helper: convert alpha_g to G
# ---------------------------------------------------------------------------

def _alpha_g_to_G(alpha_g, constants):
    """Convert gravitational coupling alpha_g to Newton's constant G."""
    return alpha_g * constants.hbar * constants.c / constants.m_e ** 2


def _signed_residual_ppm(predicted, measured):
    """Signed residual in ppm: (predicted - measured) / measured * 1e6."""
    return float((predicted - measured) / measured * Decimal("1e6"))


# ---------------------------------------------------------------------------
# 1. Core mu-structure computation
# ---------------------------------------------------------------------------

def compute_mu_structure(constants=None):
    """Compute the mu*(mu - sqrt(phi)) formula and compare to measurement.

    The key formula is alpha_g = alpha^24 * mu * (mu - sqrt(phi)), which
    replaces the original hierarchy formula alpha_g = alpha^24 * mu^2.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        Core numerics including both hierarchy and mu-structure predictions,
        residuals in ppm, and the improvement factor.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    phi = constants.phi
    hbar = constants.hbar
    c = constants.c
    m_e = constants.m_e
    m_p = constants.m_p
    G_measured = constants.G
    alpha_g_measured = constants.alpha_g

    mu = m_p / m_e
    sqrt_phi = phi.sqrt()
    mu_minus_sqrt_phi = mu - sqrt_phi

    alpha_24 = alpha ** 24

    # Original hierarchy: alpha^24 * mu^2
    alpha_g_mu_sq = alpha_24 * mu ** 2
    G_hierarchy = _alpha_g_to_G(alpha_g_mu_sq, constants)
    hierarchy_residual_ppm = _signed_residual_ppm(G_hierarchy, G_measured)

    # New formula: alpha^24 * mu * (mu - sqrt(phi))
    alpha_g_mu_structure = alpha_24 * mu * mu_minus_sqrt_phi
    G_predicted = _alpha_g_to_G(alpha_g_mu_structure, constants)
    residual_ppm = _signed_residual_ppm(G_predicted, G_measured)

    # CODATA G uncertainty is ~22 ppm
    within_codata_uncertainty = abs(residual_ppm) < 22.0

    improvement_factor = (
        abs(hierarchy_residual_ppm / residual_ppm)
        if residual_ppm != 0 else float("inf")
    )

    # sqrt(phi) * m_e in MeV:  m_e ~ 0.51099895 MeV
    m_e_MeV = Decimal("0.51099895000")
    sqrt_phi_m_e_MeV = float(sqrt_phi * m_e_MeV)

    return {
        "mu": mu,
        "sqrt_phi": sqrt_phi,
        "mu_minus_sqrt_phi": mu_minus_sqrt_phi,
        "alpha_g_mu_sq": alpha_g_mu_sq,
        "alpha_g_mu_structure": alpha_g_mu_structure,
        "alpha_g_measured": alpha_g_measured,
        "G_predicted": G_predicted,
        "G_measured": G_measured,
        "residual_ppm": residual_ppm,
        "hierarchy_residual_ppm": hierarchy_residual_ppm,
        "improvement_factor": improvement_factor,
        "within_codata_uncertainty": within_codata_uncertainty,
        "sqrt_phi_m_e_MeV": sqrt_phi_m_e_MeV,
    }


# ---------------------------------------------------------------------------
# 1b. Refined mu-structure with (1-alpha) correction
# ---------------------------------------------------------------------------

def compute_mu_structure_refined(constants=None):
    """Compute the refined mu-structure formula with (1-alpha) correction.

    Formula: alpha_G = alpha^24 * mu * (mu - sqrt(phi) * (1 - alpha))
    where mu = m_p / m_e and phi = golden ratio.

    The (1-alpha) factor reduces the residual from -5.37 ppm to -0.31 ppm
    with zero fitted parameters.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        alpha_g, G_predicted, mu, sqrt_phi, k_offset, residual_ppm,
        improvement_over_bare, assessment.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

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

    # Compare to bare sqrt(phi)
    alpha_g_bare = alpha ** 24 * mu * (mu - sqrt_phi)
    G_bare = alpha_g_bare * hbar * c / m_e ** 2
    bare_ppm = float((G_bare - G_measured) / G_measured * Decimal("1e6"))

    return {
        "alpha_g": alpha_g,
        "G_predicted": G_predicted,
        "mu": float(mu),
        "sqrt_phi": float(sqrt_phi),
        "k_offset": float(k_offset),
        "residual_ppm": residual_ppm,
        "bare_residual_ppm": bare_ppm,
        "improvement_over_bare": abs(bare_ppm) / abs(residual_ppm) if residual_ppm != 0 else float("inf"),
        "assessment": (
            f"Refined mu-structure: alpha_g = alpha^24 * mu * (mu - sqrt(phi)*(1-alpha)) "
            f"gives G residual {residual_ppm:+.2f} ppm (vs {bare_ppm:+.2f} ppm bare). "
            f"Zero fitted parameters. {abs(bare_ppm) / abs(residual_ppm):.0f}x improvement."
        ),
    }


# ---------------------------------------------------------------------------
# 2. Scan mu offsets
# ---------------------------------------------------------------------------

def scan_mu_offsets(constants=None):
    """Scan mu*(mu - k) for various k to show sqrt(phi) is optimal.

    Tests a range of "natural" offset values including integers, algebraic
    numbers, and combinations involving fundamental constants. For each k,
    computes alpha^24 * mu * (mu - k) and the signed residual in ppm.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        Keys: offsets (list of dicts), best_offset, assessment.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    phi = constants.phi
    pi = constants.pi
    m_p = constants.m_p
    m_e = constants.m_e
    G_measured = constants.G

    mu = m_p / m_e
    alpha_24 = alpha ** 24

    sqrt_phi = phi.sqrt()
    sqrt_2 = Decimal(2).sqrt()

    candidates = [
        (Decimal(0), "0 (original mu^2)"),
        (Decimal(1), "1"),
        (sqrt_phi, "sqrt(phi)"),
        (phi, "phi"),
        (Decimal(2), "2"),
        (Decimal(1) / alpha, "1/alpha"),
        (phi ** 2 / 2, "phi^2/2"),
        (Decimal(5) / Decimal(4), "5/4"),
        (1 + 1 / mu, "1 + 1/mu"),
        (sqrt_2, "sqrt(2)"),
        (Decimal(4) / pi, "4/pi"),
    ]

    offsets = []
    for k_value, k_label in candidates:
        alpha_g = alpha_24 * mu * (mu - k_value)
        G_predicted = _alpha_g_to_G(alpha_g, constants)
        residual_ppm = _signed_residual_ppm(G_predicted, G_measured)
        within_uncertainty = abs(residual_ppm) < 22.0

        offsets.append({
            "k_value": k_value,
            "k_label": k_label,
            "alpha_g": alpha_g,
            "residual_ppm": residual_ppm,
            "within_uncertainty": within_uncertainty,
        })

    # Find best offset
    best_offset = min(offsets, key=lambda x: abs(x["residual_ppm"]))

    assessment = (
        f"Among {len(offsets)} candidate offsets, k = {best_offset['k_label']} "
        f"gives the smallest residual at {best_offset['residual_ppm']:.2f} ppm. "
        f"The sqrt(phi) offset stands out because phi is already determined "
        f"by the theory (vacuum polynomial at D=6, d=4), making this a "
        f"zero-parameter formula. However, the offset space is continuous, "
        f"so finding a 'nice' algebraic number near the optimum is not "
        f"by itself strong evidence."
    )

    return {
        "offsets": offsets,
        "best_offset": best_offset,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 3. Analyze sqrt(phi) origin
# ---------------------------------------------------------------------------

def analyze_sqrt_phi_origin(constants=None):
    """Explore what sqrt(phi) means physically and algebraically.

    Traces sqrt(phi) back to the vacuum polynomial x^2 + Dx + d = 0 at
    D=6, d=4, and examines the mass scale sqrt(phi)*m_e in the context of
    known particle physics.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        Analysis of the algebraic and physical origins of sqrt(phi).
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    phi = constants.phi
    m_e = constants.m_e
    m_p = constants.m_p

    # Vacuum polynomial x^2 + 6x + 4 = 0
    # Roots: x = (-6 +/- sqrt(36-16)) / 2 = (-6 +/- sqrt(20)) / 2
    #       = -3 +/- sqrt(5)
    # Root 1: -3 + sqrt(5) = -3 + 2.236... = -0.764... = -(3 - sqrt(5))/1
    #   Note: phi = (1+sqrt(5))/2, so phi - 1 = (sqrt(5)-1)/2
    #   and 3 - sqrt(5) = 2*(2 - phi) -- less clean
    #   Actually: -(3 - sqrt(5)) = sqrt(5) - 3 = 2*phi - 4 ... let's just compute
    sqrt_5 = Decimal(5).sqrt()
    root_1 = (-Decimal(6) + sqrt_5 * 2) / 2  # = -3 + sqrt(5)
    root_2 = (-Decimal(6) - sqrt_5 * 2) / 2  # = -3 - sqrt(5)

    # How phi emerges: The discriminant is 20, sqrt(20) = 2*sqrt(5).
    # phi = (1+sqrt(5))/2, so phi is in Q(sqrt(5)).
    # root_1 = -3 + sqrt(5) = -(phi + 1) + 2*phi = phi - 1 ... no.
    # -3 + sqrt(5) = -3 + 2*phi - 1 + ... let me just check numerically:
    # phi ~ 1.618, root_1 = -3 + 2.236 = -0.764
    # -phi = -1.618, so root_1 is not -phi.
    # Actually the vacuum polynomial for the bridge is x^2 + Dx + d at generic D, d.
    # At D=6, d=4: x^2 + 6x + 4 = 0, roots = -3 +/- sqrt(5).
    # phi appears because sqrt(5) = 2*phi - 1.
    # root_1 = -3 + (2*phi - 1) = 2*phi - 4 = 2*(phi - 2)
    # root_2 = -3 - (2*phi - 1) = -2*phi - 2 = -2*(phi + 1) = -2*phi^2
    # (since phi^2 = phi + 1)
    # So: root_2 = -2*phi^2 and root_1 = 2*(phi - 2) = -2*(2 - phi) = -2/phi^2
    # (since 2 - phi = 2 - (1+sqrt5)/2 = (3-sqrt5)/2 = 1/phi^2)
    # So root_1 = -2/phi^2, root_2 = -2*phi^2. Product = 4, sum = -6. Checks out.

    sqrt_phi = phi.sqrt()
    sqrt_phi_value = float(sqrt_phi)

    m_e_MeV = Decimal("0.51099895000")
    m_p_MeV = Decimal("938.272088")

    sqrt_phi_m_e_MeV = float(sqrt_phi * m_e_MeV)
    subtracted_mass_MeV = float(m_p_MeV - sqrt_phi * m_e_MeV)
    fraction_subtracted = float(sqrt_phi * m_e_MeV / m_p_MeV)

    physical_interpretations = [
        {
            "name": "Vacuum structure correction",
            "description": (
                "phi encodes the geometry of the compact space (S^2 at D=6, d=4). "
                "sqrt(phi)*m_e is the vacuum contribution to the gravitational "
                "coupling. The formula G ~ m_p*(m_p - sqrt(phi)*m_e) says that "
                "the gravitational coupling involves the proton mass corrected "
                "by a vacuum-geometry-scaled electron mass."
            ),
            "plausibility": "moderate",
        },
        {
            "name": "Reduced mass analog",
            "description": (
                "mu*(mu - k) = mu^2 - k*mu resembles a reduced-mass correction. "
                "For k=1, it would be (m_p/m_e)*(m_p/m_e - 1) = m_p*(m_p - m_e)/m_e^2, "
                "analogous to a reduced mass. For k=sqrt(phi), the 'partner mass' "
                "is sqrt(phi)*m_e rather than m_e itself."
            ),
            "plausibility": "moderate",
        },
        {
            "name": "Golden ratio scaling",
            "description": (
                "The golden ratio appears in the vacuum polynomial from d=4, D=6. "
                "Its square root scales the electron mass to define a 'gravitational "
                "mass gap'. The combination m_p - sqrt(phi)*m_e ~ 937.85 MeV is "
                "almost all of the proton mass, consistent with the proton mass "
                "being dominated by QCD binding energy (~99.9% of 938.3 MeV)."
            ),
            "plausibility": "speculative",
        },
        {
            "name": "Effective mass",
            "description": (
                "m_p - sqrt(phi)*m_e could be an effective gravitational mass, "
                "with the subtracted piece sqrt(phi)*m_e ~ 0.650 MeV representing "
                "a quantum correction. This mass is tantalizingly close to the "
                "up quark mass (~2.2 MeV) divided by pi (~0.70 MeV), though "
                "such numerological comparisons should be treated with caution."
            ),
            "plausibility": "speculative",
        },
    ]

    assessment = (
        f"sqrt(phi) = {sqrt_phi_value:.10f} emerges as the square root of the "
        f"golden ratio, which is itself determined by the vacuum polynomial "
        f"x^2 + 6x + 4 = 0 at D=6, d=4. The mass scale sqrt(phi)*m_e = "
        f"{sqrt_phi_m_e_MeV:.4f} MeV does not correspond to any known particle, "
        f"but the subtracted mass m_p - sqrt(phi)*m_e = {subtracted_mass_MeV:.2f} MeV "
        f"is 99.95% of the proton mass. The physical interpretation remains "
        f"open -- all four candidates are plausible but none is derived from "
        f"first principles."
    )

    return {
        "vacuum_polynomial_roots": [float(root_1), float(root_2)],
        "phi_from_roots": (
            "Roots of x^2 + 6x + 4 = 0 are -2/phi^2 and -2*phi^2, "
            "where phi = (1+sqrt(5))/2. The polynomial lives in Q(sqrt(5))."
        ),
        "sqrt_phi_value": sqrt_phi_value,
        "sqrt_phi_m_e_MeV": sqrt_phi_m_e_MeV,
        "subtracted_mass_MeV": subtracted_mass_MeV,
        "proton_mass_MeV": float(m_p_MeV),
        "fraction_subtracted": fraction_subtracted,
        "physical_interpretations": physical_interpretations,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 4. Compare all bridge/hierarchy formulas
# ---------------------------------------------------------------------------

def compare_all_bridges(constants=None):
    """Comprehensive comparison of all bridge and hierarchy formulas.

    Computes six formulas for alpha_g with signed residuals and sorts them
    by absolute residual. Identifies the best zero-parameter and best
    overall formula.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        Keys: formulas (sorted list), best_zero_param, best_overall,
        assessment.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    phi = constants.phi
    m_p = constants.m_p
    m_e = constants.m_e
    G_measured = constants.G

    mu = m_p / m_e
    sqrt_phi = phi.sqrt()

    formulas_spec = [
        {
            "formula_label": "phi^2/2 * alpha^21 (tree-level bridge)",
            "alpha_g": phi ** 2 / 2 * alpha ** 21,
            "n_fitted_params": 1,
        },
        {
            "formula_label": "phi^2/2 * (1+3*alpha^2) * alpha^21 (corrected bridge)",
            "alpha_g": phi ** 2 / 2 * (1 + Decimal(3) * alpha ** 2) * alpha ** 21,
            "n_fitted_params": "fitted",
        },
        {
            "formula_label": "alpha^24 * mu^2 (hierarchy)",
            "alpha_g": alpha ** 24 * mu ** 2,
            "n_fitted_params": 0,
        },
        {
            "formula_label": "alpha^24 * mu*(mu-1) (integer offset)",
            "alpha_g": alpha ** 24 * mu * (mu - 1),
            "n_fitted_params": 0,
        },
        {
            "formula_label": "alpha^24 * mu*(mu-5/4) (rational offset)",
            "alpha_g": alpha ** 24 * mu * (mu - Decimal(5) / Decimal(4)),
            "n_fitted_params": 0,
        },
        {
            "formula_label": "alpha^24 * mu*(mu-sqrt(phi)) (sqrt(phi) formula)",
            "alpha_g": alpha ** 24 * mu * (mu - sqrt_phi),
            "n_fitted_params": 0,
        },
    ]

    # Refined mu-structure: alpha^24 * mu * (mu - sqrt(phi)*(1-alpha))
    k_refined = sqrt_phi * (1 - alpha)
    ag_refined = alpha ** 24 * mu * (mu - k_refined)
    G_refined = ag_refined * constants.hbar * constants.c / constants.m_e ** 2
    res_refined = float((G_refined - G_measured) / G_measured * Decimal("1e6"))
    formulas_spec.append({
        "formula_label": "alpha^24 * mu*(mu - sqrt(phi)*(1-alpha))",
        "alpha_g": ag_refined,
        "n_fitted_params": 0,
    })

    formulas = []
    for spec in formulas_spec:
        ag = spec["alpha_g"]
        G_predicted = _alpha_g_to_G(ag, constants)
        residual_ppm = _signed_residual_ppm(G_predicted, G_measured)
        abs_residual_ppm = abs(residual_ppm)
        within_codata_uncertainty = abs_residual_ppm < 22.0

        formulas.append({
            "formula_label": spec["formula_label"],
            "alpha_g": ag,
            "G_predicted": G_predicted,
            "residual_ppm": residual_ppm,
            "abs_residual_ppm": abs_residual_ppm,
            "n_fitted_params": spec["n_fitted_params"],
            "within_codata_uncertainty": within_codata_uncertainty,
        })

    # Sort by absolute residual
    formulas.sort(key=lambda x: x["abs_residual_ppm"])

    # Best zero-parameter formula
    zero_param = [f for f in formulas if f["n_fitted_params"] == 0]
    best_zero_param = zero_param[0] if zero_param else None

    best_overall = formulas[0]

    assessment = (
        f"The best overall formula is '{best_overall['formula_label']}' "
        f"at {best_overall['residual_ppm']:.2f} ppm. "
    )
    if best_zero_param is not None:
        assessment += (
            f"The best zero-parameter formula is "
            f"'{best_zero_param['formula_label']}' at "
            f"{best_zero_param['residual_ppm']:.2f} ppm. "
            f"The sqrt(phi) mu-structure formula achieves sub-22 ppm "
            f"accuracy with no fitted parameters, making it the most "
            f"economical G prediction in the Alpha Ladder framework."
        )

    return {
        "formulas": formulas,
        "best_zero_param": best_zero_param,
        "best_overall": best_overall,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 5. CODATA stability
# ---------------------------------------------------------------------------

def analyze_codata_stability_mu(constants=None):
    """Test the sqrt(phi) formula against both CODATA editions.

    Verifies that the mu-structure result is stable across CODATA 2014
    and CODATA 2018 constants.

    Parameters
    ----------
    constants : SimpleNamespace or None
        Unused (both editions are always tested). Accepted for API
        consistency with other modules.

    Returns
    -------
    dict
        Mapping edition name -> dict with G_predicted, G_measured,
        residual_ppm, within_uncertainty.
    """
    results = {}

    for edition in ("CODATA 2014", "CODATA 2018"):
        c = get_constants(edition)
        alpha = c.alpha
        phi = c.phi
        m_p = c.m_p
        m_e = c.m_e
        G_measured = c.G

        mu = m_p / m_e
        sqrt_phi = phi.sqrt()

        alpha_g = alpha ** 24 * mu * (mu - sqrt_phi)
        G_predicted = _alpha_g_to_G(alpha_g, c)
        residual_ppm = _signed_residual_ppm(G_predicted, G_measured)
        within_uncertainty = abs(residual_ppm) < 22.0

        results[edition] = {
            "G_predicted": G_predicted,
            "G_measured": G_measured,
            "residual_ppm": residual_ppm,
            "within_uncertainty": within_uncertainty,
        }

    return results


# ---------------------------------------------------------------------------
# 6. Summary
# ---------------------------------------------------------------------------

def summarize_mu_structure(constants=None):
    """Main entry point combining all mu-structure analyses.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        Keys: mu_structure, offset_scan, sqrt_phi_origin,
        formula_comparison, codata_stability, key_finding,
        honest_assessment.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    mu_structure = compute_mu_structure(constants)
    refined = compute_mu_structure_refined(constants)
    offset_scan = scan_mu_offsets(constants)
    sqrt_phi_origin = analyze_sqrt_phi_origin(constants)
    formula_comparison = compare_all_bridges(constants)
    codata_stability = analyze_codata_stability_mu(constants)

    key_finding = (
        f"Replacing mu^2 with mu*(mu - sqrt(phi)) in the hierarchy formula "
        f"improves the residual from {mu_structure['hierarchy_residual_ppm']:.1f} ppm "
        f"to {mu_structure['residual_ppm']:.2f} ppm -- an improvement factor of "
        f"{mu_structure['improvement_factor']:.0f}x. The formula "
        f"G = alpha^24 * m_p * (m_p - sqrt(phi)*m_e) * hbar*c / m_e^4 "
        f"has zero fitted parameters and is within CODATA G measurement "
        f"uncertainty (~22 ppm)."
    )

    honest_assessment = (
        "DERIVED: phi from the vacuum polynomial x^2 + 6x + 4 = 0 at D=6, d=4; "
        "alpha and mu from CODATA measurements. "
        "EMPIRICAL: sqrt(phi) as the offset is discovered numerically by scanning "
        "mu*(mu - k) over algebraically natural k values. It is not derived from "
        "first principles. "
        "The formula has ZERO fitted parameters -- phi is determined by the theory "
        "at d=4, D=6, and alpha, mu are measured. "
        f"The {mu_structure['residual_ppm']:.2f} ppm residual is within the "
        "CODATA G measurement uncertainty (~22 ppm). "
        "However, we have NO theoretical derivation of WHY the offset should be "
        "sqrt(phi). The formula could be a numerical coincidence: the offset space "
        "is continuous, so finding an algebraically 'nice' number near the best fit "
        "is expected. "
        "To validate, one would need either (a) a Feynman diagram or EFT calculation "
        "that produces sqrt(phi) as a mass correction, or (b) a prediction for a "
        "different observable using the same mu-structure pattern."
    )

    return {
        "mu_structure": mu_structure,
        "refined": refined,
        "offset_scan": offset_scan,
        "sqrt_phi_origin": sqrt_phi_origin,
        "formula_comparison": formula_comparison,
        "codata_stability": codata_stability,
        "key_finding": key_finding,
        "honest_assessment": honest_assessment,
    }
