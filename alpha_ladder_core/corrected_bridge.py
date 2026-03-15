"""
Radiative correction to the phi^2/2 bridge coefficient.

Discovery
---------
The bridge coefficient C in alpha_g = C * alpha^21 was found to be:

    C_exact = alpha_g / alpha^21 = 1.3092269314... (CODATA 2018)
    phi^2/2 = 1.3090169944...                      (160 ppm residual)
    phi^2/2 * (1 + 3*alpha^2) = 1.3092261152...    (0.62 ppm residual)

The multiplicative correction (1 + 3*alpha^2) closes 99.6% of the 160 ppm
gap between phi^2/2 and the measured bridge coefficient. The factor 3 has
a remarkable degeneracy: for the specific case d=4 spacetime dimensions and
n=2 extra compact dimensions (S^2), the number 3 simultaneously equals:

    d - 1   = 3   (spatial dimensions)
    n(n+1)/2 = 3  (SO(3) isometry generators of S^2)
    n + 1    = 3  (extra-dimensional + 1)

This triple coincidence is specific to n=2 and does not persist for other
compactification dimensions.

What is derived vs empirical
-----------------------------
- C_exact is measured (from G, m_e, hbar, c, alpha via CODATA).
- phi^2/2 as the tree-level bridge is motivated by the vacuum polynomial
  living in Q(sqrt(5)), but ultimately selected by numerical search.
- The correction 3*alpha^2 is empirically observed to close the gap.
- The interpretation of "3" as (d-1), n(n+1)/2, or (n+1) is speculation;
  all three happen to coincide for d=4, n=2.
- Higher-order coefficients (c_3, c_4, ...) are purely fitted.
"""

from decimal import Decimal, getcontext

getcontext().prec = 50

from alpha_ladder_core.constants import get_constants


# ---------------------------------------------------------------------------
# 1. Compute corrected bridge
# ---------------------------------------------------------------------------

def compute_corrected_bridge(constants=None):
    """Compute the corrected bridge coefficient and compare to C_exact.

    The corrected bridge is phi^2/2 * (1 + 3*alpha^2), which closes 99.6%
    of the gap between phi^2/2 and the measured C_exact.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        Keys: C_exact, C_uncorrected, C_corrected_leading, C_corrected_nlo,
        uncorrected_ppm, corrected_leading_ppm, corrected_nlo_ppm,
        G_uncorrected, G_corrected, G_measured, correction_factor,
        fraction_of_gap_explained.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    phi = constants.phi
    hbar = constants.hbar
    c = constants.c
    m_e = constants.m_e
    alpha_g = constants.alpha_g
    G_measured = constants.G

    alpha_21 = alpha ** 21

    # Exact bridge from measurement
    C_exact = alpha_g / alpha_21

    # Uncorrected: phi^2/2
    C_uncorrected = phi ** 2 / 2

    # Leading correction: phi^2/2 * (1 + 3*alpha^2)
    correction_factor = Decimal(3) * alpha ** 2
    C_corrected_leading = C_uncorrected * (1 + correction_factor)

    # NLO correction: phi^2/2 * (1 + alpha^2*(3 + 8/5*alpha))
    nlo_factor = alpha ** 2 * (Decimal(3) + Decimal(8) / Decimal(5) * alpha)
    C_corrected_nlo = C_uncorrected * (1 + nlo_factor)

    # Residuals in ppm
    uncorrected_ppm = float(abs(C_uncorrected - C_exact) / C_exact * Decimal("1e6"))
    corrected_leading_ppm = float(abs(C_corrected_leading - C_exact) / C_exact * Decimal("1e6"))
    corrected_nlo_ppm = float(abs(C_corrected_nlo - C_exact) / C_exact * Decimal("1e6"))

    # G predictions
    G_uncorrected = C_uncorrected * alpha_21 * hbar * c / m_e ** 2
    G_corrected = C_corrected_leading * alpha_21 * hbar * c / m_e ** 2

    # Fraction of gap explained
    gap_original = float(abs(C_uncorrected - C_exact))
    gap_remaining = float(abs(C_corrected_leading - C_exact))
    fraction_of_gap_explained = (
        1.0 - gap_remaining / gap_original if gap_original > 0 else 0.0
    )

    return {
        "C_exact": C_exact,
        "C_uncorrected": C_uncorrected,
        "C_corrected_leading": C_corrected_leading,
        "C_corrected_nlo": C_corrected_nlo,
        "uncorrected_ppm": uncorrected_ppm,
        "corrected_leading_ppm": corrected_leading_ppm,
        "corrected_nlo_ppm": corrected_nlo_ppm,
        "G_uncorrected": G_uncorrected,
        "G_corrected": G_corrected,
        "G_measured": G_measured,
        "correction_factor": correction_factor,
        "fraction_of_gap_explained": fraction_of_gap_explained,
    }


# ---------------------------------------------------------------------------
# 2. Analyze correction origin
# ---------------------------------------------------------------------------

def analyze_correction_origin(d=4, n=2, constants=None):
    """Explore what the factor 3 in 3*alpha^2 could mean physically.

    Tests several interpretations of the integer factor and checks which
    ones match for d=4, n=2 and how they generalize.

    Parameters
    ----------
    d : int
        Number of spacetime dimensions (default 4).
    n : int
        Number of extra compact dimensions (default 2).
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        Keys: factor_value, interpretations, degeneracy_note, assessment.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    alpha_f = float(alpha)

    D = d + n  # total dimensions

    interpretations = []

    # 1. Spatial dimensions: d - 1
    val_spatial = d - 1
    interpretations.append({
        "name": "Spatial dimensions (d-1)",
        "formula": f"d - 1 = {d} - 1",
        "value_for_d4_n2": val_spatial,
        "matches": val_spatial == 3,
        "generalization_note": (
            f"For general d: d-1. At d=4: 3. At d=5: 4. "
            f"Tied to spacetime dimension, independent of n."
        ),
    })

    # 2. SU(2) generators / SO(3) isometry of S^n
    val_isometry = n * (n + 1) // 2
    interpretations.append({
        "name": "SO(n+1) isometry generators of S^n",
        "formula": f"n(n+1)/2 = {n}*{n+1}/2",
        "value_for_d4_n2": val_isometry,
        "matches": val_isometry == 3,
        "generalization_note": (
            f"For general n: n(n+1)/2 (dimension of SO(n+1)). "
            f"At n=2: 3. At n=3: 6. At n=4: 10. "
            f"Counts gauge bosons from KK reduction on S^n."
        ),
    })

    # 3. SU(3) color
    interpretations.append({
        "name": "SU(3) color (N_c)",
        "formula": "N_c = 3",
        "value_for_d4_n2": 3,
        "matches": True,
        "generalization_note": (
            "Fixed at 3 for QCD. Does not depend on d or n. "
            "Would need a mechanism linking color to the gravitational "
            "bridge coefficient."
        ),
    })

    # 4. D - d + 1 = n + 1
    val_extra = n + 1
    interpretations.append({
        "name": "Extra dimensions + 1 (n+1 = D-d+1)",
        "formula": f"n + 1 = {n} + 1 = D - d + 1 = {D} - {d} + 1",
        "value_for_d4_n2": val_extra,
        "matches": val_extra == 3,
        "generalization_note": (
            f"For general n: n+1. At n=2: 3. At n=3: 4. At n=4: 5. "
            f"Could relate to KK tower structure."
        ),
    })

    # 5. QED vertex correction alpha/(2*pi)
    val_qed = alpha_f / (2 * 3.14159265358979)
    interpretations.append({
        "name": "One-loop QED vertex correction",
        "formula": "alpha/(2*pi)",
        "value_for_d4_n2": round(val_qed, 8),
        "matches": False,
        "generalization_note": (
            f"alpha/(2*pi) = {val_qed:.6e}, which is ~0.00116. "
            f"Much larger than 3*alpha^2 = {3 * alpha_f**2:.6e}. "
            f"Wrong functional form (linear in alpha, not quadratic)."
        ),
    })

    # 6. Graviton polarizations
    # Massless graviton in d dims: d(d-3)/2 polarizations
    # Massive KK graviton: d(d-1)/2 - 1 polarizations
    pol_massless = d * (d - 3) // 2
    pol_massive = d * (d - 1) // 2 - 1
    interpretations.append({
        "name": "Graviton polarizations",
        "formula": f"massless: d(d-3)/2 = {pol_massless}; massive KK: d(d-1)/2 - 1 = {pol_massive}",
        "value_for_d4_n2": pol_massless,
        "matches": pol_massless == 3 or pol_massive == 3,
        "generalization_note": (
            f"In d=4: massless graviton has {pol_massless} polarizations, "
            f"massive KK graviton has {pol_massive}. Neither equals 3 cleanly "
            f"(massless = 2, massive = 5)."
        ),
    })

    # 7. d*n KK vectors
    val_kk_vectors = d * n
    interpretations.append({
        "name": "KK vector fields (d*n)",
        "formula": f"d * n = {d} * {n}",
        "value_for_d4_n2": val_kk_vectors,
        "matches": val_kk_vectors == 3,
        "generalization_note": (
            f"d*n = {val_kk_vectors} KK vectors. Does not equal 3 for d=4, n=2 "
            f"(gives 8). Not the right counting."
        ),
    })

    # Check divergence for n=3
    divergence_examples = {
        "n=2": {"d-1": d - 1, "n(n+1)/2": 2 * 3 // 2, "n+1": 3},
        "n=3": {"d-1": d - 1, "n(n+1)/2": 3 * 4 // 2, "n+1": 4},
        "n=4": {"d-1": d - 1, "n(n+1)/2": 4 * 5 // 2, "n+1": 5},
    }

    degeneracy_note = (
        "For n=2, three distinct physical interpretations -- (d-1)=3 spatial "
        "dimensions, n(n+1)/2=3 SO(3) isometry generators, and (n+1)=3 -- all "
        "yield the same integer 3. This is a coincidence specific to n=2. "
        f"For n=3 they give {divergence_examples['n=3']}, and for n=4 they "
        f"give {divergence_examples['n=4']}. Without an independent way to "
        "determine n (it is fixed at 2 by the theory), we cannot distinguish "
        "these interpretations from the bridge coefficient alone."
    )

    assessment = (
        "The most physically motivated interpretation depends on the "
        "theoretical framework. In the Alpha Ladder's 6D KK context, "
        "n(n+1)/2 = dim(SO(3)) is natural because the S^2 isometry group "
        "generates the KK gauge bosons that contribute to loop corrections. "
        "However, (d-1) = 3 spatial dimensions is the simplest explanation "
        "and (n+1) = 3 also has a natural KK interpretation. The three-way "
        "degeneracy at n=2 means the data cannot distinguish these options."
    )

    return {
        "factor_value": 3,
        "interpretations": interpretations,
        "degeneracy_note": degeneracy_note,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 3. Scan correction series
# ---------------------------------------------------------------------------

def scan_correction_series(constants=None, max_order=5):
    """Expand the correction as a power series in alpha.

    Determines the coefficients c_k in:
        C = phi^2/2 * (1 + sum_{k=2}^{max_order} c_k * alpha^k)

    At each order, c_k is determined by dividing the remaining fractional
    residual by alpha^k, with all lower-order coefficients fixed.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.
    max_order : int
        Maximum power of alpha to include (default 5).

    Returns
    -------
    dict
        Keys: coefficients, series_converging, best_order, assessment.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    phi = constants.phi
    alpha_g = constants.alpha_g
    G_measured = constants.G

    alpha_21 = alpha ** 21
    C_exact = alpha_g / alpha_21
    C_tree = phi ** 2 / 2

    # Fractional correction: C_exact = C_tree * (1 + delta)
    delta = C_exact / C_tree - 1

    # G measurement uncertainty in ppm (CODATA 2018: 6.67430(15)e-11)
    # Relative uncertainty = 0.00015 / 6.67430 ~ 22 ppm
    G_unc_ppm = 22.0

    coefficients = []
    remainder = delta

    for k in range(2, max_order + 1):
        alpha_k = alpha ** k

        # Determine c_k
        c_k_decimal = remainder / alpha_k
        c_k_float = float(c_k_decimal)

        # Round to nearest integer to check if it is close to a simple number
        c_k_rounded = round(c_k_float)

        # For k=2, we expect c_2 = 3; use the rounded value
        # For higher orders, keep the float
        if k == 2:
            c_k_used = Decimal(c_k_rounded)
        else:
            c_k_used = c_k_decimal

        # Update remainder
        remainder = remainder - c_k_used * alpha_k

        # Residual after this order
        C_approx = C_tree * (1 + delta - remainder)
        residual_after_ppm = float(abs(C_approx - C_exact) / C_exact * Decimal("1e6"))

        coefficients.append({
            "order": k,
            "coefficient": float(c_k_used),
            "coefficient_exact": c_k_float,
            "nearest_integer": c_k_rounded,
            "residual_after_ppm": residual_after_ppm,
        })

    # Check convergence: |c_{k+1} * alpha^{k+1}| < |c_k * alpha^k|
    alpha_f = float(alpha)
    series_converging = True
    for i in range(len(coefficients) - 1):
        term_k = abs(coefficients[i]["coefficient"]) * alpha_f ** coefficients[i]["order"]
        term_k1 = abs(coefficients[i + 1]["coefficient"]) * alpha_f ** coefficients[i + 1]["order"]
        if term_k1 > term_k:
            series_converging = False
            break

    # Find best order (where residual drops below G uncertainty)
    best_order = max_order
    for entry in coefficients:
        if entry["residual_after_ppm"] < G_unc_ppm:
            best_order = entry["order"]
            break

    c2 = coefficients[0]["coefficient"] if coefficients else 0
    c2_exact = coefficients[0]["coefficient_exact"] if coefficients else 0

    assessment = (
        f"The leading correction coefficient c_2 = {c2_exact:.6f}, which rounds "
        f"to {coefficients[0]['nearest_integer']} (used as {c2}). "
        f"After the leading correction, the residual drops from "
        f"~160 ppm to ~{coefficients[0]['residual_after_ppm']:.2f} ppm. "
        f"The series {'converges' if series_converging else 'does not clearly converge'} "
        f"in the sense that successive terms decrease. "
        f"The residual drops below the CODATA G uncertainty ({G_unc_ppm} ppm) "
        f"at order {best_order}."
    )

    return {
        "coefficients": coefficients,
        "series_converging": series_converging,
        "best_order": best_order,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 4. Compare bridge to hierarchy
# ---------------------------------------------------------------------------

def compare_bridge_hierarchy(constants=None):
    """Compare the corrected bridge formula to the hierarchy formula.

    Corrected bridge: alpha_g = phi^2/2 * (1 + 3*alpha^2) * alpha^21
    Hierarchy:        alpha_g = alpha^24 * mu^2

    If both are valid, then phi^2/2 * (1 + 3*alpha^2) = alpha^3 * mu^2.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        Keys: bridge_alpha_g, hierarchy_alpha_g, ratio, difference_ppm,
        epsilon, interpretation.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    phi = constants.phi
    m_p = constants.m_p
    m_e = constants.m_e

    mu = m_p / m_e

    # Bridge formula
    C_corrected = phi ** 2 / 2 * (1 + Decimal(3) * alpha ** 2)
    bridge_alpha_g = C_corrected * alpha ** 21

    # Hierarchy formula
    hierarchy_alpha_g = alpha ** 24 * mu ** 2

    # Ratio and difference
    ratio = float(bridge_alpha_g / hierarchy_alpha_g)
    difference_ppm = float(
        abs(bridge_alpha_g - hierarchy_alpha_g) / hierarchy_alpha_g * Decimal("1e6")
    )

    # Epsilon: phi^2/2 * (1 + 3*alpha^2) = alpha^3 * mu^2 * (1 + epsilon)
    lhs = C_corrected
    rhs = alpha ** 3 * mu ** 2
    epsilon = float(lhs / rhs - 1)

    if abs(difference_ppm) < 1000:
        interpretation = (
            f"The corrected bridge and hierarchy formulae agree to "
            f"{difference_ppm:.1f} ppm. The relationship "
            f"phi^2/2 * (1 + 3*alpha^2) = alpha^3 * mu^2 * (1 + epsilon) "
            f"with epsilon = {epsilon:.6e} suggests these are two "
            f"representations of the same underlying structure. The bridge "
            f"formula (with phi) may be the 'tree-level + radiative correction' "
            f"form, while the hierarchy formula (with mu) may be the exact "
            f"all-orders result expressed in terms of mass ratios."
        )
    else:
        interpretation = (
            f"The corrected bridge and hierarchy formulae differ by "
            f"{difference_ppm:.1f} ppm. The residual epsilon = {epsilon:.6e} "
            f"is non-negligible, suggesting higher-order corrections in the "
            f"bridge formula or that the two formulae capture different physics."
        )

    return {
        "bridge_alpha_g": bridge_alpha_g,
        "hierarchy_alpha_g": hierarchy_alpha_g,
        "ratio": ratio,
        "difference_ppm": difference_ppm,
        "epsilon": epsilon,
        "interpretation": interpretation,
    }


# ---------------------------------------------------------------------------
# 5. Analyze CODATA editions
# ---------------------------------------------------------------------------

def analyze_codata_editions():
    """Run the corrected bridge against both CODATA 2014 and 2018.

    Tests whether the correction is stable across editions.

    Returns
    -------
    dict
        Mapping edition name -> dict with C_exact, uncorrected_ppm,
        corrected_ppm, G_predicted, G_measured.
    """
    results = {}

    for edition in ("CODATA 2014", "CODATA 2018"):
        c = get_constants(edition)
        alpha = c.alpha
        phi = c.phi
        hbar = c.hbar
        c_light = c.c
        m_e = c.m_e
        alpha_g = c.alpha_g
        G_measured = c.G

        alpha_21 = alpha ** 21
        C_exact = alpha_g / alpha_21
        C_tree = phi ** 2 / 2
        C_corrected = C_tree * (1 + Decimal(3) * alpha ** 2)

        uncorrected_ppm = float(abs(C_tree - C_exact) / C_exact * Decimal("1e6"))
        corrected_ppm = float(abs(C_corrected - C_exact) / C_exact * Decimal("1e6"))

        G_predicted = C_corrected * alpha_21 * hbar * c_light / m_e ** 2

        results[edition] = {
            "C_exact": C_exact,
            "uncorrected_ppm": uncorrected_ppm,
            "corrected_ppm": corrected_ppm,
            "G_predicted": G_predicted,
            "G_measured": G_measured,
        }

    return results


# ---------------------------------------------------------------------------
# 6. Summary
# ---------------------------------------------------------------------------

def summarize_corrected_bridge(constants=None):
    """Main entry point combining all analyses.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants. Defaults to CODATA 2018.

    Returns
    -------
    dict
        Keys: corrected_bridge, correction_origin, series,
        bridge_vs_hierarchy, editions, key_finding, honest_assessment.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    corrected = compute_corrected_bridge(constants)
    origin = analyze_correction_origin(d=4, n=2, constants=constants)
    series = scan_correction_series(constants)
    hierarchy = compare_bridge_hierarchy(constants)
    editions = analyze_codata_editions()

    key_finding = (
        "The bridge coefficient C in alpha_g = C * alpha^21 is within 160 ppm "
        "of phi^2/2, but applying a multiplicative correction (1 + 3*alpha^2) "
        f"closes 99.6% of this gap, reducing the residual to "
        f"{corrected['corrected_leading_ppm']:.2f} ppm. The factor 3 has a "
        "remarkable triple degeneracy at n=2: it simultaneously equals (d-1) "
        "spatial dimensions, n(n+1)/2 SO(3) isometry generators of S^2, and "
        "(n+1). This correction has the form of a one-loop radiative correction "
        "at order alpha^2, consistent with the bridge coefficient receiving "
        "quantum corrections from the KK gauge sector."
    )

    honest_assessment = (
        "DERIVED: C_exact from CODATA measurements; phi^2/2 as tree-level "
        "bridge from algebraic search in Q(sqrt(5)). "
        "EMPIRICAL: The correction 3*alpha^2 is numerically observed to close "
        "the gap. The coefficient 3 is not derived from first principles. "
        "SPECULATIVE: The interpretation of 3 as (d-1), n(n+1)/2, or (n+1) "
        "is post-hoc. The three-way degeneracy at n=2 prevents discrimination. "
        "Higher-order coefficients (c_3, c_4, ...) are purely fitted and "
        "their physical meaning is unknown. The correction looks like a "
        "radiative correction, but no Feynman diagram calculation has been "
        "performed to confirm this."
    )

    return {
        "corrected_bridge": corrected,
        "correction_origin": origin,
        "series": series,
        "bridge_vs_hierarchy": hierarchy,
        "editions": editions,
        "key_finding": key_finding,
        "honest_assessment": honest_assessment,
    }
