"""
Second predictions beyond Newton's G from the Alpha Ladder framework.

The Alpha Ladder formula:
    G = alpha^24 * mu * (mu - sqrt(phi)*(1-alpha)) * hbar*c / m_e^2

predicts Newton's G to 0.33 ppm with zero fitted parameters.  For scientific
credibility the framework needs SECOND predictions -- quantities other than G
that can be independently tested.

This module computes three classes of second predictions:

1. **Time variation of G**: If fundamental constants drift, the formula
   predicts dG/G = A*(dalpha/alpha) + B*(dmu/mu) with A ~ 24.0007.
   The coefficient 24 = d*D is unique to this framework.

2. **M_6 (6D Planck mass) constraints**: The KK relation
   M_Pl^2 = M_6^4 * 4*pi*a_0^2 combined with the fifth-force survival
   window a_0 in [30, 71] um places M_6 in [3.1, 4.8] TeV (LHC scale).

3. **Correlated variations**: If alpha varies, the framework predicts
   correlated G variation with a specific ratio distinguishable from
   competing frameworks.

All calculations use pure Python + math + decimal (no numpy/scipy).
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50

from alpha_ladder_core.constants import get_constants


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_float_constants(constants):
    """Extract float versions of the constants needed by this module."""
    if constants is None:
        constants = get_constants("CODATA 2018")
    alpha = float(constants.alpha)
    mu = float(constants.m_p / constants.m_e)
    phi = float(constants.phi)
    sqrt_phi = math.sqrt(phi)
    hbar = float(constants.hbar)
    cc = float(constants.c)
    m_e = float(constants.m_e)
    G = float(constants.G)
    M_Pl = math.sqrt(hbar * cc / G)       # kg
    l_Pl = math.sqrt(hbar * G / cc ** 3)   # meters
    return {
        "alpha": alpha, "mu": mu, "phi": phi, "sqrt_phi": sqrt_phi,
        "hbar": hbar, "c": cc, "m_e": m_e, "G": G,
        "M_Pl": M_Pl, "l_Pl": l_Pl, "constants": constants,
    }


def _log_spaced(x_min, x_max, n):
    """Return a list of n log-spaced values from x_min to x_max."""
    log_min = math.log10(x_min)
    log_max = math.log10(x_max)
    return [10 ** (log_min + i * (log_max - log_min) / (n - 1))
            for i in range(n)]


# ---------------------------------------------------------------------------
# 1. Time variation coefficients
# ---------------------------------------------------------------------------

def compute_time_variation_coefficients(constants=None):
    """Compute the log-derivative coefficients dG/G = A*(dalpha/alpha) + B*(dmu/mu).

    From G = alpha^24 * mu * (mu - k) * hbar*c / m_e^2 where k = sqrt(phi)*(1-alpha):
        ln G = 24*ln(alpha) + ln(mu) + ln(mu - k) + const
    Differentiating:
        A = d(lnG)/d(ln alpha) = 24 + alpha*dk/dalpha / (mu - k)
          = 24 + sqrt(phi)*alpha / (mu - k)     [since dk/dalpha = -sqrt(phi)]
          but the sign: k = sqrt(phi)*(1-alpha), dk/dalpha = -sqrt(phi)
          d(ln(mu-k))/d(ln alpha) = alpha * sqrt(phi) / (mu - k)
          (positive, because increasing alpha decreases k, increasing mu-k)
        B = d(lnG)/d(ln mu) = 1 + mu/(mu - k) = (2*mu - k) / (mu - k)

    Returns
    -------
    dict
        A_coefficient, B_coefficient, formula_description,
        comparison_with_other_frameworks, assessment.
    """
    p = _get_float_constants(constants)
    alpha = p["alpha"]
    mu = p["mu"]
    sqrt_phi = p["sqrt_phi"]

    k = sqrt_phi * (1.0 - alpha)
    mu_minus_k = mu - k

    # A = 24 + alpha * sqrt(phi) / (mu - k)
    A = 24.0 + alpha * sqrt_phi / mu_minus_k

    # B = (2*mu - k) / (mu - k)
    B = (2.0 * mu - k) / mu_minus_k

    comparison = {
        "Alpha Ladder (d*D=24)": A,
        "Standard KK (5D, G~1/R)": 1.0,
        "Bekenstein varying-alpha": 0.0,
        "Dirac Large Numbers": 1.0,
        "Randall-Sundrum": 2.0,
    }

    return {
        "A_coefficient": A,
        "B_coefficient": B,
        "formula_description": (
            "dG/G = {:.4f} * (dalpha/alpha) + {:.4f} * (dmu/mu).  "
            "The dominant coefficient 24 = d*D comes from alpha^24 in the "
            "hierarchy formula and is unique to the Alpha Ladder framework."
        ).format(A, B),
        "comparison_with_other_frameworks": comparison,
        "assessment": (
            "The coefficient A ~ 24 is unique among all known gravity "
            "frameworks.  Standard KK and Dirac give A=1; Randall-Sundrum "
            "gives A=2; Bekenstein gives A=0.  A measurement of dG/G "
            "correlated with dalpha/alpha at the 24:1 ratio would be strong "
            "evidence for the Alpha Ladder."
        ),
    }


# ---------------------------------------------------------------------------
# 2. Current experimental bounds
# ---------------------------------------------------------------------------

def compute_current_bounds(constants=None):
    """Compare framework predictions with current experimental bounds.

    Current experimental limits:
        |dalpha/alpha| < 1e-17 /yr  (atomic clocks)
        |dmu/mu|      < 1e-16 /yr  (molecular spectroscopy)
        |dG/G|        < 4e-14 /yr  (lunar laser ranging)

    Returns
    -------
    dict
        alpha_bound, mu_bound, G_bound, predicted_dG_from_alpha_bound,
        ratio_to_current_sensitivity, testable, years_to_testability,
        assessment.
    """
    tv = compute_time_variation_coefficients(constants)
    A = tv["A_coefficient"]

    alpha_bound = 1e-17    # per year
    mu_bound = 1e-16       # per year
    G_bound = 4e-14        # per year

    # If dalpha/alpha = alpha_bound and dmu/mu = 0:
    predicted_dG = A * alpha_bound   # ~ 2.4e-16 /yr

    ratio = predicted_dG / G_bound   # ~ 0.006 => well below

    return {
        "alpha_bound": alpha_bound,
        "mu_bound": mu_bound,
        "G_bound": G_bound,
        "predicted_dG_from_alpha_bound": predicted_dG,
        "ratio_to_current_sensitivity": ratio,
        "testable": False,
        "years_to_testability": 30,
        "assessment": (
            "The predicted dG/G ~ {:.1e} /yr from the alpha clock bound is "
            "{:.0f}x below current LLR sensitivity ({:.0e} /yr).  "
            "Testing requires roughly two orders of magnitude improvement in "
            "LLR or equivalent G-monitoring experiments, plausibly achievable "
            "in ~30 years with space-based laser ranging."
        ).format(predicted_dG, 1.0 / ratio, G_bound),
    }


# ---------------------------------------------------------------------------
# 3. M_6 from a_0
# ---------------------------------------------------------------------------

def compute_m6_from_a0(a0_meters, constants=None):
    """Compute the 6D Planck mass M_6 from the internal radius a_0.

    KK relation (natural units): M_Pl^2 = M_6^4 * 4*pi*a_0^2
    => M_6 = (M_Pl^2 / (4*pi*a_0^2))^{1/4}

    Parameters
    ----------
    a0_meters : float
        Internal radius in meters.

    Returns
    -------
    dict
        a0_meters, M_6_kg, M_6_eV, M_6_TeV, ratio_to_electron_mass,
        ratio_to_proton_mass, nearest_known_scale, alpha_power, assessment.
    """
    p = _get_float_constants(constants)
    hbar = p["hbar"]
    cc = p["c"]
    M_Pl = p["M_Pl"]

    # Convert M_Pl to natural units (1/meters)
    M_Pl_nat = M_Pl * cc / hbar

    # KK relation in natural units
    M_6_nat4 = M_Pl_nat ** 2 / (4.0 * math.pi * a0_meters ** 2)
    M_6_nat = M_6_nat4 ** 0.25

    # Convert back to kg and eV
    M_6_kg = M_6_nat * hbar / cc
    eV_per_J = 1.0 / 1.602176634e-19
    M_6_eV = M_6_kg * cc ** 2 * eV_per_J
    M_6_TeV = M_6_eV / 1e12

    # Reference scales
    m_e_eV = 0.51099895e6
    m_p_eV = 938.272046e6
    ratio_electron = M_6_eV / m_e_eV
    ratio_proton = M_6_eV / m_p_eV

    # Known particle masses in eV for comparison
    known_scales = {
        "electron": 0.511e6,
        "proton": 938.3e6,
        "W boson": 80.4e9,
        "Higgs boson": 125.1e9,
        "top quark": 173.0e9,
        "LHC reach": 14.0e12,
        "Planck mass": 1.22089e28,
    }

    # Find nearest known scale
    nearest_name = None
    nearest_ratio = float("inf")
    for name, scale_eV in known_scales.items():
        r = abs(math.log10(M_6_eV / scale_eV))
        if r < nearest_ratio:
            nearest_ratio = r
            nearest_name = name

    # alpha_power: what power of alpha gives M_6 from m_e?
    alpha = p["alpha"]
    if M_6_eV > 0 and m_e_eV > 0:
        alpha_power = math.log(M_6_eV / m_e_eV) / math.log(alpha)
    else:
        alpha_power = float("nan")

    return {
        "a0_meters": a0_meters,
        "M_6_kg": M_6_kg,
        "M_6_eV": M_6_eV,
        "M_6_TeV": M_6_TeV,
        "ratio_to_electron_mass": ratio_electron,
        "ratio_to_proton_mass": ratio_proton,
        "nearest_known_scale": nearest_name,
        "alpha_power": alpha_power,
        "assessment": (
            "For a_0 = {:.2e} m, M_6 = {:.2f} TeV.  Nearest known scale: {}.  "
            "M_6 = m_e * alpha^({:.2f})."
        ).format(a0_meters, M_6_TeV, nearest_name, alpha_power),
    }


# ---------------------------------------------------------------------------
# 4. Scan a_0 - M_6 landscape
# ---------------------------------------------------------------------------

def scan_a0_m6_landscape(constants=None):
    """Scan a_0 from l_Pl to 1 mm and compute M_6 at each point.

    Marks the survival window [30, 71] um and excluded regions.

    Returns
    -------
    dict
        scan_results (list of dicts), survival_window, assessment.
    """
    p = _get_float_constants(constants)
    l_Pl = p["l_Pl"]

    a0_values = _log_spaced(l_Pl, 1e-3, 20)

    # Survival window bounds
    a0_min_surv = 30e-6   # 30 um
    a0_max_surv = 71e-6   # 71 um

    scan_results = []
    survival_entries = []

    for a0 in a0_values:
        m6 = compute_m6_from_a0(a0, constants)
        in_survival = a0_min_surv <= a0 <= a0_max_surv
        excluded_eot_wash = a0 > a0_max_surv
        unphysical = a0 < l_Pl

        entry = {
            "a0_meters": a0,
            "M_6_TeV": m6["M_6_TeV"],
            "M_6_eV": m6["M_6_eV"],
            "in_survival_window": in_survival,
            "excluded_eot_wash": excluded_eot_wash,
            "unphysical": unphysical,
        }
        scan_results.append(entry)
        if in_survival:
            survival_entries.append(entry)

    # Compute survival window M_6 bounds from the actual window edges
    m6_at_min = compute_m6_from_a0(a0_min_surv, constants)
    m6_at_max = compute_m6_from_a0(a0_max_surv, constants)

    # Larger a_0 => smaller M_6
    survival_window = {
        "a0_min": a0_min_surv,
        "a0_max": a0_max_surv,
        "M6_min_TeV": m6_at_max["M_6_TeV"],   # larger a0 -> smaller M6
        "M6_max_TeV": m6_at_min["M_6_TeV"],   # smaller a0 -> larger M6
    }

    return {
        "scan_results": scan_results,
        "survival_window": survival_window,
        "assessment": (
            "The survival window a_0 in [30, 71] um maps to M_6 in "
            "[{:.1f}, {:.1f}] TeV.  This is at the LHC energy frontier, "
            "making M_6 a genuinely testable prediction.  Values of a_0 "
            "above 71 um are excluded by Eot-Wash; values below l_Pl are "
            "unphysical."
        ).format(survival_window["M6_min_TeV"], survival_window["M6_max_TeV"]),
    }


# ---------------------------------------------------------------------------
# 5. Correlated variations
# ---------------------------------------------------------------------------

def analyze_correlated_variations(dalpha_alpha=1e-17, constants=None):
    """Predict correlated G variation given a dalpha/alpha drift rate.

    If QCD-dominated coupling relates mu to alpha via dmu/mu ~ R * dalpha/alpha
    with R ~ 38 (from QCD Lambda dependence), then:
        dG/G = (A + B*R) * dalpha/alpha

    Parameters
    ----------
    dalpha_alpha : float
        Rate of alpha variation in per-year units.

    Returns
    -------
    dict
        dalpha_alpha, predicted_dG_G, predicted_dG_G_with_qcd,
        R_qcd, violates_llr_bound, assessment.
    """
    tv = compute_time_variation_coefficients(constants)
    A = tv["A_coefficient"]
    B = tv["B_coefficient"]

    R_qcd = 38.0

    predicted_dG_G = A * dalpha_alpha
    predicted_dG_G_with_qcd = (A + B * R_qcd) * dalpha_alpha

    G_bound = 4e-14  # LLR bound per year

    violates_llr = abs(predicted_dG_G_with_qcd) > G_bound

    return {
        "dalpha_alpha": dalpha_alpha,
        "predicted_dG_G": predicted_dG_G,
        "predicted_dG_G_with_qcd": predicted_dG_G_with_qcd,
        "R_qcd": R_qcd,
        "violates_llr_bound": violates_llr,
        "assessment": (
            "For dalpha/alpha = {:.1e} /yr: dG/G = {:.1e} /yr (alpha only) "
            "or {:.1e} /yr (with QCD coupling R={}).  The QCD-enhanced "
            "prediction is {:.0f}x larger.  {}"
        ).format(
            dalpha_alpha,
            predicted_dG_G,
            predicted_dG_G_with_qcd,
            R_qcd,
            abs(predicted_dG_G_with_qcd / predicted_dG_G) if predicted_dG_G != 0 else 0,
            "This violates current LLR bounds." if violates_llr
            else "This is still below current LLR sensitivity."
        ),
    }


# ---------------------------------------------------------------------------
# 6. Mu tension analysis
# ---------------------------------------------------------------------------

def analyze_mu_tension(constants=None):
    """Analyze the tension between bridge and mu-structure formulas.

    Bridge:       alpha_g = phi^2/2 * (1+3*alpha^2) * alpha^21
    Mu-structure: alpha_g = alpha^24 * mu * (mu - sqrt(phi)*(1-alpha))

    If both are exact:
        phi^2/2 * (1+3*alpha^2) = alpha^3 * mu * (mu - sqrt(phi)*(1-alpha))

    Solve for mu_predicted: quadratic mu^2 - sqrt(phi)*(1-alpha)*mu - C = 0
    where C = phi^2/2 * (1+3*alpha^2) / alpha^3.

    Returns
    -------
    dict
        mu_measured, mu_predicted_from_bridge, tension_ppm,
        tension_sigma, interpretation, possible_resolutions, assessment.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = float(constants.alpha)
    phi = float(constants.phi)
    sqrt_phi = math.sqrt(phi)
    mu_measured = float(constants.m_p / constants.m_e)

    # C = phi^2/2 * (1+3*alpha^2) / alpha^3
    C = phi ** 2 / 2.0 * (1.0 + 3.0 * alpha ** 2) / alpha ** 3

    # Quadratic: mu^2 - sqrt(phi)*(1-alpha)*mu - C = 0
    k = sqrt_phi * (1.0 - alpha)
    discriminant = k ** 2 + 4.0 * C
    mu_predicted = (k + math.sqrt(discriminant)) / 2.0

    # Tension
    tension_ppm = (mu_predicted - mu_measured) / mu_measured * 1e6
    # CODATA mu uncertainty is ~0.30 ppb = 3e-10
    mu_rel_unc = 3.0e-10
    mu_abs_unc = mu_measured * mu_rel_unc
    tension_sigma = abs(mu_predicted - mu_measured) / mu_abs_unc

    possible_resolutions = [
        "Higher-order radiative corrections to the bridge coefficient "
        "(e.g., alpha^3 or alpha^4 terms beyond 3*alpha^2)",
        "The bridge formula phi^2/2*(1+3*alpha^2) may need NLO QCD "
        "or electroweak corrections",
        "The mu-structure k = sqrt(phi)*(1-alpha) may receive threshold "
        "corrections from the KK tower",
        "Both formulas may be leading-order approximations to a single "
        "exact formula not yet discovered",
    ]

    return {
        "mu_measured": mu_measured,
        "mu_predicted_from_bridge": mu_predicted,
        "tension_ppm": tension_ppm,
        "tension_sigma": tension_sigma,
        "interpretation": (
            "The bridge and mu-structure formulas predict slightly different "
            "values of mu.  The {:.1f}-sigma tension ({:.2f} ppm) is too large "
            "to be measurement error and indicates that at least one formula "
            "needs higher-order corrections."
        ).format(tension_sigma, abs(tension_ppm)),
        "possible_resolutions": possible_resolutions,
        "assessment": (
            "Mu tension is {:.1f} sigma ({:.2f} ppm).  This is a genuine "
            "discrepancy that must be resolved for the framework to be "
            "self-consistent.  The most likely resolution is higher-order "
            "corrections to the bridge coefficient."
        ).format(tension_sigma, abs(tension_ppm)),
    }


# ---------------------------------------------------------------------------
# 7. Summary
# ---------------------------------------------------------------------------

def summarize_second_predictions(constants=None):
    """Main entry point: compute all second predictions and summarize.

    Returns
    -------
    dict
        time_variation, current_bounds, m6_landscape, correlated_variations,
        mu_tension, strongest_prediction, honest_assessment.
    """
    tv = compute_time_variation_coefficients(constants)
    cb = compute_current_bounds(constants)
    ml = scan_a0_m6_landscape(constants)
    cv = analyze_correlated_variations(dalpha_alpha=1e-17, constants=constants)
    mt = analyze_mu_tension(constants)

    return {
        "time_variation": tv,
        "current_bounds": cb,
        "m6_landscape": ml,
        "correlated_variations": cv,
        "mu_tension": mt,
        "strongest_prediction": (
            "M_6 in [3-5] TeV from the survival window is the most testable "
            "second prediction.  It is at the LHC energy frontier and could "
            "be probed by graviton emission searches at the HL-LHC or a "
            "future 100 TeV collider."
        ),
        "honest_assessment": (
            "The Alpha Ladder currently has ONE confirmed prediction (G to "
            "0.33 ppm) and ZERO confirmed second predictions.  The time "
            "variation signature (A~24) is unique but ~170x below current "
            "sensitivity.  The M_6 prediction is testable at the LHC but "
            "requires specific search strategies for n=2 extra dimensions.  "
            "The mu tension (~{:.0f} sigma) is a genuine internal inconsistency "
            "that needs resolution.  The framework is at the stage where it "
            "makes falsifiable predictions but has not yet been independently "
            "tested."
        ).format(mt["tension_sigma"]),
    }
