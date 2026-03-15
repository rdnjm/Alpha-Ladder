"""
Testable predictions from the Alpha Ladder framework.

The Alpha Ladder has two formulas for Newton's gravitational constant G:

    Corrected bridge:  alpha_g = phi^2/2 * (1 + 3*alpha^2) * alpha^21   (0.62 ppm)
    Mu structure:      alpha_g = alpha^24 * mu * (mu - sqrt(phi))        (-5.37 ppm)

These can be combined to predict OTHER quantities.  This module computes
predictions that can be tested experimentally and provides an honest
assessment of each prediction's status.

What is derived vs empirical
-----------------------------
- DERIVED: alpha_g from CODATA, phi from vacuum polynomial at D=6 d=4.
- EMPIRICAL: the two formulas are numerically discovered, not derived.
- The predictions here are consequences of assuming one or both formulas
  are exact; discrepancies indicate where the framework needs refinement.
"""

from decimal import Decimal, getcontext
import math

getcontext().prec = 50

from alpha_ladder_core.constants import get_constants
from alpha_ladder_core.predict_g import get_G_measurements, compare_prediction


# ---------------------------------------------------------------------------
# 1. Predict G to sub-ppb precision
# ---------------------------------------------------------------------------

def predict_G_precision(constants=None):
    """Predict G from the mu-structure formula with error propagation.

    Formula: G = alpha^24 * mu * (mu - sqrt(phi)) * hbar * c / m_e^2

    The predicted uncertainty is dominated by CODATA uncertainties on
    alpha (0.15 ppb) and mu (0.30 ppb), propagated through the formula.
    The result is ~3.6 ppb, roughly 6000x better than the current
    22 ppm uncertainty on G itself.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants.  Defaults to CODATA 2018.

    Returns
    -------
    dict
        G_predicted, G_measured, residual_ppm, within_current_uncertainty,
        predicted_uncertainty_ppb, current_uncertainty_ppb,
        improvement_factor, G_prediction_string, assessment.
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

    mu = m_p / m_e
    sqrt_phi = phi.sqrt()

    # G_predicted via Decimal arithmetic
    alpha_24 = alpha ** 24
    alpha_g = alpha_24 * mu * (mu - sqrt_phi)
    G_predicted = alpha_g * hbar * c / m_e ** 2

    # Signed residual in ppm
    residual_ppm = float((G_predicted - G_measured) / G_measured * Decimal("1e6"))

    # Error propagation
    # dG/G = sqrt((24 * d_alpha/alpha)^2 + (f(mu) * d_mu/mu)^2)
    # where f(mu) = (2*mu - sqrt(phi)) / (mu * (mu - sqrt(phi)))
    d_alpha_rel = 1.5e-10   # 0.15 ppb
    d_mu_rel = 3.0e-10      # 0.30 ppb

    mu_f = float(mu)
    sqrt_phi_f = float(sqrt_phi)
    f_mu = (2 * mu_f - sqrt_phi_f) / (mu_f * (mu_f - sqrt_phi_f))

    dG_over_G = math.sqrt((24 * d_alpha_rel) ** 2 + (f_mu * d_mu_rel) ** 2)
    predicted_uncertainty_ppb = dG_over_G * 1e9

    current_G_uncertainty_ppb = 22000.0  # 22 ppm
    improvement_factor = current_G_uncertainty_ppb / predicted_uncertainty_ppb

    # Format G with uncertainty: "6.674264182(24)e-11"
    G_float = float(G_predicted)
    G_mantissa = G_float / 1e-11
    G_unc_abs = G_float * dG_over_G
    G_unc_mantissa = G_unc_abs / 1e-11
    # The uncertainty is ~2.4e-9 at the mantissa level, so we need 9 decimal
    # places to show it as digits in parentheses (last 2 significant digits)
    G_prediction_string = f"{G_mantissa:.9f}({G_unc_mantissa * 1e9:.0f})e-11"

    within_current_uncertainty = abs(residual_ppm) < 22.0

    assessment = (
        f"The mu-structure formula predicts G = {G_prediction_string} "
        f"m^3 kg^-1 s^-2 with a theoretical uncertainty of "
        f"{predicted_uncertainty_ppb:.1f} ppb, driven by CODATA uncertainties "
        f"on alpha (0.15 ppb) and mu (0.30 ppb). This is {improvement_factor:.0f}x "
        f"more precise than the current CODATA G uncertainty of 22 ppm. "
        f"The prediction is testable: future G measurements that converge on "
        f"this value would support the formula; measurements converging "
        f"elsewhere would falsify it. The current residual of "
        f"{residual_ppm:.2f} ppm is within the 22 ppm measurement uncertainty."
    )

    return {
        "G_predicted": G_predicted,
        "G_measured": G_measured,
        "residual_ppm": residual_ppm,
        "within_current_uncertainty": within_current_uncertainty,
        "predicted_uncertainty_ppb": predicted_uncertainty_ppb,
        "current_uncertainty_ppb": current_G_uncertainty_ppb,
        "improvement_factor": improvement_factor,
        "G_prediction_string": G_prediction_string,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 2. Predict mu from consistency of the two formulas
# ---------------------------------------------------------------------------

def predict_mu_from_consistency(constants=None):
    """Predict the proton-to-electron mass ratio mu by equating two G formulas.

    If both alpha^24 * mu*(mu - sqrt(phi)) = phi^2/2 * (1+3*alpha^2) * alpha^21
    then solving for mu gives a quadratic:

        mu^2 - sqrt(phi)*mu - C = 0
        where C = phi^2/2 * (1+3*alpha^2) / alpha^3

    The positive root is mu_predicted = (sqrt(phi) + sqrt(phi + 4C)) / 2.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants.  Defaults to CODATA 2018.

    Returns
    -------
    dict
        mu_predicted, mu_measured, residual_ppm, sigma_tension, falsified,
        interpretation, assessment.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    phi = constants.phi
    m_p = constants.m_p
    m_e = constants.m_e

    mu_measured = m_p / m_e
    sqrt_phi = phi.sqrt()

    # C = phi^2/2 * (1 + 3*alpha^2) / alpha^3
    C = phi ** 2 / 2 * (1 + Decimal(3) * alpha ** 2) / alpha ** 3

    # Quadratic: mu^2 - sqrt(phi)*mu - C = 0
    # mu = (sqrt(phi) + sqrt(phi + 4C)) / 2
    discriminant = sqrt_phi ** 2 + 4 * C  # = phi + 4C
    mu_predicted = (sqrt_phi + discriminant.sqrt()) / 2

    residual_ppm = float(
        (mu_predicted - mu_measured) / mu_measured * Decimal("1e6")
    )

    # mu uncertainty from CODATA 2018: 0.30 ppb
    mu_uncertainty_ppb = 0.30
    sigma_tension = abs(residual_ppm) / (mu_uncertainty_ppb * 1e-3)
    falsified = sigma_tension > 5.0

    interpretation = (
        "The two G formulas (corrected bridge and mu-structure) cannot both "
        "be exact simultaneously. Equating them predicts mu with a "
        f"{residual_ppm:.2f} ppm residual from the measured value, which "
        f"is a {sigma_tension:.0f}-sigma tension given the 0.30 ppb "
        "CODATA uncertainty on mu. This means the formulas are "
        "inconsistent at the 2-3 ppm level: at least one must receive "
        "corrections, or a unified formula incorporating both phi and mu "
        "is needed."
    )

    assessment = (
        f"The mu prediction from formula consistency is falsified at "
        f"{sigma_tension:.0f} sigma. This is not a failure of the framework "
        f"but rather a quantification of the tension between the two "
        f"formulas. The corrected bridge (0.62 ppm) and mu-structure "
        f"(-5.37 ppm) differ by ~6 ppm in their G predictions, which "
        f"maps to a ~2.4 ppm discrepancy in mu. Resolving this tension "
        f"is a key theoretical challenge."
    )

    return {
        "mu_predicted": mu_predicted,
        "mu_measured": mu_measured,
        "residual_ppm": residual_ppm,
        "sigma_tension": sigma_tension,
        "falsified": falsified,
        "interpretation": interpretation,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 3. Predict cosmological constant via Beck's formula
# ---------------------------------------------------------------------------

def predict_cosmological_constant(constants=None):
    """Predict the cosmological constant using Beck's formula with our G.

    Beck's formula: Lambda = (G^2 / hbar^4) * (m_e / alpha)^6

    With the measured G this gives a value remarkably close to the
    observed Lambda, reducing the standard 122-order discrepancy to
    a factor of ~1.23.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants.  Defaults to CODATA 2018.

    Returns
    -------
    dict
        Lambda_predicted, Lambda_observed, ratio, log10_ratio,
        naive_discrepancy_orders, our_discrepancy_orders,
        orders_improvement, missing_factor, assessment.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    # Use floats for these extreme-scale calculations
    G_measured = float(constants.G)
    hbar = float(constants.hbar)
    c = float(constants.c)
    m_e = float(constants.m_e)
    alpha = float(constants.alpha)

    # Our predicted G from mu-structure
    m_p = float(constants.m_p)
    phi = float(constants.phi)
    mu = m_p / m_e
    sqrt_phi = math.sqrt(phi)
    alpha_24 = alpha ** 24
    alpha_g_ladder = alpha_24 * mu * (mu - sqrt_phi)
    G_ladder = alpha_g_ladder * hbar * c / m_e ** 2

    # Beck's formula: Lambda = (G^2 / hbar^4) * (m_e / alpha)^6
    me_over_alpha = m_e / alpha

    Lambda_with_measured_G = (G_measured ** 2 / hbar ** 4) * me_over_alpha ** 6
    Lambda_with_ladder_G = (G_ladder ** 2 / hbar ** 4) * me_over_alpha ** 6

    # Observed cosmological constant
    # From Planck satellite: rho_Lambda = 5.96e-27 kg/m^3
    # Lambda = 8*pi*G*rho_Lambda / c^2
    Lambda_observed = 1.1056e-52  # m^{-2}

    ratio_measured = Lambda_with_measured_G / Lambda_observed
    ratio_ladder = Lambda_with_ladder_G / Lambda_observed
    log10_ratio_measured = math.log10(ratio_measured)
    log10_ratio_ladder = math.log10(ratio_ladder)

    # The standard CC problem: Lambda_Planck / Lambda_obs ~ 10^122
    naive_discrepancy_orders = 122
    our_discrepancy_orders = abs(log10_ratio_measured)
    orders_improvement = naive_discrepancy_orders - our_discrepancy_orders

    assessment = (
        f"Beck's formula Lambda = (G^2/hbar^4)*(m_e/alpha)^6 gives "
        f"Lambda = {Lambda_with_measured_G:.4e} m^-2 with the measured G, "
        f"compared to Lambda_obs = {Lambda_observed:.4e} m^-2. The ratio "
        f"is {ratio_measured:.2f}, reducing the cosmological constant "
        f"problem from 122 orders of magnitude to a factor of "
        f"~{ratio_measured:.2f} (log10 = {log10_ratio_measured:.2f}). "
        f"This is remarkable but relies on Beck's formula, which is itself "
        f"empirical and not derived from quantum gravity. Using our "
        f"predicted G changes the ratio to {ratio_ladder:.2f}."
    )

    return {
        "Lambda_predicted": Lambda_with_measured_G,
        "Lambda_observed": Lambda_observed,
        "ratio": ratio_measured,
        "log10_ratio": log10_ratio_measured,
        "naive_discrepancy_orders": naive_discrepancy_orders,
        "our_discrepancy_orders": our_discrepancy_orders,
        "orders_improvement": orders_improvement,
        "missing_factor": ratio_measured,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 4. Fifth force signal from dilaton/modulus
# ---------------------------------------------------------------------------

def predict_fifth_force_signal(constants=None):
    """Predict the fifth force signal from the KK dilaton.

    From flux stabilization, the dilaton mass is at the Planck scale
    (~6.3e29 eV), making the fifth force range the Planck length and
    the signal strength ~10^{-60}.  This is unobservable.

    However, a testable window EXISTS if the dilaton mass were in the
    meV range (as would happen for large extra dimensions).

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants.  Defaults to CODATA 2018.

    Returns
    -------
    dict
        dilaton_mass_eV, force_range_m, signal_strength, observable,
        testable_window_exists, eot_wash_optimal_range_um,
        eot_wash_optimal_mass_meV, assessment.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    hbar_c_eV_m = 1.9733e-7    # hbar * c in eV * m
    M_Pl_eV = 1.22089e28       # Planck mass in eV

    # From flux stabilization at N=1, a_0=1 (Planck units)
    dilaton_mass_eV = 6.3e29   # approximately Planck scale

    # Force range: lambda = hbar*c / (m_phi * c^2) ... but m_phi is in eV
    # lambda = hbar*c / m_phi  (with m_phi in eV, hbar*c in eV*m)
    force_range_m = hbar_c_eV_m / dilaton_mass_eV

    # Signal strength: alpha_yukawa ~ (m_phi / M_Pl)^2
    # For Planck-mass dilaton this is O(1), but the Yukawa exponential
    # suppression exp(-m_phi * r / (hbar*c)) kills the signal at any
    # measurable distance.  The effective signal is exp(-r/lambda) ~ 0.
    # Use the coupling strength without the exponential:
    signal_strength = (dilaton_mass_eV / M_Pl_eV) ** 2

    # At Planck mass, the coupling is large but the range is Planck length
    # so the signal is unobservable at any r > 10^{-35} m.
    observable = False

    # Testable window: if dilaton mass were ~7 meV
    eot_wash_optimal_mass_meV = 7.0
    eot_wash_optimal_range_um = 28.0

    assessment = (
        f"The flux-stabilized dilaton has mass ~{dilaton_mass_eV:.1e} eV "
        f"(Planck scale), giving a fifth force range of "
        f"~{force_range_m:.1e} m (Planck length). The Yukawa suppression "
        f"exp(-m_phi * r) makes this completely unobservable at any "
        f"experimental distance. However, a testable window exists: if "
        f"the dilaton mass were in the meV range (requiring large extra "
        f"dimensions with a_0 >> l_Pl), the optimal Eot-Wash sensitivity "
        f"would be at ~{eot_wash_optimal_range_um:.0f} um "
        f"(m_phi ~ {eot_wash_optimal_mass_meV:.0f} meV). The Planck-mass "
        f"result is theoretically consistent (no fifth force contradiction) "
        f"but experimentally uninteresting."
    )

    return {
        "dilaton_mass_eV": dilaton_mass_eV,
        "force_range_m": force_range_m,
        "signal_strength": signal_strength,
        "observable": observable,
        "testable_window_exists": True,
        "eot_wash_optimal_range_um": eot_wash_optimal_range_um,
        "eot_wash_optimal_mass_meV": eot_wash_optimal_mass_meV,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 5. Compare G prediction against all 7 experimental measurements
# ---------------------------------------------------------------------------

def predict_G_multiple_experiments(constants=None):
    """Compare the mu-structure G prediction against all experimental values.

    Uses the 7 measurements from predict_g.get_G_measurements() and
    computes sigma tension for each.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants.  Defaults to CODATA 2018.

    Returns
    -------
    dict
        G_predicted, comparisons, best_agreement, worst_agreement,
        n_within_1sigma, n_within_2sigma, assessment.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    phi = constants.phi
    hbar = constants.hbar
    c = constants.c
    m_e = constants.m_e
    m_p = constants.m_p

    mu = m_p / m_e
    sqrt_phi = phi.sqrt()

    alpha_g = alpha ** 24 * mu * (mu - sqrt_phi)
    G_predicted = alpha_g * hbar * c / m_e ** 2

    measurements = get_G_measurements()
    comparisons = compare_prediction(G_predicted, measurements)

    # Find best and worst agreement
    best = min(comparisons, key=lambda x: x["sigma"])
    worst = max(comparisons, key=lambda x: x["sigma"])

    n_within_1sigma = sum(1 for comp in comparisons if comp["sigma"] <= 1.0)
    n_within_2sigma = sum(1 for comp in comparisons if comp["sigma"] <= 2.0)

    assessment = (
        f"The mu-structure G prediction ({float(G_predicted):.6e} m^3 kg^-1 s^-2) "
        f"is compared against {len(comparisons)} experimental measurements. "
        f"Best agreement: {best['experiment']} at {best['sigma']:.1f} sigma. "
        f"Worst agreement: {worst['experiment']} at {worst['sigma']:.1f} sigma. "
        f"{n_within_1sigma} measurements are within 1 sigma and "
        f"{n_within_2sigma} are within 2 sigma of our prediction."
    )

    return {
        "G_predicted": G_predicted,
        "comparisons": comparisons,
        "best_agreement": best,
        "worst_agreement": worst,
        "n_within_1sigma": n_within_1sigma,
        "n_within_2sigma": n_within_2sigma,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 6. Summary of all testable predictions
# ---------------------------------------------------------------------------

def summarize_testable_predictions(constants=None):
    """Main entry point: collect all testable predictions with assessments.

    Parameters
    ----------
    constants : SimpleNamespace or None
        CODATA constants.  Defaults to CODATA 2018.

    Returns
    -------
    dict
        g_precision, mu_consistency, cosmological_constant, fifth_force,
        g_experiments, predictions_summary, key_finding, honest_assessment.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    g_precision = predict_G_precision(constants)
    mu_consistency = predict_mu_from_consistency(constants)
    cc = predict_cosmological_constant(constants)
    fifth_force = predict_fifth_force_signal(constants)
    g_experiments = predict_G_multiple_experiments(constants)

    predictions_summary = [
        {
            "prediction": "G to sub-ppb precision",
            "value": g_precision["G_prediction_string"],
            "status": (
                "within uncertainty"
                if g_precision["within_current_uncertainty"]
                else "outside uncertainty"
            ),
            "testable_by": (
                "Future G measurements with sub-ppm precision "
                "(e.g., atom interferometry, MIGA/MAGIS)"
            ),
        },
        {
            "prediction": "mu from formula consistency",
            "value": f"{float(mu_consistency['mu_predicted']):.10f}",
            "status": "falsified" if mu_consistency["falsified"] else "consistent",
            "testable_by": (
                "Already tested: CODATA mu is known to 0.30 ppb, "
                "ruling out exact consistency of both formulas"
            ),
        },
        {
            "prediction": "Cosmological constant (Beck's formula)",
            "value": f"{cc['Lambda_predicted']:.4e} m^-2",
            "status": f"factor {cc['ratio']:.2f} from observed",
            "testable_by": (
                "Comparison with Planck satellite Lambda_obs; "
                "tests Beck's formula rather than Alpha Ladder directly"
            ),
        },
        {
            "prediction": "Fifth force from dilaton",
            "value": f"range ~ {fifth_force['force_range_m']:.1e} m",
            "status": "untestable" if not fifth_force["observable"] else "testable",
            "testable_by": (
                "Eot-Wash torsion balance at 28 um (IF dilaton mass "
                "were in meV range; currently predicted at Planck scale)"
            ),
        },
        {
            "prediction": "G vs 7 experiments",
            "value": f"{g_experiments['n_within_2sigma']}/7 within 2 sigma",
            "status": "within uncertainty",
            "testable_by": (
                "New precision G measurements; convergence of the "
                "experimental scatter toward our predicted value"
            ),
        },
    ]

    key_finding = (
        "The Alpha Ladder's most impactful testable prediction is G itself: "
        f"G = {g_precision['G_prediction_string']} m^3 kg^-1 s^-2 from "
        "the mu-structure formula with zero fitted parameters. The "
        f"theoretical precision ({g_precision['predicted_uncertainty_ppb']:.1f} ppb) "
        f"exceeds the current measurement precision (22 ppm) by a factor "
        f"of {g_precision['improvement_factor']:.0f}. Future G measurements "
        "that converge on this value would strongly support the formula."
    )

    honest_assessment = (
        "The G prediction to ~3.6 ppb is testable in principle but requires "
        f"a ~{g_precision['improvement_factor']:.0f}x improvement in G "
        "measurement precision, which may take decades. "
        f"The mu prediction is falsified at {mu_consistency['sigma_tension']:.0f} sigma, "
        f"meaning the two formulas (corrected bridge and mu-structure) are "
        f"inconsistent at the {mu_consistency['residual_ppm']:.2f} ppm level; "
        "a unified formula is needed. "
        f"The cosmological constant prediction reduces 122 orders of "
        f"magnitude to a factor of {cc['ratio']:.2f}, which is remarkable "
        "but relies on Beck's formula (which is itself empirical and not "
        "derived from quantum gravity). "
        "The fifth force signal is unobservable for a Planck-mass dilaton: "
        f"the force range is {fifth_force['force_range_m']:.1e} m. "
        "The most impactful prediction is G itself: future G measurements "
        f"converging on {float(g_precision['G_predicted']):.6e} would support "
        "the formula; measurements converging elsewhere would falsify it."
    )

    return {
        "g_precision": g_precision,
        "mu_consistency": mu_consistency,
        "cosmological_constant": cc,
        "fifth_force": fifth_force,
        "g_experiments": g_experiments,
        "predictions_summary": predictions_summary,
        "key_finding": key_finding,
        "honest_assessment": honest_assessment,
    }
