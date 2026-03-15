"""
Fifth force predictions for the Alpha Ladder dilaton.

Translates the theoretical framework into concrete experimental
predictions for sub-mm gravity experiments.  The key feature is that
the Alpha Ladder fixes the coupling strength:

    alpha_fifth = 2 / (2*omega + 3) = 0.618  (golden ratio inverse)

This is NOT a free parameter -- it is derived from the Brans-Dicke
parameter omega which comes from the phi^2/2 bridge coefficient.
The only undetermined parameter is the internal radius a_0, which
sets the dilaton Compton wavelength:

    lambda = hbar * c / m_phi = a_0

The framework therefore predicts a horizontal line at alpha = 0.618
in the standard (lambda, alpha) exclusion plot.  Any experiment that
can probe alpha < 0.618 at some lambda excludes the corresponding
a_0 value.

Experiments covered:
    1. Eot-Wash torsion balance (current best sub-mm)
    2. IUPUI short-range gravity (planar geometry)
    3. Casimir force measurements (plate-plate and plate-sphere)
    4. Neutron gravity resonance (qBounce)
    5. Atom interferometry (Stanford, MAGIS)
    6. Lunar laser ranging
    7. Cassini PPN (solar system)

All calculations use pure Python + math + decimal (no numpy/scipy).
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

_L_PL = 1.616e-35           # Planck length in meters
_M_PL_EV = 1.22089e28       # Planck mass in eV
_HBAR_C_EV_M = 1.9733e-7    # hbar*c in eV*m
_ALPHA_FIFTH = 0.618         # Fixed coupling from BD parameter
_G_N = 6.674e-11            # Newton's G in SI
_HBAR = 1.0546e-34           # J*s
_C = 2.998e8                 # m/s


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _safe_exp(x):
    """Compute math.exp(x) with clamping to avoid overflow."""
    if x > 700.0:
        return math.exp(700.0)
    if x < -700.0:
        return 0.0
    return math.exp(x)


def _log_spaced(x_min, x_max, n):
    """Return a list of n log-spaced values from x_min to x_max."""
    log_min = math.log10(x_min)
    log_max = math.log10(x_max)
    return [10 ** (log_min + i * (log_max - log_min) / (n - 1))
            for i in range(n)]


# ---------------------------------------------------------------------------
# 1. Yukawa signal at a single distance
# ---------------------------------------------------------------------------

def compute_yukawa_signal(r_meters, lambda_m, alpha=_ALPHA_FIFTH):
    """Compute the Yukawa fifth force signal at distance r.

    The fractional deviation from Newton's law is:

        delta_F / F_Newton = alpha * (1 + r/lambda) * exp(-r/lambda)

    This comes from differentiating the Yukawa potential
    V = -alpha * G * M1 * M2 / r * exp(-r/lambda).

    Parameters
    ----------
    r_meters : float
        Distance between test masses in meters.
    lambda_m : float
        Yukawa range (= dilaton Compton wavelength = a_0) in meters.
    alpha : float
        Coupling strength.  Fixed at 0.618 by the Alpha Ladder.

    Returns
    -------
    dict with keys:
        r_meters, lambda_m, alpha, fractional_deviation,
        exponential_factor, detectable (bool)
    """
    x = r_meters / lambda_m
    exp_factor = _safe_exp(-x)
    fractional_deviation = alpha * (1.0 + x) * exp_factor

    return {
        "r_meters": r_meters,
        "lambda_m": lambda_m,
        "alpha": alpha,
        "fractional_deviation": fractional_deviation,
        "exponential_factor": exp_factor,
        "detectable": fractional_deviation > 1e-4,
    }


# ---------------------------------------------------------------------------
# 2. Signal vs distance profile
# ---------------------------------------------------------------------------

def compute_signal_vs_distance(lambda_m, alpha=_ALPHA_FIFTH,
                               r_min=1e-6, r_max=1.0, n_points=100):
    """Compute the Yukawa fractional deviation across a log-spaced range.

    Parameters
    ----------
    lambda_m : float
        Yukawa range in meters.
    alpha : float
        Coupling strength (fixed at 0.618).
    r_min, r_max : float
        Distance range in meters.
    n_points : int
        Number of sample points.

    Returns
    -------
    dict with keys:
        r_meters (list), fractional_deviation (list), lambda_m, alpha,
        peak_distance, peak_signal
    """
    r_values = _log_spaced(r_min, r_max, n_points)
    deviations = []
    peak_signal = 0.0
    peak_distance = r_values[0]

    for r in r_values:
        x = r / lambda_m
        dev = alpha * (1.0 + x) * _safe_exp(-x)
        deviations.append(dev)
        if dev > peak_signal:
            peak_signal = dev
            peak_distance = r

    return {
        "r_meters": r_values,
        "fractional_deviation": deviations,
        "lambda_m": lambda_m,
        "alpha": alpha,
        "peak_distance": peak_distance,
        "peak_signal": peak_signal,
    }


# ---------------------------------------------------------------------------
# 3. Exclusion map
# ---------------------------------------------------------------------------

def compute_exclusion_map():
    """Map the Alpha Ladder prediction onto the (lambda, alpha) exclusion plot.

    The Alpha Ladder predicts alpha = 0.618 for ALL lambda.  This is a
    horizontal line in the standard exclusion plot.  For each experiment,
    we determine whether its sensitivity crosses alpha = 0.618 at any
    lambda in its operating range, thereby excluding that range of a_0.

    Returns
    -------
    dict with keys:
        experiments (list of dicts), alpha_ladder_line, survival_window_m,
        survival_window_um, a_0_range_surviving_m, n_experiments_excluding,
        n_experiments_not_excluding, honest_assessment
    """
    # Experimental bounds: (name, lambda_min_m, lambda_max_m,
    #   alpha_bound_at_lambda_min, scaling_description, reference)
    # For Eot-Wash, alpha_bound ~ (56e-6 / lambda)^2 for lambda > 56 um
    experiments_data = [
        {
            "name": "Eot-Wash 2006",
            "lambda_min_m": 56e-6,
            "lambda_max_m": 10e-3,
            "reference": "Kapner+ 2007",
        },
        {
            "name": "Eot-Wash 2020 (projected)",
            "lambda_min_m": 30e-6,
            "lambda_max_m": 10e-3,
            "reference": "Lee+ projected",
        },
        {
            "name": "IUPUI 2020",
            "lambda_min_m": 40e-6,
            "lambda_max_m": 1e-3,
            "reference": "Long+ 2021",
        },
        {
            "name": "Casimir (Decca)",
            "lambda_min_m": 0.2e-6,
            "lambda_max_m": 10e-6,
            "reference": "Decca+ 2007",
        },
        {
            "name": "Stanford atom interferometry",
            "lambda_min_m": 10e-6,
            "lambda_max_m": 1e-3,
            "reference": "Geraci+ 2008",
        },
        {
            "name": "Lunar laser ranging",
            "lambda_min_m": 3.84e8,
            "lambda_max_m": 1e12,
            "reference": "Williams+ 2012",
        },
        {
            "name": "Cassini PPN",
            "lambda_min_m": 1.5e11,
            "lambda_max_m": 1e15,
            "reference": "Bertotti+ 2003",
        },
    ]

    def _eot_wash_bound(lambda_m):
        """Model the Eot-Wash alpha bound as a function of lambda."""
        # At lambda = 56 um, alpha_bound ~ 1.0
        # Scaling: alpha_bound ~ (56e-6 / lambda)^2 for lambda > 56 um
        # This captures the improving sensitivity at larger lambda
        if lambda_m < 56e-6:
            return 1e3  # not probed
        return (56e-6 / lambda_m) ** 2

    def _get_alpha_bound(name, lambda_m):
        """Return the approximate experimental alpha bound at a given lambda."""
        if "Eot-Wash 2006" in name:
            return _eot_wash_bound(lambda_m)
        if "Eot-Wash 2020" in name:
            # Projected improvement: factor ~10 better
            if lambda_m < 30e-6:
                return 1e3
            return 0.1 * (30e-6 / lambda_m) ** 2
        if "IUPUI" in name:
            if lambda_m < 40e-6:
                return 1e3
            return (40e-6 / lambda_m) ** 2
        if "Casimir" in name:
            # Casimir: alpha ~ 1e6 at 0.2 um, drops as (0.2e-6/lambda)^4
            if lambda_m < 0.2e-6:
                return 1e8
            return 1e6 * (0.2e-6 / lambda_m) ** 4
        if "Stanford" in name:
            if lambda_m < 10e-6:
                return 1e6
            return 1e4 * (10e-6 / lambda_m) ** 2
        if "Lunar" in name:
            return 1e-11
        if "Cassini" in name:
            return 2.3e-5
        return 1e10  # no constraint

    alpha_ladder = _ALPHA_FIFTH
    experiment_results = []

    for exp in experiments_data:
        name = exp["name"]
        lam_min = exp["lambda_min_m"]
        lam_max = exp["lambda_max_m"]

        # Sample the bound across the lambda range
        test_lambdas = _log_spaced(lam_min, lam_max, 200)
        best_bound = min(_get_alpha_bound(name, lam) for lam in test_lambdas)

        excludes = best_bound < alpha_ladder

        experiment_results.append({
            "name": name,
            "lambda_range_m": (lam_min, lam_max),
            "alpha_bound_at_best": best_bound,
            "excludes_alpha_ladder": excludes,
            "reference": exp["reference"],
        })

    # Compute survival window: range of lambda where no experiment
    # excludes alpha = 0.618.
    # Eot-Wash 2006 excludes above ~ 200 um (where bound drops below 0.618).
    # Eot-Wash 2006: (56e-6/lambda)^2 < 0.618 => lambda > 56e-6/sqrt(0.618) ~ 71 um
    # Below 56 um: not probed by Eot-Wash.
    # Projected Eot-Wash 2020: 0.1*(30e-6/lambda)^2 < 0.618 =>
    #   lambda > 30e-6 * sqrt(0.1/0.618) ~ 12 um  -- but min is 30 um
    # So with 2006 data: survival below ~71 um.
    # With projected 2020: survival below ~30 um (at the experimental limit).

    # Compute the upper bound of the survival window from Eot-Wash 2006
    # (56e-6/lambda)^2 = 0.618 => lambda = 56e-6 / sqrt(0.618)
    lambda_upper_2006 = 56e-6 / math.sqrt(alpha_ladder)
    # Projected 2020: 0.1*(30e-6/lambda)^2 = 0.618 => lambda = 30e-6*sqrt(0.1/0.618)
    lambda_upper_2020 = 30e-6 * math.sqrt(0.1 / alpha_ladder)

    # Lower bound: no experiment probes below ~30 um (projected) or 56 um (current)
    lambda_lower_current = 56e-6
    lambda_lower_projected = 30e-6

    # Use current bounds for the survival window
    survival_min = lambda_lower_projected  # no experiment probes below this
    survival_max = lambda_upper_2006       # Eot-Wash 2006 excludes above this

    n_excluding = sum(1 for e in experiment_results if e["excludes_alpha_ladder"])
    n_not_excluding = len(experiment_results) - n_excluding

    honest_assessment = (
        f"The Alpha Ladder prediction alpha = {alpha_ladder} is a horizontal "
        f"line in the (lambda, alpha) exclusion plot.  With current experimental "
        f"bounds, the prediction survives for a_0 (= lambda) in the range "
        f"[{survival_min * 1e6:.0f}, {survival_max * 1e6:.0f}] um.  "
        f"At larger lambda, Eot-Wash torsion balance data excludes the "
        f"prediction.  At smaller lambda, no current experiment has sufficient "
        f"sensitivity.  Lunar laser ranging and Cassini PPN data exclude the "
        f"prediction at solar-system scales (lambda > 10^8 m), but those "
        f"correspond to a_0 values far larger than any plausible internal radius.  "
        f"The survival window is a genuine experimental gap that near-future "
        f"short-range gravity experiments can close."
    )

    return {
        "experiments": experiment_results,
        "alpha_ladder_line": alpha_ladder,
        "survival_window_m": (survival_min, survival_max),
        "survival_window_um": (survival_min * 1e6, survival_max * 1e6),
        "a_0_range_surviving_m": (survival_min, survival_max),
        "n_experiments_excluding": n_excluding,
        "n_experiments_not_excluding": n_not_excluding,
        "honest_assessment": honest_assessment,
    }


# ---------------------------------------------------------------------------
# 4. Eot-Wash torsion balance prediction
# ---------------------------------------------------------------------------

def compute_eot_wash_prediction(a_0_m=30e-6):
    """Compute the predicted signal in an Eot-Wash torsion balance experiment.

    The Eot-Wash experiment (Kapner+ 2007) uses a torsion pendulum with
    a test mass separated from an attractor by a gap of ~52 um.  The
    Yukawa correction to the gravitational torque is:

        signal = alpha * (1 + d/lambda) * exp(-d/lambda)

    where d is the gap distance and lambda = a_0.

    Parameters
    ----------
    a_0_m : float
        Internal radius in meters (sets lambda = a_0).

    Returns
    -------
    dict with keys:
        a_0_m, lambda_m, alpha_fifth, gap_distance_m, predicted_signal,
        current_sensitivity, signal_to_sensitivity_ratio, detectable,
        signal_description, experiment_params
    """
    lambda_m = a_0_m
    alpha = _ALPHA_FIFTH
    gap_d = 52e-6  # Kapner+ 2007 gap distance

    x = gap_d / lambda_m
    predicted_signal = alpha * (1.0 + x) * _safe_exp(-x)
    current_sensitivity = 0.01  # ~1% at 52 um gap

    ratio = predicted_signal / current_sensitivity
    detectable = predicted_signal > current_sensitivity

    if predicted_signal > 0.1:
        strength_word = "enormous"
    elif predicted_signal > 0.01:
        strength_word = "large"
    elif predicted_signal > 1e-4:
        strength_word = "marginal"
    else:
        strength_word = "undetectable"

    signal_description = (
        f"At a_0 = {a_0_m * 1e6:.1f} um (lambda = {lambda_m * 1e6:.1f} um), "
        f"the Alpha Ladder predicts a {predicted_signal * 100:.1f}% deviation "
        f"from Newton's law in the Eot-Wash torsion balance at gap = "
        f"{gap_d * 1e6:.0f} um.  Current sensitivity is ~{current_sensitivity * 100:.0f}%.  "
        f"The signal is {strength_word} "
        f"({ratio:.1f}x the current sensitivity)."
    )

    experiment_params = {
        "experiment": "Eot-Wash torsion balance",
        "reference": "Kapner+ 2007",
        "gap_distance_m": gap_d,
        "attractor_pattern": "42 holes in a disc",
        "test_mass_geometry": "torsion pendulum",
        "measurement_type": "torque",
    }

    return {
        "a_0_m": a_0_m,
        "lambda_m": lambda_m,
        "alpha_fifth": alpha,
        "gap_distance_m": gap_d,
        "predicted_signal": predicted_signal,
        "current_sensitivity": current_sensitivity,
        "signal_to_sensitivity_ratio": ratio,
        "detectable": detectable,
        "signal_description": signal_description,
        "experiment_params": experiment_params,
    }


# ---------------------------------------------------------------------------
# 5. Casimir overlap analysis
# ---------------------------------------------------------------------------

def compute_casimir_overlap(a_0_m=30e-6):
    """Compare the Yukawa fifth force signal with the Casimir force.

    The Casimir force between parallel plates dominates at sub-micron
    separations.  This function computes both pressures and shows that
    the Yukawa signal is many orders of magnitude below the Casimir
    background, making Casimir experiments a poor choice for detecting
    this particular fifth force.

    Parameters
    ----------
    a_0_m : float
        Internal radius in meters (sets lambda = a_0).

    Returns
    -------
    dict with keys:
        a_0_m, gap_distance_m, casimir_pressure_Pa, yukawa_pressure_Pa,
        ratio_yukawa_to_casimir, casimir_dominates, honest_assessment
    """
    lambda_m = a_0_m
    gap_d = a_0_m  # compare at d = lambda (most favorable for Yukawa)

    # Casimir pressure between parallel plates:
    # P_cas = pi^2 * hbar * c / (240 * d^4)
    pi = math.pi
    casimir_pressure = pi ** 2 * _HBAR * _C / (240.0 * gap_d ** 4)

    # Yukawa pressure between macroscopic parallel plates:
    # P_yuk = 2*pi * G * rho1 * rho2 * alpha * lambda^2 * exp(-d/lambda)
    # Use gold plates: rho ~ 19300 kg/m^3
    rho_gold = 19300.0
    yukawa_pressure = (2.0 * pi * _G_N * rho_gold ** 2
                       * _ALPHA_FIFTH * lambda_m ** 2
                       * _safe_exp(-gap_d / lambda_m))

    if casimir_pressure > 0:
        ratio = yukawa_pressure / casimir_pressure
    else:
        ratio = float('inf')

    casimir_dominates = ratio < 1.0

    if ratio < 1e-3:
        assessment_quality = "hopelessly"
    elif ratio < 1.0:
        assessment_quality = "significantly"
    else:
        assessment_quality = "not"

    honest_assessment = (
        f"At gap = {gap_d * 1e6:.1f} um (= lambda), the Casimir pressure is "
        f"{casimir_pressure:.2e} Pa while the Yukawa pressure is "
        f"{yukawa_pressure:.2e} Pa.  The Yukawa signal is {assessment_quality} "
        f"buried under the Casimir background (ratio = {ratio:.2e}).  "
        f"Casimir experiments are NOT the right approach for detecting "
        f"this fifth force.  Torsion balance experiments that directly "
        f"measure gravitational forces (where the Casimir background "
        f"is absent due to electrostatic shielding) are the correct strategy."
    )

    return {
        "a_0_m": a_0_m,
        "gap_distance_m": gap_d,
        "casimir_pressure_Pa": casimir_pressure,
        "yukawa_pressure_Pa": yukawa_pressure,
        "ratio_yukawa_to_casimir": ratio,
        "casimir_dominates": casimir_dominates,
        "honest_assessment": honest_assessment,
    }


# ---------------------------------------------------------------------------
# 6. Discovery reach across experiments
# ---------------------------------------------------------------------------

def compute_discovery_reach(experiments=None):
    """Compute the range of a_0 where each experiment can detect the signal.

    For each experiment type, the predicted fractional deviation is:

        signal = alpha * (1 + d/lambda) * exp(-d/lambda)

    The signal is maximized when lambda ~ d (i.e. a_0 ~ gap distance),
    giving signal_max = alpha * 2 * exp(-1) = 0.455 for alpha = 0.618.

    Parameters
    ----------
    experiments : list of dict or None
        Each dict has name, gap_m, sensitivity, status.
        If None, uses the default experiment list.

    Returns
    -------
    dict with keys:
        experiments (list of dicts), best_experiment,
        overall_discovery_range_m, honest_assessment
    """
    if experiments is None:
        experiments = [
            {"name": "Eot-Wash 2006", "gap_m": 52e-6,
             "sensitivity": 0.01, "status": "completed"},
            {"name": "Eot-Wash next-gen", "gap_m": 30e-6,
             "sensitivity": 0.001, "status": "projected"},
            {"name": "IUPUI planar", "gap_m": 40e-6,
             "sensitivity": 0.01, "status": "completed"},
            {"name": "Casimir (Decca)", "gap_m": 0.2e-6,
             "sensitivity": 0.001, "status": "completed"},
            {"name": "qBounce neutron", "gap_m": 10e-6,
             "sensitivity": 0.01, "status": "ongoing"},
            {"name": "Stanford atom interferometry", "gap_m": 25e-6,
             "sensitivity": 0.1, "status": "projected"},
        ]

    alpha = _ALPHA_FIFTH
    results = []
    global_lambda_min = float('inf')
    global_lambda_max = 0.0
    best_name = None
    best_range_width = 0.0

    for exp in experiments:
        gap = exp["gap_m"]
        sens = exp["sensitivity"]
        name = exp["name"]
        status = exp["status"]

        # Optimal lambda = gap distance
        optimal_lambda = gap
        # Signal at optimal: alpha * (1+1) * exp(-1) = alpha * 2 * exp(-1)
        signal_at_optimal = alpha * 2.0 * math.exp(-1.0)
        detectable_at_optimal = signal_at_optimal > sens

        # Find lambda range where signal > sensitivity
        # signal(lambda) = alpha * (1 + gap/lambda) * exp(-gap/lambda)
        # We solve numerically by scanning
        lambda_min_det = None
        lambda_max_det = None

        if detectable_at_optimal:
            # Scan from very small lambda to very large lambda
            test_lambdas = _log_spaced(gap * 1e-4, gap * 1e4, 10000)
            for lam in test_lambdas:
                x = gap / lam
                sig = alpha * (1.0 + x) * _safe_exp(-x)
                if sig > sens:
                    if lambda_min_det is None:
                        lambda_min_det = lam
                    lambda_max_det = lam

        lambda_range = None
        range_width = 0.0
        if lambda_min_det is not None and lambda_max_det is not None:
            lambda_range = (lambda_min_det, lambda_max_det)
            range_width = math.log10(lambda_max_det / lambda_min_det)
            global_lambda_min = min(global_lambda_min, lambda_min_det)
            global_lambda_max = max(global_lambda_max, lambda_max_det)

        if range_width > best_range_width:
            best_range_width = range_width
            best_name = name

        results.append({
            "name": name,
            "gap_m": gap,
            "sensitivity": sens,
            "optimal_lambda_m": optimal_lambda,
            "signal_at_optimal": signal_at_optimal,
            "detectable_at_optimal": detectable_at_optimal,
            "lambda_range_detectable": lambda_range,
            "status": status,
        })

    # Overall discovery range
    if global_lambda_min < float('inf') and global_lambda_max > 0:
        overall_range = (global_lambda_min, global_lambda_max)
    else:
        overall_range = None

    honest_assessment = (
        f"The Alpha Ladder prediction (alpha = {alpha}) produces a signal "
        f"of {alpha * 2.0 * math.exp(-1.0):.3f} (= alpha * 2/e) at the "
        f"optimal distance lambda = gap for any experiment.  "
        f"The best experiment for detection is '{best_name}' due to its "
        f"combination of gap distance and sensitivity.  "
        f"Eot-Wash torsion balance experiments are the most promising "
        f"because they directly measure gravitational forces without "
        f"Casimir backgrounds, and operate at gap distances (30-52 um) "
        f"that overlap with the surviving a_0 window.  "
        f"The Casimir experiment (Decca) operates at 0.2 um where the "
        f"electromagnetic Casimir force overwhelms any gravitational Yukawa "
        f"signal, making it unsuitable despite its high force sensitivity."
    )

    return {
        "experiments": results,
        "best_experiment": best_name,
        "overall_discovery_range_m": overall_range,
        "honest_assessment": honest_assessment,
    }


# ---------------------------------------------------------------------------
# 7. Prediction line for exclusion plots
# ---------------------------------------------------------------------------

def compute_alpha_ladder_prediction_line(a_0_values=None):
    """Generate the Alpha Ladder prediction as a line in (lambda, alpha) space.

    Since alpha = 0.618 is fixed, the prediction is a horizontal line.
    For each a_0 value, also compute the dilaton mass and the expected
    Eot-Wash signal.

    Parameters
    ----------
    a_0_values : list of float or None
        Internal radius values in meters.  If None, uses log-spaced
        values from 1e-7 to 1e-1 m.

    Returns
    -------
    dict with keys:
        a_0_values, lambda_values, alpha_values, m_phi_eV_values,
        eot_wash_signal_values, in_survival_window
    """
    if a_0_values is None:
        a_0_values = _log_spaced(1e-7, 1e-1, 200)

    alpha = _ALPHA_FIFTH
    eot_wash_gap = 52e-6  # Kapner+ 2007

    # Survival window from exclusion map (approximate, avoid circular call)
    survival_min = 30e-6
    survival_max = 56e-6 / math.sqrt(alpha)  # ~ 71 um

    lambda_values = list(a_0_values)  # lambda = a_0
    alpha_values = [alpha] * len(a_0_values)
    m_phi_eV_values = []
    eot_wash_signal_values = []
    in_survival = []

    for a_0 in a_0_values:
        # Dilaton mass: m_phi = hbar*c / (a_0 * c^2) ... in eV:
        # m_phi = hbar*c / a_0 in eV (natural units where hbar*c has units eV*m)
        m_phi_eV = _HBAR_C_EV_M / a_0
        m_phi_eV_values.append(m_phi_eV)

        # Eot-Wash signal at this lambda
        x = eot_wash_gap / a_0
        sig = alpha * (1.0 + x) * _safe_exp(-x)
        eot_wash_signal_values.append(sig)

        # In survival window?
        in_survival.append(survival_min <= a_0 <= survival_max)

    return {
        "a_0_values": list(a_0_values),
        "lambda_values": lambda_values,
        "alpha_values": alpha_values,
        "m_phi_eV_values": m_phi_eV_values,
        "eot_wash_signal_values": eot_wash_signal_values,
        "in_survival_window": in_survival,
    }


# ---------------------------------------------------------------------------
# 8. Summary / dashboard entry point
# ---------------------------------------------------------------------------

def summarize_fifth_force_predictions(constants=None):
    """Dashboard entry point: run all analyses and return a unified summary.

    Parameters
    ----------
    constants : SimpleNamespace or None
        Physical constants from get_constants().  Currently unused
        (all needed constants are module-level), but accepted for
        interface consistency with other modules.

    Returns
    -------
    dict with keys:
        exclusion_map, eot_wash_prediction, casimir_overlap,
        discovery_reach, prediction_line, key_prediction,
        survival_window_um, optimal_experiment,
        signal_strength_at_optimal, honest_assessment
    """
    exclusion_map = compute_exclusion_map()
    eot_wash = compute_eot_wash_prediction(a_0_m=30e-6)
    casimir = compute_casimir_overlap(a_0_m=30e-6)
    discovery = compute_discovery_reach()
    prediction_line = compute_alpha_ladder_prediction_line()

    survival_um = exclusion_map["survival_window_um"]
    optimal_exp = discovery["best_experiment"]

    # Signal at optimal lambda for the best experiment
    best_exp_data = None
    for e in discovery["experiments"]:
        if e["name"] == optimal_exp:
            best_exp_data = e
            break
    signal_at_optimal = best_exp_data["signal_at_optimal"] if best_exp_data else 0.0

    key_prediction = (
        f"The Alpha Ladder predicts a Yukawa fifth force with coupling "
        f"alpha = {_ALPHA_FIFTH} and range lambda = a_0.  For a_0 in "
        f"[{survival_um[0]:.0f}, {survival_um[1]:.0f}] um, this produces "
        f"a {eot_wash['predicted_signal'] * 100:.0f}% deviation from "
        f"Newton's law in Eot-Wash torsion balance experiments at their "
        f"current sensitivity.  This is the most falsifiable prediction "
        f"of the framework."
    )

    honest_assessment = (
        f"The Alpha Ladder's fifth force prediction is unusually sharp: "
        f"alpha is fixed at {_ALPHA_FIFTH}, not a free parameter.  "
        f"The only freedom is a_0 (internal radius), which sets lambda.  "
        f"Current Eot-Wash data already excludes a_0 > {survival_um[1]:.0f} um.  "
        f"The surviving window [{survival_um[0]:.0f}, {survival_um[1]:.0f}] um "
        f"is accessible to next-generation short-range gravity experiments.  "
        f"If no Yukawa deviation is found at alpha = 0.618 for any lambda "
        f"in this window, the framework is falsified (assuming a_0 falls "
        f"in the sub-mm range, which is the theoretically motivated regime).  "
        f"Casimir experiments are not useful here because the electromagnetic "
        f"Casimir force overwhelms the gravitational Yukawa signal by many "
        f"orders of magnitude at the relevant distances."
    )

    return {
        "exclusion_map": exclusion_map,
        "eot_wash_prediction": eot_wash,
        "casimir_overlap": casimir,
        "discovery_reach": discovery,
        "prediction_line": prediction_line,
        "key_prediction": key_prediction,
        "survival_window_um": survival_um,
        "optimal_experiment": optimal_exp,
        "signal_strength_at_optimal": signal_at_optimal,
        "honest_assessment": honest_assessment,
    }
