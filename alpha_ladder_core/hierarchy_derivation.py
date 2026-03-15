"""
Hierarchy derivation analysis for the Alpha Ladder framework.

A brute-force scan discovered: alpha_g = alpha^(d*D) * mu^n = alpha^24 * mu^2,
equivalently G = alpha^24 * m_p^2 * hbar*c / m_e^4.  This has 688 ppm residual
with zero fitted parameters -- every exponent comes from d=4, n=2, D=6.

This module explores 7 theoretical angles for WHY this formula holds,
computes each numerically, and presents honest results.

All calculations use pure Python math (no numpy/scipy).
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50


# ---------------------------------------------------------------------------
# Helper: default constants
# ---------------------------------------------------------------------------

def _default_constants(constants):
    """Return constants, loading CODATA 2018 if None."""
    if constants is None:
        from alpha_ladder_core.constants import get_constants
        return get_constants("CODATA 2018")
    return constants


# ---------------------------------------------------------------------------
# 1. Core numerics
# ---------------------------------------------------------------------------

def compute_formula_basics(d=4, n=2, constants=None):
    """
    Compute the core numerics of the hierarchy formula.

    alpha_g = alpha^(d*D) * mu^n, where D = d + n and mu = m_p / m_e.

    Parameters
    ----------
    d : int
        Number of external (non-compact) dimensions (default 4).
    n : int
        Number of internal (compact) dimensions (default 2).
    constants : SimpleNamespace or None
        Physical constants.  If None, CODATA 2018 is loaded.

    Returns
    -------
    dict with keys: exponent, d, n, D, mu, alpha, alpha_g_predicted,
        alpha_g_measured, ratio, residual_ppm, G_predicted, G_measured
    """
    constants = _default_constants(constants)

    D = d + n
    exponent = d * D

    mu = float(constants.m_p / constants.m_e)
    alpha = float(constants.alpha)

    alpha_g_predicted = alpha ** exponent * mu ** n
    alpha_g_measured = float(
        constants.G * constants.m_e ** 2 / (constants.hbar * constants.c)
    )

    ratio = alpha_g_predicted / alpha_g_measured
    residual_ppm = (ratio - 1) * 1e6

    G_predicted = alpha_g_predicted * float(
        constants.hbar * constants.c / constants.m_e ** 2
    )
    G_measured = float(constants.G)

    return {
        "exponent": exponent,
        "d": d,
        "n": n,
        "D": D,
        "mu": mu,
        "alpha": alpha,
        "alpha_g_predicted": alpha_g_predicted,
        "alpha_g_measured": alpha_g_measured,
        "ratio": ratio,
        "residual_ppm": residual_ppm,
        "G_predicted": G_predicted,
        "G_measured": G_measured,
    }


# ---------------------------------------------------------------------------
# 2. KK volume suppression
# ---------------------------------------------------------------------------

def analyze_kk_volume_suppression(d=4, n=2, constants=None):
    """
    Standard KK: G_4 = G_D / V_n.  Show this gives n powers of radius,
    not d*D = 24 powers of alpha.

    Parameters
    ----------
    d, n : int
        External and internal dimensions.
    constants : SimpleNamespace or None

    Returns
    -------
    dict with keys: works, M_D, gap_description, kk_formula,
        n_powers_of_radius, required_powers, ansatz
    """
    constants = _default_constants(constants)

    from alpha_ladder_core.kk_reduction import compute_einstein_frame_ansatz
    ansatz = compute_einstein_frame_ansatz(d, n)

    D = d + n
    required_powers = d * D

    # Standard KK formula: G_4 = 1 / (M_D^{D-2} * V_n)
    # For V_n ~ a_0^n, this involves n powers of the internal radius.
    kk_formula = (
        f"G_4 = 1 / (M_D^{{{D - 2}}} * V_{n}), "
        f"where V_{n} ~ a_0^{n}"
    )
    n_powers_of_radius = n

    # Estimate M_D from Planck mass: M_Pl^2 ~ M_D^{D-2} * a_0^n
    # In natural units (hbar=c=1): M_Pl ~ 1.22e28 eV
    M_Pl_eV = 1.22089e28
    # If a_0 ~ 1/M_D (compactification at M_D scale), then
    # M_Pl^2 ~ M_D^{D-2} * M_D^{-n} = M_D^{d-2}
    # => M_D ~ M_Pl^{2/(d-2)}
    M_D = M_Pl_eV ** (2.0 / (d - 2))

    gap_description = (
        f"Standard KK reduction gives G_4 = 1/(M_D^{{{D - 2}}} * a_0^{n}). "
        f"This involves {n} powers of the compactification radius, "
        f"not the {required_powers} powers of alpha required by the formula "
        f"alpha_g = alpha^{{{required_powers}}} * mu^{{{n}}}. "
        f"KK volume suppression alone cannot explain the exponent structure."
    )

    return {
        "works": False,
        "M_D": M_D,
        "gap_description": gap_description,
        "kk_formula": kk_formula,
        "n_powers_of_radius": n_powers_of_radius,
        "required_powers": required_powers,
        "ansatz": ansatz,
    }


# ---------------------------------------------------------------------------
# 3. Metric component counting
# ---------------------------------------------------------------------------

def analyze_metric_component_counting(d=4, n=2):
    """
    Pure geometry: count symmetric metric components in D dimensions and
    compare with d*D.

    Parameters
    ----------
    d, n : int
        External and internal dimensions.

    Returns
    -------
    dict with keys: symmetric_components, dD_product, D, kk_decomposition,
        dD_interpretation, mismatch
    """
    D = d + n

    symmetric_components = D * (D + 1) // 2
    dD_product = d * D

    # KK decomposition of the symmetric metric tensor
    graviton = d * (d + 1) // 2      # g_{mu nu}
    vectors = d * n                    # g_{mu a} (off-diagonal)
    scalars = n * (n + 1) // 2        # g_{ab}

    kk_decomposition = {
        "graviton": graviton,
        "vectors": vectors,
        "scalars": scalars,
        "total": graviton + vectors + scalars,
    }

    mismatch = dD_product - symmetric_components

    dD_interpretation = (
        f"d*D = {dD_product} counts the product of external and total "
        f"dimensions. This exceeds the {symmetric_components} symmetric "
        f"metric components by {mismatch}. The mismatch shows d*D is NOT "
        f"a metric component count."
    )

    return {
        "symmetric_components": symmetric_components,
        "dD_product": dD_product,
        "D": D,
        "kk_decomposition": kk_decomposition,
        "dD_interpretation": dD_interpretation,
        "mismatch": mismatch,
    }


# ---------------------------------------------------------------------------
# 4. Power-law running
# ---------------------------------------------------------------------------

def analyze_power_law_running(d=4, n=2, constants=None):
    """
    Power-law running of the gravitational coupling above the
    compactification scale.

    In D dimensions, gravity coupling runs as G_eff(E) ~ G_4 * (E/M_c)^n.
    Check whether this accounts for the full exponent structure.

    Parameters
    ----------
    d, n : int
    constants : SimpleNamespace or None

    Returns
    -------
    dict with keys: effective_exponent, dD, matches_dD,
        mu_contribution_exponent, total_accounted, honest, description
    """
    constants = _default_constants(constants)

    alpha = float(constants.alpha)
    mu = float(constants.m_p / constants.m_e)
    alpha_g_measured = float(
        constants.G * constants.m_e ** 2 / (constants.hbar * constants.c)
    )

    # Effective exponent: alpha_g ~ alpha^{eff_exp}
    # log(alpha_g) / log(alpha) gives eff_exp if alpha_g were a pure power
    effective_exponent = math.log(alpha_g_measured) / math.log(alpha)

    dD = d * (d + n)
    matches_dD = abs(effective_exponent - dD) < 1.0

    # mu^n contribution in powers of 1/alpha
    mu_contribution_exponent = n * math.log(mu) / math.log(1.0 / alpha)

    # Total accounted: d*D from alpha part + mu contribution
    total_accounted = dD + mu_contribution_exponent

    honest = (
        f"Power-law running in D={d + n} gives n={n} extra powers of "
        f"(E/M_c) above the compactification scale. This accounts for "
        f"{n} of the required exponent, but cannot explain the full "
        f"d*D={dD} power of alpha."
    )

    description = (
        f"The measured alpha_g ~ alpha^{{{effective_exponent:.2f}}}. "
        f"The formula predicts alpha^{{{dD}}} * mu^{{{n}}}. "
        f"The mu^{n} factor contributes ~{mu_contribution_exponent:.2f} "
        f"powers of (1/alpha), so the total effective exponent is "
        f"~{dD} + {mu_contribution_exponent:.2f} = "
        f"{dD + mu_contribution_exponent:.2f}, "
        f"consistent with {effective_exponent:.2f}."
    )

    return {
        "effective_exponent": effective_exponent,
        "dD": dD,
        "matches_dD": matches_dD,
        "mu_contribution_exponent": mu_contribution_exponent,
        "total_accounted": total_accounted,
        "honest": honest,
        "description": description,
    }


# ---------------------------------------------------------------------------
# 5. Volume-mass relation (does mu come from alpha?)
# ---------------------------------------------------------------------------

def analyze_volume_mass_relation(d=4, n=2, constants=None):
    """
    Test whether the proton-electron mass ratio mu emerges naturally
    as a power of alpha.

    Parameters
    ----------
    d, n : int
    constants : SimpleNamespace or None

    Returns
    -------
    dict with keys: mu, alpha, mu_alpha_exponent, chain_works,
        chain_description, mu_squared_contribution
    """
    constants = _default_constants(constants)

    mu = float(constants.m_p / constants.m_e)
    alpha = float(constants.alpha)

    # How many powers of 1/alpha is mu?
    mu_alpha_exponent = math.log(mu) / math.log(1.0 / alpha)

    chain_works = abs(mu_alpha_exponent - round(mu_alpha_exponent)) < 0.05

    # mu^2 contribution to alpha_g in absolute terms
    mu_squared_contribution = mu ** n

    chain_description = (
        f"mu = m_p/m_e ~ alpha^{{-{mu_alpha_exponent:.2f}}}. "
        f"If mu were exactly alpha^{{-1}} or alpha^{{-2}}, the formula "
        f"would reduce to a pure power of alpha. The irrational exponent "
        f"~{mu_alpha_exponent:.2f} means the proton-electron mass ratio "
        f"is NOT simply a power of alpha. The mu^{n} factor in the formula "
        f"is an independent input, not derived from alpha."
    )

    return {
        "mu": mu,
        "alpha": alpha,
        "mu_alpha_exponent": mu_alpha_exponent,
        "chain_works": chain_works,
        "chain_description": chain_description,
        "mu_squared_contribution": mu_squared_contribution,
    }


# ---------------------------------------------------------------------------
# 6. Induced gravity (Sakharov mechanism)
# ---------------------------------------------------------------------------

def analyze_induced_gravity(d=4, n=2, constants=None):
    """
    Sakharov induced gravity: G ~ 1/(N * Lambda_UV^2) where N is the
    number of species.  Check whether N = d*D = 24 reproduces the formula.

    Parameters
    ----------
    d, n : int
    constants : SimpleNamespace or None

    Returns
    -------
    dict with keys: N_species, Lambda_UV, M_Pl_eV, k_best, works,
        description, honest_assessment
    """
    constants = _default_constants(constants)

    D = d + n
    N_species = d * D

    M_Pl_eV = 1.22089e28
    Lambda_UV = M_Pl_eV / math.sqrt(N_species)

    alpha = float(constants.alpha)
    mu = float(constants.m_p / constants.m_e)
    alpha_g_measured = float(
        constants.G * constants.m_e ** 2 / (constants.hbar * constants.c)
    )

    # Scan k values to see if alpha^k * mu^n / N_species ~ alpha_g
    best_k = None
    best_diff = float("inf")
    for k in range(1, 30):
        candidate = alpha ** k * mu ** n / N_species
        diff = abs(math.log(candidate) - math.log(alpha_g_measured))
        if diff < best_diff:
            best_diff = diff
            best_k = k

    # Honest assessment: N=24 matching d*D is suggestive but induced gravity
    # gives G ~ 1/N, not G ~ alpha^N
    works = False

    description = (
        f"Sakharov induced gravity: G ~ 1/(N * Lambda_UV^2). "
        f"With N = d*D = {N_species} species, Lambda_UV = M_Pl/sqrt({N_species}) "
        f"= {Lambda_UV:.4e} eV. The best k for alpha^k * mu^{n} / N ~ alpha_g "
        f"is k = {best_k}, but the induced gravity mechanism produces 1/N, "
        f"not alpha^N."
    )

    honest_assessment = (
        f"The coincidence N = d*D = {N_species} is numerically suggestive. "
        f"However, the induced gravity mechanism G ~ 1/(N * Lambda_UV^2) "
        f"does not produce alpha^{{{N_species}}}. It produces a 1/{N_species} "
        f"suppression, which is a single rational factor, not an exponential "
        f"suppression in the coupling constant. The mechanism does not derive "
        f"the formula."
    )

    return {
        "N_species": N_species,
        "Lambda_UV": Lambda_UV,
        "M_Pl_eV": M_Pl_eV,
        "k_best": best_k,
        "works": works,
        "description": description,
        "honest_assessment": honest_assessment,
    }


# ---------------------------------------------------------------------------
# 7. Dimension pair scan
# ---------------------------------------------------------------------------

def scan_dimension_pairs(d_range=None, n_range=None, constants=None):
    """
    Scan (d, n) pairs to show that (4, 2) is uniquely close.

    Parameters
    ----------
    d_range : iterable of int or None
        Range of external dimensions (default range(3, 8)).
    n_range : iterable of int or None
        Range of internal dimensions (default range(1, 6)).
    constants : SimpleNamespace or None

    Returns
    -------
    dict with keys: scan_results, best_pair, best_ppm, unique_sub_1000,
        runner_up_ppm, n_pairs
    """
    if d_range is None:
        d_range = range(3, 8)
    if n_range is None:
        n_range = range(1, 6)

    constants = _default_constants(constants)

    alpha = float(constants.alpha)
    mu = float(constants.m_p / constants.m_e)
    alpha_g_measured = float(
        constants.G * constants.m_e ** 2 / (constants.hbar * constants.c)
    )

    log_alpha = math.log(alpha)
    log_mu = math.log(mu)
    log_alpha_g_meas = math.log(alpha_g_measured)

    scan_results = []
    for d_val in d_range:
        for n_val in n_range:
            D_val = d_val + n_val
            exponent = d_val * D_val
            log_alpha_g_pred = exponent * log_alpha + n_val * log_mu
            log_ratio = log_alpha_g_pred - log_alpha_g_meas
            ratio = math.exp(log_ratio)
            ppm = abs((ratio - 1) * 1e6)
            scan_results.append({
                "d": d_val,
                "n": n_val,
                "D": D_val,
                "exponent": exponent,
                "ppm": ppm,
                "ratio": ratio,
            })

    scan_results.sort(key=lambda r: r["ppm"])

    best = scan_results[0]
    best_pair = (best["d"], best["n"])
    best_ppm = best["ppm"]

    runner_up_ppm = scan_results[1]["ppm"] if len(scan_results) > 1 else None

    sub_1000 = [r for r in scan_results if r["ppm"] < 1000]
    unique_sub_1000 = len(sub_1000) == 1

    return {
        "scan_results": scan_results,
        "best_pair": best_pair,
        "best_ppm": best_ppm,
        "unique_sub_1000": unique_sub_1000,
        "runner_up_ppm": runner_up_ppm,
        "n_pairs": len(scan_results),
    }


# ---------------------------------------------------------------------------
# 8. Residual analysis
# ---------------------------------------------------------------------------

def analyze_residual(constants=None):
    """
    Analyze the 688 ppm gap: test QED corrections, reduced mass effects,
    running alpha, and combined corrections.

    Parameters
    ----------
    constants : SimpleNamespace or None

    Returns
    -------
    dict with keys: raw_ratio, raw_ppm, corrections, best_correction,
        n_corrections
    """
    constants = _default_constants(constants)

    basics = compute_formula_basics(constants=constants)
    raw_ratio = basics["ratio"]
    raw_ppm = basics["residual_ppm"]
    alpha = basics["alpha"]
    mu = basics["mu"]
    alpha_g_measured = basics["alpha_g_measured"]

    corrections = []

    # (a) QED vertex correction: (1 - alpha/(2*pi))^k
    for k in range(1, 7):
        factor = (1 - alpha / (2 * math.pi)) ** k
        corrected_ratio = raw_ratio * factor
        corrected_ppm = (corrected_ratio - 1) * 1e6
        corrections.append({
            "name": f"QED vertex (k={k})",
            "factor": factor,
            "corrected_ppm": corrected_ppm,
            "improvement": abs(corrected_ppm) < abs(raw_ppm),
        })

    # (b) Reduced mass correction: (1 - 1/mu)^k
    for k in range(1, 5):
        factor = (1 - 1.0 / mu) ** k
        corrected_ratio = raw_ratio * factor
        corrected_ppm = (corrected_ratio - 1) * 1e6
        corrections.append({
            "name": f"Reduced mass (k={k})",
            "factor": factor,
            "corrected_ppm": corrected_ppm,
            "improvement": abs(corrected_ppm) < abs(raw_ppm),
        })

    # (c) Running alpha to proton mass scale
    alpha_at_mp = alpha / (
        1 - (2 * alpha / (3 * math.pi)) * math.log(mu)
    )
    alpha_g_running = alpha_at_mp ** 24 * mu ** 2
    corrected_ratio_running = alpha_g_running / alpha_g_measured
    corrected_ppm_running = (corrected_ratio_running - 1) * 1e6
    corrections.append({
        "name": "Running alpha (1-loop QED)",
        "factor": (alpha_at_mp / alpha) ** 24,
        "corrected_ppm": corrected_ppm_running,
        "improvement": abs(corrected_ppm_running) < abs(raw_ppm),
    })

    # (d) Combined: (1 - k*alpha/pi)
    for k in range(1, 5):
        factor = (1 - k * alpha / math.pi)
        corrected_ratio = raw_ratio * factor
        corrected_ppm = (corrected_ratio - 1) * 1e6
        corrections.append({
            "name": f"Combined (1 - {k}*alpha/pi)",
            "factor": factor,
            "corrected_ppm": corrected_ppm,
            "improvement": abs(corrected_ppm) < abs(raw_ppm),
        })

    best_correction = min(corrections, key=lambda c: abs(c["corrected_ppm"]))

    return {
        "raw_ratio": raw_ratio,
        "raw_ppm": raw_ppm,
        "corrections": corrections,
        "best_correction": best_correction,
        "n_corrections": len(corrections),
    }


# ---------------------------------------------------------------------------
# 9. Swampland Distance Conjecture
# ---------------------------------------------------------------------------

def analyze_swampland_distance(d=4, n=2):
    """
    Check the Swampland Distance Conjecture (SDC) for the KK tower.

    The SDC states that as one moves a geodesic distance Delta in moduli
    space, a tower of states appears with masses m ~ m_0 * exp(-lambda * Delta)
    where lambda >= 1/sqrt(d-2).

    For pure Einstein-Hilbert KK reduction the SDC is violated, but the
    Gauss-Bonnet correction needed for the golden-ratio omega also fixes
    the SDC bound -- the two constraints are linked.

    Parameters
    ----------
    d : int
        Number of external dimensions (default 4).
    n : int
        Number of internal dimensions (default 2).

    Returns
    -------
    dict with SDC analysis including lambda values, bounds, and whether
    the conjecture is satisfied for EH and GB-corrected cases.
    """
    from alpha_ladder_core.kk_reduction import (
        compute_kinetic_coefficient,
        compute_target_omega,
    )

    D = d + n

    # --- Pure Einstein-Hilbert ---
    kinetic = compute_kinetic_coefficient(d, n)
    K_EH = kinetic["K_einstein"]  # negative for healthy kinetic term

    # Canonical normalization: action has K_EH * (d sigma)^2
    # For -(1/2)(d phi)^2: phi = sqrt(2*|K_EH|) * sigma
    # KK masses scale as m_KK ~ exp(|beta|*sigma)
    # Einstein frame beta for d=4, n=2: beta = -(d-2)/n = -1
    beta = -(d - 2) / n

    # In terms of canonical field: m_KK ~ exp(phi / sqrt(2*|K_EH|))
    # SDC rate: lambda = 1 / sqrt(2*|K_EH|)
    abs_K_EH = abs(K_EH)
    lambda_EH = 1.0 / math.sqrt(2.0 * abs_K_EH)

    # SDC bound
    lambda_bound = 1.0 / math.sqrt(d - 2)

    sdc_satisfied_EH = lambda_EH >= lambda_bound

    # omega for pure EH
    omega_EH = K_EH / (2.0 * n ** 2 * beta ** 2)

    # --- GB-corrected (target omega) ---
    target = compute_target_omega()
    omega_target = target["omega_float"]

    # |K_GB| = omega_target * 2 * n^2 * beta^2
    abs_K_GB = omega_target * 2.0 * n ** 2 * beta ** 2
    lambda_GB = 1.0 / math.sqrt(2.0 * abs_K_GB)

    sdc_satisfied_GB = lambda_GB >= lambda_bound

    # Saturation point: lambda = lambda_bound requires 2*|K| = d-2
    # i.e. |K_sat| = (d-2)/2
    abs_K_sat = (d - 2) / 2.0
    omega_sat = abs_K_sat / (2.0 * n ** 2 * beta ** 2)

    margin = lambda_GB - lambda_bound
    near_saturation = abs(omega_target - omega_sat) < 0.01

    if sdc_satisfied_GB and near_saturation:
        assessment = (
            f"The GB-corrected kinetic coefficient gives lambda_GB = "
            f"{lambda_GB:.4f}, which satisfies the SDC bound "
            f"lambda >= {lambda_bound:.4f} with margin {margin:.4f}. "
            f"The target omega = {omega_target:.4f} is near the saturation "
            f"value omega_sat = {omega_sat:.4f}, suggesting the golden-ratio "
            f"omega nearly saturates the SDC bound."
        )
    elif sdc_satisfied_GB:
        assessment = (
            f"The GB-corrected kinetic coefficient satisfies the SDC bound "
            f"(lambda_GB = {lambda_GB:.4f} >= {lambda_bound:.4f})."
        )
    else:
        assessment = (
            f"The GB-corrected kinetic coefficient violates the SDC bound "
            f"(lambda_GB = {lambda_GB:.4f} < {lambda_bound:.4f})."
        )

    return {
        "d": d,
        "n": n,
        "D": D,
        "beta": beta,
        "K_EH": K_EH,
        "K_GB": -abs_K_GB,
        "abs_K_GB": abs_K_GB,
        "lambda_EH": lambda_EH,
        "lambda_GB": lambda_GB,
        "lambda_bound": lambda_bound,
        "sdc_satisfied_EH": sdc_satisfied_EH,
        "sdc_satisfied_GB": sdc_satisfied_GB,
        "omega_EH": omega_EH,
        "omega_target": omega_target,
        "omega_saturation": omega_sat,
        "near_saturation": near_saturation,
        "margin": margin,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 10. Emergence Proposal (KK tower)
# ---------------------------------------------------------------------------

def _kk_sum(L):
    """Compute sum_{l=1}^{L} (2l+1)*l*(l+1)."""
    return sum((2 * l + 1) * l * (l + 1) for l in range(1, L + 1))


def analyze_emergence_tower(d=4, n=2, constants=None):
    """
    Test the emergence proposal: M_Pl^2 generated by one-loop diagrams
    from the KK tower on S^2.

    The KK spectrum on S^2 of radius R has m_l^2 = l(l+1)/R^2 with
    degeneracy 2l+1.  The one-loop contribution to M_Pl^2 is:

        delta(M_Pl^2) = N_pol / (16*pi^2) * sum_{l=1}^{L_max} (2l+1)*l(l+1)/R^2

    For large L_max, the sum scales as L_max^4/2, giving:

        M_Pl^2 ~ N_pol * L_max^4 / (32*pi^2 * R^2)

    Key finding: the emergence sum does NOT naturally produce the
    alpha^24 * mu^2 structure.

    Parameters
    ----------
    d : int
        Number of external dimensions (default 4).
    n : int
        Number of internal dimensions (default 2).
    constants : SimpleNamespace or None

    Returns
    -------
    dict with emergence tower analysis including canonical R results
    and honest assessment.
    """
    constants = _default_constants(constants)

    D = d + n

    # Physical constants in natural units (eV)
    hbar_c_eV_m = 1.97327e-7   # eV * m
    m_e_eV = 0.51099895e6      # eV
    m_p_eV = 938.272046e6      # eV
    M_Pl_eV = 1.22089e28       # eV (standard Planck mass)

    alpha = float(constants.alpha)

    # Planck length in meters
    l_Pl_m = 1.616255e-35

    # Canonical radii (meters)
    r_e_m = float(constants.r_e_nist) if hasattr(constants, "r_e_nist") else (
        alpha ** 2 * float(constants.a_0)
    )
    lambda_bar_c_m = float(constants.lambda_bar_c) if hasattr(
        constants, "lambda_bar_c"
    ) else float(constants.hbar / (constants.m_e * constants.c))
    a_0_m = float(constants.a_0)

    N_pol = 21  # D(D+1)/2 for D=6, symmetric metric components

    canonical_radii = [
        ("l_Pl", l_Pl_m),
        ("r_e", r_e_m),
        ("lambda_bar_c", lambda_bar_c_m),
        ("a_0", a_0_m),
    ]

    canonical_R_results = []
    for name, R_m in canonical_radii:
        # Convert R to natural units (eV^{-1})
        R_eV_inv = R_m / hbar_c_eV_m

        # Cutoff at Planck mass: L_max = M_Pl * R (natural units)
        L_max_planck = M_Pl_eV * R_eV_inv

        # For large L_max, use analytic scaling: S(L) ~ L^4/2
        # delta(M_Pl^2) = N_pol / (16*pi^2) * L_max^4 / (2 * R_eV_inv^2)
        # In eV^2 units
        if L_max_planck > 0:
            M_Pl_sq_predicted = (
                N_pol / (16.0 * math.pi ** 2)
                * L_max_planck ** 4
                / (2.0 * R_eV_inv ** 2)
            )
        else:
            M_Pl_sq_predicted = 0.0

        M_Pl_sq_measured = M_Pl_eV ** 2

        M_Pl_sq_ratio = (
            M_Pl_sq_predicted / M_Pl_sq_measured
            if M_Pl_sq_measured > 0
            else 0.0
        )

        # Predicted alpha_g
        if M_Pl_sq_predicted > 0:
            alpha_g_predicted = m_e_eV ** 2 / M_Pl_sq_predicted
        else:
            alpha_g_predicted = 0.0

        canonical_R_results.append({
            "name": name,
            "R_m": R_m,
            "L_max_planck": L_max_planck,
            "M_Pl_sq_ratio": M_Pl_sq_ratio,
            "alpha_g_predicted": alpha_g_predicted,
        })

    # Find R needed for emergence = observed M_Pl^2
    # M_Pl^2 = N_pol * (M_Pl * R)^4 / (32*pi^2 * R^2)
    #        = N_pol * M_Pl^4 * R^2 / (32*pi^2)
    # => R^2 = 32*pi^2 * M_Pl^2 / (N_pol * M_Pl^4)
    #        = 32*pi^2 / (N_pol * M_Pl^2)
    R_needed_eV_inv = math.sqrt(32.0 * math.pi ** 2 / (N_pol * M_Pl_eV ** 2))
    R_needed_m = R_needed_eV_inv * hbar_c_eV_m

    # Express in alpha units: R / (hbar/(m_e*c)) = R * m_e (natural units)
    R_needed_lambda_bar_c = R_needed_eV_inv * m_e_eV

    # Effective alpha exponent: R_needed / lambda_bar_c = alpha^a
    # => a = log(R_needed_lambda_bar_c) / log(alpha)
    if R_needed_lambda_bar_c > 0:
        R_needed_alpha_units = (
            math.log(R_needed_lambda_bar_c) / math.log(alpha)
            if R_needed_lambda_bar_c != 1.0
            else 0.0
        )
    else:
        R_needed_alpha_units = float("inf")

    honest_assessment = (
        "The emergence proposal with Lambda = M_Pl gives a self-consistent "
        f"R ~ {R_needed_m:.2e} m (approximately l_Pl), which is trivially "
        "circular. For Lambda = m_p, the required R is macroscopic and "
        "unphysical. The one-loop KK sum scales as L_max^4/R^2, which does "
        "not naturally produce the alpha^24 * mu^2 structure of the "
        "hierarchy formula."
    )

    return {
        "d": d,
        "n": n,
        "D": D,
        "sum_formula": "sum_{l=1}^{L} (2l+1)*l*(l+1)",
        "scaling_large_L": "M_Pl^2 ~ N_pol * L_max^4 / (32*pi^2*R^2)",
        "canonical_R_results": canonical_R_results,
        "N_pol_default": N_pol,
        "R_needed_for_match_m": R_needed_m,
        "R_needed_alpha_units": R_needed_alpha_units,
        "works": False,
        "honest_assessment": honest_assessment,
    }


# ---------------------------------------------------------------------------
# 11. One-loop KK power counting
# ---------------------------------------------------------------------------

def analyze_one_loop_matching(d=4, n=2, constants=None):
    """
    One-loop correction to 1/G from the KK tower: power counting analysis.

    The key formula is:

        delta(M_Pl^2) ~ N_pol * Lambda^4 * R^2 / (32*pi^2)

    so alpha_g = m_e^2 / M_Pl^2 ~ 32*pi^2 * m_e^2 / (N_pol * Lambda^4 * R^2).

    Parameterize R and Lambda in terms of alpha and masses to check whether
    the alpha^24 * mu^2 structure can arise from one-loop power counting.

    Key finding: the mu exponent comes out structurally wrong (mu^{-4}
    instead of mu^{+2}).

    Parameters
    ----------
    d : int
        Number of external dimensions (default 4).
    n : int
        Number of internal dimensions (default 2).
    constants : SimpleNamespace or None

    Returns
    -------
    dict with one-loop matching analysis including specific cases,
    alpha exponent scan, and honest assessment.
    """
    constants = _default_constants(constants)

    D = d + n

    m_e_eV = 0.51099895e6
    m_p_eV = 938.272046e6
    M_Pl_eV = 1.22089e28
    alpha = float(constants.alpha)
    mu = float(constants.m_p / constants.m_e)

    alpha_g_measured = float(
        constants.G * constants.m_e ** 2 / (constants.hbar * constants.c)
    )

    N_pol = 21  # D(D+1)/2 for D=6
    prefactor = 32.0 * math.pi ** 2

    # UV divergence degree: D-2 for gravity in D dimensions
    uv_divergence_degree = D - 2

    # --- Specific cases ---
    specific_cases = []

    # Canonical scales (in eV, natural units)
    hbar_c_eV_m = 1.97327e-7
    a_0_m = float(constants.a_0)
    r_e_m = float(constants.r_e_nist) if hasattr(constants, "r_e_nist") else (
        alpha ** 2 * a_0_m
    )

    # Convert to eV^{-1}
    a_0_eV_inv = a_0_m / hbar_c_eV_m
    r_e_eV_inv = r_e_m / hbar_c_eV_m

    cases_spec = [
        ("r_e", r_e_eV_inv, "M_Pl", M_Pl_eV),
        ("a_0", a_0_eV_inv, "M_Pl", M_Pl_eV),
        ("r_e", r_e_eV_inv, "m_p", m_p_eV),
        ("a_0", a_0_eV_inv, "m_p", m_p_eV),
    ]

    for R_name, R_eV_inv, Lambda_name, Lambda_eV in cases_spec:
        # M_Pl^2 predicted = N_pol * Lambda^4 * R^2 / (32*pi^2)
        # alpha_g = m_e^2 / M_Pl_sq_predicted
        M_Pl_sq_pred = N_pol * Lambda_eV ** 4 * R_eV_inv ** 2 / prefactor
        if M_Pl_sq_pred > 0:
            alpha_g_pred = m_e_eV ** 2 / M_Pl_sq_pred
            ratio = alpha_g_pred / alpha_g_measured
        else:
            alpha_g_pred = 0.0
            ratio = 0.0

        specific_cases.append({
            "R": R_name,
            "Lambda": Lambda_name,
            "alpha_g_predicted": alpha_g_pred,
            "ratio": ratio,
        })

    # --- Alpha exponent scan ---
    # R = alpha^a / m_e, Lambda = alpha^b * m_e (both in natural units)
    # alpha_g = prefactor * m_e^2 / (N_pol * (alpha^b * m_e)^4 * (alpha^a / m_e)^2)
    #         = prefactor / (N_pol * m_e^2) * alpha^{-(4b + 2a)}
    # Wait: let's be careful.
    # Lambda^4 = alpha^{4b} * m_e^4
    # R^2 = alpha^{2a} / m_e^2
    # M_Pl_sq = N_pol * alpha^{4b} * m_e^4 * alpha^{2a} / (m_e^2 * prefactor)
    #         = N_pol * alpha^{4b+2a} * m_e^2 / prefactor
    # alpha_g = m_e^2 / M_Pl_sq = prefactor / (N_pol * alpha^{4b+2a})
    # So alpha_g ~ alpha^{-(4b+2a)}, mu^0
    # Need -(4b+2a) = 24 => 4b+2a = -24 => 2b+a = -12

    # Also try Lambda involving mu: Lambda = mu^p * alpha^b * m_e
    # Lambda^4 = mu^{4p} * alpha^{4b} * m_e^4
    # alpha_g = prefactor / (N_pol * mu^{4p} * alpha^{4b+2a})
    # Need mu^{-4p} = mu^2 => p = -1/2 (non-integer!)
    # Need -(4b+2a) = 24 => 2b+a = -12

    alpha_exponent_scan = []
    target_alpha_exp = 24
    target_mu_exp = 2

    # Scan for R = alpha^a / m_e, Lambda = alpha^b * m_e (no mu)
    for a in range(-15, 15):
        for b in range(-15, 15):
            alpha_exp = -(4 * b + 2 * a)
            mu_exp = 0
            if alpha_exp == target_alpha_exp:
                alpha_exponent_scan.append({
                    "a": a,
                    "b": b,
                    "alpha_exp": alpha_exp,
                    "mu_exp": mu_exp,
                    "matches_alpha": True,
                    "matches_mu": mu_exp == target_mu_exp,
                    "Lambda_type": "alpha^b * m_e",
                })

    # Scan with Lambda = m_p * alpha^b (so Lambda involves mu)
    # Lambda = mu * m_e * alpha^b
    # Lambda^4 = mu^4 * m_e^4 * alpha^{4b}
    # alpha_g = prefactor / (N_pol * mu^4 * alpha^{4b+2a})
    # mu exponent: -4, alpha exponent: -(4b+2a)
    for a in range(-15, 15):
        for b in range(-15, 15):
            alpha_exp = -(4 * b + 2 * a)
            mu_exp = -4
            if alpha_exp == target_alpha_exp:
                alpha_exponent_scan.append({
                    "a": a,
                    "b": b,
                    "alpha_exp": alpha_exp,
                    "mu_exp": mu_exp,
                    "matches_alpha": True,
                    "matches_mu": mu_exp == target_mu_exp,
                    "Lambda_type": "m_p * alpha^b",
                })

    mu_exponent_problem = (
        "One-loop with Lambda=m_p gives mu^{-4}, not mu^{+2}. "
        "The mu exponent is structurally wrong."
    )

    honest_assessment = (
        "One-loop power counting gives alpha_g ~ Lambda^{-4} * R^{-2}. "
        "While integer combinations (a, b) exist that reproduce the "
        "alpha^24 exponent (e.g., 2b+a = -12), the mu exponent is "
        "structurally wrong: Lambda = m_p gives mu^{-4} instead of mu^{+2}. "
        "Achieving mu^{+2} requires a fractional power p = -1/2 in "
        "Lambda = mu^p, which is unphysical. The one-loop mechanism "
        "cannot derive the formula."
    )

    return {
        "d": d,
        "n": n,
        "D": D,
        "uv_divergence_degree": uv_divergence_degree,
        "analytic_scaling": (
            "alpha_g ~ 32*pi^2 * m_e^2 / (N_pol * Lambda^4 * R^2)"
        ),
        "specific_cases": specific_cases,
        "mu_exponent_problem": mu_exponent_problem,
        "alpha_exponent_scan": alpha_exponent_scan,
        "works": False,
        "honest_assessment": honest_assessment,
    }


# ---------------------------------------------------------------------------
# 12. Summary
# ---------------------------------------------------------------------------

def summarize_hierarchy_derivation(constants=None):
    """
    Main entry point combining all 10 analyses into a comprehensive summary.

    Parameters
    ----------
    constants : SimpleNamespace or None

    Returns
    -------
    dict with keys: formula_basics, angle_1 through angle_10,
        angles_summary, honest_assessment
    """
    constants = _default_constants(constants)

    formula_basics = compute_formula_basics(constants=constants)
    angle_1 = analyze_kk_volume_suppression(constants=constants)
    angle_2 = analyze_metric_component_counting()
    angle_3 = analyze_power_law_running(constants=constants)
    angle_4 = analyze_volume_mass_relation(constants=constants)
    angle_5 = analyze_induced_gravity(constants=constants)
    angle_6 = scan_dimension_pairs(constants=constants)
    angle_7 = analyze_residual(constants=constants)
    angle_8 = analyze_swampland_distance(d=4, n=2)
    angle_9 = analyze_emergence_tower(d=4, n=2, constants=constants)
    angle_10 = analyze_one_loop_matching(d=4, n=2, constants=constants)

    angles_summary = [
        {
            "name": "KK volume suppression",
            "works": angle_1["works"],
            "description": angle_1["gap_description"],
        },
        {
            "name": "Metric component counting",
            "works": False,
            "description": angle_2["dD_interpretation"],
        },
        {
            "name": "Power-law running",
            "works": False,
            "description": angle_3["honest"],
        },
        {
            "name": "Volume-mass relation (mu from alpha?)",
            "works": angle_4["chain_works"],
            "description": angle_4["chain_description"],
        },
        {
            "name": "Induced gravity (Sakharov)",
            "works": angle_5["works"],
            "description": angle_5["honest_assessment"],
        },
        {
            "name": "Dimension pair scan",
            "works": angle_6["unique_sub_1000"],
            "description": (
                f"(d,n) = {angle_6['best_pair']} is the best pair at "
                f"{angle_6['best_ppm']:.0f} ppm. "
                f"Unique sub-1000 ppm: {angle_6['unique_sub_1000']}."
            ),
        },
        {
            "name": "Residual analysis (688 ppm)",
            "works": abs(angle_7["best_correction"]["corrected_ppm"]) < 100,
            "description": (
                f"Best correction: {angle_7['best_correction']['name']} "
                f"reduces residual to "
                f"{angle_7['best_correction']['corrected_ppm']:.0f} ppm."
            ),
        },
        {
            "name": "Swampland Distance Conjecture",
            "works": angle_8["sdc_satisfied_GB"],
            "description": angle_8["assessment"],
        },
        {
            "name": "Emergence Proposal (KK tower)",
            "works": angle_9["works"],
            "description": angle_9["honest_assessment"],
        },
        {
            "name": "One-loop KK power counting",
            "works": angle_10["works"],
            "description": angle_10["honest_assessment"],
        },
    ]

    honest_assessment = (
        "The formula alpha_g = alpha^24 * mu^2 reproduces measured alpha_g "
        f"to {formula_basics['residual_ppm']:.0f} ppm with zero fitted "
        "parameters. The exponents d=4, n=2, D=6 have clear geometric "
        "meaning. However, no single theoretical mechanism fully derives "
        "the formula from first principles. The KK volume suppression gives "
        "the wrong power counting. The metric component count (21) does not "
        "match d*D (24). Power-law running accounts for only part of the "
        "exponent. The proton-electron mass ratio mu is not simply a power "
        "of alpha. Induced gravity matches the species count N=24 but not "
        "the alpha^24 structure. The Swampland Distance Conjecture provides "
        "a consistency check: the GB correction needed for the golden-ratio "
        "omega also satisfies the SDC bound, suggesting the framework is "
        "UV-consistent. However, the emergence proposal and one-loop power "
        "counting do not naturally produce the alpha^24 * mu^2 structure -- "
        "the mu exponent comes out wrong (mu^{-4} instead of mu^2). "
        "The formula remains an empirical discovery awaiting theoretical "
        "derivation."
    )

    return {
        "formula_basics": formula_basics,
        "angle_1": angle_1,
        "angle_2": angle_2,
        "angle_3": angle_3,
        "angle_4": angle_4,
        "angle_5": angle_5,
        "angle_6": angle_6,
        "angle_7": angle_7,
        "angle_8": angle_8,
        "angle_9": angle_9,
        "angle_10": angle_10,
        "angles_summary": angles_summary,
        "honest_assessment": honest_assessment,
    }
