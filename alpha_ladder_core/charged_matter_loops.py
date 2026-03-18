"""
Charged matter one-loop corrections on S^2 with Dirac monopole background.

The existing oneloop_g_correction.py computes the one-loop correction to G
from PURE GRAVITY KK modes (graviton tower on S^2).  The feynman_diagram.py
computes the sigma -> KK photon -> sigma loop giving the (1-alpha) factor.

This module extends the loop calculation to include CHARGED MATTER fields
(scalars, fermions, vectors) propagating on S^2 in the presence of a Dirac
monopole background with magnetic charge n.

Monopole harmonics
------------------
On a round S^2 of radius R with a U(1) monopole of charge n, a scalar
field with U(1) charge q has KK modes labelled by j = j_min, j_min+1, ...
where j_min = |q*n| / 2.  The eigenvalues of the covariant Laplacian are:

    lambda_j = [j(j+1) - q^2*n^2/4] / R^2

with degeneracy d_j = 2j + 1.  For q*n = 0 (neutral fields), this reduces
to the standard spherical harmonics l(l+1)/R^2 starting from l = 0.

For charged fields (q*n != 0), the spectrum starts from j_min > 0 and the
zero mode is lifted.  This modifies the spectral zeta function and hence
the one-loop effective potential.

One-loop effective potential (Coleman-Weinberg)
-----------------------------------------------
For a field of mass m_j at KK level j, the CW contribution to the 4D
effective potential is:

    V_CW = (+/-) sum_j (2j+1) * m_j^4 / (64*pi^2) * [ln(m_j^2/mu^2) - c_s]

where +/- is for bosons/fermions and c_s is a scheme-dependent constant
(3/2 for scalars, 5/6 for vectors in MS-bar).

The key question: does the monopole-modified spectrum produce a different
spectral zeta value that could contribute to Lambda_6 or the G correction?

Matter content
--------------
From anomaly_cancellation.py, the anomaly-free groups provide:
    E8 x E8:    740 hypers, 496 vectors
    SO(32):     740 hypers, 496 vectors
    E7 x E7:    132 hypers, 266 vectors
    E6 x E7:    147 hypers, 211 vectors

6D multiplet field counting (on S^2):
    1 hypermultiplet = 4 real scalars + 2 Weyl fermions
    1 vector multiplet = 1 gauge boson + 2 Weyl fermions (+ 1 real scalar in some conventions)

All calculations use pure Python math (no numpy/scipy).
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50

from alpha_ladder_core.constants import get_constants


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

_R_GAUGE_MATCHED = 0.550293  # Planck lengths
_M_PL_EV = 1.22089e28


# ---------------------------------------------------------------------------
# 1. Monopole harmonic KK spectrum
# ---------------------------------------------------------------------------

def compute_monopole_kk_spectrum(charge_q, n_monopole=1, l_max=100):
    """
    Compute KK spectrum for a field with U(1) charge q in monopole
    background n on S^2.

    The covariant Laplacian eigenvalues are:

        lambda_j = j(j+1) - q^2*n^2/4

    for j = j_min, j_min+1, ..., j_min+l_max, where j_min = |q*n|/2.
    The mass-squared in Planck units is lambda_j / R^2.

    For neutral fields (q=0), this reduces to the standard spectrum
    l(l+1) starting from l=0.

    Parameters
    ----------
    charge_q : int or float
        U(1) charge of the field.
    n_monopole : int
        Monopole number (default 1, Dirac quantization).
    l_max : int
        Number of KK levels to include (default 100).

    Returns
    -------
    dict with keys:
        charge_q : float
        n_monopole : int
        j_min : float -- minimum angular momentum
        levels : list of dict -- each with j, eigenvalue, degeneracy
        n_levels : int
        spectrum_shift : str -- description of how monopole modifies spectrum
    """
    qn = abs(charge_q * n_monopole)
    j_min = qn / 2.0

    # For half-integer j_min, j takes half-integer values
    # For integer j_min, j takes integer values
    levels = []
    for k in range(l_max):
        j = j_min + k
        eigenvalue = j * (j + 1) - (qn / 2.0) ** 2
        degeneracy = int(2 * j + 1)

        levels.append({
            "j": j,
            "eigenvalue": eigenvalue,
            "degeneracy": degeneracy,
            "mass_sq_factor": eigenvalue,  # multiply by 1/R^2 for physical mass
        })

    shift_desc = "standard (neutral)" if qn == 0 else (
        f"monopole-shifted: j starts at {j_min}, zero mode lifted"
    )

    return {
        "charge_q": charge_q,
        "n_monopole": n_monopole,
        "j_min": j_min,
        "levels": levels,
        "n_levels": len(levels),
        "spectrum_shift": shift_desc,
    }


# ---------------------------------------------------------------------------
# 2. Spectral zeta for monopole harmonics
# ---------------------------------------------------------------------------

def _bernoulli_polynomial(n, x):
    """Bernoulli polynomial B_n(x) for n = 0..6."""
    if n == 0:
        return 1.0
    elif n == 1:
        return x - 0.5
    elif n == 2:
        return x ** 2 - x + 1.0 / 6.0
    elif n == 3:
        return x ** 3 - 1.5 * x ** 2 + 0.5 * x
    elif n == 4:
        return x ** 4 - 2.0 * x ** 3 + x ** 2 - 1.0 / 30.0
    elif n == 5:
        return x ** 5 - 2.5 * x ** 4 + (5.0 / 3.0) * x ** 3 - x / 6.0
    elif n == 6:
        return (x ** 6 - 3.0 * x ** 5 + 2.5 * x ** 4
                - 0.5 * x ** 2 + 1.0 / 42.0)
    else:
        raise ValueError(f"B_{n} not implemented (need n <= 6)")


def _hurwitz_zeta_neg_int(n, a):
    """Compute zeta_H(-n, a) = -B_{n+1}(a) / (n+1) for n >= 0."""
    return -_bernoulli_polynomial(n + 1, a) / (n + 1)


def compute_monopole_spectral_zeta(charge_q, n_monopole=1, s=-1,
                                    l_max=10000):
    """
    Compute the spectral zeta function for a charged field on S^2 with
    monopole background:

        zeta(s) = sum_{j >= j_min} (2j+1) * [j(j+1) - (qn/2)^2]^{-s}

    For s > 1, the sum converges and is computed directly.

    For s = -1, we use the Hurwitz zeta analytic continuation.  The
    substitution u = j + 1/2 gives:

        zeta(-1) = sum_{u >= u_min} 2u * [u^2 - (1/4 + c^2)]

    where u_min = j_min + 1/2 and c = qn/2.  Expanding:

        = 2 * [zeta_H(-3, u_min) - (1/4 + c^2) * zeta_H(-1, u_min)]

    with zeta_H(-n, a) = -B_{n+1}(a)/(n+1).  This is EXACT.

    Key results:
        - Neutral (q=0): zeta = -17/480 (matches oneloop_g_correction.py)
        - Charged (q=1, n=1): zeta = 1/10 exactly (SIGN FLIP!)

    For NEUTRAL fields (q=0), the j=0 term has eigenvalue 0 and is
    excluded (matching the convention where l starts from 1).

    Parameters
    ----------
    charge_q : int or float
        U(1) charge.
    n_monopole : int
        Monopole number.
    s : float
        Zeta function argument.
    l_max : int
        Cutoff for partial sums (used for s > 1).

    Returns
    -------
    dict with keys:
        s : float
        charge_q : float
        n_monopole : int
        j_min : float
        value : float -- the zeta-regularized value
        method : str
        converged : bool
        comparison_neutral : float or None
    """
    qn = abs(charge_q * n_monopole)
    j_min = qn / 2.0
    c_sq = (qn / 2.0) ** 2

    # For neutral fields, skip j=0 (eigenvalue = 0)
    start_j = j_min if j_min > 0 else 1.0

    if s > 1.0:
        # Direct summation (converges)
        total = 0.0
        j = start_j
        for _ in range(l_max):
            eig = j * (j + 1) - c_sq
            if eig > 0:
                total += (2 * j + 1) * eig ** (-s)
            j += 1.0

        return {
            "s": s,
            "charge_q": charge_q,
            "n_monopole": n_monopole,
            "j_min": j_min,
            "value": total,
            "method": "direct",
            "converged": True,
            "comparison_neutral": None,
        }

    if abs(s - (-1.0)) < 1e-12:
        # s = -1: EXACT analytic continuation via Hurwitz zeta
        #
        # Substitution: u = j + 1/2, so j(j+1) = u^2 - 1/4, 2j+1 = 2u.
        # u_min = j_min + 1/2 = start_j + 1/2
        #
        # zeta(-1) = sum_{u >= u_min, step 1} 2u * [u^2 - 1/4 - c^2]
        #          = 2 * [sum u^3 - (1/4 + c^2) * sum u]
        #          = 2 * [zeta_H(-3, u_min) - (1/4 + c^2) * zeta_H(-1, u_min)]
        #
        # where zeta_H(-n, a) = -B_{n+1}(a)/(n+1)

        u_min = start_j + 0.5
        coeff = 0.25 + c_sq  # 1/4 + (qn/2)^2

        zh_m3 = _hurwitz_zeta_neg_int(3, u_min)  # zeta_H(-3, u_min)
        zh_m1 = _hurwitz_zeta_neg_int(1, u_min)  # zeta_H(-1, u_min)

        value = 2.0 * (zh_m3 - coeff * zh_m1)

        # Neutral comparison: u_min = 1.5 (j starts from 1), c_sq = 0
        zh_m3_neutral = _hurwitz_zeta_neg_int(3, 1.5)
        zh_m1_neutral = _hurwitz_zeta_neg_int(1, 1.5)
        comparison = 2.0 * (zh_m3_neutral - 0.25 * zh_m1_neutral)

        return {
            "s": s,
            "charge_q": charge_q,
            "n_monopole": n_monopole,
            "j_min": j_min,
            "value": value,
            "method": "hurwitz_analytic_continuation",
            "converged": True,
            "comparison_neutral": comparison,
        }

    # General s <= 1: use direct partial sum (not analytically continued)
    total = 0.0
    j = start_j
    for _ in range(l_max):
        eig = j * (j + 1) - c_sq
        if eig > 0:
            if s == 0:
                total += (2 * j + 1)
            else:
                total += (2 * j + 1) * eig ** (-s)
        j += 1.0

    return {
        "s": s,
        "charge_q": charge_q,
        "n_monopole": n_monopole,
        "j_min": j_min,
        "value": total,
        "method": "direct_partial_sum",
        "converged": False,
        "comparison_neutral": None,
    }


# ---------------------------------------------------------------------------
# 3. Coleman-Weinberg effective potential from charged KK tower
# ---------------------------------------------------------------------------

def compute_charged_cw_potential(n_scalars, n_fermions, n_vectors,
                                  charge_q=1, n_monopole=1,
                                  R=None, l_max=100, mu_R_sq=1.0):
    """
    Compute the one-loop Coleman-Weinberg effective potential from charged
    matter fields on S^2 with monopole background.

    V_CW = sum_j (2j+1) * m_j^4 * [ln(m_j^2/mu^2) - c_s] / (64*pi^2)

    with boson/fermion sign and c_s = 3/2 (scalar), 5/6 (vector), 3/2 (fermion).

    The masses are m_j^2 = [j(j+1) - (qn/2)^2] / R^2.

    Parameters
    ----------
    n_scalars : int
        Number of real scalar DOF.
    n_fermions : int
        Number of Weyl fermion DOF.
    n_vectors : int
        Number of vector DOF.
    charge_q : int or float
        U(1) charge under the monopole gauge field (default 1).
    n_monopole : int
        Monopole number (default 1).
    R : float or None
        Physical radius in Planck units.  If None, uses gauge-matched value.
    l_max : int
        KK level cutoff.
    mu_R_sq : float
        Renormalization scale ratio mu^2*R^2 (default 1.0).

    Returns
    -------
    dict with keys:
        V_scalars : float -- scalar contribution
        V_fermions : float -- fermion contribution (negative sign)
        V_vectors : float -- vector contribution
        V_total : float -- sum
        V_total_planck : float -- in M_Pl^4 units
        n_scalars, n_fermions, n_vectors : int
        charge_q, n_monopole : int
        R : float
        j_min : float
        n_modes_summed : int
        effective_lambda6 : float or None -- V_CW interpreted as Lambda_6 * R^2
        assessment : str
    """
    if R is None:
        R = _R_GAUGE_MATCHED

    qn = abs(charge_q * n_monopole)
    j_min = qn / 2.0
    c_sq = (qn / 2.0) ** 2
    R2 = R * R
    R4 = R2 * R2
    prefactor = 1.0 / (64.0 * math.pi ** 2)

    # Scheme-dependent constants
    c_s_scalar = 3.0 / 2.0
    c_s_fermion = 3.0 / 2.0
    c_s_vector = 5.0 / 6.0

    start_j = j_min if j_min > 0 else 1.0

    V_scalars = 0.0
    V_fermions = 0.0
    V_vectors = 0.0
    n_modes = 0

    j = start_j
    for _ in range(l_max):
        eig = j * (j + 1) - c_sq
        if eig <= 0:
            j += 1.0
            continue

        m_sq = eig / R2
        m4 = m_sq * m_sq
        deg = 2 * j + 1
        ln_term_s = math.log(m_sq / mu_R_sq) - c_s_scalar if m_sq > 0 else 0
        ln_term_f = math.log(m_sq / mu_R_sq) - c_s_fermion if m_sq > 0 else 0
        ln_term_v = math.log(m_sq / mu_R_sq) - c_s_vector if m_sq > 0 else 0

        V_scalars += deg * m4 * ln_term_s * prefactor
        V_fermions -= deg * m4 * ln_term_f * prefactor  # fermion sign
        V_vectors += deg * m4 * ln_term_v * prefactor

        n_modes += int(deg)
        j += 1.0

    V_s_total = n_scalars * V_scalars
    V_f_total = n_fermions * V_fermions
    V_v_total = n_vectors * V_vectors
    V_total = V_s_total + V_f_total + V_v_total

    # Interpret V_CW as effective Lambda_6 * R^2 contribution
    # V_CW acts as an R-dependent correction to the potential
    # If V_CW ~ Lambda_eff * R^2 at the gauge-matched R, then
    # Lambda_eff = V_CW / R^2
    effective_lambda6 = V_total / R2 if R > 0 else None

    assessment = (
        f"Charged CW potential from {n_scalars}S + {n_fermions}F + {n_vectors}V "
        f"with charge q={charge_q} in monopole n={n_monopole}.  "
        f"j_min = {j_min}, {n_modes} modes summed.  "
        f"V_total = {V_total:.6e} M_Pl^4.  "
    )
    if effective_lambda6 is not None:
        assessment += f"Effective Lambda_6 contribution: {effective_lambda6:.6e}.  "
        if abs(effective_lambda6) > 1.0:
            assessment += "This is O(1) or larger -- potentially significant.  "
        else:
            assessment += "This is << 1 -- negligible compared to Lambda_6 = 14.2.  "

    return {
        "V_scalars": V_s_total,
        "V_fermions": V_f_total,
        "V_vectors": V_v_total,
        "V_total": V_total,
        "V_total_planck": V_total,
        "n_scalars": n_scalars,
        "n_fermions": n_fermions,
        "n_vectors": n_vectors,
        "charge_q": charge_q,
        "n_monopole": n_monopole,
        "R": R,
        "j_min": j_min,
        "n_modes_summed": n_modes,
        "effective_lambda6": effective_lambda6,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 4. One-loop correction to G from charged matter
# ---------------------------------------------------------------------------

def compute_charged_oneloop_g_correction(n_scalars, n_fermions, n_vectors,
                                          charge_q=1, n_monopole=1,
                                          R=None, constants=None):
    """
    Compute the one-loop correction to Newton's constant from charged
    matter KK modes in monopole background.

    This extends oneloop_g_correction.py (pure gravity) to include
    charged matter.  The correction formula is:

        delta(1/16piG) = sum_s c_s * N_s * zeta_charged(-1) / (16*pi^2*R^2)

    where zeta_charged(-1) is the monopole-modified spectral zeta.

    Parameters
    ----------
    n_scalars : int
    n_fermions : int
    n_vectors : int
    charge_q : int or float
    n_monopole : int
    R : float or None
    constants : SimpleNamespace or None

    Returns
    -------
    dict with keys:
        zeta_neutral : float -- standard zeta_{S^2}(-1) = -17/480
        zeta_charged : float -- monopole-modified zeta
        delta_zeta : float -- difference
        delta_G_over_G_gravity : float -- pure gravity contribution
        delta_G_over_G_matter : float -- charged matter contribution
        delta_G_over_G_total : float -- combined
        target_3_alpha_sq : float -- target value 3*alpha^2
        ratio_to_target : float -- total / target
        n_scalars, n_fermions, n_vectors : int
        assessment : str
    """
    if R is None:
        R = _R_GAUGE_MATCHED

    alpha = 1.0 / 137.036
    if constants is not None:
        try:
            alpha = float(constants.alpha)
        except (AttributeError, TypeError):
            pass

    target = 3.0 * alpha ** 2

    # Pure gravity contribution (from oneloop_g_correction.py)
    zeta_neutral = -17.0 / 480.0
    c_grav_total = (
        (1.0 / 6.0) * 3  # scalars: c_0 * N_scalar
        + 1.0 * 3          # vectors: c_1 * N_vector
        + (5.0 / 6.0) * 1  # tensors: c_2 * N_tensor
    )
    sector_factor_grav = zeta_neutral / (16.0 * math.pi ** 2 * R ** 2)
    delta_G_gravity = c_grav_total * sector_factor_grav

    # Charged matter contribution
    zeta_result = compute_monopole_spectral_zeta(
        charge_q, n_monopole, s=-1, l_max=10000
    )
    zeta_charged = zeta_result["value"]
    delta_zeta = zeta_charged - zeta_neutral

    # Matter Seeley-DeWitt coefficients
    c_matter = (
        (1.0 / 6.0) * n_scalars
        + 1.0 * n_vectors
        - (1.0 / 12.0) * n_fermions  # fermions contribute with opposite sign
    )
    sector_factor_matter = zeta_charged / (16.0 * math.pi ** 2 * R ** 2)
    delta_G_matter = c_matter * sector_factor_matter

    delta_G_total = delta_G_gravity + delta_G_matter

    ratio = delta_G_total / target if target != 0 else 0.0

    assessment = (
        f"Pure gravity: delta(G)/G = {delta_G_gravity:.6e} "
        f"(zeta_neutral = {zeta_neutral:.6f}).  "
        f"Charged matter ({n_scalars}S + {n_fermions}F + {n_vectors}V, "
        f"q={charge_q}): delta(G)/G = {delta_G_matter:.6e} "
        f"(zeta_charged = {zeta_charged:.6e}).  "
        f"Total: {delta_G_total:.6e}.  "
        f"Target 3*alpha^2 = {target:.6e}.  "
        f"Ratio: {ratio:.4e}.  "
    )
    if abs(ratio - 1.0) < 0.01:
        assessment += "MATCHES target to < 1%."
    elif abs(ratio) > 0.1:
        assessment += f"Within an order of magnitude (ratio {ratio:.2e})."
    else:
        assessment += f"Far from target (ratio {ratio:.2e})."

    return {
        "zeta_neutral": zeta_neutral,
        "zeta_charged": zeta_charged,
        "delta_zeta": delta_zeta,
        "delta_G_over_G_gravity": delta_G_gravity,
        "delta_G_over_G_matter": delta_G_matter,
        "delta_G_over_G_total": delta_G_total,
        "target_3_alpha_sq": target,
        "ratio_to_target": ratio,
        "n_scalars": n_scalars,
        "n_fermions": n_fermions,
        "n_vectors": n_vectors,
        "charge_q": charge_q,
        "n_monopole": n_monopole,
        "R": R,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 5. Scan over anomaly-free matter content
# ---------------------------------------------------------------------------

def scan_anomaly_free_matter_loops(n_monopole=1, constants=None):
    """
    For each anomaly-free gauge group from anomaly_cancellation.py,
    compute the charged one-loop correction to G and the CW potential.

    6D multiplet field counting:
        1 hypermultiplet = 4 real scalars + 2 Weyl fermions
        1 vector multiplet = 1 gauge boson + 2 Weyl fermions + 1 scalar

    Parameters
    ----------
    n_monopole : int
        Monopole number (default 1).
    constants : SimpleNamespace or None

    Returns
    -------
    dict with keys:
        results : list of dict -- per-group results
        n_groups : int
        any_match_3alpha2 : bool
        any_significant_lambda6 : bool
        best_g_correction : dict or None
        best_lambda6 : dict or None
        assessment : str
    """
    try:
        from alpha_ladder_core.anomaly_cancellation import _ANOMALY_FREE_GROUPS
    except ImportError:
        return {
            "results": [],
            "n_groups": 0,
            "any_match_3alpha2": False,
            "any_significant_lambda6": False,
            "best_g_correction": None,
            "best_lambda6": None,
            "assessment": "FAILED: anomaly_cancellation module not available.",
        }

    results = []
    any_match = False
    any_sig_lambda6 = False
    best_g = None
    best_l6 = None

    for key, data in _ANOMALY_FREE_GROUPS.items():
        n_hyper = data.get("n_hyper")
        n_vector = data.get("n_vector")

        if n_hyper is None or n_vector is None:
            continue

        # 6D multiplet field counting
        n_scalars = 4 * n_hyper + n_vector  # hyper scalars + vector scalars
        n_fermions = 2 * n_hyper + 2 * n_vector  # hyper + vector fermions (Weyl)
        n_vectors_gauge = n_vector  # gauge bosons

        # G correction from charged matter
        g_corr = compute_charged_oneloop_g_correction(
            n_scalars=n_scalars,
            n_fermions=n_fermions,
            n_vectors=n_vectors_gauge,
            charge_q=1,
            n_monopole=n_monopole,
            constants=constants,
        )

        # CW potential from charged matter
        cw = compute_charged_cw_potential(
            n_scalars=n_scalars,
            n_fermions=n_fermions,
            n_vectors=n_vectors_gauge,
            charge_q=1,
            n_monopole=n_monopole,
        )

        entry = {
            "group": data["group"],
            "key": key,
            "n_hyper": n_hyper,
            "n_vector": n_vector,
            "n_scalars": n_scalars,
            "n_fermions": n_fermions,
            "n_vectors_gauge": n_vectors_gauge,
            "g_correction": g_corr,
            "cw_potential": cw,
        }
        results.append(entry)

        if abs(g_corr["ratio_to_target"] - 1.0) < 0.1:
            any_match = True

        if cw["effective_lambda6"] is not None and abs(cw["effective_lambda6"]) > 1.0:
            any_sig_lambda6 = True

        # Track best
        if best_g is None or abs(g_corr["ratio_to_target"] - 1.0) < abs(best_g["g_correction"]["ratio_to_target"] - 1.0):
            best_g = entry
        if best_l6 is None or (cw["effective_lambda6"] is not None and
                                abs(cw["effective_lambda6"] - 14.2) <
                                abs((best_l6["cw_potential"]["effective_lambda6"] or 0) - 14.2)):
            best_l6 = entry

    assessment = (
        f"Scanned {len(results)} anomaly-free groups with monopole n={n_monopole}.  "
    )
    if any_match:
        assessment += "At least one group matches 3*alpha^2 to < 10%.  "
    else:
        assessment += "No group matches 3*alpha^2.  "
    if any_sig_lambda6:
        assessment += "At least one group produces O(1) effective Lambda_6.  "
    else:
        assessment += "No group produces significant effective Lambda_6.  "

    return {
        "results": results,
        "n_groups": len(results),
        "any_match_3alpha2": any_match,
        "any_significant_lambda6": any_sig_lambda6,
        "best_g_correction": best_g,
        "best_lambda6": best_l6,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 6. Summary
# ---------------------------------------------------------------------------

def summarize_charged_matter_loops(constants=None):
    """
    Run the full charged matter loop analysis and return a summary.

    This is the main entry point.  It computes:
        1. Monopole KK spectrum for q=1
        2. Monopole spectral zeta at s=-1
        3. G correction from E8xE8 matter content
        4. CW potential from E8xE8 matter content
        5. Scan over all anomaly-free groups
        6. Comparison with pure gravity result

    Parameters
    ----------
    constants : SimpleNamespace or None

    Returns
    -------
    dict with keys:
        monopole_spectrum : dict
        spectral_zeta : dict
        e8_g_correction : dict
        e8_cw_potential : dict
        group_scan : dict
        pure_gravity_zeta : float
        monopole_zeta : float
        zeta_shift : float
        gap2_impact : str
        honest_assessment : str
    """
    # 1. Monopole spectrum for q=1, n=1
    spectrum = compute_monopole_kk_spectrum(charge_q=1, n_monopole=1, l_max=20)

    # 2. Spectral zeta at s=-1
    zeta_result = compute_monopole_spectral_zeta(charge_q=1, n_monopole=1, s=-1)

    # 3. E8xE8 matter content: 740 hypers + 496 vectors
    n_scalars_e8 = 4 * 740 + 496  # = 3456
    n_fermions_e8 = 2 * 740 + 2 * 496  # = 2472
    n_vectors_e8 = 496

    e8_g_corr = compute_charged_oneloop_g_correction(
        n_scalars=n_scalars_e8,
        n_fermions=n_fermions_e8,
        n_vectors=n_vectors_e8,
        charge_q=1,
        n_monopole=1,
        constants=constants,
    )

    e8_cw = compute_charged_cw_potential(
        n_scalars=n_scalars_e8,
        n_fermions=n_fermions_e8,
        n_vectors=n_vectors_e8,
        charge_q=1,
        n_monopole=1,
    )

    # 4. Full scan over anomaly-free groups
    group_scan = scan_anomaly_free_matter_loops(n_monopole=1, constants=constants)

    # 5. Comparison
    zeta_neutral = -17.0 / 480.0
    zeta_charged = zeta_result["value"]
    zeta_shift = zeta_charged - zeta_neutral

    # Gap 2 impact: does the charged matter loop change the screening story?
    # The Salam-Sezgin dilaton is already Planck-mass, so screening is trivially
    # resolved.  But does the matter loop modify the dilaton mass?
    gap2_impact = (
        "The charged matter loop modifies the spectral zeta from "
        f"{zeta_neutral:.6f} (neutral) to {zeta_charged:.6e} (charged, q=1).  "
    )
    if abs(zeta_shift) < 1e-3:
        gap2_impact += (
            "The shift is small.  The Salam-Sezgin dilaton remains Planck-mass "
            "and screening is still trivially resolved."
        )
    else:
        gap2_impact += (
            f"The shift delta_zeta = {zeta_shift:.6e} is significant.  "
            "This could modify the effective potential and dilaton mass."
        )

    # Honest assessment
    honest = (
        "Extending the one-loop calculation to include charged matter fields "
        "in the Dirac monopole background on S^2.  "
        f"The monopole (n=1) shifts the KK spectrum: j starts from j_min = 1/2 "
        f"(for unit charge) instead of l=1, modifying the spectral zeta.  "
        f"Neutral zeta_{'{'}S^2{'}'}(-1) = {zeta_neutral:.6f}.  "
        f"Charged zeta = {zeta_charged:.6e}.  "
    )
    if group_scan["any_match_3alpha2"]:
        honest += (
            "At least one anomaly-free group produces a G correction "
            "matching 3*alpha^2 -- this would close the loop on the "
            "bridge coefficient derivation.  "
        )
    else:
        honest += (
            "No anomaly-free group produces a G correction matching "
            "3*alpha^2.  The bridge coefficient remains empirical.  "
        )
    if group_scan["any_significant_lambda6"]:
        honest += (
            "The CW potential from charged matter produces an O(1) "
            "effective Lambda_6 contribution.  "
        )
    else:
        honest += (
            "The CW potential from charged matter is too small to "
            "account for Lambda_6 = 14.2.  "
        )

    return {
        "monopole_spectrum": spectrum,
        "spectral_zeta": zeta_result,
        "e8_g_correction": e8_g_corr,
        "e8_cw_potential": e8_cw,
        "group_scan": group_scan,
        "pure_gravity_zeta": zeta_neutral,
        "monopole_zeta": zeta_charged,
        "zeta_shift": zeta_shift,
        "gap2_impact": gap2_impact,
        "honest_assessment": honest,
    }
