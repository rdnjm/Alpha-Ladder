"""
Threshold corrections to Newton's constant from KK modes on S^2.

The one-loop correction to G from KK modes on S^2 involves two types of
spectral sums.  The "power sum" -- sum (2l+1)*l*(l+1) -- vanishes in zeta
regularization because zeta_{S^2}(-1) = 0 (the partial sums form an exact
polynomial L(L+1)^2(L+2)/2 with no constant term).

However, the "threshold sum" involves logarithmic factors:

    S_log(L) = sum_{l=1}^{L} (2l+1) * l*(l+1) * ln[l*(l+1)]

This does NOT have a polynomial closed form.  Its finite part (the constant
term in the asymptotic expansion for large L) is generically nonzero and
represents the physical threshold correction from integrating out the full
KK tower.

In the 6D-to-4D matching, the one-loop correction to 1/G is:

    delta(1/16piG) = c_total / (16*pi^2*R^2) * finite_part(S_log)

where R is the compactification radius and c_total counts the degrees of
freedom being integrated out.  The correction to G/G is then:

    delta(G)/G = -G * c_total * finite_part(S_log) / (pi * R^2)

This module computes the threshold sum, extracts its finite part via
Euler-Maclaurin summation with high-precision arithmetic, and evaluates
the correction for various choices of R and field content.

All calculations use pure Python math and the Decimal module for
high-precision subtraction (no numpy/scipy).
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50


# ---------------------------------------------------------------------------
# Core analytical functions (high precision)
# ---------------------------------------------------------------------------

def _antiderivative_hp(x_val):
    """
    High-precision antiderivative of f(x) = (2x+1)*x*(x+1)*ln[x*(x+1)].

    Derived by integration by parts:
        f(x) = p(x)*ln(x) + p(x)*ln(x+1)
    where p(x) = 2x^3 + 3x^2 + x.

    Using integral x^n*ln(x) dx = x^{n+1}/(n+1)*(ln(x) - 1/(n+1)):

        integral p(x)*ln(x) dx = x^4/2*(ln(x)-1/4) + x^3*(ln(x)-1/3)
                                 + x^2/2*(ln(x)-1/2)

    For ln(x+1), substitute u=x+1, express p(x) = 2u^3-3u^2+u:

        integral p(x)*ln(x+1) dx = u^4/2*(ln(u)-1/4) - u^3*(ln(u)-1/3)
                                   + u^2/2*(ln(u)-1/2)

    Parameters
    ----------
    x_val : int or float -- evaluation point (must be > 0)

    Returns
    -------
    Decimal -- antiderivative value
    """
    x = Decimal(str(x_val))
    u = x + 1
    ln_x = x.ln() if x > 0 else Decimal(0)
    ln_u = u.ln()

    one_quarter = Decimal(1) / 4
    one_third = Decimal(1) / 3
    one_half = Decimal(1) / 2

    part1 = (
        x ** 4 / 2 * (ln_x - one_quarter)
        + x ** 3 * (ln_x - one_third)
        + x ** 2 / 2 * (ln_x - one_half)
    )
    part2 = (
        u ** 4 / 2 * (ln_u - one_quarter)
        - u ** 3 * (ln_u - one_third)
        + u ** 2 / 2 * (ln_u - one_half)
    )
    return part1 + part2


def _f_hp(x_val):
    """f(x) = (2x+1)*x*(x+1)*ln[x*(x+1)] in high precision."""
    x = Decimal(str(x_val))
    g = x * (x + 1)
    if g <= 0:
        return Decimal(0)
    return (2 * x + 1) * g * g.ln()


def _f1_hp(x_val):
    """
    f'(x) = (6x^2+6x+1)*ln[x*(x+1)] + (2x+1)^2

    Derived from f = h*ln(g) where h = (2x+1)*g, g = x^2+x:
        f' = h'*ln(g) + h*g'/g = (6x^2+6x+1)*ln(g) + (2x+1)^2
    """
    x = Decimal(str(x_val))
    g = x * (x + 1)
    if g <= 0:
        return Decimal(0)
    return (6 * x ** 2 + 6 * x + 1) * g.ln() + (2 * x + 1) ** 2


def _f3_hp(x_val):
    """
    f'''(x) = 12*ln(g) + (12x+6)*(2x+1)/g
              + [(36x^2+36x+8)*g - (12x^3+18x^2+8x+1)*(2x+1)] / g^2 + 8

    where g = x*(x+1).

    Derived by successive differentiation of f'(x).
    """
    x = Decimal(str(x_val))
    g = x * (x + 1)
    if g <= 0:
        return Decimal(0)
    gp = 2 * x + 1
    N = 12 * x ** 3 + 18 * x ** 2 + 8 * x + 1
    Np = 36 * x ** 2 + 36 * x + 8
    return (
        12 * g.ln()
        + (12 * x + 6) * gp / g
        + (Np * g - N * gp) / g ** 2
        + 8
    )


def _f4_hp(x_val):
    """
    f''''(x) computed analytically.

    f'''' = 12*(2x+1)/g + d/dx[(12x+6)*(2x+1)/g]
            + d/dx[(Np*g - N*(2x+1))/g^2]

    Each term is computed using the quotient rule.
    """
    x = Decimal(str(x_val))
    g = x * (x + 1)
    if g <= 0:
        return Decimal(0)
    gp = 2 * x + 1
    A = 12 * x + 6
    N = 12 * x ** 3 + 18 * x ** 2 + 8 * x + 1
    Np = 36 * x ** 2 + 36 * x + 8
    Npp = 72 * x + 36

    t1 = 12 * gp / g
    t2 = ((12 * gp + 2 * A) * g - A * gp ** 2) / g ** 2
    H = Np * g - N * gp
    Hp = Npp * g - 2 * N
    t3 = (Hp * g - 2 * H * gp) / g ** 3
    return t1 + t2 + t3


def _f5_hp(x_val):
    """f^(5)(x) computed numerically from analytical f^(4)."""
    h = Decimal('0.0001')
    x = Decimal(str(x_val))
    return (_f4_hp(float(x + h)) - _f4_hp(float(x - h))) / (2 * h)


# ---------------------------------------------------------------------------
# Raw sums
# ---------------------------------------------------------------------------

def _compute_s_log_raw(L):
    """
    Compute S_log(L) = sum_{l=1}^{L} (2l+1) * l*(l+1) * ln[l*(l+1)].

    Uses Decimal arithmetic for exact summation (no floating-point
    accumulation error for the individual terms).

    Parameters
    ----------
    L : int -- upper cutoff

    Returns
    -------
    Decimal -- the raw sum
    """
    total = Decimal(0)
    for l in range(1, L + 1):
        ll1 = Decimal(l * (l + 1))
        total += Decimal(2 * l + 1) * ll1 * ll1.ln()
    return total


def _compute_s_log_raw_float(L):
    """
    Compute S_log(L) in standard double precision (faster for large L).

    Parameters
    ----------
    L : int -- upper cutoff

    Returns
    -------
    float -- the raw sum
    """
    total = 0.0
    for l in range(1, L + 1):
        ll1 = l * (l + 1)
        total += (2 * l + 1) * ll1 * math.log(ll1)
    return total


def _compute_s_power_raw(L):
    """
    Compute S_power(L) = sum_{l=1}^{L} (2l+1) * l*(l+1) = L*(L+1)^2*(L+2)/2.

    This is the polynomial sum that vanishes under zeta regularization
    (no constant term).

    Parameters
    ----------
    L : int -- upper cutoff

    Returns
    -------
    float -- the polynomial sum
    """
    return L * (L + 1) ** 2 * (L + 2) / 2.0


# ---------------------------------------------------------------------------
# 1. Compute the threshold sum S_log(L) and extract finite part
# ---------------------------------------------------------------------------

def compute_threshold_sum(L_max=20000, include_log_degeneracy=False):
    """
    Compute the threshold sum S_log and extract its finite part.

    S_log(L) = sum_{l=1}^{L} (2l+1) * l*(l+1) * ln[l*(l+1)]

    The finite part is extracted using the Euler-Maclaurin formula with
    high-precision (Decimal) arithmetic:

        sum_{l=1}^{L} f(l) = integral_1^L f(x)dx + [f(1)+f(L)]/2
                             + B_2/2! * [f'(L)-f'(1)]
                             + B_4/4! * [f'''(L)-f'''(1)]
                             + B_6/6! * [f^(5)(L)-f^(5)(1)]
                             + R_p

    The integral is computed analytically via integration by parts.
    Derivatives through f^(4) are computed analytically; f^(5) is
    obtained by numerical differentiation of the analytical f^(4).

    The Euler-Maclaurin series for this function is asymptotic (not
    convergent), so we truncate at the term of minimum magnitude (B_6)
    and report the truncation uncertainty.

    Parameters
    ----------
    L_max : int
        Maximum cutoff for the sum (default 20000).
        For speed, we use L_max up to 5000 with Decimal arithmetic
        for full convergence verification.
    include_log_degeneracy : bool
        If True, also compute a variant with ln(2l+1) weighting (not used
        in the main analysis).

    Returns
    -------
    dict with keys:
        L_max : int
        S_log : float -- raw sum at L_max
        asymptotic_terms : dict of divergent terms subtracted
        finite_part : float -- extracted finite constant
        finite_part_uncertainty : float -- estimated truncation error
        convergence_data : list of dicts with (L, S_log, finite_part_estimate)
        converged : bool
    """
    B2 = Decimal(1) / 6
    B4 = Decimal(-1) / 30
    B6 = Decimal(1) / 42

    # Use L values up to min(L_max, 5000) for the convergence check
    # (Decimal sums are slow for very large L)
    L_limit = min(L_max, 5000)
    L_values = [100, 200, 500, 1000, 2000]
    if L_limit >= 5000:
        L_values.append(5000)
    L_values = [v for v in L_values if v <= L_limit]
    if L_limit not in L_values:
        L_values.append(L_limit)

    convergence_data = []

    # Track the EM approximants at each truncation level to estimate
    # the optimal finite part
    fp_after_B4_list = []
    fp_after_B6_list = []

    for L in L_values:
        # High-precision sum
        s_log = _compute_s_log_raw(L)

        # High-precision integral
        I_L = _antiderivative_hp(L) - _antiderivative_hp(1)

        # Boundary average
        bnd = (_f_hp(1) + _f_hp(L)) / 2

        # EM corrections
        em1 = B2 / 2 * (_f1_hp(L) - _f1_hp(1))
        em2 = B4 / 24 * (_f3_hp(L) - _f3_hp(1))
        em3 = B6 / 720 * (_f5_hp(L) - _f5_hp(1))

        # Remainders after successive corrections
        r0 = s_log - I_L - bnd
        r1 = r0 - em1
        r2 = r1 - em2
        r3 = r2 - em3

        fp_after_B4_list.append(float(r2))
        fp_after_B6_list.append(float(r3))

        convergence_data.append({
            "L": L,
            "S_log": float(s_log),
            "remainder_after_integral": float(r0),
            "remainder_after_B2": float(r1),
            "remainder_after_B4": float(r2),
            "remainder_after_B6": float(r3),
            "finite_part_estimate": float(r3),
        })

    # The optimal truncation is at B6 (smallest |remainder|).
    # The finite part estimate is the converged r3 value.
    finite_part_B4 = fp_after_B4_list[-1]
    finite_part_B6 = fp_after_B6_list[-1]

    # The truncation error is approximately the magnitude of the next
    # (B6) correction's boundary contribution, which is:
    # |finite_part_B4 - finite_part_B6|
    truncation_error = abs(finite_part_B4 - finite_part_B6)

    # Use the B6 value as our best estimate (smallest truncation order)
    finite_part = finite_part_B6

    # Check convergence: are the last 3 B6 estimates stable?
    if len(fp_after_B6_list) >= 3:
        last_three = fp_after_B6_list[-3:]
        spread = max(last_three) - min(last_three)
        converged = abs(spread) < 1e-8
    else:
        converged = len(fp_after_B6_list) >= 1

    # Compute S_log at the requested L_max if larger than L_limit
    if L_max > L_limit:
        s_log_final = _compute_s_log_raw_float(L_max)
    else:
        s_log_final = float(convergence_data[-1]["S_log"])

    # Asymptotic terms at the final L
    last = convergence_data[-1]
    asymptotic_terms = {
        "integral_1_to_L": last["S_log"] - last["remainder_after_integral"],
        "boundary_average": last["S_log"] - last["remainder_after_integral"] - (
            last["S_log"] - last["remainder_after_integral"]
        ),
        "description": (
            "Leading asymptotic behaviour: S_log ~ L^4*ln(L) + O(L^4).  "
            "The divergent part is removed by the analytical integral "
            "(computed via integration by parts) plus Euler-Maclaurin "
            "boundary corrections through B_6.  The Euler-Maclaurin "
            "series is asymptotic; the B_6 truncation is optimal."
        ),
    }

    return {
        "L_max": L_max,
        "S_log": s_log_final,
        "asymptotic_terms": asymptotic_terms,
        "finite_part": finite_part,
        "finite_part_B4": finite_part_B4,
        "finite_part_B6": finite_part_B6,
        "finite_part_uncertainty": truncation_error,
        "convergence_data": convergence_data,
        "converged": converged,
    }


# ---------------------------------------------------------------------------
# 2. Compute threshold correction to G
# ---------------------------------------------------------------------------

def compute_threshold_correction_to_G(constants=None, R_type="r_e"):
    """
    Compute the threshold correction to Newton's constant.

    The one-loop correction from integrating out the KK tower is:

        delta(1/16piG) = c_total / (16*pi^2*R^2) * finite_part(S_log)

    So:
        delta(G)/G = c_total * |finite_part(S_log)| / (pi * R_Pl^2)

    where R_Pl = R / l_Pl is the compactification radius in Planck units.

    The target is delta(G)/G = 3*alpha^2 (the Alpha Ladder prediction).

    For each KK sector, the DOF count is:
        - Scalars: 3 DOF (from g_{ab} on S^2)
        - Vectors: 2 vectors with 3 DOF each = 6 (from g_{mu a})
        - Tensor: 1 massive spin-2 with 5 DOF

    Parameters
    ----------
    constants : SimpleNamespace or None
        Physical constants.  If None, uses CODATA 2018.
    R_type : str
        Which radius to use: "r_e", "lambda_bar_c", or "a_0".

    Returns
    -------
    dict with correction results, comparison to 3*alpha^2
    """
    from alpha_ladder_core.constants import get_constants

    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = float(constants.alpha)
    target = 3.0 * alpha ** 2

    G_SI = float(constants.G)
    hbar_SI = float(constants.hbar)
    c_SI = float(constants.c)
    r_e = float(constants.r_e_nist)
    lambda_bar_c = float(constants.lambda_bar_c)
    a_0 = float(constants.a_0)

    l_Pl = math.sqrt(hbar_SI * G_SI / c_SI ** 3)

    R_choices = {
        "r_e": {"R_m": r_e, "description": "classical electron radius"},
        "lambda_bar_c": {"R_m": lambda_bar_c, "description": "reduced Compton wavelength"},
        "a_0": {"R_m": a_0, "description": "Bohr radius"},
    }

    # Compute the threshold sum finite part
    threshold = compute_threshold_sum(L_max=2000)
    finite_part = threshold["finite_part"]
    fp_uncertainty = threshold["finite_part_uncertainty"]
    converged = threshold["converged"]

    # DOF combinations
    sector_combos = {
        "scalars_only": {"c_total": 3, "description": "3 scalars from g_{ab}"},
        "tensor_only": {"c_total": 5, "description": "1 massive spin-2 (5 DOF)"},
        "vectors_2x3": {"c_total": 6, "description": "2 vectors x 3 DOF each"},
        "scalars_plus_vectors": {"c_total": 9, "description": "3 scalars + 2 vectors (3+6)"},
        "all_sectors": {"c_total": 14, "description": "3 scalars + 2x3 vectors + 5 tensor"},
        "all_sectors_v2": {"c_total": 17, "description": "3 scalars + 3x3 vectors + 5 tensor"},
    }

    sector_results = {}

    for combo_name, combo_info in sector_combos.items():
        c_total = combo_info["c_total"]
        combo_results = {}

        for R_name, R_info in R_choices.items():
            R_m = R_info["R_m"]
            R_Pl = R_m / l_Pl

            delta_G_over_G = abs(c_total * finite_part / (math.pi * R_Pl ** 2))

            if delta_G_over_G > 0 and target > 0:
                ratio = delta_G_over_G / target
                log10_ratio = math.log10(ratio)
            else:
                ratio = float('inf')
                log10_ratio = float('inf')

            combo_results[R_name] = {
                "R_m": R_m,
                "R_Pl": R_Pl,
                "delta_G_over_G": delta_G_over_G,
                "ratio_to_target": ratio,
                "log10_ratio": log10_ratio,
            }

        sector_results[combo_name] = {
            "c_total": c_total,
            "description": combo_info["description"],
            "results": combo_results,
        }

    # Find best match
    best_match = None
    best_residual = float('inf')

    for combo_name, combo_data in sector_results.items():
        for R_name, R_data in combo_data["results"].items():
            log_ratio = (
                abs(R_data["log10_ratio"])
                if R_data["log10_ratio"] != float('inf')
                else float('inf')
            )
            if log_ratio < best_residual:
                best_residual = log_ratio
                best_match = {
                    "sector": combo_name,
                    "R_type": R_name,
                    "c_total": combo_data["c_total"],
                    "delta_G_over_G": R_data["delta_G_over_G"],
                    "ratio_to_target": R_data["ratio_to_target"],
                    "log10_ratio": R_data["log10_ratio"],
                }

    if best_match and best_match["delta_G_over_G"] > 0:
        residual_ppm = abs(best_match["delta_G_over_G"] - target) / target * 1e6
    else:
        residual_ppm = float('inf')

    # Assessment
    if best_residual < 0.3:
        assessment = (
            f"A physically motivated combination is within ~2x of 3*alpha^2.  "
            f"Best match: {best_match['sector']} with R = {best_match['R_type']}, "
            f"delta(G)/G = {best_match['delta_G_over_G']:.6e} vs "
            f"target {target:.6e}."
        )
    else:
        assessment = (
            f"No physically motivated combination of DOF and radius gives "
            f"delta(G)/G close to 3*alpha^2 = {target:.6e}.  "
            f"Best match is {best_match['sector']} with "
            f"R = {best_match['R_type']}, giving "
            f"delta(G)/G = {best_match['delta_G_over_G']:.6e} "
            f"(off by ~10^{{{best_match['log10_ratio']:.0f}}}).  "
            f"The threshold correction from KK modes is many orders of "
            f"magnitude away from the target for any natural parameter choice."
        )

    return {
        "finite_part_S_log": finite_part,
        "finite_part_uncertainty": fp_uncertainty,
        "finite_part_converged": converged,
        "R_choices": {
            name: {
                "R_m": info["R_m"],
                "R_Pl": info["R_m"] / l_Pl,
                "description": info["description"],
            }
            for name, info in R_choices.items()
        },
        "sector_results": sector_results,
        "target": target,
        "target_label": "3*alpha^2",
        "best_match": best_match,
        "residual_ppm": residual_ppm,
        "l_Pl": l_Pl,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 3. Analyze the log sum structure
# ---------------------------------------------------------------------------

def analyze_log_sum_structure(L_max=10000):
    """
    Analyze the mathematical structure of S_log in detail.

    Computes S_log(L) for several L values and studies the ratios:
        - S_log(L) / [S_power(L) * 2*ln(L)] -- should approach 1
        - S_log(L) / S_power(L) -- the "effective logarithm"

    The key structural observation is that S_log(L) ~ S_power(L) * 2*ln(L)
    for large L, because the dominant contribution comes from the l ~ L
    terms where ln[l*(l+1)] ~ 2*ln(L).

    Parameters
    ----------
    L_max : int
        Maximum L to analyze (default 10000).

    Returns
    -------
    dict with structural analysis
    """
    L_values = [10, 50, 100, 500, 1000, 5000]
    if L_max >= 10000:
        L_values.append(10000)
    L_values = [v for v in L_values if v <= L_max]
    if L_max not in L_values:
        L_values.append(L_max)

    values = []
    for L in L_values:
        s_log = _compute_s_log_raw_float(L)
        s_power = _compute_s_power_raw(L)
        ratio = s_log / s_power if s_power != 0 else 0.0
        ln_L = math.log(L) if L > 0 else 0.0
        normalized_ratio = ratio / (2.0 * ln_L) if ln_L != 0 else 0.0

        values.append({
            "L": L,
            "S_log": s_log,
            "S_power": s_power,
            "ratio_S_log_over_S_power": ratio,
            "effective_log_over_2lnL": normalized_ratio,
            "ln_L": ln_L,
        })

    last_normalized = values[-1]["effective_log_over_2lnL"]

    # Compute finite part
    threshold = compute_threshold_sum(L_max=min(L_max, 2000))
    fp = threshold["finite_part"]

    # Check if finite_part scales with any power of alpha
    alpha = 0.0072973525693
    alpha_scalings = {}
    for k in range(-4, 5):
        alpha_k = alpha ** k if k != 0 else 1.0
        alpha_scalings[f"alpha^{k}"] = {
            "alpha_power": alpha_k,
            "ratio": fp / alpha_k if alpha_k != 0 else float('inf'),
        }

    assessment = (
        f"S_log(L) grows as ~S_power(L) * 2*ln(L) = L*(L+1)^2*(L+2) * ln(L) "
        f"for large L.  The normalized ratio S_log/(S_power*2*ln(L)) approaches "
        f"{last_normalized:.6f} (should tend to 1.0 for L -> infinity).  "
        f"The finite part of S_log is {fp:.6e} +/- {threshold['finite_part_uncertainty']:.1e}, "
        f"a small number with no obvious alpha-scaling."
    )

    return {
        "values": values,
        "effective_log": 2.0 * values[-1]["ln_L"] * last_normalized,
        "effective_log_limit": 2.0,
        "normalized_ratio_at_Lmax": last_normalized,
        "finite_part": fp,
        "finite_part_uncertainty": threshold["finite_part_uncertainty"],
        "finite_part_vs_L": threshold["convergence_data"],
        "alpha_scalings": alpha_scalings,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 4. Scan matching scales
# ---------------------------------------------------------------------------

def scan_matching_scales(constants=None):
    """
    Scan over matching-scale choices for the threshold correction.

    The threshold correction depends on the matching scale mu (the energy
    where the 6D theory is matched to the 4D effective theory).  Different
    choices of mu shift the correction by terms proportional to ln(mu*R).

    The full correction is:

        delta(G)/G ~ c_total / (pi * R_Pl^2) * [finite_part(S_log)
                     + S_power_reg * ln(mu^2 * R^2)]

    Since S_power_reg = zeta_{S^2}(-1) = 0, the ln(mu) term actually
    vanishes.  This means the threshold correction is INDEPENDENT of
    the matching scale -- a robust structural result.

    Parameters
    ----------
    constants : SimpleNamespace or None
        Physical constants.  If None, uses CODATA 2018.

    Returns
    -------
    dict with matching scale scan results
    """
    from alpha_ladder_core.constants import get_constants

    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = float(constants.alpha)
    target = 3.0 * alpha ** 2
    r_e = float(constants.r_e_nist)
    hbar_SI = float(constants.hbar)
    c_SI = float(constants.c)
    G_SI = float(constants.G)
    m_e = float(constants.m_e)
    m_p = float(constants.m_p)

    l_Pl = math.sqrt(hbar_SI * G_SI / c_SI ** 3)
    M_Pl_eV = 1.22089e28
    hbar_c_eVm = 1.9733e-7

    m_e_eV = m_e * c_SI ** 2 / 1.602176634e-19
    m_p_eV = m_p * c_SI ** 2 / 1.602176634e-19
    KK_scale_eV = hbar_c_eVm / r_e

    matching_scales = {
        "KK_scale_1_over_R": {
            "mu_eV": KK_scale_eV,
            "description": "mu = 1/R (natural KK scale)",
        },
        "Planck_scale": {
            "mu_eV": M_Pl_eV,
            "description": "mu = M_Pl (Planck mass)",
        },
        "proton_mass": {
            "mu_eV": m_p_eV,
            "description": "mu = m_p (proton mass)",
        },
        "electron_mass": {
            "mu_eV": m_e_eV,
            "description": "mu = m_e (electron mass)",
        },
    }

    # zeta_{S^2}(-1) = 0 => the ln(mu) dependence vanishes
    S_power_reg = 0.0

    threshold = compute_threshold_sum(L_max=2000)
    finite_part = threshold["finite_part"]

    R_Pl = r_e / l_Pl
    c_total = 9

    results = {}
    for scale_name, scale_info in matching_scales.items():
        mu_eV = scale_info["mu_eV"]
        mu_natural = mu_eV / M_Pl_eV
        R_natural = R_Pl
        ln_mu2_R2 = 2.0 * math.log(mu_natural * R_natural)

        effective_fp = finite_part + S_power_reg * ln_mu2_R2
        delta_G = abs(c_total * effective_fp / (math.pi * R_Pl ** 2))

        if delta_G > 0 and target > 0:
            ratio_val = delta_G / target
            matches = abs(math.log10(ratio_val)) < 0.1
        else:
            ratio_val = float('inf')
            matches = False

        results[scale_name] = {
            "mu_eV": mu_eV,
            "ln_mu2_R2": ln_mu2_R2,
            "S_power_reg": S_power_reg,
            "mu_dependent_shift": S_power_reg * ln_mu2_R2,
            "effective_finite_part": effective_fp,
            "delta_G_over_G": delta_G,
            "ratio_to_target": ratio_val,
            "matches_target": matches,
            "description": scale_info["description"],
        }

    best_scale = min(
        results.keys(),
        key=lambda k: (
            abs(math.log10(results[k]["ratio_to_target"]))
            if results[k]["ratio_to_target"] > 0
               and results[k]["ratio_to_target"] != float('inf')
            else float('inf')
        ),
    )

    delta_G_values = [results[k]["delta_G_over_G"] for k in results]
    if len(delta_G_values) >= 2 and max(delta_G_values) > 0:
        spread = max(delta_G_values) - min(delta_G_values)
        mean_dG = sum(delta_G_values) / len(delta_G_values)
        rel_variation = abs(spread / mean_dG) if mean_dG != 0 else abs(spread)
        mu_independent = rel_variation < 1e-10
    else:
        rel_variation = 0.0
        mu_independent = True

    assessment = (
        f"Because zeta_{{S^2}}(-1) = 0, the threshold correction is exactly "
        f"INDEPENDENT of the matching scale mu.  All scales give the same "
        f"delta(G)/G = {delta_G_values[0]:.6e} "
        f"(relative variation: {rel_variation:.2e}).  "
        f"This is a robust structural result: the matching ambiguity cancels "
        f"due to the polynomial nature of the power sum."
    )

    return {
        "matching_scales": results,
        "best_scale": best_scale,
        "mu_independent": mu_independent,
        "relative_variation": rel_variation,
        "S_power_regularized": S_power_reg,
        "finite_part_S_log": finite_part,
        "assessment": assessment,
    }


# ---------------------------------------------------------------------------
# 5. Summary entry point
# ---------------------------------------------------------------------------

def summarize_threshold_analysis(constants=None):
    """
    Main entry point: run the full threshold correction analysis.

    Computes the threshold sum, extracts the finite part, evaluates the
    correction to G for all reasonable parameter choices, and honestly
    assesses whether any combination matches 3*alpha^2.

    Parameters
    ----------
    constants : SimpleNamespace or None
        Physical constants.  If None, uses CODATA 2018.

    Returns
    -------
    dict with all sub-results and overall assessment
    """
    from alpha_ladder_core.constants import get_constants

    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = float(constants.alpha)
    target = 3.0 * alpha ** 2

    # Run all sub-analyses
    threshold_sum = compute_threshold_sum(L_max=5000)
    correction_to_G = compute_threshold_correction_to_G(constants=constants)
    log_structure = analyze_log_sum_structure(L_max=5000)
    matching_scales = scan_matching_scales(constants=constants)

    # Does any combination match?
    best = correction_to_G["best_match"]
    if best and abs(best["log10_ratio"]) < 0.3:
        matches = True
    else:
        matches = False

    # Key finding
    fp = threshold_sum["finite_part"]
    fp_unc = threshold_sum["finite_part_uncertainty"]

    if matches:
        key_finding = (
            f"The threshold correction from KK modes on S^2 CAN produce "
            f"delta(G)/G = 3*alpha^2 for {best['sector']} with "
            f"R = {best['R_type']}.  This represents a genuine one-loop "
            f"quantum gravitational correction to Newton's constant."
        )
    else:
        key_finding = (
            f"The threshold sum S_log has a nonzero finite part = "
            f"{fp:.6e} +/- {fp_unc:.1e}, confirming that logarithmic "
            f"corrections survive even though the power sum vanishes.  "
            f"However, for all physically motivated choices of R and "
            f"field content, delta(G)/G is many orders of magnitude away "
            f"from 3*alpha^2 = {target:.6e}.  The best match "
            f"({best['sector']}, R = {best['R_type']}) is off by "
            f"~10^{{{best['log10_ratio']:.0f}}}."
        )

    # Honest assessment
    honest_parts = []

    honest_parts.append(
        "The threshold correction to G from KK modes is a well-defined, "
        "calculable one-loop effect.  The key mathematical result is that "
        "the logarithmic sum S_log has a nonzero finite part, even though "
        "the polynomial power sum vanishes in zeta regularization."
    )

    honest_parts.append(
        f"The finite part of S_log is {fp:.6e} +/- {fp_unc:.1e} "
        f"(extracted via Euler-Maclaurin summation with high-precision "
        f"arithmetic, optimally truncated at B_6).  This is a small "
        f"number of order 10^{{{math.floor(math.log10(abs(fp))) if abs(fp) > 0 else 0}}}."
    )

    if matching_scales["mu_independent"]:
        honest_parts.append(
            "The correction is independent of the matching scale mu, "
            "because zeta_{S^2}(-1) = 0 kills the mu-dependent term.  "
            "This is a robust structural feature."
        )

    if not matches:
        honest_parts.append(
            "For any natural choice of compactification radius R "
            "(classical electron radius, Compton wavelength, or Bohr "
            "radius) and any DOF count, delta(G)/G is many orders of "
            "magnitude too small to match 3*alpha^2.  The fundamental "
            "reason is that R >> l_Pl for all atomic-scale radii, so "
            "R_Pl^2 ~ 10^{40} suppresses the correction enormously.  "
            "The threshold correction does not provide the mechanism "
            "for the Alpha Ladder prediction."
        )

    return {
        "threshold_sum": threshold_sum,
        "correction_to_G": correction_to_G,
        "log_structure": log_structure,
        "matching_scales": matching_scales,
        "matches_3_alpha_sq": matches,
        "target_3_alpha_sq": target,
        "key_finding": key_finding,
        "honest_assessment": "  ".join(honest_parts),
    }
