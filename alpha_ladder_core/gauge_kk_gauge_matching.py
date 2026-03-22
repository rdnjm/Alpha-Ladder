"""
KK Gauge Matching and Coleman-Weinberg Potential
=================================================

Tree-level matching of Kaluza-Klein gauge fields to electromagnetism
and the Coleman-Weinberg effective potential for the dilaton in the
Alpha Ladder S^2 compactification framework.

Physics basis: Dereli-Senikoglu (arXiv:2601.08443) 6D Einstein-Hilbert
action reduced on M_4 x S^2.  The coset SO(3)/SO(2) yields two physical
gauge fields.  Identifying one with the photon fixes the dilaton vev
through alpha_KK = alpha_EM.

Pure Python -- only ``import math``, no numpy/scipy.

Functions
---------
1. tree_level_matching      -- dilaton vev from alpha_KK = alpha_EM
2. kk_mass_spectrum         -- KK mode masses on S^2
3. coleman_weinberg_potential -- one-loop CW potential V(phi)
4. find_cw_minimum          -- locate CW minimum by bisection
5. scan_n_charged           -- sweep n_charged for coincidence
6. alpha_running            -- one-loop running of gauge coupling
7. consistency_with_G       -- check dilaton vev vs golden ratio
8. summarize_gauge_matching -- full report
"""

import math

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
ALPHA_EM = 1.0 / 137.035999084          # CODATA 2018
G_N = 6.674298e-11                       # m^3 kg^-1 s^-2
C_LIGHT = 2.99792458e8                   # m s^-1
HBAR = 1.054571817e-34                   # J s
M_E_KG = 9.1093837015e-31               # electron mass (kg)
M_PL_EV = 1.22089e28                     # Planck mass (eV)
PHI_GOLDEN = (1.0 + math.sqrt(5.0)) / 2.0  # golden ratio


# ===================================================================
# 1. tree_level_matching
# ===================================================================
def tree_level_matching(alpha_em=None):
    """Compute the dilaton vev that matches alpha_KK = alpha_EM.

    In the S^2 reduction the KK gauge coupling is g_KK = e^{2 phi},
    so alpha_KK = g_KK^2 / (4 pi) = e^{4 phi} / (4 pi).
    Setting alpha_KK = alpha_EM gives
        phi_vev = (1/4) ln(4 pi alpha_EM).

    Returns a dict with all derived quantities.
    """
    if alpha_em is None:
        alpha_em = ALPHA_EM

    phi_vev = 0.25 * math.log(4.0 * math.pi * alpha_em)

    e_to_phi = math.exp(phi_vev)
    e_to_2phi = math.exp(2.0 * phi_vev)
    e_to_4phi = math.exp(4.0 * phi_vev)

    g_kk = math.exp(2.0 * phi_vev)          # = e^{2 phi}
    alpha_kk = e_to_4phi / (4.0 * math.pi)  # should equal alpha_em

    verification = abs(alpha_kk - alpha_em) / alpha_em < 1.0e-12

    return {
        "phi_vev": phi_vev,
        "g_KK": g_kk,
        "alpha_KK": alpha_kk,
        "alpha_EM": alpha_em,
        "e_to_phi": e_to_phi,
        "e_to_2phi": e_to_2phi,
        "e_to_4phi": e_to_4phi,
        "verification": verification,
        "description": (
            "Tree-level matching of KK gauge coupling to alpha_EM.  "
            "phi_vev = (1/4) ln(4 pi alpha_EM) = {:.6f}.  "
            "alpha_KK at this vev = {:.10e}, alpha_EM = {:.10e}.  "
            "Match verified: {}.".format(
                phi_vev, alpha_kk, alpha_em, verification
            )
        ),
    }


# ===================================================================
# 2. kk_mass_spectrum
# ===================================================================
def kk_mass_spectrum(phi_vev, a_0=1.0, l_max=20):
    """KK mode masses on S^2 with dilaton vev.

    The physical radius is R_phys = a_0 * e^{phi_vev}.
    Mode masses: m_l = sqrt(l(l+1)) / R_phys,  l = 1 .. l_max.
    Degeneracy of each mode: 2l+1.
    """
    R_phys = a_0 * math.exp(phi_vev)
    m1 = math.sqrt(2.0) / R_phys  # l=1

    modes = []
    for l in range(1, l_max + 1):
        m_l = math.sqrt(l * (l + 1)) / R_phys
        modes.append({
            "l": l,
            "m_l": m_l,
            "m_l_over_m1": m_l / m1,
            "degeneracy": 2 * l + 1,
        })

    return {
        "modes": modes,
        "R_phys": R_phys,
        "phi_vev": phi_vev,
        "a_0": a_0,
    }


# ===================================================================
# helpers for Coleman-Weinberg
# ===================================================================
def _m_l_sq(l, phi, a_0):
    """KK mass squared for mode l at dilaton value phi."""
    return l * (l + 1) / (a_0 ** 2 * math.exp(2.0 * phi))


def _m_l_fourth(l, phi, a_0):
    """m_l^4."""
    return (l * (l + 1)) ** 2 / (a_0 ** 4 * math.exp(4.0 * phi))


# ===================================================================
# 3. coleman_weinberg_potential
# ===================================================================
def coleman_weinberg_potential(phi, a_0=1.0, n_charged=1, l_max=50,
                                mu_renorm=1.0):
    """One-loop Coleman-Weinberg potential at a single phi value.

    Contributions
    -------------
    * n_charged complex scalars  (2 n_charged real dof each):
        V_scalar = 2 n_charged / (64 pi^2) sum_{l=1}^{l_max}
                   (2l+1) m_l^4 [ln(m_l^2/mu^2) - 3/2]

    * 2 physical massive vectors (3 polarisations each):
        V_gauge  = 2 * 3 / (64 pi^2) sum_{l=1}^{l_max}
                   (2l+1) m_l^4 [ln(m_l^2/mu^2) - 5/6]

    The l=0 scalar mode is massless (constant on S^2); we start the
    sum at l=1 for both scalars and vectors.

    Returns V_CW, dV_dphi at the given phi.
    """
    prefactor = 1.0 / (64.0 * math.pi ** 2)
    mu_sq = mu_renorm ** 2

    V_total = 0.0
    dV_total = 0.0

    for l in range(1, l_max + 1):
        deg = 2 * l + 1
        m4 = _m_l_fourth(l, phi, a_0)
        m2 = _m_l_sq(l, phi, a_0)
        if m2 <= 0:
            continue
        ln_ratio = math.log(m2 / mu_sq)

        # scalar contribution
        V_s = 2.0 * n_charged * deg * m4 * (ln_ratio - 1.5)
        # d/dphi [m^4 (ln-C)] = m^4 (-4 ln + 4C - 2)
        dV_s = 2.0 * n_charged * deg * m4 * (-4.0 * ln_ratio + 4.0 * 1.5 - 2.0)

        # gauge contribution (2 vectors x 3 polarisations)
        V_g = 6.0 * deg * m4 * (ln_ratio - 5.0 / 6.0)
        dV_g = 6.0 * deg * m4 * (-4.0 * ln_ratio + 4.0 * 5.0 / 6.0 - 2.0)

        V_total += V_s + V_g
        dV_total += dV_s + dV_g

    V_total *= prefactor
    dV_total *= prefactor

    return {"V_CW": V_total, "dV_dphi": dV_total, "phi": phi}


def coleman_weinberg_scan(a_0=1.0, n_charged=1, l_max=50, mu_renorm=1.0,
                           phi_min=-3.0, phi_max=1.0, n_points=100):
    """Evaluate V_CW and dV/dphi on an evenly-spaced phi grid."""
    step = (phi_max - phi_min) / (n_points - 1)
    phi_grid = [phi_min + i * step for i in range(n_points)]
    V_list = []
    dV_list = []

    for p in phi_grid:
        res = coleman_weinberg_potential(p, a_0, n_charged, l_max, mu_renorm)
        V_list.append(res["V_CW"])
        dV_list.append(res["dV_dphi"])

    # detect minima: dV changes sign from negative to positive
    has_minimum = False
    phi_minimum = None
    V_at_minimum = None
    for i in range(len(dV_list) - 1):
        if dV_list[i] < 0 and dV_list[i + 1] > 0:
            # linear interpolation
            frac = -dV_list[i] / (dV_list[i + 1] - dV_list[i])
            phi_min_est = phi_grid[i] + frac * step
            V_min_est = V_list[i] + frac * (V_list[i + 1] - V_list[i])
            has_minimum = True
            phi_minimum = phi_min_est
            V_at_minimum = V_min_est
            break

    return {
        "phi_grid": phi_grid,
        "V_CW": V_list,
        "dV_dphi": dV_list,
        "has_minimum": has_minimum,
        "phi_minimum": phi_minimum,
        "V_at_minimum": V_at_minimum,
    }


# ===================================================================
# 4. find_cw_minimum
# ===================================================================
def find_cw_minimum(a_0=1.0, n_charged=1, l_max=50, mu_renorm=1.0):
    """Scan phi in [-3, 1] for a CW minimum; refine by bisection.

    Returns phi_min, comparison with tree-level vev, dilaton mass.
    """
    tree = tree_level_matching()
    phi_vev = tree["phi_vev"]

    # coarse scan
    scan = coleman_weinberg_scan(a_0, n_charged, l_max, mu_renorm,
                                  phi_min=-3.0, phi_max=1.0, n_points=200)

    if not scan["has_minimum"]:
        return {
            "phi_min": None,
            "phi_vev": phi_vev,
            "ratio": None,
            "m_phi_squared": None,
            "m_phi": None,
            "coincidence_with_alpha_matching": None,
            "description": "No Coleman-Weinberg minimum found in [-3, 1].",
            "honest_assessment": (
                "The CW potential from charged scalars + gauge fields "
                "does not develop a minimum in the scanned range.  "
                "The flat direction is not lifted toward a local minimum "
                "by one-loop effects with these field contents."
            ),
        }

    # bisection refinement on dV/dphi = 0
    phi_lo = scan["phi_minimum"] - 0.05
    phi_hi = scan["phi_minimum"] + 0.05

    dV_lo = coleman_weinberg_potential(phi_lo, a_0, n_charged, l_max,
                                       mu_renorm)["dV_dphi"]
    dV_hi = coleman_weinberg_potential(phi_hi, a_0, n_charged, l_max,
                                       mu_renorm)["dV_dphi"]

    # make sure we bracket a zero
    if dV_lo * dV_hi > 0:
        # widen
        phi_lo = scan["phi_minimum"] - 0.5
        phi_hi = scan["phi_minimum"] + 0.5
        dV_lo = coleman_weinberg_potential(phi_lo, a_0, n_charged, l_max,
                                           mu_renorm)["dV_dphi"]
        dV_hi = coleman_weinberg_potential(phi_hi, a_0, n_charged, l_max,
                                           mu_renorm)["dV_dphi"]

    if dV_lo * dV_hi > 0:
        # cannot bracket -- fall back to coarse estimate
        phi_min_refined = scan["phi_minimum"]
    else:
        for _ in range(80):
            phi_mid = 0.5 * (phi_lo + phi_hi)
            dV_mid = coleman_weinberg_potential(phi_mid, a_0, n_charged,
                                                l_max, mu_renorm)["dV_dphi"]
            if dV_mid == 0.0:
                break
            if dV_lo * dV_mid < 0:
                phi_hi = phi_mid
                dV_hi = dV_mid
            else:
                phi_lo = phi_mid
                dV_lo = dV_mid
        phi_min_refined = 0.5 * (phi_lo + phi_hi)

    # dilaton mass: d^2V/dphi^2 by finite difference
    eps = 1.0e-6
    dV_plus = coleman_weinberg_potential(phi_min_refined + eps, a_0,
                                         n_charged, l_max, mu_renorm)["dV_dphi"]
    dV_minus = coleman_weinberg_potential(phi_min_refined - eps, a_0,
                                          n_charged, l_max, mu_renorm)["dV_dphi"]
    m_phi_sq = (dV_plus - dV_minus) / (2.0 * eps)

    ratio = phi_min_refined / phi_vev if phi_vev != 0 else None
    dev_pct = abs(phi_min_refined - phi_vev) / abs(phi_vev) * 100.0 if phi_vev != 0 else None

    m_phi = math.sqrt(m_phi_sq) if m_phi_sq > 0 else None

    return {
        "phi_min": phi_min_refined,
        "phi_vev": phi_vev,
        "ratio": ratio,
        "m_phi_squared": m_phi_sq,
        "m_phi": m_phi,
        "coincidence_with_alpha_matching": dev_pct,
        "description": (
            "CW minimum at phi = {:.6f}.  Tree-level matched vev = {:.6f}.  "
            "Ratio = {:.4f}.  Deviation = {:.2f}%.  "
            "Dilaton mass^2 = {:.6e}.".format(
                phi_min_refined, phi_vev,
                ratio if ratio else 0.0,
                dev_pct if dev_pct else 0.0,
                m_phi_sq,
            )
        ),
        "honest_assessment": (
            "The CW minimum location depends on the renormalisation scale mu "
            "and the field content (n_charged).  A coincidence with the "
            "tree-level alpha-matching vev would be striking but is not "
            "guaranteed.  The CW potential is one-loop and scheme-dependent; "
            "a true dynamical determination of alpha would require a "
            "UV-complete calculation."
        ),
    }


# ===================================================================
# 5. scan_n_charged
# ===================================================================
def scan_n_charged(n_values=None, a_0=1.0, l_max=50):
    """Scan over n_charged values and find CW minimum for each."""
    if n_values is None:
        n_values = [1, 2, 3, 4, 5, 10, 20, 50, 100]

    tree = tree_level_matching()
    phi_vev = tree["phi_vev"]

    results = []
    for n in n_values:
        res = find_cw_minimum(a_0=a_0, n_charged=n, l_max=l_max)
        results.append({
            "n_charged": n,
            "phi_min": res["phi_min"],
            "ratio_to_vev": res["ratio"],
            "m_phi_squared": res["m_phi_squared"],
            "deviation_percent": res["coincidence_with_alpha_matching"],
        })

    return {
        "phi_vev": phi_vev,
        "results": results,
    }


# ===================================================================
# 6. alpha_running
# ===================================================================
def alpha_running(phi_vev, a_0=1.0, mu_low=None, mu_high=None, n_steps=50):
    """One-loop running of the gauge coupling with KK tower thresholds.

    At each scale mu, all KK modes with m_l < mu are active.
    beta coefficients:
      - charged scalar:  b_l = +1/3  per real dof  (2 per complex)
      - gauge vector:    b_l = -11/3 per adjoint vector

    For 1 complex scalar + 1 vector at each KK level l (degeneracy 2l+1):
      b_l_total = (2l+1) * [2*(1/3) + (-11/3)]  = (2l+1) * (-3)

    Running:  1/alpha(mu) = 1/alpha(M_KK) - b_eff/(2 pi) ln(mu / M_KK)
    where b_eff accumulates as modes are crossed.
    """
    R_phys = a_0 * math.exp(phi_vev)

    if mu_low is None:
        mu_low = 1.0e-3   # eV scale (lab)
    if mu_high is None:
        mu_high = 1.0e28   # Planck scale (eV)

    # build KK mass table up to mu_high
    kk_masses = []
    l = 1
    while True:
        m_l = math.sqrt(l * (l + 1)) / R_phys
        if m_l > mu_high * 10.0:
            break
        kk_masses.append((l, m_l))
        l += 1
        if l > 10000:
            break

    # logarithmic mu grid
    log_mu_lo = math.log10(mu_low)
    log_mu_hi = math.log10(mu_high)
    step = (log_mu_hi - log_mu_lo) / (n_steps - 1)

    alpha_compactification = ALPHA_EM  # matched at the compactification scale

    # the reference scale is the lightest KK mass
    if not kk_masses:
        return {
            "running": [],
            "alpha_at_lab": ALPHA_EM,
            "alpha_at_compactification": ALPHA_EM,
            "description": "No KK modes below mu_high; no running.",
        }

    m_1 = kk_masses[0][1]  # lightest KK mass

    running = []
    for i in range(n_steps):
        log_mu = log_mu_lo + i * step
        mu = 10.0 ** log_mu

        # count active modes and accumulate beta
        n_active = 0
        b_eff = 0.0
        for (l_val, m_l) in kk_masses:
            if m_l < mu:
                deg = 2 * l_val + 1
                # scalar: 2/3 per complex scalar, gauge: -11/3 per vector
                b_eff += deg * (2.0 / 3.0 - 11.0 / 3.0)  # = deg * (-3)
                n_active += 1

        # running from m_1
        if mu > m_1 and m_1 > 0:
            inv_alpha = 1.0 / alpha_compactification - b_eff / (2.0 * math.pi) * math.log(mu / m_1)
        else:
            inv_alpha = 1.0 / alpha_compactification

        alpha_mu = 1.0 / inv_alpha if inv_alpha > 0 else float("inf")

        running.append({
            "mu": mu,
            "alpha_at_mu": alpha_mu,
            "n_active_modes": n_active,
        })

    alpha_lab = running[0]["alpha_at_mu"]

    return {
        "running": running,
        "alpha_at_lab": alpha_lab,
        "alpha_at_compactification": alpha_compactification,
        "description": (
            "One-loop running with KK thresholds.  "
            "alpha at lab ({:.1e} eV) = {:.6e}.  "
            "alpha at compactification = {:.6e}.  "
            "{} KK modes below Planck scale.".format(
                mu_low, alpha_lab, alpha_compactification, len(kk_masses)
            )
        ),
    }


# ===================================================================
# 7. consistency_with_G
# ===================================================================
def consistency_with_G(phi_vev):
    """Check whether the dilaton vev relates to the golden ratio or
    other Alpha Ladder quantities.

    The dilaton vev phi_dilaton = -0.597 is unrelated a priori to the
    golden ratio phi_golden = 1.618.  We check systematically for
    any clean relationship.
    """
    phi_g = PHI_GOLDEN
    e_phi = math.exp(phi_vev)
    e_2phi = math.exp(2.0 * phi_vev)

    checks = []

    # check exp(k * phi_vev) vs phi_golden for integer and half-integer k
    for k_num, k_den in [
        (1, 1), (-1, 1), (2, 1), (-2, 1), (3, 1), (-3, 1), (4, 1), (-4, 1),
        (1, 2), (-1, 2), (3, 2), (-3, 2),
    ]:
        k = k_num / k_den
        val = math.exp(k * phi_vev)
        ratio = val / phi_g
        is_close = abs(ratio - round(ratio)) < 0.01
        checks.append({
            "test": "exp({:.1f} * phi_vev) vs phi_golden".format(k),
            "value": val,
            "ratio_to_phi_golden": ratio,
            "close_to_integer_ratio": is_close,
        })

    # check e^{2 phi_vev} vs phi_golden^2 / 2 (the bridge coefficient)
    bridge = phi_g ** 2 / 2.0
    ratio_bridge = e_2phi / bridge
    checks.append({
        "test": "e^{2 phi_vev} / (phi_golden^2/2)",
        "value": ratio_bridge,
        "ratio_to_phi_golden": ratio_bridge,
        "close_to_integer_ratio": abs(ratio_bridge - round(ratio_bridge)) < 0.01,
    })

    # check e^{phi_vev} vs alpha_EM
    ratio_alpha = e_phi / ALPHA_EM
    checks.append({
        "test": "e^{phi_vev} / alpha_EM",
        "value": ratio_alpha,
        "ratio_to_phi_golden": None,
        "close_to_integer_ratio": abs(ratio_alpha - round(ratio_alpha)) < 0.05,
    })

    # check if ln(phi_golden) / phi_vev is a simple fraction
    ln_phi_g = math.log(phi_g)
    ratio_ln = ln_phi_g / phi_vev
    checks.append({
        "test": "ln(phi_golden) / phi_vev",
        "value": ratio_ln,
        "ratio_to_phi_golden": None,
        "close_to_integer_ratio": abs(ratio_ln - round(ratio_ln)) < 0.05,
    })

    any_match = any(c["close_to_integer_ratio"] for c in checks
                    if c["close_to_integer_ratio"] is not None)

    return {
        "phi_vev": phi_vev,
        "phi_golden": phi_g,
        "relationships_checked": checks,
        "any_match": any_match,
        "honest_assessment": (
            "The dilaton vev phi_dilaton = {:.6f} and the golden ratio "
            "phi_golden = {:.6f} are numerically unrelated.  The former "
            "is fixed by alpha_EM via the gauge matching condition; the "
            "latter appears in the bridge coefficient through an independent "
            "algebraic structure (vacuum polynomial).  No clean relationship "
            "exp(k * phi_dilaton) = phi_golden^n was found for simple k, n.  "
            "This is expected: the gauge matching and the bridge formula "
            "address different aspects of the framework.".format(
                phi_vev, phi_g
            )
        ),
    }


# ===================================================================
# 8. summarize_gauge_matching
# ===================================================================
def summarize_gauge_matching():
    """Full summary of all gauge-matching and CW results."""
    results = {}

    # 1. tree-level matching
    tree = tree_level_matching()
    results["tree_level"] = tree

    # 2. KK spectrum
    spec = kk_mass_spectrum(tree["phi_vev"], a_0=1.0, l_max=20)
    results["kk_spectrum"] = spec

    # 3. CW potential scan
    cw_scan = coleman_weinberg_scan(a_0=1.0, n_charged=1, l_max=50)
    results["cw_scan"] = {
        "has_minimum": cw_scan["has_minimum"],
        "phi_minimum": cw_scan["phi_minimum"],
        "V_at_minimum": cw_scan["V_at_minimum"],
    }

    # 4. CW minimum (refined)
    cw_min = find_cw_minimum(a_0=1.0, n_charged=1, l_max=50)
    results["cw_minimum"] = cw_min

    # 5. n_charged scan
    n_scan = scan_n_charged()
    results["n_charged_scan"] = n_scan

    # 6. alpha running
    ar = alpha_running(tree["phi_vev"])
    results["alpha_running"] = {
        "alpha_at_lab": ar["alpha_at_lab"],
        "alpha_at_compactification": ar["alpha_at_compactification"],
        "description": ar["description"],
    }

    # 7. consistency with G
    g_check = consistency_with_G(tree["phi_vev"])
    results["G_consistency"] = {
        "any_match": g_check["any_match"],
        "honest_assessment": g_check["honest_assessment"],
    }

    # 8. overall summary
    results["summary"] = {
        "key_result_1": (
            "Tree-level matching works exactly: setting alpha_KK = alpha_EM "
            "determines phi_vev = {:.6f}.".format(tree["phi_vev"])
        ),
        "key_result_2": (
            "CW minimum: {} (phi_min = {}).".format(
                "found" if cw_min["phi_min"] is not None else "not found",
                "{:.6f}".format(cw_min["phi_min"]) if cw_min["phi_min"] is not None else "N/A",
            )
        ),
        "key_result_3": (
            "Coincidence with alpha matching: {}%.".format(
                "{:.2f}".format(cw_min["coincidence_with_alpha_matching"])
                if cw_min["coincidence_with_alpha_matching"] is not None else "N/A"
            )
        ),
        "gap_1_status": (
            "Gap 1 (vacuum polynomial / radius determination) remains open "
            "regardless of CW results."
        ),
        "gap_2_status": (
            "Gap 2 (correction coefficients c_2, c_3) remains open -- "
            "the one-loop KK contribution is too small by ~40 orders."
        ),
        "honest_caveat": (
            "Adding charged matter is an extension beyond the minimal "
            "pure-gravity framework.  The CW potential is one-loop, "
            "scheme-dependent, and sensitive to the renormalisation scale.  "
            "A coincidence between phi_CW_min and phi_alpha_match, while "
            "interesting, would not constitute a derivation of alpha_EM "
            "from geometry without a UV-complete treatment."
        ),
    }

    return results


# ===================================================================
# __main__ -- print full report
# ===================================================================
if __name__ == "__main__":

    print("=" * 72)
    print("KK Gauge Matching & Coleman-Weinberg Potential")
    print("Alpha Ladder S^2 Framework")
    print("=" * 72)

    # --- 1. Tree-level matching ---
    print("\n--- 1. Tree-level matching ---")
    tree = tree_level_matching()
    print("  phi_vev           = {:.8f}".format(tree["phi_vev"]))
    print("  g_KK              = {:.8e}".format(tree["g_KK"]))
    print("  alpha_KK          = {:.12e}".format(tree["alpha_KK"]))
    print("  alpha_EM          = {:.12e}".format(tree["alpha_EM"]))
    print("  e^(phi_vev)       = {:.8f}".format(tree["e_to_phi"]))
    print("  e^(2 phi_vev)     = {:.8f}".format(tree["e_to_2phi"]))
    print("  e^(4 phi_vev)     = {:.8f}".format(tree["e_to_4phi"]))
    print("  Verification      = {}".format(tree["verification"]))

    # --- 2. KK mass spectrum ---
    print("\n--- 2. KK mass spectrum (a_0 = 1, Planck units) ---")
    spec = kk_mass_spectrum(tree["phi_vev"], a_0=1.0, l_max=10)
    print("  R_phys = {:.6f} (Planck lengths)".format(spec["R_phys"]))
    print("  {:>4s}  {:>14s}  {:>10s}  {:>5s}".format(
        "l", "m_l", "m_l/m_1", "deg"))
    for m in spec["modes"]:
        print("  {:4d}  {:14.6e}  {:10.6f}  {:5d}".format(
            m["l"], m["m_l"], m["m_l_over_m1"], m["degeneracy"]))

    # --- 3. Coleman-Weinberg scan ---
    print("\n--- 3. Coleman-Weinberg potential scan ---")
    cw = coleman_weinberg_scan(a_0=1.0, n_charged=1, l_max=50)
    print("  Minimum found: {}".format(cw["has_minimum"]))
    if cw["has_minimum"]:
        print("  phi_minimum   = {:.6f}".format(cw["phi_minimum"]))
        print("  V_at_minimum  = {:.6e}".format(cw["V_at_minimum"]))

    # --- 4. Refined CW minimum ---
    print("\n--- 4. CW minimum (bisection-refined) ---")
    cw_min = find_cw_minimum(a_0=1.0, n_charged=1, l_max=50)
    print("  " + cw_min["description"])
    print("  " + cw_min["honest_assessment"])

    # --- 5. n_charged scan ---
    print("\n--- 5. Scan over n_charged ---")
    n_scan = scan_n_charged()
    print("  phi_vev (tree) = {:.6f}".format(n_scan["phi_vev"]))
    print("  {:>5s}  {:>12s}  {:>10s}  {:>14s}  {:>10s}".format(
        "n", "phi_min", "ratio", "m_phi^2", "dev %"))
    for r in n_scan["results"]:
        phi_str = "{:.6f}".format(r["phi_min"]) if r["phi_min"] is not None else "N/A"
        ratio_str = "{:.4f}".format(r["ratio_to_vev"]) if r["ratio_to_vev"] is not None else "N/A"
        m_str = "{:.4e}".format(r["m_phi_squared"]) if r["m_phi_squared"] is not None else "N/A"
        d_str = "{:.2f}".format(r["deviation_percent"]) if r["deviation_percent"] is not None else "N/A"
        print("  {:5d}  {:>12s}  {:>10s}  {:>14s}  {:>10s}".format(
            r["n_charged"], phi_str, ratio_str, m_str, d_str))

    # --- 6. Alpha running ---
    print("\n--- 6. Alpha running ---")
    ar = alpha_running(tree["phi_vev"])
    print("  " + ar["description"])
    print("  Sample points:")
    for pt in ar["running"][::max(1, len(ar["running"]) // 8)]:
        print("    mu = {:.2e}  alpha = {:.6e}  active modes = {}".format(
            pt["mu"], pt["alpha_at_mu"], pt["n_active_modes"]))

    # --- 7. Consistency with G ---
    print("\n--- 7. Consistency with G / golden ratio ---")
    g_check = consistency_with_G(tree["phi_vev"])
    print("  " + g_check["honest_assessment"])
    for c in g_check["relationships_checked"]:
        flag = " <--" if c.get("close_to_integer_ratio") else ""
        print("    {}: value = {:.6f}{}".format(
            c["test"], c["value"], flag))

    # --- 8. Summary ---
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    summary = summarize_gauge_matching()
    for key, val in summary["summary"].items():
        print("\n  [{}]".format(key))
        print("  " + str(val))

    print("\n" + "=" * 72)
    print("END OF REPORT")
    print("=" * 72)
