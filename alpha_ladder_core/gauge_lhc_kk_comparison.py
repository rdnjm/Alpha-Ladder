"""
LHC KK Graviton Cross-Section Ratio: S^2 vs ADD
=================================================

Computes the ratio of KK graviton production cross-sections between the
Alpha Ladder S^2 compactification and the standard ADD torus (T^2)
compactification used in LHC exclusion analyses.

The LHC excludes M_6 > 5 TeV for ADD with n=2 extra dimensions via
mono-jet + missing energy searches (gg -> G_KK + jet).  The signal
strength is proportional to the effective number of accessible KK modes
N_eff(sqrt(s)), which differs between T^2 and S^2 due to:

  1. Different spectra: T^2 has m_{n1,n2} = sqrt(n1^2+n2^2)/R,
     S^2 has m_l = sqrt(l(l+1))/R_phys with degeneracy 2l+1.
  2. Different lightest modes: T^2 lightest is 1/R, S^2 is sqrt(2)/R_phys.
  3. Nearly identical radii at the same M_6 (within ~2%).

If the S^2 spectrum produces fewer accessible modes (ratio < 1), the
LHC exclusion is loosened and the M_6 bound is weaker.

The mode counts are enormous (~ 10^30 at M_6 = 5 TeV, E_cm = 14 TeV),
so we use:
  - Analytic formulas for N_eff and phase-space weighted sums
  - Exact lattice enumeration only for the first ~20 mass levels
    in the spectrum comparison table

Pure Python -- only ``import math``, no numpy/scipy.

Functions
---------
1. add_kk_spectrum           -- T^2 first N_levels KK mass levels (exact)
2. s2_kk_spectrum            -- S^2 KK modes below m_max (exact, fast)
3. effective_modes           -- N_eff for both models at given M_6 (analytic)
4. cross_section_ratio       -- scan M_6 for the ratio N_eff(S2)/N_eff(ADD)
5. rescaled_exclusion        -- rescale the ADD M_6 > 5 TeV bound to S^2
6. spectrum_comparison_table -- side-by-side first 20 KK levels
7. summarize_lhc_comparison  -- full report
"""

import math

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
ALPHA_EM = 1.0 / 137.035999084            # CODATA 2018
M_PL_EV = 1.22089e28                       # Planck mass (eV)
HBAR_C_EVM = 1.9733e-7                     # hbar*c in eV m
EV_PER_TEV = 1.0e12                        # eV per TeV

# Derived
PHI_VEV = 0.25 * math.log(4.0 * math.pi * ALPHA_EM)   # dilaton vev ~ -0.597
E_PHI_VEV = math.exp(PHI_VEV)                          # e^{phi_vev} ~ 0.5503


# ===================================================================
# Helper: radii from M_6
# ===================================================================
def _R_add(M_6_eV):
    """ADD torus radius in eV^{-1}: R = M_Pl / (2*pi*M_6^2)."""
    return M_PL_EV / (2.0 * math.pi * M_6_eV ** 2)


def _R_s2(M_6_eV):
    """S^2 physical radius in eV^{-1}: R_phys = a_0*e^{phi_vev},
    where a_0 = M_Pl / (sqrt(4*pi)*M_6^2)."""
    a_0 = M_PL_EV / (math.sqrt(4.0 * math.pi) * M_6_eV ** 2)
    return a_0 * E_PHI_VEV


def _a0_from_M6(M_6_eV):
    """Bare radius a_0 in eV^{-1} from M_6 in eV."""
    return M_PL_EV / (math.sqrt(4.0 * math.pi) * M_6_eV ** 2)


# ===================================================================
# Analytic mode counts
# ===================================================================
def _N_eff_add_analytic(R, E):
    """Analytic N_eff for ADD T^2: pi*(E*R)^2 (area of disk, minus origin).

    For n=2 flat extra dimensions, the number of lattice points (n1,n2)
    with n1^2+n2^2 <= (E*R)^2 is pi*(E*R)^2 to leading order (Gauss
    circle problem).  We subtract 1 for the (0,0) mode.
    """
    x = E * R
    return math.pi * x * x - 1.0


def _N_eff_s2_analytic(R_phys, E):
    """Analytic N_eff for S^2: sum_{l=1}^{l_max} (2l+1) where l_max
    is the largest l with sqrt(l(l+1))/R_phys < E.

    l(l+1) < (E*R_phys)^2
    l_max ~ E*R_phys - 1/2  (for large E*R)

    Exact: l_max = floor((-1 + sqrt(1 + 4*(E*R_phys)^2)) / 2)
    Sum_{l=1}^{l_max} (2l+1) = (l_max+1)^2 - 1
    """
    x_sq = (E * R_phys) ** 2
    discriminant = 1.0 + 4.0 * x_sq
    l_max_float = (-1.0 + math.sqrt(discriminant)) / 2.0
    l_max = int(l_max_float)
    if l_max < 1:
        return 0.0
    return float((l_max + 1) ** 2 - 1)


def _sigma_weighted_add_analytic(R, E):
    """Phase-space weighted sum for ADD T^2 (analytic).

    sigma ~ sum_{modes with m<E} (1 - m^2/E^2)

    For a continuous approximation with n=2:
    integral_0^E rho(m) * (1 - m^2/E^2) dm
    where rho(m) = 2*pi*m*R^2

    = 2*pi*R^2 * integral_0^E m*(1 - m^2/E^2) dm
    = 2*pi*R^2 * [E^2/2 - E^2/4]
    = 2*pi*R^2 * E^2/4
    = pi*R^2*E^2 / 2
    = N_eff / 2  (since N_eff = pi*R^2*E^2)
    """
    x = E * R
    return math.pi * x * x / 2.0


def _sigma_weighted_s2_analytic(R_phys, E):
    """Phase-space weighted sum for S^2 (analytic).

    sigma ~ sum_{l=1}^{l_max} (2l+1) * (1 - l(l+1)/(E*R)^2)

    = sum (2l+1) - (1/(E*R)^2) * sum (2l+1)*l(l+1)

    sum_{l=1}^L (2l+1) = (L+1)^2 - 1

    sum_{l=1}^L (2l+1)*l(l+1) = sum (2l^3 + 3l^2 + l)
      = 2 * [L(L+1)/2]^2 + 3*L(L+1)(2L+1)/6 + L(L+1)/2
      = L^2(L+1)^2/2 + L(L+1)(2L+1)/2 + L(L+1)/2
      = L(L+1)/2 * [L(L+1) + (2L+1) + 1]
      = L(L+1)/2 * [(L+1)^2 + L]  ... let me just compute directly.
    """
    x_sq = (E * R_phys) ** 2
    discriminant = 1.0 + 4.0 * x_sq
    l_max_float = (-1.0 + math.sqrt(discriminant)) / 2.0
    L = int(l_max_float)
    if L < 1:
        return 0.0

    # sum_{l=1}^L (2l+1) = (L+1)^2 - 1
    S1 = float((L + 1) ** 2 - 1)

    # sum_{l=1}^L (2l+1)*l*(l+1)
    # = sum 2l^3 + 3l^2 + l
    # sum l = L(L+1)/2
    # sum l^2 = L(L+1)(2L+1)/6
    # sum l^3 = [L(L+1)/2]^2
    sum_l = L * (L + 1) / 2.0
    sum_l2 = L * (L + 1) * (2 * L + 1) / 6.0
    sum_l3 = (L * (L + 1) / 2.0) ** 2
    S2 = 2.0 * sum_l3 + 3.0 * sum_l2 + sum_l

    sigma = S1 - S2 / x_sq
    return sigma


# ===================================================================
# 1. ADD (torus) KK spectrum -- first N levels (exact enumeration)
# ===================================================================
def add_kk_spectrum(R, m_max, n_extra=2):
    """Compute the ADD T^2 KK spectrum: first mass levels by exact enumeration.

    For the spectrum table we only need the first ~20 distinct mass levels.
    The masses are sqrt(k)/R where k = n1^2 + n2^2 for integers n1, n2.
    We enumerate k values (sums of two squares) in order.

    Parameters
    ----------
    R : float
        Compactification radius in eV^{-1} (natural units).
    m_max : float
        Maximum mass in eV.  Only modes with m < m_max are returned.
    n_extra : int
        Number of extra dimensions (must be 2).

    Returns
    -------
    dict with keys:
        N_modes : int
            Total number of KK modes (analytic) with m < m_max.
        masses : list of float
            Sorted list of distinct mass values in eV (first levels only).
        mass_degeneracies : list of (float, int)
            (mass, degeneracy) pairs sorted by mass (first levels only).
        rho_analytic : float
            Analytic density of states dN/dm at m_max: 2*pi*m_max*R^2.
        N_analytic : float
            Analytic count pi*(m_max*R)^2.
    """
    if n_extra != 2:
        raise ValueError("Only n_extra=2 (T^2) is implemented")

    # We enumerate sums of two squares k = n1^2 + n2^2 up to k_max.
    # For the first ~30 distinct k values, we only need to scan a small grid.
    # k_max for 30 levels is at most ~50 (plenty of margin).
    k_max_limit = int((m_max * R) ** 2)
    # Cap the enumeration to keep it fast -- only need first levels
    k_max_enum = min(k_max_limit, 200)

    l_scan = int(math.sqrt(k_max_enum)) + 1

    mass_sq_counts = {}  # k -> count of (n1,n2) pairs with n1^2+n2^2 = k

    for n1 in range(-l_scan, l_scan + 1):
        for n2 in range(-l_scan, l_scan + 1):
            k = n1 * n1 + n2 * n2
            if k == 0:
                continue
            if k <= k_max_enum:
                mass_sq_counts[k] = mass_sq_counts.get(k, 0) + 1

    # Build sorted list of (mass, degeneracy)
    mass_degens = []
    for k, count in sorted(mass_sq_counts.items()):
        m = math.sqrt(k) / R
        if m < m_max:
            mass_degens.append((m, count))

    masses = [md[0] for md in mass_degens]

    rho_analytic = 2.0 * math.pi * m_max * R * R
    N_analytic = math.pi * (m_max * R) ** 2

    # Use analytic total for N_modes (exact enumeration only for first levels)
    N_modes = int(round(N_analytic))

    return {
        "N_modes": N_modes,
        "masses": masses,
        "mass_degeneracies": mass_degens,
        "rho_analytic": rho_analytic,
        "N_analytic": N_analytic,
    }


# ===================================================================
# 2. S^2 KK spectrum
# ===================================================================
def s2_kk_spectrum(R_phys, m_max):
    """Compute the S^2 KK spectrum.

    Parameters
    ----------
    R_phys : float
        Physical radius of S^2 in eV^{-1} (natural units).
        R_phys = a_0 * e^{phi_vev} where a_0 is the bare radius.
    m_max : float
        Maximum mass in eV.

    Returns
    -------
    dict with keys:
        N_modes : int
            Total number of modes (counting degeneracy 2l+1) with m_l < m_max.
        mode_list : list of (int, float, int)
            (l, m_l, degeneracy) for each angular momentum level.
        N_analytic : float
            Asymptotic count (m_max*R_phys)^2.
    """
    mode_list = []
    N_modes = 0

    l = 1
    while True:
        m_l = math.sqrt(l * (l + 1)) / R_phys
        if m_l >= m_max:
            break
        deg = 2 * l + 1
        mode_list.append((l, m_l, deg))
        N_modes += deg
        l += 1

    N_analytic = (m_max * R_phys) ** 2

    return {
        "N_modes": N_modes,
        "mode_list": mode_list,
        "N_analytic": N_analytic,
    }


# ===================================================================
# 3. Effective modes at given M_6 (analytic)
# ===================================================================
def effective_modes(model, M_6_TeV, E_cm_TeV=14.0):
    """Compute N_eff = accessible KK modes with mass < E_cm for a given M_6.

    Uses analytic formulas since mode counts are ~ 10^30 at LHC energies.

    Parameters
    ----------
    model : str
        "ADD" for T^2 torus, "S2" for S^2 sphere.
    M_6_TeV : float
        6D Planck mass in TeV.
    E_cm_TeV : float
        Centre-of-mass energy in TeV (default: 14 TeV LHC).

    Returns
    -------
    dict with keys:
        N_eff : float
            Total accessible KK modes (analytic).
        R : float
            Compactification radius in eV^{-1}.
        R_meters : float
            Compactification radius in meters.
        m_lightest_eV : float
            Lightest KK mode mass in eV.
        sigma_weighted : float
            Phase-space weighted effective cross-section sum.
    """
    M_6_eV = M_6_TeV * EV_PER_TEV
    E_cm_eV = E_cm_TeV * EV_PER_TEV

    if model.upper() == "ADD":
        R = _R_add(M_6_eV)
        N_eff = _N_eff_add_analytic(R, E_cm_eV)
        m_lightest = 1.0 / R   # lightest mode: (1,0) or (0,1)
        sigma_w = _sigma_weighted_add_analytic(R, E_cm_eV)
    elif model.upper() == "S2":
        R = _R_s2(M_6_eV)
        N_eff = _N_eff_s2_analytic(R, E_cm_eV)
        m_lightest = math.sqrt(2.0) / R   # l=1 mode
        sigma_w = _sigma_weighted_s2_analytic(R, E_cm_eV)
    else:
        raise ValueError("model must be 'ADD' or 'S2'")

    R_meters = R * HBAR_C_EVM

    return {
        "N_eff": N_eff,
        "R": R,
        "R_meters": R_meters,
        "m_lightest_eV": m_lightest,
        "sigma_weighted": sigma_w,
    }


def _compute_both_models(M_6_TeV, E_cm_TeV=14.0):
    """Compute N_eff for both models and their ratio.

    Returns
    -------
    dict with keys:
        N_eff_add, N_eff_s2, ratio, R_add, R_s2,
        R_add_m, R_s2_m, m_lightest_add, m_lightest_s2,
        sigma_ratio.
    """
    add = effective_modes("ADD", M_6_TeV, E_cm_TeV)
    s2 = effective_modes("S2", M_6_TeV, E_cm_TeV)

    ratio = s2["N_eff"] / add["N_eff"] if add["N_eff"] > 0 else float("inf")
    sig_ratio = (s2["sigma_weighted"] / add["sigma_weighted"]
                 if add["sigma_weighted"] > 0 else float("inf"))

    return {
        "N_eff_add": add["N_eff"],
        "N_eff_s2": s2["N_eff"],
        "ratio": ratio,
        "sigma_ratio": sig_ratio,
        "R_add": add["R"],
        "R_s2": s2["R"],
        "R_add_m": add["R_meters"],
        "R_s2_m": s2["R_meters"],
        "m_lightest_add": add["m_lightest_eV"],
        "m_lightest_s2": s2["m_lightest_eV"],
    }


# ===================================================================
# 4. Cross-section ratio scan
# ===================================================================
def cross_section_ratio(M_6_TeV_values=None, E_cm_TeV=14.0):
    """Scan M_6 and compute the ratio N_eff(S^2)/N_eff(ADD) at each.

    Also computes a phase-space-weighted ratio:
        sigma_ratio = sigma_weighted(S^2) / sigma_weighted(ADD)

    Parameters
    ----------
    M_6_TeV_values : list of float or None
        M_6 values to scan in TeV.  Default: 1 to 20 in steps of 1.
    E_cm_TeV : float
        Centre-of-mass energy in TeV.

    Returns
    -------
    list of dict, each with keys:
        M_6, N_ratio, sigma_ratio, N_add, N_s2, sigma_add, sigma_s2.
    """
    if M_6_TeV_values is None:
        M_6_TeV_values = [float(i) for i in range(1, 21)]

    results = []

    for M_6 in M_6_TeV_values:
        both = _compute_both_models(M_6, E_cm_TeV)

        add_data = effective_modes("ADD", M_6, E_cm_TeV)
        s2_data = effective_modes("S2", M_6, E_cm_TeV)

        results.append({
            "M_6": M_6,
            "N_ratio": both["ratio"],
            "sigma_ratio": both["sigma_ratio"],
            "N_add": both["N_eff_add"],
            "N_s2": both["N_eff_s2"],
            "sigma_add": add_data["sigma_weighted"],
            "sigma_s2": s2_data["sigma_weighted"],
        })

    return results


# ===================================================================
# 5. Rescaled exclusion
# ===================================================================
def rescaled_exclusion(M_6_excluded_ADD=5.0, E_cm_TeV=14.0):
    """Rescale the ADD LHC exclusion M_6 > 5 TeV to S^2 geometry.

    The LHC excludes ADD models where the signal exceeds sigma_limit.
    At M_6 = 5 TeV, sigma(ADD) = sigma_limit.  For S^2, we find the
    M_6 where sigma(S^2, M_6) = sigma(ADD, 5 TeV), i.e. where the
    S^2 model produces the same signal strength.

    Uses analytic formulas.  The cross-section is proportional to
    N_eff / M_Pl^4, and since M_Pl is the 4D Planck mass (fixed), the
    signal scales directly with N_eff.  We solve:
        N_eff(M_6_excl, S2) = N_eff(5 TeV, ADD)

    Also solves using the phase-space weighted sigma.

    Returns
    -------
    dict with keys:
        M_6_excl_add : float
            ADD exclusion (input), TeV.
        M_6_excl_s2 : float
            Rescaled S^2 exclusion (from N_eff matching), TeV.
        M_6_excl_s2_sigma : float
            Rescaled S^2 exclusion (from sigma matching), TeV.
        ratio_at_exclusion : float
            N_eff(S2)/N_eff(ADD) at the ADD exclusion point.
        sigma_ratio_at_exclusion : float
            sigma(S2)/sigma(ADD) at the ADD exclusion point.
        N_eff_target : float
            N_eff(ADD) at the exclusion point.
        sigma_target : float
            Phase-space weighted sigma(ADD) at the exclusion point.
        a_0_threshold_s2_m : float
            Bare radius a_0 in meters at the S^2 exclusion.
        R_phys_threshold_s2_m : float
            Physical radius in meters at the S^2 exclusion.
    """
    E_cm_eV = E_cm_TeV * EV_PER_TEV

    # Target: values at ADD exclusion point
    add_ref = effective_modes("ADD", M_6_excluded_ADD, E_cm_TeV)
    s2_ref = effective_modes("S2", M_6_excluded_ADD, E_cm_TeV)
    N_target = add_ref["N_eff"]
    sigma_target = add_ref["sigma_weighted"]

    ratio_at_excl = s2_ref["N_eff"] / N_target if N_target > 0 else float("inf")
    sigma_ratio_at_excl = (s2_ref["sigma_weighted"] / sigma_target
                           if sigma_target > 0 else float("inf"))

    # Binary search for M_6 where N_eff(S2, M_6) = N_target
    # Lower M_6 => larger R => more modes.  N_eff is monotonically
    # decreasing in M_6.
    M_lo, M_hi = 0.5, 50.0

    for _ in range(200):
        M_mid = (M_lo + M_hi) / 2.0
        s2_mid = effective_modes("S2", M_mid, E_cm_TeV)
        if s2_mid["N_eff"] > N_target:
            M_lo = M_mid
        else:
            M_hi = M_mid
        if M_hi - M_lo < 1e-6:
            break

    M_6_excl_s2 = (M_lo + M_hi) / 2.0

    # Phase-space weighted binary search
    M_lo_s, M_hi_s = 0.5, 50.0
    for _ in range(200):
        M_mid = (M_lo_s + M_hi_s) / 2.0
        s2_mid = effective_modes("S2", M_mid, E_cm_TeV)
        if s2_mid["sigma_weighted"] > sigma_target:
            M_lo_s = M_mid
        else:
            M_hi_s = M_mid
        if M_hi_s - M_lo_s < 1e-6:
            break

    M_6_excl_s2_sigma = (M_lo_s + M_hi_s) / 2.0

    # Compute a_0 and R_phys at the S^2 N_eff exclusion
    M_excl_eV = M_6_excl_s2 * EV_PER_TEV
    a_0_nat = _a0_from_M6(M_excl_eV)
    a_0_m = a_0_nat * HBAR_C_EVM
    R_phys_m = a_0_m * E_PHI_VEV

    return {
        "M_6_excl_add": M_6_excluded_ADD,
        "M_6_excl_s2": M_6_excl_s2,
        "M_6_excl_s2_sigma": M_6_excl_s2_sigma,
        "ratio_at_exclusion": ratio_at_excl,
        "sigma_ratio_at_exclusion": sigma_ratio_at_excl,
        "N_eff_target": N_target,
        "sigma_target": sigma_target,
        "a_0_threshold_s2_m": a_0_m,
        "R_phys_threshold_s2_m": R_phys_m,
    }


# ===================================================================
# 6. Spectrum comparison table
# ===================================================================
def spectrum_comparison_table(M_6_TeV=5.0):
    """Side-by-side comparison of the first 20 KK mass levels for ADD and S^2.

    Uses exact lattice enumeration for ADD (first ~20 distinct k values)
    and exact angular momentum enumeration for S^2.

    Parameters
    ----------
    M_6_TeV : float
        6D Planck mass in TeV.

    Returns
    -------
    dict with keys:
        M_6_TeV : float
        R_add_m, R_s2_m : float
            Radii in meters.
        add_levels : list of (float, int)
            First 20 (mass_eV, degeneracy) for ADD.
        s2_levels : list of (int, float, int)
            First 20 (l, mass_eV, degeneracy) for S^2.
        add_cumulative : list of int
            Cumulative mode count after each ADD level.
        s2_cumulative : list of int
            Cumulative mode count after each S^2 level.
    """
    M_6_eV = M_6_TeV * EV_PER_TEV

    # ADD radius and spectrum (first levels via capped enumeration)
    R_add = _R_add(M_6_eV)
    # Use m_max = 30/R to get first ~20 distinct mass levels
    m_max_add = 30.0 / R_add
    add_spec = add_kk_spectrum(R_add, m_max_add)
    add_levels = add_spec["mass_degeneracies"][:20]
    add_cum = []
    running = 0
    for _, deg in add_levels:
        running += deg
        add_cum.append(running)

    # S2 radius and spectrum
    R_s2 = _R_s2(M_6_eV)
    m_max_s2 = 30.0 / R_s2
    s2_spec = s2_kk_spectrum(R_s2, m_max_s2)
    s2_levels = s2_spec["mode_list"][:20]
    s2_cum = []
    running = 0
    for _, _, deg in s2_levels:
        running += deg
        s2_cum.append(running)

    return {
        "M_6_TeV": M_6_TeV,
        "R_add_m": R_add * HBAR_C_EVM,
        "R_s2_m": R_s2 * HBAR_C_EVM,
        "add_levels": add_levels,
        "s2_levels": s2_levels,
        "add_cumulative": add_cum,
        "s2_cumulative": s2_cum,
    }


# ===================================================================
# 7. Summary report
# ===================================================================
def summarize_lhc_comparison():
    """Print a comprehensive report comparing LHC KK signals for S^2 vs ADD.

    Returns
    -------
    dict with all computed results.
    """
    lines = []
    lines.append("=" * 72)
    lines.append("LHC KK Graviton Cross-Section: S^2 (Alpha Ladder) vs ADD (T^2)")
    lines.append("=" * 72)

    # --- Radius comparison ---
    lines.append("")
    lines.append("1. RADIUS COMPARISON AT SAME M_6")
    lines.append("-" * 40)
    for M_6 in [3.0, 5.0, 10.0, 15.0]:
        both = _compute_both_models(M_6)
        lines.append(
            "  M_6 = {:5.1f} TeV:  R_ADD = {:.3e} m,  R_S2 = {:.3e} m,  "
            "R_S2/R_ADD = {:.4f}".format(
                M_6, both["R_add_m"], both["R_s2_m"],
                both["R_s2"] / both["R_add"] if both["R_add"] > 0 else 0
            )
        )

    # Analytic ratio
    R_ratio_analytic = E_PHI_VEV * 2.0 * math.pi / math.sqrt(4.0 * math.pi)
    lines.append("")
    lines.append("  Analytic R_S2/R_ADD = e^{{phi_vev}} * 2*pi / sqrt(4*pi)")
    lines.append("                      = {:.6f} * {:.6f} / {:.6f}".format(
        E_PHI_VEV, 2.0 * math.pi, math.sqrt(4.0 * math.pi)))
    lines.append("                      = {:.6f}".format(R_ratio_analytic))

    # --- Spectrum table ---
    lines.append("")
    lines.append("2. SPECTRUM COMPARISON (M_6 = 5 TeV)")
    lines.append("-" * 40)
    table = spectrum_comparison_table(5.0)
    lines.append("  R_ADD = {:.3e} m,  R_S2 = {:.3e} m".format(
        table["R_add_m"], table["R_s2_m"]))
    lines.append("")
    lines.append("  {:<5s} {:>14s} {:>6s} {:>8s}   {:<5s} {:>14s} {:>6s} {:>8s}".format(
        "Level", "m_ADD (eV)", "deg", "cumul",
        "l", "m_S2 (eV)", "deg", "cumul"))
    lines.append("  " + "-" * 68)
    n_show = min(20, len(table["add_levels"]), len(table["s2_levels"]))
    for i in range(n_show):
        m_a, d_a = table["add_levels"][i]
        l_s, m_s, d_s = table["s2_levels"][i]
        lines.append(
            "  {:>5d} {:>14.4e} {:>6d} {:>8d}   {:>5d} {:>14.4e} {:>6d} {:>8d}".format(
                i + 1, m_a, d_a, table["add_cumulative"][i],
                l_s, m_s, d_s, table["s2_cumulative"][i]))

    # --- N_eff scan ---
    lines.append("")
    lines.append("3. CROSS-SECTION RATIO SCAN")
    lines.append("-" * 40)
    lines.append("  E_cm = 14 TeV (LHC design energy)")
    lines.append("  (Mode counts are ~ 10^30; analytic formulas used)")
    lines.append("")
    lines.append("  {:>6s} {:>14s} {:>14s} {:>10s} {:>10s}".format(
        "M_6", "N_ADD", "N_S2", "N_ratio", "sig_ratio"))
    lines.append("  " + "-" * 56)

    scan = cross_section_ratio()
    for row in scan:
        lines.append("  {:>5.1f}T {:>14.6e} {:>14.6e} {:>10.6f} {:>10.6f}".format(
            row["M_6"], row["N_add"], row["N_s2"],
            row["N_ratio"], row["sigma_ratio"]))

    # --- Parton-level scan ---
    lines.append("")
    lines.append("4. PARTON-LEVEL EFFECTIVE ENERGY SCAN (M_6 = 5 TeV)")
    lines.append("-" * 40)
    lines.append("  Typical parton sqrt(s_hat) is 1-5 TeV at the LHC.")
    lines.append("")
    lines.append("  {:>8s} {:>14s} {:>14s} {:>10s} {:>10s}".format(
        "E_part", "N_ADD", "N_S2", "N_ratio", "sig_ratio"))
    lines.append("  " + "-" * 58)

    for E_part in [1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 14.0]:
        both = _compute_both_models(5.0, E_part)
        lines.append("  {:>7.1f}T {:>14.6e} {:>14.6e} {:>10.6f} {:>10.6f}".format(
            E_part, both["N_eff_add"], both["N_eff_s2"],
            both["ratio"], both["sigma_ratio"]))

    # --- Rescaled exclusion ---
    lines.append("")
    lines.append("5. RESCALED LHC EXCLUSION")
    lines.append("-" * 40)
    excl = rescaled_exclusion()
    lines.append("  ADD exclusion:           M_6 > {:.1f} TeV".format(
        excl["M_6_excl_add"]))
    lines.append("  S^2 exclusion (N_eff):   M_6 > {:.3f} TeV".format(
        excl["M_6_excl_s2"]))
    lines.append("  S^2 exclusion (sigma):   M_6 > {:.3f} TeV".format(
        excl["M_6_excl_s2_sigma"]))
    lines.append("  N_ratio at ADD excl:     {:.6f}".format(
        excl["ratio_at_exclusion"]))
    lines.append("  sigma_ratio at ADD excl: {:.6f}".format(
        excl["sigma_ratio_at_exclusion"]))
    lines.append("  N_eff target (ADD):      {:.6e}".format(excl["N_eff_target"]))
    lines.append("  a_0 threshold (S^2):     {:.3e} m".format(
        excl["a_0_threshold_s2_m"]))
    lines.append("  R_phys threshold (S^2):  {:.3e} m".format(
        excl["R_phys_threshold_s2_m"]))

    # --- Eot-Wash window ---
    lines.append("")
    lines.append("6. IMPLICATIONS FOR EOT-WASH WINDOW")
    lines.append("-" * 40)
    # Eot-Wash optimal a_0 ~ 28 um from previous analyses
    a_0_eotwash = 28e-6  # meters
    a_0_nat = a_0_eotwash / HBAR_C_EVM
    M_6_eotwash_eV = (M_PL_EV ** 2 / (4.0 * math.pi * a_0_nat ** 2)) ** 0.25
    M_6_eotwash_TeV = M_6_eotwash_eV / EV_PER_TEV

    lines.append("  Eot-Wash optimal a_0:    {:.0f} um".format(
        a_0_eotwash * 1e6))
    lines.append("  Corresponding M_6:       {:.4f} TeV".format(M_6_eotwash_TeV))
    lines.append("  S^2 LHC exclusion M_6:   {:.3f} TeV".format(
        excl["M_6_excl_s2"]))

    if M_6_eotwash_TeV < excl["M_6_excl_s2"]:
        lines.append("  STATUS: Eot-Wash window is EXCLUDED by LHC "
                      "(M_6 = {:.4f} < {:.3f} TeV)".format(
                          M_6_eotwash_TeV, excl["M_6_excl_s2"]))
    else:
        lines.append("  STATUS: Eot-Wash window is OPEN "
                      "(M_6 = {:.4f} > {:.3f} TeV)".format(
                          M_6_eotwash_TeV, excl["M_6_excl_s2"]))

    # --- Analytic understanding ---
    lines.append("")
    lines.append("7. ANALYTIC UNDERSTANDING OF THE RATIO")
    lines.append("-" * 40)
    lines.append("  For large E*R (continuum limit):")
    lines.append("    N_ADD = pi * (E*R_ADD)^2       [Gauss circle problem]")
    lines.append("    N_S2  = (E*R_S2)^2             [sum of 2l+1 up to l_max ~ E*R]")
    lines.append("")
    lines.append("  Therefore N_S2/N_ADD = (1/pi) * (R_S2/R_ADD)^2")
    lines.append("                       = (1/pi) * ({:.6f})^2".format(R_ratio_analytic))
    lines.append("                       = {:.6f}".format(
        R_ratio_analytic ** 2 / math.pi))
    lines.append("")
    lines.append("  The factor of 1/pi arises because:")
    lines.append("    - T^2 modes fill a 2D disk: N ~ pi*r^2 (area of circle)")
    lines.append("    - S^2 modes fill a 1D tower: N ~ l_max^2 = r^2 (no pi factor)")
    lines.append("  The S^2 spectrum lacks the pi prefactor of the lattice count.")
    lines.append("")
    lines.append("  Phase-space weighted ratio (sigma):")
    lines.append("    sigma_ADD = pi*R^2*E^2 / 2 = N_ADD/2")
    lines.append("    sigma_S2  ~ N_S2/2         [same 1/2 factor in continuum limit]")
    lines.append("    sigma_ratio = N_ratio       [identical in continuum limit]")
    lines.append("  Both spectra give the same phase-space suppression (1/2) for n=2,")
    lines.append("  so the N_eff ratio fully determines the cross-section ratio.")

    # --- Key messages ---
    lines.append("")
    lines.append("8. KEY MESSAGES")
    lines.append("-" * 40)
    lines.append("  a) The S^2 KK spectrum is sparser than ADD for light modes.")
    lines.append("     Lightest S^2 mode: sqrt(2)/R vs ADD: 1/R  (ratio sqrt(2) ~ 1.414).")
    lines.append("  b) Total mode count scales as ~R^2*E^2 for both, but S^2 has")
    lines.append("     a smaller prefactor: N_S2/N_ADD ~ (R_S2/R_ADD)^2/pi ~ {:.4f}.".format(
        R_ratio_analytic ** 2 / math.pi))
    lines.append("  c) At M_6 = 5 TeV, N_eff(S^2)/N_eff(ADD) ~ {:.4f}.".format(
        excl["ratio_at_exclusion"]))
    lines.append("  d) The LHC exclusion for S^2 is WEAKER than ADD.")
    lines.append("     ADD: M_6 > {:.1f} TeV  =>  S^2: M_6 > {:.3f} TeV.".format(
        excl["M_6_excl_add"], excl["M_6_excl_s2"]))

    # Assess survival
    if excl["M_6_excl_s2"] < excl["M_6_excl_add"]:
        loosening = excl["M_6_excl_add"] - excl["M_6_excl_s2"]
        pct = 100.0 * loosening / excl["M_6_excl_add"]
        lines.append("  e) The S^2 bound is loosened by {:.3f} TeV ({:.1f}%) "
                      "relative to ADD.".format(loosening, pct))
    else:
        lines.append("  e) The S^2 bound is tightened relative to ADD.")

    lines.append("  f) The 5 TeV ADD exclusion does NOT apply directly to Alpha Ladder.")
    lines.append("     The correct S^2 exclusion is M_6 > {:.3f} TeV.".format(
        excl["M_6_excl_s2"]))

    lines.append("")
    lines.append("=" * 72)

    report = "\n".join(lines)
    print(report)

    return {
        "report": report,
        "scan": scan,
        "exclusion": excl,
        "spectrum_table": table,
        "R_ratio_analytic": R_ratio_analytic,
    }


# ===================================================================
# __main__
# ===================================================================
if __name__ == "__main__":
    summarize_lhc_comparison()
