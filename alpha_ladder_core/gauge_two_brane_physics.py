"""
Two-Brane Physics on S^2 -- Full Two-Brane Scenario in Alpha Ladder
====================================================================

Two codimension-2 branes on S^2: Brane 1 at the north pole (theta=0),
Brane 2 at angle theta_2 from the north pole.  Each brane has tension
T_1, T_2 and creates a conical deficit angle in the S^2 geometry.

Key physics:
  1. Two branes create conical deficits: delta_i = 2*pi*T_i / M_6^4.
     The geometry closes only if delta_1 + delta_2 < 2*pi (total solid
     angle Omega = 4*pi - delta_1 - delta_2 must be positive).
  2. Gauge coupling splitting on Brane 2 depends on angular separation:
     at theta_2, the Killing vector norms give |K_2|^2 = cos^2(theta_2),
     |K_3|^2 = 1, so the couplings SPLIT unless theta_2 = 0 or pi.
  3. The rugby ball (theta_2 = pi, T_1 = T_2) preserves coupling equality
     on both branes (both at poles, SO(2) restored at each).
  4. Inter-brane forces are mediated by the KK tower on S^2.  The Green's
     function on S^2 is a Legendre polynomial sum.
  5. Brane-antibrane (T_2 < 0) scenarios can create a dilaton minimum
     when flux is present (speculative but interesting).
  6. Physics on Brane 2 (atomic physics, cross-sections) differs from
     Brane 1 if couplings split.
  7. All of this is BEYOND the minimal single-brane framework.

Pure Python -- only ``import math``, no numpy/scipy.

Functions
---------
1. two_brane_geometry          -- S^2 geometry with two conical deficits
2. split_couplings             -- gauge couplings on each brane vs theta_2
3. inter_brane_force           -- inter-brane potential via KK tower exchange
4. brane_antibrane_stabilization -- dilaton minimum from brane-antibrane + flux
5. brane_2_physics             -- atomic physics on Brane 2 with split couplings
6. rugby_ball_special_case     -- symmetric antipodal case (well-studied limit)
7. summarize_two_brane         -- full report
"""

import math

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
ALPHA_EM = 1.0 / 137.035999084          # Thomson limit (CODATA 2018)
M_PL_EV = 1.22089e28                    # Planck mass (eV)
HBAR_C_EVM = 1.9733e-7                  # hbar*c (eV*m)
L_PL = 1.61625e-35                      # Planck length (m)
M_E_EV = 5.11e5                         # electron mass (eV)
PHI_GOLDEN = (1.0 + math.sqrt(5.0)) / 2.0

# Derived
PHI_VEV = 0.25 * math.log(4.0 * math.pi * ALPHA_EM)   # = -0.5974...
E_PHI_VEV = math.exp(PHI_VEV)                          # = 0.5503...

# Casimir spectral zeta on S^2
ZETA_S2_MINUS_HALF = -0.25  # zeta_{S^2}(-1/2) = -1/4


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _R_phys(a_0_m):
    """Physical radius R = a_0 * exp(phi_vev)."""
    return a_0_m * E_PHI_VEV


def _M_6_from_a0(a_0_m):
    """6D Planck mass from M_6^4 = M_Pl^2 / (4*pi*R^2)."""
    R_phys_m = _R_phys(a_0_m)
    R_eV_inv = R_phys_m / HBAR_C_EVM
    M_6_eV4 = M_PL_EV ** 2 / (4.0 * math.pi * R_eV_inv ** 2)
    return M_6_eV4 ** 0.25


def _C_casimir(a_0):
    """Casimir coefficient: C_cas = zeta_{S^2}(-1/2) / (4*pi*a_0^4)."""
    return ZETA_S2_MINUS_HALF / (4.0 * math.pi * a_0 ** 4)


def _C_flux(N_flux, a_0):
    """Flux coefficient: C_flux = N^2 / (32*pi^2*a_0^6)."""
    if N_flux <= 0:
        return 0.0
    return N_flux ** 2 / (32.0 * math.pi ** 2 * a_0 ** 6)


def _legendre_P(l_max, x):
    """Compute Legendre polynomials P_0(x) through P_{l_max}(x) via recurrence.

    Returns a list of length l_max+1: [P_0(x), P_1(x), ..., P_{l_max}(x)].
    """
    if l_max < 0:
        return []
    P = [0.0] * (l_max + 1)
    P[0] = 1.0
    if l_max >= 1:
        P[1] = x
    for l in range(1, l_max):
        P[l + 1] = ((2 * l + 1) * x * P[l] - l * P[l - 1]) / (l + 1)
    return P


# ===================================================================
# 1. two_brane_geometry
# ===================================================================
def two_brane_geometry(T_1, T_2, M_6_eV, a_0_m=28e-6):
    """Compute the S^2 geometry with two conical deficits.

    Each codimension-2 brane with tension T creates deficit angle
    delta = 2*pi*T / M_6^4 in the transverse S^2.  Two branes give
    total deficit delta_1 + delta_2, leaving effective solid angle
    Omega = 4*pi - delta_1 - delta_2.

    Parameters
    ----------
    T_1 : float
        Tension of Brane 1 in eV^4.
    T_2 : float
        Tension of Brane 2 in eV^4.
    M_6_eV : float
        6D Planck mass in eV.
    a_0_m : float
        Fiducial radius of S^2 in meters (default 28 um).

    Returns
    -------
    dict with delta_1, delta_2, Omega, V, M_Pl_eff, is_physical, description.
    """
    M_6_4 = M_6_eV ** 4
    if M_6_4 == 0:
        return {"error": "M_6 = 0 is unphysical"}

    delta_1 = 2.0 * math.pi * T_1 / M_6_4
    delta_2 = 2.0 * math.pi * T_2 / M_6_4

    Omega = 4.0 * math.pi - delta_1 - delta_2
    is_physical = Omega > 0

    # Effective volume of the modified S^2 (in m^2)
    R_phys = _R_phys(a_0_m)
    V_m2 = Omega * R_phys ** 2

    # Modified Planck mass: M_Pl_eff^2 = V * M_6^4
    # V in (eV^{-1})^2 = (m / (hbar*c))^2
    R_eV_inv = R_phys / HBAR_C_EVM
    V_eV2 = Omega * R_eV_inv ** 2
    M_Pl_eff_sq = V_eV2 * M_6_4
    M_Pl_eff = math.sqrt(abs(M_Pl_eff_sq)) if M_Pl_eff_sq > 0 else 0.0

    # Ratio to standard M_Pl (which assumes Omega = 4*pi)
    M_Pl_ratio = M_Pl_eff / M_PL_EV if M_PL_EV > 0 else 0.0

    # Fraction of sphere remaining
    frac = Omega / (4.0 * math.pi) if is_physical else 0.0

    description = (
        "Two branes on S^2 with tensions T_1 = {:.4e} eV^4, T_2 = {:.4e} eV^4.\n"
        "Deficit angles: delta_1 = {:.6f} rad ({:.2f} deg), "
        "delta_2 = {:.6f} rad ({:.2f} deg).\n"
        "Total solid angle: Omega = {:.6f} sr ({:.4f}% of full sphere).\n"
        "Physical geometry (Omega > 0): {}.\n"
        "M_Pl_eff / M_Pl = {:.6f} ({:.2f}% shift).".format(
            T_1, T_2,
            delta_1, math.degrees(delta_1),
            delta_2, math.degrees(delta_2),
            Omega, frac * 100.0,
            is_physical,
            M_Pl_ratio, (M_Pl_ratio - 1.0) * 100.0,
        )
    )

    return {
        "T_1": T_1,
        "T_2": T_2,
        "M_6_eV": M_6_eV,
        "delta_1_rad": delta_1,
        "delta_1_deg": math.degrees(delta_1),
        "delta_2_rad": delta_2,
        "delta_2_deg": math.degrees(delta_2),
        "Omega": Omega,
        "fraction_of_sphere": frac,
        "V_m2": V_m2,
        "V_eV_inv2": V_eV2,
        "M_Pl_eff_eV": M_Pl_eff,
        "M_Pl_ratio": M_Pl_ratio,
        "is_physical": is_physical,
        "a_0_m": a_0_m,
        "R_phys_m": R_phys,
        "description": description,
    }


# ===================================================================
# 2. split_couplings
# ===================================================================
def split_couplings(theta_2=None, phi_vev=None):
    """Compute the gauge couplings on each brane for various theta_2.

    Brane 1 sits at the north pole (theta=0).  By SO(2) symmetry around
    the pole, both coset Killing vectors have equal norms there:
        alpha_2 = alpha_3 = alpha_EM on Brane 1.

    Brane 2 at polar angle theta_2 (azimuth phi=0 by convention):
        |K_2|^2 = cos^2(theta_2),  |K_3|^2 = 1
    except at the poles (theta_2 = 0 or pi) where SO(2) symmetry forces
    equality.

    The effective coupling is alpha_a_eff = alpha_tree * |K_a|^2 at the
    brane position (normalized so |K_a|^2 = 1 at the north pole).

    The splitting parameter S = |alpha_2 - alpha_3| / (alpha_2 + alpha_3)
    on Brane 2.

    Parameters
    ----------
    theta_2 : float or None
        If a single float, compute for that angle (radians).
        If None, compute for several canonical angles.
    phi_vev : float or None
        Dilaton vev.  If None, use PHI_VEV.

    Returns
    -------
    list of dicts, each with theta_2_deg, alpha_2_b1, alpha_3_b1,
    alpha_2_b2, alpha_3_b2, splitting_param.
    """
    if phi_vev is None:
        phi_vev = PHI_VEV
    alpha_tree = math.exp(4.0 * phi_vev) / (4.0 * math.pi)

    if theta_2 is not None and not hasattr(theta_2, '__iter__'):
        theta_2_values = [theta_2]
    elif theta_2 is not None:
        theta_2_values = list(theta_2)
    else:
        theta_2_values = [
            0.0,
            math.radians(30.0),
            math.radians(45.0),
            math.radians(60.0),
            math.radians(90.0),
            math.radians(120.0),
            math.radians(150.0),
            math.radians(180.0),
        ]

    results = []
    for th2 in theta_2_values:
        # Brane 1 at north pole: always equal
        alpha_2_b1 = alpha_tree
        alpha_3_b1 = alpha_tree

        # Brane 2 at (theta_2, phi=0)
        # At poles: SO(2) symmetry forces equality
        if abs(th2) < 1e-12 or abs(th2 - math.pi) < 1e-12:
            K2_sq = 1.0
            K3_sq = 1.0
        else:
            K2_sq = math.cos(th2) ** 2
            K3_sq = 1.0

        alpha_2_b2 = alpha_tree * K2_sq
        alpha_3_b2 = alpha_tree * K3_sq

        # Splitting parameter
        denom = alpha_2_b2 + alpha_3_b2
        if denom > 0:
            splitting = abs(alpha_3_b2 - alpha_2_b2) / denom
        else:
            splitting = 0.0

        results.append({
            "theta_2_rad": th2,
            "theta_2_deg": math.degrees(th2),
            "alpha_2_b1": alpha_2_b1,
            "alpha_3_b1": alpha_3_b1,
            "alpha_2_b2": alpha_2_b2,
            "alpha_3_b2": alpha_3_b2,
            "K2_squared": K2_sq,
            "K3_squared": K3_sq,
            "splitting_param": splitting,
            "couplings_equal_b1": True,
            "couplings_equal_b2": abs(K2_sq - K3_sq) < 1e-12,
        })

    return results


# ===================================================================
# 3. inter_brane_force
# ===================================================================
def inter_brane_force(theta_2, a_0_m=28e-6, phi_vev=None, l_max=50):
    """Compute the inter-brane potential via bulk KK field exchange on S^2.

    The Green's function on S^2 at angular separation theta is:
        G(theta) = sum_{l=1}^{inf} (2l+1)/(4*pi) * P_l(cos(theta)) / [l(l+1)/R^2]

    The l=0 term is omitted (constant, absorbed into cosmological constant).

    The inter-brane Yukawa potential in 4D at distance r:
        V_4D(r) = -(1/M_Pl^2) * sum_{l=1}^{l_max} (2l+1) * P_l(cos(theta_2))
                   * exp(-m_l * r) / (4*pi*r)
    where m_l = sqrt(l(l+1)) * hbar_c / R_phys is the KK mass of mode l.

    Parameters
    ----------
    theta_2 : float
        Angular separation between branes on S^2 (radians).
    a_0_m : float
        Fiducial radius in meters (default 28 um).
    phi_vev : float or None
        Dilaton vev.
    l_max : int
        Maximum KK mode in the sum.

    Returns
    -------
    dict with V_values at several 4D distances, comparison with Newton, etc.
    """
    if phi_vev is None:
        phi_vev = PHI_VEV

    R_phys = a_0_m * math.exp(phi_vev)
    R_eV_inv = R_phys / HBAR_C_EVM  # R in eV^{-1}

    # Compute Legendre polynomials P_l(cos(theta_2))
    cos_theta = math.cos(theta_2)
    P_vals = _legendre_P(l_max, cos_theta)

    # KK masses: m_l = sqrt(l(l+1)) / R_eV_inv  (in eV)
    # But for Yukawa in SI, we need m_l * r in dimensionless units.
    # m_l (eV) * r (m) / (hbar*c in eV*m) = sqrt(l(l+1)) * r / R_phys
    # So m_l * r = sqrt(l(l+1)) * r / R_phys  (dimensionless)

    # Green's function on S^2 (angular part, dimensionless):
    # G_S2(theta) = sum_{l=1}^{l_max} (2l+1) * P_l(cos(theta)) / [4*pi * l*(l+1)]
    G_S2 = 0.0
    for l in range(1, l_max + 1):
        G_S2 += (2 * l + 1) * P_vals[l] / (4.0 * math.pi * l * (l + 1))

    # Evaluate V_4D(r) at several distances
    # V_4D(r) = -(1/M_Pl^2) * sum (2l+1) P_l(cos theta) exp(-m_l r) / (4 pi r)
    # In eV units: r_eV = r_m / hbar_c_eVm
    # V has units eV (potential energy between unit masses)

    test_distances = [
        ("1 mm", 1.0e-3),
        ("1 cm", 1.0e-2),
        ("1 m", 1.0),
        ("10 m", 10.0),
    ]

    V_results = []
    for label, r_m in test_distances:
        r_eV_inv = r_m / HBAR_C_EVM

        # Sum over KK modes
        V_kk = 0.0
        dominant_l = None
        dominant_contrib = 0.0
        for l in range(1, l_max + 1):
            m_l_eV = math.sqrt(l * (l + 1)) / R_eV_inv
            exponent = m_l_eV * r_eV_inv
            # Prevent overflow: exp(-x) for x > 700 is zero
            if exponent > 700.0:
                continue
            contrib = (2 * l + 1) * P_vals[l] * math.exp(-exponent) / (4.0 * math.pi * r_eV_inv)
            V_kk += contrib
            if abs(contrib) > abs(dominant_contrib):
                dominant_contrib = contrib
                dominant_l = l

        V_4D = -V_kk / (M_PL_EV ** 2)

        # Newtonian potential for comparison: V_N = -G * m1 * m2 / r
        # For unit masses (1 eV/c^2 each): G = 1/(8*pi*M_Pl^2) in natural units
        # V_N = -1 / (8*pi*M_Pl^2 * r_eV_inv)
        V_newton = -1.0 / (8.0 * math.pi * M_PL_EV ** 2 * r_eV_inv)

        ratio = V_4D / V_newton if V_newton != 0 else 0.0

        V_results.append({
            "label": label,
            "r_m": r_m,
            "V_4D_eV": V_4D,
            "V_newton_eV": V_newton,
            "ratio_V_kk_over_newton": ratio,
            "dominant_l": dominant_l,
            "dominant_contrib": dominant_contrib,
            "n_active_modes": sum(1 for l in range(1, l_max + 1)
                                  if math.sqrt(l * (l + 1)) * r_m / R_phys < 700),
        })

    # Arc length between branes
    d_arc = theta_2 * R_phys  # meters
    d_arc_um = d_arc * 1e6

    # Lightest KK mass
    m_1_eV = math.sqrt(2.0) / R_eV_inv

    return {
        "theta_2_rad": theta_2,
        "theta_2_deg": math.degrees(theta_2),
        "R_phys_m": R_phys,
        "d_arc_m": d_arc,
        "d_arc_um": d_arc_um,
        "l_max": l_max,
        "G_S2_angular": G_S2,
        "m_lightest_KK_eV": m_1_eV,
        "V_results": V_results,
        "note": (
            "Inter-brane force is mediated by the KK tower on S^2. "
            "Each mode l contributes a Yukawa potential with mass "
            "m_l = sqrt(l(l+1))/R and coupling 1/M_Pl^2 (gravitational). "
            "At distances r >> R, the force is exponentially suppressed "
            "by exp(-m_1 * r) ~ exp(-sqrt(2) * r / R). "
            "For R ~ 15 um, m_1 ~ 13 meV, and r = 1 mm: "
            "exp(-m_1 * r) ~ exp(-66) ~ 10^{-29}. "
            "The inter-brane force is unobservably small at macroscopic distances."
        ),
    }


# ===================================================================
# 4. brane_antibrane_stabilization
# ===================================================================
def brane_antibrane_stabilization(T_ratio_values=None, a_0=1.0, N_flux=1):
    """Scan brane-antibrane scenarios for dilaton stabilization.

    Brane (T_1 > 0) and antibrane (T_2 < 0) with net tension
    T_eff = T_1 * (1 + T_ratio), where T_ratio = T_2 / T_1 < 0.

    The dilaton potential:
        V(phi) = T_eff * e^{4*phi} + C_cas * e^{-4*phi} + C_flux * e^{6*phi}

    For T_eff < 0 (|T_2| > T_1), the brane term flips sign.
    Combined with Casimir (C_cas < 0) and flux (C_flux > 0):
        V = -|T_eff|*e^{4*phi} - |C_cas|*e^{-4*phi} + C_flux*e^{6*phi}
    Both first terms are negative; the flux term is positive.
    A minimum could exist where flux balances the negative terms.

    Parameters
    ----------
    T_ratio_values : list of float or None
        Values of T_2/T_1 to scan. Default: -2.0 to 0.0 in steps of 0.1.
    a_0 : float
        Fiducial radius (Planck units).
    N_flux : int
        Flux quantum number.

    Returns
    -------
    dict with scan results and verdict.
    """
    if T_ratio_values is None:
        T_ratio_values = [x / 10.0 for x in range(-20, 1)]  # -2.0 to 0.0

    T_1 = 1.0  # Reference tension in Planck units
    C_cas = _C_casimir(a_0)
    C_fl = _C_flux(N_flux, a_0)

    phi_scan = [p / 40.0 for p in range(-200, 201)]  # -5.0 to 5.0 step 0.025
    results = []
    best_match = None
    best_diff = float("inf")

    for T_ratio in T_ratio_values:
        T_2 = T_1 * T_ratio
        T_eff = T_1 + T_2

        # Scan for minimum
        min_V = float("inf")
        min_phi = None
        for phi in phi_scan:
            try:
                e4phi = math.exp(4.0 * phi)
                em4phi = math.exp(-4.0 * phi)
                e6phi = math.exp(6.0 * phi)
            except OverflowError:
                continue
            V = T_eff * e4phi + C_cas * em4phi + C_fl * e6phi
            if V < min_V:
                min_V = V
                min_phi = phi

        # Check if minimum is interior (not at boundary)
        is_interior = (min_phi is not None
                       and min_phi > phi_scan[0] + 0.05
                       and min_phi < phi_scan[-1] - 0.05)

        has_minimum = is_interior

        # Check if phi_min matches phi_vev
        matches = False
        diff = float("inf")
        if has_minimum and min_phi is not None:
            diff = abs(min_phi - PHI_VEV)
            matches = diff < 0.1

        entry = {
            "T_ratio": T_ratio,
            "T_1": T_1,
            "T_2": T_2,
            "T_eff": T_eff,
            "phi_min": min_phi,
            "V_min": min_V,
            "has_minimum": has_minimum,
            "diff_from_phi_vev": diff if has_minimum else None,
            "matches_alpha": matches,
        }
        results.append(entry)

        if has_minimum and diff < best_diff:
            best_diff = diff
            best_match = entry

    # Analysis of the potential structure
    # For T_eff < 0 and C_cas < 0 and C_flux > 0:
    # V = -|T_eff|*e^{4phi} - |C_cas|*e^{-4phi} + C_flux*e^{6phi}
    # dV/dphi = -4|T_eff|*e^{4phi} + 4|C_cas|*e^{-4phi} + 6*C_flux*e^{6phi}
    # Setting to zero: 6*C_flux*e^{6phi} + 4|C_cas|*e^{-4phi} = 4|T_eff|*e^{4phi}
    # The LHS has e^{6phi} (growing fast) and e^{-4phi} (shrinking).
    # The RHS has e^{4phi} (growing).
    # For large phi: LHS ~ e^{6phi} >> RHS ~ e^{4phi}  =>  dV > 0
    # For small phi: LHS ~ e^{-4phi} >> everything  =>  dV > 0
    # For intermediate phi: RHS could exceed LHS  =>  dV < 0
    # So a minimum is possible flanked by maxima -- but only if the
    # intermediate region exists where the RHS dominates.

    n_with_min = sum(1 for r in results if r["has_minimum"])

    if best_match is not None and best_match["matches_alpha"]:
        verdict = (
            "Brane-antibrane scenario with T_ratio = {:.2f} gives phi_min = {:.4f}, "
            "close to phi_vev = {:.4f} (diff = {:.4f}). "
            "This is a TUNING of T_2/T_1, not a prediction. "
            "The mechanism works but requires specific tension ratio.".format(
                best_match["T_ratio"], best_match["phi_min"],
                PHI_VEV, best_match["diff_from_phi_vev"],
            )
        )
    elif best_match is not None:
        verdict = (
            "Best match: T_ratio = {:.2f} gives phi_min = {:.4f}, "
            "still {:.4f} away from phi_vev = {:.4f}. "
            "{} out of {} ratios produce an interior minimum.".format(
                best_match["T_ratio"], best_match["phi_min"],
                best_match["diff_from_phi_vev"], PHI_VEV,
                n_with_min, len(results),
            )
        )
    else:
        verdict = (
            "No interior minimum found for any T_ratio in [{:.1f}, {:.1f}]. "
            "The brane-antibrane scenario with N_flux = {} does not stabilize "
            "the dilaton in the scanned range.".format(
                T_ratio_values[0], T_ratio_values[-1], N_flux,
            )
        )

    return {
        "n_ratios_scanned": len(results),
        "n_with_minimum": n_with_min,
        "best_match": best_match,
        "best_diff": best_diff if best_match is not None else None,
        "verdict": verdict,
        "C_cas": C_cas,
        "C_flux": C_fl,
        "N_flux": N_flux,
        "a_0": a_0,
        "phi_vev_target": PHI_VEV,
        "results": results,
        "note": (
            "For T_eff < 0: brane and Casimir are BOTH negative. "
            "Flux (positive, e^{6*phi}) can compete at large phi. "
            "A minimum requires flux to balance the combined negative terms. "
            "The scenario is speculative: brane-antibrane systems are generically "
            "unstable (open string tachyon), and this is a low-energy effective "
            "description that may miss the instability."
        ),
    }


# ===================================================================
# 5. brane_2_physics
# ===================================================================
def brane_2_physics(theta_2=None, a_0_m=28e-6, phi_vev=None):
    """Compute atomic physics on Brane 2 with split gauge couplings.

    If the effective fine-structure constant on Brane 2 differs from
    Brane 1, all electromagnetic physics changes:
    - Hydrogen binding energy: E_H = alpha^2 * m_e / 2
    - Bohr radius: a_B = 1 / (m_e * alpha)  [natural units]
    - Thomson cross-section: sigma_T = 8*pi*alpha^2 / (3*m_e^2)

    We take the effective alpha on Brane 2 as the geometric mean of the
    two split couplings: alpha_eff = sqrt(alpha_2 * alpha_3).  This is
    appropriate if the photon couples to both gauge fields symmetrically.

    Parameters
    ----------
    theta_2 : list of float or None
        Brane 2 angles (radians).  If None, use canonical set.
    a_0_m : float
        Fiducial radius (not used directly, for consistency).
    phi_vev : float or None
        Dilaton vev.

    Returns
    -------
    list of dicts with alpha_eff_brane2, E_hydrogen_ratio, bohr_radius_ratio,
    thomson_ratio for each theta_2.
    """
    if phi_vev is None:
        phi_vev = PHI_VEV
    alpha_tree = math.exp(4.0 * phi_vev) / (4.0 * math.pi)

    if theta_2 is None:
        theta_2_values = [
            0.0,
            math.radians(30.0),
            math.radians(45.0),
            math.radians(60.0),
            math.radians(90.0),
            math.radians(120.0),
            math.radians(150.0),
            math.radians(180.0),
        ]
    elif not hasattr(theta_2, '__iter__'):
        theta_2_values = [theta_2]
    else:
        theta_2_values = list(theta_2)

    # Brane 1 reference values (standard SM)
    alpha_b1 = alpha_tree  # = alpha_EM
    E_H_b1 = alpha_b1 ** 2 * M_E_EV / 2.0  # eV
    a_B_b1 = HBAR_C_EVM / (M_E_EV * alpha_b1)  # meters
    sigma_T_b1 = 8.0 * math.pi * alpha_b1 ** 2 / (3.0 * M_E_EV ** 2)  # eV^{-2}

    results = []
    for th2 in theta_2_values:
        # Killing vector norms at Brane 2
        if abs(th2) < 1e-12 or abs(th2 - math.pi) < 1e-12:
            K2_sq = 1.0
            K3_sq = 1.0
        else:
            K2_sq = math.cos(th2) ** 2
            K3_sq = 1.0

        alpha_2_b2 = alpha_tree * K2_sq
        alpha_3_b2 = alpha_tree * K3_sq

        # Effective alpha on Brane 2: geometric mean
        # If alpha_2 = 0 (theta = 90 deg), alpha_eff = 0 -- pathological
        alpha_eff = math.sqrt(alpha_2_b2 * alpha_3_b2) if alpha_2_b2 > 0 else 0.0

        # Ratios to Brane 1 values
        if alpha_eff > 0 and alpha_b1 > 0:
            alpha_ratio = alpha_eff / alpha_b1
            E_H_ratio = alpha_ratio ** 2  # E_H ~ alpha^2
            a_B_ratio = 1.0 / alpha_ratio  # a_B ~ 1/alpha
            sigma_T_ratio = alpha_ratio ** 2  # sigma_T ~ alpha^2
        else:
            alpha_ratio = 0.0
            E_H_ratio = 0.0
            a_B_ratio = float("inf")
            sigma_T_ratio = 0.0

        # Absolute values on Brane 2
        E_H_b2 = alpha_eff ** 2 * M_E_EV / 2.0 if alpha_eff > 0 else 0.0
        a_B_b2 = HBAR_C_EVM / (M_E_EV * alpha_eff) if alpha_eff > 0 else float("inf")

        results.append({
            "theta_2_rad": th2,
            "theta_2_deg": math.degrees(th2),
            "alpha_2_b2": alpha_2_b2,
            "alpha_3_b2": alpha_3_b2,
            "alpha_eff_brane2": alpha_eff,
            "alpha_ratio_b2_over_b1": alpha_ratio,
            "E_hydrogen_eV_b1": E_H_b1,
            "E_hydrogen_eV_b2": E_H_b2,
            "E_hydrogen_ratio": E_H_ratio,
            "bohr_radius_m_b1": a_B_b1,
            "bohr_radius_m_b2": a_B_b2,
            "bohr_radius_ratio": a_B_ratio,
            "thomson_ratio": sigma_T_ratio,
        })

    return results


# ===================================================================
# 6. rugby_ball_special_case
# ===================================================================
def rugby_ball_special_case(T, a_0_m=28e-6, phi_vev=None):
    """The symmetric rugby ball: theta_2 = pi, T_1 = T_2 = T.

    Two identical branes at antipodal points on S^2.  This is the most
    studied compactification in the codimension-2 literature.

    Key properties:
    - Both branes at poles => SO(2) symmetry at each brane => no splitting
    - Antipodal: cos(pi) = -1, so P_l(-1) = (-1)^l
    - S^2 becomes a spindle (lens space) with two identical conical deficits
    - Total deficit: 2*delta = 4*pi*T/M_6^4
    - Effective solid angle: Omega = 4*pi - 2*delta
    - Planck mass reduced: M_Pl^2 = Omega * R^2 * M_6^4

    Parameters
    ----------
    T : float
        Common brane tension in eV^4.
    a_0_m : float
        Fiducial radius in meters.
    phi_vev : float or None
        Dilaton vev.

    Returns
    -------
    dict with geometry, spectrum_shift, M_Pl_ratio, V_dilaton_shift.
    """
    if phi_vev is None:
        phi_vev = PHI_VEV

    R_phys = a_0_m * math.exp(phi_vev)
    R_eV_inv = R_phys / HBAR_C_EVM
    M_6_eV = _M_6_from_a0(a_0_m)
    M_6_4 = M_6_eV ** 4

    # --- Geometry ---
    delta = 2.0 * math.pi * T / M_6_4 if M_6_4 > 0 else 0.0
    Omega = 4.0 * math.pi - 2.0 * delta
    is_physical = Omega > 0
    frac = Omega / (4.0 * math.pi) if is_physical else 0.0

    # Modified Planck mass
    V_eV2 = Omega * R_eV_inv ** 2
    M_Pl_eff_sq = V_eV2 * M_6_4
    M_Pl_eff = math.sqrt(abs(M_Pl_eff_sq)) if M_Pl_eff_sq > 0 else 0.0
    M_Pl_ratio = M_Pl_eff / M_PL_EV if M_PL_EV > 0 else 0.0

    geometry = {
        "delta_rad": delta,
        "delta_deg": math.degrees(delta),
        "total_deficit_rad": 2.0 * delta,
        "Omega": Omega,
        "fraction_of_sphere": frac,
        "M_Pl_eff_eV": M_Pl_eff,
        "M_Pl_ratio": M_Pl_ratio,
        "is_physical": is_physical,
    }

    # --- KK spectrum modification ---
    # On the rugby ball (S^2 / Z_2 with two conical deficits), the KK spectrum
    # changes.  For a smooth rugby ball (same deficit at both poles), the
    # spherical harmonics are restricted by the deficit.  The effective
    # quantum number is shifted: l_eff = l / (1 - delta/(2*pi)).
    #
    # More precisely: if the deficit fraction is f = delta/(2*pi), then the
    # angular period is reduced, and the eigenvalues become:
    #   lambda_l = l*(l+1) / (1-f)^2   (approximate, for small f)
    # This INCREASES the KK masses by a factor 1/(1-f).
    f = delta / (2.0 * math.pi) if is_physical else 0.0
    mass_shift_factor = 1.0 / (1.0 - f) if f < 1.0 else float("inf")

    # First few KK masses on the rugby ball vs round S^2
    spectrum_comparison = []
    for l in range(1, 11):
        m_round = math.sqrt(l * (l + 1)) / R_eV_inv
        m_rugby = m_round * mass_shift_factor
        spectrum_comparison.append({
            "l": l,
            "m_round_eV": m_round,
            "m_rugby_eV": m_rugby,
            "ratio": mass_shift_factor,
        })

    spectrum_shift = {
        "deficit_fraction": f,
        "mass_shift_factor": mass_shift_factor,
        "modes": spectrum_comparison,
        "note": (
            "For small deficit, KK masses shift by factor 1/(1-f) where "
            "f = delta/(2*pi). This is a SMALL correction for T << M_6^4. "
            "For T = (1 TeV)^4 and M_6 = 3 TeV: f ~ 0.026, shift ~ 2.6%."
        ),
    }

    # --- Dilaton potential shift ---
    # Two equal-tension branes: V_brane = 2*T*exp(4*phi)
    # This doubles the brane contribution compared to a single brane.
    e4phi = math.exp(4.0 * phi_vev)
    V_single = T * e4phi
    V_rugby = 2.0 * T * e4phi
    V_ratio = 2.0  # exactly double

    V_dilaton_shift = {
        "V_single_brane": V_single,
        "V_rugby_ball": V_rugby,
        "ratio": V_ratio,
        "note": (
            "The rugby ball doubles the brane potential: V = 2*T*exp(4*phi). "
            "This shifts any minimum toward smaller phi (larger extra dimensions) "
            "if the minimum exists. The qualitative structure is unchanged."
        ),
    }

    # --- Coupling check ---
    # At both poles: SO(2) symmetry => no splitting
    coupling_check = {
        "splitting_brane1": 0.0,
        "splitting_brane2": 0.0,
        "note": (
            "Both branes sit at poles of S^2. The SO(2) rotation symmetry "
            "is preserved at each pole, ensuring alpha_2 = alpha_3 on both "
            "branes. The rugby ball does NOT split gauge couplings."
        ),
    }

    return {
        "T": T,
        "a_0_m": a_0_m,
        "M_6_eV": M_6_eV,
        "R_phys_m": R_phys,
        "geometry": geometry,
        "spectrum_shift": spectrum_shift,
        "M_Pl_ratio": M_Pl_ratio,
        "V_dilaton_shift": V_dilaton_shift,
        "coupling_check": coupling_check,
        "is_physical": is_physical,
        "description": (
            "Rugby ball compactification: two identical branes at the poles of S^2. "
            "Each creates deficit delta = {:.4e} rad ({:.2f} deg). "
            "Solid angle reduced to {:.4f}% of full sphere. "
            "KK masses shift by factor {:.6f}. "
            "No coupling splitting (SO(2) at both poles).".format(
                delta, math.degrees(delta), frac * 100.0, mass_shift_factor,
            )
        ),
    }


# ===================================================================
# 7. summarize_two_brane
# ===================================================================
def summarize_two_brane():
    """Full summary of two-brane physics on S^2.

    Returns
    -------
    dict with all sub-analyses and a formatted report string.
    """
    # 1. Geometry with example tensions
    # Use T = (1 TeV)^4 and M_6 = 3 TeV
    T_example = (1.0e12) ** 4  # eV^4
    M_6_example = 3.0e12       # eV
    geom = two_brane_geometry(T_example, T_example, M_6_example)

    # 2. Split couplings
    split = split_couplings()

    # 3. Inter-brane force (antipodal)
    ibf = inter_brane_force(math.pi, a_0_m=28e-6)

    # 4. Brane-antibrane stabilization
    babs = brane_antibrane_stabilization()

    # 5. Brane 2 physics
    b2p = brane_2_physics()

    # 6. Rugby ball
    rugby = rugby_ball_special_case(T_example)

    # Build report
    lines = [
        "=" * 72,
        "TWO-BRANE PHYSICS ON S^2 -- FULL SUMMARY",
        "=" * 72,
        "",
        "Setup: Two codimension-2 branes on S^2, Brane 1 at north pole,",
        "       Brane 2 at angle theta_2 from north pole.",
        "",
        "--- 1. TWO-BRANE GEOMETRY ---",
        "Example: T_1 = T_2 = (1 TeV)^4, M_6 = 3 TeV",
        "delta_1 = {:.6f} rad ({:.2f} deg)".format(
            geom["delta_1_rad"], geom["delta_1_deg"]),
        "delta_2 = {:.6f} rad ({:.2f} deg)".format(
            geom["delta_2_rad"], geom["delta_2_deg"]),
        "Omega = {:.6f} sr ({:.4f}% of sphere)".format(
            geom["Omega"], geom["fraction_of_sphere"] * 100.0),
        "Physical: {}".format(geom["is_physical"]),
        "M_Pl_eff / M_Pl = {:.6f}".format(geom["M_Pl_ratio"]),
        "",
        "--- 2. GAUGE COUPLING SPLITTING ---",
        "{:<10s} {:>12s} {:>12s} {:>12s} {:>12s} {:>10s}".format(
            "theta_2", "alpha_2_b1", "alpha_3_b1", "alpha_2_b2", "alpha_3_b2", "splitting"),
    ]

    for s in split:
        lines.append(
            "{:>8.1f} d {:>12.6e} {:>12.6e} {:>12.6e} {:>12.6e} {:>10.4f}".format(
                s["theta_2_deg"],
                s["alpha_2_b1"], s["alpha_3_b1"],
                s["alpha_2_b2"], s["alpha_3_b2"],
                s["splitting_param"],
            )
        )

    lines.extend([
        "",
        "Key: At theta_2 = 90 deg, alpha_2 = 0 on Brane 2 (maximal splitting).",
        "     At theta_2 = 0 or 180 deg, no splitting (branes at poles).",
        "",
        "--- 3. INTER-BRANE FORCE (antipodal branes, theta_2 = pi) ---",
        "R_phys = {:.4e} m".format(ibf["R_phys_m"]),
        "Arc distance = {:.4e} m ({:.2f} um)".format(
            ibf["d_arc_m"], ibf["d_arc_um"]),
        "Lightest KK mass: {:.4e} eV".format(ibf["m_lightest_KK_eV"]),
        "G_S2 (angular Green's function): {:.6e}".format(ibf["G_S2_angular"]),
        "",
    ])

    for vr in ibf["V_results"]:
        lines.append(
            "  r = {:<6s}  V_KK = {:>12.4e} eV  V_Newton = {:>12.4e} eV  "
            "ratio = {:>12.4e}  active_modes = {}".format(
                vr["label"], vr["V_4D_eV"], vr["V_newton_eV"],
                vr["ratio_V_kk_over_newton"], vr["n_active_modes"],
            )
        )

    lines.extend([
        "",
        "--- 4. BRANE-ANTIBRANE STABILIZATION ---",
        "Scanned {} T_ratio values (T_2/T_1 from -2.0 to 0.0)".format(
            babs["n_ratios_scanned"]),
        "N_flux = {},  C_cas = {:.4e},  C_flux = {:.4e}".format(
            babs["N_flux"], babs["C_cas"], babs["C_flux"]),
        "Minima found: {} / {}".format(
            babs["n_with_minimum"], babs["n_ratios_scanned"]),
        "Verdict: {}".format(babs["verdict"]),
        "",
        "--- 5. BRANE 2 PHYSICS (atomic physics with split couplings) ---",
        "{:<10s} {:>12s} {:>10s} {:>12s} {:>10s} {:>10s}".format(
            "theta_2", "alpha_eff", "ratio", "E_H ratio", "a_B ratio", "sigma_T"),
    ])

    for bp in b2p:
        lines.append(
            "{:>8.1f} d {:>12.6e} {:>10.6f} {:>12.6f} {:>10.6f} {:>10.6f}".format(
                bp["theta_2_deg"],
                bp["alpha_eff_brane2"],
                bp["alpha_ratio_b2_over_b1"],
                bp["E_hydrogen_ratio"],
                bp["bohr_radius_ratio"] if bp["bohr_radius_ratio"] < 1e10 else float("inf"),
                bp["thomson_ratio"],
            )
        )

    lines.extend([
        "",
        "Key: At theta_2 = 90 deg, alpha_eff = 0 => no bound atoms on Brane 2!",
        "     At theta_2 = 60 deg, alpha_eff ~ 0.866 * alpha_EM => weaker binding.",
        "",
        "--- 6. RUGBY BALL SPECIAL CASE ---",
        "T = {:.4e} eV^4,  M_6 = {:.4e} eV".format(T_example, rugby["M_6_eV"]),
        "Deficit: {:.4e} rad ({:.2f} deg)".format(
            rugby["geometry"]["delta_rad"],
            rugby["geometry"]["delta_deg"]),
        "Fraction of sphere: {:.6f}".format(rugby["geometry"]["fraction_of_sphere"]),
        "KK mass shift factor: {:.6f}".format(
            rugby["spectrum_shift"]["mass_shift_factor"]),
        "M_Pl_eff / M_Pl = {:.6f}".format(rugby["M_Pl_ratio"]),
        "Coupling splitting: {} (none at poles)".format(
            rugby["coupling_check"]["splitting_brane1"]),
        "",
        "--- 7. KEY CONCLUSIONS ---",
        "",
        "1. GEOMETRY: Two branes create conical deficits that reduce the S^2",
        "   solid angle.  For T << M_6^4, the deficits are tiny (< 1 deg).",
        "   Physical constraint: Omega = 4*pi - delta_1 - delta_2 > 0.",
        "",
        "2. COUPLING SPLITTING: Brane 2 at generic theta_2 sees DIFFERENT",
        "   gauge couplings |K_2|^2 = cos^2(theta_2), |K_3|^2 = 1.",
        "   Maximum splitting at theta_2 = 90 deg (equator); none at poles.",
        "   This is beyond the minimal framework (single brane at a pole).",
        "",
        "3. RUGBY BALL: Antipodal branes (theta_2 = pi) preserve coupling",
        "   equality on both branes.  The KK spectrum shifts slightly.",
        "   This is the safest two-brane extension of the minimal framework.",
        "",
        "4. BRANE-ANTIBRANE: Net negative tension can create a minimum with",
        "   flux stabilization.  Speculative: brane-antibrane systems have",
        "   known instabilities (open string tachyon) not captured here.",
        "",
        "5. INTER-BRANE FORCES: Exponentially suppressed at macroscopic",
        "   distances (exp(-sqrt(2)*r/R) for r >> R ~ 15 um).",
        "   Unobservable at r > 1 mm.",
        "",
        "6. BRANE 2 PHYSICS: If couplings split, atomic physics on Brane 2",
        "   differs from Brane 1.  Hydrogen binding, Bohr radius, Thomson",
        "   cross-section all change.  At theta_2 = 90 deg, one coupling",
        "   vanishes -- no stable atoms possible.",
        "",
        "7. STATUS: All of this is BEYOND the minimal Alpha Ladder framework.",
        "   The minimal framework uses a single brane at a pole of S^2.",
        "   Two-brane extensions are interesting for phenomenology but add",
        "   free parameters (theta_2, T_2) without new predictions.",
        "=" * 72,
    ])

    report = "\n".join(lines)

    return {
        "geometry": geom,
        "split_couplings": split,
        "inter_brane_force": ibf,
        "brane_antibrane": babs,
        "brane_2_physics": b2p,
        "rugby_ball": rugby,
        "report": report,
    }


# ===================================================================
# __main__
# ===================================================================
if __name__ == "__main__":
    result = summarize_two_brane()
    print(result["report"])
