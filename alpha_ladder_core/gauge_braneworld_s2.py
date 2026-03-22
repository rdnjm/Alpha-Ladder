"""
Braneworld Physics on S^2 -- SM Matter on a 4D Brane in 6D Bulk
================================================================

The setup:
  - 6D bulk: Einstein-Hilbert + Gauss-Bonnet gravity, dilaton phi,
    gauge fields from SO(3)/SO(2) isometry of S^2
  - 4D brane: codimension-2 brane located at a point on S^2,
    carrying all Standard Model matter
  - Brane tension T sources the bulk equations and creates a
    conical deficit angle in the transverse S^2

Key physics:
  - Codimension-2 branes create conical deficits (Gibbons-Hawking-Horowitz
    "rugby ball" compactification)
  - Brane tension generates a dilaton potential V_brane = T * exp(4*phi)
  - Combined with Casimir (exp(-4*phi)) and flux (exp(6*phi)), a minimum
    can exist for the dilaton
  - Gap 1 (vacuum polynomial x^2+6x+4=0) is NOT closed: the brane
    potential gives a single phi_min, not a quadratic with sqrt(5) roots

Pure Python -- only ``import math``, no numpy/scipy.

Functions
---------
1. brane_tension_potential  -- V_brane = n_branes * T * exp(4*phi)
2. casimir_potential        -- V_cas = C_cas * exp(-4*phi), C_cas < 0
3. combined_potential       -- V_brane + V_Casimir + V_flux with minimum search
4. find_brane_stabilization -- scan T to find phi_min = phi_vev
5. deficit_angle            -- conical deficit from brane tension
6. brane_localized_coupling -- effective alpha and G on the brane
7. brane_matter_spectrum    -- KK graviton and gauge towers seen from brane
8. does_brane_close_gap1    -- honest assessment: no, it does not
9. summarize_braneworld     -- full report
"""

import math

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
ALPHA_EM = 1.0 / 137.035999084          # Thomson limit (CODATA 2018)
M_PL_EV = 1.22089e28                    # Planck mass (eV)
HBAR_C_EVM = 1.9733e-7                  # hbar*c (eV*m)
L_PL = 1.61625e-35                      # Planck length (m)
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


def _M_6_from_R(R_phys_m):
    """6D Planck mass from M_6^4 = M_Pl^2 / (4*pi*R^2)."""
    R_eV_inv = R_phys_m / HBAR_C_EVM
    M_6_eV4 = M_PL_EV ** 2 / (4.0 * math.pi * R_eV_inv ** 2)
    return M_6_eV4 ** 0.25


def _C_casimir(a_0):
    """Casimir coefficient: C_cas = zeta_{S^2}(-1/2) / (4*pi*a_0^4).

    For S^2: zeta(-1/2) = -1/4, so C_cas = -1/(16*pi*a_0^4) < 0.
    """
    return ZETA_S2_MINUS_HALF / (4.0 * math.pi * a_0 ** 4)


def _C_flux(N_flux, a_0):
    """Flux coefficient: C_flux = N^2 / (32*pi^2*a_0^6)."""
    if N_flux <= 0:
        return 0.0
    return N_flux ** 2 / (32.0 * math.pi ** 2 * a_0 ** 6)


# ===================================================================
# 1. brane_tension_potential
# ===================================================================
def brane_tension_potential(T, phi, n_branes=1):
    """Compute the brane tension contribution to the dilaton potential.

    For a codimension-2 brane in 6D with induced metric e^{2*sigma} g_{4D}:
      S_brane = -T * int d^4x * e^{4*phi} * sqrt(-g_4)
    so the effective potential contribution is:
      V_brane(phi) = n_branes * T * exp(4*phi)

    Parameters
    ----------
    T : float
        Brane tension in Planck units (M_Pl^4).
    phi : float
        Dilaton field value.
    n_branes : int
        Number of branes (1 for single brane, 2 for rugby ball).

    Returns
    -------
    dict with V, dV_dphi, n_branes, T, phi.
    """
    e4phi = math.exp(4.0 * phi)
    V = n_branes * T * e4phi
    dV = 4.0 * n_branes * T * e4phi

    return {
        "V": V,
        "dV_dphi": dV,
        "n_branes": n_branes,
        "T": T,
        "phi": phi,
        "note": ("Brane tension generates V = n*T*exp(4*phi). "
                 "For T > 0, this is monotonically increasing: "
                 "drives phi toward -infinity (shrinking extra dimensions).")
    }


# ===================================================================
# 2. casimir_potential
# ===================================================================
def casimir_potential(phi, a_0=1.0):
    """Bulk Casimir energy as a function of the dilaton.

    The Casimir energy on S^2 depends on the physical radius R = a_0*e^phi:
      V_Casimir ~ 1/R^4 ~ e^{-4*phi}/a_0^4
    The coefficient includes the spectral zeta: zeta_{S^2}(-1/2) = -1/4.
      C_cas = zeta(-1/2) / (4*pi*a_0^4) = -1/(16*pi*a_0^4)

    So V_Casimir = C_cas * exp(-4*phi) < 0 for all phi.

    Parameters
    ----------
    phi : float
        Dilaton field value.
    a_0 : float
        Fiducial radius of S^2 (Planck units).

    Returns
    -------
    dict with V_cas, dV_cas_dphi, C_cas, zeta_value.
    """
    C_cas = _C_casimir(a_0)
    em4phi = math.exp(-4.0 * phi)
    V_cas = C_cas * em4phi
    dV_cas = -4.0 * C_cas * em4phi

    return {
        "V_cas": V_cas,
        "dV_cas_dphi": dV_cas,
        "C_cas": C_cas,
        "zeta_value": ZETA_S2_MINUS_HALF,
        "phi": phi,
        "a_0": a_0,
        "sign": "negative (attractive, drives phi toward +infinity)",
        "note": ("C_cas < 0 because zeta_{S^2}(-1/2) = -1/4. "
                 "The Casimir force favors EXPANDING extra dimensions.")
    }


# ===================================================================
# 3. combined_potential
# ===================================================================
def combined_potential(phi, T, a_0=1.0, N_flux=0):
    """Combined dilaton potential: brane + Casimir + flux.

    V_total = T*exp(4*phi) + C_cas*exp(-4*phi) + C_flux*exp(6*phi)

    The brane drives phi -> -inf, Casimir drives phi -> +inf.
    Flux (if present) drives phi -> -inf more strongly (exp(6*phi)).

    A minimum requires competing terms with different phi-dependence.

    Parameters
    ----------
    phi : float
        Dilaton field value.
    T : float
        Brane tension (Planck units).
    a_0 : float
        Fiducial radius.
    N_flux : int
        Flux quantum number (0 = no flux).

    Returns
    -------
    dict with V, dV_dphi, components, has_minimum info.
    """
    C_cas = _C_casimir(a_0)
    C_fl = _C_flux(N_flux, a_0)

    e4phi = math.exp(4.0 * phi)
    em4phi = math.exp(-4.0 * phi)
    e6phi = math.exp(6.0 * phi)

    V_brane = T * e4phi
    V_casimir = C_cas * em4phi
    V_flux = C_fl * e6phi

    V = V_brane + V_casimir + V_flux
    dV = 4.0 * T * e4phi - 4.0 * C_cas * em4phi + 6.0 * C_fl * e6phi

    # Analytic minimum for brane+Casimir only (no flux):
    # dV/dphi = 4*T*e^{4phi} - 4*C_cas*e^{-4phi} = 0
    # => e^{8phi} = C_cas / T
    # Requires C_cas/T > 0. Since C_cas < 0, need T < 0 for minimum.
    has_minimum_no_flux = False
    phi_min_no_flux = None
    if N_flux == 0 and T != 0.0:
        ratio = C_cas / T
        if ratio > 0:
            has_minimum_no_flux = True
            phi_min_no_flux = math.log(ratio) / 8.0

    return {
        "V": V,
        "dV_dphi": dV,
        "V_brane": V_brane,
        "V_casimir": V_casimir,
        "V_flux": V_flux,
        "C_cas": C_cas,
        "C_flux": C_fl,
        "T": T,
        "phi": phi,
        "N_flux": N_flux,
        "has_minimum_analytic_no_flux": has_minimum_no_flux,
        "phi_min_analytic_no_flux": phi_min_no_flux,
        "note": ("For T > 0 and C_cas < 0 (physical values): "
                 "brane and Casimir push phi in the SAME direction (toward -inf). "
                 "No minimum without flux or negative tension.")
    }


# ===================================================================
# 4. find_brane_stabilization
# ===================================================================
def find_brane_stabilization(a_0=1.0, T_values=None, N_flux=1):
    """Scan brane tension to find phi_min matching phi_vev.

    For each T, compute the combined potential (brane + Casimir + flux)
    and search numerically for a minimum by scanning phi in [-5, 5].

    Parameters
    ----------
    a_0 : float
        Fiducial radius (Planck units).
    T_values : list or None
        Tension values to scan. If None, log-space from 1e-10 to 1e10.
    N_flux : int
        Flux quantum number.

    Returns
    -------
    dict with scan results, best match, and verdict.
    """
    if T_values is None:
        T_values = [10.0 ** (x / 5.0 - 10.0) for x in range(101)]
        # Also include negative tensions
        T_values = ([-t for t in reversed(T_values)] + T_values)

    C_cas = _C_casimir(a_0)
    C_fl = _C_flux(N_flux, a_0)

    phi_scan = [p / 20.0 for p in range(-100, 101)]  # -5.0 to 5.0 step 0.05
    results = []
    best = None
    best_diff = float("inf")

    for T in T_values:
        # Evaluate V at each phi, find minimum
        min_V = float("inf")
        min_phi = None
        for phi in phi_scan:
            e4phi = math.exp(4.0 * phi)
            em4phi = math.exp(-4.0 * phi)
            e6phi = math.exp(6.0 * phi)
            V = T * e4phi + C_cas * em4phi + C_fl * e6phi
            if V < min_V:
                min_V = V
                min_phi = phi

        # Check if minimum is interior (not at boundary)
        is_interior = (min_phi is not None
                       and min_phi > phi_scan[0]
                       and min_phi < phi_scan[-1])

        if is_interior and min_phi is not None:
            diff = abs(min_phi - PHI_VEV)
            matches = diff < 0.05  # within 0.05 of phi_vev
            entry = {
                "T": T,
                "phi_min": min_phi,
                "V_min": min_V,
                "diff_from_phi_vev": diff,
                "matches_alpha": matches,
            }
            results.append(entry)
            if diff < best_diff:
                best_diff = diff
                best = entry

    # Determine verdict
    if best is not None and best["matches_alpha"]:
        verdict = (
            "A brane tension T = {:.4e} gives phi_min = {:.4f}, "
            "close to phi_vev = {:.4f} (diff = {:.4f}). "
            "The brane tension CAN tune the dilaton vev to match alpha_EM, "
            "but this is a tuning, not a prediction.".format(
                best["T"], best["phi_min"], PHI_VEV, best["diff_from_phi_vev"]
            )
        )
    elif best is not None:
        verdict = (
            "Best match: T = {:.4e} gives phi_min = {:.4f}, "
            "still {:.4f} away from phi_vev = {:.4f}. "
            "The brane tension alone does not dynamically select phi_vev.".format(
                best["T"], best["phi_min"], best["diff_from_phi_vev"], PHI_VEV
            )
        )
    else:
        verdict = (
            "No interior minimum found for any T in the scan. "
            "For T > 0 and C_cas < 0, brane and Casimir push in the same direction. "
            "A minimum requires flux (N > 0) or negative tension."
        )

    return {
        "n_T_scanned": len(T_values),
        "n_minima_found": len(results),
        "best_match": best,
        "best_diff": best_diff if best is not None else None,
        "verdict": verdict,
        "phi_vev_target": PHI_VEV,
        "N_flux": N_flux,
        "a_0": a_0,
        "sample_results": results[:10] if len(results) > 10 else results,
    }


# ===================================================================
# 5. deficit_angle
# ===================================================================
def deficit_angle(T, M_6_eV):
    """Conical deficit angle from a codimension-2 brane.

    For a codimension-2 brane, the transverse space develops a conical
    singularity with deficit angle:
      delta = 2*pi * T / M_6^4

    The internal space changes from S^2 (solid angle 4*pi) to a
    "football" with effective solid angle 4*pi - n_branes*delta.

    Parameters
    ----------
    T : float
        Brane tension in units of M_6^4 (dimensionless ratio T/M_6^4
        is what matters; here T is in eV^4 and M_6 in eV).
    M_6_eV : float
        6D Planck mass in eV.

    Returns
    -------
    dict with delta, delta_degrees, physical, effective_solid_angle.
    """
    M_6_4 = M_6_eV ** 4
    if M_6_4 == 0:
        return {"error": "M_6 = 0 is unphysical"}

    delta = 2.0 * math.pi * T / M_6_4
    delta_deg = math.degrees(delta)
    physical = 0 < delta < 2.0 * math.pi

    # For rugby ball (2 branes at poles)
    eff_solid_1 = 4.0 * math.pi - delta          # 1 brane
    eff_solid_2 = 4.0 * math.pi - 2.0 * delta    # 2 branes (rugby ball)
    frac_1 = eff_solid_1 / (4.0 * math.pi)
    frac_2 = eff_solid_2 / (4.0 * math.pi)
    rugby_physical = 0 < 2.0 * delta < 4.0 * math.pi

    return {
        "delta_rad": delta,
        "delta_deg": delta_deg,
        "physical_single": physical,
        "effective_solid_angle_1brane": eff_solid_1,
        "effective_solid_angle_2brane": eff_solid_2,
        "fraction_of_sphere_1brane": frac_1,
        "fraction_of_sphere_2brane": frac_2,
        "rugby_ball_physical": rugby_physical,
        "T": T,
        "M_6_eV": M_6_eV,
        "T_over_M6_4": T / M_6_4,
        "note": ("Deficit angle delta = 2*pi*T/M_6^4. "
                 "Physical range: 0 < delta < 2*pi (single brane), "
                 "0 < 2*delta < 4*pi (rugby ball with equal tensions).")
    }


# ===================================================================
# 6. brane_localized_coupling
# ===================================================================
def brane_localized_coupling(phi_vev=None, a_0=1.0):
    """Effective gauge and gravitational couplings on the brane.

    The brane-localized gauge coupling:
      alpha_brane = exp(4*phi_vev) / (4*pi)
    This is the same as alpha_KK from bulk gauge kinetic term evaluated
    at the brane location. For a round S^2, this is position-independent.

    The effective Newton's constant on the brane:
      G_4D = 1 / (M_Pl^2)  [from dimensional reduction]
    The leading correction from a deficit angle delta is O(delta^2):
      G_brane = G_4D * (1 + c * delta^2 + ...)

    Parameters
    ----------
    phi_vev : float or None
        Dilaton vev. If None, use PHI_VEV (matched to alpha_EM).
    a_0 : float
        Fiducial radius (Planck units).

    Returns
    -------
    dict with alpha_brane, G_brane, corrections.
    """
    if phi_vev is None:
        phi_vev = PHI_VEV

    e4phi = math.exp(4.0 * phi_vev)
    alpha_brane = e4phi / (4.0 * math.pi)

    # Compare with alpha_EM
    alpha_ratio = alpha_brane / ALPHA_EM
    alpha_ppm = (alpha_ratio - 1.0) * 1e6

    # Gravitational coupling (leading order, no deficit correction)
    # G_4D = 1/(8*pi*M_Pl^2) in natural units
    # For small deficit delta, correction is O(delta^2)
    # We cannot compute the exact coefficient without solving the full
    # linearized equations, so we parametrize:
    # G_brane = G_4D * (1 + c_2 * delta^2) with c_2 ~ O(1)
    G_4D_natural = 1.0 / (8.0 * math.pi * M_PL_EV ** 2)

    return {
        "phi_vev": phi_vev,
        "alpha_brane": alpha_brane,
        "alpha_EM": ALPHA_EM,
        "alpha_match_ppm": alpha_ppm,
        "G_4D": G_4D_natural,
        "G_brane_leading": G_4D_natural,
        "deficit_correction": "O(delta^2), parametrically small for small T/M_6^4",
        "note": ("alpha_brane = exp(4*phi_vev)/(4*pi). "
                 "At phi_vev = -0.597, alpha_brane = alpha_EM to < 1 ppm. "
                 "This matching is exact at tree level by construction. "
                 "The gravitational coupling receives O(delta^2) corrections "
                 "from the brane bending, negligible for small deficit angles.")
    }


# ===================================================================
# 7. brane_matter_spectrum
# ===================================================================
def brane_matter_spectrum(phi_vev=None, a_0_m=28e-6):
    """KK spectrum as seen from the brane.

    SM particles live on the brane and have standard 4D masses.
    Bulk KK modes (gravitons, gauge bosons) propagate in the full 6D
    and appear as a tower of massive 4D particles coupling to brane matter.

    Graviton KK tower on S^2:
      m_l = sqrt(l(l+1) - 2) * hbar_c / R_phys,  l = 2, 3, ...
    Gauge boson KK tower:
      m_l = sqrt(l(l+1)) * hbar_c / R_phys,  l = 1, 2, ...

    Parameters
    ----------
    phi_vev : float or None
        Dilaton vev. If None, use PHI_VEV.
    a_0_m : float
        Fiducial radius in meters (default 28 um for Eot-Wash window).

    Returns
    -------
    dict with graviton_modes, gauge_modes, lightest_mode, coupling_strength.
    """
    if phi_vev is None:
        phi_vev = PHI_VEV

    R_phys = a_0_m * math.exp(phi_vev)
    R_phys_eV_inv = R_phys / HBAR_C_EVM
    M_6 = _M_6_from_R(R_phys)

    # Graviton KK: l = 2, 3, ..., 11 (first 10 modes)
    graviton_modes = []
    for l in range(2, 12):
        arg = l * (l + 1) - 2
        if arg > 0:
            m_eV = math.sqrt(arg) * HBAR_C_EVM / R_phys
            # Coupling to brane matter: 1/M_Pl per mode (gravitational strength)
            coupling = 1.0 / M_PL_EV
            graviton_modes.append({
                "l": l,
                "m_eV": m_eV,
                "degeneracy": 2 * l + 1,
                "coupling_eV_inv": coupling,
                "coupling_type": "gravitational (1/M_Pl)",
            })

    # Gauge KK: l = 1, 2, ..., 10 (first 10 modes)
    gauge_modes = []
    for l in range(1, 11):
        m_eV = math.sqrt(l * (l + 1)) * HBAR_C_EVM / R_phys
        # Gauge coupling: g_KK = sqrt(4*pi*alpha_brane)
        alpha_br = math.exp(4.0 * phi_vev) / (4.0 * math.pi)
        g_kk = math.sqrt(4.0 * math.pi * alpha_br)
        gauge_modes.append({
            "l": l,
            "m_eV": m_eV,
            "degeneracy": 2 * l + 1,
            "coupling_g": g_kk,
            "coupling_type": "gauge (g_KK = sqrt(4*pi*alpha))",
        })

    lightest_graviton = graviton_modes[0] if graviton_modes else None
    lightest_gauge = gauge_modes[0] if gauge_modes else None

    # Identify overall lightest bulk mode
    lightest = None
    if lightest_graviton and lightest_gauge:
        if lightest_graviton["m_eV"] < lightest_gauge["m_eV"]:
            lightest = lightest_graviton
        else:
            lightest = lightest_gauge
    elif lightest_graviton:
        lightest = lightest_graviton
    elif lightest_gauge:
        lightest = lightest_gauge

    return {
        "R_phys_m": R_phys,
        "a_0_m": a_0_m,
        "M_6_eV": M_6,
        "graviton_modes": graviton_modes,
        "gauge_modes": gauge_modes,
        "lightest_bulk_mode": lightest,
        "n_graviton_modes": len(graviton_modes),
        "n_gauge_modes": len(gauge_modes),
        "note": ("Graviton KK modes couple at 1/M_Pl (gravitational strength). "
                 "Gauge KK modes couple at g_KK ~ 0.3 (much stronger). "
                 "But gauge KK modes only couple to charged brane matter. "
                 "Lightest graviton: l=2 (spin-2, m ~ sqrt(4)*hbar_c/R). "
                 "Lightest gauge: l=1 (spin-1, m ~ sqrt(2)*hbar_c/R).")
    }


# ===================================================================
# 8. does_brane_close_gap1
# ===================================================================
def does_brane_close_gap1(a_0=1.0, N_flux=1):
    """Assess whether the brane potential closes Gap 1 (vacuum polynomial).

    Gap 1 asks: why does the dilaton minimization produce the polynomial
    x^2 + 6x + 4 = 0, whose root x = -3 + sqrt(5) gives phi_vev?

    The brane potential V = T*exp(4*phi) + C_cas*exp(-4*phi) + C_flux*exp(6*phi).

    Substituting x = exp(2*phi):
      V(x) = T*x^2 + C_cas/x^2 + C_flux*x^3

    dV/dx = 2*T*x - 2*C_cas/x^3 + 3*C_flux*x^2 = 0

    Multiplying by x^3:
      2*T*x^4 + 3*C_flux*x^5 - 2*C_cas = 0

    This is a QUINTIC in x, not a quadratic. The vacuum polynomial
    x^2 + 6x + 4 = 0 does not emerge.

    Without flux (N=0):
      dV/dx = 2*T*x - 2*C_cas/x^3 = 0
      => x^4 = C_cas/T
    A single equation with a unique solution. No polynomial.

    The brane potential CANNOT produce the quadratic because it gives
    at most a single algebraic condition for phi_min, not a polynomial
    with roots involving sqrt(5).

    Returns
    -------
    dict with honest assessment.
    """
    C_cas = _C_casimir(a_0)
    C_fl = _C_flux(N_flux, a_0)

    # The target polynomial
    # x^2 + 6x + 4 = 0 => x = (-6 +/- sqrt(36-16))/2 = -3 +/- sqrt(5)
    x_plus = -3.0 + math.sqrt(5.0)   # = -0.7639...
    x_minus = -3.0 - math.sqrt(5.0)  # = -5.2360...
    phi_from_x_plus = 0.5 * math.log(abs(x_plus)) if x_plus > 0 else None

    # Note: x = exp(2*phi) > 0 always, but x_plus < 0!
    # So the vacuum polynomial root x = -3+sqrt(5) < 0 cannot correspond
    # to x = exp(2*phi) which is strictly positive.
    # This is already known: the identification uses a DIFFERENT substitution.

    # Without flux: the minimization gives x^4 = C_cas/T (single value)
    no_flux_equation = "x^4 = C_cas / T  =>  phi_min = (1/8)*ln(C_cas/T)"

    # With flux: quintic 3*C_flux*x^5 + 2*T*x^4 - 2*C_cas = 0
    with_flux_equation = "3*C_flux*x^5 + 2*T*x^4 - 2*C_cas = 0  (quintic)"

    return {
        "polynomial_emerges": False,
        "target_polynomial": "x^2 + 6x + 4 = 0",
        "target_roots": {"x_plus": x_plus, "x_minus": x_minus},
        "actual_equation_no_flux": no_flux_equation,
        "actual_equation_with_flux": with_flux_equation,
        "gap1_status": "OPEN",
        "reason": (
            "The brane tension generates a dilaton potential whose minimization "
            "gives a single algebraic condition for phi_min (x^4 = C_cas/T without "
            "flux, or a quintic with flux). Neither reduces to the quadratic "
            "x^2 + 6x + 4 = 0. The vacuum polynomial requires a different mechanism "
            "or a deeper principle connecting T, C_cas, and the algebraic structure "
            "of the golden ratio."
        ),
        "additional_obstacle": (
            "The vacuum polynomial root x = -3+sqrt(5) = -0.764 is NEGATIVE, "
            "while x = exp(2*phi) > 0 always. The original identification of the "
            "polynomial uses a different variable substitution, not x = exp(2*phi). "
            "The brane potential naturally uses exponentials of phi, making the "
            "connection to the quadratic polynomial even more indirect."
        ),
        "C_cas": C_cas,
        "C_flux": C_fl,
    }


# ===================================================================
# 9. summarize_braneworld
# ===================================================================
def summarize_braneworld():
    """Full summary of braneworld physics on S^2.

    Returns
    -------
    dict with all key results and a formatted report string.
    """
    # 1. Brane potential at phi_vev
    bp = brane_tension_potential(T=1.0, phi=PHI_VEV)

    # 2. Casimir at phi_vev
    cp = casimir_potential(PHI_VEV)

    # 3. Combined potential (example)
    comb = combined_potential(PHI_VEV, T=1.0, N_flux=1)

    # 4. Stabilization scan
    stab = find_brane_stabilization(a_0=1.0, N_flux=1)

    # 5. Deficit angle for typical M_6
    # M_6 ~ 3 TeV = 3e12 eV, T ~ (1 TeV)^4 = (1e12)^4
    T_example = (1e12) ** 4  # eV^4
    M_6_example = 3e12        # eV
    da = deficit_angle(T_example, M_6_example)

    # 6. Brane coupling
    bc = brane_localized_coupling()

    # 7. Matter spectrum
    ms = brane_matter_spectrum()

    # 8. Gap 1
    g1 = does_brane_close_gap1()

    report_lines = [
        "=" * 72,
        "BRANEWORLD PHYSICS ON S^2 -- SUMMARY",
        "=" * 72,
        "",
        "Setup: 6D bulk (EH + GB gravity, dilaton, gauge fields from S^2 isometry)",
        "       4D brane at a point on S^2 carrying all SM matter",
        "       Codimension-2 brane => conical deficit angle",
        "",
        "--- 1. BRANE TENSION POTENTIAL ---",
        "V_brane(phi) = T * exp(4*phi)",
        "At phi_vev = {:.4f}: V_brane = T * {:.6f}".format(
            PHI_VEV, math.exp(4.0 * PHI_VEV)),
        "dV/dphi = 4*T*exp(4*phi): always positive for T > 0",
        "=> Brane tension drives phi toward -infinity (shrinking S^2)",
        "",
        "--- 2. CASIMIR POTENTIAL ---",
        "V_Casimir(phi) = C_cas * exp(-4*phi)",
        "C_cas = zeta_{S^2}(-1/2) / (4*pi*a_0^4) = -1/(16*pi*a_0^4)",
        "zeta_{S^2}(-1/2) = -1/4  =>  C_cas < 0",
        "=> Casimir drives phi toward +infinity (expanding S^2)",
        "",
        "CRITICAL: For T > 0 and C_cas < 0, BOTH terms drive phi -> -inf!",
        "  V_brane ~ T*exp(4*phi) -> 0 as phi -> -inf",
        "  V_Casimir ~ C_cas*exp(-4*phi) -> -inf as phi -> -inf",
        "No minimum exists for brane + Casimir alone with T > 0.",
        "",
        "--- 3. COMBINED POTENTIAL (with flux) ---",
        "V = T*exp(4*phi) + C_cas*exp(-4*phi) + C_flux*exp(6*phi)",
        "Flux adds exp(6*phi) term which can compete with Casimir exp(-4*phi)",
        "",
        "--- 4. STABILIZATION SCAN ---",
        "Scanned {} T values with N_flux = {}".format(
            stab["n_T_scanned"], stab["N_flux"]),
        "Minima found: {}".format(stab["n_minima_found"]),
        "Verdict: {}".format(stab["verdict"]),
        "",
        "--- 5. DEFICIT ANGLE ---",
        "Example: T = (1 TeV)^4, M_6 = 3 TeV",
        "delta = {:.4f} rad = {:.2f} deg".format(
            da["delta_rad"], da["delta_deg"]),
        "Physical (delta < 2*pi): {}".format(da["physical_single"]),
        "Fraction of sphere remaining (1 brane): {:.6f}".format(
            da["fraction_of_sphere_1brane"]),
        "",
        "--- 6. BRANE-LOCALIZED COUPLINGS ---",
        "alpha_brane = exp(4*phi_vev)/(4*pi) = {:.10f}".format(
            bc["alpha_brane"]),
        "alpha_EM = {:.10f}".format(ALPHA_EM),
        "Match: {:.2f} ppm".format(bc["alpha_match_ppm"]),
        "",
        "--- 7. BRANE MATTER SPECTRUM ---",
        "R_phys = {:.4e} m".format(ms["R_phys_m"]),
        "M_6 = {:.4e} eV".format(ms["M_6_eV"]),
        "Lightest graviton KK: l=2, m = {:.4e} eV".format(
            ms["graviton_modes"][0]["m_eV"]),
        "Lightest gauge KK: l=1, m = {:.4e} eV".format(
            ms["gauge_modes"][0]["m_eV"]),
        "",
        "--- 8. GAP 1 (VACUUM POLYNOMIAL) ---",
        "Status: {}".format(g1["gap1_status"]),
        "Target: {}".format(g1["target_polynomial"]),
        "Actual equation (no flux): {}".format(g1["actual_equation_no_flux"]),
        "Actual equation (with flux): {}".format(g1["actual_equation_with_flux"]),
        "",
        g1["reason"],
        "",
        g1["additional_obstacle"],
        "",
        "--- 9. KEY CONCLUSIONS ---",
        "1. Brane tension generates V_brane = T*exp(4*phi) -- a real potential",
        "2. For T > 0 and C_cas < 0, no minimum from brane+Casimir alone",
        "3. Adding flux can create a minimum; T tunes phi_min",
        "4. Gap 1 (vacuum polynomial) NOT closed by the brane mechanism",
        "5. The deficit angle constrains T < M_6^4 for physical compactification",
        "6. alpha matching is exact at tree level (by construction of phi_vev)",
        "7. The braneworld interpretation is physically consistent and standard",
        "8. Paper 3 should adopt this: gravity+gauge bulk, matter on brane",
        "=" * 72,
    ]

    report = "\n".join(report_lines)

    return {
        "brane_potential": bp,
        "casimir_potential": cp,
        "combined_potential": comb,
        "stabilization_scan": stab,
        "deficit_angle": da,
        "brane_coupling": bc,
        "matter_spectrum": ms,
        "gap1_analysis": g1,
        "report": report,
    }


# ===================================================================
# __main__
# ===================================================================
if __name__ == "__main__":
    result = summarize_braneworld()
    print(result["report"])
