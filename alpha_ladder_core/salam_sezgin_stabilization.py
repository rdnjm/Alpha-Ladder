"""
Salam-Sezgin / Freund-Rubin stabilization of the S^2 radius.

The existing flux_stabilization module stabilizes the breathing mode sigma
(volume modulus) but NOT the bare radius a_0.  The potential V_min(a_0) is
monotonically decreasing -- no finite a_0 is selected.

This module adds a positive 6D cosmological constant Lambda_6 > 0 to the
action, following the Salam-Sezgin (1984) and Randjbar-Daemi-Salam-Strathdee
(1983) construction.  Lambda_6 contributes a term proportional to the internal
volume that grows with R, providing the repulsive large-R force needed to
create a genuine minimum in the radial potential.

The effective 4D potential for the physical S^2 radius R is:

    V(R) = Lambda_6 * R^2 - 2 / R^2 + n^2 / (2 * R^4)

where:
    Lambda_6 > 0   : 6D cosmological constant (Planck units, 4*pi*M_6^4 = 1)
    n >= 1          : monopole flux quantum (integer, Dirac quantization)
    R               : physical radius of S^2

The three terms:
    1. Lambda_6 * R^2  -- cosmological constant (grows with volume, repulsive)
    2. -2 / R^2        -- internal curvature (Ricci scalar R_2 = 2/R^2)
    3. n^2 / (2*R^4)   -- Dirac-quantized monopole flux energy

The stationarity condition dV/dR = 0 gives:

    Lambda_6 * R^6 + 2 * R^2 - n^2 = 0

This is a cubic in R^2.  For Lambda_6 > 0 and n >= 1, there is exactly one
positive real root R_0^2, giving the stabilized radius.

Key results:
    1. The Salam-Sezgin mechanism DOES select a finite radius R_0 -- unlike
       the pure flux potential in flux_stabilization.py.

    2. For the gauge-matched radius R_phys = 0.550293 l_Pl (from the RDSS
       monopole result), the required Lambda_6 is computed below.

    3. The dilaton mass at the minimum is V''(R_0), giving a physical mass
       scale for the breathing mode fluctuations.

    4. The 4D cosmological constant is V(R_0), which must be compared to
       the observed Lambda_4 ~ 10^{-122} M_Pl^4.

    5. Whether Lambda_6 is "natural" (O(1) in Planck units) or fine-tuned
       is a key question this module answers.

All calculations use pure Python math (no numpy/scipy).
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50


# ---------------------------------------------------------------------------
# Physical constants used in this module
# ---------------------------------------------------------------------------

# Gauge-matched radius from KK gauge matching condition:
# phi_vev = (1/4)*ln(4*pi*alpha_em), R_phys = exp(phi_vev) with a_0 = 1
# For alpha_em = 1/137.036: R_phys ~ 0.550293 Planck lengths
_R_GAUGE_MATCHED = 0.550293  # Planck lengths

# Planck mass in eV
_M_PL_EV = 1.22089e28

# hbar * c in eV * m
_HBAR_C_EV_M = 1.9733e-7

# Observed 4D cosmological constant in Planck units (Lambda_4 / M_Pl^4)
_LAMBDA_4_OBS_PLANCK = 2.888e-122


# ---------------------------------------------------------------------------
# 1. Salam-Sezgin potential V(R)
# ---------------------------------------------------------------------------

def compute_salam_sezgin_potential(R_grid, Lambda_6, n=1):
    """
    Compute the Salam-Sezgin effective potential V(R) on a grid of radii.

    The potential for the physical S^2 radius R is:

        V(R) = Lambda_6 * R^2 - 2 / R^2 + n^2 / (2 * R^4)

    where we work in units with 4*pi*M_6^4 = 1.

    The three terms represent:
        - Lambda_6 * R^2 : 6D cosmological constant (volume energy)
        - -2 / R^2       : internal curvature (S^2 Ricci scalar = 2/R^2)
        - n^2 / (2*R^4)  : quantized monopole flux energy

    Parameters
    ----------
    R_grid : list of float
        Grid of physical radius values (Planck units).  Must be > 0.
    Lambda_6 : float
        6D cosmological constant in Planck units.  Must be > 0.
    n : int
        Monopole flux quantum (default 1).

    Returns
    -------
    dict with keys:
        R_grid : list of float -- input radii
        V_cc : list of float -- Lambda_6 * R^2 at each point
        V_curv : list of float -- -2 / R^2 at each point
        V_flux : list of float -- n^2 / (2*R^4) at each point
        V_total : list of float -- sum at each point
        Lambda_6 : float
        n : int
        formula : str -- description of the potential
    """
    V_cc = []
    V_curv = []
    V_flux = []
    V_total = []

    n_sq_half = n * n / 2.0

    for R in R_grid:
        if R <= 0:
            V_cc.append(float('inf'))
            V_curv.append(float('-inf'))
            V_flux.append(float('inf'))
            V_total.append(float('nan'))
            continue

        R2 = R * R
        R4 = R2 * R2

        v_cc = Lambda_6 * R2
        v_curv = -2.0 / R2
        v_flux = n_sq_half / R4

        V_cc.append(v_cc)
        V_curv.append(v_curv)
        V_flux.append(v_flux)
        V_total.append(v_cc + v_curv + v_flux)

    return {
        "R_grid": list(R_grid),
        "V_cc": V_cc,
        "V_curv": V_curv,
        "V_flux": V_flux,
        "V_total": V_total,
        "Lambda_6": Lambda_6,
        "n": n,
        "formula": (
            f"V(R) = Lambda_6 * R^2 - 2/R^2 + n^2/(2*R^4) "
            f"with Lambda_6 = {Lambda_6:.6e}, n = {n}"
        ),
    }


# ---------------------------------------------------------------------------
# 2. Find the Salam-Sezgin minimum
# ---------------------------------------------------------------------------

def find_ss_minimum(Lambda_6, n=1):
    """
    Find the stabilized radius R_0 by solving the stationarity condition.

    Setting dV/dR = 0:

        2*Lambda_6*R + 4/R^3 - 2*n^2/R^5 = 0

    Multiplying through by R^5 / 2:

        Lambda_6 * R^6 + 2 * R^2 - n^2 = 0       ... (*)

    This is a cubic in w = R^2:

        Lambda_6 * w^3 + 2 * w - n^2 = 0

    For Lambda_6 > 0 and n >= 1, the LHS is strictly increasing in w > 0
    (derivative = 3*Lambda_6*w^2 + 2 > 0), and it goes from -n^2 < 0 at
    w = 0 to +infinity as w -> infinity.  Therefore there is EXACTLY ONE
    positive real root.

    We find it by bisection (robust, no cubic formula needed).

    The second derivative of V at R_0 determines whether this is a minimum:

        V''(R) = 2*Lambda_6 - 12/R^4 + 10*n^2/R^6

    At the minimum, V''(R_0) > 0 is verified.

    Parameters
    ----------
    Lambda_6 : float
        6D cosmological constant (must be > 0).
    n : int
        Monopole flux quantum (default 1).

    Returns
    -------
    dict with keys:
        minimum_exists : bool
        R_0 : float or None -- stabilized radius in Planck units
        R_0_squared : float or None -- w = R_0^2
        V_at_minimum : float or None -- V(R_0)
        V_double_prime : float or None -- V''(R_0)
        is_minimum : bool
        Lambda_6 : float
        n : int
        cubic_equation : str -- the equation solved
        status : str -- DERIVED / FAILED
    """
    if Lambda_6 <= 0:
        return {
            "minimum_exists": False,
            "R_0": None,
            "R_0_squared": None,
            "V_at_minimum": None,
            "V_double_prime": None,
            "is_minimum": False,
            "Lambda_6": Lambda_6,
            "n": n,
            "cubic_equation": "Lambda_6 * w^3 + 2*w - n^2 = 0",
            "status": "FAILED: Lambda_6 must be positive",
        }

    n_sq = float(n * n)

    # Solve Lambda_6 * w^3 + 2*w - n^2 = 0 by bisection on w > 0
    # f(w) = Lambda_6 * w^3 + 2*w - n^2
    # f(0) = -n^2 < 0
    # f(w_max) > 0 for large enough w_max

    def f(w):
        return Lambda_6 * w * w * w + 2.0 * w - n_sq

    # Find upper bound: start at w = n^2 and double until f > 0
    w_lo = 0.0
    w_hi = max(n_sq, 1.0)
    safety = 0
    while f(w_hi) < 0 and safety < 200:
        w_hi *= 2.0
        safety += 1

    if f(w_hi) < 0:
        return {
            "minimum_exists": False,
            "R_0": None,
            "R_0_squared": None,
            "V_at_minimum": None,
            "V_double_prime": None,
            "is_minimum": False,
            "Lambda_6": Lambda_6,
            "n": n,
            "cubic_equation": "Lambda_6 * w^3 + 2*w - n^2 = 0",
            "status": "FAILED: bisection could not bracket root",
        }

    # Bisection: 100 iterations gives ~30 digits of precision
    for _ in range(100):
        w_mid = (w_lo + w_hi) / 2.0
        if f(w_mid) < 0:
            w_lo = w_mid
        else:
            w_hi = w_mid

    w_root = (w_lo + w_hi) / 2.0
    R_0 = math.sqrt(w_root)

    # Evaluate V at minimum
    R2 = w_root
    R4 = R2 * R2
    V_min = Lambda_6 * R2 - 2.0 / R2 + n_sq / (2.0 * R4)

    # Second derivative: V''(R) = 2*Lambda_6 - 12/R^4 + 10*n^2/R^6
    # Derivation:
    #   V'(R) = 2*L*R + 4*R^{-3} - 2*n^2*R^{-5}
    #   V''(R) = 2*L - 12*R^{-4} + 10*n^2*R^{-6}
    R6 = R4 * R2
    V_pp = 2.0 * Lambda_6 - 12.0 / R4 + 10.0 * n_sq / R6

    is_min = V_pp > 0

    return {
        "minimum_exists": is_min,
        "R_0": R_0,
        "R_0_squared": w_root,
        "V_at_minimum": V_min,
        "V_double_prime": V_pp,
        "is_minimum": is_min,
        "Lambda_6": Lambda_6,
        "n": n,
        "cubic_equation": (
            f"Lambda_6 * w^3 + 2*w - n^2 = 0 => "
            f"{Lambda_6:.6e} * w^3 + 2*w - {n_sq:.1f} = 0, "
            f"root w = {w_root:.10e}, R_0 = {R_0:.10e}"
        ),
        "status": "DERIVED" if is_min else "FAILED: stationary point is not a minimum",
    }


# ---------------------------------------------------------------------------
# 3. Dilaton mass from the Salam-Sezgin minimum
# ---------------------------------------------------------------------------

def compute_ss_dilaton_mass(Lambda_6, n=1):
    """
    Compute the dilaton (breathing mode) mass from V''(R_0).

    The mass-squared of the radial fluctuation about R_0 is:

        m_R^2 = V''(R_0)

    in Planck units.  The physical mass in eV:

        m_R_eV = M_Pl * sqrt(V''(R_0))

    and the Compton wavelength:

        lambda_c = hbar*c / m_R

    Parameters
    ----------
    Lambda_6 : float
        6D cosmological constant (Planck units).
    n : int
        Monopole flux quantum (default 1).

    Returns
    -------
    dict with keys:
        m_phi_squared : float or None -- V''(R_0) in Planck units
        m_phi_planck : float or None -- sqrt(V''(R_0))
        m_phi_eV : float or None -- mass in eV
        lambda_compton_m : float or None -- Compton wavelength in meters
        R_0 : float or None -- stabilized radius
        Lambda_6 : float
        n : int
        mass_scale : str -- classification of mass scale
        status : str
    """
    minimum = find_ss_minimum(Lambda_6, n=n)

    if not minimum["minimum_exists"]:
        return {
            "m_phi_squared": None,
            "m_phi_planck": None,
            "m_phi_eV": None,
            "lambda_compton_m": None,
            "R_0": None,
            "Lambda_6": Lambda_6,
            "n": n,
            "mass_scale": None,
            "status": "FAILED: no stable minimum",
        }

    V_pp = minimum["V_double_prime"]
    R_0 = minimum["R_0"]

    if V_pp <= 0:
        return {
            "m_phi_squared": V_pp,
            "m_phi_planck": None,
            "m_phi_eV": None,
            "lambda_compton_m": None,
            "R_0": R_0,
            "Lambda_6": Lambda_6,
            "n": n,
            "mass_scale": "tachyonic_or_massless",
            "status": "FAILED: V''(R_0) <= 0",
        }

    m_planck = math.sqrt(V_pp)
    m_eV = _M_PL_EV * m_planck
    lambda_c = _HBAR_C_EV_M / m_eV if m_eV > 0 else None

    # Classify mass scale
    if m_planck > 0.01:
        mass_scale = "planck_scale"
    elif m_planck > 1e-15:
        mass_scale = "intermediate"
    elif m_planck > 1e-30:
        mass_scale = "sub_eV"
    else:
        mass_scale = "ultra_light"

    return {
        "m_phi_squared": V_pp,
        "m_phi_planck": m_planck,
        "m_phi_eV": m_eV,
        "lambda_compton_m": lambda_c,
        "R_0": R_0,
        "Lambda_6": Lambda_6,
        "n": n,
        "mass_scale": mass_scale,
        "status": "DERIVED",
    }


# ---------------------------------------------------------------------------
# 4. 4D cosmological constant from V(R_0)
# ---------------------------------------------------------------------------

def compute_4d_cosmological_constant(Lambda_6, n=1):
    """
    Compute the effective 4D cosmological constant from the potential at
    the Salam-Sezgin minimum.

    The 4D CC is Lambda_4 = V(R_0), the value of the potential at the
    stabilized radius.  This must be compared to the observed value:

        Lambda_4^obs ~ 2.888e-122 M_Pl^4

    The 122-order discrepancy (or lack thereof) is a key diagnostic.

    Parameters
    ----------
    Lambda_6 : float
        6D cosmological constant (Planck units).
    n : int
        Monopole flux quantum (default 1).

    Returns
    -------
    dict with keys:
        Lambda_4 : float or None -- V(R_0) in Planck units
        Lambda_4_obs : float -- observed value
        ratio : float or None -- Lambda_4 / Lambda_4_obs
        log10_ratio : float or None -- log10(|ratio|)
        discrepancy_orders : float or None -- orders of magnitude discrepancy
        sign : str or None -- "positive", "negative", or "zero"
        Lambda_6 : float
        n : int
        honest_assessment : str
        status : str
    """
    minimum = find_ss_minimum(Lambda_6, n=n)

    if not minimum["minimum_exists"]:
        return {
            "Lambda_4": None,
            "Lambda_4_obs": _LAMBDA_4_OBS_PLANCK,
            "ratio": None,
            "log10_ratio": None,
            "discrepancy_orders": None,
            "sign": None,
            "Lambda_6": Lambda_6,
            "n": n,
            "honest_assessment": "No minimum exists; 4D CC cannot be computed.",
            "status": "FAILED",
        }

    Lambda_4 = minimum["V_at_minimum"]

    # Sign
    if Lambda_4 > 1e-300:
        sign = "positive"
    elif Lambda_4 < -1e-300:
        sign = "negative"
    else:
        sign = "zero"

    # Ratio to observed
    ratio = None
    log10_ratio = None
    discrepancy_orders = None

    if abs(Lambda_4) > 0 and _LAMBDA_4_OBS_PLANCK > 0:
        ratio = Lambda_4 / _LAMBDA_4_OBS_PLANCK
        if abs(ratio) > 0:
            log10_ratio = math.log10(abs(ratio))
            discrepancy_orders = abs(log10_ratio)

    # Honest assessment
    if Lambda_4 > 0 and discrepancy_orders is not None:
        if discrepancy_orders < 1:
            honest = (
                "V(R_0) is positive and within an order of magnitude of "
                "Lambda_4^obs.  This would be remarkable and likely indicates "
                "a coincidence or fine-tuning of Lambda_6."
            )
        elif discrepancy_orders < 10:
            honest = (
                f"V(R_0) is positive but {discrepancy_orders:.1f} orders of "
                "magnitude from Lambda_4^obs.  Significant discrepancy but "
                "potentially addressable by additional contributions."
            )
        else:
            honest = (
                f"V(R_0) is positive and {discrepancy_orders:.1f} orders of "
                "magnitude from Lambda_4^obs.  This is the standard cosmological "
                "constant problem: a Planck-scale potential generically gives "
                "O(1) in Planck units, not O(10^{-122}).  No mechanism in this "
                "minimal framework resolves this."
            )
    elif Lambda_4 < 0:
        honest = (
            "V(R_0) is NEGATIVE (anti-de Sitter vacuum).  The observed universe "
            "has positive Lambda_4, so this vacuum is phenomenologically excluded "
            "without additional uplifting."
        )
    else:
        honest = (
            "V(R_0) is approximately zero.  While tantalizing, achieving exact "
            "cancellation requires fine-tuning of Lambda_6 against the curvature "
            "and flux contributions."
        )

    return {
        "Lambda_4": Lambda_4,
        "Lambda_4_obs": _LAMBDA_4_OBS_PLANCK,
        "ratio": ratio,
        "log10_ratio": log10_ratio,
        "discrepancy_orders": discrepancy_orders,
        "sign": sign,
        "Lambda_6": Lambda_6,
        "n": n,
        "honest_assessment": honest,
        "status": "DERIVED",
    }


# ---------------------------------------------------------------------------
# 5. Scan over Lambda_6 values
# ---------------------------------------------------------------------------

def scan_lambda6(lambda6_range, n=1):
    """
    Scan over a range of Lambda_6 values, computing the stabilized radius,
    dilaton mass, and 4D cosmological constant for each.

    This reveals how the stabilization depends on the 6D cosmological
    constant and whether any "special" values of Lambda_6 produce
    phenomenologically interesting results.

    Parameters
    ----------
    lambda6_range : list of float
        List of Lambda_6 values to scan (all must be > 0).
    n : int
        Monopole flux quantum (default 1).

    Returns
    -------
    dict with keys:
        results : list of dict -- one entry per Lambda_6
        n : int
        n_scanned : int -- number of Lambda_6 values tested
        all_stable : bool -- True if every Lambda_6 gives a minimum
        R_0_range : tuple of (min, max) or None
        mass_range_eV : tuple of (min, max) or None
        description : str
    """
    results = []
    all_stable = True
    radii = []
    masses_eV = []

    for L6 in lambda6_range:
        minimum = find_ss_minimum(L6, n=n)
        mass = compute_ss_dilaton_mass(L6, n=n)
        cc = compute_4d_cosmological_constant(L6, n=n)

        entry = {
            "Lambda_6": L6,
            "minimum_exists": minimum["minimum_exists"],
            "R_0": minimum["R_0"],
            "V_at_minimum": minimum["V_at_minimum"],
            "V_double_prime": minimum["V_double_prime"],
            "m_phi_eV": mass["m_phi_eV"],
            "m_phi_planck": mass.get("m_phi_planck"),
            "mass_scale": mass["mass_scale"],
            "Lambda_4": cc["Lambda_4"],
            "discrepancy_orders": cc["discrepancy_orders"],
        }
        results.append(entry)

        if not minimum["minimum_exists"]:
            all_stable = False
        else:
            if minimum["R_0"] is not None:
                radii.append(minimum["R_0"])
            if mass["m_phi_eV"] is not None:
                masses_eV.append(mass["m_phi_eV"])

    R_0_range = (min(radii), max(radii)) if radii else None
    mass_range_eV = (min(masses_eV), max(masses_eV)) if masses_eV else None

    description = (
        f"Scanned {len(lambda6_range)} values of Lambda_6 with n = {n}.  "
        f"{'All' if all_stable else 'Not all'} values produce stable minima."
    )
    if R_0_range is not None:
        description += (
            f"  R_0 range: [{R_0_range[0]:.6e}, {R_0_range[1]:.6e}] l_Pl."
        )
    if mass_range_eV is not None:
        description += (
            f"  Mass range: [{mass_range_eV[0]:.4e}, {mass_range_eV[1]:.4e}] eV."
        )

    return {
        "results": results,
        "n": n,
        "n_scanned": len(lambda6_range),
        "all_stable": all_stable,
        "R_0_range": R_0_range,
        "mass_range_eV": mass_range_eV,
        "description": description,
    }


# ---------------------------------------------------------------------------
# 6. Invert: find Lambda_6 for a target radius
# ---------------------------------------------------------------------------

def find_lambda6_for_target_radius(R_target, n=1):
    """
    Given a desired stabilized radius R_target, find the Lambda_6 that
    produces it.

    From the stationarity condition:

        Lambda_6 * R^6 + 2*R^2 - n^2 = 0

    Solving for Lambda_6:

        Lambda_6 = (n^2 - 2*R^2) / R^6

    For Lambda_6 > 0 we need:

        n^2 > 2*R^2  =>  R < n / sqrt(2)

    If R >= n/sqrt(2), Lambda_6 <= 0 and the Salam-Sezgin mechanism does
    not apply (the curvature term alone can balance the flux without a
    cosmological constant).

    DERIVED: this is an exact algebraic inversion, not numerical.

    Parameters
    ----------
    R_target : float
        Desired stabilized radius in Planck units.
    n : int
        Monopole flux quantum (default 1).

    Returns
    -------
    dict with keys:
        Lambda_6 : float or None -- required 6D cosmological constant
        R_target : float
        n : int
        R_max : float -- maximum R for which Lambda_6 > 0
        is_natural : bool or None -- True if |Lambda_6| ~ O(1)
        naturalness_ratio : float or None -- Lambda_6 / 1.0
        V_at_minimum : float or None -- 4D CC at this minimum
        m_phi_eV : float or None -- dilaton mass
        formula : str
        honest_assessment : str
        status : str
    """
    n_sq = float(n * n)
    R_max = n / math.sqrt(2.0)

    if R_target <= 0:
        return {
            "Lambda_6": None,
            "R_target": R_target,
            "n": n,
            "R_max": R_max,
            "is_natural": None,
            "naturalness_ratio": None,
            "V_at_minimum": None,
            "m_phi_eV": None,
            "formula": "Lambda_6 = (n^2 - 2*R^2) / R^6",
            "honest_assessment": "R_target must be positive.",
            "status": "FAILED",
        }

    R2 = R_target * R_target
    R6 = R2 * R2 * R2

    Lambda_6 = (n_sq - 2.0 * R2) / R6

    if Lambda_6 <= 0:
        return {
            "Lambda_6": Lambda_6,
            "R_target": R_target,
            "n": n,
            "R_max": R_max,
            "is_natural": None,
            "naturalness_ratio": None,
            "V_at_minimum": None,
            "m_phi_eV": None,
            "formula": (
                f"Lambda_6 = (n^2 - 2*R^2) / R^6 = "
                f"({n_sq} - {2*R2:.6e}) / {R6:.6e} = {Lambda_6:.6e}"
            ),
            "honest_assessment": (
                f"R_target = {R_target:.6f} exceeds R_max = {R_max:.6f} "
                f"for n = {n}.  Lambda_6 = {Lambda_6:.6e} <= 0.  "
                "The Salam-Sezgin mechanism does not apply; use a different "
                "stabilization mechanism or increase the flux quantum n."
            ),
            "status": "FAILED: Lambda_6 <= 0",
        }

    # Verify by computing the minimum and mass
    minimum = find_ss_minimum(Lambda_6, n=n)
    mass = compute_ss_dilaton_mass(Lambda_6, n=n)

    V_min = minimum["V_at_minimum"]
    m_phi_eV = mass["m_phi_eV"]

    # Naturalness: Lambda_6 is "natural" if O(1) in Planck units
    naturalness_ratio = Lambda_6 / 1.0  # ratio to M_Pl^6
    is_natural = 0.01 < abs(Lambda_6) < 100.0

    # Verification: R_0 from find_ss_minimum should match R_target
    R_0_check = minimum["R_0"]
    if R_0_check is not None:
        relative_error = abs(R_0_check - R_target) / R_target
        verification = f"Verification: R_0 = {R_0_check:.10f}, error = {relative_error:.2e}"
    else:
        verification = "Verification failed: no minimum found"
        relative_error = None

    honest = (
        f"For R_target = {R_target:.6f} l_Pl with n = {n}: "
        f"Lambda_6 = {Lambda_6:.6e}.  "
    )
    if is_natural:
        honest += (
            "Lambda_6 is O(1) in Planck units -- this is NATURAL, "
            "requiring no fine-tuning.  "
        )
    elif Lambda_6 > 100:
        honest += (
            f"Lambda_6 = {Lambda_6:.2e} >> 1 in Planck units -- this is "
            "unnaturally LARGE, suggesting the target radius is too small "
            "for this mechanism.  "
        )
    else:
        honest += (
            f"Lambda_6 = {Lambda_6:.2e} << 1 in Planck units -- this is "
            "unnaturally SMALL, requiring fine-tuning.  "
        )
    honest += verification

    return {
        "Lambda_6": Lambda_6,
        "R_target": R_target,
        "n": n,
        "R_max": R_max,
        "is_natural": is_natural,
        "naturalness_ratio": naturalness_ratio,
        "V_at_minimum": V_min,
        "m_phi_eV": m_phi_eV,
        "formula": (
            f"Lambda_6 = (n^2 - 2*R^2) / R^6 = "
            f"({n_sq} - {2*R2:.6e}) / {R6:.6e} = {Lambda_6:.6e}"
        ),
        "honest_assessment": honest,
        "status": "DERIVED",
    }


# ---------------------------------------------------------------------------
# 7. Gauge matching constraint
# ---------------------------------------------------------------------------

def gauge_matching_constraint(alpha_em=None):
    """
    From the KK gauge matching condition, the dilaton vev is:

        phi_vev = (1/4) * ln(4 * pi * alpha_em)

    and the physical S^2 radius (with bare a_0 = 1 in Planck units) is:

        R_phys = a_0 * exp(phi_vev) = exp(phi_vev)

    The KK gauge coupling is g_KK = exp(2 * phi_vev), so R_phys = sqrt(g_KK).

    This function computes what Lambda_6 is needed to stabilize the S^2
    at this gauge-matched radius, and what the resulting physical
    predictions are.

    EMPIRICAL: R_phys depends on alpha_em (measured).
    DERIVED: Lambda_6 is computed from the stationarity condition.

    Parameters
    ----------
    alpha_em : float or None
        Fine-structure constant.  If None, uses 1/137.036 (CODATA 2018).

    Returns
    -------
    dict with keys:
        R_phys : float -- gauge-matched radius in Planck units
        alpha_em : float -- fine-structure constant used
        phi_vev : float -- dilaton vev from gauge matching
        g_KK : float -- KK gauge coupling exp(2*phi_vev)
        Lambda_6 : float or None -- required 6D CC
        is_natural : bool or None
        m_phi_eV : float or None -- dilaton mass
        Lambda_4 : float or None -- 4D CC at minimum
        Lambda_4_obs : float -- observed 4D CC
        discrepancy_orders : float or None
        honest_assessment : str
        status : str
    """
    if alpha_em is None:
        alpha_em = 1.0 / 137.036

    # Gauge matching: alpha_KK = exp(4*phi_vev) / (4*pi) = alpha_em
    # => phi_vev = (1/4) * ln(4*pi*alpha_em)
    phi_vev = 0.25 * math.log(4.0 * math.pi * alpha_em)
    g_KK = math.exp(2.0 * phi_vev)
    R_phys = math.exp(phi_vev)  # = sqrt(g_KK), with a_0 = 1

    # Find Lambda_6 for this radius
    result = find_lambda6_for_target_radius(R_phys, n=1)

    Lambda_6 = result["Lambda_6"]
    m_phi_eV = result["m_phi_eV"]
    V_min = result["V_at_minimum"]

    # 4D CC analysis
    cc = compute_4d_cosmological_constant(Lambda_6, n=1) if Lambda_6 is not None and Lambda_6 > 0 else None

    discrepancy = cc["discrepancy_orders"] if cc is not None else None

    # Build honest assessment
    honest = (
        f"Gauge matching: alpha_em = {alpha_em:.6e}, "
        f"phi_vev = {phi_vev:.6f}, g_KK = {g_KK:.6f}, "
        f"R_phys = exp(phi_vev) = {R_phys:.6f} l_Pl.  "
    )
    if Lambda_6 is not None and Lambda_6 > 0:
        honest += (
            f"Required Lambda_6 = {Lambda_6:.6e} (Planck units).  "
        )
        if result["is_natural"]:
            honest += "This is NATURAL (O(1) in Planck units).  "
        elif Lambda_6 > 100:
            honest += "This is unnaturally LARGE.  "
        else:
            honest += "This is unnaturally SMALL (fine-tuned).  "

        if m_phi_eV is not None:
            honest += f"Dilaton mass: {m_phi_eV:.4e} eV.  "
        if discrepancy is not None:
            honest += (
                f"4D CC discrepancy: {discrepancy:.1f} orders of magnitude "
                "from observed Lambda_4.  "
            )
    else:
        honest += (
            "Lambda_6 <= 0: the gauge-matched radius is too large for "
            "the Salam-Sezgin mechanism with n = 1.  Consider larger n."
        )

    return {
        "R_phys": R_phys,
        "alpha_em": alpha_em,
        "phi_vev": phi_vev,
        "g_KK": g_KK,
        "Lambda_6": Lambda_6,
        "is_natural": result.get("is_natural"),
        "m_phi_eV": m_phi_eV,
        "Lambda_4": V_min,
        "Lambda_4_obs": _LAMBDA_4_OBS_PLANCK,
        "discrepancy_orders": discrepancy,
        "honest_assessment": honest,
        "status": result["status"],
    }


# ---------------------------------------------------------------------------
# 8. Full summary
# ---------------------------------------------------------------------------

def summarize_salam_sezgin(constants=None):
    """
    Run the full Salam-Sezgin stabilization analysis and return a summary.

    This is the main entry point.  It computes:
        1. The gauge-matching constraint (what Lambda_6 is needed)
        2. A scan over Lambda_6 from 1e-4 to 1e4
        3. The potential on a radius grid for the gauge-matched Lambda_6
        4. The dilaton mass and 4D CC
        5. Gap closure analysis

    Parameters
    ----------
    constants : SimpleNamespace or None
        If provided, uses constants.alpha for alpha_em.  Otherwise
        uses CODATA 2018 default.

    Returns
    -------
    dict with keys:
        gauge_match : dict -- gauge matching constraint results
        scan : dict -- scan over Lambda_6 values
        potential : dict or None -- V(R) on grid for gauge-matched Lambda_6
        dilaton_mass : dict or None -- mass at gauge-matched minimum
        cc_analysis : dict or None -- 4D CC analysis
        gap1_status : str -- whether screening gap is resolved
        radius_fixed : bool -- whether a finite R_0 is selected
        overall_assessment : str
        honest_assessment : str
        first_principles : bool
    """
    # Extract alpha_em from constants if available
    alpha_em = None
    if constants is not None:
        try:
            alpha_em = float(constants.alpha)
        except (AttributeError, TypeError):
            alpha_em = None

    # 1. Gauge matching
    gauge = gauge_matching_constraint(alpha_em=alpha_em)

    # 2. Scan Lambda_6 from 1e-4 to 1e4 (logarithmic, 50 points)
    log_min = -4.0
    log_max = 4.0
    n_scan = 50
    step = (log_max - log_min) / (n_scan - 1)
    lambda6_values = [10.0 ** (log_min + i * step) for i in range(n_scan)]
    scan = scan_lambda6(lambda6_values, n=1)

    # 3. Potential on R grid for gauge-matched Lambda_6
    potential = None
    dilaton = None
    cc_analysis = None
    radius_fixed = False

    Lambda_6_gm = gauge["Lambda_6"]
    if Lambda_6_gm is not None and Lambda_6_gm > 0:
        # Build R grid around the gauge-matched radius
        R_center = gauge["R_phys"]
        R_min_grid = max(R_center * 0.1, 0.01)
        R_max_grid = R_center * 5.0
        n_grid = 100
        r_step = (R_max_grid - R_min_grid) / (n_grid - 1)
        R_grid = [R_min_grid + i * r_step for i in range(n_grid)]

        potential = compute_salam_sezgin_potential(R_grid, Lambda_6_gm, n=1)
        dilaton = compute_ss_dilaton_mass(Lambda_6_gm, n=1)
        cc_analysis = compute_4d_cosmological_constant(Lambda_6_gm, n=1)
        radius_fixed = True

    # 4. Gap closure analysis
    # The old flux_stabilization stabilized sigma but NOT a_0.
    # The Salam-Sezgin mechanism stabilizes R = a_0 * exp(sigma) directly.
    # Does this close Gap #1 (radius determination)?

    if radius_fixed and dilaton is not None:
        m_eV = dilaton.get("m_phi_eV")
        if m_eV is not None and m_eV > 0:
            # Screening: exp(-m*r) at 1 meter
            m_over_hbarc = m_eV / _HBAR_C_EV_M
            exponent_1m = m_over_hbarc * 1.0
            dilaton_decoupled = exponent_1m > 230  # exp(-230) < 10^{-100}

            if dilaton_decoupled:
                gap1_status = (
                    "RESOLVED (trivially): dilaton mass is Planck-scale, "
                    f"m_phi = {m_eV:.4e} eV.  "
                    f"exp(-m_phi * 1m) < 10^(-{exponent_1m / math.log(10):.0f}).  "
                    "The fifth force vanishes identically at all experimental "
                    "scales.  The 3854x screening gap is irrelevant."
                )
            else:
                gap1_status = (
                    f"OPEN: dilaton mass m_phi = {m_eV:.4e} eV is sub-Planck.  "
                    "Screening analysis from screening.py still applies."
                )
        else:
            gap1_status = "UNKNOWN: dilaton mass could not be computed."
    else:
        gap1_status = (
            "OPEN: Salam-Sezgin mechanism does not apply for the "
            "gauge-matched radius (Lambda_6 <= 0 or no minimum)."
        )

    # 5. Overall assessment
    if radius_fixed:
        Lambda_6_val = Lambda_6_gm
        R_0_val = gauge["R_phys"]
        m_eV_val = dilaton["m_phi_eV"] if dilaton else None
        cc_val = cc_analysis["Lambda_4"] if cc_analysis else None
        cc_disc = cc_analysis["discrepancy_orders"] if cc_analysis else None

        overall = (
            "The Salam-Sezgin mechanism (Lambda_6 > 0 on S^2) successfully "
            f"stabilizes the S^2 radius at R_0 = {R_0_val:.6f} l_Pl "
            f"(gauge-matched) with Lambda_6 = {Lambda_6_val:.6e}.  "
        )
        if m_eV_val is not None:
            overall += f"Dilaton mass: {m_eV_val:.4e} eV.  "
        if cc_disc is not None:
            overall += (
                f"4D CC: {cc_val:.4e} M_Pl^4, {cc_disc:.1f} orders from observed.  "
            )
        overall += (
            "This is the first mechanism in the Alpha Ladder framework that "
            "selects a FINITE radius for the compact S^2, resolving the "
            "radius determination problem identified in radius_determination.py."
        )
    else:
        overall = (
            "The Salam-Sezgin mechanism does NOT stabilize the S^2 at the "
            "gauge-matched radius for n = 1.  The gauge-matched radius "
            f"R = {gauge['R_phys']:.6f} l_Pl exceeds R_max = {1.0/math.sqrt(2.0):.6f} "
            "for n = 1.  Consider n >= 2 or an alternative mechanism."
        )

    honest = (
        "The Salam-Sezgin mechanism adds a single new parameter (Lambda_6) "
        "to the 6D action.  The stabilized radius depends on this parameter, "
        "so the mechanism trades the radius determination problem for the "
        "question: what sets Lambda_6?  In a UV-complete theory, Lambda_6 "
        "would be determined by the fundamental dynamics.  In the minimal "
        "framework, it is a free parameter.  "
    )
    if cc_analysis is not None and cc_analysis["Lambda_4"] is not None:
        if cc_analysis["discrepancy_orders"] is not None and cc_analysis["discrepancy_orders"] > 10:
            honest += (
                "The 4D cosmological constant at the minimum is generically "
                "O(1) in Planck units, reproducing the standard 122-order "
                "discrepancy.  The CC problem is NOT resolved by this mechanism."
            )
        else:
            honest += (
                "The 4D cosmological constant is surprisingly close to the "
                "observed value.  This warrants careful investigation for "
                "potential fine-tuning."
            )
    else:
        honest += (
            "The 4D CC could not be computed (mechanism did not produce a minimum)."
        )

    return {
        "gauge_match": gauge,
        "scan": scan,
        "potential": potential,
        "dilaton_mass": dilaton,
        "cc_analysis": cc_analysis,
        "gap1_status": gap1_status,
        "radius_fixed": radius_fixed,
        "overall_assessment": overall,
        "honest_assessment": honest,
        "first_principles": radius_fixed,
    }
