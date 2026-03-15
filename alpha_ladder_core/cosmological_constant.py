"""
Cosmological constant analysis for the Alpha Ladder framework.

The flux-stabilized vacuum energy V_min is O(1) in Planck units.
The observed cosmological constant Lambda ~ 2.888e-122 M_Pl^4.
The discrepancy is ~122 orders of magnitude -- the worst fine-tuning
problem in physics.

This module:
1. Extracts V_min from flux stabilization and converts to physical units.
2. Compares with the observed Lambda.
3. Catalogs known CC mechanisms and why none resolve it here.
4. Scans V_min vs N to show the problem is universal.
5. Applies Weinberg's 1989 no-go theorem.

Honest conclusion: the CC problem is NOT solved by this framework.
This is not a failure specific to Alpha Ladder -- it is the universal
cosmological constant problem that affects ALL known theories of gravity.

All calculations use pure Python math (no numpy/scipy).
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

_LAMBDA_OBS_PLANCK = 2.888e-122    # Observed Lambda in Planck units (M_Pl^4)
_M_PL_EV = 1.22089e28             # Planck mass in eV
_M_PL_KG = 2.176434e-8            # Planck mass in kg
_L_PL = 1.61625e-35               # Planck length in m
_RHO_LAMBDA_SI = 5.96e-27         # Lambda energy density in kg/m^3


# ---------------------------------------------------------------------------
# 1. Extract vacuum energy
# ---------------------------------------------------------------------------

def extract_vacuum_energy(N=1, a_0=1.0, constants=None):
    """
    Extract the vacuum energy V_min from flux stabilization and convert
    to physical units.

    Calls flux_stabilization.find_flux_minimum() to get V_at_minimum
    in Planck units, then converts to:
    - eV^4 (particle physics units)
    - kg/m^3 (cosmological energy density)
    - Effective Lambda in s^{-2}

    Parameters
    ----------
    N : int
        Flux quantum number (default 1).
    a_0 : float
        Internal radius in Planck units (default 1.0).
    constants : ignored
        Planck units used throughout.

    Returns
    -------
    dict with keys:
        V_min_planck : float or None -- V_min in Planck units (M_Pl^4)
        V_min_eV4 : float or None -- V_min in eV^4
        rho_vacuum_kg_m3 : float or None -- energy density in kg/m^3
        N : int
        a_0 : float
        minimum_exists : bool
        sign : str -- "positive", "negative", or "zero"
        magnitude_planck : float or None -- |V_min| in Planck units
        log10_magnitude : float or None -- log10(|V_min|) in Planck units
    """
    from alpha_ladder_core.flux_stabilization import find_flux_minimum

    result = find_flux_minimum(N=N, a_0=a_0, constants=constants)
    V_min = result.get("V_at_minimum")
    minimum_exists = result["minimum_exists"]

    V_min_eV4 = None
    rho_vacuum = None
    sign = "zero"
    magnitude = None
    log10_mag = None

    if V_min is not None and minimum_exists:
        # V_min is in Planck units (M_Pl^4 = 1)
        # Convert to eV^4: multiply by M_Pl_eV^4
        V_min_eV4 = V_min * (_M_PL_EV ** 4)

        # Convert to kg/m^3 using rho = V * M_Pl^4 / (hbar^3 * c^5)
        # In Planck units, rho (kg/m^3) = V_min * M_Pl^4 * c / (hbar^3 * G^2)
        # Simpler: rho_Pl = M_Pl * c^2 / L_Pl^3 = 5.155e96 kg/m^3
        rho_planck = 5.155e96  # Planck density in kg/m^3
        rho_vacuum = V_min * rho_planck

        if V_min > 0:
            sign = "positive"
        elif V_min < 0:
            sign = "negative"

        magnitude = abs(V_min)
        if magnitude > 0:
            log10_mag = math.log10(magnitude)

    return {
        "V_min_planck": V_min,
        "V_min_eV4": V_min_eV4,
        "rho_vacuum_kg_m3": rho_vacuum,
        "N": N,
        "a_0": a_0,
        "minimum_exists": minimum_exists,
        "sign": sign,
        "magnitude_planck": magnitude,
        "log10_magnitude": log10_mag,
    }


# ---------------------------------------------------------------------------
# 2. Compare with observation
# ---------------------------------------------------------------------------

def compare_with_observation(N=1, a_0=1.0, constants=None):
    """
    Compare the flux-stabilized vacuum energy with the observed
    cosmological constant.

    The observed Lambda:
        Lambda_obs ~ 2.888e-122 M_Pl^4
        rho_Lambda ~ 5.96e-27 kg/m^3

    Parameters
    ----------
    N : int
        Flux quantum number.
    a_0 : float
        Internal radius in Planck units.
    constants : ignored

    Returns
    -------
    dict with keys:
        V_min_planck : float or None
        Lambda_obs_planck : float
        ratio : float or None -- |V_min| / Lambda_obs
        log10_ratio : float or None
        discrepancy_orders : int or None -- number of orders of magnitude
        is_fine_tuning : bool -- True if ratio >> 1
        Lambda_obs_eV4 : float
        Lambda_obs_kg_m3 : float
        honest_assessment : str
    """
    vacuum = extract_vacuum_energy(N=N, a_0=a_0, constants=constants)

    Lambda_obs = _LAMBDA_OBS_PLANCK
    Lambda_obs_eV4 = Lambda_obs * (_M_PL_EV ** 4)
    Lambda_obs_kg_m3 = _RHO_LAMBDA_SI

    ratio = None
    log10_ratio = None
    discrepancy_orders = None
    is_fine_tuning = False

    V_min = vacuum["V_min_planck"]

    if V_min is not None and V_min != 0:
        ratio = abs(V_min) / Lambda_obs
        if ratio > 0:
            log10_ratio = math.log10(ratio)
            discrepancy_orders = int(round(log10_ratio))
            is_fine_tuning = log10_ratio > 1.0  # More than 10x discrepancy

    honest_assessment = (
        "The flux-stabilized vacuum energy is O(1) in Planck units "
        f"(V_min ~ {V_min:.4e} M_Pl^4 for N={N}), while the observed "
        f"cosmological constant is Lambda ~ {Lambda_obs:.3e} M_Pl^4. "
    )

    if discrepancy_orders is not None:
        honest_assessment += (
            f"The discrepancy is ~{discrepancy_orders} orders of magnitude. "
            f"This is the cosmological constant problem -- the worst "
            f"fine-tuning problem in physics. It is NOT specific to the "
            f"Alpha Ladder framework; it affects ALL known theories that "
            f"include gravity and quantum fields."
        )
    else:
        honest_assessment += (
            "Unable to compute the ratio (V_min may be zero or undefined)."
        )

    return {
        "V_min_planck": V_min,
        "Lambda_obs_planck": Lambda_obs,
        "ratio": ratio,
        "log10_ratio": log10_ratio,
        "discrepancy_orders": discrepancy_orders,
        "is_fine_tuning": is_fine_tuning,
        "Lambda_obs_eV4": Lambda_obs_eV4,
        "Lambda_obs_kg_m3": Lambda_obs_kg_m3,
        "honest_assessment": honest_assessment,
    }


# ---------------------------------------------------------------------------
# 3. CC mechanisms catalog
# ---------------------------------------------------------------------------

def analyze_cc_mechanisms():
    """
    Catalog known approaches to the cosmological constant problem
    and assess their applicability to the Alpha Ladder framework.

    Returns
    -------
    dict with keys:
        mechanisms : list of dict
        n_mechanisms : int
        any_resolution : bool -- False (none resolve it here)
        universal_problem : bool -- True
        honest_assessment : str
    """
    mechanisms = [
        {
            "name": "Supersymmetry",
            "description": (
                "SUSY relates boson and fermion vacuum energies, providing "
                "cancellation. Exact SUSY gives Lambda = 0."
            ),
            "applicable_here": False,
            "reason": (
                "The Alpha Ladder framework is non-supersymmetric. Even if "
                "SUSY were added, broken SUSY gives Lambda ~ M_SUSY^4 ~ "
                "(1 TeV)^4, still 60 orders too large."
            ),
            "resolves_problem": False,
        },
        {
            "name": "Anthropic / Landscape",
            "description": (
                "In a vast landscape of vacua (e.g., string landscape with "
                "~10^{500} vacua), anthropic selection picks vacua with "
                "small Lambda compatible with structure formation."
            ),
            "applicable_here": False,
            "reason": (
                "The Alpha Ladder framework has a discrete landscape from "
                "flux quanta N, but this gives O(10) vacua, not 10^{500}. "
                "Insufficient for the anthropic argument."
            ),
            "resolves_problem": False,
        },
        {
            "name": "Sequestering",
            "description": (
                "Sequestering mechanisms isolate the vacuum energy from "
                "the gravitational sector, preventing large V_min from "
                "contributing to the cosmological constant."
            ),
            "applicable_here": False,
            "reason": (
                "Requires specific global symmetries in the action that "
                "are not present in the minimal 6D EH+GB framework."
            ),
            "resolves_problem": False,
        },
        {
            "name": "Self-tuning",
            "description": (
                "Scalar fields (like the dilaton) dynamically adjust to "
                "cancel the vacuum energy. The dilaton sigma could in "
                "principle play this role."
            ),
            "applicable_here": False,
            "reason": (
                "Weinberg's 1989 no-go theorem shows that no scalar "
                "adjustment mechanism can dynamically relax Lambda without "
                "fine-tuning, given a smooth, Lorentz-invariant potential. "
                "The dilaton sigma has exactly such a potential."
            ),
            "resolves_problem": False,
        },
        {
            "name": "Unimodular gravity",
            "description": (
                "In unimodular gravity, the cosmological constant is an "
                "integration constant rather than a vacuum energy contribution."
            ),
            "applicable_here": False,
            "reason": (
                "Does not explain why the integration constant is small. "
                "Shifts the problem from 'why is V_min small?' to 'why is "
                "the integration constant small?'"
            ),
            "resolves_problem": False,
        },
    ]

    return {
        "mechanisms": mechanisms,
        "n_mechanisms": len(mechanisms),
        "any_resolution": False,
        "universal_problem": True,
        "honest_assessment": (
            "None of the known approaches to the cosmological constant "
            "problem resolve it within the Alpha Ladder framework. This is "
            "not surprising: the CC problem is unsolved in ALL known theories "
            "of gravity. The framework is honest about this limitation -- "
            "it does not manufacture a false resolution."
        ),
    }


# ---------------------------------------------------------------------------
# 4. CC scan over flux quanta
# ---------------------------------------------------------------------------

def compute_cc_scan(N_max=10, a_0=1.0, constants=None):
    """
    Compute V_min for each flux quantum N = 1..N_max and compare
    with the observed Lambda.

    Parameters
    ----------
    N_max : int
        Maximum flux quantum (default 10).
    a_0 : float
        Internal radius in Planck units (default 1.0).
    constants : ignored

    Returns
    -------
    dict with keys:
        scan_results : list of dict -- one per N with V_min, log10_ratio
        N_max : int
        a_0 : float
        all_O1_planck : bool -- True if all |V_min| are O(1)
        any_close_to_observation : bool -- False (expected)
        min_log10_ratio : float -- smallest discrepancy
        max_log10_ratio : float -- largest discrepancy
        honest_assessment : str
    """
    from alpha_ladder_core.flux_stabilization import find_flux_minimum

    scan_results = []
    log10_ratios = []

    for N in range(1, N_max + 1):
        result = find_flux_minimum(N=N, a_0=a_0, constants=constants)
        V_min = result.get("V_at_minimum")

        entry = {
            "N": N,
            "V_min_planck": V_min,
            "minimum_exists": result["minimum_exists"],
        }

        if V_min is not None and V_min != 0:
            ratio = abs(V_min) / _LAMBDA_OBS_PLANCK
            log10_r = math.log10(ratio) if ratio > 0 else 0.0
            entry["log10_ratio"] = log10_r
            entry["ratio"] = ratio
            log10_ratios.append(log10_r)

            # Classify magnitude
            abs_v = abs(V_min)
            if abs_v > 0.01:
                entry["magnitude_class"] = "O(1) Planck"
            elif abs_v > 1e-10:
                entry["magnitude_class"] = "sub-Planck"
            else:
                entry["magnitude_class"] = "small"
        else:
            entry["log10_ratio"] = None
            entry["ratio"] = None
            entry["magnitude_class"] = "N/A"

        scan_results.append(entry)

    all_O1 = all(
        e.get("magnitude_class") == "O(1) Planck"
        for e in scan_results
        if e.get("magnitude_class") != "N/A"
    )

    any_close = any(
        e.get("log10_ratio") is not None and e["log10_ratio"] < 10
        for e in scan_results
    )

    min_ratio = min(log10_ratios) if log10_ratios else 0.0
    max_ratio = max(log10_ratios) if log10_ratios else 0.0

    return {
        "scan_results": scan_results,
        "N_max": N_max,
        "a_0": a_0,
        "all_O1_planck": all_O1,
        "any_close_to_observation": any_close,
        "min_log10_ratio": min_ratio,
        "max_log10_ratio": max_ratio,
        "honest_assessment": (
            f"Scanned N = 1..{N_max} at a_0 = {a_0}. ALL vacuum energies "
            f"are O(1) in Planck units. The smallest discrepancy with "
            f"observation is ~10^{{{min_ratio:.0f}}} and the largest is "
            f"~10^{{{max_ratio:.0f}}}. No flux quantum produces a vacuum "
            f"energy anywhere near the observed Lambda ~ 10^{{-122}}. "
            f"This confirms that the CC problem is not an artifact of a "
            f"specific N value but a universal feature of the framework."
        ),
    }


# ---------------------------------------------------------------------------
# 5. Weinberg no-go theorem
# ---------------------------------------------------------------------------

def assess_no_go_theorem():
    """
    Apply Weinberg's 1989 no-go theorem to the Alpha Ladder framework.

    Weinberg proved that no scalar field adjustment mechanism can
    dynamically relax the cosmological constant to zero (or a small
    value) if:
    1. The potential is a smooth function of the scalar fields.
    2. The vacuum is Lorentz-invariant.
    3. The scalar fields have standard kinetic terms.

    All three conditions hold for the dilaton sigma in the Alpha Ladder
    framework.

    Returns
    -------
    dict with keys:
        applies_to_framework : bool -- True
        conditions : list of dict -- each with name, met, reason
        all_conditions_met : bool
        smooth_potential : bool
        lorentz_invariant : bool
        standard_kinetic : bool
        conclusion : str
        reference : str
    """
    conditions = [
        {
            "name": "Smooth potential",
            "met": True,
            "reason": (
                "V(sigma) = A*e^{4s} + B*e^{2s} + C*e^{6s} is an analytic "
                "function of sigma. No cusps, discontinuities, or singular points."
            ),
        },
        {
            "name": "Lorentz-invariant vacuum",
            "met": True,
            "reason": (
                "The framework assumes a 4D Lorentz-invariant vacuum with "
                "sigma = sigma_0 = const. The KK reduction preserves 4D "
                "Poincare invariance."
            ),
        },
        {
            "name": "Standard kinetic terms",
            "met": True,
            "reason": (
                "The dilaton sigma has a canonical kinetic term (1/2)(d sigma)^2 "
                "in the 4D Einstein frame (after Weyl rescaling from the "
                "Jordan frame)."
            ),
        },
    ]

    all_met = all(c["met"] for c in conditions)

    return {
        "applies_to_framework": all_met,
        "conditions": conditions,
        "all_conditions_met": all_met,
        "smooth_potential": conditions[0]["met"],
        "lorentz_invariant": conditions[1]["met"],
        "standard_kinetic": conditions[2]["met"],
        "conclusion": (
            "Weinberg's no-go theorem applies: the dilaton sigma cannot "
            "dynamically relax the cosmological constant to a small value. "
            "The vacuum energy V_min is determined by the shape of V(sigma), "
            "which is O(1) in Planck units for the flux-stabilized potential. "
            "No smooth deformation of the potential can change this without "
            "fine-tuning the coefficients A, B, C."
        ),
        "reference": "Weinberg, Rev. Mod. Phys. 61 (1989) 1",
    }


# ---------------------------------------------------------------------------
# 6. Summary / dashboard entry point
# ---------------------------------------------------------------------------

def summarize_cosmological_constant(constants=None):
    """
    Run the full cosmological constant analysis and return a summary
    for the Streamlit dashboard.

    Parameters
    ----------
    constants : ignored
        Planck units used throughout.

    Returns
    -------
    dict with keys:
        vacuum_energy : dict from extract_vacuum_energy
        comparison : dict from compare_with_observation
        mechanisms : dict from analyze_cc_mechanisms
        cc_scan : dict from compute_cc_scan
        no_go : dict from assess_no_go_theorem
        Lambda_obs_planck : float
        discrepancy_orders : int or None
        is_resolved : bool -- False
        overall_assessment : str
        honest_assessment : str
    """
    vacuum = extract_vacuum_energy(N=1, a_0=1.0, constants=constants)
    comparison = compare_with_observation(N=1, a_0=1.0, constants=constants)
    mechanisms = analyze_cc_mechanisms()
    scan = compute_cc_scan(N_max=10, a_0=1.0, constants=constants)
    no_go = assess_no_go_theorem()

    discrepancy = comparison.get("discrepancy_orders")

    overall_assessment = (
        "The cosmological constant problem is NOT resolved by the Alpha "
        "Ladder framework. The flux-stabilized vacuum energy is O(1) in "
        f"Planck units, while observation requires ~10^{{-122}}. "
    )

    if discrepancy is not None:
        overall_assessment += (
            f"The discrepancy is ~{discrepancy} orders of magnitude. "
        )

    overall_assessment += (
        "This is the universal CC problem that afflicts ALL known theories "
        "of gravity and quantum fields. The framework is honest about this "
        "limitation. Weinberg's no-go theorem (1989) proves that no scalar "
        "adjustment mechanism (including the dilaton) can resolve it without "
        "fine-tuning."
    )

    honest_assessment = (
        "The cosmological constant problem is the hardest unsolved problem "
        "in theoretical physics. No framework -- not string theory, not loop "
        "quantum gravity, not any other approach -- has a satisfactory "
        "resolution. The Alpha Ladder framework joins this universal failure "
        "honestly: V_min is O(1) Planck, and no mechanism within the "
        "framework can reduce it to 10^{-122}. This is not a specific "
        "weakness of the framework but a reflection of a problem that lies "
        "beyond all current theoretical tools."
    )

    return {
        "vacuum_energy": vacuum,
        "comparison": comparison,
        "mechanisms": mechanisms,
        "cc_scan": scan,
        "no_go": no_go,
        "Lambda_obs_planck": _LAMBDA_OBS_PLANCK,
        "discrepancy_orders": discrepancy,
        "is_resolved": False,
        "overall_assessment": overall_assessment,
        "honest_assessment": honest_assessment,
    }
