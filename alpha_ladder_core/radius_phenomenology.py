"""
Radius phenomenology for the Alpha Ladder framework.

Maps observable consequences as a function of the internal radius a_0
of the compact 2-manifold.  The dilaton mass scales as:

    m_phi ~ M_Pl * (l_Pl / a_0)

This single unfixed parameter creates a sharp phenomenological landscape:
at a_0 ~ l_Pl the dilaton is invisible (pure GR); at a_0 ~ 30 um the
dilaton sits in the sub-mm testable window; at a_0 >> 0.1 mm the theory
is excluded by fifth-force searches.

All calculations use pure Python + math + decimal (no numpy/scipy).
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50

# ---------------------------------------------------------------------------
# Physical constants (Decimal, 50-digit precision)
# ---------------------------------------------------------------------------

_HBAR = Decimal("1.054571817e-34")       # J*s
_C = Decimal("2.99792458e8")             # m/s
_G = Decimal("6.67430e-11")              # m^3 kg^-1 s^-2
_EV_J = Decimal("1.602176634e-19")       # 1 eV in Joules
_HBAR_C_EV_M = Decimal("1.9733e-7")     # hbar*c in eV*m
_L_PL = Decimal("1.61625e-35")          # Planck length in m
_M_PL_EV = Decimal("1.22089e28")        # Planck mass in eV
_ALPHA_FIFTH = Decimal("0.618")         # 2 / (2*omega + 3), BD coupling


# ---------------------------------------------------------------------------
# 1. Dilaton mass vs internal radius
# ---------------------------------------------------------------------------

def compute_dilaton_mass_vs_radius(a_0_values=None):
    """Compute dilaton mass m_phi = M_Pl * l_Pl / a_0 over a grid of radii.

    Parameters
    ----------
    a_0_values : list of float or None
        Internal radii in metres.  If None, a log-spaced grid of 50 points
        from 1e-35 m to 1e-2 m is generated.

    Returns
    -------
    dict with keys:
        a_0_values      : list of float  -- internal radii (m)
        m_phi_eV        : list of float  -- dilaton mass (eV)
        lambda_compton_m: list of float  -- Compton wavelength (m)
        scaling_law     : str
        l_Pl            : float
        M_Pl_eV         : float
    """
    if a_0_values is None:
        log10_min = Decimal("-35")
        log10_max = Decimal("-2")
        n = 50
        a_0_values = []
        for i in range(n):
            exponent = float(log10_min + i * (log10_max - log10_min) / (n - 1))
            a_0_values.append(10.0 ** exponent)

    m_phi_eV = []
    lambda_compton_m = []

    l_pl_f = float(_L_PL)
    m_pl_f = float(_M_PL_EV)
    hbar_c_f = float(_HBAR_C_EV_M)

    for a0 in a_0_values:
        m = m_pl_f * (l_pl_f / a0)
        m_phi_eV.append(m)
        if m > 0.0:
            lam = hbar_c_f / m
        else:
            lam = float("inf")
        lambda_compton_m.append(lam)

    return {
        "a_0_values": list(a_0_values),
        "m_phi_eV": m_phi_eV,
        "lambda_compton_m": lambda_compton_m,
        "scaling_law": "m_phi = M_Pl * l_Pl / a_0",
        "l_Pl": l_pl_f,
        "M_Pl_eV": m_pl_f,
    }


# ---------------------------------------------------------------------------
# 2. Screening amplitude vs internal radius
# ---------------------------------------------------------------------------

def compute_screening_amplitude_vs_radius(r_test=0.1, a_0_values=None):
    """Compute the screening amplitude alpha_eff at a fixed test distance.

    alpha_eff = alpha_fifth * exp(-r_test / lambda_compton)
    where lambda_compton = hbar*c / m_phi(a_0).

    Parameters
    ----------
    r_test : float
        Test distance in metres (default 0.1 m = 10 cm).
    a_0_values : list of float or None
        Internal radii in metres.  If None, uses the same default grid
        as compute_dilaton_mass_vs_radius.

    Returns
    -------
    dict with keys:
        a_0_values       : list of float
        alpha_eff        : list of float
        r_test           : float
        alpha_fifth      : float  (0.618)
        log10_alpha_eff  : list of float  (clamped at -300 for underflow)
    """
    mass_data = compute_dilaton_mass_vs_radius(a_0_values)
    a_0_list = mass_data["a_0_values"]
    lambdas = mass_data["lambda_compton_m"]

    alpha_fifth_f = float(_ALPHA_FIFTH)
    alpha_eff = []
    log10_alpha_eff = []

    for lam in lambdas:
        if lam <= 0.0 or lam == float("inf"):
            alpha_eff.append(0.0)
            log10_alpha_eff.append(-300.0)
            continue

        exponent = -r_test / lam
        # Guard against extreme underflow
        if exponent < -690.0:
            alpha_eff.append(0.0)
            log10_alpha_eff.append(-300.0)
        else:
            val = alpha_fifth_f * math.exp(exponent)
            alpha_eff.append(val)
            if val > 0.0:
                log10_val = math.log10(val)
                log10_alpha_eff.append(max(log10_val, -300.0))
            else:
                log10_alpha_eff.append(-300.0)

    return {
        "a_0_values": a_0_list,
        "alpha_eff": alpha_eff,
        "r_test": r_test,
        "alpha_fifth": alpha_fifth_f,
        "log10_alpha_eff": log10_alpha_eff,
    }


# ---------------------------------------------------------------------------
# 3. Experimental bounds
# ---------------------------------------------------------------------------

def compute_experimental_bounds():
    """Return a catalogue of experimental bounds on (alpha, lambda) from fifth-force searches.

    These are stored as data, not computed.

    Returns
    -------
    dict with keys:
        bounds      : dict of dicts, keyed by experiment identifier
        n_bounds    : int
        description : str
    """
    bounds = {
        "eot_wash_2006": {
            "name": "Eot-Wash 2006",
            "lambda_min_m": 56e-6,
            "alpha_max": 1.0,
            "reference": "Kapner et al., PRL 98, 021101 (2007)",
            "description": "Sub-mm torsion balance at U. Washington",
        },
        "eot_wash_2020": {
            "name": "Eot-Wash 2020 (projected)",
            "lambda_min_m": 30e-6,
            "alpha_max": 0.1,
            "reference": "Lee et al., projected",
            "description": "Next-generation sub-mm test",
        },
        "cassini": {
            "name": "Cassini",
            "lambda_min_m": 1.496e11,
            "alpha_max": 2.3e-5,
            "reference": "Bertotti et al., Nature 425, 374 (2003)",
            "description": "Solar conjunction PPN gamma",
        },
        "microscope": {
            "name": "MICROSCOPE",
            "lambda_min_m": 1e7,
            "alpha_max": 1e-15,
            "reference": "Touboul et al., PRL 129, 121102 (2022)",
            "description": "Equivalence principle test in orbit",
        },
        "lunar_laser_ranging": {
            "name": "Lunar Laser Ranging",
            "lambda_min_m": 3.84e8,
            "alpha_max": 1e-11,
            "reference": "Williams et al., CQG 29, 184004 (2012)",
            "description": "Earth-Moon distance precision tracking",
        },
    }

    return {
        "bounds": bounds,
        "n_bounds": len(bounds),
        "description": (
            "Fifth-force experimental bounds on Yukawa coupling (alpha, lambda).  "
            "Each entry gives the minimum Compton wavelength probed and the maximum "
            "coupling strength allowed at that scale."
        ),
    }


# ---------------------------------------------------------------------------
# 4. Testable window
# ---------------------------------------------------------------------------

def compute_testable_window():
    """Find the range of a_0 where the framework is experimentally testable.

    The lower bound is a_0 > l_Pl (physical requirement).
    The upper bound is a_0 < a_0_max where m_phi(a_0_max) equals
    the meV threshold (~2e-3 eV).

    Returns
    -------
    dict with keys:
        a_0_min            : float  -- Planck length (m)
        a_0_max_testable   : float  -- upper radius for meV threshold (m)
        a_0_eot_wash       : float  -- radius giving m_phi ~ 7 meV (m)
        a_0_cassini        : float  -- radius where Cassini bound applies (m)
        m_phi_at_eot_wash  : float  -- dilaton mass at a_0_eot_wash (eV)
        m_phi_at_cassini_radius : float  -- dilaton mass at a_0_cassini (eV)
        window_exists      : bool
        window_range_m     : tuple of float  (a_0_min, a_0_max_testable)
        window_description : str
        experiments_in_window : list of str
    """
    l_pl = float(_L_PL)
    m_pl = float(_M_PL_EV)

    # meV threshold: m_phi = 2e-3 eV
    m_threshold_eV = 2e-3
    a_0_max = m_pl * l_pl / m_threshold_eV

    # Eot-Wash optimal radius: m_phi ~ 7 meV  => a_0 ~ M_Pl * l_Pl / 7e-3
    m_eot_wash_eV = 7e-3
    a_0_eot_wash = m_pl * l_pl / m_eot_wash_eV

    # Cassini bound mass: at solar scales the screening amplitude
    # must satisfy alpha < 2.3e-5.  The effective Cassini mass threshold
    # is roughly m_cassini ~ 1.35e-17 eV.
    m_cassini_min_eV = 1.35e-17
    a_0_cassini = m_pl * l_pl / m_cassini_min_eV

    m_phi_at_eot_wash = m_pl * l_pl / a_0_eot_wash
    m_phi_at_cassini = m_pl * l_pl / a_0_cassini

    window_exists = a_0_max > l_pl

    return {
        "a_0_min": l_pl,
        "a_0_max_testable": a_0_max,
        "a_0_eot_wash": a_0_eot_wash,
        "a_0_cassini": a_0_cassini,
        "m_phi_at_eot_wash": m_phi_at_eot_wash,
        "m_phi_at_cassini_radius": m_phi_at_cassini,
        "window_exists": window_exists,
        "window_range_m": (l_pl, a_0_max),
        "window_description": (
            "The testable window spans a_0 from the Planck length "
            "({:.3e} m) to {:.3e} m, where the dilaton mass drops below "
            "the meV threshold.  Within this range, sub-mm torsion balance "
            "experiments (Eot-Wash class) can probe dilaton-mediated forces."
        ).format(l_pl, a_0_max),
        "experiments_in_window": [
            "Eot-Wash 2006",
            "Eot-Wash 2020 (projected)",
        ],
    }


# ---------------------------------------------------------------------------
# 5. Full phenomenology at a specific radius
# ---------------------------------------------------------------------------

def compute_phenomenology_at_radius(a_0, constants=None):
    """Compute ALL phenomenological consequences at a specific internal radius.

    Parameters
    ----------
    a_0 : float
        Internal radius of the compact 2-manifold (metres).
    constants : SimpleNamespace or None
        Physical constants from get_constants().  Currently unused but
        accepted for API consistency with other modules.

    Returns
    -------
    dict with keys:
        a_0               : float  -- input radius (m)
        m_phi_eV          : float  -- dilaton mass (eV)
        lambda_compton_m  : float  -- Compton wavelength (m)
        alpha_at_10cm     : float  -- screening amplitude at 0.1 m
        alpha_at_1m       : float  -- screening amplitude at 1.0 m
        alpha_at_1AU      : float  -- screening amplitude at 1 AU
        ppn_gamma_minus_1 : float  -- |gamma - 1| at 1 AU
        passes_cassini    : bool
        passes_eot_wash   : bool
        classification    : str    -- one of "invisible", "sub-mm",
                                      "marginal", "excluded"
    """
    l_pl = float(_L_PL)
    m_pl = float(_M_PL_EV)
    hbar_c = float(_HBAR_C_EV_M)
    alpha_fifth = float(_ALPHA_FIFTH)

    # Dilaton mass and Compton wavelength
    m_phi = m_pl * (l_pl / a_0)
    if m_phi > 0.0:
        lam_c = hbar_c / m_phi
    else:
        lam_c = float("inf")

    # Screening amplitudes at representative distances
    r_10cm = 0.1
    r_1m = 1.0
    r_1AU = 1.496e11  # metres

    def _alpha_eff(r):
        if lam_c <= 0.0 or lam_c == float("inf"):
            return 0.0
        exponent = -r / lam_c
        if exponent < -690.0:
            return 0.0
        return alpha_fifth * math.exp(exponent)

    alpha_10cm = _alpha_eff(r_10cm)
    alpha_1m = _alpha_eff(r_1m)
    alpha_1AU = _alpha_eff(r_1AU)

    # PPN |gamma - 1| at solar-system scale
    # In BD theory: |gamma - 1| = 2 * alpha_eff at the relevant scale
    ppn_gamma_minus_1 = 2.0 * alpha_1AU

    # Cassini bound: |gamma - 1| < 2.3e-5
    passes_cassini = ppn_gamma_minus_1 < 2.3e-5

    # Eot-Wash bound: at sub-mm scales, alpha < 1.0 for lambda > 56 um
    # We check if the Compton wavelength is probed and the coupling exceeds bound
    bounds = compute_experimental_bounds()["bounds"]
    eot_wash = bounds["eot_wash_2006"]
    if lam_c >= eot_wash["lambda_min_m"]:
        passes_eot_wash = alpha_10cm <= eot_wash["alpha_max"]
    else:
        # Compton wavelength below experimental reach -- not probed, passes by default
        passes_eot_wash = True

    # Classification
    if m_phi >= 1e9:
        classification = "invisible"
    elif m_phi >= 2e-3:
        classification = "sub-mm"
    elif m_phi >= 1e-3:
        classification = "marginal"
    else:
        # Below meV: check if excluded
        if not passes_cassini or not passes_eot_wash:
            classification = "excluded"
        else:
            classification = "excluded"  # sub-meV generically excluded by fifth-force

    return {
        "a_0": a_0,
        "m_phi_eV": m_phi,
        "lambda_compton_m": lam_c,
        "alpha_at_10cm": alpha_10cm,
        "alpha_at_1m": alpha_1m,
        "alpha_at_1AU": alpha_1AU,
        "ppn_gamma_minus_1": ppn_gamma_minus_1,
        "passes_cassini": passes_cassini,
        "passes_eot_wash": passes_eot_wash,
        "classification": classification,
    }


# ---------------------------------------------------------------------------
# 6. Dashboard entry point
# ---------------------------------------------------------------------------

def summarize_radius_phenomenology(constants=None):
    """Compute the full phenomenological landscape as a function of a_0.

    This is the dashboard entry point.  It assembles:
    - mass vs radius curve
    - testable window
    - screening at lab and solar scales
    - key radii (Planck, Eot-Wash threshold, Cassini threshold, meV threshold)
    - an honest overall assessment

    Parameters
    ----------
    constants : SimpleNamespace or None
        Physical constants from get_constants().

    Returns
    -------
    dict with keys:
        mass_curve           : dict from compute_dilaton_mass_vs_radius
        testable_window      : dict from compute_testable_window
        key_radii            : dict of notable radii (m) and their physics
        experimental_bounds  : dict from compute_experimental_bounds
        overall_assessment   : str
        framework_prediction : str
    """
    mass_curve = compute_dilaton_mass_vs_radius()
    testable_window = compute_testable_window()
    experimental_bounds = compute_experimental_bounds()

    l_pl = float(_L_PL)
    m_pl = float(_M_PL_EV)

    key_radii = {
        "planck": {
            "a_0_m": l_pl,
            "m_phi_eV": m_pl,
            "description": "Dilaton at Planck mass, invisible, pure GR recovered",
        },
        "eot_wash_optimal": {
            "a_0_m": testable_window["a_0_eot_wash"],
            "m_phi_eV": testable_window["m_phi_at_eot_wash"],
            "description": "Dilaton at ~7 meV, optimal for Eot-Wash class experiments",
        },
        "mev_threshold": {
            "a_0_m": testable_window["a_0_max_testable"],
            "m_phi_eV": 2e-3,
            "description": "Dilaton at 2 meV, below this the 3854x gap returns",
        },
        "cassini_threshold": {
            "a_0_m": testable_window["a_0_cassini"],
            "m_phi_eV": testable_window["m_phi_at_cassini_radius"],
            "description": "Dilaton mass where Cassini PPN bound applies",
        },
    }

    overall_assessment = (
        "The framework does not predict a_0 from first principles.  "
        "With flux stabilization, a_0 is determined by the minimum of a "
        "potential that involves Planck-scale ingredients, giving a_0 ~ l_Pl.  "
        "To get a_0 in the testable window (30 um to 0.1 mm) requires a "
        "hierarchy mechanism not present in the minimal 6D EH+GB framework.  "
        "However, the framework makes a conditional prediction: IF a_0 is "
        "in the sub-mm range, THEN specific signatures should appear in "
        "next-generation gravity experiments.  The framework tells "
        "experimentalists WHERE to look, parametrized by a_0."
    )

    framework_prediction = (
        "Conditional prediction: for a_0 in [{:.1e}, {:.1e}] m, "
        "a Yukawa-type fifth force with coupling alpha = {:.3f} and "
        "Compton wavelength lambda = hbar*c / (M_Pl * l_Pl / a_0) should "
        "be detectable by sub-mm torsion balance experiments.  The strongest "
        "signal appears at a_0 ~ 30 um, giving m_phi ~ 7 meV and "
        "lambda ~ 28 um.  No signal is predicted if a_0 ~ l_Pl (pure GR)."
    ).format(
        l_pl,
        testable_window["a_0_max_testable"],
        float(_ALPHA_FIFTH),
    )

    return {
        "mass_curve": mass_curve,
        "testable_window": testable_window,
        "key_radii": key_radii,
        "experimental_bounds": experimental_bounds,
        "overall_assessment": overall_assessment,
        "framework_prediction": framework_prediction,
    }
