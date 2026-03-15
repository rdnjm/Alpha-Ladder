"""
Flux stabilization of the volume modulus via quantized 2-form flux.

The casimir_stabilization module showed that the pure graviton Casimir
energy on S^2 has a NEGATIVE coefficient (A < 0), which cannot produce
a stable minimum for the volume modulus sigma.  This module adds a
quantized 2-form flux F_2 threading the S^2, which provides the missing
positive-definite contribution to the effective potential.

The effective 4D potential with flux is:

    V(sigma) = A * e^{4*sigma} + B * e^{2*sigma} + C * e^{6*sigma}

where:
    A = Casimir coefficient (negative, from casimir_stabilization module)
    B = -chi / (2 * a_0^2) = curvature coefficient (B = -1.0 for S^2)
    C = N^2 / (32 * pi^2 * a_0^6) = flux coefficient (POSITIVE for N != 0)

The flux quantum N is an integer (Dirac quantization condition for a
2-form on S^2).  C > 0 always, providing the repulsive term at small
radii that was missing from the pure Casimir potential.

Key results:

    1. Flux stabilization WORKS -- a stable minimum exists for all N >= 1.

    2. The minimum is at Planck-scale radius (a_0 ~ l_Pl) with
       Planck-scale dilaton mass (m_phi ~ M_Pl).

    3. A Planck-mass dilaton is completely invisible at all experimental
       scales: exp(-M_Pl * r) ~ 0 for ANY r > 10^{-35} m.

    4. This trivially dissolves the 3854x screening gap (Gap #1): the
       Yukawa suppression exp(-m_phi * r) is so extreme that the dilaton
       fifth force vanishes identically at all measurable distances.

    5. But it also eliminates all dilaton phenomenology -- the framework
       reduces to pure GR at observable scales, which is consistent but
       uninteresting.

    6. To get sub-Planck masses (and observable dilaton effects) requires
       a_0 >> l_Pl (large extra dimensions), which introduces a new
       hierarchy problem.

    7. The flux quantum N is the only new discrete parameter (not
       continuous), so this is a genuine stabilization mechanism with
       no fine-tuning.

All calculations use pure Python math (no numpy/scipy).
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50


# ---------------------------------------------------------------------------
# 1. Flux coefficient
# ---------------------------------------------------------------------------

def compute_flux_coefficient(N, a_0=1.0):
    """
    Compute the flux contribution coefficient for N quanta of 2-form flux
    threading an S^2 of radius a_0.

    The quantized 2-form flux on S^2 contributes an energy density that
    scales as e^{6*sigma} in the effective 4D potential.  The coefficient
    is determined by the Dirac quantization condition:

        C = N^2 / (32 * pi^2 * a_0^6)

    This arises from:
        - The flux energy (1/2) |F_2|^2 integrated over S^2
        - The quantization condition: integral of F_2 over S^2 = 2*pi*N
        - The volume element: 4*pi*a_0^2 for S^2 of radius a_0
        - Combining: C = N^2 / (2 * (4*pi*a_0^2)^2 * a_0^2)
                       = N^2 / (32 * pi^2 * a_0^6)

    Parameters
    ----------
    N : int
        Flux quantum number (integer by Dirac quantization).
    a_0 : float
        Radius of the S^2 in Planck units (default 1.0).

    Returns
    -------
    dict with keys:
        C : float -- flux coefficient
        N : int -- flux quantum number
        a_0 : float -- S^2 radius
        formula : str -- description of the formula
        positive : bool -- always True for N != 0
    """
    C = N ** 2 / (32.0 * math.pi ** 2 * a_0 ** 6)

    return {
        "C": C,
        "N": N,
        "a_0": a_0,
        "formula": (
            f"C = N^2 / (32 * pi^2 * a_0^6) = {N}^2 / (32 * pi^2 * "
            f"{a_0}^6) = {C:.10e}"
        ),
        "positive": N != 0,
    }


# ---------------------------------------------------------------------------
# 2. Flux potential on a sigma grid
# ---------------------------------------------------------------------------

def compute_flux_potential(sigma_grid=None, N=1, a_0=1.0, constants=None):
    """
    Compute the full effective potential V(sigma) with Casimir, curvature,
    and flux contributions on a grid of sigma values.

    The potential is:

        V(sigma) = A * e^{4*sigma} + B * e^{2*sigma} + C * e^{6*sigma}

    where:
        A = Casimir coefficient from casimir_stabilization module (negative)
        B = -chi / (2 * a_0^2) with chi = 2 for S^2 (negative)
        C = N^2 / (32 * pi^2 * a_0^6) from flux quantization (positive)

    Parameters
    ----------
    sigma_grid : list of float or None
        Grid of sigma values.  If None, uses 81 points from -3 to 5.
    N : int
        Flux quantum number (default 1).
    a_0 : float
        Radius of S^2 in Planck units (default 1.0).
    constants : ignored
        Planck units used throughout.

    Returns
    -------
    dict with keys:
        sigma_grid : list of float
        V_casimir : list of float -- A * e^{4*sigma} at each point
        V_curvature : list of float -- B * e^{2*sigma} at each point
        V_flux : list of float -- C * e^{6*sigma} at each point
        V_total : list of float -- sum at each point
        A : float -- Casimir coefficient
        B : float -- curvature coefficient
        C : float -- flux coefficient
        N : int -- flux quantum
        a_0 : float -- S^2 radius
        chi : int -- Euler characteristic (2 for S^2)
    """
    from alpha_ladder_core.casimir_stabilization import compute_casimir_energy_s2

    # Build sigma grid if not provided
    if sigma_grid is None:
        n_points = 81
        sigma_min = -3.0
        sigma_max = 5.0
        step = (sigma_max - sigma_min) / (n_points - 1)
        sigma_grid = [sigma_min + i * step for i in range(n_points)]

    # Casimir coefficient A
    casimir = compute_casimir_energy_s2(a_radius=a_0)
    A = casimir["coefficient"] / a_0 ** 4

    # Curvature coefficient B (S^2: chi = 2)
    chi = 2
    B = -chi / (2.0 * a_0 ** 2)

    # Flux coefficient C
    flux = compute_flux_coefficient(N, a_0)
    C = flux["C"]

    # Evaluate on grid
    V_casimir = []
    V_curvature = []
    V_flux = []
    V_total = []

    for sig in sigma_grid:
        e2s = math.exp(2.0 * sig)
        e4s = e2s * e2s
        e6s = e4s * e2s

        v_cas = A * e4s
        v_curv = B * e2s
        v_flx = C * e6s

        V_casimir.append(v_cas)
        V_curvature.append(v_curv)
        V_flux.append(v_flx)
        V_total.append(v_cas + v_curv + v_flx)

    return {
        "sigma_grid": sigma_grid,
        "V_casimir": V_casimir,
        "V_curvature": V_curvature,
        "V_flux": V_flux,
        "V_total": V_total,
        "A": A,
        "B": B,
        "C": C,
        "N": N,
        "a_0": a_0,
        "chi": chi,
    }


# ---------------------------------------------------------------------------
# 3. Find flux-stabilized minimum
# ---------------------------------------------------------------------------

def find_flux_minimum(N=1, a_0=1.0, constants=None):
    """
    Find the stable minimum of the flux-corrected potential.

    Setting V'(sigma) = 0 with the substitution u = e^{2*sigma} gives
    the quadratic equation:

        6C * u^2 + 4A * u + 2B = 0

    with solutions:

        u = [-4A +/- sqrt(16A^2 - 48CB)] / (12C)

    For S^2 with B < 0 and C > 0:
        discriminant = 16A^2 - 48CB = 16A^2 + 48C|B| > 0 always

    The positive root u_+ always exists, giving a real sigma_0.

    The second derivative at the stationary point is:

        V''(sigma_0) = 16A * u^2 + 4B * u + 36C * u^3

    where u = e^{2*sigma_0}.  A positive V'' confirms a local minimum.

    Parameters
    ----------
    N : int
        Flux quantum number (default 1).
    a_0 : float
        Radius of S^2 in Planck units (default 1.0).
    constants : ignored
        Planck units used throughout.

    Returns
    -------
    dict with keys:
        minimum_exists : bool
        sigma_0 : float or None
        a_stabilized : float or None -- stabilized radius a_0 * exp(-sigma_0)
        V_at_minimum : float or None
        V_double_prime : float or None
        is_minimum : bool
        u_plus : float or None
        u_minus : float or None
        discriminant : float
        A : float -- Casimir coefficient
        B : float -- curvature coefficient
        C : float -- flux coefficient
        N : int
        a_0 : float
    """
    from alpha_ladder_core.casimir_stabilization import compute_casimir_energy_s2

    # Coefficients
    casimir = compute_casimir_energy_s2(a_radius=a_0)
    A = casimir["coefficient"] / a_0 ** 4

    chi = 2
    B = -chi / (2.0 * a_0 ** 2)

    flux = compute_flux_coefficient(N, a_0)
    C = flux["C"]

    # Solve 6C u^2 + 4A u + 2B = 0
    # Quadratic: a_q = 6C, b_q = 4A, c_q = 2B
    a_q = 6.0 * C
    b_q = 4.0 * A
    c_q = 2.0 * B

    discriminant = b_q ** 2 - 4.0 * a_q * c_q  # 16A^2 - 48CB

    minimum_exists = False
    sigma_0 = None
    a_stabilized = None
    V_at_minimum = None
    V_double_prime = None
    is_minimum = False
    u_plus = None
    u_minus = None

    if discriminant >= 0 and a_q != 0:
        sqrt_disc = math.sqrt(discriminant)
        u_plus = (-b_q + sqrt_disc) / (2.0 * a_q)
        u_minus = (-b_q - sqrt_disc) / (2.0 * a_q)

        # We need u > 0 (since u = e^{2*sigma})
        # Try u_plus first (the larger root)
        u_stable = None
        if u_plus is not None and u_plus > 0:
            u_stable = u_plus
        elif u_minus is not None and u_minus > 0:
            u_stable = u_minus

        if u_stable is not None:
            sigma_0 = 0.5 * math.log(u_stable)
            a_stabilized = a_0 * math.exp(-sigma_0)

            # V at minimum
            e2s = u_stable
            e4s = e2s * e2s
            e6s = e4s * e2s
            V_at_minimum = A * e4s + B * e2s + C * e6s

            # V''(sigma_0) = 16A * u^2 + 4B * u + 36C * u^3
            V_double_prime = (
                16.0 * A * e4s
                + 4.0 * B * e2s
                + 36.0 * C * e6s
            )

            is_minimum = V_double_prime > 0
            minimum_exists = is_minimum

    return {
        "minimum_exists": minimum_exists,
        "sigma_0": sigma_0,
        "a_stabilized": a_stabilized,
        "V_at_minimum": V_at_minimum,
        "V_double_prime": V_double_prime,
        "is_minimum": is_minimum,
        "u_plus": u_plus,
        "u_minus": u_minus,
        "discriminant": discriminant,
        "A": A,
        "B": B,
        "C": C,
        "N": N,
        "a_0": a_0,
    }


# ---------------------------------------------------------------------------
# 4. Dilaton mass from flux stabilization
# ---------------------------------------------------------------------------

def compute_flux_dilaton_mass(N=1, a_0=1.0, constants=None):
    """
    Compute the dilaton mass from flux stabilization of the volume modulus.

    If a stable minimum exists at sigma_0, the dilaton mass-squared in
    Planck units is:

        m_phi^2 = V''(sigma_0)

    The physical mass in eV is:

        m_phi_eV = M_Pl * sqrt(m_phi^2)

    where M_Pl = 1.22089e28 eV is the Planck mass.

    The Compton wavelength is:

        lambda = hbar * c / m_phi

    The experimental threshold for fifth-force searches is 2 meV.

    Parameters
    ----------
    N : int
        Flux quantum number (default 1).
    a_0 : float
        Radius of S^2 in Planck units (default 1.0).
    constants : ignored
        Planck units used throughout.

    Returns
    -------
    dict with keys:
        m_phi_squared : float or None -- V''(sigma_0) in Planck units
        m_phi_eV : float or None -- mass in eV
        lambda_compton_m : float or None -- Compton wavelength in meters
        exceeds_threshold : bool or None -- whether m > 2 meV
        N : int
        a_0 : float
        first_principles : bool -- True if minimum exists
        comparison_with_threshold : str
        mass_scale : str -- "planck_scale" or "sub_planck"
    """
    M_Pl_eV = 1.22089e28       # Planck mass in eV
    hbar_c_eV_m = 1.9733e-7    # hbar * c in eV * m
    threshold_eV = 2e-3         # 2 meV

    minimum = find_flux_minimum(N=N, a_0=a_0, constants=constants)

    m_phi_squared = None
    m_phi_eV = None
    lambda_compton_m = None
    exceeds_threshold = None
    first_principles = False
    mass_scale = None

    if minimum["minimum_exists"]:
        m_phi_squared = minimum["V_double_prime"]
        first_principles = True

        if m_phi_squared > 0:
            m_phi_planck = math.sqrt(m_phi_squared)
            m_phi_eV = M_Pl_eV * m_phi_planck
            lambda_compton_m = hbar_c_eV_m / m_phi_eV
            exceeds_threshold = m_phi_eV > threshold_eV

            # Classify the mass scale
            if m_phi_planck > 0.01:
                mass_scale = "planck_scale"
            else:
                mass_scale = "sub_planck"

            comparison = (
                f"m_phi = {m_phi_eV:.4e} eV.  "
                f"Threshold: {threshold_eV:.0e} eV.  "
                f"{'Exceeds' if exceeds_threshold else 'Below'} threshold.  "
                f"Compton wavelength: {lambda_compton_m:.4e} m."
            )
        else:
            comparison = (
                "m_phi^2 <= 0: tachyonic or massless.  "
                "Flux stabilization produces an unstable direction."
            )
    else:
        comparison = (
            "No stable minimum exists.  The dilaton mass cannot be "
            "computed from flux stabilization with N = {N}, a_0 = {a_0}."
        )

    return {
        "m_phi_squared": m_phi_squared,
        "m_phi_eV": m_phi_eV,
        "lambda_compton_m": lambda_compton_m,
        "exceeds_threshold": exceeds_threshold,
        "N": N,
        "a_0": a_0,
        "first_principles": first_principles,
        "comparison_with_threshold": comparison,
        "mass_scale": mass_scale,
    }


# ---------------------------------------------------------------------------
# 5. Scan over flux quanta
# ---------------------------------------------------------------------------

def scan_flux_quanta(N_max=10, a_0=1.0, constants=None):
    """
    Scan over flux quantum numbers N = 1 to N_max and compute the
    stabilized minimum and dilaton mass for each.

    This reveals the discrete landscape of flux vacua: each integer N
    gives a distinct vacuum with its own sigma_0, stabilized radius,
    and dilaton mass.

    Parameters
    ----------
    N_max : int
        Maximum flux quantum to scan (default 10).
    a_0 : float
        Radius of S^2 in Planck units (default 1.0).
    constants : ignored
        Planck units used throughout.

    Returns
    -------
    dict with keys:
        results : list of dict -- one entry per N with minimum and mass data
        N_max : int
        a_0 : float
        all_stable : bool -- True if every N has a stable minimum
        mass_range_eV : tuple of (min, max) mass in eV, or None
        description : str
    """
    results = []
    all_stable = True
    masses_eV = []

    for N in range(1, N_max + 1):
        minimum = find_flux_minimum(N=N, a_0=a_0, constants=constants)
        mass = compute_flux_dilaton_mass(N=N, a_0=a_0, constants=constants)

        entry = {
            "N": N,
            "minimum_exists": minimum["minimum_exists"],
            "sigma_0": minimum["sigma_0"],
            "a_stabilized": minimum["a_stabilized"],
            "V_at_minimum": minimum["V_at_minimum"],
            "V_double_prime": minimum["V_double_prime"],
            "is_minimum": minimum["is_minimum"],
            "C": minimum["C"],
            "m_phi_squared": mass["m_phi_squared"],
            "m_phi_eV": mass["m_phi_eV"],
            "lambda_compton_m": mass["lambda_compton_m"],
            "mass_scale": mass["mass_scale"],
        }
        results.append(entry)

        if not minimum["minimum_exists"]:
            all_stable = False

        if mass["m_phi_eV"] is not None:
            masses_eV.append(mass["m_phi_eV"])

    mass_range_eV = None
    if masses_eV:
        mass_range_eV = (min(masses_eV), max(masses_eV))

    description = (
        f"Scanned N = 1..{N_max} with a_0 = {a_0}.  "
        f"{'All' if all_stable else 'Not all'} flux quanta produce stable minima.  "
    )
    if mass_range_eV is not None:
        description += (
            f"Mass range: {mass_range_eV[0]:.4e} to {mass_range_eV[1]:.4e} eV."
        )

    return {
        "results": results,
        "N_max": N_max,
        "a_0": a_0,
        "all_stable": all_stable,
        "mass_range_eV": mass_range_eV,
        "description": description,
    }


# ---------------------------------------------------------------------------
# 6. Gap closure analysis
# ---------------------------------------------------------------------------

def compute_flux_gap_closure(constants=None):
    """
    Determine whether flux stabilization closes theoretical gaps #1 and #3.

    Gap #3 (dilaton mass from first principles):
        Closed if a stable minimum exists for any N >= 1.  The flux
        potential provides the positive e^{6*sigma} term that the pure
        Casimir potential lacked, creating a genuine minimum.

    Gap #1 (screening sufficiency -- the 3854x ratio):
        Closed if the dilaton mass is large enough that screening is
        irrelevant.  For a Planck-mass dilaton:
            exp(-M_Pl * r) ~ 0 for any r > 10^{-35} m
        The fifth force vanishes identically at all experimental scales.
        The 3854x ratio between predicted and required screening becomes
        irrelevant because the dilaton is effectively decoupled.

    Parameters
    ----------
    constants : ignored
        Planck units used throughout.

    Returns
    -------
    dict with keys:
        gap1_resolved : bool
        gap3_resolved : bool
        gap1_mechanism : str
        gap3_mechanism : str
        screening_suppression_at_1m : float
        screening_suppression_at_1AU : float
        dilaton_effectively_decoupled : bool
        honest_assessment : str
    """
    # Use N=1 as the minimal flux quantum
    minimum = find_flux_minimum(N=1, a_0=1.0, constants=constants)
    mass = compute_flux_dilaton_mass(N=1, a_0=1.0, constants=constants)

    gap3_resolved = minimum["minimum_exists"]

    # Screening suppression: exp(-m_phi * r) in natural units
    # m_phi in eV, r in meters: need m_phi / (hbar*c) * r
    # 1 / (hbar*c) = 1 / (1.9733e-7 eV*m) = 5.068e6 eV^{-1} m^{-1}
    hbar_c_eV_m = 1.9733e-7
    AU_m = 1.496e11

    screening_suppression_1m = 0.0
    screening_suppression_1AU = 0.0
    dilaton_decoupled = False
    gap1_resolved = False

    if mass["m_phi_eV"] is not None and mass["m_phi_eV"] > 0:
        m_over_hbarc = mass["m_phi_eV"] / hbar_c_eV_m  # in m^{-1}

        # exp(-m * r) for r = 1 m
        exponent_1m = m_over_hbarc * 1.0
        # For Planck-scale masses this exponent is ~ 10^{34}, so
        # exp(-exponent) is effectively zero.  Clamp to avoid overflow.
        if exponent_1m > 700:
            screening_suppression_1m = 0.0
        else:
            screening_suppression_1m = math.exp(-exponent_1m)

        # exp(-m * r) for r = 1 AU
        exponent_1AU = m_over_hbarc * AU_m
        if exponent_1AU > 700:
            screening_suppression_1AU = 0.0
        else:
            screening_suppression_1AU = math.exp(-exponent_1AU)

        # Dilaton is effectively decoupled if suppression < 10^{-100}
        # at 1 meter (which is true for any Planck-scale mass)
        dilaton_decoupled = exponent_1m > 230  # e^{-230} < 10^{-100}

        # Gap #1 is resolved if dilaton is decoupled
        gap1_resolved = dilaton_decoupled

    gap3_mechanism = (
        "The quantized 2-form flux on S^2 provides a positive e^{6*sigma} "
        "term in the effective potential.  Combined with the negative Casimir "
        "(e^{4*sigma}) and curvature (e^{2*sigma}) terms, this creates a "
        "genuine minimum.  The minimum exists for all integer N >= 1 with "
        "no fine-tuning required."
    )

    gap1_mechanism = (
        "The flux-stabilized dilaton has Planck-scale mass (m_phi ~ M_Pl "
        "for a_0 ~ l_Pl).  The Yukawa suppression factor exp(-m_phi * r) "
        "is effectively zero at ALL experimental scales: "
        f"exp(-m_phi * 1m) ~ {screening_suppression_1m:.1e}.  "
        "The 3854x ratio between predicted and required screening amplitude "
        "becomes irrelevant because the dilaton fifth force is identically "
        "zero.  The framework reduces to pure GR at observable distances."
    )

    honest_assessment = (
        "Flux stabilization provides a honest, first-principles resolution "
        "of both Gap #1 and Gap #3.  However, the resolution is in some "
        "sense trivial: a Planck-mass dilaton is completely invisible, which "
        "means the alpha ladder framework makes no testable dilaton predictions "
        "at observable scales.  The predicted G_vacuum from the phi^2/2 bridge "
        "remains valid, but the Yukawa screening profile becomes unobservable.  "
        "To get phenomenologically interesting (sub-Planck) dilaton masses "
        "requires a_0 >> l_Pl, i.e. large extra dimensions, which introduces "
        "a new hierarchy problem.  The flux quantum N is discrete and natural "
        "(no fine-tuning), but the overall scale is set by a_0."
    )

    return {
        "gap1_resolved": gap1_resolved,
        "gap3_resolved": gap3_resolved,
        "gap1_mechanism": gap1_mechanism,
        "gap3_mechanism": gap3_mechanism,
        "screening_suppression_at_1m": screening_suppression_1m,
        "screening_suppression_at_1AU": screening_suppression_1AU,
        "dilaton_effectively_decoupled": dilaton_decoupled,
        "honest_assessment": honest_assessment,
    }


# ---------------------------------------------------------------------------
# 7. Summary / dashboard entry point
# ---------------------------------------------------------------------------

def summarize_flux_stabilization(constants=None):
    """
    Run the full flux stabilization pipeline and return a summary for
    the Streamlit dashboard.

    This is the main entry point.  It computes:
        1. The flux potential with N = 1 (default)
        2. The stabilized minimum
        3. The dilaton mass
        4. A scan over N = 1..5
        5. Gap closure analysis

    Parameters
    ----------
    constants : ignored
        Planck units used throughout.

    Returns
    -------
    dict with keys:
        potential : dict -- flux potential on sigma grid
        minimum : dict -- stabilized minimum data
        mass : dict -- dilaton mass from flux stabilization
        scan : dict -- scan over N = 1..5
        gap_closure : dict -- gap #1 and #3 closure status
        overall_assessment : str
        first_principles : bool
        flux_quantum_N : int
        stabilized_radius : float or None
        dilaton_mass_eV : float or None
    """
    potential = compute_flux_potential(N=1, a_0=1.0, constants=constants)
    minimum = find_flux_minimum(N=1, a_0=1.0, constants=constants)
    mass = compute_flux_dilaton_mass(N=1, a_0=1.0, constants=constants)
    scan = scan_flux_quanta(N_max=5, a_0=1.0, constants=constants)
    gap_closure = compute_flux_gap_closure(constants=constants)

    first_principles = minimum["minimum_exists"]

    if first_principles:
        overall = (
            "Flux stabilization succeeds: the quantized 2-form flux on S^2 "
            "provides a positive e^{6*sigma} contribution that, combined with "
            "the negative Casimir and curvature terms, creates a stable minimum "
            "for the volume modulus.  The minimum exists for all integer flux "
            f"quanta N >= 1.  For N = 1, a_0 = 1 (Planck units): "
            f"sigma_0 = {minimum['sigma_0']:.6f}, "
            f"stabilized radius = {minimum['a_stabilized']:.6e} l_Pl, "
            f"dilaton mass = {mass['m_phi_eV']:.4e} eV "
            f"({mass['mass_scale']}).  "
            "The Planck-scale mass means the dilaton is completely decoupled "
            "from observable physics.  The 3854x screening gap (Gap #1) is "
            "trivially resolved because exp(-M_Pl * r) = 0 at all experimental "
            "scales.  Gap #3 (dilaton mass from first principles) is also closed."
        )
    else:
        overall = (
            "Flux stabilization failed to produce a minimum for N = 1.  "
            "This is unexpected and suggests a computational error."
        )

    return {
        "potential": potential,
        "minimum": minimum,
        "mass": mass,
        "scan": scan,
        "gap_closure": gap_closure,
        "overall_assessment": overall,
        "first_principles": first_principles,
        "flux_quantum_N": 1,
        "stabilized_radius": minimum["a_stabilized"],
        "dilaton_mass_eV": mass["m_phi_eV"],
    }
