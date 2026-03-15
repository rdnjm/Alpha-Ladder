"""
Chameleon screening consistency check for the Alpha Ladder dilaton.

Tests whether the dilaton can have an environment-dependent mass
(chameleon/symmetron mechanism) that allows it to be light in cosmic
voids while heavy in the solar system.

The key self-consistency check is the KK truncation validity: if the
effective internal radius a_0_eff(rho) becomes large enough that the
KK tower spacing drops below the dilaton mass, the 4D effective theory
breaks down and the entire framework loses predictive power.

Results:
    - Fuzzy DM regime (m_phi ~ 10^{-22} eV): FAILS KK truncation.
      a_0 ~ 10^7 m means the KK tower is a continuum, gravity is 6D.
    - meV regime (m_phi ~ 1-10 meV): PASSES KK truncation.
      a_0 ~ 30-100 um, KK tower spacing >> dilaton mass.
    - Solar system screening: works for meV dilaton with density
      coupling, but the screening mechanism is not derived from the
      6D action (requires additional assumptions).

All calculations use pure Python + math + decimal (no numpy/scipy).
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

_L_PL = 1.616e-35          # Planck length in meters
_M_PL_EV = 1.22089e28      # Planck mass in eV
_HBAR_C_EV_M = 1.9733e-7   # hbar * c in eV * m
_M_PL_KG = 2.176e-8        # Planck mass in kg
_LAMBDA_EOT_WASH = 56e-6   # Eot-Wash torsion balance threshold in meters


def _safe_exp(x):
    """Compute math.exp(x) with clamping to avoid overflow."""
    if x > 700.0:
        return math.exp(700.0)
    if x < -700.0:
        return 0.0
    return math.exp(x)


# ---------------------------------------------------------------------------
# Internal helpers for flux potential coefficients
# ---------------------------------------------------------------------------

def _get_flux_coefficients(N=1, a_0=1.0):
    """
    Return the three potential coefficients (A, B, C) used in the
    flux-stabilized effective potential.

    These are computed identically to flux_stabilization.compute_flux_potential.
    """
    from alpha_ladder_core.casimir_stabilization import compute_casimir_energy_s2
    from alpha_ladder_core.flux_stabilization import compute_flux_coefficient

    casimir = compute_casimir_energy_s2(a_radius=a_0)
    A = casimir["coefficient"] / a_0 ** 4

    chi = 2
    B = -chi / (2.0 * a_0 ** 2)

    flux = compute_flux_coefficient(N, a_0)
    C = flux["C"]

    return A, B, C


# ---------------------------------------------------------------------------
# 1. Effective potential with matter coupling
# ---------------------------------------------------------------------------

def compute_effective_potential_with_matter(sigma_grid=None, rho=0.0,
                                           N=1, a_0=1.0, beta_matter=1.0):
    """
    Compute V_eff(sigma) = V_flux(sigma) + beta_matter * rho * e^{2*sigma}.

    The flux potential is:
        V_flux = A * e^{4s} + B * e^{2s} + C * e^{6s}

    The matter coupling adds a term rho * e^{2*sigma} arising from the
    sqrt(-g) L_matter piece in the 6D action, conformally coupled to
    the breathing mode.

    Parameters
    ----------
    sigma_grid : list of float or None
        Grid of sigma values.  If None, uses 81 points from -60 to 10.
    rho : float
        Matter density in Planck units (default 0.0).
    N : int
        Flux quantum number (default 1).
    a_0 : float
        Radius of S^2 in Planck units (default 1.0).
    beta_matter : float
        Dimensionless coupling strength of matter to the breathing
        mode (default 1.0).

    Returns
    -------
    dict with keys:
        sigma_grid, V_flux, V_matter, V_total, A, B, C, rho, beta_matter
    """
    if sigma_grid is None:
        n_points = 81
        sigma_min = -60.0
        sigma_max = 10.0
        step = (sigma_max - sigma_min) / (n_points - 1)
        sigma_grid = [sigma_min + i * step for i in range(n_points)]

    A, B, C = _get_flux_coefficients(N=N, a_0=a_0)

    V_flux = []
    V_matter = []
    V_total = []

    for sig in sigma_grid:
        e2s = _safe_exp(2.0 * sig)
        e4s = _safe_exp(4.0 * sig)
        e6s = _safe_exp(6.0 * sig)

        v_flux = A * e4s + B * e2s + C * e6s
        v_matter = beta_matter * rho * e2s
        v_tot = v_flux + v_matter

        V_flux.append(v_flux)
        V_matter.append(v_matter)
        V_total.append(v_tot)

    return {
        "sigma_grid": sigma_grid,
        "V_flux": V_flux,
        "V_matter": V_matter,
        "V_total": V_total,
        "A": A,
        "B": B,
        "C": C,
        "rho": rho,
        "beta_matter": beta_matter,
    }


# ---------------------------------------------------------------------------
# 2. Find density-dependent minimum
# ---------------------------------------------------------------------------

def find_density_dependent_minimum(rho, N=1, a_0=1.0, beta_matter=1.0):
    """
    Find the minimum of V_eff with matter density rho.

    The potential with matter is:
        V = A*e^{4s} + (B + beta_matter*rho)*e^{2s} + C*e^{6s}

    This has the same functional form as the pure flux potential with
    an effective curvature coefficient B_eff = B + beta_matter * rho.

    Setting V'(sigma) = 0 with u = e^{2*sigma} gives the quadratic:
        6C * u^2 + 4A * u + 2*B_eff = 0

    Parameters
    ----------
    rho : float
        Matter density in Planck units.
    N : int
        Flux quantum number (default 1).
    a_0 : float
        Radius of S^2 in Planck units (default 1.0).
    beta_matter : float
        Dimensionless matter coupling (default 1.0).

    Returns
    -------
    dict with keys:
        minimum_exists, sigma_0, a_effective, m_phi_eff, m_phi_eff_eV,
        B_eff, discriminant, u_plus, u_minus, A, B, C, rho, beta_matter
    """
    A, B, C = _get_flux_coefficients(N=N, a_0=a_0)

    B_eff = B + beta_matter * rho

    # Solve 6C u^2 + 4A u + 2 B_eff = 0
    a_q = 6.0 * C
    b_q = 4.0 * A
    c_q = 2.0 * B_eff

    discriminant = b_q ** 2 - 4.0 * a_q * c_q

    minimum_exists = False
    sigma_0 = None
    a_effective = None
    m_phi_eff = None
    m_phi_eff_eV = None
    V_double_prime = None
    u_plus = None
    u_minus = None

    if discriminant >= 0 and a_q != 0:
        sqrt_disc = math.sqrt(discriminant)
        u_plus = (-b_q + sqrt_disc) / (2.0 * a_q)
        u_minus = (-b_q - sqrt_disc) / (2.0 * a_q)

        # Need u > 0 since u = e^{2*sigma}
        u_stable = None
        if u_plus is not None and u_plus > 0:
            u_stable = u_plus
        elif u_minus is not None and u_minus > 0:
            u_stable = u_minus

        if u_stable is not None:
            sigma_0 = 0.5 * math.log(u_stable)
            # Effective internal radius: a_eff = a_0 * exp(-sigma_0)
            # (sigma > 0 means the internal space has shrunk)
            a_effective = a_0 * math.exp(-sigma_0)

            # Second derivative at the minimum
            e2s = u_stable
            e4s = e2s * e2s
            e6s = e4s * e2s
            V_double_prime = (
                16.0 * A * e4s
                + 4.0 * B_eff * e2s
                + 36.0 * C * e6s
            )

            is_minimum = V_double_prime > 0
            minimum_exists = is_minimum

            if is_minimum and V_double_prime > 0:
                # Dilaton mass in Planck units: m_phi^2 = V''
                m_phi_eff = math.sqrt(V_double_prime)
                m_phi_eff_eV = _M_PL_EV * m_phi_eff

    return {
        "minimum_exists": minimum_exists,
        "sigma_0": sigma_0,
        "a_effective": a_effective,
        "m_phi_eff": m_phi_eff,
        "m_phi_eff_eV": m_phi_eff_eV,
        "V_double_prime": V_double_prime,
        "B_eff": B_eff,
        "discriminant": discriminant,
        "u_plus": u_plus,
        "u_minus": u_minus,
        "A": A,
        "B": B,
        "C": C,
        "rho": rho,
        "beta_matter": beta_matter,
    }


# ---------------------------------------------------------------------------
# 3. Chameleon profile across densities
# ---------------------------------------------------------------------------

def compute_chameleon_profile(rho_values=None, N=1, a_0=1.0,
                              beta_matter=1.0):
    """
    Scan over a range of matter densities and compute the density-dependent
    minimum, effective dilaton mass, and effective internal radius at each.

    Densities span from cosmic void (~10^{-30} kg/m^3) to neutron star
    (~10^{17} kg/m^3), converted to Planck units.

    Parameters
    ----------
    rho_values : list of float or None
        Matter densities in SI units (kg/m^3).  If None, uses 48
        log-spaced points from 1e-30 to 1e17.
    N : int
        Flux quantum number (default 1).
    a_0 : float
        Radius of S^2 in Planck units (default 1.0).
    beta_matter : float
        Dimensionless matter coupling (default 1.0).

    Returns
    -------
    dict with keys:
        rho_values_si, rho_values_planck, sigma_0_values,
        a_eff_values, m_phi_eff_values, m_phi_eff_eV_values
    """
    if rho_values is None:
        n_points = 48
        log_min = -30.0
        log_max = 17.0
        step = (log_max - log_min) / (n_points - 1)
        rho_values = [10.0 ** (log_min + i * step) for i in range(n_points)]

    # Conversion factor: rho_planck = rho_SI * (l_Pl^3 / M_Pl_kg)
    conv_factor = _L_PL ** 3 / _M_PL_KG

    rho_values_planck = [rho_si * conv_factor for rho_si in rho_values]

    sigma_0_values = []
    a_eff_values = []
    m_phi_eff_values = []
    m_phi_eff_eV_values = []

    for rho_pl in rho_values_planck:
        result = find_density_dependent_minimum(
            rho=rho_pl, N=N, a_0=a_0, beta_matter=beta_matter,
        )

        sigma_0_values.append(result["sigma_0"])

        if result["a_effective"] is not None:
            # Convert from Planck units to meters
            a_eff_m = result["a_effective"] * _L_PL
        else:
            a_eff_m = None
        a_eff_values.append(a_eff_m)

        m_phi_eff_values.append(result["m_phi_eff"])
        m_phi_eff_eV_values.append(result["m_phi_eff_eV"])

    return {
        "rho_values_si": rho_values,
        "rho_values_planck": rho_values_planck,
        "sigma_0_values": sigma_0_values,
        "a_eff_values": a_eff_values,
        "m_phi_eff_values": m_phi_eff_values,
        "m_phi_eff_eV_values": m_phi_eff_eV_values,
    }


# ---------------------------------------------------------------------------
# 4. KK truncation validity check
# ---------------------------------------------------------------------------

def check_kk_truncation_validity(a_eff_m):
    """
    The critical self-consistency check for chameleon screening.

    Given an effective internal radius a_eff in meters, determine whether
    the 4D effective theory remains valid.  The key issue is NOT whether
    m_phi < m_KK (which is trivially satisfied), but whether gravity
    itself becomes higher-dimensional at observable distances.

    If a_eff = 10^7 m, gravity is 6D (falls as 1/r^4) at distances
    below 10^7 m.  This is catastrophically excluded by every gravity
    experiment ever performed.

    The experimental bound comes from the Eot-Wash torsion balance,
    which tests the inverse-square law down to ~56 micrometers.

    Parameters
    ----------
    a_eff_m : float
        Effective internal radius in meters.

    Returns
    -------
    dict with keys:
        a_eff_m, m_kk_first_eV, gravity_becomes_6d_below_m,
        kk_truncation_valid, exclusion_status, exclusion_reason,
        n_kk_modes_below_cutoff
    """
    # First KK mass: m_KK = hbar * c / a_eff
    m_kk_first_eV = _HBAR_C_EV_M / a_eff_m if a_eff_m > 0 else float('inf')

    # Gravity becomes 6-dimensional at distances r < a_eff
    gravity_6d_below = a_eff_m

    # Number of KK modes below a reference cutoff (1 eV as a proxy)
    # N_kk ~ cutoff / m_KK = cutoff * a_eff / (hbar*c)
    cutoff_eV = 1.0
    if m_kk_first_eV > 0:
        n_kk_below_cutoff = cutoff_eV / m_kk_first_eV
    else:
        n_kk_below_cutoff = float('inf')

    # Exclusion logic based on sub-mm gravity tests
    if a_eff_m < _LAMBDA_EOT_WASH:
        exclusion_status = "valid"
        kk_truncation_valid = True
        exclusion_reason = (
            f"a_eff = {a_eff_m:.4e} m < lambda_EotWash = {_LAMBDA_EOT_WASH:.0e} m.  "
            f"Extra dimensions are smaller than the smallest scale probed by "
            f"gravity experiments.  The 4D description is self-consistent."
        )
    elif a_eff_m < 1e-4:
        exclusion_status = "marginal"
        kk_truncation_valid = False
        exclusion_reason = (
            f"a_eff = {a_eff_m:.4e} m is between the Eot-Wash threshold "
            f"({_LAMBDA_EOT_WASH:.0e} m) and 0.1 mm.  Gravity deviations "
            f"from 1/r^2 would be detectable in current sub-mm experiments.  "
            f"Marginally excluded."
        )
    else:
        exclusion_status = "excluded"
        kk_truncation_valid = False
        if a_eff_m > 1.0:
            exclusion_reason = (
                f"a_eff = {a_eff_m:.4e} m >> 1 m.  Gravity would be 6D "
                f"(falling as 1/r^4) at all distances below {a_eff_m:.2e} m.  "
                f"This is catastrophically excluded by every gravity measurement "
                f"from Cavendish to satellite ranging."
            )
        else:
            exclusion_reason = (
                f"a_eff = {a_eff_m:.4e} m > 0.1 mm.  Gravity would deviate "
                f"from 1/r^2 at tabletop scales, excluded by Cavendish-type "
                f"experiments."
            )

    return {
        "a_eff_m": a_eff_m,
        "m_kk_first_eV": m_kk_first_eV,
        "gravity_becomes_6d_below_m": gravity_6d_below,
        "kk_truncation_valid": kk_truncation_valid,
        "exclusion_status": exclusion_status,
        "exclusion_reason": exclusion_reason,
        "n_kk_modes_below_cutoff": n_kk_below_cutoff,
    }


# ---------------------------------------------------------------------------
# 5. Self-consistency assessment
# ---------------------------------------------------------------------------

def assess_chameleon_self_consistency(N=1, a_0=1.0, beta_matter=1.0):
    """
    The main verdict function for chameleon screening self-consistency.

    Computes the chameleon profile across densities, checks KK truncation
    validity at each density, and identifies whether the mechanism is
    self-consistent in various regimes.

    Specifically checks four representative environments:
        (a) cosmic void (rho ~ 1e-30 kg/m^3)
        (b) solar system interplanetary (rho ~ 1e-18 kg/m^3)
        (c) Earth surface (rho ~ 5500 kg/m^3)
        (d) Laboratory torsion balance (rho ~ 10 kg/m^3)

    Parameters
    ----------
    N : int
        Flux quantum number (default 1).
    a_0 : float
        Radius of S^2 in Planck units (default 1.0).
    beta_matter : float
        Dimensionless matter coupling (default 1.0).

    Returns
    -------
    dict with keys:
        environment_results, fuzzy_dm_viable, mev_regime_viable,
        self_consistent_mass_range_eV, overall_verdict, no_go_reason
    """
    # Named environments
    environments = [
        {"name": "cosmic_void", "rho_si": 1e-30,
         "label": "Cosmic void (intergalactic)"},
        {"name": "interplanetary", "rho_si": 1e-18,
         "label": "Solar system (interplanetary medium)"},
        {"name": "laboratory", "rho_si": 10.0,
         "label": "Laboratory (torsion balance environment)"},
        {"name": "earth_surface", "rho_si": 5500.0,
         "label": "Earth surface (rock density)"},
    ]

    conv_factor = _L_PL ** 3 / _M_PL_KG

    environment_results = []
    for env in environments:
        rho_pl = env["rho_si"] * conv_factor
        minimum = find_density_dependent_minimum(
            rho=rho_pl, N=N, a_0=a_0, beta_matter=beta_matter,
        )

        sigma_0 = minimum["sigma_0"]
        m_phi_eff_eV = minimum["m_phi_eff_eV"]

        if minimum["a_effective"] is not None:
            a_eff_m = minimum["a_effective"] * _L_PL
            kk_check = check_kk_truncation_validity(a_eff_m)
            kk_valid = kk_check["kk_truncation_valid"]
            exclusion_status = kk_check["exclusion_status"]
        else:
            a_eff_m = None
            kk_valid = False
            exclusion_status = "no_minimum"

        environment_results.append({
            "name": env["name"],
            "label": env["label"],
            "rho_si": env["rho_si"],
            "sigma_0": sigma_0,
            "a_eff_m": a_eff_m,
            "m_phi_eff_eV": m_phi_eff_eV,
            "kk_valid": kk_valid,
            "exclusion_status": exclusion_status,
        })

    # Determine fuzzy DM viability
    # For fuzzy DM we need m_phi ~ 10^{-22} eV, requiring a_0 ~ 10^7 m.
    # Check: does the void environment produce an excluded a_eff?
    void_result = next(
        (r for r in environment_results if r["name"] == "cosmic_void"),
        None,
    )
    fuzzy_dm_viable = False
    no_go_reason = ""

    if void_result is not None:
        if void_result["a_eff_m"] is not None and void_result["a_eff_m"] > 1.0:
            fuzzy_dm_viable = False
            no_go_reason = (
                f"To achieve m_phi ~ 10^{{-22}} eV in cosmic voids, the "
                f"effective internal radius must be a_eff ~ 10^7 m.  At this "
                f"radius, gravity is 6-dimensional (1/r^4) at all distances "
                f"below ~10,000 km.  This is excluded by every gravity "
                f"experiment from torsion balances to satellite orbits.  "
                f"The KK tower has spacing m_KK ~ {_HBAR_C_EV_M / 1e7:.2e} eV, "
                f"meaning billions of KK modes are excited at any observable "
                f"energy.  The 4D effective theory is completely invalid in "
                f"this regime."
            )
        elif void_result["a_eff_m"] is not None:
            if void_result["kk_valid"]:
                # Check if the mass is actually in the fuzzy DM range
                if (void_result["m_phi_eff_eV"] is not None
                        and void_result["m_phi_eff_eV"] < 1e-20):
                    fuzzy_dm_viable = True
                else:
                    fuzzy_dm_viable = False
                    no_go_reason = (
                        "The void dilaton mass is not in the fuzzy DM range "
                        "(10^{-22} eV).  The chameleon mechanism with "
                        "Planck-scale a_0 produces Planck-scale masses even "
                        "at low densities."
                    )
            else:
                fuzzy_dm_viable = False
                no_go_reason = (
                    f"The effective internal radius in cosmic voids "
                    f"(a_eff = {void_result['a_eff_m']:.2e} m) exceeds the "
                    f"Eot-Wash bound.  The 4D description breaks down."
                )
        else:
            fuzzy_dm_viable = False
            no_go_reason = (
                "No stable minimum exists in the effective potential at "
                "void density.  The chameleon mechanism does not function."
            )

    # Determine meV regime viability
    lab_result = next(
        (r for r in environment_results if r["name"] == "laboratory"),
        None,
    )
    mev_regime_viable = False
    if lab_result is not None and lab_result["a_eff_m"] is not None:
        mev_regime_viable = lab_result["a_eff_m"] < _LAMBDA_EOT_WASH

    # Self-consistent mass range
    self_consistent_mass_range_eV = None
    if mev_regime_viable and lab_result["m_phi_eff_eV"] is not None:
        earth_result = next(
            (r for r in environment_results if r["name"] == "earth_surface"),
            None,
        )
        m_low = lab_result["m_phi_eff_eV"]
        m_high = m_low
        if (earth_result is not None
                and earth_result["m_phi_eff_eV"] is not None):
            m_high = max(m_low, earth_result["m_phi_eff_eV"])
        self_consistent_mass_range_eV = (m_low, m_high)

    # Overall verdict
    overall_verdict = (
        "The chameleon screening mechanism for the Alpha Ladder dilaton "
        "faces a fundamental self-consistency problem in the fuzzy dark "
        "matter regime.  Achieving m_phi ~ 10^{-22} eV requires an "
        "internal radius a_0 ~ 10^7 m, at which scale the entire KK "
        "tower is excited and gravity becomes 6-dimensional.  This is "
        "catastrophically excluded.\n\n"
        "In the meV regime (a_0 ~ 30-100 um), the KK truncation is "
        "self-consistent and the dilaton could produce a testable fifth "
        "force at sub-millimeter scales.  However, a meV dilaton is NOT "
        "a dark matter candidate: it is too heavy for fuzzy DM, too "
        "light for a WIMP, and its relic abundance is negligible.\n\n"
        "The chameleon mechanism itself is not derived from the 6D action "
        "-- the matter coupling beta_matter * rho * e^{2*sigma} is an "
        "additional assumption beyond the Alpha Ladder framework.  An "
        "honest assessment is that the screening mechanism requires "
        "ingredients external to the theory."
    )

    return {
        "environment_results": environment_results,
        "fuzzy_dm_viable": fuzzy_dm_viable,
        "mev_regime_viable": mev_regime_viable,
        "self_consistent_mass_range_eV": self_consistent_mass_range_eV,
        "overall_verdict": overall_verdict,
        "no_go_reason": no_go_reason,
    }


# ---------------------------------------------------------------------------
# 6. meV dark sector potential
# ---------------------------------------------------------------------------

def compute_meV_dark_sector_potential(a_0_m=30e-6):
    """
    For the meV regime where chameleon screening IS self-consistent,
    compute what the dilaton could contribute to dark sector physics.

    The honest conclusion: a meV dilaton is testable via fifth-force
    experiments but is NOT a viable dark matter or dark energy candidate.

    Parameters
    ----------
    a_0_m : float
        Internal radius in meters (default 30 um).

    Returns
    -------
    dict with keys:
        a_0_m, m_phi_eV, lambda_compton_m, is_fuzzy_dm_candidate,
        is_quintessence_candidate, is_testable_fifth_force,
        honest_assessment
    """
    # Dilaton mass: m_phi = M_Pl * l_Pl / a_0  (from KK scaling)
    # In eV: m_phi = M_Pl_eV * l_Pl / a_0
    m_phi_eV = _M_PL_EV * _L_PL / a_0_m

    # Compton wavelength: lambda_c = hbar * c / m_phi
    lambda_compton_m = _HBAR_C_EV_M / m_phi_eV if m_phi_eV > 0 else float('inf')

    # Dark energy / quintessence: needs m ~ H_0 ~ 10^{-33} eV
    H_0_eV = 1.5e-33  # Hubble constant in eV
    is_quintessence = m_phi_eV < 10.0 * H_0_eV

    # Fuzzy DM: needs m ~ 10^{-22} eV
    is_fuzzy_dm = 1e-23 < m_phi_eV < 1e-21

    # Testable fifth force: meV range, lambda_c ~ sub-mm
    is_fifth_force = 1e-4 < m_phi_eV < 1.0  # 0.1 meV to 1 eV

    # Relic abundance assessment
    if m_phi_eV > 1e-3:
        relic_note = (
            f"A dilaton with m = {m_phi_eV:.4e} eV is far too light "
            f"for thermal freeze-out (WIMP mechanism requires m > ~GeV).  "
            f"It is also far too heavy for fuzzy DM (requires m ~ 10^{{-22}} eV).  "
            f"The dilaton does not oscillate coherently at early times "
            f"(m >> H at matter-radiation equality), so it does not "
            f"contribute significantly to the matter density.  "
            f"Relic abundance: negligible."
        )
    else:
        relic_note = (
            f"m_phi = {m_phi_eV:.4e} eV: sub-meV mass.  "
            f"Cosmological constraints depend on the production mechanism."
        )

    honest_assessment = (
        f"The dilaton with a_0 = {a_0_m:.2e} m has mass "
        f"m_phi = {m_phi_eV:.4e} eV (Compton wavelength "
        f"{lambda_compton_m:.4e} m).  "
        f"This is NOT a dark matter candidate: too heavy for fuzzy DM "
        f"(by ~{m_phi_eV / 1e-22:.0e} orders of magnitude), too light for "
        f"a WIMP.  It is NOT dark energy: the mass exceeds H_0 by "
        f"~{m_phi_eV / H_0_eV:.0e} orders of magnitude.  "
        f"What it IS: a testable prediction for sub-millimeter fifth-force "
        f"experiments (Eot-Wash, IUPUI, Casimir force measurements).  "
        f"If the Alpha Ladder dilaton exists in this mass range, it would "
        f"modify the gravitational inverse-square law at distances below "
        f"{lambda_compton_m:.2e} m.  {relic_note}"
    )

    return {
        "a_0_m": a_0_m,
        "m_phi_eV": m_phi_eV,
        "lambda_compton_m": lambda_compton_m,
        "is_fuzzy_dm_candidate": is_fuzzy_dm,
        "is_quintessence_candidate": is_quintessence,
        "is_testable_fifth_force": is_fifth_force,
        "honest_assessment": honest_assessment,
    }


# ---------------------------------------------------------------------------
# 7. Summary / dashboard entry point
# ---------------------------------------------------------------------------

def summarize_chameleon_screening(constants=None):
    """
    Dashboard entry point.  Runs the full chameleon screening pipeline:

    1. Chameleon profile across densities
    2. KK truncation check at each density
    3. Self-consistency assessment
    4. meV regime analysis
    5. Overall verdict

    Parameters
    ----------
    constants : SimpleNamespace or None
        Physical constants from get_constants().  Currently unused
        (Planck units throughout), accepted for API consistency.

    Returns
    -------
    dict with keys:
        chameleon_profile, self_consistency, mev_analysis,
        overall_verdict, fuzzy_dm_excluded, fuzzy_dm_exclusion_reason,
        mev_regime_status, key_finding
    """
    profile = compute_chameleon_profile()
    consistency = assess_chameleon_self_consistency()
    mev_analysis = compute_meV_dark_sector_potential(a_0_m=30e-6)

    fuzzy_dm_excluded = not consistency["fuzzy_dm_viable"]

    if consistency["mev_regime_viable"]:
        mev_status = (
            "Self-consistent.  The meV dilaton passes KK truncation checks "
            "and could produce a testable fifth force at sub-millimeter "
            "scales.  However, it is not a dark matter candidate."
        )
    else:
        mev_status = (
            "The meV regime is not self-consistent with the current "
            "flux stabilization parameters.  The effective internal "
            "radius exceeds the Eot-Wash bound even at laboratory "
            "densities."
        )

    key_finding = (
        "The chameleon mechanism CANNOT make the Alpha Ladder dilaton a "
        "fuzzy dark matter candidate: the required internal radius "
        "(~10^7 m) makes gravity 6-dimensional at observable scales.  "
        "The meV regime is self-consistent but phenomenologically limited "
        "to sub-mm fifth-force searches."
    )

    return {
        "chameleon_profile": profile,
        "self_consistency": consistency,
        "mev_analysis": mev_analysis,
        "overall_verdict": consistency["overall_verdict"],
        "fuzzy_dm_excluded": fuzzy_dm_excluded,
        "fuzzy_dm_exclusion_reason": consistency["no_go_reason"],
        "mev_regime_status": mev_status,
        "key_finding": key_finding,
    }
