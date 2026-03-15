"""
Solar system compatibility checks for the massive dilaton.

Computes PPN parameters, fifth-force exclusion bounds, and minimum
dilaton mass required to pass the Cassini constraint on |gamma - 1|.
"""

import math


def compute_minimum_dilaton_mass(omega_BD):
    """Compute the minimum dilaton mass to pass the Cassini bound."""
    hbar = 1.054571817e-34
    c = 2.99792458e8
    eV_to_J = 1.602176634e-19
    AU = 1.496e11

    denom = 2.0 * omega_BD + 3.0
    alpha_fifth = 2.0 / denom if denom > 0 else 0.0
    cassini_limit = 2.3e-5

    if alpha_fifth > cassini_limit:
        suppression_factor = math.log(alpha_fifth / cassini_limit)
        lambda_max = AU / suppression_factor
    else:
        suppression_factor = 0.0
        lambda_max = AU

    m_phi_min_kg = hbar / (c * lambda_max)
    m_phi_min_eV = m_phi_min_kg * c**2 / eV_to_J

    return {
        "m_phi_min_eV": m_phi_min_eV,
        "m_phi_min_kg": m_phi_min_kg,
        "lambda_max_m": lambda_max,
        "lambda_max_au": lambda_max / AU,
        "suppression_factor": suppression_factor,
    }


def compute_ppn_parameters(omega_BD, m_phi_eV):
    """Compute PPN parameters for a massive Brans-Dicke dilaton."""
    hbar = 1.054571817e-34
    c = 2.99792458e8
    eV_to_J = 1.602176634e-19
    AU = 1.496e11

    denom = 2.0 * omega_BD + 3.0
    beta_coupling = 1.0 / math.sqrt(denom) if denom > 0 else float('inf')
    alpha_fifth = 2.0 * beta_coupling**2  # = 2/(2w+3)

    gamma_PPN_massless = (1.0 + omega_BD) / (2.0 + omega_BD)

    # Compton wavelength
    m_phi_kg = m_phi_eV * eV_to_J / c**2
    lambda_compton_m = hbar / (m_phi_kg * c)

    # gamma_PPN(r) = 1 - alpha_fifth * exp(-r / lambda)
    r_cassini = AU
    gamma_PPN_at_cassini = 1.0 - alpha_fifth * math.exp(-r_cassini / lambda_compton_m)
    gamma_deviation = abs(gamma_PPN_at_cassini - 1.0)
    cassini_bound = 2.3e-5
    passes_cassini = gamma_deviation < cassini_bound

    # Required minimum mass
    min_mass_result = compute_minimum_dilaton_mass(omega_BD)

    return {
        "omega_BD": omega_BD,
        "m_phi_eV": m_phi_eV,
        "beta_coupling": beta_coupling,
        "alpha_fifth": alpha_fifth,
        "lambda_compton_m": lambda_compton_m,
        "gamma_PPN_massless": gamma_PPN_massless,
        "gamma_PPN_at_cassini": gamma_PPN_at_cassini,
        "gamma_deviation_at_cassini": gamma_deviation,
        "cassini_bound": cassini_bound,
        "passes_cassini": passes_cassini,
        "required_min_mass_eV": min_mass_result["m_phi_min_eV"],
    }


def compute_ppn_profile(omega_BD, m_phi_eV, n_points=500):
    """Compute PPN gamma as a function of distance."""
    hbar = 1.054571817e-34
    c = 2.99792458e8
    eV_to_J = 1.602176634e-19
    AU = 1.496e11

    denom = 2.0 * omega_BD + 3.0
    alpha_fifth = 2.0 / denom if denom > 0 else 0.0
    m_phi_kg = m_phi_eV * eV_to_J / c**2
    lambda_compton = hbar / (m_phi_kg * c)

    r_min = 0.01
    r_max = 10.0 * AU
    log_min = math.log10(r_min)
    log_max = math.log10(r_max)
    r_meters = [10 ** (log_min + i * (log_max - log_min) / (n_points - 1)) for i in range(n_points)]

    gamma_PPN = [1.0 - alpha_fifth * math.exp(-r / lambda_compton) for r in r_meters]
    gamma_deviation = [abs(g - 1.0) for g in gamma_PPN]

    landmarks = {
        "lab": {"r_meters": 0.1, "label": "Lab (0.1 m)"},
        "Earth_radius": {"r_meters": 6.371e6, "label": "Earth radius"},
        "Moon_orbit": {"r_meters": 3.844e8, "label": "Moon orbit"},
        "Mercury_orbit": {"r_meters": 5.791e10, "label": "Mercury orbit"},
        "1_AU": {"r_meters": AU, "label": "1 AU (Cassini)"},
    }
    for key in landmarks:
        r = landmarks[key]["r_meters"]
        landmarks[key]["gamma_PPN"] = 1.0 - alpha_fifth * math.exp(-r / lambda_compton)
        landmarks[key]["gamma_deviation"] = abs(landmarks[key]["gamma_PPN"] - 1.0)

    cassini_bound = 2.3e-5
    transition_radius_m = None
    for r, dev in zip(r_meters, gamma_deviation):
        if dev < cassini_bound:
            transition_radius_m = r
            break

    return {
        "r_meters": r_meters,
        "gamma_PPN": gamma_PPN,
        "gamma_deviation": gamma_deviation,
        "landmarks": landmarks,
        "transition_radius_m": transition_radius_m,
    }


def compute_fifth_force_point(omega_BD, m_phi_eV):
    """Compute the dilaton's position in (alpha, lambda) exclusion space."""
    hbar = 1.054571817e-34
    c = 2.99792458e8
    eV_to_J = 1.602176634e-19

    denom = 2.0 * omega_BD + 3.0
    alpha_fifth = 2.0 / denom if denom > 0 else 0.0
    m_phi_kg = m_phi_eV * eV_to_J / c**2
    lambda_fifth = hbar / (m_phi_kg * c)

    return {
        "alpha_fifth": alpha_fifth,
        "lambda_fifth_m": lambda_fifth,
        "log10_alpha": math.log10(alpha_fifth) if alpha_fifth > 0 else float('-inf'),
        "log10_lambda_m": math.log10(lambda_fifth) if lambda_fifth > 0 else float('-inf'),
    }


def generate_exclusion_bounds():
    """Return approximate exclusion boundaries from 5 fifth-force experiments.

    Each boundary is a list of [log10(lambda/m), log10(alpha)] points.
    The excluded region is above the boundary curve.

    For Cassini, the boundary is computed from the actual constraint:
        alpha < |gamma-1|_max * exp(r_Cassini / lambda)
    which rises steeply for lambda << AU, making short-range forces allowed.
    """
    AU = 1.496e11
    cassini_limit = 2.3e-5

    # Cassini boundary from the actual Yukawa constraint formula
    # alpha_excluded(lambda) = cassini_limit * exp(AU / lambda)
    cassini_boundary = []
    for log_lam in [9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 13.0, 14.0]:
        lam = 10**log_lam
        log_alpha = math.log10(cassini_limit) + AU / lam / math.log(10)
        cassini_boundary.append([log_lam, log_alpha])

    bounds = [
        {
            "name": "Eot-Wash",
            "boundary": [
                [-4.0, 0.0], [-3.5, -1.0], [-3.0, -2.0],
                [-2.5, -2.5], [-2.0, -2.5], [-1.5, -2.0],
            ],
            "excluded_side": "above",
            "reference": "Adelberger et al., Prog. Part. Nucl. Phys. 62 (2009) 102",
        },
        {
            "name": "MICROSCOPE",
            "boundary": [
                [5.0, -12.0], [6.0, -12.5], [7.0, -13.0],
                [8.0, -13.0], [9.0, -12.5],
            ],
            "excluded_side": "above",
            "reference": "Touboul et al., Phys. Rev. Lett. 129 (2022) 121102",
        },
        {
            "name": "Lunar Laser Ranging",
            "boundary": [
                [8.0, -10.0], [8.5, -12.0], [9.0, -13.0],
                [9.5, -13.5],
            ],
            "excluded_side": "above",
            "reference": "Williams et al., Class. Quantum Grav. 29 (2012) 184004",
        },
        {
            "name": "Cassini",
            "boundary": cassini_boundary,
            "excluded_side": "above",
            "reference": "Bertotti et al., Nature 425 (2003) 374",
        },
        {
            "name": "Planetary ephemerides",
            "boundary": [
                [12.0, -11.0], [13.0, -11.5], [14.0, -11.0],
                [15.0, -10.0],
            ],
            "excluded_side": "above",
            "reference": "Fienga et al., Celest. Mech. Dyn. Astron. 123 (2015) 325",
        },
    ]
    return bounds


def check_dilaton_exclusion(omega_BD, m_phi_eV):
    """Check whether the dilaton point lies in the allowed region."""
    point = compute_fifth_force_point(omega_BD, m_phi_eV)
    bounds = generate_exclusion_bounds()

    log_alpha = point["log10_alpha"]
    log_lambda = point["log10_lambda_m"]

    excluded_by = []
    closest_bound = None
    closest_margin = float('inf')

    for bound in bounds:
        boundary = bound["boundary"]
        # Check if the dilaton lambda falls within this experiment's range
        lam_values = [b[0] for b in boundary]
        lam_min, lam_max = min(lam_values), max(lam_values)

        if lam_min <= log_lambda <= lam_max:
            # Interpolate the boundary alpha at this lambda
            for i in range(len(boundary) - 1):
                l0, a0 = boundary[i]
                l1, a1 = boundary[i + 1]
                if l0 <= log_lambda <= l1:
                    t = (log_lambda - l0) / (l1 - l0) if l1 != l0 else 0
                    alpha_bound = a0 + t * (a1 - a0)
                    margin = alpha_bound - log_alpha
                    if log_alpha > alpha_bound:
                        excluded_by.append(bound["name"])
                    if abs(margin) < abs(closest_margin):
                        closest_margin = margin
                        closest_bound = bound["name"]
                    break

    return {
        "dilaton_point": point,
        "bounds": bounds,
        "excluded_by": excluded_by,
        "allowed": len(excluded_by) == 0,
        "closest_bound": {"name": closest_bound, "margin_dex": closest_margin} if closest_bound else None,
    }


def summarize_solar_system(constants):
    """Top-level entry point: assemble all solar system compatibility results."""
    from alpha_ladder_core.dilaton import compute_bd_parameter

    bd = compute_bd_parameter(constants)
    omega = bd["omega"]
    m_phi_eV = bd["dilaton_mass_min_eV"]

    ppn = compute_ppn_parameters(omega, m_phi_eV)
    fifth_force = check_dilaton_exclusion(omega, m_phi_eV)
    minimum_mass = compute_minimum_dilaton_mass(omega)

    passes_all = ppn["passes_cassini"] and fifth_force["allowed"]

    summary = (
        f"Solar system compatibility: "
        f"gamma deviation at 1 AU = {ppn['gamma_deviation_at_cassini']:.2e} "
        f"({'PASSES' if ppn['passes_cassini'] else 'FAILS'} Cassini bound of {ppn['cassini_bound']:.1e}). "
        f"Fifth-force exclusion: {'allowed' if fifth_force['allowed'] else 'excluded by ' + ', '.join(fifth_force['excluded_by'])}."
    )

    return {
        "ppn": ppn,
        "fifth_force": fifth_force,
        "minimum_mass": minimum_mass,
        "passes_all": passes_all,
        "summary": summary,
    }
