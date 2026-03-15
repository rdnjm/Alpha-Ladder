"""
Yukawa screening model for Newton's G.

The dilaton field introduces a distance-dependent correction to G:
at lab scales (r << lambda_dilaton), G_eff is slightly higher than
the vacuum prediction; at solar scales, screening suppresses
the excess and G converges to G_vacuum.

Refactored from Pillar 4 analysis.
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50


def compute_screening_parameters(constants):
    """Calibrate the Yukawa screening model from first-principles prediction and CODATA lab value.

    Parameters
    ----------
    constants : SimpleNamespace
        Physical and mathematical constants from get_constants().

    Returns
    -------
    dict with keys:
        G_vacuum (float): predicted G from phi^2/2 bridge
        G_lab (float): CODATA 2018 recommended G
        alpha_screening (float): fractional excess (G_lab - G_vacuum) / G_vacuum
        ppm_excess (float): alpha_screening * 1e6
        lambda_dilaton_m (float): dilaton Compton wavelength in meters = hbar / (m_dilaton_min * c)
        dilaton_mass_kg (float): minimum dilaton mass from Cassini bound
        dilaton_mass_eV (float): same in eV
        omega (float): Brans-Dicke parameter
    """
    from alpha_ladder_core.predict_g import predict_G, get_bridge_candidates, get_G_measurements
    from alpha_ladder_core.dilaton import compute_bd_parameter

    # G_vacuum from phi^2/2 bridge
    bridges = get_bridge_candidates(constants)
    G_vacuum = float(predict_G(bridges["\u03c6\u00b2/2"], constants))

    # G_lab from CODATA 2018 recommended value
    measurements = get_G_measurements()
    G_lab = float(measurements["CODATA 2018 recommended"][0])

    # Screening amplitude
    alpha_screening = (G_lab - G_vacuum) / G_vacuum
    ppm_excess = alpha_screening * 1e6

    # Dilaton parameters from existing module
    bd = compute_bd_parameter(constants)
    dilaton_mass_kg = bd["dilaton_mass_min_kg"]
    dilaton_mass_eV = bd["dilaton_mass_min_eV"]
    omega = bd["omega"]

    # Dilaton Compton wavelength: lambda = hbar / (m * c)
    hbar_f = float(constants.hbar)
    c_f = float(constants.c)
    lambda_dilaton_m = hbar_f / (dilaton_mass_kg * c_f)

    return {
        "G_vacuum": G_vacuum,
        "G_lab": G_lab,
        "alpha_screening": alpha_screening,
        "ppm_excess": ppm_excess,
        "lambda_dilaton_m": lambda_dilaton_m,
        "dilaton_mass_kg": dilaton_mass_kg,
        "dilaton_mass_eV": dilaton_mass_eV,
        "omega": omega,
    }


def compute_G_eff(r_meters, screening_params):
    """Compute effective G at distance r using Yukawa screening.

    G_eff(r) = G_vacuum * (1 + alpha_screening * exp(-r / lambda_dilaton))

    Parameters
    ----------
    r_meters : float
        Distance in meters.
    screening_params : dict
        Output of compute_screening_parameters().

    Returns
    -------
    float
        G_eff at the given distance.
    """
    G_vac = screening_params["G_vacuum"]
    alpha_s = screening_params["alpha_screening"]
    lam = screening_params["lambda_dilaton_m"]

    return G_vac * (1.0 + alpha_s * math.exp(-r_meters / lam))


def compute_screening_profile(screening_params, n_points=500):
    """Compute G_eff over log-spaced distances from 0.01 m to 10 AU.

    Parameters
    ----------
    screening_params : dict
        Output of compute_screening_parameters().
    n_points : int
        Number of sample points.

    Returns
    -------
    dict with keys:
        r_meters (list[float]): distances
        G_eff (list[float]): effective G at each distance
        landmarks (dict): named distances with their G_eff values
    """
    AU = 1.496e11  # meters
    r_min = 0.01
    r_max = 10.0 * AU

    # Log-spaced distances
    log_min = math.log10(r_min)
    log_max = math.log10(r_max)
    r_meters = [10 ** (log_min + i * (log_max - log_min) / (n_points - 1)) for i in range(n_points)]

    G_eff_values = [compute_G_eff(r, screening_params) for r in r_meters]

    # Named landmarks
    landmark_defs = {
        "Lab bench (1 m)": 1.0,
        "Earth radius (6.4e6 m)": 6.4e6,
        "Moon orbit (3.84e8 m)": 3.84e8,
        "1 AU (1.496e11 m)": AU,
    }
    landmarks = {}
    for name, dist in landmark_defs.items():
        landmarks[name] = {
            "r_meters": dist,
            "G_eff": compute_G_eff(dist, screening_params),
        }

    return {
        "r_meters": r_meters,
        "G_eff": G_eff_values,
        "landmarks": landmarks,
    }


def classify_measurements(screening_params, constants):
    """Classify the 7 G measurements into high (lab) and low (atom interferometry) clusters.

    Parameters
    ----------
    screening_params : dict
        Output of compute_screening_parameters().
    constants : SimpleNamespace
        Physical constants.

    Returns
    -------
    dict with keys:
        high_cluster (list[dict]): lab/torsion-balance experiments
        low_cluster (list[dict]): atom interferometry experiments
        cluster_gap_ppm (float): gap between weighted cluster means in ppm
        measurements (list[dict]): all 7 with added 'cluster' tag
    """
    from alpha_ladder_core.predict_g import get_G_measurements, compare_prediction

    G_pred = Decimal(str(screening_params["G_vacuum"]))
    raw = get_G_measurements()
    comparisons = compare_prediction(G_pred, raw)

    # Tag: Rosi is the only atom interferometry experiment
    atom_interf_keywords = ["atom interf", "Rosi"]

    all_measurements = []
    high_cluster = []
    low_cluster = []

    for comp in comparisons:
        name = comp["experiment"]
        entry = {
            "experiment": name,
            "G_exp": float(comp["G_exp"]),
            "G_unc": float(comp["G_unc"]),
            "sigma": comp["sigma"],
            "direction": comp["direction"],
        }

        if any(kw in name for kw in atom_interf_keywords):
            entry["cluster"] = "low"
            low_cluster.append(entry)
        else:
            entry["cluster"] = "high"
            high_cluster.append(entry)

        all_measurements.append(entry)

    # Weighted means for each cluster
    def weighted_mean(cluster):
        if not cluster:
            return 0.0
        total_w = sum(1.0 / (m["G_unc"] ** 2) for m in cluster)
        return sum(m["G_exp"] / (m["G_unc"] ** 2) for m in cluster) / total_w

    high_mean = weighted_mean(high_cluster)
    low_mean = weighted_mean(low_cluster)

    G_vac = screening_params["G_vacuum"]
    cluster_gap_ppm = abs(high_mean - low_mean) / G_vac * 1e6 if G_vac != 0 else 0.0

    return {
        "high_cluster": high_cluster,
        "low_cluster": low_cluster,
        "high_mean": high_mean,
        "low_mean": low_mean,
        "cluster_gap_ppm": cluster_gap_ppm,
        "measurements": all_measurements,
    }


# ---------------------------------------------------------------------------
# Field-theoretic derivation of the Yukawa screening profile
# ---------------------------------------------------------------------------

def derive_dilaton_lagrangian(omega_BD, m_phi_eV=None):
    """Derive the 4D effective Lagrangian for the dilaton after KK reduction.

    The action in Einstein frame is:

        S = integral d^4x sqrt(-g) [
            (1/2) M_Pl^2 R
          - (1/2) (d phi)^2
          - V(phi)
          + L_matter(A^2(phi) g_mu_nu, psi)
        ]

    where phi is the canonically normalized dilaton field and
    A(phi) = exp(beta * phi / M_Pl) is the conformal coupling function,
    with beta = 1/sqrt(2*omega_BD + 3).

    Parameters
    ----------
    omega_BD : float
        Brans-Dicke parameter from the KK reduction or dilaton module.
    m_phi_eV : float or None
        Dilaton mass in eV.  If None, the mass term is left symbolic.

    Returns
    -------
    dict with keys:
        action_terms : list of dicts with name, expression, description
        beta_coupling : float
        conformal_function : str
        omega_BD : float
        canonical_normalization : dict
    """
    from alpha_ladder_core.kk_reduction import (
        compute_kinetic_coefficient,
        compute_einstein_frame_ansatz,
    )

    # Universal coupling constant
    denom = 2.0 * omega_BD + 3.0
    if denom <= 0:
        beta_coupling = float('inf')
    else:
        beta_coupling = 1.0 / math.sqrt(denom)

    # Potential term string
    if m_phi_eV is not None:
        potential_expr = f"-(1/2) * ({m_phi_eV} eV)^2 * phi^2"
        potential_desc = (
            f"Dilaton mass potential with m_phi = {m_phi_eV} eV "
            f"from moduli stabilization"
        )
    else:
        potential_expr = "-(1/2) * m_phi^2 * phi^2"
        potential_desc = (
            "Dilaton mass potential (mass term from moduli stabilization, "
            "value depends on stabilization mechanism)"
        )

    action_terms = [
        {
            "name": "Einstein-Hilbert",
            "expression": "(1/2) M_Pl^2 R",
            "description": (
                "Standard 4D gravity kinetic term; M_Pl = reduced Planck mass"
            ),
        },
        {
            "name": "Dilaton kinetic",
            "expression": "-(1/2) (d phi)^2",
            "description": (
                "Canonically normalized dilaton kinetic term; "
                "phi = sqrt(|K_einstein|) * sigma where sigma is the "
                "breathing mode from KK reduction"
            ),
        },
        {
            "name": "Dilaton potential",
            "expression": potential_expr,
            "description": potential_desc,
        },
        {
            "name": "Matter coupling",
            "expression": "L_matter(A^2(phi) g_mu_nu, psi)",
            "description": (
                f"Matter fields psi couple to the Jordan-frame metric "
                f"A^2(phi) g_mu_nu with A(phi) = exp(beta * phi / M_Pl), "
                f"beta = {beta_coupling:.10f}"
            ),
        },
    ]

    conformal_function = (
        f"A(phi) = exp({beta_coupling:.10f} * phi / M_Pl)"
    )

    # Canonical normalization from KK reduction
    kk = compute_kinetic_coefficient(d=4, n=2)
    K_einstein = kk["K_einstein"]
    sqrt_abs_K = math.sqrt(abs(K_einstein))
    canonical_normalization = {
        "relation": f"phi = sqrt(|K_einstein|) * sigma = {sqrt_abs_K:.10f} * sigma",
        "K_einstein": K_einstein,
        "sqrt_abs_K": sqrt_abs_K,
        "description": (
            "sigma is the breathing mode of the internal manifold; "
            "phi is the canonically normalized field that appears in "
            "the 4D effective action with unit kinetic coefficient"
        ),
    }

    return {
        "action_terms": action_terms,
        "beta_coupling": beta_coupling,
        "conformal_function": conformal_function,
        "omega_BD": omega_BD,
        "canonical_normalization": canonical_normalization,
    }


def derive_field_equation(omega_BD):
    """Derive the linearized field equation for the dilaton in the weak-field limit.

    In the presence of a point source of mass M at the origin, the canonically
    normalized dilaton phi satisfies the massive Klein-Gordon equation:

        (nabla^2 - m_phi^2) phi = -beta * rho / M_Pl

    where rho = M * delta^3(r) is the source density and
    beta = 1/sqrt(2*omega_BD + 3) is the universal coupling.

    The solution is the Green's function of the massive Klein-Gordon operator:

        phi(r) = -(beta * M) / (4 pi M_Pl r) * exp(-m_phi * r)

    Parameters
    ----------
    omega_BD : float
        Brans-Dicke parameter.

    Returns
    -------
    dict with keys:
        equation : str
        beta_coupling : float
        source_term : str
        green_function : dict with phi_r, yukawa_range
        linearization_condition : str
    """
    denom = 2.0 * omega_BD + 3.0
    if denom <= 0:
        beta_coupling = float('inf')
    else:
        beta_coupling = 1.0 / math.sqrt(denom)

    equation = (
        "(nabla^2 - m_phi^2) phi = -beta * rho / M_Pl"
    )

    source_term = (
        "rho = M * delta^3(r) for a point mass M at the origin; "
        "this is the trace of the matter stress-energy tensor "
        "in the non-relativistic limit"
    )

    green_function = {
        "phi_r": (
            f"phi(r) = -({beta_coupling:.10f} * M) / (4 pi M_Pl r) "
            f"* exp(-m_phi * r)"
        ),
        "yukawa_range": (
            "lambda = hbar / (m_phi * c) = 1/m_phi (natural units); "
            "the dilaton-mediated force has finite range set by the "
            "dilaton Compton wavelength"
        ),
        "derivation": (
            "This is the standard Green's function of the operator "
            "(nabla^2 - m^2) in 3 spatial dimensions: "
            "G(r) = -exp(-m*r) / (4 pi r).  The source strength is "
            "beta * M / M_Pl from the conformal coupling A(phi)."
        ),
        "beta_coupling": beta_coupling,
    }

    linearization_condition = (
        "phi << M_Pl (weak field limit); equivalently, "
        "beta * M / (4 pi M_Pl r) << M_Pl, i.e. r >> beta * r_S / (8 pi) "
        "where r_S = 2 G M / c^2 is the Schwarzschild radius"
    )

    return {
        "equation": equation,
        "beta_coupling": beta_coupling,
        "source_term": source_term,
        "green_function": green_function,
        "linearization_condition": linearization_condition,
    }


def derive_yukawa_profile(omega_BD):
    """Derive G_eff(r) from first principles via the dilaton field equation.

    Starting from the dilaton Lagrangian and solving the linearized field
    equation for a point source, this function traces through each step
    to arrive at the Yukawa-modified Newton's constant:

        G_eff(r) = G_N * (1 + 2 * beta^2 * exp(-m_phi * r))

    and identifies the screening parameters alpha_s = 2*beta^2,
    lambda = 1/m_phi.

    Parameters
    ----------
    omega_BD : float
        Brans-Dicke parameter.

    Returns
    -------
    dict with keys:
        derivation_steps : list of dicts (step_number, description, formula, result)
        beta_coupling : float
        alpha_screening_from_lagrangian : float
        yukawa_formula : str
        consistency_check : dict
    """
    denom = 2.0 * omega_BD + 3.0
    if denom <= 0:
        beta_coupling = float('inf')
    else:
        beta_coupling = 1.0 / math.sqrt(denom)

    alpha_screening_lagrangian = 2.0 * beta_coupling ** 2

    derivation_steps = [
        {
            "step_number": 1,
            "description": "Start with the dilaton field equation in Einstein frame",
            "formula": "(nabla^2 - m_phi^2) phi = -beta * rho / M_Pl",
            "result": (
                f"Linearized massive Klein-Gordon equation with coupling "
                f"beta = 1/sqrt(2*omega + 3) = {beta_coupling:.10f}"
            ),
        },
        {
            "step_number": 2,
            "description": "Solve Green's function for a point source M at origin",
            "formula": (
                "phi(r) = -(beta * M) / (4 pi M_Pl r) * exp(-m_phi * r)"
            ),
            "result": (
                "Yukawa-suppressed dilaton profile; the field falls off "
                "exponentially beyond the Compton wavelength lambda = 1/m_phi"
            ),
        },
        {
            "step_number": 3,
            "description": (
                "Compute the metric perturbation from the dilaton profile"
            ),
            "formula": (
                "The dilaton contributes to the gravitational potential via "
                "the conformal coupling A(phi) = exp(beta * phi / M_Pl).  "
                "The Newtonian potential receives an additional term: "
                "delta Phi_N = -beta * phi / M_Pl = "
                "+(beta^2 * M) / (4 pi M_Pl^2 r) * exp(-m_phi * r)"
            ),
            "result": (
                "The time-time metric perturbation h_00 = -2 Phi_N / c^2 "
                "acquires a Yukawa correction proportional to beta^2"
            ),
        },
        {
            "step_number": 4,
            "description": "Read off the effective Newton's constant",
            "formula": (
                "G_eff(r) = G_N * (1 + 2 * beta^2 * exp(-m_phi * r))"
            ),
            "result": (
                f"Comparing with the standard Newtonian potential "
                f"Phi_N = -G_eff * M / r, we identify the Yukawa correction.  "
                f"The factor of 2 arises because the dilaton couples to "
                f"both the source and the test mass (two vertices in the "
                f"tree-level exchange diagram), each contributing a factor "
                f"of beta.  With omega_BD = {omega_BD:.6f}: "
                f"2*beta^2 = {alpha_screening_lagrangian:.10f}"
            ),
        },
        {
            "step_number": 5,
            "description": "Identify the screening parameters",
            "formula": (
                "alpha_screening = 2 * beta^2 = "
                f"2 / (2*omega + 3) = {alpha_screening_lagrangian:.10f}; "
                "lambda_dilaton = 1/m_phi = hbar / (m_phi * c)"
            ),
            "result": (
                f"The Yukawa profile G_eff = G_vacuum * "
                f"(1 + {alpha_screening_lagrangian:.10f} * exp(-r/lambda)) "
                f"is derived from the dilaton Lagrangian, not assumed.  "
                f"G_vacuum = G_N is the bare Newton's constant from the "
                f"alpha ladder prediction."
            ),
        },
    ]

    yukawa_formula = (
        f"G_eff(r) = G_vacuum * (1 + 2*beta^2 * exp(-r/lambda)) "
        f"= G_vacuum * (1 + {alpha_screening_lagrangian:.10f} * exp(-r/lambda))"
    )

    consistency_check = {
        "omega_BD": omega_BD,
        "beta_coupling": beta_coupling,
        "alpha_from_lagrangian": alpha_screening_lagrangian,
        "formula_matches_phenomenological": True,
        "note": (
            "The Yukawa profile derived from the Lagrangian has the same "
            "functional form as the phenomenological model in "
            "compute_screening_parameters.  The amplitude alpha_screening "
            "is predicted by the coupling beta, while the range lambda "
            "is set by the dilaton mass m_phi."
        ),
    }

    return {
        "derivation_steps": derivation_steps,
        "beta_coupling": beta_coupling,
        "alpha_screening_from_lagrangian": alpha_screening_lagrangian,
        "yukawa_formula": yukawa_formula,
        "consistency_check": consistency_check,
    }


def verify_screening_consistency(constants):
    """Cross-check the Lagrangian derivation against the phenomenological model.

    Computes the dilaton coupling beta from the Brans-Dicke parameter omega_BD
    (obtained from the dilaton module), derives the predicted screening amplitude
    alpha_screening_lagrangian = 2 * beta^2, and compares it with the empirical
    value alpha_screening_empirical = (G_lab - G_vacuum) / G_vacuum.

    Parameters
    ----------
    constants : SimpleNamespace
        Physical and mathematical constants from get_constants().

    Returns
    -------
    dict with keys:
        omega_BD : float
        beta_coupling : float
        alpha_screening_lagrangian : float
        alpha_screening_empirical : float
        ratio : float
        consistent : bool
        interpretation : str
    """
    from alpha_ladder_core.dilaton import compute_bd_parameter

    # Get omega from the dilaton module
    bd = compute_bd_parameter(constants)
    omega_BD = bd["omega"]

    # Coupling constant
    denom = 2.0 * omega_BD + 3.0
    if denom <= 0:
        beta_coupling = float('inf')
    else:
        beta_coupling = 1.0 / math.sqrt(denom)

    # Lagrangian prediction
    alpha_screening_lagrangian = 2.0 * beta_coupling ** 2

    # Empirical value from the phenomenological model
    sp = compute_screening_parameters(constants)
    alpha_screening_empirical = sp["alpha_screening"]

    # Ratio and consistency
    if alpha_screening_empirical != 0:
        ratio = alpha_screening_lagrangian / alpha_screening_empirical
    else:
        ratio = float('inf')

    # Tolerance: within a factor of 10 is "broadly consistent" for this
    # kind of order-of-magnitude cross-check between a tree-level
    # Lagrangian prediction and an empirical extraction
    tolerance = 10.0
    consistent = (1.0 / tolerance) < ratio < tolerance

    # Physical interpretation
    if abs(ratio - 1.0) < 0.01:
        interpretation = (
            "The Lagrangian prediction and empirical extraction agree to "
            "better than 1%.  The dilaton coupling beta derived from the "
            "Brans-Dicke parameter omega_BD fully accounts for the observed "
            "excess of G_lab over G_vacuum."
        )
    elif consistent:
        interpretation = (
            f"The Lagrangian prediction (alpha = {alpha_screening_lagrangian:.6e}) "
            f"and empirical value (alpha = {alpha_screening_empirical:.6e}) "
            f"differ by a factor of {ratio:.4f}.  "
            f"This discrepancy may indicate: "
            f"(1) loop corrections to the tree-level dilaton exchange; "
            f"(2) additional scalar fields beyond the single dilaton; "
            f"(3) non-linear effects in the conformal coupling A(phi); or "
            f"(4) the need for a more precise determination of omega_BD "
            f"from the KK reduction (e.g. including Gauss-Bonnet corrections)."
        )
    else:
        interpretation = (
            f"Large discrepancy: Lagrangian alpha = {alpha_screening_lagrangian:.6e}, "
            f"empirical alpha = {alpha_screening_empirical:.6e} "
            f"(ratio = {ratio:.4e}).  "
            f"The tree-level single-dilaton exchange does not account for "
            f"the observed screening amplitude.  This suggests that the "
            f"effective omega_BD from the KK reduction (omega = {omega_BD:.6f}) "
            f"does not directly set the phenomenological screening strength.  "
            f"Possible resolutions: (a) the screening is dominated by a "
            f"different modulus field with a different coupling; (b) the "
            f"relationship between the KK breathing mode and the physical "
            f"dilaton involves mixing with other scalars; (c) non-perturbative "
            f"effects modify the coupling at the relevant energy scale."
        )

    return {
        "omega_BD": omega_BD,
        "beta_coupling": beta_coupling,
        "alpha_screening_lagrangian": alpha_screening_lagrangian,
        "alpha_screening_empirical": alpha_screening_empirical,
        "ratio": ratio,
        "consistent": consistent,
        "interpretation": interpretation,
    }
