"""
Universe Slider: explore what happens to physics when alpha takes
a different value.

Given a hypothetical fine-structure constant, recompute all derived
quantities and flag whether the resulting universe supports chemistry,
stellar fusion, and perturbative QED.
"""

from decimal import Decimal, getcontext

getcontext().prec = 50


def recompute_physics(alpha_val, constants):
    """Given a hypothetical alpha value, recompute all derived quantities.

    Args:
        alpha_val: Decimal or float - the hypothetical fine-structure constant.
        constants: SimpleNamespace from ``get_constants()``.

    Returns:
        dict with keys:
        - alpha: the input alpha (as Decimal)
        - inv_alpha: 1 / alpha
        - alpha_powers: dict of {n: alpha^n for n in 1..24}
        - alpha_g: predicted gravitational coupling (using phi^2/2 bridge)
        - G_predicted: predicted Newton's G (phi^2/2 bridge)
        - alpha_g_hierarchy: gravitational coupling from alpha^24 * mu^2
        - G_hierarchy: predicted Newton's G (zero free parameters)
        - r_e: classical electron radius (alpha * lambda_bar_c)
        - a_0_predicted: predicted Bohr radius (lambda_bar_c / alpha)
        - binding_energy_ratio: alpha^2 / 2 (hydrogen ground-state binding
          energy as a fraction of the electron rest-mass energy)
        - dark_coupling: alpha^10
        - gap: alpha / alpha_g
        - chemistry_viable: bool - atoms can form stable bonds
        - stars_viable: bool - stellar fusion is possible
        - notes: list of string observations about this universe
    """
    alpha = Decimal(str(alpha_val))

    # Inverse
    inv_alpha = Decimal(1) / alpha

    # Powers 1..24
    alpha_powers = {}
    power = Decimal(1)
    for n in range(1, 25):
        power *= alpha
        alpha_powers[n] = power

    # Bridge coefficient: phi^2 / 2
    phi = constants.phi
    bridge_coeff = phi ** 2 / 2

    # Gravitational coupling: alpha_g = bridge_coeff * alpha^21
    alpha_g = bridge_coeff * alpha_powers[21]

    # G_predicted = alpha_g * hbar * c / m_e^2
    hbar = constants.hbar
    c = constants.c
    m_e = constants.m_e
    G_predicted = alpha_g * hbar * c / m_e ** 2

    # --- Second prediction: alpha^24 * mu^2 (zero free parameters) ---
    m_p = constants.m_p
    mu = m_p / m_e
    alpha_g_hierarchy = alpha_powers[24] * mu ** 2
    G_hierarchy = alpha_g_hierarchy * hbar * c / m_e ** 2

    # Classical electron radius: r_e = alpha * lambda_bar_c
    lambda_bar_c = constants.lambda_bar_c
    r_e = alpha * lambda_bar_c

    # Bohr radius: a_0 = lambda_bar_c / alpha = hbar / (m_e * c * alpha)
    a_0_predicted = lambda_bar_c / alpha

    # Binding energy ratio: hydrogen ground state / rest mass ~ alpha^2 / 2
    binding_energy_ratio = alpha_powers[2] / 2

    # Dark coupling proxy
    dark_coupling = alpha_powers[10]

    # Hierarchy gap between EM and gravity
    gap = alpha / alpha_g

    # Viability flags
    alpha_upper_chem = Decimal(1) / Decimal(85)
    alpha_lower_chem = Decimal(1) / Decimal(200)
    alpha_upper_stars = Decimal(1) / Decimal(85)
    alpha_lower_stars = Decimal(1) / Decimal(180)

    chemistry_viable = alpha_lower_chem <= alpha <= alpha_upper_chem
    stars_viable = alpha_lower_stars <= alpha <= alpha_upper_stars

    # Observations
    notes = []

    if alpha > Decimal("0.1"):
        notes.append(
            "QED perturbation theory breaks down; higher-order corrections "
            "are no longer small."
        )

    if alpha > alpha_upper_chem:
        notes.append(
            "Atoms become unstable for Z > ~85; chemistry is severely limited."
        )

    if alpha < alpha_lower_chem:
        notes.append(
            "Molecular bonds are too weak for complex chemistry."
        )

    if not stars_viable:
        direction = "high" if alpha > alpha_upper_stars else "low"
        notes.append(
            f"Stellar fusion is not viable; alpha is too {direction} for "
            f"stable stellar structure."
        )

    real_alpha = constants.alpha
    if alpha != real_alpha:
        ratio = float(alpha / real_alpha)
        notes.append(
            f"Alpha is {ratio:.4f}x the real-universe value "
            f"(1/{float(inv_alpha):.2f} vs 1/{float(Decimal(1)/real_alpha):.2f})."
        )

    if chemistry_viable and stars_viable:
        notes.append(
            "This universe supports both complex chemistry and stellar fusion."
        )

    # Report hierarchy shift
    real_alpha_g = constants.alpha_g
    real_gap = float(real_alpha / real_alpha_g)
    new_gap = float(gap)
    notes.append(
        f"EM-gravity hierarchy gap: {new_gap:.4e} "
        f"(real universe: {real_gap:.4e})."
    )

    return {
        "alpha": alpha,
        "inv_alpha": inv_alpha,
        "alpha_powers": alpha_powers,
        "alpha_g": alpha_g,
        "G_predicted": G_predicted,
        "alpha_g_hierarchy": alpha_g_hierarchy,
        "G_hierarchy": G_hierarchy,
        "r_e": r_e,
        "a_0_predicted": a_0_predicted,
        "binding_energy_ratio": binding_energy_ratio,
        "dark_coupling": dark_coupling,
        "gap": gap,
        "chemistry_viable": chemistry_viable,
        "stars_viable": stars_viable,
        "notes": notes,
    }


def compute_sensitivity_curve(constants, n_points=100):
    """Compute physics across a range of alpha values from 1/200 to 1/85.

    Args:
        constants: SimpleNamespace from ``get_constants()``.
        n_points: number of sample points (default 100).

    Returns:
        list of dicts from ``recompute_physics``, one per alpha value,
        sorted from smallest to largest alpha.
    """
    alpha_min = Decimal(1) / Decimal(200)
    alpha_max = Decimal(1) / Decimal(85)
    step = (alpha_max - alpha_min) / (n_points - 1) if n_points > 1 else Decimal(0)

    results = []
    for i in range(n_points):
        alpha_val = alpha_min + step * i
        results.append(recompute_physics(alpha_val, constants))
    return results
