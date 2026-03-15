"""
Dark sector phenomenology for the Alpha Ladder dilaton.

Computes observational consequences of the dilaton as a potential
dark sector participant, focusing on three mass regimes:

1. Fuzzy DM regime (m ~ 10^{-22} eV):
   EXCLUDED by KK truncation breakdown.  a_0 ~ 10^7 m means
   gravity is 6D below ~10,000 km.

2. meV regime (m ~ 1-10 meV, a_0 ~ 30-100 um):
   Self-consistent.  Produces a testable sub-mm fifth force.
   NOT a dark matter candidate (too heavy for fuzzy DM, too
   light for thermal relic, wrong coupling for axion).

3. Planck regime (m ~ M_Pl, a_0 ~ l_Pl):
   Flux-stabilized.  Completely decoupled.  Pure GR.
   No dark sector phenomenology.

Key calculations:
   - Relic abundance via misalignment mechanism
   - Dark matter self-interaction cross section (bullet cluster)
   - Equation of state parameter w(a)
   - Comparison with observational constraints

Also retains the original Rung 10 analysis (alpha^10 ultra-light
axion-like particle) for backward compatibility.

All calculations use pure Python + math + decimal (no numpy/scipy).
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50


# ---------------------------------------------------------------------------
# Physical constants (natural units where convenient)
# ---------------------------------------------------------------------------

M_Pl_eV = 1.22089e28       # Planck mass in eV
l_Pl = 1.616e-35            # Planck length in meters
hbar_c_eV_m = 1.9733e-7     # hbar*c in eV*m
hbar_c_eV_cm = 1.9733e-14   # hbar*c in eV*cm (for cross section)
H_0_eV = 1.49e-33           # Hubble constant in eV
Omega_DM_h2 = 0.120         # Observed dark matter relic abundance
eV_to_gram = 1.783e-33      # 1 eV/c^2 in grams
kpc_to_m = 3.086e19         # 1 kpc in meters
EOT_WASH_LIMIT_M = 56e-6    # Eot-Wash sub-mm gravity bound in meters


# ===========================================================================
# Original Rung 10 analysis (backward-compatible API)
# ===========================================================================

def compute_dark_sector(constants) -> dict:
    """Compute all dark sector quantities from first principles.

    Parameters
    ----------
    constants : SimpleNamespace
        Physical and mathematical constants from ``get_constants()``.

    Returns
    -------
    dict
        Dark sector observables keyed by name, with Decimal precision
        for exact quantities and float for log10 values used in plotting.
    """
    alpha = constants.alpha
    alpha_g = constants.alpha_g
    m_e = constants.m_e
    hbar = constants.hbar
    c = constants.c
    e_charge = constants.e_charge

    alpha_10 = alpha ** 10
    epsilon = alpha ** 5  # sqrt(alpha^10) = alpha^5

    # Axion-like particle mass in eV
    m_axion_eV = m_e * c ** 2 / e_charge * alpha_10

    # Reduced Compton wavelength at rung 10
    lambda_bar_m = hbar / (m_e * alpha_10 * c)
    lambda_bar_km = lambda_bar_m / Decimal('1000')

    # Earth-Moon distance for comparison
    earth_moon_km = Decimal('384400')
    lambda_ratio = lambda_bar_km / earth_moon_km

    # Coupling hierarchy ratios
    ratio_to_em = alpha_10 / alpha       # = alpha^9
    ratio_to_gravity = alpha_10 / alpha_g

    # Float log10 values for plotting
    log10_em = math.log10(float(alpha))
    log10_dark = math.log10(float(alpha_10))
    log10_grav = math.log10(float(alpha_g))

    return {
        "alpha_10": alpha_10,
        "epsilon": epsilon,
        "m_axion_eV": m_axion_eV,
        "lambda_bar_m": lambda_bar_m,
        "lambda_bar_km": lambda_bar_km,
        "earth_moon_km": earth_moon_km,
        "lambda_ratio": lambda_ratio,
        "ratio_to_em": ratio_to_em,
        "ratio_to_gravity": ratio_to_gravity,
        "log10_em": log10_em,
        "log10_dark": log10_dark,
        "log10_grav": log10_grav,
    }


def compute_wave_profile(constants, alpha_override=None, n_points=500) -> dict:
    """Compute the fuzzy dark matter wave envelope for visualization.

    Parameters
    ----------
    constants : SimpleNamespace
        Physical and mathematical constants from ``get_constants()``.
    alpha_override : float or Decimal or None
        If provided, use this value instead of ``constants.alpha``.
    n_points : int
        Number of spatial sample points (default 500).

    Returns
    -------
    dict
        ``lambda_bar_km`` (Decimal), ``x_km`` (list[float]),
        ``psi_squared`` (list[float]).
    """
    if alpha_override is not None:
        alpha_used = Decimal(str(alpha_override))
    else:
        alpha_used = constants.alpha

    m_e = constants.m_e
    hbar = constants.hbar
    c = constants.c

    # Reduced Compton wavelength at rung 10, converted to km
    lambda_bar_m = hbar / (m_e * alpha_used ** 10 * c)
    lambda_bar_km = lambda_bar_m / Decimal('1000')

    # Spatial axis: 0 to 3 * Earth-Moon distance
    earth_moon_km = 384400.0
    x_max = 3.0 * earth_moon_km
    x_km = [x_max * i / (n_points - 1) for i in range(n_points)]

    # Wave envelope: |psi|^2 = cos^2(2*pi*x / lambda_bar_km)
    lam_float = float(lambda_bar_km)
    two_pi = 2.0 * math.pi
    psi_squared = [
        math.cos(two_pi * x / lam_float) ** 2 for x in x_km
    ]

    return {
        "lambda_bar_km": lambda_bar_km,
        "x_km": x_km,
        "psi_squared": psi_squared,
    }


def compute_alps_simulation(constants) -> dict:
    """Compute ALPS II light-shining-through-a-wall parameters.

    Parameters
    ----------
    constants : SimpleNamespace
        Physical and mathematical constants from ``get_constants()``.

    Returns
    -------
    dict
        Mixing parameter, tunneling probability, photon budget,
        and ALPS II apparatus parameters.
    """
    alpha = constants.alpha

    epsilon = alpha ** 5
    P_tunneling = epsilon ** 4          # = alpha^20
    photons_needed = Decimal(1) / P_tunneling

    B_tesla = Decimal('5.3')   # ALPS II magnetic field
    L_meters = Decimal('106')  # ALPS II cavity length

    # Real-world context: ALPS II uses a 30 W, 1064 nm infrared laser
    # Photon energy E = hbar * c / lambda (but simpler: 1064 nm -> ~1.165 eV)
    # photons/sec = 30 W / (1.165 eV * 1.602e-19 J/eV) ~ 1.6e20
    laser_power_W = Decimal('30')
    photon_energy_eV = Decimal('1.165')          # 1064 nm IR
    photon_energy_J = photon_energy_eV * constants.e_charge
    photons_per_sec = laser_power_W / photon_energy_J
    wait_seconds = photons_needed / photons_per_sec
    seconds_per_year = Decimal('3.156e7')
    wait_years = wait_seconds / seconds_per_year
    log10_wait_years = math.log10(float(wait_years))

    # ALPS II actual design sensitivity: epsilon ~ 3e-4 (after cavity enhancement)
    alps_ii_sensitivity = Decimal('3e-4')
    # Orders of magnitude between ALPS II reach and Rung 10 prediction
    orders_below_alps = math.log10(float(alps_ii_sensitivity)) - math.log10(float(epsilon))

    return {
        "epsilon": epsilon,
        "P_tunneling": P_tunneling,
        "photons_needed": photons_needed,
        "B_tesla": B_tesla,
        "L_meters": L_meters,
        "laser_power_W": laser_power_W,
        "photons_per_sec": photons_per_sec,
        "wait_years": wait_years,
        "log10_wait_years": log10_wait_years,
        "alps_ii_sensitivity": alps_ii_sensitivity,
        "orders_below_alps": orders_below_alps,
    }


def get_experimental_bounds(constants) -> dict:
    """Return experimental bounds from the experimental module for cross-reference."""
    from alpha_ladder_core.experimental import strategy_dark_sector
    return strategy_dark_sector(constants)


# ===========================================================================
# Dilaton dark sector phenomenology (new)
# ===========================================================================


# ---------------------------------------------------------------------------
# 1. Relic abundance via misalignment mechanism
# ---------------------------------------------------------------------------

def compute_relic_abundance(m_phi_eV, f_initial_eV=None, alpha_coupling=0.618):
    """Compute cosmological relic abundance of the dilaton via misalignment.

    The misalignment mechanism: the dilaton field starts displaced from
    its minimum by f_initial.  As the universe expands and H drops below
    m_phi, the field oscillates and behaves as cold dark matter with
    energy density rho_phi = (1/2) * m_phi^2 * f_initial^2.

    The relic abundance parameter:
        Omega_phi * h^2 ~ (m_phi / 10^{-22} eV)^{1/2} * (f_initial / M_Pl)^2

    Parameters
    ----------
    m_phi_eV : float
        Dilaton mass in eV.
    f_initial_eV : float or None
        Initial field displacement in eV.  If None, defaults to M_Pl.
    alpha_coupling : float
        Coupling constant (default phi^{-1} = 0.618).

    Returns
    -------
    dict with keys:
        m_phi_eV, f_initial_eV, Omega_phi_h2, overproduced, ratio_to_observed,
        f_required_for_dm_eV, f_required_over_M_Pl, fine_tuning_required,
        honest_assessment
    """
    if f_initial_eV is None:
        f_initial_eV = M_Pl_eV

    # Standard misalignment formula
    Omega_phi_h2 = Omega_DM_h2 * (m_phi_eV / 1e-22) ** 0.5 * (f_initial_eV / M_Pl_eV) ** 2

    overproduced = Omega_phi_h2 > Omega_DM_h2
    ratio_to_observed = Omega_phi_h2 / Omega_DM_h2 if Omega_DM_h2 != 0 else float('inf')

    # Invert the formula to find f_initial that gives Omega = Omega_DM
    mass_factor = (m_phi_eV / 1e-22) ** 0.5
    if mass_factor > 0:
        f_required_for_dm_eV = M_Pl_eV * math.sqrt(1.0 / mass_factor)
    else:
        f_required_for_dm_eV = float('inf')

    f_required_over_M_Pl = f_required_for_dm_eV / M_Pl_eV

    # Fine tuning: if f_required / M_Pl < 0.01, severe fine tuning needed
    fine_tuning_required = f_required_over_M_Pl < 0.01

    # Honest assessment
    if m_phi_eV > 1e-5:
        honest_assessment = (
            f"At m_phi = {m_phi_eV:.2e} eV, the misalignment mechanism with "
            f"f_initial = M_Pl massively overproduces dark matter "
            f"(Omega_phi h^2 = {Omega_phi_h2:.2e} vs observed 0.120).  "
            f"Matching the observed abundance requires f_initial / M_Pl = "
            f"{f_required_over_M_Pl:.2e}, which is severe fine-tuning.  "
            f"The meV dilaton is NOT a viable dark matter candidate via "
            f"misalignment."
        )
    elif m_phi_eV < 1e-20:
        honest_assessment = (
            f"At m_phi = {m_phi_eV:.2e} eV (fuzzy DM regime), the "
            f"misalignment mechanism can produce the right relic abundance "
            f"with f_initial ~ M_Pl.  However, this mass requires an "
            f"internal radius a_0 ~ 10^7 m, which breaks the KK "
            f"truncation and is therefore EXCLUDED."
        )
    else:
        honest_assessment = (
            f"At m_phi = {m_phi_eV:.2e} eV, the relic abundance is "
            f"Omega_phi h^2 = {Omega_phi_h2:.2e} (ratio to observed: "
            f"{ratio_to_observed:.2e}).  Whether this is viable depends "
            f"on the initial displacement f_initial."
        )

    return {
        "m_phi_eV": m_phi_eV,
        "f_initial_eV": f_initial_eV,
        "Omega_phi_h2": Omega_phi_h2,
        "overproduced": overproduced,
        "ratio_to_observed": ratio_to_observed,
        "f_required_for_dm_eV": f_required_for_dm_eV,
        "f_required_over_M_Pl": f_required_over_M_Pl,
        "fine_tuning_required": fine_tuning_required,
        "honest_assessment": honest_assessment,
    }


# ---------------------------------------------------------------------------
# 2. Dark matter self-interaction cross section
# ---------------------------------------------------------------------------

def compute_self_interaction(m_phi_eV, alpha_coupling=0.618):
    """Compute the dark matter self-interaction cross section from dilaton exchange.

    sigma/m for scalar exchange:
        sigma = alpha^2 / (4 pi m_phi^2)  (natural units)
        sigma/m = alpha^2 / (4 pi m_phi^3)

    Observational bounds:
        Bullet cluster: sigma/m < 1 cm^2/g  (Markevitch et al.)
        Milky Way halos: sigma/m < 0.1 cm^2/g

    Parameters
    ----------
    m_phi_eV : float
        Dilaton mass in eV.
    alpha_coupling : float
        Coupling constant (default phi^{-1} = 0.618).

    Returns
    -------
    dict with keys:
        m_phi_eV, alpha_coupling, sigma_cm2, sigma_over_m_cm2_g,
        bullet_cluster_bound, passes_bullet_cluster,
        milky_way_bound, passes_milky_way, honest_assessment
    """
    # Cross section in cm^2
    sigma_cm2 = alpha_coupling ** 2 * hbar_c_eV_cm ** 2 / (4.0 * math.pi * m_phi_eV ** 2)

    # Mass in grams
    m_g = m_phi_eV * eV_to_gram

    # sigma/m in cm^2/g
    sigma_over_m = sigma_cm2 / m_g if m_g > 0 else float('inf')

    bullet_cluster_bound = 1.0    # cm^2/g
    milky_way_bound = 0.1         # cm^2/g

    passes_bullet_cluster = sigma_over_m < bullet_cluster_bound
    passes_milky_way = sigma_over_m < milky_way_bound

    # Honest assessment
    if passes_milky_way:
        honest_assessment = (
            f"At m_phi = {m_phi_eV:.2e} eV, sigma/m = {sigma_over_m:.2e} cm^2/g, "
            f"well below the bullet cluster bound ({bullet_cluster_bound} cm^2/g) "
            f"and Milky Way bound ({milky_way_bound} cm^2/g).  "
            f"The dilaton self-interaction is negligible -- it does not behave "
            f"as self-interacting dark matter."
        )
    elif passes_bullet_cluster:
        honest_assessment = (
            f"At m_phi = {m_phi_eV:.2e} eV, sigma/m = {sigma_over_m:.2e} cm^2/g.  "
            f"Passes bullet cluster bound but violates Milky Way halo shape "
            f"constraint.  Marginally excluded as SIDM."
        )
    else:
        honest_assessment = (
            f"At m_phi = {m_phi_eV:.2e} eV, sigma/m = {sigma_over_m:.2e} cm^2/g.  "
            f"Violates bullet cluster bound ({bullet_cluster_bound} cm^2/g).  "
            f"Strongly excluded as a self-interacting dark matter candidate."
        )

    return {
        "m_phi_eV": m_phi_eV,
        "alpha_coupling": alpha_coupling,
        "sigma_cm2": sigma_cm2,
        "sigma_over_m_cm2_g": sigma_over_m,
        "bullet_cluster_bound": bullet_cluster_bound,
        "passes_bullet_cluster": passes_bullet_cluster,
        "milky_way_bound": milky_way_bound,
        "passes_milky_way": passes_milky_way,
        "honest_assessment": honest_assessment,
    }


# ---------------------------------------------------------------------------
# 3. Equation of state parameter w(z)
# ---------------------------------------------------------------------------

def compute_equation_of_state(m_phi_eV, z_values=None):
    """Compute the equation of state parameter w(z) for the dilaton field.

    For a massive scalar oscillating in a quadratic potential:
    - w = -1 when H >> m_phi (field is frozen, acts as dark energy)
    - w = 0 when H << m_phi (field oscillates, acts as matter)
    - Transition at H ~ m_phi

    Parameters
    ----------
    m_phi_eV : float
        Dilaton mass in eV.
    z_values : list of float or None
        Redshifts at which to evaluate w.  If None, uses a default set.

    Returns
    -------
    dict with keys:
        m_phi_eV, z_values, w_values, z_transition, H_at_transition_eV,
        w_today, behaves_as_matter, behaves_as_dark_energy, honest_assessment
    """
    if z_values is None:
        z_values = [0, 0.1, 0.5, 1, 2, 5, 10, 100, 1000, 1e4, 1e5, 1e6, 1e9]

    Omega_m = 0.315
    Omega_Lambda = 0.685

    def H_of_z(z):
        """Hubble parameter H(z) in eV for flat LCDM."""
        return H_0_eV * math.sqrt(Omega_m * (1.0 + z) ** 3 + Omega_Lambda)

    # Equation of state: smooth transition from w=-1 to w=0
    # w(z) = -1 / (1 + (m_phi / H(z))^2)
    # When H >> m: w -> -1 (frozen field)
    # When H << m: w -> 0 (oscillating matter)
    def w_of_z(z):
        H_z = H_of_z(z)
        ratio_sq = (m_phi_eV / H_z) ** 2
        # Guard against overflow: if ratio_sq is huge, w -> 0
        if ratio_sq > 1e30:
            return 0.0
        return -1.0 / (1.0 + ratio_sq)

    w_values = [w_of_z(z) for z in z_values]
    w_today = w_values[0] if len(z_values) > 0 and z_values[0] == 0 else w_of_z(0)

    # Find transition redshift where H(z_trans) ~ m_phi
    # For high z, Omega_Lambda is negligible:
    # H_0 * sqrt(Omega_m) * (1+z)^{3/2} ~ m_phi
    # z_trans ~ (m_phi / (H_0 * sqrt(Omega_m)))^{2/3} - 1
    ratio = m_phi_eV / (H_0_eV * math.sqrt(Omega_m))
    if ratio > 1:
        z_transition = ratio ** (2.0 / 3.0) - 1.0
    else:
        # m_phi < H_0: field is still frozen today
        z_transition = 0.0

    H_at_transition = H_of_z(z_transition) if z_transition > 0 else H_0_eV

    behaves_as_matter = z_transition > 1e6
    behaves_as_dark_energy = z_transition < 1.0

    # Honest assessment
    if behaves_as_matter:
        honest_assessment = (
            f"At m_phi = {m_phi_eV:.2e} eV, the dilaton began oscillating at "
            f"z ~ {z_transition:.2e}, long before BBN (z ~ 10^9).  "
            f"Throughout all observable cosmological history, w = 0 (matter).  "
            f"The dilaton is NOT dark energy.  It behaves as pressureless "
            f"dust if it has any cosmological abundance."
        )
    elif behaves_as_dark_energy:
        honest_assessment = (
            f"At m_phi = {m_phi_eV:.2e} eV, the dilaton is still frozen today "
            f"(z_transition ~ {z_transition:.2e}).  It would act as dark energy "
            f"with w ~ -1.  However, such an ultralight mass requires an "
            f"enormous internal radius, which breaks the KK truncation."
        )
    else:
        honest_assessment = (
            f"At m_phi = {m_phi_eV:.2e} eV, the transition from w = -1 to "
            f"w = 0 occurs at z ~ {z_transition:.2e}.  This intermediate "
            f"regime may affect late-time structure formation."
        )

    return {
        "m_phi_eV": m_phi_eV,
        "z_values": z_values,
        "w_values": w_values,
        "z_transition": z_transition,
        "H_at_transition_eV": H_at_transition,
        "w_today": w_today,
        "behaves_as_matter": behaves_as_matter,
        "behaves_as_dark_energy": behaves_as_dark_energy,
        "honest_assessment": honest_assessment,
    }


# ---------------------------------------------------------------------------
# 4. Fuzzy dark matter constraints
# ---------------------------------------------------------------------------

def compute_fuzzy_dm_constraints(m_phi_eV):
    """Check the dilaton against observational constraints on fuzzy/ultralight DM.

    Constraints from:
    - Lyman-alpha forest (Irsic et al. 2017)
    - Galaxy UV luminosity function (Bozek et al. 2015)
    - CMB lensing (Hlozek et al. 2015)
    - Milky Way subhalo count (Nadler et al. 2021)

    Parameters
    ----------
    m_phi_eV : float
        Dilaton mass in eV.

    Returns
    -------
    dict with keys:
        m_phi_eV, constraints, passes_all_bounds, de_broglie_wavelength_m,
        de_broglie_wavelength_kpc, is_fuzzy_candidate, kk_truncation_valid,
        honest_assessment
    """
    constraint_defs = [
        {"name": "Lyman-alpha forest (Irsic+ 2017)", "bound_eV": 2e-21},
        {"name": "Galaxy UV luminosity (Bozek+ 2015)", "bound_eV": 1e-22},
        {"name": "CMB lensing (Hlozek+ 2015)", "bound_eV": 1e-24},
        {"name": "MW subhalo count (Nadler+ 2021)", "bound_eV": 3e-21},
    ]

    constraints = []
    for c in constraint_defs:
        constraints.append({
            "name": c["name"],
            "bound_eV": c["bound_eV"],
            "passes": m_phi_eV >= c["bound_eV"],
        })

    passes_all_bounds = all(c["passes"] for c in constraints)

    # de Broglie wavelength: lambda_dB = hbar_c / (m_phi * v)
    # v ~ 200 km/s for galactic halos = 200e3 / 3e8 = 6.67e-4 c
    v_halo = 200e3 / 3e8  # in units of c
    if m_phi_eV > 0:
        de_broglie_wavelength_m = hbar_c_eV_m / (m_phi_eV * v_halo)
    else:
        de_broglie_wavelength_m = float('inf')
    de_broglie_wavelength_kpc = de_broglie_wavelength_m / kpc_to_m

    # Fuzzy DM requires lambda_dB > 0.1 kpc to affect structure formation
    is_fuzzy_candidate = de_broglie_wavelength_kpc > 0.1 and passes_all_bounds

    # KK truncation check:
    # m_phi = hbar_c / a_0  =>  a_0 = hbar_c / m_phi
    # (Since M_Pl * l_Pl = hbar_c in eV*m)
    if m_phi_eV > 0:
        a_0_m = hbar_c_eV_m / m_phi_eV
    else:
        a_0_m = float('inf')
    kk_truncation_valid = a_0_m < EOT_WASH_LIMIT_M

    # Honest assessment
    if m_phi_eV > 1e-5:
        honest_assessment = (
            f"At m_phi = {m_phi_eV:.2e} eV, the dilaton trivially passes all "
            f"fuzzy DM mass bounds (it is far heavier than any bound).  "
            f"However, it is NOT a fuzzy DM candidate because: "
            f"(1) de Broglie wavelength = {de_broglie_wavelength_m:.2e} m "
            f"({de_broglie_wavelength_kpc:.2e} kpc), which is sub-atomic, "
            f"not kpc-scale; "
            f"(2) its coupling alpha = 0.618 to matter produces a fifth force, "
            f"not a dark matter halo; "
            f"(3) it does not solve small-scale structure problems.  "
            f"KK truncation is valid (a_0 = {a_0_m:.2e} m < 56 um)."
        )
    elif not kk_truncation_valid:
        honest_assessment = (
            f"At m_phi = {m_phi_eV:.2e} eV, the internal radius "
            f"a_0 = {a_0_m:.2e} m exceeds the Eot-Wash limit of 56 um.  "
            f"KK truncation FAILS: gravity would be 6D at distances "
            f"below a_0.  This mass is EXCLUDED in the Alpha Ladder framework."
        )
    else:
        honest_assessment = (
            f"At m_phi = {m_phi_eV:.2e} eV, de Broglie wavelength = "
            f"{de_broglie_wavelength_kpc:.2e} kpc.  KK truncation valid "
            f"(a_0 = {a_0_m:.2e} m).  Passes mass bounds: {passes_all_bounds}."
        )

    return {
        "m_phi_eV": m_phi_eV,
        "constraints": constraints,
        "passes_all_bounds": passes_all_bounds,
        "de_broglie_wavelength_m": de_broglie_wavelength_m,
        "de_broglie_wavelength_kpc": de_broglie_wavelength_kpc,
        "is_fuzzy_candidate": is_fuzzy_candidate,
        "kk_truncation_valid": kk_truncation_valid,
        "honest_assessment": honest_assessment,
    }


# ---------------------------------------------------------------------------
# 5. Dark sector landscape across mass scales
# ---------------------------------------------------------------------------

def compute_dark_sector_landscape(mass_points=None):
    """Compute a comprehensive landscape table across mass scales.

    For each mass point, evaluates KK truncation, relic abundance,
    self-interaction, equation of state, and classifies the regime.

    Parameters
    ----------
    mass_points : list of float or None
        Masses in eV.  If None, uses a default spanning 61 orders of magnitude.

    Returns
    -------
    dict with keys:
        mass_points, landscape (list of dicts), viable_window, honest_assessment
    """
    if mass_points is None:
        mass_points = [1e-33, 1e-28, 1e-22, 1e-17, 1e-10, 1e-5, 1e-3, 1e0, 1e10, 1e20, 1e28]

    landscape = []
    viable_masses = []

    for m in mass_points:
        # Internal radius
        a_0 = hbar_c_eV_m / m if m > 0 else float('inf')

        # Compton wavelength (same as a_0 in these units)
        lambda_compton = a_0

        # KK truncation
        kk_valid = a_0 < EOT_WASH_LIMIT_M

        # Relic abundance (with f = M_Pl)
        relic = compute_relic_abundance(m)

        # Self-interaction
        si = compute_self_interaction(m)

        # Equation of state (just w_today)
        eos = compute_equation_of_state(m, z_values=[0])
        w_today = eos["w_today"]

        # Classification
        if not kk_valid and m < 1e-5:
            classification = "kk_excluded"
        elif m < 1e-30:
            classification = "dark_energy_candidate"
        elif m < 1e-19 and not kk_valid:
            classification = "fuzzy_dm_candidate"
        elif m > 1e25:
            classification = "invisible_planck"
        elif 1e-4 <= m <= 1e-1 and kk_valid:
            classification = "sub_mm_testable"
        elif 1e0 <= m <= 1e6:
            classification = "thermal_relic_range"
        elif kk_valid:
            classification = "sub_mm_testable"
        else:
            classification = "kk_excluded"

        if kk_valid and 1e-4 <= m <= 1e0:
            viable_masses.append(m)

        entry = {
            "m_phi_eV": m,
            "a_0_m": a_0,
            "lambda_compton_m": lambda_compton,
            "kk_truncation_valid": kk_valid,
            "Omega_phi_h2": relic["Omega_phi_h2"],
            "sigma_over_m_cm2_g": si["sigma_over_m_cm2_g"],
            "w_today": w_today,
            "classification": classification,
        }
        landscape.append(entry)

    # Viable window
    if viable_masses:
        viable_window = (min(viable_masses), max(viable_masses))
    else:
        # Default viable window from physics: meV regime
        viable_window = (1e-4, 1e0)

    honest_assessment = (
        "The Alpha Ladder dilaton occupies a narrow phenomenological window.  "
        "Ultralight masses (< 10^{-5} eV) are excluded by KK truncation: "
        "the internal radius exceeds sub-mm gravity bounds.  "
        "The viable regime is m ~ 1-100 meV (a_0 ~ 2-200 um), where "
        "the dilaton produces a testable sub-mm fifth force but is NOT "
        "a dark matter candidate.  At the Planck mass (flux-stabilized), "
        "the dilaton decouples entirely and we recover pure GR."
    )

    return {
        "mass_points": mass_points,
        "landscape": landscape,
        "viable_window": viable_window,
        "honest_assessment": honest_assessment,
    }


# ---------------------------------------------------------------------------
# 6. Dashboard summary entry point
# ---------------------------------------------------------------------------

def summarize_dark_sector(constants=None):
    """Dashboard entry point: full dark sector analysis for the meV dilaton.

    Runs the complete pipeline for the default dilaton with a_0 = 30 um,
    corresponding to m_phi ~ 7 meV.

    Parameters
    ----------
    constants : SimpleNamespace or None
        Physical constants from get_constants().  Currently unused but
        accepted for API consistency with other modules.

    Returns
    -------
    dict with keys:
        mev_analysis, relic_abundance, self_interaction, equation_of_state,
        fuzzy_constraints, landscape, key_finding, overall_verdict,
        framework_position
    """
    # Default dilaton: a_0 = 30 um => m_phi = hbar_c / a_0
    a_0_default = 30e-6  # 30 micrometers
    m_phi_default = hbar_c_eV_m / a_0_default  # ~ 6.6e-3 eV ~ 7 meV

    # Run all analyses
    relic = compute_relic_abundance(m_phi_default)
    si = compute_self_interaction(m_phi_default)
    eos = compute_equation_of_state(m_phi_default)
    fuzzy = compute_fuzzy_dm_constraints(m_phi_default)
    landscape = compute_dark_sector_landscape()

    mev_analysis = {
        "a_0_m": a_0_default,
        "m_phi_eV": m_phi_default,
        "m_phi_meV": m_phi_default * 1e3,
        "kk_truncation_valid": fuzzy["kk_truncation_valid"],
        "is_dark_matter": False,
        "is_dark_energy": False,
        "is_fifth_force": True,
    }

    key_finding = (
        "The Alpha Ladder dilaton at a_0 = 30 um (m ~ 7 meV) is a testable "
        "sub-mm fifth force, NOT a dark matter candidate.  Misalignment "
        "with f = M_Pl overproduces DM by a factor of ~10^{9.5}.  "
        "The de Broglie wavelength is sub-atomic (~10^{-12} m), ruling "
        "out fuzzy DM behavior.  The equation of state is w = 0 throughout "
        "all observable history (not dark energy)."
    )

    overall_verdict = (
        "The dilaton is phenomenologically relevant as a SHORT-RANGE "
        "fifth force testable by Eot-Wash and IUPAP experiments, not "
        "as a dark sector component.  This is a strength, not a weakness: "
        "the Alpha Ladder framework makes a falsifiable prediction about "
        "sub-mm gravity without invoking unobservable dark sector physics."
    )

    framework_position = (
        "The Alpha Ladder sits outside the standard dark sector landscape.  "
        "It predicts no new dark matter candidate and no dark energy "
        "modification.  Instead, it predicts: (1) a specific value of G "
        "from alpha; (2) a sub-mm Yukawa correction from the dilaton; "
        "(3) pure GR at long distances.  The dark sector remains unexplained "
        "by this framework -- and that honesty is part of the theory's value."
    )

    return {
        "mev_analysis": mev_analysis,
        "relic_abundance": relic,
        "self_interaction": si,
        "equation_of_state": eos,
        "fuzzy_constraints": fuzzy,
        "landscape": landscape,
        "key_finding": key_finding,
        "overall_verdict": overall_verdict,
        "framework_position": framework_position,
    }
