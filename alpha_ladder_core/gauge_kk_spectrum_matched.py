"""
Kaluza-Klein Spectrum on S^2 with Dilaton Vev Matched to alpha_EM
=================================================================

Computes the full KK spectrum for scalar, vector, graviton, and fermion
fields on an S^2 internal space, with the dilaton vev determined by
matching the KK gauge coupling to alpha_EM:

    phi_vev = (1/4) ln(4 pi alpha_EM)  ~  -0.597

The physical radius of the internal sphere is R_phys = a_0 * e^{phi_vev},
making all KK modes 1/e^{phi_vev} ~ 1.817x heavier than the naive
unmatched case.

KK mass formulas on S^2:
    Scalar:    m_l^2 = l(l+1)/R^2,         l = 0, 1, 2, ...   deg = 2l+1
    Vector:    m_l^2 = l(l+1)/R^2,         l = 1, 2, 3, ...   deg = 2(2l+1)
    Graviton:  m_l^2 = [l(l+1)-2]/R^2,     l = 2, 3, 4, ...   deg = 2l+1
    Fermion:   m_l^2 = (l+1/2)^2/R^2,      l = 0, 1, 2, ...   deg = 2(2l+1)

Pure Python -- only ``import math``, no numpy/scipy.

Functions
---------
1. scalar_spectrum           -- scalar KK masses on S^2
2. vector_spectrum           -- vector (gauge) KK masses on S^2
3. graviton_spectrum         -- graviton KK masses on S^2
4. fermion_spectrum          -- Dirac fermion KK masses on S^2
5. full_spectrum             -- all spins combined, sorted by mass
6. spectrum_at_notable_a0    -- lightest modes at notable a_0 values
7. matched_vs_unmatched      -- side-by-side comparison
8. cumulative_modes_vs_energy-- N(<E) vs energy, comparison with ADD
9. summarize_kk_spectrum     -- full report
"""

import math

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
ALPHA_EM = 1.0 / 137.035999084            # CODATA 2018
M_PL_EV = 1.22089e28                       # Planck mass (eV)
L_PL = 1.61625e-35                         # Planck length (m)
HBAR_C_EVM = 1.9733e-7                     # hbar*c in eV m

# Derived
PHI_VEV = 0.25 * math.log(4.0 * math.pi * ALPHA_EM)   # dilaton vev ~ -0.597
E_PHI_VEV = math.exp(PHI_VEV)                          # e^{phi_vev} ~ 0.5503
MASS_ENHANCEMENT = 1.0 / E_PHI_VEV                     # ~ 1.817
EV_PER_TEV = 1.0e12


# ===================================================================
# Helpers
# ===================================================================

def _r_phys(a_0_m, matched):
    """Physical radius of S^2 in meters."""
    if matched:
        return a_0_m * E_PHI_VEV
    return a_0_m


def _m6_from_a0(a_0_m):
    """M_6 in eV from a_0 in meters.  M_Pl^2 = 4 pi a_0^2 M_6^4."""
    a_0_nat = a_0_m / HBAR_C_EVM
    m6_fourth = M_PL_EV ** 2 / (4.0 * math.pi * a_0_nat ** 2)
    return m6_fourth ** 0.25


def _mass_eV(eigenvalue, r_phys_m):
    """Convert dimensionless eigenvalue to mass in eV.

    m = eigenvalue * hbar_c / R_phys_m
    """
    return eigenvalue * HBAR_C_EVM / r_phys_m


# ===================================================================
# 1. scalar_spectrum
# ===================================================================

def scalar_spectrum(a_0_m, l_max=50, matched=True):
    """Compute scalar KK masses on S^2.

    m_l = sqrt(l(l+1)) / R_phys,  l = 0, 1, 2, ...
    Degeneracy = 2l + 1.
    l = 0 is the zero mode (massless 4D dilaton).

    Returns dict with keys:
        modes       -- list of dicts {l, m_eV, m_TeV, degeneracy, cumulative_modes}
        R_phys_m    -- physical radius in meters
        phi_vev     -- dilaton vev used (0 if unmatched)
        a_0_m       -- bare radius
    """
    r = _r_phys(a_0_m, matched)
    vev = PHI_VEV if matched else 0.0
    modes = []
    cumulative = 0
    for l in range(0, l_max + 1):
        eigenvalue = math.sqrt(l * (l + 1))
        m_eV = _mass_eV(eigenvalue, r)
        deg = 2 * l + 1
        cumulative += deg
        modes.append({
            "l": l,
            "m_eV": m_eV,
            "m_TeV": m_eV / EV_PER_TEV,
            "degeneracy": deg,
            "cumulative_modes": cumulative,
        })
    return {
        "modes": modes,
        "R_phys_m": r,
        "phi_vev": vev,
        "a_0_m": a_0_m,
    }


# ===================================================================
# 2. vector_spectrum
# ===================================================================

def vector_spectrum(a_0_m, l_max=50, matched=True):
    """Compute vector (gauge) KK masses on S^2.

    From SO(3)/SO(2) coset: 2 physical gauge fields, each with
    m_l = sqrt(l(l+1)) / R_phys,  l = 1, 2, 3, ...
    Degeneracy per field = 2l + 1.  Total degeneracy = 2(2l+1).

    No l=0 mode (the massless gauge boson is the zero mode of one
    linear combination; the other has no zero mode due to SO(2) stabilizer).

    Returns dict with keys:
        modes       -- list of dicts
        R_phys_m    -- physical radius
        phi_vev     -- dilaton vev used
        a_0_m       -- bare radius
    """
    r = _r_phys(a_0_m, matched)
    vev = PHI_VEV if matched else 0.0
    modes = []
    cumulative = 0
    for l in range(1, l_max + 1):
        eigenvalue = math.sqrt(l * (l + 1))
        m_eV = _mass_eV(eigenvalue, r)
        deg_per_field = 2 * l + 1
        total_deg = 2 * deg_per_field
        cumulative += total_deg
        modes.append({
            "l": l,
            "m_eV": m_eV,
            "m_TeV": m_eV / EV_PER_TEV,
            "degeneracy_per_field": deg_per_field,
            "total_degeneracy": total_deg,
            "cumulative_modes": cumulative,
        })
    return {
        "modes": modes,
        "R_phys_m": r,
        "phi_vev": vev,
        "a_0_m": a_0_m,
    }


# ===================================================================
# 3. graviton_spectrum
# ===================================================================

def graviton_spectrum(a_0_m, l_max=50, matched=True):
    """Compute graviton KK masses on S^2.

    m_l^2 = [l(l+1) - 2] / R_phys^2,  l = 2, 3, 4, ...
    (l=0: zero mode graviton.  l=1: diffeomorphisms, removed.)
    For l=2: m_2 = sqrt(4)/R = 2/R.
    Degeneracy = 2l + 1.

    Returns dict with keys:
        modes       -- list of dicts
        R_phys_m    -- physical radius
        phi_vev     -- dilaton vev used
        a_0_m       -- bare radius
    """
    r = _r_phys(a_0_m, matched)
    vev = PHI_VEV if matched else 0.0
    modes = []
    cumulative = 0
    for l in range(2, l_max + 1):
        eigenvalue = math.sqrt(l * (l + 1) - 2)
        m_eV = _mass_eV(eigenvalue, r)
        deg = 2 * l + 1
        cumulative += deg
        modes.append({
            "l": l,
            "m_eV": m_eV,
            "m_TeV": m_eV / EV_PER_TEV,
            "degeneracy": deg,
            "cumulative_modes": cumulative,
        })
    return {
        "modes": modes,
        "R_phys_m": r,
        "phi_vev": vev,
        "a_0_m": a_0_m,
    }


# ===================================================================
# 4. fermion_spectrum
# ===================================================================

def fermion_spectrum(a_0_m, l_max=50, matched=True):
    """Compute Dirac fermion KK masses on S^2.

    m_l = (l + 1/2) / R_phys,  l = 0, 1, 2, ...
    Degeneracy = 2(2l+1) (spinor on S^2 has twice the scalar degeneracy).

    Note: fermions are NOT in the minimal pure-gravity framework.
    This is for the extended framework with charged matter (Paper 3).

    Returns dict with keys:
        modes       -- list of dicts
        R_phys_m    -- physical radius
        phi_vev     -- dilaton vev used
        a_0_m       -- bare radius
    """
    r = _r_phys(a_0_m, matched)
    vev = PHI_VEV if matched else 0.0
    modes = []
    cumulative = 0
    for l in range(0, l_max + 1):
        eigenvalue = l + 0.5
        m_eV = _mass_eV(eigenvalue, r)
        deg = 2 * (2 * l + 1)
        cumulative += deg
        modes.append({
            "l": l,
            "m_eV": m_eV,
            "m_TeV": m_eV / EV_PER_TEV,
            "degeneracy": deg,
            "cumulative_modes": cumulative,
        })
    return {
        "modes": modes,
        "R_phys_m": r,
        "phi_vev": vev,
        "a_0_m": a_0_m,
    }


# ===================================================================
# 5. full_spectrum
# ===================================================================

def full_spectrum(a_0_m, l_max=30, matched=True):
    """Combine all field types into one sorted list of KK modes.

    Includes scalar (l>=1), vector (l>=1), graviton (l>=2).
    Fermions excluded from minimal framework but included for
    completeness with a flag.

    Returns dict with keys:
        modes               -- sorted list of {spin, l, m_eV, m_TeV, degeneracy}
        count_below_1TeV    -- total modes with m < 1 TeV
        count_below_10TeV   -- total modes with m < 10 TeV
        count_below_100TeV  -- total modes with m < 100 TeV
        lightest_massive    -- the lightest massive mode
        mass_gap_eV         -- mass of lightest massive mode
        density_at_1TeV     -- approximate dN/dE at 1 TeV (modes per eV)
        R_phys_m            -- physical radius
    """
    r = _r_phys(a_0_m, matched)
    all_modes = []

    # Scalar: l >= 1 (l=0 is zero mode, exclude from massive spectrum)
    for l in range(1, l_max + 1):
        eigenvalue = math.sqrt(l * (l + 1))
        m_eV = _mass_eV(eigenvalue, r)
        deg = 2 * l + 1
        all_modes.append({
            "spin": "scalar",
            "l": l,
            "m_eV": m_eV,
            "m_TeV": m_eV / EV_PER_TEV,
            "degeneracy": deg,
        })

    # Vector: l >= 1
    for l in range(1, l_max + 1):
        eigenvalue = math.sqrt(l * (l + 1))
        m_eV = _mass_eV(eigenvalue, r)
        total_deg = 2 * (2 * l + 1)
        all_modes.append({
            "spin": "vector",
            "l": l,
            "m_eV": m_eV,
            "m_TeV": m_eV / EV_PER_TEV,
            "degeneracy": total_deg,
        })

    # Graviton: l >= 2
    for l in range(2, l_max + 1):
        eigenvalue = math.sqrt(l * (l + 1) - 2)
        m_eV = _mass_eV(eigenvalue, r)
        deg = 2 * l + 1
        all_modes.append({
            "spin": "graviton",
            "l": l,
            "m_eV": m_eV,
            "m_TeV": m_eV / EV_PER_TEV,
            "degeneracy": deg,
        })

    # Sort by mass
    all_modes.sort(key=lambda m: m["m_eV"])

    # Counts at thresholds
    count_1 = sum(m["degeneracy"] for m in all_modes if m["m_eV"] < 1.0e12)
    count_10 = sum(m["degeneracy"] for m in all_modes if m["m_eV"] < 1.0e13)
    count_100 = sum(m["degeneracy"] for m in all_modes if m["m_eV"] < 1.0e14)

    # Lightest massive mode
    lightest = all_modes[0] if all_modes else None

    # Mass gap
    mass_gap = lightest["m_eV"] if lightest else 0.0

    # Density of states at 1 TeV: approximate as Delta_N / Delta_E
    # around E = 1 TeV (within +/- 10%)
    e_center = 1.0e12  # 1 TeV in eV
    e_lo = 0.9e12
    e_hi = 1.1e12
    n_in_bin = sum(m["degeneracy"] for m in all_modes
                   if e_lo <= m["m_eV"] <= e_hi)
    delta_e = e_hi - e_lo
    density_1TeV = n_in_bin / delta_e if delta_e > 0 else 0.0

    return {
        "modes": all_modes,
        "count_below_1TeV": count_1,
        "count_below_10TeV": count_10,
        "count_below_100TeV": count_100,
        "lightest_massive": lightest,
        "mass_gap_eV": mass_gap,
        "density_at_1TeV": density_1TeV,
        "R_phys_m": r,
    }


# ===================================================================
# 6. spectrum_at_notable_a0
# ===================================================================

def spectrum_at_notable_a0(matched=True):
    """Compute the lightest 5 KK modes for each notable a_0 value.

    Notable values:
        l_Pl = 1.616e-35 m
        1 um  = 1e-6 m
        28 um = 2.8e-5 m  (Eot-Wash optimal)
        30 um = 3.0e-5 m
        0.1 mm = 1e-4 m

    Returns list of dicts, each with:
        a_0_m, a_0_label, M_6_eV, M_6_TeV, R_phys_m, lightest_5
    """
    notable = [
        (L_PL, "l_Pl (1.616e-35 m)"),
        (1.0e-6, "1 um"),
        (2.8e-5, "28 um (Eot-Wash optimal)"),
        (3.0e-5, "30 um"),
        (1.0e-4, "0.1 mm"),
    ]
    results = []
    for a_0_m, label in notable:
        spec = full_spectrum(a_0_m, l_max=10, matched=matched)
        m6 = _m6_from_a0(a_0_m)
        lightest_5 = spec["modes"][:5]
        results.append({
            "a_0_m": a_0_m,
            "a_0_label": label,
            "M_6_eV": m6,
            "M_6_TeV": m6 / EV_PER_TEV,
            "R_phys_m": spec["R_phys_m"],
            "lightest_5": lightest_5,
        })
    return results


# ===================================================================
# 7. matched_vs_unmatched
# ===================================================================

def matched_vs_unmatched(a_0_m):
    """Side-by-side comparison of matched vs unmatched scalar spectrum.

    Shows first 10 scalar modes and the constant mass enhancement factor.

    Returns dict with keys:
        matched_modes       -- first 10 scalar modes (matched)
        unmatched_modes     -- first 10 scalar modes (unmatched)
        mass_ratio          -- m_matched / m_unmatched = 1/e^{phi_vev}
        enhancement_factor  -- same as mass_ratio
        phi_vev             -- dilaton vev
        R_matched_m         -- matched physical radius
        R_unmatched_m       -- unmatched physical radius
        detectability_note  -- string summarizing impact
    """
    spec_m = scalar_spectrum(a_0_m, l_max=10, matched=True)
    spec_u = scalar_spectrum(a_0_m, l_max=10, matched=False)

    matched_modes = spec_m["modes"][:10]
    unmatched_modes = spec_u["modes"][:10]

    # Compute actual ratios for verification
    ratios = []
    for mm, um in zip(matched_modes[1:], unmatched_modes[1:]):
        if um["m_eV"] > 0:
            ratios.append(mm["m_eV"] / um["m_eV"])

    avg_ratio = sum(ratios) / len(ratios) if ratios else MASS_ENHANCEMENT

    note = (
        "Matching the dilaton vev to alpha_EM makes all KK modes {:.3f}x "
        "heavier. This shrinks the physical radius by factor {:.4f}, "
        "pushing the lightest modes to higher energies. Detection thresholds "
        "increase, but this matching is physically required for consistency "
        "with the observed electromagnetic coupling."
    ).format(MASS_ENHANCEMENT, E_PHI_VEV)

    return {
        "matched_modes": matched_modes,
        "unmatched_modes": unmatched_modes,
        "mass_ratio": avg_ratio,
        "enhancement_factor": MASS_ENHANCEMENT,
        "phi_vev": PHI_VEV,
        "R_matched_m": spec_m["R_phys_m"],
        "R_unmatched_m": spec_u["R_phys_m"],
        "detectability_note": note,
    }


# ===================================================================
# 8. cumulative_modes_vs_energy
# ===================================================================

def cumulative_modes_vs_energy(a_0_m, E_max_eV=1.0e13, n_points=100,
                               matched=True):
    """Compute N(<E) = total KK modes with mass below E.

    Sums over scalar (l>=1), vector (l>=1), graviton (l>=2).
    Compares with ADD prediction N_ADD = pi*(E*R_ADD)^2
    where R_ADD = a_0 for the torus case (single extra-dimension radius
    mapped to S^2 radius).

    Returns dict with keys:
        E_grid          -- list of energy values (eV)
        N_cumulative    -- list of cumulative mode counts
        N_add           -- list of ADD comparison counts
        R_phys_m        -- physical radius used
        R_add_m         -- ADD radius for comparison
    """
    r = _r_phys(a_0_m, matched)

    # Build full spectrum at high enough l_max
    # For S^2: m_l ~ l/R for large l, so l_max ~ E_max * R / hbar_c
    l_needed = int(E_max_eV * r / HBAR_C_EVM) + 5
    l_needed = max(l_needed, 50)
    l_needed = min(l_needed, 5000)  # cap to avoid excessive computation

    # Precompute all massive modes with mass < E_max
    all_masses_deg = []  # (mass_eV, degeneracy)

    # Scalar: l >= 1
    for l in range(1, l_needed + 1):
        m_eV = math.sqrt(l * (l + 1)) * HBAR_C_EVM / r
        if m_eV > E_max_eV:
            break
        all_masses_deg.append((m_eV, 2 * l + 1))

    # Vector: l >= 1
    for l in range(1, l_needed + 1):
        m_eV = math.sqrt(l * (l + 1)) * HBAR_C_EVM / r
        if m_eV > E_max_eV:
            break
        all_masses_deg.append((m_eV, 2 * (2 * l + 1)))

    # Graviton: l >= 2
    for l in range(2, l_needed + 1):
        m_eV = math.sqrt(l * (l + 1) - 2) * HBAR_C_EVM / r
        if m_eV > E_max_eV:
            break
        all_masses_deg.append((m_eV, 2 * l + 1))

    # Sort by mass
    all_masses_deg.sort(key=lambda x: x[0])

    # Energy grid (log-spaced)
    if E_max_eV <= 0:
        return {"E_grid": [], "N_cumulative": [], "N_add": [],
                "R_phys_m": r, "R_add_m": a_0_m}

    # Find minimum mass for log grid start
    e_min = all_masses_deg[0][0] * 0.5 if all_masses_deg else 1.0
    log_min = math.log10(max(e_min, 1.0))
    log_max = math.log10(E_max_eV)

    E_grid = []
    for i in range(n_points):
        log_e = log_min + (log_max - log_min) * i / (n_points - 1)
        E_grid.append(10.0 ** log_e)

    # Cumulative count at each energy
    N_cumulative = []
    idx = 0
    running_count = 0
    for E in E_grid:
        while idx < len(all_masses_deg) and all_masses_deg[idx][0] <= E:
            running_count += all_masses_deg[idx][1]
            idx += 1
        N_cumulative.append(running_count)

    # ADD comparison: N_ADD = pi * (E * R)^2 for 2 extra dimensions (torus)
    # Using the same a_0 as radius
    r_add_nat = a_0_m / HBAR_C_EVM  # in eV^{-1}
    N_add = []
    for E in E_grid:
        n_add = math.pi * (E * r_add_nat) ** 2
        N_add.append(n_add)

    return {
        "E_grid": E_grid,
        "N_cumulative": N_cumulative,
        "N_add": N_add,
        "R_phys_m": r,
        "R_add_m": a_0_m,
    }


# ===================================================================
# 9. summarize_kk_spectrum
# ===================================================================

def summarize_kk_spectrum():
    """Full summary report of KK spectrum with matched dilaton vev.

    Returns dict with keys:
        constants   -- key physical constants used
        key_results -- list of key physics results
        report      -- formatted string report
    """
    lines = []
    sep = "=" * 72

    lines.append(sep)
    lines.append("KK Spectrum on S^2 with Matched Dilaton Vev")
    lines.append(sep)
    lines.append("")

    # Constants
    lines.append("Physical Constants:")
    lines.append("  alpha_EM         = {:.12f}".format(ALPHA_EM))
    lines.append("  phi_vev          = {:.6f}".format(PHI_VEV))
    lines.append("  e^{{phi_vev}}      = {:.6f}".format(E_PHI_VEV))
    lines.append("  Mass enhancement = {:.6f}".format(MASS_ENHANCEMENT))
    lines.append("")

    # Key message 1: matching effect
    lines.append("-" * 72)
    lines.append("1. MATCHING EFFECT")
    lines.append("-" * 72)
    lines.append("  The matched dilaton vev phi_vev = {:.6f} shifts all".format(
        PHI_VEV))
    lines.append("  KK masses up by factor 1/e^{{phi_vev}} = {:.6f}".format(
        MASS_ENHANCEMENT))
    lines.append("  Physical radius: R_phys = a_0 * {:.6f}".format(E_PHI_VEV))
    lines.append("")

    # Key message 2: lightest modes by spin
    lines.append("-" * 72)
    lines.append("2. LIGHTEST MODES BY SPIN")
    lines.append("-" * 72)
    lines.append("  Scalar  l=1: m = sqrt(2)/R   = {:.6f}/R".format(math.sqrt(2)))
    lines.append("  Vector  l=1: m = sqrt(2)/R   = {:.6f}/R".format(math.sqrt(2)))
    lines.append("  Graviton l=2: m = sqrt(4)/R  = 2.000000/R")
    lines.append("  Fermion l=0: m = 0.5/R       = 0.500000/R")
    lines.append("")
    lines.append("  Ordering (lightest first): fermion(l=0) < scalar(l=1)")
    lines.append("    = vector(l=1) < graviton(l=2)")
    lines.append("")

    # Key message 3: Eot-Wash case study
    a_0_ew = 2.8e-5  # 28 um
    spec_matched = full_spectrum(a_0_ew, l_max=10, matched=True)
    spec_unmatched = full_spectrum(a_0_ew, l_max=10, matched=False)

    r_m = _r_phys(a_0_ew, True)
    r_u = _r_phys(a_0_ew, False)
    m6 = _m6_from_a0(a_0_ew)

    lines.append("-" * 72)
    lines.append("3. EOT-WASH CASE STUDY (a_0 = 28 um)")
    lines.append("-" * 72)
    lines.append("  M_6             = {:.4e} eV = {:.4f} TeV".format(
        m6, m6 / EV_PER_TEV))
    lines.append("  R_phys (matched)  = {:.4e} m".format(r_m))
    lines.append("  R_phys (unmatched) = {:.4e} m".format(r_u))
    lines.append("")
    lines.append("  Lightest massive mode (matched):")
    if spec_matched["lightest_massive"]:
        lm = spec_matched["lightest_massive"]
        lines.append("    spin={}, l={}, m = {:.4e} eV = {:.4f} meV".format(
            lm["spin"], lm["l"], lm["m_eV"], lm["m_eV"] * 1e3))
    lines.append("  Lightest massive mode (unmatched):")
    if spec_unmatched["lightest_massive"]:
        lu = spec_unmatched["lightest_massive"]
        lines.append("    spin={}, l={}, m = {:.4e} eV = {:.4f} meV".format(
            lu["spin"], lu["l"], lu["m_eV"], lu["m_eV"] * 1e3))
    lines.append("")

    # Key message 4: discrete spectrum
    lines.append("-" * 72)
    lines.append("4. DISCRETE SPECTRUM WITH GROWING DEGENERACY")
    lines.append("-" * 72)
    lines.append("  Unlike ADD's quasi-continuous KK lattice on T^n,")
    lines.append("  the S^2 spectrum is discrete with degeneracy (2l+1).")
    lines.append("  This gives a distinctive pattern in cross-sections:")
    lines.append("  resonance peaks at specific energies, not a continuum.")
    lines.append("")

    # Key message 5: detection implications
    lines.append("-" * 72)
    lines.append("5. DETECTION IMPLICATIONS")
    lines.append("-" * 72)
    lines.append("  Matching alpha_EM makes all modes 1.817x heavier.")
    lines.append("  Detection thresholds increase, but this matching is")
    lines.append("  physically required for consistency with observed")
    lines.append("  electromagnetic coupling.")
    lines.append("")

    # Notable a_0 table
    lines.append("-" * 72)
    lines.append("6. SPECTRUM AT NOTABLE a_0 VALUES")
    lines.append("-" * 72)
    notable = spectrum_at_notable_a0(matched=True)
    for entry in notable:
        lines.append("")
        lines.append("  a_0 = {}".format(entry["a_0_label"]))
        lines.append("    M_6     = {:.4e} eV ({:.4e} TeV)".format(
            entry["M_6_eV"], entry["M_6_TeV"]))
        lines.append("    R_phys  = {:.4e} m".format(entry["R_phys_m"]))
        lines.append("    Lightest 5 massive modes:")
        for mode in entry["lightest_5"]:
            lines.append("      spin={:8s}  l={}  m = {:.4e} eV  deg = {}".format(
                mode["spin"], mode["l"], mode["m_eV"], mode["degeneracy"]))
    lines.append("")

    # Matched vs unmatched table
    lines.append("-" * 72)
    lines.append("7. MATCHED vs UNMATCHED (a_0 = 28 um)")
    lines.append("-" * 72)
    comp = matched_vs_unmatched(a_0_ew)
    lines.append("  Mass ratio: {:.6f}".format(comp["mass_ratio"]))
    lines.append("  R_matched   = {:.4e} m".format(comp["R_matched_m"]))
    lines.append("  R_unmatched = {:.4e} m".format(comp["R_unmatched_m"]))
    lines.append("")
    lines.append("  {:>3s}  {:>14s}  {:>14s}  {:>8s}".format(
        "l", "m_matched(eV)", "m_unmatched(eV)", "ratio"))
    for mm, um in zip(comp["matched_modes"], comp["unmatched_modes"]):
        ratio = mm["m_eV"] / um["m_eV"] if um["m_eV"] > 0 else 0.0
        lines.append("  {:3d}  {:14.6e}  {:14.6e}  {:8.4f}".format(
            mm["l"], mm["m_eV"], um["m_eV"], ratio))
    lines.append("")

    lines.append(sep)
    lines.append("END OF REPORT")
    lines.append(sep)

    report = "\n".join(lines)

    constants = {
        "alpha_EM": ALPHA_EM,
        "phi_vev": PHI_VEV,
        "e_phi_vev": E_PHI_VEV,
        "mass_enhancement": MASS_ENHANCEMENT,
        "M_Pl_eV": M_PL_EV,
        "l_Pl": L_PL,
        "hbar_c_eVm": HBAR_C_EVM,
    }

    key_results = [
        "Matched dilaton vev shifts all KK masses up by factor {:.6f}".format(
            MASS_ENHANCEMENT),
        "Lightest modes: fermion(l=0) at 0.5/R, scalar/vector(l=1) at sqrt(2)/R, graviton(l=2) at 2/R",
        "For a_0=28um: lightest mode ~ {:.1f} meV (matched) vs {:.1f} meV (unmatched)".format(
            spec_matched["lightest_massive"]["m_eV"] * 1e3 if spec_matched["lightest_massive"] else 0,
            spec_unmatched["lightest_massive"]["m_eV"] * 1e3 if spec_unmatched["lightest_massive"] else 0),
        "S^2 spectrum is discrete with (2l+1) degeneracy, unlike ADD's quasi-continuous lattice",
        "Matching alpha_EM is physically required despite making detection harder",
    ]

    return {
        "constants": constants,
        "key_results": key_results,
        "report": report,
    }


# ===================================================================
# __main__
# ===================================================================

if __name__ == "__main__":
    result = summarize_kk_spectrum()
    print(result["report"])
