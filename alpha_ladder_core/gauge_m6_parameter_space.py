"""
(a_0, M_6) Parameter Space with Dilaton Vev Constraint
======================================================

Maps the one-dimensional parameter space of the Alpha Ladder S^2
compactification.  The dilaton vev phi_vev = (1/4) ln(4 pi alpha_EM)
fixes the gauge coupling to alpha_EM = 1/137 but leaves the bare
compactification radius a_0 as a free parameter.

Given a_0, the 6D Planck mass M_6 is determined by:
    M_Pl^2 = 4 pi a_0^2 M_6^4

The physical radius is R_phys = a_0 e^{phi_vev} (45% smaller than a_0),
and the KK mass scale is m_KK ~ hbar c / R_phys.

Pure Python -- only ``import math``, no numpy/scipy.

Functions
---------
1. compute_parameter_point  -- full physics at a single a_0
2. parameter_space_scan     -- log-spaced scan from l_Pl to 1 mm
3. notable_points           -- 11 physically interesting a_0 values
4. experimental_reach       -- which experiments probe which (a_0, M_6)
5. dilaton_vev_effects      -- quantify vev vs unmatched (phi=0) case
6. can_anything_fix_a0      -- survey of mechanisms to fix a_0
7. summarize_parameter_space -- full report
"""

import math

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
ALPHA_EM = 1.0 / 137.035999084            # CODATA 2018
G_N = 6.674298e-11                         # m^3 kg^-1 s^-2
C_LIGHT = 2.99792458e8                     # m s^-1
HBAR = 1.054571817e-34                     # J s
M_PL_EV = 1.22089e28                       # Planck mass (eV)
L_PL = 1.61625e-35                         # Planck length (m)
HBAR_C_EVM = 1.9733e-7                     # hbar*c in eV m
PHI_GOLDEN = (1.0 + math.sqrt(5.0)) / 2.0 # golden ratio

# Derived
PHI_VEV = 0.25 * math.log(4.0 * math.pi * ALPHA_EM)   # dilaton vev
E_PHI_VEV = math.exp(PHI_VEV)                          # e^{phi_vev} ~ 0.5503
EV_PER_TEV = 1.0e12                                    # eV per TeV
EV_PER_GEV = 1.0e9                                     # eV per GeV
EV_PER_MEV = 1.0e6                                     # eV per MeV


# ===================================================================
# helpers
# ===================================================================
def _m6_from_a0(a_0_m):
    """Compute M_6 in eV from a_0 in meters.

    M_Pl^2 = 4 pi a_0^2 M_6^4
    => M_6 = (M_Pl^2 / (4 pi a_0^2))^{1/4}

    a_0 must be in natural units (eV^{-1}), so convert:
        a_0_nat = a_0_m / hbar_c_eVm   (meters -> eV^{-1})
    """
    a_0_nat = a_0_m / HBAR_C_EVM  # eV^{-1}
    m6_fourth = M_PL_EV ** 2 / (4.0 * math.pi * a_0_nat ** 2)
    m6_eV = m6_fourth ** 0.25
    return m6_eV


def _classify_regime(m6_eV):
    """Classify the M_6 energy regime."""
    m6_TeV = m6_eV / EV_PER_TEV
    if m6_TeV > 1.0e16:
        return "Planck"
    elif m6_TeV > 100.0:
        return "above_LHC"
    elif m6_TeV > 1.0:
        return "TeV"
    elif m6_eV > EV_PER_GEV:
        return "GeV"
    elif m6_eV > EV_PER_MEV:
        return "MeV"
    elif m6_eV > 1.0:
        return "eV"
    elif m6_eV > 1.0e-3:
        return "meV"
    else:
        return "sub_meV"


def _format_eV(val_eV):
    """Format an eV value in the most natural unit."""
    if val_eV >= EV_PER_TEV:
        return "{:.4e} TeV".format(val_eV / EV_PER_TEV)
    elif val_eV >= EV_PER_GEV:
        return "{:.4e} GeV".format(val_eV / EV_PER_GEV)
    elif val_eV >= EV_PER_MEV:
        return "{:.4e} MeV".format(val_eV / EV_PER_MEV)
    elif val_eV >= 1.0:
        return "{:.4e} eV".format(val_eV)
    elif val_eV >= 1.0e-3:
        return "{:.4e} meV".format(val_eV * 1.0e3)
    else:
        return "{:.4e} eV".format(val_eV)


def _format_m(val_m):
    """Format a length in the most natural unit."""
    if val_m >= 1.0e-3:
        return "{:.4e} mm".format(val_m * 1.0e3)
    elif val_m >= 1.0e-6:
        return "{:.4e} um".format(val_m * 1.0e6)
    elif val_m >= 1.0e-9:
        return "{:.4e} nm".format(val_m * 1.0e9)
    elif val_m >= 1.0e-15:
        return "{:.4e} fm".format(val_m * 1.0e15)
    else:
        return "{:.4e} m".format(val_m)


# ===================================================================
# 1. compute_parameter_point
# ===================================================================
def compute_parameter_point(a_0_m):
    """Compute full physics for a given a_0 in meters.

    Returns a dict with M_6, R_phys, m_KK, regime classification,
    and comparison with the unmatched (phi=0) case.
    """
    # M_6 from Planck mass relation
    m6_eV = _m6_from_a0(a_0_m)
    m6_TeV = m6_eV / EV_PER_TEV

    # Physical radius with dilaton vev
    R_phys = a_0_m * E_PHI_VEV
    R_bare = a_0_m  # unmatched case

    # KK mass scale: m_KK = hbar*c / R_phys
    m_KK_eV = HBAR_C_EVM / R_phys
    m_KK_bare_eV = HBAR_C_EVM / R_bare

    # l=1 mode: m_{l=1} = sqrt(2) / R_phys (in natural units)
    m_KK_l1_eV = math.sqrt(2.0) * m_KK_eV

    # regime classification
    regime = _classify_regime(m6_eV)

    # ratio of a_0 to Planck length
    a0_over_lPl = a_0_m / L_PL

    return {
        "a_0_m": a_0_m,
        "a_0_over_l_Pl": a0_over_lPl,
        "M_6_eV": m6_eV,
        "M_6_TeV": m6_TeV,
        "R_phys_m": R_phys,
        "R_bare_m": R_bare,
        "m_KK_eV": m_KK_eV,
        "m_KK_l1_eV": m_KK_l1_eV,
        "m_KK_bare_eV": m_KK_bare_eV,
        "regime": regime,
        "R_phys_over_R_bare": E_PHI_VEV,
        "m_KK_over_m_KK_bare": 1.0 / E_PHI_VEV,
    }


# ===================================================================
# 2. parameter_space_scan
# ===================================================================
def parameter_space_scan(a_0_values=None):
    """Scan over a range of a_0 values (log-spaced from l_Pl to 1 mm).

    Returns list of parameter point dicts plus summary statistics.
    """
    if a_0_values is None:
        # 30 log-spaced points from l_Pl to 1 mm
        log_min = math.log10(L_PL)
        log_max = math.log10(1.0e-3)
        n_points = 30
        step = (log_max - log_min) / (n_points - 1)
        a_0_values = [10.0 ** (log_min + i * step) for i in range(n_points)]

    points = []
    for a0 in a_0_values:
        pt = compute_parameter_point(a0)
        points.append(pt)

    # summary statistics
    m6_values = [p["M_6_eV"] for p in points]
    mkk_values = [p["m_KK_eV"] for p in points]

    summary = {
        "n_points": len(points),
        "a_0_min_m": min(p["a_0_m"] for p in points),
        "a_0_max_m": max(p["a_0_m"] for p in points),
        "M_6_min_eV": min(m6_values),
        "M_6_max_eV": max(m6_values),
        "M_6_min_TeV": min(m6_values) / EV_PER_TEV,
        "M_6_max_TeV": max(m6_values) / EV_PER_TEV,
        "m_KK_min_eV": min(mkk_values),
        "m_KK_max_eV": max(mkk_values),
    }

    return {
        "points": points,
        "summary": summary,
    }


# ===================================================================
# 3. notable_points
# ===================================================================
def notable_points():
    """Compute parameter space at 11 physically interesting a_0 values."""
    entries = [
        (L_PL, "Planck scale (minimal radius)"),
        (10.0 * L_PL, "10 x Planck length"),
        (100.0 * L_PL, "100 x Planck length"),
        (1.0e-9, "1 nm (atomic scale)"),
        (1.0e-6, "1 um (sub-mm experiments begin)"),
        (10.0e-6, "10 um (torsion pendulum range)"),
        (28.2e-6, "28.2 um (Eot-Wash optimal, m_phi ~ 7 meV)"),
        (30.0e-6, "30 um (round number near Eot-Wash)"),
        (56.0e-6, "56 um (Eot-Wash 2006 reach)"),
        (0.1e-3, "0.1 mm (sub-mm gravity frontier)"),
        (1.0e-3, "1 mm (macroscopic scale)"),
    ]

    results = []
    for a0, desc in entries:
        pt = compute_parameter_point(a0)
        pt["description"] = desc
        results.append(pt)

    return results


# ===================================================================
# 4. experimental_reach
# ===================================================================
def experimental_reach():
    """Map which experiments probe which part of the (a_0, M_6) space.

    For each experiment, compute the a_0 threshold from the observable
    constraint and the corresponding M_6 range.
    """
    experiments = []

    # Helper: given M_6 in eV, find a_0 in meters
    # M_6 = (M_Pl^2 / (4 pi a_0_nat^2))^{1/4}
    # => a_0_nat = M_Pl / (sqrt(4 pi) M_6^2)^{1/2} ... let's invert properly
    # M_6^4 = M_Pl^2 / (4 pi a_0_nat^2)
    # a_0_nat^2 = M_Pl^2 / (4 pi M_6^4)
    # a_0_nat = M_Pl / (2 sqrt(pi) M_6^2)
    # a_0_m = a_0_nat * hbar_c_eVm

    def a0_from_m6(m6_eV):
        a0_nat = M_PL_EV / (2.0 * math.sqrt(math.pi) * m6_eV ** 2)
        return a0_nat * HBAR_C_EVM

    def m_kk_from_a0(a0_m):
        R_phys = a0_m * E_PHI_VEV
        return HBAR_C_EVM / R_phys

    # 1. LHC (14 TeV): can produce KK modes if M_6 < 7 TeV
    m6_lhc = 7.0e12  # 7 TeV in eV
    a0_lhc = a0_from_m6(m6_lhc)
    mkk_lhc = m_kk_from_a0(a0_lhc)
    experiments.append({
        "experiment": "LHC (14 TeV)",
        "observable": "Missing energy from KK graviton production",
        "constraint": "M_6 < 7 TeV",
        "a_0_threshold_m": a0_lhc,
        "M_6_at_threshold_eV": m6_lhc,
        "M_6_at_threshold_TeV": m6_lhc / EV_PER_TEV,
        "m_KK_at_threshold_eV": mkk_lhc,
        "a_0_direction": "a_0 > {:.2e} m".format(a0_lhc),
        "status": "current",
        "verdict": (
            "Requires a_0 > {:.2e} m ({:.1f} um).  "
            "This is in the sub-mm range accessible to torsion "
            "balance experiments.".format(a0_lhc, a0_lhc * 1.0e6)
        ),
    })

    # 2. FCC (100 TeV): M_6 < 50 TeV
    m6_fcc = 50.0e12  # 50 TeV in eV
    a0_fcc = a0_from_m6(m6_fcc)
    mkk_fcc = m_kk_from_a0(a0_fcc)
    experiments.append({
        "experiment": "FCC (100 TeV)",
        "observable": "Missing energy / resonances from KK tower",
        "constraint": "M_6 < 50 TeV",
        "a_0_threshold_m": a0_fcc,
        "M_6_at_threshold_eV": m6_fcc,
        "M_6_at_threshold_TeV": m6_fcc / EV_PER_TEV,
        "m_KK_at_threshold_eV": mkk_fcc,
        "a_0_direction": "a_0 > {:.2e} m".format(a0_fcc),
        "status": "planned",
        "verdict": (
            "Requires a_0 > {:.2e} m ({:.1f} nm).  "
            "Much smaller radius than torsion balance reach; "
            "complementary probe.".format(a0_fcc, a0_fcc * 1.0e9)
        ),
    })

    # 3. Eot-Wash torsion balance: sensitive to m_KK < 7 meV
    #    m_KK = hbar*c / R_phys < 7e-3 eV
    #    R_phys > hbar*c / 7e-3 = 2.82e-5 m
    #    a_0 > R_phys / e^{phi_vev} = 2.82e-5 / 0.5503
    m_kk_eotwash = 7.0e-3  # 7 meV in eV
    R_phys_eotwash = HBAR_C_EVM / m_kk_eotwash
    a0_eotwash = R_phys_eotwash / E_PHI_VEV
    m6_eotwash = _m6_from_a0(a0_eotwash)
    experiments.append({
        "experiment": "Eot-Wash torsion balance",
        "observable": "Yukawa deviation from 1/r^2 at sub-mm distances",
        "constraint": "m_KK < 7 meV (lambda > 28 um)",
        "a_0_threshold_m": a0_eotwash,
        "M_6_at_threshold_eV": m6_eotwash,
        "M_6_at_threshold_TeV": m6_eotwash / EV_PER_TEV,
        "m_KK_at_threshold_eV": m_kk_eotwash,
        "a_0_direction": "a_0 > {:.2e} m".format(a0_eotwash),
        "status": "current",
        "verdict": (
            "Requires a_0 > {:.2e} m ({:.1f} um).  "
            "This is the most sensitive current probe of the "
            "parameter space.  M_6 ~ {:.1f} TeV at threshold.".format(
                a0_eotwash, a0_eotwash * 1.0e6,
                m6_eotwash / EV_PER_TEV
            )
        ),
    })

    # 4. Cassini PPN: constrains dilaton coupling via Shapiro delay
    #    Requires massive dilaton; constrains coupling strength
    #    Not a direct a_0 threshold but mass-dependent
    experiments.append({
        "experiment": "Cassini PPN (Shapiro delay)",
        "observable": "Post-Newtonian parameter gamma deviation",
        "constraint": "|gamma - 1| < 2.3e-5",
        "a_0_threshold_m": None,
        "M_6_at_threshold_eV": None,
        "M_6_at_threshold_TeV": None,
        "m_KK_at_threshold_eV": None,
        "a_0_direction": "Constrains dilaton coupling, not a_0 directly",
        "status": "current",
        "verdict": (
            "Cassini constrains the dilaton-matter coupling.  "
            "For Planck-mass dilaton (flux-stabilized scenario), "
            "the dilaton decouples and Cassini is automatically "
            "satisfied regardless of a_0."
        ),
    })

    # 5. Neutron interferometry: sub-mm Yukawa forces
    #    Sensitive to m_KK in the eV range (nm distances)
    m_kk_neutron = 1.0  # ~1 eV (nm scale)
    R_phys_neutron = HBAR_C_EVM / m_kk_neutron
    a0_neutron = R_phys_neutron / E_PHI_VEV
    m6_neutron = _m6_from_a0(a0_neutron)
    experiments.append({
        "experiment": "Neutron interferometry",
        "observable": "Short-range Yukawa forces at nm-um distances",
        "constraint": "m_KK ~ 0.1-100 eV (nm-um range)",
        "a_0_threshold_m": a0_neutron,
        "M_6_at_threshold_eV": m6_neutron,
        "M_6_at_threshold_TeV": m6_neutron / EV_PER_TEV,
        "m_KK_at_threshold_eV": m_kk_neutron,
        "a_0_direction": "a_0 ~ {:.2e} m".format(a0_neutron),
        "status": "current",
        "verdict": (
            "Probes a_0 ~ {:.2e} m (nm scale), corresponding "
            "to M_6 ~ {:.2e} TeV.  Well above LHC reach in M_6 "
            "but accessible to precision force measurements.".format(
                a0_neutron, m6_neutron / EV_PER_TEV
            )
        ),
    })

    # 6. ATLAS/CMS missing energy: direct KK graviton production
    #    Similar to LHC entry but for graviton-specific channels
    m6_atlas = 5.0e12  # 5 TeV current limit
    a0_atlas = a0_from_m6(m6_atlas)
    mkk_atlas = m_kk_from_a0(a0_atlas)
    experiments.append({
        "experiment": "ATLAS/CMS missing energy (graviton)",
        "observable": "Mono-jet + MET from KK graviton emission",
        "constraint": "M_6 < 5 TeV (current 95% CL exclusion, n=2)",
        "a_0_threshold_m": a0_atlas,
        "M_6_at_threshold_eV": m6_atlas,
        "M_6_at_threshold_TeV": m6_atlas / EV_PER_TEV,
        "m_KK_at_threshold_eV": mkk_atlas,
        "a_0_direction": "a_0 > {:.2e} m".format(a0_atlas),
        "status": "current",
        "verdict": (
            "Current ATLAS/CMS limits exclude M_6 < 5 TeV for n=2 "
            "extra dimensions, corresponding to a_0 > {:.2e} m "
            "({:.1f} um).".format(a0_atlas, a0_atlas * 1.0e6)
        ),
    })

    return experiments


# ===================================================================
# 5. dilaton_vev_effects
# ===================================================================
def dilaton_vev_effects():
    """Quantify what phi_vev = -0.597 changes vs the unmatched phi=0 case."""
    e_phi = E_PHI_VEV
    e_inv = 1.0 / e_phi

    # Unmatched gauge coupling: alpha_unmatched = e^{4*0} / (4 pi) = 1/(4 pi)
    alpha_unmatched = 1.0 / (4.0 * math.pi)
    alpha_matched = ALPHA_EM

    return {
        "phi_vev": PHI_VEV,
        "e_to_phi_vev": e_phi,
        "R_phys_over_a_0": e_phi,
        "radius_reduction_percent": (1.0 - e_phi) * 100.0,
        "m_KK_enhancement_factor": e_inv,
        "m_KK_enhancement_percent": (e_inv - 1.0) * 100.0,
        "M_6_unchanged": True,
        "M_6_depends_on": "a_0 only (M_Pl^2 = 4 pi a_0^2 M_6^4, no phi dependence)",
        "alpha_matched": alpha_matched,
        "alpha_unmatched": alpha_unmatched,
        "alpha_ratio": alpha_matched / alpha_unmatched,
        "detection_impact": (
            "The dilaton vev makes KK modes 82% heavier than the "
            "unmatched case.  This pushes them further from experimental "
            "reach at any given a_0.  Detection is HARDER with the vev, "
            "not easier.  This is an honest consequence of matching "
            "alpha_EM -- the physical radius shrinks."
        ),
        "gauge_coupling_impact": (
            "Without the vev (phi=0), the KK gauge coupling gives "
            "alpha = 1/(4 pi) = {:.6f}, roughly 11x larger than alpha_EM.  "
            "The vev suppresses the coupling to match the observed "
            "fine-structure constant.".format(alpha_unmatched)
        ),
    }


# ===================================================================
# 6. can_anything_fix_a0
# ===================================================================
def can_anything_fix_a0():
    """Survey of mechanisms that might provide a second equation to fix a_0."""
    mechanisms = [
        {
            "mechanism": "Alpha running (SM RG)",
            "can_fix_a0": False,
            "verdict": (
                "SM running accounts for the observed alpha(m_Z) = 1/128.  "
                "It does not provide an independent constraint on a_0."
            ),
            "reference": None,
        },
        {
            "mechanism": "Coleman-Weinberg potential",
            "can_fix_a0": False,
            "verdict": (
                "CW potential from charged KK scalars + vectors does not "
                "develop a minimum that fixes a_0.  See kk_gauge_matching.py."
            ),
            "reference": "kk_gauge_matching.py",
        },
        {
            "mechanism": "Flux quantization",
            "can_fix_a0": False,
            "verdict": (
                "Flux quantization N = integer provides a discrete label "
                "but does not fix the continuous parameter a_0.  The flux "
                "potential has a minimum for each N but that minimum's "
                "location depends on unknown coefficients."
            ),
            "reference": "flux_stabilization.py (Alpha-Ladder)",
        },
        {
            "mechanism": "Casimir energy",
            "can_fix_a0": False,
            "verdict": (
                "Casimir energy on S^2 gives a dilaton mass but the "
                "potential V_min(a_0) is monotonic -- no minimum that "
                "selects a preferred a_0.  Casimir alone yields a no-go."
            ),
            "reference": "casimir_stabilization.py (Alpha-Ladder)",
        },
        {
            "mechanism": "Anthropic selection",
            "can_fix_a0": False,
            "verdict": (
                "Anthropic reasoning constrains a_0 to the range where "
                "stable atoms exist, but this is not a dynamical mechanism.  "
                "It gives a broad window, not a unique prediction."
            ),
            "reference": None,
        },
        {
            "mechanism": "Higher-derivative terms (R^2, GB)",
            "can_fix_a0": False,
            "verdict": (
                "Higher-derivative corrections to the 6D action could in "
                "principle generate a potential for a_0, but the coefficients "
                "are unknown.  The Gauss-Bonnet term is topological for S^2 "
                "and does not help."
            ),
            "reference": "kk_reduction.py (Alpha-Ladder)",
        },
        {
            "mechanism": "Vacuum polynomial / scaling symmetry",
            "can_fix_a0": False,
            "verdict": (
                "The vacuum polynomial x^2 + Dx + d = 0 selects D=6 but "
                "does not fix a_0.  Scaling symmetry of the EH action "
                "makes the exponent of a_0 exactly zero."
            ),
            "reference": "radius_determination.py (Alpha-Ladder)",
        },
        {
            "mechanism": "Flux + Casimir combined",
            "can_fix_a0": False,
            "verdict": (
                "Combining flux and Casimir terms gives a balance condition "
                "that determines the dilaton mass but still leaves a_0 as "
                "a flat direction due to scaling symmetry."
            ),
            "reference": "radius_determination.py (Alpha-Ladder)",
        },
    ]

    n_tried = len(mechanisms)
    n_succeed = sum(1 for m in mechanisms if m["can_fix_a0"])

    return {
        "mechanisms": mechanisms,
        "n_tried": n_tried,
        "n_succeed": n_succeed,
        "overall_assessment": (
            "{} mechanisms examined, {} succeed in fixing a_0.  "
            "The compactification radius remains a free parameter "
            "of the framework.  This is an honest gap: the theory "
            "predicts G as a function of a_0 (and alpha, mu, phi) "
            "but does not determine a_0 from first principles.".format(
                n_tried, n_succeed
            )
        ),
    }


# ===================================================================
# 7. summarize_parameter_space
# ===================================================================
def summarize_parameter_space():
    """Full summary of the (a_0, M_6) parameter space."""
    vev = dilaton_vev_effects()
    notable = notable_points()
    expt = experimental_reach()
    fix = can_anything_fix_a0()

    # Find Eot-Wash point
    eotwash_pt = None
    for pt in notable:
        if "28.2" in pt.get("description", ""):
            eotwash_pt = pt
            break

    return {
        "key_messages": [
            (
                "phi_vev = {:.6f} fixes alpha_EM = 1/137 but leaves a_0 "
                "as a free parameter.".format(PHI_VEV)
            ),
            (
                "The parameter space is one-dimensional: choose a_0, "
                "get M_6 = (M_Pl^2 / (4 pi a_0^2))^{{1/4}}."
            ),
            (
                "The Eot-Wash window (a_0 ~ 30 um, M_6 ~ {:.1f} TeV) "
                "is the most experimentally interesting point.".format(
                    eotwash_pt["M_6_TeV"] if eotwash_pt else 5.0
                )
            ),
            (
                "The dilaton vev makes KK modes {:.0f}% heavier "
                "(harder to detect).".format(vev["m_KK_enhancement_percent"])
            ),
            (
                "No known mechanism fixes a_0 within the framework "
                "({} mechanisms tried, all fail).".format(fix["n_tried"])
            ),
        ],
        "dilaton_vev_effects": vev,
        "notable_points": notable,
        "experimental_reach": expt,
        "mechanisms_to_fix_a0": fix,
    }


# ===================================================================
# __main__ -- print full report
# ===================================================================
if __name__ == "__main__":

    W = 72

    print("=" * W)
    print("(a_0, M_6) Parameter Space with Dilaton Vev Constraint")
    print("Alpha Ladder S^2 Framework")
    print("=" * W)

    # --- Dilaton vev ---
    print("\n--- Dilaton Vev Effects ---")
    vev = dilaton_vev_effects()
    print("  phi_vev                = {:.8f}".format(vev["phi_vev"]))
    print("  e^(phi_vev)            = {:.8f}".format(vev["e_to_phi_vev"]))
    print("  R_phys / a_0           = {:.8f}  ({:.1f}% reduction)".format(
        vev["R_phys_over_a_0"], vev["radius_reduction_percent"]))
    print("  m_KK enhancement       = {:.4f}x  ({:.1f}% heavier)".format(
        vev["m_KK_enhancement_factor"], vev["m_KK_enhancement_percent"]))
    print("  alpha (matched)        = {:.10e}".format(vev["alpha_matched"]))
    print("  alpha (unmatched phi=0)= {:.10e}".format(vev["alpha_unmatched"]))
    print("  Detection impact: {}".format(vev["detection_impact"]))

    # --- Notable points ---
    print("\n" + "=" * W)
    print("Notable Points in Parameter Space")
    print("=" * W)
    notable = notable_points()

    header = "  {:>12s}  {:>12s}  {:>12s}  {:>14s}  {:>14s}  {:>10s}".format(
        "a_0", "a_0/l_Pl", "M_6", "R_phys", "m_KK", "regime")
    print(header)
    print("  " + "-" * (len(header) - 2))

    for pt in notable:
        print("  {:>12s}  {:>12.2e}  {:>12s}  {:>14s}  {:>14s}  {:>10s}".format(
            _format_m(pt["a_0_m"]),
            pt["a_0_over_l_Pl"],
            _format_eV(pt["M_6_eV"]),
            _format_m(pt["R_phys_m"]),
            _format_eV(pt["m_KK_eV"]),
            pt["regime"],
        ))

    print("\n  Descriptions:")
    for pt in notable:
        print("    a_0 = {:>12s}  -- {}".format(
            _format_m(pt["a_0_m"]), pt["description"]))

    # --- Parameter space scan ---
    print("\n" + "=" * W)
    print("Full Parameter Space Scan (l_Pl to 1 mm)")
    print("=" * W)
    scan = parameter_space_scan()
    s = scan["summary"]
    print("  {} points scanned".format(s["n_points"]))
    print("  a_0 range : {:.2e} m  to  {:.2e} m".format(
        s["a_0_min_m"], s["a_0_max_m"]))
    print("  M_6 range : {:.2e} eV ({:.2e} TeV)  to  {:.2e} eV ({:.2e} TeV)".format(
        s["M_6_min_eV"], s["M_6_min_TeV"],
        s["M_6_max_eV"], s["M_6_max_TeV"]))
    print("  m_KK range: {:.2e} eV  to  {:.2e} eV".format(
        s["m_KK_min_eV"], s["m_KK_max_eV"]))

    print("\n  {:>12s}  {:>14s}  {:>14s}  {:>14s}  {:>10s}".format(
        "a_0 (m)", "M_6 (eV)", "R_phys (m)", "m_KK (eV)", "regime"))
    print("  " + "-" * 68)
    for pt in scan["points"]:
        print("  {:>12.3e}  {:>14.4e}  {:>14.4e}  {:>14.4e}  {:>10s}".format(
            pt["a_0_m"], pt["M_6_eV"], pt["R_phys_m"],
            pt["m_KK_eV"], pt["regime"]))

    # --- Experimental reach ---
    print("\n" + "=" * W)
    print("Experimental Reach")
    print("=" * W)
    expts = experimental_reach()
    for ex in expts:
        print("\n  [{}]  ({})".format(ex["experiment"], ex["status"]))
        print("    Observable  : {}".format(ex["observable"]))
        print("    Constraint  : {}".format(ex["constraint"]))
        if ex["a_0_threshold_m"] is not None:
            print("    a_0 threshold: {:.3e} m ({})".format(
                ex["a_0_threshold_m"], ex["a_0_direction"]))
            print("    M_6 at threshold: {} ({:.2f} TeV)".format(
                _format_eV(ex["M_6_at_threshold_eV"]),
                ex["M_6_at_threshold_TeV"]))
            print("    m_KK at threshold: {}".format(
                _format_eV(ex["m_KK_at_threshold_eV"])))
        else:
            print("    a_0 threshold: {}".format(ex["a_0_direction"]))
        print("    Verdict: {}".format(ex["verdict"]))

    # --- Can anything fix a_0? ---
    print("\n" + "=" * W)
    print("Can Anything Fix a_0?")
    print("=" * W)
    fix = can_anything_fix_a0()
    print("\n  {} mechanisms examined, {} succeed.\n".format(
        fix["n_tried"], fix["n_succeed"]))

    for m in fix["mechanisms"]:
        status = "YES" if m["can_fix_a0"] else "NO"
        ref = "  [{}]".format(m["reference"]) if m["reference"] else ""
        print("  [{:3s}] {}{}".format(status, m["mechanism"], ref))
        print("        {}".format(m["verdict"]))

    print("\n  Overall: {}".format(fix["overall_assessment"]))

    # --- Summary ---
    print("\n" + "=" * W)
    print("SUMMARY")
    print("=" * W)
    summary = summarize_parameter_space()
    for i, msg in enumerate(summary["key_messages"], 1):
        print("\n  {}. {}".format(i, msg))

    # --- Matched vs unmatched comparison table ---
    print("\n" + "=" * W)
    print("Matched (phi_vev) vs Unmatched (phi=0) Comparison")
    print("=" * W)
    print("\n  At a_0 = 28.2 um (Eot-Wash optimal):")
    pt_ew = compute_parameter_point(28.2e-6)
    print("    R_phys (matched)   = {:.4e} m".format(pt_ew["R_phys_m"]))
    print("    R_bare (unmatched) = {:.4e} m".format(pt_ew["R_bare_m"]))
    print("    m_KK (matched)     = {} (heavier)".format(
        _format_eV(pt_ew["m_KK_eV"])))
    print("    m_KK (unmatched)   = {} (lighter)".format(
        _format_eV(pt_ew["m_KK_bare_eV"])))
    print("    M_6                = {} (same either way)".format(
        _format_eV(pt_ew["M_6_eV"])))
    print("    regime             = {}".format(pt_ew["regime"]))

    print("\n" + "=" * W)
    print("HONEST ASSESSMENT")
    print("=" * W)
    print("""
  The (a_0, M_6) parameter space is one-dimensional.  The dilaton vev
  phi_vev = {:.6f} is determined by alpha_EM matching, but a_0 remains
  free.  This is not a failure of the framework -- it is a consequence
  of the scaling symmetry of the Einstein-Hilbert action on S^2.

  The most experimentally interesting region is a_0 ~ 30 um, where
  M_6 ~ 5 TeV (accessible to colliders) and m_KK ~ 7 meV (accessible
  to torsion balance experiments).  However, the dilaton vev makes KK
  modes {:.0f}% heavier than the unmatched case, pushing detection
  thresholds higher.

  Eight mechanisms have been examined to fix a_0; all fail.  This
  remains an open theoretical gap.
""".format(PHI_VEV, vev["m_KK_enhancement_percent"]))

    print("=" * W)
    print("END OF REPORT")
    print("=" * W)
