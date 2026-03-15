"""
Theoretical context for the alpha ladder framework.

Structures the tree-level vs empirical screening discrepancy, positions
the framework in the literature, and honestly assesses anomaly status.
"""

import math


def compute_screening_discrepancy(constants):
    """Structure the tree-level vs empirical amplitude discrepancy."""
    from alpha_ladder_core.dilaton import compute_bd_parameter
    from alpha_ladder_core.screening import compute_screening_parameters

    bd = compute_bd_parameter(constants)
    omega = bd["omega"]
    denom = 2.0 * omega + 3.0
    alpha_tree = 2.0 / denom if denom > 0 else 0.0

    sp = compute_screening_parameters(constants)
    alpha_empirical = sp["alpha_screening"]

    ratio = alpha_tree / alpha_empirical if alpha_empirical != 0 else float('inf')
    log10_ratio = math.log10(abs(ratio)) if ratio != 0 and not math.isinf(ratio) else float('inf')

    resolutions = [
        {
            "name": "Heavy dilaton (lambda < 0.1 mm)",
            "mechanism": "If m_phi > ~2e-3 eV (Compton wavelength < 0.1 mm), the Yukawa fifth force is exponentially suppressed at all lab and solar system scales. Eot-Wash allows alpha = 0.618 for lambda < 0.1 mm. The discrepancy dissolves.",
            "plausibility": "high",
            "description": "The 3854x ratio assumes the 160 ppm G excess (G_CODATA vs G_vacuum) IS dilaton screening. But the 7 best G measurements span 530 ppm, and 160 ppm is only 30% of that spread -- well within experimental scatter. If the dilaton is heavy enough that its Compton wavelength is below ~0.1 mm, the Yukawa force is undetectable at all experimental scales. The alpha_screening in screening.py was a phenomenological fit, not a measurement of actual dilaton screening. Testable: next-generation sub-mm gravity experiments probing below 50 microns.",
        },
        {
            "name": "Environment-dependent mass (chameleon)",
            "mechanism": "Chameleon or symmetron mechanism: the dilaton mass increases in dense environments, further suppressing the lab-scale fifth force beyond the Yukawa cutoff",
            "plausibility": "moderate",
            "description": "If the dilaton potential has a density-dependent minimum, the effective mass in the lab (rho ~ 1 g/cm^3) could be much larger than in vacuum or space. This could provide additional suppression even if the vacuum mass alone is marginal. However, this requires a specific form of the potential that is not derived from the 6D action.",
        },
        {
            "name": "Moduli mixing",
            "mechanism": "The 4 scalar moduli (1 dilaton + 3 from symmetric 2x2 internal metric on S^2) have mass eigenstates that differ from coupling eigenstates, giving partial cancellation in the net matter coupling",
            "plausibility": "moderate",
            "description": "The physical dilaton (mass eigenstate) is a superposition of coupling eigenstates. With 4 moduli, the matter-coupled combination can have a suppressed projection onto the mass eigenstate, contributing a natural factor of 2-4x suppression without fine-tuning. Not sufficient alone, but relevant if the dilaton mass is near the boundary of detectability.",
        },
        {
            "name": "Open problem",
            "mechanism": "The framework does not derive the dilaton mass from first principles. Pure 6D EH+GB on S^2 gives V(sigma) = const (no mass term). Whether the mass is large enough for the heavy dilaton resolution depends on this unresolved question.",
            "plausibility": "status",
            "description": "The dilaton mass must come from an additional mechanism (flux stabilization, matter loops, or non-perturbative effects). If the mass is large (m_phi > 2e-3 eV), the heavy dilaton resolution dissolves the discrepancy entirely. If the mass is small (m_phi << 2e-3 eV), the full 3854x gap between tree-level and empirical screening amplitudes remains an open problem requiring either environment-dependent mass or moduli mixing effects.",
        },
    ]

    casimir_available = False
    try:
        from alpha_ladder_core.casimir_stabilization import summarize_casimir_stabilization
        cas = summarize_casimir_stabilization()
        casimir_available = True
        if cas.get("gaps_status", {}).get("gap3_resolved"):
            # Update heavy dilaton resolution with computed mass
            resolutions[0]["mechanism"] += f" Casimir-derived mass: {cas.get('dilaton_mass', {}).get('m_phi_eV', 'N/A')} eV."
            resolutions[0]["plausibility"] = "confirmed"
        else:
            # Update open problem resolution to note Casimir was attempted
            resolutions[3]["mechanism"] += " Casimir stabilization from the graviton tower was computed but does not produce a stable minimum (negative Casimir coefficient). The dilaton mass remains an input parameter."
    except ImportError:
        pass

    return {
        "alpha_tree": alpha_tree,
        "alpha_empirical": alpha_empirical,
        "ratio": ratio,
        "log10_ratio": log10_ratio,
        "resolutions": resolutions,
    }


def position_in_literature():
    """Compare this framework to 4 known frameworks in the literature."""
    comparisons = [
        {
            "framework": "Salam-Sezgin 6D gauged SUGRA",
            "similarities": [
                "6D starting point with S^2 compactification",
                "Single scalar (dilaton) in 4D effective theory",
                "Brans-Dicke-like gravity sector",
            ],
            "differences": [
                "Requires N=1 SUSY and monopole flux on S^2",
                "Gauge group constrained by anomaly cancellation",
                "Our framework uses pure gravity + GB, no SUSY",
            ],
            "key_reference": "Salam & Sezgin, Phys. Lett. B 147 (1984) 47",
        },
        {
            "framework": "KKLT (Kachru-Kallosh-Linde-Trivedi)",
            "similarities": [
                "Extra dimensions compactified to produce 4D gravity",
                "Moduli stabilization is a central concern",
                "Scalar potential from higher-dimensional physics",
            ],
            "differences": [
                "10D string theory starting point, not 6D",
                "Requires fluxes, anti-D3 branes, ~100 moduli",
                "Our framework: minimal 6D, 2 extra dims, single dilaton",
            ],
            "key_reference": "Kachru et al., Phys. Rev. D 68 (2003) 046005",
        },
        {
            "framework": "Nishino-Sezgin anomaly-free 6D",
            "similarities": [
                "6D supergravity with specific gauge group",
                "Green-Schwarz anomaly cancellation mechanism",
                "Compactification on smooth internal manifold",
            ],
            "differences": [
                "Requires specific gauge group for anomaly freedom",
                "Supersymmetric matter content prescribed",
                "Our framework is non-supersymmetric, pure gravity + GB",
            ],
            "key_reference": "Nishino & Sezgin, Nucl. Phys. B 278 (1986) 353",
        },
        {
            "framework": "ADD (Arkani-Hamed-Dimopoulos-Dvali) large extra dimensions",
            "similarities": [
                "Extra dimensions modify gravity at observable scales",
                "Hierarchy problem motivation",
                "Testable predictions for sub-mm gravity experiments",
            ],
            "differences": [
                "ADD uses flat extra dimensions, not curved S^2",
                "ADD targets the hierarchy problem, not alpha_g",
                "Our framework derives G from alpha via specific bridge coefficient",
            ],
            "key_reference": "Arkani-Hamed et al., Phys. Lett. B 429 (1998) 263",
        },
    ]

    return {"comparisons": comparisons}


def analyze_anomaly_status():
    """Assess the gravitational anomaly status of the 6D framework."""
    anomaly_types = [
        {
            "type": "Pure gravitational",
            "status": "absent",
            "description": "Pure gravity in 6D has no gravitational anomaly (the relevant index theorem vanishes for the graviton alone in d=6).",
        },
        {
            "type": "Mixed gauge-gravitational",
            "status": "not applicable",
            "description": "The minimal framework has no gauge fields beyond the metric and dilaton, so mixed anomalies do not arise.",
        },
        {
            "type": "Matter sector",
            "status": "unaddressed",
            "description": "A complete theory coupling to the Standard Model would require anomaly cancellation via Green-Schwarz mechanism or specific matter content. This is not addressed in the current framework.",
        },
    ]

    # Try to enrich with computed anomaly cancellation results
    try:
        from alpha_ladder_core.anomaly_cancellation import summarize_anomaly_cancellation
        ac = summarize_anomaly_cancellation()
        if ac.get("pure_gravity_safe"):
            anomaly_types[0]["status"] = "proven_safe"
            anomaly_types[0]["description"] = (
                "Proven: pure gravity in 6D is anomaly-free (graviton is non-chiral, "
                "I_8 coefficient = 0). Computed in anomaly_cancellation module."
            )
        if ac.get("group_scan", {}).get("any_contain_sm"):
            anomaly_types[2]["status"] = "constrained"
            anomaly_types[2]["description"] = (
                "Green-Schwarz anomaly cancellation requires embedding the SM in a "
                f"larger group. {ac['group_scan']['n_contain_sm']} anomaly-free groups "
                f"contain the SM; minimal is {ac['group_scan']['minimal_group']}. "
                "The G prediction is unaffected (gravity sector only)."
            )
    except ImportError:
        pass

    return {
        "anomaly_types": anomaly_types,
        "green_schwarz_required": False,
        "interpretation": (
            "The minimal 6D pure gravity + Gauss-Bonnet framework is anomaly-free "
            "because it contains no chiral fields. However, coupling to chiral matter "
            "(as needed for a complete Standard Model embedding) would require anomaly "
            "cancellation, likely via a Green-Schwarz mechanism with specific gauge groups. "
            "This is acknowledged as an open problem for the full theory."
        ),
    }


def summarize_theoretical_status(constants):
    """Top-level: screening discrepancy + literature + anomalies + honest assessment."""
    discrepancy = compute_screening_discrepancy(constants)
    literature = position_in_literature()
    anomalies = analyze_anomaly_status()

    strengths = [
        "Derives G from alpha with no free parameters (bridge = phi^2/2)",
        "Predicts G within 160 ppm of CODATA value",
        "6D origin provides a natural explanation for the 20+1 = 21 structure",
        "Massive dilaton passes solar system tests (Cassini bound)",
        "Vacuum polynomial uniquely selects n=2 extra dimensions",
    ]

    open_problems = [
        f"Tree-level screening amplitude is {discrepancy['ratio']:.0f}x larger than empirical value (if dilaton is light)",
        "No Green-Schwarz anomaly cancellation for matter sector coupling",
        "Moduli stabilization mechanism not derived from first principles",
        "Connection to string theory landscape unclear",
        "No prediction for cosmological constant",
    ]

    casimir_analysis = None
    try:
        from alpha_ladder_core.casimir_stabilization import summarize_casimir_stabilization
        casimir_analysis = summarize_casimir_stabilization()
    except ImportError:
        pass

    flux_available = False
    flux_result = None
    try:
        from alpha_ladder_core.flux_stabilization import summarize_flux_stabilization
        flux_result = summarize_flux_stabilization()
        flux_available = True
    except ImportError:
        pass

    if flux_available and flux_result and flux_result.get("gap_closure", {}).get("gap3_resolved"):
        open_problems[2] = (
            "Moduli stabilization: volume modulus stabilized by flux "
            "(Planck-scale mass); shape moduli by GB. Internal radius "
            "a_0 not determined."
        )
        open_problems[0] = (
            "Screening gap dissolved: Planck-mass dilaton decouples "
            "entirely. But dilaton phenomenology also lost."
        )

    # Try to integrate radius determination
    radius_result = None
    try:
        from alpha_ladder_core.radius_determination import summarize_radius_determination
        radius_result = summarize_radius_determination()
    except ImportError:
        pass

    # Try to integrate anomaly cancellation
    anomaly_result = None
    try:
        from alpha_ladder_core.anomaly_cancellation import summarize_anomaly_cancellation
        anomaly_result = summarize_anomaly_cancellation()
    except ImportError:
        pass

    # Try to integrate cosmological constant
    cc_result = None
    try:
        from alpha_ladder_core.cosmological_constant import summarize_cosmological_constant
        cc_result = summarize_cosmological_constant()
    except ImportError:
        pass

    # Update open problems with new module conclusions
    if radius_result:
        # Replace or update the a_0 open problem
        open_problems.append(
            "Internal radius a_0 undetermined: scaling symmetry proven "
            "(EH exponent = 0, GB topological for n=2). Flux+Casimir "
            "balance does not fix a_0."
        )

    if anomaly_result and anomaly_result.get("g_prediction_safe"):
        open_problems[1] = (
            "Anomaly cancellation: pure gravity safe; SM embedding requires "
            f"Green-Schwarz with group like {anomaly_result.get('group_scan', {}).get('minimal_group', 'E8xE8')}. "
            "G prediction unaffected."
        )

    if cc_result:
        open_problems[4] = (
            f"Cosmological constant: V_min ~ O(1) Planck, observed ~10^{{-122}}. "
            f"Discrepancy ~{cc_result.get('discrepancy_orders', 122)} orders. "
            f"Weinberg no-go applies. Universal unsolved problem."
        )

    return {
        "discrepancy": discrepancy,
        "literature": literature,
        "anomalies": anomalies,
        "strengths": strengths,
        "open_problems": open_problems,
        "casimir_analysis": casimir_analysis,
        "flux_analysis": flux_result,
        "radius_analysis": radius_result,
        "anomaly_analysis": anomaly_result,
        "cc_analysis": cc_result,
    }
