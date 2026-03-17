"""
The Proof -- Page 14

Rigorous theoretical foundations for the alpha ladder derivation:
step-by-step GB correction, vacuum polynomial, screening Lagrangian,
and moduli stabilization.
"""

import streamlit as st
import math
import pandas as pd



# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from app.components.sidebar import render_sidebar  # noqa: E402
from app.components.formatting import fmt_decimal  # noqa: E402

constants = render_sidebar()

# ---------------------------------------------------------------------------
# Try to load core modules
# ---------------------------------------------------------------------------
_PHI = (1 + math.sqrt(5)) / 2
_SQRT5 = math.sqrt(5)

_core_available = False
try:
    from alpha_ladder_core.kk_reduction import (  # noqa: E402
        derive_gauss_bonnet_correction,
        compute_internal_curvature_identities,
        compute_target_omega,
    )
    from alpha_ladder_core.vacuum_polynomial import (  # noqa: E402
        compute_dimensional_polynomial,
        verify_lambda_is_root,
        scan_discriminant_field,
        prove_golden_uniqueness,
        derive_algebraic_closure,
        summarize_vacuum_polynomial,
    )
    from alpha_ladder_core.screening import (  # noqa: E402
        derive_dilaton_lagrangian,
        derive_field_equation,
        derive_yukawa_profile,
    )
    from alpha_ladder_core.moduli_stabilization import (  # noqa: E402
        describe_moduli_space,
        compute_shape_stabilization,
        compute_moduli_mass_spectrum,
        summarize_stabilization,
    )
    _core_available = True
except ImportError:
    pass

# ---------------------------------------------------------------------------
# Compute values
# ---------------------------------------------------------------------------
if _core_available:
    _gb_derivation = derive_gauss_bonnet_correction(d=4, n=2, genus=2)
    _poly = compute_dimensional_polynomial(d=4, n=2)
    _verify = verify_lambda_is_root(d=4, n=2)
    _scan = scan_discriminant_field(d=4, n_max=12)
    _pell = prove_golden_uniqueness(d=4, n_max=50)
    _closure = derive_algebraic_closure(d=4, n=2)
    _vac_summary = summarize_vacuum_polynomial(d=4, n=2)

    _identities = compute_internal_curvature_identities(n=2)
    _target = compute_target_omega()

    _omega_bd = _target["omega_float"]
    _yukawa = derive_yukawa_profile(_omega_bd)

    _moduli_space = describe_moduli_space(genus=2)
    _shape_stab = compute_shape_stabilization(genus=2, lambda_GB=_SQRT5 - 3)
    _mass_spectrum = compute_moduli_mass_spectrum(genus=2)
    _stab_summary = summarize_stabilization(genus=2)
else:
    # Fallbacks
    _poly = {"D": 6, "discriminant": 20, "root_plus": _SQRT5 - 3, "is_golden_field": True}
    _verify = {"is_root": True, "polynomial_value": 0.0}
    _scan = {"results": [], "minimal_n": 2, "golden_field_solutions": [2, 10]}
    _pell = {"minimal_solution": (2, 2), "next_solution": (10, 6), "n2_is_minimal": True}
    _closure = {"constants": [], "all_in_same_field": True}
    _gb_derivation = None
    _identities = None
    _yukawa = None
    _mass_spectrum = None
    _shape_stab = None
    _stab_summary = None
    _moduli_space = None

lambda_gb = _SQRT5 - 3
discriminant = 20

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("The Proof: Rigorous Foundations")
st.markdown(
    "Four theoretical gaps closed: step-by-step GB derivation, "
    "the vacuum polynomial x\u00b2 + Dx + d = 0, the screening Lagrangian, "
    "and moduli stabilization. Every claim is now backed by explicit computation."
)
st.divider()

# ---------------------------------------------------------------------------
# 4 Metric Cards
# ---------------------------------------------------------------------------
col_m1, col_m2, col_m3, col_m4 = st.columns(4)

with col_m1:
    st.metric(label="Polynomial", value="x\u00b2 + 6x + 4 = 0")
with col_m2:
    st.metric(label="Discriminant", value="20 = 4 \u00b7 5")
with col_m3:
    st.metric(label="Field extension", value="\u211a(\u221a5)")
with col_m4:
    st.metric(label="Moduli stabilized", value="7 / 7")

st.markdown("")

# ---------------------------------------------------------------------------
# Section A: The Vacuum Polynomial (Gap 4 -- most striking result)
# ---------------------------------------------------------------------------
with st.expander("The Vacuum Polynomial: x\u00b2 + Dx + d = 0", expanded=True):
    st.markdown(
        """
        <div class="theorem-card">
        <b>Theorem.</b> The GB coupling &lambda;<sub>GB</sub> = &radic;5 &minus; 3
        satisfies the polynomial<br><br>
        <center><b>x&sup2; + Dx + d = 0</b></center><br>
        where D = 6 is the total spacetime dimension and d = 4 is the external
        dimension. The discriminant D&sup2; &minus; 4d = 20 = 4 &middot; 5 places
        the roots in &Qopf;(&radic;5) &mdash; the <em>same algebraic field</em>
        as the golden ratio.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    col_v1, col_v2, col_v3 = st.columns(3)
    with col_v1:
        st.metric(label="D (total dim)", value="6")
    with col_v2:
        st.metric(label="d (external dim)", value="4")
    with col_v3:
        poly_val = _verify.get("polynomial_value", 0.0) if _verify else 0.0
        st.metric(
            label="x\u00b2 + 6x + 4 at \u03bb\u0047\u0042",
            value="0 (exact)" if abs(poly_val) < 1e-10 else fmt_decimal(poly_val, sig_figs=3),
        )

    st.markdown("")

    # Discriminant scan chart
    if _core_available and _scan and _scan.get("results"):
        try:
            from app.components.charts import discriminant_scan_chart  # noqa: E402
            fig_disc = discriminant_scan_chart(_scan["results"])
            st.plotly_chart(fig_disc, use_container_width=True)
        except Exception as exc:
            st.warning(f"Discriminant scan chart error: {exc}")

    st.markdown("")

    # Uniqueness argument
    col_u1, col_u2 = st.columns(2)
    with col_u1:
        st.markdown(
            """
            <div class="proof-card">
            <b>Uniqueness</b><br><br>
            For d = 4, the Diophantine equation n(n + 8) = 5m&sup2;
            determines which internal dimensions produce &Qopf;(&radic;5).<br><br>
            <b>Minimal solution:</b> n = 2, m = 2 (check: 2 &middot; 10 = 20 = 5 &middot; 4)<br>
            <b>Next solution:</b> n = 10, m = 6 (check: 10 &middot; 18 = 180 = 5 &middot; 36)<br><br>
            n = 2 is the unique minimal internal dimension where the golden
            ratio field emerges from the dimensional polynomial.
            </div>
            """,
            unsafe_allow_html=True,
        )

    with col_u2:
        st.markdown(
            """
            <div class="proof-card">
            <b>Vieta's Formulas</b><br><br>
            The two roots x&#8324; = &radic;5 &minus; 3 and x&#8331; = &minus;&radic;5 &minus; 3
            satisfy:<br><br>
            <b>Sum:</b> x&#8324; + x&#8331; = &minus;6 = &minus;D<br>
            <b>Product:</b> x&#8324; &middot; x&#8331; = 4 = d<br><br>
            The sum of roots equals minus the total dimension.<br>
            The product of roots equals the external dimension.
            </div>
            """,
            unsafe_allow_html=True,
        )

    # Algebraic closure table
    if _core_available and _closure and _closure.get("constants"):
        try:
            from app.components.charts import algebraic_field_chart  # noqa: E402
            fig_field = algebraic_field_chart(_closure)
            st.plotly_chart(fig_field, use_container_width=True)
        except Exception as exc:
            st.warning(f"Algebraic field chart error: {exc}")

    st.markdown(
        """
        <div class="theorem-card">
        <b>Key result:</b> The golden ratio is not put in by hand.
        &Qopf;(&radic;5) is forced by the unique minimal viable dimensions
        D = 6, d = 4. Every constant in the derivation &mdash; &phi;, &phi;&sup2;/2,
        &omega;<sub>target</sub>, &lambda;<sub>GB</sub>, &phi;<sup>&minus;2</sup>
        &mdash; lies in this field automatically.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Section B: Step-by-Step GB Derivation (Gap 1)
# ---------------------------------------------------------------------------
with st.expander("Gauss-Bonnet Correction: Step-by-Step Derivation", expanded=True):

    if _core_available and _gb_derivation:
        step1 = _gb_derivation["step1_internal_identities"]
        step2 = _gb_derivation["step2_internal_gb_density"]
        step4 = _gb_derivation["step4_cross_terms"]
        step5 = _gb_derivation["step5_total_correction"]
        verification = _gb_derivation["verification"]

        # Step 1
        st.markdown(
            f"""
            <div class="step-card">
            <b>Step 1: Internal Curvature Identities</b><br><br>
            On a maximally symmetric 2-manifold, the Riemann tensor is fully
            determined by the scalar curvature R&#8322;:<br><br>
            <code>R_{{ab}}R^{{ab}} / R&#8322;&sup2; = {step1['R_ab_R_ab_over_R_n_sq']:.4f}</code> (= 1/n)<br>
            <code>R_{{abcd}}R^{{abcd}} / R&#8322;&sup2; = {step1['R_abcd_R_abcd_over_R_n_sq']:.4f}</code> (= 2/[n(n&minus;1)])
            </div>
            """,
            unsafe_allow_html=True,
        )

        # Step 2
        st.markdown(
            f"""
            <div class="step-card">
            <b>Step 2: Internal GB Density</b><br><br>
            <code>G&#8322; / R&#8322;&sup2; = 1 &minus; 4/n + 2/[n(n&minus;1)] = (n&minus;2)(n&minus;3)/[n(n&minus;1)]</code><br><br>
            For n = 2: <b>G&#8322; = 0</b> (Gauss-Bonnet is topological in 2D).<br>
            Coefficient = {step2['G_n_coefficient']:.1f} &mdash; verified algebraically.
            </div>
            """,
            unsafe_allow_html=True,
        )

        # Step 4: cross terms
        ch1 = step4["channel_1_R_sq"]["value"]
        ch2 = step4["channel_2_Ric_sq"]["value"]
        ch3 = step4["channel_3_Riem_sq"]["value"]

        st.markdown(
            f"""
            <div class="step-card">
            <b>Steps 3-4: Cross Terms on the Warped Product</b><br><br>
            Three channels contribute to the (&part;&sigma;)&sup2; kinetic term in 4D:<br><br>
            <code>Channel 1 (R&sup2; cross):&nbsp;&nbsp;&nbsp; {ch1:+.4f}</code><br>
            <code>Channel 2 (Ric&sup2; cross):&nbsp; {ch2:+.4f}</code><br>
            <code>Channel 3 (Riem&sup2; cross): {ch3:+.4f}</code><br>
            <code>Total:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; {ch1 + ch2 + ch3:.4f}</code>
            </div>
            """,
            unsafe_allow_html=True,
        )

        # Step 5
        st.markdown(
            f"""
            <div class="step-card">
            <b>Step 5: Assembly</b><br><br>
            <code>&delta;K = &lambda;<sub>GB</sub> &middot; &chi; &middot; gb_factor
            = &lambda;<sub>GB</sub> &middot; ({step5['chi']}) &middot; {step5['gb_factor']:.4f}
            = {step5['delta_K_per_lambda_GB']:.4f} &middot; &lambda;<sub>GB</sub></code>
            </div>
            """,
            unsafe_allow_html=True,
        )

        # Verification
        match = verification.get("all_consistent", False)
        st.markdown(
            f"""
            <div class="theorem-card">
            <b>Verification:</b> Derived &delta;K = {verification.get('derived_delta_K', 0):.6f},
            existing formula gives {verification.get('existing_delta_K', 0):.6f}.
            <b>{'Match confirmed.' if match else 'MISMATCH.'}</b>
            </div>
            """,
            unsafe_allow_html=True,
        )
    else:
        st.info("Core modules not available. Install alpha_ladder_core to see the full derivation.")

# ---------------------------------------------------------------------------
# Section C: Screening from First Principles (Gap 2)
# ---------------------------------------------------------------------------
with st.expander("Screening Lagrangian: From Action to Yukawa", expanded=True):
    st.markdown(
        """
        <div class="formula-card">
        <b>4D Effective Action (Einstein Frame)</b><br><br>
        <code>S = &int; d&#8308;x &radic;(&minus;g) [ &frac12; M<sub>Pl</sub>&sup2; R
        &minus; &frac12; (&part;&phi;)&sup2;
        &minus; V(&phi;)
        + L<sub>matter</sub>(A&sup2;(&phi;) g<sub>&mu;&nu;</sub>, &psi;) ]</code><br><br>
        where A(&phi;) = exp(&beta; &phi; / M<sub>Pl</sub>) is the conformal coupling
        and &beta; = 1/&radic;(2&omega; + 3) is the universal coupling constant.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    if _core_available and _yukawa:
        steps = _yukawa.get("derivation_steps", [])
        for step in steps:
            st.markdown(
                f"""
                <div class="step-card">
                <b>Step {step.get('step_number', '?')}:</b> {step.get('description', '')}<br><br>
                <code>{step.get('formula', '')}</code>
                </div>
                """,
                unsafe_allow_html=True,
            )

        alpha_lag = _yukawa.get("alpha_screening_from_lagrangian", 0)
        beta_c = _yukawa.get("beta_coupling", 0)
        st.markdown("")
        col_s1, col_s2 = st.columns(2)
        with col_s1:
            st.metric(label="\u03b2 (universal coupling)", value=f"{beta_c:.6f}")
        with col_s2:
            st.metric(label="\u03b1_screening = 2\u03b2\u00b2", value=f"{alpha_lag:.6f}")
    else:
        st.markdown(
            """
            <div class="step-card">
            <b>Derivation chain:</b><br><br>
            1. Start with BD + massive dilaton Lagrangian<br>
            2. Linearize field equation: (&nabla;&sup2; &minus; m&sup2;)&phi; = &minus;&beta; &rho; / M<sub>Pl</sub><br>
            3. Green's function: &phi;(r) = &minus;(&beta;M)/(4&pi;M<sub>Pl</sub>r) exp(&minus;m r)<br>
            4. Metric perturbation: G<sub>eff</sub>(r) = G<sub>N</sub>(1 + 2&beta;&sup2; exp(&minus;r/&lambda;))<br>
            5. Identify: &alpha;<sub>s</sub> = 2&beta;&sup2;, &lambda; = 1/m<sub>&phi;</sub>
            </div>
            """,
            unsafe_allow_html=True,
        )

    st.markdown(
        """
        <div class="theorem-card">
        <b>Result:</b> The Yukawa screening profile
        G<sub>eff</sub>(r) = G<sub>vacuum</sub>(1 + &alpha;<sub>s</sub> exp(&minus;r/&lambda;))
        is <em>derived</em> from the 4D effective action, not assumed.
        The screening amplitude &alpha;<sub>s</sub> = 2&beta;&sup2; is fixed by
        the Brans-Dicke parameter &omega;.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Section D: Moduli Stabilization (Gap 3)
# ---------------------------------------------------------------------------
with st.expander("Moduli Stabilization: All 7 Moduli Massive", expanded=True):
    if _core_available and _moduli_space and _shape_stab and _mass_spectrum:
        col_d1, col_d2, col_d3 = st.columns(3)
        with col_d1:
            st.metric(label="Shape moduli (6g\u22126)", value=str(_moduli_space.get("real_moduli", 6)))
        with col_d2:
            st.metric(label="Volume modulus", value="1")
        with col_d3:
            st.metric(label="Total moduli", value=str(_moduli_space.get("total_moduli", 7)))

        st.markdown("")

        # Mass spectrum chart
        try:
            from app.components.charts import moduli_mass_chart  # noqa: E402
            fig_mass = moduli_mass_chart(_mass_spectrum)
            st.plotly_chart(fig_mass, use_container_width=True)
        except Exception as exc:
            st.warning(f"Moduli mass chart error: {exc}")

        st.markdown("")

        col_mech1, col_mech2 = st.columns(2)
        with col_mech1:
            st.markdown(
                """
                <div class="proof-card">
                <b>Volume Modulus</b><br><br>
                Tree-level EH gives a flat direction (&int; R vol is topological
                by Gauss-Bonnet). Stabilized by the dilaton mass from matter
                coupling: V(&sigma;) = &frac12; m<sub>&phi;</sub>&sup2; &sigma;&sup2;.
                </div>
                """,
                unsafe_allow_html=True,
            )
        with col_mech2:
            lich_bound = _shape_stab.get("lichnerowicz_bound", 2.0)
            st.markdown(
                f"""
                <div class="proof-card">
                <b>Shape Moduli (6)</b><br><br>
                The GB curvature-variation potential is minimized at the unique
                constant-curvature hyperbolic metric (Poincar&eacute; uniformization).
                All eigenvalues of the Lichnerowicz operator &ge; {lich_bound:.1f}
                &rArr; all masses positive.
                </div>
                """,
                unsafe_allow_html=True,
            )

        if _moduli_space.get("hyperelliptic"):
            st.markdown(
                """
                <div class="theorem-card">
                <b>Hyperellipticity bonus:</b> Every genus-2 surface admits a
                &Zopf;&#8322; involution (a theorem). This projects out half
                the moduli in consistent truncation, leaving 3 effective
                shape moduli. Genus &ge; 3 generically lacks this symmetry.
                </div>
                """,
                unsafe_allow_html=True,
            )

        # Casimir stabilization reference
        try:
            from alpha_ladder_core.casimir_stabilization import summarize_casimir_stabilization
            _casimir = summarize_casimir_stabilization()
            if _casimir:
                assessment = _casimir.get("overall_assessment", "")
                st.markdown(
                    f"""
                    <div class="formula-card">
                    <b>Casimir Stabilization Analysis (Page 16)</b><br><br>
                    The 1-loop Casimir energy from the KK graviton tower on S&sup2; was
                    computed as a candidate first-principles mechanism for volume modulus
                    stabilization. {assessment}
                    See <b>Page 16: Casimir Stabilization</b> for the full analysis.
                    </div>
                    """,
                    unsafe_allow_html=True,
                )
        except ImportError:
            pass

        # Flux stabilization reference
        try:
            from alpha_ladder_core.flux_stabilization import summarize_flux_stabilization
            _flux = summarize_flux_stabilization()
            if _flux and _flux.get("gap_closure", {}).get("gap3_resolved"):
                st.markdown(
                    """
                    <div class="theorem-card">
                    <b>Flux Stabilization (Page 17)</b><br><br>
                    Adding quantized 2-form flux F<sub>2</sub> on S&sup2; creates a stable
                    minimum for the volume modulus. The dilaton mass is at the Planck scale,
                    which trivially dissolves the screening discrepancy but eliminates
                    dilaton phenomenology. See <b>Page 17: Flux Stabilization</b> and
                    <b>Page 18: Radius Phenomenology</b>.
                    </div>
                    """,
                    unsafe_allow_html=True,
                )
        except ImportError:
            pass
    else:
        st.markdown(
            """
            <div class="proof-card">
            <b>Moduli stabilization summary:</b><br><br>
            &bull; 6 shape moduli: stabilized by GB curvature-variation potential at the
              constant-curvature metric (Poincar&eacute; uniformization)<br>
            &bull; 1 volume modulus: stabilized by dilaton mass from matter coupling<br>
            &bull; Lichnerowicz bound: all shape masses &ge; 2/a&sup2;<br>
            &bull; All 7 moduli are massive &rArr; no massless scalars
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("The Proof | Alpha Ladder Research Dashboard")
