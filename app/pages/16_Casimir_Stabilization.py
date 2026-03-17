"""
Casimir Stabilization -- Page 16

Analysis of 1-loop Casimir energy from the KK graviton tower on S^2
as a candidate mechanism for volume modulus stabilization.
"""

import streamlit as st
import math



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
    from alpha_ladder_core.casimir_stabilization import (  # noqa: E402
        compute_kk_spectrum_s2,
        compute_spectral_zeta_s2,
        compute_casimir_energy_s2,
        compute_effective_potential,
        find_casimir_minimum,
        compute_dilaton_mass_casimir,
        summarize_casimir_stabilization,
        compute_fermion_spectral_zeta_s2,
        compute_vector_spectral_zeta_s2,
        compute_matter_casimir_coefficient,
        scan_anomaly_free_matter_casimir,
    )
    _core_available = True
except ImportError:
    pass

# ---------------------------------------------------------------------------
# Compute values (guarded)
# ---------------------------------------------------------------------------
_spectrum = None
_zeta = None
_casimir = None
_potential = None
_minimum = None
_mass = None
_summary = None

if _core_available:
    try:
        _spectrum = compute_kk_spectrum_s2(l_max=20)
    except Exception:
        pass
    try:
        _zeta = compute_spectral_zeta_s2()
    except Exception:
        pass
    try:
        _casimir = compute_casimir_energy_s2()
    except Exception:
        pass
    try:
        _potential = compute_effective_potential(sigma=0.0)
    except Exception:
        pass
    try:
        _minimum = find_casimir_minimum()
    except Exception:
        pass
    try:
        _mass = compute_dilaton_mass_casimir()
    except Exception:
        pass
    try:
        _summary = summarize_casimir_stabilization()
    except Exception:
        pass

_matter_scan = None
if _core_available:
    try:
        _matter_scan = scan_anomaly_free_matter_casimir()
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Casimir Stabilization: 1-Loop Analysis")
st.markdown(
    "Can the 1-loop Casimir energy from the Kaluza-Klein graviton tower on "
    "S\u00b2 stabilize the volume modulus? This page presents the honest "
    "computation and its negative result -- the Casimir coefficient from the "
    "graviton tower alone does not produce a stable minimum."
)
st.divider()

# ---------------------------------------------------------------------------
# 3 Metric Cards
# ---------------------------------------------------------------------------
col_m1, col_m2, col_m3 = st.columns(3)

with col_m1:
    total_modes = _spectrum.get("total_modes", "N/A") if _spectrum else "N/A"
    st.metric(label="Total KK modes (l ≤ 20)", value=str(total_modes))
with col_m2:
    st.metric(label="Bosonic dof", value="9")
with col_m3:
    st.metric(label="Field decomposition", value="2 + 4 + 3")

st.markdown("")

# ---------------------------------------------------------------------------
# Section A: The Volume Modulus Problem
# ---------------------------------------------------------------------------
with st.expander("A. The Volume Modulus Problem", expanded=True):
    st.markdown(
        """
        <div class="theorem-card">
        <b>Theorem (Gauss-Bonnet).</b> On a closed 2-surface &Sigma;, the
        integral of the scalar curvature is topological:<br><br>
        <center><b>&int;<sub>&Sigma;</sub> R<sub>2</sub> vol = 4&pi; &chi;(&Sigma;)</b></center><br>
        where &chi; is the Euler characteristic. For genus g: &chi; = 2 &minus; 2g.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="step-card">
        <b>Why this matters:</b> The classical Einstein-Hilbert action on
        the internal space reduces to a topological invariant. The volume
        modulus &sigma; (breathing mode of the extra dimensions) has <em>no
        potential</em> at tree level &mdash; V<sub>classical</sub>(&sigma;)
        is flat.<br><br>
        To stabilize &sigma;, we need quantum corrections. The leading
        candidate is the 1-loop Casimir energy from the tower of
        Kaluza-Klein modes propagating on the internal manifold.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    col_a1, col_a2 = st.columns(2)
    with col_a1:
        st.markdown(
            """
            <div class="proof-card">
            <b>Tree-level potential</b><br><br>
            V<sub>tree</sub>(&sigma;) = &Lambda;<sub>6</sub> e<sup>2&sigma;</sup>
            &minus; &frac12; R<sub>2</sub> e<sup>0</sup><br><br>
            But R<sub>2</sub> integrates to 4&pi;&chi; (topological), so the
            &sigma;-dependence drops out.  The volume modulus is a flat
            direction of the classical potential.
            </div>
            """,
            unsafe_allow_html=True,
        )
    with col_a2:
        st.markdown(
            """
            <div class="proof-card">
            <b>Quantum correction strategy</b><br><br>
            V<sub>1-loop</sub>(&sigma;) = V<sub>curv</sub>(&sigma;)
            + V<sub>Cas</sub>(&sigma;)<br><br>
            The curvature contribution scales as e<sup>&minus;2&sigma;</sup>
            while the Casimir contribution scales as e<sup>&minus;4&sigma;</sup>.
            A stable minimum requires these to balance with the correct signs.
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section B: KK Spectrum on S^2
# ---------------------------------------------------------------------------
with st.expander("B. KK Spectrum on S\u00b2", expanded=True):
    st.markdown(
        """
        <div class="step-card">
        <b>Field content from 6D gravity on M<sub>4</sub> &times; S&sup2;</b><br><br>
        The 6D graviton decomposes into 4D fields at each KK level l:<br><br>
        &bull; <b>Spin-2 (metric):</b> 2 dof &mdash; the 4D graviton tower<br>
        &bull; <b>Spin-1 (gravi-photon):</b> 4 dof &mdash; KK gauge bosons<br>
        &bull; <b>Spin-0 (scalar):</b> 3 dof &mdash; breathing + shape modes<br><br>
        Total: <b>9 bosonic degrees of freedom</b> per KK level.
        Each level l has degeneracy (2l + 1) and mass&sup2;
        = l(l + 1) / a&sup2;.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    if _core_available and _spectrum:
        # KK spectrum chart
        try:
            from app.components.charts import kk_spectrum_chart  # noqa: E402
            fig_kk = kk_spectrum_chart(_spectrum)
            st.plotly_chart(fig_kk, use_container_width=True)
        except ImportError:
            st.info("Chart function kk_spectrum_chart not yet available in charts module.")
        except Exception as exc:
            st.warning(f"KK spectrum chart error: {exc}")

        st.markdown("")

        # Spectrum summary metrics
        col_b1, col_b2, col_b3 = st.columns(3)
        with col_b1:
            l_max = _spectrum.get("l_max", 20)
            st.metric(label="l max", value=str(l_max))
        with col_b2:
            total_degen = _spectrum.get("total_degeneracy", "N/A")
            st.metric(label="Total degeneracy", value=str(total_degen))
        with col_b3:
            st.metric(label="Mass gap (l=1)", value="2 / a₀²")
    else:
        st.markdown(
            """
            <div class="formula-card">
            <b>KK mass spectrum on S&sup2;</b><br><br>
            <code>m&sup2;(l) = l(l + 1) / a&sup2;,&nbsp;&nbsp;
            degeneracy = 2l + 1,&nbsp;&nbsp; l = 0, 1, 2, ...</code><br><br>
            Total modes up to l<sub>max</sub> = 20:
            &sum; (2l+1) = (l<sub>max</sub>+1)&sup2; = 441 per field type.
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section C: Effective Potential
# ---------------------------------------------------------------------------
with st.expander("C. Effective Potential", expanded=True):

    st.markdown(
        """
        <div class="formula-card">
        <b>Curvature contribution</b><br><br>
        <code>V<sub>curv</sub>(&sigma;) = A &middot; e<sup>&minus;2&sigma;</sup></code><br><br>
        where A encodes the internal curvature R<sub>2</sub> = 2/a&sup2;
        (positive for S&sup2;). This provides a repulsive (positive)
        contribution that grows as the extra dimensions shrink.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="formula-card">
        <b>Casimir contribution</b><br><br>
        <code>V<sub>Cas</sub>(&sigma;) = B &middot; e<sup>&minus;4&sigma;</sup></code><br><br>
        where B is the Casimir coefficient computed from the spectral zeta
        function &zeta;(s) of the Laplacian on S&sup2;, summed over all 9
        bosonic dof with appropriate signs.<br><br>
        <b>The sign of B determines stability.</b> B &gt; 0 (attractive)
        could balance V<sub>curv</sub>; B &lt; 0 (repulsive) cannot.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    if _core_available and _potential:
        # Potential chart
        try:
            from app.components.charts import casimir_potential_chart  # noqa: E402
            fig_pot = casimir_potential_chart(_potential)
            st.plotly_chart(fig_pot, use_container_width=True)
        except ImportError:
            st.info("Chart function casimir_potential_chart not yet available in charts module.")
        except Exception as exc:
            st.warning(f"Casimir potential chart error: {exc}")

        st.markdown("")

        # Coefficient values
        col_c1, col_c2 = st.columns(2)
        with col_c1:
            A_val = _potential.get("A_coefficient", None)
            st.metric(
                label="A (curvature coefficient)",
                value=fmt_decimal(A_val, sig_figs=6) if A_val is not None else "N/A",
            )
        with col_c2:
            B_val = _potential.get("B_coefficient", None)
            st.metric(
                label="B (Casimir coefficient)",
                value=fmt_decimal(B_val, sig_figs=6) if B_val is not None else "N/A",
            )

        # Genus-2 comparison if available
        genus2_data = _potential.get("genus2_data")
        if genus2_data:
            st.markdown("")
            st.markdown(
                f"""
                <div class="step-card">
                <b>Genus-2 comparison:</b> On a genus-2 surface (&chi; = &minus;2),
                the curvature contribution flips sign (R &lt; 0 for hyperbolic
                metrics). The Casimir coefficient becomes
                B<sub>g=2</sub> = {fmt_decimal(genus2_data.get('B_coefficient'), sig_figs=6)}.
                </div>
                """,
                unsafe_allow_html=True,
            )
    else:
        st.markdown(
            """
            <div class="step-card">
            <b>Effective potential structure:</b><br><br>
            V<sub>eff</sub>(&sigma;) = A e<sup>&minus;2&sigma;</sup>
            + B e<sup>&minus;4&sigma;</sup><br><br>
            For a minimum at &sigma;<sub>0</sub>, we need:<br>
            &bull; dV/d&sigma; = 0 &rArr; &sigma;<sub>0</sub> =
            &frac12; ln(&minus;2B/A)<br>
            &bull; d&sup2;V/d&sigma;&sup2; &gt; 0 at &sigma;<sub>0</sub><br><br>
            This requires A &gt; 0 and B &lt; 0 (or equivalently, B/A &lt; 0
            with the appropriate sign convention). The graviton tower on
            S&sup2; produces the wrong sign for B.
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section D: Dilaton Mass Prediction -- THE KEY RESULT
# ---------------------------------------------------------------------------
with st.expander("D. Dilaton Mass: The Negative Result", expanded=True):

    # The central negative result
    st.markdown(
        """
        <div class="warning-card">
        <b>Key Result: No Stable Minimum</b><br><br>
        The Casimir coefficient B from the graviton tower on S&sup2; is
        <b>negative</b> (repulsive). Combined with the positive curvature
        contribution A, both terms push the volume modulus toward
        decompactification (&sigma; &rarr; +&infin;).<br><br>
        <b>No stable minimum exists for the volume modulus from Casimir
        stabilization alone.</b><br><br>
        The dilaton mass <em>cannot</em> be derived from this mechanism.
        Gap #3 (volume modulus stabilization) remains open for the Casimir
        channel.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    # What WAS computed
    st.markdown("#### What was computed")

    if _core_available and _casimir:
        col_d1, col_d2, col_d3 = st.columns(3)
        with col_d1:
            zeta_val = _casimir.get("spectral_zeta_value", None)
            st.metric(
                label="Spectral zeta value",
                value=fmt_decimal(zeta_val, sig_figs=6) if zeta_val is not None else "N/A",
            )
        with col_d2:
            casimir_coeff = _casimir.get("casimir_coefficient", None)
            st.metric(
                label="Casimir coefficient B",
                value=fmt_decimal(casimir_coeff, sig_figs=6) if casimir_coeff is not None else "N/A",
            )
        with col_d3:
            sign_str = _casimir.get("sign", "negative")
            st.metric(label="Sign of B", value=sign_str.upper())
    elif _core_available and _zeta:
        col_d1, col_d2 = st.columns(2)
        with col_d1:
            zeta_val = _zeta.get("value", None)
            st.metric(
                label="Spectral zeta value",
                value=fmt_decimal(zeta_val, sig_figs=6) if zeta_val is not None else "N/A",
            )
        with col_d2:
            st.metric(label="Sign of B", value="NEGATIVE")
    else:
        st.markdown(
            """
            <div class="step-card">
            <b>Computed quantities:</b><br><br>
            &bull; Spectral zeta function &zeta;<sub>S&sup2;</sub>(s)
            = &sum;<sub>l=0</sub><sup>&infin;</sup> (2l+1)
            [l(l+1)]<sup>&minus;s</sup> (regularized)<br>
            &bull; Casimir energy E<sub>Cas</sub> = &frac12; &zeta;(&minus;1/2)
            summed over 9 bosonic dof<br>
            &bull; The coefficient B &lt; 0 for the graviton tower on
            S&sup2;
            </div>
            """,
            unsafe_allow_html=True,
        )

    st.markdown("")

    # What REMAINS OPEN
    st.markdown("#### What remains open")
    st.markdown(
        """
        <div class="warning-card">
        <b>The dilaton mass is undetermined.</b><br><br>
        Since V<sub>eff</sub>(&sigma;) has no minimum, we cannot extract
        m<sub>&phi;</sub>&sup2; = d&sup2;V/d&sigma;&sup2;|<sub>&sigma;=&sigma;<sub>0</sub></sub>.
        The volume modulus remains a runaway direction. The mass spectrum
        computed on Page 14 (moduli stabilization via curvature-variation)
        remains the only available stabilization mechanism.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    # Possible resolutions
    st.markdown("#### Possible resolutions")
    st.markdown(
        """
        <div class="proof-card">
        <b>Candidate mechanisms that could produce a stable minimum:</b><br><br>
        1. <b>Flux stabilization:</b> Adding form-field fluxes threading
        the internal manifold generates a positive potential
        V<sub>flux</sub> ~ e<sup>&minus;6&sigma;</sup> that could
        overwhelm the Casimir repulsion.<br><br>
        2. <b>Higher-genus topology:</b> On a genus-2 surface, the negative
        curvature changes the sign structure of V<sub>curv</sub>.
        Combined with the Casimir term, this may admit a minimum.<br><br>
        3. <b>Fermionic Casimir contributions:</b> The graviton tower
        is purely bosonic. Including gravitino or bulk fermion modes
        flips the sign of the Casimir energy (fermions contribute with
        opposite sign), potentially stabilizing the modulus.<br><br>
        4. <b>Non-perturbative effects:</b> Instantons or gaugino
        condensation in the internal space could generate an
        exponentially suppressed but stabilizing superpotential.<br><br>
        5. <b>Mixed mechanism:</b> The curvature-variation potential from
        Page 14 (which <em>does</em> stabilize shape moduli) combined
        with a small Casimir correction may produce an effective
        minimum for the volume modulus at a shifted location.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Section E: Impact on Open Problems
# ---------------------------------------------------------------------------
with st.expander("E. Impact on Open Problems", expanded=True):

    if _core_available and _summary:
        gap_status = _summary.get("gap_status", {})

        gap3 = gap_status.get("gap3", {})
        gap1 = gap_status.get("gap1", {})

        gap3_status = gap3.get("status", "open")
        gap1_status = gap1.get("status", "dependent")

        col_e1, col_e2 = st.columns(2)
        with col_e1:
            st.metric(
                label="Gap #3: Volume stabilization",
                value=gap3_status.upper(),
            )
        with col_e2:
            st.metric(
                label="Gap #1: Dilaton mass",
                value=gap1_status.upper(),
            )

        st.markdown("")

        assessment = _summary.get("assessment", "")
        if assessment:
            st.markdown(
                f"""
                <div class="theorem-card">
                <b>Assessment:</b> {assessment}
                </div>
                """,
                unsafe_allow_html=True,
            )
    else:
        col_e1, col_e2 = st.columns(2)
        with col_e1:
            st.metric(label="Gap #3: Volume stabilization", value="OPEN")
        with col_e2:
            st.metric(label="Gap #1: Dilaton mass", value="DEPENDENT")

    st.markdown("")

    st.markdown(
        """
        <div class="theorem-card">
        <b>Honest assessment:</b> The Casimir energy from the pure graviton
        KK tower on S&sup2; does not stabilize the volume modulus. This is
        a genuine negative result, not a failure of computation &mdash; the
        physics simply does not cooperate.<br><br>
        This rules out the simplest Casimir mechanism and narrows the space
        of viable stabilization strategies. The curvature-variation
        stabilization from Page 14 (Gap #3) remains the primary mechanism
        for shape moduli, but the volume modulus requires additional
        ingredients beyond the pure graviton Casimir effect.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="step-card">
        <b>What this result contributes:</b><br><br>
        &bull; Explicit computation of the spectral zeta function on
        S&sup2; for the full graviton tower<br>
        &bull; Proof that B &lt; 0 for bosonic modes &mdash; a
        no-go for Casimir stabilization alone<br>
        &bull; Enumeration of alternative mechanisms that could
        circumvent the no-go<br>
        &bull; Clarification that Gap #3 from Page 14 remains open
        for the volume direction<br><br>
        The negative result is itself a contribution: it rules out one
        candidate mechanism and sharpens the requirements on any
        successful stabilization proposal.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Section F: Matter Loop Corrections
# ---------------------------------------------------------------------------
with st.expander("F. Matter Loop Corrections", expanded=True):

    st.markdown(
        """
        <div class="step-card">
        <b>Can matter fields flip the Casimir sign?</b><br><br>
        The pure graviton tower gives a negative Casimir coefficient
        (&zeta;<sub>S&sup2;</sub>(&minus;1/2) &lt; 0). Fermions contribute
        with <em>opposite sign</em> to the Casimir energy. If the anomaly-free
        matter content includes enough fermionic degrees of freedom, the total
        coefficient A could flip positive, enabling stabilization.<br><br>
        <b>Key result:</b> The fermion spectral zeta &zeta;<sub>F</sub>(&minus;1/2) = 0
        exactly (because B<sub>3</sub>(1/2) = 0), so individual fermion species
        contribute <em>zero</em> at leading order. The sign flip depends entirely
        on the vector boson spectrum.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    if _core_available and _matter_scan:
        results = _matter_scan.get("results", [])
        any_flip = _matter_scan.get("any_sign_flip", False)

        if results:
            st.markdown("#### Anomaly-Free Gauge Groups")

            # Build table data
            table_rows = []
            for r in results:
                cas = r.get("casimir", {})
                table_rows.append({
                    "Gauge Group": r.get("group", ""),
                    "n_scalars": r.get("n_scalars_total", 0),
                    "n_fermions": r.get("n_fermions_total", 0),
                    "n_vectors": r.get("n_vectors_total", 0),
                    "A_total": f"{cas.get('A_total', 0):.6f}",
                    "Sign Flipped": "Yes" if cas.get("sign_flipped", False) else "No",
                })

            st.table(table_rows)

            st.markdown("")

        # Assessment card
        if any_flip:
            st.markdown(
                """
                <div class="theorem-card">
                <b>Result:</b> At least one anomaly-free gauge group produces
                a positive total Casimir coefficient. Matter loop corrections
                CAN flip the sign, opening the door to Casimir stabilization
                with the full matter content.
                </div>
                """,
                unsafe_allow_html=True,
            )
        else:
            st.markdown(
                """
                <div class="warning-card">
                <b>Result:</b> No anomaly-free gauge group flips the sign of
                the total Casimir coefficient. The fermion spectral zeta
                vanishes at s = &minus;1/2 (B<sub>3</sub>(1/2) = 0 exactly),
                so fermions contribute zero at leading order. The vector boson
                contribution alone is insufficient to overcome the negative
                graviton Casimir energy.<br><br>
                <b>The Casimir no-go persists even with full matter content.</b>
                Stabilization requires flux contributions (Page 17) or
                non-perturbative effects.
                </div>
                """,
                unsafe_allow_html=True,
            )
    else:
        st.markdown(
            """
            <div class="step-card">
            <b>Matter loop computation not available.</b><br>
            The matter Casimir correction requires the anomaly cancellation
            module. Run the full pipeline to see results.
            </div>
            """,
            unsafe_allow_html=True,
        )

    st.markdown("")

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Casimir Stabilization | Alpha Ladder Research Dashboard")
