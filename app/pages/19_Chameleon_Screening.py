"""
Chameleon Screening Consistency Check -- Page 19

Tests whether the Alpha Ladder dilaton can have an environment-dependent
mass via the chameleon/symmetron mechanism.  The key self-consistency check
is KK truncation: if the effective internal radius becomes macroscopic,
gravity turns 6-dimensional and the 4D effective theory breaks down.
"""

import streamlit as st
import math
import pandas as pd


# ---------------------------------------------------------------------------
# Custom CSS (matches Pages 16-18)
# ---------------------------------------------------------------------------
st.markdown(
    """
    <style>
    .proof-card {
        background-color: #1a1d23;
        border: 1px solid #2e3440;
        border-radius: 8px;
        padding: 1.2rem;
        margin: 0.5rem 0;
    }
    .formula-card {
        background-color: #1a1d23;
        border-left: 3px solid #f59e0b;
        padding: 0.8rem 1rem;
        margin: 0.4rem 0;
        border-radius: 0 8px 8px 0;
    }
    .theorem-card {
        background-color: #1a1d23;
        border-left: 3px solid #34d399;
        padding: 0.8rem 1rem;
        margin: 0.4rem 0;
        border-radius: 0 8px 8px 0;
    }
    .step-card {
        background-color: #1a1d23;
        border-left: 3px solid #60a5fa;
        padding: 0.8rem 1rem;
        margin: 0.4rem 0;
        border-radius: 0 8px 8px 0;
    }
    .warning-card {
        background-color: #1a1d23;
        border-left: 3px solid #f87171;
        padding: 0.8rem 1rem;
        margin: 0.4rem 0;
        border-radius: 0 8px 8px 0;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

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
_core_available = False
try:
    from alpha_ladder_core.chameleon_screening import (  # noqa: E402
        summarize_chameleon_screening,
        check_kk_truncation_validity,
        compute_meV_dark_sector_potential,
    )
    _core_available = True
except ImportError:
    pass

# ---------------------------------------------------------------------------
# Compute values (guarded, cached)
# ---------------------------------------------------------------------------
_summary = None

if _core_available:
    @st.cache_data
    def _cached_summarize():
        return summarize_chameleon_screening()

    try:
        _summary = _cached_summarize()
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Chameleon Screening Consistency Check")
st.markdown(
    "Can the dilaton mass depend on the local matter density, becoming "
    "heavy in the solar system (screened) and light in cosmic voids "
    "(unscreened)?  This page tests the self-consistency of the "
    "chameleon/symmetron mechanism within the Alpha Ladder framework."
)
st.divider()

# ---------------------------------------------------------------------------
# 4 Metric Cards
# ---------------------------------------------------------------------------
col_m1, col_m2, col_m3, col_m4 = st.columns(4)

with col_m1:
    st.metric(label="Mechanism", value="Chameleon")
with col_m2:
    st.metric(label="Critical test", value="KK truncation")
with col_m3:
    st.metric(label="Fuzzy DM", value="EXCLUDED")
with col_m4:
    st.metric(label="meV regime", value="PASSES")

st.markdown("")

# ---------------------------------------------------------------------------
# Section A: The Question
# ---------------------------------------------------------------------------
with st.expander("A. The Question", expanded=True):

    st.markdown(
        """
        <div class="step-card">
        <b>Can the dilaton mass depend on the environment?</b><br><br>
        The chameleon mechanism is a well-known screening strategy in
        scalar-tensor gravity: the effective mass of a scalar field
        increases in regions of high matter density, suppressing fifth-force
        effects in the solar system while allowing the field to be light
        (and cosmologically active) in low-density voids.<br><br>
        The symmetron mechanism achieves a similar effect through a
        symmetry-breaking potential whose vacuum expectation value depends
        on the local density.  In both cases, the field couples to matter
        density so that its effective potential -- and thus its mass --
        varies with environment.<br><br>
        For the Alpha Ladder dilaton, this would require the effective
        internal radius a<sub>0,eff</sub> to change with density.
        The self-consistency question is whether this change pushes
        a<sub>0,eff</sub> into a regime where the 4D effective theory
        (the KK truncation) breaks down.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Section B: KK Truncation Check
# ---------------------------------------------------------------------------
with st.expander("B. KK Truncation Check", expanded=True):

    st.markdown(
        """
        <div class="formula-card">
        <b>The critical self-consistency test</b><br><br>
        If the effective internal radius a<sub>0,eff</sub>(&rho;) exceeds
        the Eot-Wash bound (~56 &mu;m), gravity becomes higher-dimensional
        at distances below a<sub>0,eff</sub>.  This is catastrophically
        excluded by every gravity experiment ever performed.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    if _core_available and _summary:
        consistency = _summary.get("self_consistency", {})
        env_results = consistency.get("environment_results", [])

        if env_results:
            table_rows = []
            for env in env_results:
                rho_si = env.get("rho_si", 0)
                a_eff = env.get("a_eff_m")
                m_phi = env.get("m_phi_eff_eV")
                kk_valid = env.get("kk_valid", False)
                status = env.get("exclusion_status", "unknown")

                table_rows.append({
                    "Environment": env.get("label", env.get("name", "?")),
                    "rho (kg/m^3)": f"{rho_si:.2e}",
                    "a_eff (m)": fmt_decimal(a_eff, sig_figs=3) if a_eff is not None else "N/A",
                    "m_phi (eV)": fmt_decimal(m_phi, sig_figs=3) if m_phi is not None else "N/A",
                    "KK Valid": "Yes" if kk_valid else "No",
                    "Status": status,
                })

            st.dataframe(
                pd.DataFrame(table_rows),
                use_container_width=True,
                hide_index=True,
            )
        else:
            st.info("No environment results available from the core module.")

        st.markdown("")

        st.markdown(
            """
            <div class="warning-card">
            <b>Key finding:</b> For fuzzy DM (m<sub>&phi;</sub> ~
            10<sup>&minus;22</sup> eV), the required internal radius is
            a<sub>0</sub> ~ 10<sup>7</sup> m.  At this scale, the KK tower
            spacing collapses and gravity becomes 6-dimensional below
            ~10,000 km.  This is excluded by every sub-mm and solar-system
            gravity test.
            </div>
            """,
            unsafe_allow_html=True,
        )

        # Chameleon profile chart
        profile = _summary.get("chameleon_profile")
        if profile:
            try:
                from app.components.charts import chameleon_profile_chart  # noqa: E402
                fig_profile = chameleon_profile_chart(profile)
                st.plotly_chart(fig_profile, use_container_width=True)
            except ImportError:
                st.info("Chart function chameleon_profile_chart not yet available in charts module.")
            except Exception as exc:
                st.warning(f"Chameleon profile chart error: {exc}")

            st.markdown("")

            try:
                from app.components.charts import kk_truncation_chart  # noqa: E402
                fig_kk = kk_truncation_chart(profile)
                st.plotly_chart(fig_kk, use_container_width=True)
            except ImportError:
                st.info("Chart function kk_truncation_chart not yet available in charts module.")
            except Exception as exc:
                st.warning(f"KK truncation chart error: {exc}")
    else:
        st.markdown(
            """
            <div class="step-card">
            <b>Environment scan:</b><br><br>
            The chameleon profile is computed across four representative
            environments:<br><br>
            &bull; Cosmic void (&rho; ~ 10<sup>&minus;30</sup> kg/m&sup3;)<br>
            &bull; Interplanetary medium (&rho; ~ 10<sup>&minus;18</sup> kg/m&sup3;)<br>
            &bull; Laboratory (&rho; ~ 10 kg/m&sup3;)<br>
            &bull; Earth surface (&rho; ~ 5500 kg/m&sup3;)<br><br>
            At each density, the effective internal radius a<sub>0,eff</sub>
            and dilaton mass m<sub>&phi;,eff</sub> are computed from the
            density-dependent minimum of the flux-stabilized potential.
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section C: The No-Go Result
# ---------------------------------------------------------------------------
with st.expander("C. The No-Go Result", expanded=True):

    if _core_available and _summary:
        exclusion_reason = _summary.get("fuzzy_dm_exclusion_reason", "")
        if exclusion_reason:
            st.error(exclusion_reason)

        st.markdown("")

        st.markdown(
            """
            <div class="warning-card">
            <b>Why chameleon screening cannot rescue fuzzy DM:</b><br><br>
            The chameleon mechanism couples the dilaton potential to the
            local matter density.  To achieve m<sub>&phi;</sub> ~
            10<sup>&minus;22</sup> eV in voids, the field must sit at a
            minimum where a<sub>0,eff</sub> ~ 10<sup>7</sup> m.  Even if
            the field shifts to a different minimum in high-density regions,
            the low-density regime is where fuzzy DM phenomenology operates
            -- and that regime is precisely where the KK truncation
            breaks down.<br><br>
            This is not a failure of the screening mechanism itself.  It is a
            structural incompatibility between the extra-dimensional origin
            of the dilaton and the enormous Compton wavelength required for
            fuzzy DM behavior.
            </div>
            """,
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            """
            <div class="warning-card">
            <b>No-go result (summary):</b><br><br>
            The fuzzy DM regime requires a<sub>0</sub> ~ 10<sup>7</sup> m.
            At this scale, gravity is 6-dimensional below ~10,000 km.
            The entire 4D effective theory -- including the chameleon
            mechanism -- loses self-consistency.
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section D: The meV Regime
# ---------------------------------------------------------------------------
with st.expander("D. The meV Regime", expanded=True):

    if _core_available and _summary:
        mev = _summary.get("mev_analysis", {})

        m_phi_eV = mev.get("m_phi_eV")
        lambda_c = mev.get("lambda_compton_m")
        is_fifth_force = mev.get("is_testable_fifth_force", False)

        col_d1, col_d2, col_d3 = st.columns(3)

        with col_d1:
            st.metric(
                label="m_phi",
                value=fmt_decimal(m_phi_eV, sig_figs=3) + " eV" if m_phi_eV is not None else "N/A",
            )
        with col_d2:
            st.metric(
                label="Compton wavelength",
                value=fmt_decimal(lambda_c, sig_figs=3) + " m" if lambda_c is not None else "N/A",
            )
        with col_d3:
            st.metric(
                label="Fifth force testable",
                value="YES" if is_fifth_force else "NO",
            )

        st.markdown("")

        if is_fifth_force:
            st.success(
                "The meV dilaton passes all KK truncation checks and produces "
                "a testable Yukawa fifth force at sub-millimeter scales.  "
                "This is the self-consistent regime of the Alpha Ladder framework."
            )
        else:
            st.warning(
                "The meV regime does not produce a testable fifth force with "
                "the current parameters."
            )

        st.markdown("")

        mev_status = _summary.get("mev_regime_status", "")
        if mev_status:
            st.markdown(
                f"""
                <div class="theorem-card">
                <b>meV regime assessment:</b><br><br>
                {mev_status}
                </div>
                """,
                unsafe_allow_html=True,
            )

        honest = mev.get("honest_assessment", "")
        if honest:
            st.markdown("")
            st.markdown(
                f"""
                <div class="step-card">
                <b>Honest assessment:</b><br><br>
                {honest}
                </div>
                """,
                unsafe_allow_html=True,
            )
    else:
        st.markdown(
            """
            <div class="theorem-card">
            <b>meV regime (a<sub>0</sub> ~ 30 &mu;m):</b><br><br>
            &bull; m<sub>&phi;</sub> ~ 7 meV<br>
            &bull; Compton wavelength ~ 30 &mu;m<br>
            &bull; KK truncation: VALID (a<sub>0</sub> &lt; 56 &mu;m)<br>
            &bull; Produces testable Yukawa fifth force<br>
            &bull; NOT a dark matter candidate
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section E: Overall Verdict
# ---------------------------------------------------------------------------
with st.expander("E. Overall Verdict", expanded=True):

    if _core_available and _summary:
        key_finding = _summary.get("key_finding", "")
        if key_finding:
            st.info(key_finding)

        st.markdown("")

        overall = _summary.get("overall_verdict", "")
        if overall:
            st.markdown(
                f"""
                <div class="proof-card">
                <b>Overall verdict:</b><br><br>
                {overall}
                </div>
                """,
                unsafe_allow_html=True,
            )

        st.markdown("")

        st.markdown(
            """
            <div class="theorem-card">
            <b>Summary of results:</b><br><br>
            &bull; Fuzzy DM (m ~ 10<sup>&minus;22</sup> eV): <b>EXCLUDED</b>
            by KK truncation breakdown<br>
            &bull; meV regime (m ~ 1&ndash;10 meV): <b>SELF-CONSISTENT</b>,
            testable fifth force<br>
            &bull; Planck regime (m ~ M<sub>Pl</sub>): pure GR, no
            observable chameleon effects<br><br>
            The chameleon mechanism is structurally incompatible with using
            the Alpha Ladder dilaton as fuzzy dark matter.  The meV regime
            remains the only self-consistent phenomenological window.
            </div>
            """,
            unsafe_allow_html=True,
        )
    else:
        st.info(
            "The chameleon mechanism CANNOT make the Alpha Ladder dilaton a "
            "fuzzy dark matter candidate: the required internal radius "
            "(~10^7 m) makes gravity 6-dimensional at observable scales.  "
            "The meV regime is self-consistent but phenomenologically limited "
            "to sub-mm fifth-force searches."
        )

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Chameleon Screening | Alpha Ladder Research Dashboard")
