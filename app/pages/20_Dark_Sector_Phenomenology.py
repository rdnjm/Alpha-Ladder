"""
Dark Sector Phenomenology -- Page 20

Comprehensive analysis of the Alpha Ladder dilaton as a potential dark
sector participant.  The key result: the meV dilaton is a testable
sub-mm fifth force, NOT a dark matter or dark energy candidate.
"""

import streamlit as st
import math
import pandas as pd

st.set_page_config(page_title="Dark Sector | Alpha Ladder", layout="wide")

# ---------------------------------------------------------------------------
# Custom CSS (matches Pages 16-19)
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
    from alpha_ladder_core.dark_sector import (  # noqa: E402
        summarize_dark_sector,
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
        return summarize_dark_sector()

    try:
        _summary = _cached_summarize()
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Dark Sector Phenomenology")
st.markdown(
    "What role does the Alpha Ladder dilaton play in dark sector physics?  "
    "This page evaluates the dilaton against observational constraints on "
    "dark matter, dark energy, and self-interacting dark matter.  The "
    "honest result: the dilaton is a testable sub-mm fifth force, not a "
    "dark sector component."
)
st.divider()

# ---------------------------------------------------------------------------
# 4 Metric Cards
# ---------------------------------------------------------------------------
col_m1, col_m2, col_m3, col_m4 = st.columns(4)

with col_m1:
    if _summary:
        mev = _summary.get("mev_analysis", {})
        m_meV = mev.get("m_phi_meV")
        st.metric(label="m_phi", value=f"{m_meV:.1f} meV" if m_meV else "~7 meV")
    else:
        st.metric(label="m_phi", value="~7 meV")
with col_m2:
    st.metric(label="Dark matter", value="NO")
with col_m3:
    st.metric(label="Dark energy", value="NO")
with col_m4:
    st.metric(label="Fifth force", value="YES")

st.markdown("")

# ---------------------------------------------------------------------------
# Section A: Dark Sector Landscape
# ---------------------------------------------------------------------------
with st.expander("A. Dark Sector Landscape", expanded=True):

    st.markdown(
        """
        <div class="step-card">
        <b>Mass landscape overview:</b> The dilaton's phenomenology depends
        entirely on its mass.  We survey eleven representative mass scales
        spanning 61 orders of magnitude, from dark-energy scales
        (10<sup>&minus;33</sup> eV) to the Planck mass (10<sup>28</sup> eV),
        checking KK truncation validity and classifying each regime.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    if _core_available and _summary:
        landscape = _summary.get("landscape")
        if landscape and isinstance(landscape, dict):
            entries = landscape.get("landscape", [])
        elif isinstance(landscape, list):
            entries = landscape
        else:
            entries = []

        if entries:
            table_rows = []
            for entry in entries:
                m_eV = entry.get("m_phi_eV", 0)
                a_0 = entry.get("a_0_m")
                kk_valid = entry.get("kk_truncation_valid", False)
                classification = entry.get("classification", "unknown")

                table_rows.append({
                    "m_phi (eV)": f"{m_eV:.2e}",
                    "a_0 (m)": f"{a_0:.2e}" if a_0 is not None else "N/A",
                    "KK Valid": "Yes" if kk_valid else "No",
                    "Classification": classification.replace("_", " ").title(),
                })

            st.dataframe(
                pd.DataFrame(table_rows),
                use_container_width=True,
                hide_index=True,
            )

            st.markdown("")

            # Landscape chart
            try:
                from app.components.charts import dark_sector_landscape_chart  # noqa: E402
                fig_landscape = dark_sector_landscape_chart(
                    landscape if isinstance(landscape, dict) else {"landscape": entries}
                )
                st.plotly_chart(fig_landscape, use_container_width=True)
            except ImportError:
                st.info(
                    "Chart function dark_sector_landscape_chart not yet "
                    "available in charts module."
                )
            except Exception as exc:
                st.warning(f"Landscape chart error: {exc}")
        else:
            st.info("No landscape data available.")
    else:
        st.markdown(
            """
            <div class="formula-card">
            <b>Classification scheme:</b><br><br>
            &bull; <b>Dark energy candidate:</b> m &lt; 10<sup>&minus;30</sup> eV<br>
            &bull; <b>Fuzzy DM candidate:</b> m ~ 10<sup>&minus;22</sup> eV
            (excluded by KK)<br>
            &bull; <b>Sub-mm testable:</b> m ~ 10<sup>&minus;4</sup>&ndash;10<sup>&minus;1</sup> eV<br>
            &bull; <b>Thermal relic range:</b> m ~ 1&ndash;10<sup>6</sup> eV<br>
            &bull; <b>Invisible Planck:</b> m &gt; 10<sup>25</sup> eV
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section B: Relic Abundance
# ---------------------------------------------------------------------------
with st.expander("B. Relic Abundance", expanded=True):

    st.markdown(
        """
        <div class="formula-card">
        <b>Misalignment mechanism:</b><br><br>
        <code>&Omega;<sub>&phi;</sub> h&sup2; ~
        (m<sub>&phi;</sub> / 10<sup>&minus;22</sup> eV)<sup>1/2</sup>
        &middot; (f / M<sub>Pl</sub>)&sup2;</code><br><br>
        With f = M<sub>Pl</sub> and m<sub>&phi;</sub> ~ 7 meV, the
        dilaton overproduces dark matter by many orders of magnitude.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    if _core_available and _summary:
        relic = _summary.get("relic_abundance", {})

        Omega_h2 = relic.get("Omega_phi_h2")
        ratio = relic.get("ratio_to_observed")
        f_over_mpl = relic.get("f_required_over_M_Pl")

        col_b1, col_b2, col_b3 = st.columns(3)

        with col_b1:
            st.metric(
                label="Omega_phi h^2",
                value=fmt_decimal(Omega_h2, sig_figs=3) if Omega_h2 is not None else "N/A",
            )
        with col_b2:
            st.metric(
                label="Ratio to observed",
                value=fmt_decimal(ratio, sig_figs=3) if ratio is not None else "N/A",
            )
        with col_b3:
            st.metric(
                label="f_required / M_Pl",
                value=fmt_decimal(f_over_mpl, sig_figs=3) if f_over_mpl is not None else "N/A",
            )

        st.markdown("")

        overproduced = relic.get("overproduced", True)
        if overproduced:
            st.warning(
                "The dilaton with f = M_Pl OVERPRODUCES dark matter.  "
                "The relic abundance exceeds the observed value by many "
                "orders of magnitude.  The dilaton is not a viable dark "
                "matter candidate via the misalignment mechanism unless "
                "f is severely fine-tuned."
            )

        honest = relic.get("honest_assessment", "")
        if honest:
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
            <div class="step-card">
            <b>Expected result:</b><br><br>
            With m ~ 7 meV and f = M<sub>Pl</sub>, the misalignment
            mechanism produces &Omega;<sub>&phi;</sub> h&sup2; many
            orders of magnitude above the observed 0.120.  The dilaton
            cannot be all of dark matter.
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section C: Self-Interaction
# ---------------------------------------------------------------------------
with st.expander("C. Self-Interaction", expanded=True):

    st.markdown(
        """
        <div class="formula-card">
        <b>Self-interaction cross section:</b><br><br>
        <code>&sigma;/m = &alpha;&sup2; / (4&pi; m<sub>&phi;</sub>&sup3;)</code><br><br>
        Observational bounds: &sigma;/m &lt; 1 cm&sup2;/g (bullet cluster),
        &sigma;/m &lt; 0.1 cm&sup2;/g (Milky Way halos).
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    if _core_available and _summary:
        si = _summary.get("self_interaction", {})

        sigma_over_m = si.get("sigma_over_m_cm2_g")
        bc_bound = si.get("bullet_cluster_bound", 1.0)
        mw_bound = si.get("milky_way_bound", 0.1)
        passes_bc = si.get("passes_bullet_cluster", False)
        passes_mw = si.get("passes_milky_way", False)

        col_c1, col_c2, col_c3 = st.columns(3)

        with col_c1:
            st.metric(
                label="sigma/m (cm^2/g)",
                value=fmt_decimal(sigma_over_m, sig_figs=3) if sigma_over_m is not None else "N/A",
            )
        with col_c2:
            st.metric(
                label="Bullet cluster",
                value="PASS" if passes_bc else "FAIL",
            )
        with col_c3:
            st.metric(
                label="Milky Way",
                value="PASS" if passes_mw else "FAIL",
            )

        honest = si.get("honest_assessment", "")
        if honest:
            st.markdown("")
            st.markdown(
                f"""
                <div class="step-card">
                <b>Assessment:</b><br><br>
                {honest}
                </div>
                """,
                unsafe_allow_html=True,
            )
    else:
        st.markdown(
            """
            <div class="step-card">
            <b>Expected result:</b><br><br>
            The meV dilaton self-interaction cross section is negligible
            compared to bullet cluster and Milky Way bounds.  The dilaton
            does not behave as self-interacting dark matter.
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section D: Equation of State
# ---------------------------------------------------------------------------
with st.expander("D. Equation of State", expanded=True):

    st.markdown(
        """
        <div class="formula-card">
        <b>Equation of state:</b><br><br>
        <code>w(z) = &minus;1 / (1 + (m<sub>&phi;</sub> / H(z))&sup2;)</code><br><br>
        When H &gt;&gt; m<sub>&phi;</sub>: w &rarr; &minus;1 (frozen, dark-energy-like).<br>
        When H &lt;&lt; m<sub>&phi;</sub>: w &rarr; 0 (oscillating, matter-like).
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    if _core_available and _summary:
        eos = _summary.get("equation_of_state", {})

        w_today = eos.get("w_today")
        z_trans = eos.get("z_transition")
        behaves_matter = eos.get("behaves_as_matter", False)

        col_d1, col_d2, col_d3 = st.columns(3)

        with col_d1:
            st.metric(
                label="w(today)",
                value=f"{w_today:.4f}" if w_today is not None else "N/A",
            )
        with col_d2:
            st.metric(
                label="z_transition",
                value=fmt_decimal(z_trans, sig_figs=3) if z_trans is not None else "N/A",
            )
        with col_d3:
            st.metric(
                label="Behaves as matter",
                value="YES" if behaves_matter else "NO",
            )

        honest = eos.get("honest_assessment", "")
        if honest:
            st.markdown("")
            st.markdown(
                f"""
                <div class="step-card">
                <b>Assessment:</b><br><br>
                {honest}
                </div>
                """,
                unsafe_allow_html=True,
            )
    else:
        st.markdown(
            """
            <div class="step-card">
            <b>Expected result:</b><br><br>
            For m ~ 7 meV, H<sub>0</sub> / m &lt;&lt; 1, so w &approx; 0
            today.  The dilaton behaves as pressureless matter (not dark
            energy) throughout all observable history.  The transition from
            w = &minus;1 to w = 0 occurred at extremely high redshift.
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section E: Fuzzy DM Constraints
# ---------------------------------------------------------------------------
with st.expander("E. Fuzzy DM Constraints", expanded=True):

    if _core_available and _summary:
        fuzzy = _summary.get("fuzzy_constraints", {})

        constraints = fuzzy.get("constraints", [])
        if constraints:
            constraint_rows = []
            for c in constraints:
                constraint_rows.append({
                    "Constraint": c.get("name", "?"),
                    "Bound (eV)": f"{c.get('bound_eV', 0):.2e}",
                    "Passes": "Yes" if c.get("passes", False) else "No",
                })

            st.dataframe(
                pd.DataFrame(constraint_rows),
                use_container_width=True,
                hide_index=True,
            )

        st.markdown("")

        db_wave_m = fuzzy.get("de_broglie_wavelength_m")
        db_wave_kpc = fuzzy.get("de_broglie_wavelength_kpc")
        is_fuzzy = fuzzy.get("is_fuzzy_candidate", False)

        col_e1, col_e2, col_e3 = st.columns(3)

        with col_e1:
            st.metric(
                label="de Broglie wavelength",
                value=fmt_decimal(db_wave_m, sig_figs=3) + " m" if db_wave_m is not None else "N/A",
            )
        with col_e2:
            st.metric(
                label="de Broglie wavelength",
                value=fmt_decimal(db_wave_kpc, sig_figs=3) + " kpc" if db_wave_kpc is not None else "N/A",
            )
        with col_e3:
            st.metric(
                label="Fuzzy DM candidate",
                value="YES" if is_fuzzy else "NO",
            )

        honest = fuzzy.get("honest_assessment", "")
        if honest:
            st.markdown("")
            st.markdown(
                f"""
                <div class="warning-card">
                <b>Assessment:</b><br><br>
                {honest}
                </div>
                """,
                unsafe_allow_html=True,
            )
    else:
        st.markdown(
            """
            <div class="warning-card">
            <b>Fuzzy DM assessment:</b><br><br>
            The meV dilaton has a de Broglie wavelength many orders of
            magnitude below the 0.1 kpc threshold required for fuzzy DM
            behavior.  It passes all observational mass bounds (Lyman-alpha,
            UV luminosity, CMB lensing, subhalo count) trivially because
            it is far too heavy, but it cannot produce the wave-like
            interference patterns characteristic of fuzzy dark matter.
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section F: Framework Position
# ---------------------------------------------------------------------------
with st.expander("F. Framework Position", expanded=True):

    if _core_available and _summary:
        framework_position = _summary.get("framework_position", "")
        if framework_position:
            st.markdown(
                f"""
                <div class="proof-card">
                <b>Framework position:</b><br><br>
                {framework_position}
                </div>
                """,
                unsafe_allow_html=True,
            )

        st.markdown("")

        key_finding = _summary.get("key_finding", "")
        if key_finding:
            st.info(key_finding)

        st.markdown("")

        overall = _summary.get("overall_verdict", "")
        if overall:
            st.markdown(
                f"""
                <div class="theorem-card">
                <b>Overall verdict:</b><br><br>
                {overall}
                </div>
                """,
                unsafe_allow_html=True,
            )
    else:
        st.info(
            "The Alpha Ladder dilaton at a_0 = 30 um (m ~ 7 meV) is a testable "
            "sub-mm fifth force, NOT a dark matter candidate.  The dark sector "
            "remains unexplained by this framework -- and that honesty is part "
            "of the theory's value."
        )

        st.markdown("")

        st.markdown(
            """
            <div class="theorem-card">
            <b>What the framework DOES predict:</b><br><br>
            &bull; A specific value of G from &alpha;<br>
            &bull; A sub-mm Yukawa correction from the dilaton<br>
            &bull; Pure GR at long distances<br><br>
            <b>What the framework does NOT explain:</b><br><br>
            &bull; Dark matter<br>
            &bull; Dark energy<br>
            &bull; The cosmological constant
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Dark Sector Phenomenology | Alpha Ladder Research Dashboard")
