"""
Page 40 -- Charged Matter Loops
Monopole background modifies KK spectrum, spectral zeta sign flip, anomaly-free group scan.
"""

import sys
import os
import streamlit as st
import pandas as pd


sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from app.components.sidebar import render_sidebar

constants = render_sidebar()

_core_available = False
try:
    from alpha_ladder_core.charged_matter_loops import summarize_charged_matter_loops
    _core_available = True
except ImportError:
    pass

_charts_available = False
try:
    from app.components.charts import (
        zeta_comparison_chart,
        group_scan_chart,
    )
    _charts_available = True
except ImportError:
    pass


@st.cache_data
def _get_cml_summary(_constants):
    if _core_available:
        return summarize_charged_matter_loops(_constants)
    return None


# ---------------------------------------------------------------------------
# Title and Introduction
# ---------------------------------------------------------------------------
st.title("Charged Matter Loops")
st.markdown("*Monopole background and spectral zeta sign flip on S^2*")

st.markdown("""
<div class="theorem-card">
<strong>Key idea:</strong> A monopole background on S^2 modifies the KK spectrum
from j(j+1) to (j+n)(j+n+1) where n is the monopole number. This shifts the
spectral zeta function from the pure-gravity value zeta(-1/2) = -17/480
to the charged value zeta(-1/2) = +1/10, flipping the sign and enabling
Casimir-driven radius stabilization.
</div>
""", unsafe_allow_html=True)

if not _core_available:
    st.error(
        "Core module `alpha_ladder_core.charged_matter_loops` not found. "
        "Install alpha_ladder_core to enable computations on this page."
    )

summary = _get_cml_summary(constants)

# ---------------------------------------------------------------------------
# Section 1: Monopole Spectrum
# ---------------------------------------------------------------------------
st.header("1. Monopole Spectrum")

if summary:
    mono = summary.get("monopole_spectrum", {})
    if mono:
        c1, c2, c3 = st.columns(3)
        charge_q = mono.get("charge_q")
        c1.metric("Charge q", f"{charge_q}" if charge_q is not None else "N/A")
        j_min = mono.get("j_min")
        c2.metric("j_min", f"{j_min}" if j_min is not None else "N/A")
        n_levels = mono.get("n_levels")
        c3.metric("Number of levels", f"{n_levels}" if n_levels is not None else "N/A")

        levels = mono.get("levels", [])
        if levels:
            display_levels = levels[:10]
            rows = []
            for lvl in display_levels:
                rows.append({
                    "j": f"{lvl.get('j', 'N/A')}",
                    "Eigenvalue": f"{lvl.get('eigenvalue', 0):.4f}",
                    "Degeneracy": f"{lvl.get('degeneracy', 0)}",
                })
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
        else:
            st.info("No level data available.")
    else:
        c1, c2, c3 = st.columns(3)
        c1.metric("Charge q", "N/A")
        c2.metric("j_min", "N/A")
        c3.metric("Number of levels", "N/A")
        st.markdown("""
<div class="step-card">
<strong>Monopole spectrum:</strong> Data not available.
The monopole background shifts the angular momentum quantum number j by the monopole number n.
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from charged matter loops module.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 2: Spectral Zeta Sign Flip
# ---------------------------------------------------------------------------
st.header("2. Spectral Zeta Sign Flip")

if summary:
    zeta = summary.get("spectral_zeta", {})
    if zeta:
        if _charts_available:
            st.plotly_chart(zeta_comparison_chart(zeta), use_container_width=True)

        pure_grav = summary.get("pure_gravity_zeta")
        mono_zeta = summary.get("monopole_zeta")
        zeta_shift = summary.get("zeta_shift")

        st.markdown(f"""
<div class="proof-card">
<strong>Spectral zeta sign flip:</strong><br><br>
Neutral (pure gravity): zeta_S2(-1/2) = {pure_grav if pure_grav is not None else 'N/A'}
= -17/480<br>
Charged (monopole): zeta_S2(-1/2) = {mono_zeta if mono_zeta is not None else 'N/A'}
= +1/10<br><br>
Shift: {zeta_shift if zeta_shift is not None else 'N/A'}<br>
The sign flip from negative to positive is the key result: it means the
Casimir energy now has the correct sign to compete with flux energy and
produce a stable minimum.
</div>
""", unsafe_allow_html=True)

        comp = zeta.get("comparison_neutral")
        if comp is not None:
            st.markdown(f"""
<div class="step-card">
<strong>Comparison:</strong> Charged zeta value = {zeta.get('value', 'N/A')},
method = {zeta.get('method', 'N/A')}, converged = {zeta.get('converged', 'N/A')},
neutral comparison = {comp}
</div>
""", unsafe_allow_html=True)
    else:
        st.markdown("""
<div class="step-card">
<strong>Spectral zeta:</strong> Data not available.
The neutral zeta is -17/480 (negative, no stabilization) while the
charged zeta is +1/10 (positive, enabling stabilization).
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from charged matter loops module.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 3: E8xE8 G Correction
# ---------------------------------------------------------------------------
st.header("3. E8xE8 G Correction")

if summary:
    e8 = summary.get("e8_g_correction", {})
    if e8:
        c1, c2, c3 = st.columns(3)
        dg_grav = e8.get("delta_G_over_G_gravity")
        c1.metric("dG/G (gravity)", f"{dg_grav:.4e}" if dg_grav is not None else "N/A")
        dg_matt = e8.get("delta_G_over_G_matter")
        c2.metric("dG/G (matter)", f"{dg_matt:.4e}" if dg_matt is not None else "N/A")
        ratio = e8.get("ratio_to_target")
        c3.metric("Ratio to 3*alpha^2", f"{ratio:.4f}" if ratio is not None else "N/A")

        c1, c2, c3 = st.columns(3)
        dg_total = e8.get("delta_G_over_G_total")
        c1.metric("dG/G (total)", f"{dg_total:.4e}" if dg_total is not None else "N/A")
        target = e8.get("target_3_alpha_sq")
        c2.metric("Target (3*alpha^2)", f"{target:.4e}" if target is not None else "N/A")
        delta_zeta = e8.get("delta_zeta")
        c3.metric("Delta zeta", f"{delta_zeta:.6f}" if delta_zeta is not None else "N/A")

        assessment = e8.get("assessment", "")
        if assessment:
            st.markdown(f"""
<div class="step-card">
<strong>E8xE8 Assessment:</strong> {assessment}
</div>
""", unsafe_allow_html=True)
    else:
        c1, c2, c3 = st.columns(3)
        c1.metric("dG/G (gravity)", "N/A")
        c2.metric("dG/G (matter)", "N/A")
        c3.metric("Ratio to 3*alpha^2", "N/A")
        st.markdown("""
<div class="step-card">
<strong>E8xE8 G correction:</strong> Data not available.
The one-loop correction to G from charged matter fields in the E8xE8 spectrum.
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from charged matter loops module.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 4: Anomaly-Free Group Scan
# ---------------------------------------------------------------------------
st.header("4. Anomaly-Free Group Scan")

if summary:
    gscan = summary.get("group_scan", {})
    if gscan:
        if _charts_available:
            st.plotly_chart(group_scan_chart(gscan), use_container_width=True)

        results = gscan.get("results", [])
        if results:
            rows = []
            for entry in results:
                g_corr = entry.get("g_correction", {})
                cw_pot = entry.get("cw_potential", {})
                rows.append({
                    "Group": entry.get("group", "N/A"),
                    "n_hyper": f"{entry.get('n_hyper', 'N/A')}",
                    "n_vector": f"{entry.get('n_vector', 'N/A')}",
                    "Ratio to 3*alpha^2": f"{g_corr.get('ratio_to_target', 0):.4f}",
                    "eff. Lambda_6": f"{cw_pot.get('effective_lambda6', 0):.4e}",
                })
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

        st.markdown(f"""
<div class="step-card">
<strong>Group scan:</strong> {gscan.get('n_groups', 'N/A')} anomaly-free groups tested.
Any match 3*alpha^2: {gscan.get('any_match_3alpha2', 'N/A')}.
Any significant Lambda_6: {gscan.get('any_significant_lambda6', 'N/A')}.
</div>
""", unsafe_allow_html=True)

        assessment = gscan.get("assessment", "")
        if assessment:
            st.markdown(f"""
<div class="theorem-card">
<strong>Group Scan Assessment:</strong> {assessment}
</div>
""", unsafe_allow_html=True)
    else:
        st.info("Group scan data not available.")
elif _core_available:
    st.info("No data returned from charged matter loops module.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 5: Gap Impact
# ---------------------------------------------------------------------------
st.header("5. Gap Impact")

if summary:
    gap2 = summary.get("gap2_impact", "")
    honest = summary.get("honest_assessment", "")

    if gap2:
        st.markdown(f"""
<div class="step-card">
<strong>Gap 2 Impact:</strong><br><br>
{gap2}
</div>
""", unsafe_allow_html=True)
    else:
        st.markdown("""
<div class="step-card">
<strong>Gap 2 impact:</strong> Data not available.
</div>
""", unsafe_allow_html=True)

    if honest:
        st.markdown(f"""
<div class="warning-card">
<strong>Honest Assessment:</strong><br><br>
{honest}
</div>
""", unsafe_allow_html=True)
    else:
        st.markdown("""
<div class="step-card">
<strong>Assessment:</strong> No honest assessment available from the module.
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from charged matter loops module.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Charged Matter Loops | Alpha Ladder Research Dashboard")
