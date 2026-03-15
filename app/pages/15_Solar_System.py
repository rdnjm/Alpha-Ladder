"""
Solar System Compatibility -- Page 15

Demonstrates that the massive dilaton passes Cassini and fifth-force
exclusion bounds, structures the screening discrepancy honestly,
and positions the framework in the literature.
"""

import streamlit as st
import math


# ---------------------------------------------------------------------------
# Custom CSS
# ---------------------------------------------------------------------------
st.markdown(
    """
    <style>
    .solar-card {
        background-color: #1a1d23;
        border: 1px solid #2e3440;
        border-radius: 8px;
        padding: 1.2rem;
        margin: 0.5rem 0;
    }
    .pass-card {
        background-color: #1a1d23;
        border-left: 3px solid #34d399;
        padding: 0.8rem 1rem;
        margin: 0.4rem 0;
        border-radius: 0 8px 8px 0;
    }
    .warn-card {
        background-color: #1a1d23;
        border-left: 3px solid #f59e0b;
        padding: 0.8rem 1rem;
        margin: 0.4rem 0;
        border-radius: 0 8px 8px 0;
    }
    .fail-card {
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
    from alpha_ladder_core.solar_system import (  # noqa: E402
        compute_ppn_parameters,
        compute_ppn_profile,
        compute_fifth_force_point,
        check_dilaton_exclusion,
        compute_minimum_dilaton_mass,
        summarize_solar_system,
    )
    from alpha_ladder_core.theoretical_context import (  # noqa: E402
        compute_screening_discrepancy,
        position_in_literature,
        analyze_anomaly_status,
        summarize_theoretical_status,
    )
    from alpha_ladder_core.dilaton import compute_bd_parameter  # noqa: E402
    _core_available = True
except ImportError:
    pass

# ---------------------------------------------------------------------------
# Chart imports
# ---------------------------------------------------------------------------
try:
    from app.components.charts import (  # noqa: E402
        ppn_profile_chart,
        fifth_force_exclusion_chart,
        screening_discrepancy_chart,
        literature_comparison_chart,
    )
    _charts_available = True
except ImportError:
    _charts_available = False

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Solar System Compatibility")
st.markdown(
    "Demonstrating that the massive dilaton passes Cassini tracking, "
    "fifth-force exclusion bounds, and honestly assessing the screening "
    "amplitude discrepancy."
)

if not _core_available or constants is None:
    st.warning("Core modules not available. Please check your installation.")
    st.stop()

# ---------------------------------------------------------------------------
# Compute values
# ---------------------------------------------------------------------------
bd = compute_bd_parameter(constants)
omega = bd["omega"]
m_phi_eV = bd["dilaton_mass_min_eV"]

ppn = compute_ppn_parameters(omega, m_phi_eV)
ppn_profile = compute_ppn_profile(omega, m_phi_eV)
fifth_force_point = compute_fifth_force_point(omega, m_phi_eV)
exclusion = check_dilaton_exclusion(omega, m_phi_eV)
min_mass = compute_minimum_dilaton_mass(omega)
theoretical = summarize_theoretical_status(constants)

# ---------------------------------------------------------------------------
# 4 Metric Cards
# ---------------------------------------------------------------------------
c1, c2, c3, c4 = st.columns(4)
with c1:
    st.metric(
        label="gamma_PPN at 1 AU",
        value=f"{ppn['gamma_PPN_at_cassini']:.10f}",
    )
with c2:
    st.metric(
        label="Cassini bound |gamma-1|",
        value=f"{ppn['cassini_bound']:.1e}",
    )
with c3:
    st.metric(
        label="Min dilaton mass",
        value=fmt_decimal(min_mass["m_phi_min_eV"]),
    )
with c4:
    st.metric(
        label="Compton wavelength",
        value=fmt_decimal(min_mass["lambda_max_m"]),
    )

st.divider()

# ---------------------------------------------------------------------------
# Section A: PPN Parameters
# ---------------------------------------------------------------------------
with st.expander("**A. PPN Parameters**", expanded=True):
    if _charts_available:
        fig = ppn_profile_chart(ppn_profile)
        st.plotly_chart(fig, use_container_width=True)

    # Landmark table
    st.subheader("Landmark Values")
    landmarks = ppn_profile.get("landmarks", {})
    if landmarks:
        import pandas as pd
        rows = []
        for key, lm in landmarks.items():
            passes = lm["gamma_deviation"] < ppn["cassini_bound"]
            rows.append({
                "Location": lm.get("label", key),
                "Distance (m)": f"{lm['r_meters']:.3e}",
                "gamma_PPN": f"{lm['gamma_PPN']:.12f}",
                "|gamma - 1|": f"{lm['gamma_deviation']:.3e}",
                "Passes Cassini": "Yes" if passes else "No",
            })
        df = pd.DataFrame(rows)
        st.dataframe(df, use_container_width=True, hide_index=True)

    # Transition radius
    if ppn_profile.get("transition_radius_m") is not None:
        tr = ppn_profile["transition_radius_m"]
        AU = 1.496e11
        st.markdown(
            f'<div class="pass-card">'
            f"<b>Result:</b> The massive dilaton passes the Cassini bound at distances "
            f"beyond {tr:.3e} m ({tr/AU:.4f} AU). At 1 AU, |gamma - 1| = "
            f"{ppn['gamma_deviation_at_cassini']:.2e}, safely below the Cassini limit of "
            f"{ppn['cassini_bound']:.1e}. The dilaton mass provides "
            f"{min_mass['suppression_factor']:.1f} e-foldings of Yukawa suppression."
            f"</div>",
            unsafe_allow_html=True,
        )
    elif ppn["passes_cassini"]:
        st.markdown(
            '<div class="pass-card"><b>Result:</b> Massive dilaton passes Cassini bound.</div>',
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section B: Fifth-Force Exclusion
# ---------------------------------------------------------------------------
with st.expander("**B. Fifth-Force Exclusion**", expanded=True):
    if _charts_available:
        fig = fifth_force_exclusion_chart(
            exclusion["dilaton_point"],
            exclusion["bounds"],
        )
        st.plotly_chart(fig, use_container_width=True)

    if exclusion["allowed"]:
        closest = exclusion.get("closest_bound")
        margin_text = ""
        if closest:
            margin_text = f" Closest bound: {closest['name']} (margin: {closest['margin_dex']:.1f} dex)."
        st.markdown(
            f'<div class="pass-card">'
            f"<b>Result:</b> The dilaton point lies in the allowed region, "
            f"below all experimental exclusion curves.{margin_text}"
            f"</div>",
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            f'<div class="fail-card">'
            f"<b>Result:</b> The dilaton is excluded by: "
            f"{', '.join(exclusion['excluded_by'])}."
            f"</div>",
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section C: Screening Amplitude
# ---------------------------------------------------------------------------
with st.expander("**C. Screening Amplitude Discrepancy**", expanded=True):
    discrepancy = theoretical["discrepancy"]

    if _charts_available:
        fig = screening_discrepancy_chart(discrepancy)
        st.plotly_chart(fig, use_container_width=True)

    st.markdown(
        f"**Tree-level prediction:** alpha = 2 beta^2 = {discrepancy['alpha_tree']:.4e}"
    )
    st.markdown(
        f"**Empirical value:** alpha = (G_lab - G_vac) / G_vac = {discrepancy['alpha_empirical']:.4e}"
    )
    st.markdown(
        f"**Ratio:** {discrepancy['ratio']:.0f}x (log10 = {discrepancy['log10_ratio']:.2f})"
    )

    st.subheader("Possible Resolutions")
    for res in discrepancy["resolutions"]:
        plaus = res["plausibility"]
        if plaus == "high":
            card_class = "pass-card"
        elif plaus == "moderate":
            card_class = "warn-card"
        else:
            card_class = "solar-card"
        st.markdown(
            f'<div class="{card_class}">'
            f"<b>{res['name']}</b> (plausibility: {plaus})<br>"
            f"<i>{res['mechanism']}</i><br>"
            f"{res['description']}"
            f"</div>",
            unsafe_allow_html=True,
        )

    # Amber warning
    st.markdown(
        '<div class="warn-card">'
        "<b>Status:</b> The ~3800x ratio is a <i>conditional</i> discrepancy. "
        "IF the dilaton is light (Compton wavelength &gt; 0.1 mm), the tree-level "
        "coupling is 3854x too strong. IF the dilaton is heavy (wavelength &lt; 0.1 mm), "
        "no discrepancy exists -- the Yukawa force is undetectable at all lab scales, "
        "and the sub-ppm G residual (corrected bridge) is well within experimental scatter. "
        "The dilaton mass origin remains an open problem."
        "</div>",
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Section D: Theoretical Context
# ---------------------------------------------------------------------------
with st.expander("**D. Theoretical Context**", expanded=True):
    literature = theoretical["literature"]
    anomalies = theoretical["anomalies"]

    if _charts_available:
        fig = literature_comparison_chart(literature)
        st.plotly_chart(fig, use_container_width=True)

    # Anomaly status
    st.subheader("Anomaly Status")
    for anom in anomalies["anomaly_types"]:
        status = anom["status"]
        if status == "absent":
            card_class = "pass-card"
        elif status == "not applicable":
            card_class = "solar-card"
        else:
            card_class = "warn-card"
        st.markdown(
            f'<div class="{card_class}">'
            f"<b>{anom['type']}</b> -- {status}<br>"
            f"{anom['description']}"
            f"</div>",
            unsafe_allow_html=True,
        )

    st.markdown(f"*{anomalies['interpretation']}*")

    # Strengths and open problems
    col_s, col_o = st.columns(2)
    with col_s:
        st.subheader("Strengths")
        for s in theoretical["strengths"]:
            st.markdown(f"- {s}")

    with col_o:
        st.subheader("Open Problems")
        for p in theoretical["open_problems"]:
            st.markdown(f"- {p}")

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Solar System Compatibility | Alpha Ladder Research Dashboard")
