"""
Page 32 -- Residual Mapping
Maps delta = sqrt(phi) - k_exact against SM constants.
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
    from alpha_ladder_core.residual_mapping import summarize_residual_mapping
    _core_available = True
except ImportError:
    pass

_charts_available = False
try:
    from app.components.charts import residual_scan_chart, k_closed_form_chart
    _charts_available = True
except ImportError:
    pass


@st.cache_data
def _get_summary(_constants):
    if _core_available:
        return summarize_residual_mapping(_constants)
    return None


st.title("Residual Mapping")
st.markdown("*Mapping δ = √φ − k_exact against Standard Model constants*")

summary = _get_summary(constants)

# --- A. The Residual ---
st.header("A. The Residual")

if summary:
    di = summary["delta_info"]
    c1, c2, c3 = st.columns(3)
    c1.metric("k_exact", f"{di['k_exact']:.6f}")
    c2.metric("√φ", f"{di['sqrt_phi']:.6f}")
    c3.metric("δ", f"{di['delta']:.6e}")

    st.markdown(f"""
<div class="formula-card">
<strong>The gap:</strong> δ = √φ − k_exact = {di['delta']:.6e}<br>
k_exact is the offset that makes α_g = α²⁴ × μ × (μ − k) match G exactly.<br>
√φ = {di['sqrt_phi']:.6f} overshoots by δ, giving the -0.31 ppm residual (with (1−α) correction).
</div>
""", unsafe_allow_html=True)
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# --- B. Delta Structure ---
st.header("B. Delta Structure Analysis")

if summary:
    ds = summary["delta_structure"]
    c1, c2, c3 = st.columns(3)
    c1.metric("δ / α", f"{ds['delta_over_alpha']:.4f}")
    c2.metric("α power", f"{ds['delta_as_alpha_power']:.2f}")
    c3.metric("δ × μ", f"{ds['delta_times_mu']:.2f}")

    st.markdown(f"""
<div class="step-card">
{ds['interpretation']}
</div>
""", unsafe_allow_html=True)

    rows = [
        {"Quantity": "δ / α", "Value": f"{ds['delta_over_alpha']:.6f}"},
        {"Quantity": "δ / α²", "Value": f"{ds['delta_over_alpha_sq']:.2f}"},
        {"Quantity": "δ × μ", "Value": f"{ds['delta_times_mu']:.4f}"},
        {"Quantity": "δ / Schwinger", "Value": f"{ds['delta_over_schwinger']:.4f}"},
        {"Quantity": "δ / √φ", "Value": f"{ds['delta_over_sqrt_phi']:.6f}"},
        {"Quantity": "δ as αⁿ", "Value": f"n = {ds['delta_as_alpha_power']:.4f}"},
    ]
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

# --- C. Anomalous Magnetic Moment Scan ---
st.header("C. Anomalous Magnetic Moment Terms")

if summary:
    g2 = summary["g2_scan"]
    st.markdown(f"""
<div class="step-card">
{g2['assessment']}
</div>
""", unsafe_allow_html=True)

    rows = []
    for c in g2["candidates"]:
        rows.append({
            "Expression": c["name"],
            "Value": f"{c['value']:.6e}",
            "δ / expr": f"{c['ratio']:.4f}",
            "Nearest int": c.get("nearest_integer", ""),
            "Frac. part": f"{c.get('closeness_to_integer', ''):.4f}" if isinstance(c.get("closeness_to_integer"), (int, float)) and c["closeness_to_integer"] != float("inf") else "N/A",
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

# --- D. Zeta Function Scan ---
st.header("D. Zeta Function and Mathematical Constants")

if summary:
    zs = summary["zeta_scan"]
    st.markdown(f"""
<div class="step-card">
{zs['assessment']}
</div>
""", unsafe_allow_html=True)

    rows = []
    for c in zs["candidates"]:
        rows.append({
            "Expression": c["name"],
            "Value": f"{c['value']:.6e}",
            "δ / expr": f"{c['ratio']:.4f}",
            "Nearest int": c.get("nearest_integer", ""),
            "Frac. part": f"{c.get('closeness_to_integer', ''):.4f}" if isinstance(c.get("closeness_to_integer"), (int, float)) and c["closeness_to_integer"] != float("inf") else "N/A",
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

# --- E. Composite Expression Search ---
st.header("E. Composite Expression Search")

if summary:
    cs = summary["composite_scan"]

    c1, c2 = st.columns(2)
    c1.metric("Within 100 ppm", cs["n_within_100ppm"])
    c2.metric("Within 1000 ppm", cs["n_within_1000ppm"])

    st.markdown(f"""
<div class="theorem-card">
{cs['assessment']}
</div>
""", unsafe_allow_html=True)

    if _charts_available and cs["best_matches"]:
        st.plotly_chart(residual_scan_chart(cs), use_container_width=True)

    if cs["best_matches"]:
        rows = []
        for c in cs["best_matches"]:
            rows.append({
                "Expression": c["name"],
                "Value": f"{c['value']:.6e}",
                "Residual (ppm)": f"{c['residual_ppm']:.1f}",
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

# --- F. Closed-Form Search for k ---
st.header("F. Closed-Form Candidates for k_exact")

if summary:
    kf = summary["k_closed_forms"]

    st.markdown(f"""
<div class="theorem-card">
{kf['assessment']}
</div>
""", unsafe_allow_html=True)

    if _charts_available and kf["best_matches"]:
        st.plotly_chart(k_closed_form_chart(kf), use_container_width=True)

    if kf["best_matches"]:
        rows = []
        for c in kf["best_matches"]:
            rows.append({
                "Expression": c["name"],
                "Value": f"{c['value']:.6f}",
                "k_exact": f"{c['k_exact']:.6f}",
                "Residual (ppm)": f"{c['residual_ppm']:.1f}",
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

# --- G. Honest Assessment ---
st.header("G. Honest Assessment")

if summary:
    st.markdown(f"""
<div class="warning-card">
<strong>Key Finding:</strong> {summary['key_finding']}<br><br>
{summary['honest_assessment']}
</div>
""", unsafe_allow_html=True)

st.divider()
st.caption("Residual Mapping | Alpha Ladder Research Dashboard")
