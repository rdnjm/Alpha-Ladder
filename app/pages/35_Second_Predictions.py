"""
Page 35 -- Second Predictions
Testable predictions beyond G: time variation, M_6 scale, correlated variations.
"""

import sys
import os
import streamlit as st
import pandas as pd

st.set_page_config(page_title="Second Predictions | Alpha Ladder", layout="wide")

st.markdown("""
<style>
.proof-card {
    background: #1a1a2e;
    border-radius: 10px;
    padding: 20px;
    margin: 10px 0;
}
.formula-card {
    background: #1a1a2e;
    border-left: 4px solid #f59e0b;
    border-radius: 10px;
    padding: 20px;
    margin: 10px 0;
}
.theorem-card {
    background: #1a1a2e;
    border-left: 4px solid #34d399;
    border-radius: 10px;
    padding: 20px;
    margin: 10px 0;
}
.step-card {
    background: #1a1a2e;
    border-left: 4px solid #60a5fa;
    border-radius: 10px;
    padding: 20px;
    margin: 10px 0;
}
.warning-card {
    background: #1a1a2e;
    border-left: 4px solid #f87171;
    border-radius: 10px;
    padding: 20px;
    margin: 10px 0;
}
</style>
""", unsafe_allow_html=True)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from app.components.sidebar import render_sidebar

constants = render_sidebar()

_core_available = False
try:
    from alpha_ladder_core.second_predictions import summarize_second_predictions
    _core_available = True
except ImportError:
    pass


@st.cache_data
def _get_summary(_constants):
    if _core_available:
        return summarize_second_predictions(_constants)
    return None


st.title("Second Predictions")
st.markdown("*Testable predictions beyond G: time variation, M_6 scale, mu tension*")

summary = _get_summary(constants)

# --- A. Time Variation of G ---
st.header("A. Time Variation: dG/G = 24 * (dalpha/alpha)")

if summary:
    tv = summary["time_variation"]
    c1, c2 = st.columns(2)
    c1.metric("A coefficient (alpha)", tv["A_coefficient"])
    c2.metric("B coefficient (mu)", tv["B_coefficient"])

    st.markdown(f"""
<div class="theorem-card">
<strong>Master formula:</strong><br><br>
dG/G = {tv['A_coefficient']} * (dalpha/alpha) + {tv['B_coefficient']} * (dmu/mu)
</div>
""", unsafe_allow_html=True)

    comp = tv.get("comparison_with_other_frameworks", {})
    if comp:
        rows = []
        for framework, coefficient in comp.items():
            rows.append({
                "Framework": framework,
                "Predicted coefficient": coefficient,
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.markdown(f"""
<div class="step-card">
{tv['assessment']}
</div>
""", unsafe_allow_html=True)
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# --- B. Current Experimental Bounds ---
st.header("B. Current Experimental Bounds")

if summary:
    cb = summary["current_bounds"]
    c1, c2, c3 = st.columns(3)
    c1.metric("alpha bound (dalpha/alpha)", f"{cb['alpha_bound']:.2e}")
    c2.metric("mu bound (dmu/mu)", f"{cb['mu_bound']:.2e}")
    c3.metric("G bound (dG/G)", f"{cb['G_bound']:.2e}")

    st.markdown(f"""
<div class="step-card">
<strong>Predicted dG from alpha bound:</strong> {cb['predicted_dG_from_alpha_bound']:.2e}<br>
<strong>Ratio to current G sensitivity:</strong> {cb['ratio_to_current_sensitivity']:.1f}x
</div>
""", unsafe_allow_html=True)

    ratio = cb["ratio_to_current_sensitivity"]
    st.markdown(f"""
<div class="formula-card">
Currently {ratio:.0f}x below detection threshold
</div>
""", unsafe_allow_html=True)

    if not cb.get("testable", True):
        st.markdown("""
<div class="warning-card">
<strong>Not yet testable with current experiments.</strong>
The predicted time variation of G from the Alpha Ladder framework falls below
the sensitivity of current atomic clock and lunar laser ranging experiments.
</div>
""", unsafe_allow_html=True)

# --- C. M_6 Landscape ---
st.header("C. 6D Planck Mass from Compactification Radius")

if summary:
    m6 = summary["m6_landscape"]
    sw = m6["survival_window"]
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("a0 min (um)", f"{sw['a0_min']*1e6:.1f}")
    c2.metric("a0 max (um)", f"{sw['a0_max']*1e6:.1f}")
    c3.metric("M6 min (TeV)", f"{sw['M6_min_TeV']:.1f}")
    c4.metric("M6 max (TeV)", f"{sw['M6_max_TeV']:.1f}")

    st.markdown(f"""
<div class="theorem-card">
<strong>Key prediction:</strong><br><br>
If a_0 is in the survival window [{sw['a0_min']*1e6:.0f}-{sw['a0_max']*1e6:.0f} um],
then M_6 = {sw['M6_min_TeV']:.0f}-{sw['M6_max_TeV']:.0f} TeV (LHC scale).
</div>
""", unsafe_allow_html=True)

    scan = m6.get("scan_results", [])
    if scan:
        rows = []
        for r in scan:
            if r.get("unphysical"):
                status = "Unphysical"
            elif r.get("excluded_eot_wash"):
                status = "Excluded"
            elif r.get("in_survival_window"):
                status = "Survival window"
            else:
                status = "Below reach"
            rows.append({
                "a0 (m)": f"{r['a0_meters']:.2e}",
                "M_6 (TeV)": f"{r['M_6_TeV']:.2f}",
                "In survival window": "YES" if r.get("in_survival_window") else "",
                "Status": status,
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.markdown(f"""
<div class="step-card">
{m6['assessment']}
</div>
""", unsafe_allow_html=True)

# --- D. Correlated Variations ---
st.header("D. Correlated Constant Variations")

if summary:
    cv = summary["correlated_variations"]
    c1, c2 = st.columns(2)
    c1.metric("Predicted dG/G (EM only)", f"{cv['predicted_dG_G']:.2e}")
    c2.metric("Predicted dG/G (with QCD)", f"{cv['predicted_dG_G_with_qcd']:.2e}")

    st.markdown(f"""
<div class="step-card">
<strong>QCD coupling R ~ 38:</strong> Nuclear masses depend on alpha_s, which runs with alpha.
The QCD contribution amplifies the correlated variation by a factor of ~38.
</div>
""", unsafe_allow_html=True)

    st.markdown("""
<div class="formula-card">
<strong>With QCD coupling:</strong> dG/G = (24 + 2*38) * dalpha/alpha = 100 * dalpha/alpha
</div>
""", unsafe_allow_html=True)

    st.markdown(f"""
<div class="step-card">
{cv['assessment']}
</div>
""", unsafe_allow_html=True)

# --- E. Mu Tension ---
st.header("E. Bridge vs Mu-Structure Tension")

if summary:
    mt = summary["mu_tension"]
    c1, c2, c3 = st.columns(3)
    c1.metric("mu measured", f"{mt['mu_measured']:.6f}")
    c2.metric("mu from bridge", f"{mt['mu_predicted_from_bridge']:.6f}")
    c3.metric("Tension (sigma)", f"{mt['tension_sigma']:.1f}")

    st.markdown(f"""
<div class="warning-card">
<strong>Interpretation:</strong> {mt['interpretation']}
</div>
""", unsafe_allow_html=True)

    resolutions = mt.get("possible_resolutions", [])
    if resolutions:
        rows = []
        for i, res in enumerate(resolutions, 1):
            rows.append({
                "#": i,
                "Possible resolution": res,
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.markdown(f"""
<div class="step-card">
{mt['assessment']}
</div>
""", unsafe_allow_html=True)

# --- F. Strongest Prediction ---
st.header("F. Strongest Testable Prediction")

if summary:
    sp = summary["strongest_prediction"]
    st.markdown(f"""
<div class="theorem-card">
{sp}
</div>
""", unsafe_allow_html=True)

# --- G. Honest Assessment ---
st.header("G. Honest Assessment")

if summary:
    st.markdown(f"""
<div class="warning-card">
{summary['honest_assessment']}
</div>
""", unsafe_allow_html=True)

st.divider()
st.caption("Second Predictions | Alpha Ladder Research Dashboard")
