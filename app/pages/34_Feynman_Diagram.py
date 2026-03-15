"""
Page 34 -- Feynman Diagram Calculation
Explicit one-loop diagram: sigma -> KK photon loop -> sigma.
"""

import sys
import os
import streamlit as st
import pandas as pd


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
    from alpha_ladder_core.feynman_diagram import summarize_feynman_diagram
    _core_available = True
except ImportError:
    pass


@st.cache_data
def _get_summary(_constants):
    if _core_available:
        return summarize_feynman_diagram(_constants)
    return None


st.title("Feynman Diagram Calculation")
st.markdown("*Explicit one-loop self-energy: sigma -> KK photon loop -> sigma on S^2*")

summary = _get_summary(constants)

# --- A. The Diagram ---
st.header("A. The Diagram: sigma -> KK Photon -> sigma")

if summary:
    diag = summary.get("diagram_description", {})

    st.markdown("""
<div class="theorem-card">
<pre style="color: #e2e8f0; font-family: monospace; font-size: 14px; margin: 10px 0;">
sigma ----[vertex]---- KK photon (l, m) ----[vertex]---- sigma
                            |                    |
                            ------[loop]----------
</pre>
<br>
The volume modulus sigma couples to the U(1) gauge field through the gauge kinetic
function f(sigma) = Vol(S<sup>2</sup>)(sigma). One-loop self-energy from KK photon
modes on S<sup>2</sup>.
</div>
""", unsafe_allow_html=True)

    gc = summary.get("gauge_coupling")
    if gc:
        c1, c2 = st.columns(2)
        c1.metric("Vertex factor", f"{gc.get('vertex_factor', 'N/A')}")
        c2.metric("Gauge coupling index n", f"{gc.get('n', 'N/A')}")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# --- B. Gauge Kinetic Coupling ---
st.header("B. Gauge Kinetic Coupling")

if summary:
    gc = summary.get("gauge_coupling")
    if gc:
        c1, c2 = st.columns(2)
        c1.metric("Vertex factor", f"{gc.get('vertex_factor', 'N/A')}")
        c2.metric("Coupling index n", f"{gc.get('n', 'N/A')}")

        coupling_desc = gc.get("coupling_description", "")
        if coupling_desc:
            st.markdown(f"""
<div class="step-card">
{coupling_desc}
</div>
""", unsafe_allow_html=True)

        scan = gc.get("scan_results")
        if scan:
            rows = []
            for entry in scan:
                rows.append({
                    "n": entry.get("n", ""),
                    "Vertex factor": entry.get("vertex_factor", ""),
                    "Coupling form": entry.get("coupling_form", ""),
                })
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

# --- C. KK Photon Spectrum ---
st.header("C. KK Photon Spectrum")

if summary:
    kk = summary.get("kk_spectrum")
    if kk:
        modes = kk.get("modes", [])
        if modes:
            rows = []
            for mode in modes:
                rows.append({
                    "l": mode.get("l", ""),
                    "mass^2 / R^2": mode.get("mass_sq_over_R2", ""),
                    "Degeneracy (2l+1)": mode.get("degeneracy", ""),
                })
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

        explanation = kk.get("explanation", "")
        if explanation:
            st.markdown(f"""
<div class="step-card">
{explanation}
</div>
""", unsafe_allow_html=True)

# --- D. One-Loop Self-Energy ---
st.header("D. One-Loop Self-Energy")

if summary:
    se = summary.get("self_energy")
    if se:
        c1, c2, c3 = st.columns(3)
        c1.metric("Coupling factor", f"{se.get('coupling_factor', 'N/A')}")
        c2.metric("Loop factor", f"{se.get('loop_factor', 'N/A')}")
        c3.metric("Volume factor", f"{se.get('volume_factor', 'N/A')}")

        product = se.get("product", "")
        st.markdown(f"""
<div class="formula-card">
<strong>Cancellation:</strong><br><br>
coupling_factor * loop_factor * volume_factor = <strong>{product}</strong><br><br>
The 4*pi from Vol(S<sup>2</sup>) cancels the 1/(4*pi) from the loop measure.
</div>
""", unsafe_allow_html=True)

        eff = se.get("effective_correction")
        if eff is not None:
            st.metric("Effective correction", f"{eff}")

        assessment = se.get("assessment", "")
        if assessment:
            st.markdown(f"""
<div class="step-card">
{assessment}
</div>
""", unsafe_allow_html=True)

# --- E. Mass Correction ---
st.header("E. Mass Correction")

if summary:
    mc = summary.get("mass_correction")
    if mc:
        c1, c2, c3 = st.columns(3)
        c1.metric("Correction to mass", f"{mc.get('correction_to_mass', 'N/A')}")
        c2.metric("k_bare", f"{mc.get('k_bare', 'N/A')}")
        c3.metric("k_corrected", f"{mc.get('k_corrected', 'N/A')}")

        st.markdown(f"""
<div class="theorem-card">
<strong>Mass correction result:</strong><br><br>
m_eff = m_bare * (1 - alpha)<br><br>
k = sqrt(phi) * (1 - alpha)<br><br>
The one-loop diagram reproduces the (1 - alpha) factor that appears in the
Alpha Ladder prediction for G.
</div>
""", unsafe_allow_html=True)

# --- F. Consistency Check ---
st.header("F. Consistency Check")

if summary:
    con = summary.get("consistency")
    if con:
        c1, c2 = st.columns(2)
        c1.metric("G predicted", f"{float(con.get('G_predicted', 0)):.6e}")
        c2.metric("Residual", f"{con.get('residual_ppm', 0):+.2f} ppm")

        consistent = con.get("consistent", False)
        if consistent:
            st.markdown(f"""
<div class="theorem-card">
<strong>Consistency verified:</strong> The Feynman diagram calculation reproduces the
same G prediction as the algebraic derivation. The one-loop correction is exactly
(1 - alpha), confirming the S<sup>2</sup> volume cancellation mechanism.
</div>
""", unsafe_allow_html=True)
        else:
            assessment = con.get("assessment", "Consistency check failed.")
            st.markdown(f"""
<div class="warning-card">
<strong>Consistency issue:</strong> {assessment}
</div>
""", unsafe_allow_html=True)

# --- G. Scheme Dependence ---
st.header("G. Scheme Dependence")

if summary:
    sa = summary.get("scheme_analysis")
    if sa:
        indep = sa.get("scheme_independent_part", "")
        if indep:
            st.markdown(f"""
<div class="step-card">
<strong>Scheme-independent part:</strong><br><br>
{indep}
</div>
""", unsafe_allow_html=True)

        dep = sa.get("scheme_dependent_part", "")
        if dep:
            st.markdown(f"""
<div class="warning-card">
<strong>Scheme-dependent part:</strong><br><br>
{dep}
</div>
""", unsafe_allow_html=True)

        conclusion = sa.get("conclusion", "")
        if conclusion:
            st.markdown(f"""
<div class="step-card">
<strong>Conclusion:</strong><br><br>
{conclusion}
</div>
""", unsafe_allow_html=True)

# --- H. Honest Assessment ---
st.header("H. Honest Assessment")

if summary:
    ha = summary.get("honest_assessment", "")
    if ha:
        st.markdown(f"""
<div class="warning-card">
{ha}
</div>
""", unsafe_allow_html=True)

st.divider()
st.caption("Feynman Diagram Calculation | Alpha Ladder Research Dashboard")
