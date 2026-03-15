"""
Page 33 -- One-Alpha Derivation
Why (1-alpha) appears: S^2 volume cancellation.
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
    from alpha_ladder_core.one_alpha_derivation import summarize_one_alpha_derivation
    _core_available = True
except ImportError:
    pass


@st.cache_data
def _get_summary(_constants):
    if _core_available:
        return summarize_one_alpha_derivation(_constants)
    return None


st.title("The (1 - alpha) Derivation")
st.markdown("*Why bare alpha appears without pi factors: S^2 volume cancellation*")

summary = _get_summary(constants)

# --- A. The Cancellation ---
st.header("A. The Volume Cancellation")

if summary:
    vc = summary["volume_cancellation"]
    c1, c2, c3 = st.columns(3)
    c1.metric("Loop factor 1/(4pi)", f"{vc['loop_factor']:.6f}")
    c2.metric("Vol(S^2)/R^2", f"{vc['volume_factor']:.4f}")
    c3.metric("Product", f"{vc['product']:.1f}")

    st.markdown(f"""
<div class="theorem-card">
<strong>The key identity:</strong><br><br>
One-loop mass correction = alpha / (4*pi) * Vol(S^2) / R^2<br>
= alpha / (4*pi) * 4*pi<br>
= <strong>alpha</strong><br><br>
The 4*pi from the S^2 area exactly cancels the 1/(4*pi) from the 4D loop measure.
This is why (1 - alpha) appears with bare alpha, not alpha/pi or alpha/(2*pi).
</div>
""", unsafe_allow_html=True)
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# --- B. Candidate Mechanisms ---
st.header("B. Candidate Mechanisms")

if summary:
    mech = summary["mechanisms"]

    st.markdown(f"""
<div class="step-card">
{mech['assessment']}
</div>
""", unsafe_allow_html=True)

    rows = []
    for m in mech["mechanisms"]:
        rows.append({
            "Mechanism": m["name"],
            "Expression": m["expression"],
            "Ratio to alpha": f"{m['ratio_to_alpha']:.6f}",
            "Match": "YES" if m["matches_alpha"] else "",
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

# --- C. Uniqueness of S^2 ---
st.header("C. Why S^2 is Unique")

if summary:
    uniq = summary["uniqueness"]

    st.markdown(f"""
<div class="formula-card">
{uniq['assessment']}
</div>
""", unsafe_allow_html=True)

    rows = []
    for r in uniq["scan_results"]:
        rows.append({
            "Internal space": r["sphere"],
            "Vol/R^n": f"{r['volume_over_Rn']:.6f}",
            "Product with 1/(4pi)": f"{r['product_with_loop']:.6f}",
            "Gives bare alpha": "YES" if r["gives_bare_alpha"] else "",
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.markdown("""
<div class="step-card">
S^2 appears in three independent roles:<br>
1. <strong>Vacuum polynomial:</strong> x^2 + 6x + 4 = 0 requires n=2 for phi<br>
2. <strong>Volume cancellation:</strong> Vol(S^2)/R^2 = 4*pi cancels loop factor<br>
3. <strong>KK spectrum:</strong> S^2 gives the right degeneracies for the hierarchy<br><br>
All three point to the same geometry.
</div>
""", unsafe_allow_html=True)

# --- D. Sign Derivation ---
st.header("D. Why (1 - alpha), Not (1 + alpha)")

if summary:
    sign = summary.get("sign_analysis")
    if sign:
        c1, c2 = st.columns(2)
        c1.metric("Sign", sign["correction_form"])
        c2.metric("Moves toward data", "YES" if sign["numerical_cross_check"]["minus_moves_toward_data"] else "NO")

        st.markdown(f"""
<div class="step-card">
{sign['physical_argument']}
</div>
""", unsafe_allow_html=True)

        st.markdown(f"""
<div class="step-card">
<strong>Graviton propagator:</strong> {sign['graviton_propagator_sign']}
</div>
""", unsafe_allow_html=True)

# --- E. Degeneracy ---
st.header("E. Why Coefficient 1, Not 3")

if summary:
    deg = summary.get("degeneracy_analysis")
    if deg:
        c1, c2 = st.columns(2)
        c1.metric("l=0 degeneracy", deg["l0_degeneracy"])
        c2.metric("l=1 degeneracy", deg["l1_degeneracy"])

        st.markdown(f"""
<div class="step-card">
{deg['physical_argument']}
</div>
""", unsafe_allow_html=True)

        st.markdown(f"""
<div class="formula-card">
<strong>Why not 3?</strong> {deg['why_not_3']}
</div>
""", unsafe_allow_html=True)

        rows = []
        for entry in deg["kk_spectrum"]:
            rows.append({
                "l": entry["l"],
                "Degeneracy (2l+1)": entry["degeneracy"],
                "m^2/R^2": entry["mass_sq_over_R2"],
                "Sector": entry["sector"],
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

# --- F. Mode Sum & EM vs Gravitational ---
st.header("F. Mode Sum & Why Electromagnetic")

if summary:
    ms = summary.get("mode_sum")
    if ms:
        c1, c2, c3 = st.columns(3)
        c1.metric("Zeta(-1)", ms["spectral_zeta_minus1"])
        c2.metric("Polynomial identity", "Verified" if ms["polynomial_identity_verified"] else "Failed")
        em_gap = ms["gravitational_vs_em"]["orders_of_magnitude_gap"]
        c3.metric("EM/grav gap", f"{em_gap:.0f} orders")

        st.markdown(f"""
<div class="warning-card">
<strong>Critical finding:</strong> {ms['gravitational_vs_em']['conclusion']}
</div>
""", unsafe_allow_html=True)

        st.markdown(f"""
<div class="step-card">
{ms['assessment']}
</div>
""", unsafe_allow_html=True)

# --- G. The Derivation Chain ---
st.header("G. Complete Derivation Chain")

if summary:
    pred = summary["prediction"]

    c1, c2 = st.columns(2)
    c1.metric("G predicted", f"{float(pred['G_predicted']):.6e}")
    c2.metric("Residual", f"{pred['residual_ppm']:+.2f} ppm")

    for i, step in enumerate(pred["derivation_chain"]):
        st.markdown(f"""
<div class="step-card">
{step}
</div>
""", unsafe_allow_html=True)

# --- H. Key Finding ---
st.header("H. Key Finding")

if summary:
    st.markdown(f"""
<div class="theorem-card">
{summary['key_finding']}
</div>
""", unsafe_allow_html=True)

# --- I. Honest Assessment ---
st.header("I. Honest Assessment")

if summary:
    st.markdown(f"""
<div class="warning-card">
{summary['honest_assessment']}
</div>
""", unsafe_allow_html=True)

st.divider()
st.caption("One-Alpha Derivation | Alpha Ladder Research Dashboard")
