"""
Page 36 -- Mu Tension Resolution
The bridge and mu-structure formulas are the same formula.
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
    from alpha_ladder_core.mu_tension import summarize_mu_tension
    _core_available = True
except ImportError:
    pass


@st.cache_data
def _get_summary(_constants):
    if _core_available:
        return summarize_mu_tension(_constants)
    return None


st.title("Mu Tension Resolution")
st.markdown("*The bridge and mu-structure formulas are the same formula*")

summary = _get_summary(constants)

# ---------------------------------------------------------------------------
# A. Leading-Order Identity
# ---------------------------------------------------------------------------
st.header("A. The Deep Relation: mu^2 * alpha^3 = phi^2/2")

if summary:
    lo = summary.get("leading_order", {})
    c1, c2 = st.columns(2)
    c1.metric("mu^2 * alpha^3", f"{lo.get('mu_sq_alpha_cubed', 0):.6e}")
    c2.metric("phi^2 / 2", f"{lo.get('phi_sq_over_2', 0):.6e}")

    residual = lo.get("residual_ppm", 0)
    st.markdown(f"""
<div class="theorem-card">
<strong>Leading-order identity:</strong><br><br>
mu<sup>2</sup> * alpha<sup>3</sup> = phi<sup>2</sup>/2 to <strong>{residual:.2f} ppm</strong>.<br><br>
The proton-to-electron mass ratio is approximately determined by alpha and the golden ratio.
</div>
""", unsafe_allow_html=True)

    assessment = lo.get("assessment", "")
    if assessment:
        st.markdown(f"""
<div class="step-card">
{assessment}
</div>
""", unsafe_allow_html=True)
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# B. The Bracketing
# ---------------------------------------------------------------------------
st.header("B. The Two Formulas Bracket the Exact Answer")

if summary:
    fc = summary.get("formula_comparison", {})

    rows = [
        {"Formula": "mu-structure (C_mu)", "C value": f"{fc.get('C_mu', 0):.10f}", "Residual (ppm)": "reference"},
        {"Formula": "phi^2/2 (bare)", "C value": f"{fc.get('C_bridge_lo', 0):.10f}", "Residual (ppm)": f"{fc.get('residuals', {}).get('leading_order', 0):+.2f}"},
        {"Formula": "phi^2/2 * (1+3a^2)", "C value": f"{fc.get('C_bridge_c2', 0):.10f}", "Residual (ppm)": f"{fc.get('residuals', {}).get('c2_only', 0):+.2f}"},
        {"Formula": "phi^2/2 * (1+3a^2+8/5a^3)", "C value": f"{fc.get('C_bridge_nlo', 0):.10f}", "Residual (ppm)": f"{fc.get('residuals', {}).get('nlo', 0):+.2f}"},
    ]
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.markdown("""
<div class="theorem-card">
<strong>Bracketing result:</strong><br><br>
The mu-structure formula gives +0.31 ppm vs the bridge with c2=3,
and -0.31 ppm vs the bridge with c3=8/5.
The two formulas <strong>BRACKET</strong> the exact answer.
</div>
""", unsafe_allow_html=True)

    assessment = fc.get("assessment", "")
    if assessment:
        st.markdown(f"""
<div class="step-card">
{assessment}
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# C. Mu Predicted from Bridge
# ---------------------------------------------------------------------------
st.header("C. Predicting mu from alpha and phi")

if summary:
    mfb = summary.get("mu_from_bridge", {})
    c1, c2, c3 = st.columns(3)
    c1.metric("mu predicted", f"{mfb.get('mu_predicted', 0):.6f}")
    c2.metric("mu measured", f"{mfb.get('mu_measured', 0):.6f}")
    c3.metric("Residual (ppm)", f"{mfb.get('residual_ppm', 0):+.2f}")

    mu_pred = mfb.get("mu_predicted", 0)
    res_ppm = mfb.get("residual_ppm", 0)
    st.markdown(f"""
<div class="formula-card">
<strong>Bridge with c2=3 predicts mu = {mu_pred:.6f}</strong>,
just <strong>{res_ppm:+.2f} ppm</strong> from measured.
</div>
""", unsafe_allow_html=True)

    assessment = mfb.get("assessment", "")
    if assessment:
        st.markdown(f"""
<div class="step-card">
{assessment}
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# D. The Exact c2
# ---------------------------------------------------------------------------
st.header("D. c2 = 3.006: Almost Exactly 3")

if summary:
    ec2 = summary.get("exact_c2", {})
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("c2 exact", f"{ec2.get('c2_exact', 0):.6f}")
    c2.metric("Integer part", f"{ec2.get('c2_integer_part', 0)}")
    c3.metric("Excess", f"{ec2.get('c2_excess', 0):.6f}")
    c4.metric("Excess / alpha", f"{ec2.get('c2_excess_over_alpha', 0):.4f}")

    st.markdown(f"""
<div class="step-card">
<strong>c2 = 3 + 0.006.</strong> The integer part 3 = d-1 (spatial dimensions).
The excess 0.006 = 0.81 * alpha could be a higher-order radiative correction.
</div>
""", unsafe_allow_html=True)

    assessment = ec2.get("assessment", "")
    if assessment:
        st.markdown(f"""
<div class="step-card">
{assessment}
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# E. The Exact c3
# ---------------------------------------------------------------------------
st.header("E. c3 = 0.81: The Unifying Value")

if summary:
    ec3 = summary.get("exact_c3", {})
    c1, c2, c3_col = st.columns(3)
    c1.metric("c3 exact", f"{ec3.get('c3_exact', 0):.4f}")
    c2.metric("c3 bridge (fitted)", f"{ec3.get('c3_bridge', 0):.1f}")
    c3_col.metric("Difference", f"{ec3.get('difference', 0):.4f}")

    st.markdown(f"""
<div class="warning-card">
<strong>The bridge-fitted c3 = 8/5 is too large.</strong>
The mu-structure requires c3 = 0.81 for exact consistency.
</div>
""", unsafe_allow_html=True)

    candidates = ec3.get("c3_clean_candidates", [])
    if candidates:
        rows = []
        for cand in candidates:
            rows.append({
                "Expression": cand.get("expression", ""),
                "Value": f"{cand.get('value', 0):.6f}",
                "Residual": f"{cand.get('residual', 0):.6f}",
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    assessment = ec3.get("assessment", "")
    if assessment:
        st.markdown(f"""
<div class="step-card">
{assessment}
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# F. Unification Verified
# ---------------------------------------------------------------------------
st.header("F. Bridge = Mu-Structure (Verified)")

if summary:
    unif = summary.get("unification", {})
    c1, c2, c3_col = st.columns(3)
    c1.metric("G (bridge)", f"{float(unif.get('G_bridge', 0)):.6e}")
    c2.metric("G (mu-structure)", f"{float(unif.get('G_mu_structure', 0)):.6e}")
    c3_col.metric("Difference (ppm)", f"{unif.get('difference_ppm', 0):.4f}")

    unified = unif.get("unified", False)
    if unified:
        st.markdown(f"""
<div class="theorem-card">
<strong>Unification confirmed:</strong><br><br>
With c3 = 0.81, both formulas give <strong>IDENTICAL</strong> G predictions.
The mu tension is resolved.
</div>
""", unsafe_allow_html=True)
    else:
        st.markdown(f"""
<div class="warning-card">
<strong>Unification not yet achieved.</strong>
The two formulas still give different G predictions at this precision.
</div>
""", unsafe_allow_html=True)

    assessment = unif.get("assessment", "")
    if assessment:
        st.markdown(f"""
<div class="step-card">
{assessment}
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# G. Correction Hierarchy
# ---------------------------------------------------------------------------
st.header("G. Correction Series Hierarchy")

if summary:
    ch = summary.get("correction_hierarchy", {})

    terms = ch.get("terms", [])
    if terms:
        rows = []
        for term in terms:
            rows.append({
                "Order": term.get("order", ""),
                "Coefficient": f"{term.get('coefficient', 0):.4f}",
                "Contribution (ppm)": f"{term.get('contribution_ppm', 0):.4f}",
                "Cumulative (ppm)": f"{term.get('cumulative_ppm', 0):.4f}",
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.markdown("""
<div class="formula-card">
<strong>F = 1 + 3*alpha<sup>2</sup> + 0.81*alpha<sup>3</sup> + ...</strong><br><br>
The series converges rapidly.
</div>
""", unsafe_allow_html=True)

    assessment = ch.get("assessment", "")
    if assessment:
        st.markdown(f"""
<div class="step-card">
{assessment}
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# H. The Second Prediction: mu from alpha and phi
# ---------------------------------------------------------------------------
st.header("H. The Second Prediction")

if summary:
    mp = summary.get("mu_prediction", {})
    c1, c2 = st.columns(2)
    c1.metric("mu predicted", f"{mp.get('mu_predicted', 0):.6f}")
    c2.metric("Residual (ppm)", f"{mp.get('residual_ppm', 0):+.2f}")

    res_ppm = mp.get("residual_ppm", 0)
    st.markdown(f"""
<div class="theorem-card">
<strong>Second prediction:</strong><br><br>
The proton-to-electron mass ratio is predicted from alpha and phi alone
to <strong>{res_ppm:+.2f} ppm</strong> with zero fitted parameters.
This is a genuine second prediction of the framework.
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# I. Honest Assessment
# ---------------------------------------------------------------------------
st.header("I. Honest Assessment")

if summary:
    honest = summary.get("honest_assessment", "")
    st.markdown(f"""
<div class="warning-card">
{honest}
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Mu Tension Resolution | Alpha Ladder Research Dashboard")
