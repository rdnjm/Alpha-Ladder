"""
Page 37 -- Dimension Uniqueness
Why d=4, D=6: three independent constraints select a unique geometry.
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
    from alpha_ladder_core.dimension_uniqueness import summarize_dimension_uniqueness
    _core_available = True
except ImportError:
    pass


@st.cache_data
def _get_summary(_constants):
    if _core_available:
        return summarize_dimension_uniqueness(_constants)
    return None


st.title("Dimension Uniqueness")
st.markdown("*Why d=4, D=6: three constraints select a unique geometry*")

summary = _get_summary(constants)


def _fmt_ppm(val):
    """Format residual ppm: scientific notation for large values."""
    if val is None:
        return "N/A"
    if abs(val) > 10000:
        return f"{val:.1e}"
    return f"{val:+.1f}"


# ---------------------------------------------------------------------------
# A. Exponent Constraint
# ---------------------------------------------------------------------------
st.header("A. Exponent d*D = 24")

if summary:
    exp_data = summary.get("exponent_scan", {})
    scan_results = exp_data.get("scan_results", [])

    if scan_results:
        rows = []
        raw_residuals = []
        for row in scan_results:
            res = row.get("residual_ppm", None)
            raw_residuals.append(res)
            rows.append({
                "d": row.get("d", ""),
                "n": row.get("n", ""),
                "D": row.get("D", ""),
                "Exponent": row.get("exponent", ""),
                "Residual (ppm)": _fmt_ppm(res),
            })
        df = pd.DataFrame(rows)

        def _highlight_sub1000(row_series):
            idx = row_series.name
            res = raw_residuals[idx] if idx < len(raw_residuals) else None
            if res is not None and abs(res) < 1000:
                return ["background-color: #1a3a2e"] * len(row_series)
            return [""] * len(row_series)

        st.dataframe(
            df.style.apply(_highlight_sub1000, axis=1),
            use_container_width=True,
            hide_index=True,
        )

    assessment = exp_data.get("assessment", "")

    st.markdown("""
<div class="theorem-card">
<strong>Only two (d,n) pairs give exponent 24 and sub-1000 ppm residual: (3,5) and (4,2).</strong>
</div>
""", unsafe_allow_html=True)

    if assessment:
        st.markdown(f"""
<div class="step-card">
{assessment}
</div>
""", unsafe_allow_html=True)
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# B. Volume Cancellation Constraint
# ---------------------------------------------------------------------------
st.header("B. S^n Volume Cancellation: Only n=2")

if summary:
    vol_data = summary.get("volume_scan", {})
    scan_results = vol_data.get("scan_results", [])

    if scan_results:
        rows = []
        for row in scan_results:
            n_val = row.get("n", "")
            product = row.get("product", 0)
            rows.append({
                "n": n_val,
                "Sphere": f"S^{n_val}",
                "Volume / R^n": f"{row.get('vol_over_rn', 0):.6f}",
                "Product with 1/(4pi)": f"{product:.6f}",
                "Gives Bare Alpha": "YES" if abs(product - 1.0) < 1e-10 else "",
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.markdown("""
<div class="theorem-card">
<strong>Only S<sup>2</sup> produces exact volume cancellation. This selects n=2, eliminating (3,5).</strong>
</div>
""", unsafe_allow_html=True)

    assessment = vol_data.get("assessment", "")
    if assessment:
        st.markdown(f"""
<div class="step-card">
{assessment}
</div>
""", unsafe_allow_html=True)
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# C. Vacuum Polynomial Constraint
# ---------------------------------------------------------------------------
st.header("C. Vacuum Polynomial: Only d=4, D=6 Gives phi")

if summary:
    poly_data = summary.get("polynomial_scan", {})
    results = poly_data.get("scan_results", [])

    if results:
        rows = []
        for row in results:
            rows.append({
                "d": row.get("d", ""),
                "n": row.get("n", ""),
                "D": row.get("D", ""),
                "Discriminant": row.get("discriminant", ""),
                "Involves phi": row.get("involves_phi", ""),
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.markdown("""
<div class="theorem-card">
<strong>x<sup>2</sup>+6x+4=0 gives roots -3+/-sqrt(5), involving the golden ratio.
x<sup>2</sup>+8x+3=0 gives roots involving sqrt(13), not phi.</strong>
</div>
""", unsafe_allow_html=True)

    assessment = poly_data.get("assessment", "")
    if assessment:
        st.markdown(f"""
<div class="step-card">
{assessment}
</div>
""", unsafe_allow_html=True)
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# D. Uniqueness Proof
# ---------------------------------------------------------------------------
st.header("D. The Intersection: d=4, n=2, D=6")

if summary:
    proof = summary.get("uniqueness_proof", {})

    c1, c2, c3 = st.columns(3)
    c1.metric("Exponent pairs", str(proof.get("exponent_constraint", "")))
    c2.metric("Volume n", str(proof.get("volume_constraint", "")))
    c3.metric("Polynomial pairs", str(proof.get("polynomial_constraint", "")))

    unique_pair = proof.get("unique_pair", {})
    st.markdown(f"""
<div class="proof-card">
<strong>The intersection of three independent constraints is exactly ONE point:
d={unique_pair.get('d', 4)}, n={unique_pair.get('n', 2)}, D={unique_pair.get('D', 6)}.</strong>
</div>
""", unsafe_allow_html=True)
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# E. c3 = phi/2
# ---------------------------------------------------------------------------
st.header("E. The Third Coefficient: c3 = phi/2")

if summary:
    c3_data = summary.get("c3_derivation", {})

    c1, c2, c3_col = st.columns(3)
    c1.metric("c3 (phi/2)", f"{c3_data.get('c3_phi_half', 0):.6f}")
    c2.metric("c3 exact", f"{c3_data.get('c3_exact', 0):.6f}")
    c3_col.metric("Residual (%)", f"{c3_data.get('c3_residual_percent', 0):.4f}")

    st.markdown("""
<div class="formula-card">
<strong>F = 1 + (d-1)*alpha<sup>2</sup> + (phi/2)*alpha<sup>3</sup>
= 1 + 3*alpha<sup>2</sup> + 0.809*alpha<sup>3</sup></strong>
</div>
""", unsafe_allow_html=True)

    st.markdown("""
<div class="step-card">
<strong>Both coefficients derive from geometry:</strong> c2 = d-1 = 3, c3 = phi/2.
</div>
""", unsafe_allow_html=True)
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# F. Mu Prediction
# ---------------------------------------------------------------------------
st.header("F. Second Prediction: mu from alpha and phi")

if summary:
    mu_data = summary.get("mu_prediction", {})

    c1, c2, c3_col, c4 = st.columns(4)
    c1.metric("mu predicted", f"{mu_data.get('mu_predicted', 0):.6f}")
    c2.metric("mu measured", f"{mu_data.get('mu_measured', 0):.6f}")
    residual_ppm = mu_data.get("residual_ppm", 0)
    c3_col.metric("Residual (ppm)", f"{residual_ppm:+.2f}")
    sigma = mu_data.get("sigma_tension", 0)
    c4.metric("Sigma tension", f"{sigma:.2f}")

    st.markdown(f"""
<div class="theorem-card">
<strong>mu predicted to {residual_ppm:+.2f} ppm ({sigma:.2f} sigma). At measurement precision.</strong>
</div>
""", unsafe_allow_html=True)
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# G. G Prediction
# ---------------------------------------------------------------------------
st.header("G. Primary Prediction: G from alpha, mu, phi")

if summary:
    g_data = summary.get("G_prediction", {})

    c1, c2, c3_col = st.columns(3)
    c1.metric("G unified", f"{float(g_data.get('G_unified', 0)):.6e}")
    c2.metric("G measured", f"{float(g_data.get('G_measured', 0)):.6e}")
    g_residual = g_data.get("residual_unified_ppm", 0)
    c3_col.metric("Residual (ppm)", f"{g_residual:+.1f}")

    st.markdown(f"""
<div class="theorem-card">
<strong>G predicted to {g_residual:+.1f} ppm with zero fitted parameters.</strong>
</div>
""", unsafe_allow_html=True)
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# H. The Complete Formula
# ---------------------------------------------------------------------------
st.header("H. The Complete Derived Formula")

if summary:
    formula_data = summary.get("complete_formula", {})

    components = formula_data.get("formula_components", [])
    if components:
        rows = []
        for comp in components:
            rows.append({
                "Component": comp.get("component", ""),
                "Value": comp.get("value", ""),
                "Origin": comp.get("origin", ""),
                "Derived": comp.get("derived", ""),
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    n_fitted = formula_data.get("n_fitted_params", "N/A")
    st.metric("Fitted parameters", str(n_fitted))

    g_pred = formula_data.get("G_predicted", 0)
    f_residual = formula_data.get("residual_ppm", 0)
    st.markdown(f"""
<div class="proof-card">
<strong>Complete formula:</strong><br><br>
G = phi<sup>2</sup>/2 * (1 + 3*alpha<sup>2</sup> + (phi/2)*alpha<sup>3</sup>) * alpha<sup>21</sup> * hbar*c / m_e<sup>2</sup><br><br>
G = {g_pred:.6e} ({f_residual:+.2f} ppm from measured)
</div>
""", unsafe_allow_html=True)
else:
    st.info("Core module not available. Install alpha_ladder_core.")

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
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Dimension Uniqueness | Alpha Ladder Research Dashboard")
