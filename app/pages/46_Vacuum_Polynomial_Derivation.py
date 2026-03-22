"""
Page 46 -- Vacuum Polynomial Derivation: Can x^2+6x+4=0 Be Derived?
Documents the systematic FAILURE to derive the vacuum polynomial from first principles.
All results are negative. This is an honest assessment.
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
    from alpha_ladder_core.vacuum_polynomial_derivation import (
        summarize_vacuum_polynomial_derivation,
    )
    _core_available = True
except ImportError:
    pass


@st.cache_data
def _get_summary(_constants):
    if _core_available:
        return summarize_vacuum_polynomial_derivation(_constants)
    return None


# ---------------------------------------------------------------------------
# Title and Introduction
# ---------------------------------------------------------------------------
st.title("Vacuum Polynomial: Can x\u00b2+6x+4=0 Be Derived?")

st.markdown("""
<div class="warning-card">
<strong>The vacuum polynomial is the WEAKEST LINK in the Alpha Ladder.</strong>
Three systematic approaches were tried. All failed.
The polynomial x<sup>2</sup>+Dx+d=0 with (D=6, d=4) produces the golden ratio
and hence C<sub>0</sub> = phi<sup>2</sup>/2, but it is postulated, not derived.
</div>
""", unsafe_allow_html=True)

if not _core_available:
    st.error(
        "Core module `alpha_ladder_core.vacuum_polynomial_derivation` not found. "
        "Install alpha_ladder_core to enable computations on this page."
    )

summary = _get_summary(constants)

# ---------------------------------------------------------------------------
# Section 1: Algebraic Coefficient Search
# ---------------------------------------------------------------------------
st.header("1. Algebraic Coefficient Search")

if summary:
    coeff = summary.get("coefficient_scan", {})
    phi_producing = coeff.get("phi_producing", [])

    c1, c2, c3 = st.columns(3)
    c1.metric("Polynomials scanned", str(coeff.get("polynomials_scanned", "N/A")))
    c2.metric("Phi-producing", str(len(phi_producing)))
    c3.metric("Unique to (4,2)", str(coeff.get("unique_to_d4_n2", "N/A")))

    if phi_producing:
        rows = []
        for entry in phi_producing:
            roots = entry.get("roots", [0, 0])
            rows.append({
                "Polynomial": entry.get("polynomial", ""),
                "Source A": entry.get("source_A", ""),
                "Source B": entry.get("source_B", ""),
                "Disc": str(entry.get("disc", "")),
                "Root+": f"{roots[0]:.6f}",
                "Root-": f"{roots[1]:.6f}",
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.markdown("""
<div class="step-card">
<strong>Result:</strong> Only x<sup>2</sup>+6x+4 (A=D=6, B=d=4) produces phi from
KK dimension numbers. The selection is numerological: coefficients are the dimension
numbers plugged into a quadratic template.
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from coefficient scan.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 2: Dimension Pair Scan
# ---------------------------------------------------------------------------
st.header("2. Dimension Pair Scan")

if summary:
    dim_scan = summary.get("dimension_scan", {})
    phi_pairs = dim_scan.get("phi_pairs", [])

    c1, c2, c3 = st.columns(3)
    c1.metric("Pairs scanned", str(dim_scan.get("pairs_scanned", "N/A")))
    c2.metric("Minimal pair", str(dim_scan.get("minimal_pair", "N/A")))
    c3.metric("Minimal physical", str(dim_scan.get("minimal_physical_pair", "N/A")))

    if phi_pairs:
        rows = []
        for p in phi_pairs:
            rows.append({
                "d": str(p.get("d", "")),
                "n": str(p.get("n", "")),
                "D = d+n": str(p.get("D", "")),
                "Disc": str(p.get("disc", "")),
                "k": str(p.get("k", "")),
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.markdown("""
<div class="theorem-card">
<strong>Finding:</strong> (4,2) gives disc=20 (minimal physical case). 7+ other pairs
also give sqrt(5). The dimension pair is not unique.
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from dimension pair scan.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 3: Physical Quantities as Roots
# ---------------------------------------------------------------------------
st.header("3. Physical Quantities as Roots")

if summary:
    kk = summary.get("kk_root_check", {})
    quantities = kk.get("quantities_tested", [])
    roots = kk.get("roots", {})

    c1, c2, c3 = st.columns(3)
    c1.metric("Quantities tested", str(len(quantities)))
    c2.metric("Root+", f"{roots.get('r_plus', 0):.6f}")
    c3.metric("Root-", f"{roots.get('r_minus', 0):.6f}")

    if quantities:
        rows = []
        for q in quantities:
            rows.append({
                "Quantity": q.get("name", ""),
                "Value": f"{q.get('value', 0):.6e}",
                "P(x)": f"{q.get('polynomial_value', 0):.6e}",
                "Is root": str(q.get("is_root", False)),
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.markdown("""
<div class="warning-card">
<strong>Both roots are negative</strong> (-0.764, -5.236). All physical KK quantities
are non-negative. No physical quantity satisfies the polynomial.
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from KK root check.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 4: Moduli Space Geometry
# ---------------------------------------------------------------------------
st.header("4. Moduli Space Geometry")

if summary:
    moduli = summary.get("moduli_geometry", {})

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Moduli dimension", str(moduli.get("moduli_dimension", "N/A")))
    c2.metric("Curvature", str(moduli.get("curvature", "N/A")))
    mc = moduli.get("metric_constant")
    c3.metric("Metric constant", f"{mc:.6f}" if mc is not None else "N/A")
    c4.metric("SS stationarity deg.", str(moduli.get("ss_stationarity_degree", "N/A")))

    st.markdown("""
<div class="step-card">
<strong>Result:</strong> The moduli space is 1D and flat. No polynomial emerges from
its geometry. The Salam-Sezgin stationarity condition is cubic, not quadratic.
The mass matrix is 1x1 so its characteristic polynomial is degree 1.
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from moduli space geometry check.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 5: Swampland Constraints
# ---------------------------------------------------------------------------
st.header("5. Swampland Constraints")

if summary:
    swamp = summary.get("swampland_constraints", {})
    results = swamp.get("results", [])

    c1, c2, c3 = st.columns(3)
    c1.metric("Conjectures checked", str(swamp.get("conjectures_checked", "N/A")))
    c2.metric("Any constrains C_0", str(swamp.get("any_constrains_C0", "N/A")))
    c3.metric("dS tension", str(swamp.get("ds_tension", "N/A")))

    if results:
        rows = []
        for r in results:
            rows.append({
                "Conjecture": r.get("name", ""),
                "Statement": r.get("statement", ""),
                "Constrains C_0": str(r.get("constrains_C0", False)),
                "Status": r.get("status", ""),
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.markdown("""
<div class="warning-card">
<strong>No Swampland conjecture constrains C<sub>0</sub>.</strong> The dS conjecture
conflicts with the Salam-Sezgin vacuum (V'=0 at a dS minimum), but this tension
concerns the vacuum structure, not C<sub>0</sub> specifically.
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from Swampland constraint check.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 6: Alternative Polynomial x^2+3x+1
# ---------------------------------------------------------------------------
st.header("6. Alternative Polynomial x\u00b2+3x+1")

if summary:
    alt = summary.get("alternative_polynomial", {})
    alt_roots = alt.get("roots", {})

    c1, c2, c3 = st.columns(3)
    c1.metric("Disc", str(alt.get("disc", "N/A")))
    c2.metric("Root r_1 = -1/phi", f"{alt_roots.get('r1', 0):.6f}")
    c3.metric("Root r_2 = -phi^2", f"{alt_roots.get('r2', 0):.6f}")

    coeffs = alt.get("coefficients_from_kk", {})
    st.markdown(f"""
<div class="formula-card">
<strong>x<sup>2</sup>+3x+1=0</strong> has roots -1/phi and -phi<sup>2</sup>.<br><br>
Coefficients: A = {coeffs.get('A', 'N/A')}, B = {coeffs.get('B', 'N/A')}<br>
C<sub>0</sub> from roots: {alt.get('relation_to_bridge', 'N/A')}<br>
C<sub>0</sub> match: {alt.get('C0_match', 'N/A')}<br>
More natural: {alt.get('is_more_natural', 'N/A')} (arguable -- A=dim(SO(3)) is geometric)
</div>
""", unsafe_allow_html=True)

    honest = alt.get("honest_assessment", "")
    if honest:
        st.markdown(f"""
<div class="step-card">
<strong>Assessment:</strong> {honest}
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from alternative polynomial analysis.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 7: Honest Assessment
# ---------------------------------------------------------------------------
st.header("7. Honest Assessment")

if summary:
    status = summary.get("honest_status", {})

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Derived", str(status.get("is_derived", "N/A")))
    c2.metric("Ansatz", str(status.get("is_ansatz", "N/A")))
    c3.metric("Alternatives exist", str(status.get("alternatives_exist", "N/A")))
    c4.metric("Best alternative", str(status.get("best_alternative", "N/A")))

    st.markdown(f"""
<div class="proof-card">
<strong>is_derived = {status.get('is_derived', 'N/A')}</strong><br>
<strong>is_ansatz = {status.get('is_ansatz', 'N/A')}</strong><br><br>
Status: {status.get('status', 'N/A')}
</div>
""", unsafe_allow_html=True)

    honest = status.get("honest_assessment", "")
    if honest:
        st.markdown(f"""
<div class="warning-card">
<strong>The vacuum polynomial remains a phenomenological ansatz.</strong> It is the
weakest link in the derivation chain. Everything golden-ratio-related
(C<sub>0</sub>=phi<sup>2</sup>/2, c<sub>3</sub>=phi/2) depends on it.<br><br>
{honest}
</div>
""", unsafe_allow_html=True)

    # What works / what fails summary
    what_works = summary.get("what_works", [])
    what_fails = summary.get("what_fails", [])

    if what_works or what_fails:
        st.subheader("Summary")
        col_w, col_f = st.columns(2)
        with col_w:
            st.markdown("**What works:**")
            for item in what_works:
                st.markdown(f"- {item}")
        with col_f:
            st.markdown("**What fails:**")
            for item in what_fails:
                st.markdown(f"- {item}")

    remaining = summary.get("remaining_gap", "")
    if remaining:
        st.markdown(f"""
<div class="theorem-card">
<strong>Remaining gap:</strong> {remaining}
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from honest status assessment.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Vacuum Polynomial Derivation | Alpha Ladder Research Dashboard")
