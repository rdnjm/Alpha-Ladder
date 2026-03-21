"""
Page 41 -- Coefficient Derivation: c_2 = 3 and c_3 = phi/2
Systematic search for first-principles origins of the bridge correction coefficients.
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
    from alpha_ladder_core.c2_derivation import summarize_c2_derivation
    _core_available = True
except ImportError:
    pass


@st.cache_data
def _get_summary(_constants):
    if _core_available:
        return summarize_c2_derivation(_constants)
    return None


# ---------------------------------------------------------------------------
# Title and Introduction
# ---------------------------------------------------------------------------
st.title("Coefficient Derivation: c_2 = 3 and c_3 = phi/2")
st.markdown("*Systematic search for first-principles origins of the bridge correction coefficients*")

st.markdown("""
<div class="warning-card">
<strong>Gap 3 (c_2=3) status: Partial.</strong> Gap 4 (c_3=phi/2) status: Open.
No complete derivation exists.
</div>
""", unsafe_allow_html=True)

if not _core_available:
    st.error(
        "Core module `alpha_ladder_core.c2_derivation` not found. "
        "Install alpha_ladder_core to enable computations on this page."
    )

summary = _get_summary(constants)

# ---------------------------------------------------------------------------
# Section 1: Heat Kernel: 1/a_1 = 3
# ---------------------------------------------------------------------------
st.header("1. Heat Kernel: 1/a_1 = 3")

if summary and summary.get("heat_kernel"):
    hk = summary["heat_kernel"]
    c1, c2, c3 = st.columns(3)
    c1.metric("a_1 density", f"{hk.get('a1_density', 'N/A')}")
    c2.metric("1 / a_1", f"{hk.get('inverse_a1', 'N/A')}")
    c3.metric("Matches c_2", str(hk.get("matches_c2", "N/A")))

    st.markdown("""
<div class="formula-card">
<strong>Heat kernel coefficient:</strong> a_1 = R/6 = (2/a_0^2) / 6 = 1/3 on the unit S^2
(where R = 2 for the round unit sphere).<br><br>
The inverse 1/a_1 = 3 matches c_2 exactly. This is a spectral invariant of the
Laplacian on S^2 and requires no fitted parameters.
</div>
""", unsafe_allow_html=True)
elif summary:
    st.info("Heat kernel data not available in summary.")
else:
    c1, c2, c3 = st.columns(3)
    c1.metric("a_1 density", "1/3")
    c2.metric("1 / a_1", "3")
    c3.metric("Matches c_2", "True")
    st.markdown("""
<div class="formula-card">
<strong>Heat kernel coefficient:</strong> a_1 = R/6 = (2/a_0^2) / 6 = 1/3 on the unit S^2
(where R = 2 for the round unit sphere).<br><br>
The inverse 1/a_1 = 3 matches c_2 exactly. This is a spectral invariant of the
Laplacian on S^2 and requires no fitted parameters.
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Section 2: Holomorphic Euler Characteristic
# ---------------------------------------------------------------------------
st.header("2. Holomorphic Euler Characteristic: chi(T_{S^2}) = 3")

if summary and summary.get("chi_tangent_bundle"):
    chi = summary["chi_tangent_bundle"]
    c1, c2, c3 = st.columns(3)
    c1.metric("chi(T)", f"{chi.get('chi_T', 'N/A')}")
    c2.metric("h^0", f"{chi.get('h0', 'N/A')}")
    c3.metric("h^1", f"{chi.get('h1', 'N/A')}")
elif summary:
    st.info("Holomorphic Euler characteristic data not available in summary.")
else:
    c1, c2, c3 = st.columns(3)
    c1.metric("chi(T)", "3")
    c2.metric("h^0", "3")
    c3.metric("h^1", "0")

st.markdown("""
<div class="theorem-card">
<strong>Theorem:</strong> chi(O(2)) = 3 via Hirzebruch-Riemann-Roch. The 3 holomorphic
sections {1, z, z^2} of O(2) = T_{CP^1} are the 3 Killing vectors generating SO(3).
S^2 is the only sphere that is a complex manifold.
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Section 3: Gauge Matching Identity
# ---------------------------------------------------------------------------
st.header("3. Gauge Matching Identity")

if summary and summary.get("gauge_matching"):
    gm = summary["gauge_matching"]
    c1, c2, c3 = st.columns(3)
    c1.metric("g_KK^2", f"{gm.get('g_KK_squared', 'N/A')}")
    c2.metric("Ratio", f"{gm.get('ratio', 'N/A')}")
    c3.metric("Residual (ppm)", f"{gm.get('residual_ppm', 'N/A')}")
elif summary:
    st.info("Gauge matching data not available in summary.")
else:
    c1, c2, c3 = st.columns(3)
    c1.metric("g_KK^2", "N/A")
    c2.metric("Ratio", "N/A")
    c3.metric("Residual (ppm)", "N/A")

st.markdown("""
<div class="proof-card">
<strong>Proof:</strong> g_KK^4 / (16 pi^2) = alpha^2 exactly. All factors of pi cancel
via the gauge matching condition g_KK^2 = 4 pi alpha.
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Section 4: Sphere Scan: Uniqueness of n=2
# ---------------------------------------------------------------------------
st.header("4. Sphere Scan: Uniqueness of n=2")

if summary and summary.get("sphere_scan"):
    ss = summary["sphere_scan"]
    scan_results = ss.get("scan_results", [])
    if scan_results:
        rows = []
        for entry in scan_results:
            rows.append({
                "n": entry.get("n", ""),
                "1/a_1": entry.get("inverse_a1", ""),
                "dim(SO(n+1))": entry.get("dim_SO", ""),
                "chi+1": entry.get("chi_plus_1", ""),
                "n+1": entry.get("n_plus_1", ""),
                "d-1": entry.get("d_minus_1", ""),
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
    else:
        st.info("No sphere scan results in summary.")
elif summary:
    st.info("Sphere scan data not available in summary.")
else:
    # Fallback: show known results for n=1..5
    fallback_rows = []
    for n in range(1, 6):
        inv_a1 = round(n * (n + 1) / 2, 1) if n > 0 else "N/A"
        dim_so = n * (n + 1) // 2
        fallback_rows.append({
            "n": n,
            "1/a_1": inv_a1,
            "dim(SO(n+1))": dim_so,
            "chi+1": n + 1 if n % 2 == 0 else 1,
            "n+1": n + 1,
            "d-1": n + 2 - 1,
        })
    st.dataframe(pd.DataFrame(fallback_rows), use_container_width=True, hide_index=True)

st.markdown("""
<div class="theorem-card">
<strong>Uniqueness:</strong> Five quantities all equal 3 ONLY at n=2.
This degeneracy is unique to S^2.
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Section 5: One-Loop Attempts
# ---------------------------------------------------------------------------
st.header("5. One-Loop Attempts")

if summary and summary.get("loop_attempts"):
    la = summary["loop_attempts"]
    attempts = la.get("attempts", [])
    if attempts:
        rows = []
        for entry in attempts:
            rows.append({
                "Name": entry.get("name", ""),
                "Result": entry.get("result", ""),
                "Target": entry.get("target", ""),
                "Status": entry.get("status", ""),
                "Gap": entry.get("gap", ""),
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
    else:
        st.info("No loop attempt data in summary.")
elif summary:
    st.info("Loop attempts data not available in summary.")
else:
    st.dataframe(pd.DataFrame([{
        "Name": "N/A", "Result": "N/A", "Target": "3",
        "Status": "N/A", "Gap": "N/A",
    }]), use_container_width=True, hide_index=True)

st.markdown("""
<div class="warning-card">
<strong>No explicit one-loop calculation reproduces c_2 = 3 exactly.</strong>
The mechanism is identified but the coefficient is not derived.
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Section 6: c_3 = phi/2 from Vacuum Polynomial
# ---------------------------------------------------------------------------
st.header("6. c_3 = phi/2 from Vacuum Polynomial")

if summary and summary.get("c3_vacuum_polynomial"):
    vp = summary["c3_vacuum_polynomial"]
    c1, c2, c3 = st.columns(3)
    c1.metric("c_3 value", f"{vp.get('c3_value', 'N/A')}")
    c2.metric("C_0", f"{vp.get('C_0', 'N/A')}")
    c3.metric("c_3 / C_0", f"{vp.get('ratio_c3_to_C0', 'N/A')}")
elif summary:
    st.info("Vacuum polynomial data not available in summary.")
else:
    c1, c2, c3 = st.columns(3)
    c1.metric("c_3 value", "phi/2")
    c2.metric("C_0", "phi^2/2")
    c3.metric("c_3 / C_0", "1/phi")

st.markdown("""
<div class="formula-card">
<strong>Vacuum polynomial origin:</strong> c_3 = (r_+ + d) / d where r_+ = -3 + sqrt(5)
is the root of x^2 + 6x + 4 = 0 and d = 4. This is the same vacuum polynomial
that produces C_0 = phi^2/2.
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Section 7: Pentagonal Polynomial
# ---------------------------------------------------------------------------
st.header("7. Pentagonal Polynomial")

if summary and summary.get("c3_pentagonal"):
    pent = summary["c3_pentagonal"]

    st.markdown("""
<div class="formula-card">
<strong>Pentagonal connection:</strong> The substitution f = (r+d)/d maps x^2+6x+4=0
to 4f^2-2f-1=0, the minimal polynomial of cos(pi/5). This is a mathematical
consequence of the vacuum polynomial ansatz, not new numerology.
</div>
""", unsafe_allow_html=True)
else:
    st.markdown("""
<div class="formula-card">
<strong>Pentagonal connection:</strong> The substitution f = (r+d)/d maps x^2+6x+4=0
to 4f^2-2f-1=0, the minimal polynomial of cos(pi/5). This is a mathematical
consequence of the vacuum polynomial ansatz, not new numerology.
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="step-card">
<strong>WARNING:</strong> Although c_3 = C_0/phi provides a natural ratio, the geometric
resummation that extended this into an infinite 1/phi series has been retracted
(>14 sigma tension with Alighanbari et al. 2025). The identity does NOT license
higher-order extrapolation.
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Section 8: Rationality No-Go
# ---------------------------------------------------------------------------
st.header("8. Rationality No-Go")

if summary and summary.get("c3_rationality_nogo"):
    nogo = summary["c3_rationality_nogo"]

st.markdown("""
<div class="proof-card">
<strong>No-go theorem:</strong> phi/2 is irrational. All S^2 spectral data (eigenvalues,
degeneracies, Seeley-DeWitt coefficients, zeta values) are rational. Therefore phi/2
CANNOT come from S^2 spectral geometry. It must enter from the vacuum polynomial.
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Section 9: Honest Assessment
# ---------------------------------------------------------------------------
st.header("9. Honest Assessment")

if summary and summary.get("honest_assessment"):
    ha = summary["honest_assessment"]

    # Warning card with full text
    assessment_text = ha.get("text", "")
    if assessment_text:
        st.markdown(f"""
<div class="warning-card">
<strong>Honest Assessment:</strong><br><br>
{assessment_text}
</div>
""", unsafe_allow_html=True)

    # Key metrics
    c1, c2 = st.columns(2)
    c1.metric("Derivation achieved", str(ha.get("derivation_achieved", False)))
    c2.metric("Best candidate", ha.get("best_candidate", "N/A"))

    # Two columns: What Works vs What Fails
    col_works, col_fails = st.columns(2)

    with col_works:
        st.subheader("What Works")
        what_works = ha.get("what_works", [])
        if what_works:
            for item in what_works:
                st.markdown(f"- {item}")
        else:
            st.markdown("- No data available")

    with col_fails:
        st.subheader("What Fails")
        what_fails = ha.get("what_fails", [])
        if what_fails:
            for item in what_fails:
                st.markdown(f"- {item}")
        else:
            st.markdown("- No data available")
elif summary:
    st.info("Honest assessment data not available in summary.")
else:
    st.markdown("""
<div class="warning-card">
<strong>Honest Assessment:</strong><br><br>
Module not available. Install alpha_ladder_core to see the full honest assessment
of the coefficient derivation status.
</div>
""", unsafe_allow_html=True)

    c1, c2 = st.columns(2)
    c1.metric("Derivation achieved", "False")
    c2.metric("Best candidate", "N/A")

    col_works, col_fails = st.columns(2)

    with col_works:
        st.subheader("What Works")
        st.markdown("- Heat kernel 1/a_1 = 3 identifies c_2")
        st.markdown("- chi(T_{S^2}) = 3 confirms algebraic origin")
        st.markdown("- Vacuum polynomial produces c_3 = phi/2")

    with col_fails:
        st.subheader("What Fails")
        st.markdown("- No one-loop calculation derives c_2 = 3")
        st.markdown("- c_3 = phi/2 has no spectral origin (rationality no-go)")
        st.markdown("- No complete first-principles derivation exists")

# ---------------------------------------------------------------------------
# 10. One-Loop Result: Delta_tot = 43/15
# ---------------------------------------------------------------------------
st.header("10. One-Loop Result: 43/15")

st.markdown("""
<div class="formula-card">
<strong>Full one-loop spectral zeta calculation (zeta/phase3):</strong><br><br>
The complete quadratic fluctuation action on M4 x S2 with Dirac monopole
background gives 13 sectors (graviton, trace, dilaton, neutral vectors,
TT gravitons, charged vectors, FP ghosts, Lorenz ghosts). The total
one-loop correction to 1/G is:<br><br>
Delta_tot = 43/15 = 2.867
</div>
""", unsafe_allow_html=True)

st.markdown("")

c1, c2, c3 = st.columns(3)
c1.metric("Delta_tot (one-loop)", "43/15 = 2.867")
c2.metric("c_2 (topological)", "3")
c3.metric("Difference", "2/15 = 0.133")

st.markdown("")

st.markdown("""
<div class="warning-card">
<strong>The one-loop output (43/15) is NOT c_2.</strong><br><br>
c_2 = 3 is the topological/geometric identification (five-fold coincidence on S2).
The one-loop calculation gives 43/15 = 2.867, which differs by 2/15.
These are conceptually distinct: the loop identifies the mechanism (correct sign,
correct magnitude), but does not produce the exact integer.<br><br>
<strong>Effect on G:</strong> Using c_2 = 43/15 instead of 3 shifts G by
7.1 ppm. This is below the 22 ppm CODATA uncertainty but would be clearly
detectable by next-generation experiments aiming for 1-10 ppm.<br><br>
A sub-10 ppm measurement of G would distinguish between c_2 = 3 (topological)
and Delta_tot = 43/15 (one-loop).
</div>
""", unsafe_allow_html=True)

st.markdown("")

st.markdown("""
<div class="proof-card">
<strong>SUSY denominator no-go theorem (zeta/phase5):</strong><br><br>
No combination of 6D N=(1,0) SUSY multiplets can produce the residual
c_matter = 2/15. All SUSY matter contributions live in (1/12)*Z,
but 2/15 = 8/60 and 8 is not divisible by 5. The exact integer c_2 = 3
must have non-perturbative or topological origin.
</div>
""", unsafe_allow_html=True)

st.markdown("")

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Coefficient Derivation | Alpha Ladder Research Dashboard")
