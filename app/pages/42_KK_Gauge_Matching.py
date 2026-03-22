"""
Page 42 -- KK Gauge Matching
Tree-level matching of KK gauge fields to electromagnetism
and the Coleman-Weinberg effective potential for the dilaton.
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
    from alpha_ladder_core.gauge_kk_gauge_matching import summarize_gauge_matching
    _core_available = True
except ImportError:
    pass


@st.cache_data
def _get_summary():
    if _core_available:
        return summarize_gauge_matching()
    return None


# ---------------------------------------------------------------------------
# Title and Introduction
# ---------------------------------------------------------------------------
st.title("KK Gauge Matching")
st.markdown("*Tree-level matching of KK gauge coupling to alpha_EM and Coleman-Weinberg potential*")

st.markdown("""
<div class="formula-card">
<strong>Gauge matching condition:</strong><br><br>
phi_vev = (1/4) ln(4 pi alpha_EM) &approx; -0.597<br>
alpha_KK = e<sup>4 phi_vev</sup> / (4 pi) = alpha_EM exactly<br><br>
The S&sup2; reduction yields two gauge fields from the SO(3)/SO(2) coset.
Identifying one with the photon fixes the dilaton vev through alpha_KK = alpha_EM.
</div>
""", unsafe_allow_html=True)

if not _core_available:
    st.error(
        "Core module `alpha_ladder_core.gauge_kk_gauge_matching` not found. "
        "Install alpha_ladder_core to enable computations on this page."
    )

summary = _get_summary()

# ---------------------------------------------------------------------------
# Section 1: Tree-Level Matching
# ---------------------------------------------------------------------------
st.header("1. Tree-Level Matching")

if summary:
    tree = summary.get("tree_level", {})
    c1, c2, c3 = st.columns(3)
    phi_vev = tree.get("phi_vev")
    c1.metric("phi_vev", f"{phi_vev:.6f}" if phi_vev is not None else "N/A")
    alpha_kk = tree.get("alpha_KK")
    c2.metric("alpha_KK", f"{alpha_kk:.12e}" if alpha_kk is not None else "N/A")
    verification = tree.get("verification")
    c3.metric("Verification", str(verification) if verification is not None else "N/A")

    c1b, c2b, c3b = st.columns(3)
    g_kk = tree.get("g_KK")
    c1b.metric("g_KK", f"{g_kk:.8e}" if g_kk is not None else "N/A")
    e_phi = tree.get("e_to_phi")
    c2b.metric("e^(phi_vev)", f"{e_phi:.8f}" if e_phi is not None else "N/A")
    e_4phi = tree.get("e_to_4phi")
    c3b.metric("e^(4 phi_vev)", f"{e_4phi:.8f}" if e_4phi is not None else "N/A")
elif _core_available:
    st.info("No data returned from tree-level matching.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 2: KK Mass Spectrum
# ---------------------------------------------------------------------------
st.header("2. KK Mass Spectrum")

if summary:
    spec = summary.get("kk_spectrum", {})
    modes = spec.get("modes", [])

    if modes:
        r_phys = spec.get("R_phys")
        if r_phys is not None:
            st.metric("R_phys (Planck lengths)", f"{r_phys:.6f}")

        rows = []
        for m in modes:
            rows.append({
                "l": m.get("l", ""),
                "m_l": f"{m.get('m_l', 0):.6e}",
                "m_l / m_1": f"{m.get('m_l_over_m1', 0):.6f}",
                "Degeneracy": m.get("degeneracy", ""),
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
    else:
        st.info("No spectrum data available.")
elif _core_available:
    st.info("No spectrum data returned.")
else:
    st.info("Core module not available.")

# ---------------------------------------------------------------------------
# Section 3: Coleman-Weinberg Potential
# ---------------------------------------------------------------------------
st.header("3. Coleman-Weinberg Potential")

if summary:
    cw_scan = summary.get("cw_scan", {})
    cw_min = summary.get("cw_minimum", {})

    st.markdown("""
<div class="step-card">
<strong>Coleman-Weinberg one-loop effective potential:</strong><br>
V_CW = Sum m^4(phi) log(m^2(phi) / mu^2) / (64 pi^2)<br><br>
The CW potential generates a phi-dependent quantum correction.
We scan for coincidence between phi_CW_min and phi_alpha_match.
</div>
""", unsafe_allow_html=True)

    c1, c2, c3 = st.columns(3)
    has_min = cw_scan.get("has_minimum", False)
    c1.metric("CW minimum found", str(has_min))
    phi_min = cw_scan.get("phi_minimum")
    c2.metric("phi at minimum", f"{phi_min:.6f}" if phi_min is not None else "None")
    v_min = cw_scan.get("V_at_minimum")
    c3.metric("V at minimum", f"{v_min:.6e}" if v_min is not None else "None")

    # Refined CW minimum
    if cw_min:
        desc = cw_min.get("description", "")
        if desc:
            st.markdown(f"""
<div class="formula-card">
<strong>Refined CW minimum (bisection):</strong><br>
{desc}
</div>
""", unsafe_allow_html=True)

        honest = cw_min.get("honest_assessment", "")
        if honest:
            st.markdown(f"""
<div class="step-card">
{honest}
</div>
""", unsafe_allow_html=True)

    # n_charged scan
    n_scan = summary.get("n_charged_scan", {})
    if n_scan:
        results = n_scan.get("results", [])
        if results:
            # Check if any entry has a minimum
            any_minimum = any(r.get("phi_min") is not None for r in results)
            st.subheader("n_charged Scan")
            if any_minimum:
                rows = []
                for r in results:
                    phi_str = f"{r['phi_min']:.6f}" if r.get("phi_min") is not None else "No min"
                    ratio_str = f"{r['ratio_to_vev']:.4f}" if r.get("ratio_to_vev") is not None else "--"
                    dev_str = f"{r['deviation_percent']:.2f}%" if r.get("deviation_percent") is not None else "--"
                    rows.append({
                        "n_charged": r.get("n_charged", ""),
                        "phi_min": phi_str,
                        "ratio_to_vev": ratio_str,
                        "deviation": dev_str,
                    })
                st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
            else:
                st.markdown("""
<div class="warning-card">
<strong>Result:</strong> The Coleman-Weinberg potential has NO minimum for any
n_charged value tested (1 to 100). The one-loop CW potential does not
dynamically determine phi_vev. The tree-level gauge matching
phi_vev = (1/4)*ln(4*pi*alpha) remains the only mechanism that fixes
the dilaton vev.  This is an honest negative result.
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No CW data returned.")
else:
    st.info("Core module not available.")

# ---------------------------------------------------------------------------
# Section 4: Alpha Running and G Consistency
# ---------------------------------------------------------------------------
st.header("4. Alpha Running and G Consistency")

if summary:
    ar = summary.get("alpha_running", {})
    if ar:
        c1, c2 = st.columns(2)
        c1.metric("alpha at lab", f"{ar.get('alpha_at_lab', 0):.8e}")
        c2.metric("alpha at compactification", f"{ar.get('alpha_at_compactification', 0):.8e}")

        desc = ar.get("description", "")
        if desc:
            st.markdown(f"""
<div class="step-card">
{desc}
</div>
""", unsafe_allow_html=True)

    g_check = summary.get("G_consistency", {})
    if g_check:
        any_match = g_check.get("any_match", False)
        st.metric("Any G / golden ratio match", str(any_match))

        honest = g_check.get("honest_assessment", "")
        if honest:
            st.markdown(f"""
<div class="step-card">
{honest}
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No running / consistency data returned.")
else:
    st.info("Core module not available.")

# ---------------------------------------------------------------------------
# Section 5: Honest Assessment
# ---------------------------------------------------------------------------
st.header("5. Honest Assessment")

if summary:
    s = summary.get("summary", {})
    for key in ["key_result_1", "key_result_2", "key_result_3"]:
        val = s.get(key, "")
        if val:
            st.markdown(f"""
<div class="step-card">
<strong>{key.replace('_', ' ').title()}:</strong> {val}
</div>
""", unsafe_allow_html=True)

    for key in ["gap_1_status", "gap_2_status"]:
        val = s.get(key, "")
        if val:
            st.markdown(f"""
<div class="theorem-card">
{val}
</div>
""", unsafe_allow_html=True)

    caveat = s.get("honest_caveat", "")
    if caveat:
        st.markdown(f"""
<div class="warning-card">
<strong>Honest Caveat:</strong><br><br>
{caveat}
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No summary data returned.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("KK Gauge Matching | Alpha Ladder Research Dashboard")
