"""
Page 43 -- Braneworld: Why SM Must Be on a Brane
Four escape routes for KK gauge matching vs alpha running.
Braneworld is the only viable option.
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
    from alpha_ladder_core.gauge_escape_routes import summarize_escape_routes
    _core_available = True
except ImportError:
    pass


@st.cache_data
def _get_summary():
    if _core_available:
        return summarize_escape_routes()
    return None


# ---------------------------------------------------------------------------
# Title and Introduction
# ---------------------------------------------------------------------------
st.title("Braneworld: Why SM Must Be on a Brane")
st.markdown("*Charged bulk matter on S^2 destroys alpha running -- the braneworld is mandatory*")

st.markdown("""
<div class="theorem-card">
<strong>The theorem:</strong> Charged bulk matter on S&sup2; shifts 1/alpha by ~10<sup>17</sup> ppm.
The only consistent picture is SM matter on a 4D brane, with gravity and gauge fields
propagating in the 6D bulk. This is the standard braneworld approach
(Randall-Sundrum, ADD, Horava-Witten) and is not ad hoc.
</div>
""", unsafe_allow_html=True)

if not _core_available:
    st.error(
        "Core module `alpha_ladder_core.gauge_escape_routes` not found. "
        "Install alpha_ladder_core to enable computations on this page."
    )

summary = _get_summary()

# ---------------------------------------------------------------------------
# Section 1: The Problem
# ---------------------------------------------------------------------------
st.header("1. The Problem")

st.markdown("""
<div class="warning-card">
<strong>Alpha running disaster:</strong><br><br>
Tree-level matching gives phi_vev = -0.597, so alpha_KK = alpha_EM = 1/137.036 exactly.
But if charged matter propagates in the 6D bulk on S&sup2; (a_0 = 28 um):<br>
-- KK tower has modes at ~18 meV spacing<br>
-- Billions of charged modes below the electron mass<br>
-- Shifts 1/alpha by ~10<sup>17</sup> ppm (observed precision: 0.37 ppb)<br>
-- DISASTROUS for alpha running consistency
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Section 2: Four Escape Routes
# ---------------------------------------------------------------------------
st.header("2. Four Escape Routes")

if summary:
    comp = summary.get("comparison", {})
    table = comp.get("table", [])

    if table:
        rows = []
        for r in table:
            rs = "Yes" if r.get("running_safe") is True else (
                "No" if r.get("running_safe") is False else str(r.get("running_safe", "N/A"))
            )
            ms = "Yes" if r.get("matching_works") is True else (
                "No" if r.get("matching_works") is False else str(r.get("matching_works", "N/A"))
            )
            eo = "Yes" if r.get("eotwash_open") else "No"
            rows.append({
                "Route": r.get("route", ""),
                "Running Safe": rs,
                "Matching Works": ms,
                "Eot-Wash Open": eo,
                "Viable": str(r.get("viable", "")),
                "Summary": r.get("summary", ""),
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
    else:
        st.info("No comparison table available.")
elif _core_available:
    st.info("No data returned from escape routes analysis.")
else:
    st.info("Core module not available.")

# ---------------------------------------------------------------------------
# Section 3: Route 2 -- Braneworld (the winner)
# ---------------------------------------------------------------------------
st.header("3. Route 2: Braneworld (the Winner)")

if summary:
    brane = summary.get("route_2_brane", {})
    brane_deep = summary.get("braneworld", {})

    st.markdown(f"""
<div class="proof-card">
<strong>Braneworld setup:</strong><br><br>
SM charged fields confined to a 4D brane at a point on S&sup2;.
Gravity + gauge fields propagate in the 6D bulk.<br><br>
<strong>Running safe:</strong> {brane.get('running_safe', 'N/A')}<br>
<strong>Matching works:</strong> {brane.get('matching_works', 'N/A')}<br>
<strong>Eot-Wash viable:</strong> {brane.get('eotwash_viable', 'N/A')}<br>
<strong>Verdict:</strong> {brane.get('verdict', 'N/A')}
</div>
""", unsafe_allow_html=True)

    # Bulk vs brane fields
    bulk_fields = brane.get("bulk_fields", [])
    brane_fields = brane.get("brane_fields", [])

    if bulk_fields or brane_fields:
        c1, c2 = st.columns(2)
        with c1:
            st.subheader("Bulk Fields")
            for f in bulk_fields:
                st.markdown(f"- {f}")
        with c2:
            st.subheader("Brane Fields")
            for f in brane_fields:
                st.markdown(f"- {f}")

    # Braneworld deep dive
    if brane_deep:
        c1, c2, c3 = st.columns(3)
        c1.metric("Position dependent", str(brane_deep.get("position_dependent", "N/A")))
        alpha_b = brane_deep.get("alpha_at_brane")
        c2.metric("alpha at brane", f"{alpha_b:.12e}" if alpha_b is not None else "N/A")
        c3.metric("Positions checked", str(brane_deep.get("brane_positions_checked", "N/A")))

        all_agree = brane_deep.get("all_positions_agree", False)
        st.metric("All positions agree", str(all_agree))

        squashed = brane_deep.get("squashed_note", "")
        if squashed:
            st.markdown(f"""
<div class="step-card">
{squashed}
</div>
""", unsafe_allow_html=True)

        # Bulk content details
        bulk_content = brane_deep.get("bulk_content", {})
        if bulk_content:
            st.subheader("Bulk Content Detail")
            rows = []
            for key, val in bulk_content.items():
                rows.append({
                    "Component": key.replace("_", " ").title(),
                    "Field": val.get("field", ""),
                    "4D Modes": val.get("4D_modes", ""),
                    "Role": val.get("role", ""),
                })
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
elif _core_available:
    st.info("No braneworld data returned.")
else:
    st.info("Core module not available.")

# ---------------------------------------------------------------------------
# Section 4: Routes That Fail
# ---------------------------------------------------------------------------
st.header("4. Routes That Fail")

if summary:
    # Route 1: Planck scale
    r1 = summary.get("route_1_planck", {})
    if r1:
        verdict = r1.get("verdict", "")
        st.markdown(f"""
<div class="step-card">
<strong>Route 1: Planck-scale radius</strong><br><br>
a_0 = l_Pl: all KK modes at Planck scale. Running is safe but no observables.<br>
M_6 = {r1.get('M_6', 0):.4e} eV<br>
Running safe: {r1.get('running_safe', 'N/A')}<br>
Eot-Wash viable: {r1.get('eotwash_viable', 'N/A')}<br>
<strong>Verdict:</strong> {verdict}
</div>
""", unsafe_allow_html=True)

    # Route 3: Tiny charge
    r3 = summary.get("route_3_tiny_q", {})
    if r3:
        verdict = r3.get("verdict", "")
        q_max = r3.get("Q_max_for_safe_running")
        alpha_at_q = r3.get("alpha_KK_at_Q_max")
        st.markdown(f"""
<div class="step-card">
<strong>Route 3: Tiny charge Q</strong><br><br>
Reduce the KK charge Q to suppress running contribution.<br>
Q_max for safe running: {f'{q_max:.4e}' if q_max is not None else 'N/A'}<br>
alpha_KK at Q_max: {f'{alpha_at_q:.4e}' if alpha_at_q is not None else 'N/A'}<br>
Shift at Q=1: {f'{r3.get("shift_ppm_Q1", 0):.4e}'} ppm<br>
<strong>Verdict:</strong> {verdict}
</div>
""", unsafe_allow_html=True)

    # Route 4: Heavy bulk matter
    r4 = summary.get("route_4_heavy_bulk", {})
    if r4:
        verdict = r4.get("verdict", "")
        st.markdown(f"""
<div class="step-card">
<strong>Route 4: Heavy bulk matter</strong><br><br>
Give the 6D scalar a bulk mass M_bulk to decouple KK modes.<br>
Need M_bulk > M_6 to fully decouple.<br>
KK spacing: {f'{r4.get("kk_spacing_eV", 0):.4e}'} eV<br>
<strong>Verdict:</strong> {verdict}
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No route data returned.")
else:
    st.info("Core module not available.")

# ---------------------------------------------------------------------------
# Section 5: Honest Assessment
# ---------------------------------------------------------------------------
st.header("5. Honest Assessment")

if summary:
    key_messages = summary.get("key_messages", [])
    if key_messages:
        for i, msg in enumerate(key_messages, 1):
            st.markdown(f"""
<div class="step-card">
<strong>{i}.</strong> {msg}
</div>
""", unsafe_allow_html=True)

    st.markdown("""
<div class="warning-card">
<strong>Bottom line:</strong> The braneworld interpretation is the ONLY option that
preserves all three requirements: alpha matching, safe running, and the Eot-Wash window.
This is the standard approach in extra-dimension physics and is not ad hoc.
Paper 3 should adopt the braneworld interpretation.
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
st.caption("Braneworld | Alpha Ladder Research Dashboard")
