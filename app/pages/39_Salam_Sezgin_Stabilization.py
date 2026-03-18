"""
Page 39 -- Salam-Sezgin Radius Stabilization
Lambda_6 > 0 cosmological constant fixes the internal radius a_0.
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
    from alpha_ladder_core.salam_sezgin_stabilization import summarize_salam_sezgin
    _core_available = True
except ImportError:
    pass

_charts_available = False
try:
    from app.components.charts import (
        ss_potential_chart,
        ss_lambda6_scan_chart,
    )
    _charts_available = True
except ImportError:
    pass


@st.cache_data
def _get_ss_summary(_constants):
    if _core_available:
        return summarize_salam_sezgin(_constants)
    return None


# ---------------------------------------------------------------------------
# Title and Introduction
# ---------------------------------------------------------------------------
st.title("Salam-Sezgin Radius Stabilization")
st.markdown("*Fixing the internal radius via 6D cosmological constant*")

st.markdown("""
<div class="warning-card">
<strong>The Problem:</strong> The internal radius a_0 was unfixed by all three mechanisms
in the previous analysis (Coleman-Weinberg, warped S^2, orbifold S^2/Z_2).
The Salam-Sezgin mechanism adds a positive 6D cosmological constant Lambda_6 > 0
to the action, which competes with curvature and flux to produce a stable minimum
for the radius modulus.
</div>
""", unsafe_allow_html=True)

if not _core_available:
    st.error(
        "Core module `alpha_ladder_core.salam_sezgin_stabilization` not found. "
        "Install alpha_ladder_core to enable computations on this page."
    )

summary = _get_ss_summary(constants)

# ---------------------------------------------------------------------------
# Section 1: Gauge Matching
# ---------------------------------------------------------------------------
st.header("1. Gauge Matching")

st.markdown("""
<div class="formula-card">
<strong>Salam-Sezgin effective potential:</strong><br><br>
V(R) = Lambda_6 * R^2 - 2 / R^2 + n^2 / (2 * R^4)<br><br>
The positive Lambda_6 term grows at large R while curvature and flux dominate at small R,
producing a stable minimum at finite R_0.
</div>
""", unsafe_allow_html=True)

if summary:
    gauge = summary.get("gauge_match", {})
    if gauge:
        c1, c2, c3 = st.columns(3)
        r_phys = gauge.get("R_phys")
        c1.metric("R_phys", f"{r_phys:.4e}" if r_phys is not None else "N/A")
        lambda_6 = gauge.get("Lambda_6")
        c2.metric("Lambda_6", f"{lambda_6:.4e}" if lambda_6 is not None else "N/A")
        c3.metric("Is natural", str(gauge.get("is_natural", "N/A")))

        c1, c2, c3 = st.columns(3)
        phi_vev = gauge.get("phi_vev")
        c1.metric("phi_vev", f"{phi_vev:.4e}" if phi_vev is not None else "N/A")
        g_kk = gauge.get("g_KK")
        c2.metric("g_KK", f"{g_kk:.4e}" if g_kk is not None else "N/A")
        status = gauge.get("status", "N/A")
        c3.metric("Status", str(status))
    else:
        st.info("Gauge matching data not available.")
elif _core_available:
    st.info("No data returned from Salam-Sezgin module.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 2: Effective Potential V(R)
# ---------------------------------------------------------------------------
st.header("2. Effective Potential V(R)")

if summary:
    potential = summary.get("potential")
    if potential:
        if _charts_available:
            st.plotly_chart(ss_potential_chart(potential), use_container_width=True)
        else:
            st.info("Charts module not available.")

        gauge = summary.get("gauge_match", {})
        dilaton = summary.get("dilaton_mass", {})

        c1, c2, c3 = st.columns(3)
        r_0 = gauge.get("R_phys") if gauge else None
        c1.metric("R_0 (minimum)", f"{r_0:.4e}" if r_0 is not None else "N/A")

        # V at minimum from potential data
        v_total = potential.get("V_total")
        if v_total is not None and hasattr(v_total, '__len__'):
            v_min = min(v_total)
            c2.metric("V at minimum", f"{v_min:.4e}")
        else:
            c2.metric("V at minimum", "N/A")

        v_double_prime = dilaton.get("m_phi_squared") if dilaton else None
        c3.metric("V''(R_0)", f"{v_double_prime:.4e}" if v_double_prime is not None else "N/A")
    else:
        st.info("Potential data not available.")
elif _core_available:
    st.info("No data returned from Salam-Sezgin module.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 3: Dilaton Mass
# ---------------------------------------------------------------------------
st.header("3. Dilaton Mass")

if summary:
    dilaton = summary.get("dilaton_mass")
    if dilaton:
        c1, c2, c3 = st.columns(3)
        m_phi_eV = dilaton.get("m_phi_eV")
        c1.metric("m_phi (eV)", f"{m_phi_eV:.4e}" if m_phi_eV is not None else "N/A")
        lambda_c = dilaton.get("lambda_compton_m")
        c2.metric("Compton wavelength (m)", f"{lambda_c:.4e}" if lambda_c is not None else "N/A")
        mass_scale = dilaton.get("mass_scale", "N/A")
        c3.metric("Mass scale", str(mass_scale))

        st.markdown(f"""
<div class="formula-card">
<strong>Dilaton mass from V''(R_0):</strong><br><br>
m_phi^2 = V''(R_0) evaluated at the Salam-Sezgin minimum.<br>
m_phi = {m_phi_eV:.4e} eV, Compton wavelength = {lambda_c:.4e} m<br>
Mass scale: {mass_scale}
</div>
""", unsafe_allow_html=True)
    else:
        c1, c2, c3 = st.columns(3)
        c1.metric("m_phi (eV)", "N/A")
        c2.metric("Compton wavelength (m)", "N/A")
        c3.metric("Mass scale", "N/A")
        st.markdown("""
<div class="step-card">
<strong>Dilaton mass:</strong> Data not available. The dilaton mass is determined
by the curvature of the potential at the minimum.
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from Salam-Sezgin module.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 4: Lambda_6 Scan
# ---------------------------------------------------------------------------
st.header("4. Lambda_6 Scan")

if summary:
    scan = summary.get("scan")
    if scan:
        if _charts_available:
            st.plotly_chart(ss_lambda6_scan_chart(scan), use_container_width=True)

        scan_results = scan.get("results", [])
        if scan_results:
            rows = []
            for entry in scan_results:
                rows.append({
                    "Lambda_6": f"{entry.get('Lambda_6', 0):.4e}",
                    "R_0": f"{entry.get('R_0', 0):.4e}",
                    "m_phi (eV)": f"{entry.get('m_phi_eV', 0):.4e}",
                    "Discrepancy (orders)": f"{entry.get('discrepancy_orders', 0):.1f}",
                })
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

        st.markdown(f"""
<div class="step-card">
<strong>Scan summary:</strong> {scan.get('n_scanned', 'N/A')} values of Lambda_6 tested.
All stable: {scan.get('all_stable', 'N/A')}.
R_0 range: {scan.get('R_0_range', 'N/A')}.
Mass range: {scan.get('mass_range_eV', 'N/A')} eV.
{scan.get('description', '')}
</div>
""", unsafe_allow_html=True)
    else:
        st.info("Scan data not available.")
elif _core_available:
    st.info("No data returned from Salam-Sezgin module.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 5: Cosmological Constant
# ---------------------------------------------------------------------------
st.header("5. Cosmological Constant")

if summary:
    cc = summary.get("cc_analysis")
    if cc:
        c1, c2 = st.columns(2)
        lambda_4 = cc.get("Lambda_4")
        c1.metric("Lambda_4", f"{lambda_4:.4e}" if lambda_4 is not None else "N/A")
        disc = cc.get("discrepancy_orders")
        c2.metric("Discrepancy (orders)", f"{disc:.1f}" if disc is not None else "N/A")

        c1, c2 = st.columns(2)
        lambda_4_obs = cc.get("Lambda_4_obs")
        c1.metric("Lambda_4 (observed)", f"{lambda_4_obs:.4e}" if lambda_4_obs is not None else "N/A")
        sign = cc.get("sign", "N/A")
        c2.metric("Sign", str(sign))

        honest = cc.get("honest_assessment", "")
        if honest:
            st.markdown(f"""
<div class="warning-card">
<strong>CC Honest Assessment:</strong><br><br>
{honest}
</div>
""", unsafe_allow_html=True)
    else:
        c1, c2 = st.columns(2)
        c1.metric("Lambda_4", "N/A")
        c2.metric("Discrepancy (orders)", "N/A")
        st.markdown("""
<div class="step-card">
<strong>Cosmological constant analysis:</strong> Data not available.
The CC discrepancy is expected to be ~122 orders of magnitude.
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from Salam-Sezgin module.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 6: Gap Closure
# ---------------------------------------------------------------------------
st.header("6. Gap Closure")

if summary:
    gap1 = summary.get("gap1_status", "")
    radius_fixed = summary.get("radius_fixed", False)
    overall = summary.get("overall_assessment", "")

    if gap1:
        st.markdown(f"""
<div class="proof-card">
<strong>Gap 1 Status:</strong><br><br>
{gap1}<br><br>
Radius fixed: {'Yes' if radius_fixed else 'No'}<br>
First principles: {'Yes' if summary.get('first_principles', False) else 'No'}
</div>
""", unsafe_allow_html=True)
    else:
        st.markdown("""
<div class="step-card">
<strong>Gap 1:</strong> Status data not available.
</div>
""", unsafe_allow_html=True)

    if overall:
        st.markdown(f"""
<div class="theorem-card">
<strong>Overall Assessment:</strong><br><br>
{overall}
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from Salam-Sezgin module.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 7: Honest Assessment
# ---------------------------------------------------------------------------
st.header("7. Honest Assessment")

if summary:
    honest = summary.get("honest_assessment", "")
    if honest:
        st.markdown(f"""
<div class="warning-card">
<strong>Honest Assessment:</strong><br><br>
{honest}
</div>
""", unsafe_allow_html=True)
    else:
        st.markdown("""
<div class="step-card">
<strong>Assessment:</strong> No honest assessment available from the module.
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from Salam-Sezgin module.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Salam-Sezgin Stabilization | Alpha Ladder Research Dashboard")
