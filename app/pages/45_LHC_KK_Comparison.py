"""
Page 45 -- LHC KK Comparison: S^2 vs ADD
Cross-section ratio and rescaled exclusion bounds.
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
    from alpha_ladder_core.gauge_lhc_kk_comparison import summarize_lhc_comparison
    _core_available = True
except ImportError:
    pass


@st.cache_data
def _get_summary():
    if _core_available:
        return summarize_lhc_comparison()
    return None


# ---------------------------------------------------------------------------
# Title and Introduction
# ---------------------------------------------------------------------------
st.title("LHC KK Comparison: S\u00b2 vs ADD")
st.markdown("*Cross-section ratio and rescaled exclusion bounds for the Alpha Ladder S^2 compactification*")

if not _core_available:
    st.error(
        "Core module `alpha_ladder_core.gauge_lhc_kk_comparison` not found. "
        "Install alpha_ladder_core to enable computations on this page."
    )

summary = _get_summary()

# ---------------------------------------------------------------------------
# Section 1: Mode Counting
# ---------------------------------------------------------------------------
st.header("1. Mode Counting")

st.markdown("""
<div class="formula-card">
<strong>Lightest KK modes:</strong><br><br>
ADD (T&sup2;): m_1 = 1/R<br>
S&sup2;: m_1 = sqrt(2)/R_phys (scalar/vector l=1)<br><br>
<strong>Mode count scaling (continuum limit):</strong><br>
N_ADD = pi * (E*R_ADD)&sup2; [Gauss circle problem]<br>
N_S2 = (E*R_S2)&sup2; [sum of 2l+1 up to l_max ~ E*R]<br>
N_S2/N_ADD ~ (1/pi) * (R_S2/R_ADD)&sup2; ~ 0.26
</div>
""", unsafe_allow_html=True)

if summary:
    r_ratio = summary.get("R_ratio_analytic")

    c1, c2, c3 = st.columns(3)
    c1.metric("Lightest S\u00b2 mode", "sqrt(2)/R ~ 1.414/R")
    c2.metric("S\u00b2/ADD mode ratio", f"{r_ratio**2 / 3.14159:.4f}" if r_ratio is not None else "~0.26")
    c3.metric("R_S2/R_ADD", f"{r_ratio:.6f}" if r_ratio is not None else "N/A")
elif _core_available:
    st.info("No mode counting data returned.")
else:
    st.info("Core module not available.")

# ---------------------------------------------------------------------------
# Section 2: LHC Exclusion Rescaling
# ---------------------------------------------------------------------------
st.header("2. LHC Exclusion Rescaling")

if summary:
    excl = summary.get("exclusion", {})

    c1, c2, c3 = st.columns(3)
    m6_add = excl.get("M_6_excl_add")
    c1.metric("ADD M_6 exclusion", f"{m6_add:.1f} TeV" if m6_add is not None else "5.0 TeV")
    m6_s2 = excl.get("M_6_excl_s2")
    c2.metric("S\u00b2 rescaled M_6", f"{m6_s2:.3f} TeV" if m6_s2 is not None else "N/A")
    ratio = excl.get("ratio_at_exclusion")
    c3.metric("N_eff ratio at exclusion", f"{ratio:.6f}" if ratio is not None else "N/A")

    c1b, c2b, c3b = st.columns(3)
    m6_s2_sig = excl.get("M_6_excl_s2_sigma")
    c1b.metric("S\u00b2 exclusion (sigma)", f"{m6_s2_sig:.3f} TeV" if m6_s2_sig is not None else "N/A")
    sig_ratio = excl.get("sigma_ratio_at_exclusion")
    c2b.metric("sigma ratio", f"{sig_ratio:.6f}" if sig_ratio is not None else "N/A")
    a0_thresh = excl.get("a_0_threshold_s2_m")
    c3b.metric("a_0 threshold (S\u00b2)", f"{a0_thresh:.3e} m" if a0_thresh is not None else "N/A")

    # Cross-section scan table
    scan = summary.get("scan", [])
    if scan:
        st.subheader("Cross-Section Ratio Scan (E_cm = 14 TeV)")
        rows = []
        for row in scan:
            rows.append({
                "M_6 (TeV)": f"{row.get('M_6', 0):.1f}",
                "N_ADD": f"{row.get('N_add', 0):.6e}",
                "N_S2": f"{row.get('N_s2', 0):.6e}",
                "N ratio": f"{row.get('N_ratio', 0):.6f}",
                "sigma ratio": f"{row.get('sigma_ratio', 0):.6f}",
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    # Spectrum comparison
    spec_table = summary.get("spectrum_table", {})
    add_levels = spec_table.get("add_levels", [])
    s2_levels = spec_table.get("s2_levels", [])
    if add_levels and s2_levels:
        st.subheader("Spectrum Comparison (M_6 = 5 TeV, first 20 levels)")
        n_show = min(20, len(add_levels), len(s2_levels))
        add_cum = spec_table.get("add_cumulative", [])
        s2_cum = spec_table.get("s2_cumulative", [])
        rows = []
        for i in range(n_show):
            m_a, d_a = add_levels[i]
            l_s, m_s, d_s = s2_levels[i]
            rows.append({
                "ADD level": i + 1,
                "m_ADD (eV)": f"{m_a:.4e}",
                "deg_ADD": d_a,
                "cum_ADD": add_cum[i] if i < len(add_cum) else "",
                "l (S\u00b2)": l_s,
                "m_S2 (eV)": f"{m_s:.4e}",
                "deg_S2": d_s,
                "cum_S2": s2_cum[i] if i < len(s2_cum) else "",
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
elif _core_available:
    st.info("No exclusion data returned.")
else:
    st.info("Core module not available.")

# ---------------------------------------------------------------------------
# Section 3: Eot-Wash Window
# ---------------------------------------------------------------------------
st.header("3. Eot-Wash Window")

if summary:
    excl = summary.get("exclusion", {})
    m6_s2 = excl.get("M_6_excl_s2", 0)

    st.markdown(f"""
<div class="step-card">
<strong>Eot-Wash window status:</strong><br><br>
The Eot-Wash optimal point is a_0 = 28 um, corresponding to M_6 ~ 5 TeV.<br>
The S&sup2; LHC exclusion is M_6 > {m6_s2:.3f} TeV (weaker than ADD's 5.0 TeV).<br><br>
Since the S&sup2; bound is looser than ADD, the Eot-Wash window is marginally
at the boundary. The S&sup2; spectrum produces fewer modes than ADD at the same
M_6, so the LHC exclusion is relaxed.
</div>
""", unsafe_allow_html=True)

    st.markdown("""
<div class="theorem-card">
<strong>Key insight:</strong> The 5 TeV ADD exclusion does NOT apply directly to Alpha Ladder.
The correct exclusion uses the S&sup2; mode count, which is ~26% of ADD.
This loosens the bound significantly.
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No Eot-Wash data returned.")
else:
    st.info("Core module not available.")

# ---------------------------------------------------------------------------
# Section 4: Honest Assessment
# ---------------------------------------------------------------------------
st.header("4. Honest Assessment")

if summary:
    st.markdown("""
<div class="warning-card">
<strong>Honest Assessment:</strong><br><br>
The S&sup2; KK spectrum is sparser than ADD for light modes
(lightest mode sqrt(2)/R vs 1/R, ratio sqrt(2) ~ 1.414).
Total mode count scales as R&sup2; E&sup2; for both, but S&sup2; has a smaller
prefactor: N_S2/N_ADD ~ (R_S2/R_ADD)&sup2;/pi ~ 0.26.
The LHC exclusion for S&sup2; is WEAKER than ADD, meaning
the Eot-Wash window remains open (marginally).
The phase-space weighted ratio equals the mode count ratio
in the continuum limit (both give the same 1/2 suppression for n=2).
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
st.caption("LHC KK Comparison | Alpha Ladder Research Dashboard")
