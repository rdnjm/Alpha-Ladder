"""
Page 44 -- KK Spectrum on S^2
Full Kaluza-Klein spectrum with dilaton vev matched to alpha_EM.
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
    from alpha_ladder_core.gauge_kk_spectrum_matched import (
        summarize_kk_spectrum,
        spectrum_at_notable_a0,
        matched_vs_unmatched,
    )
    _core_available = True
except ImportError:
    pass


@st.cache_data
def _get_summary():
    if _core_available:
        return summarize_kk_spectrum()
    return None


@st.cache_data
def _get_notable_a0():
    if _core_available:
        return spectrum_at_notable_a0(matched=True)
    return None


@st.cache_data
def _get_matched_vs_unmatched():
    if _core_available:
        return matched_vs_unmatched(2.8e-5)  # 28 um (Eot-Wash)
    return None


# ---------------------------------------------------------------------------
# Title and Introduction
# ---------------------------------------------------------------------------
st.title("KK Spectrum on S\u00b2")
st.markdown("*Full Kaluza-Klein spectrum with dilaton vev matched to alpha_EM*")

st.markdown("""
<div class="formula-card">
<strong>KK mass formulas on S&sup2;:</strong><br><br>
Scalar: m_l&sup2; = l(l+1)/R&sup2;, l = 0, 1, 2, ... deg = 2l+1<br>
Vector: m_l&sup2; = l(l+1)/R&sup2;, l = 1, 2, 3, ... deg = 2(2l+1)<br>
Graviton: m_l&sup2; = [l(l+1)-2]/R&sup2;, l = 2, 3, 4, ... deg = 2l+1<br>
Fermion: m_l&sup2; = (l+1/2)&sup2;/R&sup2;, l = 0, 1, 2, ... deg = 2(2l+1)<br><br>
Physical radius: R_phys = a_0 * e^(phi_vev), with phi_vev = (1/4) ln(4 pi alpha_EM).
</div>
""", unsafe_allow_html=True)

if not _core_available:
    st.error(
        "Core module `alpha_ladder_core.gauge_kk_spectrum_matched` not found. "
        "Install alpha_ladder_core to enable computations on this page."
    )

summary = _get_summary()

# ---------------------------------------------------------------------------
# Section 1: Matched vs Unmatched
# ---------------------------------------------------------------------------
st.header("1. Matched vs Unmatched")

if summary:
    consts = summary.get("constants", {})

    c1, c2, c3 = st.columns(3)
    phi_vev = consts.get("phi_vev")
    c1.metric("phi_vev", f"{phi_vev:.6f}" if phi_vev is not None else "N/A")
    e_phi = consts.get("e_phi_vev")
    c2.metric("e^(phi_vev)", f"{e_phi:.6f}" if e_phi is not None else "N/A")
    enhancement = consts.get("mass_enhancement")
    c3.metric("Mass enhancement", f"{enhancement:.6f}x" if enhancement is not None else "N/A")

    e_phi_str = f"{e_phi:.4f}" if e_phi else "N/A"
    enh_str = f"{enhancement:.3f}" if enhancement else "N/A"
    st.markdown(f"""
<div class="step-card">
<strong>R_phys = a_0 * e^(phi_vev) = a_0 * {e_phi_str}</strong><br>
All KK modes are {enh_str}x heavier than naive (unmatched) case.
This is physically required for consistency with observed alpha_EM.
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No constants data returned.")
else:
    st.info("Core module not available.")

# ---------------------------------------------------------------------------
# Section 2: Key Results
# ---------------------------------------------------------------------------
st.header("2. Key Results")

if summary:
    key_results = summary.get("key_results", [])
    for i, result in enumerate(key_results, 1):
        st.markdown(f"""
<div class="step-card">
<strong>{i}.</strong> {result}
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No key results returned.")
else:
    st.info("Core module not available.")

# ---------------------------------------------------------------------------
# Section 3: Spectrum at Notable Radii
# ---------------------------------------------------------------------------
st.header("3. Spectrum at Notable Radii")

notable = _get_notable_a0()

if notable:
    rows = []
    for entry in notable:
        lightest = entry.get("lightest_5", [])
        lightest_str = ""
        if lightest:
            first = lightest[0]
            lightest_str = f"{first.get('spin', '')} l={first.get('l', '')} m={first.get('m_eV', 0):.4e} eV"
        rows.append({
            "a_0": entry.get("a_0_label", ""),
            "M_6 (TeV)": f"{entry.get('M_6_TeV', 0):.4e}",
            "R_phys (m)": f"{entry.get('R_phys_m', 0):.4e}",
            "Lightest mode": lightest_str,
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
elif _core_available:
    st.info("No notable a_0 data returned.")
else:
    st.info("Core module not available.")

# ---------------------------------------------------------------------------
# Section 4: Matched vs Unmatched (a_0 = 28 um)
# ---------------------------------------------------------------------------
st.header("4. Matched vs Unmatched (a_0 = 28 um)")

mvu = _get_matched_vs_unmatched()

if mvu:
    c1, c2, c3 = st.columns(3)
    mass_ratio = mvu.get("mass_ratio")
    c1.metric("Mass ratio (matched/unmatched)", f"{mass_ratio:.6f}" if mass_ratio is not None else "N/A")
    r_m = mvu.get("R_matched_m")
    c2.metric("R matched (m)", f"{r_m:.4e}" if r_m is not None else "N/A")
    r_u = mvu.get("R_unmatched_m")
    c3.metric("R unmatched (m)", f"{r_u:.4e}" if r_u is not None else "N/A")

    matched_modes = mvu.get("matched_modes", [])
    unmatched_modes = mvu.get("unmatched_modes", [])
    if matched_modes and unmatched_modes:
        rows = []
        for mm, um in zip(matched_modes, unmatched_modes):
            ratio = mm["m_eV"] / um["m_eV"] if um.get("m_eV", 0) > 0 else 0.0
            rows.append({
                "l": mm.get("l", ""),
                "m_matched (eV)": f"{mm.get('m_eV', 0):.6e}",
                "m_unmatched (eV)": f"{um.get('m_eV', 0):.6e}",
                "ratio": f"{ratio:.4f}",
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    note = mvu.get("detectability_note", "")
    if note:
        st.markdown(f"""
<div class="step-card">
{note}
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No matched-vs-unmatched data returned.")
else:
    st.info("Core module not available.")

# ---------------------------------------------------------------------------
# Section 5: S^2 vs ADD
# ---------------------------------------------------------------------------
st.header("5. S\u00b2 vs ADD")

st.markdown("""
<div class="step-card">
<strong>Spectrum comparison:</strong><br><br>
Unlike ADD's quasi-continuous KK lattice on T^n, the S&sup2; spectrum is discrete
with degeneracy (2l+1). This gives a distinctive pattern in cross-sections:
resonance peaks at specific energies, not a continuum.<br><br>
The S&sup2; spectrum is sparser than ADD -- approximately 26% of the ADD mode count
at the same energy cutoff. The factor of 1/pi arises because T&sup2; modes fill a
2D disk (N ~ pi r&sup2;) while S&sup2; modes fill a 1D tower (N ~ l_max&sup2;, no pi factor).
</div>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Section 4: Honest Assessment
# ---------------------------------------------------------------------------
st.header("6. Honest Assessment")

if summary:
    st.markdown("""
<div class="warning-card">
<strong>Honest Assessment:</strong><br><br>
The matched dilaton vev is physically required but makes detection harder.
All KK modes shift up by 1.817x.
The S&sup2; spectrum is qualitatively different from ADD (discrete vs lattice),
leading to distinct experimental signatures.
The spectrum is fully determined by a single parameter (a_0)
once alpha_EM is matched.
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
st.caption("KK Spectrum on S\u00b2 | Alpha Ladder Research Dashboard")
