"""
Alpha Ladder -- Home Content Page

Overview of the Alpha Ladder theory and its current results.
"""

import streamlit as st
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from app.components.sidebar import render_sidebar  # noqa: E402

constants = render_sidebar()

# ---------------------------------------------------------------------------
# Title
# ---------------------------------------------------------------------------
st.title("Alpha Ladder Research Dashboard")
st.markdown("**Predicting Newton's Constant from the Fine-Structure Constant**")
st.divider()

# ---------------------------------------------------------------------------
# Key Results -- metric cards
# ---------------------------------------------------------------------------
col1, col2, col3 = st.columns(3)

with col1:
    st.metric(
        label="G residual",
        value="-0.33 ppm",
        delta="CODATA 2022",
        delta_color="normal",
    )

with col2:
    st.metric(
        label="One-loop corroboration",
        value="43/15 = 2.867",
        delta="7.1 ppm from c_2=3",
        delta_color="normal",
    )

with col3:
    st.metric(
        label="Fitted parameters",
        value="0",
        delta="fully derived",
        delta_color="normal",
    )

st.divider()

# ---------------------------------------------------------------------------
# The Complete Formula
# ---------------------------------------------------------------------------
st.subheader("The Complete Formula")

st.markdown(
    """
    <div class="result-card">
    The <b>Alpha Ladder</b> derives Newton's gravitational constant G
    from the fine-structure constant α, the golden ratio φ,
    and the proton-to-electron mass ratio μ -- with <b>zero fitted
    parameters</b>:<br><br>
    <center><code>G = (φ²/2) × (1 + 3α² + (φ/2)α³) × α²¹ × ℏc / mₑ²</code></center><br>
    <b>Every component is derived:</b><br>
    &bull; <b>φ²/2</b> -- bridge coefficient from the vacuum polynomial x² + 6x + 4 = 0<br>
    &bull; <b>3α²</b> -- topological identification: 3 = χ(T<sub>S²</sub>) via HRR = 1/a<sub>1</sub> = d-1 = dim SO(3), all degenerate at n=2. One-loop gives 43/15 ≈ 2.867 (7.1 ppm from 3, testable by next-gen G experiments)<br>
    &bull; <b>(φ/2)α³</b> -- from vacuum polynomial root shift: φ/2 = (r<sub>+</sub>+d)/d, same polynomial as C<sub>0</sub><br>
    &bull; <b>α²¹</b> -- exponent from 6D metric components (d×D - 3 = 21)<br>
    &bull; <b>d=4, D=6</b> -- uniquely selected by 3 independent constraints
    </div>
    """,
    unsafe_allow_html=True,
)

# ---------------------------------------------------------------------------
# Two equivalent forms
# ---------------------------------------------------------------------------
st.subheader("Two Equivalent Forms")

col_a, col_b = st.columns(2)

with col_a:
    st.markdown(
        """
        <div class="theory-box">
        <b>Bridge form (corrected):</b><br><br>
        <code>α_G = φ²/2 × F × α²¹</code><br><br>
        where F = 1 + 3α² + (φ/2)α³<br><br>
        Residual: <b>-0.33 ppm</b> (CODATA 2022)
        </div>
        """,
        unsafe_allow_html=True,
    )

with col_b:
    st.markdown(
        """
        <div class="theory-box">
        <b>Mu-structure form:</b><br><br>
        <code>α_G = α²⁴ × μ × (μ - √φ × (1-α))</code><br><br>
        Same prediction, different factorization<br><br>
        Residual: <b>-0.31 ppm</b> (CODATA 2018)
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Retraction notice
# ---------------------------------------------------------------------------
st.markdown(
    """
    <div class="warning-card">
    <b>v5 Retraction (March 2026):</b> A geometric resummation and derived
    mu prediction proposed in v4 have been retracted -- inconsistent with
    Alighanbari et al. (2025) H2+ spectroscopy (Nature 644, 69) at >14 sigma.
    The G prediction (-0.31 ppm) is unaffected. The framework makes one
    testable prediction (G to sub-ppm), not two.
    </div>
    """,
    unsafe_allow_html=True,
)

# ---------------------------------------------------------------------------
# Navigation hints
# ---------------------------------------------------------------------------
st.subheader("Explore")

_PAGES_DIR = str(Path(__file__).resolve().parent)

nav_col1, nav_col2 = st.columns(2)

with nav_col1:
    st.markdown("**Start Here**")
    st.page_link(_PAGES_DIR + "/13_The_Derivation.py", label="The Derivation -- complete 6D-to-G chain")
    st.page_link(_PAGES_DIR + "/37_Dimension_Uniqueness.py", label="Dimension Uniqueness -- d=4, D=6 from 3 constraints")
    st.page_link(_PAGES_DIR + "/41_Coefficient_Derivation.py", label="Coefficient Derivation -- c_2=3 (HRR + one-loop 43/15), c_3=phi/2")
    st.page_link(_PAGES_DIR + "/46_Vacuum_Polynomial_Derivation.py", label="Vacuum Polynomial -- weakest link, 3 approaches all failed")

    st.markdown("")
    st.markdown("**Gap Analysis**")
    st.page_link(_PAGES_DIR + "/39_Salam_Sezgin_Stabilization.py", label="Salam-Sezgin -- radius stabilized (Gap 1 conditional)")
    st.page_link(_PAGES_DIR + "/40_Charged_Matter_Loops.py", label="Charged Matter Loops -- monopole zeta sign flip")
    st.page_link(_PAGES_DIR + "/27_Corrected_Bridge.py", label="Corrected Bridge -- phi^2/2*(1+3a^2+(phi/2)a^3)")

with nav_col2:
    st.markdown("**Experimental**")
    st.page_link(_PAGES_DIR + "/30_Testable_Predictions.py", label="Testable Predictions -- G to -0.31 ppm, c_2=3 vs 43/15 at 7.1 ppm")
    st.page_link(_PAGES_DIR + "/21_Fifth_Force_Predictions.py", label="Fifth Force -- unobservable (Planck-mass dilaton)")
    st.page_link(_PAGES_DIR + "/15_Solar_System.py", label="Solar System -- Cassini bound, screening resolved")

    st.markdown("")
    st.markdown("**Gauge Physics**")
    st.page_link(_PAGES_DIR + "/43_Braneworld.py", label="Braneworld -- SM must be on 4D brane")
    st.page_link(_PAGES_DIR + "/42_KK_Gauge_Matching.py", label="KK Gauge Matching -- alpha_KK = alpha_EM exactly")
    st.page_link(_PAGES_DIR + "/45_LHC_KK_Comparison.py", label="LHC KK Comparison -- S^2 vs ADD")

    st.markdown("")
    st.markdown("**Legacy**")
    st.page_link(_PAGES_DIR + "/01_Constant_Core.py", label="Constant Core -- original exploratory pages (1-11)")

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption(
    "Alpha Ladder Research Dashboard | "
    "Data sourced from CODATA recommended values | "
    "Computations use Python Decimal with 50-digit precision"
)
