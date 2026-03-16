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
# Custom CSS
# ---------------------------------------------------------------------------
st.markdown(
    """
    <style>
    .theory-box {
        background-color: #1a1d23;
        border-left: 4px solid #60a5fa;
        padding: 1.2rem;
        border-radius: 0 8px 8px 0;
        margin: 1rem 0;
    }
    .result-card {
        background-color: #1a1d23;
        border-left: 4px solid #34d399;
        padding: 1.2rem;
        border-radius: 0 8px 8px 0;
        margin: 1rem 0;
    }
    .nav-hint {
        background-color: #1e2128;
        border: 1px solid #3b4252;
        border-radius: 8px;
        padding: 1rem;
        margin: 0.5rem 0;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

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
        label="mu residual",
        value="+0.0004 ppm",
        delta="measurement precision",
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
    from the fine-structure constant alpha, the golden ratio phi,
    and the proton-to-electron mass ratio mu -- with <b>zero fitted
    parameters</b>:<br><br>
    <center><code>G = (phi^2/2) * (1 + 3*alpha^2 + (phi/2)*alpha^3) * alpha^21 * hbar*c / m_e^2</code></center><br>
    <b>Every component is derived:</b><br>
    &bull; <b>phi^2/2</b> -- bridge coefficient from the vacuum polynomial x^2 + 6x + 4 = 0<br>
    &bull; <b>3*alpha^2</b> -- radiative correction, 3 = d-1 = n(n+1)/2 = n+1 at n=2<br>
    &bull; <b>(phi/2)*alpha^3</b> -- NLO correction from S^2 volume cancellation<br>
    &bull; <b>alpha^21</b> -- exponent from 6D metric components (d*D - 3 = 21)<br>
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
        <code>alpha_G = phi^2/2 * F * alpha^21</code><br><br>
        where F = 1 + 3*alpha^2 + (phi/2)*alpha^3<br><br>
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
        <code>alpha_G = alpha^24 * mu * (mu - sqrt(phi)*(1-alpha))</code><br><br>
        Same prediction, different factorization<br><br>
        Residual: <b>-0.31 ppm</b> (CODATA 2018)
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Navigation hints
# ---------------------------------------------------------------------------
st.subheader("Explore")

nav_col1, nav_col2 = st.columns(2)

with nav_col1:
    st.markdown(
        """
        <div class="nav-hint">
        <b>Start Here:</b> <em>The Derivation</em> (Core Results)<br>
        The complete 6D-to-G derivation chain with all steps traced.
        </div>
        """,
        unsafe_allow_html=True,
    )
    st.markdown(
        """
        <div class="nav-hint">
        <b>Key Result:</b> <em>Dimension Uniqueness</em> (Core Results)<br>
        Three constraints uniquely select d=4, D=6.
        </div>
        """,
        unsafe_allow_html=True,
    )

with nav_col2:
    st.markdown(
        """
        <div class="nav-hint">
        <b>Experimental:</b> <em>Fifth Force Predictions</em><br>
        Testable predictions for Eot-Wash and other experiments.
        </div>
        """,
        unsafe_allow_html=True,
    )
    st.markdown(
        """
        <div class="nav-hint">
        <b>Legacy:</b> <em>Constant Core</em> (pages 1-11)<br>
        The original exploratory dashboard pages.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption(
    "Alpha Ladder Research Dashboard | "
    "Data sourced from CODATA recommended values | "
    "Computations use Python Decimal with 50-digit precision"
)
