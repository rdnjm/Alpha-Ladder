"""
Alpha Ladder Research Dashboard -- Home Page

Main entry point for the Streamlit application. Displays the 4 Pillars
overview, a brief explanation of the theory, and navigation hints.
"""

import streamlit as st

st.set_page_config(
    page_title="Alpha Ladder",
    page_icon=":ladder:",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ---------------------------------------------------------------------------
# Custom CSS for a clean, scientific aesthetic
# ---------------------------------------------------------------------------
st.markdown(
    """
    <style>
    /* Page-specific styles for Home (global styles are in sidebar.py) */
    .theory-box {
        background-color: #1a1d23;
        border-left: 4px solid #60a5fa;
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
# Sidebar
# ---------------------------------------------------------------------------
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.components.sidebar import render_sidebar  # noqa: E402

constants = render_sidebar()

# ---------------------------------------------------------------------------
# Title
# ---------------------------------------------------------------------------
st.title("Alpha Ladder Research Dashboard")
st.markdown("**Coupling Constants as Geometric Powers of Alpha**")
st.divider()

# ---------------------------------------------------------------------------
# 4 Pillars as metric cards
# ---------------------------------------------------------------------------
col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric(
        label="Pillar 1: Geometric Ladder",
        value="Verified",
        delta="αⁿ structure",
        delta_color="normal",
    )

with col2:
    st.metric(
        label="Pillar 2: The Bridge",
        value="Verified",
        delta="α_G = C · α²¹",
        delta_color="normal",
    )

with col3:
    st.metric(
        label="Pillar 3: Dark Sector",
        value="Active",
        delta="α¹⁰ candidate",
        delta_color="normal",
    )

with col4:
    st.metric(
        label="Pillar 4: The Prediction",
        value="Active",
        delta="G from first principles",
        delta_color="normal",
    )

st.divider()

# ---------------------------------------------------------------------------
# Theory overview
# ---------------------------------------------------------------------------
st.subheader("The Core Idea")

st.markdown(
    """
    <div class="theory-box">
    The <b>Alpha Ladder</b> hypothesis proposes that the fundamental coupling
    constants of nature are <em>geometric powers of the fine-structure
    constant</em>, α = e<sup>2</sup> / (4π ε₀ ℏ c)
    &approx; 1/137.036.
    <br><br>
    In particular, the gravitational coupling constant α_G is connected
    to α through:
    <br><br>
    <code>α_G = C · α²¹</code>
    <br><br>
    where <b>C</b> is a dimensionless bridge coefficient expressible in
    terms of fundamental mathematical constants (φ, π, e). This
    relationship spans the <b>42 orders of magnitude</b> between
    electromagnetism and gravity -- the so-called "Great Desert" --
    and opens the door to predicting Newton's G from first principles.
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
        <b>Constant Core</b> (Page 1)<br>
        Full-precision CODATA constants, electron geometry verification,
        and the 42-order gap between α and α_G.
        </div>
        """,
        unsafe_allow_html=True,
    )

with nav_col2:
    st.markdown(
        """
        <div class="nav-hint">
        <b>Geometric Ladder</b> (Page 2)<br>
        Interactive log-scale chart of αⁿ, detailed rung table,
        bridge candidates, and residual analysis.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
st.subheader("Export")

st.markdown(
    """
    <div class="nav-hint">
    <b>Download Report</b><br>
    Generate a single PDF containing key results from all 24 dashboard pages,
    including computed metrics, tables, gap status, and honest assessments.
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown("")

_pdf_available = False
try:
    from app.components.pdf_export import generate_pdf  # noqa: E402
    _pdf_available = True
except ImportError:
    pass

if _pdf_available:
    if st.button("Generate PDF Report"):
        with st.spinner("Generating PDF..."):
            try:
                pdf_bytes = generate_pdf(constants=constants)
                st.download_button(
                    label="Download PDF",
                    data=pdf_bytes,
                    file_name="alpha_ladder_report.pdf",
                    mime="application/pdf",
                )
            except Exception as exc:
                st.error(f"PDF generation failed: {exc}")
else:
    st.info("PDF export requires fpdf2: pip install fpdf2")

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption(
    "Alpha Ladder Research Dashboard | "
    "Data sourced from CODATA recommended values | "
    "Computations use Python Decimal with 50-digit precision"
)
