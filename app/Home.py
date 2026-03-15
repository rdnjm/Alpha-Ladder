"""
Alpha Ladder Research Dashboard -- Entry Point

Uses st.navigation() for programmatic page grouping.
"""

import streamlit as st
from pathlib import Path

st.set_page_config(
    page_title="Alpha Ladder",
    page_icon=":ladder:",
    layout="wide",
    initial_sidebar_state="expanded",
)

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

_DIR = Path(__file__).resolve().parent / "pages"


def _page(filename, title):
    """Create a st.Page reference."""
    return st.Page(str(_DIR / filename), title=title)


pages = {
    "Home": [
        _page("00_Home_Content.py", "Home"),
    ],
    "Core Results": [
        _page("13_The_Derivation.py", "The Derivation"),
        _page("37_Dimension_Uniqueness.py", "Dimension Uniqueness"),
        _page("34_Feynman_Diagram.py", "Feynman Diagram"),
        _page("36_Mu_Tension.py", "Mu Tension"),
        _page("35_Second_Predictions.py", "Second Predictions"),
        _page("12_The_Prediction.py", "The Prediction"),
    ],
    "Experimental": [
        _page("21_Fifth_Force_Predictions.py", "Fifth Force Predictions"),
        _page("30_Testable_Predictions.py", "Testable Predictions"),
        _page("15_Solar_System.py", "Solar System"),
    ],
    "Theory": [
        _page("14_The_Proof.py", "The Proof"),
        _page("33_One_Alpha_Derivation.py", "One-Alpha Derivation"),
        _page("27_Corrected_Bridge.py", "Corrected Bridge"),
        _page("28_Mu_Structure.py", "Mu Structure"),
        _page("31_Unified_Formula.py", "Unified Formula"),
        _page("25_Bridge_Significance.py", "Bridge Significance"),
        _page("26_Hierarchy_Derivation.py", "Hierarchy Derivation"),
        _page("29_Literature_Comparison.py", "Literature Comparison"),
        _page("32_Residual_Mapping.py", "Residual Mapping"),
    ],
    "Gap Analysis": [
        _page("16_Casimir_Stabilization.py", "Casimir Stabilization"),
        _page("17_Flux_Stabilization.py", "Flux Stabilization"),
        _page("18_Radius_Phenomenology.py", "Radius Phenomenology"),
        _page("19_Chameleon_Screening.py", "Chameleon Screening"),
        _page("20_Dark_Sector_Phenomenology.py", "Dark Sector Phenomenology"),
        _page("22_Radius_Determination.py", "Radius Determination"),
        _page("23_Anomaly_Cancellation.py", "Anomaly Cancellation"),
        _page("24_Cosmological_Constant.py", "Cosmological Constant"),
    ],
    "Legacy Dashboard": [
        _page("01_Constant_Core.py", "Constant Core"),
        _page("02_Geometric_Ladder.py", "Geometric Ladder"),
        _page("03_Bridge_Lab.py", "Bridge Lab"),
        _page("04_Universe_Slider.py", "Universe Slider"),
        _page("05_Phi_Scanner.py", "Phi Scanner"),
        _page("06_Particle_Harmonics.py", "Particle Harmonics"),
        _page("07_Rung_Spacing.py", "Rung Spacing"),
        _page("08_Dilaton_Lab.py", "Dilaton Lab"),
        _page("09_Experimental.py", "Experimental"),
        _page("10_Alpha_Units.py", "Alpha Units"),
        _page("11_Dark_Sector.py", "Dark Sector"),
    ],
}

pg = st.navigation(pages)
pg.run()
