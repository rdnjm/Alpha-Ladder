"""
Consolidated CSS for the Alpha Ladder dashboard.

Single source of truth for all card styles, typography, and visual polish.
All pages inherit these styles via inject_global_css() called from render_sidebar().
"""

import streamlit as st


def inject_global_css():
    """Inject the global CSS stylesheet into the current Streamlit page."""
    st.markdown(_GLOBAL_CSS, unsafe_allow_html=True)


_GLOBAL_CSS = """
<style>
/* ===== Font imports ===== */
@import url('https://fonts.googleapis.com/css2?family=Crimson+Pro:wght@400;600;700&family=Source+Sans+3:wght@300;400;500;600&family=Fira+Mono:wght@400;500&display=swap');

/* ===== Base typography ===== */
.stApp, .stMarkdown, .stMarkdown p, .stText,
div[data-testid="stExpander"] p,
.stAlert p {
    font-family: 'Source Sans 3', 'Segoe UI', sans-serif;
    font-size: 1.1rem !important;
    line-height: 1.7 !important;
}

/* Markdown body text */
div[data-testid="stMarkdownContainer"] p,
div[data-testid="stMarkdownContainer"] li,
div[data-testid="stMarkdownContainer"] span:not(h1 span):not(h2 span):not(h3 span),
div[data-testid="stMarkdownContainer"] b,
div[data-testid="stMarkdownContainer"] strong,
div[data-testid="stMarkdownContainer"] em {
    font-family: 'Source Sans 3', 'Segoe UI', sans-serif;
    font-size: 1.15rem !important;
    line-height: 1.75 !important;
}

/* ===== Headings -- Crimson Pro serif ===== */
h1, div[data-testid="stMarkdownContainer"] h1,
div[data-testid="stHeading"] h1 {
    font-family: 'Crimson Pro', Georgia, serif !important;
    font-size: 2.2rem !important;
    font-weight: 700 !important;
    letter-spacing: 0.02em !important;
    text-shadow: 0 1px 2px rgba(0,0,0,0.3);
}
h2, div[data-testid="stMarkdownContainer"] h2,
div[data-testid="stHeading"] h2 {
    font-family: 'Crimson Pro', Georgia, serif !important;
    font-size: 1.75rem !important;
    font-weight: 600 !important;
    color: #93c5fd !important;
}
h3, div[data-testid="stMarkdownContainer"] h3,
div[data-testid="stHeading"] h3 {
    font-family: 'Crimson Pro', Georgia, serif !important;
    font-size: 1.4rem !important;
    font-weight: 600 !important;
}
h4, div[data-testid="stMarkdownContainer"] h4 {
    font-family: 'Crimson Pro', Georgia, serif !important;
    font-size: 1.2rem !important;
}

/* ===== Metric styling ===== */
.stMetric .metric-container {
    background-color: #1a1d23;
    border: 1px solid #2e3440;
    border-radius: 10px;
    padding: 1rem;
}
.stMetric label {
    font-family: 'Fira Mono', monospace;
    font-size: 1.05rem !important;
}
div[data-testid="stMetricLabel"],
div[data-testid="stMetricLabel"] p,
div[data-testid="stMetricLabel"] label,
div[data-testid="stMetricLabel"] div {
    font-size: 1.05rem !important;
    line-height: 1.5 !important;
}
div[data-testid="stMetricValue"] {
    font-family: 'Fira Mono', monospace;
    font-size: 2rem !important;
}

/* Metric top accent bar */
div[data-testid="stMetric"] > div {
    border-top: 2px solid #60a5fa;
    border-radius: 10px;
    padding: 0.75rem 1rem;
}

/* ===== Code / monospace ===== */
code, .stCode, pre,
div[data-testid="stMarkdownContainer"] code {
    font-family: 'Fira Mono', Consolas, monospace !important;
    font-size: 1.15rem !important;
    line-height: 1.6 !important;
    padding: 0.15em 0.4em !important;
}
.formula, .equation, .geom-verify code,
.target-box code, .best-card code,
.decomp-card code, .converter-result code {
    font-family: 'Fira Mono', Consolas, monospace !important;
    font-size: 1.15rem !important;
}
sup { font-size: 0.75em !important; }
sub { font-size: 0.75em !important; }

/* ===== Custom HTML containers ===== */
.target-box, .best-card, .decomp-card,
.theory-box, .nav-hint, .verdict-box,
.connection-card, .gap-display,
.const-table, .phi-table, .unit-table,
.rung-table, .spacing-result, .harmonics-summary {
    font-size: 1.1rem !important;
    line-height: 1.7 !important;
}

/* ===== Sidebar text ===== */
section[data-testid="stSidebar"] div[data-testid="stMarkdownContainer"] p,
section[data-testid="stSidebar"] div[data-testid="stMarkdownContainer"] span,
section[data-testid="stSidebar"] .stMarkdown p {
    font-size: 1.05rem !important;
    line-height: 1.6 !important;
}

/* ===== Card classes (consolidated from all pages) ===== */

/* Generic cards */
.proof-card, .deriv-card, .pred-card, .gap-card,
.mechanism-card, .anomaly-card, .cc-card, .hierarchy-card,
.bridge-card, .mu-card, .lit-card, .unified-card,
.residual-card, .alpha-card, .feynman-card, .second-card,
.tension-card, .dim-card, .fix-card, .solar-card,
.casimir-card, .flux-card, .radius-card, .chameleon-card,
.dark-card, .fifth-card, .reconcile-box, .best-fit-card,
.comparison-highlight, .viability-card, .bridge-lab-target,
.best-bridge, .weinberg-detail, .alps-step,
.converter-result, .harmonics-summary, .falsify-card {
    background-color: #1a1d23;
    border: 1px solid #2e3440;
    border-radius: 10px;
    padding: 1.2rem;
    margin: 0.5rem 0;
    box-shadow: 0 2px 8px rgba(0,0,0,0.3);
    transition: border-color 0.2s;
}

/* Formula cards -- amber left border */
.formula-card {
    background: linear-gradient(135deg, #1a1d23 0%, #1e2128 100%);
    border-left: 3px solid #f59e0b;
    border-top: none;
    border-right: none;
    border-bottom: none;
    padding: 0.8rem 1rem;
    margin: 0.4rem 0;
    border-radius: 0 10px 10px 0;
    box-shadow: 0 2px 8px rgba(0,0,0,0.3);
    transition: border-color 0.2s;
}

/* Takeaway cards -- emerald left border */
.takeaway-card {
    background: linear-gradient(135deg, #1a1d23 0%, #1e2128 100%);
    border-left: 3px solid #34d399;
    border-top: none;
    border-right: none;
    border-bottom: none;
    padding: 0.8rem 1rem;
    margin: 0.4rem 0;
    border-radius: 0 10px 10px 0;
    box-shadow: 0 2px 8px rgba(0,0,0,0.3);
    transition: border-color 0.2s;
}

/* Step cards */
.step-card {
    background-color: #1a1d23;
    border: 1px solid #3b4252;
    border-radius: 10px;
    padding: 0.8rem 1rem;
    margin: 0.3rem 0;
    box-shadow: 0 2px 8px rgba(0,0,0,0.3);
    transition: border-color 0.2s;
}

/* Warning cards -- red left border */
.warning-card {
    background-color: #1a1d23;
    border-left: 3px solid #f87171;
    border-top: none;
    border-right: none;
    border-bottom: none;
    padding: 0.8rem 1rem;
    margin: 0.4rem 0;
    border-radius: 0 10px 10px 0;
    box-shadow: 0 2px 8px rgba(0,0,0,0.3);
}

/* Theorem cards -- blue left border */
.theorem-card {
    background: linear-gradient(135deg, #1a1d23 0%, #1e2128 100%);
    border-left: 3px solid #60a5fa;
    border-top: none;
    border-right: none;
    border-bottom: none;
    padding: 0.8rem 1rem;
    margin: 0.4rem 0;
    border-radius: 0 10px 10px 0;
    box-shadow: 0 2px 8px rgba(0,0,0,0.3);
}

/* Navigation hints */
.nav-hint {
    background-color: #1e2128;
    border: 1px solid #3b4252;
    border-radius: 10px;
    padding: 0.6rem 1rem;
    margin: 0.3rem 0;
    font-size: 0.95rem !important;
    box-shadow: 0 2px 8px rgba(0,0,0,0.3);
}

/* Result card (home page) */
.result-card {
    background: linear-gradient(135deg, #1a1d23 0%, #1e2128 100%);
    border: 1px solid #60a5fa;
    border-radius: 10px;
    padding: 1.5rem;
    margin: 1rem 0;
    box-shadow: 0 2px 8px rgba(0,0,0,0.3);
}

/* Wrap table */
.wrap-table table {
    table-layout: fixed;
    width: 100%;
}
.wrap-table td, .wrap-table th {
    word-wrap: break-word;
    overflow-wrap: break-word;
    white-space: normal;
}

/* Viability states */
.viability-card.viable {
    border-left: 3px solid #34d399;
}
.viability-card.marginal {
    border-left: 3px solid #fbbf24;
}
.viability-card.non-viable {
    border-left: 3px solid #f87171;
}

/* ===== Styled dividers ===== */
hr {
    border: none;
    height: 1px;
    background: linear-gradient(90deg, transparent, #3b4252 20%, #60a5fa 50%, #3b4252 80%, transparent);
    margin: 1.5rem 0;
}
</style>
"""
