"""
Shared sidebar component for the Alpha Ladder dashboard.

Renders the CODATA edition selector and the 4 Pillars status panel.
All pages should call render_sidebar() to get the active constants namespace.
"""

import streamlit as st

# ---------------------------------------------------------------------------
# Graceful imports from the core library (may not exist yet during
# concurrent development). Fall back to stub data so the UI remains
# functional even before the backend modules are merged.
# ---------------------------------------------------------------------------
try:
    from alpha_ladder_core.constants import (
        get_constants,
        available_editions,
        DEFAULT_EDITION,
    )
    _CORE_AVAILABLE = True
except ImportError:
    _CORE_AVAILABLE = False

    DEFAULT_EDITION = "CODATA 2014"

    def available_editions():
        """Stub: return placeholder editions until the core module is ready."""
        return ["CODATA 2018", "CODATA 2014"]

    def get_constants(edition):
        """Stub: return None; pages must guard against this."""
        return None


# ---------------------------------------------------------------------------
# Pillar definitions
# ---------------------------------------------------------------------------
_PILLARS = [
    {
        "number": 1,
        "title": "Geometric Ladder",
        "description": "αⁿ structure verified",
        "icon_done": "check",
        "icon_pending": "hourglass_flowing_sand",
    },
    {
        "number": 2,
        "title": "The Bridge",
        "description": "α_G = C · α²¹",
        "icon_done": "check",
        "icon_pending": "hourglass_flowing_sand",
    },
    {
        "number": 3,
        "title": "The Dark Sector",
        "description": "α¹⁰ candidate",
        "icon_done": "check",
        "icon_pending": "hourglass_flowing_sand",
    },
    {
        "number": 4,
        "title": "The Prediction",
        "description": "G from first principles",
        "icon_done": "check",
        "icon_pending": "hourglass_flowing_sand",
    },
]


def _pillar_status(pillar_number):
    """Determine whether a pillar should show as verified or pending.

    Currently pillars 1 and 2 are considered verified (the geometric
    ladder and bridge coefficient are computed); pillars 3 and 4 are
    still pending further analysis.

    Returns True (verified) or False (pending).
    """
    return pillar_number in (1, 2, 3, 4)


_GLOBAL_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Fira+Mono:wght@400;500&display=swap');

/* Global base font size */
.stApp, .stMarkdown, .stMarkdown p, .stText,
div[data-testid="stExpander"] p,
.stAlert p {
    font-size: 1.1rem !important;
    line-height: 1.7 !important;
}

/* Markdown containers -- body text, lists, spans (NOT headings) */
div[data-testid="stMarkdownContainer"] p,
div[data-testid="stMarkdownContainer"] li,
div[data-testid="stMarkdownContainer"] span:not(h1 span):not(h2 span):not(h3 span),
div[data-testid="stMarkdownContainer"] b,
div[data-testid="stMarkdownContainer"] strong,
div[data-testid="stMarkdownContainer"] em {
    font-size: 1.15rem !important;
    line-height: 1.75 !important;
}

/* Headings -- preserve native sizes */
h1, div[data-testid="stMarkdownContainer"] h1 { font-size: 2.2rem !important; }
h2, div[data-testid="stMarkdownContainer"] h2 { font-size: 1.75rem !important; }
h3, div[data-testid="stMarkdownContainer"] h3 { font-size: 1.4rem !important; }
h4, div[data-testid="stMarkdownContainer"] h4 { font-size: 1.2rem !important; }

/* Streamlit heading containers */
div[data-testid="stHeading"] h1 { font-size: 2.2rem !important; }
div[data-testid="stHeading"] h2 { font-size: 1.75rem !important; }
div[data-testid="stHeading"] h3 { font-size: 1.4rem !important; }

/* Metric labels (often contain equations with Unicode) */
div[data-testid="stMetricLabel"],
div[data-testid="stMetricLabel"] p,
div[data-testid="stMetricLabel"] label,
div[data-testid="stMetricLabel"] div {
    font-size: 1.05rem !important;
    line-height: 1.5 !important;
}

.stMetric .metric-container {
    background-color: #1a1d23;
    border: 1px solid #2e3440;
    border-radius: 8px;
    padding: 1rem;
}
.stMetric label {
    font-family: 'Fira Mono', monospace;
    font-size: 1.05rem !important;
}
div[data-testid="stMetricValue"] {
    font-family: 'Fira Mono', monospace;
    font-size: 1.8rem !important;
}

/* Inline code / equation blocks */
code, .stCode, pre,
div[data-testid="stMarkdownContainer"] code {
    font-size: 1.15rem !important;
    line-height: 1.6 !important;
    padding: 0.15em 0.4em !important;
}

/* Monospace formula displays in custom HTML */
.formula, .equation, .geom-verify code,
.target-box code, .best-card code,
.decomp-card code, .converter-result code {
    font-size: 1.15rem !important;
}

/* Superscripts and subscripts in equations */
sup { font-size: 0.75em !important; }
sub { font-size: 0.75em !important; }

/* Custom HTML equation/formula containers */
.target-box, .best-card, .decomp-card,
.theory-box, .nav-hint, .verdict-box,
.decomp-card, .connection-card, .gap-display,
.const-table, .phi-table, .unit-table,
.rung-table, .spacing-result, .harmonics-summary {
    font-size: 1.1rem !important;
    line-height: 1.7 !important;
}

/* Sidebar text */
section[data-testid="stSidebar"] div[data-testid="stMarkdownContainer"] p,
section[data-testid="stSidebar"] div[data-testid="stMarkdownContainer"] span,
section[data-testid="stSidebar"] .stMarkdown p {
    font-size: 1.05rem !important;
    line-height: 1.6 !important;
}
</style>
"""


def render_sidebar():
    """Render the shared sidebar and return the selected constants namespace.

    Returns
    -------
    object or None
        The constants namespace for the chosen CODATA edition, or None
        if the core module is not yet available.
    """
    # Inject global styles on every page
    st.markdown(_GLOBAL_CSS, unsafe_allow_html=True)

    with st.sidebar:
        st.header("Alpha Ladder")

        # -- Core module availability notice --
        if not _CORE_AVAILABLE:
            st.warning(
                "Core library not yet available. "
                "Displaying placeholder data."
            )

        # -- CODATA Edition selector --
        editions = available_editions()
        default_idx = 0
        if DEFAULT_EDITION in editions:
            default_idx = editions.index(DEFAULT_EDITION)

        selected_edition = st.selectbox(
            "CODATA Edition",
            editions,
            index=default_idx,
            key="codata_edition",
        )

        st.divider()

        # -- 4 Pillars status panel --
        st.subheader("4 Pillars")

        for pillar in _PILLARS:
            verified = _pillar_status(pillar["number"])
            if verified:
                status_icon = ":white_check_mark:"
                status_label = "Verified"
            else:
                status_icon = ":hourglass_flowing_sand:"
                status_label = "Pending"

            st.markdown(
                f"**{status_icon} Pillar {pillar['number']}: {pillar['title']}**"
            )
            st.caption(f"{pillar['description']} -- {status_label}")

        st.divider()

        # -- PDF Export --
        _pdf_ok = False
        try:
            from app.components.pdf_export import generate_pdf as _gen_pdf
            _pdf_ok = True
        except ImportError:
            try:
                import importlib
                _mod = importlib.import_module("app.components.pdf_export")
                _gen_pdf = _mod.generate_pdf
                _pdf_ok = True
            except (ImportError, AttributeError):
                pass

        if _pdf_ok:
            _c = get_constants(selected_edition)
            if st.button("Download PDF Report", key="sidebar_pdf_btn"):
                with st.spinner("Generating..."):
                    try:
                        _pdf_bytes = _gen_pdf(constants=_c)
                        st.download_button(
                            label="Save PDF",
                            data=_pdf_bytes,
                            file_name="alpha_ladder_report.pdf",
                            mime="application/pdf",
                            key="sidebar_pdf_dl",
                        )
                    except Exception as _exc:
                        st.error(f"PDF failed: {_exc}")

        st.divider()
        st.caption("Alpha Ladder Research Dashboard")

    # -- Fetch and return constants for the selected edition --
    constants = get_constants(selected_edition)
    return constants
