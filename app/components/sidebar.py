"""
Shared sidebar component for the Alpha Ladder dashboard.

Renders the CODATA edition selector, global CSS, and PDF export.
All pages should call render_sidebar() to get the active constants namespace.
"""

import streamlit as st

from app.components.styles import inject_global_css

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
        return ["CODATA 2022", "CODATA 2018", "CODATA 2014"]

    def get_constants(edition):
        """Stub: return None; pages must guard against this."""
        return None


def render_sidebar():
    """Render the shared sidebar and return the selected constants namespace.

    Returns
    -------
    object or None
        The constants namespace for the chosen CODATA edition, or None
        if the core module is not yet available.
    """
    # Inject global styles on every page
    inject_global_css()

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
        st.caption("Alpha Ladder Research Dashboard")

    # -- Fetch and return constants for the selected edition --
    constants = get_constants(selected_edition)
    return constants
