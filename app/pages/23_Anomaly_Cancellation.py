"""
Anomaly Cancellation -- Page 23

Pure 6D gravity is anomaly-free. Matter coupling requires Green-Schwarz.
The G prediction is UNAFFECTED by anomaly constraints.
"""

import streamlit as st
import pandas as pd


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from app.components.sidebar import render_sidebar  # noqa: E402

constants = render_sidebar()

# ---------------------------------------------------------------------------
# Try to load core module
# ---------------------------------------------------------------------------
_core_available = False
try:
    from alpha_ladder_core.anomaly_cancellation import (  # noqa: E402
        summarize_anomaly_cancellation,
    )
    _core_available = True
except ImportError:
    pass

# ---------------------------------------------------------------------------
# Compute values (guarded, cached)
# ---------------------------------------------------------------------------
summary = None

if _core_available:
    @st.cache_data
    def _get_summary():
        return summarize_anomaly_cancellation()

    try:
        summary = _get_summary()
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Anomaly Cancellation")
st.markdown(
    "Pure 6D gravity is anomaly-free. Coupling to Standard Model matter "
    "requires Green-Schwarz anomaly cancellation. The G prediction is unaffected."
)
st.divider()

# ---------------------------------------------------------------------------
# A. Pure Gravity is Safe
# ---------------------------------------------------------------------------
st.header("A. Pure Gravity is Safe")

if summary:
    anomaly_poly = summary["anomaly_polynomial"]

    col_a1, col_a2, col_a3 = st.columns(3)

    with col_a1:
        st.metric(
            label="Graviton I_8",
            value="0",
            delta="Non-chiral, anomaly-free",
        )

    with col_a2:
        st.metric(
            label="Gravitino I_8",
            value=f"1/{5760}",
            delta="Chiral, anomalous",
        )

    with col_a3:
        st.metric(
            label="Pure Gravity",
            value="SAFE",
        )

    st.markdown("")

    st.markdown(
        f"""
        <div class="theorem-card">
        <b>Result:</b> {anomaly_poly["interpretation"]}
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.markdown(
        """
        <div class="theorem-card">
        <b>Pure gravity in 6D:</b><br><br>
        The graviton is non-chiral in 6D, contributing zero to the anomaly
        polynomial I_8.  Pure gravity (EH + GB) is anomaly-free in any
        dimension.  No Green-Schwarz mechanism is needed for the minimal
        framework.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# B. Matter Requires Green-Schwarz
# ---------------------------------------------------------------------------
st.header("B. Matter Requires Green-Schwarz")

st.markdown(
    """
    <div class="step-card">
    <b>The anomaly cancellation condition:</b><br><br>
    For (1,0) SUGRA in 6D with n_H hypermultiplets, n_V vector multiplets,
    and n_T tensor multiplets:<br>
    <code>n_H - n_V + 29 * n_T = 273</code><br><br>
    The anomaly polynomial I_8 must factorize as X_4 . X_4' for the
    Green-Schwarz mechanism to cancel the anomaly via a tree-level
    B ^ X_4 coupling.
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown("")

if summary:
    gs_checks = summary["green_schwarz_checks"]

    table_rows = []
    for key, check in gs_checks.items():
        table_rows.append({
            "Gauge Group": check["gauge_group"],
            "GS Factorizes": "Yes" if check["factorizes"] else "No",
            "GS Applies": "Yes" if check["gs_mechanism_applies"] else "No",
            "Anomaly Condition": "Met" if check.get("anomaly_condition_met") else "Not met" if check.get("anomaly_condition_met") is False else "N/A",
        })

    st.dataframe(
        pd.DataFrame(table_rows),
        use_container_width=True,
        hide_index=True,
    )
else:
    st.markdown(
        """
        <div class="formula-card">
        <b>Green-Schwarz checks:</b><br><br>
        E8 x E8: Factorizes (anomaly-free)<br>
        SO(32): Factorizes (anomaly-free)<br>
        SU(3) x SU(2) x U(1): Does NOT factorize (anomalous alone)
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# C. Which Groups Work
# ---------------------------------------------------------------------------
st.header("C. Which Groups Work")

if summary:
    group_scan = summary["group_scan"]

    col_c1, col_c2, col_c3 = st.columns(3)

    with col_c1:
        st.metric(
            label="Anomaly-Free Groups",
            value=str(group_scan["n_anomaly_free"]),
        )

    with col_c2:
        st.metric(
            label="Contain SM",
            value=str(group_scan["n_contain_sm"]),
        )

    with col_c3:
        st.metric(
            label="Minimal Group",
            value=group_scan["minimal_group"] or "N/A",
        )

    st.markdown("")

    # Groups table
    table_rows = []
    for g in group_scan["groups"]:
        chain = " \u2192 ".join(g["sm_embedding_chain"]) if g["sm_embedding_chain"] else "N/A"
        table_rows.append({
            "Group": g["name"],
            "Dim": g["dim"],
            "Rank": g["rank"],
            "Anomaly-Free": "Yes" if g["gs_factorizes"] else "No",
            "Contains SM": "Yes" if g["contains_sm"] else "No",
            "Origin": g["origin"],
        })

    st.dataframe(
        pd.DataFrame(table_rows),
        use_container_width=True,
        hide_index=True,
    )

    st.markdown("")

    # Group chart
    try:
        from app.components.charts import anomaly_group_chart  # noqa: E402
        fig = anomaly_group_chart(group_scan)
        st.plotly_chart(fig, use_container_width=True)
    except ImportError:
        st.info("Chart function anomaly_group_chart not yet available.")
    except Exception as exc:
        st.warning(f"Group chart error: {exc}")

    st.markdown("")
    st.info(group_scan["honest_assessment"])
else:
    st.markdown(
        """
        <div class="formula-card">
        <b>Anomaly-free 6D gauge groups containing the SM:</b><br><br>
        E8 × E8 (dim 496, from heterotic string)<br>
        SO(32) (dim 496, from Type I / heterotic)<br>
        E7 × E7 (dim 266, from 6D SUGRA)<br>
        E6 × E7 (dim 211, from 6D SUGRA)<br><br>
        The raw SM group SU(3) × SU(2) × U(1) is NOT anomaly-free in 6D.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# D. Constraints on Alpha Ladder
# ---------------------------------------------------------------------------
st.header("D. Constraints on Alpha Ladder")

if summary:
    constraints = summary["alpha_ladder_constraints"]

    col_d1, col_d2, col_d3 = st.columns(3)

    with col_d1:
        st.metric(
            label="G Prediction",
            value="UNAFFECTED",
        )

    with col_d2:
        st.metric(
            label="Gravity Sector",
            value="Anomaly-Free",
        )

    with col_d3:
        st.metric(
            label="Matter Sector",
            value="Constrained",
        )

    st.markdown("")

    for step in constraints["detailed_reasoning"]:
        st.markdown(
            f"""
            <div class="step-card">{step}</div>
            """,
            unsafe_allow_html=True,
        )

    st.markdown("")

    st.markdown(
        f"""
        <div class="proof-card">
        <b>Honest assessment:</b><br><br>
        {constraints["honest_assessment"]}
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.markdown(
        """
        <div class="proof-card">
        <b>Key result:</b><br><br>
        The G prediction from the φ²/2 bridge derives entirely from the
        gravity sector (EH + GB → KK → Brans-Dicke → G). No gauge fields
        or chiral matter enter the derivation. Anomaly cancellation constrains
        the matter sector but leaves the gravity sector -- and therefore the
        G prediction -- completely unchanged.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Anomaly Cancellation | Alpha Ladder Research Dashboard")
