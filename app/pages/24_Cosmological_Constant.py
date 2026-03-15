"""
Cosmological Constant -- Page 24

The flux-stabilized vacuum energy is O(1) Planck. The observed Lambda
is 10^{-122}. This is the universal CC problem. Honest assessment:
no resolution exists within this or any other known framework.
"""

import streamlit as st
import pandas as pd


# ---------------------------------------------------------------------------
# Custom CSS (matches Pages 16-23)
# ---------------------------------------------------------------------------
st.markdown(
    """
    <style>
    .proof-card {
        background-color: #1a1d23;
        border: 1px solid #2e3440;
        border-radius: 8px;
        padding: 1.2rem;
        margin: 0.5rem 0;
    }
    .formula-card {
        background-color: #1a1d23;
        border-left: 3px solid #f59e0b;
        padding: 0.8rem 1rem;
        margin: 0.4rem 0;
        border-radius: 0 8px 8px 0;
    }
    .theorem-card {
        background-color: #1a1d23;
        border-left: 3px solid #34d399;
        padding: 0.8rem 1rem;
        margin: 0.4rem 0;
        border-radius: 0 8px 8px 0;
    }
    .step-card {
        background-color: #1a1d23;
        border-left: 3px solid #60a5fa;
        padding: 0.8rem 1rem;
        margin: 0.4rem 0;
        border-radius: 0 8px 8px 0;
    }
    .warning-card {
        background-color: #1a1d23;
        border-left: 3px solid #f87171;
        padding: 0.8rem 1rem;
        margin: 0.4rem 0;
        border-radius: 0 8px 8px 0;
    }
    .wrap-table {
        width: 100%;
        border-collapse: collapse;
    }
    .wrap-table th {
        background-color: #2e3440;
        color: #e0e0e0;
        padding: 8px 12px;
        text-align: left;
        font-size: 0.85rem;
    }
    .wrap-table td {
        padding: 8px 12px;
        border-top: 1px solid #2e3440;
        word-wrap: break-word;
        white-space: normal;
        font-size: 0.85rem;
        color: #d8dee9;
    }
    .wrap-table tr:hover {
        background-color: #1e2230;
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
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from app.components.sidebar import render_sidebar  # noqa: E402
from app.components.formatting import fmt_decimal  # noqa: E402

constants = render_sidebar()

# ---------------------------------------------------------------------------
# Try to load core module
# ---------------------------------------------------------------------------
_core_available = False
try:
    from alpha_ladder_core.cosmological_constant import (  # noqa: E402
        summarize_cosmological_constant,
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
        return summarize_cosmological_constant()

    try:
        summary = _get_summary()
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Cosmological Constant")
st.markdown(
    "The flux-stabilized vacuum energy is O(1) in Planck units. The observed "
    "Lambda is ~10^{-122}. This ~122-order discrepancy is the universal "
    "cosmological constant problem."
)
st.divider()

# ---------------------------------------------------------------------------
# A. Vacuum Energy from Flux
# ---------------------------------------------------------------------------
st.header("A. Vacuum Energy from Flux")

if summary:
    vacuum = summary["vacuum_energy"]

    col_a1, col_a2, col_a3 = st.columns(3)

    with col_a1:
        v_str = f"{vacuum['V_min_planck']:.4e}" if vacuum["V_min_planck"] is not None else "N/A"
        st.metric(
            label="V_min (Planck units)",
            value=v_str,
        )

    with col_a2:
        st.metric(
            label="Sign",
            value=vacuum["sign"].title(),
        )

    with col_a3:
        log_str = f"{vacuum['log10_magnitude']:.1f}" if vacuum["log10_magnitude"] is not None else "N/A"
        st.metric(
            label="log10(|V_min|)",
            value=log_str,
        )
else:
    st.markdown(
        """
        <div class="formula-card">
        <b>Vacuum energy:</b><br><br>
        V_min ~ O(1) M_Pl^4 for flux-stabilized potential with N=1, a_0 = l_Pl.
        The vacuum energy is at the Planck scale -- not small.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# B. The 122-Order Discrepancy
# ---------------------------------------------------------------------------
st.header("B. The 122-Order Discrepancy")

if summary:
    comparison = summary["comparison"]

    col_b1, col_b2, col_b3 = st.columns(3)

    with col_b1:
        st.metric(
            label="Lambda_obs (Planck)",
            value="~10^{-122}",
        )

    with col_b2:
        orders = comparison.get("discrepancy_orders")
        st.metric(
            label="Discrepancy",
            value=f"~{orders} orders" if orders else "N/A",
        )

    with col_b3:
        st.metric(
            label="Fine-Tuning Required",
            value="Yes" if comparison["is_fine_tuning"] else "No",
        )

    st.markdown("")

    st.markdown(
        f"""
        <div class="warning-card">
        <b>The cosmological constant problem:</b><br><br>
        {comparison["honest_assessment"]}
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.markdown(
        """
        <div class="warning-card">
        <b>The 122-order discrepancy:</b><br><br>
        V_min ~ O(1) M_Pl^4 from flux stabilization.<br>
        Lambda_obs ~ 2.888 x 10^{-122} M_Pl^4.<br>
        Ratio: ~10^{122}.<br><br>
        This is the worst fine-tuning problem in physics.
        It affects ALL known theories of gravity.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# C. Why This Is Universal
# ---------------------------------------------------------------------------
st.header("C. Why This Is Universal")

if summary:
    mechanisms = summary["mechanisms"]

    table_rows = []
    for m in mechanisms["mechanisms"]:
        table_rows.append({
            "Approach": m["name"],
            "Applicable Here": "Yes" if m["applicable_here"] else "No",
            "Resolves Problem": "Yes" if m["resolves_problem"] else "No",
            "Reason": m["reason"],
        })

    df = pd.DataFrame(table_rows)
    st.markdown(
        df.to_html(index=False, classes="wrap-table", escape=False),
        unsafe_allow_html=True,
    )

    st.markdown("")
    st.info(mechanisms["honest_assessment"])
else:
    st.markdown(
        """
        <div class="step-card">
        <b>Known approaches (all fail):</b><br><br>
        1. SUSY: Not applicable (framework is non-SUSY)<br>
        2. Anthropic/Landscape: Insufficient vacua (~10, not 10^500)<br>
        3. Sequestering: Requires absent global symmetries<br>
        4. Self-tuning: Blocked by Weinberg no-go (1989)<br>
        5. Unimodular gravity: Shifts but doesn't solve the problem
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# D. No-Go Theorem (Weinberg 1989)
# ---------------------------------------------------------------------------
st.header("D. No-Go Theorem (Weinberg 1989)")

if summary:
    no_go = summary["no_go"]

    for cond in no_go["conditions"]:
        icon = "Met" if cond["met"] else "Not met"
        st.markdown(
            f"""
            <div class="step-card">
            <b>{cond["name"]}:</b> {icon}<br>
            {cond["reason"]}
            </div>
            """,
            unsafe_allow_html=True,
        )

    st.markdown("")

    st.markdown(
        f"""
        <div class="theorem-card">
        <b>Conclusion:</b><br><br>
        {no_go["conclusion"]}<br><br>
        <i>Reference: {no_go["reference"]}</i>
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.markdown(
        """
        <div class="theorem-card">
        <b>Weinberg's no-go theorem (1989):</b><br><br>
        No scalar field adjustment mechanism can dynamically relax Lambda
        to a small value if: (1) the potential is smooth, (2) the vacuum is
        Lorentz-invariant, (3) kinetic terms are standard.  All three
        conditions hold for the dilaton sigma in the Alpha Ladder framework.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# E. Honest Assessment
# ---------------------------------------------------------------------------
st.header("E. Honest Assessment")

if summary:
    # CC scan chart
    try:
        from app.components.charts import cc_scan_chart  # noqa: E402
        fig = cc_scan_chart(summary["cc_scan"])
        st.plotly_chart(fig, use_container_width=True)
    except ImportError:
        st.info("Chart function cc_scan_chart not yet available.")
    except Exception as exc:
        st.warning(f"CC scan chart error: {exc}")

    st.markdown("")

    st.markdown(
        f"""
        <div class="proof-card">
        <b>Overall assessment:</b><br><br>
        {summary["overall_assessment"]}
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        f"""
        <div class="theorem-card">
        <b>Why this is honest:</b><br><br>
        {summary["honest_assessment"]}
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.markdown(
        """
        <div class="proof-card">
        <b>Honest assessment:</b><br><br>
        The cosmological constant problem is the hardest unsolved problem in
        theoretical physics.  No framework -- not string theory, not loop
        quantum gravity, not any other approach -- has a satisfactory resolution.
        The Alpha Ladder framework joins this universal failure honestly:
        V_min is O(1) Planck, and no mechanism within the framework can
        reduce it to 10^{-122}.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Cosmological Constant | Alpha Ladder Research Dashboard")
