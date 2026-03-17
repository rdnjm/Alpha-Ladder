"""
Cosmological Constant -- Page 24

The flux-stabilized vacuum energy is O(1) Planck. The observed Lambda
is 10^{-122}. This is the universal CC problem. Honest assessment:
no resolution exists within this or any other known framework.
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
    "\u039b is ≈10\u207b\u00b9\u00b2\u00b2. This ≈122-order discrepancy is the universal "
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
        V_min ≈ O(1) M_Pl⁴ for flux-stabilized potential with N=1, a₀ = l_Pl.
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
            label="Λ_obs (Planck)",
            value="≈10⁻\u00b9\u00b2\u00b2",
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
        V_min ≈ O(1) M_Pl⁴ from flux stabilization.<br>
        Λ_obs ≈ 2.888 × 10⁻\u00b9\u00b2\u00b2 M_Pl⁴.<br>
        Ratio: ≈10\u00b9\u00b2\u00b2.<br><br>
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
        2. Anthropic/Landscape: Insufficient vacua (≈10, not 10⁵⁰⁰)<br>
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
        No scalar field adjustment mechanism can dynamically relax Λ
        to a small value if: (1) the potential is smooth, (2) the vacuum is
        Lorentz-invariant, (3) kinetic terms are standard. All three
        conditions hold for the dilaton σ in the Alpha Ladder framework.
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
        reduce it to 10⁻\u00b9\u00b2\u00b2.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Cosmological Constant | Alpha Ladder Research Dashboard")
