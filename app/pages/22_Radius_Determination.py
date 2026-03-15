"""
Radius Determination -- Page 22

Proves WHY the internal radius a_0 is a free parameter (scaling symmetry
of the 6D EH+GB action for n=2) and catalogs mechanisms that could fix it.
"""

import streamlit as st
import pandas as pd

st.set_page_config(page_title="Radius Determination | Alpha Ladder", layout="wide")

# ---------------------------------------------------------------------------
# Custom CSS (matches Pages 16-21)
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

constants = render_sidebar()

# ---------------------------------------------------------------------------
# Try to load core module
# ---------------------------------------------------------------------------
_core_available = False
try:
    from alpha_ladder_core.radius_determination import (  # noqa: E402
        summarize_radius_determination,
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
        return summarize_radius_determination()

    try:
        summary = _get_summary()
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Radius Determination")
st.markdown(
    "Why the internal radius a_0 is NOT determined by the minimal 6D EH+GB "
    "framework, and what mechanisms could fix it."
)
st.divider()

# ---------------------------------------------------------------------------
# A. The Scaling Symmetry
# ---------------------------------------------------------------------------
st.header("A. The Scaling Symmetry")

if summary:
    scaling = summary["scaling_symmetry"]
    proof_steps = scaling["proof_steps"]

    col_a1, col_a2, col_a3 = st.columns(3)

    with col_a1:
        st.metric(
            label="EH Exponent (n-2)",
            value=str(scaling["eh_scaling_exponent"]),
            delta="Scale-invariant" if scaling["eh_is_scale_invariant"] else "Breaks symmetry",
        )

    with col_a2:
        st.metric(
            label="GB Term",
            value="Topological" if scaling["gb_is_topological"] else "Metric-dependent",
        )

    with col_a3:
        st.metric(
            label="a_0 Status",
            value="FREE" if scaling["is_scale_invariant"] else "FIXED",
        )

    st.markdown("")

    for step in proof_steps:
        st.markdown(
            f"""
            <div class="step-card">{step}</div>
            """,
            unsafe_allow_html=True,
        )

    st.markdown("")
    st.info(scaling["interpretation"])
else:
    st.markdown(
        """
        <div class="theorem-card">
        <b>Scaling symmetry proof:</b><br><br>
        Under g_ab -> t^2 g_ab, the EH action scales as t^{n-2}.
        For n=2, the exponent is 0: scale-invariant.  The GB term
        in 2D equals 4*pi*chi (Gauss-Bonnet theorem): topological.
        Therefore a_0 is a flat direction of the tree-level action.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# B. What Could Fix a_0
# ---------------------------------------------------------------------------
st.header("B. What Could Fix a_0")

if summary:
    mechanisms = summary["mechanisms"]["mechanisms"]

    table_rows = []
    for m in mechanisms:
        table_rows.append({
            "Mechanism": m["name"],
            "Status": m["status"].replace("_", " ").title(),
            "Computed": "Yes" if m["computed"] else "No",
            "Summary": m["result_summary"][:100] + "..." if len(m["result_summary"]) > 100 else m["result_summary"],
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
        <b>Mechanisms that could fix a_0:</b><br><br>
        1. Flux + Casimir balance<br>
        2. Matter loop corrections<br>
        3. Non-perturbative effects<br>
        4. Anthropic selection<br>
        5. Higher-derivative corrections<br><br>
        None are available within the minimal EH+GB framework.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# C. Flux-Casimir Balance
# ---------------------------------------------------------------------------
st.header("C. Flux-Casimir Balance")

if summary:
    balance = summary["flux_casimir_balance"]

    col_c1, col_c2 = st.columns(2)

    with col_c1:
        st.metric(
            label="a_0 Solution",
            value="None found" if balance["a_0_solution"] is None else f"{balance['a_0_solution']:.4f}",
        )

    with col_c2:
        st.metric(
            label="V_min(a_0) Monotonic",
            value="Yes" if balance["is_monotonic"] else "No",
        )

    st.markdown("")

    st.markdown(
        f"""
        <div class="warning-card">
        <b>Result:</b> {balance["honest_result"]}
        </div>
        """,
        unsafe_allow_html=True,
    )

    # Show V_min values table
    v_rows = []
    for entry in balance["V_min_values"]:
        v_rows.append({
            "a_0 (Planck)": f"{entry['a_0']:.1f}",
            "V_min": f"{entry['V_min']:.6e}" if entry["V_min"] is not None else "N/A",
            "sigma_0": f"{entry['sigma_0']:.4f}" if entry["sigma_0"] is not None else "N/A",
        })

    st.dataframe(
        pd.DataFrame(v_rows),
        use_container_width=True,
        hide_index=True,
    )
else:
    st.markdown(
        """
        <div class="warning-card">
        <b>Flux-Casimir balance:</b><br><br>
        V_min(a_0) is monotonically decreasing (scales as a_0^{-4}).
        No finite a_0 minimizes the vacuum energy.  The internal radius
        remains undetermined by flux + Casimir alone.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# D. Radius Landscape
# ---------------------------------------------------------------------------
st.header("D. Radius Landscape")

if summary:
    landscape = summary["landscape"]

    table_rows = []
    for entry in landscape["landscape"]:
        table_rows.append({
            "N": entry["N"],
            "sigma_0": f"{entry['sigma_0']:.4f}" if entry["sigma_0"] is not None else "N/A",
            "m_phi (eV)": f"{entry['m_phi_eV']:.4e}" if entry["m_phi_eV"] is not None else "N/A",
            "Mass Scale": entry.get("mass_scale", "N/A"),
            "Classification": entry["classification"].replace("_", " ").title(),
        })

    st.dataframe(
        pd.DataFrame(table_rows),
        use_container_width=True,
        hide_index=True,
    )

    st.markdown("")

    # Landscape chart
    try:
        from app.components.charts import radius_landscape_chart  # noqa: E402
        fig = radius_landscape_chart(landscape)
        st.plotly_chart(fig, use_container_width=True)
    except ImportError:
        st.info("Chart function radius_landscape_chart not yet available.")
    except Exception as exc:
        st.warning(f"Landscape chart error: {exc}")
else:
    st.markdown(
        """
        <div class="formula-card">
        <b>Landscape (a_0 = l_Pl):</b><br><br>
        At Planck-scale a_0, ALL flux quanta N = 1..10 produce Planck-scale
        dilaton masses.  The dilaton is invisible at all experimental scales.
        To get sub-Planck masses requires a_0 >> l_Pl.
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
    st.markdown(
        f"""
        <div class="proof-card">
        <b>Overall:</b><br><br>
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
        The internal radius a_0 is NOT determined by the minimal 6D EH+GB
        framework.  This is proven by a scaling symmetry: the EH action is
        scale-invariant for n=2 and the GB term is topological.  Five
        mechanisms are cataloged that could break this symmetry, but none
        are available within the minimal framework.  The framework makes
        conditional predictions: IF a_0 is in the sub-mm range, specific
        signatures appear in fifth-force experiments.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Radius Determination | Alpha Ladder Research Dashboard")
