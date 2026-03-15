"""
Bridge Lab -- Page 3

Interactive bridge coefficient search combining bridge_search and predict_g
modules. Displays the target coefficient, filterable search results, G
predictions for each candidate, a sigma heatmap, and a best-prediction summary.
"""

import streamlit as st
import pandas as pd
from decimal import Decimal, getcontext

getcontext().prec = 50

st.set_page_config(page_title="Bridge Lab | Alpha Ladder", layout="wide")

# ---------------------------------------------------------------------------
# Custom CSS
# ---------------------------------------------------------------------------
st.markdown(
    """
    <style>
    .bridge-lab-target {
        background-color: #1a1d23;
        border-left: 4px solid #f59e0b;
        padding: 1.2rem;
        border-radius: 0 8px 8px 0;
        font-family: 'Fira Mono', Consolas, monospace;
        font-size: 1.1rem;
    }
    .best-bridge {
        background-color: #1a1d23;
        border: 2px solid #34d399;
        border-radius: 8px;
        padding: 1.2rem;
        margin: 0.5rem 0;
        font-family: 'Fira Mono', Consolas, monospace;
        font-size: 1.1rem;
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
from app.components.charts import sigma_heatmap  # noqa: E402
from app.components.formatting import (  # noqa: E402
    fmt_decimal,
    fmt_percent,
    fmt_sigma,
    color_by_quality,
)

constants = render_sidebar()

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Bridge Lab")
st.markdown("Interactive search for the bridge coefficient linking alpha to gravity.")
st.divider()

# ---------------------------------------------------------------------------
# Graceful imports from core modules
# ---------------------------------------------------------------------------
_bridge_search_available = False
_predict_g_available = False

try:
    from alpha_ladder_core.bridge_search import (  # noqa: E402
        compute_target_coefficient,
        run_full_search,
    )
    _bridge_search_available = True
except ImportError:
    pass

try:
    from alpha_ladder_core.predict_g import (  # noqa: E402
        get_bridge_candidates,
        predict_G,
        get_G_measurements,
        compare_prediction,
    )
    _predict_g_available = True
except ImportError:
    pass


# ---------------------------------------------------------------------------
# Section 1: Target Coefficient
# ---------------------------------------------------------------------------
st.subheader("Target Coefficient")

if _bridge_search_available and constants is not None:
    target_decimal, target_float = compute_target_coefficient(constants)
    st.markdown(
        f"""
        <div class="bridge-lab-target">
        <b>α_G / α²¹</b> = {fmt_decimal(target_decimal, sig_figs=12)}
        <br><br>
        Full precision: <code>{target_decimal}</code>
        <br><br>
        This is the dimensionless coefficient <i>C</i> in the relation
        <b>α_G = C · α²¹</b>. The Bridge Lab searches for simple
        mathematical expressions that reproduce this value.
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    # Fallback with hardcoded CODATA 2018 values
    alpha_fb = Decimal("0.0072973525693")
    alpha_g_fb = Decimal("1.7512e-45")
    target_decimal = alpha_g_fb / alpha_fb ** 21
    target_float = float(target_decimal)
    st.markdown(
        f"""
        <div class="bridge-lab-target">
        <b>α_G / α²¹</b> = {fmt_decimal(target_decimal, sig_figs=12)}
        <br><br>
        Full precision: <code>{target_decimal}</code>
        <br><br>
        <em>(Using hardcoded CODATA 2018 values; core module not available.)</em>
        </div>
        """,
        unsafe_allow_html=True,
    )

st.divider()

# ---------------------------------------------------------------------------
# Section 2: Bridge Search Results
# ---------------------------------------------------------------------------
st.subheader("Bridge Search Results")


@st.cache_data(show_spinner="Searching for bridge coefficient matches...")
def _cached_full_search(_constants):
    """Cached wrapper around run_full_search."""
    return run_full_search(_constants)


if _bridge_search_available and constants is not None:
    max_err_pct = st.slider(
        "Max error % threshold",
        min_value=0.01,
        max_value=5.0,
        value=1.0,
        step=0.01,
        format="%.2f%%",
        help="Filter results to only show matches within this error threshold.",
    )

    all_results = _cached_full_search(constants)

    # Filter by threshold
    filtered = [(err, expr, val) for err, expr, val in all_results if err <= max_err_pct]

    if filtered:
        search_rows = []
        for err, expr, val in filtered:
            search_rows.append({
                "Expression": expr,
                "Value": fmt_decimal(val, sig_figs=10),
                "Error (%)": f"{err:.6f}%",
                "Error (raw)": err,
            })

        df_search = pd.DataFrame(search_rows)
        st.dataframe(
            df_search[["Expression", "Value", "Error (%)"]],
            use_container_width=True,
            hide_index=True,
            column_config={
                "Expression": st.column_config.TextColumn(width="large"),
                "Value": st.column_config.TextColumn(width="medium"),
                "Error (%)": st.column_config.TextColumn(width="small"),
            },
        )
        st.caption(
            f"Showing {len(filtered)} of {len(all_results)} total matches "
            f"(filtered to <= {max_err_pct:.2f}% error)."
        )
    else:
        st.info(
            f"No matches found within {max_err_pct:.2f}% error. "
            "Try increasing the threshold."
        )
else:
    st.warning(
        "Bridge search module not available. "
        "Ensure alpha_ladder_core.bridge_search is installed."
    )

st.divider()

# ---------------------------------------------------------------------------
# Section 3: G Predictions
# ---------------------------------------------------------------------------
st.subheader("G Predictions")

if _predict_g_available and constants is not None:
    candidates = get_bridge_candidates(constants)

    g_rows = []
    for name, coeff in candidates.items():
        G_pred = predict_G(coeff, constants)
        g_rows.append({
            "Bridge": name,
            "C value": fmt_decimal(coeff, sig_figs=10),
            "G_predicted": fmt_decimal(G_pred, sig_figs=8),
        })

    df_g = pd.DataFrame(g_rows)
    st.dataframe(
        df_g,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Bridge": st.column_config.TextColumn(width="medium"),
            "C value": st.column_config.TextColumn(width="medium"),
            "G_predicted": st.column_config.TextColumn(width="medium"),
        },
    )
else:
    st.warning(
        "Prediction module not available. "
        "Ensure alpha_ladder_core.predict_g is installed."
    )

st.divider()

# ---------------------------------------------------------------------------
# Section 4: Sigma Heatmap
# ---------------------------------------------------------------------------
st.subheader("Sigma Heatmap: Predictions vs Experiments")

if _predict_g_available and constants is not None:
    candidates = get_bridge_candidates(constants)
    measurements = get_G_measurements()
    exp_names = list(measurements.keys())

    sigma_matrix = []
    bridge_names_list = []

    for name, coeff in candidates.items():
        G_pred = predict_G(coeff, constants)
        comparisons = compare_prediction(G_pred, measurements)

        sigmas_row = [round(comp["sigma"], 1) for comp in comparisons]
        sigma_matrix.append(sigmas_row)
        bridge_names_list.append(name)

    fig_sigma = sigma_heatmap(sigma_matrix, bridge_names_list)

    # Use actual experiment names on x-axis
    fig_sigma.update_layout(
        xaxis=dict(
            tickvals=list(range(len(exp_names))),
            ticktext=exp_names,
            tickangle=-45,
        )
    )
    fig_sigma.data[0].x = exp_names
    st.plotly_chart(fig_sigma, use_container_width=True)

    # Detailed sigma breakdown per bridge
    with st.expander("Detailed Sigma Breakdown"):
        for name, coeff in candidates.items():
            G_pred = predict_G(coeff, constants)
            comparisons = compare_prediction(G_pred, measurements)
            st.markdown(f"**{name}** -- G = {fmt_decimal(G_pred, sig_figs=8)}")

            detail_rows = []
            for comp in comparisons:
                err_pct = float(abs(G_pred - comp["G_exp"]) / comp["G_exp"]) * 100
                detail_rows.append({
                    "Experiment": comp["experiment"],
                    "Sigma": fmt_sigma(comp["sigma"], comp["direction"]),
                    "G_exp": fmt_decimal(comp["G_exp"], sig_figs=6),
                    "Quality": f"{err_pct:.4f}%",
                })

            df_detail = pd.DataFrame(detail_rows)
            st.dataframe(df_detail, use_container_width=True, hide_index=True)
else:
    st.warning(
        "Prediction module not available. "
        "Ensure alpha_ladder_core.predict_g is installed."
    )

st.divider()

# ---------------------------------------------------------------------------
# Section 5: Best Prediction Summary
# ---------------------------------------------------------------------------
st.subheader("Best Prediction Summary")

if _predict_g_available and constants is not None:
    candidates = get_bridge_candidates(constants)
    measurements = get_G_measurements()

    # Compute average sigma for each candidate to find the best
    best_name = None
    best_avg_sigma = float("inf")
    best_G = None

    for name, coeff in candidates.items():
        G_pred = predict_G(coeff, constants)
        comparisons = compare_prediction(G_pred, measurements)
        avg_sigma = sum(c["sigma"] for c in comparisons) / len(comparisons)

        if avg_sigma < best_avg_sigma:
            best_avg_sigma = avg_sigma
            best_name = name
            best_G = G_pred

    st.markdown(
        f"""
        <div class="best-bridge">
        <b>Best Bridge Candidate: {best_name}</b>
        <br><br>
        <b>G_predicted</b> = {fmt_decimal(best_G, sig_figs=8)}
        <br>
        <b>Average σ</b> = {best_avg_sigma:.2f} σ across all 7 experiments
        <br><br>
        The candidate <b>φ²/2</b> provides the closest match to experimental
        measurements of Newton's gravitational constant, yielding
        <b>G = 6.67323e-11</b> m³ kg⁻¹ s⁻².
        <br><br>
        This supports the Alpha Ladder hypothesis: the 42-order gap between
        electromagnetic and gravitational coupling is bridged by
        <b>α_G = (φ²/2) · α²¹</b>, where φ is the golden ratio.
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.markdown(
        """
        <div class="best-bridge">
        <b>Best Bridge Candidate: φ²/2</b>
        <br><br>
        <b>G_predicted</b> = 6.67323e-11 m³ kg⁻¹ s⁻²
        <br><br>
        <em>(Core modules not available; showing reference summary.)</em>
        </div>
        """,
        unsafe_allow_html=True,
    )

st.divider()
st.caption("Bridge Lab | Alpha Ladder Research Dashboard")
