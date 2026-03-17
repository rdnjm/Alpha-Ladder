"""
Phi Scanner -- Page 5

Phi coincidence scanner for Standard Model parameters.
Shows phi-expression matches, Weinberg angle deep dive,
the 21 connection, and summary metrics.
"""

import streamlit as st


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from app.components.sidebar import render_sidebar  # noqa: E402
from app.components.formatting import (  # noqa: E402
    fmt_decimal,
    fmt_percent,
    color_by_quality,
)

constants = render_sidebar()

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Phi Scanner")
st.markdown(
    "Scanning Standard Model parameters for φ-expression coincidences."
)
st.divider()

# ---------------------------------------------------------------------------
# Graceful imports from core modules
# ---------------------------------------------------------------------------
_core_available = False
try:
    from alpha_ladder_core.phi_coincidence import (  # noqa: E402
        run_full_scan,
        deep_dive_weinberg,
        compute_21_connection,
    )
    _core_available = True
except ImportError:
    st.error(
        "Could not import alpha_ladder_core.phi_coincidence. "
        "Ensure the core library is installed."
    )
    st.stop()

import pandas as pd  # noqa: E402

# ---------------------------------------------------------------------------
# Build constants dict for core functions
# ---------------------------------------------------------------------------
if constants is not None and hasattr(constants, "alpha"):
    _constants_arg = constants
else:
    # Fallback: use a dict with CODATA 2018 alpha
    _constants_arg = {"alpha": 0.0072973525693}


# ---------------------------------------------------------------------------
# Section 1: SM Parameter Scan
# ---------------------------------------------------------------------------
st.subheader("SM Parameter Scan")
st.markdown(
    "Each Standard Model parameter is tested against φ-expressions. "
    "Rows are color-coded by match quality: "
    "green (< 0.5%), yellow (< 2%), red (worse)."
)


@st.cache_data
def _cached_full_scan(_alpha_val):
    """Cache the full scan keyed on alpha value."""
    c = {"alpha": _alpha_val}
    return run_full_scan(c)


_alpha_float = (
    float(constants.alpha)
    if constants is not None and hasattr(constants, "alpha")
    else 0.0072973525693
)
scan_results = _cached_full_scan(_alpha_float)

# Build styled dataframe
rows = []
for r in scan_results:
    err = r["best_err"]
    rows.append({
        "Parameter": r["param_name"],
        "Value": f"{r['param_val']:.8g}",
        "Best Phi-Expression": r["best_expr"],
        "Expression Value": f"{r['best_val']:.8g}",
        "Error %": f"{err:.4f}",
        "_err_raw": err,
    })

df_scan = pd.DataFrame(rows)


# Display the scan table without the internal _err_raw column
df_display = df_scan.drop(columns=["_err_raw"])

# Use Streamlit dataframe with column config
st.dataframe(
    df_display,
    use_container_width=True,
    hide_index=True,
    column_config={
        "Parameter": st.column_config.TextColumn(width="large"),
        "Value": st.column_config.TextColumn(width="medium"),
        "Best Phi-Expression": st.column_config.TextColumn(width="large"),
        "Expression Value": st.column_config.TextColumn(width="medium"),
        "Error %": st.column_config.TextColumn(width="small"),
    },
)

# Color-coded legend using HTML for quality indication
st.markdown(
    """
    <div style="display: flex; gap: 1.5rem; font-size: 1.0rem; margin-top: 0.25rem;">
        <span style="color: #00c853;">&#9632; &lt; 0.5% error</span>
        <span style="color: #ffd600;">&#9632; &lt; 2% error</span>
        <span style="color: #d50000;">&#9632; &ge; 2% error</span>
    </div>
    """,
    unsafe_allow_html=True,
)

# Additionally show a colored HTML table for richer visual feedback
with st.expander("Color-coded detail view"):
    html_rows = ""
    for r in scan_results:
        err = r["best_err"]
        if err < 0.5:
            row_color = "rgba(0, 200, 83, 0.18)"
        elif err < 2.0:
            row_color = "rgba(255, 214, 0, 0.18)"
        else:
            row_color = "rgba(213, 0, 0, 0.18)"

        html_rows += (
            f'<tr style="background-color: {row_color};">'
            f'<td style="padding: 6px 10px;">{r["param_name"]}</td>'
            f'<td style="padding: 6px 10px; font-family: monospace; font-size: 1.1rem;">{r["param_val"]:.8g}</td>'
            f'<td style="padding: 6px 10px;">{r["best_expr"]}</td>'
            f'<td style="padding: 6px 10px; font-family: monospace; font-size: 1.1rem;">{r["best_val"]:.8g}</td>'
            f'<td style="padding: 6px 10px; font-family: monospace; font-size: 1.1rem;">{err:.4f}%</td>'
            f"</tr>"
        )

    st.markdown(
        f"""
        <table style="width: 100%; border-collapse: collapse; font-size: 1.0rem;">
        <thead>
            <tr style="border-bottom: 2px solid #3b4252;">
                <th style="padding: 8px 10px; text-align: left;">Parameter</th>
                <th style="padding: 8px 10px; text-align: left;">Value</th>
                <th style="padding: 8px 10px; text-align: left;">Best Phi-Expression</th>
                <th style="padding: 8px 10px; text-align: left;">Expression Value</th>
                <th style="padding: 8px 10px; text-align: left;">Error %</th>
            </tr>
        </thead>
        <tbody>
            {html_rows}
        </tbody>
        </table>
        """,
        unsafe_allow_html=True,
    )

st.divider()

# ---------------------------------------------------------------------------
# Section 2: Deep Dive -- Weinberg Angle
# ---------------------------------------------------------------------------
st.subheader("Deep Dive: Weinberg Angle")

with st.expander("sin²(θ_W) candidates ranked by error", expanded=False):
    weinberg_results = deep_dive_weinberg(_constants_arg)

    wb_rows = []
    for err_pct, name, val, inside in weinberg_results:
        wb_rows.append({
            "Candidate": name,
            "Value": f"{val:.8f}",
            "Error %": f"{err_pct:.4f}",
            "Within Error Bar": "Yes" if inside else "No",
        })

    df_wb = pd.DataFrame(wb_rows)
    st.dataframe(
        df_wb,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Candidate": st.column_config.TextColumn(width="large"),
            "Value": st.column_config.TextColumn(width="medium"),
            "Error %": st.column_config.TextColumn(width="small"),
            "Within Error Bar": st.column_config.TextColumn(width="small"),
        },
    )

    st.markdown(
        "**Reference value:** sin²(θ_W) = 0.23122 (MS-bar at M_Z). "
        "Error bar containment uses Δ < 0.00004."
    )

    # Highlight the best candidate
    if weinberg_results:
        best = weinberg_results[0]
        color = color_by_quality(best[0])
        st.markdown(
            f'<div class="weinberg-detail">'
            f"<b>Best candidate:</b> {best[1]} = {best[2]:.8f}<br>"
            f'<b>Error:</b> <span style="color: {color};">{best[0]:.4f}%</span><br>'
            f'<b>Within error bar:</b> {"Yes" if best[3] else "No"}'
            f"</div>",
            unsafe_allow_html=True,
        )

st.divider()

# ---------------------------------------------------------------------------
# Section 3: The 21 Connection
# ---------------------------------------------------------------------------
st.subheader("The 21 Connection")
st.markdown(
    "Connections between the number 21, Riemann tensor components, "
    "and Standard Model gauge group dimensions."
)

connection = compute_21_connection()

col1, col2 = st.columns(2)

with col1:
    st.markdown("**Riemann Tensor & Combinatorics**")
    st.markdown(
        f"""
        <div class="connection-card">
        <table style="font-family: monospace; font-size: 1.1rem; width: 100%;">
            <tr><td>Riemann components (4D)</td><td style="text-align:right;">{connection['riemann_4d']}</td></tr>
            <tr><td>Triangular number T(6)</td><td style="text-align:right;">{connection['triangular_6']}</td></tr>
            <tr><td>C(7, 2)</td><td style="text-align:right;">{connection['C_7_2']}</td></tr>
            <tr><td>SO(7) dimension</td><td style="text-align:right;">{connection['SO7_dim']}</td></tr>
            <tr><td>Symmetric 6x6 components</td><td style="text-align:right;">{connection['symmetric_6x6']}</td></tr>
        </table>
        </div>
        """,
        unsafe_allow_html=True,
    )

with col2:
    st.markdown("**SM Gauge Group Dimensions**")
    sm = connection["sm_gauge_dims"]
    st.markdown(
        f"""
        <div class="connection-card">
        <table style="font-family: monospace; font-size: 1.1rem; width: 100%;">
            <tr><td>SU(3)</td><td style="text-align:right;">{sm['SU(3)']}</td></tr>
            <tr><td>SU(2)</td><td style="text-align:right;">{sm['SU(2)']}</td></tr>
            <tr><td>U(1)</td><td style="text-align:right;">{sm['U(1)']}</td></tr>
            <tr><td><b>Total gauge bosons</b></td><td style="text-align:right;"><b>{sm['total_gauge_bosons']}</b></td></tr>
        </table>
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("**Symmetric Tensor Components by Dimension**")
    sym = connection["symmetric_tensor_by_dim"]
    sym_html = ""
    for d, comps in sorted(sym.items()):
        highlight = ' style="color: #f59e0b; font-weight: bold;"' if comps == 21 else ""
        sym_html += f"<tr><td>D = {d}</td><td style=\"text-align:right;\"{highlight}>{comps}</td></tr>"

    st.markdown(
        f"""
        <div class="connection-card">
        <table style="font-family: monospace; font-size: 1.1rem; width: 100%;">
            {sym_html}
        </table>
        </div>
        """,
        unsafe_allow_html=True,
    )

st.divider()

# ---------------------------------------------------------------------------
# Section 4: Summary Metrics
# ---------------------------------------------------------------------------
st.subheader("Summary Metrics")

total_params = len(scan_results)
sub_05 = sum(1 for r in scan_results if r["best_err"] < 0.5)
sub_01 = sum(1 for r in scan_results if r["best_err"] < 0.1)

col_m1, col_m2, col_m3 = st.columns(3)

with col_m1:
    st.metric("Total Parameters Scanned", total_params)

with col_m2:
    st.metric("Sub-0.5% Matches", sub_05, delta=f"{sub_05/total_params*100:.0f}%")

with col_m3:
    st.metric("Sub-0.1% Matches", sub_01, delta=f"{sub_01/total_params*100:.0f}%")

st.divider()
st.caption("Phi Scanner | Alpha Ladder Research Dashboard")
