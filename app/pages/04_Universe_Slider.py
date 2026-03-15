"""
Universe Slider -- Page 4

What-if explorer: drag the alpha slider to see how physics changes in
hypothetical universes with different fine-structure constants.
"""

import streamlit as st
import pandas as pd
from decimal import Decimal, getcontext

getcontext().prec = 50


# ---------------------------------------------------------------------------
# Custom CSS
# ---------------------------------------------------------------------------
st.markdown(
    """
    <style>
    .viability-card {
        border-radius: 8px;
        padding: 1rem;
        text-align: center;
        font-weight: bold;
        font-size: 1.1rem;
    }
    .viable-yes {
        background-color: rgba(52, 211, 153, 0.15);
        border: 2px solid #34d399;
        color: #34d399;
    }
    .viable-no {
        background-color: rgba(248, 113, 113, 0.15);
        border: 2px solid #f87171;
        color: #f87171;
    }
    .notes-box {
        background-color: #1a1d23;
        border: 1px solid #2e3440;
        border-radius: 8px;
        padding: 1rem;
        margin: 0.5rem 0;
    }
    .comparison-highlight {
        background-color: #1a1d23;
        border-left: 4px solid #60a5fa;
        padding: 1rem;
        border-radius: 0 8px 8px 0;
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
from app.components.charts import ladder_chart  # noqa: E402
from app.components.formatting import fmt_decimal  # noqa: E402

constants = render_sidebar()

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Universe Slider")
st.markdown(
    "Explore hypothetical universes by adjusting the fine-structure constant. "
    "See how key quantities and viability conditions change in real time."
)
st.divider()

# ---------------------------------------------------------------------------
# Graceful import of universe_slider module
# ---------------------------------------------------------------------------
_core_available = False
try:
    from alpha_ladder_core.universe_slider import recompute_physics  # noqa: E402
    _core_available = True
except ImportError:
    pass

# ---------------------------------------------------------------------------
# Real alpha value (for default slider position and comparison)
# ---------------------------------------------------------------------------
if constants is not None and hasattr(constants, "alpha"):
    real_alpha = constants.alpha
else:
    real_alpha = Decimal("0.0072973525693")

real_alpha_float = float(real_alpha)
real_inv_alpha = float(Decimal(1) / real_alpha)

# ---------------------------------------------------------------------------
# Section 1: Alpha Slider
# ---------------------------------------------------------------------------
st.subheader("Alpha Slider")

alpha_min = 1.0 / 200.0  # 0.005
alpha_max = 1.0 / 85.0   # ~0.01176

alpha_val = st.slider(
    "Fine-structure constant (α)",
    min_value=alpha_min,
    max_value=alpha_max,
    value=real_alpha_float,
    step=0.00001,
    format="%.6f",
    help="Drag to explore hypothetical universes. Real α = 1/137.036...",
)

inv_alpha_display = 1.0 / alpha_val if alpha_val != 0 else float("inf")

col_a, col_b = st.columns(2)
with col_a:
    st.metric("α", f"{alpha_val:.6f}")
with col_b:
    st.metric("1 / α", f"{inv_alpha_display:.2f}")

st.divider()

# ---------------------------------------------------------------------------
# Cached computation
# ---------------------------------------------------------------------------
@st.cache_data(show_spinner="Recomputing physics for this universe...")
def _cached_recompute(alpha_float, _constants):
    """Cached wrapper around recompute_physics."""
    return recompute_physics(alpha_float, _constants)


# ---------------------------------------------------------------------------
# Compute physics for the selected alpha
# ---------------------------------------------------------------------------
if _core_available and constants is not None:
    physics = _cached_recompute(alpha_val, constants)

    # Also compute the real-universe baseline for comparison
    real_physics = _cached_recompute(real_alpha_float, constants)

    # -------------------------------------------------------------------
    # Section 2: Key Quantities (metric cards)
    # -------------------------------------------------------------------
    st.subheader("Key Quantities")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.metric(
            "G (phi^2/2 bridge)",
            fmt_decimal(physics["G_predicted"], sig_figs=6),
            help="G from (phi^2/2) * alpha^21 bridge.",
        )

    with col2:
        st.metric(
            "G (alpha^24 * mu^2)",
            fmt_decimal(physics["G_hierarchy"], sig_figs=6),
            help="G from alpha^24 * mu^2 (zero free parameters, 688 ppm).",
        )

    with col3:
        st.metric(
            "G (CODATA)",
            fmt_decimal(constants.G, sig_figs=6),
            help="Measured Newton's G from CODATA.",
        )

    col4, col5, col6 = st.columns(3)

    with col4:
        st.metric(
            "alpha_G (bridge)",
            fmt_decimal(physics["alpha_g"], sig_figs=6),
            help="Gravitational coupling: (phi^2/2) * alpha^21.",
        )

    with col5:
        st.metric(
            "alpha_G (hierarchy)",
            fmt_decimal(physics["alpha_g_hierarchy"], sig_figs=6),
            help="Gravitational coupling: alpha^24 * mu^2.",
        )

    with col6:
        st.metric(
            "rₑ / a₀",
            f"{fmt_decimal(physics['r_e'], sig_figs=6)} / {fmt_decimal(physics['a_0_predicted'], sig_figs=6)}",
            help="Classical electron radius / Bohr radius.",
        )

    st.divider()

    # -------------------------------------------------------------------
    # Section 3: Viability Assessment
    # -------------------------------------------------------------------
    st.subheader("Viability Assessment")

    alpha_dec = Decimal(str(alpha_val))
    qed_perturbative = alpha_dec < Decimal("0.1")

    col_v1, col_v2, col_v3 = st.columns(3)

    with col_v1:
        chem_class = "viable-yes" if physics["chemistry_viable"] else "viable-no"
        chem_text = "YES" if physics["chemistry_viable"] else "NO"
        st.markdown(
            f'<div class="viability-card {chem_class}">'
            f"Chemistry Viable?<br>{chem_text}</div>",
            unsafe_allow_html=True,
        )

    with col_v2:
        stars_class = "viable-yes" if physics["stars_viable"] else "viable-no"
        stars_text = "YES" if physics["stars_viable"] else "NO"
        st.markdown(
            f'<div class="viability-card {stars_class}">'
            f"Stars Viable?<br>{stars_text}</div>",
            unsafe_allow_html=True,
        )

    with col_v3:
        qed_class = "viable-yes" if qed_perturbative else "viable-no"
        qed_text = "YES" if qed_perturbative else "NO"
        st.markdown(
            f'<div class="viability-card {qed_class}">'
            f"QED Perturbative?<br>{qed_text}</div>",
            unsafe_allow_html=True,
        )

    st.divider()

    # -------------------------------------------------------------------
    # Section 4: Notes
    # -------------------------------------------------------------------
    st.subheader("Notes")

    if physics["notes"]:
        notes_html = "<div class='notes-box'><ul>"
        for note in physics["notes"]:
            notes_html += f"<li>{note}</li>"
        notes_html += "</ul></div>"
        st.markdown(notes_html, unsafe_allow_html=True)
    else:
        st.info("No special observations for this alpha value.")

    st.divider()

    # -------------------------------------------------------------------
    # Section 5: Comparison Table
    # -------------------------------------------------------------------
    st.subheader("Comparison: This Universe vs Our Universe")

    comparison_data = [
        {
            "Quantity": "α",
            "This Universe": fmt_decimal(physics["alpha"], sig_figs=8),
            "Our Universe": fmt_decimal(real_physics["alpha"], sig_figs=8),
        },
        {
            "Quantity": "1 / α",
            "This Universe": f"{float(physics['inv_alpha']):.4f}",
            "Our Universe": f"{float(real_physics['inv_alpha']):.4f}",
        },
        {
            "Quantity": "G (phi^2/2 bridge)",
            "This Universe": fmt_decimal(physics["G_predicted"], sig_figs=8),
            "Our Universe": fmt_decimal(real_physics["G_predicted"], sig_figs=8),
        },
        {
            "Quantity": "G (alpha^24 * mu^2)",
            "This Universe": fmt_decimal(physics["G_hierarchy"], sig_figs=8),
            "Our Universe": fmt_decimal(real_physics["G_hierarchy"], sig_figs=8),
        },
        {
            "Quantity": "alpha_G (bridge)",
            "This Universe": fmt_decimal(physics["alpha_g"], sig_figs=6),
            "Our Universe": fmt_decimal(real_physics["alpha_g"], sig_figs=6),
        },
        {
            "Quantity": "alpha_G (hierarchy)",
            "This Universe": fmt_decimal(physics["alpha_g_hierarchy"], sig_figs=6),
            "Our Universe": fmt_decimal(real_physics["alpha_g_hierarchy"], sig_figs=6),
        },
        {
            "Quantity": "rₑ (electron radius)",
            "This Universe": fmt_decimal(physics["r_e"], sig_figs=6),
            "Our Universe": fmt_decimal(real_physics["r_e"], sig_figs=6),
        },
        {
            "Quantity": "Bohr radius (a₀)",
            "This Universe": fmt_decimal(physics["a_0_predicted"], sig_figs=6),
            "Our Universe": fmt_decimal(real_physics["a_0_predicted"], sig_figs=6),
        },
        {
            "Quantity": "Binding energy ratio",
            "This Universe": fmt_decimal(physics["binding_energy_ratio"], sig_figs=6),
            "Our Universe": fmt_decimal(real_physics["binding_energy_ratio"], sig_figs=6),
        },
        {
            "Quantity": "Dark coupling (α¹⁰)",
            "This Universe": fmt_decimal(physics["dark_coupling"], sig_figs=6),
            "Our Universe": fmt_decimal(real_physics["dark_coupling"], sig_figs=6),
        },
        {
            "Quantity": "EM-gravity gap",
            "This Universe": fmt_decimal(physics["gap"], sig_figs=4),
            "Our Universe": fmt_decimal(real_physics["gap"], sig_figs=4),
        },
        {
            "Quantity": "Chemistry viable?",
            "This Universe": "Yes" if physics["chemistry_viable"] else "No",
            "Our Universe": "Yes" if real_physics["chemistry_viable"] else "No",
        },
        {
            "Quantity": "Stars viable?",
            "This Universe": "Yes" if physics["stars_viable"] else "No",
            "Our Universe": "Yes" if real_physics["stars_viable"] else "No",
        },
    ]

    df_comparison = pd.DataFrame(comparison_data)
    st.dataframe(
        df_comparison,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Quantity": st.column_config.TextColumn(width="medium"),
            "This Universe": st.column_config.TextColumn(width="large"),
            "Our Universe": st.column_config.TextColumn(width="large"),
        },
    )

    st.divider()

    # -------------------------------------------------------------------
    # Section 6: Ladder Visualization
    # -------------------------------------------------------------------
    st.subheader("Geometric Ladder for Current Alpha")

    rungs_data = []
    for n in range(1, 25):
        val = physics["alpha_powers"][n]
        label = ""
        if n == 10:
            label = "Dark Matter Candidate"
        elif n == 21:
            label = "Gravity Floor (Bridge)"
        rungs_data.append({
            "power": n,
            "value": val,
            "label": label,
        })

    fig_ladder = ladder_chart(rungs_data)
    fig_ladder.update_layout(
        title=f"Geometric Ladder: αⁿ  (α = {alpha_val:.6f}, 1/α = {inv_alpha_display:.2f})",
    )
    st.plotly_chart(fig_ladder, use_container_width=True)

else:
    st.warning(
        "Universe slider module not available. "
        "Ensure alpha_ladder_core.universe_slider is installed."
    )

    st.markdown(
        """
        <div class="comparison-highlight">
        The Universe Slider requires the <code>alpha_ladder_core.universe_slider</code>
        module to recompute physics for hypothetical alpha values.
        <br><br>
        When available, this page will display:
        <ul>
        <li>Key quantities (G, α_G, rₑ, Bohr radius) as metric cards</li>
        <li>Viability indicators (chemistry, stars, QED) with color coding</li>
        <li>Detailed notes about the hypothetical universe</li>
        <li>Side-by-side comparison with our universe</li>
        <li>Geometric ladder visualization for the selected alpha</li>
        </ul>
        </div>
        """,
        unsafe_allow_html=True,
    )

st.divider()
st.caption("Universe Slider | Alpha Ladder Research Dashboard")
