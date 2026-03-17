"""
Rung Spacing -- Page 7

Three-tab layout for rung spacing search: rational spacings,
irrational spacings, and continuous optimization.
"""

import streamlit as st
import math


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from app.components.sidebar import render_sidebar  # noqa: E402
from app.components.formatting import fmt_decimal, color_by_quality  # noqa: E402
from app.components.charts import spacing_score_chart  # noqa: E402

constants = render_sidebar()

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Rung Spacing Search")
st.markdown(
    "Brute force search for the best rung spacing "
    "in the Standard Model mass spectrum."
)
st.divider()

# ---------------------------------------------------------------------------
# Graceful imports
# ---------------------------------------------------------------------------
_core_available = False
try:
    from alpha_ladder_core.rung_spacing import (  # noqa: E402
        compute_rungs,
        search_rational_spacings,
        search_irrational_spacings,
        search_continuous_optimum,
        get_best_fit_details,
        score_spacing,
    )
    _core_available = True
except ImportError:
    st.error(
        "Could not import alpha_ladder_core.rung_spacing. "
        "Ensure the core library is installed."
    )
    st.stop()

import pandas as pd  # noqa: E402

# ---------------------------------------------------------------------------
# Build constants argument and compute rungs
# ---------------------------------------------------------------------------
if constants is not None and hasattr(constants, "alpha"):
    _constants_arg = constants
else:
    _constants_arg = {"alpha": 0.0072973525693}

rung_dict = compute_rungs(_constants_arg)
rung_values = list(rung_dict.values())
n_particles = len(rung_dict)


# ---------------------------------------------------------------------------
# Cached search functions
# ---------------------------------------------------------------------------
@st.cache_data
def _cached_rational(rung_dict_items):
    """Cache rational spacing search. Convert items back to dict internally."""
    rd = dict(rung_dict_items)
    return search_rational_spacings(rd)


@st.cache_data
def _cached_irrational(rung_dict_items):
    """Cache irrational spacing search."""
    rd = dict(rung_dict_items)
    return search_irrational_spacings(rd)


@st.cache_data
def _cached_continuous(rung_vals_tuple):
    """Cache continuous optimization search."""
    return search_continuous_optimum(list(rung_vals_tuple))


# Convert to hashable types for caching
_rung_items = tuple(sorted(rung_dict.items()))
_rung_vals_tuple = tuple(rung_values)

rational_results = _cached_rational(_rung_items)
irrational_results = _cached_irrational(_rung_items)
best_spacing, best_score, local_mins = _cached_continuous(_rung_vals_tuple)

# ---------------------------------------------------------------------------
# Three-tab layout
# ---------------------------------------------------------------------------
tab1, tab2, tab3 = st.tabs([
    "Rational Spacings",
    "Irrational Spacings",
    "Continuous Optimization",
])

# ---------------------------------------------------------------------------
# Tab 1: Rational Spacings
# ---------------------------------------------------------------------------
with tab1:
    st.subheader("Rational Spacings (1/k)")
    st.markdown(
        f"Testing spacings 1/k for k = 1 to 24 against {n_particles} particles."
    )

    # Random expectation: for tolerance_frac=0.15, probability = 2*0.15 = 0.30
    random_prob = 0.30
    random_expected = n_particles * random_prob

    rat_rows = []
    for k, spacing, avg, matches, details in rational_results:
        excess = matches - random_expected
        rat_rows.append({
            "k": k,
            "Spacing (1/k)": f"{spacing:.6f}",
            "Matches": matches,
            "Expected (Random)": f"{random_expected:.1f}",
            "Excess": f"{excess:+.1f}",
            "Avg Error": f"{avg:.4f}",
        })

    df_rat = pd.DataFrame(rat_rows)

    # Highlight the best row(s)
    def _highlight_best_rational(row):
        if int(row["Matches"]) == max(int(r["Matches"]) for r in rat_rows):
            return ["background-color: rgba(0, 200, 83, 0.15)"] * len(row)
        return [""] * len(row)

    styled_rat = df_rat.style.apply(_highlight_best_rational, axis=1)

    st.dataframe(
        styled_rat,
        use_container_width=True,
        hide_index=True,
        column_config={
            "k": st.column_config.NumberColumn(width="small"),
            "Spacing (1/k)": st.column_config.TextColumn(width="medium"),
            "Matches": st.column_config.NumberColumn(width="small"),
            "Expected (Random)": st.column_config.TextColumn(width="small"),
            "Excess": st.column_config.TextColumn(width="small"),
            "Avg Error": st.column_config.TextColumn(width="small"),
        },
    )

# ---------------------------------------------------------------------------
# Tab 2: Irrational Spacings
# ---------------------------------------------------------------------------
with tab2:
    st.subheader("Irrational Spacings")
    st.markdown(
        f"Testing 15 irrational spacings from the legacy analysis "
        f"against {n_particles} particles."
    )

    irr_rows = []
    for name, spacing, avg, matches, details in irrational_results:
        excess = matches - random_expected
        irr_rows.append({
            "Spacing Name": name,
            "Spacing Value": f"{spacing:.6f}",
            "Matches": matches,
            "Expected (Random)": f"{random_expected:.1f}",
            "Excess": f"{excess:+.1f}",
            "Avg Error": f"{avg:.4f}",
        })

    df_irr = pd.DataFrame(irr_rows)

    def _highlight_best_irrational(row):
        if int(row["Matches"]) == max(int(r["Matches"]) for r in irr_rows):
            return ["background-color: rgba(0, 200, 83, 0.15)"] * len(row)
        return [""] * len(row)

    styled_irr = df_irr.style.apply(_highlight_best_irrational, axis=1)

    st.dataframe(
        styled_irr,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Spacing Name": st.column_config.TextColumn(width="large"),
            "Spacing Value": st.column_config.TextColumn(width="medium"),
            "Matches": st.column_config.NumberColumn(width="small"),
            "Expected (Random)": st.column_config.TextColumn(width="small"),
            "Excess": st.column_config.TextColumn(width="small"),
            "Avg Error": st.column_config.TextColumn(width="small"),
        },
    )

# ---------------------------------------------------------------------------
# Tab 3: Continuous Optimization
# ---------------------------------------------------------------------------
with tab3:
    st.subheader("Continuous Optimization")
    st.markdown(
        "Brute force scan of spacings from 0.01 to 2.0 in steps of 0.001. "
        "The global best and top local minima are shown."
    )

    # Global best
    st.markdown(
        f"""
        <div class="best-fit-card">
        <b>Global Best Spacing:</b> {best_spacing:.4f}<br>
        <b>RMS Score:</b> {best_score:.6f}
        </div>
        """,
        unsafe_allow_html=True,
    )

    # Known constants for identification
    phi = (1 + math.sqrt(5)) / 2
    _known_constants = {
        "1/2": 0.5,
        "1/3": 1 / 3,
        "1/4": 0.25,
        "1/5": 0.2,
        "1/6": 1 / 6,
        "1/φ": 1 / phi,
        "1/π": 1 / math.pi,
        "1/e": 1 / math.e,
        "ln2": math.log(2),
        "1/φ²": 1 / phi**2,
        "φ/π": phi / math.pi,
        "2/π": 2 / math.pi,
        "1/√(2)": 1 / math.sqrt(2),
        "1/√(3)": 1 / math.sqrt(3),
        "φ - 1": phi - 1,
        "ln(φ)": math.log(phi),
        "π/6": math.pi / 6,
        "1/(2·φ)": 1 / (2 * phi),
        "1/(3·φ)": 1 / (3 * phi),
        "φ/4": phi / 4,
        "1": 1.0,
        "1/7": 1 / 7,
        "1/8": 0.125,
    }

    def _identify_nearest_constant(spacing_val):
        """Find the nearest known constant to a spacing value."""
        best_name = "?"
        best_diff = float("inf")
        for name, val in _known_constants.items():
            diff = abs(spacing_val - val)
            if diff < best_diff:
                best_diff = diff
                best_name = name
        return best_name, best_diff

    # Top 10 local minima
    st.markdown("**Top 10 Local Minima**")
    top_n = min(10, len(local_mins))
    min_rows = []
    for spacing_val, rms in local_mins[:top_n]:
        const_name, const_diff = _identify_nearest_constant(spacing_val)
        min_rows.append({
            "Spacing": f"{spacing_val:.4f}",
            "RMS Score": f"{rms:.6f}",
            "Nearest Constant": const_name,
            "Distance": f"{const_diff:.4f}",
        })

    df_mins = pd.DataFrame(min_rows)

    def _highlight_global_best(row):
        if row["Spacing"] == f"{best_spacing:.4f}":
            return ["background-color: rgba(245, 158, 11, 0.2)"] * len(row)
        return [""] * len(row)

    styled_mins = df_mins.style.apply(_highlight_global_best, axis=1)

    st.dataframe(
        styled_mins,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Spacing": st.column_config.TextColumn(width="small"),
            "RMS Score": st.column_config.TextColumn(width="medium"),
            "Nearest Constant": st.column_config.TextColumn(width="medium"),
            "Distance": st.column_config.TextColumn(width="small"),
        },
    )

    # Spacing score chart
    st.markdown("**Score Distribution**")
    chart_data = []
    for spacing_val, rms in local_mins[:top_n]:
        const_name, _ = _identify_nearest_constant(spacing_val)
        chart_data.append({
            "spacing": spacing_val,
            "score": rms,
            "label": f"{spacing_val:.3f} ({const_name})",
        })

    fig = spacing_score_chart(chart_data)
    st.plotly_chart(fig, use_container_width=True)

st.divider()

# ---------------------------------------------------------------------------
# Best Fit Detail (below tabs)
# ---------------------------------------------------------------------------
st.subheader("Best Fit Detail")
st.markdown(
    f"The winning spacing **{best_spacing:.4f}** applied to all particles."
)

best_details = get_best_fit_details(best_spacing, rung_dict)

detail_rows = []
for name, n, nearest_rung, k, delta, frac, match in best_details:
    detail_rows.append({
        "Particle": name,
        "Rung n": f"{n:.4f}",
        "Nearest Rung": f"{nearest_rung:.4f}",
        "k": k,
        "Delta": f"{delta:.4f}",
        "Match": "Yes" if match else "No",
    })

df_details = pd.DataFrame(detail_rows)


def _highlight_match_detail(row):
    if row["Match"] == "Yes":
        return ["background-color: rgba(0, 200, 83, 0.15)"] * len(row)
    return ["background-color: rgba(213, 0, 0, 0.08)"] * len(row)


styled_details = df_details.style.apply(_highlight_match_detail, axis=1)

st.dataframe(
    styled_details,
    use_container_width=True,
    hide_index=True,
    column_config={
        "Particle": st.column_config.TextColumn(width="medium"),
        "Rung n": st.column_config.TextColumn(width="small"),
        "Nearest Rung": st.column_config.TextColumn(width="small"),
        "k": st.column_config.NumberColumn(width="small"),
        "Delta": st.column_config.TextColumn(width="small"),
        "Match": st.column_config.TextColumn(width="small"),
    },
)

n_detail_matches = sum(1 for d in best_details if d[6])
st.markdown(
    f"**{n_detail_matches} of {len(best_details)}** particles match "
    f"with spacing {best_spacing:.4f} (tolerance < 15% of spacing)."
)

st.divider()
st.caption("Rung Spacing | Alpha Ladder Research Dashboard")
