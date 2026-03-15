"""
Particle Harmonics -- Page 6

Particles plotted on the alpha-rung number line.
Shows mass spectrum table, number line visualization, and summary statistics.
"""

import streamlit as st

st.set_page_config(page_title="Particle Harmonics | Alpha Ladder", layout="wide")

# ---------------------------------------------------------------------------
# Custom CSS
# ---------------------------------------------------------------------------
st.markdown(
    """
    <style>
    .harmonics-summary {
        background-color: #1a1d23;
        border: 1px solid #2e3440;
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
from app.components.formatting import fmt_decimal  # noqa: E402
from app.components.charts import particle_number_line  # noqa: E402

constants = render_sidebar()

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Particle Harmonics")
st.markdown("Particles plotted on the alpha-rung number line.")
st.divider()

# ---------------------------------------------------------------------------
# Graceful imports
# ---------------------------------------------------------------------------
_core_available = False
try:
    from alpha_ladder_core.particle_harmonics import compute_harmonics  # noqa: E402
    _core_available = True
except ImportError:
    st.error(
        "Could not import alpha_ladder_core.particle_harmonics. "
        "Ensure the core library is installed."
    )
    st.stop()

import pandas as pd  # noqa: E402

# ---------------------------------------------------------------------------
# Build constants argument
# ---------------------------------------------------------------------------
if constants is not None and hasattr(constants, "alpha"):
    _constants_arg = constants
else:
    _constants_arg = {"alpha": 0.0072973525693}

# ---------------------------------------------------------------------------
# Compute harmonics
# ---------------------------------------------------------------------------
harmonics = compute_harmonics(_constants_arg)

# ---------------------------------------------------------------------------
# Section 1: Mass Spectrum Table
# ---------------------------------------------------------------------------
st.subheader("Mass Spectrum Table")
st.markdown(
    "All particles with mass, mass ratio to electron, α rung number (n), "
    "nearest half-integer, closeness, and match status."
)

rows = []
for h in harmonics:
    rows.append({
        "Particle": h["name"],
        "Mass (MeV)": f"{h['mass']:.4f}",
        "Ratio (m/mₑ)": f"{h['ratio']:.6f}",
        "Rung n": f"{h['rung_n']:.4f}",
        "Nearest Half-Int": f"{h['nearest_half']:.1f}",
        "Closeness": f"{h['closeness']:.4f}",
        "Match": "Yes" if h["is_match"] else "No",
    })

df_harmonics = pd.DataFrame(rows)

# Apply conditional formatting via Pandas Styler
def _highlight_match(row):
    """Return background styles based on match status."""
    if row["Match"] == "Yes":
        return ["background-color: rgba(0, 200, 83, 0.15)"] * len(row)
    else:
        return [""] * len(row)


styled_df = df_harmonics.style.apply(_highlight_match, axis=1)

st.dataframe(
    styled_df,
    use_container_width=True,
    hide_index=True,
    column_config={
        "Particle": st.column_config.TextColumn(width="medium"),
        "Mass (MeV)": st.column_config.TextColumn(width="small"),
        "Ratio (m/mₑ)": st.column_config.TextColumn(width="medium"),
        "Rung n": st.column_config.TextColumn(width="small"),
        "Nearest Half-Int": st.column_config.TextColumn(width="small"),
        "Closeness": st.column_config.TextColumn(width="small"),
        "Match": st.column_config.TextColumn(width="small"),
    },
)

st.divider()

# ---------------------------------------------------------------------------
# Section 2: Number Line Plot
# ---------------------------------------------------------------------------
st.subheader("Number Line Plot")
st.markdown(
    "Particles as dots on a horizontal axis where x = rung number. "
    "Green diamonds indicate matches (closeness < 0.05), red diamonds indicate non-matches."
)

# Build data for the chart function
chart_data = []
for h in harmonics:
    chart_data.append({
        "name": h["name"],
        "rung": h["rung_n"],
        "match": h["is_match"],
    })

fig = particle_number_line(chart_data)
st.plotly_chart(fig, use_container_width=True)

st.divider()

# ---------------------------------------------------------------------------
# Section 3: Summary
# ---------------------------------------------------------------------------
st.subheader("Summary")

n_particles = len(harmonics)
n_matches = sum(1 for h in harmonics if h["is_match"])
n_non_matches = n_particles - n_matches

# Random expectation: with tolerance 0.05 on each side of a half-integer,
# we cover 0.10 out of every 0.50 interval, so probability ~20% per particle.
random_expected = n_particles * 0.10 / 0.50

col1, col2, col3 = st.columns(3)

with col1:
    st.metric("Total Particles", n_particles)

with col2:
    st.metric(
        "Matches (near rungs)",
        n_matches,
        delta=f"{n_matches / n_particles * 100:.0f}%" if n_particles > 0 else "0%",
    )

with col3:
    st.metric(
        "Random Expectation",
        f"{random_expected:.1f}",
        delta=f"{n_matches - random_expected:+.1f} excess",
        delta_color="normal" if n_matches > random_expected else "inverse",
    )

st.markdown(
    f"""
    <div class="harmonics-summary">
    <b>Result:</b> {n_matches} of {n_particles} particles land within 0.05 of a
    half-integer rung, compared to a random expectation of {random_expected:.1f}.
    {"This represents a statistically significant excess." if n_matches > random_expected + 1 else "The observed count is close to random expectation."}
    </div>
    """,
    unsafe_allow_html=True,
)

st.divider()
st.caption("Particle Harmonics | Alpha Ladder Research Dashboard")
