"""
Dilaton Lab -- Page 8

Dilaton interpretation of the gravity rung (alpha^21).
Displays 20+1 decomposition, alternative splits, Brans-Dicke analysis,
and 6D reconciliation.
"""

import streamlit as st
from decimal import Decimal, getcontext

getcontext().prec = 50

st.set_page_config(page_title="Dilaton Lab | Alpha Ladder", layout="wide")

# ---------------------------------------------------------------------------
# Custom CSS
# ---------------------------------------------------------------------------
st.markdown(
    """
    <style>
    .decomp-card {
        background-color: #1a1d23;
        border: 1px solid #2e3440;
        border-radius: 8px;
        padding: 1rem;
        margin: 0.5rem 0;
        font-family: 'Fira Mono', Consolas, monospace;
        font-size: 1.1rem;
    }
    .reconcile-box {
        background-color: #1a1d23;
        border-left: 4px solid #60a5fa;
        padding: 1.2rem;
        border-radius: 0 8px 8px 0;
        font-family: 'Fira Mono', monospace;
        font-size: 1.1rem;
    }
    .highlight-row {
        background-color: #1e293b;
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
from app.components.formatting import fmt_decimal, fmt_percent  # noqa: E402

constants = render_sidebar()

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Dilaton Lab")
st.markdown(
    "Dilaton interpretation of the gravity rung: decomposing α²¹ "
    "into geometric and scalar sectors."
)
st.divider()

# ---------------------------------------------------------------------------
# Try to load core modules
# ---------------------------------------------------------------------------
_core_available = False
try:
    from alpha_ladder_core.dilaton import (  # noqa: E402
        decompose_21,
        compute_bd_parameter,
        reconcile_6d,
    )
    _core_available = True
except ImportError:
    pass

# ---------------------------------------------------------------------------
# Section 1: 20+1 Decomposition
# ---------------------------------------------------------------------------
st.subheader("20+1 Decomposition")
st.markdown(
    "The exponent 21 splits as **20 + 1**: the geometric sector "
    "(Riemann tensor, 20 independent components in 4D) plus one "
    "scalar degree of freedom (the dilaton)."
)

if _core_available and constants is not None:
    decomp = decompose_21(constants)
else:
    # Fallback values
    alpha_val = 0.0072973525693
    phi = (1 + 5 ** 0.5) / 2
    decomp = {
        "alpha_20": alpha_val ** 20,
        "alpha_1": alpha_val,
        "bridge": phi ** 2 / 2,
        "alpha_21": alpha_val ** 21,
        "alternative_splits": [],
    }

col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric(
        label="α²⁰ (Geometric Sector)",
        value=fmt_decimal(decomp["alpha_20"], sig_figs=6),
    )

with col2:
    st.metric(
        label="α¹ (Dilaton Sector)",
        value=fmt_decimal(decomp["alpha_1"], sig_figs=10),
    )

with col3:
    st.metric(
        label="Bridge Coefficient (φ²/2)",
        value=fmt_decimal(decomp["bridge"], sig_figs=8),
    )

with col4:
    st.metric(
        label="α²¹ (Full Product)",
        value=fmt_decimal(decomp["alpha_21"], sig_figs=6),
    )

st.divider()

# ---------------------------------------------------------------------------
# Section 2: Alternative Splits
# ---------------------------------------------------------------------------
st.subheader("Alternative Decompositions of 21")
st.markdown(
    "Any factorization a + b = 21 is mathematically valid, but only "
    "**20 + 1** has a clear physical interpretation: Riemann tensor "
    "components plus one scalar."
)

if decomp.get("alternative_splits"):
    import pandas as pd

    rows = []
    for sp in decomp["alternative_splits"]:
        rows.append({
            "Split": f"{sp['a']} + {sp['b']}",
            "αᵃ": fmt_decimal(sp["alpha_a"], sig_figs=6),
            "αᵇ": fmt_decimal(sp["alpha_b"], sig_figs=6),
            "Geometric Meaning": sp["meaning"],
        })

    df = pd.DataFrame(rows)
    st.dataframe(
        df,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Split": st.column_config.TextColumn(width="small"),
            "αᵃ": st.column_config.TextColumn(width="medium"),
            "αᵇ": st.column_config.TextColumn(width="medium"),
            "Geometric Meaning": st.column_config.TextColumn(width="large"),
        },
    )

    st.info(
        "The 20+1 split is highlighted as the physically motivated decomposition. "
        "The Riemann curvature tensor in 4D has exactly 20 independent components, "
        "and the remaining scalar couples as a dilaton field."
    )
else:
    st.warning(
        "Alternative split data not available. "
        "Ensure the alpha_ladder_core.dilaton module is installed."
    )

st.divider()

# ---------------------------------------------------------------------------
# Section 3: Brans-Dicke Analysis
# ---------------------------------------------------------------------------
st.subheader("Brans-Dicke Analysis")
st.markdown(
    "If the bridge coefficient φ²/2 arises from a scalar-tensor theory, "
    "it implies a specific Brans-Dicke parameter omega."
)

if _core_available and constants is not None:
    bd = compute_bd_parameter(constants)
else:
    # Fallback values
    phi_f = (1 + 5 ** 0.5) / 2
    phi2 = phi_f ** 2
    bridge = phi2 / 2.0
    numerator = 3 * phi2 / 2 - 4
    denominator = 2 - phi2
    omega = numerator / denominator if denominator != 0 else None

    import math as _math
    hbar_f = 1.054571817e-34
    c_f = 2.99792458e8
    AU = 1.496e11
    cassini_limit = 2.3e-5
    denom_fifth = 2.0 * omega + 3.0 if omega is not None else 3.0
    alpha_fifth = 2.0 / denom_fifth if denom_fifth != 0 else 0.0
    if alpha_fifth > cassini_limit:
        suppression = _math.log(alpha_fifth / cassini_limit)
        lambda_max = AU / suppression
    else:
        suppression = 0.0
        lambda_max = AU
    m_d_min = hbar_f / (c_f * lambda_max)
    m_d_min_eV = m_d_min * c_f ** 2 / 1.602176634e-19

    bd = {
        "omega": omega,
        "phi_squared_over_2": bridge,
        "cassini_bound": 40000,
        "omega_excluded_massless": omega is not None and abs(omega) < 40000,
        "omega_negative": omega is not None and omega < 0,
        "dilaton_mass_min_kg": m_d_min,
        "dilaton_mass_min_eV": m_d_min_eV,
        "cassini_lambda_max_m": lambda_max,
        "cassini_suppression_factor": suppression,
        "dark_scale_eV": 0.0072973525693 ** 5 * 1.22089e28,
    }

col_bd1, col_bd2 = st.columns(2)

with col_bd1:
    st.markdown("**Brans-Dicke Parameter**")
    if bd["omega"] is not None:
        st.metric(label="omega (BD)", value=f"{bd['omega']:.4f}")
    else:
        st.metric(label="omega (BD)", value="Undefined")

    st.metric(
        label="φ²/2 (bare ratio)",
        value=f"{bd['phi_squared_over_2']:.6f}",
    )

with col_bd2:
    st.markdown("**Cassini Bound Comparison**")
    st.metric(label="Cassini lower bound on |omega|", value=f"{bd['cassini_bound']:,}")

    if bd["omega"] is not None:
        st.metric(
            label="|omega| from ladder",
            value=f"{abs(bd['omega']):.4f}",
        )

# Key conclusions
if bd.get("omega_negative"):
    st.warning(
        "omega is **negative**, placing this theory outside standard "
        "Brans-Dicke territory. Negative omega corresponds to a ghost-like "
        "scalar (wrong-sign kinetic term) in the Jordan frame."
    )

if bd.get("omega_excluded_massless"):
    st.error(
        "For a **massless** dilaton, the Cassini bound requires |omega| > 40,000. "
        "The ladder value is far below this threshold. A massless dilaton with "
        "this coupling is **excluded** by Solar System tests."
    )

st.info(
    "However, if the dilaton is **massive** (acquires a mass through some mechanism), "
    "the Cassini bound does not apply at distances shorter than the Compton wavelength. "
    "This opens the door for a massive dilaton interpretation."
)

with st.expander("Dilaton Mass Requirements"):
    st.markdown(f"- Minimum dilaton mass to evade Cassini: "
                f"**{fmt_decimal(bd['dilaton_mass_min_kg'], sig_figs=4)} kg**")
    st.markdown(f"- Minimum dilaton mass in eV: "
                f"**{fmt_decimal(bd['dilaton_mass_min_eV'], sig_figs=4)} eV**")
    if bd.get("dark_scale_eV") is not None:
        st.markdown(f"- Dark/dilaton mass scale from ladder (α⁵ · m_Pl): "
                    f"**{fmt_decimal(bd['dark_scale_eV'], sig_figs=4)} eV**")
        st.markdown(
            "The ladder's own mass scale is far above the Cassini minimum, "
            "suggesting self-consistency: the dilaton can be massive enough "
            "to evade Solar System bounds."
        )

st.divider()

# ---------------------------------------------------------------------------
# Section 4: 6D Reconciliation
# ---------------------------------------------------------------------------
st.subheader("6D Reconciliation")
st.markdown(
    "Two complementary ways to understand why the exponent is 21: "
    "bottom-up from 4D gravity, or top-down from 6D geometry."
)

if _core_available and constants is not None:
    recon = reconcile_6d(constants)
else:
    recon = {
        "g_uv_4d": 10,
        "g_ab_2d": 3,
        "g_ua_mixed": 8,
        "total_6d_metric": 21,
        "riemann_4d": 20,
        "scalar_dilaton": 1,
        "riemann_plus_scalar": 21,
        "views": {
            "bottom_up": "4D gravity needs 20+1 to couple to the SM",
            "top_down": "6D gravity has 21 metric components naturally",
        },
    }

col_6d_left, col_6d_right = st.columns(2)

with col_6d_left:
    st.markdown("**Bottom-Up: 4D Riemann + Scalar**")
    st.markdown(
        """
        <div class="reconcile-box">
        <b>21 = 20 + 1</b><br><br>
        Riemann tensor (4D): <b>20</b> independent components<br>
        Scalar (dilaton): <b>1</b> degree of freedom<br>
        <br>
        Total: <b>"""
        + str(recon["riemann_plus_scalar"])
        + """</b>
        </div>
        """,
        unsafe_allow_html=True,
    )
    st.caption(recon["views"]["bottom_up"])

with col_6d_right:
    st.markdown("**Top-Down: 6D Metric Decomposition**")
    st.markdown(
        """
        <div class="reconcile-box">
        <b>21 = 10 + 3 + 8</b><br><br>
        g_uv (4D metric): <b>"""
        + str(recon["g_uv_4d"])
        + """</b> components (symmetric 4x4)<br>
        g_ab (2D internal): <b>"""
        + str(recon["g_ab_2d"])
        + """</b> components (symmetric 2x2)<br>
        g_ua (mixed): <b>"""
        + str(recon["g_ua_mixed"])
        + """</b> components (4 x 2)<br>
        <br>
        Total: <b>"""
        + str(recon["total_6d_metric"])
        + """</b>
        </div>
        """,
        unsafe_allow_html=True,
    )
    st.caption(recon["views"]["top_down"])

st.markdown("---")
st.markdown(
    "Both decompositions yield **21** -- the same exponent that appears in the "
    "Alpha Ladder bridge. This is either a deep structural coincidence or "
    "evidence that gravity's coupling strength is set by the geometry of "
    "a higher-dimensional theory compactified to 4D."
)

st.divider()
st.caption("Dilaton Lab | Alpha Ladder Research Dashboard")
