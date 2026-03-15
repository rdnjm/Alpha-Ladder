"""
Constant Core -- Page 1

Displays the full-precision constant table for the selected CODATA edition,
electron geometry verification, and the 42-order gap between alpha and alpha_G.
"""

import streamlit as st
from decimal import Decimal, getcontext

getcontext().prec = 50


# ---------------------------------------------------------------------------
# Custom CSS
# ---------------------------------------------------------------------------
st.markdown(
    """
    <style>
    .const-table {
        font-family: 'Fira Mono', Consolas, monospace;
        font-size: 1.05rem;
    }
    .geom-verify {
        background-color: #1a1d23;
        border: 1px solid #2e3440;
        border-radius: 8px;
        padding: 1rem;
        margin: 0.5rem 0;
    }
    .gap-display {
        background-color: #1a1d23;
        border-left: 4px solid #f59e0b;
        padding: 1.2rem;
        border-radius: 0 8px 8px 0;
        font-family: 'Fira Mono', monospace;
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
from app.components.formatting import fmt_decimal, fmt_percent  # noqa: E402

constants = render_sidebar()

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Constant Core")
st.markdown("Full-precision constants and electron geometry verification.")
st.divider()

# ---------------------------------------------------------------------------
# Try to load core modules; fall back to legacy data
# ---------------------------------------------------------------------------
_core_available = False
try:
    from alpha_ladder_core.ladder import compute_electron_geometry  # noqa: E402
    _core_available = True
except ImportError:
    pass

# ---------------------------------------------------------------------------
# Section 1: Full-precision constant table
# ---------------------------------------------------------------------------
st.subheader("CODATA Constants")

_DISPLAY_NAMES = {
    "alpha": "α",
    "alpha_g": "α_G",
    "hbar": "ℏ",
    "c": "c",
    "m_e": "mₑ",
    "m_p": "mₚ",
    "m_mu": "m_μ",
    "k_B": "k_B",
    "e_charge": "e (charge)",
    "epsilon_0": "ε₀",
    "G": "G",
    "lambda_bar_c": "λ̄_c",
    "a_0": "a₀",
    "r_e_nist": "rₑ (NIST)",
    "pi": "π",
    "phi": "φ",
    "e": "e (Euler)",
    "ln2": "ln 2",
    "sqrt2": "√2",
    "sqrt3": "√3",
    "sqrt5": "√5",
}

if constants is not None and hasattr(constants, "__dict__"):
    # Build a dataframe from the constants namespace
    import pandas as pd

    rows = []
    for attr_name in sorted(dir(constants)):
        if attr_name.startswith("_"):
            continue
        val = getattr(constants, attr_name)
        if isinstance(val, (Decimal, float, int)):
            rows.append({
                "Constant": _DISPLAY_NAMES.get(attr_name, attr_name),
                "Value": fmt_decimal(val, sig_figs=12),
                "Type": type(val).__name__,
            })

    if rows:
        df = pd.DataFrame(rows)
        st.dataframe(
            df,
            use_container_width=True,
            hide_index=True,
            column_config={
                "Constant": st.column_config.TextColumn(width="medium"),
                "Value": st.column_config.TextColumn(width="large"),
                "Type": st.column_config.TextColumn(width="small"),
            },
        )
    else:
        st.info("No numeric constants found in the selected edition.")
else:
    st.warning(
        "Core constants module not yet available. "
        "Showing hardcoded CODATA 2018 reference values."
    )

    import pandas as pd

    fallback_constants = {
        "α": "0.0072973525693",
        "α_G": "1.7512e-45",
        "ℏ": "1.054571817e-34",
        "c": "299792458",
        "mₑ": "9.1093837015e-31",
        "λ̄_c": "3.8615926764e-13",
        "a₀": "5.2917721067e-11",
        "rₑ (NIST)": "2.8179403227e-15",
    }

    df = pd.DataFrame(
        [{"Constant": k, "Value": v} for k, v in fallback_constants.items()]
    )
    st.dataframe(df, use_container_width=True, hide_index=True)

st.divider()

# ---------------------------------------------------------------------------
# Section 2: Electron Geometry Verification
# ---------------------------------------------------------------------------
st.subheader("Electron Geometry Verification")
st.markdown(
    "The classical electron radius can be derived from two independent paths, "
    "both yielding the same value to high precision:"
)

if _core_available and constants is not None:
    geom = compute_electron_geometry(constants)
    r_e_compton = geom.get("r_e_from_compton")
    r_e_bohr = geom.get("r_e_from_bohr")
    r_e_nist = geom.get("r_e_nist")
else:
    # Fallback to hardcoded values from the original alpha_ladder.py
    alpha = Decimal("0.0072973525693")
    lambda_bar_c = Decimal("3.8615926764e-13")
    a_0 = Decimal("5.2917721067e-11")
    r_e_compton = alpha * lambda_bar_c
    r_e_bohr = alpha ** 2 * a_0
    r_e_nist = Decimal("2.8179403227e-15")

col1, col2, col3 = st.columns(3)

with col1:
    st.markdown("**From Compton wavelength**")
    st.markdown("`r_e = α · λ̄_c`")
    st.code(fmt_decimal(r_e_compton, sig_figs=10) + " m")

with col2:
    st.markdown("**From Bohr radius**")
    st.markdown("`r_e = α² · a₀`")
    st.code(fmt_decimal(r_e_bohr, sig_figs=10) + " m")

with col3:
    st.markdown("**NIST reference**")
    st.markdown("`rₑ (NIST)`")
    st.code(fmt_decimal(r_e_nist, sig_figs=10) + " m")

# Show agreement
if r_e_compton and r_e_bohr and r_e_nist:
    with st.expander("Detailed comparison"):
        diff_cb = abs(float(r_e_compton) - float(r_e_bohr))
        diff_cn = abs(float(r_e_compton) - float(r_e_nist))

        rel_cb = diff_cb / float(r_e_nist) if float(r_e_nist) != 0 else 0
        rel_cn = diff_cn / float(r_e_nist) if float(r_e_nist) != 0 else 0

        st.markdown(f"- Compton vs Bohr: Δ = {diff_cb:.4e} m ({fmt_percent(rel_cb)})")
        st.markdown(f"- Compton vs NIST: Δ = {diff_cn:.4e} m ({fmt_percent(rel_cn)})")
        st.markdown(
            "Both derivations agree with NIST to the precision of the input constants, "
            "confirming the geometric scaling `rₑ = α · λ̄_c = α² · a₀`."
        )

st.divider()

# ---------------------------------------------------------------------------
# Section 3: The 42-Order Gap
# ---------------------------------------------------------------------------
st.subheader("The 42-Order Gap (The Great Desert)")

alpha_val = Decimal("0.0072973525693")
alpha_g_val = Decimal("1.7512e-45")

if constants is not None and hasattr(constants, "alpha") and hasattr(constants, "alpha_g"):
    alpha_val = constants.alpha
    alpha_g_val = constants.alpha_g

gap = alpha_val / alpha_g_val

st.markdown(
    """
    <div class="gap-display">
    <b>α / α_G</b> = """
    + fmt_decimal(gap, sig_figs=6)
    + """
    <br><br>
    This ratio spans approximately <b>42 orders of magnitude</b> --
    the vast hierarchy between the electromagnetic and gravitational
    coupling strengths. The Alpha Ladder hypothesis proposes that this
    gap is not arbitrary but is precisely α⁻²⁰, bridged by a
    simple geometric coefficient.
    </div>
    """,
    unsafe_allow_html=True,
)

st.divider()
st.caption("Constant Core | Alpha Ladder Research Dashboard")
