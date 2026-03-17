"""
Alpha Units -- Page 10

SI to alpha-unit converter and notable conversions table.
"""

import streamlit as st
import math
from decimal import Decimal, getcontext

getcontext().prec = 50



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
# Header
# ---------------------------------------------------------------------------
st.title("Alpha Units")
st.markdown(
    "Natural unit system built on the fine-structure constant α, "
    "the electron mass mₑ, and the speed of light c."
)
st.divider()

# ---------------------------------------------------------------------------
# Try to load core modules
# ---------------------------------------------------------------------------
_core_available = False
try:
    from alpha_ladder_core.alpha_units import (  # noqa: E402
        get_alpha_units,
        si_to_alpha,
        alpha_to_si,
        express_in_alpha_powers,
    )
    _core_available = True
except ImportError:
    pass

# ---------------------------------------------------------------------------
# Fallback: build alpha units from hardcoded constants
# ---------------------------------------------------------------------------
_FALLBACK_UNITS = None

if not _core_available or constants is None:
    _alpha = 0.0072973525693
    _a_0 = 5.2917721067e-11
    _c = 2.99792458e8
    _m_e = 9.1093837015e-31
    _e_charge = 1.602176634e-19
    _k_B = 1.380649e-23

    _r_e = _alpha ** 2 * _a_0
    _t_unit = _r_e / _c
    _energy = _m_e * _c ** 2
    _force = _energy / _r_e

    _FALLBACK_UNITS = {
        "length": _r_e,
        "time": _t_unit,
        "mass": _m_e,
        "energy": _energy,
        "charge": _e_charge,
        "force": _force,
        "temperature": _energy / _k_B,
    }


def _get_units():
    """Get alpha unit definitions, using core or fallback."""
    if _core_available and constants is not None:
        return get_alpha_units(constants)
    return _FALLBACK_UNITS


def _si_to_alpha_safe(value, unit_type):
    """Convert SI to alpha units, using core or fallback."""
    if _core_available and constants is not None:
        return si_to_alpha(value, unit_type, constants)
    units = _get_units()
    if unit_type not in units:
        raise ValueError(f"Unknown unit_type {unit_type!r}")
    return value / units[unit_type]


def _alpha_to_si_safe(value, unit_type):
    """Convert alpha units to SI, using core or fallback."""
    if _core_available and constants is not None:
        return alpha_to_si(value, unit_type, constants)
    units = _get_units()
    if unit_type not in units:
        raise ValueError(f"Unknown unit_type {unit_type!r}")
    return value * units[unit_type]


def _express_powers_safe(si_value, unit_type):
    """Express SI value as alpha powers, using core or fallback."""
    if _core_available and constants is not None:
        return express_in_alpha_powers(si_value, unit_type, constants)
    alpha = 0.0072973525693
    alpha_val = _si_to_alpha_safe(si_value, unit_type)
    if alpha_val <= 0:
        return None
    exact_power = math.log(alpha_val) / math.log(alpha)
    nearest_power = round(exact_power)
    approx = alpha ** nearest_power
    residual_pct = abs(approx - alpha_val) / abs(alpha_val) * 100.0
    return {
        "alpha_unit_value": alpha_val,
        "nearest_power": nearest_power,
        "exact_power": exact_power,
        "residual_pct": residual_pct,
    }


# ---------------------------------------------------------------------------
# Section 1: Alpha Unit System
# ---------------------------------------------------------------------------
st.subheader("Alpha Unit System")
st.markdown(
    "Base units defined in terms of α and electron properties. "
    "Each unit is a combination of the classical electron radius rₑ, "
    "electron mass mₑ, and the speed of light c."
)

units = _get_units()

if units is not None:
    import pandas as pd

    _SI_LABELS = {
        "length": "metres (m)",
        "time": "seconds (s)",
        "mass": "kilograms (kg)",
        "energy": "Joules (J)",
        "charge": "Coulombs (C)",
        "force": "Newtons (N)",
        "temperature": "Kelvin (K)",
    }

    _DEFINITIONS = {
        "length": "rₑ = α² · a₀",
        "time": "rₑ / c",
        "mass": "mₑ",
        "energy": "mₑ · c²",
        "charge": "e (elementary charge)",
        "force": "mₑ · c² / rₑ",
        "temperature": "mₑ · c² / k_B",
    }

    rows = []
    for utype in ["length", "time", "mass", "energy", "charge", "force", "temperature"]:
        if utype in units:
            rows.append({
                "Unit Type": utype.capitalize(),
                "Definition": _DEFINITIONS.get(utype, ""),
                "SI Equivalent": fmt_decimal(units[utype], sig_figs=8),
                "SI Unit": _SI_LABELS.get(utype, ""),
            })

    df_units = pd.DataFrame(rows)
    st.dataframe(
        df_units,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Unit Type": st.column_config.TextColumn(width="small"),
            "Definition": st.column_config.TextColumn(width="medium"),
            "SI Equivalent": st.column_config.TextColumn(width="medium"),
            "SI Unit": st.column_config.TextColumn(width="small"),
        },
    )
else:
    st.warning(
        "Alpha unit system not available. "
        "Ensure the alpha_ladder_core.alpha_units module is installed."
    )

st.divider()

# ---------------------------------------------------------------------------
# Section 2: SI to Alpha-Units Converter
# ---------------------------------------------------------------------------
st.subheader("SI to Alpha-Units Converter")

available_types = list(units.keys()) if units else [
    "length", "time", "mass", "energy", "charge", "force", "temperature"
]

col_conv1, col_conv2 = st.columns([1, 2])

with col_conv1:
    unit_type_fwd = st.selectbox(
        "Unit type",
        options=available_types,
        format_func=lambda x: x.capitalize(),
        key="fwd_unit_type",
    )

    si_value = st.number_input(
        "SI value",
        value=1.0e-15,
        format="%e",
        key="fwd_si_value",
        help="Enter the value in SI units. Use scientific notation (e.g. 1.6e-19).",
    )

with col_conv2:
    if si_value != 0 and units is not None:
        try:
            alpha_result = _si_to_alpha_safe(si_value, unit_type_fwd)

            st.markdown(
                f"""
                <div class="converter-result">
                <b>{fmt_decimal(si_value, sig_figs=6)} {_SI_LABELS.get(unit_type_fwd, '')}</b><br>
                = <b>{fmt_decimal(alpha_result, sig_figs=8)}</b> alpha-units ({unit_type_fwd})
                </div>
                """,
                unsafe_allow_html=True,
            )

            # Alpha-power decomposition
            if si_value > 0:
                powers = _express_powers_safe(si_value, unit_type_fwd)
                if powers is not None:
                    st.markdown("**Alpha-Power Decomposition**")
                    col_p1, col_p2, col_p3 = st.columns(3)
                    with col_p1:
                        st.metric(
                            label="Nearest αⁿ",
                            value=f"α{powers['nearest_power']}",
                        )
                    with col_p2:
                        st.metric(
                            label="Exact power",
                            value=f"{powers['exact_power']:.4f}",
                        )
                    with col_p3:
                        st.metric(
                            label="Residual",
                            value=f"{powers['residual_pct']:.2f}%",
                        )

        except (ValueError, ZeroDivisionError) as exc:
            st.error(f"Conversion error: {exc}")
    else:
        st.info("Enter a non-zero SI value to see the conversion.")

st.divider()

# ---------------------------------------------------------------------------
# Section 3: Reverse Converter (Alpha-Units to SI)
# ---------------------------------------------------------------------------
st.subheader("Alpha-Units to SI Converter")

col_rev1, col_rev2 = st.columns([1, 2])

with col_rev1:
    unit_type_rev = st.selectbox(
        "Unit type",
        options=available_types,
        format_func=lambda x: x.capitalize(),
        key="rev_unit_type",
    )

    alpha_value = st.number_input(
        "Alpha-unit value",
        value=1.0,
        format="%e",
        key="rev_alpha_value",
        help="Enter the value in alpha-units.",
    )

with col_rev2:
    if alpha_value != 0 and units is not None:
        try:
            si_result = _alpha_to_si_safe(alpha_value, unit_type_rev)

            _SI_LABELS_REV = _SI_LABELS if '_SI_LABELS' in dir() else {
                "length": "metres (m)", "time": "seconds (s)",
                "mass": "kilograms (kg)", "energy": "Joules (J)",
                "charge": "Coulombs (C)", "force": "Newtons (N)",
                "temperature": "Kelvin (K)",
            }

            st.markdown(
                f"""
                <div class="converter-result">
                <b>{fmt_decimal(alpha_value, sig_figs=6)}</b> alpha-units ({unit_type_rev})<br>
                = <b>{fmt_decimal(si_result, sig_figs=8)} {_SI_LABELS_REV.get(unit_type_rev, '')}</b>
                </div>
                """,
                unsafe_allow_html=True,
            )

        except (ValueError, ZeroDivisionError) as exc:
            st.error(f"Conversion error: {exc}")
    else:
        st.info("Enter a non-zero alpha-unit value to see the conversion.")

st.divider()

# ---------------------------------------------------------------------------
# Section 4: Notable Conversions
# ---------------------------------------------------------------------------
st.subheader("Notable Conversions")
st.markdown(
    "Pre-computed conversions of interesting physical quantities into alpha-units."
)

def _build_notable(consts):
    """Build notable conversions list from constants namespace."""
    _f = lambda v: float(v) if v is not None else None

    hbar = _f(consts.hbar)
    c = _f(consts.c)
    G = _f(consts.G)
    m_e = _f(consts.m_e)

    # Planck quantities (derived)
    l_Pl = (hbar * G / c ** 3) ** 0.5
    t_Pl = l_Pl / c
    E_Pl = (hbar * c ** 5 / G) ** 0.5

    return [
        {"name": "Planck length", "si_value": l_Pl, "unit_type": "length"},
        {"name": "Classical electron radius (rₑ)", "si_value": _f(consts.r_e_nist), "unit_type": "length"},
        {"name": "Bohr radius (a₀)", "si_value": _f(consts.a_0), "unit_type": "length"},
        {"name": "Compton wavelength (reduced)", "si_value": _f(consts.lambda_bar_c), "unit_type": "length"},
        {"name": "Proton mass", "si_value": _f(consts.m_p), "unit_type": "mass"},
        {"name": "Electron mass", "si_value": m_e, "unit_type": "mass"},
        {"name": "Muon mass", "si_value": _f(consts.m_mu), "unit_type": "mass"},
        {"name": "Electron rest energy", "si_value": m_e * c ** 2, "unit_type": "energy"},
        {"name": "Planck energy", "si_value": E_Pl, "unit_type": "energy"},
        {"name": "Planck time", "si_value": t_Pl, "unit_type": "time"},
        {"name": "Elementary charge", "si_value": _f(consts.e_charge), "unit_type": "charge"},
    ]


if constants is not None:
    _NOTABLE = _build_notable(constants)
else:
    # Fallback with CODATA 2018 hardcoded
    _NOTABLE = [
        {"name": "Planck length", "si_value": 1.616255e-35, "unit_type": "length"},
        {"name": "Classical electron radius (rₑ)", "si_value": 2.8179403227e-15, "unit_type": "length"},
        {"name": "Bohr radius (a₀)", "si_value": 5.2917721067e-11, "unit_type": "length"},
        {"name": "Compton wavelength (reduced)", "si_value": 3.8615926764e-13, "unit_type": "length"},
        {"name": "Proton mass", "si_value": 1.67262192369e-27, "unit_type": "mass"},
        {"name": "Electron mass", "si_value": 9.1093837015e-31, "unit_type": "mass"},
        {"name": "Muon mass", "si_value": 1.883531627e-28, "unit_type": "mass"},
        {"name": "Electron rest energy", "si_value": 8.1871057769e-14, "unit_type": "energy"},
        {"name": "Planck energy", "si_value": 1.9561e9, "unit_type": "energy"},
        {"name": "Planck time", "si_value": 5.39124e-44, "unit_type": "time"},
        {"name": "Elementary charge", "si_value": 1.602176634e-19, "unit_type": "charge"},
    ]

if units is not None:
    import pandas as pd

    notable_rows = []
    for item in _NOTABLE:
        try:
            alpha_val = _si_to_alpha_safe(item["si_value"], item["unit_type"])
            powers = _express_powers_safe(item["si_value"], item["unit_type"])

            notable_rows.append({
                "Quantity": item["name"],
                "SI Value": fmt_decimal(item["si_value"], sig_figs=6),
                "Unit Type": item["unit_type"].capitalize(),
                "Alpha-Unit Value": fmt_decimal(alpha_val, sig_figs=6),
                "Nearest αⁿ": f"α{powers['nearest_power']}" if powers else "N/A",
                "Exact Power": f"{powers['exact_power']:.3f}" if powers else "N/A",
                "Residual (%)": f"{powers['residual_pct']:.2f}" if powers else "N/A",
            })
        except (ValueError, ZeroDivisionError):
            notable_rows.append({
                "Quantity": item["name"],
                "SI Value": fmt_decimal(item["si_value"], sig_figs=6),
                "Unit Type": item["unit_type"].capitalize(),
                "Alpha-Unit Value": "Error",
                "Nearest αⁿ": "N/A",
                "Exact Power": "N/A",
                "Residual (%)": "N/A",
            })

    df_notable = pd.DataFrame(notable_rows)
    st.dataframe(
        df_notable,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Quantity": st.column_config.TextColumn(width="medium"),
            "SI Value": st.column_config.TextColumn(width="medium"),
            "Unit Type": st.column_config.TextColumn(width="small"),
            "Alpha-Unit Value": st.column_config.TextColumn(width="medium"),
            "Nearest αⁿ": st.column_config.TextColumn(width="small"),
            "Exact Power": st.column_config.TextColumn(width="small"),
            "Residual (%)": st.column_config.TextColumn(width="small"),
        },
    )

    st.info(
        "Quantities with low residual percentages sit near exact α-power rungs, "
        "suggesting they may be naturally expressible in the α-unit system. "
        "The electron mass mₑ and charge e are exact by construction (residual = 0%)."
    )
else:
    st.warning("Cannot compute notable conversions without the alpha unit system.")

st.divider()
st.caption("Alpha Units | Alpha Ladder Research Dashboard")
