"""
Experimental Strategies -- Page 9

Five testable strategies for the Alpha Ladder hypothesis,
each displayed as an expandable card with analysis tables and plots.
"""

import streamlit as st
import math
from decimal import Decimal, getcontext

getcontext().prec = 50


# ---------------------------------------------------------------------------
# Custom CSS
# ---------------------------------------------------------------------------
st.markdown(
    """
    <style>
    .strategy-header {
        font-family: 'Fira Mono', Consolas, monospace;
        font-size: 1rem;
    }
    .verdict-box {
        background-color: #1a1d23;
        border: 1px solid #2e3440;
        border-radius: 8px;
        padding: 1.2rem;
        margin: 0.5rem 0;
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
st.title("Experimental Strategies")
st.markdown("Five testable pathways for validating (or falsifying) the Alpha Ladder.")
st.divider()

# ---------------------------------------------------------------------------
# Try to load core modules
# ---------------------------------------------------------------------------
_core_available = False
try:
    from alpha_ladder_core.experimental import (  # noqa: E402
        strategy_mass_ratios,
        strategy_multiple_paths,
        strategy_dark_sector,
        strategy_muon_g2,
        strategy_experimental_approaches,
    )
    _core_available = True
except ImportError:
    pass

# ---------------------------------------------------------------------------
# Fallback constants for computations when core is unavailable
# ---------------------------------------------------------------------------
_FALLBACK_ALPHA = 0.0072973525693

if constants is not None and hasattr(constants, "alpha"):
    _alpha = float(constants.alpha)
else:
    _alpha = _FALLBACK_ALPHA


# ---------------------------------------------------------------------------
# Strategy 1: Mass Ratios
# ---------------------------------------------------------------------------
with st.expander("Strategy 1: Mass Ratios -- Do particle masses fall on α rungs?", expanded=False):
    st.markdown(
        "If the Alpha Ladder is fundamental, particle mass ratios should "
        "cluster near integer powers of 1/α."
    )

    if _core_available and constants is not None:
        ratios = strategy_mass_ratios(constants)
    else:
        # Minimal fallback
        m_e = 9.1093837015e-31
        m_p = 1.67262192369e-27
        m_mu = 1.883531627e-28
        raw_ratios = {
            "mₚ / mₑ": m_p / m_e,
            "mμ / mₑ": m_mu / m_e,
            "mₚ / mμ": m_p / m_mu,
        }
        ratios = []
        for name, ratio in raw_ratios.items():
            n_exact = math.log(ratio) / math.log(1 / _alpha)
            n_round = round(n_exact)
            coeff = ratio * _alpha ** n_round
            ratios.append({
                "name": name,
                "ratio": ratio,
                "n_exact": n_exact,
                "n_round": n_round,
                "coeff_at_round": coeff,
                "nearby_rungs": [],
            })

    import pandas as pd

    rows = []
    for r in ratios:
        rows.append({
            "Mass Ratio": r["name"],
            "Value": f"{r['ratio']:.6f}",
            "Exact rung n": f"{r['n_exact']:.4f}",
            "Nearest integer n": r["n_round"],
            "Coefficient at nearest": f"{r['coeff_at_round']:.6f}",
        })

    df = pd.DataFrame(rows)
    st.dataframe(df, use_container_width=True, hide_index=True)

    st.caption(
        "A coefficient near 1 at the nearest integer rung indicates a clean "
        "\u03b1-power relationship."
    )


# ---------------------------------------------------------------------------
# Strategy 2: Multiple Paths to \u03b1\u1d33
# ---------------------------------------------------------------------------
with st.expander("Strategy 2: Multiple Paths to \u03b1\u1d33 -- Convergence test", expanded=False):
    st.markdown(
        "Three independent routes to the gravitational coupling constant "
        "should converge if the ladder is correct."
    )

    if _core_available and constants is not None:
        paths = strategy_multiple_paths(constants)
    else:
        alpha_d = Decimal(str(_alpha))
        phi_d = (1 + Decimal(5).sqrt()) / Decimal(2)
        hbar = Decimal('1.054571817e-34')
        c = Decimal('2.99792458e8')
        m_e = Decimal('9.1093837015e-31')
        m_p = Decimal('1.67262192369e-27')
        G = Decimal('6.67430e-11')

        alpha_21 = alpha_d ** 21
        alpha_g_A = (phi_d ** 2 / Decimal(2)) * alpha_21
        G_A = alpha_g_A * hbar * c / m_e ** 2

        alpha_g_proton = G * m_p ** 2 / (hbar * c)
        n_proton = math.log(float(alpha_g_proton)) / math.log(float(alpha_d))

        hbar_c_over_G = hbar * c / G
        m_Pl = hbar_c_over_G.sqrt()
        ratio_e_Pl = m_e / m_Pl
        n_planck = math.log(float(ratio_e_Pl)) / math.log(float(alpha_d))

        mu = m_p / m_e
        alpha_g_D = alpha_d ** 24 * mu ** 2
        G_D = alpha_g_D * hbar * c / m_e ** 2

        paths = {
            "path_A": {
                "description": "\u03c6\u00b2/2 \u00b7 \u03b1\u00b2\u00b9",
                "alpha_G": float(alpha_g_A),
                "G_predicted": float(G_A),
            },
            "path_B": {
                "description": "proton-based \u03b1_G",
                "alpha_G_proton": float(alpha_g_proton),
                "n_proton": n_proton,
                "nearest_rung": round(n_proton),
                "coeff_at_rung": float(alpha_g_proton) / float(alpha_d) ** round(n_proton),
            },
            "path_C": {
                "description": "Planck mass relationship",
                "m_e_over_m_Pl": float(ratio_e_Pl),
                "alpha_G_from_planck": float(ratio_e_Pl ** 2),
                "n_planck": n_planck,
                "n_alpha_G": 2 * n_planck,
            },
            "path_D": {
                "description": "alpha^24 * mu^2 (zero free parameters)",
                "alpha_G": float(alpha_g_D),
                "G_predicted": float(G_D),
                "residual_ppm": float((G_D - G) / G * Decimal('1000000')),
            },
        }

    col_a, col_b, col_c, col_d = st.columns(4)

    with col_a:
        st.markdown("**Path A: Bridge Formula**")
        st.markdown("\u03b1_G = \u03c6\u00b2/2 \u00b7 \u03b1\u00b2\u00b9")
        st.metric(
            label="\u03b1_G (Path A)",
            value=fmt_decimal(paths["path_A"]["alpha_G"], sig_figs=6),
        )
        st.metric(
            label="G predicted (m\u00b3 kg\u207b\u00b9 s\u207b\u00b2)",
            value=fmt_decimal(paths["path_A"]["G_predicted"], sig_figs=6),
        )

    with col_b:
        st.markdown("**Path B: Proton-Based**")
        st.markdown("\u03b1_G(proton) = G \u00b7 m\u209a\u00b2 / (\u210f \u00b7 c)")
        st.metric(
            label="\u03b1_G (proton)",
            value=fmt_decimal(paths["path_B"]["alpha_G_proton"], sig_figs=6),
        )
        st.metric(
            label="Nearest rung",
            value=str(paths["path_B"]["nearest_rung"]),
        )
        st.metric(
            label="Coefficient at rung",
            value=f"{paths['path_B']['coeff_at_rung']:.4f}",
        )

    with col_c:
        st.markdown("**Path C: Planck Mass**")
        st.markdown("m\u2091 / m\u209a\u2097 \u2192 \u03b1\u207f")
        st.metric(
            label="m\u2091 / m\u209a\u2097",
            value=fmt_decimal(paths["path_C"]["m_e_over_m_Pl"], sig_figs=6),
        )
        st.metric(
            label="n (from Planck)",
            value=f"{paths['path_C']['n_planck']:.4f}",
        )
        st.metric(
            label="n for \u03b1_G = (m\u2091/m\u209a\u2097)\u00b2",
            value=f"{paths['path_C']['n_alpha_G']:.4f}",
        )

    with col_d:
        st.markdown("**Path D: Hierarchy Formula**")
        st.markdown("alpha_G = alpha^24 * mu^2")
        st.metric(
            label="\u03b1_G (Path D)",
            value=fmt_decimal(paths["path_D"]["alpha_G"], sig_figs=6),
        )
        st.metric(
            label="G predicted",
            value=fmt_decimal(paths["path_D"]["G_predicted"], sig_figs=6),
        )
        st.metric(
            label="Residual (ppm)",
            value=f"{paths['path_D']['residual_ppm']:.0f}",
        )

    st.info(
        "All four paths should yield \u03b1_G values consistent with "
        "the ladder exponent. Path D (alpha^24 * mu^2) uses zero free "
        "parameters. Convergence strengthens the hypothesis; divergence "
        "would falsify it."
    )


# ---------------------------------------------------------------------------
# Strategy 3: Dark Sector
# ---------------------------------------------------------------------------
with st.expander("Strategy 3: Dark Sector -- Dark photon coupling prediction", expanded=False):
    st.markdown(
        "The ladder predicts a dark photon kinetic mixing parameter "
        "\u03b5 = \u221a(\u03b1\u00b9\u2070)."
    )

    if _core_available and constants is not None:
        dark = strategy_dark_sector(constants)
    else:
        alpha_10 = _alpha ** 10
        epsilon = math.sqrt(alpha_10)
        below_current = epsilon < 1e-3
        dark = {
            "alpha_10": alpha_10,
            "epsilon_predicted": epsilon,
            "experimental_bounds": {
                "BaBar_2017": {"epsilon_max": 1e-3, "mass_range": "1-10 GeV"},
                "NA64_2019": {"epsilon_max": 1e-4, "mass_range": "1-100 MeV"},
                "LDMX_projected": {"epsilon_max": 1e-6, "mass_range": "1-100 MeV"},
            },
            "below_current_bounds": below_current,
            "orders_below_sensitivity": math.log10(1e-3 / epsilon) if below_current else None,
        }

    col_d1, col_d2 = st.columns(2)

    with col_d1:
        st.metric(label="\u03b1\u00b9\u2070", value=fmt_decimal(dark["alpha_10"], sig_figs=6))
        st.metric(
            label="\u03b5 (predicted mixing)",
            value=fmt_decimal(dark["epsilon_predicted"], sig_figs=6),
        )

        if dark["below_current_bounds"]:
            st.success(
                f"Prediction is below current experimental bounds "
                f"by ~{dark['orders_below_sensitivity']:.1f} orders of magnitude."
            )
        else:
            st.error("Prediction is above current bounds -- potentially already excluded.")

    with col_d2:
        st.markdown("**Experimental Bounds**")

        import pandas as pd

        bound_rows = []
        for exp_name, bdata in dark["experimental_bounds"].items():
            bound_rows.append({
                "Experiment": exp_name.replace("_", " "),
                "ε_max": fmt_decimal(bdata["epsilon_max"], sig_figs=2),
                "Mass Range": bdata["mass_range"],
            })

        df_bounds = pd.DataFrame(bound_rows)
        st.dataframe(df_bounds, use_container_width=True, hide_index=True)

    # Simple plot: prediction vs bounds
    st.markdown("**Prediction vs Experimental Bounds**")

    try:
        import plotly.graph_objects as go
        from app.components.charts import _apply_theme, _COLOR_GREEN, _COLOR_RED, _COLOR_HIGHLIGHT_21

        experiments = ["BaBar 2017", "NA64 2019", "LDMX (proj.)"]
        eps_bounds = [1e-3, 1e-4, 1e-6]
        prediction_val = dark["epsilon_predicted"]

        fig = go.Figure()

        fig.add_trace(go.Bar(
            x=experiments,
            y=eps_bounds,
            name="Upper Bound",
            marker_color=_COLOR_RED,
            opacity=0.7,
        ))

        fig.add_hline(
            y=prediction_val,
            line_dash="dash",
            line_color=_COLOR_GREEN,
            line_width=2,
            annotation_text=f"Ladder prediction: {prediction_val:.2e}",
            annotation_position="top left",
            annotation_font_color=_COLOR_GREEN,
        )

        fig.update_layout(
            yaxis_type="log",
            yaxis_title="ε (kinetic mixing)",
            xaxis_title="Experiment",
            title="Dark Photon Mixing: Prediction vs Bounds",
            height=400,
        )
        fig = _apply_theme(fig)
        st.plotly_chart(fig, use_container_width=True)

    except Exception:
        st.warning("Could not render dark sector plot (Plotly not available).")


# ---------------------------------------------------------------------------
# Strategy 4: Muon g-2
# ---------------------------------------------------------------------------
with st.expander("Strategy 4: Muon g-2 -- Anomaly analysis", expanded=False):
    st.markdown(
        "The muon anomalous magnetic moment discrepancy may align with "
        "an alpha-power rung."
    )

    if _core_available and constants is not None:
        g2 = strategy_muon_g2(constants)
    else:
        alpha_d = Decimal(str(_alpha))
        delta_a_mu = 2.51e-9
        a_mu_exp = 1.16592061e-3
        n_anomaly = math.log(delta_a_mu) / math.log(float(alpha_d))
        nearest_rung_anomaly = round(n_anomaly)
        anomaly_coeff = delta_a_mu / float(alpha_d) ** nearest_rung_anomaly
        frac_anomaly = delta_a_mu / a_mu_exp
        n_frac = math.log(frac_anomaly) / math.log(float(alpha_d))
        nearest_rung_frac = round(n_frac)

        g2 = {
            "delta_a_mu": delta_a_mu,
            "a_mu_exp": a_mu_exp,
            "fractional_anomaly": frac_anomaly,
            "n_anomaly": n_anomaly,
            "nearest_rung_anomaly": nearest_rung_anomaly,
            "anomaly_coeff": anomaly_coeff,
            "n_fractional": n_frac,
            "nearest_rung_fractional": nearest_rung_frac,
        }

    col_g1, col_g2 = st.columns(2)

    with col_g1:
        st.markdown("**Anomaly Magnitude**")
        st.metric(
            label="\u0394a\u03bc (experimental discrepancy)",
            value=fmt_decimal(g2["delta_a_mu"], sig_figs=4),
        )
        st.metric(
            label="Exact \u03b1-power",
            value=f"{g2['n_anomaly']:.4f}",
        )
        st.metric(
            label="Nearest integer rung",
            value=str(g2["nearest_rung_anomaly"]),
        )
        st.metric(
            label="Coefficient at nearest rung",
            value=f"{g2['anomaly_coeff']:.4f}",
        )

    with col_g2:
        st.markdown("**Fractional Anomaly**")
        st.metric(
            label="\u0394a\u03bc / a\u03bc(exp)",
            value=fmt_decimal(g2["fractional_anomaly"], sig_figs=4),
        )
        st.metric(
            label="Exact \u03b1-power (fractional)",
            value=f"{g2['n_fractional']:.4f}",
        )
        st.metric(
            label="Nearest integer rung (fractional)",
            value=str(g2["nearest_rung_fractional"]),
        )

    st.info(
        "If the anomaly magnitude or its fractional value sits near an "
        "integer \u03b1-power rung, it suggests the new physics responsible "
        "for the discrepancy may share the same \u03b1-scaling structure."
    )


# ---------------------------------------------------------------------------
# Strategy 5: Experimental Approaches
# ---------------------------------------------------------------------------
with st.expander("Strategy 5: Experimental Approaches to G -- Precision tests", expanded=False):
    st.markdown(
        "The Alpha Ladder predicts a specific value of Newton's constant G. "
        "Which experiments can test this?"
    )

    if _core_available and constants is not None:
        approaches = strategy_experimental_approaches(constants)
    else:
        approaches = {
            "ladder_prediction": 6.67323e-11,
            "hierarchy_prediction": 6.67889e-11,
            "hierarchy_difference_ppm": 688,
            "codata_2018": {"value": 6.67430e-11, "uncertainty": 0.00015e-11},
            "difference": 1.07e-14,
            "difference_ppm": 16,
            "required_precision_ppm": 5,
            "approaches": {
                "A_atom_interferometry": {
                    "description": "Atom interferometry (most promising)",
                    "rosi_2014_value": 6.67191e-11,
                    "technique": "Cold atom clouds as test masses",
                    "projected_precision_ppm": "5-10",
                    "groups": ["Stanford (Kasevich)", "Florence (Tino)", "Wuhan (Zhu)"],
                },
                "B_MEMS_oscillators": {
                    "description": "MEMS oscillators",
                    "technique": "Microscale torsion oscillators on silicon chips",
                    "current_precision_ppm": 100,
                },
                "C_targeted_reanalysis": {
                    "description": "Targeted re-analysis of Rosi/Tino data",
                    "rosi_uncertainty_ppm": 99,
                    "prediction_sigma_from_rosi": 1.3,
                    "note": "Only existing measurement consistent with the ladder",
                },
            },
        }

    col_e1, col_e2, col_e3, col_e4 = st.columns(4)

    with col_e1:
        st.metric(
            label="G (bridge)",
            value=fmt_decimal(approaches["ladder_prediction"], sig_figs=6),
            delta=f"{approaches['difference_ppm']} ppm",
            delta_color="off",
        )

    with col_e2:
        st.metric(
            label="G (alpha^24 * mu^2)",
            value=fmt_decimal(approaches["hierarchy_prediction"], sig_figs=6),
            delta=f"{approaches['hierarchy_difference_ppm']} ppm",
            delta_color="off",
        )

    with col_e3:
        st.metric(
            label="G (CODATA 2018)",
            value=fmt_decimal(approaches["codata_2018"]["value"], sig_figs=6),
        )

    with col_e4:
        st.metric(
            label="Required precision",
            value=f"{approaches['required_precision_ppm']} ppm",
            help="Precision needed to distinguish bridge prediction from CODATA",
        )

    st.markdown("---")

    for key, approach in approaches["approaches"].items():
        st.markdown(f"**{approach['description']}**")

        if "technique" in approach:
            st.markdown(f"- Technique: {approach['technique']}")

        if "projected_precision_ppm" in approach:
            st.markdown(f"- Projected precision: {approach['projected_precision_ppm']} ppm")

        if "current_precision_ppm" in approach:
            st.markdown(f"- Current precision: {approach['current_precision_ppm']} ppm")

        if "groups" in approach:
            st.markdown(f"- Active groups: {', '.join(approach['groups'])}")

        if "rosi_2014_value" in approach:
            st.markdown(
                f"- Rosi 2014 value: {fmt_decimal(approach['rosi_2014_value'], sig_figs=6)} "
                f"m\u00b3 kg\u207b\u00b9 s\u207b\u00b2"
            )

        if "prediction_sigma_from_rosi" in approach:
            st.markdown(
                f"- Ladder prediction is {approach['prediction_sigma_from_rosi']}\u03c3 "
                f"from Rosi measurement"
            )

        if "rosi_uncertainty_ppm" in approach:
            st.markdown(f"- Rosi published uncertainty: \u00b1{approach['rosi_uncertainty_ppm']} ppm")

        if "note" in approach:
            st.caption(approach["note"])

        st.markdown("")

st.divider()

# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------
st.subheader("Verdict: What To Do Next")

col_v1, col_v2, col_v3, col_v4 = st.columns(4)

with col_v1:
    st.markdown(
        """
        <div class="verdict-box">
        <b>Cheapest</b><br><br>
        Re-analyze existing Rosi/Tino atom interferometry data with the
        ladder prediction as the null hypothesis. Requires no new
        equipment -- only computational effort.
        </div>
        """,
        unsafe_allow_html=True,
    )

with col_v2:
    st.markdown(
        """
        <div class="verdict-box">
        <b>Fastest</b><br><br>
        Check particle mass ratios against alpha-power rungs using
        existing high-precision mass measurements from Penning trap
        experiments. Can be done immediately.
        </div>
        """,
        unsafe_allow_html=True,
    )

with col_v3:
    st.markdown(
        """
        <div class="verdict-box">
        <b>Strongest</b><br><br>
        Next-generation atom interferometry at 5 ppm precision
        (Stanford/Florence/Wuhan). Would provide a definitive test
        of the G prediction within the decade.
        </div>
        """,
        unsafe_allow_html=True,
    )

with col_v4:
    st.markdown(
        """
        <div class="verdict-box">
        <b>Boldest</b><br><br>
        Search for dark photon at epsilon = \u03b1\u2075 mixing. Would
        require next-generation intensity-frontier experiments
        (LDMX, Belle II). Discovery would be transformative.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.divider()
st.caption("Experimental Strategies | Alpha Ladder Research Dashboard")
