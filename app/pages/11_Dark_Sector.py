"""
Dark Sector (Rung 10) -- Page 11

Interactive exploration of alpha^10 as an ultra-light axion-like particle:
Fuzzy Wave Visualizer, Interaction Gap, and ALPS II photon-mixing simulator.
"""

import streamlit as st
import math
from decimal import Decimal, getcontext

getcontext().prec = 50

st.set_page_config(page_title="Dark Sector | Alpha Ladder", layout="wide")

# ---------------------------------------------------------------------------
# Custom CSS
# ---------------------------------------------------------------------------
st.markdown(
    """
    <style>
    .dark-header {
        font-family: 'Fira Mono', Consolas, monospace;
        font-size: 1rem;
    }
    .dark-card {
        background-color: #1a1d23;
        border: 1px solid #2e3440;
        border-radius: 8px;
        padding: 1.2rem;
        margin: 0.5rem 0;
    }
    .alps-step {
        background-color: #1a1d23;
        border-left: 3px solid #a78bfa;
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
from app.components.formatting import fmt_decimal  # noqa: E402

constants = render_sidebar()

# ---------------------------------------------------------------------------
# Try to load core modules
# ---------------------------------------------------------------------------
_core_available = False
try:
    from alpha_ladder_core.dark_sector import (  # noqa: E402
        compute_dark_sector,
        compute_wave_profile,
        compute_alps_simulation,
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
# Header
# ---------------------------------------------------------------------------
st.title("The Dark Sector: Rung 10")
st.divider()

# ---------------------------------------------------------------------------
# Compute dark sector values
# ---------------------------------------------------------------------------
if _core_available and constants is not None:
    dark = compute_dark_sector(constants)
else:
    alpha_10 = _alpha ** 10
    epsilon = _alpha ** 5
    m_e_kg = 9.1093837015e-31
    hbar_Js = 1.054571817e-34
    c_ms = 2.99792458e8
    e_charge = 1.602176634e-19
    m_axion_eV = m_e_kg * c_ms ** 2 / e_charge * alpha_10
    lambda_bar_m = hbar_Js / (m_e_kg * alpha_10 * c_ms)
    lambda_bar_km = lambda_bar_m / 1000.0
    earth_moon_km = 384400.0
    lambda_ratio = lambda_bar_km / earth_moon_km
    ratio_to_em = alpha_10 / _alpha  # alpha^9
    # Fallback alpha_g estimate
    phi = (1 + math.sqrt(5)) / 2
    alpha_g_fallback = (phi ** 2 / 2) * _alpha ** 21
    ratio_to_gravity = alpha_10 / alpha_g_fallback
    log10_em = math.log10(_alpha)
    log10_dark = math.log10(alpha_10)
    log10_grav = math.log10(alpha_g_fallback)

    dark = {
        "alpha_10": alpha_10,
        "epsilon": epsilon,
        "m_axion_eV": m_axion_eV,
        "lambda_bar_km": lambda_bar_km,
        "earth_moon_km": earth_moon_km,
        "lambda_ratio": lambda_ratio,
        "ratio_to_em": ratio_to_em,
        "ratio_to_gravity": ratio_to_gravity,
        "log10_em": log10_em,
        "log10_dark": log10_dark,
        "log10_grav": log10_grav,
    }

# ---------------------------------------------------------------------------
# 4 Metric Cards
# ---------------------------------------------------------------------------
col_m1, col_m2, col_m3, col_m4 = st.columns(4)

with col_m1:
    st.metric(
        label="\u03b1\u00b9\u2070",
        value=fmt_decimal(dark["alpha_10"], sig_figs=6),
    )

with col_m2:
    st.metric(
        label="\u03b5 (kinetic mixing)",
        value=fmt_decimal(dark["epsilon"], sig_figs=6),
    )

with col_m3:
    st.metric(
        label="Axion mass (eV)",
        value=fmt_decimal(dark["m_axion_eV"], sig_figs=6),
    )

with col_m4:
    st.metric(
        label="Wavelength",
        value=f"{float(dark['lambda_bar_km']):,.0f} km",
    )

st.markdown("")

# ---------------------------------------------------------------------------
# Section A: Fuzzy Wave Visualizer
# ---------------------------------------------------------------------------
with st.expander("Fuzzy Wave Visualizer -- Rung 10 Density Profile", expanded=True):

    col_chart, col_controls = st.columns([2, 1])

    with col_controls:
        _alpha_real = float(constants.alpha) if (constants is not None and hasattr(constants, "alpha")) else _FALLBACK_ALPHA

        _alpha_options = {
            f"0.0050000  (weak coupling)": 0.005,
            f"{_alpha_real:.7f}  (real \u03b1)": _alpha_real,
            f"0.0100000  (strong coupling)": 0.01,
        }
        _labels = list(_alpha_options.keys())

        selected_label = st.select_slider(
            "Fine-structure constant (\u03b1)",
            options=_labels,
            value=_labels[1],
            key="alpha_slider",
        )
        alpha_val = _alpha_options[selected_label]

        if _core_available and constants is not None:
            wave = compute_wave_profile(constants, alpha_override=Decimal(str(alpha_val)))
        else:
            # Fallback wave computation
            alpha_used = alpha_val
            lam_m = hbar_Js / (m_e_kg * alpha_used ** 10 * c_ms)
            lam_km = lam_m / 1000.0
            x_max = 3.0 * 384400.0
            n_pts = 500
            x_km_list = [x_max * i / (n_pts - 1) for i in range(n_pts)]
            two_pi = 2.0 * math.pi
            psi_sq = [math.cos(two_pi * x / lam_km) ** 2 for x in x_km_list]
            wave = {
                "lambda_bar_km": Decimal(str(lam_km)),
                "x_km": x_km_list,
                "psi_squared": psi_sq,
            }

        earth_moon_ratio = float(wave["lambda_bar_km"]) / 384400.0

        st.metric(
            label="Wavelength (km)",
            value=f"{float(wave['lambda_bar_km']):,.0f}",
        )
        st.metric(
            label="Earth-Moon ratio",
            value=f"{earth_moon_ratio:.2f}",
        )

    with col_chart:
        try:
            import plotly.graph_objects as go
            from app.components.charts import _apply_theme, _COLOR_HIGHLIGHT_10, _COLOR_DEFAULT

            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=wave["x_km"], y=wave["psi_squared"],
                mode="lines", fill="tozeroy",
                line=dict(color="#a78bfa", width=2),
                fillcolor="rgba(167, 139, 250, 0.2)",
                name="|psi|^2",
            ))
            # Vertical dashed line at Earth-Moon distance
            fig.add_vline(
                x=384400, line_dash="dash", line_color="#60a5fa", line_width=2,
                annotation_text="Earth-Moon", annotation_position="top right",
                annotation_font_color="#60a5fa",
            )
            fig.update_layout(
                title="Fuzzy Dark Matter Wave Envelope",
                xaxis_title="Distance (km)",
                yaxis_title="|psi|^2",
                height=450,
            )
            fig = _apply_theme(fig)
            st.plotly_chart(fig, use_container_width=True)

        except Exception:
            st.warning("Could not render wave profile chart (Plotly not available).")

# ---------------------------------------------------------------------------
# Section B: Interaction Gap
# ---------------------------------------------------------------------------
with st.expander("Interaction Gap -- Where Does Rung 10 Sit?", expanded=True):

    try:
        import plotly.graph_objects as go
        from app.components.charts import _apply_theme, _COLOR_HIGHLIGHT_10, _COLOR_DEFAULT

        fig2 = go.Figure()
        categories = ["EM (\u03b1)", "Rung 10 (\u03b1\u00b9\u2070)", "Gravity (\u03b1_g)"]
        log_values = [dark["log10_em"], dark["log10_dark"], dark["log10_grav"]]
        colors = ["#60a5fa", "#a78bfa", "#f59e0b"]

        fig2.add_trace(go.Bar(
            y=categories, x=log_values,
            orientation="h",
            marker_color=colors,
            text=[f"{v:.1f}" for v in log_values],
            textposition="outside",
        ))

        # Add annotations for gaps
        gap_em_dark = abs(dark["log10_dark"] - dark["log10_em"])
        gap_dark_grav = abs(dark["log10_grav"] - dark["log10_dark"])

        fig2.add_annotation(
            x=(dark["log10_em"] + dark["log10_dark"]) / 2, y=0.5,
            text=f"~10^{gap_em_dark:.0f} gap", showarrow=False,
            font=dict(color="#e0e0e0", size=12),
        )
        fig2.add_annotation(
            x=(dark["log10_dark"] + dark["log10_grav"]) / 2, y=1.5,
            text=f"~10^{gap_dark_grav:.0f} gap", showarrow=False,
            font=dict(color="#e0e0e0", size=12),
        )

        fig2.update_layout(
            title="Coupling Strength Hierarchy (log10 scale)",
            xaxis_title="log10(coupling)",
            height=350,
        )
        fig2 = _apply_theme(fig2)
        st.plotly_chart(fig2, use_container_width=True)

    except Exception:
        st.warning("Could not render interaction gap chart (Plotly not available).")

    # Two-column metrics for exact ratio values
    col_r1, col_r2 = st.columns(2)

    with col_r1:
        st.metric(
            label="\u03b1\u00b9\u2070 / \u03b1 = \u03b1\u2079",
            value=fmt_decimal(dark["ratio_to_em"], sig_figs=6),
        )

    with col_r2:
        st.metric(
            label="\u03b1\u00b9\u2070 / \u03b1_g",
            value=fmt_decimal(dark["ratio_to_gravity"], sig_figs=6),
        )

# ---------------------------------------------------------------------------
# Section C: Axion-Photon Mixer (ALPS II Lab)
# ---------------------------------------------------------------------------
with st.expander("Axion-Photon Mixer -- ALPS II Simulation", expanded=True):

    st.markdown(
        "The ALPS II (Any Light Particle Search) experiment at DESY uses "
        "a **Light Shining through a Wall** (LSW) technique: a high-power "
        "laser beam passes through a strong magnetic field where photons "
        "may convert to axion-like particles, traverse an opaque barrier, "
        "then reconvert to detectable photons in a second magnetic field region."
    )

    st.markdown("")
    col_p1, col_p2 = st.columns(2)
    with col_p1:
        st.metric(label="Magnetic field (B)", value="5.3 T")
    with col_p2:
        st.metric(label="Cavity length (L)", value="106 m")

    st.markdown("")

    if st.button("Fire Photon", key="fire_photon", type="primary"):
        st.session_state["photon_fired"] = True

    if st.session_state.get("photon_fired", False):
        if _core_available and constants is not None:
            alps = compute_alps_simulation(constants)
        else:
            eps_fb = _alpha ** 5
            P_fb = eps_fb ** 4
            photons_fb = 1.0 / P_fb
            pps_fb = 30.0 / (1.165 * 1.602176634e-19)
            wait_s = photons_fb / pps_fb
            wait_y = wait_s / 3.156e7
            alps = {
                "epsilon": eps_fb,
                "P_tunneling": P_fb,
                "photons_needed": photons_fb,
                "B_tesla": 5.3,
                "L_meters": 106,
                "photons_per_sec": pps_fb,
                "wait_years": wait_y,
                "log10_wait_years": math.log10(wait_y),
                "alps_ii_sensitivity": 3e-4,
                "orders_below_alps": math.log10(3e-4) - math.log10(eps_fb),
            }

        _eps_str = fmt_decimal(float(alps["epsilon"]), sig_figs=6)
        _p_str = fmt_decimal(float(alps["P_tunneling"]), sig_figs=6)
        _n_str = fmt_decimal(float(alps["photons_needed"]), sig_figs=4)

        # Derive laser context from base ALPS values (resilient to module version)
        _photons_needed_f = float(alps["photons_needed"])
        _laser_pps = 30.0 / (1.165 * 1.602176634e-19)  # 30 W, 1064 nm
        _wait_yrs = _photons_needed_f / _laser_pps / 3.156e7
        _log_wait = f"{math.log10(_wait_yrs):.0f}"
        _pps_str = fmt_decimal(_laser_pps, sig_figs=3)
        _eps_f = float(alps["epsilon"])
        _orders_below = f"{math.log10(3e-4) - math.log10(_eps_f):.1f}"

        st.markdown(
            f"""
            <div class="alps-step">
            <b>Step 1:</b> Photon enters 5.3 T magnetic field (ALPS II: 30 W laser at 1064 nm
            producing ~{_pps_str} photons/sec)
            </div>
            <div class="alps-step">
            <b>Step 2:</b> Mixing probability &epsilon; = &alpha;<sup>5</sup> = {_eps_str}
            &mdash; the photon-to-axion conversion amplitude at Rung 5
            </div>
            <div class="alps-step">
            <b>Step 3:</b> Round-trip tunneling P = &epsilon;<sup>4</sup> = &alpha;<sup>20</sup> = {_p_str}
            &mdash; two conversions (photon &rarr; axion &rarr; photon), each at &epsilon;<sup>2</sup>
            </div>
            <div class="alps-step">
            <b>Result:</b> At {_pps_str} photons/sec, you would need to run for ~10<sup>{_log_wait}</sup> years
            to expect a single detection ({_n_str} photons required). The age of the universe
            is ~10<sup>10</sup> years.
            </div>
            """,
            unsafe_allow_html=True,
        )

        st.markdown("")

        st.markdown(
            f"""
            <div class="dark-card">
            <b>Why no detection?</b> Rung 10 sits {_orders_below} orders of magnitude below
            ALPS II's design sensitivity (&epsilon; ~ 3 &times; 10<sup>-4</sup>). This is not a
            failure of the experiment -- it means the alpha-ladder dark photon at &epsilon; = &alpha;<sup>5</sup>
            ~ 2 &times; 10<sup>-11</sup> lives in a regime that <em>no</em> current LSW experiment can reach.
            Next-generation proposals (JURA, ALPS III) aim for &epsilon; ~ 10<sup>-6</sup>, still
            five orders short. Indirect astrophysical probes (stellar cooling, CMB distortions)
            may be the first to test this prediction.
            </div>
            """,
            unsafe_allow_html=True,
        )

    st.markdown("")

    # Prediction vs experimental bounds chart
    try:
        import plotly.graph_objects as go
        from app.components.charts import _apply_theme, _COLOR_GREEN, _COLOR_RED

        # Get experimental bounds
        dark_exp = None
        try:
            from alpha_ladder_core.experimental import strategy_dark_sector as _strat_dark
            if _core_available and constants is not None:
                dark_exp = _strat_dark(constants)
        except ImportError:
            pass

        if dark_exp is None:
            dark_exp = {
                "epsilon_predicted": _alpha ** 5,
                "experimental_bounds": {
                    "BaBar_2017": {"epsilon_max": 1e-3, "mass_range": "1-10 GeV"},
                    "NA64_2019": {"epsilon_max": 1e-4, "mass_range": "1-100 MeV"},
                    "LDMX_projected": {"epsilon_max": 1e-6, "mass_range": "1-100 MeV"},
                },
            }

        experiments = ["BaBar 2017", "NA64 2019", "LDMX (proj.)"]
        eps_bounds = [1e-3, 1e-4, 1e-6]
        prediction_val = float(dark["epsilon"])

        fig3 = go.Figure()
        fig3.add_trace(go.Bar(
            x=experiments, y=eps_bounds, name="Upper Bound",
            marker_color=_COLOR_RED, opacity=0.7,
        ))
        fig3.add_hline(
            y=prediction_val, line_dash="dash", line_color=_COLOR_GREEN, line_width=2,
            annotation_text=f"Rung 5 prediction: {prediction_val:.2e}",
            annotation_position="top left", annotation_font_color=_COLOR_GREEN,
        )
        fig3.update_layout(
            yaxis_type="log",
            yaxis_title="\u03b5 (kinetic mixing)",
            xaxis_title="Experiment",
            title="ALPS II Sensitivity vs Ladder Prediction",
            height=400,
        )
        fig3 = _apply_theme(fig3)
        st.plotly_chart(fig3, use_container_width=True)

    except Exception:
        st.warning("Could not render ALPS II sensitivity chart (Plotly not available).")

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Dark Sector (Rung 10) | Alpha Ladder Research Dashboard")
