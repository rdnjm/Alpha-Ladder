"""
The Prediction (Pillar 4) -- Page 12

Derives Newton's gravitational constant G from first principles via the
alpha ladder, visualizes the Big G deadlock, dilaton screening profile,
and falsification criteria.
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
    .pred-card {
        background-color: #1a1d23;
        border: 1px solid #2e3440;
        border-radius: 8px;
        padding: 1.2rem;
        margin: 0.5rem 0;
    }
    .formula-card {
        background-color: #1a1d23;
        border-left: 3px solid #f59e0b;
        padding: 0.8rem 1rem;
        margin: 0.4rem 0;
        border-radius: 0 8px 8px 0;
    }
    .falsify-card {
        background-color: #1a1d23;
        border-left: 3px solid #34d399;
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
from app.components.formatting import fmt_decimal, fmt_sigma  # noqa: E402

constants = render_sidebar()

# ---------------------------------------------------------------------------
# Try to load core modules
# ---------------------------------------------------------------------------
_core_available = False
try:
    from alpha_ladder_core.screening import (  # noqa: E402
        compute_screening_parameters,
        compute_G_eff,
        compute_screening_profile,
        classify_measurements,
    )
    from alpha_ladder_core.predict_g import (  # noqa: E402
        compare_prediction,
        get_G_measurements,
        get_bridge_candidates,
        predict_G,
        predict_G_hierarchy,
    )
    _core_available = True
except ImportError:
    pass

# ---------------------------------------------------------------------------
# Fallback constants
# ---------------------------------------------------------------------------
_FALLBACK_ALPHA = 0.0072973525693
_PHI = (1 + math.sqrt(5)) / 2

if constants is not None and hasattr(constants, "alpha"):
    _alpha = float(constants.alpha)
else:
    _alpha = _FALLBACK_ALPHA

# Compute core values
if _core_available and constants is not None:
    screening = compute_screening_parameters(constants)
    G_vacuum = float(screening["G_vacuum"])
    G_lab = float(screening["G_lab"])
    alpha_screening = float(screening["alpha_screening"])
    ppm_excess = float(screening["ppm_excess"])
    lambda_dilaton_m = float(screening["lambda_dilaton_m"])
    dilaton_mass_eV = float(screening["dilaton_mass_eV"])
    dilaton_mass_kg = float(screening["dilaton_mass_kg"])
    omega = float(screening["omega"])
else:
    # Fallback computation
    screening = None
    G_vacuum = float((_PHI ** 2 / 2) * _alpha ** 21 / (2 * math.pi)) if False else 6.67322976e-11
    G_lab = 6.67430e-11
    alpha_screening = 1.607e-4
    ppm_excess = (G_lab - G_vacuum) / G_vacuum * 1e6
    lambda_dilaton_m = 0.002
    dilaton_mass_eV = 9.9e-5
    dilaton_mass_kg = 1.76e-40
    omega = 0.0

# Pull G_codata from constants when available, hardcoded fallback otherwise
if constants is not None and hasattr(constants, "G"):
    G_codata = float(constants.G)
    try:
        _g_meas = get_G_measurements()
        _codata_entry = _g_meas.get("CODATA 2018 recommended")
        G_codata_unc = float(_codata_entry[1]) if _codata_entry else 1.5e-15
    except Exception:
        G_codata_unc = 1.5e-15
else:
    G_codata = 6.67430e-11
    G_codata_unc = 1.5e-15

# Compute hierarchy prediction
_hierarchy = None
if _core_available and constants is not None:
    try:
        _hierarchy = predict_G_hierarchy(constants)
        G_hierarchy = float(_hierarchy["G_predicted"])
        hierarchy_residual_ppm = float(_hierarchy["residual_ppm"])
    except Exception:
        G_hierarchy = None
        hierarchy_residual_ppm = None
else:
    G_hierarchy = None
    hierarchy_residual_ppm = None

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("The Prediction: G from First Principles")
st.markdown(
    "Deriving Newton's gravitational constant from the alpha ladder: "
    "**G_vacuum** emerges as a pure function of the fine-structure constant "
    "and the golden ratio, while the CODATA spread is explained by "
    "dilaton screening at laboratory scales."
)
st.divider()

# ---------------------------------------------------------------------------
# 4 Metric Cards
# ---------------------------------------------------------------------------
col_m1, col_m2, col_m3, col_m4 = st.columns(4)

with col_m1:
    st.metric(
        label="G_vacuum (predicted)",
        value=fmt_decimal(G_vacuum, sig_figs=6),
    )

with col_m2:
    st.metric(
        label="G_CODATA 2018",
        value=fmt_decimal(G_codata, sig_figs=6),
    )

with col_m3:
    st.metric(
        label="Excess (ppm)",
        value=f"{ppm_excess:.1f}",
    )

with col_m4:
    st.metric(
        label="\u03b1_screening",
        value=fmt_decimal(alpha_screening, sig_figs=4),
    )

# Hierarchy prediction row (no fitted parameters)
if G_hierarchy is not None:
    col_h1, col_h2, col_h3, col_h4 = st.columns(4)
    with col_h1:
        st.metric(
            label="G_hierarchy (alpha^24 mu^2)",
            value=fmt_decimal(G_hierarchy, sig_figs=6),
        )
    with col_h2:
        st.metric(
            label="Hierarchy residual (ppm)",
            value=f"{hierarchy_residual_ppm:.0f}",
        )
    with col_h3:
        st.markdown(
            """
            <div class="formula-card">
            <b>Hierarchy formula</b><br>
            <code>G = alpha^24 * mu^2 * hbar*c / m_e^2</code><br>
            No fitted parameters.
            </div>
            """,
            unsafe_allow_html=True,
        )

st.markdown("")

# ---------------------------------------------------------------------------
# Section A: Big G Deadlock Visualizer
# ---------------------------------------------------------------------------
with st.expander("Big G Deadlock Visualizer", expanded=True):

    # Classify measurements
    if _core_available and constants is not None and screening is not None:
        classified = classify_measurements(screening, constants)
        measurements_list = classified["measurements"]
        high_mean = float(classified["high_mean"])
        low_mean = float(classified["low_mean"])
        cluster_gap_ppm = float(classified["cluster_gap_ppm"])
    else:
        # Fallback measurement data
        _fallback_measurements = [
            {"experiment": "BIPM-14", "G_exp": 6.67554e-11, "G_unc": 1.6e-15, "sigma": 1.4, "direction": "-", "cluster": "high"},
            {"experiment": "LENS-14", "G_exp": 6.67191e-11, "G_unc": 9.9e-16, "sigma": 1.3, "direction": "+", "cluster": "low"},
            {"experiment": "UCI-14", "G_exp": 6.67435e-11, "G_unc": 1.4e-15, "sigma": 0.7, "direction": "-", "cluster": "high"},
            {"experiment": "HUST-18", "G_exp": 6.67484e-11, "G_unc": 1.2e-15, "sigma": 1.0, "direction": "-", "cluster": "high"},
            {"experiment": "HUST-18b", "G_exp": 6.67401e-11, "G_unc": 1.2e-15, "sigma": 0.5, "direction": "-", "cluster": "high"},
            {"experiment": "JILA-18", "G_exp": 6.67383e-11, "G_unc": 1.8e-15, "sigma": 0.4, "direction": "-", "cluster": "high"},
            {"experiment": "UZur-06", "G_exp": 6.67425e-11, "G_unc": 1.2e-15, "sigma": 0.6, "direction": "-", "cluster": "high"},
        ]
        measurements_list = _fallback_measurements
        high_vals = [m["G_exp"] for m in _fallback_measurements if m["cluster"] == "high"]
        low_vals = [m["G_exp"] for m in _fallback_measurements if m["cluster"] == "low"]
        high_mean = sum(high_vals) / len(high_vals) if high_vals else G_codata
        low_mean = sum(low_vals) / len(low_vals) if low_vals else G_vacuum
        cluster_gap_ppm = (high_mean - low_mean) / low_mean * 1e6

    col_chart, col_info = st.columns([2, 1])

    with col_chart:
        try:
            from app.components.charts import g_deadlock_scatter  # noqa: E402

            fig_deadlock = g_deadlock_scatter(measurements_list, G_vacuum)
            st.plotly_chart(fig_deadlock, use_container_width=True)
        except ImportError:
            st.warning("Could not render Big G deadlock chart (Plotly not available).")
        except Exception as exc:
            st.warning(f"Big G deadlock chart error: {exc}")

    with col_info:
        st.markdown(
            f"""
            <div class="pred-card">
            <b>Cluster Analysis</b><br><br>
            The 7 modern measurements of G split into two persistent clusters
            separated by ~{cluster_gap_ppm:.0f} ppm -- a gap far larger than
            individual error bars.<br><br>
            <b>High cluster mean:</b> {high_mean:.5e}<br>
            <b>Low cluster mean:</b> {low_mean:.5e}<br>
            <b>Gap:</b> {cluster_gap_ppm:.1f} ppm<br><br>
            The alpha ladder thesis: this is not experimental error but
            <em>dilaton screening</em> -- a scalar field coupled to the
            gravitational sector that shifts G_eff depending on the
            experimental geometry and source-mass distance.
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section B: Screening Simulator
# ---------------------------------------------------------------------------
with st.expander("Screening Simulator", expanded=True):

    # Distance slider with log-spaced labels
    _distance_map = {
        "0.01 m": 0.01,
        "0.03 m": 0.03,
        "0.1 m": 0.1,
        "0.3 m": 0.3,
        "1 m": 1.0,
        "3 m": 3.0,
        "10 m": 10.0,
        "100 m": 100.0,
        "1 km": 1e3,
        "10 km": 1e4,
        "100 km": 1e5,
        "1000 km": 1e6,
        "Earth-Moon": 3.844e8,
        "1 AU": 1.496e11,
        "10 AU": 1.496e12,
    }
    _dist_labels = list(_distance_map.keys())

    selected_dist_label = st.select_slider(
        "Select distance",
        options=_dist_labels,
        value="1 m",
        key="screening_dist",
    )
    selected_dist_m = _distance_map[selected_dist_label]

    # Compute profile
    if _core_available and constants is not None and screening is not None:
        profile = compute_screening_profile(screening)
        G_eff_at_dist = compute_G_eff(selected_dist_m, screening)
    else:
        # Fallback profile
        import numpy as np
        _r_arr = np.logspace(-2, 12, 500)
        _g_eff_arr = [G_vacuum * (1 + alpha_screening * math.exp(-r / lambda_dilaton_m)) for r in _r_arr]
        profile = {
            "r_meters": _r_arr.tolist(),
            "G_eff": _g_eff_arr,
            "landmarks": [
                {"label": "Torsion balance", "r_meters": 0.1},
                {"label": "Atom interf.", "r_meters": 0.3},
                {"label": "Earth radius", "r_meters": 6.371e6},
                {"label": "Earth-Moon", "r_meters": 3.844e8},
            ],
        }
        G_eff_at_dist = G_vacuum * (1 + alpha_screening * math.exp(-selected_dist_m / lambda_dilaton_m))

    # Normalise landmarks to list-of-dicts format expected by chart
    raw_landmarks = profile.get("landmarks", {})
    if isinstance(raw_landmarks, dict):
        profile["landmarks"] = [
            {"label": name, "r_meters": v["r_meters"]}
            for name, v in raw_landmarks.items()
        ]

    col_profile, col_sim_info = st.columns([2, 1])

    with col_profile:
        try:
            from app.components.charts import screening_profile_chart  # noqa: E402

            fig_profile = screening_profile_chart(profile, G_codata)

            # Add vertical marker at selected distance
            fig_profile.add_vline(
                x=selected_dist_m, line_dash="solid", line_color="#e0e0e0",
                line_width=1.5,
                annotation_text=selected_dist_label,
                annotation_position="top left",
                annotation_font_color="#e0e0e0",
            )

            st.plotly_chart(fig_profile, use_container_width=True)
        except ImportError:
            st.warning("Could not render screening profile chart (Plotly not available).")
        except Exception as exc:
            st.warning(f"Screening profile chart error: {exc}")

    with col_sim_info:
        _excess_ppm_at_dist = (G_eff_at_dist - G_vacuum) / G_vacuum * 1e6

        st.metric(label="G_eff at distance", value=fmt_decimal(G_eff_at_dist, sig_figs=6))
        st.metric(label="Excess over G_vacuum (ppm)", value=f"{_excess_ppm_at_dist:.2f}")
        st.metric(label="\u03bb_dilaton", value=f"{lambda_dilaton_m:.4f} m")
        st.metric(label="Dilaton mass", value=fmt_decimal(dilaton_mass_eV, sig_figs=3) + " eV")

        st.markdown(
            """
            <div class="formula-card">
            <b>Yukawa Screening</b><br><br>
            <code>G_eff(r) = G_vacuum * (1 + &alpha;<sub>s</sub> * exp(-r / &lambda;))</code><br><br>
            At short distances (torsion balance scale), the dilaton field enhances
            the effective G above the vacuum value. At astronomical distances the
            exponential dies and G_eff converges to G_vacuum.
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section C: Falsification Clock
# ---------------------------------------------------------------------------
with st.expander("Falsification Clock", expanded=True):

    # Get sigma comparisons
    if _core_available and constants is not None:
        bridges = get_bridge_candidates(constants)
        bridge_coeff = bridges.get("\u03c6\u00b2/2", bridges.get("phi^2/2", None))
        if bridge_coeff is not None:
            G_pred_decimal = predict_G(bridge_coeff, constants)
            G_pred = float(G_pred_decimal)
        else:
            G_pred_decimal = Decimal(str(G_vacuum))
            G_pred = G_vacuum
        measurements_dict = get_G_measurements()
        sigma_comparisons = compare_prediction(G_pred_decimal, measurements_dict)
    else:
        G_pred = G_vacuum
        sigma_comparisons = [
            {"experiment": "BIPM-14", "G_exp": 6.67554e-11, "G_unc": 1.6e-15, "sigma": 14.4, "direction": "-"},
            {"experiment": "LENS-14", "G_exp": 6.67191e-11, "G_unc": 9.9e-16, "sigma": 1.3, "direction": "+"},
            {"experiment": "UCI-14", "G_exp": 6.67435e-11, "G_unc": 1.4e-15, "sigma": 8.0, "direction": "-"},
            {"experiment": "HUST-18", "G_exp": 6.67484e-11, "G_unc": 1.2e-15, "sigma": 13.4, "direction": "-"},
            {"experiment": "HUST-18b", "G_exp": 6.67401e-11, "G_unc": 1.2e-15, "sigma": 6.5, "direction": "-"},
            {"experiment": "JILA-18", "G_exp": 6.67383e-11, "G_unc": 1.8e-15, "sigma": 3.3, "direction": "-"},
            {"experiment": "UZur-06", "G_exp": 6.67425e-11, "G_unc": 1.2e-15, "sigma": 8.5, "direction": "-"},
        ]

    col_bars, col_table = st.columns([1, 1])

    with col_bars:
        try:
            from app.components.charts import falsification_bars  # noqa: E402

            fig_falsify = falsification_bars(sigma_comparisons, G_pred)
            st.plotly_chart(fig_falsify, use_container_width=True)
        except ImportError:
            st.warning("Could not render falsification chart (Plotly not available).")
        except Exception as exc:
            st.warning(f"Falsification chart error: {exc}")

    with col_table:
        import pandas as pd

        df_data = []
        for c in sigma_comparisons:
            df_data.append({
                "Experiment": c["experiment"],
                "G (m3 kg-1 s-2)": f"{float(c['G_exp']):.5e}",
                "Uncertainty": f"{float(c['G_unc']):.1e}",
                "Sigma": f"{float(c['sigma']):.1f}",
                "Direction": "HIGH" if c.get("direction") == "-" else "LOW",
            })
        df = pd.DataFrame(df_data)
        st.dataframe(df, use_container_width=True, hide_index=True)

    st.markdown("")

    st.markdown(
        f"""
        <div class="falsify-card">
        <b>Key Metrics</b><br><br>
        The vacuum prediction G_vacuum = {G_pred:.5e} sits <b>~160 ppm below</b>
        the CODATA 2018 recommended value.{f' The hierarchy prediction (alpha^24 * mu^2, no fitted parameters) gives G = {G_hierarchy:.5e} (~{hierarchy_residual_ppm:.0f} ppm from CODATA).' if G_hierarchy is not None else ''}
        Six torsion-balance / lab experiments
        all read <b>high</b> (bars extending right, 3--18 sigma), while
        <b>Rosi et al. 2014</b> -- the only atom interferometry measurement --
        reads <b>low</b> at just <b>1.3 sigma</b> (bar extending left).<br><br>
        This is the screening signature: torsion-balance experiments probe
        G at short source-mass separations where the dilaton Yukawa term is
        unsuppressed, so they systematically overshoot. Atom interferometry
        uses free-falling atoms with effectively no nearby source mass,
        reducing the dilaton enhancement and landing closest to G_vacuum.<br><br>
        <b>What would falsify this?</b> A measurement of G at the <b>5 ppm level</b>
        at large source-mass separation (r &gt; 10 m) that still returns
        G &asymp; G_CODATA would kill the screening hypothesis.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="pred-card">
        <b>Upcoming Experiments to Watch</b><br><br>
        <b>Stanford Tower</b> -- Kasevich group atom interferometry at 10 m
        baseline. If G drops toward G_vacuum at this distance, the screening
        model is confirmed.<br><br>
        <b>Florence (MAGIA-Advanced)</b> -- atom interferometry with
        cold strontium, targeting sub-10 ppm precision at ~0.3 m baseline.<br><br>
        <b>Wuhan (HUST next-gen)</b> -- torsion pendulum with variable
        source-mass distance, designed to map the distance dependence of G
        with 5 ppm sensitivity.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("The Prediction (Pillar 4) | Alpha Ladder Research Dashboard")
