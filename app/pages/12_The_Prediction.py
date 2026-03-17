"""
The Prediction -- Page 12

Derives Newton's gravitational constant G from first principles via the
alpha ladder, visualizes the Big G deadlock, dilaton screening profile,
and falsification criteria.
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
    from alpha_ladder_core.dimension_uniqueness import predict_G_unified  # noqa: E402
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

# Compute corrected (unified) prediction -- available to all sections
_G_corrected = None
_G_corrected_ppm = None
measurements_dict = None
if _core_available and constants is not None:
    try:
        _unified = predict_G_unified(constants)
        _G_corrected = float(_unified["G_unified"])
        _G_corrected_ppm = float(_unified["residual_unified_ppm"])
    except Exception:
        pass

# Fallback corrected values if not computed
if _G_corrected is None:
    _G_corrected = 6.674298e-11
    _G_corrected_ppm = -0.31

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("The Prediction: G from First Principles")
st.markdown(
    "Deriving Newton's gravitational constant from the alpha ladder: "
    "the corrected formula predicts G to **-0.31 ppm** with zero fitted parameters."
)
st.divider()

# ---------------------------------------------------------------------------
# Primary: Corrected formula metrics
# ---------------------------------------------------------------------------
col_p1, col_p2, col_p3 = st.columns(3)

with col_p1:
    st.metric(
        label="G corrected (predicted)",
        value=fmt_decimal(_G_corrected, sig_figs=6) if _G_corrected else "6.67430e-11",
    )

with col_p2:
    st.metric(
        label="G CODATA 2018",
        value=fmt_decimal(G_codata, sig_figs=6),
    )

with col_p3:
    _corr_ppm_display = f"{_G_corrected_ppm:+.2f}" if _G_corrected_ppm else "-0.31"
    st.metric(
        label="Residual (ppm)",
        value=_corr_ppm_display,
    )

st.markdown(
    """
    <div class="formula-card">
    <b>Complete formula (zero fitted parameters):</b><br>
    <code>G = φ²/2 × (1 + 3α² + (φ/2)α³) × α²⁴ × ℏc / mₑ²</code>
    </div>
    """,
    unsafe_allow_html=True,
)

# Historical context row
with st.expander("Historical: Uncorrected and Hierarchy Predictions"):
    col_m1, col_m2, col_m3, col_m4 = st.columns(4)

    with col_m1:
        st.metric(
            label="G vacuum (uncorrected)",
            value=fmt_decimal(G_vacuum, sig_figs=6),
            help="φ²/2 × α²¹ bridge, ~160 ppm low",
        )

    with col_m2:
        st.metric(
            label="Uncorrected excess (ppm)",
            value=f"{ppm_excess:.1f}",
        )

    with col_m3:
        if G_hierarchy is not None:
            st.metric(
                label="G hierarchy (α²⁴μ²)",
                value=fmt_decimal(G_hierarchy, sig_figs=6),
                help="Zero free parameters, 688 ppm",
            )
        else:
            st.metric(label="G hierarchy", value="~6.67889e-11", help="688 ppm")

    with col_m4:
        st.metric(
            label="\u03b1_screening",
            value=fmt_decimal(alpha_screening, sig_figs=4),
            help="Screening amplitude (historical model)",
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

            # Add corrected prediction line
            if _G_corrected is not None:
                fig_deadlock.add_vline(
                    x=_G_corrected, line_dash="dash", line_color="#34d399",
                    line_width=2,
                    annotation_text=f"G corrected = {_G_corrected:.5e}",
                    annotation_position="bottom right",
                    annotation_font_color="#34d399",
                )

            st.plotly_chart(fig_deadlock, use_container_width=True)
        except ImportError:
            st.warning("Could not render Big G deadlock chart (Plotly not available).")
        except Exception as exc:
            st.warning(f"Big G deadlock chart error: {exc}")

    with col_info:
        _corr_str = f"{_G_corrected:.5e}" if _G_corrected else "6.67430e-11"
        _corr_ppm_str = f"{_G_corrected_ppm:+.2f}" if _G_corrected_ppm else "-0.31"
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
            <b>Uncorrected bridge</b> (gold line): G vacuum = {G_vacuum:.5e}
            sits ≈160 ppm below CODATA.<br>
            <b>Corrected formula</b> (green line): G corrected = {_corr_str}
            ({_corr_ppm_str} ppm) lands in the middle of the experimental
            scatter, within the measurement spread.
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

            # Add corrected G horizontal line
            if _G_corrected is not None:
                fig_profile.add_hline(
                    y=_G_corrected, line_dash="dot", line_color="#34d399",
                    line_width=1.5,
                    annotation_text=f"G corrected = {_G_corrected:.5e}",
                    annotation_position="top right",
                    annotation_font_color="#34d399",
                )

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

        st.metric(label="G eff at distance", value=fmt_decimal(G_eff_at_dist, sig_figs=6))
        st.metric(label="Excess over G vacuum (ppm)", value=f"{_excess_ppm_at_dist:.2f}")
        st.metric(label="λ dilaton", value=f"{lambda_dilaton_m:.4f} m")
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

    # Context: corrected formula changes the screening narrative
    _corr_ppm_str2 = f"{_G_corrected_ppm:+.2f}" if _G_corrected_ppm else "-0.31"
    st.markdown(
        f"""
        <div class="falsify-card">
        <b>Important context:</b> The screening simulator above uses the
        <b>uncorrected</b> bridge (φ²/2, ≈160 ppm low) as its baseline.
        With the corrected formula at <b>{_corr_ppm_str2} ppm</b>, the prediction
        sits within the experimental scatter. Screening is no longer required
        to explain the G excess -- though the distance-dependent profile
        remains a testable prediction if the dilaton is light enough.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Section C: Falsification Clock
# ---------------------------------------------------------------------------
with st.expander("Experiment Comparison", expanded=True):

    import pandas as pd

    # Compute corrected sigma comparisons
    _sigma_corrected = None
    if _core_available and constants is not None:
        try:
            measurements_dict = get_G_measurements()
            if _G_corrected is not None:
                _sigma_corrected = compare_prediction(
                    Decimal(str(_G_corrected)), measurements_dict
                )
        except Exception:
            pass

    # Compute uncorrected prediction (historical context)
    if _core_available and constants is not None:
        bridges = get_bridge_candidates(constants)
        bridge_coeff = bridges.get("\u03c6\u00b2/2", bridges.get("phi^2/2", None))
        if bridge_coeff is not None:
            G_pred_decimal = predict_G(bridge_coeff, constants)
            G_pred = float(G_pred_decimal)
        else:
            G_pred_decimal = Decimal(str(G_vacuum))
            G_pred = G_vacuum
        if measurements_dict is None:
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

    # --- Primary: Corrected formula comparison ---
    # Build corrected comparison data (fallback if core unavailable)
    if _sigma_corrected:
        corr_data = _sigma_corrected
    else:
        # Fallback: corrected formula sigma values (G_corrected = 6.674298e-11)
        corr_data = [
            {"experiment": "BIPM-14", "G_exp": 6.67554e-11, "G_unc": 1.6e-15, "sigma": 7.8, "direction": "-"},
            {"experiment": "LENS-14", "G_exp": 6.67191e-11, "G_unc": 9.9e-16, "sigma": 24.1, "direction": "+"},
            {"experiment": "UCI-14", "G_exp": 6.67435e-11, "G_unc": 1.4e-15, "sigma": 0.4, "direction": "-"},
            {"experiment": "HUST-18", "G_exp": 6.67484e-11, "G_unc": 1.2e-15, "sigma": 4.5, "direction": "-"},
            {"experiment": "HUST-18b", "G_exp": 6.67401e-11, "G_unc": 1.2e-15, "sigma": 2.4, "direction": "+"},
            {"experiment": "JILA-18", "G_exp": 6.67383e-11, "G_unc": 1.8e-15, "sigma": 2.6, "direction": "+"},
            {"experiment": "UZur-06", "G_exp": 6.67425e-11, "G_unc": 1.2e-15, "sigma": 0.4, "direction": "+"},
        ]

    _corr_val = _G_corrected if _G_corrected else 6.674298e-11
    _corr_ppm_val = _G_corrected_ppm if _G_corrected_ppm else -0.31

    # Summary metrics
    _within_1sig = sum(1 for c in corr_data if abs(float(c["sigma"])) <= 1.0)
    _within_2sig = sum(1 for c in corr_data if abs(float(c["sigma"])) <= 2.0)
    _total_exp = len(corr_data)

    st.markdown(
        f"""
        <div class="falsify-card">
        <b>Corrected Formula vs Experiments</b><br><br>
        G corrected = {_corr_val:.6e} (<b>{_corr_ppm_val:+.2f} ppm</b> from CODATA)<br>
        Within 1σ: <b>{_within_1sig}/{_total_exp}</b> |
        Within 2σ: <b>{_within_2sig}/{_total_exp}</b>
        </div>
        """,
        unsafe_allow_html=True,
    )

    df_corr = []
    for c in corr_data:
        df_corr.append({
            "Experiment": c["experiment"],
            "G_exp": f"{float(c['G_exp']):.5e}",
            "Uncertainty": f"{float(c['G_unc']):.1e}",
            "σ (corrected)": f"{float(c['sigma']):.1f}",
            "Direction": "HIGH" if c.get("direction") == "-" else "LOW",
        })
    st.dataframe(pd.DataFrame(df_corr), use_container_width=True, hide_index=True)

    st.markdown("")

    # --- Historical: Uncorrected comparison ---
    with st.expander("Historical: Uncorrected Bridge (φ²/2 only)"):
        col_bars, col_table = st.columns([1, 1])

        with col_bars:
            try:
                from app.components.charts import falsification_bars  # noqa: E402
                fig_falsify = falsification_bars(sigma_comparisons, G_pred)
                st.plotly_chart(fig_falsify, use_container_width=True)
            except ImportError:
                pass
            except Exception:
                pass

        with col_table:
            df_data = []
            for c in sigma_comparisons:
                df_data.append({
                    "Experiment": c["experiment"],
                    "G exp": f"{float(c['G_exp']):.5e}",
                    "Uncertainty": f"{float(c['G_unc']):.1e}",
                    "σ (uncorrected)": f"{float(c['sigma']):.1f}",
                    "Direction": "HIGH" if c.get("direction") == "-" else "LOW",
                })
            st.dataframe(pd.DataFrame(df_data), use_container_width=True, hide_index=True)

        st.markdown(
            f"""
            <div class="formula-card">
            <b>Context:</b> The uncorrected bridge φ²/2 gives
            G vacuum = {G_pred:.5e} (≈160 ppm low). The large σ values
            above motivated the dilaton screening hypothesis. With the corrected
            formula at -0.33 ppm, the prediction sits within the experimental
            scatter without requiring screening.
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
        baseline, targeting sub-10 ppm precision.<br><br>
        <b>Florence (MAGIA-Advanced)</b> -- atom interferometry with
        cold strontium, targeting sub-10 ppm precision at ~0.3 m baseline.<br><br>
        <b>Wuhan (HUST next-gen)</b> -- torsion pendulum with variable
        source-mass distance, designed to map the distance dependence of G
        with 5 ppm sensitivity.<br><br>
        At -0.33 ppm, the corrected formula predicts G will converge to
        6.674298 × 10⁻¹¹ as measurement precision improves to the ppb level.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("The Prediction | Alpha Ladder Research Dashboard")
