"""
Fifth Force Predictions -- Page 21

The most important dashboard page: shows experimentalists exactly where
to look for the Alpha Ladder's falsifiable prediction.  The coupling
alpha = 0.618 is fixed by theory; the only freedom is a_0 (internal
radius), which sets the Yukawa range lambda.
"""

import streamlit as st
import pandas as pd


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
# Try to load core module
# ---------------------------------------------------------------------------
_core_available = False
try:
    from alpha_ladder_core.fifth_force_predictions import (  # noqa: E402
        summarize_fifth_force_predictions,
    )
    _core_available = True
except ImportError:
    pass

# ---------------------------------------------------------------------------
# Compute values (guarded, cached)
# ---------------------------------------------------------------------------
summary = None

if _core_available:
    @st.cache_data
    def _get_summary():
        return summarize_fifth_force_predictions()

    try:
        summary = _get_summary()
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Fifth Force Predictions")
st.markdown(
    "Where experimentalists should look for the Alpha Ladder's falsifiable "
    "prediction: a Yukawa fifth force with coupling α = 0.618 "
    "(fixed by theory) and range λ = a₀ (internal radius)."
)
st.divider()

# ---------------------------------------------------------------------------
# A. The Prediction
# ---------------------------------------------------------------------------
st.header("A. The Prediction")

if summary:
    st.info(summary["key_prediction"])

    survival_um = summary["survival_window_um"]
    optimal_exp = summary["optimal_experiment"]

    col_a1, col_a2, col_a3 = st.columns(3)

    with col_a1:
        st.metric(
            label="Coupling α",
            value="0.618",
            delta="Fixed by theory",
        )

    with col_a2:
        st.metric(
            label="Survival Window",
            value=f"{survival_um[0]:.0f} - {survival_um[1]:.0f} um",
        )

    with col_a3:
        st.metric(
            label="Best Experiment",
            value=optimal_exp,
        )
else:
    st.info(
        "The Alpha Ladder predicts a Yukawa fifth force with coupling "
        "α = 0.618 (fixed by theory). The range λ = a₀ is "
        "the only free parameter. The surviving window is approximately "
        "30 - 71 um."
    )
    col_a1, col_a2, col_a3 = st.columns(3)
    with col_a1:
        st.metric(label="Coupling α", value="0.618", delta="Fixed by theory")
    with col_a2:
        st.metric(label="Survival Window", value="30 - 71 um")
    with col_a3:
        st.metric(label="Best Experiment", value="Eot-Wash next-gen")

st.markdown("")

# ---------------------------------------------------------------------------
# B. Exclusion Map
# ---------------------------------------------------------------------------
st.header("B. Exclusion Map")

st.markdown(
    """
    <div class="step-card">
    <b>Interpretation:</b> The Alpha Ladder predicts a horizontal line
    at α = 0.618 in the standard exclusion plot. Experiments that
    reach below this line exclude the corresponding λ range.
    Any experiment with α_bound &lt; 0.618 at some λ rules out
    that value of a₀.
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown("")

if summary:
    exclusion_map = summary["exclusion_map"]
    experiments = exclusion_map["experiments"]

    table_rows = []
    for exp in experiments:
        lam_range = exp["lambda_range_m"]
        # Format lambda range for display
        if lam_range[0] < 1e-3:
            lam_str = f"{lam_range[0]*1e6:.1f} - {lam_range[1]*1e6:.0f} um"
        elif lam_range[0] < 1.0:
            lam_str = f"{lam_range[0]*1e3:.1f} - {lam_range[1]*1e3:.0f} mm"
        else:
            lam_str = f"{lam_range[0]:.2e} - {lam_range[1]:.2e} m"

        alpha_bound = exp["alpha_bound_at_best"]
        if alpha_bound > 1e6:
            alpha_str = f"{alpha_bound:.0e}"
        elif alpha_bound > 1.0:
            alpha_str = f"{alpha_bound:.1f}"
        else:
            alpha_str = fmt_decimal(alpha_bound, sig_figs=3)

        excludes = exp["excludes_alpha_ladder"]

        table_rows.append({
            "Experiment": exp["name"],
            "Lambda Range": lam_str,
            "Best Alpha Bound": alpha_str,
            "Excludes AL": "Yes" if excludes else "No",
            "Reference": exp["reference"],
        })

    st.dataframe(
        pd.DataFrame(table_rows),
        use_container_width=True,
        hide_index=True,
    )

    st.markdown("")

    # Exclusion plot chart
    try:
        from app.components.charts import exclusion_plot_chart  # noqa: E402
        fig_excl = exclusion_plot_chart(exclusion_map)
        st.plotly_chart(fig_excl, use_container_width=True)
    except ImportError:
        st.info("Chart function exclusion_plot_chart not yet available.")
    except Exception as exc:
        st.warning(f"Exclusion plot chart error: {exc}")
else:
    st.markdown(
        """
        <div class="formula-card">
        <b>Exclusion map (summary):</b><br><br>
        Eot-Wash 2006 excludes a₀ &gt; ~71 um. Lunar laser ranging and
        Cassini PPN exclude a₀ at solar-system scales. Casimir experiments
        do not reach α = 0.618 at any λ due to their much higher
        α bounds. The survival window is approximately 30 - 71 um.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# C. Eot-Wash Signal
# ---------------------------------------------------------------------------
st.header("C. Eot-Wash Signal")

if summary:
    eot_wash = summary["eot_wash_prediction"]

    predicted_signal = eot_wash["predicted_signal"]
    current_sensitivity = eot_wash["current_sensitivity"]
    ratio = eot_wash["signal_to_sensitivity_ratio"]

    col_c1, col_c2, col_c3 = st.columns(3)

    with col_c1:
        st.metric(
            label="Predicted Signal",
            value=f"{predicted_signal * 100:.1f}%",
        )

    with col_c2:
        st.metric(
            label="Current Sensitivity",
            value=f"{current_sensitivity * 100:.0f}%",
        )

    with col_c3:
        st.metric(
            label="Signal / Sensitivity",
            value=f"{ratio:.1f}x",
        )

    st.markdown("")

    st.success(eot_wash["signal_description"])

    st.markdown("")

    # Signal vs distance chart
    try:
        from alpha_ladder_core.fifth_force_predictions import (  # noqa: E402
            compute_signal_vs_distance,
        )
        from app.components.charts import signal_vs_distance_chart  # noqa: E402

        signal_data = compute_signal_vs_distance(lambda_m=30e-6)
        fig_signal = signal_vs_distance_chart(signal_data)
        st.plotly_chart(fig_signal, use_container_width=True)
    except ImportError:
        st.info("Signal vs distance chart not yet available.")
    except Exception as exc:
        st.warning(f"Signal vs distance chart error: {exc}")
else:
    st.markdown(
        """
        <div class="theorem-card">
        <b>Eot-Wash signal (a₀ = 30 um):</b><br><br>
        At a₀ = 30 um, the predicted Yukawa deviation is ≈11% at the
        52 um gap distance. Current Eot-Wash sensitivity is ≈1%.
        The signal is ≈11x the current sensitivity -- easily detectable
        if a₀ falls in this range.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# D. Casimir vs Yukawa
# ---------------------------------------------------------------------------
st.header("D. Casimir vs Yukawa")

st.markdown(
    """
    <div class="step-card">
    <b>Why Casimir experiments are not useful here:</b><br><br>
    At sub-micron separations, the electromagnetic Casimir force between
    conducting plates overwhelms any gravitational Yukawa signal by many
    orders of magnitude.  Torsion balance experiments (which use
    electrostatic shielding to eliminate the Casimir background) are the
    correct strategy.
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown("")

if summary:
    casimir = summary["casimir_overlap"]

    col_d1, col_d2, col_d3 = st.columns(3)

    with col_d1:
        st.metric(
            label="Casimir Pressure",
            value=fmt_decimal(casimir["casimir_pressure_Pa"], sig_figs=3) + " Pa",
        )

    with col_d2:
        st.metric(
            label="Yukawa Pressure",
            value=fmt_decimal(casimir["yukawa_pressure_Pa"], sig_figs=3) + " Pa",
        )

    with col_d3:
        st.metric(
            label="Yukawa / Casimir Ratio",
            value=fmt_decimal(casimir["ratio_yukawa_to_casimir"], sig_figs=3),
        )

    st.markdown("")

    st.warning(casimir["honest_assessment"])
else:
    st.markdown(
        """
        <div class="warning-card">
        <b>Casimir vs Yukawa (a₀ = 30 um):</b><br><br>
        The Casimir pressure dominates the Yukawa pressure by many orders
        of magnitude.  Casimir experiments are NOT the right approach for
        detecting this fifth force.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# E. Discovery Reach
# ---------------------------------------------------------------------------
st.header("E. Discovery Reach")

if summary:
    discovery = summary["discovery_reach"]
    reach_experiments = discovery["experiments"]

    table_rows = []
    for exp in reach_experiments:
        lam_range = exp.get("lambda_range_detectable")
        if lam_range:
            lam_str = f"{lam_range[0]*1e6:.1f} - {lam_range[1]*1e6:.0f} um"
        else:
            lam_str = "N/A"

        table_rows.append({
            "Experiment": exp["name"],
            "Gap (um)": f"{exp['gap_m']*1e6:.1f}",
            "Sensitivity": f"{exp['sensitivity']*100:.1f}%",
            "Signal at Optimal": f"{exp['signal_at_optimal']*100:.1f}%",
            "Detectable": "Yes" if exp["detectable_at_optimal"] else "No",
            "Lambda Range": lam_str,
            "Status": exp["status"],
        })

    st.dataframe(
        pd.DataFrame(table_rows),
        use_container_width=True,
        hide_index=True,
    )

    st.markdown("")

    # Discovery reach chart
    try:
        from app.components.charts import discovery_reach_chart  # noqa: E402
        fig_reach = discovery_reach_chart(discovery)
        st.plotly_chart(fig_reach, use_container_width=True)
    except ImportError:
        st.info("Discovery reach chart not yet available.")
    except Exception as exc:
        st.warning(f"Discovery reach chart error: {exc}")

    st.markdown("")

    st.info(discovery["honest_assessment"])
else:
    st.markdown(
        """
        <div class="theorem-card">
        <b>Discovery reach (summary):</b><br><br>
        The optimal signal for any experiment is α × 2/e = 45.5% when
        λ = gap distance. Eot-Wash torsion balance experiments are
        the most promising because they directly measure gravitational
        forces without Casimir backgrounds and operate at gap distances
        (30-52 um) overlapping the surviving a_0 window.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# F. Salam-Sezgin Update
# ---------------------------------------------------------------------------
st.header("F. Salam-Sezgin Update")

st.markdown(
    """
    <div class="warning-card">
    <b>UPDATE (Page 39):</b> The Salam-Sezgin mechanism stabilizes the
    S^2 radius with Lambda_6 = 14.2 (O(1) in Planck units). The dilaton
    acquires a Planck-scale mass m_phi ~ 2 x 10^29 eV, giving a Yukawa
    range of ~3 x 10^-37 m -- far below any foreseeable experiment.<br><br>
    With the Salam-Sezgin stabilization, the fifth force is
    <b>unobservable at all experimental scales</b>. The 3854x screening
    gap dissolves trivially (Planck-mass dilaton decouples identically).
    The analysis above (Eot-Wash window, survival region) applies only
    if an alternative stabilization mechanism produces a lighter dilaton.
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown("")

# ---------------------------------------------------------------------------
# G. The Bottom Line
# ---------------------------------------------------------------------------
st.header("G. The Bottom Line")

if summary:
    st.markdown(
        f"""
        <div class="proof-card">
        <b>Honest assessment:</b><br><br>
        {summary["honest_assessment"]}
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.markdown(
        """
        <div class="proof-card">
        <b>Honest assessment:</b><br><br>
        The fifth force prediction depends on the dilaton mass, which
        depends on the stabilization mechanism. With the Salam-Sezgin
        mechanism (Lambda_6 = 14.2), the dilaton is Planck-mass and the
        fifth force is unobservable. If a different mechanism produces a
        lighter dilaton, the Eot-Wash window [30, 71] um becomes relevant.
        The framework's primary testable prediction is G to -0.31 ppm,
        not the fifth force.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

st.markdown(
    """
    <div class="theorem-card">
    <b>Conclusion:</b><br><br>
    The fifth force coupling alpha = 0.618 is fixed by the framework, but
    the Yukawa range depends on the dilaton mass. The Salam-Sezgin
    stabilization gives a Planck-mass dilaton, making the fifth force
    unobservable. The framework's falsifiable prediction is G to sub-ppm
    precision, not a fifth force signal.
    </div>
    """,
    unsafe_allow_html=True,
)

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Fifth Force Predictions | Alpha Ladder Research Dashboard")
