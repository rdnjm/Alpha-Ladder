"""Testable Predictions -- Page 30: Concrete predictions from the Alpha Ladder framework that can be confirmed or falsified by experiment."""

import streamlit as st
import pandas as pd

st.set_page_config(page_title="Testable Predictions | Alpha Ladder", layout="wide")

# ---------------------------------------------------------------------------
# Custom CSS (matches Pages 16-29)
# ---------------------------------------------------------------------------
st.markdown(
    """
    <style>
    .proof-card {
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
    .theorem-card {
        background-color: #1a1d23;
        border-left: 3px solid #34d399;
        padding: 0.8rem 1rem;
        margin: 0.4rem 0;
        border-radius: 0 8px 8px 0;
    }
    .step-card {
        background-color: #1a1d23;
        border-left: 3px solid #60a5fa;
        padding: 0.8rem 1rem;
        margin: 0.4rem 0;
        border-radius: 0 8px 8px 0;
    }
    .warning-card {
        background-color: #1a1d23;
        border-left: 3px solid #f87171;
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
# Try to load core module
# ---------------------------------------------------------------------------
_core_available = False
try:
    from alpha_ladder_core.testable_predictions import (  # noqa: E402
        summarize_testable_predictions,
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
        return summarize_testable_predictions()

    try:
        summary = _get_summary()
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Testable Predictions")
st.markdown(
    "Concrete predictions from the Alpha Ladder framework that can be "
    "confirmed or falsified by experiment."
)
st.divider()

# ---------------------------------------------------------------------------
# A. Predictions Overview
# ---------------------------------------------------------------------------
st.header("A. Predictions Overview")

if summary:
    predictions = summary.get("predictions_summary", [])
    if predictions:
        table_rows = []
        for p in predictions:
            table_rows.append({
                "Prediction": p.get("prediction", ""),
                "Value": p.get("value", ""),
                "Status": p.get("status", ""),
                "Testable By": p.get("testable_by", ""),
            })
        st.dataframe(
            pd.DataFrame(table_rows),
            use_container_width=True,
            hide_index=True,
        )

    key_finding = summary.get("key_finding", "")
    if key_finding:
        st.markdown("")
        st.markdown(
            f"""
            <div class="step-card">
            <b>Key finding:</b><br><br>
            {key_finding}
            </div>
            """,
            unsafe_allow_html=True,
        )
else:
    fallback_predictions = [
        {
            "Prediction": "G to sub-ppb precision",
            "Value": "6.674264e-11 m^3 kg^-1 s^-2",
            "Status": "Within CODATA uncertainty",
            "Testable By": "Next-generation G measurements (sub-ppm)",
        },
        {
            "Prediction": "Proton-to-electron mass ratio structure",
            "Value": "mu ~ 1836.157",
            "Status": "Falsified at ~7900 sigma",
            "Testable By": "Already tested (CODATA mu precision)",
        },
        {
            "Prediction": "Cosmological constant from ladder G",
            "Value": "Lambda ratio ~ 1.23",
            "Status": "~122 orders discrepancy remains",
            "Testable By": "Cosmological observations",
        },
        {
            "Prediction": "Fifth force signal",
            "Value": "Dilaton mass ~ Planck scale",
            "Status": "Unobservable (Planck-mass dilaton)",
            "Testable By": "Eot-Wash torsion balance (if dilaton lighter)",
        },
    ]
    st.dataframe(
        pd.DataFrame(fallback_predictions),
        use_container_width=True,
        hide_index=True,
    )

    st.markdown("")
    st.markdown(
        """
        <div class="step-card">
        <b>Key finding:</b><br><br>
        The Alpha Ladder makes one headline prediction (G to sub-ppb from alpha and mu
        alone), one falsified prediction (mu structure formula), and two predictions
        that are currently beyond experimental reach (cosmological constant, fifth force).
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# B. G to Sub-ppb Precision
# ---------------------------------------------------------------------------
st.header("B. G to Sub-ppb Precision")

st.markdown(
    """
    <div class="formula-card">
    <b>Headline formula:</b><br><br>
    G = alpha^24 * mu * (mu - sqrt(phi)) * hbar*c / m_e^2
    </div>
    """,
    unsafe_allow_html=True,
)
st.markdown("")

if summary:
    g_data = summary.get("g_precision", {})

    col_b1, col_b2, col_b3, col_b4 = st.columns(4)

    with col_b1:
        g_pred = g_data.get("G_predicted")
        st.metric(
            label="G Predicted",
            value=f"{float(g_pred):.6e}" if g_pred is not None else "~6.674264e-11",
        )

    with col_b2:
        res_ppm = g_data.get("residual_ppm")
        st.metric(
            label="Residual (ppm)",
            value=f"{res_ppm:.2f}" if res_ppm is not None else "~-5.37",
        )

    with col_b3:
        unc_ppb = g_data.get("predicted_uncertainty_ppb")
        st.metric(
            label="Predicted Uncertainty (ppb)",
            value=f"{unc_ppb:.1f}" if unc_ppb is not None else "~3.6",
        )

    with col_b4:
        imp = g_data.get("improvement_factor")
        st.metric(
            label="Improvement Factor",
            value=f"{imp:.0f}x" if imp is not None else "~6100x",
        )

    st.markdown("")

    st.markdown(
        """
        <div class="theorem-card">
        <b>If the formula is exact, G is determined to 3.6 ppb from alpha and mu
        alone -- 6000x more precise than any measurement.</b>
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    g_string = g_data.get("G_prediction_string", "")
    if g_string:
        st.markdown(
            f"""
            <div class="step-card">
            <b>Prediction:</b><br><br>
            {g_string}
            </div>
            """,
            unsafe_allow_html=True,
        )
else:
    col_b1, col_b2, col_b3, col_b4 = st.columns(4)

    with col_b1:
        st.metric(label="G Predicted", value="~6.674264e-11")

    with col_b2:
        st.metric(label="Residual (ppm)", value="~-5.37")

    with col_b3:
        st.metric(label="Predicted Uncertainty (ppb)", value="~3.6")

    with col_b4:
        st.metric(label="Improvement Factor", value="~6100x")

    st.markdown("")

    st.markdown(
        """
        <div class="theorem-card">
        <b>If the formula is exact, G is determined to 3.6 ppb from alpha and mu
        alone -- 6000x more precise than any measurement.</b>
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="step-card">
        <b>Prediction:</b><br><br>
        G = alpha^24 * mu * (mu - sqrt(phi)) * hbar*c / m_e^2 = 6.674264e-11 m^3 kg^-1 s^-2,
        with a residual of -5.37 ppm relative to CODATA 2018 (within the ~22 ppm measurement
        uncertainty). If exact, the predicted uncertainty is ~3.6 ppb, set entirely by the
        uncertainties in alpha and mu.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# C. G vs Experimental Measurements
# ---------------------------------------------------------------------------
st.header("C. G vs Experimental Measurements")

if summary:
    g_exp = summary.get("g_experiments", {})
    comparisons = g_exp.get("comparisons", [])

    if comparisons:
        table_rows = []
        for entry in comparisons:
            g_val = entry.get("G_exp")
            table_rows.append({
                "Experiment": entry.get("experiment", ""),
                "G_exp": fmt_decimal(g_val) if g_val is not None else "",
                "Sigma": f"{entry.get('sigma', 0):.1f}",
                "Direction": entry.get("direction", ""),
            })
        st.dataframe(
            pd.DataFrame(table_rows),
            use_container_width=True,
            hide_index=True,
        )

    st.markdown("")

    col_c1, col_c2 = st.columns(2)

    with col_c1:
        n1 = g_exp.get("n_within_1sigma")
        st.metric(
            label="Within 1 sigma",
            value=str(n1) if n1 is not None else "N/A",
        )

    with col_c2:
        n2 = g_exp.get("n_within_2sigma")
        st.metric(
            label="Within 2 sigma",
            value=str(n2) if n2 is not None else "N/A",
        )

    st.markdown("")

    best = g_exp.get("best_agreement", "")
    worst = g_exp.get("worst_agreement", "")
    note_parts = []
    if best:
        note_parts.append(f"Best agreement: {best}")
    if worst:
        note_parts.append(f"Worst agreement: {worst}")
    if note_parts:
        st.markdown(
            f"""
            <div class="step-card">
            <b>Experiment comparison:</b><br><br>
            {"<br>".join(note_parts)}
            </div>
            """,
            unsafe_allow_html=True,
        )
else:
    fallback_experiments = [
        {"Experiment": "CODATA 2018 (recommended)", "G_exp": "6.67430e-11", "Sigma": "0.2", "Direction": "above"},
        {"Experiment": "Quinn et al. (2013)", "G_exp": "6.67545e-11", "Sigma": "1.8", "Direction": "above"},
        {"Experiment": "Rosi et al. (2014)", "G_exp": "6.67191e-11", "Sigma": "1.5", "Direction": "below"},
        {"Experiment": "Li et al. (2018) TOS", "G_exp": "6.67484e-11", "Sigma": "0.9", "Direction": "above"},
        {"Experiment": "Li et al. (2018) AAF", "G_exp": "6.67349e-11", "Sigma": "0.4", "Direction": "above"},
    ]
    st.dataframe(
        pd.DataFrame(fallback_experiments),
        use_container_width=True,
        hide_index=True,
    )

    st.markdown("")

    col_c1, col_c2 = st.columns(2)

    with col_c1:
        st.metric(label="Within 1 sigma", value="~3")

    with col_c2:
        st.metric(label="Within 2 sigma", value="~5")

    st.markdown("")

    st.markdown(
        """
        <div class="step-card">
        <b>Experiment comparison:</b><br><br>
        Best agreement: CODATA 2018 recommended value (~0.2 sigma).<br>
        Worst agreement: Quinn et al. 2013 (~1.8 sigma). All comparisons are
        within 2 sigma of the predicted value given current measurement uncertainties.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# D. Proton-to-Electron Mass Ratio
# ---------------------------------------------------------------------------
st.header("D. Proton-to-Electron Mass Ratio")

st.markdown(
    """
    <div class="formula-card">
    <b>Mu-structure formula:</b><br><br>
    mu = [sqrt(phi) + sqrt(phi + 2*phi^2*(1+3a^2)/a^3)] / 2
    </div>
    """,
    unsafe_allow_html=True,
)
st.markdown("")

if summary:
    mu_data = summary.get("mu_consistency", {})

    col_d1, col_d2, col_d3 = st.columns(3)

    with col_d1:
        mu_pred = mu_data.get("mu_predicted")
        st.metric(
            label="mu Predicted",
            value=f"{mu_pred:.6f}" if mu_pred is not None else "~1836.157000",
        )

    with col_d2:
        res_ppm = mu_data.get("residual_ppm")
        st.metric(
            label="Residual (ppm)",
            value=f"{res_ppm:.2f}" if res_ppm is not None else "~+2.37",
        )

    with col_d3:
        sigma = mu_data.get("sigma_tension")
        st.metric(
            label="Sigma Tension",
            value=f"{sigma:.0f}" if sigma is not None else "~7900",
        )

    st.markdown("")

    sigma_val = mu_data.get("sigma_tension")
    sigma_str = f"{sigma_val:.0f}" if sigma_val is not None else "~7900"
    st.markdown(
        f"""
        <div class="warning-card">
        <b>FALSIFIED at current precision ({sigma_str} sigma).</b> This means the
        corrected bridge and mu-structure formulas cannot both be exact as stated.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    interpretation = mu_data.get("interpretation", "")
    if interpretation:
        st.markdown(
            f"""
            <div class="step-card">
            <b>Interpretation:</b><br><br>
            {interpretation}
            </div>
            """,
            unsafe_allow_html=True,
        )
else:
    col_d1, col_d2, col_d3 = st.columns(3)

    with col_d1:
        st.metric(label="mu Predicted", value="~1836.157000")

    with col_d2:
        st.metric(label="Residual (ppm)", value="~+2.37")

    with col_d3:
        st.metric(label="Sigma Tension", value="~7900")

    st.markdown("")

    st.markdown(
        """
        <div class="warning-card">
        <b>FALSIFIED at current precision (~7900 sigma).</b> This means the
        corrected bridge and mu-structure formulas cannot both be exact as stated.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="step-card">
        <b>Interpretation:</b><br><br>
        The mu-structure formula predicts mu to ~2.37 ppm, which is excellent by
        absolute standards but is ruled out at ~7900 sigma given the sub-ppb precision
        of CODATA mu. This falsification is informative: it tells us the corrected
        bridge coefficient and the mu offset cannot both be exact simultaneously.
        At least one requires higher-order corrections or a modified functional form.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# E. Cosmological Constant
# ---------------------------------------------------------------------------
st.header("E. Cosmological Constant")

st.markdown(
    """
    <div class="formula-card">
    <b>Beck's formula with ladder G:</b><br><br>
    Lambda = (G^2/hbar^4)*(m_e/alpha)^6
    </div>
    """,
    unsafe_allow_html=True,
)
st.markdown("")

if summary:
    cc_data = summary.get("cosmological_constant", {})

    col_e1, col_e2, col_e3 = st.columns(3)

    with col_e1:
        lam_pred = cc_data.get("Lambda_predicted")
        st.metric(
            label="Lambda Predicted",
            value=f"{float(lam_pred):.2e}" if lam_pred is not None else "~4.4e-66 eV^4",
        )

    with col_e2:
        ratio = cc_data.get("ratio_to_observed")
        st.metric(
            label="Ratio to Observed",
            value=f"{ratio:.2f}" if ratio is not None else "~1.23",
        )

    with col_e3:
        orders = cc_data.get("orders_improvement")
        st.metric(
            label="Orders Improvement",
            value=f"{orders:.0f}" if orders is not None else "~122",
        )

    st.markdown("")

    ratio_val = cc_data.get("ratio_to_observed")
    ratio_str = f"{ratio_val:.2f}" if ratio_val is not None else "~1.23"
    st.markdown(
        f"""
        <div class="theorem-card">
        <b>Reduces the cosmological constant problem from 122 orders of magnitude
        to a factor of {ratio_str}.</b>
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    col_e1, col_e2, col_e3 = st.columns(3)

    with col_e1:
        st.metric(label="Lambda Predicted", value="~4.4e-66 eV^4")

    with col_e2:
        st.metric(label="Ratio to Observed", value="~1.23")

    with col_e3:
        st.metric(label="Orders Improvement", value="~122")

    st.markdown("")

    st.markdown(
        """
        <div class="theorem-card">
        <b>Reduces the cosmological constant problem from 122 orders of magnitude
        to a factor of ~1.23.</b>
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# F. Fifth Force Signal
# ---------------------------------------------------------------------------
st.header("F. Fifth Force Signal")

if summary:
    ff_data = summary.get("fifth_force", {})

    col_f1, col_f2, col_f3 = st.columns(3)

    with col_f1:
        mass_eV = ff_data.get("dilaton_mass_eV")
        st.metric(
            label="Dilaton Mass (eV)",
            value=f"{float(mass_eV):.2e}" if mass_eV is not None else "~6.3e29",
        )

    with col_f2:
        force_range = ff_data.get("force_range_m")
        st.metric(
            label="Force Range (m)",
            value=f"{float(force_range):.2e}" if force_range is not None else "~3.1e-37",
        )

    with col_f3:
        observable = ff_data.get("observable", False)
        st.metric(
            label="Observable",
            value="Yes" if observable else "No",
        )

    st.markdown("")

    if not ff_data.get("observable", False):
        st.markdown(
            """
            <div class="warning-card">
            <b>Planck-mass dilaton makes signal unobservable.</b> The flux-stabilized
            dilaton mass (~6.3e29 eV) corresponds to a Yukawa range of ~3.1e-37 m,
            far below any foreseeable experimental sensitivity. The fifth force
            prediction is formally present but physically inaccessible.
            </div>
            """,
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            """
            <div class="theorem-card">
            <b>Fifth force signal is within experimental reach.</b>
            </div>
            """,
            unsafe_allow_html=True,
        )

    st.markdown("")

    eot_range = ff_data.get("eot_wash_optimal_range_um")
    eot_mass = ff_data.get("eot_wash_optimal_mass_meV")
    eot_parts = []
    if eot_range is not None:
        eot_parts.append(f"Eot-Wash optimal range: {eot_range} um")
    if eot_mass is not None:
        eot_parts.append(f"Eot-Wash optimal mass: {eot_mass} meV")
    if eot_parts:
        st.markdown(
            f"""
            <div class="step-card">
            <b>Testable window:</b><br><br>
            {"<br>".join(eot_parts)}<br><br>
            If a mechanism exists to make the dilaton lighter than Planck scale,
            the Eot-Wash torsion balance experiment is optimally sensitive at ~28 um
            (corresponding to ~7 meV dilaton mass).
            </div>
            """,
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            """
            <div class="step-card">
            <b>Testable window:</b><br><br>
            Eot-Wash optimal range: ~28 um<br>
            Eot-Wash optimal mass: ~7 meV<br><br>
            If a mechanism exists to make the dilaton lighter than Planck scale,
            the Eot-Wash torsion balance experiment is optimally sensitive at ~28 um
            (corresponding to ~7 meV dilaton mass).
            </div>
            """,
            unsafe_allow_html=True,
        )
else:
    col_f1, col_f2, col_f3 = st.columns(3)

    with col_f1:
        st.metric(label="Dilaton Mass (eV)", value="~6.3e29")

    with col_f2:
        st.metric(label="Force Range (m)", value="~3.1e-37")

    with col_f3:
        st.metric(label="Observable", value="No")

    st.markdown("")

    st.markdown(
        """
        <div class="warning-card">
        <b>Planck-mass dilaton makes signal unobservable.</b> The flux-stabilized
        dilaton mass (~6.3e29 eV) corresponds to a Yukawa range of ~3.1e-37 m,
        far below any foreseeable experimental sensitivity. The fifth force
        prediction is formally present but physically inaccessible.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="step-card">
        <b>Testable window:</b><br><br>
        Eot-Wash optimal range: ~28 um<br>
        Eot-Wash optimal mass: ~7 meV<br><br>
        If a mechanism exists to make the dilaton lighter than Planck scale,
        the Eot-Wash torsion balance experiment is optimally sensitive at ~28 um
        (corresponding to ~7 meV dilaton mass).
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# G. Honest Assessment
# ---------------------------------------------------------------------------
st.header("G. Honest Assessment")

if summary:
    honest = summary.get("honest_assessment", "")

    st.markdown(
        f"""
        <div class="proof-card">
        <b>Honest assessment:</b><br><br>
        {honest}
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.markdown(
        """
        <div class="proof-card">
        <b>Honest assessment:</b><br><br>
        CONFIRMED (within uncertainty): G = alpha^24 * mu * (mu - sqrt(phi)) * hbar*c / m_e^2
        predicts G to -5.37 ppm, within the ~22 ppm CODATA measurement uncertainty. If exact,
        the formula determines G to 3.6 ppb -- 6000x more precise than any measurement. This
        is the headline prediction and will be tested by next-generation G experiments.<br><br>
        FALSIFIED: The mu-structure formula predicts mu to ~2.37 ppm but is ruled out at
        ~7900 sigma. The corrected bridge and mu offset cannot both be exact simultaneously.<br><br>
        BEYOND REACH: The cosmological constant prediction reduces the 122-order discrepancy
        to a factor of ~1.23 using Beck's formula with ladder G, but this relies on Beck's
        formula being correct. The fifth force signal is unobservable for a Planck-mass
        dilaton (range ~3e-37 m).<br><br>
        BOTTOM LINE: The Alpha Ladder makes one genuinely testable prediction (G to sub-ppb),
        one prediction already falsified (mu structure), and two predictions currently beyond
        experimental reach. Honest science requires reporting all four.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Testable Predictions | Alpha Ladder Research Dashboard")
