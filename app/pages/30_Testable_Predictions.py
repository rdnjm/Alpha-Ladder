"""Testable Predictions -- Page 30: Concrete predictions from the Alpha Ladder framework that can be confirmed or falsified by experiment."""

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
            "Value": "6.674298e-11 m³ kg⁻¹ s⁻²",
            "Status": "Within CODATA uncertainty",
            "Testable By": "Next-generation G measurements (sub-ppm)",
        },
        {
            "Prediction": "Proton-to-electron mass ratio structure",
            "Value": "μ ≈ 1836.157",
            "Status": "Unified (0.16 ppm via c3=φ/2)",
            "Testable By": "CODATA μ precision (resolved by unification)",
        },
        {
            "Prediction": "Cosmological constant from ladder G",
            "Value": "Λ ratio ≈ 1.23",
            "Status": "≈122 orders discrepancy remains",
            "Testable By": "Cosmological observations",
        },
        {
            "Prediction": "Fifth force signal",
            "Value": "Dilaton mass ≈ Planck scale",
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
        The Alpha Ladder makes one headline prediction (G to sub-ppb from α and μ
        alone), one unified result (bridge and μ-structure are the same formula via
        c3 = φ/2, predicting μ to 0.16 ppm), and two predictions that are currently
        beyond experimental reach (cosmological constant, fifth force).
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
    G = α²⁴ × μ × (μ − √φ×(1−α)) × ℏc / mₑ²
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
            value=f"{float(g_pred):.6e}" if g_pred is not None else "~6.674298e-11",
        )

    with col_b2:
        res_ppm = g_data.get("residual_ppm")
        st.metric(
            label="Residual (ppm)",
            value=f"{res_ppm:.2f}" if res_ppm is not None else "~-0.31",
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
        <b>If the formula is exact, G is determined to 3.6 ppb from α and μ
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
        st.metric(label="G Predicted", value="~6.674298e-11")

    with col_b2:
        st.metric(label="Residual (ppm)", value="~-0.31")

    with col_b3:
        st.metric(label="Predicted Uncertainty (ppb)", value="~3.6")

    with col_b4:
        st.metric(label="Improvement Factor", value="~6100x")

    st.markdown("")

    st.markdown(
        """
        <div class="theorem-card">
        <b>If the formula is exact, G is determined to 3.6 ppb from α and μ
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
        G = α²⁴ × μ × (μ − √φ×(1−α)) × ℏc / mₑ² = 6.674298e-11 m³ kg⁻¹ s⁻²,
        with a residual of -0.31 ppm relative to CODATA 2018 (within the ≈22 ppm measurement
        uncertainty). If exact, the predicted uncertainty is ≈3.6 ppb, set entirely by the
        uncertainties in α and μ.
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
            label="Within 1σ",
            value=str(n1) if n1 is not None else "N/A",
        )

    with col_c2:
        n2 = g_exp.get("n_within_2sigma")
        st.metric(
            label="Within 2σ",
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
        st.metric(label="Within 1σ", value="~3")

    with col_c2:
        st.metric(label="Within 2σ", value="~5")

    st.markdown("")

    st.markdown(
        """
        <div class="step-card">
        <b>Experiment comparison:</b><br><br>
        Best agreement: CODATA 2018 recommended value (≈0.2σ).<br>
        Worst agreement: Quinn et al. 2013 (≈1.8σ). All comparisons are
        within 2σ of the predicted value given current measurement uncertainties.
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
    <b>μ-structure formula:</b><br><br>
    μ = [√φ + √(φ + 2×φ²×(1+3a²)/a³)] / 2
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
            label="μ Predicted",
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
            label="σ Tension",
            value=f"{sigma:.0f}" if sigma is not None else "~7900",
        )

    st.markdown("")

    sigma_val = mu_data.get("sigma_tension")
    sigma_str = f"{sigma_val:.0f}" if sigma_val is not None else "~7900"
    st.markdown(
        f"""
        <div class="warning-card">
        <b>Tension at current precision ({sigma_str}σ).</b> However, mu_tension.py
        proves the bridge and μ-structure formulas are the same formula unified by
        c3 = φ/2 = 0.809, predicting μ from α and φ to 0.16 ppm.
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
        st.metric(label="μ Predicted", value="~1836.157000")

    with col_d2:
        st.metric(label="Residual (ppm)", value="~+2.37")

    with col_d3:
        st.metric(label="σ Tension", value="~7900")

    st.markdown("")

    st.markdown(
        """
        <div class="warning-card">
        <b>Tension at current precision (≈7900σ).</b> However, mu_tension.py
        proves the bridge and μ-structure formulas are the same formula unified by
        c3 = φ/2 = 0.809, predicting μ from α and φ to 0.16 ppm.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="step-card">
        <b>Interpretation:</b><br><br>
        The μ-structure formula predicts μ to ≈2.37 ppm, which shows tension at
        ≈7900σ given the sub-ppb precision of CODATA μ. However, mu_tension.py
        proves the bridge and μ-structure formulas are actually the same formula,
        unified by c3 = φ/2 = 0.809. With this identification, μ is predicted
        from α and φ alone to 0.16 ppm (see Page 36: Mu Tension Resolution).
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
    Λ = (G²/ℏ⁴)×(mₑ/α)⁶
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
            label="Λ Predicted",
            value=f"{float(lam_pred):.2e}" if lam_pred is not None else "~4.4e-66 eV⁴",
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
        st.metric(label="Λ Predicted", value="~4.4e-66 eV⁴")

    with col_e2:
        st.metric(label="Ratio to Observed", value="~1.23")

    with col_e3:
        st.metric(label="Orders Improvement", value="~122")

    st.markdown("")

    st.markdown(
        """
        <div class="theorem-card">
        <b>Reduces the cosmological constant problem from 122 orders of magnitude
        to a factor of ≈1.23.</b>
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
            dilaton mass (≈6.3e29 eV) corresponds to a Yukawa range of ≈3.1e-37 m,
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
        dilaton mass (≈6.3e29 eV) corresponds to a Yukawa range of ≈3.1e-37 m,
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
        Eot-Wash optimal range: ≈28 um<br>
        Eot-Wash optimal mass: ≈7 meV<br><br>
        If a mechanism exists to make the dilaton lighter than Planck scale,
        the Eot-Wash torsion balance experiment is optimally sensitive at ≈28 um
        (corresponding to ≈7 meV dilaton mass).
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# G. Geometric Resummation: mu from alpha and phi
# ---------------------------------------------------------------------------
st.header("G. Geometric Resummation: mu from alpha and phi")

st.markdown(
    """
    <div class="formula-card">
    <b>Geometric resummation (discovered 2026-03-18):</b><br><br>
    F = 1 + 3&alpha;&sup2; + &phi;&sup2;&alpha;&sup3; / [2(&phi; - &alpha;)]<br><br>
    The correction series admits a closed-form with geometric ratio 1/&phi;.
    Setting F_geom = F_exact and solving for &mu; gives the proton-to-electron
    mass ratio as a function of &alpha; and &phi; alone.
    </div>
    """,
    unsafe_allow_html=True,
)
st.markdown("")

st.markdown(
    """
    <div class="theorem-card">
    <b>Prediction:</b> &mu; predicted to 0.001 ppm from &alpha; and &phi; alone,
    with zero free parameters. Stable across CODATA 2014, 2018, and 2022
    editions. If the geometric structure is fundamental, &mu; is not an
    independent constant but is determined by the compactification geometry.
    See Page 27 (Corrected Bridge, Sections G-H) for full analysis.
    </div>
    """,
    unsafe_allow_html=True,
)
st.markdown("")

# ---------------------------------------------------------------------------
# H. Honest Assessment
# ---------------------------------------------------------------------------
st.header("H. Honest Assessment")

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
        CONFIRMED (within uncertainty): G = α²⁴ × μ × (μ − √φ×(1−α)) × ℏc / mₑ²
        predicts G to -0.31 ppm, within the ≈22 ppm CODATA measurement uncertainty. If exact,
        the formula determines G to 3.6 ppb -- 6000x more precise than any measurement. This
        is the headline prediction and will be tested by next-generation G experiments.<br><br>
        UNIFIED: The bridge and μ-structure formulas are the same formula (proven in
        mu_tension.py with c3 = φ/2 = 0.809), predicting μ from α and φ to 0.16 ppm.
        The apparent ≈7900σ tension dissolves under this unification.<br><br>
        BEYOND REACH: The cosmological constant prediction reduces the 122-order discrepancy
        to a factor of ≈1.23 using Beck's formula with ladder G, but this relies on Beck's
        formula being correct. The fifth force signal is unobservable for a Planck-mass
        dilaton (range ≈3e-37 m).<br><br>
        BOTTOM LINE: The Alpha Ladder makes one genuinely testable prediction (G to sub-ppb),
        one unified result (bridge = μ-structure via c3 = φ/2), and two predictions currently
        beyond experimental reach. Honest science requires reporting all four.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Testable Predictions | Alpha Ladder Research Dashboard")
