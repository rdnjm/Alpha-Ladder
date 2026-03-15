"""
Hierarchy Derivation -- Page 26

The formula alpha_g = alpha^24 * mu^2 reproduces the gravitational coupling
to 688 ppm with zero fitted parameters. Every exponent comes from d=4, n=2, D=6.
This page explores 10 theoretical angles for why this formula holds.
"""

import streamlit as st
import pandas as pd


# ---------------------------------------------------------------------------
# Custom CSS (matches Pages 16-25)
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
    from alpha_ladder_core.hierarchy_derivation import (  # noqa: E402
        summarize_hierarchy_derivation,
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
        return summarize_hierarchy_derivation()

    try:
        summary = _get_summary()
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Hierarchy Derivation")
st.markdown(
    "The formula G = alpha^24 * m_p^2 * hbar*c / m_e^4 reproduces the "
    "gravitational constant to 688 ppm with zero fitted parameters. Every "
    "exponent comes from d=4, n=2, D=6."
)
st.divider()

# ---------------------------------------------------------------------------
# A. The Formula
# ---------------------------------------------------------------------------
st.header("A. The Formula")

if summary:
    basics = summary["formula_basics"]

    col_a1, col_a2, col_a3 = st.columns(3)

    with col_a1:
        st.metric(
            label="Exponent (d*D)",
            value=str(basics.get("exponent", 24)),
        )

    with col_a2:
        ppm_val = basics.get("residual_ppm")
        ppm_str = f"{ppm_val:.0f} ppm" if ppm_val is not None else "~688 ppm"
        st.metric(
            label="Residual",
            value=ppm_str,
        )

    with col_a3:
        mu_val = basics.get("mu")
        mu_str = f"{mu_val:.2f}" if mu_val is not None else "~1836"
        st.metric(
            label="mu = m_p / m_e",
            value=mu_str,
        )

    st.markdown("")

    formula_text = basics.get(
        "formula_latex",
        "alpha_g = alpha^{d*D} * mu^n = alpha^24 * mu^2"
    )
    st.markdown(
        f"""
        <div class="formula-card">
        <b>The hierarchy formula:</b><br><br>
        {formula_text}<br><br>
        With d=4 spacetime dimensions, n=2 extra dimensions, D=6 total dimensions.
        Zero fitted parameters. Residual ~688 ppm.
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    col_a1, col_a2, col_a3 = st.columns(3)

    with col_a1:
        st.metric(label="Exponent (d*D)", value="24")

    with col_a2:
        st.metric(label="Residual", value="~688 ppm")

    with col_a3:
        st.metric(label="mu = m_p / m_e", value="~1836")

    st.markdown("")

    st.markdown(
        """
        <div class="formula-card">
        <b>The hierarchy formula:</b><br><br>
        alpha_g = alpha^{d*D} * mu^n = alpha^24 * mu^2<br><br>
        With d=4 spacetime dimensions, n=2 extra dimensions, D=6 total dimensions.
        Zero fitted parameters. Residual ~688 ppm.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# B. Uniqueness -- Dimension Scan
# ---------------------------------------------------------------------------
st.header("B. Uniqueness -- Dimension Scan")

if summary:
    angle_6 = summary["angle_6"]

    # Dimension scan chart
    try:
        from app.components.charts import dimension_scan_chart  # noqa: E402
        fig = dimension_scan_chart(angle_6)
        st.plotly_chart(fig, use_container_width=True)
    except ImportError:
        st.info("Chart function dimension_scan_chart not yet available.")
    except Exception as exc:
        st.warning(f"Dimension scan chart error: {exc}")

    st.markdown("")

    # Scan results table
    scan_results = angle_6.get("scan_results", [])
    if scan_results:
        table_rows = []
        for r in scan_results:
            table_rows.append({
                "d": r.get("d", ""),
                "n": r.get("n", ""),
                "D = d + n": r.get("D", ""),
                "Exponent (d*D)": r.get("exponent", ""),
                "ppm": r.get("ppm", ""),
            })
        df = pd.DataFrame(table_rows)
        if "ppm" in df.columns:
            df = df.sort_values("ppm", ascending=True)
        st.dataframe(df, use_container_width=True, hide_index=True)

    st.markdown("")

    uniqueness = angle_6.get("assessment", "")
    st.markdown(
        f"""
        <div class="theorem-card">
        <b>Uniqueness:</b><br><br>
        {uniqueness if uniqueness else "(4,2) is the unique pair with sub-1000 ppm agreement. No other (d,n) comes close."}
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.markdown(
        """
        <div class="theorem-card">
        <b>Dimension scan:</b><br><br>
        Scanning all integer pairs (d, n) with d in [3..11] and n in [1..7]:<br>
        only (d=4, n=2) yields sub-1000 ppm agreement with G_observed.<br><br>
        The exponent d*D = 4*6 = 24 is uniquely selected by the data.
        No fitted parameters are involved -- the integers come from
        spacetime and compactification geometry alone.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# C. Standard KK -- Volume Suppression
# ---------------------------------------------------------------------------
st.header("C. Standard KK Volume Suppression")

if summary:
    angle_1 = summary["angle_1"]

    status = angle_1.get("status", "fails")
    reason = angle_1.get("reason", "")

    st.markdown(
        f"""
        <div class="warning-card">
        <b>Standard KK volume suppression: {status}</b><br><br>
        {reason if reason else (
            "In standard Kaluza-Klein, G_4 ~ G_D / V_n where V_n is the "
            "volume of the compact space. This gives G ~ (M_Pl^(2-D)) * R^n. "
            "The volume suppression produces n powers of the compactification "
            "radius, not d*D = 24 powers of alpha. The exponent structure "
            "does not match."
        )}
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    detail = angle_1.get("detail", "")
    if detail:
        st.markdown(
            f"""
            <div class="step-card">
            <b>Detail:</b><br><br>
            {detail}
            </div>
            """,
            unsafe_allow_html=True,
        )
else:
    st.markdown(
        """
        <div class="warning-card">
        <b>Standard KK volume suppression: fails</b><br><br>
        In standard Kaluza-Klein, G_4 ~ G_D / V_n where V_n is the
        volume of the compact space. This gives G ~ (M_Pl^{2-D}) * R^n.
        The volume suppression produces n powers of the compactification
        radius, not d*D = 24 powers of alpha. The exponent structure
        does not match.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# D. Metric Counting
# ---------------------------------------------------------------------------
st.header("D. Metric Counting")

if summary:
    angle_2 = summary["angle_2"]

    components = angle_2.get("components", {})
    graviton = components.get("graviton", 10)
    vectors = components.get("vectors", 8)
    scalars = components.get("scalars", 3)
    total = components.get("total", 21)
    target = components.get("target_exponent", 24)

    st.markdown(
        f"""
        <div class="step-card">
        <b>KK decomposition of the D=6 metric:</b><br><br>
        Symmetric metric g_MN in D=6 has D(D+1)/2 = 21 independent components.<br><br>
        KK decomposition:<br>
        - Graviton g_mu,nu: {graviton} components (4D symmetric tensor)<br>
        - Vectors g_mu,a: {vectors} components (4D vectors from mixed indices)<br>
        - Scalars g_a,b: {scalars} components (moduli from internal metric)
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        f"""
        <div class="step-card">
        <b>Counting:</b><br><br>
        Total independent components: {graviton} + {vectors} + {scalars} = {total}<br>
        Target exponent d*D: {target}<br>
        Mismatch: {total} vs {target}
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    assessment = angle_2.get("assessment", "")
    if assessment:
        st.markdown(
            f"""
            <div class="theorem-card">
            <b>Assessment:</b><br><br>
            {assessment}
            </div>
            """,
            unsafe_allow_html=True,
        )
else:
    st.markdown(
        """
        <div class="step-card">
        <b>KK decomposition of the D=6 metric:</b><br><br>
        Symmetric metric g_MN in D=6 has D(D+1)/2 = 21 independent components.<br><br>
        KK decomposition:<br>
        - Graviton g_mu,nu: 10 components (4D symmetric tensor)<br>
        - Vectors g_mu,a: 8 components (4D vectors from mixed indices)<br>
        - Scalars g_a,b: 3 components (moduli from internal metric)
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="step-card">
        <b>Counting:</b><br><br>
        Total independent components: 10 + 8 + 3 = 21<br>
        Target exponent d*D: 24<br>
        Mismatch: 21 vs 24. The metric counting does not directly
        produce the exponent 24.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# E. Power-Law Running
# ---------------------------------------------------------------------------
st.header("E. Power-Law Running")

if summary:
    angle_3 = summary["angle_3"]

    eff_exp = angle_3.get("effective_exponent")
    exp_str = f"{eff_exp}" if eff_exp is not None else "N/A"

    st.markdown(
        f"""
        <div class="step-card">
        <b>Power-law running of gauge couplings:</b><br><br>
        In D=6, gauge couplings run as power laws above the compactification
        scale due to KK tower contributions. The effective exponent from
        RG running analysis: {exp_str}<br><br>
        {angle_3.get("mechanism", "")}
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    assessment = angle_3.get("assessment", "")
    if assessment:
        st.markdown(
            f"""
            <div class="theorem-card">
            <b>Assessment:</b><br><br>
            {assessment}
            </div>
            """,
            unsafe_allow_html=True,
        )
else:
    st.markdown(
        """
        <div class="step-card">
        <b>Power-law running of gauge couplings:</b><br><br>
        In D=6, gauge couplings run as power laws above the compactification
        scale due to KK tower contributions. The beta function receives
        contributions from all KK modes, producing an enhanced running
        that could generate large hierarchies.<br><br>
        The question is whether this running naturally produces the
        exponent 24 = d*D from first principles.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# F. Induced Gravity (Species Scale)
# ---------------------------------------------------------------------------
st.header("F. Induced Gravity (Species Scale)")

if summary:
    angle_5 = summary["angle_5"]

    n_species = angle_5.get("N_species")
    lambda_uv = angle_5.get("Lambda_UV")

    n_str = str(n_species) if n_species is not None else "N/A"
    uv_str = f"{lambda_uv}" if lambda_uv is not None else "N/A"

    st.markdown(
        f"""
        <div class="step-card">
        <b>Induced gravity from N species:</b><br><br>
        N species: {n_str}<br>
        Lambda_UV: {uv_str}<br><br>
        {angle_5.get("mechanism", "")}
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    assessment = angle_5.get("assessment", "")
    if assessment:
        st.markdown(
            f"""
            <div class="theorem-card">
            <b>Assessment:</b><br><br>
            {assessment}
            </div>
            """,
            unsafe_allow_html=True,
        )
else:
    st.markdown(
        """
        <div class="step-card">
        <b>Induced gravity from N species:</b><br><br>
        With N=24 species (one per metric component or degree of freedom),
        the species scale Lambda_species ~ M_Pl / sqrt(N) modifies the
        effective Planck mass. Induced gravity contributions from N
        species loops generate M_Pl^2 ~ N * Lambda_UV^2.<br><br>
        This could explain why 24 appears as the exponent, if each
        "species" contributes one power of alpha to the hierarchy.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# G. Residual Analysis
# ---------------------------------------------------------------------------
st.header("G. Residual Analysis")

if summary:
    angle_7 = summary["angle_7"]

    # Residual analysis chart
    try:
        from app.components.charts import residual_analysis_chart  # noqa: E402
        fig = residual_analysis_chart(angle_7)
        st.plotly_chart(fig, use_container_width=True)
    except ImportError:
        st.info("Chart function residual_analysis_chart not yet available.")
    except Exception as exc:
        st.warning(f"Residual analysis chart error: {exc}")

    st.markdown("")

    # Correction candidates table
    candidates = angle_7.get("correction_candidates", [])
    if candidates:
        table_rows = []
        for c in candidates:
            table_rows.append({
                "Correction": c.get("name", ""),
                "Type": c.get("type", ""),
                "Magnitude": c.get("magnitude", ""),
                "Could Explain 688 ppm": c.get("viable", ""),
            })
        st.dataframe(
            pd.DataFrame(table_rows),
            use_container_width=True,
            hide_index=True,
        )

    st.markdown("")

    residual_assessment = angle_7.get("assessment", "")
    if residual_assessment:
        st.markdown(
            f"""
            <div class="formula-card">
            <b>Residual assessment:</b><br><br>
            {residual_assessment}
            </div>
            """,
            unsafe_allow_html=True,
        )
else:
    st.markdown(
        """
        <div class="formula-card">
        <b>Residual analysis:</b><br><br>
        The 688 ppm residual between alpha^24 * mu^2 and G_observed / G_natural
        is small but nonzero. Possible sources of the residual:<br><br>
        - Radiative corrections (loop effects at the compactification scale)<br>
        - Threshold corrections from massive KK modes<br>
        - Higher-order terms in the alpha expansion<br>
        - Gauss-Bonnet or curvature-squared corrections to the action<br>
        - Running of alpha and mu between laboratory and compactification scales
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# H. Swampland Distance Conjecture
# ---------------------------------------------------------------------------
st.header("H. Swampland Distance Conjecture")

if summary:
    angle_8 = summary["angle_8"]

    # SDC lambda chart
    try:
        from app.components.charts import sdc_lambda_chart  # noqa: E402
        fig = sdc_lambda_chart(angle_8)
        st.plotly_chart(fig, use_container_width=True)
    except ImportError:
        st.info("Chart function sdc_lambda_chart not yet available.")
    except Exception as exc:
        st.warning(f"SDC chart error: {exc}")

    st.markdown("")

    col_h1, col_h2, col_h3 = st.columns(3)

    with col_h1:
        lam_eh = angle_8.get("lambda_EH")
        st.metric(
            label="lambda_EH (pure EH)",
            value=fmt_decimal(lam_eh) if lam_eh is not None else "~0.354",
        )

    with col_h2:
        lam_gb = angle_8.get("lambda_GB")
        st.metric(
            label="lambda_GB (with GB)",
            value=fmt_decimal(lam_gb) if lam_gb is not None else "~0.728",
        )

    with col_h3:
        lam_bound = angle_8.get("lambda_bound")
        st.metric(
            label="SDC bound 1/sqrt(d-2)",
            value=fmt_decimal(lam_bound) if lam_bound is not None else "~0.707",
        )

    st.markdown("")

    sdc_eh = angle_8.get("sdc_satisfied_EH", False)
    sdc_gb = angle_8.get("sdc_satisfied_GB", True)
    st.markdown(
        f"""
        <div class="theorem-card">
        <b>Swampland Distance Conjecture check:</b><br><br>
        Pure Einstein-Hilbert: lambda_EH = {fmt_decimal(lam_eh) if lam_eh is not None else '~0.354'}
        -- SDC {'satisfied' if sdc_eh else 'NOT satisfied'}<br>
        With Gauss-Bonnet correction: lambda_GB = {fmt_decimal(lam_gb) if lam_gb is not None else '~0.728'}
        -- SDC {'satisfied' if sdc_gb else 'NOT satisfied'}<br><br>
        The GB correction is required for consistency with the Swampland Distance Conjecture.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    col_h4, col_h5, col_h6 = st.columns(3)

    with col_h4:
        omega_t = angle_8.get("omega_target")
        st.metric(
            label="omega_target",
            value=fmt_decimal(omega_t) if omega_t is not None else "~0.118",
        )

    with col_h5:
        omega_s = angle_8.get("omega_saturation")
        st.metric(
            label="omega_saturation",
            value=fmt_decimal(omega_s) if omega_s is not None else "~0.125",
        )

    with col_h6:
        margin = angle_8.get("margin")
        st.metric(
            label="Margin (lambda_GB - bound)",
            value=fmt_decimal(margin) if margin is not None else "~0.021",
        )

    st.markdown("")

    sdc_assessment = angle_8.get("assessment", "")
    if sdc_assessment:
        st.markdown(
            f"""
            <div class="step-card">
            <b>Assessment:</b><br><br>
            {sdc_assessment}
            </div>
            """,
            unsafe_allow_html=True,
        )
else:
    col_h1, col_h2, col_h3 = st.columns(3)

    with col_h1:
        st.metric(label="lambda_EH (pure EH)", value="~0.354")

    with col_h2:
        st.metric(label="lambda_GB (with GB)", value="~0.728")

    with col_h3:
        st.metric(label="SDC bound 1/sqrt(d-2)", value="~0.707")

    st.markdown("")

    st.markdown(
        """
        <div class="theorem-card">
        <b>Swampland Distance Conjecture check:</b><br><br>
        Pure Einstein-Hilbert: lambda_EH ~ 0.354 -- SDC NOT satisfied<br>
        With Gauss-Bonnet correction: lambda_GB ~ 0.728 -- SDC satisfied<br><br>
        The GB correction is required for consistency with the Swampland Distance Conjecture.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    col_h4, col_h5, col_h6 = st.columns(3)

    with col_h4:
        st.metric(label="omega_target", value="~0.118")

    with col_h5:
        st.metric(label="omega_saturation", value="~0.125")

    with col_h6:
        st.metric(label="Margin (lambda_GB - bound)", value="~0.021")

    st.markdown("")

    st.markdown(
        """
        <div class="step-card">
        <b>Assessment:</b><br><br>
        The SDC decay rate lambda must satisfy lambda >= 1/sqrt(d-2) for
        consistency with quantum gravity. Pure Einstein-Hilbert in D=6
        gives lambda_EH ~ 0.354, violating the bound. The Gauss-Bonnet
        correction raises this to lambda_GB ~ 0.728, satisfying the SDC.
        Near-saturation of the bound at the golden-ratio omega hints at
        a deeper structure.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# I. Emergence Proposal
# ---------------------------------------------------------------------------
st.header("I. Emergence Proposal")

if summary:
    angle_9 = summary["angle_9"]

    # Emergence R scan chart
    try:
        from app.components.charts import emergence_R_scan_chart  # noqa: E402
        fig = emergence_R_scan_chart(angle_9)
        st.plotly_chart(fig, use_container_width=True)
    except ImportError:
        st.info("Chart function emergence_R_scan_chart not yet available.")
    except Exception as exc:
        st.warning(f"Emergence chart error: {exc}")

    st.markdown("")

    sum_formula = angle_9.get("sum_formula", "sum_{l=1}^{L} (2l+1)*l*(l+1)")
    scaling = angle_9.get("scaling_large_L", "M_Pl^2 ~ N_pol * L_max^4 / (32*pi^2*R^2)")

    st.markdown(
        f"""
        <div class="step-card">
        <b>KK tower sum and Planck mass scaling:</b><br><br>
        KK angular momentum sum: {sum_formula}<br><br>
        Large-L scaling: {scaling}
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    # Canonical R results table
    canonical = angle_9.get("canonical_R_results", [])
    if canonical:
        table_rows = []
        for c in canonical:
            table_rows.append({
                "Name": c.get("name", ""),
                "R (m)": fmt_decimal(c.get("R_m")) if c.get("R_m") is not None else "",
                "L_max (Planck)": c.get("L_max_planck", ""),
                "M_Pl^2 ratio": fmt_decimal(c.get("M_Pl_sq_ratio")) if c.get("M_Pl_sq_ratio") is not None else "",
            })
        st.dataframe(
            pd.DataFrame(table_rows),
            use_container_width=True,
            hide_index=True,
        )

    st.markdown("")

    r_needed = angle_9.get("R_needed_for_match_m")
    st.metric(
        label="R needed for self-consistency (m)",
        value=fmt_decimal(r_needed) if r_needed is not None else "N/A",
    )

    st.markdown("")

    honest_emerg = angle_9.get("honest_assessment", "")
    if honest_emerg:
        st.markdown(
            f"""
            <div class="warning-card">
            <b>Honest assessment (works = {angle_9.get('works', False)}):</b><br><br>
            {honest_emerg}
            </div>
            """,
            unsafe_allow_html=True,
        )
else:
    st.markdown(
        """
        <div class="step-card">
        <b>KK tower sum and Planck mass scaling:</b><br><br>
        KK angular momentum sum: sum_{l=1}^{L} (2l+1)*l*(l+1)<br><br>
        Large-L scaling: M_Pl^2 ~ N_pol * L_max^4 / (32*pi^2*R^2)
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="warning-card">
        <b>Honest assessment (works = False):</b><br><br>
        The emergence proposal attempts to generate M_Pl^2 from summing
        one-loop contributions of KK modes. While the scaling is correct
        in form, the required compactification radius for self-consistency
        does not match any physically motivated scale. The approach does
        not naturally reproduce alpha^24 * mu^2.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# J. One-Loop KK Power Counting
# ---------------------------------------------------------------------------
st.header("J. One-Loop KK Power Counting")

if summary:
    angle_10 = summary["angle_10"]

    analytic = angle_10.get("analytic_scaling", "alpha_g ~ 32*pi^2 * m_e^2 / (N_pol * Lambda^4 * R^2)")

    st.markdown(
        f"""
        <div class="step-card">
        <b>Analytic scaling from one-loop power counting:</b><br><br>
        UV divergence degree: {angle_10.get('uv_divergence_degree', 4)} (= D - 2)<br><br>
        {analytic}
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    # Specific cases table
    cases = angle_10.get("specific_cases", [])
    if cases:
        table_rows = []
        for c in cases:
            table_rows.append({
                "R": c.get("R", ""),
                "Lambda": c.get("Lambda", ""),
                "alpha_g predicted": fmt_decimal(c.get("alpha_g_predicted")) if c.get("alpha_g_predicted") is not None else "",
                "Ratio to measured": fmt_decimal(c.get("ratio")) if c.get("ratio") is not None else "",
            })
        st.dataframe(
            pd.DataFrame(table_rows),
            use_container_width=True,
            hide_index=True,
        )

    st.markdown("")

    mu_problem = angle_10.get("mu_exponent_problem", "")
    if mu_problem:
        st.markdown(
            f"""
            <div class="warning-card">
            <b>mu exponent problem:</b><br><br>
            {mu_problem}
            </div>
            """,
            unsafe_allow_html=True,
        )

    st.markdown("")

    honest_loop = angle_10.get("honest_assessment", "")
    if honest_loop:
        st.markdown(
            f"""
            <div class="warning-card">
            <b>Honest assessment (works = {angle_10.get('works', False)}):</b><br><br>
            {honest_loop}
            </div>
            """,
            unsafe_allow_html=True,
        )
else:
    st.markdown(
        """
        <div class="step-card">
        <b>Analytic scaling from one-loop power counting:</b><br><br>
        UV divergence degree: 4 (= D - 2)<br><br>
        alpha_g ~ 32*pi^2 * m_e^2 / (N_pol * Lambda^4 * R^2)
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="warning-card">
        <b>mu exponent problem:</b><br><br>
        One-loop power counting produces alpha_g proportional to m_e^2,
        which gives mu^(-2) rather than the required mu^(+2). The sign
        of the mu exponent is wrong, indicating that naive one-loop
        KK power counting cannot reproduce the empirical formula.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="warning-card">
        <b>Honest assessment (works = False):</b><br><br>
        One-loop KK power counting gets the parametric scaling wrong:
        it produces mu^(-2) instead of mu^(+2). No choice of cutoff
        or compactification radius fixes this sign. The approach fails
        to reproduce alpha^24 * mu^2.
        </div>
        """,
        unsafe_allow_html=True,
    )

st.markdown("")

# ---------------------------------------------------------------------------
# K. Honest Assessment
# ---------------------------------------------------------------------------
st.header("K. Honest Assessment")

if summary:
    honest = summary.get("honest_assessment", "")

    st.markdown(
        f"""
        <div class="proof-card">
        <b>Overall assessment:</b><br><br>
        {honest}
        </div>
        """,
        unsafe_allow_html=True,
    )
else:
    st.markdown(
        """
        <div class="proof-card">
        <b>Overall assessment:</b><br><br>
        The formula alpha_g = alpha^24 * mu^2 is a striking numerical fact:
        it reproduces G to 688 ppm with zero fitted parameters, and (d=4, n=2)
        is uniquely selected among all integer dimension pairs.<br><br>
        However, no first-principles derivation currently exists. Standard KK
        volume suppression does not produce the d*D exponent structure. Metric
        counting gives 21, not 24. Power-law running and induced gravity are
        suggestive but incomplete.<br><br>
        The formula remains an empirical observation awaiting theoretical
        explanation. Its uniqueness in the dimension scan is a strong hint
        that something real underlies it, but intellectual honesty requires
        acknowledging that the mechanism is not yet understood.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Hierarchy Derivation | Alpha Ladder Research Dashboard")
