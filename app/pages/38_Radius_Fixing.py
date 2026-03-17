"""
Page 38 -- Radius Fixing: Three New Mechanisms
Coleman-Weinberg, Warped S^2, and Orbifold S^2/Z_2 approaches to fixing a_0.
"""

import sys
import os
import streamlit as st
import pandas as pd


sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from app.components.sidebar import render_sidebar

constants = render_sidebar()

_core_available = False
try:
    from alpha_ladder_core.radius_fixing import (
        summarize_radius_fixing,
        summarize_cw_mechanism,
        summarize_warp_mechanism,
        summarize_orbifold_mechanism,
    )
    _core_available = True
except ImportError:
    pass

_charts_available = False
try:
    from app.components.charts import (
        cw_potential_chart,
        warped_potential_chart,
        orbifold_potential_chart,
    )
    _charts_available = True
except ImportError:
    pass


@st.cache_data
def _get_cw_summary(_constants):
    if _core_available:
        return summarize_cw_mechanism(_constants)
    return None


@st.cache_data
def _get_warp_summary(_constants):
    if _core_available:
        return summarize_warp_mechanism(_constants)
    return None


@st.cache_data
def _get_orbifold_summary(_constants):
    if _core_available:
        return summarize_orbifold_mechanism(_constants)
    return None


@st.cache_data
def _get_overall_summary(_constants):
    if _core_available:
        return summarize_radius_fixing(_constants)
    return None


# ---------------------------------------------------------------------------
# Title and Introduction
# ---------------------------------------------------------------------------
st.title("Radius Fixing: Three New Mechanisms")
st.markdown("*Attempting to fix the internal radius a₀ via quantum and geometric effects*")

st.markdown("""
<div class="warning-card">
<strong>The Problem:</strong> The internal radius a₀ is unfixed due to a classical scaling symmetry
(the Einstein-Hilbert exponent vanishes for n=2, and the Gauss-Bonnet term is topological).
The flux + Casimir balance is monotonic and does not produce a stable minimum.
Here we explore three new mechanisms that might break this degeneracy.
</div>
""", unsafe_allow_html=True)

if not _core_available:
    st.error(
        "Core module `alpha_ladder_core.radius_fixing` not found. "
        "Install alpha_ladder_core to enable computations on this page."
    )

# ---------------------------------------------------------------------------
# Section 1: Coleman-Weinberg Mechanism
# ---------------------------------------------------------------------------
st.header("1. Coleman-Weinberg Mechanism")

cw_summary = _get_cw_summary(constants)

st.markdown("""
<div class="formula-card">
<strong>Coleman-Weinberg one-loop effective potential:</strong><br><br>
V_CW = &Sigma; m⁴(a₀) log(m²(a₀) / &mu;²) / (64 π²)<br><br>
The KK mass spectrum m_n ≈ n / a₀ generates a radius-dependent quantum correction.
Matter content from anomaly-constrained groups (E8×E8: 740 hypers, 496 vectors).
The total potential is V_total = V_CW + V_flux.
</div>
""", unsafe_allow_html=True)

if cw_summary:
    # Key metrics
    c1, c2, c3 = st.columns(3)
    c1.metric("Minimum found", str(cw_summary.get("minimum_found", False)))
    min_a0 = cw_summary.get("minimum_a_0")
    c2.metric("a₀ at min", f"{min_a0:.4e}" if min_a0 is not None else "None")
    c3.metric("Gauge group", cw_summary.get("group_used", "N/A"))

    # CW potential details at a_0 = 1
    cw_pot = cw_summary.get("cw_potential_at_a0_1", {})
    if cw_pot:
        st.subheader("CW Potential Components at a₀ = 1")
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("V_CW", f"{cw_pot.get('V_cw', 0):.4e}")
        c2.metric("Scalar sum", f"{cw_pot.get('scalar_sum', 0):.4e}")
        c3.metric("Fermion sum", f"{cw_pot.get('fermion_sum', 0):.4e}")
        c4.metric("Vector sum", f"{cw_pot.get('vector_sum', 0):.4e}")

    # Table of scan data
    scan_result = cw_summary.get("scan_result", {})
    scan_data = scan_result.get("scan_data", [])
    if scan_data:
        rows = []
        for entry in scan_data:
            rows.append({
                "a₀ (Planck)": f"{entry.get('a_0', 0):.4e}",
                "V_CW": f"{entry.get('V_cw', 0):.4e}",
                "V_flux": f"{entry.get('V_flux', 0):.4e}",
                "V_total": f"{entry.get('V_total', 0):.4e}",
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    # Chart
    if scan_data and _charts_available:
        chart_data = {
            "a_0_values": [e["a_0"] for e in scan_data],
            "V_cw": [e["V_cw"] for e in scan_data],
            "V_flux_min": [e["V_flux"] for e in scan_data],
            "V_total": [e["V_total"] for e in scan_data],
        }
        st.plotly_chart(cw_potential_chart(chart_data), use_container_width=True)

    # Physics summary
    physics = cw_summary.get("physics_summary", "")
    if physics:
        st.markdown(f"""
<div class="step-card">
<strong>Physics:</strong> {physics}
</div>
""", unsafe_allow_html=True)

    # Honest assessment
    honest = cw_summary.get("honest_assessment", "")
    if honest:
        st.markdown(f"""
<div class="theorem-card">
<strong>Result:</strong> {honest}
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from Coleman-Weinberg mechanism.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 2: Warped S^2 Compactification
# ---------------------------------------------------------------------------
st.header("2. Warped S² Compactification")

warp_summary = _get_warp_summary(constants)

st.markdown("""
<div class="formula-card">
<strong>Warp ansatz:</strong> A(&theta;) = &epsilon; cos(&theta;)<br><br>
A small warp factor on S² introduces gradient energy ≈ &epsilon;² / a₀²
and modifies the curvature integral. The warp parameter &epsilon; controls the deviation
from a round sphere. For &epsilon; = 0, this reduces to the unwarped case.
</div>
""", unsafe_allow_html=True)

if warp_summary:
    # Key metrics
    scan_result = warp_summary.get("scan_result", {})
    c1, c2, c3 = st.columns(3)
    c1.metric("Epsilon values tested", str(scan_result.get("n_epsilon", "N/A")))
    c2.metric("Minimum found", str(warp_summary.get("minimum_found", False)))
    best_eps = scan_result.get("best_epsilon")
    c3.metric("Best epsilon", f"{best_eps}" if best_eps is not None else "None")

    # Single-point gradient details
    sp = warp_summary.get("single_point", {})
    if sp:
        st.markdown(f"""
<div class="formula-card">
<strong>Gradient energy at &epsilon;={sp.get('epsilon', 'N/A')}, a₀={sp.get('a_0', 'N/A')}:</strong><br>
V_grad = {sp.get('V_grad', 0):.6e}, V_curv = {sp.get('V_curv', 0):.6e}, V_warp = {sp.get('V_warp', 0):.6e}<br>
Gradient scaling: {sp.get('gradient_scaling', 'N/A')}, Curvature scaling: {sp.get('curvature_scaling', 'N/A')}<br>
Perturbative: {sp.get('is_perturbative', 'N/A')}
</div>
""", unsafe_allow_html=True)

    # Scan data table -- group by epsilon, show a few a_0 samples
    scan_data = scan_result.get("scan_data", [])
    if scan_data:
        # Collect unique epsilons
        epsilons = sorted(set(e["epsilon"] for e in scan_data))
        # Build chart data
        if _charts_available:
            chart_data = {
                "a_0_values": sorted(set(e["a_0"] for e in scan_data)),
                "epsilon_values": epsilons,
                "V_total_by_epsilon": {},
            }
            for eps in epsilons:
                eps_entries = sorted(
                    [e for e in scan_data if e["epsilon"] == eps],
                    key=lambda x: x["a_0"],
                )
                chart_data["V_total_by_epsilon"][eps] = [e["V_total"] for e in eps_entries]
            st.plotly_chart(warped_potential_chart(chart_data), use_container_width=True)

        # Summary table: one row per epsilon (sample at a_0=1.0)
        rows = []
        for eps in epsilons:
            eps_at_1 = [e for e in scan_data if e["epsilon"] == eps and abs(e["a_0"] - 1.0) < 0.5]
            if eps_at_1:
                e = eps_at_1[0]
                rows.append({
                    "epsilon": f"{eps}",
                    "V_flux (a₀≈1)": f"{e.get('V_flux', 0):.4e}",
                    "V_warp (a₀≈1)": f"{e.get('V_warp', 0):.4e}",
                    "V_total (a₀≈1)": f"{e.get('V_total', 0):.4e}",
                })
        if rows:
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    # Physics summary and assessment
    physics = warp_summary.get("physics_summary", "")
    if physics:
        st.markdown(f"""
<div class="step-card">
<strong>Physics:</strong> {physics}
</div>
""", unsafe_allow_html=True)

    honest = warp_summary.get("honest_assessment", "")
    if honest:
        st.markdown(f"""
<div class="theorem-card">
<strong>Result:</strong> {honest}
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from warped mechanism.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 3: Orbifold S^2/Z_2
# ---------------------------------------------------------------------------
st.header("3. Orbifold S²/Z₂")

orbifold_summary = _get_orbifold_summary(constants)

st.markdown("""
<div class="theorem-card">
<strong>Z₂ orbifold construction:</strong> Identify antipodal points on S²,
producing S²/Z₂ = RP². Only even-l KK modes survive (l=0,2,4,...),
modifying the Casimir coefficient. Brane tension at the 2 fixed points contributes
V_brane = 2 T₀ / a₀², which scales differently from V_casimir ≈ a₀⁻⁴.
</div>
""", unsafe_allow_html=True)

if orbifold_summary:
    # Key metrics
    orb_cas = orbifold_summary.get("orbifold_casimir", {})
    crit = orbifold_summary.get("critical_tension", {})
    c1, c2, c3 = st.columns(3)
    c1.metric("Minimum found", str(orbifold_summary.get("minimum_found", False)))
    t0_crit = crit.get("T_0_critical")
    c2.metric("Critical T_0", f"{t0_crit:.4e}" if t0_crit is not None else "None found")
    zeta_val = orb_cas.get("zeta_orb_value")
    c3.metric("ζ_orb(-1/2)", f"{zeta_val:.6f}" if zeta_val is not None else "N/A")

    # Orbifold Casimir details
    if orb_cas:
        st.markdown(f"""
<div class="formula-card">
<strong>Orbifold Casimir:</strong> ζ_orb(-1/2) = {orb_cas.get('zeta_orb_value', 'N/A'):.6f}
(cf. full S²: ζ(-1/2) ≈ -0.25).<br>
Coefficient: {orb_cas.get('casimir_coefficient', 0):.6e}, Modes included: {orb_cas.get('n_modes_included', 'N/A')}<br>
Sign differs from full S²: {orb_cas.get('differs_from_s2', 'N/A')}
</div>
""", unsafe_allow_html=True)

    # Radius scan table
    rscan = orbifold_summary.get("radius_scan", {})
    scan_data = rscan.get("scan_data", [])
    if scan_data:
        # Build chart from first T_0 value
        t0_values = sorted(set(e["T_0"] for e in scan_data))
        if t0_values and _charts_available:
            first_t0 = t0_values[len(t0_values) // 2]  # pick middle T_0
            t0_entries = sorted(
                [e for e in scan_data if e["T_0"] == first_t0],
                key=lambda x: x["a_0"],
            )
            chart_data = {
                "a_0_values": [e["a_0"] for e in t0_entries],
                "V_casimir": [e["V_casimir_orb"] for e in t0_entries],
                "V_brane": [e["V_brane"] for e in t0_entries],
                "V_total": [e["V_total"] for e in t0_entries],
                "T_0": first_t0,
            }
            st.plotly_chart(orbifold_potential_chart(chart_data), use_container_width=True)

        # T_0 scan summary table
        rows = []
        for t0 in t0_values:
            t0_entries = [e for e in scan_data if e["T_0"] == t0]
            # Check if V_total has a minimum (non-monotonic)
            v_vals = [e["V_total"] for e in sorted(t0_entries, key=lambda x: x["a_0"])]
            has_min = False
            min_a0 = "N/A"
            min_v = "N/A"
            for i in range(1, len(v_vals) - 1):
                if v_vals[i] < v_vals[i - 1] and v_vals[i] < v_vals[i + 1]:
                    has_min = True
                    sorted_entries = sorted(t0_entries, key=lambda x: x["a_0"])
                    min_a0 = f"{sorted_entries[i]['a_0']:.4e}"
                    min_v = f"{sorted_entries[i]['V_total']:.4e}"
                    break
            rows.append({
                "T_0": f"{t0:.4e}",
                "Min exists": str(has_min),
                "a₀ at min": min_a0,
                "V at min": min_v,
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    # Critical tension info
    if crit:
        st.markdown(f"""
<div class="step-card">
<strong>Critical brane tension search:</strong> {crit.get('description', 'N/A')}
</div>
""", unsafe_allow_html=True)

    # Physics summary and assessment
    physics = orbifold_summary.get("physics_summary", "")
    if physics:
        st.markdown(f"""
<div class="step-card">
<strong>Physics:</strong> {physics}
</div>
""", unsafe_allow_html=True)

    honest = orbifold_summary.get("honest_assessment", "")
    if honest:
        st.markdown(f"""
<div class="theorem-card">
<strong>Result:</strong> {honest}
</div>
""", unsafe_allow_html=True)
elif _core_available:
    st.info("No data returned from orbifold mechanism.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Section 4: Summary -- What Could Fix a_0
# ---------------------------------------------------------------------------
st.header("4. What Could Fix a₀")

overall = _get_overall_summary(constants)

if overall:
    # Mechanism summary cards
    mech_data = [
        ("Coleman-Weinberg", overall.get("mechanism_1_cw", {})),
        ("Warped S²", overall.get("mechanism_2_warp", {})),
        ("Orbifold S²/Z₂", overall.get("mechanism_3_orbifold", {})),
    ]

    rows = []
    for name, mech in mech_data:
        found = mech.get("minimum_found", False)
        icon = "PASS" if found else "FAIL"
        phys = mech.get("physics_summary", "")
        # Truncate physics summary for card display
        short_phys = phys
        st.markdown(f"""
<div class="step-card">
<strong>[{icon}] {mech.get('mechanism_name', name)}</strong><br>
{short_phys}
</div>
""", unsafe_allow_html=True)
        rows.append({
            "Mechanism": name,
            "Fixes a₀": "Yes" if found else "No",
            "Key physics": mech.get("honest_assessment", ""),
        })

    # Summary table
    st.subheader("Mechanism Comparison")
    df = pd.DataFrame(rows)
    st.markdown(
        df.to_html(index=False, classes="wrap-table", escape=False),
        unsafe_allow_html=True,
    )

    # Overall verdict
    verdict = overall.get("overall_assessment", "")
    if verdict:
        st.markdown(f"""
<div class="proof-card">
<strong>Overall Verdict:</strong><br><br>
{verdict}
</div>
""", unsafe_allow_html=True)

    # Honest assessment
    honest = overall.get("honest_assessment", "")
    if honest:
        st.markdown(f"""
<div class="warning-card">
<strong>Honest Assessment:</strong><br><br>
{honest}
</div>
""", unsafe_allow_html=True)

    # Key metrics row
    any_fixes = overall.get("any_mechanism_fixes_radius", False)
    best = overall.get("best_mechanism")
    c1, c2, c3 = st.columns(3)
    c1.metric("Mechanisms tested", "3")
    c2.metric("Any fix a₀", "Yes" if any_fixes else "No")
    c3.metric("Best mechanism", best if best else "None")
elif _core_available:
    st.info("No summary data returned.")
else:
    st.info("Core module not available. Install alpha_ladder_core.")

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Radius Fixing | Alpha Ladder Research Dashboard")
