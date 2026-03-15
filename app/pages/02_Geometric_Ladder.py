"""
Geometric Ladder -- Page 2

The alpha ladder maps powers of the fine-structure constant to physical
scales. Starting from r_e = alpha * lambda_bar_c = alpha^2 * a_0, each
rung multiplies by alpha (~1/137) in length or 1/alpha (~137) in energy.

The ladder connects atomic physics (rung 0) through nuclear/QCD scales
(rungs 1-2), electroweak physics (rung 2.5), the GUT scale (rung 9),
the Planck scale (rung 11.5), and gravity (rung 21).
"""

import streamlit as st
import pandas as pd
from decimal import Decimal, getcontext

getcontext().prec = 50


# ---------------------------------------------------------------------------
# Custom CSS
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
    .step-card {
        background-color: #1a1d23;
        border-left: 3px solid #60a5fa;
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
    .warning-card {
        background-color: #1a1d23;
        border-left: 3px solid #f87171;
        padding: 0.8rem 1rem;
        margin: 0.4rem 0;
        border-radius: 0 8px 8px 0;
    }
    .rung-table {
        font-family: 'Fira Mono', Consolas, monospace;
        font-size: 1.1rem;
    }
    .bridge-card {
        background-color: #1a1d23;
        border: 1px solid #2e3440;
        border-radius: 8px;
        padding: 1rem;
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
from app.components.charts import ladder_chart, sigma_heatmap  # noqa: E402
from app.components.formatting import (  # noqa: E402
    fmt_decimal,
    fmt_percent,
    fmt_sigma,
    color_by_quality,
)

constants = render_sidebar()

# ---------------------------------------------------------------------------
# Core module imports (guarded)
# ---------------------------------------------------------------------------
_core_available = False
_rung_map_available = False
try:
    from alpha_ladder_core.ladder import calculate_geometric_rungs, compute_physical_rung_map
    _core_available = True
    _rung_map_available = True
except ImportError:
    try:
        from alpha_ladder_core.ladder import calculate_geometric_rungs
        _core_available = True
    except ImportError:
        pass

# ---------------------------------------------------------------------------
# Cached computation for physical rung map
# ---------------------------------------------------------------------------
rung_map = None
if _rung_map_available and constants is not None:
    @st.cache_data(show_spinner="Computing physical rung map...")
    def _get_rung_map(_constants):
        return compute_physical_rung_map(_constants)

    try:
        rung_map = _get_rung_map(constants)
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Cached computation for ladder data (bridge candidates)
# ---------------------------------------------------------------------------
@st.cache_data(show_spinner="Computing geometric ladder...")
def _compute_ladder_fallback():
    """Compute the geometric ladder using hardcoded CODATA 2018 values."""
    alpha = Decimal("0.0072973525693")
    alpha_g = Decimal("1.7512e-45")
    pi = Decimal("3.14159265358979323846264338327950288419716939937510")
    phi = (1 + Decimal(5).sqrt()) / 2
    e = Decimal("2.71828182845904523536028747135266249775724709369995")

    rungs = []
    for n in range(1, 25):
        val = alpha ** n
        ratio = val / alpha_g
        label = ""
        if n == 10:
            label = "Dark Matter Candidate (The Blank)"
        elif n == 21:
            label = "Gravity Floor (The Bridge)"
        rungs.append({
            "power": n,
            "value": val,
            "ratio_to_gravity": ratio,
            "label": label,
        })

    alpha_21 = alpha ** 21

    bridges = {
        "phi * alpha^21": phi * alpha_21,
        "(phi^2/2) * alpha^21": (phi ** 2 / 2) * alpha_21,
        "(5/12) * pi * alpha^21": (Decimal(5) / Decimal(12)) * pi * alpha_21,
        "sqrt(e) / cbrt(2) * alpha^21": (
            e.sqrt() / Decimal(2) ** (Decimal(1) / Decimal(3))
        ) * alpha_21,
    }

    return rungs, bridges, alpha_g


@st.cache_data(show_spinner="Computing geometric ladder...")
def _compute_ladder_from_core(_constants):
    """Compute ladder using the core module."""
    result = calculate_geometric_rungs(_constants)
    return result


def _get_ladder_data():
    """Get ladder data from core module or fallback."""
    if _core_available and constants is not None:
        try:
            result = _compute_ladder_from_core(constants)
            if isinstance(result, tuple) and len(result) >= 2:
                return result[0], result[1], getattr(constants, "alpha_g", Decimal("1.7512e-45"))
            return result, {}, Decimal("1.7512e-45")
        except Exception as exc:
            st.warning(f"Core module error: {exc}. Using fallback.")
    return _compute_ladder_fallback()


rungs_data, bridges, alpha_g = _get_ladder_data()

# ===================================================================
# Header
# ===================================================================
st.title("Geometric Ladder")
st.markdown(
    "The alpha ladder maps powers of the fine-structure constant to physical "
    "scales. Starting from r_e = alpha * lambda_bar_c = alpha^2 * a_0, each "
    "rung multiplies by 1/alpha (~137) in energy."
)
st.divider()

# ===================================================================
# Section A: The Base Identity
# ===================================================================
st.subheader("A. The Base Identity")

if rung_map is not None:
    bi = rung_map["base_identity"]
    st.markdown(
        '<div class="formula-card">'
        "<strong>r_e = alpha * lambda_bar_c = alpha^2 * a_0</strong><br>"
        f"Relative error between the two expressions: {bi['relative_error']:.2e}"
        "</div>",
        unsafe_allow_html=True,
    )
    col_a, col_b, col_c = st.columns(3)
    col_a.metric("a_0 (Bohr radius)", f"{bi['a_0']:.6e} m")
    col_b.metric("lambda_bar_c (Compton)", f"{bi['lambda_bar_c']:.6e} m")
    col_c.metric("r_e (classical radius)", f"{bi['r_e_from_compton']:.6e} m")
else:
    st.markdown(
        '<div class="formula-card">'
        "<strong>r_e = alpha * lambda_bar_c = alpha^2 * a_0</strong><br>"
        "a_0 = 5.29177e-11 m (Bohr radius)<br>"
        "lambda_bar_c = 3.86159e-13 m (reduced Compton wavelength)<br>"
        "r_e = 2.81794e-15 m (classical electron radius)<br>"
        "Each step down by alpha reduces the length by a factor of ~137."
        "</div>",
        unsafe_allow_html=True,
    )

st.divider()

# ===================================================================
# Section B: Physical Rung Map
# ===================================================================
st.subheader("B. Physical Rung Map")

st.markdown(
    "Each rung k corresponds to a length scale alpha^k * a_0 and an energy "
    "scale E_k = hbar*c / (alpha^k * a_0). The table shows which known "
    "physical scales land near each rung."
)


def _dex_color(dex):
    """Return a CSS color based on match quality in dex."""
    if dex < 0.2:
        return "#34d399"  # green
    elif dex < 0.5:
        return "#f59e0b"  # amber
    elif dex <= 1.0:
        return "#fb923c"  # orange
    else:
        return "#f87171"  # red


if rung_map is not None:
    all_rungs = rung_map["rungs"]

    # Show rungs 0-12 by default
    def _build_rung_table(rung_list):
        rows = []
        for r in rung_list:
            k = r["k"]
            np_ = r["nearest_particle"]
            ns_ = r["nearest_scale"]
            # Pick whichever match is closer
            if np_["offset_dex"] <= ns_["offset_dex"]:
                match_name = np_["name"]
                match_dex = np_["offset_dex"]
            else:
                match_name = ns_["name"]
                match_dex = ns_["offset_dex"]

            rows.append({
                "Rung k": k,
                "Length (m)": f"{r['length_m']:.3e}",
                "Energy (eV)": f"{r['energy_eV']:.3e}",
                "Nearest Physical Scale": match_name if match_name else "--",
                "Match (dex)": f"{match_dex:.2f}" if match_name else "--",
            })
        return rows

    # Primary table: rungs 0-12
    primary_rungs = [r for r in all_rungs if r["k"] <= 12]
    primary_rows = _build_rung_table(primary_rungs)
    df_primary = pd.DataFrame(primary_rows)

    st.dataframe(
        df_primary,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Rung k": st.column_config.NumberColumn(width="small"),
            "Length (m)": st.column_config.TextColumn(width="medium"),
            "Energy (eV)": st.column_config.TextColumn(width="medium"),
            "Nearest Physical Scale": st.column_config.TextColumn(width="large"),
            "Match (dex)": st.column_config.TextColumn(width="small"),
        },
    )

    # Full table in expander
    with st.expander("Full rung table (k = 0 to 30)"):
        full_rows = _build_rung_table(all_rungs)
        df_full = pd.DataFrame(full_rows)
        st.dataframe(df_full, use_container_width=True, hide_index=True)

else:
    # Fallback: static content
    st.markdown(
        '<div class="proof-card">'
        "<strong>Physical rung map (representative values, CODATA 2018)</strong><br><br>"
        "<table class='rung-table' style='width:100%; border-collapse:collapse;'>"
        "<tr style='border-bottom:1px solid #2e3440;'>"
        "<th>Rung k</th><th>Length (m)</th><th>Energy (eV)</th><th>Nearest Scale</th><th>dex</th></tr>"
        "<tr><td>0</td><td>5.29e-11</td><td>3.73e+03</td><td>--</td><td>--</td></tr>"
        "<tr><td>1</td><td>3.86e-13</td><td>5.11e+05</td><td>Electron</td><td>0.00</td></tr>"
        "<tr><td>2</td><td>2.82e-15</td><td>7.00e+07</td><td>Muon / Strange</td><td>~0.2</td></tr>"
        "<tr><td>3</td><td>2.06e-17</td><td>9.59e+09</td><td>W / Z / Higgs</td><td>~0.5</td></tr>"
        "<tr><td>4</td><td>1.50e-19</td><td>1.31e+12</td><td>--</td><td>--</td></tr>"
        "<tr><td>9</td><td>~1e-32</td><td>~1e+25</td><td>GUT scale</td><td>~0.0</td></tr>"
        "<tr><td>12</td><td>~1e-38</td><td>~1e+31</td><td>near Planck</td><td>~0.5</td></tr>"
        "</table>"
        "</div>",
        unsafe_allow_html=True,
    )

st.divider()

# ===================================================================
# Section C: Where Particles Live
# ===================================================================
st.subheader("C. Where Particles Live")

st.markdown(
    "Each Standard Model particle has a mass that corresponds to a "
    "fractional rung position k = log(m / m_e) / log(1/alpha). "
    "Half-integer rungs turn out to be physically significant."
)

if rung_map is not None:
    pr = rung_map["particle_rungs"]
    pr_rows = []
    for p in pr:
        pr_rows.append({
            "Particle": p["name"],
            "Mass (eV)": f"{p['mass_eV']:.4e}",
            "Rung Position": f"{p['rung_position']:.3f}",
            "Nearest Half-Rung": f"{p['nearest_half']:.1f}",
            "Residual": f"{p['residual']:+.3f}",
        })

    df_particles = pd.DataFrame(pr_rows)
    st.dataframe(
        df_particles,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Particle": st.column_config.TextColumn(width="medium"),
            "Mass (eV)": st.column_config.TextColumn(width="medium"),
            "Rung Position": st.column_config.TextColumn(width="small"),
            "Nearest Half-Rung": st.column_config.TextColumn(width="small"),
            "Residual": st.column_config.TextColumn(width="small"),
        },
    )

    st.markdown(
        '<div class="step-card">'
        "The electroweak sector (W, Z, Higgs, Top) clusters at rung 2.5"
        "</div>",
        unsafe_allow_html=True,
    )
    st.markdown(
        '<div class="step-card">'
        "The proton sits at rung 1.53 -- set by QCD, not QED"
        "</div>",
        unsafe_allow_html=True,
    )
    st.markdown(
        '<div class="step-card">'
        "The GUT scale lands at rung 9.0 -- almost exactly integer"
        "</div>",
        unsafe_allow_html=True,
    )

else:
    st.markdown(
        '<div class="proof-card">'
        "<strong>SM particle rung positions (approximate)</strong><br><br>"
        "Electron: k = 0.00 (definition)<br>"
        "Muon: k = 1.09<br>"
        "Proton: k = 1.53<br>"
        "W boson: k = 2.47<br>"
        "Z boson: k = 2.50<br>"
        "Higgs: k = 2.57<br>"
        "Top quark: k = 2.64<br>"
        "</div>",
        unsafe_allow_html=True,
    )
    st.markdown(
        '<div class="step-card">'
        "The electroweak sector (W, Z, Higgs, Top) clusters at rung 2.5"
        "</div>",
        unsafe_allow_html=True,
    )
    st.markdown(
        '<div class="step-card">'
        "The proton sits at rung 1.53 -- set by QCD, not QED"
        "</div>",
        unsafe_allow_html=True,
    )
    st.markdown(
        '<div class="step-card">'
        "The GUT scale lands at rung 9.0 -- almost exactly integer"
        "</div>",
        unsafe_allow_html=True,
    )

st.divider()

# ===================================================================
# Section D: Special Scales
# ===================================================================
st.subheader("D. Special Scales")

st.markdown(
    "Beyond the Standard Model particles, several landmark energy scales "
    "also map onto the alpha ladder."
)

if rung_map is not None:
    sr = rung_map["scale_rungs"]
    sr_rows = []
    for s in sr:
        rp = s["rung_position"]
        rp_str = f"{rp:.2f}" if rp == rp else "N/A"  # NaN check
        sr_rows.append({
            "Scale": s["name"],
            "Energy (eV)": f"{s['energy_eV']:.3e}",
            "Rung Position": rp_str,
        })

    df_scales = pd.DataFrame(sr_rows)
    st.dataframe(
        df_scales,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Scale": st.column_config.TextColumn(width="large"),
            "Energy (eV)": st.column_config.TextColumn(width="medium"),
            "Rung Position": st.column_config.TextColumn(width="small"),
        },
    )

    k_pl = rung_map["k_planck"]
    k_gr = rung_map["k_gravity"]
    gap = k_gr - k_pl

    col_d1, col_d2, col_d3 = st.columns(3)
    col_d1.metric("k_planck", f"{k_pl:.2f}")
    col_d2.metric("k_gravity", f"{k_gr:.2f}")
    col_d3.metric("Gap (rungs)", f"{gap:.2f}")

    st.markdown(
        '<div class="formula-card">'
        f"Gravity at effective rung {k_gr:.2f}. "
        f"Planck length at rung {k_pl:.2f}. "
        f"The gap of ~{gap:.1f} rungs between Planck and gravity is unexplained."
        "</div>",
        unsafe_allow_html=True,
    )

else:
    st.markdown(
        '<div class="proof-card">'
        "<strong>Special scales (approximate rung positions)</strong><br><br>"
        "QCD Lambda (~200 MeV): k ~ 1.2<br>"
        "Weak scale v (246 GeV): k ~ 2.6<br>"
        "GUT scale (~1e16 GeV): k ~ 9.0<br>"
        "Planck mass: k ~ 11.5<br>"
        "Neutrino mass (~0.05 eV): k ~ -1.9<br>"
        "Hubble scale: k ~ -17<br>"
        "</div>",
        unsafe_allow_html=True,
    )

    st.markdown(
        '<div class="formula-card">'
        "Gravity at effective rung ~20.95. "
        "Planck length at rung ~11.47. "
        "The gap of ~9.5 rungs between Planck and gravity is unexplained."
        "</div>",
        unsafe_allow_html=True,
    )

st.divider()

# ===================================================================
# Section E: The Ladder Chart
# ===================================================================
st.subheader("E. Alpha^n Ladder (Log Scale)")

fig = ladder_chart(rungs_data)
st.plotly_chart(fig, use_container_width=True)

st.divider()

# ===================================================================
# Section F: Bridge Candidates
# ===================================================================
st.subheader("F. Bridge Candidates: alpha_G = C * alpha^21")

st.markdown(
    f"**Measured alpha_G** = {fmt_decimal(alpha_g, sig_figs=6)}"
)

# Experimental G measurements for sigma comparison
G_measurements = {
    "CODATA 2018": (Decimal("6.67430e-11"), Decimal("0.00015e-11")),
    "Quinn 2013 (BIPM)": (Decimal("6.67545e-11"), Decimal("0.00018e-11")),
    "Rosi 2014": (Decimal("6.67191e-11"), Decimal("0.00099e-11")),
    "Newman 2014": (Decimal("6.67435e-11"), Decimal("0.00013e-11")),
    "Li 2018 (HUST-A)": (Decimal("6.67418e-11"), Decimal("0.00009e-11")),
    "Li 2018 (HUST-B)": (Decimal("6.67484e-11"), Decimal("0.00009e-11")),
    "CODATA 2014": (Decimal("6.67408e-11"), Decimal("0.00031e-11")),
}

# G prediction helper
hbar = Decimal("1.054571817e-34")
c_light = Decimal("299792458")
m_e = Decimal("9.1093837015e-31")


def _predict_G(bridge_val):
    """Predict G from a bridge candidate value (alpha_G = bridge_val)."""
    return bridge_val * hbar * c_light / m_e ** 2


bridge_rows = []
sigma_matrix = []
bridge_names_list = []

for name, pred_alpha_g in bridges.items():
    residual = abs(float(alpha_g) - float(pred_alpha_g)) / float(alpha_g) * 100
    G_pred = _predict_G(pred_alpha_g)

    bridge_rows.append({
        "Bridge": name,
        "Predicted alpha_G": fmt_decimal(pred_alpha_g, sig_figs=8),
        "Residual (%)": f"{residual:.4f}%",
        "Predicted G": fmt_decimal(G_pred, sig_figs=6),
    })

    # Compute sigma against each experiment
    sigmas_row = []
    for exp_name, (G_exp, G_unc) in G_measurements.items():
        diff = G_pred - G_exp
        sigma_val = float(abs(diff) / G_unc) if G_unc != 0 else float("inf")
        sigmas_row.append(round(sigma_val, 1))

    sigma_matrix.append(sigmas_row)
    bridge_names_list.append(name.split(" * ")[0].strip() if " * " in name else name)

# Display bridge table
df_bridges = pd.DataFrame(bridge_rows)
st.dataframe(df_bridges, use_container_width=True, hide_index=True)

# Display sigma heatmap
with st.expander("Sigma Deviations vs Experimental Measurements", expanded=True):
    fig_sigma = sigma_heatmap(sigma_matrix, bridge_names_list)
    # Override experiment labels to use actual names
    exp_names = list(G_measurements.keys())
    fig_sigma.update_layout(
        xaxis=dict(
            tickvals=list(range(len(exp_names))),
            ticktext=exp_names,
            tickangle=-45,
        )
    )
    fig_sigma.data[0].x = exp_names
    st.plotly_chart(fig_sigma, use_container_width=True)

# ---------------------------------------------------------------------------
# Per-bridge detailed sigma breakdown
# ---------------------------------------------------------------------------
with st.expander("Detailed Sigma Breakdown per Bridge"):
    for name, pred_alpha_g in bridges.items():
        G_pred = _predict_G(pred_alpha_g)
        st.markdown(f"**{name}** -- G = {fmt_decimal(G_pred, sig_figs=8)}")

        detail_rows = []
        for exp_name, (G_exp, G_unc) in G_measurements.items():
            diff = G_pred - G_exp
            sigma_val = float(abs(diff) / G_unc) if G_unc != 0 else float("inf")
            direction = "+" if diff > 0 else "-"
            err_pct = float(abs(diff) / G_exp) * 100

            detail_rows.append({
                "Experiment": exp_name,
                "Sigma": fmt_sigma(sigma_val, direction),
                "Delta": f"{float(diff):.3e}",
                "Quality": err_pct,
            })

        df_detail = pd.DataFrame(detail_rows)
        st.dataframe(df_detail, use_container_width=True, hide_index=True)

st.divider()
st.caption("Geometric Ladder | Alpha Ladder Research Dashboard")
