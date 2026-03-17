"""
Reusable Plotly chart builders for the Alpha Ladder dashboard.

All charts use a consistent dark theme with scientific styling:
dark background, visible grid lines, and monospace fonts for numeric labels.
"""

import plotly.graph_objects as go

# ---------------------------------------------------------------------------
# Shared layout configuration
# ---------------------------------------------------------------------------
_DARK_THEME = dict(
    paper_bgcolor="#0e1117",
    plot_bgcolor="#1a1d23",
    font=dict(family="Source Sans 3, Segoe UI, sans-serif", color="#e0e0e0", size=12),
    title_font=dict(family="Crimson Pro, Georgia, serif", size=16, color="#ffffff"),
    xaxis=dict(
        gridcolor="#262b33",
        zerolinecolor="#3b4252",
        tickfont=dict(family="Fira Mono, Consolas, monospace"),
        showline=True,
        linecolor="#3b4252",
        linewidth=1,
    ),
    yaxis=dict(
        gridcolor="#262b33",
        zerolinecolor="#3b4252",
        tickfont=dict(family="Fira Mono, Consolas, monospace"),
        showline=True,
        linecolor="#3b4252",
        linewidth=1,
    ),
    margin=dict(l=60, r=30, t=50, b=50),
)

# Accent palette (Tailwind-inspired)
_COLOR_DEFAULT = "#60a5fa"      # blue-400
_COLOR_HIGHLIGHT_10 = "#a78bfa"  # violet-400
_COLOR_HIGHLIGHT_21 = "#f59e0b"  # amber-500
_COLOR_GREEN = "#34d399"         # emerald-400
_COLOR_RED = "#f87171"           # red-400


def _apply_theme(fig):
    """Apply the shared dark theme to a Plotly figure."""
    fig.update_layout(**_DARK_THEME)
    return fig


# ---------------------------------------------------------------------------
# Chart: Geometric Ladder (log-scale bar chart of alpha^n)
# ---------------------------------------------------------------------------
def ladder_chart(rungs_data):
    """Create a log-scale bar chart of alpha^n values.

    Parameters
    ----------
    rungs_data : list[dict]
        Each dict should have keys: "power" (int), "value" (Decimal or float),
        and optionally "label" (str).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    powers = [r["power"] for r in rungs_data]
    values = [float(r["value"]) for r in rungs_data]
    labels = [r.get("label", "") for r in rungs_data]

    colors = []
    for p in powers:
        if p == 10:
            colors.append(_COLOR_HIGHLIGHT_10)
        elif p == 21:
            colors.append(_COLOR_HIGHLIGHT_21)
        else:
            colors.append(_COLOR_DEFAULT)

    hover_texts = []
    for p, v, lbl in zip(powers, values, labels):
        text = f"α<sup>{p}</sup> = {v:.6e}"
        if lbl:
            text += f"<br>{lbl}"
        hover_texts.append(text)

    fig = go.Figure(
        data=go.Bar(
            x=[f"n={p}" for p in powers],
            y=values,
            marker_color=colors,
            hovertext=hover_texts,
            hoverinfo="text",
        )
    )

    fig.update_layout(
        title="Geometric Ladder: αⁿ",
        xaxis_title="Rung (n)",
        yaxis_title="Value",
        yaxis_type="log",
    )

    # Annotate highlighted rungs
    for p, v, color, lbl in zip(powers, values, colors, labels):
        if lbl:
            fig.add_annotation(
                x=f"n={p}",
                y=v,
                text=lbl,
                showarrow=True,
                arrowhead=2,
                arrowcolor=color,
                font=dict(color=color, size=10),
                yshift=15,
            )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Sigma Heatmap (G predictions vs experiments)
# ---------------------------------------------------------------------------
def sigma_heatmap(comparisons, bridge_names):
    """Create a heatmap of sigma deviations for G predictions.

    Parameters
    ----------
    comparisons : list[list[float]]
        2D array where comparisons[i][j] is the sigma deviation for
        bridge i against experiment j.
    bridge_names : list[str]
        Names of bridge candidates (rows).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not comparisons or not bridge_names:
        fig = go.Figure()
        fig.add_annotation(text="No data available", showarrow=False)
        return _apply_theme(fig)

    # Derive experiment names from the first dimension
    n_experiments = len(comparisons[0]) if comparisons else 0
    experiment_labels = [f"Exp {i+1}" for i in range(n_experiments)]

    fig = go.Figure(
        data=go.Heatmap(
            z=comparisons,
            x=experiment_labels,
            y=bridge_names,
            colorscale=[
                [0.0, "#00c853"],
                [0.25, "#ffd600"],
                [0.5, "#ff6d00"],
                [1.0, "#d50000"],
            ],
            colorbar=dict(title="| σ |"),
            hovertemplate="Bridge: %{y}<br>Experiment: %{x}<br>Sigma: %{z:.1f}<extra></extra>",
        )
    )

    fig.update_layout(
        title="Sigma Deviations: G Predictions vs Measurements",
        xaxis_title="Experiment",
        yaxis_title="Bridge Candidate",
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Particle Number Line (particles on alpha-rung number line)
# ---------------------------------------------------------------------------
def particle_number_line(harmonics_data):
    """Plot particles on an alpha-rung number line.

    Parameters
    ----------
    harmonics_data : list[dict]
        Each dict should have keys: "name" (str), "rung" (float),
        and optionally "match" (bool).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not harmonics_data:
        fig = go.Figure()
        fig.add_annotation(text="No harmonics data available", showarrow=False)
        return _apply_theme(fig)

    names = [h["name"] for h in harmonics_data]
    rungs = [float(h["rung"]) for h in harmonics_data]
    is_match = [h.get("match", False) for h in harmonics_data]

    colors = [_COLOR_GREEN if m else _COLOR_RED for m in is_match]
    sizes = [14 if m else 10 for m in is_match]

    fig = go.Figure(
        data=go.Scatter(
            x=rungs,
            y=[0] * len(rungs),
            mode="markers+text",
            marker=dict(color=colors, size=sizes, symbol="diamond"),
            text=names,
            textposition="top center",
            textfont=dict(size=9, color="#e0e0e0"),
            hovertemplate="%{text}<br>Rung: %{x:.4f}<extra></extra>",
        )
    )

    # Add integer rung lines
    max_rung = max(rungs) if rungs else 5
    for i in range(int(max_rung) + 2):
        fig.add_vline(
            x=i, line_dash="dot", line_color="#3b4252", line_width=1
        )

    fig.update_layout(
        title="Particle Masses on the Alpha Rung Number Line",
        xaxis_title="Rung n  (mass ratio = (1/α)ⁿ)",
        yaxis=dict(visible=False, range=[-1, 1]),
        showlegend=False,
        height=350,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Spacing Score (bar chart of spacing scores)
# ---------------------------------------------------------------------------
def spacing_score_chart(spacing_results):
    """Create a bar chart of rung-spacing scores.

    Parameters
    ----------
    spacing_results : list[dict]
        Each dict should have keys: "spacing" (float or str),
        "score" (float), and optionally "label" (str).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not spacing_results:
        fig = go.Figure()
        fig.add_annotation(text="No spacing data available", showarrow=False)
        return _apply_theme(fig)

    labels = [
        str(r.get("label", r.get("spacing", "?"))) for r in spacing_results
    ]
    scores = [float(r["score"]) for r in spacing_results]

    # Color by score quality (lower is better)
    max_score = max(scores) if scores else 1
    colors = []
    for s in scores:
        ratio = s / max_score if max_score > 0 else 0
        if ratio < 0.3:
            colors.append(_COLOR_GREEN)
        elif ratio < 0.6:
            colors.append(_COLOR_HIGHLIGHT_21)
        else:
            colors.append(_COLOR_RED)

    fig = go.Figure(
        data=go.Bar(
            x=labels,
            y=scores,
            marker_color=colors,
            hovertemplate="Spacing: %{x}<br>Score: %{y:.4f}<extra></extra>",
        )
    )

    fig.update_layout(
        title="Rung Spacing Scores (lower is better)",
        xaxis_title="Spacing",
        yaxis_title="RMS Error",
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Big G Deadlock Scatter (horizontal scatter with error bars)
# ---------------------------------------------------------------------------
def g_deadlock_scatter(measurements_data, G_prediction):
    """Horizontal scatter with error bars showing the Big G deadlock.

    Parameters
    ----------
    measurements_data : list[dict]
        Each dict has keys: experiment (str), G_exp (float), G_unc (float),
        cluster ("high" or "low").
    G_prediction : float
        The vacuum prediction for G.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not measurements_data:
        fig = go.Figure()
        fig.add_annotation(text="No measurement data available", showarrow=False)
        return _apply_theme(fig)

    experiments = [m["experiment"] for m in measurements_data]
    g_values = [float(m["G_exp"]) for m in measurements_data]
    g_unc = [float(m["G_unc"]) for m in measurements_data]
    clusters = [m.get("cluster", "high") for m in measurements_data]

    colors = [
        _COLOR_DEFAULT if c == "high" else _COLOR_HIGHLIGHT_10
        for c in clusters
    ]

    fig = go.Figure()

    # Individual measurement points with error bars
    for i, (exp, gv, gu, col, cl) in enumerate(
        zip(experiments, g_values, g_unc, colors, clusters)
    ):
        fig.add_trace(go.Scatter(
            x=[gv],
            y=[exp],
            error_x=dict(type="data", array=[gu], visible=True, color=col),
            mode="markers",
            marker=dict(color=col, size=10, symbol="circle"),
            name=f"{cl.title()} cluster" if i == 0 or clusters[i] != clusters[i - 1] else None,
            showlegend=(i == 0 or clusters[i] != clusters[i - 1]),
            legendgroup=cl,
            hovertemplate=f"{exp}<br>G = {{x:.5e}} +/- {gu:.2e}<extra></extra>",
        ))

    # Prediction vertical line
    fig.add_vline(
        x=G_prediction, line_dash="solid", line_color=_COLOR_HIGHLIGHT_21,
        line_width=2,
        annotation_text=f"G_vac = {G_prediction:.5e}",
        annotation_position="top right",
        annotation_font_color=_COLOR_HIGHLIGHT_21,
    )

    fig.update_layout(
        title="Big G Deadlock: Experimental Measurements",
        xaxis_title="G (m\u00b3 kg\u207b\u00b9 s\u207b\u00b2)",
        yaxis_title="",
        height=450,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Screening Profile (G_eff vs distance on log x-axis)
# ---------------------------------------------------------------------------
def screening_profile_chart(profile_data, G_codata):
    """Line chart on log x-axis showing G_eff vs distance.

    Parameters
    ----------
    profile_data : dict
        Keys: r_meters (list), G_eff (list), landmarks (list of dicts with
        label (str) and r_meters (float)).
    G_codata : float
        CODATA recommended G value.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not profile_data or not profile_data.get("r_meters"):
        fig = go.Figure()
        fig.add_annotation(text="No screening data available", showarrow=False)
        return _apply_theme(fig)

    r_meters = profile_data["r_meters"]
    g_eff = profile_data["G_eff"]
    landmarks = profile_data.get("landmarks", [])

    G_vacuum = g_eff[-1] if g_eff else 6.67323e-11

    fig = go.Figure()

    # Invisible baseline at G_vacuum for fill reference
    fig.add_trace(go.Scatter(
        x=r_meters, y=[G_vacuum] * len(r_meters),
        mode="lines",
        line=dict(color="rgba(0,0,0,0)", width=0),
        showlegend=False,
        hoverinfo="skip",
    ))

    # G_eff(r) curve, filled down to the G_vacuum baseline
    fig.add_trace(go.Scatter(
        x=r_meters, y=g_eff,
        mode="lines",
        line=dict(color=_COLOR_HIGHLIGHT_21, width=2.5),
        name="G_eff(r)",
        fill="tonexty",
        fillcolor="rgba(245, 158, 11, 0.15)",
    ))

    # G_vacuum horizontal reference
    fig.add_hline(
        y=G_vacuum, line_dash="dash", line_color=_COLOR_GREEN, line_width=1.5,
        annotation_text=f"G_vac = {G_vacuum:.5e}",
        annotation_position="bottom left",
        annotation_font_color=_COLOR_GREEN,
    )

    # G_CODATA horizontal reference
    fig.add_hline(
        y=G_codata, line_dash="dash", line_color=_COLOR_DEFAULT, line_width=1.5,
        annotation_text=f"G_CODATA = {G_codata:.5e}",
        annotation_position="top left",
        annotation_font_color=_COLOR_DEFAULT,
    )

    # Landmark vertical lines
    for lm in landmarks:
        fig.add_vline(
            x=lm["r_meters"], line_dash="dot", line_color="#4c566a", line_width=1,
            annotation_text=lm.get("label", ""),
            annotation_position="top right",
            annotation_font_color="#d8dee9",
            annotation_font_size=10,
        )

    fig.update_layout(
        title="Screening Profile: G_eff vs Distance",
        xaxis_title="Distance (m)",
        yaxis_title="G_eff (m\u00b3 kg\u207b\u00b9 s\u207b\u00b2)",
        xaxis_type="log",
        height=500,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Falsification Bars (sigma deviations per experiment)
# ---------------------------------------------------------------------------
def falsification_bars(comparisons, G_pred):
    """Horizontal bar chart of sigma deviations from prediction.

    Bars are *signed*: experiments that read HIGH (G_exp > G_pred)
    extend to the right (+), while experiments reading LOW extend
    left (-). This makes the directional clustering visible.

    Parameters
    ----------
    comparisons : list[dict]
        Each dict has keys: experiment (str), sigma (float), direction (str).
        direction "+" means G_pred > G_exp (experiment reads LOW).
        direction "-" means G_pred < G_exp (experiment reads HIGH).
    G_pred : float
        The predicted G value (used in title).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not comparisons:
        fig = go.Figure()
        fig.add_annotation(text="No comparison data available", showarrow=False)
        return _apply_theme(fig)

    experiments = [c["experiment"] for c in comparisons]
    raw_sigmas = [abs(float(c["sigma"])) for c in comparisons]
    directions = [c.get("direction", "+") for c in comparisons]

    # Signed sigma: positive = experiment reads HIGH, negative = reads LOW
    signed_sigmas = []
    for s, d in zip(raw_sigmas, directions):
        # direction "-" means G_pred < G_exp => experiment reads HIGH => positive bar
        # direction "+" means G_pred > G_exp => experiment reads LOW  => negative bar
        signed_sigmas.append(s if d == "-" else -s)

    colors = []
    for s, d in zip(raw_sigmas, directions):
        if s < 2:
            colors.append(_COLOR_GREEN)
        elif s <= 5:
            colors.append(_COLOR_HIGHLIGHT_21)
        else:
            colors.append(_COLOR_RED)

    hover_labels = []
    for exp, s, d in zip(experiments, raw_sigmas, directions):
        tag = "HIGH" if d == "-" else "LOW"
        hover_labels.append(f"{exp}<br>{s:.1f} sigma ({tag})<extra></extra>")

    fig = go.Figure()

    fig.add_trace(go.Bar(
        y=experiments,
        x=signed_sigmas,
        orientation="h",
        marker_color=colors,
        hovertemplate=hover_labels,
    ))

    # Zero line (= prediction)
    fig.add_vline(
        x=0, line_color="#e0e0e0", line_width=1.5,
    )

    # Reference lines at +/- 2-sigma and 5-sigma
    for ref in [2.0, 5.0]:
        fig.add_vline(
            x=ref, line_dash="dash", line_color="#4c566a", line_width=1,
            annotation_text=f"{ref:.0f}\u03c3", annotation_position="top right",
            annotation_font_color="#9e9e9e", annotation_font_size=10,
        )
        fig.add_vline(
            x=-ref, line_dash="dash", line_color="#4c566a", line_width=1,
        )

    fig.update_layout(
        title=f"Falsification Clock: Sigma Deviations from G_pred = {G_pred:.5e}",
        xaxis_title="\u03c3 (+ = experiment reads HIGH, \u2212 = reads LOW)",
        yaxis_title="",
        height=400,
        showlegend=False,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Genus Scan (required GB coupling per genus)
# ---------------------------------------------------------------------------
def genus_scan_chart(scan_data, omega_target):
    """Vertical bar chart of required |lambda_GB| per genus.

    Parameters
    ----------
    scan_data : list[dict]
        Each dict from scan_golden_point() with keys: genus, required_gb_coupling.
    omega_target : float
        Target omega value (for annotation).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not scan_data:
        fig = go.Figure()
        fig.add_annotation(text="No scan data available", showarrow=False)
        return _apply_theme(fig)

    labels = [f"g={r['genus']}" for r in scan_data]
    values = [abs(r["required_gb_coupling"]) for r in scan_data]
    colors = [
        _COLOR_HIGHLIGHT_21 if r["genus"] == 2 else _COLOR_GREEN
        for r in scan_data
    ]

    fig = go.Figure(
        data=go.Bar(
            x=labels,
            y=values,
            marker_color=colors,
            hovertemplate="Genus %{x}<br>|λ_GB| = %{y:.6f}<extra></extra>",
        )
    )

    fig.add_hline(
        y=1.0, line_dash="dash", line_color="#e0e0e0", line_width=1,
        annotation_text="Natural coupling threshold",
        annotation_position="top right",
        annotation_font_color="#e0e0e0",
    )

    fig.update_layout(
        title="Required GB Coupling by Genus",
        xaxis_title="Genus",
        yaxis_title="|λ_GB|",
        height=400,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Algebraic Chain (waterfall from baseline to target omega)
# ---------------------------------------------------------------------------
def algebraic_chain_chart(chain_data):
    """Waterfall chart showing omega derivation from baseline to target.

    Parameters
    ----------
    chain_data : dict
        Keys: omega_baseline, gap, omega_target.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not chain_data:
        fig = go.Figure()
        fig.add_annotation(text="No chain data available", showarrow=False)
        return _apply_theme(fig)

    fig = go.Figure(
        data=go.Waterfall(
            x=["EH Baseline", "GB Shift", "Target omega"],
            y=[chain_data["omega_baseline"], -chain_data["gap"], chain_data["omega_target"]],
            measure=["absolute", "relative", "total"],
            connector=dict(line=dict(color="#4c566a")),
            increasing=dict(marker=dict(color=_COLOR_DEFAULT)),
            decreasing=dict(marker=dict(color=_COLOR_HIGHLIGHT_21)),
            totals=dict(marker=dict(color=_COLOR_GREEN)),
            hovertemplate="%{x}<br>%{y:.8f}<extra></extra>",
        )
    )

    fig.add_hline(
        y=chain_data["omega_target"], line_dash="dash",
        line_color=_COLOR_GREEN, line_width=1,
        annotation_text=f"ω = {chain_data['omega_target']:.8f}",
        annotation_position="top right",
        annotation_font_color=_COLOR_GREEN,
    )

    fig.update_layout(
        title="From Baseline ω to Target",
        yaxis_title="ω",
        height=400,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Discriminant Scan (which dimensions produce Q(sqrt(5)))
# ---------------------------------------------------------------------------
def discriminant_scan_chart(scan_results):
    """Horizontal bar chart showing which internal dimensions n produce Q(sqrt(5)).

    Parameters
    ----------
    scan_results : list[dict]
        Each dict from scan_discriminant_field() with keys: n (int),
        D (int), discriminant (int), squarefree_part (int),
        is_golden_field (bool).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not scan_results:
        fig = go.Figure()
        fig.add_annotation(text="No scan results available", showarrow=False)
        return _apply_theme(fig)

    y_labels = [f"n={r['n']} (D={r['D']})" for r in scan_results]
    x_values = [r["squarefree_part"] for r in scan_results]
    colors = [
        _COLOR_HIGHLIGHT_21 if r["is_golden_field"] else _COLOR_DEFAULT
        for r in scan_results
    ]

    fig = go.Figure(
        data=go.Bar(
            x=x_values,
            y=y_labels,
            orientation="h",
            marker_color=colors,
            hovertemplate=(
                "n=%{customdata[0]}, D=%{customdata[1]}<br>"
                "Discriminant: %{customdata[2]}<br>"
                "Squarefree part: %{x}<extra></extra>"
            ),
            customdata=[
                [r["n"], r["D"], r["discriminant"]] for r in scan_results
            ],
        )
    )

    # Vertical reference line at squarefree_part = 5
    fig.add_vline(
        x=5, line_dash="dash", line_color=_COLOR_HIGHLIGHT_21, line_width=1.5,
        annotation_text="√5 field",
        annotation_position="top right",
        annotation_font_color=_COLOR_HIGHLIGHT_21,
    )

    fig.update_layout(
        title="Discriminant Scan: Which Dimensions Produce ℚ(√5)?",
        xaxis_title="Squarefree Part of Discriminant",
        yaxis_title="",
        height=400,
        showlegend=False,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Moduli Mass Spectrum (vertical bar chart of moduli masses)
# ---------------------------------------------------------------------------
def moduli_mass_chart(mass_spectrum):
    """Vertical bar chart of moduli masses.

    Parameters
    ----------
    mass_spectrum : dict
        From compute_moduli_mass_spectrum() with keys: volume_modulus
        (dict with mass_squared or None), shape_moduli (list of dicts
        with mass_squared, index).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not mass_spectrum:
        fig = go.Figure()
        fig.add_annotation(text="No mass spectrum available", showarrow=False)
        return _apply_theme(fig)

    x_labels = []
    y_values = []
    colors = []

    volume = mass_spectrum.get("volume_modulus")
    if volume is not None and volume.get("mass_squared") is not None:
        x_labels.append("Volume (\u03c3)")
        y_values.append(float(volume["mass_squared"]))
        colors.append(_COLOR_HIGHLIGHT_21)

    for sm in mass_spectrum.get("shape_moduli", []):
        idx = sm.get("index", len(x_labels))
        x_labels.append(f"Shape {idx}")
        y_values.append(float(sm["mass_squared"]))
        colors.append(_COLOR_GREEN)

    if not x_labels:
        fig = go.Figure()
        fig.add_annotation(text="No mass data to display", showarrow=False)
        return _apply_theme(fig)

    fig = go.Figure(
        data=go.Bar(
            x=x_labels,
            y=y_values,
            marker_color=colors,
            hovertemplate="Modulus: %{x}<br>m\u00b2 = %{y:.6e}<extra></extra>",
        )
    )

    # Stability threshold at y = 0
    fig.add_hline(
        y=0, line_dash="dash", line_color="#e0e0e0", line_width=1,
        annotation_text="Stability threshold",
        annotation_position="top right",
        annotation_font_color="#e0e0e0",
    )

    fig.update_layout(
        title="Moduli Mass Spectrum",
        xaxis_title="Modulus",
        yaxis_title="m\u00b2 (mass squared)",
        height=400,
        showlegend=False,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Algebraic Field Table (constants in Q(sqrt(5)))
# ---------------------------------------------------------------------------
def algebraic_field_chart(closure_data):
    """Table visualization showing all constants in Q(sqrt(5)).

    Parameters
    ----------
    closure_data : dict
        From derive_algebraic_closure() with key constants (list of dicts
        with name, value, exact_form, minimal_polynomial, in_field).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not closure_data or not closure_data.get("constants"):
        fig = go.Figure()
        fig.add_annotation(text="No algebraic closure data available", showarrow=False)
        return _apply_theme(fig)

    constants = closure_data["constants"]

    names = [c["name"] for c in constants]
    exact_forms = [c["exact_form"] for c in constants]
    min_polys = [c["minimal_polynomial"] for c in constants]
    values = [f"{float(c['value']):.10g}" for c in constants]

    header_bg = "#1a1d23"
    cell_bg = "#0e1117"

    fig = go.Figure(
        data=go.Table(
            header=dict(
                values=["Name", "Exact Form", "Minimal Polynomial", "Value"],
                fill_color=header_bg,
                font=dict(
                    family="Fira Mono, Consolas, monospace",
                    color="#ffffff",
                    size=13,
                ),
                align="left",
                line_color="#262b33",
            ),
            cells=dict(
                values=[names, exact_forms, min_polys, values],
                fill_color=cell_bg,
                font=dict(
                    family="Fira Mono, Consolas, monospace",
                    color=[
                        "#ffffff",       # Name
                        "#e0e0e0",       # Exact Form
                        "#e0e0e0",       # Minimal Polynomial
                        _COLOR_GREEN,    # Value
                    ],
                    size=12,
                ),
                align="left",
                line_color="#262b33",
            ),
        )
    )

    fig.update_layout(
        title="Algebraic Closure: All Constants in ℚ(√5)",
        title_font=dict(family="Crimson Pro, Georgia, serif", size=16, color="#ffffff"),
        paper_bgcolor="#0e1117",
        margin=dict(l=20, r=20, t=50, b=20),
        height=350,
    )

    return fig


# ---------------------------------------------------------------------------
# Chart: PPN Profile (|gamma_PPN - 1| vs distance)
# ---------------------------------------------------------------------------
def ppn_profile_chart(ppn_profile_data):
    """Line chart of |gamma_PPN - 1| vs distance on log-x axis.

    Parameters
    ----------
    ppn_profile_data : dict
        From compute_ppn_profile() with keys: r_meters, gamma_deviation, landmarks.
    """
    if not ppn_profile_data or not ppn_profile_data.get("r_meters"):
        fig = go.Figure()
        fig.add_annotation(text="No PPN profile data available", showarrow=False)
        return _apply_theme(fig)

    r_meters = ppn_profile_data["r_meters"]
    gamma_dev = ppn_profile_data["gamma_deviation"]
    landmarks = ppn_profile_data.get("landmarks", {})
    cassini_bound = 2.3e-5

    fig = go.Figure()

    # Main curve
    fig.add_trace(go.Scatter(
        x=r_meters, y=gamma_dev,
        mode="lines",
        line=dict(color=_COLOR_HIGHLIGHT_21, width=2.5),
        name="|γ − 1|",
        hovertemplate="r = %{x:.2e} m<br>|γ−1| = %{y:.2e}<extra></extra>",
    ))

    # Cassini bound horizontal line
    fig.add_hline(
        y=cassini_bound, line_dash="dash", line_color=_COLOR_RED, line_width=1.5,
        annotation_text=f"Cassini bound ({cassini_bound:.1e})",
        annotation_position="top right",
        annotation_font_color=_COLOR_RED,
    )

    # Shaded danger zone
    fig.add_hrect(
        y0=cassini_bound, y1=1.0,
        fillcolor="rgba(248, 113, 113, 0.08)",
        line_width=0,
    )

    # Landmark vertical lines
    for key, lm in landmarks.items():
        fig.add_vline(
            x=lm["r_meters"], line_dash="dot", line_color="#4c566a", line_width=1,
            annotation_text=lm.get("label", key),
            annotation_position="top right",
            annotation_font_color="#d8dee9",
            annotation_font_size=10,
        )

    fig.update_layout(
        title="PPN Profile: |γ − 1| vs Distance",
        xaxis_title="Distance (m)",
        yaxis_title="|γ_PPN − 1|",
        xaxis_type="log",
        yaxis_type="log",
        height=400,
        showlegend=True,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Fifth-Force Exclusion Plot
# ---------------------------------------------------------------------------
def fifth_force_exclusion_chart(dilaton_point, bounds):
    """Log-log plot of (lambda, alpha) exclusion space.

    Parameters
    ----------
    dilaton_point : dict
        With keys: log10_alpha, log10_lambda_m.
    bounds : list[dict]
        Each with keys: name, boundary (list of [log10_lambda, log10_alpha]).
    """
    fig = go.Figure()

    # Exclusion region colors
    region_colors = [
        "rgba(96, 165, 250, 0.15)",   # blue
        "rgba(167, 139, 250, 0.15)",   # violet
        "rgba(52, 211, 153, 0.15)",    # green
        "rgba(248, 113, 113, 0.15)",   # red
        "rgba(245, 158, 11, 0.15)",    # amber
    ]
    line_colors = [
        _COLOR_DEFAULT, _COLOR_HIGHLIGHT_10, _COLOR_GREEN,
        _COLOR_RED, _COLOR_HIGHLIGHT_21,
    ]

    for i, bound in enumerate(bounds):
        boundary = bound["boundary"]
        x_vals = [b[0] for b in boundary]
        y_vals = [b[1] for b in boundary]
        color_idx = i % len(region_colors)

        fig.add_trace(go.Scatter(
            x=x_vals, y=y_vals,
            mode="lines",
            line=dict(color=line_colors[color_idx], width=1.5),
            fill="tozeroy",
            fillcolor=region_colors[color_idx],
            name=bound["name"],
            hovertemplate=f"{bound['name']}<br>log₁₀(λ) = %{{x:.1f}}<br>log₁₀(α) = %{{y:.1f}}<extra></extra>",
        ))

    # Dilaton point
    if dilaton_point:
        fig.add_trace(go.Scatter(
            x=[dilaton_point["log10_lambda_m"]],
            y=[dilaton_point["log10_alpha"]],
            mode="markers",
            marker=dict(
                color=_COLOR_HIGHLIGHT_21, size=16,
                symbol="diamond", line=dict(color="#ffffff", width=2),
            ),
            name="Dilaton",
            hovertemplate=(
                f"Dilaton<br>"
                f"log₁₀(λ) = {dilaton_point['log10_lambda_m']:.2f}<br>"
                f"log₁₀(α) = {dilaton_point['log10_alpha']:.2f}<extra></extra>"
            ),
        ))

    fig.update_layout(
        title="Fifth-Force Exclusion Plot",
        xaxis_title="log₁₀(λ / m)",
        yaxis_title="log₁₀(α)",
        height=450,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Screening Discrepancy (tree-level vs empirical)
# ---------------------------------------------------------------------------
def screening_discrepancy_chart(discrepancy_data):
    """Horizontal bar chart comparing tree-level vs empirical screening amplitude.

    Parameters
    ----------
    discrepancy_data : dict
        With keys: alpha_tree, alpha_empirical, ratio.
    """
    if not discrepancy_data:
        fig = go.Figure()
        fig.add_annotation(text="No discrepancy data available", showarrow=False)
        return _apply_theme(fig)

    import math as _m
    alpha_tree = discrepancy_data["alpha_tree"]
    alpha_empirical = discrepancy_data["alpha_empirical"]
    ratio = discrepancy_data["ratio"]

    fig = go.Figure()

    fig.add_trace(go.Bar(
        y=["Tree-level (2 beta^2)"],
        x=[_m.log10(alpha_tree) if alpha_tree > 0 else 0],
        orientation="h",
        marker_color=_COLOR_HIGHLIGHT_21,
        name="Tree-level",
        hovertemplate=f"α_tree = {alpha_tree:.4e}<extra></extra>",
    ))

    fig.add_trace(go.Bar(
        y=["Empirical (G_lab - G_vac)/G_vac"],
        x=[_m.log10(abs(alpha_empirical)) if alpha_empirical != 0 else 0],
        orientation="h",
        marker_color=_COLOR_GREEN,
        name="Empirical",
        hovertemplate=f"α_empirical = {alpha_empirical:.4e}<extra></extra>",
    ))

    fig.add_annotation(
        x=0, y=1.15,
        xref="paper", yref="paper",
        text=f"Ratio: {ratio:.0f}x",
        showarrow=False,
        font=dict(size=14, color=_COLOR_HIGHLIGHT_21),
    )

    fig.update_layout(
        title="Screening Amplitude: Tree-level vs Empirical",
        xaxis_title="log₁₀(α)",
        yaxis_title="",
        height=300,
        showlegend=False,
        barmode="group",
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Literature Comparison Table
# ---------------------------------------------------------------------------
def literature_comparison_chart(literature_data):
    """Plotly table comparing the framework to 4 literature models.

    Parameters
    ----------
    literature_data : dict
        With key 'comparisons', a list of dicts with keys:
        framework, similarities, differences, key_reference.
    """
    if not literature_data or not literature_data.get("comparisons"):
        fig = go.Figure()
        fig.add_annotation(text="No literature data available", showarrow=False)
        return _apply_theme(fig)

    comparisons = literature_data["comparisons"]

    frameworks = [c["framework"] for c in comparisons]
    similarities = ["<br>".join(f"- {s}" for s in c["similarities"]) for c in comparisons]
    differences = ["<br>".join(f"- {d}" for d in c["differences"]) for c in comparisons]
    references = [c["key_reference"] for c in comparisons]

    header_bg = "#1a1d23"
    cell_bg = "#0e1117"

    fig = go.Figure(
        data=go.Table(
            header=dict(
                values=["Framework", "Similarities", "Differences", "Key Reference"],
                fill_color=header_bg,
                font=dict(family="Fira Mono, Consolas, monospace", color="#ffffff", size=12),
                align="left",
                line_color="#262b33",
            ),
            cells=dict(
                values=[frameworks, similarities, differences, references],
                fill_color=cell_bg,
                font=dict(family="Fira Mono, Consolas, monospace", color="#e0e0e0", size=11),
                align="left",
                line_color="#262b33",
                height=40,
            ),
        )
    )

    fig.update_layout(
        title="Literature Comparison",
        title_font=dict(family="Crimson Pro, Georgia, serif", size=16, color="#ffffff"),
        paper_bgcolor="#0e1117",
        margin=dict(l=20, r=20, t=50, b=20),
        height=350,
    )

    return fig


# ---------------------------------------------------------------------------
# Chart: Casimir Effective Potential (V_eff vs sigma)
# ---------------------------------------------------------------------------
def casimir_potential_chart(potential_data):
    """Line chart of V_eff(sigma) with classical, Casimir, and total curves.

    Parameters
    ----------
    potential_data : dict
        From compute_effective_potential() with keys: sigma_grid,
        V_classical, V_casimir, V_total.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not potential_data or not potential_data.get("sigma_grid"):
        fig = go.Figure()
        fig.add_annotation(text="No potential data available", showarrow=False)
        return _apply_theme(fig)

    sigma = potential_data["sigma_grid"]
    V_class = potential_data["V_classical"]
    V_cas = potential_data["V_casimir"]
    V_total = potential_data["V_total"]

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=sigma, y=V_class,
        mode="lines",
        line=dict(color=_COLOR_DEFAULT, width=2, dash="dash"),
        name="V_classical (curvature)",
        hovertemplate="σ = %{x:.2f}<br>V_class = %{y:.4e}<extra></extra>",
    ))

    fig.add_trace(go.Scatter(
        x=sigma, y=V_cas,
        mode="lines",
        line=dict(color=_COLOR_HIGHLIGHT_10, width=2, dash="dot"),
        name="V_Casimir (1-loop)",
        hovertemplate="σ = %{x:.2f}<br>V_Casimir = %{y:.4e}<extra></extra>",
    ))

    fig.add_trace(go.Scatter(
        x=sigma, y=V_total,
        mode="lines",
        line=dict(color=_COLOR_HIGHLIGHT_21, width=3),
        name="V_total",
        hovertemplate="σ = %{x:.2f}<br>V_total = %{y:.4e}<extra></extra>",
    ))

    # Zero line
    fig.add_hline(y=0, line_dash="solid", line_color="#4c566a", line_width=1)

    fig.update_layout(
        title="Effective Potential: Casimir Stabilization",
        xaxis_title="σ (breathing mode)",
        yaxis_title="V_eff (Planck units)",
        height=450,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: KK Spectrum (bar chart of KK masses with degeneracies)
# ---------------------------------------------------------------------------
def kk_spectrum_chart(spectrum_data):
    """Bar chart of KK eigenvalues with degeneracy as bar width/color.

    Parameters
    ----------
    spectrum_data : dict
        From compute_kk_spectrum_s2() with keys: eigenvalues (list),
        degeneracies (list).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not spectrum_data or not spectrum_data.get("eigenvalues"):
        fig = go.Figure()
        fig.add_annotation(text="No spectrum data available", showarrow=False)
        return _apply_theme(fig)

    eigenvalues = spectrum_data["eigenvalues"]
    degeneracies = spectrum_data["degeneracies"]

    # Show first 20 modes for clarity
    n_show = min(20, len(eigenvalues))
    l_values = list(range(1, n_show + 1))
    ev_show = eigenvalues[:n_show]
    deg_show = degeneracies[:n_show]

    # Color by degeneracy
    max_deg = max(deg_show) if deg_show else 1
    colors = [
        f"rgba(96, 165, 250, {0.3 + 0.7 * d / max_deg})"
        for d in deg_show
    ]

    fig = go.Figure(
        data=go.Bar(
            x=[f"l={l}" for l in l_values],
            y=ev_show,
            marker_color=colors,
            text=[f"d={d}" for d in deg_show],
            textposition="auto",
            hovertemplate="l=%{customdata}<br>m² = l(l+1) = %{y}<br>degeneracy = %{text}<extra></extra>",
            customdata=l_values,
        )
    )

    fig.update_layout(
        title="KK Spectrum on S²: Eigenvalues and Degeneracies",
        xaxis_title="Angular momentum l",
        yaxis_title="m² = l(l+1) / a²",
        height=400,
        showlegend=False,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Flux Potential (three-term potential: Casimir + Curvature + Flux)
# ---------------------------------------------------------------------------
def flux_potential_chart(potential_data, min_location=None):
    """Line chart of the three-term potential V(sigma).

    Parameters
    ----------
    potential_data : dict
        Keys: sigma_grid (list), V_casimir (list), V_curvature (list),
        V_flux (list), V_total (list).
    min_location : float or None
        If provided, marks the minimum with a diamond marker.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not potential_data or not potential_data.get("sigma_grid"):
        fig = go.Figure()
        fig.add_annotation(text="No potential data available", showarrow=False)
        return _apply_theme(fig)

    sigma = potential_data["sigma_grid"]
    V_casimir = potential_data["V_casimir"]
    V_curvature = potential_data["V_curvature"]
    V_flux = potential_data["V_flux"]
    V_total = potential_data["V_total"]

    fig = go.Figure()

    # Casimir term (blue dashed)
    fig.add_trace(go.Scatter(
        x=sigma, y=V_casimir,
        mode="lines",
        line=dict(color=_COLOR_DEFAULT, width=2, dash="dash"),
        name="V_Casimir (A e^4s)",
        hovertemplate="σ = %{x:.2f}<br>V_Casimir = %{y:.4e}<extra></extra>",
    ))

    # Curvature term (violet dotted)
    fig.add_trace(go.Scatter(
        x=sigma, y=V_curvature,
        mode="lines",
        line=dict(color=_COLOR_HIGHLIGHT_10, width=2, dash="dot"),
        name="V_curvature (B e^2s)",
        hovertemplate="σ = %{x:.2f}<br>V_curvature = %{y:.4e}<extra></extra>",
    ))

    # Flux term (green dash-dot)
    fig.add_trace(go.Scatter(
        x=sigma, y=V_flux,
        mode="lines",
        line=dict(color=_COLOR_GREEN, width=2, dash="dashdot"),
        name="V_flux (C e^6s)",
        hovertemplate="σ = %{x:.2f}<br>V_flux = %{y:.4e}<extra></extra>",
    ))

    # Total potential (amber solid thick)
    fig.add_trace(go.Scatter(
        x=sigma, y=V_total,
        mode="lines",
        line=dict(color=_COLOR_HIGHLIGHT_21, width=3),
        name="V_total",
        hovertemplate="σ = %{x:.2f}<br>V_total = %{y:.4e}<extra></extra>",
    ))

    # Mark minimum with diamond marker
    if min_location is not None:
        # Find closest sigma index to the minimum location
        import bisect
        idx = bisect.bisect_left(sigma, min_location)
        if idx >= len(sigma):
            idx = len(sigma) - 1
        V_at_min = V_total[idx]

        fig.add_trace(go.Scatter(
            x=[min_location],
            y=[V_at_min],
            mode="markers",
            marker=dict(
                color=_COLOR_HIGHLIGHT_21, size=14,
                symbol="diamond", line=dict(color="#ffffff", width=2),
            ),
            name="Minimum",
            hovertemplate=(
                f"σ₀ = {min_location:.4f}<br>"
                f"V_min = {V_at_min:.4e}<extra></extra>"
            ),
        ))

    # Zero line
    fig.add_hline(y=0, line_dash="solid", line_color="#4c566a", line_width=1)

    fig.update_layout(
        title="Three-Term Potential: Casimir + Curvature + Flux",
        xaxis_title="σ (breathing mode)",
        yaxis_title="V_eff (Planck units)",
        height=450,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Radius Mass (dilaton mass vs internal radius, log-log)
# ---------------------------------------------------------------------------
def radius_mass_chart(mass_data, testable_window=None):
    """Log-log line chart of m_phi (eV) vs a_0 (m).

    Parameters
    ----------
    mass_data : dict
        Keys: a_0_values (list of float, meters), m_phi_eV (list of float, eV).
    testable_window : dict or None
        If provided, dict with keys a_0_min, a_0_max_testable (in meters).
        A shaded vertical band marks the testable region.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not mass_data or not mass_data.get("a_0_values"):
        fig = go.Figure()
        fig.add_annotation(text="No mass data available", showarrow=False)
        return _apply_theme(fig)

    a_0 = mass_data["a_0_values"]
    m_phi = mass_data["m_phi_eV"]

    fig = go.Figure()

    # Main curve
    fig.add_trace(go.Scatter(
        x=a_0, y=m_phi,
        mode="lines",
        line=dict(color=_COLOR_HIGHLIGHT_21, width=2.5),
        name="m_φ(a₀)",
        hovertemplate="a₀ = %{x:.2e} m<br>m_φ = %{y:.2e} eV<extra></extra>",
    ))

    # Testable window shaded band
    if testable_window:
        a_min = testable_window.get("a_0_min")
        a_max = testable_window.get("a_0_max_testable")
        if a_min is not None and a_max is not None:
            fig.add_vrect(
                x0=a_min, x1=a_max,
                fillcolor="rgba(52, 211, 153, 0.12)",
                line_width=0,
                annotation_text="Testable window",
                annotation_position="top left",
                annotation_font_color=_COLOR_GREEN,
            )

    # Horizontal line: 2 meV threshold (amber)
    fig.add_hline(
        y=2e-3, line_dash="dash", line_color=_COLOR_HIGHLIGHT_21, line_width=1.5,
        annotation_text="2 meV threshold",
        annotation_position="bottom right",
        annotation_font_color=_COLOR_HIGHLIGHT_21,
    )

    # Horizontal line: M_Pl (red dashed)
    fig.add_hline(
        y=1.22e28, line_dash="dash", line_color=_COLOR_RED, line_width=1.5,
        annotation_text="M_Pl",
        annotation_position="top right",
        annotation_font_color=_COLOR_RED,
    )

    # Vertical line: l_Pl (red)
    fig.add_vline(
        x=1.616e-35, line_dash="solid", line_color=_COLOR_RED, line_width=1.5,
        annotation_text="l_Pl",
        annotation_position="bottom right",
        annotation_font_color=_COLOR_RED,
    )

    # Vertical line: Eot-Wash 30 um bound (green)
    fig.add_vline(
        x=30e-6, line_dash="solid", line_color=_COLOR_GREEN, line_width=1.5,
        annotation_text="Eot-Wash 30 µm",
        annotation_position="top left",
        annotation_font_color=_COLOR_GREEN,
    )

    fig.update_layout(
        title="Dilaton Mass vs Internal Radius",
        xaxis_title="a₀ (m)",
        yaxis_title="m_φ (eV)",
        xaxis_type="log",
        yaxis_type="log",
        height=500,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Flux Scan (dilaton mass vs flux quantum N)
# ---------------------------------------------------------------------------
def flux_scan_chart(scan_data):
    """Bar chart of m_phi_eV vs flux quantum N.

    Parameters
    ----------
    scan_data : dict
        Key: results (list of dicts with N (int) and m_phi_eV (float)).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not scan_data or not scan_data.get("results"):
        fig = go.Figure()
        fig.add_annotation(text="No scan data available", showarrow=False)
        return _apply_theme(fig)

    results = scan_data["results"]
    n_values = [r["N"] for r in results]
    m_values = [float(r["m_phi_eV"]) for r in results]

    fig = go.Figure(
        data=go.Bar(
            x=[f"N={n}" for n in n_values],
            y=m_values,
            marker_color=_COLOR_HIGHLIGHT_21,
            hovertemplate="N = %{customdata}<br>m_φ = %{y:.4e} eV<extra></extra>",
            customdata=n_values,
        )
    )

    fig.update_layout(
        title="Dilaton Mass vs Flux Quantum N",
        xaxis_title="Flux Quantum N",
        yaxis_title="m_φ (eV)",
        yaxis_type="log",
        height=400,
        showlegend=False,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Chameleon Profile (dilaton mass vs matter density)
# ---------------------------------------------------------------------------
def chameleon_profile_chart(profile_data):
    """Log-log line chart of effective dilaton mass vs matter density.

    Parameters
    ----------
    profile_data : dict
        From compute_chameleon_profile() with keys: rho_values_si,
        m_phi_eff_eV_values.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import math as _m

    rho_si = profile_data.get("rho_values_si", [])
    m_phi = profile_data.get("m_phi_eff_eV_values", [])

    # Filter out None / non-positive values
    valid = [
        (r, m) for r, m in zip(rho_si, m_phi)
        if m is not None and r is not None and m > 0 and r > 0
    ]

    fig = go.Figure()

    if not valid:
        fig.add_annotation(text="No valid profile data available", showarrow=False)
        return _apply_theme(fig)

    rho_valid, m_valid = zip(*valid)

    log_rho = [_m.log10(r) for r in rho_valid]
    log_m = [_m.log10(m) for m in m_valid]

    fig.add_trace(go.Scatter(
        x=log_rho, y=log_m,
        mode="lines",
        line=dict(color=_COLOR_RED, width=2.5),
        name="m_φ(ρ)",
        hovertemplate=(
            "log₁₀(ρ) = %{x:.1f} kg/m³<br>"
            "log₁₀(m_φ) = %{y:.1f} eV<extra></extra>"
        ),
    ))

    fig.update_layout(
        title="Chameleon Profile: Dilaton Mass vs Density",
        xaxis_title="log₁₀(ρ / kg m⁻³)",
        yaxis_title="log₁₀(m_φ / eV)",
        height=400,
        showlegend=True,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: KK Truncation (effective internal radius vs density)
# ---------------------------------------------------------------------------
def kk_truncation_chart(profile_data):
    """Log-log line chart of effective internal radius vs density
    with Eot-Wash exclusion line.

    Parameters
    ----------
    profile_data : dict
        From compute_chameleon_profile() with keys: rho_values_si,
        a_eff_values.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import math as _m

    rho_si = profile_data.get("rho_values_si", [])
    a_eff = profile_data.get("a_eff_values", [])

    valid = [
        (r, a) for r, a in zip(rho_si, a_eff)
        if a is not None and r is not None and a > 0 and r > 0
    ]

    fig = go.Figure()

    if not valid:
        fig.add_annotation(
            text="No valid KK truncation data available", showarrow=False,
        )
        return _apply_theme(fig)

    rho_valid, a_valid = zip(*valid)

    log_rho = [_m.log10(r) for r in rho_valid]
    log_a = [_m.log10(a) for a in a_valid]

    # Effective radius curve
    fig.add_trace(go.Scatter(
        x=log_rho, y=log_a,
        mode="lines",
        line=dict(color=_COLOR_DEFAULT, width=2.5),
        name="a_eff(ρ)",
        hovertemplate=(
            "log₁₀(ρ) = %{x:.1f} kg/m³<br>"
            "log₁₀(a_eff) = %{y:.2f} m<extra></extra>"
        ),
    ))

    # Eot-Wash exclusion line at 56 um
    eot_wash_log = _m.log10(56e-6)
    fig.add_hline(
        y=eot_wash_log,
        line_dash="dash",
        line_color=_COLOR_RED,
        line_width=2,
        annotation_text="Eot-Wash 56 µm",
        annotation_position="bottom right",
        annotation_font_color=_COLOR_RED,
    )

    fig.update_layout(
        title="Effective Internal Radius vs Density (red: Eot-Wash limit)",
        xaxis_title="log₁₀(ρ / kg m⁻³)",
        yaxis_title="log₁₀(a_eff / m)",
        height=400,
        showlegend=True,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Dark Sector Landscape (bar chart by mass scale)
# ---------------------------------------------------------------------------
def dark_sector_landscape_chart(landscape_data):
    """Horizontal bar chart of the dark sector landscape classifications.

    Parameters
    ----------
    landscape_data : dict
        With key 'landscape', a list of dicts each containing
        m_phi_eV, classification, kk_truncation_valid.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import math as _m

    entries = landscape_data.get("landscape", [])

    if not entries:
        fig = go.Figure()
        fig.add_annotation(text="No landscape data available", showarrow=False)
        return _apply_theme(fig)

    log_masses = [
        _m.log10(e["m_phi_eV"]) if e.get("m_phi_eV", 0) > 0 else 0
        for e in entries
    ]
    classifications = [e.get("classification", "unknown") for e in entries]

    color_map = {
        "dark_energy_candidate": "#a78bfa",   # violet
        "fuzzy_dm_candidate": _COLOR_DEFAULT,  # blue
        "kk_excluded": _COLOR_RED,             # red
        "sub_mm_testable": _COLOR_GREEN,       # green
        "thermal_relic_range": "#f59e0b",      # amber
        "invisible_planck": "#6b7280",         # gray
    }

    colors = [color_map.get(c, "#6b7280") for c in classifications]
    labels = [c.replace("_", " ").title() for c in classifications]

    fig = go.Figure(
        data=go.Bar(
            x=log_masses,
            y=labels,
            orientation="h",
            marker_color=colors,
            hovertemplate=(
                "log₁₀(m_φ) = %{x:.0f} eV<br>"
                "%{y}<extra></extra>"
            ),
        )
    )

    fig.update_layout(
        title="Dark Sector Landscape by Mass Scale",
        xaxis_title="log₁₀(m_φ / eV)",
        yaxis_title="",
        height=max(300, 40 * len(entries)),
        showlegend=False,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Exclusion Plot (Alpha Ladder prediction vs experimental bounds)
# ---------------------------------------------------------------------------
def exclusion_plot_chart(exclusion_data):
    """Plot the (lambda, alpha) exclusion map with the Alpha Ladder prediction line.

    Parameters
    ----------
    exclusion_data : dict
        From compute_exclusion_map() with keys: experiments, alpha_ladder_line,
        survival_window_m.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import math as _m

    fig = go.Figure()

    # Alpha Ladder prediction line (horizontal at 0.618)
    al_alpha = exclusion_data.get("alpha_ladder_line", 0.618)
    log_alpha = _m.log10(al_alpha)

    fig.add_trace(go.Scatter(
        x=[-7, -1],
        y=[log_alpha, log_alpha],
        mode="lines",
        name=f"Alpha Ladder (α = {al_alpha})",
        line=dict(color=_COLOR_RED, width=3, dash="solid"),
        hovertemplate=f"α = {al_alpha}<extra></extra>",
    ))

    # Experiment bounds (approximate regions)
    experiments = exclusion_data.get("experiments", [])
    exp_colors = [
        _COLOR_DEFAULT, _COLOR_GREEN, _COLOR_HIGHLIGHT_10,
        _COLOR_HIGHLIGHT_21, "#1abc9c", "#e67e22", "#95a5a6",
    ]

    for i, exp in enumerate(experiments):
        lam_range = exp["lambda_range_m"]
        alpha_bound = exp["alpha_bound_at_best"]
        if alpha_bound > 0 and alpha_bound < 1e10:
            log_lam_min = _m.log10(lam_range[0])
            log_lam_max = _m.log10(lam_range[1])
            log_ab = _m.log10(alpha_bound) if alpha_bound > 0 else -15
            color = exp_colors[i % len(exp_colors)]
            fig.add_trace(go.Scatter(
                x=[(log_lam_min + log_lam_max) / 2],
                y=[log_ab],
                mode="markers+text",
                name=exp["name"],
                text=[exp["name"]],
                textposition="top center",
                marker=dict(size=10, color=color),
                textfont=dict(size=9),
                hovertemplate=(
                    f"{exp['name']}<br>"
                    f"log₁₀(α_bound) = {log_ab:.1f}<extra></extra>"
                ),
            ))

    # Survival window shading
    sw = exclusion_data.get("survival_window_m")
    if sw:
        fig.add_vrect(
            x0=_m.log10(sw[0]), x1=_m.log10(sw[1]),
            fillcolor="rgba(52, 211, 153, 0.15)", line_width=0,
            annotation_text="Survival Window",
            annotation_position="top left",
            annotation_font_color=_COLOR_GREEN,
        )

    fig.update_layout(
        title="Exclusion Plot: Alpha Ladder Prediction vs Experimental Bounds",
        xaxis_title="log₁₀(λ / m)",
        yaxis_title="log₁₀(α)",
        height=500,
        showlegend=True,
        legend=dict(font=dict(size=10)),
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Signal vs Distance (Yukawa fractional deviation)
# ---------------------------------------------------------------------------
def signal_vs_distance_chart(signal_data):
    """Plot the Yukawa fractional deviation as a function of distance.

    Parameters
    ----------
    signal_data : dict
        From compute_signal_vs_distance() with keys: r_meters,
        fractional_deviation, lambda_m, peak_distance, peak_signal.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import math as _m

    r_values = signal_data.get("r_meters", [])
    deviations = signal_data.get("fractional_deviation", [])

    if not r_values or not deviations:
        fig = go.Figure()
        fig.add_annotation(text="No signal data available", showarrow=False)
        return _apply_theme(fig)

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=[_m.log10(r) for r in r_values],
        y=[d * 100 for d in deviations],
        mode="lines",
        name="Yukawa signal",
        line=dict(color=_COLOR_RED, width=2),
        hovertemplate="log₁₀(r) = %{x:.2f}<br>signal = %{y:.2f}%<extra></extra>",
    ))

    # Sensitivity line at 1%
    fig.add_hline(
        y=1.0, line_dash="dash", line_color="#ffffff", line_width=1.5,
        annotation_text="Eot-Wash sensitivity (1%)",
        annotation_position="top right",
        annotation_font_color="#ffffff",
    )

    # Mark the peak
    peak_r = signal_data.get("peak_distance")
    peak_s = signal_data.get("peak_signal")
    if peak_r and peak_s:
        fig.add_trace(go.Scatter(
            x=[_m.log10(peak_r)],
            y=[peak_s * 100],
            mode="markers",
            name=f"Peak: {peak_s * 100:.1f}%",
            marker=dict(size=12, color=_COLOR_GREEN, symbol="star"),
            hovertemplate=f"Peak: {peak_s * 100:.1f}%<extra></extra>",
        ))

    lambda_um = signal_data.get("lambda_m", 0) * 1e6
    fig.update_layout(
        title=f"Yukawa Signal vs Distance (λ = {lambda_um:.0f} µm)",
        xaxis_title="log₁₀(r / m)",
        yaxis_title="Fractional deviation (%)",
        height=400,
        showlegend=True,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Discovery Reach (bar chart by experiment)
# ---------------------------------------------------------------------------
def discovery_reach_chart(reach_data):
    """Bar chart showing discovery reach for each experiment.

    Parameters
    ----------
    reach_data : dict
        From compute_discovery_reach() with key 'experiments', each having
        name, signal_at_optimal, detectable_at_optimal.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    experiments = reach_data.get("experiments", [])

    if not experiments:
        fig = go.Figure()
        fig.add_annotation(text="No discovery reach data available", showarrow=False)
        return _apply_theme(fig)

    names = [e["name"] for e in experiments]
    signals = [e["signal_at_optimal"] * 100 for e in experiments]
    detectable = [e["detectable_at_optimal"] for e in experiments]
    colors = [_COLOR_GREEN if d else _COLOR_RED for d in detectable]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=names,
        y=signals,
        marker_color=colors,
        text=[f"{s:.1f}%" for s in signals],
        textposition="auto",
        hovertemplate="%{x}<br>signal = %{y:.1f}%<extra></extra>",
    ))

    fig.update_layout(
        title="Signal at Optimal Lambda for Each Experiment",
        xaxis_title="Experiment",
        yaxis_title="Signal at optimal (%)",
        height=400,
        showlegend=False,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Radius Landscape (N vs a_0 classification)
# ---------------------------------------------------------------------------
def radius_landscape_chart(landscape_data):
    """Scatter chart: flux quantum N vs log10(a_0 / l_Pl) colored by classification.

    Parameters
    ----------
    landscape_data : dict
        From compute_radius_landscape() with key 'landscape', each entry
        having N, a_stabilized, classification.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import math as _m

    entries = landscape_data.get("landscape", [])
    if not entries:
        fig = go.Figure()
        fig.add_annotation(text="No landscape data available", showarrow=False)
        return _apply_theme(fig)

    l_pl = 1.61625e-35  # Planck length in metres

    color_map = {
        "invisible_planck": "#6b7280",   # gray
        "sub_mm_testable": _COLOR_GREEN,  # green
        "marginal": _COLOR_HIGHLIGHT_21,  # amber
        "excluded": _COLOR_RED,           # red
    }

    ns = []
    log_a0s = []
    colors = []
    labels = []

    for e in entries:
        N = e.get("N", 0)
        a_stab = e.get("a_stabilized")
        cls = e.get("classification", "unknown")

        if a_stab is not None and a_stab > 0:
            # a_stabilized is in Planck units (l_Pl = 1), so log10(a_0/l_Pl) = log10(a_stab)
            log_a0 = _m.log10(a_stab) if a_stab > 0 else 0
        else:
            log_a0 = 0

        ns.append(N)
        log_a0s.append(log_a0)
        colors.append(color_map.get(cls, "#6b7280"))
        labels.append(cls.replace("_", " ").title())

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=ns,
        y=log_a0s,
        mode="markers+text",
        text=labels,
        textposition="top center",
        textfont=dict(size=9),
        marker=dict(size=12, color=colors),
        hovertemplate="N = %{x}<br>log₁₀(a₀/l_Pl) = %{y:.2f}<extra></extra>",
    ))

    fig.update_layout(
        title="Radius Landscape: Flux Quantum vs Stabilized Radius",
        xaxis_title="Flux Quantum N",
        yaxis_title="log₁₀(a₀ / l_Pl)",
        height=400,
        showlegend=False,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Anomaly Group Chart (gauge groups and their status)
# ---------------------------------------------------------------------------
def anomaly_group_chart(groups_data):
    """Horizontal bar chart of gauge groups colored by anomaly-free + SM status.

    Parameters
    ----------
    groups_data : dict
        From scan_anomaly_free_groups() with key 'groups', each having
        name, gs_factorizes, contains_sm, dim.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    groups = groups_data.get("groups", [])
    if not groups:
        fig = go.Figure()
        fig.add_annotation(text="No group data available", showarrow=False)
        return _apply_theme(fig)

    names = []
    dims = []
    colors = []

    for g in groups:
        names.append(g["name"])
        dims.append(g.get("dim", 0))

        if g["gs_factorizes"] and g["contains_sm"]:
            colors.append(_COLOR_GREEN)       # anomaly-free + contains SM
        elif g["gs_factorizes"]:
            colors.append(_COLOR_HIGHLIGHT_21)  # anomaly-free but no SM
        else:
            colors.append(_COLOR_RED)           # not anomaly-free

    fig = go.Figure()
    fig.add_trace(go.Bar(
        y=names,
        x=dims,
        orientation="h",
        marker_color=colors,
        text=[str(d) for d in dims],
        textposition="auto",
        hovertemplate="%{y}<br>dim = %{x}<extra></extra>",
    ))

    fig.update_layout(
        title="6D Gauge Groups: Anomaly Status",
        xaxis_title="Group Dimension",
        yaxis_title="",
        height=max(300, 50 * len(groups)),
        showlegend=False,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Cosmological Constant Scan (V_min vs N)
# ---------------------------------------------------------------------------
def cc_scan_chart(scan_data):
    """Line chart: flux quantum N vs log10(|V_min|), with observed Lambda line.

    Parameters
    ----------
    scan_data : dict
        From compute_cc_scan() with key 'scan_results', each having
        N, V_min_planck, log10_ratio.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import math as _m

    results = scan_data.get("scan_results", [])
    if not results:
        fig = go.Figure()
        fig.add_annotation(text="No CC scan data available", showarrow=False)
        return _apply_theme(fig)

    ns = []
    log_vmins = []

    for entry in results:
        V_min = entry.get("V_min_planck")
        if V_min is not None and V_min != 0:
            ns.append(entry["N"])
            log_vmins.append(_m.log10(abs(V_min)))

    if not ns:
        fig = go.Figure()
        fig.add_annotation(text="No valid V_min values", showarrow=False)
        return _apply_theme(fig)

    fig = go.Figure()

    # V_min line
    fig.add_trace(go.Scatter(
        x=ns,
        y=log_vmins,
        mode="lines+markers",
        name="log₁₀(|V_min|)",
        line=dict(color=_COLOR_DEFAULT, width=2),
        marker=dict(size=8),
        hovertemplate="N = %{x}<br>log₁₀(|V_min|) = %{y:.2f}<extra></extra>",
    ))

    # Observed Lambda line at -122
    fig.add_hline(
        y=-122, line_dash="dash", line_color=_COLOR_RED, line_width=2,
        annotation_text="Observed Λ ≈ 10⁻¹²²",
        annotation_position="bottom right",
        annotation_font_color=_COLOR_RED,
    )

    fig.update_layout(
        title="Vacuum Energy vs Flux Quantum",
        xaxis_title="Flux Quantum N",
        yaxis_title="log₁₀(|V_min| / M_Pl⁴)",
        height=400,
        showlegend=True,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Expression Density Histogram (Bridge Significance)
# ---------------------------------------------------------------------------
def expression_density_chart(histogram_data, target_value, coverage_data):
    """Bar histogram of expression values in [0.5, 2.0].

    Parameters
    ----------
    histogram_data : dict
        Keys: bin_centers (list of float), counts (list of int).
    target_value : float
        The target coefficient (~1.309).
    coverage_data : dict
        Keys: coverage_fraction (float), values_in_interval (int).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    bin_centers = histogram_data.get("bin_centers", [])
    counts = histogram_data.get("counts", [])

    if not bin_centers or not counts:
        fig = go.Figure()
        fig.add_annotation(text="No histogram data available", showarrow=False,
                           xref="paper", yref="paper", x=0.5, y=0.5)
        fig.update_layout(height=450)
        return _apply_theme(fig)

    fig = go.Figure()

    # Bar chart
    fig.add_trace(go.Bar(
        x=bin_centers,
        y=counts,
        marker_color=_COLOR_DEFAULT,
        opacity=0.6,
        name="Expressions",
        hovertemplate="Value: %{x:.3f}<br>Count: %{y}<extra></extra>",
    ))

    # Target vertical line
    fig.add_vline(
        x=target_value, line_color=_COLOR_RED, line_width=2.5, line_dash="solid",
    )
    fig.add_annotation(
        x=target_value, y=max(counts) * 0.95,
        text=f"Target: {target_value:.4f}",
        showarrow=False,
        xanchor="left", xshift=8,
        font=dict(color=_COLOR_RED, size=12),
    )

    # Epsilon band
    epsilon = 0.00016
    fig.add_vrect(
        x0=target_value * (1 - epsilon),
        x1=target_value * (1 + epsilon),
        fillcolor="rgba(248, 113, 113, 0.25)",
        line_width=0,
    )

    # Coverage annotation
    cov_frac = coverage_data.get("coverage_fraction", 0.0)
    vals_in = coverage_data.get("values_in_interval", 0)
    fig.add_annotation(
        text=f"Coverage: {cov_frac * 100:.1f}% | {vals_in} values in interval",
        showarrow=False,
        xref="paper", yref="paper",
        x=0.02, y=0.98,
        xanchor="left", yanchor="top",
        font=dict(color="#e0e0e0", size=12),
    )

    fig.update_layout(
        title="Expression Value Distribution in [0.5, 2.0]",
        xaxis_title="Expression Value",
        yaxis_title="Count",
        height=450,
        showlegend=False,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Significance Summary (Bridge Significance)
# ---------------------------------------------------------------------------
def significance_summary_chart(summary_data):
    """Side-by-side coverage fraction and p-value comparison.

    Parameters
    ----------
    summary_data : dict
        Keys: coverage (dict with coverage_fraction),
              bonferroni (dict with p_corrected),
              empirical_pvalue (dict with p_empirical),
              monte_carlo (dict with p_monte_carlo).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    from plotly.subplots import make_subplots

    coverage_frac = summary_data.get("coverage", {}).get("coverage_fraction", 0.0)
    p_bonf = summary_data.get("bonferroni", {}).get("p_corrected", 0.0)
    p_emp = summary_data.get("empirical_pvalue", {}).get("p_empirical", 0.0)
    p_mc = summary_data.get("monte_carlo", {}).get("p_monte_carlo", 0.0)

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("Interval Coverage Fraction", "P-values (lower = more significant)"),
    )

    # Left: coverage bar
    fig.add_trace(
        go.Bar(
            x=["Coverage"], y=[coverage_frac],
            marker_color=_COLOR_HIGHLIGHT_21,
            showlegend=False,
            hovertemplate="Coverage: %{y:.4f}<extra></extra>",
        ),
        row=1, col=1,
    )
    fig.add_hline(
        y=0.05, line_dash="dash", line_color=_COLOR_RED, line_width=1.5,
        annotation_text="5% significance",
        annotation_position="top right",
        annotation_font_color=_COLOR_RED,
        row=1, col=1,
    )

    # Right: p-value horizontal bars
    p_names = ["Bonferroni", "Empirical", "Monte Carlo"]
    p_values = [p_bonf, p_emp, p_mc]
    p_colors = [_COLOR_RED, _COLOR_HIGHLIGHT_21, _COLOR_GREEN]

    fig.add_trace(
        go.Bar(
            y=p_names, x=p_values,
            orientation="h",
            marker_color=p_colors,
            showlegend=False,
            hovertemplate="%{y}: %{x:.4f}<extra></extra>",
        ),
        row=1, col=2,
    )
    fig.add_vline(
        x=0.05, line_dash="dash", line_color="#ffffff", line_width=1.5,
        annotation_text="p = 0.05",
        annotation_position="top right",
        annotation_font_color="#ffffff",
        row=1, col=2,
    )

    fig.update_layout(height=400)

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Dimension Scan (scatter of (d,n) pairs colored by residual ppm)
# ---------------------------------------------------------------------------
def dimension_scan_chart(scan_data):
    """Scatter plot of (d,n) dimension pairs colored by residual ppm.

    Parameters
    ----------
    scan_data : dict
        From scan_dimension_pairs() with key 'scan_results', each having
        d, n, D, exponent, ppm. Also 'best_pair'.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import math as _m

    results = scan_data.get("scan_results", [])
    if not results:
        fig = go.Figure()
        fig.add_annotation(text="No dimension scan data", showarrow=False)
        return _apply_theme(fig)

    best = scan_data.get("best_pair", (4, 2))

    ds = [e["d"] for e in results]
    ns = [e["n"] for e in results]
    ppms = [e["ppm"] for e in results]
    exponents = [e["exponent"] for e in results]

    # Color by log10(ppm), marker size inversely related to log(ppm)
    log_ppms = [_m.log10(max(p, 1)) for p in ppms]
    max_log = max(log_ppms) if log_ppms else 1
    sizes = [max(8, 40 - 5 * lp) for lp in log_ppms]

    # Separate best point
    colors = []
    for e in results:
        if (e["d"], e["n"]) == tuple(best):
            colors.append(_COLOR_GREEN)
        else:
            colors.append(_COLOR_DEFAULT)

    hover_texts = [
        f"d={e['d']}, n={e['n']}<br>D={e['D']}, exp={e['exponent']}<br>ppm={e['ppm']:.0f}"
        for e in results
    ]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=ds,
        y=ns,
        mode="markers+text",
        marker=dict(size=sizes, color=colors, line=dict(width=1, color="#ffffff")),
        text=[f"{p:.0f}" for p in ppms],
        textposition="top center",
        textfont=dict(size=8),
        hovertext=hover_texts,
        hoverinfo="text",
    ))

    # Annotate best point
    best_entry = next((e for e in results if (e["d"], e["n"]) == tuple(best)), None)
    if best_entry:
        fig.add_annotation(
            x=best_entry["d"], y=best_entry["n"],
            text=f"BEST: ({best[0]},{best[1]})<br>{best_entry['ppm']:.0f} ppm",
            showarrow=True, arrowhead=2, arrowcolor=_COLOR_GREEN,
            font=dict(color=_COLOR_GREEN, size=12),
            bgcolor="rgba(26, 29, 35, 0.9)",
            bordercolor=_COLOR_GREEN,
        )

    fig.update_layout(
        title="Dimension Scan: Residual ppm by (d, n) Pair",
        xaxis_title="External dimensions d",
        yaxis_title="Internal dimensions n",
        height=500,
        showlegend=False,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Residual Analysis (horizontal bar of correction candidates)
# ---------------------------------------------------------------------------
def residual_analysis_chart(residual_data):
    """Horizontal bar chart of correction candidates for the 688 ppm residual.

    Parameters
    ----------
    residual_data : dict
        From analyze_residual() with keys 'corrections' (list of dicts
        with name, corrected_ppm, improvement) and 'raw_ppm'.

    Returns
    -------
    plotly.graph_objects.Figure
    """
    corrections = residual_data.get("corrections", [])
    raw_ppm = residual_data.get("raw_ppm", 0)

    if not corrections:
        fig = go.Figure()
        fig.add_annotation(text="No correction data available", showarrow=False)
        return _apply_theme(fig)

    # Sort by absolute corrected_ppm
    sorted_corr = sorted(corrections, key=lambda c: abs(c["corrected_ppm"]))

    names = [c["name"] for c in sorted_corr]
    abs_ppms = [abs(c["corrected_ppm"]) for c in sorted_corr]
    colors = [_COLOR_GREEN if c.get("improvement", False) else _COLOR_RED for c in sorted_corr]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        y=names,
        x=abs_ppms,
        orientation="h",
        marker_color=colors,
        text=[f"{p:.0f}" for p in abs_ppms],
        textposition="auto",
        hovertemplate="%{y}<br>|corrected ppm| = %{x:.0f}<extra></extra>",
    ))

    # Reference line at raw residual
    fig.add_vline(
        x=abs(raw_ppm), line_dash="dash", line_color="#ffffff", line_width=1.5,
        annotation_text=f"Raw: {abs(raw_ppm):.0f} ppm",
        annotation_position="top right",
        annotation_font_color="#ffffff",
    )

    fig.update_layout(
        title="Correction Candidates: |Residual| After Each Correction",
        xaxis_title="|Corrected ppm|",
        yaxis_title="",
        height=max(350, 30 * len(corrections)),
        showlegend=False,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Swampland Distance Conjecture Lambda Rates
# ---------------------------------------------------------------------------
def sdc_lambda_chart(sdc_data):
    """Horizontal bar chart comparing lambda values against the SDC bound.

    Parameters
    ----------
    sdc_data : dict
        Keys: lambda_EH (float), lambda_GB (float), lambda_bound (float),
        sdc_satisfied_EH (bool), sdc_satisfied_GB (bool).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    lambda_EH = sdc_data.get("lambda_EH", 0.0)
    lambda_GB = sdc_data.get("lambda_GB", 0.0)
    lambda_bound = sdc_data.get("lambda_bound", 0.0)
    satisfied_EH = sdc_data.get("sdc_satisfied_EH", False)
    satisfied_GB = sdc_data.get("sdc_satisfied_GB", False)

    names = ["Pure EH", "GB-corrected"]
    values = [lambda_EH, lambda_GB]
    colors = [
        _COLOR_GREEN if satisfied_EH else _COLOR_RED,
        _COLOR_GREEN if satisfied_GB else _COLOR_RED,
    ]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        y=names,
        x=values,
        orientation="h",
        marker_color=colors,
        text=[f"{v:.4f}" for v in values],
        textposition="auto",
        hovertemplate="%{y}<br>λ = %{x:.4f}<extra></extra>",
    ))

    # SDC bound vertical line
    fig.add_vline(
        x=lambda_bound, line_dash="dash", line_color="#ffffff", line_width=2,
        annotation_text=f"SDC bound ({lambda_bound:.3f})",
        annotation_position="top right",
        annotation_font_color="#ffffff",
    )

    fig.update_layout(
        title="Swampland Distance Conjecture: λ Rates",
        xaxis_title="λ (SDC rate)",
        yaxis_title="",
        height=300,
        showlegend=False,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Emergence M_Pl^2 at Canonical Radii
# ---------------------------------------------------------------------------
def emergence_R_scan_chart(emergence_data):
    """Bar chart of M_Pl^2 ratio for different canonical radii.

    Parameters
    ----------
    emergence_data : dict
        Key 'canonical_R_results': list of dicts each with
        name (str), R_m (float), M_Pl_sq_ratio (float).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    import math as _m

    results = emergence_data.get("canonical_R_results", [])
    if not results:
        fig = go.Figure()
        fig.add_annotation(text="No emergence data available", showarrow=False)
        return _apply_theme(fig)

    names = [r["name"] for r in results]
    ratios = [r["M_Pl_sq_ratio"] for r in results]
    colors = [
        _COLOR_GREEN if 0.1 <= r["M_Pl_sq_ratio"] <= 10.0 else _COLOR_DEFAULT
        for r in results
    ]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=names,
        y=ratios,
        marker_color=colors,
        text=[f"{r:.2e}" for r in ratios],
        textposition="auto",
        hovertemplate="%{x}<br>M_Pl² ratio = %{y:.3e}<extra></extra>",
    ))

    # Match line at y=1
    fig.add_hline(
        y=1.0, line_dash="dash", line_color=_COLOR_RED, line_width=2,
        annotation_text="Match",
        annotation_position="top right",
        annotation_font_color=_COLOR_RED,
    )

    fig.update_layout(
        title="Emergence: M_Pl² from KK Tower at Canonical Radii",
        xaxis_title="",
        yaxis_title="M_Pl² (predicted) / M_Pl² (measured)",
        yaxis_type="log",
        height=450,
        showlegend=False,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Correction Series Convergence (dual-axis: coefficients + residual)
# ---------------------------------------------------------------------------
def correction_series_chart(series_data):
    """Dual-axis chart: bar coefficients (left) and residual ppm line (right).

    Parameters
    ----------
    series_data : dict
        Key 'coefficients': list of dicts each with
        order (int), coefficient (float), residual_after_ppm (float).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    from plotly.subplots import make_subplots

    coefficients = series_data.get("coefficients", [])
    if not coefficients:
        fig = go.Figure()
        fig.add_annotation(text="No correction series data available",
                           showarrow=False)
        return _apply_theme(fig)

    orders = [c["order"] for c in coefficients]
    coeff_abs = [abs(c["coefficient"]) for c in coefficients]
    residuals = [c["residual_after_ppm"] for c in coefficients]
    x_labels = [f"α^{o}" for o in orders]

    # Decide whether right axis should be log scale
    r_pos = [r for r in residuals if r > 0]
    use_log = (max(r_pos) / min(r_pos) > 50) if len(r_pos) >= 2 else False

    # Green palette for bars, scaled by magnitude
    max_coeff = max(coeff_abs) if coeff_abs else 1.0
    green_colors = []
    for ca in coeff_abs:
        frac = ca / max_coeff if max_coeff > 0 else 0.5
        r_val = int(20 + 30 * (1 - frac))
        g_val = int(140 + 115 * frac)
        b_val = int(80 + 70 * (1 - frac))
        green_colors.append(f"rgb({r_val}, {g_val}, {b_val})")

    fig = make_subplots(specs=[[{"secondary_y": True}]])

    # Left axis: coefficient bars
    fig.add_trace(
        go.Bar(
            x=x_labels,
            y=coeff_abs,
            name="|Coefficient|",
            marker_color=green_colors,
            text=[f"{ca:.2f}" for ca in coeff_abs],
            textposition="auto",
            hovertemplate="%{x}<br>|c| = %{y:.4f}<extra></extra>",
        ),
        secondary_y=False,
    )

    # Right axis: residual line
    fig.add_trace(
        go.Scatter(
            x=x_labels,
            y=residuals,
            name="Residual (ppm)",
            mode="lines+markers",
            line=dict(color="#fb923c", width=2.5),   # orange-400
            marker=dict(size=10, color="#f87171"),     # red-400
            hovertemplate="%{x}<br>residual = %{y:.2f} ppm<extra></extra>",
        ),
        secondary_y=True,
    )

    # Reference lines on right axis
    fig.add_hline(
        y=0.62, line_dash="dash", line_color="#94a3b8", line_width=1.5,
        annotation_text="Leading-order residual (0.62 ppm)",
        annotation_position="top left",
        annotation_font_color="#94a3b8",
        secondary_y=True,
    )
    fig.add_hline(
        y=22, line_dash="dash", line_color="#facc15", line_width=1.5,
        annotation_text="CODATA uncertainty",
        annotation_position="bottom left",
        annotation_font_color="#facc15",
        secondary_y=True,
    )

    # Axis labels
    fig.update_yaxes(title_text="|Coefficient|", secondary_y=False)
    fig.update_yaxes(
        title_text="Residual after correction (ppm)",
        type="log" if use_log else "linear",
        secondary_y=True,
    )

    fig.update_layout(
        title="Correction Series: Coefficients and Residual Convergence",
        xaxis_title="Correction Order (αᵏ)",
        height=480,
        showlegend=True,
        legend=dict(x=0.01, y=0.99, xanchor="left", yanchor="top",
                    bgcolor="rgba(14, 17, 23, 0.7)"),
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Mu Offset Scan (horizontal bar chart of residual_ppm vs k)
# ---------------------------------------------------------------------------
def mu_offset_scan_chart(offset_data):
    """Horizontal bar chart showing residual_ppm for each k value tested.

    Parameters
    ----------
    offset_data : dict
        Dict with key ``offsets``, a list of dicts each containing
        ``k_label`` (str), ``k_value`` (float), and ``residual_ppm`` (float).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    offsets = offset_data.get("offsets", [])

    k_labels = [o["k_label"] for o in offsets]
    residuals = [o["residual_ppm"] for o in offsets]

    # Color logic: green for sqrt(phi), else red/blue by sign
    colors = []
    for o in offsets:
        if "sqrt(phi)" in o["k_label"]:
            colors.append("#34d399")
        elif o["residual_ppm"] >= 0:
            colors.append("#f87171")
        else:
            colors.append("#60a5fa")

    fig = go.Figure()

    fig.add_trace(go.Bar(
        y=k_labels,
        x=residuals,
        orientation="h",
        marker_color=colors,
        hovertemplate="k = %{y}<br>Residual = %{x:.1f} ppm<extra></extra>",
    ))

    # Zero line
    fig.add_vline(x=0, line_dash="dash", line_color="#fafafa", line_width=1)

    # CODATA uncertainty band at +/-22 ppm
    fig.add_vline(x=22, line_dash="dash", line_color="gray", line_width=1)
    fig.add_vline(x=-22, line_dash="dash", line_color="gray", line_width=1)

    fig.update_layout(
        title="Residual (ppm) for α²⁴ · μ · (μ − k)",
        xaxis_title="Residual (ppm)",
        height=400,
        showlegend=False,
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Formula Comparison (horizontal bar chart of bridge/hierarchy residuals)
# ---------------------------------------------------------------------------
def formula_comparison_chart(comparison_data):
    """Horizontal bar chart comparing all bridge/hierarchy formulas.

    Parameters
    ----------
    comparison_data : dict
        Dict with key ``formulas``, a list of dicts each containing
        ``formula_label`` (str), ``residual_ppm`` (float),
        ``abs_residual_ppm`` (float), and ``n_fitted_params`` (int or str).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    formulas = comparison_data.get("formulas", [])

    labels = [f["formula_label"] for f in formulas]
    residuals = [f["residual_ppm"] for f in formulas]

    # Color by abs residual magnitude
    colors = []
    for f in formulas:
        abs_r = f["abs_residual_ppm"]
        if abs_r < 22:
            colors.append("#34d399")
        elif abs_r < 200:
            colors.append("#f59e0b")
        else:
            colors.append("#f87171")

    fig = go.Figure()

    fig.add_trace(go.Bar(
        y=labels,
        x=residuals,
        orientation="h",
        marker_color=colors,
        hovertemplate="%{y}<br>Residual = %{x:.1f} ppm<extra></extra>",
    ))

    # Zero line
    fig.add_vline(x=0, line_dash="dash", line_color="#fafafa", line_width=1)

    # CODATA uncertainty band at +/-22 ppm
    fig.add_vline(x=22, line_dash="dash", line_color="gray", line_width=1)
    fig.add_vline(x=-22, line_dash="dash", line_color="gray", line_width=1)

    fig.update_layout(
        title="Bridge and Hierarchy Formula Residuals (ppm)",
        xaxis_title="Residual (ppm)",
        height=350,
        showlegend=False,
    )

    return _apply_theme(fig)


def residual_scan_chart(scan_data):
    """Horizontal bar chart of candidate expressions and their closeness to delta.

    Parameters
    ----------
    scan_data : dict
        Must have "best_matches" key with list of dicts having "name" and "residual_ppm".

    Returns
    -------
    plotly.graph_objects.Figure
    """
    entries = scan_data.get("best_matches", [])[:10]
    if not entries:
        fig = go.Figure()
        fig.add_annotation(text="No data", xref="paper", yref="paper", x=0.5, y=0.5)
        return _apply_theme(fig)

    names = [e["name"] for e in entries]
    ppms = [e["residual_ppm"] for e in entries]

    colors = ["#34d399" if p < 100 else "#f59e0b" if p < 1000 else "#f87171" for p in ppms]

    fig = go.Figure(go.Bar(
        x=ppms,
        y=names,
        orientation="h",
        marker_color=colors,
        text=[f"{p:.0f}" for p in ppms],
        textposition="outside",
    ))

    fig.update_layout(
        title="Closest Composite Expressions to δ",
        xaxis_title="Residual (ppm)",
        yaxis_title="",
        height=400,
        margin=dict(l=250, r=50, t=50, b=50),
        xaxis=dict(type="log"),
    )

    return _apply_theme(fig)


def k_closed_form_chart(k_data):
    """Horizontal bar chart of closed-form candidates for k_exact.

    Parameters
    ----------
    k_data : dict
        Must have "best_matches" key with list of dicts having "name" and "residual_ppm".

    Returns
    -------
    plotly.graph_objects.Figure
    """
    entries = k_data.get("best_matches", [])[:10]
    if not entries:
        fig = go.Figure()
        fig.add_annotation(text="No data", xref="paper", yref="paper", x=0.5, y=0.5)
        return _apply_theme(fig)

    names = [e["name"] for e in entries]
    ppms = [e["residual_ppm"] for e in entries]

    colors = ["#34d399" if p < 100 else "#f59e0b" if p < 1000 else "#60a5fa" for p in ppms]

    fig = go.Figure(go.Bar(
        x=ppms,
        y=names,
        orientation="h",
        marker_color=colors,
        text=[f"{p:.0f}" for p in ppms],
        textposition="outside",
    ))

    fig.update_layout(
        title="Closed-Form Candidates for k_exact",
        xaxis_title="Residual vs k_exact (ppm)",
        yaxis_title="",
        height=400,
        margin=dict(l=280, r=50, t=50, b=50),
        xaxis=dict(type="log"),
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Coleman-Weinberg + Flux Potential vs Radius
# ---------------------------------------------------------------------------
def cw_potential_chart(scan_data):
    """Line chart of V_CW, V_flux_min, and V_total vs a_0.

    Parameters
    ----------
    scan_data : dict
        Keys: a_0_values (list), V_cw (list), V_flux_min (list), V_total (list).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not scan_data or not scan_data.get("a_0_values"):
        fig = go.Figure()
        fig.add_annotation(text="No data available", showarrow=False,
                           xref="paper", yref="paper", x=0.5, y=0.5)
        return _apply_theme(fig)

    a0 = scan_data["a_0_values"]
    V_cw = scan_data["V_cw"]
    V_flux = scan_data["V_flux_min"]
    V_total = scan_data["V_total"]

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=a0, y=V_cw,
        mode="lines",
        line=dict(color=_COLOR_DEFAULT, width=2, dash="dash"),
        name="V_CW",
        hovertemplate="a₀ = %{x:.4e}<br>V_CW = %{y:.4e}<extra></extra>",
    ))

    fig.add_trace(go.Scatter(
        x=a0, y=V_flux,
        mode="lines",
        line=dict(color=_COLOR_HIGHLIGHT_21, width=2, dash="dot"),
        name="V_flux_min",
        hovertemplate="a₀ = %{x:.4e}<br>V_flux = %{y:.4e}<extra></extra>",
    ))

    fig.add_trace(go.Scatter(
        x=a0, y=V_total,
        mode="lines",
        line=dict(color=_COLOR_GREEN, width=3),
        name="V_total",
        hovertemplate="a₀ = %{x:.4e}<br>V_total = %{y:.4e}<extra></extra>",
    ))

    # Mark minimum of V_total if it exists (not at boundary)
    if len(V_total) > 2:
        min_idx = min(range(len(V_total)), key=lambda i: V_total[i])
        if 0 < min_idx < len(V_total) - 1:
            fig.add_trace(go.Scatter(
                x=[a0[min_idx]],
                y=[V_total[min_idx]],
                mode="markers",
                marker=dict(
                    color=_COLOR_RED, size=14,
                    symbol="star", line=dict(color="#ffffff", width=2),
                ),
                name="Minimum",
                hovertemplate=(
                    f"a₀ = {a0[min_idx]:.4e}<br>"
                    f"V_min = {V_total[min_idx]:.4e}<extra></extra>"
                ),
            ))

    fig.add_hline(y=0, line_dash="solid", line_color="#4c566a", line_width=1)

    fig.update_layout(
        title="Coleman-Weinberg + Flux Potential vs Radius",
        xaxis_title="a₀ (Planck units)",
        yaxis_title="Effective Potential (Planck units)",
        xaxis=dict(type="log"),
        height=400,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Warped S^2 Potential vs Radius
# ---------------------------------------------------------------------------
def warped_potential_chart(warp_data):
    """Line chart of V_total(a_0) for different warp parameters epsilon.

    Parameters
    ----------
    warp_data : dict
        Keys: a_0_values (list), epsilon_values (list),
        V_total_by_epsilon (dict mapping epsilon -> list of V values).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not warp_data or not warp_data.get("a_0_values"):
        fig = go.Figure()
        fig.add_annotation(text="No data available", showarrow=False,
                           xref="paper", yref="paper", x=0.5, y=0.5)
        return _apply_theme(fig)

    a0 = warp_data["a_0_values"]
    epsilon_values = warp_data.get("epsilon_values", [])
    V_by_eps = warp_data.get("V_total_by_epsilon", {})

    _PALETTE = [
        _COLOR_DEFAULT, _COLOR_HIGHLIGHT_21, _COLOR_GREEN,
        _COLOR_HIGHLIGHT_10, _COLOR_RED,
        "#38bdf8", "#fb923c", "#a3e635", "#e879f9", "#fbbf24",
    ]

    fig = go.Figure()

    for i, eps in enumerate(epsilon_values):
        key = eps if eps in V_by_eps else str(eps)
        V_vals = V_by_eps.get(key, V_by_eps.get(eps, []))
        if not V_vals:
            continue
        color = _PALETTE[i % len(_PALETTE)]
        fig.add_trace(go.Scatter(
            x=a0, y=V_vals,
            mode="lines",
            line=dict(color=color, width=2),
            name=f"eps = {eps}",
            hovertemplate=f"eps = {eps}<br>" + "a₀ = %{x:.4e}<br>V = %{y:.4e}<extra></extra>",
        ))

    fig.add_hline(y=0, line_dash="solid", line_color="#4c566a", line_width=1)

    fig.update_layout(
        title="Warped S² Potential vs Radius",
        xaxis_title="a₀ (Planck units)",
        yaxis_title="V_total (Planck units)",
        xaxis=dict(type="log"),
        height=400,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return _apply_theme(fig)


# ---------------------------------------------------------------------------
# Chart: Orbifold S^2/Z_2 Potential vs Radius
# ---------------------------------------------------------------------------
def orbifold_potential_chart(orbifold_data):
    """Line chart of V_casimir, V_brane, and V_total for the orbifold mechanism.

    Parameters
    ----------
    orbifold_data : dict
        Keys: a_0_values (list), V_casimir (list), V_brane (list),
        V_total (list), T_0 (float).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    if not orbifold_data or not orbifold_data.get("a_0_values"):
        fig = go.Figure()
        fig.add_annotation(text="No data available", showarrow=False,
                           xref="paper", yref="paper", x=0.5, y=0.5)
        return _apply_theme(fig)

    a0 = orbifold_data["a_0_values"]
    V_casimir = orbifold_data["V_casimir"]
    V_brane = orbifold_data["V_brane"]
    V_total = orbifold_data["V_total"]
    T_0 = orbifold_data.get("T_0", None)

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=a0, y=V_casimir,
        mode="lines",
        line=dict(color=_COLOR_DEFAULT, width=2, dash="dash"),
        name="V_casimir",
        hovertemplate="a₀ = %{x:.4e}<br>V_casimir = %{y:.4e}<extra></extra>",
    ))

    fig.add_trace(go.Scatter(
        x=a0, y=V_brane,
        mode="lines",
        line=dict(color=_COLOR_HIGHLIGHT_21, width=2, dash="dot"),
        name="V_brane",
        hovertemplate="a₀ = %{x:.4e}<br>V_brane = %{y:.4e}<extra></extra>",
    ))

    fig.add_trace(go.Scatter(
        x=a0, y=V_total,
        mode="lines",
        line=dict(color=_COLOR_GREEN, width=3),
        name="V_total",
        hovertemplate="a₀ = %{x:.4e}<br>V_total = %{y:.4e}<extra></extra>",
    ))

    # Mark minimum of V_total if it exists (not at boundary)
    if len(V_total) > 2:
        min_idx = min(range(len(V_total)), key=lambda i: V_total[i])
        if 0 < min_idx < len(V_total) - 1:
            fig.add_trace(go.Scatter(
                x=[a0[min_idx]],
                y=[V_total[min_idx]],
                mode="markers",
                marker=dict(
                    color=_COLOR_RED, size=14,
                    symbol="star", line=dict(color="#ffffff", width=2),
                ),
                name="Minimum",
                hovertemplate=(
                    f"a₀ = {a0[min_idx]:.4e}<br>"
                    f"V_min = {V_total[min_idx]:.4e}<extra></extra>"
                ),
            ))

    fig.add_hline(y=0, line_dash="solid", line_color="#4c566a", line_width=1)

    title = "Orbifold S²/ℤ₂ Potential vs Radius"
    if T_0 is not None:
        title += f" (T_0 = {T_0:.4e})"

    fig.update_layout(
        title=title,
        xaxis_title="a₀ (Planck units)",
        yaxis_title="Effective Potential (Planck units)",
        xaxis=dict(type="log"),
        height=400,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return _apply_theme(fig)