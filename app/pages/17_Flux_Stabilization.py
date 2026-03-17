"""
Flux Stabilization -- Page 17

Adding quantized 2-form flux on S^2 resolves the Casimir no-go result
from Page 16, creating a stable minimum for the volume modulus.
"""

import streamlit as st
import math
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
# Try to load core modules
# ---------------------------------------------------------------------------
_PHI = (1 + math.sqrt(5)) / 2
_SQRT5 = math.sqrt(5)

_core_available = False
try:
    from alpha_ladder_core.flux_stabilization import (  # noqa: E402
        compute_flux_potential,
        find_flux_minimum,
        scan_flux_quanta,
    )
    _core_available = True
except ImportError:
    pass

# ---------------------------------------------------------------------------
# Compute values (guarded)
# ---------------------------------------------------------------------------
_potential = None
_minimum = None
_scan = None

if _core_available:
    try:
        _potential = compute_flux_potential()
    except Exception:
        pass
    try:
        _minimum = find_flux_minimum()
    except Exception:
        pass
    try:
        _scan = scan_flux_quanta()
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Flux Stabilization: Resolving the Casimir No-Go")
st.markdown(
    "Adding quantized 2-form flux creates a stable minimum for the volume "
    "modulus. This resolves the negative result from Page 16 and closes "
    "Gap #3 (volume stabilization) while simultaneously dissolving Gap #1 "
    "(dilaton mass/screening)."
)
st.divider()

# ---------------------------------------------------------------------------
# 3 Metric Cards
# ---------------------------------------------------------------------------
col_m1, col_m2, col_m3 = st.columns(3)

with col_m1:
    st.metric(label="Gap #3 Status", value="Closed")
with col_m2:
    st.metric(label="Gap #1 Status", value="Closed")
with col_m3:
    st.metric(label="Mechanism", value="F₂ flux on S²")

st.markdown("")

# ---------------------------------------------------------------------------
# Section A: The Three-Term Potential
# ---------------------------------------------------------------------------
with st.expander("A. The Three-Term Potential", expanded=True):
    st.markdown(
        """
        <div class="formula-card">
        <b>Three-Term Effective Potential</b><br><br>
        <center><b>V(&sigma;) = A e<sup>4&sigma;</sup>
        + B e<sup>2&sigma;</sup>
        + C e<sup>6&sigma;</sup></b></center><br>
        where &sigma; is the volume modulus (breathing mode) of the internal
        S&sup2;. Each term has a distinct physical origin and scaling.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="step-card">
        <b>A e<sup>4&sigma;</sup> &mdash; Casimir (1-loop quantum)</b><br><br>
        The 1-loop Casimir energy from the KK graviton tower on S&sup2;.
        As computed on Page 16, <b>A &lt; 0</b> for the bosonic graviton
        tower. This term drives decompactification on its own.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown(
        """
        <div class="step-card">
        <b>B e<sup>2&sigma;</sup> &mdash; Curvature (tree-level)</b><br><br>
        The classical curvature contribution from the internal manifold.
        For S&sup2;, R<sub>2</sub> = 2/a&sup2; &gt; 0 but after
        dimensional reduction <b>B &lt; 0</b> in the effective potential
        (the curvature acts as a negative cosmological constant in the
        modulus direction).
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown(
        """
        <div class="step-card">
        <b>C e<sup>6&sigma;</sup> &mdash; Flux (quantized, the new ingredient)</b><br><br>
        A quantized 2-form field strength F<sub>2</sub> threading the
        internal S&sup2;. Dirac quantization requires
        &int;<sub>S&sup2;</sub> F<sub>2</sub> = 2&pi;N for integer N.
        The flux energy density scales as |F<sub>2</sub>|&sup2; &sim;
        N&sup2;/a<sup>4</sup>, giving <b>C &gt; 0</b>. This is the
        stabilizing term.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="theorem-card">
        <b>Key insight:</b> The flux term e<sup>6&sigma;</sup> grows
        faster than the Casimir term e<sup>4&sigma;</sup> as &sigma;
        increases. This creates a potential wall at large &sigma; that
        prevents the volume modulus from running away to
        decompactification. Combined with the curvature term at small
        &sigma;, a stable minimum is guaranteed for any N &ge; 1.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Section B: Stabilized Minimum
# ---------------------------------------------------------------------------
with st.expander("B. Stabilized Minimum", expanded=True):

    if _core_available and _potential:
        # Potential chart
        try:
            from app.components.charts import flux_potential_chart  # noqa: E402
            min_location = None
            if _minimum:
                min_location = _minimum.get("sigma_0")
            fig_pot = flux_potential_chart(_potential, min_location=min_location)
            st.plotly_chart(fig_pot, use_container_width=True)
        except ImportError:
            st.info(
                "Chart function flux_potential_chart not yet available "
                "in charts module."
            )
        except Exception as exc:
            st.warning(f"Flux potential chart error: {exc}")

        st.markdown("")

        if _minimum:
            col_b1, col_b2, col_b3 = st.columns(3)
            with col_b1:
                sigma_0 = _minimum.get("sigma_0")
                st.metric(
                    label="σ₀ (minimum)",
                    value=fmt_decimal(sigma_0, sig_figs=6) if sigma_0 is not None else "N/A",
                )
            with col_b2:
                V_min = _minimum.get("V_min")
                st.metric(
                    label="V min",
                    value=fmt_decimal(V_min, sig_figs=6) if V_min is not None else "N/A",
                )
            with col_b3:
                V_pp = _minimum.get("V_double_prime")
                stability = "Stable" if V_pp is not None and V_pp > 0 else "N/A"
                st.metric(label="V'' > 0?", value=stability)

        st.markdown("")

        if _scan and _scan.get("results"):
            st.markdown("#### Flux Quantum Scan (N = 1 .. 5)")

            scan_df = pd.DataFrame(_scan["results"])
            display_cols = [c for c in ["N", "sigma_0", "V_min", "m_phi_eV"] if c in scan_df.columns]
            if display_cols:
                st.dataframe(scan_df[display_cols], use_container_width=True)

            # Flux scan chart
            try:
                from app.components.charts import flux_scan_chart  # noqa: E402
                fig_scan = flux_scan_chart(_scan)
                st.plotly_chart(fig_scan, use_container_width=True)
            except ImportError:
                pass
            except Exception as exc:
                st.warning(f"Flux scan chart error: {exc}")

    else:
        st.markdown(
            """
            <div class="step-card">
            <b>Stabilization mechanism:</b><br><br>
            Setting dV/d&sigma; = 0 for the three-term potential:<br><br>
            4A e<sup>4&sigma;</sup> + 2B e<sup>2&sigma;</sup>
            + 6C e<sup>6&sigma;</sup> = 0<br><br>
            With A &lt; 0, B &lt; 0, and C &gt; 0, the equation
            admits a solution &sigma;<sub>0</sub> where
            d&sup2;V/d&sigma;&sup2; &gt; 0 (stable minimum).<br><br>
            The minimum exists for <em>every</em> integer flux quantum
            N &ge; 1, with &sigma;<sub>0</sub> shifting slightly as
            N increases.
            </div>
            """,
            unsafe_allow_html=True,
        )

        st.markdown("")

        st.markdown(
            """
            <div class="formula-card">
            <b>Stability condition:</b><br><br>
            <code>d&sup2;V/d&sigma;&sup2;|<sub>&sigma;=&sigma;<sub>0</sub></sub>
            = 16A e<sup>4&sigma;<sub>0</sub></sup>
            + 4B e<sup>2&sigma;<sub>0</sub></sup>
            + 36C e<sup>6&sigma;<sub>0</sub></sup> &gt; 0</code><br><br>
            The dominant C e<sup>6&sigma;</sup> term ensures V'' &gt; 0
            at the minimum.
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section C: The Planck Mass Consequence
# ---------------------------------------------------------------------------
with st.expander("C. The Planck Mass Consequence", expanded=True):

    st.markdown(
        """
        <div class="warning-card">
        <b>Critical result: The stabilized dilaton mass is at the Planck
        scale (~10<sup>28</sup> eV)</b><br><br>
        The flux stabilization mechanism works <em>too well</em>. Because
        all three potential ingredients are Planck-scale quantities,
        the resulting dilaton mass is enormous &mdash; far beyond any
        observable energy scale.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="proof-card">
        <b>All potential ingredients are Planck-scale:</b><br><br>
        &bull; Casimir: V<sub>Cas</sub> ~ 1/a<sup>4</sup> with
        a ~ l<sub>Pl</sub><br>
        &bull; Curvature: V<sub>curv</sub> ~ 1/a&sup2; with
        a ~ l<sub>Pl</sub><br>
        &bull; Flux: V<sub>flux</sub> ~ N&sup2;/a<sup>6</sup> with
        a ~ l<sub>Pl</sub><br><br>
        The minimum occurs at a ~ l<sub>Pl</sub>, giving
        m<sub>&phi;</sub> ~ M<sub>Pl</sub>.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown(
        """
        <div class="proof-card">
        <b>m<sub>&phi;</sub> ~ M<sub>Pl</sub> means complete decoupling:</b><br><br>
        The Yukawa profile of the dilaton is
        exp(&minus;m<sub>&phi;</sub> r). With m<sub>&phi;</sub> ~
        M<sub>Pl</sub> ~ 10<sup>28</sup> eV:<br><br>
        &bull; exp(&minus;M<sub>Pl</sub> &middot; r) = 0 for any
        r &gt; 10<sup>&minus;35</sup> m<br>
        &bull; The dilaton field vanishes at <em>all</em> observable
        distances<br>
        &bull; No fifth force, no PPN deviation, no screening
        phenomenology
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown(
        """
        <div class="proof-card">
        <b>The screening gap dissolves:</b><br><br>
        On Page 15, the tree-level screening amplitude was
        &alpha;<sub>eff</sub> = 0.618 &times; exp(&minus;m<sub>&phi;</sub> r).
        With m<sub>&phi;</sub> ~ M<sub>Pl</sub>:<br><br>
        &bull; &alpha;<sub>eff</sub> = 0.618 &times; 0 = 0 at all
        laboratory and astrophysical distances<br>
        &bull; The 3854x screening gap between tree-level and empirical
        amplitudes is dissolved &mdash; not by reducing the tree-level
        prediction, but by killing the dilaton field entirely<br>
        &bull; G<sub>eff</sub>(r) = G<sub>Newton</sub> everywhere, as
        observed
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="theorem-card">
        <b>Summary:</b> Gap #3 is closed (a stable minimum exists for
        the volume modulus). Gap #1 is dissolved (the dilaton is invisible
        at all observable scales). The framework reduces to pure GR at
        all observable distances, consistent with all experimental
        bounds.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Section D: What Flux Adds to the Framework
# ---------------------------------------------------------------------------
with st.expander("D. What Flux Adds to the Framework", expanded=True):

    st.markdown(
        """
        <div class="formula-card">
        <b>Honest assessment of the flux stabilization result:</b><br><br>
        &bull; The flux quantum N is a <em>discrete</em> parameter
        (integer), not a continuous tuning knob. Dirac quantization
        fixes N &in; {1, 2, 3, ...}.<br><br>
        &bull; The stable minimum exists for <em>all</em> N &ge; 1.
        There is no fine-tuning: any nonzero flux stabilizes the
        volume modulus. This is a robust structural result.<br><br>
        &bull; However, the result is "too good" &mdash; it kills
        <em>all</em> fifth-force phenomenology. The dilaton becomes
        completely invisible, and the framework reduces to pure GR
        with no observable signatures beyond standard gravity.<br><br>
        &bull; The framework's prediction for G remains
        (G = &phi;&sup2;/2 bridge), but the screening mechanism
        becomes irrelevant &mdash; there is nothing to screen.<br><br>
        &bull; For observable consequences to emerge, the internal
        radius a<sub>0</sub> must be <em>much larger</em> than
        l<sub>Pl</sub>. See <b>Page 18</b> for the experimental
        landscape parametrized by a<sub>0</sub>.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Flux Stabilization | Alpha Ladder Research Dashboard")
