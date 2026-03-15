"""
Radius Phenomenology -- Page 18

Observable consequences of the Alpha Ladder framework parametrized by
the internal radius a_0. The dilaton mass scales as m_phi ~ M_Pl * l_Pl / a_0,
so larger extra dimensions produce a lighter, potentially detectable dilaton.
"""

import streamlit as st
import math
import pandas as pd


# ---------------------------------------------------------------------------
# Custom CSS (matches Page 16)
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
# Try to load core modules
# ---------------------------------------------------------------------------
_PHI = (1 + math.sqrt(5)) / 2

_core_available = False
try:
    from alpha_ladder_core.radius_phenomenology import (  # noqa: E402
        compute_dilaton_mass_vs_radius,
        compute_testable_window,
    )
    _core_available = True
except ImportError:
    pass

_bounds_available = False
try:
    from alpha_ladder_core.radius_phenomenology import (  # noqa: E402
        compute_experimental_bounds,
    )
    _bounds_available = True
except ImportError:
    pass

# ---------------------------------------------------------------------------
# Compute values (guarded)
# ---------------------------------------------------------------------------
_mass_data = None
_window = None
_bounds = None

if _core_available:
    try:
        _mass_data = compute_dilaton_mass_vs_radius()
    except Exception:
        pass
    try:
        _window = compute_testable_window()
    except Exception:
        pass

if _bounds_available:
    try:
        _bounds = compute_experimental_bounds()
    except Exception:
        pass

# ---------------------------------------------------------------------------
# Physical constants for fallback calculations
# ---------------------------------------------------------------------------
_M_PL_EV = 1.22e28       # Planck mass in eV
_L_PL = 1.616e-35        # Planck length in meters
_HC_EV_M = 1.9733e-7     # hbar*c in eV*m

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("Radius Phenomenology: The Experimental Landscape")
st.markdown(
    "Observable consequences of the Alpha Ladder framework parametrized "
    "by the internal radius a\u2080. The dilaton mass scales inversely "
    "with a\u2080, creating a direct link between extra-dimension size "
    "and experimental detectability."
)
st.divider()

# ---------------------------------------------------------------------------
# 4 Metric Cards
# ---------------------------------------------------------------------------
col_m1, col_m2, col_m3, col_m4 = st.columns(4)

with col_m1:
    st.metric(label="a\u2080 range", value="l_Pl .. 1 mm")
with col_m2:
    st.metric(label="Testable window", value="10 \u03bcm .. 100 \u03bcm")
with col_m3:
    st.metric(label="Key threshold", value="2 meV")
with col_m4:
    st.metric(label="Best probe", value="Eot-Wash")

st.markdown("")

# ---------------------------------------------------------------------------
# Section A: The Scaling Law
# ---------------------------------------------------------------------------
with st.expander("A. The Scaling Law", expanded=True):

    st.markdown(
        """
        <div class="formula-card">
        <b>Dilaton Mass Scaling</b><br><br>
        <center><b>m<sub>&phi;</sub> = M<sub>Pl</sub> &middot;
        l<sub>Pl</sub> / a<sub>0</sub></b></center><br>
        The dilaton mass is set by the ratio of the Planck length to
        the internal radius. Larger extra dimensions produce a lighter
        dilaton with a longer Yukawa range &lambda; = 1/m<sub>&phi;</sub>.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="step-card">
        <b>The scaling explained:</b><br><br>
        &bull; At a<sub>0</sub> = l<sub>Pl</sub>, the dilaton has
        m<sub>&phi;</sub> = M<sub>Pl</sub> ~ 10<sup>28</sup> eV
        (invisible, as shown on Page 17)<br>
        &bull; As a<sub>0</sub> increases, m<sub>&phi;</sub> drops
        inversely<br>
        &bull; At a<sub>0</sub> ~ 30 &mu;m, the dilaton enters the
        sub-meV range where short-range gravity experiments become
        sensitive<br>
        &bull; The Yukawa range &lambda; = &hbar;c / m<sub>&phi;</sub>
        equals a<sub>0</sub> (the Compton wavelength matches the
        extra-dimension size)
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown("#### Key Radii")

    radii_data = [
        {
            "a_0": "l_Pl (1.6e-35 m)",
            "m_phi": "M_Pl (1.2e28 eV)",
            "lambda": "1.6e-35 m",
            "Classification": "Invisible",
        },
        {
            "a_0": "30 um (3e-5 m)",
            "m_phi": "~7 meV",
            "lambda": "~30 um",
            "Classification": "Sub-mm testable",
        },
        {
            "a_0": "0.1 mm (1e-4 m)",
            "m_phi": "~2 meV",
            "lambda": "~0.1 mm",
            "Classification": "Threshold",
        },
        {
            "a_0": "1 mm (1e-3 m)",
            "m_phi": "~0.2 meV",
            "lambda": "~1 mm",
            "Classification": "Excluded",
        },
    ]

    st.dataframe(
        pd.DataFrame(radii_data),
        use_container_width=True,
        hide_index=True,
    )

# ---------------------------------------------------------------------------
# Section B: Mass vs Radius
# ---------------------------------------------------------------------------
with st.expander("B. Mass vs Radius", expanded=True):

    if _core_available and _mass_data:
        try:
            from app.components.charts import radius_mass_chart  # noqa: E402
            testable_window = None
            if _window:
                testable_window = _window
            fig_mass = radius_mass_chart(_mass_data, testable_window=testable_window)
            st.plotly_chart(fig_mass, use_container_width=True)
        except ImportError:
            st.info(
                "Chart function radius_mass_chart not yet available "
                "in charts module."
            )
        except Exception as exc:
            st.warning(f"Mass vs radius chart error: {exc}")

        st.markdown("")

        if _window:
            col_bw1, col_bw2 = st.columns(2)
            with col_bw1:
                a_min = _window.get("a_0_min")
                st.metric(
                    label="Testable window lower bound",
                    value=fmt_decimal(a_min, sig_figs=3) + " m" if a_min is not None else "N/A",
                )
            with col_bw2:
                a_max = _window.get("a_0_max_testable")
                st.metric(
                    label="Testable window upper bound",
                    value=fmt_decimal(a_max, sig_figs=3) + " m" if a_max is not None else "N/A",
                )
    else:
        st.markdown(
            """
            <div class="step-card">
            <b>Mass-radius relationship:</b><br><br>
            On a log-log plot, m<sub>&phi;</sub> vs a<sub>0</sub> is a
            straight line with slope &minus;1:<br><br>
            <code>log(m<sub>&phi;</sub>) = log(M<sub>Pl</sub>
            l<sub>Pl</sub>) &minus; log(a<sub>0</sub>)</code><br><br>
            The testable window sits where this line crosses the
            sensitivity bands of sub-mm gravity experiments (Eot-Wash,
            Irvine, Indiana) at a<sub>0</sub> ~ 10&ndash;100 &mu;m.
            </div>
            """,
            unsafe_allow_html=True,
        )

# ---------------------------------------------------------------------------
# Section C: Screening Landscape
# ---------------------------------------------------------------------------
with st.expander("C. Screening Landscape", expanded=True):

    st.markdown(
        """
        <div class="step-card">
        <b>Screening amplitude vs a<sub>0</sub>:</b><br><br>
        The effective fifth-force coupling at distance r is:<br><br>
        <center>&alpha;<sub>eff</sub>(r, a<sub>0</sub>) = 2&beta;&sup2;
        &middot; exp(&minus;r / &lambda;)</center><br>
        where &beta; = 1/&radic;6 and &lambda; = a<sub>0</sub>
        (the Compton wavelength equals the internal radius).<br><br>
        &bull; At r = 10 cm (Eot-Wash): detectable if
        a<sub>0</sub> &gt; ~30 &mu;m<br>
        &bull; At r = 1 m (laboratory): detectable if
        a<sub>0</sub> &gt; ~1 mm (already excluded)<br>
        &bull; At r = 1 AU (Cassini): detectable if
        a<sub>0</sub> &gt; ~10<sup>11</sup> m (impossible for
        compact extra dimensions)
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="theorem-card">
        <b>Conditional prediction:</b> The framework makes a
        conditional prediction: <em>IF</em> a<sub>0</sub> is in the
        sub-mm range (10&ndash;100 &mu;m), <em>THEN</em> specific
        signatures appear in next-generation short-range gravity
        experiments. The signal is a Yukawa deviation from Newton's
        law with coupling strength &alpha; = 2/6 &approx; 0.333 and
        range &lambda; = a<sub>0</sub>.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Section D: Experimental Bounds
# ---------------------------------------------------------------------------
with st.expander("D. Experimental Bounds", expanded=True):

    if _bounds_available and _bounds:
        if isinstance(_bounds, list):
            bounds_df = pd.DataFrame(_bounds)
        elif isinstance(_bounds, dict) and "experiments" in _bounds:
            bounds_df = pd.DataFrame(_bounds["experiments"])
        else:
            bounds_df = None

        if bounds_df is not None and not bounds_df.empty:
            st.dataframe(bounds_df, use_container_width=True, hide_index=True)
        else:
            st.info("Experimental bounds data format not recognized.")
    else:
        st.markdown("#### Current Experimental Landscape")

        experiments_data = [
            {
                "Experiment": "Eot-Wash (2020)",
                "Range": "> 52 um",
                "Bound on alpha": "< 0.01 at 100 um",
                "Probes": "Yukawa deviation at sub-mm",
            },
            {
                "Experiment": "Irvine (2003)",
                "Range": "> 200 um",
                "Bound on alpha": "< 1 at 200 um",
                "Probes": "Newton's law at sub-mm",
            },
            {
                "Experiment": "Indiana (2019)",
                "Range": "> 40 um",
                "Bound on alpha": "< 10 at 40 um",
                "Probes": "Planar geometry ISL test",
            },
            {
                "Experiment": "Cassini (2003)",
                "Range": "1 AU",
                "Bound on alpha": "|gamma-1| < 2.3e-5",
                "Probes": "PPN parameter at solar-system scale",
            },
            {
                "Experiment": "MICROSCOPE (2022)",
                "Range": "LEO orbit",
                "Bound on alpha": "eta < 1e-15",
                "Probes": "Equivalence principle violation",
            },
        ]

        st.dataframe(
            pd.DataFrame(experiments_data),
            use_container_width=True,
            hide_index=True,
        )

    st.markdown("")

    st.markdown(
        """
        <div class="formula-card">
        <b>What each experiment probes:</b><br><br>
        &bull; <b>Eot-Wash:</b> Torsion balance measuring torque from
        short-range Yukawa forces. Sensitive to &lambda; &gt; 52 &mu;m
        with &alpha; &lt; 0.01. This is the most sensitive probe of
        the framework's testable window.<br><br>
        &bull; <b>Cassini:</b> Shapiro time delay bounding the PPN
        parameter &gamma;. Constrains long-range modifications but
        insensitive to sub-mm physics.<br><br>
        &bull; <b>MICROSCOPE:</b> Equivalence principle test.
        Constrains composition-dependent fifth forces but the
        dilaton couples universally (no EP violation in this
        framework).
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Section E: The Conditional Prediction
# ---------------------------------------------------------------------------
with st.expander("E. The Conditional Prediction", expanded=True):

    st.markdown(
        """
        <div class="theorem-card">
        <b>Main result: a conditional, falsifiable prediction</b><br><br>
        The Alpha Ladder framework does <em>not</em> predict the value
        of a<sub>0</sub>. The internal radius is a free parameter that
        must be fixed by a UV-complete theory or measured
        experimentally.<br><br>
        However, the framework makes a sharp conditional statement:<br><br>
        &bull; <b>IF</b> a<sub>0</sub> is in the range 30&ndash;100
        &mu;m, <b>THEN</b> the dilaton has mass 2&ndash;7 meV and
        Yukawa range 30&ndash;100 &mu;m<br><br>
        &bull; The coupling strength is predicted: &alpha; = 2&beta;&sup2;
        = 1/3, independent of a<sub>0</sub><br><br>
        &bull; Next-generation Eot-Wash experiments (targeting
        &lambda; ~ 10 &mu;m) probe <em>exactly</em> this window<br><br>
        &bull; A detection at &alpha; &approx; 1/3 with &lambda; in
        the sub-mm range would be strong evidence for the framework
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("")

    st.markdown(
        """
        <div class="proof-card">
        <b>Detection vs null result:</b><br><br>
        &bull; <b>Detection scenario:</b> If a Yukawa signal is found
        at &lambda; ~ 30&ndash;100 &mu;m with &alpha; ~ 0.3, this
        would (i) confirm the extra-dimension interpretation,
        (ii) measure a<sub>0</sub> directly from &lambda;, and
        (iii) validate the predicted coupling from the dilaton-matter
        vertex.<br><br>
        &bull; <b>Null result scenario:</b> If next-generation
        experiments push the bound to &lambda; &lt; 10 &mu;m with
        &alpha; &lt; 0.01, this would exclude a<sub>0</sub> &gt;
        10 &mu;m. The framework would remain consistent (with
        a<sub>0</sub> &lt; 10 &mu;m, possibly down to l<sub>Pl</sub>),
        but the dilaton would be undetectable with current technology.
        The framework would not be falsified, but would lose its
        near-term testability.<br><br>
        &bull; <b>Falsification:</b> A detection of a Yukawa force
        with &alpha; &ne; 1/3 (at the measured &lambda;) would
        falsify the specific dilaton coupling predicted by the
        framework, even if extra dimensions exist.
        </div>
        """,
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------
st.divider()
st.caption("Radius Phenomenology | Alpha Ladder Research Dashboard")
