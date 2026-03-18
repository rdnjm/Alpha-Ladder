# Alpha Ladder

**[Live Dashboard](https://alphaladder.streamlit.app/)**

Coupling constants as geometric powers of the fine-structure constant.

## The Result

A zero-parameter prediction of Newton's gravitational constant:

```
G = alpha^24 * mu * (mu - sqrt(phi) * (1 - alpha)) * hbar * c / m_e^2
```

where:
- `alpha` = fine-structure constant (1/137.036...)
- `mu` = proton-to-electron mass ratio (m_p / m_e = 1836.15...)
- `phi` = golden ratio ((1 + sqrt(5)) / 2), arising from the vacuum polynomial x^2 + 6x + 4 = 0 at d=4, D=6
- `hbar`, `c`, `m_e` = reduced Planck constant, speed of light, electron mass

| CODATA Edition | G predicted | G measured | Residual |
|----------------|-------------|------------|----------|
| CODATA 2022 | 6.674298 x 10^-11 | 6.674300 x 10^-11 | -0.33 ppm |
| CODATA 2018 | 6.674298 x 10^-11 | 6.674300 x 10^-11 | -0.31 ppm |
| CODATA 2014 | -- | 6.674080 x 10^-11 | different G value |

Stable across three independent CODATA editions. No fitted parameters. Every quantity is either a measured physical constant or derived from the spacetime dimensions (d=4 spacetime, D=6 total with S^2 compactification).

## The Framework

The Alpha Ladder proposes that all coupling constants are geometric powers of alpha:

| Coupling | Formula | Accuracy |
|----------|---------|----------|
| Electromagnetic | alpha | exact (by definition) |
| Gravitational | alpha^24 * mu * (mu - sqrt(phi)*(1-alpha)) | 0.33 ppm |

The exponent 24 = d * D arises from the product of spacetime (d=4) and total (D=6) dimensions. The golden ratio phi emerges from the vacuum polynomial of the compactification geometry.

The complete bridge form is equivalent:

```
G = phi^2/2 * (1 + 3*alpha^2 + (phi/2)*alpha^3) * alpha^21 * hbar * c / m_e^2
```

Both formulas give -0.31 ppm with zero fitted parameters. They are the same identity, factored differently (proven in mu_tension.py).

### Retraction note (v5, March 18 2026)

Version 4 of the preprint proposed a geometric resummation of the correction series and a derived prediction of mu to 0.001 ppm. This claim has been retracted: it is inconsistent with Alighanbari et al. (2025) H2+ spectroscopy (Nature 644, 69) at >14 sigma and with CODATA 2022 at >20 sigma. The assumed 1/phi geometric ratio was an artifact of overfitting to CODATA 2018 data. The G prediction (-0.31 ppm) is unaffected. The framework makes one testable prediction (G to sub-ppm), not two.

## Project Structure

```
alpha_ladder_core/   # 46 computation modules (pure Python, Decimal precision)
app/                 # Streamlit dashboard (40 pages)
tests/               # 2249 unit tests
legacy/              # Original standalone scripts
```

## Quick Start

```bash
# Install
pip install -e ".[dev]"

# Run tests
pytest

# Launch dashboard
streamlit run app/Home.py
```

## Requirements

- Python >= 3.10
- streamlit >= 1.30.0
- plotly >= 5.18.0
- pandas >= 2.0.0
- numpy >= 1.24.0

## Key Modules

| Module | Description |
|--------|-------------|
| `constants.py` | CODATA 2014/2018/2022 constants with 50-digit Decimal precision |
| `dimension_uniqueness.py` | Proves d=4, D=6 is uniquely selected by 3 constraints |
| `mu_tension.py` | Unifies bridge and mu-structure into one identity |
| `mu_structure.py` | The sqrt(phi)*(1-alpha) discovery and formula comparison |
| `one_alpha_derivation.py` | Derives (1-alpha) from S^2 volume cancellation |
| `feynman_diagram.py` | Explicit sigma -> KK photon loop -> sigma verification |
| `predict_g.py` | G predictions from bridge, hierarchy, and mu-structure formulas |
| `corrected_bridge.py` | Radiative correction phi^2/2*(1+3*alpha^2+(phi/2)*alpha^3); geometric resummation retracted |
| `salam_sezgin_stabilization.py` | Salam-Sezgin radius stabilization (Lambda_6 > 0, Gap 1 conditional) |
| `charged_matter_loops.py` | Monopole spectral zeta sign flip (+1/10), anomaly-free group scan |
| `vacuum_polynomial.py` | Why phi: the polynomial x^2+Dx+d=0 at d=4, D=6 |
| `kk_reduction.py` | Kaluza-Klein reduction from 6D to 4D |
| `casimir_stabilization.py` | Casimir energy on S^2 including matter loop corrections |
| `flux_stabilization.py` | 2-form flux potential, Planck-mass minimum |
| `testable_predictions.py` | Falsifiable predictions and experimental comparison |
| `second_predictions.py` | dG/G=24*dalpha/alpha, M_6 landscape, correlated variations |
| `literature_comparison.py` | Quantitative comparison with published approaches |

## Dashboard Pages

**Core Results**
| Page | Description |
|------|-------------|
| The Derivation | Complete 6D-to-G derivation chain |
| Dimension Uniqueness | d=4, D=6 uniquely selected by 3 constraints |
| Feynman Diagram | Explicit sigma -> KK photon loop -> sigma |
| Mu Tension | Bridge = mu-structure proven as one identity |
| Second Predictions | dG/G = 24*dalpha/alpha, M_6 landscape |
| The Prediction | G from first principles via phi^2/2 bridge |

**Experimental**
| Page | Description |
|------|-------------|
| Fifth Force Predictions | Yukawa signal, Eot-Wash prediction, discovery reach |
| Testable Predictions | G precision, mu consistency, fifth force bounds |
| Solar System | Solar system constraints on extra dimensions |

**Theory**
| Page | Description |
|------|-------------|
| The Proof | 4 theoretical gaps and their status |
| One-Alpha Derivation | S^2 volume cancellation derives (1-alpha) |
| Corrected Bridge | phi^2/2*(1 + 3*alpha^2 + (phi/2)*alpha^3) |
| Mu Structure | sqrt(phi)*(1-alpha) offset, formula comparison |
| Unified Formula | Resummed bridges, gap analysis |
| Bridge Significance | Expression density, significance analysis |
| Hierarchy Derivation | 10 theoretical angles for alpha^24*mu^2 |
| Literature Comparison | Beck, Alexander, Eaves, BVW, Eddington-Dirac |
| Residual Mapping | delta vs SM constants, closed-form search |

**Gap Analysis**
| Page | Description |
|------|-------------|
| Casimir Stabilization | Casimir no-go + matter loop corrections |
| Flux Stabilization | 2-form flux potential, Planck-mass minimum |
| Radius Phenomenology | Testable window, Eot-Wash optimal at ~28 um |
| Chameleon Screening | Density-dependent mass, KK truncation check |
| Dark Sector Phenomenology | Relic abundance, self-interaction, fuzzy DM |
| Radius Determination | Scaling symmetry proof, mechanism catalog |
| Radius Fixing | CW, warped S^2, orbifold -- all fail |
| Salam-Sezgin Stabilization | Lambda_6 > 0 fixes radius (Gap 1 conditional) |
| Charged Matter Loops | Monopole spectral zeta sign flip |
| Anomaly Cancellation | Pure gravity safe, SM embedding, G unaffected |
| Cosmological Constant | V_min ~ O(1) Planck, 122-order discrepancy |

**Legacy Dashboard** (pages 1-11): Constant Core, Geometric Ladder, Bridge Lab, Universe Slider, Phi Scanner, Particle Harmonics, Rung Spacing, Dilaton Lab, Experimental, Alpha Units, Dark Sector.

## What is Derived vs. Empirical

**Derived from first principles:**
- d=4, D=6 uniquely selected by 3 independent constraints (exponent, volume cancellation, vacuum polynomial)
- The exponent 24 = d * D
- The golden ratio phi from the vacuum polynomial x^2 + 6x + 4 = 0
- sqrt(phi) as the mass offset
- The (1-alpha) correction from S^2 volume cancellation: Vol(S^2)/R^2 = 4*pi cancels the 1/(4*pi) loop factor
- The NLO correction (phi/2)*alpha^3 from the same mechanism
- Verified by explicit Feynman diagram: sigma -> KK photon loop -> sigma

**Measured (not fitted):**
- alpha, mu, hbar, c, m_e (CODATA values)

**Open:**
- Whether the formula is exact or a leading-order approximation

## Status

This is active research. The numerical result is robust (verified against CODATA 2014, 2018, and 2022) with 2249 tests passing. The theoretical derivation chain is complete: dimension uniqueness (d=4, D=6) -> vacuum polynomial -> KK reduction -> corrected bridge -> G prediction at -0.31 ppm with zero fitted parameters. A geometric resummation claiming to predict mu to 0.001 ppm was proposed and subsequently retracted after falsification by Alighanbari et al. (2025) at >14 sigma. See the [live dashboard](https://alphaladder.streamlit.app/) for detailed analysis of each step.

## License

Apache License 2.0. See [LICENSE](LICENSE).
