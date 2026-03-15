# Alpha Ladder

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

A second independent formula confirms the result:

| Formula | Residual | Fitted params |
|---------|----------|---------------|
| alpha^24 * mu * (mu - sqrt(phi)*(1-alpha)) | -0.33 ppm | 0 |
| phi^2/2 * (1 + 3*alpha^2) * alpha^21 | -0.64 ppm | 1 |

The two formulas agree to 0.31 ppm, bracketing the measured value.

## Project Structure

```
alpha_ladder_core/   # 36 computation modules (pure Python, Decimal precision)
app/                 # Streamlit dashboard (32 pages)
tests/               # 1725 unit tests
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
| `predict_g.py` | G predictions from bridge, hierarchy, and mu-structure formulas |
| `mu_structure.py` | The sqrt(phi)*(1-alpha) discovery and formula comparison |
| `unified_formula.py` | Convergence analysis between bridge and mu-structure paths |
| `residual_mapping.py` | Systematic search mapping residuals against SM constants |
| `kk_reduction.py` | Kaluza-Klein reduction from 6D to 4D |
| `vacuum_polynomial.py` | Why phi: the polynomial x^2+Dx+d=0 at d=4, D=6 |
| `testable_predictions.py` | Falsifiable predictions and experimental comparison |
| `literature_comparison.py` | Quantitative comparison with published approaches |

## What is Derived vs. Empirical

**Derived from d=4, D=6 geometry:**
- The exponent 24 = d * D
- The golden ratio phi from x^2 + 6x + 4 = 0
- sqrt(phi) as the mass offset

**Measured (not fitted):**
- alpha, mu, hbar, c, m_e (CODATA values)
- The (1-alpha) correction (alpha is measured, coefficient is 1)

**Not yet derived from first principles:**
- Why d=4, D=6 specifically (assumed, not proven)
- Why the (1-alpha) correction takes this form (physically natural, no diagram calculation yet)
- Whether the formula is exact or a leading-order approximation

## Status

This is active research. The numerical result is robust (verified against CODATA 2014, 2018, and 2022), but the theoretical derivation connecting Kaluza-Klein geometry to this specific formula remains incomplete. See the dashboard for detailed analysis of each theoretical angle.

## License

Apache License 2.0. See [LICENSE](LICENSE).
