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

**Prediction:** G = 6.674298 x 10^-11 m^3 kg^-1 s^-2
**CODATA 2018:** G = 6.674300 x 10^-11 m^3 kg^-1 s^-2
**Residual:** -0.31 ppm (within the 22 ppm measurement uncertainty)

No fitted parameters. Every quantity is either a measured physical constant or derived from the spacetime dimensions (d=4 spacetime, D=6 total with S^2 compactification).

## The Framework

The Alpha Ladder proposes that all coupling constants are geometric powers of alpha:

| Coupling | Formula | Accuracy |
|----------|---------|----------|
| Electromagnetic | alpha | exact (by definition) |
| Gravitational | alpha^24 * mu * (mu - sqrt(phi)*(1-alpha)) | 0.31 ppm |

The exponent 24 = d * D arises from the product of spacetime (d=4) and total (D=6) dimensions. The golden ratio phi emerges from the vacuum polynomial of the compactification geometry.

## Project Structure

```
alpha_ladder_core/   # 39 computation modules (pure Python, Decimal precision)
app/                 # Streamlit dashboard (33 pages)
tests/               # 1713 unit tests
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
| `constants.py` | CODATA 2014/2018 constants with 50-digit Decimal precision |
| `predict_g.py` | G predictions from bridge, hierarchy, and mu-structure formulas |
| `mu_structure.py` | The sqrt(phi)*(1-alpha) discovery and formula comparison |
| `unified_formula.py` | Convergence analysis between bridge and mu-structure paths |
| `residual_mapping.py` | Systematic search mapping residuals against SM constants |
| `kk_reduction.py` | Kaluza-Klein reduction from 6D to 4D |
| `vacuum_polynomial.py` | Why phi: the polynomial x^2+Dx+d=0 at d=4, D=6 |
| `testable_predictions.py` | Falsifiable predictions and experimental comparison |

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

This is active research. The numerical result is robust (verified against CODATA 2014 and 2018), but the theoretical derivation connecting Kaluza-Klein geometry to this specific formula remains incomplete. See the dashboard for detailed analysis of each theoretical angle.

## License

Apache License 2.0. See [LICENSE](LICENSE).
