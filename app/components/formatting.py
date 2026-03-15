"""
Decimal display helpers for the Alpha Ladder dashboard.

Provides consistent formatting for scientific values, percentages,
sigma deviations, and quality-based color coding throughout the UI.
"""

from decimal import Decimal, InvalidOperation


def fmt_decimal(val, sig_figs=6):
    """Format a Decimal or float value in scientific notation with given significant figures.

    Parameters
    ----------
    val : Decimal, float, or str
        The numeric value to format.
    sig_figs : int
        Number of significant figures to display (default 6).

    Returns
    -------
    str
        Formatted string in scientific notation, e.g. "7.29735e-03".
    """
    if val is None:
        return "N/A"
    try:
        f = float(val)
    except (TypeError, ValueError, InvalidOperation):
        return str(val)
    if f == 0:
        return "0.000000e+00"
    fmt_str = f".{sig_figs - 1}e"
    return format(f, fmt_str)


def fmt_percent(val, decimals=4):
    """Format a value as a percentage string.

    Parameters
    ----------
    val : Decimal, float, or str
        The numeric value (already as a fraction, e.g. 0.05 for 5%).
    decimals : int
        Number of decimal places in the percentage (default 4).

    Returns
    -------
    str
        Formatted percentage string, e.g. "5.0000%".
    """
    if val is None:
        return "N/A"
    try:
        f = float(val) * 100
    except (TypeError, ValueError, InvalidOperation):
        return str(val)
    return f"{f:.{decimals}f}%"


def fmt_sigma(sigma, direction=None):
    """Format a sigma deviation string.

    Parameters
    ----------
    sigma : float or Decimal
        The magnitude of the sigma deviation.
    direction : str or None
        "+" or "-" to indicate direction. If None, uses the sign of sigma.

    Returns
    -------
    str
        Formatted string like "+2.3 sigma" or "-1.1 sigma".
    """
    if sigma is None:
        return "N/A"
    try:
        f = float(sigma)
    except (TypeError, ValueError, InvalidOperation):
        return str(sigma)
    if direction is None:
        direction = "+" if f >= 0 else "-"
        f = abs(f)
    return f"{direction}{f:.1f} σ"


def color_by_quality(err_pct):
    """Return a CSS color based on error percentage quality thresholds.

    Parameters
    ----------
    err_pct : float or Decimal
        The error percentage (e.g. 0.05 means 0.05%).

    Returns
    -------
    str
        CSS color string:
        - green (#00c853) for < 0.1%
        - yellow (#ffd600) for < 1%
        - orange (#ff6d00) for < 5%
        - red (#d50000) for >= 5%
    """
    if err_pct is None:
        return "#9e9e9e"
    try:
        f = float(err_pct)
    except (TypeError, ValueError, InvalidOperation):
        return "#9e9e9e"
    if f < 0.1:
        return "#00c853"
    elif f < 1.0:
        return "#ffd600"
    elif f < 5.0:
        return "#ff6d00"
    else:
        return "#d50000"
