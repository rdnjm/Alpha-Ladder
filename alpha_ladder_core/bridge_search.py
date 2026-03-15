"""
Search for simple mathematical expressions that match the bridge coefficient.

Refactored from legacy/bridge_search.py.
"""

from decimal import Decimal, getcontext

getcontext().prec = 50


# ---------------------------------------------------------------------------
# Unicode display helpers for expression labels
# ---------------------------------------------------------------------------

_DISPLAY_NAMES = {
    "pi": "π",
    "phi": "φ",
    "e": "e",
    "ln2": "ln2",
    "sqrt2": "√2",
    "sqrt3": "√3",
    "sqrt5": "√5",
    "1": "1",
    "2": "2",
    "3": "3",
    "4": "4",
    "5": "5",
    "6": "6",
    "7": "7",
}

_SUPERSCRIPT_MAP = str.maketrans("0123456789+-/()", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻ᐟ⁽⁾")


def _fmt_power(name, p, q):
    """Format a constant raised to a rational power using Unicode display names
    and superscript characters.

    Examples:
        _fmt_power("phi", 2, 1)  -> "φ²"
        _fmt_power("pi", -1, 3)  -> "π⁻¹ᐟ³"
        _fmt_power("sqrt2", 3, 1) -> "√2³"
    """
    display = _DISPLAY_NAMES.get(name, name)
    if q == 1:
        power_str = str(p).translate(_SUPERSCRIPT_MAP)
    else:
        power_str = f"({p}/{q})".translate(_SUPERSCRIPT_MAP)
    return f"{display}{power_str}"


def compute_target_coefficient(constants):
    """Compute the target bridge coefficient alpha_g / alpha^21.

    Returns (target_decimal, target_float) tuple.
    """
    alpha = constants.alpha
    alpha_g = constants.alpha_g
    alpha_21 = alpha ** 21
    target = alpha_g / alpha_21
    return target, float(target)


def _get_math_constants():
    """Return the dict of math constants used in searches, matching legacy script."""
    pi = Decimal("3.14159265358979323846264338327950288419716939937510")
    phi = (1 + Decimal(5).sqrt()) / 2
    e = Decimal("2.71828182845904523536028747135266249775724709369995")
    ln2 = Decimal("0.69314718055994530941723212145817656807550013436026")
    sqrt2 = Decimal(2).sqrt()
    sqrt3 = Decimal(3).sqrt()
    sqrt5 = Decimal(5).sqrt()

    return {
        "pi": pi,
        "phi": phi,
        "e": e,
        "ln2": ln2,
        "sqrt2": sqrt2,
        "sqrt3": sqrt3,
        "sqrt5": sqrt5,
        "2": Decimal(2),
        "3": Decimal(3),
        "4": Decimal(4),
        "5": Decimal(5),
        "6": Decimal(6),
        "7": Decimal(7),
        "1": Decimal(1),
    }


def search_single_constant(target_f, constants_dict=None):
    """Phase 1: Single constant with rational power p/q.

    Returns list of (err, expr, val) tuples where err < 2.0%.
    """
    if constants_dict is None:
        constants_dict = _get_math_constants()

    results = []
    for name, c in constants_dict.items():
        c_f = float(c)
        if c_f <= 0:
            continue
        for p in range(-6, 7):
            for q in range(1, 7):
                if p == 0:
                    continue
                try:
                    val = float(c) ** (p / q)
                    if val > 0:
                        err = abs(val - target_f) / target_f * 100
                        if err < 2.0:
                            expr = _fmt_power(name, p, q)
                            results.append((err, expr, val))
                except Exception:
                    pass
    return results


def search_two_constants(target_f, constants_dict=None):
    """Phase 2: A^(p/q) * B^(r/s) with rational powers.

    Returns list of (err, expr, val) tuples where err < 0.5%.
    """
    if constants_dict is None:
        constants_dict = _get_math_constants()

    results = []
    const_list = list(constants_dict.items())
    powers = [(p, q) for p in range(-4, 5) for q in range(1, 5) if p != 0]

    for i, (n1, c1) in enumerate(const_list):
        for j, (n2, c2) in enumerate(const_list):
            if j <= i:
                continue
            c1_f, c2_f = float(c1), float(c2)
            if c1_f <= 0 or c2_f <= 0:
                continue
            for (p1, q1) in powers:
                v1 = c1_f ** (p1 / q1)
                for (p2, q2) in powers:
                    try:
                        v2 = c2_f ** (p2 / q2)
                        val = v1 * v2
                        if val > 0:
                            err = abs(val - target_f) / target_f * 100
                            if err < 0.5:
                                e1 = _fmt_power(n1, p1, q1)
                                e2 = _fmt_power(n2, p2, q2)
                                expr = f"{e1} · {e2}"
                                results.append((err, expr, val))
                    except Exception:
                        pass
    return results


def search_fraction_times_constant(target_f, constants_dict=None):
    """Phase 3: (a/b) * constant^(p/q) with simple fractions.

    Returns list of (err, expr, val) tuples where err < 0.1%.
    """
    if constants_dict is None:
        constants_dict = _get_math_constants()

    results = []
    for a in range(1, 13):
        for b in range(1, 13):
            if a == b:
                continue
            frac = a / b
            for name, c in constants_dict.items():
                c_f = float(c)
                if c_f <= 0:
                    continue
                for p in range(-4, 5):
                    for q in range(1, 5):
                        if p == 0:
                            continue
                        try:
                            val = frac * (c_f ** (p / q))
                            if val > 0:
                                err = abs(val - target_f) / target_f * 100
                                if err < 0.1:
                                    ep = _fmt_power(name, p, q)
                                    expr = f"({a}/{b}) · {ep}"
                                    results.append((err, expr, val))
                        except Exception:
                            pass
    return results


def run_full_search(constants):
    """Run all three search phases, combine, sort, and deduplicate.

    Returns sorted list of top matches as (err, expr, val) tuples.
    """
    _, target_f = compute_target_coefficient(constants)
    math_consts = _get_math_constants()

    results = []
    results.extend(search_single_constant(target_f, math_consts))
    results.extend(search_two_constants(target_f, math_consts))
    results.extend(search_fraction_times_constant(target_f, math_consts))

    results.sort(key=lambda x: x[0])

    # Deduplicate by expression
    seen = set()
    deduped = []
    for err, expr, val in results:
        if expr not in seen:
            seen.add(expr)
            deduped.append((err, expr, val))

    return deduped
