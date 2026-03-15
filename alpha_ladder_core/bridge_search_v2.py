"""
Systematic search for the bridge coefficient C in alpha_g = C * alpha^21.

Given the Alpha Ladder relation:
    alpha_g = C * alpha^21
    G = C * alpha^21 * hbar * c / m_e^2

the exact C from measurement is:
    C_exact = G * m_e^2 / (alpha^21 * hbar * c) = alpha_g / alpha^21

This module performs an exhaustive search for mathematical and physical
expressions that reproduce C_exact to sub-ppm precision, testing:

1. Pure math expressions built from pi, phi, e, sqrt(N), ln(2), and
   small integer ratios.
2. Physical expressions using alpha and mu = m_p/m_e with integer powers.
3. Hybrid expressions combining both.

The central question: is there a closed-form expression for C that
(a) matches to sub-ppm, (b) has theoretical motivation, and (c) reveals
the structural relationship between the phi^2/2 bridge and the
alpha^24 * mu^2 hierarchy formula?
"""

import math
from decimal import Decimal, getcontext
from math import gcd

from alpha_ladder_core.constants import get_constants, DEFAULT_EDITION

getcontext().prec = 50


# ---------------------------------------------------------------------------
# 1. Exact bridge coefficient
# ---------------------------------------------------------------------------

def compute_exact_bridge(constants=None):
    """Compute the exact bridge coefficient C = alpha_g / alpha^21
    using Decimal arithmetic.

    Parameters
    ----------
    constants : SimpleNamespace, optional
        Physical constants from get_constants(). Uses CODATA 2018 if None.

    Returns
    -------
    dict
        C_exact: Decimal value, C_float: float value.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    alpha = constants.alpha
    alpha_g = constants.alpha_g
    alpha_21 = alpha ** 21
    C_exact = alpha_g / alpha_21

    return {
        "C_exact": C_exact,
        "C_float": float(C_exact),
    }


# ---------------------------------------------------------------------------
# Mathematical constant pool
# ---------------------------------------------------------------------------

def _math_pool():
    """Return dict of named mathematical constants as floats."""
    phi = (1 + math.sqrt(5)) / 2
    return {
        "pi": math.pi,
        "phi": phi,
        "e": math.e,
        "sqrt2": math.sqrt(2),
        "sqrt3": math.sqrt(3),
        "sqrt5": math.sqrt(5),
        "ln2": math.log(2),
    }


# ---------------------------------------------------------------------------
# 2. Simple expression search
# ---------------------------------------------------------------------------

def search_simple_expressions(constants=None, max_ppm=1000):
    """Search for pure-math and physical expressions matching C_exact.

    Generates candidates from:
    - p/q * math_const^n  (p,q in 1..12, n in -3..3)
    - math_const_A^a * math_const_B^b  (a,b in -3..3)
    - p/q * math_const_A * math_const_B
    - alpha^k * mu^j  (integer k,j)
    - alpha^k * mu^j * (simple math expression)

    Parameters
    ----------
    constants : SimpleNamespace, optional
        Uses CODATA 2018 if None.
    max_ppm : float
        Maximum |residual_ppm| to include in results.

    Returns
    -------
    list of dict
        Each dict: name, value, residual_ppm, has_fitted_params, category.
        Sorted by |residual_ppm|.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    bridge = compute_exact_bridge(constants)
    C_exact = bridge["C_float"]

    alpha = float(constants.alpha)
    mu = float(constants.m_p / constants.m_e)
    mc = _math_pool()

    results = []
    seen_names = set()

    def _add(name, value, has_fitted, category):
        if name in seen_names:
            return
        if value is None or value <= 0 or not math.isfinite(value):
            return
        ppm = (value - C_exact) / C_exact * 1e6
        if abs(ppm) < max_ppm:
            seen_names.add(name)
            results.append({
                "name": name,
                "value": value,
                "residual_ppm": ppm,
                "has_fitted_params": has_fitted,
                "category": category,
            })

    # --- Phase A: p/q * math_const^n ---
    for cname, cval in mc.items():
        for n in range(-3, 4):
            if n == 0:
                continue
            try:
                powered = cval ** n
            except Exception:
                continue
            if not math.isfinite(powered) or powered <= 0:
                continue
            # bare constant^n
            _add(f"{cname}^{n}", powered, True, "pure_math")
            # p/q * constant^n (reduce to lowest terms)
            seen_ratios = set()
            for p in range(1, 13):
                for q in range(1, 13):
                    if p == q:
                        continue
                    g = gcd(p, q)
                    pr, qr = p // g, q // g
                    if (pr, qr) in seen_ratios:
                        continue
                    seen_ratios.add((pr, qr))
                    val = (p / q) * powered
                    _add(f"({pr}/{qr})*{cname}^{n}", val, True, "pure_math")

    # --- Phase B: math_const_A^a * math_const_B^b ---
    mc_items = list(mc.items())
    for i, (n1, v1) in enumerate(mc_items):
        for j, (n2, v2) in enumerate(mc_items):
            if j <= i:
                continue
            for a in range(-3, 4):
                if a == 0:
                    continue
                try:
                    va = v1 ** a
                except Exception:
                    continue
                if not math.isfinite(va) or va <= 0:
                    continue
                for b in range(-3, 4):
                    if b == 0:
                        continue
                    try:
                        vb = v2 ** b
                        val = va * vb
                    except Exception:
                        continue
                    _add(f"{n1}^{a}*{n2}^{b}", val, True, "pure_math")

    # --- Phase C: p/q * math_A * math_B (power 1 only, quick) ---
    for i, (n1, v1) in enumerate(mc_items):
        for j, (n2, v2) in enumerate(mc_items):
            if j <= i:
                continue
            base = v1 * v2
            seen_ratios_c = set()
            for p in range(1, 13):
                for q in range(1, 13):
                    if p == q:
                        continue
                    g = gcd(p, q)
                    pr, qr = p // g, q // g
                    if (pr, qr) in seen_ratios_c:
                        continue
                    seen_ratios_c.add((pr, qr))
                    _add(f"({pr}/{qr})*{n1}*{n2}", (p / q) * base, True, "pure_math")

    # --- Phase D: alpha^k * mu^j ---
    for k in range(-5, 6):
        for j in range(-3, 4):
            if k == 0 and j == 0:
                continue
            try:
                val = (alpha ** k) * (mu ** j) if j != 0 else alpha ** k
                if k == 0:
                    val = mu ** j
            except Exception:
                continue
            label = ""
            if k != 0:
                label += f"alpha^{k}"
            if j != 0:
                if label:
                    label += "*"
                label += f"mu^{j}"
            _add(label, val, False, "physical")

    # --- Phase E: alpha^k * mu^j * simple_math ---
    simple_math = {
        "phi": mc["phi"],
        "phi^2": mc["phi"] ** 2,
        "1/phi": 1 / mc["phi"],
        "phi/2": mc["phi"] / 2,
        "phi^2/2": mc["phi"] ** 2 / 2,
        "pi": mc["pi"],
        "pi/2": mc["pi"] / 2,
        "pi/3": mc["pi"] / 3,
        "pi/4": mc["pi"] / 4,
        "2*pi": 2 * mc["pi"],
        "(5/12)*pi": (5 / 12) * mc["pi"],
        "e": mc["e"],
        "sqrt(e)": math.sqrt(mc["e"]),
        "sqrt(e)/cbrt(2)": math.sqrt(mc["e"]) / (2 ** (1 / 3)),
        "e/2": mc["e"] / 2,
        "1/e": 1 / mc["e"],
        "sqrt2": mc["sqrt2"],
        "sqrt3": mc["sqrt3"],
        "sqrt5": mc["sqrt5"],
        "ln2": mc["ln2"],
        "1/ln2": 1 / mc["ln2"],
    }
    for mname, mval in simple_math.items():
        for k in range(-5, 6):
            for j in range(-3, 4):
                if k == 0 and j == 0:
                    # pure math -- already covered above
                    _add(mname, mval, True, "pure_math")
                    continue
                try:
                    phys = 1.0
                    if k != 0:
                        phys *= alpha ** k
                    if j != 0:
                        phys *= mu ** j
                    val = mval * phys
                except Exception:
                    continue
                parts = [mname]
                if k != 0:
                    parts.append(f"alpha^{k}")
                if j != 0:
                    parts.append(f"mu^{j}")
                _add("*".join(parts), val, True, "mixed")

    # --- Phase F: p/q * alpha^k * mu^j for small ratios ---
    for p in range(1, 7):
        for q in range(1, 7):
            if p == q:
                continue
            for k in range(-5, 6):
                for j in range(-3, 4):
                    if k == 0 and j == 0:
                        continue
                    try:
                        val = (p / q)
                        if k != 0:
                            val *= alpha ** k
                        if j != 0:
                            val *= mu ** j
                    except Exception:
                        continue
                    parts = [f"({p}/{q})"]
                    if k != 0:
                        parts.append(f"alpha^{k}")
                    if j != 0:
                        parts.append(f"mu^{j}")
                    _add("*".join(parts), val, False, "physical")

    results.sort(key=lambda r: abs(r["residual_ppm"]))
    return results


# ---------------------------------------------------------------------------
# 3. Alpha-mu bridge search
# ---------------------------------------------------------------------------

def search_alpha_mu_bridges(constants=None, max_ppm=1000):
    """Search C = alpha^a * mu^b * (p/q) systematically.

    The total formula becomes:
        alpha_g = (p/q) * alpha^(21+a) * mu^b

    Parameters
    ----------
    constants : SimpleNamespace, optional
        Uses CODATA 2018 if None.
    max_ppm : float
        Maximum |residual_ppm| to include.

    Returns
    -------
    list of dict
        Each: a, b, p, q, total_exponent, mu_exponent, residual_ppm,
        formula_str. Sorted by |residual_ppm|.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    bridge = compute_exact_bridge(constants)
    C_exact = bridge["C_float"]

    alpha = float(constants.alpha)
    mu = float(constants.m_p / constants.m_e)

    # Build reduced ratio set: all p/q with p,q in 1..6, reduced to lowest terms
    ratios = set()
    for p in range(1, 7):
        for q in range(1, 7):
            g = gcd(p, q)
            ratios.add((p // g, q // g))
    # Also include 1/q for q up to 12
    for q in range(1, 13):
        ratios.add((1, q))

    results = []
    seen_formulas = set()
    for a in range(-5, 6):
        for b in range(-3, 4):
            for (p, q) in ratios:
                try:
                    C_cand = (p / q)
                    if a != 0:
                        C_cand *= alpha ** a
                    if b != 0:
                        C_cand *= mu ** b
                except Exception:
                    continue
                if C_cand <= 0 or not math.isfinite(C_cand):
                    continue
                ppm = (C_cand - C_exact) / C_exact * 1e6
                if abs(ppm) > max_ppm:
                    continue

                total_exp = 21 + a
                parts = []
                if p != q:
                    parts.append(f"({p}/{q})")
                parts.append(f"alpha^{total_exp}")
                if b != 0:
                    parts.append(f"mu^{b}")
                formula = " * ".join(parts)

                if formula in seen_formulas:
                    continue
                seen_formulas.add(formula)

                results.append({
                    "a": a,
                    "b": b,
                    "p": p,
                    "q": q,
                    "total_exponent": total_exp,
                    "mu_exponent": b,
                    "residual_ppm": ppm,
                    "formula_str": formula,
                })

    results.sort(key=lambda r: abs(r["residual_ppm"]))
    return results


# ---------------------------------------------------------------------------
# 4. Hybrid bridge search
# ---------------------------------------------------------------------------

def search_hybrid_bridges(constants=None, max_ppm=50):
    """Search C = math_expr * alpha^a * mu^b.

    Tests known mathematical constants multiplied by alpha/mu powers.

    Parameters
    ----------
    constants : SimpleNamespace, optional
        Uses CODATA 2018 if None.
    max_ppm : float
        Maximum |residual_ppm| to include.

    Returns
    -------
    list of dict
        Each: math_expr, a, b, total_exponent, mu_exponent,
        residual_ppm, formula_str, value. Sorted by |residual_ppm|.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    bridge = compute_exact_bridge(constants)
    C_exact = bridge["C_float"]

    alpha = float(constants.alpha)
    mu = float(constants.m_p / constants.m_e)
    mc = _math_pool()
    phi = mc["phi"]

    math_exprs = {
        "phi": phi,
        "phi^2": phi ** 2,
        "1/phi": 1 / phi,
        "phi/2": phi / 2,
        "phi^2/2": phi ** 2 / 2,
        "pi": mc["pi"],
        "pi/2": mc["pi"] / 2,
        "pi/3": mc["pi"] / 3,
        "pi/4": mc["pi"] / 4,
        "2*pi": 2 * mc["pi"],
        "(5/12)*pi": (5 / 12) * mc["pi"],
        "e": mc["e"],
        "sqrt(e)": math.sqrt(mc["e"]),
        "sqrt(e)/cbrt(2)": math.sqrt(mc["e"]) / (2 ** (1 / 3)),
        "e/2": mc["e"] / 2,
        "1/e": 1 / mc["e"],
        "sqrt2": mc["sqrt2"],
        "sqrt3": mc["sqrt3"],
        "sqrt5": mc["sqrt5"],
        "ln2": mc["ln2"],
        "1/ln2": 1 / mc["ln2"],
    }

    results = []
    seen = set()

    for mname, mval in math_exprs.items():
        for a in range(-5, 6):
            for b in range(-3, 4):
                try:
                    C_cand = mval
                    if a != 0:
                        C_cand *= alpha ** a
                    if b != 0:
                        C_cand *= mu ** b
                except Exception:
                    continue
                if C_cand <= 0 or not math.isfinite(C_cand):
                    continue
                ppm = (C_cand - C_exact) / C_exact * 1e6
                if abs(ppm) > max_ppm:
                    continue

                total_exp = 21 + a
                parts = [mname]
                if a != 0:
                    parts.append(f"alpha^{a}")
                if b != 0:
                    parts.append(f"mu^{b}")
                formula = " * ".join(parts)

                if formula in seen:
                    continue
                seen.add(formula)

                results.append({
                    "math_expr": mname,
                    "a": a,
                    "b": b,
                    "total_exponent": total_exp,
                    "mu_exponent": b,
                    "residual_ppm": ppm,
                    "formula_str": formula,
                    "value": C_cand,
                })

    results.sort(key=lambda r: abs(r["residual_ppm"]))
    return results


# ---------------------------------------------------------------------------
# 5. Summary
# ---------------------------------------------------------------------------

def summarize_bridge_search(constants=None):
    """Run all searches and summarize results.

    Parameters
    ----------
    constants : SimpleNamespace, optional
        Uses CODATA 2018 if None.

    Returns
    -------
    dict
        C_exact, C_float, best_pure_math, best_physical, best_hybrid,
        top_20, key_insight.
    """
    if constants is None:
        constants = get_constants("CODATA 2018")

    bridge = compute_exact_bridge(constants)

    simple = search_simple_expressions(constants, max_ppm=1000)
    alpha_mu = search_alpha_mu_bridges(constants, max_ppm=1000)
    hybrid = search_hybrid_bridges(constants, max_ppm=1000)

    # Categorize simple results
    pure_math = [r for r in simple if r["category"] == "pure_math"]
    physical = [r for r in simple if r["category"] == "physical"]

    best_pure_math = pure_math[0] if pure_math else None
    best_physical = physical[0] if physical else None

    # Best hybrid: must actually mix math + physics (a != 0 or b != 0)
    hybrid_with_physics = [r for r in hybrid if r["a"] != 0 or r["b"] != 0]
    best_hybrid = hybrid_with_physics[0] if hybrid_with_physics else None

    # Merge all results for top-20
    all_results = []
    seen_names = set()

    for r in simple:
        key = r["name"]
        if key not in seen_names:
            seen_names.add(key)
            all_results.append({
                "name": r["name"],
                "value": r["value"],
                "residual_ppm": r["residual_ppm"],
                "category": r["category"],
                "has_fitted_params": r["has_fitted_params"],
            })

    for r in alpha_mu:
        key = r["formula_str"]
        if key not in seen_names:
            seen_names.add(key)
            all_results.append({
                "name": r["formula_str"],
                "value": None,  # not stored in alpha_mu results
                "residual_ppm": r["residual_ppm"],
                "category": "physical",
                "has_fitted_params": False,
            })

    for r in hybrid:
        key = r["formula_str"]
        if key not in seen_names:
            seen_names.add(key)
            all_results.append({
                "name": r["formula_str"],
                "value": r["value"],
                "residual_ppm": r["residual_ppm"],
                "category": "mixed",
                "has_fitted_params": True,
            })

    all_results.sort(key=lambda r: abs(r["residual_ppm"]))
    top_20 = all_results[:20]

    # Generate insight
    winners = []
    if best_pure_math:
        winners.append(
            f"Best pure math: {best_pure_math['name']} "
            f"({best_pure_math['residual_ppm']:.2f} ppm)"
        )
    if best_physical:
        winners.append(
            f"Best physical (no fitted params): {best_physical['name']} "
            f"({best_physical['residual_ppm']:.2f} ppm)"
        )
    if best_hybrid:
        winners.append(
            f"Best hybrid: {best_hybrid['formula_str']} "
            f"({best_hybrid['residual_ppm']:.2f} ppm)"
        )

    # Check for sub-ppm results
    sub_ppm = [r for r in all_results if abs(r["residual_ppm"]) < 1.0]
    if sub_ppm:
        insight = (
            f"Found {len(sub_ppm)} sub-ppm expression(s). "
            + "; ".join(winners)
        )
    else:
        insight = (
            "No sub-ppm expression found in search space. "
            + "; ".join(winners)
        )

    return {
        "C_exact": bridge["C_exact"],
        "C_float": bridge["C_float"],
        "best_pure_math": best_pure_math,
        "best_physical": best_physical,
        "best_hybrid": best_hybrid,
        "top_20": top_20,
        "key_insight": insight,
    }
