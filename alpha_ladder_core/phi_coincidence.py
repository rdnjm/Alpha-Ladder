"""
Scan for phi coincidences in precision-measured Standard Model parameters.
Refactored from legacy/phi_coincidence_scan.py.
"""

import math

phi = (1 + math.sqrt(5)) / 2
pi = math.pi


def get_sm_parameters() -> dict:
    """Return dict of parameter name to value (all 14 precision SM parameters)."""
    return {
        "sin\u00b2(\u03b8_W) (MS-bar, MZ)":  0.23122,
        "sin\u00b2(\u03b8_W) (on-shell)":    0.22337,
        "Cabibbo angle sin(θ_C)":       0.22500,
        "Weinberg angle θ_W (rad)":     math.asin(math.sqrt(0.23122)),
        "mₚ / mₑ":                     1836.15267343,
        "m_μ / mₑ":                     206.7682830,
        "m_τ / mₑ":                     3477.23,
        "m_W / m_Z":                    0.88147,
        "m_H / m_Z":                    1.37186,
        "m_H / m_W":                    1.55636,
        "\u03b1\u209b(MZ)":                  0.1180,
        "\u03b1(MZ) (EM at Z pole)":     1 / 127.951,
        "G_F (Fermi, reduced)":         1.1663788e-5,
        "Euler-Mascheroni gamma":       0.5772156649,
    }


def scan_parameter(pval: float, alpha_em: float) -> tuple:
    """
    Run all 3 search phases for a single parameter value.
    Returns (best_err, best_expr, best_val).

    Phase 1: phi^(p/q)
    Phase 2: (a/b) * phi^(p/q)
    Phase 3: alpha^n * phi^(p/q)

    Ranges are IDENTICAL to legacy.
    """
    if pval <= 0:
        return (float('inf'), "", 0)

    best_err = float('inf')
    best_expr = ""
    best_val = 0

    # Phase 1: phi^(p/q)
    for p in range(-8, 9):
        for q in range(1, 7):
            if p == 0:
                continue
            candidate = phi ** (p / q)
            err = abs(candidate - pval) / abs(pval) * 100
            if err < best_err:
                best_err = err
                best_val = candidate
                best_expr = f"phi^({p}/{q})" if q != 1 else f"phi^{p}"

    # Phase 2: (a/b) * phi^(p/q)
    for a in range(1, 13):
        for b in range(1, 13):
            if a == b:
                continue
            frac = a / b
            for p in range(-6, 7):
                for q in range(1, 5):
                    if p == 0:
                        candidate = frac
                    else:
                        candidate = frac * phi ** (p / q)
                    err = abs(candidate - pval) / abs(pval) * 100
                    if err < best_err:
                        best_err = err
                        best_val = candidate
                        ep = f"phi^({p}/{q})" if q != 1 else f"phi^{p}"
                        best_expr = f"({a}/{b}) * {ep}" if p != 0 else f"{a}/{b}"

    # Phase 3: alpha^n * phi^(p/q)
    for n in range(-3, 4):
        if n == 0:
            continue
        for p in range(-4, 5):
            for q in range(1, 4):
                candidate = (alpha_em ** n) * (phi ** (p / q))
                err = abs(candidate - pval) / abs(pval) * 100
                if err < best_err:
                    best_err = err
                    best_val = candidate
                    ep = (f"phi^({p}/{q})" if (q != 1 and p != 0)
                          else (f"phi^{p}" if p != 0 else ""))
                    an = f"alpha^{n}"
                    best_expr = f"{an} * {ep}" if ep else an

    return (best_err, best_expr, best_val)


def run_full_scan(constants: dict) -> list:
    """
    Run the phi coincidence scan over all 14 SM parameters.
    Returns list of dicts with keys: param_name, param_val, best_err, best_expr, best_val.
    Sorted by parameter order (insertion order of get_sm_parameters).
    """
    alpha_em = float(constants.alpha) if hasattr(constants, 'alpha') else constants["alpha"]
    params = get_sm_parameters()
    results = []
    for pname, pval in params.items():
        best_err, best_expr, best_val = scan_parameter(pval, alpha_em)
        results.append({
            "param_name": pname,
            "param_val": pval,
            "best_err": best_err,
            "best_expr": best_expr,
            "best_val": best_val,
        })
    return results


def deep_dive_weinberg(constants: dict) -> list:
    """
    Deep dive into sin\u00b2(\u03b8_W) candidates.
    Returns list of (err_pct, name, value, inside_bool) sorted by err_pct ascending.
    IDENTICAL candidate list to legacy.
    """
    alpha_em = float(constants.alpha) if hasattr(constants, 'alpha') else constants["alpha"]
    sw2 = 0.23122

    candidates_sw2 = {
        "1/\u03c6\u00b3":                      1 / phi**3,
        "3/(4\u03c6\u00b3)":                  3 / (4 * phi**3),
        "(\u03c6\u22121)/\u03c6\u00b2":                (phi - 1) / phi**2,
        "1/(\u03c6\u00b2 + \u03c6)":              1 / (phi**2 + phi),
        "\u03c6\u00b2 / (2\u03c0 + \u03c6)":         phi**2 / (2 * pi + phi),
        "1/(2\u03c6 + 1)":                1 / (2 * phi + 1),
        "1/\u03c6\u00b3 \u00b7 (1 \u2212 \u03b1)":        (1 / phi**3) * (1 - alpha_em),
        "(3 \u2212 \u03c6) / (4 + \u03c6)":        (3 - phi) / (4 + phi),
        "(\u03c6 \u2212 1)\u00b3":                  (phi - 1)**3,
        "1/(\u03c6 + \u03c0)":                 1 / (phi + pi),
        "\u03c0/(4\u03c6\u00b2 \u2212 1)":            pi / (4 * phi**2 - 1),
        "\u03c6/(\u03c6\u00b2 + \u03c0)":            phi / (phi**2 + pi),
    }

    ranked = []
    for name, val in candidates_sw2.items():
        err = abs(val - sw2) / sw2 * 100
        inside = abs(val - sw2) < 0.00004
        ranked.append((err, name, val, inside))

    ranked.sort()
    return ranked


def compute_21_connection() -> dict:
    """
    Compute the '21 connection' -- Riemann component counts and SM group dimensions.
    Pure computation, no constants needed.
    """
    # Riemann tensor in 4D
    n = 4
    riemann_4d = n**2 * (n**2 - 1) // 12  # = 20

    # Triangular number T(6)
    triangular_6 = 6 * 7 // 2  # = 21

    # C(7,2)
    c_7_2 = 7 * 6 // 2  # = 21

    # SO(7) dimension
    so7_dim = 7 * 6 // 2  # = 21

    # Symmetric 6x6 matrix components
    sym_6x6 = 6 * 7 // 2  # = 21

    # SM gauge group dimensions
    sm_dims = {
        "SU(3)": 8,
        "SU(2)": 3,
        "U(1)": 1,
        "total_gauge_bosons": 12,
    }

    # Symmetric rank-2 tensor components for D=3..7
    sym_tensor = {}
    for d in range(3, 8):
        comps = d * (d + 1) // 2
        sym_tensor[d] = comps

    return {
        "riemann_4d": riemann_4d,
        "triangular_6": triangular_6,
        "C_7_2": c_7_2,
        "SO7_dim": so7_dim,
        "symmetric_6x6": sym_6x6,
        "sm_gauge_dims": sm_dims,
        "symmetric_tensor_by_dim": sym_tensor,
    }
