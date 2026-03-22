"""
Dilaton vev vs golden ratio relationship check.

Systematic search for any clean numerical relationship between the
dilaton vev phi_vev = -0.597 (from matching alpha_KK = alpha_EM) and
the golden ratio phi_golden = (1+sqrt(5))/2 = 1.618 (from the vacuum
polynomial x^2+6x+4=0).

Result: NO clean relationship found. These are mathematically independent:
- phi_vev is transcendental (involves ln, pi, alpha)
- phi_golden is algebraic (root of x^2-x-1=0)
A connection would require a deep identity linking alpha to the golden
ratio, which is what the Alpha Ladder formula claims at the level of G
(via the bridge coefficient), but not via the dilaton vev.
"""

import math


# Physical constants
ALPHA_EM = 1 / 137.035999084
PHI_VEV = (1 / 4) * math.log(4 * math.pi * ALPHA_EM)  # -0.59730485
PHI_GOLDEN = (1 + math.sqrt(5)) / 2  # 1.61803399
ALPHA = 0.0072973525693
C_0 = PHI_GOLDEN ** 2 / 2  # 1.30901699 (bridge coefficient)


def check_exponential_power_relations():
    """Check exp(k * phi_vev) = phi_golden^n for rational k and n.

    Scans k in [-12, 12] with denominators 1,2,3,4 and
    n in [-6, 6] with denominators 1,2,3.

    Returns
    -------
    dict with keys:
        matches : list of dict -- all hits within 1% (10000 ppm)
        closest : dict -- single best match
        any_clean_match : bool -- True if any hit within 100 ppm
    """
    matches = []
    best = None
    best_ppm = float("inf")

    k_values = sorted(set(
        num / den
        for num in range(-12, 13)
        for den in [1, 2, 3, 4]
        if num != 0
    ))

    n_values = sorted(set(
        num / den
        for num in range(-6, 7)
        for den in [1, 2, 3]
        if num != 0
    ))

    for k in k_values:
        val = math.exp(k * PHI_VEV)
        for n in n_values:
            target = PHI_GOLDEN ** n
            if target <= 0:
                continue
            ratio = val / target
            ppm = abs(ratio - 1) * 1e6
            if ppm < 10000:
                entry = {
                    "k": k,
                    "n": n,
                    "exp_val": val,
                    "phi_g_n": target,
                    "ratio": ratio,
                    "ppm_off": ppm,
                }
                matches.append(entry)
                if ppm < best_ppm:
                    best_ppm = ppm
                    best = entry

    return {
        "matches": sorted(matches, key=lambda x: x["ppm_off"]),
        "closest": best,
        "any_clean_match": best_ppm < 100 if best else False,
        "n_matches_within_1pct": len(matches),
        "description": (
            f"Scanned {len(k_values)} k values x {len(n_values)} n values. "
            f"Found {len(matches)} hits within 1%. "
            f"Closest: k={best['k']}, n={best['n']}, "
            f"{best['ppm_off']:.0f} ppm off."
            if best else "No matches found."
        ),
    }


def check_algebraic_combinations():
    """Check algebraic combinations of phi_vev and phi_golden.

    Tests formulas like e^(phi_vev) + 1 =? phi_golden, etc.

    Returns
    -------
    dict with keys:
        checks : list of dict
        any_match : bool
        closest : dict
    """
    e_phi = math.exp(PHI_VEV)
    e_neg = math.exp(-PHI_VEV)

    formulas = [
        ("e^(phi_vev) + 1", e_phi + 1, PHI_GOLDEN),
        ("e^(-phi_vev) - 1/phi_g", e_neg - 1 / PHI_GOLDEN, 1.0),
        ("e^(phi_vev) * phi_g", e_phi * PHI_GOLDEN, 1.0),
        ("e^(-phi_vev) / phi_g", e_neg / PHI_GOLDEN, 1.0),
        ("e^(2*phi_vev) + e^(phi_vev)", math.exp(2 * PHI_VEV) + e_phi, 1.0),
        ("phi_vev * phi_g", PHI_VEV * PHI_GOLDEN, -1.0),
        ("phi_vev^2", PHI_VEV ** 2, 1 / PHI_GOLDEN),
        ("phi_vev^2 * phi_g", PHI_VEV ** 2 * PHI_GOLDEN, 1.0),
        ("-phi_vev * 2", -PHI_VEV * 2, PHI_GOLDEN - 1 / PHI_GOLDEN),
        ("phi_vev + ln(phi_g)", PHI_VEV + math.log(PHI_GOLDEN), 0.0),
        ("1/(2*phi_vev)", 1 / (2 * PHI_VEV), -1 / PHI_GOLDEN),
        ("e^(phi_vev) + e^(2*phi_vev)", e_phi + e_phi ** 2, 1.0),
        ("-4*phi_vev / phi_g", -4 * PHI_VEV / PHI_GOLDEN, math.log(4 * math.pi * ALPHA) / PHI_GOLDEN),
    ]

    checks = []
    best = None
    best_ppm = float("inf")

    for name, val, target in formulas:
        if target == 0:
            ppm = abs(val) * 1e6
        else:
            ppm = abs(val / target - 1) * 1e6
        entry = {"formula": name, "value": val, "target": target, "ppm_off": ppm}
        checks.append(entry)
        if ppm < best_ppm:
            best_ppm = ppm
            best = entry

    return {
        "checks": checks,
        "any_match": best_ppm < 1000 if best else False,
        "closest": best,
    }


def check_bridge_vs_gauge_matching():
    """Check if the bridge coefficient C_0 relates to exp(4*phi_vev).

    C_0 = phi_golden^2/2 = 1.30902 (from vacuum polynomial)
    exp(4*phi_vev) = 4*pi*alpha = 0.09170 (from gauge matching)

    Returns
    -------
    dict with ratio, interpretation
    """
    c0 = C_0
    e4phi = math.exp(4 * PHI_VEV)

    ratio = c0 / e4phi
    product = c0 * e4phi
    ratio_to_14 = ratio / 14
    ratio_to_4pi = ratio / (4 * math.pi)

    return {
        "C_0": c0,
        "exp_4phi": e4phi,
        "C0_over_e4phi": ratio,
        "C0_times_e4phi": product,
        "ratio_close_to_14": abs(ratio_to_14 - 1) * 100,
        "ratio_close_to_4pi": abs(ratio_to_4pi - 1) * 100,
        "description": (
            f"C_0 / exp(4*phi_vev) = {ratio:.6f}. "
            f"Close to 14 (off by {abs(ratio_to_14 - 1) * 100:.1f}%). "
            f"Close to 4*pi (off by {abs(ratio_to_4pi - 1) * 100:.1f}%). "
            f"Neither is a clean match."
        ),
    }


def mathematical_incompatibility():
    """Explain why phi_vev and phi_golden cannot be simply related.

    Returns
    -------
    dict with the mathematical argument
    """
    return {
        "phi_vev_type": "transcendental",
        "phi_vev_formula": "(1/4) * ln(4 * pi * alpha_EM)",
        "phi_vev_involves": ["ln", "pi", "alpha (measured, likely transcendental)"],
        "phi_golden_type": "algebraic",
        "phi_golden_formula": "(1 + sqrt(5)) / 2",
        "phi_golden_satisfies": "x^2 - x - 1 = 0",
        "connection_would_require": (
            "A deep identity linking alpha_EM to the golden ratio. "
            "This is essentially what the Alpha Ladder formula claims "
            "at the level of G (via the bridge coefficient C_0 = phi^2/2), "
            "but not via the dilaton vev. The bridge uses phi_golden "
            "algebraically; the vev uses alpha transcendentally. "
            "These are different mathematical structures."
        ),
        "lindemann_weierstrass": (
            "By the Lindemann-Weierstrass theorem, exp(algebraic) is "
            "transcendental for nonzero algebraic arguments. Since "
            "phi_golden is algebraic, exp(k*phi_golden) is transcendental "
            "for any nonzero rational k. For phi_vev to equal phi_golden "
            "(or a power thereof), alpha would need to be an algebraic "
            "function of pi and phi_golden, which is not known and "
            "likely false."
        ),
    }


def summarize_phi_vev_golden_ratio():
    """Full summary of the phi_vev vs phi_golden search.

    Returns
    -------
    dict with all results and honest assessment
    """
    exp_power = check_exponential_power_relations()
    algebraic = check_algebraic_combinations()
    bridge_gauge = check_bridge_vs_gauge_matching()
    math_incompat = mathematical_incompatibility()

    return {
        "exponential_power_search": exp_power,
        "algebraic_combinations": algebraic,
        "bridge_vs_gauge": bridge_gauge,
        "mathematical_argument": math_incompat,
        "phi_vev": PHI_VEV,
        "phi_golden": PHI_GOLDEN,
        "relationship_found": False,
        "honest_assessment": (
            "No clean relationship exists between the dilaton vev "
            f"phi_vev = {PHI_VEV:.6f} and the golden ratio "
            f"phi_golden = {PHI_GOLDEN:.6f}. The closest exponential "
            f"match is exp({exp_power['closest']['k']}*phi_vev) vs "
            f"phi_golden^{exp_power['closest']['n']} at "
            f"{exp_power['closest']['ppm_off']:.0f} ppm -- not meaningful. "
            "These quantities are mathematically incompatible: phi_vev is "
            "transcendental (involves ln, pi, alpha) while phi_golden is "
            "algebraic (root of x^2-x-1=0). The bridge coefficient C_0 "
            "and the gauge matching are independent structures in the "
            "framework. A connection would require an unknown identity "
            "linking alpha to the golden ratio."
        ),
        "key_messages": [
            "phi_vev = transcendental (from ln, pi, alpha)",
            "phi_golden = algebraic (from x^2 - x - 1 = 0)",
            "No exp(k*phi_vev) = phi_golden^n match within 5000 ppm",
            "C_0 / exp(4*phi_vev) = 14.27 (2% off from 14, not clean)",
            "The bridge and gauge matching are independent structures",
            "Lindemann-Weierstrass theorem makes a simple link unlikely",
        ],
    }


if __name__ == "__main__":
    print("=" * 72)
    print("DILATON VEV vs GOLDEN RATIO: SYSTEMATIC SEARCH")
    print("=" * 72)
    print()
    print(f"phi_vev (dilaton)  = {PHI_VEV:.10f}")
    print(f"phi_golden         = {PHI_GOLDEN:.10f}")
    print(f"C_0 = phi_g^2/2    = {C_0:.10f}")
    print(f"e^(4*phi_vev)      = {math.exp(4 * PHI_VEV):.10f}")
    print(f"4*pi*alpha         = {4 * math.pi * ALPHA:.10f}")
    print()

    # 1. Exponential-power search
    print("-" * 72)
    print("1. EXPONENTIAL-POWER SEARCH: exp(k*phi_vev) =? phi_golden^n")
    print("-" * 72)
    ep = check_exponential_power_relations()
    if ep["matches"]:
        print(f"Found {len(ep['matches'])} hits within 1%:")
        print(f"{'k':>6s} {'n':>6s}  {'exp(k*phi_vev)':>16s}  {'phi_g^n':>16s}  {'ppm off':>10s}")
        print("-" * 60)
        for m in ep["matches"][:10]:
            print(f"{m['k']:6.2f} {m['n']:6.2f}  {m['exp_val']:16.8f}  {m['phi_g_n']:16.8f}  {m['ppm_off']:10.0f}")
        print()
        print(f"Closest: k={ep['closest']['k']}, n={ep['closest']['n']}, "
              f"{ep['closest']['ppm_off']:.0f} ppm")
    else:
        print("No matches within 1% found.")
    print(f"Clean match (< 100 ppm): {ep['any_clean_match']}")

    # 2. Algebraic combinations
    print()
    print("-" * 72)
    print("2. ALGEBRAIC COMBINATIONS")
    print("-" * 72)
    alg = check_algebraic_combinations()
    print(f"{'Formula':>35s}  {'Value':>14s}  {'Target':>14s}  {'ppm':>10s}")
    print("-" * 80)
    for c in alg["checks"]:
        marker = " <--" if c["ppm_off"] < 50000 else ""
        print(f"{c['formula']:>35s}  {c['value']:14.8f}  {c['target']:14.8f}  "
              f"{c['ppm_off']:10.0f}{marker}")

    # 3. Bridge vs gauge matching
    print()
    print("-" * 72)
    print("3. BRIDGE COEFFICIENT vs GAUGE MATCHING")
    print("-" * 72)
    bg = check_bridge_vs_gauge_matching()
    print(f"C_0 = phi_g^2/2          = {bg['C_0']:.8f}")
    print(f"exp(4*phi_vev) = 4*pi*a  = {bg['exp_4phi']:.8f}")
    print(f"C_0 / exp(4*phi_vev)     = {bg['C0_over_e4phi']:.6f}")
    print(f"  vs 14:   off by {bg['ratio_close_to_14']:.1f}%")
    print(f"  vs 4*pi: off by {bg['ratio_close_to_4pi']:.1f}%")

    # 4. Mathematical argument
    print()
    print("-" * 72)
    print("4. MATHEMATICAL INCOMPATIBILITY")
    print("-" * 72)
    mi = mathematical_incompatibility()
    print(f"phi_vev type:    {mi['phi_vev_type']}")
    print(f"phi_golden type: {mi['phi_golden_type']}")
    print()
    print("Lindemann-Weierstrass:")
    print(f"  {mi['lindemann_weierstrass']}")

    # Conclusion
    print()
    print("=" * 72)
    print("CONCLUSION")
    print("=" * 72)
    summary = summarize_phi_vev_golden_ratio()
    for msg in summary["key_messages"]:
        print(f"  - {msg}")
    print()
    print("No clean relationship found. The dilaton vev and the golden")
    print("ratio are independent mathematical structures in the framework.")
    print("=" * 72)
