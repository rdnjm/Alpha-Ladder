"""
Proton-to-electron mass ratio vs dilaton vev relation check.

Investigates whether the dilaton vev phi_vev = -0.597 (from matching
alpha_KK = alpha_EM) has any clean numerical relationship to the
proton-to-electron mass ratio mu = m_p/m_e = 1836.15.

Result: NO clean relationship found. The dilaton vev and mu are
independent structures in the Alpha Ladder framework.
- phi_vev comes from gauge physics (alpha matching)
- mu comes from mass physics (Higgs/Yukawa sector)
- In the braneworld picture, fermion masses are brane-localized
  and do not depend on the S^2 geometry or dilaton vev.
"""

import math


# Physical constants
ALPHA_EM = 1 / 137.035999084
PHI_VEV = (1 / 4) * math.log(4 * math.pi * ALPHA_EM)  # -0.59730485
PHI_GOLDEN = (1 + math.sqrt(5)) / 2  # 1.61803399
MU = 1836.15267344  # m_p / m_e
ALPHA = 0.0072973525693


def check_exponential_relations():
    """Check exp(k * phi_vev) against mu for various k values.

    Searches for a simple k (integer, half-integer, third-integer)
    such that exp(k * phi_vev) = mu or a simple function of mu.

    Returns
    -------
    dict with keys:
        relations : list of dict -- each with k, exp_val, ratio_to_mu
        best_k : float -- k value giving closest ratio to 1
        best_ratio : float -- closest ratio to 1
        exact_k : float -- k that would give exp(k*phi_vev) = mu exactly
        clean_match : bool -- True if any ratio within 1% of 1
    """
    k_values = []
    for num in range(-8, 9):
        for den in [1, 2, 3]:
            if num == 0:
                continue
            k_values.append(num / den)
    k_values = sorted(set(k_values))

    relations = []
    best_k = None
    best_ratio = float("inf")

    for k in k_values:
        val = math.exp(k * PHI_VEV)
        ratio = val / MU
        relations.append({
            "k": k,
            "exp_val": val,
            "ratio_to_mu": ratio,
            "inverse_ratio": MU / val,
        })
        if abs(ratio - 1) < abs(best_ratio - 1):
            best_ratio = ratio
            best_k = k

    exact_k = math.log(MU) / PHI_VEV  # k such that exp(k*phi_vev) = mu

    return {
        "relations": relations,
        "best_k": best_k,
        "best_ratio": best_ratio,
        "exact_k": exact_k,
        "clean_match": abs(best_ratio - 1) < 0.01,
        "description": (
            f"Exact k for exp(k*phi_vev) = mu: k = {exact_k:.6f} "
            f"(not a simple number). Best simple k = {best_k} gives "
            f"ratio = {best_ratio:.6e}. No clean match found."
        ),
    }


def check_combination_relations():
    """Check combinations of alpha, phi_golden, and phi_vev against mu.

    Tests various algebraic combinations to see if mu can be expressed
    in terms of the framework quantities.

    Returns
    -------
    dict with keys:
        combinations : list of dict -- each with name, value, ratio_to_mu
        any_match : bool
        closest_name : str
        closest_ratio : float
    """
    combos = [
        ("mu * alpha", MU * ALPHA),
        ("mu * alpha^2", MU * ALPHA ** 2),
        ("mu / phi_golden", MU / PHI_GOLDEN),
        ("mu / phi_golden^2", MU / PHI_GOLDEN ** 2),
        ("sqrt(mu)", math.sqrt(MU)),
        ("mu^(1/3)", MU ** (1 / 3)),
        ("ln(mu)", math.log(MU)),
        ("ln(mu)/ln(1/alpha)", math.log(MU) / math.log(1 / ALPHA)),
        ("exp(-phi_vev)", math.exp(-PHI_VEV)),
        ("exp(-2*phi_vev)", math.exp(-2 * PHI_VEV)),
        ("exp(-phi_vev)/alpha", math.exp(-PHI_VEV) / ALPHA),
        ("1/(4*pi*alpha^2)", 1 / (4 * math.pi * ALPHA ** 2)),
        ("4*pi/alpha", 4 * math.pi / ALPHA),
        ("sqrt(4*pi/alpha)", math.sqrt(4 * math.pi / ALPHA)),
        ("mu * exp(phi_vev)", MU * math.exp(PHI_VEV)),
        ("mu^2 * 4*pi*alpha", MU ** 2 * 4 * math.pi * ALPHA),
    ]

    results = []
    best_name = None
    best_ratio = float("inf")

    for name, val in combos:
        ratio = val / MU
        results.append({
            "name": name,
            "value": val,
            "ratio_to_mu": ratio,
        })
        if 0 < abs(ratio - 1) < abs(best_ratio - 1):
            best_ratio = ratio
            best_name = name

    return {
        "combinations": results,
        "any_match": abs(best_ratio - 1) < 0.01 if best_name else False,
        "closest_name": best_name,
        "closest_ratio": best_ratio,
    }


def check_bridge_consistency():
    """Verify that the mu-structure and corrected bridge agree.

    The mu-structure formula:
        alpha_G = alpha^24 * mu * (mu - sqrt(phi_golden) * (1 - alpha))

    The corrected bridge:
        alpha_G = (phi_golden^2/2) * (1 + 3*alpha^2 + (phi/2)*alpha^3) * alpha^21

    These should agree to ~0.001 ppm.

    Returns
    -------
    dict with keys:
        bridge_lhs, bridge_rhs, match_ppm, consistent
    """
    F = 1 + 3 * ALPHA ** 2 + (PHI_GOLDEN / 2) * ALPHA ** 3
    lhs = (PHI_GOLDEN ** 2 / 2) * F
    rhs = ALPHA ** 3 * MU * (MU - math.sqrt(PHI_GOLDEN) * (1 - ALPHA))
    match_ppm = abs(lhs - rhs) / lhs * 1e6

    return {
        "bridge_lhs": lhs,
        "bridge_rhs": rhs,
        "match_ppm": match_ppm,
        "consistent": match_ppm < 0.01,
        "description": (
            f"Bridge LHS (phi^2/2 * F) = {lhs:.10f}, "
            f"RHS (alpha^3 * mu-structure) = {rhs:.10f}, "
            f"match = {match_ppm:.4f} ppm."
        ),
    }


def check_specific_relations():
    """Check specific interesting numerical relations.

    Returns
    -------
    dict with detailed checks
    """
    e_phi = math.exp(PHI_VEV)
    e_neg_phi = math.exp(-PHI_VEV)

    return {
        "mu_times_e_phi": MU * e_phi,
        "mu_over_e_neg_phi": MU / e_neg_phi,
        "mu_sq_times_4pi_alpha": MU ** 2 * 4 * math.pi * ALPHA,
        "mu_sq_times_e4phi": MU ** 2 * math.exp(4 * PHI_VEV),
        "four_pi_over_alpha": 4 * math.pi / ALPHA,
        "ratio_4pi_alpha_to_mu": (4 * math.pi / ALPHA) / MU,
        "description": (
            f"mu * e^(phi_vev) = {MU * e_phi:.4f} (not a clean number). "
            f"mu^2 * 4*pi*alpha = {MU ** 2 * 4 * math.pi * ALPHA:.4f} "
            f"(= mu^2 * e^(4*phi_vev), a tautology via matching). "
            f"4*pi/alpha = {4 * math.pi / ALPHA:.2f} vs mu = {MU:.2f} "
            f"(within 6%, closest hit but not a match)."
        ),
    }


def summarize_mu_vev_analysis():
    """Full summary of the mu vs phi_vev relationship search.

    Returns
    -------
    dict with all results and honest assessment
    """
    exp_check = check_exponential_relations()
    combo_check = check_combination_relations()
    bridge_check = check_bridge_consistency()
    specific_check = check_specific_relations()

    return {
        "exponential_relations": exp_check,
        "combination_relations": combo_check,
        "bridge_consistency": bridge_check,
        "specific_relations": specific_check,
        "phi_vev": PHI_VEV,
        "mu": MU,
        "alpha": ALPHA,
        "phi_golden": PHI_GOLDEN,
        "any_relation_found": False,
        "honest_assessment": (
            "The dilaton vev phi_vev = -0.597 and the proton-to-electron "
            "mass ratio mu = 1836.15 are numerically independent in the "
            "Alpha Ladder framework. phi_vev is determined by alpha_EM "
            "alone (exp(4*phi_vev) = 4*pi*alpha). mu enters through the "
            "mu-structure formula, which is algebraically equivalent to "
            "the corrected bridge but involves different physics (mass "
            "ratios vs gauge couplings). In the braneworld picture, "
            "fermion masses are brane-localized 4D parameters from the "
            "Higgs mechanism, independent of the S^2 geometry and dilaton. "
            "No simple exp(k*phi_vev) = f(mu) relation exists."
        ),
        "key_messages": [
            "phi_vev = gauge physics (alpha matching via KK reduction)",
            "mu = mass physics (Higgs/Yukawa sector, brane-localized)",
            "These are independent structures in the current framework",
            "The mu-structure formula's success is NOT explained by gauge matching",
            "Closest numerical hit: 4*pi/alpha = 1722 vs mu = 1836 (6% off, not a match)",
            "The bridge formula agreement (0.00 ppm) confirms mu-structure equivalence",
        ],
    }


if __name__ == "__main__":
    print("=" * 72)
    print("PROTON-TO-ELECTRON MASS RATIO vs DILATON VEV")
    print("=" * 72)
    print()

    print(f"phi_vev (dilaton)  = {PHI_VEV:.10f}")
    print(f"e^(phi_vev)        = {math.exp(PHI_VEV):.10f}")
    print(f"phi_golden         = {PHI_GOLDEN:.10f}")
    print(f"mu = m_p/m_e       = {MU:.10f}")
    print(f"alpha              = {ALPHA:.13f}")
    print()

    # Exponential relations
    print("-" * 72)
    print("1. EXPONENTIAL RELATIONS: exp(k * phi_vev) vs mu")
    print("-" * 72)
    exp_r = check_exponential_relations()
    print(f"{'k':>8s}  {'exp(k*phi_vev)':>16s}  {'ratio to mu':>14s}")
    print("-" * 44)
    for r in exp_r["relations"]:
        if abs(r["k"]) <= 4 and r["k"] == int(r["k"]):
            print(f"{r['k']:8.1f}  {r['exp_val']:16.8f}  {r['ratio_to_mu']:14.8e}")
    print()
    print(f"Exact k for match: {exp_r['exact_k']:.6f} (not simple)")
    print(f"Clean match found: {exp_r['clean_match']}")

    # Combination relations
    print()
    print("-" * 72)
    print("2. ALGEBRAIC COMBINATIONS vs mu")
    print("-" * 72)
    combo_r = check_combination_relations()
    print(f"{'Expression':>30s}  {'Value':>16s}  {'ratio to mu':>14s}")
    print("-" * 65)
    for c in combo_r["combinations"]:
        marker = ""
        if 0.9 < c["ratio_to_mu"] < 1.1:
            marker = " <-- CLOSE"
        print(f"{c['name']:>30s}  {c['value']:16.6f}  {c['ratio_to_mu']:14.8e}{marker}")

    # Bridge consistency
    print()
    print("-" * 72)
    print("3. BRIDGE FORMULA CONSISTENCY")
    print("-" * 72)
    bridge_r = check_bridge_consistency()
    print(bridge_r["description"])

    # Specific checks
    print()
    print("-" * 72)
    print("4. SPECIFIC NUMERICAL CHECKS")
    print("-" * 72)
    spec = check_specific_relations()
    print(f"mu * exp(phi_vev)      = {spec['mu_times_e_phi']:.4f}")
    print(f"mu^2 * 4*pi*alpha      = {spec['mu_sq_times_4pi_alpha']:.4f}")
    print(f"4*pi/alpha             = {spec['four_pi_over_alpha']:.2f}  vs  mu = {MU:.2f}")
    print(f"  ratio                = {spec['ratio_4pi_alpha_to_mu']:.4f}  (6% off)")

    # Conclusion
    print()
    print("=" * 72)
    print("CONCLUSION")
    print("=" * 72)
    summary = summarize_mu_vev_analysis()
    for msg in summary["key_messages"]:
        print(f"  - {msg}")
    print()
    print(summary["honest_assessment"])
    print()
    print("=" * 72)
