"""
Brute force search for the best rung spacing in the Standard Model mass spectrum.
Refactored from legacy/rung_spacing_search.py.
"""

import math

from alpha_ladder_core.constants import get_particle_masses


phi_f = (1 + math.sqrt(5)) / 2


def compute_rungs(constants: dict) -> dict:
    """
    Compute rung values (log ratios) for all particles.
    Returns dict mapping particle name to rung value.
    """
    alpha = float(constants.alpha) if hasattr(constants, 'alpha') else constants["alpha"]
    inv_alpha = 1.0 / alpha
    masses = get_particle_masses()
    m_e = masses["Electron"]

    rungs = {}
    for name, m in masses.items():
        if name == "Electron":
            continue
        rungs[name] = math.log(m / m_e) / math.log(inv_alpha)
    return rungs


def score_spacing(spacing: float, rung_dict: dict, tolerance_frac: float = 0.15) -> tuple:
    """
    Score how well particles land on rungs with given spacing.
    Returns (avg_closeness, n_matches, details).
    Closeness = distance to nearest multiple of spacing, as fraction of spacing.
    Exactly as legacy.
    """
    details = []
    for name, n in rung_dict.items():
        nearest_k = round(n / spacing)
        nearest_rung = nearest_k * spacing
        delta = abs(n - nearest_rung)
        frac = delta / spacing  # 0 = perfect, 0.5 = worst
        match = frac < tolerance_frac
        details.append((name, n, nearest_rung, nearest_k, delta, frac, match))

    n_matches = sum(1 for d in details if d[6])
    avg_frac = sum(d[5] for d in details) / len(details)
    return avg_frac, n_matches, details


def search_rational_spacings(rung_dict: dict, k_max: int = 24) -> list:
    """
    Test rational spacings 1/k for k=1..k_max.
    Returns list of (k, spacing, avg, matches, details) sorted by matches descending.
    """
    results = []
    for k in range(1, k_max + 1):
        spacing = 1.0 / k
        avg, matches, details = score_spacing(spacing, rung_dict)
        results.append((k, spacing, avg, matches, details))
    results.sort(key=lambda x: -x[3])
    return results


def search_irrational_spacings(rung_dict: dict) -> list:
    """
    Test all 15 irrational spacings from legacy.
    Returns list of (name, spacing, avg, matches, details) sorted by matches descending.
    """
    irrational_spacings = {
        "1/φ":              1.0 / phi_f,
        "1/π":              1.0 / math.pi,
        "1/e":              1.0 / math.e,
        "ln2":              math.log(2),
        "1/φ²":             1.0 / phi_f**2,
        "φ/π":              phi_f / math.pi,
        "2/π":              2.0 / math.pi,
        "1/√2":             1.0 / math.sqrt(2),
        "1/√3":             1.0 / math.sqrt(3),
        "φ − 1 = 1/φ":     phi_f - 1,
        "ln(φ)":            math.log(phi_f),
        "π/6":              math.pi / 6,
        "1/(2φ)":           1.0 / (2 * phi_f),
        "1/(3φ)":           1.0 / (3 * phi_f),
        "φ/4":              phi_f / 4,
    }

    results = []
    for name, spacing in irrational_spacings.items():
        avg, matches, details = score_spacing(spacing, rung_dict)
        results.append((name, spacing, avg, matches, details))
    results.sort(key=lambda x: -x[3])
    return results


def search_continuous_optimum(rung_values: list) -> tuple:
    """
    Brute force continuous optimization: scan spacings from 0.01 to 2.0.
    Returns (best_spacing, best_score, local_mins) exactly as legacy.
    local_mins is sorted by RMS ascending.
    """
    best_score = float('inf')
    best_spacing = 0
    scores = []

    # Scan from 0.01 to 2.0 in steps of 0.001
    for i in range(10, 2001):
        s = i / 1000.0
        total = 0
        for n in rung_values:
            nearest = round(n / s) * s
            total += (n - nearest) ** 2
        rms = math.sqrt(total / len(rung_values))
        scores.append((s, rms))
        if rms < best_score:
            best_score = rms
            best_spacing = s

    # Find local minima
    local_mins = []
    for i in range(1, len(scores) - 1):
        if scores[i][1] < scores[i-1][1] and scores[i][1] < scores[i+1][1]:
            local_mins.append(scores[i])

    local_mins.sort(key=lambda x: x[1])
    return best_spacing, best_score, local_mins


def get_best_fit_details(spacing: float, rung_dict: dict) -> list:
    """
    Return details for a given spacing applied to all particles.
    Returns details list sorted by rung value, each entry is:
      (name, n, nearest_rung, k, delta, frac, match)
    """
    _, _, details = score_spacing(spacing, rung_dict)
    details.sort(key=lambda x: x[1])
    return details
