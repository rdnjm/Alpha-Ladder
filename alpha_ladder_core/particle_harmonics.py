"""
Scan particle masses for alpha-power harmonics.
Refactored from legacy/particle_harmonics.py.
"""

import math

from alpha_ladder_core.constants import get_particle_masses


def compute_harmonics(constants: dict) -> list:
    """
    Compute alpha-power rung for each particle mass relative to the electron.
    Returns list of dicts with keys:
      name, mass, ratio, rung_n, nearest_half, closeness, is_match
    Uses alpha from constants dict.
    """
    alpha = float(constants.alpha) if hasattr(constants, 'alpha') else constants["alpha"]
    inv_alpha = 1.0 / alpha
    masses = get_particle_masses()
    m_electron = masses["Electron"]

    results = []
    for name, m in masses.items():
        ratio = m / m_electron
        # n = log(ratio) / log(1/alpha)
        n = math.log(ratio) / math.log(inv_alpha)
        nearest_half = round(n * 2) / 2
        closeness = abs(n - nearest_half)
        is_match = closeness < 0.05

        results.append({
            "name": name,
            "mass": m,
            "ratio": ratio,
            "rung_n": n,
            "nearest_half": nearest_half,
            "closeness": closeness,
            "is_match": is_match,
        })

    return results
