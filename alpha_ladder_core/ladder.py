"""
Core Alpha Ladder computations, refactored from legacy/alpha_ladder.py.

Every function accepts a *constants* namespace (produced by
``get_constants()``) so that the caller controls which CODATA edition
is used.
"""

import math
from decimal import Decimal

from alpha_ladder_core.constants import get_particle_masses


def compute_electron_geometry(constants):
    """Return a dict of electron-radius estimates.

    Keys
    ----
    r_e_from_compton : Decimal
        alpha * lambda_bar_c
    r_e_from_bohr : Decimal
        alpha^2 * a_0
    r_e_nist : Decimal
        The NIST reference value carried in the constants namespace.
    """
    return {
        "r_e_from_compton": constants.alpha * constants.lambda_bar_c,
        "r_e_from_bohr": constants.alpha ** 2 * constants.a_0,
        "r_e_nist": constants.r_e_nist,
    }


def calculate_geometric_rungs(constants):
    """Build the geometric ladder of alpha powers.

    Returns
    -------
    tuple
        (rungs_list, bridges_dict, alpha_dm, pi_bridge)

        *rungs_list* -- list of dicts with keys ``power``, ``value``,
        ``ratio_to_gravity``, ``label`` for n = 1 .. 24.

        *bridges_dict* -- candidate expressions for alpha_G = C * alpha^21.

        *alpha_dm* -- Decimal, the dark-sector candidate alpha^10.

        *pi_bridge* -- Decimal, alpha^21 * pi^2.
    """
    alpha = constants.alpha
    alpha_g = constants.alpha_g
    pi = constants.pi
    phi = constants.phi
    e = constants.e

    rungs = []
    for n in range(1, 25):
        val = alpha ** n
        ratio = val / alpha_g
        label = ""
        if n == 10:
            label = "Dark Matter Candidate (The Blank)"
        elif n == 21:
            label = "Gravity Floor (The Bridge)"
        rungs.append({
            "power": n,
            "value": val,
            "ratio_to_gravity": ratio,
            "label": label,
        })

    alpha_21 = alpha ** 21

    # Bridge candidates for alpha_G = C * alpha^21
    bridges = {
        "φ · α²¹  (original)":      phi * alpha_21,
        "(φ²/2) · α²¹":          (phi ** 2 / 2) * alpha_21,
        "(5/12) · π · α²¹":          (Decimal(5) / Decimal(12)) * pi * alpha_21,
        "√e / ³√2 · α²¹":    (e.sqrt() / Decimal(2) ** (Decimal(1) / Decimal(3))) * alpha_21,
    }

    # Dark sector candidate: alpha_DM ~ alpha^10
    alpha_dm = alpha ** 10

    # Pi-Bridge
    pi_bridge = alpha_21 * (pi ** 2)

    return rungs, bridges, alpha_dm, pi_bridge


def compute_gap(constants):
    """Return alpha / alpha_g  (the ~42-order-of-magnitude gap)."""
    return constants.alpha / constants.alpha_g


# ---------------------------------------------------------------------------
# Physical rung map
# ---------------------------------------------------------------------------

def _build_sm_particles():
    """Build SM particle list in eV from the shared particle mass source."""
    masses = get_particle_masses()
    return [(name, m * 1e6) for name, m in masses.items()]


_SPECIAL_SCALES = [
    ("QCD Lambda (~200 MeV)", 200e6),
    ("Weak scale v (246 GeV)", 246e9),
    ("GUT scale (~1e16 GeV)", 1e25),
    ("Planck mass", 1.22089e28),
    ("Neutrino mass (~0.05 eV)", 0.05),
    ("CMB temperature (2.725 K)", 2.35e-4),
    ("Hubble scale (H0)", 1.49e-33),
]


def compute_physical_rung_map(constants, k_max=30):
    """Map each rung k=0..k_max to physical length/energy scales.

    The base identity r_e = alpha * lambda_bar_c = alpha^2 * a_0 defines
    the first two rungs. Each subsequent rung multiplies by alpha (~1/137)
    in length, or by 1/alpha (~137) in energy.

    Parameters
    ----------
    constants : SimpleNamespace
        Physical constants from get_constants().
    k_max : int
        Maximum rung index (default 30).

    Returns
    -------
    dict with keys:
        rungs : list of dict -- one per k=0..k_max, each with:
            k, alpha_k, length_m, energy_eV,
            nearest_particle (dict with name, energy_eV, rung_position, offset_dex),
            nearest_scale (dict with name, energy_eV, rung_position, offset_dex)
        particle_rungs : list of dict -- each SM particle with:
            name, mass_eV, rung_position (float), nearest_half (float), residual
        scale_rungs : list of dict -- each special scale with:
            name, energy_eV, rung_position (float)
        key_matches : list of dict -- all matches within 0.5 dex
        k_planck : float -- rung position of Planck length
        k_gravity : float -- effective rung of alpha_g
        base_identity : dict -- r_e = alpha * lambda_bar_c = alpha^2 * a_0 verification
    """
    alpha_f = float(constants.alpha)
    hbar_f = float(constants.hbar)
    c_f = float(constants.c)
    m_e_f = float(constants.m_e)
    e_charge_f = float(constants.e_charge)
    G_f = float(constants.G)

    sm_particles = _build_sm_particles()

    log_inv_alpha = math.log(1.0 / alpha_f)

    # Bohr radius: a_0 = hbar / (m_e * c * alpha)
    a_0_float = hbar_f / (m_e_f * c_f * alpha_f)

    # Electron mass in eV
    m_e_eV = m_e_f * c_f ** 2 / e_charge_f

    # ------------------------------------------------------------------
    # Helper: find nearest entry from a list of (name, energy_eV)
    # ------------------------------------------------------------------
    def _nearest(energy_eV, catalog):
        best_name = None
        best_eV = None
        best_offset = float("inf")
        best_rp = 0.0
        for name, cat_eV in catalog:
            if cat_eV <= 0:
                continue
            offset_dex = abs(math.log10(energy_eV / cat_eV))
            if offset_dex < best_offset:
                best_offset = offset_dex
                best_name = name
                best_eV = cat_eV
                best_rp = math.log(cat_eV / m_e_eV) / log_inv_alpha
        return {
            "name": best_name,
            "energy_eV": best_eV,
            "rung_position": best_rp,
            "offset_dex": best_offset,
        }

    # ------------------------------------------------------------------
    # Build rungs k = 0 .. k_max
    # ------------------------------------------------------------------
    rungs = []
    for k in range(k_max + 1):
        alpha_k = alpha_f ** k
        length_m = alpha_k * a_0_float
        energy_eV = hbar_f * c_f / (length_m * e_charge_f)

        nearest_particle = _nearest(energy_eV, sm_particles)
        nearest_scale = _nearest(energy_eV, _SPECIAL_SCALES)

        rungs.append({
            "k": k,
            "alpha_k": alpha_k,
            "length_m": length_m,
            "energy_eV": energy_eV,
            "nearest_particle": nearest_particle,
            "nearest_scale": nearest_scale,
        })

    # ------------------------------------------------------------------
    # Particle rungs
    # ------------------------------------------------------------------
    particle_rungs = []
    for name, mass_eV in sm_particles:
        rung_position = math.log(mass_eV / m_e_eV) / log_inv_alpha
        nearest_half = round(rung_position * 2) / 2
        residual = rung_position - nearest_half
        particle_rungs.append({
            "name": name,
            "mass_eV": mass_eV,
            "rung_position": rung_position,
            "nearest_half": nearest_half,
            "residual": residual,
        })
    particle_rungs.sort(key=lambda d: d["rung_position"])

    # ------------------------------------------------------------------
    # Scale rungs
    # ------------------------------------------------------------------
    scale_rungs = []
    for name, energy_eV in _SPECIAL_SCALES:
        if energy_eV > 0:
            rung_position = math.log(energy_eV / m_e_eV) / log_inv_alpha
        else:
            rung_position = float("nan")
        scale_rungs.append({
            "name": name,
            "energy_eV": energy_eV,
            "rung_position": rung_position,
        })

    # ------------------------------------------------------------------
    # Key matches: rung k whose energy is within 0.5 dex of a known scale
    # ------------------------------------------------------------------
    all_scales = sm_particles + _SPECIAL_SCALES
    key_matches = []
    for rung in rungs:
        k = rung["k"]
        e_k = rung["energy_eV"]
        for name, scale_eV in all_scales:
            if scale_eV <= 0:
                continue
            dex = abs(math.log10(e_k / scale_eV))
            if dex < 0.5:
                key_matches.append({
                    "k": k,
                    "scale_name": name,
                    "rung_energy_eV": e_k,
                    "scale_energy_eV": scale_eV,
                    "offset_dex": dex,
                })
    key_matches.sort(key=lambda d: d["k"])

    # ------------------------------------------------------------------
    # Planck rung position
    # ------------------------------------------------------------------
    l_Pl = math.sqrt(hbar_f * G_f / c_f ** 3)
    k_planck = math.log(l_Pl / a_0_float) / math.log(alpha_f)

    # ------------------------------------------------------------------
    # Gravity rung position
    # ------------------------------------------------------------------
    alpha_g = G_f * m_e_f ** 2 / (hbar_f * c_f)
    k_gravity = math.log(alpha_g) / math.log(alpha_f)

    # ------------------------------------------------------------------
    # Base identity verification
    # ------------------------------------------------------------------
    lambda_bar_c = hbar_f / (m_e_f * c_f)
    r_e_from_compton = alpha_f * lambda_bar_c
    r_e_from_bohr = alpha_f ** 2 * a_0_float
    base_identity = {
        "r_e_from_compton": r_e_from_compton,
        "r_e_from_bohr": r_e_from_bohr,
        "a_0": a_0_float,
        "lambda_bar_c": lambda_bar_c,
        "relative_error": abs(r_e_from_compton - r_e_from_bohr) / r_e_from_compton,
    }

    return {
        "rungs": rungs,
        "particle_rungs": particle_rungs,
        "scale_rungs": scale_rungs,
        "key_matches": key_matches,
        "k_planck": k_planck,
        "k_gravity": k_gravity,
        "base_identity": base_identity,
    }
