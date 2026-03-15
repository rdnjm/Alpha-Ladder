"""
Radius fixing mechanisms for the S^2 compactification.

The radius_determination module proved that the 6D EH+GB action has an exact
scaling symmetry for n=2 internal dimensions, leaving the internal radius a_0
as a flat direction.  The flux_stabilization module showed that adding quantized
2-form flux stabilizes the volume modulus sigma but V_min(a_0) is monotonic --
no finite a_0 minimizes the vacuum energy.

This module implements THREE mechanisms that attempt to fix a_0 by introducing
new a_0-dependent contributions to the effective potential:

    1. Coleman-Weinberg potential from anomaly-constrained matter
       (1-loop log(a_0) terms from KK towers of hypers and vectors)

    2. Warped S^2 compactification
       (gradient energy from warp factor A(theta) = epsilon*cos(theta))

    3. Orbifold S^2/Z_2
       (modified Casimir from even-l modes + brane tension at fixed points)

Honest expectations:
    - Mechanism 1 (CW): the log(a_0) factor is qualitatively new but the
      overall CW scale is suppressed by 1/(64*pi^2); may or may not compete
      with the tree-level flux potential.
    - Mechanism 2 (warp): gradient energy scales as a_0^{-2}, still monotonic;
      unlikely to produce a minimum.
    - Mechanism 3 (orbifold): brane tension ~ a_0^{-2} competes with Casimir
      ~ a_0^{-4} at different scales; BEST candidate for creating a minimum.

All calculations use pure Python math (no numpy/scipy).
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50

# ---------------------------------------------------------------------------
# Graceful imports
# ---------------------------------------------------------------------------

try:
    from alpha_ladder_core.casimir_stabilization import compute_casimir_energy_s2
    _CASIMIR_AVAILABLE = True
except ImportError:
    _CASIMIR_AVAILABLE = False

try:
    from alpha_ladder_core.flux_stabilization import find_flux_minimum
    _FLUX_AVAILABLE = True
except ImportError:
    _FLUX_AVAILABLE = False

try:
    from alpha_ladder_core.anomaly_cancellation import _ANOMALY_FREE_GROUPS
    _ANOMALY_AVAILABLE = True
except ImportError:
    _ANOMALY_AVAILABLE = False
    _ANOMALY_FREE_GROUPS = {}


# ===========================================================================
# MECHANISM 1: Coleman-Weinberg from anomaly-constrained matter
# ===========================================================================

# ---------------------------------------------------------------------------
# 1a. CW sum for real scalar on S^2
# ---------------------------------------------------------------------------

def _cw_sum_scalar(a_0, l_max=100, mu=1.0):
    """
    Coleman-Weinberg 1-loop sum for a real scalar field on S^2.

    The KK spectrum of a real scalar on S^2 of radius a_0 is:

        m_l^2 = l(l+1) / a_0^2,    l = 1, 2, ..., l_max
        degeneracy = 2l + 1

    The CW contribution (before the 1/(64*pi^2) prefactor) is:

        Sigma = sum_{l=1}^{l_max} (2l+1) * m_l^2 * log(m_l^2 / mu^2)

    The sum diverges as l_max -> infinity.  We regularize using
    Euler-Maclaurin subtraction: subtract the leading divergent piece
    (the l^3 * log(l) asymptotic behavior) to obtain a finite result.

    Parameters
    ----------
    a_0 : float
        Radius of the S^2 in Planck units.
    l_max : int
        Truncation of the KK tower (default 100).
    mu : float
        Renormalization scale in Planck units (default 1.0).

    Returns
    -------
    dict with keys:
        raw_sum : float -- unregularized sum
        subtracted_divergence : float -- leading divergent piece
        regularized_sum : float -- raw_sum - subtracted_divergence
        a_0 : float
        l_max : int
        mu : float
        field_type : str -- "real_scalar"
    """
    raw_sum = 0.0
    divergent_piece = 0.0
    mu_sq = mu * mu

    for l in range(1, l_max + 1):
        m_sq = l * (l + 1) / (a_0 * a_0)
        deg = 2 * l + 1

        if m_sq > 0 and mu_sq > 0:
            raw_sum += deg * m_sq * math.log(m_sq / mu_sq)

        # Leading asymptotic: (2l+1)*l(l+1)/a_0^2 ~ 2*l^3/a_0^2 for large l
        # log(l(l+1)/a_0^2 / mu^2) ~ 2*log(l) + log(1/a_0^2/mu^2)
        # Divergent piece ~ 2*l^3/a_0^2 * (2*log(l) + log(1/(a_0^2*mu^2)))
        if l >= 2:
            asymp_m_sq = l * l / (a_0 * a_0)
            asymp_deg = 2.0 * l
            log_piece = 2.0 * math.log(l) + math.log(1.0 / (a_0 * a_0 * mu_sq)) if (a_0 * a_0 * mu_sq) > 0 else 0.0
            divergent_piece += asymp_deg * asymp_m_sq * log_piece

    regularized = raw_sum - divergent_piece

    return {
        "raw_sum": raw_sum,
        "subtracted_divergence": divergent_piece,
        "regularized_sum": regularized,
        "a_0": a_0,
        "l_max": l_max,
        "mu": mu,
        "field_type": "real_scalar",
    }


# ---------------------------------------------------------------------------
# 1b. CW sum for Weyl fermion on S^2
# ---------------------------------------------------------------------------

def _cw_sum_fermion(a_0, l_max=100, mu=1.0):
    """
    Coleman-Weinberg 1-loop sum for a Weyl fermion on S^2.

    The KK spectrum of a Weyl fermion on S^2 of radius a_0 is:

        m_l^2 = (l + 1/2)^2 / a_0^2,    l = 0, 1, ..., l_max
        degeneracy = 2(2l + 1)

    The extra factor of 2 relative to scalars comes from the two spin
    components of the Weyl fermion on S^2.  The half-integer shift arises
    from the spin connection on S^2.

    The CW contribution (before the 1/(64*pi^2) prefactor, with fermion
    sign -1) is:

        Sigma = -sum_{l=0}^{l_max} 2(2l+1) * m_l^2 * log(m_l^2 / mu^2)

    The minus sign is the standard fermion contribution to the CW potential.

    Parameters
    ----------
    a_0 : float
        Radius of the S^2 in Planck units.
    l_max : int
        Truncation of the KK tower (default 100).
    mu : float
        Renormalization scale in Planck units (default 1.0).

    Returns
    -------
    dict with keys:
        raw_sum : float -- unregularized sum (includes fermion sign)
        subtracted_divergence : float -- leading divergent piece
        regularized_sum : float -- raw_sum - subtracted_divergence
        a_0 : float
        l_max : int
        mu : float
        field_type : str -- "weyl_fermion"
    """
    raw_sum = 0.0
    divergent_piece = 0.0
    mu_sq = mu * mu

    for l in range(0, l_max + 1):
        m_sq = (l + 0.5) ** 2 / (a_0 * a_0)
        deg = 2 * (2 * l + 1)

        if m_sq > 0 and mu_sq > 0:
            # Fermion sign: negative contribution
            raw_sum += -deg * m_sq * math.log(m_sq / mu_sq)

        # Leading asymptotic subtraction (same structure, with fermion sign)
        if l >= 2:
            asymp_m_sq = l * l / (a_0 * a_0)
            asymp_deg = 4.0 * l
            log_piece = 2.0 * math.log(l) + math.log(1.0 / (a_0 * a_0 * mu_sq)) if (a_0 * a_0 * mu_sq) > 0 else 0.0
            divergent_piece += -asymp_deg * asymp_m_sq * log_piece

    regularized = raw_sum - divergent_piece

    return {
        "raw_sum": raw_sum,
        "subtracted_divergence": divergent_piece,
        "regularized_sum": regularized,
        "a_0": a_0,
        "l_max": l_max,
        "mu": mu,
        "field_type": "weyl_fermion",
    }


# ---------------------------------------------------------------------------
# 1c. CW sum for gauge vector on S^2
# ---------------------------------------------------------------------------

def _cw_sum_vector(a_0, l_max=100, mu=1.0):
    """
    Coleman-Weinberg 1-loop sum for a gauge vector boson on S^2.

    The KK spectrum of a gauge vector on S^2 of radius a_0 is:

        m_l^2 = l(l+1) / a_0^2,    l = 1, 2, ..., l_max
        degeneracy = 2(2l + 1)

    The factor of 2 comes from the two transverse polarizations of the
    gauge field in 4D.  The l=0 mode is excluded (pure gauge for abelian,
    or the 4D gauge zero mode for non-abelian).

    The CW contribution (before the 1/(64*pi^2) prefactor) is:

        Sigma = sum_{l=1}^{l_max} 2(2l+1) * m_l^2 * log(m_l^2 / mu^2)

    Parameters
    ----------
    a_0 : float
        Radius of the S^2 in Planck units.
    l_max : int
        Truncation of the KK tower (default 100).
    mu : float
        Renormalization scale in Planck units (default 1.0).

    Returns
    -------
    dict with keys:
        raw_sum : float -- unregularized sum
        subtracted_divergence : float -- leading divergent piece
        regularized_sum : float -- raw_sum - subtracted_divergence
        a_0 : float
        l_max : int
        mu : float
        field_type : str -- "gauge_vector"
    """
    raw_sum = 0.0
    divergent_piece = 0.0
    mu_sq = mu * mu

    for l in range(1, l_max + 1):
        m_sq = l * (l + 1) / (a_0 * a_0)
        deg = 2 * (2 * l + 1)

        if m_sq > 0 and mu_sq > 0:
            raw_sum += deg * m_sq * math.log(m_sq / mu_sq)

        # Leading asymptotic subtraction
        if l >= 2:
            asymp_m_sq = l * l / (a_0 * a_0)
            asymp_deg = 4.0 * l
            log_piece = 2.0 * math.log(l) + math.log(1.0 / (a_0 * a_0 * mu_sq)) if (a_0 * a_0 * mu_sq) > 0 else 0.0
            divergent_piece += asymp_deg * asymp_m_sq * log_piece

    regularized = raw_sum - divergent_piece

    return {
        "raw_sum": raw_sum,
        "subtracted_divergence": divergent_piece,
        "regularized_sum": regularized,
        "a_0": a_0,
        "l_max": l_max,
        "mu": mu,
        "field_type": "gauge_vector",
    }


# ---------------------------------------------------------------------------
# 1d. Total CW potential for anomaly-constrained matter
# ---------------------------------------------------------------------------

def compute_cw_potential(a_0, group_key="E8_x_E8", l_max=100, mu=1.0):
    """
    Compute the total Coleman-Weinberg potential at a given a_0 for matter
    content determined by anomaly cancellation.

    The matter content for each anomaly-free group is:
        - Each hypermultiplet contributes: 4 real scalars + 2 Weyl fermions
        - Each vector multiplet contributes: 1 gauge vector + 1 Weyl gaugino

    The total CW potential is:

        V_CW = (1 / (64 * pi^2)) * (
            n_hyper * (4 * Sigma_scalar + 2 * Sigma_fermion)
          + n_vector * (Sigma_vector + Sigma_fermion)
        )

    Parameters
    ----------
    a_0 : float
        Radius of the S^2 in Planck units.
    group_key : str
        Key in _ANOMALY_FREE_GROUPS (default "E8_x_E8").
    l_max : int
        Truncation of the KK tower (default 100).
    mu : float
        Renormalization scale in Planck units (default 1.0).

    Returns
    -------
    dict with keys:
        V_cw : float -- total CW potential
        scalar_sum : float -- regularized scalar sum
        fermion_sum : float -- regularized fermion sum
        vector_sum : float -- regularized vector sum
        n_hyper : int or None
        n_vector : int or None
        n_scalar_total : int -- 4 * n_hyper
        n_fermion_total : int -- 2 * n_hyper + n_vector (gauginos)
        n_vector_total : int -- n_vector
        group_key : str
        group_name : str
        a_0 : float
        l_max : int
        mu : float
        prefactor : float -- 1 / (64 * pi^2)
        description : str
    """
    if not _ANOMALY_AVAILABLE or group_key not in _ANOMALY_FREE_GROUPS:
        return {
            "V_cw": 0.0,
            "scalar_sum": 0.0,
            "fermion_sum": 0.0,
            "vector_sum": 0.0,
            "n_hyper": None,
            "n_vector": None,
            "n_scalar_total": 0,
            "n_fermion_total": 0,
            "n_vector_total": 0,
            "group_key": group_key,
            "group_name": "unknown",
            "a_0": a_0,
            "l_max": l_max,
            "mu": mu,
            "prefactor": 1.0 / (64.0 * math.pi ** 2),
            "description": (
                f"Group '{group_key}' not found in anomaly-free database. "
                f"CW potential set to zero."
            ),
        }

    group_data = _ANOMALY_FREE_GROUPS[group_key]
    n_hyper = group_data.get("n_hyper")
    n_vector = group_data.get("n_vector")
    group_name = group_data.get("group", group_key)

    # Handle None values (some groups lack n_hyper data)
    if n_hyper is None:
        n_hyper = 0
    if n_vector is None:
        n_vector = 0

    # Compute individual KK sums
    scalar_result = _cw_sum_scalar(a_0, l_max=l_max, mu=mu)
    fermion_result = _cw_sum_fermion(a_0, l_max=l_max, mu=mu)
    vector_result = _cw_sum_vector(a_0, l_max=l_max, mu=mu)

    scalar_sum = scalar_result["regularized_sum"]
    fermion_sum = fermion_result["regularized_sum"]
    vector_sum = vector_result["regularized_sum"]

    # Multiplicities
    # Hypers: 4 real scalars + 2 Weyl fermions per hyper
    # Vectors: 1 gauge vector + 1 Weyl gaugino per vector multiplet
    n_scalar_total = 4 * n_hyper
    n_fermion_total = 2 * n_hyper + n_vector  # hyper fermions + gauginos
    n_vector_total = n_vector

    prefactor = 1.0 / (64.0 * math.pi ** 2)

    V_cw = prefactor * (
        n_hyper * (4 * scalar_sum + 2 * fermion_sum)
        + n_vector * (vector_sum + fermion_sum)
    )

    description = (
        f"Coleman-Weinberg potential for {group_name} matter on S^2 "
        f"(a_0 = {a_0:.4f}): {n_hyper} hypers ({n_scalar_total} scalars, "
        f"{2 * n_hyper} fermions) + {n_vector} vectors ({n_vector} gauge "
        f"bosons, {n_vector} gauginos). V_CW = {V_cw:.6e}."
    )

    return {
        "V_cw": V_cw,
        "scalar_sum": scalar_sum,
        "fermion_sum": fermion_sum,
        "vector_sum": vector_sum,
        "n_hyper": n_hyper,
        "n_vector": n_vector,
        "n_scalar_total": n_scalar_total,
        "n_fermion_total": n_fermion_total,
        "n_vector_total": n_vector_total,
        "group_key": group_key,
        "group_name": group_name,
        "a_0": a_0,
        "l_max": l_max,
        "mu": mu,
        "prefactor": prefactor,
        "description": description,
    }


# ---------------------------------------------------------------------------
# 1e. CW radius scan
# ---------------------------------------------------------------------------

def compute_cw_radius_scan(a_0_values=None, group_key="E8_x_E8", N_flux=1,
                           l_max=100, constants=None):
    """
    Scan a_0 values and compute V_flux(a_0) + V_CW(a_0) to search for a
    minimum that fixes the internal radius.

    The key physics: V_flux(a_0) is a monotonic power law (~a_0^{-4}),
    while V_CW(a_0) contains log(a_0) factors from the KK sums.  The
    logarithm is qualitatively different and could create a turning point.

    For each a_0, we:
        1. Compute V_flux_min(a_0) from find_flux_minimum(N=N_flux, a_0=a_0)
        2. Compute V_CW(a_0) from the anomaly-constrained matter content
        3. Form V_total = V_flux_min + V_CW and look for a minimum

    Parameters
    ----------
    a_0_values : list of float or None
        Values of a_0 to scan.  If None, uses 20 log-spaced points
        from 0.1 to 100 Planck units.
    group_key : str
        Key in _ANOMALY_FREE_GROUPS (default "E8_x_E8").
    N_flux : int
        Flux quantum number (default 1).
    l_max : int
        Truncation of the KK tower (default 100).
    constants : ignored
        Planck units used throughout.

    Returns
    -------
    dict with keys:
        scan_data : list of dict -- a_0, V_flux, V_cw, V_total for each point
        minimum_found : bool
        minimum_a_0 : float or None
        minimum_V_total : float or None
        group_key : str
        N_flux : int
        l_max : int
        n_points : int
        description : str
    """
    if a_0_values is None:
        # 20 log-spaced points from 0.1 to 100
        log_min = math.log10(0.1)
        log_max = math.log10(100.0)
        n_pts = 20
        step = (log_max - log_min) / (n_pts - 1) if n_pts > 1 else 0
        a_0_values = [10.0 ** (log_min + i * step) for i in range(n_pts)]

    scan_data = []

    for a_0 in a_0_values:
        # Flux contribution
        v_flux = None
        if _FLUX_AVAILABLE:
            flux_result = find_flux_minimum(N=N_flux, a_0=a_0, constants=constants)
            if flux_result["minimum_exists"]:
                v_flux = flux_result["V_at_minimum"]

        # CW contribution
        cw_result = compute_cw_potential(a_0, group_key=group_key,
                                         l_max=l_max, mu=1.0)
        v_cw = cw_result["V_cw"]

        # Total
        v_total = None
        if v_flux is not None:
            v_total = v_flux + v_cw

        scan_data.append({
            "a_0": a_0,
            "V_flux": v_flux,
            "V_cw": v_cw,
            "V_total": v_total,
        })

    # Search for minimum in V_total
    minimum_found = False
    minimum_a_0 = None
    minimum_V_total = None

    valid_points = [(d["a_0"], d["V_total"]) for d in scan_data
                    if d["V_total"] is not None]

    if len(valid_points) >= 3:
        for i in range(1, len(valid_points) - 1):
            a_prev, v_prev = valid_points[i - 1]
            a_curr, v_curr = valid_points[i]
            a_next, v_next = valid_points[i + 1]

            # Local minimum: V decreasing then increasing
            if v_curr < v_prev and v_curr < v_next:
                minimum_found = True
                minimum_a_0 = a_curr
                minimum_V_total = v_curr
                break

    if minimum_found:
        description = (
            f"A local minimum in V_total(a_0) was found at a_0 = {minimum_a_0:.4f} "
            f"(Planck units) with V_total = {minimum_V_total:.6e}.  The Coleman-"
            f"Weinberg log(a_0) contribution from {group_key} matter creates a "
            f"turning point in the otherwise monotonic flux potential."
        )
    else:
        description = (
            f"No local minimum found in V_total(a_0) = V_flux + V_CW for "
            f"{group_key} matter with N_flux = {N_flux}.  The CW contribution, "
            f"while containing log(a_0) factors, is suppressed by the 1/(64*pi^2) "
            f"loop factor and does not overcome the monotonic power-law behavior "
            f"of the flux potential.  The internal radius remains unfixed by this "
            f"mechanism."
        )

    return {
        "scan_data": scan_data,
        "minimum_found": minimum_found,
        "minimum_a_0": minimum_a_0,
        "minimum_V_total": minimum_V_total,
        "group_key": group_key,
        "N_flux": N_flux,
        "l_max": l_max,
        "n_points": len(a_0_values),
        "description": description,
    }


# ---------------------------------------------------------------------------
# 1f. CW mechanism summary
# ---------------------------------------------------------------------------

def summarize_cw_mechanism(constants=None):
    """
    Run the Coleman-Weinberg radius-fixing mechanism and return a full summary.

    This computes the CW potential from anomaly-constrained matter (E8 x E8
    by default) combined with the flux potential, scanning over a_0 to
    search for a minimum.

    Parameters
    ----------
    constants : ignored
        Planck units used throughout.

    Returns
    -------
    dict with keys:
        mechanism_name : str
        scan_result : dict from compute_cw_radius_scan
        minimum_found : bool
        minimum_a_0 : float or None
        cw_potential_at_a0_1 : dict from compute_cw_potential at a_0=1
        group_used : str
        honest_assessment : str
        physics_summary : str
    """
    # Single-point evaluation at a_0 = 1 for diagnostics
    cw_at_1 = compute_cw_potential(1.0, group_key="E8_x_E8", l_max=100, mu=1.0)

    # Full scan
    scan = compute_cw_radius_scan(group_key="E8_x_E8", N_flux=1,
                                  l_max=100, constants=constants)

    minimum_found = scan["minimum_found"]
    minimum_a_0 = scan["minimum_a_0"]

    if minimum_found:
        honest_assessment = (
            f"The Coleman-Weinberg mechanism from E8 x E8 matter content "
            f"produces a minimum at a_0 = {minimum_a_0:.4f} Planck units.  "
            f"However, this result depends on the specific gauge group choice "
            f"and the regularization scheme for the KK sums.  The log(a_0) "
            f"dependence is genuine but the numerical coefficient is scheme-"
            f"dependent.  This should be treated as suggestive rather than "
            f"definitive."
        )
    else:
        honest_assessment = (
            "The Coleman-Weinberg mechanism from anomaly-constrained matter "
            "does NOT fix a_0.  While the log(a_0) factors from 1-loop KK "
            "sums are qualitatively new compared to the tree-level power laws, "
            "the overall CW contribution is suppressed by 1/(64*pi^2) ~ 1.6e-3.  "
            "This loop suppression factor means the CW potential is subdominant "
            "to the tree-level flux potential at all a_0 values.  The monotonic "
            "behavior of V_flux(a_0) dominates, and the internal radius remains "
            "unfixed."
        )

    physics_summary = (
        "The Coleman-Weinberg potential arises from 1-loop corrections of "
        "KK modes on S^2.  Each mode contributes m^2 * log(m^2/mu^2), where "
        "m^2 ~ l(l+1)/a_0^2.  The crucial feature is the log(a_0) dependence, "
        "which is absent from tree-level terms.  For E8 x E8: 740 hypers "
        "(2960 scalars, 1480 fermions) + 496 vectors (496 gauge bosons, "
        "496 gauginos).  Despite the large multiplicities, the 1/(64*pi^2) "
        "prefactor suppresses the CW potential relative to the flux potential."
    )

    return {
        "mechanism_name": "Coleman-Weinberg from anomaly-constrained matter",
        "scan_result": scan,
        "minimum_found": minimum_found,
        "minimum_a_0": minimum_a_0,
        "cw_potential_at_a0_1": cw_at_1,
        "group_used": "E8_x_E8",
        "honest_assessment": honest_assessment,
        "physics_summary": physics_summary,
    }


# ===========================================================================
# MECHANISM 2: Warped S^2 compactification
# ===========================================================================

# ---------------------------------------------------------------------------
# 2a. Warp backreaction energy
# ---------------------------------------------------------------------------

def compute_warp_backreaction(epsilon, a_0=1.0):
    """
    Compute the backreaction energy from a warped S^2 compactification.

    Warp ansatz:
        ds^2_6 = e^{2A(theta)} g_{mu nu} dx^mu dx^nu + a_0^2 d_Omega_2^2

    with A(theta) = epsilon * cos(theta).

    The warp factor contributes two terms to the effective 4D potential:

    1. Gradient energy:
        V_grad = integral over S^2 of (grad A)^2 * sqrt(g_S2)
        For A = eps*cos(theta):
            |grad A|^2 = eps^2 * sin^2(theta) / a_0^2
            sqrt(g_S2) = a_0^2 * sin(theta)
            integral = a_0^2 * eps^2 * 2*pi * integral_0^pi sin^3(theta) d_theta
                     = a_0^2 * eps^2 * 2*pi * (4/3)
                     = (8*pi/3) * eps^2 * a_0^2

        After dividing by the 4D volume normalization (~ a_0^4):
            V_grad = (8*pi/3) * eps^2 / a_0^2

    2. Curvature modification:
        The warp factor modifies the effective curvature integral.
        For small epsilon, the S^2 Ricci scalar integral becomes:
            integral R * e^{2A} sqrt(g) = (4*pi / a_0^2) * sinh(2*eps) / (2*eps)
        where the sinh(2*eps)/(2*eps) factor arises from averaging e^{2A}
        over S^2.

    Total: V_warp = V_grad + V_curv.

    Parameters
    ----------
    epsilon : float
        Warp amplitude (dimensionless).  epsilon << 1 is perturbative.
    a_0 : float
        Radius of the S^2 in Planck units (default 1.0).

    Returns
    -------
    dict with keys:
        V_grad : float -- gradient energy contribution
        V_curv : float -- curvature modification
        V_warp : float -- total warp potential
        epsilon : float
        a_0 : float
        gradient_scaling : str -- "a_0^{-2}"
        curvature_scaling : str -- "a_0^{-2}"
        is_perturbative : bool -- True if epsilon < 0.5
        description : str
    """
    # Gradient energy
    V_grad = (8.0 * math.pi / 3.0) * epsilon ** 2 / (a_0 ** 2)

    # Curvature modification
    # Unwarped curvature integral = 4*pi/a_0^2 (from chi=2 Gauss-Bonnet)
    # Warped modification factor: sinh(2*eps)/(2*eps)
    if abs(epsilon) < 1e-12:
        warp_factor = 1.0
    else:
        warp_factor = math.sinh(2.0 * epsilon) / (2.0 * epsilon)

    V_curv = -4.0 * math.pi * warp_factor / (a_0 ** 2)

    V_warp = V_grad + V_curv
    is_perturbative = abs(epsilon) < 0.5

    description = (
        f"Warped S^2 with A(theta) = {epsilon:.4f} * cos(theta), "
        f"a_0 = {a_0:.4f}: V_grad = {V_grad:.6e}, V_curv = {V_curv:.6e}, "
        f"V_warp = {V_warp:.6e}.  Both terms scale as a_0^{{-2}}."
    )

    return {
        "V_grad": V_grad,
        "V_curv": V_curv,
        "V_warp": V_warp,
        "epsilon": epsilon,
        "a_0": a_0,
        "gradient_scaling": "a_0^{-2}",
        "curvature_scaling": "a_0^{-2}",
        "is_perturbative": is_perturbative,
        "description": description,
    }


# ---------------------------------------------------------------------------
# 2b. Warped potential scan
# ---------------------------------------------------------------------------

def compute_warped_potential_scan(epsilon_values=None, a_0_values=None,
                                  N_flux=1, constants=None):
    """
    Scan over (epsilon, a_0) pairs to search for a minimum of
    V_flux(a_0) + V_warp(epsilon, a_0).

    The warp gradient energy scales as a_0^{-2}, which is one power
    slower than the flux potential (~a_0^{-4}).  In principle, the
    different scaling could create a minimum.  In practice, both terms
    are monotonically decreasing in a_0, so no minimum is expected.

    Parameters
    ----------
    epsilon_values : list of float or None
        Warp amplitudes to scan.  Default: [0.01, 0.1, 0.3, 0.5, 1.0].
    a_0_values : list of float or None
        Radius values to scan.  Default: 20 log-spaced from 0.1 to 100.
    N_flux : int
        Flux quantum number (default 1).
    constants : ignored
        Planck units used throughout.

    Returns
    -------
    dict with keys:
        scan_data : list of dict -- epsilon, a_0, V_flux, V_warp, V_total
        any_minimum_found : bool
        best_epsilon : float or None
        best_a_0 : float or None
        best_V_total : float or None
        n_epsilon : int
        n_a_0 : int
        description : str
    """
    if epsilon_values is None:
        epsilon_values = [0.01, 0.1, 0.3, 0.5, 1.0]

    if a_0_values is None:
        log_min = math.log10(0.1)
        log_max = math.log10(100.0)
        n_pts = 20
        step = (log_max - log_min) / (n_pts - 1) if n_pts > 1 else 0
        a_0_values = [10.0 ** (log_min + i * step) for i in range(n_pts)]

    scan_data = []
    any_minimum_found = False
    best_epsilon = None
    best_a_0 = None
    best_V_total = None

    for eps in epsilon_values:
        row_data = []
        for a_0 in a_0_values:
            # Flux contribution
            v_flux = None
            if _FLUX_AVAILABLE:
                flux_result = find_flux_minimum(N=N_flux, a_0=a_0,
                                                constants=constants)
                if flux_result["minimum_exists"]:
                    v_flux = flux_result["V_at_minimum"]

            # Warp contribution
            warp_result = compute_warp_backreaction(eps, a_0=a_0)
            v_warp = warp_result["V_warp"]

            v_total = None
            if v_flux is not None:
                v_total = v_flux + v_warp

            entry = {
                "epsilon": eps,
                "a_0": a_0,
                "V_flux": v_flux,
                "V_warp": v_warp,
                "V_total": v_total,
            }
            scan_data.append(entry)
            row_data.append((a_0, v_total))

        # Check for minimum in this epsilon slice
        valid = [(a, v) for a, v in row_data if v is not None]
        if len(valid) >= 3:
            for i in range(1, len(valid) - 1):
                a_prev, v_prev = valid[i - 1]
                a_curr, v_curr = valid[i]
                a_next, v_next = valid[i + 1]
                if v_curr < v_prev and v_curr < v_next:
                    any_minimum_found = True
                    if best_V_total is None or v_curr < best_V_total:
                        best_epsilon = eps
                        best_a_0 = a_curr
                        best_V_total = v_curr
                    break

    if any_minimum_found:
        description = (
            f"A minimum was found at epsilon = {best_epsilon:.4f}, "
            f"a_0 = {best_a_0:.4f} with V_total = {best_V_total:.6e}.  "
            f"The warp gradient energy creates a turning point."
        )
    else:
        description = (
            "No minimum found in V_flux + V_warp for any (epsilon, a_0) pair.  "
            "The warp gradient energy scales as a_0^{-2} and the curvature "
            "modification also scales as a_0^{-2}.  While this is a slower "
            "falloff than V_flux (~a_0^{-4}), both contributions are still "
            "monotonically decreasing in a_0.  The total potential has no "
            "turning point."
        )

    return {
        "scan_data": scan_data,
        "any_minimum_found": any_minimum_found,
        "best_epsilon": best_epsilon,
        "best_a_0": best_a_0,
        "best_V_total": best_V_total,
        "n_epsilon": len(epsilon_values),
        "n_a_0": len(a_0_values),
        "description": description,
    }


# ---------------------------------------------------------------------------
# 2c. Warp mechanism summary
# ---------------------------------------------------------------------------

def summarize_warp_mechanism(constants=None):
    """
    Run the warped S^2 radius-fixing mechanism and return a full summary.

    Parameters
    ----------
    constants : ignored
        Planck units used throughout.

    Returns
    -------
    dict with keys:
        mechanism_name : str
        scan_result : dict from compute_warped_potential_scan
        minimum_found : bool
        single_point : dict -- warp backreaction at epsilon=0.1, a_0=1
        honest_assessment : str
        physics_summary : str
    """
    single_point = compute_warp_backreaction(epsilon=0.1, a_0=1.0)
    scan = compute_warped_potential_scan(N_flux=1, constants=constants)

    minimum_found = scan["any_minimum_found"]

    if minimum_found:
        honest_assessment = (
            f"The warped S^2 mechanism produces a minimum at "
            f"epsilon = {scan['best_epsilon']:.4f}, "
            f"a_0 = {scan['best_a_0']:.4f}.  This is an unexpected result: "
            f"the a_0^{{-2}} warp energy was expected to remain monotonic.  "
            f"The result should be checked with finer grids and different "
            f"epsilon values."
        )
    else:
        honest_assessment = (
            "The warped S^2 mechanism does NOT fix a_0, as expected.  "
            "The warp ansatz A(theta) = epsilon*cos(theta) introduces "
            "gradient energy (~a_0^{-2}) and modifies the curvature integral "
            "(also ~a_0^{-2}).  Both terms scale with the SAME power of a_0 "
            "and both decrease monotonically.  Adding a monotonic a_0^{-2} "
            "term to the monotonic a_0^{-4} flux potential cannot create a "
            "minimum.  This is a genuine negative result: warping the S^2 "
            "does not break the flat direction."
        )

    physics_summary = (
        "A warped compactification ds^2 = e^{2A(theta)} g_{mu nu} dx^mu dx^nu "
        "+ a_0^2 d_Omega_2^2 with A(theta) = epsilon*cos(theta) adds gradient "
        "energy from (nabla A)^2 and modifies the curvature integral.  "
        "The gradient energy = (8*pi/3)*eps^2/a_0^2 and the curvature = "
        "-4*pi*sinh(2*eps)/(2*eps*a_0^2).  Both scale as a_0^{-2}, which is "
        "one power slower than the a_0^{-4} flux potential.  Despite the "
        "different scaling, the combination remains monotonic in a_0."
    )

    return {
        "mechanism_name": "Warped S^2 compactification",
        "scan_result": scan,
        "minimum_found": minimum_found,
        "single_point": single_point,
        "honest_assessment": honest_assessment,
        "physics_summary": physics_summary,
    }


# ===========================================================================
# MECHANISM 3: Orbifold S^2/Z_2
# ===========================================================================

# ---------------------------------------------------------------------------
# 3a. Orbifold Casimir energy
# ---------------------------------------------------------------------------

def compute_orbifold_casimir(a_0=1.0, l_max=100):
    """
    Compute the Casimir energy on the orbifold S^2/Z_2.

    The Z_2 identification (antipodal map) on S^2 projects out odd-l
    KK modes, keeping only even-l modes: l = 0, 2, 4, 6, ...

    The spectral zeta function on S^2/Z_2 is:

        zeta_{orb}(s) = sum_{l=0,2,4,...}^{l_max} (2l+1) * [l(l+1)]^{-s}

    (skipping l=0 for the zeta function, since l(l+1)=0 at l=0).

    The Casimir coefficient is determined by zeta_{orb}(-1/2), which may
    differ in sign from the full S^2 zeta function zeta_S2(-1/2) < 0.

    We compute this numerically using Euler-Maclaurin regularization
    on the even-l subseries.

    Parameters
    ----------
    a_0 : float
        Radius of the S^2/Z_2 in Planck units (default 1.0).
    l_max : int
        Maximum l in the sum (default 100).  Only even l are included.

    Returns
    -------
    dict with keys:
        zeta_orb_value : float -- regularized zeta_{orb}(-1/2)
        casimir_coefficient : float -- N_eff * zeta_orb / (4*pi)
        V_casimir_orb : float -- Casimir energy at this a_0
        n_modes_included : int -- number of even-l modes summed
        differs_from_s2 : bool -- True if sign differs from S^2
        a_0 : float
        l_max : int
        scaling_power : int -- -4
        description : str
    """
    # Compute zeta_{orb}(-1/2) via direct summation + Euler-Maclaurin
    # subtraction on the even-l subseries.
    #
    # For even l = 2k with k = 1, 2, ..., l_max//2:
    #   eigenvalue = l(l+1) = 2k(2k+1)
    #   degeneracy = 2l+1 = 4k+1
    #   term = (4k+1) * [2k(2k+1)]^{1/2}  (for s = -1/2)

    raw_sum = 0.0
    divergent_piece = 0.0
    n_modes = 0
    k_max = l_max // 2

    for k in range(1, k_max + 1):
        l = 2 * k
        eigenval = l * (l + 1)
        deg = 2 * l + 1
        term = deg * math.sqrt(eigenval)
        raw_sum += term
        n_modes += 1

        # Leading asymptotic: (4k+1)*sqrt(2k(2k+1)) ~ (4k)*(2k) = 8*k^2
        if k >= 2:
            asymp = 8.0 * k * k
            divergent_piece += asymp

    zeta_orb = raw_sum - divergent_piece

    # Casimir coefficient: C_cas = N_eff * zeta_orb / (4*pi)
    # N_eff = 9 for 6D graviton (same as on S^2)
    N_eff = 9
    casimir_coeff = N_eff * zeta_orb / (4.0 * math.pi)

    # Casimir energy: V_cas = C_cas / a_0^4
    V_casimir_orb = casimir_coeff / (a_0 ** 4)

    # Compare with full S^2: zeta_S2(-1/2) ~ -0.25 (negative)
    # The orbifold sum has fewer terms, so the sign could change
    s2_zeta_negative = True  # known from casimir_stabilization module
    differs = (zeta_orb > 0) != (not s2_zeta_negative)

    description = (
        f"Orbifold Casimir on S^2/Z_2 (a_0 = {a_0:.4f}): "
        f"zeta_orb(-1/2) = {zeta_orb:.6f} "
        f"({'positive' if zeta_orb > 0 else 'negative'}), "
        f"C_cas = {casimir_coeff:.6f}, V_cas = {V_casimir_orb:.6e}.  "
        f"{n_modes} even-l modes summed (l = 2, 4, ..., {2 * k_max})."
    )

    return {
        "zeta_orb_value": zeta_orb,
        "casimir_coefficient": casimir_coeff,
        "V_casimir_orb": V_casimir_orb,
        "n_modes_included": n_modes,
        "differs_from_s2": differs,
        "a_0": a_0,
        "l_max": l_max,
        "scaling_power": -4,
        "description": description,
    }


# ---------------------------------------------------------------------------
# 3b. Brane tension potential
# ---------------------------------------------------------------------------

def compute_brane_tension_potential(T_0, a_0=1.0):
    """
    Compute the brane tension contribution from orbifold fixed points.

    The S^2/Z_2 orbifold has 2 fixed points (north and south poles).
    Each fixed point can support a codimension-2 brane with tension T_0.
    The brane tension contributes to the 4D effective potential as:

        V_brane = 2 * T_0 / a_0^2

    The factor of 2 comes from the two fixed points.  The a_0^{-2}
    scaling arises from the brane codimension: a codimension-2 brane
    in the internal space contributes energy ~ T / Vol(M_2) ~ T / a_0^2.

    This is the KEY new term: it scales as a_0^{-2}, which is SLOWER
    than the Casimir a_0^{-4}.  The competition between these two
    scalings can create a minimum.

    Parameters
    ----------
    T_0 : float
        Brane tension at each fixed point (in Planck units).
    a_0 : float
        Radius of the S^2/Z_2 in Planck units (default 1.0).

    Returns
    -------
    dict with keys:
        V_brane : float -- brane tension potential
        T_0 : float -- brane tension parameter
        a_0 : float
        n_fixed_points : int -- 2
        scaling_power : int -- -2
        description : str
    """
    V_brane = 2.0 * T_0 / (a_0 ** 2)

    description = (
        f"Brane tension potential: V_brane = 2 * T_0 / a_0^2 = "
        f"2 * {T_0:.6e} / {a_0:.4f}^2 = {V_brane:.6e}.  "
        f"Two fixed points on S^2/Z_2 (north and south poles), "
        f"each with tension T_0 = {T_0:.6e}."
    )

    return {
        "V_brane": V_brane,
        "T_0": T_0,
        "a_0": a_0,
        "n_fixed_points": 2,
        "scaling_power": -2,
        "description": description,
    }


# ---------------------------------------------------------------------------
# 3c. Total orbifold potential
# ---------------------------------------------------------------------------

def compute_orbifold_total_potential(a_0, T_0=0.01, N_flux=1, l_max=100):
    """
    Compute the total effective potential on the orbifold S^2/Z_2.

    V_total = V_casimir_orb(a_0) + V_flux(a_0) + V_brane(T_0, a_0)

    The key competition:
        - V_casimir ~ a_0^{-4}  (from KK tower on S^2/Z_2, even-l only)
        - V_flux ~ a_0^{-4}     (from V_min of flux-stabilized sigma)
        - V_brane ~ a_0^{-2}    (from brane tension at fixed points)

    At small a_0: V_casimir + V_flux dominate (faster falloff -> larger)
    At large a_0: V_brane dominates (slower falloff -> larger relative)

    If V_casimir + V_flux is negative and V_brane is positive, then
    V_total can change sign, creating a minimum at the crossover.

    Parameters
    ----------
    a_0 : float
        Radius of the S^2/Z_2 in Planck units.
    T_0 : float
        Brane tension (default 0.01 Planck units).
    N_flux : int
        Flux quantum number (default 1).
    l_max : int
        Truncation of KK tower (default 100).

    Returns
    -------
    dict with keys:
        V_total : float -- total potential
        V_casimir_orb : float -- orbifold Casimir contribution
        V_flux : float or None -- flux contribution
        V_brane : float -- brane tension contribution
        a_0 : float
        T_0 : float
        N_flux : int
        description : str
    """
    # Orbifold Casimir
    cas_result = compute_orbifold_casimir(a_0=a_0, l_max=l_max)
    v_casimir = cas_result["V_casimir_orb"]

    # Flux contribution (from find_flux_minimum at this a_0)
    v_flux = None
    if _FLUX_AVAILABLE:
        flux_result = find_flux_minimum(N=N_flux, a_0=a_0)
        if flux_result["minimum_exists"]:
            v_flux = flux_result["V_at_minimum"]

    # Brane tension
    brane_result = compute_brane_tension_potential(T_0, a_0=a_0)
    v_brane = brane_result["V_brane"]

    # Total
    v_total = v_casimir + v_brane
    if v_flux is not None:
        v_total += v_flux

    description = (
        f"Orbifold potential at a_0 = {a_0:.4f}: "
        f"V_cas = {v_casimir:.6e}, V_flux = {v_flux}, "
        f"V_brane = {v_brane:.6e}, V_total = {v_total:.6e}."
    )

    return {
        "V_total": v_total,
        "V_casimir_orb": v_casimir,
        "V_flux": v_flux,
        "V_brane": v_brane,
        "a_0": a_0,
        "T_0": T_0,
        "N_flux": N_flux,
        "description": description,
    }


# ---------------------------------------------------------------------------
# 3d. Orbifold radius scan
# ---------------------------------------------------------------------------

def compute_orbifold_radius_scan(T_0_values=None, a_0_values=None,
                                  N_flux=1, l_max=100):
    """
    Scan over brane tensions and a_0 values to find the orbifold mechanism
    sweet spot.

    For each T_0 value, scan a_0 and look for a minimum in V_total(a_0).
    The a_0^{-2} brane tension competes with the a_0^{-4} Casimir+flux
    at different scales, potentially creating a minimum.

    Parameters
    ----------
    T_0_values : list of float or None
        Brane tensions to scan.  Default: [0.001, 0.01, 0.1, 1.0, 10.0].
    a_0_values : list of float or None
        Radius values to scan.  Default: 20 log-spaced from 0.1 to 100.
    N_flux : int
        Flux quantum number (default 1).
    l_max : int
        Truncation of KK tower (default 100).

    Returns
    -------
    dict with keys:
        scan_data : list of dict -- T_0, a_0, V_total for each point
        minima_found : list of dict -- T_0, a_0, V_total where minimum exists
        any_minimum_found : bool
        best_T_0 : float or None
        best_a_0 : float or None
        best_V_total : float or None
        n_T_0 : int
        n_a_0 : int
        description : str
    """
    if T_0_values is None:
        T_0_values = [0.001, 0.01, 0.1, 1.0, 10.0]

    if a_0_values is None:
        log_min = math.log10(0.1)
        log_max = math.log10(100.0)
        n_pts = 20
        step = (log_max - log_min) / (n_pts - 1) if n_pts > 1 else 0
        a_0_values = [10.0 ** (log_min + i * step) for i in range(n_pts)]

    scan_data = []
    minima_found = []
    best_T_0 = None
    best_a_0 = None
    best_V_total = None

    for T_0 in T_0_values:
        row_data = []
        for a_0 in a_0_values:
            result = compute_orbifold_total_potential(
                a_0=a_0, T_0=T_0, N_flux=N_flux, l_max=l_max
            )
            v_total = result["V_total"]

            scan_data.append({
                "T_0": T_0,
                "a_0": a_0,
                "V_total": v_total,
                "V_casimir_orb": result["V_casimir_orb"],
                "V_flux": result["V_flux"],
                "V_brane": result["V_brane"],
            })
            row_data.append((a_0, v_total))

        # Check for minimum in this T_0 slice
        if len(row_data) >= 3:
            for i in range(1, len(row_data) - 1):
                a_prev, v_prev = row_data[i - 1]
                a_curr, v_curr = row_data[i]
                a_next, v_next = row_data[i + 1]

                if v_curr < v_prev and v_curr < v_next:
                    minima_found.append({
                        "T_0": T_0,
                        "a_0": a_curr,
                        "V_total": v_curr,
                    })
                    if best_V_total is None or v_curr < best_V_total:
                        best_T_0 = T_0
                        best_a_0 = a_curr
                        best_V_total = v_curr
                    break

    any_minimum_found = len(minima_found) > 0

    if any_minimum_found:
        description = (
            f"Orbifold mechanism: {len(minima_found)} minimum/minima found "
            f"across {len(T_0_values)} brane tensions.  Best minimum at "
            f"T_0 = {best_T_0:.4e}, a_0 = {best_a_0:.4f} with "
            f"V_total = {best_V_total:.6e}.  The brane tension a_0^{{-2}} "
            f"competes with Casimir+flux a_0^{{-4}} to create a turning point."
        )
    else:
        description = (
            f"Orbifold mechanism: no minimum found across "
            f"{len(T_0_values)} brane tensions and "
            f"{len(a_0_values)} a_0 values.  The brane tension term "
            f"(~a_0^{{-2}}) does not overcome the Casimir+flux terms "
            f"(~a_0^{{-4}}) in the scanned range to create a turning point."
        )

    return {
        "scan_data": scan_data,
        "minima_found": minima_found,
        "any_minimum_found": any_minimum_found,
        "best_T_0": best_T_0,
        "best_a_0": best_a_0,
        "best_V_total": best_V_total,
        "n_T_0": len(T_0_values),
        "n_a_0": len(a_0_values),
        "description": description,
    }


# ---------------------------------------------------------------------------
# 3e. Find critical brane tension
# ---------------------------------------------------------------------------

def find_critical_brane_tension(N_flux=1, l_max=100):
    """
    Binary search for the critical brane tension T_0* where a minimum
    in V_total(a_0) first appears on the orbifold S^2/Z_2.

    The idea: for T_0 = 0, V_total ~ a_0^{-4} (monotonic, no minimum).
    As T_0 increases, the a_0^{-2} brane term grows.  At some critical
    T_0*, the potential develops an inflection point, and for T_0 > T_0*
    a genuine minimum appears.

    We search by checking whether a minimum exists at a coarse grid of
    a_0 values, and use binary search on T_0.

    Parameters
    ----------
    N_flux : int
        Flux quantum number (default 1).
    l_max : int
        Truncation of KK tower (default 100).

    Returns
    -------
    dict with keys:
        T_0_critical : float or None -- critical brane tension
        minimum_exists_above : bool
        minimum_exists_below : bool
        search_converged : bool
        n_iterations : int
        a_0_at_inflection : float or None
        description : str
    """
    # Helper: does a minimum exist for given T_0?
    def _has_minimum(T_0):
        log_min = math.log10(0.1)
        log_max = math.log10(100.0)
        n_pts = 30
        step = (log_max - log_min) / (n_pts - 1) if n_pts > 1 else 0
        a_0_values = [10.0 ** (log_min + i * step) for i in range(n_pts)]

        v_values = []
        for a_0 in a_0_values:
            result = compute_orbifold_total_potential(
                a_0=a_0, T_0=T_0, N_flux=N_flux, l_max=l_max
            )
            v_values.append((a_0, result["V_total"]))

        # Check for local minimum
        for i in range(1, len(v_values) - 1):
            _, v_prev = v_values[i - 1]
            a_curr, v_curr = v_values[i]
            _, v_next = v_values[i + 1]
            if v_curr < v_prev and v_curr < v_next:
                return True, a_curr
        return False, None

    # Binary search
    T_lo = 0.0
    T_hi = 100.0
    max_iter = 40
    tol = 1e-6

    # First check: does T_hi have a minimum?
    has_hi, a_hi = _has_minimum(T_hi)
    if not has_hi:
        # Even large T_0 doesn't produce a minimum
        return {
            "T_0_critical": None,
            "minimum_exists_above": False,
            "minimum_exists_below": False,
            "search_converged": False,
            "n_iterations": 0,
            "a_0_at_inflection": None,
            "description": (
                "No minimum found even at T_0 = 100.  The orbifold brane "
                "tension does not create a turning point in the scanned "
                "a_0 range [0.1, 100] Planck units."
            ),
        }

    # T_lo = 0 has no minimum (monotonic), T_hi has minimum
    # Binary search for transition
    n_iter = 0
    a_0_inflection = None

    for _ in range(max_iter):
        n_iter += 1
        T_mid = (T_lo + T_hi) / 2.0
        has_mid, a_mid = _has_minimum(T_mid)

        if has_mid:
            T_hi = T_mid
            a_0_inflection = a_mid
        else:
            T_lo = T_mid

        if (T_hi - T_lo) < tol:
            break

    T_critical = (T_lo + T_hi) / 2.0
    converged = (T_hi - T_lo) < tol

    description = (
        f"Critical brane tension T_0* = {T_critical:.6e} (Planck units).  "
        f"For T_0 > T_0*, the orbifold potential V_total(a_0) develops a "
        f"local minimum.  For T_0 < T_0*, the potential is monotonic.  "
        f"Binary search converged in {n_iter} iterations "
        f"({'converged' if converged else 'not fully converged'}).  "
        f"The inflection point appears near a_0 ~ {a_0_inflection:.4f} "
        f"Planck units."
        if a_0_inflection is not None else
        f"Critical brane tension T_0* = {T_critical:.6e}.  "
        f"Search used {n_iter} iterations."
    )

    return {
        "T_0_critical": T_critical,
        "minimum_exists_above": True,
        "minimum_exists_below": False,
        "search_converged": converged,
        "n_iterations": n_iter,
        "a_0_at_inflection": a_0_inflection,
        "description": description,
    }


# ---------------------------------------------------------------------------
# 3f. Orbifold mechanism summary
# ---------------------------------------------------------------------------

def summarize_orbifold_mechanism(constants=None):
    """
    Run the orbifold S^2/Z_2 radius-fixing mechanism and return a full
    summary.

    This is the BEST candidate mechanism because the brane tension
    (a_0^{-2}) has a genuinely different scaling from the Casimir+flux
    (a_0^{-4}).  The competition between these scalings can create a
    minimum at a_0 where the two contributions balance.

    Parameters
    ----------
    constants : ignored
        Planck units used throughout.

    Returns
    -------
    dict with keys:
        mechanism_name : str
        orbifold_casimir : dict from compute_orbifold_casimir at a_0=1
        brane_tension_example : dict from compute_brane_tension_potential
        radius_scan : dict from compute_orbifold_radius_scan
        critical_tension : dict from find_critical_brane_tension
        minimum_found : bool
        best_a_0 : float or None
        best_T_0 : float or None
        honest_assessment : str
        physics_summary : str
    """
    orb_casimir = compute_orbifold_casimir(a_0=1.0, l_max=100)
    brane_example = compute_brane_tension_potential(T_0=0.01, a_0=1.0)
    scan = compute_orbifold_radius_scan(N_flux=1, l_max=100)
    critical = find_critical_brane_tension(N_flux=1, l_max=100)

    minimum_found = scan["any_minimum_found"]
    best_a_0 = scan["best_a_0"]
    best_T_0 = scan["best_T_0"]

    if minimum_found:
        honest_assessment = (
            f"The orbifold S^2/Z_2 mechanism DOES produce a minimum at "
            f"a_0 = {best_a_0:.4f} Planck units for brane tension "
            f"T_0 = {best_T_0:.4e}.  The critical brane tension is "
            f"T_0* = {critical['T_0_critical']:.4e} (Planck units).  "
            f"This is a genuine radius-fixing mechanism, but the required "
            f"brane tension is a NEW free parameter.  The value of a_0 "
            f"depends on T_0, trading one unknown for another.  To make "
            f"this predictive would require deriving T_0 from first "
            f"principles (e.g., from anomaly cancellation or flux "
            f"quantization on the brane)."
        )
    else:
        honest_assessment = (
            "The orbifold S^2/Z_2 mechanism does NOT produce a minimum "
            "in the scanned range.  The brane tension (a_0^{-2}) and "
            "Casimir+flux (a_0^{-4}) terms have different scalings, which "
            "in principle allows a minimum.  However, in the numerical scan "
            "with T_0 in [0.001, 10] and a_0 in [0.1, 100], no turning "
            "point was found.  The critical brane tension "
            + (f"T_0* = {critical['T_0_critical']:.4e} "
               if critical['T_0_critical'] is not None
               else "was not found in the search range ")
            + "may require a_0 values outside the scanned range, or the "
            "orbifold Casimir coefficient has the wrong sign for this "
            "mechanism to work."
        )

    physics_summary = (
        "The orbifold S^2/Z_2 projects out odd-l KK modes, modifying the "
        "Casimir coefficient.  The two Z_2 fixed points (poles of S^2) "
        "support codimension-2 branes with tension T_0, contributing "
        "V_brane = 2*T_0/a_0^2 to the potential.  The crucial feature is "
        "the a_0^{-2} scaling of V_brane vs a_0^{-4} for Casimir+flux.  "
        "At small a_0, the a_0^{-4} terms dominate; at large a_0, the "
        "a_0^{-2} brane tension dominates.  If these have opposite signs, "
        "a crossover produces a minimum.  The orbifold Casimir coefficient "
        f"zeta_orb(-1/2) = {orb_casimir['zeta_orb_value']:.6f}."
    )

    return {
        "mechanism_name": "Orbifold S^2/Z_2 with brane tension",
        "orbifold_casimir": orb_casimir,
        "brane_tension_example": brane_example,
        "radius_scan": scan,
        "critical_tension": critical,
        "minimum_found": minimum_found,
        "best_a_0": best_a_0,
        "best_T_0": best_T_0,
        "honest_assessment": honest_assessment,
        "physics_summary": physics_summary,
    }


# ===========================================================================
# OVERALL SUMMARY
# ===========================================================================

def summarize_radius_fixing(constants=None):
    """
    Run all three radius-fixing mechanisms and return a combined summary.

    This is the main entry point for the radius fixing analysis.  It
    evaluates:
        1. Coleman-Weinberg from anomaly-constrained matter
        2. Warped S^2 compactification
        3. Orbifold S^2/Z_2 with brane tension

    and determines which (if any) mechanism successfully fixes a_0.

    Parameters
    ----------
    constants : ignored
        Planck units used throughout.

    Returns
    -------
    dict with keys:
        mechanism_1_cw : dict from summarize_cw_mechanism
        mechanism_2_warp : dict from summarize_warp_mechanism
        mechanism_3_orbifold : dict from summarize_orbifold_mechanism
        any_mechanism_fixes_radius : bool
        best_mechanism : str or None
        overall_assessment : str
        honest_assessment : str
    """
    mech1 = summarize_cw_mechanism(constants=constants)
    mech2 = summarize_warp_mechanism(constants=constants)
    mech3 = summarize_orbifold_mechanism(constants=constants)

    results = [
        ("Coleman-Weinberg (anomaly matter)", mech1["minimum_found"]),
        ("Warped S^2", mech2["minimum_found"]),
        ("Orbifold S^2/Z_2", mech3["minimum_found"]),
    ]

    any_fixes = any(success for _, success in results)
    best_mechanism = None
    for name, success in results:
        if success:
            best_mechanism = name
            break

    # Build overall assessment
    n_success = sum(1 for _, s in results if s)
    mechanism_status = []
    for name, success in results:
        status = "SUCCEEDS" if success else "FAILS"
        mechanism_status.append(f"  - {name}: {status}")
    status_lines = "\n".join(mechanism_status)

    overall_assessment = (
        f"Three mechanisms tested to fix the internal radius a_0 of S^2:\n"
        f"{status_lines}\n"
        f"{n_success} of 3 mechanisms produce a minimum.  "
        + (f"Best mechanism: {best_mechanism}."
           if best_mechanism else
           "No mechanism fixes a_0 within the scanned parameter ranges.")
    )

    honest_assessment = (
        "The internal radius a_0 remains the central open problem of the "
        "6D Alpha Ladder framework.  The scaling symmetry (n=2 makes EH "
        "scale-invariant and GB topological) is robust and cannot be "
        "broken by simple modifications within the minimal framework.  "
        "The three mechanisms tested here represent increasing levels of "
        "beyond-minimal physics: (1) matter loops require specifying the "
        "full gauge group and are loop-suppressed; (2) warping adds "
        "gradient energy but with the same monotonic behavior; (3) the "
        "orbifold introduces brane tension as a new parameter.  Even if "
        "mechanism (3) succeeds, it trades one unknown (a_0) for another "
        "(T_0).  A fully predictive determination of a_0 likely requires "
        "UV-complete information (e.g., from string theory) that is "
        "beyond the scope of the effective field theory framework."
    )

    return {
        "mechanism_1_cw": mech1,
        "mechanism_2_warp": mech2,
        "mechanism_3_orbifold": mech3,
        "any_mechanism_fixes_radius": any_fixes,
        "best_mechanism": best_mechanism,
        "overall_assessment": overall_assessment,
        "honest_assessment": honest_assessment,
    }
