"""
Casimir stabilization of the dilaton mass.

Computes the 1-loop Casimir energy from the KK tower of the 6D graviton
on a compact 2-manifold (S^2 or genus-2 surface) and determines whether
it can stabilize the volume modulus sigma.

The tree-level KK reduction of the 6D Einstein-Hilbert action on a compact
2-manifold gives a flat direction for the volume modulus sigma because the
integrated internal curvature is topological (Gauss-Bonnet theorem).  The
Casimir energy from the KK tower provides a 1-loop correction that could
break this flat direction.

The effective 4D potential is:
    V_eff(sigma) = A * e^{4*sigma} + B * e^{2*sigma}

where:
    B = -(1/2) * chi / a_0^2   (curvature contribution, M_Pl = 1)
    A = N_eff * Z_reg / (4*pi*a_0^4)   (Casimir contribution)

The spectral zeta function on S^2 is:
    zeta_S2(s) = sum_{l>=1} (2l+1) * [l(l+1)]^{-s}

The Casimir energy involves zeta_S2(-1/2), which requires analytic
continuation.  The honest result is that zeta_S2(-1/2) < 0 (negative
Casimir coefficient), so Casimir stabilization from the graviton tower
alone does NOT produce a stable minimum.

All calculations use pure Python math (no numpy/scipy).
"""

import math
from decimal import Decimal, getcontext

try:
    from alpha_ladder_core.anomaly_cancellation import _ANOMALY_FREE_GROUPS
    _ANOMALY_AVAILABLE = True
except ImportError:
    _ANOMALY_AVAILABLE = False
    _ANOMALY_FREE_GROUPS = {}

getcontext().prec = 50


# ---------------------------------------------------------------------------
# Helper: Bernoulli polynomials B_n(x) for n = 0..6
# ---------------------------------------------------------------------------

def _bernoulli_polynomial(n, x):
    """
    Compute the Bernoulli polynomial B_n(x) for n = 0, 1, ..., 6.

    Uses the explicit formulas for the first 7 Bernoulli polynomials.
    These are needed for the Hurwitz zeta function at negative integers.

    Parameters
    ----------
    n : int -- order of the Bernoulli polynomial (0 <= n <= 6)
    x : float -- argument

    Returns
    -------
    float -- B_n(x)
    """
    if n == 0:
        return 1.0
    elif n == 1:
        return x - 0.5
    elif n == 2:
        return x ** 2 - x + 1.0 / 6.0
    elif n == 3:
        return x ** 3 - 1.5 * x ** 2 + 0.5 * x
    elif n == 4:
        return x ** 4 - 2.0 * x ** 3 + x ** 2 - 1.0 / 30.0
    elif n == 5:
        return x ** 5 - 2.5 * x ** 4 + (5.0 / 3.0) * x ** 3 - x / 6.0
    elif n == 6:
        return (x ** 6 - 3.0 * x ** 5 + 2.5 * x ** 4
                - 0.5 * x ** 2 + 1.0 / 42.0)
    else:
        raise ValueError(f"Bernoulli polynomial B_{n} not implemented (n must be 0..6)")


# ---------------------------------------------------------------------------
# Helper: Hurwitz zeta at negative integers via Bernoulli polynomials
# ---------------------------------------------------------------------------

def _hurwitz_zeta_negative_int(n, a):
    """
    Compute zeta_H(-n, a) = -B_{n+1}(a) / (n+1) for non-negative integer n.

    This is the standard identity relating the Hurwitz zeta function at
    negative integers to the Bernoulli polynomials.

    Parameters
    ----------
    n : int -- non-negative integer; evaluates zeta_H(-n, a)
    a : float -- the shift parameter of the Hurwitz zeta

    Returns
    -------
    float -- zeta_H(-n, a)
    """
    if n < 0:
        raise ValueError(f"n must be non-negative, got {n}")
    if n + 1 > 6:
        raise ValueError(
            f"Need B_{n + 1}(a) but only B_0..B_6 are implemented"
        )
    return -_bernoulli_polynomial(n + 1, a) / (n + 1)


# ---------------------------------------------------------------------------
# Helper: Hurwitz zeta for Re(s) > 1 via partial sums
# ---------------------------------------------------------------------------

def _hurwitz_zeta_positive(s, a, terms=200):
    """
    Compute zeta_H(s, a) = sum_{k=0}^{terms-1} (k + a)^{-s} for Re(s) > 1.

    This is a direct partial sum approximation, valid when s > 1 and a > 0.
    The truncation error is O(terms^{1-s}).

    Parameters
    ----------
    s : float -- the exponent (must be > 1 for convergence)
    a : float -- the shift parameter (must be > 0)
    terms : int -- number of terms in the partial sum

    Returns
    -------
    float -- approximate zeta_H(s, a)
    """
    if s <= 1.0:
        raise ValueError(
            f"Hurwitz zeta partial sum requires s > 1, got s = {s}"
        )
    if a <= 0:
        raise ValueError(f"Hurwitz zeta requires a > 0, got a = {a}")

    total = 0.0
    for k in range(terms):
        total += (k + a) ** (-s)
    return total


# ---------------------------------------------------------------------------
# Helper: generalized binomial coefficient C(alpha, k)
# ---------------------------------------------------------------------------

def _generalized_binomial(alpha, k):
    """
    Compute the generalized binomial coefficient C(alpha, k) for real alpha
    and non-negative integer k.

    C(alpha, k) = alpha * (alpha - 1) * ... * (alpha - k + 1) / k!

    Parameters
    ----------
    alpha : float -- the upper index (can be any real number)
    k : int -- the lower index (non-negative integer)

    Returns
    -------
    float -- C(alpha, k)
    """
    if k < 0:
        return 0.0
    if k == 0:
        return 1.0
    result = 1.0
    for i in range(k):
        result *= (alpha - i)
    for i in range(1, k + 1):
        result /= i
    return result


# ---------------------------------------------------------------------------
# 1. KK spectrum on S^2
# ---------------------------------------------------------------------------

def compute_kk_spectrum_s2(l_max=100):
    """
    Compute KK eigenvalues and field content for the 6D graviton on S^2.

    The eigenvalues of the scalar Laplacian on S^2 of radius a are:
        m_l^2 = l(l+1) / a^2,    degeneracy = 2l + 1

    The 6D graviton has D(D-3)/2 = 6*3/2 = 9 physical polarizations.
    Upon KK reduction on S^2, these decompose into 4D fields:
        - 1 graviton (2 dof): the massless 4D graviton (l=0 mode)
        - 2 vectors (4 dof): gravi-photons from mixed components
        - 3 scalars (3 dof): from internal components of the 6D metric

    For l >= 1, each KK level contains massive spin-2, spin-1, and spin-0
    fields with total 9 polarizations per level.

    Parameters
    ----------
    l_max : int -- maximum angular momentum quantum number (default 100)

    Returns
    -------
    dict with keys:
        eigenvalues : list of float -- m_l^2 = l(l+1) for l = 0..l_max
        degeneracies : list of int -- 2l+1 for l = 0..l_max
        field_content : dict describing the 6D graviton decomposition
        l_max : int
        n_massive_modes : int -- total massive degrees of freedom
        description : str
    """
    eigenvalues = []
    degeneracies = []
    n_massive_dof = 0

    for l in range(l_max + 1):
        eigenvalues.append(l * (l + 1))
        degeneracies.append(2 * l + 1)
        if l >= 1:
            n_massive_dof += (2 * l + 1)

    # 6D graviton polarizations: D(D-3)/2 with D = 6
    D = 6
    total_polarizations = D * (D - 3) // 2  # = 9

    field_content = {
        "D": D,
        "total_polarizations": total_polarizations,
        "decomposition": {
            "graviton_4D": {
                "spin": 2,
                "dof": 2,
                "description": "massless 4D graviton from l=0 mode",
            },
            "vectors": {
                "spin": 1,
                "count": 2,
                "dof_per_vector": 2,
                "dof_total": 4,
                "description": (
                    "gravi-photons from g_{mu a} components "
                    "(2 internal indices on S^2)"
                ),
            },
            "scalars": {
                "spin": 0,
                "count": 3,
                "dof_total": 3,
                "description": (
                    "scalars from g_{ab} components: 1 trace (breathing mode) "
                    "+ 2 from symmetric traceless tensor on S^2"
                ),
            },
        },
        "check_total": 2 + 4 + 3,  # = 9
        "consistent": (2 + 4 + 3) == total_polarizations,
    }

    return {
        "eigenvalues": eigenvalues,
        "degeneracies": degeneracies,
        "field_content": field_content,
        "l_max": l_max,
        "n_massive_modes": n_massive_dof,
        "description": (
            f"KK spectrum on S^2: {l_max + 1} levels (l = 0..{l_max}).  "
            f"6D graviton has {total_polarizations} polarizations decomposing "
            f"into 4D graviton (2) + vectors (4) + scalars (3).  "
            f"Total massive dof summed over l >= 1: {n_massive_dof}."
        ),
    }


# ---------------------------------------------------------------------------
# 2. Spectral zeta function on S^2
# ---------------------------------------------------------------------------

def compute_spectral_zeta_s2(s, l_max=1000):
    """
    Compute the spectral zeta function on S^2:

        zeta_S2(s) = sum_{l >= 1} (2l+1) * [l(l+1)]^{-s}

    For Re(s) > 1, this is computed as a direct partial sum.

    For s = -1/2, analytic continuation is required.  We use the
    substitution u = l + 1/2, so l(l+1) = u^2 - 1/4, and expand
    [u^2 - 1/4]^{-s} = u^{-2s} * [1 - 1/(4u^2)]^{-s} in a binomial
    series.  This gives:

        zeta_S2(s) = 2 * sum_{j=0}^J C(s+j-1, j) * (1/4)^j
                       * zeta_H(2s+2j-1, 3/2)

    where C is the generalized binomial coefficient and zeta_H is the
    Hurwitz zeta function.  The factor of 2 comes from the degeneracy
    2l+1 = 2u in the new variable.

    For s = -1/2, the argument of the Hurwitz zeta is:
        2*(-1/2) + 2j - 1 = 2j - 2

    So j=0 gives zeta_H(-2, 3/2), j=1 gives zeta_H(0, 3/2), etc.
    These are evaluated using the Bernoulli polynomial identity for
    non-positive integers and partial sums for positive arguments.

    Parameters
    ----------
    s : float -- the exponent
    l_max : int -- maximum l for partial sums (default 1000)

    Returns
    -------
    dict with keys:
        zeta_value : float -- the computed value
        s : float
        l_max : int
        method : str -- "partial_sum" or "analytic_continuation"
        partial_sum : float or None -- partial sum if also computed
        convergence_diagnostics : dict
    """
    # Direct partial sum (for any s, but only converges for Re(s) > 1)
    partial_sum = 0.0
    for l in range(1, l_max + 1):
        partial_sum += (2 * l + 1) * (l * (l + 1)) ** (-s)

    if s > 1.0:
        # Partial sum converges; also compute tail estimate
        # Last term as fraction of total gives convergence diagnostic
        last_term = (2 * l_max + 1) * (l_max * (l_max + 1)) ** (-s)
        relative_last = abs(last_term / partial_sum) if partial_sum != 0 else 0.0

        return {
            "zeta_value": partial_sum,
            "s": s,
            "l_max": l_max,
            "method": "partial_sum",
            "partial_sum": partial_sum,
            "convergence_diagnostics": {
                "last_term": last_term,
                "relative_last_term": relative_last,
                "converged": relative_last < 1e-10,
                "description": (
                    f"Direct sum of {l_max} terms.  "
                    f"Last term / total = {relative_last:.2e}."
                ),
            },
        }

    # --- Analytic continuation for s <= 1 ---
    # Use the Hurwitz zeta expansion with J+1 terms
    J = 7  # use 8 terms (j = 0..7)

    # zeta_S2(s) = 2 * sum_{j=0}^J c_j(s) * zeta_H(2s + 2j - 1, 3/2)
    # where c_j(s) = C(s+j-1, j) * (1/4)^j

    a = 1.5  # = 3/2
    total = 0.0
    term_details = []

    for j in range(J + 1):
        # Coefficient: C(s+j-1, j) * (1/4)^j
        c_j = _generalized_binomial(s + j - 1, j) * (0.25 ** j)

        # Hurwitz zeta argument: 2s + 2j - 1
        hz_arg = 2.0 * s + 2.0 * j - 1.0

        # Evaluate the Hurwitz zeta
        if hz_arg < 0 or abs(hz_arg - round(hz_arg)) < 1e-12:
            # Non-positive integer or close to it: use Bernoulli formula
            # zeta_H(-n, a) = -B_{n+1}(a) / (n+1)
            n_int = -int(round(hz_arg))
            if n_int >= 0 and n_int + 1 <= 6:
                hz_val = _hurwitz_zeta_negative_int(n_int, a)
            else:
                # For zeta_H(0, a) = 1/2 - a (this is B_1(a) = a - 1/2,
                # and zeta_H(0,a) = -B_1(a) = 1/2 - a)
                # Actually _hurwitz_zeta_negative_int(0, a) handles n=0
                # which gives zeta_H(0, a) = -B_1(a)/1 = -(a - 1/2) = 1/2 - a
                # This is correct.  But if n_int+1 > 6, we cannot compute.
                hz_val = 0.0  # fallback; should not be reached for J <= 7
        elif hz_arg > 1.0:
            hz_val = _hurwitz_zeta_positive(hz_arg, a, terms=500)
        elif abs(hz_arg - 1.0) < 1e-12:
            # zeta_H(1, a) diverges (pole)
            # This should not arise for our parameter range
            hz_val = float('inf')
        else:
            # 0 < hz_arg < 1, not an integer: cannot compute with our tools
            # For the specific case s = -1/2, hz_arg = 2j - 2, always integer
            # So this branch should not be reached.
            hz_val = 0.0  # fallback

        contribution = 2.0 * c_j * hz_val
        total += contribution

        term_details.append({
            "j": j,
            "c_j": c_j,
            "hurwitz_arg": hz_arg,
            "hurwitz_value": hz_val,
            "contribution": contribution,
        })

    # Convergence diagnostic: compare last few terms
    if len(term_details) >= 2:
        last_contrib = abs(term_details[-1]["contribution"])
        second_last = abs(term_details[-2]["contribution"])
        ratio = last_contrib / second_last if second_last != 0 else 0.0
    else:
        ratio = 1.0

    return {
        "zeta_value": total,
        "s": s,
        "l_max": l_max,
        "method": "analytic_continuation",
        "partial_sum": partial_sum,
        "convergence_diagnostics": {
            "n_terms": J + 1,
            "term_details": term_details,
            "last_term_ratio": ratio,
            "converged": ratio < 0.1,
            "description": (
                f"Analytic continuation via Hurwitz zeta expansion "
                f"with {J + 1} terms.  "
                f"|last term / second-to-last| = {ratio:.4e}.  "
                f"Partial sum (non-convergent) = {partial_sum:.6f}."
            ),
        },
    }


# ---------------------------------------------------------------------------
# 3. Casimir energy on S^2
# ---------------------------------------------------------------------------

def compute_casimir_energy_s2(a_radius, field_content=None):
    """
    Compute the Casimir energy of the 6D graviton tower on S^2.

    The Casimir energy is:
        V_Cas = N_eff * zeta_reg / (4 * pi * a^4)

    where:
        N_eff = 9 (total bosonic degrees of freedom of the 6D graviton)
        zeta_reg = zeta_S2(-1/2) (regularized spectral sum)
        a = radius of the S^2

    The sign of V_Cas is determined by the sign of zeta_S2(-1/2).
    For the graviton tower on S^2, zeta_S2(-1/2) < 0, so V_Cas < 0
    (attractive Casimir force).

    Parameters
    ----------
    a_radius : float -- radius of the S^2 in Planck units
    field_content : dict or None -- if None, uses 6D graviton (N_eff = 9)

    Returns
    -------
    dict with keys:
        V_casimir : float -- Casimir energy
        coefficient : float -- C_cas = N_eff * zeta_reg / (4*pi)
        N_eff : int -- total bosonic degrees of freedom
        zeta_value : float -- zeta_S2(-1/2)
        a_radius : float
        scaling_power : int -- -4 (V ~ a^{-4})
        sign : str -- "negative" or "positive"
        field_content : dict
        description : str
    """
    # Default field content: 6D graviton
    if field_content is None:
        N_eff = 9
        field_content = {
            "source": "6D graviton",
            "total_polarizations": 9,
            "bosonic": True,
        }
    else:
        N_eff = field_content.get("total_polarizations", 9)

    # Compute the spectral zeta function at s = -1/2
    zeta_result = compute_spectral_zeta_s2(s=-0.5)
    zeta_value = zeta_result["zeta_value"]

    # Casimir coefficient
    C_cas = N_eff * zeta_value / (4.0 * math.pi)

    # Casimir energy
    V_casimir = C_cas / (a_radius ** 4)

    sign = "negative" if V_casimir < 0 else "positive"

    return {
        "V_casimir": V_casimir,
        "coefficient": C_cas,
        "N_eff": N_eff,
        "zeta_value": zeta_value,
        "a_radius": a_radius,
        "scaling_power": -4,
        "sign": sign,
        "field_content": field_content,
        "description": (
            f"Casimir energy V_Cas = C_cas / a^4 with C_cas = N_eff * zeta_reg "
            f"/ (4*pi).  N_eff = {N_eff}, zeta_S2(-1/2) = {zeta_value:.6f}, "
            f"C_cas = {C_cas:.6f}.  V_Cas = {V_casimir:.6e} (Planck units).  "
            f"Sign: {sign} -- Casimir energy is {'attractive' if sign == 'negative' else 'repulsive'}."
        ),
    }


# ---------------------------------------------------------------------------
# 4. Effective potential
# ---------------------------------------------------------------------------

def compute_effective_potential(sigma=None, constants=None):
    """
    Compute the effective 4D potential for the volume modulus sigma.

    The potential combines the classical curvature contribution and the
    1-loop Casimir energy:

        V_total(sigma) = V_curv(sigma) + V_cas(sigma)
        V_curv(sigma) = B * exp(2*sigma)
        V_cas(sigma) = A * exp(4*sigma)

    where:
        B = -chi / (2 * a_0^2)    (classical curvature, M_Pl = 1)
        A = C_cas / a_0^4         (Casimir energy coefficient)

    For S^2 (chi = 2):  B = -1/a_0^2 < 0
    For genus 2 (chi = -2):  B = 1/a_0^2 > 0

    Uses a_0 = 1.0 (Planck units) and M_Pl = 1.

    The potential is evaluated on a grid of 41 sigma values from -3 to 3.

    Parameters
    ----------
    sigma : ignored (grid is computed internally)
    constants : ignored (Planck units used throughout)

    Returns
    -------
    dict with keys:
        sigma_grid : list of float
        V_classical : list of float
        V_casimir : list of float
        V_total : list of float
        A_coefficient : float
        B_coefficient : float
        chi : int
        genus : int
        a_0 : float
        description : str
        genus2_data : dict with same structure for genus 2
    """
    a_0 = 1.0

    # Compute the Casimir coefficient
    casimir = compute_casimir_energy_s2(a_radius=a_0)
    C_cas = casimir["coefficient"]

    # --- S^2 (genus 0, chi = 2) ---
    chi_s2 = 2
    B_s2 = -chi_s2 / (2.0 * a_0 ** 2)      # = -1.0
    A_s2 = C_cas / a_0 ** 4                  # Casimir coefficient

    # Sigma grid
    n_points = 41
    sigma_min = -3.0
    sigma_max = 3.0
    step = (sigma_max - sigma_min) / (n_points - 1)
    sigma_grid = [sigma_min + i * step for i in range(n_points)]

    V_classical_s2 = []
    V_casimir_s2 = []
    V_total_s2 = []

    for sig in sigma_grid:
        v_cl = B_s2 * math.exp(2.0 * sig)
        v_cas = A_s2 * math.exp(4.0 * sig)
        V_classical_s2.append(v_cl)
        V_casimir_s2.append(v_cas)
        V_total_s2.append(v_cl + v_cas)

    # --- Genus 2 (chi = -2) ---
    chi_g2 = -2
    B_g2 = -chi_g2 / (2.0 * a_0 ** 2)       # = +1.0
    A_g2 = C_cas / a_0 ** 4                   # same Casimir coefficient

    V_classical_g2 = []
    V_casimir_g2 = []
    V_total_g2 = []

    for sig in sigma_grid:
        v_cl = B_g2 * math.exp(2.0 * sig)
        v_cas = A_g2 * math.exp(4.0 * sig)
        V_classical_g2.append(v_cl)
        V_casimir_g2.append(v_cas)
        V_total_g2.append(v_cl + v_cas)

    genus2_data = {
        "sigma_grid": sigma_grid,
        "V_classical": V_classical_g2,
        "V_casimir": V_casimir_g2,
        "V_total": V_total_g2,
        "A_coefficient": A_g2,
        "B_coefficient": B_g2,
        "chi": chi_g2,
        "genus": 2,
        "a_0": a_0,
        "description": (
            f"Genus-2 potential: B = {B_g2:.6f} (positive, chi = {chi_g2}), "
            f"A = {A_g2:.6f} (Casimir).  "
            f"Since A < 0 and B > 0, V = A*e^{{4s}} + B*e^{{2s}} has a "
            f"stationary point but it is a MAXIMUM (V'' < 0), not a minimum."
        ),
    }

    return {
        "sigma_grid": sigma_grid,
        "V_classical": V_classical_s2,
        "V_casimir": V_casimir_s2,
        "V_total": V_total_s2,
        "A_coefficient": A_s2,
        "B_coefficient": B_s2,
        "chi": chi_s2,
        "genus": 0,
        "a_0": a_0,
        "description": (
            f"S^2 potential: B = {B_s2:.6f} (negative, chi = {chi_s2}), "
            f"A = {A_s2:.6f} (Casimir).  "
            f"Both A and B are negative, so V = A*e^{{4s}} + B*e^{{2s}} is "
            f"monotonically negative and has no minimum."
        ),
        "genus2_data": genus2_data,
    }


# ---------------------------------------------------------------------------
# 5. Find Casimir minimum
# ---------------------------------------------------------------------------

def find_casimir_minimum(constants=None):
    """
    Analyse whether the Casimir-corrected potential has a stable minimum.

    For V(sigma) = A * e^{4*sigma} + B * e^{2*sigma}:
        V'(sigma) = 4A * e^{4s} + 2B * e^{2s} = 0
        => e^{2*s_0} = -B / (2A)

    For a real solution we need -B/(2A) > 0.

    The second derivative at the stationary point:
        V''(s_0) = 16A * e^{4s_0} + 4B * e^{2s_0}
                 = 16A * [B/(2A)]^2 + 4B * [-B/(2A)]
                 = 16A * B^2/(4A^2) + 4B * [-B/(2A)]
                 = 4B^2/A - 2B^2/A
                 = 2B^2/A

    For a MINIMUM we need V'' > 0, which requires A > 0.

    For pure graviton Casimir: A < 0, so:
    - S^2 (B < 0, A < 0): -B/(2A) < 0, no real stationary point
    - Genus 2 (B > 0, A < 0): -B/(2A) > 0, stationary point exists
      but V'' = 2B^2/A < 0, so it is a MAXIMUM

    Parameters
    ----------
    constants : ignored (Planck units used throughout)

    Returns
    -------
    dict with analysis for both S^2 and genus 2
    """
    a_0 = 1.0

    # Get Casimir coefficient
    casimir = compute_casimir_energy_s2(a_radius=a_0)
    C_cas = casimir["coefficient"]

    results = {}

    for label, chi, genus in [("s2", 2, 0), ("genus2", -2, 2)]:
        B = -chi / (2.0 * a_0 ** 2)
        A = C_cas / a_0 ** 4

        # Check for stationary point: need -B/(2A) > 0
        if A != 0:
            ratio = -B / (2.0 * A)
        else:
            ratio = None

        stationary_exists = ratio is not None and ratio > 0

        sigma_0 = None
        a_stabilized = None
        V_min = None
        V_double_prime = None
        is_minimum = False
        is_maximum = False

        if stationary_exists:
            sigma_0 = 0.5 * math.log(ratio)
            # a_stabilized = a_0 * exp(-sigma_0) (since beta = -1 for d=4, n=2)
            a_stabilized = a_0 * math.exp(-sigma_0)

            # V at stationary point
            e2s = math.exp(2.0 * sigma_0)
            e4s = e2s ** 2
            V_min = A * e4s + B * e2s

            # V'' = 2 B^2 / A
            V_double_prime = 2.0 * B ** 2 / A

            is_minimum = V_double_prime > 0
            is_maximum = V_double_prime < 0

        # Honest assessment
        if not stationary_exists:
            honest = (
                f"No stationary point exists for {label} "
                f"(chi = {chi}, A = {A:.6f}, B = {B:.6f}).  "
                f"-B/(2A) = {ratio:.6f} <= 0, so e^{{2*sigma_0}} has no "
                f"real positive solution.  The potential is monotonic and "
                f"Casimir stabilization does not work."
            )
        elif is_maximum:
            honest = (
                f"Stationary point exists at sigma_0 = {sigma_0:.6f} "
                f"for {label} (chi = {chi}), but it is a MAXIMUM "
                f"(V'' = 2B^2/A = {V_double_prime:.6f} < 0 because A < 0).  "
                f"Casimir stabilization does NOT produce a stable minimum.  "
                f"The volume modulus rolls away from this point."
            )
        else:
            honest = (
                f"Stable minimum exists at sigma_0 = {sigma_0:.6f} "
                f"for {label} (chi = {chi}) with V'' = {V_double_prime:.6f} > 0."
            )

        results[label] = {
            "minimum_exists": is_minimum,
            "stationary_point_exists": stationary_exists,
            "sigma_0": sigma_0,
            "a_stabilized": a_stabilized,
            "V_min": V_min,
            "V_double_prime": V_double_prime,
            "is_minimum": is_minimum,
            "is_maximum": is_maximum,
            "A": A,
            "B": B,
            "chi": chi,
            "genus": genus,
            "honest_assessment": honest,
        }

    # Primary result is S^2; genus 2 is comparison
    primary = results["s2"]

    return {
        "minimum_exists": primary["is_minimum"],
        "stationary_point_exists": primary["stationary_point_exists"],
        "sigma_0": primary["sigma_0"],
        "a_stabilized": primary["a_stabilized"],
        "V_min": primary["V_min"],
        "V_double_prime": primary["V_double_prime"],
        "is_minimum": primary["is_minimum"],
        "is_maximum": primary["is_maximum"],
        "A": primary["A"],
        "B": primary["B"],
        "chi": primary["chi"],
        "genus": primary["genus"],
        "honest_assessment": primary["honest_assessment"],
        "genus2_result": results["genus2"],
    }


# ---------------------------------------------------------------------------
# 6. Dilaton mass from Casimir stabilization
# ---------------------------------------------------------------------------

def compute_dilaton_mass_casimir(constants=None):
    """
    Compute the dilaton mass from Casimir stabilization, if a minimum exists.

    If a stable minimum exists in the Casimir potential:
        m_phi^2 = V''(sigma_0) = 2 * B^2 / A

    The mass in eV is obtained by converting from Planck units:
        m_eV = m_Planck * sqrt(m_phi^2)
    where m_Planck = 1.22089e28 eV.

    The dilaton must satisfy the experimental bound m_phi > 2e-3 eV
    (from fifth-force searches at sub-millimeter scales).

    If no minimum exists (which is the honest result for the pure graviton
    Casimir energy), this function reports the failure and explains why.

    Parameters
    ----------
    constants : ignored (Planck units used throughout)

    Returns
    -------
    dict with keys:
        m_phi_squared : float or None
        m_phi_eV : float or None
        lambda_compton_m : float or None
        exceeds_threshold : bool or None
        minimum_exists : bool
        first_principles : bool
        comparison_with_threshold : str
        honest_assessment : str
    """
    minimum_result = find_casimir_minimum()

    m_Planck_eV = 1.22089e28  # eV
    hbar_c_eV_m = 1.9733e-7   # hbar * c in eV * m
    threshold_eV = 2e-3        # 2 meV

    m_phi_squared = None
    m_phi_eV = None
    lambda_compton_m = None
    exceeds_threshold = None
    first_principles = False

    if minimum_result["minimum_exists"]:
        # Stable minimum exists
        m_phi_squared = minimum_result["V_double_prime"]
        if m_phi_squared > 0:
            m_phi_planck = math.sqrt(m_phi_squared)
            m_phi_eV = m_Planck_eV * m_phi_planck
            lambda_compton_m = hbar_c_eV_m / m_phi_eV
            exceeds_threshold = m_phi_eV > threshold_eV
            first_principles = True
            comparison = (
                f"m_phi = {m_phi_eV:.4e} eV.  "
                f"Threshold: {threshold_eV:.0e} eV.  "
                f"{'Exceeds' if exceeds_threshold else 'Below'} threshold."
            )
        else:
            comparison = (
                "m_phi^2 <= 0: tachyonic or massless.  "
                "Casimir stabilization produces an unstable direction."
            )
    else:
        comparison = (
            "No stable minimum exists in the Casimir potential.  "
            "The dilaton mass cannot be computed from Casimir stabilization "
            "alone because the pure graviton tower on S^2 produces a "
            "negative Casimir coefficient A < 0, which means V''(s_0) < 0 "
            "(maximum, not minimum) or no stationary point at all."
        )

    # Honest assessment
    honest_parts = []

    if not minimum_result["minimum_exists"]:
        honest_parts.append(
            "The Casimir energy from the pure 6D graviton KK tower on S^2 "
            "does NOT stabilize the volume modulus."
        )
        honest_parts.append(
            "The spectral zeta function zeta_S2(-1/2) is NEGATIVE, making "
            "the Casimir coefficient A negative."
        )
        honest_parts.append(
            "For S^2 (chi = 2): both A and B are negative, so the potential "
            "is monotonically decreasing -- no stationary point."
        )
        honest_parts.append(
            "For genus 2 (chi = -2): B > 0, A < 0, so a stationary point "
            "exists but V'' = 2B^2/A < 0 -- it is a MAXIMUM."
        )
        honest_parts.append(
            "Possible resolutions: (1) adding fermion contributions (opposite "
            "Casimir sign could flip A to positive); (2) flux contributions; "
            "(3) higher-loop corrections; (4) different regularization scheme."
        )
    else:
        honest_parts.append(
            f"Casimir stabilization produces a minimum with "
            f"m_phi = {m_phi_eV:.4e} eV."
        )

    return {
        "m_phi_squared": m_phi_squared,
        "m_phi_eV": m_phi_eV,
        "lambda_compton_m": lambda_compton_m,
        "exceeds_threshold": exceeds_threshold,
        "minimum_exists": minimum_result["minimum_exists"],
        "first_principles": first_principles,
        "comparison_with_threshold": comparison,
        "honest_assessment": "  ".join(honest_parts),
    }


# ---------------------------------------------------------------------------
# 7a. Fermion spectral zeta on S^2
# ---------------------------------------------------------------------------

def compute_fermion_spectral_zeta_s2(s, l_max=1000):
    """
    Compute the fermion spectral zeta function on S^2.

    Fermion eigenvalues on S^2 are (l + 1/2)^2 for l >= 0, with
    degeneracy 2(2l + 1).  The spectral zeta function is:

        zeta_F(s) = sum_{l >= 0} 2(2l+1) * [(l+1/2)^2]^{-s}
                  = 2 * sum_{l >= 0} (2l+1) * (l+1/2)^{-2s}
                  = 4 * sum_{l >= 0} (l+1/2)^{1-2s}
                  = 4 * zeta_H(2s-1, 1/2)

    For s = -1/2:
        zeta_F(-1/2) = 4 * zeta_H(-2, 1/2) = 4 * (-B_3(1/2) / 3)

    B_3(1/2) = (1/2)^3 - 1.5*(1/2)^2 + 0.5*(1/2)
             = 0.125 - 0.375 + 0.25 = 0   (exactly)

    So zeta_F(-1/2) = 0 for fermions on S^2.

    Parameters
    ----------
    s : float -- the exponent
    l_max : int -- maximum l for partial sums (default 1000)

    Returns
    -------
    dict with keys:
        zeta_value : float
        s : float
        method : str
        description : str
    """
    if s > 1.0:
        # Direct partial sum for convergent regime
        total = 0.0
        for l in range(0, l_max + 1):
            eigenvalue = (l + 0.5) ** 2
            degeneracy = 2 * (2 * l + 1)
            total += degeneracy * eigenvalue ** (-s)
        return {
            "zeta_value": total,
            "s": s,
            "method": "partial_sum",
            "description": (
                f"Fermion spectral zeta on S^2 via partial sum with "
                f"{l_max + 1} terms.  zeta_F({s}) = {total:.10f}."
            ),
        }

    # Analytic continuation: zeta_F(s) = 4 * zeta_H(2s - 1, 1/2)
    hz_arg = 2.0 * s - 1.0
    n_int = -int(round(hz_arg))

    if n_int >= 0 and abs(hz_arg + n_int) < 1e-12 and n_int + 1 <= 6:
        # Use Bernoulli polynomial identity: zeta_H(-n, a) = -B_{n+1}(a)/(n+1)
        hz_val = _hurwitz_zeta_negative_int(n_int, 0.5)
        zeta_value = 4.0 * hz_val
        method = "analytic_continuation"
    elif hz_arg > 1.0:
        hz_val = _hurwitz_zeta_positive(hz_arg, 0.5, terms=l_max)
        zeta_value = 4.0 * hz_val
        method = "partial_sum_hurwitz"
    else:
        # Fallback: direct partial sum (non-convergent, for diagnostics)
        total = 0.0
        for l in range(0, l_max + 1):
            eigenvalue = (l + 0.5) ** 2
            degeneracy = 2 * (2 * l + 1)
            total += degeneracy * eigenvalue ** (-s)
        zeta_value = total
        method = "partial_sum_fallback"

    return {
        "zeta_value": zeta_value,
        "s": s,
        "method": method,
        "description": (
            f"Fermion spectral zeta on S^2 via {method}.  "
            f"zeta_F({s}) = {zeta_value:.10f}.  "
            f"For s=-1/2, B_3(1/2) = 0 exactly, so zeta_F(-1/2) = 0."
        ),
    }


# ---------------------------------------------------------------------------
# 7b. Vector boson spectral zeta on S^2
# ---------------------------------------------------------------------------

def compute_vector_spectral_zeta_s2(s, l_max=1000):
    """
    Compute the vector boson spectral zeta function on S^2.

    Vector eigenvalues on S^2 are l(l+1) - 1 = (l-1)(l+2) for l >= 1,
    with degeneracy 2(2l+1).  The l=1 mode gives eigenvalue 0 (gauge
    zero mode) and is excluded from the massive spectrum.

    For Re(s) > 1, compute via direct partial sum from l = 2 to l_max.

    For s = -1/2, use analytic continuation via the substitution
    u = l + 1/2, giving eigenvalue u^2 - 9/4.  This is expressed in
    terms of Hurwitz zeta functions.

    Parameters
    ----------
    s : float -- the exponent
    l_max : int -- maximum l for partial sums (default 1000)

    Returns
    -------
    dict with keys:
        zeta_value : float
        s : float
        l_max : int
        method : str
        description : str
    """
    if s > 1.0:
        # Direct partial sum (convergent for s > 1)
        total = 0.0
        for l in range(2, l_max + 1):
            eigenvalue = (l - 1) * (l + 2)  # = l(l+1) - 2  ... actually l^2+l-2
            degeneracy = 2 * (2 * l + 1)
            total += degeneracy * eigenvalue ** (-s)
        return {
            "zeta_value": total,
            "s": s,
            "l_max": l_max,
            "method": "partial_sum",
            "description": (
                f"Vector spectral zeta on S^2 via partial sum (l=2..{l_max}).  "
                f"zeta_V({s}) = {total:.10f}."
            ),
        }

    # Analytic continuation for s <= 1
    # Eigenvalue = (l-1)(l+2) = l^2 + l - 2.  Substituting u = l + 1/2:
    #   l^2 + l - 2 = u^2 - 1/4 - 2 = u^2 - 9/4
    # Degeneracy 2(2l+1) = 4u.
    #
    # zeta_V(s) = sum_{l>=2} 4u * (u^2 - 9/4)^{-s},  u = l+1/2 >= 5/2
    #
    # For s = -1/2:
    #   zeta_V(-1/2) = sum_{l>=2} 4u * (u^2 - 9/4)^{1/2}
    #
    # We use binomial expansion: (u^2 - 9/4)^{-s} = u^{-2s} * (1 - 9/(4u^2))^{-s}
    #   = u^{-2s} * sum_{j>=0} C(-s, j) * (-9/4)^j * u^{-2j}
    #   = sum_{j>=0} C(-s, j) * (-9/4)^j * u^{-2s-2j}
    #
    # Then: zeta_V(s) = 4 * sum_j C(-s,j) * (-9/4)^j * sum_{l>=2} u^{1-2s-2j}
    #                  = 4 * sum_j C(-s,j) * (-9/4)^j * zeta_H(2s+2j-1, 5/2)
    #                                                   (shifted to start at u=5/2)
    # Actually: sum_{l>=2} u^{-alpha} = sum_{k>=0} (k + 5/2)^{-alpha} = zeta_H(alpha, 5/2)

    J = 7  # 8 terms in the binomial expansion
    a = 2.5  # = 5/2, since u starts at l=2 -> u = 2.5
    total = 0.0

    for j in range(J + 1):
        # Coefficient: C(-s, j) * (-9/4)^j
        c_j = _generalized_binomial(-s, j) * ((-9.0 / 4.0) ** j)

        # Hurwitz zeta argument: 2s + 2j - 1
        hz_arg = 2.0 * s + 2.0 * j - 1.0

        if hz_arg < 0 or abs(hz_arg - round(hz_arg)) < 1e-12:
            n_int = -int(round(hz_arg))
            if n_int >= 0 and n_int + 1 <= 6:
                hz_val = _hurwitz_zeta_negative_int(n_int, a)
            else:
                hz_val = 0.0
        elif hz_arg > 1.0:
            hz_val = _hurwitz_zeta_positive(hz_arg, a, terms=500)
        elif abs(hz_arg - 1.0) < 1e-12:
            hz_val = float('inf')
        else:
            hz_val = 0.0

        total += 4.0 * c_j * hz_val

    return {
        "zeta_value": total,
        "s": s,
        "l_max": l_max,
        "method": "analytic_continuation",
        "description": (
            f"Vector spectral zeta on S^2 via analytic continuation "
            f"(Hurwitz zeta expansion with {J + 1} terms).  "
            f"zeta_V({s}) = {total:.10f}.  "
            f"l=1 gauge zero mode excluded; sum starts from l=2."
        ),
    }


# ---------------------------------------------------------------------------
# 7c. Matter Casimir coefficient
# ---------------------------------------------------------------------------

def compute_matter_casimir_coefficient(n_scalars=0, n_fermions=0,
                                       n_vectors=0, a_radius=1.0):
    """
    Compute the total Casimir coefficient combining all spin species.

    The Casimir energy on S^2 receives contributions from different
    spin fields with different signs:
      - Bosonic scalars: +zeta_scalar(-1/2) per degree of freedom
      - Fermions: -zeta_fermion(-1/2) per degree of freedom
      - Vectors: +zeta_vector(-1/2) per degree of freedom

    The graviton tower contributes 9 bosonic scalar-type degrees of
    freedom (the existing pure graviton computation).

    For 6D SUSY multiplets on S^2:
      - Each hypermultiplet: 4 real scalars + 2 fermion dof
      - Each vector multiplet: 1 vector + 2 fermion dof (gaugino)

    The caller provides the raw counts (n_scalars, n_fermions, n_vectors)
    for the matter sector.  The function adds the graviton contribution
    (9 scalar dof) automatically.

    Total:
        A_total = [N_grav * zeta_S + n_scalars * zeta_S
                   - n_fermions * zeta_F + n_vectors * zeta_V]
                  / (4 * pi * a^4)

    Parameters
    ----------
    n_scalars : int -- number of matter scalar dof
    n_fermions : int -- number of matter fermion dof
    n_vectors : int -- number of matter vector dof
    a_radius : float -- radius of S^2 (default 1.0)

    Returns
    -------
    dict with keys:
        A_total : float
        A_graviton : float
        A_matter : float
        sign_flipped : bool
        n_scalars : int
        n_fermions : int
        n_vectors : int
        description : str
    """
    N_grav = 9  # graviton tower scalar-type dof

    zeta_S = compute_spectral_zeta_s2(s=-0.5)["zeta_value"]
    zeta_F = compute_fermion_spectral_zeta_s2(s=-0.5)["zeta_value"]
    zeta_V = compute_vector_spectral_zeta_s2(s=-0.5)["zeta_value"]

    prefactor = 1.0 / (4.0 * math.pi * a_radius ** 4)

    A_graviton = N_grav * zeta_S * prefactor
    A_matter = (n_scalars * zeta_S - n_fermions * zeta_F
                + n_vectors * zeta_V) * prefactor
    A_total = A_graviton + A_matter

    # Sign flip: graviton alone gives A < 0; check if matter flips it
    sign_flipped = (A_total > 0) and (A_graviton < 0)

    return {
        "A_total": A_total,
        "A_graviton": A_graviton,
        "A_matter": A_matter,
        "sign_flipped": sign_flipped,
        "n_scalars": n_scalars,
        "n_fermions": n_fermions,
        "n_vectors": n_vectors,
        "description": (
            f"Total Casimir coefficient A = {A_total:.8e}.  "
            f"Graviton (9 scalar dof): A_grav = {A_graviton:.8e}.  "
            f"Matter ({n_scalars}S + {n_fermions}F + {n_vectors}V): "
            f"A_matter = {A_matter:.8e}.  "
            f"Sign flipped by matter: {sign_flipped}."
        ),
    }


# ---------------------------------------------------------------------------
# 7d. Scan anomaly-free groups for matter Casimir
# ---------------------------------------------------------------------------

def scan_anomaly_free_matter_casimir():
    """
    Loop over anomaly-free gauge groups and compute the matter Casimir
    coefficient for each.

    For each group in _ANOMALY_FREE_GROUPS that has both n_hyper and
    n_vector defined (skip groups with n_hyper = None), compute the
    matter Casimir coefficient using the 6D multiplet field counting:
      - n_scalars_total = 4 * n_hyper
      - n_fermions_total = 2 * n_hyper + 2 * n_vector
      - n_vectors_total = n_vector

    Returns
    -------
    dict with keys:
        results : list of dict -- per-group results
        any_sign_flip : bool -- whether any group flips the Casimir sign
        description : str
    """
    results = []
    any_sign_flip = False

    for key, data in _ANOMALY_FREE_GROUPS.items():
        n_hyper = data.get("n_hyper")
        n_vector = data.get("n_vector")

        if n_hyper is None or n_vector is None:
            continue

        # 6D multiplet field counting on S^2
        n_scalars_total = 4 * n_hyper
        n_fermions_total = 2 * n_hyper + 2 * n_vector
        n_vectors_total = n_vector

        casimir = compute_matter_casimir_coefficient(
            n_scalars=n_scalars_total,
            n_fermions=n_fermions_total,
            n_vectors=n_vectors_total,
        )

        entry = {
            "group": data["group"],
            "key": key,
            "n_hyper": n_hyper,
            "n_vector": n_vector,
            "n_scalars_total": n_scalars_total,
            "n_fermions_total": n_fermions_total,
            "n_vectors_total": n_vectors_total,
            "casimir": casimir,
        }
        results.append(entry)

        if casimir["sign_flipped"]:
            any_sign_flip = True

    return {
        "results": results,
        "any_sign_flip": any_sign_flip,
        "description": (
            f"Scanned {len(results)} anomaly-free groups for matter Casimir.  "
            f"Sign flip found: {any_sign_flip}."
        ),
    }


# ---------------------------------------------------------------------------
# 7. Summary / dashboard entry point
# ---------------------------------------------------------------------------

def summarize_casimir_stabilization(constants=None):
    """
    Run the full Casimir stabilization pipeline and return a summary.

    This is the main entry point for the Streamlit dashboard.  It calls
    all sub-functions and assembles the results into a single dict with
    an overall assessment and gap status.

    The pipeline:
        1. Compute KK spectrum on S^2
        2. Compute spectral zeta function zeta_S2(-1/2)
        3. Compute Casimir energy
        4. Compute effective potential on sigma grid
        5. Search for stable minimum
        6. Compute dilaton mass (if minimum exists)

    Parameters
    ----------
    constants : ignored (Planck units used throughout)

    Returns
    -------
    dict with all sub-results and overall assessment
    """
    # Run the pipeline
    spectrum = compute_kk_spectrum_s2(l_max=100)
    zeta = compute_spectral_zeta_s2(s=-0.5)
    casimir = compute_casimir_energy_s2(a_radius=1.0)
    potential = compute_effective_potential()
    minimum = find_casimir_minimum()
    mass = compute_dilaton_mass_casimir()

    # Matter loop Casimir corrections
    matter_loop_results = None
    if _ANOMALY_AVAILABLE:
        try:
            matter_loop_results = scan_anomaly_free_matter_casimir()
        except Exception:
            pass

    # Overall assessment
    zeta_negative = zeta["zeta_value"] < 0
    has_minimum = minimum["minimum_exists"]

    if has_minimum and mass["first_principles"]:
        overall = (
            "Casimir stabilization succeeds: the 1-loop Casimir energy "
            "from the KK graviton tower produces a stable minimum for the "
            "volume modulus, yielding a dilaton mass from first principles."
        )
    else:
        overall = (
            "Casimir stabilization from the pure graviton tower FAILS to "
            "produce a stable minimum.  The spectral zeta function "
            f"zeta_S2(-1/2) = {zeta['zeta_value']:.6f} is negative, giving "
            f"a negative Casimir coefficient.  This means the 1-loop potential "
            f"does not have the right shape to trap the volume modulus.  "
            f"The dilaton mass cannot be derived from Casimir stabilization "
            f"alone.  This is an honest result: the graviton tower by itself "
            f"is insufficient, and additional ingredients (fermions, fluxes, "
            f"or non-perturbative effects) are needed."
        )

    # Gap status
    # Gap 1: vacuum polynomial (handled elsewhere)
    # Gap 3: dilaton mass from first principles
    gap1_resolved = False  # not addressed by this module
    gap3_resolved = has_minimum and mass["first_principles"]

    gaps_status = {
        "gap1_resolved": gap1_resolved,
        "gap3_resolved": gap3_resolved,
        "gap3_detail": (
            "Resolved: dilaton mass computed from Casimir stabilization."
            if gap3_resolved else
            "NOT resolved: Casimir stabilization alone does not produce a "
            "stable minimum.  The dilaton mass requires additional physics "
            "beyond the pure graviton Casimir energy."
        ),
    }

    # What IS computed from first principles
    first_principles_results = {
        "kk_spectrum": True,
        "spectral_zeta_function": True,
        "casimir_coefficient_sign": True,
        "effective_potential_shape": True,
        "stable_minimum": False,
        "dilaton_mass": False,
        "description": (
            "From first principles, we compute: the KK spectrum, the "
            "spectral zeta function (via analytic continuation), the sign "
            "of the Casimir coefficient (negative), and the shape of the "
            "effective potential.  What we CANNOT derive is a stable minimum "
            "or a dilaton mass, because the coefficient has the wrong sign."
        ),
    }

    return {
        "spectrum": spectrum,
        "zeta": zeta,
        "casimir": casimir,
        "potential": potential,
        "minimum": minimum,
        "mass": mass,
        "overall_assessment": overall,
        "gaps_status": gaps_status,
        "first_principles_results": first_principles_results,
        "zeta_s2_value": zeta["zeta_value"],
        "zeta_s2_negative": zeta_negative,
        "casimir_coefficient": casimir["coefficient"],
        "has_stable_minimum": has_minimum,
        "matter_loop_results": matter_loop_results,
    }
