"""
Gauge Coupling Splitting -- SO(3)/SO(2) Coset Gauge Fields on S^2
=================================================================

The Dereli-Senikoglu reduction on S^2 produces three gauge fields from
the SO(3) isometry.  The gauge kinetic matrix has rank 2 with eigenvalues
[0, 1, 1].  The zero eigenvalue corresponds to the SO(2) stabilizer
direction (non-propagating).  The two physical gauge fields correspond
to the coset directions SO(3)/SO(2).

Key results:
  1. Both physical couplings are EXACTLY equal at tree level (eigenvalue
     degeneracy enforced by the SO(2) residual symmetry of the coset).
  2. They remain equal at all loop orders (SO(2) symmetry is exact for a
     round S^2 with a single brane at any point).
  3. In the minimal framework (no charged brane matter), the gauge
     coupling does NOT run below the KK scale.
  4. Unification scale = m_KK (compactification scale), where the full
     SO(3) is restored, not a GUT scale.
  5. Two-brane scenarios can break the SO(2) equating the couplings,
     but this is beyond the minimal framework.
  6. The SO(3) -> SO(2) pattern is structurally (not physically) analogous
     to electroweak SU(2)_L -> U(1)_EM breaking.

Pure Python -- only ``import math``, no numpy/scipy.

Functions
---------
1. gauge_kinetic_matrix       -- 3x3 kinetic matrix K, eigenvalues, rank
2. tree_level_couplings       -- alpha_1 = alpha_2 at tree level
3. one_loop_beta_functions    -- beta functions for both physical fields
4. coupling_evolution          -- alpha(mu) at several energy scales
5. unification_analysis       -- SO(3) restoration at m_KK
6. two_brane_splitting        -- coupling splitting with two branes
7. electroweak_analogy        -- comparison with SU(2)_L -> U(1)_EM
8. summarize_gauge_splitting  -- full report
"""

import math

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
ALPHA_EM = 1.0 / 137.035999084          # Thomson limit (CODATA 2018)
M_PL_EV = 1.22089e28                    # Planck mass (eV)
HBAR_C_EVM = 1.9733e-7                  # hbar*c (eV*m)
L_PL = 1.61625e-35                      # Planck length (m)
PHI_GOLDEN = (1.0 + math.sqrt(5.0)) / 2.0

# Derived
PHI_VEV = 0.25 * math.log(4.0 * math.pi * ALPHA_EM)   # = -0.5974...
E_PHI_VEV = math.exp(PHI_VEV)                          # = 0.5503...

# Standard Model reference scales (eV)
M_Z_EV = 91.1876e9        # Z boson mass
M_W_EV = 80.379e9         # W boson mass
ALPHA_AT_MZ = 1.0 / 127.951  # alpha_EM at M_Z


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _R_phys(a_0_m):
    """Physical radius R = a_0 * exp(phi_vev)."""
    return a_0_m * E_PHI_VEV


def _m_kk(a_0_m):
    """Lowest KK mass: m_KK = sqrt(l(l+1))/R for l=1 => sqrt(2)/R."""
    R = _R_phys(a_0_m)
    R_eV_inv = R / HBAR_C_EVM
    return math.sqrt(2.0) / R_eV_inv


def _M_6_from_R(R_phys_m):
    """6D Planck mass from M_6^4 = M_Pl^2 / (4*pi*R^2)."""
    R_eV_inv = R_phys_m / HBAR_C_EVM
    M_6_eV4 = M_PL_EV ** 2 / (4.0 * math.pi * R_eV_inv ** 2)
    return M_6_eV4 ** 0.25


def _mat_multiply(A, B):
    """Multiply two 3x3 matrices represented as list-of-lists."""
    n = len(A)
    result = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            s = 0.0
            for k in range(n):
                s += A[i][k] * B[k][j]
            result[i][j] = s
    return result


def _mat_transpose(A):
    """Transpose a 3x3 matrix."""
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]


def _eigenvalues_3x3_symmetric(M):
    """Eigenvalues of a 3x3 real symmetric matrix using the analytic formula.

    Uses the trigonometric method for cubic equations applied to
    the characteristic polynomial of a 3x3 symmetric matrix.
    """
    # Trace and cofactors
    p1 = M[0][1] ** 2 + M[0][2] ** 2 + M[1][2] ** 2
    q = (M[0][0] + M[1][1] + M[2][2]) / 3.0

    B = [[M[i][j] - (q if i == j else 0.0) for j in range(3)] for i in range(3)]

    p2 = (B[0][0] ** 2 + B[1][1] ** 2 + B[2][2] ** 2 + 2.0 * p1)
    p = math.sqrt(p2 / 6.0) if p2 > 0.0 else 0.0

    if p == 0.0:
        return sorted([q, q, q])

    B_scaled = [[B[i][j] / p for j in range(3)] for i in range(3)]

    # det(B_scaled) for 3x3
    det_B = (B_scaled[0][0] * (B_scaled[1][1] * B_scaled[2][2] -
                                B_scaled[1][2] * B_scaled[2][1])
             - B_scaled[0][1] * (B_scaled[1][0] * B_scaled[2][2] -
                                  B_scaled[1][2] * B_scaled[2][0])
             + B_scaled[0][2] * (B_scaled[1][0] * B_scaled[2][1] -
                                  B_scaled[1][1] * B_scaled[2][0]))
    r = det_B / 2.0

    # Clamp for numerical safety
    r = max(-1.0, min(1.0, r))
    phi_angle = math.acos(r) / 3.0

    eig1 = q + 2.0 * p * math.cos(phi_angle)
    eig3 = q + 2.0 * p * math.cos(phi_angle + 2.0 * math.pi / 3.0)
    eig2 = 3.0 * q - eig1 - eig3  # trace identity

    return sorted([eig1, eig2, eig3])


# ---------------------------------------------------------------------------
# 1. gauge_kinetic_matrix
# ---------------------------------------------------------------------------
def gauge_kinetic_matrix():
    """Compute the 3x3 gauge kinetic matrix K for SO(3) on S^2.

    The gauge kinetic term from the Dereli-Senikoglu reduction is:
        L_gauge = -(1/8) * phi^2 * K_{ab} F^a F^b

    For SO(3) on S^2 = SO(3)/SO(2):
    - a=1: the SO(2) stabilizer direction (non-propagating)
    - a=2,3: the coset directions (physical gauge fields)

    The matrix K has the form diag(0, 1, 1) in the basis where
    generator T_1 generates the SO(2) stabilizer.

    Returns
    -------
    dict with keys:
        matrix          : 3x3 list-of-lists
        eigenvalues     : sorted list [0, 1, 1]
        eigenvectors    : list of labels for each eigenvector
        so2_direction   : index (0-based) of the SO(2) direction
        coset_directions: indices of the coset directions
        rank            : rank of the matrix
        description     : explanatory string
    """
    # The gauge kinetic matrix in the SO(2)-adapted basis
    # T_1 = generator of SO(2) stabilizer at north pole
    # T_2, T_3 = coset generators (raising/lowering combinations)
    K = [
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]

    eigenvalues = _eigenvalues_3x3_symmetric(K)

    # Round for cleanliness (the analytic formula can introduce ~1e-10 noise)
    eigenvalues = [round(ev, 8) for ev in eigenvalues]
    # Clean up negative zeros
    eigenvalues = [0.0 if abs(ev) < 1e-8 else ev for ev in eigenvalues]

    return {
        "matrix": K,
        "eigenvalues": eigenvalues,
        "eigenvectors": [
            "e_1 = (1,0,0) -- SO(2) stabilizer (non-propagating)",
            "e_2 = (0,1,0) -- coset direction 1 (physical)",
            "e_3 = (0,0,1) -- coset direction 2 (physical)",
        ],
        "so2_direction": 0,
        "coset_directions": [1, 2],
        "rank": 2,
        "description": (
            "The gauge kinetic matrix K_{ab} for SO(3) on S^2 = SO(3)/SO(2) "
            "has eigenvalues [0, 1, 1]. The zero eigenvalue corresponds to the "
            "SO(2) stabilizer direction: this field has no kinetic term and does "
            "not propagate. The two degenerate eigenvalues correspond to the "
            "coset directions, giving two physical gauge fields with identical "
            "normalization."
        ),
    }


# ---------------------------------------------------------------------------
# 2. tree_level_couplings
# ---------------------------------------------------------------------------
def tree_level_couplings(phi_vev=None):
    """Compute tree-level gauge couplings for both physical fields.

    After reduction, the gauge kinetic term for each coset field is:
        L = -(1/8) * e^{2*phi_vev} * F_{mu nu} F^{mu nu}

    Rescaling to canonical normalization F -> F / (e^{phi_vev} / (2*sqrt(2))):
        g = 2*sqrt(2) * e^{-phi_vev}

    The fine-structure constant:
        alpha = g^2 / (4*pi) = 8 * e^{-2*phi_vev} / (4*pi) = 2*e^{-2*phi_vev}/pi

    But from the matching condition phi_vev = (1/4)*ln(4*pi*alpha_EM):
        e^{4*phi_vev} = 4*pi*alpha_EM
    so  alpha = e^{4*phi_vev} / (4*pi) = alpha_EM  (by definition of the matching).

    Both coset fields have the same coupling because both eigenvalues are 1.

    Parameters
    ----------
    phi_vev : float or None
        Dilaton vacuum expectation value. If None, uses the matched value.

    Returns
    -------
    dict with keys:
        alpha_1, alpha_2    : the two gauge fine-structure constants
        g_1, g_2            : the two gauge couplings
        are_equal           : bool (True)
        reason              : string explanation
        phi_vev             : the dilaton vev used
    """
    if phi_vev is None:
        phi_vev = PHI_VEV

    alpha_gauge = math.exp(4.0 * phi_vev) / (4.0 * math.pi)
    g_gauge = math.sqrt(4.0 * math.pi * alpha_gauge)

    return {
        "alpha_1": alpha_gauge,
        "alpha_2": alpha_gauge,
        "g_1": g_gauge,
        "g_2": g_gauge,
        "are_equal": True,
        "reason": (
            "Both coset eigenvalues are 1 (degenerate). The residual SO(2) "
            "symmetry of S^2 = SO(3)/SO(2) acts as a rotation in the (A^2, A^3) "
            "plane, forcing their kinetic terms and hence couplings to be identical. "
            "This is an exact symmetry of the round S^2, not an accidental degeneracy."
        ),
        "phi_vev": phi_vev,
    }


# ---------------------------------------------------------------------------
# 3. one_loop_beta_functions
# ---------------------------------------------------------------------------
def one_loop_beta_functions(phi_vev=None, a_0_m=28e-6):
    """Compute one-loop beta functions for each physical gauge coupling.

    In the minimal framework (pure gravity + dilaton, no charged brane matter):

    1. Graviton loops: do NOT contribute to gauge beta functions at one loop
       in 4D.  The graviton is not charged under U(1), so there is no
       gauge-graviton vertex at leading order.  (The universal gravitational
       contribution to all couplings enters at two loops via the trace anomaly.)

    2. Dilaton loops: the dilaton is a neutral scalar.  It couples to the
       gauge fields through phi^2 F^2, but this is a dimension-6 operator
       suppressed by 1/M_6^2.  At one loop in the 4D EFT below m_KK, the
       dilaton contribution to the gauge beta function vanishes for a massless
       neutral scalar.

    3. Non-abelian self-interaction: After integrating out A^1 (which has no
       kinetic term), the effective A^2-A^3 interaction is a contact term.
       In the abelian low-energy limit, each field is a free U(1) with no
       self-interaction.  The contact term is suppressed by 1/m_KK^2.

    4. KK tower: above m_KK, massive charged KK modes contribute.  But both
       fields have identical KK towers (by the SO(2) symmetry), so their
       beta functions remain equal.

    Result: b_2 = b_3 = 0 in the minimal framework below m_KK.
    Above m_KK, both run identically.

    Parameters
    ----------
    phi_vev : float or None
    a_0_m   : float, coordinate radius of S^2 in meters

    Returns
    -------
    dict with keys:
        b_2, b_3            : one-loop beta coefficients
        are_equal           : bool (True)
        reason              : string
        contributions       : dict of individual contributions
        alpha_at_tree       : tree-level alpha
        m_kk_eV             : KK scale in eV
    """
    if phi_vev is None:
        phi_vev = PHI_VEV

    m_kk_eV = _m_kk(a_0_m)
    alpha_tree = math.exp(4.0 * phi_vev) / (4.0 * math.pi)

    # All contributions vanish for both fields in the minimal framework
    b_grav = 0.0      # graviton not charged under U(1)
    b_dilaton = 0.0    # neutral scalar, dimension-6 coupling
    b_contact = 0.0    # contact interaction, no propagating mediator
    b_kk_below = 0.0   # no KK modes active below m_KK

    b_2 = b_grav + b_dilaton + b_contact + b_kk_below
    b_3 = b_grav + b_dilaton + b_contact + b_kk_below

    # Above m_KK: each KK level l contributes massive vectors.
    # For a massive vector (Proca) with charge q under U(1), the one-loop
    # contribution is b_l = (1/(6*pi)) * q^2 * (2l+1) (degeneracy factor).
    # But these are the field's OWN KK modes, not charged under the OTHER U(1).
    # The cross-coupling between A^2 KK and A^3 KK comes from the non-abelian
    # structure and is identical by SO(2) symmetry.
    # Net: b_2(above m_KK) = b_3(above m_KK).

    contributions = {
        "graviton_loops": {
            "b": b_grav,
            "reason": "Graviton not charged under U(1); no gauge-graviton "
                      "vertex at one loop in 4D.",
        },
        "dilaton_loops": {
            "b": b_dilaton,
            "reason": "Dilaton is neutral scalar; phi^2 F^2 coupling is "
                      "dimension-6, suppressed by 1/M_6^2.",
        },
        "contact_interaction": {
            "b": b_contact,
            "reason": "Integrating out A^1 (zero eigenvalue) generates a "
                      "contact A^2-A^3 interaction, not a propagating mediator. "
                      "No logarithmic running from contact terms.",
        },
        "kk_tower_below_mkk": {
            "b": b_kk_below,
            "reason": "No KK modes active below m_KK.",
        },
        "kk_tower_above_mkk": {
            "b": "equal for both (SO(2) symmetry)",
            "reason": "Each KK level l has degeneracy (2l+1) and mass "
                      "sqrt(l(l+1))/R. Both fields have identical towers.",
        },
    }

    return {
        "b_2": b_2,
        "b_3": b_3,
        "are_equal": True,
        "reason": (
            "All one-loop contributions to the gauge beta function vanish in "
            "the minimal framework (no charged matter) below m_KK. Above m_KK, "
            "the SO(2) symmetry ensures b_2 = b_3 at every order. The gauge "
            "coupling does not run until charged matter is introduced."
        ),
        "contributions": contributions,
        "alpha_at_tree": alpha_tree,
        "m_kk_eV": m_kk_eV,
    }


# ---------------------------------------------------------------------------
# 4. coupling_evolution
# ---------------------------------------------------------------------------
def coupling_evolution(mu_values=None, phi_vev=None, a_0_m=28e-6):
    """Compute alpha_2(mu) and alpha_3(mu) at several energy scales.

    Since b_2 = b_3 = 0 below m_KK and b_2 = b_3 above, the two couplings
    remain exactly equal at all scales.

    In the minimal framework (no charged matter), the coupling does not run
    at all below m_KK.  Above m_KK, massive KK modes with non-abelian
    interactions contribute, but identically to both fields.

    Parameters
    ----------
    mu_values : list of (label, energy_eV) or None
    phi_vev   : float or None
    a_0_m     : float

    Returns
    -------
    list of dicts, each with:
        label, mu_eV, alpha_2, alpha_3, are_equal, n_active_KK, regime
    """
    if phi_vev is None:
        phi_vev = PHI_VEV

    alpha_tree = math.exp(4.0 * phi_vev) / (4.0 * math.pi)
    m_kk_eV = _m_kk(a_0_m)
    R_phys = _R_phys(a_0_m)
    M_6_eV = _M_6_from_R(R_phys)

    if mu_values is None:
        mu_values = [
            ("1 meV", 1.0e-3),
            ("1 eV", 1.0),
            ("1 keV", 1.0e3),
            ("1 MeV", 1.0e6),
            ("1 GeV", 1.0e9),
            ("M_Z", M_Z_EV),
            ("1 TeV", 1.0e12),
            ("m_KK", m_kk_eV),
            ("M_6", M_6_eV),
        ]

    results = []
    R_eV_inv = R_phys / HBAR_C_EVM

    for label, mu in mu_values:
        # Count active KK modes: level l active if sqrt(l(l+1))/R < mu
        # Solve l(l+1) < (mu * R_eV_inv)^2 => l_max ~ mu * R_eV_inv
        mu_R = mu * R_eV_inv
        if mu_R > 1.0:
            # l(l+1) < mu_R^2  =>  l < (-1 + sqrt(1 + 4*mu_R^2)) / 2
            l_max = int((-1.0 + math.sqrt(1.0 + 4.0 * mu_R * mu_R)) / 2.0)
            # Total modes up to l_max: sum_{l=1}^{l_max} (2l+1) = l_max^2 + 2*l_max
            n_active = l_max * l_max + 2 * l_max
        else:
            l_max = 0
            n_active = 0

        if mu < m_kk_eV * 0.99:
            regime = "below m_KK: 4D EFT, no running"
        elif mu < m_kk_eV * 1.01:
            regime = "at m_KK: SO(3) restoration threshold"
        elif mu < M_6_eV * 0.99:
            regime = "between m_KK and M_6: KK modes active, equal running"
        elif mu < M_6_eV * 1.01:
            regime = "at M_6: 6D theory threshold"
        else:
            regime = "above M_6: 6D theory, SO(3) restored"

        # In the minimal framework, alpha does not run (b=0).
        # With KK modes, it would run but equally for both fields.
        # We report the tree-level value for honesty.
        alpha_val = alpha_tree

        results.append({
            "label": label,
            "mu_eV": mu,
            "alpha_2": alpha_val,
            "alpha_3": alpha_val,
            "are_equal": True,
            "n_active_KK": n_active,
            "l_max": l_max,
            "regime": regime,
        })

    return results


# ---------------------------------------------------------------------------
# 5. unification_analysis
# ---------------------------------------------------------------------------
def unification_analysis(a_0_m=28e-6, phi_vev=None):
    """Analyze the SO(3) gauge unification at the KK scale.

    Below m_KK: two separate U(1) gauge fields from the coset, equal coupling.
    At m_KK: SO(3) gauge symmetry is restored -- all three directions active.
    Above m_KK: SO(3) gauge theory in 4D (plus massive KK modes).
    Above M_6: full 6D theory.

    This is NOT like GUT unification (10^16 GeV).  The "unification" here is
    just the restoration of the compactification symmetry at the KK scale,
    and there is nothing to unify since the couplings never split.

    Parameters
    ----------
    a_0_m   : float
    phi_vev : float or None

    Returns
    -------
    dict
    """
    if phi_vev is None:
        phi_vev = PHI_VEV

    alpha_tree = math.exp(4.0 * phi_vev) / (4.0 * math.pi)
    R_phys = _R_phys(a_0_m)
    m_kk_eV = _m_kk(a_0_m)
    M_6_eV = _M_6_from_R(R_phys)

    # GUT scale for comparison
    M_GUT_eV = 2.0e25  # ~2 x 10^16 GeV

    return {
        "unification_scale_eV": m_kk_eV,
        "unification_scale_GeV": m_kk_eV / 1.0e9,
        "coupling_at_unification": alpha_tree,
        "M_6_eV": M_6_eV,
        "ratio_to_GUT": M_GUT_eV / m_kk_eV,
        "comparison_with_GUT": {
            "GUT_scale_GeV": M_GUT_eV / 1.0e9,
            "KK_scale_GeV": m_kk_eV / 1.0e9,
            "orders_of_magnitude_apart": math.log10(M_GUT_eV / m_kk_eV)
            if m_kk_eV > 0 else float("inf"),
            "mechanism": {
                "GUT": "Three SM gauge couplings (alpha_1, alpha_2, alpha_3) run "
                       "from different values at M_Z and converge at ~10^16 GeV due "
                       "to different beta functions from different matter content.",
                "S2": "Two coset gauge couplings are ALWAYS equal (SO(2) symmetry). "
                      "There is nothing to converge. The SO(3) restoration at m_KK "
                      "is a geometric effect, not a dynamical unification.",
            },
        },
        "description": (
            "The SO(3) gauge symmetry is restored at the KK scale m_KK = "
            "sqrt(2)/R_phys. This is not a dynamical unification like a GUT: "
            "the two coset couplings were never different, and the 'unification' "
            "is simply the geometric restoration of the full isometry group. "
            "At energies above m_KK, the non-propagating SO(2) direction "
            "acquires a kinetic term from the KK tower, and all three SO(3) "
            "gauge fields become dynamical."
        ),
    }


# ---------------------------------------------------------------------------
# 6. two_brane_splitting
# ---------------------------------------------------------------------------
def two_brane_splitting(theta_1=0.0, theta_2_values=None):
    """Analyze gauge coupling splitting with two branes on S^2.

    For a SINGLE brane at any point on S^2, one can always choose coordinates
    so the brane sits at the north pole.  The residual SO(2) (rotation around
    the pole) ensures A^2 and A^3 have identical couplings to brane matter.

    For TWO branes at different positions, the coordinate freedom is exhausted
    by placing Brane 1 at the north pole.  Brane 2 sits at polar angle theta_2.
    The SO(2) rotation around the north pole maps:
        A^2 -> cos(chi) A^2 + sin(chi) A^3
        A^3 -> -sin(chi) A^2 + cos(chi) A^3
    for azimuthal rotation by chi.

    On Brane 1 (north pole): the coupling is SO(2)-symmetric, so both fields
    couple equally.  No splitting on Brane 1.

    On Brane 2 (at theta_2, phi_azimuth=0 by convention): the brane position
    breaks the SO(2) freedom only if we demand BOTH branes to be special.
    But the coupling of the coset fields to a brane at (theta, 0) is:

        The coset Killing vectors at (theta, phi) are:
        K_2 = -sin(phi) d/d(theta) - cos(phi)*cot(theta) d/d(phi)
        K_3 = cos(phi) d/d(theta) - sin(phi)*cot(theta) d/d(phi)

    At phi_azimuth=0:
        |K_2|^2 = sin^2(0) + cos^2(0)*cot^2(theta) = cot^2(theta)
        |K_3|^2 = cos^2(0) + sin^2(0)*cot^2(theta) = 1

    Wait -- the norms depend on how we decompose. The gauge fields couple
    to brane currents via A^a_mu * J^mu_a. The coupling strength depends on
    the Killing vector overlap with the brane position.

    For a brane at the north pole (theta=0), both K_2 and K_3 have equal
    magnitude by SO(2) symmetry.

    For a brane at general theta_2 with phi_azimuth=0:
    The effective coupling to each field depends on the projection of the
    Killing vectors onto the brane.  However, the KEY POINT is:

    Even with two branes, the BULK gauge kinetic matrix is still diag(0,1,1).
    The brane-localized terms add corrections to the effective 4D couplings.
    These corrections are proportional to the brane tensions and the Killing
    vector norms at the brane positions.

    For a round S^2, the Killing vectors of SO(3) satisfy:
        sum_a K^a(x) K^a(x) = 2   (constant, independent of position)

    This is a consequence of the homogeneity of S^2.  Therefore, for ANY
    single point, the sum |K_2|^2 + |K_3|^2 is the same (it equals 2 minus
    |K_1|^2, where K_1 = d/d(phi) at the north pole).

    The individual norms |K_2|^2 and |K_3|^2 CAN differ at a generic point.

    At (theta, 0):
        K_1 = d/d(phi) => |K_1|^2 = sin^2(theta) (metric on S^2)
        |K_2|^2 + |K_3|^2 = 2 - sin^2(theta)

    But |K_2|^2 and |K_3|^2 individually depend on the azimuthal angle phi
    of the brane, which we set to 0.  At phi=0:
        K_2 propto d/d(theta)  => |K_2|^2 = 1  (from the g_{theta theta}=1 on unit S^2)
        K_3 propto cot(theta) d/d(phi) => |K_3|^2 = cos^2(theta)/sin^2(theta) * sin^2(theta) = cos^2(theta)

    Hmm, this requires more care.  The Killing vectors of SO(3) on the unit
    S^2 with metric ds^2 = d(theta)^2 + sin^2(theta) d(phi)^2 are:

        K_1 = d/d(phi)
        K_2 = -sin(phi) d/d(theta) - cos(phi) cos(theta)/sin(theta) d/d(phi)
        K_3 = cos(phi) d/d(theta) - sin(phi) cos(theta)/sin(theta) d/d(phi)

    At phi=0:
        K_2 = -cos(0)*cos(theta)/sin(theta) d/d(phi) = -cot(theta) d/d(phi)
              => |K_2|^2 = cot^2(theta) * sin^2(theta) = cos^2(theta)
        K_3 = cos(0) d/d(theta) = d/d(theta)
              => |K_3|^2 = 1

    So at (theta, phi=0): |K_2|^2 = cos^2(theta), |K_3|^2 = 1.

    The effective coupling of the gauge field A^a to brane matter is
    proportional to |K_a|^2 at the brane position (the Killing vector
    norm determines how much the gauge transformation acts on the brane).

    alpha_a_eff = alpha_tree * |K_a(brane)|^2 / |K_a(pole)|^2

    At the north pole (theta->0): |K_2|^2 -> 1, |K_3|^2 -> 1 (both equal).
    At general theta: |K_2|^2 = cos^2(theta), |K_3|^2 = 1.

    The splitting on Brane 2:
        delta_alpha / alpha = (|K_3|^2 - |K_2|^2) / (|K_3|^2 + |K_2|^2)
                            = (1 - cos^2(theta)) / (1 + cos^2(theta))
                            = sin^2(theta) / (1 + cos^2(theta))

    Parameters
    ----------
    theta_1         : float, polar angle of Brane 1 (default 0 = north pole)
    theta_2_values  : list of float or None

    Returns
    -------
    list of dicts
    """
    if theta_2_values is None:
        theta_2_values = [
            0.0,
            math.pi / 6.0,
            math.pi / 4.0,
            math.pi / 3.0,
            math.pi / 2.0,
            2.0 * math.pi / 3.0,
            3.0 * math.pi / 4.0,
            math.pi,
        ]

    alpha_tree = math.exp(4.0 * PHI_VEV) / (4.0 * math.pi)

    results = []

    for theta_2 in theta_2_values:
        # Brane 1 at north pole: |K_2|^2 = |K_3|^2 = 1 (by SO(2) symmetry)
        alpha_2_brane1 = alpha_tree
        alpha_3_brane1 = alpha_tree

        # Brane 2 at (theta_2, phi=0):
        # |K_2|^2 = cos^2(theta_2), |K_3|^2 = 1
        # Handle theta_2 = 0 (north pole) and theta_2 = pi (south pole)
        if abs(theta_2) < 1e-12 or abs(theta_2 - math.pi) < 1e-12:
            # At poles: both Killing vectors have equal norms by symmetry
            K2_sq = 1.0
            K3_sq = 1.0
        else:
            K2_sq = math.cos(theta_2) ** 2
            K3_sq = 1.0

        # Effective couplings on Brane 2
        # The brane-localized gauge coupling gets weighted by the
        # Killing vector norm at the brane position.
        alpha_2_brane2 = alpha_tree * K2_sq
        alpha_3_brane2 = alpha_tree * K3_sq

        # Splitting measure
        if (K2_sq + K3_sq) > 0:
            splitting = abs(K3_sq - K2_sq) / (K3_sq + K2_sq)
        else:
            splitting = 0.0

        results.append({
            "theta_2": theta_2,
            "theta_2_deg": math.degrees(theta_2),
            "K2_squared": K2_sq,
            "K3_squared": K3_sq,
            "alpha_2_brane1": alpha_2_brane1,
            "alpha_3_brane1": alpha_3_brane1,
            "alpha_2_brane2": alpha_2_brane2,
            "alpha_3_brane2": alpha_3_brane2,
            "splitting": splitting,
            "are_equal_brane1": True,
            "are_equal_brane2": abs(K2_sq - K3_sq) < 1e-12,
        })

    return results


# ---------------------------------------------------------------------------
# 7. electroweak_analogy
# ---------------------------------------------------------------------------
def electroweak_analogy():
    """Compare SO(3) -> SO(2) on S^2 with SU(2)_L -> U(1)_EM in the SM.

    Structural similarities and key differences between:
    - Geometric breaking: SO(3) -> SO(2) by compactification on S^2
    - Electroweak breaking: SU(2)_L x U(1)_Y -> U(1)_EM by the Higgs

    Returns
    -------
    dict
    """
    # Group theory data
    so3_data = {
        "group": "SO(3)",
        "dimension": 3,
        "rank": 1,
        "generators": 3,
        "unbroken_subgroup": "SO(2) ~ U(1)",
        "unbroken_generators": 1,
        "broken_generators": 2,
        "coset": "SO(3)/SO(2) ~ S^2",
        "coset_dimension": 2,
        "casimir_so3": "l(l+1) for spin-l representation",
        "branching_rule": "spin-l -> m = -l, ..., +l under SO(2)",
    }

    su2_data = {
        "group": "SU(2)_L x U(1)_Y",
        "dimension": 4,
        "rank": 2,
        "generators": 4,
        "unbroken_subgroup": "U(1)_EM",
        "unbroken_generators": 1,
        "broken_generators": 3,
        "coset": "SU(2)xU(1) / U(1)_EM ~ S^3",
        "coset_dimension": 3,
        "casimir_su2": "T(T+1) for isospin-T representation",
        "branching_rule": "isospin-T -> Q = T_3 + Y/2 under U(1)_EM",
    }

    comparison = [
        {
            "aspect": "Parent group",
            "S2": "SO(3)",
            "EW": "SU(2)_L x U(1)_Y",
            "similar": False,
            "note": "SO(3) ~ SU(2)/Z_2 locally, but S^2 has no hypercharge",
        },
        {
            "aspect": "Unbroken subgroup",
            "S2": "SO(2) ~ U(1)",
            "EW": "U(1)_EM",
            "similar": True,
            "note": "Both leave a U(1) unbroken",
        },
        {
            "aspect": "Breaking mechanism",
            "S2": "Geometric (topology of S^2)",
            "EW": "Spontaneous (Higgs VEV)",
            "similar": False,
            "note": "This is the KEY difference: S^2 is explicit, Higgs is spontaneous",
        },
        {
            "aspect": "Goldstone bosons",
            "S2": "None (explicit breaking)",
            "EW": "3 eaten by W+, W-, Z",
            "similar": False,
            "note": "Geometric breaking does not produce Goldstones",
        },
        {
            "aspect": "Massive vectors",
            "S2": "KK tower (massive at m_l = sqrt(l(l+1))/R)",
            "EW": "W+, W-, Z (massive by Higgs mechanism)",
            "similar": True,
            "note": "Both have massive vector bosons, but by different mechanisms",
        },
        {
            "aspect": "Massless vector fate",
            "S2": "SO(2) direction: NON-PROPAGATING (zero eigenvalue)",
            "EW": "U(1)_EM: massless photon (PROPAGATING)",
            "similar": False,
            "note": "OPPOSITE fate: in S^2, the unbroken direction has no kinetic "
                    "term; in EW, the unbroken direction gives the photon",
        },
        {
            "aspect": "Number of physical gauge fields",
            "S2": "2 (coset directions)",
            "EW": "4 (W+, W-, Z, gamma)",
            "similar": False,
            "note": "S^2 has fewer physical fields",
        },
        {
            "aspect": "Coupling unification",
            "S2": "Trivial (couplings always equal)",
            "EW": "Non-trivial (g and g' differ; sin^2(theta_W) ~ 0.23)",
            "similar": False,
            "note": "S^2 has no analog of the Weinberg angle",
        },
    ]

    return {
        "so3_group_theory": so3_data,
        "su2_group_theory": su2_data,
        "comparison": comparison,
        "structural_analogy": True,
        "physical_analogy": False,
        "key_difference": (
            "In electroweak breaking, the unbroken U(1)_EM gives the massless, "
            "propagating photon. In the S^2 compactification, the unbroken SO(2) "
            "direction has ZERO eigenvalue in the gauge kinetic matrix and does "
            "NOT propagate. The fates of the unbroken direction are exactly "
            "opposite. This is because EW breaking is spontaneous (Goldstones "
            "eaten to make W/Z massive, photon remains), while S^2 breaking is "
            "explicit (geometry dictates which directions get kinetic terms)."
        ),
    }


# ---------------------------------------------------------------------------
# 8. summarize_gauge_splitting
# ---------------------------------------------------------------------------
def summarize_gauge_splitting(a_0_m=28e-6):
    """Full summary of gauge coupling splitting analysis.

    Parameters
    ----------
    a_0_m : float, coordinate radius in meters

    Returns
    -------
    dict with all sub-analyses and summary text
    """
    kinetic = gauge_kinetic_matrix()
    tree = tree_level_couplings()
    beta = one_loop_beta_functions(a_0_m=a_0_m)
    evolution = coupling_evolution(a_0_m=a_0_m)
    unification = unification_analysis(a_0_m=a_0_m)
    two_brane = two_brane_splitting()
    ew = electroweak_analogy()

    # Check: is there any splitting anywhere?
    any_splitting_single_brane = False  # Never, by SO(2) symmetry
    max_two_brane_splitting = max(r["splitting"] for r in two_brane)

    summary_lines = [
        "GAUGE COUPLING SPLITTING ANALYSIS: SO(3)/SO(2) ON S^2",
        "=" * 56,
        "",
        "1. GAUGE KINETIC MATRIX",
        "   Eigenvalues: {} (rank {})".format(kinetic["eigenvalues"], kinetic["rank"]),
        "   SO(2) direction: index {} (non-propagating)".format(
            kinetic["so2_direction"]
        ),
        "   Coset directions: indices {} (physical)".format(
            kinetic["coset_directions"]
        ),
        "",
        "2. TREE-LEVEL COUPLINGS",
        "   alpha_2 = alpha_3 = {:.6e}  (= alpha_EM)".format(tree["alpha_1"]),
        "   Equal: {} -- {}".format(tree["are_equal"], tree["reason"][:80]),
        "",
        "3. ONE-LOOP BETA FUNCTIONS",
        "   b_2 = {}, b_3 = {}".format(beta["b_2"], beta["b_3"]),
        "   Equal: {}".format(beta["are_equal"]),
        "   In minimal framework: gauge coupling does NOT run below m_KK",
        "",
        "4. COUPLING EVOLUTION",
    ]

    for ev in evolution:
        kk_str = "{:d}".format(ev["n_active_KK"]) if ev["n_active_KK"] < 100000 \
            else "{:.2e}".format(ev["n_active_KK"])
        summary_lines.append(
            "   {:<8s}  alpha_2 = {:.6e}  alpha_3 = {:.6e}  "
            "equal={!s:<5}  KK={:>10s}  {}".format(
                ev["label"],
                ev["alpha_2"],
                ev["alpha_3"],
                ev["are_equal"],
                kk_str,
                ev["regime"],
            )
        )

    summary_lines.extend([
        "",
        "5. UNIFICATION ANALYSIS",
        "   KK unification scale: {:.4e} eV ({:.4e} GeV)".format(
            unification["unification_scale_eV"],
            unification["unification_scale_GeV"],
        ),
        "   Coupling at unification: {:.6e}".format(
            unification["coupling_at_unification"]
        ),
        "   GUT scale / KK scale: {:.2e} ({:.1f} orders of magnitude)".format(
            unification["ratio_to_GUT"],
            unification["comparison_with_GUT"]["orders_of_magnitude_apart"],
        ),
        "   Nature: geometric restoration, NOT dynamical unification",
        "",
        "6. TWO-BRANE SPLITTING",
    ])

    for tb in two_brane:
        summary_lines.append(
            "   theta_2 = {:>6.1f} deg  |K_2|^2 = {:.4f}  |K_3|^2 = {:.4f}  "
            "splitting = {:.4f}  equal_brane2 = {}".format(
                tb["theta_2_deg"],
                tb["K2_squared"],
                tb["K3_squared"],
                tb["splitting"],
                tb["are_equal_brane2"],
            )
        )

    summary_lines.extend([
        "   Max splitting (two-brane): {:.4f}".format(max_two_brane_splitting),
        "",
        "7. ELECTROWEAK ANALOGY",
        "   Structural analogy: {}".format(ew["structural_analogy"]),
        "   Physical analogy: {}".format(ew["physical_analogy"]),
        "   Key difference: {}".format(ew["key_difference"][:100]),
        "",
        "8. CONCLUSIONS",
        "   - Single brane: NO splitting at any scale (exact SO(2) symmetry)",
        "   - Minimal framework: NO running (no charged matter)",
        "   - Two branes: splitting possible, up to {:.0f}% for equatorial brane".format(
            max_two_brane_splitting * 100
        ),
        "   - Unification is geometric (m_KK), not dynamical (GUT)",
        "   - The framework makes NO prediction for coupling unification",
        "   - The SO(3)->SO(2) pattern is structurally, not physically, like EW",
    ])

    summary_text = "\n".join(summary_lines)

    return {
        "kinetic_matrix": kinetic,
        "tree_level": tree,
        "beta_functions": beta,
        "evolution": evolution,
        "unification": unification,
        "two_brane": two_brane,
        "electroweak": ew,
        "any_splitting_single_brane": any_splitting_single_brane,
        "max_two_brane_splitting": max_two_brane_splitting,
        "summary_text": summary_text,
    }


# ---------------------------------------------------------------------------
# __main__
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    result = summarize_gauge_splitting()
    print(result["summary_text"])
