"""
Running of the Electromagnetic Fine-Structure Constant with KK Thresholds
=========================================================================

Computes the one-loop running of alpha_EM from the Thomson limit up to
arbitrary energy scales, including Standard Model fermion thresholds and
(optionally) Kaluza-Klein charged scalar modes from S^2 compactification.

Key results
-----------
* Below the KK scale, only SM fermion loops drive the running.
  alpha(m_e) = 1/137.036, alpha(m_Z) = 1/127.9 -- fully accounted for.
* In the MINIMAL Alpha Ladder framework (pure gravity on S^2), the KK
  modes are uncharged under U(1)_EM, so there is NO KK correction to
  alpha running.
* If charged matter is added (Paper 3 extension), each 6D charged scalar
  produces a tower of massive charged scalars at m_l = sqrt(l(l+1)) / R
  with degeneracy (2l+1).  For a_0 ~ 28 um the lightest KK mode is at
  ~18 meV -- far below the electron mass -- so these modes modify alpha
  at ALL observable scales.  This is a stringent consistency constraint.

Pure Python -- only ``import math``, no numpy/scipy.

Functions
---------
1. sm_alpha_running        -- alpha_EM(mu) from SM fermion thresholds
2. kk_alpha_running        -- alpha_EM(mu) with KK charged scalar tower
3. compare_sm_vs_kk        -- side-by-side table at key energy scales
4. does_it_match_observation -- consistency check against alpha(m_e), alpha(m_Z)
5. running_profile         -- plot-ready grid of 1/alpha(mu) vs mu
6. summarize_alpha_running -- full physics report
"""

import math

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
ALPHA_EM_0 = 1.0 / 137.035999084       # Thomson limit (CODATA 2018)
ALPHA_EM_MZ = 1.0 / 127.944            # at Z pole (PDG)
M_E = 5.11e5                           # electron mass (eV)
M_MU = 1.057e8                         # muon mass (eV)
M_TAU = 1.7769e9                       # tau mass (eV)
M_U = 2.2e6                            # up quark (eV)
M_D = 4.7e6                            # down quark (eV)
M_S = 9.5e7                            # strange quark (eV)
M_C = 1.27e9                           # charm quark (eV)
M_B = 4.18e9                           # bottom quark (eV)
M_T = 1.73e11                          # top quark (eV)
M_Z = 9.12e10                          # Z boson mass (eV)
M_PL_EV = 1.22089e28                   # Planck mass (eV)
HBAR_C_EVM = 1.9733e-7                 # hbar*c (eV*m)
PHI_GOLDEN = (1.0 + math.sqrt(5.0)) / 2.0

# SM charged fermions: (name, mass_eV, |Q|, N_c)
SM_FERMIONS = [
    ("electron", M_E,   1.0,   1),
    ("up",       M_U,   2.0/3, 3),
    ("down",     M_D,   1.0/3, 3),
    ("muon",     M_MU,  1.0,   1),
    ("strange",  M_S,   1.0/3, 3),
    ("charm",    M_C,   2.0/3, 3),
    ("tau",      M_TAU, 1.0,   1),
    ("bottom",   M_B,   1.0/3, 3),
    ("top",      M_T,   2.0/3, 3),
]
# Sort by mass for threshold stepping
SM_FERMIONS.sort(key=lambda f: f[1])


def _b_coefficient(Q, Nc):
    """One-loop QED beta-function coefficient for a Dirac fermion.

    d(1/alpha)/d(ln mu) = -(2/(3*pi)) * Nc * Q^2
    per active fermion.  Returns the coefficient Nc*Q^2.
    """
    return Nc * Q * Q


# ===================================================================
# 1. sm_alpha_running
# ===================================================================
def sm_alpha_running(mu_eV):
    """Compute alpha_EM(mu) using one-loop QED running with SM thresholds.

    Parameters
    ----------
    mu_eV : float
        Energy scale in eV.

    Returns
    -------
    dict with keys:
        alpha       -- alpha_EM at scale mu
        inv_alpha   -- 1/alpha_EM at scale mu
        active      -- list of active fermion names at scale mu
        mu_eV       -- the input scale
        description -- human-readable summary
    """
    # Reference point: alpha at the electron mass
    # We define alpha(m_e) = alpha_Thomson = 1/137.036 (the difference
    # between alpha(0) and alpha(m_e) is < 0.01 and irrelevant here).
    inv_alpha_ref = 1.0 / ALPHA_EM_0
    mu_ref = M_E  # reference scale = electron mass

    if mu_eV <= 0:
        return {
            "alpha": ALPHA_EM_0,
            "inv_alpha": inv_alpha_ref,
            "active": [],
            "mu_eV": mu_eV,
            "description": "Scale <= 0: returning Thomson limit.",
        }

    if mu_eV < mu_ref:
        # Below electron mass: no charged fermions active, alpha frozen.
        return {
            "alpha": ALPHA_EM_0,
            "inv_alpha": inv_alpha_ref,
            "active": [],
            "mu_eV": mu_eV,
            "description": "Below m_e: alpha frozen at Thomson limit.",
        }

    # Step through thresholds from m_e upward
    inv_alpha = inv_alpha_ref
    current_mu = mu_ref
    active = []

    for name, mass, Q, Nc in SM_FERMIONS:
        if mu_eV <= mass:
            # Haven't reached this threshold yet; run from current_mu to mu_eV
            # with currently active fermions, then stop.
            break
        # This fermion becomes active at its mass threshold.
        # First, run from current_mu to this threshold with existing active set.
        if mass > current_mu:
            b_sum = sum(_b_coefficient(Q_f, Nc_f)
                        for _, _, Q_f, Nc_f in SM_FERMIONS
                        if (_, _, Q_f, Nc_f) in _dummy_never_matches(active))
            # Actually: sum over already-active fermions
            b_active = 0.0
            for aname in active:
                for fn, fm, fQ, fNc in SM_FERMIONS:
                    if fn == aname:
                        b_active += _b_coefficient(fQ, fNc)
                        break
            delta = -(2.0 / (3.0 * math.pi)) * b_active * math.log(mass / current_mu)
            inv_alpha += delta
            current_mu = mass
        active.append(name)
    else:
        # All fermions activated; we exited the loop normally.
        pass

    # Final run from current_mu (last activated threshold) to mu_eV
    b_active = 0.0
    for aname in active:
        for fn, fm, fQ, fNc in SM_FERMIONS:
            if fn == aname:
                b_active += _b_coefficient(fQ, fNc)
                break
    if mu_eV > current_mu:
        delta = -(2.0 / (3.0 * math.pi)) * b_active * math.log(mu_eV / current_mu)
        inv_alpha += delta

    alpha = 1.0 / inv_alpha if inv_alpha > 0 else float('inf')

    return {
        "alpha": alpha,
        "inv_alpha": inv_alpha,
        "active": list(active),
        "mu_eV": mu_eV,
        "description": (
            "1/alpha({:.3e} eV) = {:.4f},  alpha = {:.8e},  "
            "{} active fermions.".format(mu_eV, inv_alpha, alpha, len(active))
        ),
    }


def _dummy_never_matches(active):
    """Placeholder that returns empty -- used only to satisfy a dead branch."""
    return []


# ===================================================================
# 2. kk_alpha_running
# ===================================================================
def kk_alpha_running(mu_eV, a_0_m, n_charged=1, Q_charge=1.0, matched=True,
                     l_max=None):
    """Compute alpha_EM(mu) including KK charged scalar contributions.

    Above the KK threshold, each level-l mode of a 6D charged scalar on S^2
    contributes to the running.  Mass: m_l = sqrt(l(l+1)) / R_phys,
    degeneracy: 2l+1, spin-0 (complex scalar).

    Parameters
    ----------
    mu_eV     : float  -- energy scale in eV
    a_0_m     : float  -- compactification radius in meters
    n_charged : int    -- number of distinct 6D charged scalars (0 = minimal)
    Q_charge  : float  -- electric charge of the 6D scalar (in units of e)
    matched   : bool   -- if True, use dilaton-matched R_phys
    l_max     : int    -- max angular momentum to include (None = auto)

    Returns
    -------
    dict with keys:
        alpha, inv_alpha, n_kk_active, m_kk_lightest_eV, active_sm,
        mu_eV, description
    """
    # Physical radius
    if matched:
        phi_vev = 0.25 * math.log(4.0 * math.pi * ALPHA_EM_0)
        R_phys = a_0_m * math.exp(phi_vev)
    else:
        R_phys = a_0_m

    # Lightest KK mass: l=1, m_1 = sqrt(2) / R_phys  (in natural units)
    # Convert: m (eV) = hbar*c / R (m) * sqrt(l(l+1))
    m_scale_eV = HBAR_C_EVM / R_phys  # = hbar*c / R in eV

    m_kk_lightest = m_scale_eV * math.sqrt(2.0)  # l=1

    # Start with SM running up to mu
    sm = sm_alpha_running(mu_eV)
    inv_alpha = sm["inv_alpha"]
    active_sm = sm["active"]

    if n_charged == 0:
        return {
            "alpha": sm["alpha"],
            "inv_alpha": inv_alpha,
            "n_kk_active": 0,
            "m_kk_lightest_eV": m_kk_lightest,
            "active_sm": active_sm,
            "mu_eV": mu_eV,
            "description": (
                "Minimal framework (n_charged=0): no KK correction to alpha.  "
                + sm["description"]
            ),
        }

    # Determine l_max: include modes with m_l < mu_eV
    if l_max is None:
        # m_l = m_scale * sqrt(l(l+1)), so l(l+1) < (mu/m_scale)^2
        ratio_sq = (mu_eV / m_scale_eV) ** 2 if m_scale_eV > 0 else 0
        # l^2 + l - ratio_sq < 0  =>  l < (-1 + sqrt(1 + 4*ratio_sq)) / 2
        if ratio_sq > 0:
            l_max_float = (-1.0 + math.sqrt(1.0 + 4.0 * ratio_sq)) / 2.0
            l_max = int(math.floor(l_max_float))
        else:
            l_max = 0

    if l_max < 1:
        # No KK modes active
        return {
            "alpha": sm["alpha"],
            "inv_alpha": inv_alpha,
            "n_kk_active": 0,
            "m_kk_lightest_eV": m_kk_lightest,
            "active_sm": active_sm,
            "mu_eV": mu_eV,
            "description": (
                "Below KK threshold: no active KK modes.  " + sm["description"]
            ),
        }

    # KK contribution: each level l (l >= 1, up to l_max) contributes
    # delta(1/alpha) = -(1/(3*pi)) * Q^2 * (2l+1) * n_charged * ln(mu/m_l)
    # where m_l = m_scale * sqrt(l(l+1)).
    #
    # Total: -(n_charged*Q^2/(3*pi)) * Sum_{l=1}^{L} (2l+1)*ln(mu/(m_scale*sqrt(l(l+1))))
    #      = -(n_charged*Q^2/(3*pi)) * [ln(mu/m_scale)*S1 - 0.5*S2]
    # where S1 = Sum (2l+1) = L(L+2)
    #       S2 = Sum (2l+1)*ln(l(l+1))
    #
    # For S2, use the identity:
    #   Sum_{l=1}^{L} (2l+1)*ln(l(l+1)) = Sum (2l+1)*ln(l) + Sum (2l+1)*ln(l+1)
    # and the key relation: Sum_{l=1}^{L} (2l+1)*ln(l) can be written as
    #   2*Sum l*ln(l) + Sum ln(l) = 2*Sum l*ln(l) + ln(L!)
    # We compute S2 via Stirling for large L, or exactly for small L.

    Q2 = Q_charge * Q_charge
    L = l_max

    # Total number of active KK modes (counting degeneracy)
    n_kk_active = L * (L + 2)  # Sum_{l=1}^{L} (2l+1) = L^2+2L

    # For moderate L (< 50000), compute S2 exactly by direct summation.
    # For very large L, use the analytical approximation.
    USE_EXACT_THRESHOLD = 50000

    if L <= USE_EXACT_THRESHOLD:
        S2 = 0.0
        for l in range(1, L + 1):
            S2 += (2 * l + 1) * math.log(l * (l + 1))
    else:
        # Analytical approximation for S2 = Sum_{l=1}^{L} (2l+1)*ln(l(l+1))
        # Using Euler-Maclaurin / Stirling:
        # Sum (2l+1)*ln(l(l+1)) ~ integral_1^L (2x+1)*ln(x(x+1)) dx + corrections
        #
        # Integral of (2x+1)*ln(x(x+1)) dx:
        # Let u = x(x+1) = x^2+x, du = (2x+1)dx
        # So integral = integral ln(u) du = u*ln(u) - u
        # = x(x+1)*ln(x(x+1)) - x(x+1) evaluated from 1 to L
        #
        uL = L * (L + 1)
        u1 = 1 * 2  # = 2
        integral = (uL * math.log(uL) - uL) - (u1 * math.log(u1) - u1)

        # Euler-Maclaurin corrections (first two terms):
        # f(L) = (2L+1)*ln(L(L+1)), f(1) = 3*ln(2)
        fL = (2 * L + 1) * math.log(L * (L + 1))
        f1 = 3.0 * math.log(2.0)
        em_correction = 0.5 * (fL + f1)  # trapezoidal correction

        S2 = integral + em_correction

    S1 = L * (L + 2)
    ln_mu_over_mscale = math.log(mu_eV / m_scale_eV)

    delta_kk = -(n_charged * Q2 / (3.0 * math.pi)) * (
        ln_mu_over_mscale * S1 - 0.5 * S2
    )

    # But wait: the SM running was computed with 1/alpha(m_e) = 137.036 as
    # input.  If KK modes are active below m_e, they shift what 1/alpha(m_e)
    # SHOULD be (the observed value already includes their effect if they exist).
    # For a self-consistent treatment, we add the KK contribution on top.
    # The consistency check is done in does_it_match_observation().

    inv_alpha_kk = inv_alpha + delta_kk
    alpha_kk = 1.0 / inv_alpha_kk if inv_alpha_kk > 0 else float('inf')

    return {
        "alpha": alpha_kk,
        "inv_alpha": inv_alpha_kk,
        "n_kk_active": n_kk_active,
        "m_kk_lightest_eV": m_kk_lightest,
        "active_sm": active_sm,
        "mu_eV": mu_eV,
        "description": (
            "1/alpha({:.3e} eV) = {:.4f} (SM+KK),  "
            "{} active KK modes (n_charged={}, Q={}).  "
            "Lightest KK: {:.3e} eV."
            .format(mu_eV, inv_alpha_kk, n_kk_active, n_charged,
                    Q_charge, m_kk_lightest)
        ),
    }


# ===================================================================
# 3. compare_sm_vs_kk
# ===================================================================
def compare_sm_vs_kk(a_0_m, n_charged=1, matched=True):
    """Compare alpha at key energy scales with and without KK modes.

    Parameters
    ----------
    a_0_m     : float -- compactification radius in meters
    n_charged : int   -- number of 6D charged scalars
    matched   : bool  -- use dilaton-matched R_phys

    Returns
    -------
    dict with keys:
        table : list of dicts, each with mu_eV, label, alpha_SM, alpha_KK,
                delta_alpha_ppm, n_kk_active, inv_alpha_SM, inv_alpha_KK
        a_0_m, m_kk_eV, description
    """
    # Compute lightest KK mass
    if matched:
        phi_vev = 0.25 * math.log(4.0 * math.pi * ALPHA_EM_0)
        R_phys = a_0_m * math.exp(phi_vev)
    else:
        R_phys = a_0_m

    m_kk = HBAR_C_EVM / R_phys * math.sqrt(2.0)

    # M_6 from the formula M_6^4 = M_Pl^2 / (4*pi*R^2)
    # in eV: M_6 = (M_Pl^2 / (4*pi*R^2))^(1/4) ... but R in eV^{-1}
    R_eV_inv = R_phys / HBAR_C_EVM  # R in eV^{-1}
    M_6_eV4 = M_PL_EV ** 2 / (4.0 * math.pi * R_eV_inv ** 2)
    M_6_eV = M_6_eV4 ** 0.25

    scales = [
        ("m_e",         M_E),
        ("1 GeV",       1.0e9),
        ("m_Z",         M_Z),
        ("1 TeV",       1.0e12),
        ("m_KK",        m_kk),
        ("10 m_KK",     10.0 * m_kk),
        ("100 m_KK",    100.0 * m_kk),
        ("M_6",         M_6_eV),
    ]

    # Sort by energy
    scales.sort(key=lambda s: s[1])

    table = []
    for label, mu in scales:
        sm = sm_alpha_running(mu)
        kk = kk_alpha_running(mu, a_0_m, n_charged=n_charged, matched=matched)
        delta_ppm = 0.0
        if sm["alpha"] > 0 and kk["alpha"] > 0:
            delta_ppm = (kk["alpha"] - sm["alpha"]) / sm["alpha"] * 1.0e6
        table.append({
            "label": label,
            "mu_eV": mu,
            "alpha_SM": sm["alpha"],
            "inv_alpha_SM": sm["inv_alpha"],
            "alpha_KK": kk["alpha"],
            "inv_alpha_KK": kk["inv_alpha"],
            "delta_alpha_ppm": delta_ppm,
            "n_kk_active": kk["n_kk_active"],
        })

    return {
        "table": table,
        "a_0_m": a_0_m,
        "m_kk_eV": m_kk,
        "M_6_eV": M_6_eV,
        "description": (
            "Comparison of alpha_EM running: SM-only vs SM+KK.  "
            "a_0 = {:.3e} m,  m_KK = {:.3e} eV,  M_6 = {:.3e} eV."
            .format(a_0_m, m_kk, M_6_eV)
        ),
    }


# ===================================================================
# 4. does_it_match_observation
# ===================================================================
def does_it_match_observation(a_0_m, n_charged=1, Q_charge=1.0,
                              matched=True):
    """Check whether the KK tower is consistent with observed alpha values.

    The observed alpha(m_e) = 1/137.036 and alpha(m_Z) = 1/127.944 are
    fully accounted for by SM running.  Any additional KK charged modes
    active below m_Z would shift these values.

    For the minimal framework (n_charged=0), there is no issue.
    For n_charged >= 1 with a_0 ~ 28 um, the KK scale is ~18 meV,
    so the entire tower is active at m_e and above.

    Parameters
    ----------
    a_0_m     : float -- compactification radius in meters
    n_charged : int   -- number of 6D charged scalars
    Q_charge  : float -- charge of the 6D scalar
    matched   : bool  -- dilaton-matched radius

    Returns
    -------
    dict with keys:
        alpha_predicted_at_me  -- alpha(m_e) with SM + KK
        alpha_predicted_at_mZ  -- alpha(m_Z) with SM + KK
        alpha_observed_at_me   -- 1/137.036
        alpha_observed_at_mZ   -- 1/127.944
        shift_at_me_ppm        -- fractional shift at m_e in ppm
        shift_at_mZ_ppm        -- fractional shift at m_Z in ppm
        consistent             -- True if shifts < 100 ppm
        m_kk_eV                -- lightest KK mass
        description            -- summary string
    """
    kk_me = kk_alpha_running(M_E, a_0_m, n_charged=n_charged,
                              Q_charge=Q_charge, matched=matched)
    kk_mz = kk_alpha_running(M_Z, a_0_m, n_charged=n_charged,
                              Q_charge=Q_charge, matched=matched)
    sm_me = sm_alpha_running(M_E)
    sm_mz = sm_alpha_running(M_Z)

    # The "predicted" values if KK modes exist: the SM running already
    # gives the right answer, and KK adds on top.  The shift is the
    # KK-only part.
    shift_me = (kk_me["inv_alpha"] - sm_me["inv_alpha"])
    shift_mz = (kk_mz["inv_alpha"] - sm_mz["inv_alpha"])

    # In ppm of 1/alpha
    shift_me_ppm = shift_me / sm_me["inv_alpha"] * 1.0e6
    shift_mz_ppm = shift_mz / sm_mz["inv_alpha"] * 1.0e6

    # alpha(m_e) is measured to ~0.37 ppb, alpha(m_Z) to ~150 ppm.
    # Any shift > 100 ppm at m_e is inconsistent.
    threshold_ppm = 100.0  # generous threshold
    consistent = (abs(shift_me_ppm) < threshold_ppm and
                  abs(shift_mz_ppm) < threshold_ppm)

    if n_charged == 0:
        consistent = True
        desc = (
            "Minimal framework (n_charged=0): no charged KK modes.  "
            "SM running accounts for alpha exactly.  Fully consistent."
        )
    else:
        desc = (
            "With {} charged scalar(s) (Q={}):  ".format(n_charged, Q_charge)
            + "delta(1/alpha) at m_e = {:.2f} (shift {:.1f} ppm).  ".format(
                shift_me, shift_me_ppm)
            + "delta(1/alpha) at m_Z = {:.2f} (shift {:.1f} ppm).  ".format(
                shift_mz, shift_mz_ppm)
            + "Lightest KK mode: {:.3e} eV.  ".format(kk_me["m_kk_lightest_eV"])
            + ("CONSISTENT." if consistent else
               "INCONSISTENT: KK tower over-runs alpha.")
        )

    return {
        "alpha_predicted_at_me": kk_me["alpha"],
        "alpha_predicted_at_mZ": kk_mz["alpha"],
        "alpha_observed_at_me": ALPHA_EM_0,
        "alpha_observed_at_mZ": ALPHA_EM_MZ,
        "inv_alpha_shift_at_me": shift_me,
        "inv_alpha_shift_at_mZ": shift_mz,
        "shift_at_me_ppm": shift_me_ppm,
        "shift_at_mZ_ppm": shift_mz_ppm,
        "consistent": consistent,
        "m_kk_eV": kk_me["m_kk_lightest_eV"],
        "n_charged": n_charged,
        "Q_charge": Q_charge,
        "description": desc,
    }


# ===================================================================
# 5. running_profile
# ===================================================================
def running_profile(a_0_m, n_charged=1, mu_min_eV=1e-3, mu_max_eV=1e16,
                    n_points=200, matched=True):
    """Compute 1/alpha(mu) on a log-spaced grid for plotting.

    Parameters
    ----------
    a_0_m      : float -- compactification radius (m)
    n_charged  : int   -- number of 6D charged scalars
    mu_min_eV  : float -- minimum energy scale (eV)
    mu_max_eV  : float -- maximum energy scale (eV)
    n_points   : int   -- number of grid points
    matched    : bool  -- dilaton-matched radius

    Returns
    -------
    dict with keys:
        mu_eV          : list of energy scales
        inv_alpha_SM   : list of 1/alpha from SM only
        inv_alpha_KK   : list of 1/alpha from SM + KK
        m_kk_eV        : lightest KK mass
        description    : summary
    """
    log_min = math.log10(mu_min_eV)
    log_max = math.log10(mu_max_eV)
    step = (log_max - log_min) / (n_points - 1) if n_points > 1 else 0

    mu_list = []
    inv_alpha_sm = []
    inv_alpha_kk = []

    for i in range(n_points):
        log_mu = log_min + i * step
        mu = 10.0 ** log_mu
        mu_list.append(mu)

        sm = sm_alpha_running(mu)
        inv_alpha_sm.append(sm["inv_alpha"])

        kk = kk_alpha_running(mu, a_0_m, n_charged=n_charged, matched=matched)
        inv_alpha_kk.append(kk["inv_alpha"])

    # KK lightest mass
    if matched:
        phi_vev = 0.25 * math.log(4.0 * math.pi * ALPHA_EM_0)
        R_phys = a_0_m * math.exp(phi_vev)
    else:
        R_phys = a_0_m
    m_kk = HBAR_C_EVM / R_phys * math.sqrt(2.0)

    return {
        "mu_eV": mu_list,
        "inv_alpha_SM": inv_alpha_sm,
        "inv_alpha_KK": inv_alpha_kk,
        "m_kk_eV": m_kk,
        "n_points": n_points,
        "description": (
            "Running profile: {} points from {:.1e} to {:.1e} eV.  "
            "m_KK = {:.3e} eV."
            .format(n_points, mu_min_eV, mu_max_eV, m_kk)
        ),
    }


# ===================================================================
# 6. summarize_alpha_running
# ===================================================================
def summarize_alpha_running():
    """Generate a comprehensive summary of alpha running physics.

    Returns
    -------
    dict with keys:
        sm_check        -- SM running check at m_Z
        minimal_check   -- minimal framework consistency
        charged_check   -- charged matter consistency (a_0=28um)
        comparison      -- comparison table at a_0=28um
        key_messages    -- list of key physics takeaways
    """
    lines = []
    lines.append("=" * 72)
    lines.append("RUNNING OF ALPHA_EM WITH KK THRESHOLD CORRECTIONS")
    lines.append("=" * 72)

    # --- SM running check ---
    lines.append("")
    lines.append("1. STANDARD MODEL RUNNING")
    lines.append("-" * 40)
    sm_me = sm_alpha_running(M_E)
    sm_mz = sm_alpha_running(M_Z)
    lines.append("  alpha(m_e)  = 1/{:.4f}  (input: 1/137.036)"
                 .format(sm_me["inv_alpha"]))
    lines.append("  alpha(m_Z)  = 1/{:.4f}  (observed: 1/127.944)"
                 .format(sm_mz["inv_alpha"]))
    lines.append("  Active fermions at m_Z: {}".format(
        ", ".join(sm_mz["active"])))

    # Discrepancy at m_Z
    disc_mz = sm_mz["inv_alpha"] - 127.944
    lines.append("  Discrepancy at m_Z: delta(1/alpha) = {:.2f}"
                 .format(disc_mz))
    lines.append("  (Expected: one-loop QED misses EW corrections, "
                 "hadronic vacuum polarization, etc.)")

    sm_check = {
        "inv_alpha_at_me": sm_me["inv_alpha"],
        "inv_alpha_at_mZ": sm_mz["inv_alpha"],
        "inv_alpha_mZ_observed": 127.944,
        "discrepancy": disc_mz,
    }

    # --- Minimal framework ---
    lines.append("")
    lines.append("2. MINIMAL FRAMEWORK (no charged matter)")
    lines.append("-" * 40)
    a_0 = 28.0e-6  # 28 um
    minimal = does_it_match_observation(a_0, n_charged=0, matched=True)
    lines.append("  " + minimal["description"])
    lines.append("  Result: KK modes from pure gravity are uncharged.")
    lines.append("  No correction to alpha_EM.  Fully consistent.")

    # --- Charged matter extension ---
    lines.append("")
    lines.append("3. CHARGED MATTER EXTENSION (Paper 3)")
    lines.append("-" * 40)
    for nc in [1, 2, 5, 10]:
        check = does_it_match_observation(a_0, n_charged=nc, Q_charge=1.0,
                                          matched=True)
        status = "OK" if check["consistent"] else "EXCLUDED"
        lines.append("  n_charged={:2d}, Q=1:  shift(m_e) = {:+.2e} ppm,  "
                     "shift(m_Z) = {:+.2e} ppm  [{}]"
                     .format(nc, check["shift_at_me_ppm"],
                             check["shift_at_mZ_ppm"], status))

    charged_check_1 = does_it_match_observation(a_0, n_charged=1, matched=True)
    lines.append("")
    lines.append("  Lightest KK mass (a_0=28um): {:.3e} eV = {:.3e} meV"
                 .format(charged_check_1["m_kk_eV"],
                         charged_check_1["m_kk_eV"] * 1e3))

    # --- Comparison table ---
    lines.append("")
    lines.append("4. COMPARISON TABLE (a_0 = 28 um, n_charged=1)")
    lines.append("-" * 40)
    comp = compare_sm_vs_kk(a_0, n_charged=1, matched=True)
    lines.append("  {:>12s}  {:>12s}  {:>14s}  {:>14s}  {:>12s}"
                 .format("Scale", "1/alpha_SM", "1/alpha_KK",
                         "delta(ppm)", "n_KK"))
    for row in comp["table"]:
        # Format 1/alpha_KK: use scientific notation if extreme
        inv_kk = row["inv_alpha_KK"]
        if abs(inv_kk) > 1e6:
            inv_kk_str = "{:+.3e}".format(inv_kk)
        else:
            inv_kk_str = "{:+.4f}".format(inv_kk)
        # Format delta ppm
        dppm = row["delta_alpha_ppm"]
        if inv_kk < 0:
            dppm_str = "Landau pole"
        elif abs(dppm) > 1e9:
            dppm_str = "{:+.2e}".format(dppm)
        else:
            dppm_str = "{:+.1f}".format(dppm)
        # Format n_KK
        nkk = row["n_kk_active"]
        if nkk > 1e6:
            nkk_str = "{:.2e}".format(nkk)
        else:
            nkk_str = str(nkk)
        lines.append("  {:>12s}  {:12.4f}  {:>14s}  {:>14s}  {:>12s}"
                     .format(row["label"], row["inv_alpha_SM"],
                             inv_kk_str, dppm_str, nkk_str))

    # --- Key messages ---
    lines.append("")
    lines.append("5. KEY MESSAGES")
    lines.append("-" * 40)
    key_messages = [
        ("SM running from alpha(m_e)=1/137 to alpha(m_Z)=1/128 is "
         "well-understood.  Our one-loop QED gives 1/{:.1f} at m_Z "
         "(misses ~{:.1f} from EW/hadronic effects)."
         .format(sm_mz["inv_alpha"], abs(disc_mz))),

        ("In the MINIMAL Alpha Ladder framework, KK modes are gravitational "
         "(uncharged under U(1)_EM).  No correction to alpha running."),

        ("With charged matter (Paper 3 extension): KK modes at ~{:.0f} meV "
         "(for a_0=28um) are active below m_e.  They modify alpha at ALL "
         "observable scales."
         .format(charged_check_1["m_kk_eV"] * 1e3)),

        ("This is a NEW consistency constraint: the charged KK tower must "
         "not shift 1/alpha(m_e) by more than ~100 ppm (current precision "
         "is 0.37 ppb)."),

        ("For Q=1: each charged scalar shifts 1/alpha(m_e) by ~{:.2e} ppm.  "
         "Even ONE charged 6D scalar is wildly excluded."
         .format(abs(charged_check_1["shift_at_me_ppm"]))),
    ]
    for i, msg in enumerate(key_messages, 1):
        lines.append("  {}. {}".format(i, msg))

    report = "\n".join(lines)

    return {
        "sm_check": sm_check,
        "minimal_check": minimal,
        "charged_check": charged_check_1,
        "comparison": comp,
        "key_messages": key_messages,
        "report": report,
    }


# ===================================================================
# __main__ -- print full report
# ===================================================================
if __name__ == "__main__":
    result = summarize_alpha_running()
    print(result["report"])
