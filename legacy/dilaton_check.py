"""
Check: does the dilaton interpretation of rung 21 hold quantitatively?

In Brans-Dicke theory:  G_eff = 1/Phi_BD  (Phi_BD is the BD scalar field)
In string theory:       G ~ g_s^2 * alpha'^(D-2)/2  where g_s = e^<dilaton>
In Kaluza-Klein:        G_4D ~ G_5D / R  where R is the compactification radius

Question: if we write alpha_G = alpha^(20+1) * C, does the "20" connect
to the Riemann tensor and the "+1" to a scalar (dilaton) contribution?
"""

import math
from decimal import Decimal, getcontext
getcontext().prec = 50

alpha = Decimal('0.0072973525693')
phi = (1 + Decimal(5).sqrt()) / Decimal(2)

# -----------------------------------------------------------------------
# Decompose alpha^21 = alpha^20 * alpha^1
# The bridge: alpha_G = (phi^2/2) * alpha^21
#           = (phi^2/2) * alpha^20 * alpha
#
# Interpretation:
#   alpha^20  <-- 20 Riemann components (geometric/local)
#   alpha^1   <-- 1 scalar component (dilaton/global)
#   phi^2/2   <-- the bridge coefficient (vacuum geometry)
# -----------------------------------------------------------------------

alpha_20 = alpha ** 20
alpha_1 = alpha
bridge = phi ** 2 / Decimal(2)

print("=" * 70)
print("  DILATON DECOMPOSITION OF THE GRAVITY RUNG")
print("=" * 70)
print()
print(f"  alpha_G = (phi^2/2) * alpha^21")
print(f"         = (phi^2/2) * alpha^20 * alpha^1")
print()
print(f"  Geometric sector (alpha^20):  {alpha_20:.6e}")
print(f"  Dilaton sector   (alpha^1):   {alpha_1:.6e}")
print(f"  Bridge coeff     (phi^2/2):   {float(bridge):.10f}")
print()

# -----------------------------------------------------------------------
# Alternative decompositions: is 20+1 special vs other splits?
# If 20+1 is physically meaningful, the "1" part should be simpler
# or more natural than other splits like 19+2, 18+3, etc.
# -----------------------------------------------------------------------
print("  Alternative decompositions of 21 = a + b:")
print(f"  {'Split':<10s}  {'alpha^a':<14s}  {'alpha^b':<14s}  {'Geometric meaning of a'}")
print(f"  {'-'*10}  {'-'*14}  {'-'*14}  {'-'*40}")

splits = {
    (20, 1): "Riemann(4D) + scalar",
    (19, 2): "???",
    (18, 3): "???",
    (15, 6): "Riemann(3D)=6 + 15 ???",
    (14, 7): "G2 Lie algebra + 7 ???",
    (12, 9): "SM gauge bosons + 9 ???",
    (10, 11): "Ricci(4D) or Weyl(4D) + 11",
    (6, 15):  "metric(3D) + Riemann(3D)=6 ???",
}

for (a, b), meaning in splits.items():
    print(f"  {a}+{b:<7d}  {float(alpha**a):<14.6e}  {float(alpha**b):<14.6e}  {meaning}")

# -----------------------------------------------------------------------
# Key test: in Brans-Dicke, the dilaton couples with strength 1/(2*omega+3)
# where omega is the BD parameter. GR is omega -> infinity.
# If phi governs the dilaton, what omega does that imply?
# -----------------------------------------------------------------------
print()
print("=" * 70)
print("  BRANS-DICKE PARAMETER FROM PHI")
print("=" * 70)
print()

# In BD theory: G_measured = G_bare * (2*omega + 4) / (2*omega + 3)
# The deviation from GR is of order 1/omega
# Solar system tests require omega > 40,000 (Cassini 2003)
#
# But if phi^2/2 IS the dilaton contribution, then the BD-like correction is:
# G_eff / G_bare = phi^2/2 ~ 1.309
# That's a 30% deviation from "bare" gravity, which would mean omega ~ 1
# This is EXCLUDED by solar system tests unless the dilaton is massive
# (short-range, screened at solar system scales)

bare_ratio = float(bridge)
# (2w+4)/(2w+3) = phi^2/2 => 2w+4 = (phi^2/2)(2w+3) => ...
# 2w + 4 = phi^2*w + 3*phi^2/2
# w(2 - phi^2) = 3*phi^2/2 - 4
# w = (3*phi^2/2 - 4) / (2 - phi^2)
phi_f = float(phi)
phi2 = phi_f ** 2
numerator = 3 * phi2 / 2 - 4
denominator = 2 - phi2
if denominator != 0:
    omega = numerator / denominator
    print(f"  If G_eff/G_bare = phi^2/2 = {bare_ratio:.6f}")
    print(f"  Then BD parameter omega = {omega:.4f}")
    print()
    if omega < 0:
        print(f"  omega is NEGATIVE ({omega:.2f})")
        print(f"  This is allowed in some scalar-tensor theories")
        print(f"  (e.g., string-inspired models with phantom dilaton)")
    if abs(omega) < 40000:
        print(f"  |omega| = {abs(omega):.0f} << 40,000 (Cassini bound)")
        print(f"  This is EXCLUDED for a massless dilaton.")
        print(f"  BUT: if the dilaton has mass m_d >> 1/AU,")
        print(f"  the solar system constraint doesn't apply.")
        print()
        # What dilaton mass would screen at solar system scales?
        # Compton wavelength < 1 AU = 1.496e11 m
        # m_d > hbar / (c * 1 AU)
        hbar_f = 1.054571817e-34
        c_f = 2.99792458e8
        AU = 1.496e11
        m_d_min = hbar_f / (c_f * AU)
        m_d_min_eV = m_d_min * c_f**2 / 1.602176634e-19
        print(f"  Minimum dilaton mass to evade Cassini:")
        print(f"    m_d > {m_d_min:.2e} kg")
        print(f"    m_d > {m_d_min_eV:.2e} eV")
        print(f"    Compton wavelength < 1 AU")
        print()

        # What scale does alpha^10 (dark matter rung) correspond to as a mass?
        # If alpha_DM ~ alpha^10 sets a coupling, the associated mass scale
        # in natural units would be m ~ alpha^5 * m_Pl (geometric mean)
        m_Pl_eV = 1.22089e28  # eV
        m_dark_scale = float(alpha) ** 5 * m_Pl_eV
        print(f"  For comparison, alpha^5 * m_Planck = {m_dark_scale:.2e} eV")
        print(f"  (a possible dark/dilaton mass scale from the ladder)")

# -----------------------------------------------------------------------
# The 20+1 vs 6D metric
# -----------------------------------------------------------------------
print()
print("=" * 70)
print("  RECONCILING 20+1 WITH 6D")
print("=" * 70)
print()
print("  Earlier we found: 21 = symmetric rank-2 tensor components in 6D")
print("  Now we find:      21 = 20 (Riemann in 4D) + 1 (dilaton)")
print()
print("  These are COMPATIBLE if the 6D metric decomposes as:")
print("    g_MN (6D, 21 components) = g_uv (4D, 10) + g_ab (2D, 3) + g_ua (mixed, 8)")
print(f"    Total: 10 + 3 + 8 = 21")
print()
print("  The Riemann tensor of 4D (20 components) + 1 scalar is a DIFFERENT")
print("  counting, but both give 21. This suggests two complementary views:")
print()
print("    View A (bottom-up):  4D gravity needs 20+1 to couple to the SM")
print("    View B (top-down):   6D gravity has 21 metric components naturally")
print()
print("  Both views point to the same rung.")
