"""
Scan for phi coincidences in precision-measured Standard Model parameters.
Looking for a second independent appearance of phi to corroborate the bridge.
"""

import math

phi = (1 + math.sqrt(5)) / 2
pi = math.pi

# -----------------------------------------------------------------------
# Precision SM parameters (PDG/CODATA values)
# -----------------------------------------------------------------------
params = {
    "sin^2(theta_W) (MS-bar, MZ)":  0.23122,       # +/- 0.00004
    "sin^2(theta_W) (on-shell)":    0.22337,        # +/- 0.00010
    "Cabibbo angle sin(theta_C)":   0.22500,        # Vus
    "Weinberg angle theta_W (rad)": math.asin(math.sqrt(0.23122)),
    "m_p / m_e":                    1836.15267343,
    "m_mu / m_e":                   206.7682830,
    "m_tau / m_e":                  3477.23,
    "m_W / m_Z":                    0.88147,
    "m_H / m_Z":                    1.37186,
    "m_H / m_W":                    1.55636,
    "alpha_s(MZ)":                  0.1180,         # strong coupling
    "alpha(MZ) (EM at Z pole)":     1/127.951,
    "G_F (Fermi, reduced)":         1.1663788e-5,   # GeV^-2
    "Euler-Mascheroni gamma":       0.5772156649,
}

# -----------------------------------------------------------------------
# Search: is each parameter close to phi^(p/q), pi^(p/q), or
# simple combinations (a/b) * phi^(p/q)?
# -----------------------------------------------------------------------
print("=" * 85)
print("  PHI COINCIDENCE SCANNER")
print("=" * 85)
print()

alpha_em = 0.0072973525693

for pname, pval in params.items():
    if pval <= 0:
        continue

    best_err = float('inf')
    best_expr = ""
    best_val = 0

    # phi^(p/q)
    for p in range(-8, 9):
        for q in range(1, 7):
            if p == 0:
                continue
            candidate = phi ** (p / q)
            err = abs(candidate - pval) / abs(pval) * 100
            if err < best_err:
                best_err = err
                best_val = candidate
                best_expr = f"phi^({p}/{q})" if q != 1 else f"phi^{p}"

    # (a/b) * phi^(p/q)
    for a in range(1, 13):
        for b in range(1, 13):
            if a == b:
                continue
            frac = a / b
            for p in range(-6, 7):
                for q in range(1, 5):
                    if p == 0:
                        candidate = frac
                    else:
                        candidate = frac * phi ** (p / q)
                    err = abs(candidate - pval) / abs(pval) * 100
                    if err < best_err:
                        best_err = err
                        best_val = candidate
                        ep = f"phi^({p}/{q})" if q != 1 else f"phi^{p}"
                        best_expr = f"({a}/{b}) * {ep}" if p != 0 else f"{a}/{b}"

    # alpha^n * phi^(p/q)
    for n in range(-3, 4):
        if n == 0:
            continue
        for p in range(-4, 5):
            for q in range(1, 4):
                candidate = (alpha_em ** n) * (phi ** (p / q))
                err = abs(candidate - pval) / abs(pval) * 100
                if err < best_err:
                    best_err = err
                    best_val = candidate
                    ep = f"phi^({p}/{q})" if (q != 1 and p != 0) else (f"phi^{p}" if p != 0 else "")
                    an = f"alpha^{n}"
                    best_expr = f"{an} * {ep}" if ep else an

    marker = " <<<" if best_err < 0.5 else (" <" if best_err < 2.0 else "")
    print(f"  {pname:<38s} = {pval:<14.8g}")
    print(f"    best: {best_expr:<30s} = {best_val:<14.8g}  err = {best_err:.4f}%{marker}")
    print()

# -----------------------------------------------------------------------
# Deep dive: sin^2(theta_W) vs phi
# -----------------------------------------------------------------------
print("=" * 85)
print("  DEEP DIVE: sin^2(theta_W) and phi")
print("=" * 85)
print()

sw2 = 0.23122
candidates_sw2 = {
    "1/phi^3":                      1 / phi**3,
    "3/(4*phi^3)":                  3 / (4 * phi**3),
    "(phi-1)/phi^2":                (phi - 1) / phi**2,
    "1/(phi^2 + phi)":              1 / (phi**2 + phi),
    "phi^2 / (2*pi + phi)":         phi**2 / (2*pi + phi),
    "1/(2*phi + 1)":                1 / (2*phi + 1),
    "1/phi^3 * (1 - alpha)":        (1/phi**3) * (1 - alpha_em),
    "(3 - phi) / (4 + phi)":        (3 - phi) / (4 + phi),
    "(phi - 1)^3":                  (phi - 1)**3,
    "1/(phi + pi)":                 1 / (phi + pi),
    "pi/(4*phi^2 - 1)":            pi / (4*phi**2 - 1),
    "phi/(phi^2 + pi)":            phi / (phi**2 + pi),
}

print(f"  sin^2(theta_W) = {sw2} +/- 0.00004")
print(f"  Experimental range: [{sw2 - 0.00004:.5f}, {sw2 + 0.00004:.5f}]")
print()

ranked = []
for name, val in candidates_sw2.items():
    err = abs(val - sw2) / sw2 * 100
    inside = abs(val - sw2) < 0.00004
    ranked.append((err, name, val, inside))

ranked.sort()
for err, name, val, inside in ranked:
    status = "INSIDE ERROR BAR" if inside else ""
    print(f"  {err:8.4f}%  |  {val:.8f}  |  {name:<28s}  {status}")

# -----------------------------------------------------------------------
# Riemann tensor connection
# -----------------------------------------------------------------------
print()
print("=" * 85)
print("  THE 21 CONNECTION")
print("=" * 85)
print()
print(f"  Riemann curvature tensor in 4D: {4*4*(4*4-1)//12} independent components = 20")
print(f"  ... but Ricci tensor has 10, Weyl tensor has 10, total Riemann = 20")
print()
print(f"  Actually: n*(n^2-1)/12 for n=4 dimensions:")
print(f"    4 * (16-1) / 12 = 4 * 15 / 12 = {4*15//12} ... that's 5, not 21")
print()
print(f"  Correct formula: n^2 * (n^2 - 1) / 12 for n=4:")
print(f"    16 * 15 / 12 = {16*15//12}")
print()
print(f"  So Riemann tensor has 20 independent components in 4D, not 21.")
print()
print(f"  But: 21 = triangular number T(6) = 6*7/2")
print(f"  And: 21 = C(7,2) = ways to choose 2 items from 7")
print(f"  And: 21 = dimension of SO(7) Lie algebra")
print(f"  And: 21 = number of components in a symmetric 6x6 matrix")
print(f"  And: 21 = dimension of SU(3) x SU(2) x U(1) adjoint rep (8+3+1=12... no)")
print()

# What IS 21 in terms of SM gauge group?
print(f"  SM gauge group dimensions:")
print(f"    SU(3): dim = 8  (gluons)")
print(f"    SU(2): dim = 3  (W+, W-, Z precursors)")
print(f"    U(1):  dim = 1  (hypercharge)")
print(f"    Total: 8 + 3 + 1 = 12 gauge bosons")
print(f"    Plus Higgs doublet (4 real) + gravity (1?) = 12 + 4 + 1 = 17... no")
print()
print(f"  However: independent components of a SYMMETRIC RANK-2 TENSOR in D dims:")
for d in range(3, 8):
    comps = d * (d + 1) // 2
    print(f"    D={d}: {comps} components" + ("  <-- 21!" if comps == 21 else ""))
