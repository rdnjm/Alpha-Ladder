"""
Brute force: predict Newton's gravitational constant G from the Alpha Ladder.

If alpha_G = C * alpha^21, and alpha_G = G * m_e^2 / (hbar * c), then:
    G = C * alpha^21 * hbar * c / m_e^2

alpha, hbar, c, m_e are all known to sub-ppm precision.
C is the bridge coefficient -- the only unknown.
Each candidate for C predicts a different value of G.
"""

from decimal import Decimal, getcontext

getcontext().prec = 50

# -----------------------------------------------------------------------
# CODATA 2018 constants (high precision, NOT dependent on G)
# -----------------------------------------------------------------------
alpha = Decimal('0.0072973525693')       # 1.5e-10 relative uncertainty
hbar = Decimal('1.054571817e-34')        # J*s, exact in 2019 SI
c = Decimal('299792458')                 # m/s, exact by definition
m_e = Decimal('9.1093837015e-31')        # kg, 3e-10 relative uncertainty

alpha_21 = alpha ** 21

# -----------------------------------------------------------------------
# Bridge candidates
# -----------------------------------------------------------------------
pi = Decimal('3.14159265358979323846264338327950288419716939937510')
phi = (1 + Decimal(5).sqrt()) / 2
e = Decimal('2.71828182845904523536028747135266249775724709369995')

candidates = {
    "phi^2 / 2":            phi ** 2 / 2,
    "(5/12) * pi":          Decimal(5) / Decimal(12) * pi,
    "sqrt(e) / cbrt(2)":    e.sqrt() / Decimal(2) ** (Decimal(1) / Decimal(3)),
}

# -----------------------------------------------------------------------
# Predict G from each candidate
# -----------------------------------------------------------------------
def predict_G(bridge_coeff):
    alpha_g = bridge_coeff * alpha_21
    return alpha_g * hbar * c / m_e ** 2

# -----------------------------------------------------------------------
# Experimental G measurements (various groups, chronological)
# -----------------------------------------------------------------------
G_measurements = {
    "CODATA 2018 recommended":          (Decimal('6.67430e-11'), Decimal('0.00015e-11')),
    "Quinn et al. 2013 (BIPM)":         (Decimal('6.67545e-11'), Decimal('0.00018e-11')),
    "Rosi et al. 2014 (atom interf.)":  (Decimal('6.67191e-11'), Decimal('0.00099e-11')),
    "Newman et al. 2014":               (Decimal('6.67435e-11'), Decimal('0.00013e-11')),
    "Li et al. 2018 (HUST-A)":          (Decimal('6.67418e-11'), Decimal('0.00009e-11')),
    "Li et al. 2018 (HUST-B)":          (Decimal('6.67484e-11'), Decimal('0.00009e-11')),
    "CODATA 2014 recommended":          (Decimal('6.67408e-11'), Decimal('0.00031e-11')),
}

print("=" * 80)
print("  ALPHA LADDER -- PREDICTED VALUES OF NEWTON'S G")
print("=" * 80)
print()
print("  Formula:  G = C * alpha^21 * hbar * c / m_e^2")
print("  where C is the bridge coefficient alpha_G / alpha^21")
print()

for name, coeff in candidates.items():
    G_pred = predict_G(coeff)
    print(f"  Bridge: {name}")
    print(f"  C = {float(coeff):.15f}")
    print(f"  G = {G_pred:.5e} m^3 kg^-1 s^-2")
    print()

    # Compare against every experimental measurement
    for exp_name, (G_exp, G_unc) in G_measurements.items():
        diff = G_pred - G_exp
        sigma = float(abs(diff) / G_unc) if G_unc != 0 else float('inf')
        direction = "+" if diff > 0 else "-"
        print(f"    vs {exp_name:<36s}  {direction}{abs(sigma):5.1f} sigma  "
              f"(delta = {float(diff):.3e})")
    print()

# -----------------------------------------------------------------------
# Summary: which candidate is most consistent across ALL measurements?
# -----------------------------------------------------------------------
print("=" * 80)
print("  SUMMARY: Average |sigma| across all measurements")
print("=" * 80)
print()

for name, coeff in candidates.items():
    G_pred = predict_G(coeff)
    total_sigma = 0
    count = 0
    for exp_name, (G_exp, G_unc) in G_measurements.items():
        diff = G_pred - G_exp
        sigma = float(abs(diff) / G_unc)
        total_sigma += sigma
        count += 1
    avg = total_sigma / count
    print(f"  {name:<28s}  avg |sigma| = {avg:.1f}")

print()

# -----------------------------------------------------------------------
# The prediction
# -----------------------------------------------------------------------
G_best = predict_G(phi ** 2 / 2)
print("=" * 80)
print("  THE PREDICTION (assuming phi^2/2 bridge)")
print("=" * 80)
print()
print(f"  G_predicted = {G_best:.8e} m^3 kg^-1 s^-2")
print()
print(f"  This value is derived ENTIRELY from:")
print(f"    alpha   = {alpha}  (QED, sub-ppb)")
print(f"    hbar    = {hbar}  (exact, 2019 SI)")
print(f"    c       = {c}  (exact, by definition)")
print(f"    m_e     = {m_e}  (sub-ppb)")
print(f"    phi     = (1+sqrt(5))/2  (exact, mathematical)")
print()
print(f"  Predicted precision: limited only by m_e uncertainty (~3e-10)")
print(f"  That is ~70,000x more precise than current G measurements.")
print("=" * 80)
