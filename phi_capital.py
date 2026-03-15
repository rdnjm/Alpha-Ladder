"""
Investigate the bridge coefficient Phi = alpha_G / alpha^21.
Is it a fundamental constant in its own right?
"""

from decimal import Decimal, getcontext
getcontext().prec = 50

alpha = Decimal('0.0072973525664')
alpha_g = Decimal('1.7512e-45')
pi = Decimal('3.14159265358979323846264338327950288419716939937510')
phi = (1 + Decimal(5).sqrt()) / 2  # golden ratio
e = Decimal('2.71828182845904523536028747135266249775724709369995')

alpha_21 = alpha ** 21
Phi = alpha_g / alpha_21

print(f"Phi = alpha_G / alpha^21")
print(f"Phi = {Phi}")
print(f"Phi (float) = {float(Phi):.15f}")
print()

# -------------------------------------------------------------------
# How sensitive is Phi to the uncertainty in G?
# G = 6.67430(15) x 10^-11 (CODATA 2018, ~22 ppm uncertainty)
# alpha_G = G * m_e^2 / (hbar * c)
# So alpha_G inherits G's ~22 ppm relative uncertainty
# -------------------------------------------------------------------
m_e = Decimal('9.10938370150e-31')   # kg
hbar = Decimal('1.054571817e-34')    # J*s
c = Decimal('299792458')              # m/s

G_central = Decimal('6.67430e-11')
G_low = Decimal('6.67415e-11')       # -1 sigma
G_high = Decimal('6.67445e-11')      # +1 sigma

def alpha_g_from_G(G):
    return G * m_e ** 2 / (hbar * c)

alpha_g_calc = alpha_g_from_G(G_central)
Phi_central = alpha_g_calc / alpha_21
Phi_low = alpha_g_from_G(G_low) / alpha_21
Phi_high = alpha_g_from_G(G_high) / alpha_21

print("--- Sensitivity to G uncertainty ---")
print(f"  G = {G_central}  =>  alpha_G = {alpha_g_calc:.6e}  =>  Phi = {float(Phi_central):.10f}")
print(f"  G = {G_low} (-1s) =>  alpha_G = {alpha_g_from_G(G_low):.6e}  =>  Phi = {float(Phi_low):.10f}")
print(f"  G = {G_high} (+1s) =>  alpha_G = {alpha_g_from_G(G_high):.6e}  =>  Phi = {float(Phi_high):.10f}")
print(f"  Phi range: [{float(Phi_low):.10f}, {float(Phi_high):.10f}]")
print(f"  Phi uncertainty: +/- {float((Phi_high - Phi_low) / 2):.2e}")
print()

# -------------------------------------------------------------------
# What known constants fall in the Phi range?
# -------------------------------------------------------------------
candidates = {
    "phi^2 / 2  = (phi+1)/2":          float(phi ** 2 / 2),
    "(5/12) * pi":                      5.0 / 12.0 * float(pi),
    "sqrt(e) / cbrt(2)":               float(e.sqrt() / Decimal(2) ** (Decimal(1)/Decimal(3))),
    "4/3  (exact)":                     4.0 / 3.0,
    "ln(phi) + 1":                      float(phi.ln() + 1),
    "1 + 1/pi":                         1.0 + 1.0 / float(pi),
    "sqrt(phi) * phi^(-1/4)":          float(phi) ** 0.5 * float(phi) ** (-0.25),
    "2*phi - 2":                        2.0 * float(phi) - 2.0,
    "phi / sqrt(phi+1)":               float(phi / (phi + 1).sqrt()),
    "e / (e - phi/pi)":                float(e) / (float(e) - float(phi)/float(pi)),
    "7/4 - phi/4":                     7.0/4.0 - float(phi)/4.0,
    "1 + pi/e^2":                       1.0 + float(pi) / float(e)**2,
    "cos(pi/5) * 2  = phi":            float(phi),  # for reference
    "cos(1)^(-1)  = sec(1)":           1.0 / float(Decimal(1).copy_abs()),  # need math
}

# Use math for trig
import math
candidates["sec(1 radian)"] = 1.0 / math.cos(1.0)
candidates["2*cos(pi/5)"] = 2.0 * math.cos(math.pi / 5)
candidates["sqrt(phi) + 1/phi^2"] = math.sqrt(float(phi)) + 1.0/float(phi)**2
candidates["1 + Catalan/pi"] = 1.0 + 0.9159655941 / math.pi  # Catalan's constant

# Clean up bad ones and add more
del candidates["cos(1)^(-1)  = sec(1)"]
del candidates["cos(pi/5) * 2  = phi"]

Phi_f = float(Phi)
print("--- Candidate expressions for Phi ---")
print(f"    Target Phi = {Phi_f:.15f}")
print()

ranked = []
for name, val in candidates.items():
    err = abs(val - Phi_f) / Phi_f * 100
    ranked.append((err, name, val))

ranked.sort()
for err, name, val in ranked:
    marker = " <<<" if err < 0.05 else ""
    print(f"  {err:10.4f}%  |  {val:.15f}  |  {name}{marker}")

# -------------------------------------------------------------------
# Key question: does Phi fall within G's error bar for any expression?
# -------------------------------------------------------------------
print()
print("--- Do any expressions fall within G's uncertainty window? ---")
print(f"    Phi window: [{float(Phi_low):.10f}, {float(Phi_high):.10f}]")
print()
for err, name, val in ranked[:10]:
    inside = float(Phi_low) <= val <= float(Phi_high)
    status = "INSIDE" if inside else "outside"
    print(f"  {status:7s}  |  {val:.15f}  |  {name}")
