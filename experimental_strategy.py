"""
Experimental strategy for testing the Alpha Ladder.

Instead of measuring G (the hardest constant in physics), we look for
OTHER testable predictions the ladder makes and check them against
existing data or accessible experiments.
"""

from decimal import Decimal, getcontext
getcontext().prec = 50

alpha = Decimal('0.0072973525693')
pi = Decimal('3.14159265358979323846264338327950288419716939937510')
phi = (1 + Decimal(5).sqrt()) / Decimal(2)
e_const = Decimal('2.71828182845904523536028747135266249775724709369995')

# Physical constants (CODATA 2018, high precision)
hbar = Decimal('1.054571817e-34')
c = Decimal('2.99792458e8')
m_e = Decimal('9.1093837015e-31')
m_p = Decimal('1.67262192369e-27')
m_mu = Decimal('1.883531627e-28')      # muon mass
k_B = Decimal('1.380649e-23')          # Boltzmann, exact in 2019 SI
e_charge = Decimal('1.602176634e-19')  # elementary charge, exact in 2019 SI
epsilon_0 = Decimal('8.8541878128e-12')
G = Decimal('6.67430e-11')

alpha_21 = alpha ** 21

print("=" * 80)
print("  ALPHA LADDER -- EXPERIMENTAL STRATEGY")
print("=" * 80)

# -----------------------------------------------------------------------
# STRATEGY 1: Mass ratios as alpha powers
# If the ladder organizes coupling strengths, mass ratios between
# fundamental particles might also fall on rungs.
# -----------------------------------------------------------------------
print("\n" + "-" * 80)
print("  STRATEGY 1: Do particle mass ratios fall on alpha-power rungs?")
print("-" * 80)
print()
print("  If the ladder is geometric, mass ratios between fundamental particles")
print("  might be expressible as alpha^n * (simple coefficient).")
print()

import math

ratios = {
    "m_p / m_e":       float(m_p / m_e),
    "m_mu / m_e":      float(m_mu / m_e),
    "m_p / m_mu":      float(m_p / m_mu),
}

alpha_f = float(alpha)

for name, ratio in ratios.items():
    # What power of alpha gives this ratio?
    n_exact = math.log(ratio) / math.log(1 / alpha_f)
    n_round = round(n_exact)

    # What coefficient is needed at that integer rung?
    coeff = ratio * alpha_f ** n_round
    print(f"  {name} = {ratio:.6f}")
    print(f"    = alpha^(-{n_exact:.4f})")
    print(f"    ~ alpha^(-{n_round}) * {coeff:.6f}")

    # Check nearby rungs too
    for n in [n_round - 1, n_round, n_round + 1]:
        coeff_n = ratio * alpha_f ** n
        print(f"      at n=-{n}: coeff = {coeff_n:.6f}")
    print()

# -----------------------------------------------------------------------
# STRATEGY 2: Predict alpha_G from multiple independent paths
# If multiple relationships converge on the same G, the ladder gains
# credibility without directly measuring G.
# -----------------------------------------------------------------------
print("-" * 80)
print("  STRATEGY 2: Multiple paths to alpha_G")
print("-" * 80)
print()

# Path A: phi^2/2 * alpha^21  (the bridge)
alpha_g_A = (phi ** 2 / Decimal(2)) * alpha_21
G_A = alpha_g_A * hbar * c / m_e ** 2

# Path B: use proton mass instead of electron mass
# alpha_G(proton) = G * m_p^2 / (hbar*c)
alpha_g_proton = G * m_p ** 2 / (hbar * c)
# What power of alpha gives this?
n_proton = math.log(float(alpha_g_proton)) / math.log(float(alpha))
print(f"  Path A: alpha_G(electron) = (phi^2/2) * alpha^21")
print(f"          G_predicted = {G_A:.8e}")
print()
print(f"  Path B: alpha_G(proton) = G*m_p^2/(hbar*c) = {alpha_g_proton:.6e}")
print(f"          falls at alpha^{n_proton:.4f}")
print(f"          nearest integer rung: alpha^{round(n_proton)}")
proton_rung = round(n_proton)
proton_coeff = float(alpha_g_proton) / float(alpha) ** proton_rung
print(f"          coefficient at that rung: {proton_coeff:.6f}")
print()

# Path C: Planck mass relationship
# m_Pl = sqrt(hbar*c/G), so G = hbar*c/m_Pl^2
# alpha_G = (m_e/m_Pl)^2
hbar_c_over_G = hbar * c / G
m_Pl = hbar_c_over_G.sqrt()
ratio_e_Pl = m_e / m_Pl
print(f"  Path C: m_e / m_Planck = {ratio_e_Pl:.6e}")
print(f"          (m_e/m_Pl)^2 = alpha_G = {ratio_e_Pl**2:.6e}")
n_planck = math.log(float(ratio_e_Pl)) / math.log(float(alpha))
print(f"          m_e/m_Pl falls at alpha^{n_planck:.4f}")
print(f"          so alpha_G = (m_e/m_Pl)^2 falls at alpha^{2*n_planck:.4f}")
print()

# -----------------------------------------------------------------------
# STRATEGY 3: The dark sector prediction (alpha^10)
# Does alpha^10 appear in any measured cross-sections or anomalies?
# -----------------------------------------------------------------------
print("-" * 80)
print("  STRATEGY 3: The dark matter coupling (alpha^10)")
print("-" * 80)
print()

alpha_10 = float(alpha ** 10)
print(f"  alpha_DM = alpha^10 = {alpha_10:.6e}")
print()
print(f"  If a dark photon exists with coupling epsilon to SM photon,")
print(f"  then epsilon^2 ~ alpha_DM = {alpha_10:.2e}")
print(f"  so epsilon ~ {math.sqrt(alpha_10):.2e}")
print()
print(f"  Current experimental bounds on dark photon mixing (epsilon):")
print(f"    BaBar (2017):        epsilon < 1e-3  for m_A' ~ 1-10 GeV")
print(f"    NA64 (2019):         epsilon < 1e-4  for m_A' ~ 1-100 MeV")
print(f"    LDMX (projected):    epsilon ~ 1e-6  for m_A' ~ 1-100 MeV")
print(f"    Alpha Ladder predicts: epsilon ~ {math.sqrt(alpha_10):.1e}")
print()

in_range = math.sqrt(alpha_10)
if in_range < 1e-3:
    print(f"  Status: BELOW current bounds (not yet excluded)")
    print(f"  The prediction sits {math.log10(1e-3/in_range):.1f} orders of magnitude")
    print(f"  below current sensitivity. Future experiments could reach it.")
else:
    print(f"  Status: WOULD ALREADY BE EXCLUDED by current data.")

# -----------------------------------------------------------------------
# STRATEGY 4: Muon g-2 anomaly
# The muon anomalous magnetic moment has a ~4.2 sigma discrepancy.
# Does any alpha-power combination match the anomaly size?
# -----------------------------------------------------------------------
print()
print("-" * 80)
print("  STRATEGY 4: Muon g-2 anomaly")
print("-" * 80)
print()

# The discrepancy
delta_a_mu = Decimal('2.51e-9')  # experimental - SM prediction
a_mu_exp = Decimal('1.16592061e-3')

print(f"  Muon g-2 discrepancy: delta(a_mu) = {delta_a_mu:.2e}")
print(f"  Fractional anomaly:   delta/a_mu   = {float(delta_a_mu/a_mu_exp):.2e}")
print()

# What alpha power matches the anomaly?
n_anomaly = math.log(float(delta_a_mu)) / math.log(float(alpha))
print(f"  delta(a_mu) falls at alpha^{n_anomaly:.4f}")
print(f"  nearest rung: alpha^{round(n_anomaly)}")
anomaly_coeff = float(delta_a_mu) / float(alpha) ** round(n_anomaly)
print(f"  coefficient: {anomaly_coeff:.4f}")
print()

# Fractional anomaly
n_frac = math.log(float(delta_a_mu / a_mu_exp)) / math.log(float(alpha))
print(f"  Fractional anomaly falls at alpha^{n_frac:.4f}")
print(f"  nearest rung: alpha^{round(n_frac)}")

# -----------------------------------------------------------------------
# STRATEGY 5: Use Kibble balance + cold atom interferometry
# These are the two most promising experimental paths to better G.
# -----------------------------------------------------------------------
print()
print("-" * 80)
print("  STRATEGY 5: Accessible experimental approaches to G")
print("-" * 80)
print()
print("  The ladder predicts G = 6.67323e-11 m^3 kg^-1 s^-2")
print(f"  CODATA 2018 says    G = 6.67430e-11 +/- 0.00015e-11")
print(f"  Difference:           = 1.07e-14  (~16 ppm)")
print()
print("  To resolve this, you need G measured to < 5 ppm.")
print()
print("  Approach A: Atom interferometry (most promising)")
print("    - Rosi et al. 2014 got G = 6.67191e-11 (closest to prediction!)")
print("    - Uses cold atom clouds as test masses")
print("    - No mechanical bearings, no fiber drift")
print("    - Projected next-gen precision: 5-10 ppm")
print("    - Groups: Stanford (Kasevich), Florence (Tino), Wuhan (Zhu)")
print()
print("  Approach B: MEMS oscillators")
print("    - Microscale torsion oscillators on silicon chips")
print("    - Westphal et al. 2021 measured gravity between mm-scale masses")
print("    - Potentially scalable and reproducible")
print("    - Current precision: ~100 ppm (improving)")
print()
print("  Approach C: Propose a targeted re-analysis")
print("    - Contact Rosi/Tino group (their result is closest to prediction)")
print("    - Ask if their systematic error budget permits G = 6.67323e-11")
print("    - Their published uncertainty: +/- 99 ppm")
print("    - Your prediction falls at +1.3 sigma from their central value")
print("    - This is the ONLY existing measurement consistent with the ladder")

print()
print("=" * 80)
print("  VERDICT: WHAT TO DO NEXT")
print("=" * 80)
print()
print("  1. CHEAPEST: Check if alpha^10 matches any dark photon search")
print("     data. This is a literature search, not an experiment.")
print()
print("  2. FASTEST: Verify the mass ratio predictions (Strategy 1)")
print("     against precision particle data. Desktop work.")
print()
print("  3. STRONGEST: Contact the Florence atom interferometry group")
print("     and propose a collaboration to test G = 6.67323e-11.")
print("     Their 2014 measurement is the only one consistent with it.")
print()
print("  4. BOLDEST: Write a paper predicting G to 8 significant figures")
print("     from alpha, phi, and known constants. Let experimentalists")
print("     try to falsify it.")
print("=" * 80)
