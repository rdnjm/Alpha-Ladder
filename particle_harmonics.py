import math
from decimal import Decimal, getcontext

getcontext().prec = 50

# The Alpha Constant from your notes
alpha = Decimal('0.0072973525693')
inv_alpha = 1 / float(alpha)

# Standard Model Masses (in MeV/c^2)
# Using electron as the base unit
masses = {
    "Electron (Base)": 0.511,
    "Muon": 105.66,
    "Tau": 1776.86,
    "Up Quark": 2.2,
    "Down Quark": 4.7,
    "Strange Quark": 95.0,
    "Charm Quark": 1275.0,
    "Bottom Quark": 4180.0,
    "Top Quark": 172760.0,
    "W Boson": 80377.0,
    "Z Boson": 91187.0,
    "Higgs Boson": 125100.0,
    "Proton": 938.27
}

def scan_particle_harmonics():
    print(f"{'Particle':<18} | {'Mass Ratio':<12} | {'Alpha Rung (n)':<15} | {'Closeness'}")
    print("-" * 65)

    m_electron = masses["Electron (Base)"]

    for name, m in masses.items():
        ratio = m / m_electron
        # Equation: Ratio = (1/alpha)^n  => n = log(Ratio) / log(1/alpha)
        n = math.log(ratio) / math.log(inv_alpha)

        # Check for near-integer or half-integer rungs
        nearest_half = round(n * 2) / 2
        closeness = abs(n - nearest_half)

        note = "MATCH!" if closeness < 0.05 else ""

        print(f"{name:<18} | {ratio:<12.2f} | {n:<15.4f} | {note}")

scan_particle_harmonics()
