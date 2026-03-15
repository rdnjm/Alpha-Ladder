"""
Brute force search: what rung spacing best fits the Standard Model mass spectrum?

Instead of assuming integer or half-integer rungs, test every rational
spacing 1/k for k=1..24 and score how well the particles land on those rungs.
Also test irrational spacings involving phi, pi, e, ln2.
"""

import math
from decimal import Decimal, getcontext

getcontext().prec = 50

alpha = Decimal('0.0072973525693')
inv_alpha = 1 / float(alpha)
phi_f = (1 + math.sqrt(5)) / 2

# Standard Model masses in MeV/c^2 (PDG 2023 central values)
# Exclude electron (trivial base) and quarks with >10% uncertainty
masses_precise = {
    "Muon":         105.6583755,    # +/- 0.0000023
    "Tau":          1776.86,        # +/- 0.12
    "Proton":       938.272088,     # +/- 0.000058
    "W Boson":      80377.0,        # +/- 12
    "Z Boson":      91187.6,        # +/- 2.1
    "Higgs Boson":  125100.0,       # +/- 140  (using 125.10 GeV)
    "Top Quark":    172760.0,       # +/- 300
}

masses_uncertain = {
    "Up Quark":     2.16,           # +/- 0.49  (23% uncertainty!)
    "Down Quark":   4.67,           # +/- 0.48  (10%)
    "Strange":      93.4,           # +/- 8.6   (9%)
    "Charm":        1270.0,         # +/- 20    (1.6%)
    "Bottom":       4180.0,         # +/- 30    (0.7%)
}

all_masses = {**masses_precise, **masses_uncertain}
m_e = 0.51099895000  # electron mass MeV

# Compute all rungs (log ratios)
rungs = {}
for name, m in all_masses.items():
    rungs[name] = math.log(m / m_e) / math.log(inv_alpha)


# -----------------------------------------------------------------------
# SEARCH 1: Rational spacings 1/k
# -----------------------------------------------------------------------
def score_spacing(spacing, rung_dict, tolerance_frac=0.15):
    """Score how well particles land on rungs with given spacing.
    Returns (avg_closeness, n_matches, details).
    Closeness = distance to nearest multiple of spacing, as fraction of spacing.
    """
    details = []
    for name, n in rung_dict.items():
        nearest_k = round(n / spacing)
        nearest_rung = nearest_k * spacing
        delta = abs(n - nearest_rung)
        frac = delta / spacing  # 0 = perfect, 0.5 = worst
        match = frac < tolerance_frac
        details.append((name, n, nearest_rung, nearest_k, delta, frac, match))

    n_matches = sum(1 for d in details if d[6])
    avg_frac = sum(d[5] for d in details) / len(details)
    return avg_frac, n_matches, details


print("=" * 80)
print("  SEARCH 1: Best rational rung spacing (1/k)")
print("=" * 80)
print()
print("  Testing spacings 1/1, 1/2, 1/3, ..., 1/24")
print("  Scoring: what fraction of particles land within 15% of a rung?")
print()

results = []
for k in range(1, 25):
    spacing = 1.0 / k
    avg, matches, details = score_spacing(spacing, rungs)
    results.append((k, spacing, avg, matches, details))

# Also compute expected random matches for each k
# If rungs are spaced by s, a random value has probability 2*0.15*s/s = 0.30
# of landing within 15% of a rung. With N particles, expected = 0.30 * N.
N = len(rungs)
expected_random = 0.30 * N

results.sort(key=lambda x: -x[3])  # sort by number of matches

print(f"  {'k':>3s}  {'Spacing':>8s}  {'Matches':>7s}  {'Expected':>8s}  {'Excess':>6s}  {'Avg Err':>7s}")
print(f"  {'-'*3}  {'-'*8}  {'-'*7}  {'-'*8}  {'-'*6}  {'-'*7}")
for k, spacing, avg, matches, details in results[:15]:
    excess = matches - expected_random
    print(f"  {k:3d}  {spacing:8.4f}  {matches:7d}  {expected_random:8.1f}  {excess:+6.1f}  {avg:7.4f}")

# -----------------------------------------------------------------------
# SEARCH 2: Irrational spacings
# -----------------------------------------------------------------------
print()
print("=" * 80)
print("  SEARCH 2: Irrational rung spacings")
print("=" * 80)
print()

irrational_spacings = {
    "1/phi":            1.0 / phi_f,
    "1/pi":             1.0 / math.pi,
    "1/e":              1.0 / math.e,
    "ln2":              math.log(2),
    "1/phi^2":          1.0 / phi_f**2,
    "phi/pi":           phi_f / math.pi,
    "2/pi":             2.0 / math.pi,
    "1/sqrt(2)":        1.0 / math.sqrt(2),
    "1/sqrt(3)":        1.0 / math.sqrt(3),
    "phi - 1 = 1/phi":  phi_f - 1,
    "ln(phi)":          math.log(phi_f),
    "pi/6":             math.pi / 6,
    "1/(2*phi)":        1.0 / (2*phi_f),
    "1/(3*phi)":        1.0 / (3*phi_f),
    "phi/4":            phi_f / 4,
}

irr_results = []
for name, spacing in irrational_spacings.items():
    avg, matches, details = score_spacing(spacing, rungs)
    irr_results.append((name, spacing, avg, matches, details))

irr_results.sort(key=lambda x: -x[3])

print(f"  {'Spacing':<20s}  {'Value':>8s}  {'Matches':>7s}  {'Expected':>8s}  {'Excess':>6s}  {'Avg Err':>7s}")
print(f"  {'-'*20}  {'-'*8}  {'-'*7}  {'-'*8}  {'-'*6}  {'-'*7}")
for name, spacing, avg, matches, details in irr_results:
    excess = matches - expected_random
    print(f"  {name:<20s}  {spacing:8.4f}  {matches:7d}  {expected_random:8.1f}  {excess:+6.1f}  {avg:7.4f}")


# -----------------------------------------------------------------------
# SEARCH 3: Continuous optimization
# What spacing minimizes the total distance-to-nearest-rung?
# -----------------------------------------------------------------------
print()
print("=" * 80)
print("  SEARCH 3: Continuous optimization (brute force)")
print("=" * 80)
print()

rung_values = list(rungs.values())

best_score = float('inf')
best_spacing = 0
scores = []

# Scan from 0.01 to 2.0 in steps of 0.001
for i in range(10, 2001):
    s = i / 1000.0
    total = 0
    for n in rung_values:
        nearest = round(n / s) * s
        total += (n - nearest) ** 2
    rms = math.sqrt(total / len(rung_values))
    scores.append((s, rms))
    if rms < best_score:
        best_score = rms
        best_spacing = s

# Find top 10 local minima
from itertools import islice
local_mins = []
for i in range(1, len(scores) - 1):
    if scores[i][1] < scores[i-1][1] and scores[i][1] < scores[i+1][1]:
        local_mins.append(scores[i])

local_mins.sort(key=lambda x: x[1])

print(f"  Global best: spacing = {best_spacing:.4f}, RMS error = {best_score:.6f}")
print()
print(f"  Top 10 local minima:")
print(f"  {'Spacing':>8s}  {'RMS Err':>8s}  {'Nearest known constant':<30s}")
print(f"  {'-'*8}  {'-'*8}  {'-'*30}")

known = {
    1/phi_f: "1/phi",
    1/math.pi: "1/pi",
    1/math.e: "1/e",
    0.5: "1/2",
    1/3: "1/3",
    0.25: "1/4",
    1/6: "1/6",
    math.log(2): "ln(2)",
    phi_f - 1: "phi - 1 = 1/phi",
    1/math.sqrt(2): "1/sqrt(2)",
    math.log(phi_f): "ln(phi)",
    math.pi/6: "pi/6",
    2/math.pi: "2/pi",
    phi_f/4: "phi/4",
    1/(2*phi_f): "1/(2*phi)",
}

for s, rms in local_mins[:10]:
    nearest_const = min(known.items(), key=lambda kv: abs(kv[0] - s))
    delta_pct = abs(nearest_const[0] - s) / s * 100
    label = f"{nearest_const[1]} ({delta_pct:.1f}% off)" if delta_pct < 5 else "---"
    print(f"  {s:8.4f}  {rms:8.6f}  {label}")


# -----------------------------------------------------------------------
# Show the winning spacing applied to all particles
# -----------------------------------------------------------------------
print()
print("=" * 80)
print(f"  BEST FIT: spacing = {best_spacing:.4f}")
print("=" * 80)
print()

avg, matches, details = score_spacing(best_spacing, rungs)
details.sort(key=lambda x: x[1])

print(f"  {'Particle':<15s}  {'n':>7s}  {'Nearest Rung':>12s}  {'k':>4s}  {'Delta':>7s}  {'Match?':>6s}")
print(f"  {'-'*15}  {'-'*7}  {'-'*12}  {'-'*4}  {'-'*7}  {'-'*6}")
for name, n, nearest, k, delta, frac, match in details:
    tag = "YES" if match else ""
    print(f"  {name:<15s}  {n:7.4f}  {nearest:12.4f}  {k:4.0f}  {delta:7.4f}  {tag:>6s}")
