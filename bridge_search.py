from decimal import Decimal, getcontext
from itertools import product

getcontext().prec = 50

alpha = Decimal('0.0072973525664')
alpha_g = Decimal('1.7512e-45')
alpha_21 = alpha ** 21

# The exact coefficient we need
target = alpha_g / alpha_21
print(f"Target coefficient: {target}")
print(f"Target (float):     {float(target):.10f}")
print()

# ---------------------------------------------------------------------------
# Search space: combinations of pi, phi, e, 2, 3, ln2 with small powers
# ---------------------------------------------------------------------------
pi = Decimal('3.14159265358979323846264338327950288419716939937510')
phi = (1 + Decimal(5).sqrt()) / 2
e = Decimal('2.71828182845904523536028747135266249775724709369995')
ln2 = Decimal('0.69314718055994530941723212145817656807550013436026')
sqrt2 = Decimal(2).sqrt()
sqrt3 = Decimal(3).sqrt()
sqrt5 = Decimal(5).sqrt()

constants = {
    "pi": pi,
    "phi": phi,
    "e": e,
    "ln2": ln2,
    "sqrt2": sqrt2,
    "sqrt3": sqrt3,
    "sqrt5": sqrt5,
    "2": Decimal(2),
    "3": Decimal(3),
    "4": Decimal(4),
    "5": Decimal(5),
    "6": Decimal(6),
    "7": Decimal(7),
    "1": Decimal(1),
}

target_f = float(target)
results = []

# --- Phase 1: Single constant with rational power p/q ---
print("=== Phase 1: Single constant^(p/q) ===")
for name, c in constants.items():
    c_f = float(c)
    if c_f <= 0:
        continue
    for p in range(-6, 7):
        for q in range(1, 7):
            if p == 0:
                continue
            try:
                val = float(c) ** (p / q)
                if val > 0:
                    err = abs(val - target_f) / target_f * 100
                    if err < 2.0:
                        expr = f"{name}^({p}/{q})" if q != 1 else f"{name}^{p}"
                        results.append((err, expr, val))
            except:
                pass

# --- Phase 2: A * B with rational powers ---
print("=== Phase 2: A^(p/q) * B^(r/s) ===")
const_list = list(constants.items())
powers = [(p, q) for p in range(-4, 5) for q in range(1, 5) if p != 0]

for i, (n1, c1) in enumerate(const_list):
    for j, (n2, c2) in enumerate(const_list):
        if j <= i:
            continue
        c1_f, c2_f = float(c1), float(c2)
        if c1_f <= 0 or c2_f <= 0:
            continue
        for (p1, q1) in powers:
            v1 = c1_f ** (p1 / q1)
            for (p2, q2) in powers:
                try:
                    v2 = c2_f ** (p2 / q2)
                    val = v1 * v2
                    if val > 0:
                        err = abs(val - target_f) / target_f * 100
                        if err < 0.5:
                            e1 = f"{n1}^({p1}/{q1})" if q1 != 1 else f"{n1}^{p1}"
                            e2 = f"{n2}^({p2}/{q2})" if q2 != 1 else f"{n2}^{p2}"
                            expr = f"{e1} * {e2}"
                            results.append((err, expr, val))
                except:
                    pass

# --- Phase 3: Simple fractions a/b * constant^(p/q) ---
print("=== Phase 3: (a/b) * constant^(p/q) ===")
for a in range(1, 13):
    for b in range(1, 13):
        if a == b:
            continue
        frac = a / b
        for name, c in constants.items():
            c_f = float(c)
            if c_f <= 0:
                continue
            for p in range(-4, 5):
                for q in range(1, 5):
                    if p == 0:
                        continue
                    try:
                        val = frac * (c_f ** (p / q))
                        if val > 0:
                            err = abs(val - target_f) / target_f * 100
                            if err < 0.1:
                                ep = f"{name}^({p}/{q})" if q != 1 else f"{name}^{p}"
                                expr = f"({a}/{b}) * {ep}"
                                results.append((err, expr, val))
                    except:
                        pass

# --- Sort and display ---
results.sort(key=lambda x: x[0])
seen = set()
print(f"\n{'='*72}")
print(f"  TOP MATCHES  (target = {target_f:.10f})")
print(f"{'='*72}")
for err, expr, val in results[:40]:
    if expr not in seen:
        seen.add(expr)
        print(f"  {err:8.4f}%  |  {val:.10f}  |  {expr}")
