from decimal import Decimal, getcontext

# Set high precision for deep physics calculations
getcontext().prec = 50

# ---------------------------------------------------------------------------
# Foundational Constants
# ---------------------------------------------------------------------------
# Fine-structure constant: alpha = e^2 / (4*pi*epsilon_0*hbar*c)
alpha = Decimal('0.0072973525664')

# Gravitational coupling constant (electron)
alpha_g = Decimal('1.7512e-45')

# High-precision pi
pi = Decimal('3.14159265358979323846264338327950288419716939937510')

# Golden ratio phi = (1 + sqrt(5)) / 2
phi = (1 + Decimal(5).sqrt()) / 2

# Euler's number
e = Decimal('2.71828182845904523536028747135266249775724709369995')

# ---------------------------------------------------------------------------
# Electron Geometry Scaling
#   r_e = alpha * lambda_bar_c = alpha^2 * a_0
# ---------------------------------------------------------------------------
# Reduced Compton wavelength (m)
lambda_bar_c = Decimal('3.8615926764e-13')
# Bohr radius (m)
a_0 = Decimal('5.2917721067e-11')

r_e_from_compton = alpha * lambda_bar_c        # alpha * lambda_bar_c
r_e_from_bohr = alpha ** 2 * a_0               # alpha^2 * a_0
r_e_nist = Decimal('2.8179403227e-15')          # NIST value


# ---------------------------------------------------------------------------
# Geometric Ladder
# ---------------------------------------------------------------------------
def calculate_geometric_rungs():
    rungs = []
    for n in range(1, 25):
        val = alpha ** n
        ratio = val / alpha_g
        label = ""
        if n == 10:
            label = "Dark Matter Candidate (The Blank)"
        elif n == 21:
            label = "Gravity Floor (The Bridge)"
        rungs.append({
            "power": n,
            "value": val,
            "ratio_to_gravity": ratio,
            "label": label,
        })

    alpha_21 = alpha ** 21

    # Bridge candidates for alpha_G = C * alpha^21
    bridges = {
        "phi * alpha^21  (original)":      phi * alpha_21,
        "(phi^2 / 2) * alpha^21":          (phi ** 2 / 2) * alpha_21,
        "(5/12) * pi * alpha^21":          (Decimal(5) / Decimal(12)) * pi * alpha_21,
        "sqrt(e) / cbrt(2) * alpha^21":    (e.sqrt() / Decimal(2) ** (Decimal(1) / Decimal(3))) * alpha_21,
    }

    # Dark sector candidate: alpha_DM ~ alpha^10
    alpha_dm = alpha ** 10

    # Pi-Bridge
    pi_bridge = alpha_21 * (pi ** 2)

    return rungs, bridges, alpha_dm, pi_bridge


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    rungs_data, bridges, alpha_dm, bridge_val = calculate_geometric_rungs()

    print("=" * 72)
    print("  ALPHA LADDER  --  Geometric Powers of the Fine-Structure Constant")
    print("=" * 72)

    # --- Electron geometry verification ---
    print("\n--- Electron Geometry Scaling ---")
    print(f"  r_e = alpha * lambda_bar_c  = {r_e_from_compton:.6e} m")
    print(f"  r_e = alpha^2 * a_0         = {r_e_from_bohr:.6e} m")
    print(f"  r_e (NIST)                  = {r_e_nist:.6e} m")

    # --- The 42-order gap ---
    gap = alpha / alpha_g
    print("\n--- The 42-Order Gap (The Great Desert) ---")
    print(f"  alpha / alpha_G = {gap:.2e}")

    # --- Full ladder ---
    print("\n--- Geometric Rungs ---")
    for r in rungs_data:
        tag = f"  << {r['label']}" if r['label'] else ""
        print(
            f"  Alpha^{r['power']:2d}: {r['value']:.6e}  |  "
            f"Ratio to alpha_G: {float(r['ratio_to_gravity']):.6e}{tag}"
        )

    # --- Dark sector ---
    print("\n--- The Dark Sector ---")
    print(f"  alpha_DM  = alpha^10        = {alpha_dm:.6e}")

    # --- Bridge candidates ---
    print("\n--- The Geometric Bridge (Candidates) ---")
    print(f"  alpha_G (measured)          = {alpha_g:.6e}")
    print()
    for name, pred in bridges.items():
        residual = abs(alpha_g - pred) / alpha_g * 100
        print(f"  {name:<36s} = {pred:.6e}  |  residual: {residual:.4f}%")

    # --- Pi-Bridge ---
    print(f"\n  Pi-Bridge (alpha^21 * pi^2) = {bridge_val:.6e}")
    print("=" * 72)
