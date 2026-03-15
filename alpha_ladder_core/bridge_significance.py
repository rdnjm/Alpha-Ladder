"""
Bridge significance analysis -- look-elsewhere effect for phi^2/2.

Quantifies whether the match between phi^2/2 and alpha_g/alpha^21
(~160 ppm residual found among 135,816 candidates) is statistically
significant, or could happen by chance in a large search.

Two complementary approaches:
1. Coverage analysis: what fraction of the target interval is "close"
   to some expression? (This is the empirical p-value.)
2. Algebraic analysis: phi^2/2 lives in Q(sqrt(5)), the same field
   forced by the vacuum polynomial.
"""

import bisect
import math
import random

from alpha_ladder_core.bridge_search import _get_math_constants, compute_target_coefficient


# ---------------------------------------------------------------------------
# Simple label formatting (no Unicode, plain ASCII)
# ---------------------------------------------------------------------------

def _plain_power(name, p, q):
    """Format a constant raised to a rational power as a plain string."""
    if q == 1:
        return f"{name}^({p})"
    return f"{name}^({p}/{q})"


# ---------------------------------------------------------------------------
# Enumeration
# ---------------------------------------------------------------------------

def enumerate_all_expression_values():
    """Evaluate all candidates from the 3-phase search space.

    No filtering by error -- every finite positive evaluation is collected.

    Returns dict with keys:
        values, unique_values, labels, phase_counts, total_count,
        phi_expressions, phi_count
    """
    constants_dict = _get_math_constants()
    const_list = list(constants_dict.items())
    powers_phase1 = [(p, q) for p in range(-6, 7) for q in range(1, 7) if p != 0]
    powers_phase23 = [(p, q) for p in range(-4, 5) for q in range(1, 5) if p != 0]

    values = []
    labels = []
    phase_counts = {"phase1": 0, "phase2": 0, "phase3": 0}

    # Phase 1: single constant c^(p/q)
    for name, c in const_list:
        c_f = float(c)
        if c_f <= 0:
            continue
        for p, q in powers_phase1:
            try:
                val = c_f ** (p / q)
                if math.isfinite(val) and val > 0:
                    label = _plain_power(name, p, q)
                    values.append(val)
                    labels.append(label)
                    phase_counts["phase1"] += 1
            except Exception:
                pass

    # Phase 2: A^(p/q) * B^(r/s), i < j
    for i, (n1, c1) in enumerate(const_list):
        for j, (n2, c2) in enumerate(const_list):
            if j <= i:
                continue
            c1_f, c2_f = float(c1), float(c2)
            if c1_f <= 0 or c2_f <= 0:
                continue
            for p1, q1 in powers_phase23:
                try:
                    v1 = c1_f ** (p1 / q1)
                except Exception:
                    continue
                for p2, q2 in powers_phase23:
                    try:
                        v2 = c2_f ** (p2 / q2)
                        val = v1 * v2
                        if math.isfinite(val) and val > 0:
                            l1 = _plain_power(n1, p1, q1)
                            l2 = _plain_power(n2, p2, q2)
                            label = f"{l1}*{l2}"
                            values.append(val)
                            labels.append(label)
                            phase_counts["phase2"] += 1
                    except Exception:
                        pass

    # Phase 3: (a/b) * c^(p/q)
    for a in range(1, 13):
        for b in range(1, 13):
            if a == b:
                continue
            frac = a / b
            for name, c in const_list:
                c_f = float(c)
                if c_f <= 0:
                    continue
                for p, q in powers_phase23:
                    try:
                        val = frac * (c_f ** (p / q))
                        if math.isfinite(val) and val > 0:
                            label = f"({a}/{b})*{_plain_power(name, p, q)}"
                            values.append(val)
                            labels.append(label)
                            phase_counts["phase3"] += 1
                    except Exception:
                        pass

    total_count = len(values)

    # Phi expressions
    phi_expressions = [
        (values[i], labels[i])
        for i in range(total_count)
        if "phi" in labels[i]
    ]
    phi_count = len(phi_expressions)

    # Unique values (sorted)
    unique_values = sorted(set(values))

    return {
        "values": values,
        "unique_values": unique_values,
        "labels": labels,
        "phase_counts": phase_counts,
        "total_count": total_count,
        "phi_expressions": phi_expressions,
        "phi_count": phi_count,
    }


# ---------------------------------------------------------------------------
# Coverage analysis
# ---------------------------------------------------------------------------

def compute_coverage(expression_values, interval_min=0.5, interval_max=2.0, epsilon=0.00016):
    """Compute what fraction of [interval_min, interval_max] is within
    relative epsilon of any expression value.

    Parameters
    ----------
    expression_values : list of float
        The unique_values list from enumerate_all_expression_values().
    interval_min, interval_max : float
        Bounds of the target interval.
    epsilon : float
        Relative tolerance (160 ppm = 0.00016).

    Returns dict with coverage_fraction, covered_length, values_in_interval,
    interval_length.
    """
    interval_length = interval_max - interval_min

    # Filter to values in the interval
    in_interval = [v for v in expression_values if interval_min <= v <= interval_max]
    in_interval.sort()

    if not in_interval:
        return {
            "coverage_fraction": 0.0,
            "covered_length": 0.0,
            "values_in_interval": 0,
            "interval_length": interval_length,
        }

    # Build bands and clamp
    bands = []
    for v in in_interval:
        lo = max(v * (1 - epsilon), interval_min)
        hi = min(v * (1 + epsilon), interval_max)
        bands.append((lo, hi))

    # Sort by start
    bands.sort(key=lambda x: x[0])

    # Merge overlapping intervals
    merged = [bands[0]]
    for lo, hi in bands[1:]:
        if lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))

    covered_length = sum(hi - lo for lo, hi in merged)
    coverage_fraction = covered_length / interval_length

    return {
        "coverage_fraction": coverage_fraction,
        "covered_length": covered_length,
        "values_in_interval": len(in_interval),
        "interval_length": interval_length,
    }


# ---------------------------------------------------------------------------
# Bonferroni correction
# ---------------------------------------------------------------------------

def compute_bonferroni_pvalue(n_trials, epsilon=0.00016):
    """Bonferroni upper bound on the p-value.

    p_single = 2 * epsilon (band width relative to unit interval),
    p_corrected = min(1.0, n_trials * p_single).

    Returns dict with p_single, p_corrected, significant_at_005.
    """
    p_single = 2 * epsilon
    p_corrected = min(1.0, n_trials * p_single)
    return {
        "p_single": p_single,
        "p_corrected": p_corrected,
        "significant_at_005": p_corrected < 0.05,
    }


# ---------------------------------------------------------------------------
# Empirical p-value
# ---------------------------------------------------------------------------

def compute_empirical_pvalue(expression_values, target_value=None,
                             interval_min=0.5, interval_max=2.0, epsilon=0.00016):
    """Empirical p-value equals the coverage fraction: probability that a
    uniformly random point in the interval is within epsilon of some expression.

    Parameters
    ----------
    expression_values : list of float
        The unique_values list.
    target_value : float or None
        Approximate target (default 1.3).
    interval_min, interval_max : float
        Interval bounds.
    epsilon : float
        Relative tolerance.

    Returns dict with p_empirical, p_bonferroni, interpretation.
    """
    if target_value is None:
        target_value = 1.3

    coverage = compute_coverage(expression_values, interval_min, interval_max, epsilon)
    p_empirical = coverage["coverage_fraction"]

    bonferroni = compute_bonferroni_pvalue(len(expression_values), epsilon)
    p_bonferroni = bonferroni["p_corrected"]

    if p_empirical >= 0.05:
        interpretation = (
            "The numerical match alone is not statistically significant "
            "after accounting for the look-elsewhere effect."
        )
    else:
        interpretation = (
            "The numerical match is statistically significant "
            "even after the look-elsewhere correction."
        )

    return {
        "p_empirical": p_empirical,
        "p_bonferroni": p_bonferroni,
        "interpretation": interpretation,
    }


# ---------------------------------------------------------------------------
# Monte Carlo validation
# ---------------------------------------------------------------------------

def run_monte_carlo_validation(expression_values, interval_min=0.5, interval_max=2.0,
                               epsilon=0.00016, n_samples=10000, seed=42):
    """Monte Carlo estimate of the coverage fraction.

    Draws n_samples uniform random targets and checks whether each falls
    within relative epsilon of any expression value using bisect for
    O(log N) lookup.

    Returns dict with p_monte_carlo, p_monte_carlo_uncertainty, consistent,
    n_samples, n_hits.
    """
    sorted_vals = sorted(v for v in expression_values
                         if interval_min <= v <= interval_max)

    rng = random.Random(seed)
    n_hits = 0

    for _ in range(n_samples):
        target = rng.uniform(interval_min, interval_max)
        idx = bisect.bisect_left(sorted_vals, target)

        hit = False
        # Check the nearest candidates (idx-1 and idx)
        for candidate_idx in (idx - 1, idx):
            if 0 <= candidate_idx < len(sorted_vals):
                v = sorted_vals[candidate_idx]
                if abs(v - target) / target <= epsilon:
                    hit = True
                    break
        if hit:
            n_hits += 1

    p_mc = n_hits / n_samples
    uncertainty = math.sqrt(p_mc * (1 - p_mc) / n_samples) if n_samples > 0 else 0.0

    # Compute coverage for consistency check
    coverage = compute_coverage(expression_values, interval_min, interval_max, epsilon)
    coverage_frac = coverage["coverage_fraction"]

    consistent = abs(p_mc - coverage_frac) < 3 * uncertainty if uncertainty > 0 else True

    return {
        "p_monte_carlo": p_mc,
        "p_monte_carlo_uncertainty": uncertainty,
        "consistent": consistent,
        "n_samples": n_samples,
        "n_hits": n_hits,
    }


# ---------------------------------------------------------------------------
# Phi subset analysis
# ---------------------------------------------------------------------------

def analyze_phi_subset(expression_values_data, target_value=None):
    """Restrict coverage analysis to phi-involving expressions only.

    Parameters
    ----------
    expression_values_data : dict
        Full output from enumerate_all_expression_values().
    target_value : float or None
        Approximate target (default 1.3).

    Returns dict with phi_count, phi_coverage_fraction, phi_p_empirical,
    best_phi_match, algebraic_field_note.
    """
    if target_value is None:
        target_value = 1.3

    phi_expressions = expression_values_data["phi_expressions"]
    phi_count = len(phi_expressions)

    phi_values = sorted(set(v for v, _ in phi_expressions))

    coverage = compute_coverage(phi_values)
    phi_coverage_fraction = coverage["coverage_fraction"]

    # Best match to target
    best_match = None
    best_residual = float("inf")
    for val, label in phi_expressions:
        residual = abs(val - target_value) / target_value * 1e6
        if residual < best_residual:
            best_residual = residual
            best_match = {"value": val, "label": label, "residual_ppm": residual}

    if best_match is None:
        best_match = {"value": 0.0, "label": "", "residual_ppm": float("inf")}

    algebraic_field_note = (
        "phi^2/2 lives in Q(sqrt(5)), the same algebraic field forced by "
        "the vacuum polynomial x^2 + Dx + d = 0 with D = -(1+sqrt(5))/2."
    )

    return {
        "phi_count": phi_count,
        "phi_coverage_fraction": phi_coverage_fraction,
        "phi_p_empirical": phi_coverage_fraction,
        "best_phi_match": best_match,
        "algebraic_field_note": algebraic_field_note,
    }


# ---------------------------------------------------------------------------
# Expression density histogram
# ---------------------------------------------------------------------------

def compute_expression_density(expression_values, n_bins=100,
                               interval_min=0.5, interval_max=2.0):
    """Histogram of expression values in the interval.

    Returns dict with bin_centers, counts, density.
    """
    interval_length = interval_max - interval_min
    bin_width = interval_length / n_bins

    counts = [0] * n_bins
    total_in_interval = 0

    for v in expression_values:
        if interval_min <= v <= interval_max:
            idx = int((v - interval_min) / bin_width)
            if idx == n_bins:
                idx = n_bins - 1
            counts[idx] += 1
            total_in_interval += 1

    bin_centers = [interval_min + (i + 0.5) * bin_width for i in range(n_bins)]

    if total_in_interval > 0 and bin_width > 0:
        density = [c / total_in_interval / bin_width for c in counts]
    else:
        density = [0.0] * n_bins

    return {
        "bin_centers": bin_centers,
        "counts": counts,
        "density": density,
    }


# ---------------------------------------------------------------------------
# Main summary
# ---------------------------------------------------------------------------

def summarize_bridge_significance(constants):
    """Main entry point -- run the full significance analysis pipeline.

    Parameters
    ----------
    constants : SimpleNamespace
        CODATA constants (must have alpha and alpha_g attributes).

    Returns dict with search_space, target_value, residual_ppm, coverage,
    bonferroni, empirical_pvalue, monte_carlo, phi_analysis, histogram,
    honest_assessment.
    """
    _, target_f = compute_target_coefficient(constants)

    expr_data = enumerate_all_expression_values()
    unique_values = expr_data["unique_values"]
    total_count = expr_data["total_count"]
    unique_count = len(unique_values)

    # Find best match residual
    best_residual = float("inf")
    for v in unique_values:
        residual = abs(v - target_f) / target_f * 1e6
        if residual < best_residual:
            best_residual = residual

    coverage = compute_coverage(unique_values)
    bonferroni = compute_bonferroni_pvalue(total_count)
    empirical = compute_empirical_pvalue(unique_values, target_value=target_f)
    mc = run_monte_carlo_validation(unique_values)
    phi_analysis = analyze_phi_subset(expr_data, target_value=target_f)
    histogram = compute_expression_density(unique_values)

    coverage_pct = coverage["coverage_fraction"] * 100
    p_emp = empirical["p_empirical"]
    is_or_not = "is not" if p_emp >= 0.05 else "is"

    honest_assessment = (
        f"The bridge search evaluated {total_count} expressions, of which "
        f"{unique_count} unique values cover {coverage_pct:.1f}% of the "
        f"target interval [0.5, 2.0] within 160 ppm. "
        f"The empirical p-value (probability of a random target matching "
        f"by chance) is {p_emp:.4f}, meaning the numerical coincidence "
        f"alone {is_or_not} surprising. "
        f"The deeper significance lies in the algebraic connection: "
        f"phi^2/2 belongs to Q(sqrt(5)), the same field forced by the "
        f"vacuum polynomial. This algebraic link cannot be captured by "
        f"the look-elsewhere correction."
    )

    return {
        "search_space": {
            "total_count": total_count,
            "unique_count": unique_count,
            "phase_counts": expr_data["phase_counts"],
        },
        "target_value": target_f,
        "residual_ppm": best_residual,
        "coverage": coverage,
        "bonferroni": bonferroni,
        "empirical_pvalue": empirical,
        "monte_carlo": mc,
        "phi_analysis": phi_analysis,
        "histogram": histogram,
        "honest_assessment": honest_assessment,
    }
