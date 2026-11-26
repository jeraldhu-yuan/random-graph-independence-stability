"""Exact checks for the follow-up note. Run beside verify_random_stability.py.

python3 -B verify_stability_followup.py
python3 -B verify_stability_followup.py --sample-86 8192
"""

import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import shake_256
from itertools import combinations
from math import comb, factorial

from verify_random_stability import coefficients_by_vertex_deletion, enumerate_independent_sets


def triangular(j):
    return j * (j - 1) // 2


def odd_bounds(v):
    # On [8.5, 8.6] the alternating tail decreases from the first power.
    upper = sum((-v)**m / F(2**triangular(2 * m + 1)) for m in range(13))
    return upper - v**13 / 2**triangular(27), upper


def crossing_bounds():
    lo, hi = F(17, 2), F(43, 5)
    assert odd_bounds(lo)[0] > 0 and odd_bounds(hi)[1] < 0
    for _ in range(48):
        mid = (lo + hi) / 2
        neg, pos = odd_bounds(mid)
        if neg > 0:
            lo = mid
        elif pos < 0:
            hi = mid
        else:
            raise ArithmeticError("Refine the alternating-series enclosure")

    # For T(xi), even partial sums bound from below and odd ones from above.
    # Evaluate positive and negative terms at opposite interval endpoints.
    tau_lo = tau_hi = F(0)
    for m in range(1, 12):
        denom = 2**triangular(2 * m)
        if m % 2:
            if m <= 10:
                tau_lo += lo**m / denom
            tau_hi += hi**m / denom
        else:
            tau_lo -= hi**m / denom
            tau_hi -= lo**m / denom

    def exp_complement_sum(x, terms):
        return sum((-1)**(m + 1) * x**m / factorial(m) for m in range(1, terms + 1))

    prob_lo = exp_complement_sum(1 / tau_hi, 14)
    prob_hi = exp_complement_sum(1 / tau_lo, 15)
    assert F("8.57160089651") < lo < hi < F("8.57160089652")
    assert F("3.15699422295") < tau_lo < tau_hi < F("3.15699422296")
    assert F("0.27149223459") < prob_lo < prob_hi < F("0.27149223460")
    rho = F(99, 50)
    assert 9 * F(99, 100) * (1 - 1 / rho - 1 / rho**3) > F(16, 5)
    for v in (9, 81):
        assert 1 - F(v, 8) + F(v * v, 1024) < 0
    return lo, tau_lo, prob_lo


def hurwitz_minors(descending):
    """Yield exact leading determinants by fraction-free elimination."""
    d = len(descending) - 1
    matrix = [[descending[2 * col - row + 1] if 0 <= 2 * col - row + 1 <= d else 0
               for col in range(d)] for row in range(d)]
    previous = 1
    for col in range(d):
        pivot = matrix[col][col]
        yield pivot
        if pivot == 0:
            return
        for i in range(col + 1, d):
            for j in range(col + 1, d):
                numerator = pivot * matrix[i][j] - matrix[i][col] * matrix[col][j]
                quotient, remainder = divmod(numerator, previous)
                if remainder:
                    raise ArithmeticError("Inexact Bareiss division")
                matrix[i][j] = quotient
        previous = pivot


def right_half_plane_count(descending):
    """Count roots by Routh's rule. A zero pivot needs a separate limit."""
    degree = len(descending) - 1
    upper = [F(a) for a in descending[::2]]
    lower = [F(a) for a in descending[1::2]]
    lower += [F(0)] * (len(upper) - len(lower))
    column = [upper[0], lower[0]]
    for _ in range(2, degree + 1):
        if lower[0] == 0:
            return None
        ratio = upper[0] / lower[0]
        row = [a - ratio * b for a, b in zip(upper[1:], lower[1:])] + [F(0)]
        upper, lower = lower, row
        column.append(row[0])
    if 0 in column:
        return None
    return sum((x < 0) != (y < 0) for x, y in zip(column, column[1:]))


def check_truncated_series():
    d = 16
    cases = [("3.156994", None), ("3.156995", 6), ("3.17", 4),
             ("3.25", 3), ("4.1", 2)]
    for value, expected in cases:
        t = F(value)
        scale = t.denominator * 2**triangular(d)
        # Read as descending powers to test the reciprocal polynomial.
        coeffs = [int(t * scale)] + [scale // 2**triangular(j) for j in range(1, d + 1)]
        minors = list(hurwitz_minors(coeffs))
        assert len(minors) == d and 0 not in minors
        first_negative = next((j for j, x in enumerate(minors, 1) if x < 0), None)
        assert first_negative == expected, (value, first_negative)
    print("Full degree-16 Hurwitz test passed on both sides of the crossing.")


def finite_upper_bound():
    n, d = 65_600_000, 44
    upper = F(comb(n, d + 1), 2**triangular(d + 1))
    for s in range(2, d + 1):
        weighted = sum(comb(s, i) * comb(n - s, s - i) * (2**triangular(i) - 1)
                       for i in range(2, s + 1))
        upper += 256 * F(weighted, comb(n, s))
    assert 2 * F(15, 17)**2 > F(3, 2)
    assert upper < F("0.009215")
    return upper


PAIRS_86 = tuple(combinations(range(86), 2))


def graph_86(seed):
    bits = shake_256(f"BC:86:{seed}".encode("ascii")).digest((len(PAIRS_86) + 7) // 8)
    adjacency = [0] * 86
    for index, (u, v) in enumerate(PAIRS_86):
        if bits[index // 8] & (1 << (index % 8)):
            adjacency[u] |= 1 << v
            adjacency[v] |= 1 << u
    return adjacency


def check_graph():
    graph = graph_86(6)
    counts = enumerate_independent_sets(graph)
    assert counts == coefficients_by_vertex_deletion(graph)
    assert counts == [1, 86, 1839, 13105, 35078, 37685, 17134, 3495, 385, 39, 3]
    minors = list(hurwitz_minors(counts[::-1]))
    assert minors[1:3] == [4530, -5819319]
    assert right_half_plane_count(counts[::-1]) == 2
    print(f"86-vertex seed 6: second determinant {minors[1]}, third {minors[2]}.")


def sample_86(num_graphs):
    first_negative = Counter()
    root_counts = Counter()
    undetermined = 0
    for seed in range(num_graphs):
        counts = enumerate_independent_sets(graph_86(seed))
        descending = counts[::-1]
        for index, determinant in enumerate(hurwitz_minors(descending), 1):
            if determinant < 0:
                first_negative[index] += 1
                count = right_half_plane_count(descending)
                if count == 0:
                    raise ArithmeticError("Hurwitz and Routh classifications disagree")
                root_counts[count] += 1
                break
            if determinant == 0:
                undetermined += 1
                break
    print(f"86 vertices, SHAKE256 seeds 0 through {num_graphs - 1}: "
          f"{sum(first_negative.values())} unstable, {first_negative[2]} flagged by the four-coefficient test, "
          f"{undetermined} undetermined at a zero pivot.")
    for index, count in sorted(first_negative.items()):
        print(f"First negative Hurwitz determinant at order {index}: {count} graphs.")
    for count in sorted(k for k in root_counts if k is not None):
        print(f"Unstable graphs with {count} right-half-plane roots: {root_counts[count]}.")
    if root_counts[None]:
        print(f"Unstable graphs with a singular Routh table: {root_counts[None]}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-86", type=int, default=0, metavar="COUNT")
    args = parser.parse_args()
    if args.sample_86 < 0:
        parser.error("COUNT must be nonnegative")
    xi, tau, probability = crossing_bounds()
    print(f"Certified enclosures give xi = {float(xi):.10f}, tau = {float(tau):.10f}, "
          f"limsup = {float(probability):.10f} (rounded).")
    check_truncated_series()
    print(f"n=65,600,000: certified instability upper bound {float(finite_upper_bound()):.10f} (rounded).")
    check_graph()
    if args.sample_86:
        sample_86(args.sample_86)
