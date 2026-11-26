"""Reproduce the finite checks in the random-graph stability note. Python 3.8+.

Run: python3 verify_random_stability.py
All decisions use integers or fractions. Decimal output is for display only.
"""

from fractions import Fraction
from functools import lru_cache
from hashlib import shake_256
from itertools import combinations
from math import comb


def moments(n, s, q):
    """First and second moments of the number of independent s-sets; q = 1 - p."""
    m = comb(n, s)
    edges = s * (s - 1) // 2
    mean = m * q**edges
    second = sum(
        m * comb(s, i) * comb(n - s, s - i) * q**(2 * edges - i * (i - 1) // 2)
        for i in range(max(0, 2 * s - n), s + 1)
    )
    return mean, second


def check_moments_by_enumeration(n, q):
    """Average over every labelled graph, weighting each by its exact probability."""
    pairs = list(combinations(range(n), 2))
    required_absences = []
    for size in range(n + 1):
        masks = []
        for vertices in combinations(range(n), size):
            masks.append(sum(1 << i for i, (u, v) in enumerate(pairs)
                             if u in vertices and v in vertices))
        required_absences.append(masks)

    e = len(pairs)
    denom = q.denominator**e
    weight = [(q.denominator - q.numerator)**j * q.numerator**(e - j)
              for j in range(e + 1)]
    first = [0] * (n + 1)
    second = [0] * (n + 1)
    for graph in range(1 << e):
        graph_weight = weight[bin(graph).count("1")]
        for s, masks in enumerate(required_absences):
            x = sum((graph & mask) == 0 for mask in masks)
            first[s] += graph_weight * x
            second[s] += graph_weight * x * x
    for s in range(n + 1):
        observed = (Fraction(first[s], denom), Fraction(second[s], denom))
        if observed != moments(n, s, q):
            raise ArithmeticError(f"Moment mismatch at n={n}, s={s}, q={q}")


def example_graph():
    n = 86
    size = (comb(n, 2) + 7) // 8
    bits = shake_256(b"BC:86:622").digest(size)
    adjacency = [0] * n
    for index, (u, v) in enumerate(combinations(range(n), 2)):
        if bits[index // 8] & (1 << (index % 8)):
            adjacency[u] |= 1 << v
            adjacency[v] |= 1 << u
    return adjacency


def enumerate_independent_sets(adjacency):
    """List each independent set in increasing vertex order."""
    counts = [0] * (len(adjacency) + 1)

    def visit(candidates, size):
        counts[size] += 1
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            vertex = bit.bit_length() - 1
            visit(candidates & ~adjacency[vertex], size + 1)

    visit((1 << len(adjacency)) - 1, 0)
    while counts[-1] == 0:
        counts.pop()
    return counts


def coefficients_by_vertex_deletion(adjacency):
    """Evaluate I(B) by deleting a vertex or its closed neighborhood, then decode."""
    n = len(adjacency)
    base = 1 << (n + 1)  # Every coefficient is at most 2**n.

    @lru_cache(None)
    def evaluate(vertices):
        if not vertices:
            return 1
        v = vertices.bit_length() - 1
        rest = vertices ^ (1 << v)
        return evaluate(rest) + base * evaluate(rest & ~adjacency[v])

    value = evaluate((1 << n) - 1)
    coeffs = []
    while value:
        value, digit = divmod(value, base)
        coeffs.append(digit)
    return coeffs


def check_root_disk(coeffs):
    """Certify a root within 1/100 of 11/5 + 9i by Rouche's theorem."""
    # Write z = (11 + 45i + w)/5. The disk in w has radius 1/20.
    d = len(coeffs) - 1
    powers = [(1, 0)]
    for _ in range(d):
        x, y = powers[-1]
        powers.append((11 * x - 45 * y, 45 * x + 11 * y))

    # Taylor coefficients b_j of 5**d I((11 + 45i + w)/5).
    taylor = []
    for j in range(d + 1):
        real = imag = 0
        for i in range(j, d + 1):
            factor = coeffs[i] * 5**(d - i) * comb(i, j)
            x, y = powers[i - j]
            real += factor * x
            imag += factor * y
        taylor.append((real, imag))
    radius = Fraction(1, 20)
    linear_lower = radius * max(map(abs, taylor[1]))
    rest_upper = sum((abs(x) + abs(y)) * radius**j
                     for j, (x, y) in enumerate(taylor) if j != 1)
    if not rest_upper < Fraction(87, 100) * linear_lower:
        raise ArithmeticError("Root disk certificate failed")
    return rest_upper / linear_lower


def finite_probability_certificate():
    n, k, q = 260_000, 30, Fraction(1, 2)
    delta = Fraction(1, 8)
    stats = {s: moments(n, s, q) for s in range(k - 3, k + 2)}
    mean = {s: x for s, (x, _) in stats.items()}
    ratio = mean[k - 1] * mean[k - 2] / mean[k - 3]
    ratio_bound = ratio * (1 + delta)**2 / (1 - delta)
    if not ratio_bound < Fraction(74, 100):
        raise ArithmeticError("Four-coefficient test failed")

    # Bonferroni for X_k > 0, followed by Markov and three Chebyshev bounds.
    hit = mean[k] - (stats[k][1] - mean[k]) / 2
    error = sum((stats[s][1] / mean[s]**2 - 1) / delta**2
                for s in range(k - 3, k))
    lower_bound = hit - mean[k + 1] - error
    if not lower_bound > Fraction(1016, 10000):
        raise ArithmeticError("Finite probability certificate failed")
    return ratio_bound, lower_bound


if __name__ == "__main__":
    for n, q in [(5, Fraction(1, 3)), (5, Fraction(2, 3)), (6, Fraction(1, 2))]:
        check_moments_by_enumeration(n, q)
        print(f"Verified all first and second moments on {2**comb(n, 2):,} "
              f"labelled graphs for n={n}, q={q}.")

    graph = example_graph()
    coeffs = enumerate_independent_sets(graph)
    expected = [1, 86, 1810, 12478, 31706, 31559, 12792, 2225, 205, 21, 2]
    if coeffs != expected or coefficients_by_vertex_deletion(graph) != expected:
        raise ArithmeticError("The two graph counts must reproduce the certificate")
    determinant = coeffs[-2] * coeffs[-3] - coeffs[-1] * coeffs[-4]
    if determinant != -145:
        raise ArithmeticError("Unexpected Hurwitz determinant")
    print(f"86-vertex graph, coefficients in ascending order: {coeffs}")
    print(f"Second Hurwitz determinant: {determinant}")
    disk_ratio = check_root_disk(coeffs)
    print(f"Root in |z - (11/5 + 9i)| < 1/100; comparison ratio {float(disk_ratio):.8f}")
    ratio_bound, probability = finite_probability_certificate()
    print(f"n=260,000: coefficient ratio upper bound {float(ratio_bound):.8f}; "
          f"instability probability lower bound {float(probability):.8f}")
