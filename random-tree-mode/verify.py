#!/usr/bin/env python3
"""Finite certificates for the mode of the independence sequence of a Cayley tree.

Run with python3 verify.py. No third-party packages, network, or input files.
The fixed-point certificate and small-tree checks are rigorous. The final
seeded sample is descriptive. The limiting normal law has an analytic proof.
"""

from collections import Counter
from fractions import Fraction
from heapq import heapify, heappop, heappush
from itertools import product
from math import factorial, prod
from random import Random
from time import monotonic


GRID = 1024
DEPTH = 40
POISSON_CUTOFF = 17
MASS = 1 << 60
PRECISION = 1 << 112
LOG_TERMS = 44
WORD = 16


def check(condition, message):
    if not condition:
        raise ArithmeticError(message)


def ceildiv(a, b):
    return -(-a // b)


def power_bounds(limit):
    """Enclose 2**(-i/GRID) in units of PRECISION, all integer arithmetic."""
    q = PRECISION
    target = q ** GRID
    lo, hi = 0, q
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if 2 * mid ** GRID <= target:
            lo = mid
        else:
            hi = mid
    low, high = [q], [q]
    for _ in range(limit):
        low.append(low[-1] * lo // q)
        high.append(ceildiv(high[-1] * hi, q))
    return low, high


def log1p_bound(x, upper):
    """One directed bound for ln(1+x/Q) in units of Q, for 0 <= x <= Q."""
    q = PRECISION
    divide = ceildiv if upper else int.__floordiv__
    w = divide(q * x, 2 * q + x)
    w2 = divide(w * w, q)
    power = w
    value = 0
    for r in range(LOG_TERMS):
        value += divide(2 * power, 2 * r + 1)
        power = divide(power * w2, q)
    # w_exact <= 1/3.  Tail <= 2 (1/3)^(2R+1)/((2R+1)(1-1/9)).
    tail_units = ceildiv(18 * q, 8 * (2 * LOG_TERMS + 1) * 3 ** (2 * LOG_TERMS + 1))
    return value + (tail_units if upper else 0)


def grid_table(low, high):
    """Integer approximations to y_i=J log2(1+2^(-i/J)), scaled by MASS."""
    ln2_low = log1p_bound(PRECISION, False)
    ln2_high = log1p_bound(PRECISION, True)
    table = []
    for a, b in zip(low, high):
        lower = GRID * MASS * log1p_bound(a, False) // ln2_high
        upper = ceildiv(GRID * MASS * log1p_bound(b, True), ln2_low)
        check(0 <= lower <= GRID * MASS and 0 <= upper - lower <= 2,
              "grid coordinate enclosure")
        table.append(divmod(lower, MASS))
    return table


def normalize(weights, denominator):
    """Floor nonnegative masses on a common denominator, return residual to 0."""
    result = [w // denominator for w in weights]
    remainder = MASS - sum(result)
    check(0 <= remainder < len(result), "normalization or total mass")
    result[0] += remainder
    return result


def convolution(a, b):
    # No carries between digits: each raw coefficient is at most MASS**2.
    check(MASS * MASS < 1 << (8 * WORD), "packed polynomial word too small")
    def pack(values):
        return int.from_bytes(b"".join(x.to_bytes(WORD, "little") for x in values), "little")
    count = len(a) + len(b) - 1
    blob = (pack(a) * pack(b)).to_bytes(WORD * count, "little")
    raw = [int.from_bytes(blob[i:i+WORD], "little") for i in range(0, len(blob), WORD)]
    return normalize(raw, MASS)


def poisson_sum(ordinary):
    """Law of sum of d ordinary messages, mixed over Poisson(1) given d<=K."""
    k = POISSON_CUTOFF
    weights = [factorial(k) // factorial(d) for d in range(k + 1)]
    mixture = [0] * (k * GRID + 1)
    current = [MASS]
    for d, weight in enumerate(weights):
        if d:
            current = convolution(current, ordinary)
        for i, mass in enumerate(current):
            mixture[i] += weight * mass
    return normalize(mixture, sum(weights))


def update_message(sums, table):
    split = [0] * (GRID + 1)
    for mass, (left, frac) in zip(sums, table):
        split[left] += mass * (MASS - frac)
        if frac:
            split[left + 1] += mass * frac
    return normalize(split, MASS)


def decimal_endpoints(lower, upper, places):
    scale = 10 ** places
    left = lower.numerator * scale // lower.denominator
    right = ceildiv(upper.numerator * scale, upper.denominator)
    def format_point(a):
        return f"{a // scale}.{a % scale:0{places}d}"
    return format_point(left), format_point(right)


def density_certificate():
    j, h, k = GRID, DEPTH, POISSON_CUTOFF
    powers_lo, powers_hi = power_bounds((k + 1) * j)
    table = grid_table(powers_lo, powers_hi)
    ordinary = [MASS] + [0] * j
    spine = ordinary.copy()
    for _ in range(h):
        sums = poisson_sum(ordinary)
        spine = update_message(convolution(sums, spine), table)
        ordinary = update_message(sums, table)
    def root_average(powers):
        return Fraction(sum(m * (PRECISION - x) for m, x in zip(spine, powers)),
                        MASS * PRECISION)
    boundary = Fraction(h + 1, 2 ** h)
    interpolation = Fraction(1, 2 * j * j)
    # One update includes at most k convolutions, a mixture, an optional
    # spine convolution, interpolation using a width-2/MASS table, and a
    # final normalization. Each normalization loses <length/MASS in TV.
    update_tv = (Fraction(2, factorial(k + 1))
                 + Fraction((k + 4) * ((k + 1) * j + 2), MASS))
    error = boundary + interpolation + h * (h + 1) // 2 * update_tv
    lower = root_average(powers_hi) - error
    upper = root_average(powers_lo) + error
    check(Fraction(30612549, 10**8) <= lower <= upper <= Fraction(30612646, 10**8),
          "published density enclosure")
    limits = decimal_endpoints(lower, upper, 8)
    print(f"Rigorous c enclosure: [{limits[0]}, {limits[1]}]")
    return lower, upper


def independence_fraction(lower, upper):
    # The unique zero of exp(-x)-x lies between these test points.
    a, b = Fraction(56714329, 10**8), Fraction(56714330, 10**8)
    def partial(x, order):
        return sum((-x)**i / factorial(i) for i in range(order + 1))
    check(partial(a, 17) > a and partial(b, 18) < b, "omega brackets")
    lo, hi = decimal_endpoints(lower / b, upper / a, 6)
    print(f"Rigorous c/omega enclosure: [{lo}, {hi}].")


LEAF = ()
TREE_A = (LEAF, (LEAF, LEAF), (LEAF, ((LEAF,),)))
TREE_B = ((LEAF,), (LEAF,), ((LEAF,), (LEAF,)))


def rooted_adjacency(shape):
    graph = []
    def add(branch, parent):
        v = len(graph)
        graph.append([])
        if parent is not None:
            graph[parent].append(v)
            graph[v].append(parent)
        for child in branch:
            add(child, v)
    add(shape, None)
    return graph


def automorphisms(shape):
    return prod(automorphisms(x) for x in shape) * prod(factorial(c) for c in Counter(shape).values())


def brute_polynomials(graph):
    n = len(graph)
    edges = [(u, v) for u, row in enumerate(graph) for v in row if u < v]
    result = [[0] * (n + 1) for _ in range(2)]
    for mask in range(1 << n):
        if all((mask & (1 << u)) == 0 or (mask & (1 << v)) == 0 for u, v in edges):
            result[mask & 1][bin(mask).count("1")] += 1
    return result


def polynomial(graph):
    """Exact tree recurrence, multiplying polynomials by integer packing."""
    n = len(graph)
    # Every coefficient in a partial product counts subsets of at most n vertices.
    bits = n + 1
    unit = 1 << bits
    parent = [-1] * n
    order = [0]
    for u in order:
        for v in graph[u]:
            if v != parent[u]:
                parent[v] = u
                order.append(v)
    absent, present = [1] * n, [unit] * n
    for v in reversed(order[1:]):
        u = parent[v]
        absent[u] *= absent[v] + present[v]
        present[u] *= absent[v]
    encoded = absent[0] + present[0]
    coeff = []
    while encoded:
        coeff.append(encoded & (unit - 1))
        encoded >>= bits
    return coeff


def gadget_certificate():
    pairs = []
    for shape, expected in [(TREE_A, (130, 40, 395, 140)), (TREE_B, (117, 36, 336, 120))]:
        graph = rooted_adjacency(shape)
        counts = brute_polynomials(graph)
        z0, z1 = map(sum, counts)
        d0, d1 = (sum(i * a for i, a in enumerate(row)) for row in counts)
        check(len(graph) == 10 and (z0, z1, d0, d1) == expected, "fringe gadget")
        pairs.append((Fraction(z1, z0), Fraction(d0, z0), Fraction(d1, z1)))
    a, b = pairs
    check(a[0] == b[0] and a[1]-b[1] == a[2]-b[2] == Fraction(1, 6), "slot identity")
    aut = automorphisms(TREE_A), automorphisms(TREE_B)
    check(aut == (2, 4), "automorphism orders")
    # e < 3 since 1/j! <= 1/2^(j-1) for j>=2, strictly for j>=4.
    lower = Fraction(1, 36 * sum(aut) * 3**10)
    upper = 4 * 6 * 4**2
    print(f"Rigorous variance bounds: {lower} < v <= {upper}. Gadget shift: 1/6.")


def decode_pruefer(code):
    n = len(code) + 2
    degrees = [1] * n
    for v in code:
        degrees[v] += 1
    leaves = [i for i, d in enumerate(degrees) if d == 1]
    heapify(leaves)
    graph = [[] for _ in range(n)]
    def edge(u, v):
        graph[u].append(v)
        graph[v].append(u)
    for v in code:
        u = heappop(leaves)
        edge(u, v)
        degrees[v] -= 1
        if degrees[v] == 1:
            heappush(leaves, v)
    edge(*leaves)
    return graph


def ball_shape(graph, v, parent, depth):
    if depth == 0:
        return ()
    return tuple(sorted(ball_shape(graph, u, v, depth - 1) for u in graph[v] if u != parent))


def ball_parameters(shape, depth):
    if depth == 0:
        return 1, 1
    child_data = [ball_parameters(s, depth - 1) for s in shape]
    return 1 + sum(m for m, _ in child_data), sum(b for _, b in child_data)


def check_ball_law(counts, n, depth):
    for shape, count in counts.items():
        m, b = ball_parameters(shape, depth)
        a = automorphisms(shape)
        labellings = factorial(n-1) // factorial(n-m)
        completions = 1 if m == n else b * (n-m+b)**(n-m-1)
        check(count * a == labellings * completions, "rooted ball law")


def finite_check():
    total = 0
    for n in range(2, 7):
        balls = {h: Counter() for h in (1, 2)}
        for code in product(range(n), repeat=n-2):
            tree = decode_pruefer(code)
            first, second = brute_polynomials(tree)
            coeff = polynomial(tree)
            check(coeff + [0] * (n+1-len(coeff)) == [a+b for a, b in zip(first, second)],
                  "polynomial versus subset enumeration")
            for h, counts in balls.items():
                counts[ball_shape(tree, 0, -1, h)] += 1
            total += 1
        for h, counts in balls.items():
            check_ball_law(counts, n, h)
    check(total == 1441, "number of labelled trees tested")
    print(f"Rigorous finite checks: polynomials and rooted ball laws on all {total} labelled trees of orders 2 through 6.")


def sample():
    rng = Random(260920961)
    for n in (64, 128, 256):
        modes, means = [], []
        for _ in range(48):
            tree = decode_pruefer([rng.randrange(n) for _ in range(n-2)])
            coeff = polynomial(tree)
            peak = max(coeff)
            modes.append(max(i for i, a in enumerate(coeff) if a == peak))
            means.append(Fraction(sum(i * a for i, a in enumerate(coeff)), sum(coeff)))
        mbar = Fraction(sum(modes), len(modes))
        variance = sum((m-mbar)**2 for m in modes) / (len(modes)-1)
        gap = max(abs(m-u) for m, u in zip(modes, means))
        print(f"Non-rigorous sample n={n}, 48 trees: mean(M/n)={float(mbar/n):.5f}, "
              f"sample Var(M)/n={float(variance/n):.5f}, max |M-mu|={float(gap):.4f}")


if __name__ == "__main__":
    start = monotonic()
    gadget_certificate()
    finite_check()
    independence_fraction(*density_certificate())
    sample()
    print(f"Elapsed: {monotonic()-start:.1f} seconds.")
