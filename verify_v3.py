"""Certificates and a small labelled simulation for the v3 manuscript. Python 3.8+."""

import argparse
from fractions import Fraction as F
from hashlib import shake_256
from math import comb, factorial, isqrt


# Integer interval endpoints in units of 10**-34. All divisions round outwards.
SCALE = 10**34


def ceildiv(n, d):
    return -((-n) // d)


class Box:
    def __init__(self, lo, hi):
        assert lo <= hi
        self.lo, self.hi = lo, hi

    @classmethod
    def rational(cls, x, y=None):
        x, y = F(x), F(x if y is None else y)
        return cls(
            x.numerator * SCALE // x.denominator,
            ceildiv(y.numerator * SCALE, y.denominator),
        )

    def __add__(self, other):
        return Box(self.lo + other.lo, self.hi + other.hi)

    def __sub__(self, other):
        return Box(self.lo - other.hi, self.hi - other.lo)

    def __mul__(self, other):
        p = [a * b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Box(min(p) // SCALE, ceildiv(max(p), SCALE))

    def __pow__(self, k):
        assert self.lo >= 0 and k >= 0
        if k == 0:
            return Box(SCALE, SCALE)
        d = SCALE ** (k - 1)
        return Box(self.lo**k // d, ceildiv(self.hi**k, d))


def cells(lo, hi, n):
    lo, hi = F(lo), F(hi)
    return [
        Box.rational(lo + (hi - lo) * j / n, lo + (hi - lo) * (j + 1) / n)
        for j in range(n)
    ]


def b(j):
    return j * (j - 1) // 2


def normalization_checks():
    """Check the reciprocal coefficient products and the padding correction exactly."""
    q = F(7, 10)
    for n in range(3, 18):
        means = [comb(n, j) * q ** b(j) for j in range(n + 1)]
        for k in range(2, n + 1):
            ratio = means[k] / means[k - 1]
            product = F(1)
            for j in range(1, k + 1):
                h = j - 1
                product *= q**h * F(k - h, k) * F(n - k + 1, n - k + 1 + h)
                assert ratio**j * means[k - j] / means[k] == product
                assert 0 < product <= q ** b(j)
    # K4 has I(w)=1+4w. Padding at k=2 supplies a single extra z.
    mean_two = 6 * q
    ratio = mean_two / 4
    for z in (F(1, 2), F(-2), -4 / ratio):
        direct = (ratio * z) ** 2 * (1 + 4 / (ratio * z)) / mean_two
        assert direct == z * (4 * ratio + ratio**2 * z) / mean_two
    print(
        "RIGOROUS: finite reciprocal identities, coefficient bounds and padded origin."
    )


def make_coeff(s, order, degree=15):
    return [
        (s ** b(j)) * Box.rational((-1) ** j * factorial(j) // factorial(j - order))
        for j in range(order, degree + 1)
    ]


def enclosure(coeff, x):
    value = coeff[-1]
    for term in reversed(coeff[:-1]):
        value = value * x + term
    return value + Box(-1, 1)


def check_tail():
    # For x<=7, s<=5/16 and orders r<=2, the first omitted term at j=16
    # is bounded by 16^2 * 7^16 * (5/16)^120. Subsequent ratios are <1/2.
    s = F(5, 16)
    assert 2 * 16**2 * 7**16 * s**120 < F(1, SCALE)
    assert 7 * s**16 * F(17, 15) < F(1, 2)


def rectangle_checks():
    check_tail()
    q0, q1 = F(17, 25), F(18, 25)
    for q in cells(q0, q1, 120):
        c = make_coeff(q**4, 0)
        a = Box.rational(F(10, 3))
        assert enclosure(c, a).hi < 0 and enclosure(c, a * q**2).hi < 0
    z = Box.rational
    for s in cells(q1**4, F(31, 100), 80):
        f, df, ddf = (make_coeff(s, r) for r in range(3))
        left, right = z("1.6"), z("3.5")
        assert enclosure(f, left).lo > 0 and enclosure(f, right).lo > 0
        assert enclosure(df, left).hi < 0 and enclosure(df, right).lo > 0
        for x in cells("1.6", "3.5", 100):
            d2 = enclosure(ddf, x)
            assert d2.lo > SCALE // 25
            assert (x * d2 + enclosure(df, x)).lo > 0
    for s in cells(F(29, 40) ** 4, F(31, 100), 100):
        df = make_coeff(s, 1)
        for x in cells(3, F(13, 2), 130):
            assert enclosure(df, x).lo > 0
    # Constants in the all-degree annulus estimates, paper (9)–(10).
    s = F(5, 16)
    root_upper = F(5591, 10000)
    assert root_upper**2 > s
    v = F(20, 19)
    # A polynomial lower bound for the central Fourier terms, decreasing in s.
    assert 4 * s * root_upper < 1
    assert (
        1
        - (v + 1 / v) * root_upper
        + 2 * s * s
        - 2 * v**3 * s**4 * root_upper / (1 - v * s**3 * root_upper)
    ) > F(1, 20)
    assert F(1, 20) * (1 - s * (1 + v + 1 / v)) > s**6 * v**3 / (1 - s**4 * v)
    assert (1 + v + 1 / v) * s**2 < 1
    assert F(1, 20) * F(19, 20) ** 3 > s**4 * (2 + root_upper)
    assert F(19, 20) > s**3
    assert 1 - 2 + 4 * F(27, 100) - 8 * F(27, 100) ** 3 + 16 * F(27, 100) ** 6 < 0
    assert 4 - 24 * F(27, 100) ** 2 > 0
    print(
        "RIGOROUS: compact parameter boxes, derivative signs and all-degree tail margins."
    )


def exact_series(s, x, derivative=0, tol=F(1, 10**24)):
    """Sum f_s or a derivative at a nonnegative rational x, with a tail enclosure."""
    total = F(0)
    term = F((-1) ** derivative * factorial(derivative)) * s ** b(derivative)
    j = derivative
    qs = s**j
    while True:
        total += term
        ratio = x * qs * F(j + 1, j + 1 - derivative)
        nxt = -term * ratio
        if ratio <= F(1, 2) and abs(nxt) < tol:
            return min(total, total + nxt), max(total, total + nxt)
        j += 1
        qs *= s
        term = nxt


def sign(s, x, derivative=0):
    tol = F(1, 10**22)
    for _ in range(8):
        lower, upper = exact_series(s, x, derivative, tol)
        if lower > 0:
            return 1
        if upper < 0:
            return -1
        tol *= tol
    raise ArithmeticError("Sign needs more precision")


def bisect(s, lo, hi, derivative=0, steps=45):
    lo, hi = F(lo), F(hi)
    direction = sign(s, lo, derivative)
    assert sign(s, hi, derivative) == -direction
    for _ in range(steps):
        mid = (lo + hi) / 2
        if sign(s, mid, derivative) == direction:
            lo = mid
        else:
            hi = mid
    return lo, hi


def interval_series(slo, shi, lo, hi):
    """Interval f_s(x) by signed monomials; the entire omitted tail is bounded."""
    a = d = F(0)
    for j in range(2000):
        lower = slo ** b(j) * lo**j
        upper = shi ** b(j) * hi**j
        if j % 2:
            a -= upper
            d -= lower
        else:
            a += lower
            d += upper
        ratio = shi**j * hi
        if ratio <= F(1, 2) and upper * ratio < F(1, 10**20):
            return a - 2 * upper * ratio, d + 2 * upper * ratio
    raise ArithmeticError("Increase the term budget")


def roots_pair(q):
    s = q**4
    if s <= F(27, 100):
        mid = F(2)
    else:
        lo, hi = bisect(s, F("1.6"), F("3.5"), 1, 52)
        mid = (lo + hi) / 2
    assert sign(s, mid) < 0
    left = bisect(s, 0, mid)
    # When mid=2, f_s(2/s)=1-(2/s)f_s(2)>0. Stay inside the two-zero disk.
    right_end = min(2 / s, 1 / q**6) if mid == 2 else 1 / q**6
    right = bisect(s, mid, right_end)
    return left, right


def heights(q):
    roots = roots_pair(q)
    crossings = []
    for left, right in roots:
        lower, upper = interval_series(q**4, q**4, left / q**2, right / q**2)
        crossings.append((1 - upper, 1 - lower))
    return crossings


def exp_negative(x):
    """Rational enclosure of exp(-x), by exp(-x/2^h) and repeated squaring."""
    x = F(x)
    h = 0
    while x > 1:
        x /= 2
        h += 1
    v = term = F(1)
    for j in range(1, 31):
        term *= -x / j
        v += term
    lo = v + term * (-x) / 31
    lo, hi = round_out((lo, v), 120)
    for _ in range(h):
        lo, hi = round_out((lo * lo, hi * hi), 120)
    return max(lo, F(0)), hi


def inside(interval, lo, hi):
    assert F(lo) < interval[0] <= interval[1] < F(hi), (lo, hi)


def transition_checks():
    # The first two odd zeros coalesce at Q. A tangent bound encloses the minimum.
    ql, qh = F(".74572240"), F(".74572242")
    center = F("2.32038")
    assert sign(ql**4, center) < 0
    val = exact_series(qh**4, center)[0]
    sl, sh = exact_series(qh**4, center, 1)
    assert val - max(abs(sl), abs(sh)) ** 2 / (2 * F(1, 25)) > 0
    assert F(".3092493") < ql**4 < qh**4 < F(".3092494")
    assert F(".25427758") <= 1 - qh < 1 - ql <= F(".25427760")
    for q in (ql, qh):
        assert sign(q**4, F("2.32036"), 1) < 0 and sign(q**4, F("2.32040"), 1) > 0
    for x, expected in ((F("2.32036"), -1), (F("2.32040"), 1)):
        enc = enclosure(make_coeff(Box.rational(ql**4, qh**4), 1), Box.rational(x))
        assert enc.hi < 0 if expected < 0 else enc.lo > 0
    a, d = interval_series(ql**4, qh**4, F("2.32036") / qh**2, F("2.32040") / ql**2)
    endpoint = (1 - d, 1 - a)
    inside(endpoint, ".6816", ".6820")
    # At q_c the second low odd zero is exactly the first divided by q^2.
    for q, side in ((F(".726"), -1), (F(".727"), 1)):
        first, second = roots_pair(q)
        assert first[1] / q**2 < second[0] if side < 0 else first[0] / q**2 > second[1]
    assert F("1.68") / qh**2 > 3 and F("3.3") / F(".725") ** 2 < F(13, 2)
    r1, r2 = roots_pair(F(".725"))
    assert r1[0] > F("1.68") and r2[1] < F("3.3")
    dlow, dhigh = F(".6823"), F(".6824")
    assert 1 - dlow - dlow**3 > 0 > 1 - dhigh - dhigh**3
    print("RIGOROUS: s_* in (.3092493,.3092494); 1-Q in (.25427758,.25427760).")
    print("RIGOROUS: 1-q_c in (.273,.274); root of 1-q-q^3 in (.6823,.6824).")
    print(
        "RIGOROUS: q_c in (.726,.727), Q in (.74572240,.74572242), t_* in (.6816,.6820)."
    )
    return endpoint


def round_out(bounds, digits=12):
    scale = 10**digits
    lo, hi = bounds
    return F((lo * scale).__floor__(), scale), F((hi * scale).__ceil__(), scale)


def probability_below_band(a, bounds, limit=140):
    """Enclose the largest stable mass in the mixed phase, including tied atoms."""
    al, au = round_out(a)
    bl, bu = round_out(bounds)
    if not 0 < al <= au < bl <= bu < 1:
        # The outward-rounded parameter intervals can overlap extremely close to Q.
        return F(0), F(1)
    # c(b) is increasing on (0,1).
    c_upper = round_out((F(0), exp_negative(1 / bu - 1)[1] / bu))[1]
    best_lo = best_hi = F(0)
    for m in range(1, limit + 1):
        # Poisson masses decrease with the mean throughout [m/bu,m/bl].
        def tailpart(mu, first):
            first = max(first, 0)
            subtotal = sum((mu**j / factorial(j) for j in range(first, m + 1)), F(0))
            e = exp_negative(mu)
            return round_out((e[0] * subtotal, e[1] * subtotal))

        first_lo = (al * m / bu).__ceil__()
        first_hi = (au * m / bl).__ceil__()
        best_lo = max(best_lo, tailpart(F(m) / bl, first_hi)[0])
        best_hi = max(best_hi, tailpart(F(m) / bu, first_lo)[1])
        omitted = min(c_upper ** (m + 1), tailpart(F(m) / bu, 0)[1])
        if omitted < best_lo:
            return best_lo, best_hi
    return best_lo, min(F(1), max(best_hi, omitted))


def show(label, bounds, digits=7):
    lo, hi = round_out(bounds, digits)
    scale = 10**digits

    def fmt(x):
        z = int(x * scale)
        return (
            ("-" if z < 0 else "")
            + str(abs(z) // scale)
            + "."
            + str(abs(z) % scale).zfill(digits)
        )

    print(f"RIGOROUS: {label} in [{fmt(lo)}, {fmt(hi)}].")


def root_frequency(pair, qmin, qmax):
    """Outward square roots in 1e-10 units for y = sqrt(x/q^3)."""
    scale = 10**10
    lo = pair[0] / qmax**3
    hi = pair[1] / qmin**3
    return (
        F(isqrt((lo * scale**2).__floor__()), scale),
        F(1 + isqrt((hi * scale**2).__floor__()), scale),
    )


def upper_probability(band):
    bl, bu = round_out(band)
    assert bl > 0
    return 1 - exp_negative(1 / bu)[1], 1 - exp_negative(1 / bl)[0]


def examples(endpoint):
    from verify_stability_followup import crossing_bounds

    crossing_bounds()
    print(
        "RIGOROUS: inherited q=1/2 crossing and peak reproduced with rational arithmetic."
    )
    ql, qh = F(".726"), F(".727")
    low, high = roots_pair(ql), roots_pair(qh)
    for j, bracket in enumerate(((low[0][0], high[0][1]), (high[1][0], low[1][1])), 1):
        show(f"x_{j} at q_c", bracket, 5)
        show(f"imaginary frequency {j} at q_c", root_frequency(bracket, ql, qh), 5)
    show(
        "nonzero imaginary frequency at Q",
        root_frequency((F("2.32036"), F("2.32040")), F(".74572240"), F(".74572242")),
        5,
    )
    show("limsup at q_c", upper_probability((F(1), F(1))))
    for q in (F(".7"), F(".74")):
        hi, lo = heights(q)
        print(f"q = {q}")
        for j, pair in enumerate(roots_pair(q), 1):
            show(f"x_{j}", pair)
            show(f"imaginary frequency {j}", root_frequency(pair, q, q))
        show("first and upper crossing", hi)
        show("second crossing", lo)
        if lo[1] < 0:
            inside(hi, "1.1946219", "1.1946220")
            prob = upper_probability(hi)
            inside(prob, ".5670291", ".5670292")
            show("limsup", prob)
        else:
            inside(hi, ".8600525", ".8600526")
            inside(lo, ".4164890", ".4164891")
            m, M = probability_below_band(lo, hi)
            inside((1 - M, 1 - m), ".5084363", ".5084364")
            show("liminf", (1 - M, 1 - m))
    tl, tu = round_out(endpoint)
    x, y = 1 / tu, 1 / tl
    show("eta at Q", (x, y), 5)
    prob = (1 - x * exp_negative(x)[1], 1 - y * exp_negative(y)[0])
    inside(prob, ".6616", ".6617")
    show("liminf at Q", prob, 4)


def requested_parameter(q):
    """Additional point certificates at an exact rational q in (0,1)."""
    if not 0 < q < 1:
        raise ValueError("Expected 0 < q < 1")
    print(f"Additional rational parameter q = {q}")
    if q >= F(".74572242"):
        print("RIGOROUS: above Q; first-instability constant 0, probability limit 1.")
        return
    if q >= F(".74572240"):
        s = q**4
        left, right = bisect(s, F("1.6"), F("3.5"), 1, 65)
        mid = (left + right) / 2
        val = exact_series(s, mid)
        dl, du = exact_series(s, mid, 1)
        if val[0] - max(abs(dl), abs(du)) ** 2 / (2 * F(1, 25)) > 0:
            print("RIGOROUS: above Q by strict convexity; probability limit 1.")
            return
        if val[1] >= 0:
            print(
                "The rational input is unresolved at the current spectral-bracket precision."
            )
            return
    hi, lo = heights(q)
    for j, pair in enumerate(roots_pair(q), 1):
        show(f"x_{j}", pair)
        show(f"imaginary frequency {j}", root_frequency(pair, q, q))
    show("upper crossing B", hi)
    show("lower crossing A before restriction to t >= 0", lo)
    if lo[1] < 0 and hi[0] > 1:
        show("first-instability constant", hi)
        show("upper probability endpoint, lower endpoint 0", upper_probability(hi))
    elif lo[0] > 0 and hi[1] < 1:
        ml, mu = probability_below_band(lo, hi)
        show("lower probability endpoint, upper endpoint 1", (1 - mu, 1 - ml))
        print("RIGOROUS: first-instability constant 0.")
    else:
        print("The crossing intervals straddle a regime boundary at this precision.")


def monte_carlo():
    from verify_random_stability import (
        coefficients_by_vertex_deletion,
        enumerate_independent_sets,
    )
    from verify_stability_followup import hurwitz_minors

    def sample(n, threshold, first_seed, count_graph):
        size = 128
        unstable = singular = 0
        for seed in range(first_seed, first_seed + size):
            raw = shake_256(f"v3:{n}:{threshold}:{seed}".encode()).digest(comb(n, 2))
            adj = [0] * n
            ix = 0
            for u in range(n):
                for v in range(u + 1, n):
                    if raw[ix] < threshold:
                        adj[u] |= 1 << v
                        adj[v] |= 1 << u
                    ix += 1
            coeff = count_graph(adj)
            for minor in hurwitz_minors(coeff[::-1]):
                if minor <= 0:
                    singular += minor == 0
                    unstable += minor < 0
                    break
        print(
            f"NON-RIGOROUS Monte Carlo: p={F(threshold, 256)}, n={n}, {size} samples, "
            f"{unstable} unstable, {singular} inconclusive, {size - unstable - singular} strictly stable."
        )
        return unstable, singular

    for threshold, observed in ((68, 2), (76, 6)):
        assert sample(36, threshold, 0, enumerate_independent_sets) == (0, 0)
        assert sample(72, threshold, 1, coefficients_by_vertex_deletion) == (
            observed,
            0,
        )


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit(
            "Run without -O so all rational certificate assertions execute."
        )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-simulation", action="store_true")
    parser.add_argument(
        "--q", type=F, help="also enclose the crossings at an exact rational q"
    )
    args = parser.parse_args()
    normalization_checks()
    rectangle_checks()
    endpoint = transition_checks()
    examples(endpoint)
    if args.q is not None:
        requested_parameter(args.q)
    if not args.skip_simulation:
        monte_carlo()
