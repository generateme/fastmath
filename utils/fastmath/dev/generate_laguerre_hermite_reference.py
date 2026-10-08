"""One-off generator of independent reference values for the generalized Laguerre and the Hermite
polynomials of `fastmath.polynomials`: `eval-laguerre-L`, `eval-hermite-H`, `eval-hermite-He` and their
`*-ratio` and polynomial object forms.

Produces `test/resources/polynomials/laguerre_hermite_reference.edn`. Everything is exact rational
arithmetic (`fractions.Fraction`); the order and the arguments are doubles taken at their exact binary
value. The polynomials come from the explicit sums (not from the recurrences the library uses):

  L(n; a)(x) = sum_k (-1)^k C(n+a, n-k) x^k / k!            (generalized binomial coefficient)
  H(n)(x)    = n! sum_m (-1)^m (2x)^(n-2m) / (m! (n-2m)!)
  He(n)(x)   = n! sum_m (-1)^m x^(n-2m) / (m! (n-2m)! 2^m)

and the generator asserts that they satisfy the three term recurrences. The values and the first
derivatives are evaluated at the exact binary value of every double argument, then rounded once to a
double. Rows whose value or derivative is not a finite double are left out.

Layout (ratios are `[numerator denominator]` pairs, coefficients ascending in the power of x):
  - `:laguerre`    `[{:order a :decimal-exact? bool :coefficients [[ratio ...] ...] :grid [[n x value derivative] ...]} ...]`
  - `:hermite-h`   `{:coefficients [[ratio ...] ...] :grid [...]}`
  - `:hermite-he`  `{:coefficients [[ratio ...] ...] :grid [...]}`
`:decimal-exact?` is true when `rationalize` of the order gives its exact binary value (for example 0.5
or 2.5, not 1 + 2^-30): only then the `*-ratio` forms are compared with `:coefficients`.

Run with (from repo root, using the `uv`-managed Python env mentioned in AGENTS.md):
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_laguerre_hermite_reference.py
"""
from fractions import Fraction
from math import factorial

OUT = "/home/ts/clojure/fastmath/test/resources/polynomials/laguerre_hermite_reference.edn"
DBL_MAX = 1.7976931348623157e308


def to_double(q):
    try:
        v = float(Fraction(q))
    except OverflowError:
        return None
    return v if abs(v) <= DBL_MAX else None


def fmt(v):
    return repr(float(v))


def ratio(q):
    q = Fraction(q)
    return "[%d %d]" % (q.numerator, q.denominator)


def ratios(qs):
    return "[" + " ".join(ratio(q) for q in qs) + "]"


def horner(cs, x):
    acc = Fraction(0)
    for c in reversed(cs):
        acc = acc * x + c
    return acc


def derivative(cs):
    return [i * c for i, c in enumerate(cs)][1:] or [Fraction(0)]


def gbinom(top, k):
    acc = Fraction(1)
    for j in range(k):
        acc = acc * (top - j) / (j + 1)
    return acc


def laguerre(n, a):
    return [(-1) ** k * gbinom(n + a, n - k) / factorial(k) for k in range(n + 1)]


def hermite_h(n):
    cs = [Fraction(0)] * (n + 1)
    for m in range(n // 2 + 1):
        cs[n - 2 * m] = Fraction((-1) ** m * factorial(n) * 2 ** (n - 2 * m), factorial(m) * factorial(n - 2 * m))
    return cs


def hermite_he(n):
    cs = [Fraction(0)] * (n + 1)
    for m in range(n // 2 + 1):
        cs[n - 2 * m] = Fraction((-1) ** m * factorial(n), factorial(m) * factorial(n - 2 * m) * 2 ** m)
    return cs


def shift(cs):
    return [Fraction(0)] + cs


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def scale(cs, s):
    return [c * s for c in cs]


def check_recurrences():
    """The explicit sums satisfy the three term recurrences (a consistency check of the generator)."""
    for n in range(2, 25):
        for a in (Fraction(0), Fraction(1, 2), Fraction(-3, 2), Fraction(5, 2)):
            # n L_n = (2n - 1 + a - x) L_(n-1) - (n - 1 + a) L_(n-2)
            lhs = scale(laguerre(n, a), n)
            rhs = padd(padd(scale(laguerre(n - 1, a), 2 * n - 1 + a), scale(shift(laguerre(n - 1, a)), -1)),
                       scale(laguerre(n - 2, a), -(n - 1 + a)))
            assert padd(lhs, scale(rhs, -1)) == [0] * len(padd(lhs, scale(rhs, -1))), ("laguerre", n, a)
        # H_n = 2 x H_(n-1) - 2 (n - 1) H_(n-2);  He_n = x He_(n-1) - (n - 1) He_(n-2)
        d = padd(hermite_h(n), scale(padd(scale(shift(hermite_h(n - 1)), 2), scale(hermite_h(n - 2), -2 * (n - 1))), -1))
        assert all(c == 0 for c in d), ("H", n)
        d = padd(hermite_he(n), scale(padd(shift(hermite_he(n - 1)), scale(hermite_he(n - 2), -(n - 1))), -1))
        assert all(c == 0 for c in d), ("He", n)


def decimal_exact(p):
    return Fraction(str(p)) == Fraction(p)


def grid_rows(polys, xs):
    rows = []
    for n, poly in polys:
        ds = derivative(poly)
        for x in xs:
            xf = Fraction(x)
            v, d = to_double(horner(poly, xf)), to_double(horner(ds, xf))
            if v is None or d is None:
                continue
            rows.append("[%d %s %s %s]" % (n, fmt(x), fmt(v), fmt(d)))
    return rows


def symmetric(xs):
    return sorted(set(xs) | {-x for x in xs})


LAGUERRE_XS = symmetric([0.0, 1e-9, 0.1, 0.5, 1.0, 2.0, 3.5, 7.0, 10.0, 25.0, 50.0, 100.0, 700.0])
HERMITE_XS = symmetric([0.0, 1e-9, 0.1, 0.5, 0.7071067811865476, 1.0, 1.5, 3.0, 5.0, 10.0, 30.0, 100.0])
LAGUERRE_DEGREES = [0, 1, 2, 3, 4, 5, 6, 8, 10, 20, 30, 50, 100]
HERMITE_DEGREES = [0, 1, 2, 3, 4, 5, 6, 8, 10, 20, 30, 50, 100, 200]
ORDERS = [0.0, 0.5, 1.0, 2.5, 7.0, 0.3, 1e-9, -0.5, -1.0, -2.0, -2.5, -7.0, 1.0 + 2.0 ** -30]

check_recurrences()

with open(OUT, "w") as f:
    f.write(";; Generated by utils/fastmath/dev/generate_laguerre_hermite_reference.py (exact rational arithmetic)\n")
    f.write("{:laguerre [\n")
    total = 0
    for a in ORDERS:
        q = Fraction(a)
        rows = grid_rows([(n, laguerre(n, q)) for n in LAGUERRE_DEGREES], LAGUERRE_XS)
        total += len(rows)
        f.write("  {:order %s :decimal-exact? %s\n   :coefficients [%s]\n   :grid [%s]}\n" % (
            fmt(a), "true" if decimal_exact(a) else "false",
            " ".join(ratios(laguerre(n, q)) for n in range(0, 21)), "\n    ".join(rows)))
    f.write(" ]\n")
    for key, make in ((":hermite-h", hermite_h), (":hermite-he", hermite_he)):
        rows = grid_rows([(n, make(n)) for n in HERMITE_DEGREES], HERMITE_XS)
        total += len(rows)
        f.write(" %s {:coefficients [%s]\n   :grid [%s]}\n" % (
            key, " ".join(ratios(make(n)) for n in range(0, 41)), "\n    ".join(rows)))
    f.write("}\n")

print("grid rows:", total)
print("wrote", OUT)
