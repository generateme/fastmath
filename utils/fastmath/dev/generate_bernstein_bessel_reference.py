"""One-off generator of independent reference values for the Bernstein basis polynomials and the Bessel
polynomials of `fastmath.polynomials`: `eval-bernstein`, `bernstein`, `eval-bessel-y`, `eval-bessel-t` and
their `*-ratio` and polynomial object forms.

Produces `test/resources/polynomials/bernstein_bessel_reference.edn`. Everything is exact rational
arithmetic (`fractions.Fraction`); the arguments are doubles taken at their exact binary value. The
polynomials come from the explicit sums (not from the recurrences the library uses):

  b(k, n)(x) = C(n, k) x^k (1-x)^(n-k)                                  (0 outside 0 <= k <= n)
  y_n(x)     = sum_k (n+k)! / ((n-k)! k!) (x/2)^k                       (Bessel polynomial)
  theta_n(x) = x^n y_n(1/x) = sum_k (n+k)! / ((n-k)! k! 2^k) x^(n-k)    (reverse Bessel polynomial)

and the generator asserts the three term recurrences of y and theta. The values and the first derivatives
are evaluated at the exact binary value of every double argument, then rounded once to a double. A row is
left out when the value or the derivative is not a finite double, or when the value is non-zero and below
1e-290 in magnitude (subnormal results are not a target of the tests).

Layout (ratios are `[numerator denominator]` pairs, coefficients ascending in the power of x):
  - `:bernstein`   `{:grid [[degree order x value] ...] :coefficients [[degree order [ratio ...]] ...]}`
  - `:bessel-y`    `{:coefficients [[ratio ...] ...] :grid [[n x value derivative] ...]}`
  - `:bessel-t`    `{:coefficients [[ratio ...] ...] :grid [[n x value derivative] ...]}`

Run with (from repo root, using the `uv`-managed Python env mentioned in AGENTS.md):
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_bernstein_bessel_reference.py
"""
from fractions import Fraction
from math import comb, factorial

OUT = "/home/ts/clojure/fastmath/test/resources/polynomials/bernstein_bessel_reference.edn"
DBL_MAX = 1.7976931348623157e308


def to_double(q):
    try:
        v = float(Fraction(q))
    except OverflowError:
        return None
    return v if abs(v) <= DBL_MAX else None


def keep(v):
    return v is not None and (v == 0.0 or abs(v) >= 1e-290)


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


def bessel_y(n):
    return [Fraction(factorial(n + k), factorial(n - k) * factorial(k) * 2 ** k) for k in range(n + 1)]


def bessel_t(n):
    return [Fraction(factorial(n + (n - j)), factorial(j) * factorial(n - j) * 2 ** (n - j)) for j in range(n + 1)]


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def check_recurrences():
    """y_n = (2n-1) x y_(n-1) + y_(n-2);  theta_n = (2n-1) theta_(n-1) + x^2 theta_(n-2)."""
    for n in range(2, 30):
        lhs = bessel_y(n)
        rhs = padd([Fraction(0)] + [(2 * n - 1) * c for c in bessel_y(n - 1)], bessel_y(n - 2))
        assert lhs == rhs, ("y", n)
        lhs = bessel_t(n)
        rhs = padd([(2 * n - 1) * c for c in bessel_t(n - 1)], [Fraction(0), Fraction(0)] + bessel_t(n - 2))
        assert lhs == rhs, ("theta", n)


def bernstein_value(n, k, x):
    if k < 0 or k > n:
        return Fraction(0)
    return comb(n, k) * x ** k * (1 - x) ** (n - k)


def bernstein_coefficients(n, k):
    """Coefficients of b(k, n) in powers of x: (-1)^(l-k) C(n, l) C(l, k)."""
    if k < 0 or k > n:
        return [Fraction(0)] * (n + 1)
    return [Fraction(0) if l < k else Fraction((-1) ** (l - k) * comb(n, l) * comb(l, k)) for l in range(n + 1)]


def symmetric(xs):
    return sorted(set(xs) | {-x for x in xs})


BERNSTEIN_X = [0.0, 1e-9, 1e-4, 0.01, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.99, 1 - 1e-9, 1.0, 1.5, 2.0, -0.5, -2.0, 1e-3, 0.999]
BERNSTEIN_DEGREES = [0, 1, 2, 3, 4, 5, 10, 30, 100, 200, 500, 1000, 1100, 2000, 5000]
BESSEL_X = symmetric([0.0, 1e-9, 0.1, 0.5, 1.0, 2.0, 3.5, 7.0, 10.0, 30.0, 100.0])
BESSEL_DEGREES = [0, 1, 2, 3, 4, 5, 6, 8, 10, 20, 30, 50, 100]


def orders(n):
    return sorted({o for o in (-1, 0, 1, n // 4, n // 2, n - 1, n, n + 1, n + 3) if True})


check_recurrences()

bern_rows = []
for n in BERNSTEIN_DEGREES:
    for k in orders(n):
        for x in BERNSTEIN_X:
            v = to_double(bernstein_value(n, k, Fraction(x)))
            if keep(v):
                bern_rows.append("[%d %d %s %s]" % (n, k, fmt(x), fmt(v)))

bern_coefs = []
for n in (0, 1, 2, 3, 4, 5, 8, 12, 20, 40):
    for k in orders(n):
        bern_coefs.append("[%d %d %s]" % (n, k, ratios(bernstein_coefficients(n, k))))


def grid(make):
    rows = []
    for n in BESSEL_DEGREES:
        poly = make(n)
        ds = derivative(poly)
        for x in BESSEL_X:
            xf = Fraction(x)
            v, d = to_double(horner(poly, xf)), to_double(horner(ds, xf))
            if keep(v) and d is not None:
                rows.append("[%d %s %s %s]" % (n, fmt(x), fmt(v), fmt(d)))
    return rows


with open(OUT, "w") as f:
    f.write(";; Generated by utils/fastmath/dev/generate_bernstein_bessel_reference.py (exact rational arithmetic)\n")
    f.write("{:bernstein {:grid [%s]\n  :coefficients [%s]}\n" % ("\n    ".join(bern_rows), "\n    ".join(bern_coefs)))
    for key, make in ((":bessel-y", bessel_y), (":bessel-t", bessel_t)):
        f.write(" %s {:coefficients [%s]\n   :grid [%s]}\n" % (
            key, " ".join(ratios(make(n)) for n in range(0, 41)), "\n    ".join(grid(make))))
    f.write("}\n")

print("bernstein rows:", len(bern_rows), "coefficient sets:", len(bern_coefs))
print("wrote", OUT)
