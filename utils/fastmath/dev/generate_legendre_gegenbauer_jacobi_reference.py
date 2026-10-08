"""One-off generator of independent reference values for the Legendre, Gegenbauer
(ultraspherical) and Jacobi polynomials of `fastmath.polynomials`: `eval-legendre-P`,
`eval-gegenbauer-C`, `eval-jacobi-P` and their `*-ratio` and polynomial object forms.

Produces `test/resources/polynomials/legendre_gegenbauer_jacobi_reference.edn`. Everything is exact
rational arithmetic (`fractions.Fraction`); the parameters are doubles taken at their exact binary
value. Legendre and Gegenbauer coefficients come from the three term recurrences, the Jacobi ones from
the explicit sum, which is valid for every real parameters, also where the recurrence divides by zero:

  P(n; a, b)(x) = sum_s C(n+a, n-s) C(n+b, s) ((x-1)/2)^s ((x+1)/2)^(n-s)

with the generalized binomial coefficient C(t, k) = t (t-1) ... (t-k+1) / k!. The values and the first
derivatives are evaluated at the exact binary value of every double argument, then rounded once to a
double.

Layout (ratios are `[numerator denominator]` pairs, coefficients ascending in the power of x):
  - `:legendre`    `{:coefficients [[ratio ...] ...] :grid [[n x value derivative] ...]}`
  - `:gegenbauer`  `[{:alpha a :decimal-exact? bool :coefficients ... :grid ...} ...]`
  - `:jacobi`      `[{:alpha a :beta b :decimal-exact? bool :coefficients ... :grid ...} ...]`
`:decimal-exact?` is true when `rationalize` of the parameters gives their exact binary value (for
example 0.25 or 5.5, not 1 + 2^-30): only then the `*-ratio` forms are compared with `:coefficients`.
Rows whose value or derivative is not a finite double are left out.

Run with (from repo root, using the `uv`-managed Python env mentioned in AGENTS.md):
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_legendre_gegenbauer_jacobi_reference.py
"""
from fractions import Fraction

OUT = "/home/ts/clojure/fastmath/test/resources/polynomials/legendre_gegenbauer_jacobi_reference.edn"
DBL_MAX = 1.7976931348623157e308
EPS30 = 2.0 ** -30


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


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def pscale(a, s):
    return [c * s for c in a]


def pmul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def horner(cs, x):
    acc = Fraction(0)
    for c in reversed(cs):
        acc = acc * x + c
    return acc


def derivative(cs):
    return [i * c for i, c in enumerate(cs)][1:] or [Fraction(0)]


def legendre(n):
    a, b = [Fraction(1)], [Fraction(0), Fraction(1)]
    if n == 0:
        return a
    for i in range(2, n + 1):
        a, b = b, pscale(padd(pscale(pmul(b, [Fraction(0), Fraction(1)]), 2 * i - 1), pscale(a, -(i - 1))), Fraction(1, i))
    return b


def gegenbauer(n, al):
    a, b = [Fraction(1)], [Fraction(0), 2 * al]
    if n == 0:
        return a
    for i in range(2, n + 1):
        t = padd(pscale(pmul(b, [Fraction(0), Fraction(1)]), 2 * (al + i - 1)), pscale(a, -(i + 2 * al - 2)))
        a, b = b, pscale(t, Fraction(1, i))
    return b


def gbinom(top, k):
    acc = Fraction(1)
    for j in range(k):
        acc = acc * (top - j) / (j + 1)
    return acc


def jacobi(n, al, be):
    m, p = [Fraction(-1), Fraction(1)], [Fraction(1), Fraction(1)]
    pm, pp = [[Fraction(1)]], [[Fraction(1)]]
    for _ in range(n):
        pm.append(pmul(pm[-1], m))
        pp.append(pmul(pp[-1], p))
    total = [Fraction(0)]
    for s in range(n + 1):
        c = gbinom(n + al, n - s) * gbinom(n + be, s)
        total = padd(total, pscale(pmul(pm[s], pp[n - s]), c))
    return pscale(total, Fraction(1, 2 ** n))


def x_grid():
    xs = {0.0, 1e-9, 0.1, 0.5, 0.7071067811865476, 0.9, 1.5, 3.0, 10.0, 1.0, 1.0 + 1e-6}
    for k in (3, 6, 9, 12):
        xs.add(1.0 - 10.0 ** -k)
    return sorted(xs | {-x for x in xs})


def decimal_exact(*params):
    return all(Fraction(str(p)) == Fraction(p) for p in params)


def grid_rows(cs, degrees_and_polys, xs):
    rows = []
    for n, poly in degrees_and_polys:
        ds = derivative(poly)
        for x in xs:
            xf = Fraction(x)
            v, d = to_double(horner(poly, xf)), to_double(horner(ds, xf))
            if v is None or d is None:
                continue
            rows.append("[%d %s %s %s]" % (n, fmt(x), fmt(v), fmt(d)))
    return rows


XS = x_grid()
LEGENDRE_DEGREES = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 20, 50, 100, 200, 500]
DEGREES = [0, 1, 2, 3, 4, 5, 6, 8, 10, 20, 30]

GEGENBAUER_ALPHAS = [0.0, 0.25, 0.5, 0.5 + EPS30, 0.5 - EPS30, 0.75, 1.0, 1.0 + EPS30, 1.0 - EPS30,
                     1.5, 2.0, 5.5, -0.25, -0.5, -1.0, -1.5, -2.0]
JACOBI_PAIRS = [(0.0, 0.0), (0.5, -0.5), (1.5, 2.5), (-0.5, 0.25), (2.0, 0.0), (-1.0, -1.0), (-1.0, 0.0),
                (0.0, -1.0), (-1.0, -2.0), (-2.0, 0.0), (-0.5, -1.5), (-2.0, -3.0), (-1.0, 3.0), (3.0, -1.0),
                (-3.0, 0.5), (-0.75, -1.25), (-1.0 + 2.0 ** -10, -1.0), (-1.0 + 2.0 ** -20, -1.0),
                (-1.0 + 2.0 ** -30, -1.0), (-1.5 + 2.0 ** -20, -1.5)]

with open(OUT, "w") as f:
    f.write(";; Generated by utils/fastmath/dev/generate_legendre_gegenbauer_jacobi_reference.py (exact rational arithmetic)\n")
    f.write("{:legendre {:coefficients [%s]\n  :grid [%s]}\n" % (
        " ".join(ratios(legendre(n)) for n in range(0, 41)),
        "\n    ".join(grid_rows(None, [(n, legendre(n)) for n in LEGENDRE_DEGREES], XS))))
    f.write(" :gegenbauer [\n")
    total = 0
    for al in GEGENBAUER_ALPHAS:
        a = Fraction(al)
        rows = grid_rows(None, [(n, gegenbauer(n, a)) for n in DEGREES], XS)
        total += len(rows)
        f.write("  {:alpha %s :decimal-exact? %s\n   :coefficients [%s]\n   :grid [%s]}\n" % (
            fmt(al), "true" if decimal_exact(al) else "false",
            " ".join(ratios(gegenbauer(n, a)) for n in range(0, 21)), "\n    ".join(rows)))
    f.write(" ]\n :jacobi [\n")
    for al, be in JACOBI_PAIRS:
        a, b = Fraction(al), Fraction(be)
        rows = grid_rows(None, [(n, jacobi(n, a, b)) for n in DEGREES], XS)
        total += len(rows)
        f.write("  {:alpha %s :beta %s :decimal-exact? %s\n   :coefficients [%s]\n   :grid [%s]}\n" % (
            fmt(al), fmt(be), "true" if decimal_exact(al, be) else "false",
            " ".join(ratios(jacobi(n, a, b)) for n in range(0, 15)), "\n    ".join(rows)))
    f.write(" ]}\n")

print("x points:", len(XS), "gegenbauer+jacobi grid rows:", total)
print("wrote", OUT)
