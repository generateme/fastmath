"""Generator of reference values for the one dimensional interpolators and the grid interpolators of
`fastmath.interpolation.{linear,cubic,monotone,step,sprague,barycentric}` and for isotonic regression
(`fastmath.ml.regression/pava`).

Oracles (independent implementations, never the fastmath code):

    :cubic        scipy `CubicSpline(bc_type='natural')`, extrapolate=True (continues the end cubic)
    :linear       `numpy.interp` (inside) and scipy `interp1d(fill_value='extrapolate')` (all queries)
    :pchip        scipy `PchipInterpolator(extrapolate=True)`
    :monotone     an independent Fritsch-Carlson clone written below (`fc_*`); the clone has switches
                  that remove the last-interval limiting and the sign test, used for the
                  `:fc-needed` counters (the two repairs are demonstrably needed)
    :polynomial   mpmath Lagrange polynomial, 60 digits, from the exact double inputs; scipy
                  `BarycentricInterpolator` as a second (double precision) oracle
    :barycentric  mpmath Floater-Hormann rational interpolant (paper formula for the weights), 50 digits
    :sprague      colour-science `SpragueInterpolator`; outside the domain the same formula with the
                  index of the first / last interval forced (colour raises ValueError there)
    :step-*       scipy `interp1d` kinds 'next' (step-before), 'previous' (step-after),
                  'nearest-up' (step with point 0.5, ties go to the right knot)
    :isotonic     scikit-learn `IsotonicRegression(increasing=True / False)` with sample weights
    :bilinear     scipy `RegularGridInterpolator(method='linear')`, extrapolating
    :cubic-2d     two passes of scipy `CubicSpline(bc_type='natural')` (along ys, then along xs)

Domain: small hand-made datasets (uniform, nonuniform, negative, tiny and huge scale, two and three
points, plateaus, sign changes, one hundred Chebyshev nodes) and seeded random datasets (numpy
`default_rng`, fixed seeds). Queries: all knots, midpoints, seeded random points inside the domain
and points outside on both sides. Why: the previous implementations had no oracle tests; `monotone`
lost monotonicity (the last interval was not limited and there was no sign test), `barycentric` found
exact node hits with a binary search, `:desc` of isotonic regression meant ascending.

Produces `test/resources/interpolation/interp1d_reference.edn`: a map of sections, every case a map
with `:xs`, `:ys`, `:qs` (queries) and the oracle values (`:expected`, ...). The output is
deterministic; run it twice and diff.

Run from anywhere:
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_interp1d_reference.py
"""
from pathlib import Path

import colour
import mpmath
import numpy as np
import scipy
import sklearn
from colour.algebra import SpragueInterpolator
from scipy.interpolate import (BarycentricInterpolator, CubicSpline, PchipInterpolator,
                               RegularGridInterpolator, interp1d)
from sklearn.isotonic import IsotonicRegression

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "test" / "resources" / "interpolation" / "interp1d_reference.edn"

mpmath.mp.dps = 60


# ---------------------------------------------------------------------------------------------
# EDN writing

def edn(o):
    if isinstance(o, dict):
        return "{" + " ".join(":%s %s" % (k, edn(v)) for k, v in o.items()) + "}"
    if isinstance(o, (list, tuple, np.ndarray)):
        return "[" + " ".join(edn(v) for v in o) + "]"
    if isinstance(o, (bool, np.bool_)):
        return "true" if o else "false"
    if isinstance(o, (int, np.integer)):
        return str(int(o))
    if isinstance(o, (float, np.floating)):
        s = repr(float(o))
        assert "inf" not in s and "nan" not in s, s
        return s
    if isinstance(o, str):
        return '"%s"' % o
    raise TypeError(type(o))


def section(name, cases):
    lines = [" :%s [" % name]
    for c in cases:
        lines.append("  " + edn(c))
    lines.append(" ]")
    return "\n".join(lines)


# ---------------------------------------------------------------------------------------------
# queries

def queries(xs, rng, n_random=12, outside=(0.01, 0.5, 2.0)):
    xs = np.asarray(xs, dtype=float)
    lo, hi = xs[0], xs[-1]
    w = hi - lo
    q = list(xs)
    q += list(0.5 * (xs[1:] + xs[:-1]))
    q += list(rng.uniform(lo, hi, n_random))
    for f in outside:
        q += [lo - f * w, hi + f * w]
    return np.array(q)


def case(name, xs, ys, rng, **kw):
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    return {"name": name, "xs": xs, "ys": ys, "qs": queries(xs, rng, **kw)}


def raw(s):
    return s


# ---------------------------------------------------------------------------------------------
# data sets

def datasets(rng):
    ds = []
    ds.append(case("two", [0.0, 2.0], [1.0, 5.0], rng))
    ds.append(case("three", [0.0, 1.0, 3.0], [1.0, -2.0, 4.0], rng))
    x = np.arange(6.0)
    ds.append(case("uni6", x, np.sin(x), rng))
    x = np.array([0.0, 0.3, 0.35, 1.1, 1.6, 2.9, 3.0, 4.4])
    ds.append(case("nonuni8", x, np.cos(2.0 * x) + 0.3 * x, rng))
    x = np.array([-3.0, -2.5, -1.0, 0.5, 4.0, 4.25])
    ds.append(case("negative", x, np.array([2.0, -1.0, 0.5, 3.0, -2.0, 1.0]), rng))
    x = 1.0e-6 * np.arange(7.0)
    ds.append(case("tiny-x", x, 1.0e8 * np.sin(1.0e6 * x), rng))
    x = 1.0e5 * np.arange(1.0, 8.0)
    ds.append(case("huge-x", x, 1.0e-3 * np.cos(x / 3.0e5), rng))
    x = np.sort(rng.uniform(-5.0, 5.0, 60))
    ds.append(case("many60", x, np.sin(1.3 * x) + 0.1 * rng.standard_normal(60), rng))
    x = np.linspace(0.0, 10.0, 200)
    ds.append(case("many200", x, np.exp(-0.2 * x) * np.cos(x), rng))
    return ds


def shape_datasets(rng):
    """Datasets that stress shape preservation (plateaus, sign changes, steep jumps)."""
    ds = []
    ds.append(case("plateau", np.arange(9.0), np.array([0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 5.0, 5.0, 5.0]), rng))
    ds.append(case("sign-changes", np.arange(7.0), np.array([0.0, 1.0, 0.0, -1.0, 0.0, 2.0, -3.0]), rng))
    ds.append(case("last-interval", np.arange(4.0), np.array([0.0, 1.0, 2.0, 2.01]), rng))
    ds.append(case("extremum", np.arange(4.0), np.array([0.0, 3.0, 2.0, 0.0]), rng))
    ds.append(case("steep", np.array([0.0, 1.0, 1.1, 2.0, 5.0, 5.2]), np.array([0.0, 0.01, 10.0, 10.5, 10.6, 30.0]), rng))
    ds.append(case("decreasing", np.arange(6.0), np.array([9.0, 7.0, 7.0, 3.0, 2.9, -4.0]), rng))
    ds.append(case("two", np.array([0.0, 2.0]), np.array([1.0, 5.0]), rng))
    ds.append(case("three-flat", np.arange(3.0), np.array([1.0, 1.0, 1.0]), rng))
    x = np.sort(rng.uniform(0.0, 10.0, 25))
    ds.append(case("random-sorted", x, rng.standard_normal(25), rng))
    ds.append(case("tiny-y", np.arange(5.0), 1.0e-200 * np.array([0.0, 1.0, 1.0, 3.0, 3.5]), rng))
    ds.append(case("huge-y", np.arange(5.0), 1.0e150 * np.array([0.0, 1.0, 1.0, 3.0, 3.5]), rng))
    return ds


# ---------------------------------------------------------------------------------------------
# cubic, linear, pchip (scipy)

def make_cubic(ds):
    out = []
    for d in ds:
        f = CubicSpline(d["xs"], d["ys"], bc_type="natural", extrapolate=True)
        out.append(dict(d, expected=f(d["qs"])))
    return out


def make_linear(ds):
    out = []
    for d in ds:
        f = interp1d(d["xs"], d["ys"], kind="linear", fill_value="extrapolate", assume_sorted=True)
        e = f(d["qs"])
        inside = (d["qs"] >= d["xs"][0]) & (d["qs"] <= d["xs"][-1])
        assert np.allclose(e[inside], np.interp(d["qs"][inside], d["xs"], d["ys"]), rtol=1e-12, atol=0)
        out.append(dict(d, expected=e))
    return out


def make_pchip(ds):
    out = []
    for d in ds:
        f = PchipInterpolator(d["xs"], d["ys"], extrapolate=True)
        out.append(dict(d, expected=f(d["qs"])))
    return out


# ---------------------------------------------------------------------------------------------
# independent Fritsch-Carlson clone

def fc_tangents(xs, ys, last_interval=True, sign_test=True):
    """Fritsch and Carlson (1980): arithmetic mean interior tangents, one sided end tangents, tangents set
    to zero where the secants change sign or vanish, then every interval is limited to the circle
    alpha^2 + beta^2 <= 9."""
    n = len(xs)
    d = [(ys[k + 1] - ys[k]) / (xs[k + 1] - xs[k]) for k in range(n - 1)]
    m = [0.0] * n
    m[0] = d[0]
    m[n - 1] = d[n - 2]
    for k in range(1, n - 1):
        a, b = d[k - 1], d[k]
        if sign_test and (a == 0.0 or b == 0.0 or (a > 0.0) != (b > 0.0)):
            m[k] = 0.0
        else:
            m[k] = 0.5 * (a + b)
    stop = n - 1 if last_interval else n - 2
    for k in range(stop):
        if d[k] == 0.0:
            m[k] = 0.0
            m[k + 1] = 0.0
            continue
        alpha = m[k] / d[k]
        beta = m[k + 1] / d[k]
        s = alpha * alpha + beta * beta
        if s > 9.0:
            tau = 3.0 / s ** 0.5
            m[k] = tau * alpha * d[k]
            m[k + 1] = tau * beta * d[k]
    return m


def hermite(xs, ys, m, x):
    n = len(xs)
    i = int(np.searchsorted(xs, x, side="right")) - 1
    i = min(max(i, 0), n - 2)
    h = xs[i + 1] - xs[i]
    t = (x - xs[i]) / h
    h00 = (1 + 2 * t) * (1 - t) ** 2
    h10 = t * (1 - t) ** 2
    h01 = t * t * (3 - 2 * t)
    h11 = t * t * (t - 1)
    return h00 * ys[i] + h10 * h * m[i] + h01 * ys[i + 1] + h11 * h * m[i + 1]


def fc_eval(xs, ys, qs, **kw):
    xs = [float(v) for v in xs]
    ys = [float(v) for v in ys]
    m = fc_tangents(xs, ys, **kw)
    out = []
    for x in qs:
        # exact knot: the value of the data (the formula gives it up to rounding)
        hit = [k for k in range(len(xs)) if xs[k] == x]
        out.append(ys[hit[0]] if hit else hermite(xs, ys, m, float(x)))
    return np.array(out)


def make_monotone(ds):
    out = []
    for d in ds:
        out.append(dict(d, expected=fc_eval(d["xs"], d["ys"], d["qs"])))
    return out


def violations(xs, ys, last_interval, sign_test, grid_points=10000):
    """Number of intervals where the clone leaves [min(y_i, y_i+1), max(y_i, y_i+1)] by more than 1e-12."""
    xs = [float(v) for v in xs]
    ys = [float(v) for v in ys]
    m = fc_tangents(xs, ys, last_interval=last_interval, sign_test=sign_test)
    bad = 0
    for i in range(len(xs) - 1):
        lo, hi = min(ys[i], ys[i + 1]), max(ys[i], ys[i + 1])
        tol = 1e-12 * max(1.0, hi - lo)
        for t in np.linspace(0.0, 1.0, grid_points // (len(xs) - 1) + 1):
            x = xs[i] + t * (xs[i + 1] - xs[i])
            v = hermite(xs, ys, m, x)
            if v < lo - tol or v > hi + tol:
                bad += 1
                break
    return bad


def random_monotone_sets(rng, count=200):
    sets = []
    for k in range(count):
        n = int(rng.integers(5, 16))
        steps = rng.uniform(0.05, 2.0, n - 1)
        steps[rng.random(n - 1) < 0.15] = 0.01
        xs = np.concatenate([[0.0], np.cumsum(steps)])
        inc = np.zeros(n - 1)
        u = rng.random(n - 1)
        for j in range(n - 1):
            if u[j] < 0.3:
                inc[j] = 0.0
            elif u[j] < 0.5:
                inc[j] = rng.uniform(5.0, 50.0)
            else:
                inc[j] = rng.uniform(0.001, 1.0)
        ys = np.concatenate([[0.0], np.cumsum(inc)])
        if k % 2 == 1:
            ys = -ys
        sets.append({"xs": xs, "ys": ys})
    return sets


def random_general_sets(rng, count=100):
    sets = []
    for k in range(count):
        n = int(rng.integers(3, 13))
        xs = np.concatenate([[rng.uniform(-3.0, 3.0)], np.zeros(n - 1)])
        xs[1:] = xs[0] + np.cumsum(rng.uniform(0.05, 2.0, n - 1))
        ys = rng.standard_normal(n) * rng.choice([1.0, 10.0, 0.01])
        flat = rng.random(n) < 0.2
        for j in range(1, n):
            if flat[j]:
                ys[j] = ys[j - 1]
        sets.append({"xs": xs, "ys": ys})
    return sets


# ---------------------------------------------------------------------------------------------
# polynomial (mpmath Lagrange) and Floater-Hormann (mpmath)

def mp_lagrange(xs, ys, x):
    xs = [mpmath.mpf(float(v)) for v in xs]
    ys = [mpmath.mpf(float(v)) for v in ys]
    x = mpmath.mpf(float(x))
    total = mpmath.mpf(0)
    for k in range(len(xs)):
        term = ys[k]
        for j in range(len(xs)):
            if j != k:
                term *= (x - xs[j]) / (xs[k] - xs[j])
        total += term
    return float(total)


def mp_floater_hormann(xs, ys, d, x):
    """Floater and Hormann (2007), eq. (18): w_k = (-1)^(k-d) sum_{i in J_k} prod_{j=i, j!=k}^{i+d} 1/|x_k-x_j|,
    J_k = {i : k-d <= i <= k, 0 <= i <= n-d} (n+1 nodes in the paper, here len(xs) - 1 = n)."""
    xm = [mpmath.mpf(float(v)) for v in xs]
    ym = [mpmath.mpf(float(v)) for v in ys]
    xq = mpmath.mpf(float(x))
    n = len(xm) - 1
    for k in range(n + 1):
        if xm[k] == xq:
            return float(ym[k])
    num = mpmath.mpf(0)
    den = mpmath.mpf(0)
    for k in range(n + 1):
        s = mpmath.mpf(0)
        for i in range(max(0, k - d), min(k, n - d) + 1):
            p = mpmath.mpf(1)
            for j in range(i, i + d + 1):
                if j != k:
                    p /= abs(xm[k] - xm[j])
            s += p
        w = (-1) ** (k - d) * s
        num += w * ym[k] / (xq - xm[k])
        den += w / (xq - xm[k])
    return float(num / den)


def chebyshev(n, a, b):
    k = np.arange(n)
    c = np.cos(np.pi * (2 * k + 1) / (2 * n))[::-1]
    return 0.5 * (a + b) + 0.5 * (b - a) * c


def make_polynomial(rng):
    out = []
    specs = [
        ("two", np.array([0.0, 2.0]), lambda x: 3.0 * x - 1.0),
        ("three", np.array([0.0, 1.0, 3.0]), lambda x: x * x - 2.0),
        ("equi5", np.linspace(-1.0, 1.0, 5), lambda x: 1.0 / (1.0 + 25.0 * x * x)),
        ("equi9", np.linspace(-1.0, 1.0, 9), lambda x: np.exp(x) * np.sin(3.0 * x)),
        ("cheb15", chebyshev(15, -1.0, 1.0), lambda x: 1.0 / (1.0 + 25.0 * x * x)),
        ("equi15", np.linspace(0.0, 1.0, 15), lambda x: np.cos(5.0 * x)),
        ("nonuni6", np.array([-2.0, -1.5, 0.1, 0.9, 2.5, 7.0]), lambda x: np.sin(x) + 0.1 * x),
        ("shift-scale", 1.0e4 + 1.0e-3 * np.arange(8.0), lambda x: np.sin(1.0e3 * (x - 1.0e4))),
    ]
    for name, xs, f in specs:
        ys = f(xs)
        qs = queries(xs, rng, n_random=8, outside=(0.01, 0.3))
        out.append({"name": name, "xs": xs, "ys": ys, "qs": qs,
                    "expected": np.array([mp_lagrange(xs, ys, q) for q in qs]),
                    "scipy": BarycentricInterpolator(xs, ys, random_state=0)(qs)})
    # many nodes on a tiny interval: the node polynomial 1 / prod(x_k - x_j) underflows without scaling
    xs = chebyshev(100, 0.0, 1.0e-4)
    ys = np.exp(1.0e4 * xs)
    # inside the domain only: the second barycentric form (like Neville) loses accuracy when it extrapolates a degree 99 polynomial
    qs = np.concatenate([xs[::17], rng.uniform(0.0, 1.0e-4, 12)])
    out.append({"name": "cheb100-tiny", "xs": xs, "ys": ys, "qs": qs,
                "expected": np.array([mp_lagrange(xs, ys, q) for q in qs]),
                "scipy": BarycentricInterpolator(xs, ys, random_state=0)(qs)})
    return out


def make_barycentric(rng):
    out = []
    specs = [
        ("nonuni8", np.array([0.0, 0.3, 0.35, 1.1, 1.6, 2.9, 3.0, 4.4]),
         lambda x: np.cos(2.0 * x) + 0.3 * x, (0, 1, 2, 3, 5)),
        ("uni10", np.linspace(-2.0, 2.0, 10), lambda x: np.abs(x) + 0.2 * x, (0, 1, 2, 3, 4)),
        ("negative6", np.array([-3.0, -2.5, -1.0, 0.5, 4.0, 4.25]),
         lambda x: x * np.sin(x), (0, 1, 2, 4)),
        ("two", np.array([0.0, 2.0]), lambda x: 3.0 * x - 1.0, (0, 1)),
    ]
    for name, xs, f, orders in specs:
        ys = f(xs)
        qs = queries(xs, rng, n_random=8, outside=(0.01, 0.3))
        for d in orders:
            out.append({"name": name, "order": d, "xs": xs, "ys": ys, "qs": qs,
                        "expected": np.array([mp_floater_hormann(xs, ys, d, q) for q in qs])})
    return out


# ---------------------------------------------------------------------------------------------
# sprague (colour-science)

class SpragueForced(SpragueInterpolator):
    """colour's formula with a forced interval index: the same evaluation outside the domain."""

    def _validate_interpolation_range(self, x):
        pass

    def at(self, x, i):
        X = (x - self._xp[i]) / (self._xp[i + 1] - self._xp[i])
        r = self._yp
        r_s = np.asarray((r[i - 2], r[i - 1], r[i], r[i + 1], r[i + 2], r[i + 3]))
        w_s = np.asarray(((2, -16, 0, 16, -2, 0), (-1, 16, -30, 16, -1, 0), (-9, 39, -70, 66, -33, 7),
                          (13, -64, 126, -124, 61, -12), (-5, 25, -50, 50, -25, 5)))
        a = np.dot(w_s, r_s) / 24
        return float(r[i] + sum(a[k] * X ** (k + 1) for k in range(5)))


def make_sprague(rng):
    out = []
    specs = [
        ("min6", -1.0 + 0.5 * np.arange(6.0), lambda x: np.sin(x)),
        ("n7", np.arange(7.0), lambda x: x ** 3 - 2.0 * x),
        ("n8-offset", 3.25 + 0.25 * np.arange(8.0), lambda x: np.exp(-x) * np.cos(4.0 * x)),
        ("n25", -4.0 + 0.5 * np.arange(25.0), lambda x: np.sin(x) + 0.1 * x),
        ("n12-scaled", 1.0e6 * (1.0 + np.arange(12.0)), lambda x: 1.0e-3 * np.sin(x / 3.0e6)),
    ]
    for name, xs, f in specs:
        ys = f(xs)
        h0 = xs[1] - xs[0]
        inside = queries(xs, rng, n_random=12, outside=())
        interp = SpragueForced(xs, ys)
        pub = SpragueInterpolator(xs, ys)(inside)
        # sanity: the forced-index formula equals the colour public call inside the domain
        e_in = np.array([interp.at(q, max(int(np.searchsorted(interp._xp, q)) - 1, 2)) for q in inside])
        assert np.allclose(pub, e_in, rtol=1e-12, atol=1e-12 * max(1.0, np.abs(ys).max()))
        n = len(xs)
        left = [xs[0] - 0.3 * h0, xs[0] - 1.0 * h0, xs[0] - 5.0 * h0]
        right = [xs[-1] + 0.3 * h0, xs[-1] + 1.0 * h0, xs[-1] + 5.0 * h0]
        e_l = [interp.at(q, 2) for q in left]
        e_r = [interp.at(q, n) for q in right]
        out.append({"name": name, "xs": xs, "ys": ys, "qs": inside, "expected": pub,
                    "qs-outside": np.array(left + right), "expected-outside": np.array(e_l + e_r)})
    return out


# ---------------------------------------------------------------------------------------------
# step functions (scipy interp1d)

def make_step(rng):
    out = []
    specs = [
        ("irregular5", np.array([0.0, 1.0, 3.0, 4.0, 6.0]), np.array([2.0, -1.0, 5.0, 0.5, 7.0])),
        ("two", np.array([-1.0, 1.0]), np.array([4.0, -4.0])),
        ("uni8", np.arange(8.0), np.array([3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0])),
        ("negative", np.array([-8.0, -4.0, -2.0, 0.0, 1.0]), np.array([0.5, 0.25, -3.0, 1.5, 9.0])),
    ]
    for name, xs, ys in specs:
        # midpoints of dyadic knots are exact: ties are tested exactly
        qs = list(xs) + list(0.5 * (xs[1:] + xs[:-1])) + list(rng.uniform(xs[0], xs[-1], 12)) \
            + [xs[0] - 3.0, xs[0] - 1e-9, xs[-1] + 1e-9, xs[-1] + 3.0]
        qs = np.array(qs)
        fill = (ys[0], ys[-1])
        mk = lambda kind: interp1d(xs, ys, kind=kind, bounds_error=False, fill_value=fill, assume_sorted=True)(qs)
        out.append({"name": name, "xs": xs, "ys": ys, "qs": qs,
                    "step-before": mk("next"), "step-after": mk("previous"), "step": mk("nearest-up")})
    return out


# ---------------------------------------------------------------------------------------------
# isotonic regression (scikit-learn)

def make_isotonic(rng):
    out = []
    for k in range(14):
        n = [2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 60][k]
        ys = np.round(rng.standard_normal(n) * 3.0, 1) if k % 2 == 0 else rng.standard_normal(n)
        if k % 3 == 0:
            ys = np.round(ys)  # many ties
        ws = rng.uniform(0.1, 5.0, n)
        if k % 4 == 1:
            ws = np.round(ws)
            ws[ws == 0] = 1.0
        if k == 13:
            ws = 10.0 ** rng.uniform(-3, 3, n)
        x = np.arange(n, dtype=float)
        inc = IsotonicRegression(increasing=True).fit_transform(x, ys, sample_weight=ws)
        dec = IsotonicRegression(increasing=False).fit_transform(x, ys, sample_weight=ws)
        out.append({"ys": ys, "ws": ws, "increasing": inc, "decreasing": dec})
    return out


# ---------------------------------------------------------------------------------------------
# grids

def make_grids(rng):
    bil, cub = [], []
    specs = [
        ("2x2", np.array([0.0, 1.0]), np.array([0.0, 2.0]), lambda x, y: 1.0 + 2.0 * x - y + 0.5 * x * y),
        ("3x2", np.array([0.0, 1.0, 3.0]), np.array([-1.0, 1.0]), lambda x, y: np.sin(x) + y * y),
        ("uni5x4", np.arange(5.0), np.arange(4.0), lambda x, y: np.sin(x) * np.cos(y)),
        ("nonuni7x6", np.array([0.0, 0.4, 0.5, 1.5, 2.0, 3.5, 4.0]), np.array([-1.0, -0.2, 0.3, 1.1, 2.0, 2.2]),
         lambda x, y: np.exp(-0.3 * x) * np.sin(y) + 0.1 * x * y),
        ("uni9x9", np.linspace(0.0, 1.0, 9), np.linspace(0.0, 1.0, 9), lambda x, y: np.exp(x) * np.sin(3.0 * y)),
    ]
    for name, xs, ys, f in specs:
        X, Y = np.meshgrid(xs, ys, indexing="ij")
        vss = f(X, Y)
        qx = queries(xs, rng, n_random=6, outside=(0.05, 0.4))
        qy = queries(ys, rng, n_random=6, outside=(0.05, 0.4))
        m = min(len(qx), len(qy))
        pts = np.stack([rng.permutation(qx)[:m], rng.permutation(qy)[:m]], axis=1)
        # exact nodes and cell corners / outside corners
        nodes = np.array([[xs[0], ys[0]], [xs[-1], ys[-1]], [xs[0], ys[-1]], [xs[-1] + 0.3, ys[0] - 0.3],
                          [xs[0] - 0.2, ys[-1] + 0.2], [xs[len(xs) // 2], ys[len(ys) // 2]]])
        pts = np.vstack([pts, nodes])
        rgi = RegularGridInterpolator((xs, ys), vss, method="linear", bounds_error=False, fill_value=None)
        bil.append({"name": name, "xs": xs, "ys": ys, "vss": [list(r) for r in vss],
                    "pts": [list(p) for p in pts], "expected": rgi(pts)})

        def natural_2d(x, y):
            rows = np.array([CubicSpline(ys, vss[i], bc_type="natural")(y) for i in range(len(xs))])
            return float(CubicSpline(xs, rows, bc_type="natural")(x))

        cub.append({"name": name, "xs": xs, "ys": ys, "vss": [list(r) for r in vss],
                    "pts": [list(p) for p in pts], "expected": np.array([natural_2d(p[0], p[1]) for p in pts])})
    return bil, cub


# ---------------------------------------------------------------------------------------------

def main():
    rng = np.random.default_rng(20261014)
    sections = []

    base = datasets(rng)
    sections.append(section("cubic", make_cubic(base)))
    sections.append(section("linear", make_linear(base)))
    shape = shape_datasets(rng)
    sections.append(section("pchip", make_pchip(base + shape)))
    sections.append(section("monotone", make_monotone(shape + base[:6])))

    mono = random_monotone_sets(rng)
    general = random_general_sets(rng)
    sections.append(section("random-monotone", mono))
    sections.append(section("random-general", general))

    # the two repairs are demonstrably needed: count the intervals that leave their data range
    sets = mono + general + [{"xs": np.arange(4.0), "ys": np.array([0.0, 1.0, 2.0, 2.01])},
                             {"xs": np.arange(4.0), "ys": np.array([0.0, 3.0, 2.0, 0.0])}]
    needed = {}
    for label, kw in (("full", dict(last_interval=True, sign_test=True)),
                      ("no-last-interval", dict(last_interval=False, sign_test=True)),
                      ("no-sign-test", dict(last_interval=True, sign_test=False)),
                      ("neither", dict(last_interval=False, sign_test=False))):
        needed[label] = sum(1 for s in sets if violations(s["xs"], s["ys"], **kw) > 0)
    for label in ("no-last-interval", "no-sign-test", "neither"):
        assert needed[label] > 0, needed
    assert needed["full"] == 0, needed
    sections.append(" :fc-needed " + edn(needed))
    print("datasets with an interval outside its data range:", needed)

    sections.append(section("polynomial", make_polynomial(rng)))
    sections.append(section("barycentric", make_barycentric(rng)))
    sections.append(section("sprague", make_sprague(rng)))
    sections.append(section("step", make_step(rng)))
    sections.append(section("isotonic", make_isotonic(rng)))
    bil, cub = make_grids(rng)
    sections.append(section("bilinear", bil))
    sections.append(section("cubic-2d", cub))

    header = (";; Reference values for the 1D and grid interpolators and isotonic regression.\n"
              ";; Generated by utils/fastmath/dev/generate_interp1d_reference.py; do not edit.\n"
              ";; scipy %s, numpy %s, scikit-learn %s, mpmath %s, colour-science %s\n"
              % (scipy.__version__, np.__version__, sklearn.__version__, mpmath.__version__, colour.__version__))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(header + "{\n" + "\n".join(sections) + "\n}\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
