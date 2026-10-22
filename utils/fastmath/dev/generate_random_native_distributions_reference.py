"""Generator of reference data for the native distributions of `fastmath.random.distributions`
(Random Audit WP 5.2, 5.5, 7.2, 7.5, 7.6, 7.7).

Tools: mpmath 1.3.0 (oracle for `logarithmic` and `inverse-gaussian`), scipy 1.18.1 (`stats.dirichlet`,
`stats.multinomial`, `stats.multivariate_normal`, `stats.invgauss` as a cross-check of the mpmath values), numpy 2.5.2
(`default_rng(seed)` for the generated parameters and points only).  No reference value depends on an unseeded draw;
every random choice is made with a fixed seed.

Precision: mpmath with 60 digits for `logarithmic` and `inverse-gaussian`, rounded once to double; scipy values are
double precision (the closed forms: dirichlet, multinomial, multivariate normal; relative error ~1e-15).  The `x`
and `p` arguments are exact doubles (printed with `repr`), the reference is a function of that exact double.

Domain:
    :logarithmic        p in .1 .5 .9 .99 .999; cdf at every k in 1..200 and at 78 more k up to 1e5 (k >= 150
                        included); the quantile (smallest k with cdf(k) >= p) for 750 seeded uniform p, the tail
                        levels 1e-6, 1e-3, 1-1e-4, 1-1e-6, 1-1e-9, 1-1e-12.
    :dirichlet          alpha [2 3 4], [.5 .5 .5], [1 1], [5], [.3 1.2 2.5 4 .9]; 5 seeded interior points each,
                        pdf and log-pdf (scipy).
    :multinomial        (n, p): (10, [.2 .3 .5]), (1, [.25 .25 .5]), (20, [.1 .2 .3 .4]), (6, [.5 0 .5]) pmf at the
                        points (zero-probability category included; scipy), cdf (exact rational enumeration) at
                        a few points.
    :multi-normal       dimensions 1..5 (seeded SPD covariance, means), 4 interior points each, pdf and log-pdf
                        (scipy); the far tail: Sigma [[1 .5] [.5 2]] at (40, -40), a 3-d case at (30, 30, 30), a
                        1-d case at 60 (log-pdf only, the pdf underflows).
    :inverse-gaussian   (mu, lambda): lambda/mu in 1 100 354 355 400 1000 1e4 (mu = 1), plus mu = .01 lambda = 100,
                        mu = 5 lambda = 1775 (lambda/mu = 355) and mu = 2 lambda = 200: cdf (mpmath, 60 digits) at
                        x = mu (1 + k / sqrt(lambda/mu)) for k in -4..8 and at x = .2 mu, .5 mu.  Cross-checked
                        against scipy `invgauss(mu/lambda, scale=lambda)` where the scipy value is finite
                        (relative 1e-8, asserted).

Run (writes test/resources/random/native_distributions_reference.edn, byte-identical on rerun):
    cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_random_native_distributions_reference.py
"""
import itertools
import math
import os
import random
import warnings
from fractions import Fraction

import mpmath as mp
import numpy as np
import scipy
import scipy.stats as st

mp.mp.dps = 60
warnings.filterwarnings("ignore")

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                   "test", "resources", "random", "native_distributions_reference.edn")


# ---------------------------------------------------------------- EDN writer

def edn(v, indent=0):
    if isinstance(v, bool):
        return "true" if v else "false"
    if v is None:
        return "nil"
    if isinstance(v, str):
        return v if v.startswith(":") else '"' + v + '"'
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating, mp.mpf)):
        f = float(v)
        if f != f:
            return "##NaN"
        if f == float("inf"):
            return "##Inf"
        if f == float("-inf"):
            return "##-Inf"
        return repr(f)
    if isinstance(v, (list, tuple, np.ndarray)):
        return "[" + " ".join(edn(x) for x in v) + "]"
    if isinstance(v, dict):
        pad = " " * (indent + 1)
        items = [edn(k) + " " + edn(x, indent + 1) for k, x in v.items()]
        return "{" + ("\n" + pad).join(items) + "}"
    raise TypeError(type(v))


# ---------------------------------------------------------------- logarithmic

LOG_P = [0.1, 0.5, 0.9, 0.99, 0.999]
LOG_KMAX = 100000


def logarithmic_reference(p):
    """cdf(k) for k in 1..LOG_KMAX as a list of mpf, via the cumulative sum of p^j/j / (-log(1-p))."""
    mpp = mp.mpf(p)
    norm = -mp.log(1 - mpp)
    cum = []
    pw = mp.mpf(1)
    s = mp.mpf(0)
    for j in range(1, LOG_KMAX + 1):
        pw *= mpp
        s += pw / j
        cum.append(s / norm)
    return cum


def logarithmic():
    rng = random.Random(20251021)
    ks = list(range(1, 201)) + [230, 250, 300, 400, 500, 700, 1000, 1500, 2000, 3000, 5000, 7000, 10000,
                                 20000, 30000, 50000, 70000, 100000]
    while len(set(ks)) < 278:
        ks.append(int(math.exp(rng.uniform(math.log(200), math.log(LOG_KMAX)))))
    ks = sorted(set(ks))
    tails = [1e-6, 1e-3, 1 - 1e-4, 1 - 1e-6, 1 - 1e-9, 1 - 1e-12]
    out = []
    for p in LOG_P:
        cum = logarithmic_reference(p)
        cdf = [cum[k - 1] for k in ks]
        ps = [rng.random() for _ in range(750)] + tails
        icdf = []
        for q in ps:
            # smallest k with cdf(k) >= q
            lo, hi = 0, len(cum)
            mq = mp.mpf(q)
            while lo < hi:
                mid = (lo + hi) // 2
                if cum[mid] >= mq:
                    hi = mid
                else:
                    lo = mid + 1
            assert lo < len(cum), (p, q)
            icdf.append(lo + 1)
        out.append({":p": p, ":ks": ks, ":cdf": cdf, ":icdf-p": ps, ":icdf": icdf})
    return out


# ---------------------------------------------------------------- dirichlet

def dirichlet():
    rng = np.random.default_rng(101)
    out = []
    for alpha in ([2.0, 3.0, 4.0], [0.5, 0.5, 0.5], [1.0, 1.0], [5.0], [0.3, 1.2, 2.5, 4.0, 0.9]):
        pts = []
        for _ in range(5):
            x = rng.dirichlet(np.maximum(alpha, 0.5) * 2.0)
            x = np.array([float(v) for v in x])
            x[-1] = 1.0 - float(np.sum(x[:-1]))  # simplex up to the rounding of the last coordinate
            pts.append(x)
        d = st.dirichlet(alpha)
        out.append({":alpha": alpha,
                    ":points": [list(map(float, x)) for x in pts],
                    ":pdf": [float(d.pdf(x)) for x in pts],
                    ":lpdf": [float(d.logpdf(x)) for x in pts]})
    return out


# ---------------------------------------------------------------- multinomial

def multinomial_cdf(n, p, x):
    """P(X <= x componentwise), exact rational enumeration of the outcomes."""
    ps = [Fraction(v) for v in p]
    tot = sum(ps)
    ps = [v / tot for v in ps]
    total = Fraction(0)
    for c in itertools.product(*[range(0, min(n, xi) + 1) for xi in x]):
        if sum(c) != n:
            continue
        coef = math.factorial(n)
        for ci in c:
            coef //= math.factorial(ci)
        term = Fraction(coef)
        for ci, pi in zip(c, ps):
            term *= pi ** ci
        total += term
    return float(total)


def multinomial():
    out = []
    specs = [
        (10, [0.2, 0.3, 0.5], [[2, 3, 5], [10, 0, 0], [0, 5, 5], [3, 3, 4], [1, 1, 8]], [[2, 3, 5], [4, 4, 6], [10, 10, 10]]),
        (1, [0.25, 0.25, 0.5], [[1, 0, 0], [0, 0, 1]], [[1, 0, 1], [0, 0, 0]]),
        (20, [0.1, 0.2, 0.3, 0.4], [[2, 4, 6, 8], [20, 0, 0, 0], [5, 5, 5, 5]], [[5, 5, 5, 5]]),
        (6, [0.5, 0.0, 0.5], [[3, 0, 3], [6, 0, 0], [0, 0, 6], [3, 1, 2]], [[3, 0, 6], [2, 0, 4]]),
    ]
    for n, p, pts, cdf_pts in specs:
        d = st.multinomial(n, p)
        out.append({":n": n, ":ps": p, ":points": pts,
                    ":pmf": [float(d.pmf(x)) for x in pts],
                    ":cdf-points": cdf_pts,
                    ":cdf": [multinomial_cdf(n, p, x) for x in cdf_pts]})
    return out


# ---------------------------------------------------------------- multi-normal

def multi_normal():
    rng = np.random.default_rng(7)
    out = []
    for dim in range(1, 6):
        a = rng.standard_normal((dim, dim))
        s = a @ a.T + np.eye(dim)
        s = (s + s.T) / 2.0
        mu = rng.standard_normal(dim)
        pts = [mu + rng.standard_normal(dim) * 1.5 for _ in range(4)]
        d = st.multivariate_normal(mu, s)
        out.append({":means": mu, ":covariances": [list(map(float, r)) for r in s],
                    ":points": [list(map(float, x)) for x in pts],
                    ":pdf": [float(d.pdf(x)) for x in pts],
                    ":lpdf": [float(d.logpdf(x)) for x in pts]})
    far = []
    for mu, s, x in (([0.0, 0.0], [[1.0, 0.5], [0.5, 2.0]], [40.0, -40.0]),
                     ([1.0, -1.0, 0.5], [[2.0, 0.3, 0.1], [0.3, 1.0, 0.2], [0.1, 0.2, 0.5]], [30.0, 30.0, 30.0]),
                     ([0.0], [[1.0]], [60.0])):
        d = st.multivariate_normal(mu, s)
        far.append({":means": mu, ":covariances": s, ":point": x, ":lpdf": float(d.logpdf(x))})
    return {":interior": out, ":far": far}


# ---------------------------------------------------------------- inverse gaussian

def ig_cdf(mu, lam, x):
    mu, lam, x = mp.mpf(mu), mp.mpf(lam), mp.mpf(x)
    s = mp.sqrt(lam / x)
    return mp.ncdf(s * (x / mu - 1)) + mp.exp(2 * lam / mu) * mp.ncdf(-s * (x / mu + 1))


def inverse_gaussian():
    out = []
    for mu, lam in ((1.0, 1.0), (1.0, 100.0), (1.0, 354.0), (1.0, 355.0), (1.0, 400.0), (1.0, 1000.0),
                    (1.0, 10000.0), (0.01, 100.0), (5.0, 1775.0), (2.0, 200.0)):
        ratio = lam / mu
        xs = [float(mu * (1.0 + k / math.sqrt(ratio))) for k in (-4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 5, 8)]
        xs += [0.2 * mu, 0.5 * mu]
        xs = sorted(set(x for x in xs if x > 0.0))
        cdf = [ig_cdf(mu, lam, x) for x in xs]
        sci = st.invgauss(mu / lam, scale=lam)
        for x, c in zip(xs, cdf):
            v = float(sci.cdf(x))
            if math.isfinite(v) and v > 1e-300 and c > 1e-300:
                assert abs(v - float(c)) <= 1e-8 * float(c) + 1e-300, (mu, lam, x, v, float(c))
        out.append({":mu": mu, ":lambda": lam, ":xs": xs, ":cdf": cdf})
    return out


# ---------------------------------------------------------------- output

def main():
    data = {
        ":logarithmic": logarithmic(),
        ":dirichlet": dirichlet(),
        ":multinomial": multinomial(),
        ":multi-normal": multi_normal(),
        ":inverse-gaussian": inverse_gaussian(),
    }
    header = (";; Reference data for the native distributions (logarithmic, dirichlet, multinomial, multi-normal, inverse-gaussian).\n"
              ";; Generated by utils/fastmath/dev/generate_random_native_distributions_reference.py\n"
              ";; (mpmath %s, scipy %s, numpy %s). Do not edit by hand.\n" % (mp.__version__, scipy.__version__, np.__version__))
    body = "{" + "\n ".join(k + " " + edn(v, 1) for k, v in data.items()) + "}\n"
    with open(OUT, "w") as f:
        f.write(header + body)
    print("written", OUT)


if __name__ == "__main__":
    main()
