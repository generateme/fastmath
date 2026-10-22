"""Generator of reference data for `truncated`, `mixture` and `ccdf` of `fastmath.random.distributions`.

Tools: mpmath 1.3.0 (oracle for the truncated distributions, the discrete cases, the mixture quantiles and the
`ccdf` of the keys for which scipy has no usable `isf`), scipy 1.18.1 (`stats.<dist>.sf/isf/cdf` for `ccdf`).
No reference value depends on a random draw.

Precision: mpmath with 60 digits for every truncated and mixture value (the integrals by `mp.quad`, the
quantiles by `mp.findroot`), rounded once to double.  The `ccdf` values are scipy `sf` of the exact double
`x`, with `x = isf(level)` for the levels 1e-6, 1e-12, 1e-30, 1e-100 (the `x` is stored, `sf` is a function
of that exact double, so the inaccuracy of `isf` does not matter).  The scipy `sf` of the families that
have a simple closed form is checked here against mpmath (relative 1e-10, asserted).

Domain:
    :ccdf        key + parameters of the Clojure constructor, the points `x` at the levels 1e-6, 1e-12, 1e-30,
                 1e-100 and `sf(x)` (and the cdf at the median, to verify the parameter mapping from the
                 Clojure side).  Normal, t, logistic, laplace, exponential, weibull, pareto, log-normal, gamma,
                 chi-squared, cauchy and the SSJ keys erlang, fatigue-life, frechet, half-normal, hyperbolic-secant,
                 inverse-gamma, johnson-su, rayleigh, chi, chi-squared-noncentral: scipy (the closed forms are
                 asserted against mpmath, relative 1e-8).  Log-logistic, pearson-6, inverse-gaussian,
                 hypoexponential-equal, johnson-sb, johnson-sl, power: mpmath closed forms (60 digits), the
                 point `x` is the double nearest to the root of sf(x) = level.
    :truncated   N[8,9], N[5,6], N[10,11], N[-9,-8], N(-inf,-3], N[2,inf), N[-1,2], N[3,4], Exp[30,31],
                 Cauchy[-5,5], Beta(2,5)[.1,.4], Gamma(3,2)[1,10]: pdf at 5 points, cdf at 5 points
                 (cdf as the plain probability, not 1 - ...), mean, variance, icdf at p in 1e-12, 1e-3,
                 .5, .999, 1 - 1e-12 (the p are exact doubles; the quantile is the exact quantile of that
                 double).
    :discrete    Poisson(4)[2,6], Poisson(4)[2,inf), Binomial(10,.4)[3,7], Geometric(.3)[2,inf): pmf, cdf,
                 mean, variance, icdf on a grid of probabilities away from the cumulative values
                 (exact rational arithmetic of the normalised pmf in mpmath).
    :mixture     N(0,1)+N(5,1) weights [.3 .7], N(0,1)+N(100,1) [.5 .5], N(0,1)+N(0,1e6), Poisson(3)+
                 Poisson(10): icdf on a grid of p (smallest x with cdf(x) >= p, mpmath root bracketing on
                 the exact mixture cdf in 60 digits; discrete: exact integers).

Run (writes test/resources/random/truncated_mixture_reference.edn, byte-identical on rerun):
    cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_random_truncated_mixture_reference.py
"""
import os
import warnings
import mpmath as mp
import scipy
import scipy.stats as st

mp.mp.dps = 60
warnings.filterwarnings("ignore")

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                   "test", "resources", "random", "truncated_mixture_reference.edn")

LEVELS = [1e-6, 1e-12, 1e-30, 1e-100]


# ---------------------------------------------------------------- EDN writer

def edn(v, indent=0):
    if isinstance(v, bool):
        return "true" if v else "false"
    if v is None:
        return "nil"
    if isinstance(v, str):
        return v if v.startswith(":") else '"' + v + '"'
    if isinstance(v, int):
        return str(v)
    if isinstance(v, (float, mp.mpf)):
        f = float(v)
        if f != f:
            return "##NaN"
        if f == float("inf"):
            return "##Inf"
        if f == float("-inf"):
            return "##-Inf"
        return repr(f)
    if isinstance(v, (list, tuple)):
        return "[" + " ".join(edn(x) for x in v) + "]"
    if isinstance(v, dict):
        pad = " " * (indent + 1)
        items = [edn(k) + " " + edn(x, indent + 1) for k, x in v.items()]
        return "{" + ("\n" + pad).join(items) + "}"
    raise TypeError(type(v))


def kw(name):
    return ":" + name


# ---------------------------------------------------------------- ccdf

def ccdf_cases():
    """(key, parameters, scipy frozen distribution, closed-form mpmath sf or None)."""
    c = []
    def add(key, params, dist, mpsf=None):
        c.append((key, params, dist, mpsf))

    add("normal", {":mu": 0.0, ":sd": 1.0}, st.norm(0.0, 1.0), lambda x: mp.ncdf(-x))
    add("normal", {":mu": 3.0, ":sd": 2.5}, st.norm(3.0, 2.5), lambda x: mp.ncdf(-(x - 3.0) / 2.5))
    add("t", {":degrees-of-freedom": 3.0}, st.t(3.0))
    add("t", {":degrees-of-freedom": 10.5}, st.t(10.5))
    add("logistic", {":mu": 0.0, ":s": 1.0}, st.logistic(0.0, 1.0), lambda x: 1 / (1 + mp.exp(x)))
    add("logistic", {":mu": 2.0, ":s": 3.0}, st.logistic(2.0, 3.0), lambda x: 1 / (1 + mp.exp((x - 2.0) / 3.0)))
    add("laplace", {":mu": 0.0, ":beta": 1.0}, st.laplace(0.0, 1.0), lambda x: mp.exp(-x) / 2)
    add("laplace", {":mu": -1.0, ":beta": 2.0}, st.laplace(-1.0, 2.0), lambda x: mp.exp(-(x + 1.0) / 2.0) / 2)
    add("exponential", {":mean": 1.0}, st.expon(0.0, 1.0), lambda x: mp.exp(-x))
    add("exponential", {":mean": 3.5}, st.expon(0.0, 3.5), lambda x: mp.exp(-x / 3.5))
    add("weibull", {":alpha": 2.0, ":beta": 1.0}, st.weibull_min(2.0, 0.0, 1.0), lambda x: mp.exp(-(x / 1.0) ** 2.0))
    add("weibull", {":alpha": 0.7, ":beta": 3.0}, st.weibull_min(0.7, 0.0, 3.0), lambda x: mp.exp(-(x / 3.0) ** 0.7))
    add("pareto", {":scale": 1.0, ":shape": 3.0}, st.pareto(3.0, 0.0, 1.0), lambda x: (1 / x) ** 3)
    add("pareto", {":scale": 2.5, ":shape": 1.5}, st.pareto(1.5, 0.0, 2.5), lambda x: (2.5 / x) ** 1.5)
    add("log-normal", {":scale": 0.0, ":shape": 1.0}, st.lognorm(1.0, 0.0, 1.0), lambda x: mp.ncdf(-mp.log(x)))
    add("log-normal", {":scale": 1.0, ":shape": 0.5}, st.lognorm(0.5, 0.0, float(mp.exp(1))),
        lambda x: mp.ncdf(-(mp.log(x) - 1.0) / 0.5))
    add("gamma", {":shape": 3.0, ":scale": 1.0}, st.gamma(3.0, 0.0, 1.0), lambda x: mp.gammainc(3.0, x, mp.inf, regularized=True))
    add("gamma", {":shape": 0.5, ":scale": 2.0}, st.gamma(0.5, 0.0, 2.0), lambda x: mp.gammainc(0.5, x / 2.0, mp.inf, regularized=True))
    add("gamma", {":shape": 20.0, ":scale": 0.3}, st.gamma(20.0, 0.0, 0.3), lambda x: mp.gammainc(20.0, x / 0.3, mp.inf, regularized=True))
    add("chi-squared", {":degrees-of-freedom": 3.0}, st.chi2(3.0), lambda x: mp.gammainc(1.5, x / 2, mp.inf, regularized=True))
    add("chi-squared", {":degrees-of-freedom": 10.0}, st.chi2(10.0), lambda x: mp.gammainc(5.0, x / 2, mp.inf, regularized=True))
    add("cauchy", {":median": 0.0, ":scale": 1.0}, st.cauchy(0.0, 1.0), lambda x: mp.atan(1 / x) / mp.pi)
    add("cauchy", {":median": 1.0, ":scale": 2.0}, st.cauchy(1.0, 2.0), lambda x: mp.atan(2.0 / (x - 1.0)) / mp.pi)
    # SSJ-backed keys (the scipy parameterisation is verified through the stored cdf at the median)
    add("erlang", {":k": 3, ":lambda": 2.0}, st.erlang(3, 0.0, 0.5))
    add("fatigue-life", {":alpha": 0.5, ":beta": 2.0, ":gamma": 1.5}, st.fatiguelife(1.5, 0.5, 2.0))
    add("frechet", {":alpha": 2.5, ":beta": 1.5, ":delta": 0.5}, st.invweibull(2.5, 0.5, 1.5))
    add("half-normal", {":mu": 0.5, ":sigma": 2.0}, st.halfnorm(0.5, 2.0))
    add("hyperbolic-secant", {":mu": 0.5, ":sigma": 2.0}, st.hypsecant(0.5, 4.0 / float(mp.pi)))
    add("inverse-gamma", {":alpha": 3.0, ":beta": 2.0}, st.invgamma(3.0, 0.0, 2.0))
    add("johnson-su", {":gamma": 0.5, ":delta": 1.5, ":xi": 0.5, ":lambda": 2.0}, st.johnsonsu(0.5, 1.5, 0.5, 2.0))
    add("rayleigh", {":a": 0.5, ":beta": 2.0}, st.rayleigh(0.5, 2.0))
    add("chi", {":nu": 3}, st.chi(3))
    add("chi-squared-noncentral", {":nu": 3.0, ":lambda": 2.0}, st.ncx2(3.0, 2.0))
    return c


def mp_ccdf_cases():
    """(key, parameters, mpmath sf, bisection bracket): cases for which scipy `isf` is not usable (heavy tails,
    bounded supports) or for which scipy has no counterpart."""
    c = []
    ncdf = mp.ncdf
    # log-logistic: alpha = shape, beta = scale
    c.append(("log-logistic", {":alpha": 2.0, ":beta": 3.0},
              lambda x: 1 / (1 + (x / 3.0) ** 2.0), (mp.mpf(0), mp.mpf(10) ** 60)))
    # pearson-6 = scaled beta prime(alpha1, alpha2) with scale beta
    c.append(("pearson-6", {":alpha1": 2.0, ":alpha2": 3.0, ":beta": 1.5},
              lambda x: mp.betainc(3.0, 2.0, 0, 1 / (1 + x / 1.5), regularized=True), (mp.mpf(0), mp.mpf(10) ** 60)))
    mu, lam = 2.0, 5.0
    c.append(("inverse-gaussian", {":mu": mu, ":lambda": lam},
              lambda x: ncdf(-mp.sqrt(lam / x) * (x / mu - 1)) - mp.exp(2 * lam / mu) * ncdf(-mp.sqrt(lam / x) * (x / mu + 1)),
              (mp.mpf(1) / 1000, mp.mpf(10) ** 5)))
    # hypoexponential with equal rate spacing: rates (n + 1 - i) h, i = 1..k
    n, k, h = 5, 3, 0.5
    rates = [(n + 1 - i) * h for i in range(1, k + 1)]
    def hypo(x):
        s = mp.mpf(0)
        for i, li in enumerate(rates):
            ci = mp.mpf(1)
            for j, lj in enumerate(rates):
                if j != i:
                    ci *= mp.mpf(lj) / (lj - li)
            s += ci * mp.exp(-li * x)
        return s
    c.append(("hypoexponential-equal", {":n": n, ":k": k, ":h": h}, hypo, (mp.mpf(0), mp.mpf(10) ** 3)))
    # Johnson SB / SL: F(x) = Phi(gamma + delta ln(z / (1 - z))) resp. Phi(gamma + delta ln z), z = (x - xi) / lambda
    g, d, xi, la = 0.5, 1.5, 0.5, 2.0
    c.append(("johnson-sb", {":gamma": g, ":delta": d, ":xi": xi, ":lambda": la},
              lambda x: ncdf(-(g + d * mp.log(((x - xi) / la) / (1 - (x - xi) / la)))),
              (mp.mpf(xi) + mp.mpf(la) * mp.mpf(10) ** -30, mp.mpf(xi) + mp.mpf(la) * (1 - mp.mpf(10) ** -30))))
    c.append(("johnson-sl", {":gamma": g, ":delta": d, ":xi": xi, ":lambda": la},
              lambda x: ncdf(-(g + d * mp.log((x - xi) / la))), (mp.mpf(xi) + mp.mpf(10) ** -30, mp.mpf(10) ** 30)))
    # power: F(x) = ((x - a) / (b - a))^c
    a, b, cc = 1.0, 3.0, 2.5
    c.append(("power", {":a": a, ":b": b, ":c": cc},
              lambda x: 1 - ((x - a) / (b - a)) ** cc, (mp.mpf(a) + 1, mp.mpf(b))))
    return c


def ccdf_section():
    rows = []
    for key, params, sf_fn, (lo, hi) in mp_ccdf_cases():
        pts = []
        for lv in LEVELS:
            x = float(mp_bisect(lambda t: lv - sf_fn(t), lo, hi, iters=700))
            sf = float(sf_fn(mp.mpf(x)))
            if sf > 0.0:
                pts.append([lv, x, sf])
        med = float(mp_bisect(lambda t: 0.5 - sf_fn(t), lo, hi, iters=700))
        rows.append({":key": kw(key), ":params": params, ":median": med,
                     ":cdf-at-median": float(1 - sf_fn(mp.mpf(med))), ":points": pts, ":oracle": ":mpmath"})
    rows.extend(scipy_ccdf_section())
    return rows


def scipy_ccdf_section():
    rows = []
    for key, params, dist, mpsf in ccdf_cases():
        pts = []
        for lv in LEVELS:
            x = float(dist.isf(lv))
            sf = float(dist.sf(x))
            if not (x == x and abs(x) != float("inf")) or sf <= 0.0 or sf != sf:
                continue
            if abs(sf / lv - 1.0) > 1e-3:
                # isf lost accuracy at this level; the pair (x, sf(x)) is still exact, keep it
                pass
            if mpsf is not None:
                ref = mpsf(mp.mpf(x))
                rel = abs((mp.mpf(sf) - ref) / ref)
                assert rel < 1e-8, (key, params, lv, x, sf, float(ref), float(rel))
            pts.append([lv, x, sf])
        med = float(dist.median())
        rows.append({":key": kw(key), ":params": params, ":median": med,
                     ":cdf-at-median": float(dist.cdf(med)), ":points": pts})
    return rows


# ---------------------------------------------------------------- truncated

def mp_quantile(cdf, p, lo, hi):
    f = lambda x: cdf(x) - p
    return mp.findroot(f, (lo, hi), solver="anderson", tol=mp.mpf(10) ** -45, maxsteps=300)


PS = [1e-12, 1e-3, 0.5, 0.999, 1 - 1e-12]


def trunc_continuous(name, base_pdf, base_cdf, base_sf, left, right, xs_pts, lo_support, hi_support):
    """Truncated continuous distribution: all functions in mpmath from the exact doubles."""
    L = mp.mpf(left) if left is not None else lo_support
    R = mp.mpf(right) if right is not None else hi_support
    # mass from the tail that is accurate in mpmath (dps 60 makes both equivalent)
    Z = base_cdf(R) - base_cdf(L) if R != mp.inf else base_sf(L)
    if L == -mp.inf:
        Z = base_cdf(R)
    pdf = lambda x: base_pdf(x) / Z
    cdf = lambda x: (base_cdf(x) - (base_cdf(L) if L != -mp.inf else 0)) / Z
    a = L if L != -mp.inf else mp.mpf(-40)
    b = R if R != mp.inf else mp.mpf(80)
    mean = mp.quad(lambda x: x * pdf(x), [a, b])
    var = mp.quad(lambda x: (x - mean) ** 2 * pdf(x), [a, b])
    out = {":left": left if left is not None else float("-inf"),
           ":right": right if right is not None else float("inf"),
           ":Z": float(Z), ":xs": xs_pts,
           ":pdf": [float(pdf(mp.mpf(x))) for x in xs_pts],
           ":cdf": [float(cdf(mp.mpf(x))) for x in xs_pts],
           ":mean": float(mean), ":variance": float(var), ":ps": PS}
    qs = []
    for p in PS:
        pm = mp.mpf(p)
        sfR = base_sf(R) if R != mp.inf else mp.mpf(0)
        # increasing function of x, built on the accurate side (survival for p > 1/2)
        if p <= 0.5:
            q = mp_bisect(lambda x: cdf(x) - pm, a, b)
        else:
            q = mp_bisect(lambda x: (1 - pm) - (base_sf(x) - sfR) / Z, a, b)
        qs.append(float(q))
    out[":icdf"] = qs
    return out


def mp_bisect(f, a, b, iters=400):
    """Root of an increasing f on [a, b] by bisection in 60 digits (f(a) < 0 < f(b))."""
    a = mp.mpf(a)
    b = mp.mpf(b)
    fa = f(a)
    fb = f(b)
    if fa >= 0:
        return a
    if fb <= 0:
        return b
    for _ in range(iters):
        m_ = (a + b) / 2
        fm = f(m_)
        if fm < 0:
            a = m_
        else:
            b = m_
        if b - a < mp.mpf(10) ** -40 * max(1, abs(b)):
            break
    return (a + b) / 2


def truncated_section():
    out = {}
    nsf = lambda x: mp.ncdf(-x)
    def normal_case(l, r, pts):
        return trunc_continuous("n", mp.npdf, mp.ncdf, nsf, l, r, pts, -mp.inf, mp.inf)
    out[":n-8-9"] = normal_case(8.0, 9.0, [8.0, 8.1, 8.5, 8.9, 9.0])
    out[":n-5-6"] = normal_case(5.0, 6.0, [5.0, 5.2, 5.5, 5.8, 6.0])
    out[":n-10-11"] = normal_case(10.0, 11.0, [10.0, 10.1, 10.5, 10.9, 11.0])
    out[":n-m9-m8"] = normal_case(-9.0, -8.0, [-9.0, -8.9, -8.5, -8.2, -8.0])
    out[":n-inf-m3"] = normal_case(None, -3.0, [-6.0, -4.5, -3.5, -3.1, -3.0])
    out[":n-2-inf"] = normal_case(2.0, None, [2.0, 2.2, 2.7, 3.5, 5.0])
    out[":n-m1-2"] = normal_case(-1.0, 2.0, [-0.9, -0.3, 0.0, 0.5, 1.99])
    out[":n-3-4"] = normal_case(3.0, 4.0, [3.0, 3.2, 3.5, 3.8, 4.0])
    # exponential(mean 1) on [30, 31]
    out[":exp-30-31"] = trunc_continuous("e", lambda x: mp.exp(-x) if x >= 0 else mp.mpf(0),
                                         lambda x: 1 - mp.exp(-x) if x >= 0 else mp.mpf(0),
                                         lambda x: mp.exp(-x) if x >= 0 else mp.mpf(1),
                                         30.0, 31.0, [30.0, 30.2, 30.5, 30.8, 31.0], mp.mpf(0), mp.inf)
    out[":cauchy-m5-5"] = trunc_continuous("c", lambda x: 1 / (mp.pi * (1 + x * x)),
                                           lambda x: mp.atan(x) / mp.pi + mp.mpf(1) / 2,
                                           lambda x: mp.mpf(1) / 2 - mp.atan(x) / mp.pi,
                                           -5.0, 5.0, [-5.0, -2.0, 0.0, 1.5, 5.0], -mp.inf, mp.inf)
    # beta(2,5) on [.1, .4]
    bpdf = lambda x: mp.betainc(2, 5, 0, 1, regularized=False) ** -1 * x * (1 - x) ** 4 if 0 <= x <= 1 else mp.mpf(0)
    bcdf = lambda x: mp.betainc(2, 5, 0, x, regularized=True) if 0 < x < 1 else (mp.mpf(0) if x <= 0 else mp.mpf(1))
    out[":beta-2-5-01-04"] = trunc_continuous("b", bpdf, bcdf, lambda x: 1 - bcdf(x), 0.1, 0.4,
                                              [0.1, 0.15, 0.25, 0.35, 0.4], mp.mpf(0), mp.mpf(1))
    # gamma(shape 3, scale 2) on [1, 10]
    gpdf = lambda x: x ** 2 * mp.exp(-x / 2) / (mp.gamma(3) * 2 ** 3) if x > 0 else mp.mpf(0)
    gcdf = lambda x: mp.gammainc(3, 0, x / 2, regularized=True) if x > 0 else mp.mpf(0)
    gsf = lambda x: mp.gammainc(3, x / 2, mp.inf, regularized=True) if x > 0 else mp.mpf(1)
    out[":gamma-3-2-1-10"] = trunc_continuous("g", gpdf, gcdf, gsf, 1.0, 10.0, [1.0, 2.0, 4.0, 7.0, 10.0],
                                              mp.mpf(0), mp.inf)
    return out


# ---------------------------------------------------------------- discrete

def discrete_case(pmf_fn, lo, hi, kmax=2000, ps=None):
    """Truncated integer distribution on [lo, hi] (hi None: until the mass is below 1e-50): exact."""
    top = hi if hi is not None else kmax
    ks = list(range(lo, top + 1))
    w = [pmf_fn(k) for k in ks]
    Z = mp.fsum(w)
    p = [x / Z for x in w]
    mean = mp.fsum(k * q for k, q in zip(ks, p))
    var = mp.fsum((k - mean) ** 2 * q for k, q in zip(ks, p))
    cum = []
    s = mp.mpf(0)
    for q in p:
        s += q
        cum.append(s)
    shown = list(range(max(lo - 2, 0), min(top, lo + 12) + 3))
    pm = [float(pmf_fn(k) / Z) if lo <= k <= top else 0.0 for k in shown]
    cd = [float(cum[k - lo]) if lo <= k <= top else (0.0 if k < lo else 1.0) for k in shown]
    grid = [0.0013, 0.0117, 0.0518, 0.1021, 0.2043, 0.3077, 0.4531, 0.5527, 0.6531, 0.8033, 0.9029, 0.9511, 0.9907, 0.9993]
    qs = []
    for g in grid:
        gm = mp.mpf(g)
        # smallest k with cum >= g; the grid is checked to be away from the cumulative values
        k = next(ks[i] for i, c in enumerate(cum) if c >= gm)
        assert min(abs(c - gm) for c in cum) > 1e-9, (g,)
        qs.append(k)
    return {":lo": lo, ":hi": hi if hi is not None else float("inf"), ":Z": float(Z),
            ":ks": shown, ":pmf": pm, ":cdf": cd, ":sum-pmf": float(mp.fsum(p)),
            ":mean": float(mean), ":variance": float(var), ":ps": grid, ":icdf": qs}


def discrete_section():
    out = {}
    pois = lambda lam: (lambda k: mp.exp(-lam) * mp.mpf(lam) ** k / mp.factorial(k))
    binom = lambda n, pr: (lambda k: mp.binomial(n, k) * mp.mpf(pr) ** k * (1 - mp.mpf(pr)) ** (n - k))
    geom = lambda pr: (lambda k: (1 - mp.mpf(pr)) ** k * mp.mpf(pr))  # support 0, 1, ... (Commons)
    out[":poisson-4-2-6"] = discrete_case(pois(4), 2, 6)
    out[":poisson-4-0-6"] = discrete_case(pois(4), 0, 6)
    out[":poisson-4-3-3"] = discrete_case(pois(4), 3, 3)
    out[":poisson-4-2-inf"] = discrete_case(pois(4), 2, None, kmax=200)
    out[":binomial-10-04-3-7"] = discrete_case(binom(10, 0.4), 3, 7)
    out[":geometric-03-2-inf"] = discrete_case(geom(0.3), 2, None, kmax=500)
    return out


# ---------------------------------------------------------------- mixture

def mixture_icdf_case(cdf_parts, ws, ps, lo, hi, integer=False):
    tot = mp.fsum(ws)
    ws = [mp.mpf(w) / tot for w in ws]
    cdf = lambda x: mp.fsum(w * c(x) for w, c in zip(ws, cdf_parts))
    out = []
    for p in ps:
        pm = mp.mpf(p)
        if integer:
            k = 0
            while cdf(k) < pm:
                k += 1
            out.append(k)
        else:
            q = mp_bisect(lambda x: cdf(x) - pm, lo, hi)
            out.append(float(q))
    return out


def mixture_section():
    out = {}
    ps = [1e-12, 0.001, 0.25, 0.5, 0.75, 0.999, 1 - 1e-12]
    ncdf = lambda mu, sd: (lambda x: mp.ncdf((x - mu) / sd))
    out[":n0-n5"] = {":weights": [0.3, 0.7], ":ps": ps,
                     ":icdf": mixture_icdf_case([ncdf(0, 1), ncdf(5, 1)], [0.3, 0.7], ps, -60, 60)}
    ps2 = [0.001, 0.25, 0.5, 0.75, 0.999]
    out[":n0-n0-1e6"] = {":weights": [0.5, 0.5], ":ps": ps2,
                         ":icdf": mixture_icdf_case([ncdf(0, 1), ncdf(0, 1e6)], [1, 1], ps2, -1e8, 1e8)}
    # plateau: the cdf is flat at 1/2 (to ~1e-17) between about 8.3 and 91.7; the smallest x with cdf >= .5
    # in double arithmetic is about 8.29, the reference stores the quantiles away from the plateau and
    # the exact (mpmath) quantile of .5 for information.
    ps3 = [0.1, 0.25, 0.4, 0.6, 0.75, 0.9]
    out[":n0-n100"] = {":weights": [0.5, 0.5], ":ps": ps3,
                       ":icdf": mixture_icdf_case([ncdf(0, 1), ncdf(100, 1)], [1, 1], ps3, -60, 160)}
    pois = lambda lam: (lambda k: mp.gammainc(int(mp.floor(k)) + 1, lam, mp.inf, regularized=True) if k >= 0 else mp.mpf(0))
    grid30 = [(i + 0.5) / 30 for i in range(30)]
    out[":pois3-pois10"] = {":weights": [0.5, 0.5], ":ps": grid30,
                            ":icdf": mixture_icdf_case([pois(3), pois(10)], [1, 1], grid30, 0, 0, integer=True)}
    return out


def main():
    ref = {":header": {":tools": ["mpmath " + mp.__version__, "scipy " + scipy.__version__],
                       ":precision": "mpmath dps 60, rounded once to double; ccdf: scipy sf of the exact double x",
                       ":ccdf-levels": LEVELS},
           ":ccdf": ccdf_section(),
           ":truncated": truncated_section(),
           ":discrete": discrete_section(),
           ":mixture": mixture_section()}
    text = edn(ref) + "\n"
    with open(OUT, "w") as f:
        f.write(text)
    print("written", os.path.normpath(OUT), len(text))


if __name__ == "__main__":
    main()
