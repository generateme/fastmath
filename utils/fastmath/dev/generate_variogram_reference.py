"""Generator of reference values for `fastmath.kernel.variogram` (models and estimators).

Tools: mpmath 1.3.0 (80 digits; textbook formulas written independently of the Clojure code), scikit-gstat
1.0.24 (estimators `matheron`, `cressie`, `dowd`), numpy 2.5 (RandomState for the data, the Weibull quantile
`method='weibull'` for the Genton estimator, a re-implementation of the Gunst-Hartfield M-estimator), gstools 1.7.0
(cross-check of the TPL stable model, asserted in the generator, not stored).

Domain and why:
  * :matern - `matern` with nugget 0, psill 1, range 1: 31 orders (0.1 .. 200) x 23 lags (1e-300 .. 100). Before the fix
    the model gave NaN at small lags for large orders (0 * Inf) and wrong values from order 150.
    The semivariance is 1 - 2^(1-nu)/Gamma(nu) d^nu K_nu(d) with d = sqrt(2 nu) x. Values below 1e-80 are not
    resolved by 80 digits (stored as the rounded difference, within 1e-12 of the truth).
  * :bessel - `bessel` (BesselJ model) with nugget 0, psill 1, range 1: 16 orders (0 .. 300) x 16 lags (1e-300 .. 1000);
    the semivariance is 1 - Gamma(nu + 1) (2/x)^nu J_nu(x).
  * :hyperspherical and :superspherical - the orders inside the documented accuracy limits (error of the 2F1
    evaluation: hyperspherical dimension up to 60 for non-integer parameter d and up to 31 for odd dimensions,
    superspherical order up to 15 for integer orders), 16 values of h in [0.001, 1].
  * :models - the other models (exponential, gaussian, expower, rational, cauchy, power, tpower, hole, spherical,
    cubic, pentaspherical, circular, linear) for two parameter sets and 17 relative lags (tiny, fractions of the range,
    the range, beyond it). The gstat models are in the R generator.
  * :tplstable - `tplstable` (E_n exponential integral) at parameters where `fastmath.special/En` is exact:
    integer orders 1 + 2H/beta and H 0.8 / beta 1.9; the fractional orders are in :tplstable-fractional
    (they follow the fix of `En`).
  * :estimators - the estimators on 12 difference vectors (N 1 .. 40, one with an outlier, one with ties):
    gamma = half of the estimate of the variogram. Dowd: scikit-gstat uses the constant 2.198, the value stored is
    `dowd(|d|) / 2.198 * 2.198109338317728` (the constant of fastmath, the reciprocal of the median of chi-squared(1)).

Produces `test/resources/kernel/variogram_reference.edn`.

Run from anywhere:
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_variogram_reference.py
"""
from math import comb
from pathlib import Path

import mpmath as mp
import numpy as np
from skgstat import estimators as skest

mp.mp.dps = 80

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "test" / "resources" / "kernel" / "variogram_reference.edn"

DOWD_FASTMATH = 2.198109338317728


def edn_double(x):
    s = repr(float(x))
    assert "inf" not in s and "nan" not in s, s
    return s


def vec(xs):
    return "[" + " ".join(edn_double(x) for x in xs) + "]"


def mpf(x):
    return mp.mpf(x)


# ------------------------------------------------------------------ matern, bessel
MATERN_ORDERS = [0.1, 0.3, 0.5, 0.9, 1.0, 1.1, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0, 12.0, 15.0, 17.5, 19.9,
                 20.0, 20.1, 22.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 170.0, 171.5, 200.0]
MATERN_LAGS = [1e-300, 1e-100, 1e-30, 1e-15, 1e-8, 1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.3, 1.7,
               2.0, 2.5, 3.0, 5.0, 10.0, 100.0]
BESSEL_ORDERS = [0.0, 0.3, 0.5, 1.0, 1.5, 2.0, 3.5, 5.0, 10.0, 20.0, 50.0, 100.0, 150.0, 170.0, 200.0, 300.0]
BESSEL_LAGS = [1e-300, 1e-50, 1e-10, 1e-4, 1e-2, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 40.0, 100.0, 1000.0]


def matern_gamma(nu, x):
    nu = mpf(nu)
    d = mp.sqrt(2 * nu) * mpf(x)
    f = mpf(2) ** (1 - nu) / mp.gamma(nu) * d ** nu * mp.besselk(nu, d)
    return 1 - f


def bessel_gamma(nu, x):
    nu = mpf(nu)
    d = mpf(x)
    return 1 - mp.gamma(nu + 1) * (2 / d) ** nu * mp.besselj(nu, d)


def grid_block(key, orders, lags, fn):
    o, l, g = [], [], []
    for nu in orders:
        for x in lags:
            o.append(nu)
            l.append(x)
            g.append(fn(nu, x))
    return " :%s {:order %s :lag %s :gamma %s}" % (key, vec(o), vec(l), vec(g))


# ------------------------------------------------------------------ hyperspherical, superspherical
HYPER_ORDERS = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0, 11.0, 13.0, 15.0, 20.0, 21.0, 25.0, 30.0, 31.0, 40.0, 60.0]
SUPER_ORDERS = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0, 15.0]
SPH_H = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999, 1.0]


def spherical_family(d_of_order):
    def fn(order, h):
        d = d_of_order(mpf(order))
        den = mp.hyp2f1(mpf(1) / 2, d, mpf(3) / 2, 1)
        return mpf(h) * mp.hyp2f1(mpf(1) / 2, d, mpf(3) / 2, mpf(h) ** 2) / den
    return fn


# ------------------------------------------------------------------ simple models
PARAM_SETS = [(0.0, 1.0, 1.0), (0.3, 2.5, 7.5)]
REL_LAGS = [1e-9, 1e-3, 0.05, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 0.99, 1.0, 1.0000001, 1.2, 1.5, 3.0, 10.0]
BETA_SETS = {"power": [0.5, 1.0, 1.5, 2.0], "expower": [0.5, 1.0, 1.5, 2.0], "rational": [0.5, 2.0, 7.0],
             "cauchy": [0.5, 1.0, 4.0], "tpower": [0.5, 1.5, 3.0]}


def unit_model(name, h, beta):
    """Textbook correlation-type function g(h), h = x / range, the semivariance is nugget + psill * g."""
    h = mpf(h)
    one = mpf(1)
    if name == "exponential":
        return one - mp.exp(-h)
    if name == "gaussian":
        return one - mp.exp(-h ** 2)
    if name == "expower":
        return one - mp.exp(-h ** mpf(beta))
    if name == "rational":
        b = mpf(beta)
        return one - (one + h ** 2 / b) ** (-b)
    if name == "cauchy":
        return one - (one + h ** 2) ** (-mpf(beta))
    if name == "power":
        return h ** mpf(beta)
    if name == "tpower":
        return one if h > 1 else one - (one - h) ** mpf(beta)
    if name == "hole":
        return one - mp.sinc(mp.pi * h)
    if name == "spherical":
        return one if h > 1 else mpf(3) / 2 * h - h ** 3 / 2
    if name == "cubic":
        return one if h > 1 else 7 * h ** 2 - mpf(35) / 4 * h ** 3 + mpf(7) / 2 * h ** 5 - mpf(3) / 4 * h ** 7
    if name == "pentaspherical":
        return one if h > 1 else mpf(15) / 8 * h - mpf(5) / 4 * h ** 3 + mpf(3) / 8 * h ** 5
    if name == "circular":
        return one if h > 1 else 2 / mp.pi * (h * mp.sqrt(1 - h ** 2) + mp.asin(h))
    if name == "linear":
        return one if h > 1 else h
    raise ValueError(name)


MODELS = ["exponential", "gaussian", "expower", "rational", "cauchy", "power", "tpower", "hole", "spherical", "cubic",
          "pentaspherical", "circular", "linear"]


def models_block():
    lines = []
    for name in MODELS:
        for beta in BETA_SETS.get(name, [None]):
            for (nugget, psill, rng) in PARAM_SETS:
                lags = [r * rng for r in REL_LAGS]
                g = [nugget + psill * unit_model(name, mpf(x) / mpf(rng), beta) for x in lags]
                bs = "" if beta is None else " :beta %s" % edn_double(beta)
                lines.append("  {:model :%s :nugget %s :psill %s :range %s%s :lag %s :gamma %s}"
                             % (name, edn_double(nugget), edn_double(psill), edn_double(rng), bs, vec(lags), vec(g)))
    return " :models [\n" + "\n".join(lines) + "\n ]"


# ------------------------------------------------------------------ tplstable
TPL_REL = [1e-9, 1e-4, 0.01, 0.1, 0.3, 0.5, 0.8, 1.0, 2.0, 5.0, 20.0, 100.0]
TPL_EXACT = [(0.8, 1.9), (0.5, 1.0), (0.5, 0.5), (0.75, 1.5), (0.9, 0.6)]
TPL_FRACTIONAL = [(0.5, 1.2), (0.3, 0.8), (0.25, 1.0), (0.5, 2.0), (0.6, 1.7), (0.1, 1.3)]


def tpl_gamma(h_exp, beta, nugget, psill, rng, x):
    ha = 2 * mpf(h_exp) / mpf(beta)
    h = mpf(x) / mpf(rng)
    return nugget + psill * (1 - ha * mp.expint(1 + ha, h ** mpf(beta)))


def tpl_cross_check():
    import gstools as gs
    for (hurst, beta) in TPL_EXACT:
        model = gs.TPLStable(dim=2, var=1.0, len_scale=1.0, hurst=hurst, alpha=beta)
        for x in (0.05, 0.3, 1.0, 4.0):
            ref = float(tpl_gamma(hurst, beta, 0.0, 1.0, 1.0, x))
            got = float(model.variogram(x))
            assert abs(ref - got) < 1e-10, (hurst, beta, x, ref, got)


def tpl_block(key, sets):
    lines = []
    for (hurst, beta) in sets:
        for (nugget, psill, rng) in PARAM_SETS:
            lags = [r * rng for r in TPL_REL]
            g = [tpl_gamma(hurst, beta, nugget, psill, rng, x) for x in lags]
            lines.append("  {:H %s :beta %s :nugget %s :psill %s :range %s :lag %s :gamma %s}"
                         % (edn_double(hurst), edn_double(beta), edn_double(nugget), edn_double(psill),
                            edn_double(rng), vec(lags), vec(g)))
    return " :%s [\n" % key + "\n".join(lines) + "\n ]"


# ------------------------------------------------------------------ estimators
def gunst_hartfield(d):
    a = np.sqrt(np.abs(d))
    m = np.median(a)
    s = np.median(np.abs(a - m)) / 0.6745
    t = (a - m) / s
    inside = np.abs(t) < 4.0
    psi = np.where(inside, t * (16.0 - t ** 2) ** 2, 0.0)
    dpsi = np.where(inside, (16.0 - t ** 2) * (16.0 - 5.0 * t ** 2), 0.0)
    n = len(d)
    return (m + s * psi.sum() / dpsi.sum()) ** 4 / (0.457 + 0.494 / n + 0.045 / n ** 2)


def genton_weibull(d):
    n = len(d)
    ds = np.array([abs(d[i] - d[j]) for i in range(n) for j in range(i + 1, n)])
    q = comb(n // 2 + 1, 2) / len(ds)
    return (2.2191 * np.quantile(ds, q, method="weibull")) ** 2


def estimators_block():
    rs = np.random.RandomState(2024)
    cases = []
    for n in (1, 2, 3, 4, 5, 8, 10, 11, 18, 25, 40):
        cases.append(rs.normal(0.0, 1.5, n))
    outlier = rs.normal(0.0, 1.0, 12)
    outlier[3] = 25.0
    cases.append(outlier)
    cases.append(np.array([1.0, -1.0, 1.0, -1.0, 2.0, 2.0, -2.0]))  # ties
    lines = []
    for d in cases:
        n = len(d)
        a = np.abs(d)
        matheron = skest.matheron(d)
        cressie = skest.cressie(a)
        dowd = skest.dowd(a) / 2.198 * DOWD_FASTMATH
        parts = ["{:n %d :diffs %s :matheron %s :cressie %s :dowd %s" % (n, vec(d), edn_double(matheron),
                                                                         edn_double(cressie), edn_double(dowd))]
        if n >= 2:
            parts.append(" :genton-weibull %s" % edn_double(0.5 * genton_weibull(d)))
        if n >= 4:
            with np.errstate(all="ignore"):
                mr = 0.5 * gunst_hartfield(d)
            if np.isfinite(mr):  # undefined (NaN) when the median absolute deviation is zero (ties)
                parts.append(" :m-robust %s" % edn_double(mr))
        lines.append("  " + "".join(parts) + "}")
    return " :estimators [\n" + "\n".join(lines) + "\n ]"


def main():
    tpl_cross_check()
    blocks = [
        grid_block("matern", MATERN_ORDERS, MATERN_LAGS, matern_gamma),
        grid_block("bessel", BESSEL_ORDERS, BESSEL_LAGS, bessel_gamma),
        grid_block("hyperspherical", HYPER_ORDERS, SPH_H, spherical_family(lambda o: -(o - 1) / 2)),
        grid_block("superspherical", SUPER_ORDERS, SPH_H, spherical_family(lambda o: -o)),
        models_block(),
        tpl_block("tplstable", TPL_EXACT),
        tpl_block("tplstable-fractional", TPL_FRACTIONAL),
        estimators_block(),
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("{\n" + "\n".join(blocks) + "\n}\n")
    print("mpmath", mp.__version__, "wrote", OUT)


if __name__ == "__main__":
    main()
