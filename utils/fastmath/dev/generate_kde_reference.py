"""Generator of reference data for the kernels and the exact evaluation of `fastmath.kernel.density`.

Tools: mpmath 1.3.0 (oracle, dps 40 for the sums, dps 30 for the quadratures), numpy 2.5.2 (second
implementation, double precision), scipy 1.18.1 `stats.gaussian_kde`, scikit-learn 1.9.0
`KernelDensity`, statsmodels 0.15.0 `KDEUnivariate` (`fft=False`, exact sum). R `stats::density` is
not used (it bins).

Precision: every density is the direct sum  f(x) = 1/(n h) * sum_i K((x - x_i)/h)  over ALL points (no
cut-off, no binning) evaluated in mpmath with 40 digits from the exact binary values of x, x_i and h,
rounded once to double.  Densities that underflow are stored as 0.0.  Kernel values, `k2`, `x2k` and the
kernel radii are mpmath values (dps 30 quadrature) rounded to double.

Kernels (the 14 standard forms of `fastmath.kernel.density`; compact kernels on [-1, 1] and integrating
to 1, gaussian with sd 1; the integral of every kernel is stored as `:integral` and is 1 to 1e-29):
    uniform 1/2 (closed support: K(+-1) = 1/2),  triangular 1 - |x|,  epanechnikov 3/4 (1 - x^2),
    quartic 15/16 (1 - x^2)^2,  triweight 35/32 (1 - x^2)^3,  tricube 70/81 (1 - |x|^3)^3,
    gaussian exp(-x^2/2)/sqrt(2 pi),  cosine pi/4 cos(pi x/2),  logistic 1/(e^x + 2 + e^-x),
    sigmoid 2/pi / (e^x + e^-x),  silverman 1/2 exp(-|x|/sqrt 2) sin(|x|/sqrt 2 + pi/4) (4th order, K(+-Inf) = 0),
    wigner 2/pi sqrt(1 - x^2),  cauchy 2/(pi (1 + 4 x^2)) (the code's 1/(pi/2 (1 + (2x)^2))),
    laplace exp(-|x|)/2.

Domain: 13 datasets x 14 kernels = 182 cases, each with an explicit bandwidth `h` and 240 evaluation
points: the data points, points inside the gaps (the KE-01 failure region: within the kernel radius to the
LEFT of the right neighbour but beyond it from the left one), at +-R h from the extreme points and +-1e-6
relatively beside them (+-1e-3 for wigner: the square root at the edge amplifies the rounding of 1 - u^2
in double arithmetic, 2.5e-9 of the maximum at the exact edge), tails, far outside (1e3 h and 1e5 h) and seeded random points around the data.  `R`
is the support radius (1) for compact kernels and the radius where |K| <= 1e-10 for the others (30 for
cauchy, whose true radius is 39894).  For `uniform` (jump) and `wigner` (square root singularity) the points
within 1e-9 of |x - x_i| = h are dropped (not decidable in double arithmetic).  Datasets: single point, the two KE-01 anchors
(epanechnikov [0 10] h 1 at 9.5 = 0.28125; gaussian [0 100] h 1 at 98 = 0.026995...), 5/30/200 normal points,
a bimodal sample with a wide gap, a sample with many ties, a tight cluster with a far outlier, equal spacing,
a heavy tail (t with 2 df), a large offset (1e6) and two clusters 100 h apart.  Section `:anchors` repeats the
KE-01 values in exact form.

Why: KE-01 (the forward scan of `kde-` returned 0 or too small densities in every gap wider than the cut-off
span), KE-06/07 (cauchy cut-off), KE-08, and the work packages 1.2, 1.3, 1.6 of the Kernel Density audit.

Cross-checks (written in the header of the produced file with the measured maxima): numpy double
direct sum for 14 kernels; scipy `gaussian_kde` (bw_method = h / sd, gaussian); scikit-learn `KernelDensity`
(gaussian, tophat = uniform, epanechnikov, cosine, linear = triangular, exponential = laplace); statsmodels
`KDEUnivariate` where the kernel is defined there.  Differences are |other - oracle| / max_x oracle(x).

Produces `test/resources/kernel/kde_reference.edn`:
    {:kernels [...]                  keyword per kernel, in this order
     :anchors [{:id :kernel :data :h :x :density}]
     :kernel-values {:x [...] :values {kernel [...]}}      K at each x of :x (also negative and +-Inf)
     :kernel-constants {kernel {:integral :k2 :x2k :delta0 :nrd-scale :efficiency :support-radius
                                :radius}}
     :cases [{:id :dataset :kernel :data :h :x :density}]}
`:radius` is the largest |x| with |K(x)| > 1e-10 (for compact kernels 1.0); `:support-radius` 1.0 or
##Inf; `:delta0 = k2^0.2 x2k^-0.4`, `:nrd-scale = delta0 (0.5/sqrt(pi))^-0.2`, `:efficiency = k2 sqrt(x2k)`
(Inf/0 for silverman (x2k = 0) and cauchy (x2k = Inf) as in the old code).

Run from anywhere (a rerun gives a byte-identical file):
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_kde_reference.py
"""
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
import scipy.stats
import sklearn
import sklearn.neighbors
import statsmodels
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "test" / "resources" / "kernel" / "kde_reference.edn"

mp.mp.dps = 40
PI = mp.pi

KERNELS = ["uniform", "triangular", "epanechnikov", "quartic", "triweight", "tricube", "gaussian",
           "cosine", "logistic", "sigmoid", "silverman", "wigner", "cauchy", "laplace"]
COMPACT = {"uniform", "triangular", "epanechnikov", "quartic", "triweight", "tricube", "cosine", "wigner"}


# ------------------------------------------------------------------ EDN writer
class Kw(str):
    pass


def fmt(v):
    if v is None:
        return "nil"
    if isinstance(v, (bool, np.bool_)):
        return "true" if v else "false"
    if isinstance(v, Kw):
        return ":" + str(v)
    if isinstance(v, str):
        return json.dumps(v)
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        v = float(v)
        if math.isnan(v):
            return "##NaN"
        if math.isinf(v):
            return "##Inf" if v > 0 else "##-Inf"
        return repr(v)
    if isinstance(v, (list, tuple, np.ndarray)):
        return "[" + " ".join(fmt(e) for e in v) + "]"
    if isinstance(v, dict):
        return "{" + " ".join(fmt(k) + " " + fmt(e) for k, e in v.items()) + "}"
    raise TypeError(type(v))


def emit(v, depth=0):
    """Top-level map: one key per line; a list of maps: one map per line; a map of maps: one entry per line."""
    if isinstance(v, dict) and depth <= 1:
        ind = "  " * (depth + 1)
        lines = []
        for k, e in v.items():
            lines.append(ind + fmt(k) + " " + emit(e, depth + 1))
        return "{\n" + "\n".join(lines) + "}" if depth == 0 else "{\n" + "\n".join(lines) + "}"
    if isinstance(v, list) and v and all(isinstance(e, dict) for e in v):
        ind = "  " * (depth + 1)
        return "[\n" + "\n".join(ind + fmt(e) for e in v) + "]"
    return fmt(v)


# ------------------------------------------------------------------ oracle kernels (mpmath)
def K_mp(name, u):
    u = mp.mpf(u)
    a = abs(u)
    if name == "uniform":
        return mp.mpf(1) / 2 if a <= 1 else mp.mpf(0)
    if name == "triangular":
        return 1 - a if a <= 1 else mp.mpf(0)
    if name == "epanechnikov":
        return mp.mpf(3) / 4 * (1 - u * u) if a <= 1 else mp.mpf(0)
    if name == "quartic":
        return mp.mpf(15) / 16 * (1 - u * u) ** 2 if a <= 1 else mp.mpf(0)
    if name == "triweight":
        return mp.mpf(35) / 32 * (1 - u * u) ** 3 if a <= 1 else mp.mpf(0)
    if name == "tricube":
        return mp.mpf(70) / 81 * (1 - a ** 3) ** 3 if a <= 1 else mp.mpf(0)
    if name == "gaussian":
        return mp.exp(-u * u / 2) / mp.sqrt(2 * PI)
    if name == "cosine":
        return PI / 4 * mp.cos(PI / 2 * u) if a <= 1 else mp.mpf(0)
    if name == "logistic":
        if a == mp.inf:
            return mp.mpf(0)
        return 1 / (mp.exp(a) + 2 + mp.exp(-a))
    if name == "sigmoid":
        if a == mp.inf:
            return mp.mpf(0)
        return (2 / PI) / (mp.exp(a) + mp.exp(-a))
    if name == "silverman":
        if a == mp.inf:
            return mp.mpf(0)
        t = a / mp.sqrt(2)
        return mp.exp(-t) * mp.sin(t + PI / 4) / 2
    if name == "wigner":
        return 2 / PI * mp.sqrt(1 - u * u) if a <= 1 else mp.mpf(0)
    if name == "cauchy":
        if a == mp.inf:
            return mp.mpf(0)
        return 2 / (PI * (1 + 4 * u * u))
    if name == "laplace":
        if a == mp.inf:
            return mp.mpf(0)
        return mp.exp(-a) / 2
    raise KeyError(name)


def K_np(name, u):
    """Independent double precision (numpy) textbook formulas."""
    a = np.abs(u)
    inb = a <= 1
    pi = np.pi
    if name == "uniform":
        return np.where(inb, 0.5, 0.0)
    if name == "triangular":
        return np.where(inb, 1 - a, 0.0)
    if name == "epanechnikov":
        return np.where(inb, 0.75 * (1 - u * u), 0.0)
    if name == "quartic":
        return np.where(inb, 15 / 16 * (1 - u * u) ** 2, 0.0)
    if name == "triweight":
        return np.where(inb, 35 / 32 * (1 - u * u) ** 3, 0.0)
    if name == "tricube":
        return np.where(inb, 70 / 81 * (1 - a ** 3) ** 3, 0.0)
    if name == "gaussian":
        return np.exp(-0.5 * u * u) / np.sqrt(2 * pi)
    if name == "cosine":
        return np.where(inb, pi / 4 * np.cos(pi / 2 * u), 0.0)
    with np.errstate(over="ignore"):
        if name == "logistic":
            return 1 / (np.exp(a) + 2 + np.exp(-a))
        if name == "sigmoid":
            return (2 / pi) / (np.exp(a) + np.exp(-a))
    if name == "silverman":
        t = a / np.sqrt(2)
        return 0.5 * np.exp(-t) * np.sin(t + pi / 4)
    if name == "wigner":
        return np.where(inb, 2 / pi * np.sqrt(np.clip(1 - u * u, 0, None)), 0.0)
    if name == "cauchy":
        return 2 / (pi * (1 + 4 * u * u))
    if name == "laplace":
        return 0.5 * np.exp(-a)
    raise KeyError(name)


def density_mp(name, data, h, xs):
    """Exact (no cut-off) direct sum, mpmath 40 digits, one rounding to double."""
    n = len(data)
    d_mp = [mp.mpf(float(d)) for d in data]
    h_mp = mp.mpf(float(h))
    compact = name in COMPACT
    out = []
    for x in xs:
        xm = mp.mpf(float(x))
        terms = []
        for d in d_mp:
            u = (xm - d) / h_mp
            if compact and abs(u) > 1:
                continue
            terms.append(K_mp(name, u))
        out.append(float(mp.fsum(terms) / (n * h_mp)) if terms else 0.0)
    return out


def density_np(name, data, h, xs):
    data = np.asarray(data, dtype=float)
    xs = np.asarray(xs, dtype=float)
    u = (xs[:, None] - data[None, :]) / h
    return K_np(name, u).sum(1) / (len(data) * h)


# ------------------------------------------------------------------ constants
def tail_radius(name):
    """Largest x with |K(x)| > 1e-10 (compact kernels: the support radius 1)."""
    if name in COMPACT:
        return mp.mpf(1)
    thr = mp.mpf("1e-10")
    f = lambda x: abs(K_mp(name, x)) - thr
    if name == "silverman":
        # oscillating tail: scan downward from the envelope bound 0.5 exp(-x/sqrt 2) <= 1e-10
        dx = mp.mpf("0.001")
        x = mp.mpf(40)
        while f(x) <= 0:
            x -= dx
        lo, hi = x, x + dx
    else:
        lo, hi = mp.mpf(1), mp.mpf(10) ** 7
        assert f(lo) > 0 and f(hi) < 0
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def kernel_constants(name):
    mp.mp.dps = 30
    pts = [0, 1] if name in COMPACT else [0, 5, 10, 20, 40, 80, 160, mp.inf]
    if name == "cauchy":
        pts = [0, 1, 10, 100, 1000, mp.inf]
    F = lambda u: K_mp(name, u)
    integral = 2 * mp.quad(F, pts)
    k2 = 2 * mp.quad(lambda u: F(u) ** 2, pts)
    if name == "cauchy":
        x2k = mp.inf
        assert mp.quad(lambda u: u * u * F(u), [0, 1, 10, 100, 1000]) > mp.mpf(10) ** 2
    else:
        x2k = 2 * mp.quad(lambda u: u * u * F(u), pts)
        if abs(x2k) < mp.mpf("1e-25"):
            assert name == "silverman"
            x2k = mp.mpf(0)
    mp.mp.dps = 40
    radius = tail_radius(name)
    if x2k == 0:
        delta0, nscale, eff = mp.inf, mp.inf, mp.mpf(0)
    elif x2k == mp.inf:
        delta0, nscale, eff = mp.mpf(0), mp.mpf(0), mp.inf
    else:
        delta0 = k2 ** mp.mpf("0.2") * x2k ** mp.mpf("-0.4")
        nscale = delta0 / (mp.mpf("0.5") / mp.sqrt(PI)) ** mp.mpf("0.2")
        eff = k2 * mp.sqrt(x2k)

    def fl(v):
        return math.inf if v == mp.inf else float(v)

    return {"integral": float(integral), "k2": float(k2), "x2k": fl(x2k), "delta0": fl(delta0),
            "nrd-scale": fl(nscale), "efficiency": fl(eff),
            "support-radius": 1.0 if name in COMPACT else math.inf, "radius": float(radius)}


# ------------------------------------------------------------------ datasets
def datasets():
    ds = []

    def add(name, data, h, extra=()):
        ds.append({"name": name, "data": [float(v) for v in data], "h": float(h), "extra": list(extra)})

    add("single", [3.7], 0.8)
    add("anchor-gap10", [0.0, 10.0], 1.0, [9.5, 9.2, 9.0, 8.99, 5.0, 0.5, 1.0, 10.5, 0.25, 9.75])
    add("anchor-gap100", [0.0, 100.0], 1.0, [98.0, 99.0, 96.0, 97.5, 50.0, 2.0, 1.5, 100.5, 95.0, 90.0])
    add("n5-normal", np.random.default_rng(2026101901).normal(0, 1, 5), 0.7)
    add("n30-normal", np.random.default_rng(2026101902).normal(0, 1, 30), 0.45)
    add("n200-normal", np.random.default_rng(2026101903).normal(0, 1, 200), 0.3)
    g = np.random.default_rng(2026101904)
    add("bimodal-wide-gap", np.concatenate([g.normal(-6, 0.5, 30), g.normal(6, 0.5, 30)]), 0.4)
    g = np.random.default_rng(2026101905)
    add("ties200", np.round(g.normal(0, 3, 200) * 2) / 2, 0.26)
    g = np.random.default_rng(2026101906)
    add("cluster-outlier", np.concatenate([g.normal(0, 0.05, 30), [50.0]]), 0.1)
    add("equally-spaced", np.arange(20) * 0.5, 0.7)
    add("heavy-tail-t2", np.random.default_rng(2026101907).standard_t(2, 100), 0.5)
    add("offset-1e6", 1e6 + np.random.default_rng(2026101908).normal(0, 1, 50), 0.3)
    g = np.random.default_rng(2026101909)
    add("two-clusters-100h", np.concatenate([g.normal(0, 0.1, 15), g.normal(25, 0.1, 15)]), 0.25)
    return ds


def eval_points(data, h, R, rng, extra, uniform_filter, total=240, eps=1e-6):
    u = np.unique(np.asarray(data, dtype=float))
    mn, mx = float(u[0]), float(u[-1])
    pts = list(extra)
    for f in [1.0, 1 - eps, 1 + eps, 0.5, 0.25, 0.75, 0.9, 0.99, 1.01, 1.05, 1.5, 2.0]:
        pts += [mn - f * R * h, mx + f * R * h]
    for k in [0.5, 1, 2, 3, 5, 10, 50, 1e3, 1e5]:
        pts += [mn - k * h, mx + k * h]
    if len(u) > 1:
        g = np.diff(u)
        for gi in np.argsort(-g, kind="stable")[:5]:
            a, b = float(u[gi]), float(u[gi + 1])
            cand = [a + t * (b - a) for t in (0.1, 0.25, 0.5, 0.75, 0.9)]
            for f in (0.25, 0.5, 0.9, 1 - eps, 1 + eps, 1.01, 1.5):
                cand += [a + f * R * h, b - f * R * h]
            pts += [c for c in cand if a < c < b]
    sel = sorted({0, len(u) - 1, len(u) // 2, len(u) // 4, (3 * len(u)) // 4})
    for i in sel:
        d = float(u[i])
        for f in (1 - eps, 1 + eps, 0.5, 0.01):
            pts += [d - f * h, d + f * h]
    step = max(1, len(u) // 30)
    pts += [float(v) for v in u[::step]]
    seen, uniq = set(), []
    for p in pts:
        if p not in seen:
            seen.add(p)
            uniq.append(p)
    if uniform_filter:
        uniq = [p for p in uniq if not np.any(np.abs(np.abs((p - u) / h) - 1.0) < 1e-9)]
    uniq = uniq[:total - 8]
    seen = set(uniq)
    span = h * min(R, 4.0) * 1.3
    while len(uniq) < total:
        c = float(u[rng.integers(len(u))])
        p = float(c + rng.uniform(-1, 1) * span)
        if p in seen:
            continue
        if uniform_filter and np.any(np.abs(np.abs((p - u) / h) - 1.0) < 1e-9):
            continue
        seen.add(p)
        uniq.append(p)
    return sorted(uniq)


# ------------------------------------------------------------------ cross-checks
def other_implementations(name, data, h, xs):
    """Returns {impl: density array} for the implementations that define this kernel."""
    out = {"numpy": density_np(name, data, h, xs)}
    data = np.asarray(data, dtype=float)
    xs = np.asarray(xs, dtype=float)
    n = len(data)
    if name == "gaussian" and n > 1 and np.std(data, ddof=1) > 0:
        kde = scipy.stats.gaussian_kde(data, bw_method=h / np.std(data, ddof=1))
        out["scipy.gaussian_kde"] = kde(xs)
    skl = {"gaussian": "gaussian", "uniform": "tophat", "epanechnikov": "epanechnikov",
           "cosine": "cosine", "triangular": "linear", "laplace": "exponential"}
    if name in skl:
        kd = sklearn.neighbors.KernelDensity(kernel=skl[name], bandwidth=h, atol=0.0, rtol=0.0,
                                             algorithm="kd_tree")
        kd.fit(data[:, None])
        with np.errstate(divide="ignore"):
            out["sklearn.KernelDensity"] = np.exp(kd.score_samples(xs[:, None]))
    smk = {"gaussian": "gau", "uniform": "uni", "epanechnikov": "epa", "triangular": "tri",
           "quartic": "biw", "triweight": "triw", "cosine": "cos"}
    if name in smk:
        try:
            kd = sm.nonparametric.KDEUnivariate(data)
            kd.fit(kernel=smk[name], bw=h, fft=False, cut=0)
            out["statsmodels.KDEUnivariate"] = np.array([float(np.squeeze(kd.evaluate(x))) for x in xs])
        except Exception as e:  # pragma: no cover - recorded in the header
            out["statsmodels.KDEUnivariate"] = None
    return out


def main():
    X = [0.0, 1e-9, 1e-6, 0.1, 0.25, 0.5, 0.75, 0.9, 1 - 1e-6, 1 - 1e-9, 1.0, 1 + 1e-9, 1 + 1e-6, 1.5, 2.0,
         3.0, 5.0, 6.65, 10.0, 23.1, 30.0, 50.0, 100.0, 1000.0, 1e6]
    KX = X + [-v for v in X[1:]] + [math.inf, -math.inf]
    kvals = {k: [float(K_mp(k, mp.mpf(v)) if not math.isinf(v) else mp.mpf(0)) for v in KX] for k in KERNELS}
    consts = {k: kernel_constants(k) for k in KERNELS}

    anchors = []
    for kid, kern, data, h, x in [("ke-01-epanechnikov", "epanechnikov", [0.0, 10.0], 1.0, 9.5),
                                  ("ke-01-epanechnikov-9.2", "epanechnikov", [0.0, 10.0], 1.0, 9.2),
                                  ("ke-01-gaussian", "gaussian", [0.0, 100.0], 1.0, 98.0),
                                  ("ke-01-uniform-9.5", "uniform", [0.0, 10.0], 1.0, 9.5)]:
        anchors.append({"id": Kw(kid), "kernel": Kw(kern), "data": data, "h": h, "x": x,
                        "density": density_mp(kern, data, h, [x])[0]})
    assert anchors[0]["density"] == 0.28125
    assert abs(anchors[2]["density"] - 0.026995483256594) < 1e-14

    cases = []
    maxdiff = {}
    dss = datasets()
    for di, ds in enumerate(dss):
        for ki, kern in enumerate(KERNELS):
            rng = np.random.default_rng([2026101910, di, ki])
            R = float(consts[kern]["radius"]) if kern != "cauchy" else 30.0
            xs = eval_points(ds["data"], ds["h"], R, rng, ds["extra"], kern in ("uniform", "wigner"),
                             eps=1e-3 if kern == "wigner" else 1e-6)
            assert len(xs) == 240
            dens = density_mp(kern, ds["data"], ds["h"], xs)
            fmax = max(dens)
            for impl, arr in other_implementations(kern, ds["data"], ds["h"], xs).items():
                if arr is None:
                    continue
                err = float(np.max(np.abs(np.asarray(arr) - np.asarray(dens))) / fmax)
                key = (impl, kern)
                maxdiff[key] = max(maxdiff.get(key, 0.0), err)
            cases.append({"id": Kw(f"{ds['name']}-{kern}"), "dataset": Kw(ds["name"]), "kernel": Kw(kern),
                          "data": ds["data"], "h": ds["h"], "x": xs, "density": dens})

    for (impl, kern), v in maxdiff.items():
        # scipy whitens the data (x / sd): at the offset 1e6 dataset that costs ~1e-10 in u
        lim = 1e-9 if impl.startswith(("statsmodels", "scipy")) else 1e-12
        assert v < lim, (impl, kern, v)

    by_impl = {}
    for (impl, kern), v in sorted(maxdiff.items()):
        by_impl.setdefault(impl, []).append(f"{kern} {v:.1e}")
    header = [
        ";; kde_reference.edn - generated by utils/fastmath/dev/generate_kde_reference.py (do not edit)",
        ";; Oracle: exact direct-sum KDE in mpmath %s (40 digits, no cut-off, no binning), 13 datasets x 14 kernels,"
        " 240 points per case; kernel constants by mpmath quad (30 digits)." % mp.__version__,
        ";; Second implementations (max |other - oracle| / max oracle per kernel over all datasets):",
    ]
    for impl, items in by_impl.items():
        header.append(";;   %s: %s" % (impl, ", ".join(items)))
    header.append(";; numpy %s, scipy %s, scikit-learn %s, statsmodels %s. See the generator docstring for the schema."
                  % (np.__version__, scipy.__version__, sklearn.__version__, statsmodels.__version__))
    body = {
        Kw("kernels"): [Kw(k) for k in KERNELS],
        Kw("anchors"): [{Kw(k): v for k, v in a.items()} for a in anchors],
        Kw("kernel-values"): {Kw("x"): KX, Kw("values"): {Kw(k): v for k, v in kvals.items()}},
        Kw("kernel-constants"): {Kw(k): {Kw(a): b for a, b in c.items()} for k, c in consts.items()},
        Kw("cases"): [{Kw(k): v for k, v in c.items()} for c in cases],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(header) + "\n" + emit(body) + "\n")
    print("wrote", OUT, OUT.stat().st_size, "bytes;", len(cases), "cases")
    for impl, items in by_impl.items():
        print(impl, "max", max(float(i.split()[-1]) for i in items))


if __name__ == "__main__":
    main()
