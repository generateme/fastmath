"""Generator of the Python part of test/resources/kernel/bandwidth_reference.edn  (step 2 of 2).

Regeneration (byte-identical; run in this order, from anywhere):
    Rscript /home/ts/clojure/fastmath/utils/fastmath/dev/generate_bandwidth_reference.R
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_bandwidth_reference.py
Step 1 writes the 32 datasets and the R outputs; this script reads that file (mini EDN reader below), drops
the keys it owns (`:py` in every case, `:cv-info`, `:cv-gaussian-battery`, and the header lines marked
`[py]`) and writes the whole file again (idempotent: a second run of this script alone gives the same bytes).
`generate_rlcv_reference.py` imports the criteria of this file.

Tools: numpy 2.5.2, scipy 1.18.1 (only `minimize_scalar` for the refinement); no Clojure code is used;
kernels are written from textbook formulas.  Double precision.

1. Bandwidth rules (per case, key `:py`).  Type 7 quantiles (`numpy.percentile`, linear), sample sd (n - 1):
       :nrd-chain    1.06 lo n^-0.2   lo = min(sd, IQR/1.34), else sd, else |x[0]|, else 1   (the R bw.nrd0
                     fallback chain with the 1.06 constant: the proposal D6/F2 for the Clojure `:nrd`)
       :nrd0-chain   0.9 lo n^-0.2   (equals R bw.nrd0 to rounding; independent check)
       :nrd-plain    1.06 min(sd, IQR/1.34) n^-0.2  without fallback (equals R bw.nrd, 0 when IQR = 0)
       :sj-ste :sj-dpi :sj-ste-tight   port of R `bw.SJ` (R 4.5.2 sources: `bw_pair_cnts`, C `bw_den`,
                     `bw_den_binned`, `bw_phi4`, `bw_phi6` with DELMAX 1000, `uniroot` = `R_zeroin2` ported line
                     by line, `tol = 0.1 * lower` evaluated lazily AFTER the interval widening, as R does)
   with the same outcome schema as the R part ({:value v} or {:error true :message s}).  `n < 2` gives
   "need at least 2 data points" for all of them.  The port is compared with R on all 32 datasets
   (agreement is written in the header of the produced file).
   Observed on the committed data: max relative difference to R 1.7e-16 (nrd0), 2.2e-16 (nrd), 4.3e-15
   (sj-ste), 4.9e-15 (sj-dpi), 2.7e-15 (sj-ste-tight), all error messages equal.

2. Cross-validation criteria (the independent oracle used by `generate_rlcv_reference.py` too).  z =
   (x - mean) / sd (sample sd), h in standardised units, K the kernel, n points:
       f(t)   = 1/(n h) sum_j K((t - z_j)/h)
       f_i    = 1/((n-1) h) sum_{j != i} K((z_i - z_j)/h)          (leave one out; no floor)
       LSCV   = int f^2 dt - 2 mean(f_i)                           (minimised)
       LCV    = mean(log f_i)                                      (maximised; NaN if some f_i < 0, -Inf if 0)
       RLCV   = mean(g_a(f_i)) - b*(h)                             (maximised)
       a = (1/sqrt 2) sqrt(ln n) / n ;  g_a(x) = ln x (x >= a), ln a - 1 + x/a (x < a)
       b*  printed:     int B(f) dt,  B(v) = v (v >= a), v^2/(2a) (v < a)     (Wu 2019 as printed, KB-12)
       b*  continuous:  B(v) = v - a/2 (v >= a), v^2/(2a) (v < a)
   Integrals: trapezoid on a uniform grid over [min z - e h, max z + e h], e = the extent (radius where
   |K| <= 1e-10: 6.65 gaussian, 23.1 logistic, 22.6 sigmoid, 29.7 silverman, 22.4 laplace, 1 compact kernels,
   100 cauchy), spacing h/50 (or the stated divisor).  f on the grid is the exact sum; for the infinite
   kernels except cauchy the sum over data farther than the window 10 (gaussian), 47, 47, 65, 46 h (logistic,
   sigmoid, silverman, laplace) is skipped, where |K| < 1e-20; cauchy sums every datum.
   Cauchy tail: its extent 100 cuts a mass of 2 * (1/pi) * (pi/2 - atan(200)) = 3.18e-3 of every kernel, so
   the cauchy integrals are lower than the integrals over the whole line (the tail error of the criterion
   is of that order; with the specified extent 1000 it would be 3.18e-4).
   The scale of the search: s = nrd(z) * nrd-scale with nrd(z) = 1.06 lo n^-0.2 (chain above) and nrd-scale =
   delta0 * (0.5/sqrt(pi))^-0.2, delta0 = k2^0.2 x2k^-0.4 (s = nrd(z) for silverman, cauchy); the h grid is
   4000 points log-uniform over [0.01, 5] s, refined by bounded Brent in log h (xatol 1e-9) in the bracket
   of the best grid point; `:at-bound` when the optimum is within 1e-3 (relative) of a grid end.

3. Gaussian battery `:cv-gaussian-battery`: 32 datasets (n 20 ... 500; normal, scaled sd 0.01 / 10 / 1000, offset
   1e6, lognormal, bimodal, trimodal, t(3), t(2), exponential, ties, uniform), seeded numpy
   `default_rng`.  Per dataset the optimum h (ORIGINAL units, = standardised optimum * sd) of :lscv :lcv
   :rlcv-printed :rlcv-continuous with the value (of the standardised criterion), `:h-std`, `:at-bound`.
   The closed form of the gaussian LSCV integral (K*K = N(0, 2)) is compared with the grid integral at
   h = s (value in the header).

Why: KB-01/02 (units and range), KB-05/06 (bound and local optimum), KB-07/08 (nrd), KB-11 (floor), KB-12
(RLCV), D4/D5/D6/D7 of the Kernel Density audit; work packages 2.1, 2.3, 2.5, 2.6.
"""
import math
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "test" / "resources" / "kernel" / "bandwidth_reference.edn"

KERNELS = ["uniform", "triangular", "epanechnikov", "quartic", "triweight", "tricube", "gaussian",
           "cosine", "logistic", "sigmoid", "silverman", "wigner", "cauchy", "laplace"]
PI = math.pi

# (k2, x2k) closed forms (verified against mpmath in the audit and in kde_reference.edn)
CONSTS = {
    "uniform": (0.5, 1 / 3), "triangular": (2 / 3, 1 / 6), "epanechnikov": (0.6, 0.2),
    "quartic": (5 / 7, 1 / 7), "triweight": (350 / 429, 1 / 9), "tricube": (175 / 247, 35 / 243),
    "gaussian": (0.5 / math.sqrt(PI), 1.0), "cosine": (PI ** 2 / 16, 1 - 8 / PI ** 2),
    "logistic": (1 / 6, PI ** 2 / 3), "sigmoid": (2 / PI ** 2, PI ** 2 / 4),
    "silverman": (3 * math.sqrt(2) / 16, 0.0), "wigner": (16 / (3 * PI ** 2), 0.25),
    "cauchy": (1 / PI, math.inf), "laplace": (0.25, 2.0),
}
EXT = {"uniform": 1.0, "triangular": 1.0, "epanechnikov": 1.0, "quartic": 1.0, "triweight": 1.0,
       "tricube": 1.0, "cosine": 1.0, "wigner": 1.0, "gaussian": 6.65, "logistic": 23.1, "sigmoid": 22.6,
       "silverman": 29.7, "laplace": 22.4, "cauchy": 100.0}
# window radius of the grid sums (|K| < 1e-20 beyond); None = every datum
WIN = {"uniform": 1.0, "triangular": 1.0, "epanechnikov": 1.0, "quartic": 1.0, "triweight": 1.0,
       "tricube": 1.0, "cosine": 1.0, "wigner": 1.0, "gaussian": 10.0, "logistic": 47.0, "sigmoid": 47.0,
       "silverman": 65.0, "laplace": 46.0, "cauchy": None}


# ------------------------------------------------------------------ EDN
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
        return '"' + v.replace("\\", "\\\\").replace('"', '\\"') + '"'
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
        return "{\n" + "\n".join(ind + fmt(k) + " " + emit(e, depth + 1) for k, e in v.items()) + "}"
    if isinstance(v, list) and v and all(isinstance(e, dict) for e in v):
        ind = "  " * (depth + 1)
        return "[\n" + "\n".join(ind + fmt(e) for e in v) + "]"
    return fmt(v)


_TOK = re.compile(r'\s+|,|;[^\n]*|(\{|\}|\[|\]|"(?:[^"\\]|\\.)*"|##-?(?:Inf|NaN)|[^\s\[\]{}",;]+)')


def read_edn(text):
    toks = [m.group(1) for m in _TOK.finditer(text) if m.group(1) is not None]
    pos = [0]

    def val():
        t = toks[pos[0]]
        pos[0] += 1
        if t == "{":
            d = {}
            while toks[pos[0]] != "}":
                k = val()
                d[k] = val()
            pos[0] += 1
            return d
        if t == "[":
            a = []
            while toks[pos[0]] != "]":
                a.append(val())
            pos[0] += 1
            return a
        if t.startswith('"'):
            return re.sub(r'\\(.)', r'\1', t[1:-1])
        if t == "##Inf":
            return math.inf
        if t == "##-Inf":
            return -math.inf
        if t == "##NaN":
            return math.nan
        if t.startswith(":"):
            return Kw(t[1:])
        if t == "true":
            return True
        if t == "false":
            return False
        if t == "nil":
            return None
        if re.fullmatch(r"[-+]?\d+", t):
            return int(t)
        return float(t)

    v = val()
    assert pos[0] == len(toks)
    return v


# ------------------------------------------------------------------ nrd family
def nrd_lo(x):
    """lo of R bw.nrd0: min(sd, IQR/1.34), else sd, else |x[0]|, else 1 (x in the original order)."""
    x = np.asarray(x, dtype=float)
    sd = float(np.std(x, ddof=1))
    q = np.percentile(x, [25, 75])
    lo = min(sd, float((q[1] - q[0]) / 1.34))
    if not lo:
        lo = sd or abs(float(x[0])) or 1.0
    return lo


def nrd_plain(x, c):
    x = np.asarray(x, dtype=float)
    q = np.percentile(x, [25, 75])
    return c * min(float(np.std(x, ddof=1)), float((q[1] - q[0]) / 1.34)) * len(x) ** -0.2


def nrd_chain(x, c):
    return c * nrd_lo(x) * len(x) ** -0.2


# ------------------------------------------------------------------ port of R bw.SJ
class RError(Exception):
    pass


def zeroin2(f, ax, bx, fa, fb, tol, maxit=1000):
    """R_zeroin2 (src/library/stats/src/zeroin.c)."""
    eps = np.finfo(float).eps
    a, b = ax, bx
    c, fc = a, fa
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    for _ in range(maxit + 1):
        prev_step = b - a
        if abs(fc) < abs(fb):
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
        tol_act = 2 * eps * abs(b) + tol / 2
        new_step = (c - b) / 2
        if abs(new_step) <= tol_act or fb == 0.0:
            return b
        if abs(prev_step) >= tol_act and abs(fa) > abs(fb):
            cb = c - b
            if a == c:
                t1 = fb / fa
                p = cb * t1
                q = 1.0 - t1
            else:
                q = fa / fc
                t1 = fb / fc
                t2 = fb / fa
                p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0))
                q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0)
            if p > 0:
                q = -q
            else:
                p = -p
            if p < (0.75 * cb * q - abs(tol_act * q) / 2) and p < abs(prev_step * q / 2):
                new_step = p / q
        if abs(new_step) < tol_act:
            new_step = tol_act if new_step > 0 else -tol_act
        a, fa = b, fb
        b += new_step
        fb = f(b)
        if (fb > 0 and fc > 0) or (fb < 0 and fc < 0):
            c, fc = a, fa
    raise RError("zeroin did not converge")


def bw_pair_cnts(x, nb, binned):
    """stats:::bw_pair_cnts: returns (d, cnt[0..nb-1])."""
    n = len(x)
    if binned:
        d = (x.max() - x.min()) * 1.01 / nb
        if not d > 0:
            raise RError("data are constant in bandwidth calculation")
        xx = np.trunc(np.abs(x) / d) * np.sign(x)
        xx = xx - xx.min() + 1
        w = np.bincount(xx.astype(np.int64), minlength=nb + 1)[1:nb + 1].astype(float)
        conv = np.convolve(w, w[::-1])
        cnt = np.zeros(nb)
        cnt[1:] = conv[nb:2 * nb - 1]
        cnt[0] = 0.5 * float(np.sum(w * (w - 1.0)))
        return d, cnt
    rang = (x.max() - x.min()) * 1.01
    dd = rang / nb
    if not dd > 0:
        raise RError("data are constant in bandwidth calculation")
    ii = np.trunc(x / dd).astype(np.int64)
    iu = np.triu_indices(n, 1)
    cnt = np.bincount(np.abs(ii[iu[0]] - ii[iu[1]]), minlength=nb).astype(float)
    return dd, cnt


def pw(x, p):
    """R's x^p: NaN for negative x (Python would return a complex number)."""
    return float(np.power(np.float64(x), p))


def div(a, b):
    with np.errstate(all="ignore"):
        return float(np.float64(a) / np.float64(b))


def _phi(cnt, n, d, h, order):
    with np.errstate(all="ignore"):
        return _phi0(cnt, n, d, h, order)


def _phi0(cnt, n, d, h, order):
    i = np.arange(len(cnt))
    delta = (i * d) / h
    delta = delta * delta
    keep = delta < 1000.0
    dl = delta[keep]
    if order == 4:
        term = np.exp(-dl / 2) * (dl * dl - 6 * dl + 3)
        s = 2 * float(np.sum(term * cnt[keep])) + n * 3
        return float(np.float64(s) / np.float64(n * (n - 1) * h ** 5.0 * math.sqrt(2 * PI)))
    term = np.exp(-dl / 2) * (dl * dl * dl - 15 * dl * dl + 45 * dl - 15)
    s = 2 * float(np.sum(term * cnt[keep])) - 15 * n
    return float(np.float64(s) / np.float64(n * (n - 1) * h ** 7.0 * math.sqrt(2 * PI)))


def bw_sj(x, method="ste", tol=None, nb=1000):
    """Port of stats::bw.SJ with the default lower/upper (tol = None means 0.1 * lower after widening)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 2:
        raise RError("need at least 2 data points")
    d, cnt = bw_pair_cnts(x, nb, n > nb / 2)
    q = np.percentile(x, [25, 75])
    scale = min(float(np.std(x, ddof=1)), float(q[1] - q[0]) / 1.349)
    a = 1.24 * scale * n ** (-1 / 7)
    b = 1.23 * scale * n ** (-1 / 9)
    c1 = 1 / (2 * math.sqrt(PI) * n)
    TD = -_phi(cnt, n, d, b, 6)
    if not math.isfinite(TD) or TD <= 0:
        raise RError("sample is too sparse to find TD")
    if method == "dpi":
        return pw(div(c1, _phi(cnt, n, d, pw(2.394 / (n * TD), 1 / 7), 4)), 1 / 5)
    hmax = 1.144 * scale * n ** (-1 / 5)
    lower, upper = 0.1 * hmax, hmax
    alph2 = 1.357 * pw(div(_phi(cnt, n, d, a, 4), TD), 1 / 7)
    if not math.isfinite(alph2):
        raise RError("sample is too sparse to find alph2")

    def fSD(h):
        s = _phi(cnt, n, d, alph2 * h ** (5 / 7), 4)
        return pw(div(c1, s), 1 / 5) - h

    itry = 1
    while True:
        p = fSD(lower) * fSD(upper)
        if math.isnan(p):
            raise RError("missing value where TRUE/FALSE needed")
        if not p > 0:
            break
        if itry > 99:
            raise RError("no solution in the specified range of bandwidths")
        if itry % 2:
            upper = upper * 1.2
        else:
            lower = lower / 1.2
        itry += 1
    t = 0.1 * lower if tol is None else tol
    return zeroin2(fSD, lower, upper, fSD(lower), fSD(upper), t)


def outcome(fn):
    try:
        v = fn()
    except RError as e:
        return {Kw("error"): True, Kw("message"): str(e)}
    if not math.isfinite(v):
        return {Kw("error"): True, Kw("message"): "NA"}
    return {Kw("value"): float(v)}


# ------------------------------------------------------------------ kernels (textbook formulas)
def K_np(name, u):
    a = np.abs(u)
    inb = a <= 1
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
        return np.exp(-0.5 * u * u) / math.sqrt(2 * PI)
    if name == "cosine":
        return np.where(inb, PI / 4 * np.cos(PI / 2 * u), 0.0)
    with np.errstate(over="ignore"):
        if name == "logistic":
            return 1 / (np.exp(a) + 2 + np.exp(-a))
        if name == "sigmoid":
            return (2 / PI) / (np.exp(a) + np.exp(-a))
    if name == "silverman":
        t = a / math.sqrt(2)
        return 0.5 * np.exp(-t) * np.sin(t + PI / 4)
    if name == "wigner":
        return np.where(inb, 2 / PI * np.sqrt(np.clip(1 - u * u, 0, None)), 0.0)
    if name == "cauchy":
        return 2 / (PI * (1 + 4 * u * u))
    if name == "laplace":
        return 0.5 * np.exp(-a)
    raise KeyError(name)


def nrd_scale(name):
    k2, x2k = CONSTS[name]
    if x2k == 0.0 or math.isinf(x2k):
        return None
    return k2 ** 0.2 * x2k ** -0.4 * (1 / ((0.5 / math.sqrt(PI)) ** 0.2))


def search_scale(z, name):
    """s = nrd(z) * nrd-scale (1.06 chain), nrd(z) alone where the scale is undefined."""
    base = nrd_chain(z, 1.06)
    sc = nrd_scale(name)
    return base if sc is None else base * sc


# ------------------------------------------------------------------ criteria
def rlcv_a(n):
    return math.sqrt(math.log(n)) / (math.sqrt(2) * n)


def loo_density(z, h, name):
    n = len(z)
    M = K_np(name, (z[:, None] - z[None, :]) / h)
    np.fill_diagonal(M, 0.0)
    return M.sum(1) / ((n - 1) * h)


def grid_f(z, h, name, div):
    """Exact kernel sum on the uniform grid [min - e h, max + e h] with spacing <= h/div."""
    e = EXT[name]
    n = len(z)
    lo = float(z.min()) - e * h
    hi = float(z.max()) + e * h
    N = int(math.ceil((hi - lo) * div / h))
    delta = (hi - lo) / N
    wr = WIN[name]
    if wr is None:
        g = lo + delta * np.arange(N + 1)
        f = np.empty(N + 1)
        blk = max(1, int(2e6 // n))
        for s in range(0, N + 1, blk):
            f[s:s + blk] = K_np(name, (g[s:s + blk, None] - z[None, :]) / h).sum(1)
    else:
        half = int(wr * div) + 2
        off = np.arange(-half, half + 1)
        c = np.floor((z - lo) / delta).astype(np.int64)
        f = np.zeros(N + 1)
        blk = max(1, int(2e6 // len(off)))
        for s in range(0, n, blk):
            idx = c[s:s + blk, None] + off[None, :]
            ok = (idx >= 0) & (idx <= N)
            idc = np.where(ok, idx, 0)
            t = (lo + idc * delta - z[s:s + blk, None]) / h
            kv = np.where(ok, K_np(name, t), 0.0)
            f += np.bincount(idc.ravel(), weights=kv.ravel(), minlength=N + 1)
    f /= n * h
    return delta, f


def trap(delta, y):
    return float(delta * (y.sum() - 0.5 * (y[0] + y[-1])))


def criteria(z, h, name, div=50, need_grid=True):
    """All components and the four criteria (as defined in the module docstring)."""
    n = len(z)
    a = rlcv_a(n)
    loo = loo_density(z, h, name)
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_log = float(np.mean(np.log(loo)))
    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.where(loo >= a, np.log(np.where(loo >= a, loo, 1.0)), math.log(a) - 1 + loo / a)
    mean_g = float(np.mean(g))
    mean_loo = float(np.mean(loo))
    out = {"mean-loo": mean_loo, "mean-log-loo": mean_log, "mean-g": mean_g, "lcv": mean_log}
    if need_grid:
        delta, f = grid_f(z, h, name, div)
        i2 = trap(delta, f * f)
        low = f * f / (2 * a)
        bp = trap(delta, np.where(f >= a, f, low))
        bc = trap(delta, np.where(f >= a, f - a / 2, low))
        out.update({"int-f2": i2, "b-printed": bp, "b-continuous": bc,
                    "measure-above-a": trap(delta, (f >= a).astype(float)), "lscv": i2 - 2 * mean_loo,
                    "rlcv-printed": mean_g - bp, "rlcv-continuous": mean_g - bc})
    return out


CRIT = [("lscv", -1.0), ("lcv", 1.0), ("rlcv-printed", 1.0), ("rlcv-continuous", 1.0)]  # (name, sign: +1 maximise)


def scan_optimum(z, name, ngrid=4000, div=50, lo_f=0.01, hi_f=5.0):
    """Log-grid scan + bounded Brent refinement for the four criteria. Returns dict criterion -> result."""
    s = search_scale(z, name)
    hs = np.exp(np.linspace(math.log(lo_f * s), math.log(hi_f * s), ngrid))
    vals = {c: np.empty(ngrid) for c, _ in CRIT}
    for i, h in enumerate(hs):
        r = criteria(z, float(h), name, div)
        for c, _ in CRIT:
            vals[c][i] = r[c]
    res = {}
    for c, sign in CRIT:
        v = vals[c] * sign
        v = np.where(np.isnan(v), -np.inf, v)
        i = int(np.argmax(v))
        j0, j1 = max(i - 1, 0), min(i + 1, ngrid - 1)

        def obj(lh, c=c, sign=sign):
            r = criteria(z, math.exp(lh), name, div)[c] * sign
            return -r if math.isfinite(r) else math.inf

        best_h, best_v = float(hs[i]), float(vals[c][i])
        if math.isfinite(v[i]):
            m = minimize_scalar(obj, bounds=(math.log(hs[j0]), math.log(hs[j1])), method="bounded",
                                options={"xatol": 1e-9})
            if -m.fun >= v[i]:
                best_h, best_v = math.exp(m.x), float(-m.fun * sign)
        at_bound = bool(abs(best_h / hs[0] - 1) < 1e-3 or abs(best_h / hs[-1] - 1) < 1e-3)
        res[c] = {"h": best_h, "value": best_v, "at-bound": at_bound, "grid-index": i}
    return res, float(s), float(hs[0]), float(hs[-1])


def job_scan(args):
    z, name, ngrid, div = args
    return scan_optimum(np.asarray(z), name, ngrid, div)


def standardise(x):
    x = np.asarray(x, dtype=float)
    return (x - x.mean()) / x.std(ddof=1)


# ------------------------------------------------------------------ gaussian closed form (check)
def int_f2_gaussian_exact(z, h):
    n = len(z)
    d = (z[:, None] - z[None, :]) / h
    return float(np.exp(-d * d / 4).sum() / (2 * math.sqrt(PI)) / (n * n * h))


# ------------------------------------------------------------------ battery datasets
def battery():
    out = []

    def add(name, x):
        out.append((name, np.asarray(x, dtype=float)))

    seed = [2026102000]

    def rng():
        seed[0] += 1
        return np.random.default_rng(seed[0])

    for n in (20, 50, 100, 200, 500):
        add(f"normal-n{n}", rng().normal(0, 1, n))
    for sd in (0.01, 10.0, 1000.0):
        add(f"normal-sd{sd:g}-n100", sd * rng().normal(0, 1, 100))
    for n in (30, 100, 500):
        add(f"normal-offset1e6-n{n}", 1e6 + rng().normal(0, 1, n))
    add("normal-offset1e6-sd1000-n200", 1e6 + 1000 * rng().normal(0, 1, 200))
    for n in (30, 100, 300, 500):
        add(f"lognormal-n{n}", rng().lognormal(0, 1, n))
    for n in (50, 100, 200, 500):
        g = rng()
        c = g.random(n) < 0.5
        add(f"bimodal-n{n}", np.where(c, g.normal(-3, 0.5, n), g.normal(3, 0.5, n)))
    for n in (300, 500):
        g = rng()
        u = g.random(n)
        add(f"trimix-n{n}", np.where(u < 0.5, g.normal(0, 1, n), np.where(u < 0.8, g.normal(5, 0.1, n),
                                                                           g.laplace(-4, 0.3, n))))
    for n in (100, 500):
        add(f"t3-n{n}", rng().standard_t(3, n))
    add("t2-n200", rng().standard_t(2, 200))
    for n in (100, 300):
        add(f"exponential-n{n}", rng().exponential(1.0, n))
    add("ties-integers-n100", np.round(rng().normal(0, 2, 100)))
    add("ties-integers-n200", np.round(rng().normal(0, 3, 200)))
    add("ties-halves-n200", np.round(rng().normal(0, 3, 200) * 2) / 2)
    add("uniform-n100", rng().uniform(0, 1, 100))
    add("normal-negative-offset-sd0.1-n60", -50 + 0.1 * rng().normal(0, 1, 60))
    assert len(out) == 32, len(out)
    return out


def main():
    text = OUT.read_text()
    data = read_edn(text)
    cases = data[Kw("cases")] if isinstance(data, dict) else None
    assert cases, "run generate_bandwidth_reference.R first"
    data = {k: v for k, v in data.items()}
    hdr = [l for l in text.splitlines() if l.startswith(";;") and "[py]" not in l]

    # ---- 1. nrd family and SJ port
    rel = {"nrd0": 0.0, "nrd": 0.0, "ste": 0.0, "dpi": 0.0, "ste-tight": 0.0}
    status_mismatch = []
    new_cases = []
    for c in cases:
        x = np.asarray(c[Kw("data")], dtype=float)
        n = len(x)

        def guard(fn):
            if n < 2:
                return {Kw("error"): True, Kw("message"): "need at least 2 data points"}
            return {Kw("value"): float(fn())}

        py = {
            Kw("nrd-chain"): guard(lambda: nrd_chain(x, 1.06)),
            Kw("nrd0-chain"): guard(lambda: nrd_chain(x, 0.9)),
            Kw("nrd-plain"): guard(lambda: nrd_plain(x, 1.06)),
            Kw("sj-ste"): outcome(lambda: bw_sj(x, "ste")),
            Kw("sj-dpi"): outcome(lambda: bw_sj(x, "dpi")),
            Kw("sj-ste-tight"): outcome(lambda: bw_sj(x, "ste", tol=_tight(x))),
        }
        for rk, pk, tag in [("nrd0", "nrd0-chain", "nrd0"), ("nrd", "nrd-plain", "nrd"), ("sj-ste", "sj-ste", "ste"),
                            ("sj-dpi", "sj-dpi", "dpi"), ("sj-ste-tight", "sj-ste-tight", "ste-tight")]:
            r, p = c[Kw(rk)], py[Kw(pk)]
            re_, pe = r.get(Kw("error"), False), p.get(Kw("error"), False)
            if re_ or pe:
                if re_ != pe or (rk.startswith("sj") and r[Kw("message")] != p[Kw("message")]):
                    status_mismatch.append((c[Kw("name")], rk, r, p))
                continue
            rv, pv = r[Kw("value")], p[Kw("value")]
            rel[tag] = max(rel[tag], abs(pv - rv) / max(abs(rv), 1e-300))
        nc = {k: v for k, v in c.items() if k != Kw("data")}
        nc[Kw("py")] = py
        nc[Kw("data")] = c[Kw("data")]
        new_cases.append(nc)
    assert not status_mismatch, status_mismatch
    for tag, v in rel.items():
        assert v < 1e-9, (tag, v)

    # ---- 2. gaussian battery
    jobs, ds = [], battery()
    for name, x in ds:
        jobs.append((list(standardise(x)), "gaussian", 4000, 50))
    with ProcessPoolExecutor(max_workers=14) as ex:
        results = list(ex.map(job_scan, jobs, chunksize=1))
    cv = []
    max_exact_rel = 0.0
    for (name, x), (res, s, hlo, hhi) in zip(ds, results):
        sd = float(np.std(x, ddof=1))
        z = standardise(x)
        ex_ = int_f2_gaussian_exact(z, s)
        gr = criteria(z, s, "gaussian")["int-f2"]
        max_exact_rel = max(max_exact_rel, abs(gr - ex_) / ex_)
        entry = {Kw("name"): Kw(name), Kw("n"): len(x), Kw("mean"): float(x.mean()), Kw("sd"): sd,
                 Kw("s-std"): s, Kw("h-std-range"): [hlo, hhi], Kw("data"): [float(v) for v in x]}
        for c, _ in CRIT:
            r = res[c]
            entry[Kw(c)] = {Kw("h"): r["h"] * sd, Kw("h-std"): r["h"], Kw("value"): r["value"],
                            Kw("at-bound"): r["at-bound"], Kw("grid-index"): r["grid-index"]}
        # sort data to the end of the entry
        d = entry.pop(Kw("data"))
        entry[Kw("data")] = d
        cv.append(entry)

    # ---- write
    nbound = {c: sum(1 for e in cv if e[Kw(c)][Kw("at-bound")]) for c, _ in CRIT}
    hdr += [
        ";; [py] numpy %s, scipy %s. Python part: :py in every case (nrd with the bw.nrd0 fallback chain, bw.SJ port), :cv-info, :cv-gaussian-battery."
        % (np.__version__, scipy.__version__),
        ";; [py] Port of R bw.SJ vs R on the 32 datasets (max relative difference where both give a value; error status and message equal everywhere):"
        " nrd0 %.1e, nrd %.1e, sj-ste %.1e, sj-dpi %.1e, sj-ste-tight %.1e."
        % (rel["nrd0"], rel["nrd"], rel["ste"], rel["dpi"], rel["ste-tight"]),
        ";; [py] Gaussian battery: grid integral of f^2 (spacing h/50) vs the closed form at h = s: max relative difference %.1e;"
        " optima on a bound of the search interval: lscv %d, lcv %d, rlcv-printed %d, rlcv-continuous %d of %d."
        % (max_exact_rel, nbound["lscv"], nbound["lcv"], nbound["rlcv-printed"], nbound["rlcv-continuous"], len(cv)),
    ]
    info = {
        Kw("description"): "Gaussian battery: optimum h (original units = standardised optimum * sd) of the four CV criteria "
                           "on a 4000 point log grid over [0.01, 5] s refined by bounded Brent; s = 1.06 lo n^-0.2 on standardised data; "
                           "criteria as in the generator docstring; at-bound = within 1e-3 of a grid end.",
        Kw("grid-points"): 4000, Kw("grid-factors"): [0.01, 5.0], Kw("integral-divisor"): 50,
        Kw("rlcv-a"): "sqrt(ln n) / (sqrt 2 n)",
    }
    body = {Kw("cases"): new_cases, Kw("cv-info"): info, Kw("cv-gaussian-battery"): cv}
    OUT.write_text("\n".join(hdr) + "\n" + emit(body) + "\n")
    print("wrote", OUT, OUT.stat().st_size)
    print("rel diff vs R:", rel)
    print("battery int-f2 grid vs exact:", max_exact_rel, "at-bound:", nbound)


def _tight(x):
    if len(x) < 2:
        return math.nan
    q = np.percentile(x, [25, 75])
    sc = min(float(np.std(x, ddof=1)), float(q[1] - q[0]) / 1.349)
    hm = 1.144 * sc * len(x) ** -0.2
    return 1e-10 * hm if hm > 0 else math.nan


if __name__ == "__main__":
    main()
