"""Generator of reference values for `fastmath.kernel.window-test` (kernel audit, window fixes).

Produces `test/resources/kernel/window_reference.edn`.

Tools: Python 3, mpmath 1.3 (60 digits, definitions of the windows written independently of the
fastmath code from the published formulas: Doerry, Catalog of Window Taper Functions for Sidelobe
Control, SAND2017-4042; Harris 1978; Singla and Singh; Heinzel et al. 2002), scipy 1.18
(`scipy.signal.windows`).

Why: the audited windows had wrong raw scales (b-spline, tukey, webster-hamming), unstable constants
(singla-singh, b-spline order > 25, quadrature constants of avci-nacaroglu / parzen-exponential /
parzen-geometric), zero ends (knab), clipped negative lobes (blackman-harris-family), a non-standard
vorbis, N = 1 NaN (mottaghi-kashtiban-shayesteh).  Each `:cases` entry is
`{:id :group :window :opts :N :normalize? :expected}` where `:normalize?` is the value passed to the
window (`false` raw scale: continuous windows have integral 1, discrete ones mean 1; `true` the
default normalisation: division by f(0) for continuous model windows, by the maximum for discrete
ones).  `:interior true` marks cases where the two end values are not compared (singular ends).
`:periodic true` marks cases where `:expected` is a periodic window of `:N` samples (the symmetric window of
N + 1 samples without the last one, which is what `fastmath.kernel/window` returns for `:symmetric? false`).

Window additions (WP 5.1 - 5.4): `cosine` (scipy `cosine`), `kaiser` (scipy `kaiser`, mpmath raw scale),
`general-gaussian` (scipy `general_gaussian`, mpmath raw scale), the twelve flat-top sets (mpmath cosine sums and
scipy `general_cosine`).

Run from anywhere (output is byte-identical on every run):
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_window_reference.py
"""
import math
import warnings
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.signal import windows as W

warnings.filterwarnings("ignore")  # scipy chebwin warns for attenuation < 45 dB
mp.mp.dps = 60

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "test" / "resources" / "kernel" / "window_reference.edn"

PI = mp.pi
CASES = []


def fl(x):
    return repr(float(x))


def vec(xs):
    return "[" + " ".join(fl(x) for x in xs) + "]"


def opt_edn(v):
    if isinstance(v, list):
        return "[" + " ".join(fl(c) for c in v) + "]"
    if isinstance(v, int) and not isinstance(v, bool):
        return str(v)
    return fl(v)


def opts_edn(opts):
    return "{" + " ".join(":%s %s" % (k, opt_edn(v)) for k, v in opts.items()) + "}"


def emit(group, window, opts, N, normalize, vals, interior=False, tag="", periodic=False):
    """`periodic`: the expected values are the periodic window of N samples, which fastmath builds as the
    symmetric window of N + 1 samples without the last one (the test does this, `:periodic true`)."""
    ident = "%s %s N=%d %s%s" % (window, opts_edn(opts), N, "raw" if normalize is False else "peak",
                                 " periodic" if periodic else "")
    ident = ident.replace('"', "'")
    CASES.append("{:id \"%s\" :group :%s :window :%s :opts %s :N %d :normalize? %s%s%s :expected %s}"
                 % (ident, group, window, opts_edn(opts), N, "false" if normalize is False else "true",
                    " :interior true" if interior else "", " :periodic true" if periodic else "", vec(vals)))


def grid(N):
    if N == 1:
        return [mp.mpf(0)]
    return [mp.mpf(-1) / 2 + mp.mpf(i) / (N - 1) for i in range(N)]


def integral(shape, breaks=()):
    """Integral of the even `shape` over [-1/2, 1/2]."""
    pts = [mp.mpf(0)] + sorted(set(mp.mpf(b) for b in breaks if 0 < b < 0.5)) + [mp.mpf(1) / 2]
    return 2 * mp.quad(shape, pts)


def continuous_cases(group, window, opts, shape, Ns, breaks=(), clip=False, interior=False, tag="", raw_only=False,
                     scale=None):
    """Cases `:normalize? false` (shape / integral) and `true` (shape / shape(0))."""
    z = scale if scale is not None else integral(shape, breaks)
    f0 = shape(mp.mpf(0))
    for N in Ns:
        xs = grid(N)
        v = [shape(x) for x in xs]
        if clip:
            v = [max(mp.mpf(0), a) for a in v]
        emit(group, window, opts, N, False, [a / z for a in v], interior, tag)
        if not raw_only:
            emit(group, window, opts, N, True, [a / f0 for a in v], interior, tag)


# ----------------------------------------------------------------------------------------------
# b-spline: order * B_n(order * x), B_n the centred cardinal B-spline (truncated power sum)
# ----------------------------------------------------------------------------------------------
def bspline_raw(n):
    def f(x):
        u = n * x
        with mp.workdps(120):
            s = mp.mpf(0)
            for k in range(n + 1):
                d = u + mp.mpf(n) / 2 - k
                if d > 0:
                    s += (-1) ** k * mp.binomial(n, k) * d ** (n - 1)
            return n * s / mp.factorial(n - 1)
    return f


for order in (2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 45, 60):
    f = bspline_raw(order)
    f0 = f(mp.mpf(0))
    for N in (9, 16):
        v = [f(x) if abs(x) < 0.5 else mp.mpf(0) for x in grid(N)]
        emit("b-spline", "b-spline", {"order": order}, N, False, v)
        emit("b-spline", "b-spline", {"order": order}, N, True, [a / f0 for a in v])


# ----------------------------------------------------------------------------------------------
# tukey: plateau + cosine tapers, alpha > 1 is hann
# ----------------------------------------------------------------------------------------------
def tukey_shape(a):
    a = mp.mpf(min(a, 1.0))
    if a == 0:
        return (lambda x: mp.mpf(1)), []
    w = a / 2

    def f(x):
        ax = abs(x)
        if ax <= mp.mpf(1) / 2 - w:
            return mp.mpf(1)
        return (1 + mp.cos(PI * (ax - (mp.mpf(1) / 2 - w)) * 2 / a)) / 2
    return f, [mp.mpf(1) / 2 - w]


for a in (0.0, 0.1, 0.5, 0.9, 1.0, 1.5):
    f, br = tukey_shape(a)
    continuous_cases("tukey", "tukey", {"alpha": a}, f, (9, 10), br)


# ----------------------------------------------------------------------------------------------
# webster-hamming (Doerry eq. 147): alpha(v) cos^v + (1 - alpha(v)) cos^(v+2)
# ----------------------------------------------------------------------------------------------
def webster_shape(v):
    v = mp.mpf(v)
    al = (2 + 3 * v + v * v) / (23 + 9 * v + v * v)
    return lambda x: al * mp.cos(PI * x) ** v + (1 - al) * mp.cos(PI * x) ** (v + 2)


for v in (-0.4, 0.0, 1.0, 3.0, 10.0, 100.0):
    continuous_cases("webster-hamming", "webster-hamming", {"v": v}, webster_shape(v), (9,), interior=(v < 0))


# ----------------------------------------------------------------------------------------------
# constants by quadrature (singular points): avci-nacaroglu, parzen-exponential, parzen-geometric
# ----------------------------------------------------------------------------------------------
for alpha in (0.5, 1.0, 2.0, 4.0):
    pa = PI * alpha
    continuous_cases("quadrature", "avci-nacaroglu", {"alpha": alpha},
                     lambda x, pa=pa: mp.exp(pa * (mp.sqrt(max(mp.mpf(0), 1 - 4 * x * x)) - 1)), (9,))

for alpha, r in ((3.0, 0.5), (1.0, 0.5), (2.0, 0.7), (1.0, 1.5), (0.5, 0.25)):
    continuous_cases("quadrature", "parzen-exponential", {"alpha": alpha, "r": r},
                     lambda x, a=alpha, r=r: mp.exp(-abs(2 * a * x) ** r), (9,), tag="-a%s-r%s" % (alpha, r))

for alpha, r in ((1.0, 0.5), (2.0, 0.8), (1.0, 3.0), (0.5, 0.25)):
    continuous_cases("quadrature", "parzen-geometric", {"alpha": alpha, "r": r},
                     lambda x, a=alpha, r=r: 1 / (1 + abs(2 * a * x) ** r), (9,), tag="-a%s-r%s" % (alpha, r))


# ----------------------------------------------------------------------------------------------
# knab: true limit at the ends (s = 0)
# ----------------------------------------------------------------------------------------------
def knab_shape(alpha):
    pa = PI * mp.mpf(alpha)

    def f(x):
        s = mp.sqrt(max(mp.mpf(0), 1 - 4 * x * x))
        return pa / mp.sinh(pa) if s == 0 else mp.sinh(pa * s) / (mp.sinh(pa) * s)
    return f


for alpha in (0.5, 1.0, 2.0, 3.0):
    continuous_cases("knab", "knab", {"alpha": alpha}, knab_shape(alpha), (2, 5, 8))


# ----------------------------------------------------------------------------------------------
# sum of cosines (positive magnitudes, centred form) without clipping
# ----------------------------------------------------------------------------------------------
def cossum(a):
    a = [mp.mpf(c) for c in a]
    return lambda x: sum(c * mp.cos(2 * PI * l * x) for l, c in enumerate(a)) / a[0]


FLATTOP = [0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368]
FTNI = [0.2810639, 0.5208972, 0.1980399]
HANN4 = [0.5, 0.5]
for coeffs in (FLATTOP, FTNI, [1.0, 1.5, 0.7, -0.3], HANN4):
    tag = "-%d-%s" % (len(coeffs), "%g" % coeffs[1])
    shape = cossum(coeffs)
    # `blackman-harris-family`: raw scale 1 / a0, peak = shape / shape(0)
    for N in (5, 8, 9):
        xs = grid(N)
        v = [shape(x) for x in xs]
        emit("family", "blackman-harris-family", {"coeffs": coeffs}, N, False, v, tag=tag)
        emit("family", "blackman-harris-family", {"coeffs": coeffs}, N, True, [a / shape(mp.mpf(0)) for a in v], tag=tag)
shape = cossum(FTNI)
for N in (2, 3, 8, 9, 64):
    xs = grid(N)
    v = [shape(x) for x in xs]
    emit("flat-top-3", "flat-top-3", {}, N, False, v)
    emit("flat-top-3", "flat-top-3", {}, N, True, [a / shape(mp.mpf(0)) for a in v])
shape = cossum(FLATTOP)
for N in (2, 5, 8, 9, 64):
    xs = grid(N)
    v = [shape(x) for x in xs]
    emit("flat-top-3", "flat-top", {}, N, False, v)
    emit("flat-top-3", "flat-top", {}, N, True, [a / shape(mp.mpf(0)) for a in v])

# scipy flattop (general_cosine) on the same coefficients: shape / sum(a)
for N in (5, 8, 9, 65, 100):
    sc = W.flattop(N, sym=True) / sum(FLATTOP)
    emit("family-scipy", "blackman-harris-family", {"coeffs": FLATTOP}, N, True, sc, tag="-scipy-flattop")


# ----------------------------------------------------------------------------------------------
# mottaghi-kashtiban-shayesteh: Doerry table 1, unclipped
# ----------------------------------------------------------------------------------------------
def mks_shape(N):
    a0 = mp.mpf("0.5363") - mp.mpf("0.14") / (N - 1)
    a = [a0, mp.mpf("0.996") - a0, mp.mpf(0), mp.mpf("0.004")]
    return a, (lambda x: sum(c * mp.cos(2 * PI * l * x) for l, c in enumerate(a)))


for N in (1, 2, 3, 4, 5, 6, 8, 9, 16, 65):
    if N == 1:
        emit("mks", "mottaghi-kashtiban-shayesteh", {}, 1, False, [mp.mpf(1)])
        emit("mks", "mottaghi-kashtiban-shayesteh", {}, 1, True, [mp.mpf(1)])
        continue
    a, shape = mks_shape(N)
    v = [shape(x) for x in grid(N)]
    emit("mks", "mottaghi-kashtiban-shayesteh", {}, N, False, [t / a[0] for t in v])
    emit("mks", "mottaghi-kashtiban-shayesteh", {}, N, True, [t / shape(mp.mpf(0)) for t in v])

# ----------------------------------------------------------------------------------------------
# vorbis: standard discrete window sin(pi/2 sin^2(pi (n + 1/2) / N))
# ----------------------------------------------------------------------------------------------
for N in (2, 3, 4, 5, 8, 9, 16, 64):
    v = [mp.sin(PI / 2 * mp.sin(PI * (n + mp.mpf(1) / 2) / N) ** 2) for n in range(N)]
    emit("vorbis", "vorbis", {}, N, False, v)
    mx = max(v)
    emit("vorbis", "vorbis", {}, N, True, [a / mx for a in v])


# ----------------------------------------------------------------------------------------------
# singla-singh: w = 2 (1 - I_t(m+1, m+1)), t = abs(2x)  (Singla and Singh eq. 8, 11)
# ----------------------------------------------------------------------------------------------
def singla_shape(m):
    return lambda x: 1 - mp.betainc(m + 1, m + 1, 0, abs(2 * x), regularized=True)


for order in (0, 1, 2, 3, 5, 8, 10, 13, 15, 19, 21, 30, 45, 60, 100, 1000):
    shape = singla_shape(order)
    for N in (17, 16):
        xs = grid(N)
        v = [shape(x) for x in xs]
        emit("singla-singh", "singla-singh", {"order": order}, N, False, [2 * a for a in v])
        emit("singla-singh", "singla-singh", {"order": order}, N, True, [a / shape(mp.mpf(0)) for a in v])


# ----------------------------------------------------------------------------------------------
# taylor: Doerry eq. 208-212 on the fastmath grid (endpoints included)
# ----------------------------------------------------------------------------------------------
def taylor_F(sll, nbar):
    eta = mp.mpf(10) ** (mp.mpf(abs(sll)) / 20)
    A = mp.acosh(eta) / mp.pi
    A2 = A * A
    s2 = mp.mpf(nbar) ** 2 / (A2 + (nbar - mp.mpf(1) / 2) ** 2)
    F = []
    for m in range(1, nbar):
        num = mp.mpf(1)
        den = mp.mpf(1)
        for i in range(1, nbar):
            num *= 1 - mp.mpf(m) ** 2 / (s2 * (A2 + (i - mp.mpf(1) / 2) ** 2))
            if i != m:
                den *= 1 - mp.mpf(m) ** 2 / mp.mpf(i) ** 2
        F.append((-1) ** (m + 1) * num / (2 * den))
    return F


for level, n in ((-50.0, 4), (-30.0, 2), (-80.0, 8), (-35.0, 3), (-100.0, 4), (-50.0, 2)):
    F = taylor_F(level, n)
    shape = lambda x, F=F: 1 + 2 * sum(Fm * mp.cos(2 * PI * (m + 1) * x) for m, Fm in enumerate(F))
    continuous_cases("taylor", "taylor", {"level": level, "n": n}, shape, (2, 3, 8, 9, 64), scale=mp.mpf(1),
                     tag="-%s-n%d" % (level, n))


# ----------------------------------------------------------------------------------------------
# discrete windows by the spectrum definition (mean 1):
#   w[n] = 1/N sum_{k=0}^{N-1} W_k cos(2 pi k (n - (N-1)/2) / N),  W_k = N P(x0 cos(pi k / N)) / P(x0)
# ----------------------------------------------------------------------------------------------
def spectrum_window(N, mu_db, P):
    eta = mp.mpf(10) ** (mp.mpf(mu_db) / 20)
    x0 = mp.cosh(mp.acosh(eta) / (N - 1))
    P0 = P(N - 1, x0)
    spec = [N * P(N - 1, x0 * mp.cos(PI * k / N)) / P0 for k in range(N)]
    c = mp.mpf(N - 1) / 2
    return [sum(spec[k] * mp.cos(2 * PI * k * (n - c) / N) for k in range(N)) / N for n in range(N)]


def emit_discrete(group, window, opts, N, w, tag=""):
    emit(group, window, opts, N, False, w, tag=tag)
    mx = max(abs(a) for a in w)
    emit(group, window, opts, N, True, [a / mx for a in w], tag=tag)


for level in (-40.0, -60.0, -80.0):
    for N in (8, 9, 33):
        emit_discrete("discrete", "dolph-chebyshev", {"level": level}, N,
                      spectrum_window(N, abs(level), lambda n, x: mp.chebyt(n, x)), tag="-%s" % level)
for level, alpha in ((-50.0, 0.5), (-50.0, 2.0), (-70.0, 3.0), (-60.0, -0.5)):
    for N in (8, 9, 33):
        emit_discrete("discrete", "ultraspherical", {"level": level, "alpha": alpha}, N,
                      spectrum_window(N, abs(level), lambda n, x, a=alpha: mp.gegenbauer(n, a, x)),
                      tag="-%s-a%s" % (level, alpha))
for level in (-50.0, -80.0):
    for N in (8, 9, 33):
        emit_discrete("discrete", "saramaki", {"level": level}, N,
                      spectrum_window(N, abs(level), lambda n, x: mp.gegenbauer(n, 1, x)), tag="-%s" % level)
        # Doerry eq. 266: mu = 1.0754 |L| + 1.7388
        emit_discrete("discrete", "legendre", {"level": level}, N,
                      spectrum_window(N, 1.0754 * abs(level) + 1.7388, lambda n, x: mp.legendre(n, x)),
                      tag="-%s" % level)


def sk_window(N):
    Nm = N - 1
    ep = mp.mpf("0.02") + mp.mpf("0.001") * Nm + 1 / (2 * mp.mpf(Nm) + 50)

    def sinc(z):
        return mp.mpf(1) if z == 0 else mp.sin(PI * z) / (PI * z)
    v = [ep if (n == 0 or n == Nm) else sinc((n - mp.mpf(Nm) / 2) / (mp.mpf("0.654") * Nm)) ** mp.mpf("2.5")
         for n in range(N)]
    s = sum(v)
    return [N * a / s for a in v]


for N in (2, 3, 8, 9, 33, 64):
    emit_discrete("discrete", "shayesteh-kashtiban", {}, N, sk_window(N))

# ----------------------------------------------------------------------------------------------
# scipy: dolph-chebyshev (chebwin), kaiser-bessel (kaiser), kaiser-bessel-derived, gaussian, tukey
# ----------------------------------------------------------------------------------------------
for level in (40.0, 60.0, 80.0, 100.0):
    for N in (8, 9, 64, 65):
        sc = W.chebwin(N, level)
        emit("scipy", "dolph-chebyshev", {"level": -level}, N, True, sc, tag="-scipy-%s" % level)

for alpha in (0.5, 1.0, 2.0, 3.0):
    for N in (8, 9, 64):
        sc = W.kaiser(N, math.pi * alpha)
        emit("scipy", "kaiser-bessel", {"alpha": alpha}, N, True, sc, tag="-scipy-%s" % alpha)

for alpha in (0.5, 1.0, 2.0):
    for N in (2, 4, 8, 16, 64):
        sc = W.kaiser_bessel_derived(N, math.pi * alpha)
        emit("scipy", "kaiser-bessel-derived", {"alpha": alpha}, N, False, sc, tag="-scipy-%s" % alpha)

for alpha in (1.0, 2.5):
    for N in (8, 9, 65):
        sc = W.gaussian(N, (N - 1) / (2 * alpha))
        emit("scipy", "gaussian", {"alpha": alpha}, N, True, sc, tag="-scipy-%s" % alpha)

for a in (0.1, 0.5, 0.9, 1.0, 1.5):
    for N in (8, 9, 64):
        emit("scipy", "tukey", {"alpha": a}, N, True, W.tukey(N, min(a, 1.0)), tag="-scipy-%s" % a)

# ----------------------------------------------------------------------------------------------
# Window additions (WP 5.1 - 5.4): cosine, kaiser, general-gaussian, the flat-top sets
# ----------------------------------------------------------------------------------------------
# cosine: sin(pi (n + 1/2) / N), discrete.  Default :normalize? is false (the formula); true is the maximum 1.
for N in (1, 2, 3, 4, 5, 8, 9, 16, 64, 65, 256):
    v = [mp.sin(PI * (n + mp.mpf(1) / 2) / N) for n in range(N)]
    emit("cosine", "cosine", {}, N, False, v)
    mx = max(v)
    emit("cosine", "cosine", {}, N, True, [a / mx for a in v])
# scipy.signal.windows.cosine: symmetric (sym=True) and periodic (sym=False)
for N in (1, 2, 3, 4, 5, 8, 9, 16, 64, 65, 256):
    emit("cosine-scipy", "cosine", {}, N, False, W.cosine(N, sym=True))
    if N >= 2:
        emit("cosine-scipy", "cosine", {}, N, False, W.cosine(N, sym=False), periodic=True)

# kaiser: scipy.signal.windows.kaiser(N, beta), symmetric and periodic; default normalisation (division by the
# middle value, the maximum 1 for odd N)
for beta in (0.0, 0.5, math.pi, 5.0, 8.6, 14.0, 30.0, 100.0):
    for N in (2, 3, 8, 9, 64, 65):
        emit("kaiser-scipy", "kaiser", {"beta": beta}, N, True, W.kaiser(N, beta, sym=True))
        emit("kaiser-scipy", "kaiser", {"beta": beta}, N, True, W.kaiser(N, beta, sym=False), periodic=True)
# kaiser: the raw scale (integral 1) from the definition I0(beta sqrt(1 - 4 x^2)) in mpmath
for beta in (0.0, 0.5, math.pi, 8.0, 20.0, 60.0):
    continuous_cases("kaiser-mp", "kaiser", {"beta": beta},
                     lambda x, b=mp.mpf(beta): mp.besseli(0, b * mp.sqrt(max(mp.mpf(0), 1 - 4 * x * x))),
                     (2, 5, 8, 9))

# general-gaussian: exp(-1/2 abs(x / sig)^(2 p)); scipy.signal.windows.general_gaussian(N, p, sig * (N - 1)),
# the periodic window of N samples is the symmetric window of N + 1 samples, so scipy gets sig * N for sym=False
for p in (0.5, 1.0, 1.5, 3.0, 8.0):
    for sig in (0.1, 0.25, 0.5):
        for N in (8, 9, 65):
            emit("general-gaussian-scipy", "general-gaussian", {"p": p, "sig": sig}, N, True,
                 W.general_gaussian(N, p, sig * (N - 1)))
            emit("general-gaussian-scipy", "general-gaussian", {"p": p, "sig": sig}, N, True,
                 W.general_gaussian(N, p, sig * N, sym=False), periodic=True)
# general-gaussian: the raw scale (integral 1) from the definition in mpmath
for p in (0.1, 0.25, 0.5, 1.0, 1.5, 4.0, 10.0):
    for sig in (0.05, 0.25, 0.5, 2.0):
        continuous_cases("general-gaussian-mp", "general-gaussian", {"p": p, "sig": sig},
                         lambda x, p=mp.mpf(p), s=mp.mpf(sig): mp.exp(-abs(x / s) ** (2 * p) / 2),
                         (2, 9), breaks=(sig,))

# exponential (scipy mapping in the docstring): scipy.signal.windows.exponential(N, tau=(N - 1) / (2 alpha)) with the
# center at the middle; the periodic window of N samples is the symmetric window of N + 1 samples, tau = N / (2 alpha)
for alpha in (0.5, 1.0, 3.0):
    for N in (8, 9, 65):
        emit("exponential-scipy", "exponential", {"alpha": alpha}, N, True,
             W.exponential(N, tau=(N - 1) / (2 * alpha), sym=True))
        emit("exponential-scipy", "exponential", {"alpha": alpha}, N, True,
             W.exponential(N, tau=N / (2 * alpha), sym=False), periodic=True)

# flat-top sets: sums of cosines with the positive magnitudes (the same lists as in generate_flattop_psl.py;
# hft169d with a_4 = 0.33673418, see the comment there)
FLATTOP_SETS = [
    ("flat-top-hp", ["1", "1.912510941", "1.079173272", "0.1832630879"]),
    ("flat-top-srs", ["1", "1.93", "1.29", "0.388", "0.028"]),
    ("flat-top-bk", ["1", "1.93", "1.29", "0.388", "0.0322"]),
    ("hft70", ["1", "1.90796", "1.07349", "0.18199"]),
    ("hft90d", ["1", "1.942604", "1.340318", "0.440811", "0.043097"]),
    ("hft95", ["1", "1.9383379", "1.3045202", "0.4028270", "0.0350665"]),
    ("hft116d", ["1", "1.9575375", "1.4780705", "0.6367431", "0.1228389", "0.0066288"]),
    ("hft144d", ["1", "1.96760033", "1.57983607", "0.81123644", "0.22583558", "0.02773848", "0.00090360"]),
    ("hft169d", ["1", "1.97441842", "1.65409888", "0.95788186", "0.33673418", "0.06364621", "0.00521942",
                 "0.00010599"]),
    ("hft196d", ["1", "1.979280420", "1.710288951", "1.081629853", "0.448734314", "0.112376628", "0.015122992",
                 "0.000871252", "0.000011896"]),
    ("hft223d", ["1", "1.98298997309", "1.75556083063", "1.19037717712", "0.56155440797", "0.17296769663",
                 "0.03233247087", "0.00324954578", "0.00013801040", "0.00000132725"]),
    ("hft248d", ["1", "1.985844164102", "1.791176438506", "1.282075284005", "0.667777530266", "0.240160796576",
                 "0.056656381764", "0.008134974479", "0.000624544650", "0.000019808998", "0.000000132974"]),
]
for name, cs in FLATTOP_SETS:
    shape = cossum(cs)    # a_0 = 1, so the raw scale is the sum itself (integral 1)
    for N in (2, 3, 8, 9, 64):
        v = [shape(x) for x in grid(N)]
        emit("flat-top-sets", name, {}, N, False, v)
        emit("flat-top-sets", name, {}, N, True, [a / shape(mp.mpf(0)) for a in v])
    # scipy.signal.windows.general_cosine(N, a) is the same sum (a = positive magnitudes); peak value sum(a)
    af = [float(c) for c in cs]
    for N in (8, 9, 65):
        emit("flat-top-sets-scipy", name, {}, N, True, W.general_cosine(N, af, sym=True) / sum(af))
        emit("flat-top-sets-scipy", name, {}, N, True, W.general_cosine(N, af, sym=False) / sum(af), periodic=True)

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w") as fh:
    fh.write(";; Generated by utils/fastmath/dev/generate_window_reference.py (mpmath %s, scipy %s, numpy %s).\n"
             % (mp.__version__, __import__("scipy").__version__, np.__version__))
    fh.write("{:cases [\n")
    for c in CASES:
        fh.write(" " + c + "\n")
    fh.write("]}\n")
print("wrote", OUT, len(CASES), "cases")
