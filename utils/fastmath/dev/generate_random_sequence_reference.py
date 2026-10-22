"""Generator of reference data for the sequence generators of `fastmath.random` (Random Audit, work package 9.3).

Oracles (none of them reads fastmath or Commons Math):

    :sobol     `scipy.stats.qmc.Sobol(d, scramble=False).random(1000)` (Joe-Kuo direction numbers, index 0 is the origin).
               Dimensions 1 2 5 20 (all columns) and 1000 (first 3 columns).
    :halton    independent Python implementation of the Kocis-Whiten Halton sequence as defined by Commons Math
               `HaltonSequenceGenerator`: bases are the first primes, the digit `a` of base `b` and weight `w` is permuted to
               `(w * a) mod b`, weights 1 2 3 3 8 11 12 14 7 18 12 13 17 18 29 14 18 43 41 44 40 30 47 65 71 28 40 60 79 89 56 50
               52 61 108 56 66 63 60 66, index starts at 0 (the origin is the first point). It differs from the plain Halton
               (`scipy.stats.qmc.Halton`) for d >= 2. Dimensions 1 2 3 5 10 40, first 200 points.
    :r2        mpmath (50 digits): the positive root phi of x^(d+1) = x + 1, alpha_i = phi^-i (i = 1..d), for d = 1..15;
               the points are the closed form x_n = frac(0.5 + n * alpha), n = 1..100 (first point n = 1), for d = 1 2 5 15,
               rounded to double from the 50 digit value.
    :jitter    scale table of the jitter of Roberts, "A simple method to construct isotropic quasirandom blue-noise point
               sequences" (extremelearning.com.au): jitter offset of point i (0 based) is `scale(i) * u`, u uniform in [0, 1),
               `scale(i) = jitter * sqrt(pi) * f * d0 * (i + 1 - i0)^p` with article constants
                   R2      d0 = 0.76  i0 = 0.70  f = 0.25  p = -0.5
                   Halton  d0 = 0.90  i0 = 0.70  f = 0.25  p = -0.5
                   Sobol   d0 = 0.16  i0 = 0.58  f = 0.40  p = -0.2
               (article: k = lambda * delta * sqrt(pi) / (4 sqrt(i + 0.3)) for R2, i 0 based). Values for i = 0..4 and 1000 at
               jitter = 1.0, evaluated with mpmath (50 digits).

Run:  cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_random_sequence_reference.py
Output: test/resources/random/sequence_reference.edn (deterministic, byte-identical on rerun).
"""
import math
import os
import platform
import warnings

import mpmath as mp
import numpy as np
import scipy
from scipy.stats import qmc

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                   "test", "resources", "random", "sequence_reference.edn")

mp.mp.dps = 50
warnings.simplefilter("ignore")

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109,
          113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173]
WEIGHTS = [1, 2, 3, 3, 8, 11, 12, 14, 7, 18, 12, 13, 17, 18, 29, 14, 18, 43, 41, 44, 40, 30, 47, 65, 71, 28, 40, 60, 79, 89,
           56, 50, 52, 61, 108, 56, 66, 63, 60, 66]


def num(x):
    x = float(x)
    if math.isnan(x):
        return "##NaN"
    if math.isinf(x):
        return "##Inf" if x > 0 else "##-Inf"
    return repr(x)


def vec(xs):
    return "[" + " ".join(num(x) for x in xs) + "]"


def rows(a):
    return "[" + "\n   ".join(vec(r) for r in a) + "]"


def radical_inverse(k, base, weight):
    x, f = 0.0, 1.0 / base
    while k > 0:
        x += ((k % base) * weight % base) * f
        k //= base
        f /= base
    return x


def halton_kw(d, n):
    return [[radical_inverse(k, PRIMES[j], WEIGHTS[j]) for j in range(d)] for k in range(n)]


def r2_alpha(d):
    phi = mp.findroot(lambda x: x ** (d + 1) - x - 1, 1.5)
    return [phi ** (-(i + 1)) for i in range(d)]


def r2_points(d, n):
    al = r2_alpha(d)
    return [[float(mp.frac(mp.mpf("0.5") + k * a)) for a in al] for k in range(1, n + 1)]


def jitter_scales(d0, i0, f, p, idx):
    return [float(mp.sqrt(mp.pi) * mp.mpf(f) * mp.mpf(d0) * (i + 1 - mp.mpf(i0)) ** mp.mpf(p)) for i in idx]


def main():
    parts = []
    parts.append(':header {:tool "python %s numpy %s scipy %s mpmath %s"\n'
                 '  :precision "mpmath 50 digits (R2, jitter); IEEE double (Sobol, Halton)"\n'
                 '  :domain "Sobol d 1 2 5 20 and 1000 (3 columns), 1000 points; Halton-Kocis-Whiten d 1 2 3 5 10 40, 200 points; '
                 'R2 alpha d 1..15, points d 1 2 5 15 n 1..100; jitter scale i 0..4 and 1000"}'
                 % (platform.python_version(), np.__version__, scipy.__version__, mp.__version__))

    sobol = []
    for d, cols in [(1, 1), (2, 2), (5, 5), (20, 20), (1000, 3)]:
        pts = qmc.Sobol(d, scramble=False).random(1000)[:, :cols]
        sobol.append("  %d {:columns %d :points\n   %s}" % (d, cols, rows(pts)))
    parts.append(":sobol {\n%s}" % "\n".join(sobol))

    halton = []
    for d in [1, 2, 3, 5, 10, 40]:
        halton.append("  %d {:points\n   %s}" % (d, rows(halton_kw(d, 200))))
    parts.append(":halton {\n%s}" % "\n".join(halton))

    alphas = []
    for d in range(1, 16):
        al = r2_alpha(d)
        alphas.append("  %d {:alpha %s\n     :alpha-digits [%s]}" % (
            d, vec([float(a) for a in al]), " ".join('"%s"' % mp.nstr(a, 40) for a in al)))
    parts.append(":r2-alpha {\n%s}" % "\n".join(alphas))

    r2 = []
    for d in [1, 2, 5, 15]:
        r2.append("  %d {:points\n   %s}" % (d, rows(r2_points(d, 100))))
    parts.append(":r2 {\n%s}" % "\n".join(r2))

    idx = [0, 1, 2, 3, 4, 1000]
    jit = []
    for name, (d0, i0, f, p) in [(":r2", (0.76, 0.70, 0.25, -0.5)), (":halton", (0.90, 0.70, 0.25, -0.5)),
                                 (":sobol", (0.16, 0.58, 0.40, -0.2))]:
        jit.append("  %s {:constants {:d0 %s :i0 %s :f %s :p %s}\n   :indices %s\n   :scale %s}" % (
            name, num(d0), num(i0), num(f), num(p), "[" + " ".join(str(i) for i in idx) + "]",
            vec(jitter_scales(d0, i0, f, p, idx))))
    parts.append(":jitter-scale {\n%s}" % "\n".join(jit))

    with open(OUT, "w") as fh:
        fh.write("{" + "\n".join(parts) + "}\n")


main()
