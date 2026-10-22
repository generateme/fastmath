"""Generator of reference data for the noise functions of `fastmath.random` (Random Audit, work package 9.4), step 2 of 2.

Step 1 is `generate_random_noise_reference.clj`: it dumps the seed dependent lookup tables of the library (permutation `perm`,
`perm12` = perm mod 12 and the value table `valueLUT`, all produced by `java.util.Random` inside `NoiseConfig`) to a JSON file.
This script reads them (they are the only thing taken from the library) and evaluates an INDEPENDENT Python re-implementation of
the algorithms from their published definitions:

    hash        lattice hash of the FastNoise / "coherent noise" family: `v = (1619 x) ^ (31337 y) ^ (6971 z) ^ seed ^ (1013 offset)`,
                `v = v * v * v * 60493`, `v ^ (v >> 13)`, 32-bit wrapping arithmetic.
    gradient    fade/lerp gradient noise (Perlin): 1-D gradients +-2 from `perm`, 2-D four diagonal gradients (hash & 3), 3-D the 12 edge
                gradients of Perlin 2002 chosen by `perm12[hash & 511]`; corner dot products blended with the fade curve
                (none = nearest cell corner * 0.5, linear t, hermite 3t^2 - 2t^3 (Perlin 1985), quintic 6t^5 - 15t^4 + 10t^3 (Perlin 2002)).
    value       lattice values from `valueLUT` blended with the same fade curves.
    simplex     Gustavson, "Simplex noise demystified" (2005): skew F2 = (sqrt 3 - 1)/2, unskew G2 = (3 - sqrt 3)/6, F3 = 1/3, G3 = 1/6,
                kernel radius^2 0.5 (2-D) and 0.6 (3-D), scale 70 / 32, 1-D kernel (1 - x^2)^4 scale 0.395 (the interpolation is ignored).
    fractals    libnoise style: single = one octave with offset 0; fbm = sum gain^k v_k / sum gain^k; billow = same with |2 v| - 1 per octave;
                ridged multifractal = sum lacunarity^-k (1 - |v_k|)^2 w_k with the feedback weight w_{k+1} = clamp(signal_k / gain, 0, 1),
                normalised by the sum of the spectral weights lacunarity^-k (Spec WP 3.3, NEW normaliser). Octave k uses the lattice offset
                `perm[k]` and the coordinates multiplied k times by the lacunarity (repeated multiplication, as in the library, so the
                coordinates are bit-identical). `normalize? true` maps to (s + 1)/2, ridged to s; raw ridged = 2 s - 1.

Cases: generators single fbm billow ridgemulti x (value and gradient with 4 interpolations, simplex with :hermite only) x dims 1 2 3 x
normalize? true/false = 216 cases, 12 points per case (the 4 points of the baseline `noise-points` of that dimension from
utils/fastmath/dev/generate_random_baseline.clj + 8 points of numpy PCG64: 4 in [-30, 30]^d, 2 in [-300, 300]^d, 2 in [-1000, 1000]^d; larger coordinates are avoided on purpose: the rounding of
the lattice coordinates (skew, 2.1^k scaling) amplifies differences of operation order, 1e-12 at |x| = 1e4).
Two configurations (as in the baseline):
    :normalize? true   {:seed 12345 :octaves 5 :lacunarity 2.1 :gain 0.45}
    :normalize? false  {:seed 12345 :octaves 6 :lacunarity 2.0 :gain 0.5}
No simplex case uses a non-hermite interpolation (it is ignored by the algorithm; the baseline file covers it).

Run (two steps, from the repository root):
    java -cp "$(cat /tmp/fm_cp.txt)" clojure.main utils/fastmath/dev/generate_random_noise_reference.clj [tables.json]
    cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_random_noise_reference.py [tables.json]
Default tables file: /tmp/fastmath_noise_tables.json.
Output: test/resources/random/noise_reference.edn (deterministic, byte-identical on rerun).
"""
import json
import math
import os
import platform
import sys

import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                   "test", "resources", "random", "noise_reference.edn")
TABLES = sys.argv[1] if len(sys.argv) > 1 else "/tmp/fastmath_noise_tables.json"

D = json.load(open(TABLES))
SEED = D["seed"]
PERM = D["perm"] + D["perm"]
PERM12 = D["perm12"] + D["perm12"]
LUT = D["valueLUT"]

CONFIGS = {True: dict(seed=SEED, octaves=5, lacunarity=2.1, gain=0.45),
           False: dict(seed=SEED, octaves=6, lacunarity=2.0, gain=0.5)}

BASELINE_POINTS = {1: [[0.3], [-4.7], [123.456], [-0.0001]],
                   2: [[0.3, 0.4], [-4.7, 2.2], [123.456, -78.9], [10.0, 10.0]],
                   3: [[0.3, 0.4, 0.5], [-4.7, 2.2, 9.9], [123.4, -78.9, 1.1], [0.5, 0.5, 0.5]]}


def i32(x):
    x &= 0xFFFFFFFF
    return x - (1 << 32) if x >= (1 << 31) else x


def hsh(off, *c):
    v = 0
    for m, x in zip((1619, 31337, 6971), c):
        v ^= i32(m * x)
    v ^= SEED
    v ^= i32(1013 * off)
    v = i32(v)
    v = i32(i32(i32(v * v) * v) * 60493)
    return i32(v ^ (v >> 13))  # python >> on negative ints is arithmetic, as in Java


def fade(t, it):
    if it == "quintic":
        return t * t * t * (t * (t * 6 - 15) + 10)
    if it == "hermite":
        return t * t * (3 - 2 * t)
    return t


def lerp(a, b, t):
    return a + t * (b - a)


fl = math.floor

GR3 = [(1, 1, 0), (-1, 1, 0), (1, -1, 0), (-1, -1, 0), (1, 0, 1), (-1, 0, 1), (1, 0, -1), (-1, 0, -1),
       (0, 1, 1), (0, -1, 1), (0, 1, -1), (0, -1, -1)]


# ---- gradient noise
def g1(off, x, it):
    x0 = fl(x)
    xd = x - x0

    def G(i, d):
        return d * (2.0, -2.0)[PERM[(i & 255) + off] & 1]
    if it == "none":
        return 0.5 * G(x0, xd)
    return lerp(G(x0, xd), G(x0 + 1, xd - 1), fade(xd, it))


def g2(off, x, y, it):
    x0, y0 = fl(x), fl(y)
    xd, yd = x - x0, y - y0
    gx = [1, -1, 1, -1]
    gy = [1, -1, -1, 1]

    def G(i, j, a, b):
        k = hsh(off, i, j) & 3
        return a * gx[k] + b * gy[k]
    if it == "none":
        return 0.5 * G(x0, y0, xd, yd)
    xs, ys = fade(xd, it), fade(yd, it)
    return lerp(lerp(G(x0, y0, xd, yd), G(x0 + 1, y0, xd - 1, yd), xs),
                lerp(G(x0, y0 + 1, xd, yd - 1), G(x0 + 1, y0 + 1, xd - 1, yd - 1), xs), ys)


def G3f(off, i, j, k, a, b, c):
    g = GR3[PERM12[hsh(off, i, j, k) & 511]]
    return a * g[0] + b * g[1] + c * g[2]


def g3(off, x, y, z, it):
    x0, y0, z0 = fl(x), fl(y), fl(z)
    xd, yd, zd = x - x0, y - y0, z - z0
    if it == "none":
        return 0.5 * G3f(off, x0, y0, z0, xd, yd, zd)
    xs, ys, zs = fade(xd, it), fade(yd, it), fade(zd, it)

    def L(dz):
        z1 = z0 + dz
        zz = zd - dz
        return lerp(lerp(G3f(off, x0, y0, z1, xd, yd, zz), G3f(off, x0 + 1, y0, z1, xd - 1, yd, zz), xs),
                    lerp(G3f(off, x0, y0 + 1, z1, xd, yd - 1, zz), G3f(off, x0 + 1, y0 + 1, z1, xd - 1, yd - 1, zz), xs), ys)
    return lerp(L(0), L(1), zs)


# ---- value noise
def v1(off, x, it):
    x0 = fl(x)

    def L(i):
        return LUT[PERM[(i & 255) + off]]
    return L(x0) if it == "none" else lerp(L(x0), L(x0 + 1), fade(x - x0, it))


def v2(off, x, y, it):
    x0, y0 = fl(x), fl(y)

    def L(i, j):
        return LUT[hsh(off, i, j) & 255]
    if it == "none":
        return L(x0, y0)
    xs, ys = fade(x - x0, it), fade(y - y0, it)
    return lerp(lerp(L(x0, y0), L(x0 + 1, y0), xs), lerp(L(x0, y0 + 1), L(x0 + 1, y0 + 1), xs), ys)


def v3(off, x, y, z, it):
    x0, y0, z0 = fl(x), fl(y), fl(z)

    def L(i, j, k):
        return LUT[hsh(off, i, j, k) & 255]
    if it == "none":
        return L(x0, y0, z0)
    xs, ys, zs = fade(x - x0, it), fade(y - y0, it), fade(z - z0, it)

    def P(dz):
        return lerp(lerp(L(x0, y0, z0 + dz), L(x0 + 1, y0, z0 + dz), xs),
                    lerp(L(x0, y0 + 1, z0 + dz), L(x0 + 1, y0 + 1, z0 + dz), xs), ys)
    return lerp(P(0), P(1), zs)


# ---- simplex noise (Gustavson 2005)
F2 = 0.5 * (math.sqrt(3) - 1)
G2 = (3 - math.sqrt(3)) / 6
F3 = 1 / 3
G3 = 1 / 6


def s1(off, x):
    i0 = fl(x)
    i1 = i0 + 1
    x0 = x - i0
    x1 = x0 - 1
    S = [1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8]

    def c(i, xd):
        t = 1 - xd * xd
        if t < 0:
            return 0
        t *= t
        return t * t * xd * S[PERM[(i & 255) + off] & 15]
    return 0.395 * (c(i0, x0) + c(i1, x1))


def s2(off, x, y):
    s = (x + y) * F2
    i = fl(x + s)
    j = fl(y + s)
    t = (i + j) * G2
    X0 = i - t
    Y0 = j - t
    x0 = x - X0
    y0 = y - Y0
    i1, j1 = (1, 0) if x0 > y0 else (0, 1)
    x1 = x0 - i1 + G2
    y1 = y0 - j1 + G2
    x2 = x0 - 1 + 2 * G2
    y2 = y0 - 1 + 2 * G2
    n = 0
    for (ci, cj, xx, yy) in ((i, j, x0, y0), (i + i1, j + j1, x1, y1), (i + 1, j + 1, x2, y2)):
        t = 0.5 - xx * xx - yy * yy
        if t >= 0:
            g = GR3[PERM12[hsh(off, ci, cj) & 255]]
            t *= t
            n += t * t * (g[0] * xx + g[1] * yy)
    return 70 * n


def s3(off, x, y, z):
    s = (x + y + z) * F3
    i = fl(x + s)
    j = fl(y + s)
    k = fl(z + s)
    t = (i + j + k) * G3
    x0 = x - (i - t)
    y0 = y - (j - t)
    z0 = z - (k - t)
    if x0 >= y0:
        if y0 >= z0:
            o = (1, 0, 0, 1, 1, 0)
        elif x0 >= z0:
            o = (1, 0, 0, 1, 0, 1)
        else:
            o = (0, 0, 1, 1, 0, 1)
    else:
        if y0 < z0:
            o = (0, 0, 1, 0, 1, 1)
        elif x0 < z0:
            o = (0, 1, 0, 0, 1, 1)
        else:
            o = (0, 1, 0, 1, 1, 0)
    i1, j1, k1, i2, j2, k2 = o
    cs = ((i, j, k, x0, y0, z0),
          (i + i1, j + j1, k + k1, x0 - i1 + G3, y0 - j1 + G3, z0 - k1 + G3),
          (i + i2, j + j2, k + k2, x0 - i2 + 2 * G3, y0 - j2 + 2 * G3, z0 - k2 + 2 * G3),
          (i + 1, j + 1, k + 1, x0 - 1 + 3 * G3, y0 - 1 + 3 * G3, z0 - 1 + 3 * G3))
    n = 0
    for (a, b, c, xx, yy, zz) in cs:
        t = 0.6 - xx * xx - yy * yy - zz * zz
        if t >= 0:
            g = GR3[PERM12[hsh(off, a, b, c) & 511]]
            t *= t
            n += t * t * (g[0] * xx + g[1] * yy + g[2] * zz)
    return 32 * n


def base(nt, it, d, off, p):
    if nt == "gradient":
        return [None, g1, g2, g3][d](off, *p, it)
    if nt == "value":
        return [None, v1, v2, v3][d](off, *p, it)
    return [None, s1, s2, s3][d](off, *p)


# ---- fractals
def fractal(gen, nt, it, d, p, oct_, lac, gain, norm):
    if gen == "single":
        v = base(nt, it, d, 0, p)
        return (v + 1) / 2 if norm else v
    if gen in ("fbm", "billow"):
        s = 0.0
        amp = 1.0
        tot = 0.0
        q = list(p)
        for k in range(oct_):
            v = base(nt, it, d, PERM[k], q)
            if gen == "billow":
                v = abs(2 * v) - 1
            s += amp * v
            tot += amp
            amp *= gain
            q = [c * lac for c in q]
        s /= tot
        return (s + 1) / 2 if norm else s
    # ridged multifractal, normaliser 1 / sum of spectral weights lacunarity^-k
    s = 0.0
    w = 1.0
    freq = 1.0
    wsum = 0.0
    q = list(p)
    for k in range(oct_):
        v = base(nt, it, d, PERM[k], q)
        sig = (1 - abs(v)) ** 2 * w
        sw = 1.0 / freq
        s += sig * sw
        wsum += sw
        w = min(1.0, max(0.0, sig / gain))
        freq *= lac
        q = [c * lac for c in q]
    s /= wsum
    return s if norm else 2 * s - 1


def num(x):
    return repr(float(x))


def main():
    rng = np.random.Generator(np.random.PCG64(20241021))
    points = {}
    for d in (1, 2, 3):
        extra = np.vstack([rng.uniform(-30, 30, (4, d)), rng.uniform(-300, 300, (2, d)), rng.uniform(-1000, 1000, (2, d))])
        points[d] = BASELINE_POINTS[d] + [[float(c) for c in r] for r in extra]

    cases = []
    for gen in ("single", "fbm", "billow", "ridgemulti"):
        for nt in ("value", "gradient", "simplex"):
            for it in (("none", "linear", "hermite", "quintic") if nt != "simplex" else ("hermite",)):
                for d in (1, 2, 3):
                    for norm in (True, False):
                        c = CONFIGS[norm]
                        vals = [fractal(gen, nt, it, d, p, c["octaves"], c["lacunarity"], c["gain"], norm) for p in points[d]]
                        cases.append((gen, nt, it, d, norm, vals))
    assert len(cases) == 216

    out = []
    out.append(':header {:tool "python %s numpy %s (PCG64 seed 20241021 for the extra points)"\n'
               '  :precision "IEEE double, independent re-implementation; tables (perm perm12 valueLUT) of seed 12345 dumped from NoiseConfig"\n'
               '  :domain "216 cases = 4 generators x (value, gradient with 4 interpolations; simplex :hermite) x dims 1 2 3 x normalize? true/false; '
               '12 points per case: the 4 baseline points of the dimension, 4 in [-30 30], 2 in [-300 300], 2 in [-1000 1000]; ridged with the new '
               'normaliser 1/sum(lacunarity^-k)"}' % (platform.python_version(), np.__version__))
    out.append(':configs {true {:seed %d :octaves %d :lacunarity %s :gain %s :normalize? true}\n'
               '  false {:seed %d :octaves %d :lacunarity %s :gain %s :normalize? false}}' % (
                   SEED, CONFIGS[True]["octaves"], num(CONFIGS[True]["lacunarity"]), num(CONFIGS[True]["gain"]),
                   SEED, CONFIGS[False]["octaves"], num(CONFIGS[False]["lacunarity"]), num(CONFIGS[False]["gain"])))
    out.append(":points {\n%s}" % "\n".join(
        " %d [%s]" % (d, "\n  ".join("[" + " ".join(num(c) for c in p) + "]" for p in points[d])) for d in (1, 2, 3)))
    cs = []
    for gen, nt, it, d, norm, vals in cases:
        cs.append(" {:generator :%s :noise-type :%s :interpolation :%s :dim %d :normalize? %s\n  :values [%s]}" % (
            gen, nt, it, d, "true" if norm else "false", " ".join(num(v) for v in vals)))
    out.append(":cases [\n%s]" % "\n".join(cs))
    with open(OUT, "w") as fh:
        fh.write("{" + "\n".join(out) + "}\n")


main()
