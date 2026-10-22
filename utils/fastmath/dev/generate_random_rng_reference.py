"""Generator of reference data for the RNG layer of `fastmath.random` (Random Audit, work package 9.1).

Oracles (all outside the library, no value is read from fastmath or Commons Math):

    :mersenne  numpy `MT19937._legacy_seeding(np.array([hi, lo], uint32))` + `random_raw()`.
               This is the reference `init_by_array` + `genrand_int32` of Matsumoto & Nishimura; the 64 bit seed of
               Commons Math `MersenneTwister.setSeed(long)` is `init_by_array({seed >>> 32, seed & 0xffffffff})`.
    :jdk       `java.util.Random` re-implemented in Python from its specification (48 bit LCG
               `s = (s * 0x5DEECE66D + 0xB) mod 2^48`, scrambled seed `(seed ^ 0x5DEECE66D) mod 2^48`, legacy `nextInt(bound)`,
               `nextDouble = ((next(26) << 27) + next(27)) * 2^-53`, polar `nextGaussian`).
    :well512a  the WELL512a recurrence of Panneton, L'Ecuyer & Matsumoto (ACM TOMS 2006) including the published
               errata mask 0xda442d24, written from the paper. The paper has no seeding; the state fill below is the one of
               Commons Math `AbstractWell.setSeed(int[])` (written from the documented algorithm, not read from the
               library): `v[0..1] = {seed >>> 32, seed & 0xffffffff}` (signed 32 bit), then for `i = 2..15`
               `l = v[i - 2]` (sign extended to 64 bit), `v[i] = (int) ((1812433253 * (l ^ (l >> 30)) + i) & 0xffffffff)`,
               start index 0. The sign extension of `l` is the Commons quirk (the MT reference uses an unsigned word); if the
               library fill differed from this, the stream test fails.

Derived values follow the Commons Math `BitsStreamGenerator` formulas (MT, WELL) and the `java.util.Random` ones (JDK):

    nextInt()        next(32)
    nextLong()       (next(32) << 32) | (next(32) & 0xffffffff)         (BitsStreamGenerator)
                     (next(32) << 32) + next(32)                         (java.util.Random)
    nextDouble()     ((next(26) << 26) | next(26)) * 2^-52               (BitsStreamGenerator)
                     ((next(26) << 27) + next(27)) * 2^-53               (java.util.Random)
    nextFloat()      next(23) * 2^-23  (BitsStreamGenerator),  next(24) / 2^24 (java.util.Random)
    nextBoolean()    next(1) != 0
    nextInt(n)       power of two: (n * next(31)) >> 31, otherwise rejection on next(31)
    nextGaussian()   Box-Muller (BitsStreamGenerator), polar (java.util.Random)

fastmath shapes derived from them (Spec 1.1-1.5):
    (irandom r n)          nextInt(n)                       (irandom r mn mx)   mn + nextInt(mx - mn)
    (lrandom r n)          unbiased bounded long (rule below)
    (lrandom r mn mx)      mn + bounded(mx - mn)            wide path for spans >= 2^63 (rule below)
    (drandom r mn mx)      mn + (mx - mn) * nextDouble      (grandom r mean std)   mean + std * nextGaussian
    (brandom r)            nextBoolean                      (brandom r p)          nextDouble < p

Unbiased bounded long (Spec 1.1, an independent port of the rule, applied to the raw MT nextLong stream):
    mx > 0: rem = 2^64 mod mx; draw x (signed) until (x - MIN) <= 2^64 - 1 - rem; value = x mod mx (floor mod).
    Spec 1.2/1.3 wide path (span = mx - mn >= 2^63, i.e. the long difference is negative while mx > mn):
    draw x until unsigned(x) < unsigned(mx - mn); value = mn + x (wrapping).
Each stream is 20000 accepted draws from MT seed 42. The file stores the first 200 values, the number of raw draws
consumed and the SHA-256 (hex) of all 20000 values printed in decimal and joined by a newline (no trailing newline).

Run:  cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_random_rng_reference.py
Output: test/resources/random/rng_reference.edn (deterministic, byte-identical on rerun).
"""
import hashlib
import math
import os
import platform

import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                   "test", "resources", "random", "rng_reference.edn")

M32 = 0xFFFFFFFF
M64 = (1 << 64) - 1
MIN = -(1 << 63)
MAX = (1 << 63) - 1

SEEDS = [0, 1, 42, 123456789, 4294967303, -1, MIN, MAX, 1337, 281474976710657]


def s32(x):
    x &= M32
    return x - (1 << 32) if x >= (1 << 31) else x


def s64(x):
    x &= M64
    return x - (1 << 64) if x >= (1 << 63) else x


def seed_words(seed):
    u = seed & M64
    return (u >> 32) & M32, u & M32


# ---------------------------------------------------------------- raw streams, `next(bits)`

class MTStream:
    def __init__(self, seed):
        hi, lo = seed_words(seed)
        bg = np.random.MT19937()
        bg._legacy_seeding(np.array([hi, lo], dtype=np.uint32))
        self.bg = bg

    def next(self, bits):
        return int(self.bg.random_raw()) >> (32 - bits)


class WellStream:
    """WELL512a: recurrence of the paper, state fill of Commons Math (see the module docstring)."""

    def __init__(self, seed):
        hi, lo = seed_words(seed)
        v = [s32(hi), s32(lo)] + [0] * 14
        for i in range(2, 16):
            l = v[i - 2]  # signed, python >> is arithmetic like Java long >>
            v[i] = s32((1812433253 * (l ^ (l >> 30)) + i) & M32)
        self.v = [x & M32 for x in v]
        self.i = 0

    def next(self, bits):
        v, i = self.v, self.i
        z0 = v[(i + 15) & 15]
        vi, vi1, vi2 = v[i], v[(i + 13) & 15], v[(i + 9) & 15]
        z1 = (vi ^ ((vi << 16) & M32)) ^ (vi1 ^ ((vi1 << 15) & M32))
        z2 = vi2 ^ (vi2 >> 11)
        z3 = z1 ^ z2
        z4 = ((z0 ^ ((z0 << 2) & M32)) ^ (z1 ^ ((z1 << 18) & M32)) ^ ((z2 << 28) & M32)
              ^ (z3 ^ (((z3 << 5) & M32) & 0xDA442D24)))
        v[i] = z3
        v[(i + 15) & 15] = z4
        self.i = (i + 15) & 15
        return z4 >> (32 - bits)


class JDKStream:
    def __init__(self, seed):
        self.s = (seed ^ 0x5DEECE66D) & ((1 << 48) - 1)

    def next(self, bits):
        self.s = (self.s * 0x5DEECE66D + 0xB) & ((1 << 48) - 1)
        return self.s >> (48 - bits)


class Gen:
    """`kind` is 'bits' (Commons BitsStreamGenerator semantics: MT, WELL) or 'jdk' (java.util.Random)."""

    def __init__(self, stream, kind):
        self.st = stream
        self.kind = kind
        self.cached = None

    def nx(self, bits):
        v = self.st.next(bits)
        return s32(v) if bits == 32 else v

    def next_int(self):
        return self.nx(32)

    def next_long(self):
        hi = self.nx(32)
        lo = self.nx(32)
        if self.kind == 'bits':
            return s64((hi << 32) | (lo & M32))
        return s64((hi << 32) + lo)

    def next_double(self):
        if self.kind == 'bits':
            hi = self.nx(26) << 26
            lo = self.nx(26)
            return (hi | lo) * 2.0 ** -52
        return ((self.nx(26) << 27) + self.nx(27)) * 2.0 ** -53

    def next_float(self):
        if self.kind == 'bits':
            return self.nx(23) * 2.0 ** -23
        return self.nx(24) / float(1 << 24)

    def next_boolean(self):
        return self.nx(1) != 0

    def next_int_n(self, n):
        if (n & -n) == n:
            return (n * self.nx(31)) >> 31
        while True:
            bits = self.nx(31)
            val = bits % n
            if bits - val + (n - 1) < (1 << 31):
                return val

    def next_gaussian(self):
        if self.kind == 'bits':
            if self.cached is not None:
                v, self.cached = self.cached, None
                return v
            x = self.next_double()
            y = self.next_double()
            alpha = 2 * math.pi * x
            r = math.sqrt(-2 * math.log(y))
            self.cached = r * math.sin(alpha)
            return r * math.cos(alpha)
        if self.cached is not None:
            v, self.cached = self.cached, None
            return v
        while True:
            v1 = 2 * self.next_double() - 1
            v2 = 2 * self.next_double() - 1
            s = v1 * v1 + v2 * v2
            if 0 < s < 1:
                break
        mul = math.sqrt(-2 * math.log(s) / s)
        self.cached = v2 * mul
        return v1 * mul


ALGOS = [(':mersenne', MTStream, 'bits'), (':jdk', JDKStream, 'jdk'), (':well512a', WellStream, 'bits')]


def derive(make):
    """All derived shapes, each computed from a fresh generator `make()`."""

    def run(f, n):
        g = make()
        return [f(g) for _ in range(n)]

    return {
        ':ints': run(lambda g: g.next_int(), 24),
        ':longs': run(lambda g: g.next_long(), 6),
        ':doubles': run(lambda g: g.next_double(), 6),
        ':floats': run(lambda g: g.next_float(), 6),
        ':booleans': run(lambda g: g.next_boolean(), 40),
        ':gaussians': run(lambda g: g.next_gaussian(), 6),
        ':irandom-10': run(lambda g: g.next_int_n(10), 30),
        ':irandom-16': run(lambda g: g.next_int_n(16), 10),
        ':irandom-1e9': run(lambda g: g.next_int_n(1000000000), 10),
        ':irandom-range-m3-4': run(lambda g: -3 + g.next_int_n(7), 10),
        ':lrandom-1000': run(lambda g: g.next_long() % 1000, 10),
        ':lrandom-range-m5-5': run(lambda g: -5 + g.next_long() % 10, 10),
        ':drandom-range-m2-6': run(lambda g: -2.0 + (6.0 - -2.0) * g.next_double(), 6),
        ':grandom-10-2': run(lambda g: 10.0 + 2.0 * g.next_gaussian(), 4),
        ':brandom-p-0.3': run(lambda g: g.next_double() < 0.3, 20),
    }


# ---------------------------------------------------------------- unbiased bounded long

def bounded_long(g, mx):
    """Spec 1.1 for mx > 0. Returns (value, raw draws consumed)."""
    rem = (1 << 64) % mx
    limit = (1 << 64) - 1 - rem
    n = 0
    while True:
        x = g.next_long()
        n += 1
        if (x - MIN) <= limit:
            return x % mx, n


def wide_long(g, mn, mx):
    """Spec 1.2 wide path (span >= 2^63). Returns (value, raw draws consumed)."""
    diff = (mx - mn) & M64
    n = 0
    while True:
        x = g.next_long()
        n += 1
        if (x & M64) < diff:
            return s64(mn + x), n


def stream(f, count=20000):
    g = Gen(MTStream(42), 'bits')
    vals, consumed = [], 0
    for _ in range(count):
        v, n = f(g)
        vals.append(v)
        consumed += n
    return vals, consumed


# ---------------------------------------------------------------- printing

def num(x):
    if isinstance(x, bool):
        return "true" if x else "false"
    if isinstance(x, int):
        return str(x)
    x = float(x)
    if math.isnan(x):
        return "##NaN"
    if math.isinf(x):
        return "##Inf" if x > 0 else "##-Inf"
    return repr(x)


def digest(xs):
    return hashlib.sha256("\n".join(str(x) for x in xs).encode("ascii")).hexdigest()


def vec(xs):
    return "[" + " ".join(num(x) for x in xs) + "]"


def main():
    parts = []
    header = ('{:tool "python %s numpy %s"\n'
              '  :oracles {:mersenne "numpy MT19937 legacy init_by_array + random_raw"\n'
              '            :jdk "java.util.Random specification re-implemented (48 bit LCG)"\n'
              '            :well512a "WELL512a recurrence (Panneton et al. 2006, errata mask) + Commons AbstractWell state fill"}\n'
              '  :precision "IEEE double, 32/64 bit integers exact"\n'
              '  :domain "10 seeds x 24 raw ints plus derived shapes; unbiased/wide bounded long streams of 20000 draws (MT seed 42)"\n'
              '  :seeds %s}') % (platform.python_version(), np.__version__, vec(SEEDS))
    parts.append(":header " + header)

    for name, cls, kind in ALGOS:
        per_shape = {}
        for seed in SEEDS:
            d = derive(lambda: Gen(cls(seed), kind))
            for k, v in d.items():
                per_shape.setdefault(k, []).append((seed, v))
        lines = []
        for k, items in per_shape.items():
            lines.append("  %s {%s}" % (k, "\n   ".join("%d %s" % (s, vec(v)) for s, v in items)))
        parts.append("%s {\n%s}" % (name, "\n".join(lines)))

    unbiased = []
    for mx in [3 * 2 ** 61, 2 ** 62 + 1, 10 ** 18, 2 ** 32 + 5]:
        vals, consumed = stream(lambda g, mx=mx: bounded_long(g, mx))
        unbiased.append("  %d {:draws 20000 :raw-consumed %d\n   :sha256 \"%s\"\n   :head %s}" % (mx, consumed, digest(vals), vec(vals[:200])))
    parts.append(":lrandom-unbiased {\n%s}" % "\n".join(unbiased))

    wide = []
    for mn, mx in [(MIN, MAX), (-(2 ** 62), 2 ** 62), (MIN, 0), (-5 * 10 ** 18, 5 * 10 ** 18), (-1, MAX)]:
        vals, consumed = stream(lambda g, mn=mn, mx=mx: wide_long(g, mn, mx))
        wide.append("  [%d %d] {:draws 20000 :raw-consumed %d\n   :sha256 \"%s\"\n   :head %s}" % (mn, mx, consumed, digest(vals), vec(vals[:200])))
    parts.append(":lrandom-wide {\n%s}" % "\n".join(wide))

    with open(OUT, "w") as f:
        f.write("{" + "\n".join(parts) + "}\n")


main()
