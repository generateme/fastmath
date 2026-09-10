"""One-off generator of independent reference values for the `Si` (sine
integral) test in `fastmath.special-test`.

Produces `test/resources/special/si_reference.edn`, computed with `mpmath`
(`mpmath.si`, mpmath's own name for the standard sine integral
Si(x) = integral_0^x sin(t)/t dt -- not to be confused with fastmath's
lowercase `si`, which is the shifted quantity `Si(x) - pi/2`; `si` is tested
in `fastmath.special-test` purely via that relation, no separate external
reference is needed for it).

Domain was chosen to exercise every branch of fastmath's piecewise Cephes-style
rational-polynomial implementation (`t = x*x`):
  - `t <= 36`  (|x| <= 6)   : power-series-like rational approximation
  - `t <= 144` (6 < |x| <= 12) : first asymptotic-form rational approximation
  - `t < Inf`  (|x| > 12)   : second asymptotic-form rational approximation,
    all the way up to where `x*x` itself overflows to double `Inf`
    (`x` around `+-1.3407807929942596e154`, i.e. `sqrt(Double.MAX_VALUE)`)
  - beyond that overflow threshold: exact `+-pi/2` (no series needed)

The grid below is dense near 0 and near each branch boundary (6.0, 12.0, and
the overflow threshold) to catch any discontinuity at the switches, plus a
wide, sparser sweep out to very large |x|.

Run with (from repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_si_reference.py
"""
import mpmath as mp

mp.mp.dps = 50

OUT = "test/resources/special/si_reference.edn"

OVERFLOW_X = 1.3407807929942596e154  # sqrt(Double.MAX_VALUE), where x*x -> Inf


def to_double(x):
    d = float(x)
    if d != d or d in (float('inf'), float('-inf')):
        return None
    return d


def fmt(x):
    return repr(x)


def edn_vec(xs):
    return "[" + " ".join(fmt(x) for x in xs) + "]"


xs = []

# dense near zero (both signs, incl. exact 0.0)
xs += [round(0.01 * i, 4) for i in range(-50, 51)]

# dense sweep across the first two branches and their boundary (0..14, both signs)
xs += [round(0.05 * i, 4) for i in range(0, 281)]      # 0 .. 14.0
xs += [round(-0.05 * i, 4) for i in range(1, 281)]     # -0.05 .. -14.0

# tight bracketing right around the two branch switches (t=36 -> x=6, t=144 -> x=12)
for edge in (6.0, 12.0):
    for eps in (-1e-6, -1e-9, 0.0, 1e-9, 1e-6):
        xs.append(edge + eps)
        xs.append(-(edge + eps))

# medium-to-large sweep (third branch), both signs
xs += [round(1.0 * i, 4) for i in range(14, 201)]      # 14 .. 200
xs += [round(-1.0 * i, 4) for i in range(14, 201)]

# wide log-spaced sweep out to the overflow threshold and just beyond it
import math
for e in [i * 0.5 for i in range(6, 620)]:  # 10^3 .. 10^309ish, but we clip at OVERFLOW_X*10
    x = 10.0 ** e
    if x > OVERFLOW_X * 10:
        break
    xs.append(x)
    xs.append(-x)

# explicit values bracketing the overflow threshold itself
for factor in (0.999999, 0.9999999999, 1.0, 1.0000000001, 1.000001, 10.0, 1e10):
    x = OVERFLOW_X * factor
    xs.append(x)
    xs.append(-x)

# de-duplicate while preserving order, drop anything mpmath/double can't represent sanely
seen = set()
xs_dedup = []
for x in xs:
    xf = float(x)
    if xf not in seen and xf == xf:  # skip accidental NaN
        seen.add(xf)
        xs_dedup.append(xf)
xs = xs_dedup

ref = []
xs_out = []
for x in xs:
    v = to_double(mp.si(mp.mpf(x)))
    if v is not None:
        xs_out.append(x)
        ref.append(v)

with open(OUT, "w") as f:
    f.write("{:Si {:arg " + edn_vec(xs_out) + "\n")
    f.write("      :ref " + edn_vec(ref) + "}}\n")

print("Si points:", len(xs_out))
print("wrote", OUT)
