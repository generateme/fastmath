"""One-off generator of independent reference values for the `Ci`/`Cin`
(cosine integral / entire cosine integral) tests in `fastmath.special-test`.

Produces `test/resources/special/ci_reference.edn`, computed with `mpmath`:
  - `Ci` : `mpmath.ci` directly, restricted to `x >= 0` (fastmath's `Ci` only
    accepts non-negative `x` -- for negative `x`, `mpmath.ci` itself returns a
    complex value `Ci(|x|) + i*pi`, and Julia's `SpecialFunctions.cosint`
    raises a `DomainError`; fastmath follows the latter convention).
  - `Cin`: NOT taken from `mpmath.ci` composed with fastmath's own formula --
    computed independently here from mpmath's own Euler-Mascheroni constant,
    log and ci: `Cin(x) = gamma + log(|x|) - ci(|x|)` for `x != 0`, and `0`
    at `x = 0` (removable singularity). `Cin` is entire and even, so its
    domain is all of the reals, including negative `x` and `0`.

Domains mirror `generate_si_reference.py`'s density strategy: dense near
zero, dense brackets around every internal branch switch of the piecewise
implementation (`Ci`'s at `x = 3, 6, 12`; both share the `x*x` double-overflow
branch around `+-1.3407807929942596e154`), and a wide sparse sweep beyond
that out to very large `x`.

Run with (from repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_ci_reference.py
"""
import mpmath as mp

mp.mp.dps = 50

OUT = "test/resources/special/ci_reference.edn"

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


def dedup(xs):
    seen = set()
    out = []
    for x in xs:
        xf = float(x)
        if xf not in seen and xf == xf:
            seen.add(xf)
            out.append(xf)
    return out


# ---------------- Ci domain: x >= 0 only ----------------
ci_xs = []
ci_xs += [round(0.001 * i, 5) for i in range(1, 51)]      # 0.001 .. 0.05, near the x=0 singularity (excl. 0 itself)
ci_xs += [round(0.02 * i, 4) for i in range(1, 701)]      # 0.02 .. 14.0, dense across branches 1-3

for edge in (3.0, 6.0, 12.0):
    for eps in (-1e-6, -1e-9, 0.0, 1e-9, 1e-6):
        ci_xs.append(edge + eps)

ci_xs += [round(1.0 * i, 4) for i in range(14, 201)]      # 14 .. 200

for e in [i * 0.5 for i in range(6, 620)]:
    x = 10.0 ** e
    if x > OVERFLOW_X * 10:
        break
    ci_xs.append(x)

for factor in (0.999999, 0.9999999999, 1.0, 1.0000000001, 1.000001, 10.0, 1e10):
    ci_xs.append(OVERFLOW_X * factor)

ci_xs = dedup(ci_xs)

ci_ref = []
ci_xs_out = []
for x in ci_xs:
    v = to_double(mp.ci(mp.mpf(x)))
    if v is not None:
        ci_xs_out.append(x)
        ci_ref.append(v)

# ---------------- Cin domain: all reals (even, entire) ----------------
cin_xs = [0.0]
cin_xs += [round(0.001 * i, 5) for i in range(-50, 51)]
cin_xs += [round(0.02 * i, 4) for i in range(-700, 701)]

for edge in (3.0, 6.0, 12.0):
    for eps in (-1e-6, -1e-9, 0.0, 1e-9, 1e-6):
        cin_xs.append(edge + eps)
        cin_xs.append(-(edge + eps))

cin_xs += [round(1.0 * i, 4) for i in range(-200, 201)]

for e in [i * 0.5 for i in range(6, 620)]:
    x = 10.0 ** e
    if x > OVERFLOW_X * 10:
        break
    cin_xs.append(x)
    cin_xs.append(-x)

for factor in (0.999999, 0.9999999999, 1.0, 1.0000000001, 1.000001, 10.0, 1e10):
    cin_xs.append(OVERFLOW_X * factor)
    cin_xs.append(-OVERFLOW_X * factor)

cin_xs = dedup(cin_xs)

GAMMA = mp.euler


def cin_ref_fn(x):
    if x == 0.0:
        return mp.mpf(0)
    ax = abs(mp.mpf(x))
    return GAMMA + mp.log(ax) - mp.ci(ax)


cin_ref = []
cin_xs_out = []
for x in cin_xs:
    v = to_double(cin_ref_fn(x))
    if v is not None:
        cin_xs_out.append(x)
        cin_ref.append(v)

with open(OUT, "w") as f:
    f.write("{:Ci {:arg " + edn_vec(ci_xs_out) + "\n")
    f.write("      :ref " + edn_vec(ci_ref) + "}\n")
    f.write(" :Cin {:arg " + edn_vec(cin_xs_out) + "\n")
    f.write("       :ref " + edn_vec(cin_ref) + "}}\n")

print("Ci points:", len(ci_xs_out))
print("Cin points:", len(cin_xs_out))
print("wrote", OUT)
