"""One-off generator of independent reference values for the exponential- and
logarithmic-integral family (`E0`, `E1`, `Ein`, `En`, `Ei`, `li`, `Li`) tests
in `fastmath.special-test`.

Produces `test/resources/special/ei_reference.edn`, computed with `mpmath`:
  - `E0`   : `mpmath.expint(0, x)`, all real `x` (fastmath's `E0` is defined
    for every real `x`, including `+-Inf`).
  - `E1`   : `mpmath.expint(1, x)` / `mpmath.e1(x)`, restricted to `x >= 0`
    (for negative `x`, mpmath itself returns a complex value -- fastmath,
    like Julia's `SpecialFunctions.expint`, only accepts non-negative `x`).
  - `Ein`  : computed independently of fastmath's own formula --
    `gamma + log(x) - expint(1,x)`... note the SIGN: `Ein(x) = E1(x) +
    log(x) + gamma`, evaluated directly via `mpmath.expint(1,x) + log(x) +
    mpmath.euler` for `x > 0`, and via the alternating entire power series
    `sum_{k=1}^inf (-1)^(k+1) x^k / (k*k!)` (summed directly in mpmath, not
    ported from fastmath's Taylor-coefficient code) for `x <= 0` down to the
    fixed domain boundary `-2.15` (fastmath's `Ein` is undefined -- `NaN` --
    below that, a known, documented limitation from reusing `E1`'s Taylor
    coefficients, which are only accurate up to `|x| = 2.15`).
  - `En`   : `mpmath.expint(n, x)`, for a grid of integer orders (positive,
    negative, and zero -- real for `x < 0` too, since these reduce to
    elementary functions of `x` and `exp(-x)`) and non-integer orders
    (positive `x` only -- complex for `x < 0` per mpmath itself).
  - `Ei`   : `mpmath.ei(x)`, all real `x != 0`.
  - `li`   : `mpmath.li(x)`, `x > 0`, `x != 1` (singularity there).
  - `Li`   : `mpmath.li(x, offset=True)`, `x > 0`.

Run with (from repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/clojure/fastmath && (cd /home/ts/penv && uv run python \
        utils/fastmath/dev/generate_ei_reference.py)
"""
import mpmath as mp

mp.mp.dps = 50

OUT = "test/resources/special/ei_reference.edn"


def to_double(x):
    d = float(x)
    if d != d or d in (float('inf'), float('-inf')):
        return None
    return d


def edn_vec(xs):
    return "[" + " ".join(repr(x) for x in xs) + "]"


def dedup(xs):
    seen = set()
    out = []
    for x in xs:
        xf = float(x)
        if xf not in seen and xf == xf:
            seen.add(xf)
            out.append(xf)
    return out


def build_1arg(fn, xs):
    xs = dedup(xs)
    out_x, out_r = [], []
    for x in xs:
        v = to_double(fn(x))
        if v is not None:
            out_x.append(x)
            out_r.append(v)
    return {"arg": out_x, "ref": out_r}


# ------------- shared branch-threshold grids (mirrors E1's own thresholds) --
POS_THRESH_DENSE = []
for lo, hi, step in [(0.0001, 0.0044, 0.0002), (0.0044, 0.053, 0.002),
                     (0.053, 0.6, 0.01), (0.6, 2.15, 0.03),
                     (2.15, 4.0, 0.05), (4.0, 10.0, 0.15),
                     (10.0, 20.0, 0.3), (20.0, 200.0, 4.0),
                     (200.0, 800.0, 15.0)]:
    n = int(round((hi - lo) / step))
    POS_THRESH_DENSE += [round(lo + i * step, 6) for i in range(n + 1)]

for edge in (0.0044, 0.053, 0.6, 2.15, 4.0, 10.0, 20.0, 200.0):
    for eps in (-1e-6, -1e-9, 1e-9, 1e-6):
        POS_THRESH_DENSE.append(edge + eps)

POS_THRESH_DENSE = dedup([x for x in POS_THRESH_DENSE if x > 0])

# ---------------- E0 : all real x ----------------
e0_xs = [-x for x in POS_THRESH_DENSE] + POS_THRESH_DENSE + [0.0, -0.0]
e0_block = build_1arg(lambda x: mp.expint(0, mp.mpf(x)), e0_xs)

# ---------------- E1 : x >= 0 ----------------
e1_block = build_1arg(lambda x: mp.expint(1, mp.mpf(x)), POS_THRESH_DENSE)

# ---------------- Ein : x >= -2.15 ----------------
GAMMA = mp.euler


def ein_series(x):
    x = mp.mpf(x)
    s = mp.mpf(0)
    term = mp.mpf(1)
    for k in range(1, 300):
        term = term * (-x) / k
        contrib = -term / k
        s += contrib
        if k > 5 and abs(contrib) < abs(s) * mp.mpf('1e-45'):
            break
    return s


def ein_ref_fn(x):
    if x == 0.0:
        return mp.mpf(0)
    if x > 0:
        return mp.expint(1, mp.mpf(x)) + mp.log(mp.mpf(x)) + GAMMA
    return ein_series(x)


ein_neg_xs = [-x for x in POS_THRESH_DENSE if x <= 2.15]
ein_xs = ein_neg_xs + POS_THRESH_DENSE + [0.0]
ein_block = build_1arg(ein_ref_fn, ein_xs)

# ---------------- En : (n, x) pairs ----------------
int_orders = [-3, -2, -1, 0, 1, 2, 3, 5, 10]
frac_orders = [0.5, 1.5, 2.5, -0.5, -1.5]
en_xs_pos = [0.001, 0.01, 0.1, 0.5, 0.999, 1.0, 1.001, 1.5, 2.0, 2.15, 3.0,
            5.0, 10.0, 20.0, 50.0, 100.0, 300.0]
en_xs_neg = [-0.001, -0.01, -0.1, -0.5, -1.0, -1.5, -2.0, -3.0, -5.0, -10.0]

en_n, en_x, en_ref = [], [], []
for n in int_orders:
    for x in en_xs_pos:
        v = to_double(mp.expint(n, mp.mpf(x)))
        if v is not None:
            en_n.append(float(n)); en_x.append(x); en_ref.append(v)
    if n <= 0:
        for x in en_xs_neg:
            v = to_double(mp.expint(n, mp.mpf(x)))
            if v is not None:
                en_n.append(float(n)); en_x.append(x); en_ref.append(v)
for n in frac_orders:
    for x in en_xs_pos:
        v = to_double(mp.expint(mp.mpf(n), mp.mpf(x)))
        if v is not None:
            en_n.append(n); en_x.append(x); en_ref.append(v)

en_block = {"order": en_n, "arg": en_x, "ref": en_ref}

# ---------------- Ei : all real x != 0 ----------------
ei_xs = [-x for x in POS_THRESH_DENSE] + POS_THRESH_DENSE
ei_block = build_1arg(lambda x: mp.ei(mp.mpf(x)), ei_xs)

# ---------------- li : x > 0, x != 1 ----------------
li_xs = [0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999, 0.9999,
        1.0001, 1.001, 1.01, 1.1, 1.3, 1.5, 2.0, 3.0, 5.0, 10.0, 50.0,
        100.0, 1000.0, 1.0e6, 1.0e10, 1.0e20, 1.0e50, 1.0e100]
li_block = build_1arg(lambda x: mp.li(mp.mpf(x)), li_xs)

# ---------------- Li : x > 0 ----------------
li_offset_block = build_1arg(lambda x: mp.li(mp.mpf(x), offset=True), li_xs)

with open(OUT, "w") as f:
    f.write("{:E0 {:arg " + edn_vec(e0_block["arg"]) + "\n")
    f.write("      :ref " + edn_vec(e0_block["ref"]) + "}\n")
    f.write(" :E1 {:arg " + edn_vec(e1_block["arg"]) + "\n")
    f.write("      :ref " + edn_vec(e1_block["ref"]) + "}\n")
    f.write(" :Ein {:arg " + edn_vec(ein_block["arg"]) + "\n")
    f.write("       :ref " + edn_vec(ein_block["ref"]) + "}\n")
    f.write(" :En {:order " + edn_vec(en_block["order"]) + "\n")
    f.write("      :arg " + edn_vec(en_block["arg"]) + "\n")
    f.write("      :ref " + edn_vec(en_block["ref"]) + "}\n")
    f.write(" :Ei {:arg " + edn_vec(ei_block["arg"]) + "\n")
    f.write("      :ref " + edn_vec(ei_block["ref"]) + "}\n")
    f.write(" :li {:arg " + edn_vec(li_block["arg"]) + "\n")
    f.write("      :ref " + edn_vec(li_block["ref"]) + "}\n")
    f.write(" :Li {:arg " + edn_vec(li_offset_block["arg"]) + "\n")
    f.write("      :ref " + edn_vec(li_offset_block["ref"]) + "}}\n")

print("E0 :", len(e0_block["arg"]))
print("E1 :", len(e1_block["arg"]))
print("Ein:", len(ein_block["arg"]))
print("En :", len(en_block["arg"]))
print("Ei :", len(ei_block["arg"]))
print("li :", len(li_block["arg"]))
print("Li :", len(li_offset_block["arg"]))
print("wrote", OUT)
