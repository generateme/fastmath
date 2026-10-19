"""One-off generator of independent reference values for the generalized
exponential integral `fastmath.special/En` (`fastmath.special.en-test`).

Produces `test/resources/special/en_reference.edn`, computed with `mpmath`
(`mpmath.expint(n, x)`, 50 digits; every value is cross-checked against a
70-digit evaluation and the generator aborts if the two differ by more than
1e-30 relative -- near-integer orders are the delicate part).  Arguments are
the exact doubles that are written to the file, so the Clojure test evaluates
`En` at exactly the point the oracle used.

Blocks (maps keyed by role, parallel vectors `:order`, `:arg`, `:ref`, except
`:grid` which is a row-major table to keep the file small):
  - `:grid`          orders 0.5 .. 3.5 step 0.05 (61), x 0.01 .. 2.99 step 0.01
                     (299); `:ref` is row-major, `(+ (* i 299) j)` for order i,
                     x j.  Contains the isolated points (1.55, 2.88),
                     (1.6, 2.66), (1.65, 2.46) where the old near-integer
                     Taylor branch fired accidentally (error up to 5.3).
  - `:near-integer`  orders `k +- d` for k in 1, 2, 3, 6, 11, d in 0.4 .. 1e-12,
                     plus `k` exactly, several x < 3 and x just below 3.
  - `:tiny-x`        x from 1e-300 to 1e-5 (finite references only).
  - `:branch-switch` x in [2.9, 3.1] step 0.01 and 3 +- 1e-9, 1e-6 (series for
                     x < 3, continued fraction otherwise), several orders.
  - `:large-x`       x from 3 to 700, positive, negative, fractional and
                     integer orders.
  - `:large-order`   fractional orders 4.5 .. 1200.5, x from 1e-5 to 10.
  - `:negative-order` negative fractional and negative integer orders, x > 0
                     and (negative integer orders only) x < 0.
  - `:posint`        positive integer orders (unchanged series), x > 0.
  - `:dispatch`      orders 0 and 1 (dispatch to E0 / E1), incl. x < 0 for 0.

Run with (from repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_en_reference.py
Regeneration is byte-identical (fixed point sets, `repr` of doubles).
"""
import os

import mpmath as mp

# repo root = four directory levels above this file (utils/fastmath/dev/)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
OUT = os.path.join(ROOT, "test/resources/special/en_reference.edn")


def to_double(v):
    d = float(v)
    if d != d or d in (float('inf'), float('-inf')):
        return None
    return d


def en_ref(n, x):
    """Reference `En(n, x)` as a double, or None when not finite/real."""
    vals = []
    for dps in (50, 70):
        mp.mp.dps = dps
        vals.append(mp.expint(mp.mpf(n), mp.mpf(x)))
    a, b = vals
    if isinstance(a, mp.mpc) or isinstance(b, mp.mpc):
        return None
    if abs(a - b) > mp.mpf('1e-30') * abs(b):
        raise RuntimeError("oracle unstable at n=%r x=%r: %s vs %s" % (n, x, a, b))
    return to_double(b)


def edn_vec(xs):
    return "[" + " ".join(repr(x) for x in xs) + "]"


def dedup_pairs(pairs):
    seen = set()
    out = []
    for p in pairs:
        key = (float(p[0]), float(p[1]))
        if key not in seen:
            seen.add(key)
            out.append(key)
    return out


def build(pairs):
    pairs = dedup_pairs(pairs)
    ns, xs, rs = [], [], []
    for n, x in pairs:
        r = en_ref(n, x)
        if r is not None:
            ns.append(n)
            xs.append(x)
            rs.append(r)
    return {"order": ns, "arg": xs, "ref": rs}


def cross(orders, xs):
    return [(float(n), float(x)) for n in orders for x in xs]


# ---------------- grid: 61 orders x 299 x ----------------
grid_orders = [round(0.5 + 0.05 * i, 2) for i in range(61)]
grid_xs = [round(0.01 + 0.01 * j, 2) for j in range(299)]
grid_ref = []
for n in grid_orders:
    for x in grid_xs:
        r = en_ref(n, x)
        assert r is not None
        grid_ref.append(r)

# ---------------- near-integer orders ----------------
near_ds = [0.4, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 1e-3, 1e-4,
           1e-6, 1e-9, 1e-12]
near_orders = []
for k in (1, 2, 3, 6, 11):
    near_orders.append(float(k))
    for d in near_ds:
        near_orders.append(k + d)
        near_orders.append(k - d)
# near the zero-order boundary as well (n = d and n = -d)
for d in near_ds:
    near_orders.append(d)
    near_orders.append(-d)
near_xs = [0.01, 0.1, 0.5, 1.0, 2.0, 2.9, 2.999999]
near_block = build(cross(near_orders, near_xs))

# ---------------- tiny x ----------------
tiny_orders = [0.5, 1.5, 2.0, 2.5, 3.0, 5.5, -0.5, -1.5, -2.0, -3.5]
tiny_xs = [1e-300, 1e-200, 1e-100, 1e-50, 1e-20, 1e-10, 1e-5]
tiny_block = build(cross(tiny_orders, tiny_xs))

# ---------------- branch switch at x = 3 ----------------
sw_xs = [round(2.9 + 0.01 * i, 2) for i in range(21)]
sw_xs += [3.0 - 1e-6, 3.0 - 1e-9, 3.0 + 1e-9, 3.0 + 1e-6]
sw_orders = [0.5, 1.5, 2.0, 2.5, 3.7, 6.0, 10.5, -0.5, -1.0, -2.5, 0.25, 0.75]
sw_block = build(cross(sw_orders, sw_xs))

# ---------------- large x ----------------
large_xs = [3.0, 3.5, 4.0, 5.0, 7.5, 10.0, 20.0, 50.0, 100.0, 300.0, 700.0]
large_orders = [0.5, 1.5, 2.0, 2.5, 5.0, 5.5, 10.5, 20.5, -0.5, -2.5, -10.5,
                -1.0, -3.0]
large_block = build(cross(large_orders, large_xs))

# ---------------- large fractional orders ----------------
big_orders = [4.5, 7.001, 20.5, 50.2, 99.7, 100.5, 150.3, 500.5, 1200.5]
big_xs = [1e-5, 0.01, 0.5, 1.0, 2.0, 2.9, 2.999, 3.0, 10.0]
big_block = build(cross(big_orders, big_xs))

# ---------------- negative orders ----------------
neg_orders_frac = [-0.25, -0.5, -1.5, -2.5, -5.5, -10.5]
neg_xs = [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 2.9, 3.0, 10.0, 50.0]
neg_pairs = cross(neg_orders_frac, neg_xs)
neg_int_orders = [-1.0, -2.0, -3.0, -5.0, -10.0]
neg_pairs += cross(neg_int_orders, neg_xs)
# negative integer orders are real for x < 0
neg_pairs += cross([0.0] + neg_int_orders,
                   [-0.001, -0.01, -0.1, -0.5, -1.0, -2.0, -3.0, -5.0, -10.0])
neg_block = build(neg_pairs)

# ---------------- positive integer orders ----------------
posint_orders = [2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0]
posint_xs = [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 2.9, 3.0, 5.0, 10.0, 50.0]
posint_block = build(cross(posint_orders, posint_xs))

# ---------------- dispatch: n = 0 and n = 1 ----------------
disp_xs = [0.001, 0.1, 0.5, 1.0, 2.0, 2.9, 3.0, 10.0, 100.0]
disp_pairs = cross([0.0, 1.0], disp_xs)
disp_pairs += cross([0.0], [-0.1, -1.0, -5.0])
disp_block = build(disp_pairs)

BLOCKS = [("near-integer", near_block), ("tiny-x", tiny_block),
          ("branch-switch", sw_block), ("large-x", large_block),
          ("large-order", big_block),
          ("negative-order", neg_block), ("posint", posint_block),
          ("dispatch", disp_block)]

with open(OUT, "w") as f:
    f.write("{:grid {:orders " + edn_vec(grid_orders) + "\n")
    f.write("       :xs " + edn_vec(grid_xs) + "\n")
    f.write("       :ref " + edn_vec(grid_ref) + "}")
    for name, b in BLOCKS:
        f.write("\n :" + name + " {:order " + edn_vec(b["order"]) + "\n")
        f.write("   :arg " + edn_vec(b["arg"]) + "\n")
        f.write("   :ref " + edn_vec(b["ref"]) + "}")
    f.write("}\n")

print("grid:", len(grid_orders), "x", len(grid_xs), "=", len(grid_ref))
for name, b in BLOCKS:
    print(name + ":", len(b["ref"]))
print("wrote", OUT)
