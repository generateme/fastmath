"""One-off generator of independent reference values for the three low-order
hypergeometric functions (`hypergeometric-0F0`, `hypergeometric-1F0`,
`hypergeometric-0F1`) tests in `fastmath.special-test`.

Produces `test/resources/special/hyp_low_order_reference.edn`, computed with
`mpmath`:

  - `0F0` : `mpmath.hyper([], [], x)` (== `exp(x)`), a grid of real `x`
    including values close to the double overflow/underflow boundaries.
  - `1F0` : `mpmath.hyper([a], [], x)` (== `(1-x)**(-a)` on the real line),
    for a grid of `a` (negative integers, where the series truncates to a
    polynomial and so is real for every `x`; and other real `a`, restricted
    to `x < 1` since `(1-x)**(-a)` for `x >= 1` and non-integer `-a` is
    complex -- fastmath, using `Math.pow`, returns `NaN` there, matching
    mpmath's own complex result on the real axis).
  - `0F1` : `mpmath.hyp0f1(a, x)`, for a grid of `a` (positive, and negative
    non-integers -- `a` a non-positive integer is a pole of the series,
    confirmed to raise `ZeroDivisionError` in mpmath itself, matching
    fastmath's `NaN` there) and a wide range of `x`, both signs.

Run with (from repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/clojure/fastmath && (cd /home/ts/penv && uv run python \
        utils/fastmath/dev/generate_hyp_low_order_reference.py)
"""
import mpmath as mp

mp.mp.dps = 50

OUT = "/home/ts/clojure/fastmath/test/resources/special/hyp_low_order_reference.edn"


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
        try:
            v = to_double(fn(x))
        except Exception:
            v = None
        if v is not None:
            out_x.append(x)
            out_r.append(v)
    return {"arg": out_x, "ref": out_r}


def build_2arg(fn, pairs_a, xs_per_a):
    out_a, out_x, out_r = [], [], []
    for a in pairs_a:
        for x in dedup(xs_per_a(a)):
            try:
                v = to_double(fn(a, x))
            except Exception:
                v = None
            if v is not None:
                out_a.append(a)
                out_x.append(x)
                out_r.append(v)
    return {"order": out_a, "arg": out_x, "ref": out_r}


# ---------------- 0F0 : all real x, incl. near over/underflow -------------
f00_xs = ([round(-50.0 + 0.5 * i, 6) for i in range(201)] +
         [700.0, 709.0, 709.5, 709.78, 709.79, 710.0, -700.0, -745.0,
          -745.13, -746.0, 0.0, -0.0])
f00_block = build_1arg(lambda x: mp.hyper([], [], mp.mpf(x)), f00_xs)

# ---------------- 1F0 : (a, x) grid, x < 1 for non-integer a --------------
f10_neg_int_orders = [-5.0, -4.0, -3.0, -2.0, -1.0, 0.0]
f10_other_orders = [-2.5, -1.5, -0.5, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]

f10_xs_wide = [-20.0, -10.0, -5.0, -2.0, -1.0, -0.5, -0.1, 0.0, 0.1, 0.5,
              0.9, 0.99, 0.999, 1.5, 2.0, 5.0, 10.0]  # safe for polynomial a
f10_xs_lt1 = [-20.0, -10.0, -5.0, -2.0, -1.0, -0.5, -0.1, 0.0, 0.1, 0.5,
             0.9, 0.99, 0.999, 0.9999]


def f10_fn(a, x):
    return mp.hyper([mp.mpf(a)], [], mp.mpf(x))


f10_block = build_2arg(f10_fn, f10_neg_int_orders, lambda a: f10_xs_wide)
f10_block_other = build_2arg(f10_fn, f10_other_orders, lambda a: f10_xs_lt1)
for k in ("order", "arg", "ref"):
    f10_block[k] = f10_block[k] + f10_block_other[k]

# ---------------- 0F1 : (a, x) grid, wide real domain ---------------------
f01_pos_orders = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0]
f01_neg_noninteger_orders = [-0.5, -1.5, -2.5, -3.5, -5.5]

f01_xs = ([round(-200.0 + 2.0 * i, 6) for i in range(0, 201)] +
         [-0.5, -0.1, -0.01, -0.001, 0.001, 0.01, 0.1, 0.5, 1.0, 300.0,
          400.0])


def f01_fn(a, x):
    return mp.hyp0f1(mp.mpf(a), mp.mpf(x))


f01_block = build_2arg(f01_fn, f01_pos_orders + f01_neg_noninteger_orders,
                       lambda a: f01_xs)

with open(OUT, "w") as f:
    f.write("{:0F0 {:arg " + edn_vec(f00_block["arg"]) + "\n")
    f.write("       :ref " + edn_vec(f00_block["ref"]) + "}\n")
    f.write(" :1F0 {:order " + edn_vec(f10_block["order"]) + "\n")
    f.write("       :arg " + edn_vec(f10_block["arg"]) + "\n")
    f.write("       :ref " + edn_vec(f10_block["ref"]) + "}\n")
    f.write(" :0F1 {:order " + edn_vec(f01_block["order"]) + "\n")
    f.write("       :arg " + edn_vec(f01_block["arg"]) + "\n")
    f.write("       :ref " + edn_vec(f01_block["ref"]) + "}}\n")

print("0F0:", len(f00_block["arg"]))
print("1F0:", len(f10_block["arg"]))
print("0F1:", len(f01_block["arg"]))
print("wrote", OUT)
