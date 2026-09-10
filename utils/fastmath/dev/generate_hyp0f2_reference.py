"""One-off generator of independent reference values for `hypergeometric-0F2`
tests in `fastmath.special-test`.

Produces `test/resources/special/hyp0f2_reference.edn`, computed with
`mpmath.hyper([], [a, b], x)`:

  - `generic`        : `a`, `b` chosen so neither is a non-positive integer,
    `x` ranging over a domain empirically confirmed reliable for BOTH the
    positive-x (`maclaurin-0F2`) and negative-x (`weniger-0F2`) code paths
    (`|x| <= 5000`; see the note on `weniger-0F2`'s breakdown below).
  - `wide-positive`   : same `a`, `b` grid, `x` positive only, extended much
    further (up to 1e7) -- the positive-x path (`maclaurin-0F2`) stays
    accurate to double precision far beyond where the negative-x path
    breaks down.

Poles (`a` or `b` a non-positive integer, `x != 0`) are NOT included (mpmath
raises `ZeroDivisionError`/`pole in hypergeometric series` there, matching
fastmath's `##NaN`, added this session); they are covered by dedicated
edge-case assertions instead.

IMPORTANT, discovered while building this reference: `weniger-0F2` (used for
negative `x`) progressively loses accuracy as `|x|` grows, and the rate of
degradation depends heavily on `a`/`b`:
  - when `a` and `b` are BOTH positive, it stays accurate (~1e-10 relative
    or better) out to at least `|x| = 5000`, and only breaks down
    catastrophically (wrong order of magnitude, even wrong sign) beyond
    roughly `|x| ~ 20000-50000` (e.g. at `a=1.5, b=2.5`: `x=-20000` is
    already ~1e-5 relative error, `x=-50000` gives a positive ~7.6e16 where
    the true value is a negative ~-4.0e17);
  - when `a` or `b` is negative (even non-integer, i.e. nowhere near an
    actual pole), the safe range shrinks dramatically -- e.g. at
    `a=-1.7, b=-2.7`, `x=-5000` is already off by more than 100% (wrong
    sign); empirically, `|x| <= 200` stays accurate (~1e-8 relative or
    better) across the WHOLE `a`, `b` grid used here, including the
    negative cases.
This is a genuine numerical-stability limitation of the acceleration
algorithm itself, left undocumented/unguarded in the code (by user's
explicit choice this session -- treated as a known limitation, not fixed).
This generator therefore uses two different domains: `|x| <= 200` for the
full `a`, `b` grid (`generic`), and an extended `|x| <= 5000` negative-`x`
domain restricted to `a`, `b` BOTH positive (`negative-x-wide`), where that
range is confirmed safe.

Run with (from repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/clojure/fastmath && (cd /home/ts/penv && uv run python \
        utils/fastmath/dev/generate_hyp0f2_reference.py)
"""
import mpmath as mp

mp.mp.dps = 50

OUT = "/home/ts/clojure/fastmath/test/resources/special/hyp0f2_reference.edn"


def to_double(x):
    d = float(x)
    if d != d or d in (float('inf'), float('-inf')):
        return None
    return d


def edn_vec(xs):
    return "[" + " ".join(repr(x) for x in xs) + "]"


def build_3arg(fn, triples):
    out_a, out_b, out_x, out_r = [], [], [], []
    for a, b, x in triples:
        try:
            v = to_double(fn(a, b, x))
        except Exception:
            v = None
        if v is not None:
            out_a.append(a)
            out_b.append(b)
            out_x.append(x)
            out_r.append(v)
    return {"a": out_a, "b": out_b, "x": out_x, "ref": out_r}


def h(a, b, x):
    return mp.hyper([], [mp.mpf(a), mp.mpf(b)], mp.mpf(x))


a_grid = [-4.5, -2.3, -1.7, -0.5, 0.3, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]
b_grid = [-3.5, -2.7, 0.5, 1.0, 1.5, 2.5, 3.0, 5.0, 8.0]
a_grid_pos = [0.3, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]
b_grid_pos = [0.5, 1.0, 1.5, 2.5, 3.0, 5.0, 8.0]

X_SAFE_ALL = [-200.0, -100.0, -50.0, -20.0, -10.0, -5.0, -2.0, -1.0, -0.5,
             -0.1, -0.01, 0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0,
             50.0, 100.0, 200.0]
X_NEG_WIDE = [-5000.0, -2000.0, -1000.0, -500.0]
X_WIDE_POS = [10000.0, 50000.0, 100000.0, 500000.0, 1000000.0, 5000000.0,
             10000000.0]

generic_triples = []
for a in a_grid:
    for b in b_grid:
        for x in X_SAFE_ALL:
            generic_triples.append((a, b, x))
generic_block = build_3arg(h, generic_triples)

neg_wide_triples = []
for a in a_grid_pos:
    for b in b_grid_pos:
        for x in X_NEG_WIDE:
            neg_wide_triples.append((a, b, x))
neg_wide_block = build_3arg(h, neg_wide_triples)

wide_pos_triples = []
for a in a_grid:
    for b in b_grid:
        for x in X_WIDE_POS:
            wide_pos_triples.append((a, b, x))
wide_pos_block = build_3arg(h, wide_pos_triples)

with open(OUT, "w") as f:
    f.write("{:generic {:a " + edn_vec(generic_block["a"]) + "\n")
    f.write("           :b " + edn_vec(generic_block["b"]) + "\n")
    f.write("           :x " + edn_vec(generic_block["x"]) + "\n")
    f.write("           :ref " + edn_vec(generic_block["ref"]) + "}\n")
    f.write(" :negative-x-wide {:a " + edn_vec(neg_wide_block["a"]) + "\n")
    f.write("                   :b " + edn_vec(neg_wide_block["b"]) + "\n")
    f.write("                   :x " + edn_vec(neg_wide_block["x"]) + "\n")
    f.write("                   :ref " + edn_vec(neg_wide_block["ref"]) + "}\n")
    f.write(" :wide-positive {:a " + edn_vec(wide_pos_block["a"]) + "\n")
    f.write("                 :b " + edn_vec(wide_pos_block["b"]) + "\n")
    f.write("                 :x " + edn_vec(wide_pos_block["x"]) + "\n")
    f.write("                 :ref " + edn_vec(wide_pos_block["ref"]) + "}}\n")

print("generic          :", len(generic_block["a"]))
print("negative-x-wide  :", len(neg_wide_block["a"]))
print("wide-positive    :", len(wide_pos_block["a"]))
print("wrote", OUT)
