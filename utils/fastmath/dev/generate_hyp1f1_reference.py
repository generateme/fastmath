"""One-off generator of independent reference values for `hypergeometric-1F1`
(Kummer's confluent hypergeometric function `M`, `kummers-M` in
`fastmath.special`) tests in `fastmath.special-test`.

Produces `test/resources/special/hyp1f1_reference.edn`, computed with
`mpmath.hyp1f1(a, b, x)`:

  - `generic`  : `a`, `b` chosen so neither is a non-positive integer and
    `a != b` -- the plain, pole-free regime, wide `x` range (`b` is never a
    non-positive integer here, so there is no pole regardless of `a`).
  - `equal`    : `a == b`, split into two mathematically distinct
    sub-regimes: `a` a non-positive integer (the series has a removable
    `0/0` coincidence at the term where both Pochhammer symbols vanish
    together; the correct value is the truncated exponential series, *not*
    `exp(x)`) and `a` anything else (where the series is genuinely,
    unambiguously `exp(x)`).
  - `neg-int-a`: `a` a non-positive integer, `b` not (any sign/magnitude);
    the series always terminates to a finite polynomial for any `x`
    (numerator vanishes identically before the denominator ever could).
  - `neg-int-both-safe` : both `a` and `b` non-positive integers with
    `a >= b` (`a` no more negative than `b`, `a != b`); still terminates
    safely (numerator vanishes at or before the point the denominator
    would).

`neg-int-both-pole` (`a < b`, both non-positive integers -- a genuine pole,
confirmed to raise `ZeroDivisionError`/`pole in hypergeometric series` in
mpmath itself, matching fastmath's `NaN`) and the positive-integer-`a`
with non-positive-integer-`b` pole case are NOT included in the bulk
reference (mpmath raises there); they are covered by dedicated edge-case
assertions instead.

Run with (from repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/clojure/fastmath && (cd /home/ts/penv && uv run python \
        utils/fastmath/dev/generate_hyp1f1_reference.py)
"""
import mpmath as mp

mp.mp.dps = 50

OUT = "/home/ts/clojure/fastmath/test/resources/special/hyp1f1_reference.edn"


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
    return mp.hyp1f1(mp.mpf(a), mp.mpf(b), mp.mpf(x))


X_WIDE = [-50.0, -20.0, -10.0, -5.0, -2.0, -1.0, -0.5, -0.1, -0.01, 0.0,
         0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
X_OVERFLOW = [100.0, 200.0, 500.0, 700.0, 709.0, 710.0, 720.0, 1000.0,
             2000.0, -700.0, -1000.0]

# ---------------- generic : a, b non-integer/positive, b never a
# ---------------- non-positive integer, a != b -----------------------------
generic_a = [-5.5, -3.5, -2.2, -1.7, -0.5, 0.3, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0,
            10.0]
generic_b = [-4.5, -2.3, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]

generic_triples = []
for a in generic_a:
    for b in generic_b:
        if a == b:
            continue
        for x in X_WIDE:
            generic_triples.append((a, b, x))
generic_block = build_3arg(h, generic_triples)

# ---------------- overflow probe: a,b generic small, x huge ----------------
overflow_triples = []
for a, b in [(2.0, 3.0), (0.5, 3.0), (5.0, 0.5), (1.0, 1.5), (-2.0, 3.0)]:
    for x in X_OVERFLOW:
        overflow_triples.append((a, b, x))
overflow_block = build_3arg(h, overflow_triples)

# ---------------- equal : a == b -------------------------------------------
equal_exp_orders = [0.5, 1.0, 1.5, 2.0, 3.5, 5.0, -0.5, -1.5, -2.5]  # exp(x)
equal_trunc_orders = [-1.0, -2.0, -3.0, -4.0, -5.0, -7.0, -10.0]  # truncated

equal_triples = []
for a in equal_exp_orders + equal_trunc_orders:
    for x in X_WIDE:
        equal_triples.append((a, a, x))
equal_block = build_3arg(h, equal_triples)

# ---------------- neg-int-a : a non-positive integer, b generic ------------
neg_int_a_orders = [-1.0, -2.0, -3.0, -5.0, -8.0]
neg_int_a_bs = [0.5, 1.5, 3.0, -2.5, -7.5]

neg_int_a_triples = []
for a in neg_int_a_orders:
    for b in neg_int_a_bs:
        for x in X_WIDE:
            neg_int_a_triples.append((a, b, x))
neg_int_a_block = build_3arg(h, neg_int_a_triples)

# ---------------- neg-int-both-safe : both non-positive integers, a >= b,
# ---------------- a != b ----------------------------------------------------
safe_pairs = [(-1.0, -3.0), (-1.0, -5.0), (-2.0, -5.0), (-2.0, -7.0),
             (-3.0, -8.0), (-1.0, -2.0), (-4.0, -9.0)]

safe_triples = []
for a, b in safe_pairs:
    for x in X_WIDE:
        safe_triples.append((a, b, x))
safe_block = build_3arg(h, safe_triples)

with open(OUT, "w") as f:
    f.write("{:generic {:a " + edn_vec(generic_block["a"]) + "\n")
    f.write("           :b " + edn_vec(generic_block["b"]) + "\n")
    f.write("           :x " + edn_vec(generic_block["x"]) + "\n")
    f.write("           :ref " + edn_vec(generic_block["ref"]) + "}\n")
    f.write(" :overflow {:a " + edn_vec(overflow_block["a"]) + "\n")
    f.write("            :b " + edn_vec(overflow_block["b"]) + "\n")
    f.write("            :x " + edn_vec(overflow_block["x"]) + "\n")
    f.write("            :ref " + edn_vec(overflow_block["ref"]) + "}\n")
    f.write(" :equal {:a " + edn_vec(equal_block["a"]) + "\n")
    f.write("         :b " + edn_vec(equal_block["b"]) + "\n")
    f.write("         :x " + edn_vec(equal_block["x"]) + "\n")
    f.write("         :ref " + edn_vec(equal_block["ref"]) + "}\n")
    f.write(" :neg-int-a {:a " + edn_vec(neg_int_a_block["a"]) + "\n")
    f.write("             :b " + edn_vec(neg_int_a_block["b"]) + "\n")
    f.write("             :x " + edn_vec(neg_int_a_block["x"]) + "\n")
    f.write("             :ref " + edn_vec(neg_int_a_block["ref"]) + "}\n")
    f.write(" :neg-int-both-safe {:a " + edn_vec(safe_block["a"]) + "\n")
    f.write("                     :b " + edn_vec(safe_block["b"]) + "\n")
    f.write("                     :x " + edn_vec(safe_block["x"]) + "\n")
    f.write("                     :ref " + edn_vec(safe_block["ref"]) + "}}\n")

print("generic           :", len(generic_block["a"]))
print("overflow          :", len(overflow_block["a"]))
print("equal             :", len(equal_block["a"]))
print("neg-int-a         :", len(neg_int_a_block["a"]))
print("neg-int-both-safe :", len(safe_block["a"]))
print("wrote", OUT)
