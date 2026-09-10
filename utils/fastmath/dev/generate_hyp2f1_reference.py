"""One-off generator of independent reference values for `hypergeometric-2F1`
(Gauss's hypergeometric function) tests in `fastmath.special-test`.

Produces `test/resources/special/hyp2f1_reference.edn`, computed with
`mpmath.hyp2f1(a, b, c, x)`:

  - `generic`     : a, b, c chosen so neither a nor b is a non-positive
    integer (the "normal", non-terminating series), x in [-5, 1) (the
    natural |x|<1 convergence domain plus the analytically-continued
    negative-x range, both confirmed real by mpmath here; x >= 1 is a
    genuine branch cut for generic parameters, see below).
  - `terminating` : a or b a non-positive integer -- an exact finite
    polynomial, valid for ANY real x (fixed this session, previously only
    handled for both a, b non-positive integers with |x| < 0.72); x
    ranges over both signs and past x=1.
  - `x1-boundary`  : x = 1.0 exactly, a grid split by the sign of c-a-b
    (Gauss's summation theorem: finite for c-a-b>0, diverges for
    c-a-b<=0).
  - `near-integer-diff` : (a, b) chosen so b-a is very close to an integer
    (both positive and negative), x on both sides of 1 -- this specific
    region contained TWO confirmed bugs this session (see below), so is
    exercised extra thoroughly.

For x > 1 with GENERIC (non-terminating) a, b, c, mpmath returns a
genuinely complex value (confirmed directly, e.g. at a=1.5, b=2.5, c=3.5,
x=5.0); this is a real branch cut, not a bug, and is not included in the
bulk reference (`build_3or4arg` drops complex results automatically).

TWO genuine bugs were found and fixed in `hypergeometric-2F1` (in
`fastmath.special.hypergeometric`) while building this reference, both
confirmed against mpmath:
  - A genuine INFINITE LOOP (not just a wrong value -- confirmed hung for
    multiple minutes before being force-killed) for real x that routes to
    the internal `inf-2F1` helper (roughly `|x| >= 1.39`) whenever it is
    invoked with `b < a` (root cause: an internal loop in a helper named
    `P`, bounded by `(== n m)` with `n` only ever increasing from `0`,
    never terminates for negative `m = round(b-a)`). This can happen even
    though the top-level function already normalizes `a <= b` once,
    because `general-2F1`'s own `c-a-b < 0` transformation can produce an
    internal recursive call with `b < a` again. Fixed by swapping `a`, `b`
    (exact, by `hypergeometric-2F1`'s own a<->b symmetry) at the start of
    `inf-2F1` whenever `b < a`, restoring `m >= 0`.
  - A first attempt at a blanket "b-a near an integer -> NaN" guard (before
    finding the true root cause above) was tried and reverted after it was
    found to cause 69 test regressions in `regularized-beta` (which
    legitimately routes through `inf-2F1` with near-integer differences
    that converge FINE) -- i.e. "b-a near an integer" alone does not
    predict a hang; the real condition is specifically `b < a` (fixed
    above).
  - `a` or `b` a non-positive integer (an exact finite polynomial, valid
    for any x) was only handled when BOTH were non-positive integers AND
    `|x| < 0.72`, giving `##NaN` for every other combination (e.g. only
    `a` non-positive, or `|x| >= 0.72`, or `x > 1`) even though the true
    value is perfectly finite and real. Fixed with a direct polynomial
    evaluation, added ahead of every other special-cased branch.

Run with (from repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/clojure/fastmath && (cd /home/ts/penv && uv run python \
        utils/fastmath/dev/generate_hyp2f1_reference.py)
"""
import mpmath as mp

mp.mp.dps = 50

OUT = "/home/ts/clojure/fastmath/test/resources/special/hyp2f1_reference.edn"


def to_double(x):
    d = float(x)
    if d != d or d in (float('inf'), float('-inf')):
        return None
    return d


def edn_vec(xs):
    return "[" + " ".join(repr(x) for x in xs) + "]"


def build_4arg(fn, quads):
    out_a, out_b, out_c, out_x, out_r = [], [], [], [], []
    for a, b, c, x in quads:
        try:
            v = fn(a, b, c, x)
            if hasattr(v, "imag") and v.imag != 0:
                v = None
            else:
                v = to_double(v)
        except Exception:
            v = None
        if v is not None:
            out_a.append(a)
            out_b.append(b)
            out_c.append(c)
            out_x.append(x)
            out_r.append(v)
    return {"a": out_a, "b": out_b, "c": out_c, "x": out_x, "ref": out_r}


def f(a, b, c, x):
    return mp.hyp2f1(mp.mpf(a), mp.mpf(b), mp.mpf(c), mp.mpf(x))


X_GENERIC = [-5.0, -2.0, -1.0, -0.5, -0.1, -0.01, 0.0, 0.01, 0.1, 0.3, 0.5,
            0.7, 0.9, 0.99]
generic_abc = [(0.5, 1.5, 3.0), (1.0, 2.0, 4.0), (1.5, 2.5, 3.5),
              (2.0, 3.0, 0.5), (0.3, 0.7, 1.5), (-0.5, 1.5, 2.5),
              (1.5, -0.5, 2.0), (1.0, 1.0, 0.5), (1.0, 1.5, 1.5),
              (1.0, 2.0, 2.0), (1.0, 1.5, 2.5), (1.0, 1.5, 2.0),
              (2.0, 2.0, 4.0)]

generic_quads = []
for a, b, c in generic_abc:
    for x in X_GENERIC:
        generic_quads.append((a, b, c, x))
generic_block = build_4arg(f, generic_quads)

# ---------------- terminating: a or b a non-positive integer --------------
X_TERM = [-5.0, -2.0, -1.0, -0.5, 0.0, 0.5, 0.9, 1.0, 1.5, 2.0, 5.0, 10.0]
term_abc = [(-3.0, 2.5, 3.5), (-2.0, -3.0, 3.5), (-1.0, 5.0, 2.0),
           (-4.0, 1.5, -3.0), (-2.0, 2.5, -3.0), (0.0, 2.5, 3.5)]
term_quads = []
for a, b, c in term_abc:
    for x in X_TERM:
        term_quads.append((a, b, c, x))
term_block = build_4arg(f, term_quads)

# ---------------- x=1 boundary, split by sign of c-a-b ---------------------
x1_abc = [(1.0, 1.5, 3.5), (0.5, 0.5, 3.0), (2.0, 1.0, 5.0),  # c-a-b>0
         (1.0, 1.5, 2.5), (0.5, 1.5, 2.0), (2.0, 3.0, 5.0)]   # c-a-b=0
x1_quads = [(a, b, c, 1.0) for a, b, c in x1_abc]
x1_block = build_4arg(f, x1_quads)

# ---------------- near-integer b-a, both signs, x incl. > 1 ---------------
nid_pairs = [(1.5, 2.5), (1.45, 2.5), (1.55, 2.5), (2.5, 1.5), (2.5, 1.45),
            (2.5, 1.55), (0.9, 3.1), (3.1, 0.9), (1.999, 3.0),
            (3.0, 1.999)]
nid_c = 3.5
X_NID = [-2.0, -0.5, 0.5, 0.9, 1.01, 1.1, 1.4, 1.5, 2.0, 5.0, 10.0]
nid_quads = []
for a, b in nid_pairs:
    for x in X_NID:
        nid_quads.append((a, b, nid_c, x))
nid_block = build_4arg(f, nid_quads)

with open(OUT, "w") as f2:
    f2.write("{:generic {:a " + edn_vec(generic_block["a"]) + "\n")
    f2.write("           :b " + edn_vec(generic_block["b"]) + "\n")
    f2.write("           :c " + edn_vec(generic_block["c"]) + "\n")
    f2.write("           :x " + edn_vec(generic_block["x"]) + "\n")
    f2.write("           :ref " + edn_vec(generic_block["ref"]) + "}\n")
    f2.write(" :terminating {:a " + edn_vec(term_block["a"]) + "\n")
    f2.write("               :b " + edn_vec(term_block["b"]) + "\n")
    f2.write("               :c " + edn_vec(term_block["c"]) + "\n")
    f2.write("               :x " + edn_vec(term_block["x"]) + "\n")
    f2.write("               :ref " + edn_vec(term_block["ref"]) + "}\n")
    f2.write(" :x1-boundary {:a " + edn_vec(x1_block["a"]) + "\n")
    f2.write("               :b " + edn_vec(x1_block["b"]) + "\n")
    f2.write("               :c " + edn_vec(x1_block["c"]) + "\n")
    f2.write("               :x " + edn_vec(x1_block["x"]) + "\n")
    f2.write("               :ref " + edn_vec(x1_block["ref"]) + "}\n")
    f2.write(" :near-integer-diff {:a " + edn_vec(nid_block["a"]) + "\n")
    f2.write("                     :b " + edn_vec(nid_block["b"]) + "\n")
    f2.write("                     :c " + edn_vec(nid_block["c"]) + "\n")
    f2.write("                     :x " + edn_vec(nid_block["x"]) + "\n")
    f2.write("                     :ref " + edn_vec(nid_block["ref"]) + "}}\n")

print("generic            :", len(generic_block["a"]))
print("terminating        :", len(term_block["a"]))
print("x1-boundary        :", len(x1_block["a"]))
print("near-integer-diff  :", len(nid_block["a"]))
print("wrote", OUT)
