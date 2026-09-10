"""One-off generator of independent reference values for `hypergeometric-2F0`
tests in `fastmath.special-test`.

Produces `test/resources/special/hyp2f0_reference.edn`, computed with
`mpmath.hyper([a, b], [], x)`:

  - `terminating` : `a` or `b` a non-positive integer -- the series always
    terminates to a finite polynomial in `x` (no denominator parameter to
    ever cause a pole), well-defined and REAL for every `x`, both signs.
  - `generic-negative-safe` : neither `a` nor `b` a non-positive integer,
    `|x| <= 0.5` -- high precision confirmed.
  - `generic-negative-moderate` : same `a`, `b` grid, `0.5 < |x| <= 5.0` --
    looser tolerance expected (see note below).

IMPORTANT, discovered while building this reference: `hypergeometric-2F0`'s
underlying series `sum_n (a)_n (b)_n x^n / n!` diverges for every `x != 0`
(radius of convergence `0`) UNLESS `a` or `b` is a non-positive integer (the
`terminating` case above, fixed this session to use a direct, exact
polynomial evaluation -- it previously had isolated `##NaN` glitches at
specific coincidental `(a, b, x)`, e.g. `a=-2.0, b=1.0, x=-1.0`). For the
generic case it is necessarily computed via resummation of that divergent
asymptotic series (`weniger-2F0`, untouched this session). Confirmed via
mpmath:
  - for `x < 0`, mpmath's own resummation (`mpmath.hyper`) agrees with
    fastmath's value to double precision for small `|x|`, but the agreement
    degrades progressively as `|x|` grows -- ~1e-15 relative at `|x|<=0.1`,
    ~1e-6 at `|x|<=0.5`, ~1e-2 at `|x|<=5`, over 100% (wrong order of
    magnitude) by `|x|~1000`. This is a genuine numerical-algorithm
    limitation (left undocumented/unguarded in the code, by user's explicit
    choice this session -- treated as known, not fixed), so the reference
    grid below does not extend past `|x| = 5.0`;
  - for `x > 0`, mpmath's resummation instead returns a genuinely COMPLEX
    number (confirmed directly, e.g. at `a=1.5, b=2.5, x=0.5`), and its
    real part does NOT match fastmath's real result there -- this part is
    expected, not a bug: resumming a divergent series has a
    branch-cut-like ambiguity, and different resummation conventions
    (Borel summation vs. the Pade/Weniger-type acceleration used here) can
    legitimately disagree once you cross to the other side of it. This
    generator therefore does not attempt to validate `x > 0` (generic case)
    at all.

Run with (from repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/clojure/fastmath && (cd /home/ts/penv && uv run python \
        utils/fastmath/dev/generate_hyp2f0_reference.py)
"""
import mpmath as mp

mp.mp.dps = 50

OUT = "/home/ts/clojure/fastmath/test/resources/special/hyp2f0_reference.edn"


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
            v = fn(a, b, x)
            if hasattr(v, "imag") and v.imag != 0:
                v = None
            else:
                v = to_double(v)
        except Exception:
            v = None
        if v is not None:
            out_a.append(a)
            out_b.append(b)
            out_x.append(x)
            out_r.append(v)
    return {"a": out_a, "b": out_b, "x": out_x, "ref": out_r}


def h(a, b, x):
    return mp.hyper([mp.mpf(a), mp.mpf(b)], [], mp.mpf(x))


X_BOTH_SIGNS = [-100.0, -50.0, -20.0, -10.0, -5.0, -2.0, -1.0, -0.5, -0.1,
               -0.01, 0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0,
               100.0]
X_NEG_SAFE = [-0.5, -0.3, -0.2, -0.1, -0.05, -0.01, -0.001, 0.0]
X_NEG_MODERATE = [-5.0, -3.0, -2.0, -1.5, -1.0, -0.7, -0.5]

term_orders = [0.0, -1.0, -2.0, -3.0, -5.0, -8.0]
other_pos = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
other_neg = [-4.5, -2.3, -1.7, -0.5]
other_grid = other_pos + other_neg

term_triples = []
for a in term_orders:
    for b in other_grid + term_orders:
        if a == b:
            continue
        for x in X_BOTH_SIGNS:
            term_triples.append((a, b, x))
term_block = build_3arg(h, term_triples)

generic_a = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, -0.5, -1.7, -2.3, -4.5]
generic_b = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, -0.5, -1.7, -2.3, -4.5]

safe_triples = []
mod_triples = []
for a in generic_a:
    for b in generic_b:
        for x in X_NEG_SAFE:
            safe_triples.append((a, b, x))
        for x in X_NEG_MODERATE:
            mod_triples.append((a, b, x))
safe_block = build_3arg(h, safe_triples)
mod_block = build_3arg(h, mod_triples)

with open(OUT, "w") as f:
    f.write("{:terminating {:a " + edn_vec(term_block["a"]) + "\n")
    f.write("               :b " + edn_vec(term_block["b"]) + "\n")
    f.write("               :x " + edn_vec(term_block["x"]) + "\n")
    f.write("               :ref " + edn_vec(term_block["ref"]) + "}\n")
    f.write(" :generic-negative-safe {:a " + edn_vec(safe_block["a"]) + "\n")
    f.write("                         :b " + edn_vec(safe_block["b"]) + "\n")
    f.write("                         :x " + edn_vec(safe_block["x"]) + "\n")
    f.write("                         :ref " + edn_vec(safe_block["ref"]) + "}\n")
    f.write(" :generic-negative-moderate {:a " + edn_vec(mod_block["a"]) + "\n")
    f.write("                             :b " + edn_vec(mod_block["b"]) + "\n")
    f.write("                             :x " + edn_vec(mod_block["x"]) + "\n")
    f.write("                             :ref " + edn_vec(mod_block["ref"]) + "}}\n")

print("terminating               :", len(term_block["a"]))
print("generic-negative-safe     :", len(safe_block["a"]))
print("generic-negative-moderate :", len(mod_block["a"]))
print("wrote", OUT)
