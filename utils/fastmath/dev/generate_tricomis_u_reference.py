"""One-off generator of independent reference values for `tricomis-U` tests
in `fastmath.special-test`.

Produces `test/resources/special/tricomis_u_reference.edn`, computed with
`mpmath.hyperu(a, b, x)`, for a grid of `a`, `b` (including the two internal
special-cases `a = b` and `a = b - 1`) and `x > 0.0` only.

IMPORTANT: `x = 0.0` is deliberately NOT included in this bulk reference.
mpmath's own direct evaluation of `hyperu(a, b, 0)` for `b >= 1` was found
to be UNRELIABLE there: for `a=-0.5, b=3.0` it literally returns `+inf`,
but evaluating at a sequence of very small positive `x` (1e-3, 1e-6, 1e-10)
shows the true limit is clearly `-inf` (a consistently, increasingly
negative trend) -- i.e. mpmath's `+inf` at the exact singular point does
not reliably track the correct sign. The `x=0` boundary values used in this
session's fix and tests were instead derived from these small-x limit
trends by hand (and, for `a` a non-positive integer, from an independently
confirmed closed form `(-1)^n (b)_n`, matching such limits exactly at
several `(a, b)` pairs), not from a bulk mpmath `x=0` sweep.

A SEPARATE precision issue, also discovered while building this reference
(NOT the x=0 unreliability above): for the general branch (`a != b`,
`a != b-1`), `tricomis-U(a,b,x) = x^(-a) * hypergeometric-2F0(a, 1+a-b,
-1/x)`. For small `x`, `-1/x` is a LARGE-magnitude negative argument fed to
`hypergeometric-2F0`, which -- as documented separately in that function's
own test file section -- has a known, already-accepted precision
limitation there (accurate to ~1e-15 relative for small `|arg|`, degrading
to ~1e-2 by `|arg|~5`, catastrophically wrong by `|arg|~1000`). This is
therefore not a new bug in `tricomis-U` itself, just that limitation
surfacing through the `1/x` transform; this generator's smallest `x` is
`0.01` (`|arg| <= 100`) to stay clear of the worst of it, and the resulting
reference is compared with a moderately loosened tolerance in the tests
(observed max relative error ~1.9e-3 in the `x >= 0.01` domain used here).

Run with (from repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/clojure/fastmath && (cd /home/ts/penv && uv run python \
        utils/fastmath/dev/generate_tricomis_u_reference.py)
"""
import mpmath as mp

mp.mp.dps = 50

OUT = "/home/ts/clojure/fastmath/test/resources/special/tricomis_u_reference.edn"


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


def u(a, b, x):
    return mp.hyperu(mp.mpf(a), mp.mpf(b), mp.mpf(x))


a_grid = [-3.0, -2.0, -1.0, -0.5, 0.3, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0]
b_grid = [-2.0, -0.5, 0.3, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0]
x_grid = [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 300.0]

triples = []
for a in a_grid:
    for b in b_grid:
        for x in x_grid:
            triples.append((a, b, x))

# explicit a=b and a=b-1 coincidences (internal special-cased branches)
for a in [0.5, 1.0, 1.5, 2.0, 2.5, -1.5]:
    for x in x_grid:
        triples.append((a, a, x))          # a = b
        triples.append((a, a + 1.0, x))    # a = b - 1

block = build_3arg(u, triples)

with open(OUT, "w") as f:
    f.write("{:a " + edn_vec(block["a"]) + "\n")
    f.write(" :b " + edn_vec(block["b"]) + "\n")
    f.write(" :x " + edn_vec(block["x"]) + "\n")
    f.write(" :ref " + edn_vec(block["ref"]) + "}\n")

print("points:", len(block["a"]))
print("wrote", OUT)
