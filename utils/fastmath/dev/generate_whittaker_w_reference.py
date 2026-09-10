"""One-off generator of independent reference values for `whittaker-W` tests
in `fastmath.special-test`.

Produces `test/resources/special/whittaker_w_reference.edn`, computed with
`mpmath.whitw(kappa, mu, x)` (mpmath's own Whittaker W function, matching
fastmath's `(kappa, mu, x)` argument order exactly), over a grid of `kappa`,
`mu` (including negative, to exercise `tricomis-U`'s domain/pole structure
through the transformed parameters `mu-kappa+1/2`, `2mu+1`) and `x > 0.0`
(`x = 0.0` is a separate boundary case, covered by dedicated edge-case
assertions instead, since its correct value depends discontinuously on the
sign of `mu+0.5`, mirroring `whittaker-M`). Points landing exactly on a
genuine pole are skipped (mpmath raises there; caught and dropped by
`build_3arg` below).

Unlike `whittaker-M` (which showed a confirmed intermediate-overflow
limitation via `kummers-M`), no analogous OVERFLOW issue was found for
`whittaker-W` here: `tricomis-U` is the decaying/bounded solution of
Kummer's equation (its defining asymptotic property is `U(a,b,x) ~ x^(-a)`
as `x -> infinity`, a power law, not exponential growth like `kummers-M`),
so it does not tend to overflow the way `kummers-M` does. (A SEPARATE,
genuine `##NaN`-producing overflow WAS found and fixed this session in
`tricomis-U`'s own `a=b` branch for very large `x` -- see that function's
own test section; the `x` grid below already benefits from that fix.)

A different, NOT-fixed (documented as a known limitation) precision issue
was found instead: for SMALL `x`, `tricomis-U`'s general branch feeds a
LARGE-magnitude negative argument (`-1/x`) to `hypergeometric-2F0`, which
has its own already-documented precision limitation there (see that
function's test section). Through the `mu`, `kappa` -> `a`, `b` transform
used by `whittaker-W`, this surfaces at a LARGER `x` threshold than for
`tricomis-U` tested directly, because the `mu`, `kappa` grid used here
reaches larger `|a|`, `|b|` (e.g. `kappa=5.0, mu=-3.0` gives `a=-7.5,
b=-5.0`) -- confirmed: `x=0.01` still shows ~95% relative error for some
points in this grid, only dropping to ~1e-7 by `x=0.1` and to machine
precision by `x=0.2`. This generator therefore starts at `x=0.1`.

Run with (from repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/clojure/fastmath && (cd /home/ts/penv && uv run python \
        utils/fastmath/dev/generate_whittaker_w_reference.py)
"""
import mpmath as mp

mp.mp.dps = 50

OUT = "/home/ts/clojure/fastmath/test/resources/special/whittaker_w_reference.edn"


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


def w(kappa, mu, x):
    return mp.whitw(mp.mpf(kappa), mp.mpf(mu), mp.mpf(x))


kappa_grid = [-3.0, -1.5, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0, 5.0]
mu_grid = [-3.0, -2.3, -1.7, -0.9, -0.5, -0.1, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0]
X = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 300.0, 700.0, 1000.0]

triples = []
for kappa in kappa_grid:
    for mu in mu_grid:
        for x in X:
            triples.append((kappa, mu, x))
block = build_3arg(w, triples)

with open(OUT, "w") as f:
    f.write("{:kappa " + edn_vec(block["a"]) + "\n")
    f.write(" :mu " + edn_vec(block["b"]) + "\n")
    f.write(" :x " + edn_vec(block["x"]) + "\n")
    f.write(" :ref " + edn_vec(block["ref"]) + "}\n")

print("points:", len(block["a"]))
print("wrote", OUT)
