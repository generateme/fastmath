"""One-off generator of independent reference values for `whittaker-M` tests
in `fastmath.special-test`.

Produces `test/resources/special/whittaker_m_reference.edn`, computed with
`mpmath.whitm(kappa, mu, x)` (mpmath's own Whittaker M function, matching
fastmath's `(kappa, mu, x)` argument order exactly), over a grid of `kappa`,
`mu` (including negative, to exercise `kummers-M`'s pole structure through
the transformed parameters `mu-kappa+1/2`, `2mu+1`) and `x > 0.0` (`x = 0.0`
is a separate boundary case, covered by dedicated edge-case assertions
instead, since its correct value depends discontinuously on the sign of
`mu+0.5`). Points landing exactly on a genuine pole are skipped (mpmath
raises there; caught and dropped by `build_3arg` below).

IMPORTANT, discovered while building this reference: `whittaker-M` can
incorrectly return `##Inf` when `kummers-M`'s own intermediate value
individually overflows double range, even though the final SCALED product
(kummers-M's huge value times the tiny `z^2` prefactor) would be a valid
finite double -- e.g. `whittaker-M(-3.0, -2.3, 700.0)` returns `##Inf`, but
the true value is `~9.3e159`, well within double range (confirmed: 93 of
1315 points at `x` in `{700, 1000}` in an earlier, wider version of this
grid were affected; none below `x = 300`). This is a genuine
numerical-algorithm limitation (a real fix would need a log-space/rescaled
variant of `kummers-M`), left undocumented/unguarded in the code (by user's
explicit choice this session -- treated as a known limitation, not fixed).
This generator therefore stays within `x <= 300`, where the issue was not
observed.

Run with (from repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/clojure/fastmath && (cd /home/ts/penv && uv run python \
        utils/fastmath/dev/generate_whittaker_m_reference.py)
"""
import mpmath as mp

mp.mp.dps = 50

OUT = "/home/ts/clojure/fastmath/test/resources/special/whittaker_m_reference.edn"


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
    return mp.whitm(mp.mpf(kappa), mp.mpf(mu), mp.mpf(x))


kappa_grid = [-3.0, -1.5, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0, 5.0]
mu_grid = [-3.0, -2.3, -1.7, -0.9, -0.5, -0.1, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0]
X = [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 300.0]

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
