"""One-off generator of independent reference values for `digamma`/`trigamma`
tests in `fastmath.special-test`.

Produces `test/resources/special/digamma_trigamma_reference.edn`: a map with
`:arg`, `:digamma` and `:trigamma` parallel vectors, computed with
`scipy.special.digamma` / `scipy.special.polygamma(1, x)`.

Grids are chosen to exercise both the shifting-loop branch (x < 8 for
digamma, x < 10 for trigamma) and the asymptotic-series branch, on both the
positive domain and the negative (reflection-formula) domain, while avoiding
exact integer poles.

Run with (from repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_digamma_trigamma_reference.py
"""
import scipy.special as sp

OUT = "test/resources/special/digamma_trigamma_reference.edn"


def fmt(x):
    return repr(float(x))


def edn_vec(xs):
    return "[" + " ".join(fmt(x) for x in xs) + "]"


# positive grids
pos_small = [round(0.05 + 0.05 * i, 10) for i in range(0, 159)]           # 0.05 .. 7.95, step 0.05 (shifting-loop branch)
pos_med = [8.0 + i for i in range(0, 42)]                                 # 8 .. 49 (asymptotic branch, near threshold)
pos_large = [50.0, 100.0, 500.0, 1000.0, 1e4, 1e5, 1e6, 1e8, 1e10, 1e12]  # far asymptotic tail

# negative grids (offset to avoid exact integer poles)
neg_small = [-(0.05 + 0.1 * i) for i in range(0, 80)]                     # -0.05 .. -7.95 (step .1, .x5 offsets)
neg_med = [-(k + 0.37) for k in range(1, 50)]                             # -1.37 .. -49.37
neg_large = [-(500.37), -(1000.37), -(1e4 + 0.37), -(1e6 + 0.37), -(1e8 + 0.37)]

xs = pos_small + pos_med + pos_large + neg_small + neg_med + neg_large

dg = [sp.digamma(x) for x in xs]
tg = [sp.polygamma(1, x) for x in xs]

with open(OUT, "w") as f:
    f.write("{:arg " + edn_vec(xs) + "\n")
    f.write(" :digamma " + edn_vec(dg) + "\n")
    f.write(" :trigamma " + edn_vec(tg) + "}\n")

print("count:", len(xs))
print("min/max x:", min(xs), max(xs))
print("wrote", OUT)
