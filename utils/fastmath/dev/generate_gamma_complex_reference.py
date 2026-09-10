"""
Generates reference values for `fastmath.special/log-gamma-complex` and
`fastmath.special/gamma-complex` using mpmath's `loggamma`/`gamma` with
complex arguments. Writes an EDN file with two blocks (:log-gamma, :gamma),
each {:z [...] :ref [...]} with [re im] pairs throughout.

Requires: mpmath (uv run python in /home/ts/penv).
"""
import mpmath as mp
import random
import math

mp.mp.dps = 30
random.seed(2468)


def edn_c(c):
    return "[" + repr(c.real) + " " + repr(c.imag) + "]"


def rand_z():
    mode = random.choice(["small", "medium", "large", "near-real-axis"])
    if mode == "small":
        return complex(round(random.uniform(-5, 8), 3), round(random.uniform(-5, 5), 3))
    if mode == "medium":
        return complex(round(random.uniform(-3, 15), 2), round(random.uniform(-15, 15), 2))
    if mode == "large":
        return complex(round(random.uniform(-10, 40), 1), round(random.uniform(-40, 40), 1))
    # near real axis but not exactly (avoid poles), small imaginary offset
    re = round(random.uniform(-8, 8), 2)
    im = round(random.choice([1, -1]) * random.uniform(0.001, 0.2), 4)
    return complex(re, im)


points_lg = []
points_g = []
tries = 0
while len(points_lg) < 90 and tries < 3000:
    tries += 1
    z = rand_z()
    try:
        lg = mp.loggamma(z)
        g = mp.gamma(z)
        if not (mp.isfinite(lg) and mp.isfinite(g)):
            continue
        points_lg.append((z, complex(lg)))
        points_g.append((z, complex(g)))
    except Exception:
        continue

print("log-gamma", len(points_lg), "gamma", len(points_g))

with open("/home/ts/clojure/fastmath/test/resources/special/gamma_complex_reference.edn", "w") as f:
    f.write("{\n")
    f.write("  :log-gamma\n")
    f.write("  {:z [" + " ".join(edn_c(z) for z, _ in points_lg) + "]\n")
    f.write("   :ref [" + " ".join(edn_c(v) for _, v in points_lg) + "]}\n\n")
    f.write("  :gamma\n")
    f.write("  {:z [" + " ".join(edn_c(z) for z, _ in points_g) + "]\n")
    f.write("   :ref [" + " ".join(edn_c(v) for _, v in points_g) + "]}\n\n")
    f.write("}\n")

print("done")
