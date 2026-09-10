"""
Generates reference values for `fastmath.special/hypergeometric-pFq` using
mpmath's `hyper`. Writes an EDN file with several blocks, each a map of
{:ps [...] :qs [...] :x [...] :ref [...]} where :ps/:qs are vectors of
vectors (each entry is the numerator/denominator parameter vector for one
test point) so the test side can zip them with `hypergeometric-pFq`
directly.

Requires: mpmath (uv run python in /home/ts/penv).
"""
import mpmath as mp
import random
import json

mp.mp.dps = 40
random.seed(2024)


def edn_vec(xs):
    return "[" + " ".join(repr(float(x)) for x in xs) + "]"


def write_block(f, name, points):
    # points: list of (ps, qs, x, ref)
    f.write(f"  :{name}\n")
    f.write("  {:ps [" + " ".join(edn_vec(p[0]) for p in points) + "]\n")
    f.write("   :qs [" + " ".join(edn_vec(p[1]) for p in points) + "]\n")
    f.write("   :x [" + " ".join(repr(float(p[2])) for p in points) + "]\n")
    f.write("   :ref [" + " ".join(repr(float(p[3])) for p in points) + "]}\n\n")


def rand_generic(lo=-4.0, hi=4.0):
    return round(random.uniform(lo, hi), 2)


def mk_points(gen_fn, n, dps=None):
    out = []
    tries = 0
    while len(out) < n and tries < n * 20:
        tries += 1
        try:
            ps, qs, x = gen_fn()
            v = mp.hyper(ps, qs, x)
            if mp.im(v) != 0 or not mp.isfinite(v):
                continue
            out.append((ps, qs, x, float(v)))
        except Exception:
            continue
    return out


# 1) entire (p<=q), x>=0, wide magnitude range: MacLaurin path
def gen_entire_pos():
    p = random.randint(0, 3)
    q = random.randint(p, p + 2)
    ps = [rand_generic() for _ in range(p)]
    qs = [round(random.uniform(0.5, 6.0), 2) for _ in range(q)]  # avoid accidental poles
    x = round(random.uniform(0.0, 100.0), 2)
    return ps, qs, x


# 2) entire (p<=q), x<0, moderate magnitude: Weniger path
def gen_entire_neg():
    p = random.randint(0, 3)
    q = random.randint(p, p + 2)
    ps = [rand_generic() for _ in range(p)]
    qs = [round(random.uniform(0.5, 6.0), 2) for _ in range(q)]
    x = round(random.uniform(-50.0, -0.01), 2)
    return ps, qs, x


# 3) radius-1 (p=q+1), |x|<1: MacLaurin path
def gen_radius1_inside():
    q = random.randint(0, 2)
    p = q + 1
    ps = [rand_generic() for _ in range(p)]
    qs = [round(random.uniform(0.5, 6.0), 2) for _ in range(q)]
    x = round(random.uniform(-0.95, 0.7), 3)
    return ps, qs, x


# 4) radius-1 (p=q+1), |x|>=0.72 (incl. > 1 for x<0 side only, since x>1
#    generic leaves the real line for non-terminating params): Weniger path
def gen_radius1_outside():
    q = random.randint(0, 2)
    p = q + 1
    ps = [rand_generic() for _ in range(p)]
    qs = [round(random.uniform(0.5, 6.0), 2) for _ in range(q)]
    x = round(random.uniform(-5.0, -0.72), 3)
    return ps, qs, x


# 5) terminating: a numerator param a non-positive integer, wide x incl. huge
def gen_terminating():
    p = random.randint(1, 3)
    q = random.randint(0, 2)
    ps = [rand_generic() for _ in range(p)]
    ps[random.randrange(p)] = float(-random.randint(0, 6))
    qs = [round(random.uniform(0.5, 6.0), 2) for _ in range(q)]
    x = random.choice([
        round(random.uniform(-5, 5), 2),
        round(random.uniform(-1e5, 1e5), 1),
    ])
    return ps, qs, x


# 6) cancellation: one numerator exactly equals one denominator
def gen_cancellation():
    p = random.randint(1, 3)
    q = random.randint(1, 2)
    ps = [rand_generic() for _ in range(p)]
    qs = [round(random.uniform(0.5, 6.0), 2) for _ in range(q)]
    shared = rand_generic()
    ps[random.randrange(p)] = shared
    qs[random.randrange(q)] = shared
    x = round(random.uniform(-5, 5), 2)
    return ps, qs, x


blocks = {
    "entire-positive": mk_points(gen_entire_pos, 60),
    "entire-negative": mk_points(gen_entire_neg, 60),
    "radius1-inside": mk_points(gen_radius1_inside, 60),
    "radius1-outside": mk_points(gen_radius1_outside, 60),
    "terminating": mk_points(gen_terminating, 60),
    "cancellation": mk_points(gen_cancellation, 40),
}

for k, v in blocks.items():
    print(k, len(v))

with open("/home/ts/clojure/fastmath/test/resources/special/pfq_reference.edn", "w") as f:
    f.write("{\n")
    for name, points in blocks.items():
        write_block(f, name, points)
    f.write("}\n")

print("done")
