"""
Generates reference values for `fastmath.special/hypergeometric-pFq-complex`
using mpmath's `hyper` with complex arguments. Writes an EDN file with
several blocks, each a map of {:ps [...] :qs [...] :z [...] :ref [...]}
where :ps/:qs are vectors of vectors of [re im] pairs (one numerator/
denominator parameter vector per test point) and :z/:ref are vectors of
[re im] pairs, so the test side can zip them with
`hypergeometric-pFq-complex` directly.

Requires: mpmath (uv run python in /home/ts/penv).
"""
import mpmath as mp
import random

mp.mp.dps = 40
random.seed(4242)


def edn_c(c):
    return "[" + repr(c.real) + " " + repr(c.imag) + "]"


def edn_vec_c(xs):
    return "[" + " ".join(edn_c(x) for x in xs) + "]"


def write_block(f, name, points):
    f.write(f"  :{name}\n")
    f.write("  {:ps [" + " ".join(edn_vec_c(p[0]) for p in points) + "]\n")
    f.write("   :qs [" + " ".join(edn_vec_c(p[1]) for p in points) + "]\n")
    f.write("   :z [" + " ".join(edn_c(p[2]) for p in points) + "]\n")
    f.write("   :ref [" + " ".join(edn_c(p[3]) for p in points) + "]}\n\n")


def rand_generic(lo=-4.0, hi=4.0):
    return complex(round(random.uniform(lo, hi), 2), round(random.uniform(lo, hi), 2))


def rand_pos_denom():
    # avoid accidental poles/near-zero in denominator params
    return complex(round(random.uniform(0.5, 6.0), 2), round(random.uniform(-3.0, 3.0), 2))


def mk_points(gen_fn, n):
    out = []
    tries = 0
    while len(out) < n and tries < n * 20:
        tries += 1
        try:
            ps, qs, z = gen_fn()
            v = mp.hyper(ps, qs, z)
            if not mp.isfinite(v):
                continue
            out.append((ps, qs, z, complex(v)))
        except Exception:
            continue
    return out


# 1) entire (p<=q), Re(z)>0, wide magnitude: MacLaurin path
def gen_entire_pos():
    p = random.randint(0, 3)
    q = random.randint(p, p + 2)
    ps = [rand_generic() for _ in range(p)]
    qs = [rand_pos_denom() for _ in range(q)]
    ang = random.uniform(-1.3, 1.3)  # roughly toward positive real axis
    mag = random.uniform(0.0, 80.0)
    z = mag * complex(round(abs(__import__("math").cos(ang)), 4),
                       round(__import__("math").sin(ang), 4))
    return ps, qs, z


# 2) entire (p<=q), Re(z)<0, wide magnitude: Weniger path (confirmed safe to ~60)
def gen_entire_neg():
    p = random.randint(0, 3)
    q = random.randint(p, p + 2)
    ps = [rand_generic() for _ in range(p)]
    qs = [rand_pos_denom() for _ in range(q)]
    ang = random.uniform(-1.3, 1.3)
    mag = random.uniform(0.0, 55.0)
    z = mag * complex(-round(abs(__import__("math").cos(ang)), 4),
                       round(__import__("math").sin(ang), 4))
    return ps, qs, z


# 3) radius-1 (p=q+1), |z|<0.72: MacLaurin path
def gen_radius1_inside():
    q = random.randint(0, 2)
    p = q + 1
    ps = [rand_generic() for _ in range(p)]
    qs = [rand_pos_denom() for _ in range(q)]
    mag = random.uniform(0.0, 0.65)
    ang = random.uniform(0, 2 * __import__("math").pi)
    z = mag * complex(round(__import__("math").cos(ang), 4), round(__import__("math").sin(ang), 4))
    return ps, qs, z


# 4) radius-1 (p=q+1), |z|>=0.72: Weniger path, confirmed reliable on both
#    sides of the real axis at moderate magnitude
def gen_radius1_outside():
    q = random.randint(0, 1)  # keep p,q small: p=3,q=2 showed more failures
    p = q + 1
    ps = [rand_generic() for _ in range(p)]
    qs = [rand_pos_denom() for _ in range(q)]
    mag = random.uniform(0.72, 2.5)
    ang = random.uniform(0, 2 * __import__("math").pi)
    z = mag * complex(round(__import__("math").cos(ang), 4), round(__import__("math").sin(ang), 4))
    return ps, qs, z


# 5) terminating: a numerator param a non-positive integer, wide |z| incl. huge
def gen_terminating():
    p = random.randint(1, 3)
    q = random.randint(0, 2)
    ps = [rand_generic() for _ in range(p)]
    ps[random.randrange(p)] = complex(-float(random.randint(0, 6)), 0.0)
    qs = [rand_pos_denom() for _ in range(q)]
    mag = random.choice([random.uniform(0, 5), random.uniform(0, 1e5)])
    ang = random.uniform(0, 2 * __import__("math").pi)
    z = mag * complex(round(__import__("math").cos(ang), 4), round(__import__("math").sin(ang), 4))
    return ps, qs, z


# 6) cancellation: one numerator exactly equals one denominator (non-integer)
def gen_cancellation():
    p = random.randint(1, 2)
    q = random.randint(1, 2)
    ps = [rand_generic() for _ in range(p)]
    qs = [rand_pos_denom() for _ in range(q)]
    shared = rand_generic()
    ps[random.randrange(p)] = shared
    qs[random.randrange(q)] = shared
    z = rand_generic(-2.5, 2.5)
    return ps, qs, z


blocks = {
    "entire-positive": mk_points(gen_entire_pos, 50),
    "entire-negative": mk_points(gen_entire_neg, 50),
    "radius1-inside": mk_points(gen_radius1_inside, 50),
    "radius1-outside": mk_points(gen_radius1_outside, 70),
    "terminating": mk_points(gen_terminating, 50),
    "cancellation": mk_points(gen_cancellation, 55),
}

for k, v in blocks.items():
    print(k, len(v))

with open("/home/ts/clojure/fastmath/test/resources/special/pfq_complex_reference.edn", "w") as f:
    f.write("{\n")
    for name, points in blocks.items():
        write_block(f, name, points)
    f.write("}\n")

print("done")
