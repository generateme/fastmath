"""
Generates reference values for `fastmath.special/tricomis-U-complex` using
mpmath's `hyperu` with complex arguments. Writes an EDN file with blocks:
:generic (Re(z)>=0, wide |z|), :negative-re-small (Re(z)<0, small |z|,
the reliable sub-region there), :integer-b (b an integer, Re(z)>=0 so the
asymptotic path is used), :zero-z (z=0 exactly, a,b real). Each block is
{:a [...] :b [...] :z [...] :ref [...]} with [re im] pairs throughout.

The reference is restricted to the confirmed-reliable domains documented in
tricomis-U-complex's own docstring (Re(z)<0 large |z| is a known,
undocumented-boundary, sporadic-instability gap, not included here).

Requires: mpmath (uv run python in /home/ts/penv).
"""
import mpmath as mp
import random

mp.mp.dps = 30
random.seed(13579)


def edn_c(c):
    return "[" + repr(c.real) + " " + repr(c.imag) + "]"


def rand_c(lo=-4.0, hi=4.0):
    return complex(round(random.uniform(lo, hi), 2), round(random.uniform(lo, hi), 2))


def mk(gen_fn, n):
    out = []
    tries = 0
    while len(out) < n and tries < n * 25:
        tries += 1
        try:
            a, b, z = gen_fn()
            v = mp.hyperu(a, b, z)
            if not mp.isfinite(v):
                continue
            out.append((a, b, z, complex(v)))
        except Exception:
            continue
    return out


def gen_generic():
    a = rand_c()
    b = rand_c()
    mag = random.uniform(0.05, 250.0)
    ang = random.uniform(-1.5, 1.5)  # keep Re(z) >= 0-ish
    z = mag * complex(round(abs(__import__("math").cos(ang)), 4),
                       round(__import__("math").sin(ang), 4))
    return a, b, z


def gen_negative_re_small():
    a = rand_c()
    b = rand_c()
    mag = random.uniform(0.05, 4.0)
    ang = random.uniform(-1.5, 1.5)
    z = mag * complex(-round(abs(__import__("math").cos(ang)), 4),
                       round(__import__("math").sin(ang), 4))
    return a, b, z


def gen_integer_b():
    a = rand_c()
    b = complex(float(random.randint(-4, 6)), 0.0)
    mag = random.uniform(0.05, 100.0)
    ang = random.uniform(-1.5, 1.5)
    z = mag * complex(round(abs(__import__("math").cos(ang)), 4),
                       round(__import__("math").sin(ang), 4))
    return a, b, z


def gen_zero_z():
    a = complex(round(random.uniform(-4, 4), 2), 0.0)
    b = complex(round(random.uniform(-4, 4), 2), 0.0)
    return a, b, complex(0.0, 0.0)


def write_block(f, name, points):
    f.write(f"  :{name}\n")
    f.write("  {:a [" + " ".join(edn_c(p[0]) for p in points) + "]\n")
    f.write("   :b [" + " ".join(edn_c(p[1]) for p in points) + "]\n")
    f.write("   :z [" + " ".join(edn_c(p[2]) for p in points) + "]\n")
    f.write("   :ref [" + " ".join(edn_c(p[3]) for p in points) + "]}\n\n")


blocks = {
    "generic": mk(gen_generic, 60),
    "negative-re-small": mk(gen_negative_re_small, 65),
    "integer-b": mk(gen_integer_b, 40),
    "zero-z": mk(gen_zero_z, 30),
}

for k, v in blocks.items():
    print(k, len(v))

with open("/home/ts/clojure/fastmath/test/resources/special/tricomis_u_complex_reference.edn", "w") as f:
    f.write("{\n")
    for name, points in blocks.items():
        write_block(f, name, points)
    f.write("}\n")

print("done")
