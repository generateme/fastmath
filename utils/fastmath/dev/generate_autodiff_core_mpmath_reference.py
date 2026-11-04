# Reference values for derivative rules of fastmath.calculus.autodiff which have no rule in Julia (DiffRules.jl,
# LogExpFunctions.jl): first and second derivatives of the exact functions.
#
# - sinc and exprel around 0, where the rules switch between a Taylor series and the direct forms (the direct forms
#   cancel near 0): 0 and +/- 10^k for k from -10 to 1 in steps of 1/4, both sides of the thresholds 0.5;
# - exprel at large |x|, expexp, cexpexp, cloglog, loglog and log1psq inside their domains.
#
# Tool: Python mpmath 1.3.0, 50 decimal digits, derivatives by mpmath.diff, run in the /home/ts/penv environment:
#
#   cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_autodiff_core_mpmath_reference.py
#
# Output: test/resources/calculus/autodiff_core_mpmath_reference.edn
# {<id> [[x d1 d2] ..]}, the ids are the names of fastmath.calculus.autodiff.

import os
import mpmath as mp

mp.mp.dps = 50


def sinc(x):
    return mp.sinc(mp.pi * x) if x != 0 else mp.mpf(1)


def exprel(x):
    return mp.expm1(x) / x if x != 0 else mp.mpf(1)


def around_zero():
    xs = [0.0]
    for i in range(0, 45):
        x = float(mp.mpf(10) ** (mp.mpf(-10) + mp.mpf(i) / 4))
        xs += [x, -x]
    return xs


general = [-30.0, -3.7, -1.2, -0.3, 0.2, 0.9, 2.5, 30.0]
probability = [1.0e-10, 0.01, 0.2, 0.5, 0.7, 0.99, 0.999999]

functions = [
    ("sinc", sinc, around_zero()),
    ("exprel", exprel, around_zero() + [-700.0, -30.0, 30.0, 700.0]),
    ("expexp", lambda x: mp.exp(-mp.exp(-x)), general),
    ("cexpexp", lambda x: -mp.expm1(-mp.exp(x)), general),
    ("cloglog", lambda x: mp.log(-mp.log1p(-x)), probability),
    ("loglog", lambda x: -mp.log(-mp.log(x)), probability),
    ("log1psq", lambda x: mp.log1p(x * x), general + [1.0e100, -1.0e100])]


def row(f, x):
    return "[%r %r %r]" % (x, float(mp.diff(f, mp.mpf(x), 1)), float(mp.diff(f, mp.mpf(x), 2)))


out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                   "test", "resources", "calculus", "autodiff_core_mpmath_reference.edn")
with open(out, "w") as f:
    f.write("{")
    f.write("\n ".join(":%s [\n  %s]" % (name, "\n  ".join(row(g, x) for x in xs)) for name, g, xs in functions))
    f.write("}\n")
