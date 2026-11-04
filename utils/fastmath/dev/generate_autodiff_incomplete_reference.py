# Reference values for the derivatives in x of the incomplete gamma and beta functions of
# fastmath.calculus.autodiff (ForwardDiff.jl has no rule for them: a MethodError for gamma_inc in a and beta_inc,
# NaN for the upper gamma(a, x)).
#
# Tool: Python mpmath 1.3.0, 40 significant digits, numerical differentiation mpmath.diff of the functions
#       themselves (independent of the closed forms used by fastmath), rounded to doubles; run in /home/ts/penv:
#
#   uv run python utils/fastmath/dev/generate_autodiff_incomplete_reference.py
#
# Output: test/resources/calculus/autodiff_incomplete_reference.edn
# {<id> [[shape.. x d1 d2] ..]}, the ids are the names of fastmath.calculus.autodiff, the shape parameters in the
# order of the arguments of fastmath (the beta functions take x first, its row is [a b x d1 d2]).
# The shape parameters include negative values, where gamma(a) or beta(a, b) is negative.

import os
import mpmath as mp

mp.mp.dps = 40

def lower(s, x): return mp.gammainc(s, 0, x)
def upper(s, x): return mp.gammainc(s, x, mp.inf)
def gp(a, x): return mp.gammainc(a, 0, x, regularized=True)
def gq(a, x): return mp.gammainc(a, x, mp.inf, regularized=True)
def ib(a, b, x): return mp.betainc(a, b, 0, x)
def rb(a, b, x): return mp.betainc(a, b, 0, x, regularized=True)

gamma_points = [(s, x) for s in ("0.5", "1", "2.5", "10", "-0.5", "-1.5") for x in ("0.1", "1", "3.5")]
gamma_points += [("200", "200"), ("200", "180")]
beta_points = [(a, b, x) for (a, b) in (("0.5", "0.5"), ("2", "3"), ("1", "1"), ("7.5", "1.5"), ("-0.5", "2"), ("2", "-0.5"))
               for x in ("0.05", "0.3", "0.9")]
beta_points += [("300", "200", "0.6")]

def num(v):
    d = float(v)
    if d == float("inf"):
        return "##Inf"
    if d == float("-inf"):
        return "##-Inf"
    return repr(d)

def rows(f, points):
    out = []
    for p in points:
        args = [mp.mpf(a) for a in p]
        g = lambda t: f(*args[:-1], t)
        d1 = mp.diff(g, args[-1], 1)
        d2 = mp.diff(g, args[-1], 2)
        out.append("[" + " ".join([num(a) for a in args] + [num(d1), num(d2)]) + "]")
    return out

functions = [("lower-incomplete-gamma", lower, gamma_points),
             ("upper-incomplete-gamma", upper, gamma_points),
             ("regularized-gamma-p", gp, gamma_points),
             ("regularized-gamma-q", gq, gamma_points),
             ("incomplete-beta", ib, beta_points),
             ("regularized-beta", rb, beta_points)]

path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "test", "resources", "calculus",
                    "autodiff_incomplete_reference.edn")
with open(path, "w") as io:
    io.write("{\n")
    for name, f, points in functions:
        io.write(" :" + name + " [")
        for r in rows(f, points):
            io.write("\n  " + r)
        io.write("]\n")
    io.write("}\n")
