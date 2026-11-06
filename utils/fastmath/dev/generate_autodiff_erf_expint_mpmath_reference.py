# Reference values for the derivative rules of fastmath.calculus.autodiff (phase 5b) which have no rule in Julia:
# Lambert W (both branches), E0, Ein, Cin, si, li and Li, and Ei (the ForwardDiff derivative of expinti loses
# digits); first and second derivatives of the exact functions.
#
# Tool: Python mpmath 1.3.0, 50 decimal digits, derivatives by mpmath.diff (Lambert W: the closed forms
# W' = W/(x(1+W)), W'' = -W^2 (W+2)/(x^2 (1+W)^3)), run in the /home/ts/penv environment:
#
#   cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_autodiff_erf_expint_mpmath_reference.py
#
# Output: test/resources/calculus/autodiff_erf_expint_mpmath_reference.edn
# {<id> [[x d1 d2] ..]                the ids are the names of fastmath.calculus.autodiff
#  :En-real-branch [[n x d1 d2] ..]}   En of a non-positive integer order at a negative x (real there), derivatives in x
# Points inside the domains, away from the singularities (Lambert W branch point -1/e, li at 1); Ein also below -2.15
# (-3, -10).

import os
import mpmath as mp

mp.mp.dps = 50


def ein(x):
    # Ein(x) = sum_{k>=1} (-1)^(k+1) x^k / (k k!), entire
    return mp.nsum(lambda k: (-1) ** (k + 1) * x ** k / (k * mp.factorial(k)), [1, mp.inf])


def cin(x):
    return mp.quad(lambda t: (1 - mp.cos(t)) / t if t != 0 else mp.mpf(0), [0, x])


def lambert_w0(x):
    return mp.lambertw(x, 0).real


def lambert_wm1(x):
    return mp.lambertw(x, -1).real


lambert = {lambert_w0: 0, lambert_wm1: -1}

functions = [
    ("lambert-W", lambert_w0, [-0.3, -0.1, 0.0, 0.5, 1.0, 2.5, 2.8, 10.0, 1.0e6, 1.0e300]),
    ("lambert-W-1", lambert_wm1, [-0.36, -0.3, -0.1, -1.0e-3, -1.0e-10, -1.0e-300]),
    ("E0", lambda x: mp.exp(-x) / x, [-3.7, -0.3, 0.2, 2.5, 30.0]),
    # -3 and -10 appended (not sorted): the earlier rows keep their places in the output
    ("Ein", ein, [-2.1, -1.2, -0.3, 0.0, 0.2, 0.9, 2.5, 7.0, -3.0, -10.0]),
    ("Cin", cin, [-7.0, -1.2, 0.0, 1.0e-3, 0.2, 2.5, 30.0]),
    ("si", lambda x: mp.si(x) - mp.pi / 2, [-3.7, -0.3, 0.2, 2.5, 30.0]),
    ("Ei", mp.ei, [-30.0, -2.5, -0.3, 0.2, 0.9, 2.5, 7.0, 30.0, 300.0, 709.5, 710.0, 716.0]),
    ("li", mp.li, [0.1, 0.5, 1.5, 2.0, 10.0, 1.0e6]),
    ("Li", lambda x: mp.li(x) - mp.li(2), [0.1, 0.5, 1.5, 2.0, 10.0, 1.0e6])]


def edn(v):
    v = float(v)
    if v != v:
        return "##NaN"
    if v in (float("inf"), float("-inf")):
        return "##Inf" if v > 0 else "##-Inf"
    return repr(v)


def lambert_row(k, x):
    # closed forms with W at 50 digits: mpmath.diff steps across 0 next to it (x = -1e-300); 1 and -2 at x = 0
    if x == 0:
        return "[%r %r %r]" % (x, 1.0, -2.0)
    x = mp.mpf(x)
    w = mp.lambertw(x, k).real
    return "[%s %s %s]" % (edn(x), edn(w / (x * (1 + w))), edn(-w * w * (w + 2) / (x * x * (1 + w) ** 3)))


def row(f, x):
    if f in lambert:
        return lambert_row(lambert[f], x)
    return "[%s %s %s]" % (edn(x), edn(mp.diff(f, mp.mpf(x), 1)), edn(mp.diff(f, mp.mpf(x), 2)))


out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                   "test", "resources", "calculus", "autodiff_erf_expint_mpmath_reference.edn")
def en_row(n, x):
    g = lambda y: mp.re(mp.expint(n, y))
    return "[%s %s %s %s]" % (edn(n), edn(x), edn(mp.diff(g, mp.mpf(x), 1)), edn(mp.diff(g, mp.mpf(x), 2)))


with open(out, "w") as f:
    f.write("{")
    f.write("\n ".join(":%s [\n  %s]" % (name, "\n  ".join(row(g, x) for x in xs)) for name, g, xs in functions))
    f.write("\n :En-real-branch [\n  %s]" % "\n  ".join(en_row(n, x) for n in (0, -1, -2, -3) for x in (-2.5, -1.0, -0.3)))
    f.write("}\n")
