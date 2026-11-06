# Reference values for the derivative rules of the Bessel, spherical Bessel and Airy functions of
# fastmath.calculus.autodiff (phase 5c): first and second derivatives in x, the order constant. Also the order -1/2 of
# the spherical functions (where the recurrence divides by 2n+1 = 0).
#
# mpmath for all of them: the ForwardDiff derivatives of SpecialFunctions.jl lose digits at some points
# (besseli at order -1.5, x = 20: 2.5e-14; airyai at -7: 1.7e-14; sphericalbesselj(0, 0.3), which differentiates
# sin(x)/x: 8e-15), where the autodiff values agree with mpmath to all digits.
#
# Tool: Python mpmath 1.3.0, 50 decimal digits, derivatives by mpmath.diff, run in the /home/ts/penv environment:
#
#   cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_autodiff_bessel_mpmath_reference.py
#
# Output: test/resources/calculus/autodiff_bessel_mpmath_reference.edn
# {<id> [[order x d1 d2] ..]}, the ids are the names of fastmath.calculus.autodiff (the Airy functions with order 0.0,
# unused); for the half odd order K the order is the odd integer (the order of K is order/2). The fixed-order functions
# (bessel-J0, spherical-bessel-j1, ...) are tested against the rows of their order. Definitions (DLMF 10.47): j_n = sqrt(pi/(2x)) J_(n+1/2),
# y_n = sqrt(pi/(2x)) Y_(n+1/2), i1_n = sqrt(pi/(2x)) I_(n+1/2), i2_n = sqrt(pi/(2x)) I_(-n-1/2),
# k_n = sqrt(pi/(2x)) K_(n+1/2).

import os
import mpmath as mp

mp.mp.dps = 50


def edn(v):
    v = float(v)
    if v != v:
        return "##NaN"
    if v in (float("inf"), float("-inf")):
        return "##Inf" if v > 0 else "##-Inf"
    return repr(v)


def sph(f, shift):
    # at a negative x (integer orders, j only) the reflection j_n(-x) = (-1)^n j_n(x): the product form at a negative
    # argument is on the other branch of sqrt (the opposite sign)
    def at(n, x):
        if x < 0:
            return (-1) ** int(n) * mp.sqrt(mp.pi / (-2 * x)) * f(shift(n), -x)
        return mp.sqrt(mp.pi / (2 * x)) * f(shift(n), x)
    return lambda n: (lambda x: at(n, x))


up = lambda n: n + mp.mpf(1) / 2
down = lambda n: -n - mp.mpf(1) / 2

positive = [0.3, 1.0, 5.0, 20.0]
orders = [0.0, 1.0, 2.0, 2.5, -1.5, 7.0]
with_negative = positive + [-0.7, -4.0]


def bessel(f):
    return lambda n: (lambda x: f(n, x))


def airy(f):
    return lambda n: f


# orders next to -1/2 (where the recurrence through both neighbouring orders cancels), for every spherical family except y: the
# value of fastmath.special/spherical-bessel-y there has the documented error about eps/d (d the distance from -1/2), the
# derivatives follow it (1.6e-4 at d = 1e-12); k: 6.4e-14 (orders -0.499, -0.501, x <= 1), 1.4e-12 between the grid points
# (1e-5 < d < 1e-3, x <= 1.5: k at the neighbouring order, bessel-K next to an integer order)
near_minus_half = [-0.5 - 1.0e-12, -0.5 + 1.0e-12, -0.500001, -0.499999, -0.501, -0.499, -0.53, -0.47, -0.2]

functions = [
    ("bessel-J", bessel(mp.besselj), orders, with_negative),
    ("bessel-Y", bessel(mp.bessely), orders, positive),
    ("bessel-I", bessel(mp.besseli), orders, with_negative),
    ("bessel-K", bessel(mp.besselk), orders, positive),
    ("airy-Ai", airy(mp.airyai), [0.0], [-7.0, -1.2, 0.0, 0.7, 3.0, 12.0]),
    ("airy-Ai'", airy(lambda x: mp.airyai(x, derivative=1)), [0.0], [-7.0, -1.2, 0.0, 0.7, 3.0, 12.0]),
    ("airy-Bi", airy(mp.airybi), [0.0], [-7.0, -1.2, 0.0, 0.7, 3.0, 12.0]),
    ("airy-Bi'", airy(lambda x: mp.airybi(x, derivative=1)), [0.0], [-7.0, -1.2, 0.0, 0.7, 3.0, 12.0]),
    ("spherical-bessel-1-i", sph(mp.besseli, up), [0.0, 1.0, 2.0, 3.5, -0.5] + near_minus_half, positive),
    ("spherical-bessel-2-i", sph(mp.besseli, down), [0.0, 1.0, 2.0, 3.5, -0.5] + near_minus_half, positive),
    ("spherical-bessel-k", sph(mp.besselk, up), [0.0, 1.0, 2.0, 3.5, -0.5] + near_minus_half, positive),
    ("spherical-bessel-j", sph(mp.besselj, up), [0.0, 1.0, 2.0, 3.5, -0.5] + near_minus_half, positive + [-0.7, -4.0]),
    ("spherical-bessel-y", sph(mp.bessely, up), [0.0, 1.0, 2.0, 3.5, -0.5], positive),
    ("bessel-K-half-odd", lambda n: (lambda x: mp.besselk(mp.mpf(n) / 2, x)), [1, 3, 5, -3], positive),
    ("bessel-K-half-odd-scaled", lambda n: (lambda x: mp.exp(x) * mp.besselk(mp.mpf(n) / 2, x)), [1, 3, 5, -3], positive)]


def row(f, n, x):
    # negative x: integer orders only (the value is complex otherwise)
    if x < 0 and n != int(n):
        return None
    g = f(n)
    return "[%s %s %s %s]" % (edn(n), edn(x), edn(mp.diff(g, mp.mpf(x), 1)), edn(mp.diff(g, mp.mpf(x), 2)))


out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                   "test", "resources", "calculus", "autodiff_bessel_mpmath_reference.edn")
with open(out, "w") as f:
    f.write("{")
    f.write("\n ".join(":%s [\n  %s]" % (name, "\n  ".join(r for r in (row(g, n, x) for n in ns for x in xs) if r))
                       for name, g, ns, xs in functions))
    f.write("}\n")
