# Reference values for the derivative rules of the Jacobi elliptic functions, their inverses, the amplitude and the
# incomplete elliptic integrals of fastmath.calculus.autodiff (phase 5d): first and second derivatives in the
# differentiated argument (u, x or phi), the parameter m = k^2 (and the characteristic n) constant.
#
# Tool: Python mpmath 1.3.0, 50 decimal digits, derivatives by mpmath.diff, run in the /home/ts/penv environment:
#
#   cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_autodiff_elliptic_mpmath_reference.py
#
# Output: test/resources/calculus/autodiff_elliptic_mpmath_reference.edn
# {<id> [[a m d1 d2] ..]} (a = u, x or phi), :elliptic-PI [[n phi m d1 d2] ..]; the ids are the names of
# fastmath.calculus.autodiff.
# - Jacobi functions: mpmath.ellipfun (parameter m, also m = 1.5 > 1), the quotients of sn, cn, dn;
# - jacobi-am: the root phi of ellipf(phi, m) = u;
# - inverses: at x = f(u0) with u0 = 0.3 and 0.7 quarter periods (K(m); K(1/m)/sqrt(m) for m > 1), in the principal range of
#   fastmath.special), the derivatives of the inverse 1/f'(u0) and -f''(u0)/f'(u0)^3;
# - integrals: ellipf, ellipe, ellippi, and D = (F - E)/m.

import os
import mpmath as mp

mp.mp.dps = 50


def edn(v):
    if isinstance(v, mp.mpc):
        # mpmath returns complex numbers with a zero imaginary part for some real arguments (negative m)
        assert abs(v.imag) <= mp.mpf(10) ** -30 * max(1, abs(v.real)), v
        v = v.real
    v = float(v)
    if v != v:
        return "##NaN"
    if v in (float("inf"), float("-inf")):
        return "##Inf" if v > 0 else "##-Inf"
    return repr(v)


def sn(u, m): return mp.ellipfun("sn", u, m=m)
def cn(u, m): return mp.ellipfun("cn", u, m=m)
def dn(u, m): return mp.ellipfun("dn", u, m=m)


jacobi = {
    "sn": sn, "cn": cn, "dn": dn,
    "sc": lambda u, m: sn(u, m) / cn(u, m), "sd": lambda u, m: sn(u, m) / dn(u, m),
    "cs": lambda u, m: cn(u, m) / sn(u, m), "cd": lambda u, m: cn(u, m) / dn(u, m),
    "ds": lambda u, m: dn(u, m) / sn(u, m), "dc": lambda u, m: dn(u, m) / cn(u, m),
    "ns": lambda u, m: 1 / sn(u, m), "nc": lambda u, m: 1 / cn(u, m), "nd": lambda u, m: 1 / dn(u, m)}

ms = [0.3, 0.8, -0.5]
us = [-0.9, 0.3, 1.0, 3.0]
phis = [-1.2, 0.3, 1.0, 2.5]


def am(u, m):
    return mp.findroot(lambda p: mp.ellipf(p, m) - u, u)


def quarter_period(m):
    # K(m) for m <= 1; for m > 1 the real quarter period K(1/m)/sqrt(m) (sn(u|m) = sn(sqrt(m) u|1/m)/sqrt(m))
    return mp.ellipk(m) if m <= 1 else mp.ellipk(1 / mp.mpf(m)) / mp.sqrt(m)


# m = 1.5 (> 1, where fastmath.special uses the real analytic continuation): u inside the first quarter period only
# (the quotients have poles where dn = 0 beyond it)
ms_jacobi = [(m, us) for m in ms] + [(1.5, [-0.9, 0.3, 1.0])]

rows = {}
for name, f in jacobi.items():
    rows["jacobi-" + name] = [(u, m, mp.diff(lambda t: f(t, m), mp.mpf(u), 1), mp.diff(lambda t: f(t, m), mp.mpf(u), 2))
                              for m, uss in ms_jacobi for u in uss]
    inv = []
    for m in ms + [1.5]:
        for frac in (0.3, 0.7):
            u0 = frac * quarter_period(m)
            d1 = mp.diff(lambda t: f(t, m), u0, 1)
            d2 = mp.diff(lambda t: f(t, m), u0, 2)
            inv.append((f(u0, m), m, 1 / d1, -d2 / d1 ** 3))
    rows["jacobi-a" + name] = inv
# jacobi-am for m >= 0 only: fastmath.special/jacobi-am is NaN for a negative m
rows["jacobi-am"] = [(u, m, mp.diff(lambda t: am(t, m), mp.mpf(u), 1), mp.diff(lambda t: am(t, m), mp.mpf(u), 2))
                     for m in ms if m >= 0 for u in us]
for name, f in [("elliptic-F", mp.ellipf), ("elliptic-E", mp.ellipe),
                ("elliptic-D", lambda p, m: (mp.ellipf(p, m) - mp.ellipe(p, m)) / m)]:
    rows[name] = [(p, m, mp.diff(lambda t: f(t, m), mp.mpf(p), 1), mp.diff(lambda t: f(t, m), mp.mpf(p), 2))
                  for m in ms for p in phis]
pi_rows = [(n, p, m, mp.diff(lambda t: mp.ellippi(n, t, m), mp.mpf(p), 1), mp.diff(lambda t: mp.ellippi(n, t, m), mp.mpf(p), 2))
           for n in (0.4, -0.7) for m in ms for p in phis]

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                   "test", "resources", "calculus", "autodiff_elliptic_mpmath_reference.edn")
with open(out, "w") as fo:
    fo.write("{")
    fo.write("\n ".join(":%s [\n  %s]" % (name, "\n  ".join("[%s]" % " ".join(edn(v) for v in r) for r in rs))
                        for name, rs in rows.items()))
    fo.write("\n :elliptic-PI [\n  %s]" % "\n  ".join("[%s]" % " ".join(edn(v) for v in r) for r in pi_rows))
    fo.write("}\n")
