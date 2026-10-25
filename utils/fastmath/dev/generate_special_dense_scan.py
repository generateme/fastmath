"""Dense scan grids for the branch-boundary check of `fastmath.special` (audit WP 9.2).

Writes `/tmp/special_dense_scan.edn`, a map `{function-name [[x reference]...]}` with `mpmath` (30 digits)
values on fine uniform grids (a branch boundary of the piecewise approximations falls between two grid
points) and on log-spaced grids for large arguments.  Not stored in the repository (about 10 MB);
`scan_special_dense.clj` compares the functions with it and prints the worst error per interval.

Run with:
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_special_dense_scan.py
"""
import math

import mpmath as mp

mp.mp.dps = 30
OUT = "/tmp/special_dense_scan.edn"


def frange(a, b, step):
    n = int(round((b - a) / step))
    return [a + i * step for i in range(n + 1)]


def logrange(a, b, n):
    la, lb = math.log10(a), math.log10(b)
    return [10 ** (la + (lb - la) * i / (n - 1)) for i in range(n)]


FUNCS = {
    'bessel-J0': (lambda x: mp.besselj(0, x), frange(0.013, 200, 0.013) + logrange(200, 1e6, 1500)),
    'bessel-J1': (lambda x: mp.besselj(1, x), frange(0.013, 200, 0.013) + logrange(200, 1e6, 1500)),
    'bessel-Y0': (lambda x: mp.bessely(0, x), frange(0.013, 200, 0.013) + logrange(200, 1e6, 1500)),
    'bessel-Y1': (lambda x: mp.bessely(1, x), frange(0.013, 200, 0.013) + logrange(200, 1e6, 1500)),
    'bessel-I0': (lambda x: mp.besseli(0, x), frange(0.05, 700, 0.05)),
    'bessel-I1': (lambda x: mp.besseli(1, x), frange(0.05, 700, 0.05)),
    'bessel-K0': (lambda x: mp.besselk(0, x), frange(0.05, 700, 0.05)),
    'bessel-K1': (lambda x: mp.besselk(1, x), frange(0.05, 700, 0.05)),
    'digamma': (mp.digamma, [x for x in frange(-20.0127, 20, 0.0127) if abs(x - round(x)) > 1e-3] + logrange(20, 1e8, 1500)),
    'trigamma': (lambda x: mp.psi(1, x), [x for x in frange(-20.0127, 20, 0.0127) if abs(x - round(x)) > 1e-3] + logrange(20, 1e8, 1500)),
    'zeta': (mp.zeta, [x for x in frange(-60.0013, 60, 0.0107) if abs(x - 1) > 1e-3]),
    'airy-Ai': (mp.airyai, frange(-100.0013, 100, 0.0213)),
    "airy-Ai'": (lambda x: mp.airyai(x, derivative=1), frange(-100.0013, 100, 0.0213)),
    'airy-Bi': (mp.airybi, frange(-100.0013, 100, 0.0213)),
    "airy-Bi'": (lambda x: mp.airybi(x, derivative=1), frange(-100.0013, 100, 0.0213)),
    'Si': (mp.si, frange(0.0107, 100, 0.0107) + logrange(100, 1e5, 1000)),
    'Ci': (mp.ci, frange(0.0107, 100, 0.0107) + logrange(100, 1e5, 1000)),
    'E1': (mp.e1, frange(0.0107, 100, 0.0107)),
    'Ei': (mp.ei, frange(0.0107, 100, 0.0107)),
    'si': (lambda x: mp.si(x) - mp.pi / 2, frange(0.0107, 100, 0.0107)),
    'Cin': (lambda x: mp.euler + mp.log(x) - mp.ci(x), frange(0.0107, 100, 0.0107)),
    'E0': (lambda x: mp.exp(-x) / x, frange(0.0107, 100, 0.0107)),
    'Ein': (lambda x: mp.e1(x) + mp.log(x) + mp.euler, frange(0.0107, 100, 0.0107)),
    'li': (mp.li, [x for x in frange(0.0107, 100, 0.0107) if abs(x - 1) > 1e-2]),
    'Li': (lambda x: mp.li(x) - mp.li(2), frange(0.0107, 100, 0.0107)),
}


def edn(v):
    return repr(float(v))


with open(OUT, 'w') as f:
    f.write('{\n')
    for name, (fn, xs) in FUNCS.items():
        f.write(' "%s" [\n' % name)
        for x in xs:
            try:
                v = fn(mp.mpf(x))
                if isinstance(v, mp.mpc) or not mp.isfinite(v):
                    continue
                f.write('  [%s %s]\n' % (edn(x), edn(v)))
            except Exception:
                pass
        f.write(' ]\n')
    f.write('}\n')
print('written', OUT)


# ----------------------------------------------------------------------------- real orders

OUT2 = "/tmp/special_dense_scan_orders.edn"
ORDERS = [0.3, 1.5, 2.7, 5.2, 10.5, 17.0, 25.5, 40.0, 75.5, 120.0]
XS = logrange(0.01, 1000, 220)
FUNCS2 = {
    'bessel-J': mp.besselj,
    'bessel-Y': mp.bessely,
    'bessel-I': mp.besseli,
    'bessel-K': mp.besselk,
}
with open(OUT2, 'w') as f:
    f.write('{\n')
    for name, fn in FUNCS2.items():
        f.write(' "%s" [\n' % name)
        for v in ORDERS:
            for x in XS:
                try:
                    r = fn(mp.mpf(v), mp.mpf(x))
                    if isinstance(r, mp.mpc) or not mp.isfinite(r) or abs(r) > 1e300 or (r != 0 and abs(r) < 1e-300):
                        continue
                    f.write('  [%s %s %s]\n' % (edn(v), edn(x), edn(r)))
                except Exception:
                    pass
        f.write(' ]\n')
    f.write('}\n')
print('written', OUT2)
