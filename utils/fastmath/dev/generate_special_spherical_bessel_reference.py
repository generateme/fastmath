# Reference values for the spherical Bessel fixes of fastmath.special (upstream group 4):
#
# - :zero        the values at x = +0.0 and -0.0 of spherical-bessel-j, -y, -1-i, -2-i, -k: the one-sided limits (the
#                limit from the left for integer orders at -0.0, x < 0 is in the domain there; the limit from the right
#                otherwise). A limit is classified by f(+-1e-50) / f(+-1e-25): growth gives +-Inf, decay a signed zero,
#                otherwise the finite value. The modified functions are summed from their power series with 1/Gamma
#                (0 at the poles), mpmath.besseli returns 0 at a negative integer order and a tiny x;
# - :K           bessel-K next to integer orders, m + mu, mu from 1e-5 to 0.3 (x <= 1.5: the Temme series and the power
#                series);
# - :I           bessel-I at negative orders next to a negative integer, -(m + d);
# - :Y           bessel-Y next to integer orders, m + d (x < 7 loses about eps/d; documented, not fixed);
# - :spherical   the five spherical functions next to half-integer orders, c + d.
#
# Tool: Python mpmath 1.3.0, 50 decimal digits. Run in the /home/ts/penv environment:
#
#   cd /home/ts/clojure/fastmath && /home/ts/penv/.venv/bin/python utils/fastmath/dev/generate_special_spherical_bessel_reference.py
#
# Output: test/resources/special/spherical_bessel_reference.edn

import os
import mpmath as mp

mp.mp.dps = 50


def fmt(x):
    x = float(x)
    if x == float('inf'):
        return '##Inf'
    if x == float('-inf'):
        return '##-Inf'
    return repr(x)


def vec(xs):
    return '[' + ' '.join(fmt(x) for x in xs) + ']'


def strvec(xs):
    return '[' + ' '.join('"%s"' % x for x in xs) + ']'


def iser(v, x):
    return mp.fsum((x / 2) ** (2 * k + v) * mp.rgamma(k + 1) * mp.rgamma(k + v + 1) for k in range(40))


def sph(nm, n, x):
    x = mp.mpf(x)
    n = mp.mpf(n)
    c = mp.sqrt(mp.pi / (2 * x))
    if nm == 'j':
        return c * mp.besselj(n + 0.5, x)
    if nm == 'y':
        return c * mp.bessely(n + 0.5, x)
    if nm == 'i1':
        return c * iser(n + 0.5, x)
    if nm == 'i2':
        return c * iser(-n - 0.5, x)
    return c * mp.besselk(n + 0.5, x)


PARITY = {'j': 0, 'i1': 0, 'y': 1, 'i2': 1}


def right_limit(nm, n):
    a = sph(nm, n, mp.mpf('1e-25'))
    b = sph(nm, n, mp.mpf('1e-50'))
    r = abs(b / a)
    if r > 10:
        return mp.inf if b > 0 else -mp.inf
    if r < 0.1:
        return 0.0 if b > 0 else -0.0
    return b


def at_zero(nm, n, negative_zero):
    v = right_limit(nm, n)
    if negative_zero and nm != 'k' and n == int(n) and (int(n) + PARITY[nm]) % 2 == 1:
        return -v
    return v


# x = 0
zf, zn, zx, zr = [], [], [], []
for n in [0.0, 1.0, 2.0, 3.0, 4.0, 0.5, 1.5, -0.3, -0.5, -0.7, -1.0, -1.5, -2.0, -2.5, -3.0, -3.3, -4.0, -5.5, -6.0,
          7.3]:
    for nm in ['j', 'y', 'i1', 'i2', 'k']:
        for neg in (False, True):
            zf.append(nm); zn.append(n); zx.append(-0.0 if neg else 0.0); zr.append(at_zero(nm, n, neg))

# K next to integer orders
kn, kx, kr = [], [], []
for m in [0, 1, 2, 5, 10, 20]:
    for mu in [1e-5, 5e-5, 1e-4, 1.4e-4, 1.6e-4, 3e-4, 1e-3, 1e-2, 0.1, 0.3]:
        for s in (1, -1):
            v = m + s * mu
            if v <= 0:
                continue
            for x in [0.01, 0.1, 0.5, 1.0, 1.4, 1.5]:
                kn.append(v); kx.append(x); kr.append(mp.besselk(mp.mpf(v), mp.mpf(x)))

# I at negative orders next to a negative integer
i_n, i_x, i_r = [], [], []
for m in [1, 2, 3, 5]:
    for d in [1e-12, 1e-9, 1e-6, 1e-3]:
        for s in (1, -1):
            v = -(m + s * d)
            for x in [0.01, 0.1, 1.0, 10.0]:
                i_n.append(v); i_x.append(x); i_r.append(mp.besseli(mp.mpf(v), mp.mpf(x)))

# Y next to integer orders
yn, yx, yr = [], [], []
for m in [0, 1, 2, 3]:
    for d in [1e-12, 1e-9, 1e-6, 1e-3]:
        for s in (1, -1):
            v = m + s * d
            for x in [0.1, 1.0, 3.0, 6.9, 7.5, 10.0]:
                yn.append(v); yx.append(x); yr.append(mp.bessely(mp.mpf(v), mp.mpf(x)))

# spherical next to half-integer orders
sf, sn, sx, sr = [], [], [], []
for c in [-1.5, -0.5, 0.5, 1.5, 2.5]:
    for d in [1e-12, 1e-9, 1e-6, 1e-3]:
        for s in (1, -1):
            n = c + s * d
            for x in [0.1, 1.0, 10.0]:
                for nm in ['j', 'y', 'i1', 'i2', 'k']:
                    sf.append(nm); sn.append(n); sx.append(x)
                    sr.append(sph(nm, n, x) if nm not in ('i1', 'i2') else
                              mp.sqrt(mp.pi / (2 * mp.mpf(x))) * mp.besseli((n + 0.5) if nm == 'i1' else (-n - 0.5),
                                                                            mp.mpf(x)))

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..',
                   'test', 'resources', 'special', 'spherical_bessel_reference.edn')
with open(out, 'w') as f:
    f.write('{:zero {:family %s :order %s :x %s :ref %s}\n' % (strvec(zf), vec(zn), vec(zx), vec(zr)))
    f.write(' :K {:order %s :x %s :ref %s}\n' % (vec(kn), vec(kx), vec(kr)))
    f.write(' :I {:order %s :x %s :ref %s}\n' % (vec(i_n), vec(i_x), vec(i_r)))
    f.write(' :Y {:order %s :x %s :ref %s}\n' % (vec(yn), vec(yx), vec(yr)))
    f.write(' :spherical {:family %s :order %s :x %s :ref %s}}\n' % (strvec(sf), vec(sn), vec(sx), vec(sr)))
print('zero', len(zf), 'K', len(kn), 'I', len(i_n), 'Y', len(yn), 'spherical', len(sf))
