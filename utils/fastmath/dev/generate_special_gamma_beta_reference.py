# Reference values for the gamma and beta fixes of fastmath.special (upstream group 2):
#
# - :gamma-small  P(a,x), Q(a,x) and the upper incomplete gamma Gamma(a,x) for small shapes a (1e-300 to 0.9) and
#                 x from 1e-300 to 700 (both sides of the x = 1.5 switch of the small-shape series);
# - :gamma-large  P(a,x), Q(a,x) for large shapes a (1e4 to 1e300; 1.5e4 to 5e4 around the switch to Temme's expansion), x = a + r sqrt(a) and x = lambda a; computed by
#                 Gauss-Legendre quadrature of the density in v = (t/a - 1) sqrt(a) (mpmath.gammainc does not
#                 converge there, R pgamma is off by 3e-9 at a = 1e20); only values >= 1e-300;
# - :beta-poles   B(p, q) where p + q is a non-positive integer and p a positive integer: (-1)^p B(p, 1 - p - q)
#                 (the pole of Gamma(q) cancels the pole of Gamma(p + q)); doubles above 2^53 are even;
# - :beta-tiny    I_x(a, b) and B(x; a, b) with a shape below 1e-3 (mpmath.betainc).
#
# Tool: Python mpmath 1.3.0, 60 decimal digits (60 + log10(a) for the large-shape quadrature). Run in the /home/ts/penv environment:
#
#   cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_special_gamma_beta_reference.py
#
# Output: test/resources/special/gamma_beta_reference.edn

import os
import mpmath as mp

mp.mp.dps = 60


def fmt(x):
    x = float(x)
    if x == float('inf'):
        return '##Inf'
    if x == float('-inf'):
        return '##-Inf'
    return repr(x)


def vec(xs):
    return '[' + ' '.join(fmt(x) for x in xs) + ']'


def block(name, cols):
    return ':%s {%s}' % (name, ' '.join(':%s %s' % (k, vec(v)) for k, v in cols))


# gamma, small shapes
gs = {'a': [], 'x': [], 'p': [], 'q': [], 'upper': []}
for a in [1e-300, 1e-200, 1e-100, 1e-50, 1e-20, 1e-10, 1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.0999, 0.1, 0.3, 0.9]:
    for x in [1e-300, 1e-100, 1e-20, 1e-5, 0.01, 0.1, 0.5, 1.0, 1.4999999999999998, 1.5, 1.5000000000000002,
              2.0, 5.0, 20.0, 100.0, 700.0]:
        A, X = mp.mpf(a), mp.mpf(x)
        q = mp.gammainc(A, X, regularized=True)
        p = mp.gammainc(A, 0, X, regularized=True)
        u = mp.gammainc(A, X)
        gs['a'].append(a); gs['x'].append(x); gs['p'].append(p); gs['q'].append(q); gs['upper'].append(u)


# gamma, large shapes, by quadrature
def tail(a, x, lower):
    # the exponent cancels terms of the size a log(a): the precision grows with log10(a)
    mp.mp.dps = 60 + int(mp.log10(a))
    a, x = mp.mpf(a), mp.mpf(x)
    sa = mp.sqrt(a)
    v0 = (x / a - 1) * sa
    c = a * mp.log(a) - mp.loggamma(a)
    f = lambda v: mp.exp(c + (a - 1) * mp.log1p(v / sa) - a * (1 + v / sa)) / sa
    slope = abs((a - 1) / (sa + v0) - sa)
    w = 60 / max(1, slope)
    if lower:
        lo = max(-sa, v0 - w)
        return mp.quad(f, mp.linspace(lo, v0, 201), method='gauss-legendre')
    return mp.quad(f, mp.linspace(v0, v0 + w, 201), method='gauss-legendre')


gl = {'a': [], 'x': [], 'p': [], 'q': []}
for a in [1e4, 1.5e4, 3e4, 5e4, 1e5, 1e6, 1e8, 1e10, 1e15, 1e20, 1e30, 1e100, 1e300]:
    xs = [a + r * a ** 0.5 for r in (-30.0, -5.0, -1.0, 0.0, 1.0, 5.0, 30.0)]
    xs += [a * l for l in (0.9, 0.99, 0.999, 1.001, 1.01, 1.1)]
    for x in xs:
        p, q = tail(a, x, True), tail(a, x, False)
        if p < q:
            q = 1 - p
        else:
            p = 1 - q
        if min(p, q) >= mp.mpf('1e-300'):
            gl['a'].append(a); gl['x'].append(x); gl['p'].append(p); gl['q'].append(q)

mp.mp.dps = 60

# beta where p + q is a non-positive integer
bp = {'p': [], 'q': [], 'ref': []}
for p in [1.0, 2.0, 3.0, 10.0, 1e3, 1e6, 1e9, 1e12, 1e15, 2.0 ** 52, 2.0 ** 53 + 2, 1e20, 1e100, 1e300]:
    for k in (0, 3):
        q = -(p + k)
        if q == -p - k and q != -p - k + 1:
            P = mp.mpf(p)
            even = (p >= 2 ** 53) or int(p) % 2 == 0
            ref = (1 if even else -1) * mp.beta(P, 1 - P - mp.mpf(q))
            bp['p'].append(p); bp['q'].append(q); bp['ref'].append(ref)

# incomplete beta with a tiny shape
bt = {'x': [], 'a': [], 'b': [], 'reg': [], 'inc': []}
for a in [1e-300, 1e-100, 1e-20, 1e-8, 1e-4]:
    for b in [1e-300, 1e-20, 1e-3, 0.5, 1.0, 3.0, 100.0]:
        for x in [1e-10, 0.001, 0.3, 0.5, 0.7, 0.999]:
            for aa, bb in ((a, b), (b, a)):
                A, B, X = mp.mpf(aa), mp.mpf(bb), mp.mpf(x)
                bt['x'].append(x); bt['a'].append(aa); bt['b'].append(bb)
                bt['reg'].append(mp.betainc(A, B, 0, X, regularized=True))
                bt['inc'].append(mp.betainc(A, B, 0, X))

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..',
                   'test', 'resources', 'special', 'gamma_beta_reference.edn')
with open(out, 'w') as f:
    f.write('{' + '\n '.join([
        block('gamma-small', [(k, gs[k]) for k in ('a', 'x', 'p', 'q', 'upper')]),
        block('gamma-large', [(k, gl[k]) for k in ('a', 'x', 'p', 'q')]),
        block('beta-poles', [(k, bp[k]) for k in ('p', 'q', 'ref')]),
        block('beta-tiny', [(k, bt[k]) for k in ('x', 'a', 'b', 'reg', 'inc')])]) + '}\n')
