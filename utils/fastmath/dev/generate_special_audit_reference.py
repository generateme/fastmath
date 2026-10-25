"""Generator of the independent reference values of the `fastmath.special` audit
(`fastmath.special.audit-test`).

Produces `test/resources/special/audit_reference.edn`, computed with `mpmath`
(60 digits): a map `{function-name [[arg... reference]...]}` for the functions of
`fastmath.special` that had no (or too few) direct tests, plus a `"hankel"` block
`[[name order x re im]...]` for the complex Hankel functions, and a `"complex"` block
`[[name [args...] [re im]]...]` for the functions of a complex argument (complex numbers are `[re im]`).  The arguments are the
exact doubles written to the file; infinities are `##Inf`, `##-Inf`.

Blocks (one function per section below):
  - erf family, inverse erf/erfc (the inverse of erfc in the far tail by
    a root finder on `log erfc`, `erfinv(1 - x)` cannot be used there), gamma
    (747 points on [-185, 171.7]; subnormal results are left out), log-gamma,
    `log-gamma-1p`, `inv-gamma-1pm1`, incomplete gamma (s > 0 and s <= 0),
    regularized gamma;
  - Bessel J, Y, I, K of orders 0, 1 and real orders, spherical Bessel, jinc
    (normalised, 2 J1(pi x) / (pi x)), Hankel (principal branch, x < 0 too);
  - beta (also negative arguments), log-beta, regularized/incomplete beta,
    digamma, trigamma, polygamma, zeta, eta, Dirichlet beta, harmonic numbers;
  - Airy (up to -1e6), Lambert W, Minkowski, Owen's T;
  - elliptic K, E, F, D, Pi (incomplete, all amplitudes), Carlson R_*, Jacobi
    functions (m in [-2, 3]);
  - confluent and Gauss hypergeometric: `hypergeometric-2F1` (1500 random moderate
    parameters, seed 5), `hypergeometric-0F1`, `0F2` (x >= -50), `kummers-M` (x >= -10),
    `whittaker-M`, `tricomis-U` (a > 0, a an integer, or 0 < b <= 8.5),
    `whittaker-W` (random parameters, seed 11).

Run with (from the repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_special_audit_reference.py
Regeneration is byte-identical (fixed point sets and seeds, `repr` of doubles).
"""
import fractions
import math
import random

import mpmath as mp

mp.mp.dps = 60

OUT = "/home/ts/clojure/fastmath/test/resources/special/audit_reference.edn"

CASES = {}
HANKEL = []
COMPLEX = []


def put(fn, args, ref):
    CASES.setdefault(fn, []).append([float(a) for a in args] + [ref])


def add(fn, args, f):
    """Appends `[args... f(args...)]` to the block `fn`; skips complex and non-finite values."""
    try:
        v = f(*[mp.mpf(a) for a in args])
        if isinstance(v, mp.mpc) or not mp.isfinite(v):
            return
        put(fn, args, float(v))
    except Exception:
        pass


def add_overflow(fn, args, f):
    """As `add`, but a value outside the double range becomes an infinity."""
    try:
        v = f(*[mp.mpf(a) for a in args])
        if isinstance(v, mp.mpc):
            return
        put(fn, args, float(v) if abs(v) <= mp.mpf('1.7976931348623157e308') else (math.inf if v > 0 else -math.inf))
    except Exception:
        pass


# ----------------------------------------------------------------------------- erf, gamma

def erf_gamma():
    xs = [-30, -10, -5, -2.5, -1, -0.5, -1e-3, 1e-8, 1e-3, 0.1, 0.5, 1, 2, 3.5, 6, 10, 20, 27, 40, 100]
    for x in xs:
        add('erf', [x], mp.erf)
        add('erfc', [x], mp.erfc)
    for x in [-30, -10, -5, -1, -0.5, 0.5, 1, 5, 10, 26, 27.5, 30]:
        add('erfc', [x], mp.erfc)
    for x in [-0.999999, -0.9, -0.5, -1e-8, 1e-8, 0.1, 0.5, 0.9, 0.999, 0.999999, 1 - 1e-12]:
        add('inv-erf', [x], mp.erfinv)

    def inverfc(x):
        if x >= 1:
            return mp.erfinv(1 - x) if x <= 1.5 else -mp.erfinv(x - 1)
        if x > mp.mpf('1e-3'):
            return mp.erfinv(1 - x)
        t0 = mp.sqrt(-mp.log(x) - mp.log(mp.pi * (-mp.log(x))) / 2)
        return mp.findroot(lambda t: mp.log(mp.erfc(t)) - mp.log(x), t0, tol=mp.mpf(10) ** (-40))

    for x in [1e-300, 1e-200, 1e-100, 1e-50, 1e-20, 1e-15, 1e-10, 1e-5, 1e-3, 0.01, 0.1, 0.5, 1, 1.5, 1.9,
              2 - 1e-10, 2 - 1e-3]:
        add('inv-erfc', [x], inverfc)

    for x in [-100.5, -10.3, -2.5, -0.5, 1e-5, 0.1, 0.5, 1, 2, 3, 10, 100, 1e5, 1e10, 1e100, 1e300]:
        add('log-gamma', [x], mp.loggamma)
    for x in [-0.5, -0.4, -0.1, -1e-6, 1e-6, 0.1, 0.5, 1, 1.5]:
        add('log-gamma-1p', [x], lambda x: mp.loggamma(1 + x))
        add('inv-gamma-1pm1', [x], lambda x: 1 / mp.gamma(1 + x) - 1)

    rnd = random.Random(7)
    xs = [rnd.uniform(120, 171.6) for _ in range(150)] + [140.0, 140.5, 141, 142, 143, 150, 160, 171, 171.5,
                                                          171.6, 171.62, 171.7]
    xs += [-rnd.uniform(100, 184) for _ in range(250)] + [-170.5, -171.5, -172.5, -180.5, -183.2, -185.5,
                                                          -140.3, -141.5, -139.5]
    xs += [rnd.uniform(-100, 140) for _ in range(300)] + [1e-310, 5e-324, 1e-300, -1e-300, 0.5, 1, 2, 10, 20.5, 100]
    for x in xs:
        v = mp.gamma(mp.mpf(x))
        if 0 < abs(v) < 2.3e-308:  # subnormal results have fewer digits
            continue
        add_overflow('gamma', [x], mp.gamma)

    pos = [1e-10, 1e-5, 1e-3, 0.01, 0.1, 0.5, 1, 2, 5, 10, 30, 100, 500]
    for s in [0.1, 0.5, 1, 2, 5, 20, 100]:
        for x in pos:
            add('lower-incomplete-gamma', [s, x], lambda s, x: mp.gammainc(s, 0, x))
            add('upper-incomplete-gamma', [s, x], lambda s, x: mp.gammainc(s, x))
            add('regularized-gamma-p', [s, x], lambda a, x: mp.gammainc(a, 0, x, regularized=True))
            add('regularized-gamma-q', [s, x], lambda a, x: mp.gammainc(a, x, mp.inf, regularized=True))
    for s in [-0.5, -1.5, 0, -1]:
        for x in [0.1, 1, 5]:
            add('upper-incomplete-gamma', [s, x], lambda s, x: mp.gammainc(s, x))
    for s in [-0.5, -1.5, -2.5, -0.1, -7.3]:  # analytic continuation to s <= 0
        for x in [0.01, 0.1, 0.5, 1, 2, 5, 20]:
            add('lower-incomplete-gamma', [s, x], lambda s, x: mp.gammainc(s, 0, x))
            add('regularized-gamma-p', [s, x], lambda a, x: mp.gammainc(a, 0, x, regularized=True))
            add('regularized-gamma-q', [s, x], lambda a, x: mp.gammainc(a, x, mp.inf, regularized=True))


# ----------------------------------------------------------------------------- Bessel

def bessel():
    for x in [-50, -12.3, -5, -1.5, -0.5, -1e-3, 1e-6, 0.01, 0.3, 1, 2, 5, 9.5, 20, 60, 150, 700]:
        add('bessel-J0', [x], lambda x: mp.besselj(0, x))
        add('bessel-J1', [x], lambda x: mp.besselj(1, x))
        add('bessel-Y0', [x], lambda x: mp.bessely(0, x))
        add('bessel-Y1', [x], lambda x: mp.bessely(1, x))
        add('bessel-I0', [x], lambda x: mp.besseli(0, x))
        add('bessel-I1', [x], lambda x: mp.besseli(1, x))
        add('bessel-K0', [x], lambda x: mp.besselk(0, x))
        add('bessel-K1', [x], lambda x: mp.besselk(1, x))
        for v in [0, 2, 5, 17, 50, 0.5, 1.5, 3.7, 25.3, 100, -1.5, -3.2]:
            add_overflow('bessel-J', [v, x], mp.besselj)
            add_overflow('bessel-Y', [v, x], mp.bessely)
            add_overflow('bessel-I', [v, x], mp.besseli)
            add_overflow('bessel-K', [v, x], mp.besselk)

    def sph(f, sign):
        return lambda n, x: sign * mp.sqrt(mp.pi / (2 * x)) * f(n + 0.5, x)

    for n in range(0, 6):
        for x in [0.01, 0.5, 1, 3, 10, 50, 200]:
            add('spherical-bessel-j', [n, x], sph(mp.besselj, 1))
            add('spherical-bessel-y', [n, x], sph(mp.bessely, 1))
            add('spherical-bessel-1-i', [n, x], sph(mp.besseli, 1))
            add('spherical-bessel-2-i', [n, x], lambda n, x: mp.sqrt(mp.pi / (2 * x)) * mp.besseli(-n - 0.5, x))
            add('spherical-bessel-k', [n, x], sph(mp.besselk, 1))
        for x in [-0.01, -0.5, -3, -20]:  # parity for integer orders
            xx = mp.mpf(x)
            put('spherical-bessel-j', [n, x], float((-1) ** n * mp.sqrt(mp.pi / (2 * -xx)) * mp.besselj(n + 0.5, -xx)))
            put('spherical-bessel-y', [n, x], float((-1) ** (n + 1) * mp.sqrt(mp.pi / (2 * -xx)) * mp.bessely(n + 0.5, -xx)))

    for x in [-30, -10, -3, -1.5, -1, -1e-6, 1e-8, 1e-4, 0.01, 0.1, 0.5, 1, 1.2197, 2, 3, 5, 10, 30, 100, 1e3, 1e5]:
        xx = mp.mpf(x)
        put('jinc', [x], float(2 * mp.besselj(1, mp.pi * xx) / (mp.pi * xx)))

    for nu in [0, 0.5, 1, 2.3, -1.5, 3]:
        for x in [-5, -1, -0.1, 0.1, 1, 5, 20]:
            h1 = mp.hankel1(mp.mpf(nu), mp.mpf(x))
            h2 = mp.hankel2(mp.mpf(nu), mp.mpf(x))
            HANKEL.append(['hankel-1', nu, x, float(h1.real), float(h1.imag)])
            HANKEL.append(['hankel-2', nu, x, float(h2.real), float(h2.imag)])


# ----------------------------------------------------------------------------- beta, zeta, polygamma

def mink(x):
    """Minkowski ?(x) on [0, 1) from the continued fraction of the exact rational x."""
    f = fractions.Fraction(x)
    a = math.floor(f)
    fr = f - a
    cf, k = [], 0
    while fr != 0 and k < 60:
        fr = 1 / fr
        ai = math.floor(fr)
        fr -= ai
        cf.append(ai)
        k += 1
    sign, e, tot = 1, 0, mp.mpf(0)
    for ai in cf:
        e += ai
        tot += sign * mp.mpf(2) ** (-e)
        sign = -sign
    return a + 2 * tot


def beta_zeta():
    for x in [-50, -12.3, -5, -1.5, -0.5, -1e-3, 1e-6, 0.01, 0.3, 1, 2, 5, 9.5, 20, 60, 150, 700]:
        add('digamma', [x], mp.digamma)
        add('trigamma', [x], lambda x: mp.psi(1, x))
        for m in [2, 3, 5, 10]:
            add('polygamma', [m, x], lambda m, x: mp.psi(int(m), x))
        add('zeta', [x], mp.zeta)
        add('eta', [x], mp.altzeta)
        add('dirichlet-beta', [x], lambda x: mp.dirichlet(x, [0, 1, 0, -1]))
        for p in [0.3, 2, 40, 1e-3, 100.5]:
            # beta(0.3, -12.3): the double sum is exactly -12 (a pole cancellation), not the reference
            if (p, x) != (0.3, -12.3):
                add('beta', [p, x], mp.beta)
                add('log-beta', [p, x], lambda p, q: mp.log(abs(mp.beta(p, q))))
            add('zeta', [x, p], mp.zeta)
    CASES['dirichlet-beta'] = [c for c in CASES['dirichlet-beta'] if c[0] != -5.0]  # trivial zero: absolute 2e-15
    for x in [0.001, 0.2, 0.5, 0.9, 0.999]:
        for a in [0.5, 2, 30, 500]:
            for b in [0.5, 3, 70, 900]:
                add('regularized-beta', [x, a, b], lambda x, a, b: mp.betainc(a, b, 0, x, regularized=True))
                add('incomplete-beta', [x, a, b], lambda x, a, b: mp.betainc(a, b, 0, x))
    for p, q in [(-1, 1), (1, -1), (2, -5), (-5, 2), (3, -4), (4, -10), (1, -30)]:  # poles cancelling
        v = mp.beta(mp.mpf(p), mp.mpf(q))
        put('beta', [p, q], float(v) if mp.isfinite(v) else math.inf)
    for n in [0, 1, 2, 5, 10, 100, 1000, 1e6]:
        add('harmonic-number', [n], lambda n: mp.harmonic(n))
    for n in [1, 2, 5, 10, 100]:
        for m in [2, 3, 0.5, 10]:
            add('harmonic-number', [n, m], lambda n, m: mp.zeta(m) - mp.zeta(m, n + 1))
    for x in [0.1, 0.25, 1 / 3, 0.5, 0.618033988749895, 0.7, 0.9, 1.5, -0.3, 3.14159, 10.25, 1e-5]:
        put('minkowski', [x], float(mink(x)))


# ----------------------------------------------------------------------------- Airy, Lambert, Owen

def airy_lambert_owen():
    # entire cosine integral Cin(x) = gamma + ln x - Ci(x): small x (the difference cancels), both sides of the switch x = 2
    for x in [1e-12, 1e-8, 1e-6, 1e-4, 1e-3, 0.01, 0.1, 0.5, 1, 1.5, 1.999, 2.001, 3, 5, 10, 50]:
        add('Cin', [x], lambda x: mp.euler + mp.log(x) - mp.ci(x))
    for x in [-1e6, -1e5, -1e4, -1e3, -100, -30, -10, -5, -2, -1, -0.1, 0, 0.1, 1, 2, 5, 10, 30, 100]:
        X = mp.mpf(x)
        put('airy-Ai', [x], float(mp.airyai(X)))
        put("airy-Ai'", [x], float(mp.airyai(X, derivative=1)))
        put('airy-Bi', [x], float(mp.airybi(X)))
        put("airy-Bi'", [x], float(mp.airybi(X, derivative=1)))
    for x in [-0.36787, -0.3, -0.1, -1e-8, 1e-8, 0.1, 1, 2.718281828, 10, 100, 1e4, 1e10, 1e100, 1e300]:
        add('lambert-W', [x], lambda x: mp.lambertw(x, 0).real)
    for x in [-0.367879, -0.36, -0.3, -0.1, -1e-3, -1e-8, -1e-100, -1e-300]:
        add('lambert-W-1', [x], lambda x: mp.lambertw(x, -1).real)

    def owen(h, a):
        return mp.quad(lambda t: mp.exp(-h * h * (1 + t * t) / 2) / (1 + t * t), [0, a]) / (2 * mp.pi)

    for h in [0, 0.01, 0.5, 1, 2, 5, 10, -1]:
        for a in [0, 0.01, 0.5, 1, 2, 10, 100, -0.5, 1e6]:
            add('owens-t', [h, a], owen)
    for h in [0.01, 1, 3]:
        for a in [-3, 0.1, 1, 4]:
            add('owens-t', [h, a], owen)


# ----------------------------------------------------------------------------- elliptic

def elliptic():
    def D(phi, m):
        if m != 0:
            return (mp.ellipf(phi, m) - mp.ellipe(phi, m)) / m
        return mp.quad(lambda t: mp.sin(t) ** 2 / mp.sqrt(1 - m * mp.sin(t) ** 2), [0, phi])

    for m in [-5, -1, -0.1, 0, 0.1, 0.5, 0.9, 0.999, 0.999999, 1]:
        add('elliptic-K', [m], mp.ellipk)
        add('elliptic-E', [m], mp.ellipe)
        for phi in [0.1, 0.5, 1, 1.5, 2, 3, 5, -1]:
            add('elliptic-F', [phi, m], mp.ellipf)
            add('elliptic-E', [phi, m], mp.ellipe)
            add('elliptic-PI', [0.3, phi, m], lambda n, p, m: mp.ellippi(n, p, m))
            add('elliptic-D', [phi, m], D)
        for n in [-2, -0.5, 0.3, 0.9]:  # n > 1 is a principal value (NaN)
            add('elliptic-PI', [n, m], mp.ellippi)
        add('elliptic-D', [m], lambda m: (mp.ellipk(m) - mp.ellipe(m)) / m if m != 0 else mp.pi / 4)
    for x, y, z in [(1, 2, 3), (0, 1, 2), (1e-10, 1, 1e10), (0.5, 0.5, 0.5), (5, 5, 1), (1e-5, 2, 3), (100, 1, 1e-3)]:
        add('elliptic-Rf', [x, y, z], mp.elliprf)
        add('elliptic-Rd', [x, y, z], mp.elliprd)
        add('elliptic-Rg', [x, y, z], mp.elliprg)
        for p in [0.5, 2, 10]:
            add('elliptic-Rj', [x, y, z, p], mp.elliprj)
    for x, y in [(0, 1), (1, 2), (2, 1), (1e-8, 1), (5, 5), (0.5, -1), (3, -0.5), (1, 1e10), (1e10, 1)]:
        add('elliptic-Rc', [x, y], mp.elliprc)
    names = ['sn', 'cn', 'dn', 'sc', 'sd', 'cs', 'cd', 'ds', 'dc', 'ns', 'nc', 'nd']
    for m in [0, 0.2, 0.7, 0.99]:
        for u in [-3, -0.5, 0.1, 0.5, 1, 2, 5, 10, 50]:
            if abs(u) < mp.ellipk(m):
                add('jacobi-am', [u, m], lambda u, m: mp.asin(mp.ellipfun('sn', u, m=m)))
    for m in [0, 0.2, 0.7, 0.99, -2, -0.5, 1.5, 3, 0.5]:  # m outside [0, 1]: real analytic continuation
        for u in ([-3, -0.5, 0.1, 0.5, 1, 2, 5, 10, 50] if 0 <= m < 1 else [-2.5, -0.7, 0.3, 0.9, 1.5, 2.5, 4]):
            for nm in names:
                add('jacobi-' + nm, [u, m], (lambda nm: lambda u, m: mp.ellipfun(nm, u, m=m))(nm))


# ----------------------------------------------------------------------------- hypergeometric

def hypergeometric():
    rnd = random.Random(5)
    for _ in range(1500):
        a = rnd.choice([0.5, 1, 2, 2.5, 3.7, 5, 8, 0.1, -0.5, -2.5])
        b = rnd.choice([0.5, 1.5, 3, 4, 6, 7.3, 9, -1.5])
        c = rnd.choice([0.5, 1, 2, 3, 4.5, 6, 8, 10, 12.5, -0.5, -3.5, 15])
        x = rnd.choice([-20, -5, -2, -0.9, -0.5, -0.1, 0.1, 0.3, 0.5, 0.72, 0.8, 0.9, 0.99])
        add('hypergeometric-2F1', [a, b, c, x], mp.hyp2f1)
    CASES['hypergeometric-2F1'] = [c for c in CASES['hypergeometric-2F1'] if c[-1] != 0.0]
    # regressions: c-a-b < 0 with c-b a non-positive integer (the Euler transformation terminates), integer c-a-b,
    # a negative non-integer parameter
    for a, b, c, x in [(0.1, 20, 2, 0.5), (30, 0.5, 2, 0.5), (0.1, 5, 2, 0.5), (2.5, 20, -4.5, -5), (10, 0.5, -4.5, -5),
                       (10, 20, 30, 0.5), (1, 50, 60, 0.5), (30, 20, 60, 0.5), (75, -2.5, 30, 0.5), (-2.5, 0.5, -0.5, 0.5),
                       (-2.5, 1.5, 0.5, 0.99), (2.5, 1.5, 15, 0.5), (5, 6, 15, 0.5), (8, 4, 15, 0.5)]:
        add('hypergeometric-2F1', [a, b, c, x], mp.hyp2f1)

    # the variogram models: 2F1(1/2, -k/2; 3/2; x), k up to 180 (terminating for odd k): direct sums cancel
    for k2 in range(1, 181, 3):
        for x in [0.01, 0.3, 0.5, 0.7, 0.95, 0.999]:
            add('hypergeometric-2F1', [0.5, -k2 / 2.0, 1.5, x], mp.hyp2f1)

    rnd = random.Random(11)
    xs = [-200, -50, -10, -3, -1, -0.1, 0.01, 0.3, 1, 2.5, 6, 15, 40, 120]
    A = [-6.5, -3, -0.5, 0.1, 0.5, 1, 2.5, 7, 20, 60]
    B = [-4.5, -1.5, 0.5, 1, 1.5, 3, 8.5, 25, 100]
    for _ in range(500):
        a, b, c, x = rnd.choice(A), rnd.choice(B), rnd.choice(B), rnd.choice(xs)
        add('hypergeometric-0F1', [b, x], lambda b, x: mp.hyp0f1(b, x))
        if x >= -10:  # large negative x: the value is below the absolute rounding noise
            add('kummers-M', [a, b, x], lambda a, b, x: mp.hyp1f1(a, b, x))
        if x >= -50:  # large negative x: alternating series, cancellation
            add('hypergeometric-0F2', [b, c, x], lambda b, c, x: mp.hyper([], [b, c], x))
        xp = rnd.choice([0.01, 0.1, 0.5, 1, 3, 10, 30, 100])
        if a > 0 or a == math.floor(a) or 0 < b <= 8.5:  # the documented range of tricomis-U for a < 0
            add_overflow('tricomis-U', [a, b, xp], lambda a, b, x: mp.hyperu(a, b, x))
        k = rnd.choice([-2.5, -1, 0, 0.5, 1, 3])
        mm = rnd.choice([0.25, 0.5, 1, 2, 3.5, 6])
        add('whittaker-M', [k, mm, xp], lambda k, m, x: mp.whitm(k, m, x))
        add('whittaker-W', [k, mm, xp], lambda k, m, x: mp.whitw(k, m, x))


# ----------------------------------------------------------------------------- complex arguments

def complex_functions():
    rnd = random.Random(21)

    def c(z):
        return [float(z.real), float(z.imag)]

    for _ in range(300):
        z = mp.mpc(rnd.choice([-30, -8.5, -3, -1.2, -0.4, 0.1, 0.6, 1.5, 4, 12, 40]),
                   rnd.choice([-25, -6, -2, -0.7, 0, 0.3, 1.1, 3.5, 9, 30]))
        if not (z.imag == 0 and z.real <= 0 and z.real == int(z.real)):   # the poles are covered by the edge tests
            COMPLEX.append(['gamma-complex', [c(z)], c(mp.gamma(z))])
            COMPLEX.append(['log-gamma-complex', [c(z)], c(mp.loggamma(z))])  # the analytic branch (Im not wrapped)
        a = mp.mpc(rnd.choice([0.3, 1.5, -0.7, 3, 7, -2.5]), rnd.choice([0, 0.4, -1.2, 2]))
        b = mp.mpc(rnd.choice([0.5, 1.5, 2.5, 3.2, -1.5, 6.5]), rnd.choice([0, 0.3, -0.8]))
        zu = mp.mpc(rnd.choice([0.2, 1, 3, 10, 40, -2, -8]), rnd.choice([0, 0.5, -3, 6, 20]))
        u = mp.hyperu(a, b, zu)
        if not (zu.imag == 0 and zu.real < 0):   # on the cut, and the documented limits of tricomis-U-complex
            if abs(u) < 1e300:
                COMPLEX.append(['tricomis-U-complex', [c(a), c(b), c(zu)], c(u)])
        n = rnd.choice([1, 3, 5, 7, 9, 13])
        k = mp.besselk(mp.mpf(n) / 2, z)
        COMPLEX.append(['bessel-K-half-odd-complex', [n, c(z)], c(k)])
        COMPLEX.append(['bessel-K-half-odd-scaled-complex', [n, c(z)], c(k * mp.exp(z))])
    ps = [[mp.mpc(0.5, 0.2), mp.mpc(1.3, -0.4)], [mp.mpc(-2.5, 0), mp.mpc(0.7, 0.3)], [mp.mpc(2, 1), mp.mpc(0.5, 0)]]
    qs = [[mp.mpc(2.5, 0.3)], [mp.mpc(1.5, -0.2)], [mp.mpc(3.5, 0.5)]]
    for p, q in zip(ps, qs):
        for _ in range(40):
            z = mp.mpc(rnd.choice([-0.9, -0.4, 0.1, 0.5, 0.8, -3, -1.5]), rnd.choice([-0.5, 0, 0.3, 0.9, 2]))
            COMPLEX.append(['hypergeometric-pFq-complex', [[c(x) for x in p], [c(x) for x in q], c(z)], c(mp.hyper(p, q, z))])


# ----------------------------------------------------------------------------- output

def edn(v):
    if isinstance(v, float):
        if math.isnan(v):
            return '##NaN'
        if math.isinf(v):
            return '##Inf' if v > 0 else '##-Inf'
        return repr(v)
    return repr(float(v))


def edn_nested(x):
    if isinstance(x, list):
        return '[' + ' '.join(edn_nested(e) for e in x) + ']'
    return edn(float(x)) if not isinstance(x, int) else repr(x)


def main():
    erf_gamma()
    bessel()
    beta_zeta()
    airy_lambert_owen()
    elliptic()
    hypergeometric()
    complex_functions()
    with open(OUT, 'w') as f:
        f.write(';; mpmath %s reference values for fastmath.special.audit-test: {function-name [[args... ref]...]}\n'
                ';; generated by utils/fastmath/dev/generate_special_audit_reference.py\n{\n' % mp.__version__)
        for fn in sorted(CASES):
            seen, rows = set(), []
            for c in CASES[fn]:
                key = tuple(repr(x) for x in c)
                if key not in seen:
                    seen.add(key)
                    rows.append(c)
            f.write(' "%s"\n [\n' % fn)
            for c in rows:
                f.write('  [' + ' '.join(edn(x) for x in c) + ']\n')
            f.write(' ]\n')
        f.write(' "hankel"\n [\n')
        for h in HANKEL:
            f.write('  ["%s" ' % h[0] + ' '.join(edn(float(x)) for x in h[1:]) + ']\n')
        f.write(' ]\n "complex"\n [\n')
        for name, args, ref in COMPLEX:
            f.write('  ["%s" %s %s]\n' % (name, edn_nested(args), edn_nested(ref)))
        f.write(' ]\n}\n')
    print({k: len(v) for k, v in sorted(CASES.items())}, len(HANKEL))


if __name__ == '__main__':
    main()
