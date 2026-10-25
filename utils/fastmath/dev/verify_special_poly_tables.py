"""Regenerates the coefficient tables of `fastmath.special.poly` from their closed forms
and compares them with the values in the source, coefficient by coefficient (a value test
at the usual 1e-9 tolerance cannot find a slip in the last retained term).

Input: the JSON written by `dump_special_poly_tables.clj`.  Prints one line per table with the
largest deviation and the list of coefficients outside the tolerance; exit status 1 if any.

Tables and their definitions:
  - `uk-poly-10`, `uk-poly-20`, `uk-poly-jn`: Debye polynomials `u_k(p) = p^k * sum_j c_j p^(2j)` of the
    uniform asymptotic expansions of Bessel functions, `u_0 = 1`,
    `u_(k+1)(p) = p^2 (1 - p^2) u_k'(p) / 2 + 1/8 int_0^p (1 - 5 t^2) u_k(t) dt` (exact rationals;
    `uk-poly-jn` is stored as integers over a common denominator and must match exactly).
  - `a-ap-poly-10/20/30`: `s_k(mu)`, `mu = 4 nu^2`, the coefficients of `x^(-2k)` in the derivative of the
    phase of the Hankel function, `theta'(x) = 1 / (P^2 + Q^2)` with `P + iQ = sum_k i^k a_k(nu) x^-k`,
    `a_k = prod_(j<=k)(mu - (2j-1)^2) / (k! 8^k)` (exact rationals).
  - `j0-roots`, `j1-roots` (hi + lo double-double): the zeros of `J0 * J1` and `J1 * J1'`;
    `j0-polys`, `j1-polys`: Taylor coefficients of `J0`, `J1` about those roots (mpmath).
  - `abi-1`, `abi-2`, `abip-1`, `abip-2`: Chebyshev coefficients of `Bi(x) exp(-2/3 x^1.5)` and
    `Bi'(x) exp(-2/3 x^1.5)` on [2, 6] and [6, 10]; `hg-2F1-poly`: of `2F1(1, 3/2; 5/2; t/5)` on [-1, 1].
  - `bessel-y-chebyshev-weights`: the 2-D Chebyshev coefficients of `Y_(t+1)(x)`, `t` in [-1, 1] (orders
    0..2), `x` in [6, 19] (`x1 = 2 (x - 6) / 13 - 1`).

Run with (from the repo root):
    clojure -M utils/fastmath/dev/dump_special_poly_tables.clj /tmp/special_poly_tables.json
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/verify_special_poly_tables.py /tmp/special_poly_tables.json
"""
import json
import sys
from fractions import Fraction
from math import comb

import mpmath as mp

mp.mp.dps = 50

FAILED = []


def report(name, worst, bad, count, tol):
    status = 'ok' if not bad else 'FAIL'
    print('%-28s %-5s n=%-4d worst=%.2e (tol %.0e)%s' % (
        name, status, count, worst, tol, '' if not bad else '  bad: ' + ', '.join(bad[:8])))
    if bad:
        FAILED.append(name)


def rel(a, b):
    """Relative difference of the mp numbers a (reference) and b (table)."""
    a, b = mp.mpf(a), mp.mpf(b)
    if a == 0:
        return abs(b)
    return abs(a - b) / abs(a)


# ----------------------------------------------------------------------------- Debye polynomials

def debye_polys(K):
    """Exact Debye polynomials u_0..u_K as {power of p: Fraction}."""
    u = [{0: Fraction(1)}]
    for k in range(K):
        cur = u[-1]
        nxt = {}
        for e, c in cur.items():
            # p^2 (1 - p^2) u' / 2
            if e > 0:
                nxt[e + 1] = nxt.get(e + 1, 0) + c * e / 2
                nxt[e + 3] = nxt.get(e + 3, 0) - c * e / 2
            # 1/8 int (1 - 5 t^2) t^e dt
            nxt[e + 1] = nxt.get(e + 1, 0) + c / (8 * (e + 1))
            nxt[e + 3] = nxt.get(e + 3, 0) - 5 * c / (8 * (e + 3))
        u.append(nxt)
    return u


def check_uk(data, name, K):
    u = debye_polys(K)
    worst, bad, count = 0.0, [], 0
    exact = name == 'uk-poly-jn'
    for key, entry in data[name].items():
        k = int(key[1:])
        ref = [u[k][k + 2 * j] for j in range(k + 1)]
        coeffs = entry['coeffs']
        if len(coeffs) != k + 1:
            bad.append('%s length %d != %d' % (key, len(coeffs), k + 1))
            continue
        for j, (r, c) in enumerate(zip(ref, coeffs)):
            count += 1
            if exact:
                table = Fraction(int(c), int(entry['denominator']))
                if table != r:
                    bad.append('%s[%d]' % (key, j))
                    worst = max(worst, float(abs(table - r) / abs(r)))
            else:
                d = float(rel(mp.mpf(r.numerator) / r.denominator, c))
                worst = max(worst, d)
                if d > 4e-16:
                    bad.append('%s[%d] %.1e' % (key, j, d))
    report(name, worst, bad, count, 0 if exact else 4e-16)


# ----------------------------------------------------------------------------- phase of the Hankel function

def theta_prime_coeffs(K, mus):
    """For every integer mu in `mus`: exact series coefficients s_1..s_K of 1/(P^2 + Q^2) in z^2, z = 1/x."""
    out = {}
    for mu in mus:
        a = [Fraction(1)]
        for k in range(1, 2 * K + 1):
            a.append(a[-1] * (mu - (2 * k - 1) ** 2) / (8 * k))
        # P = sum (-1)^j a_(2j) z^(2j); Q = sum (-1)^j a_(2j+1) z^(2j+1)
        n = 2 * K
        pq = [Fraction(0)] * (n + 1)   # coefficients of z^i in P^2 + Q^2 (only even i survive)
        P = [(-1) ** (i // 2) * a[i] if i % 2 == 0 else Fraction(0) for i in range(n + 1)]
        Q = [(-1) ** (i // 2) * a[i] if i % 2 == 1 else Fraction(0) for i in range(n + 1)]
        for i in range(n + 1):
            for j in range(n + 1 - i):
                pq[i + j] += P[i] * P[j] + Q[i] * Q[j]
        inv = [Fraction(0)] * (n + 1)
        inv[0] = Fraction(1)
        for i in range(1, n + 1):
            inv[i] = -sum(pq[j] * inv[i - j] for j in range(1, i + 1))
        out[mu] = [inv[2 * k] for k in range(1, K + 1)]
    return out


def interpolate(xs, ys):
    """Exact Newton interpolation: the coefficients (low to high) of the polynomial through (xs, ys)."""
    n = len(xs)
    coef = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (xs[i] - xs[i - j])
    poly = [Fraction(0)]
    for i in range(n - 1, -1, -1):
        # poly = poly * (x - xs[i]) + coef[i]
        new = [Fraction(0)] * (len(poly) + 1)
        for d, c in enumerate(poly):
            new[d + 1] += c
            new[d] -= c * xs[i]
        new[0] += coef[i]
        poly = new
    return poly[:n]


def check_aap(data):
    K = 30
    mus = list(range(0, K + 1))
    series = theta_prime_coeffs(K, mus)
    for name in ['a-ap-poly-10', 'a-ap-poly-20', 'a-ap-poly-30']:
        worst, bad, count = 0.0, [], 0
        for key, entry in data[name].items():
            k = int(key[1:])
            ys = [series[mu][k - 1] for mu in mus[:k + 1]]
            ref = interpolate([Fraction(m) for m in mus[:k + 1]], ys)
            coeffs = entry['coeffs']
            if len(coeffs) != k + 1:
                bad.append('%s length %d != %d' % (key, len(coeffs), k + 1))
                continue
            for j, (r, c) in enumerate(zip(ref, coeffs)):
                count += 1
                rv = mp.mpf(r.numerator) / r.denominator
                d = float(rel(rv, c))
                worst = max(worst, d)
                if d > 4e-16:
                    bad.append('%s[%d] %.1e' % (key, j, d))
        report(name, worst, bad, count, 4e-16)


# ----------------------------------------------------------------------------- Bessel zeros and Taylor series

def bessel_taylor(order, r0, n):
    """Taylor coefficients c_0..c_(n-1) of J_order(r0 + h) in h from the ODE x^2 y'' + x y' + (x^2 - order^2) y = 0."""
    c = [mp.besselj(order, r0), mp.besselj(order, r0, derivative=1)]
    nu2 = mp.mpf(order) ** 2
    for k in range(0, n - 2):
        cm1 = c[k - 1] if k >= 1 else 0
        cm2 = c[k - 2] if k >= 2 else 0
        c.append(-((2 * k + 1) * (k + 1) * r0 * c[k + 1] + (k * k + r0 * r0 - nu2) * c[k] + 2 * r0 * cm1 + cm2)
                 / (r0 * r0 * (k + 2) * (k + 1)))
    return c[:n]


def check_j_roots_and_polys(data):
    for fn, order in [('j0', 0), ('j1', 1)]:
        roots = data[fn + '-roots']
        worst, bad = 0.0, []
        for k, (hi, lo) in enumerate(roots):
            r = mp.mpf(hi) + mp.mpf(lo)
            # nearest stationary point / zero: Newton on the product f * f'
            f = (lambda x: mp.besselj(order, x) * mp.besselj(order, x, derivative=1)) if order == 1 else \
                (lambda x: mp.besselj(0, x) * mp.besselj(1, x))
            ex = mp.findroot(f, r, tol=mp.mpf(10) ** -40)
            d = float(abs(ex - r) / ex)
            worst = max(worst, d)
            if d > 2e-16:   # the low part of the double-double is only accurate to about 1e-16
                bad.append('root %d %.1e' % (k, d))
        report(fn + '-roots', worst, bad, len(roots), 2e-16)

        # The polynomials are fits (not Taylor series) about the root: they are checked as functions on the
        # window of `x` where `bessel-J0` / `bessel-J1` use them, [(n + 1) pi / 2, min((n + 2) pi / 2, 26)]
        # (see the callers in fastmath.special).
        worst, bad, count = 0.0, [], 0
        for key, entry in data[fn + '-polys'].items():
            k = int(key)
            hi, lo = roots[k]
            r0 = mp.mpf(hi) + mp.mpf(lo)
            lo_x, hi_x = (k + 1) * mp.pi / 2, min((k + 2) * mp.pi / 2, mp.mpf(26))   # used for x < 26
            coeffs = [mp.mpf(c) for c in entry['coeffs']]
            err = mp.mpf(0)
            N = 400
            for i in range(N + 1):
                x = lo_x + (hi_x - lo_x) * i / N
                h = x - r0
                val = mp.mpf(0)
                for c in reversed(coeffs):
                    val = val * h + c
                err = max(err, abs(val - mp.besselj(order, x)))
                count += 1
            e = float(err)
            worst = max(worst, e)
            if e > 5e-16:
                bad.append('%s %.1e' % (key, e))
        report(fn + '-polys (as functions)', worst, bad, count, 5e-16)


# ----------------------------------------------------------------------------- Chebyshev tables

def chebyshev_coeffs(f, n, oversample=4):
    """Coefficients c_0..c_(n-1) with f ~ sum c_k T_k on [-1, 1] (c_0 is the mean, not halved)."""
    M = oversample * n
    nodes = [mp.cos(mp.pi * (j + mp.mpf(1) / 2) / M) for j in range(M)]
    vals = [f(x) for x in nodes]
    out = []
    for k in range(n):
        s = sum(v * mp.cos(k * mp.pi * (j + mp.mpf(1) / 2) / M) for j, v in enumerate(vals))
        out.append(s / M if k == 0 else 2 * s / M)
    return out


def check_cheb(name, table, f, tol_abs=3e-16):
    ref = chebyshev_coeffs(f, len(table))
    worst, bad = 0.0, []
    for j, (r, c) in enumerate(zip(ref, table)):
        d = float(abs(r - mp.mpf(c)))
        scale = float(abs(ref[0]))
        worst = max(worst, d)
        if d > tol_abs * max(1.0, scale):
            bad.append('[%d] %.1e' % (j, d))
    report(name, worst, bad, len(table), tol_abs)


def check_chebyshev(data):
    def Bi_scaled(lo, hi, deriv):
        def f(t):
            x = lo + (hi - lo) * (t + 1) / 2
            e = mp.exp(-mp.mpf(2) / 3 * x ** mp.mpf('1.5'))
            return (mp.airybi(x, derivative=1) if deriv else mp.airybi(x)) * e
        return f
    check_cheb('abi-1', data['abi-1'], Bi_scaled(2, 6, False))
    check_cheb('abi-2', data['abi-2'], Bi_scaled(6, 10, False))
    check_cheb('abip-1', data['abip-1'], Bi_scaled(2, 6, True))
    check_cheb('abip-2', data['abip-2'], Bi_scaled(6, 10, True))
    check_cheb('hg-2F1-poly', data['hg-2F1-poly'],
               lambda t: mp.hyp2f1(1, mp.mpf(3) / 2, mp.mpf(5) / 2, t / 5))


def check_bessel_y_weights(data):
    """2-D Chebyshev coefficients: weights[i][j] multiplies T_i(nu) T_j(x1), x1 = 2 (x - 6) / 13 - 1."""
    W = data['bessel-y-chebyshev-weights']
    ni, nj = len(W), len(W[0])
    M_nu, M_x = 2 * ni, 2 * nj
    nus = [mp.cos(mp.pi * (a + mp.mpf(1) / 2) / M_nu) for a in range(M_nu)]
    xs = [mp.cos(mp.pi * (b + mp.mpf(1) / 2) / M_x) for b in range(M_x)]
    mp.mp.dps = 30
    vals = [[mp.bessely(nu + 1, 6 + 13 * (x1 + 1) / 2) for x1 in xs] for nu in nus]
    worst, bad = 0.0, []
    for i in range(ni):
        for j in range(nj):
            s = mp.mpf(0)
            for a in range(M_nu):
                ci = mp.cos(i * mp.pi * (a + mp.mpf(1) / 2) / M_nu)
                for b in range(M_x):
                    s += vals[a][b] * ci * mp.cos(j * mp.pi * (b + mp.mpf(1) / 2) / M_x)
            c = s / (M_nu * M_x) * (1 if i == 0 else 2) * (1 if j == 0 else 2)
            d = float(abs(c - mp.mpf(W[i][j])))
            worst = max(worst, d)
            if d > 3e-16:
                bad.append('[%d][%d] %.1e' % (i, j, d))
    mp.mp.dps = 50
    report('bessel-y-chebyshev-weights', worst, bad, ni * nj, 3e-16)


def main(path):
    data = json.load(open(path))
    check_uk(data, 'uk-poly-10', 10)
    check_uk(data, 'uk-poly-20', 20)
    check_uk(data, 'uk-poly-jn', 21)
    check_aap(data)
    check_j_roots_and_polys(data)
    check_chebyshev(data)
    check_bessel_y_weights(data)
    print('FAILED: %s' % FAILED if FAILED else 'all tables match')
    sys.exit(1 if FAILED else 0)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '/tmp/special_poly_tables.json')
