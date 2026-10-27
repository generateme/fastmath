# Reference real roots of quadratics and cubics for test/fastmath/solver/closed_form_test.clj
# The coefficients are doubles, the roots are computed from them exactly (mpmath, 120 digits) and rounded to doubles.
# Cases have well separated roots (relative gap >= 1e-3, complex pairs at Im/|Re| >= 1e-3) so that the forward error of a backward stable method is small
# and the number of real roots does not depend on the rounding of the coefficients.
# run: cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_solver_closed_form_reference.py
import random, mpmath as mp
mp.mp.dps = 120
rnd = random.Random(20261027)

def sign_exp():
    return rnd.choice([-1, 1]) * 10.0 ** rnd.uniform(-9, 9)

def sep_ok(rs):
    rs = sorted(abs(x) for x in rs)
    return True

def quad_case():
    r1, r2 = sign_exp(), sign_exp()
    if abs(r1 - r2) < 1e-3 * max(abs(r1), abs(r2)):
        return None
    a = rnd.choice([-1, 1]) * 10.0 ** rnd.uniform(-3, 3)
    b = -a * (r1 + r2)
    c = a * r1 * r2
    return [a, b, c]

def real_roots(coefs):
    mc = [mp.mpf(x) for x in coefs]
    n = len(mc) - 1
    # scale x = s y to tame the dynamic range
    s = max([abs(mc[i] / mc[0]) ** (mp.mpf(1) / i) for i in range(1, n + 1) if mc[i] != 0] or [mp.mpf(1)])
    sc = [mc[i] / (mc[0] * s ** i) for i in range(n + 1)]
    rs = mp.polyroots(sc, maxsteps=20000, extraprec=4000)
    rs = [r * s for r in rs]
    out = [mp.re(r) for r in rs if abs(mp.im(r)) <= mp.mpf(10) ** -60 * max(abs(r), mp.mpf(10) ** -300)]
    return sorted(out)

def cubic_case(kind):
    if kind == 'three':
        rs = [sign_exp() for _ in range(3)]
        if min(abs(x - y) / max(abs(x), abs(y)) for i, x in enumerate(rs) for y in rs[i + 1:]) < 1e-3:
            return None
        cx = None
    else:
        r = sign_exp()
        re_ = sign_exp()
        im = abs(re_) * 10.0 ** rnd.uniform(-3, 1)  # Im/|Re| >= 1e-3: a pair this far from the real axis is not a double root within the rounding error
        # (x - r)(x^2 - 2 re x + re^2 + im^2)
        cx = (r, re_, im)
    a = rnd.choice([-1, 1]) * 10.0 ** rnd.uniform(-3, 3)
    if cx is None:
        r1, r2, r3 = rs
        co = [a, -a * (r1 + r2 + r3), a * (r1 * r2 + r1 * r3 + r2 * r3), -a * r1 * r2 * r3]
    else:
        r, re_, im = cx
        p, q = -2 * re_, re_ * re_ + im * im
        co = [a, a * (p - r), a * (q - r * p), -a * r * q]
    return [float(x) for x in co]

def fmt(xs):
    return '[' + ' '.join(repr(float(x)) for x in xs) + ']'

quads, cubics = [], []
while len(quads) < 300:
    c = quad_case()
    if c:
        c = [float(x) for x in c]
        quads.append((c, real_roots(c)))
while len(cubics) < 250:
    c = cubic_case('three')
    if c:
        rr = real_roots(c)
        if len(rr) == 3:
            cubics.append((c, rr))
while len(cubics) < 400:
    c = cubic_case('one')
    cubics.append((c, real_roots(c)))
with open('/home/ts/clojure/fastmath/test/resources/solver/closed_form_reference.edn', 'w') as f:
    f.write('{:quadratic [\n')
    for c, r in quads:
        f.write(' {:coefs %s :roots %s}\n' % (fmt(c), fmt(r)))
    f.write(' ]\n :cubic [\n')
    for c, r in cubics:
        f.write(' {:coefs %s :roots %s}\n' % (fmt(c), fmt(r)))
    f.write(' ]}\n')
print(len(quads), len(cubics), sum(len(r) == 3 for _, r in cubics), sum(len(r) == 1 for _, r in cubics))
