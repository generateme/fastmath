# Reference complex roots of random polynomials for test/fastmath/solver/polyroot_test.clj
# The coefficients are doubles (real: [c0 c1 ...]; complex: [[re im] ...]), the roots are computed from them (mpmath, 150 digits)
# and rounded to doubles. The roots are well separated (relative gap >= 1e-2) so the problem is well conditioned.
# run: cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_solver_polyroot_reference.py
import random, cmath, math, mpmath as mp
mp.mp.dps = 150
rnd = random.Random(20261027)

def separated(rs):
    for i, x in enumerate(rs):
        for y in rs[i + 1:]:
            if abs(x - y) < 1e-2 * max(abs(x), abs(y)):
                return False
    return True

def poly_from_roots(rs, lead):
    co = [mp.mpc(lead)]
    for r in rs:
        # multiply by (x - r): coefficients in descending order
        new = co + [mp.mpc(0)]
        for i in range(len(co)):
            new[i + 1] -= r * co[i]
        co = new
    return co[::-1]  # ascending

def rand_root():
    return cmath.rect(10.0 ** rnd.uniform(-3, 3), rnd.uniform(0, 2 * math.pi))

def true_roots(co):
    n = len(co) - 1
    sc = [mp.mpc(c) for c in co]
    s = max([abs(sc[0] / sc[n]) ** (mp.mpf(1) / n)] + [abs(sc[i] / sc[n]) ** (mp.mpf(1) / (n - i)) for i in range(1, n) if sc[i] != 0])
    cs = [sc[i] * s ** i / (sc[n] * s ** n) for i in range(n + 1)]
    rs = mp.polyroots(cs[::-1], maxsteps=20000, extraprec=6000)
    return [r * s for r in rs]

def fl(x): return repr(float(x))

cases = []
while len(cases) < 100:
    n = rnd.randint(3, 20)
    rs = []
    while len(rs) < n:
        if n - len(rs) >= 2 and rnd.random() < 0.6:
            z = rand_root()
            if abs(z.imag) < 0.2 * abs(z):
                continue
            rs += [z, z.conjugate()]
        else:
            rs.append((1 if rnd.random() < .5 else -1) * 10.0 ** rnd.uniform(-3, 3) + 0j)
    if not separated(rs):
        continue
    co = poly_from_roots([mp.mpc(r) for r in rs], 10.0 ** rnd.uniform(-2, 2) * rnd.choice([-1, 1]))
    c = [float(mp.re(x)) for x in co]
    tr = true_roots(c)
    cases.append(('real', c, tr))
while len(cases) < 160:
    n = rnd.randint(3, 15)
    rs = [rand_root() for _ in range(n)]
    if not separated(rs):
        continue
    co = poly_from_roots([mp.mpc(r) for r in rs], complex(rnd.uniform(-2, 2), rnd.uniform(-2, 2)))
    c = [(float(mp.re(x)), float(mp.im(x))) for x in co]
    tr = true_roots([mp.mpc(a, b) for a, b in c])
    cases.append(('complex', c, tr))
with open('/home/ts/clojure/fastmath/test/resources/solver/polyroot_reference.edn', 'w') as f:
    f.write('[\n')
    for kind, c, tr in cases:
        cs = '[' + ' '.join(fl(x) if kind == 'real' else '[%s %s]' % (fl(x[0]), fl(x[1])) for x in c) + ']'
        rs = '[' + ' '.join('[%s %s]' % (fl(mp.re(r)), fl(mp.im(r))) for r in tr) + ']'
        f.write(' {:kind :%s :coefs %s :roots %s}\n' % (kind, cs, rs))
    f.write(']\n')
print(len(cases))
