# Reference values for the log-exp family of fastmath.core: log2, logexpm1, logsumexp and every smooth-max family,
# with emphasis on the regions of the upstream group 1 fixes: powers of two, subnormals and the range-reduction
# boundaries near 1, sqrt(2) and the smallest normal (log2), overflow of
# expm1 (logexpm1), overflow and underflow of exp(alpha x) and |x|^alpha, ties and near ties (sigmoid kernel,
# mellowmax cancellation). Infinite and NaN arguments are asserted by hand in the tests.
#
# Tool: Python mpmath 1.3.0; 60 decimal digits, 700 digits when the elements are below 1e-100 or above 1e100
# (wide dynamic range) or form near ties (log(mean(exp)) cancels). The pairwise families use the kernels and constants of fastmath.core/smooth-max.
# Run in the /home/ts/penv environment:
#
#   cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_core_logexp_reference.py
#
# Output: test/resources/core/logexp_reference.edn
# {:log2 [[x y] ..], :logexpm1 [[x y] ..], :logsumexp [[xs y] ..], :smooth-max {<family> [[xs alpha y] ..]}}

import math
import os
import random
import mpmath as mp


def fmt(x):
    x = float(x)
    if x == float('inf'):
        return '##Inf'
    if x == float('-inf'):
        return '##-Inf'
    return repr(x)


def vec(xs):
    return '[' + ' '.join(fmt(x) for x in xs) + ']'


def precision(xs):
    small = max(abs(x) for x in xs) < 1e-100
    wide = max(abs(x) for x in xs) > 1e100
    spread = max(xs) - min(xs)
    near_tie = len(xs) > 1 and 0 < spread < 1e-6 * max(1.0, max(abs(x) for x in xs))
    return 700 if small or wide or near_tie else 60


rnd = random.Random(20261105)

# log2: random magnitudes, neighbours of powers of two, subnormals
log2_xs = [mp.e ** (rnd.uniform(-744, 744)) for _ in range(300)]
log2_xs = [float(x) for x in log2_xs]
for k in range(-1070, 1024, 37):
    p = 2.0 ** k
    log2_xs += [p * (1 + 2 ** -52), p * (1 - 2 ** -53), p * 1.5]
log2_xs += [rnd.uniform(0, 1) * 1e-310 for _ in range(30)] + [5e-324, 1e-320, 0.1, 0.5, 3.0, 10.0, 1e10]


def neighbours(x, k):
    out = []
    for i in range(-k, k + 1):
        y = x
        for _ in range(abs(i)):
            y = math.nextafter(y, math.inf if i > 0 else -math.inf)
        out.append(y)
    return out


# range-reduction boundaries: neighbours of 1 (cancellation), of sqrt(2) (x or x/2), of the smallest normal, the largest
log2_xs += neighbours(1.0, 40) + [1 + 1e-15, 1 + 1e-12, 1 + 1e-9, 1 - 1e-15, 1 - 1e-12, 1 - 1e-9]
log2_xs += [y * 2.0 ** e for e in (0, -3, 10, 500) for y in neighbours(math.sqrt(2), 5)]
log2_xs += neighbours(2.2250738585072014e-308, 5) + neighbours(1.7976931348623157e308, 0) + [math.nextafter(1.7976931348623157e308, 0)]


def log2_row(x):
    mp.mp.dps = 60
    return '[%s %s]' % (fmt(x), fmt(mp.log(mp.mpf(x), 2)))


# logexpm1: tiny to huge arguments, around the 18 / 33.3 / 36.8 / 709.78 points
logexpm1_xs = [1e-300, 1e-100, 1e-10, 1e-5, 0.01, 0.5, 1.0, 5.0, 17.9, 18.0, 18.1, 30.0, 33.3, 33.4, 36.7, 36.8, 36.9,
               100.0, 709.0, 709.78, 709.79, 710.0, 800.0, 1e5, 1e100, 1e300]


def logexpm1_row(x):
    mp.mp.dps = 60
    return '[%s %s]' % (fmt(x), fmt(mp.log(mp.expm1(mp.mpf(x)))))


# logsumexp: random vectors, ties, wide ranges
lse_cases = [[rnd.gauss(0, s) for _ in range(n)] for n in (1, 2, 3, 10, 100) for s in (1.0, 30.0, 300.0) for _ in range(2)]
lse_cases += [[1.0, 1.0], [3.0, 3.0, 3.0], [-1000.0, 5.0, 5.0], [1e308, 1e308], [-1e308, 1.0], [-745.0, -745.0],
              [700.0, 710.0, 720.0]]


def lse_row(xs):
    mp.mp.dps = 60
    return '[%s %s]' % (vec(xs), fmt(mp.log(mp.fsum(mp.exp(mp.mpf(x)) for x in xs))))


# smooth-max definitions (alpha != 0)
def hard(a, b, al):
    return max(a, b) if al > 0 else min(a, b)


def f_lse(xs, al):
    return mp.log(mp.fsum(mp.exp(al * x) for x in xs)) / al


def f_mellowmax(xs, al):
    m = max(xs) if al > 0 else min(xs)
    return m + mp.log(mp.fsum(mp.exp(al * (x - m)) for x in xs) / len(xs)) / al


def f_boltzmann(xs, al):
    m = max(al * x for x in xs)
    ws = [mp.exp(al * x - m) for x in xs]
    return mp.fsum(x * w for x, w in zip(xs, ws)) / mp.fsum(ws)


def f_p_norm(xs, al):
    return mp.fsum(abs(x) ** al for x in xs) ** (1 / al)


def fold(xs, pair):
    r = xs[0]
    for b in xs[1:]:
        r = pair(r, b)
    return r


def f_smu(xs, al):
    eps = abs(1 / al)
    s = 1 if al > 0 else -1
    return fold(xs, lambda a, b: (a + b + s * mp.sqrt((a - b) ** 2 + eps)) / 2)


def kernel_fold(g, c):
    def f(xs, al):
        k = mp.mpf(c) * (-1 / al)
        return fold(xs, lambda a, b: b - k * g((b - a) / k))
    return f


def g_sigmoid(x):
    return 1 / mp.ln2 if x == 0 else x / (1 - mp.mpf(2) ** (-x))


def g_circular(x):
    return x if x > 1 else (0 if x < -1 else 1 + (x - mp.sqrt(2 - x * x)) / 2)


def g_quadratic(x):
    return x if x > 1 else (0 if x < -1 else (1 + x * (2 + x)) / 4)


def g_cubic(x):
    return x if x > 1 else (0 if x < -1 else (1 + 3 * x * (1 + x) - abs(x ** 3)) / 6)


def g_quartic(x):
    return x if x > 1 else (0 if x < -1 else (1 + x) ** 2 * (3 - x * (x - 2)) / 16)


def f_exponential(xs, al):
    return fold(xs, lambda a, b: hard(a, b, al) + mp.log(1 + mp.mpf(2) ** (-abs(al * (a - b))), 2) / al)


families = [('lse', f_lse), ('boltzmann', f_boltzmann), ('mellowmax', f_mellowmax), ('p-norm', f_p_norm),
            ('smu', f_smu), ('sigmoid', kernel_fold(g_sigmoid, float(mp.log(2)))),
            ('circular', kernel_fold(g_circular, 3.414213562373096)), ('quadratic', kernel_fold(g_quadratic, 4.0)),
            ('cubic', kernel_fold(g_cubic, 6.0)), ('quartic', kernel_fold(g_quartic, 5.333333333333333)),
            ('exponential', f_exponential)]

sm_cases = [[rnd.gauss(0, s) for _ in range(n)] for n in (1, 2, 3, 10) for s in (0.01, 1.0, 30.0) for _ in range(2)]
sm_cases += [[1.0, 2.0, 3.0, 0.5], [1.0, 1.0], [3.0, 3.0, 1.0], [-2.0, -2.0], [400.0, 1.0], [-400.0, -500.0],
             [600.0, 1.0], [1e200, 1.0], [1e308, 1.0], [-1e300, 1.0], [1e-200, 1e-200], [1e-200, 2e-200],
             [1e-160, 1e-160], [1.0, 1.0000000001], [5.0, 5.0 + 1e-12, 5.0 - 1e-12],
             [1e308, 5e307], [-1e308, -2e307]]
alphas = [2.0, -2.0, 0.01, 50.0]
# smu with a tiny alpha: epsilon = 1/|alpha| is comparable to (a - b)^2 beyond 1e150
# mellowmax around the cumulant-expansion switch (|alpha (x - m)| = 1e-8) and at tiny and subnormal alpha
sm_extra = {'smu': [([1e153, 0.0], 1e-305), ([0.0, 1e153], 1e-305), ([1e153, 0.0], -1e-305), ([3e152, -1e152], 2e-305)],
            'mellowmax': [([1.0, 1.0 + 0.9e-8], 1.0), ([1.0, 1.0 + 1.1e-8], 1.0), ([0.3, 1.7, -2.2], 1e-9),
                          ([0.3, 1.7, -2.2], 5e-9), ([0.3, 1.7, -2.2], -5e-9), ([0.3, 1.7, -2.2], 1e-300),
                          ([1.0, 2.0, 3.0, 0.5], 5e-324), ([1.0, 2.0, 3.0, 0.5], -5e-324)]}


def sm_rows(name, f):
    rows = []
    for xs in sm_cases:
        for al in alphas:
            if name == 'p-norm' and any(x == 0 for x in xs):
                continue
            mp.mp.dps = precision(xs)
            y = f([mp.mpf(x) for x in xs], mp.mpf(al))
            rows.append('[%s %s %s]' % (vec(xs), fmt(al), fmt(y)))
    for xs, al in sm_extra.get(name, []):
        mp.mp.dps = 700
        y = f([mp.mpf(x) for x in xs], mp.mpf(al))
        rows.append('[%s %s %s]' % (vec(xs), fmt(al), fmt(y)))
    return rows


out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..',
                   'test', 'resources', 'core', 'logexp_reference.edn')
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, 'w') as f:
    f.write('{:log2 [\n  %s]\n' % '\n  '.join(log2_row(x) for x in log2_xs))
    f.write(' :logexpm1 [\n  %s]\n' % '\n  '.join(logexpm1_row(x) for x in logexpm1_xs))
    f.write(' :logsumexp [\n  %s]\n' % '\n  '.join(lse_row(xs) for xs in lse_cases))
    f.write(' :smooth-max {%s}}\n' % '\n  '.join(':%s [\n   %s]' % (name, '\n   '.join(sm_rows(name, g)))
                                                    for name, g in families))
