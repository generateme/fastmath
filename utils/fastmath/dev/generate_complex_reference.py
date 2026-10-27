# Reference values for fastmath.complex (test/resources/complex/complex_reference.edn).
#
# Tool: Python mpmath 1.3.0, 60 digits plus twice the number of leading zeros of the smallest part of z and 1/z (/home/ts/penv, `uv run python`).
# Domain: points a+bi with a, b from a list of magnitudes (1e-300 .. 1e300, around 1, around the
# overflow of cosh at 710) and both signs, with a != 0 and b != 0 (the branch cuts and signed zeros are
# covered by the Julia grid, generate_complex_cuts_reference.jl, because mpmath has no signed zero).
# Output: one map, :points [[re im] ...] and for every function a vector of [re im] (nil when mpmath
# fails), for `pow` pairs of arguments. Infinities are written as ##Inf.
# Regeneration is byte-identical (fixed order, repr of doubles).
import mpmath as mp
import math
import sys

mp.mp.dps = 60

mags = [1e-300, 1e-200, 1e-100, 1e-10, 1e-5, 0.1, 0.5, 0.9, 0.99999, 1.0, 1.00001, 2.0, 10.0, 100.0,
        710.0, 800.0, 1e5, 1e10, 1e100, 1e200, 1e300]
vals = sorted(set([-m for m in mags] + mags))
points = [(a, b) for a in vals for b in vals]
# huge arguments of trigonometric functions are expensive and overflow anyway: thin them out
points = [(a, b) for (a, b) in points if not (abs(a) > 1e10 and abs(b) > 1e10 and (abs(a) < 1e200 or abs(b) < 1e200))]

fns = {
    'sqrt': mp.sqrt, 'exp': mp.exp, 'log': mp.log,
    'sin': mp.sin, 'cos': mp.cos, 'tan': mp.tan,
    'sec': mp.sec, 'csc': mp.csc, 'cot': mp.cot,
    'sinh': mp.sinh, 'cosh': mp.cosh, 'tanh': mp.tanh,
    'sech': mp.sech, 'csch': mp.csch, 'coth': mp.coth,
    'asin': mp.asin, 'acos': mp.acos, 'atan': mp.atan,
    'asec': mp.asec, 'acsc': mp.acsc, 'acot': mp.acot,
    'asinh': mp.asinh, 'acosh': mp.acosh, 'atanh': mp.atanh,
    'asech': lambda z: mp.acosh(1 / z), 'acsch': lambda z: mp.asinh(1 / z), 'acoth': lambda z: mp.atanh(1 / z),
    'log1p': mp.log1p, 'expm1': mp.expm1, 'cbrt': mp.cbrt,
    'nth-root-5': lambda z: mp.root(z, 5),
}

pow_bases = [(a, b) for a in [-100.0, -2.0, -0.5, -1e-5, 1e-5, 0.5, 2.0, 100.0, 1e50]
             for b in [-3.0, -0.5, 1e-5, 0.5, 3.0, 1e50]]
pow_exps = [(2.0, 0.0), (3.0, 0.0), (-2.0, 0.0), (0.5, 0.0), (-1.5, 0.0), (10.3, 0.0), (1e-5, 0.0), (0.0, 1.0),
            (0.0, -2.5), (1.5, -2.0), (-0.3, 0.7), (40.0, 0.0), (-33.0, 0.0), (1.0, 1e-5)]


def fmt(x):
    if math.isinf(x):
        return '##Inf' if x > 0 else '##-Inf'
    if math.isnan(x):
        return '##NaN'
    return repr(x)


def tofloat(x):
    try:
        return float(x)
    except OverflowError:
        return math.inf if x > 0 else -math.inf


def row(w):
    if w is None:
        return 'nil'
    return '[' + fmt(tofloat(w.real)) + ' ' + fmt(tofloat(w.imag)) + ']'


def digits_for(*zs):
    # mpmath loses the small component of asin, acos, atan... near zero, near the branch points and for
    # reciprocals of huge numbers when the working precision is too low: the precision follows the smallest of
    # the parts, 1/(the largest part) and the smallest part of 1/z (min/max^2)
    parts = [abs(float(t)) for z in zs for t in (z.real, z.imag) if t != 0]
    mn = min(parts + [1.0])
    mx = max(parts + [1.0])
    small = min(mn, 1.0 / mx, mn / (mx * mx), 1.0)
    return 60 + 2 * int(math.ceil(-math.log10(small)))


def apply(f, *args):
    try:
        with mp.workdps(digits_for(*args)):
            w = f(*args)
            if isinstance(w, mp.mpf):
                w = mp.mpc(w, 0)
            return w
    except Exception as e:
        return None


out = sys.stdout
out.write('{\n')
out.write(':points [' + ' '.join('[' + fmt(a) + ' ' + fmt(b) + ']' for a, b in points) + ']\n')
for name, f in fns.items():
    out.write(':' + name + ' [' + ' '.join(row(apply(f, mp.mpc(a, b))) for a, b in points) + ']\n')
out.write(':pow-args [' + ' '.join('[[' + fmt(a) + ' ' + fmt(b) + '] [' + fmt(c) + ' ' + fmt(d) + ']]'
                                    for (a, b) in pow_bases for (c, d) in pow_exps) + ']\n')


def pow_row(a, b, c, d):
    # an overflowing power has no meaningful direction (the sign of the infinite parts follows sin and cos of
    # the huge argument): such points are skipped
    w = apply(lambda z, w: z ** w, mp.mpc(a, b), mp.mpc(c, d))
    if w is not None and (abs(w.real) > 1e308 or abs(w.imag) > 1e308):
        return None
    return w


out.write(':pow [' + ' '.join(row(pow_row(a, b, c, d))
                              for (a, b) in pow_bases for (c, d) in pow_exps) + ']\n')
out.write('}\n')
