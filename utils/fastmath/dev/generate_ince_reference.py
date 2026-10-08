"""One-off generator of independent reference values for the Ince polynomials of `fastmath.polynomials`:
`ince-C-coeffs`, `ince-S-coeffs`, `ince-C`, `ince-S`, `ince-C-radial`, `ince-S-radial`.

Produces `test/resources/polynomials/ince_reference.edn`. The eigenproblem is built from the differential
equation itself, not from the coefficient recurrences the library uses (which are the DLMF 28.31 matrices):

  Ince's equation   w'' + e sin(2 x) w' + (a - p e cos(2 x)) w = 0

The Fourier-type basis of each of the four cases (C or S, even or odd p) is

  C, p even: cos(2 r x),        r = 0 .. p/2            S, p even: sin(2 (r+1) x),  r = 0 .. p/2 - 1
  C, p odd:  cos((2 r + 1) x),  r = 0 .. (p-1)/2        S, p odd:  sin((2 r + 1) x), r = 0 .. (p-1)/2

For every basis function b_r the function L[b_r] = b_r'' + e sin(2x) b_r' - p e cos(2x) b_r is a combination of
the basis; its coefficients (the matrix of -L) are found by solving the interpolation system at distinct
points in (0, pi/2) in `mpmath` (60 digits), then the eigenvalues `a` of that matrix are sorted increasingly
and the one with index `m/2` (C, p even), `(m-1)/2` (C, p odd), `m/2 - 1` (S, p even) or `(m-1)/2` (S, p odd)
is the one of the polynomial of degree `m`. The eigenvectors are normalised to the Euclidean norm 1 and the
sign is chosen so that C(0) > 0 and S'(0) > 0. The generator asserts that the series satisfies the equation at
random points to 1e-45 and that `a = m^2` for `e = 0`.

Layout: `[{:kind :C|:S :p p :m m :e e :a a :coefficients [c0 c1 ...]} ...]`, coefficients in the basis order above.

A second file, `ince_large_e_reference.edn`, has the same layout plus `:zero-value` (C(0) for C, S'(0) for S, the
quantity whose sign fixes the sign of the vector) for large `|e|` (up to 1e8) and `p` up to 30, where a
non-symmetric eigen solver loses the eigenvectors; the sign of a vector is only meaningful where `|:zero-value|`
is above rounding.

Run with (from repo root, using the `uv`-managed Python env mentioned in AGENTS.md):
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_ince_reference.py
"""
import random
import mpmath as mp

mp.mp.dps = 60
OUT = "/home/ts/clojure/fastmath/test/resources/polynomials/ince_reference.edn"


def basis(kind, p):
    """List of (function, first derivative, second derivative) of the basis, as lambdas of mp numbers."""
    if kind == "C":
        freqs = [2 * r for r in range(p // 2 + 1)] if p % 2 == 0 else [2 * r + 1 for r in range((p - 1) // 2 + 1)]
        return [((lambda x, k=k: mp.cos(k * x)), (lambda x, k=k: -k * mp.sin(k * x)), (lambda x, k=k: -k * k * mp.cos(k * x)))
                for k in freqs]
    freqs = [2 * (r + 1) for r in range(p // 2)] if p % 2 == 0 else [2 * r + 1 for r in range((p - 1) // 2 + 1)]
    return [((lambda x, k=k: mp.sin(k * x)), (lambda x, k=k: k * mp.cos(k * x)), (lambda x, k=k: -k * k * mp.sin(k * x)))
            for k in freqs]


def operator(kind, p, e, b, x):
    """L[b](x) = b'' + e sin(2x) b' - p e cos(2x) b."""
    f, d1, d2 = b
    return d2(x) + e * mp.sin(2 * x) * d1(x) - p * e * mp.cos(2 * x) * f(x)


def matrix(kind, p, e):
    bs = basis(kind, p)
    size = len(bs)
    points = [(j + mp.mpf(1) / 2) * mp.pi / (2 * size) for j in range(size)]
    b_matrix = mp.matrix(size, size)
    for j, x in enumerate(points):
        for k, b in enumerate(bs):
            b_matrix[j, k] = b[0](x)
    m_matrix = mp.matrix(size, size)
    for r, b in enumerate(bs):
        rhs = mp.matrix([-operator(kind, p, e, b, x) for x in points])  # -L[b_r]
        column = mp.lu_solve(b_matrix, rhs)
        for k in range(size):
            m_matrix[k, r] = column[k]
    return m_matrix, bs


def eigen(kind, p, e):
    m_matrix, bs = matrix(kind, p, e)
    values, vectors = mp.eig(m_matrix)
    order = sorted(range(len(values)), key=lambda i: mp.re(values[i]))
    result = []
    for i in order:
        assert abs(mp.im(values[i])) < mp.mpf(10) ** -40, ("complex eigenvalue", kind, p, e, values[i])
        v = [mp.re(vectors[k, i]) for k in range(len(values))]
        norm = mp.sqrt(sum(c * c for c in v))
        v = [c / norm for c in v]
        # sign: C(0) > 0, S'(0) > 0
        if kind == "C":
            at_zero = sum(c * b[0](mp.mpf(0)) for c, b in zip(v, bs))
        else:
            at_zero = sum(c * b[1](mp.mpf(0)) for c, b in zip(v, bs))
        if at_zero < 0:
            v = [-c for c in v]
        result.append((mp.re(values[i]), v))
    return result, bs


def index_of(kind, p, m):
    if kind == "C":
        return m // 2 if p % 2 == 0 else (m - 1) // 2
    return m // 2 - 1 if p % 2 == 0 else (m - 1) // 2


def fmt(v):
    return repr(float(v))


random.seed(7)
entries = []
checked = 0
for e_value in (0.0, 0.3, 0.6, -0.6, 2.0, 5.0, 10.0):
    e = mp.mpf(e_value)
    for p in list(range(0, 9)) + [12]:
        for kind in ("C", "S"):
            if kind == "S" and p == 0:
                continue
            eigs, bs = eigen(kind, p, e)
            for m in range(0 if kind == "C" else 1, p + 1):
                if (p - m) % 2:
                    continue
                a, v = eigs[index_of(kind, p, m)]
                if e_value == 0.0:
                    assert abs(a - m * m) < mp.mpf(10) ** -45, (kind, p, m, a)
                for _ in range(3):
                    x = mp.mpf(random.random()) * 3
                    w = sum(c * b[0](x) for c, b in zip(v, bs))
                    w1 = sum(c * b[1](x) for c, b in zip(v, bs))
                    w2 = sum(c * b[2](x) for c, b in zip(v, bs))
                    residual = w2 + e * mp.sin(2 * x) * w1 + (a - p * e * mp.cos(2 * x)) * w
                    assert abs(residual) < mp.mpf(10) ** -45, (kind, p, m, e_value, residual)
                    checked += 1
                entries.append("  {:kind :%s :p %d :m %d :e %s :a %s\n   :coefficients [%s]}" % (
                    kind, p, m, fmt(e), fmt(a), " ".join(fmt(c) for c in v)))

with open(OUT, "w") as f:
    f.write(";; Generated by utils/fastmath/dev/generate_ince_reference.py (mpmath, 60 digits, built from the differential equation)\n")
    f.write("[\n" + "\n".join(entries) + "\n]\n")

print("entries:", len(entries), "equation checks:", checked)
print("wrote", OUT)


# ---- large |e| ---------------------------------------------------------------------------------------------
LARGE_OUT = "/home/ts/clojure/fastmath/test/resources/polynomials/ince_large_e_reference.edn"
large_entries = []
large_checked = 0
for e_value in (-2500.0, -300.0, 300.0, 1000.0, 2500.0, 1e5, 1e8):
    e = mp.mpf(e_value)
    for p in (4, 9, 12, 20, 30):
        for kind in ("C", "S"):
            eigs, bs = eigen(kind, p, e)
            for m in range(0 if kind == "C" else 1, p + 1):
                if (p - m) % 2:
                    continue
                a, v = eigs[index_of(kind, p, m)]
                for _ in range(2):
                    x = mp.mpf(random.random()) * 3
                    w = sum(c * b[0](x) for c, b in zip(v, bs))
                    w1 = sum(c * b[1](x) for c, b in zip(v, bs))
                    w2 = sum(c * b[2](x) for c, b in zip(v, bs))
                    residual = w2 + e * mp.sin(2 * x) * w1 + (a - p * e * mp.cos(2 * x)) * w
                    assert abs(residual) < mp.mpf(10) ** -35 * (1 + abs(e) * p), (kind, p, m, e_value, residual)
                    large_checked += 1
                index = 0 if kind == "C" else 1
                zero_value = sum(c * b[index](mp.mpf(0)) for c, b in zip(v, bs))
                large_entries.append("  {:kind :%s :p %d :m %d :e %s :a %s :zero-value %s\n   :coefficients [%s]}" % (
                    kind, p, m, fmt(e), fmt(a), fmt(zero_value), " ".join(fmt(c) for c in v)))

with open(LARGE_OUT, "w") as f:
    f.write(";; Generated by utils/fastmath/dev/generate_ince_reference.py (mpmath, 60 digits): large |e|\n")
    f.write("[\n" + "\n".join(large_entries) + "\n]\n")
print("large e entries:", len(large_entries), "equation checks:", large_checked)
print("wrote", LARGE_OUT)
