"""Generator of reference values for `fastmath.kernel.rbf` (radial basis function families).

Tools: sympy (exact rational arithmetic for the compactly supported families), mpmath (50 digits,
for the Bessel, Whittaker-type and shifted surface spline families) and, through
`generate_rbf_kernel_reference.R`, the R package fields (`Wendland`). Every number is converted to the
nearest double (`float` of an exact rational or of an mpmath number, correctly rounded) and written
with `repr`. The distances are doubles; the oracle is evaluated at exactly that double.

Domain and derivation (each family is derived independently from the literature, not from the code):

  :wu (Wu 1995, Fasshauer MATH 590 chapter 11)
      psi_l(r) = integral over t from 2r - 1 to 1 of (1 - t^2)^l (1 - (2r - t)^2)^l (the convolution
      (1 - .^2)^l * (1 - .^2)^l at 2r), psi_{l,k} = D^k psi_l with D f(r) = -f'(r) / r, normalised by
      psi_{l,k}(0) = 1. Exact for integer l (psi_l is a polynomial of degree 4l + 1 on [0, 1]): l 0..4,
      k 0..l (k <= l makes D^k a polynomial). A non-integer l (1.5, k 0) uses mpmath quadrature.
  :wendland (Wendland 1995, Fasshauer Definition 4.4)
      phi_{s,k} = I^k (1 - r)_+^l with l = floor(s / 2) + k + 1 and I f(r) = integral from r to 1 of
      t f(t) dt, normalised by phi(0) = 1; s 1..8, k 0..4, exact rationals.
  :gneiting (Gneiting 2002, Fasshauer chapter 11)
      the turning bands operator phi_s = phi_{s+2} + r phi'_{s+2} / s applied to the Wendland function
      (1 - r)_+^{l+1} ((l + 1) r + 1); s 1..4, l in (3.5, 4, 5.5, 7), exact rationals.
  :shifted-surface-splines (Luh, arXiv:1011.0124, equations (1) and (2) with beta in place of lambda)
      t = beta - s / 2; odd s: (-1)^ceil(t) (1 + r^2)^t; even s (t a positive integer):
      (-1)^(1 + t) (1 + r^2)^t log sqrt(1 + r^2); s 1..5.
  :poisson
      Gamma(p + 1) 2^p J_p(r) / r^p with p = d / 2 - 1 (value 1 at r = 0), d in (1, 1.5, 2 .. 7).
  :matern
      2^(1 - v) / Gamma(v) r^v K_v(r) with v = order / 2 (value 1 at r = 0), order 1, 3, 5, 7, 9.
  :whittaker
      the integral from 0 to 1 / r of (1 - r t)^(k - 1) t^alpha exp(-beta t) dt normalised by its value at
      r = 0 (Gamma(alpha + 1) / beta^(alpha + 1)); alpha 0, 1, 0.5; k 2, 3; beta 0.5, 1, 2.5.
  :gaussians-laguerre
      exp(-r^2) L_n^(d/2)(r^2) (unnormalised, L_n^(d/2)(0) at r = 0); d 1, 2, 3, 2.5; n 0..4 and 20 (the recurrence must not suffer from the cancellation of the monomial coefficients).
  :fields-wendland
      the R package fields, see the .R generator (k 1..4, dimension 1..8).

Why: the old code compared boxed vectors in `wu`/`whittaker` (the closed forms were never reached, wu-33
at r = 1e-3 gave -20958 instead of 0.9978), used a truncated exponent in the shifted surface splines, did
not normalise poisson, and the Wendland exponent differed from the textbook for odd s.

Produces `test/resources/kernel/rbf_kernel_reference.edn`: a map keyed by family; each family is a vector of
cases `{<parameters> :r [...] :v [...]}` with parallel vectors (r = distance, v = value, scale 1).
Deterministic (no randomness), regenerates byte-identically.

Run from anywhere:
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_rbf_kernel_reference.py
"""
import subprocess
from fractions import Fraction
from pathlib import Path

import mpmath
import sympy
from sympy import Rational, integrate, symbols, cancel, diff

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = ROOT / "test" / "resources" / "kernel" / "rbf_kernel_reference.edn"
RSCRIPT = HERE / "generate_rbf_kernel_reference.R"

mpmath.mp.dps = 50
mp = mpmath.mp

r_sym, t_sym = symbols("r t")

# distances (doubles) shared by the compactly supported families
R_COMPACT = [0.0, 1e-3, 1e-2, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 1.0, 1.5]
# distances for the globally supported families
R_GLOBAL = [0.0, 1e-9, 1e-3, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 25.0]
R_SPECIAL = [0.0, 1e-9, 1e-5, 1e-3, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]


def edn_double(x):
    s = repr(float(x))
    assert "inf" not in s and "nan" not in s, s
    return s


def vec(xs):
    return "[" + " ".join(edn_double(x) for x in xs) + "]"


def to_double(r):
    return float(Fraction(int(r.p), int(r.q)))


def exact(x):
    """The double `x` as an exact sympy rational."""
    f = Fraction(x)
    return Rational(f.numerator, f.denominator)


def case(params, rs, vs):
    ps = " ".join(":%s %s" % (k, v) for k, v in params)
    return " {%s :r %s :v %s}" % (ps, vec(rs), vec(vs))


def poly_value(poly, r):
    v = poly.subs(r_sym, exact(r))
    if v.is_Rational:
        return to_double(v)
    return float(v.evalf(50))


# ---------------------------------------------------------------- wu
def wu_poly(l, k):
    base = integrate(
        (1 - t_sym**2) ** l * (1 - (2 * r_sym - t_sym) ** 2) ** l,
        (t_sym, 2 * r_sym - 1, 1),
    )
    f = sympy.expand(base)
    for _ in range(k):
        f = cancel(-diff(f, r_sym) / r_sym)
    f = sympy.expand(f)
    f0 = f.subs(r_sym, 0)
    assert f0 != 0
    return sympy.expand(f / f0)


def wu_cases():
    rows = []
    for l in range(0, 5):
        for k in range(0, l + 1):
            poly = wu_poly(l, k)
            vs = [poly_value(poly, r) if r <= 1.0 else 0.0 for r in R_COMPACT]
            rows.append(case([("l", l), ("k", k)], R_COMPACT, vs))
    # a non-integer l: base function by quadrature
    def psi(l, r):
        if r >= 1:
            return mp.mpf(0)
        r = mp.mpf(r)
        f = lambda t: (1 - t * t) ** l * (1 - (2 * r - t) ** 2) ** l
        return mp.quad(f, [2 * r - 1, r, 1])

    l = mp.mpf("1.5")
    p0 = psi(l, 0.0)
    vs = [float(psi(l, r) / p0) for r in R_COMPACT]
    rows.append(case([("l", "1.5"), ("k", 0)], R_COMPACT, vs))
    return rows


# ---------------------------------------------------------------- wendland
def wendland_poly(s, k):
    l = s // 2 + k + 1
    f = (1 - r_sym) ** l
    for _ in range(k):
        f = integrate(t_sym * f.subs(r_sym, t_sym), (t_sym, r_sym, 1))
    f = sympy.expand(f)
    return sympy.expand(f / f.subs(r_sym, 0))


def wendland_cases():
    rows = []
    for s in range(1, 9):
        for k in range(0, 5):
            poly = wendland_poly(s, k)
            vs = [poly_value(poly, r) if r <= 1.0 else 0.0 for r in R_COMPACT]
            rows.append(case([("s", s), ("k", k)], R_COMPACT, vs))
    return rows


# ---------------------------------------------------------------- gneiting
def gneiting_cases():
    rows = []
    for s in range(1, 5):
        for l_str in ("7/2", "4", "11/2", "7"):
            l = Rational(l_str)
            phi = (1 - r_sym) ** (l + 1) * ((l + 1) * r_sym + 1)
            tau = sympy.expand(sympy.simplify(phi + r_sym * diff(phi, r_sym) / s))
            vs = [poly_value(tau, r) if r <= 1.0 else 0.0 for r in R_COMPACT]
            rows.append(case([("s", s), ("l", edn_double(float(l)))], R_COMPACT, vs))
    return rows


# ---------------------------------------------------------------- shifted surface splines
def shifted_cases():
    rows = []
    rs = [0.0, 1e-3, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    spec = {
        1: [1.0, 1.5, 2.0, 3.5],
        2: [2.0, 3.0, 4.0],
        3: [2.0, 2.5, 3.0, 4.5],
        4: [3.0, 4.0, 5.0],
        5: [3.0, 3.5, 4.0, 5.5],
    }
    for s, betas in spec.items():
        for beta in betas:
            t = mp.mpf(beta) - mp.mpf(s) / 2
            assert t > 0
            vs = []
            for r in rs:
                u = 1 + mp.mpf(r) ** 2
                if s % 2 == 1:
                    sign = -1 if int(mp.ceil(t)) % 2 else 1
                    vs.append(float(sign * u**t))
                else:
                    assert t == int(t)
                    sign = -1 if (1 + int(t)) % 2 else 1
                    vs.append(float(sign * u**t * mp.log(mp.sqrt(u))))
            rows.append(case([("s", s), ("beta", edn_double(beta))], rs, vs))
    return rows


# ---------------------------------------------------------------- poisson
def poisson_cases():
    rows = []
    for d in (1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0):
        p = mp.mpf(d) / 2 - 1
        vs = []
        for r in R_SPECIAL:
            if r == 0.0:
                vs.append(1.0)
            else:
                x = mp.mpf(r)
                vs.append(float(mp.gamma(p + 1) * 2**p * mp.besselj(p, x) / x**p))
        rows.append(case([("d", edn_double(d))], R_SPECIAL, vs))
    return rows


# ---------------------------------------------------------------- matern
def matern_cases():
    rows = []
    rs = R_SPECIAL + [200.0, 700.0, 745.0, 800.0]
    for order in (1, 3, 5, 7, 9):
        v = mp.mpf(order) / 2
        vs = []
        for r in rs:
            if r == 0.0:
                vs.append(1.0)
            else:
                x = mp.mpf(r)
                vs.append(float(2 ** (1 - v) / mp.gamma(v) * x**v * mp.besselk(v, x)))
        rows.append(case([("order", order)], rs, vs))
    return rows


# ---------------------------------------------------------------- whittaker
def whittaker_cases():
    rows = []
    rs = [0.0, 1e-3, 0.1, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0, 1000.0, 1e5]
    for alpha in ("0", "1", "0.5"):
        for k in (2, 3):
            for beta in ("0.5", "1", "2.5"):
                a, b = mp.mpf(alpha), mp.mpf(beta)
                norm = mp.gamma(a + 1) / b ** (a + 1)
                vs = []
                for r in rs:
                    if r == 0.0:
                        vs.append(1.0)
                    else:
                        x = mp.mpf(r)
                        f = lambda t: (1 - x * t) ** (k - 1) * t**a * mp.exp(-b * t)
                        vs.append(float(mp.quad(f, [0, 1 / (2 * x), 1 / x]) / norm))
                rows.append(case([("alpha", alpha), ("k", k), ("beta", beta)], rs, vs))
    return rows


# ---------------------------------------------------------------- gaussians-laguerre
def laguerre_cases():
    rows = []
    for d in ("1", "2", "3", "2.5"):
        for n in (0, 1, 2, 3, 4, 20):
            alpha = mp.mpf(d) / 2
            vs = [float(mp.exp(-mp.mpf(r) ** 2) * mp.laguerre(n, alpha, mp.mpf(r) ** 2)) for r in R_GLOBAL]
            rows.append(case([("dimension", d), ("degree", n)], R_GLOBAL, vs))
    return rows


def fields_cases():
    res = subprocess.run(["Rscript", str(RSCRIPT)], capture_output=True, text=True, check=True)
    lines = [ln for ln in res.stdout.splitlines() if ln.startswith(" {")]
    assert len(lines) == 32, len(lines)
    return lines


def block(name, rows):
    return " :%s [\n%s]" % (name, "\n".join(rows) + "\n")


def main():
    header = ";; rbf_kernel_reference.edn - generated by utils/fastmath/dev/generate_rbf_kernel_reference.py\n"
    header += ";; sympy %s, mpmath %s, R fields (see the generator)\n" % (sympy.__version__, mpmath.__version__)
    parts = [
        block("wu", wu_cases()),
        block("wendland", wendland_cases()),
        block("gneiting", gneiting_cases()),
        block("shifted-surface-splines", shifted_cases()),
        block("poisson", poisson_cases()),
        block("matern", matern_cases()),
        block("whittaker", whittaker_cases()),
        block("gaussians-laguerre", laguerre_cases()),
        block("fields-wendland", fields_cases()),
    ]
    text = header + "{\n" + "\n".join(parts) + "}\n"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(text)
    print("wrote", OUT, len(text), "bytes")


if __name__ == "__main__":
    main()
