"""Generator of reference values for the vector kernels of `fastmath.kernel` / `fastmath.kernel.vector`.

Tools: scikit-learn (`rbf_kernel`, `polynomial_kernel`, `sigmoid_kernel`, `laplacian_kernel`,
`additive_chi2_kernel`, `gaussian_process.kernels.{Matern,RationalQuadratic,ExpSineSquared}`),
scipy (`BSpline.basis_element`), mpmath (80 digits: quadrature, Bessel J, sech, sin), numpy. The
R package kernlab is used by `generate_kernel_vector_kernlab_reference.R` (another file).

Domain: one case per (kernel, parameter set) in the registry, each with deterministic points
(`numpy.random.RandomState(20260514)`, fixed order of draws): 6 random pairs, an identical pair, a pair
at distance about 1e-9 and a far pair (offset 20). Nothing is stored for which the oracle is not
independent of the formula in `vector.clj`: a mapping to a library kernel is used where one exists
(gaussian, exponential and laplacian via Matern nu 1/2, laplacian L1, polynomial, sigmoid,
rational-quadratic, cauchy, matern odd orders, periodic, chi-square mappings), a direct textbook formula
in mpmath otherwise (wave, bessel, bessel2, dirichlet, hyperbolic-secant, pearson as in Weka `Puk`).
The mappings:
    gaussian(sigma)              rbf_kernel(gamma = 1 / (2 sigma^2))
    exponential(sigma)           Matern(length_scale = 2 sigma^2, nu = 1/2)
    laplacian(sigma)             Matern(length_scale = sigma, nu = 1/2);  :abs distance: laplacian_kernel(gamma = 1/sigma)
    rational-quadratic(c)        RationalQuadratic(length_scale = sqrt(c/2), alpha = 1)
    cauchy(sigma)                RationalQuadratic(length_scale = sigma/sqrt(2), alpha = 1)
    matern(order, theta)         Matern(length_scale = theta / sqrt(2), nu = order / 2)
    periodic(sigma, periodicity) ExpSineSquared(length_scale = sigma, periodicity)
    chi-square                   sum (x + y) / 2 + additive_chi2_kernel / 2
    chi-square2                  1 + 2 additive_chi2_kernel
The `:geometric` cases compare with the volume of the intersection of two n-balls of radius r at
centre distance d divided by the ball volume, computed by mpmath quadrature of
2 V_{n-1} / V_n * integral_{d/2r}^{1} (1 - t^2)^((n-1)/2) dt (and cross-checked with the regularised
incomplete beta function), n = 1..7, 9, 12 and 20; the `:b-spline` cases with scipy's cardinal B-spline of degree
2n+1 (zero outside the support).
Why: the recursion seed of the geometric kernel was twice the right value (31x off for n = 4 near the
support edge), pearson used a square root where the Ustun / Weka form has none, chi-square gave NaN at 0/0,
b-spline returned noise outside its support.

Produces `test/resources/kernel/kernel_vector_reference.edn`: a map `{:cases [...]}` where a case is
`{:kernel :params :oracle :x [[...]] :y [[...]] :expected [...]}` (point i of `:x` pairs with point i of
`:y`; `:distance :abs` in `:params` is the keyword of `fastmath.vector/distances`).

Run from anywhere:
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_kernel_vector_reference.py
"""
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
import sklearn
from scipy.interpolate import BSpline
from sklearn.gaussian_process.kernels import ExpSineSquared, Matern, RationalQuadratic
from sklearn.metrics.pairwise import (
    additive_chi2_kernel,
    laplacian_kernel,
    polynomial_kernel,
    rbf_kernel,
    sigmoid_kernel,
)

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "test" / "resources" / "kernel" / "kernel_vector_reference.edn"

mp.mp.dps = 80
rng = np.random.RandomState(20260514)
CASES = []


def pairs(d, lo=-1.5, hi=1.5, nonneg=False, zeros=False):
    """6 random pairs, an identical pair, a pair at distance ~1e-9, a far pair."""
    if nonneg:
        lo, hi = 0.0, 2.0
    x = rng.uniform(lo, hi, (6, d))
    y = rng.uniform(lo, hi, (6, d))
    if zeros:
        x[1, 0] = 0.0
        y[1, 0] = 0.0
        x[2, :] = 0.0
        y[2, :] = 0.0
        y[3, 0] = 0.0
    x = np.vstack([x, x[0:1], x[0:1], x[1:2]])
    y = np.vstack([y, x[0:1], x[0:1] + 1e-9, x[1:2] + 20.0])
    return x, y


def diag(kfun, x, y):
    return [float(kfun(x[i : i + 1], y[i : i + 1])[0, 0]) for i in range(len(x))]


def add(kernel, params, x, y, expected, oracle):
    CASES.append(
        dict(
            kernel=kernel,
            params=params,
            oracle=oracle,
            x=x.tolist(),
            y=y.tolist(),
            expected=[float(e) for e in expected],
        )
    )


def euclid(x, y):
    return np.sqrt(np.sum((x - y) ** 2, axis=1))


def mpf_dist(x, y):
    return mp.sqrt(sum((mp.mpf(float(a)) - mp.mpf(float(b))) ** 2 for a, b in zip(x, y)))


# --- sklearn mappings -----------------------------------------------------------------------------
x, y = pairs(3)
for sigma in (0.5, 1.0, 2.3):
    add("gaussian", dict(sigma=sigma), x, y,
        diag(lambda a, b: rbf_kernel(a, b, gamma=1.0 / (2 * sigma**2)), x, y), "sklearn rbf_kernel")
x, y = pairs(2)
for sigma in (0.4, 1.0, 1.7):
    add("exponential", dict(sigma=sigma), x, y,
        diag(Matern(length_scale=2 * sigma**2, nu=0.5), x, y), "sklearn Matern nu=1/2, length_scale=2 sigma^2")
for sigma in (0.3, 1.0, 2.5):
    add("laplacian", dict(sigma=sigma), x, y,
        diag(Matern(length_scale=sigma, nu=0.5), x, y), "sklearn Matern nu=1/2, length_scale=sigma")
x, y = pairs(3)
for sigma in (0.3, 1.0, 2.5):
    add("laplacian", dict(sigma=sigma, distance="abs"), x, y,
        diag(lambda a, b: laplacian_kernel(a, b, gamma=1.0 / sigma), x, y), "sklearn laplacian_kernel (L1)")
for alpha, c, p in ((1.0, 0.0, 2.0), (0.7, 1.0, 3.0), (2.0, 0.5, 2.0)):
    add("polynomial", dict(alpha=alpha, c=c, p=p), x, y,
        diag(lambda a, b: polynomial_kernel(a, b, degree=p, gamma=alpha, coef0=c), x, y), "sklearn polynomial_kernel")
x, y = pairs(3, lo=-0.8, hi=0.8)
for alpha, c in ((1.0, 0.0), (0.5, 0.3), (2.0, -0.2)):
    add("hyperbolic-tangent", dict(alpha=alpha, c=c), x, y,
        diag(lambda a, b: sigmoid_kernel(a, b, gamma=alpha, coef0=c), x, y), "sklearn sigmoid_kernel")
x, y = pairs(3)
add("linear", {}, x, y, [float(np.dot(a, b)) for a, b in zip(x, y)], "numpy dot")
for c in (0.5, 1.0, 4.0):
    add("rational-quadratic", dict(c=c), x, y,
        diag(RationalQuadratic(length_scale=np.sqrt(c / 2), alpha=1.0), x, y),
        "sklearn RationalQuadratic alpha=1, length_scale=sqrt(c/2)")
for sigma in (0.5, 1.0, 3.0):
    add("cauchy", dict(sigma=sigma), x, y,
        diag(RationalQuadratic(length_scale=sigma / np.sqrt(2), alpha=1.0), x, y),
        "sklearn RationalQuadratic alpha=1, length_scale=sigma/sqrt(2)")
x, y = pairs(2)
for kw, order, theta in (("matern-12", 1, 1.0), ("matern-32", 3, 1.0), ("matern-52", 5, 1.0),
                         ("matern", 1, 2.0), ("matern", 3, 0.7), ("matern", 5, 3.1),
                         ("matern", 7, 1.3), ("matern", 9, 2.0), ("matern", 11, 0.9)):
    params = {} if kw != "matern" else dict(order=order, theta=theta)
    add(kw, params, x, y,
        diag(Matern(length_scale=theta / np.sqrt(2), nu=order / 2.0), x, y),
        "sklearn Matern nu=order/2, length_scale=theta/sqrt(2)")
x, y = pairs(1, lo=-3.0, hi=3.0)
for sigma, per in ((1.0, 1.0), (0.6, 2.5), (2.0, 0.7)):
    add("periodic", dict(sigma=sigma, periodicity=per), x, y,
        diag(ExpSineSquared(length_scale=sigma, periodicity=per), x, y), "sklearn ExpSineSquared")

# chi-square: sklearn additive_chi2_kernel A = -sum (x-y)^2 / (x+y), 0/0 = 0
x, y = pairs(4, nonneg=True, zeros=True)
A = np.array(diag(additive_chi2_kernel, x, y))
add("chi-square", {}, x, y, np.sum(x + y, axis=1) / 2 + A / 2,
    "sum(x+y)/2 + sklearn additive_chi2_kernel / 2")
add("chi-square2", {}, x, y, 1 + 2 * A, "1 + 2 sklearn additive_chi2_kernel")

# --- formulas in numpy / mpmath ---------------------------------------------------------------------
x, y = pairs(3)
for c in (0.5, 1.0, 2.0):
    add("multiquadratic", dict(c=c), x, y, np.sqrt(euclid(x, y) ** 2 + c * c), "numpy sqrt(d^2 + c^2)")
    add("inverse-multiquadratic", dict(c=c), x, y, 1 / np.sqrt(euclid(x, y) ** 2 + c * c), "numpy 1/sqrt(d^2 + c^2)")
for p in (1.0, 2.0, 1.5):
    add("power", dict(p=p), x, y, -(euclid(x, y) ** p), "numpy -d^p")
    add("log", dict(p=p), x, y, -np.log1p(euclid(x, y) ** p), "numpy -log(1 + d^p)")
    add("generalized-t-student", dict(p=p), x, y, 1 / (1 + euclid(x, y) ** p), "numpy 1/(1 + d^p)")

# pearson VII: Ustun 2006 / Weka Puk form
for sigma, omega in ((1.0, 2.0), (1.0, 1.0), (2.0, 3.0), (0.5, 0.7), (1.5, 8.0)):
    d = euclid(x, y)
    add("pearson", dict(sigma=sigma, omega=omega), x, y,
        1.0 / (1.0 + (2.0 * d * np.sqrt(2.0 ** (1.0 / omega) - 1.0) / sigma) ** 2) ** omega,
        "Ustun 2006 / Weka Puk")

xn, yn = pairs(3, nonneg=True, zeros=True)
add("histogram", {}, xn, yn, np.sum(np.minimum(xn, yn), axis=1), "numpy sum(min)")
add("generalized-histogram", dict(p=1.5), xn, yn, np.sum(np.minimum(xn**1.5, yn**1.5), axis=1), "numpy sum(min(x^p, y^p))")
add("hellinger", {}, xn, yn, np.sum(np.sqrt(xn * yn), axis=1), "numpy sum(sqrt(x y))")

for sigma in (0.5, 1.0, 2.0):
    out = []
    for i in range(len(x)):
        d = mpf_dist(x[i], y[i])
        out.append(mp.mpf(1) if d == 0 else sigma / d * mp.sin(d / sigma))
    add("wave", dict(sigma=sigma), x, y, out, "mpmath (sigma/d) sin(d/sigma)")
for a in (0.5, 1.0, 3.0):
    add("hyperbolic-secant", dict(a=a), x, y, [mp.sech(a * mpf_dist(x[i], y[i])) for i in range(len(x))],
        "mpmath sech(a d)")
for sigma, n, v in ((1.0, 2.0, -1.0), (1.3, 2.0, -0.5), (0.8, 1.0, 1.0), (1.0, 1.5, 0.5)):
    out = []
    for i in range(len(x)):
        d = mpf_dist(x[i], y[i])
        vp = mp.mpf(v) + 1
        out.append(mp.besselj(vp, sigma * d) * d ** (n * vp) if d != 0 else (mp.mpf(1) if vp == 0 else mp.mpf(0)))
    add("bessel", dict(sigma=sigma, n=n, v=v), x, y, out, "mpmath besselj(v+1, sigma d) d^(n (v+1))")
for sigma, order, degree in ((1.0, 0.0, 1.0), (0.7, 1.0, 1.0), (1.2, 1.0, 2.0), (2.0, 0.5, 1.0), (1.0, 2.0, 3.0)):
    out = []
    lim = 1 / (mp.gamma(order + 1) * mp.mpf(2) ** order)
    for i in range(len(x)):
        d = mpf_dist(x[i], y[i])
        b = sigma * d
        if b < mp.mpf("1e-5"):
            out.append(mp.mpf(1))
        else:
            out.append((mp.besselj(order, b) * b ** (-order) / lim) ** degree)
    add("bessel2", dict(sigma=sigma, order=order, degree=degree), x, y, out,
        "mpmath (besselj(order, b) b^-order / lim)^degree")
# dirichlet: product over the components of sin((n + 1/2) delta) / (2 sin(delta / 2)), n + 1/2 at delta = 0
for n in (1.0, 2.0, 4.0):
    xd, yd = pairs(2)
    out = []
    for i in range(len(xd)):
        dirich = mp.mpf(1)
        for a, b in zip(xd[i], yd[i]):
            delta = mp.mpf(float(a)) - mp.mpf(float(b))
            dirich *= (n + mp.mpf(1) / 2) if delta == 0 else mp.sin((n + mp.mpf(1) / 2) * delta) / (2 * mp.sin(delta / 2))
        out.append(dirich)
    add("dirichlet", dict(n=n), xd, yd, out, "mpmath sin((n+1/2) delta) / (2 sin(delta/2)) per component")

x2, y2 = pairs(2)
for sigma, k, d_ in ((1.0, 1.0, 1.0), (0.5, 1.0, 2.0), (1.5, 2.0, 1.0), (0.8, 3.0, 0.5)):
    xx = np.abs(x2) if (k != int(k)) else x2
    yy = np.abs(y2) if (k != int(k)) else y2
    add("anova", dict(sigma=sigma, k=k, d=d_), xx, yy,
        np.sum(np.exp(-sigma * (xx**k - yy**k) ** 2) ** d_, axis=1), "numpy sum_i exp(-sigma (x_i^k - y_i^k)^2)^d")

# spline (kernlab splinedot formula, written here in the product form of Genton / crsouza)
xs, ys = pairs(3, nonneg=True, zeros=True)
mn = np.minimum(xs, ys)
add("spline", {}, xs, ys,
    np.prod(1 + xs * ys * (1 + mn) - (xs + ys) / 2 * mn**2 + mn**3 / 3, axis=1), "numpy (kernlab splinedot formula)")

# b-spline: cardinal B-spline of degree 2n+1 (scipy), zero outside the support
xb = np.array([[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0],
               [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
yb = np.array([[0.0, 0.0], [0.3, -0.2], [1.0, 0.5], [1.4, 1.4], [-2.2, 0.1], [2.9, 0.0], [2.0, 2.0],
               [3.0, 0.0], [4.5, 0.2], [30.0, 0.0], [100.0, 0.0]])
for n in (0, 1, 2, 3, 4, 5):
    deg = 2 * n + 1
    basis = BSpline.basis_element(np.arange(deg + 2) - (deg + 1) / 2.0, extrapolate=False)
    out = []
    for i in range(len(xb)):
        val = 1.0
        for a, b in zip(xb[i], yb[i]):
            t = basis(abs(a - b))
            val *= 0.0 if np.isnan(t) else float(t)
        out.append(val)
    add("b-spline", dict(n=n), xb, yb, out, "scipy BSpline.basis_element of degree 2n+1")


# geometric: intersection volume of two n-balls
def ball_volume(k):
    return mp.pi ** (mp.mpf(k) / 2) / mp.gamma(mp.mpf(k) / 2 + 1)


def lens_ratio(n, s):
    """Intersection volume of two unit n-balls at centre distance 2 s over the ball volume."""
    if s >= 1:
        return mp.mpf(0)
    if s == 0:
        return mp.mpf(1)
    integral = mp.quad(lambda t: (1 - t * t) ** (mp.mpf(n - 1) / 2), [s, (s + 1) / 2, 1])
    ratio = 2 * ball_volume(n - 1) / ball_volume(n) * integral
    beta = mp.betainc((n + 1) / mp.mpf(2), mp.mpf(1) / 2, 0, 1 - s * s, regularized=True)
    assert abs(ratio - beta) < mp.mpf("1e-40"), (n, s, ratio, beta)
    return ratio


svals = [0.0, 1e-9, 0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 0.95, 0.99, 0.999, 0.999999, 1.0, 1.3]
for n in (1, 2, 3, 4, 5, 6, 7, 9, 12, 20):
    for r in (1.0, 0.7):
        xg = np.zeros((len(svals), 1))
        yg = np.array([[2 * r * s] for s in svals])
        out = [lens_ratio(n, mp.mpf(float(yg[i, 0])) / (2 * mp.mpf(r))) for i in range(len(svals))]
        add("geometric", dict(n=n, r=r), xg, yg, out, "mpmath quadrature of the n-ball intersection volume")


# --- EDN output -------------------------------------------------------------------------------------
def edn(v):
    if isinstance(v, dict):
        return "{" + " ".join(":%s %s" % (k, edn(x)) for k, x in v.items()) + "}"
    if isinstance(v, (list, tuple)):
        return "[" + " ".join(edn(x) for x in v) + "]"
    if isinstance(v, str):
        return ":" + v if v in ("abs",) else '"%s"' % v
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    f = float(v)
    s = repr(f)
    assert "inf" not in s and "nan" not in s, s
    return s


def main():
    lines = ["{:cases ["]
    for c in CASES:
        lines.append(
            " {:kernel :%s :params %s :oracle \"%s\" :x %s :y %s :expected %s}"
            % (c["kernel"], edn(c["params"]), c["oracle"], edn(c["x"]), edn(c["y"]), edn(c["expected"]))
        )
    lines.append("]}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n")
    print("sklearn", sklearn.__version__, "scipy", scipy.__version__, "mpmath", mp.__version__, "numpy", np.__version__)
    print("wrote", len(CASES), "cases to", OUT)


if __name__ == "__main__":
    main()
