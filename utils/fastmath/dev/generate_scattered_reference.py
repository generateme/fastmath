"""Generator of reference values for the scattered data interpolators of `fastmath.interpolation`.

Tools: numpy 2.5 (inverse distance weighting and the normal equations of the regularized RBF system,
written out in float64), scipy 1.18 (`scipy.interpolate.RBFInterpolator`), scikit-learn 1.9
(`GaussianProcessRegressor`, `ConstantKernel * RBF`, fixed hyperparameters, `optimizer=None`).
Doubles are printed with `repr` (shortest round trip representation).

Domain (all data sets are generated with `numpy.random.default_rng` and fixed seeds, the order is fixed):
  * :shepard - 1-d (7 nodes), 2-d (15) and 3-d (20) nodes with 30 queries each (never at a node) and
    p in 0, 1, 2, 3.5; value `sum(y / d^p) / sum(1 / d^p)` (`fastmath.interpolation.shepard/shepard`).
  * :rbf - 2-d (30 nodes) and 3-d (40 nodes) data, kernels gaussian (epsilon = 1 / scale), thin plate spline
    (degree 1), linear (`-r`, degree 0), cubic (degree 1), multiquadric (`-sqrt(1 + (epsilon r)^2)`, degree 0)
    with the smoothing 0 (interpolation, 500 queries and the values at the nodes), 1e-3 and 1e-1 (100 queries).
    `scipy` `smoothing` is the ridge `(K + smoothing I) w = y` with an unregularized polynomial block, this
    is `:regularization :ridge` of `fastmath.interpolation.rbf/rbf`. `:kernel`, `:params` name the matching
    `fastmath.kernel/rbf` kernel (`:radial-powers` with beta 1 has the sign `-r` as scipy `linear`).
  * :rbf-normal - the former behaviour of `rbf`, `(A^T A + lambda I) w = A^T y` with A the (augmented)
    system matrix, solved with numpy for the gaussian kernel (kscale 1 and 10, lambda 1e-4 and 1e-2, without
    polynomial and with a constant + linear polynomial); `:regularization :normal` of `rbf`.
  * :gp - `gaussian-process` (gaussian kernel with sigma, `kscale`, `noise` = alpha) against sklearn:
    mean, standard deviation, full covariance of the latent function, and the log marginal likelihood `L` (of the
    normalised `ys`). With `:normalize? true` the `ys` are normalised here with the sample standard deviation
    (ddof 1, as fastmath does; sklearn `normalize_y` uses ddof 0) and the prediction is mapped back.
    Includes 400 uniform points in [0, 1] (sigma 0.3, noise 1e-8; the determinant of the covariance underflows
    to zero) and 40 points with sigma 0.01 and kscale 1e20.
Why: the audit found shepard NaN at exact and near hits, a `:lambda` different from the usual smoothing, and a gp `L`
which was +Inf/-Inf from `log(det)`.

Produces `test/resources/interpolation/scattered_reference.edn`.

Run from anywhere:
    cd /home/ts/penv && uv run python \\
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_scattered_reference.py
"""
from pathlib import Path

import numpy as np
import scipy
import sklearn
from scipy.interpolate import RBFInterpolator
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "test" / "resources" / "interpolation" / "scattered_reference.edn"


def d(x):
    s = repr(float(x))
    assert "inf" not in s and "nan" not in s, s
    return s


def vec(xs):
    return "[" + " ".join(d(x) for x in xs) + "]"


def mat(rows):
    return "[" + " ".join(vec(r) for r in rows) + "]"


def smooth_fn(x):
    # a smooth test function of any dimension
    x = np.atleast_2d(x)
    out = np.sin(3.0 * x[:, 0]) + 0.5 * np.cos(2.0 * x[:, 1 % x.shape[1]])
    if x.shape[1] > 1:
        out = out + x[:, 0] * x[:, -1]
    return out


# ---------------------------------------------------------------- shepard
def idw(xs, ys, q, p):
    dist = np.sqrt(((xs - q) ** 2).sum(axis=1))
    w = 1.0 / dist ** p
    return float((w * ys).sum() / w.sum())


def shepard_cases():
    lines = []
    for name, dim, n, seed in (("d1", 1, 7, 101), ("d2", 2, 15, 102), ("d3", 3, 20, 103)):
        rng = np.random.default_rng(seed)
        xs = rng.random((n, dim)) * 4.0 - 2.0
        ys = rng.normal(size=n) * 3.0 + 1.0
        qs = rng.random((30, dim)) * 6.0 - 3.0
        for p in (0.0, 1.0, 2.0, 3.5):
            vals = [idw(xs, ys, q, p) for q in qs]
            lines.append("  {:name :%s :p %s :xs %s :ys %s :queries %s :values %s}"
                         % (name, d(p), mat(xs), vec(ys), mat(qs), vec(vals)))
    return lines


# ---------------------------------------------------------------- rbf
RBF_KERNELS = (
    # name, scipy kernel, scale (epsilon = 1 / scale), degree, fastmath kernel, fastmath params
    ("gaussian", "gaussian", 0.4, -1, "gaussian", "{:scale 0.4}"),
    ("thin-plate", "thin_plate_spline", 1.0, 1, "thin-plate", "{}"),
    ("linear", "linear", 1.0, 0, "radial-powers", "{:beta 1.0}"),
    ("cubic", "cubic", 1.0, 1, "radial-powers-3", "{}"),
    ("multiquadric", "multiquadric", 0.5, 0, "generalized-multiquadratic", "{:beta 0.5 :scale 0.5}"),
)
SMOOTHINGS = (0.0, 1e-3, 1e-1)


def rbf_cases():
    lines = []
    for dim, n, seed in ((2, 30, 201), (3, 40, 202)):
        rng = np.random.default_rng(seed)
        xs = rng.random((n, dim))
        ys = smooth_fn(xs)
        qs = rng.random((500, dim))
        for name, sk, scale, degree, fk, fparams in RBF_KERNELS:
            cases = []
            for smoothing in SMOOTHINGS:
                it = RBFInterpolator(xs, ys, kernel=sk, epsilon=1.0 / scale, degree=degree, smoothing=smoothing)
                nq = 500 if smoothing == 0.0 else 100
                vals = it(qs[:nq])
                extra = ""
                if smoothing == 0.0:
                    extra = " :node-values %s" % vec(it(xs))
                cases.append("{:smoothing %s :values %s%s}" % (d(smoothing), vec(vals), extra))
            lines.append("  {:name :%s-%dd :kernel :%s :params %s :scale %s :degree %d :xs %s :ys %s :queries %s :cases [%s]}"
                         % (name, dim, fk, fparams, d(scale), degree, mat(xs), vec(ys), mat(qs), " ".join(cases)))
    return lines


def gauss_kernel_matrix(a, b, scale):
    d2 = ((a[:, None, :] - b[None, :, :]) ** 2).sum(axis=2)
    return np.exp(-d2 / scale ** 2)


def rbf_normal_cases():
    lines = []
    rng = np.random.default_rng(301)
    n, dim, scale = 20, 2, 0.5
    xs = rng.random((n, dim))
    ys = smooth_fn(xs)
    qs = rng.random((25, dim))
    for degree in (-1, 1):
        for kscale in (1.0, 10.0):
            for lam in (1e-4, 1e-2):
                K = kscale * gauss_kernel_matrix(xs, xs, scale)
                Kq = kscale * gauss_kernel_matrix(qs, xs, scale)
                if degree == 1:
                    P = np.hstack([np.ones((n, 1)), xs])
                    Pq = np.hstack([np.ones((len(qs), 1)), qs])
                    A = np.block([[K, P], [P.T, np.zeros((3, 3))]])
                    b = np.concatenate([ys, np.zeros(3)])
                else:
                    A, b = K, ys
                w = np.linalg.solve(A.T @ A + lam * np.eye(A.shape[0]), A.T @ b)
                vals = Kq @ w[:n] + (Pq @ w[n:] if degree == 1 else 0.0)
                lines.append("  {:scale %s :kscale %s :lambda %s :degree %d :xs %s :ys %s :queries %s :values %s}"
                             % (d(scale), d(kscale), d(lam), degree, mat(xs), vec(ys), mat(qs), vec(vals)))
    return lines


# ---------------------------------------------------------------- gp
def gp_case(name, xs, ys, qs, sigma, kscale, noise, normalize):
    ys = np.asarray(ys, dtype=float)
    ymean, ystd = (ys.mean(), ys.std(ddof=1)) if normalize else (0.0, 1.0)
    yn = (ys - ymean) / ystd
    kernel = ConstantKernel(kscale, constant_value_bounds="fixed") * RBF(sigma, length_scale_bounds="fixed")
    gpr = GaussianProcessRegressor(kernel=kernel, alpha=noise, optimizer=None, normalize_y=False).fit(xs, yn)
    mu, cov = gpr.predict(qs, return_cov=True)
    _, sd = gpr.predict(qs, return_std=True)
    L = gpr.log_marginal_likelihood(kernel.theta)
    return ("  {:name :%s :sigma %s :kscale %s :noise %s :normalize? %s :xs %s :ys %s :queries %s :mean %s :std %s :cov %s :L %s}"
            % (name, d(sigma), d(kscale), d(noise), "true" if normalize else "false", mat(xs), vec(ys), mat(qs),
               vec(mu * ystd + ymean), vec(sd * ystd), mat(cov * ystd ** 2), d(L)))


def gp_cases():
    lines = []
    rng = np.random.default_rng(401)
    xs = np.sort(rng.random(12)).reshape(-1, 1) * 4.0
    ys = np.sin(xs[:, 0]) + 0.1 * rng.normal(size=12)
    qs = np.linspace(-0.5, 4.5, 15).reshape(-1, 1)
    lines.append(gp_case("plain-1d", xs, ys, qs, 0.6, 1.5, 1e-3, False))

    rng = np.random.default_rng(402)
    xs = rng.random((25, 2)) * 3.0
    ys = 100.0 + 20.0 * smooth_fn(xs) + rng.normal(size=25)
    qs = rng.random((12, 2)) * 3.0
    lines.append(gp_case("normalized-2d", xs, ys, qs, 0.7, 2.0, 1e-4, True))

    rng = np.random.default_rng(403)
    xs = rng.random((18, 3))
    ys = smooth_fn(xs) * 5.0 + 3.0
    qs = rng.random((10, 3))
    lines.append(gp_case("normalized-3d-noise", xs, ys, qs, 0.6, 1.0, 1e-2, True))

    rng = np.random.default_rng(404)
    xs = np.sort(rng.random(400)).reshape(-1, 1)
    ys = np.sin(6.0 * xs[:, 0])
    qs = np.linspace(0.0, 1.0, 9).reshape(-1, 1)
    lines.append(gp_case("uniform-400", xs, ys, qs, 0.3, 1.0, 1e-8, False))

    rng = np.random.default_rng(405)
    xs = np.sort(rng.random(40)).reshape(-1, 1)
    ys = np.cos(10.0 * xs[:, 0]) + 0.1 * rng.normal(size=40)
    qs = np.linspace(0.0, 1.0, 7).reshape(-1, 1)
    lines.append(gp_case("narrow-40", xs, ys, qs, 0.01, 1e20, 1e-8, False))
    return lines


def main():
    sections = [
        (":shepard", shepard_cases()),
        (":rbf", rbf_cases()),
        (":rbf-normal", rbf_normal_cases()),
        (":gp", gp_cases()),
    ]
    out = ["{"]
    for key, lines in sections:
        out.append(" %s [" % key)
        out.extend(lines)
        out.append(" ]")
    out.append("}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(out) + "\n")
    print("scipy", scipy.__version__, "sklearn", sklearn.__version__, "numpy", np.__version__,
          "wrote", sum(len(l) for _, l in sections), "cases to", OUT)


if __name__ == "__main__":
    main()
