"""Generator of reference values for `fastmath.interpolation.akima/akima` and `fastmath.interpolation.loess/loess`.

Tools: scipy (`Akima1DInterpolator` with method 'akima' and 'makima', `CubicSpline` natural) and
statsmodels (`nonparametric.lowess` with delta = 0 and, for user weights, the Cython routine
`_smoothers_lowess.lowess` with `resid_weights`), numpy for the seeded data. The version numbers are
written into the header of the output file.

Domain:
  * akima: 12 datasets (smooth uniform and nonuniform, plateau, extrema, noisy 15 points, 6 and 5 points,
    exactly linear data, a unit step, tiny and huge scales, a dataset with a repeated slope pattern) and
    3 tiny datasets with 2, 3 and 4 points (scipy accepts them; fastmath requires 5 and the algorithm is
    tested through a private function). Queries: 40 points on a grid reaching 25 percent of the range
    beyond both ends, all knots and all interval midpoints. The interpolants are evaluated with
    extrapolate = True (the end cubics continue).
  * loess: 4 datasets (uniform 40 points, nonuniform 50 points with 2 outliers, noisy 15 points, 7 points)
    with bandwidth in (0.15, 0.3, 0.5, 0.8, 1.0), robustness iterations 0 and 2, smoothed values at the
    knots with delta = 0.0; plus the natural cubic spline (scipy CubicSpline, bc_type natural) through the
    smoothed values at queries (10 percent beyond both ends, knots excluded); plus user weights (iters 0).
Why: the former implementation (Apache Commons Math) differed from scipy up to 0.25 on the end intervals of
Akima and from statsmodels up to 1.45 in LOESS (window not the true k nearest neighbours, upper median).

Produces `test/resources/interpolation/akima_loess_reference.edn`.

Run from anywhere:
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_akima_loess_reference.py
"""
from pathlib import Path

import numpy as np
import scipy
import statsmodels
from scipy.interpolate import Akima1DInterpolator, CubicSpline
from statsmodels.nonparametric._smoothers_lowess import lowess as lowess_core
from statsmodels.nonparametric.smoothers_lowess import lowess

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "test" / "resources" / "interpolation" / "akima_loess_reference.edn"


def edn_double(x):
    s = repr(float(x))
    assert "inf" not in s and "nan" not in s, s
    return s


def edn_vec(a):
    return "[" + " ".join(edn_double(v) for v in a) + "]"


def sorted_uniform(rng, n, lo, hi):
    return np.sort(rng.uniform(lo, hi, n))


def akima_datasets():
    rng = np.random.default_rng(20261014)
    ds = []
    x = np.linspace(0.0, 2.0 * np.pi, 12)
    ds.append(("smooth-uniform", x, np.sin(x)))
    x = sorted_uniform(rng, 14, -3.0, 5.0)
    ds.append(("smooth-nonuniform", x, np.exp(-x / 4.0) * np.cos(1.7 * x)))
    x = np.arange(12.0)
    ds.append(("plateau", x, np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 3, 3, 3], dtype=float)))
    x = np.array([0.0, 0.7, 1.1, 2.0, 2.4, 3.9, 4.0, 5.5, 6.1, 7.3])
    ds.append(("extrema", x, np.array([0.0, 3.0, -2.0, 4.0, -1.0, 5.0, -3.0, 2.0, 0.0, 1.0])))
    x = sorted_uniform(rng, 15, 0.0, 10.0)
    ds.append(("noisy15", x, np.sin(x) + 0.3 * rng.standard_normal(15)))
    x = np.array([0.0, 1.0, 2.5, 3.0, 4.5, 6.0])
    ds.append(("small6", x, np.array([1.0, 2.0, 0.5, 3.0, 2.0, 2.5])))
    x = np.array([0.0, 1.0, 2.0, 3.5, 4.0])
    ds.append(("minimal5", x, np.array([2.0, -1.0, 0.5, 3.0, 1.0])))
    x = np.array([-2.0, -1.5, -0.3, 0.2, 1.1, 2.0, 3.7])
    ds.append(("linear", x, 2.0 * x + 1.0))
    x = np.arange(10.0)
    ds.append(("step", x, np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1], dtype=float)))
    x = np.linspace(0.0, 1.0, 9) * 1.0e-8
    ds.append(("tiny-scale", x, 1.0e-9 * np.array([1.0, 3.0, 2.0, 5.0, 4.0, 4.5, 7.0, 6.0, 8.0])))
    x = np.linspace(0.0, 1.0, 9) * 1.0e8
    ds.append(("huge-scale", x, 1.0e9 * np.array([1.0, 3.0, 2.0, 5.0, 4.0, 4.5, 7.0, 6.0, 8.0])))
    x = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    ds.append(("repeated-slopes", x, np.array([0.0, 1.0, 2.0, 3.0, 3.0, 3.0, 4.0, 5.0])))
    return ds


def akima_small():
    return [
        ("two", np.array([0.0, 1.5]), np.array([1.0, 4.0])),
        ("three", np.array([0.0, 1.0, 3.0]), np.array([1.0, 0.0, 2.0])),
        ("four", np.array([0.0, 1.0, 2.0, 4.0]), np.array([1.0, 3.0, 2.0, 5.0])),
    ]


def queries(x, extend):
    lo, hi = x[0], x[-1]
    r = hi - lo
    grid = np.linspace(lo - extend * r, hi + extend * r, 40)
    mids = 0.5 * (x[1:] + x[:-1])
    return np.concatenate([grid, x, mids])


def akima_entry(name, x, y):
    q = queries(x, 0.25)
    out = {}
    for method in ("akima", "makima"):
        f = Akima1DInterpolator(x, y, method=method, extrapolate=True)
        out[method] = f(q)
    return (name, x, y, q, out["akima"], out["makima"])


def loess_datasets():
    rng = np.random.default_rng(8675309)
    ds = []
    x = np.linspace(0.0, 10.0, 40)
    ds.append(("uniform40", x, np.sin(x) + 0.25 * rng.standard_normal(40)))
    x = sorted_uniform(rng, 50, 0.0, 10.0)
    y = np.sin(x) + 0.25 * rng.standard_normal(50)
    y[5] += 8.0
    y[30] -= 6.0
    ds.append(("nonuniform50", x, y))
    x = sorted_uniform(rng, 15, 0.0, 10.0)
    ds.append(("noisy15", x, np.sin(x) + 0.3 * rng.standard_normal(15)))
    x = np.array([0.0, 0.5, 1.5, 2.0, 3.5, 4.0, 6.0])
    ds.append(("seven", x, np.array([1.0, 2.0, 0.0, 3.0, 2.5, 4.0, 1.0])))
    return ds


FRACS = (0.15, 0.3, 0.5, 0.8, 1.0)
ITERS = (0, 2)


def loess_entries(ds):
    rows = []
    for name, x, y in ds:
        n = len(x)
        for frac in FRACS:
            k = int(frac * n + 1e-10)
            if k < 2:
                continue
            for it in ITERS:
                smooth = lowess(y, x, frac=frac, it=it, delta=0.0, return_sorted=True)[:, 1]
                cs = CubicSpline(x, smooth, bc_type="natural", extrapolate=True)
                r = x[-1] - x[0]
                q = np.linspace(x[0] - 0.1 * r, x[-1] + 0.1 * r, 37)
                rows.append((name, frac, it, k, smooth, q, cs(q)))
    return rows


def weighted_entries(ds):
    rng = np.random.default_rng(424242)
    rows = []
    for name, x, y in ds:
        if name not in ("noisy15", "nonuniform50"):
            continue
        n = len(x)
        w = rng.uniform(0.2, 3.0, n)
        w[3] = 0.0
        for frac in (0.3, 0.5):
            res, _ = lowess_core(y, x, x, w, frac, 0, 0.0, 0)
            rows.append((name, frac, w, res[:, 1]))
    return rows


def main():
    lines = [
        ";; Generated by utils/fastmath/dev/generate_akima_loess_reference.py (do not edit)",
        ";; scipy %s, statsmodels %s, numpy %s" % (scipy.__version__, statsmodels.__version__, np.__version__),
        "{",
    ]
    lines.append(" :akima [")
    for name, x, y, q, ak, mk in [akima_entry(*d) for d in akima_datasets()]:
        lines.append(
            "  {:name \"%s\" :xs %s :ys %s :queries %s :akima %s :makima %s}"
            % (name, edn_vec(x), edn_vec(y), edn_vec(q), edn_vec(ak), edn_vec(mk))
        )
    lines.append(" ]")
    lines.append(" :akima-small [")
    for name, x, y, q, ak, mk in [akima_entry(*d) for d in akima_small()]:
        lines.append(
            "  {:name \"%s\" :xs %s :ys %s :queries %s :akima %s :makima %s}"
            % (name, edn_vec(x), edn_vec(y), edn_vec(q), edn_vec(ak), edn_vec(mk))
        )
    lines.append(" ]")
    ds = loess_datasets()
    lines.append(" :loess-data [")
    for name, x, y in ds:
        lines.append("  {:name \"%s\" :xs %s :ys %s}" % (name, edn_vec(x), edn_vec(y)))
    lines.append(" ]")
    lines.append(" :loess [")
    for name, frac, it, k, smooth, q, spl in loess_entries(ds):
        lines.append(
            "  {:name \"%s\" :bandwidth %s :iters %d :k %d :smooth %s :queries %s :spline %s}"
            % (name, edn_double(frac), it, k, edn_vec(smooth), edn_vec(q), edn_vec(spl))
        )
    lines.append(" ]")
    lines.append(" :loess-weighted [")
    for name, frac, w, smooth in weighted_entries(ds):
        lines.append(
            "  {:name \"%s\" :bandwidth %s :weights %s :smooth %s}"
            % (name, edn_double(frac), edn_vec(w), edn_vec(smooth))
        )
    lines.append(" ]")
    lines.append("}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n")
    print("scipy", scipy.__version__, "statsmodels", statsmodels.__version__, "wrote", OUT)


if __name__ == "__main__":
    main()
