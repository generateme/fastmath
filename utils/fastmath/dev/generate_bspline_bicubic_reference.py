"""Generator of reference values for `fastmath.interpolation.bspline/b-spline`,
`fastmath.interpolation.bicubic/bicubic` and `fastmath.interpolation.ssj/cubic-smoothing`.

Tools: scipy (`make_interp_spline`, `make_smoothing_spline`), numpy (`polyfit` for the rho = 0 limit
of the smoothing spline) and an independent bicubic Hermite patch clone written here in pure Python.
Versions are printed on the run and stored in the `:versions` key of the output.

Domain:
  :bspline    9 datasets (uniform, nonuniform, noisy, minimal 6 points, plateau, step, huge ys, tiny xs,
              shifted xs) for degrees 1..5 with default knots, query points inside, on every knot,
              and 10 percent of the span outside on both sides (scipy extrapolate=True continues the
              end polynomial); two cases with user knots (clamped with interior knots chosen by hand,
              uniform knots with spacing 1.4 spanning exactly the data).
  :smoothing  `make_smoothing_spline(x, y, w, lam)` with lam = (1 - rho) / rho for rho in
              1, 0.9, 0.5, 0.1, 0.001 (weights none and seeded random; measured equivalence
              of the SSJ weights of `cubic-smoothing` to scipy: w = sqrt(weights)),
              rho = 0 as the least squares
              line (numpy polyfit); end derivatives of scipy to check the linear extrapolation.
  :bicubic    grids uni9x9, uni5x5, nonuni7x6, 2x2, 3x4 and the function
              f(x, y) = sin(1.3 x) cos(0.9 y) + 0.2 x y; patch values at seeded random points inside
              and outside the grid from a Hermite patch built in pure Python with the derivative rule
              of `bicubic`: central differences (f[i+1] - f[i-1]) / (x[i+1] - x[i-1]) at interior
              nodes (the rule of Apache Commons Math BicubicInterpolator), second order one-sided
              three point differences on nonuniform spacing at boundary nodes, the secant for a
              two point axis; the mixed derivative is the tensor product of both stencils;
              the patch of the nearest edge cell is continued outside of the grid.
Why: the SSJ based b-spline was a control polygon (error up to 2.7 against make_interp_spline);
the Apache bicubic took zero derivatives at the boundary (outer ring error 20 times larger).

Produces `test/resources/interpolation/bspline_bicubic_reference.edn`. Deterministic: fixed seeds and
fixed order, doubles are printed with `repr` (round trip). Rerun and diff to check byte identity.

Run from anywhere:
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_bspline_bicubic_reference.py
"""
import math
from pathlib import Path

import numpy as np
import scipy
from scipy.interpolate import make_interp_spline, make_smoothing_spline

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "test" / "resources" / "interpolation" / "bspline_bicubic_reference.edn"


def d(x):
    s = repr(float(x))
    assert "inf" not in s and "nan" not in s, s
    return s


def dv(xs):
    return "[" + " ".join(d(x) for x in xs) + "]"


def dvv(rows):
    return "[" + " ".join(dv(r) for r in rows) + "]"


# ---------------------------------------------------------------- b-spline


def bspline_datasets():
    rng = np.random.default_rng(20261014)
    ds = {}
    x = np.linspace(0.0, 9.0, 10)
    ds["uni10"] = (x, np.sin(x) + 0.1 * x)
    x = np.cumsum(rng.uniform(0.2, 1.8, 12))
    ds["nonuni12"] = (x, np.cos(0.7 * x) * np.exp(-0.05 * x))
    x = np.linspace(-3.0, 3.0, 15)
    ds["noisy15"] = (x, np.exp(-x * x) + 0.1 * rng.standard_normal(15))
    x = np.array([0.0, 0.4, 1.1, 1.5, 2.6, 3.0])
    ds["min6"] = (x, np.array([1.0, -0.5, 2.0, 0.3, -1.2, 0.8]))
    x = np.arange(8.0)
    ds["plateau8"] = (x, np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0]))
    x = np.array([0.0, 1.0, 2.0, 2.01, 3.0, 5.0, 8.0, 8.5, 9.0])
    ds["jump9"] = (x, np.array([0.0, 1.0, 2.0, 2.01, 10.0, 10.5, 3.0, -4.0, 0.0]))
    x = np.linspace(0.0, 4.0, 11)
    ds["huge_y"] = (x, 1.0e8 * np.sin(1.5 * x) + 3.0e7)
    x = np.linspace(0.0, 1.0e-6, 9) * 1.0 + 2.0e-7
    ds["tiny_x"] = (x, np.sin(4.0e6 * x) + 0.5)
    x = np.linspace(1.0e3, 1.0e3 + 12.0, 13) + np.array([0.0, 0.1, -0.1, 0.2, 0.0, 0.3, -0.2, 0.1, 0.0, 0.2, -0.1, 0.0, 0.0])
    ds["shifted13"] = (x, np.log1p(x - 999.0))
    return ds


def queries(x, rng):
    span = x[-1] - x[0]
    inside = list(np.sort(rng.uniform(x[0], x[-1], 24)))
    mids = list(0.5 * (x[1:] + x[:-1]))
    outside = [x[0] - 0.1 * span, x[0] - 0.01 * span, x[-1] + 0.01 * span, x[-1] + 0.1 * span]
    return list(x) + mids + inside + outside


def gen_bspline():
    rng = np.random.default_rng(7)
    ds = bspline_datasets()
    out_ds = []
    cases = []
    for name, (x, y) in ds.items():
        q = queries(x, rng)
        out_ds.append((name, x, y, q))
        for k in (1, 2, 3, 4, 5):
            spl = make_interp_spline(x, y, k=k)
            cases.append((name, k, [float(v) for v in spl(np.array(q))], [float(v) for v in spl.t]))
    # user supplied knots
    custom = []
    x = np.linspace(0.0, 7.0, 8)
    y = np.array([0.5, 1.5, -1.0, 2.0, 0.0, 1.0, -2.0, 0.7])
    t = np.r_[(0.0,) * 4, [1.5, 3.0, 4.0, 5.5], (7.0,) * 4]
    spl = make_interp_spline(x, y, k=3, t=t)
    q = list(np.linspace(-0.7, 7.7, 31))
    custom.append(("clamped-interior", 3, x, y, t, q, [float(v) for v in spl(np.array(q))]))
    t = np.arange(-3, 9) * 1.4
    spl = make_interp_spline(x, y, k=3, t=t)
    q = list(np.linspace(-0.7, 7.7, 31))
    custom.append(("uniform-wide", 3, x, y, t, q, [float(v) for v in spl(np.array(q))]))
    return out_ds, cases, custom


# ---------------------------------------------------------------- smoothing


def gen_smoothing():
    rng = np.random.default_rng(99)
    out = []
    data = {}
    x = np.linspace(0.0, 6.0, 13)
    data["uni13"] = (x, np.sin(x) + 0.2 * rng.standard_normal(13))
    x = np.sort(rng.uniform(0.0, 5.0, 17))
    data["nonuni17"] = (x, np.cos(2.0 * x) + 0.15 * rng.standard_normal(17))
    x = np.array([0.0, 1.0, 2.5, 3.0, 4.2])
    data["five"] = (x, np.array([1.0, 3.0, 2.0, -1.0, 0.5]))
    wts = {}
    for name, (x, y) in data.items():
        wts[name] = rng.uniform(0.5, 2.0, len(x))
    for name, (x, y) in data.items():
        span = x[-1] - x[0]
        q = list(np.linspace(x[0] - 0.3 * span, x[-1] + 0.3 * span, 25)) + list(x)
        for wname, w in (("none", None), ("random", wts[name])):
            for rho in (1.0, 0.9, 0.5, 0.1, 0.001):
                lam = (1.0 - rho) / rho
                if rho == 1.0:
                    lam = 0.0
                # measured equivalence: scipy w = sqrt(weights of cubic-smoothing)
                spl = make_smoothing_spline(x, y, w=None if w is None else np.sqrt(w), lam=lam)
                dspl = spl.derivative()
                out.append(dict(name=name, weights=w, rho=rho, lam=lam, q=q,
                                v=[float(v) for v in spl(np.array(q))],
                                d0=float(dspl(x[0])), d1=float(dspl(x[-1])),
                                v0=float(spl(x[0])), v1=float(spl(x[-1]))))
        # rho = 0: least squares line
        a, b = np.polyfit(x, y, 1)
        out.append(dict(name=name, weights=None, rho=0.0, lam=None, q=q,
                        v=[float(a * v + b) for v in q],
                        d0=float(a), d1=float(a), v0=float(a * x[0] + b), v1=float(a * x[-1] + b)))
    return data, out


# ---------------------------------------------------------------- bicubic clone


def stencil(xs, i):
    """Three point derivative stencil [(index, weight)] of the axis xs at node i."""
    n = len(xs)
    if n == 2:
        h = xs[1] - xs[0]
        return [(0, -1.0 / h), (1, 1.0 / h)]
    if i == 0:
        h1 = xs[1] - xs[0]
        h2 = xs[2] - xs[1]
        return [(0, -(2.0 * h1 + h2) / (h1 * (h1 + h2))),
                (1, (h1 + h2) / (h1 * h2)),
                (2, -h1 / (h2 * (h1 + h2)))]
    if i == n - 1:
        h1 = xs[n - 1] - xs[n - 2]
        h2 = xs[n - 2] - xs[n - 3]
        return [(n - 3, h1 / (h2 * (h1 + h2))),
                (n - 2, -(h1 + h2) / (h1 * h2)),
                (n - 1, (2.0 * h1 + h2) / (h1 * (h1 + h2)))]
    dx = xs[i + 1] - xs[i - 1]
    return [(i - 1, -1.0 / dx), (i + 1, 1.0 / dx)]


def derivatives(xs, ys, f):
    nx, ny = len(xs), len(ys)
    fx = [[sum(w * f[a][j] for a, w in stencil(xs, i)) for j in range(ny)] for i in range(nx)]
    fy = [[sum(w * f[i][b] for b, w in stencil(ys, j)) for j in range(ny)] for i in range(nx)]
    fxy = [[sum(wa * wb * f[a][b] for a, wa in stencil(xs, i) for b, wb in stencil(ys, j))
            for j in range(ny)] for i in range(nx)]
    return fx, fy, fxy


def hermite_basis(u):
    """Return (h00, h01, h10, h11) at u: value at 0, value at 1, derivative at 0, derivative at 1."""
    u2 = u * u
    u3 = u2 * u
    return (2 * u3 - 3 * u2 + 1, -2 * u3 + 3 * u2, u3 - 2 * u2 + u, u3 - u2)


def cell(axis, v):
    n = len(axis)
    i = 0
    for j in range(n - 1):
        if axis[j] <= v:
            i = j
    return min(i, n - 2)


def patch(xs, ys, f, fx, fy, fxy, x, y):
    i = cell(xs, x)
    j = cell(ys, y)
    hx = xs[i + 1] - xs[i]
    hy = ys[j + 1] - ys[j]
    u = (x - xs[i]) / hx
    v = (y - ys[j]) / hy
    bu = hermite_basis(u)
    bv = hermite_basis(v)
    # value basis and (scaled) derivative basis at the corner a (0 or 1)
    vu = (bu[0], bu[1])
    du = (bu[2] * hx, bu[3] * hx)
    vv = (bv[0], bv[1])
    dvv_ = (bv[2] * hy, bv[3] * hy)
    s = 0.0
    for a in (0, 1):
        for b in (0, 1):
            ii, jj = i + a, j + b
            s += f[ii][jj] * vu[a] * vv[b]
            s += fx[ii][jj] * du[a] * vv[b]
            s += fy[ii][jj] * vu[a] * dvv_[b]
            s += fxy[ii][jj] * du[a] * dvv_[b]
    return s


def target(x, y):
    return math.sin(1.3 * x) * math.cos(0.9 * y) + 0.2 * x * y


def gen_bicubic():
    rng = np.random.default_rng(31415)
    grids = {
        "uni9x9": (list(np.linspace(0.0, 4.0, 9)), list(np.linspace(0.0, 4.0, 9))),
        "uni5x5": (list(np.linspace(0.0, 4.0, 5)), list(np.linspace(0.0, 4.0, 5))),
        "nonuni7x6": ([0.0, 0.35, 0.9, 1.7, 2.4, 3.3, 4.0], [0.0, 0.6, 1.3, 2.5, 3.1, 4.0]),
        "grid2x2": ([0.0, 1.5], [0.2, 2.0]),
        "grid3x4": ([0.0, 0.7, 2.0], [0.0, 0.5, 1.5, 2.2]),
    }
    out = []
    for name, (xs, ys) in grids.items():
        f = [[target(x, y) for y in ys] for x in xs]
        fx, fy, fxy = derivatives(xs, ys, f)
        lo_x, hi_x, lo_y, hi_y = xs[0], xs[-1], ys[0], ys[-1]
        pts = []
        for _ in range(40):
            pts.append((float(rng.uniform(lo_x, hi_x)), float(rng.uniform(lo_y, hi_y))))
        # outside of the grid, one point per side and the corners
        spx, spy = hi_x - lo_x, hi_y - lo_y
        for px, py in ((lo_x - 0.3 * spx, 0.5 * (lo_y + hi_y)), (hi_x + 0.3 * spx, 0.4 * (lo_y + hi_y)),
                       (0.45 * (lo_x + hi_x), lo_y - 0.3 * spy), (0.55 * (lo_x + hi_x), hi_y + 0.3 * spy),
                       (lo_x - 0.2 * spx, lo_y - 0.2 * spy), (hi_x + 0.2 * spx, hi_y + 0.2 * spy)):
            pts.append((px, py))
        # exactly on the nodes
        for i in range(len(xs)):
            for j in range(len(ys)):
                pts.append((xs[i], ys[j]))
        vals = [patch(xs, ys, f, fx, fy, fxy, x, y) for x, y in pts]
        out.append((name, xs, ys, f, pts, vals))
    return out


# ---------------------------------------------------------------- output


def main():
    ds, cases, custom = gen_bspline()
    sm_data, sm = gen_smoothing()
    bc = gen_bicubic()
    L = []
    L.append("{")
    L.append(" :versions {:scipy \"%s\" :numpy \"%s\"}" % (scipy.__version__, np.__version__))
    L.append(" :bspline")
    L.append(" {:datasets")
    L.append("  {")
    for name, x, y, q in ds:
        L.append("   \"%s\" {:xs %s :ys %s :queries %s}" % (name, dv(x), dv(y), dv(q)))
    L.append("  }")
    L.append("  :cases")
    L.append("  [")
    for name, k, vals, t in cases:
        L.append("   {:dataset \"%s\" :degree %d :expected %s :knots %s}" % (name, k, dv(vals), dv(t)))
    L.append("  ]")
    L.append("  :custom-knots")
    L.append("  [")
    for name, k, x, y, t, q, vals in custom:
        L.append("   {:name \"%s\" :degree %d :xs %s :ys %s :knots %s :queries %s :expected %s}"
                 % (name, k, dv(x), dv(y), dv(t), dv(q), dv(vals)))
    L.append("  ]}")
    L.append(" :smoothing")
    L.append(" {:datasets")
    L.append("  {")
    for name, (x, y) in sm_data.items():
        L.append("   \"%s\" {:xs %s :ys %s}" % (name, dv(x), dv(y)))
    L.append("  }")
    L.append("  :cases")
    L.append("  [")
    for c in sm:
        L.append("   {:dataset \"%s\" :weights %s :rho %s :lam %s :queries %s :expected %s :d0 %s :d1 %s :v0 %s :v1 %s}"
                 % (c["name"], dv(c["weights"]) if c["weights"] is not None else "nil", d(c["rho"]),
                    d(c["lam"]) if c["lam"] is not None else "nil", dv(c["q"]), dv(c["v"]),
                    d(c["d0"]), d(c["d1"]), d(c["v0"]), d(c["v1"])))
    L.append("  ]}")
    L.append(" :bicubic")
    L.append(" [")
    for name, xs, ys, f, pts, vals in bc:
        L.append("  {:name \"%s\" :xs %s :ys %s :vss %s" % (name, dv(xs), dv(ys), dvv(f)))
        L.append("   :points %s" % dvv([list(p) for p in pts]))
        L.append("   :expected %s}" % dv(vals))
    L.append(" ]")
    L.append("}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(L) + "\n")
    print("scipy", scipy.__version__, "numpy", np.__version__, "wrote", OUT)


if __name__ == "__main__":
    main()
