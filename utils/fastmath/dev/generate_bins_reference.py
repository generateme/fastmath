"""Reference data for the bin estimators in test/fastmath/stats/bins_test.clj

Tools: Python 3.14, numpy 2.5.2.
Run: cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_bins_reference.py
Output: test/resources/stats/bins_reference.edn (regeneration is byte-identical: legacy
RandomState streams, fixed order, repr() of doubles).

Why: fastmath.stats.bins implements the 2d normal-reference rules (scott-2d, scott-correlated,
freedman-diaconis-2d), leave-one-out cross-validation (Rudemo 1982, numpy `stone`) and the
Shimazaki-Shinomoto cost in 1d and 2d. This script is an independent numpy implementation of
the same specification:
- percentiles: numpy method 'weibull' (R type 6), which equals Commons Math
  StatUtils/percentile used by fastmath (checked on a sample);
- 1d candidates: log-spaced bin counts between max(1, floor(fd/4)) and min(20 fd, n,
  max(1, floor(range / (2 * smallest gap)))), rounded half up, duplicates removed; fd is the
  Freedman-Diaconis count ceil(range / (2 IQR / cbrt(n)));
- 1d bins: index floor((x - origin) / h) clamped to the bins; shifted origins min - j h / s with
  one extra bin for j > 0;
- 1d costs: cross-validation (2 - (n+1) sum p^2) / h (numpy `_hist_bin_stone` jhat), SS
  (2 mean - biased var) / h^2 over all bins; ties to the smaller bin count;
- 2d candidates (square cells only, other cell types have no external oracle): sizes
  log-spaced between s/20 and 4 s, s = sqrt(freedman-diaconis-2d area), not below
  sqrt(4 gap_x gap_y); cells floor((x - sv_x - shift) / size), shifts j size / s along the
  diagonal; CV over occupied cells, SS with the cells of the index range of occupied cells
  (counted as doubles); ties to the larger size; the 'outlier' set (one point at
  (1e10, -1e10)) has an index range beyond a long key;
- numpy `stone` bin count of every 1d set, for the sets where its optimum lies inside the
  candidate range.
"""
import math
import os
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                   "test", "resources", "stats", "bins_reference.edn")
C2 = (48 * math.pi) ** 0.25

rs = np.random.RandomState(20261116)
sets1 = {
    "bimodal": np.concatenate([rs.normal(-2, 0.5, 600), rs.normal(1.5, 1.0, 400)]),
    "normal": rs.normal(10, 2, 500),
    "t3": rs.standard_t(3, 500),
    "integers": np.round(rs.normal(0, 3, 1000)),
    "three": np.array([1.0, 2.5, 7.0]),
}
x = rs.normal(0, 1, 500)
sets2 = {
    "correlated": (x, 0.9 * x + math.sqrt(1 - 0.81) * rs.normal(0, 1, 500)),
    "mixture": tuple(np.concatenate(p) for p in zip(
        (rs.normal(-2, 0.5, 200), rs.normal(0, 0.5, 200)),
        (rs.normal(2, 0.7, 200), rs.normal(1, 0.7, 200)),
        (rs.normal(0, 1.0, 200), rs.normal(-2, 0.3, 200)))),
    "anisotropic": (rs.normal(0, 1, 400), rs.normal(0, 4, 400)),
    "integers": (np.round(rs.normal(0, 3, 400)), np.round(rs.normal(0, 2, 400))),
    "outlier": (np.append(rs.normal(0, 1, 500), 1.0e10), np.append(rs.normal(0, 1, 500), -1.0e10)),
}


def pct(a, q):
    return float(np.percentile(a, q, method="weibull"))


def gap(a):
    d = np.diff(np.sort(a))
    d = d[d > 0]
    return float(d.min()) if d.size else math.inf


def log_spaced(lo, hi, k):
    if k < 2 or lo == hi:
        return [lo]
    step = math.log(hi / lo) / (k - 1)
    return [lo * math.exp(step * i) for i in range(k)]


def fd_count(a):
    n = len(a)
    h = 2 * (pct(a, 75) - pct(a, 25)) / n ** (1 / 3)
    return max(1, math.ceil(np.ptp(a) / h))


def candidates_1d(a, k=200):
    n = len(a)
    rng = float(np.ptp(a))
    fd = fd_count(a)
    gl = max(1.0, math.floor(rng / (2 * gap(a))))
    hi = int(min(20.0 * fd, min(float(n), gl)))
    lo = min(hi, max(1, math.floor(fd / 4)))
    out = []
    for v in log_spaced(float(lo), float(hi), k):
        r = math.floor(v + 0.5)
        if r not in out:
            out.append(r)
    return out


def cost_1d(a, method, nbins, shifts):
    n = len(a)
    mn, mx = float(a.min()), float(a.max())
    h = (mx - mn) / nbins
    total = 0.0
    for j in range(shifts):
        nb = nbins if j == 0 else nbins + 1
        origin = mn - j * h / shifts
        idx = np.clip(np.floor((a - origin) / h).astype(np.int64), 0, nb - 1)
        k = np.bincount(idx, minlength=nb).astype(float)
        ss = float((k * k).sum())
        if method == "cv":
            total += (2 - (n + 1) * ss / (n * n)) / h
        else:
            kbar = n / nb
            total += (2 * kbar - (ss / nb - kbar * kbar)) / (h * h)
    return total / shifts


def best_1d(a, method, shifts):
    if len(a) < 2 or np.ptp(a) == 0:
        return 1
    best, bc = 1, math.inf
    for nb in candidates_1d(a):
        c = cost_1d(a, method, nb, shifts)
        if c < bc:
            best, bc = nb, c
    return best


def sd(a):
    return float(np.std(a, ddof=1))


def scott2d(xs, ys):
    f = len(xs) ** -0.25
    return C2 * sd(xs) * f * C2 * sd(ys) * f


def line_area(c, spread, n):
    h = c * spread * n ** (-1 / 3)
    return h * h


def scott_corr(xs, ys):
    # points on a line (one constant axis, or 1 - rho^2 below 4096 ulp(1)): 1d Scott width along it
    n = len(xs)
    vx, vy = np.var(xs, ddof=1), np.var(ys, ddof=1)
    if (vx == 0) != (vy == 0):
        return line_area(3.5, math.sqrt(vx + vy), n)
    cov = np.cov(xs, ys, ddof=1)[0, 1]
    one = 1 - cov * cov / (vx * vy)
    if one < 4096 * 2.0 ** -52:
        return line_area(3.5, math.sqrt(vx + vy), n)
    return float(math.sqrt(96 * math.pi * math.sqrt(vx * vy) * one ** 1.5 / (n * (1 / vx + 1 / vy))))


def fd2d(xs, ys):
    f = len(xs) ** -0.25
    sx = (pct(xs, 75) - pct(xs, 25)) / 1.349
    sy = (pct(ys, 75) - pct(ys, 25)) / 1.349
    return C2 * sx * f * C2 * sy * f


def best_2d_square(xs, ys, method, sv, shifts, k=200):
    n = len(xs)
    base = math.sqrt(fd2d(xs, ys))
    min_size = math.sqrt(4 * gap(xs) * gap(ys))
    sizes = [s for s in log_spaced(base / 20, base * 4, k) if s >= min_size] or [min_size]
    best, bc = None, math.inf
    for s in reversed(sizes):
        total = 0.0
        for j in range(shifts):
            sh = j * s / shifts
            q = np.floor((xs - sv[0] - sh) / s).astype(np.int64)
            r = np.floor((ys - sv[1] - sh) / s).astype(np.int64)
            _, cnt = np.unique(np.stack([q, r], axis=1), axis=0, return_counts=True)
            ss = float((cnt.astype(float) ** 2).sum())
            a = s * s
            if method == "cv":
                total += (2 - (n + 1) * ss / (n * n)) / a
            else:
                m = float(q.max() - q.min() + 1) * float(r.max() - r.min() + 1)
                kbar = n / m
                total += (2 * kbar - (ss / m - kbar * kbar)) / (a * a)
        c = total / shifts
        if c < bc:
            best, bc = s, c
    return best * best


def num(v):
    return "##NaN" if math.isnan(v) else repr(float(v))


def vec(vs):
    return "[" + " ".join(num(v) for v in vs) + "]"


def ivec(vs):
    return "[" + " ".join(str(int(v)) for v in vs) + "]"


lines = [";; generated by utils/fastmath/dev/generate_bins_reference.py, do not edit", "{:one-d {"]
for name, a in sets1.items():
    stone = len(np.histogram_bin_edges(a, "stone")) - 1
    cands = candidates_1d(a) if np.ptp(a) > 0 else [1]
    lines.append(f" :{name} {{:data {vec(a)}")
    lines.append(f"  :candidates {ivec(cands)}")
    lines.append(f"  :cv-costs {vec([cost_1d(a, 'cv', nb, 1) for nb in cands])}")
    lines.append(f"  :ss-costs {vec([cost_1d(a, 'ss', nb, 1) for nb in cands])}")
    lines.append(f"  :cv {best_1d(a, 'cv', 1)} :ss {best_1d(a, 'ss', 1)}")
    lines.append(f"  :cv-shifts-5 {best_1d(a, 'cv', 5)} :ss-shifts-5 {best_1d(a, 'ss', 5)}")
    lines.append(f"  :stone {stone} :stone-in-range {'true' if cands[0] <= stone <= cands[-1] else 'false'}}}")
lines.append("}")
lines.append(" :two-d {")
for name, (xs, ys) in sets2.items():
    lines.append(f" :{name} {{:xs {vec(xs)}")
    lines.append(f"  :ys {vec(ys)}")
    lines.append(f"  :scott {num(scott2d(xs, ys))} :scott-correlated {num(scott_corr(xs, ys))}"
                 f" :freedman-diaconis {num(fd2d(xs, ys))}")
    for method in ("cv", "ss"):
        for tag, sv, shifts in (("", (0.0, 0.0), 1), ("-sv", (0.3, 0.7), 1), ("-shifts-3", (0.0, 0.0), 3)):
            lines.append(f"  :{method}{tag} {num(best_2d_square(xs, ys, method, sv, shifts))}")
    lines[-1] += "}"
lines.append("}}")

with open(OUT, "w") as f:
    f.write("\n".join(lines) + "\n")
