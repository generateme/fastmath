"""Generator of test/resources/kernel/rlcv_reference.edn: LSCV / LCV / RLCV criteria and their optima for the
14 KDE kernels of `fastmath.kernel.density`.

Tools: numpy 2.5.2, scipy 1.18.1; the criteria are the independent implementation of
`generate_bandwidth_reference.py` (imported; read its docstring for the exact definitions of f, the leave one
out f_i, LSCV, LCV, RLCV with the `:printed` and `:continuous` bias term, the extents, the integration grid
and the search scale s = nrd(z) * nrd-scale).  Wu, "Robust Likelihood Cross Validation for Kernel Density
Estimation", JBES 37(4), 2019 (a_n = sqrt(ln n) / (sqrt 2 n) for d = 1 on standardised data).  Double precision.

Domain: 6 datasets (normal n 30 / 100 / 300, lognormal n 100 (skewed), bimodal n 200, exponential n 150), seeded
numpy `default_rng`, x 14 kernels = 84 cases.  For each case, on z = (x - mean) / sd (sample sd), h in
standardised units:
  * the global optimum of each criterion: 4000 point log grid over [0.01, 5] s refined by bounded Brent in
    log h (xatol 1e-9), with the criterion value and `:at-bound` (within 1e-3 relative of a grid end);
  * five fixed bandwidths h = s * [0.25 0.5 1 2 4] with all components (int-f2, mean-loo, mean-log-loo, mean-g,
    b-printed, b-continuous, measure-above-a = int I(f >= a) dt, so that b-printed - b-continuous =
    a/2 * measure-above-a) and the criteria lscv, lcv, rlcv-printed, rlcv-continuous, with the integrals on
    the grid of spacing h/50 (`:fixed`) and h/1000 (`:fine`; the trapezoid error of the kinked kernels is
    O((h/50)^2), so a consumer integrating with another spacing should compare with `:fine` and a tolerance
    that covers `:fixed` - `:fine`).
LCV is undefined (NaN) when some leave-one-out density is negative (silverman) and -Inf when it is 0 (compact
kernels with an isolated point): the stored value is then ##NaN or ##-Inf, and the optimum is searched among the
finite values.  Cauchy: extent 100 (the specified 1000 is too wide for a direct sum), the mass cut is
3.18e-3 per kernel (see the docstring of `generate_bandwidth_reference.py`).

Schema (the file is one map):
    {:cv-info {...}
     :datasets [{:name :n :mean :sd :a :data}]            original data (unstandardised); :a = rlcv a_n
     :cases [{:dataset :kernel :s :h-range [lo hi]
              :lscv {:h :value :at-bound :grid-index} :lcv {...} :rlcv-printed {...} :rlcv-continuous {...}
              :fixed-h [h1 ... h5]
              :fixed {:int-f2 [..5] :mean-loo [] :mean-log-loo [] :mean-g [] :b-printed [] :b-continuous []
                      :measure-above-a [] :lscv [] :lcv [] :rlcv-printed [] :rlcv-continuous []}
              :fine {same keys}}]}
`:h` is in standardised units.  Why: KB-01 ... KB-06, KB-11, KB-12; work packages 2.3 and 2.5 of the Kernel
Density audit.

Run from anywhere (a rerun gives a byte-identical file; ~10 min with 14 processes):
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_rlcv_reference.py
"""
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import scipy

sys.path.insert(0, str(Path(__file__).resolve().parent))
import generate_bandwidth_reference as gb  # noqa: E402
from generate_bandwidth_reference import Kw  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "test" / "resources" / "kernel" / "rlcv_reference.edn"
FACTORS = [0.25, 0.5, 1.0, 2.0, 4.0]
COMP = ["int-f2", "mean-loo", "mean-log-loo", "mean-g", "b-printed", "b-continuous", "measure-above-a",
        "lscv", "lcv", "rlcv-printed", "rlcv-continuous"]


def datasets():
    out = []
    seed = [2026102100]

    def rng():
        seed[0] += 1
        return np.random.default_rng(seed[0])

    for n in (30, 100, 300):
        out.append((f"normal-n{n}", rng().normal(0, 1, n)))
    out.append(("lognormal-n100", rng().lognormal(0, 1, 100)))
    g = rng()
    c = g.random(200) < 0.5
    out.append(("bimodal-n200", np.where(c, g.normal(-3, 0.5, 200), g.normal(3, 0.5, 200))))
    out.append(("exponential-n150", rng().exponential(1.0, 150)))
    return out


def job(args):
    z, kern = args
    z = np.asarray(z)
    res, s, hlo, hhi = gb.scan_optimum(z, kern, 4000, 50)
    hs = [s * f for f in FACTORS]
    fixed = {c: [] for c in COMP}
    fine = {c: [] for c in COMP}
    for h in hs:
        r = gb.criteria(z, h, kern, 50)
        rf = gb.criteria(z, h, kern, 1000)
        for c in COMP:
            fixed[c].append(r[c])
            fine[c].append(rf[c])
    return res, s, hlo, hhi, hs, fixed, fine


def main():
    ds = datasets()
    jobs = [(list(gb.standardise(x)), k) for _, x in ds for k in gb.KERNELS]
    with ProcessPoolExecutor(max_workers=14) as ex:
        results = list(ex.map(job, jobs, chunksize=1))
    cases = []
    it = iter(results)
    for name, x in ds:
        for kern in gb.KERNELS:
            res, s, hlo, hhi, hs, fixed, fine = next(it)
            case = {Kw("dataset"): Kw(name), Kw("kernel"): Kw(kern), Kw("s"): s, Kw("h-range"): [hlo, hhi]}
            for c, _ in gb.CRIT:
                r = res[c]
                case[Kw(c)] = {Kw("h"): r["h"], Kw("value"): r["value"], Kw("at-bound"): r["at-bound"],
                               Kw("grid-index"): r["grid-index"]}
            case[Kw("fixed-h")] = hs
            case[Kw("fixed")] = {Kw(c): v for c, v in fixed.items()}
            case[Kw("fine")] = {Kw(c): v for c, v in fine.items()}
            cases.append(case)
    dsets = [{Kw("name"): Kw(n), Kw("n"): len(x), Kw("mean"): float(x.mean()), Kw("sd"): float(x.std(ddof=1)),
              Kw("a"): gb.rlcv_a(len(x)), Kw("data"): [float(v) for v in x]} for n, x in ds]
    nb = {c: sum(1 for k in cases if k[Kw(c)][Kw("at-bound")]) for c, _ in gb.CRIT}
    # max |fixed - fine| / |fine| over finite values, per criterion (the trapezoid error of the h/50 grid)
    gap = {}
    for k in cases:
        for c in ("lscv", "rlcv-printed", "rlcv-continuous"):
            for a, b in zip(k[Kw("fixed")][Kw(c)], k[Kw("fine")][Kw(c)]):
                if math.isfinite(a) and math.isfinite(b) and b != 0:
                    gap[c] = max(gap.get(c, 0.0), abs(a - b) / abs(b))
    hdr = [
        ";; rlcv_reference.edn - generated by utils/fastmath/dev/generate_rlcv_reference.py (do not edit)",
        ";; numpy %s, scipy %s. Independent numpy criteria (LSCV, LCV, RLCV printed and continuous) for 14 kernels x 6 datasets; schema in the generator docstring."
        % (np.__version__, scipy.__version__),
        ";; Optima on a bound of the search interval [0.01, 5] s: lscv %d, lcv %d, rlcv-printed %d, rlcv-continuous %d of %d cases."
        % (nb["lscv"], nb["lcv"], nb["rlcv-printed"], nb["rlcv-continuous"], len(cases)),
        ";; Max relative difference of the criteria between the integration spacings h/50 (:fixed) and h/1000 (:fine): lscv %.1e, rlcv-printed %.1e, rlcv-continuous %.1e."
        % (gap["lscv"], gap["rlcv-printed"], gap["rlcv-continuous"]),
    ]
    info = {Kw("grid-points"): 4000, Kw("grid-factors"): [0.01, 5.0], Kw("fixed-factors"): FACTORS,
            Kw("integral-divisors"): [50, 1000], Kw("extent"): {Kw(k): v for k, v in gb.EXT.items()},
            Kw("rlcv-a"): "sqrt(ln n) / (sqrt 2 n)",
            Kw("scale"): "s = 1.06 lo n^-0.2 * nrd-scale (delta0 * (0.5/sqrt(pi))^-0.2), no scale for silverman and cauchy"}
    body = {Kw("cv-info"): info, Kw("datasets"): dsets, Kw("cases"): cases}
    OUT.write_text("\n".join(hdr) + "\n" + gb.emit(body) + "\n")
    print("wrote", OUT, OUT.stat().st_size, "bytes", len(cases), "cases")
    print("at-bound", nb, "h/50 vs h/1000:", gap)


if __name__ == "__main__":
    main()
