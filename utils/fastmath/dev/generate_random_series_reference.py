"""Generator of reference data for the time-series functions of `fastmath.random` (Random Audit, work package 9.5):
`ma ar arma fi arfima` (and their private helpers `fi-diffs`). None of the oracles reads fastmath.

Oracles
    :lfilter       `scipy.signal.lfilter` with zero initial conditions: ma `b = [1 theta]`, ar `a = [1 -phi]`, arma both. Orders 1-4.
    :emulation     independent Python emulation of the library recursion (Spec 4.1/4.6): "the history is seeded with the first cnt
                   signal values, the recursion runs from index cnt (ar: Y_t = sum phi_i Y_{t-i} + e_t on the raw seeded history;
                   ma: Y_t = e_t + sum theta_i e_{t-i}), the first `burn` outputs after the seed are dropped" (default burn = cnt):
                   output j = process value at signal index cnt + burn + j. For ma it equals lfilter exactly; for ar of order >= 2 and for
                   arma it is NOT the zero-initial-condition filter (the raw-seeded start), the difference to :lfilter decays geometrically.
                   arma = ma stage then ar stage, each with its own burn.
    :arima-sim     R `stats::arima.sim` with explicit `innov`, `start.innov` and `n.start = 2 (p + q)` (generate_random_series_reference.R,
                   called through Rscript).
    :fi-exact      psi weights psi_k = Gamma(k + d) / (Gamma(d) k!) by mpmath (50 digits); y_t = sum_{k<=t} psi_k e_{t-k} (zero initial
                   conditions, mpmath arithmetic, rounded to double); the library output j is y_{j+1} (e_0 is the seed).
    :fi-diffs      coefficients c_k = -(-1)^k binom(d, k), k = 1..50 (mpmath `binomial`), the AR coefficients of the truncated filter.
    :fi-trunc      the raw-seeded AR(L) emulation with the coefficients c_k.
    :fi-int        numpy: d >= 1 `cumsum` applied d times; d <= -1 `np.diff(x, prepend=0)` applied |d| times (full length series).
    :arfima        3-stage emulation ma -> fi (exact or truncated) -> ar.
    :theory        statsmodels `ArmaProcess` (acovf, acf up to lag 4000): variance, acf lags 0..5, Bartlett standard errors
                   sqrt(sum_{j>=1} (rho_{j+k} + rho_{j-k} - 2 rho_j rho_k)^2) (per sqrt(N)) of the sample acf and
                   sqrt(2 sum_{j in Z} rho_j^2) of the relative variance.

Signal: numpy PCG64(20241021) standard normal, N = 120 (stored in the file; the tests use prefixes, outputs are causal).

Run:  cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_random_series_reference.py
      (needs Rscript on the path; writes /tmp/fastmath_series_signal.txt and /tmp/fastmath_series_r.txt)
Output: test/resources/random/series_reference.edn (deterministic, byte-identical on rerun).
"""
import os
import platform
import subprocess

import mpmath as mp
import numpy as np
import scipy
import statsmodels
from scipy.signal import lfilter
from statsmodels.tsa.arima_process import ArmaProcess

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "..", "test", "resources", "random", "series_reference.edn")
SIGNAL_TXT = "/tmp/fastmath_series_signal.txt"
R_TXT = "/tmp/fastmath_series_r.txt"

mp.mp.dps = 50
N = 120

MODELS = [  # name, phi, theta (same table as in the R script)
    ("ma-1", [], [0.5]),
    ("ma-2", [], [0.5, -0.3]),
    ("ma-3", [], [0.4, -0.3, 0.2]),
    ("ma-4", [], [0.5, -0.3, 0.2, 0.1]),
    ("ar-1", [0.6], []),
    ("ar-2", [0.5, -0.3], []),
    ("ar-3", [0.2, 0.3, 0.1], []),
    ("ar-4", [0.3, 0.2, 0.1, 0.05], []),
    ("arma-11", [0.7], [0.4]),
    ("arma-22", [0.5, -0.2], [0.3, 0.1]),
    ("arma-33", [0.2, 0.3, 0.1], [0.4, -0.3, 0.2]),
    ("arma-44", [0.3, 0.2, 0.1, 0.05], [0.5, -0.3, 0.2, 0.1]),
    ("arma-21", [0.5, -0.3], [0.4]),
    ("arma-12", [0.5], [0.4, 0.2]),
]

BURN_MODELS = ["ma-2", "ar-2", "ar-3", "arma-22", "arma-12"]
BURNS = [0, 1, 5, 20]

FI_D = [0.1, 0.3, 0.45, -0.3, 0.7, 1.3]
FI_TRUNC = [(0.2, 1), (0.3, 5), (0.4, 20), (-0.3, 10)]
FI_INT = [1, 2, 3, -1, -2]

THEORY = [  # name, phi, theta
    ("ar-1", [0.7], []),
    ("ar-2", [0.5, -0.3], []),
    ("ma-2", [], [0.5, -0.4]),
    ("arma-11", [0.7], [0.4]),
]


def num(x):
    return repr(float(x))


def vec(xs):
    return "[" + " ".join(num(x) for x in xs) + "]"


def ma_emul(x, theta, burn=None):
    q = len(theta)
    burn = q if burn is None else burn
    y = [x[t] + sum(theta[i] * x[t - 1 - i] for i in range(q)) for t in range(q, len(x))]
    return np.array(y[burn:])


def ar_emul(x, phi, burn=None):
    p = len(phi)
    burn = p if burn is None else burn
    y = list(x[:p])
    for t in range(p, len(x)):
        y.append(sum(phi[i] * y[t - 1 - i] for i in range(p)) + x[t])
    return np.array(y[p + burn:])


def arma_emul(x, phi, theta, burn_ma=None, burn_ar=None):
    z = ma_emul(x, theta, burn_ma) if theta else x
    return ar_emul(z, phi, burn_ar) if phi else z


def psi(d, n):
    d = mp.mpf(d)
    return [mp.gamma(k + d) / (mp.gamma(d) * mp.factorial(k)) for k in range(n)]


def fi_exact(x, d):
    """y_t = sum_k psi_k x_{t-k} in mpmath, rounded to double; full length."""
    ps = psi(d, len(x))
    xm = [mp.mpf(float(v)) for v in x]
    return np.array([float(sum(ps[k] * xm[t - k] for k in range(t + 1))) for t in range(len(x))])


def fi_coeffs(d, n):
    d = mp.mpf(d)
    return [float(-(-1) ** k * mp.binomial(d, k)) for k in range(1, n + 1)]


def read_r(path):
    echo, ver, models = None, None, {}
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("E "):
            echo = [float(v) for v in line[2:].split()]
        elif line.startswith("V "):
            ver = line[2:]
        elif line.startswith("M "):
            name, rest = line[2:].split(" ", 1)
            parts = rest.split(";")
            models[name] = ([float(v) for v in parts[0].split()], [float(v) for v in parts[1].split()],
                            [float(v) for v in parts[2].split()])
    return echo, ver, models


def theory(phi, theta):
    ap = ArmaProcess(np.r_[1.0, -np.array(phi)], np.r_[1.0, np.array(theta)])
    k = 4000
    acov = ap.acovf(k)
    rho = acov / acov[0]
    se = []
    for lag in range(1, 6):
        s = 0.0
        for j in range(1, k - lag - 1):
            s += (rho[j + lag] + rho[abs(j - lag)] - 2 * rho[j] * rho[lag]) ** 2
        se.append(float(np.sqrt(s)))
    var_rel = float(np.sqrt(2 * (1 + 2 * np.sum(rho[1:k] ** 2))))
    return float(acov[0]), [float(r) for r in rho[:6]], se, var_rel


def main():
    rng = np.random.Generator(np.random.PCG64(20241021))
    e = rng.standard_normal(N)
    with open(SIGNAL_TXT, "w") as fh:
        fh.write("\n".join(repr(float(v)) for v in e) + "\n")
    subprocess.run(["Rscript", os.path.join(HERE, "generate_random_series_reference.R"), SIGNAL_TXT, R_TXT], check=True)
    echo, rver, rmodels = read_r(R_TXT)
    assert echo == [float(v) for v in e], "R text round trip of the signal is not exact"

    out = []
    out.append(':header {:tool "python %s numpy %s scipy %s mpmath %s statsmodels %s; %s"\n'
               '  :precision "lfilter/numpy/R: IEEE double (R text %%.17g round trip checked); psi weights and fi-diffs: mpmath 50 digits; '
               'theory: statsmodels acovf to lag 4000, Bartlett sums to lag 4000"\n'
               '  :domain "signal PCG64(20241021) standard normal N = %d; ma ar arma orders 1-4 (stable ar parts); fi exact d 0.1 0.3 0.45 -0.3 0.7 1.3 '
               '(N = 120), fi-diffs k <= 50, truncated fi, integer fi d 1 2 3 -1 -2 (N = 60), burn-in 0 1 5 20 (arfima 0 3), 4 ARMA theory models"}'
               % (platform.python_version(), np.__version__, scipy.__version__, mp.__version__, statsmodels.__version__, rver, N))
    out.append(":signal " + vec(e))

    cases = []
    for name, phi, theta in MODELS:
        b = [1.0] + list(theta)
        a = [1.0] + [-v for v in phi]
        lf = lfilter(b, a, e)
        em = arma_emul(e, phi, theta)
        rm = rmodels[name]
        assert rm[0] == phi and rm[1] == theta, name
        cases.append(" {:id :%s :phi %s :theta %s\n  :lfilter %s\n  :emulation %s\n  :arima-sim %s}" % (
            name, vec(phi), vec(theta), vec(lf), vec(em), vec(rm[2])))
    out.append(":filters [\n%s]" % "\n".join(cases))

    bs = []
    for name, phi, theta in MODELS:
        if name not in BURN_MODELS:
            continue
        for burn in BURNS:
            if not theta:
                em = ar_emul(e, phi, burn)
            elif not phi:
                em = ma_emul(e, theta, burn)
            else:
                em = arma_emul(e, phi, theta, burn, burn)
            bs.append(" {:id :%s :burn-in %d :emulation %s}" % (name, burn, vec(em)))
    out.append(":burn-in [\n%s]" % "\n".join(bs))

    fe = []
    for d in FI_D:
        y = fi_exact(e, d)
        fe.append(" %s {:psi %s\n  :out %s}" % (num(d), vec([float(p) for p in psi(d, 12)]), vec(y[1:])))
    out.append(":fi-exact {\n%s}" % "\n".join(fe))

    out.append(":fi-diffs {\n%s}" % "\n".join(" %s %s" % (num(d), vec(fi_coeffs(d, 50))) for d in FI_D))

    ft = []
    for d, lim in FI_TRUNC:
        c = fi_coeffs(d, lim)
        ft.append(" {:d %s :dlimit %d :emulation %s :emulation-burn-0 %s}" % (num(d), lim, vec(ar_emul(e, c)), vec(ar_emul(e, c, 0))))
    out.append(":fi-trunc [\n%s]" % "\n".join(ft))

    fi_int = []
    e60 = e[:60]
    for d in FI_INT:
        y = e60.copy()
        for _ in range(abs(d)):
            y = np.cumsum(y) if d > 0 else np.diff(y, prepend=0.0)
        fi_int.append(" %d %s" % (d, vec(y)))
    out.append(":fi-int {\n%s}" % "\n".join(fi_int))

    af = []
    # exact: phi [0.5], d 0.2, theta [0.3]
    x = ma_emul(e, [0.3])
    z = fi_exact(x, 0.2)[1:]
    af.append(" {:phi [0.5] :d 0.2 :dlimit 0 :theta [0.3] :emulation %s}" % vec(ar_emul(z, [0.5])))
    # truncated: phi [0.5 0.1], d 0.3, dlimit 20, theta [0.2]
    x = ma_emul(e, [0.2])
    z = ar_emul(x, fi_coeffs(0.3, 20))
    af.append(" {:phi [0.5 0.1] :d 0.3 :dlimit 20 :theta [0.2] :emulation %s}" % vec(ar_emul(z, [0.5, 0.1])))
    # truncated with burn-in b in every stage: phi [0.5 0.1], d 0.3, dlimit 20, theta [0.2]
    for burn in (0, 3):
        x = ma_emul(e, [0.2], burn)
        z = ar_emul(x, fi_coeffs(0.3, 20), burn)
        af.append(" {:phi [0.5 0.1] :d 0.3 :dlimit 20 :theta [0.2] :burn-in %d :emulation %s}" % (burn, vec(ar_emul(z, [0.5, 0.1], burn))))
    out.append(":arfima [\n%s]" % "\n".join(af))

    th = []
    for name, phi, theta in THEORY:
        var, acf, se, vse = theory(phi, theta)
        th.append(" {:id :%s :phi %s :theta %s :variance %s\n  :acf %s\n  :acf-se-sqrt-n %s :variance-rel-se-sqrt-n %s}" % (
            name, vec(phi), vec(theta), num(var), vec(acf), vec(se), num(vse)))
    out.append(":theory [\n%s]" % "\n".join(th))

    with open(OUT, "w") as fh:
        fh.write("{" + "\n".join(out) + "}\n")


main()
