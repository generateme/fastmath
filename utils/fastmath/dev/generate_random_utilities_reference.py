"""Generator of reference data for `fastmath.random.utilities-test` (Random Audit WP 7.3, 9.6).

Tools: mpmath 1.3.0 (exact negative binomial pmf/cdf for the `nbi`/`nbii` fallback region), R gamlss.dist
(`dNBI/pNBI/qNBI/dNBII/pNBII/qNBII`, run through `generate_random_utilities_reference.R`, version recorded in the
output under `:versions`).  No random numbers are used: the output is deterministic.

Precision: mpmath with 60 digits, rounded once to double.  The `mu` and `sigma` arguments are exact doubles (printed with
`repr`) and the reference is the exact law for exactly these doubles: `nbi` has r = 1/sigma, p = r/(r+mu); `nbii` has
r = mu/sigma, p = 1/(1+sigma) (so that mean = mu and variance = mu + sigma mu^2, resp. mu (1+sigma)).  The recursion
pmf(k+1) = pmf(k) (k+r)/(k+1) (1-p), pmf(0) = p^r is evaluated in 60 digits, error far below 1e-40 for k <= 2e4.

Domain:
    :fallback    sigma in 1e-9 1e-6 1e-4 5e-5 1e-3 x mu in 1 100 1e4, both families; k = round(mu + j sd) for
                 j in -6 -4 -3 -2 -1 0 1 2 3 4 6 (non-negative, pmf > 1e-200) and 0 1 2; exact pmf, log-pmf, cdf, mean,
                 variance.
    :gamlss      see the R script: pmf, log-pmf, cdf, sf, quantiles for six (mu, sigma) per family.

Run (writes test/resources/random/utilities_reference.edn, byte-identical on rerun):
    cd /home/ts/penv && uv run python /home/ts/clojure/fastmath/utils/fastmath/dev/generate_random_utilities_reference.py
"""
import json
import os
import subprocess
import tempfile

import mpmath as mp

mp.mp.dps = 60

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "..", "test", "resources", "random", "utilities_reference.edn")
RSCRIPT = os.path.join(HERE, "generate_random_utilities_reference.R")


# ---------------------------------------------------------------- EDN writer

def edn(v, indent=0):
    if isinstance(v, bool):
        return "true" if v else "false"
    if v is None:
        return "nil"
    if isinstance(v, str):
        return v if v.startswith(":") else '"' + v + '"'
    if isinstance(v, int):
        return str(v)
    if isinstance(v, (float, mp.mpf)):
        f = float(v)
        if f != f:
            return "##NaN"
        if f == float("inf"):
            return "##Inf"
        if f == float("-inf"):
            return "##-Inf"
        return repr(f)
    if isinstance(v, (list, tuple)):
        return "[" + " ".join(edn(x) for x in v) + "]"
    if isinstance(v, dict):
        pad = " " * (indent + 1)
        items = [edn(k) + " " + edn(x, indent + 1) for k, x in v.items()]
        return "{" + ("\n" + pad).join(items) + "}"
    raise TypeError(type(v))


# ---------------------------------------------------------------- exact negative binomial

def nb_exact(r, p, ks):
    """pmf, log pmf and cdf of NB(r, p) (failures before the r-th success) at every k of `ks`, 60 digits."""
    q = 1 - p
    kmax = max(ks)
    want = set(ks)
    pmf_k = mp.power(p, r)
    cdf_k = mp.mpf(0)
    res = {}
    for k in range(kmax + 1):
        cdf_k += pmf_k
        if k in want:
            res[k] = (pmf_k, cdf_k)
        pmf_k = pmf_k * (k + r) / (k + 1) * q
    return res


def fallback_case(family, mu, sigma):
    mu_m, sigma_m = mp.mpf(mu), mp.mpf(sigma)
    if family == "nbi":
        r = 1 / sigma_m
        p = r / (r + mu_m)
        var = mu_m + sigma_m * mu_m ** 2
    else:
        r = mu_m / sigma_m
        p = 1 / (1 + sigma_m)
        var = mu_m * (1 + sigma_m)
    sd = mp.sqrt(var)
    ks = {0, 1, 2}
    for j in (-6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6):
        k = int(mp.nint(mu_m + j * sd))
        if k >= 0:
            ks.add(k)
    ks = sorted(ks)
    exact = nb_exact(r, p, ks)
    kk = [k for k in ks if exact[k][0] > mp.mpf("1e-200")]
    return {"mu": mu, "sigma": sigma, "ks": kk,
            "pmf": [float(exact[k][0]) for k in kk],
            "lpmf": [float(mp.log(exact[k][0])) for k in kk],
            "cdf": [float(exact[k][1]) for k in kk],
            "mean": mu, "var": float(var)}


def main():
    fallback = {"nbi": [], "nbii": []}
    for family in ("nbi", "nbii"):
        for sigma in (1e-9, 1e-6, 1e-4, 5e-5, 1e-3):
            for mu in (1.0, 100.0, 1e4):
                fallback[family].append(fallback_case(family, mu, sigma))

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "gamlss.json")
        subprocess.run(["Rscript", RSCRIPT, path], check=True, capture_output=True)
        with open(path) as f:
            rj = json.load(f)

    def num(x):
        return float(x)

    def conv(case):
        return {k: ([num(x) for x in v] if isinstance(v, list) else num(v)) for k, v in case.items()}

    gamlss = {"nbi": [conv(c) for c in rj["nbi"]], "nbii": [conv(c) for c in rj["nbii"]]}
    out = {":versions": {":mpmath": mp.__version__, ":R": rj["version"]["R"], ":gamlss.dist": rj["version"]["gamlss.dist"]},
           ":fallback": {":" + k: v for k, v in fallback.items()},
           ":gamlss": {":" + k: v for k, v in gamlss.items()}}

    def keyed(v):
        """Turn the dict keys of the nested structure into EDN keywords (the values keep their types)."""
        if isinstance(v, dict):
            return {(k if k.startswith(":") else ":" + k): keyed(x) for k, x in v.items()}
        if isinstance(v, list):
            return [keyed(x) for x in v]
        return v

    text = ";; Reference data for fastmath.random.utilities-test. Generated by utils/fastmath/dev/generate_random_utilities_reference.py\n"
    text += edn(keyed(out)) + "\n"
    with open(OUT, "w") as f:
        f.write(text)
    print("written", os.path.normpath(OUT))


if __name__ == "__main__":
    main()
