"""Generator of reference values for `fastmath.calculus.samples/integrate-samples`.

Tool: scipy 1.18.1 (numpy 2.x), Python 3.14.
- `scipy.integrate.trapezoid`          (method :trapezoid, not cumulative),
- `scipy.integrate.simpson`            (method :simpson, not cumulative; the modern implementation, even N
                                        uses the Cartwright correction for the last interval, x may be non-uniform),
- `scipy.integrate.cumulative_trapezoid(initial=0)` (method :trapezoid, cumulative),
- `scipy.integrate.cumulative_simpson(initial=0)`   (method :simpson, cumulative; N < 3 falls back to the
                                        trapezoid rule).
`cumulative_simpson` rejects a non-increasing `x` ("Input x must be strictly increasing") though its formulas are
valid for a strictly decreasing `x`; for the decreasing-x cases the same scipy internals
(`_cumulatively_sum_simpson_integrals` with `_cumulative_simpson_unequal_intervals`) are called directly,
skipping only the check. `fastmath` does not impose that restriction.

Domain: sample counts N in {2, 3, 4, 5, 6, 7, 10, 11, 50, 101} (odd and even, both smallest sizes for
Simpson), random standard-normal samples on a uniform grid (spacing as a number `dx`, and the same grid as an
explicit `x` sequence), on non-uniform grids (random positive steps), with a negative `dx` and a decreasing `x`.
Additionally: polynomials where Simpson is exact (a cubic on a uniform odd grid, a quadratic on a non-uniform
grid; `:exact` is computed with exact rational arithmetic), and a few hand cases. All random data comes
from fixed seeds, the output is deterministic.
Why: the rules differ for odd/even N, for uniform/non-uniform grids and for the cumulative variants, each of
which is a separate code path.

N = 0 and N = 1 are not here (scipy raises or returns 0.0 depending on the function; fastmath defines them and
tests them by hand).

Produces `test/resources/calculus/samples_reference.edn`: a map of groups, each a vector of cases
`{:id ... :ys [...] :x [...] | :dx d, :trapezoid r, :trapezoid-cumulative [...], :simpson r,
:simpson-cumulative [...] (:exact r)}`.

Run from anywhere:
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_samples_reference.py
"""
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy
from scipy import integrate
from scipy.integrate import _quadrature as _q

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "test" / "resources" / "calculus" / "samples_reference.edn"

NS = [2, 3, 4, 5, 6, 7, 10, 11, 50, 101]


def fl(v):
    return repr(float(v))


def vec(vs):
    return "[" + " ".join(fl(v) for v in vs) + "]"


def cumulative_simpson(y, x=None, dx=1.0):
    """Public scipy function, or the same internals without the increasing-x check."""
    try:
        return integrate.cumulative_simpson(y, x=x, dx=dx, initial=0)
    except ValueError as e:
        assert "strictly increasing" in str(e) and x[0] > x[-1], e
        d = np.diff(x)
        assert np.all(d < 0)
        res = _q._cumulatively_sum_simpson_integrals(
            y, d, _q._cumulative_simpson_unequal_intervals, np)
        return np.concatenate(([0.0], res))


def refs(y, x=None, dx=None):
    y = np.asarray(y, dtype=float)
    out = {}
    if x is not None:
        x = np.asarray(x, dtype=float)
        out[":trapezoid"] = integrate.trapezoid(y, x=x)
        out[":trapezoid-cumulative"] = integrate.cumulative_trapezoid(y, x=x, initial=0)
        out[":simpson"] = integrate.simpson(y, x=x)
        out[":simpson-cumulative"] = cumulative_simpson(y, x=x)
    else:
        out[":trapezoid"] = integrate.trapezoid(y, dx=dx)
        out[":trapezoid-cumulative"] = integrate.cumulative_trapezoid(y, dx=dx, initial=0)
        out[":simpson"] = integrate.simpson(y, dx=dx)
        out[":simpson-cumulative"] = cumulative_simpson(y, dx=dx)
    return out


def case(cid, y, x=None, dx=None, exact=None):
    r = refs(y, x, dx)
    parts = [f':id "{cid}"', f":ys {vec(y)}"]
    parts.append(f":x {vec(x)}" if x is not None else f":dx {fl(dx)}")
    for k in (":trapezoid", ":simpson"):
        parts.append(f"{k} {fl(r[k])}")
    for k in (":trapezoid-cumulative", ":simpson-cumulative"):
        assert len(r[k]) == len(y)
        parts.append(f"{k} {vec(r[k])}")
    if exact is not None:
        parts.append(f":exact {fl(exact)}")
    return "{" + "\n   ".join(parts) + "}"


def group(name, cases):
    return f" {name}\n [" + "\n  ".join(cases) + "]"


def poly_exact(coeffs, a, b):
    """Exact integral of sum c_k x^k on [a, b] (Fractions)."""
    return sum(Fraction(c) * (Fraction(b) ** (k + 1) - Fraction(a) ** (k + 1)) / (k + 1)
               for k, c in enumerate(coeffs))


def poly_val(coeffs, xs):
    return [float(sum(Fraction(c) * Fraction(float(x)) ** k for k, c in enumerate(coeffs))) for x in xs]


def main():
    assert scipy.__version__ == "1.18.1", scipy.__version__
    rng = np.random.default_rng(20240607)

    uniform_dx, uniform_x, nonuniform, negative = [], [], [], []
    for n in NS:
        y = rng.standard_normal(n)
        dx = float(rng.uniform(0.05, 0.5))
        uniform_dx.append(case(f"dx-n{n}", y, dx=dx))
        y = rng.standard_normal(n)
        x = np.linspace(-1.5, 2.5, n)
        uniform_x.append(case(f"linspace-n{n}", y, x=x))
        y = rng.standard_normal(n)
        x = 0.3 + np.cumsum(np.concatenate(([0.0], rng.uniform(0.05, 1.0, n - 1))))
        nonuniform.append(case(f"random-x-n{n}", y, x=x))
    for n in NS:
        y = rng.standard_normal(n)
        negative.append(case(f"negative-dx-n{n}", y, dx=-float(rng.uniform(0.05, 0.5))))
        y = rng.standard_normal(n)
        x = 5.0 - np.cumsum(np.concatenate(([0.0], rng.uniform(0.05, 1.0, n - 1))))
        negative.append(case(f"decreasing-x-n{n}", y, x=x))

    known = []
    # cubic, uniform grid, odd N: Simpson is exact
    cubic = [Fraction(1, 2), Fraction(-3), Fraction(2), Fraction(1, 4)]
    for n in (3, 5, 11, 101):
        x = np.linspace(-1.0, 2.0, n)
        known.append(case(f"cubic-uniform-n{n}", poly_val(cubic, x), x=x,
                          exact=poly_exact(cubic, -1.0, 2.0)))
    x = np.linspace(-1.0, 2.0, 11)
    known.append(case("cubic-dx-n11", poly_val(cubic, x), dx=float(x[1] - x[0]),
                      exact=poly_exact(cubic, -1.0, 2.0)))
    # quadratic, non-uniform grid, odd N: Simpson is exact
    quad = [Fraction(2), Fraction(-1, 3), Fraction(5, 2)]
    for n in (3, 7, 11):
        x = np.sort(np.concatenate(([0.0, 3.0], rng.uniform(0.1, 2.9, n - 2))))
        known.append(case(f"quadratic-nonuniform-n{n}", poly_val(quad, x), x=x,
                          exact=poly_exact(quad, 0.0, 3.0)))
    # constant and linear data (every rule exact, even N included)
    for n in (2, 4, 6):
        x = np.linspace(0.0, 1.0, n)
        known.append(case(f"linear-n{n}", 3.0 * x + 1.0, x=x, exact=2.5))
    # hand cases
    known.append(case("hand-x-eq-y-0-9", np.arange(10.0), x=np.arange(10.0)))        # simpson 40.5
    known.append(case("hand-cube-0-9", np.arange(10.0) ** 3, x=np.arange(10.0)))     # simpson 1640.5
    known.append(case("hand-two-points", [1.0, 5.0], dx=2.0))
    known.append(case("hand-three-points", [1.0, 4.0, 9.0], x=[0.0, 1.0, 3.0]))
    known.append(case("hand-four-points", [1.0, 2.0, 4.0, 8.0], x=[0.0, 1.0, 2.0, 4.0]))
    known.append(case("hand-zeros", [0.0] * 5, dx=1.0))

    text = "{" + "\n".join([
        group(":uniform-dx", uniform_dx),
        group(":uniform-x", uniform_x),
        group(":nonuniform-x", nonuniform),
        group(":negative", negative),
        group(":known", known),
    ]) + "}\n"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(text)
    print("wrote", OUT, len(text), "bytes")


if __name__ == "__main__":
    main()
