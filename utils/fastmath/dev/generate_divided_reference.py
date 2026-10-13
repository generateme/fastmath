"""Generator of reference values for `fastmath.calculus.divided` (Newton and Hermite divided differences).

Tool: sympy 1.14 with 60-digit floating point arithmetic (mpmath backend). The derivatives
f^(k) of every test function are obtained symbolically (`sympy.diff`) and evaluated at the nodes
with 60 digits; the table is then built with the standard confluent recurrence:
    f[x_i .. x_(i+l)] = (f[x_(i+1) .. x_(i+l)] - f[x_i .. x_(i+l-1)]) / (x_(i+l) - x_i)   if x_i != x_(i+l)
    f[x_i .. x_(i+l)] = f^(l)(x_i) / l!                                                   otherwise
(repeated nodes are adjacent, hence equal end points mean that all the nodes between are equal).
The values are converted to the nearest double. Nodes are dyadic rationals (exact doubles).

Domain: the functions exp, sin, cos, log (positive nodes only) and the polynomial
x^5 - 2 x^3 + x on node lists with repeated nodes (a pair, a triple, several runs, one node
repeated six times, runs at both ends) and distinct nodes (two Newton cases: a plain table and nodes of
mixed sign). Why: Hermite values with exact derivatives must agree to 1e-14, and the fallback with
finite-difference derivatives is checked against the same table.

Produces `test/resources/calculus/divided_reference.edn`: a vector of maps
`{:f <id> :xs [...] :levels [[...] ...]}` where `:levels` is the whole table, the level 0 (function
values) first and the single highest level last, as `[(f x0 ... xn)]`. The ids of the functions are
defined identically in `test/fastmath/calculus/divided_test.clj`.

Run from anywhere:
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_divided_reference.py
"""
from pathlib import Path

import sympy
from sympy import Float, Rational, factorial

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "test" / "resources" / "calculus" / "divided_reference.edn"

DIGITS = 60
X = sympy.Symbol("x")

FUNCTIONS = {
    "exp": sympy.exp(X),
    "sin": sympy.sin(X),
    "cos": sympy.cos(X),
    "log": sympy.log(X),
    "poly": X**5 - 2 * X**3 + X,
}

CASES = [
    ("exp", ["0", "0", "1", "1"]),
    ("sin", ["0", "0", "0", "1", "1"]),
    ("exp", ["1/2", "1/2", "1/2", "1/2"]),
    ("cos", ["0", "1", "1", "2"]),
    ("cos", ["0", "0", "1", "1", "1", "2", "2"]),
    ("log", ["1", "1", "2", "3", "3", "3"]),
    ("poly", ["-1", "-1", "0", "0", "2"]),
    ("exp", ["0", "0"]),
    ("exp", ["3/2", "3/2", "3/2", "3/2", "3/2", "3/2"]),
    ("sin", ["-1", "-1", "-1/2", "-1/2", "-1/2", "1", "1"]),
    ("exp", ["2"]),
    ("exp", ["1", "2", "3"]),
    ("sin", ["-1", "-1/2", "0", "1", "3/2"]),
    ("log", ["1", "2", "4", "8"]),
    ("poly", ["-2", "-1", "0", "1", "2", "3"]),
]


def edn_double(x):
    s = repr(float(x))
    assert "inf" not in s and "nan" not in s, s
    return s


def table(fid, xs):
    f = FUNCTIONS[fid]
    xv = [Rational(s) for s in xs]
    n = len(xv)
    cache = {}

    def deriv_over_fact(k, x):
        key = (k, x)
        if key not in cache:
            cache[key] = Float(sympy.diff(f, X, k).subs(X, x).evalf(DIGITS) / factorial(k), DIGITS)
        return cache[key]

    levels = [[deriv_over_fact(0, xv[i]) for i in range(n)]]
    for lvl in range(1, n):
        prev = levels[-1]
        cur = []
        for i in range(n - lvl):
            if xv[i] == xv[i + lvl]:
                cur.append(deriv_over_fact(lvl, xv[i]))
            else:
                cur.append((prev[i + 1] - prev[i]) / Float(xv[i + lvl] - xv[i], DIGITS))
        levels.append(cur)
    return levels


def main():
    lines = ["["]
    for fid, xs in CASES:
        levels = table(fid, xs)
        xs_edn = " ".join(edn_double(Rational(s)) for s in xs)
        lv_edn = " ".join("[" + " ".join(edn_double(v) for v in lvl) + "]" for lvl in levels)
        lines.append(" {:f :%s :xs [%s] :levels [%s]}" % (fid, xs_edn, lv_edn))
    lines.append("]")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n")
    print("sympy", sympy.__version__, "wrote", len(CASES), "cases to", OUT)


if __name__ == "__main__":
    main()
