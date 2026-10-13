"""Generator of reference values for `fastmath.calculus.finite/fd-coeffs` and `fd-coeffs-for-accuracy`.

Tool: sympy (`sympy.calculus.finite_diff.finite_diff_weights`, Fornberg's algorithm in exact
rational arithmetic); the exact rational weights are converted to the nearest double
(`float(Fraction(p, q))`, correctly rounded).

Domain: all 168 combinations of the derivative order n = 1..7, the order of accuracy acc = 1..8 and
the stencil kind in (central, forward, backward). The offsets of the stencil are built with the same
rule as `fd-coeffs-for-accuracy`:
    num-central = 2 * floor((n + 1) / 2) - 1 + acc,   num-side = floor(num-central / 2),
    num-coeffs  = num-central + 1 if n is even else num-central
    forward  : 0 .. num-coeffs - 1
    backward : -(num-coeffs - 1) .. 0
    central  : -num-side .. num-side
Why: the former implementation solved a Vandermonde system by LU, which is ill-conditioned (51 of
these 168 combinations differed by more than 1e-9, worst 202 for n = 7, acc = 8, forward); the
whole grid is the range of the orders the library handles with an inline stencil (n <= 6,
acc <= 6) and with the general path.

Produces `test/resources/calculus/fd_coeffs_reference.edn`: a vector of maps
`{:n :acc :kind :offsets [...] :coeffs [...]}` where `:coeffs` are the weights of the n-th
derivative at 0 (to be divided by h^n).

Run from anywhere:
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_fd_coeffs_reference.py
"""
from fractions import Fraction
from pathlib import Path

import sympy
from sympy.calculus.finite_diff import finite_diff_weights

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "test" / "resources" / "calculus" / "fd_coeffs_reference.edn"

ORDERS = range(1, 8)
ACCURACIES = range(1, 9)
KINDS = ("central", "forward", "backward")


def offsets_for(n, acc, kind):
    num_central = 2 * ((n + 1) // 2) - 1 + acc
    num_side = num_central // 2
    num_coeffs = num_central + 1 if n % 2 == 0 else num_central
    if kind == "forward":
        return list(range(num_coeffs))
    if kind == "backward":
        return list(range(-(num_coeffs - 1), 1))
    return list(range(-num_side, num_side + 1))


def to_double(r):
    return float(Fraction(int(r.p), int(r.q)))


def edn_double(x):
    s = repr(float(x))
    assert "inf" not in s and "nan" not in s, s
    return s


def main():
    rows = []
    for n in ORDERS:
        for acc in ACCURACIES:
            for kind in KINDS:
                offsets = offsets_for(n, acc, kind)
                weights = finite_diff_weights(n, offsets, 0)[n][-1]
                assert len(weights) == len(offsets)
                coeffs = [to_double(sympy.Rational(w)) for w in weights]
                rows.append((n, acc, kind, offsets, coeffs))
    lines = ["["]
    for n, acc, kind, offsets, coeffs in rows:
        lines.append(
            " {:n %d :acc %d :kind :%s :offsets [%s] :coeffs [%s]}"
            % (n, acc, kind, " ".join(str(o) for o in offsets), " ".join(edn_double(c) for c in coeffs))
        )
    lines.append("]")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n")
    print("sympy", sympy.__version__, "wrote", len(rows), "cases to", OUT)


if __name__ == "__main__":
    main()
