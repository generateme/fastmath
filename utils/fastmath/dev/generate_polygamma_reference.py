"""One-off generator of independent reference values for the `polygamma`
test in `fastmath.special-test`.

Produces `test/resources/special/polygamma_reference.edn`, computed with
`mpmath.polygamma` (an entirely independent implementation from fastmath's).

Domains were chosen empirically (see investigation notes in the PR/session):
  - `x>0`: safe and accurate for the *entire* order range tested here (up to
    m=300), after fixing a premature-`Gamma`-overflow bug in `polygamma`
    (same root cause, and same log-space fix, as an earlier `xi` fix).
  - `x<=0`: restricted to small orders (`m<=6`) and moderate `|x|` (<=10).
    `polygamma` at nonpositive `x` goes through an internal `cotderiv`
    helper (the m-th derivative of `cot(pi*z)`) which, for even `m`, is
    mathematically exactly zero at half-integer `z` (an odd-function
    symmetry) but is computed from a merely-double-precision `pi`; the
    resulting tiny nonzero floating-point residual gets amplified by
    `pi^(m+1)`, corrupting results starting around `m~12` near
    half-integer `x` and spreading to a wider domain for larger `m` -- a
    known, documented, NOT fixed in this session, limitation (see
    `polygamma`'s docstring). The chosen (m<=6, |x|<=10) box was verified
    point-by-point against this reference generator to be free of the
    issue.

Run with (from repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_polygamma_reference.py
"""
import mpmath as mp

mp.mp.dps = 40

OUT = "test/resources/special/polygamma_reference.edn"


def to_double(x):
    d = float(x)
    if d != d or d in (float('inf'), float('-inf')):
        return None
    return d


def fmt(x):
    return repr(x)


def edn_vec(xs):
    return "[" + " ".join(fmt(x) for x in xs) + "]"


# ---- x>0: wide order range, safe & accurate after the overflow fix ----
pos_ms = list(range(2, 11)) + [20, 30, 50, 80, 100, 120, 130, 140, 150, 160, 170, 200, 250, 300]
pos_xs = [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 5.0, 10.0, 50.0, 100.0, 1e4, 1e6, 1e10]

pos_m_out, pos_x_out, pos_ref = [], [], []
for m in pos_ms:
    for x in pos_xs:
        # fastmath's zeta(s,x) internally computes a leading term x^-s in
        # double precision; if that itself underflows to exactly 0.0 (or
        # into the reduced-precision subnormal range), regardless of how
        # huge the true polygamma value ends up being once multiplied by
        # Gamma(m+1), fastmath returns 0.0 or a low-precision result -- a
        # known, documented double-precision limitation, not a bug. Skip
        # such points here rather than asserting a value fastmath cannot
        # produce accurately this way.
        if x ** (-(m + 1)) < 2.2250738585072014e-308:  # smallest normal double
            continue
        try:
            v = to_double(mp.polygamma(m, mp.mpf(x)))
        except Exception as e:
            print("skip", m, x, e)
            continue
        if v is not None:
            pos_m_out.append(m)
            pos_x_out.append(x)
            pos_ref.append(v)

# ---- x<=0: restricted to small orders / moderate |x| (see module docstring) ----
neg_ms = [2, 3, 4, 5, 6]
neg_xs = [-0.1, -0.25, -0.5, -0.75, -0.9, -1.1, -1.5, -2.5, -3.5, -5.5, -8.5, -10.0]

neg_m_out, neg_x_out, neg_ref = [], [], []
for m in neg_ms:
    for x in neg_xs:
        try:
            v = to_double(mp.polygamma(m, mp.mpf(x)))
        except Exception as e:
            print("skip", m, x, e)
            continue
        if v is not None:
            neg_m_out.append(m)
            neg_x_out.append(x)
            neg_ref.append(v)

with open(OUT, "w") as f:
    f.write("{:pos {:order " + edn_vec(pos_m_out) + "\n")
    f.write("       :arg " + edn_vec(pos_x_out) + "\n")
    f.write("       :ref " + edn_vec(pos_ref) + "}\n")
    f.write(" :neg {:order " + edn_vec(neg_m_out) + "\n")
    f.write("       :arg " + edn_vec(neg_x_out) + "\n")
    f.write("       :ref " + edn_vec(neg_ref) + "}}\n")

print("pos points:", len(pos_m_out))
print("neg points:", len(neg_m_out))
print("wrote", OUT)
