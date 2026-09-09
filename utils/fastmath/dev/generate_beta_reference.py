"""One-off generator of independent reference values for `beta`,
`regularized-beta` and `incomplete-beta` negative-argument tests in
`fastmath.special-test`.

Produces `test/resources/special/beta_reference.edn`, computed with `mpmath`
(analytic continuation via `mpmath.beta` / `mpmath.betainc`).

Two datasets:
  - `:beta`     {:p :q :ref}                for `beta(p, q)`
  - `:incbeta`  {:x :a :b :reg :inc}         for `regularized-beta`/`incomplete-beta`

Grids intentionally avoid true singularities:
  - `beta`: p, q, and p+q must not be non-positive integers.
  - `incbeta`: a, b, and a+b must not be non-positive integers (a+b being a
    non-positive integer is a genuine, non-removable singularity of the
    *regularized* incomplete beta function, unlike for `beta` itself).

Run with (from repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_beta_reference.py
"""
import mpmath as mp

mp.mp.dps = 30

OUT = "test/resources/special/beta_reference.edn"


def fmt(x):
    return repr(float(x))


def edn_vec(xs):
    return "[" + " ".join(fmt(x) for x in xs) + "]"


def is_nonpos_int(x, tol=1e-6):
    # tolerance-based: also excludes points landing *near* a pole, where
    # Gamma(x) blows up and reference/implementation both lose precision
    return x <= tol and abs(x - round(x)) < tol


# ---- beta(p, q) grid ----
p_vals = [-(0.3 + 0.5 * i) for i in range(0, 20)] + [1.7, 3.2]   # negative + a couple positive
q_vals = [-(0.2 + 0.7 * i) for i in range(0, 20)] + [2.1, 4.4]

beta_p, beta_q, beta_ref = [], [], []
for p in p_vals:
    for q in q_vals:
        s = p + q
        if is_nonpos_int(p) or is_nonpos_int(q) or is_nonpos_int(s):
            continue
        beta_p.append(p)
        beta_q.append(q)
        beta_ref.append(mp.beta(p, q))

# ---- regularized-beta / incomplete-beta grid ----
x_vals = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95]
a_vals = [-(0.3 + 0.6 * i) for i in range(0, 8)] + [1.3, 2.7]
b_vals = [-(0.4 + 0.5 * i) for i in range(0, 8)] + [1.9, 3.1]

ib_x, ib_a, ib_b, ib_reg, ib_inc = [], [], [], [], []
for x in x_vals:
    for a in a_vals:
        for b in b_vals:
            s = a + b
            if is_nonpos_int(a) or is_nonpos_int(b) or is_nonpos_int(s):
                continue
            ib_x.append(x)
            ib_a.append(a)
            ib_b.append(b)
            ib_reg.append(mp.betainc(a, b, 0, x, regularized=True))
            ib_inc.append(mp.betainc(a, b, 0, x, regularized=False))

with open(OUT, "w") as f:
    f.write("{:beta {:p " + edn_vec(beta_p) + "\n")
    f.write("        :q " + edn_vec(beta_q) + "\n")
    f.write("        :ref " + edn_vec(beta_ref) + "}\n")
    f.write(" :incbeta {:x " + edn_vec(ib_x) + "\n")
    f.write("           :a " + edn_vec(ib_a) + "\n")
    f.write("           :b " + edn_vec(ib_b) + "\n")
    f.write("           :reg " + edn_vec(ib_reg) + "\n")
    f.write("           :inc " + edn_vec(ib_inc) + "}}\n")

print("beta grid points:", len(beta_p))
print("incbeta grid points:", len(ib_x))
print("wrote", OUT)
