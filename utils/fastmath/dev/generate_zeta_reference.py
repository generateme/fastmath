"""One-off generator of independent reference values for `zeta`, `eta`,
`dirichlet-beta` and `xi` tests in `fastmath.special-test`.

Produces `test/resources/special/zeta_reference.edn`, computed with `mpmath`
(`mpmath.zeta`, `mpmath.altzeta`; `dirichlet-beta`/`xi` via their standard
defining formulas using mpmath's own independent gamma/zeta, not fastmath's).

Domains were chosen empirically (see investigation notes in the PR/session):
  - 1-arg `zeta` and `eta`: wide, safe across the whole double range after the
    log-space reflection-branch fix (saturating to +-Inf only where the true
    value itself exceeds double precision, which is not exercised here -
    the negative grid stays comfortably within double range).
  - 2-arg (Hurwitz) `zeta` and negative-domain `dirichlet-beta`: restricted to
    a conservatively safe range (s in [-3, 50]), since the Euler-Maclaurin
    direct-sum/asymptotic-tail branch selection in `zeta`'s 2-arity becomes
    unreliable for more negative `s` (a known, documented, NOT fixed in this
    session, limitation - see `zeta`'s docstring).
  - `xi`: wide (|s| up to ~400), safe after the log-space rewrite.

Run with (from repo root, using the `uv`-managed Python env mentioned in
AGENTS.md):
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_zeta_reference.py
"""
import mpmath as mp

mp.mp.dps = 40

OUT = "test/resources/special/zeta_reference.edn"


def to_double(x):
    """mpf -> python float, or None if it over/underflows double precision
    (kept out of the reference grids entirely: those are exercised via a
    handful of hand-verified explicit assertions in the test file instead,
    not the bulk data-driven comparison)."""
    d = float(x)
    if d != d or d in (float('inf'), float('-inf')):
        return None
    return d


def fmt(x):
    return repr(x)


def edn_vec(xs):
    return "[" + " ".join(fmt(x) for x in xs) + "]"


def dbeta(s):
    return mp.mpf(4) ** (-s) * (mp.zeta(s, mp.mpf('0.25')) - mp.zeta(s, mp.mpf('0.75')))


def xi(s):
    # mirrors fastmath's own reflection to avoid raw gamma poles at s<0
    if s < 0:
        s = 1 - s
    if s == 0 or s == 1:
        return mp.mpf('0.5')
    hs = s / 2
    return hs * (s - 1) * mp.pi ** (-hs) * mp.gamma(hs) * mp.zeta(s)


# ---- 1-arg zeta & eta: wide domain ----
z1_s = ([round(-0.25 * i, 4) for i in range(0, 4)] +      # 0 .. -0.75
        [-(1 + 0.37 * i) for i in range(0, 40)] +          # -1.37 .. -15.8ish, non-integer
        list(range(-1, -220, -1)) +                        # all integers -1..-219 (odd = nonzero, even = trivial zero)
        [0.001, 0.01, 0.1, 0.3, 0.7, 0.9, 0.999, 1.001, 1.01, 1.1,
         1.5, 2.0, 3.0, 5.0, 10.0, 50.0, 100.0, 1000.0, 100000.0])
z1_s_f, z1_ref, eta_ref = [], [], []
for s in z1_s:
    zv = to_double(mp.zeta(mp.mpf(s)))
    ev = to_double(mp.altzeta(mp.mpf(s)))
    if zv is not None and ev is not None:
        z1_s_f.append(s)
        z1_ref.append(zv)
        eta_ref.append(ev)
z1_s = z1_s_f

# ---- 2-arg (Hurwitz) zeta: restricted safe domain ----
z2_s = ([round(-0.3 * i, 4) for i in range(0, 11)] +       # 0 .. -3.0
        list(range(0, 51, 5)) + [0.5, 1.5, 2.5, 10.5, 30.5])
z2_z = [0.1, 0.3, 0.5, 0.9, 1.3, 2.0, 2.7, 5.1, 10.3]
z2_s_out, z2_z_out, z2_ref = [], [], []
for s in z2_s:
    for z in z2_z:
        v = to_double(mp.zeta(mp.mpf(s), mp.mpf(z)))
        if v is not None:
            z2_s_out.append(s)
            z2_z_out.append(z)
            z2_ref.append(v)

# ---- dirichlet-beta: wide positive, restricted negative ----
db_x_raw = ([round(-0.3 * i, 4) for i in range(0, 11) if round(-0.3 * i, 4) != -1] +  # avoid landing exactly on -1 (still fine, but keep grid non-integer-focused)
            [0.001, 0.01, 0.1, 0.3, 0.7, 0.9, 0.999, 1.001, 1.01, 1.1,
             1.5, 2.0, 3.0, 5.0, 10.0, 50.0, 100.0, 1000.0])
db_x, db_ref = [], []
for x in db_x_raw:
    v = to_double(dbeta(mp.mpf(x)))
    if v is not None:
        db_x.append(x)
        db_ref.append(v)

# ---- xi: wide domain (safe up to ~430 given the fix) ----
xi_s = ([round(-0.37 * i, 4) for i in range(0, 40)] +      # 0 .. ~-14.4, non-integer
        [round(0.37 * i, 4) for i in range(0, 40)] +       # 0 .. ~14.4
        list(range(-400, 401, 25)) +
        [0.5, -0.5, 1.5, -1.5, 100.3, -100.3, 300.7, -300.7])
xi_s_f, xi_ref = [], []
for s in xi_s:
    v = to_double(xi(mp.mpf(s)))
    if v is not None:
        xi_s_f.append(s)
        xi_ref.append(v)
xi_s = xi_s_f

with open(OUT, "w") as f:
    f.write("{:zeta1 {:arg " + edn_vec(z1_s) + "\n")
    f.write("         :ref " + edn_vec(z1_ref) + "}\n")
    f.write(" :eta {:arg " + edn_vec(z1_s) + "\n")
    f.write("       :ref " + edn_vec(eta_ref) + "}\n")
    f.write(" :zeta2 {:s " + edn_vec(z2_s_out) + "\n")
    f.write("         :z " + edn_vec(z2_z_out) + "\n")
    f.write("         :ref " + edn_vec(z2_ref) + "}\n")
    f.write(" :dbeta {:arg " + edn_vec(db_x) + "\n")
    f.write("         :ref " + edn_vec(db_ref) + "}\n")
    f.write(" :xi {:arg " + edn_vec(xi_s) + "\n")
    f.write("      :ref " + edn_vec(xi_ref) + "}}\n")

print("zeta1/eta points:", len(z1_s))
print("zeta2 points:", len(z2_s_out))
print("dbeta points:", len(db_x))
print("xi points:", len(xi_s))
print("wrote", OUT)
