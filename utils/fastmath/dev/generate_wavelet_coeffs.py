"""Generator of `resources/wavelets/coeffs.edn`, the wavelet filter table of
`fastmath.transform.wavelets`.

Every entry is `"name" [dec_lo dec_hi rec_lo rec_hi]` in PyWavelets orientation
(`pywt.Wavelet(name).filter_bank`). Sources:

- PyWavelets (exact copy of `filter_bank`): haar, db1-db38, sym2-sym20,
  coif1-coif5, bior*, rbio*, dmey (62 taps).
- Spectral factorisation in `mpmath` (120 digits) of the Daubechies polynomial
  P(y) = sum_k C(N-1+k, k) y^k, y = (2 - z - 1/z)/4:
  - db39-db45: all roots inside the unit circle (minimum phase); the same
    code reproduces PyWavelets db2-db38 exactly (checked on every run),
  - sym21-sym45: the root of every reciprocal pair is the one MATLAB chose,
    identified from the MATLAB table (roots of the MATLAB filter, ambiguous
    pairs resolved by the smallest coefficient distance; the best choice must
    beat the second best by a factor of at least 20). The same procedure
    reproduces PyWavelets sym2-sym20 to the precision of their table.
- MATLAB values copied by hand (`wavelet_coeffs_matlab.edn`), no precise source
  found: mb*, fk*, han*, bl7, bl9, bl10, beyl, vaid.
- Derived: battle3, the orthonormal cubic-spline Battle-Lemarie scaling filter
  (Mallat 1989) sampled from its frequency response and truncated to 161
  symmetric taps (n = -80..80), followed by one zero so the filter length is
  even (the quadrature mirror construction needs it); orthonormality error
  4.6e-13, multilevel reconstruction error 2.1e-11 (pywt, N 2048 and 601).

`la8`-`la20` are not stored: they equal sym4-sym10 (aliases in Clojure).

Run from anywhere:
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_wavelet_coeffs.py
"""
import itertools
import re
from pathlib import Path

import mpmath as mp
import numpy as np
import pywt

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "resources" / "wavelets" / "coeffs.edn"
MATLAB = Path(__file__).resolve().parent / "wavelet_coeffs_matlab.edn"

mp.mp.dps = 120


# ---------------------------------------------------------------- helpers

def qmf(filt):
    """PyWavelets `qmf`: reversed filter with every odd element negated."""
    q = list(filt[::-1])
    for i in range(1, len(q), 2):
        q[i] = -q[i]
    return q


def orthogonal_bank(rec_lo):
    """[dec_lo dec_hi rec_lo rec_hi] from a scaling filter, PyWavelets convention."""
    rec_lo = [float(x) for x in rec_lo]
    rec_hi = qmf(rec_lo)
    return [rec_lo[::-1], rec_hi[::-1], rec_lo, rec_hi]


def read_matlab():
    text = MATLAB.read_text()
    out = {}
    for name, body in re.findall(r'\["([^"]+)" \[((?:\[[^\[\]]*\]\s*)+)\]\]', text):
        dl, dh = [np.array([float(x) for x in v.split()]) for v in re.findall(r'\[([^\[\]]*)\]', body)]
        out[name] = (dl, dh)
    return out


# ------------------------------------------------- Daubechies factorisation

def y_roots(n):
    coeffs = [mp.binomial(n - 1 + k, k) for k in range(n)]
    return mp.polyroots(coeffs[::-1], maxsteps=500, extraprec=400) if n > 1 else []


def z_pair(y):
    """Reciprocal pair (inside, outside) of z with (2 - z - 1/z)/4 = y."""
    b = 2 - 4 * y
    d = mp.sqrt(b * b - 4)
    z1, z2 = (b + d) / 2, (b - d) / 2
    return (z1, z2) if abs(z1) < 1 else (z2, z1)


def root_groups(n):
    """Choice groups: a real y root gives one pair, a complex conjugate pair of y roots gives two."""
    ys = y_roots(n)
    used, groups = set(), []
    for i, y in enumerate(ys):
        if i in used:
            continue
        used.add(i)
        if abs(mp.im(y)) < mp.mpf(10) ** -60:
            groups.append([z_pair(y)])
        else:
            j = min((k for k in range(len(ys)) if k not in used), key=lambda k: abs(ys[k] - mp.conj(y)))
            used.add(j)
            groups.append([z_pair(y), z_pair(ys[j])])
    return groups


def scaling_filter(n, zs):
    """Coefficients (ascending powers) of ((1+x)/2)^n prod (1 - z x), normalised to sum sqrt(2)."""
    poly = [mp.mpc(1)]
    for _ in range(n):
        poly = np.convolve(poly, [mp.mpc(1), mp.mpc(1)]).tolist()
    for z in zs:
        poly = np.convolve(poly, [mp.mpc(1), -z]).tolist()
    h = [mp.re(c) for c in poly]
    s = mp.fsum(h)
    return [c * mp.sqrt(2) / s for c in h]


def db_scaling(n):
    return scaling_filter(n, [z_pair(y)[0] for y in y_roots(n)])


def float_filter(n, zs):
    p = np.array([1.0 + 0j])
    for _ in range(n):
        p = np.convolve(p, [1, 1])
    for z in zs:
        p = np.convolve(p, [1, -z])
    h = p.real
    return h * np.sqrt(2) / h.sum()


def sym_scaling(n, reference_rec_lo):
    """Symlet scaling filter whose root choice reproduces `reference_rec_lo`."""
    groups = root_groups(n)
    fgroups = [[(complex(a), complex(b)) for a, b in g] for g in groups]
    roots = np.roots(reference_rec_lo[::-1])
    inverse = 1 / roots[np.abs(roots + 1) > 0.05]   # drop the zeros at -1
    base, ambiguity = [], []
    for g in fgroups:
        d = [sum(np.min(np.abs(inverse - p[c])) for p in g) for c in (0, 1)]
        base.append(0 if d[0] < d[1] else 1)
        ambiguity.append(min(d) / max(d))
    open_groups = [i for i, a in enumerate(ambiguity) if a > 0.05]
    assert len(open_groups) <= 16, (n, len(open_groups))
    scored = []
    for bits in itertools.product([0, 1], repeat=len(open_groups)):
        choice = list(base)
        for i, b in zip(open_groups, bits):
            choice[i] = b
        h = float_filter(n, [p[c] for g, c in zip(fgroups, choice) for p in g])
        scored.append((float(np.max(np.abs(h - reference_rec_lo))), choice))
    scored.sort(key=lambda s: s[0])
    best = scored[0]
    # the runner-up must be clearly worse; sym45 has the smallest observed ratio (48)
    runner_up_ratio = scored[1][0] / best[0] if len(scored) > 1 else float("inf")
    assert runner_up_ratio > 20, (n, best[0], scored[1][0])
    zs = [p[c] for g, c in zip(groups, best[1]) for p in g]
    return scaling_filter(n, zs), best[0], runner_up_ratio


# ---------------------------------------------------------- Battle-Lemarie

def battle_lemarie(degree, taps, grid=1 << 14):
    w = 2 * np.pi * np.arange(grid) / grid

    def phi_hat(w):
        p = degree + 1
        x = w / 2
        sinc = np.where(x == 0, 1.0, np.sin(x) / np.where(x == 0, 1, x))
        k = np.arange(-2000, 2001)[:, None]
        t = w[None, :] / 2 + np.pi * k
        with np.errstate(divide='ignore', invalid='ignore'):
            terms = np.where(t == 0, 1.0, (np.sin(w[None, :] / 2) / np.where(t == 0, 1, t)) ** (2 * p))
        return sinc ** p / np.sqrt(terms.sum(0))

    h = np.real(np.fft.ifft(np.sqrt(2) * phi_hat(2 * w) / phi_hat(w)))
    half = taps // 2
    return np.concatenate([h[-half:], h[:half + 1]])


# ------------------------------------------------------------- the table

matlab = read_matlab()
table = {}
sources = {}

existing = (["haar"] + ["db%d" % n for n in range(1, 39)] + ["sym%d" % n for n in range(2, 21)]
            + ["coif%d" % n for n in range(1, 6)] + ["dmey"]
            + [n for n in pywt.wavelist(kind="discrete") if n.startswith(("bior", "rbio"))])
for name in existing:
    table[name] = [list(map(float, f)) for f in pywt.Wavelet(name).filter_bank]
    sources[name] = "pywt"

for n in range(2, 39):
    gen = orthogonal_bank(db_scaling(n))
    assert gen == table["db%d" % n], "db%d differs from pywt" % n

for n in range(2, 21):
    ref = np.array(pywt.Wavelet("sym%d" % n).rec_lo)
    h, dist, _ = sym_scaling(n, ref)
    assert dist < 1e-10, ("sym%d" % n, dist)

report = []
for n in range(39, 46):
    name = "db%d" % n
    table[name] = orthogonal_bank(db_scaling(n))
    sources[name] = "generated"
    report.append("%s: |generated - MATLAB| %.1e" % (name, np.max(np.abs(np.array(table[name][0]) - matlab[name][0]))))

for n in range(21, 46):
    name = "sym%d" % n
    h, dist, ratio = sym_scaling(n, matlab[name][0][::-1])
    table[name] = orthogonal_bank(h)
    sources[name] = "generated"
    report.append("%s: |generated - MATLAB| %.1e, runner-up %.0fx worse" % (name, dist, ratio))

for name, (dl, dh) in matlab.items():
    if name.startswith(("db", "sym")):
        continue
    assert np.array_equal(np.array(qmf(list(dl[::-1])))[::-1], dh), name   # MATLAB dec_hi follows the pywt convention
    table[name] = orthogonal_bank(list(dl[::-1]))
    sources[name] = "matlab"

table["battle3"] = orthogonal_bank(list(battle_lemarie(3, 161)) + [0.0])
sources["battle3"] = "derived"


def fmt(xs):
    return "[" + " ".join(repr(float(x)) for x in xs) + "]"


def sort_key(name):
    m = re.match(r"([a-z]+)(.*)", name)
    rest = m.group(2)
    nums = tuple(float(x) for x in re.findall(r"\d+", rest)) if rest else ()
    return (m.group(1), nums)


lines = [";; Generated by utils/fastmath/dev/generate_wavelet_coeffs.py (pywt %s). Do not edit." % pywt.__version__,
         ";; name -> [dec_lo dec_hi rec_lo rec_hi], PyWavelets orientation.",
         ";; Sources: pywt - copied from PyWavelets; generated - mpmath spectral factorisation;",
         ";; matlab - MATLAB values copied by hand; derived - computed from the definition.",
         "{"]
for name in sorted(table, key=sort_key):
    lines.append(' ;; %s' % sources[name])
    lines.append(' "%s" [%s]' % (name, " ".join(fmt(f) for f in table[name])))
lines.append("}")
OUT.write_text("\n".join(lines) + "\n")

print("\n".join(report))
print("entries:", len(table))
print("wrote", OUT)
