"""Generator of reference data for the flat-top windows of `fastmath.kernel.window` (kernel audit, WP 5.4).

Produces `test/resources/kernel/flattop_psl_reference.edn`.

Tools: Python 3, mpmath 1.3 (50 digits), numpy 2.x (float64 FFT).

Why: the flat-top sets `flat-top-hp`, `flat-top-srs`, `flat-top-bk`, `hft70`, `hft90d`, `hft95`, `hft116d`,
`hft144d`, `hft169d`, `hft196d`, `hft223d`, `hft248d` are sums of cosines.  Their coefficients were copied from
Heinzel, Ruediger and Schilling, Spectrum and spectral density estimation by the Discrete Fourier Transform
(DFT), including a comprehensive list of window functions and some new flat-top windows, 2002, section D.2,
D.3 (flat-top-bk: Gade and Herlufsen 1987 = Octave `flattopwin`).  A wrong digit changes the peak sidelobe
level (PSL), so the test compares the PSL of the fastmath window with the PSL of the exact transfer function.

For a window `w(x) = sum_k a_k cos(2 pi k x)` on `[-1/2, 1/2]` (`a_k` the positive magnitudes of the
fastmath coefficient list) the transfer function in units of the bin `1/N` is (`sinc(f) = sin(pi f) / (pi f)`)

    W(f) = a_0 sinc(f) + sum_{k >= 1} a_k / 2 (sinc(f - k) + sinc(f + k))
         = sin(pi f) / pi * (a_0 / f + sum_{k >= 1} (-1)^k a_k f / (f^2 - k^2)),     W(0) = a_0.

Every case records
  :coeffs       the positive magnitudes
  :end-value    sum (-1)^k a_k / sum a_k, the value of the default normalised window at both ends
  :nenbw        a0^-2 (a_0^2 + sum_{k>=1} a_k^2 / 2) in bins - sum a_k, closed form, exact
  :first-zero   first zero of W(f) (mpmath, bisection), bins
  :psl-db       max_{f >= first zero, f <= 100} 20 log10 |W(f) / W(0)| (mpmath, scan step 0.01 + golden section)
  :peak-f       position of this sidelobe, bins
  :dft-psl-db   the same measured as the fastmath test does: periodic window N = 8192 in float64, FFT with
                zero padding to 2^18, max over f in [first zero, 100] (shows where double precision limits)
  :heinzel      Table 2 / section D values {:psl-db :first-zero :peak-f :nenbw} (nil where Heinzel has none)
  :printed      only hft169d: the coefficients exactly as printed in the report (see the comment at SETS)

The sidelobes of the equiripple designs (hft90d, hft116d, hft144d) have nearly equal heights, so the position
of the highest one in the report (:peak-f) can differ from the position found here; the level agrees.

Run from anywhere (output is byte-identical on every run):
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_flattop_psl.py
"""
from pathlib import Path

import mpmath as mp
import numpy as np

mp.mp.dps = 50

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "test" / "resources" / "kernel" / "flattop_psl_reference.edn"

N_DFT = 8192
PAD = 1 << 18
F_MAX = 100

# name -> (coefficient strings (positive magnitudes), Heinzel values or None)
# Heinzel: PSLL from Table 2 (positive dB in the paper, negative here), first zero, highest sidelobe location,
# NENBW from sections D.2, D.3 (text; first zero of FTNI 3.00 / FTSRS 4.72 are quoted, FTHP 4.00)
SETS = [
    ("flat-top-hp", ["1", "1.912510941", "1.079173272", "0.1832630879"],
     dict(psl=-70.4, zero=4.0, peak=4.65, nenbw=3.4279)),
    ("flat-top-srs", ["1", "1.93", "1.29", "0.388", "0.028"],
     dict(psl=-76.6, zero=4.72, peak=5.37, nenbw=3.7702)),
    ("flat-top-bk", ["1", "1.93", "1.29", "0.388", "0.0322"], None),
    ("hft70", ["1", "1.90796", "1.07349", "0.18199"],
     dict(psl=-70.4, zero=4.0, peak=4.65, nenbw=3.4129)),
    ("hft90d", ["1", "1.942604", "1.340318", "0.440811", "0.043097"],
     dict(psl=-90.2, zero=5.0, peak=5.58, nenbw=3.8832)),
    ("hft95", ["1", "1.9383379", "1.3045202", "0.4028270", "0.0350665"],
     dict(psl=-95.0, zero=5.0, peak=7.49, nenbw=3.8112)),
    ("hft116d", ["1", "1.9575375", "1.4780705", "0.6367431", "0.1228389", "0.0066288"],
     dict(psl=-116.8, zero=6.0, peak=7.52, nenbw=4.2186)),
    ("hft144d", ["1", "1.96760033", "1.57983607", "0.81123644", "0.22583558", "0.02773848", "0.00090360"],
     dict(psl=-144.1, zero=7.0, peak=7.07, nenbw=4.5386)),
    # a_4 is 0.33673418, not the printed 0.33673420 (Heinzel eq. 182): the printed digits give
    # 1 + sum c_k = 2.0e-8 (not 0, condition (161) of the report), PSL -167.88 dB at 10.42 bins; with a_4 lowered
    # by 2e-8 the end value is 0 and the PSL is -169.52 dB at 10.413 bins (report: -169.5 dB at 10.41 bins).
    # Lowering any other single coefficient gives -167.4 ... -169.45 dB (see PRINTED below).
    ("hft169d", ["1", "1.97441842", "1.65409888", "0.95788186", "0.33673418", "0.06364621", "0.00521942",
                 "0.00010599"],
     dict(psl=-169.5, zero=8.0, peak=10.41, nenbw=4.8347)),
    ("hft196d", ["1", "1.979280420", "1.710288951", "1.081629853", "0.448734314", "0.112376628", "0.015122992",
                 "0.000871252", "0.000011896"],
     dict(psl=-196.2, zero=9.0, peak=9.06, nenbw=5.1134)),
    ("hft223d", ["1", "1.98298997309", "1.75556083063", "1.19037717712", "0.56155440797", "0.17296769663",
                 "0.03233247087", "0.00324954578", "0.00013801040", "0.00000132725"],
     dict(psl=-223.0, zero=10.0, peak=11.38, nenbw=5.3888)),
    ("hft248d", ["1", "1.985844164102", "1.791176438506", "1.282075284005", "0.667777530266", "0.240160796576",
                 "0.056656381764", "0.008134974479", "0.000624544650", "0.000019808998", "0.000000132974"],
     dict(psl=-248.4, zero=11.0, peak=13.37, nenbw=5.6512)),
]


# The coefficients of HFT169D as printed in the report (end value 2.0e-8, not 0)
PRINTED = {"hft169d": ["1", "1.97441842", "1.65409888", "0.95788186", "0.33673420", "0.06364621", "0.00521942",
                       "0.00010599"]}


def transfer(a, f):
    """W(f) for the positive magnitudes `a` (mpf list), f an mpf."""
    s = a[0] / f
    for k in range(1, len(a)):
        s += (-1) ** k * a[k] * f / (f * f - k * k)
    return mp.sin(mp.pi * f) / mp.pi * s


def first_zero(a):
    h = mp.mpf("0.01")
    f = mp.mpf("0.005")
    prev = transfer(a, f)
    while True:
        g = f + h
        cur = transfer(a, g)
        if prev * cur < 0:
            lo, hi = f, g
            for _ in range(170):
                mid = (lo + hi) / 2
                if transfer(a, mid) * prev > 0:
                    lo = mid
                else:
                    hi = mid
            return (lo + hi) / 2
        f, prev = g, cur


def peak_sidelobe(a, f0):
    """Max |W(f)| / a0 for f in [f0, F_MAX]: scan step 0.01, golden section around the best grid point."""
    h = mp.mpf("0.01")
    n = int((F_MAX - f0) / h)
    vals = [(abs(transfer(a, f0 + i * h)), f0 + i * h) for i in range(n + 1)]
    best = max(vals)
    f = best[1]
    lo, hi = f - h, f + h
    gr = (mp.sqrt(5) - 1) / 2
    c = hi - gr * (hi - lo)
    d = lo + gr * (hi - lo)
    fc, fd = abs(transfer(a, c)), abs(transfer(a, d))
    for _ in range(120):
        if fc > fd:
            hi, d, fd = d, c, fc
            c = hi - gr * (hi - lo)
            fc = abs(transfer(a, c))
        else:
            lo, c, fc = c, d, fd
            d = lo + gr * (hi - lo)
            fd = abs(transfer(a, d))
    fm = (lo + hi) / 2
    return abs(transfer(a, fm)) / a[0], fm


def dft_psl(a, f0):
    """float64 periodic window of N_DFT samples, FFT with zero padding to PAD, max |W| / |W(0)| for f in [f0, F_MAX]."""
    af = [float(c) for c in a]
    j = np.arange(N_DFT)
    x = j / N_DFT - 0.5
    w = np.zeros(N_DFT)
    for k, c in enumerate(af):
        w += c * np.cos(2.0 * np.pi * k * x)
    spec = np.abs(np.fft.rfft(w, PAD))
    f = np.arange(spec.size) * (N_DFT / PAD)
    sel = (f >= float(f0)) & (f <= F_MAX)
    return float(20 * np.log10(spec[sel].max() / spec[0]))


def fl(x):
    return repr(float(x))


def vec(xs):
    return "[" + " ".join(fl(x) for x in xs) + "]"


rows = []
for name, cs, hz in SETS:
    a = [mp.mpf(c) for c in cs]
    f0 = first_zero(a)
    pk, fpk = peak_sidelobe(a, f0)
    psl = 20 * mp.log10(pk)
    nenbw = (a[0] ** 2 + sum(c * c for c in a[1:]) / 2) / a[0] ** 2
    ssum = sum(a)
    end = sum((-1) ** k * c for k, c in enumerate(a)) / ssum
    dpsl = dft_psl(a, f0)
    heinzel = "nil" if hz is None else "{:psl-db %s :first-zero %s :peak-f %s :nenbw %s}" % (
        fl(hz["psl"]), fl(hz["zero"]), fl(hz["peak"]), fl(hz["nenbw"]))
    printed = ""
    if name in PRINTED:
        pa = [mp.mpf(c) for c in PRINTED[name]]
        pf0 = first_zero(pa)
        ppk, _ = peak_sidelobe(pa, pf0)
        pend = sum((-1) ** k * c for k, c in enumerate(pa)) / sum(pa)
        printed = ' :printed {:coeffs %s :end-value %s :psl-db %s}' % (vec(pa), fl(pend), fl(20 * mp.log10(ppk)))
    rows.append(
        '{:window :%s :coeffs %s :end-value %s :nenbw %s :first-zero %s :psl-db %s :peak-f %s :dft-psl-db %s '
        ':heinzel %s%s}' % (name, vec(a), fl(end), fl(nenbw), fl(f0), fl(psl), fl(fpk), fl(dpsl), heinzel, printed))
    print("%-13s zero %8.5f  psl %9.3f at %8.4f  dft %9.3f  nenbw %.5f  end %.3e  (heinzel %s)" % (
        name, float(f0), float(psl), float(fpk), dpsl, float(nenbw), float(end), None if hz is None else hz["psl"]))

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w") as fh:
    fh.write(";; Generated by utils/fastmath/dev/generate_flattop_psl.py (mpmath %s, numpy %s).\n"
             % (mp.__version__, np.__version__))
    fh.write(";; N-DFT %d, zero padding %d, sidelobe search range [first zero, %d] bins.\n" % (N_DFT, PAD, F_MAX))
    fh.write("{:n-dft %d :pad %d :f-max %d\n :cases [\n" % (N_DFT, PAD, F_MAX))
    for r in rows:
        fh.write(" " + r + "\n")
    fh.write("]}\n")
print("wrote", OUT)
