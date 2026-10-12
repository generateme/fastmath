"""Generator of reference values for `fastmath.signal-test`.

Produces `test/resources/signal/signal_reference.edn`:

- `:effects` - outputs of Python ports of the original effect sources for the
  4000-sample input `:x` (0.5 sin(2 pi 440 n/44100) + 0.3 sin(2 pi 3000 n/44100))
  and the unit impulse. Ports: RBJ Audio EQ Cookbook biquads with
  `scipy.signal.lfilter`, swh LADSPA (`svf_1214`, `decimator_1202`,
  `divider_1186`, `foverdrive_1196`), stmlib `Svf`, Common-DSP
  `SvfLinearTrapOptimised2` (Cytomic), mda-vst3 `mdaThruZeroProcessor`, CMT
  (`vcf303`, `lofi` bandwidth limiter), Audacity (`BassTreble` 2.4.2, `Echo`),
  omins `slew_limiter`. Each entry: `[effect-name params input-key output]`.
- `:spectra` - `scipy.signal.periodogram` (boxcar, detrend=False, fs=2.5) for
  N in {7, 8, 15, 16}, one- and two-sided, plus numpy magnitudes and phases.
- `:welch` - `scipy.signal.welch` (Hann window 16, noverlap 8, detrend=False)
  with 6 and 7 segments, mean and median averages, density and spectrum scaling.
- `:spectrogram` - `scipy.signal.spectrogram` for an odd window (9, noverlap 6).
- `:savgol` - `scipy.signal.savgol_filter(mode='constant')`, derivatives 0..3.
- `:convolve` - `scipy.signal.convolve` / `correlate` in all modes.
- `:hilbert` - `scipy.signal.hilbert`; `:resample` - `scipy.signal.resample`.

Run from anywhere:
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_signal_reference.py
"""
import math
from pathlib import Path

import numpy as np
import scipy.signal as ss

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "test" / "resources" / "signal" / "signal_reference.edn"

FS = 44100.0
N = 4000
n = np.arange(N)
X = 0.5 * np.sin(2 * np.pi * 440 * n / FS) + 0.3 * np.sin(2 * np.pi * 3000 * n / FS)
IMP = np.zeros(N)
IMP[0] = 1.0
INPUTS = {"x": X, "impulse": IMP}

# ---------------------------------------------------------------- biquads (RBJ)


def band_alpha(w0, bw):
    return np.sin(w0) * np.sinh(np.log(2) / 2 * bw * w0 / np.sin(w0))


def rbj(kind, fc, fs=FS, bw=1.0, gain=0.0, slope=1.5):
    w0 = 2 * np.pi * fc / fs
    c, s = np.cos(w0), np.sin(w0)
    A = 10 ** (gain / 40)
    if kind in ("lp", "hp", "bp", "eq"):
        al = band_alpha(w0, bw)
    if kind == "lp":
        b, a = [(1 - c) / 2, 1 - c, (1 - c) / 2], [1 + al, -2 * c, 1 - al]
    elif kind == "hp":
        b, a = [(1 + c) / 2, -(1 + c), (1 + c) / 2], [1 + al, -2 * c, 1 - al]
    elif kind == "bp":
        b, a = [al, 0.0, -al], [1 + al, -2 * c, 1 - al]
    elif kind == "eq":
        b, a = [1 + al * A, -2 * c, 1 - al * A], [1 + al / A, -2 * c, 1 - al / A]
    else:
        al = s / 2 * np.sqrt((A + 1 / A) * (1 / slope - 1) + 2)
        sa = 2 * np.sqrt(A) * al
        if kind == "ls":
            b = [A * ((A + 1) - (A - 1) * c + sa), 2 * A * ((A - 1) - (A + 1) * c), A * ((A + 1) - (A - 1) * c - sa)]
            a = [(A + 1) + (A - 1) * c + sa, -2 * ((A - 1) + (A + 1) * c), (A + 1) + (A - 1) * c - sa]
        else:
            b = [A * ((A + 1) + (A - 1) * c + sa), -2 * A * ((A - 1) + (A + 1) * c), A * ((A + 1) + (A - 1) * c - sa)]
            a = [(A + 1) - (A - 1) * c + sa, 2 * ((A - 1) - (A + 1) * c), (A + 1) - (A - 1) * c - sa]
    return np.array(b) / a[0], np.array(a) / a[0]


def biquad(x, *args, **kw):
    b, a = rbj(*args, **kw)
    return ss.lfilter(b, a, x)


def dj_eq(x, hi, mid, low, shelf_slope=1.5, peak_bw=1.0):
    y = biquad(x, "eq", 100.0, bw=peak_bw, gain=low)
    y = biquad(y, "eq", 1000.0, bw=peak_bw, gain=mid)
    return biquad(y, "hs", 10000.0, gain=hi, slope=shelf_slope)

# ---------------------------------------------------------------- simple filters


def simple_lowpass(x, rate=FS, cutoff=2000.0):
    dt = 1 / rate
    tau = 1 / (cutoff * 2 * np.pi)
    alpha = dt / (tau + dt)
    y = 0.0
    out = []
    for v in x:
        y = v * alpha + (y - y * alpha)
        out.append(y)
    return np.array(out)


def phaser(x, delay=0.5):
    a1 = (1 - delay) / (1 + delay)
    z = 0.0
    out = []
    for v in x:
        y = z - a1 * v
        z = v + a1 * y
        out.append(y)
    return np.array(out)

# ---------------------------------------------------------------- mda ThruZero (vst3)


def thru(x, speed=0.3, depth=0.43, mix=0.47, feedback=0.3, depth_mod=1.0, sr=FS):
    rat = 10 ** (3 * speed - 2) * 2 / sr
    if speed < 0.01:
        rat = 0.0
    dep = 2000 * depth * depth
    dem = dep - dep * depth_mod
    dep -= dem
    wet, dry, fb = mix, 1 - mix, 1.9 * feedback - 0.95
    buf = np.zeros(2048)
    ph, bp, f1 = 0.0, 0, 0.0
    out = []
    for a in x:
        ph += rat
        if ph > 1:
            ph -= 2
        bp = (bp - 1) & 0x7FF
        buf[bp] = a + fb * f1
        tmpf = dem + dep * (1 - ph * ph)
        tmp = int(tmpf)
        tmpf -= tmp
        tmp = (tmp + bp) & 0x7FF
        tmpi = (tmp + 1) & 0x7FF
        f1 = buf[tmp]
        f1 = tmpf * (buf[tmpi] - f1) + f1
        out.append(a * dry - f1 * wet)
    return np.array(out)

# ---------------------------------------------------------------- swh svf_1214


def svf1(x, fc=1000.0, q=0.5, res=0.0, t="l", R=3, sr=FS):
    f = 2 * math.sin(math.pi * fc / (sr * R))
    qq = 2 * math.cos(q ** 0.1 * math.pi * 0.5)
    qn = math.sqrt(qq / 2 + 0.01)
    h = b = l = 0.0
    out = []
    for s in x:
        inp = qn * (s + b * res)
        for _ in range(R):
            b = b - b * b * b * 0.001
            h = inp - l - qq * b
            b = b + f * h
            l = l + f * b
            o = {"l": l, "h": h, "b": b, "n": l + h, "p": l - h}[t]
            inp = o
        out.append(o)
    return np.array(out)

# ---------------------------------------------------------------- stmlib Svf


def tan_st(f, kind):
    P = math.pi
    if kind == "exact":
        return math.tan(P * min(f, 0.497))
    if kind == "dirty":
        return f * (P + 3.736e-01 * P ** 3 * f * f)
    f2 = f * f
    if kind == "fast":
        return f * (P + f2 * (3.260e-01 * P ** 3 + 1.823e-01 * P ** 5 * f2))
    a, b, c, d, e = (3.333314036e-01 * P ** 3, 1.333923995e-01 * P ** 5, 5.33740603e-02 * P ** 7,
                     2.900525e-03 * P ** 9, 9.5168091e-03 * P ** 11)
    return f * (P + f2 * (a + f2 * (b + f2 * (c + f2 * (d + f2 * e)))))


def svf2(x, fc=1000.0, Q=0.5, mode="lp", tan="exact", sr=FS):
    g = tan_st(fc / sr, tan)
    r = 1 / Q
    h = 1 / (1 + r * g + g * g)
    s1 = s2 = 0.0
    out = []
    for v in x:
        hp = (v - r * s1 - g * s1 - s2) * h
        bp = g * hp + s1
        s1 = g * hp + bp
        lp = g * bp + s2
        s2 = g * bp + lp
        out.append({"lp": lp, "bp": bp, "bpn": bp * r, "hp": hp, "notch": lp + hp}[mode])
    return np.array(out)

# ---------------------------------------------------------------- Cytomic / Common-DSP


def svf3(x, fc=1000.0, Q=0.5, typ="lp", gaindb=0.0, sr=FS):
    A = 10 ** (gaindb / 40)
    g = math.tan(math.pi * fc / sr)
    k = 1 / (Q * A) if typ == "bell" else 1 / Q
    if typ == "ls":
        g /= math.sqrt(A)
    if typ == "hs":
        g *= math.sqrt(A)
    a1 = 1 / (1 + g * (g + k))
    a2 = g * a1
    a3 = g * a2
    m = {"lp": (0, 0, 1), "bp": (0, 1, 0), "hp": (1, -k, -1), "notch": (1, -k, 0), "peak": (1, -k, -2),
         "ap": (1, -2 * k, 0), "bell": (1, k * (A * A - 1), 0), "ls": (1, k * (A - 1), A * A - 1),
         "hs": (A * A, k * (1 - A) * A, 1 - A * A)}[typ]
    ic1 = ic2 = 0.0
    out = []
    for v0 in x:
        v3 = v0 - ic2
        v1 = a1 * ic1 + a2 * v3
        v2 = ic2 + a2 * ic1 + a3 * v3
        ic1 = 2 * v1 - ic1
        ic2 = 2 * v2 - ic2
        out.append(m[0] * v0 + m[1] * v1 + m[2] * v2)
    return np.array(out)

# ---------------------------------------------------------------- swh decimator_1202


def decimator(x, bits=2.0, fsd=4410.0, sr=FS):
    if bits >= 31 or bits < 1:
        step, stepr = 0.0, 1.0
    else:
        step = 0.5 ** (bits - 0.999)
        stepr = 1 / step
    ratio = 1.0 if fsd >= sr else fsd / sr
    count, last = 0.0, 0.0
    out = []
    for v in x:
        count += ratio
        if count >= 1.0:
            count -= 1.0
            y = (v + (-1.0 if v < 0 else 1.0) * step * 0.5) * stepr
            last = v - math.fmod(y, 1.0) * step  # modf fractional part
        out.append(last)
    return np.array(out)

# ---------------------------------------------------------------- CMT vcf303 (trigger at start)


def vcf303(x, cutoff=0.8, reso=0.8, env=0.5, decay=1.0, trigger=True, sr=FS):
    e0 = math.exp(5.613 - 0.8 * env + 2.1553 * cutoff - 0.7696 * (1 - reso)) * math.pi / sr
    e1 = math.exp(6.109 + 1.5876 * env + 2.1553 * cutoff - 1.2 * (1 - reso)) * math.pi / sr
    c0 = e1 - e0 if trigger else 0.0
    d = 0.1 ** (1 / ((0.2 + 2.3 * decay) * sr))
    dec = d ** 64
    res = math.exp(-1.2 + 3.455 * reso)

    def abc(c0):
        w = e0 + c0
        k = math.exp(-w / res)
        a = 2 * math.cos(2 * w) * k
        b = -k * k
        return a, b, (1 - a - b) * 0.2
    a, b, c = abc(c0)
    d1 = d2 = 0.0
    ep = 0
    out = []
    for v in x:
        s = a * d1 + b * d2 + c * v
        out.append(s)
        d2, d1 = d1, s
        ep += 1
        if ep >= 64:
            ep = 0
            c0 *= dec
            a, b, c = abc(c0)
    return np.array(out)

# ---------------------------------------------------------------- swh divider_1186


def divider(x, den=2):
    out_, amp, count, lamp, last, z = 1.0, 0.0, 0.0, 0.0, 0.0, 0
    o = []
    for v in x:
        count += 1
        if (v > 0 and last <= 0) or (v < 0 and last >= 0):
            z += 1
            if den == 1:
                out_ = -1.0 if out_ > 0 else 1.0
                lamp = amp / count
                z, count, amp = 0, 0, 0
        amp += abs(v)
        if den > 1 and z % den == den - 1:
            out_ = -1.0 if out_ > 0 else 1.0
            lamp = amp / count
            z, count, amp = 0, 0, 0
        last = v
        o.append(out_ * lamp)
    return np.array(o)

# ---------------------------------------------------------------- Audacity BassTreble 2.4.2


def basstreble(x, bass=1.0, treble=1.0, gain=0.0, slope=0.4, bf=250.0, tf=4000.0, sr=FS):
    def coef(hz, g, typ):
        w = 2 * math.pi * hz / sr
        a = math.exp(math.log(10) * g / 40)
        b = math.sqrt((a * a + 1) / slope - (a - 1) ** 2)
        c, s = math.cos(w), math.sin(w)
        if typ == "b":
            return (a * ((a + 1) - (a - 1) * c + b * s), 2 * a * ((a - 1) - (a + 1) * c),
                    a * ((a + 1) - (a - 1) * c - b * s), (a + 1) + (a - 1) * c + b * s,
                    -2 * ((a - 1) + (a + 1) * c), (a + 1) + (a - 1) * c - b * s)
        return (a * ((a + 1) + (a - 1) * c + b * s), -2 * a * ((a - 1) + (a + 1) * c),
                a * ((a + 1) + (a - 1) * c - b * s), (a + 1) - (a - 1) * c + b * s,
                2 * ((a - 1) - (a + 1) * c), (a + 1) - (a - 1) * c - b * s)
    B, T, G = coef(bf, bass, "b"), coef(tf, treble, "t"), 10 ** (gain / 20)
    yb = ss.lfilter(B[:3], B[3:], x)
    return ss.lfilter(T[:3], T[3:], yb) * G

# ---------------------------------------------------------------- others


def bandwidth_limit(x, freq=1000.0, sr=FS):
    dx = freq / sr
    s = 0.0
    o = []
    for v in x:
        s = min(s + dx, v) if v >= s else max(s - dx, v)
        o.append(s)
    return np.array(o)


def echo(x, delay=0.5, decay=0.5, sr=FS):
    hl = int(sr * delay)
    h = np.zeros(hl)
    p = 0
    o = []
    for v in x:
        if p == hl:
            p = 0
        h[p] = v + h[p] * decay
        o.append(h[p])
        p += 1
    return np.array(o)


def slew(x, maxrise=500.0, maxfall=500.0, sr=FS):
    y = 0.0
    o = []
    for v in x:
        y = y + min(max(v - y, -maxfall / sr), maxrise / sr)
        o.append(y)
    return np.array(o)


EFFECTS = [
    ["simple-lowpass", {}, "x", simple_lowpass(X)],
    ["simple-highpass", {"cutoff": 500.0}, "x", X - simple_lowpass(X, cutoff=500.0)],
    ["biquad-lp", {"fc": 1500.0, "bw": 0.5}, "x", biquad(X, "lp", 1500.0, bw=0.5)],
    ["biquad-lp", {}, "impulse", biquad(IMP, "lp", 1000.0)],
    ["biquad-hp", {"fc": 1500.0}, "x", biquad(X, "hp", 1500.0)],
    ["biquad-bp", {"fc": 3000.0, "bw": 2.0}, "x", biquad(X, "bp", 3000.0, bw=2.0)],
    ["biquad-eq", {"fc": 3000.0, "gain": 6.0, "bw": 0.5}, "x", biquad(X, "eq", 3000.0, gain=6.0, bw=0.5)],
    ["biquad-hs", {"fc": 2000.0, "gain": -6.0}, "x", biquad(X, "hs", 2000.0, gain=-6.0)],
    ["biquad-ls", {"fc": 500.0, "gain": 6.0, "slope": 1.0}, "x", biquad(X, "ls", 500.0, gain=6.0, slope=1.0)],
    ["dj-eq", {"hi": 6.0, "mid": -3.0, "low": 4.0}, "x", dj_eq(X, 6.0, -3.0, 4.0)],
    ["phaser-allpass", {"delay": 0.3}, "x", phaser(X, 0.3)],
    ["divider", {}, "x", divider(X, 2)],
    ["divider", {"denom": 1}, "x", divider(X, 1)],
    ["bandwidth-limit", {}, "x", bandwidth_limit(X)],
    ["distort", {"factor": 2.0}, "x", 3.0 * X / (2.0 + np.abs(X))],
    ["foverdrive", {"drive": 2.5}, "x", X * (np.abs(X) + 2.5) / (X * X + 1.5 * np.abs(X) + 1)],
    ["decimator", {"bits": 4.0}, "x", decimator(X, bits=4.0)],
    ["decimator", {}, "x", decimator(X)],
    ["decimator", {"bits": 0.5, "fs": 50000.0}, "x", decimator(X, bits=0.5, fsd=50000.0)],
    ["basstreble", {}, "x", basstreble(X)],
    ["basstreble", {"bass": 6.0, "treble": -6.0, "gain": 2.0}, "x", basstreble(X, 6.0, -6.0, 2.0)],
    ["echo", {"delay": 0.01}, "x", echo(X, 0.01)],
    ["vcf303", {"trigger": True}, "x", vcf303(X)],
    ["vcf303", {"cutoff": 0.3, "resonance": 0.5}, "x", vcf303(X, cutoff=0.3, reso=0.5, trigger=False)],
    ["slew-limit", {"maxrise": 2000.0, "maxfall": 800.0}, "x", slew(X, 2000.0, 800.0)],
    ["mda-thru-zero", {}, "x", thru(X)],
    ["mda-thru-zero", {"speed": 0.6, "depth": 0.7, "mix": 0.3, "feedback": 0.8, "depth-mod": 0.5}, "x",
     thru(X, speed=0.6, depth=0.7, mix=0.3, feedback=0.8, depth_mod=0.5)],
    ["mda-thru-zero", {"speed": 0.005}, "x", thru(X, speed=0.005)],
    ["gain", {"gain": 0.5}, "x", 0.5 * X],
    ["clipping", {"pregain": 3.0}, "x", np.clip(3 * X, -1, 1)],
    ["clipping", {"pregain": 3.0, "method": "hyperbolic"}, "x", np.tanh(3 * X)],
    ["clipping", {"pregain": 3.0, "method": "soft"}, "x",
     np.where(3 * X < -1, -2 / 3, np.where(3 * X > 1, 2 / 3, 3 * X - (3 * X) ** 3 / 3))],
    ["svf1", {"oversamples": 3, "resonance": 0.0}, "x", svf1(X)],
    ["svf1", {"oversamples": 3, "resonance": 0.0, "kind": "allpass"}, "x", svf1(X, t="p")],
    ["svf1", {"Q": 0.9, "cutoff": 3000.0, "kind": "bandpass", "resonance": 0.2}, "x",
     svf1(X, fc=3000.0, q=0.9, res=0.2, t="b", R=1)],
    ["svf1", {"kind": "highpass", "oversamples": 2}, "x", svf1(X, res=0.01, t="h", R=2)],
    ["svf1", {"kind": "notch"}, "x", svf1(X, res=0.01, t="n", R=1)],
]
for kind, mode in [("lowpass", "lp"), ("highpass", "hp"), ("bandpass", "bp"),
                   ("bandpass-normalized", "bpn"), ("notch", "notch")]:
    EFFECTS.append(["svf2", {"Q": 2.0, "kind": kind}, "x", svf2(X, Q=2.0, mode=mode)])
for tan in ["dirty", "fast", "accurate"]:
    EFFECTS.append(["svf2", {"Q": 2.0, "cutoff": 5000.0, "tan": tan}, "impulse",
                    svf2(IMP, fc=5000.0, Q=2.0, tan=tan)])
for kind, typ in [("lowpass", "lp"), ("bandpass", "bp"), ("highpass", "hp"), ("notch", "notch"), ("peak", "peak"),
                  ("allpass", "ap"), ("bell", "bell"), ("lowshelf", "ls"), ("highshelf", "hs")]:
    EFFECTS.append(["svf3", {"Q": 2.0, "kind": kind, "gaindb": 6.0, "cutoff": 2000.0}, "x",
                    svf3(X, fc=2000.0, Q=2.0, typ=typ, gaindb=6.0)])

# ---------------------------------------------------------------- spectra

rng = np.random.default_rng(20261012)
SPECTRA = []
for size in (7, 8, 15, 16):
    x = rng.normal(size=size)
    for onesided in (True, False):
        f, psd = ss.periodogram(x, fs=2.5, window="boxcar", detrend=False, return_onesided=onesided,
                                scaling="density")
        _, power = ss.periodogram(x, fs=2.5, window="boxcar", detrend=False, return_onesided=onesided,
                                  scaling="spectrum")
        X_ = np.fft.rfft(x) if onesided else np.fft.fft(x)
        w = np.ones(len(X_))
        if onesided:
            w[1:] = 2.0
            if size % 2 == 0:
                w[-1] = 1.0
        mag = w * np.abs(X_)
        SPECTRA.append({"x": x, "one-sided?": onesided, "freqs": np.arange(len(X_)) * 2.5 / size,
                        "psd": psd, "power": power, "energy": power * size * size,
                        "magnitude": mag, "amplitude": mag / size, "asd": np.sqrt(psd),
                        "phase": np.angle(X_), "real": X_.real, "imag": X_.imag})

WIN16 = ss.get_window("hann", 16)
WELCH = []
for size in (56, 64):  # 6 and 7 segments
    x = rng.normal(size=size)
    entry = {"x": x, "window": WIN16}
    for avg in ("mean", "median"):
        for scaling in ("density", "spectrum"):
            f, p = ss.welch(x, fs=3.0, window=WIN16, noverlap=8, detrend=False, scaling=scaling, average=avg)
            entry[avg + "-" + scaling] = p
    entry["freqs"] = f
    WELCH.append(entry)

WIN9 = ss.get_window("hamming", 9)
xs = rng.normal(size=40)
SPEC = {"x": xs, "window": WIN9}
for scaling, key in (("density", "psd"), ("spectrum", "power")):
    f, t, S = ss.spectrogram(xs, fs=2.0, window=WIN9, noverlap=6, detrend=False, scaling=scaling, mode="psd")
    SPEC[key] = S.T
SPEC["freqs"], SPEC["times"] = f, t

# ---------------------------------------------------------------- savgol, convolution, hilbert, resample

xsg = rng.normal(size=30)
SAVGOL = {"x": xsg, "default": ss.savgol_filter(xsg, 5, 2, mode="constant")}
for d in range(4):
    SAVGOL["deriv-%d" % d] = ss.savgol_filter(xsg, 7, 3, deriv=d, delta=0.5, mode="constant")

a7 = rng.normal(size=7)
b4 = rng.normal(size=4)
CONV = {"a": a7, "b": b4}
for mode in ("full", "same", "valid"):
    CONV["convolve-ab-" + mode] = ss.convolve(a7, b4, mode=mode, method="direct")
    CONV["convolve-ba-" + mode] = ss.convolve(b4, a7, mode=mode, method="direct")
    CONV["correlate-ab-" + mode] = ss.correlate(a7, b4, mode=mode, method="direct")
    CONV["correlate-ba-" + mode] = ss.correlate(b4, a7, mode=mode, method="direct")

HILBERT = []
for size in (9, 10):
    x = rng.normal(size=size)
    h = ss.hilbert(x)
    HILBERT.append({"x": x, "real": h.real, "imag": h.imag})

RESAMPLE = []
for size, targets in ((10, (1, 4, 7, 15)), (9, (2, 5, 12))):
    x = rng.normal(size=size)
    for t in targets:
        RESAMPLE.append({"x": x, "n": t, "y": ss.resample(x, t)})

# ---------------------------------------------------------------- EDN


def edn(v):
    if isinstance(v, (bool, np.bool_)):
        return "true" if v else "false"
    if isinstance(v, str):
        return ":" + v
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        return repr(float(v))
    if isinstance(v, dict):
        return "{" + " ".join(":%s %s" % (k, edn(x)) for k, x in v.items()) + "}"
    if isinstance(v, np.ndarray) and v.ndim == 2:
        return "[" + "\n  ".join(edn(r) for r in v) + "]"
    return "[" + " ".join(edn(x) for x in v) + "]"


def effect_edn(e):
    name, params, inp, out = e
    return "[:%s %s :%s %s]" % (name, edn(params), inp, edn(np.asarray(out, float)))


OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w") as fo:
    fo.write(";; generated by utils/fastmath/dev/generate_signal_reference.py\n")
    fo.write("{:x " + edn(X) + "\n")
    fo.write(" :effects [" + "\n  ".join(effect_edn(e) for e in EFFECTS) + "]\n")
    fo.write(" :spectra " + edn(SPECTRA) + "\n")
    fo.write(" :welch " + edn(WELCH) + "\n")
    fo.write(" :spectrogram " + edn(SPEC) + "\n")
    fo.write(" :savgol " + edn(SAVGOL) + "\n")
    fo.write(" :convolve " + edn(CONV) + "\n")
    fo.write(" :hilbert " + edn(HILBERT) + "\n")
    fo.write(" :resample " + edn(RESAMPLE) + "}\n")

print("wrote", OUT)
