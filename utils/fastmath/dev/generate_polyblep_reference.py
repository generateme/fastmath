"""Generator of reference values for the `polyblep-*` waveforms of `fastmath.signal.waveform`.

A literal transliteration of Tale's `poly_blep.jsfx-inc` (https://github.com/TaleTN/JSFX, WTFPL):
every waveform function is copied from the JSFX source, with the polynomials expanded inline
where the source does so. `t - (t|0)` is `frac1`, which asserts a non-negative argument
(`|0` truncates, the Clojure port uses `floor`). No JSFX interpreter is involved.

Differences from the JSFX, as in the Clojure port:
- `pw` is clamped to [0, 1] and, where a waveform divides by it (`half2`, `trip`, `hwsaw`), to >= 1e-4,
- `abs(dt)` under `sqrt` in `arc` (the term is 0 for a negative `dt`).

Produces `test/resources/signal/polyblep_reference.edn`: for every waveform the samples of the
generator with fs = 1000, phase = 0.37, amplitude = 1.5, 80 samples, for the phase steps
dt in {0.013, 0.1, 0.23} (f = dt * fs) and for the extra parameter:
pw in {0.25, 0.7}, n in {1, 2, 5}, fc = f2 / f in {0.6, 1.5, 3.4}.
`:time-varying` has three waveforms with f(t) = 100 + 80 t, pw(t) = 0.3 + 0.4 t, f2 = 250.

Run from anywhere:
    cd /home/ts/penv && uv run python \
        /home/ts/clojure/fastmath/utils/fastmath/dev/generate_polyblep_reference.py
"""
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "test" / "resources" / "signal" / "polyblep_reference.edn"

PI = math.pi
TWO_PI = 6.2831853071795865


def sqr(x):
    return x * x


def frac1(x):
    """`x - (x|0)`"""
    assert x >= 0, x
    return x - int(x)


def blep(t, dt):
    if t < dt:
        return -sqr(t / dt - 1)
    if t > 1 - dt:
        return sqr((t - 1) / dt + 1)
    return 0.0


def blamp(t, dt):
    if t < dt:
        t = t / dt - 1
        return sqr(t) * t * -1 / 3
    if t > 1 - dt:
        t = (t - 1) / dt + 1
        return sqr(t) * t * 1 / 3
    return 0.0


def bluh(t, dt):
    if t < dt:
        t = sqr(t / dt - 1)
        return (sqr(t) - t) * -4
    if t > 1 - dt:
        t = sqr((t - 1) / dt + 1)
        return (sqr(t) - t) * 4
    return 0.0


def gain(dt):
    """`_poly_setdt`"""
    if dt <= 0.2:
        return 1.0
    if dt < 0.25:
        return 1 - sqr((dt - 0.2) * 20)
    return 0.0


# ---- waveforms: (this.t, dt, p) -> sample without the gain `a`

def sinn(t1, dt, n):
    n = int(n)
    if n >= 2:
        y1 = math.sin(int(t1 * n * 2 + 0.5) * PI / n)
        dt2 = 0.5 / n
        t1 += dt2 * 0.5
    else:
        y1 = 0.70710678118654752 if t1 < 0.5 else -0.70710678118654752
        dt2 = 0.25
    t2 = t1 + 0.5
    dt3 = dt2 * TWO_PI
    t3 = dt3
    y2 = 0.0
    s0 = 0.0
    for _ in range(n):
        t1 = frac1(t1)
        t2 = frac1(t2)
        s1 = s0
        s0 = math.sin(t3)
        y2 += (blep(t1, dt) - blep(t2, dt)) * (s0 - s1)
        t1 += dt2
        t2 += dt2
        t3 += dt3
    return (0.5 if n >= 2 else 0.70710678118654752) * y2 + y1


def half(t, dt, _):
    return ((blamp(t, dt) + blamp(t + 0.5 - int(t + 0.5), dt)) * TWO_PI * dt
            + (math.sin(t * TWO_PI) * 2 if t < 0.5 else 0) - 0.63661977236758134)


def half2(this_t, dt, pw):
    pw = max(pw, 1e-4)
    t1 = pw * 0.5 + this_t + 0.75
    t2 = t1 - pw
    t1 = frac1(t1)
    return ((blamp(t2 - int(t2), dt) + blamp(t1, dt)) / pw * PI * dt
            + (math.sin(t1 / pw * PI) * 2 if t1 < pw else 0) - pw * 1.2732395447351627)


def full(this_t, dt, _):
    t = this_t + 0.25
    t = frac1(t)
    return (blamp(t, dt) * TWO_PI * dt
            + math.sin(t * PI) * 2 - 1.2732395447351627)


def sinp(t, dt, _):
    return (-blep(t + 0.5 - int(t + 0.5), dt)
            + blamp(t, dt) * PI * dt
            + (math.sin(t * PI) * 2 if t < 0.5 else 0) - 0.63661977236758134)


def sawp(this_t, dt, _):
    t = this_t + 0.5
    t = frac1(t)
    return (-blep(t, dt)
            + (math.cos(t * TWO_PI) + 1 if t >= 0.5 else 0) - 0.5)


def alt(t, dt, _):
    return ((blamp(t + 0.5 - int(t + 0.5), dt) - blamp(t, dt)) * -TWO_PI * dt
            + (math.sin(t * 12.566370614359173) if t < 0.5 else 0))


def camel(t1, dt, _):
    t2 = t1 + 0.5
    t3 = t2 + 0.25
    return ((blamp(t3 - int(t3), dt) * 2
             + blamp(t2 - int(t2), dt)
             + blamp(t1, dt)) * 12.566370614359173 * dt
            + (abs(math.sin(t1 * 12.566370614359173)) * 2 if t1 < 0.5 else 0) - 0.63661977236758134)


def camela(t1, dt, _):
    t2 = t1 + 0.5
    t3 = t2 + 0.25
    y = math.sin(t1 * 12.566370614359173)
    if t1 < 0.5:
        y = abs(y)
    return ((blamp(t2 - int(t2), dt)
             + blamp(t3 - int(t3), dt)) * 12.566370614359173 * dt
            + y - 0.31830988618379067)


def camel2(this_t, dt, _):
    t1 = this_t + 0.25
    t2 = t1 + 0.5
    t1 = frac1(t1)
    return ((blamp(t2 - int(t2), dt)
             - blamp(t1, dt)) * 12.566370614359173 * dt
            + math.sin(abs(t1 - 0.5) * 12.566370614359173))


def sin2(t, dt, _):
    return ((bluh(t + 0.5 - int(t + 0.5), dt)
             - bluh(t, dt)) * sqr(dt) * -1.7724538509055160
            + (1 if t < 0.5 else -1) * sqr(math.sin(t * TWO_PI)))


def para(t, dt, _):
    return ((bluh(t + 0.5 - int(t + 0.5), dt)
             - bluh(t, dt)) * 0.66666666666666667 * sqr(dt)
            + (t * 16 - 8) * (abs(t - 0.5) - 0.5))


def arc(this_t, dt, _):
    t = this_t + 0.25
    t = frac1(t)
    return (blamp(t, dt) * math.sqrt(abs(dt)) * TWO_PI
            + math.sqrt(t - sqr(t)) * 4 - 1)


def tri(t0, dt, _):
    t1 = t0 + 0.25
    t2 = t1 + 0.5
    y = abs((int(t0 + 0.75) - t0) * 4 - 1) - 1
    t1 -= int(t1)
    if t1 < dt:
        p = t1 / dt - 1
        y -= sqr(p) * p * 1.3333333333333333 * dt
    elif t1 > 1 - dt:
        p = (t1 - 1) / dt + 1
        y += sqr(p) * p * 1.3333333333333333 * dt
    t2 -= int(t2)
    if t2 < dt:
        p = t2 / dt - 1
        y += sqr(p) * p * 1.3333333333333333 * dt
    elif t2 > 1 - dt:
        p = (t2 - 1) / dt + 1
        y -= sqr(p) * p * 1.3333333333333333 * dt
    return y


def hyptri(t1, dt, _):
    t2 = t1 + 0.5
    return ((blamp(t1, dt) * 1.2732395447351627
             - blamp(t2 - int(t2), dt) * 9.273) * dt
            + math.exp((t1 if t1 < 0.5 else 1 - t1) * 4) * 0.31303528549933130 - 1)


def trin(this_t, dt, n):
    n = int(n)
    dt2 = 0.5 / n
    t1 = dt2 * 0.5 + 0.25 + this_t
    t1 -= int(t1)
    y1 = int(t1 * n * 2) * 2 / n
    if t1 < 0.5:
        y1 -= 1
    else:
        y1 = 3 - y1
    t2 = t1 + 0.5
    y2 = 0.0
    for _ in range(n):
        if t1 < dt:
            y2 += sqr(t1 / dt - 1)
        elif t1 > 1 - dt:
            y2 -= sqr((t1 - 1) / dt + 1)
        t1 += dt2
        t1 -= int(t1)
        t2 -= int(t2)
        if t2 < dt:
            y2 -= sqr(t2 / dt - 1)
        elif t2 > 1 - dt:
            y2 += sqr((t2 - 1) / dt + 1)
        t2 += dt2
    return y2 / n + y1


def tri2(this_t, dt, pw):
    pw = min(max(0.0001, pw), 0.9999)
    t1 = pw * 0.5 + this_t
    t2 = t1 - pw + 1
    t1 -= int(t1)
    return ((blamp(t1, dt)
             - blamp(t2 - int(t2), dt)) / (pw - sqr(pw)) * dt
            + (t1 / pw if t1 < pw else (t1 - pw) / (pw - 1) + 1) * 2 - 1)


def trip(this_t, dt, pw):
    pw = max(pw, 1e-4)
    t3 = this_t + 0.75
    t1 = t3 + pw * 0.5
    t2 = t1 - pw
    t1 -= int(t1)
    return ((blamp(t3 - int(t3), dt) * -2
             + blamp(t2 - int(t2), dt)
             + blamp(t1, dt)) / pw * 2 * dt
            + (2 - abs(t1 * 4 / pw - 2) if t1 < pw else 0) - pw)


def trap(t0, dt, _):
    t1 = t0 + 0.125
    t3 = t1 + 0.25
    t2 = t3 + 0.25
    t4 = t2 + 0.25
    return (min(max(-1, (abs((int(t0 + 0.75) - t0) * -4 + 1) - 1) * 2), 1)
            + (blamp(t1 - int(t1), dt)
               - blamp(t2 - int(t2), dt)
               + blamp(t3 - int(t3), dt)
               - blamp(t4 - int(t4), dt)) * 4 * dt)


def trap2(t0, dt, pw):
    pw = min(0.9999, pw)
    t1 = t0 + pw * 0.25 + 0.25
    t2 = t1 + 0.5
    t3 = t2 - pw * 0.5
    t4 = t3 + 0.5
    scale = 1 / (1 - pw)
    return (min(max(-1, (abs((int(t0 + 0.75) - t0) * 4 - 1) - 1) * scale), 1)
            + (blamp(t1 - int(t1), dt)
               - blamp(t2 - int(t2), dt)
               - blamp(t3 - int(t3), dt)
               + blamp(t4 - int(t4), dt)) * scale * 2 * dt)


def sqr_(t1, dt, _):
    t2 = t1 + 0.5
    y = (1 if t1 < 0.5 else 0) * 2 - 1
    if t1 < dt:
        y -= sqr(t1 / dt - 1)
    elif t1 > 1 - dt:
        y += sqr((t1 - 1) / dt + 1)
    t2 -= int(t2)
    if t2 < dt:
        y += sqr(t2 / dt - 1)
    elif t2 > 1 - dt:
        y -= sqr((t2 - 1) / dt + 1)
    return y


def lpsqr2(t1, dt, _):
    t2 = t1 + 0.5
    y = (0.5 - sqr(t1 * 2 - 1) if t1 < 0.5 else sqr((1 - t1) * 2) - 0.5) * 2
    if t1 < dt:
        p = t1 / dt - 1
        y -= sqr(p) * p * 1.3333333333333333 * dt
    elif t1 > 1 - dt:
        p = (t1 - 1) / dt + 1
        y += sqr(p) * p * 1.3333333333333333 * dt
    t2 -= int(t2)
    if t2 < dt:
        p = t2 / dt - 1
        y += sqr(p) * p * 1.3333333333333333 * dt
    elif t2 > 1 - dt:
        p = (t2 - 1) / dt + 1
        y -= sqr(p) * p * 1.3333333333333333 * dt
    return y


def lpsqr3(t1, dt, _):
    t2 = t1 + 0.5
    if t1 < 0.5:
        p = t1 * 2 - 1
        y = sqr(p) * p * 2 + 1
    else:
        p = (1 - t1) * 2
        y = sqr(p) * p * 2 - 1
    if t1 < dt:
        p = t1 / dt - 1
        y -= sqr(p) * p * 2 * dt
    elif t1 > 1 - dt:
        p = (t1 - 1) / dt + 1
        y += sqr(p) * p * 2 * dt
    t2 -= int(t2)
    if t2 < dt:
        p = t2 / dt - 1
        y += sqr(p) * p * 2 * dt
    elif t2 > 1 - dt:
        p = (t2 - 1) / dt + 1
        y -= sqr(p) * p * 2 * dt
    return y


def intlpsqr2(this_t, dt, _):
    t0 = this_t + 0.41503077834420436
    t1 = t0 + 0.5
    t2 = t1 + 0.5
    t0 -= int(t0)
    if t0 < 0.5:
        y = 0.5
    else:
        y = -0.5
        t0 += y
    return ((((64 * t0 - 96) * t0 + 24) * t0 + 2) * y
            + (bluh(t1 - int(t1), dt)
               - bluh(t2 - int(t2), dt)) * sqr(dt)) * -0.54691816067802716


def intlpsqr3(this_t, dt, _):
    t0 = this_t + 0.37399153968783991
    t1 = t0 + 0.5
    t2 = t1 + 0.5
    t0 -= int(t0)
    if t0 < 0.5:
        y = -0.5
    else:
        y = 0.5
        t0 -= y
    return (((((64 * t0 - 128) * t0 + 96) * t0 - 16) * t0 - 2) * y
            + (bluh(t1 - int(t1), dt)
               - bluh(t2 - int(t2), dt)) * sqr(dt)) * -0.72405970419833828


def bpsqr2(t1, dt, _):
    t2 = t1 + 0.5
    t2 -= int(t2)
    return (((blamp(t1, dt)
              - blamp(t2, dt)) * 3.2
             + (bluh(t2, dt)
                - bluh(t1, dt)) * 0.53333333333333333 * dt) * dt
            + ((0.75 - t1) * t1 - 0.0625 if t1 < 0.5 else (t1 - 1.75) * t1 + 0.6875) * 12.8)


def sqrm1(t1, dt, _):
    t2 = t1 + 0.5
    return ((math.sin(t1 * TWO_PI) * (sqr(dt) * 2.0929829946242788 - 0.63661977236758134)
             + (1 if t1 < 0.5 else 0)) * 2 - 1
            + blep(t1, dt)
            - blep(t2 - int(t2), dt))


def sqr2(this_t, dt, pw):
    t1 = pw * 0.25 + this_t + 0.75
    t2 = t1 - 0.5
    t3 = t2 - pw * 0.5
    t4 = t3 + 0.5
    t1 -= int(t1)
    t3 -= int(t3)
    return ((blep(t1, dt)
             - blep(t2 - int(t2), dt)
             + blep(t3, dt)
             - blep(t4 - int(t4), dt)) * 0.5
            + (1 if t1 < 0.5 else 0) + (1 if t3 < 0.5 else 0) - 1)


def rect(t1, dt, pw):
    t2 = t1 - pw + 1
    y = ((1 if t1 < pw else 0) - pw) * 2
    if t1 < dt:
        y -= sqr(t1 / dt - 1)
    elif t1 > 1 - dt:
        y += sqr((t1 - 1) / dt + 1)
    t2 -= int(t2)
    if t2 < dt:
        y += sqr(t2 / dt - 1)
    elif t2 > 1 - dt:
        y -= sqr((t2 - 1) / dt + 1)
    return y


def saw(this_t, dt, _):
    t = this_t + 0.5
    t = frac1(t)
    return -blep(t, dt) + t * 2 - 1


def ramp(t, dt, _):
    return blep(t, dt) + t * -2 + 1


def intsaw(t, dt, _):
    return (blamp(t + 0.5 - int(t + 0.5), dt) * dt
            - ((t * -2 + 1) * (1 if t >= 0.5 else 0) + sqr(t) - 0.083333333333333333)) * -8


def cubsaw(this_t, dt, _):
    t1 = this_t + 0.5
    t0 = (1 if t1 >= 1 else 0) - t1 + 0.5
    return (bluh(t1 - int(t1), dt) * sqr(dt) * 0.0625
            + (sqr(t0) - 0.25) * t0) * 20.784609690826528


def sinsaw(t, dt, _):
    return blep(t, dt) + math.sin((0.5 - t) * PI)


def hpsaw6(this_t, dt, _):
    t = 1.5 - this_t
    t -= int(t)
    y = sqr(sqr(t) * t) * -2 + 0.28571428571428571
    if t < dt:
        p = t / dt - 1
        y -= (sqr(p) * 2 * dt + p) * p
    elif t > 1 - dt:
        p = (t - 1) / dt + 1
        y += (sqr(p) * 2 * dt + p) * p
    return y


def sawm1(t1, dt, _):
    t2 = t1 + 0.5
    t2 -= int(t2)
    return ((math.sin(t1 * TWO_PI) * (sqr(dt) * 1.0464914973121394 - 0.31830988618379067)
             + t2) * 2 - 1
            - blep(t2, dt))


def sawn(this_t, dt, n):
    n = int(n)
    t = this_t + 0.5
    t -= int(t)
    y1 = int((n + 1) * t) / n * 2 - 1
    if t < dt:
        y1 += sqr(t / dt - 1)
    elif t > 1 - dt:
        y1 -= sqr((t - 1) / dt + 1)
    dt2 = 1 / (n + 1)
    y2 = 0.0
    for _ in range(n):
        t += dt2
        t -= int(t)
        if t < dt:
            y2 -= sqr(t / dt - 1)
        elif t > 1 - dt:
            y2 += sqr((t - 1) / dt + 1)
    return y2 / n + y1


def hwsaw(this_t, dt, pw):
    pw = max(pw, 1e-4)
    t2 = this_t + 0.5
    t1 = t2 + pw
    t1 -= int(t1)
    y = (t1 * 2 / pw if t1 < pw else 0) - pw
    if t1 < dt:
        p = t1 / dt - 1
        y -= sqr(p) * p / (pw * 3) * dt
    elif t1 > 1 - dt:
        p = (t1 - 1) / dt + 1
        y += sqr(p) * p / (pw * 3) * dt
    t2 -= int(t2)
    if t2 < dt:
        p = t2 / dt - 1
        y += (sqr(p) / (pw * 3) * dt + p) * p
    elif t2 > 1 - dt:
        p = (t2 - 1) / dt + 1
        y -= (sqr(p) / (pw * 3) * dt + p) * p
    return y


def logit3(this_t, dt, _):
    t = this_t + 0.5
    t -= int(t)
    y = ((t * 8 - 12) * t + 6) * t - 1
    if t < dt:
        y += sqr(t / dt - 1)
    elif t > 1 - dt:
        y -= sqr((t - 1) / dt + 1)
    return y


def stairs(t1, dt, _):
    t2 = t1 + 0.25
    t3 = t2 + 0.25
    t4 = t3 + 0.25
    y = int(t1 * 4) * -0.66666666666666667 + 1
    if t1 < dt:
        y -= sqr(t1 / dt - 1)
    elif t1 > 1 - dt:
        y += sqr((t1 - 1) / dt + 1)
    for tk in (t2, t3, t4):
        tk -= int(tk)
        if tk < dt:
            y += sqr(tk / dt - 1) * 0.33333333333333333
        elif tk > 1 - dt:
            y -= sqr((tk - 1) / dt + 1) * 0.33333333333333333
    return y


def stairs3(t1, dt, pw):
    t2 = t1 + 0.5
    t3 = t2 - pw * 0.5
    t4 = t3 + 0.5
    if t1 < 0.5:
        y = 1 if t1 * 2 < pw else 0.33333333333333333
    else:
        y = -0.33333333333333333 if t1 * 2 - 1 < pw else -1
    y = y - pw * 0.66666666666666667 + 0.33333333333333333
    if t1 < dt:
        y -= sqr(t1 / dt - 1)
    elif t1 > 1 - dt:
        y += sqr((t1 - 1) / dt + 1)
    for tk in (t2, t3, t4):
        tk -= int(tk)
        if tk < dt:
            y += sqr(tk / dt - 1) * 0.33333333333333333
        elif tk > 1 - dt:
            y -= sqr((tk - 1) / dt + 1) * 0.33333333333333333
    return y


def stairs2(t1, dt, pw):
    t2 = t1 - pw * 0.5 + 0.5
    t3 = t2 + 0.5
    y = ((1 if t1 * 2 < pw else 0) if t1 * 2 - 1 < pw else -1) - pw + 0.5
    if t1 < dt:
        y -= sqr(t1 / dt - 1)
    elif t1 > 1 - dt:
        y += sqr((t1 - 1) / dt + 1)
    for tk in (t2, t3):
        tk -= int(tk)
        if tk < dt:
            y += sqr(tk / dt - 1) * 0.5
        elif tk > 1 - dt:
            y -= sqr((tk - 1) / dt + 1) * 0.5
    return y


def hstri(this_t, dt, fc):
    t = this_t + 0.5
    if fc > 0:
        t -= int(t)
        x = t * fc + 0.75
        y = abs((int(x + 0.75) - x) * 4 - 1) - 1
        dc = -1 / fc
        pw = dc * 0.5 + 1
        x = fc - int(fc)
        scale = abs(x * 2 - 1) - 1
        if x > 0.5:
            x -= 0.5
        else:
            dc = -dc
        y -= (sqr(x) * 2 - x) * dc
        y += blep(t, dt) * scale
        n = math.ceil(fc * 2)
        scale = fc * 4 * dt
        if n & 1:
            n -= 1
            t += pw
            scale = -scale
        for _ in range(n):
            t -= int(t)
            y += blamp(t, dt) * scale
            t += pw
            scale = -scale
        return y
    return 0.0


def hssqr(this_t, dt, fc):
    t = this_t + 0.5
    t -= int(t)
    x = t * fc
    y = (1 if x - int(x) < 0.5 else 0) * -2
    x = int(fc)
    y += 2 - x / max(fc, 0.5) if x + 0.5 > fc else (x + 1) / fc
    n = math.ceil(fc - 0.5)
    pw = 1 - 0.5 / fc if fc != 0 else -math.inf
    if x < n:
        scale = -1
    else:
        t += pw
        scale = 1
    for _ in range(max(int(n + n), 0)):
        t = frac1(t)  # the loop runs only for fc > 0.5, then pw > 0
        y += blep(t, dt) * scale
        t += pw
        scale = -scale
    return y


def hssaw(this_t, dt, fc):
    n = max(math.ceil(fc) - 1, 0)
    scale = fc - n
    pw = scale / fc if fc > 0 else 1
    t = this_t + 0.5
    t -= int(t)
    y = t * fc
    y = (y - int(y)) * 2 - sqr(pw) * fc + pw - 1
    if t < dt:
        y += sqr(t / dt - 1) * scale
    elif t > 1 - dt:
        y -= sqr((t - 1) / dt + 1) * scale
    if n > 0:
        pw = (pw - 1) / n + 1
        for _ in range(n):
            t += pw
            t -= int(t)
            if t < dt:
                y += sqr(t / dt - 1)
            elif t > 1 - dt:
                y -= sqr((t - 1) / dt + 1)
    return y


class Circ:
    """`poly_circ`, the only waveform with a state (`lp`)."""

    def __init__(self):
        self.lp = 0.0

    def __call__(self, this_t, dt, _):
        t1 = this_t
        t2 = t1 + 0.015625
        t3 = t2 + 0.5
        t2 -= int(t2)
        self.lp += (blep(t2, dt)
                    - blep(t3 - int(t3), dt)
                    + (1 if t2 < 0.5 else 0) * 2 - 1
                    - self.lp) / (dt + 0.028134884879909565) * dt
        return (self.lp + math.sin(t1 * TWO_PI)) * 0.5


# waveform name -> (function or factory, kind of the extra parameter)
WAVES = {
    "sinn": (sinn, "n"), "half": (half, None), "half2": (half2, "pw"), "full": (full, None),
    "sinp": (sinp, None), "alt": (alt, None), "camel": (camel, None), "camela": (camela, None),
    "camel2": (camel2, None), "sin2": (sin2, None), "para": (para, None), "arc": (arc, None),
    "circ": (Circ, "state"), "tri": (tri, None), "hyptri": (hyptri, None), "trin": (trin, "n"),
    "tri2": (tri2, "pw"), "trip": (trip, "pw"), "trap": (trap, None), "trap2": (trap2, "pw"),
    "sqr": (sqr_, None), "lpsqr2": (lpsqr2, None), "lpsqr3": (lpsqr3, None),
    "intlpsqr2": (intlpsqr2, None), "intlpsqr3": (intlpsqr3, None), "bpsqr2": (bpsqr2, None),
    "sqrm1": (sqrm1, None), "sqr2": (sqr2, "pw"), "rect": (rect, "pw"),
    "saw": (saw, None), "ramp": (ramp, None), "sawp": (sawp, None), "intsaw": (intsaw, None),
    "cubsaw": (cubsaw, None), "sinsaw": (sinsaw, None), "hpsaw6": (hpsaw6, None),
    "sawm1": (sawm1, None), "sawn": (sawn, "n"), "hwsaw": (hwsaw, "pw"), "logit3": (logit3, None),
    "stairs": (stairs, None), "stairs3": (stairs3, "pw"), "stairs2": (stairs2, "pw"),
    "hstri": (hstri, "f2"), "hssqr": (hssqr, "f2"), "hssaw": (hssaw, "f2"),
}

FS = 1000.0
PHASE = 0.37
AMPLITUDE = 1.5
N = 80


def wrap(x):
    r = x - math.floor(x)
    return r if r < 1.0 else 0.0


def samples(wave, f, p=None):
    """The generator: `f` (Hz) and `p` (pw / n / f2) are numbers or functions of time."""
    fn, kind = WAVES[wave]
    body = fn() if kind == "state" else fn
    fof = f if callable(f) else (lambda t: f)
    pof = p if callable(p) else (lambda t: p)
    acc = 0.0
    out = []
    for i in range(N):
        tm = i / FS
        dt = fof(tm) / FS
        ph = wrap(acc + PHASE)
        if kind == "pw":
            q = min(1.0, max(0.0, pof(tm)))
        elif kind == "n":
            q = float(p)
        elif kind == "f2":
            q = pof(tm) / (dt * FS) if dt > 0 else 0.0
        else:
            q = 0.0
        y = AMPLITUDE * gain(dt) * body(ph, dt, q)
        assert math.isfinite(y), (wave, i, y)
        out.append(y)
        acc = wrap(acc + dt)
    return out


def vec(xs):
    return "[" + " ".join(repr(float(x)) for x in xs) + "]"


def fmt(x):
    return repr(float(x)) if isinstance(x, float) else str(x)


lines = []
for wave, (_, kind) in WAVES.items():
    for dt in (0.013, 0.1, 0.23):
        f = dt * FS
        if kind in (None, "state"):
            variants = [{}]
        elif kind == "pw":
            variants = [{"pw": 0.25}, {"pw": 0.7}]
        elif kind == "n":
            variants = [{"n": 1}, {"n": 2}, {"n": 5}]
        else:
            variants = [{"f2": fc * f} for fc in (0.6, 1.5, 3.4)]
        for v in variants:
            p = v.get(kind) if kind not in (None, "state") else None
            head = "{:wave :%s :f %s" % (wave, fmt(f))
            head += "".join(" :%s %s" % (k, fmt(x)) for k, x in v.items())
            lines.append(head + " :samples " + vec(samples(wave, f, p)) + "}")

# time-varying f(t), pw(t)
f_tv = lambda t: 100.0 + 80.0 * t
pw_tv = lambda t: 0.3 + 0.4 * t
tv = [
    "{:wave :saw :samples " + vec(samples("saw", f_tv)) + "}",
    "{:wave :rect :samples " + vec(samples("rect", f_tv, pw_tv)) + "}",
    "{:wave :hssaw :f2 250.0 :samples " + vec(samples("hssaw", f_tv, 250.0)) + "}",
]

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w") as fo:
    fo.write(";; generated by utils/fastmath/dev/generate_polyblep_reference.py\n")
    fo.write("{:fs %s :phase %s :amplitude %s :n-samples %d\n" % (fmt(FS), fmt(PHASE), fmt(AMPLITUDE), N))
    fo.write(" :cases [\n  " + "\n  ".join(lines) + "]\n")
    fo.write(" :time-varying [\n  " + "\n  ".join(tv) + "]}\n")

print("wrote", OUT, len(lines), "cases")
