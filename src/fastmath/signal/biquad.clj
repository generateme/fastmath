(ns fastmath.signal.biquad
  "Biquad (second order IIR) filter coefficients from the Audio EQ Cookbook by Robert Bristow-Johnson.

  Every function returns a [[BiquadConf]] for the difference equation

  `y[n] = b0 x[n] + b1 x[n-1] + b2 x[n-2] - a1 y[n-1] - a2 y[n-2]`

  with `b0 b1 b2 a1 a2` already divided by `a0` (the cookbook sign convention, as `scipy.signal.lfilter` with `a = [1 a1 a2]`); `a0` is kept unnormalized.

  Common arguments:

  * `fc` - center or corner frequency, `0 < fc < fs/2`
  * `fs` - sample rate
  * `bw` - bandwidth in octaves (between -3 dB points for band filters, between midpoint gain points for the peaking equalizer), converted with the digital cookbook formula, see [[bandwidth->Q]]
  * `gain` - gain in dB (peaking and shelving filters)
  * `Q` - shelf quality factor; [[shelf-slope->Q]] converts the cookbook shelf slope `S`

  Throws `ex-info` when `fc` is not in `(0, fs/2)` or `bw`, `Q` are not positive.

  Reference: https://webaudio.github.io/Audio-EQ-Cookbook/Audio-EQ-Cookbook.txt"
  (:require [fastmath.core :as m]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(deftype ^{:doc "Biquad coefficients: `b0 b1 b2 a1 a2` divided by `a0`, `a0` unnormalized."}
    BiquadConf [^double b0 ^double b1 ^double b2 ^double a0 ^double a1 ^double a2])

(defn- omega
  "Normalized angular frequency `2 pi fc / fs`, checked to be in `(0, pi)`."
  ^double [^double fc ^double fs]
  (when-not (and (m/pos? fs) (m/< 0.0 fc (m/* 0.5 fs)))
    (throw (ex-info "Frequency should be in (0, fs/2)." {:fc fc :fs fs})))
  (m/* m/TWO_PI (m// fc fs)))

(defn- positive!
  ^double [what ^double v]
  (when-not (m/pos? v)
    (throw (ex-info (str what " should be positive.") {what v})))
  v)

(defn bandwidth->Q
  "Converts bandwidth `bw` in octaves to the quality factor `Q`.

  * `[bw]` - analog prototype: `1 / (2 sinh(ln(2)/2 bw))`
  * `[bw fc fs]` - digital filter (Audio EQ Cookbook, bilinear transform warping): `1 / (2 sinh(ln(2)/2 bw w0/sin(w0)))`, `w0 = 2 pi fc / fs`

  The biquad functions of this namespace use the digital form."
  (^double [^double bw] (m// (m/* 2.0 (m/sinh (m/* bw m/LN2_2)))))
  (^double [^double bw ^double fc ^double fs]
   (let [w0 (omega fc fs)]
     (m// (m/* 2.0 (m/sinh (m/* bw m/LN2_2 (m// w0 (m/sin w0)))))))))

(defn Q->bandwidth
  "Converts the quality factor `Q` to bandwidth in octaves, the inverse of [[bandwidth->Q]] (analog `[Q]` or digital `[Q fc fs]`)."
  (^double [^double Q] (m// (m/asinh (m// (m// Q) 2.0)) m/LN2_2))
  (^double [^double Q ^double fc ^double fs]
   (let [w0 (omega fc fs)]
     (m/* (Q->bandwidth Q) (m// (m/sin w0) w0)))))

(defn shelf-slope->Q
  "Converts the cookbook shelf slope `S` (1 is the steepest monotonic slope) for a shelf of `gain` dB to the quality factor `Q`.

  `1/Q = sqrt((A + 1/A)(1/S - 1) + 2)`, `A = 10^(gain/40)`. Throws `ex-info` when `S` is not positive or too large for the gain (negative square root)."
  ^double [^double gain ^double S]
  (positive! :slope S)
  (let [A (m/exp10 (m/* gain 0.025))
        r (m/+ (m/* (m/+ A (m// A)) (m/- (m// S) 1.0)) 2.0)]
    (when-not (m/pos? r)
      (throw (ex-info "Shelf slope too large for this gain." {:slope S :gain gain})))
    (m// (m/sqrt r))))

(defn- band-alpha
  "`sin(w0) / (2 Q)` for bandwidth `bw` in octaves (digital form)."
  ^double [^double w0 ^double bw]
  (positive! :bandwidth bw)
  (m/* (m/sin w0) (m/sinh (m/* bw m/LN2_2 (m// w0 (m/sin w0))))))

(defn equalizer
  "Peaking equalizer: `gain` dB at `fc`, bandwidth `bw` octaves, sample rate `fs`."
  [^double fc ^double gain ^double bw ^double fs]
  (let [w0 (omega fc fs)
        cw (m/cos w0)
        A (m/exp10 (m/* gain 0.025))
        alpha (band-alpha w0 bw)
        adivA (m// alpha A)
        amulA (m/* alpha A)

        a0 (m/inc adivA)
        a0r (m// a0)

        b0 (m/* a0r (m/inc amulA))
        b1 (m/* a0r -2.0 cw)
        b2 (m/* a0r (m/- 1.0 amulA))
        a1 b1
        a2 (m/* a0r (m/- 1.0 adivA))]
    (BiquadConf. b0 b1 b2 a0 a1 a2)))

(defn highshelf
  "High shelf: `gain` dB above `fc`, quality factor `Q` (see [[shelf-slope->Q]]), sample rate `fs`."
  [^double fc ^double gain ^double Q ^double fs]
  (let [w0 (omega fc fs)
        cw (m/cos w0)
        sw (m/sin w0)
        A (m/exp10 (m/* gain 0.025))
        A+ (m/inc A)
        A- (m/dec A)
        alpha (m/* sw (m// (m/sqrt A) (positive! :Q Q))) ;; 2 sqrt(A) sin(w0)/(2Q)
        cA+ (m/* cw A+)
        cA- (m/* cw A-)
        
        a0 (m/- (m/+ A+ alpha) cA-)
        a0r (m// a0)
        
        b0 (m/* a0r A (m/+ A+ cA- alpha))
        b1 (m/* a0r -2.0 A (m/+ A- cA+))
        b2 (m/* a0r A (m/- (m/+ A+ cA-) alpha))
        a1 (m/* a0r 2.0 (m/- A- cA+))
        a2 (m/* a0r (m/- A+ cA- alpha))]
    (BiquadConf. b0 b1 b2 a0 a1 a2)))

(defn lowshelf
  "Low shelf: `gain` dB below `fc`, quality factor `Q` (see [[shelf-slope->Q]]), sample rate `fs`."
  [^double fc ^double gain ^double Q ^double fs]
  (let [w0 (omega fc fs)
        cw (m/cos w0)
        sw (m/sin w0)
        A (m/exp10 (m/* gain 0.025))
        A+ (m/inc A)
        A- (m/dec A)
        alpha (m/* sw (m// (m/sqrt A) (positive! :Q Q))) ;; 2 sqrt(A) sin(w0)/(2Q)
        cA+ (m/* cw A+)
        cA- (m/* cw A-)
        
        a0 (m/+ A+ cA- alpha)
        a0r (m// a0)
        
        b0 (m/* a0r A (m/- (m/+ A+ alpha) cA-))
        b1 (m/* a0r 2.0 A (m/- A- cA+))
        b2 (m/* a0r A (m/- A+ cA- alpha))
        a1 (m/* a0r -2.0 (m/+ A- cA+))
        a2 (m/* a0r (m/- (m/+ A+ cA-) alpha))]
    (BiquadConf. b0 b1 b2 a0 a1 a2)))

(defn lowpass
  "Low pass filter: corner frequency `fc`, bandwidth `bw` octaves, sample rate `fs`."
  [^double fc ^double bw ^double fs]
  (let [w0 (omega fc fs)
        cw (m/cos w0)
        alpha (band-alpha w0 bw)

        a0 (m/inc alpha)
        a0r (m// a0)
        
        cw- (m/- 1.0 cw)
        
        b0 (m/* a0r 0.5 cw-)
        b1 (m/* a0r cw-)
        b2 b0
        a1 (m/* a0r -2.0 cw)
        a2 (m/* a0r (m/- 1.0 alpha))]
    (BiquadConf. b0 b1 b2 a0 a1 a2)))

(defn highpass
  "High pass filter: corner frequency `fc`, bandwidth `bw` octaves, sample rate `fs`."
  [^double fc ^double bw ^double fs]
  (let [w0 (omega fc fs)
        cw (m/cos w0)
        alpha (band-alpha w0 bw)

        a0 (m/inc alpha)
        a0r (m// a0)
        cw+ (m/inc cw)
        
        b0 (m/* a0r 0.5 cw+)
        b1 (m/* a0r (m/- cw+))
        b2 b0
        a1 (m/* a0r -2.0 cw)
        a2 (m/* a0r (m/- 1.0 alpha))]
    (BiquadConf. b0 b1 b2 a0 a1 a2)))

(defn bandpass
  "Band pass filter with 0 dB peak gain: center frequency `fc`, bandwidth `bw` octaves, sample rate `fs`."
  [^double fc ^double bw ^double fs]
  (let [w0 (omega fc fs)
        cw (m/cos w0)
        alpha (band-alpha w0 bw)
        
        a0 (m/inc alpha)
        a0r (m// a0)
        
        b0 (m/* a0r alpha)
        b1 0.0
        b2 (m/* a0r (m/- alpha))
        a1 (m/* a0r -2.0 cw)
        a2 (m/* a0r (m/- 1.0 alpha))]
    (BiquadConf. b0 b1 b2 a0 a1 a2)))

(defn notch
  "Notch (band stop) filter: center frequency `fc`, bandwidth `bw` octaves, sample rate `fs`."
  [^double fc ^double bw ^double fs]
  (let [w0 (omega fc fs)
        cw (m/cos w0)
        alpha (band-alpha w0 bw)
        
        a0 (m/inc alpha)
        a0r (m// a0)
        cw2 (m/* -2.0 cw)
        
        b0 a0r
        b1 (m/* a0r cw2)
        b2 a0r
        a1 b1
        a2 (m/* a0r (m/- 1.0 alpha))]
    (BiquadConf. b0 b1 b2 a0 a1 a2)))

(defn allpass
  "All pass filter (phase shift, unit magnitude): center frequency `fc`, bandwidth `bw` octaves, sample rate `fs`."
  [^double fc ^double bw ^double fs]
  (let [w0 (omega fc fs)
        cw (m/cos w0)
        alpha (band-alpha w0 bw)
        
        a0 (m/inc alpha)
        a0r (m// a0)
        cw2 (m/* -2.0 cw)
        
        b0 (m/* a0r (m/- 1.0 alpha))
        b1 (m/* a0r cw2)
        b2 (m/* a0r (m/inc alpha))
        a1 b1
        a2 b0]
    (BiquadConf. b0 b1 b2 a0 a1 a2)))
