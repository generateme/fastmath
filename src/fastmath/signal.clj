(ns fastmath.signal
  "Signal processing: effects (filters), spectral analysis, convolution, filtering, smoothing and signal generation.

  A signal is any sequence of doubles.

  ## Effects

  An effect is a stateful sample-by-sample filter created by the [[effect]] multimethod. Apply effects with [[apply-effects]] (any sequence) or [[apply-effects-raw]] (`double-array` only), compose them with [[compose-effects]], reset them with [[reset-effects]]. An effect can be called as a function: with a sample it returns the effect with the next state, without arguments it returns the last output.

  Every effect is created with `(effect name)` (default parameters) or `(effect name params)`. All names are listed in [[effects-list]]. Unknown keyword options (`:kind`, `:tan`, `:method`) throw `ex-info`. Effects which keep a mutable buffer (`:echo`, `:mda-thru-zero`, `:iir`) should not be shared between threads.

  #### :simple-lowpass, :simple-highpass

  One pole RC filter (highpass is `x - lowpass(x)`).

  * `:rate` - sample rate (default 44100.0)
  * `:cutoff` - cutoff frequency (default 2000.0)

  #### :biquad-eq, :biquad-hs, :biquad-ls, :biquad-lp, :biquad-hp, :biquad-bp

  Audio EQ Cookbook biquads (see `fastmath.signal.biquad`): peaking equalizer, high and low shelf, lowpass, highpass, bandpass.

  * `:fc` - center or corner frequency (default 1000.0)
  * `:gain` - gain in dB for `:biquad-eq`, `:biquad-hs`, `:biquad-ls` (default 0.0)
  * `:bw` - bandwidth in octaves for `:biquad-eq`, `:biquad-lp`, `:biquad-hp`, `:biquad-bp` (default 1.0)
  * `:slope` - shelf slope `S` for `:biquad-hs`, `:biquad-ls` (default 1.5, `S = 1` is the steepest monotonic shelf), converted with `biquad/shelf-slope->Q`
  * `:fs` - sample rate (default 44100.0)

  #### :dj-eq

  High shelf at 10 kHz, peaking equalizers at 1 kHz and 100 Hz (frequencies scaled by `rate/44100`).

  * `:hi` - high frequency gain in dB (default 0.0)
  * `:mid` - mid frequency gain in dB (default 0.0)
  * `:low` - low frequency gain in dB (default 0.0)
  * `:shelf-slope` - shelf slope for the high shelf (default 1.5)
  * `:peak-bw` - bandwidth in octaves for mid and low (default 1.0)
  * `:rate` - sample rate (default 44100.0)

  #### :phaser-allpass

  First order allpass, coefficient `(1 - delay) / (1 + delay)`.

  * `:delay` - delay factor (default 0.5)

  #### :divider

  Octave divider, swh `divider_1186`.

  * `:denom` - divisor, an integer at least 1 (default 2, the original 1)

  #### :fm

  Modulate the signal on a cosine carrier, quantize, demodulate.

  * `:quant` - quantization levels (0.0 - no quantization, default 10.0)
  * `:omega` - carrier factor (default 0.014)
  * `:phase` - deviation factor (default 0.00822)

  #### :bandwidth-limit

  CMT `lofi` bandwidth limiter (slew of `freq/rate` per sample).

  * `:rate` - sample rate (default 44100.0)
  * `:freq` - cutoff frequency (default 1000.0)

  #### :distort

  `(1 + factor) x / (factor + |x|)`

  * `:factor` - distortion factor (default 1.0)

  #### :foverdrive

  Fast overdrive, swh `foverdrive_1196`.
  
  * `:drive` - drive (default 2.0)

  #### :decimator

  Sample rate and bit depth reduction, swh `decimator_1202` (`bits` outside `[1, 31)` disable bit reduction, `fs >= rate` disables rate reduction).

  * `:bits` - bit depth (default 2.0)
  * `:fs` - decimator sample rate (default 4410.0)
  * `:rate` - input sample rate (default 44100.0)

  #### :basstreble

  Audacity `BassTreble` (2.4.2) shelves.

  * `:bass` - bass gain in dB (default 1.0, the original 0.0)
  * `:treble` - treble gain in dB (default 1.0, the original 0.0)
  * `:gain` - output gain in dB (default 0.0)
  * `:rate` - sample rate (default 44100.0)
  * `:slope` - shelf slope for both (default 0.4)
  * `:bass-freq` - bass frequency (default 250.0)
  * `:treble-freq` - treble frequency (default 4000.0)

  #### :echo

  Audacity `Echo`.

  * `:delay` - delay time in seconds, from one sample to `1e7` samples (default 0.5)
  * `:decay` - amount of echo in the signal (default 0.5)
  * `:rate` - sample rate (default 44100.0)

  #### :vcf303

  CMT `vcf303` (TB-303 like resonant filter with envelope).

  * `:rate` - sample rate (default 44100.0)
  * `:trigger` - when `true`, the envelope fires when the effect is created or reset (default `false`); the original fires on a rising edge of a trigger input, use the `reset` argument of [[apply-effects]] to retrigger
  * `:cutoff` - cutoff (0-1, default 0.8)
  * `:resonance` - resonance (0-1, default 0.8)
  * `:env-mod` - envelope modulation (0-1, default 0.5)
  * `:decay` - decay (0-1, default 1.0)
  * `:gain` - output gain (default 1.0)
  
  #### :slew-limit

  omins `slew_limiter`.

  * `:rate` - sample rate (default 44100.0)
  * `:maxrise` - maximum rise per second (default 500.0)
  * `:maxfall` - maximum fall per second (default 500.0)

  #### :mda-thru-zero

  mda ThruZero flanger (mda-vst3 processor, one channel, without per-block parameter smoothing and denormal guard). All parameters are in `[0, 1]`.

  * `:rate` - sample rate (default 44100.0)
  * `:speed` - modulation rate, `10^(3 speed - 2)` seconds per period, below 0.01 the modulation stops (default 0.3)
  * `:depth` - depth (default 0.43)
  * `:mix` - wet amount, the wet signal is subtracted (default 0.47)
  * `:depth-mod` - depth modulation (default 1.0)
  * `:feedback` - feedback (default 0.3)

  #### :iir

  `(effect :iir filter)` wraps an IIR filter created in `fastmath.signal.iir`. The filter is reset when the effect is created or reset. There is no 1-arity (throws `ex-info`).

  #### :gain

  `(effect :gain g)` or `(effect :gain {:gain g})` multiplies the signal by `g` (default 1.0).

  #### :clipping

  * `:method` - `:hard` (clamp to `[-1, 1]`), `:hyperbolic` (`tanh`), `:soft` (cubic `x - x^3/3`, `±2/3` outside `[-1, 1]`), default `:hard`
  * `:pregain` - gain applied before clipping (default 1.0)

  #### :svf1

  State variable filter, swh `svf_1214` (LADSPA). Faithful to the original: `:allpass` is `low - high` (not a unit gain allpass), when oversampling the output of the previous pass is the input of the next one, the input is scaled by `sqrt(q/2 + 0.01)` (so the lowpass DC gain is below 1), `Q` above about 1.0655 gives `NaN` in the original and throws `ex-info` here.

  * `:rate` - sample rate (default 44100.0)
  * `:cutoff` - cutoff frequency (default 1000.0)
  * `:kind` - `:lowpass` (default), `:highpass`, `:bandpass`, `:notch`, `:allpass`
  * `:Q` - resonance in `[0, 1.0655]` (default 0.5)
  * `:oversamples` - oversampling, at least 1 (default 1, the original 3)
  * `:resonance` - feedback of the band output (default 0.01, the original 0.0)

  #### :svf2

  State variable filter, Mutable Instruments stmlib `Svf`.

  * `:rate` - sample rate (default 44100.0)
  * `:cutoff` - cutoff frequency (default 1000.0)
  * `:Q` - quality factor, positive (default 0.5)
  * `:kind` - `:lowpass` (default), `:highpass`, `:bandpass`, `:bandpass-normalized` (bandpass divided by `Q`, unit gain at `cutoff`), `:notch` (lowpass + highpass)
  * `:tan` - `tan` approximation: `:exact` (default), `:dirty`, `:fast`, `:accurate`

  #### :svf3

  Linear trapezoidal state variable filter by Andrew Simper (Cytomic), Common-DSP `SvfLinearTrapOptimised2`.

  * `:rate` - sample rate (default 44100.0)
  * `:cutoff` - cutoff frequency (default 1000.0)
  * `:Q` - quality factor, positive (default 0.5)
  * `:kind` - `:lowpass` (default), `:bandpass`, `:highpass`, `:notch`, `:peak` (`high - low`, as Common-DSP), `:allpass`, `:bell`, `:lowshelf`, `:highshelf`
  * `:gaindb` - gain in dB for `:bell`, `:lowshelf`, `:highshelf` (default 0.0)
  * `:tan` - as in `:svf2`

  ## Spectral analysis

  [[fft-energy]], [[fft-magnitude]], [[fft-amplitude]], [[fft-power]] scale coefficients from `fastmath.transform/fft`; [[spectrum]], [[stft]] and [[periodogram]] analyse a signal.

  ## Convolution and filtering

  [[convolve]], [[correlate]] (direct), [[fft-convolve]], [[fft-correlate]] (FFT based), [[filter-signal]], [[filter-filter-signal]] apply FIR coefficients (`fastmath.signal.fir`) or IIR filters (`fastmath.signal.iir`).

  ## Smoothing

  Savitzky-Golay [[savgol-filter]], [[moving-average-filter]], [[kernel-smoothing-filter]].

  ## Generation

  [[waveform]], [[chirp]], [[add-waveforms]], [[gain-waveform]], [[sample-waveform]].

  ## File operations

  [[save-signal]], [[load-signal]]: 16 bit signed, big endian. Use Audacity or SoX to convert to/from audio files."
  (:require [fastmath.core :as m]
            [clojure.java.io :refer [file make-parents output-stream input-stream]]
            [fastmath.vector :as v]
            [fastmath.transform :as trans]
            [fastmath.kernel :as ker]
            [fastmath.stats :as stats]
            [fastmath.complex :as cplx]

            [fastmath.signal.waveform :as wv]
            [fastmath.signal.chirp :as chirp]
            [fastmath.signal.pad :as pad]
            [fastmath.signal.biquad :as biquad])
  (:import [fastmath.vector Vec3]
           [clojure.lang IFn IFn$OD IFn$DDDD]
           [org.apache.commons.math3.linear Array2DRowRealMatrix SingularValueDecomposition]
           [org.apache.commons.math3.util MathArrays]
           [fastmath.signal.biquad BiquadConf]
           [uk.me.berndporr.iirj Cascade]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; ## Signal processing
;;
;; Looks like a fastest implementation (reduce and vector types are 8-10 times slower).

(declare single-pass)

(deftype SampleAndState [^double sample state]
  Object
  (toString [_] (str sample))) ; sample and effect state, StateWithF or vector of StateWithF

;; Type representing list node consisting resulting sample, effect functions, state value and link to next node (or nil if last node)
(deftype EffectsList [effect-name ^double sample effect state next]
  Object
  (toString [_] (str (if next 
                       (str next " -> " (name effect-name))
                       (name effect-name)) " (" sample ")"))
  IFn
  (invoke [_] sample)
  (invoke [e n]
    (single-pass e n)))

(defn effect-node
  "Creates a single effect from the step function `f` of a stateful sample-by-sample filter.

  This is the constructor used by the methods of [[effect]]; to create an effect use [[effect]].

  Parameters:

  - `effect-name` (keyword): Name of the effect, shown when the effect is printed.
  - `f` (function): The step function. Called without arguments it returns the initial state, called with a sample and a state it returns the output sample together with the next state (`SampleAndState`).

  Returns the effect with the initial state `(f)` and the last output sample `0.0`.

  See also [[effect]], [[compose-effects]]."
  [effect-name f]
  (EffectsList. effect-name 0.0 f (f) nil))

(defn- compose-two-effects
  "Join two effect lists."
  [^EffectsList e1 ^EffectsList e2]
  (EffectsList. (.effect-name e1) (.sample e1) (.effect e1) (.state e1) (if (.next e1)
                                                                          (compose-two-effects (.next e1) e2)
                                                                          e2)))

(defn compose-effects
  "Composes effects into a single effect, as `comp` does for functions.

  `(compose-effects e1 e2 ... en)` passes a sample through `en` first and through `e1` last; the output of `e1` is the output of the composition.

  Parameters:

  - `e`, `es` (effects): One or more effects created with [[effect]].

  Returns the composed effect, which is applied, called and reset as a single effect while each part keeps its own state. With one effect the effect itself is returned. The arguments are not changed.

  See also [[effect]], [[apply-effects]], [[reset-effects]]."
  [^EffectsList e & es]
  (reduce compose-two-effects e es))

(defn reset-effects
  "Resets the state of an effect (single or composed) to the initial one.

  Parameters:

  - `e` (effect): Effect created with [[effect]] or [[compose-effects]].

  Returns a new effect with the initial state of every part and the last output sample `0.0`. The argument is not changed.

  See also [[apply-effects]] (its `reset` argument), [[compose-effects]]."
  [^EffectsList e]
  (EffectsList. (.effect-name e) 0.0 (.effect e) ((.effect e)) (when (.next e)
                                                                 (reset-effects (.next e)))))

(defn single-pass
  "Processes one sample with an effect.

  Parameters:

  - `e` (effect): Effect created with [[effect]] or [[compose-effects]].
  - `sample` (double): Input sample.

  Returns a new effect with the advanced state. The output sample is the value of the returned effect called without arguments: `((single-pass e x))`. The argument is not changed. Calling an effect with a sample, `(e x)`, is the same as `(single-pass e x)`.

  See also [[apply-effects]], [[reset-effects]]."
  [^EffectsList e ^double sample]
  (if-not (.next e)
    (let [^SampleAndState r ((.effect e) sample (.state e))]
      (EffectsList. (.effect-name e) (.sample r) (.effect e) (.state r) nil))
    (let [^EffectsList prev (single-pass (.next e) sample)
          ^SampleAndState r ((.effect e) (.sample prev) (.state e))]
      (EffectsList. (.effect-name e) (.sample r) (.effect e) (.state r) prev))))

(defn apply-effects-raw
  "Applies an effect to a signal given as `double-array`.

  Parameters:

  - `in` (`double-array`): Input signal.
  - `effects` (effect): Effect or composition of effects, see [[effect]] and [[compose-effects]].
  - `reset` (integer, optional): When positive, the state of the effect is reset (see [[reset-effects]]) after every `reset` samples. Default `0`, no reset.

  Returns a new `double-array` of the same length; the input array is not changed.

  The `effects` value is not advanced, so the same effect can be applied to several signals. The exceptions are the effects with a mutable buffer (`:echo`, `:mda-thru-zero`, `:iir`): do not reuse them for another signal.

  See also [[apply-effects]], [[single-pass]]."
  (^doubles [^doubles in effects ^long reset]
   (let [len (alength in)
         ^doubles out (double-array len)]
     (loop [idx (int 0)
            effects-and-state effects]
       (when (m/< idx len)
         (let [sample (aget in idx)
               ^EffectsList res (single-pass effects-and-state sample)
               idx+ (m/inc idx)] 
           (aset out idx ^double (.sample res))
           (recur idx+
                  (if (and (m/pos? reset) (m/zero? ^long (m/mod idx+ reset)))
                    (reset-effects effects)
                    res)))))
     out))
  (^doubles [^doubles in effects] (apply-effects-raw in effects 0)))

(defn apply-effects
  "Applies an effect to a signal given as any sequence of numbers.

  Parameters:

  - `in` (sequence of numbers): Input signal.
  - `effects` (effect): Effect or composition of effects, see [[effect]] and [[compose-effects]].
  - `reset` (integer, optional): When positive, the state of the effect is reset (see [[reset-effects]]) after every `reset` samples. Default `0`, no reset.

  Returns the filtered signal, a sequence of doubles of the same length (`nil` for an empty signal).

  See [[apply-effects-raw]] for the notes on the reuse of effects with a mutable buffer."
  ([in effects ^long reset] (m/double-array->seq (apply-effects-raw (m/seq->double-array in) effects reset)))
  ([in effects] (apply-effects in effects 0)))

;; ## Helper functions

(defn db->linear
  "Converts a gain in decibels to the linear amplitude ratio `10^(x/20)`.

  `0` dB gives `1.0`, `20` dB gives `10.0`, `-20` dB gives `0.1`, about `6.02` dB gives `2.0`. For a power ratio the exponent is `x/10` instead.

  Returns a double. See also [[linear->db]]."
  ^double [^double x]
  (m/exp10 (m// x 20.0)))

(defn linear->db
  "Converts a linear amplitude ratio to decibels: `20 log10(x)`.

  `1.0` gives `0`, `10.0` gives `20`, `2.0` gives about `6.02`. `0.0` gives `-Infinity` and a negative `x` gives `NaN`. For a power ratio the factor is `10` instead of `20`.

  Returns a double. See also [[db->linear]]."
  ^double [^double x]
  (m/* 20.0 (m/log10 x)))

;; ## Effects / Filters

(defmulti effect
  "Creates an effect: a stateful sample-by-sample filter, given its name.

  Parameters:

  - `name` (keyword): Name of the effect, one of [[effects-list]].
  - `params` (map, optional): Parameters of the effect, missing keys take the defaults. All parameters are listed for each effect in the docstring of the namespace. `:gain` also accepts the gain number and `:iir` requires an IIR filter from `fastmath.signal.iir`.

  Returns the effect. An effect is a function: called with a sample it returns the effect with the next state, called without arguments it returns the last output sample (`0.0` for a new effect). It prints as the names of its parts and the last sample.

  Throws `IllegalArgumentException` for an unknown name and `ex-info` for an invalid parameter (for example an unknown `:kind`).

  See also [[apply-effects]], [[compose-effects]], [[reset-effects]], [[single-pass]]."
  (fn [m & _] m))

;; ### Simple Low/High pass filters
(defn- calc-filter-alpha
  "Calculate alpha factor"
  ^double [^double rate ^double cutoff]
  (let [tinterval (m// rate)
        tau (m// (m/* cutoff m/TWO_PI))]
    (m// tinterval (m/+ tau tinterval))))

(defmethod effect :simple-lowpass
  ([m] (effect m {}))
  ([m {:keys [rate cutoff]
       :or {rate 44100.0 cutoff 2000.0}}]
   (let [alpha (calc-filter-alpha rate cutoff)] 
     (effect-node m (fn
                      ([^double sample ^double prev]
                       (let [s1 (m/* sample alpha)
                             s2 (m/- prev (m/* prev alpha))
                             nprev (m/+ s1 s2)]
                         (SampleAndState. nprev nprev)))
                      ([] 0.0))))))

(defmethod effect :simple-highpass
  ([m] (effect m {}))
  ([m conf]
   (let [lpfilter (effect :simple-lowpass conf)]
     (effect-node m (fn
                      ([^double sample lp]
                       (let [^EffectsList res (single-pass lp sample)]
                         (SampleAndState. (m/- sample (.sample res)) res)))
                      ([] lpfilter))))))


;; Store state in `StateBiquad` type.
(deftype StateBiquad [^double x2 ^double x1 ^double y2 ^double y1])

(defn- make-biquad-filter
  "Create biquad effect based on passed configuration: `y = b0 x + b1 x1 + b2 x2 - a1 y1 - a2 y2`."
  [m ^BiquadConf c]
  (effect-node m (fn
                   ([^double sample ^StateBiquad state]
                    (let [y (-> (m/* (.b0 c) sample)
                                (m/+ (m/* (.b1 c) (.x1 state)))
                                (m/+ (m/* (.b2 c) (.x2 state)))
                                (m/- (m/* (.a1 c) (.y1 state)))
                                (m/- (m/* (.a2 c) (.y2 state))))]
                      (SampleAndState. y (StateBiquad. (.x1 state) sample (.y1 state) y))))
                   ([] (StateBiquad. 0.0 0.0 0.0 0.0)))))

;; ### Biquad equalizer

(defmethod effect :biquad-eq
  ([m] (effect m {}))
  ([m {:keys [fc gain bw fs]
       :or {fc 1000.0 gain 0.0 bw 1.0 fs 44100.0}}]
   (make-biquad-filter m (biquad/equalizer fc gain bw fs))))

;; ### Biquad high/low shelf, `:slope` is the cookbook shelf slope S

(defmethod effect :biquad-hs
  ([m] (effect m {}))
  ([m {:keys [fc gain slope fs]
       :or {fc 1000.0 gain 0.0 slope 1.5 fs 44100.0}}]
   (make-biquad-filter m (biquad/highshelf fc gain (biquad/shelf-slope->Q gain slope) fs))))

(defmethod effect :biquad-ls
  ([m] (effect m {}))
  ([m {:keys [fc gain slope fs]
       :or {fc 1000.0 gain 0.0 slope 1.5 fs 44100.0}}]
   (make-biquad-filter m (biquad/lowshelf fc gain (biquad/shelf-slope->Q gain slope) fs))))

;; ### Biquad lowpass/highpass/bandpass

(defn- lhb-params
  "Create parameters for lp hp and bp biquad filters."
  [f {:keys [fc bw fs]
      :or {fc 1000.0 bw 1.0 fs 44100.0}}]
  (f fc bw fs))

(defmethod effect :biquad-lp ([m] (effect m {})) ([m conf] (make-biquad-filter m (lhb-params biquad/lowpass conf))))
(defmethod effect :biquad-hp ([m] (effect m {})) ([m conf] (make-biquad-filter m (lhb-params biquad/highpass conf))))
(defmethod effect :biquad-bp ([m] (effect m {})) ([m conf] (make-biquad-filter m (lhb-params biquad/bandpass conf))))

;; ### DJ Equalizer

(defmethod effect :dj-eq
  ([m] (effect m {}))
  ([m {:keys [hi mid low shelf-slope peak-bw ^double rate]
       :or {hi 0.0 mid 0.0 low 0.0 shelf-slope 1.5 peak-bw 1.0 rate 44100.0}}]
   (let [b (compose-effects
            (effect :biquad-hs {:fc (m/* rate (m// 10000.0 44100.0)) :gain hi :slope shelf-slope :fs rate})
            (effect :biquad-eq {:fc (m/* rate (m// 1000.0 44100.0)) :gain mid :bw peak-bw :fs rate})
            (effect :biquad-eq {:fc (m/* rate (m// 100.0 44100.0)) :gain low :bw peak-bw :fs rate}))]
     (effect-node m (fn
                      ([sample state]
                       (let [^EffectsList res (single-pass state sample)]
                         (SampleAndState. (.sample res) res)))
                      ([] b))))))

;; ### Phaser

(defmethod effect :phaser-allpass
  ([m] (effect m {}))
  ([m {:keys [^double delay]
       :or {delay 0.5}}]
   (let [a1 (m// (m/- 1.0 delay) (m/inc delay))]
     (effect-node m (fn
                      ([^double sample ^double zm1]
                       (let [y (m/+ zm1 (m/* sample (m/- a1)))
                             new-zm1 (m/+ sample (m/* y a1))]
                         (SampleAndState. y new-zm1)))
                      ([] 0.0))))))

;; ### Divider

(deftype StateDivider [^double out ^double amp ^double count ^double lamp ^double last ^long zeroxs])

(defmethod effect :divider
  ([m] (effect m {}))
  ([m {:keys [denom]
       :or {denom 2}}]
   (when-not (and (number? denom) (m/== (double denom) (m/floor (double denom))) (m/>= (double denom) 1.0))
     (throw (ex-info "Divider denominator should be an integer at least 1." {:denom denom})))
   (let [denom (long denom)]
     (effect-node m
                  (fn
                    ([^double sample ^StateDivider state]
                     (let [count (m/inc (.count state))
                           ^StateDivider s1 (if (or (and (m/> sample 0.0) (m/<= (.last state) 0.0))
                                                    (and (m/neg? sample) (m/>= (.last state) 0.0)))
                                              (if (m/== denom 1)
                                                (StateDivider. (if (m/pos? (.out state)) -1.0 1.0) 0.0 0.0 (m// (.amp state) count) (.last state) 0)
                                                (StateDivider. (.out state) (.amp state) count (.lamp state) (.last state) (m/inc (.zeroxs state))))
                                              (StateDivider. (.out state) (.amp state) count (.lamp state) (.last state) (.zeroxs state)))
                           amp (m/+ (.amp s1) (m/abs sample))
                           ^StateDivider s2 (if (and (m/> denom 1)
                                                     (m/== (m/rem (.zeroxs s1) denom) (m/dec denom)))
                                              (StateDivider. (if (m/pos? (.out s1)) -1.0 1.0) 0.0 0.0 (m// amp (.count s1)) (.last s1) 0)
                                              (StateDivider. (.out s1) amp (.count s1) (.lamp s1) (.last s1) (.zeroxs s1)))]
                       (SampleAndState. (m/* (.out s2) (.lamp s2)) (StateDivider. (.out s2) (.amp s2) (.count s2) (.lamp s2) sample (.zeroxs s2)))))
                    ([] (StateDivider. 1.0 0.0 0.0 0.0 0.0 0)))))))

;; ### FM filter

(deftype StateFm [^double pre ^double integral ^double t lp])

(defmethod effect :fm
  ([m] (effect m {}))
  ([m {:keys [^double quant ^double omega ^double phase]
       :or {quant 10.0 omega 0.014 phase 0.00822}}]
   (let [lp-chain (compose-effects (effect :simple-lowpass {:rate 100000.0 :cutoff 25000.0})
                                   (effect :simple-lowpass {:rate 100000.0 :cutoff 10000.0})
                                   (effect :simple-lowpass {:rate 100000.0 :cutoff 1000.0}))]
     (effect-node m
                  (fn
                    ([^double sample ^StateFm state]
                     (let [sig (m/* sample phase)
                           new-integral (m/+ (.integral state) sig)
                           m (m/cos (m/+ new-integral (m/* omega (.t state))))
                           m (if (m/pos? quant)
                               (m/norm (unchecked-int (m/norm m -1.0 1.0 0.0 quant)) 0.0 quant -1.0 1.0)
                               m)
                           dem (m/abs (m/- m (.pre state)))
                           ^EffectsList res (single-pass (.lp state) dem)
                           demf (m// (m/* 2.0 (m/- (.sample res) omega)) phase)]
                       (SampleAndState. (m/constrain demf -1.0 1.0) (StateFm. m new-integral (m/inc (.t state)) res))))
                    ([] (StateFm. 0.0 0.0 0.0 lp-chain)))))))


;; ### Bandwidth limit

(defmethod effect :bandwidth-limit
  ([m] (effect m {}))
  ([m {:keys [^double freq ^double rate]
       :or {freq 1000.0 rate 44100.0}}]
   (let [dx (m// freq rate)]
     (effect-node m (fn
                      ([^double sample ^double state]
                       (let [res (if (m/>= sample state)
                                   (m/min (m/+ state dx) sample)
                                   (m/max (m/- state dx) sample))]
                         (SampleAndState. res res)))
                      ([] 0.0))))))

;; ### Distortion

(defmethod effect :distort
  ([m] (effect m {}))
  ([m {:keys [^double factor]
       :or {factor 1.0}}]
   (let [nfact (m/inc factor)]
     (effect-node m (fn 
                      ([^double sample state]
                       (let [div (m/+ factor (m/abs sample))
                             res (m/* nfact (m// sample div))]
                         (SampleAndState. res state)))
                      ([]))))))

;; ### Fast overdrive

(defmethod effect :foverdrive
  ([m] (effect m {}))
  ([m {:keys [^double drive]
       :or {drive 2.0}}]
   (let [drivem1 (m/dec drive)]
     (effect-node m (fn
                      ([^double sample state]
                       (let [fx (m/abs sample)
                             res (m// (m/* sample (m/+ fx drive)) (m/inc (m/+ (m/* sample sample) (m/* fx drivem1))))]
                         (SampleAndState. res state)))
                      ([]))))))

;; ### Decimator (swh `decimator_1202`)

(deftype StateDecimator [^double count ^double last])

(defmethod effect :decimator
  ([m] (effect m {}))
  ([m {:keys [^double bits ^double fs ^double rate]
       :or {bits 2.0 fs 4410.0 rate 44100.0}}]
   (let [step (if (or (m/>= bits 31.0) (m/< bits 1.0)) 0.0 (m/pow 0.5 (m/- bits 0.999)))
         stepr (if (m/zero? step) 1.0 (m// step))
         ratio (if (m/>= fs rate) 1.0 (m// fs rate))]
     (effect-node m (fn
                      ([^double sample ^StateDecimator state]
                       (let [ncount (m/+ (.count state) ratio)]
                         (if (m/>= ncount 1.0)
                           (let [y (m/* stepr (m/+ sample (m/* (if (m/neg? sample) -1.0 1.0) step 0.5)))
                                 delta (m/* step (m/- y (m/trunc y))) ;; C `modf` fractional part
                                 last (m/- sample delta)]
                             (SampleAndState. last (StateDecimator. (m/dec ncount) last)))
                           (SampleAndState. (.last state) (StateDecimator. ncount (.last state))))))
                      ([] (StateDecimator. 0.0 0.0)))))))

;; ### BassTreble

(deftype StateBassTreble [^double xn1Bass ^double xn2Bass ^double yn1Bass ^double yn2Bass
                          ^double xn1Treble ^double xn2Treble ^double yn1Treble ^double yn2Treble])

(defmethod effect :basstreble
  ([m] (effect m {}))
  ([m {:keys [^double bass ^double treble ^double gain ^double rate ^double slope ^double bass-freq ^double treble-freq]
       :or {bass 1.0 treble 1.0 gain 0.0 rate 44100.0 slope 0.4 bass-freq 250.0 treble-freq 4000.0}}]
   (let [data-gain (db->linear gain)
         wb (m// (m/* m/TWO_PI bass-freq) rate)
         wt (m// (m/* m/TWO_PI treble-freq) rate) 
         cwb (m/cos wb)
         cwt (m/cos wt)
         ab (m/exp (m// (m/* 2.302585092994046 bass) 40.0))
         ab+ (m/inc ab)
         ab- (m/dec ab)
         at (m/exp (m// (m/* 2.302585092994046 treble) 40.0))
         at+ (m/inc at)
         at- (m/dec at)
         bb (m/sqrt (m/- (m// (m/inc (m/sq ab)) slope) (m/sq (m/dec ab))))
         bt (m/sqrt (m/- (m// (m/inc (m/sq at)) slope) (m/sq (m/dec at))))
         bswb (m/* bb (m/sin wb))
         bswt (m/* bt (m/sin wt))

         b0b (m/* ab (m/+ (m/- ab+ (m/* ab- cwb)) bswb))
         b1b (m/* 2.0 ab (m/- ab- (m/* ab+ cwb)))
         b2b (m/* ab (m/- (m/- ab+ (m/* ab- cwb)) bswb))
         a0b (m/+ ab+ (m/* ab- cwb) bswb)
         a1b (m/* -2.0 (m/+ ab- (m/* ab+ cwb)))
         a2b (m/- (m/+ ab+ (m/* ab- cwb)) bswb)

         b0t (m/* at (m/+ at+ (m/* at- cwt) bswt))
         b1t (m/* -2.0 at (m/+ at- (m/* at+ cwt)))
         b2t (m/* at (m/- (m/+ at+ (m/* at- cwt)) bswt))
         a0t (m/+ (m/- at+ (m/* at- cwt)) bswt)
         a1t (m/* 2.0 (m/- at- (m/* at+ cwt)))
         a2t (m/- (m/- at+ (m/* at- cwt)) bswt)]
     (effect-node m (fn
                      ([^double sample ^StateBassTreble state]
                       (let [outb (m// (-> (m/* b0b sample)
                                         (m/+ (m/* b1b (.xn1Bass state)))
                                         (m/+ (m/* b2b (.xn2Bass state)))
                                         (m/- (m/* a1b (.yn1Bass state)))
                                         (m/- (m/* a2b (.yn2Bass state)))) a0b)
                             outt (m// (-> (m/* b0t outb)
                                         (m/+ (m/* b1t (.xn1Treble state)))
                                         (m/+ (m/* b2t (.xn2Treble state)))
                                         (m/- (m/* a1t (.yn1Treble state)))
                                         (m/- (m/* a2t (.yn2Treble state)))) a0t)]
                         (SampleAndState. (m/* outt data-gain)
                                          (StateBassTreble. sample (.xn1Bass state) outb (.yn1Bass state)
                                                            outb (.xn1Treble state) outt (.yn1Treble state)))))
                      ([] (StateBassTreble. 0.0 0.0 0.0 0.0
                                            0.0 0.0 0.0 0.0)))))))

;; ### Echo (audacity)

(deftype StateEcho [^doubles buffer ^long position])

(defmethod effect :echo
  ([m] (effect m {}))
  ([m {:keys [^double delay ^double decay ^double rate]
       :or {delay 0.5 decay 0.5 rate 44100.0}}]
   (let [samples (m/* delay rate)
         buffer-len (long samples)]
     (when-not (m/<= 1.0 samples 1.0e7)
       (throw (ex-info "Echo delay should be from one sample to 1e7 samples long." {:delay delay :rate rate})))
     (effect-node m (fn
                      ([^double sample ^StateEcho state]
                       (let [^doubles buffer (.buffer state)
                             position (.position state)
                             result (m/+ sample (m/* decay (aget buffer position)))]
                         (aset buffer position result)
                         (SampleAndState. result (StateEcho. buffer (m/rem (m/inc position) buffer-len)))))
                      ([]
                       (StateEcho. (double-array buffer-len 0.0) 0)))))))

;; ### Vcf303
;;
(deftype StateVcf303 [^double d1 ^double d2 ^double c0 ^long env-pos ^Vec3 abc])

(defmethod effect :vcf303
  ([m] (effect m {}))
  ([m {:keys [^double rate trigger ^double cutoff ^double resonance ^double env-mod ^double decay ^double gain]
       :or {rate 44100.0 trigger false cutoff 0.8 resonance 0.8 env-mod 0.5 decay 1.0 gain 1.0}}]
   (let [scale (m// m/PI rate)
         e0 (m/* scale
               (m/exp (-> (m/- 5.613 (m/* 0.8 env-mod))
                          (m/+ (m/* 2.1553 cutoff))
                          (m/- (m/* 0.7696 (m/- 1.0 resonance))))))        
         d (m/pow (->> decay
                       (m/* 2.3)
                       (m/+ 0.2)
                       (m/* rate)
                       (m// 1.0)
                       (m/pow 0.1)) 64.0)
         r (m/exp (m/- (m/* 3.455 resonance) 1.20))
         recalc-abc (fn [^double vc0]
                      (let [whopping (m/+ e0 vc0)
                            k (m/exp (m// (m/- whopping) r))
                            a (m/* (m/+ k k) (m/cos (m/+ whopping whopping)))
                            b (m/* (m/- k) k)
                            c (m/* 0.2 (m/- (m/- 1.0 a) b))]
                        (Vec3. a b c)))
         init-c0       (if trigger
                         (m/- (m/* scale
                               (m/exp (-> (m/+ 6.109 (m/* 1.5876 env-mod))
                                          (m/+ (m/* 2.1553 cutoff))
                                          (m/- (m/* 1.2 (m/- 1.0 resonance)))))) e0)
                         0.0)]    
     (effect-node m (fn
                      ([^double sample ^StateVcf303 state]
                       (let [^Vec3 abc (.abc state)
                             result (-> (m/* (.x abc) (.d1 state))
                                        (m/+ (m/* (.y abc) (.d2 state)))
                                        (m/+ (m/* (.z abc) sample)))
                             d2 (.d1 state)
                             d1 result
                             env-pos (m/inc (.env-pos state))]
                         (if (m/>= env-pos 64)
                           (let [c0 (m/* d (.c0 state))]
                             (SampleAndState. (m/* gain result)
                                              (StateVcf303. d1 d2 c0 0 (recalc-abc c0))))
                           (SampleAndState. (m/* gain result)
                                            (StateVcf303. d1 d2 (.c0 state) env-pos abc)))))
                      ([] (StateVcf303. 0.0 0.0 init-c0 0 (recalc-abc init-c0))))))))

;; ### Slew limiter

(defmethod effect :slew-limit
  ([m] (effect m {}))
  ([m {:keys [^double rate ^double maxrise ^double maxfall]
       :or {rate 44100.0 maxrise 500.0 maxfall 500.0}}]
   (let [maxinc (m// maxrise rate)
         maxdec (m/- (m// maxfall rate))]
     (effect-node m (fn
                      ([^double sample ^double prev]
                       (let [increment (m/- sample prev) 
                             nsample (m/+ prev (m/constrain increment maxdec maxinc))]
                         (SampleAndState. nsample nsample)))
                      ([] 0.0))))))


;; ### mda ThruZero (mda-vst3 `mdaThruZeroProcessor`), one channel

(deftype StateMdaThruZero [^doubles buffer ^double ph ^long bp ^double f])

(defmethod effect :mda-thru-zero
  ([m] (effect m {}))
  ([m {:keys [^double rate ^double speed ^double depth ^double mix ^double depth-mod ^double feedback]
       :or {rate 44100.0 speed 0.3 depth 0.43 mix 0.47 feedback 0.3 depth-mod 1.0}}]
   (let [rat (if (m/< speed 0.01) 0.0 (m// (m/* (m/pow 10.0 (m/- (m/* 3.0 speed) 2.0)) 2.0) rate))
         dep (m/* 2000.0 (m/sq depth))
         dem (m/- dep (m/* dep depth-mod))
         dep (m/- dep dem)
         wet mix
         dry (m/- 1.0 wet)
         fb (m/- (m/* 1.9 feedback) 0.95)]
     (effect-node m (fn
                      ([^double sample ^StateMdaThruZero state]
                       (let [^doubles buffer (.buffer state)
                             ph (m/+ (.ph state) rat)
                             ph (if (m/> ph 1.0) (m/- ph 2.0) ph)
                             bp (m/bit-and (m/long-dec (.bp state)) 0x7ff)]
                         (aset buffer bp (m/+ sample (m/* fb (.f state))))
                         (let [tmpf (m/+ dem (m/* dep (m/- 1.0 (m/sq ph))))
                               tmp (long tmpf)
                               tmpf (m/- tmpf tmp)
                               tmp (m/bit-and (m/long-add tmp bp) 0x7ff)
                               tmpi (m/bit-and (m/long-inc tmp) 0x7ff)
                               f (aget buffer tmp)
                               f (m/+ (m/* tmpf (m/- (aget buffer tmpi) f)) f)
                               result (m/- (m/* sample dry) (m/* f wet))]
                           (SampleAndState. result (StateMdaThruZero. buffer ph bp f)))))
                      ([] (StateMdaThruZero. (double-array 2048) 0.0 0 0.0)))))))

;; ### IIR filter, gain and clipping

(defmethod effect :iir
  ([m] (throw (ex-info "The :iir effect requires an IIR filter, see `fastmath.signal.iir`." {:effect m})))
  ([m iir]
   (when-not (instance? Cascade iir)
     (throw (ex-info "The :iir effect requires an IIR filter, see `fastmath.signal.iir`." {:effect m :filter iir})))
   (let [^Cascade iir iir]
     (effect-node m (fn ([^double sample ^Cascade state]
                        (SampleAndState. (.filter state sample) state))
                      ([] (.reset iir) iir))))))

(defmethod effect :gain
  ([m] (effect m 1.0))
  ([m gain]
   (let [gain (double (if (map? gain) (get gain :gain 1.0) gain))]
     (effect-node m (fn ([^double sample _] (SampleAndState. (m/* gain sample) nil))
                      ([] nil))))))

(defmethod effect :clipping
  ([m] (effect m {}))
  ([m {:keys [method ^double pregain]
       :or {method :hard pregain 1.0}}]
   (let [f (case method
             :hard (fn ^double [^double v] (m/constrain v -1.0 1.0))
             :hyperbolic (fn ^double [^double v] (m/tanh v))
             :soft (fn ^double [^double v] (cond
                                            (m/< v -1.0) -0.6666666666666666
                                            (m/> v 1.0) 0.6666666666666666
                                            :else (m/- v (m/* m/THIRD (m/cb v)))))
             (throw (ex-info "Unknown clipping method." {:method method :allowed [:hard :hyperbolic :soft]})))]
     (effect-node m (fn ([^double sample _]
                        (SampleAndState. (f (m/* pregain sample)) nil))
                      ([] nil))))))

;; ### State variable filters

(defn- positive-Q!
  [Q]
  (when-not (and (number? Q) (m/pos? (double Q)))
    (throw (ex-info "Q should be positive." {:Q Q})))
  Q)

(defn- svf-unknown
  [what v allowed]
  (throw (ex-info (str "Unknown " what ".") {what v :allowed allowed})))

;; swh `svf_1214` (LADSPA)

(deftype SVF1State [^double h ^double b ^double l ^double p ^double n])

(defmethod effect :svf1
  ([m] (effect m {}))
  ([m {:keys [rate cutoff kind Q oversamples resonance]
       :or {rate 44100.0 cutoff 1000.0 kind :lowpass Q 0.5 oversamples 1 resonance 0.01}}]
   (let [rate (double rate)
         cutoff (double cutoff)
         Q (double Q)
         oversamples (long oversamples)
         resonance (double resonance)
         ^IFn$OD fun (case kind
               :lowpass (fn ^double [^SVF1State s] (.l s))
               :highpass (fn ^double [^SVF1State s] (.h s))
               :bandpass (fn ^double [^SVF1State s] (.b s))
               :notch (fn ^double [^SVF1State s] (.n s))
               :allpass (fn ^double [^SVF1State s] (.p s))
               (svf-unknown :kind kind [:lowpass :highpass :bandpass :notch :allpass]))
         f (m/* 2.0 (m/sinpi (m// cutoff (m/* rate oversamples))))
         q (m/* 2.0 (m/cos (m/* m/HALF_PI (m/pow Q 0.1))))
         qn2 (m/+ (m/* 0.5 q) 0.01)]
     (when (m/neg? qn2)
       (throw (ex-info "Q is too large for :svf1 (the original filter returns NaN), maximum is about 1.0655." {:Q Q})))
     (when (m/< oversamples 1)
       (throw (ex-info "Oversampling should be at least 1." {:oversamples oversamples})))
     (let [qn (m/sqrt qn2)]
       (effect-node m (fn ([^double sample ^SVF1State state]
                          (loop [i 0
                                 in (m/* qn (m/+ sample (m/* (.b state) resonance)))
                                 ^SVF1State state state]
                            (if (m/== i oversamples)
                              (SampleAndState. (.invokePrim fun state) state)
                              (let [b (m/- (.b state) (m/* 0.001 (m/cb (.b state))))
                                    h (m/- in (.l state) (m/* b q))
                                    b (m/+ b (m/* f h))
                                    l (m/+ (.l state) (m/* f b))
                                    n (m/+ l h)
                                    p (m/- l h)
                                    nstate (SVF1State. h b l p n)]
                                (recur (m/inc i) (.invokePrim fun nstate) nstate)))))
                        ([] (SVF1State. 0.0 0.0 0.0 0.0 0.0))))))))

;; stmlib `Svf` (Mutable Instruments), https://github.com/pichenettes/stmlib/blob/master/dsp/filter.h

(def ^:const ^:private SVF2-M_PI_POW_2  (m/* m/M_PI m/M_PI))
(def ^:const ^:private SVF2-M_PI_POW_3  (m/* SVF2-M_PI_POW_2 m/M_PI))
(def ^:const ^:private SVF2-M_PI_POW_5  (m/* SVF2-M_PI_POW_3 SVF2-M_PI_POW_2))
(def ^:const ^:private SVF2-M_PI_POW_7  (m/* SVF2-M_PI_POW_5 SVF2-M_PI_POW_2))
(def ^:const ^:private SVF2-M_PI_POW_9  (m/* SVF2-M_PI_POW_7 SVF2-M_PI_POW_2))
(def ^:const ^:private SVF2-M_PI_POW_11 (m/* SVF2-M_PI_POW_9 SVF2-M_PI_POW_2))

(def ^:const ^:private SVF2-DIRTY_A    (m/* 3.736e-01 SVF2-M_PI_POW_3))
(def ^:const ^:private SVF2-FAST_A     (m/* 3.260e-01 SVF2-M_PI_POW_3))
(def ^:const ^:private SVF2-FAST_B     (m/* 1.823e-01 SVF2-M_PI_POW_5))
(def ^:const ^:private SVF2-ACCURATE_A (m/* 3.333314036e-01 SVF2-M_PI_POW_3))
(def ^:const ^:private SVF2-ACCURATE_B (m/* 1.333923995e-01 SVF2-M_PI_POW_5))
(def ^:const ^:private SVF2-ACCURATE_C (m/* 5.33740603e-02 SVF2-M_PI_POW_7))
(def ^:const ^:private SVF2-ACCURATE_D (m/* 2.900525e-03 SVF2-M_PI_POW_9))
(def ^:const ^:private SVF2-ACCURATE_E (m/* 9.5168091e-03 SVF2-M_PI_POW_11))

(defn- tan-exact
  ^double [^double f]
  (m/tan (m/* m/M_PI (m/min f 0.497))))

(defn- tan-dirty
  ^double [^double f]
  (m/* f (m/+ m/M_PI (m/* f f SVF2-DIRTY_A))))

(defn- tan-fast
  ^double [^double f]
  (let [f2 (m/* f f)]
    (->> SVF2-FAST_B
         (m/* f2)
         (m/+ SVF2-FAST_A)
         (m/* f2)
         (m/+ m/M_PI)
         (m/* f))))

(defn- tan-accurate
  ^double [^double f]
  (let [f2 (m/* f f)]
    (->> SVF2-ACCURATE_E
         (m/* f2)
         (m/+ SVF2-ACCURATE_D)
         (m/* f2)
         (m/+ SVF2-ACCURATE_C)
         (m/* f2)
         (m/+ SVF2-ACCURATE_B)
         (m/* f2)
         (m/+ SVF2-ACCURATE_A)
         (m/* f2)
         (m/+ m/M_PI)
         (m/* f))))

(defn- svf-g
  "stmlib `tan` approximations of `tan(pi f)`."
  ^double [tan ^double ratio]
  (case tan
    :exact (tan-exact ratio)
    :dirty (tan-dirty ratio)
    :fast (tan-fast ratio)
    :accurate (tan-accurate ratio)
    (svf-unknown :tan tan [:exact :dirty :fast :accurate])))

(deftype SVF2State [^double state1 ^double state2])

(defmethod effect :svf2
  ([m] (effect m {}))
  ([m {:keys [rate cutoff Q kind tan]
       :or {rate 44100.0 cutoff 1000.0 kind :lowpass Q 0.5 tan :exact}}]
   (let [g (svf-g tan (m// (double cutoff) (double rate)))
         r (m// (double (positive-Q! Q)))
         h (m// (m/+ 1.0 (m/* g (m/+ r g))))
         ^IFn$DDDD out (case kind
               :lowpass (fn ^double [^double lp ^double _bp ^double _hp] lp)
               :highpass (fn ^double [^double _lp ^double _bp ^double hp] hp)
               :bandpass (fn ^double [^double _lp ^double bp ^double _hp] bp)
               :bandpass-normalized (fn ^double [^double _lp ^double bp ^double _hp] (m/* r bp))
               :notch (fn ^double [^double lp ^double _bp ^double hp] (m/+ lp hp))
               (svf-unknown :kind kind [:lowpass :highpass :bandpass :bandpass-normalized :notch]))]
     (effect-node m (fn ([^double sample ^SVF2State state]
                        (let [s1 (.state1 state)
                              s2 (.state2 state)
                              hp (m/* h (m/- sample (m/* r s1) (m/* g s1) s2))
                              bp (m/+ (m/* g hp) s1)
                              lp (m/+ (m/* g bp) s2)]
                          (SampleAndState. (.invokePrim out lp bp hp)
                                           (SVF2State. (m/+ (m/* g hp) bp) (m/+ (m/* g bp) lp)))))
                      ([] (SVF2State. 0.0 0.0)))))))


;; Cytomic (Andrew Simper) linear trapezoidal SVF,
;; https://github.com/FredAntonCorvest/Common-DSP/blob/master/Filter/SvfLinearTrapOptimised2.hpp

(deftype SVF3State [^double ic1 ^double ic2])

(defn- compute-a
  ^Vec3 [^double g ^double k]
  (let [a1 (m// (m/+ 1.0 (m/* g (m/+ g k))))
        a2 (m/* g a1)]
    (Vec3. a1 a2 (m/* g a2))))

(defn- compute-m
  ^Vec3 [kind ^double k ^double A]
  (case kind
    :lowpass (Vec3. 0.0 0.0 1.0)
    :bandpass (Vec3. 0.0 1.0 0.0)
    :highpass (Vec3. 1.0 (m/- k) -1.0)
    :notch (Vec3. 1.0 (m/- k) 0.0)
    :peak (Vec3. 1.0 (m/- k) -2.0)
    :allpass (Vec3. 1.0 (m/* -2.0 k) 0.0)
    :bell (Vec3. 1.0 (m/* k (m/dec (m/* A A))) 0.0)
    :lowshelf (Vec3. 1.0 (m/* k (m/dec A)) (m/dec (m/* A A)))
    :highshelf (Vec3. (m/* A A) (m/* k (m/- 1.0 A) A) (m/- 1.0 (m/* A A)))
    (svf-unknown :kind kind [:lowpass :bandpass :highpass :notch :peak :allpass :bell :lowshelf :highshelf])))

(defmethod effect :svf3
  ([m] (effect m {}))
  ([m {:keys [rate cutoff Q kind tan gaindb]
       :or {rate 44100.0 cutoff 1000.0 kind :lowpass Q 0.5 tan :exact gaindb 0.0}}]
   (let [Q (double (positive-Q! Q))
         A (m/exp10 (m// (double gaindb) 40.0))
         Asqrt (m/sqrt A)
         g (svf-g tan (m// (double cutoff) (double rate)))
         k (if (= kind :bell) (m// (m/* A Q)) (m// Q))
         ^Vec3 mv (compute-m kind k A)
         ^Vec3 av (case kind
                    :lowshelf (compute-a (m// g Asqrt) k)
                    :highshelf (compute-a (m/* g Asqrt) k)
                    (compute-a g k))]
     (effect-node m (fn ([^double sample ^SVF3State state]
                        (let [v3 (m/- sample (.ic2 state))
                              v1 (m/+ (m/* (.x av) (.ic1 state))
                                    (m/* (.y av) v3))
                              v2 (m/+ (.ic2 state)
                                    (m/* (.y av) (.ic1 state))
                                    (m/* (.z av) v3))]
                          (SampleAndState. (m/+ (m/* sample (.x mv))
                                              (m/* v1 (.y mv))
                                              (m/* v2 (.z mv)))
                                           (SVF3State. (m/- (m/* 2.0 v1) (.ic1 state))
                                                       (m/- (m/* 2.0 v2) (.ic2 state))))))
                      ([] (SVF3State. 0.0 0.0)))))))


(def ^{:doc "Sorted sequence of the names (keywords) of all effects which can be created with [[effect]]."}
  effects-list (sort (keys (methods effect))))

;; ## File operations

(defn save-signal
  "Saves a signal to a file as 16 bit signed samples, big endian.

  Values are clipped to `[-1, 1]` and mapped linearly onto `[-32768, 32767]`, so `0.0` reads back as about `1.5e-5` (see [[load-signal]]). Use Audacity or SoX to convert the file to or from an audio format.

  Parameters:

  - `sig` (sequence of numbers): Signal to save.
  - `filename` (string or `java.io.File`): Path of the file. Missing parent directories are created, an existing file is overwritten.

  Returns the signal as a `double-array` (not clipped)."
  [sig filename]
  (make-parents filename)
  (let [^java.io.DataOutputStream out (java.io.DataOutputStream. (output-stream filename))
        s (m/seq->double-array sig)]
    (try
      (dotimes [i (alength s)]
        (.writeShort out (short (m/cnorm (aget s i) -1.0 1.0 Short/MIN_VALUE Short/MAX_VALUE))))
      (.flush out)
      (finally (. out clojure.core/close)))
    s))

(defn load-signal
  "Reads a signal from a file of 16 bit signed samples, big endian.

  Samples are mapped linearly from `[-32768, 32767]` onto `[-1, 1]` (see [[save-signal]]). Use Audacity or SoX to convert an audio file to this format.

  Parameters:

  - `filename` (string or `java.io.File`): Path of the file.

  Returns the signal as a `double-array` of length `floor(bytes / 2)`, a trailing odd byte is ignored. Throws `java.io.FileNotFoundException` when the file does not exist."
  [filename]
  (let [^java.io.File f (file filename)
        len (m/long-quot (.length f) 2)
        ^java.io.DataInputStream in (java.io.DataInputStream. (input-stream filename))
        ^doubles buffer (double-array len)]
    (try
      (dotimes [i len]
        (aset buffer i (m/cnorm (.readShort in) Short/MIN_VALUE Short/MAX_VALUE -1.0 1.0)))
      (finally (. in clojure.core/close)))
    buffer))

;; convolution / correlation

(defn- non-empty-signals!
  [sig1 sig2]
  (when (or (empty? sig1) (empty? sig2))
    (throw (ex-info "Signals should not be empty." {:size1 (count sig1) :size2 (count sig2)}))))

(defn- convolution-slice
  "Part of the full convolution `full` for `mode`, as `scipy.signal.convolve`."
  [mode ^doubles full ^long size1 ^long size2]
  (let [len (alength full)
        [^long start ^long n] (case mode
                                (:full :all) [0 len]
                                (:same :first) [(m/long-quot (m/long-sub len size1) 2) size1]
                                :valid (let [s (m/long-inc (m/long-sub (m/long-max size1 size2) (m/long-min size1 size2)))]
                                         [(m/long-quot (m/long-sub len s) 2) s])
                                (throw (ex-info "Unknown convolution mode." {:mode mode :allowed [:full :same :valid]})))]
    (vec (java.util.Arrays/copyOfRange full start (m/long-add start n)))))

(defn convolve
  "Computes the direct convolution of two signals, zero padded outside.

  Needs time proportional to `N1 N2`, for long signals see [[fft-convolve]].

  Parameters:

  - `sig1`, `sig2` (sequences of numbers): Signals of lengths `N1` and `N2`, not empty.
  - `mode` (keyword, optional): The part of the result to return, as in `scipy.signal.convolve`:
    - `:full` (default, alias `:all`) - the full convolution, length `N1 + N2 - 1`
    - `:same` (alias `:first`) - length `N1`, centred with respect to `:full`
    - `:valid` - length `|N1 - N2| + 1`, values computed without zero padding

  Returns a vector. Throws `ex-info` for an empty signal or an unknown mode.

  See also [[fft-convolve]], [[correlate]], [[filter-signal]]."
  ([sig1 sig2] (convolve sig1 sig2 :full))
  ([sig1 sig2 mode]
   (non-empty-signals! sig1 sig2)
   (let [a1 (m/seq->double-array sig1)
         a2 (m/seq->double-array sig2)]
     (convolution-slice mode (MathArrays/convolve a1 a2) (alength a1) (alength a2)))))

(defn fft-convolve
  "Computes the convolution of two signals using FFT (zero padded to a power of two).

  The result is the same as of [[convolve]] up to rounding errors, in time proportional to `N log N`, `N = N1 + N2 - 1`, so it is faster for long signals.

  Parameters:

  - `sig1`, `sig2` (sequences of numbers): Signals, not empty.
  - `mode` (keyword, optional): The part of the result to return, `:full` (default), `:same` or `:valid`, see [[convolve]].

  Returns a vector. Throws `ex-info` for an empty signal or an unknown mode.

  See also [[convolve]], [[fft-correlate]]."
  ([sig1 sig2] (fft-convolve sig1 sig2 :full))
  ([sig1 sig2 mode]
   (non-empty-signals! sig1 sig2)
   (let [s1 (count sig1)
         s2 (count sig2)
         tsize (m/long-dec (m/long-add s1 s2))
         size (m/round-up-pow2 tsize)
         fft1 (trans/fft (pad/zero (m/seq->double-array sig1) size :right))
         fft2 (trans/fft (pad/zero (m/seq->double-array sig2) size :right))
         full (trans/ifft (with-meta (mapv cplx/mult fft1 fft2) (meta fft1)))]
     (convolution-slice mode (m/seq->double-array (subvec full 0 tsize)) s1 s2))))

(defn correlate
  "Computes the direct cross-correlation of two real signals, zero padded outside.

  It is the convolution of `sig1` with the reversed `sig2`, as `scipy.signal.correlate`.

  Parameters:

  - `sig1`, `sig2` (sequences of numbers): Real signals, not empty.
  - `mode` (keyword, optional): The part of the result to return, `:full` (default), `:same` or `:valid`, see [[convolve]].

  Returns a vector. Throws `ex-info` for an empty signal or an unknown mode.

  See also [[fft-correlate]], [[convolve]]."
  ([sig1 sig2] (correlate sig1 sig2 :full))
  ([sig1 sig2 mode]
   (convolve sig1 (reverse sig2) mode)))

(defn fft-correlate
  "Computes the cross-correlation of two real signals using FFT, the result of [[correlate]] up to rounding errors.

  Parameters:

  - `sig1`, `sig2` (sequences of numbers): Real signals, not empty.
  - `mode` (keyword, optional): The part of the result to return, `:full` (default), `:same` or `:valid`, see [[convolve]].

  Returns a vector. Throws `ex-info` for an empty signal or an unknown mode.

  See also [[correlate]], [[fft-convolve]]."
  ([sig1 sig2] (fft-correlate sig1 sig2 :full))
  ([sig1 sig2 mode]
   (fft-convolve sig1 (reverse sig2) mode)))

;; filtering

(defn filter-signal
  "Applies FIR coefficients or an IIR filter to a signal.

  Parameters:

  - `FIR-or-IIR`: One of:
    - FIR coefficients (a sequence of numbers, see `fastmath.signal.fir`) - centred convolution, mode `:same` of [[convolve]] (as `numpy.convolve(xs, fir, 'same')` for `N >= M`, not a causal `lfilter`); computed directly for short signals and with FFT for long ones.
    - an IIR filter (see `fastmath.signal.iir`) - causal filtering, sample by sample. The filter is mutable and keeps its state between calls, reset it with [[reset-IIR!]] before filtering another signal.
  - `xs` (sequence of numbers): Signal to filter.

  Returns the filtered signal as a vector, of the same length as `xs`.

  See also [[filter-filter-signal]] (zero-phase), [[filter-signal-1]], [[convolve]]."
  [FIR-or-IIR xs]
  (if (instance? Cascade FIR-or-IIR)
    (let [^Cascade iir FIR-or-IIR]
      (mapv (fn [^double x] (.filter iir x)) xs))
    (let [cf (if (m/< (m/+ (count xs) (count FIR-or-IIR)) 2048)
               convolve fft-convolve)]
      (cf xs FIR-or-IIR :same))))

(defn filter-signal-1
  "Applies an IIR filter to a single sample or creates the function doing it.

  Parameters:

  - `IIR`: Filter from `fastmath.signal.iir`.
  - `x` (double, optional): Sample to filter.

  With `x` returns the filtered sample as a double. Without `x` returns a function of one sample which does the same.

  The filter is mutable: every call advances its internal state, use [[reset-IIR!]] to start again.

  See also [[filter-signal]], [[reset-IIR!]]."
  ([^Cascade IIR]
   (fn ^double [^double x] (.filter IIR x)))
  (^double [^Cascade IIR ^double x]
   (.filter IIR x)))

(defn reset-IIR!
  "Resets the internal state of an IIR filter.

  The filter keeps its state between calls of [[filter-signal]] and [[filter-signal-1]], reset it before filtering another signal.

  Parameters:

  - `IIR`: Filter from `fastmath.signal.iir`.

  Returns the same (mutated) filter."
  [^Cascade IIR]
  (.reset IIR)
  IIR)

(defn- maybe-reset-filter!
  [FIR-or-IIR]
  (when (instance? Cascade FIR-or-IIR)
    (.reset ^Cascade FIR-or-IIR)))

(defn filter-filter-signal
  "Applies a filter twice, forward and backward, to cancel the phase shift (zero-phase filtering).

  Parameters:

  - `FIR-or-IIR`: FIR coefficients or an IIR filter, see [[filter-signal]].
  - `xs` (sequence of numbers): Signal to filter.

  An IIR filter is reset between the passes but not before the forward pass. There is no edge padding, unlike `scipy.signal.filtfilt`. The magnitude response is squared.

  Returns the filtered signal as a vector, of the same length as `xs`.

  See also [[filter-signal]], [[zero-phase]]."
  [FIR-or-IIR xs]
  (let [forward (filter-signal FIR-or-IIR xs)]
    (maybe-reset-filter! FIR-or-IIR)
    (vec (rseq (filter-signal FIR-or-IIR (rseq forward))))))

;; signal smoothing

(defn- perform-convolution
  "Correlation of `signal` with odd length `coeffs` centred at sample, zero padding."
  [^doubles coeffs ^long fc signal]
  (let [full (convolve signal (reverse coeffs))]
    (subvec full fc (m/long-add fc (count signal)))))

(defn savgol-filter
  "Creates Savitzky-Golay smoothing (or differentiating) filter.

  Arguments:

  * length - odd length of the kernel (default: 5)
  * order - polynomial order, less than `length` (default: 2)
  * derivative - order of the derivative, `0 <= derivative <= order` (default: 0)
  * delta - spacing of the samples, used for derivatives (default: 1.0)

  Boundary rule is to pad with zeros (as `scipy.signal.savgol_filter` with `mode='constant'`).

  Returns filtering function which accepts collection of numbers and returns filtered signal (a vector). Throws `ex-info` for invalid arguments."
  ([] (savgol-filter 5))
  ([^long length] (savgol-filter length 2))
  ([^long length ^long order] (savgol-filter length order 0))
  ([^long length ^long order ^long derivative] (savgol-filter length order derivative 1.0))
  ([^long length ^long order ^long derivative ^double delta]
   (when-not (and (m/pos? length) (m/odd? length))
     (throw (ex-info "Length should be positive and odd." {:length length})))
   (when-not (m/< -1 order length)
     (throw (ex-info "Order should be non-negative and less than length." {:length length :order order})))
   (when-not (m/<= 0 derivative order)
     (throw (ex-info "Derivative should be between 0 and order." {:order order :derivative derivative})))
   (when-not (m/pos? delta)
     (throw (ex-info "Delta should be positive." {:delta delta})))
   (let [fc (m/long-quot (m/long-dec length) 2)
         factor (m// (m/factorial derivative) (m/pow delta derivative))
         coeffs (-> (for [v (range (m/- fc) (m/inc fc))]
                      (map #(m/pow v %) (range (m/inc order))))
                    (m/seq->double-double-array)
                    (Array2DRowRealMatrix.)
                    (SingularValueDecomposition.)
                    (.getSolver)
                    (.getInverse)
                    (.getRow derivative))
         coeffs (double-array (map (fn [^double c] (m/* factor c)) coeffs))]
     (fn [signal]
       (perform-convolution coeffs fc signal)))))

(defn moving-average-filter
  "Creates a moving average filter: every sample is replaced by the mean of `length` samples centred at it, the signal is zero padded at both ends.

  Parameters:

  - `length` (odd positive integer, optional): Size of the window. Default `5`.

  Returns a function which accepts a sequence of numbers and returns the smoothed signal as a vector of the same length. Throws `ex-info` for an even or non-positive `length`.

  This is [[savgol-filter]] of order 0. See also [[kernel-smoothing-filter]]."
  ([] (moving-average-filter 5))
  ([^long length] (savgol-filter length 0 0)))

(defn kernel-smoothing-filter
  "Creates Nadaraya-Watson kernel-weighted average.

  Arguments:

  * kernel - function of two numbers `(kernel x y)`, e.g. `(fastmath.kernel/kernel :gaussian)`
  * length - odd length of the kernel (default: 5)
  * step - distance between consecutive samples (default: 1.0)

  Weights `kernel(0, i step)`, `i = -(length-1)/2, ..., (length-1)/2`, are normalized to sum 1. Boundary rule is to pad with zeros (weights are not renormalized at edges).

  Returns filtering function which returns a vector. See also [[savgol-filter]]."
  ([kernel] (kernel-smoothing-filter kernel 5))
  ([kernel ^long length] (kernel-smoothing-filter kernel length 1.0))
  ([kernel ^long length ^double step]
   (when-not (and (m/pos? length) (m/odd? length))
     (throw (ex-info "Length should be positive and odd." {:length length})))
   (let [fc (m/long-quot (m/long-dec length) 2)
         coeffs (map (fn [^long v]
                       (kernel 0.0 (m/* step v))) (range (m/- fc) (m/inc fc)))
         ^double sum (reduce m/+ coeffs)
         coeffs (double-array (map (fn [^double v]
                                     (m// v sum)) coeffs))]
     (partial perform-convolution coeffs fc))))

;; Spectral core.
;; Coefficients are scaled back to the unnormalized DFT (`:norm :backward`); bins of a one-sided spectrum
;; which stand for two coefficients (all except DC and Nyquist) are doubled, as in `scipy.signal.periodogram`.

(defn- fft-layout
  "Layout of the coefficients `txs`: `{:one-sided? :n :scale}`.

  Read from the `:fastmath.transform/fft` metadata, `options` (`:kind`, `:n`, `:nyquist?`, `:norm`) override it. `:scale` converts coefficients to the unnormalized DFT."
  [txs options]
  (let [given (into {} (filter (comp some? val)) (select-keys options [:kind :n :nyquist? :norm]))
        {:keys [kind n nyquist? norm]
         :or {kind :real norm :backward}} (merge (::trans/fft (meta txs)) given)
        cnt (count txs)
        one-sided? (case kind
                     :real true
                     :complex false
                     (throw (ex-info "Unknown spectrum kind." {:kind kind :allowed [:real :complex]})))
        n (long (or n (cond
                        (not one-sided?) cnt
                        (nil? nyquist?) (throw (ex-info "Signal length of a one-sided spectrum is unknown: provide fft metadata, :n or :nyquist?." {:count cnt}))
                        nyquist? (m/long-mult 2 (m/long-dec cnt))
                        :else (m/long-dec (m/long-mult 2 cnt)))))]
    (when (m/zero? cnt)
      (throw (ex-info "Spectrum is empty." {})))
    (when-not (m/== cnt (if one-sided? (m/long-inc (m/long-quot n 2)) n))
      (throw (ex-info "Signal length :n does not match the number of coefficients." {:n n :count cnt :kind kind})))
    {:one-sided? one-sided?
     :n n
     :scale (case norm
              :backward 1.0
              :ortho (m/sqrt n)
              :forward (double n)
              (throw (ex-info "Unknown norm." {:norm norm :allowed [:backward :ortho :forward]})))}))

(defn- one-sided-weight
  "2 for bins of a one-sided spectrum which stand for two coefficients, 1 otherwise."
  ^double [one-sided? ^long n ^long k]
  (if (and one-sided?
           (m/pos? k)
           (not (and (m/even? n) (m/== k (m/long-quot n 2)))))
    2.0 1.0))

(defn- magnitudes
  [txs {:keys [one-sided? ^long n ^double scale]}]
  (vec (map-indexed (fn [^long k z]
                      (m/* (one-sided-weight one-sided? n k) scale (v/mag z))) txs)))

(defn- energies
  [txs {:keys [one-sided? ^long n ^double scale]}]
  (vec (map-indexed (fn [^long k z]
                      (m/* (one-sided-weight one-sided? n k) (m/sq (m/* scale (v/mag z))))) txs)))

(defn- scale-values
  [xs ^double factor]
  (mapv (fn [^double x] (m/* factor x)) xs))

(defn- spectral-values
  "Values of `method` for the coefficients `txs` of a segment windowed by a window with sum `wsum` and sum of squares `wsum2` (both `n` without window)."
  [method txs layout fs wsum wsum2]
  (let [fs (double fs)
        wsum (double wsum)
        wsum2 (double wsum2)]
    (case method
      :magnitude (magnitudes txs layout)
      :amplitude (scale-values (magnitudes txs layout) (m// wsum))
      :energy (energies txs layout)
      :power (scale-values (energies txs layout) (m// (m/sq wsum)))
      :psd (scale-values (energies txs layout) (m// (m/* fs wsum2)))
      :asd (mapv m/sqrt (scale-values (energies txs layout) (m// (m/* fs wsum2))))
      :phase (mapv v/heading txs)
      :real (mapv first txs)
      :imag (mapv second txs)
      :complex (vec txs)
      (throw (ex-info "Unknown spectrum method." {:method method
                                                  :allowed [:magnitude :amplitude :energy :power :psd :asd
                                                            :phase :real :imag :complex]})))))

(defn- ->db
  "Decibels: `10 log10` for power quantities, `20 log10` for amplitude quantities (0 gives `-Infinity`), other methods unchanged."
  [method xs]
  (case method
    (:energy :power :psd) (mapv (fn [^double x] (m/* 10.0 (m/log10 x))) xs)
    (:magnitude :amplitude :asd) (mapv (fn [^double x] (m/* 20.0 (m/log10 x))) xs)
    xs))

(defn fft-energy
  "Returns the energy spectrum `|X_k|^2` of the coefficients `txs` (from `fastmath.transform/fft`).

  For a one-sided spectrum the bins standing for two coefficients (all but DC and, for even `N`, Nyquist) are doubled, so the sum equals the sum over the two-sided spectrum (Parseval).

  Options (default: the `:fastmath.transform/fft` metadata of `txs`):

  * `:kind` - `:real` for a one-sided spectrum (`N/2+1` coefficients of a real signal), `:complex` for a two-sided spectrum (`N` coefficients)
  * `:n` - signal length `N`
  * `:nyquist?` - without `:n` and metadata, `true` when `N` is even (`N = 2(count-1)`), `false` for odd `N` (`N = 2 count - 1`); a one-sided spectrum without metadata, `:n` and `:nyquist?` throws `ex-info`
  * `:norm` - normalization used by `fft`: `:backward` (default), `:ortho`, `:forward`; values are computed for the unnormalized DFT

  Returns a vector of doubles. Throws `ex-info` when `:n` does not match the number of coefficients."
  ([txs] (fft-energy txs nil))
  ([txs options] (energies txs (fft-layout txs options))))

(defn fft-magnitude
  "Returns the magnitude spectrum `|X_k|` of the coefficients `txs` (from `fastmath.transform/fft`).

  For a one-sided spectrum the bins standing for two coefficients (all but DC and, for even `N`, Nyquist) are doubled, so a sine of amplitude `A` at bin `k` gives `A N`.

  Options as in [[fft-energy]]. Returns a vector of doubles."
  ([txs] (fft-magnitude txs nil))
  ([txs options] (magnitudes txs (fft-layout txs options))))

(defn fft-amplitude
  "Returns the amplitude spectrum, [[fft-magnitude]] divided by the signal length `N`: a one-sided spectrum gives the peak amplitude of each sinusoidal component (and the mean at DC).

  Options as in [[fft-energy]]. Returns a vector of doubles."
  ([txs] (fft-amplitude txs nil))
  ([txs options] (let [layout (fft-layout txs options)]
                   (scale-values (magnitudes txs layout) (m// (double (:n layout)))))))

(defn fft-power
  "Returns the power spectrum, [[fft-energy]] divided by `N^2`: the mean square of each component, sums to the mean square of the signal (`scipy.signal.periodogram` with `scaling='spectrum'`).

  Options as in [[fft-energy]]. Returns a vector of doubles."
  ([txs] (fft-power txs nil))
  ([txs options] (let [layout (fft-layout txs options)]
                   (scale-values (energies txs layout) (m// (m/sq (:n layout)))))))

(defn fft-frequencies
  "Returns the frequencies of the `N` DFT bins for sample rate `fs`: `k fs / N` for `k = 0, ..., N-1` (from 0 up to `fs - fs/N`; bins above `N/2` stand for negative frequencies `(k - N) fs / N`).

  For a one-sided spectrum take the first `N/2+1` values. Returns a vector of doubles."
  [^double fs ^long N]
  (let [step (m// fs N)]
    (mapv (fn [^long k] (m/* k step)) (range N))))

(defn- fft-times
  "Times of the segment centres, `(wlen/2 + i shift) / fs` (as `scipy.signal.spectrogram`)."
  [^double fs ^long shift ^long wlen]
  (map (fn [^long i] (m// (m/+ (m/* 0.5 wlen) (m/* i shift)) fs)) (range)))

(defn stft
  "Computes the Short-Time Fourier Transform (STFT) of a 1D signal.

  The signal is split into segments of the window length, shifted by `round((1 - overlap) wlen)` samples (at least 1); the last partial segment is dropped. Each segment is multiplied by the window and transformed.

  Input parameters:

  * `xs` - Input time-domain signal (sequence of doubles), not shorter than the window.
  * `options` - A map of configuration keys:
    * `:window` - A sequence of window coefficients (e.g., from `fastmath.kernel/window`). Default: Gaussian window of length `min(N, 256)` (256 is the scipy default segment length).
    * `:overlap` - The fraction of overlap between adjacent segments, in `[0, 1)`. Default is `0.5`.
    * `:fs` - The sampling frequency. Default is `1.0`.
    * `:method` - Values for each segment, as in [[spectrum]] with the window taken into account: `:power` (energy / (sum of window)^2, `scipy.signal.spectrogram` with `scaling='spectrum'`), `:psd` (energy / (fs sum of squared window), `scaling='density'`), `:asd` (`sqrt(psd)`), `:energy`, `:magnitude`, `:amplitude` (magnitude / sum of window), `:phase`, `:real`, `:imag`, `:complex` (coefficients). Default is `:power`. Like [[spectrum]], one-sided `:magnitude` and `:amplitude` are doubled on the bins standing for two coefficients (peak amplitude; scipy `mode='magnitude'` does not double), so there `amplitude^2 = 2 power` (peak versus mean square).
    * `:db?` - convert values to decibels: `10 log10` for `:energy`, `:power`, `:psd`, `20 log10` for `:magnitude`, `:amplitude`, `:asd` (0 gives `-Infinity`); other methods unchanged. Default `false`.
    * `:spectrum`, `:norm` - passed to `fastmath.transform/fft` (one-sided by default).

  Returns a map containing:

  * `:N` - The number of samples in the input signal.
  * `:spectrum` - A vector of spectra (vectors), one for each segment.
  * `:freqs` - Frequencies of the bins.
  * `:times` - Times of the segment centres, `(wlen/2 + i shift) / fs`.

  Throws `ex-info` for an empty signal or window, a signal shorter than the window, `:overlap` outside `[0, 1)` or an unknown method."
  ([xs] (stft xs nil))
  ([xs {:keys [window overlap fs method db?]
        :or {overlap 0.5 fs 1.0 method :power}
        :as options}]
   (let [N (count xs)
         fs (double fs)
         _ (when (m/zero? N)
             (throw (ex-info "Signal is empty." {})))
         _ (when-not (and (number? overlap) (m/<= 0.0 (double overlap)) (m/< (double overlap) 1.0))
             (throw (ex-info "Overlap should be in [0, 1)." {:overlap overlap})))
         window (vec (or window (ker/window :gaussian (m/long-min N 256))))
         wlen (count window)
         _ (when (m/zero? wlen)
             (throw (ex-info "Window is empty." {})))
         _ (when (m/< N wlen)
             (throw (ex-info "Signal is shorter than the window." {:N N :window-length wlen})))
         shift (m/long-max 1 (m/round (m/* (m/- 1.0 (double overlap)) wlen)))
         wsum (double (reduce m/+ 0.0 window))
         wsum2 (double (reduce (fn [^double acc ^double w] (m/+ acc (m/* w w))) 0.0 window))
         fft-options (select-keys options [:spectrum :norm])
         spectra (mapv (fn [segment]
                         (let [txs (trans/fft (mapv (fn [^double w ^double x] (m/* w x)) window segment) fft-options)
                               s (spectral-values method txs (fft-layout txs nil) fs wsum wsum2)]
                           (if db? (->db method s) s)))
                       (partition wlen shift xs))]
     {:N N
      :spectrum spectra
      :freqs (vec (take (count (first spectra)) (fft-frequencies fs wlen)))
      :times (vec (take (count spectra) (fft-times fs shift wlen)))})))

(defn spectrum
  "Returns the frequency spectrum of a signal using FFT.

  Values follow `scipy.signal.periodogram` with a boxcar window and no detrending.

  Input parameters:

  * `xs` - Input signal. Can be a sequence of doubles (time-domain) or a sequence of `Vec2` coefficients from `fastmath.transform/fft` (frequency-domain).
  * `options` - A map of configuration keys:
    * `:method` - `:magnitude` (`|X_k|`), `:amplitude` (default, magnitude / N), `:energy` (`|X_k|^2`), `:power` (energy / N^2, `scaling='spectrum'`), `:psd` (power spectral density, energy / (N fs), `scaling='density'`), `:asd` (amplitude spectral density, `sqrt(psd)`), `:phase` (angle in radians), `:real`, `:imag` or `:complex` (coefficients). For a one-sided spectrum magnitude and energy of the bins standing for two coefficients (all but DC and, for even `N`, Nyquist) are doubled.
    * `:db?` - convert to decibels: `10 log10` for `:energy`, `:power`, `:psd`, `20 log10` for `:magnitude`, `:amplitude`, `:asd` (0 gives `-Infinity`); other methods unchanged. Default `false`.
    * `:fs` - The sampling frequency. Default is `1.0`.
    * `:domain` - `:time` (default, `xs` is transformed with `fastmath.transform/fft`, `options` are passed to it, e.g. `:spectrum :double-sided`, `:norm`) or `:frequency` (`xs` are coefficients).
    * `:kind`, `:n`, `:nyquist?`, `:norm` - layout of coefficients for `:domain :frequency` when they have no `fft` metadata, see [[fft-energy]].

  Returns a map containing:

  * `:N` - The number of samples of the signal.
  * `:spectrum` - A vector of spectral values.
  * `:freqs` - Frequencies of the bins (see [[fft-frequencies]]).

  Throws `ex-info` for an empty signal, an unknown `:method` or `:domain`."
  ([xs] (spectrum xs nil))
  ([xs {:keys [method db? fs domain]
        :or {method :amplitude db? false fs 1.0 domain :time}
        :as options}]
   (let [fs (double fs)
         [txs layout-options] (case domain
                                :time [(trans/fft xs options) nil]
                                :frequency [xs options]
                                (throw (ex-info "Unknown domain." {:domain domain :allowed [:time :frequency]})))
         layout (fft-layout txs layout-options)
         n (double (:n layout))
         sp (spectral-values method txs layout fs n n)
         sp (if db? (->db method sp) sp)]
     {:N (:n layout)
      :spectrum sp
      :freqs (vec (take (count sp) (fft-frequencies fs (:n layout))))})))

(defn- median-bias
  "Bias of the median of `n` chi-squared(2) estimates, as `scipy.signal._spectral_py._median_bias`."
  ^double [^long n]
  (loop [i 1
         acc 1.0]
    (if (m/> i (m/long-quot (m/long-dec n) 2))
      acc
      (let [i2 (m/* 2.0 i)]
        (recur (m/long-inc i) (m/+ acc (m/- (m// (m/inc i2)) (m// i2))))))))

(defn periodogram
  "Estimate the spectral density of a signal with Welch's method.

  The signal is split into windowed segments (see [[stft]]), spectral values are calculated for each segment and averaged. It equals `scipy.signal.welch` with `detrend=False` and `noverlap = wlen - round((1 - overlap) wlen)`.

  Input parameters:

  * `xs` - Input time-domain signal (sequence of doubles), not shorter than the window.
  * `options` - A map of configuration keys:
    * `:method` - values calculated for segments, as in [[stft]]: `:psd` (default, `scaling='density'`), `:power` (`scaling='spectrum'`), `:energy`, `:magnitude`, `:amplitude`, `:asd`, `:phase`, `:real`, `:imag`.
    * `:average` - aggregation of segments: `:mean` (default), `:median` or `:umedian` (median divided by the bias of the median for the number of segments, `average='median'` in scipy).
    * `:db?` - convert the averaged values to decibels, as in [[stft]].
    * `:window`, `:overlap`, `:fs`, `:spectrum`, `:norm` - as in [[stft]].

  Returns a map containing:

  * `:freqs` - Frequencies of the bins.
  * `:spectrum` - Averaged spectral values.

  Throws `ex-info` for a signal shorter than the window, an unknown `:method` (also `:complex`) or `:average`."
  ([xs] (periodogram xs nil))
  ([xs {:keys [method average db?] :or {method :psd average :mean} :as options}]
   (when (= :complex method)
     (throw (ex-info "Complex coefficients can't be averaged." {:method method})))
   (let [{:keys [freqs spectrum]} (stft xs (assoc options :method method :db? false))
         columns (apply mapv vector spectrum)
         avg (case average
               :mean (mapv stats/mean columns)
               :median (mapv (fn [col] (stats/median col :r7)) columns)
               :umedian (let [bias (median-bias (count spectrum))]
                          (mapv (fn [col] (m// (stats/median col :r7) bias)) columns))
               (throw (ex-info "Unknown average." {:average average :allowed [:mean :median :umedian]})))]
     {:freqs freqs
      :spectrum (if db? (->db method avg) avg)})))

;; Waveforms

(defn waveform
  "Create a waveform function (oscillator) for signal generation.

  This function returns a stateless oscillator that maps time (as a double) to a signal value. It supports a variety of wave shapes, including standard geometric oscillators, 'analog-style' approximations, and additive band-limited synthesis to reduce aliasing. Parameters such as frequency, amplitude, and phase can be provided as constant numbers or as functions of time, enabling complex modulations (FM/AM).

  Input parameters:
  * `wave-type` - A keyword selecting the waveform shape:
      * Basic: `:sine`, `:square`, `:saw`, `:triangle`.
      * Analog-style: `:analog-sine`, `:analog-saw`, `:analog-triangle`.
      * Band-limited: `:band-limited-square`, `:band-limited-saw`, `:band-limited-triangle`.
  * `options` - A map of configuration keys:
      * `:f` - Frequency in Hz (number, default `1.0`), or a phase function `(fn [t])` returning the phase in cycles (used instead of `f t`), for frequency modulation as in [[chirp]].
      * `:phase` - Phase offset in cycles [0.0 to 1.0]. Can be a number or a function `(fn [t])`. Default: `0.0`.
      * `:amplitude` - Peak amplitude. Can be a number or a function `(fn [t])` for amplitude modulation. Default: `1.0`.
      * `:duty` - Duty cycle for `:square` waves, ranging from 0.0 to 1.0 or `(fn [t])` for modulation. Default: `0.5`.
      * `:up?` - Boolean for `:saw` and `:analog-saw` waves. Set to `true` (default) for a rising ramp, `false` for a falling ramp.
      * `:bands` - Number of harmonics for `:band-limited-*` types. Increasing this value improves shape accuracy but increases computation. Default: `5`.

  The waves are periodic also for negative time; shapes and their relation to `scipy.signal.square` / `sawtooth` are described in `fastmath.signal.waveform`.

  Returns a function `(fn [t])` that accepts time `t` as a double and returns the signal value as a double. Throws `ex-info` for an unknown `wave-type`."
  ([wave-type] (waveform wave-type nil))
  ([wave-type options]
   (case wave-type
     :sine (wv/->sine options)
     :analog-sine (wv/->analog-sine options)
     :square (wv/->square options)
     :band-limited-square (wv/->band-limited-square options)
     :saw (wv/->saw options)
     :analog-saw (wv/->analog-saw options)
     :band-limited-saw (wv/->band-limited-saw options)
     :triangle (wv/->triangle options)
     :analog-triangle (wv/->analog-triangle options)
     :band-limited-triangle (wv/->band-limited-triangle options)
     (throw (ex-info "Unknown wave type." {:wave-type wave-type
                                           :allowed [:sine :analog-sine :square :band-limited-square :saw :analog-saw
                                                     :band-limited-saw :triangle :analog-triangle :band-limited-triangle]})))))

(defn chirp
  "Create a frequency-swept sine wave (chirp) oscillator.

  Generates a signal where the instantaneous frequency changes over a specified time interval according to a defined modulation profile. Chirps are essential in radar, sonar, and system identification to analyze frequency-dependent responses and impulse behavior.

  Input parameters:
  * `chirp-type` - A keyword specifying the sweep profile (see `fastmath.signal.chirp`):
      * `:linear` - Frequency changes linearly with time.
      * `:quadratic-up` - Frequency changes quadratically, parabola vertex at $t=0$.
      * `:quadratic-down` - Frequency changes quadratically, parabola vertex at $t=time$.
      * `:logarithmic` - Frequency changes exponentially (geometric sweep), frequencies must be positive.
      * `:hyperbolic` - Frequency follows a hyperbolic curve, frequencies must be positive.
      * `:constant` - Constant frequency `:f0`.
  * `options` - A map of configuration keys:
      * `:f0` - Starting frequency in Hz at $t=0$ (default 0.0; 1e-10 for `:logarithmic` and `:hyperbolic`).
      * `:f1` - Ending frequency in Hz at $t=time$ (default 5.0).
      * `:time` - The time (in seconds) at which the frequency reaches `:f1` (default 1.0, must be positive).
      * `:phase0` - Phase offset of the sweep in cycles (default 0.0).
      * `:amplitude` - Peak amplitude. Can be a number or a modulation function `(fn [t])`.
      * `:phase` - Additional phase offset in cycles, as for [[waveform]].

  The value is `amplitude * sin(2 pi phase(t))`; it equals `scipy.signal.chirp(t, f0, time, f1, method, phi=-90 + 360 phase0)` (`:quadratic-up`/`-down` are `method='quadratic'` with `vertex_zero` true/false).

  Returns a stateless function `(fn [t])` that accepts time `t` as a double and returns the signal value as a double. Throws `ex-info` for an unknown `chirp-type` or invalid options."
  ([chirp-type] (chirp chirp-type nil))
  ([chirp-type options]
   (let [modulation (case chirp-type
                      :linear (chirp/->linear options)
                      :quadratic-up (chirp/->quadratic-up options)
                      :quadratic-down (chirp/->quadratic-down options)
                      :logarithmic (chirp/->logarithmic options)
                      :hyperbolic (chirp/->hyperbolic options)
                      :constant (chirp/->constant options)
                      (throw (ex-info "Unknown chirp type." {:chirp-type chirp-type
                                                             :allowed [:linear :quadratic-up :quadratic-down
                                                                       :logarithmic :hyperbolic :constant]})))]
     (waveform :sine (merge options {:f modulation})))))

(defn add-waveforms
  "Adds waveforms.

  Parameters:

  - `w`, `w1`, `w2`, ... (functions): Waveforms, functions `(fn [t])` of the time in seconds returning the amplitude, as created by [[waveform]] or [[chirp]].

  Returns the waveform which is the sum of the given ones at `t`. A single waveform is returned unchanged.

  See also [[gain-waveform]], [[sample-waveform]]."
  ([w] w)
  ([w1 w2] (fn [^double t] (m/+ (double (w1 t)) (double (w2 t)))))
  ([w1 w2 w3] (fn [^double t] (m/+ (double (w1 t)) (double (w2 t)) (double (w3 t)))))
  ([w1 w2 w3 w4] (fn [^double t] (m/+ (double (w1 t)) (double (w2 t)) (double (w3 t)) (double (w4 t)))))
  ([w1 w2 w3 w4 w5] (fn [^double t] (m/+ (double (w1 t)) (double (w2 t)) (double (w3 t))
                                        (double (w4 t)) (double (w5 t)))))
  ([w1 w2 w3 w4 w5 & r]
   (reduce add-waveforms (add-waveforms w1 w2 w3 w4 w5) r)))

(defn gain-waveform
  "Scales the amplitude of a waveform.

  Parameters:

  - `w` (function): Waveform, a function `(fn [t])` of the time in seconds returning the amplitude, see [[waveform]].
  - `gain` (double): Factor, a negative value inverts the waveform.

  Returns the waveform `(fn [t])` equal to `gain * (w t)`.

  See also [[add-waveforms]], [[sample-waveform]]."
  [w ^double gain]
  (fn ^double [^double t] (m/* gain (double (w t)))))

(defn sample-waveform
  "Discretize a continuous waveform function into a sequence of samples.

  Evaluates a time-dependent function at `t = i / fs` for every `i >= 0` with `t < time` (`ceil(fs time)` samples, as `numpy.arange(0, time, 1/fs)`; a product `fs time` within `1e-9` relative of an integer counts as that integer, so `fs 10, time 0.3` gives 3 samples).

  Input parameters:
  * `waveform-fn` - A function `(fn [t])` that accepts time as a double and returns the signal amplitude. Usually created via [[waveform]] or [[chirp]].
  * `fs` - Sampling frequency in Hz (samples per second), positive.
  * `time` - The total duration of the signal to generate in seconds. Defaults to `1.0`.

  Returns a vector of doubles. Throws `ex-info` when `fs` is not positive."
  ([waveform-fn ^double fs] (sample-waveform waveform-fn fs 1.0))
  ([waveform-fn ^double fs ^double time]
   (when-not (m/pos? fs)
     (throw (ex-info "Sampling frequency should be positive." {:fs fs})))
   (let [x (m/* fs time)
         r (m/round x)
         n (if (m/<= (m/abs (m/- x r)) (m/* 1.0e-9 (m/max 1.0 (m/abs x)))) r (long (m/ceil x)))]
     (mapv (fn [^long i] (double (waveform-fn (m// i fs))))
           (range (m/long-max 0 n))))))

;; padding

(defn pad
  "Pad a 1D signal to a specified length using various boundary conditions.

  Extends a signal to a target length `N`, which is a common preprocessing step for transforms that require power-of-two input sizes (like FFT or DWT) or to mitigate edge artifacts during filtering and convolution. The function supports multiple padding strategies to maintain signal continuity or satisfy specific boundary assumptions.

  Parameters:

  * `signal` - A sequence or array of real numbers representing the input signal.
  * `N` - The desired target length (must be greater than or equal to the current signal length). If not provided, defaults to the smallest power of two greater than or equal to the signal length.
  * `pad-method` - A keyword specifying the padding strategy (default is `:periodic`):
      * `:zero` - Appends zeros.
      * `:edge` - Repeats the last/first value.
      * `:linear` - Linear ramp between the edge and zero.
      * `:periodic` - Circular padding (repeats the signal).
      * `:symmetric` - Mirroring at the edge (repeats edge values).
      * `:antisymmetric` - Mirroring with sign inversion.
      * `:reflect` - Mirroring at the edge (does not repeat edge values).
      * `:antireflect` - Anti-mirroring at the edge.
  * `side` - A keyword specifying where to apply the padding: `:left`, `:right`, or `:both` (default).

  Returns a double array of length `N` containing the padded signal."
  ([signal] (pad signal nil))
  ([signal N] (pad signal N :periodic))
  ([signal N pad-method] (pad signal N pad-method :both))
  ([signal N pad-method side]
   (let [N (long (or N (m/round-up-pow2 (count signal))))]
     (if (m/> (count signal) N)
       (throw (ex-info "New length of the signal is lower than signal size."
                       {:N N :signal-length (count signal)}))
       (let [asignal (m/seq->double-array signal)]
         (case pad-method
           :zero (pad/zero asignal N side)
           :edge (pad/edge asignal N side)
           :linear (pad/linear asignal N side)
           :periodic (pad/periodic asignal N side)
           :symmetric (pad/symmetric asignal N side)
           :antisymmetric (pad/antisymmetric asignal N side)
           :reflect (pad/reflect asignal N side)
           :antireflect (pad/antireflect asignal N side)
           (throw (ex-info "Unknown padding method" {:pad-method pad-method}))))))))

;;

(defn hilbert
  "Calculates the analytic signal `x + i H(x)` of a real signal with FFT, as `scipy.signal.hilbert`.

  `H(x)` is the Hilbert transform of `x`: the signal with every component shifted in phase by 90 degrees. The analytic signal is the base of [[amplitude-envelope]], [[instantaneous-phase]] and [[instantaneous-frequency]].

  Parameters:

  - `xs` (sequence of numbers): Real signal, not empty.

  Returns a vector of complex numbers (`Vec2`: real and imaginary part) of the same length, the real part is `xs`. Throws `ex-info` for an empty signal."
  [xs]
  (let [t (trans/fft xs {:spectrum :double-sided})
        cnt (count t)
        hN (m// cnt 2)
        nyquist? (m/even? cnt)]
    (-> (map-indexed (fn [^long id z]
                       (cond
                         (m/zero? id) z
                         (and nyquist? (m/== id hN)) z ;; nyquist
                         (m/<= id hN) (cplx/scale z 2.0)
                         :else cplx/ZERO)) t)
        (trans/ifft {:kind :complex :real? false}))))

(defn amplitude-envelope
  "Calculates the amplitude envelope of a real signal: the magnitude of its analytic signal ([[hilbert]]).

  Parameters:

  - `xs` (sequence of numbers): Real signal, not empty.

  Returns a lazy sequence of doubles of the same length. Throws `ex-info` for an empty signal.

  See also [[instantaneous-phase]], [[instantaneous-frequency]]."
  [xs]
  (->> (hilbert xs) (map cplx/abs)))

(defn instantaneous-phase
  "Calculates the instantaneous phase of a real signal in radians: the argument of its analytic signal ([[hilbert]]).

  Parameters:

  - `xs` (sequence of numbers): Real signal, not empty.
  - `unwrap?` (boolean, optional): When `true` (default), jumps larger than `pi` are removed by adding multiples of `2 pi`, so the phase is continuous. With `false` the values are in `[-pi, pi]`.

  Returns a lazy sequence of doubles of the same length. Throws `ex-info` for an empty signal.

  See also [[instantaneous-frequency]], [[amplitude-envelope]]."
  ([xs] (instantaneous-phase xs true))
  ([xs unwrap?]
   (let [args (->> (hilbert xs) (map cplx/arg))]
     (if unwrap? (v/unwrap args m/TWO_PI) args))))

(defn instantaneous-frequency
  "Calculates the instantaneous frequency of a real signal: the differences of the unwrapped [[instantaneous-phase]] times `fs / (2 pi)`.

  Parameters:

  - `xs` (sequence of numbers): Real signal of `N` samples, not empty.
  - `fs` (double, optional): Sampling frequency. Default `1.0`, the frequency is then in cycles per sample; with `fs` in Hz it is in Hz.

  Returns a vector of `N - 1` doubles. Throws `ex-info` for an empty signal.

  See also [[amplitude-envelope]], [[hilbert]]."
  ([xs] (instantaneous-frequency xs 1.0))
  ([xs ^double fs]
   (-> (instantaneous-phase xs)
       (v/differences)
       (v/div m/TWO_PI)
       (v/mult fs))))

;;

(defn zero-phase
  "Sets the phase of every component of the spectrum of a signal to zero and keeps the magnitudes.

  The signal is transformed with FFT, every coefficient is replaced by its magnitude and the result is transformed back. The result is real and symmetric (the sample `k` equals the sample `N - k`), so its energy is concentrated around the first sample.

  Parameters:

  - `xs` (sequence of numbers): Real signal, not empty.

  Returns a vector of the same length. Throws `ex-info` for an empty signal.

  See also [[filter-filter-signal]] (zero-phase filtering), [[hilbert]]."
  [xs]
  (let [X (trans/fft xs)]
    (trans/ifft (with-meta (mapv (comp cplx/complex cplx/abs) X) (meta X)))))

;;

(defn resample
  "Resamples a signal to `nsize` samples using FFT interpolation, as `scipy.signal.resample`.

  The signal is assumed to be periodic, so for a signal which is not, the edges can show ringing.

  Parameters:

  - `xs` (sequence of numbers): Real signal, not empty.
  - `nsize` (integer): Number of samples of the result, at least `1`. For `nsize` equal to the length of `xs` the signal is returned as doubles.

  Returns a vector of doubles. Throws `ex-info` for an empty signal or `nsize` less than `1`.

  See also [[pad]], [[spectrum]]."
  [xs ^long nsize]
  (let [size (count xs)]
    (when (m/zero? size)
      (throw (ex-info "Signal should not be empty." {})))
    (when (m/< nsize 1)
      (throw (ex-info "New size should be positive." {:nsize nsize})))
    (if (m/== nsize size)
      (mapv double xs)
      (let [ratio (m// (double nsize) size)
            m (m/long-min size nsize)
            m2 (m/long-inc (m/long-quot m 2))
            X (as-> (vec (take m2 (trans/fft xs))) X
                (if (m/odd? m) X (update X (m/long-quot m 2) cplx/scale (if (m/< nsize size) 2.0 0.5)))
                (if (m/< nsize size) X
                    (concat X (repeat (m/long-sub (m/long-quot nsize 2) (m/long-quot size 2)) cplx/ZERO)))
                (map #(cplx/scale % ratio) X))]
        (trans/ifft X {:kind :real :n nsize})))))
