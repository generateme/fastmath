(ns fastmath.signal
  "Signal processing (effect) and generation (oscillators).

  Singal is any sequence with double values.

  ## Signal processing

  To process signal use [[apply-effects]] or [[apply-effects-raw]] (operates on `double-array` only) function.

  Effect is signal filter, created by [[effect]] multimethod. Effects can be composed with [[compose-effects]]. Effect can be treated as function and can be called for given sample.

  Each effect has it's own parametrization which should be passed during creation.
  
  List of all available effects is under [[effects-list]] value.

  ### Effects parametrization

  Each effect has its own parameters.

  #### :simple-lowpass, :simple-highpass

  * `:rate` - sample rate (default 44100.0)
  * `:cutoff` - cutoff frequency (default 2000.0)

  #### :biquad-eq

  Biquad equalizer
  
  * `:fc` - center frequency
  * `:gain` - gain
  * `:bw` - bandwidth (default: 1.0)
  * `:fs` - sampling rate (defatult: 44100.0)

  #### :biquad-hs, :biquad-ls

  Biquad highpass and lowpass shelf filters
  
  * `:fc` - center frequency
  * `:gain` - gain
  * `:slope` - shelf slope (default 1.5)
  * `:fs` - sampling rate (default 44100.0)

  #### :biquad-lp, :biquad-hp, :biquad-bp

  Biquad lowpass, highpass and bandpass filters
  
  * `:fc` - cutoff/center frequency
  * `:bw` - bandwidth (default 1.0)
  * `:fs` - sampling rate (default 44100.0)

  #### :dj-eq

  * `:high` - high frequency gain (10000Hz)
  * `:mid` - mid frequency gain (1000Hz)
  * `:low` - low frequency gain (100Hz)
  * `:shelf-slope` - shelf slope for high frequency (default 1.5)
  * `:peak-bw` - peak bandwidth for mid and low frequencies (default 1.0)
  * `:rate` - sampling rate (default 44100.0)

  #### :phaser-allpass

  * `:delay` - delay factor (default: 0.5)

  #### :divider

  * `:denom` (long, default 2.0)

  #### :fm

  Modulate and demodulate signal using frequency

   * `:quant` - quantization value (0.0 - if no quantization, default 10)
   * `:omega` - carrier factor (default 0.014)
   * `:phase` - deviation factor (default 0.00822)

  #### :bandwidth-limit

  https://searchcode.com/file/18573523/cmt/src/lofi.cpp#

  * `:rate` - sample rate (default 44100.0)
  * `:freq` - cutoff frequency (default 1000.0)

  #### :distort

  * `:factor` - distortion factor (default 1.0)

  #### :foverdrive

  Fast overdrive
  
  * `:drive` - drive (default 2.0)

  #### :decimator

  * `:bits` - bit depth (default 2)
  * `:fs` - decimator sample rate (default 4410.0)
  * `:rate` - input sample rate (default 44100.0)

  #### :basstreble

  * `:bass` - bass gain (default 1.0)
  * `:treble` - treble gain (default 1.0)
  * `:gain` - gain (default 0.0)
  * `:rate` - sample rate (default 44100.0)
  * `:slope` - slope for both (default 0.4)
  * `:bass-freq` - bass freq (default 250.0)
  * `:treble-freq` - treble freq (default 4000.0)

  #### :echo

  * `:delay` - delay time in seconds (default 0.5)
  * `:decay` - decay (amount echo in signal, default 0.5)
  * `:rate` - sample rate (default 44100.0)
  
  _Warning! Echo filter uses mutable array as a internal state, don't use the same filter in paraller processing._

  #### :vcf303

  * `:rate` - sample rate (default 44100.0)
  * `:trigger` - boolean, trigger some action (default `false`), set true when you reset filter every line
  * `:cutoff` - cutoff frequency (values 0-1, default 0.8)
  * `:resonance` - resonance (values 0-1, default 0.8)
  * `:env-mod` - envelope modulation (values 0-1, default 0.5)
  * `:decay` - decay (values 0-1, default 1.0)
  * `:gain` - gain output signal (default: 1.0)
  
  #### :slew-limit

  http://git.drobilla.net/cgit.cgi/omins.lv2.git/tree/src/slew_limiter.c

  * `:rate` - sample rate
  * `:maxrise` - maximum change for rising signal (in terms of 1/rate steps, default 500)
  * `:maxfall` - maximum change for falling singal (default 500)

  #### :mda-thru-zero

  * `:rate` - sample rate
  * `:speed` - effect rate
  * `:depth`
  * `:mix`
  * `:depth-mod`
  * `:feedback`
  
  _Warning: internal state is kept in doubles array._
  
  ## Oscillators

  [[oscillator]] creates function which generates signal value for given time.

  To sample generated wave to signal, call [[oscillator->signal]] with following parameters:

  * `f` - oscillator
  * `samplerate` - sample rate (samples per second)
  * `seconds` - duration

  To convert signal to oscillator (using interpolation) use [[signal->oscillator]] passing signal and duration.

  Add oscillators using [[oscillators-sum]].

  ## Smoothing filter

  Savitzky-Golay smoothing filter [[savgol-filter]].
  
  ## File operations

  You can [[save-signal]] or [[load-signal]]. Representation is 16 bit signed, big endian. Use Audacity or SoX to convert to/from audio files."
  (:require [fastmath.core :as m]
            [clojure.java.io :refer [file make-parents output-stream input-stream]]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [fastmath.interpolation.linear :as linear-interp])
  (:import [fastmath.vector Vec3]
           [clojure.lang IFn]
           [org.apache.commons.math3.linear Array2DRowRealMatrix SingularValueDecomposition]
           [org.apache.commons.math3.util MathArrays]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

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
  "Create `EffectsList` node from effect function and initial state."
  [effect-name f]
  (EffectsList. effect-name 0.0 f (f) nil))

(defn- compose-two-effects
  "Join two effect lists."
  [^EffectsList e1 ^EffectsList e2]
  (EffectsList. (.effect-name e1) (.sample e1) (.effect e1) (.state e1) (if (.next e1)
                                                                          (compose-two-effects (.next e1) e2)
                                                                          e2)))

(defn compose-effects
  "Compose effects."
  [^EffectsList e & es]
  (reduce compose-two-effects e es))

(defn reset-effects
  "Resets effects state to initial one."
  [^EffectsList e]
  (EffectsList. (.effect-name e) 0.0 (.effect e) ((.effect e)) (when (.next e)
                                                                 (reset-effects (.next e)))))

(defn single-pass
  "Process on sample using effects, returns `EffectsList` with result and new effect states."
  [^EffectsList e ^double sample]
  (if-not (.next e)
    (let [^SampleAndState r ((.effect e) sample (.state e))]
      (EffectsList. (.effect-name e) (.sample r) (.effect e) (.state r) nil))
    (let [^EffectsList prev (single-pass (.next e) sample)
          ^SampleAndState r ((.effect e) (.sample prev) (.state e))]
      (EffectsList. (.effect-name e) (.sample r) (.effect e) (.state r) prev))))

(defn apply-effects-raw
  "Apply effects to signal as `double-array`.

  If `reset` is positive, reinit state each `reset` number of samples.

  Returns new signal as `double-array`."
  (^doubles [^doubles in effects ^long reset]
   (let [len (alength in)
         ^doubles out (double-array len)]
     (loop [idx (int 0)
            effects-and-state effects]
       (when (< idx len)
         (let [sample (aget in idx)
               ^EffectsList res (single-pass effects-and-state sample)
               idx+ (inc idx)] 
           (aset out idx ^double (.sample res))
           (recur idx+
                  (if (and (pos? reset) (zero? ^long (mod idx+ reset)))
                    (reset-effects effects)
                    res)))))
     out))
  (^doubles [^doubles in effects] (apply-effects-raw in effects 0)))

(defn apply-effects
  "Apply effects to signal as any sequence.

  If `reset` is positive, reinit state each `reset` number of samples.

  Returns new signal."
  ([in effects ^long reset] (m/double-array->seq (apply-effects-raw (m/seq->double-array in) effects reset)))
  ([in effects] (apply-effects in effects 0)))

;; ## Helper functions

(defn db->linear
  "DB to Linear"
  ^double [^double x]
  (m/pow 10.0 (/ x 20.0)))

(defn linear->db
  "Linear to DB"
  ^double [^double x]
  (* 20.0 (m/log10 x)))

;; ## Effects / Filters

(defmulti effect
  "Create effect for given name (as keyword) and optional parameters.

  List of all possible effects is under [[effects-list]].

  Effect is a custom type which contains: name, sample (result of last call), effect function and current state.

  Effect can be considered as function: call with sample to effect with next state or call without parameter to obtain latest result. Effects are also composable with [[compose-effects]]."
  (fn [m & _] m))

;; ### Simple Low/High pass filters
(defn- calc-filter-alpha
  "Calculate alpha factor"
  ^double [^double rate ^double cutoff]
  (let [tinterval (/ rate)
        tau (/ (* cutoff m/TWO_PI))]
    (/ tinterval (+ tau tinterval))))

(defmethod effect :simple-lowpass
  ([m] (effect m {}))
  ([m {:keys [rate cutoff]
       :or {rate 44100.0 cutoff 2000.0}}]
   (let [alpha (calc-filter-alpha rate cutoff)] 
     (effect-node m (fn
                      ([^double sample ^double prev]
                       (let [s1 (* sample alpha)
                             s2 (- prev (* prev alpha))
                             nprev (+ s1 s2)]
                         (SampleAndState. nprev nprev)))
                      ([] 0.0))))))

(defmethod effect :simple-highpass
  ([m] (effect m {}))
  ([m conf]
   (let [lpfilter (effect :simple-lowpass conf)]
     (effect-node m (fn
                      ([^double sample lp]
                       (let [^EffectsList res (single-pass lp sample)]
                         (SampleAndState. (- sample (.sample res)) res)))
                      ([] lpfilter))))))

;; ### Biquad filters

;; Store biquad effect configuration in `BiquadConf` type
(deftype BiquadConf [^double b0 ^double b1 ^double b2 ^double a1 ^double a2])

(defn- biquad-eq-params
  "Calculate configuration for biquad equalizer
   fc - center frequency
   gain
   bw - bandwidth
   fs - sample rate"
  [^double fc ^double gain ^double bw ^double fs]
  (let [w (/ (* m/TWO_PI (m/constrain fc 1.0 (* 0.5 fs))) fs)
        cw (m/cos w)
        sw (m/sin w)
        J (m/pow 10.0 (* gain 0.025))
        g (-> bw
              (m/constrain 0.0001 4.0)
              (* m/LN2_2 w)
              (/ sw)
              (m/sinh)
              (* sw))
        a0r (/ (+ 1.0 (/ g J)))

        b0 (* a0r (+ 1.0 (* g J)))
        b1 (* a0r -2.0 cw)
        b2 (* a0r (- 1.0 (* g J)))
        a1 (- b1)
        a2 (* a0r (- (/ g J) 1.0))]
    (BiquadConf. b0 b1 b2 a1 a2)))

(defn- biquad-hs-params 
  "Calculate configuration for biquad high shelf
   fc - center frequency
   gain
   slope - shelf slope
   fs - sample rate"
  [^double fc ^double gain ^double slope ^double fs]
  (let [w (/ (* m/TWO_PI (m/constrain fc 1.0 (* 0.5 fs))) fs)
        cw (m/cos w)
        sw (m/sin w)
        A (m/pow 10.0 (* gain 0.025))
        iA (inc A)
        dA (dec A)
        b (m/sqrt (- (/ (inc (* A A)) (m/constrain slope 0.0001 1.0)) (* dA dA)))
        apc (* cw iA)
        amc (* cw dA)
        bs (* b sw)
        a0r (->> amc
                 (- iA)
                 (+ bs)
                 (/ 1.0))

        b0 (* a0r A (+ iA amc bs))
        b1 (* a0r A -2.0 (+ dA apc))
        b2 (* a0r A (- (+ iA amc) bs))
        a1 (* a0r -2.0 (- dA apc))
        a2 (* a0r (+ (dec (- A)) amc bs))]
    (BiquadConf. b0 b1 b2 a1 a2)))

(defn- biquad-ls-params 
  "Calculate configuration for biquad low shelf
   fc - center frequency
   gain
   slope - shelf slope
   fs - sample rate"
  [^double fc ^double gain ^double slope ^double fs]
  (let [w (/ (* m/TWO_PI (m/constrain fc 1.0 (* 0.5 fs))) fs)
        cw (m/cos w)
        sw (m/sin w)
        A (m/pow 10.0 (* gain 0.025))
        iA (inc A)
        dA (dec A)
        b (m/sqrt (- (/ (inc (* A A)) (m/constrain slope 0.0001 1.0)) (* dA dA)))
        apc (* cw iA)
        amc (* cw dA)
        bs (* b sw)
        a0r (->> amc
                 (+ iA)
                 (+ bs)
                 (/ 1.0))

        b0 (* a0r A (- (+ iA bs) amc))
        b1 (* a0r A 2.0 (- dA apc))
        b2 (* a0r A (- iA amc bs))
        a1 (* a0r 2.0 (+ dA apc))
        a2 (* a0r (+ bs (- (dec (- A)) amc)))]
    (BiquadConf. b0 b1 b2 a1 a2)))

(defn- biquad-lp-params
  "Calculate configuration for biquad low pass"
  [^double fc ^double bw ^double fs]
  (let [omega (* m/TWO_PI (/ fc fs))
        sn (m/sin omega)
        cs (m/cos omega)
        alpha (* sn (m/sinh (* m/LN2_2 bw (/ omega sn))))
        a0r (/ (inc alpha))
        cs- (- 1.0 cs)
        
        b0 (* a0r 0.5 cs-)
        b1 (* a0r cs-)
        b2 b0
        a1 (* a0r 2.0 cs)
        a2 (* a0r (dec alpha))]
    (BiquadConf. b0 b1 b2 a1 a2)))

(defn- biquad-hp-params
  "Calculate configuration for biquad high pass"
  [^double fc ^double bw ^double fs]
  (let [omega (* m/TWO_PI (/ fc fs))
        sn (m/sin omega)
        cs (m/cos omega)
        alpha (* sn (m/sinh (* m/LN2_2 bw (/ omega sn))))
        a0r (/ (inc alpha))
        cs+ (inc cs)
        
        b0 (* a0r 0.5 cs+)
        b1 (* a0r (- cs+))
        b2 b0
        a1 (* a0r 2.0 cs)
        a2 (* a0r (dec alpha))]
    (BiquadConf. b0 b1 b2 a1 a2)))

(defn- biquad-bp-params
  "Calculate configuration for biquad band pass"
  [^double fc ^double bw ^double fs]
  (let [omega (* m/TWO_PI (/ fc fs))
        sn (m/sin omega)
        cs (m/cos omega)
        alpha (if (zero? sn) 0.0 (* sn (m/sinh (* m/LN2_2 bw (/ omega sn)))))
        a0r (/ (inc alpha))
        
        b0 (* a0r alpha)
        b1 0.0
        b2 (* a0r (- alpha))
        a1 (* a0r 2.0 cs)
        a2 (* a0r (dec alpha))]
    (BiquadConf. b0 b1 b2 a1 a2)))

;; Store state in `StateBiquad` type.
(deftype StateBiquad [^double x2 ^double x1 ^double y2 ^double y1])

(defn- make-biquad-filter
  "Create biquad effect based on passed configuration"
  [m ^BiquadConf c]
  (effect-node m (fn
                   ([^double sample ^StateBiquad state]
                    (let [y (-> (* (.b0 c) sample)
                                (+ (* (.b1 c) (.x1 state)))
                                (+ (* (.b2 c) (.x2 state)))
                                (+ (* (.a1 c) (.y1 state)))
                                (+ (* (.a2 c) (.y2 state))))]
                      (SampleAndState. y (StateBiquad. (.x1 state) sample (.y1 state) y))))
                   ([] (StateBiquad. 0.0 0.0 0.0 0.0)))))

;; ### Biquad equalizer

(defmethod effect :biquad-eq
  ([m] (effect m {}))
  ([m {:keys [fc gain bw fs]
       :or {fc 1000.0 gain 0.0 bw 1.0 fs 44100.0}}]
   (make-biquad-filter m (biquad-eq-params fc gain bw fs))))

;; ### Biquad high/low shelf

(defmethod effect :biquad-hs
  ([m] (effect m {}))
  ([m {:keys [fc gain slope fs]
       :or {fc 1000.0 gain 0.0 slope 1.5 fs 44100.0}}]
   (make-biquad-filter m (biquad-hs-params fc gain slope fs))))

(defmethod effect :biquad-ls
  ([m] (effect m {}))
  ([m {:keys [fc gain slope fs]
       :or {fc 1000.0 gain 0.0 slope 1.5 fs 44100.0}}]
   (make-biquad-filter m (biquad-ls-params fc gain slope fs))))

;; ### Biquad lowpass/highpass/bandpass

(defn- lhb-params
  "Create parameters for lp hp and bp biquad filters."
  [f {:keys [fc bw fs]
      :or {fc 1000.0 bw 1.0 fs 44100.0}}]
  (f fc bw fs))

(defmethod effect :biquad-lp ([m] (effect m {})) ([m conf] (make-biquad-filter m (lhb-params biquad-lp-params conf))))
(defmethod effect :biquad-hp ([m] (effect m {})) ([m conf] (make-biquad-filter m (lhb-params biquad-hp-params conf))))
(defmethod effect :biquad-bp ([m] (effect m {})) ([m conf] (make-biquad-filter m (lhb-params biquad-bp-params conf))))

;; ### DJ Equalizer

(defmethod effect :dj-eq
  ([m] (effect m {}))
  ([m {:keys [hi mid low shelf-slope peak-bw ^double rate]
       :or {hi 0.0 mid 0.0 low 0.0 shelf-slope 1.5 peak-bw 1.0 rate 44100.0}}]
   (let [b (compose-effects
            (effect :biquad-hs {:fc (* rate (/ 10000.0 44100.0)) :gain hi :slope shelf-slope :fs rate})
            (effect :biquad-eq {:fc (* rate (/ 1000.0 44100.0)) :gain mid :bw peak-bw :fs rate})
            (effect :biquad-eq {:fc (* rate (/ 100.0 44100.0)) :gain low :bw peak-bw :fs rate}))]
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
   (let [a1 (/ (- 1.0 delay) (inc delay))]
     (effect-node m (fn
                      ([^double sample ^double zm1]
                       (let [y (+ zm1 (* sample (- a1)))
                             new-zm1 (+ sample (* y a1))]
                         (SampleAndState. y new-zm1)))
                      ([] 0.0))))))

;; ### Divider

(deftype StateDivider [^double out ^double amp ^double count ^double lamp ^double last ^int zeroxs])

(defmethod effect :divider
  ([m] (effect m {}))
  ([m {:keys [^long denom]
       :or {denom 2.0}}]
   (effect-node m
                (fn
                  ([^double sample ^StateDivider state]
                   (let [count (inc (.count state))
                         ^StateDivider s1 (if (or (and (> sample 0.0) (<= (.last state) 0.0))
                                                  (and (neg? sample) (>= (.last state) 0.0)))
                                            (if (== denom 1)
                                              (StateDivider. (if (pos? (.out state)) -1.0 1.0) 0.0 0.0 (/ (.amp state) count) (.last state) 0)
                                              (StateDivider. (.out state) (.amp state) count (.lamp state) (.last state) (inc (.zeroxs state))))
                                            (StateDivider. (.out state) (.amp state) count (.lamp state) (.last state) (.zeroxs state)))
                         amp (+ (.amp s1) (m/abs sample))
                         ^StateDivider s2 (if (and (> denom 1)
                                                   (== ^long (rem (.zeroxs s1) denom) (dec denom)))
                                            (StateDivider. (if (pos? (.out s1)) -1.0 1.0) 0.0 0 (/ amp (.count s1)) (.last s1) 0)
                                            (StateDivider. (.out s1) amp (.count s1) (.lamp s1) (.last s1) (.zeroxs s1)))]
                     (SampleAndState. (* (.out s2) (.lamp s2)) (StateDivider. (.out s2) (.amp s2) (.count s2) (.lamp s2) sample (.zeroxs s2)))))
                  ([] (StateDivider. 1.0 0.0 0.0 0.0 0.0 0.0))))))

;; ### FM filter

(deftype StateFm [^double pre ^double integral ^double t lp])

(defmethod effect :fm
  ([m] (effect m {}))
  ([m {:keys [^double quant ^double omega ^double phase]
       :or {quant 10.0 omega 0.014 phase 0.00822}}]
   (let [lp-chain (compose-effects (effect :simple-lowpass {:rate 100000 :cutoff 25000})
                                   (effect :simple-lowpass {:rate 100000 :cutoff 10000})
                                   (effect :simple-lowpass {:rate 100000 :cutoff 1000}))]
     (effect-node m
                  (fn
                    ([^double sample ^StateFm state]
                     (let [sig (* sample phase)
                           new-integral (+ (.integral state) sig)
                           m (m/cos (+ new-integral (* omega (.t state))))
                           m (if (pos? quant)
                               (m/norm (unchecked-int (m/norm m -1.0 1.0 0.0 quant)) 0.0 quant -1.0 1.0)
                               m)
                           dem (m/abs (- m (.pre state)))
                           ^EffectsList res (single-pass (.lp state) dem)
                           demf (/ (* 2.0 (- (.sample res) omega)) phase)]
                       (SampleAndState. (m/constrain demf -1.0 1.0) (StateFm. m new-integral (inc (.t state)) res))))
                    ([] (StateFm. 0.0 0.0 0.0 lp-chain)))))))


;; ### Bandwidth limit

(defmethod effect :bandwidth-limit
  ([m] (effect m {}))
  ([m {:keys [^double freq ^double rate]
       :or {freq 1000.0 rate 44100.0}}]
   (let [dx (/ freq rate)]
     (effect-node m (fn
                      ([^double sample ^double state]
                       (let [res (if (>= sample state)
                                   (min (+ state dx) sample)
                                   (max (- state dx) sample))]
                         (SampleAndState. res res)))
                      ([] 0.0))))))

;; ### Distortion

(defmethod effect :distort
  ([m] (effect m {}))
  ([m {:keys [^double factor]
       :or {factor 1.0}}]
   (let [nfact (inc factor)]
     (effect-node m (fn 
                      ([^double sample state]
                       (let [div (+ factor (m/abs sample))
                             res (* nfact (/ sample div))]
                         (SampleAndState. res state)))
                      ([]))))))

;; ### Fast overdrive

(defmethod effect :foverdrive
  ([m] (effect m {}))
  ([m {:keys [^double drive]
       :or {drive 2.0}}]
   (let [drivem1 (dec drive)]
     (effect-node m (fn
                      ([^double sample state]
                       (let [fx (m/abs sample)
                             res (/ (* sample (+ fx drive)) (inc (+ (* sample sample) (* fx drivem1))))]
                         (SampleAndState. res state)))
                      ([]))))))

;; ### Decimator

(deftype StateDecimator [^double count ^double last])

(defmethod effect :decimator
  ([m] (effect m {}))
  ([m {:keys [^double bits ^double fs ^double rate]
       :or {bits 2.0 fs 4410.0 rate 44100.0}}]
   (let [step (m/pow 0.5 (- bits 0.9999))
         stepr (/ step)
         ratio (/ fs rate)]
     (effect-node m (fn
                      ([^double sample ^StateDecimator state]
                       (let [ncount (+ (.count state) ratio)]
                         (if (>= ncount 1.0)
                           (let [delta (* step ^double (m/remainder (->> sample
                                                                         m/sgn
                                                                         (* step 0.5)
                                                                         (+ sample)
                                                                         (* stepr)) 1.0))
                                 last (- sample delta)]
                             (SampleAndState. last (StateDecimator. (dec ncount) last)))
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
         wb (/ (* m/TWO_PI bass-freq) rate)
         wt (/ (* m/TWO_PI treble-freq) rate) 
         cwb (m/cos wb)
         cwt (m/cos wt)
         ab (m/exp (/ (* 2.302585092994046 bass) 40.0))
         ab+ (inc ab)
         ab- (dec ab)
         at (m/exp (/ (* 2.302585092994046 treble) 40.0))
         at+ (inc at)
         at- (dec at)
         bb (m/sqrt (- (/ (inc (m/sq ab)) slope) (m/sq (dec ab))))
         bt (m/sqrt (- (/ (inc (m/sq at)) slope) (m/sq (dec at))))
         bswb (* bb (m/sin wb))
         bswt (* bt (m/sin wt))

         b0b (* ab (+ (- ab+ (* ab- cwb)) bswb))
         b1b (* 2.0 ab (- ab- (* ab+ cwb)))
         b2b (* ab (- (- ab+ (* ab- cwb)) bswb))
         a0b (+ (+ ab+ (* ab- cwb)) bswb)
         a1b (* -2.0 (+ ab- (* ab+ cwb)))
         a2b (- (+ ab+ (* ab- cwb)) bswb)

         b0t (* at (+ (+ at+ (* at- cwt)) bswt))
         b1t (* -2.0 at (+ at- (* at+ cwt)))
         b2t (* at (- (+ at+ (* at- cwt)) bswt))
         a0t (+ (- at+ (* at- cwt)) bswt)
         a1t (* 2.0 (- at- (* at+ cwt)))
         a2t (- (- at+ (* at- cwt)) bswt)]
     (effect-node m (fn
                      ([^double sample ^StateBassTreble state]
                       (let [outb (/ (-> (* b0b sample)
                                         (+ (* b1b (.xn1Bass state)))
                                         (+ (* b2b (.xn2Bass state)))
                                         (- (* a1b (.yn1Bass state)))
                                         (- (* a2b (.yn2Bass state)))) a0b)
                             outt (/ (-> (* b0t outb)
                                         (+ (* b1t (.xn1Treble state)))
                                         (+ (* b2t (.xn2Treble state)))
                                         (- (* a1t (.yn1Treble state)))
                                         (- (* a2t (.yn2Treble state)))) a0t)]
                         (SampleAndState. (* outt data-gain)
                                          (StateBassTreble. sample (.xn1Bass state) outb (.yn1Bass state)
                                                            outb (.xn1Treble state) outt (.yn1Treble state)))))
                      ([] (StateBassTreble. 0.0 0.0 0.0 0.0
                                            0.0 0.0 0.0 0.0)))))))

;; ### Echo (audacity)

(deftype StateEcho [^doubles buffer ^int position])

(defmethod effect :echo
  ([m] (effect m {}))
  ([m {:keys [^double delay ^double decay ^double rate]
       :or {delay 0.5 decay 0.5 rate 44100.0}}]
   (let [buffer-len (int (min 10000000 (* delay rate)))]
     (effect-node m (fn
                      ([^double sample ^StateEcho state]
                       (let [result (+ sample (* decay (aget ^doubles (.buffer state) (.position state))))]
                         (aset ^doubles (.buffer state) (.position state) result)
                         (SampleAndState. result (StateEcho. (.buffer state) (rem (inc (.position state)) buffer-len)))))
                      ([]
                       (StateEcho. (double-array buffer-len 0.0) 0)))))))

;; ### Vcf303
;;
(deftype StateVcf303 [^double d1 ^double d2 ^double c0 ^int env-pos ^Vec3 abc])

(defmethod effect :vcf303
  ([m] (effect m {}))
  ([m {:keys [^double rate trigger ^double cutoff ^double resonance ^double env-mod ^double decay ^double gain]
       :or {rate 44100.0 trigger false cutoff 0.8 resonance 0.8 env-mod 0.5 decay 1.0 gain 1.0}}]
   (let [scale (/ m/PI rate)
         e0 (* scale
               (m/exp (-> (- 5.613 (* 0.8 env-mod))
                          (+ (* 2.1553 cutoff))
                          (- (* 0.7696 (- 1.0 resonance))))))        
         d (m/pow (->> decay
                       (* 2.3)
                       (+ 0.2)
                       (* rate)
                       (/ 1.0)
                       (m/pow 0.1)) 64.0)
         r (m/exp (- (* 3.455 resonance) 1.20))
         recalc-abc (fn [^double vc0]
                      (let [whopping (+ e0 vc0)
                            k (m/exp (/ (- whopping) r))
                            a (* (+ k k) (m/cos (+ whopping whopping)))
                            b (* (- k) k)
                            c (* 0.2 (- (- 1.0 a) b))]
                        (Vec3. a b c)))
         init-c0       (if trigger
                         (- (* scale
                               (m/exp (-> (+ 6.109 (* 1.5876 env-mod))
                                          (+ (* 2.1553 cutoff))
                                          (- (* 1.2 (- 1.0 resonance)))))) e0)
                         0.0)]    
     (effect-node m (fn
                      ([^double sample ^StateVcf303 state]
                       (let [^Vec3 abc (.abc state)
                             result (-> (* (.x abc) (.d1 state))
                                        (+ (* (.y abc) (.d2 state)))
                                        (+ (* (.z abc) sample)))
                             d2 (.d1 state)
                             d1 result
                             env-pos (inc (.env-pos state))]
                         (if (>= env-pos 64)
                           (let [c0 (* d (.c0 state))]
                             (SampleAndState. (* gain result)
                                              (StateVcf303. d1 d2 c0 0 (recalc-abc c0))))
                           (SampleAndState. (* gain result)
                                            (StateVcf303. d1 d2 (.c0 state) env-pos abc)))))
                      ([] (StateVcf303. 0.0 0.0 init-c0 0 (recalc-abc init-c0))))))))

;; ### Slew limiter

(defmethod effect :slew-limit
  ([m] (effect m {}))
  ([m {:keys [^double rate ^double maxrise ^double maxfall]
       :or {rate 44100.0 maxrise 500.0 maxfall 500.0}}]
   (let [maxinc (/ maxrise rate)
         maxdec (- (/ maxfall rate))]
     (effect-node m (fn
                      ([^double sample ^double prev]
                       (let [increment (- sample prev) 
                             nsample (+ prev (m/constrain increment maxdec maxinc))]
                         (SampleAndState. nsample nsample)))
                      ([] 0.0))))))


;; 
(deftype StateMdaThruZero [^doubles buffer ^double ph ^long bp ^double f])

(defmethod effect :mda-thru-zero
  ([m] (effect m {}))
  ([m {:keys [^double rate ^double speed ^double depth ^double mix ^double depth-mod ^double feedback]
       :or {rate 44100.0 speed 0.3 depth 0.43 mix 0.47 feedback 0.3 depth-mod 1.0}}]
   (let [rat (/ (* (m/pow 10.0 (- 2.0 (* 3.0 speed))) 2.0) rate)
         dep (* 2000.0 (m/sq depth))
         dem (- dep (* dep depth-mod))
         dep (- dep dem)
         wet mix
         dry (- 1.0 wet)
         fb (- (* 1.9 feedback) 0.95)]
     (effect-node m (fn
                      ([^double sample ^StateMdaThruZero state]
                       (let [ph (+ (.ph state) rat)
                             ph (if (> ph 1.0) (- ph 2.0) ph)
                             bp (bit-and (dec (.bp state)) 0x7ff)]
                         (aset ^doubles (.buffer state) bp (+ sample (* fb (.f state))))
                         (let [tmpf (+ dem (* dep (- 1.0 (m/sq ph))))
                               tmp (unchecked-int tmpf)
                               tmpf (- tmpf tmp)
                               tmp (bit-and (+ tmp bp) 0x7ff)
                               tmpi (bit-and (inc tmp) 0x7ff)
                               f (aget ^doubles (.buffer state) tmp)
                               f (+ (* tmpf (- (aget ^doubles (.buffer state) tmpi) f)) f)
                               result (+ (* sample dry) (* f wet))]
                           (SampleAndState. result (StateMdaThruZero. (.buffer state) ph bp f)))))
                      ([] (StateMdaThruZero. (double-array 2048) 0.0 0 0.0)))))))

(def ^{:doc "List of effects."}
  effects-list (sort (keys (methods effect))))

;; ## File operations

(defn save-signal
  "Save signal to file.

  Representation is: 16 bit signed, big endian file
  You can use Audacity/SOX utilities to convert files to audio."
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
  "Read signal from file

  Expected representation is 16 bit signed, big endian file."
  [filename]
  (let [^java.io.File f (file filename)
        len (/ (.length f) 2)
        ^java.io.DataInputStream in (java.io.DataInputStream. (input-stream filename))
        ^doubles buffer (double-array len)]
    (try
      (dotimes [i len]
        (aset ^doubles buffer (int i) (double (m/cnorm (.readShort in) Short/MIN_VALUE Short/MAX_VALUE -1.0 1.0))))
      (finally (. in clojure.core/close)))
    buffer))

;; ## Signal generators
;;
;; Here you have defined multimethods to create waves from various oscilators
;;
;; Parameters are:
;;
;; * oscilator name (see `oscillators` variable)
;; * frequency
;; * amplitude
;; * phase (0-1)
;;
;; Multimethod creates oscillator function accepting `double` (time) and resulting `double` from [-1.0 1.0] range.

(defmulti oscillator
  "Create oscillator.

  Parameters are:

  * oscilator name (see `oscillators` variable)
  * frequency
  * amplitude
  * phase (0-1)
  
  Multimethod creates oscillator function accepting `double` (as time) and returns `double` from [-1.0 1.0] range.

  To convert `oscillator` to signal, call [[signal-from-oscillator]].

  To add oscillators, call [[sum-oscillators]]."
  (fn [f _ _ _] f))

(defmethod oscillator :sin [_ ^double f ^double a ^double p]
  (fn ^double [^double x]
    (* a
       (m/sin (+ (* p m/TWO_PI) (* x m/TWO_PI f))))))

(def ^:private snoise (r/fbm-noise {:noise-type :simplex
                                    :octaves 1                                    
                                    :normalize? false}))

(defmethod oscillator :noise [_ ^double f ^double a ^double p]
  (fn ^double [^double x]
    (* a ^double (snoise (* (+ p x) f) 1.23456789))))

(defmethod oscillator :saw [_ ^double f ^double a ^double p] 
  (fn ^double [^double x]
    (let [rp (* 2.0 a)
          p2 (* f (mod (+ (* a p) a x) 1.0))]
      (* rp (- p2 (m/floor p2) 0.5)))))

(defmethod oscillator :square [_ ^double f ^double a ^double p]
  (fn ^double [^double x]
    (if (< (mod (+ p (* x f)) 1.0) 0.5)
      a
      (- a))))

(defmethod oscillator :triangle [_ ^double f ^double a ^double p]
  (let [saw (oscillator :saw f a p)]
    (fn ^double [^double x]
      (- (* 2.0 (m/abs (double (saw x)))) a))))

(defmethod oscillator :cut-triangle [_ ^double f ^double a ^double p]
  (let [tri (oscillator :triangle f a p)]
    (fn ^double [^double x]
      (let [namp (* 0.5 a)]
        (* 2.0 (m/constrain (double (tri x)) (- namp) namp))))))

(defmethod oscillator :constant [_ _ ^double a _] (constantly a))

(def ^{:doc "List of oscillator names used with [[oscillator]]"}
  oscillators (sort (keys (methods oscillator))))

(defn oscillators-sum
  "Create oscillator which is sum of all oscillators."
  [& fs]
  (reduce #(fn ^double [^double x] (+ ^double (%1 x) ^double (%2 x))) fs))

(defn oscillator-gain
  [fs ^double gain]
  (fn ^double [^double x]
    (* gain ^double (fs x))))

(defn oscillator->signal
  "Create signal from oscillator.

  Parameters are:

  * f - oscillator
  * samplerate - in Hz
  * seconds - duration

  Returns sampled signal as double array."
  [f ^double samplerate ^double seconds]
  (let [len (* samplerate seconds)
        ^doubles buffer (double-array len)]
    (dotimes [i len]
      (aset ^doubles buffer i (m/constrain ^double (f (m/norm i 0 len 0 seconds)) -1.0 1.0)))
    buffer))

(defn signal->oscillator
  "Create oscillator from signal.

  Parameters:

  * sig - signal as sequence
  * seconds - duration
  * interpolator - interpolation (see [[fastmath.interpolation]]). Default: [[linear-smile]]."
  ([sig ^double seconds] (signal->oscillator sig seconds linear-interp/linear))
  ([sig ^double seconds interpolator]
   (let [c (count sig)
         step (/ seconds c)] 
     (interpolator (for [^long i (range c)]
                     (* i step)) sig))))

;; signal smoothing

(defn- perform-convolution
  [coeffs fc signal]
  (->> (MathArrays/convolve (m/seq->double-array signal) coeffs)
       (drop fc)
       (take (count signal))))

(defn savgol-filter
  "Creates Savitzky-Golay smoothing filter.

  Arguments:

  * length - length of the kernel (default: 5)
  * order - polynomial order (default: 2)
  * derivative - signal derivative (default: 0)

  Returns filtering function which accepts collection of numbers and returns filtered signal."
  ([] (savgol-filter 5))
  ([^long length] (savgol-filter length 2))
  ([^long length ^long order] (savgol-filter length order 0))
  ([^long length ^long order ^long derivative]
   (assert (odd? length) "Length must be odd!")
   (let [fc (/ (dec length) 2)
         coeffs (-> (for [v (range (- fc) (inc fc))]
                      (map #(m/pow v %) (range (inc order))))
                    (m/seq->double-double-array)
                    (Array2DRowRealMatrix.)
                    (SingularValueDecomposition.)
                    (.getSolver)
                    (.getInverse)
                    (.getRow derivative))]
     (fn [signal]
       (let [ns (perform-convolution coeffs fc signal)]
         (if (even? derivative)
           ns
           (map (fn [^double v] (* -1.0 v)) ns)))))))

(defn moving-average-filter
  "Creates moving average filter.

  Arguments:

  * length - length of the kernel (default: 5)

  Returns filtering function. See also [[savgol-filter]]."
  ([] (moving-average-filter 5))
  ([^long length] (savgol-filter length 1 0)))

(defn kernel-smoothing-filter
  "Creates Nadaraya-Watson kernel-weighted average

  Arguments:

  * kernel - [[kernel]] function (default `gaussian`)
  * length - length of the kernel (default: 5)
  * step - distance between consecutive samples (default: 1.0)

  Returns filtering function. See also [[savgol-filter]]."
  ([kernel] (kernel-smoothing-filter kernel 5))
  ([kernel ^long length] (kernel-smoothing-filter kernel length 1.0))
  ([kernel ^long length ^double step]
   (assert (odd? length) "Length must be odd!")
   (let [fc (/ (dec length) 2)
         coeffs (map (fn [^long v]
                       (kernel 0 (* step v))) (range (- fc) (inc fc)))
         ^double sum (reduce m/+ coeffs)
         coeffs (double-array (map (fn [^double v]
                                     (/ v sum)) coeffs))]
     (partial perform-convolution coeffs fc))))

(m/unuse-primitive-operators)
