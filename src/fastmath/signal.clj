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
            [fastmath.interpolation.linear :as linear-interp]
            [fastmath.transform :as trans]
            [fastmath.kernel :as ker]
            [fastmath.stats :as stats]
            [fastmath.special :as special]
            [fastmath.complex :as cplx]

            [fastmath.signal.waveform :as wv]
            [fastmath.signal.chirp :as chirp]
            [fastmath.signal.pad :as pad]
            [fastmath.signal.biquad :as biquad])
  (:import [fastmath.vector Vec3]
           [clojure.lang IFn]
           [org.apache.commons.math3.linear Array2DRowRealMatrix SingularValueDecomposition]
           [org.apache.commons.math3.util MathArrays]
           [fastmath.signal.biquad BiquadConf]
           [uk.me.berndporr.iirj Cascade]))

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
  "DB to Linear (power)"
  ^double [^double x]
  (m/exp10 (/ x 20.0)))

(defn linear->db
  "Linear (power) to DB"
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
   (make-biquad-filter m (biquad/equalizer fc gain bw fs))))

;; ### Biquad high/low shelf

(defmethod effect :biquad-hs
  ([m] (effect m {}))
  ([m {:keys [fc gain slope fs]
       :or {fc 1000.0 gain 0.0 slope 1.5 fs 44100.0}}]
   (make-biquad-filter m (biquad/highshelf fc gain slope fs))))

(defmethod effect :biquad-ls
  ([m] (effect m {}))
  ([m {:keys [fc gain slope fs]
       :or {fc 1000.0 gain 0.0 slope 1.5 fs 44100.0}}]
   (make-biquad-filter m (biquad/lowshelf fc gain slope fs))))

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


(defmethod effect :iir
  [m ^Cascade iir]
  (effect-node m (fn ([^double sample ^Cascade state]
                     (SampleAndState. (.filter state sample) state))
                   ([] (.reset iir) iir))))

(defmethod effect :gain
  ([m] (effect m 1.0))
  ([m ^double gain] (effect-node m (fn ([^double sample _] (SampleAndState. (m/* gain sample) nil))
                                     ([] nil)))))

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
                                            :else (m/- v (m/* m/THIRD (m/cb v))))))]
     (effect-node m (fn ([^double sample _]
                        (let [gsample (m/* pregain sample)] 
                          (SampleAndState. (f gsample) nil)))
                      ([] nil))))))

;; LADSPA version

(deftype SVF1State [^double h ^double b ^double l ^double p ^double n])

(defmethod effect :svf1
  ([m] (effect m {}))
  ([m {:keys [^double rate ^double cutoff kind ^double Q ^long oversamples ^double resonance]
       :or {rate 44100 cutoff 1000 kind :lowpass Q 0.5 oversamples 1 resonance 0.01}}]
   (let [fun (case kind
               :lowpass (fn ^double [^SVF1State s] (.l s))
               :highpass (fn ^double [^SVF1State s] (.h s))
               :bandpass (fn ^double [^SVF1State s] (.b s))
               :notch (fn ^double [^SVF1State s] (.n s))
               :allpass (fn ^double [^SVF1State s] (.p s)))
         f (m/* 2.0 (m/sinpi (m// cutoff (m/* rate oversamples))))
         q (m/* 2.0 (m/cos (m/* m/HALF_PI (m/pow Q 0.1))))
         qn (m/sqrt (m/+ (m/* 0.5 q) 0.01))]
     (effect-node m (fn ([^double sample ^SVF1State state]
                        (loop [i (long 0)
                               in (m/* qn (m/+ sample (m/* (.b state) resonance)))
                               ^SVF1State state state]
                          (if (m/== i oversamples)
                            (SampleAndState. (fun state) state)
                            (let [b (m/- (.b state) (m/* 0.001 (m/cb (.b state))))
                                  h (m/- in (.l state) (m/* b q))
                                  b (m/+ b (m/* f h))
                                  l (m/+ (.l state) (m/* f b))
                                  n (m/+ l h)
                                  p (m/- l h)
                                  nstate (SVF1State. h b l p n)]
                              (recur (m/inc i) (fun nstate) nstate)))))
                      ([] (SVF1State. 0.0 0.0 0.0 0.0 0.0)))))))

;; https://github.com/genmeblog/soundsynth/blob/master/src/sound/filter.clj

(def ^:const ^:private SVF2-M_PI_POW_2  (m/* m/M_PI m/M_PI))
(def ^:const ^:private SVF2-M_PI_POW_3  (m/* SVF2-M_PI_POW_2 m/M_PI))
(def ^:const ^:private SVF2-M_PI_POW_5  (m/* SVF2-M_PI_POW_3 SVF2-M_PI_POW_2))
(def ^:const ^:private SVF2-M_PI_POW_7  (m/* SVF2-M_PI_POW_5 SVF2-M_PI_POW_2))
(def ^:const ^:private SVF2-M_PI_POW_9  (m/* SVF2-M_PI_POW_7 SVF2-M_PI_POW_2))
(def ^:const ^:private SVF2-M_PI_POW_11 (m/* SVF2-M_PI_POW_9 SVF2-M_PI_POW_2))

(def ^:const ^:private SVF2-DIRTY_A    (m/* 3.739e-01 SVF2-M_PI_POW_3))
(def ^:const ^:private SVF2-FAST_A     (m/* 3.26e-01 SVF2-M_PI_POW_3))
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

(deftype SVF2State [^double state1 ^double state2 ^double lp ^double bp ^double bpn ^double hp])

(defmethod effect :svf2
  ([m] (effect m {}))
  ([m {:keys [^double rate ^double cutoff ^double Q kind tan]
       :or {rate 44100 cutoff 1000 kind :lowpass Q 0.5 tan :exact}}]
   (let [fun (case kind
               :lowpass (fn ^double [^SVF2State s] (.lp s))
               :highpass (fn ^double [^SVF2State s] (.hp s))
               :bandpass (fn ^double [^SVF2State s] (.bp s))
               :notch (fn ^double [^SVF2State s] (.bpn s)))
         ratio (m// cutoff rate)
         g (case tan
             :exact (tan-exact ratio)
             :dirty (tan-dirty ratio)
             :fast (tan-fast ratio)
             :accurate (tan-accurate ratio))
         r (m// Q)
         h (m// (m/+ 1.0 (m/* g (m/+ r g))))]
     (effect-node m (fn ([^double sample ^SVF2State state]
                        (let [hp (m/* h
                                      (m/- sample
                                           (m/* r (.state1 state))
                                           (m/* g (.state1 state))
                                           (.state2 state)))
                              bp (m/+ (m/* g hp)
                                      (.state1 state))
                              lp (m/+ (m/* g bp)
                                      (.state2 state))
                              nstate (SVF2State. (m/+ (m/* g hp) bp)
                                                 (m/+ (m/* g bp) lp)
                                                 lp bp (m/* r bp) hp)]
                          (SampleAndState. (fun nstate) nstate)))
                      ([] (SVF2State. 0.0 0.0 0.0 0.0 0.0 0.0)))))))


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
    :bandpass (Vec3. 0.0 1.0 1.0)
    :highpass (Vec3. 1.0 (m/- k) -1.0)
    :notch (Vec3. 1.0 (m/- k) 0.0)
    :peak (Vec3. 1.0 (m/- k) -2.0)
    :allpass (Vec3. 1.0 (m/* -2.0 k) 0.0)
    :bell (Vec3. 1.0 (m/* k (m/dec (m/* A A))) 0.0)
    :lowshelf (Vec3. 1.0 (m/* k (m/dec A)) (m/dec (m/* A A)))
    :highshelf (Vec3. (m/* A A) (m/* k (m/- 1.0 A) A) (m/- 1.0 (m/* A A)))))

(defmethod effect :svf3
  ([m] (effect m {}))
  ([m {:keys [^double rate ^double cutoff ^double Q kind tan ^double gaindb]
       :or {rate 44100 cutoff 1000 kind :lowpass Q 0.5 tan :exact gaindb 0.0}}]
   (let [A (m/exp10 (m// gaindb 40.0))
         Asqrt (m/sqrt A)
         ratio (m// cutoff rate)
         g (case tan
             :exact (tan-exact ratio)
             :dirty (tan-dirty ratio)
             :fast (tan-fast ratio)
             :accurate (tan-accurate ratio))
         k (if (= kind :BELL) (m// (m/* A Q)) (m// Q))
         ^Vec3 av (case kind
                    :lowshelf (compute-a (m// g Asqrt) k)
                    :highshelf (compute-a (m/* g Asqrt) k)
                    (compute-a g k))
         ^Vec3 mv (compute-m kind k A)]
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
                      ([] (SVF3State. 0.0 0.0 )))))))


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

;; convolution / correlation

(defn- get-slice
  [mode xs ^long size1 ^long size2]
  (cond
    (#{:all :full} mode) xs
    (#{:same :first} mode) (->> xs (drop (m// (m/- (count xs) size1) 2)) (take size1))
    :else (let [s (m/inc (m/abs (m/- size1 size2)))]
            (->> xs (drop (m// (m/- (count xs) s) 2)) (take s)))))

(defn convolve
  "Perform direct convolution with zero padding."
  ([sig1 sig2] (convolve sig1 sig2 :full))
  ([sig1 sig2 mode]
   (let [a1 (m/seq->double-array sig1)
         a2 (m/seq->double-array sig2)]
     (get-slice mode (MathArrays/convolve a1 a2) (alength a1) (alength a2)))))

(defn fft-convolve
  "Perform convolution using fft method."
  ([sig1 sig2] (fft-convolve sig1 sig2 :full))
  ([sig1 sig2 mode]
   (let [s1 (count sig1)
         s2 (count sig2)
         tsize (m/dec (m/+ s1 s2))
         size (m/round-up-pow2 tsize)
         fft1 (trans/fft (pad/zero (m/seq->double-array sig1) size :right))
         fft2 (trans/fft (pad/zero (m/seq->double-array sig2) size :right))]
     (get-slice mode (take tsize (trans/ifft (with-meta (mapv cplx/mult fft1 fft2) (meta fft1)))) s1 s2))))

(defn correlate
  "Perform direct correlation with zero padding"
  ([sig1 sig2] (correlate sig1 sig2 :full))
  ([sig1 sig2 mode]
   (convolve sig1 (reverse sig2) mode)))

(defn fft-correlate
  "Perform correlation using fft method"
  ([sig1 sig2] (fft-correlate sig1 sig2 :full))
  ([sig1 sig2 mode]
   (fft-convolve sig1 (reverse sig2) mode)))

;; filtering

(defn filter-signal
  "Apply FIR coefficients or IIR filter to a signal."
  [FIR-or-IIR xs]
  (if (instance? Cascade FIR-or-IIR)
    (map (fn [^double x] (.filter ^Cascade FIR-or-IIR x)) xs)
    (let [cf (if (m/< 2048 (m/+ (count xs) (count FIR-or-IIR)))
               convolve fft-convolve)]
      (cf xs FIR-or-IIR :same))))

(defn filter-signal-1
  "Apply IIR filter to a value"
  ([^Cascade IIR]
   (fn ^double [^double x] (.filter IIR x)))
  (^double [^Cascade IIR ^double x]
   (.filter IIR x)))

(defn reset-IIR!
  "Reset internal state of the IIR filter. Reset always to reuse on different signal."
  [^Cascade IIR]
  (.reset IIR)
  IIR)

(defn- maybe-reset-filter!
  [FIR-or-IIR xs]
  (when (instance? Cascade FIR-or-IIR) (.reset FIR-or-IIR))
  xs)

(defn filter-filter-signal
  "Apply filter twice (forward and backward) to align phase."
  [FIR-or-IIR xs]
  (->> (filter-signal FIR-or-IIR xs)
       (reverse)
       (maybe-reset-filter! FIR-or-IIR)
       (filter-signal FIR-or-IIR)
       (reverse)))

;; signal smoothing

(defn- perform-convolution
  [coeffs fc signal]
  (->> (convolve signal coeffs)
       (drop fc)
       (take (count signal))))

(defn savgol-filter
  "Creates Savitzky-Golay smoothing filter.

  Arguments:

  * length - length of the kernel (default: 5)
  * order - polynomial order (default: 2)
  * derivative - signal derivative (default: 0)

  Boundary rule is to pad with zeros.  

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

;; https://appliedacousticschalmers.github.io/scaling-of-the-dft/AES2020_eBrief/
;; https://dewesoft.com/blog/guide-to-fft-analysis
;; a book: https://brianmcfee.net/dstbook-site/content/intro.html


(defn fft-energy
  "Returns the energy spectrum (magnitude squared) of a transformed signal.

  Computes the energy for each frequency bin in a signal's Fourier representation. For real-valued signals, this function automatically accounts for energy conservation in single-sided spectra by doubling the energy of all bins except the DC component and the Nyquist frequency.

  Input parameters:

  * `txs` - A sequence of complex coefficients (frequency domain), typically `Vec2` objects produced by `fastmath.transform/fft`.
  * `options` - A map of configuration keys to override or provide metadata:
    * `:kind` - The type of the input spectrum: `:real` (default) or `:complex`.
    * `:nyquist?` - For real signals, a boolean indicating if the original time-domain signal length was even (required to correctly identify the Nyquist bin).

  Output:
  Returns a sequence of doubles representing the energy (squared magnitude) for each frequency bin."
  [txs options]
  (let [{:keys [kind nyquist?] :or {kind :real nyquist? true}} (merge (::fft (meta txs)) options)]
    (if (= :complex kind)
      (map v/magsq txs)
      (let [len- (m/dec (count txs))]
        (map-indexed (fn [^long id v]
                       (if (or (m/zero? id) (and nyquist? (m/== len- id)))
                         (v/magsq v)
                         (m/* 2.0 (v/magsq v)))) txs)))))

(defn fft-magnitude
  "Returns the magnitude spectrum of a transformed signal.

  Computes the absolute value (magnitude) for each frequency bin in a signal's Fourier representation. For real-valued signals, this function accounts for the single-sided spectrum representation by scaling coefficients (except DC and Nyquist) to ensure the magnitude correctly reflects the amplitude of the signal components.

  Input parameters:

  * `txs` - A sequence of complex coefficients (frequency domain), typically `Vec2` objects produced by `fastmath.transform/fft`.
  * `options` - A map of configuration keys to override or provide metadata:
    * `:kind` - The type of the input spectrum: `:real` (default) or `:complex`.
    * `:nyquist?` - For real signals, a boolean indicating if the original time-domain signal length was even (required to correctly identify the Nyquist bin).

  Output:
  Returns a sequence of doubles representing the magnitude (absolute value) for each frequency bin."
  [txs options]
  (-> txs (fft-energy options) (v/mult 2.0) v/sqrt))

(defn- fft-infer-N
  ^long [txs options]
  (let [options (merge (::fft (meta txs)) options)
        {:keys [kind nyquist?] :or {kind :real nyquist? true}} options
        N (count txs)]
    (if (= :real kind)
      (if nyquist? (m/* 2.0 (m/dec N)) (m/+ N (m/dec N)))
      N)))

(defn fft-amplitude
  "Returns the amplitude spectrum of a transformed signal.

  Computes the peak amplitude for each frequency bin by normalizing the magnitude spectrum. In Fourier analysis, raw FFT coefficients are proportional to the signal length $N$; this function divides the magnitudes by $N$ (and accounts for single-sided scaling in real signals) so that the resulting values correspond to the actual amplitudes of the sinusoidal components in the original time-domain signal.

  Input parameters:

  * `txs` - A sequence of complex coefficients (frequency domain), typically `Vec2` objects produced by `fastmath.transform/fft`.
  * `options` - A map of configuration keys to override or provide metadata:
    * `:kind` - The type of the input spectrum: `:real` (default) or `:complex`.
    * `:nyquist?` - For real signals, a boolean indicating if the original time-domain signal length was even (required to correctly calculate the normalization factor $N$).

  Output:
  Returns a sequence of doubles representing the normalized peak amplitude for each frequency bin."
  [txs options]
  (v/div (fft-magnitude txs options) (fft-infer-N txs options)))

(defn fft-power
  "Returns the power spectrum of a transformed signal.

  Computes the power distribution for each frequency bin by normalizing the energy spectrum by the square of the signal length $N$. This representation describes how much power (mean square amplitude) is contained in each frequency component. For real-valued signals, it automatically accounts for the single-sided spectrum scaling before normalization.

  Input parameters:

  * `txs` - A sequence of complex coefficients (frequency domain), typically `Vec2` objects produced by `fastmath.transform/fft`.
  * `options` - A map of configuration keys to override or provide metadata:
    * `:kind` - The type of the input spectrum: `:real` (default) or `:complex`.
    * `:nyquist?` - For real signals, a boolean indicating if the original time-domain signal length was even (required to correctly calculate the normalization factor $N$).

  Output:
  Returns a sequence of doubles representing the power (mean square) for each frequency bin."
  [txs options]
  (v/div (fft-energy txs options) (m/sq (fft-infer-N txs options))))

(defn fft-frequencies
  "Returns a sequence of frequency values corresponding to FFT bins.

  Generates the frequency axis for a signal's Fourier representation. It maps discrete bin indices to their physical frequency values based on the sampling rate and the signal length, facilitating the interpretation of the spectrum in Hertz (or the reciprocal units of the sampling interval).

  Input parameters:

  * `fs` - Sampling frequency (samples per second) of the original time-domain signal.
  * `N` - Total number of points in the FFT (typically the signal or window length).

  Output:
  Returns a sequence of `N` doubles representing the center frequency of each bin, starting from 0 (DC) up to the sampling frequency."
  [^double fs ^long N]
  (let [z (m// fs N)] (v/mult (range N) z)))

(defn- fft-times
  ([^double fs ^long shift] (fft-times fs shift 0))
  ([^double fs ^long shift ^long wlen]
   (let [step (m// shift fs)
         mid (m// wlen 2 fs)]
     (map (fn [^long n]
            (m/+ mid (m/* n step))) (range)))))

(defn stft
  "Computes the Short-Time Fourier Transform (STFT) of a 1D signal.

  The STFT provides a time-frequency representation of a signal by performing Fourier transforms over short, overlapping windowed segments of the data. This process captures how the frequency content of a non-stationary signal evolves over time, providing the underlying data structure for spectrograms.

  Input parameters:

  * `xs` - Input time-domain signal (sequence of doubles).
  * `options` - A map of configuration keys:
    * `:window` - A sequence of coefficients representing the window function (e.g., from `fastmath.kernel/window`). If not provided, a Gaussian window is automatically generated.
    * `:overlap` - The fraction of overlap between adjacent segments, typically between 0.0 and 1.0. Default is `0.5`.
    * `:fs` - The sampling frequency of the original signal. Default is `1.0`.
    * `:method` - The type of spectral values to return for each segment: `:magnitude`, `:amplitude`, `:energy`, `:power`, `:psd` (Power Spectral Density), `:phase`, or `nil` (returns raw `Vec2` complex coefficients). Default is `:power`.

  Returns a map containing:

  * `:N` - The total number of samples in the input signal.
  * `:spectrum` - A sequence of sequences, where each inner sequence represents the frequency spectrum of a time-localized window.
  * `:freqs` - A sequence of frequency values corresponding to the bins in the spectrum.
  * `:times` - A sequence of time values corresponding to the center of each windowed segment."
  ([xs {:keys [window ^double overlap ^double fs method db?]
        :or {overlap 0.5 fs 1.0 method :power}
        :as options}]
   (let [N (count xs)
         window (or window (ker/window :gaussian (m/max 5 (m/round (m/* 0.05 N)))))
         wlen (count window)
         shift (m/max 1 (m/round (m/* (m/- 1.0 overlap) wlen)))
         scale (v/sum window)
         scale-sq (m/sq scale)
         scale2 (m/* fs (v/sum (v/sq window)))
         xxs (->> (partition wlen shift xs)
                  (map (fn [xs]
                         (let [txs (-> (v/emult window xs)
                                       (trans/fft options))
                               s (case method
                                   :magnitude (fft-magnitude txs options)
                                   :amplitude (v/div (fft-magnitude txs options) scale)
                                   :energy (fft-energy txs options)
                                   :power (v/div (fft-energy txs options) scale-sq)
                                   :psd (v/div (fft-energy txs options) scale2)
                                   :phase (map v/heading txs)
                                   txs)]
                           (if (and db? (not (#{:phase} method)))
                             (map (fn [^double x] (m/* 10.0 (m/log10 (m/max x m/EPSILON)))) s)
                             s)))))]
     {:N N
      :spectrum xxs
      :freqs (take (count (first xxs)) (fft-frequencies fs wlen))
      :times (take (count xxs) (fft-times fs shift wlen))})))

(defn spectrum
  "Returns the frequency spectrum of a signal using FFT.

  Analyzes the frequency content of a time-domain signal or processes existing Fourier coefficients to produce various spectral representations. The function handles the necessary normalization and scaling to ensure the output reflects physical quantities such as peak amplitude, power, or spectral density.

  Input parameters:

  * `xs` - Input signal. Can be a sequence of doubles (time-domain) or a sequence of `fastmath.vector.Vec2` complex coefficients (frequency-domain).
  * `options` - A map of configuration keys:
    * `:method` - The type of spectral values to return. Options: `:magnitude`, `:amplitude` (default), `:energy`, `:power`, `:psd` (Power Spectral Density), `:asd` (Amplitude Spectral Density), `:phase` (angle in radians), `:real` (real part only), `:imag` (imaginary part only) or `:complex` (transformed data).
    * `:db?` - Boolean flag; if true, converts the resulting spectrum values to decibels using $10\\log_{10}(x)$. Default is `false`.
    * `:fs` - The sampling frequency of the signal. Default is `1.0`.
    * `:domain` - Specifies if the input `xs` is in the `:time` domain (requires performing an FFT) or already in the `:frequency` domain. Default is `:time`.

  Returns a map containing:

  * `:N` - The number of samples in the original signal.
  * `:spectrum` - A sequence of doubles representing the calculated spectral values.
  * `:freqs` - A sequence of frequency values (in the same units as `:fs`) corresponding to each bin in the spectrum."
  ([xs] (spectrum xs nil))
  ([xs {:keys [method db? ^double fs domain]
        :or {method :amplitude db? false fs 1.0 domain :time}
        :as options}]
   (let [N (count xs)
         step (m// fs N)
         sp (let [txs (if (= :time domain)
                        (trans/fft xs options)
                        xs)]
              (case method
                :magnitude (fft-magnitude txs options)
                :amplitude (fft-amplitude txs options)
                :energy (fft-energy txs options)
                :power (fft-power txs options)
                :psd (v/div (fft-power txs options) step)
                :asd (v/div (fft-amplitude txs options) step)
                :phase (map v/heading txs)
                :real (map first txs)
                :imag (map second txs)
                :complex txs))
         sp (if (and db? (not (#{:phase :real :imag :complex} method)))
              (map (fn [^double x] (m/* 10.0 (m/log10 (m/max x m/EPSILON)))) sp)
              sp)]
     {:N N
      :spectrum sp
      :freqs (take (count sp) (fft-frequencies fs N))})))

;; https://arxiv.org/pdf/gr-qc/0509116
(defn- median-bias
  [^long n]
  (let [n+ (m/inc n)]
    (- (special/digamma n+)
       (special/digamma (m/* 0.5 n+)))))

(defn periodogram
  "Estimate the spectral density of a signal using windowed averaging.

  Calculates the periodogram (an estimate of the spectral density) of a 1D signal by dividing the data into overlapping segments, computing the spectrum for each segment, and aggregating the results. This technique, based on Welch's method, reduces the noise and variance of the spectral estimate compared to a single FFT of the entire signal.

  Input parameters:
  * `xs` - Input time-domain signal (sequence of doubles).
  * `options` - A map of configuration keys:
    * `:method` - The type of spectral values to calculate for segments: `:psd` (Power Spectral Density, default), `:power`, `:magnitude`, `:amplitude`, `:energy`, or `:phase`.
    * `:average` - Aggregation method for segments: `:mean` (default), `:median` (more robust to outliers) or `:umedian` (unbiased median, like in SciPy).
    * `:window` - A sequence of coefficients for the window function.
    * `:overlap` - Fraction of overlap between segments.
    * `:fs` - Sampling frequency of the signal.

  Returns a map containing:
  * `:freqs` - A sequence of frequency values (bins).
  * `:spectrum` - A sequence of aggregated spectral values corresponding to the frequencies."
  ([xs] (periodogram xs nil))
  ([xs {:keys [method average] :or {method :psd average :mean} :as options}]
   (let [{:keys [freqs spectrum]} (stft xs (assoc options :method method))]
     {:freqs freqs
      :spectrum (case average
                  :mean (v/average-vectors spectrum)
                  :umedian (v/div (map stats/median (apply map vector spectrum)) (median-bias (count spectrum)))
                  :median (map stats/median (apply map vector spectrum)))})))

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
      * `:f` - Frequency in Hz. Can be a number or a function `(fn [t])` for frequency modulation. Default: `1.0`.
      * `:phase` - Phase offset in cycles [0.0 to 1.0]. Can be a number or a function `(fn [t])`. Default: `0.0`.
      * `:amplitude` - Peak amplitude. Can be a number or a function `(fn [t])` for amplitude modulation. Default: `1.0`.
      * `:duty` - Duty cycle for `:square` waves, ranging from 0.0 to 1.0 or `(fn [t])` for modulation. Default: `0.5`.
      * `:up?` - Boolean for `:saw` and `:analog-saw` waves. Set to `true` (default) for a rising ramp, `false` for a falling ramp.
      * `:bands` - Number of harmonics for `:band-limited-*` types. Increasing this value improves shape accuracy but increases computation. Default: `5`.

  Returns a function `(fn [t])` that accepts time `t` as a double and returns the signal value as a double."
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
     :band-limited-triangle (wv/->band-limited-triangle options))))

(defn chirp
  "Create a frequency-swept sine wave (chirp) oscillator.

  Generates a signal where the instantaneous frequency changes over a specified time interval according to a defined modulation profile. Chirps are essential in radar, sonar, and system identification to analyze frequency-dependent responses and impulse behavior.

  Input parameters:
  * `chirp-type` - A keyword specifying the sweep profile:
      * `:linear` - Frequency changes linearly with time.
      * `:quadratic-up` - Frequency increases quadratically.
      * `:quadratic-down` - Frequency decreases quadratically.
      * `:logarithmic` - Frequency changes exponentially (geometric sweep).
      * `:hyperbolic` - Frequency follows a hyperbolic curve.
  * `options` - A map of configuration keys:
      * `:f0` - Starting frequency in Hz at $t=0$.
      * `:f1` - Ending frequency in Hz at $t=time$.
      * `:time` - The time duration (in seconds) over which the sweep from `:f0` to `:f1` occurs.
      * `:amplitude` - Peak amplitude. Can be a number or a modulation function `(fn [t])`.
      * `:phase` - Initial phase offset in cycles [0.0 to 1.0].

  Returns a stateless function `(fn [t])` that accepts time `t` as a double and returns the signal value as a double."
  ([chirp-type] (chirp chirp-type nil))
  ([chirp-type options]
   (let [modulation (case chirp-type
                      :linear (chirp/->linear options)
                      :quadratic-up (chirp/->quadratic-up options)
                      :quadratic-down (chirp/->quadratic-down options)
                      :logarithmic (chirp/->logarithmic options)
                      :hyperbolic (chirp/->hyperbolic options))]
     (waveform :sine (merge options {:f modulation})))))

(defn add-waveforms
  "Adds two or more waveforms."
  ([w] w)
  ([w1 w2] (fn [^double t] (m/+ (double (w1 t)) (double (w2 t)))))
  ([w1 w2 w3] (fn [^double t] (m/+ (double (w1 t)) (double (w2 t)) (double (w3 t)))))
  ([w1 w2 w3 w4] (fn [^double t] (m/+ (double (w1 t)) (double (w2 t)) (double (w3 t)) (double (w4 t)))))
  ([w1 w2 w3 w4 w5] (fn [^double t] (m/+ (double (w1 t)) (double (w2 t)) (double (w3 t))
                                        (double (w4 t)) (double (w5 t)))))
  ([w1 w2 w3 w4 w5 & r]
   (reduce add-waveforms (add-waveforms w1 w2 w3 w4 w5) r)))

(defn gain-waveform
  "Scale amplitude with given `gain` value. Returns waveform function."
  [w ^double gain]
  (fn ^double [^double t] (m/* gain (double (w t)))))

(defn sample-waveform
  "Discretize a continuous waveform function into a sequence of samples.

  Evaluates a time-dependent function at regular intervals to produce a discrete-time signal.

  Input parameters:
  * `waveform-fn` - A function `(fn [t])` that accepts time as a double and returns the signal amplitude. Usually created via [[waveform]] or [[chirp]].
  * `fs` - Sampling frequency in Hz (samples per second).
  * `time` - The total duration of the signal to generate in seconds. Defaults to `1.0`.

  Returns a sequence of doubles representing the sampled signal in the time domain."
  ([waveform-fn ^double fs] (sample-waveform waveform-fn fs 1.0))
  ([waveform-fn ^double fs ^double time]
   (m/sample waveform-fn 0.0 time (m/* fs time))))

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
  "Calculates analytical signal."
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
  "Calculate envelope of the signal."
  [xs]
  (->> (hilbert xs) (map cplx/abs)))

(defn instantaneous-phase
  "Calculate instantaneous phase."
  ([xs] (instantaneous-phase xs true))
  ([xs unwrap?]
   (let [args (->> (hilbert xs) (map cplx/arg))]
     (if unwrap? (v/unwrap args m/TWO_PI) args))))

(defn instantaneous-frequency
  "Calculate instantaneous phase."
  ([xs] (instantaneous-frequency xs 1.0))
  ([xs ^double fs]
   (-> (instantaneous-phase xs)
       (v/differences)
       (v/div m/TWO_PI)
       (v/mult fs))))

;;

(defn zero-phase
  "Set phase to zero"
  [xs]
  (let [X (trans/fft xs)]
    (trans/ifft (with-meta (mapv (comp cplx/complex cplx/abs) X) (meta X)))))

;;

(defn resample
  "Resample signal using FFT interpolation."
  [xs ^long nsize]
  (let [size (count xs)]
    (if (m/== nsize size)
      xs
      (let [ratio (m// (double nsize) size)
            m (m/min size nsize)
            m2 (m/inc (m// m 2))
            X (as-> (vec (take m2 (trans/fft xs))) X
                (if (m/odd? m) X (update X (m// m 2) cplx/scale (if (m/< nsize size) 2.0 0.5)))
                (if (m/< nsize size) X
                    (concat X (repeat (if (m/even? size)
                                        (m// (m/- nsize size) 2)
                                        (m// (m/inc (m/- nsize size)) 2)) cplx/ZERO)))
                (map #(cplx/scale % ratio) X))]
        (trans/ifft X {:nyquist? (m/even? nsize)})))))



;; DEPRECATED

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
  {:deprecated "Use `waveform` and `chirp` functions."}
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
          p2 (* f (m/mod (+ (* a p) a x) 1.0))]
      (* rp (- p2 (m/floor p2) 0.5)))))

(defmethod oscillator :square [_ ^double f ^double a ^double p]
  (fn ^double [^double x]
    (if (< (m/mod (+ p (* x f)) 1.0) 0.5)
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

(def ^{:doc "List of oscillator names used with [[oscillator]]"
       :deprecated true}
  oscillators (sort (keys (methods oscillator))))

(defn oscillators-sum
  "Create oscillator which is sum of all oscillators."
  {:deprecated "Use `add-waveforms`."}
  [& fs]
  (reduce #(fn ^double [^double x] (+ ^double (%1 x) ^double (%2 x))) fs))

(defn oscillator-gain
  {:deprecated "Use `gain-waveform`."}
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
  {:deprecated "Use `sample-waveform`"}
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
  {:deprecated true}
  ([sig ^double seconds] (signal->oscillator sig seconds linear-interp/linear))
  ([sig ^double seconds interpolator]
   (let [c (count sig)
         step (/ seconds c)] 
     (interpolator (for [^long i (range c)]
                     (* i step)) sig))))
