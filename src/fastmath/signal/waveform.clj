(ns fastmath.signal.waveform
  "Waveform generators (oscillators), used by `fastmath.signal/waveform`.

  Every `->name` constructor takes an options map and returns a function of time `t` (seconds). The waveforms are periodic in the phase `phase + f t` (in cycles), also for negative time:

  * `:f` - frequency in Hz (default 1.0); a function `(fn [t])` is used directly as the phase in cycles (instead of `f t`), which allows frequency modulation, see `fastmath.signal.chirp`
  * `:phase` - phase offset in cycles (default 0.0), number or `(fn [t])`
  * `:amplitude` - peak amplitude (default 1.0), number or `(fn [t])`

  Shapes at phase 0 (as `scipy.signal` waveforms with a shifted phase):

  * [[->sine]] - `sin(2 pi phase)`
  * [[->square]] - 1 for the first `:duty` (default 0.5) part of the cycle, -1 after, as `scipy.signal.square(2 pi phase, duty)`
  * [[->saw]] - 0 rising to 1 at half a cycle, -1 after the jump (`:up? false` falling), as `scipy.signal.sawtooth(2 pi (phase + 1/2))`
  * [[->triangle]] - 0 rising to 1 at a quarter cycle, as `scipy.signal.sawtooth(2 pi (phase + 1/4), 0.5)`
  * `->analog-*` - cheap rational approximations of sine, saw and triangle with rounded shapes
  * `->band-limited-*` - first `:bands` (default 5) terms of the Fourier series of square, saw and triangle, without aliasing
  * [[polyblep-hyptri]] - lazy sequence of samples of a band-limited \"hyper triangle\" (PolyBLEP, Tale's `poly_blep.jsfx-inc`)"
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec3]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- cycle-frac
  "Fractional part in [0, 1), periodic also for negative values (`m/frac` is symmetric around 0)."
  ^double [^double x]
  (m/- x (m/floor x)))

(defn- ->time-fn [v]  (if (number? v) (let [v (double v)] (fn ^double [^double t] (m/* t v))) v))
(defn- ->const-fn [v]  (if (number? v) (let [v (double v)] (fn ^double [_] v)) v))

(defn- exp-curve
  ^double [^double x]
  (m// (m/+ 3.0 (m/* x (m/+ -13.0 (m/* 5.0 x))))
       (m/+ 3.0 x x)))

(defn ->sine
  "Sine wave: `amplitude sin(2 pi (phase + f t))`."
  ([] (->sine nil))
  ([{:keys [f phase amplitude]
     :or {f 1.0 phase 0.0 amplitude 1.0}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (fn ^double [^double t]
       (m/* (double (afn t))
            (m/sin (m/* m/TWO_PI (m/+ (double (pfn t)) (double (ffn t))))))))))

(defn- quadratic-sine-approx
  ^double [^double v]
  (let [phase (cycle-frac v)
        hphase? (m/< phase 0.5)
        x (m/- phase (if hphase? 0.25 0.75))
        v (m/- 1.0 (m/* 16.0 x x))]
    (if hphase? v (m/- v))))

(defn ->analog-sine
  "Piecewise quadratic approximation of the sine wave."
  ([] (->analog-sine nil))
  ([{:keys [f phase amplitude]
     :or {f 1.0 phase 0.0 amplitude 1.0}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (fn ^double [^double t]
       (m/* (double (afn t))
            (quadratic-sine-approx (m/+ (double (pfn t)) (double (ffn t)))))))))

(defn ->square
  "Square wave: `amplitude` for the first `:duty` part (0 to 1, default 0.5, number or `(fn [t])`) of every cycle, `-amplitude` for the rest."
  ([] (->square nil))
  ([{:keys [f phase amplitude duty]
     :or {f 1.0 phase 0.0 amplitude 1.0 duty 0.5}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)
         wfn (->const-fn duty)]
     (fn ^double [^double t]
       (m/* (double (afn t))
            (if (m/< (cycle-frac (m/+ (double (pfn t)) (double (ffn t))))
                     (m/constrain (double (wfn t)) 0.0 1.0))
              1.0 -1.0))))))

(defn ->saw
  "Sawtooth wave: 0 at phase 0, rising (`:up?` true, default) or falling to the jump at half a cycle."
  ([] (->saw nil))
  ([{:keys [f  phase amplitude up?]
     :or {f 1.0 phase 0.0 amplitude 1.0 up? true}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (if up?
       (fn ^double [^double t]
         (m/* (double (afn t))
              (m/dec (m/* 2.0 (cycle-frac (m/+ 0.5 (double (pfn t)) (double (ffn t))))))))
       (fn ^double [^double t]
         (m/* (double (afn t))
              (m/dec (m/* 2.0 (m/- 1.0 (cycle-frac (m/+ 0.5 (double (pfn t)) (double (ffn t)))))))))))))

(defn ->analog-saw
  "Rounded approximation of the sawtooth wave (`:up?` as in [[->saw]])."
  ([] (->analog-saw nil))
  ([{:keys [f phase amplitude up?]
     :or {f 1.0 phase 0.0 amplitude 1.0 up? true}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (if up?
       (fn ^double [^double t]
         (let [x (cycle-frac (m/+ 0.255969349108945 (double (pfn t)) (double (ffn t))))]
           (m/* -1.0 (double (afn t)) (exp-curve x))))
       (fn ^double [^double t]
         (let [x (cycle-frac (m/+ 0.255969349108945 (double (pfn t)) (double (ffn t))))]
           (m/* (double (afn t)) (exp-curve x))))))))

(defn ->triangle
  "Triangle wave: 0 at phase 0, 1 at a quarter cycle, -1 at three quarters."
  ([] (->triangle nil))
  ([{:keys [f phase amplitude]
     :or {f 1.0 phase 0.0 amplitude 1.0}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (fn ^double [^double t]
       (let [v (cycle-frac (m/+ (double (pfn t)) (double (ffn t))))]
         (m/* 4.0 (double (afn t))
              (double (cond
                        (m/< v 0.25) v
                        (m/< v 0.75) (m/- 0.5 v)
                        :else (m/dec v)))))))))

(defn- analog-triangle-approx
  ^double [^double phase]
  (let [x (-> (cycle-frac phase)
              (m/+ 0.1279846745544725)
              (cycle-frac))
        hx? (m/>= x 0.5)
        x (-> (m/* 2.0 x)
              (cycle-frac))]
    (if hx?
      (exp-curve x)
      (m/- (exp-curve x)))))

(defn ->analog-triangle
  "Rounded approximation of the triangle wave."
  ([] (->analog-triangle nil))
  ([{:keys [f phase amplitude]
     :or {f 1.0 phase 0.0 amplitude 1.0}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (fn ^double [^double t]
       (m/* (double (afn t))
            (analog-triangle-approx (m/+ (double (pfn t)) (double (ffn t)))))))))

;; band-limited

(defn- square-bl
  ^double [^long bands ^double t]
  (reduce (fn [^double s ^long b]
            (m/+ s (let [b+ (m/+ 1.0 b b)]                          
                     (m// (m/sin (m/* m/TWO_PI (m/* b+ t)))
                          b+)))) 0.0 (range bands)))

(defn ->band-limited-square
  "Square wave from the first `:bands` (default 5) odd harmonics: `4/pi sum sin(2 pi (2k+1) phase) / (2k+1)`."
  ([] (->band-limited-square nil))
  ([{:keys [f phase amplitude ^long bands]
     :or {bands 5 f 1.0 phase 0.0 amplitude 1.0}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (fn ^double [^double t]
       (m// (m/* (double (afn t))
                 (square-bl bands (m/+ (double (pfn t)) (double (ffn t)))))
            m/QUARTER_PI)))))

(defn- saw-bl
  ^double [^long bands ^double t]
  (reduce (fn [^double s ^long b]
            (m/+ s (let [b+ (m/inc b)]                          
                     (-> (m/sin (m/* m/TWO_PI (m/* b+ t)))
                         (m// b+)
                         (m/* (if (m/even? b) 1.0 -1.0)))))) 0.0 (range bands)))

(defn ->band-limited-saw
  "Sawtooth wave from the first `:bands` (default 5) harmonics: `2/pi sum (-1)^k sin(2 pi (k+1) phase) / (k+1)` (negated for `:up? false`)."
  ([] (->band-limited-saw nil))
  ([{:keys [f phase amplitude up? ^long bands]
     :or {f 1.0 phase 0.0 amplitude 1.0 up? true bands 5}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (if up?
       (fn ^double [^double t]
         (m// (m/* (double (afn t)) (saw-bl bands (m/+ (double (pfn t)) (double (ffn t)))))
              m/HALF_PI))
       (fn ^double [^double t]
         (m// (m/* (double (afn t)) (saw-bl bands (m/+ (double (pfn t)) (double (ffn t)))))
              m/-HALF_PI))))))

(defn- triangle-bl
  ^double [^long bands ^double t]
  (reduce (fn [^double s ^long b]
            (m/+ s (let [b+ (m/+ 1 b b)
                         b+2 (m/sq b+)]
                     (-> (m/sin (m/* m/TWO_PI (m/* b+ t)))
                         (m// b+2)
                         (m/* (if (m/even? b) 1.0 -1.0)))))) 0.0 (range bands)))

(defn ->band-limited-triangle
  "Triangle wave from the first `:bands` (default 5) odd harmonics: `8/pi^2 sum (-1)^k sin(2 pi (2k+1) phase) / (2k+1)^2`."
  ([] (->band-limited-triangle nil))
  ([{:keys [f phase amplitude ^long bands]
     :or {f 1.0 phase 0.0 amplitude 1.0 bands 5}}]
   (let [ffn (->time-fn f)
         afn (->const-fn amplitude)
         pfn (->const-fn phase)]
     (fn ^double [^double t]
       (m/* 0.8105694691387022 ;; 8/pi^2
            (double (afn t))
            (triangle-bl bands (m/+ (double (pfn t)) (double (ffn t)))))))))

;; poly_blep.jsfx-inc, licence: http://www.wtfpl.net/
;; blep/blamp/bluh

(defn blep
  "Band-limited step"
  ^double [^double t ^double dt]
  (cond
    (m/< t dt) (m/- (m/sq (m/dec (m// t dt))))
    (m/> t (m/- 1.0 dt)) (m/sq (m/inc (m// (m/dec t) dt)))
    :else 0.0))

(defn blamp
  "Band-limited ramp"
  ^double [^double t ^double dt]
  (cond
    (m/< t dt) (m/* -0.3333333333333333 (m/cb (m/dec (m// t dt))))
    (m/> t (m/- 1.0 dt)) (m/* 0.3333333333333333 (m/cb (m/inc (m// (m/dec t) dt))))
    :else 0.0))

(defn bluh
  "Band-limited curve"
  ^double [^double t ^double dt]
  (cond
    (m/< t dt)           (let [t (m/sq (m/dec (m// t dt)))]
                           (m/* -4.0 (m/- (m/sq t) t)))
    (m/> t (m/- 1.0 dt)) (let [t (m/sq (m/inc (m// (m/dec t) dt)))]
                           (m/* 4.0 (m/- (m/sq t) t)))
    :else                0.0))

;;

(defmacro polyblep
  "Defines `polyblep-<n>`, a sample generator from `fun` of `(phase, phase-step)` returning the sample value.

  The generated function takes options `:f` (Hz), `:phase` (cycles), `:amplitude` (numbers or functions of time in seconds) and `:fs` (sample rate, default 44100) and returns an infinite lazy sequence of samples."
  [n fun]
  `(defn ~(symbol (str "polyblep-" n))
     ~(str "Infinite lazy sequence of samples of the PolyBLEP " n " wave, see [[polyblep]] for the options.")
     [opts#]
     (let [f# (->const-fn (or (:f opts#) 1.0))
           phase# (->const-fn (or (:phase opts#) 0.0))
           amplitude# (->const-fn (or (:amplitude opts#) 1.0))
           fs# (double (or (:fs opts#) 44100.0))]
       (->> (iterate (fn [^Vec3 stp#]
                       (let [time-in-sec# (m// (.y stp#) fs#)
                             current-phase# (double (phase# time-in-sec#))
                             phase# (m/mod (m/+ (.z stp#) current-phase#) 1.0)
                             current-frequency# (double (f# time-in-sec#))
                             phase-step# (m// current-frequency# fs#)
                             current-amplitude# (double (amplitude# time-in-sec#))]
                         (Vec3. (m/* current-amplitude# (~fun phase# phase-step#))
                                (m/inc (.y stp#))
                                (m/+ (.z stp#) phase-step#)))) (Vec3. 0.0 0.0 0.0))
            (rest)
            (map first)))))
;;

(defn- pb-hyptri
  ^double [^double phase ^double phase-step]
  (m/+ -1.0
       (m/* phase-step (m/- (m// (blamp phase phase-step) m/QUARTER_PI)
                            (m/* (blamp (cycle-frac (m/+ phase 0.5)) phase-step) 9.273)))
       (m/* (m/exp (m/* 4.0 (if (m/< phase 0.5) phase (m/- 1.0 phase))))
            0.3130352854993313)))

#_{:clj-kondo/ignore [:unresolved-symbol]}
(polyblep hyptri pb-hyptri)

