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
  * `polyblep-*` - 46 quasi band-limited (PolyBLEP) waves, a port of Tale's [`poly_blep.jsfx-inc`](https://github.com/TaleTN/JSFX/blob/main/Effects/Tale/poly_blep.jsfx-inc) (WTFPL), see below

  ## PolyBLEP waves

  `(polyblep-<wave> options)` returns an infinite lazy sequence of samples (not a function of time), the sample `i` is taken at the time `i / fs`. The jumps, kinks and curvature changes of the naive wave are smoothed with the polynomials [[blep]], [[blamp]] and [[bluh]] over the two samples around them, which removes most of the aliasing. Options are described in [[polyblep]]: `:f`, `:phase`, `:amplitude`, `:fs` and the extra parameter of some waves: `:pw` (pulse width), `:n` (number of steps) or `:f2` (follower frequency of hard sync). Waves are silent for the frequency `f` of `fs / 4` and above (faded out from `0.2 fs`).

  * sine: `sinn` (`:n`), `half`, `half2` (`:pw`), `full`, `sinp`, `alt`, `camel`, `camela`, `camel2`, `sin2`
  * curves: `para`, `arc`, `circ` (the only wave with a state)
  * triangle and trapezoid: `tri`, `hyptri`, `trin` (`:n`), `tri2` (`:pw`), `trip` (`:pw`), `trap`, `trap2` (`:pw`)
  * square and rectangle: `sqr`, `lpsqr2`, `lpsqr3`, `intlpsqr2`, `intlpsqr3`, `bpsqr2`, `sqrm1`, `sqr2` (`:pw`), `rect` (`:pw`)
  * sawtooth: `saw`, `ramp`, `sawp`, `intsaw`, `cubsaw`, `sinsaw`, `hpsaw6`, `sawm1`, `sawn` (`:n`), `hwsaw` (`:pw`), `logit3`
  * staircase: `stairs`, `stairs3` (`:pw`), `stairs2` (`:pw`)
  * hard sync: `hstri`, `hssqr`, `hssaw` (`:f2`)"
  (:require [fastmath.core :as m]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- cycle-frac
  "Fractional part in [0, 1), periodic also for negative values (`m/frac` is symmetric around 0)."
  ^double [^double x]
  (let [r (m/- x (m/floor x))]
    ;; a tiny negative `x` rounds `x - floor(x)` to 1.0
    (if (m/< r 1.0) r 0.0)))

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

;; PolyBLEP waveforms, a port of Tale's `poly_blep.jsfx-inc`
;; https://github.com/TaleTN/JSFX/blob/main/Effects/Tale/poly_blep.jsfx-inc
;; licence: WTFPL, http://www.wtfpl.net/
;; Adapted from "Phaseshaping Oscillator Algorithms for Musical Sound Synthesis" by
;; Jari Kleimola, Victor Lazzarini, Joseph Timoney and Vesa Valimaki.

(defn blep
  "Band-limited step residual (PolyBLEP).

  `t` is the phase in [0, 1), `dt` the phase step per sample (frequency / sample rate).
  Returns a 2nd order polynomial in the two samples around the discontinuity (`t` just after, or just before, the wrap from 1 to 0), 0 otherwise.
  Add it to a naive waveform with a jump of height 1 (subtract for a jump down) to remove most of the aliasing."
  ^double [^double t ^double dt]
  (cond
    (m/< t dt) (m/- (m/sq (m/dec (m// t dt))))
    (m/> t (m/- 1.0 dt)) (m/sq (m/inc (m// (m/dec t) dt)))
    :else 0.0))

(defn blamp
  "Band-limited ramp residual (PolyBLAMP), the integral of [[blep]].

  `t` is the phase in [0, 1), `dt` the phase step per sample.
  Returns a 3rd order polynomial around a change of slope (of size 1 per unit phase) at the wrap from 1 to 0, 0 otherwise."
  ^double [^double t ^double dt]
  (cond
    (m/< t dt) (m/* -0.3333333333333333 (m/cb (m/dec (m// t dt))))
    (m/> t (m/- 1.0 dt)) (m/* 0.3333333333333333 (m/cb (m/inc (m// (m/dec t) dt))))
    :else 0.0))

(defn bluh
  "Band-limited curve residual, a 4th order polynomial (the integral of [[blamp]], scaled).

  `t` is the phase in [0, 1), `dt` the phase step per sample.
  Returns a polynomial around a change of curvature at the wrap from 1 to 0, 0 otherwise."
  ^double [^double t ^double dt]
  (cond
    (m/< t dt)           (let [t (m/sq (m/dec (m// t dt)))]
                           (m/* -4.0 (m/- (m/sq t) t)))
    (m/> t (m/- 1.0 dt)) (let [t (m/sq (m/inc (m// (m/dec t) dt)))]
                           (m/* 4.0 (m/- (m/sq t) t)))
    :else                0.0))

(defn- blep-at
  "[[blep]] of the phase `t + off` wrapped to [0, 1)."
  ^double [^double t ^double off ^double dt]
  (blep (cycle-frac (m/+ t off)) dt))

(defn- blamp-at
  "[[blamp]] of the phase `t + off` wrapped to [0, 1)."
  ^double [^double t ^double off ^double dt]
  (blamp (cycle-frac (m/+ t off)) dt))

(defn- bluh-at
  "[[bluh]] of the phase `t + off` wrapped to [0, 1)."
  ^double [^double t ^double off ^double dt]
  (bluh (cycle-frac (m/+ t off)) dt))

(defmacro ^:private bool->double
  [c]
  `(if ~c 1.0 0.0))

(defn- fade-gain
  "JSFX gain `a`: 1 up to the phase step 0.2, fading to 0 at 0.25 (above that the fundamental is above sample rate / 4 and the wave is silent)."
  ^double [^double dt]
  (cond
    (m/<= dt 0.2) 1.0
    (m/< dt 0.25) (m/- 1.0 (m/sq (m/* (m/- dt 0.2) 20.0)))
    :else 0.0))

;; Waveform bodies.
;; `t` is the phase in [0, 1), `dt` the phase step, `p` the extra parameter (pulse width, number of steps or hard sync frequency ratio).
;; Bodies return the sample without the amplitude and the `fade-gain`, `polyblep-samples` applies both.

;; sine family

(defn- pb-sinn
  ^double [^double t ^double dt ^double p]
  (let [n (long p)
        stepped? (m/>= n 2)
        dt2 (if stepped? (m// 0.5 p) 0.25)
        y1 (double (cond
                     stepped? (m/sin (m// (m/* (m/floor (m/+ (m/* t p 2.0) 0.5)) m/PI) p))
                     (m/< t 0.5) 0.7071067811865476
                     :else -0.7071067811865476))
        t1 (if stepped? (m/+ t (m/* dt2 0.5)) t)
        step3 (m/* dt2 m/TWO_PI)]
    (loop [i 0
           t1 t1
           t2 (m/+ t1 0.5)
           t3 step3
           s0 0.0
           y2 0.0]
      (if (m/< i n)
        (let [t1 (cycle-frac t1)
              t2 (cycle-frac t2)
              s (m/sin t3)]
          (recur (m/long-inc i) (m/+ t1 dt2) (m/+ t2 dt2) (m/+ t3 step3) s
                 (m/+ y2 (m/* (m/- (blep t1 dt) (blep t2 dt)) (m/- s s0)))))
        (m/+ (m/* (if stepped? 0.5 0.7071067811865476) y2) y1)))))

(defn- pb-half
  ^double [^double t ^double dt]
  (m/- (m/+ (m/* (m/+ (blamp t dt) (blamp-at t 0.5 dt)) m/TWO_PI dt)
            (if (m/< t 0.5) (m/* 2.0 (m/sin (m/* m/TWO_PI t))) 0.0))
       m/TWO_INV_PI))

(defn- pb-half2
  ^double [^double t ^double dt ^double pw]
  (let [pw (m/max pw 1.0e-4)
        t1u (m/+ (m/* pw 0.5) t 0.75)
        t1 (cycle-frac t1u)]
    (m/- (m/+ (m/* (m/+ (blamp (cycle-frac (m/- t1u pw)) dt) (blamp t1 dt))
                   (m// m/PI pw) dt)
              (if (m/< t1 pw) (m/* 2.0 (m/sin (m/* (m// t1 pw) m/PI))) 0.0))
         (m/* pw m/FOUR_INV_PI))))

(defn- pb-full
  ^double [^double t ^double dt]
  (let [t (cycle-frac (m/+ t 0.25))]
    (m/- (m/+ (m/* (blamp t dt) m/TWO_PI dt)
              (m/* 2.0 (m/sin (m/* t m/PI))))
         m/FOUR_INV_PI)))

(defn- pb-sinp
  ^double [^double t ^double dt]
  (m/- (m/+ (m/- (blep-at t 0.5 dt))
            (m/* (blamp t dt) m/PI dt)
            (if (m/< t 0.5) (m/* 2.0 (m/sin (m/* t m/PI))) 0.0))
       m/TWO_INV_PI))

(defn- pb-alt
  ^double [^double t ^double dt]
  (m/+ (m/* (m/- (blamp-at t 0.5 dt) (blamp t dt)) (m/- m/TWO_PI) dt)
       (if (m/< t 0.5) (m/sin (m/* 4.0 m/PI t)) 0.0)))

(defn- pb-camel
  ^double [^double t ^double dt]
  (m/- (m/+ (m/* (m/+ (m/* 2.0 (blamp-at t 0.75 dt))
                      (blamp-at t 0.5 dt)
                      (blamp t dt))
                 (m/* 4.0 m/PI) dt)
            (if (m/< t 0.5) (m/* 2.0 (m/abs (m/sin (m/* 4.0 m/PI t)))) 0.0))
       m/TWO_INV_PI))

(defn- pb-camela
  ^double [^double t ^double dt]
  (let [y (m/sin (m/* 4.0 m/PI t))
        y (if (m/< t 0.5) (m/abs y) y)]
    (m/- (m/+ (m/* (m/+ (blamp-at t 0.5 dt) (blamp-at t 0.75 dt)) (m/* 4.0 m/PI) dt)
              y)
         m/INV_PI)))

(defn- pb-camel2
  ^double [^double t ^double dt]
  (let [t1 (cycle-frac (m/+ t 0.25))]
    (m/+ (m/* (m/- (blamp-at t 0.75 dt) (blamp t1 dt)) (m/* 4.0 m/PI) dt)
         (m/sin (m/* (m/abs (m/- t1 0.5)) 4.0 m/PI)))))

(defn- pb-sin2
  ^double [^double t ^double dt]
  (m/+ (m/* (m/- (bluh-at t 0.5 dt) (bluh t dt)) (m/sq dt) -1.7724538509055160)
       (m/* (if (m/< t 0.5) 1.0 -1.0) (m/sq (m/sin (m/* m/TWO_PI t))))))

;; curves

(defn- pb-para
  ^double [^double t ^double dt]
  (m/+ (m/* (m/- (bluh-at t 0.5 dt) (bluh t dt)) 0.6666666666666666 (m/sq dt))
       (m/* (m/- (m/* t 16.0) 8.0) (m/- (m/abs (m/- t 0.5)) 0.5))))

(defn- pb-arc
  ^double [^double t ^double dt]
  (let [t (cycle-frac (m/+ t 0.25))]
    (m/+ (m/* (blamp t dt) (m/sqrt (m/abs dt)) m/TWO_PI)
         (m/* 4.0 (m/sqrt (m/- t (m/sq t))))
         -1.0)))

;; triangle family

(defn- pb-tri
  ^double [^double t ^double dt]
  (m/+ (m/- (m/abs (m/- (m/* (m/- (m/floor (m/+ t 0.75)) t) 4.0) 1.0)) 1.0)
       (m/* 4.0 dt (m/- (blamp-at t 0.25 dt) (blamp-at t 0.75 dt)))))

(defn- pb-hyptri
  ^double [^double t ^double dt]
  (m/+ -1.0
       (m/* dt (m/- (m/* (blamp t dt) m/FOUR_INV_PI)
                    (m/* (blamp-at t 0.5 dt) 9.273)))
       (m/* (m/exp (m/* 4.0 (if (m/< t 0.5) t (m/- 1.0 t))))
            0.3130352854993313)))

(defn- pb-trin
  ^double [^double t ^double dt ^double p]
  (let [n (long p)
        dt2 (m// 0.5 p)
        t1 (cycle-frac (m/+ (m/* dt2 0.5) 0.25 t))
        y1 (m// (m/* (m/floor (m/* t1 p 2.0)) 2.0) p)
        y1 (if (m/< t1 0.5) (m/dec y1) (m/- 3.0 y1))]
    (loop [i 0
           t1 t1
           t2 (m/+ t1 0.5)
           y2 0.0]
      (if (m/< i n)
        (let [y2 (m/- y2 (blep t1 dt))
              t1 (cycle-frac (m/+ t1 dt2))
              t2 (cycle-frac t2)
              y2 (m/+ y2 (blep t2 dt))]
          (recur (m/long-inc i) t1 (m/+ t2 dt2) y2))
        (m/+ (m// y2 p) y1)))))

(defn- pb-tri2
  ^double [^double t ^double dt ^double pw]
  (let [pw (m/min (m/max pw 1.0e-4) 0.9999)
        t1u (m/+ (m/* pw 0.5) t)
        t1 (cycle-frac t1u)
        t2 (cycle-frac (m/inc (m/- t1u pw)))]
    (m/+ (m/* (m// (m/- (blamp t1 dt) (blamp t2 dt)) (m/- pw (m/sq pw))) dt)
         (m/dec (m/* 2.0 (if (m/< t1 pw)
                           (m// t1 pw)
                           (m/inc (m// (m/- t1 pw) (m/dec pw)))))))))

(defn- pb-trip
  ^double [^double t ^double dt ^double pw]
  (let [pw (m/max pw 1.0e-4)
        t3 (m/+ t 0.75)
        t1u (m/+ t3 (m/* pw 0.5))
        t1 (cycle-frac t1u)]
    (m/- (m/+ (m/* (m// (m/+ (m/* -2.0 (blamp (cycle-frac t3) dt))
                             (blamp (cycle-frac (m/- t1u pw)) dt)
                             (blamp t1 dt))
                        pw)
                   2.0 dt)
              (if (m/< t1 pw) (m/- 2.0 (m/abs (m/- (m// (m/* t1 4.0) pw) 2.0))) 0.0))
         pw)))

(defn- pb-trap
  ^double [^double t ^double dt]
  (m/+ (m/min (m/max -1.0 (m/* 2.0 (m/dec (m/abs (m/inc (m/* (m/- (m/floor (m/+ t 0.75)) t) -4.0)))))) 1.0)
       (m/* 4.0 dt (m/+ (blamp-at t 0.125 dt)
                        (m/- (blamp-at t 0.625 dt))
                        (blamp-at t 0.375 dt)
                        (m/- (blamp-at t 0.875 dt))))))

(defn- pb-trap2
  ^double [^double t ^double dt ^double pw]
  (let [pw (m/min pw 0.9999)
        scale (m// 1.0 (m/- 1.0 pw))
        t1 (m/+ t (m/* pw 0.25) 0.25)
        t2 (m/+ t1 0.5)
        t3 (m/- t2 (m/* pw 0.5))
        t4 (m/+ t3 0.5)]
    (m/+ (m/min (m/max -1.0 (m/* scale (m/dec (m/abs (m/dec (m/* (m/- (m/floor (m/+ t 0.75)) t) 4.0)))))) 1.0)
         (m/* (m/+ (blamp (cycle-frac t1) dt)
                   (m/- (blamp (cycle-frac t2) dt))
                   (m/- (blamp (cycle-frac t3) dt))
                   (blamp (cycle-frac t4) dt))
              scale 2.0 dt))))

;; square family

(defn- pb-sqr
  ^double [^double t ^double dt]
  (m/+ (m/dec (m/* 2.0 (bool->double (m/< t 0.5))))
       (blep t dt)
       (m/- (blep-at t 0.5 dt))))

(defn- pb-lpsqr2
  ^double [^double t ^double dt]
  (let [y (m/* 2.0 (if (m/< t 0.5)
                     (m/- 0.5 (m/sq (m/dec (m/* t 2.0))))
                     (m/- (m/sq (m/* (m/- 1.0 t) 2.0)) 0.5)))]
    (m/+ y (m/* 4.0 dt (m/- (blamp t dt) (blamp-at t 0.5 dt))))))

(defn- pb-lpsqr3
  ^double [^double t ^double dt]
  (let [y (if (m/< t 0.5)
            (m/inc (m/* 2.0 (m/cb (m/dec (m/* t 2.0)))))
            (m/dec (m/* 2.0 (m/cb (m/* (m/- 1.0 t) 2.0)))))]
    (m/+ y (m/* 6.0 dt (m/- (blamp t dt) (blamp-at t 0.5 dt))))))

(defn- pb-intlpsqr2
  ^double [^double t ^double dt]
  (let [t0u (m/+ t 0.41503077834420436)
        t1 (m/+ t0u 0.5)
        t2 (m/+ t1 0.5)
        t0w (cycle-frac t0u)
        low? (m/< t0w 0.5)
        y (if low? 0.5 -0.5)
        t0 (if low? t0w (m/- t0w 0.5))]
    (m/* (m/+ (m/* y (m/+ (m/* (m/+ (m/* (m/- (m/* 64.0 t0) 96.0) t0) 24.0) t0) 2.0))
              (m/* (m/- (bluh (cycle-frac t1) dt) (bluh (cycle-frac t2) dt)) (m/sq dt)))
         -0.54691816067802716)))

(defn- pb-intlpsqr3
  ^double [^double t ^double dt]
  (let [t0u (m/+ t 0.37399153968783991)
        t1 (m/+ t0u 0.5)
        t2 (m/+ t1 0.5)
        t0w (cycle-frac t0u)
        low? (m/< t0w 0.5)
        y (if low? -0.5 0.5)
        t0 (if low? t0w (m/- t0w 0.5))]
    (m/* (m/+ (m/* y (m/- (m/* (m/- (m/* (m/+ (m/* (m/- (m/* 64.0 t0) 128.0) t0) 96.0) t0) 16.0) t0) 2.0))
              (m/* (m/- (bluh (cycle-frac t1) dt) (bluh (cycle-frac t2) dt)) (m/sq dt)))
         -0.72405970419833828)))

(defn- pb-bpsqr2
  ^double [^double t ^double dt]
  (let [t2 (cycle-frac (m/+ t 0.5))]
    (m/+ (m/* dt (m/+ (m/* 3.2 (m/- (blamp t dt) (blamp t2 dt)))
                      (m/* 0.5333333333333333 dt (m/- (bluh t2 dt) (bluh t dt)))))
         (m/* 12.8 (if (m/< t 0.5)
                     (m/- (m/* (m/- 0.75 t) t) 0.0625)
                     (m/+ (m/* (m/- t 1.75) t) 0.6875))))))

(defn- pb-sqrm1
  ^double [^double t ^double dt]
  (m/+ (m/dec (m/* 2.0 (m/+ (m/* (m/sin (m/* m/TWO_PI t))
                                 (m/- (m/* (m/sq dt) 2.0929829946242788) m/TWO_INV_PI))
                            (bool->double (m/< t 0.5)))))
       (blep t dt)
       (m/- (blep-at t 0.5 dt))))

(defn- pb-sqr2
  ^double [^double t ^double dt ^double pw]
  (let [t1u (m/+ (m/* pw 0.25) t 0.75)
        t2 (m/- t1u 0.5)
        t3u (m/- t2 (m/* pw 0.5))
        t4 (m/+ t3u 0.5)
        t1 (cycle-frac t1u)
        t3 (cycle-frac t3u)]
    (m/+ (m/* 0.5 (m/+ (blep t1 dt)
                       (m/- (blep (cycle-frac t2) dt))
                       (blep t3 dt)
                       (m/- (blep (cycle-frac t4) dt))))
         (bool->double (m/< t1 0.5))
         (bool->double (m/< t3 0.5))
         -1.0)))

(defn- pb-rect
  ^double [^double t ^double dt ^double pw]
  (m/+ (m/* 2.0 (m/- (bool->double (m/< t pw)) pw))
       (blep t dt)
       (m/- (blep (cycle-frac (m/inc (m/- t pw))) dt))))

;; sawtooth family

(defn- pb-saw
  ^double [^double t ^double dt]
  (let [t (cycle-frac (m/+ t 0.5))]
    (m/+ (m/- (blep t dt)) (m/dec (m/* t 2.0)))))

(defn- pb-ramp
  ^double [^double t ^double dt]
  (m/+ (blep t dt) (m/inc (m/* t -2.0))))

(defn- pb-sawp
  ^double [^double t ^double dt]
  (let [t (cycle-frac (m/+ t 0.5))]
    (m/+ (m/- (blep t dt))
         (if (m/>= t 0.5) (m/inc (m/cos (m/* m/TWO_PI t))) 0.0)
         -0.5)))

(defn- pb-intsaw
  ^double [^double t ^double dt]
  (m/* -8.0 (m/- (m/* (blamp-at t 0.5 dt) dt)
                 (m/- (m/+ (m/* (m/inc (m/* t -2.0)) (bool->double (m/>= t 0.5)))
                           (m/sq t))
                      0.08333333333333333))))

(defn- pb-cubsaw
  ^double [^double t ^double dt]
  (let [t1 (m/+ t 0.5)
        t0 (m/+ (m/- (bool->double (m/>= t1 1.0)) t1) 0.5)]
    (m/* 20.784609690826528
         (m/+ (m/* (bluh (cycle-frac t1) dt) (m/sq dt) 0.0625)
              (m/* (m/- (m/sq t0) 0.25) t0)))))

(defn- pb-sinsaw
  ^double [^double t ^double dt]
  (m/+ (blep t dt) (m/sin (m/* (m/- 0.5 t) m/PI))))

(defn- pb-hpsaw6
  ^double [^double t ^double dt]
  (let [t (cycle-frac (m/- 1.5 t))]
    (m/+ (m/* -2.0 (m/sq (m/* (m/sq t) t)))
         0.2857142857142857
         (blep t dt)
         (m/* 6.0 dt (blamp t dt)))))

(defn- pb-sawm1
  ^double [^double t ^double dt]
  (let [t2 (cycle-frac (m/+ t 0.5))]
    (m/- (m/dec (m/* 2.0 (m/+ (m/* (m/sin (m/* m/TWO_PI t))
                                   (m/- (m/* (m/sq dt) 1.0464914973121394) m/INV_PI))
                              t2)))
         (blep t2 dt))))

(defn- pb-sawn
  ^double [^double t ^double dt ^double p]
  (let [n (long p)
        t (cycle-frac (m/+ t 0.5))
        y1 (m/- (m/dec (m/* 2.0 (m// (m/floor (m/* (m/inc p) t)) p)))
                (blep t dt))
        dt2 (m// 1.0 (m/inc p))]
    (loop [i 0
           t t
           y2 0.0]
      (if (m/< i n)
        (let [t (cycle-frac (m/+ t dt2))]
          (recur (m/long-inc i) t (m/+ y2 (blep t dt))))
        (m/+ (m// y2 p) y1)))))

(defn- pb-hwsaw
  ^double [^double t ^double dt ^double pw]
  (let [pw (m/max pw 1.0e-4)
        t2 (cycle-frac (m/+ t 0.5))
        t1 (cycle-frac (m/+ t 0.5 pw))]
    (m/- (m/+ (m/- (if (m/< t1 pw) (m// (m/* t1 2.0) pw) 0.0) pw)
              (m/* (m// (blamp t1 dt) pw) dt))
         (m/+ (m/* (m// (blamp t2 dt) pw) dt) (blep t2 dt)))))

(defn- pb-logit3
  ^double [^double t ^double dt]
  (let [t (cycle-frac (m/+ t 0.5))]
    (m/- (m/dec (m/* (m/+ (m/* (m/- (m/* t 8.0) 12.0) t) 6.0) t))
         (blep t dt))))

;; staircases

(defn- pb-stairs
  ^double [^double t ^double dt]
  (m/- (m/+ (m/inc (m/* (m/floor (m/* t 4.0)) -0.6666666666666666))
            (blep t dt))
       (m/* 0.3333333333333333 (m/+ (blep-at t 0.25 dt) (blep-at t 0.5 dt) (blep-at t 0.75 dt)))))

(defn- pb-stairs3
  ^double [^double t ^double dt ^double pw]
  (let [y (m/+ (if (m/< t 0.5)
                 (if (m/< (m/* t 2.0) pw) 1.0 0.3333333333333333)
                 (if (m/< (m/dec (m/* t 2.0)) pw) -0.3333333333333333 -1.0))
               (m/* pw -0.6666666666666666)
               0.3333333333333333)
        t2 (m/+ t 0.5)
        t3 (m/- t2 (m/* pw 0.5))
        t4 (m/+ t3 0.5)]
    (m/- (m/+ y (blep t dt))
         (m/* 0.3333333333333333 (m/+ (blep (cycle-frac t2) dt)
                                      (blep (cycle-frac t3) dt)
                                      (blep (cycle-frac t4) dt))))))

(defn- pb-stairs2
  ^double [^double t ^double dt ^double pw]
  (let [y (m/+ (if (m/< (m/dec (m/* t 2.0)) pw)
                 (bool->double (m/< (m/* t 2.0) pw))
                 -1.0)
               (m/- pw)
               0.5)
        t2 (m/+ (m/- t (m/* pw 0.5)) 0.5)
        t3 (m/+ t2 0.5)]
    (m/- (m/+ y (blep t dt))
         (m/* 0.5 (m/+ (blep (cycle-frac t2) dt) (blep (cycle-frac t3) dt))))))

;; hard sync, `fc` is the follower / leader frequency ratio

(defn- pb-hstri
  ^double [^double t ^double dt ^double fc]
  (if (m/> fc 0.0)
    (let [t (cycle-frac (m/+ t 0.5))
          x (m/+ (m/* t fc) 0.75)
          y (m/dec (m/abs (m/dec (m/* 4.0 (m/- (m/floor (m/+ x 0.75)) x)))))
          dc (m// -1.0 fc)
          pw (m/inc (m/* dc 0.5))
          x (m/- fc (m/floor fc))
          scale (m/dec (m/abs (m/dec (m/* x 2.0))))
          high? (m/> x 0.5)
          x (if high? (m/- x 0.5) x)
          dc (if high? dc (m/- dc))
          y (m/+ (m/- y (m/* dc (m/- (m/* 2.0 (m/sq x)) x)))
                 (m/* (blep t dt) scale))
          n (long (m/ceil (m/* fc 2.0)))
          odd? (m/odd? n)
          n (if odd? (m/long-dec n) n)
          t (if odd? (m/+ t pw) t)
          step (m/* fc 4.0 dt)
          step (if odd? (m/- step) step)]
      (loop [i 0
             t t
             y y
             scale step]
        (if (m/< i n)
          (let [t (cycle-frac t)]
            (recur (m/long-inc i) (m/+ t pw) (m/+ y (m/* (blamp t dt) scale)) (m/- scale)))
          y)))
    0.0))

(defn- pb-hssqr
  ^double [^double t ^double dt ^double fc]
  (let [t (cycle-frac (m/+ t 0.5))
        x (m/* t fc)
        y (m/* -2.0 (bool->double (m/< (m/- x (m/floor x)) 0.5)))
        xi (m/floor fc)
        y (m/+ y (if (m/> (m/+ xi 0.5) fc)
                   (m/- 2.0 (m// xi (m/max fc 0.5)))
                   (m// (m/inc xi) fc)))
        n (m/ceil (m/- fc 0.5))
        cnt (long (m/+ n n))
        pw (m/- 1.0 (m// 0.5 fc))
        first? (m/< xi n)
        t (if first? t (m/+ t pw))
        scale (if first? -1.0 1.0)]
    (loop [i 0
           t t
           y y
           scale scale]
      (if (m/< i cnt)
        (let [t (cycle-frac t)]
          (recur (m/long-inc i) (m/+ t pw) (m/+ y (m/* (blep t dt) scale)) (m/- scale)))
        y))))

(defn- pb-hssaw
  ^double [^double t ^double dt ^double fc]
  (let [nd (m/max (m/dec (m/ceil fc)) 0.0)
        n (long nd)
        scale (m/- fc nd)
        pw (if (m/> fc 0.0) (m// scale fc) 1.0)
        t (cycle-frac (m/+ t 0.5))
        y (m/* t fc)
        y (m/+ (m/- (m/* 2.0 (m/- y (m/floor y))) (m/* (m/sq pw) fc)) pw -1.0)
        y (m/- y (m/* (blep t dt) scale))]
    (if (m/> n 0)
      (let [pw (m/inc (m// (m/dec pw) nd))]
        (loop [i 0
               t t
               y y]
          (if (m/< i n)
            (let [t (cycle-frac (m/+ t pw))]
              (recur (m/long-inc i) t (m/- y (blep t dt))))
            y)))
      y)))

;; circle, the only stateful waveform

(defn- pb-circ-body
  "Returns a new body of [[polyblep-circ]] with its own low-pass state, a sample depends on all previous ones.
  The state is private to one sequence of samples, which is realized once and in order."
  []
  (let [lp (double-array 1)]
    (fn ^double [^double t ^double dt ^double _p]
      (let [t2u (m/+ t 0.015625)
            t2 (cycle-frac t2u)
            target (m/+ (blep t2 dt)
                        (m/- (blep-at t2u 0.5 dt))
                        (m/dec (m/* 2.0 (bool->double (m/< t2 0.5)))))
            state (m/+ (aget lp 0)
                       (m/* (m// (m/- target (aget lp 0)) (m/+ dt 0.028134884879909565)) dt))]
        (aset lp 0 state)
        (m/* 0.5 (m/+ state (m/sin (m/* m/TWO_PI t))))))))

;; sample generators

(defn- polyblep-doc
  [wave param doc]
  (str doc ". Port of `poly_" wave "` from Tale's `poly_blep.jsfx-inc` (PolyBLEP: polynomial band-limited step).\n\n"
       "Returns an infinite lazy sequence of samples of the wave band-limited around its jumps, see [[polyblep]] for the options. "
       (case param
         (nil :state) ""
         :pw "Uses `:pw`, the pulse width in [0, 1] (default 0.5). "
         :n "Uses `:n`, the number of steps (integer, at least 1, default 1). "
         :f2 "Uses `:f2`, the follower frequency in Hz (required), the leader frequency is `:f`. ")
       "Silent for a phase step `f / fs` of 0.25 and above, faded out above 0.2."))

(defn- polyblep-param-fn
  "Returns a function of the time and the phase step with the extra parameter of the body."
  ^clojure.lang.IFn$DDD [param {:keys [pw n f2 fs]}]
  (case param
    nil (fn ^double [^double _t ^double _dt] 0.0)
    :pw (let [pwfn (->const-fn (if (nil? pw) 0.5 pw))]
          (fn ^double [^double t ^double _dt]
            (m/min 1.0 (m/max 0.0 (double (pwfn t))))))
    :n (let [n (double (if (nil? n) 1 n))]
         (fn ^double [^double _t ^double _dt] n))
    :f2 (let [f2fn (->const-fn f2)
              fs (double fs)]
          (fn ^double [^double t ^double dt]
            (if (m/> dt 0.0)
              (m// (double (f2fn t)) (m/* dt fs))
              0.0)))))

(defn- check-polyblep-options!
  [param {:keys [fs pw n f2]}]
  (when-not (and (number? fs) (m/> (double fs) 0.0))
    (throw (ex-info "polyblep: :fs must be a positive number" {:fs fs})))
  (case param
    :pw (when (and (number? pw) (or (m/< (double pw) 0.0) (m/> (double pw) 1.0)))
          (throw (ex-info "polyblep: :pw must be in [0, 1]" {:pw pw})))
    :n (when-not (or (nil? n) (and (number? n) (m/>= (double n) 1.0) (== (double n) (m/floor (double n)))))
         (throw (ex-info "polyblep: :n must be an integer, at least 1" {:n n})))
    :f2 (when-not (or (number? f2) (fn? f2))
          (throw (ex-info "polyblep: :f2 (follower frequency in Hz, number or function of time) is required" {:f2 f2})))
    nil))

(defn- with-ignored-parameter
  "Lifts a body `(t, dt)` to the body `(t, dt, p)` ignoring `p`."
  ^clojure.lang.IFn$DDDD [^clojure.lang.IFn$DDD body]
  (fn ^double [^double t ^double dt ^double _p]
    (.invokePrim body t dt)))

(defn- polyblep-samples
  "Lazy sequence of `amplitude fade-gain body(phase, dt, p)` for the time `i / fs`, `i = 0, 1, ...`.

  The phase is the accumulated phase step `f / fs` wrapped to [0, 1) plus the `:phase` offset."
  [^clojure.lang.IFn$DDDD body param opts]
  (let [opts (update opts :fs #(if (nil? %) 44100.0 %))
        _ (check-polyblep-options! param opts)
        {:keys [f phase amplitude fs]} opts
        fs (double fs)
        ffn (->const-fn (if (nil? f) 1.0 f))
        pfn (->const-fn (if (nil? phase) 0.0 phase))
        afn (->const-fn (if (nil? amplitude) 1.0 amplitude))
        paramfn (polyblep-param-fn param opts)]
    (letfn [(samples [^long i ^double acc]
              (lazy-seq
               (let [time (m// (double i) fs)
                     dt (m// (double (ffn time)) fs)
                     ph (cycle-frac (m/+ acc (double (pfn time))))
                     p (.invokePrim paramfn time dt)]
                 (cons (m/* (double (afn time)) (fade-gain dt) (.invokePrim body ph dt p))
                       (samples (m/long-inc i) (cycle-frac (m/+ acc dt)))))))]
      (samples 0 0.0))))

(defmacro polyblep
  "Defines `polyblep-<wave>`, a generator of samples of a PolyBLEP wave.

  Arguments:

  * `wave` - symbol, the name; the wave is `poly_<wave>` in Tale's `poly_blep.jsfx-inc`
  * `param` - the extra parameter of the wave: `nil` (none), `:state` (none, a stateful wave: `fun` is then a function without arguments returning a new function `(t, dt, p)` for every generated sequence), `:pw` (pulse width in [0, 1], default 0.5, number or `(fn [t])`), `:n` (number of steps, integer at least 1, default 1) or `:f2` (follower frequency in Hz for hard sync, required, number or `(fn [t])`; the follower / leader ratio is `:f2 / :f`)
  * `doc` - short description of the wave
  * `fun` - function of the phase in [0, 1) and the phase step `f / fs` (and the extra parameter, when `param` is set) returning the sample of amplitude 1, with primitive `double` arguments; `polyblep-samples` applies the amplitude and the gain

  The generated function takes an options map and returns an infinite lazy sequence of samples.
  The options (see below) are checked when the function is called, `ex-info` is thrown for an invalid one.

  The wave is multiplied by the amplitude and by the gain of Tale's `poly_setdt`: 1 for the phase step `f / fs` up to 0.2, fading to 0 at 0.25 (the wave is silent above).
  A negative frequency gives a negative phase step: the wave is then not corrected.

  Options:

  * `:f` - frequency in Hz (default 1.0), number or `(fn [t])` of time in seconds
  * `:phase` - phase offset in cycles (default 0.0), number or `(fn [t])`
  * `:amplitude` - amplitude (default 1.0), number or `(fn [t])`
  * `:fs` - sample rate in Hz (default 44100.0), positive; the sample `i` is taken at the time `i / fs`
  * `:pw`, `:n`, `:f2` - the extra parameter, when the wave has one"
  [wave param doc fun]
  `(defn ~(symbol (str "polyblep-" wave))
     ~(polyblep-doc wave param doc)
     [opts#]
     (polyblep-samples ~(case param
                          nil `(with-ignored-parameter ~fun)
                          :state `(~fun)
                          fun)
                       ~(when (#{:pw :n :f2} param) param)
                       opts#)))

;; sine family

#_{:clj-kondo/ignore [:unresolved-symbol]}
(do
  (polyblep sinn :n "Stepped sine" pb-sinn)
  (polyblep half nil "Half-wave rectified sine (fixed 0.5 pulse width), anti-aliasing depends on the frequency" pb-half)
  (polyblep half2 :pw "Half-wave rectified sine with variable pulse width, anti-aliasing is worse for a pulse width below 1" pb-half2)
  (polyblep full nil "Full-wave rectified sine" pb-full)
  (polyblep sinp nil "Pulse sine (fixed 0.5 pulse width)" pb-sinp)
  (polyblep alt nil "Alternating sine" pb-alt)
  (polyblep camel nil "Camel sine, anti-aliasing doesn't work well at higher frequencies" pb-camel)
  (polyblep camela nil "Alternating camel sine" pb-camela)
  (polyblep camel2 nil "Bipolar camel sine" pb-camel2)
  (polyblep sin2 nil "Bipolar squared sine" pb-sin2)

  ;; curves
  (polyblep para nil "Parabola" pb-para)
  (polyblep arc nil "Cycloid, anti-aliasing is probably not correct" pb-arc)
  (polyblep circ :state "Circle (approximation): the mix of a sine and a low-pass filtered square, the only wave with a state" pb-circ-body)

  ;; triangle family
  (polyblep tri nil "Triangle" pb-tri)
  (polyblep hyptri nil "Hyper triangular wave" pb-hyptri)
  (polyblep trin :n "Stepped triangle" pb-trin)
  (polyblep tri2 :pw "Modified triangle with variable pulse width" pb-tri2)
  (polyblep trip :pw "Triangular pulse with variable pulse width" pb-trip)
  (polyblep trap nil "Trapezoid (fixed 0.5 pulse width)" pb-trap)
  (polyblep trap2 :pw "Trapezoid with variable pulse width" pb-trap2)

  ;; square family
  (polyblep sqr nil "Square" pb-sqr)
  (polyblep lpsqr2 nil "Low-pass filtered square approximation of order 2" pb-lpsqr2)
  (polyblep lpsqr3 nil "Low-pass filtered square approximation of order 3" pb-lpsqr3)
  (polyblep intlpsqr2 nil "Integrated low-pass filtered square of order 2" pb-intlpsqr2)
  (polyblep intlpsqr3 nil "Integrated low-pass filtered square of order 3" pb-intlpsqr3)
  (polyblep bpsqr2 nil "Band-pass filtered square approximation" pb-bpsqr2)
  (polyblep sqrm1 nil "Square wave minus the fundamental" pb-sqrm1)
  (polyblep sqr2 :pw "Modified square with variable pulse width" pb-sqr2)
  (polyblep rect :pw "Rectangle (pulse) with variable pulse width" pb-rect)

  ;; sawtooth family
  (polyblep saw nil "Sawtooth" pb-saw)
  (polyblep ramp nil "Ramp" pb-ramp)
  (polyblep sawp nil "Saw pulse" pb-sawp)
  (polyblep intsaw nil "Integrated sawtooth" pb-intsaw)
  (polyblep cubsaw nil "Cubic sawtooth" pb-cubsaw)
  (polyblep sinsaw nil "Sine sawtooth" pb-sinsaw)
  (polyblep hpsaw6 nil "High-pass filtered sawtooth approximation" pb-hpsaw6)
  (polyblep sawm1 nil "Sawtooth minus the fundamental" pb-sawm1)
  (polyblep sawn :n "Stepped sawtooth" pb-sawn)
  (polyblep hwsaw :pw "Half-wave rectified sawtooth with variable pulse width" pb-hwsaw)
  (polyblep logit3 nil "Logit approximation" pb-logit3)

  ;; staircases
  (polyblep stairs nil "Staircase (fixed 0.5 pulse width)" pb-stairs)
  (polyblep stairs3 :pw "Staircase with variable pulse width" pb-stairs3)
  (polyblep stairs2 :pw "Uneven staircase with variable pulse width" pb-stairs2)

  ;; hard sync
  (polyblep hstri :f2 "Hard sync triangle: the triangle of the follower frequency `:f2`, reset by the leader `:f`" pb-hstri)
  (polyblep hssqr :f2 "Hard sync square: the square of the follower frequency `:f2`, reset by the leader `:f`" pb-hssqr)
  (polyblep hssaw :f2 "Hard sync sawtooth: the sawtooth of the follower frequency `:f2`, reset by the leader `:f`" pb-hssaw))
