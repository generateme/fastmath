(ns fastmath.signal.chirp
  "Phase functions of frequency sweeps (chirps), used by `fastmath.signal/chirp`.

  Every constructor takes an options map and returns a function of time `t` giving the phase in cycles; `sin(2 pi phase(t))` is the chirp, the instantaneous frequency is the derivative of the phase. The sweeps are the ones of `scipy.signal.chirp` (which uses cosine and a phase in degrees):

  * `:time` - time `t1` at which the frequency reaches `:f1` (default 1.0), must be positive
  * `:f0` - frequency at `t = 0`, `:f1` - frequency at `t = t1`
  * `:phase0` - phase offset in cycles (default 0.0)

  Throws `ex-info` for a non-positive `:time` and for non-positive frequencies of logarithmic and hyperbolic sweeps."
  (:require [fastmath.core :as m]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- positive!
  ^double [what ^double v]
  (when-not (m/pos? v)
    (throw (ex-info (str what " should be positive.") {what v})))
  v)

(defn ->constant
  "Constant frequency `:f0` (default 0.0): phase `phase0 + f0 t`."
  ([] (->constant nil))
  ([{:keys [^double f0 ^double phase0]
     :or {f0 0.0 phase0 0.0}}]
   (fn ^double [^double t]
     (m/+ phase0 (m/* f0 t)))))

(defn ->linear
  "Linear sweep, frequency `f0 + (f1 - f0) t / t1`. Defaults: `:f0` 0.0, `:f1` 5.0, `:time` 1.0."
  ([] (->linear nil))
  ([{:keys [^double time ^double f0 ^double f1 ^double phase0]
     :or {time 1.0 f0 0.0 f1 5.0 phase0 0.0}}]
   (let [c2 (m/* 0.5 (m// (m/- f1 f0) (positive! :time time)))]
     (fn ^double [^double t]
       (m/+ phase0 (m/+ (m/* c2 t t)
                        (m/* f0 t)))))))

(defn ->quadratic-up
  "Quadratic sweep with the parabola vertex at `t = 0` (`scipy` `vertex_zero=True`), frequency `f0 + (f1 - f0) t^2 / t1^2`. Defaults: `:f0` 0.0, `:f1` 5.0, `:time` 1.0."
  ([] (->quadratic-up nil))
  ([{:keys [^double time ^double f0 ^double f1 ^double phase0]
     :or {time 1.0 f0 0.0 f1 5.0 phase0 0.0}}]
   (let [v (m/* m/THIRD (m// (m/- f1 f0) (m/sq (positive! :time time))))]
     (fn ^double [^double t]
       (m/+ phase0 (m/+ (m/* f0 t)
                        (m/* v (m/cb t))))))))

(defn ->quadratic-down
  "Quadratic sweep with the parabola vertex at `t = t1` (`scipy` `vertex_zero=False`), frequency `f1 + (f0 - f1) (t1 - t)^2 / t1^2`. Defaults: `:f0` 0.0, `:f1` 5.0, `:time` 1.0."
  ([] (->quadratic-down nil))
  ([{:keys [^double time ^double f0 ^double f1 ^double phase0]
     :or {time 1.0 f0 0.0 f1 5.0 phase0 0.0}}]
   (let [v (m/* m/THIRD (m// (m/- f1 f0) (m/sq (positive! :time time))))
         time3 (m/cb time)]
     (fn ^double [^double t]
       (m/+ phase0 (m/+ (m/* f1 t)
                        (m/* v (m/- (m/cb (m/- time t)) time3))))))))

(defn ->logarithmic
  "Exponential (geometric) sweep, frequency `f0 (f1/f0)^(t/t1)`; constant `f0` when `f0 = f1`. Frequencies must be positive. Defaults: `:f0` 1e-10, `:f1` 5.0, `:time` 1.0."
  ([] (->logarithmic nil))
  ([{:keys [^double time ^double f0 ^double f1 ^double phase0]
     :or {time 1.0 f0 m/EPSILON f1 5.0 phase0 0.0}}]
   (positive! :time time)
   (positive! :f0 f0)
   (positive! :f1 f1)
   (if (m/== f0 f1)
     (->constant {:f0 f0 :phase0 phase0})
     (let [k (m// f1 f0)
           v (m// (m/* time f0) (m/log k))]
       (fn ^double [^double t]
         (m/+ phase0 (m/* v (m/dec (m/pow k (m// t time))))))))))

(defn ->hyperbolic
  "Hyperbolic sweep, frequency `f0 f1 t1 / ((f0 - f1) t + f1 t1)`; constant `f0` when `f0 = f1`. Beyond the singularity `t = f1 t1 / (f1 - f0)` the phase uses the absolute value of the logarithm argument, as scipy. Frequencies must be positive. Defaults: `:f0` 1e-10, `:f1` 5.0, `:time` 1.0."
  ([] (->hyperbolic nil))
  ([{:keys [^double time ^double f0 ^double f1 ^double phase0]
     :or {time 1.0 f0 m/EPSILON f1 5.0 phase0 0.0}}]
   (positive! :time time)
   (positive! :f0 f0)
   (positive! :f1 f1)
   (if (m/== f0 f1)
     (->constant {:f0 f0 :phase0 phase0})
     (let [v1 (m// (m/- f1 f0) (m/* f1 time))
           v2 (m// (m/- f0) v1)]
       (fn ^double [^double t]
         (m/+ phase0 (m/* v2 (m/log (m/abs (m/- 1.0 (m/* v1 t)))))))))))
