(ns fastmath.signal.chirp
  (:require [fastmath.core :as m]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn ->constant
  ([] (->constant nil))
  ([{:keys [^double f0 ^double phase0]
     :or {f0 0.0 phase0 0.0}}]
   (fn ^double [^double t]
     (m/+ phase0 (m/* f0 t)))))

(defn ->linear
  ([] (->linear nil))
  ([{:keys [^double time ^double f0 ^double f1 ^double phase0]
     :or {time 1.0 f0 0.0 f1 5.0 phase0 0.0}}]
   (fn ^double [^double t]
     (let [c2 (m/* 0.5 (m// (m/- f1 f0) time))]
       (m/+ phase0 (m/+ (m/* c2 t t)
                        (m/* f0 t)))))))

(defn ->quadratic-up
  ([] (->quadratic-up nil))
  ([{:keys [^double time ^double f0 ^double f1 ^double phase0]
     :or {time 1.0 f0 0.0 f1 5.0 phase0 0.0}}]
   (let [v (m/* m/THIRD (m// (m/- f1 f0) (m/sq time)))]
     (fn ^double [^double t]
       (m/+ phase0 (m/+ (m/* f0 t)
                        (m/* v (m/cb t))))))))

(defn ->quadratic-down
  ([] (->quadratic-down nil))
  ([{:keys [^double time ^double f0 ^double f1 ^double phase0]
     :or {time 1.0 f0 0.0 f1 5.0 phase0 0.0}}]
   (let [v (m/* m/THIRD (m// (m/- f1 f0) (m/sq time)))
         time3 (m/cb time)]
     (fn ^double [^double t]
       (m/+ phase0 (m/+ (m/* f1 t)
                        (m/* v (m/- (m/cb (m/- time t)) time3))))))))

(defn ->logarithmic
  ([] (->logarithmic nil))
  ([{:keys [^double time ^double f0 ^double f1 ^double phase0]
     :or {time 1.0 f0 m/EPSILON f1 5.0 phase0 0.0}}]
   (let [f0 (m/max m/EPSILON f0)
         k (m// f1 f0)
         v (m// (m/* time f0) (m/log k))]
     (fn ^double [^double t]
       (m/+ phase0 (m/* v (m/dec (m/pow k (m// t time)))))))))

(defn ->hyperbolic
  ([] (->hyperbolic nil))
  ([{:keys [^double time ^double f0 ^double f1 ^double phase0]
     :or {time 1.0 f0 m/EPSILON f1 5.0 phase0 0.0}}]
   (let [f0 (m/max m/EPSILON f0)
         v1 (m// (m/- f1 f0) (m/* f1 time))
         v2 (m// (m/- f0) v1)]
     (fn ^double [^double t]
       (m/+ phase0 (m/* v2 (m/log (m/- 1.0 (m/* v1 t)))))))))
