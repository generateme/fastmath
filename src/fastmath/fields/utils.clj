(ns fastmath.fields.utils
  (:require [fastmath.random :as r]
            [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec2]))

#_(set! *warn-on-reflection* true)

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(def unitx (Vec2. 1.0 0.0))
(def zerov (Vec2. 0.0 0.0))

(defn sdrand
  "Symetric random from [-mx -mn] and [mn mx]"
  ^double  [^double mn ^double mx]
  (let [rand (r/drand mn mx)]
    (r/randval rand (* -1.0 rand))))

(defn sirand
  "Symmetric irand"
  (^long [^long mx] (sirand 1 mx))
  (^long [^long mn ^long mx]
   (let [rand (r/irand mn mx)]
     (r/randval rand (* -1 rand)))))

;;

(defn closest
  ^long [P ^long n U]
  (loop [i (long 0)
         d2min Double/MAX_VALUE
         j (long 0)]
    (if (< i n)
      (let [d2 (v/dist-sq (P i) U)
            low? (< d2 d2min)]
        (recur (inc i) (if low? d2 d2min) (if low? i j)))
      j)))

(defn vratio
  ^double [^Vec2 P ^Vec2 Q ^Vec2 U]
  (let [PmQ (v/sub P Q)]
    (if (v/is-zero? PmQ)
      1.0
      (-> (v/sub U Q)
          (v/emult PmQ)
          (v/sum)
          (/ (v/magsq PmQ))
          (* 2.0)))))

(defn voronoi
  ^double [P ^long n ^long q U]
  (let [Pq (P q)]
    (loop [i (long 0)
           ratiomax Double/MIN_VALUE]
      (if (< i n)
        (if (== i q)
          (recur (inc i) ratiomax)
          (recur (inc i) (max ratiomax (vratio (P i) Pq U))))
        ratiomax))))
