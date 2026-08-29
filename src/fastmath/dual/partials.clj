(ns fastmath.dual.partials
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(defn mult
  ([partials1 partials2] (mult partials1 partials2 1.0 1.0))
  ([partials1 partials2 ^double v1 ^double v2]
   (cond
     (and (seq partials1) (seq partials2)) (v/add (v/mult partials1 v1)
                                                  (v/mult partials2 v2))
     (m/zero? (count partials1)) (v/mult partials2 v2)
     :else (v/mult partials1 v1))))

(defn div
  ([partials1 partials2] (mult partials1 partials2 1.0 -1.0))
  ([partials1 partials2 ^double v1 ^double v2]
   (mult partials1 partials2 (m// v2) (m/- (m// v1 (m/* v2 v2))))))

