(ns fastmath.fields.u
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn unpolar
  ([] {:type :regular})
  ([^double amount _]
   (let [vvar2 (m// (m/* amount 0.5) m/PI)]
     (fn [^Vec2 v]
       (let [r (m/* vvar2 (m/exp (.y v)))
             s (m/sin (.x v))
             c (m/cos (.x v))]
         (Vec2. (m/* r s) (m/* r c)))))))

