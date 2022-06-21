(ns fastmath.fields.u
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn unpolar
  ([] {:type :regular})
  ([^double amount _]
   (let [vvar2 (/ (* amount 0.5) m/PI)]
     (fn [^Vec2 v]
       (let [r (* vvar2 (m/exp (.y v)))
             s (m/sin (.x v))
             c (m/cos (.x v))]
         (Vec2. (* r s) (* r c)))))))
