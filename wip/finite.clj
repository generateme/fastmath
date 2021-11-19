(ns fastmath.finite
  "Finite differentiation"
  (:require [fastmath.core :as m])
  (:import [org.apache.commons.math3.linear Array2DRowRealMatrix ArrayRealVector LUDecomposition]
           [org.apache.commons.math3.util CombinatoricsUtils]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;; coeffs

(defn central-coeffs
  [^long derivative ^long accuracy]
  (assert (pos? accuracy) "Accuracy should be positive number")
  (let [p (int (m/ceil (dec (+ (m/floor (* 0.5 (inc derivative))) accuracy))))
        rp (range (- p) (inc p))
        cp (count rp)
        B (-> (map-indexed (fn [idx _]
                             (if (= idx derivative) (CombinatoricsUtils/factorialDouble derivative) 0.0)) (range cp))
              (m/seq->double-array)
              (ArrayRealVector.))
        A (-> (take cp (iterate (fn [xs]
                                  (map (fn [^long a ^long b] (* a b)) rp xs)) (repeat cp 1)))
              (m/seq->double-double-array)
              (Array2DRowRealMatrix.))]
    [rp (-> A
            (LUDecomposition.)
            (.getSolver)
            (.solve B)
            (.toArray))]))

#_(defmacro make-diff
    [derivative accuracy]
    (let [[rp coeffs] (central-coeffs derivative accuracy)]
      `(fn fname#
         ([f# ^double step#]
          (let [step# (m/pow step# ~derivative)]
            (fn [x#]
              (/ ))))
         ([f#] (fname# f# 0.01)))))

(seq (central-coeffs 1 1))
