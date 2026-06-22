(ns fastmath.calculus.divided
  (:require [fastmath.core :as m])
  (:import [fastmath.java Array]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(defn divided-all
  "Divided differences f[x0, x1, ..., xn], returns every level as a sequence of values."
  [f xs]
  (let [ys (map f xs)
        axs (m/seq->double-array xs)
        n (alength axs)]
    (reduce (fn [buff ^long level]
              (conj buff (->> (first buff)
                              (partition 2 1)
                              (map-indexed (fn [^long id [^double y0 ^double y1]]
                                             (m// (m/- y1 y0)
                                                  (m/- (Array/aget axs (m/+ id level))
                                                       (Array/aget axs id)))))))) (list ys) (range 1 n))))

(defn divided
  "Divided differences f[x0, x1, ..., xn], returns every level as a sequence of values."
  ^double [f xs]
  (let [ys (double-array (map f xs))
        axs (m/seq->double-array xs)
        n (alength axs)]
    (doseq [^long i (range 1 n)
            ^long j (range (m/- n i))]
      (Array/aset ys j (m// (m/- (Array/aget ys (m/inc j))
                                 (Array/aget ys j))
                            (m/- (Array/aget axs (m/+ j i))
                                 (Array/aget axs j)))))
    (Array/aget ys 0)))
