;; Numerical recipes, p. 128
(ns fastmath.interpolation.barycentric
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [fastmath.vector Vec2]
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
#_(set! *warn-on-reflection* true)

(defn- coeffs-array
  ^doubles [^doubles xs ^long d ^long n]
  (let [n-d (m/- n d)
        arr (double-array n)]
    (dotimes [k n]
      (let [imin (m/max 0 (m/- k d))
            imax (if (m/>= k n-d) n-d (m/inc k))
            xxk (Array/aget xs k)
            ^Vec2 res (->> (range imin imax)
                           (reduce (fn [^Vec2 b ^long i]
                                     (let [jmax (m/min (m/+ i d 1) n)
                                           ^double term (reduce (fn [^double t ^long j]
                                                                  (if (m/== j k)
                                                                    t
                                                                    (m/* t (m/- xxk (Array/aget xs j)))))
                                                                1.0 (range i jmax))
                                           term (m// (.y b) term)]
                                       (Vec2. (m/+ (.x b) term) (m/- (.y b)))))
                                   (Vec2. 0.0 (if (m/odd? imin) -1.0 1.0))))]
        (Array/aset arr k (.x res))))
    arr))

(defn barycentric
  ([xs ys] (barycentric xs ys nil))
  ([xs ys {:keys [^long order]
           :or {order 1}}]
   (assert (and (m/not-neg? order)
                (m/< order (count xs)))
           "Degree should be non-negative and lower than count of xs.")
   (let [xs (m/seq->double-array xs)
         ys (m/seq->double-array ys)
         n (alength xs)
         ^doubles coeffs (coeffs-array xs order n)]
     (fn ^double [^double x]
       (let [bsid (java.util.Arrays/binarySearch xs x)]
         (if (m/not-neg? bsid)
           (Array/aget ys bsid)
           (loop [i (int 0)
                  ^Vec2 nd (Vec2. 0.0 0.0)]
             (if (m/== i n)
               (m// (.x nd) (.y nd))
               (let [temp (m// (Array/aget coeffs i)
                               (m/- x (Array/aget xs i)))]
                 (recur (m/inc i) (v/add nd (Vec2. (m/* temp (Array/aget ys i)) temp))))))))))))
