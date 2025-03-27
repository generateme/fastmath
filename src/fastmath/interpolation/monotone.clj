(ns fastmath.interpolation.monotone
  (:require [fastmath.core :as m]
            [fastmath.interpolation.common :as ic])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)
(m/use-primitive-operators)

(defn monotone
  "Monotone interpolation

  https://gist.github.com/lecho/7627739"
  [xs ys]
  (let [n (count xs)
        d (double-array (map (fn [[^double x1 ^double x2] [^double y1 ^double y2]]
                               (/ (- y2 y1) (- x2 x1))) (partition 2 1 xs) (partition 2 1 ys)))
        m (double-array (conj (map (fn [[^double d1 ^double d2]]
                                     (if d2 (* 0.5 (+ d1 d2)) d1)) (partition-all 2 1 d)) (first d)))
        stop (- n 2)
        m (m/seq->double-array (loop [idx (int 0)
                                      mi (Array/aget m 0)
                                      mi+ (Array/aget m 1)
                                      ms []]
                                 (if (== idx stop)
                                   (conj ms mi mi+)
                                   (let [di (Array/aget d idx)
                                         mi++ (Array/aget m (+ idx 2))]
                                     (if (zero? di)
                                       (recur (inc idx) 0.0 mi++ (conj ms 0))
                                       (let [a (/ mi di)
                                             b (/ mi+ di)
                                             h (m/hypot-sqrt a b)]
                                         (if (> h 3.0)
                                           (let [t (/ 3.0 h)]
                                             (recur (inc idx) (* t b di) mi++ (conj ms (* t a di))))
                                           (recur (inc idx) mi+ mi++ (conj ms mi)))))))))
        xs (m/seq->double-array xs)
        ys (m/seq->double-array ys)]
    (fn ^double [^double v]
      (let [b (java.util.Arrays/binarySearch xs v)]
        (if-not (neg? b)
          (Array/aget ys b)
          (let [i (ic/binary-search-id b n)
                mxi (Array/aget xs i)
                h (- (Array/aget xs (inc i)) mxi)
                t (/ (- v mxi) h)]
            (+ (* (m/sq (- 1.0 t))
                  (+ (* h t (Array/aget m i))
                     (* (inc (* 2.0 t)) (Array/aget ys i))))
               (* t t (+ (* h (dec t) (Array/aget m (inc i)))
                         (* (- 3.0 t t) (Array/aget ys (inc i))))))))))))

(m/unuse-primitive-operators)

