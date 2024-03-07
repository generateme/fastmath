;; Numerical recipes p.123
(ns fastmath.interpolation.cubic
  (:require [fastmath.core :as m]
            [fastmath.interpolation.common :as ic])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(deftype Coeffs [^double x1 ^double x2
                 ^double y1 ^double y2
                 ^double y21 ^double y22
                 ^double h ^double h2])

;; not optimal, but calculated once
(defn- y''
  [^doubles xs ^doubles ys]
  (let [n (count xs)
        n- (m/dec n)
        u (double-array (repeat n- 0.0))
        y2 (double-array (repeat n 0.0))]
    (dotimes [iter (m/dec n-)]
      (let [i (inc iter)
            i- (m/dec i)
            i+ (m/inc i)
            x- (Array/aget xs i-)
            x (Array/aget xs i)
            x+ (Array/aget xs i+)
            sig (m// (m/- x x-) (m/- x+ x-))
            p (m/+ (m/* sig (Array/aget y2 i-)) 2.0)
            y2i (m// (m/dec sig) p)
            y (Array/aget ys i)
            ui (m/- (m// (m/- (Array/aget ys i+) y)
                         (m/- x+ x))
                    (m// (m/- y (Array/aget ys i-))
                         (m/- x x-)))
            ui (m// (m/- (m// (m/* 6.0 ui)
                              (m/- x+ x-))
                         (m/* sig (Array/aget u i-))) p)]
        (Array/aset u i ui)
        (Array/aset y2 i y2i)))
    (loop [k (m/dec n-)
           y2 y2]
      (if (m/neg? k)
        y2
        (recur (m/dec k)
               (Array/aset y2 k (m/+ (m/* (Array/aget y2 k)
                                          (Array/aget y2 (m/inc k ))) (Array/aget u k))))))))

(defn- coeffs-array
  [^doubles xs ^doubles ys]
  (let [n (m/dec (alength xs))
        y2 (y'' xs ys)]
    (loop [i (long 0)
           arr []]
      (if (m/== i n)
        arr
        (let [i+ (m/inc i)
              x1 (Array/aget xs i)
              x2 (Array/aget xs i+)
              h (m/- x2 x1)]
          (recur (m/inc i)
                 (conj arr (Coeffs. x1 x2
                                    (Array/aget ys i) (Array/aget ys i+)
                                    (Array/aget y2 i) (Array/aget y2 i+)
                                    h (m// (* h h) 6.0)))))))))

(defn cubic
  [xs ys]
  (let [xs (m/seq->double-array xs)
        coeffs (coeffs-array xs (m/seq->double-array ys))]
    (fn ^double [^double x]
      (let [id (ic/binary-search xs x)
            ^Coeffs segment (coeffs id)
            a (m// (m/- (.x2 segment) x) (.h segment))
            b (m// (m/- x (.x1 segment)) (.h segment))]
        (m/+ (m/* a (.y1 segment))
             (m/* b (.y2 segment))
             (m/* (m/+ (m/* (m/- (m/* a a a) a)
                            (.y21 segment))
                       (m/* (m/- (m/* b b b) b)
                            (.y22 segment))) (.h2 segment)))))))

(defn cubic-2d
  [xs ys vss]
  (let [cbs (map (partial cubic (m/seq->double-array ys)) vss)
        xs (m/seq->double-array xs)]
    (fn cubic-2d-interpolation
      (^double [[^double x ^double y]] (cubic-2d-interpolation x y))
      (^double [^double x ^double y]
       (let [nys (map (fn [f] (f y)) cbs)
             interpolator (cubic xs nys)]
         (interpolator x))))))
