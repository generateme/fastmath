(ns fastmath.interpolation.linear
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.interpolation.common :as ic])
  (:import [fastmath.vector Vec3]
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
#_(set! *warn-on-reflection* true)

(defn- coeffs-array
  [^doubles xs ^doubles ys]
  (let [n (m/dec (alength xs))]
    (loop [i (long 0)
           arr []]
      (if (m/== i n)
        arr
        (let [i+ (m/inc i)
              x1 (Array/aget xs i)
              y1 (Array/aget ys i)]
          (recur (m/inc i) (conj arr (Vec3. x1 y1 (m// (m/- (Array/aget ys i+) y1)
                                                       (m/- (Array/aget xs i+) x1))))))))))

(defn linear
  "1d linear interpolation.

  * `xs` - x coordinate
  * `ys` - y coordinate

  Function extrapolates when `x` is outside a domain."
  [xs ys]
  (let [xs (m/seq->double-array xs)
        coeffs (coeffs-array xs (m/seq->double-array ys))]
    (fn ^double [^double x]
      (let [id (ic/binary-search xs x)
            ^Vec3 segment (coeffs id)]
        (m/+ (.y segment)
             (m/* (m/- x (.x segment)) (.z segment)))))))

(defn bilinear
  [xs ys vss]
  (let [vss (m/seq->double-double-array vss)
        xs (m/seq->double-array xs)
        ys (m/seq->double-array ys)]
    (fn bilinear-interpolation
      (^double [[^double x ^double y]] (bilinear-interpolation x y))
      (^double [^double x ^double y]
       (let [i (ic/binary-search xs x)
             j (ic/binary-search ys y)
             i+ (m/inc i)
             j+ (m/inc j)
             xi (Array/aget xs i)
             yj (Array/aget ys j)
             t (m// (m/- x xi) (m/- (Array/aget xs i+) xi))
             u (m// (m/- y yj) (m/- (Array/aget ys j+) yj))
             t- (m/- 1.0 t)
             u- (m/- 1.0 u)]
         (m/+ (m/* t- u- (Array/aget2d vss i  j))
              (m/* t  u- (Array/aget2d vss i+ j))
              (m/* t- u  (Array/aget2d vss i  j+))
              (m/* t  u  (Array/aget2d vss i+ j+))))))))
