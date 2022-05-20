(ns fastmath.matrix-test
  (:require [fastmath.matrix :as sut]
            [clojure.test :as t]
            [fastmath.vector :as v]
            [fastmath.core :as m]))

(t/deftest create-matrix
  (t/are [c r] (= c r)
    (sut/mat2x2 2.0) (sut/->Mat2x2 2.0 2.0 2.0 2.0)
    (sut/mat3x3 2.0) (sut/->Mat3x3 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0)
    (sut/mat4x4 2.0) (sut/->Mat4x4 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0)
    (sut/diagonal [3 4]) (sut/->Mat2x2 3.0 0.0 0.0 4.0)
    (sut/diagonal [3 4 5]) (sut/->Mat3x3 3.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 5.0)
    (sut/diagonal [3 4 5 6]) (sut/->Mat4x4 3.0 0.0 0.0 0.0
                                           0.0 4.0 0.0 0.0
                                           0.0 0.0 5.0 0.0
                                           0.0 0.0 0.0 6.0)
    (sut/eye 2) (sut/->Mat2x2 1.0 0.0 0.0 1.0)
    (sut/eye 3) (sut/->Mat3x3 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0)
    (sut/eye 4) (sut/->Mat4x4 1.0 0.0 0.0 0.0
                              0.0 1.0 0.0 0.0
                              0.0 0.0 1.0 0.0
                              0.0 0.0 0.0 1.0)))

;; https://danceswithcode.net/engineeringnotes/linear_equations/linear_equations.html

(def d22 [2.0  3.0
        5.0 -10.0])
(def d33 [-3.0 2.0 -6.0
        5.0  7.0 -5.0
        1.0  4.0 -2.0])
(def d44 [4.0  1.0  2.0 -3.0
        -3.0 3.0 -1.0  4.0
        -1.0 2.0  5.0  1.0
        5.0  4.0  3.0 -1.0])

(def m22 (apply sut/mat2x2 d22))
(def m33 (apply sut/mat3x3 d33))
(def m44 (apply sut/mat4x4 d44))

(t/deftest solving
  (t/are [A b x] (= (v/approx (sut/solve A b)) x)
    m22 (v/vec2 -16 30) (v/vec2 -2 -4)
    m33 (v/vec3 6 6 8) (v/vec3 -2 3 1)
    m44 (v/vec4 -16 20 -4 -10) (v/vec4 -1 1 -2 3)))

(t/deftest outer-product
  (t/are [v1 v2 m] (= (sut/fmap (sut/outer v1 v2) m/approx) m)
    (v/vec2 1 2) (v/vec2 -3 4) (sut/mat2x2 -3.0 4.0 -6.0 8.0)
    (v/vec3 1 2 -1) (v/vec3 -3 4 -2) (sut/mat3x3 -3 4 -2 -6 8 -4 3 -4 2)
    (v/vec4 1 2 -1 3) (v/vec4 -3 4 -2 1) (sut/mat4x4 -3 4 -2 1 -6 8 -4 2 3 -4 2 -1 -9 12 -6 3)))

(t/deftest cols
  (t/are [m r] (= (sut/cols m) r)
    m22 [(v/vec2 2 5) (v/vec2 3 -10)]
    m33 [(v/vec3 -3 5 1) (v/vec3 2 7 4) (v/vec3 -6 -5 -2)]
    m44 [(v/vec4 4 -3 -1 5) (v/vec4 1 3 2 4) (v/vec4 2 -1 5 3) (v/vec4 -3 4 1 -1)]))

(t/deftest column
  (t/are [m c res] (= (sut/col m c) res)
    m22 0 (v/vec2 2 5) m22 1 (v/vec2 3 -10)
    m33 0 (v/vec3 -3 5 1) m33 1 (v/vec3 2 7 4) m33 2 (v/vec3 -6 -5 -2)
    m44 0 (v/vec4 4 -3 -1 5) m44 1 (v/vec4 1 3 2 4) m44 2 (v/vec4 2 -1 5 3) m44 3(v/vec4 -3 4 1 -1)))

(t/deftest rows
  (t/are [m r] (= (sut/rows m) r)
    m22 [(v/vec2 2 3) (v/vec2 5 -10)]
    m33 [(v/vec3 -3 2 -6) (v/vec3 5 7 -5) (v/vec3 1 4 -2)]
    m44 [(v/vec4 4 1 2 -3) (v/vec4 -3 3 -1 4) (v/vec4 -1 2 5 1) (v/vec4 5 4 3 -1)]))

(t/deftest row
  (t/are [m r res] (= (sut/row m r) res)
    m22 0 (v/vec2 2 3) m22 1 (v/vec2 5 -10)
    m33 0 (v/vec3 -3 2 -6) m33 1 (v/vec3 5 7 -5) m33 2 (v/vec3 1 4 -2)
    m44 0 (v/vec4 4 1 2 -3) m44 1 (v/vec4 -3 3 -1 4) m44 2 (v/vec4 -1 2 5 1) m44 3 (v/vec4 5 4 3 -1)))

(t/deftest array2d
  (t/are [m d s] (= (m/double-double-array->seq (sut/mat->array2d m)) (partition s s d))
    m22 d22 2 m33 d33 3 m44 d44 4))

(t/deftest sizes
  (t/are [m s] (= s (sut/nrow m) (sut/ncol m))
    m22 2 m33 3 m44 4))

(t/deftest symmetry
  (t/are [m s] (= s (boolean (sut/symmetric? m)))
    m22 false m33 false m44 false
    (sut/add m22 (sut/transpose m22)) true
    (sut/add m33 (sut/transpose m33)) true
    (sut/add m44 (sut/transpose m44)) true
    (sut/eye 2) true
    (sut/eye 3) true
    (sut/eye 4) true
    (sut/diagonal [9 1]) true
    (sut/diagonal [9 1 2]) true
    (sut/diagonal [1 2 3 4]) true))

(t/deftest transpose
  (t/are [m1 m2] (and (= (sut/transpose m1) m2)
                      (= m1 (sut/transpose m2)))
    m22 (sut/mat2x2 2.0 5.0 3.0 -10.0)
    m33 (sut/mat3x3 -3 5 1 2 7 4 -6 -5 -2)
    m44 (sut/mat4x4 4 -3 -1 5 1 3 2 4 2 -1 5 3 -3 4 1 -1)))

(t/deftest inversion
  (t/are [m] (= (v/approx (seq (sut/inverse m)) 6)
                (v/approx (flatten (m/double-double-array->seq (.getData (org.apache.commons.math3.linear.MatrixUtils/inverse (sut/mat->RealMatrix m))))) 6))
    m22 m33 m44
    (sut/transpose m22) (sut/transpose m33) (sut/transpose m44))
  (t/are [s m] (= (sut/eye s) (sut/fmap (sut/mulm m (sut/inverse m)) m/approx))
    2 m22 3 m33 4 m44))
