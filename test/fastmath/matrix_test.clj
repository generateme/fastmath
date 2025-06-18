(ns fastmath.matrix-test
  (:require [fastmath.matrix :as sut]
            [clojure.test :as t]
            [fastmath.vector :as v]
            [fastmath.core :as m]
            [fastmath.random :as r]))

(t/deftest create-matrix
  (t/are [c r] (= c r)
    (sut/mat2x2 2.0) (sut/->Mat2x2 2.0 2.0 2.0 2.0)
    (sut/mat3x3 2.0) (sut/->Mat3x3 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0)
    (sut/mat4x4 2.0) (sut/->Mat4x4 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0)
    (sut/diagonal 3 4) (sut/->Mat2x2 3.0 0.0 0.0 4.0)
    (sut/diagonal 3 4 5) (sut/->Mat3x3 3.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 5.0)
    (sut/diagonal 3 4 5 6) (sut/->Mat4x4 3.0 0.0 0.0 0.0
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
(def m44a (sut/mat->array2d m44))
(def m44ra (sut/mat->RealMatrix m44))

(defn- creator
  [m]
  (condp instance? m
    fastmath.matrix.Mat2x2 sut/mat2x2
    fastmath.matrix.Mat3x3 sut/mat3x3
    fastmath.matrix.Mat4x4 sut/mat4x4))


(t/deftest solving
  (t/are [A b x] (= (v/approx (sut/solve A b)) x)
    m22 (v/vec2 -16 30) (v/vec2 -2.0 -4.0)
    m33 (v/vec3 6 6 8) (v/vec3 -2.0 3.0 1.0)
    m44 (v/vec4 -16 20 -4 -10) (v/vec4 -1.0 1.0 -2.0 3.0))
  (t/are [A b x] (= (seq (v/approx (sut/solve A b))) (seq x))
    m44a (double-array [-16 20 -4 -10]) [-1.0 1.0 -2.0 3.0])
  (t/are [A b x] (= (seq (.getDataRef (v/approx (sut/solve A b)))) (seq x))
    m44ra (v/vec->RealVector [-16 20 -4 -10]) [-1.0 1.0 -2.0 3.0]))

(t/deftest outer-product
  (t/are [v1 v2 m] (= (sut/fmap (sut/outer v1 v2) m/approx) m)
    (v/vec2 1 2) (v/vec2 -3 4) (sut/mat2x2 -3.0 4.0 -6.0 8.0)
    (v/vec3 1 2 -1) (v/vec3 -3 4 -2) (sut/mat3x3 -3 4 -2 -6 8 -4 3 -4 2)
    (v/vec4 1 2 -1 3) (v/vec4 -3 4 -2 1) (sut/mat4x4 -3 4 -2 1 -6 8 -4 2 3 -4 2 -1 -9 12 -6 3))
  (t/is (= (m/double-double-array->seq (sut/mat->array2d (sut/outer [1 2 3 4 5] [4 3 2 2 1])))
           [[4.0 3.0 2.0 2.0 1.0]
            [8.0 6.0 4.0 4.0 2.0]
            [12.0 9.0 6.0 6.0 3.0]
            [16.0 12.0 8.0 8.0 4.0]
            [20.0 15.0 10.0 10.0 5.0]]))
  (t/is (= (sut/outer (v/vec->RealVector [1 2])
                      (v/vec->RealVector [-3 4]))
           (sut/mat->RealMatrix (sut/mat2x2 -3.0 4.0 -6.0 8.0)))))

(t/deftest cols
  (t/are [m r] (= (sut/cols m) r)
    m22 [(v/vec2 2 5) (v/vec2 3 -10)]
    m33 [(v/vec3 -3 5 1) (v/vec3 2 7 4) (v/vec3 -6 -5 -2)]
    m44 [(v/vec4 4 -3 -1 5) (v/vec4 1 3 2 4) (v/vec4 2 -1 5 3) (v/vec4 -3 4 1 -1)])
  (t/is (= (map seq (sut/cols m44a)) '((4.0 -3.0 -1.0 5.0) (1.0 3.0 2.0 4.0)
                                       (2.0 -1.0 5.0 3.0) (-3.0 4.0 1.0 -1.0))))
  (t/is (= (sut/cols m44ra)
           (map v/vec->RealVector '((4.0 -3.0 -1.0 5.0) (1.0 3.0 2.0 4.0)
                                    (2.0 -1.0 5.0 3.0) (-3.0 4.0 1.0 -1.0))))))

(t/deftest column
  (t/are [m c res] (= (sut/col m c) res)
    m22 0 (v/vec2 2 5) m22 1 (v/vec2 3 -10)
    m33 0 (v/vec3 -3 5 1) m33 1 (v/vec3 2 7 4) m33 2 (v/vec3 -6 -5 -2)
    m44 0 (v/vec4 4 -3 -1 5) m44 1 (v/vec4 1 3 2 4) m44 2 (v/vec4 2 -1 5 3) m44 3(v/vec4 -3 4 1 -1))
  (t/is (= (map (comp seq (partial sut/col m44a)) (range 4))
           (map seq [(v/vec4 4 -3 -1 5) (v/vec4 1 3 2 4) (v/vec4 2 -1 5 3) (v/vec4 -3 4 1 -1)])))
  (t/is (= (map (partial sut/col m44ra) (range 4))
           (map v/vec->RealVector [(v/vec4 4 -3 -1 5) (v/vec4 1 3 2 4)
                                   (v/vec4 2 -1 5 3) (v/vec4 -3 4 1 -1)]))))

(t/deftest rows
  (t/are [m r] (= (sut/rows m) r)
    m22 [(v/vec2 2 3) (v/vec2 5 -10)]
    m33 [(v/vec3 -3 2 -6) (v/vec3 5 7 -5) (v/vec3 1 4 -2)]
    m44 [(v/vec4 4 1 2 -3) (v/vec4 -3 3 -1 4) (v/vec4 -1 2 5 1) (v/vec4 5 4 3 -1)])
  (t/is (= (map seq (sut/rows m44a)) '((4.0 1.0 2.0 -3.0) (-3.0 3.0 -1.0 4.0)
                                       (-1.0 2.0 5.0 1.0) (5.0 4.0 3.0 -1.0))))
  (t/is (= (sut/rows m44ra)
           (map v/vec->RealVector '((4.0 1.0 2.0 -3.0) (-3.0 3.0 -1.0 4.0)
                                    (-1.0 2.0 5.0 1.0) (5.0 4.0 3.0 -1.0))))))

(t/deftest row
  (t/are [m r res] (= (sut/row m r) res)
    m22 0 (v/vec2 2 3) m22 1 (v/vec2 5 -10)
    m33 0 (v/vec3 -3 2 -6) m33 1 (v/vec3 5 7 -5) m33 2 (v/vec3 1 4 -2)
    m44 0 (v/vec4 4 1 2 -3) m44 1 (v/vec4 -3 3 -1 4) m44 2 (v/vec4 -1 2 5 1) m44 3 (v/vec4 5 4 3 -1))
  (t/is (= (map (comp seq (partial sut/row m44a)) (range 4))
           (map seq [(v/vec4 4 1 2 -3) (v/vec4 -3 3 -1 4) (v/vec4 -1 2 5 1) (v/vec4 5 4 3 -1)])))
  (t/is (= (map (partial sut/row m44ra) (range 4))
           (map v/vec->RealVector [(v/vec4 4 1 2 -3) (v/vec4 -3 3 -1 4)
                                   (v/vec4 -1 2 5 1) (v/vec4 5 4 3 -1)]))))

(t/deftest array2d
  (t/are [m d s] (= (m/double-double-array->seq (sut/mat->array2d m)) (partition s s d))
    m22 d22 2 m33 d33 3 m44 d44 4 m44a d44 4 m44ra d44 4))

(t/deftest sizes
  (t/are [m s] (= s (sut/nrow m) (sut/ncol m))
    m22 2 m33 3 m44 4 m44a 4 m44ra 4))

(t/deftest symmetry
  (t/are [m s] (= s (boolean (sut/symmetric? m)))
    m22 false m33 false m44 false m44a false m44ra false
    (sut/add m22 (sut/transpose m22)) true
    (sut/add m33 (sut/transpose m33)) true
    (sut/add m44 (sut/transpose m44)) true
    (sut/add m44a (sut/transpose m44a)) true
    (sut/add m44ra (sut/transpose m44ra)) true
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
    m44 (sut/mat4x4 4 -3 -1 5 1 3 2 4 2 -1 5 3 -3 4 1 -1)
    m44ra (sut/mat->RealMatrix (sut/mat4x4 4 -3 -1 5 1 3 2 4 2 -1 5 3 -3 4 1 -1)))
  (t/are [m1 m2] (and (= (m/double-double-array->seq (sut/transpose m1))
                         (m/double-double-array->seq m2))
                      (= (m/double-double-array->seq m1)
                         (m/double-double-array->seq (sut/transpose m2))))
    m44a (sut/mat->array2d (sut/mat4x4 4 -3 -1 5 1 3 2 4 2 -1 5 3 -3 4 1 -1))))

(t/deftest inversion
  (t/are [m] (= (v/approx (seq (sut/mat->array (sut/inverse m))) 6)
                (v/approx (flatten (m/double-double-array->seq (.getData (org.apache.commons.math3.linear.MatrixUtils/inverse (sut/mat->RealMatrix m))))) 6))
    m22 m33 m44 m44a m44ra
    (sut/transpose m22) (sut/transpose m33) (sut/transpose m44)
    (sut/transpose m44a) (sut/transpose m44ra))
  (t/are [s m] (= (seq (sut/eye s))
                  (seq (sut/mat->array (sut/fmap (sut/mulm m (sut/inverse m)) m/approx))))
    2 m22 3 m33 4 m44 4 m44a 4 m44ra))

(t/deftest diag-and-trace
  (t/are [m v s] (and (= (v/vec->Vec (sut/diag m)) v)
                      (= (sut/trace m) s))
    m22 (v/vec2 2.0 -10.0) -8.0
    m33 (v/vec3 -3.0 7.0 -2.0) 2.0
    m44 (v/vec4 4.0 3.0 5.0 -1.0) 11.0
    m44a [4.0 3.0 5.0 -1.0] 11.0
    m44ra [4.0 3.0 5.0 -1.0] 11.0))

(t/deftest determinant
  (t/are [m] (m/approx-eq
              (sut/det m)
              (.getDeterminant (org.apache.commons.math3.linear.LUDecomposition. (sut/mat->RealMatrix m)))
              8)
    m22 m33 m44 m44a m44ra))

(t/deftest add-a-scalar
  (t/are [m res s] (= res (sut/adds m s))
    m22 (apply sut/mat2x2 (map inc d22)) 1.0
    m33 (apply sut/mat3x3 (map dec d33)) -1.0
    m44 (apply sut/mat4x4 (map inc d44)) 1.0
    m44ra (sut/mat->RealMatrix (apply sut/mat4x4 (map inc d44))) 1.0)
  (t/is (= (seq (sut/mat->array (sut/adds m44a 1.0)))
           (seq (sut/mat->array (sut/mat->array2d (apply sut/mat4x4 (map inc d44))))))))

(t/deftest add-sub
  (t/are [m s] (and (= s (sut/add m m))
                    (= ((creator m) 0.0) (sut/sub m m)))
    m22 (apply sut/mat2x2 (v/add d22 d22))
    m33 (apply sut/mat3x3 (v/add d33 d33))
    m44 (apply sut/mat4x4 (v/add d44 d44)))
  (t/are [m s] (and (= s (seq (sut/mat->array (sut/add m m))))
                    (= (repeat 16 0.0) (seq (sut/mat->array (sut/sub m m)))))
    m44a (v/add d44 d44)
    m44ra (v/add d44 d44)))

(t/deftest negation
  (t/are [m] (= ((creator m) 0.0) (sut/add (sut/sub m) m) (sut/add (sut/negate m) m))
    m22 m33 m44)
  (t/are [m] (= (repeat 16 0.0)
                (seq (sut/mat->array (sut/add (sut/sub m) m)))
                (seq (sut/mat->array (sut/add (sut/negate m) m))))
    m44a m44ra))

(t/deftest multiplication
  (t/are [m t1 t2 r] (= r (seq (sut/mat->array (sut/mulm m t1 m t2))))
    m22 false false [19.0 -24.0 -40.0 115.0]
    m22 true false [29.0 -44.0 -44.0 109.0]
    m22 false true [13.0 -20.0 -20.0 125.0]
    m22 true true [19.0 -40.0 -24.0 115.0]
    m33 false false [13.0 -16.0 20.0 15.0 39.0 -55.0 15.0 22.0 -22.0]
    m33 true false [35.0 33.0 -9.0 33.0 69.0 -55.0 -9.0 -55.0 65.0]
    m33 false true [49.0 29.0 17.0 29.0 99.0 43.0 17.0 43.0 21.0]
    m33 true true [13.0 15.0 15.0 -16.0 39.0 22.0 20.0 -55.0 -22.0]
    m44 false false [-4.0 -1.0 8.0 -3.0 0.0 20.0 -2.0 16.0 -10.0 19.0 24.0 15.0 0.0 19.0 18.0 5.0]
    m44 true false [51.0 13.0 21.0 -30.0 13.0 30.0 21.0 7.0 21.0 21.0 39.0 -8.0 -30.0 7.0 -8.0 27.0]
    m44 false true [30.0 -23.0 5.0 33.0 -23.0 35.0 8.0 -10.0 5.0 8.0 31.0 17.0 33.0 -10.0 17.0 51.0]
    m44 true true [-4.0 0.0 -10.0 0.0 -1.0 20.0 19.0 19.0 8.0 -2.0 24.0 18.0 -3.0 16.0 15.0 5.0]
    m44a false false [-4.0 -1.0 8.0 -3.0 0.0 20.0 -2.0 16.0 -10.0 19.0 24.0 15.0 0.0 19.0 18.0 5.0]
    m44a true false [51.0 13.0 21.0 -30.0 13.0 30.0 21.0 7.0 21.0 21.0 39.0 -8.0 -30.0 7.0 -8.0 27.0]
    m44a false true [30.0 -23.0 5.0 33.0 -23.0 35.0 8.0 -10.0 5.0 8.0 31.0 17.0 33.0 -10.0 17.0 51.0]
    m44a true true [-4.0 0.0 -10.0 0.0 -1.0 20.0 19.0 19.0 8.0 -2.0 24.0 18.0 -3.0 16.0 15.0 5.0]
    m44ra false false [-4.0 -1.0 8.0 -3.0 0.0 20.0 -2.0 16.0 -10.0 19.0 24.0 15.0 0.0 19.0 18.0 5.0]
    m44ra true false [51.0 13.0 21.0 -30.0 13.0 30.0 21.0 7.0 21.0 21.0 39.0 -8.0 -30.0 7.0 -8.0 27.0]
    m44ra false true [30.0 -23.0 5.0 33.0 -23.0 35.0 8.0 -10.0 5.0 8.0 31.0 17.0 33.0 -10.0 17.0 51.0]
    m44ra true true [-4.0 0.0 -10.0 0.0 -1.0 20.0 19.0 19.0 8.0 -2.0 24.0 18.0 -3.0 16.0 15.0 5.0])
  (t/are [m r] (= r (seq (sut/mat->array (sut/mulmt m m))))
    m22 [13.0 -20.0 -20.0 125.0]
    m33 [49.0 29.0 17.0 29.0 99.0 43.0 17.0 43.0 21.0]
    m44 [30.0 -23.0 5.0 33.0 -23.0 35.0 8.0 -10.0 5.0 8.0 31.0 17.0 33.0 -10.0 17.0 51.0]
    m44a [30.0 -23.0 5.0 33.0 -23.0 35.0 8.0 -10.0 5.0 8.0 31.0 17.0 33.0 -10.0 17.0 51.0]
    m44ra [30.0 -23.0 5.0 33.0 -23.0 35.0 8.0 -10.0 5.0 8.0 31.0 17.0 33.0 -10.0 17.0 51.0])
  (t/are [m r] (= r (seq (sut/mat->array (sut/tmulm m m))))
    m22 [29.0 -44.0 -44.0 109.0]
    m33 [35.0 33.0 -9.0 33.0 69.0 -55.0 -9.0 -55.0 65.0]
    m44 [51.0 13.0 21.0 -30.0 13.0 30.0 21.0 7.0 21.0 21.0 39.0 -8.0 -30.0 7.0 -8.0 27.0]
    m44a [51.0 13.0 21.0 -30.0 13.0 30.0 21.0 7.0 21.0 21.0 39.0 -8.0 -30.0 7.0 -8.0 27.0]
    m44ra [51.0 13.0 21.0 -30.0 13.0 30.0 21.0 7.0 21.0 21.0 39.0 -8.0 -30.0 7.0 -8.0 27.0])
  (t/are [m r] (= r (seq (sut/mat->array (sut/tmulmt m m))))
    m22 [19.0 -40.0 -24.0 115.0]
    m33 [13.0 15.0 15.0 -16.0 39.0 22.0 20.0 -55.0 -22.0]
    m44 [-4.0 0.0 -10.0 0.0 -1.0 20.0 19.0 19.0 8.0 -2.0 24.0 18.0 -3.0 16.0 15.0 5.0]
    m44a [-4.0 0.0 -10.0 0.0 -1.0 20.0 19.0 19.0 8.0 -2.0 24.0 18.0 -3.0 16.0 15.0 5.0]
    m44ra [-4.0 0.0 -10.0 0.0 -1.0 20.0 19.0 19.0 8.0 -2.0 24.0 18.0 -3.0 16.0 15.0 5.0])
  (t/are [m] (= (seq (sut/mat->array (sut/emulm m m))) (seq (sut/mat->array (sut/sq m))))
    m22 m33 m44 m44a m44ra)
  (t/are [m d] (let [r (r/drand)]
                 (= (seq (sut/mat->array (sut/muls m r)))
                    (v/mult d r)))
    m22 d22 m33 d33 m44 d44 m44a d44 m44ra d44)
  (t/is (= (v/vec2 (v/sum (sut/row m22 0))
                   (v/sum (sut/row m22 1)))
           (sut/mulv m22 (v/vec2 1 1))))
  (t/is (= (v/vec3 (v/sum (sut/row m33 0))
                   (v/sum (sut/row m33 1))
                   (v/sum (sut/row m33 2)))
           (sut/mulv m33 (v/vec3 1 1 1))))
  (t/is (= (v/vec4 (v/sum (sut/row m44 0))
                   (v/sum (sut/row m44 1))
                   (v/sum (sut/row m44 2))
                   (v/sum (sut/row m44 3)))
           (sut/mulv m44 (v/vec4 1 1 1 1))))
  (t/is (= [(v/sum (sut/row m44a 0))
            (v/sum (sut/row m44a 1))
            (v/sum (sut/row m44a 2))
            (v/sum (sut/row m44a 3))]
           (seq (sut/mulv m44a (double-array [1 1 1 1])))))
  (t/is (= [(v/sum (sut/row m44ra 0))
            (v/sum (sut/row m44ra 1))
            (v/sum (sut/row m44ra 2))
            (v/sum (sut/row m44ra 3))]
           (v/vec->Vec (sut/mulv m44ra (v/vec->RealVector [1 1 1 1])))))
  (t/is (= (v/vec2 (v/sum (sut/col m22 0))
                   (v/sum (sut/col m22 1)))
           (sut/vtmul m22 (v/vec2 1 1))))
  (t/is (= (v/vec3 (v/sum (sut/col m33 0))
                   (v/sum (sut/col m33 1))
                   (v/sum (sut/col m33 2)))
           (sut/vtmul m33 (v/vec3 1 1 1))))
  (t/is (= (v/vec4 (v/sum (sut/col m44 0))
                   (v/sum (sut/col m44 1))
                   (v/sum (sut/col m44 2))
                   (v/sum (sut/col m44 3)))
           (sut/vtmul m44 (v/vec4 1 1 1 1))))
  (t/is (= [(v/sum (sut/col m44a 0))
            (v/sum (sut/col m44a 1))
            (v/sum (sut/col m44a 2))
            (v/sum (sut/col m44a 3))]
           (seq (sut/vtmul m44a (double-array [1 1 1 1])))))
  (t/is (= [(v/sum (sut/col m44ra 0))
            (v/sum (sut/col m44ra 1))
            (v/sum (sut/col m44ra 2))
            (v/sum (sut/col m44ra 3))]
           (v/vec->Vec (sut/vtmul m44ra (v/vec->RealVector [1 1 1 1]))))))

;;

(t/deftest cholesky
  (t/are [t m] (v/delta-eq t (seq (sut/cholesky (apply sut/mat m) true)))
    [1.414213 2.121320
     0.000000 2.345208] [2 3 3 10]
    [1.414214 2.1213203 2.1213203
     0.000000 0.7071068 0.7071068
     0.000000 0.0000000 2.2360680] [2 3 3
                                    3 5 5
                                    3 5 10]
    [1.414214 2.1213203 2.1213203 2.828427e+00
     0.000000 0.7071068 0.7071068 1.256074e-15
     0.000000 0.0000000 2.2360680 2.683282e+00
     0.000000 0.0000000 0.0000000 2.190890e+00] [2 3 3  4
                                                 3 5 5  6
                                                 3 5 10 12
                                                 4 6 12 20])
  (t/is (v/delta-eq (seq (sut/mat->array (sut/cholesky (sut/mat->array2d (sut/mat2x2 2 3 3 10)) true)))
                    '(1.4142135623730951 2.1213203435596424 0.0 2.345207879911715)))
  (t/is (v/delta-eq (seq (sut/mat->array (sut/cholesky (sut/mat->RealMatrix (sut/mat2x2 2 3 3 10)) true)))
                    '(1.4142135623730951 2.1213203435596424 0.0 2.345207879911715))))

;;

(t/deftest norm
  (t/are [t m r] (m/delta-eq r (sut/norm m t) 1.0e-5)
    1 m22 13
    2 m22 11.33421
    :inf m22 15
    :frobenius m22 11.74734
    :max m22 10
    [1] m22 14.422205101855956

    1 m33 13
    2 m33 11.4708
    :inf m33 17
    :frobenius m33 13
    :max m33 7
    [1] m33 18.71064301229419
    
    1 m44 13
    2 m44 9.351205
    :inf m44 13
    :frobenius m44 12.12436
    :max m44 5
    [1] m44 20.32997715038671

    1 m44a 13
    2 m44a 9.351205
    :inf m44a 13
    :frobenius m44a 12.12436
    :max m44a 5
    [1] m44a 20.32997715038671

    1 m44ra 13
    2 m44ra 9.351205
    :inf m44ra 13
    :frobenius m44ra 12.12436
    :max m44ra 5
    [1] m44ra 20.32997715038671))
