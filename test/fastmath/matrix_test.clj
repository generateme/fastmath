(ns fastmath.matrix-test
  (:require [fastmath.matrix :as sut]
            [clojure.test :as t]
            [fastmath.vector :as v]
            [fastmath.core :as m]
            [fastmath.complex :as cplx]
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
(def m33s (sut/mat3x3 1 2 3 2 -3 -2 3 -2 5))

(defn- creator
  [m]
  (condp instance? m
    fastmath.matrix.Mat2x2 sut/mat2x2
    fastmath.matrix.Mat3x3 sut/mat3x3
    fastmath.matrix.Mat4x4 sut/mat4x4))

;; `rows->matNxN`/`cols->matNxN` take N separate row/column vectors as
;; N distinct arguments (not one collection of rows/columns).
(t/deftest rows-cols-diag-constructors
  (t/is (= m22 (apply sut/rows->mat2x2 (partition 2 d22))))
  (t/is (= m22 (apply sut/cols->mat2x2 (sut/cols m22))))
  (t/is (= m33 (apply sut/rows->mat3x3 (partition 3 d33))))
  (t/is (= m33 (apply sut/cols->mat3x3 (sut/cols m33))))
  (t/is (= m44 (apply sut/rows->mat4x4 (partition 4 d44))))
  (t/is (= m44 (apply sut/cols->mat4x4 (sut/cols m44))))
  (t/are [c r] (= c r)
    (sut/diag->mat2x2 5.0) (sut/mat2x2 5.0 0.0 0.0 5.0)
    (sut/diag->mat2x2 3.0 4.0) (sut/mat2x2 3.0 0.0 0.0 4.0)
    (sut/diag->mat3x3 5.0) (sut/mat3x3 5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0)
    (sut/diag->mat3x3 3.0 4.0 5.0) (sut/mat3x3 3.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 5.0)
    (sut/diag->mat4x4 5.0) (sut/mat4x4 5.0 0.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 0.0 5.0)
    (sut/diag->mat4x4 3.0 4.0 5.0 6.0) (sut/diagonal 3.0 4.0 5.0 6.0)))

;; `real-matrix` and `rows->RealMatrix` are the same implementation under two
;; names; `cols->RealMatrix` transposes the `rows->RealMatrix` result.
(t/deftest real-matrix-constructors
  (let [rows [[1.0 2.0 3.0] [4.0 5.0 6.0]]
        cols [[1.0 4.0] [2.0 5.0] [3.0 6.0]]
        flat [1.0 2.0 3.0 4.0 5.0 6.0]]
    (t/is (instance? org.apache.commons.math3.linear.RealMatrix (sut/real-matrix rows)))
    (t/is (= flat (seq (sut/mat->array (sut/real-matrix rows)))
             (seq (sut/mat->array (sut/rows->RealMatrix rows)))))
    (t/is (= flat (seq (sut/mat->array (sut/cols->RealMatrix cols)))))
    ;; already-typed double[][] input path
    (t/is (= flat (seq (sut/mat->array (sut/real-matrix (m/seq->double-double-array rows))))))))

;; `mat`/`rows->mat`/`cols->mat`: 1-arg collection form always builds a
;; `RealMatrix`; 4/9/16-arg forms build the corresponding fixed type.
(t/deftest generic-mat-constructors
  (t/is (= org.apache.commons.math3.linear.Array2DRowRealMatrix (class (sut/mat [[1 2 3] [4 5 6]]))))
  (t/is (= [1.0 2.0 3.0 4.0 5.0 6.0] (seq (sut/mat->array (sut/mat [[1 2 3] [4 5 6]])))))
  (t/is (= m22 (apply sut/mat d22)))
  (t/is (= m33 (apply sut/mat d33)))
  (t/is (= m44 (apply sut/mat d44)))
  (t/is (= org.apache.commons.math3.linear.Array2DRowRealMatrix (class (sut/rows->mat [[1 2 3] [4 5 6]]))))
  (t/is (= m22 (apply sut/rows->mat (partition 2 d22))))
  (t/is (= m33 (apply sut/rows->mat (partition 3 d33))))
  (t/is (= m44 (apply sut/rows->mat (partition 4 d44))))
  (t/is (= org.apache.commons.math3.linear.Array2DRowRealMatrix (class (sut/cols->mat [[1 2 3] [4 5 6]]))))
  (t/is (= m22 (apply sut/cols->mat (sut/cols m22))))
  (t/is (= m33 (apply sut/cols->mat (sut/cols m33))))
  (t/is (= m44 (apply sut/cols->mat (sut/cols m44)))))

(t/deftest array2d-constructors
  (t/is (= m22 (sut/array2d->mat2x2 (m/seq->double-double-array (partition 2 d22)))))
  (t/is (= m33 (sut/array2d->mat3x3 (m/seq->double-double-array (partition 3 d33)))))
  (t/is (= m44 (sut/array2d->mat4x4 (m/seq->double-double-array (partition 4 d44)))))
  (t/is (instance? org.apache.commons.math3.linear.RealMatrix
                   (sut/array2d->RealMatrix (m/seq->double-double-array (partition 2 d22)))))
  (t/is (= d22 (seq (sut/mat->array (sut/array2d->RealMatrix (m/seq->double-double-array (partition 2 d22))))))))

;; `eye`/`zero`: size 2/3/4 (1-arg, or 2-arg with `real-matrix?` false) builds
;; the fixed type, any other size (or `real-matrix?` true) builds a `RealMatrix`.
;; `diagonal`: direct scalar args build the fixed type; a vector arg always
;; builds a `RealMatrix`, regardless of its length (documented breaking change,
;; see CHANGELOG).
(t/deftest eye-zero-diagonal-dispatch
  (t/are [s c] (= c (class (sut/eye s)))
    2 fastmath.matrix.Mat2x2
    3 fastmath.matrix.Mat3x3
    4 fastmath.matrix.Mat4x4
    1 org.apache.commons.math3.linear.Array2DRowRealMatrix
    5 org.apache.commons.math3.linear.Array2DRowRealMatrix)
  (t/is (= fastmath.matrix.Mat2x2 (class (sut/eye 2 false))))
  (t/is (= org.apache.commons.math3.linear.Array2DRowRealMatrix (class (sut/eye 2 true))))
  (t/is (= (seq (sut/mat->array (sut/eye 5))) (seq (sut/mat->array (sut/eye 5 true)))))

  (t/are [s c] (= c (class (sut/zero s)))
    2 fastmath.matrix.Mat2x2
    3 fastmath.matrix.Mat3x3
    4 fastmath.matrix.Mat4x4
    1 org.apache.commons.math3.linear.Array2DRowRealMatrix
    5 org.apache.commons.math3.linear.Array2DRowRealMatrix)
  (t/is (= fastmath.matrix.Mat2x2 (class (sut/zero 2 false))))
  (t/is (= org.apache.commons.math3.linear.Array2DRowRealMatrix (class (sut/zero 2 true))))
  (t/is (= [3 5] (sut/shape (sut/zero 3 5 false))))
  (t/is (every? zero? (seq (sut/mat->array (sut/zero 3 5 false)))))

  (t/are [c r] (= c r)
    (sut/diagonal 3.0 4.0) (sut/mat2x2 3.0 0.0 0.0 4.0)
    (sut/diagonal 3.0 4.0 5.0) (sut/mat3x3 3.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 5.0)
    (sut/diagonal 3.0 4.0 5.0 6.0) (sut/mat4x4 3.0 0.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 0.0 6.0))
  (t/is (= org.apache.commons.math3.linear.Array2DRowRealMatrix (class (sut/diagonal [3.0 4.0]))))
  (t/is (= [3.0 4.0 5.0] (v/vec->Vec (sut/diag (sut/diagonal [3.0 4.0 5.0])))))
  (t/is (m/zero? (sut/entry (sut/diagonal [3.0 4.0 5.0]) 0 1))))

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

;; independent reference: R's `kronecker()`, Rscript, 2026-09-23
(t/deftest kronecker
  (t/is (= [0.0 5.0 0.0 10.0 6.0 7.0 12.0 14.0 0.0 15.0 0.0 20.0 18.0 21.0 24.0 28.0]
           (seq (sut/mat->array (sut/kronecker (sut/mat2x2 1.0 2.0 3.0 4.0)
                                               (sut/mat2x2 0.0 5.0 6.0 7.0))))))
  (t/is (= [6 6] (sut/shape (sut/kronecker m22 m33)))))

;; independent reference: R's `Matrix::bdiag()`, `rbind()`, `cbind()`, Rscript, 2026-09-23
(t/deftest block-diagonal-and-bind
  (t/is (= [2.0 3.0 0.0 0.0 0.0
            5.0 -10.0 0.0 0.0 0.0
            0.0 0.0 -3.0 2.0 -6.0
            0.0 0.0 5.0 7.0 -5.0
            0.0 0.0 1.0 4.0 -2.0]
           (seq (sut/mat->array (sut/block-diagonal m22 m33)))))
  (t/is (= [5 5] (sut/shape (sut/block-diagonal m22 m33))))
  (t/is (= [6 6] (sut/shape (sut/block-diagonal m22 m22 m22))))
  (t/is (thrown? clojure.lang.ExceptionInfo (sut/block-diagonal m22 (sut/mat [[1.0 2.0 3.0]]))))

  (let [ones (sut/mat [[1.0 1.0] [1.0 1.0]])]
    (t/is (= [2.0 3.0 5.0 -10.0 1.0 1.0 1.0 1.0]
             (seq (sut/mat->array (sut/bind-rows m22 ones)))))
    (t/is (= [2.0 3.0 1.0 1.0 5.0 -10.0 1.0 1.0]
             (seq (sut/mat->array (sut/bind-cols m22 ones)))))
    (t/is (= [6 2] (sut/shape (sut/bind-rows m22 m22 m22))))
    (t/is (= [2 6] (sut/shape (sut/bind-cols m22 m22 m22)))))
  ;; zero-padding on mismatched dimensions
  (t/is (= [2.0 3.0 0.0 5.0 -10.0 0.0 1.0 2.0 3.0]
           (seq (sut/mat->array (sut/bind-rows m22 (sut/mat [[1.0 2.0 3.0]])))))))

;; `map-rows`/`map-cols`: applying a uniform scale-by-2 to every row/col is
;; equivalent to scaling the whole matrix by 2, for both `RealMatrix` and
;; fixed representations.
(t/deftest map-rows-and-cols
  (t/is (= (sut/muls m22 2.0) (sut/map-rows #(v/mult % 2.0) m22)))
  (t/is (= (sut/muls m22 2.0) (sut/map-cols #(v/mult % 2.0) m22)))
  (t/is (= fastmath.matrix.Mat2x2 (class (sut/map-rows #(v/mult % 2.0) m22))))
  (t/is (= (v/mult d44 2.0) (seq (sut/mat->array (sut/map-rows #(v/mult % 2.0) m44ra)))
           (seq (sut/mat->array (sut/map-cols #(v/mult % 2.0) m44ra)))))
  (t/is (instance? org.apache.commons.math3.linear.RealMatrix (sut/map-rows #(v/mult % 2.0) m44ra))))

;; independent reference: R's `diff()` applied column-wise, Rscript, 2026-09-23
(t/deftest differences
  (t/is (= [2 3] (sut/shape (sut/differences m33))))
  (t/is (= [8.0 5.0 1.0 -4.0 -3.0 3.0] (seq (sut/mat->array (sut/differences m33)))))
  (t/is (= [1 3] (sut/shape (sut/differences m33 2))))
  (t/is (= [-12.0 -8.0 2.0] (seq (sut/mat->array (sut/differences m33 2)))))
  (t/is (instance? org.apache.commons.math3.linear.RealMatrix (sut/differences m22))))

;; independent reference: R's `scale()`/`sweep()`, Rscript, 2026-09-23
(t/deftest normalize-demean-standardize
  (t/is (every? true? (map v/delta-eq (sut/cols (sut/demean m33))
                          (mapv v/vec3 [[-4.0 4.0 0.0] [-2.333333333333333 2.666666666666667 -0.33333333333333304]
                                        [-1.666666666666667 -0.666666666666667 2.333333333333333]]))))
  (t/is (every? true? (map v/delta-eq (sut/rows (sut/demean m33 true))
                          (sut/cols (sut/demean (sut/transpose m33))))))
  (t/is (every? true? (map v/delta-eq (sut/cols (sut/normalize m33))
                          (mapv v/vec3 [[-0.50709255283711 0.8451542547285165 0.1690308509457033]
                                        [0.2407717061715384 0.8427009716003844 0.4815434123430768]
                                        [-0.7442084075352509 -0.6201736729460423 -0.24806946917841693]]))))
  (t/is (every? #(m/delta-eq 1.0 (v/mag %)) (sut/cols (sut/normalize m33))))
  (t/is (every? true? (map v/delta-eq (sut/cols (sut/standardize m33))
                          (mapv v/vec3 [[-1.0 1.0 0.0] [-0.9271726499455306 1.0596258856520353 -0.13245323570650427]
                                        [-0.8006407690254359 -0.32025630761017443 1.12089707663561]])))))

;; independent reference: R's `scale(center=FALSE, scale=TRUE)`, row- and
;; column-wise, Rscript, 2026-09-23
(t/deftest scale-rows-and-cols
  (t/is (every? true? (map v/delta-eq (sut/cols (sut/scale-cols m33))
                          (mapv v/vec3 [[-0.7171371656006362 1.1952286093343936 0.23904572186687872]
                                        [0.34050261230349943 1.191759143062248 0.6810052246069989]
                                        [-1.0524696231684352 -0.8770580193070293 -0.3508232077228117]]))))
  (t/is (every? true? (map v/delta-eq (sut/rows (sut/scale-rows m33))
                          (mapv v/vec3 [[-0.6060915267313264 0.40406101782088427 -1.2121830534626528]
                                        [0.7106690545187015 0.9949366763261821 -0.7106690545187015]
                                        [0.3086066999241838 1.2344267996967353 -0.6172133998483676]]))))
  (t/is (= (sut/scale-rows m33 3.0) (sut/map-rows #(v/mult % 3.0) m33))))

;; Regression test: `shift-cols` had an extra sign inversion not present in
;; `shift-rows`, silently negating both the default (demean) and any
;; explicit constant/function shift.
(t/deftest shift-rows-and-cols
  (t/is (every? #(m/delta-eq 0.0 (v/average %) 1.0e-9) (sut/rows (sut/shift-rows m22))))
  (t/is (every? #(m/delta-eq 0.0 (v/average %) 1.0e-9) (sut/cols (sut/shift-cols m22))))
  (t/is (= [7.0 8.0 10.0 -5.0] (seq (sut/mat->array (sut/shift-rows m22 5.0)))))
  (t/is (= [7.0 8.0 10.0 -5.0] (seq (sut/mat->array (sut/shift-cols m22 5.0)))))
  ;; a positive constant shift must move every element up, for both rows and cols
  (t/is (every? #(>= % 0.0) (seq (sut/mat->array (sut/shift-rows (sut/mat2x2 -5.0 -5.0 -5.0 -5.0) 5.0)))))
  (t/is (every? #(>= % 0.0) (seq (sut/mat->array (sut/shift-cols (sut/mat2x2 -5.0 -5.0 -5.0 -5.0) 5.0))))))

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
  (t/are [m r res] (v/delta-eq (sut/row m r) res)
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

(t/deftest entry
  (t/are [m r c v] (= v (sut/entry m r c))
    m22 0 1 3.0 m22 1 0 5.0
    m33 2 1 4.0
    m44a 1 2 -1.0
    m44ra 3 3 -1.0)
  (t/is (= 7.0 (sut/entry 7.0 0 0)))
  (t/is (thrown? Exception (sut/entry 7.0 0 1)))
  (t/is (thrown? IndexOutOfBoundsException (sut/entry m22 5 5))))

(t/deftest fmap
  (t/is (= (sut/mat2x2 4.0 6.0 10.0 -20.0) (sut/fmap m22 (partial * 2.0))))
  (t/is (= (v/mult d44 2.0) (seq (sut/mat->array (sut/fmap m44a (partial * 2.0))))))
  (t/is (= (v/mult d44 2.0) (seq (sut/mat->array (sut/fmap m44ra (partial * 2.0))))))
  (t/is (= 14.0 (sut/fmap 7.0 (partial * 2.0)))))

;; `mat->seq`/`mat->array`/`mat->float-array`/`mat->float-array2d`/`mat->RealMatrix`:
;; flat/2d, double/float conversions, row order preserved across representations.
(t/deftest array-conversions
  (t/are [m d] (= d (seq (sut/mat->seq m)))
    m22 d22 m33 d33 m44 d44 m44a d44 m44ra d44)
  (t/are [m d] (= d (seq (sut/mat->array m)) (seq (sut/mat->float-array m)))
    m22 d22 m33 d33 m44 d44 m44a d44 m44ra d44)
  (t/is (= (m/double-double-array->seq (sut/mat->array2d m22))
           (m/double-double-array->seq (sut/mat->float-array2d m22))))
  (t/is (instance? org.apache.commons.math3.linear.RealMatrix (sut/mat->RealMatrix m22)))
  (t/is (= d44 (seq (sut/mat->array (sut/mat->RealMatrix m44)))))
  (t/is (= [7.0] (seq (sut/mat->seq 7.0)) (seq (sut/mat->array 7.0)) (seq (sut/mat->float-array 7.0))))
  (t/is (instance? org.apache.commons.math3.linear.RealMatrix (sut/mat->RealMatrix 7.0)))
  (t/is (= [7.0] (seq (sut/mat->array (sut/mat->RealMatrix 7.0))))))

(t/deftest sizes
  (t/are [m s] (= s (sut/nrow m) (sut/ncol m))
    m22 2 m33 3 m44 4 m44a 4 m44ra 4)
  (t/is (= 1 (sut/nrow 7.0) (sut/ncol 7.0))))

(t/deftest shape-and-square
  (t/are [m s] (= s (sut/shape m))
    m22 [2 2] m33 [3 3] m44 [4 4] m44a [4 4] m44ra [4 4] 7.0 [1 1])
  (t/are [m s] (= s (boolean (sut/square? m)))
    m22 true m33 true m44 true m44a true m44ra true 7.0 true
    (sut/zero 3 5 false) false))

;; A plain `Number` implements `MatrixProto` as a degenerate 1x1 matrix.
;; `entry`/`fmap`/`mat->seq`/array conversions/`nrow`/`ncol`/`shape` for this
;; representation are exercised above; the remaining structural ops follow here.
(t/deftest number-as-1x1-matrix
  (t/is (= [[7.0]] (sut/cols 7.0) (sut/rows 7.0)))
  (t/is (= [7.0] (sut/row 7.0 0) (sut/col 7.0 0)))
  (t/is (thrown? Exception (sut/row 7.0 1)))
  (t/is (thrown? Exception (sut/col 7.0 1)))
  (t/is (true? (sut/symmetric? 7.0)))
  (t/is (= [7.0] (sut/diag 7.0))))

(t/deftest symmetry
  (t/are [m s] (= s (boolean (sut/symmetric? m)))
    m22 false m33 false m44 false m44a false m44ra false m33s true
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
    2 m22 3 m33 4 m44 4 m44a 4 m44ra)
  ;; independent reference: R's `solve()`, Rscript, 2026-09-23
  (t/is (v/delta-eq (seq (sut/mat->array (sut/inverse m22)))
                    [0.2857143 0.08571429 0.1428571 -0.05714286] 1.0e-6))
  (t/is (v/delta-eq (seq (sut/mat->array (sut/inverse m33)))
                    [-0.06976744 0.2325581 -0.3720930
                     -0.05813953 -0.1395349 0.5232558
                     -0.15116279 -0.1627907 0.3604651] 1.0e-6))
  (t/is (v/delta-eq (seq (sut/mat->array (sut/inverse m44)))
                    [-0.7413793 -0.3965517 -0.10344828 0.5344828
                     0.7931034 0.5172414 -0.01724138 -0.3275862
                     -0.2241379 -0.1896552 0.18965517 0.1034483
                     -1.2068966 -0.4827586 -0.01724138 0.6724138] 1.0e-6)))

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
    m22 m33 m44 m44a m44ra)
  ;; independent reference: R's `det()`, Rscript, 2026-09-23
  (t/are [m d] (m/delta-eq (sut/det m) d)
    m22 -35.0 m33 -86.0 m44 116.0))

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
  ;; independent reference: R's `m + m`, Rscript, 2026-09-23
  (t/is (= [4.0 6.0 10.0 -20.0] (seq (sut/mat->array (sut/add m22 m22)))))
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
  ;; independent reference: R's `m * m` (Hadamard), Rscript, 2026-09-23
  (t/is (= [4.0 9.0 25.0 100.0] (seq (sut/mat->array (sut/emulm m22 m22)))))
  ;; independent reference: R's `m %*% m` and `m %*% t(m)`/`t(m) %*% m`, Rscript, 2026-09-23
  (t/is (= [19.0 -24.0 -40.0 115.0] (seq (sut/mat->array (sut/mulm m22 m22)))))
  (t/is (= [13.0 -20.0 -20.0 125.0] (seq (sut/mat->array (sut/mulmt m22 m22)))))
  (t/is (= [29.0 -44.0 -44.0 109.0] (seq (sut/mat->array (sut/tmulm m22 m22)))))
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

;; A plain `Number`, as a degenerate 1x1 matrix: every arithmetic op reduces
;; to the corresponding scalar op, and transposition is a no-op, so `mulmt`/
;; `tmulm`/`tmulmt` must equal plain `mulm`. `negate` and 1-arity `sub` are
;; the same operation under two names (both delegate to `prot/sub`); this
;; holds for the `Number` representation too.
(t/deftest arithmetic-number-representation
  (t/is (= 10.0 (sut/add 7.0 3.0) (sut/adds 7.0 3.0)))
  (t/is (= 4.0 (sut/sub 7.0 3.0)))
  (t/is (= -7.0 (sut/sub 7.0) (sut/negate 7.0)))
  (t/is (= 21.0 (sut/emulm 7.0 3.0) (sut/muls 7.0 3.0) (sut/mulv 7.0 3.0) (sut/vtmul 7.0 3.0)))
  (t/is (= 21.0 (sut/mulm 7.0 3.0) (sut/mulmt 7.0 3.0) (sut/tmulm 7.0 3.0) (sut/tmulmt 7.0 3.0)))
  (t/is (= 7.0 (sut/transpose 7.0) (sut/trace 7.0))))

;; `negate` and 1-arity `sub` are the same operation under two names, for
;; every representation (both delegate to `prot/sub`).
(t/deftest negate-is-sub-alias
  (t/are [m] (= (sut/negate m) (sut/sub m))
    m22 m33 m44)
  (t/are [m] (= (seq (sut/mat->array (sut/negate m))) (seq (sut/mat->array (sut/sub m))))
    m44a m44ra))

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

;; `condition A 2` (default) is `norm(A,2) x norm(inverse(A),2)`, i.e. the
;; ratio of `A`'s largest to smallest singular value -- cross-checked both
;; against `singular-values` directly, and independently against R's
;; `kappa(m, exact=TRUE, norm="2")`, Rscript, 2026-09-23.
(t/deftest condition-number
  (t/are [m] (m/delta-eq (sut/condition m)
                         (let [svs (sut/singular-values m)]
                           (/ (apply max svs) (apply min svs))))
    m22 m33 m44)
  (t/are [m k] (m/delta-eq (sut/condition m) k 1.0e-5)
    m22 3.670408 m33 9.16118 m44 19.11566)
  (t/is (= 1.0 (sut/condition 7.0)))
  (t/is (= ##Inf (sut/condition 0.0))))

;; `singular?`: zero determinant, for every representation.
(t/deftest singular
  (t/are [m s] (= s (boolean (sut/singular? m)))
    m22 false (sut/mat2x2 0.0) true
    m33 false (sut/mat3x3 0.0) true
    m44 false (sut/mat4x4 0.0) true
    m44a false (sut/mat->array2d (sut/mat4x4 0.0)) true
    m44ra false (sut/mat->RealMatrix (sut/mat4x4 0.0)) true
    7.0 false 0.0 true))

;; Regression test: `inverse` on a singular `Mat2x2`/`Mat3x3`/`Mat4x4`
;; documentedly returns `nil`; `solve` must propagate that `nil` instead of
;; crashing on a `nil` intermediate result. `RealMatrix`/`double[][]` instead
;; throw `SingularMatrixException`, a pre-existing, unchanged convention for
;; those representations.
(t/deftest singular-matrix-inverse-and-solve
  (t/are [m] (nil? (sut/inverse m))
    (sut/mat2x2 0.0) (sut/mat3x3 0.0) (sut/mat4x4 0.0))
  (t/are [m b] (nil? (sut/solve m b))
    (sut/mat2x2 0.0) (v/vec2 1 2)
    (sut/mat3x3 0.0) (v/vec3 1 2 3)
    (sut/mat4x4 0.0) (v/vec4 1 2 3 4))
  (t/is (thrown? org.apache.commons.math3.linear.SingularMatrixException
                 (sut/inverse (sut/mat->RealMatrix (sut/mat2x2 0.0)))))
  (t/is (thrown? org.apache.commons.math3.linear.SingularMatrixException
                 (sut/solve (sut/mat->RealMatrix (sut/mat2x2 0.0)) (v/vec->RealVector [1 2])))))

;; A plain `Number`, as a degenerate 1x1 matrix.
(t/deftest linalg-number-representation
  (t/is (= 7.0 (sut/det 7.0)))
  (t/is (m/delta-eq (/ 1.0 7.0) (sut/inverse 7.0)))
  (t/is (= ##Inf (sut/inverse 0.0)))
  (t/is (false? (sut/singular? 7.0)))
  (t/is (true? (sut/singular? 0.0)))
  (t/is (= 3.0 (sut/solve 7.0 21.0)))
  (t/is (= 7.0 (sut/norm 7.0) (sut/norm -7.0)))
  (t/are [t] (= 7.0 (sut/norm -7.0 t))
    1 2 :inf :max :frobenius [2 2] [1]))

;;

;; `sut/eigenvectors`'s non-`:raw` scalings return genuine (possibly complex)
;; eigenvectors, i.e. sequences of `fastmath.complex` numbers (`Vec2`). These
;; helpers let the `eigen` test below verify them via the defining relation
;; `A v = lambda v` using complex arithmetic, instead of pinning brittle,
;; sign/phase-ambiguous numeric literals.

(defn- cplx-mulv
  "`m` (a real matrix) times `v` (a seq of complex numbers), via complex arithmetic."
  [m v]
  (mapv (fn [row] (reduce cplx/add (map cplx/scale v (seq row)))) (sut/rows m)))

(defn- flatten-cplx
  "Flattens a seq of complex numbers into a flat real vector `[re im re im ...]`,
  so `v/delta-eq` (real-vector equality) can compare them."
  [v]
  (vec (mapcat (juxt cplx/re cplx/im) v)))

(defn- eigen-relation-ok?
  "True when `A v_i = lambda_i v_i` holds for every eigenvalue/eigenvector pair
  of `m`, with eigenvectors scaled per `scaling`."
  [m scaling]
  (every? true?
          (map (fn [lambda v]
                 (v/delta-eq (flatten-cplx (cplx-mulv m v))
                             (flatten-cplx (mapv #(cplx/mult lambda %) v))))
               (sut/eigenvalues m) (sut/eigenvectors m scaling))))

(defn- unit-norm?
  "True when complex vector `v` has unit Euclidean length (`cplx/norm` is the
  squared magnitude of one component)."
  [v]
  (m/delta-eq (v/sum (map cplx/norm v)) 1.0 1.0e-9))

(t/deftest eigen
  (t/testing "eigenvalues, against R's `eigen(m)$values` (`stats::eigen`), reordered to fastmath's own eigenvalue order"
    (t/are [m res] (every? identity (map v/delta-eq (map vec (sut/eigenvalues m)) res))
      m22  [[3.1414284285428486 0.0] [-11.14142842854285 0.0]]
      m33  [[-4.6874352745829464 0.0] [3.3437176372914754 2.6770267769740972] [3.3437176372914754 -2.6770267769740972]]
      m33s [[-4.8042356853662502 0.0] [1.1212661599174663 0.0] [6.6829695254487973 0.0]]
      m44  [[0.5078448033470897 2.1995609029900725] [0.5078448033470897 -2.1995609029900725] [3.5229553075271176 0.0] [6.4613550857787017 0.0]]))

  (t/testing ":lapack, against R's `eigen(m)$vectors` (`stats::eigen`), reordered to fastmath's own eigenvalue order and flattened per eigenvector. A real eigenvector may differ from R's by an overall sign -- R/LAPACK's `DTREVC` normalizes only a real eigenvector's magnitude, not its sign; a complex eigenvector is compared exactly, since both conventions force its largest-magnitude component to be real, which collapses the remaining ambiguity to the same real ±1 already tolerated below"
    (letfn [(matches-up-to-sign? [a b] (or (v/delta-eq a b) (v/delta-eq a (mapv - b))))]
      (t/are [m r-vectors] (every? true? (map matches-up-to-sign?
                                              (map flatten-cplx (sut/eigenvectors m :lapack))
                                              r-vectors))
        m22  [[0.93463573000423483 0.0 0.35560659751957779 0.0]
              [-0.22256004869781562 0.0 0.97491898367178487 0.0]]
        m33  [[-0.93057911761573098 0.0 0.33373876078968029 0.0 -0.15046908454594965 0.0]
              [-0.10605071047019334 0.2339615258207475 0.81061772319998115 0.0 0.48671874637681295 -0.20004754435847763]
              [-0.10605071047019334 -0.2339615258207475 0.81061772319998115 0.0 0.48671874637681295 0.20004754435847763]]
        m33s [[0.44845463039893491 0.0 -0.83892980769391001 0.0 -0.3083589178804701 0.0]
              [0.7787333900938993 0.0 0.53606412633984912 0.0 -0.32589808162116735 0.0]
              [-0.43870576885495544 0.0 0.093978881745602205 0.0 -0.89370309284416627 0.0]]
        m44  [[-0.49902050637769313 -0.25591965989483245 0.43640885181942368 0.030142639626960701 -0.095249142437944828 -0.11702873250912969 -0.68655245528362796 0.0]
              [-0.49902050637769313 0.25591965989483245 0.43640885181942368 -0.030142639626960701 -0.095249142437944828 0.11702873250912969 -0.68655245528362796 0.0]
              [0.31408833643093526 0.0 -0.53513564844473704 0.0 0.70568896868117836 0.0 0.34202548759321061 0.0]
              [0.090664901976302267 0.0 0.33110377697542054 0.0 0.76448126145445472 0.0 0.54563592743444245 0.0]])))
  
  (t/testing ":raw is exactly the (unscaled) columns of the V matrix, for both backends"
    (doseq [m [m22 m33 m33s m44]
            backend [:acm :colt]]
      (let [ed (sut/eigen-decomposition m {:backend backend :eigenvectors-scaling :raw})
            v-cols (mapv seq (sut/cols (sut/decomposition-component ed :V)))]
        (t/is (= (mapv vec (sut/decomposition-component ed :eigenvectors))
                 (mapv vec v-cols))))))

  (t/testing "false, true/:normalized and :lapack all satisfy A v = lambda v"
    (doseq [m [m22 m33 m33s m44]
            scaling [false true :normalized :lapack]]
      (t/is (eigen-relation-ok? m scaling))))

  (t/testing "true and :normalized are equivalent scaling requests"
    (doseq [m [m22 m33 m33s m44]]
      (t/is (= (sut/eigenvectors m true) (sut/eigenvectors m :normalized)))))

  (t/testing "true/:normalized and :lapack eigenvectors have unit length"
    (doseq [m [m22 m33 m33s m44]
            scaling [true :lapack]]
      (t/is (every? unit-norm? (sut/eigenvectors m scaling)))))

  (t/testing ":lapack leaves real eigenvalues' eigenvectors unchanged (identical to :normalized)"
    (doseq [m [m22 m33s]]
      (t/is (= (sut/eigenvectors m :normalized) (sut/eigenvectors m :lapack)))))

  (t/testing ":lapack rotates each complex eigenvector so its largest-magnitude component is real"
    (doseq [m [m33 m44]
            [lambda v] (map vector (sut/eigenvalues m) (sut/eigenvectors m :lapack))
            :when (not (m/zero? (cplx/im lambda)))]
      (t/is (m/near-zero? (cplx/im (v (v/maxdim (mapv cplx/norm v)))) 1.0e-9)))))
