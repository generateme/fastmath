(ns fastmath.interpolation.common
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.matrix :as mat])
  (:import [org.apache.commons.math3.linear LUDecomposition SingularValueDecomposition
            RealMatrix RealVector SingularMatrixException MatrixUtils]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;;                                v    n=4        
;; ------0---------1--------------2------------3----------
;;  -1   0   -2    1       -3     2   -4       3     -5
;;   0   0    0    1        1     2    2       2      2 

(defn- binary-search-id
  "Find array index after binary search snapping to the lower value. The last segment snaps to the previous segment."
  ^long [^long bsres ^long n]
  (if (m/neg? bsres)
    (cond
      (m/== -1 bsres) 0
      (m/> bsres (m/- n)) (m/- (m/- bsres) 2)
      :else (m/- n 2))
    (if (m/>= bsres (m/- n 2)) (m/- n 2) bsres)))

(defn binary-search
  ^long [^doubles arr ^double v]
  (-> (java.util.Arrays/binarySearch arr v)
      (binary-search-id (alength arr))))

;;

(defn solve-LU-or-SVD
  [^RealMatrix A ^RealVector b]
  (try
    (-> A LUDecomposition. .getSolver (.solve b) v/vec->array seq vec)
    (catch SingularMatrixException _
      (do
        (println "Warning. Matrix is singular, trying SVD")
        (-> A SingularValueDecomposition. .getSolver (.solve b) v/vec->array seq vec)))))

;; rbf / kriging

(defn pi-matrix-ys
  "Adds polynomial terms"
  [xss Phi ys polynomial-terms]
  (let [pterms (map polynomial-terms xss)
        P (mat/rows->RealMatrix pterms)
        pos (count xss)
        ptsize (count (first pterms))
        size (+ ptsize pos)
        ^RealMatrix target (MatrixUtils/createRealMatrix size size)]
    (.setSubMatrix target (mat/mat->array2d Phi) 0 0)
    (.setSubMatrix target (mat/mat->array2d P) 0 pos)
    (.setSubMatrix target (mat/mat->array2d (mat/transpose P)) pos 0)
    [target (concat ys (repeat ptsize 0.0))]))
