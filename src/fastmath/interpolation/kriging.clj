;; Numerical Recipes
(ns fastmath.interpolation.kriging
  (:require [fastmath.vector :as v]
            [fastmath.distance :as d]
            [fastmath.core :as m]
            [fastmath.matrix :as mat]
            [fastmath.kernel.variogram :as vgram]
            [fastmath.interpolation.common :as ic])
  (:import [org.apache.commons.math3.linear MatrixUtils]))

(defn- V-matrix [xss variogram distance error]
  (let [m (-> (for [xs1 xss]
                (for [xs2 xss]
                  (variogram (distance xs1 xs2))))
              (mat/rows->RealMatrix))]
    (if (or (and (number? error) (m/zero? error)) (not (seq error)))
      m
      (mat/sub m (MatrixUtils/createRealDiagonalMatrix
                  (m/seq->double-array (if (sequential? error)
                                         (v/sq error)
                                         (repeat (count xss) (m/sq error)))))))))

(defn- v-vector [xss x variogram distance]
  (mapv (fn [xs] (variogram (distance xs x))) xss))

(defn kriging
  ([xss ys] (kriging xss ys nil))
  ([xss ys variogram] (kriging xss ys variogram nil))
  ([xss ys variogram {:keys [error distance polynomial-terms]
                      :or {error 0.0 polynomial-terms (constantly [1.0])}}]
   (let [variogram (or variogram (vgram/fit (vgram/empirical xss ys) :gaussian))
         ysize (count ys)
         dist (or distance (if (number? (first xss)) d/euclidean-1d d/euclidean))
         V (V-matrix xss variogram dist error)
         [V ys] (if polynomial-terms (ic/pi-matrix-ys xss V ys polynomial-terms) [V ys])
         ysv (v/vec->RealVector ys)
         w (ic/solve-LU-or-SVD V ysv)
         [w c] (if polynomial-terms (split-at ysize w) [w nil])]
     (if polynomial-terms
       (fn ^double [x] (+ (v/dot w (v-vector xss x variogram dist))
                         (v/dot c (polynomial-terms x))))
       (fn ^double [x] (v/dot w (v-vector xss x variogram dist)))))))
