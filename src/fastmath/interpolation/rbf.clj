;; http://shihchinw.github.io/2018/10/data-interpolation-with-radial-basis-functions-rbfs.html
(ns fastmath.interpolation.rbf
  (:require [fastmath.vector :as v]
            [fastmath.distance :as d]
            [fastmath.kernel :as k]
            [fastmath.core :as m]
            [fastmath.matrix :as mat]
            [fastmath.interpolation.common :as ic])
  (:import [org.apache.commons.math3.linear MatrixUtils
            RealMatrix RealVector]))

(defn- phi-matrix [xss kernel ^double kscale distance]
  (->> (for [xs1 xss]
         (for [xs2 xss]
           (* kscale ^double (kernel (distance xs1 xs2)))))
       (mat/rows->RealMatrix)))

(defn- solve
  ([^RealMatrix A ^RealVector b]
   (ic/solve-LU-or-SVD A b))
  ([^RealMatrix A ^RealVector b ^double lambda]
   (if (m/zero? lambda)
     (solve A b)
     (let [At (mat/transpose A)]
       (solve (mat/add (mat/mulm At A)
                       (MatrixUtils/createRealDiagonalMatrix
                        (m/seq->double-array (repeat (.getDimension b) lambda))))
              (mat/mulv At b))))))

(defn- phi-vector [xss x kernel kscale distance]
  (let [s (double kscale)]
    (mapv (fn [xs] (* s ^double (kernel (distance xs x)))) xss)))

(defn rbf
  ([xss ys] (rbf xss ys nil))
  ([xss ys kernel] (rbf xss ys kernel nil))
  ([xss ys kernel {:keys [^double kscale ^double lambda distance polynomial-terms]
                   :or {kscale 1.0 lambda 0.0}}]
   (let [kernel (or kernel (k/rbf :gaussian))
         ysize (count ys)
         dist (or distance (if (number? (first xss)) d/euclidean-1d d/euclidean))
         Phi (phi-matrix xss kernel kscale dist)
         [Phi ys] (if polynomial-terms (ic/pi-matrix-ys xss Phi ys polynomial-terms) [Phi ys])
         ysv (v/vec->RealVector ys)
         w (solve Phi ysv lambda)
         [w c] (if polynomial-terms (split-at ysize w) [w nil])]
     (if polynomial-terms
       (fn ^double [x] (+ (v/dot w (phi-vector xss x kernel kscale dist))
                         (v/dot c (polynomial-terms x))))
       (fn ^double [x] (v/dot w (phi-vector xss x kernel kscale dist)))))))


