;; http://shihchinw.github.io/2018/10/data-interpolation-with-radial-basis-functions-rbfs.html
(ns fastmath.rbf
  (:require [fastmath.vector :as v]
            [fastmath.distance :as d]
            [fastmath.kernel :as k]
            [fastmath.core :as m]
            [fastmath.matrix :as mat])
  (:import [org.apache.commons.math3.linear MatrixUtils LUDecomposition RealMatrix RealVector]))

(defn- distance-1d ^double [^double x ^double y] (m/abs (m/- x y)))

(defn- phi-matrix [xss kernel ^double kscale distance]
  (->> (for [xs1 xss]
         (for [xs2 xss]
           (* kscale (kernel (distance xs1 xs2)))))
       (mat/rows->RealMatrix)))

(defn- solve
  ([^RealMatrix A ^RealVector b]
   (-> A LUDecomposition. .getSolver (.solve b) v/vec->array seq vec))
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
    (mapv (fn [xs] (* s (kernel (distance xs x)))) xss)))

(defn- pi-matrix-ys
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

(defn- split-w  [^doubles w ^long size] (split-at size w))

(defn- dot ^double [a b] (reduce m/fast+ 0.0 (map m/fast* a b)))

(defn rbf
  ([xss ys] (rbf xss ys nil))
  ([xss ys {:keys [kernel ^double kscale ^double lambda distance polynomial-terms]
            :or {kernel (k/rbf :gaussian) kscale 1.0 lambda 0.0 distance d/euclidean}}]
   (let [ysize (count ys)
         R? (every? number? xss)
         dist (if R? distance-1d distance)
         Phi (phi-matrix xss kernel kscale dist)
         [Phi ys] (if polynomial-terms (pi-matrix-ys xss Phi ys polynomial-terms) [Phi ys])
         ysv (v/vec->RealVector ys)
         w (solve Phi ysv lambda)
         [w c] (if polynomial-terms (split-w w ysize) [w nil])]
     (println w)
     (if polynomial-terms
       (fn ^double [x] (+ (dot w (phi-vector xss x kernel kscale dist))
                         (dot c (double-array (polynomial-terms x)))))
       (fn ^double [x] (dot w (phi-vector xss x kernel kscale dist)))))))

(let [f (rbf [1 2 3 4] [0.8 0.9 0.8 0.1] {:kernel (k/rbf :thin-plate)
                                          :polynomial-terms (fn [^double a] [1.0 a (* a a a)])})
      xs (m/slice-range 0.5 5 1000)
      ys (time (mapv f xs))
      ch (umontreal.ssj.charts.XYLineChart. "A" "x" "y" (m/seq->double-double-array [xs ys]) 0 1)]
  (.setManualRange ch (double-array [0 6 0 2]))
  (.view ch 600 600))
