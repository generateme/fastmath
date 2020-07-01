(ns fastmath.gp
  "Gaussian Processes"
  (:require [fastmath.core :as m]
            [fastmath.kernel :as k]
            [fastmath.stats :as stats]
            [fastmath.vector :as v]
            [fastmath.random :as r])
  (:import [clojure.lang IFn]
           [org.apache.commons.math3.linear MatrixUtils CholeskyDecomposition Array2DRowRealMatrix ArrayRealVector RealMatrix]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;; bayesian optimization

;; native gaussian processes
;; Based on: https://www.cs.ubc.ca/~nando/540-2013/lectures/gp.py

(declare predict)

(deftype GaussianProcess [^double kscale kernel noise-fn xss ^ArrayRealVector ys ^double ymean L
                          ^ArrayRealVector alpha]
  IFn
  (invoke [gp v] (predict gp v false))
  (invoke [gp v stddev?] (predict gp v stddev?)))

(defn- ensure-vectors
  [xs]
  (if (sequential? (first xs)) xs (mapv vector xs)))

(defn- kernel-cov-matrix
  (^Array2DRowRealMatrix [kernel ^double scale xss xss*]
   (MatrixUtils/createRealMatrix
    (m/seq->double-double-array (for [x xss]
                                  (map #(* scale ^double (kernel x %)) xss*)))))
  (^Array2DRowRealMatrix [kernel xss xss*] (kernel-cov-matrix kernel 1.0 xss xss*)))

(defn gaussian-process
  ([xss ys] (gaussian-process xss ys nil))
  ([xss ys {:keys [^double kscale kernel noise normalize?]
            :or {kscale 1.0 kernel (k/kernel :gaussian 1.0) normalize? false}}]

   (let [xss (ensure-vectors xss)
         ymean (if normalize? (stats/mean ys) 0.0)
         ^ArrayRealVector ys (MatrixUtils/createRealVector (m/seq->double-array (map #(- ^double % ymean) ys)))
         ^RealMatrix data-cov (kernel-cov-matrix kernel kscale xss xss)
         noise-fn #(if (sequential? noise)
                     (take (count %) (cycle noise))
                     (repeat (count %) (or noise 1.0e-6)))
         ^RealMatrix diag (MatrixUtils/createRealDiagonalMatrix (m/seq->double-array (noise-fn xss)))
         ^CholeskyDecomposition chol (CholeskyDecomposition. (.add data-cov diag))
         L (.getL chol)
         alpha (.solve (.getSolver chol) ys)]

     (MatrixUtils/solveLowerTriangularSystem L ys)

     (->GaussianProcess kscale kernel noise-fn xss ys ymean L alpha))))

(defn predict
  ([gp-object xval] (predict gp-object xval false))
  ([^GaussianProcess gp-object xval stddev?]
   (let [xtest (if (sequential? xval) xval [xval])
         cov-vector (double-array (map #(* (.kscale gp-object) ^double ((.kernel gp-object) xtest %)) (.xss gp-object)))
         mu (+ (.ymean gp-object) ^double (v/dot cov-vector (.getDataRef ^ArrayRealVector (.alpha gp-object))))]
     (if-not stddev?
       mu
       (let [cov-v (MatrixUtils/createRealVector cov-vector)]
         (MatrixUtils/solveLowerTriangularSystem (.L gp-object) cov-v)
         [mu (m/safe-sqrt (- 1.0 (.dotProduct cov-v cov-v)))])))))

(defn predict-all
  ([gp-object xvals] (predict-all gp-object xvals false))
  ([gp-object xvals stddev?]
   (map #(predict gp-object % stddev?) xvals)))

(defn prior-samples
  [^GaussianProcess gp-object xvals]
  (let [xvals (ensure-vectors xvals)
        ^RealMatrix cov (kernel-cov-matrix (.kernel gp-object) (.kscale gp-object) xvals xvals)
        ^RealMatrix diag (MatrixUtils/createRealDiagonalMatrix (m/seq->double-array ((.noise-fn gp-object))))
        ^RealMatrix Lp (.getL ^CholeskyDecomposition (CholeskyDecomposition. (.add cov diag)))]
    (seq (.operate Lp (m/seq->double-array (repeatedly (count xvals) r/grand))))))

(defn posterior-samples
  ([gp-object xvals] (posterior-samples gp-object xvals false))
  ([^GaussianProcess gp-object xvals stddev?]
   (let [xvals (ensure-vectors xvals)
         ^RealMatrix cov (kernel-cov-matrix (.kernel gp-object) (.kscale gp-object) (.xss gp-object) xvals)
         cov-v (mapv #(.getColumnVector cov %) (range (.getColumnDimension cov)))]
     (run! #(MatrixUtils/solveLowerTriangularSystem (.L gp-object) %) cov-v)
     (let [Lk (MatrixUtils/createRealMatrix (m/seq->double-double-array (map #(.getDataRef ^ArrayRealVector %) cov-v)))
           diag (MatrixUtils/createRealDiagonalMatrix (m/seq->double-array ((.noise-fn gp-object) xvals)))
           ^RealMatrix k2 (kernel-cov-matrix (.kernel gp-object) (.kscale gp-object) xvals xvals)
           ^RealMatrix Lp (.getL ^CholeskyDecomposition (CholeskyDecomposition. (.add k2 (.subtract ^RealMatrix diag (.multiply ^RealMatrix Lk (.transpose Lk)))))) ;; opposite than in source
           mu (map (fn [v ^double n]
                     (+ (.ymean gp-object) (.dotProduct ^ArrayRealVector (.ys gp-object) v) n)) cov-v (seq (.operate Lp (m/seq->double-array (repeatedly (count xvals) r/grand)))))]
       (if-not stddev?
         mu
         (let [stddev (map (fn [^long id ^double v] (m/sqrt (- ^double (.getEntry k2 id id) v)))
                           (range (count xvals))
                           (map #(.dotProduct ^ArrayRealVector % ^ArrayRealVector %) cov-v))]
           (map vector mu stddev)))))))

