(ns fastmath.gp
  "Gaussian Processes

  See more [here](https://nextjournal.com/generateme/gaussian-processes#gp%2B)"
  (:require [fastmath.core :as m]
            [fastmath.kernel :as k]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [fastmath.matrix :as mat])
  (:import [clojure.lang IFn]
           [org.apache.commons.math3.stat StatUtils]
           [org.apache.commons.math3.linear DiagonalMatrix CholeskyDecomposition RealMatrix]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;; native gaussian processes
;; Based on: https://www.cs.ubc.ca/~nando/540-2013/lectures/gp.py

(defn- ensure-vectors
  [xs]
  (let [res (if (sequential? (first xs)) xs (mapv vector xs))]
    (if (vector? res) res (vec res))))

(defn- kernel-cov-matrix
  ([xss kernel ^double scale] (kernel-cov-matrix xss xss kernel scale))
  ([xss xss* kernel ^double scale]
   (-> (for [vj xss]
         (for [vi xss*]
           (kernel vi vj)))
       (mat/rows->RealMatrix)
       (mat/muls scale))))

(defn- cov-matrix
  [xss kernel ^double scale ^double noise]
  (mat/add (kernel-cov-matrix xss kernel scale)
           (DiagonalMatrix. (m/seq->double-array (repeat (count xss) noise)))))


(declare predict)
(declare predict-all)

(deftype GaussianProcess [^double kscale kernel ^double noise xss ys ^double ymean ^double ystddev
                          chol w ^long n ^double L]
  IFn
  (invoke [this xs]
    (if (and (sequential? xs)
             (sequential? (first xs)))
      (predict-all this xs)
      (predict this xs)))
  (invoke [this xs stddev?]
    (if (and (sequential? xs)
             (sequential? (first xs)))
      (predict-all this xs stddev?)
      (predict this xs stddev?))))

(def ^{:const true :private true :tag 'double} LOG2PI (m/log m/TWO_PI))

(defn L
  ([ys w ^CholeskyDecomposition chol ^long n]
   (* -0.5 (+ (v/dot ys w)
              (m/log (.getDeterminant chol))
              (* n LOG2PI))))
  ([^GaussianProcess gp-object]
   (if (m/invalid-double? (.L gp-object))
     (L (.ys gp-object)
        (.w gp-object)
        (.chol gp-object)
        (.n gp-object))
     (.L gp-object))))

(defn gaussian-process
  ([xss ys] (gaussian-process xss ys nil))
  ([xss ys {:keys [^double kscale kernel ^double noise normalize? L?]
            :or {kscale 1.0 kernel (k/kernel :gaussian 1.0) normalize? false noise 1.0e-8 L? true}}]
   (let [xss (ensure-vectors xss)
         ys (m/seq->double-array ys)
         normalize? (if (= 1 (alength ys)) false normalize?)
         ymean (if normalize? (StatUtils/mean ys) 0.0)
         ystddev (if normalize? (m/sqrt (StatUtils/variance ys)) 1.0)
         ys (v/vec->RealVector (if normalize? (StatUtils/normalize ys) ys))
         ^CholeskyDecomposition chol (CholeskyDecomposition. (cov-matrix xss kernel kscale noise))
         w (.solve (.getSolver chol) ys)
         n (count xss)
         l (if-not L? ##NaN (L ys w chol n))]
     (->GaussianProcess kscale kernel noise xss ys ymean ystddev chol w n l))))

(defn predict
  ([gp-object xval] (predict gp-object xval false))
  ([^GaussianProcess gp-object xval stddev?]
   (let [xtest (if (sequential? xval) xval [xval])
         cov-vector (-> (map (fn [xs]
                               (* (.kscale gp-object) ^double ((.kernel gp-object) xtest xs))) (.xss gp-object))
                        (v/vec->RealVector))
         mu (->> (v/dot (.w gp-object) cov-vector)
                 (* (.ystddev gp-object))
                 (+ (.ymean gp-object)))]
     (if-not stddev?
       mu
       (v/vec2 mu (->> cov-vector
                       (.solve (.getSolver ^CholeskyDecomposition (.chol gp-object)))
                       (v/dot cov-vector)
                       (- (* (.kscale gp-object) ^double ((.kernel gp-object) xtest xtest)))
                       (m/safe-sqrt)
                       (* (.ystddev gp-object))))))))

(defn- predict-all-
  ([gp-object xvals] (predict-all- gp-object xvals false))
  ([^GaussianProcess gp-object xvals stddev?]
   (let [xvals (ensure-vectors xvals)
         cov-vectors (kernel-cov-matrix xvals (.xss gp-object) (.kernel gp-object) (.kscale gp-object))
         mus (-> cov-vectors
                 (mat/mulv (.w gp-object))
                 (v/mult (.ystddev gp-object))
                 (v/shift (.ymean gp-object))
                 (v/vec->array))]
     (if-not stddev?
       mus
       (let [^RealMatrix kv (mat/transpose cov-vectors)
             sd (->> (.solve (.getSolver ^CholeskyDecomposition (.chol gp-object)) kv)
                     (mat/mulm cov-vectors)
                     (mat/sub (kernel-cov-matrix xvals (.kernel gp-object) (.kscale gp-object))))]
         [mus sd])))))

(defn- with-sd
  [^double ystddev mus sd]
  (->> (mat/diag sd)
       (v/vec->array)
       (map (fn [^double sd] (* ystddev  (m/sqrt sd))))
       (map vector mus)))


(defn predict-all
  ([gp-object xvals] (predict-all gp-object xvals false))
  ([^GaussianProcess gp-object xvals stddev?]
   (let [res (predict-all- gp-object xvals stddev?)]
     (if-not stddev?
       (seq res)
       (apply with-sd (.ystddev gp-object) res)))))

(defn prior-samples
  [^GaussianProcess gp-object xvals]
  (let [xvals (ensure-vectors xvals)
        cov-vectors (cov-matrix xvals (.kernel gp-object) (.kscale gp-object) (.noise gp-object))
        L (.getL ^CholeskyDecomposition (CholeskyDecomposition. cov-vectors))]
    (->> (repeatedly (count xvals) r/grand)
         (v/vec->RealVector)
         (mat/mulv L)
         (v/vec->array)
         (seq))))

(defn posterior-samples
  ([gp-object xvals] (posterior-samples gp-object xvals false))
  ([^GaussianProcess gp-object xvals stddev?]
   (let [[mu sd] (predict-all- gp-object xvals true)
         L (.getL ^CholeskyDecomposition (CholeskyDecomposition. sd 1.0e-6 1.0e-10))
         post (->> (repeatedly (count xvals) r/grand)
                   (v/vec->RealVector)
                   (mat/mulv L)
                   (v/vec->array)
                   (v/add mu)
                   (seq))         ]
     (if-not stddev?
       post
       (with-sd (.ystddev gp-object) post sd)))))

(m/unuse-primitive-operators)

