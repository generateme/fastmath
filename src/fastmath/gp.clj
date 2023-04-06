(ns fastmath.gp
  "Gaussian Processes

  See more [here](https://nextjournal.com/generateme/gaussian-processes#gp%2B)"
  (:require [fastmath.core :as m]
            [fastmath.kernel :as k]
            [fastmath.vector :as v]
            [fastmath.random :as r])
  (:import [clojure.lang IFn]
           [org.apache.commons.math3.stat StatUtils]
           [smile.math.matrix Matrix Matrix$Cholesky]
           [smile.math.blas UPLO]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;; native gaussian processes
;; Based on: https://www.cs.ubc.ca/~nando/540-2013/lectures/gp.py

(defn- ensure-vectors
  [xs]
  (let [res (if (sequential? (first xs)) xs (mapv vector xs))]
    (if (vector? res) res (vec res))))

(defn- kernel-cov-matrix
  (^Matrix [kernel scale xss] (kernel-cov-matrix kernel scale xss xss))
  (^Matrix [kernel ^double scale xss xss*]
   (let [cnt (count xss)
         cnt* (count xss*)
         ^Matrix m (Matrix. cnt cnt*)]
     (dotimes [j cnt]
       (let [vj (xss j)]
         (dotimes [i cnt*]
           (let [vi (xss* i)
                 k (* scale ^double (kernel vi vj))]
             (.set m j i k)))))
     m)))

(declare predict)
(declare predict-all)

(deftype GaussianProcess [^double kscale kernel ^double noise xss ys ^double ymean ^double ystddev
                          ^Matrix$Cholesky chol w ^long n ^double L]
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
  ([ys w ^Matrix$Cholesky chol ^long n]
   (* -0.5 (+ ^double (reduce m/fast+ (map m/fast* ys w))
              (.logdet chol)
              (* n LOG2PI))))
  ([^GaussianProcess gp-object]
   (if (m/invalid-double? (.L gp-object))
     (L (.ys gp-object)
        (.w gp-object)
        (.chol gp-object)
        (.n gp-object))
     (.L gp-object))))

(defn- cov-matrix
  ^Matrix [xss ^double noise kernel kscale]
  (let [cov (kernel-cov-matrix kernel kscale xss)
        ^Matrix res (if (pos? noise)
                      (reduce (fn [^Matrix m ^long id]
                                (.add m id id noise)
                                m) cov (range (count xss)))
                      cov)]
    (.uplo res UPLO/LOWER)))

(defn gaussian-process
  ([xss ys] (gaussian-process xss ys nil))
  ([xss ys {:keys [^double kscale kernel ^double noise normalize? L?]
            :or {kscale 1.0 kernel (k/kernel :gaussian 1.0) normalize? false noise 1.0e-8 L? true}}]
   (let [xss (ensure-vectors xss)
         ys (m/seq->double-array ys)
         normalize? (if (= 1 (alength ys)) false normalize?)
         ymean (if normalize? (StatUtils/mean ys) 0.0)
         ystddev (if normalize? (m/sqrt (StatUtils/variance ys)) 1.0)
         ys (if normalize? (StatUtils/normalize ys) ys)
         ^Matrix$Cholesky chol (.cholesky ^Matrix (cov-matrix xss noise kernel kscale) true)
         w (.solve chol ys)
         n (count xss)
         l (if-not L? ##NaN (L ys w chol n))]
     (->GaussianProcess kscale kernel noise xss ys ymean ystddev chol w n l))))

(defn predict
  ([gp-object xval] (predict gp-object xval false))
  ([^GaussianProcess gp-object xval stddev?]
   (let [xtest (if (sequential? xval) xval [xval])
         cov-vector (double-array (map (fn [xs]
                                         (* (.kscale gp-object) ^double ((.kernel gp-object) xtest xs))) (.xss gp-object)))
         mu (->> cov-vector
                 (map (fn [^double w ^double cv]
                        (* w cv)) (.w gp-object))
                 (reduce m/fast+)
                 (double)
                 (* (.ystddev gp-object))
                 (+ (.ymean gp-object)))]
     (if-not stddev?
       mu
       [mu (->> cov-vector
                (.solve ^Matrix$Cholesky (.chol gp-object))
                (v/dot cov-vector)
                (- (* (.kscale gp-object) ^double ((.kernel gp-object) xtest xtest)))
                (m/safe-sqrt)
                (* (.ystddev gp-object)))]))))

(defn- predict-all-
  ([gp-object xvals] (predict-all- gp-object xvals false))
  ([^GaussianProcess gp-object xvals stddev?]
   (let [xvals (ensure-vectors xvals)
         ^Matrix cov-vectors (kernel-cov-matrix (.kernel gp-object) (.kscale gp-object) xvals (.xss gp-object))
         mus (v/add (v/mult (.mv cov-vectors (.w gp-object)) (.ystddev gp-object))
                    (double-array (repeat (count xvals) (.ymean gp-object))))]
     (if-not stddev?
       mus
       (let [kv (.clone (.transpose cov-vectors))
             sd (do (.solve ^Matrix$Cholesky (.chol gp-object) kv)
                    (.sub ^Matrix (cov-matrix xvals 0.0 (.kernel gp-object) (.kscale gp-object))
                          (.mm cov-vectors kv)))]
         [mus sd])))))

(defn- with-sd
  [^double ystddev mus ^Matrix sd]
  (map vector mus (map (fn [^double sd]
                         (* ystddev  (m/sqrt sd))) (.diag sd))))

(defn predict-all
  ([gp-object xvals] (predict-all gp-object xvals false))
  ([^GaussianProcess gp-object xvals stddev?]
   (let [res (predict-all- gp-object xvals stddev?)]
     (if-not stddev?
       (seq res)
       (apply with-sd (.ystddev gp-object) res)))))

(defn- make-L
  ^Matrix [^Matrix m]
  (let [dim (.nrows m)]
    (doseq [^long r (range (dec dim))]
      (doseq [^long c (range (inc r) dim)]
        (.set m r c 0.0))))
  (.uplo m nil))

(defn- draw-multi
  ([mu cov cnt]
   (v/add mu (draw-multi cov cnt)))
  ([cov cnt]
   (-> cov
       (make-L)
       (.mv (double-array (repeatedly cnt r/grand))))))

(defn prior-samples
  [^GaussianProcess gp-object xvals]
  (let [xvals (ensure-vectors xvals)
        cov (.lu (.cholesky ^Matrix (cov-matrix xvals 0.0 (.kernel gp-object) (.kscale gp-object)) true))]
    (seq (draw-multi cov (count xvals)))))

(defn posterior-samples
  ([gp-object xvals] (posterior-samples gp-object xvals false))
  ([^GaussianProcess gp-object xvals stddev?]
   (let [[mu ^Matrix sd] (predict-all- gp-object xvals true)
         cov (.lu (.cholesky sd false))
         post (draw-multi mu cov (count xvals))]
     (if-not stddev?
       (seq post)
       (with-sd (.ystddev gp-object) post sd)))))
