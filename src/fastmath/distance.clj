(ns fastmath.distance
  "Distance objects.

  Objects implement IFn, Smile and Apache Commons Math distance interfaces."
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [clojure.lang IFn]
           [org.apache.commons.math3.ml.distance DistanceMeasure CanberraDistance EarthMoversDistance]
           [smile.math.distance Distance Metric EuclideanDistance ManhattanDistance ChebyshevDistance
            CorrelationDistance JensenShannonDistance MahalanobisDistance]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn- gen-smile
  [^Distance dst]
  (reify
    IFn (invoke [d x y] (.d dst (m/seq->double-array x) (m/seq->double-array y)))
    Metric Distance (d [_ x y] (.d dst (m/seq->double-array x) (m/seq->double-array y)))
    DistanceMeasure (compute [_ x y] (.d dst (m/seq->double-array x) (m/seq->double-array y)))))

(defn- gen-apache
  [^DistanceMeasure dst]
  (reify
    IFn (invoke [d x y] (.compute dst (m/seq->double-array x) (m/seq->double-array y)))
    Metric Distance (d [_ x y] (.compute dst (m/seq->double-array x) (m/seq->double-array y)))
    DistanceMeasure (compute [_ x y] (.compute dst (m/seq->double-array x) (m/seq->double-array y)))))

(defn- gen-fn
  [f]
  (reify
    IFn (invoke [_ x y] (f x y))
    Metric Distance (d [_ x y] (f x y))
    DistanceMeasure (compute [_ x y] (f x y))))

;; static distances

(defonce ^{:doc "Euclidean distance"} euclidean (gen-smile (EuclideanDistance.)))
(defonce ^{:doc "Manhattan distance"} manhattan (gen-smile (ManhattanDistance.)))
(defonce ^{:doc "Chebyshev distance"} chebyshev (gen-smile (ChebyshevDistance.)))
(defonce ^{:doc "Correlation distance"} correlation (gen-smile (CorrelationDistance.)))
(defonce ^{:doc "Canberra distance"} canberra (gen-apache (CanberraDistance.)))
(defonce ^{:doc "Earth-Movers distance"} earth-movers (gen-apache (EarthMoversDistance.)))
(defonce ^{:doc "Squared Euclidean distance"} euclidean-sq (gen-fn v/dist-sq))
(defonce ^{:doc "Discrete distance"} discrete (gen-fn v/dist-discrete))
(defonce ^{:doc "Cosine similarity"} cosine (gen-fn v/sim-cos))
(defonce ^{:doc "Angular distance"} angular (gen-fn v/dist-ang))
(defonce ^{:doc "Jensen-Shannon distance"} jensen-shannon (gen-smile (JensenShannonDistance.)))

(defn make-mahalanobis
  "Create Mahalonobis distance for given covariance (seq of seqs) matrix."
  [cov] (gen-smile (MahalanobisDistance. (m/seq->double-double-array cov))))

(defn make-minkowski
  "Create Minkowski distance for order `p` and optional `weights` for each dimension."
  ([^double p weights]
   (let [rp (/ p)]
     (gen-fn (fn [x y]
               (m/pow (reduce m/fast+ 0.0 (map (fn [^double w ^double xx ^double yy]
                                                 (* w (m/pow (m/abs (- xx yy)) p))) weights x y)) rp)))))
  ([^double p]
   (let [rp (/ p)]
     (gen-fn (fn [x y]
               (m/pow (reduce m/fast+ 0.0 (map (fn [^double xx ^double yy]
                                                 (m/pow (m/abs (- xx yy)) p)) x y)) rp))))))

