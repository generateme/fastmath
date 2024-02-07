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
    IFn (invoke [_ x y] (.d dst (m/seq->double-array x) (m/seq->double-array y)))
    Metric Distance (d [_ x y] (.d dst (m/seq->double-array x) (m/seq->double-array y)))
    DistanceMeasure (compute [_ x y] (.d dst (m/seq->double-array x) (m/seq->double-array y)))))

(defn- gen-apache
  [^DistanceMeasure dst]
  (reify
    IFn (invoke [_ x y] (.compute dst (m/seq->double-array x) (m/seq->double-array y)))
    Metric Distance (d [_ x y] (.compute dst (m/seq->double-array x) (m/seq->double-array y)))
    DistanceMeasure (compute [_ x y] (.compute dst (m/seq->double-array x) (m/seq->double-array y)))))

(defn- gen-fn
  [f]
  (reify
    IFn (invoke [_ x y] (f x y))
    Metric Distance (d [_ x y] (f x y))
    DistanceMeasure (compute [_ x y] (f x y))))

;; static distances

(def ^{:doc "Euclidean distance"} euclidean (gen-smile (EuclideanDistance.)))
(def ^{:doc "Manhattan distance"} manhattan (gen-smile (ManhattanDistance.)))
(def ^{:doc "Chebyshev distance"} chebyshev (gen-smile (ChebyshevDistance.)))
(def ^{:doc "Correlation distance"} correlation (gen-smile (CorrelationDistance.)))
(def ^{:doc "Canberra distance"} canberra (gen-apache (CanberraDistance.)))
(def ^{:doc "Earth-Movers distance"} earth-movers (gen-apache (EarthMoversDistance.)))
(def ^{:doc "Squared Euclidean distance"} euclidean-sq (gen-fn v/dist-sq))
(def ^{:doc "Discrete distance"} discrete (gen-fn v/dist-discrete))
(def ^{:doc "Cosine similarity"} cosine (gen-fn v/sim-cos))
(def ^{:doc "Angular distance"} angular (gen-fn v/dist-ang))
(def ^{:doc "Jensen-Shannon distance"} jensen-shannon (gen-smile (JensenShannonDistance.)))

(def ^{:doc "Weierstrass metric"} weierstrass
  (gen-fn (fn [a b] (m/acosh (- (* (m/sqrt (inc (v/magsq a)))
                                  (m/sqrt (inc (v/magsq b))))
                               (v/dot a b))))))

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
;;

(defn bound
  "Create relative distance function to be bound by a `diameter` (default: by `1.0`)."
  ([distance] (bound distance 1.0))
  ([distance ^double diameter]
   (gen-fn (fn [a b] (let [d (double (distance a b))] (/ (* diameter d) (inc d)))))))

(defn shoenberg-transform
  "Create relative distance function by applying Shoenberg transform."
  ([distance] (shoenberg-transform distance 1.0))
  ([distance ^double lambda]
   (gen-fn (fn [a b] (- 1.0 (m/exp (* -1.0 lambda ^double (distance a b))))))))

