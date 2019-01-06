(ns fastmath.distance
  "Distance functions.

  Distances objects which work with Smile, Apache Commons, fastmath vectors and native clojure sequences."
  (:require [fastmath.core :as m]
            [fastmath.vector :as v])
  (:import [clojure.lang IFn]
           [org.apache.commons.math3.ml.distance DistanceMeasure CanberraDistance EarthMoversDistance]
           [smile.math.distance Distance EuclideanDistance ManhattanDistance ChebyshevDistance
            CorrelationDistance JensenShannonDistance MahalanobisDistance MinkowskiDistance]))

(set! *warn-on-reflection* true)

(defprotocol DistProto
  (dist [t v1 v2]))

(defn- gen-smile
  [^Distance dst]
  (reify
    IFn (invoke [d x y] (dist d x y))
    Distance (d [_ x y] (.d dst ^doubles x ^doubles y))
    DistanceMeasure (compute [_ x y] (.d dst x y))
    DistProto (dist [_ v1 v2] (.d dst (double-array v1) (double-array v2)))))

(defn- gen-apache
  [^DistanceMeasure dst]
  (reify
    IFn (invoke [d x y] (dist d x y))
    Distance (d [_ x y] (.compute dst x y))
    DistanceMeasure (compute [_ x y] (.compute dst x y))
    DistProto (dist [_ v1 v2] (.compute dst (double-array v1) (double-array v2)))))

(defn- gen-fn
  [f]
  (reify
    IFn (invoke [_ x y] (f x y))
    Distance (d [_ x y] (f x y))
    DistanceMeasure (compute [_ x y] (f x y))
    DistProto (dist [_ v1 v2] (f v1 v2))))

(defmulti distance "Create distance object."
  (fn [k & r] k))

(defmethod distance :euclidean [_ & _] (gen-smile (EuclideanDistance.)))
(defmethod distance :manhattan [_ & _] (gen-smile (ManhattanDistance.)))
(defmethod distance :chebyshev [_ & _] (gen-smile (ChebyshevDistance.)))
(defmethod distance :correlation [_ & _] (gen-smile (CorrelationDistance.)))
(defmethod distance :jensen-shannon [_ & _] (gen-smile (JensenShannonDistance.)))
(defmethod distance :mahalanobis [_ & [cov]] (gen-smile (MahalanobisDistance. (m/seq->double-double-array cov))))
(defmethod distance :minkowski [_ & [p weights]]
  (if weights
    (gen-smile (MinkowskiDistance. (int p) (m/seq->double-array weights)))
    (gen-smile (MinkowskiDistance. (int p)))))

(defmethod distance :canberra [_ & _] (gen-apache (CanberraDistance.)))
(defmethod distance :earth-movers [_ & _] (gen-apache (EarthMoversDistance.)))

(defmethod distance :euclidean-sq [_ & _] (gen-fn v/dist-sq))
(defmethod distance :discrete [_ & _] (gen-fn v/dist-discrete))
(defmethod distance :cosine [_ & _] (gen-fn v/dist-cos))

;; static distances

(def euclidean (distance :euclidean))
(def manhattan (distance :manhattan))
(def chebyshev (distance :chebyshev))
(def correlation (distance :correlation))
(def canberra (distance :canberra))
(def earth-movers (distance :earth-movers))
(def euclidean-sq (distance :euclidean-sq))
(def discrete (distance :discrete))
(def cosine (distance :cosine))

;;

(def ^{:doc "List of distance creator keys."} distance-names (sort (keys (methods distance))))
