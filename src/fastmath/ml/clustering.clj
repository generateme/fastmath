(ns fastmath.ml.clustering
  (:require [fastmath.distance :as dists]
            [fastmath.random :as r]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.stats :as stats])
  (:import [org.apache.commons.math3.ml.clustering
            KMeansPlusPlusClusterer
            KMeansPlusPlusClusterer$EmptyClusterStrategy
            MultiKMeansPlusPlusClusterer
            FuzzyKMeansClusterer
            DBSCANClusterer
            DoublePoint Clusterable Cluster Clusterer CentroidCluster]
           [org.apache.commons.math3.ml.distance
            EuclideanDistance ManhattanDistance ChebyshevDistance CanberraDistance EarthMoversDistance]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- ->acm-distance [d]
  (dists/->acm-distance (case d
                          :euclidean (EuclideanDistance.)
                          :chebyshev (ChebyshevDistance.)
                          :manhattan (ManhattanDistance.)
                          :canberra (CanberraDistance.)
                          :earth-movers (EarthMoversDistance.)
                          d)))

(def ^:private kmeans-strategies
  {:largest-variance KMeansPlusPlusClusterer$EmptyClusterStrategy/LARGEST_VARIANCE
   :largest-points-number KMeansPlusPlusClusterer$EmptyClusterStrategy/LARGEST_POINTS_NUMBER
   :farthest-point KMeansPlusPlusClusterer$EmptyClusterStrategy/FARTHEST_POINT
   :error KMeansPlusPlusClusterer$EmptyClusterStrategy/ERROR})

(defn- point->vec [^Clusterable p] (vec (.getPoint p)))

(defn- cluster->data
  [add-data? ^Cluster cluster-data]
  (let [data (mapv point->vec (.getPoints cluster-data))
        res {:representative (if (instance? CentroidCluster cluster-data)
                               (point->vec (.getCenter ^CentroidCluster cluster-data))
                               (v/average-vectors data))
             :size (count data)}]
    (if add-data? (assoc res :data data) res)))

(defn- ->clusters
  [xss ^Clusterer clusterer add-data?]
  (->> xss
       (map (fn [xs] (DoublePoint. (v/vec->array xs))))
       (.cluster clusterer)
       (map (partial cluster->data add-data?))))

(defn kmeans++
  ([xss] (kmeans++ xss nil))
  ([xss {:keys [^long clusters ^long max-iters distance rng empty-cluster-strategy ^long trials add-data?]
         :or {clusters 1 max-iters -1 distance :euclidean empty-cluster-strategy :largest-variance
              trials 1 rng (r/rng :default) add-data? true}}]
   (let [strategy (kmeans-strategies empty-cluster-strategy)
         dist (->acm-distance distance)
         clusterer (KMeansPlusPlusClusterer. clusters max-iters dist rng strategy)
         clusterer (if (== trials 1)
                     clusterer
                     (MultiKMeansPlusPlusClusterer. clusterer trials))]
     (->clusters xss clusterer add-data?))))

(defn fuzzy-kmeans
  ([xss] (fuzzy-kmeans xss nil))
  ([xss {:keys [^long clusters ^double fuzziness ^long max-iters distance  ^double epsilon rng add-data?]
         :or {clusters 1 fuzziness 2 max-iters -1 distance :euclidean epsilon 1.0e-3
              rng (r/rng :default) add-data? true}}]
   (let [dist (->acm-distance distance)
         clusterer (FuzzyKMeansClusterer. clusters fuzziness max-iters dist epsilon rng)]
     (->clusters xss clusterer add-data?))))

(defn infer-dbscan-radius
  [xss dist]
  (let [indexed (map-indexed vector xss)]
    (-> (for [[^long id1 a] indexed
              [^long id2 b] indexed
              :when (m/< id1 id2)]
          (dist a b))
        (kmeans++ {:clusters 2 :add-data? false})
        (->> (map (comp first :representative))
             (reduce m/min)))))

(defn dbscan
  ([xss] (dbscan xss nil))
  ([xss {:keys [^double eps ^long points distance add-data?]
         :or {distance :euclidean add-data? true}}]
   (let [points (or points (* 2 (alength (v/vec->array (first xss)))))
         dist (->acm-distance distance)
         radius (or eps (infer-dbscan-radius xss dist))
         clusterer (DBSCANClusterer. radius points dist)]
     (->clusters xss clusterer add-data?))))
