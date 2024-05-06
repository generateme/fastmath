(ns fastmath.distance
  "Various distance functions"
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.matrix :as mat])
  (:import [clojure.lang IFn]
           [org.apache.commons.math3.ml.distance DistanceMeasure]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;; static distances

(defn euclidean-1d "Euclidean distance for doubles" ^double [^double v1 ^double v2] (m/abs (- v1 v2)))
(defn euclidean-sq-1d "Squared Euclidean distance for doubles" ^double [^double v1 ^double v2] (m/sq (- v1 v2)))
(defn canberra-1d "Canberra distance for doubles" ^double [^double v1 ^double v2] (/ (m/abs (- v1 v2))
                                                                                  (+ (m/abs v1) (m/abs v2))))

(defn euclidean "Euclidean distance" ^double [v1 v2] (v/dist v1 v2))
(defn manhattan "Manhattan distance" ^double [v1 v2] (v/dist-abs v1 v2))
(defn chebyshev "Chebyshev distance" ^double [v1 v2] (v/dist-cheb v1 v2))
(defn canberra "Canberra distance" ^double [v1 v2] (v/dist-canberra v1 v2))
(defn earth-movers "Earth-Movers distance" ^double [v1 v2] (v/dist-emd v1 v2))
(defn euclidean-sq "Squared Euclidean distance" ^double [v1 v2] (v/dist-sq v1 v2))
(defn discrete "Discrete distance" ^double [v1 v2] (v/dist-discrete v1 v2))
(defn cosine "Cosine similarity" ^double [v1 v2] (v/sim-cos v1 v2))
(defn angular "Angular distance" ^double [v1 v2] (v/dist-ang v1 v2))
(defn weierstrass "Weierstrass metric" ^double [v1 v2] (m/acosh (- (* (m/sqrt (inc (v/magsq v1)))
                                                                   (m/sqrt (inc (v/magsq v2))))
                                                                (v/dot v1 v2))))

(defn haversine "Haversine distance" ^double [lat-lon1 lat-lon2] (m/haversine-dist lat-lon1 lat-lon2))

(defn periodic-euclidean "Euclidean distance with periodic boundaries"
  [v1 v2 w]
  (let [ws (vec (if (number? w) (repeat (count v1) w) w))
        a1 (mapv (fn [^double x ^double m] (m/mod x m)) (v/abs (v/sub v1 v2)) ws)
        a2 (v/sub ws a1)]
    (v/mag (v/emn a1 a2))))

(defn- covariance->inverse-matrix [cov] (-> (mat/mat->RealMatrix (m/seq->double-double-array cov))
                                            (mat/inverse)))
(def ^:private mem-covariance->inverse-matrix (memoize covariance->inverse-matrix))

(defn mahalanobis-sq
  "Squared Mahalonobis distance for given covariance (seq of seqs) matrix."
  ([cov] (let [S-1 (covariance->inverse-matrix cov)]
           (fn ^double [x mean]
             (let [diff (m/seq->double-array (v/sub x mean))]
               (v/dot diff (mat/mulv S-1 diff))))))
  (^double [x mean cov] (let [S-1 (mem-covariance->inverse-matrix cov)
                              diff (m/seq->double-array (v/sub x mean))]
                          (v/dot diff (mat/mulv S-1 diff)))))

(defn mahalanobis
  "Mahalonobis distance for given covariance (seq of seqs) matrix."
  ([cov] (let [m (mahalanobis-sq cov)]
           (fn ^double [x mean] (v/sqrt (m x mean)))))
  (^double [x mean cov] (m/sqrt (mahalanobis-sq x mean cov))))

(defn minkowski
  "Minkowski distance for order `p` and optional `weights` for each dimension."
  ([v1 v2 ^double p]
   (m/pow (reduce m/+ 0.0 (map (fn [^double xx ^double yy]
                                 (m/pow (m/abs (- xx yy)) p)) v1 v2)) (/ p)))
  ([v1 v2 ^double p weights]
   (m/pow (reduce m/+ 0.0 (map (fn [^double w ^double xx ^double yy]
                                 (* w (m/pow (m/abs (- xx yy)) p))) weights v1 v2)) (/ p)))
  ([^double p weights]
   (let [rp (/ p)]
     (fn ^double [v1 v2]
       (m/pow (reduce m/+ 0.0 (map (fn [^double w ^double xx ^double yy]
                                     (* w (m/pow (m/abs (- xx yy)) p))) weights v1 v2)) rp))))
  ([^double p]
   (let [rp (/ p)]
     (fn ^double [v1 v2]
       (m/pow (reduce m/+ 0.0 (map (fn [^double xx ^double yy]
                                     (m/pow (m/abs (- xx yy)) p)) v1 v2)) rp)))))
;;

(defn bound
  "Create relative distance function to be bound by a `diameter` (default: by `1.0`)."
  ([distance] (bound distance 1.0))
  ([distance ^double diameter]
   (fn ^double [v1 v2] (let [d (double (distance v1 v2))] (/ (* diameter d) (inc d))))))

(defn shoenberg-transform
  "Create relative distance function by applying Shoenberg transform."
  ([distance] (shoenberg-transform distance 1.0))
  ([distance ^double lambda]
   (fn ^double [v1 v2] (- 1.0 (m/exp (* -1.0 lambda ^double (distance v1 v2)))))))

(defn ->acm-distance
  "Create DistanceMeasure Apache Commons Math object based on given distance."
  [distance]
  (reify
    IFn (invoke [_ x y] (distance x y))
    DistanceMeasure (compute [_ x y] (distance x y))))

(m/unuse-primitive-operators)
