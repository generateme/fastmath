(ns fastmath.clustering-examples
  (:require [fastmath.clustering :refer :all]
            [metadoc.examples :refer :all]
            [fastmath.random :as r]
            [fastmath.distance :as d]))

(add-examples clustering-methods-list
  (example "List of methods" clustering-methods-list))

(add-examples k-means
  (example (k-means [1 2 3 -1 -1 2 -1 11 111] 4))
  (example "Clusters group into separate maps." (regroup (k-means [1 2 3 -1 -1 2 -1 11 111] 4)))
  (example "Use as predictor" (let [cl (k-means [1 2 3 -1 -1 2 -1 11 111] 4)]
                                [(cl -1) (cl 10) (cl 100) (cl 1) (cl -1000) (cl 1000)])))

(add-examples g-means
  (example "Expect 2 clusters, uniform distribution."
    ((juxt :clusters :sizes :representatives) (g-means (repeatedly 100 #(r/randval (r/drand) (r/drand 5 6))) 4))))

(add-examples x-means
  (example "Expect 2 clusters, gaussian distribution."
    ((juxt :clusters :sizes :representatives) (x-means (repeatedly 10000 #(r/randval (r/grand) (r/grand 5 1.0))) 4))))

(add-examples deterministic-annealing
  (example (map (fn [m] (dissoc m :data))
                (-> (repeatedly 1000 #(vector (r/randval (r/grand) (r/grand 5 1.0))
                                              (r/randval (r/grand) (r/grand 5 1.0))))
                    (deterministic-annealing 4 0.5)
                    (regroup)))))

(add-examples denclue
  (example-session "Expect 2 clusters, uniform distribution."
    ((juxt :clusters :sizes :representatives) (denclue (repeatedly 100 #(r/randval (r/drand) (r/drand 5 6))) 1 10))
    (map (fn [m] (dissoc m :data)) (regroup (denclue (repeatedly 1000 #(r/randval 0.1 (r/drand) (r/drand 5 6))) 1 10)))))

(add-examples clarans
  (example (dissoc (clarans (repeatedly 1000 #(r/randval 0.1 (r/irand -10 10) (r/irand 100 150))) d/chebyshev 2) :data :clustering :obj :predict)))

(add-examples dbscan
  (example "3d vectors" (dissoc (dbscan (repeatedly 5000 #(vector (r/randval 0.1 (r/irand -10 10) (r/irand 100 150))
                                                                  (r/randval (r/irand -10 10) (r/irand 100 150))
                                                                  (r/randval (r/irand -10 10) (r/irand 100 150)))) 10 20) :data :clustering :obj :predict)))

(add-examples mec
  (example "2d vectors" (dissoc (mec (repeatedly 5000 #(vector (r/randval 0.1 (r/irand -10 10) (r/irand 100 150))
                                                               (r/randval (r/irand -10 10) (r/irand 100 150)))) d/manhattan 8 20) :data :clustering :obj :predict)))

(add-examples spectral
  (example "2d vectors" (dissoc (spectral (repeatedly 500 #(vector (r/randval 0.1 (r/irand -10 10) (r/irand 100 150))
                                                                   (r/randval (r/irand -10 10) (r/irand 100 150)))) 4 1) :data :clustering :obj :predict)))

(add-examples spectral
  (example "2d vectors" (dissoc (lloyd (repeatedly 500 #(vector (r/randval 0.1 (r/irand -10 10) (r/irand 100 150))
                                                                (r/randval (r/irand -10 10) (r/irand 100 150)))) 4 1) :data :clustering :obj :predict)))


(add-examples regroup
  (example-session "Result of clustering with regrouping"
    (k-means [1 2 3 -1 -1 2 -1 11 111] 7)
    (regroup (k-means [1 2 3 -1 -1 2 -1 11 111] 7))
    (count (regroup (k-means [1 2 3 -1 -1 2 -1 11 111] 7)))))
