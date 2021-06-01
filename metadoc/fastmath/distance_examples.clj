(ns fastmath.distance-examples
  (:require [metadoc.examples :refer :all]
            [fastmath.distance :refer :all]
            [fastmath.stats :as stats]))

(add-examples euclidean (example (euclidean [1 2 4 4] [-1 3 4 -5])))
(add-examples manhattan (example (manhattan [1 2 4 4] [-1 3 4 -5])))
(add-examples chebyshev (example (chebyshev [1 2 4 4] [-1 3 4 -5])))
(add-examples correlation (example (correlation [1 2 4 4] [-1 3 4 -5])))
(add-examples canberra (example (canberra [1 2 4 4] [-1 3 4 -5])))
(add-examples earth-movers (example (earth-movers [1 2 4 4] [-1 3 4 -5])))
(add-examples euclidean-sq (example (euclidean-sq [1 2 4 4] [-1 3 4 -5])))
(add-examples discrete (example (discrete [1 2 4 4] [-1 3 4 -5])))
(add-examples cosine (example (cosine [1 2 4 4] [-1 3 4 -5])))
(add-examples angular (example (angular [1 2 4 4] [-1 3 4 -5])))
(add-examples jensen-shannon (example (jensen-shannon [1 2 4 4] [1 3 4 5])))

(add-examples make-mahalanobis
  (example (let [cov (stats/covariance-matrix [[1 2 4 4] [-1 3 4 -5]])
                 d (make-mahalanobis cov)]
             (d [1 2] [3 4]))))

(add-examples make-minkowski
  (example-session "Usage"
    (let [d (make-minkowski 0.5)] (d [1 2 4 4] [-1 3 4 -5]))
    (let [d (make-minkowski 3.0)] (d [1 2 4 4] [-1 3 4 -5]))
    (let [d (make-minkowski 2.0 [0.1 0.2 0.4 0.3])] (d [1 2 4 4] [-1 3 4 -5]))))
