(ns fastmath.ml.clustering.hierarchical-test
  (:require [fastmath.ml.clustering.hierarchical :as sut]
            [clojure.test :as t]))

;; The public functions are tested in fastmath.ml.clustering-test against R, SciPy and Julia,
;; here only the condensed layout shared by all algorithms.

(t/deftest condensed-layout
  (t/testing "pairs i < j in row order, as R as.vector(as.dist(m)): d(1,2) d(1,3) d(1,4) d(2,3) d(2,4) d(3,4)"
    (let [m [[0 12 13 14] [12 0 23 24] [13 23 0 34] [14 24 34 0]]]
      (t/is (= [12.0 13.0 14.0 23.0 24.0 34.0]
               (vec (sut/matrix-distances (mapv double-array m)))))))
  (t/testing "index of every pair"
    (doseq [n [2 3 7]]
      (t/is (= (range (sut/condensed-size n))
               (for [i (range n) j (range (inc i) n)] (sut/condensed-index n i j))))))
  (t/testing "the largest size fits an array"
    (t/is (= 2147450880 (sut/condensed-size 65536)))
    (t/is (= 2147450879 (sut/condensed-index 65536 65534 65535))))
  (t/testing "distance-at is symmetric with zeros on the diagonal"
    (let [d (sut/matrix-distances (mapv double-array [[0 1 2] [1 0 3] [2 3 0]]))]
      (t/is (= [[0.0 1.0 2.0] [1.0 0.0 3.0] [2.0 3.0 0.0]]
               (mapv (fn [i] (mapv #(sut/distance-at d 3 i %) (range 3))) (range 3)))))))
