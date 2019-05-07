(ns fastmath.classification-examples
  (:refer-clojure :exclude [test])
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [metadoc.examples :refer :all]
            [fastmath.classification :refer :all]
            [fastmath.distance :as dist]
            [fastmath.kernel :as k]))

(def iris-data (->> (io/resource "iris.csv")
                    (io/reader)
                    (csv/read-csv)
                    (drop 1)
                    (map (fn [v]
                           (let [[x y z w nm] (map read-string v)]
                             [[x y z w] (str nm)])))))

(def split
  (let [split-point (* 0.7 (count iris-data))
        iris-shuffled (shuffle iris-data)
        iris-v (map first iris-shuffled)
        iris-l (map second iris-shuffled)
        [dd dt] (split-at split-point iris-v)
        [ld lt] (split-at split-point iris-l)]
    {:data dd
     :labels ld
     :test-data dt
     :test-labels lt}))

(def train-data (:data split))
(def train-labels (:labels split))
(def test-data (:test-data split))
(def test-labels (:test-labels split))

(add-examples accuracy
  (example (accuracy [1 2 3 4 1 2 3 4] [1 2 4 4 2 4 4 3])))

(add-examples activation-functions-list (example "Names" activation-functions-list))

(add-examples ada-boost
  (example (let [cl (ada-boost train-data train-labels)]
             (select-keys (validate cl test-data test-labels) [:invalid :stats]))))

(add-examples bayes-models-list (example "Names" bayes-models-list))

(add-examples backend
              (example-session "Usage"
                               (backend (knn test-data test-labels))
                               #_    (backend (xgboost test-data test-labels))
                               (backend (liblinear test-data test-labels))))

(add-examples model-native
              (example-session "Usage"
                               (model-native (knn train-data train-labels))
                               #_(model-native (xgboost train-data train-labels))
                               ))

(add-examples data
  (example-session "Usage"
    (data (knn [[1 2] [3 2]] [0 1]))
    (data (knn [[1 2] [3 2]] [0 1]) true)))

(add-examples labels
  (example (labels (knn train-data train-labels))))

(add-examples predict
  (example-session "Usage"
    (predict (knn train-data train-labels) [1 2 3 4])
    (predict (ada-boost train-data train-labels) [1 2 3 4] true)))

(add-examples predict-all
  (example-session "Usage"
    (predict-all (knn train-data train-labels) (take 3 test-data))
    (predict-all (ada-boost train-data train-labels) (take 3 test-data) true)))

(add-examples test
  (example "Use data provided during training" (test (knn train-data train-labels test-data test-labels)))
  (example "No test data provided" (test (knn train-data train-labels)))
  (example "Test other data" (test (knn train-data train-labels) [[1 2 3 4] [6 4 3 2]] ["setosa" "setosa"])))

(add-examples validate
  (example "Use data provided during training" (validate (qda train-data train-labels test-data test-labels)))
  (example "No test data provided" (validate (qda train-data train-labels)))
  (example "Test other data" (validate (qda train-data train-labels) [[1 2 3 4] [6 4 3 2]] ["virginica" "setosa"])))

(add-examples train
  (example "Train new data" (train (knn train-data train-labels) test-data test-labels)))

(add-examples decision-tree
  (example (let [cl (decision-tree train-data train-labels)]
             (select-keys (validate cl test-data test-labels) [:invalid :stats]))))

(add-examples fld
  (example (let [cl (fld train-data train-labels)]
             (select-keys (validate cl test-data test-labels) [:invalid :stats]))))

(add-examples lda
  (example (let [cl (lda train-data train-labels)]
             (select-keys (validate cl test-data test-labels) [:invalid :stats]))))

(add-examples logistic-regression
  (example (let [cl (logistic-regression train-data train-labels)]
             (select-keys (validate cl test-data test-labels) [:invalid :stats]))))

(add-examples loocv
  (example (loocv (knn train-data train-labels))))

(add-examples cv
  (example-session "Usage"
    (cv (knn train-data train-labels))
    (cv (knn train-data train-labels) 100)
    (cv (knn train-data train-labels) 4)
    (cv (knn train-data train-labels) 1)))

(add-examples cv-native
              (example-session "Usage"
                               (cv-native (knn test-data test-labels))
                               (cv-native (knn test-data test-labels) {:type :loocv})
                               (cv-native (knn test-data test-labels) {:type :bootstrap})
                               #_    (cv-native (xgboost test-data test-labels))
                               (cv-native (liblinear test-data test-labels))))

(add-examples bootstrap
  (example-session "Usage"
    (select-keys (bootstrap (knn train-data train-labels)) [:accuracy])
    (select-keys (bootstrap (fld train-data train-labels)) [:accuracy])
    (bootstrap (decision-tree train-data train-labels) 4)
    (select-keys (bootstrap (knn train-data train-labels) 500) [:accuracy])
    (bootstrap (neural-net {:layers [100 50]} train-data train-labels) 5)))


(add-examples gradient-tree-boost
  (example (let [cl (gradient-tree-boost train-data train-labels)]
             (select-keys (validate cl test-data test-labels) [:invalid :stats]))))

(add-examples naive-bayes
  (example (let [cl (naive-bayes train-data train-labels)]
             (select-keys (validate cl test-data test-labels) [:invalid :stats]))))

(add-examples neural-net
  (example (let [cl (neural-net {:layers [100 100]} train-data train-labels)]
             (select-keys (validate cl test-data test-labels) [:invalid :stats])))
  (example "Bad model" (let [cl (neural-net {:layers [3 3 3] :number-of-epochs 3} train-data train-labels)]
                         (select-keys (validate cl test-data test-labels) [:invalid :stats]))))

(add-examples qda
  (example (let [cl (qda train-data train-labels)]
             (select-keys (validate cl test-data test-labels) [:invalid :stats]))))

(add-examples random-forest
  (example (let [cl (random-forest train-data train-labels)]
             (select-keys (validate cl test-data test-labels) [:invalid :stats]))))

(add-examples rbf-network
              (example (let [cl (rbf-network train-data train-labels)]
                         (select-keys (validate cl test-data test-labels) [:invalid :stats])))
              (example "Custom rbfs" (let [cl (rbf-network {:rbf (take 5 (cycle [(k/rbf :linear) (k/rbf :wendland-20)]))} train-data train-labels)]
                                       (select-keys (validate cl test-data test-labels) [:invalid :stats]))))


(add-examples rda
  (example (let [cl (rda train-data train-labels)]
             (select-keys (validate cl test-data test-labels) [:invalid :stats]))))

(add-examples error-functions-list (example "Names" error-functions-list))

(add-examples split-rules-list (example "Names" split-rules-list))

(add-examples knn
  (example (let [cl (knn train-data train-labels)]
             (select-keys (validate cl test-data test-labels) [:invalid :stats])))
  (example "Different distance" (let [cl (knn {:distance dist/cosine} train-data train-labels)]
                                  (select-keys (validate cl test-data test-labels) [:invalid :stats]))))

(add-examples svm
              (example (let [cl (svm train-data train-labels)]
                         (select-keys (validate cl test-data test-labels) [:invalid :stats])))
              (example "Different kernel" (let [cl (svm {:kernel (k/kernel :gaussian 1) :epochs 10} train-data train-labels)]
                                            (select-keys (validate cl test-data test-labels) [:invalid :stats]))))

#_(add-examples xgboost
                (example-session "Usage"
                                 (let [cl (xgboost train-data train-labels test-data test-labels)]
                                   (select-keys (validate cl test-data test-labels) [:invalid :stats]))
                                 (let [cl (xgboost {:params {:eta 0.999
                                                             :max_depth 1
                                                             :grow_policy "lossguide"
                                                             :alpha 0.5}
                                                    :rounds 2} train-data train-labels test-data test-labels)]
                                   (select-keys (validate cl test-data test-labels) [:invalid :stats]))))

(add-examples liblinear
              (example (let [cl (liblinear train-data train-labels)]
                         (select-keys (validate cl test-data test-labels) [:invalid :stats])))
              (example "Different solver" (let [cl (liblinear {:solver :l1r-lr :C 0.5} train-data train-labels)]
                                            (select-keys (validate cl test-data test-labels) [:invalid :stats]))))


(add-examples confusion-map
  (example (let [cl (knn train-data train-labels)
                 pred (predict-all cl test-data)]
             (confusion-map test-labels pred))))

(add-examples liblinear-solver-list (example "Names" liblinear-solver-list))
(add-examples multiclass-strategies-list (example "Names" multiclass-strategies-list))
