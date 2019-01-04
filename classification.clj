(ns fastmath.classification
  (:require [fastmath.random :as r]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [fastmath.stats :as stat])
  (:import [smile.classification KNN AdaBoost FLD SVM QDA LDA RandomForest LogisticRegression]))

(def data (->> (io/resource "iris.csv")
               (io/reader)
               (csv/read-csv)
               (drop 1)
               (map (fn [v]
                      (let [[x y z w nm] (map read-string v)]
                        [[x y z w] (case (str nm)
                                     "setosa"     0
                                     "versicolor" 1
                                     "virginica"  2)])))))

(defn prepare-data [in]
  {:data (m/seq->double-double-array (map first in))
   :cats (int-array (map second in))})

(defn simulate
  []
  (let [split-point (r/irand (* 0.8 (count data)) (* 0.9 (count data)))
        split-set (let [[tr te] (split-at split-point (shuffle data))]
                    {:train-data (prepare-data tr)
                     :test-data (prepare-data te)})
        train-data (:train-data split-set)
        test-data (:test-data split-set)
        cl (LogisticRegression. (:data train-data) (:cats train-data) 100)
        pred (map #(.predict cl %) (:data test-data))]
    (/ (count (filter zero? (map compare pred (:cats test-data))))
       (count (:cats test-data)))))

(org.slf4j.LoggerFactory/getLogger (org.slf4j.Logger/ROOT_LOGGER_NAME))

(simulate)

(stat/mean (repeatedly 1000 simulate))

(comment

  (def data1 (repeatedly 1000 r/grand))
  (def data2 (repeatedly 1000 #(r/grand 2 10)))

  (def data (map vector data1 data2))

  (def cats (map (fn [[v1 v2]]
                   (if (or (and (pos? v1) (pos? v2))
                           (and (neg? v1) (neg? v2)))
                     0 1)) data))

  (def test-data (map m/seq->double-array  [[2 2] [-2 -2] [2 -2] [-2 2]]))


  (def cl (KNN/learn (m/seq->double-double-array data) (int-array cats)))
  (def cl (FLD. (m/seq->double-double-array data) (int-array cats)))
  (def cl (AdaBoost. (m/seq->double-double-array data) (int-array cats) 100))
  (def cl (QDA. (m/seq->double-double-array data) (int-array cats)))
  (def cl (LDA. (m/seq->double-double-array data) (int-array cats)))
  (def cl (RandomForest. (m/seq->double-double-array data) (int-array cats) 100))
  (def cl (LogisticRegression. (m/seq->double-double-array data) (int-array cats)))

  (map #(.predict cl %) test-data)


  (.predict fld (m/seq->double-array [10 10]))


  (.predict adaboost (m/seq->double-array [1 1]))



  (.predict qda (m/seq->double-array [10 -10]))

  (let [v [1 2 1 2 12 2 2 1 2 2 1 2 1]
        f (frequencies v)
        total (reduce + (vals f))]
    (into {} (map (fn [[k v]] [k (/ v total)]) f)))

  (defn roulette-wheel
    "Perform a roulette wheel selection given a list of frequencies"
    [freqs]
    (let [nfreqs (count freqs)
          tot (reduce + freqs)]
      (if (= tot 0)
        nil
        (let [dist (map #(/ % tot) freqs)
              rval (double (rand))]
          (loop [acc 0, i 0]
            (let [lb acc, ub (+ acc (nth dist i))]
              (cond (>= (+ i 1) nfreqs) i
                    (and (>= rval lb) (< rval ub)) i
                    :else (recur ub (+ i 1)))))))))

  (roulette-wheel [1 2 3 4 3 22])

  )
