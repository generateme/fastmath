(ns codrna-classification
  (:require [fastmath.classification :as cl]
            [clj-boost.core :as xgboost]
            [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [clojure.string :as s]))


(defn split [x]  (s/split x #":"))
(defn rec->vec
  ""
  [xs features]
  (let [a (double-array features)]
    (doseq [x xs
            :let [[id v] (split x)]]
      (aset a (dec (read-string id)) ^double (read-string v)))
    (seq a)))

(defn read-data
  [name]
  (map #(s/split % #"\s+") (line-seq (io/reader (io/resource name)))))

(defn process-data
  [name features]
  (let [lst (read-data name)]
    {:labels (mapv first lst)
     :data (map #(rec->vec (rest %) features) lst)}))

(def training-data (process-data "cod-rna.txt" 8))
(def test-data (process-data "cod-rna-test.txt" 8))

(def ttt (cl/lasvm (take 10000 (:data training-data)) (take 10000 (:labels training-data))))

(.valueOf ttt (double-array (first (:data test-data))))

(def knn-cl (cl/knn (:data training-data) (:labels training-data)))
(def vknn (cl/validate knn-cl (:data test-data) (:labels test-data)))

(def xgboost-cl (cl/xgboost {:params {:objective "binary:gistic"}} (:data training-data) (:labels training-data) (:data test-data) (:labels test-data)))
(def vknn (cl/validate xgboost-cl))

(def ll-cl (cl/liblinear {:solver :l1r-lr} (:data training-data) (:labels training-data) (:data test-data) (:labels test-data)))
(def vknn (cl/validate ll-cl))

(xgboost/predict (cl/model-native xgboost-cl) (xgboost/dmatrix [(vec (first (:data test-data)))]))

(xgboost-cl (first (:data test-data)) true)



(last (:labels test-data))

(cl/labels xgboost-cl)

(cl/cv-native xgboost-cl)

(:stats vknn)

(get-in vknn [:invalid :count])

(fn [d & _] (type d))
(defmulti dmatrix (fn [d & _] (cond
                                (seqable? d) :seq
                                (map? d) :map
                                :else (type d))))
