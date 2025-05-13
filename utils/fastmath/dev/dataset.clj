(ns fastmath.dev.dataset
  (:require [clojure.java.io :as io]
            [clojure.data.csv :as csv]))

(defn transform [spec data]
  (if (sequential? spec)
    (let [[nm f] spec] [nm (f data)])
    [spec (read-string data)]))
(defn zip [header row] (into {} (map transform header row)))
(defn parse-rows  [header rows]  (map (partial zip header) rows))
(defn read-csv
  [fname header]
  (->> (io/resource fname)
       (slurp)
       (csv/read-csv)
       (rest)
       (parse-rows header)))

(defn data->fn
  [data]
  (memoize (fn
             ([] data)
             ([selector] (map #(get % selector) data))
             ([selector filter-pred] (map #(get % selector) (filter filter-pred data))))))

(defn by
  ([data f]
   (group-by f (data)))
  ([data f selector]
   (update-vals (group-by f (data)) (partial map #(get % selector)))))

(def iris (data->fn (read-csv "iris.csv"
                            [:sepal-length :sepal-width :petal-length :petal-width [:species keyword]])))
(def mtcars (data->fn (read-csv "mtcars.csv"
                              [[:name identity] :mpg :cyl :disp :hp :drat :wt :qsec :vs :am :gear :carb])))

(def winequality (data->fn (read-csv "winequality-white.csv"
                                   ["fixed acidity","volatile acidity","citric acid","residual sugar","chlorides","free sulfur dioxide","total sulfur dioxide","density","pH","sulphates","alcohol","quality"])))
