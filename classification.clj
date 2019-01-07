(ns fastmath.classification
  (:require [fastmath.random :as r]
            [fastmath.core :as m]
            [fastmath.vector :as v]
            [clojure.java.io :as io]
            [clojure.data.csv :as csv]
            [fastmath.stats :as stat]
            [fastmath.distance :as dist])
  (:import [clojure.lang IFn]
           [smile.classification SoftClassifier KNN$Trainer AdaBoost$Trainer FLD$Trainer QDA$Trainer LDA$Trainer RandomForest$Trainer LogisticRegression ClassifierTrainer]
           [smile.math.distance EuclideanDistance]
           [smile.validation Accuracy Fallout FDR FMeasure Precision Recall Sensitivity Specificity ClassificationMeasure Validation]))

(set! *warn-on-reflection* false)

(defprotocol ClassificationProto
  (predict [_ v] [_ v posteriori?])
  (predict-all [_ v] [_ v posteriori?])
  (train [_ x y])
  (data [_])
  (labels [_])
  (validate [_ tx ty] [_ tx ty measures])
  (loocv [_] [_ measures] [_ tx ty] [_ tx ty measures])
  (cv [_ k] [_ k measures] [_ k tx ty] [_ k tx ty measures])
  (bootstrap [_ k] [_ k measures] [_ k tx ty] [_ k tx ty measures]))

(defmulti ^:private prepare-data (fn [k & _] k))

(defmethod prepare-data :smile
  ([_ x y labels labels->int]
   [(m/seq->double-double-array x) (int-array (map labels->int y)) labels labels->int])
  ([_ x y]
   (let [ydata (map vector (sort (distinct y)) (range))
         labels (mapv first ydata)
         labels->int (into {} ydata)]
     (prepare-data :smile x y labels labels->int))))

(def ^:private measure-objs {:accuracy (Accuracy.)
                             :fallout (Fallout.)
                             :fdr (FDR.)
                             :fscore (FMeasure.)
                             :precision (Precision.)
                             :recall (Recall.)
                             :sensitivity (Sensitivity.)
                             :specificity (Specificity.)})

(def ^:const ^:private classification-measure-class (Class/forName "smile.validation.ClassificationMeasure"))

(defn- measures->array
  [measures]
  (into-array classification-measure-class (map measure-objs measures)))

(defn- validate-smile
  [measure y-truth y-pred]
  (when-let [m (measure-objs measure)]
    (.measure ^ClassificationMeasure m y-truth y-pred)))

(defmulti ^:private trainer (fn [k c x y] k))

(defmethod trainer :smile
  [_ ^ClassifierTrainer classifier x y]
  (let [[data int-labels labels labels->int] (prepare-data :smile x y)
        trainer (.train classifier data int-labels)
        predict-raw #(.predict trainer (m/seq->double-array %))
        predict-fn (comp labels predict-raw)
        predict-fn-posteriori (if (instance? SoftClassifier trainer)
                                #(let [posteriori (double-array (count labels))]
                                   [(labels (.predict trainer (m/seq->double-array %) posteriori)) (seq posteriori)])
                                predict-fn)]
    (reify
      
      IFn
      (invoke [_ v] (predict-fn v))
      (invoke [_ v posteriori?] (if posteriori? (predict-fn-posteriori v) (predict-fn v)))

      ClassificationProto
      (predict [_ v] (predict-fn v))
      (predict [_ v posteriori?] (if posteriori? (predict-fn-posteriori v) (predict-fn v)))
      (predict-all [_ vs] (map predict-fn vs))
      (predict-all [_ vs posteriori?] (if posteriori? (map predict-fn-posteriori vs) (map predict-fn vs)))
      (train [_ x y] (trainer :smile classifier x y))
      (data [_] x)
      (labels [_] labels)

      (validate [_ tx ty]
        (let [[test-data test-int-labels] (prepare-data :smile tx ty labels labels->int)
              result (map predict-raw test-data)
              invalid (->> (map vector x test-int-labels result)
                           (filter (fn [[_ a b]] (not= a b))))
              cnt (count invalid)]
          {:invalid {:count cnt
                     :data (map first invalid)
                     :predictions (map (comp labels #(nth % 2)) invalid)
                     :truth (map (comp labels second) invalid)}
           :raw {:predictions (int-array result)
                 :truth test-int-labels}
           :stats {:accuracy (double (- 1.0 (/ cnt (count test-data))))}}))

      (validate [cl tx ty measures]
        (let [v (validate cl tx ty)
              truth (get-in v [:raw :truth])
              predictions (get-in v [:raw :predictions])]
          (reduce (fn [c m] (update c :stats assoc m (validate-smile m truth predictions))) v measures)))

      (loocv [_] (Validation/loocv classifier data int-labels))
      (loocv [_ measures] (into {} (map vector measures (Validation/loocv classifier data int-labels (measures->array measures)))))
      (loocv [_ tx ty]
        (let [[test-data test-int-labels] (prepare-data :smile tx ty labels labels->int)]
          (Validation/loocv classifier test-data test-int-labels)))
      (loocv [_ tx ty measures]
        (let [[test-data test-int-labels] (prepare-data :smile tx ty labels labels->int)]
          (into {} (map vector measures (Validation/loocv classifier test-data test-int-labels (measures->array measures))))))

      (cv [_ k] (Validation/cv k classifier data int-labels))
      (cv [_ k measures] (into {} (map vector measures (Validation/cv k classifier data int-labels (measures->array measures)))))
      (cv [_ k tx ty]
        (let [[test-data test-int-labels] (prepare-data :smile tx ty labels labels->int)]
          (Validation/cv k classifier test-data test-int-labels)))
      (cv [_ k tx ty measures]
        (let [[test-data test-int-labels] (prepare-data :smile tx ty labels labels->int)]
          (into {} (map vector measures (Validation/cv k classifier test-data test-int-labels (measures->array measures))))))

      (bootstrap [_ k] (Validation/bootstrap k classifier data int-labels))
      (bootstrap [_ k measures] (m/double-double-array->seq (Validation/bootstrap k classifier data int-labels (measures->array measures))))
      (bootstrap [_ k tx ty]
        (let [[test-data test-int-labels] (prepare-data :smile tx ty labels labels->int)]
          (Validation/bootstrap k classifier test-data test-int-labels)))
      (bootstrap [_ k tx ty measures]
        (let [[test-data test-int-labels] (prepare-data :smile tx ty labels labels->int)]
          (m/double-double-array->seq (Validation/bootstrap k classifier test-data test-int-labels (measures->array measures)))))

      )))

(defn knn
  ([x y distance k] (trainer :smile (KNN$Trainer. distance k) x y))
  ([x y k] (knn x y (EuclideanDistance.) k))
  ([x y] (knn x y 1)))

(defn ada-boost
  ([x y] (trainer :smile (AdaBoost$Trainer.) x y))
  ([x y number-of-trees] (let [^AdaBoost$Trainer t (AdaBoost$Trainer.)]
                           (trainer :smile (.setNumTrees t number-of-trees) x y)))
  ([x y number-of-trees max-nodes] (let [^AdaBoost$Trainer t (AdaBoost$Trainer.)]
                                     (trainer :smile (-> t (.setNumTrees number-of-trees)
                                                         (.setMaxNodes max-nodes)) x y))))

(defn fld
  ([x y] (trainer :smile (FLD$Trainer.) x y))
  ([x y dimensionality] (let [^FLD$Trainer t (FLD$Trainer.)]
                          (trainer :smile (.setDimension t dimensionality) x y)))
  ([x y dimensionality tolerance] (let [^FLD$Trainer t (FLD$Trainer.)]
                                    (trainer :smile (-> t (.setDimension dimensionality)
                                                        (.setTolerance tolerance)) x y))))

(defn qda
  ([x y] (trainer :smile (QDA$Trainer.) x y))
  ([x y priori] (let [^QDA$Trainer t (QDA$Trainer.)]
                  (trainer :smile (.setPriori t (m/seq->double-array priori)) x y)))
  ([x y priori tolerance] (let [^QDA$Trainer t (QDA$Trainer.)]
                            (trainer :smile (-> t (.setPriori (m/seq->double-array priori))
                                                (.setTolerance tolerance)) x y))))


(defn lda
  ([x y] (trainer :smile (LDA$Trainer.) x y))
  ([x y priori] (let [^LDA$Trainer t (LDA$Trainer.)]
                  (trainer :smile (.setPriori t (m/seq->double-array priori)) x y)))
  ([x y priori tolerance] (let [^LDA$Trainer t (LDA$Trainer.)]
                            (trainer :smile (-> t (.setPriori (m/seq->double-array priori))
                                                (.setTolerance tolerance)) x y))))


(def iris-data (->> (io/resource "iris.csv")
                    (io/reader)
                    (csv/read-csv)
                    (drop 1)
                    (map (fn [v]
                           (let [[x y z w nm] (map read-string v)]
                             [[x y z w] (str nm)])))))

(def split
  (let [split-point (r/irand (* 0.8 (count iris-data)) (* 0.9 (count iris-data)))
        iris-shuffled (shuffle iris-data)
        iris-v (map first iris-shuffled)
        iris-l (map second iris-shuffled)
        [dd dt] (split-at split-point iris-v)
        [ld lt] (split-at split-point iris-l)]
    {:data dd
     :labels ld
     :test-data dt
     :test-labels lt}))

(def knn-cl (knn (:data split) (:labels split) 3))
(def adaboost-cl (ada-boost (:data split) (:labels split) 30 10))
(def fld-cl (fld (:data split) (:labels split)))
(def qda-cl (qda (:data split) (:labels split) [0.05 0.05 0.9]))
(def lda-cl (lda (:data split) (:labels split)))

(bootstrap knn-cl 2 [:accuracy])

(validate knn-cl (:test-data split) (:test-labels split))
(validate adaboost-cl (:test-data split) (:test-labels split))
(validate qda-cl (:test-data split) (:test-labels split))
(validate lda-cl (:test-data split) (:test-labels split))

(type (object-array [1 2]))
;; => [Ljava.lang.Object;

(defn prepare-data [in]
  {:data (m/seq->double-double-array (map first in))
   :cats (int-array (map second in))})

(defn simulate
  []
  (let [split-point (r/irand (* 0.8 (count iris-data)) (* 0.9 (count iris-data)))
        split-set (let [[tr te] (split-at split-point (shuffle iris-data))]
                    {:train-data (prepare-data tr)
                     :test-data (prepare-data te)})
        train-data (:train-data split-set)
        test-data (:test-data split-set)
        cl (LogisticRegression. (:data train-data) (:cats train-data) 100)
        pred (map #(.predict cl %) (:data test-data))]
    (/ (count (filter zero? (map compare pred (:cats test-data))))
       (count (:cats test-data)))))


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
