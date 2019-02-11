(ns fastmath.regression
  (:require [fastmath.core :as m]
            [fastmath.distance :as dist]
            [fastmath.rbf :as rbf]
            [fastmath.stats :as stat]
            [fastmath.kernel.mercer :as mercer])
  (:refer-clojure :exclude [test])
  (:import [clojure.lang IFn]
           [smile.regression Regression RegressionTrainer OLS$Trainer RLS$Trainer LASSO$Trainer RidgeRegression$Trainer
            RBFNetwork$Trainer SVR$Trainer RegressionTree$Trainer RandomForest$Trainer GradientTreeBoost$Trainer
            GradientTreeBoost$Loss GaussianProcessRegression$Trainer
            NeuralNetwork$Trainer NeuralNetwork$ActivationFunction]
           [smile.validation Validation RegressionMeasure MeanAbsoluteDeviation MSE RMSE RSS]
           [smile.math.rbf RadialBasisFunction]
           [smile.math.kernel MercerKernel]))

(defprotocol RegressionProto
  (backend [_])
  (model-native [_])
  (predict [_ v])
  (predict-all [_ vs])
  (train [_ x y])
  (data [_] [_ native?])
  (test [_] [_ tx ty])
  (cv-native [_] [_ params]))

(defmulti regression (fn [k & _] k))

(def ^:private regression-measures {:rmse (RMSE.)
                                    :mad (MeanAbsoluteDeviation.)
                                    :mse (MSE.)
                                    :rss (RSS.)})

(declare validate)

(defmethod regression :smile
  [_ ^RegressionTrainer trainer x y tx ty]
  (let [data-x (m/seq->double-double-array x)
        data-y (m/seq->double-array y)
        ^Regression regression-raw (.train trainer data-x data-y)
        predict-fn #(if (sequential? (first %))
                      (seq (.predict regression-raw (m/seq->double-double-array %)))
                      (.predict regression-raw (m/seq->double-array %)))]
    (reify
      IFn
      (invoke [_ v] (predict-fn v))

      RegressionProto
      (backend [_] :smile)
      (model-native [_] regression-raw)
      (predict [_ v] (predict-fn v))
      (predict-all [_ v] (predict-fn v))
      (train [_ x y] (regression :smile trainer x y))
      (data [_] [x y])
      (data [_ native?] (if native? [data-x data-y] [x y]))
      (test [c] (when (and tx ty) (validate c tx ty)))
      (test [c tx ty] (validate c tx ty))
      (cv-native [c] (cv-native c {}))
      (cv-native [_ {:keys [^int k type measure] :or {k 10 type :cv measure :rmse}}]
        (let [measure (if (contains? regression-measures measure) measure :rmse)
              ^RegressionMeasure measure-obj (regression-measures measure)]
          (case type
            :cv {measure (Validation/cv k trainer data-x data-y measure-obj)}
            :loocv {measure (Validation/loocv trainer data-x data-y measure-obj)}
            :bootstrap (let [b (Validation/bootstrap k trainer data-x data-y measure-obj)]
                         {measure (stat/mean b) :bootstrap (vec b)})))))))

(defmacro ^:private wrap-regression
  {:style/indent 3}
  [typ rname parameter instance & r]
  (let [[x y tx ty params] (map symbol ["x" "y" "tx" "ty" "params"])
        doc (str rname " regression. Backend library: " (name typ))]
    `(defn ~rname ~doc
       {:metadoc/categories #{:reg}}
       ([~x ~y] (~rname {} ~x ~y nil nil))
       ([~params ~x ~y] (~rname ~params ~x ~y nil nil))
       ([~x ~y ~tx ~ty] (~rname {} ~x ~y ~tx ~ty))
       ([~parameter ~x ~y ~tx ~ty] (regression ~typ ~instance ~x ~y ~tx ~ty ~@r)))))

(wrap-regression :smile ols {} (OLS$Trainer.))

(wrap-regression :smile rls {} (RLS$Trainer.))

(wrap-regression :smile lasso {:keys [lambda tolerance max-iters]
                               :or {lambda 10.0 tolerance 1.0e-3 max-iters 1000}}
  (-> (LASSO$Trainer. lambda)
      (.setTolerance tolerance)
      (.setMaxNumIteration max-iters)))

(wrap-regression :smile ridge {:keys [lambda]
                               :or {lambda 0.1}}
  (RidgeRegression$Trainer. lambda))

(wrap-regression :smile rbf-network {:keys [distance rbf number-of-basis normalize?]
                                     :or {distance dist/euclidean number-of-basis 10 normalize? false}}
  (let [cl (RBFNetwork$Trainer. distance)]
    (-> (cond
          (nil? rbf) cl
          (sequential? rbf) (.setRBF cl (into-array RadialBasisFunction rbf))
          :else (.setRBF cl rbf number-of-basis))
        (.setNormalized normalize?))))

(wrap-regression :smile svr {:keys [^MercerKernel kernel C eps tolerance]
                             :or {kernel (mercer/kernel :linear) C 1.0 eps 0.001 tolerance 1.0e-3}}
  (-> (SVR$Trainer. kernel eps C)
      (.setTolerance tolerance)))

(wrap-regression :smile regression-tree {:keys [^int max-nodes node-size]
                                         :or {max-nodes 100 node-size 2}}
  (let [features (count (first x))]
    (-> (RegressionTree$Trainer. features max-nodes)
        (.setNodeSize node-size))))

(wrap-regression :smile random-forest {:keys [number-of-trees mtry node-size max-nodes subsample]
                                       :or {number-of-trees 500 node-size 2 max-nodes 100 subsample 1.0}}
  (let [mtry (or mtry (max 1 (/ (count (first x)) 3)))]
    (-> (RandomForest$Trainer. number-of-trees)
        (.setNumRandomFeatures mtry)
        (.setMaxNodes max-nodes)
        (.setSamplingRates subsample)
        (.setNodeSize node-size))))

(def loss-strategies {:least-squares GradientTreeBoost$Loss/LeastSquares
                      :least-absolute-deviation GradientTreeBoost$Loss/LeastAbsoluteDeviation
                      :huber GradientTreeBoost$Loss/Huber})

(def ^{:doc "List of loss for Gradient Tree Boost algorithm"} loss-list (keys loss-strategies))

(wrap-regression :smile gradient-tree-boost {:keys [loss number-of-trees shrinkage max-nodes subsample]
                                             :or {loss :least-squares number-of-trees 500 shrinkage 0.005
                                                  max-nodes 6 subsample 0.7}}
  (-> (GradientTreeBoost$Trainer. number-of-trees)
      (.setLoss (or (loss-strategies loss) GradientTreeBoost$Loss/LeastSquares))
      (.setMaxNodes max-nodes)
      (.setShrinkage shrinkage)
      (.setSamplingRates subsample)))

(wrap-regression :smile gaussian-process {:keys [^MercerKernel kernel lambda]
                                          :or {kernel (mercer/kernel :linear) lambda 0.5}}
  (-> (GaussianProcessRegression$Trainer. kernel lambda)))

(def activation-functions {:logistic-sigmoid NeuralNetwork$ActivationFunction/LOGISTIC_SIGMOID
                           :tanh NeuralNetwork$ActivationFunction/TANH})

(wrap-regression :smile neural-net {:keys [activation-function layers learning-rate momentum weight-decay number-of-epochs]
                                    :or {error-function :logistic-sigmoid learning-rate 0.1 momentum 0.0 weight-decay 0.0 number-of-epochs 25}}
  (let [layers (into-array Integer/TYPE (cons (count (first x)) (conj (vec layers) 1)))]
    (-> (if-not activation-function
          (NeuralNetwork$Trainer. layers)
          (NeuralNetwork$Trainer. (or (activation-functions activation-function)
                                      NeuralNetwork$ActivationFunction/LOGISTIC_SIGMOID) layers))
        (.setLearningRate learning-rate)
        (.setMomentum momentum)
        (.setWeightDecay weight-decay)
        (.setNumEpochs number-of-epochs))))

(defn validate
  "Validate data against trained regression."
  ([model] (test model))
  ([model tx ty]
   (let [pred (predict-all model tx)
         pred-arr (m/seq->double-array pred)
         truth (m/seq->double-array ty)]
     {:prediction pred
      :stats (into {} (map (fn [[k ^RegressionMeasure v]]
                             [k (.measure v truth pred-arr)]) regression-measures))})))
