(ns fastmath.regression
  (:require [fastmath.core :as m]
            [fastmath.distance :as dist]
            [fastmath.kernel :as k]
            [fastmath.stats :as stats]
            [fastmath.random :as r]
            [fastmath.vector :as v])
  (:refer-clojure :exclude [test])
  (:import [clojure.lang IFn]
           [smile.regression Regression RegressionTrainer OLS$Trainer RLS$Trainer LASSO$Trainer RidgeRegression$Trainer
            RBFNetwork$Trainer SVR$Trainer RegressionTree$Trainer RandomForest$Trainer GradientTreeBoost$Trainer
            GradientTreeBoost$Loss GaussianProcessRegression$Trainer
            NeuralNetwork$Trainer NeuralNetwork$ActivationFunction
            ElasticNet]
           [smile.validation Validation RegressionMeasure MeanAbsoluteDeviation MSE RMSE RSS]
           [org.apache.commons.math3.linear MatrixUtils CholeskyDecomposition Array2DRowRealMatrix ArrayRealVector RealMatrix]))

(set! *warn-on-reflection* false)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defprotocol RegressionProto
  (backend [_])
  (model-native [_])
  (data-native [_])
  (predict [_ v] [_ v info?])
  (predict-all [_ vs] [_ v info?])
  (train [_] [_ x y])
  (stats [_])
  (cv [_] [_ params]))

(defmulti regression (fn [k & _] k))

(def ^:private regression-measure-objs
  (let [v {:rmse (RMSE.)
           :mad (MeanAbsoluteDeviation.)
           :mse (MSE.)
           :rss (RSS.)}]
    (assoc v :all (into-array RegressionMeasure (map v (keys v))))))
(def ^:private regression-measures [:rmse :mad :mse :rss])

(defn- ensure-vectors
  [xs]
  (if (sequential? (first xs)) xs (mapv vector xs)))

(defmethod regression :smile
  [_ ^RegressionTrainer trainer x y]
  (let [data-x (m/seq->double-double-array (ensure-vectors x))
        data-y (m/seq->double-array y)
        regression-raw (delay (.train trainer data-x data-y))]
    (reify
      IFn
      (invoke [r v] (predict r v))
      (invoke [r v _] (predict r v))

      RegressionProto
      (backend [_] :smile)
      (model-native [_] @regression-raw)
      (data-native [_] [data-x data-y])
      (predict [_ v] (.predict ^Regression @regression-raw (m/seq->double-array v)))
      (predict [r v _] (predict r v))
      (predict-all [_ vs] (seq (.predict ^Regression @regression-raw (m/seq->double-double-array vs))))
      (predict-all [r vs _] (predict-all r vs))
      (train [v] (do (deref regression-raw) v))
      (train [_ nx ny] (train (regression :smile trainer nx ny)))
      (cv [c] (cv c {}))
      (cv [_ {:keys [^int k type measure] :or {k 10 type :cv measure :rmse}}]
        (let [measure-obj (regression-measure-objs measure)
              result (case type
                       :cv (Validation/cv k trainer data-x data-y measure-obj)
                       :loocv (Validation/loocv trainer data-x data-y measure-obj)
                       :bootstrap (Validation/bootstrap k trainer data-x data-y measure-obj))]
          (if (= measure :all)
            (if-not (= type :bootstrap)
              (zipmap regression-measures result)
              (zipmap regression-measures (->> result
                                               (m/double-double-array->seq)
                                               (apply map vector)
                                               (map #(hash-map :mean (stats/mean %)
                                                               :stddev (stats/stddev %))))))
            (if-not (= type :bootstrap)
              {measure result}
              {measure {:mean (stats/mean result)
                        :stddev (stats/stddev result)}})))))))

(defmacro ^:private wrap-regression
  {:style/indent 3}
  [typ rname parameter instance & r]
  (let [[x y] (map symbol ["x" "y"])
        doc (str rname " regression. Backend library: " (name typ))]
    `(defn ~rname ~doc
       {:metadoc/categories #{:reg}}
       ([~x ~y] (~rname {} ~x ~y))
       ([~parameter ~x ~y] (regression ~typ ~instance ~x ~y ~@r)))))

(wrap-regression :smile ols {} (OLS$Trainer.))

(wrap-regression :smile rls {} (RLS$Trainer.))

(wrap-regression :smile lasso {:keys [lambda tolerance max-iters]
                               :or {lambda 10.0 tolerance 1.0e-3 max-iters 1000}}
                 (-> (LASSO$Trainer. lambda)
                     (.setTolerance tolerance)
                     (.setMaxNumIteration max-iters)))

(wrap-regression :smile elastic-net {:keys [lambda1 lambda2 tolerance max-iters]
                                     :or {lambda1 0.1 lambda2 0.1 tolerance 1.0e-4 max-iters 1000}}
                 (proxy [RegressionTrainer] []
                   (train [xx yy]
                     (ElasticNet. xx yy lambda1 lambda2 tolerance max-iters))))

(wrap-regression :smile ridge {:keys [lambda]
                               :or {lambda 0.1}}
                 (RidgeRegression$Trainer. lambda))

(wrap-regression :smile rbf-network {:keys [distance rbf number-of-basis normalize?]
                                     :or {distance dist/euclidean number-of-basis 10 normalize? false}}
  (let [cl (RBFNetwork$Trainer. distance)]
    (-> (cond
          (nil? rbf) cl
          (sequential? rbf) (.setRBF cl (into-array smile.math.rbf.RadialBasisFunction (map k/smile-rbf rbf)))
          :else (.setRBF cl (k/smile-rbf rbf) number-of-basis))
        (.setNormalized normalize?))))

(wrap-regression :smile svr {:keys [kernel C eps tolerance]
                             :or {kernel (k/kernel :linear) C 1.0 eps 0.001 tolerance 1.0e-3}}
  (-> (SVR$Trainer. (k/smile-mercer kernel) eps C)
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

(def ^:private loss-strategies {:least-squares GradientTreeBoost$Loss/LeastSquares
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

(wrap-regression :smile gaussian-process {:keys [kernel lambda]
                                          :or {kernel (k/kernel :gaussian) lambda 0.5}}
  (-> (GaussianProcessRegression$Trainer. (k/smile-mercer kernel) lambda)))

(def ^:private activation-functions {:logistic-sigmoid NeuralNetwork$ActivationFunction/LOGISTIC_SIGMOID
                                     :tanh NeuralNetwork$ActivationFunction/TANH})

(wrap-regression :smile neural-net {:keys [activation-function layers learning-rate momentum weight-decay number-of-epochs]
                                    :or {error-function :logistic-sigmoid learning-rate 0.1 momentum 0.0 weight-decay 0.0 number-of-epochs 25}}
  (let [fl (if (sequential? (first x)) (count (first x)) 1) ;; first layer - input
        mid (if (seq layers) (vec layers) [(inc fl)]) ;; mid layer, if empty, insert artificial
        layers (into-array Integer/TYPE (cons fl (conj mid 1)))]
    (-> (if-not activation-function
          (NeuralNetwork$Trainer. layers)
          (NeuralNetwork$Trainer. (or (activation-functions activation-function)
                                      NeuralNetwork$ActivationFunction/LOGISTIC_SIGMOID) layers))
        (.setLearningRate learning-rate)
        (.setMomentum momentum)
        (.setWeightDecay weight-decay)
        (.setNumEpochs number-of-epochs))))

;; native gaussian processes
;; Based on: https://www.cs.ubc.ca/~nando/540-2013/lectures/gp.py

(defn- kernel-cov-matrix
  (^Array2DRowRealMatrix [kernel ^double scale xss xss*]
   (MatrixUtils/createRealMatrix
    (m/seq->double-double-array (for [x xss]
                                  (map #(* scale ^double (kernel x %)) xss*)))))
  (^Array2DRowRealMatrix [kernel xss xss*] (kernel-cov-matrix kernel 1.0 xss xss*)))

(defprotocol GPProto
  (prior-samples [_ vs] "Draw samples from prior for given vs")
  (posterior-samples [gp vs] [gp vs stddev?] "Draw samples from posterior for given vs"))

(defn gaussian-process+
  ([x y] (gaussian-process+ {} x y))
  ([{:keys [^double kscale kernel noise normalize?]
     :or {kscale 1.0 kernel (k/kernel :gaussian 1.0) normalize? false}} xs y]

   (let [xs (ensure-vectors xs)
         ymean (if normalize? (stats/mean y) 0.0)
         ^ArrayRealVector ys (MatrixUtils/createRealVector (m/seq->double-array (map #(- ^double % ymean) y)))
         ^RealMatrix data-cov (kernel-cov-matrix kernel kscale xs xs)
         noise-fn #(if (sequential? noise)
                     (take (count %) (cycle noise))
                     (repeat (count %) (or noise 1.0e-6)))
         ^RealMatrix diag (MatrixUtils/createRealDiagonalMatrix (m/seq->double-array (noise-fn xs)))
         ^CholeskyDecomposition chol (CholeskyDecomposition. (.add data-cov diag))
         L (.getL chol)
         alpha (.solve (.getSolver chol) ys)]
     
     (MatrixUtils/solveLowerTriangularSystem L ys)
     
     (reify
       IFn
       (invoke [gp v] (predict gp v))
       (invoke [gp v stddev?] (predict gp v stddev?))

       RegressionProto
       (backend [_] :fastmath)
       (predict [gp xval] (predict gp xval false))
       (predict [_ xval stddev?]
         (let [xtest (if (sequential? xval) xval [xval])
               cov-vector (double-array (map #(* kscale ^double (kernel xtest %)) xs))
               mu (+ ymean ^double (v/dot cov-vector (.getDataRef ^ArrayRealVector alpha)))]
           (if-not stddev?
             mu
             (let [cov-v (MatrixUtils/createRealVector cov-vector)]
               (MatrixUtils/solveLowerTriangularSystem L cov-v)
               [mu (m/safe-sqrt (- 1.0 (.dotProduct cov-v cov-v)))]))))
       (predict-all [gp xtest] (predict-all gp xtest false))
       (predict-all [gp xtest stddev?]
         (map #(predict gp % stddev?) xtest))
       (train [gp] gp)
       (train [_ x y] (gaussian-process+ {:kscale kscale :kernel kernel :normalize? normalize?} x y))


       GPProto
       (prior-samples [_ xs]
         (let [xs (ensure-vectors xs)
               ^RealMatrix cov (kernel-cov-matrix kernel kscale xs xs)
               ^RealMatrix diag (MatrixUtils/createRealDiagonalMatrix (m/seq->double-array (noise-fn xs)))
               ^RealMatrix Lp (.getL ^CholeskyDecomposition (CholeskyDecomposition. (.add cov diag)))]
           (seq (.operate Lp (m/seq->double-array (repeatedly (count xs) r/grand))))))

       (posterior-samples [gp xtest] (posterior-samples gp xtest false))
       (posterior-samples [_ xtest stddev?]
         (let [xtest (ensure-vectors xtest)
               ^RealMatrix cov (kernel-cov-matrix kernel kscale xs xtest)
               cov-v (mapv #(.getColumnVector cov %) (range (.getColumnDimension cov)))]
           (run! #(MatrixUtils/solveLowerTriangularSystem L %) cov-v)
           (let [Lk (MatrixUtils/createRealMatrix (m/seq->double-double-array (map #(.getDataRef ^ArrayRealVector %) cov-v)))
                 diag (MatrixUtils/createRealDiagonalMatrix (m/seq->double-array (noise-fn xtest)))
                 ^RealMatrix k2 (kernel-cov-matrix kernel kscale xtest xtest)
                 ^RealMatrix Lp (.getL ^CholeskyDecomposition (CholeskyDecomposition. (.add k2 (.subtract ^RealMatrix diag (.multiply ^RealMatrix Lk (.transpose Lk)))))) ;; opposite than in source
                 mu (map (fn [v ^double n]
                           (+ ymean (.dotProduct ys v) n)) cov-v (seq (.operate Lp (m/seq->double-array (repeatedly (count xtest) r/grand)))))]
             (if-not stddev?
               mu
               (let [stddev (map (fn [^long id ^double v] (m/sqrt (- ^double (.getEntry k2 id id) v)))
                                 (range (count xtest))
                                 (map #(.dotProduct ^ArrayRealVector % ^ArrayRealVector %) cov-v))]
                 (map vector mu stddev))))))))))

;; validation

(defn validate
  "Validate data against trained regression."
  [model tx ty]
  (let [pred (predict-all model tx)
        pred-arr (m/seq->double-array pred)
        truth (m/seq->double-array ty)]
    {:prediction pred
     :stats (zipmap regression-measures
                    (map #(let [^RegressionMeasure m (regression-measure-objs %)]
                            (.measure m truth pred-arr)) regression-measures))}))
