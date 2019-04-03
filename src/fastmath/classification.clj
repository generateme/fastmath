(ns fastmath.classification
  "Classification algorithms.

  ### Input data

  * features - sequence of sequences of numbers
  * categories - sequence of any values

  When needed you can split input data into training data and test data. Test data are optional.
  
  ### Workflow
  
  * create and train by calling classifier with parameters, train data and optional test data.
  * validate data by calling [[test]] or [[validate]].
  * cross validate [[cv]], [[loocv]], [[bootstrap]]
  * repeat or [[predict]]

  Classifier parameters are map of values specific for given algorithm. Check documentation in backend library to find documentation. Classifier can be retrained using [[train]]. New instance will be created.
  
  ### Implementation notes

  * only doubles as input data
  * categories can be any type

  ### Cross validation

  Every classifier exposes it's own cross validation method with configuratin [[cv-native]]. Additionally three Clojure level methods are defined: [[cv]], [[loocv]] and [[bootstrap]].
  
  ### SMILE

  [documentation](https://haifengl.github.io/smile/classification.html)
  
  What is missed:

  * other types than doubles (attributes)
  * maxent
  * Online classifiers
  * General Naive Bayes

  Native cross validation config is a map with keys:

  * `:k` - number of folds (default: 10)
  * `:type` - type of cross validation, one of `:cv` (default), `:loocv` and `:bootstrap`
  
  ### XGBoost

  **TURNED OFF**
  
  XGBoost is backed by [clj-boost](https://gitlab.com/alanmarazzi/clj-boost/tree/master).

  For:

  * multiclass classification use `multi:softprob` objective (default).
  * binary classification use `binary-logistic` objective.

  See [documentation](https://xgboost.readthedocs.io/en/latest/parameter.html#parameters-for-tree-booster) for parameters.

  By default `:early-stopping` is set to `20` and `:round` is set to `10`.

  Native cross validation is a map with keys

  * `:nfold` - number of folds (default: 10)
  * `:rounds` - number of boosting iterations (default: 3)
  * `:metrics` - metrics to evaluate goodness (default: `[\"merror\", \"mlogloss\"]` or `[\"error\", \"logloss\"]`)

  Note: `\"(m)error\" = 1.0 - accuracy`
  
  see [more](https://gitlab.com/alanmarazzi/clj-boost/tree/master)
  
  ### Liblinear

  https://github.com/bwaldvogel/liblinear-java

  Native cross validation expects a number as number of folds (default: 10)
  
  ### Examples

  Iris database is used."
  {:metadoc/categories {:cl "Classification"
                        :vd "Validation"
                        :dt "Data operation"}}
  (:refer-clojure :exclude [test])
  (:require [fastmath.core :as m]
            [fastmath.distance :as dist]
            [fastmath.rbf :as rbf]
            ;; [clj-boost.core :as xgboost]
            [fastmath.stats :as stat]
            [fastmath.kernel.mercer :as mercer])
  (:import [clojure.lang IFn]
           [smile.classification Classifier SoftClassifier ClassifierTrainer KNN$Trainer AdaBoost$Trainer FLD$Trainer QDA$Trainer
            LDA$Trainer DecisionTree$Trainer DecisionTree$SplitRule GradientTreeBoost$Trainer
            LogisticRegression$Trainer Maxent$Trainer NaiveBayes$Trainer NaiveBayes$Model
            NeuralNetwork$Trainer NeuralNetwork$ErrorFunction NeuralNetwork$ActivationFunction
            RBFNetwork$Trainer RDA$Trainer RandomForest$Trainer SVM$Trainer SVM$Multiclass]
           [smile.math.distance EuclideanDistance]
           [smile.math.rbf RadialBasisFunction]
           [smile.math.kernel MercerKernel]
           [smile.validation Bootstrap Validation]

           [de.bwaldvogel.liblinear Problem Feature FeatureNode Linear Model Parameter SolverType]))

(set! *warn-on-reflection* false)

(defprotocol ClassificationProto
  (backend [_] "Return name of backend library")
  (^{:metadoc/categories #{:cl}} model-native [_] "Return trained model as a backend class.")
  (^{:metadoc/categories #{:cl}} predict [_ v] [_ v posteriori?] "Predict class for given vector. With posteriori probabilities when required.")
  (^{:metadoc/categories #{:cl}} predict-all [_ vs] [_ v posteriori?] "Predict classes for given sequence of vectors. With posteriori probabilities when required.")
  (^{:metadoc/categories #{:cl}} train [_ x y] [_ x y tx ty] "Train another set of data for given classifier. Test data are optional.")
  (^{:metadoc/categories #{:dt}} data [_] [_ native?] "Return data. Transformed data for backend libarary are returned when `native?` is true.")
  (^{:metadoc/categories #{:dt}} labels [_]  "Return labels.")
  (^{:metadoc/categories #{:vd}} test [_] [_ tx ty] "Validate data using given model and provided test data. Same as [[validate]].")
  (^{:metadoc/categories #{:vd}} cv-native [_] [_ params] "Call native implementation of cross-validation"))

(defn- labels-converters
  "Convert y into label->int and int->label functions."
  [y]
  (let [sort-fn (if (instance? Comparable (first y)) sort identity)
        ydata (mapv vector (sort-fn (distinct y)) (range))]
    [(mapv first ydata) (into {} ydata)]))

(defn- liblinear-to-features
  [fs]
  (into-array Feature (map-indexed (fn [id v] (FeatureNode. (inc id) v)) fs)))

(defmulti ^:private prepare-data (fn [k & _] k))
(defmethod prepare-data :smile [_ x y labels->int] [(m/seq->double-double-array x) (int-array (mapv labels->int y))])
#_(defmethod prepare-data :xgboost [_ x y labels->int] (xgboost/dmatrix x (mapv labels->int y)))
(defmethod prepare-data :liblinear [_ x y bias labels->int]
  (let [x (if-not (pos? bias) x
                  (map #(conj (vec %) bias) x))
        fa (into-array (map liblinear-to-features x))
        problem (Problem.)]
    (set! (.-x problem) fa)
    (set! (.-y problem) (double-array (map labels->int y)))
    (set! (.-n problem) (count (first x)))
    (set! (.-l problem) (count y))
    (set! (.-bias problem) bias)
    problem))

(declare validate)
(declare accuracy)

(defmulti ^:private classifier (fn [k & _] k))

(comment defrecord LibLinearClassifier [x y tx ty bias model labels labels->int data params backend]
         IFn
         (invoke [_ v] (labels (int (Linear/predict model (liblinear-to-features v)))))
         (invoke [c v posteriori?] (let [buff (double-array (count labels))]
                                     (if posteriori?
                                       [(labels (int (Linear/predictProbability model (liblinear-to-features v) buff))) (vec buff)]
                                       (c v))))
         ClassificationProto
         (predict [c v] (c v))
         (predict [c v posteriori?] (c v posteriori?))
         (predict-all [_ vs] (map (comp labels int #(Linear/predict model (liblinear-to-features %))) vs))
         (predict-all [c vs posteriori?] (if posteriori?
                                           (map #(c % true) vs)
                                           (predict-all c vs)))
         (cv-native [c] (cv-native c 10))
         (cv-native [_ k]
                    (let [target (double-array (count x))]
                      (Linear/crossValidation data params (or k 10) target)
                      {:accuracy (accuracy y (map (comp labels int) target))}))
         (train [c x y] (train c x y nil nil))
         (train [_ x y tx ty] (classifier :liblinear params x y tx ty bias)))

(defmethod classifier :liblinear
  [_ params x y tx ty label-map bias]
  (let [[labels labels->int] (labels-converters y)
        data (prepare-data :liblinear x y (or bias -1) labels->int)
        ^Model model (delay (Linear/train data params))]
    ;; (->LibLinearClassifier x y tx ty bias model labels labels->int data params :liblinear)
    (reify
      IFn
      (invoke [_ v] (labels (int (Linear/predict @model (liblinear-to-features v)))))
      (invoke [c v posteriori?] (let [buff (double-array (count labels))]
                                  (if posteriori?
                                    [(labels (int (Linear/predictProbability @model (liblinear-to-features v) buff))) (vec buff)]
                                    (c v))))

      ClassificationProto
      (backend [_] :liblinear)
      (model-native [_] @model)

      (predict [c v] (c v))
      (predict [c v posteriori?] (c v posteriori?))
      (predict-all [_ vs] (map (comp labels int #(Linear/predict @model (liblinear-to-features %))) vs))
      (predict-all [c vs posteriori?] (if posteriori?
                                        (map #(c % true) vs)
                                        (predict-all c vs)))

      (cv-native [c] (cv-native c 10))
      (cv-native [_ k]
        (let [target (double-array (count x))]
          (Linear/crossValidation data params (or k 10) target)
          {:accuracy (accuracy y (map (comp labels int) target))}))

      (train [c x y] (train c x y nil nil))
      (train [_ x y tx ty] (classifier :liblinear params x y tx ty label-map bias))
      (test [c] (when (and tx ty) (validate c tx ty)))
      (test [c tx ty] (validate c tx ty))
      
      (data [_] [x y])
      (data [_ native?] (if native? data [x y]))
      (labels [_] labels))))

#_(defmethod classifier :xgboost
    ([_ booster-params x y tx ty label-map] (classifier :xgboost booster-params x y tx ty label-map nil))
    ([_ booster-params x y tx ty label-map booster]
     (let [[labels labels->int] (labels-converters y)
           data (prepare-data :xgboost x y labels->int)
           test-data (when tx (prepare-data :xgboost tx ty labels->int))
           binary? (= "binary:logistic" (get-in booster-params [:params :objective]))

           ;; to clean-up
           booster-params (if-not (contains? booster-params :rounds) (assoc booster-params :rounds 10) booster-params)
           booster-params (if-not (contains? booster-params :early-stopping) (assoc booster-params :early-stopping 20) booster-params)
           booster-params (if-not binary?
                            (-> booster-params
                                (assoc-in [:params :num_class] (count labels))
                                (assoc-in [:params :objective] "multi:softprob"))
                            booster-params)
           
           booster-params (-> (dissoc booster-params :booster :watches)
                              (assoc-in [:watches :train] data))
           
           booster-params (if (and tx ty) (assoc-in booster-params [:watches :test] test-data) booster-params)
           booster-params (if booster (assoc booster-params :booster booster) booster-params)

           model (delay (xgboost/fit data booster-params))
           res-range (range (count labels))
           prob->label (if binary?
                         #(if (< (first %) 0.5) (labels 0) (labels 1))
                         #(labels (apply max-key (vec %) res-range)))
           prob-map (if binary?
                      #(vector (- 1.0 %) %)
                      identity)]
       (reify
         IFn
         (invoke [_ v]  (prob->label (xgboost/predict @model (xgboost/dmatrix [v]))))
         (invoke [c v posteriori?] (if posteriori? (let [r (-> (xgboost/predict @model (xgboost/dmatrix [v])))]
                                                     [(prob->label r) (map prob-map r)]) (c v)))

         ClassificationProto
         (backend [_] :xgboost)
         (model-native [_] @model)
         
         (predict [c v] (c v))
         (predict [c v posteriori?] (c v posteriori?))
         (predict-all [_ vs] (map prob->label (seq (.predict ^ml.dmlc.xgboost4j.java.Booster @model (xgboost/dmatrix vs)))))
         (predict-all [c vs posteriori?] (if posteriori?
                                           (mapv (comp #(vector (prob->label %) (map prob-map %)) seq)
                                                 (seq (.predict ^ml.dmlc.xgboost4j.java.Booster @model (xgboost/dmatrix vs))))
                                           (predict-all c vs)))
         (cv-native [c] (cv-native c {}))
         (cv-native [_ {:keys [nfold rounds metrics]
                        :or {nfold 10 rounds 3 metrics (if binary? ["error", "logloss"] ["merror", "mlogloss"])}}]
           (xgboost/cross-validation data {:nfold nfold :rounds rounds :metrics metrics :params (booster-params :params)}))
         
         (train [c x y] (train c x y nil nil))
         (train [_ x y tx ty] (classifier :xgboost booster-params x y tx ty label-map @model))
         (test [c] (when (and tx ty) (validate c tx ty)))
         (test [c tx ty] (validate c tx ty))
         
         (data [_] [x y])
         (data [_ native?] (if native? data [x y]))
         (labels [_] labels)))))

(defmethod classifier :smile
  [_ ^ClassifierTrainer trainer x y tx ty label-map]
  (let [[labels labels->int] (labels-converters y)
        [data int-labels] (prepare-data :smile x y labels->int)
        classifier-raw (delay (.train trainer data int-labels))
        predict-raw #(.predict ^Classifier @classifier-raw (m/seq->double-array %))
        predict-fn (comp labels predict-raw)
        predict-fn-posteriori (if (instance? SoftClassifier classifier-raw)
                                #(let [posteriori (double-array (count labels))]
                                   [(labels (.predict ^SoftClassifier @classifier-raw (m/seq->double-array %) posteriori)) (seq posteriori)])
                                predict-fn)]
    (reify
      
      IFn
      (invoke [_ v] (predict-fn v))
      (invoke [_ v posteriori?] (if posteriori? (predict-fn-posteriori v) (predict-fn v)))

      ClassificationProto
      (backend [_] :smile)
      (model-native [_] @classifier-raw)
      
      (predict [_ v] (predict-fn v))
      (predict [_ v posteriori?] (if posteriori? (predict-fn-posteriori v) (predict-fn v)))
      (predict-all [_ vs] (map predict-fn vs))
      (predict-all [_ vs posteriori?] (if posteriori? (map predict-fn-posteriori vs) (map predict-fn vs)))
      (cv-native [c] (cv-native c {}))
      (cv-native [_ {:keys [^int k type] :or {k 10 type :cv}}]
        (case type
          :cv {:accuracy (Validation/cv k trainer data int-labels)}
          :loocv {:accuracy (Validation/loocv trainer data int-labels)}
          :bootstrap (let [b (Validation/bootstrap k trainer data int-labels)]
                       {:accuracy (stat/mean b) :bootstrap (vec b)})))
      
      (train [c x y] (train c x y nil nil))
      (train [_ x y tx ty] (classifier :smile trainer x y tx ty label-map))
      (test [c] (when (and tx ty) (validate c tx ty)))
      (test [c tx ty] (validate c tx ty))
      
      (data [_] [x y])
      (data [_ native?] (if native? [data int-labels] [x y]))
      (labels [_] labels))))

(defmacro ^:private wrap-classifier
  {:style/indent 3}
  [typ clname parameter instance & r]
  (let [[x y tx ty params all] (map symbol ["x" "y" "tx" "ty" "params" "all"])
        doc (str clname " classifier. Backend library: " (name typ))
        parameter (if (map? parameter) (assoc parameter :as all) parameter)
        label-map-getter (if (map? parameter) `(:label-map ~all) `(:label-map ~parameter))]
    `(defn ~clname ~doc
       {:metadoc/categories #{:cl}}
       ([~x ~y] (~clname {} ~x ~y nil nil))
       ([~params ~x ~y] (~clname ~params ~x ~y nil nil))
       ([~x ~y ~tx ~ty] (~clname {} ~x ~y ~tx ~ty))
       ([~parameter ~x ~y ~tx ~ty] (classifier ~typ ~instance ~x ~y ~tx ~ty (or ~label-map-getter identity) ~@r)))))

(wrap-classifier :smile knn {:keys [distance k]
                             :or {distance (EuclideanDistance.) k 1}}
  (KNN$Trainer. distance k))

(wrap-classifier :smile ada-boost {:keys [number-of-trees max-nodes]
                                   :or {number-of-trees 500 max-nodes 2}}
  (-> (AdaBoost$Trainer.)
      (.setNumTrees number-of-trees)
      (.setMaxNodes max-nodes)))

(wrap-classifier :smile fld {:keys [dimensionality tolerance]
                             :or {dimensionality -1 tolerance 1.0e-4}}
  (let [^FLD$Trainer t (FLD$Trainer.)]
    (-> (if-not (pos? dimensionality) t (.setDimension t dimensionality))
        (.setTolerance tolerance))))

(wrap-classifier :smile qda {:keys [priori tolerance]
                             :or {priori nil tolerance 1.0e-4}}
  (-> (QDA$Trainer.)
      (.setPriori (m/seq->double-array priori))
      (.setTolerance tolerance)))

(wrap-classifier :smile lda {:keys [priori tolerance]
                             :or {priori nil tolerance 1.0e-4}}
  (-> (LDA$Trainer.)
      (.setPriori (m/seq->double-array priori))
      (.setTolerance tolerance)))

(def split-rules {:gini DecisionTree$SplitRule/GINI
                  :entropy DecisionTree$SplitRule/ENTROPY
                  :classification-error DecisionTree$SplitRule/CLASSIFICATION_ERROR})

(def ^{:doc "List of split rules for [[decision tree]] and [[random-forest]]"} split-rules-list (keys split-rules))

(wrap-classifier :smile decision-tree {:keys [max-nodes node-size split-rule]
                                       :or {max-nodes 100 node-size 1 split-rule :gini}}
  (let [t (-> (DecisionTree$Trainer. max-nodes)
              (.setSplitRule (or (split-rules split-rule) DecisionTree$SplitRule/GINI)))]
    (if (> node-size 1) (.setNodeSize t node-size) t)))

(wrap-classifier :smile gradient-tree-boost {:keys [number-of-trees shrinkage max-nodes subsample]
                                             :or {number-of-trees 500 shrinkage 0.005 max-nodes 6 subsample 0.7}}
  (-> (GradientTreeBoost$Trainer. number-of-trees)
      (.setMaxNodes max-nodes)
      (.setShrinkage shrinkage)
      (.setSamplingRates subsample)))

(wrap-classifier :smile logistic-regression {:keys [lambda tolerance max-iterations]
                                             :or {lambda 0.0 tolerance 1.0e-5 max-iterations 500}}
  (-> (LogisticRegression$Trainer.)
      (.setRegularizationFactor lambda)
      (.setTolerance tolerance)
      (.setMaxNumIteration max-iterations)))

(def ^:private bayes-models {:multinomial NaiveBayes$Model/MULTINOMIAL
                             :bernoulli NaiveBayes$Model/BERNOULLI
                             :polyaurn NaiveBayes$Model/POLYAURN})

(def ^{:doc "List of [[naive-bayes]] models."} bayes-models-list (keys bayes-models))

(wrap-classifier :smile naive-bayes {:keys [model priori sigma]
                                     :or {model :bernoulli sigma 1.0}}
  (let [classes (count (distinct y))
        independent-variables (count (first x))
        model (or (bayes-models model) NaiveBayes$Model/BERNOULLI)]
    (-> (NaiveBayes$Trainer. model classes independent-variables)
        (.setSmooth sigma)
        (.setPriori (m/seq->double-array priori)))))

(def ^:private error-functions {:least-mean-squares NeuralNetwork$ErrorFunction/LEAST_MEAN_SQUARES
                                :cross-entropy NeuralNetwork$ErrorFunction/CROSS_ENTROPY})

(def ^{:doc "List of error functions for [[neural-net]]."} error-functions-list (keys error-functions))

(def ^:private activation-functions {:linear NeuralNetwork$ActivationFunction/LINEAR
                                     :logistic-sigmoid NeuralNetwork$ActivationFunction/LOGISTIC_SIGMOID
                                     :soft-max NeuralNetwork$ActivationFunction/SOFTMAX})

(def ^{:doc "List of activation functions for [[neural-net]]."} activation-functions-list (keys activation-functions))

(wrap-classifier :smile neural-net
    {:keys [error-function activation-function layers learning-rate momentum weight-decay number-of-epochs]
     :or {error-function :cross-entropy learning-rate 0.1 momentum 0.0 weight-decay 0.0 number-of-epochs 25}}
  (let [layers (into-array Integer/TYPE (cons (count (first x)) (conj (vec layers) (count (distinct y)))))
        ef (or (error-functions error-function) NeuralNetwork$ErrorFunction/CROSS_ENTROPY)]
    (-> (if-not activation-function
          (NeuralNetwork$Trainer. ef layers)
          (NeuralNetwork$Trainer.
           ef (or (activation-functions activation-function)
                  (if (= ef NeuralNetwork$ErrorFunction/CROSS_ENTROPY)
                    NeuralNetwork$ActivationFunction/SOFTMAX
                    NeuralNetwork$ActivationFunction/LINEAR)) layers))
        (.setLearningRate learning-rate)
        (.setMomentum momentum)
        (.setWeightDecay weight-decay)
        (.setNumEpochs number-of-epochs))))

(wrap-classifier :smile rbf-network {:keys [distance rbf number-of-basis normalize?]
                                     :or {distance dist/euclidean number-of-basis 10 normalize? false}}
  (let [cl (RBFNetwork$Trainer. distance)]
    (-> (cond
          (nil? rbf) cl
          (sequential? rbf) (.setRBF cl (into-array RadialBasisFunction rbf))
          :else (.setRBF cl rbf number-of-basis))
        (.setNormalized normalize?))))

(wrap-classifier :smile rda {:keys [alpha priori tolerance]
                             :or {alpha 0.9 tolerance 1.0e-4}}
  (-> (RDA$Trainer. alpha)
      (.setPriori (m/seq->double-array priori))
      (.setTolerance tolerance)))

(wrap-classifier :smile random-forest {:keys [number-of-trees split-rule mtry node-size max-nodes subsample]
                                       :or {number-of-trees 500 split-rule :gini node-size 1 max-nodes 100 subsample 1.0}}
                 (let [mtry (or mtry (m/floor (m/sqrt (count (first x)))))]
                   (-> (RandomForest$Trainer. number-of-trees)
                       (.setSplitRule (or (split-rules split-rule) DecisionTree$SplitRule/GINI))
                       (.setNumRandomFeatures mtry)
                       (.setMaxNodes max-nodes)
                       (.setSamplingRates subsample)
                       (.setNodeSize node-size))))

#_(wrap-classifier :xgboost xgboost xgboost-params xgboost-params)

(def ^:private multiclass-strategies {:one-vs-one SVM$Multiclass/ONE_VS_ONE
                                      :one-vs-all SVM$Multiclass/ONE_VS_ALL})

(def ^{:doc "List of multiclass strategies for [[svm]]"} multiclass-strategies-list (keys multiclass-strategies))

(wrap-classifier :smile svm {:keys [^MercerKernel kernel ^double c-or-cp ^double cn strategy-for-multiclass class-weights tolerance ^int epochs]
                             :or {kernel (mercer/kernel :linear) c-or-cp 1.0 cn 1.0 strategy-for-multiclass :one-vs-one tolerance 1.0e-3 epochs 2}}
  (let [cn (or cn c-or-cp)
        classes-no (count (distinct y))
        ^SVM$Multiclass ms (or (multiclass-strategies strategy-for-multiclass) SVM$Multiclass/ONE_VS_ONE)
        t (if (== 2 classes-no)
            (SVM$Trainer. kernel c-or-cp cn)
            (if class-weights
              (SVM$Trainer. kernel c-or-cp (double-array class-weights) ms)
              (SVM$Trainer. kernel c-or-cp classes-no ms)))]
    (-> (.setTolerance t tolerance)
        (.setNumEpochs epochs))))

(def ^:private liblinear-solvers {:l2r-lr SolverType/L2R_LR
                                  :l2r-l2loss-svc-dual SolverType/L2R_L2LOSS_SVC_DUAL
                                  :l2r-l2loss-svc SolverType/L2R_L2LOSS_SVC
                                  :l2r-l1loss-svc-dual SolverType/L2R_L1LOSS_SVC_DUAL
                                  :mcsvm-cs SolverType/MCSVM_CS
                                  :l1r-l2loss-svc SolverType/L1R_L2LOSS_SVC
                                  :l1r-lr SolverType/L1R_LR
                                  :l2r-lr-dual SolverType/L2R_LR_DUAL})

(def ^{:doc "List of [[liblinear]] solvers."} liblinear-solver-list (keys liblinear-solvers))

(wrap-classifier :liblinear liblinear {:keys [solver bias ^double C ^double eps ^int max-iters ^double p weights]
                                       :or {solver :l2r-l2loss-svc-dual bias -1 C 1.0 eps 0.01 max-iters 1000 p 0.1}}
                 (let [par (Parameter. (or (liblinear-solvers solver) SolverType/L2R_LR) C eps max-iters p)]
                   (if weights
                     (do (.setWeight par (double-array weights) (int-array (range (count weights)))) par)
                     par)) bias)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; validation metrics

(defn accuracy
  "Calculate accuracy for real and predicted sequences."
  {:metadoc/categories #{:vd}}
  [t p] (double (/ (count (filter (partial apply =) (mapv vector t p)))
                   (count t))))

(defn validate
  "Validate data against trained classifier. Same as [[test]]."
  {:metadoc/categories #{:vd}}
  ([model] (test model))
  ([model tx ty]
   (let [pred (predict-all model tx)
         invalid (->> (mapv vector tx ty pred)
                      (filter #(apply not= (rest %))))
         invalid-cnt (count invalid)]
     {:truth ty
      :prediction pred
      :invalid {:count invalid-cnt
                :data (map first invalid)
                :prediction (map #(nth % 2) invalid)
                :truth (map second invalid)}
      :stats {:accuracy (double (- 1.0 (/ invalid-cnt (count ty))))}})))

(defn- drop-nth [coll n]
  (keep-indexed #(if (not= %1 n) %2) coll))

(defn- process-chunks
  [model predict-fn data labels mapcat? id]
  (let [tdata (if mapcat? (mapcat identity (drop-nth data id)) (drop-nth data id))
        tlabels (if mapcat? (mapcat identity (drop-nth labels id)) (drop-nth labels id))
        value (nth data id)
        mtrain (train model tdata tlabels)]
    (predict-fn mtrain value)))

(def ^:dynamic ^{:doc "When `true` provide data used to test."}
  *cross-validation-debug* false)

(defn loocv
  "Leave-one-out cross validation of a classification model."
  {:metadoc/categories #{:vd}}
  ([model]
   (let [[data labels] (data model)
         plabels (mapv (partial process-chunks model predict data labels false) (range (count data)))
         stats {:accuracy (accuracy labels plabels)}]
     (if *cross-validation-debug*
       (merge stats {:data data :truth labels :prediction plabels})
       stats))))

(defn- slice [data ids] (mapv (partial nth data) ids))

(defn cv
  "Cross validation of a classification model.

  k defaults to 10% of data count."
  {:metadoc/categories #{:vd}}
  ([model]
   (cv model  (* 0.1 (count (first (data model))))))
  ([model k]
   (let [[data labels] (data model)
         dsize (count data)
         chunksize (max 2 (min (* 0.75 dsize) (/ dsize k)))
         ids (shuffle (range dsize))
         sdata (slice data ids)
         slabels (slice labels ids)
         psdata (partition-all chunksize sdata)
         plabels (mapcat (partial process-chunks model predict-all
                                  psdata (partition-all chunksize slabels) true) (range (count psdata)))
         stats {:accuracy (accuracy slabels plabels)}]
     (if *cross-validation-debug*
       (merge stats {:data sdata :truth slabels :prediction plabels})
       stats))))

(defn bootstrap
  "Perform k-round bootstrap validation.

  k defaults to 10

  Returns map where `:accuracy` is average accuracy from every round. `:boostrap` contains every round statistics."
  {:metadoc/categories #{:vd}}
  ([model] (bootstrap model 10))
  ([model k]
   (let [[data labels] (data model)
         ^Bootstrap bootstrap (Bootstrap. (count data) k)
         all (for [[train-ids test-ids] (map vector (.-train bootstrap) (.-test bootstrap))
                   :let [sdata (slice data train-ids)
                         slabels (slice labels train-ids)
                         test-labels-true (slice labels test-ids)
                         test-data (slice data test-ids)
                         m (train model sdata slabels)
                         test-labels-pred (predict-all m test-data)
                         stats {:accuracy (accuracy test-labels-true test-labels-pred)}]]
               (if *cross-validation-debug*
                 (merge stats {:data test-data :truth test-labels-true :prediction test-labels-pred})
                 stats))
         avg (stat/mean (map :accuracy all))]
     {:accuracy avg :bootstrap (into (sorted-map) (map-indexed vector all))})))

(defn confusion-map
  "Create confusion map where keys are pairs of `[truth-label, prediction-label]`"
  {:metadoc/categories #{:vd}}
  [t p]
  (reduce (fn [m tpv]
            (if (contains? m tpv)
              (update m tpv inc)
              (assoc m tpv 1))) {} (map vector t p)))


;; binary classification statistics

(defn- binary-confusion
  [t p]
  (cond
    (and t p) :tp
    (and t (not p)) :fn
    (and (not t) p) :fp
    :else :tn))

(defn- binary-confusion-val
  [t p tv fv]
  (cond
    (and (= t tv) (= p tv)) :tp
    (and (= t tv) (= p fv)) :fn
    (and (= t fv) (= p tv)) :fp
    :else :tn))

(defn binary-measures-all
  "https://en.wikipedia.org/wiki/Precision_and_recall"
  [truth prediction]
  (let [{:keys [^double tp ^double fp ^double fn ^double tn] :as details} (frequencies (map binary-confusion truth prediction))
        cp (+ tp fn)
        cn (+ fp tn)
        total (+ cp cn)
        pcp (+ tp fp)
        pcn (+ fn tn)
        ppv (/ tp pcp)
        npv (/ tn pcn)
        tpr (/ tp cp)
        fpr (/ fp cn)
        tnr (- 1.0 fpr)
        fnr (- 1.0 tpr)
        lr+ (/ tpr fpr)
        lr- (/ fnr tnr)
        f-beta (clojure.core/fn [beta] (let [b2 (* beta beta)]
                                         (* (inc b2) (/ (* ppv tpr)
                                                        (+ ppv tpr)))))
        f1-score (f-beta 1.0)]
    (merge details {:cp cp
                    :cn cn
                    :pcp pcp
                    :pcn pcn
                    :total total
                    :tpr tpr
                    :recall tpr
                    :sensitivity tpr
                    :hit-rate tpr
                    :fnr fnr
                    :miss-rate fnr
                    :fpr fpr
                    :fall-out fpr
                    :tnr tnr
                    :specificity tnr
                    :selectivity tnr
                    :prevalence (/ cp total)
                    :accuracy (/ (+ tp tn) total)
                    :ppv ppv
                    :precision ppv
                    :fdr (- 1.0 ppv)
                    :npv npv
                    :for (- 1.0 npv)
                    :lr+ lr+
                    :lr- lr-
                    :dor (/ lr+ lr-)
                    :f-measure f1-score
                    :f1-score f1-score
                    :f-beta f-beta
                    :mcc (/ (- (* tp tn) (* fp fn))
                            (m/sqrt (* (+ tp fp)
                                       (+ tp fn)
                                       (+ tn fp)
                                       (+ tn fn))))
                    :bm (dec (+ tpr tnr))
                    :mk (dec (+ ppv npv))})))

(defn binary-measures
  [truth prediction]
  (select-keys (binary-measures-all truth prediction)
               [:tp :tn :fp :fn :accuracy :fdr :f-measure :fall-out :precision :recall :sensitivity :specificity :prevalance]))

