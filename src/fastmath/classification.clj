(ns fastmath.classification
  "Classification algorithms.

  ### Input data

  * features - sequence of sequences of numbers
  * categories - sequence of any values

  ### Workflow
  
  * create classifier with parameters
  * cross validate [[cv]]
  * repeat or [[predict]]
  * to validate model against test data, call [[validate]]

  Classifier parameters are map of values specific for given algorithm. Check documentation in backend library to find documentation. Classifier can be retrained using [[train]]. New instance will be created.

  Classifier training is delayed to the actual use. To force training, call [[train]] or [[predict]].
  
  ### Implementation notes

  * only doubles as input data
  * categories can be any type

  ### Cross validation

  Every classifier exposes it's own cross validation method with configuration [[cv]].

  ~~Additionally three Clojure level methods are defined: [[cv]], [[loocv]] and [[bootstrap]].~~
  
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
            [fastmath.kernel :as k]
            [fastmath.stats :as stats]
            [fastmath.protocols :as pr])
  (:import [clojure.lang IFn]
           [smile.classification Classifier SoftClassifier ClassifierTrainer KNN$Trainer AdaBoost$Trainer FLD$Trainer QDA$Trainer
            LDA$Trainer DecisionTree$Trainer DecisionTree$SplitRule GradientTreeBoost$Trainer
            LogisticRegression$Trainer Maxent$Trainer NaiveBayes$Trainer NaiveBayes$Model
            NeuralNetwork$Trainer NeuralNetwork$ErrorFunction NeuralNetwork$ActivationFunction
            RBFNetwork$Trainer RDA$Trainer RandomForest$Trainer SVM$Trainer SVM$Multiclass]
           [smile.math.distance EuclideanDistance]
           [smile.validation Bootstrap Validation]
           [de.bwaldvogel.liblinear Problem Feature FeatureNode Linear Model Parameter SolverType]))

(set! *warn-on-reflection* false)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn- labels-converters
  "Convert y into label->int and int->label functions."
  [y]
  (let [sort-fn (if (instance? Comparable (first y)) sort identity)
        ydata (mapv vector (sort-fn (distinct y)) (range))]
    [(mapv first ydata) (into {} ydata)]))

(defn- liblinear-to-features
  [fs]
  (into-array Feature (map-indexed (fn [^long id v] (FeatureNode. (inc id) v)) fs)))

(defn- ensure-vectors
  [xs]
  (if (sequential? (first xs)) xs (mapv vector xs)))

(defmulti ^:private prepare-data (fn [k & _] k))
(defmethod prepare-data :smile [_ x y labels->int] [(m/seq->double-double-array (ensure-vectors x))
                                                    (int-array (mapv labels->int y))])
(defmethod prepare-data :liblinear [_ x y bias labels->int]
  (let [x (ensure-vectors x)
        x (if-not (pos? (double bias)) x
                  (map #(conj (vec %) bias) x))
        fa (into-array (map liblinear-to-features x))
        problem (Problem.)]
    (set! (.-x problem) fa)
    (set! (.-y problem) (double-array (map labels->int y)))
    (set! (.-n problem) (count (first x)))
    (set! (.-l problem) (count y))
    (set! (.-bias problem) bias)
    problem))

(defmulti ^:private classifier (fn [k & _] k))

(defn accuracy
  "Calculate accuracy for real and predicted sequences."
  {:metadoc/categories #{:vd}}
  [t p] (/ (count (filter (partial apply =) (mapv vector t p)))
           (double (count t))))

(defmethod classifier :liblinear
  [_ params x y bias]
  (let [[labels labels->int] (labels-converters y)
        data (prepare-data :liblinear x y (or bias -1) labels->int)
        ^Model model (delay (Linear/train data params))]
    (reify
      IFn
      (invoke [_ v] (labels (int (Linear/predict @model (liblinear-to-features v)))))
      (invoke [c v posteriori?] (let [buff (double-array (count labels))]
                                  (if posteriori?
                                    [(labels (int (Linear/predictProbability @model (liblinear-to-features v) buff))) (vec buff)]
                                    (c v))))

      pr/PredictorProto
      (backend [_] :liblinear)
      (model-native [_] @model)
      (data-native [_] [data labels])

      (predict [c v posteriori?] (c v posteriori?))
      (predict-all [c vs posteriori?] (if posteriori?
                                        (map #(c % true) vs)
                                        (map (comp labels int #(Linear/predict @model (liblinear-to-features %))) vs)))

      (cv [c] (pr/cv c 10))
      (cv [_ k]
        (let [target (double-array (count x))]
          (Linear/crossValidation data params (or k 10) target)
          {:accuracy (accuracy y (map (comp labels int) target))}))

      (train [c] (do (deref model) c))
      (train [_ x y] (pr/train (classifier :liblinear params x y bias))))))

(defmethod classifier :smile
  [_ ^ClassifierTrainer trainer x y]
  (let [[labels labels->int] (labels-converters y)
        [data int-labels :as internal-data] (prepare-data :smile x y labels->int)
        classifier-raw (delay (.train trainer data int-labels))
        predict-raw #(.predict ^Classifier @classifier-raw (m/seq->double-array %))
        predict-fn (comp labels predict-raw)
        predict-fn-posteriori #(let [posteriori (double-array (count labels))]
                                 [(labels (.predict ^SoftClassifier @classifier-raw (m/seq->double-array %) posteriori))
                                  (zipmap labels posteriori)])]
    
    (reify
      
      IFn
      (invoke [_ v] (predict-fn v))
      (invoke [_ v posteriori?] (if posteriori? (predict-fn-posteriori v) (predict-fn v)))

      pr/PredictorProto
      (backend [_] :smile)
      (model-native [_] @classifier-raw)
      (data-native [_] internal-data)
      
      (predict [_ v posteriori?] (if posteriori? (predict-fn-posteriori v) (predict-fn v)))
      (predict-all [_ vs posteriori?] (if posteriori? (map predict-fn-posteriori vs) (map predict-fn vs)))
      (cv [c] (pr/cv c {}))
      (cv [_ {:keys [^int k type] :or {k 10 type :cv}}]
        (case type
          :cv {:accuracy (Validation/cv k trainer data int-labels)}
          :loocv {:accuracy (Validation/loocv trainer data int-labels)}
          :bootstrap (let [b (Validation/bootstrap k trainer data int-labels)]
                       {:accuracy {:mean (stats/mean b)
                                   :stddev (stats/stddev b)}})))
      
      (train [c] (do (deref classifier-raw) c))
      (train [_ x y] (pr/train (classifier :smile trainer x y))))))

(defmacro ^:private wrap-classifier
  {:style/indent 3}
  [typ clname parameter instance & r]
  (let [[x y] (map symbol ["x" "y"])
        doc (str clname " classifier. Backend library: " (name typ))]
    `(defn ~clname ~doc
       {:metadoc/categories #{:cl}}
       ([~x ~y] (~clname {} ~x ~y))
       ([~parameter ~x ~y] (classifier ~typ ~instance ~x ~y ~@r)))))

(wrap-classifier :smile knn {:keys [distance k]
                             :or {distance (EuclideanDistance.) k 1}}
  (KNN$Trainer. distance k))

(wrap-classifier :smile ada-boost {:keys [number-of-trees max-nodes]
                                   :or {number-of-trees 500 max-nodes 2}}
  (-> (AdaBoost$Trainer.)
      (.setNumTrees number-of-trees)
      (.setMaxNodes max-nodes)))

(wrap-classifier :smile fld {:keys [^long dimensionality tolerance]
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

(def ^:private split-rules {:gini DecisionTree$SplitRule/GINI
                            :entropy DecisionTree$SplitRule/ENTROPY
                            :classification-error DecisionTree$SplitRule/CLASSIFICATION_ERROR})

(def ^{:doc "List of split rules for [[decision tree]] and [[random-forest]]"} split-rules-list (keys split-rules))

(wrap-classifier :smile decision-tree {:keys [max-nodes ^long node-size split-rule]
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
        independent-variables (if (sequential? (first x)) (count (first x)) 1)
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
  (let [fl (int (if (sequential? (first x)) (count (first x)) 1)) ;; first layer - input
        mid (if (seq layers) (vec layers) [(inc fl)]) ;; mid layer, if empty, insert artificial
        layers (into-array Integer/TYPE (cons fl (conj mid (count (distinct y)))))
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
          (sequential? rbf) (.setRBF cl (into-array smile.math.rbf.RadialBasisFunction (map k/smile-rbf rbf)))
          :else (.setRBF cl (k/smile-rbf rbf) number-of-basis))
        (.setNormalized normalize?))))

(wrap-classifier :smile rda {:keys [alpha priori tolerance]
                             :or {alpha 0.9 tolerance 1.0e-4}}
  (-> (RDA$Trainer. alpha)
      (.setPriori (m/seq->double-array priori))
      (.setTolerance tolerance)))

(wrap-classifier :smile random-forest {:keys [number-of-trees split-rule mtry node-size max-nodes subsample]
                                       :or {number-of-trees 500 split-rule :gini node-size 1 max-nodes 100 subsample 1.0}}
                 (let [mtry (or mtry (if (sequential? (first x)) (m/floor (m/sqrt (count (first x)))) 1))]
                   (-> (RandomForest$Trainer. ^int mtry ^int number-of-trees)
                       (.setSplitRule (or (split-rules split-rule) DecisionTree$SplitRule/GINI))
                       (.setMaxNodes max-nodes)
                       (.setSamplingRates subsample)
                       (.setNodeSize node-size))))

#_(wrap-classifier :xgboost xgboost xgboost-params xgboost-params)

(def ^:private multiclass-strategies {:one-vs-one SVM$Multiclass/ONE_VS_ONE
                                      :one-vs-all SVM$Multiclass/ONE_VS_ALL})

(def ^{:doc "List of multiclass strategies for [[svm]]"} multiclass-strategies-list (keys multiclass-strategies))

(wrap-classifier :smile svm {:keys [kernel ^double c-or-cp ^double cn strategy-for-multiclass class-weights tolerance ^int epochs]
                             :or {kernel (k/kernel :linear) c-or-cp 1.0 cn 1.0 strategy-for-multiclass :one-vs-one tolerance 1.0e-3 epochs 2}}
                 (let [kernel (k/smile-mercer kernel)
                       cn (or cn c-or-cp)
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
                     (do (.setWeights par (double-array weights) (int-array (range (count weights)))) par)
                     par)) bias)

(defn backend
  "Return name of backend library"
  [model] (pr/backend model))

(defn model-native
  "Return trained model as a backend class."
  [model] (pr/model-native model))

(defn data-native
  "Return data transformed for backend library."
  [model] (pr/data-native model))

(defn predict
  "Predict categories for given vector. If `posteriori?` is true returns also posteriori probability (default `false`)."
  ([model v] (pr/predict model v false))
  ([model v posteriori?] (pr/predict model v posteriori?)))

(defn predict-all
  "Predict categories for given sequence of vectors. If `posteriori?` is true returns also posteriori probability (default `false`)."
  ([model v] (pr/predict-all model v false))
  ([model v posteriori?] (pr/predict-all model v posteriori?)))

(defn train
  "Train another set of data for given classifier or force training already given data."
  ([model] (pr/train model))
  ([model xs ys] (pr/train model xs ys)))

(defn cv
  "Cross validation"
  ([model] (pr/cv model))
  ([model params] (pr/cv model params)))

(defn labels
  "Return labels"
  [ys] (first (labels-converters ys)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; validation metrics

(defn validate
  "Validate data against trained classifier. Same as [[test]]."
  {:metadoc/categories #{:vd}}
  [model tx ty]
  (let [pred (predict-all model tx)
        invalid (->> (mapv vector tx ty pred)
                     (filter #(apply not= (rest %))))
        invalid-cnt (count invalid)]
    {:prediction pred
     :invalid {:count invalid-cnt
               :data (map first invalid)
               :prediction (map #(nth % 2) invalid)
               :truth (map second invalid)}
     :stats {:accuracy (double (- 1.0 (/ invalid-cnt (double (count ty)))))}}))

(comment do

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
              {:accuracy avg :bootstrap (into (sorted-map) (map-indexed vector all))}))))

(defn confusion-map
  "Create confusion map where keys are pairs of `[truth-label prediction-label]`"
  {:metadoc/categories #{:vd}}
  [t p]
  (reduce (fn [m tpv]
            (if (contains? m tpv)
              (update m tpv clojure.core/inc)
              (assoc m tpv 1))) {} (map vector t p)))
