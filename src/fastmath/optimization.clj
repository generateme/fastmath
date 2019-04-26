(ns fastmath.optimization
  (:import [org.apache.commons.math3.optim.nonlinear.scalar GoalType ObjectiveFunction]
           [org.apache.commons.math3.optim.univariate SearchInterval BrentOptimizer UnivariateObjectiveFunction UnivariatePointValuePair]
           [org.apache.commons.math3.optim BaseOptimizer OptimizationData MaxEval MaxIter SimpleBounds InitialGuess PointValuePair]
           [org.apache.commons.math3.analysis UnivariateFunction MultivariateFunction]
           [org.apache.commons.math3.optim.nonlinear.scalar.noderiv BOBYQAOptimizer PowellOptimizer])
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(def ^:private univariate-set #{:brent})
(def ^:private multivariate-set #{:bobyqa :powell})

(defn- brent
  [{:keys [^double rel ^double abs]
    :or {rel 1.0e-6  abs 1.0e-10}}]
  (BrentOptimizer. rel abs))

(defn- bobyqa
  ([{:keys [number-of-points ^int dim initial-radius stopping-radius]
     :or {initial-radius BOBYQAOptimizer/DEFAULT_INITIAL_RADIUS
          stopping-radius BOBYQAOptimizer/DEFAULT_STOPPING_RADIUS}}]
   (let [number-of-points (or number-of-points (/ (+ 3 (* 3 dim)) 2))]
     (BOBYQAOptimizer. number-of-points initial-radius stopping-radius))))

(defn- powell
  [{:keys [^double rel ^double abs]
    :or {rel 1.0e-6  abs 1.0e-10}}]
  (PowellOptimizer. rel abs))

(def ^:private optimizers
  {:brent brent
   :bobyqa bobyqa
   :powell powell})

(defn- wrap-univariate-function [f] (UnivariateObjectiveFunction. (reify UnivariateFunction
                                                                    (value [_ x] (f x)))))

(defn- wrap-multivariate-function [f] (ObjectiveFunction. (reify MultivariateFunction
                                                            (value [_ xs] (apply f xs)))))

(defn- wrap-evals [evals] (if evals (MaxEval. evals) (MaxEval/unlimited)))
(defn- wrap-iters [iters] (if iters (MaxIter. iters) (MaxIter/unlimited)))

(defn- infer-lo-high
  [lo high]
  [(or lo m/EPSILON)
   (or high (- 1.0 m/EPSILON))])

(defn- search-interval
  [[lo high] init]
  (let [[lo high] (infer-lo-high lo high)]
    (if init
      (SearchInterval. lo high init)
      (SearchInterval. lo high))))

(defn- ensure-bounds-as-vectors
  [bounds]
  (map #(if (sequential? %) (vec %) [%]) bounds))

(defn- multi-bounds
  [bounds]
  (let [[lo high] (ensure-bounds-as-vectors bounds)]
    (SimpleBounds. (double-array lo) (double-array high))))

(defn- initial-guess
  [bounds initial]
  (InitialGuess. (double-array (if initial (if (sequential? initial) initial [initial])
                                   (let [[lo high] (ensure-bounds-as-vectors bounds)]
                                     (map (fn [^double l ^double h] (* 0.5 (+ l h))) lo high))))))

(defn- wrap-bounds
  [method bounds initial]
  (when (and bounds (multivariate-set method))
    (multi-bounds bounds)))

(defn wrap-initial
  [method bounds initial]
  (if (univariate-set method)
    (search-interval bounds initial)
    (initial-guess bounds initial)))

(defn- wrap-function
  [method f config]
  (if (univariate-set method)
    (wrap-univariate-function f)
    (wrap-multivariate-function f)))

(defn- parse-result
  [res]
  (condp instance? res
    UnivariatePointValuePair (let [^UnivariatePointValuePair res res]
                               [(.getPoint res) (.getValue res)])
    PointValuePair (let [^PointValuePair res res]
                     [(seq (.getPointRef res)) (.getValue res)])
    res))

(defn optimizer
  [method f {:keys [max-evals max-iters goal initial bounds stats?] :as config}]
  (assert (or (not (nil? bounds))
              (not (nil? initial))) "Provide search bounds or initial.")
  (let [lo (if initial initial (first bounds))
        config (assoc config :dim (if (and lo (sequential? lo))
                                    (count lo) 1))
        
        b (wrap-bounds method bounds initial)
        base-opt-data [(wrap-evals max-evals)
                       (wrap-iters max-iters)
                       (if (= goal :maximize) GoalType/MAXIMIZE GoalType/MINIMIZE)
                       (wrap-function method f config)]
        base-opt-data (if (and b (not (#{:powell} method))) (conj base-opt-data b) base-opt-data) 
        
        builder (optimizers method)        
        ^BaseOptimizer optimizer (builder config)]
    
    (fn [init] (let [res (->> (conj base-opt-data (wrap-initial method bounds init))
                             (into-array OptimizationData)
                             (.optimize optimizer)
                             (parse-result))]
                (if-not stats?
                  res
                  {:result res
                   :evaluations (.getEvaluations optimizer)
                   :iterations (.getIterations optimizer)})))))

(defn minimizer [method f config] (optimizer method f (assoc config :goal :minimize)))
(defn maximizer [method f config] (optimizer method f (assoc config :goal :maximize)))

(defmacro ^:private with-optimizer [opt m f c] `((~opt ~m ~f ~c) (:initial ~c)))

(defn optimize [method f config] (with-optimizer optimizer method f config))
(defn minimize [method f config] (with-optimizer minimizer method f config))
(defn maximize [method f config] (with-optimizer maximizer method f config))

(defn- goal-comparator
  [goal]
  (if (= goal :minimize)
    #(< ^double %1 ^double %2)
    #(> ^double %1 ^double %2)))

(defn- generate-points
  [f [lo high] goal N n]
  (let [^long dim (if (sequential? lo) (count lo) 1)
        N (max 10 ^int N)
        n (max 10 (m/floor (* ^double n N)))
        gen (r/jittered-sequence-generator (if (<= dim 4) :r2 :sobol) dim 0.2)
        inter (if (== dim 1) m/lerp v/einterpolate)
        genf (if (== dim 1) f (partial apply f))]
    (->> (take N gen)
         (map #(let [p (inter lo high %)]
                 [(genf p) p]))
         (sort-by first (goal-comparator goal))
         (take n)
         (map second))))

(defn- scan-and-
  "For cheap functions, scan domain and bruteforcely search for set of minimal values, then use part of this values as initial points for parallel optimization.

  Additional parameters in config:

  * N - number of points to scan per dimension (default: 10, minumum 10)
  * n - fraction of total points N) used for optimization (default: 0.1, minimum 10)"
  [optimizer goal method f {:keys [bounds ^int N ^double n]
                            :or {N 10 n 0.1}
                            :as config}]
  (assert (not (nil? bounds)) "Provide search bounds.")
  (let [goal (or goal (get config :goal :minimize))
        samples (generate-points f bounds goal N n)
        opt (repeatedly (count samples) #(optimizer method f config))
        ^int tk (get config :take 1)
        taker (if (> tk 1) (partial take tk) first)]
    (->> (pmap #(%1 %2) opt samples)
         (sort-by second (goal-comparator goal))
         (taker))))

(def scan-and-optimize (partial scan-and- optimizer nil))
(def scan-and-minimize (partial scan-and- minimizer :minimize))
(def scan-and-maximize (partial scan-and- maximizer :maximize))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn bfn6 ^double [^double x ^double y] (+ (* 100 (m/sqrt (m/abs (- y (* 0.01 x x)))))
                                            (* 0.01 (m/abs (+ 10.0 x)))))

(defn h ^double [^double x ^double y] (+ (m/sq (+ (* x x) y -11))
                                         (m/sq (+ x (* y y) -7))))

(defn d5 [a b c d e] (reduce #(+ ^double %1 (m/sq %2)) 0.0 [a b c d e]))

(time (scan-and-maximize :bobyqa bfn6 {:bounds [[-15 -3] [15 3]] :N 100 :n 0.2}))

(minimize :powell bfn6 {:bounds [[-15 -3] [15 3]] :initial [-13.217532309719662 0.9280007414397415]})

(time (scan-and-optimize :powell #(m/cos %) {:bounds [-3 3] :initial -2 :N 100 :n 0.2 :goal :maximize}))


;; => [1.9210981963566007 -5.751481824637489 0.3304425131054902]



