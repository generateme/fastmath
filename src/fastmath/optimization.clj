(ns fastmath.optimization
  (:import [org.apache.commons.math3.optim.nonlinear.scalar GoalType ObjectiveFunction]
           [org.apache.commons.math3.optim.univariate SearchInterval BrentOptimizer UnivariateObjectiveFunction UnivariatePointValuePair]
           [org.apache.commons.math3.optim BaseOptimizer OptimizationData MaxEval MaxIter SimpleBounds InitialGuess PointValuePair]
           [org.apache.commons.math3.analysis UnivariateFunction MultivariateFunction]
           [org.apache.commons.math3.optim.nonlinear.scalar.noderiv BOBYQAOptimizer])
  (:require [fastmath.core :as m]
            [fastmath.random :as r]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(def ^:private univariate-set #{:brent})
(def ^:private multivariate-set #{:bobyqa})

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

(def ^:private optimizers
  {:brent brent
   :bobyqa bobyqa})

(defn- optimizer
  [config]
  ((optimizers (:method config)) config))


(defn- wrap-univariate-function [f] (UnivariateObjectiveFunction. (reify UnivariateFunction
                                                                    (value [_ x] (f x)))))

(defn- wrap-multivariate-function [f] (ObjectiveFunction. (reify MultivariateFunction
                                                            (value [_ xs] (apply f xs)))))

(defn- evals [evals] (if evals (MaxEval. evals) (MaxEval/unlimited)))
(defn- iters [iters] (if iters (MaxIter. iters) (MaxIter/unlimited)))

(defn- infer-lo-high
  [lo high]
  [(or lo m/EPSILON)
   (or high (- 1.0 m/EPSILON))])

(defn- search-interval
  [[lo high init]]
  (let [[lo high] (infer-lo-high lo high)]
    (if init
      (SearchInterval. lo high init)
      (SearchInterval. lo high))))

(defn- bounds-and-guess
  [bounds]
  (let [[lo high init] (map #(if (sequential? %) % [%]) bounds)
        sb (SimpleBounds. (double-array lo) (double-array high))]
    [sb (InitialGuess. (double-array (if init init (map (fn [^double l ^double h]
                                                          (* 0.5 (+ l h))) lo high))))]))

(defn- bounds
  [{:keys [bounds method]}]
  (if (univariate-set method)
    (search-interval bounds)
    (bounds-and-guess bounds)))

(defn- wrap-function
  [f config]
  (if (univariate-set (:method config))
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

(defn optimize
  ([method f] (optimize method f {}))
  ([method f {:keys [max-evals max-iters goal stats?] :as config}]
   (let [[lo] (:bounds config)
         config (assoc config :method method :dim (if (and lo (sequential? lo))
                                                    (count lo) 1))
         evals (evals max-evals)
         iters (iters max-iters)
         goal (if (= goal :maximize) GoalType/MAXIMIZE GoalType/MINIMIZE)
         bounds (bounds config)
         f (wrap-function f config)
         ^BaseOptimizer optimizer (optimizer config)
         res (parse-result (.optimize optimizer (into-array OptimizationData (flatten [evals iters goal bounds f]))))]
     (if-not stats?
       res
       {:result res
        :evaluations (.getEvaluations optimizer)
        :iterations (.getIterations optimizer)}))))

(defn minimize
  ([method f config] (optimize method f (assoc config :goal :minimize)))
  ([method f] (minimize method f {})))

(defn maximize
  ([method f config] (optimize method f (assoc config :goal :maximize)))
  ([method f] (maximize method f {})))


(time (maximize :brent #(m/cos %) {:bounds [-3 3]}))


(defn s2d ^double [^double x ^double y] (+ (* 100 (m/sqrt (m/abs (- y (* 0.01 x x)))))
                                           (* 0.01 (m/abs (+ 10.0 x)))))

(defn h ^double [^double x ^double y] (+ (m/sq (+ (* x x) y -11))
                                         (m/sq (+ x (* y y) -7))))

(defn d5 [a b c d e] (reduce #(+ %1 (m/sq %2)) 0.0 [a b c d e]))

(first (sort-by first (take 100000 (map #(vector (s2d %1 %2) %1 %2) (repeatedly #(r/drand -15 -5)) (repeatedly #(r/drand -3 3))))))
;; => [1.9210981963566007 -5.751481824637489 0.3304425131054902]

(time (minimize :bobyqa d5 {:bounds [[-5 -5 -5 -5 -5] [5 5 5 5 5]] :stats? true}))
