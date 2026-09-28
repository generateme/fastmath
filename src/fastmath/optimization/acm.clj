(ns fastmath.optimization.acm
  "Apache Commons Math optimizers"
  (:require [fastmath.core :as m])
  (:import [org.apache.commons.math3.optim.univariate BrentOptimizer UnivariatePointValuePair SearchInterval UnivariateObjectiveFunction BracketFinder]
           [org.apache.commons.math3.optim.nonlinear.scalar GoalType]
           [org.apache.commons.math3.optim MaxEval MaxIter OptimizationData]
           [org.apache.commons.math3.optim.nonlinear.scalar.noderiv BOBYQAOptimizer]
           [org.apache.commons.math3.analysis UnivariateFunction]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; common

(defn goal-type
  [goal]
  (if-not goal
    GoalType/MINIMIZE
    (case goal
      :minimize GoalType/MINIMIZE
      :maximize GoalType/MAXIMIZE)))

(defn max-eval [max-evals] (if max-evals (MaxEval. max-evals) (MaxEval/unlimited)))
(defn max-iter [max-iters] (if max-iters (MaxIter. max-iters) (MaxIter/unlimited)))

(defn optimization-data ^"[Lorg.apache.commons.math3.optim.OptimizationData;" [data] (into-array OptimizationData data))


;; univariate

(defn search-interval
  ([^double low ^double high]
   (SearchInterval. low high))
  ([^double low ^double high init]
   (if init
     (SearchInterval. low high (double init))
     (search-interval low high))))

(defn univariate-function
  [f vector-arg?]
  (if vector-arg?
    (reify UnivariateFunction (value [_ x] (f [x])))
    (reify UnivariateFunction (value [_ x] (f x)))))

(defn univariate-objective-function [^UnivariateFunction f] (UnivariateObjectiveFunction. f))

(defn process-bounds
  [bounds find-bracket ^UnivariateFunction uf goal]
  (let [bounds (if (sequential? (first bounds)) (first bounds) bounds)]
    (if find-bracket
      (let [{:keys [^double grow-limit ^long max-evals]
             :or {grow-limit 100.0 max-evals 500}} (if (map? find-bracket) find-bracket {})
            [^double x ^double y] bounds
            ^BracketFinder bf (BracketFinder. grow-limit max-evals)]
        (.search bf uf goal x y)
        [(.getLo bf) (.getHi bf)])
      bounds)))

(defn parse-result
  [^UnivariatePointValuePair res vector-arg?]
  (if vector-arg?
    [[(.getPoint res)] (.getValue res)]
    [(.getPoint res) (.getValue res)]))

(defn brent-data
  [f {:keys [bounds init max-evals max-iters ^double rel ^double abs find-bracket vector-arg? goal]
      :or {rel 1.0e-8 abs 1.0e-10}}]
  (when-not bounds (throw (ex-info "No bounds defined" nil)))
  (let [goal (goal-type goal)
        uf (univariate-function f vector-arg?)
        uof (univariate-objective-function uf)
        [^double lo ^double high] (process-bounds bounds find-bracket uf goal)
        si (search-interval lo high init)
        me (max-eval max-evals)
        mi (max-iter max-iters)]
    {:goal goal
     :univariate-objective-function uof
     :search-interval si
     :max-eval me
     :max-iter mi
     :rel rel :abs abs
     :vector-arg? vector-arg?}))

(defn update-init
  [{:keys [^SearchInterval search-interval] :as brent-data} ^double init]
  (assoc brent-data :search-interval (SearchInterval. (.getMin search-interval)
                                                      (.getMax search-interval)
                                                      init)))

(def ^:private extract-optimization-data (juxt :goal :univariate-objective-function :search-interval :max-eval :max-iter))

(defn brent
  [{:keys [^double rel ^double abs vector-arg?] :as brent-data}]
  (let [^BrentOptimizer bo (BrentOptimizer. rel abs)]
    (-> (extract-optimization-data brent-data)
        (optimization-data)
        (->> (.optimize bo))
        (parse-result vector-arg?))))

;; multivariate

(defn bobyqa
  [f {:keys [bounds init max-evals max-iters vector-arg? goal number-of-points ^double initial-radius ^double stopping-radius]
      :or {initial-radius BOBYQAOptimizer/DEFAULT_INITIAL_RADIUS
           stopping-radius BOBYQAOptimizer/DEFAULT_STOPPING_RADIUS}}]
  (when-not bounds (throw (ex-info "No bounds defined" nil)))
  (let [n (count bounds)]
    (when (m/< n 2) (throw (ex-info "Number of dimensions should be equal or greater than 2" {:n n :bounds bounds})))
    (let [number-of-points (long (or number-of-points (m/inc (m/* 2 n))))
          goal (goal-type goal)]))
  )
