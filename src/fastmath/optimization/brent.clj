(ns fastmath.optimization.brent
  (:require [fastmath.optimization.common :as common])
  (:import [org.apache.commons.math3.optim.univariate BrentOptimizer UnivariatePointValuePair SearchInterval UnivariateObjectiveFunction BracketFinder]
           [org.apache.commons.math3.analysis UnivariateFunction]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

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

(defn brent
  [f {:keys [bounds init max-evals max-iters ^double rel ^double abs find-bracket vector-arg? goal]
      :or {rel 1.0e-8 abs 1.0e-10 goal :minimize}}]
  (when-not bounds (throw (ex-info "No bounds defined" nil)))
  (let [goal (common/goal-type goal)
        uf (univariate-function f vector-arg?)
        uof (univariate-objective-function uf)
        [^double lo ^double high] (process-bounds bounds find-bracket uf goal)
        si (search-interval lo high init)
        me (common/max-eval max-evals)
        mi (common/max-iter max-iters)
        ^BrentOptimizer bo (BrentOptimizer. rel abs)]
    (parse-result (.optimize bo (common/optimization-data [goal uof si me mi])) vector-arg?)))

