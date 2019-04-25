(ns fastmath.optimization
  (:import [org.apache.commons.math3.optim.nonlinear.scalar GoalType ObjectiveFunction]
           [org.apache.commons.math3.optim.univariate SearchInterval BrentOptimizer UnivariateObjectiveFunction]
           [org.apache.commons.math3.optim OptimizationData MaxEval MaxIter]
           [org.apache.commons.math3.analysis UnivariateFunction MultivariateFunction])
  (:require [fastmath.core :as m]))

(defn- interval
  ([x y] (SearchInterval. x y))
  ([x y init] (SearchInterval. x y init)))

(defn- brent
  ([{:keys [^double rel ^double abs]
     :or {rel 1.0e-6  abs 1.0e-10}}]
   (BrentOptimizer. rel abs))
  ([] (brent nil)))

(defn- wrap-univariate-function
  [f]
  (UnivariateObjectiveFunction. (reify UnivariateFunction
                                  (value [_ x] (f x)))))

(defn- wrap-multivariate-function
  [f]
  (ObjectiveFunction. (reify MultivariateFunction
                        (value [_ xs] (apply f xs)))))

(defn- prepare-params
  [{:keys [max-evals max-iters]}]
  {:max-evals (if max-evals (MaxEval. max-evals) (MaxEval/unlimited))
   :max-iters (if max-iters (MaxIter. max-iters) (MaxIter/unlimited))})



(let [b (brent)
      i (apply interval [-4 4])
      e (MaxEval/unlimited)
      f (wrap-univariate-function #(* m/PI (m/cos (dec (m/sin %)))))
      o (.optimize b (into-array OptimizationData [GoalType/MAXIMIZE i f e]))]
  (println (.getIterations b))
  [(.getPoint o) (.getValue o)])

