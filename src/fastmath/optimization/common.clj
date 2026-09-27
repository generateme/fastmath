(ns fastmath.optimization.common
  (:import [org.apache.commons.math3.optim.nonlinear.scalar GoalType]
           [org.apache.commons.math3.optim MaxEval MaxIter OptimizationData]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)


;; Apache Commons Math stuff

;; common

(defn goal-type
  [goal]
  (case goal
    :minimize GoalType/MINIMIZE
    :maximize GoalType/MAXIMIZE))

(defn max-eval [max-evals] (if max-evals (MaxEval. max-evals) (MaxEval/unlimited)))
(defn max-iter [max-iters] (if max-iters (MaxIter. max-iters) (MaxIter/unlimited)))

(defn optimization-data ^"[Lorg.apache.commons.math3.optim.OptimizationData;" [data] (into-array OptimizationData data))
