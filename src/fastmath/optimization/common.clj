(ns fastmath.optimization.common
  "Internal helpers shared by [[fastmath.optimization]], [[fastmath.optimization.acm]] and [[fastmath.optimization.lbfgsb]].

  The namespace is not part of the public API and may change without notice. It collects the option handling that every optimizer backend needs, so that each backend validates its input in exactly the same way.

  Option parsing: [[parse-goal]], [[resolve-vector-arg?]], [[throw-unknown]].

  Input validation: [[normalize-bounds]] checks and normalizes search bounds according to the capabilities of the given method.

  Function adaptation: [[->vector-fn]] turns a multi-arity function into a function of a single sequence.

  Apache Commons Math glue used by the ACM backends and by linear optimization: [[goal-type]], [[max-eval]], [[max-iter]], [[optimization-data]], [[multivariate-optimize]]."
  (:require [fastmath.core :as m])
  (:import [org.apache.commons.math3.optim.nonlinear.scalar GoalType]
           [org.apache.commons.math3.optim BaseOptimizer OptimizationData MaxEval MaxIter PointValuePair]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; options

(defn throw-unknown
  "Throws `ex-info` reporting an option value which is not one of the allowed ones.

  Parameters:

  - `what` (keyword): name of the option, for example `:goal`.
  - `value`: the rejected value.
  - `allowed` (collection): the accepted values.

  The exception data is `{what value, :allowed allowed}`. Never returns."
  [what value allowed]
  (throw (ex-info (str "Unknown " (pr-str what) " value " (pr-str value) ", allowed: " (pr-str allowed))
                  {what value :allowed allowed})))

(defn parse-goal
  "Validates an optimization goal.

  Parameters:

  - `goal`: `:minimize`, `:maximize` or `nil`.

  Returns the goal keyword; `nil` means `:minimize`. Throws `ex-info` for any other value.

  See also [[goal-type]]."
  [goal]
  (case goal
    nil :minimize
    (:minimize :maximize) goal
    (throw-unknown :goal goal #{:minimize :maximize})))

(defn resolve-vector-arg?
  "Resolves whether the optimized function receives its arguments as a single sequence.

  Parameters:

  - `method` (keyword): optimization method.
  - `vector-arg?`: `true`, `false` or `nil` (not given).

  An explicit boolean is returned as given. When `vector-arg?` is `nil`, the default is `false` for `:brent` (the function receives a number) and `true` for every other method (the function receives a sequence).

  Returns a boolean."
  [method vector-arg?]
  (if (nil? vector-arg?)
    (not= method :brent)
    (boolean vector-arg?)))

(defn ->vector-fn
  "Adapts a function to the single sequence calling convention.

  Parameters:

  - `f` (function): the function to adapt.
  - `vector-arg?` (boolean): `true` when `f` already takes one sequence.

  Returns `f` itself when `vector-arg?` is `true`, otherwise a function of one sequence which applies `f` to the elements."
  [f vector-arg?]
  (if vector-arg?
    f
    (fn [xs] (apply f xs))))

;; bounds

(def ^:private bounds-rules
  {:brent {:required? true :one-pair? true :finite? true :strict? true}
   :bobyqa {:required? true :finite? true}
   :cmaes {:required? true :finite? true}
   :nelder-mead {:finite? true :strict? true}
   :multidirectional-simplex {:finite? true :strict? true}
   :lbfgsb {:required? true}
   :powell {}
   :gradient {}
   :non-linear-gradient {}})

(defn- pair-problem
  "Returns a description of what is wrong with a single [lo hi] pair or nil."
  [method finite? strict?]
  (fn [[^double lo ^double hi]]
    (cond
      (or (Double/isNaN lo) (Double/isNaN hi)) "bounds must not be NaN"
      (m/> lo hi) "lower bound must not be greater than upper bound"
      (m/== lo Double/POSITIVE_INFINITY) "lower bound must not be +Inf"
      (m/== hi Double/NEGATIVE_INFINITY) "upper bound must not be -Inf"
      (and finite? (or (Double/isInfinite lo) (Double/isInfinite hi))) (str method " requires finite bounds")
      (and strict? (m/== lo hi)) (str method " requires lower bound less than upper bound"))))

(defn normalize-bounds
  "Validates search bounds for the given optimization method and returns them in the canonical form.

  Parameters:

  - `method` (keyword): `:brent`, `:bobyqa`, `:cmaes`, `:nelder-mead`, `:multidirectional-simplex`, `:lbfgsb`, `:powell`, `:gradient` or `:non-linear-gradient`.
  - `bounds`: a sequence of `[lo hi]` pairs, one per dimension, or `nil`. A flat `[lo hi]` is accepted for one dimension.
  - `initial`: the initial point (a number for one dimension or a sequence) or `nil`.

  Rules for all methods: no `nil` or NaN, `lo <= hi`, `lo` is not `+Inf` and `hi` is not `-Inf`, the number of pairs equals the length of `initial` when it is given, and infinite bounds require `initial` (no midpoint can be defined).

  Rules per method:

  - `:brent` - bounds are required, exactly one pair, finite with `lo < hi`.
  - `:bobyqa` and `:cmaes` - bounds are required and finite (`lo = hi` is allowed).
  - `:nelder-mead` and `:multidirectional-simplex` - bounds are optional (they only size the initial simplex); when given they are finite with `lo < hi`.
  - `:lbfgsb` - bounds are required, infinite values are allowed.
  - `:powell`, `:gradient` and `:non-linear-gradient` - bounds are optional and only checked structurally.

  Returns a vector of `[lo hi]` pairs of doubles, or `nil` when the method accepts missing bounds and `bounds` is `nil`.

  Throws `ex-info` with `{:method :bounds :initial :reason}` when a rule is violated, and with `{:method method, :allowed ...}` for an unknown method."
  [method bounds initial]
  (let [{:keys [required? one-pair? finite? strict?]} (or (bounds-rules method)
                                                          (throw-unknown :method method (set (keys bounds-rules))))
        fail (fn [reason] (throw (ex-info (str "Invalid bounds for " method ": " reason)
                                          {:method method :bounds bounds :initial initial :reason reason})))]
    (if (nil? bounds)
      (when required? (fail "bounds are required"))
      (let [flat? (and (sequential? bounds) (m/== 2 (count bounds)) (every? number? bounds))
            pairs (if flat? [bounds] bounds)]
        (when-not (and (sequential? pairs) (seq pairs))
          (fail "bounds must be a non-empty sequence of [lo hi] pairs"))
        (when-not (every? #(and (sequential? %) (m/== 2 (count %)) (every? number? %)) pairs)
          (fail "each bound must be a [lo hi] pair of numbers"))
        (let [res (mapv (fn [[lo hi]] [(double lo) (double hi)]) pairs)]
          (when (and one-pair? (not= 1 (count res)))
            (fail (str method " accepts exactly one [lo hi] pair")))
          (when-let [reason (some (pair-problem method finite? strict?) res)]
            (fail reason))
          (if (nil? initial)
            (when (some (fn [[^double lo ^double hi]] (or (Double/isInfinite lo) (Double/isInfinite hi))) res)
              (fail "infinite bounds require an initial point"))
            (let [dims (if (number? initial) 1 (count initial))]
              (when (not= dims (count res))
                (fail (str "number of bounds (" (count res) ") differs from the length of the initial point (" dims ")")))))
          res)))))

;; Apache Commons Math

(defn goal-type
  "Converts an optimization goal to the Apache Commons Math `GoalType`.

  Parameters:

  - `goal`: `:minimize`, `:maximize` or `nil` (minimize), see [[parse-goal]].

  Throws `ex-info` for any other value."
  ^GoalType [goal]
  (if (= :maximize (parse-goal goal))
    GoalType/MAXIMIZE
    GoalType/MINIMIZE))

(defn max-eval
  "Creates the Apache Commons Math limit on the number of function evaluations.

  Parameters:

  - `max-evals`: maximum number of evaluations; `nil` means 10000.

  Exceeding the limit during optimization throws an exception."
  ^MaxEval [max-evals]
  (MaxEval. (long (or max-evals 10000))))

(defn max-iter
  "Creates the Apache Commons Math limit on the number of iterations.

  Parameters:

  - `max-iters`: maximum number of iterations; `nil` means 10000.

  Exceeding the limit during optimization throws an exception."
  ^MaxIter [max-iters]
  (MaxIter. (long (or max-iters 10000))))

(defn optimization-data
  "Converts a collection of Apache Commons Math optimization data (goal, limits, objective, bounds, initial guess, ...) to the array expected by optimizers."
  ^"[Lorg.apache.commons.math3.optim.OptimizationData;" [data]
  (into-array OptimizationData data))

(defn multivariate-optimize
  "Runs an Apache Commons Math optimizer which returns a point and a value.

  Parameters:

  - `opt`: a `BaseOptimizer` returning a `PointValuePair`.
  - `data` (collection): optimization data, see [[optimization-data]].
  - `stats?` (boolean): return statistics along with the result.

  Returns `[point value]` where `point` is a vector. When `stats?` is `true`, returns a map with `:point`, `:value`, `:evaluations` and `:iterations`."
  [^BaseOptimizer opt data stats?]
  (let [^PointValuePair res (.optimize opt (optimization-data data))
        point (vec (.getPointRef res))
        value (.getValue res)]
    (if stats?
      {:point point
       :value value
       :evaluations (.getEvaluations opt)
       :iterations (.getIterations opt)}
      [point value])))
