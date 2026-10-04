(ns fastmath.optimization
  "Minimization and maximization of functions with various optimization methods, Bayesian optimization and linear programming.

  Functions are optimized by name of a method: call [[minimize]] or [[maximize]] with the method, the function and an options map. [[minimizer]] and [[maximizer]] create a function which runs the optimization for a given initial point, [[scan-and-minimize]] and [[scan-and-maximize]] scan the search domain first and then run many optimizations in parallel from the best points. [[bayesian-optimization]] optimizes expensive functions and [[linear-optimization]] solves linear programs.

  Methods:

  - `:brent` - one dimension, local, derivative free (Apache Commons Math).
  - `:bobyqa` - box constrained, derivative free, two or more dimensions (Apache Commons Math).
  - `:cmaes` - box constrained evolution strategy, derivative free, stochastic (Apache Commons Math).
  - `:nelder-mead`, `:multidirectional-simplex`, `:powell` - unconstrained, derivative free (Apache Commons Math).
  - `:gradient` (also `:non-linear-gradient`) - unconstrained conjugate gradient, numerical or given gradient (Apache Commons Math).
  - `:lbfgsb` - box constrained quasi-Newton L-BFGS-B, numerical or given gradient (see [[fastmath.optimization.lbfgsb]]).

  The functions of the methods are documented in [[fastmath.optimization.acm]] and [[fastmath.optimization.lbfgsb]], where all their options are described.

  Options common to all methods:

  - `:bounds` - sequence of `[lo hi]` pairs, one for each dimension (`[lo hi]` for one dimension). Required by `:brent`, `:bobyqa`, `:cmaes` and `:lbfgsb`, optional for the other methods where they only set the size of the initial simplex or the initial point. Bounds are validated: no NaN, `lo <= hi`, their number matches the initial point. `:brent` needs exactly one finite interval with `lo < hi`, `:bobyqa` and `:cmaes` finite bounds, the simplex methods finite bounds with `lo < hi`. Infinite bounds are allowed by `:lbfgsb` when `:initial` is given.
  - `:initial` - the initial point, default: the middle of the bounds.
  - `:goal` - `:minimize` (default) or `:maximize`. [[minimize]] and [[maximize]] set it.
  - `:vector-arg?` - `true`: the function receives the point as one sequence, `false`: as separate arguments. Default: `true`, but `false` for `:brent`, which receives a number. Whichever the form, the point is treated as a sequence of numbers (it can be an array, a vector or a lazy sequence depending on the method).
  - `:gradient` - function of the point, always one sequence, returning the gradient of the function as a sequence of numbers. It is the gradient of the function itself, also when maximizing. Used by `:lbfgsb` and `:gradient` only and ignored by the other methods. Default: finite differences with step `:gradient-h`.
  - `:max-evals`, `:max-iters` - limits of the numbers of evaluations and iterations. Exceeding a limit throws an exception, except for `:lbfgsb` (maximum of iterations is not an error, no limit of evaluations) and `:cmaes`.
  - `:stats?` - return a map with additional information instead of `[point value]`. The map depends on the method: `:point`, `:value`, `:evaluations` and `:iterations` for the Apache Commons Math methods and for linear optimization, `:point`, `:value`, `:iterations`, `:gradient` and `:status` for `:lbfgsb`.

  The result is always `[point value]`, where `value` is the value of the function, also when maximizing.

  Unknown methods, goals and option values (for example formulas or line searches) throw `ex-info` with the allowed values in the exception data.

  Scan and optimize: the `scan-and-...` functions evaluate the function at `:N` points of a jittered low discrepancy sequence, start the optimization from the best `:n` fraction of them in parallel, and return the best result. Optimization runs which fail with an exception are skipped.

  Bayesian optimization: [[bayesian-optimization]] can be used for optimizing expensive to evaluate black box functions. Refer to this [article](http://krasserm.github.io/2018/03/21/bayesian-optimization/) or this [article](https://nextjournal.com/a/LKqpdDdxiggRyHhqDG5FH?token=Ss1Qq3MzHWN8ZyEt9UC1ZZ)

  Linear optimization: [[linear-optimization]] solves linear programs with the simplex method."
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.kernel :as k]
            [fastmath.interpolation.gp :as gp]
            [fastmath.optimization.common :as common]
            [fastmath.optimization.lbfgsb :as lbfgsb]
            [fastmath.optimization.acm :as acm])
  (:import [org.apache.commons.math3.optim BaseOptimizer]
           [org.apache.commons.math3.optim.linear LinearObjectiveFunction LinearConstraint
            Relationship LinearConstraintSet SimplexSolver NonNegativeConstraint PivotSelectionRule]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^:private optimizers
  {:lbfgsb lbfgsb/lbfgsb
   :brent acm/brent
   :bobyqa acm/bobyqa
   :cmaes acm/cmaes
   :nelder-mead acm/nelder-mead
   :multidirectional-simplex acm/multidirectional-simplex
   :powell acm/powell
   :gradient acm/non-linear-gradient
   :non-linear-gradient acm/non-linear-gradient})

(defn- optimizer
  "Returns the optimizer function of the method or throws `ex-info`."
  [method]
  (or (optimizers method)
      (common/throw-unknown :method method (set (keys optimizers)))))

(defn- initial-point-optimizer
  [method f options goal]
  (let [optimize-fn (optimizer method)
        options (assoc options :goal goal)]
    (fn [initial] (optimize-fn f (assoc options :initial initial)))))

(defn minimizer
  "Creates a function which minimizes the function `f` from a given initial point.

  Parameters:

  - `method` (keyword): optimization method, see [[fastmath.optimization]].
  - `f` (function): the function to minimize.
  - `options` (map): options of the method, see [[fastmath.optimization]]. `:goal` is overridden.

  Returns a function of one argument, the initial point (a sequence of numbers, a number for `:brent`) or `nil` for the default one. The initial point replaces `:initial` of `options`. The function returns the result of [[minimize]]. The function has no zero-arity.

  Throws `ex-info` for an unknown method. The other options are validated when the returned function is called.

  See also [[maximizer]], [[minimize]], [[scan-and-minimize]]."
  [method f options]
  (initial-point-optimizer method f options :minimize))

(defn maximizer
  "Creates a function which maximizes the function `f` from a given initial point.

  Parameters:

  - `method` (keyword): optimization method, see [[fastmath.optimization]].
  - `f` (function): the function to maximize.
  - `options` (map): options of the method, see [[fastmath.optimization]]. `:goal` is overridden.

  Returns a function of one argument, the initial point (a sequence of numbers, a number for `:brent`) or `nil` for the default one. The initial point replaces `:initial` of `options`. The function returns the result of [[maximize]]. The function has no zero-arity.

  Throws `ex-info` for an unknown method. The other options are validated when the returned function is called.

  See also [[minimizer]], [[maximize]], [[scan-and-maximize]]."
  [method f options]
  (initial-point-optimizer method f options :maximize))

;;

(defn optimize
  "Optimizes the function `f` with the given method. The goal is taken from the options.

  Parameters:

  - `method` (keyword): optimization method, see [[fastmath.optimization]].
  - `f` (function): the function to optimize.
  - `options` (map): options of the method, see [[fastmath.optimization]]. `:goal` is `:minimize` (default) or `:maximize`.

  Returns `[point value]`, or a map when `:stats?` is `true`.

  Throws `ex-info` for an unknown method, goal or option value, and for invalid bounds.

  See also [[minimize]], [[maximize]]."
  [method f options]
  ((optimizer method) f options))

(defn minimize
  "Minimizes the function `f` with the given method.

  Parameters:

  - `method` (keyword): optimization method, one of `:brent`, `:bobyqa`, `:cmaes`, `:nelder-mead`, `:multidirectional-simplex`, `:powell`, `:gradient` and `:lbfgsb`.
  - `f` (function): the function to minimize.
  - `options` (map): options of the method, see [[fastmath.optimization]]. `:goal` is overridden.

  Returns `[point value]`, or a map when `:stats?` is `true`. `point` is a vector (a number for `:brent`) and `value` is the value of `f` at the point.

  Throws `ex-info` for an unknown method or option value, and for invalid bounds. Exceeding `:max-evals` or `:max-iters` throws an exception.

  See also [[maximize]], [[minimizer]], [[scan-and-minimize]]."
  [method f options]
  (optimize method f (assoc options :goal :minimize)))

(defn maximize
  "Maximizes the function `f` with the given method.

  Parameters:

  - `method` (keyword): optimization method, one of `:brent`, `:bobyqa`, `:cmaes`, `:nelder-mead`, `:multidirectional-simplex`, `:powell`, `:gradient` and `:lbfgsb`.
  - `f` (function): the function to maximize.
  - `options` (map): options of the method, see [[fastmath.optimization]]. `:goal` is overridden.

  Returns `[point value]`, or a map when `:stats?` is `true`. `point` is a vector (a number for `:brent`) and `value` is the value of `f` at the point.

  Throws `ex-info` for an unknown method or option value, and for invalid bounds. Exceeding `:max-evals` or `:max-iters` throws an exception.

  See also [[minimize]], [[maximizer]], [[scan-and-maximize]]."
  [method f options]
  (optimize method f (assoc options :goal :maximize)))

;;

(defn- goal-comparator [goal] (if (= goal :minimize) m/< m/>))

(defn- generate-points
  "Evaluates `f` (a function of a sequence) at `N` (at least `4.5 + d log2 d` for `d` dimensions) points of the bounds and returns the points sorted from the best one."
  [f bounds goal N jitter]
  (let [dim (count bounds)
        lo (map first bounds)
        hi (map second bounds)
        N (long (m/max (m/+ 4.5 (m/* dim (m/log2 dim))) (long N)))
        gen (r/jittered-sequence-generator (if (m/< dim 15) :r2 :sobol) dim jitter)]
    (->> (if (m/one? dim) (map vector gen) gen)
         (map (fn [v] (let [p (v/einterpolate lo hi v)] [(f p) p])))
         (filter (comp m/valid-double? first))
         (take N)
         (sort-by first (goal-comparator goal))
         (map second))))

(defn- wrap-optimizer
  [optimizer]
  (fn [initial]
    (try
      (optimizer initial)
      (catch Exception _ nil))))

(defn scan-and-optimize
  "Scans the search domain with a low discrepancy sequence and optimizes the function in parallel from the best scanned points.

  The function is evaluated at `:N` points of a jittered low discrepancy sequence (see [[fastmath.random/jittered-sequence-generator]]). The best of them are the initial points of the optimization with the given method. Use it for cheap functions with many local extrema.

  Parameters:

  - `method` (keyword): optimization method, see [[fastmath.optimization]].
  - `f` (function): the function to optimize.
  - `opts` (map): all options of the method, see [[fastmath.optimization]] (`:initial` is replaced by the scanned points) and:
    - `:bounds` (required) - the domain to scan, validated for the method. Infinite bounds are not allowed.
    - `:goal` - `:minimize` (default) or `:maximize`.
    - `:N` - number of scanned points, default: `100`. At least `4.5 + d log2 d` points are used for `d` dimensions.
    - `:n` - number of optimization runs: a fraction of `:N` when not greater than `1.0` (default: `0.05`), otherwise the number itself. At least one run is made.
    - `:jitter` - jitter of the sequence generator, default: `0.25`.
    - `:parallel?` - run the optimizations in parallel, default: `true`. The function has to be thread safe.
    - `:vector-arg?` - how the function receives the point, see [[fastmath.optimization]].
    - `:take-last-n` - when greater than `1`, return the best `:take-last-n` results as a sequence, default: `0`.
    - `:stats?` - return maps with additional information instead of `[point value]`.

  Returns the best result of the form of [[minimize]] (`[point value]` or a map). With `:take-last-n` returns a sequence of results sorted from the best one. Returns `nil` (an empty sequence) when all optimization runs failed.

  Runs which throw an exception, for example when a limit is exceeded, are skipped. Exceptions of the evaluation of the scanned points, an unknown method and invalid bounds or options throw.

  See also [[scan-and-minimize]], [[scan-and-maximize]], [[minimize]]."
  [method f {:keys [bounds N n jitter parallel? vector-arg? goal take-last-n stats?]
             :or {parallel? true}
             :as opts}]
  (let [N (long (or N 100))
        n (double (or n 0.05))
        jitter (double (or jitter 0.25))
        take-last-n (long (or take-last-n 0))
        optimize-fn (optimizer method)
        goal (common/parse-goal goal)
        bounds (or (common/normalize-bounds method bounds nil)
                   (throw (ex-info "Provide search bounds." {:method method :bounds bounds})))
        vector-arg? (common/resolve-vector-arg? method vector-arg?)
        opts (assoc opts :bounds bounds :goal goal :vector-arg? vector-arg?)
        run (wrap-optimizer (fn [initial] (optimize-fn f (assoc opts :initial initial))))
        samples (generate-points (common/->vector-fn f vector-arg?) bounds goal N jitter)
        nbest (long (m/max 1 (if (m/> n 1.0) n (m/floor (m/* n N)))))
        taker (if (m/> take-last-n 1) (partial take take-last-n) first)
        mapper (if parallel? pmap map)
        sort-selector (if stats? :value second)]
    (->> (mapper run samples)
         (filter identity)
         (take nbest)
         (sort-by sort-selector (goal-comparator goal))
         (taker))))

(defn scan-and-minimize
  "Scans the search domain and minimizes the function in parallel from the best scanned points.

  Parameters:

  - `method` (keyword): optimization method, see [[fastmath.optimization]].
  - `f` (function): the function to minimize.
  - `opts` (map): `:bounds` (required), the options of the method and the scan options `:N`, `:n`, `:jitter`, `:parallel?`, `:take-last-n`, see [[scan-and-optimize]]. `:goal` is overridden.

  Returns the best result of the form of [[minimize]], or `nil` when all runs failed.

  See also [[scan-and-maximize]], [[scan-and-optimize]]."
  [method f opts]
  (scan-and-optimize method f (assoc opts :goal :minimize)))

(defn scan-and-maximize
  "Scans the search domain and maximizes the function in parallel from the best scanned points.

  Parameters:

  - `method` (keyword): optimization method, see [[fastmath.optimization]].
  - `f` (function): the function to maximize.
  - `opts` (map): `:bounds` (required), the options of the method and the scan options `:N`, `:n`, `:jitter`, `:parallel?`, `:take-last-n`, see [[scan-and-optimize]]. `:goal` is overridden.

  Returns the best result of the form of [[maximize]], or `nil` when all runs failed.

  See also [[scan-and-minimize]], [[scan-and-optimize]]."
  [method f opts]
  (scan-and-optimize method f (assoc opts :goal :maximize)))

;; bo

(defmulti ^:private utility-function (fn [t & _] t))

(defmethod utility-function :default [t _]
  (common/throw-unknown :utility-function-type t #{:ucb :ei :poi}))

(defmethod utility-function :ucb
  [_ ^double kappa]
  (fn [gp x _]
    (let [[^double mean ^double stddev] (gp/predict gp x true)]
      (m/+ mean (m/* kappa stddev)))))

(defmethod utility-function :ei
  [_ ^double xi]
  (fn [gp x ^double y-max]
    (let [[^double mean ^double stddev] (gp/predict gp x true)
          diff (m/- mean y-max xi)
          z (m// diff stddev)]
      (m/+ (m/* diff (r/cdf r/default-normal z))
           (m/* stddev (r/pdf r/default-normal z))))))

(defmethod utility-function :poi
  [_ ^double xi]
  (fn [gp x ^double y-max]
    (let [[^double mean ^double stddev] (gp/predict gp x true)]
      (r/cdf r/default-normal (m// (m/- mean y-max xi) stddev)))))

(defn- gen-sequence
  [init-points bounds jitter]
  (let [dims (count bounds)
        int-fn (if (m/one? dims)
                 #(vector (m/lerp (ffirst bounds) (second (first bounds)) %))
                 #(v/einterpolate (mapv first bounds) (mapv second bounds) %))]
    (->> (r/jittered-sequence-generator (if (m/< dims 15) :r2 :sobol) dims jitter)
         (take init-points)
         (map int-fn))))

(defn- initial-values
  [f init-points bounds jitter]
  (let [pts (if (sequential? init-points)
              init-points
              (gen-sequence init-points bounds jitter))]
    [pts (map f pts)]))

(defn- bayesian-step-fn
  [f util-fn warm-up bounds gp jitter optimizer optimizer-params]
  ;; the utility function is a function of one sequence and the result has to be a single point
  (let [params (merge optimizer-params {:N warm-up :n 0.02
                                        :bounds bounds :jitter jitter :parallel? false
                                        :vector-arg? true :stats? false :take-last-n 0})]
    (fn [{:keys [x ^double y xs ys]}]
      (let [curr-gp (gp xs ys)
            curr-util (fn [r] (util-fn curr-gp (vec r) y))
            ;; unconstrained optimizers (powell, nelder-mead, ...) can leave the bounds: the point is moved back
            bx (mapv (fn [^double x [^double lo ^double hi]] (m/constrain x lo hi))
                     (or (first (scan-and-maximize optimizer curr-util params))
                         (throw (ex-info "No maximum of the utility function found" {:optimizer optimizer :bounds bounds})))
                     bounds)
            by (double (f bx))
            nxs (conj xs bx)
            nys (conj ys by)]
        {:x (if (m/> by y) bx x)
         :y (if (m/> by y) by y)
         :util-fn curr-util
         :gp curr-gp
         :xs nxs
         :ys nys
         :util-best bx}))))

(defn bayesian-optimization
  "Maximizes an expensive to evaluate black box function with Bayesian optimization.

  A Gaussian process is fitted to the visited points. In every step the point which maximizes the utility function of the process is evaluated.

  Parameters:

  - `f` (function): the function to maximize. It receives the point as one sequence, or as separate arguments when `:vector-arg?` is `false`.
  - `opts` (map):
    - `:bounds` (required) - sequence of `[lo hi]` pairs, one for each dimension, validated for the `:optimizer`. Infinite bounds are not allowed.
    - `:vector-arg?` - how `f` receives the point, default: `true`. The utility function is always optimized with sequences, so any `:optimizer` works, including `:brent`.
    - `:warm-up` - number of scanned points used to find the maximum of the utility function, default: `1000` for every dimension.
    - `:init-points` - number of initial evaluations before the optimization starts, default: `3`. The points are selected with a jittered low discrepancy sequence generator (see [[fastmath.random/jittered-sequence-generator]]). A sequence of points can be given instead.
    - `:utility-function-type` - `:ucb` (default), `:ei` or `:poi`.
    - `:utility-param` - parameter of the utility function: `kappa` for `:ucb` (default: `2.576`), `xi` for `:ei` and `:poi` (default: `0.001`).
    - `:kernel` - kernel of the Gaussian process, a keyword or a kernel, default: `:matern-52`, see [[fastmath.kernel]].
    - `:kscale` - scaling factor of the kernel, default: `1.0`.
    - `:jitter` - jitter of the sequence generators, default: `0.25`.
    - `:noise` - noise (lambda) of the Gaussian process, default: `1.0e-8`.
    - `:normalize?` - normalize data in the Gaussian process, default: `true`.
    - `:optimizer` - method used to optimize the utility function, default: `:cmaes` for one dimension and `:lbfgsb` otherwise. A point found outside of the bounds by a method without constraints is moved to the bounds.
    - `:optimizer-params` - options of the optimizer. The scan options and `:vector-arg?`, `:stats?` and `:take-last-n` are set by the function.

  Returns a lazy sequence of consecutive steps. Every step is a map with:

  - `:x` - the best visited point, `:y` - its value,
  - `:xs` - all visited points, `:ys` - their values,
  - `:gp` - the current Gaussian process regression,
  - `:util-fn` - the current utility function,
  - `:util-best` - the maximum of the utility function, the point evaluated in the step.

  Throws `ex-info` for invalid bounds and an unknown `:utility-function-type` or `:optimizer`, and when the maximum of the utility function cannot be found (all runs of the optimizer failed).

  See also [[scan-and-maximize]], [[maximize]]."
  [f {:keys [^long warm-up init-points bounds utility-function-type utility-param kernel kscale
             jitter noise optimizer optimizer-params normalize? vector-arg?]
      :or {kscale 1.0
           kernel :matern-52
           init-points 3
           utility-function-type :ucb
           jitter 0.25
           normalize? true
           noise 1.0e-8}}]
  (let [;; the default optimizer depends on the number of dimensions, which any structural check of the bounds gives
        optimizer (or optimizer (if (m/one? (count (common/normalize-bounds :powell bounds nil))) :cmaes :lbfgsb))
        bounds (or (common/normalize-bounds optimizer bounds nil)
                   (throw (ex-info "Provide search bounds." {:optimizer optimizer :bounds bounds})))
        f (common/->vector-fn f (common/resolve-vector-arg? :bayesian-optimization vector-arg?))
        warm-up (or warm-up (m/* (count bounds) 1000))
        utility-param (double (or utility-param (if (#{:ei :poi} utility-function-type) 0.001 2.576)))
        kernel (if (keyword? kernel) (k/kernel kernel) kernel)
        [xs ys] (initial-values f init-points bounds jitter)
        [maxx maxy] (first (sort-by second m/> (map vector xs ys)))
        util-fn (utility-function utility-function-type utility-param)
        gp #(gp/gaussian-process %1 %2 {:normalize? normalize? :kernel kernel :kscale kscale :noise noise})
        step-fn (bayesian-step-fn f util-fn warm-up bounds gp jitter optimizer optimizer-params)]
    (rest (iterate step-fn {:x maxx
                            :y maxy
                            :xs xs
                            :ys ys}))))

;; linear optimization

(def ^:private relations
  {:<= Relationship/LEQ '<= Relationship/LEQ :leq Relationship/LEQ
   :>= Relationship/GEQ '>= Relationship/GEQ :geq Relationship/GEQ
   := Relationship/EQ '= Relationship/EQ :eq Relationship/EQ})

(defn- constraint-relation
  ^Relationship [relation]
  (or (relations relation)
      (common/throw-unknown :relation relation (vec (keys relations)))))

(defn- build-constraint
  ^LinearConstraint [[left relation right]]
  (let [relationship (constraint-relation relation)]
    (if (number? right)
      (LinearConstraint. (m/seq->double-array left) relationship (double right))
      (LinearConstraint. (m/seq->double-array (butlast left))
                         (double (last left))
                         relationship
                         (m/seq->double-array (butlast right))
                         (double (last right))))))

(defn- pivot-rule
  ^PivotSelectionRule [rule]
  (case rule
    :dantzig PivotSelectionRule/DANTZIG
    :bland PivotSelectionRule/BLAND
    (common/throw-unknown :rule rule #{:dantzig :bland})))

(defn linear-optimization
  "Solves a linear programming problem using the simplex method.

  The objective is a linear function of any number of variables subject to a set of linear equality or inequality constraints. This is a distinct, specialized solver and does not go through [[minimize]]/[[maximize]]/[[minimizer]]/[[maximizer]]; use [[bayesian-optimization]] or the general optimizers for nonlinear problems.

  Parameters:

  - `target` (vector of numbers): coefficients of the objective function with the constant term as the last value. `[a1 a2 a3 ... c]` represents `f(x1,x2,x3,...) = a1*x1 + a2*x2 + a3*x3 + ... + c`.
  - `constraints` (flat sequence): a concatenation of triplets `left R right`, each of one of the following forms:
      - `[a1 a2 a3 ...] R n` — means `a1*x1 + a2*x2 + a3*x3 + ... R n`, where `n` is a number.
      - `[a1 a2 a3 ... ca] R [b1 b2 b3 ... cb]` — means `a1*x1 + a2*x2 + a3*x3 + ... + ca R b1*x1 + b2*x2 + b3*x3 + ... + cb`.
      - `R` is the relationship: `<=`, `>=` or `=` (as symbols or keywords), or `:leq`, `:geq` or `:eq`. Nothing else is accepted.
  - `options` (optional map):
      - `:goal` — `:minimize` (default) or `:maximize`.
      - `:rule` — pivot selection rule, `:dantzig` (default) or `:bland`.
      - `:max-iters` — maximum number of iterations, default `10000`. Exceeding it throws an exception.
      - `:non-negative?` — when `true`, restrict all variables to non-negative values, default `false`.
      - `:epsilon` — convergence tolerance, default `1.0e-6`.
      - `:max-ulps` — allowed floating point comparison tolerance expressed in ulps, default `10`.
      - `:cut-off` — pivot elements smaller than this value are treated as zero, default `1.0e-10`.
      - `:stats?` — when `true`, return a map with the iteration count, default `false`.

  Returns a pair `[point value]`, where `point` is a vector of optimal variable values and `value` is the optimal objective function value. When `:stats?` is set to `true`, returns instead a map with `:point`, `:value`, `:evaluations` (always `0`) and `:iterations`.

  Every three consecutive values of `constraints` are treated as one triplet, so the collection must contain a multiple of three elements matching the pattern above. Throws `ex-info` otherwise, and for an unknown relationship, `:goal` or `:rule`. Infeasible and unbounded problems throw the exceptions of Apache Commons Math.

  For example `(linear-optimization [-1 4 0] [[-3 1] :<= 6 [-1 -2] :>= -4 [0 1] :>= -3])` returns `[[9.999999999999995 -3.0] -21.999999999999993]`."
  ([target constraints] (linear-optimization target constraints {}))
  ([target constraints {:keys [goal ^double epsilon ^int max-ulps ^double cut-off
                               rule non-negative? max-iters stats?]
                        :or {goal :minimize epsilon 1.0e-6 max-ulps 10 cut-off 1.0e-10
                             rule :dantzig non-negative? false}}]
   (when-not (zero? (rem (count constraints) 3))
     (throw (ex-info "Constraints should be a sequence of triplets: left, relation, right" {:count (count constraints)})))
   (let [goal (common/goal-type goal)
         rule (pivot-rule rule)
         max-iter (common/max-iter max-iters)
         non-negative? (NonNegativeConstraint. non-negative?)
         target (LinearObjectiveFunction. (m/seq->double-array (butlast target))
                                          (double (last target)))
         constraints (->> constraints
                          (partition 3)
                          ^java.util.Collection (map build-constraint)
                          (LinearConstraintSet.))
         ^BaseOptimizer solver (SimplexSolver. epsilon max-ulps cut-off)]
     (common/multivariate-optimize solver [goal rule max-iter non-negative? target constraints] stats?))))


#_(let [f (fn [^double x ^double y] (inc (- (- (* x x)) (m/sq (dec y)))))
        bounds [[-4 4] [-3 3]]
        bo (bayesian-optimization f {:bounds bounds
                                     :utility-function-type :ei
                                     :utility-param 0.1
                                     :optimizer :powell})]
    (println (f 0 1))
    (last (take 30 (map (juxt :x :y) bo)))    )

#_(let [f (fn [^double x] (- (+ (/ (m/sin (* 10 m/PI x)) (+ x x)) (m/pow (dec x) 4))))
        bounds [[-0.2 0.2]]
        bo (bayesian-optimization f {:bounds bounds
                                     :utility-function-type :ucb
                                     ;; :utility-param 0
                                     :optimizer :lbfgsb})]
    (take 3 (drop 30 (map (juxt :x :y) bo))))
;; => ([(0.14580050823621077) 2.867142408546129] [(0.14564645826782044) 2.868128502656343] [(0.14550225115172954) 2.868981775657352])
;; => ([(0.14973586005947528) 2.8164430833342324] [(0.14973586005947528) 2.8164430833342324] [(0.14973586005947528) 2.8164430833342324])

;; tests

#_(defn target-1d ^double [^double x]
    (+ (/ (m/sin (* 10 m/PI x)) (+ x x)) (m/pow (dec x) 4)))

#_(defn target-2d-schwefel
    ^double [^double x ^double y]
    (- 418.9829
       (+ (* x (m/sin (m/sqrt (m/abs x))))
          (* y (m/sin (m/sqrt (m/abs y)))))))

#_(defn target-2d-booth
    ^double [^double x ^double y]
    (+ (m/sq (+ x y y -7))
       (m/sq (+ x x y -5))))

#_(time (let [f (minimizer :lbfgsb target-2d-booth {:bounds [[-10 10]
                                                             [-10 10]]
                                                    :tolerance 1.0e-10})]
          (f (v/generate-vec2 #(r/drand -9 9)))))

#_(scan-and-minimize :lbfgsb target-2d-schwefel {:bounds [[-500 500]
                                                          [-500 500]] :N 1000 :bounded? true})

#_(scan-and-minimize :gradient target-1d {:bounds [[-0.5 2.5]] :gradient-step 0.1})

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#_(do

    (defn bfn6 ^double [^double x ^double y] (+ (* 100 (m/sqrt (m/abs (- y (* 0.01 x x)))))
                                             (* 0.01 (m/abs (+ 10.0 x)))))

    (defn h ^double [^double x ^double y] (+ (m/sq (+ (* x x) y -11))
                                          (m/sq (+ x (* y y) -7))))

    (defn d5 [a b c d e] (reduce #(+ ^double %1 (m/sq %2)) 0.0 [a b c d e]))

    (time (scan-and-maximize :bobyqa bfn6 {:bounds [[-15 -3] [15 3]] :N 100 :n 0.2}))

    (minimize :powell bfn6 {:bounds [[-15 -3] [15 3]] :initial [-13.217532309719662 0.9280007414397415]})

    (time (scan-and-optimize :powell #(m/cos %) {:bounds [-3 3] :initial -2 :N 100 :n 0.2 :goal :maximize})))


;; => [1.9210981963566007 -5.751481824637489 0.3304425131054902]

#_(defn rosenbrock
    [& vs]
    (reduce (fn [^double fx [^double xi ^double xi+1]]
              (let [t1 (- 1.0 xi)
                    t2 (* 10.0 (- xi+1 (* xi xi)))]
                (+ fx (* t1 t1) (* t2 t2)))) 0.0 (partition 2 1 vs)))


#_(minimize :lbfgsb rosenbrock {:bounds (repeat 20 [-5 10])
                                :init [2 -4 2 4 -2] :m 50 :N 10 :n 1
                                :max-iters 1000})



#_(comment (defn hump
             [^double x ^double y]
             (let [x2 (* x x)
                   x4 (* x2 x2)]
               (+ (* 2.0 x2)
                  (* -1.05 x4)
                  (* x4 x2 m/SIXTH)
                  (* x y)
                  (* y y))))

           (defn hump-grad
             [^double x ^double y]
             (let [x2 (* x x)
                   x4 (* x2 x2)]
               [(+ (* 4.0 x)
                   (* -4.2 x2 x)
                   (* x4 x)
                   y)
                (+ x (* 2.0 y))]))

           (minimize :lbfgsb hump {:bounds [[-5 5.1] [-5 5.1]]
                                   :gradient-f hump-grad :N 1000
                                   :weak-wolfe? false
                                   :stats? true}))

