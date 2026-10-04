(ns fastmath.optimization.acm
  "Optimizers from the Apache Commons Math library.

  Every optimizer is a function of the objective and an options map, `(optimizer f opts)`, and creates all its Apache Commons Math objects from scratch on every call. Nothing is shared between calls, so the same function can be called from many threads.

  Optimizers:

  - [[brent]] - one dimensional, searches an interval.
  - [[bobyqa]] - derivative free, box constrained, two or more dimensions.
  - [[cmaes]] - derivative free evolutionary strategy, box constrained.
  - [[nelder-mead]] and [[multidirectional-simplex]] - derivative free simplex methods, unconstrained. Bounds only set the size of the initial simplex.
  - [[powell]] - derivative free direction set method, unconstrained.
  - [[non-linear-gradient]] - conjugate gradient method, unconstrained, with the gradient given by the user or approximated numerically.

  Common options of all optimizers:

  - `:bounds` - sequence of `[lo hi]` pairs, one for each dimension. What is required and what is allowed depends on the optimizer, see [[fastmath.optimization.common/normalize-bounds]].
  - `:initial` - the initial point, default: the middle of the bounds.
  - `:goal` - `:minimize` (default) or `:maximize`.
  - `:vector-arg?` - when `true`, the objective receives the point as one sequence (a `double[]`), when `false` as separate arguments. Default: `true`, only [[brent]] defaults to `false` and receives a number.
  - `:max-evals`, `:max-iters` - limits on the number of evaluations and iterations, default: `10000` each. Exceeding a limit throws an exception. The exceptions are [[cmaes]], which stops after `:max-iters` iterations and returns the best point found, and ignores `:max-evals`, and [[bobyqa]], which does not count iterations.
  - `:stats?` - return a map with additional information, default: `false`.

  The result is `[point value]`, where `point` is a vector (a number for [[brent]] with `:vector-arg?` set to `false`). With `:stats?` it is a map with `:point`, `:value`, `:evaluations` and `:iterations`. The value is always the value of the objective, also when maximizing.

  Invalid options and bounds throw `ex-info`. Errors reported by Apache Commons Math (limits exceeded, convergence failures) are thrown as they are."
  (:require [fastmath.core :as m]
            [fastmath.calculus.finite :as finite]
            [fastmath.optimization.common :as common])
  (:import [org.apache.commons.math3.optim.univariate BrentOptimizer UnivariatePointValuePair SearchInterval UnivariateObjectiveFunction BracketFinder]
           [org.apache.commons.math3.optim.nonlinear.scalar GoalType ObjectiveFunctionGradient ObjectiveFunction]
           [org.apache.commons.math3.optim SimpleBounds InitialGuess SimpleValueChecker]
           [org.apache.commons.math3.optim.nonlinear.scalar.noderiv BOBYQAOptimizer CMAESOptimizer CMAESOptimizer$PopulationSize CMAESOptimizer$Sigma NelderMeadSimplex MultiDirectionalSimplex SimplexOptimizer PowellOptimizer]
           [org.apache.commons.math3.optim.nonlinear.scalar.gradient NonLinearConjugateGradientOptimizer NonLinearConjugateGradientOptimizer$Formula NonLinearConjugateGradientOptimizer$IdentityPreconditioner Preconditioner]
           [org.apache.commons.math3.analysis UnivariateFunction MultivariateFunction MultivariateVectorFunction]
           [org.apache.commons.math3.random JDKRandomGenerator]
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; univariate

(defn- univariate-function
  ^UnivariateFunction [f vector-arg?]
  (if vector-arg?
    (reify UnivariateFunction (value [_ x] (double (f [x]))))
    (reify UnivariateFunction (value [_ x] (double (f x))))))

(defn- bracket
  "Moves the interval `[lo hi]` to a bracket of an extremum found by `BracketFinder`."
  [find-bracket ^UnivariateFunction uf ^GoalType goal lo hi]
  (let [{:keys [^double grow-limit ^long max-evals]
         :or {grow-limit 100.0 max-evals 500}} (if (map? find-bracket) find-bracket {})
        ^BracketFinder bf (BracketFinder. grow-limit (int max-evals))]
    (.search bf uf goal (double lo) (double hi))
    [(.getLo bf) (.getHi bf)]))

(defn brent
  "Finds a minimum or a maximum of a one dimensional function in an interval with Brent's method.

  The method is local: it finds one extremum of the interval, which is not necessarily the global one.

  Parameters:

  - `f` (function): the objective. It receives a number, or a sequence with one number when `:vector-arg?` is `true`.
  - `opts` (map):
    - `:bounds` (required) - the interval as `[[lo hi]]` or `[lo hi]`, finite with `lo < hi`.
    - `:initial` - the starting point, a number or a sequence with one number, inside of the interval (after `:find-bracket` - inside of the bracket).
    - `:goal` - `:minimize` (default) or `:maximize`.
    - `:vector-arg?` - default: `false`.
    - `:rel`, `:abs` - relative and absolute accuracy, default: `1.0e-8` and `1.0e-10`.
    - `:find-bracket` - when truthy, the interval is first moved to a bracket of an extremum found by Apache Commons Math `BracketFinder`. A map with `:grow-limit` (default: `100.0`) and `:max-evals` (default: `500`) sets its parameters.
    - `:max-evals`, `:max-iters` - limits, default: `10000` each.
    - `:stats?` - return a map with additional information, default: `false`.

  Returns `[point value]` where `point` is a number (a vector with one number when `:vector-arg?` is `true`). With `:stats?` returns a map with `:point` (in the same form), `:value`, `:evaluations` and `:iterations`.

  Throws `ex-info` for invalid bounds (missing, more than one interval, infinite, empty) and an Apache Commons Math exception when a limit is exceeded or `:initial` is outside of the interval.

  See also [[bobyqa]], [[powell]], [[fastmath.optimization/minimize]]."
  [f {:keys [bounds initial max-evals max-iters ^double rel ^double abs find-bracket vector-arg? goal stats?]
      :or {rel 1.0e-8 abs 1.0e-10}}]
  (let [[[lo hi]] (common/normalize-bounds :brent bounds initial)
        vector-arg? (common/resolve-vector-arg? :brent vector-arg?)
        goal (common/goal-type goal)
        uf (univariate-function f vector-arg?)
        [lo hi] (if find-bracket (bracket find-bracket uf goal lo hi) [lo hi])
        start (when (some? initial) (if (number? initial) initial (first initial)))
        interval (if (some? start)
                   (SearchInterval. (double lo) (double hi) (double start))
                   (SearchInterval. (double lo) (double hi)))
        ^BrentOptimizer bo (BrentOptimizer. rel abs)
        ^UnivariatePointValuePair res (.optimize bo (common/optimization-data [goal
                                                                              (UnivariateObjectiveFunction. uf)
                                                                              interval
                                                                              (common/max-eval max-evals)
                                                                              (common/max-iter max-iters)]))
        x (.getPoint res)
        point (if vector-arg? [x] x)]
    (if stats?
      {:point point
       :value (.getValue res)
       :evaluations (.getEvaluations bo)
       :iterations (.getIterations bo)}
      [point (.getValue res)])))

;; multivariate

(defn- multivariate-function
  ^MultivariateFunction [f]
  (reify MultivariateFunction (value [_ x] (double (f x)))))

(defn- objective-function
  ^ObjectiveFunction [f]
  (ObjectiveFunction. (multivariate-function f)))

(defn- simple-bounds
  ^SimpleBounds [bounds]
  (SimpleBounds. (double-array (map first bounds))
                 (double-array (map second bounds))))

(defn- bounds->steps
  ^doubles [bounds length]
  (let [length (double length)]
    (double-array (map (fn [[^double lo ^double hi]] (m/* length (m/- hi lo))) bounds))))

(defn- multivariate-base
  "Validates the options common to all multivariate optimizers.

  Returns a map with `:f` (the objective as a function of one sequence), `:bounds` (normalized or `nil`), `:goal` (Apache Commons Math `GoalType`) and `:data` (the evaluation and iteration limits and the initial guess, as optimization data).

  Throws `ex-info` when neither bounds nor the initial point are given."
  [method f {:keys [bounds initial goal max-evals max-iters vector-arg?]}]
  (let [initial (if (number? initial) [initial] initial)
        bounds (common/normalize-bounds method bounds initial)
        initial (cond
                  initial (m/seq->double-array initial)
                  bounds (common/bounds-midpoint bounds)
                  :else (throw (ex-info "Bounds or an initial point should be provided" {:method method :bounds bounds :initial initial})))]
    {:f (common/->vector-fn f (common/resolve-vector-arg? method vector-arg?))
     :bounds bounds
     :goal (common/goal-type goal)
     :data [(common/max-eval max-evals) (common/max-iter max-iters) (InitialGuess. initial)]}))

(defn- bobyqa-radius
  ^double [bounds]
  (m/* 0.5 (double (reduce m/max (map (fn [[^double lo ^double hi]] (m/- hi lo)) bounds)))))

(defn bobyqa
  "Minimizes or maximizes a function within box constraints with the BOBYQA method (bound optimization by quadratic approximation), which does not need derivatives.

  Parameters:

  - `f` (function): the objective. It receives the point as one sequence (a `double[]`), or as separate arguments when `:vector-arg?` is `false`.
  - `opts` (map):
    - `:bounds` (required) - sequence of finite `[lo hi]` pairs, at least two dimensions.
    - `:initial` - the initial point inside of the bounds, default: the middle of the bounds.
    - `:goal` - `:minimize` (default) or `:maximize`.
    - `:vector-arg?` - default: `true`.
    - `:number-of-points` - number of interpolation points, from `n+2` to `(n+1)(n+2)/2` for `n` dimensions, default: `2n+1`.
    - `:initial-radius` - initial trust region radius: `:inferred` (default, half of the widest bound range), `:default` (the Apache Commons Math default) or a number.
    - `:stopping-radius` - the trust region radius at which the method stops, default: the Apache Commons Math default (`1.0e-8`).
    - `:max-evals` - limit on the number of evaluations, default: `10000`. `:max-iters` is accepted but BOBYQA does not count iterations.
    - `:stats?` - return a map with additional information, default: `false`.

  Returns `[point value]`. With `:stats?` returns a map with `:point`, `:value`, `:evaluations` and `:iterations` (always `0`).

  Throws `ex-info` for invalid bounds or fewer than two dimensions and an Apache Commons Math exception when the number of points is out of range or the evaluation limit is exceeded.

  See also [[cmaes]], [[powell]], [[fastmath.optimization/minimize]]."
  [f {:keys [number-of-points initial-radius ^double stopping-radius]
      :or {initial-radius :inferred
           stopping-radius BOBYQAOptimizer/DEFAULT_STOPPING_RADIUS}
      :as opts}]
  (let [{:keys [bounds goal data] f :f} (multivariate-base :bobyqa f opts)
        n (count bounds)
        _ (when (m/< n 2) (throw (ex-info "Number of dimensions should be equal or greater than 2" {:n n :bounds bounds})))
        radius (condp = initial-radius
                 :default BOBYQAOptimizer/DEFAULT_INITIAL_RADIUS
                 :inferred (bobyqa-radius bounds)
                 (double initial-radius))
        ^BOBYQAOptimizer bo (BOBYQAOptimizer. (int (or number-of-points (m/inc (m/* 2 n)))) radius stopping-radius)]
    (common/multivariate-optimize bo (into [goal (objective-function f) (simple-bounds bounds)] data) (:stats? opts))))

(defn cmaes
  "Minimizes or maximizes a function within box constraints with the covariance matrix adaptation evolution strategy (CMA-ES), which does not need derivatives.

  The method is stochastic: use `:rng` with a seeded generator for reproducible results.

  Parameters:

  - `f` (function): the objective. It receives the point as one sequence (a `double[]`), or as separate arguments when `:vector-arg?` is `false`.
  - `opts` (map):
    - `:bounds` (required) - sequence of finite `[lo hi]` pairs.
    - `:initial` - the initial point inside of the bounds, default: the middle of the bounds.
    - `:goal` - `:minimize` (default) or `:maximize`.
    - `:vector-arg?` - default: `true`.
    - `:sigma` - initial standard deviation as a fraction of the range of every dimension, default: `0.2`.
    - `:population-size` - number of candidates in every generation, default: `ceil(4 + 3 ln n)` for `n` dimensions.
    - `:stop-fitness` - stop when the objective value is below this value, default: `1.0e-10`.
    - `:active-cma?` - use the active covariance matrix update, default: `true`.
    - `:diagonal-only` - number of initial iterations with a diagonal covariance matrix only, default: `0`.
    - `:check-feasible-count` - number of resampling attempts for a candidate outside of the bounds, default: `0`.
    - `:rng` - random number generator (an Apache Commons Math `RandomGenerator`), default: a new unseeded `JDKRandomGenerator` for every call. A generator given here is shared by all runs which use it, so use it from one thread only unless it is thread safe.
    - `:rel`, `:abs` - relative and absolute tolerance of the objective value change, default: `1.0e-10` each.
    - `:max-iters` - maximum number of iterations, default: `10000`. Reaching it is not an error: the best point found is returned. `:max-evals` is ignored by this method.
    - `:stats?` - return a map with additional information, default: `false`.

  Returns `[point value]`. With `:stats?` returns a map with `:point`, `:value`, `:evaluations` and `:iterations`.

  Throws `ex-info` for invalid bounds and an Apache Commons Math exception when the method does not converge.

  See also [[bobyqa]], [[powell]], [[fastmath.optimization/minimize]]."
  [f {:keys [stop-fitness active-cma? diagonal-only check-feasible-count rng ^double rel ^double abs population-size sigma max-iters stats?]
      :or {stop-fitness 1.0e-10 active-cma? true diagonal-only 0 check-feasible-count 0
           rel 1.0e-10 abs 1.0e-10 sigma 0.2}
      :as opts}]
  (let [{:keys [bounds goal data] f :f} (multivariate-base :cmaes f opts)
        population-size (CMAESOptimizer$PopulationSize. (int (or population-size (m/ceil (m/+ 4.0 (m/* 3.0 (m/log (count bounds))))))))
        sigma (CMAESOptimizer$Sigma. (bounds->steps bounds sigma))
        ^CMAESOptimizer co (CMAESOptimizer. (int (or max-iters 10000)) (double stop-fitness) (boolean active-cma?)
                                            (int diagonal-only) (int check-feasible-count)
                                            (or rng (JDKRandomGenerator.)) false (SimpleValueChecker. rel abs))]
    (common/multivariate-optimize co (into [goal (objective-function f) (simple-bounds bounds) population-size sigma] data) stats?)))

(defn- simplex-steps
  "Edge lengths of the initial simplex: a fraction of the bounds ranges, or a constant when there are no bounds."
  ^doubles [bounds ^doubles initial length]
  (if bounds
    (bounds->steps bounds (or length 0.4))
    (double-array (repeat (alength initial) (double (or length 1.0))))))

(defn- simplex-optimize
  [method f opts make-simplex]
  (let [{:keys [bounds goal data] f :f} (multivariate-base method f opts)
        ^InitialGuess guess (peek data)
        steps (simplex-steps bounds (.getInitialGuess guess) (:length opts))
        {:keys [rel abs] :or {rel 1.0e-10 abs 1.0e-10}} opts
        ^SimplexOptimizer so (SimplexOptimizer. (double rel) (double abs))]
    (common/multivariate-optimize so (into [goal (objective-function f) (make-simplex steps)] data) (:stats? opts))))

(defn nelder-mead
  "Minimizes or maximizes a function with the Nelder-Mead simplex method, which does not need derivatives and does not use constraints.

  The bounds are not constraints: they only set the size of the initial simplex, so the optimizer may leave them.

  Parameters:

  - `f` (function): the objective. It receives the point as one sequence (a `double[]`), or as separate arguments when `:vector-arg?` is `false`.
  - `opts` (map):
    - `:bounds` - sequence of finite `[lo hi]` pairs with `lo < hi`. Optional when `:initial` is given.
    - `:initial` - the initial point, default: the middle of the bounds.
    - `:goal` - `:minimize` (default) or `:maximize`.
    - `:vector-arg?` - default: `true`.
    - `:length` - edge length of the initial simplex: a fraction of the range of every dimension when `:bounds` are given (default: `0.4`), an absolute length otherwise (default: `1.0`).
    - `:rho`, `:khi`, `:gamma`, `:sigma` - reflection, expansion, contraction and shrinkage coefficients, default: `1.0`, `2.0`, `0.5` and `0.5`.
    - `:rel`, `:abs` - relative and absolute tolerance of the objective value change, default: `1.0e-10` each.
    - `:max-evals`, `:max-iters` - limits, default: `10000` each.
    - `:stats?` - return a map with additional information, default: `false`.

  Returns `[point value]`. With `:stats?` returns a map with `:point`, `:value`, `:evaluations` and `:iterations`.

  Throws `ex-info` when neither `:bounds` nor `:initial` are given or the bounds are invalid, and an Apache Commons Math exception when a limit is exceeded.

  See also [[multidirectional-simplex]], [[powell]], [[fastmath.optimization/minimize]]."
  [f {:keys [rho khi gamma sigma]
      :or {rho 1.0 khi 2.0 gamma 0.5 sigma 0.5}
      :as opts}]
  (simplex-optimize :nelder-mead f opts
                    (fn [steps] (NelderMeadSimplex. ^doubles steps (double rho) (double khi) (double gamma) (double sigma)))))

(defn multidirectional-simplex
  "Minimizes or maximizes a function with the multi-directional simplex method of Torczon, which does not need derivatives and does not use constraints.

  The bounds are not constraints: they only set the size of the initial simplex, so the optimizer may leave them.

  Parameters:

  - `f` (function): the objective. It receives the point as one sequence (a `double[]`), or as separate arguments when `:vector-arg?` is `false`.
  - `opts` (map):
    - `:bounds` - sequence of finite `[lo hi]` pairs with `lo < hi`. Optional when `:initial` is given.
    - `:initial` - the initial point, default: the middle of the bounds.
    - `:goal` - `:minimize` (default) or `:maximize`.
    - `:vector-arg?` - default: `true`.
    - `:length` - edge length of the initial simplex: a fraction of the range of every dimension when `:bounds` are given (default: `0.4`), an absolute length otherwise (default: `1.0`).
    - `:khi`, `:gamma` - expansion and contraction coefficients, default: `2.0` and `0.5`.
    - `:rel`, `:abs` - relative and absolute tolerance of the objective value change, default: `1.0e-10` each.
    - `:max-evals`, `:max-iters` - limits, default: `10000` each.
    - `:stats?` - return a map with additional information, default: `false`.

  Returns `[point value]`. With `:stats?` returns a map with `:point`, `:value`, `:evaluations` and `:iterations`.

  Throws `ex-info` when neither `:bounds` nor `:initial` are given or the bounds are invalid, and an Apache Commons Math exception when a limit is exceeded.

  See also [[nelder-mead]], [[powell]], [[fastmath.optimization/minimize]]."
  [f {:keys [khi gamma]
      :or {khi 2.0 gamma 0.5}
      :as opts}]
  (simplex-optimize :multidirectional-simplex f opts
                    (fn [steps] (MultiDirectionalSimplex. ^doubles steps (double khi) (double gamma)))))

(defn- negate-value
  "Negates the value in the result of an optimizer, `[point value]` or a map."
  [res stats?]
  (if stats?
    (update res :value m/-)
    (update res 1 m/-)))

(defn powell
  "Minimizes or maximizes a function with Powell's conjugate direction method, which does not need derivatives and does not use constraints.

  Maximization is done by minimizing the negated function: the method of Apache Commons Math stops after the first iteration and returns a wrong result when it is asked to maximize.

  Parameters:

  - `f` (function): the objective. It receives the point as one sequence (a `double[]`), or as separate arguments when `:vector-arg?` is `false`.
  - `opts` (map):
    - `:initial` - the initial point. Required when `:bounds` are not given.
    - `:bounds` - optional, only used to find the middle as the default initial point; they are not constraints.
    - `:goal` - `:minimize` (default) or `:maximize`.
    - `:vector-arg?` - default: `true`.
    - `:rel`, `:abs` - relative and absolute tolerance of the objective value change, default: `1.0e-10` each.
    - `:line-rel`, `:line-abs` - relative and absolute tolerance of the line search, default: the square roots of `:rel` and `:abs`.
    - `:max-evals`, `:max-iters` - limits, default: `10000` each.
    - `:stats?` - return a map with additional information, default: `false`.

  Returns `[point value]`, where `value` is the value of `f`. With `:stats?` returns a map with `:point`, `:value`, `:evaluations` and `:iterations`.

  Throws `ex-info` when neither `:bounds` nor `:initial` are given or the bounds are invalid, and an Apache Commons Math exception when a limit is exceeded.

  See also [[nelder-mead]], [[bobyqa]], [[fastmath.optimization/minimize]]."
  [f {:keys [^double rel ^double abs line-rel line-abs stats?]
      :or {rel 1.0e-10 abs 1.0e-10}
      :as opts}]
  (let [{:keys [goal data] f :f} (multivariate-base :powell f opts)
        maximize? (= goal GoalType/MAXIMIZE)
        objective (if maximize? (fn [x] (m/- (double (f x)))) f)
        ^PowellOptimizer po (PowellOptimizer. rel abs (double (or line-rel (m/sqrt rel))) (double (or line-abs (m/sqrt abs))))
        res (common/multivariate-optimize po (into [GoalType/MINIMIZE (objective-function objective)] data) stats?)]
    (if maximize?
      (negate-value res stats?)
      res)))

;; gradient

(defn- gradient-function
  "Wraps the gradient (a function of one sequence returning a sequence of numbers) for Apache Commons Math."
  ^MultivariateVectorFunction [gradient]
  (reify MultivariateVectorFunction
    (value [_ x]
      (let [^doubles res (or (m/seq->double-array (gradient x)) (double-array 0))
            n (alength ^doubles x)]
        (when-not (m/== n (alength res))
          (throw (ex-info "Gradient has wrong length" {:expected n :actual (alength res)})))
        res))))

(defn- hessian-preconditioner
  [f h]
  (let [hd (finite/hessian-diagonal f {:h h})]
    (reify Preconditioner
      (precondition [_ point r]
        (let [^doubles diagonal (hd point)
              ^doubles nr (double-array r)]
          (dotimes [i (alength nr)]
            (let [d (Array/aget diagonal i)]
              (when (m/> d 1.0e-6)
                (Array/aset nr i (m// (Array/aget r i) d)))))
          nr)))))

(defn non-linear-gradient
  "Minimizes or maximizes a function with the non-linear conjugate gradient method, which does not use constraints.

  The gradient is given by the user or approximated with finite differences. The method is also available under the name `:gradient` in [[fastmath.optimization/minimize]].

  Parameters:

  - `f` (function): the objective. It receives the point as one sequence (a `double[]`), or as separate arguments when `:vector-arg?` is `false`.
  - `opts` (map):
    - `:initial` - the initial point. Required when `:bounds` are not given.
    - `:bounds` - optional, only used to find the middle as the default initial point; they are not constraints.
    - `:goal` - `:minimize` (default) or `:maximize`.
    - `:vector-arg?` - default: `true`.
    - `:gradient` - function of the point, always one sequence, returning the gradient of `f` as a sequence of numbers. It is the gradient of `f` itself, also when maximizing. Default: the gradient is approximated with finite differences.
    - `:gradient-h`, `:gradient-acc` - step (default: `1.0e-6`, must be positive) and accuracy order (`2` or `4`, default: `2`) of the finite differences.
    - `:formula` - update formula of the conjugate direction: `:polak-ribiere` (default) or `:fletcher-reeves`.
    - `:preconditioner` - `:identity` (default) or `:hessian`, which divides the gradient by the diagonal of the Hessian approximated with finite differences with step `:hessian-h` (default: `5.0e-3`).
    - `:bracketing-range` - initial bracketing range of the line search, default: `1.0e-10`.
    - `:rel`, `:abs` - relative and absolute tolerance of the objective value change, default: `1.0e-10` each.
    - `:line-rel`, `:line-abs` - relative and absolute tolerance of the line search, default: the square roots of `:rel` and `:abs`.
    - `:max-evals`, `:max-iters` - limits, default: `10000` each.
    - `:stats?` - return a map with additional information, default: `false`.

  Returns `[point value]`. With `:stats?` returns a map with `:point`, `:value`, `:evaluations` and `:iterations`.

  Throws `ex-info` when neither `:bounds` nor `:initial` are given, when the bounds, `:gradient-h`, `:gradient-acc`, `:formula` or `:preconditioner` are invalid, and when the gradient has a wrong length. Apache Commons Math exceptions are thrown when a limit is exceeded.

  See also [[fastmath.optimization.lbfgsb/lbfgsb]], [[powell]], [[fastmath.optimization/minimize]]."
  [f {:keys [gradient gradient-h gradient-acc ^double rel ^double abs line-rel line-abs formula ^double bracketing-range
             preconditioner hessian-h stats?]
      :or {gradient-h 1.0e-6 gradient-acc 2 rel 1.0e-10 abs 1.0e-10 formula :polak-ribiere bracketing-range 1.0e-10
           preconditioner :identity hessian-h 5.0e-3}
      :as opts}]
  (let [{:keys [goal data] f :f} (multivariate-base :non-linear-gradient f opts)
        _ (when-not (m/pos? (double gradient-h)) (throw (ex-info "gradient-h must be positive" {:gradient-h gradient-h})))
        _ (when-not (contains? #{2 4} gradient-acc) (common/throw-unknown :gradient-acc gradient-acc #{2 4}))
        gradient (or gradient (finite/gradient f {:h gradient-h :acc gradient-acc}))
        formula (case formula
                  :polak-ribiere NonLinearConjugateGradientOptimizer$Formula/POLAK_RIBIERE
                  :fletcher-reeves NonLinearConjugateGradientOptimizer$Formula/FLETCHER_REEVES
                  (common/throw-unknown :formula formula #{:polak-ribiere :fletcher-reeves}))
        ^Preconditioner preconditioner (case preconditioner
                                         :identity (NonLinearConjugateGradientOptimizer$IdentityPreconditioner.)
                                         :hessian (hessian-preconditioner f hessian-h)
                                         (common/throw-unknown :preconditioner preconditioner #{:identity :hessian}))
        ^NonLinearConjugateGradientOptimizer nlcgo (NonLinearConjugateGradientOptimizer.
                                                    formula (SimpleValueChecker. rel abs)
                                                    (double (or line-rel (m/sqrt rel))) (double (or line-abs (m/sqrt abs)))
                                                    bracketing-range preconditioner)]
    (common/multivariate-optimize nlcgo (into [goal (objective-function f) (ObjectiveFunctionGradient. (gradient-function gradient))] data) stats?)))
