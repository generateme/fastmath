(ns fastmath.calculus
  "Integration and differentiation of functions, divided differences and integration of sampled data.

  Integration of functions:

  * [[integrate]] - univariate functions: adaptive Gauss-Kronrod (default), Gauss-Legendre, Romberg, Simpson, midpoint and trapezoid rules; infinite bounds and break points are supported.
  * [[cubature]] - multivariate functions: h-adaptive Genz-Malik cubature.
  * [[vegas]] - multivariate functions: VEGAS and VEGAS+ Monte Carlo integration.
  * [[integrate-samples]] - sampled data: trapezoid and Simpson rules, also cumulative.

  Improper integrals (infinite bounds) are mapped to finite intervals by a substitution. Reversed bounds change the sign of the integral.

  Derivatives (finite differences):

  * [[derivative]] - derivatives of any order and any order of accuracy, forward, backward and central, optionally with Richardson extrapolation ([[extrapolate]]).
  * [[gradient]], [[jacobian]], [[hessian]], [[hessian-diagonal]], [[laplacian]] - derivatives of multivariate functions.

  Divided differences, including repeated nodes (Hermite): [[divided]], [[divided-all]]."
  (:require [fastmath.calculus.common :as cc]
            [fastmath.calculus.vegas :as vegas]
            [fastmath.calculus.cubature :as cubature]
            [fastmath.calculus.quadrature :as quadrature]
            [fastmath.calculus.finite :as finite]
            [fastmath.calculus.divided :as divided]
            [fastmath.calculus.samples :as samples])
  (:import [org.apache.commons.math3.analysis UnivariateFunction]
           [org.apache.commons.math3.analysis.integration UnivariateIntegrator
            IterativeLegendreGaussIntegrator MidPointIntegrator
            RombergIntegrator SimpsonIntegrator TrapezoidIntegrator
            BaseAbstractUnivariateIntegrator]
           [org.apache.commons.math3.exception MaxCountExceededException TooManyEvaluationsException]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(defn vegas
  "VEGAS and VEGAS+ Monte Carlo integration of a multivariate function on a box.

  The algorithm adapts an importance sampling grid in each dimension (VEGAS) and, with stratified sampling, the number of samples in each hypercube (VEGAS+). The results of the iterations are combined with weights equal to the inverse of their variances. Improper integrals (infinite bounds) are handled by a substitution. The estimate is random: seed `:rng` to reproduce it.

  Parameters:

  - `f` (function): the integrand, a function of a point returning a double. It receives a `double[]` when all the bounds are finite, which is reused for the next call, so it is valid only during the call. When any bound is infinite it receives a vector. Use `nth`, `aget` or destructuring to read the coordinates.
  - `lower`, `upper` (sequences of numbers): bounds of the box, one per dimension, of the same length (one or more dimensions). A dimension with `lower` above `upper` is reversed: the result changes the sign for each reversed dimension. When `lower` equals `upper` in any dimension the integral is `0.0`.
  - `options` (map, optional):
    - `:max-iters` - maximum number of iterations, not counting the warmup, an integer >= 1 (default: `10`).
    - `:nevals` - number of evaluations per iteration, an integer >= 2 (default: `10000`).
    - `:nintervals` - number of grid intervals per dimension, an integer >= 1 (default: `1000`). The actual size is the smallest multiple of the numbers of strata which is not below it.
    - `:nstrats` - number of strata per dimension, a number or a sequence cycled over the dimensions (default: calculated from `:nevals`). Use `1` for the original VEGAS.
    - `:warmup` - number of the initial iterations which only train the grid and the allocation and are not used for the result, an integer >= 0 (default: `0`).
    - `:alpha` - grid refinement parameter, a number >= 0 (default: `0.5`, slow, for VEGAS+; `1.5`, moderate, when there is one stratum per dimension, i.e. for the original VEGAS with `:nstrats 1`).
    - `:beta` - damping of the stratification adaptation, a number >= 0 (default: `0.75`).
    - `:rel` - relative accuracy, a number >= 0 (default: `5.0e-4`).
    - `:abs` - absolute accuracy, a number >= 0 (default: `5.0e-4`).
    - `:random-sequence` - source of the samples: `:uniform` (default) or a low-discrepancy sequence `:r2`, `:sobol` or `:halton`.
    - `:jitter` - jittering factor of a low-discrepancy sequence (default: `0.75`).
    - `:rng` - random number generator (see [[fastmath.random/rng]]) which draws the samples, or the jitter of the low-discrepancy sequences. Seed it to reproduce a result. Default: a new `JDKRandomGenerator` for every call; the shared default generator is not used. The generator is used from one thread.
    - `:info?` - return the full information instead of the result (default: `false`).
    - `:record-data?` - store samples, numbers of strata, grid `:x` and `:dx` and the hypercubes in the result; requires `:info?` (default: `false`).

  The iterations stop when the standard deviation `sd` of the combined result `I` satisfies `sd <= max(rel * |I|, abs)` or when `:max-iters` iterations were done.

  Returns the estimate of the integral as a double. With `:info?` returns a map:

  - `:result` - the estimate of the integral.
  - `:sd` - standard deviation of the estimate.
  - `:iterations` - number of iterations, not counting the warmup.
  - `:evaluations` - number of calls of `f`, including the warmup.
  - `:nintervals` - actual size of the grid.
  - `:nstrats` - number of strata per dimension.
  - `:nhcubes` - number of hypercubes.
  - `:chi2-avg` - chi squared of the iterations per degree of freedom, about `1` for consistent iterations.
  - `:dof` - degrees of freedom, the number of iterations minus one.
  - `:Q` - goodness of fit of the iterations: close to `1` is very good, below `0.25` is very poor.
  - `:data` and `:hcubes` - recorded data, only with `:record-data?`.

  Corner cases:

  - after a single iteration `:dof` is `0` and `:chi2-avg` and `:Q` are `NaN`.
  - for a constant integrand the variances are zero: `:sd` is about zero and `:chi2-avg` and `:Q` are not meaningful.
  - an equal `lower` and `upper` in any dimension returns `0.0` (with `:info?` a map with zero `:sd`, `:iterations` and `:evaluations`).
  - throws `ex-info` with `{:option :value :expected}` for an invalid option value (unknown options are ignored), and with `{:rng value}` for an `:rng` which is not a generator.
  - throws `ex-info` with `{:x point :value v}` when `f` returns `NaN` or an infinite value.
  - throws `ex-info` for empty bounds and for bounds of different lengths.

  See also [[cubature]], [[integrate]]."
  ([f lower upper] (vegas f lower upper nil))
  ([f lower upper options] (vegas/vegas f lower upper options)))

(defn cubature
  "Adaptive cubature of a multivariate function on a box.

  The Genz-Malik degree 7 rule is applied to the box with the largest error estimate, which is then divided in two along the dimension with the largest fourth difference (h-adaptive cubature, a port of HCubature.jl). A one-dimensional integral is computed by the Gauss-Kronrod rule of [[integrate]]. Improper integrals (infinite bounds) are handled by a substitution.

  Parameters:

  - `f` (function): the integrand, a function of a vector of doubles returning a double.
  - `lower`, `upper` (sequences of numbers): bounds of the box, one per dimension, of the same length (one or more dimensions). A dimension with `lower` above `upper` is reversed: the result changes the sign for each reversed dimension. When `lower` equals `upper` in any dimension the integral is `0.0`.
  - `options` (map, optional):
    - `:initdiv` - number of the equal pieces each dimension is divided in at the start, an integer >= 1 (default: `2`).
    - `:max-evals` - maximum number of evaluations of `f`, never exceeded after the first round (default: `Integer/MAX_VALUE`).
    - `:max-iters` - maximum number of rounds (evaluations of the rule), including the first one (default: `64`).
    - `:rel` - relative error (default: `1.0e-7`).
    - `:abs` - absolute error (default: `1.0e-7`).
    - `:info?` - return the full information instead of the result (default: `false`).

  The calculation stops when the error estimate `E` of the result `I` satisfies `E <= max(rel * |I|, abs)`, otherwise when `:max-iters` or `:max-evals` would be exceeded by the next division.

  Returns the estimate of the integral as a double; when a limit stopped the calculation it is the best estimate found. With `:info?` returns a map:

  - `:result` - the estimate of the integral.
  - `:error` - the estimate of the error.
  - `:iterations` - number of rounds.
  - `:evaluations` - number of evaluations of `f`, equal to the number of calls.
  - `:subdivisions` - final number of boxes.
  - `:fail?` - `false` when converged, `:max-iters` or `:max-evals` when the limit was reached without convergence.

  Corner cases:

  - an integrand with `NaN` values never converges and ends with a limit.
  - a limit is checked after the convergence: a result which converged in the last allowed round has `:fail?` `false`. When the first round alone uses more than `:max-evals`, its result is returned with `:fail?` `:max-evals`.
  - throws `ex-info` with `{:option :value :expected}` when `:initdiv` is not an integer >= 1.
  - throws `ex-info` for empty bounds and for bounds of different lengths.

  See also [[vegas]], [[integrate]]."
  ([f lower upper] (cubature f lower upper nil))
  ([f lower upper options] (cubature/cubature f lower upper options)))

;;

(defn- make-integrant
  [f]
  (reify UnivariateFunction
    (value [_ x] (f x))))

(def ^:private commons-max-evals 10000000)

(defn- integrate-commons-nonempty
  "Integration with an Apache Commons Math integrator on a non-empty interval of oriented, finite or infinite bounds."
  [f ^double lower ^double upper
   {:keys [rel abs min-iters max-iters max-evals info? integrator integration-points]}]
  (let [rel (double (or rel BaseAbstractUnivariateIntegrator/DEFAULT_RELATIVE_ACCURACY))
        abs (double (or abs BaseAbstractUnivariateIntegrator/DEFAULT_ABSOLUTE_ACCURACY))
        min-iters (int (or min-iters BaseAbstractUnivariateIntegrator/DEFAULT_MIN_ITERATIONS_COUNT))
        max-evals (int (or max-evals commons-max-evals))
        max-iters (int (or max-iters
                           (case integrator
                             :romberg RombergIntegrator/ROMBERG_MAX_ITERATIONS_COUNT
                             :trapezoid TrapezoidIntegrator/TRAPEZOID_MAX_ITERATIONS_COUNT
                             :midpoint MidPointIntegrator/MIDPOINT_MAX_ITERATIONS_COUNT
                             :simpson SimpsonIntegrator/SIMPSON_MAX_ITERATIONS_COUNT
                             :gauss-legendre 64
                             BaseAbstractUnivariateIntegrator/DEFAULT_MAX_ITERATIONS_COUNT)))
        integrator-obj (case integrator
                         :romberg (RombergIntegrator. rel abs min-iters max-iters)
                         :trapezoid (TrapezoidIntegrator. rel abs min-iters max-iters)
                         :midpoint (MidPointIntegrator. rel abs min-iters max-iters)
                         :simpson (SimpsonIntegrator. rel abs min-iters max-iters)
                         :gauss-legendre (IterativeLegendreGaussIntegrator.
                                          (int (quadrature/check-integration-points (if (nil? integration-points) 7 integration-points)))
                                          rel abs min-iters max-iters)
                         (if (instance? UnivariateIntegrator integrator)
                           integrator
                           (throw (ex-info "Unknown integrator"
                                           {:option :integrator :value integrator
                                            :expected #{:gauss-kronrod :gauss-legendre :romberg :simpson :midpoint :trapezoid}}))))
        [f ^double lower ^double upper] (cc/subst-1d f lower upper)
        result (try
                 (.integrate ^UnivariateIntegrator integrator-obj max-evals (make-integrant f) lower upper)
                 (catch TooManyEvaluationsException _ :max-evals)
                 (catch MaxCountExceededException _ :max-iters))
        fail? (keyword? result)
        result' (if fail? ##NaN result)]
    (if info?
      {:result result'
       :iterations (.getIterations ^UnivariateIntegrator integrator-obj)
       :evaluations (.getEvaluations ^UnivariateIntegrator integrator-obj)
       :integrator integrator
       :fail? (if fail? result false)}
      result')))

(defn- integrate-commons
  "Integration with an Apache Commons Math integrator on oriented, finite or infinite bounds."
  [f ^double lower ^double upper {:keys [info? integrator] :as options}]
  (if (== lower upper)
    ;; the Commons integrators reject an empty interval
    (if info?
      {:result 0.0 :iterations 0 :evaluations 0 :integrator integrator :fail? false}
      0.0)
    (integrate-commons-nonempty f lower upper options)))

(defn integrate
  "Integrates a univariate function on an interval.

  The default integrator is the adaptive Gauss-Kronrod quadrature: the interval with the largest error estimate is divided in two until the total error is small enough (a port of QuadGK.jl). The other integrators come from Apache Commons Math. Improper integrals (infinite bounds) are mapped to a finite interval by a substitution.

  Parameters:

  - `f` (function): double -> double.
  - `lower`, `upper` (numbers): bounds, finite or infinite (default: `0.0` and `1.0`). Reversed bounds change the sign of the result, also for infinite bounds.
  - `options` (map, optional):
    - `:integrator` - `:gauss-kronrod` (default), `:gauss-legendre`, `:romberg`, `:simpson`, `:midpoint`, `:trapezoid`, or an instance of `org.apache.commons.math3.analysis.integration.UnivariateIntegrator`, which keeps its own tolerances and limits (only `:max-evals` and `:info?` are used).
    - `:rel` - relative error (default: `1.0e-8` for `:gauss-kronrod`, `1.0e-6` for the other integrators).
    - `:abs` - absolute error (default: `1.0e-8` for `:gauss-kronrod`, `1.0e-15` for the other integrators). The default absolute error accepts a result which is smaller than it, e.g. a narrow peak missed by all the quadrature points; see `:peak-guard`.
    - `:max-iters` - maximum number of iterations (default: `64`; `32` for `:romberg`).
    - `:max-evals` - maximum number of evaluations (default: `Integer/MAX_VALUE` for `:gauss-kronrod`, `10000000` for the other integrators, which would otherwise run practically forever for an integral which does not converge).
    - `:min-iters` - minimum number of iterations (default: `3`); not used by `:gauss-kronrod`. `:max-iters` of the other integrators has to be larger than `:min-iters`.
    - `:integration-points` - for `:gauss-kronrod` the order `n` of the rule (an integer >= 1, `n` Gauss points inside `2n + 1` Kronrod points; default: `7`, which is the 15 point rule), for `:gauss-legendre` the number of points (default: `7`).
    - `:initdiv` - for `:gauss-kronrod`, the number of equal pieces the interval (or each piece between `:points`) is divided in at the start, an integer >= 1 (default: `1`).
    - `:points` - for `:gauss-kronrod`, a sequence of break points inside the interval, in the variable of the integrand also for infinite bounds. The interval is divided at them at the start, which is useful for integrands with a singularity or a kink at a known point. Points outside of the interval are ignored.
    - `:peak-guard` - for `:gauss-kronrod`, when true the default `:abs` is `0.0`, so that the result is accepted only when the relative error is met (default: `false`). It finds peaks which the first quadrature points miss, but an integral equal to `0` is not accepted before `:max-iters` is reached.
    - `:info?` - return the full information instead of the result (default: `false`).

  Returns the integral as a double. With `:info?` returns a map:

  - `:result` - the value of the integral.
  - `:error` - estimate of the error (`:gauss-kronrod` only).
  - `:iterations` - number of iterations.
  - `:evaluations` - number of calls of `f`.
  - `:subdivisions` - final number of intervals (`:gauss-kronrod` only).
  - `:integrator` - the integrator (the other integrators only).
  - `:fail?` - `false` when converged, `:max-iters` or `:max-evals` when a limit was reached without convergence.

  A limit reached without convergence is reported in `:fail?`, but the results differ: `:gauss-kronrod` returns the best estimate found, the other integrators return `NaN` as the result.

  Corner cases:

  - an interval of zero length gives `0.0`.
  - the integrand must be finite inside of the interval; the Gauss-Kronrod rule does not evaluate it at the bounds, but the Romberg, trapezoid and Simpson rules do, also at the end points of the substituted infinite interval.
  - slowly decaying integrands of infinite intervals (like `x^-1.5`) do not converge with the trapezoid, midpoint, Simpson and Romberg rules; the evaluation limit stops them.
  - throws `ex-info` with `{:option :value :expected}` when `:integration-points` or `:initdiv` is not an integer >= 1, when `:points` contains a non-number and for an unknown `:integrator`.
  - throws `ex-info` for NaN bounds and for two equal infinite bounds.

  See also [[cubature]], [[vegas]], [[integrate-samples]]."
  ([f] (integrate f 0.0 1.0))
  ([f lower upper] (integrate f lower upper nil))
  ([f lower upper options]
   (let [[sign lower upper] (cc/normalize-bounds lower upper)
         {:keys [integrator abs peak-guard]} options
         integrator (if (nil? integrator) :gauss-kronrod integrator)]
     (cc/orient sign
                (if (= integrator :gauss-kronrod)
                  (quadrature/gk-quadrature f lower upper (if (and peak-guard (nil? abs))
                                                            (assoc options :abs 0.0)
                                                            options))
                  (integrate-commons f lower upper (assoc options :integrator integrator)))))))

(defn integrate-samples
  "Integrates sampled data with the composite trapezoid or Simpson rule.

  The samples `ys` are taken either at equally spaced points with a known spacing, or at arbitrary (not necessarily uniform) abscissae. The result is the integral over the whole sampled range, or, in the cumulative mode, the running integral at every sample. The results match `scipy.integrate.trapezoid`, `simpson`, `cumulative_trapezoid` and `cumulative_simpson` (the last two with `initial=0`).

  Parameters:

  - `ys` (sequence or `double[]` of numbers): the samples.
  - `x-or-dx` (optional): either a sequence or `double[]` of abscissae, of the same length as `ys`, or a number `dx` being the spacing between consecutive samples. Defaults to the spacing `1.0` (also when `nil`). A negative or zero spacing is allowed and is used as is, as a decreasing `x`.
  - `options` (optional map):
    - `:method` (`:trapezoid` or `:simpson`, default `:trapezoid`): the quadrature rule.
    - `:cumulative?` (boolean, default `false`): when true the running integral is returned instead of the total.

  Returns a `double` when `:cumulative?` is false. Otherwise returns a vector of doubles of the same length as `ys`, where the first element is `0.0` and the i-th element is the integral from the first to the i-th sample.

  The trapezoid rule is exact for linear data. The Simpson rule fits a parabola through each three consecutive samples; the total is exact for polynomials up to degree 3 on a uniform grid with an odd number of samples, and up to degree 2 otherwise. For an even number of samples the last interval is treated with the three-point correction of Cartwright, so any number of samples is handled. The cumulative Simpson integral is computed over every interval from the parabola through the nearest three samples. Simpson with fewer than 3 samples is the trapezoid rule.

  Corner cases:

  - no samples: returns `0.0`, or `[]` in the cumulative mode.
  - one sample: returns `0.0`, or `[0.0]` in the cumulative mode.
  - zero spacings in `x` are not rejected: the trapezoid rule gives zero for such an interval, the total Simpson rule treats the undefined terms as zero, the cumulative Simpson rule gives NaN or infinite values.
  - throws `ex-info` with `:ys-count` and `:x-count` in the ex-data when `x-or-dx` is a sequence of a different length than `ys`.
  - throws `ex-info` with `:option`, `:value` and `:expected` in the ex-data for an unknown `:method`, and for `ys` or `x-or-dx` which are neither numbers nor sequences as required.

  See also [[integrate]] (integration of functions)."
  ([ys] (samples/integrate-samples ys))
  ([ys x-or-dx] (samples/integrate-samples ys x-or-dx))
  ([ys x-or-dx options] (samples/integrate-samples ys x-or-dx options)))

;;; finite differences

(defn fx->gx+h
  "Converts `f(x)` to the function `g(x, h) = f(x + h)`, which is the input of [[extrapolate]].

  Parameters:

  - `f` (function): double -> double.

  Returns a function of two doubles `x` and `h`.

  See also [[fx->gx-h]], [[extrapolate]]."
  [f]
  (fn [^double x ^double h]
    (f (+ x h))))

(defn fx->gx-h
  "Converts `f(x)` to the function `g(x, h) = f(x - h)`, which is the input of [[extrapolate]].

  Parameters:

  - `f` (function): double -> double.

  Returns a function of two doubles `x` and `h`.

  See also [[fx->gx+h]], [[extrapolate]]."
  [f]
  (fn [^double x ^double h]
    (f (- x h))))

(defn extrapolate
  "Richardson extrapolation of a function `g(x, h)` to `h -> 0`, returns the extrapolated function of `x`.

  The function `g` is evaluated for a geometric sequence of steps `h` (started at `:init-h` and multiplied by `:contract` each time). The values are combined with Neville's scheme, which eliminates consecutive powers of `h` from the error expansion. The result is the combination with the smallest error estimate. The method is a port of Richardson.jl.

  Parameters:

  - `g` (function): the function of `x` and `h`, both doubles, returning a double. It has to be defined (and finite) for the steps used.
  - `options` (map, optional):
    - `:contract` - contraction factor of the step, must be in the open interval `(0, 1)` (default: `0.5`).
    - `:power` - the power of `h` in the leading error term; the next terms are assumed to be in powers of `contract^power` steps, e.g. `2.0` for the central differences (default: `1.0`).
    - `:init-h` - the first step, a nonzero number (default: `0.5`).
    - `:rel` - relative tolerance of the error estimate (default: square root of the ulp of `:init-h`).
    - `:abs` - absolute tolerance of the error estimate (default: machine epsilon).
    - `:max-evals` - maximum number of evaluations of `g`, an integer >= 1 (default: `Integer/MAX_VALUE`).
    - `:tol` - the iteration stops when the best error estimate grows by more than this factor with respect to the previous one (default: `2.0`).

  The extrapolation stops when the error estimate is not larger than `max(rel * |result|, abs)`, when the error estimate starts to grow, or when `:max-evals` is reached. In the last case the best value found so far is returned.

  For an infinite `x` the function `g(0, 1/h)` is evaluated for `x = ##Inf` and `g(0, -1/h)` for `x = ##-Inf`, which allows to extrapolate limits at infinity. `NaN` gives `NaN`.

  The stopping rule differs from Richardson.jl: the error estimate of the previous iteration is tested against the tolerance (one iteration later than in the original), therefore `g` can be evaluated once more. The result is not less accurate.

  Throws `ex-info` with `{:option :value :expected}` when `:contract` is not in `(0, 1)`, `:max-evals` is not an integer >= 1 or `:init-h` is zero or not finite.

  See also [[derivative]] (the `:extrapolate?` option), [[fx->gx+h]], [[fx->gx-h]]."
  ([g] (extrapolate g nil))
  ([g options] (finite/extrapolate g options)))

(defn derivative
  "Creates the `n`-th derivative of a univariate function `f` approximated with finite differences.

  Parameters:

  - `f` (function): double -> double.
  - `n` (integer >= 0): the order of the derivative (default: `1`). For `0` the function `f` itself is returned.
  - `options` (map, optional):
    - `:acc` - order of accuracy, an integer >= 1 (default: `2`). The error of the approximation is `O(h^acc)`; an odd `acc` for the central stencil gives the next even accuracy.
    - `:h` - step, a number >= 0 (default: `0`). `0` means the automatic step `eps^(1/(n+acc)) * max(1, |x|/10)`, where `eps` is the machine epsilon and `acc` is the accuracy (for the central method rounded up to even); this balances the truncation and rounding errors and grows with the magnitude of the point above `10`. An explicit positive `:h` is an absolute step.
    - `:method` - `:central` (default), `:forward` or `:backward`.
    - `:extrapolate?` - when true or a map with options of [[extrapolate]] the derivative is Richardson-extrapolated to `h -> 0` (default: `false`); the extrapolation starts with the step given in its options (`:init-h`, default `0.5`) and `:power` is `2` for the central and `1` otherwise, unless given.

  Returns a function of `x` returning the derivative as a double. The function has also the arities `(d x h)` with the absolute step `h` (0 means automatic) and `(d x h hn)` with the step and its precomputed power `h^n`. It is the function `(d x h)` that [[extrapolate]] uses.

  Throws `ex-info` with `{:option :value :expected}` when `n` is not an integer >= 0, `:acc` is not an integer >= 1, `:method` is unknown or `:h` is negative or not finite.

  Higher orders and high accuracies need smaller relative steps and are more sensitive to rounding errors: the attainable relative error grows quickly with `n`.

  See also [[f']], [[gradient]], [[hessian]], [[extrapolate]], [[fastmath.calculus.finite/fd-coeffs-for-accuracy]]."
  ([f] (derivative f 1))
  ([f n] (derivative f n nil))
  ([f n options] (finite/derivative f n options)))

(defn f'
  "Creates the first central derivative of `f` with order of accuracy 2.

  Parameters:

  - `f` (function): double -> double.

  Returns a function of `x`, the same as `(derivative f 1)`.

  See also [[derivative]], [[f'']], [[f''']]."
  [f] (finite/derivative f 1))

(defn f''
  "Creates the second central derivative of `f` with order of accuracy 2.

  Parameters:

  - `f` (function): double -> double.

  Returns a function of `x`, the same as `(derivative f 2)`.

  See also [[derivative]], [[f']], [[f''']]."
  [f] (finite/derivative f 2))

(defn f'''
  "Creates the third central derivative of `f` with order of accuracy 2.

  Parameters:

  - `f` (function): double -> double.

  Returns a function of `x`, the same as `(derivative f 3)`.

  See also [[derivative]], [[f']], [[f'']]."
  [f] (finite/derivative f 3))

(defn gradient
  "Creates the gradient (vector of the first partial derivatives) of a multivariate scalar function `f`, approximated with central finite differences.

  Parameters:

  - `f` (function): a vector of doubles -> double.
  - `options` (map, optional):
    - `:acc` - order of accuracy, an integer >= 1 (default: `2`); an odd value gives the next even accuracy.
    - `:h` - step, a number >= 0 (default: `0`). `0` means the automatic step for each component `i`: `eps^(1/(1+acc)) * max(1, |x_i|/10)` (with `acc` rounded up to even). An explicit positive `:h` is an absolute step used for every component.

  Returns a function of a point `v` (any sequence of numbers or a `double[]`) returning the gradient as a vector of doubles. The function has also the arity `(g v h)` with an explicit absolute step (0 means automatic). The function `f` receives a vector and is called when the gradient is called, the result is not lazy. The gradient of the empty point is `[]`.

  Throws `ex-info` with `{:option :value :expected}` when `:acc` is not an integer >= 1 or `:h` is negative or not finite.

  See also [[hessian]], [[jacobian]], [[derivative]]."
  ([f] (gradient f nil))
  ([f options] (finite/gradient f options)))

(defn jacobian
  "Creates the Jacobian matrix of a vector valued function `f`, approximated with central finite differences.

  Parameters:

  - `f` (function): a vector of `n` doubles -> a sequence (or array) of `m` numbers.
  - `options` (map, optional): as in [[gradient]]
    - `:acc` - order of accuracy, an integer >= 1 (default: `2`).
    - `:h` - step, a number >= 0 (default: `0` = automatic step for each component, `eps^(1/(1+acc)) * max(1, |x_i|/10)`); an explicit positive value is an absolute step.

  Returns a function of a point `v` (any sequence of numbers or a `double[]`) returning the matrix as a vector of `m` rows, each a vector of `n` partial derivatives `d f_j / d x_i`. The function has also the arity `(j v h)` with an explicit absolute step (0 means automatic). The function `f` is called when the Jacobian is called. For the empty point the result is a vector of `m` empty rows.

  Throws `ex-info` with `{:option :value :expected}` when `:acc` or `:h` are invalid, when `f` returns a number instead of a sequence or the number of returned values changes between calls.

  See also [[gradient]], [[hessian]]."
  ([f] (jacobian f nil))
  ([f options] (finite/jacobian f options)))

(defn hessian
  "Creates the Hessian matrix of a multivariate scalar function `f`, approximated with central finite differences of the second order.

  Parameters:

  - `f` (function): a vector of doubles -> double.
  - `options` (map, optional):
    - `:h` - step, a number >= 0 (default: `0`). `0` means the automatic step for each component `i`: `eps^(1/4) * max(1, |x_i|/10)`. An explicit positive `:h` is an absolute step used for every component.

  Returns a function of a point `v` (any sequence of numbers or a `double[]`) returning the symmetric matrix as a vector of vectors of doubles. The function has also the arity `(hs v h)` with an explicit absolute step (0 means automatic). The diagonal term is `(f(x+h_i) - 2 f(x) + f(x-h_i)) / h_i^2` and the mixed term is `(f(++) - f(+-) - f(-+) + f(--)) / (4 h_i h_j)`. The function `f` is called when the Hessian is called. The Hessian of the empty point is `[]`.

  Throws `ex-info` with `{:option :value :expected}` when `:h` is negative or not finite.

  See also [[hessian-diagonal]], [[laplacian]], [[gradient]]."
  ([f] (hessian f nil))
  ([f options] (finite/hessian f options)))

(defn hessian-diagonal
  "Creates the function returning the diagonal of the Hessian matrix of a multivariate scalar function `f`, approximated with central finite differences of the second order.

  Parameters:

  - `f` (function): a vector of doubles -> double.
  - `options` (map, optional):
    - `:h` - step, a number >= 0 (default: `0`). `0` means the automatic step for each component `i`: `eps^(1/4) * max(1, |x_i|/10)`. An explicit positive `:h` is an absolute step used for every component.

  Returns a function of a point `v` (any sequence of numbers or a `double[]`) returning the diagonal `(f(x+h_i) - 2 f(x) + f(x-h_i)) / h_i^2` as a vector of doubles. The function has also the arity `(hd v h)` with an explicit absolute step (0 means automatic). The function `f` is called when the diagonal is computed. It needs `2n + 1` evaluations of `f` for `n` components, while [[hessian]] needs more.

  Throws `ex-info` with `{:option :value :expected}` when `:h` is negative or not finite.

  See also [[hessian]], [[laplacian]]."
  ([f] (hessian-diagonal f nil))
  ([f options]
   (let [hd (finite/hessian-diagonal f options)]
     (fn local-hessian-diagonal
       ([v] (vec (hd v)))
       ([v h] (vec (hd v h)))))))

(defn laplacian
  "Creates the Laplacian (the sum of the second partial derivatives) of a multivariate scalar function `f`, approximated with central finite differences.

  Parameters:

  - `f` (function): a vector of doubles -> double.
  - `options` (map, optional): as in [[hessian-diagonal]]
    - `:h` - step, a number >= 0 (default: `0` = automatic step for each component, `eps^(1/4) * max(1, |x_i|/10)`); an explicit positive value is an absolute step.

  Returns a function of a point `v` (any sequence of numbers or a `double[]`) returning the Laplacian as a double. The function has also the arity `(l v h)` with an explicit absolute step (0 means automatic). The Laplacian of the empty point is `0.0`.

  Throws `ex-info` with `{:option :value :expected}` when `:h` is negative or not finite.

  See also [[hessian-diagonal]], [[hessian]]."
  ([f] (laplacian f nil))
  ([f options] (finite/laplacian f options)))

;; divided differences

(defn divided-all
  "Divided differences of `f` at the nodes `xs`: all levels of the divided differences table.

  The level `k` consists of the divided differences `f[x_i, ..., x_(i+k)]` of `k + 1` consecutive nodes. Repeated nodes are allowed and give the Hermite (confluent) divided differences: for `k + 1` equal nodes the value is `f^(k)(x) / k!`.

  Parameters:

  - `f` (function): the function of a node. A vector or a map (lookup by the node) is also accepted as `f`, but then the nodes have to be distinct or `:derivatives` has to be given.
  - `xs` (sequence of numbers): the nodes, at least one. Equal nodes must be adjacent.
  - `options` (map, optional):
    - `:derivatives` - sequence of functions `[f' f'' ...]` giving the derivatives of `f`; the element `k - 1` is the `k`-th derivative. It has to contain the derivatives up to the order of the longest run of equal nodes minus one. When it is not given, the derivatives are approximated with [[derivative]] with Richardson extrapolation (`:extrapolate? true`) started with the step `0.5` (divided by 10, up to 8 times, when `f` is not finite so far from the node). It is accurate to about 1e-6 or better for the first and the second derivative of smooth functions, and gets worse for higher orders. `f` has to be a smooth function of a double.

  Returns a list of levels in the order from the highest to the lowest: the first element is a list with the single divided difference of all nodes, the last element is the list of the values of `f` at the nodes. Every level is a list of doubles. For a single node the result is `((f(x0)))`. The first (highest) element is the value returned by [[divided]].

  Throws `ex-info` with `{:option :value :expected}` when `xs` is empty, equal nodes are not adjacent, `f` is not callable, `f` is not a function while nodes are repeated and `:derivatives` is not given, or `:derivatives` is too short.

  See also [[divided]]."
  ([f xs] (divided/divided-all f xs))
  ([f xs options] (divided/divided-all f xs options)))

(defn divided
  "Divided difference `f[x0, x1, ..., xn]` of `f` at the nodes `xs`, the highest level of the table.

  It is the leading coefficient of the Newton interpolation polynomial through the nodes. Repeated nodes are allowed and give the Hermite (confluent) divided differences: for `k + 1` equal nodes the value is `f^(k)(x) / k!`.

  Parameters:

  - `f` (function): the function of a node. A vector or a map (lookup by the node) is also accepted as `f`, but then the nodes have to be distinct or `:derivatives` has to be given.
  - `xs` (sequence of numbers): the nodes, at least one. Equal nodes must be adjacent.
  - `options` (map, optional):
    - `:derivatives` - sequence of functions `[f' f'' ...]` giving the derivatives of `f`; the element `k - 1` is the `k`-th derivative. It has to contain the derivatives up to the order of the longest run of equal nodes minus one. When it is not given, the derivatives are approximated with [[derivative]] with Richardson extrapolation (`:extrapolate? true`) started with the step `0.5` (divided by 10, up to 8 times, when `f` is not finite so far from the node). It is accurate to about 1e-6 or better for the first and the second derivative of smooth functions, and gets worse for higher orders. `f` has to be a smooth function of a double.

  Returns the divided difference as a double. For a single node it is `f(x0)`.

  Throws `ex-info` with `{:option :value :expected}` when `xs` is empty, equal nodes are not adjacent, `f` is not callable, `f` is not a function while nodes are repeated and `:derivatives` is not given, or `:derivatives` is too short.

  See also [[divided-all]]."
  (^double [f xs] (divided/divided f xs))
  (^double [f xs options] (divided/divided f xs options)))
