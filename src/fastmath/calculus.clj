(ns fastmath.calculus
  "Integration and derivatives

  Integrate univariate and multivariate functions.

  * VEGAS / VEGAS+ - Monte Carlo integration of multivariate function
  * h-Cubature - h-adaptive integration of multivariate function
  * Guass-Kronrod and Gauss-Legendre - quadrature integration of univariate functions
  * Romberg, Simpson, MidPoint and Trapezoid

  Integrant is substituted in case of improper integration bounds.

  Derivatives (finite differences method):

  * derivatives of any degree and any order of accuracy
  * gradient and hessian for multivariate functions"
  (:require [fastmath.core :as m]

            [fastmath.calculus.common :as cc]
            [fastmath.calculus.vegas :as vegas]
            [fastmath.calculus.cubature :as cubature]
            [fastmath.calculus.quadrature :as quadrature]
            [fastmath.calculus.finite :as finite])
  (:import [org.apache.commons.math3.analysis UnivariateFunction]
           [org.apache.commons.math3.analysis.integration UnivariateIntegrator
            IterativeLegendreGaussIntegrator MidPointIntegrator
            RombergIntegrator SimpsonIntegrator TrapezoidIntegrator
            BaseAbstractUnivariateIntegrator]
           
           [org.apache.commons.math3.exception MaxCountExceededException TooManyEvaluationsException]))

#_(set! *warn-on-reflection* true)
#_(set! *unchecked-math* :warn-on-boxed)

(defn vegas
  "VEGAS+ - Monte Carlo integration of multivariate function, n>1 dimensions.

  Improper integrals with infinite bounds are handled by a substitution.

  Arguments:

  * `f` - integrant
  * `lower` - seq of lower bounds
  * `upper` - seq of upper bounds
  
  Additional options:

  * `:max-iters` - maximum number of iterations, default: 10
  * `:nevals` - number of evaluations per iteration, default: 10000
  * `:nintervals` - number of grid intervals per dimension (default: 1000)
  * `:nstrats` - number of stratifications per dimension (calculated)
  * `:warmup` - number of warmup iterations (results are used to train stratification and grid spacings, default: 0
  * `:alpha` - grid refinement parameter, 0.5 slow (default for vegas+), 1.5 moderate/fast (defatult for vegas)
  * `:beta` - stratification damping parameter for startification adaptation, default: 0.75
  * `:rel` - relative accuracy, default: 5.0e-4
  * `:abs` - absolute accuracy, default: 5.0e-4
  * `:random-sequence` - random sequence used for generating samples: `:uniform` (default), low-discrepancy sequences: `:r2`, `:sobol` and `:halton`.
  * `:jitter` - jittering factor for low-discrepancy random sequence, default: 0.75
  * `:info?` - return full information about integration, default: false
  * `:record-data?` - stores samples, number of strata, x and dx, default: false (requires, `:info?` to be set to `true`)

  For original VEGAS algorithm set `:nstrats` to `1`.

  `:nstrats` can be also a list, then each dimension is divided independently according to a given number. If list is lower then number of dimensions, then it's cycled.

  Function returns a map with following keys (if info? is true, returns result otherwise):

  * `:result` - value of integral
  * `:iterations` - number of iterations (excluding warmup)
  * `:sd` - standard deviation of results
  * `:nintervals` - actual grid size
  * `:nstrats` - number of stratitfications per dimension
  * `:nhcubes` - number of hypercubes
  * `:evaluations` - number of function calls
  * `:chi2-avg` - average of chi2
  * `:dof` - degrees of freedom
  * `:Q` - goodness of fit indicator, 1 - very good, <0.25 very poor
  * `:data` - recorded data (if available)"
  ([f lower upper] (vegas f lower upper nil))
  ([f lower upper options] (vegas/vegas f lower upper options)))

(defn cubature
  "Cubature - h-adaptive integration of multivariate function, n>1 dimensions.

  Algorithm uses Genz Malik method.

  In each iteration a box with biggest error is subdivided and reevaluated.

  Improper integrals with infinite bounds are handled by a substitution.

  Arguments:

  * `f` - integrant
  * `lower` - seq of lower bounds
  * `upper` - seq of upper bounds
  
  Options:

  * `:initvid` - initial subdivision per dimension, default: 2.
  * `:max-evals` - maximum number of evaluations, default: max integer value.
  * `:max-iters` - maximum number of iterations, default: 64.
  * `:rel` - relative error, 1.0e-7
  * `:abs` - absolute error, 1.0e-7
  * `:info?` - return full information about integration, default: false

  Function returns a map containing (if info? is true, returns result otherwise):

  * `:result` - integration value
  * `:error` - integration error
  * `:iterations` - number of iterations
  * `:evaluations` - number of evaluations
  * `:subdivisions` - final number of boxes
  * `:fail?` - set to `:max-evals` or `:max-iters` when one of the limits has been reached without the convergence."
  ([f lower upper] (cubature f lower upper nil))
  ([f lower upper options] (cubature/cubature f lower upper options)))

;;

(defn- make-integrant
  [f]
  (reify UnivariateFunction
    (value [_ x] (f x))))

(defn integrate
  "Univariate integration.

  Improper integrals with infinite bounds are handled by a substitution.

  Arguments:

  * `f` - integrant
  * `lower` - lower bound
  * `upper` - upper bound

  Options:

  * `:integrator` - integration algorithm, one of: `:romberg`, `:trapezoid`, `:midpoint`, `:simpson`, `:gauss-legendre` and `:gauss-kronrod` (default).
  * `:min-iters` - minimum number of iterations (default: 3), not used in `:gauss-kronrod`
  * `:max-iters` - maximum number of iterations (default: 32 or 64)
  * `:max-evals` - maximum number of evaluations, (default: maximum integer)
  * `:rel` - relative error
  * `:abs` - absolute error
  * `:integration-points` - number of integration (quadrature) points for `:gauss-legendre` and `:gauss-kronrod`, default 7
  * `:initdiv` - initial number of subdivisions for `:gauss-kronrod`, default: 1
  * `:info?` - return full information about integration, default: false

  `:gauss-kronrod` is h-adaptive implementation

  Function returns a map containing (if info? is true, returns result otherwise):

  * `:result` - integration value
  * `:error` - integration error (`:gauss-kronrod` only)
  * `:iterations` - number of iterations
  * `:evaluations` - number of evaluations
  * `:subdivisions` - final number of boxes (`:gauss-kronrod` only)
  * `:fail?` - set to `:max-evals` or `:max-iters` when one of the limits has been reached without the convergence."
  ([f] (integrate f 0.0 1.0))
  ([f ^double lower ^double upper] (integrate f lower upper nil))
  ([f ^double lower ^double upper
    {:keys [^double rel ^double abs max-iters ^int min-iters ^int max-evals info?
            integrator ^long integration-points]
     :or {rel BaseAbstractUnivariateIntegrator/DEFAULT_RELATIVE_ACCURACY
          abs BaseAbstractUnivariateIntegrator/DEFAULT_ABSOLUTE_ACCURACY
          min-iters BaseAbstractUnivariateIntegrator/DEFAULT_MIN_ITERATIONS_COUNT
          max-evals Integer/MAX_VALUE
          integration-points 7
          integrator :gauss-kronrod
          info? false}
     :as options}]
   (let [max-iters (unchecked-int (or max-iters
                                      (case integrator
                                        :romberg RombergIntegrator/ROMBERG_MAX_ITERATIONS_COUNT
                                        :trapezoid TrapezoidIntegrator/TRAPEZOID_MAX_ITERATIONS_COUNT
                                        :midpoint MidPointIntegrator/MIDPOINT_MAX_ITERATIONS_COUNT
                                        :simpson SimpsonIntegrator/SIMPSON_MAX_ITERATIONS_COUNT
                                        :gauss-legendre 64                                        
                                        BaseAbstractUnivariateIntegrator/DEFAULT_MAX_ITERATIONS_COUNT)))
         integrator-obj (if (keyword? integrator)
                          (case integrator
                            :romberg (RombergIntegrator. rel abs min-iters max-iters)
                            :trapezoid (TrapezoidIntegrator. rel abs min-iters max-iters)
                            :midpoint (MidPointIntegrator. rel abs min-iters max-iters)
                            :simpson (SimpsonIntegrator. rel abs min-iters max-iters)
                            :gauss-legendre (IterativeLegendreGaussIntegrator. integration-points
                                                                               rel abs min-iters max-iters)
                            integrator)
                          integrator)]
     (cond (instance? UnivariateIntegrator integrator-obj)
           (let [[f ^double lower ^double upper] (cc/subst-1d f lower upper)
                 integrant (make-integrant f)
                 result (try
                          (.integrate ^UnivariateIntegrator integrator-obj max-evals integrant lower upper)
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
               result'))

           (= integrator-obj :gauss-kronrod)
           (quadrature/gk-quadrature f lower upper options)
           
           :else (throw (Exception. (str "Unknown integrator: " integrator)))))))

;;; finite difference

(defn fx->gx+h
  "Convert f(x) to g(x,h)=f(x+h)"
  [f]
  (fn [^double x ^double h]
    (f (m/+ x h))))

(defn fx->gx-h
  "Convert f(x) to g(x,h)=f(x-h)"
  [f]
  (fn [^double x ^double h]
    (f (m/- x h))))

(defn extrapolate
  "Richardson extrapolation for given function `g=g(x,h)`. Returns extrapolated function f(x).

  Options:

  * `:contract` - shrinkage factor, default=`1/2`
  * `:power` - set to `2.0` for even functions around `x0`, default `1.0`
  * `:init-h` - initial step `h`, default=`1/2`
  * `:abs` - absolute error, default: machine epsilon
  * `:rel` - relative error, default: ulp for init-h
  * `:tol` - tolerance for error, default: `2.0` 
  * `:max-evals` - maximum evaluations, default: maximum integer"
  ([g] (extrapolate g nil))
  ([g options] (finite/extrapolate g options)))

(defn derivative
  "Create nth derivative of `f` using finite difference method for given accuracy `:acc` and step `:h`.

  Returns function.

  Arguments:

  * `n` - derivative
  * `:acc` - order of accuracy (default: 2)
  * `:h` - step, (default: 0.0, automatic)
  * `:method` - `:central` (default), `:forward` or `:backward`
  * `:extrapolate?` - creates extrapolated derivative if set to true or a map with [[extrapolate]] function options"
  ([f] (derivative f 1))
  ([f ^long n] (derivative f n nil))
  ([f ^long n options] (finite/derivative f n options)))

(defn f' "First central derivative with order of accuracy 2." [f] (finite/derivative f 1))
(defn f'' "Second central derivative with order of accuracy 2." [f] (finite/derivative f 2))
(defn f''' "Third central derivative with order of accuracy 2." [f] (finite/derivative f 3))

(defn gradient
  "Create first partial derivatives of multivariate function for given accuracy `:acc` and step `:h`.

  Returns function.

  Options:

  * `:acc` - order of accuracy, 2 (default) or 4.
  * `:h` - step, default `1.0e-6`"
  ([f] (gradient f nil))
  ([f options] (finite/gradient f options)))

(defn hessian
  "Creates function returning Hessian matrix for mulitvariate function `f` and given `:h` step (default: `5.0e-3`)."
  ([f] (hessian f nil))
  ([f options] (finite/hessian f options)))
