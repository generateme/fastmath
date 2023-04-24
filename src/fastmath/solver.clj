(ns fastmath.solver
  (:import [org.apache.commons.math3.analysis UnivariateFunction]
           [org.apache.commons.math3.analysis.solvers UnivariateSolver BrentSolver BisectionSolver IllinoisSolver MullerSolver MullerSolver2 PegasusSolver RegulaFalsiSolver RiddersSolver SecantSolver]))

(defn- wrap-univariate-function ^UnivariateFunction [f] (reify UnivariateFunction (value [_ x] (f x))))

(defmacro ^:private solver-constructor
  [clazz]
  `(cond
     ~'_relative-accuracy (new ~clazz (or ~'_absolute-accuracy 1.0e-8) ~'_relative-accuracy)
     ~'_absolute-accuracy (new ~clazz (or ~'_absolute-accuracy 1.0e-8))
     :else (new ~clazz)))

(defn- get-solver
  ^UnivariateSolver [solver-type _absolute-accuracy _relative-accuracy]
  (case solver-type
    :brent (solver-constructor BrentSolver)
    :bisection (solver-constructor BisectionSolver)
    :illinois (solver-constructor IllinoisSolver)
    :muller (solver-constructor MullerSolver)
    :muller2 (solver-constructor MullerSolver2)
    :pegasus (solver-constructor PegasusSolver)
    :regula-falsi (solver-constructor RegulaFalsiSolver)
    :ridders (solver-constructor RiddersSolver)
    :secant (solver-constructor SecantSolver)
    (solver-constructor BrentSolver)))

(defn find-root
  "Find zero (root) of a function `f` in given range [`lower-bound`, `upper-bound`].

  Optional parameters:

  * `:absolute-accuracy` - default 1.0e-8
  * `:relative-accuracy`
  * `:max-iters` - maximum iterations (default: 100)
  * `:initial-value` - algorithm starting value
  * `:solver` - one of: `:brent` (default), `:bisection`, `:illinois`, `:muller`, `:muller2`, `:pegasus`, `:regula-falsi`, `:ridders` and `:secant`."
  (^double [f lower-bound upper-bound] (find-root f lower-bound upper-bound {}))
  (^double [f lower-bound upper-bound {:keys [absolute-accuracy relative-accuracy
                                              ^int max-iters initial-value solver]
                                       :or {max-iters 100 solver :brent}}]
   (let [^UnivariateSolver solver (get-solver solver absolute-accuracy relative-accuracy)]
     (if initial-value
       (.solve solver max-iters (wrap-univariate-function f) (double lower-bound) (double upper-bound) (double initial-value))
       (.solve solver max-iters (wrap-univariate-function f) (double lower-bound) (double upper-bound))))))
