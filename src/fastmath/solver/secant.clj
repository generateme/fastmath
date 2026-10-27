(ns fastmath.solver.secant
  "The secant method for a root of a function, started from the ends of an interval.

  The new point is the root of the line through the last two points. The method has superlinear convergence (order about 1.62) but, unlike the bracketing methods, it does not keep a bracket: the iteration can leave the interval and run away from the root. It is a failure here, never a result: a result outside the interval, a flat step or a point which is not a root throws."
  (:require [fastmath.core :as m]
            [fastmath.solver.common :as c])
  (:import [fastmath.solver.common IEvaluator Problem]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- core
  ^double [^Problem p]
  (let [^IEvaluator ev (.ev p)
        atol (.atol p)
        rtol (.rtol p)
        lo (.a p)
        hi (.b p)]
    (loop [x0 lo f0 (.fa p) x1 hi f1 (.fb p)]
      (when (m/== f0 f1)
        (c/throw-status :diverged (str "The secant method stalled: equal function values at " x0 " and " x1)
                        {:evaluations (.evaluations ev) :x0 x0 :x1 x1}))
      (let [t (m// f1 (m/- f1 f0))
            x2 (m/+ (m/* (m/- 1.0 t) x1) (m/* t x0))]
        (when-not (m/valid-double? x2)
          (c/throw-status :diverged (str "The secant method ran away: " x2)
                          {:evaluations (.evaluations ev) :x x2}))
        (let [f2 (.eval ev x2)]
          (cond
            (m/zero? f2) x2
            (m/<= (m/abs (m/- x2 x1)) (c/tolerance atol rtol x2))
            (if (and (m/>= x2 lo) (m/<= x2 hi) (c/certified? ev x2 f2 p))
              x2
              (c/throw-status :diverged (str "The secant method converged to a point which is not a root of the interval [" lo ", " hi "]: " x2)
                              {:evaluations (.evaluations ev) :x x2 :fx f2}))
            :else (recur x1 f1 x2 f2)))))))

(defn secant
  "Finds a root of `f` with the secant method started from the ends of `[lower, upper]`.

  Parameters:

  - `f` (function): a function of one number returning a number; it has to change the sign on the interval.
  - `lower`, `upper` (finite numbers): ends of the interval; reversed ends are swapped. The result has to lie in the interval, the iteration itself may leave it.
  - `options` (map, optional):
      - `:absolute-accuracy` - absolute tolerance of the root, default `1.0e-300`.
      - `:relative-accuracy` - tolerance relative to the root, default four times the machine epsilon. The root is found within `absolute + relative * |root|`.
      - `:max-iters` - limit of function evaluations (an integer), default `1200`. Exceeding it throws.
      - `:initial-value` - ignored by this method.
      - `:info?` - when truthy, returns the map `{:root :evaluations :status}` instead of the root, and returns the failure status instead of throwing. Default `false`.

  Returns the root as a double. A root at an end of the interval is returned after at most two evaluations. The search stops when the step is within the tolerance and the function changes the sign within two tolerances of the point.

  Throws `ex-info` for an invalid argument, and for a failed solve: no sign change on the interval (`:status :no-bracket`), `NaN` function values (`:nan`), a non number result (`:invalid-result`), the evaluation limit (`:max-evals`) or a runaway iteration (`:diverged`): an iterate is not finite, two consecutive values are equal, or the step converged to a point which is not a root of the interval. The function is evaluated outside the interval when an iterate leaves it, which gives `:nan` when the function is not defined there. The data of the exception holds the `:status` and a few scalars.

  See also [[fastmath.solver/find-root]], [[fastmath.solver.brent/brent]], [[fastmath.solver.illinois/illinois]]."
  ([f lower upper] (secant f lower upper nil))
  ([f lower upper options]
   (c/solve-bracketed core false f lower upper options)))
