(ns fastmath.solver.bisection
  "The bisection method for a root of a function on an interval with a sign change.

  Every step halves the interval, so the number of evaluations is known in advance and the method never fails on a function with a sign change, whatever its shape, jumps included. It is the slowest of the bracketing methods and works on any scale of the root."
  (:require [fastmath.core :as m]
            [fastmath.solver.common :as c])
  (:import [fastmath.solver.common IEvaluator Problem]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- core
  ^double [^Problem p]
  (let [^IEvaluator ev (.ev p)
        atol (.atol p)
        rtol (.rtol p)]
    (loop [a (.a p) fa (.fa p) b (.b p)]
      (let [mid (m/+ (m/* 0.5 a) (m/* 0.5 b))]
        (if (or (m/<= (m/- (m/* 0.5 b) (m/* 0.5 a)) (c/tolerance atol rtol mid))
                (m/== mid a)
                (m/== mid b))
          mid
          (let [fm (.eval ev mid)]
            (cond
              (m/zero? fm) mid
              (m/== (m/signum fm) (m/signum fa)) (recur mid fm b)
              :else (recur a fa mid))))))))

(defn bisection
  "Finds a root of `f` on `[lower, upper]` by halving the interval.

  Parameters:

  - `f` (function): a function of one number returning a number; it has to change the sign on the interval.
  - `lower`, `upper` (finite numbers): ends of the interval; reversed ends are swapped.
  - `options` (map, optional):
      - `:absolute-accuracy` - absolute tolerance of the root, default `1.0e-300`.
      - `:relative-accuracy` - tolerance relative to the root, default four times the machine epsilon. The root is found within `absolute + relative * |root|`.
      - `:max-iters` - limit of function evaluations (an integer), default `1200`, enough to reach the floating point resolution on any interval. Exceeding it throws.
      - `:initial-value` - ignored by this method.
      - `:info?` - when truthy, returns the map `{:root :evaluations :status}` instead of the root, and returns the failure status instead of throwing. Default `false`.

  Returns the root as a double, the middle of the last interval. A root at an end of the interval is returned after at most two evaluations. At the default tolerance the result is within a few units in the last place; on a jump of the function it is the point of the jump.

  Throws `ex-info` for an invalid argument, and for a failed solve: no sign change on the interval (`:status :no-bracket`), `NaN` function values (`:nan`), a non number result (`:invalid-result`) or the evaluation limit (`:max-evals`). The data of the exception holds the `:status` and a few scalars.

  See also [[fastmath.solver/find-root]], [[fastmath.solver.brent/brent]]."
  ([f lower upper] (bisection f lower upper nil))
  ([f lower upper options]
   (c/solve-bracketed core false f lower upper options)))
