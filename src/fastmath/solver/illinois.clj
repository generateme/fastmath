(ns fastmath.solver.illinois
  "The Illinois method for a root of a function on an interval with a sign change.\n\n  A modification of the regula falsi: the value at the retained end of the bracket is halved when the same end would be kept twice, which restores the superlinear convergence (order about 1.44)."
  (:require [fastmath.solver.common :as c]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn illinois
  "Finds a root of `f` on `[lower, upper]` with the Illinois method.

  Parameters:

  - `f` (function): a function of one number returning a number; it has to change the sign on the interval.
  - `lower`, `upper` (finite numbers): ends of the interval; reversed ends are swapped.
  - `options` (map, optional):
      - `:absolute-accuracy` - absolute tolerance of the root, default `1.0e-300`.
      - `:relative-accuracy` - tolerance relative to the root, default four times the machine epsilon. The root is found within `absolute + relative * |root|`.
      - `:max-iters` - limit of function evaluations (an integer), default `1200`. Exceeding it throws.
      - `:initial-value` - ignored by this method.
      - `:info?` - when truthy, returns the map `{:root :evaluations :status}` instead of the root, and returns the failure status instead of throwing. Default `false`.

  Returns the root as a double. A root at an end of the interval is returned after at most two evaluations. The new point always stays inside the bracket, so the method does not leave the interval; it stops when the bracket is narrower than the tolerance, or when the step is within the tolerance and the point is confirmed by a sign change at a distance of two tolerances.

  Throws `ex-info` for an invalid argument, and for a failed solve: no sign change on the interval (`:status :no-bracket`), `NaN` function values (`:nan`), a non number result (`:invalid-result`) or the evaluation limit (`:max-evals`). The data of the exception holds the `:status` and a few scalars.

  See also [[fastmath.solver/find-root]], [[fastmath.solver.regula-falsi/regula-falsi]], [[fastmath.solver.pegasus/pegasus]], [[fastmath.solver.brent/brent]]."
  ([f lower upper] (illinois f lower upper nil))
  ([f lower upper options]
   (c/solve-bracketed (partial c/false-position :illinois) false f lower upper options)))
