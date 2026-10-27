(ns fastmath.solver.halley
  "Halley's method for a root of a function with known derivatives, safeguarded by the bracket.

  The new point uses the first and the second derivative, which gives a cubic convergence near a simple root. Here the step is replaced by a bisection step whenever it would leave the bracket or does not converge fast enough, so the method does not fail on a function with a sign change. The derivatives are the options `:df` and `:d2f`; without them central finite differences (see [[fastmath.calculus/derivative]]) are used."
  (:require [fastmath.solver.common :as c]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn halley
  "Finds a root of `f` on `[lower, upper]` with Halley's method safeguarded by the bracket.

  Parameters:

  - `f` (function): a function of one number returning a number; it has to change the sign on the interval.
  - `lower`, `upper` (finite numbers): ends of the interval; reversed ends are swapped.
  - `options` (map, optional):
      - `:df`, `:d2f` - the first and the second derivative of `f`, functions of one number returning a number. Default central finite differences of order 4 and 2 (the evaluations of `f` made by them count to `:max-iters`).
      - `:absolute-accuracy` - absolute tolerance of the root, default `1.0e-300`.
      - `:relative-accuracy` - tolerance relative to the root, default four times the machine epsilon. The root is found within `absolute + relative * |root|`.
      - `:max-iters` - limit of evaluations of `f` (an integer), default `1200`; a call of a derivative is made at most once per evaluation of `f`. Exceeding the limit throws.
      - `:initial-value` - the first point, strictly inside the interval, default the middle of the interval.
      - `:info?` - when truthy, returns the map `{:root :evaluations :status}` instead of the root, and returns the failure status instead of throwing. Default `false`.

  Returns the root as a double. A root at an end of the interval is returned after at most two evaluations. A derivative which is zero, `NaN` or infinite never stops the search: the step is replaced by a bisection step.

  Throws `ex-info` for an invalid argument (an option which is not a function included), and for a failed solve: no sign change on the interval (`:status :no-bracket`), `NaN` function values (`:nan`), a non number result of `f` or of a derivative (`:invalid-result`) or the evaluation limit (`:max-evals`). The data of the exception holds the `:status` and a few scalars.

  See also [[fastmath.solver/find-root]], [[fastmath.solver.newton/newton]], [[fastmath.solver.brent/brent]]."
  ([f lower upper] (halley f lower upper nil))
  ([f lower upper options]
   (c/solve-with-derivatives :halley f lower upper options)))
