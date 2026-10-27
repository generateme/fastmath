(ns fastmath.solver.fixed-point
  "Fixed point iteration for a solution of `g(x) = x`.

  The iteration `x -> g(x)` converges to a fixed point when `g` is a contraction near it, with a linear rate `|g'|`. The optional Steffensen acceleration (Aitken's delta squared process applied to every second step) converges quadratically when the plain iteration converges. Unlike the root finders, the method needs a start and no bracket: it fails with the status `:diverged` when the iterates stop being finite or do not converge within the limit."
  (:require [fastmath.core :as m]
            [fastmath.solver.common :as c])
  (:import [clojure.lang IFn]
           [fastmath.solver.common IEvaluator]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- plain
  ^double [^IEvaluator ev ^double x0 ^double atol ^double rtol]
  (loop [x x0]
    (let [x1 (.eval ev x)]
      (when-not (m/valid-double? x1)
        (c/throw-status :diverged (str "The fixed point iteration is not finite: g(" x ") = " x1) {:evaluations (.evaluations ev) :x x}))
      (if (m/<= (m/abs (m/- x1 x)) (c/tolerance atol rtol x1))
        x1
        (recur x1)))))

(defn- steffensen
  ^double [^IEvaluator ev ^double x0 ^double atol ^double rtol]
  (loop [x x0]
    (let [x1 (.eval ev x)
          x2 (.eval ev x1)
          _ (when-not (and (m/valid-double? x1) (m/valid-double? x2))
              (c/throw-status :diverged (str "The fixed point iteration is not finite near " x) {:evaluations (.evaluations ev) :x x}))
          denom (m/+ (m/- x2 (m/* 2.0 x1)) x)
          xn (if (m/zero? denom) x2 (m/- x (m// (m/* (m/- x1 x) (m/- x1 x)) denom)))]
      (if (m/<= (m/abs (m/- xn x)) (c/tolerance atol rtol xn))
        xn
        (recur xn)))))

(defn fixed-point
  "Finds a fixed point of `g`, a solution of `g(x) = x`, by iteration from `x0`.

  Parameters:

  - `g` (function): a function of one number returning a number.
  - `x0` (finite number): the start.
  - `options` (map, optional):
      - `:steffensen?` - when truthy, uses the Steffensen acceleration; each step costs two evaluations. Default `false`.
      - `:absolute-accuracy` - absolute tolerance, default `1.0e-300`.
      - `:relative-accuracy` - tolerance relative to the result, default four times the machine epsilon. The iteration stops when the step is within `absolute + relative * |x|`.
      - `:max-iters` - limit of evaluations of `g` (an integer), default `1200`. Exceeding it throws.
      - `:info?` - when truthy, returns the map `{:root :evaluations :status}` instead of the fixed point, and returns the failure status instead of throwing. Default `false`.

  Returns the fixed point as a double. The step size is the only stopping rule: a small step of a slowly converging iteration (a rate close to 1) does not mean a small error. The error is about `step * rate / (1 - rate)`.

  Throws `ex-info` for an invalid argument, and for a failed iteration: iterates which are not finite (`:status :diverged`), `NaN` values of `g` (`:nan`), a non number result (`:invalid-result`) or the evaluation limit (`:max-evals`, which is also the result of a diverging or oscillating iteration). The data of the exception holds the `:status` and a few scalars.

  See also [[fastmath.solver/find-root]], [[fastmath.solver.newton/newton]]."
  ([g x0] (fixed-point g x0 nil))
  ([g x0 options]
   (let [{:keys [^double atol ^double rtol ^long max-evals info?]} (c/parse-options options)]
     (when-not (and (number? x0) (m/valid-double? x0))
       (c/throw-invalid-option :x0 x0 "must be a finite number"))
     (when-not (instance? IFn g)
       (c/throw-invalid-option :g g "must be a function"))
     (let [^IEvaluator ev (c/evaluator g max-evals)
           x0 (double x0)
           steffensen? (boolean (:steffensen? options))]
       (c/call-guarded
        info? ev
        (fn []
          (let [root (if steffensen? (steffensen ev x0 atol rtol) (plain ev x0 atol rtol))]
            (if info? {:root root :evaluations (.evaluations ev) :status :converged} root))))))))
