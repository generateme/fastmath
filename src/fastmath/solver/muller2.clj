(ns fastmath.solver.muller2
  "Muller's method for a root of a function on an interval, without a bracket.

  The classic Muller iteration: a parabola through the last three points, starting from the ends and the middle of the interval. The bracket is not kept, so the iteration is faster than the bracketing one but can run away: an iterate which leaves the interval, or a step which converged to a point which is not a root, throws. A parabola without real roots moves to its vertex, a flat parabola to a bisection step of the sign change seen so far."
  (:require [fastmath.core :as m]
            [fastmath.solver.common :as c])
  (:import [fastmath.solver.common IEvaluator Problem]
           [fastmath.vector Vec3 Vec4]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- core
  ^double [^Problem p]
  (let [^IEvaluator ev (.ev p)
        atol (.atol p)
        rtol (.rtol p)
        lo (.a p)
        hi (.b p)
        fa0 (.fa p)
        fb0 (.fb p)
        mid (m/+ (m/* 0.5 lo) (m/* 0.5 hi))
        fm (.eval ev mid)]
    (if (m/zero? fm)
      mid
      (let [same? (m/== (m/signum fm) (m/signum fa0))
            ;; the tightest sign change seen so far, the fallback of a useless step
            a (if same? mid lo)
            fa (if same? fm fa0)
            b (if same? hi mid)
            fb (if same? fb0 fm)]
        (loop [a a fa fa b b fb fb pa lo fpa fa0 pb hi fpb fb0 pc mid fpc fm xprev ##Inf]
          (let [^Vec4 r (c/muller-step (Vec3. pa pc pb) (Vec3. fpa fpc fpb))
                x1 (.x r)
                x2 (.y r)
                slope (.z r)
                curv (.w r)
                cand (double (cond (and (m/valid-double? x1) (m/valid-double? x2))
                                   (if (m/< (m/abs (m/- x1 pc)) (m/abs (m/- x2 pc))) x1 x2)
                                   (m/valid-double? x1) x1
                                   (m/valid-double? x2) x2
                                   ;; no real root: the vertex of the parabola, or a bisection step when it is flat
                                   (and (m/not-zero? curv) (m/valid-double? (m// slope curv))) (m/- pc (m// slope (m/* 2.0 curv)))
                                   :else ##NaN))
                x (if (or (m/nan? cand) (m/== cand pa) (m/== cand pb) (m/== cand pc))
                    (m/+ (m/* 0.5 a) (m/* 0.5 b))
                    cand)]
            (when-not (and (m/>= x lo) (m/<= x hi))
              (c/throw-status :diverged (str "Muller's method left the interval [" lo ", " hi "]: " x)
                              {:evaluations (.evaluations ev) :x x :a lo :b hi}))
            (let [fx (.eval ev x)]
              (if (m/zero? fx)
                x
                (let [same? (m/== (m/signum fx) (m/signum fa))
                      na (if (and same? (m/> x a) (m/< x b)) x a)
                      nfa (if (and same? (m/> x a) (m/< x b)) fx fa)
                      nb (if (and (not same?) (m/> x a) (m/< x b)) x b)
                      nfb (if (and (not same?) (m/> x a) (m/< x b)) fx fb)]
                  (if (m/<= (m/abs (m/- x xprev)) (c/tolerance atol rtol x))
                    (if (c/certified? ev x fx p)
                      x
                      (c/throw-status :diverged (str "Muller's method converged to a point which is not a root: " x)
                                      {:evaluations (.evaluations ev) :x x :fx fx}))
                    (recur na nfa nb nfb pb fpb pc fpc x fx x)))))))))))

(defn muller2
  "Finds a root of `f` on `[lower, upper]` with Muller's method, without a bracket.

  Parameters:

  - `f` (function): a function of one number returning a number; it has to change the sign on the interval.
  - `lower`, `upper` (finite numbers): ends of the interval; reversed ends are swapped. The interval is also the limit of the search.
  - `options` (map, optional):
      - `:absolute-accuracy` - absolute tolerance of the root, default `1.0e-300`.
      - `:relative-accuracy` - tolerance relative to the root, default four times the machine epsilon. The root is found within `absolute + relative * |root|`.
      - `:max-iters` - limit of function evaluations (an integer), default `1200`. Exceeding it throws.
      - `:initial-value` - ignored by this method.
      - `:info?` - when truthy, returns the map `{:root :evaluations :status}` instead of the root, and returns the failure status instead of throwing. Default `false`.

  Returns the root as a double. A root at an end of the interval is returned after at most two evaluations. The search stops when the step is within the tolerance and the function changes the sign within two tolerances of the point.

  Throws `ex-info` for an invalid argument, and for a failed solve: no sign change on the interval (`:status :no-bracket`), `NaN` function values (`:nan`), a non number result (`:invalid-result`), the evaluation limit (`:max-evals`) or a runaway iteration (`:diverged`): the iterate leaves the interval, or the step converged to a point which is not a root. The data of the exception holds the `:status` and a few scalars.

  See also [[fastmath.solver/find-root]], [[fastmath.solver.muller/muller]], [[fastmath.solver.secant/secant]]."
  ([f lower upper] (muller2 f lower upper nil))
  ([f lower upper options]
   (c/solve-bracketed core false f lower upper options)))
