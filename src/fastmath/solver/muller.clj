(ns fastmath.solver.muller
  "Muller's method for a root of a function on an interval with a sign change, safeguarded by the bracket.

  Muller's method fits a parabola through the last three points and moves to its root, which gives a convergence of order about 1.84 and works for a root of even multiplicity or a pair of close roots where the secant method fails. Here the new point is always inside the bracket: a parabola without a root in the bracket, or a step which does not shrink the bracket fast enough, is replaced by a bisection step."
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
        start (.start p)
        s (if (m/nan? start) (m/+ (m/* 0.5 lo) (m/* 0.5 hi)) start)
        fs (.eval ev s)]
    (if (m/zero? fs)
      s
      (let [same? (m/== (m/signum fs) (m/signum fa0))
            a (if same? s lo)
            fa (if same? fs fa0)
            b (if same? hi s)
            fb (if same? fb0 fs)]
        ;; the window is (pa, pb, pc): pc is the newest point and the center of the parabola
        (loop [a a fa fa b b fb fb pa lo fpa fa0 pb hi fpb fb0 pc s fpc fs xprev ##Inf slow (long 0)]
          (let [inside? (fn [^double x] (and (m/> x a) (m/< x b) (not (m/== x pa)) (not (m/== x pb)) (not (m/== x pc))))
                ^Vec4 r (c/muller-step (Vec3. pa pc pb) (Vec3. fpa fpc fpb))
                x1 (.x r)
                x2 (.y r)
                ok1? (inside? x1)
                ok2? (inside? x2)
                cand (double (cond (and ok1? ok2?) (if (m/< (m/abs (m/- x1 pc)) (m/abs (m/- x2 pc))) x1 x2)
                                   ok1? x1
                                   ok2? x2
                                   :else ##NaN))
                bisect? (or (m/nan? cand) (m/>= slow 2))
                x (if bisect? (m/+ (m/* 0.5 a) (m/* 0.5 b)) cand)
                fx (.eval ev x)]
            (if (m/zero? fx)
              x
              (let [same? (m/== (m/signum fx) (m/signum fa))
                    na (if same? x a)
                    nfa (if same? fx fa)
                    nb (if same? b x)
                    nfb (if same? fb fx)
                    wold (m/- b a)
                    wnew (m/- nb na)
                    nslow (if (or bisect? (m/<= wnew (m/* 0.5 wold))) 0 (m/inc slow))
                    tol (c/tolerance atol rtol x)]
                (if (or (m/<= (m/* 0.5 wnew) tol)
                        (and (m/<= (m/abs (m/- x xprev)) tol) (c/certified? ev x fx p)))
                  x
                  (recur na nfa nb nfb pb fpb pc fpc x fx x nslow))))))))))

(defn muller
  "Finds a root of `f` on `[lower, upper]` with Muller's method safeguarded by the bracket.

  Parameters:

  - `f` (function): a function of one number returning a number; it has to change the sign on the interval.
  - `lower`, `upper` (finite numbers): ends of the interval; reversed ends are swapped.
  - `options` (map, optional):
      - `:absolute-accuracy` - absolute tolerance of the root, default `1.0e-300`.
      - `:relative-accuracy` - tolerance relative to the root, default four times the machine epsilon. The root is found within `absolute + relative * |root|`.
      - `:max-iters` - limit of function evaluations (an integer), default `1200`. Exceeding it throws.
      - `:initial-value` - the third point of the first parabola, strictly inside the interval, default the middle of the interval.
      - `:info?` - when truthy, returns the map `{:root :evaluations :status}` instead of the root, and returns the failure status instead of throwing. Default `false`.

  Returns the root as a double. A root at an end of the interval is returned after at most two evaluations. A negative discriminant or a parabola without a root inside the bracket never produces `NaN`: a bisection step is made instead.

  Throws `ex-info` for an invalid argument, and for a failed solve: no sign change on the interval (`:status :no-bracket`), `NaN` function values (`:nan`), a non number result (`:invalid-result`) or the evaluation limit (`:max-evals`). The data of the exception holds the `:status` and a few scalars.

  See also [[fastmath.solver/find-root]], [[fastmath.solver.muller2/muller2]], [[fastmath.solver.brent/brent]]."
  ([f lower upper] (muller f lower upper nil))
  ([f lower upper options]
   (c/solve-bracketed core true f lower upper options)))
