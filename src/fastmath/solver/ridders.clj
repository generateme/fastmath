(ns fastmath.solver.ridders
  "Ridders' method for a root of a function on an interval with a sign change.

  Every step evaluates the function at the middle of the bracket and fits an exponential factor to the three values, which gives a quadratic convergence per two evaluations. The new point always stays inside the bracket."
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
    (loop [xl lo fl (.fa p) xh hi fh (.fb p) ans ##Inf]
      (let [xm (m/+ (m/* 0.5 xl) (m/* 0.5 xh))
            fm (.eval ev xm)]
        (if (m/zero? fm)
          xm
          (let [scale (m/max (m/abs fl) (m/abs fh) (m/abs fm))
                rl (m// fl scale)
                rh (m// fh scale)
                rm (m// fm scale)
                s (m/sqrt (m/- (m/* rm rm) (m/* rl rh)))
                xi (m/+ xm (m/* (m/- xm xl) (m// (m/* (m/signum (m/- fl fh)) rm) s)))
                xnew (if (and (m/>= xi (m/min xl xh)) (m/<= xi (m/max xl xh))) xi xm)
                fnew (.eval ev xnew)]
            (if (m/zero? fnew)
              xnew
              (let [sf (m/signum fnew)
                    split? (not (m/== sf (m/signum fm)))
                    low? (and (not split?) (not (m/== sf (m/signum fl))))
                    nxl (double (cond split? xm low? xl :else xnew))
                    nfl (double (cond split? fm low? fl :else fnew))
                    nxh (if (or split? low?) xnew xh)
                    nfh (if (or split? low?) fnew fh)
                    tol (c/tolerance atol rtol xnew)]
                (if (or (m/<= (m/abs (m/- (m/* 0.5 nxh) (m/* 0.5 nxl))) tol)
                        (and (m/<= (m/abs (m/- xnew ans)) tol) (c/certified? ev xnew fnew p)))
                  xnew
                  (recur nxl nfl nxh nfh xnew))))))))))

(defn ridders
  "Finds a root of `f` on `[lower, upper]` with Ridders' method.

  Parameters:

  - `f` (function): a function of one number returning a number; it has to change the sign on the interval.
  - `lower`, `upper` (finite numbers): ends of the interval; reversed ends are swapped.
  - `options` (map, optional):
      - `:absolute-accuracy` - absolute tolerance of the root, default `1.0e-300`.
      - `:relative-accuracy` - tolerance relative to the root, default four times the machine epsilon. The root is found within `absolute + relative * |root|`.
      - `:max-iters` - limit of function evaluations (an integer), default `1200`. Exceeding it throws.
      - `:initial-value` - ignored by this method.
      - `:info?` - when truthy, returns the map `{:root :evaluations :status}` instead of the root, and returns the failure status instead of throwing. Default `false`.

  Returns the root as a double. A root at an end of the interval is returned after at most two evaluations. Two function evaluations are made per step.

  Throws `ex-info` for an invalid argument, and for a failed solve: no sign change on the interval (`:status :no-bracket`), `NaN` function values (`:nan`), a non number result (`:invalid-result`) or the evaluation limit (`:max-evals`). The data of the exception holds the `:status` and a few scalars.

  See also [[fastmath.solver/find-root]], [[fastmath.solver.brent/brent]], [[fastmath.solver.illinois/illinois]]."
  ([f lower upper] (ridders f lower upper nil))
  ([f lower upper options]
   (c/solve-bracketed core false f lower upper options)))
