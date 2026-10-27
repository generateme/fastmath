(ns fastmath.solver.brent
  "Brent's method for a root of a function on an interval with a sign change.

  The method combines bisection, the secant method and inverse quadratic interpolation. It always keeps a bracket of the root, converges at least as fast as bisection and usually superlinearly. It is the default method of [[fastmath.solver/find-root]]."
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
        start (.start p)
        a0 (.a p)
        b0 (.b p)
        fa0 (.fa p)
        fb0 (.fb p)
        ;; a given start replaces the end of the bracket which has the same sign
        fs (if (m/nan? start) ##NaN (.eval ev start))
        use-start? (not (m/nan? start))]
    (if (and use-start? (m/zero? fs))
      start
      (let [same? (and use-start? (m/== (m/signum fs) (m/signum fa0)))
            a (if (and use-start? same?) start a0)
            fa (if (and use-start? same?) fs fa0)
            b (if (and use-start? (not same?)) start b0)
            fb (if (and use-start? (not same?)) fs fb0)]
        ;; b is the best estimate, c the other end of the bracket, a the previous b
        (loop [a a fa fa b b fb fb c a fc fa]
          (let [prev-step (m/- b a)
                swap? (m/< (m/abs fc) (m/abs fb))
                a (if swap? b a)
                fa (if swap? fb fa)
                b (if swap? c b)
                fb (if swap? fc fb)
                c (if swap? a c)
                fc (if swap? fa fc)
                tol (c/tolerance atol rtol b)
                cb (m/- c b)
                new-step (m/- (m/* 0.5 c) (m/* 0.5 b))]
            (if (or (m/<= (m/abs new-step) tol) (m/zero? fb))
              b
              (let [step (if (and (m/>= (m/abs prev-step) tol) (m/> (m/abs fa) (m/abs fb)))
                           (let [linear? (m/== a c)
                                 t1 (m// fb (if linear? fa fc))
                                 p0 (if linear?
                                      (m/* cb t1)
                                      (let [q (m// fa fc)
                                            t2 (m// fb fa)]
                                        (m/* t2 (m/- (m/* cb q (m/- q t1)) (m/* (m/- b a) (m/- t1 1.0))))))
                                 q0 (if linear?
                                      (m/- 1.0 t1)
                                      (let [q (m// fa fc)
                                            t2 (m// fb fa)]
                                        (m/* (m/- q 1.0) (m/- t1 1.0) (m/- t2 1.0))))
                                 p (if (m/pos? p0) p0 (m/- p0))
                                 q (if (m/pos? p0) (m/- q0) q0)]
                             (if (and (m/< p (m/- (m/* 0.75 cb q) (m/* 0.5 (m/abs (m/* tol q)))))
                                      (m/< p (m/abs (m/* 0.5 prev-step q))))
                               (m// p q)
                               new-step))
                           new-step)
                    step (if (m/< (m/abs step) tol) (m/copy-sign tol step) step)
                    nb (m/+ b step)
                    fnb (.eval ev nb)
                    same-c? (m/== (m/signum fnb) (m/signum fc))]
                (if same-c?
                  (recur b fb nb fnb b fb)
                  (recur b fb nb fnb c fc))))))))))

(defn brent
  "Finds a root of `f` on `[lower, upper]` with Brent's method.

  Parameters:

  - `f` (function): a function of one number returning a number; it has to change the sign on the interval.
  - `lower`, `upper` (finite numbers): ends of the interval; reversed ends are swapped.
  - `options` (map, optional):
      - `:absolute-accuracy` - absolute tolerance of the root, default `1.0e-300`.
      - `:relative-accuracy` - tolerance relative to the root, default four times the machine epsilon. The root is found within `absolute + relative * |root|`.
      - `:max-iters` - limit of function evaluations (an integer), default `1200`. Exceeding it throws.
      - `:initial-value` - the first point of the search, strictly inside the interval, default the plain bracket. The point replaces the end of the interval which has the same sign.
      - `:info?` - when truthy, returns the map `{:root :evaluations :status}` instead of the root, and returns the failure status instead of throwing. Default `false`.

  Returns the root as a double. A root at an end of the interval is returned after at most two evaluations. Function values are never used to stop the search: only the bracket decides, so a function with tiny values is solved as accurately as any other.

  Throws `ex-info` for an invalid argument, and for a failed solve: no sign change on the interval (`:status :no-bracket`), `NaN` function values (`:nan`), a non number result (`:invalid-result`) or the evaluation limit (`:max-evals`). The data of the exception holds the `:status` and a few scalars.

  See also [[fastmath.solver/find-root]], [[fastmath.solver.bisection/bisection]], [[fastmath.solver.toms748/toms748]]."
  ([f lower upper] (brent f lower upper nil))
  ([f lower upper options]
   (c/solve-bracketed core true f lower upper options)))
