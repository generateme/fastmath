(ns fastmath.solver.common
  "Internal helpers shared by the root finders `fastmath.solver.*`.

  The namespace is not part of the public API and may change without notice. It collects what every solver needs, so that each solver namespace contains only its iteration.

  Options: [[parse-options]] validates the options of [[fastmath.solver/find-root]] and returns the setup of a solve; [[throw-invalid-option]] reports a rejected value.

  Evaluation: [[evaluator]] wraps the user function, counts the evaluations, enforces the limit `:max-iters` and rejects `NaN` and non numbers; [[tolerance]] is the rule `abs + rel * |x|` which decides convergence.

  Driver: [[solve-bracketed]] checks the bracket, builds the [[Problem]] and calls the iteration of a solver; [[certified?]] confirms that a point is a root of a sign change within the tolerance.

  Failures: a solve which does not converge throws `ex-info`; the exception data holds `:status` (one of [[statuses]]), `:evaluations` and a few scalars, never the iterates. With the option `:info?` the status is returned as a map instead."
  (:require [fastmath.core :as m]
            [fastmath.calculus.finite :as finite]
            [fastmath.vector])
  (:import [fastmath.java Array]
           [fastmath.vector Vec3 Vec4]
           [clojure.lang ExceptionInfo IFn IFn$DD]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^{:doc "Statuses of a solve which did not end with a root: `:max-evals` (the limit of evaluations is reached), `:no-bracket` (no sign change), `:nan` (the function is `NaN`), `:diverged` (the iteration left the interval, stalled or the result is not a root), `:invalid-result` (the function returned a non number). A solve which found the root has the status `:converged`."}
  statuses #{:converged :max-evals :no-bracket :nan :diverged :invalid-result})

(def ^{:const true :tag 'long :doc "Default limit of function evaluations: enough for the bisection to reach the floating point resolution on any interval."}
  default-max-iters 1200)

(def ^{:const true :tag 'double :doc "Default relative tolerance, four times the machine epsilon."}
  default-relative-accuracy (* 4.0 m/MACHINE-EPSILON2))

(def ^{:const true :tag 'double :doc "Default absolute tolerance, close to the smallest normal double so that it is relevant only for a root at zero."}
  default-absolute-accuracy 1.0e-300)

;; errors

(defn throw-invalid-option
  "Throws `ex-info` reporting an invalid value of an option or an argument.

  Parameters:

  - `option` (keyword): name of the option, for example `:max-iters`.
  - `value`: the rejected value.
  - `reason` (string): what the value has to satisfy.

  The exception data is `{:option option :value value :reason reason}`. Never returns."
  [option value reason]
  (throw (ex-info (str "Invalid option " option ": " reason) {:option option :value value :reason reason})))

(definterface SolverFailure)

(defn throw-status
  "Throws `ex-info` for a solve which ended without a root.

  Parameters:

  - `status` (keyword): one of [[statuses]].
  - `message` (string): the description.
  - `data` (map): scalars which help to locate the problem (an evaluation count, a point), never a sequence of iterates.

  The exception data is `data` with `:status` added; the exception is an `ex-info` marked as a solver failure, which [[call-guarded]] turns into an info map. Never returns."
  [status message data]
  (throw (proxy [ExceptionInfo SolverFailure] [message (assoc data :status status)])))

;; options

(defn- accuracy-option
  ^double [options option ^double default]
  (let [v (get options option)]
    (cond
      (nil? v) default
      (and (number? v) (let [d (double v)] (and (m/not-neg? d) (not (m/inf? d))))) (double v)
      :else (throw-invalid-option option v "must be a finite number not less than 0"))))

(defn long-option
  "Reads the option `option` of the map `options`: an integer not less than `1`, or `default` when the option is missing. Throws `ex-info` (see [[throw-invalid-option]]) for any other value."
  ^long [options option ^long default]
  (let [v (get options option)]
    (cond
      (nil? v) default
      (and (or (instance? Long v) (instance? Integer v) (instance? Short v) (instance? Byte v)) (>= (long v) 1)) (long v)
      :else (throw-invalid-option option v "must be an integer not less than 1"))))

(defn finite-number
  "Returns `v` as a double when it is a finite number, throws `ex-info` (see [[throw-invalid-option]]) for the option `option` otherwise."
  ^double [option v]
  (if (and (number? v) (m/valid-double? v))
    (double v)
    (throw-invalid-option option v "must be a finite number")))

(defn parse-options
  "Validates the options of a solve.

  Parameters:

  - `options` (map or `nil`): `:absolute-accuracy`, `:relative-accuracy` (finite numbers not less than `0`), `:max-iters` (integer not less than `1`, the limit of function evaluations), `:initial-value` (a finite number), `:info?` (truthy to return the info map). Unknown keys are ignored.

  Returns a map with `:atol`, `:rtol` (doubles), `:max-evals` (long), `:start` (a double, `NaN` when `:initial-value` is not given) and `:info?` (boolean). Throws `ex-info` (see [[throw-invalid-option]]) for an invalid value.

  See also [[tolerance]]."
  [options]
  (when-not (or (nil? options) (map? options))
    (throw-invalid-option :options options "must be a map or nil"))
  (let [start (get options :initial-value)]
    (when (and (some? start) (not (and (number? start) (m/valid-double? start))))
      (throw-invalid-option :initial-value start "must be a finite number"))
    {:atol (accuracy-option options :absolute-accuracy default-absolute-accuracy)
     :rtol (accuracy-option options :relative-accuracy default-relative-accuracy)
     :max-evals (long-option options :max-iters default-max-iters)
     :start (if (nil? start) ##NaN (double start))
     :info? (boolean (:info? options))}))

(defn tolerance
  "Tolerance of the root at `x`: `atol + rtol * |x|`, never less than two units in the last place of `x`, so that the iteration can always end."
  ^double [^double atol ^double rtol ^double x]
  (m/max (m/+ atol (m/* rtol (m/abs x))) (m/* 2.0 (m/ulp x))))

;; evaluation

(definterface IEvaluator
  (^double eval [^double x])
  (^long evaluations []))

(defn evaluator
  "Wraps the function `f` for a solver.

  Parameters:

  - `f` (function): the function of one number.
  - `max-evals` (long): the limit of evaluations.

  Returns an `IEvaluator` with `(.eval ev x)` and `(.evaluations ev)`. The evaluation throws `ex-info` with the status `:max-evals` when the limit is reached, `:nan` when `f` returns `NaN` and `:invalid-result` when `f` does not return a number. Infinite values are returned as they are.

  Contract: `f` is called only with finite doubles."
  ^fastmath.solver.common.IEvaluator [f ^long max-evals]
  (let [^IFn f f
        n (long-array 1)]
    (reify IEvaluator
      (evaluations [_] (Array/get n 0))
      (eval [_ x]
        (let [k (Array/get n 0)]
          (when (>= k max-evals)
            (throw-status :max-evals (str "The limit of " max-evals " function evaluations (:max-iters) is reached")
                          {:evaluations k :limit max-evals}))
          (Array/set n 0 (inc k))
          (let [y (.invoke f x)]
            (if (number? y)
              (let [d (double y)]
                (if (m/nan? d)
                  (throw-status :nan (str "The function is NaN at " x) {:evaluations (inc k) :x x})
                  d))
              (throw-status :invalid-result (str "The function returned a non number at " x ": " (pr-str y))
                            {:evaluations (inc k) :x x}))))))))

;; driver

(deftype Problem [^IEvaluator ev ^double a ^double fa ^double b ^double fb ^double atol ^double rtol ^double start])

(defn certified?
  "Checks that `x` is a root of a sign change within the tolerance.

  Returns `true` when `fx` is `0` or when the function has the opposite sign at `x - 8 tol` or at `x + 8 tol` (a probe outside the bracket `[a, b]` of the problem is skipped); costs up to two evaluations. A point which converged by the step size but is not a root (an even multiplicity, a jump) fails.

  Parameters: `ev` the evaluator, `x` and `fx` the point and the value, `p` the [[Problem]] with the bracket and the tolerance."
  [^IEvaluator ev ^double x ^double fx ^Problem p]
  (or (m/zero? fx)
      (let [d (m/* 8.0 (tolerance (.atol p) (.rtol p) x))
            xl (m/- x d)
            xr (m/+ x d)]
        (boolean
         (or (and (m/>= xl (.a p)) (m/< xl x) (not (m/== (m/signum fx) (m/signum (.eval ev xl)))))
             (and (m/<= xr (.b p)) (m/> xr x) (not (m/== (m/signum fx) (m/signum (.eval ev xr))))))))))

(defn- info-map
  [root ^IEvaluator ev status]
  {:root root :evaluations (.evaluations ev) :status status})

(defn call-guarded
  "Calls the function of no arguments `thunk` of a solve with the evaluator `ev`.

  A solver failure (see [[throw-status]]) is returned as the info map `{:root ##NaN :evaluations :status}` when `info?` is truthy and rethrown otherwise. Any other exception, also an `ex-info` thrown by the user function, is rethrown."
  [info? ^IEvaluator ev thunk]
  (try
    (thunk)
    (catch ExceptionInfo e
      (if (and info? (instance? SolverFailure e))
        (info-map ##NaN ev (:status (ex-data e)))
        (throw e)))))

(defn solve-bracketed
  "Runs the iteration `core` of a solver on a sign change of `f`.

  Parameters:

  - `core` (function): `(core ^Problem p)` returns the root as a double. It reads the evaluator, the bracket `[a, b]` with the values `fa`, `fb` (opposite signs, `a < b`), the tolerance `atol`, `rtol` and `start` (`NaN` when not given). It throws with [[throw-status]] when it does not converge.
  - `uses-start?` (boolean): when `true`, `:initial-value` has to lie strictly inside the bracket; when `false` it is ignored.
  - `f`, `lower`, `upper`, `options`: the arguments of the solver.

  Bounds are swapped when reversed. A root at an end of the interval is returned after at most two evaluations. Equal bounds return the point when `f` is zero there. Returns the root, or with `:info?` the map `{:root :evaluations :status}`; the status of a failed solve is returned instead of thrown.

  Throws `ex-info` for an invalid argument (see [[throw-invalid-option]]), and for a failed solve (see [[throw-status]]) when `:info?` is not set."
  [core uses-start? f lower upper options]
  (let [{:keys [^double atol ^double rtol ^long max-evals ^double start info?]} (parse-options options)
        lo (finite-number :lower-bound lower)
        hi (finite-number :upper-bound upper)
        a (m/min lo hi)
        b (m/max lo hi)]
    (when-not (instance? IFn f)
      (throw-invalid-option :f f "must be a function"))
    (when (and uses-start? (not (m/nan? start)) (not (and (m/> start a) (m/< start b))))
      (throw-invalid-option :initial-value start "must lie strictly inside the bracket"))
    (let [^IEvaluator ev (evaluator f max-evals)]
      (call-guarded
       info? ev
       (fn []
         (let [fa (.eval ev a)]
           (if (m/zero? fa)
             (if info? (info-map a ev :converged) a)
             (let [fb (if (m/== a b) fa (.eval ev b))]
               (cond
                 (m/zero? fb) (if info? (info-map b ev :converged) b)
                 (= (m/signum fa) (m/signum fb))
                 (throw-status :no-bracket (str "The function does not change the sign on [" a ", " b "]: f(a)=" fa ", f(b)=" fb)
                               {:evaluations (.evaluations ev) :a a :b b :fa fa :fb fb})
                 :else (let [root (double (core (Problem. ev a fa b fb atol rtol start)))]
                         (if info? (info-map root ev :converged) root)))))))))))

;; false position family

(defn false-position
  "Iteration of the false position family on the bracket of the [[Problem]] `p`; returns the root.

  `variant` selects what happens to the retained end of the bracket when the new point has the same sign as the newest one: `:regula-falsi` keeps its value, `:illinois` halves it and `:pegasus` scales it by `f1 / (f1 + f)`.

  The search ends when the bracket is narrower than twice the tolerance, or when the step is within the tolerance and the point is certified (see [[certified?]]). The new point is always strictly inside the bracket, a step which could not move is replaced by the middle."
  ^double [variant ^Problem p]
  (let [^IEvaluator ev (.ev p)
        atol (.atol p)
        rtol (.rtol p)
        lo (.a p)
        hi (.b p)
        illinois? (identical? variant :illinois)
        pegasus? (identical? variant :pegasus)]
    (loop [x0 lo f0 (.fa p) x1 hi f1 (.fb p)]
      (let [scale (m/max (m/abs f0) (m/abs f1))
            r0 (m// f0 scale)
            r1 (m// f1 scale)
            t (m// r1 (m/- r1 r0))
            xi (m/+ (m/* (m/- 1.0 t) x1) (m/* t x0))
            xl (m/min x0 x1)
            xh (m/max x0 x1)
            x (if (and (m/> xi xl) (m/< xi xh)) xi (m/+ (m/* 0.5 x0) (m/* 0.5 x1)))
            fx (.eval ev x)]
        (if (m/zero? fx)
          x
          (let [opposite? (not (m/== (m/signum fx) (m/signum f1)))
                nf0 (double (cond opposite? f1
                                  illinois? (m/* 0.5 f0)
                                  pegasus? (m/* f0 (m// f1 (m/+ f1 fx)))
                                  :else f0))
                nx0 (double (if opposite? x1 x0))
                tol (tolerance atol rtol x)]
            (cond
              (m/<= (m/abs (m/- (m/* 0.5 x) (m/* 0.5 nx0))) tol) x
              (and (m/<= (m/abs (m/- x x1)) tol) (certified? ev x fx p)) x
              :else (recur nx0 nf0 x fx))))))))

;; Muller

(defn muller-step
  "Muller's quadratic through three points, centered at the middle one.

  Parameters: `xs` (`Vec3`) the abscissae `[xa xc xb]`, the center `xc` is the middle one, and `ys` (`Vec3`) the values `[ya yc yb]`. The abscissae have to be distinct.

  Returns a `Vec4` `[x y z w]`: `x` and `y` are the two real roots of the parabola (`NaN` when the discriminant is negative, `x` is the root with the larger denominator, nearer to `xc`), `z` the slope of the parabola at `xc` and `w` its curvature."
  ^Vec4 [^Vec3 xs ^Vec3 ys]
  (let [xa (.x xs) xc (.y xs) xb (.z xs)
        ya (.x ys) yc (.y ys) yb (.z ys)
        dac (m// (m/- yc ya) (m/- xc xa))
        dcb (m// (m/- yb yc) (m/- xb xc))
        d (m// (m/- dcb dac) (m/- xb xa))
        c (m/+ dac (m/* (m/- xc xa) d))
        disc (m/- (m/* c c) (m/* 4.0 yc d))
        sq (m/sqrt disc)
        dp (m/+ c sq)
        dm (m/- c sq)
        big (if (m/> (m/abs dp) (m/abs dm)) dp dm)
        small (if (m/> (m/abs dp) (m/abs dm)) dm dp)]
    (Vec4. (m/- xc (m// (m/* 2.0 yc) big)) (m/- xc (m// (m/* 2.0 yc) small)) c d)))

;; Newton and Halley

(defn derivative-evaluator
  "Wraps the derivative `df` (a function of one number) for a solver.

  Returns a function of `x` returning the derivative as a double; `NaN` and infinite derivatives are returned as they are (the caller replaces the step by a bisection step). Throws `ex-info` with the status `:invalid-result` when `df` does not return a number."
  [df]
  (let [^IFn df df]
    (fn ^double [^double x]
      (let [y (.invoke df x)]
        (if (number? y)
          (double y)
          (throw-status :invalid-result (str "The derivative returned a non number at " x ": " (pr-str y)) {:x x}))))))

(defn newton-iteration
  "Iteration of the Newton (`:newton`) or Halley (`:halley`) method on the bracket of the [[Problem]] `p`; returns the root.

  `df` and `d2f` are functions of `x` returning the first and the second derivative as doubles (`d2f` is ignored by Newton).

  The iteration is safeguarded by the bracket: a step which would leave the bracket, is not finite or does not converge at least as fast as halving (the safeguard of the bisection) is replaced by the bisection step, so the iteration cannot fail on a function with a sign change. The start is the `start` of the problem or the middle of the bracket. The search ends when the bracket is narrower than twice the tolerance, or when the step is within the tolerance and the point is certified (see [[certified?]])."
  ^double [variant df d2f ^Problem p]
  (let [^IEvaluator ev (.ev p)
        atol (.atol p)
        rtol (.rtol p)
        lo (.a p)
        hi (.b p)
        start (.start p)
        halley? (identical? variant :halley)
        ^IFn$DD df df
        ^IFn$DD d2f d2f]
    (loop [a lo fa (.fa p) b hi fb (.fb p)
           x (if (m/nan? start) (m/+ (m/* 0.5 lo) (m/* 0.5 hi)) start)
           dx-prev ##Inf
           dx-old ##Inf]
      (let [fx (.eval ev x)]
        (if (m/zero? fx)
          x
          (let [same? (m/== (m/signum fx) (m/signum fa))
                na (if same? x a)
                nfa (if same? fx fa)
                nb (if same? b x)
                nfb (if same? fb fx)
                d1 (.invokePrim df x)
                dx (if halley?
                     (let [d2 (.invokePrim d2f x)]
                       (m// (m/* 2.0 fx d1) (m/- (m/* 2.0 d1 d1) (m/* fx d2))))
                     (m// fx d1))
                xn (m/- x dx)
                slow? (m/> (m/abs (m/* 2.0 fx)) (m/abs (m/* dx-old d1)))
                use-newton? (and (m/valid-double? dx) (m/> xn na) (m/< xn nb) (not slow?))
                x2 (if use-newton? xn (m/+ (m/* 0.5 na) (m/* 0.5 nb)))
                tol (tolerance atol rtol x2)
                step (m/- x x2)]
            (if (or (m/<= (m/* 0.5 (m/- nb na)) tol)
                    (and (m/<= (m/abs step) tol)
                         (let [f2 (.eval ev x2)]
                           (or (m/zero? f2) (certified? ev x2 f2 p)))))
              x2
              (recur na nfa nb nfb x2 step dx-prev))))))))

(defn- derivative-function
  "The derivative of the order `n` as a function `x -> double`: the user function `d`, or when `d` is `nil` the central finite difference of the accuracy order `acc` of the counted function `ev`."
  [d n acc ^IEvaluator ev]
  (if (nil? d)
    (finite/derivative (fn ^double [^double x] (.eval ev x)) n {:acc acc})
    (derivative-evaluator d)))

(defn solve-with-derivatives
  "Solves with the safeguarded Newton (`:newton`) or Halley (`:halley`) iteration, see [[newton-iteration]].

  Parameters: `variant`, the function `f`, the bounds and the `options` of the solver; the options `:df` and `:d2f` are the first and the second derivative of `f` (functions of one number) and default to central finite differences (order 4 and 2) of `f`, whose evaluations count to `:max-iters`.

  Throws `ex-info` (see [[throw-invalid-option]]) when `:df` or `:d2f` is given and is not a function, before any evaluation of `f`."
  [variant f lower upper options]
  (let [df (:df options)
        d2f (:d2f options)]
    (when-not (or (nil? df) (instance? IFn df)) (throw-invalid-option :df df "must be a function"))
    (when-not (or (nil? d2f) (instance? IFn d2f)) (throw-invalid-option :d2f d2f "must be a function"))
    (solve-bracketed (fn [^Problem p]
                       (let [^IEvaluator ev (.ev p)]
                         (newton-iteration variant (derivative-function df 1 4 ev) (derivative-function d2f 2 2 ev) p)))
                     true f lower upper options)))
