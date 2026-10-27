(ns fastmath.solver
  "Root finding of univariate functions and polynomials.

  Bracketing methods find a root of a function that changes the sign on an interval and cannot fail on it: [[brent]] (default), [[bisection]], [[illinois]], [[pegasus]], [[regula-falsi]], [[ridders]], [[muller]] and [[toms748]]. [[secant]] and [[muller2]] have no bracket: they are faster on smooth functions but fail with the status `:diverged` when the iteration runs away. [[newton]] and [[halley]] use derivatives (given or approximated by finite differences) and are safeguarded by the bracket. [[fixed-point]] solves `g(x) = x`.

  [[find-root]] is the common entry point of all of them, [[find-roots]] finds all the roots in an interval and [[find-bracket]] finds an interval with a sign change around a point. Every solver lives in its own namespace `fastmath.solver.<name>` and has the same signature `(solve f lower upper options)`; the options are described in [[find-root]].

  [[quadratic]] and [[cubic]] return the real roots of the polynomials of degree two and three, [[polyroot]] all the complex roots of a polynomial of any degree.

  Failures are reported with `ex-info` whose data holds the `:status` (`:no-bracket`, `:max-evals`, `:nan`, `:diverged`, `:invalid-result`) and a few scalars; the option `:info?` returns the status as a map instead."
  (:require [fastmath.core :as m]
            [fastmath.solver.common :as c]
            [fastmath.solver.brent :as brent-ns]
            [fastmath.solver.bisection :as bisection-ns]
            [fastmath.solver.illinois :as illinois-ns]
            [fastmath.solver.pegasus :as pegasus-ns]
            [fastmath.solver.regula-falsi :as regula-falsi-ns]
            [fastmath.solver.ridders :as ridders-ns]
            [fastmath.solver.muller :as muller-ns]
            [fastmath.solver.muller2 :as muller2-ns]
            [fastmath.solver.secant :as secant-ns]
            [fastmath.solver.toms748 :as toms748-ns]
            [fastmath.solver.newton :as newton-ns]
            [fastmath.solver.halley :as halley-ns]
            [fastmath.solver.fixed-point :as fixed-point-ns]
            [fastmath.solver.closed-form :as closed-form-ns]
            [fastmath.solver.polyroot :as polyroot-ns])
  (:import [clojure.lang IFn]
           [fastmath.solver.common IEvaluator]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defmacro ^:private export
  "Defines the var of the same name in this namespace with the value and the documentation of `v`."
  [v]
  (let [mv (meta (resolve v))]
    `(def ~(with-meta (symbol (name v)) {:doc (:doc mv) :arglists (list 'quote (:arglists mv))}) ~v)))

(export brent-ns/brent)
(export bisection-ns/bisection)
(export illinois-ns/illinois)
(export pegasus-ns/pegasus)
(export regula-falsi-ns/regula-falsi)
(export ridders-ns/ridders)
(export muller-ns/muller)
(export muller2-ns/muller2)
(export secant-ns/secant)
(export toms748-ns/toms748)
(export newton-ns/newton)
(export halley-ns/halley)
(export fixed-point-ns/fixed-point)
(export closed-form-ns/quadratic)
(export closed-form-ns/cubic)
(export polyroot-ns/polyroot)

#_{:clj-kondo/ignore [:unresolved-symbol]}
(def ^{:doc "Map from the keyword of the option `:solver` of [[find-root]] to the solver function `(solve f lower upper options)`."}
  solvers
  {:brent brent
   :bisection bisection
   :illinois illinois
   :pegasus pegasus
   :regula-falsi regula-falsi
   :ridders ridders
   :muller muller
   :muller2 muller2
   :secant secant
   :toms748 toms748
   :newton newton
   :halley halley})

(defn- solver-fn
  [options]
  (let [s (or (get options :solver) :brent)]
    (or (get solvers s)
        (c/throw-invalid-option :solver s (str "must be one of " (pr-str (vec (sort (keys solvers)))))))))

(defn find-root
  "Finds a root of the function `f` on the interval `[lower, upper]`.

  Parameters:

  - `f` (function): a function of one number returning a number; it has to change the sign on the interval.
  - `lower`, `upper` (finite numbers): ends of the interval; reversed ends are swapped.
  - `options` (map, optional):
      - `:solver` - the method, one of `:brent` (default), `:bisection`, `:illinois`, `:pegasus`, `:regula-falsi`, `:ridders`, `:muller`, `:muller2`, `:secant`, `:toms748`, `:newton`, `:halley`; see the keys of [[solvers]]. An unknown method throws.
      - `:absolute-accuracy` - absolute tolerance of the root, default `1.0e-300`.
      - `:relative-accuracy` - tolerance relative to the root, default four times the machine epsilon. The root is found within `absolute + relative * |root|`, so the default is the resolution of doubles at any scale of the root.
      - `:max-iters` - limit of function evaluations (an integer), default `1200`. Exceeding it throws.
      - `:initial-value` - the first point of the search, strictly inside the interval; used by `:brent`, `:muller`, `:newton` and `:halley`, ignored by the other methods.
      - `:df`, `:d2f` - the first and the second derivative of `f`, used by `:newton` and `:halley`; the default is a central finite difference.
      - `:info?` - when truthy, returns the map `{:root :evaluations :status}` instead of the root, and returns the failure status instead of throwing. Default `false`.

  Returns the root as a double. A root at an end of the interval is returned after at most two evaluations. The function value is never used to stop the search, so a function with tiny values is solved as accurately as any other.

  Throws `ex-info` for an invalid argument or option, and for a failed solve: no sign change on the interval (`:status :no-bracket`), `NaN` function values (`:nan`), a non number result (`:invalid-result`), the evaluation limit (`:max-evals`) or a runaway iteration of `:secant` and `:muller2` (`:diverged`). The data of the exception holds the `:status`, the number of `:evaluations` and a few scalars.

  See also [[find-roots]], [[find-bracket]], [[solvers]], [[fixed-point]]."
  ([f lower upper] (find-root f lower upper nil))
  ([f lower upper options]
   ((solver-fn options) f lower upper options)))

;; find-bracket

(def ^{:private true :const true :tag 'double} MAX-DOUBLE Double/MAX_VALUE)

(defn- clamp
  ^double [^double lo ^double hi ^double x]
  (m/min hi (m/max lo x)))

(defn- number-option
  "Reads the numeric option `option` of `options`: a number accepted by the predicate `valid?` (of a double), or `default` when missing; throws `ex-info` with `reason` otherwise. `rule` is the pair `[valid? reason]`."
  ^double [options option ^double default rule]
  (let [[valid? reason] rule
        v (get options option)]
    (cond
      (nil? v) default
      (and (number? v) (valid? (double v))) (double v)
      :else (c/throw-invalid-option option v reason))))

(defn find-bracket
  "Finds an interval with a sign change of `f` around the point `x0` by expanding it geometrically.

  Parameters:

  - `f` (function): a function of one number returning a number.
  - `x0` (finite number): the point.
  - `options` (map, optional):
      - `:step` - half-width of the first interval, a positive number, default `0.1 * max(1, |x0|)`.
      - `:factor` - the growth of the interval in every step, greater than `1`, default `1.6`. In every step the end with the smaller absolute value of `f` is moved outwards.
      - `:lower`, `:upper` - limits of the search, default none (`##-Inf`, `##Inf`). An end which reaches a limit stays there.
      - `:max-iters` - limit of function evaluations, default `200`.

  Returns the vector `[a b]` with `a < b` and a sign change (or a zero) between `f(a)` and `f(b)`, ready for [[find-root]].

  Throws `ex-info` for an invalid argument, with the status `:no-bracket` when both ends reached the limits `:lower` and `:upper` without a sign change, `:max-evals` when the number of evaluations is exhausted, `:nan` when `f` is `NaN` and `:invalid-result` when it is not a number. The data holds the `:status`, the `:evaluations` and the last interval `:a`, `:b`.

  See also [[find-root]], [[find-roots]]."
  ([f x0] (find-bracket f x0 nil))
  ([f x0 options]
   (when-not (instance? IFn f) (c/throw-invalid-option :f f "must be a function"))
   (let [x0 (c/finite-number :x0 x0)
         step (number-option options :step (m/* 0.1 (m/max 1.0 (m/abs x0)))
                             [(fn [^double d] (and (m/pos? d) (m/valid-double? d))) "must be a positive finite number"])
         factor (number-option options :factor 1.6
                               [(fn [^double d] (and (m/> d 1.0) (m/valid-double? d))) "must be a finite number greater than 1"])
         ;; the search never leaves the finite doubles, f is called only with finite points
         lo-limit (m/max (m/- MAX-DOUBLE) (number-option options :lower ##-Inf [(fn [^double d] (not (m/nan? d))) "must be a number"]))
         hi-limit (m/min MAX-DOUBLE (number-option options :upper ##Inf [(fn [^double d] (not (m/nan? d))) "must be a number"]))
         max-evals (c/long-option options :max-iters 200)
         _ (when (m/> lo-limit x0) (c/throw-invalid-option :lower lo-limit "the point has to lie between :lower and :upper"))
         _ (when (m/< hi-limit x0) (c/throw-invalid-option :upper hi-limit "the point has to lie between :lower and :upper"))
         ^IEvaluator ev (c/evaluator f max-evals)
         s1 (clamp lo-limit hi-limit (m/- x0 step))
         s2 (clamp lo-limit hi-limit (m/+ x0 step))]
     (loop [x1 s1 f1 (.eval ev s1) x2 s2 f2 (.eval ev s2)]
       (let [can1? (m/> x1 lo-limit)
             can2? (m/< x2 hi-limit)]
         (cond
           (or (m/zero? f1) (m/zero? f2) (not (m/== (m/signum f1) (m/signum f2))))
           [(m/min x1 x2) (m/max x1 x2)]
           (not (or can1? can2?))
           (c/throw-status :no-bracket "No sign change was found within the limits :lower and :upper"
                           {:evaluations (.evaluations ev) :a x1 :b x2})
           (and can1? (or (not can2?) (m/< (m/abs f1) (m/abs f2))))
           (let [nx (clamp lo-limit hi-limit (m/+ x1 (m/* factor (m/- x1 x2))))]
             (recur nx (.eval ev nx) x2 f2))
           :else
           (let [nx (clamp lo-limit hi-limit (m/+ x2 (m/* factor (m/- x2 x1))))]
             (recur x1 f1 nx (.eval ev nx)))))))))

;; find-roots

(defn find-roots
  "Finds all the roots of `f` on the interval `[lower, upper]` that are sign changes at the resolution of a grid.

  The interval is divided into `:samples` equal parts; every part with a sign change is solved by [[find-root]], and a grid point with `f = 0` is a root. A pole is a sign change too, but it is not a root: a solution where the absolute value of `f` is larger than at both ends of its part is dropped.

  Parameters:

  - `f` (function): a function of one number returning a number.
  - `lower`, `upper` (finite numbers): ends of the interval, `lower < upper`.
  - `options` (map, optional):
      - `:samples` - the number of parts, an integer not less than `1`, default `1000`.
      - all the options of [[find-root]]; `:max-iters` applies to every root separately, `:info?` and `:initial-value` are not supported and ignored. The default `:absolute-accuracy` is `1e-15 * max(|lower|, |upper|)`, so that the search of a pole, which attracts the solver, ends quickly.

  Returns a vector of the roots in ascending order, empty when there is none. A root inside a part where the function changes the sign twice (two close roots), and a root of an even multiplicity (a touch of zero without a sign change) are missed: raise `:samples` to resolve close roots. Parts where `f` is `NaN` at an end are skipped.

  Throws `ex-info` for an invalid argument or option, and the failures of [[find-root]] for a part (for example `:max-evals`).

  See also [[find-root]], [[find-bracket]]."
  ([f lower upper] (find-roots f lower upper nil))
  ([f lower upper options]
   (when-not (instance? IFn f) (c/throw-invalid-option :f f "must be a function"))
   (let [lo (c/finite-number :lower-bound lower)
         hi (c/finite-number :upper-bound upper)
         _ (when-not (m/< lo hi) (c/throw-invalid-option :upper-bound upper "must be greater than the lower bound"))
         n (c/long-option options :samples 1000)
         opts (-> (dissoc options :samples :info? :initial-value)
                  (update :absolute-accuracy #(if (nil? %) (m/* 1.0e-15 (m/max (m/abs lo) (m/abs hi))) %)))
         solve (solver-fn opts)
         value (fn ^double [^double x]
                 (let [y (f x)] (if (number? y) (double y) ##NaN)))
         point (fn ^double [^long i] (if (m/== i n) hi (m/+ lo (m/* (m// (double i) (double n)) (m/- hi lo)))))
         xs (mapv point (range (m/inc n)))
         ys (mapv value xs)]
     (->> (range n)
          (mapcat (fn [^long i]
                    (let [x0 (double (xs i)) x1 (double (xs (m/inc i)))
                          y0 (double (ys i)) y1 (double (ys (m/inc i)))]
                      (cond
                        (m/zero? y0) [x0]
                        (and (m/== i (m/dec n)) (m/zero? y1)) [x1]
                        (or (m/nan? y0) (m/nan? y1) (m/zero? y1)) nil
                        (not (m/== (m/signum y0) (m/signum y1)))
                        (let [r (double (solve f x0 x1 opts))
                              yr (m/abs (value r))]
                          (when-not (m/> yr (m/max (m/abs y0) (m/abs y1))) [r]))
                        :else nil))))
          vec))))
