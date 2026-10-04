(ns fastmath.optimization.lbfgsb
  "Box constrained minimization and maximization with the L-BFGS-B algorithm.

  The namespace wraps the Java implementation of L-BFGS-B bundled with the library. The algorithm is a limited-memory quasi-Newton method which works with lower and upper bounds on every variable and needs the gradient of the objective. The gradient can be given by the user or is approximated with finite differences.

  The single optimizer function is [[lbfgsb]]. The algorithm settings (memory size, tolerances, line search) are built by [[parameters]] from the same options map.

  Bounds are always required and are validated before the optimization starts (see [[fastmath.optimization.common/normalize-bounds]]). Infinite bounds are allowed when the initial point is given. The objective function is never evaluated outside of the bounds, also when the gradient is approximated numerically. A non-finite objective value or gradient stops the optimization with an exception."
  (:require [fastmath.core :as m]
            [fastmath.optimization.common :as common])
  (:import [org.generateme.lbfgsb Parameters Parameters$LINESEARCH LBFGSB LBFGSB$Status IGradFunction]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)

(def ^:private line-searches
  {:more-thuente Parameters$LINESEARCH/MORETHUENTE_ORIG
   :orig Parameters$LINESEARCH/MORETHUENTE_ORIG
   :more-thuente-lbfgspp Parameters$LINESEARCH/MORETHUENTE_LBFGSPP
   :lbfgsb Parameters$LINESEARCH/MORETHUENTE_LBFGSPP
   :lewis-overton Parameters$LINESEARCH/LEWISOVERTON})

(defn- line-search-method
  ^Parameters$LINESEARCH [linesearch]
  (or (line-searches linesearch)
      (common/throw-unknown :linesearch linesearch (vec (keys line-searches)))))

(defn parameters
  "Creates the settings of the L-BFGS-B algorithm.

  Parameters:

  - `opts` (map): all keys are optional.
    - `:m` - number of stored correction pairs, the memory of the method, default: `6`.
    - `:abs`, `:rel` - absolute and relative tolerance of the projected gradient norm, default: `1.0e-8` each. The optimization stops when the norm is below `:abs` or below `:rel` times the norm of the point.
    - `:past` - number of past iterations used by the stall test, default: `3`. `0` disables the test.
    - `:delta` - the optimization stops when the objective changed by less than this value (scaled by the objective magnitude, at least 1) over the last `:past` iterations, default: `1.0e-10`.
    - `:max-iters` - maximum number of iterations, default: `1000`. `0` means unlimited.
    - `:max-submin` - maximum number of iterations of the subspace minimization, default: `10`.
    - `:max-linesearch` - maximum number of line search iterations, default: `20`.
    - `:linesearch` - line search method, one of `:more-thuente` (or `:orig`, default), `:more-thuente-lbfgspp` (or `:lbfgsb`) and `:lewis-overton` (experimental).
    - `:xtol` - relative tolerance of the line search interval, default: `1.0e-8`.
    - `:min-step`, `:max-step` - bounds on the line search step, default: `1.0e-20` and `1.0e20`.
    - `:ftol` - sufficient decrease parameter of the line search, default: `1.0e-4`.
    - `:wolfe` - curvature condition parameter of the line search, default: `0.9`.
    - `:weak-wolfe?` - use the weak Wolfe condition, default: `true`.
    - `:debug?` - print the progress of the algorithm to the standard output, default: `false`. The flag is global for the Java implementation and is set on every call.

  Returns a `Parameters` object.

  Throws `ex-info` when `:m`, `:abs`, `:rel`, `:max-linesearch` or `:min-step` is not positive, when `:past`, `:delta`, `:max-iters` or `:max-submin` is negative, when `:max-step` is lower than `:min-step`, when `0 < :ftol < 0.5` and `:ftol < :wolfe < 1` do not hold, and for an unknown `:linesearch`.

  See also [[lbfgsb]]."
  [{:keys [^int m ^double rel ^double abs ^int past ^double delta ^int max-iters ^int max-submin ^int max-linesearch
           linesearch ^double xtol ^double min-step ^double max-step ^double ftol ^double wolfe ^boolean weak-wolfe? ^boolean debug?]
    :or {debug? false m 6 abs 1.0e-8 rel 1.0e-8 past 3 delta 1.0e-10 max-iters 1000 max-submin 10 max-linesearch 20
         linesearch :more-thuente xtol 1.0e-8 min-step 1.0e-20 max-step 1.0e20 ftol 1.0e-4 wolfe 0.9 weak-wolfe? true}}]
  (when-not (m/pos? m) (throw (ex-info "m must be positive" {:m m})))
  (when-not (m/pos? abs) (throw (ex-info "abs must be positive" {:abs abs})))
  (when-not (m/pos? rel) (throw (ex-info "rel must be positive" {:rel rel})))
  (when-not (m/pos? max-linesearch) (throw (ex-info "max-linesearch must be positive" {:max-linesearch max-linesearch})))
  (when-not (m/pos? min-step) (throw (ex-info "min-step must be positive" {:min-step min-step})))
  (when-not (m/not-neg? past) (throw (ex-info "past must be non-negative" {:past past})))
  (when-not (m/not-neg? delta) (throw (ex-info "delta must be non-negative" {:delta delta})))
  (when-not (m/not-neg? max-iters) (throw (ex-info "max-iters must be non-negative" {:max-iters max-iters})))
  (when-not (m/not-neg? max-submin) (throw (ex-info "max-submin must be non-negative" {:max-submin max-submin})))
  (when-not (m/>= max-step min-step) (throw (ex-info "max-step must be greater than min-step" {:min-step min-step :max-step max-step})))
  (when-not (and (m/< 0.0 ftol 0.5)
                 (m/< ftol wolfe 1.0)) (throw (ex-info "ftol and wolfe must satisfy 0<ftol<0.5 and ftol<wolfe<1.0" {:ftol ftol :wolfe wolfe})))
  (let [^Parameters p (Parameters.)]
    (set! org.generateme.lbfgsb.Debug/DEBUG debug?)
    (set! (.-m p) m)
    (set! (.-epsilon p) abs)
    (set! (.-epsilon_rel p) rel)
    (set! (.-past p) past)
    (set! (.-delta p) delta)
    (set! (.-max_iterations p) max-iters)
    (set! (.-max_submin p) max-submin)
    (set! (.-max_linesearch p) max-linesearch)
    (set! (.-linesearch p) (line-search-method linesearch))
    (set! (.-xtol p) xtol)
    (set! (.min_step p) min-step)
    (set! (.max_step p) max-step)
    (set! (.ftol p) ftol)
    (set! (.wolfe p) wolfe)
    (set! (.weak_wolfe p) weak-wolfe?)
    p))

(defn- finite-difference-gradient!
  "Fills `g` with the finite difference gradient of `evaluate` at `xs` without leaving the box given by `l` and `u`.

  The central difference with step `h` is used when both `xs[i] - h` and `xs[i] + h` are inside the box. Otherwise the one-sided difference towards the inside is used; when the box is narrower than `h` in the direction, the step is shortened to the distance to the farther bound. A coordinate with a zero-width box gets gradient `0.0`.

  `xs` is changed temporarily and restored before returning."
  [evaluate xs g l u h]
  (let [^doubles xs xs
        ^doubles g g
        ^doubles l l
        ^doubles u u
        h (double h)
        f0 (delay (double (evaluate xs)))
        shifted (fn [i x shift]
                  (let [i (long i)
                        x (double x)]
                    (aset xs i (m/+ x (double shift)))
                    (let [v (double (evaluate xs))]
                      (aset xs i x)
                      v)))]
    (dotimes [i (alength xs)]
      (let [x (aget xs i)
            dl (m/- x (aget l i))
            du (m/- (aget u i) x)]
        (aset g i (double
                   (cond
                     (and (m/>= dl h) (m/>= du h)) (m// (m/- (double (shifted i x h)) (double (shifted i x (m/- h)))) (m/* 2.0 h))
                     (m/>= du h) (m// (m/- (double (shifted i x h)) (double @f0)) h)
                     (m/>= dl h) (m// (m/- (double @f0) (double (shifted i x (m/- h)))) h)
                     :else (let [s (m/max dl du)]
                             (cond
                               (m/<= s 0.0) 0.0
                               (m/>= du dl) (m// (m/- (double (shifted i x s)) (double @f0)) s)
                               :else (m// (m/- (double @f0) (double (shifted i x (m/- s)))) s))))))))))

(defn- grad-function
  "Creates the Java objective with gradient.

  The objective is minimized: for `sign` equal to `-1.0` the value and the gradient of `f` are negated. `f` receives the point as a sequence. `gradient` (a function of the point returning a sequence of `n` numbers or `nil`) is used when given, otherwise the gradient is approximated with step `h` inside the box `[l, u]`."
  ^IGradFunction [f gradient sign h l u]
  (let [sign (double sign)
        evaluate (fn ^double [^doubles xs] (m/* sign (double (f xs))))]
    (if gradient
      (reify IGradFunction
        (evaluate [_ xs] (evaluate xs))
        (gradient [_ xs g]
          (let [^doubles res (or (m/seq->double-array (gradient xs)) (double-array 0))
                n (alength ^doubles g)]
            (when-not (m/== n (alength res))
              (throw (ex-info "Gradient has wrong length" {:expected n :actual (alength res)})))
            (dotimes [i n]
              (aset ^doubles g i (m/* sign (aget res i)))))))
      (reify IGradFunction
        (evaluate [_ xs] (evaluate xs))
        (gradient [_ xs g] (finite-difference-gradient! evaluate xs g l u h))))))

(defn- status-keyword
  [^LBFGSB$Status status]
  (condp identical? status
    LBFGSB$Status/CONVERGED :converged
    LBFGSB$Status/STALLED :stalled
    LBFGSB$Status/MAX_ITERATIONS :max-iterations))

(defn lbfgsb
  "Minimizes or maximizes a function within box constraints using the L-BFGS-B algorithm.

  Parameters:

  - `f` (function): the objective. It receives the point as one sequence (a `double[]` here), or as separate arguments when `:vector-arg?` is `false`, and returns a number. It is never called outside of `:bounds`.
  - `opts` (map): all keys of [[parameters]] and:
    - `:bounds` (required) - sequence of `[lo hi]` pairs, one for each dimension. Infinite values are allowed when `:initial` is given. A flat `[lo hi]` is accepted for one dimension.
    - `:initial` - initial point, default: the middle of the bounds. A point outside of the bounds is moved to the nearest bound.
    - `:goal` - `:minimize` (default) or `:maximize`.
    - `:vector-arg?` - when `false`, `f` takes separate arguments, default: `true`.
    - `:gradient` - function of the point (always one sequence) returning the gradient of `f` as a sequence of numbers. It is the gradient of `f` itself, also when maximizing. Default: the gradient is approximated with finite differences.
    - `:gradient-h` - step of the finite differences, default: `1.0e-6`. Used only when `:gradient` is not given. Near a bound a one-sided difference is used.
    - `:stats?` - return a map with additional information, default: `false`.

  Returns `[point value]`, where `point` is a vector and `value` is the value of `f`. When `:stats?` is `true`, returns a map with:

  - `:point` and `:value` - as above,
  - `:iterations` - number of iterations,
  - `:gradient` - the gradient of `f` at the point,
  - `:status` - why the optimization stopped: `:converged` (projected gradient norm below `:abs` or `:rel`), `:stalled` (the objective stopped changing, see `:delta` and `:past`) or `:max-iterations`.

  Reaching `:max-iters` is not an error and ends with status `:max-iterations`. `:max-evals` is not supported.

  Throws `ex-info` for invalid bounds, options or gradient length, and an exception when the objective value or the gradient is not finite or when the line search fails.

  See also [[parameters]], [[fastmath.optimization/minimize]], [[fastmath.optimization/maximize]]."
  [f {:keys [bounds initial goal vector-arg? gradient gradient-h stats?]
      :as opts}]
  (let [goal (common/parse-goal goal)
        vector-arg? (common/resolve-vector-arg? :lbfgsb vector-arg?)
        gradient-h (double (or gradient-h 1.0e-6))
        _ (when-not (m/pos? gradient-h) (throw (ex-info "gradient-h must be positive" {:gradient-h gradient-h})))
        bounds (common/normalize-bounds :lbfgsb bounds initial)
        l (double-array (map first bounds))
        u (double-array (map second bounds))
        initial (if initial (m/seq->double-array initial) (common/bounds-midpoint bounds))
        sign (if (= goal :maximize) -1.0 1.0)
        gf (grad-function (common/->vector-fn f vector-arg?) gradient sign gradient-h l u)
        ^LBFGSB optimizer (LBFGSB. (parameters opts))
        x (.minimize optimizer gf initial l u)
        value (m/* sign (.-fx optimizer))]
    (if stats?
      {:point (vec x)
       :value value
       :iterations (.-k optimizer)
       :gradient (mapv #(m/* sign (double %)) (.-m_grad optimizer))
       :status (status-keyword (.-status optimizer))}
      [(vec x) value])))
