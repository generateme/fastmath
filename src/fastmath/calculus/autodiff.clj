(ns fastmath.calculus.autodiff
  "Forward-mode automatic differentiation with tagged dual numbers.

  A dual number carries a value and a vector of partial derivatives (perturbations). Arithmetic and elementary functions defined in this namespace propagate the perturbations by the rules of calculus, so a function built from them computes its derivative alongside its value, exactly up to the floating point rounding. The design follows ForwardDiff.jl.

  Every call of a derivative creates a fresh tag and marks its perturbation with it. Duals with different tags never mix their perturbations, which makes nested derivatives (a derivative of a function which itself takes a derivative) correct and avoids perturbation confusion. A dual with a newer tag is the outer layer, its value and partials can be duals with older tags.

  A function to differentiate has to use the operations of this namespace instead of `clojure.core` or `fastmath.core` arithmetic. All operations accept plain numbers and duals, a number is treated as a constant and coerced to a double. A `Dual` is not a `java.lang.Number`, so it can not be passed silently to a primitive function, which would drop the perturbation.

  Defined functions:

  - derivatives: [[derivative]],
  - arithmetic: [[+]], [[-]], [[*]], [[/]],
  - elementary functions: [[sin]], [[cos]], [[exp]], [[log]], [[pow]],
  - comparisons: [[<]], [[<=]], [[>]], [[>=]], [[=]],
  - duals: [[dual]], [[primal]], [[scalar]].

  Limitations: a function returning a function can not be differentiated (it throws). A dual which escapes its derivative call (stored in an atom, a cache or a lazy sequence realized later) is detected only when it flows into the result of an enclosing derivative, which then throws. An escaped dual used only in a comparison (a branch) is not detected and can silently give a wrong derivative, as can a dual left over from a finished, unrelated derivative call."
  (:refer-clojure :exclude [+ - * / < > <= >= =])
  (:require [fastmath.core :as m]
            [fastmath.calculus.common :as cc])
  (:import [clojure.lang IPersistentVector]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; ## Tags

(defonce ^:private tag-counter (atom 0))

(defn- fresh-tag
  "Returns a new tag, larger than all tags returned before."
  ^long []
  (long (swap! tag-counter inc)))

;; ## Protocol and the dual type

(defprotocol ADProto
  "Operations on numbers and dual numbers, dispatched on the first argument. A number argument is coerced to a double."
  (-add [x y] "Sum of `x` and `y`.")
  (-sub [x y] "Difference of `x` and `y`.")
  (-mul [x y] "Product of `x` and `y`.")
  (-div [x y] "Quotient of `x` and `y`.")
  (-pow [x y] "`x` raised to the power `y`.")
  (-cmp [x y] "Lexicographic comparison of (value, partials): `-1`, `0` or `1`, and `2` when unordered (NaN).")
  (-neg [x] "Negation of `x`.")
  (-sin [x] "Sine of `x`.")
  (-cos [x] "Cosine of `x`.")
  (-exp [x] "Exponential of `x`.")
  (-log [x] "Natural logarithm of `x`."))

(declare dual-add dual-sub dual-mul dual-div dual-pow dual-cmp
         dual-neg dual-sin dual-cos dual-exp dual-log)

;; Invariant: every tag found in `value` and in `partials` is smaller than `tag`.
(deftype Dual [^long tag value ^IPersistentVector partials]
  Object
  (toString [_] (str "#dual[" tag " " (pr-str value) " " (pr-str partials) "]"))
  ADProto
  (-add [x y] (dual-add x y))
  (-sub [x y] (dual-sub x y))
  (-mul [x y] (dual-mul x y))
  (-div [x y] (dual-div x y))
  (-pow [x y] (dual-pow x y))
  (-cmp [x y] (dual-cmp x y))
  (-neg [x] (dual-neg x))
  (-sin [x] (dual-sin x))
  (-cos [x] (dual-cos x))
  (-exp [x] (dual-exp x))
  (-log [x] (dual-log x)))

(defmethod print-method Dual [^Dual d ^java.io.Writer w]
  (.write w (.toString d)))

(defn- leaf
  "A dual unchanged, a number coerced to a double."
  [x]
  (if (instance? Dual x) x (.doubleValue ^Number x)))

(defn- cmp-doubles
  ^long [^double a ^double b]
  (cond (m/< a b) -1
        (m/> a b) 1
        (m/== a b) 0
        :else 2))

(extend-type Number
  ADProto
  (-add [x y] (if (instance? Dual y) (dual-add x y) (m/+ (.doubleValue x) (.doubleValue ^Number y))))
  (-sub [x y] (if (instance? Dual y) (dual-sub x y) (m/- (.doubleValue x) (.doubleValue ^Number y))))
  (-mul [x y] (if (instance? Dual y) (dual-mul x y) (m/* (.doubleValue x) (.doubleValue ^Number y))))
  (-div [x y] (if (instance? Dual y) (dual-div x y) (m// (.doubleValue x) (.doubleValue ^Number y))))
  (-pow [x y] (if (instance? Dual y) (dual-pow x y) (m/pow (.doubleValue x) (.doubleValue ^Number y))))
  (-cmp [x y] (if (instance? Dual y) (dual-cmp x y) (cmp-doubles (.doubleValue x) (.doubleValue ^Number y))))
  (-neg [x] (m/- (.doubleValue x)))
  (-sin [x] (m/sin (.doubleValue x)))
  (-cos [x] (m/cos (.doubleValue x)))
  (-exp [x] (m/exp (.doubleValue x)))
  (-log [x] (m/log (.doubleValue x))))

;; ## Tag layers

(defn- tag-of
  "Tag of a dual, `-1` for a number."
  ^long [x]
  (if (instance? Dual x) (.tag ^Dual x) -1))

(defn- top-tag
  ^long [x y]
  (Math/max (tag-of x) (tag-of y)))

(defn- value-at
  "Value of `x` in the layer of tag `t`: the value of a dual tagged `t`, otherwise `x` itself (a constant)."
  [x ^long t]
  (if (and (instance? Dual x) (m/== t (.tag ^Dual x))) (.value ^Dual x) x))

(defn- partials-at
  "Partials of `x` in the layer of tag `t`, `nil` (zero) for a constant."
  [x ^long t]
  (when (and (instance? Dual x) (m/== t (.tag ^Dual x))) (.partials ^Dual x)))

;; ## Partials (`nil` is a zero vector)

(declare zero-element?)

(defn- p-zero?
  "Are all partials zero (`nil` is zero)? Stops at the first nonzero element."
  [p]
  (or (nil? p)
      (let [^IPersistentVector v p
            n (.count v)]
        (loop [i 0]
          (cond
            (m/== i n) true
            (zero-element? (.nth v i)) (recur (inc i))
            :else false)))))

(defn- zero-element?
  "Is `x` zero, including all of its partials? The value is checked first, so a nonzero scalar value is found by following the values only."
  [x]
  (if (instance? Dual x)
    (and (zero-element? (.value ^Dual x)) (p-zero? (.partials ^Dual x)))
    (m/zero? (.doubleValue ^Number x))))

(defn- check-same-length!
  "Duals with the same tag must have the same number of partials (a precondition of the caller)."
  [^IPersistentVector px ^IPersistentVector py]
  (when-not (m/== (.count px) (.count py))
    (throw (ex-info "Duals with the same tag have different numbers of partials."
                    {:count-x (.count px) :count-y (.count py)}))))

(defn- p-zeros [p] (vec (repeat (count p) 0.0)))

(defn- p-add [px py]
  (cond (nil? px) py
        (nil? py) px
        :else (do (check-same-length! px py) (mapv -add px py))))

(defn- p-sub [px py]
  (cond (nil? py) px
        (nil? px) (mapv -neg py)
        :else (do (check-same-length! px py) (mapv -sub px py))))

(defn- p-scale [p s] (when p (mapv (fn [e] (-mul e s)) p)))

(defn- p-div [p s] (when p (mapv (fn [e] (-div e s)) p)))

(defn- p-cmp
  "Lexicographic comparison of partials, `nil` is a zero vector."
  ^long [px py]
  (when (and px py) (check-same-length! px py))
  (let [n (count (or px py))]
    (loop [i 0]
      (if (m/< i n)
        (let [c (long (-cmp (if px (px i) 0.0) (if py (py i) 0.0)))]
          (if (m/zero? c) (recur (inc i)) c))
        0))))

;; ## Operations on duals (at least one argument is a dual)

(defn- dual-add [x y]
  (let [t (top-tag x y)]
    (Dual. t (-add (value-at x t) (value-at y t)) (p-add (partials-at x t) (partials-at y t)))))

(defn- dual-sub [x y]
  (let [t (top-tag x y)]
    (Dual. t (-sub (value-at x t) (value-at y t)) (p-sub (partials-at x t) (partials-at y t)))))

(defn- dual-mul [x y]
  (let [t (top-tag x y)
        vx (value-at x t)
        vy (value-at y t)]
    (Dual. t (-mul vx vy) (p-add (p-scale (partials-at x t) vy) (p-scale (partials-at y t) vx)))))

(defn- dual-div [x y]
  (let [t (top-tag x y)
        vx (value-at x t)
        vy (value-at y t)
        q (-div vx vy)]
    ;; d(x/y) = dx/y - dy * (x/y)/y
    (Dual. t q (p-sub (p-div (partials-at x t) vy) (p-scale (partials-at y t) (-div q vy))))))

(defn- positive-element? [x] (m/== 1 (long (-cmp x 0.0))))

(defn- pow-dual-const
  "`x` is tagged `t`, `y` is a constant in this layer."
  [x y ^long t]
  (let [vx (value-at x t)
        px (partials-at x t)]
    (Dual. t (-pow vx y)
           (if (or (zero-element? y) (p-zero? px))
             (p-zeros px)
             (p-scale px (-mul y (-pow vx (-sub y 1.0))))))))

(defn- pow-const-dual
  "`x` is a constant in this layer, `y` is tagged `t`."
  [x y ^long t]
  (let [vy (value-at y t)
        py (partials-at y t)
        e (-pow x vy)]
    (Dual. t e
           (if (or (p-zero? py) (and (zero-element? x) (positive-element? vy)))
             (p-zeros py)
             (p-scale py (-mul e (-log x)))))))

(defn- pow-dual-dual
  "Both `x` and `y` are tagged `t`."
  [x y ^long t]
  (let [px (partials-at x t)
        py (partials-at y t)
        vx (value-at x t)
        vy (value-at y t)]
    (cond
      (p-zero? py) (pow-dual-const x vy t)
      (p-zero? px) (pow-const-dual vx y t)
      :else (let [e (-pow vx vy)
                  dx (-mul vy (-pow vx (-sub vy 1.0)))
                  dy (if (and (zero-element? vx) (positive-element? vy)) 0.0 (-mul e (-log vx)))]
              (Dual. t e (p-add (p-scale px dx) (p-scale py dy)))))))

(defn- dual-pow [x y]
  (let [t (top-tag x y)]
    (cond
      (nil? (partials-at y t)) (pow-dual-const x y t)
      (nil? (partials-at x t)) (pow-const-dual x y t)
      :else (pow-dual-dual x y t))))

(defn- dual-cmp [x y]
  (let [t (top-tag x y)
        c (long (-cmp (value-at x t) (value-at y t)))]
    (if (m/zero? c) (p-cmp (partials-at x t) (partials-at y t)) c)))

(defn- dual-neg [^Dual x]
  (Dual. (.tag x) (-neg (.value x)) (mapv -neg (.partials x))))

(defn- dual-sin [^Dual x]
  (let [v (.value x)]
    (Dual. (.tag x) (-sin v) (p-scale (.partials x) (-cos v)))))

(defn- dual-cos [^Dual x]
  (let [v (.value x)]
    (Dual. (.tag x) (-cos v) (p-scale (.partials x) (-neg (-sin v))))))

(defn- dual-exp [^Dual x]
  (let [e (-exp (.value x))]
    (Dual. (.tag x) e (p-scale (.partials x) e))))

(defn- dual-log [^Dual x]
  (let [v (.value x)]
    (Dual. (.tag x) (-log v) (p-div (.partials x) v))))

;; ## Public API: duals

(defn dual
  "Creates a dual number with the given tag, value and partial derivatives.

  Duals are created by [[derivative]], this function is meant for tests and experiments.

  Warning: the tags of this function and the tags created by [[derivative]] are the same numbers. A dual created here with a tag of a running derivative call is mixed with its perturbation and gives a wrong derivative (perturbation confusion). Do not use duals created here inside a function being differentiated.

  Parameters:

  - `tag` (integer >= 0): the tag of the perturbation. Tags created by [[derivative]] are positive and grow with each call.
  - `value` (number or dual): the value. A number is coerced to a double.
  - `partials` (non-empty sequence of numbers or duals): the partial derivatives, numbers are coerced to doubles.

  All duals found in `value` and `partials` must have tags smaller than `tag`. Duals with the same tag, combined in one computation, must have the same number of partials, otherwise the operation throws `ex-info`.

  Returns a `Dual`.

  Throws `ex-info` with `{:option :value :expected}` when `tag` is negative, `partials` is empty or a nested dual has a tag not smaller than `tag`.

  See also [[primal]], [[scalar]], [[derivative]]."
  [tag value partials]
  (let [t (long tag)
        v (leaf value)
        p (mapv leaf partials)]
    (when (m/neg? t) (cc/invalid-option! :tag tag "an integer >= 0"))
    (when (empty? p) (cc/invalid-option! :partials partials "a non-empty sequence"))
    (when-not (and (m/< (tag-of v) t) (every? (fn [e] (m/< (tag-of e) t)) p))
      (cc/invalid-option! :tag tag "a tag larger than the tags of the nested duals"))
    (Dual. t v p)))

(defn primal
  "Value of a dual without its outermost perturbation.

  Parameters:

  - `x` (number or dual).

  Returns the value of `x`, which is a dual with an older tag for a nested dual. A number is returned as a double.

  See also [[scalar]], [[dual]]."
  [x]
  (if (instance? Dual x) (.value ^Dual x) (leaf x)))

(defn scalar
  "Value of a number or a dual with all perturbations removed.

  Parameters:

  - `x` (number or dual).

  Returns a double. Use it to branch on the value of a dual when the perturbations should be ignored, see [[<]] for the comparisons which take them into account.

  See also [[primal]], [[dual]]."
  ^double [x]
  (if (instance? Dual x) (recur (.value ^Dual x)) (.doubleValue ^Number x)))

;; ## Public API: arithmetic and functions

(defn +
  "Sum of numbers and duals.

  Parameters:

  - zero or more numbers or duals.

  Returns `0.0` for no arguments, the argument for one argument (a number coerced to a double), otherwise the sum. The result is a double when all arguments are numbers, otherwise a dual.

  See also [[-]], [[*]], [[/]]."
  ([] 0.0)
  ([x] (leaf x))
  ([x y] (-add x y))
  ([x y & more] (reduce -add (-add x y) more)))

(defn -
  "Difference of numbers and duals, or negation of a single argument.

  Parameters:

  - one or more numbers or duals.

  Returns the negation of one argument, otherwise the first argument minus the rest. The result is a double when all arguments are numbers, otherwise a dual.

  See also [[+]], [[*]], [[/]]."
  ([x] (-neg x))
  ([x y] (-sub x y))
  ([x y & more] (reduce -sub (-sub x y) more)))

(defn *
  "Product of numbers and duals.

  Parameters:

  - zero or more numbers or duals.

  Returns `1.0` for no arguments, the argument for one argument (a number coerced to a double), otherwise the product. The result is a double when all arguments are numbers, otherwise a dual.

  See also [[+]], [[-]], [[/]]."
  ([] 1.0)
  ([x] (leaf x))
  ([x y] (-mul x y))
  ([x y & more] (reduce -mul (-mul x y) more)))

(defn /
  "Quotient of numbers and duals, or reciprocal of a single argument.

  Parameters:

  - one or more numbers or duals.

  Returns the reciprocal of one argument, otherwise the first argument divided by the rest. Division by zero follows the floating point rules (`Inf` or `NaN`). The result is a double when all arguments are numbers, otherwise a dual.

  See also [[+]], [[-]], [[*]]."
  ([x] (-div 1.0 x))
  ([x y] (-div x y))
  ([x y & more] (reduce -div (-div x y) more)))

(defn sin
  "Sine of a number or a dual.

  Parameters:

  - `x` (number or dual): angle in radians.

  Returns a double for a number, otherwise a dual.

  See also [[cos]]."
  [x] (-sin x))

(defn cos
  "Cosine of a number or a dual.

  Parameters:

  - `x` (number or dual): angle in radians.

  Returns a double for a number, otherwise a dual.

  See also [[sin]]."
  [x] (-cos x))

(defn exp
  "Exponential function of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual.

  See also [[log]], [[pow]]."
  [x] (-exp x))

(defn log
  "Natural logarithm of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The value is `NaN` for a negative argument and `-Inf` for zero; the derivative is `1/x`, infinite at zero.

  See also [[exp]], [[pow]]."
  [x] (-log x))

(defn pow
  "Power `x^y` of numbers and duals.

  Parameters:

  - `x` (number or dual): base.
  - `y` (number or dual): exponent.

  Returns a double when both arguments are numbers, otherwise a dual. The value follows `Math/pow`, so it is `NaN` for a negative base and a fractional exponent.

  Derivative conventions:

  - a constant exponent `y` gives the derivative `y x^(y-1)`, zero for `y = 0` (also at `x = 0`), finite for a negative base and an integer exponent; it is infinite for `0 < y < 1` at `x = 0`,
  - a constant base `x` gives the derivative `x^y log(x)`, zero for `x = 0` and `y > 0`, `NaN` for a negative base,
  - an exponent or a base which is a dual with zero partials is treated as a constant and the rules above apply,
  - when both depend on the variable, both terms are added; the result is `NaN` where the logarithm of the base is not defined.

  See also [[exp]], [[log]]."
  [x y] (-pow x y))

;; ## Public API: comparisons

(def ^:private comparison-doc
  "

  Numbers and duals are compared lexicographically: first the values, then, for equal values, the partial derivatives of the newest tag, a constant having zero partials. A dual is therefore not equal to its value when its partials are not zero, and a comparison at a tie follows the direction of the perturbation (a branch on `x < 0` at `x = 0` gives the one-sided derivative in the direction of the perturbation). Any comparison with `NaN` is false. To compare values only, use [[scalar]].

  Parameters:

  - one or more numbers or duals.

  Returns `true` for one argument, otherwise whether the relation holds for each consecutive pair.")

(defn- comparison-chain
  "Does `relation` hold for `x`, `y` and each consecutive pair of `more`?"
  [relation x y more]
  (if (relation x y)
    (if (next more)
      (recur relation y (first more) (next more))
      (relation y (first more)))
    false))

(defn- comparison-docstring [head see]
  (str head comparison-doc "\n\n  See also " see "."))

(defn ^{:doc (comparison-docstring "Are the arguments strictly increasing?" "[[<=]], [[>]], [[>=]], [[=]]")} <
  ([_] true)
  ([x y] (m/== -1 (long (-cmp x y))))
  ([x y & more] (comparison-chain < x y more)))

(defn ^{:doc (comparison-docstring "Are the arguments non-decreasing?" "[[<]], [[>]], [[>=]], [[=]]")} <=
  ([_] true)
  ([x y] (let [c (long (-cmp x y))] (or (m/== c -1) (m/zero? c))))
  ([x y & more] (comparison-chain <= x y more)))

(defn ^{:doc (comparison-docstring "Are the arguments strictly decreasing?" "[[<]], [[<=]], [[>=]], [[=]]")} >
  ([_] true)
  ([x y] (m/== 1 (long (-cmp x y))))
  ([x y & more] (comparison-chain > x y more)))

(defn ^{:doc (comparison-docstring "Are the arguments non-increasing?" "[[<]], [[<=]], [[>]], [[=]]")} >=
  ([_] true)
  ([x y] (let [c (long (-cmp x y))] (or (m/== c 1) (m/zero? c))))
  ([x y & more] (comparison-chain >= x y more)))

(defn ^{:doc (comparison-docstring "Are the arguments equal, including their partial derivatives?" "[[<]], [[<=]], [[>]], [[>=]], [[scalar]]")} =
  ([_] true)
  ([x y] (m/zero? (long (-cmp x y))))
  ([x y & more] (comparison-chain = x y more)))

;; ## Derivative

(defn- extract-derivative
  "The first partial of the layer `t` of the `result` of the function."
  [^long t result]
  (cond
    (instance? Dual result) (let [rt (.tag ^Dual result)]
                              (cond
                                (m/== rt t) ((.partials ^Dual result) 0)
                                (m/< rt t) 0.0
                                :else (throw (ex-info "A dual escaped from its derivative call into the result."
                                                      {:tag t :result-tag rt}))))
    (number? result) 0.0
    (fn? result) (throw (ex-info "Derivative of a function returning a function is not supported."
                                 {:tag t :result-type (type result)}))
    :else (throw (ex-info "Derivative of a function returning neither a number nor a dual is not supported."
                          {:tag t :result-type (type result)}))))

(defn- first-derivative
  [f]
  (fn [x]
    (let [t (fresh-tag)]
      (extract-derivative t (f (Dual. t (leaf x) [1.0]))))))

(defn derivative
  "Creates the `n`-th derivative of a univariate function `f` with forward-mode automatic differentiation.

  The function `f` has to be built from the operations of this namespace ([[+]], [[*]], [[sin]], [[pow]], ...), which accept numbers and duals; the derivative is exact up to the floating point rounding. Derivatives can be nested: `f` can call a derivative itself, every call of the returned function uses a fresh tag, so the perturbations of different derivative calls are not confused.

  Parameters:

  - `f` (function): number or dual -> number or dual.
  - `n` (integer >= 0): the order of the derivative (default: `1`). For `0` the function `f` itself is returned. A fractional `n` is truncated.

  Returns a function of `x` (a number or a dual) returning the derivative: a double, or a dual when `x` is a dual of an enclosing derivative. The derivative of a function which does not depend on `x` is `0.0`.

  Throws `ex-info` with `{:option :value :expected}` when `n` is negative. The returned function throws `ex-info` when `f` returns a function (derivatives of functions returning functions are not supported), when it returns something else than a number or a dual (for example a vector), and when a dual of a finished derivative call (stored in an atom, a cache or a lazy sequence) leaks into the result. An escaped dual which only takes part in a comparison is not detected.

  The order `n` is computed by `n` nested derivatives, every level multiplies the number of terms carried by the duals, so the cost grows geometrically with `n` (about three times per order for a typical expression).

  See also [[fastmath.calculus/derivative]] (finite differences), [[dual]]."
  ([f] (derivative f 1))
  ([f ^long n]
   (when (m/neg? n) (cc/invalid-option! :n n "an integer >= 0"))
   (loop [g f
          k n]
     (if (m/zero? k)
       g
       (recur (first-derivative g) (dec k))))))
