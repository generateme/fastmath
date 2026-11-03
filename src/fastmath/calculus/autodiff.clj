(ns fastmath.calculus.autodiff
  "Forward-mode automatic differentiation with tagged dual numbers.

  A dual number carries a value and its partial derivatives (perturbations), one for each input differentiated at once. Arithmetic and elementary functions defined in this namespace propagate the perturbations by the rules of calculus, so a function built from them computes its derivative alongside its value, exactly up to the floating point rounding. The design follows ForwardDiff.jl.

  Every call of a derivative creates a fresh tag and marks its perturbation with it. Duals with different tags never mix their perturbations, which makes nested derivatives (a derivative of a function which itself takes a derivative) correct and avoids perturbation confusion. A dual with a newer tag is the outer layer, its value and partials can be duals with older tags.

  A function to differentiate has to use the operations of this namespace instead of `clojure.core` or `fastmath.core` arithmetic. Numeric equality is [[==]]; `clojure.core/=` compares duals structurally (tag, value and partials), consistently with their hash, so duals can be used in sets and as map keys. All operations accept plain numbers and duals, a number is treated as a constant and coerced to a double. A `Dual` is not a `java.lang.Number`, so it can not be passed silently to a primitive function, which would drop the perturbation.

  Defined functions:

  - derivatives: [[derivative]], [[gradient]], [[jacobian]], [[hessian]],
  - arithmetic: [[+]], [[-]], [[*]], [[/]],
  - elementary functions: [[sin]], [[cos]], [[exp]], [[log]], [[pow]],
  - comparisons: [[<]], [[<=]], [[>]], [[>=]], [[==]],
  - duals: [[dual]], [[primal]], [[scalar]].

  Limitations: a function returning a function can not be differentiated (it throws). A dual which escapes its derivative call (stored in an atom, a cache or a lazy sequence realized later) is detected only when it flows into the result of an enclosing derivative, which then throws. An escaped dual used only in a comparison (a branch) is not detected and can silently give a wrong derivative, as can a dual left over from a finished, unrelated derivative call, or a dual of an earlier pass of the same [[gradient]], [[jacobian]] or [[hessian]] call reused in a later pass: such duals have older tags and are treated as constants."
  (:refer-clojure :exclude [+ - * / < > <= >= ==])
  (:require [fastmath.core :as m]
            [fastmath.calculus.common :as cc]
            [fastmath.calculus.autodiff.dual :as d]
            [fastmath.calculus.autodiff.ops :as ops]
            [fastmath.calculus.autodiff.derivatives :as der])
  (:import [fastmath.calculus.autodiff.dual Dual]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; The public API. The implementation lives in the internal namespaces `fastmath.calculus.autodiff.dual` (the dual
;; type, tags, partials containers), `fastmath.calculus.autodiff.ops` (arithmetic and derivative rules) and
;; `fastmath.calculus.autodiff.derivatives` (seeding, extraction and passes).

;; ## Duals

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
        v (d/leaf value)
        p (mapv d/leaf partials)]
    (when (m/neg? t) (cc/invalid-option! :tag tag "an integer >= 0"))
    (when (empty? p) (cc/invalid-option! :partials partials "a non-empty sequence"))
    (when-not (and (m/< (d/tag-of v) t) (every? (fn [e] (m/< (d/tag-of e) t)) p))
      (cc/invalid-option! :tag tag "a tag larger than the tags of the nested duals"))
    (d/new-dual t v (object-array p))))

(defn primal
  "Value of a dual without its outermost perturbation.

  Parameters:

  - `x` (number or dual).

  Returns the value of `x`, which is a dual with an older tag for a nested dual. A number is returned as a double.

  See also [[scalar]], [[dual]]."
  [x]
  (if (instance? Dual x) (.value ^Dual x) (d/leaf x)))

(defn scalar
  "Value of a number or a dual with all perturbations removed.

  Parameters:

  - `x` (number or dual).

  Returns a double. Use it to branch on the value of a dual when the perturbations should be ignored, see [[<]] for the comparisons which take them into account.

  See also [[primal]], [[dual]]."
  ^double [x]
  (d/scalar x))

;; ## Public API: arithmetic and functions

(defn +
  "Sum of numbers and duals.

  Parameters:

  - zero or more numbers or duals.

  Returns `0.0` for no arguments, the argument for one argument (a number coerced to a double), otherwise the sum. The result is a double when all arguments are numbers, otherwise a dual.

  See also [[-]], [[*]], [[/]]."
  ([] 0.0)
  ([x] (d/leaf x))
  ([x y] (ops/add x y))
  ([x y & more] (reduce ops/add (ops/add x y) more)))

(defn -
  "Difference of numbers and duals, or negation of a single argument.

  Parameters:

  - one or more numbers or duals.

  Returns the negation of one argument, otherwise the first argument minus the rest. The result is a double when all arguments are numbers, otherwise a dual.

  See also [[+]], [[*]], [[/]]."
  ([x] (ops/neg x))
  ([x y] (ops/sub x y))
  ([x y & more] (reduce ops/sub (ops/sub x y) more)))

(defn *
  "Product of numbers and duals.

  Parameters:

  - zero or more numbers or duals.

  Returns `1.0` for no arguments, the argument for one argument (a number coerced to a double), otherwise the product. The result is a double when all arguments are numbers, otherwise a dual.

  See also [[+]], [[-]], [[/]]."
  ([] 1.0)
  ([x] (d/leaf x))
  ([x y] (ops/mul x y))
  ([x y & more] (reduce ops/mul (ops/mul x y) more)))

(defn /
  "Quotient of numbers and duals, or reciprocal of a single argument.

  Parameters:

  - one or more numbers or duals.

  Returns the reciprocal of one argument, otherwise the first argument divided by the rest. Division by zero follows the floating point rules (`Inf` or `NaN`). The result is a double when all arguments are numbers, otherwise a dual.

  See also [[+]], [[-]], [[*]]."
  ([x] (ops/div 1.0 x))
  ([x y] (ops/div x y))
  ([x y & more] (reduce ops/div (ops/div x y) more)))

(defn sin
  "Sine of a number or a dual.

  Parameters:

  - `x` (number or dual): angle in radians.

  Returns a double for a number, otherwise a dual.

  See also [[cos]]."
  [x] (ops/sin x))

(defn cos
  "Cosine of a number or a dual.

  Parameters:

  - `x` (number or dual): angle in radians.

  Returns a double for a number, otherwise a dual.

  See also [[sin]]."
  [x] (ops/cos x))

(defn exp
  "Exponential function of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual.

  See also [[log]], [[pow]]."
  [x] (ops/exp x))

(defn log
  "Natural logarithm of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The value is `NaN` for a negative argument and `-Inf` for zero; the derivative is `1/x`, infinite at zero and `NaN` (as all higher derivatives) for a negative argument, outside the domain.

  See also [[exp]], [[pow]]."
  [x] (ops/log x))

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
  [x y] (ops/pow x y))

;; ## Public API: comparisons


(defn- comparison-chain
  "Does `relation` hold for `x`, `y` and each consecutive pair of `more`?"
  [relation x y more]
  (if (relation x y)
    (if (next more)
      (recur relation y (first more) (next more))
      (relation y (first more)))
    false))

(defn <
  "Are the arguments strictly increasing?

  Numbers and duals are compared lexicographically: first the values, then, for equal values, the partial derivatives of the newest tag, a constant having zero partials. A dual is therefore not equal to its value when its partials are not zero, and a comparison at a tie follows the direction of the perturbation (a branch on `x < 0` at `x = 0` gives the one-sided derivative in the direction of the perturbation). Any comparison with `NaN` is false. To compare values only, use [[scalar]].

  Parameters:

  - one or more numbers or duals.

  Returns `true` for one argument, otherwise whether the relation holds for each consecutive pair.

  See also [[<=]], [[>]], [[>=]], [[==]]."
  ([_] true)
  ([x y] (m/== -1 (ops/cmp x y)))
  ([x y & more] (comparison-chain < x y more)))

(defn <=
  "Are the arguments non-decreasing?

  Numbers and duals are compared lexicographically: first the values, then, for equal values, the partial derivatives of the newest tag, a constant having zero partials. A dual is therefore not equal to its value when its partials are not zero, and a comparison at a tie follows the direction of the perturbation (a branch on `x < 0` at `x = 0` gives the one-sided derivative in the direction of the perturbation). Any comparison with `NaN` is false. To compare values only, use [[scalar]].

  Parameters:

  - one or more numbers or duals.

  Returns `true` for one argument, otherwise whether the relation holds for each consecutive pair.

  See also [[<]], [[>]], [[>=]], [[==]]."
  ([_] true)
  ([x y] (let [c (ops/cmp x y)] (or (m/== c -1) (m/zero? c))))
  ([x y & more] (comparison-chain <= x y more)))

(defn >
  "Are the arguments strictly decreasing?

  Numbers and duals are compared lexicographically: first the values, then, for equal values, the partial derivatives of the newest tag, a constant having zero partials. A dual is therefore not equal to its value when its partials are not zero, and a comparison at a tie follows the direction of the perturbation (a branch on `x < 0` at `x = 0` gives the one-sided derivative in the direction of the perturbation). Any comparison with `NaN` is false. To compare values only, use [[scalar]].

  Parameters:

  - one or more numbers or duals.

  Returns `true` for one argument, otherwise whether the relation holds for each consecutive pair.

  See also [[<]], [[<=]], [[>=]], [[==]]."
  ([_] true)
  ([x y] (m/== 1 (ops/cmp x y)))
  ([x y & more] (comparison-chain > x y more)))

(defn >=
  "Are the arguments non-increasing?

  Numbers and duals are compared lexicographically: first the values, then, for equal values, the partial derivatives of the newest tag, a constant having zero partials. A dual is therefore not equal to its value when its partials are not zero, and a comparison at a tie follows the direction of the perturbation (a branch on `x < 0` at `x = 0` gives the one-sided derivative in the direction of the perturbation). Any comparison with `NaN` is false. To compare values only, use [[scalar]].

  Parameters:

  - one or more numbers or duals.

  Returns `true` for one argument, otherwise whether the relation holds for each consecutive pair.

  See also [[<]], [[<=]], [[>]], [[==]]."
  ([_] true)
  ([x y] (let [c (ops/cmp x y)] (or (m/== c 1) (m/zero? c))))
  ([x y & more] (comparison-chain >= x y more)))

(defn ==
  "Are the arguments numerically equal, including their partial derivatives?

  This is the numeric comparison of the tower of numbers and duals. `clojure.core/=` is the structural equality of objects instead: two duals are `=` when they have the same tag, value and partials, a dual is never `=` to a number, and the numbers inside are compared as by `clojure.core/=` (`0.0` equals `-0.0`, `NaN` is not equal to another `NaN` object).

  Numbers and duals are compared lexicographically: first the values, then, for equal values, the partial derivatives of the newest tag, a constant having zero partials. A dual is therefore not equal to its value when its partials are not zero, and a comparison at a tie follows the direction of the perturbation (a branch on `x < 0` at `x = 0` gives the one-sided derivative in the direction of the perturbation). Any comparison with `NaN` is false. To compare values only, use [[scalar]].

  Parameters:

  - one or more numbers or duals.

  Returns `true` for one argument, otherwise whether the relation holds for each consecutive pair.

  See also [[<]], [[<=]], [[>]], [[>=]], [[scalar]]."
  ([_] true)
  ([x y] (m/zero? (ops/cmp x y)))
  ([x y & more] (comparison-chain == x y more)))

;; ## Derivatives

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
  ([f] (der/derivative f 1))
  ([f n] (der/derivative f n)))

(defn gradient
  "Creates the gradient (vector of the first partial derivatives) of a multivariate scalar function `f` with forward-mode automatic differentiation.

  The function `f` receives a vector and has to be built from the operations of this namespace. The gradient is exact up to the floating point rounding. Gradients can be nested with other derivatives: every pass uses a fresh tag.

  The inputs are differentiated in passes over chunks of at most `:chunk-size` inputs, each pass calls `f` once with the inputs of the chunk replaced by duals. The default chunk size splits `n` inputs evenly into chunks of at most 8. Small chunks are faster for sparse functions (each operation depends on few inputs), larger for dense ones. The result does not depend on the chunk size, except at a tie of a comparison (see [[<]]): a tie is decided by the partial derivatives, so a branch taken at a tie can depend on the chunk.

  Parameters:

  - `f` (function): vector of numbers or duals -> number or dual.
  - `options` (map, optional):
    - `:chunk-size` - the number of inputs differentiated in one pass, an integer >= 1, clamped to the number of inputs (default: balanced chunks of at most 8).

  Returns a function of a point `v` (any sequence of numbers or duals, or a `double[]`) returning the gradient as a vector: of doubles, or of duals when `v` contains duals of an enclosing derivative. The gradient of the empty point is `[]` and `f` is not called.

  Throws `ex-info` with `{:option :value :expected}` when `:chunk-size` is less than 1. The returned function throws `ex-info` when `f` returns something else than a number or a dual, and when a dual of a finished derivative call leaks into the result of an enclosing derivative. A dual of an earlier pass of the same call, reused in a later pass, is not detected and is treated as a constant.

  See also [[jacobian]], [[hessian]], [[derivative]], [[fastmath.calculus/gradient]] (finite differences)."
  ([f] (der/gradient f nil))
  ([f options] (der/gradient f options)))

(defn jacobian
  "Creates the Jacobian matrix of a vector valued function `f` with forward-mode automatic differentiation.

  The function `f` receives a vector and has to be built from the operations of this namespace. The Jacobian is exact up to the floating point rounding.

  The inputs are differentiated in passes over chunks of at most `:chunk-size` inputs, each pass calls `f` once with the inputs of the chunk replaced by duals. The default chunk size splits `n` inputs evenly into chunks of at most 8. Small chunks are faster for sparse functions (each operation depends on few inputs), larger for dense ones. The result does not depend on the chunk size, except at a tie of a comparison (see [[<]]): a tie is decided by the partial derivatives, so a branch taken at a tie can depend on the chunk.

  Parameters:

  - `f` (function): vector of `n` numbers or duals -> sequential of `m` numbers or duals.
  - `options` (map, optional): as in [[gradient]]
    - `:chunk-size` - the number of inputs differentiated in one pass, an integer >= 1 (default: balanced chunks of at most 8).

  Returns a function of a point `v` (any sequence of numbers or duals, or a `double[]`) returning the matrix as a vector of `m` rows, each a vector of `n` partial derivatives `d f_j / d x_i`: doubles, or duals when `v` contains duals of an enclosing derivative. For the empty point `f` is called once and the result is a vector of `m` empty rows.

  Throws `ex-info` with `{:option :value :expected}` when `:chunk-size` is less than 1. The returned function throws `ex-info` when `f` returns something else than a sequential (a number, a `double[]`, a function), when an element of the result is not a number or a dual, when the number of values changes between the passes, and when a dual of a finished derivative call leaks into the result of an enclosing derivative. A dual of an earlier pass of the same call, reused in a later pass, is not detected and is treated as a constant.

  See also [[gradient]], [[hessian]], [[fastmath.calculus/jacobian]] (finite differences)."
  ([f] (der/jacobian f nil))
  ([f options] (der/jacobian f options)))

(defn hessian
  "Creates the Hessian matrix (the second partial derivatives) of a multivariate scalar function `f` with forward-mode automatic differentiation.

  The Hessian is computed as the Jacobian of the gradient (nested forward mode), it is exact up to the floating point rounding. The matrix is symmetric up to the rounding, the symmetry is not enforced.

  The inputs are differentiated in passes over chunks of at most `:chunk-size` inputs, each pass calls `f` once with the inputs of the chunk replaced by duals. The default chunk size splits `n` inputs evenly into chunks of at most 8. Small chunks are faster for sparse functions (each operation depends on few inputs), larger for dense ones. The result does not depend on the chunk size, except at a tie of a comparison (see [[<]]): a tie is decided by the partial derivatives, so a branch taken at a tie can depend on the chunk. The same chunk size is used for the gradient and for the Jacobian.

  Parameters:

  - `f` (function): vector of numbers or duals -> number or dual.
  - `options` (map, optional):
    - `:chunk-size` - the number of inputs differentiated in one pass, an integer >= 1 (default: balanced chunks of at most 8).

  Returns a function of a point `v` (any sequence of numbers or duals, or a `double[]`) returning the matrix as a vector of vectors: of doubles, or of duals when `v` contains duals of an enclosing derivative. The Hessian of the empty point is `[]` and `f` is not called.

  Throws `ex-info` with `{:option :value :expected}` when `:chunk-size` is less than 1. The returned function throws `ex-info` when `f` returns something else than a number or a dual, and when a dual of a finished derivative call leaks into the result of an enclosing derivative. A dual of an earlier pass of the same call, reused in a later pass, is not detected and is treated as a constant.

  See also [[gradient]], [[jacobian]], [[fastmath.calculus/hessian]] (finite differences)."
  ([f] (der/hessian f nil))
  ([f options] (der/hessian f options)))
