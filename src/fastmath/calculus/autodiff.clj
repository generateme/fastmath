(ns fastmath.calculus.autodiff
  "Forward-mode automatic differentiation with tagged dual numbers.

  A dual number carries a value and its partial derivatives (perturbations), one for each input differentiated at once. Arithmetic and elementary functions defined in this namespace propagate the perturbations by the rules of calculus, so a function built from them computes its derivative alongside its value, exactly up to the floating point rounding. The design follows ForwardDiff.jl.

  Every call of a derivative creates a fresh tag and marks its perturbation with it. Duals with different tags never mix their perturbations, which makes nested derivatives (a derivative of a function which itself takes a derivative) correct and avoids perturbation confusion. A dual with a newer tag is the outer layer, its value and partials can be duals with older tags.

  A function to differentiate has to use the operations of this namespace instead of `clojure.core` or `fastmath.core` arithmetic. Numeric equality is [[==]]; `clojure.core/=` compares duals structurally (tag, value and partials), consistently with their hash, so duals can be used in sets and as map keys. All operations accept plain numbers and duals, a number is treated as a constant and coerced to a double. A `Dual` is not a `java.lang.Number`, so it can not be passed silently to a primitive function, which would drop the perturbation.

  Defined functions:

  - derivatives: [[derivative]], [[gradient]], [[jacobian]], [[hessian]],
  - derivative products: [[jvp]], [[hvp]],
  - arithmetic: [[+]], [[-]], [[*]], [[/]],
  - elementary functions: [[sin]], [[cos]], [[tan]], [[asin]], [[acos]], [[atan]], [[atan2]], [[sinh]], [[cosh]], [[tanh]], [[asinh]], [[acosh]], [[atanh]], [[exp]], [[expm1]], [[log]], [[log1p]], [[log2]], [[log10]], [[pow]], [[sqrt]], [[cbrt]], [[sq]], [[hypot]],
  - piecewise functions: [[abs]], [[min]], [[max]], [[signum]], [[floor]], [[ceil]],
  - gamma and beta functions: [[gamma]], [[log-gamma]], [[log-gamma-1p]], [[inv-gamma-1pm1]], [[digamma]], [[trigamma]], [[polygamma]], [[beta]], [[log-beta]], [[lower-incomplete-gamma]], [[upper-incomplete-gamma]], [[regularized-gamma-p]], [[regularized-gamma-q]], [[incomplete-beta]], [[regularized-beta]],
  - functions with derivatives given by the user: [[rule]],
  - comparisons: [[<]], [[<=]], [[>]], [[>=]], [[==]],
  - duals: [[dual]], [[primal]], [[scalar]].

  The values of the functions are the values of `fastmath.core` and `fastmath.special`, the derivatives follow DiffRules.jl. Where a function is not defined (its value is `NaN` for arguments which are not `NaN`, for example [[acosh]] below `1`) all its derivatives are `NaN`. At the kinks of the piecewise functions the derivatives are those of ForwardDiff.jl. The shape parameters of the incomplete gamma and beta functions and the order of [[polygamma]] have no derivative: a dual depending on the variables there throws, a constant (or a dual with zero partials) is accepted.

  Limitations: a function returning a function can not be differentiated (it throws). A dual which escapes its derivative call (stored in an atom, a cache or a lazy sequence realized later) is detected only when it flows into the result of an enclosing derivative, which then throws. An escaped dual used only in a comparison (a branch) is not detected and can silently give a wrong derivative, as can a dual left over from a finished, unrelated derivative call, or a dual of an earlier pass of the same [[gradient]], [[jacobian]], [[hessian]] or [[hvp]] call reused in a later pass: such duals have older tags and are treated as constants."
  (:refer-clojure :exclude [+ - * / < > <= >= == abs min max])
  (:require [fastmath.core :as m]
            [fastmath.calculus.common :as cc]
            [fastmath.calculus.autodiff.dual :as d]
            [fastmath.calculus.autodiff.ops :as ops]
            [fastmath.calculus.autodiff.elementary :as el]
            [fastmath.calculus.autodiff.special :as sp]
            [fastmath.calculus.autodiff.rules :as rules]
            [fastmath.calculus.autodiff.derivatives :as der])
  (:import [fastmath.calculus.autodiff.dual Dual]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; The public API. The implementation lives in the internal namespaces `fastmath.calculus.autodiff.dual` (the dual
;; type, tags, partials containers), `fastmath.calculus.autodiff.ops` (arithmetic and derivative rules) and
;; `fastmath.calculus.autodiff.derivatives` (seeding, extraction, passes and derivative products).

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

;; ## Public API: elementary functions

(defn tan
  "Tangent of a number or a dual.

  Parameters:

  - `x` (number or dual): angle in radians.

  Returns a double for a number, otherwise a dual. The derivative is `1 + tan(x)^2`.

  See also [[sin]], [[cos]], [[atan]]."
  [x] (el/tan x))

(defn asin
  "Arcus sine of a number or a dual.

  Parameters:

  - `x` (number or dual): a value in `[-1, 1]`.

  Returns a double for a number, otherwise a dual. The derivative is `1/sqrt(1 - x^2)`, infinite at `-1` and `1`; outside `[-1, 1]` the value and all derivatives are `NaN`.

  See also [[acos]], [[atan]], [[sin]]."
  [x] (el/asin x))

(defn acos
  "Arcus cosine of a number or a dual.

  Parameters:

  - `x` (number or dual): a value in `[-1, 1]`.

  Returns a double for a number, otherwise a dual. The derivative is `-1/sqrt(1 - x^2)`, infinite at `-1` and `1`; outside `[-1, 1]` the value and all derivatives are `NaN`.

  See also [[asin]], [[atan]], [[cos]]."
  [x] (el/acos x))

(defn atan
  "Arcus tangent of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `1/(1 + x^2)`.

  See also [[atan2]], [[tan]]."
  [x] (el/atan x))

(defn atan2
  "Angle of the point `(x, y)` from the positive x axis, as `Math/atan2`: the first argument is `y`.

  Parameters:

  - `y` (number or dual): the second coordinate.
  - `x` (number or dual): the first coordinate.

  Returns a double when both arguments are numbers, otherwise a dual. The partial derivatives are `x/(x^2 + y^2)` in `y` and `-y/(x^2 + y^2)` in `x`, `NaN` at the origin.

  See also [[atan]], [[hypot]]."
  [y x] (el/atan2 y x))

(defn sinh
  "Hyperbolic sine of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `cosh(x)`.

  See also [[cosh]], [[tanh]], [[asinh]]."
  [x] (el/sinh x))

(defn cosh
  "Hyperbolic cosine of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `sinh(x)`.

  See also [[sinh]], [[tanh]], [[acosh]]."
  [x] (el/cosh x))

(defn tanh
  "Hyperbolic tangent of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `1 - tanh(x)^2`.

  See also [[sinh]], [[cosh]], [[atanh]]."
  [x] (el/tanh x))

(defn asinh
  "Area hyperbolic sine (inverse of [[sinh]]) of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `1/sqrt(x^2 + 1)`.

  See also [[acosh]], [[atanh]], [[sinh]]."
  [x] (el/asinh x))

(defn acosh
  "Area hyperbolic cosine (inverse of [[cosh]]) of a number or a dual.

  Parameters:

  - `x` (number or dual): a value `>= 1`.

  Returns a double for a number, otherwise a dual. The derivative is `1/sqrt(x^2 - 1)`, infinite at `1`; below `1` the value and all derivatives are `NaN`.

  See also [[asinh]], [[atanh]], [[cosh]]."
  [x] (el/acosh x))

(defn atanh
  "Area hyperbolic tangent (inverse of [[tanh]]) of a number or a dual.

  Parameters:

  - `x` (number or dual): a value in `[-1, 1]`.

  Returns a double for a number, otherwise a dual. The derivative is `1/(1 - x^2)`, infinite at `-1` and `1`; outside `[-1, 1]` the value and all derivatives are `NaN`.

  See also [[asinh]], [[acosh]], [[tanh]]."
  [x] (el/atanh x))

(defn expm1
  "`exp(x) - 1` of a number or a dual, accurate for `x` near zero.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `exp(x)`.

  See also [[exp]], [[log1p]]."
  [x] (el/expm1 x))

(defn log1p
  "`log(1 + x)` of a number or a dual, accurate for `x` near zero.

  Parameters:

  - `x` (number or dual): a value `>= -1`.

  Returns a double for a number, otherwise a dual. The derivative is `1/(1 + x)`, infinite at `-1`; below `-1` the value and all derivatives are `NaN`.

  See also [[log]], [[expm1]]."
  [x] (el/log1p x))

(defn log2
  "Logarithm base 2 of a number or a dual.

  Parameters:

  - `x` (number or dual): a value `>= 0`.

  Returns a double for a number, otherwise a dual. The derivative is `1/(x ln 2)`, infinite at zero; for a negative argument the value and all derivatives are `NaN`.

  See also [[log]], [[log10]]."
  [x] (el/log2 x))

(defn log10
  "Logarithm base 10 of a number or a dual.

  Parameters:

  - `x` (number or dual): a value `>= 0`.

  Returns a double for a number, otherwise a dual. The derivative is `1/(x ln 10)`, infinite at zero; for a negative argument the value and all derivatives are `NaN`.

  See also [[log]], [[log2]]."
  [x] (el/log10 x))

(defn sqrt
  "Square root of a number or a dual.

  Parameters:

  - `x` (number or dual): a value `>= 0`.

  Returns a double for a number, otherwise a dual. The derivative is `1/(2 sqrt(x))`, infinite at zero; for a negative argument the value and all derivatives are `NaN`.

  See also [[cbrt]], [[sq]], [[pow]]."
  [x] (el/sqrt x))

(defn cbrt
  "Cubic root of a number or a dual, defined for negative arguments too.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `1/(3 cbrt(x)^2)`, infinite at zero.

  See also [[sqrt]], [[pow]]."
  [x] (el/cbrt x))

(defn sq
  "Square `x * x` of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `2x`.

  See also [[sqrt]], [[pow]]."
  [x] (el/sq x))

(defn hypot
  "Euclidean norm of two or three numbers or duals, without an intermediate overflow.

  Parameters:

  - `x`, `y` (numbers or duals), optionally `z` (number or dual).

  Returns a double when all arguments are numbers, otherwise a dual. The partial derivative in each argument is the argument divided by the norm, `NaN` at the origin.

  See also [[sqrt]], [[atan2]]."
  ([x y] (el/hypot x y))
  ([x y z] (el/hypot x y z)))

;; ## Public API: piecewise functions

(defn abs
  "Absolute value of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is the sign of `x` given by its sign bit, as in ForwardDiff.jl: `1` at `0.0` and `-1` at `-0.0`.

  See also [[signum]], [[max]]."
  [x] (el/abs x))

(defn signum
  "Sign of a number or a dual: `1.0`, `0.0` or `-1.0`.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The value is `0.0` for `0.0`, `-0.0` and `NaN`, as `fastmath.core/signum`. The derivative is `0` everywhere, also at zero.

  See also [[abs]]."
  [x] (el/signum x))

(defn floor
  "Rounding of a number or a dual towards negative infinity.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `0` everywhere, also at the integers.

  See also [[ceil]]."
  [x] (el/floor x))

(defn ceil
  "Rounding of a number or a dual towards positive infinity.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `0` everywhere, also at the integers.

  See also [[floor]]."
  [x] (el/ceil x))

(defn max
  "The largest of numbers and duals.

  The value is the value of `fastmath.core/max`. The partial derivatives are those of the first argument when its value is strictly larger, otherwise of the second one, as in ForwardDiff.jl: the values are compared, not the partials, so at a tie the derivative of the second argument is taken (unlike [[>]], which compares the partials at a tie). With `NaN` the value is `NaN` and the derivatives are those of the second argument.

  Parameters:

  - one or more numbers or duals, folded from the left for more than two.

  Returns a double when all arguments are numbers, otherwise a dual.

  See also [[min]], [[abs]]."
  ([x] (d/leaf x))
  ([x y] (el/max x y))
  ([x y & more] (reduce el/max (el/max x y) more)))

(defn min
  "The smallest of numbers and duals.

  The value is the value of `fastmath.core/min`. The partial derivatives are those of the second argument when the value of the first one is strictly larger, otherwise of the first one, as in ForwardDiff.jl: the values are compared, not the partials, so at a tie the derivative of the first argument is taken (unlike [[<]], which compares the partials at a tie). With `NaN` the value is `NaN` and the derivatives are those of the first argument.

  Parameters:

  - one or more numbers or duals, folded from the left for more than two.

  Returns a double when all arguments are numbers, otherwise a dual.

  See also [[max]], [[abs]]."
  ([x] (d/leaf x))
  ([x y] (el/min x y))
  ([x y & more] (reduce el/min (el/min x y) more)))

;; ## Public API: gamma and beta functions

(defn gamma
  "Gamma function of a number or a dual, see [[fastmath.special/gamma]].

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `gamma(x) digamma(x)`; at the poles (zero and the negative integers) the value and all derivatives are `NaN`.

  See also [[log-gamma]], [[digamma]], [[beta]]."
  [x] (sp/gamma x))

(defn log-gamma
  "Logarithm of the absolute value of the gamma function of a number or a dual, see [[fastmath.special/log-gamma]].

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `digamma(x)`, also for negative arguments.

  See also [[gamma]], [[log-gamma-1p]], [[digamma]]."
  [x] (sp/log-gamma x))

(defn log-gamma-1p
  "`log-gamma(1 + x)` of a number or a dual, accurate for `x` near zero, see [[fastmath.special/log-gamma-1p]].

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `digamma(1 + x)`.

  See also [[log-gamma]], [[inv-gamma-1pm1]]."
  [x] (sp/log-gamma-1p x))

(defn inv-gamma-1pm1
  "`1/gamma(1 + x) - 1` of a number or a dual, accurate for `x` near zero, see [[fastmath.special/inv-gamma-1pm1]].

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `-digamma(1 + x)/gamma(1 + x)`.

  See also [[gamma]], [[log-gamma-1p]]."
  [x] (sp/inv-gamma-1pm1 x))

(defn digamma
  "Digamma function (the derivative of [[log-gamma]]) of a number or a dual, see [[fastmath.special/digamma]].

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `trigamma(x)`.

  See also [[trigamma]], [[polygamma]], [[log-gamma]]."
  [x] (sp/digamma x))

(defn trigamma
  "Trigamma function (the derivative of [[digamma]]) of a number or a dual, see [[fastmath.special/trigamma]].

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `polygamma(2, x)`.

  See also [[digamma]], [[polygamma]]."
  [x] (sp/trigamma x))

(defn polygamma
  "Polygamma function of order `m` (the `m`-th derivative of [[digamma]]) of a number or a dual, see [[fastmath.special/polygamma]].

  Parameters:

  - `m` (number or dual): the order, a nonnegative integer (any other value gives `NaN`). It has no derivative: a dual depending on the variables throws, a constant or a dual with zero partials is accepted.
  - `x` (number or dual).

  Returns a double when both arguments are numbers, otherwise a dual. The derivative in `x` is `polygamma(m + 1, x)`.

  Throws `ex-info` with `{:function :argument}` when `m` depends on the variables.

  See also [[digamma]], [[trigamma]]."
  [m x] (sp/polygamma m x))

(defn beta
  "Beta function of numbers or duals, see [[fastmath.special/beta]].

  Parameters:

  - `p`, `q` (numbers or duals).

  Returns a double when both arguments are numbers, otherwise a dual. The partial derivatives are `beta(p, q) (digamma(p) - digamma(p + q))` in `p` and `beta(p, q) (digamma(q) - digamma(p + q))` in `q`.

  See also [[log-beta]], [[gamma]]."
  [p q] (sp/beta p q))

(defn log-beta
  "Logarithm of the beta function of numbers or duals, see [[fastmath.special/log-beta]].

  Parameters:

  - `p`, `q` (numbers or duals).

  Returns a double when both arguments are numbers, otherwise a dual. The partial derivatives are `digamma(p) - digamma(p + q)` in `p` and `digamma(q) - digamma(p + q)` in `q`.

  See also [[beta]], [[log-gamma]]."
  [p q] (sp/log-beta p q))

(defn lower-incomplete-gamma
  "Lower incomplete gamma function of numbers or duals, see [[fastmath.special/lower-incomplete-gamma]].

  Parameters:

  - `s` (number or dual): the shape parameter. It has no derivative: a dual depending on the variables throws, a constant or a dual with zero partials is accepted.
  - `x` (number or dual).

  Returns a double when both arguments are numbers, otherwise a dual. The derivative in `x` is `x^(s-1) exp(-x)`.

  Throws `ex-info` with `{:function :argument}` when `s` depends on the variables.

  See also [[upper-incomplete-gamma]], [[regularized-gamma-p]]."
  [s x] (sp/lower-incomplete-gamma s x))

(defn upper-incomplete-gamma
  "Upper incomplete gamma function of numbers or duals, see [[fastmath.special/upper-incomplete-gamma]].

  Parameters:

  - `s` (number or dual): the shape parameter. It has no derivative: a dual depending on the variables throws, a constant or a dual with zero partials is accepted.
  - `x` (number or dual).

  Returns a double when both arguments are numbers, otherwise a dual. The derivative in `x` is `-x^(s-1) exp(-x)`.

  Throws `ex-info` with `{:function :argument}` when `s` depends on the variables.

  See also [[lower-incomplete-gamma]], [[regularized-gamma-q]]."
  [s x] (sp/upper-incomplete-gamma s x))

(defn regularized-gamma-p
  "Regularized lower incomplete gamma function `P(a, x)` of numbers or duals, see [[fastmath.special/regularized-gamma-p]].

  Parameters:

  - `a` (number or dual): the shape parameter. It has no derivative: a dual depending on the variables throws, a constant or a dual with zero partials is accepted.
  - `x` (number or dual).

  Returns a double when both arguments are numbers, otherwise a dual. The derivative in `x` is `x^(a-1) exp(-x) / gamma(a)`, computed without overflow for large `a`.

  Throws `ex-info` with `{:function :argument}` when `a` depends on the variables.

  See also [[regularized-gamma-q]], [[lower-incomplete-gamma]]."
  [a x] (sp/regularized-gamma-p a x))

(defn regularized-gamma-q
  "Regularized upper incomplete gamma function `Q(a, x)` of numbers or duals, see [[fastmath.special/regularized-gamma-q]].

  Parameters:

  - `a` (number or dual): the shape parameter. It has no derivative: a dual depending on the variables throws, a constant or a dual with zero partials is accepted.
  - `x` (number or dual).

  Returns a double when both arguments are numbers, otherwise a dual. The derivative in `x` is `-x^(a-1) exp(-x) / gamma(a)`, computed without overflow for large `a`.

  Throws `ex-info` with `{:function :argument}` when `a` depends on the variables.

  See also [[regularized-gamma-p]], [[upper-incomplete-gamma]]."
  [a x] (sp/regularized-gamma-q a x))

(defn incomplete-beta
  "Incomplete beta function `B(x; a, b)` of numbers or duals, see [[fastmath.special/incomplete-beta]]. The first argument is `x`.

  Parameters:

  - `x` (number or dual): a value in `[0, 1]`.
  - `a`, `b` (numbers or duals): the shape parameters. They have no derivative: a dual depending on the variables throws, a constant or a dual with zero partials is accepted.

  Returns a double when all arguments are numbers, otherwise a dual. The derivative in `x` is `x^(a-1) (1-x)^(b-1)`.

  Throws `ex-info` with `{:function :argument}` when `a` or `b` depends on the variables.

  See also [[regularized-beta]], [[beta]]."
  [x a b] (sp/incomplete-beta x a b))

(defn regularized-beta
  "Regularized incomplete beta function `I_x(a, b)` of numbers or duals, see [[fastmath.special/regularized-beta]]. The first argument is `x`.

  Parameters:

  - `x` (number or dual): a value in `[0, 1]`.
  - `a`, `b` (numbers or duals): the shape parameters. They have no derivative: a dual depending on the variables throws, a constant or a dual with zero partials is accepted.

  Returns a double when all arguments are numbers, otherwise a dual. The derivative in `x` is `x^(a-1) (1-x)^(b-1) / beta(a, b)`, computed without overflow for large `a` and `b`.

  Throws `ex-info` with `{:function :argument}` when `a` or `b` depends on the variables.

  See also [[incomplete-beta]], [[beta]]."
  [x a b] (sp/regularized-beta x a b))

;; ## Public API: rules of the user

(defn rule
  "Creates a function with a derivative given by the user, which accepts numbers and duals like the functions of this namespace.

  Use it for a function which is not built from the operations of this namespace (a primitive function, a library call), when its derivative is known. The derivative is first order only: the created function works in [[derivative]] of order `1`, [[gradient]], [[jacobian]] and [[jvp]], and throws in nested derivatives ([[derivative]] of order `2` or more, [[hessian]], [[hvp]]), which would need the derivative of the derivative.

  Parameters:

  - `f` (function): the value, a function of one or two doubles returning a number.
  - `df` (function or `nil`): the derivative of a function of one argument, a function of a double returning a number.
  - `dfx`, `dfy` (functions or `nil`): the partial derivatives of a function of two arguments in the first and the second argument, functions of two doubles returning numbers.

  A derivative given as `nil` means that the function has no derivative in this argument: a dual depending on the variables there throws, a constant (or a dual with zero partials) is accepted.

  Returns a function of one argument (for `f` and `df`) or two arguments (for `f`, `dfx` and `dfy`). Called with numbers only, it returns `(f x)` or `(f x y)` with the arguments coerced to doubles. Called with duals, the value and the derivatives are evaluated at the values of the arguments and the result is a dual; where the value is `NaN` for arguments which are not `NaN`, the derivatives are `NaN`.

  Throws `ex-info` with `{:option :value :expected}` when `f` or a derivative is neither a function nor `nil` (`f` can not be `nil`). The created function throws `ex-info` with `{:tag}` in a nested derivative, and `ex-info` with `{:function :argument}` for a dual depending on the variables in an argument without a derivative.

  See also [[derivative]], [[gradient]]."
  ([f df]
   (when-not (ifn? f) (cc/invalid-option! :f f "a function"))
   (when-not (or (nil? df) (ifn? df)) (cc/invalid-option! :df df "a function or nil"))
   (rules/rule1 f df))
  ([f dfx dfy]
   (when-not (ifn? f) (cc/invalid-option! :f f "a function"))
   (when-not (or (nil? dfx) (ifn? dfx)) (cc/invalid-option! :dfx dfx "a function or nil"))
   (when-not (or (nil? dfy) (ifn? dfy)) (cc/invalid-option! :dfy dfy "a function or nil"))
   (rules/rule2 f dfx dfy)))

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

(defn jvp
  "Creates the Jacobian-vector product (the directional derivative) of a multivariate function `f` with forward-mode automatic differentiation.

  For a scalar function the product is the directional derivative, the gradient times the direction; for a vector valued function it is the Jacobian matrix times the direction. It is computed in one pass, without the gradient or the Jacobian: `f` is called once with every input `x_i` replaced by a dual with the single partial `v_i`. The result is exact up to the floating point rounding, it can differ from the product of [[gradient]] or [[jacobian]] and the direction in the last bits (a different order of the summation). Products can be nested with other derivatives: every call uses a fresh tag.

  Parameters:

  - `f` (function): vector of `n` numbers or duals -> number or dual, or a sequential of `m` numbers or duals.

  Returns a function of a point `x` and a direction `v` (both any sequence of `n` numbers or duals, or a `double[]`) returning the product: a double for a scalar `f`, a vector of `m` elements for a vector valued `f`; duals when `x` or `v` contain duals of an enclosing derivative. For the empty point `f` is called once and the product is `0.0` (or a vector of zeros).

  The returned function throws `ex-info` with `{:count-x :count-v}` when `x` and `v` have different lengths, and `ex-info` when `f` returns something else than a number, a dual or a sequential of them (a function, a `double[]`), and when a dual of a finished derivative call leaks into the result of an enclosing derivative.

  See also [[hvp]], [[gradient]], [[jacobian]]."
  [f] (der/jvp f))

(defn hvp
  "Creates the Hessian-vector product of a multivariate scalar function `f` with forward-mode automatic differentiation.

  The product is computed as the Jacobian-vector product ([[jvp]]) of the [[gradient]] (nested forward mode), without the Hessian matrix: its cost is about one gradient of duals, the cost of the full [[hessian]] is about `n / chunk-size` of them. The result is exact up to the floating point rounding, it can differ from the product of [[hessian]] and the direction in the last bits.

  Parameters:

  - `f` (function): vector of numbers or duals -> number or dual.
  - `options` (map, optional): options of the inner [[gradient]]
    - `:chunk-size` - the number of inputs differentiated in one pass, an integer >= 1 (default: balanced chunks of at most 8).

  Returns a function of a point `x` and a direction `v` (both any sequence of `n` numbers or duals, or a `double[]`) returning the product as a vector of `n` elements: doubles, or duals when `x` or `v` contain duals of an enclosing derivative. The product for the empty point is `[]` and `f` is not called.

  Throws `ex-info` with `{:option :value :expected}` when `:chunk-size` is less than 1. The returned function throws as [[jvp]] and [[gradient]] do.

  See also [[jvp]], [[hessian]], [[gradient]]."
  ([f] (der/hvp f nil))
  ([f options] (der/hvp f options)))
