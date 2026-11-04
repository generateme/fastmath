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
  - reciprocal and inverse trigonometric and hyperbolic functions: [[sec]], [[csc]], [[cot]], [[sech]], [[csch]], [[coth]], [[asec]], [[acsc]], [[acot]], [[asech]], [[acsch]], [[acoth]],
  - pi-scaled functions: [[sinpi]], [[cospi]], [[tanpi]], [[cotpi]], [[secpi]], [[cscpi]], [[sinc]],
  - exponential and logarithmic variants: [[exp2]], [[exp10]], [[logb]], [[ln]], [[xlogx]], [[xlogy]], [[xlog1py]], [[log1pexp]], [[logexpm1]], [[log1mexp]], [[log1pmx]], [[log1psq]], [[logcosh]], [[exprel]], [[expexp]], [[cexpexp]], [[cloglog]], [[loglog]],
  - sigmoid and log-sum-exp: [[sigmoid]], [[logit]], [[logaddexp]], [[logsubexp]], [[logsumexp]], [[smooth-max]],
  - piecewise functions: [[abs]], [[min]], [[max]], [[signum]], [[floor]], [[ceil]], [[trunc]], [[rint]], [[round]], [[sgn]], [[frac]], [[sfrac]], [[mod]], [[rem]], [[copy-sign]],
  - powers and convenience functions: [[cb]], [[pow2]], [[pow3]], [[pow10]], [[spow]], [[safe-sqrt]], [[fma]], [[muladd]], [[lerp]], [[smoothstep]], [[norm]],
  - gamma and beta functions: [[gamma]], [[log-gamma]], [[log-gamma-1p]], [[inv-gamma-1pm1]], [[digamma]], [[trigamma]], [[polygamma]], [[beta]], [[log-beta]], [[lower-incomplete-gamma]], [[upper-incomplete-gamma]], [[regularized-gamma-p]], [[regularized-gamma-q]], [[incomplete-beta]], [[regularized-beta]],
  - functions with derivatives given by the user: [[rule]],
  - comparisons: [[<]], [[<=]], [[>]], [[>=]], [[==]],
  - duals: [[dual]], [[primal]], [[scalar]].

  The values of the functions are the values of `fastmath.core` and `fastmath.special`, the derivatives follow DiffRules.jl. Where a function is not defined (its value is `NaN` for arguments which are not `NaN`, for example [[acosh]] below `1`) all its derivatives are `NaN`. At the kinks of the piecewise functions the derivatives are those of ForwardDiff.jl, at the jumps (of [[frac]], [[sfrac]], [[mod]], [[rem]]) they are `NaN`. The shape parameters of the incomplete gamma and beta functions, the order of [[polygamma]] and the sharpness of [[smooth-max]] have no derivative: a dual depending on the variables there throws, a constant (or a dual with zero partials) is accepted.

  Limitations: a function returning a function can not be differentiated (it throws). A dual which escapes its derivative call (stored in an atom, a cache or a lazy sequence realized later) is detected only when it flows into the result of an enclosing derivative, which then throws. An escaped dual used only in a comparison (a branch) is not detected and can silently give a wrong derivative, as can a dual left over from a finished, unrelated derivative call, or a dual of an earlier pass of the same [[gradient]], [[jacobian]], [[hessian]] or [[hvp]] call reused in a later pass: such duals have older tags and are treated as constants."
  (:refer-clojure :exclude [+ - * / < > <= >= == abs min max mod rem])
  (:require [fastmath.core :as m]
            [fastmath.calculus.common :as cc]
            [fastmath.calculus.autodiff.dual :as d]
            [fastmath.calculus.autodiff.ops :as ops]
            [fastmath.calculus.autodiff.elementary :as el]
            [fastmath.calculus.autodiff.special :as sp]
            [fastmath.calculus.autodiff.logexp :as lx]
            [fastmath.calculus.autodiff.piecewise :as pw]
            [fastmath.calculus.autodiff.arithmetic :as ar]
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
  "Natural logarithm of a number or a dual, or the logarithm with a base.

  Parameters:

  - `x` (number or dual).
  - `base` (number or dual, optional): the base, the same as [[logb]].

  Returns a double for numbers, otherwise a dual. The value is `NaN` for a negative argument and `-Inf` for zero; the derivative is `1/x`, infinite at zero and `NaN` (as all higher derivatives) for a negative argument, outside the domain. With a base the partial derivatives are those of [[logb]].

  See also [[exp]], [[pow]], [[logb]], [[ln]]."
  ([x] (ops/log x))
  ([base x] (lx/logb base x)))

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

  Returns a double when both arguments are numbers, otherwise a dual. The partial derivatives are `x/(x^2 + y^2)` in `y` and `-y/(x^2 + y^2)` in `x`, computed without overflow or underflow of the squares; they are `NaN` at the origin and `0` (the limit) at an infinite argument.

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

  Returns a double when all arguments are numbers, otherwise a dual. The partial derivative in each argument is the argument divided by the norm, `NaN` at the origin. At an infinite argument the limits are used: the sign of the only infinite argument, `0` for the finite ones, `NaN` when several arguments are infinite.

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

  The value is the value of `fastmath.core/max`. The partial derivatives are those of the first argument when its value is strictly larger, otherwise of the second one, as in ForwardDiff.jl: the values are compared, not the partials, so at a tie the derivative of the second argument is taken (unlike [[>]], which compares the partials at a tie). With `NaN` the value is `NaN` and the derivatives are those of the second argument. Zeros of different signs are a tie: the value is the second argument, `(max 0.0 -0.0)` is `-0.0`.

  Parameters:

  - one or more numbers or duals, folded from the left for more than two.

  Returns a double when all arguments are numbers, otherwise a dual.

  See also [[min]], [[abs]]."
  ([x] (d/leaf x))
  ([x y] (el/max x y))
  ([x y & more] (reduce el/max (el/max x y) more)))

(defn min
  "The smallest of numbers and duals.

  The value is the value of `fastmath.core/min`. The partial derivatives are those of the second argument when the value of the first one is strictly larger, otherwise of the first one, as in ForwardDiff.jl: the values are compared, not the partials, so at a tie the derivative of the first argument is taken (unlike [[<]], which compares the partials at a tie). With `NaN` the value is `NaN` and the derivatives are those of the first argument. Zeros of different signs are a tie: the value is the second argument, `(min 0.0 -0.0)` is `-0.0`.

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

;; ## Public API: reciprocal and inverse trigonometric and hyperbolic functions

(defn sec
  "Secant of a number or a dual.

  Parameters:

  - `x` (number or dual): angle in radians.

  Returns a double for a number, otherwise a dual. The derivative is `sec(x) tan(x)`.

  See also [[csc]], [[cot]], [[cos]]."
  [x] (el/sec x))

(defn csc
  "Cosecant of a number or a dual.

  Parameters:

  - `x` (number or dual): angle in radians.

  Returns a double for a number, otherwise a dual. The derivative is `-csc(x) cot(x)`, infinite at the multiples of `pi`.

  See also [[sec]], [[cot]], [[sin]]."
  [x] (el/csc x))

(defn cot
  "Cotangent of a number or a dual.

  Parameters:

  - `x` (number or dual): angle in radians.

  Returns a double for a number, otherwise a dual. The derivative is `-(1 + cot(x)^2)`.

  See also [[tan]], [[csc]]."
  [x] (el/cot x))

(defn sech
  "Hyperbolic secant of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `-tanh(x) sech(x)`.

  See also [[csch]], [[coth]], [[cosh]]."
  [x] (el/sech x))

(defn csch
  "Hyperbolic cosecant of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `-coth(x) csch(x)`, infinite at `0`.

  See also [[sech]], [[coth]], [[sinh]]."
  [x] (el/csch x))

(defn coth
  "Hyperbolic cotangent of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `-csch(x)^2`, infinite at `0`.

  See also [[tanh]], [[csch]]."
  [x] (el/coth x))

(defn asec
  "Arcus secant of a number or a dual.

  Parameters:

  - `x` (number or dual): a value with `|x| >= 1`.

  Returns a double for a number, otherwise a dual. The derivative is `1/(|x| sqrt(x^2 - 1))`, infinite at `-1` and `1`; inside `(-1, 1)` the value and all derivatives are `NaN`.

  See also [[acsc]], [[acos]]."
  [x] (el/asec x))

(defn acsc
  "Arcus cosecant of a number or a dual.

  Parameters:

  - `x` (number or dual): a value with `|x| >= 1`.

  Returns a double for a number, otherwise a dual. The derivative is `-1/(|x| sqrt(x^2 - 1))`, infinite at `-1` and `1`; inside `(-1, 1)` the value and all derivatives are `NaN`.

  See also [[asec]], [[asin]]."
  [x] (el/acsc x))

(defn acot
  "Arcus cotangent of a number or a dual, `pi/2 - atan(x)`.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `-1/(1 + x^2)`.

  See also [[atan]], [[cot]]."
  [x] (el/acot x))

(defn asech
  "Area hyperbolic secant of a number or a dual.

  Parameters:

  - `x` (number or dual): a value in `(0, 1]`.

  Returns a double for a number, otherwise a dual. The derivative is `-1/(x sqrt(1 - x^2))`, infinite at `0` and `1`; outside the domain the value and all derivatives are `NaN`.

  See also [[acosh]], [[sech]]."
  [x] (el/asech x))

(defn acsch
  "Area hyperbolic cosecant of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `-1/(|x| sqrt(1 + x^2))`, infinite at `0`.

  See also [[asinh]], [[csch]]."
  [x] (el/acsch x))

(defn acoth
  "Area hyperbolic cotangent of a number or a dual.

  Parameters:

  - `x` (number or dual): a value with `|x| > 1`.

  Returns a double for a number, otherwise a dual. The derivative is `1/(1 - x^2)`; inside `[-1, 1]` the value is `NaN` or infinite.

  See also [[atanh]], [[coth]]."
  [x] (el/acoth x))

;; ## Public API: pi-scaled functions and sinc

(defn sinpi
  "Sine of `pi x` for a number or a dual.

  Parameters:

  - `x` (number or dual): angle in half-turns.

  Returns a double for a number, otherwise a dual. The value is `fastmath.core/sinpi`, `sin(pi x)`. The derivative is `pi cospi(x)`.

  See also [[cospi]], [[sin]]."
  [x] (el/sinpi x))

(defn cospi
  "Cosine of `pi x` for a number or a dual.

  Parameters:

  - `x` (number or dual): angle in half-turns.

  Returns a double for a number, otherwise a dual. The derivative is `-pi sinpi(x)`.

  See also [[sinpi]], [[cos]]."
  [x] (el/cospi x))

(defn tanpi
  "Tangent of `pi x` for a number or a dual.

  Parameters:

  - `x` (number or dual): angle in half-turns.

  Returns a double for a number, otherwise a dual. The derivative is `pi (1 + tanpi(x)^2)`.

  See also [[cotpi]], [[tan]]."
  [x] (el/tanpi x))

(defn cotpi
  "Cotangent of `pi x` for a number or a dual.

  Parameters:

  - `x` (number or dual): angle in half-turns.

  Returns a double for a number, otherwise a dual. The derivative is `-pi (1 + cotpi(x)^2)`.

  See also [[tanpi]], [[cot]]."
  [x] (el/cotpi x))

(defn secpi
  "Secant of `pi x` for a number or a dual.

  Parameters:

  - `x` (number or dual): angle in half-turns.

  Returns a double for a number, otherwise a dual. The derivative is `pi secpi(x) tanpi(x)`.

  See also [[cscpi]], [[sec]]."
  [x] (el/secpi x))

(defn cscpi
  "Cosecant of `pi x` for a number or a dual.

  Parameters:

  - `x` (number or dual): angle in half-turns.

  Returns a double for a number, otherwise a dual. The derivative is `-pi cscpi(x) cotpi(x)`.

  See also [[secpi]], [[csc]]."
  [x] (el/cscpi x))

(defn sinc
  "Normalized sinc of a number or a dual, `sin(pi x)/(pi x)`, `1` at `0`.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `(cos(pi x) - sinc(x))/x`; near `0`, where this form cancels, a Taylor series is used, so the derivatives at `0` are the limits (`0`, `-pi^2/3`, ...).

  See also [[sinpi]]."
  [x] (el/sinc x))

;; ## Public API: exponential and logarithmic variants

(defn exp2
  "Power of two, `2^x`, of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `2^x ln(2)`.

  See also [[exp10]], [[exp]], [[log2]]."
  [x] (lx/exp2 x))

(defn exp10
  "Power of ten, `10^x`, of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `10^x ln(10)`.

  See also [[exp2]], [[exp]], [[log10]]."
  [x] (lx/exp10 x))

(defn logb
  "Logarithm of `x` with the base `b`, `log(x)/log(b)`, of numbers and duals.

  Parameters:

  - `b` (number or dual): the base, positive and not `1`.
  - `x` (number or dual): a positive value.

  Returns a double when all arguments are numbers, otherwise a dual. The partial derivatives are `-log(x)/(b log(b)^2)` in `b` and `1/(x log(b))` in `x`.

  See also [[log]], [[log2]], [[log10]]."
  [b x] (lx/logb b x))

(defn ln
  "Natural logarithm of a number or a dual, the same as [[log]] with one argument.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `1/x`.

  See also [[log]], [[exp]]."
  [x] (ops/log x))

(defn xlogx
  "`x log(x)` of a number or a dual, `0` at `0`.

  Parameters:

  - `x` (number or dual): a value `>= 0`.

  Returns a double for a number, otherwise a dual. The derivative is `1 + log(x)`, `-Inf` at `0`; for a negative argument the value and all derivatives are `NaN`.

  See also [[xlogy]], [[log]]."
  [x] (lx/xlogx x))

(defn xlogy
  "`x log(y)` of numbers and duals, `0` when `x` is `0` and `y` is not `NaN`.

  Parameters:

  - `x` (number or dual).
  - `y` (number or dual): a value `>= 0`.

  Returns a double when all arguments are numbers, otherwise a dual. The partial derivatives are `log(y)` in `x` and `x/y` in `y`; the latter is `0` when `x` is zero together with all its derivatives and `y` is not `NaN`, as in LogExpFunctions.jl.

  See also [[xlogx]], [[xlog1py]]."
  [x y] (lx/xlogy x y))

(defn xlog1py
  "`x log(1 + y)` of numbers and duals, `0` when `x` is `0` and `y` is not `NaN`.

  Parameters:

  - `x` (number or dual).
  - `y` (number or dual): a value `>= -1`.

  Returns a double when all arguments are numbers, otherwise a dual. The partial derivatives are `log(1 + y)` in `x` and `x/(1 + y)` in `y`; the latter is `0` when `x` is zero together with all its derivatives and `y` is not `NaN`, as in LogExpFunctions.jl.

  See also [[xlogy]], [[log1p]]."
  [x y] (lx/xlog1py x y))

(defn log1pexp
  "`log(1 + exp(x))` (softplus) of a number or a dual, without overflow.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is [[sigmoid]] of `x`.

  See also [[logexpm1]], [[log1mexp]], [[logaddexp]]."
  [x] (lx/log1pexp x))

(defn logexpm1
  "`log(exp(x) - 1)` of a number or a dual, the inverse of [[log1pexp]].

  Parameters:

  - `x` (number or dual): a positive value.

  Returns a double for a number, otherwise a dual. The derivative is `exp(x - logexpm1(x))`; for `x <= 0` the value and all derivatives are `NaN` (or infinite at `0`).

  See also [[log1pexp]], [[expm1]]."
  [x] (lx/logexpm1 x))

(defn log1mexp
  "`log(1 - exp(x))` of a number or a dual, accurate near `0`.

  Parameters:

  - `x` (number or dual): a negative value.

  Returns a double for a number, otherwise a dual. The derivative is `-exp(x - log1mexp(x))`; for a positive argument the value and all derivatives are `NaN`.

  See also [[log1pexp]], [[logsubexp]]."
  [x] (lx/log1mexp x))

(defn log1pmx
  "`log(1 + x) - x` of a number or a dual, accurate near `0`.

  Parameters:

  - `x` (number or dual): a value `> -1`.

  Returns a double for a number, otherwise a dual. The derivative is `-x/(1 + x)`.

  See also [[log1p]]."
  [x] (lx/log1pmx x))

(defn log1psq
  "`log(1 + x^2)` of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `2x/(1 + x^2)` (`0` beyond `|x| = 1.3e154`, where `x^2` overflows).

  See also [[log1p]]."
  [x] (lx/log1psq x))

(defn logcosh
  "`log(cosh(x))` of a number or a dual, without overflow.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `tanh(x)`.

  See also [[cosh]], [[tanh]]."
  [x] (lx/logcosh x))

(defn exprel
  "`(exp(x) - 1)/x` of a number or a dual, `1` at `0`.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `(exp(x) - exprel(x))/x`; near `0`, where this form cancels, a Taylor series is used (the derivative at `0` is `1/2`), and above it for a positive `x` the form `exprel(x) (1 - 1/x) + 1/x`.

  See also [[expm1]], [[exp]]."
  [x] (lx/exprel x))

(defn expexp
  "`exp(-exp(-x))` (the Gumbel distribution function) of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `exp(-x - exp(-x))`.

  See also [[cexpexp]], [[loglog]]."
  [x] (lx/expexp x))

(defn cexpexp
  "`1 - exp(-exp(x))` of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `exp(x - exp(x))`.

  See also [[expexp]], [[cloglog]]."
  [x] (lx/cexpexp x))

(defn cloglog
  "Complementary log-log, `log(-log(1 - x))`, of a number or a dual.

  Parameters:

  - `x` (number or dual): a value in `(0, 1)`.

  Returns a double for a number, otherwise a dual. The derivative is `-1/((1 - x) log(1 - x))`; outside `[0, 1]` the value and all derivatives are `NaN`.

  See also [[loglog]], [[cexpexp]]."
  [x] (lx/cloglog x))

(defn loglog
  "Log-log, `-log(-log(x))`, of a number or a dual.

  Parameters:

  - `x` (number or dual): a value in `(0, 1)`.

  Returns a double for a number, otherwise a dual. The derivative is `-1/(x log(x))`; outside `[0, 1]` the value and all derivatives are `NaN`.

  See also [[cloglog]], [[expexp]]."
  [x] (lx/loglog x))

;; ## Public API: sigmoid and log-sum-exp

(defn sigmoid
  "Logistic sigmoid, `1/(1 + exp(-x))`, of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `sigmoid(x) (1 - sigmoid(x))`.

  See also [[logit]], [[log1pexp]]."
  [x] (lx/sigmoid x))

(defn logit
  "Logit, `log(x/(1 - x))`, of a number or a dual, the inverse of [[sigmoid]].

  Parameters:

  - `x` (number or dual): a value in `(0, 1)`.

  Returns a double for a number, otherwise a dual. The derivative is `1/(x (1 - x))`; outside `[0, 1]` the value and all derivatives are `NaN`.

  See also [[sigmoid]]."
  [x] (lx/logit x))

(defn logaddexp
  "`log(exp(x) + exp(y))` of numbers and duals, without overflow.

  Parameters:

  - `x` (number or dual).
  - `y` (number or dual).

  Returns a double when all arguments are numbers, otherwise a dual. The partial derivatives are `exp(x - z)` and `exp(y - z)`, `z` the value.

  See also [[logsubexp]], [[logsumexp]]."
  [x y] (lx/logaddexp x y))

(defn logsubexp
  "`log(|exp(x) - exp(y)|)` of numbers and duals.

  Parameters:

  - `x` (number or dual).
  - `y` (number or dual).

  Returns a double when all arguments are numbers, otherwise a dual. With `z` the value, the partial derivatives are `exp(x - z)` and `-exp(y - z)` when `x > y`, otherwise `-exp(x - z)` and `exp(y - z)`, as in LogExpFunctions.jl; at `x = y` the value is `-Inf` and the derivatives are infinite or `NaN`.

  See also [[logaddexp]], [[log1mexp]]."
  [x y] (lx/logsubexp x y))

(defn logsumexp
  "`log(sum(exp(xs)))` of a sequence of numbers and duals, without overflow.

  Parameters:

  - `xs` (non-empty sequence of numbers and duals).

  Returns a double when all elements are numbers, otherwise a dual. The partial derivatives are the softmax, `exp(x_i - z)`, `z` the value; with exactly one element `+Inf` they are its limit, `1` in that element and `0` in the others (as LogExpFunctions.jl), and `NaN` with more than one.

  Throws `ex-info` with `{:option :value :expected}` for an empty sequence (`fastmath.core/logsumexp` returns a value there).

  See also [[logaddexp]], [[smooth-max]]."
  [xs] (lx/logsumexp xs))

;; ## Public API: rounding and kinks

(defn trunc
  "Rounding of a number or a dual towards zero.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `0` everywhere, also at the integers.

  See also [[floor]], [[ceil]], [[rint]], [[frac]]."
  [x] (pw/trunc x))

(defn rint
  "Rounding of a number or a dual to the nearest integer, half to even.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `0` everywhere, also at the half-integers.

  See also [[round]], [[trunc]]."
  [x] (pw/rint x))

(defn round
  "Rounding of a number or a dual to the nearest integer, half up, as a double.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The value is `fastmath.core/round` as a double: it saturates at the `long` range (`(round 1e20)` is `9.223372036854776E18`) and is `0.0` for `NaN`. The derivative is `0` everywhere.

  See also [[rint]], [[trunc]]."
  [x] (pw/round x))

(defn sgn
  "Sign of a number or a dual: `-1.0` for a negative value, otherwise `1.0` (also for `0.0`, `-0.0` and `NaN`).

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `0` everywhere, also at zero.

  See also [[signum]], [[copy-sign]]."
  [x] (pw/sgn x))

(defn frac
  "Unsigned fractional part of a number or a dual, `|x - trunc(x)|`.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is the sign of `x` by its sign bit (`1` at `0.0`, `-1` at `-0.0`, as [[abs]]) and `NaN` at the nonzero integers (the jumps); for an infinite argument the value and the derivatives are `NaN`.

  See also [[sfrac]], [[trunc]], [[mod]]."
  [x] (pw/frac x))

(defn sfrac
  "Signed fractional part of a number or a dual, `x - trunc(x)`.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `1`, `NaN` at the nonzero integers (the jumps).

  See also [[frac]], [[rem]]."
  [x] (pw/sfrac x))

(defn mod
  "Floored modulus of numbers and duals, the sign of the divisor.

  Parameters:

  - `x` (number or dual): the dividend.
  - `y` (number or dual): the divisor.

  Returns a double when all arguments are numbers, otherwise a dual. The partial derivatives are `1` in `x` and `-floor(x/y)` in `y`; both are `NaN` when `x/y` is an integer (the jumps), as in DiffRules.jl.

  See also [[rem]], [[frac]]."
  [x y] (pw/mod x y))

(defn rem
  "Truncated remainder of numbers and duals, the sign of the dividend.

  Parameters:

  - `x` (number or dual): the dividend.
  - `y` (number or dual): the divisor.

  Returns a double when all arguments are numbers, otherwise a dual. The partial derivatives are `1` in `x` and `-trunc(x/y)` in `y`; both are `NaN` when `x/y` is an integer (the jumps), as in DiffRules.jl.

  See also [[mod]], [[sfrac]]."
  [x y] (pw/rem x y))

(defn copy-sign
  "The magnitude of the first argument with the sign of the second one, for numbers and duals.

  Parameters:

  - `magnitude` (number or dual).
  - `sign` (number or dual): only its sign bit is used.

  Returns a double when all arguments are numbers, otherwise a dual. The partial derivatives are `1` in `magnitude` when the sign bits of the arguments agree and `-1` otherwise, and `0` in `sign`.

  See also [[abs]], [[sgn]]."
  [magnitude sign] (pw/copy-sign magnitude sign))

;; ## Public API: powers and convenience functions

(defn cb
  "Cube of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `3x^2`.

  See also [[pow3]], [[sq]], [[pow]]."
  [x] (ar/cb x))

(defn pow2
  "Square of a number or a dual, the same as [[sq]].

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `2x`.

  See also [[sq]], [[pow]]."
  [x] (ar/pow2 x))

(defn pow3
  "Cube of a number or a dual, the same as [[cb]].

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `3x^2`.

  See also [[cb]], [[pow]]."
  [x] (ar/cb x))

(defn pow10
  "Tenth power of a number or a dual.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `10x^9`.

  See also [[pow]]."
  [x] (ar/pow10 x))

(defn spow
  "Symmetric power `sgn(x) |x|^e` of numbers and duals, defined for a negative base.

  Parameters:

  - `x` (number or dual): the base.
  - `e` (number or dual): the exponent.

  Returns a double when all arguments are numbers, otherwise a dual. The partial derivatives are `e |x|^(e-1)` in `x` (infinite at `0` for `e < 1`) and `spow(x, e) log|x|` in `e`, `0` at `x = 0` for `e > 0` (the limit, as [[pow]] with a constant base).

  See also [[pow]], [[abs]]."
  [x e] (pw/spow x e))

(defn safe-sqrt
  "Square root of a number or a dual, `0` for a negative value.

  Parameters:

  - `x` (number or dual).

  Returns a double for a number, otherwise a dual. The derivative is `1/(2 sqrt(x))`, and `0` for a negative argument (the value is the constant `0`).

  See also [[sqrt]]."
  [x] (pw/safe-sqrt x))

(defn fma
  "Fused multiply-add `x y + z` of numbers and duals.

  Parameters:

  - `x` (number or dual).
  - `y` (number or dual).
  - `z` (number or dual).

  Returns a double when all arguments are numbers, otherwise a dual. The partial derivatives are `y`, `x` and `1`.

  See also [[muladd]], [[*]], [[+]]."
  [x y z] (ar/fma x y z))

(defn muladd
  "Multiply-add `x y + z` of numbers and duals, the same as [[fma]].

  Parameters:

  - `x` (number or dual).
  - `y` (number or dual).
  - `z` (number or dual).

  Returns a double when all arguments are numbers, otherwise a dual. The partial derivatives are `y`, `x` and `1`.

  See also [[fma]]."
  [x y z] (ar/fma x y z))

(defn lerp
  "Linear interpolation `start + t (stop - start)` of numbers and duals.

  Parameters:

  - `start` (number or dual).
  - `stop` (number or dual).
  - `t` (number or dual): the amount, extrapolates outside `[0, 1]`.

  Returns a double when all arguments are numbers, otherwise a dual. The partial derivatives are `1 - t`, `t` and `stop - start`.

  See also [[norm]], [[smoothstep]]."
  [start stop t] (ar/lerp start stop t))

(defn smoothstep
  "Smooth Hermite step of numbers and duals: `0` up to `edge0`, `1` from `edge1`, `3t^2 - 2t^3` between, `t = (x - edge0)/(edge1 - edge0)`.

  Parameters:

  - `edge0` (number or dual).
  - `edge1` (number or dual).
  - `x` (number or dual).

  Returns a double when all arguments are numbers, otherwise a dual. Between the edges the partial derivatives are `-s (1 - t)`, `-s t` and `s` with `s = 6t(1 - t)/(edge1 - edge0)`; outside (where `t` is clamped, also at the edges) they are `0`.

  See also [[lerp]], [[norm]]."
  [edge0 edge1 x] (pw/smoothstep edge0 edge1 x))

(defn norm
  "Linear map of numbers and duals from one range to another.

  Parameters:

  - `v` (number or dual): the value to map.
  - `start`, `stop` (numbers or duals): the source range, mapped to `[0, 1]` (3 arguments).
  - `start1`, `stop1`, `start2`, `stop2` (numbers or duals): the source and the target ranges (5 arguments).

  Returns a double when all arguments are numbers, otherwise a dual. With `t = (v - start)/(stop - start)` the partial derivatives are `1/w`, `(t - 1)/w` and `-t/w` in `v`, `start` and `stop`, `w = stop - start`; with 5 arguments these are multiplied by `stop2 - start2` and the derivatives in `start2` and `stop2` are `1 - t` and `t`. For equal ends of the source range the value is a constant end (as `fastmath.core/norm`) and the derivatives in the source range and `v` are `0`.

  See also [[lerp]], [[smoothstep]]."
  ([v start stop] (ar/norm v start stop))
  ([v start1 stop1 start2 stop2] (ar/norm v start1 stop1 start2 stop2)))

(defn smooth-max
  "Smooth maximum of a sequence of numbers and duals, `fastmath.core/smooth-max`.

  Parameters:

  - `xs` (non-empty sequence of numbers and duals).
  - `alpha` (number, default `2.0`): the sharpness; it has no derivative.
  - `family` (keyword, default `:lse`): the smoothing family, see `fastmath.core/smooth-max`.

  Returns a double when all elements and `alpha` are numbers, otherwise a dual. Derivatives are available for the families `:lse`, `:mellowmax` (both the softmax of `alpha xs`), `:boltzmann`, `:p-norm` (`sign(x_i) (|x_i|/z)^(alpha-1)`, `z` the value) and `:smu` (a fold of the pairwise rule).

  Throws `ex-info` with `{:option :value :expected}` for an empty sequence (`fastmath.core/smooth-max` returns a value there). For duals throws `ex-info` with `{:function :family}` for the piecewise families (`:sigmoid`, `:circular`, `:quadratic`, `:cubic`, `:quartic`, `:exponential`), and with `{:function :argument}` when `alpha` depends on the variables.

  See also [[logsumexp]], [[max]]."
  ([xs] (lx/smooth-max xs 2.0 :lse))
  ([xs alpha] (lx/smooth-max xs alpha :lse))
  ([xs alpha family] (lx/smooth-max xs alpha family)))

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

  Throws `ex-info` with `{:option :value :expected}` when `f` or a derivative is neither a function (`fn?`: not a keyword, a map, a vector or a var) nor `nil` (`f` can not be `nil`). The created function throws `ex-info` with `{:tag}` in a nested derivative, and `ex-info` with `{:function :argument}` for a dual depending on the variables in an argument without a derivative.

  See also [[derivative]], [[gradient]]."
  ([f df]
   (when-not (fn? f) (cc/invalid-option! :f f "a function"))
   (when-not (or (nil? df) (fn? df)) (cc/invalid-option! :df df "a function or nil"))
   (rules/rule1 f df))
  ([f dfx dfy]
   (when-not (fn? f) (cc/invalid-option! :f f "a function"))
   (when-not (or (nil? dfx) (fn? dfx)) (cc/invalid-option! :dfx dfx "a function or nil"))
   (when-not (or (nil? dfy) (fn? dfy)) (cc/invalid-option! :dfy dfy "a function or nil"))
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

  The inputs are differentiated in passes over chunks of at most `:chunk-size` inputs, each pass calls `f` once with the inputs of the chunk replaced by duals. The default chunk size splits `n` inputs evenly into chunks of at most 8. Small chunks are faster for sparse functions (each operation depends on few inputs), larger for dense ones. The result does not depend on the chunk size, except at a tie of a comparison (see [[<]]): a tie is decided by the partial derivatives, so a branch taken at a tie can depend on the chunk. Also where a partial derivative is infinite: in one pass the infinite factor multiplies the zero partials of the other inputs of the chunk, `0 * Inf` is `NaN` (as in ForwardDiff.jl), so the derivatives in these inputs are `NaN`, while with `:chunk-size` `1` they are not; for example the gradient of `sqrt(x) y` at `[0 1]` is `[Inf NaN]`, or `[Inf 0]` with one input per pass.

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

  The inputs are differentiated in passes over chunks of at most `:chunk-size` inputs, each pass calls `f` once with the inputs of the chunk replaced by duals. The default chunk size splits `n` inputs evenly into chunks of at most 8. Small chunks are faster for sparse functions (each operation depends on few inputs), larger for dense ones. The result does not depend on the chunk size, except at a tie of a comparison (see [[<]]): a tie is decided by the partial derivatives, so a branch taken at a tie can depend on the chunk. Also where a partial derivative is infinite: in one pass the infinite factor multiplies the zero partials of the other inputs of the chunk, `0 * Inf` is `NaN` (as in ForwardDiff.jl), so the derivatives in these inputs are `NaN`, while with `:chunk-size` `1` they are not; for example the gradient of `sqrt(x) y` at `[0 1]` is `[Inf NaN]`, or `[Inf 0]` with one input per pass.

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

  The inputs are differentiated in passes over chunks of at most `:chunk-size` inputs, each pass calls `f` once with the inputs of the chunk replaced by duals. The default chunk size splits `n` inputs evenly into chunks of at most 8. Small chunks are faster for sparse functions (each operation depends on few inputs), larger for dense ones. The result does not depend on the chunk size, except at a tie of a comparison (see [[<]]): a tie is decided by the partial derivatives, so a branch taken at a tie can depend on the chunk. Also where a partial derivative is infinite: in one pass the infinite factor multiplies the zero partials of the other inputs of the chunk, `0 * Inf` is `NaN` (as in ForwardDiff.jl), so the derivatives in these inputs are `NaN`, while with `:chunk-size` `1` they are not; for example the gradient of `sqrt(x) y` at `[0 1]` is `[Inf NaN]`, or `[Inf 0]` with one input per pass. The same chunk size is used for the gradient and for the Jacobian.

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

  For a scalar function the product is the directional derivative, the gradient times the direction; for a vector valued function it is the Jacobian matrix times the direction. It is computed in one pass, without the gradient or the Jacobian: `f` is called once with every input `x_i` replaced by a dual with the single partial `v_i`. A partial derivative which is infinite in an input with a zero direction gives `NaN` (`0 * Inf`), as in ForwardDiff.jl. The result is exact up to the floating point rounding, it can differ from the product of [[gradient]] or [[jacobian]] and the direction in the last bits (a different order of the summation). Products can be nested with other derivatives: every call uses a fresh tag.

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
