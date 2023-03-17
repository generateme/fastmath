^{:nextjournal.clerk/visibility :hide-ns
  :nextjournal.clerk/toc true}
(ns core
  (:require [fastmath.core :as m]
            [nextjournal.clerk :as clerk]
            [clojure.walk :as walk]
            [utils :as u]))

;; # fastmath.core

;; Collection of type hinted math macros and functions. Partially backed by Java static functions and exposed as macros. They are prepared to accept primitive `long` or `double` arguments and return `long` or `double` only. There is no support for Clojure specific numeric types (like eg. Ratio or Number).

;; There is a possibility to replace `clojure.core` functions with a selection of `fastmath.core` macros. Call:

;; * `(m/use-primitive-operators)` to replace functions with macros
;; * `(m/unuse-primitive-operators)` to revert replacement.

;; Be aware that there are some differences and `fastmath.core` versions shoudn't be treated as a drop-in replacement. Also, since Clojure 1.12, always call `unuse-primitive-operators` at the end of the namespace.

;; Here is the complete list of replaced functions:

;; * `* + - /`
;; * `> < >= <= ==`
;; * `rem quot mod`
;; * `bit-or bit-and bit-xor bit-not bit-and-not`
;; * `bit-shift-left bit-shift-right unsigned-bit-shift-right`
;; * `bit-set bit-clear bit-flip bit-test`
;; * `inc dec`
;; * `zero? neg? pos? even? odd?`
;; * `min max`
;; * `abs`

(require '[fastmath.core :as m])

;; ## Primitive math ops

;; Primitive math version of some mathematical clojure functions and operations.

;; ### add,sub,mul,div

;; Multi-arity macros for primitive mathematical operations.

^{::clerk/visibility :hide}
(u/table
 [[+ true]
  [- true "negation for 1-ary"]
  [* true]
  [/ true "reciprocal for 1-ary, doesn't return Ratio"]])

;; Please note that division with mixed argument types can lead to unexpected result (comparing to Clojure). Operations are done pairwise and division acts as `quot` for integer arguments.

;; n-ary version is expanded to a pairwise execution

(clerk/code (macroexpand '(m// 1 2 3.1 -11.1)))

^{::clerk/visibility :hide}
(clerk/example
 (m/+ 1)
 (m/+ -4 2)
 (m/+ 1 2 3.132 4 5)
 (m/- 1)
 (m/- -4 2)
 (m/- 1 2 3.132 4 5)
 (m/* 1)
 (m/* -4 2)
 (m/* 1 2 3.132 4 5)
 (m// 2)
 (m// -4 2)
 (m// -2 3)
 (m// 1.0 2.0 3.132 4.0 5.0)
 (m// 1 2 4.0))

;; ### Integer division and remainders

^{::clerk/visibility :hide}
(u/table
 [[quot true "same as in Clojure"]
  [mod true "same as in Clojure"]
  [rem true "same as in Clojure"]
  [remainder true "see note below"]
  [wrap false "wraps a value to be within the range if it overflows"]])

;; * `remainder` macro returns $dividend - divisor * n$, where $n$ is the mathematical integer closest to $\frac{dividend}{divisor}$. Returned value is inside the $[\frac{-|divisor|}{2},\frac{|divisor|}{2}]$ range

^{::clerk/visibility :hide}
(clerk/example
 (m/quot 10 4)
 (m/quot -10.25 4.0)
 (m/quot 10.25 -4.0)
 (m/quot -10.25 -4.0)
 (m/mod 10 4)
 (m/mod -10.25 4.0)
 (m/mod 10.25 -4.0)
 (m/mod -10.25 -4.0)
 (m/rem 10 4)
 (m/rem -10.25 4.0)
 (m/rem 10.25 -4.0)
 (m/rem -10.25 -4.0)
 (m/remainder 10 4)
 (m/remainder -10.25 4.0)
 (m/remainder 10.25 -4.0)
 (m/remainder -10.25 -4.0)
 (m/wrap -1.25 1.25 1.0)
 (m/wrap -1.25 1.25 1.35)
 (m/wrap [-1.25 1.25] -1.35))

;; ### GCD and LCM

^{::clerk/visibility :hide}
(u/table
 [[gcd false "greatest common divisor, Stein's algorithm"]
  [lcm false "least common multiplier as (/ (* a b) (gcd a b))"]])

^{::clerk/visibility :hide}
(clerk/example
 (m/gcd 10 115)
 (m/lcm 10 115))

;; ### Incrementation and decrementation

^{::clerk/visibility :hide}
(u/table
 [[inc true "same as in Clojure"]
  [dec true "same as in Clojure"]])

^{::clerk/visibility :hide}
(clerk/example
 (m/inc 1)
 (m/inc 1.25)
 (m/dec 1)
 (m/dec 1.25))

;; ### Comparison

^{::clerk/visibility :hide}
(u/table
 [[== true "numerical equality"]
  [eq false "function, up to 4 doubles"]
  [not== true "numerical inequality"]
  [> true]
  [< true]
  [>= true]
  [<= true]
  [zero? true]
  [one? true]
  [pos? true]
  [not-pos? true]
  [neg? true]
  [not-neg? true]
  [even? true]
  [odd? true]
  [between? false "is value in closed interval?"]
  [between-? false "is value in left open interval?"]])

^{::clerk/visibility :hide}
(clerk/example
 (m/== 1 2 3 4)
 (m/eq 1.0 2.0 3.0 4.0)
 (m/not== 1 2 3 4)
 (m/< 1 2 3 4)
 (m/> 1 2 3 4)
 (m/>= 1 2 3 4)
 (m/<= 1 2 3 4)
 (m/zero? 0)
 (m/zero? 0.0)
 (m/one? 1.0)
 (m/one? -1)
 (m/pos? 0)
 (m/not-pos? 0)
 (m/neg? 0)
 (m/not-neg? 0.0)
 (m/even? 2)
 (m/odd? 2)
 (m/even? 2.9))

;; Also two additional functions if value is between ranges

;; * `between?` - closed range, $v\in[x,y]$
;; * `between-?` - left open range, $v\in(x,y]$

;; Both function can accept a two values vector as the range.

^{::clerk/visibility :hide}
(clerk/example
 (m/between? 2.0 10.0 2.0)
 (m/between? [2.0 10.0] 2.0)
 (m/between-? 2.0 10.0 2.0)
 (m/between-? [2.0 10.0] 2.0))

;; #### Approximate equality

;; There are two options:

^{::clerk/visibility :hide}
(u/table
 [[approx-eq false "see notes below"]
  [approx= false "the same as approx-eq"]
  [delta-eq false "see notes below"]
  [delta= false "the same as delta-eq"]])

;; * `approx-eq` or `approx=` - which first round numbers to specified number of decimal places (default `2`), then compare rounded numbers
;; * `delta-eq` or `delta=` - which verifies if absolute difference between two numbers is less than given accuracy (default: `1.0e-6`)

;; These functions work on doubles only.

^{::clerk/visibility :hide}
(clerk/example
 (m/approx-eq 2.230 2.231)
 (m/approx-eq 2.230 2.231 3)
 (m/delta-eq 2.230 2.231)
 (m/delta-eq 2.230 2.231 1.0e-2))

;; #### NaN and Inf

;; Set of function to test against `NAN` and `Inf` values. 

^{::clerk/visibility :hide}
(u/table
 [[nan? false "checks if arg is NaN"]
  [inf? false "checks if arg is either Inf or -Inf"]
  [pos-inf? false "check if arg is Inf"]
  [neg-inf? false "check if arg is -Inf"]
  [invalid-double? false "is arg Inf, -Inf or NaN?"]
  [valid-double? false "is not arg Inf, -Inf or NaN?"]])

^{::clerk/visibility :hide}
(clerk/example
 (m/nan? ##NaN)
 (m/inf? ##-Inf)
 (m/neg-inf? ##-Inf)
 (m/pos-inf? ##-Inf)
 (map m/invalid-double? [##Inf ##-Inf ##NaN 1.0])
 (map m/valid-double? [##Inf ##-Inf ##NaN 1.0]))

;; ### Sign

^{::clerk/visibility :hide}
(u/table
 [[sgn false "returns -1.0 for negative numbers, 1.0 otherwise"]
  [signum false "returns -1.0 for negative numbers, 0.0 for zero, 1.0 for positive numbers"]
  [copy-sign true "sets sign of second argument to the first argument"]])

^{::clerk/visibility :hide}
(clerk/example
 (m/sgn -1.0)
 (m/sgn 0.0)
 (m/sgn 1.0)
 (m/signum -1.0)
 (m/signum 0.0)
 (m/signum 1.0)
 (m/copy-sign 1.0 1.0)
 (m/copy-sign -1.0 1.0)
 (m/copy-sign 1.0 -1.0)
 (m/copy-sign -1.0 -1.0))

;; ### Other

^{::clerk/visibility :hide}
(u/table
 [[abs false "works on doubles only"]
  [iabs false "works on longs only"]
  [max true]
  [min true]])

^{::clerk/visibility :hide}
(clerk/example
 (m/abs -1.0)
 (m/iabs -1)
 (m/max 1 2 3 4.0)
 (m/min 1 2 3 4.0))

;; ### Primitive ops as functions

;; Some functions are exposed as two arity inlined functions to use for fast reduction. They work only on doubles.

^{::clerk/visibility :hide}
(u/table
 [[fast+ false]
  [fast- false]
  [fast* false]
  [fast-max false]
  [fast-min false]
  [fast-identity false "returns input as double"]])

^{::clerk/visibility :hide}
(clerk/example
 (m/fast+ 2.25 3.325)
 (m/fast- 2.25 3.325)
 (m/fast* 2.25 3.325)
 (m/fast-max 2.25 3.325)
 (m/fast-min 2.25 3.325)
 (m/fast-identity 2.25))

;; ## Bitwise operations

;; Bit manipulation

^{::clerk/visibility :hide}
(u/table
 [[bit-shift-left true "same as <<"]
  [<< true "same as bit-shift-left"]
  [bit-shift-right true "same as >>"]
  [>> true "same as bit-shift-right"]
  [unsigned-bit-shift-right true "same as >>>"]
  [>>> true "same as unsigned-bit-shift-right"]
  [bit-or true "a | b"]
  [bit-and true "a & b"]
  [bit-xor true "a ^ b"]
  [bit-not true "~a"]
  [bit-and-not true "a & ~b"]
  [bit-set true "a | (1 << n)"]
  [bit-clear true "a & ~(1 << n)"]
  [bit-flip true "a & (1 << n)"]
  [bit-test true]])

^{::clerk/visibility :hide}
(clerk/example
 (m/bit-shift-left -122 2)
 (m/<< -122 2)
 (m/bit-shift-right -122 2)
 (m/>> -122 2)
 (m/unsigned-bit-shift-right -122 2)
 (m/>>> -122 2)
 (m/bit-or 123 -243)
 (m/bit-and 123 -234)
 (m/bit-and-not 123 -234)
 (m/bit-xor 123 -234)
 (m/bit-not 123)
 (m/bit-set 123 2)
 (m/bit-clear 123 1)
 (m/bit-flip 123 5)
 (m/bit-test 123 1))

;; ## Rounding

^{::clerk/visibility :hide}
(u/table
 [[floor false "round to the lower integer, value can be scaled optionally, returns double"]
  [ceil false "round to the upper integer, value can be scaled optionally, returns double"]
  [round false "round to the nearest integer, returns long"]
  [round-even false "IEEE / IEC rounding, returns long"]
  [rint false "round to the nearest integer, returns double"]
  [qfloor true "floor by long casting with correction"]
  [qceil true "ceil by long casting with correction"]
  [qround true "round by long casting with correction"]
  [trunc false "truncate fractional part, returns double"]
  [itrunc false "truncate fractionl part, returns long"]
  [approx false "round fractional part to the given number of decimals (default: 2)"]
  [frac false "fractional part, unsigned"]
  [sfrac false "fractional part, signed"]
  [low-2-exp false "find greatest exponent for which 2^x is less than argument"]
  [high-2-exp false "find lowest exponent for which 2^x is greater than argument"]
  [low-exp false "find greatest exponent for which b^x is less than argument"]
  [high-exp false "find lowest exponent for which b^x is greater than argument"]
  [round-up-pow2 false "round long to the next (nearest) power of two"]])

^{::clerk/visibility :hide}
(clerk/example
 (map m/floor [-10.5 10.5])
 (m/floor 10.5 4.0)
 (map m/ceil [-10.5 10.5])
 (m/ceil 10.5 4.0)
 (map m/round [-10.51 -10.5 -10.49 10.49 10.5 10.51])
 (map m/round-even [-10.51 -10.5 -10.49 10.49 10.5 10.51])
 (map m/rint [-10.51 -10.5 -10.49 10.49 10.5 10.51])
 (map (fn [x] (m/qfloor x)) [-10.5 10.5])
 (map (fn [x] (m/qceil x)) [-10.5 10.5])
 (map (fn [x] (m/qround x)) [-10.51 -10.5 -10.49 10.49 10.5 10.51])
 (map m/trunc [-10.591 10.591])
 (map m/itrunc [-10.591 10.591])
 (m/approx 10.591)
 (m/approx 10.591 1)
 (m/approx 10.591 0)
 (m/approx -10.591)
 (m/approx -10.591 1)
 (m/approx -10.591 0)
 (map m/frac [-10.591 10.591])
 (map m/sfrac [-10.591 10.591])
 (m/low-2-exp 10.591)
 (m/high-2-exp 10.591)
 (m/low-exp 0.5 10.591)
 (m/high-exp 0.5 10.591)
 (map m/round-up-pow2 (range 10)))

^{::clerk/visibility :hide ::clerk/viewer u/unpaginated-table}
[["floor" "ceil" "round"]
 [(u/fgraph m/floor) (u/fgraph m/ceil) (u/fgraph m/round)]
 ["trunc" "frac" "sfrac"]
 [(u/fgraph m/trunc) (u/fgraph m/frac) (u/fgraph m/sfrac)]]

;; ## Polynomials and fma

;; Set of macros or functions which deal with polynomial evaluation plus some useful primitive ops shortcuts.

^{::clerk/visibility :hide}
(u/table
 [[muladd true "(+ (* x y) z) or Math/fma for JDK 9+"]
  [fma true "same as muladd"]
  [negmuladd true (clerk/code '(- z (* x y)))]
  [difference-of-products false (clerk/code '(- (* a b) (* c d)))]
  [sum-of-products false (clerk/code '(+ (* a b) (* c d)))]
  [mevalpoly true "evaluate polynomial in the form: coeffs[0]+coeffs[1]*x+coeffs[2]*x^2+..."]
  [evalpoly false "same as mevalpoly but function"]
  [makepoly false "creates polynomial function"]])

;; `sum-of-products` and `difference-of-products` uses some tricks (Kahan's algorithm) to protect against catastrophic cancellation. See more [here](https://pharr.org/matt/blog/2019/11/03/difference-of-floats)

^{::clerk/visibility :hide}
(clerk/example
 (m/muladd 10.01 20.02 30.03)
 (m/fma 10.01 20.02 30.03)
 (m/negmuladd 10.01 20.2 30.3)
 (m/difference-of-products 100.0 -200.1 500.0 -400.1)
 (m/sum-of-products 100.0 -200.1 500.0 -400.1))

;; Let's define a polynomial $P(x)=0.5x^4-0.2x^3-2.0x-0.5$

(def Px (m/makepoly [-0.5 -2.0 0.0 -0.2 0.5]))

^{::clerk/visibility :hide}
(u/fgraph Px [-1.5 2] nil)

^{::clerk/visibility :hide}
(clerk/example
 (Px -1.0)
 (Px 1.0)
 (m/mevalpoly -1.0 -0.5 -2.0 0.0 -0.2 0.5)
 (m/mevalpoly 1.0 -0.5 -2.0 0.0 -0.2 0.5)
 (m/evalpoly -1.0 -0.5 -2.0 0.0 -0.2 0.5)
 (m/evalpoly 1.0 -0.5 -2.0 0.0 -0.2 0.5))

;; `mevalpoly` unrolls into calls to `fma` (for JDK 9+)

(clerk/code (walk/macroexpand-all '(m/mevalpoly -2.0 2.5 1.0 0.0 7.0 3.0)))

;; ## Trigonometric 

;; ### Basic

^{::clerk/visibility :hide}
(u/table
 [[sin true (clerk/md "$\\sin(x)$")]
  [cos true (clerk/md "$\\cos(x)$")]
  [tan true (clerk/md "$\\frac{\\sin(x)}{\\cos(x)}$")]
  [cot false (clerk/md "$\\frac{1}{\\tan(x)}$")]
  [sec false (clerk/md "$\\frac{1}{\\cos(x)}$")]
  [csc false (clerk/md "$\\frac{1}{\\sin(x)}$")]
  [asin true (clerk/md "$\\arcsin(x)$")]
  [acos true (clerk/md "$\\arccos(x)$")]
  [atan true (clerk/md "$\\arctan(x)$")]
  [acot false (clerk/md "$\\frac{\\pi}{2}-\\arctan(x)$")]
  [asec false (clerk/md "$\\arccos(\\frac{1}{x})$")]
  [acsc false (clerk/md "$\\arcsin(\\frac{1}{x})$")]
  [(atan2 y x) true "the angle between the positive x-axis of a plane and the point (x,y)"]
  [qsin true "faster and less accurate sin"]
  [qcos true "faster and less accurate cos"]])

^{::clerk/visibility :hide ::clerk/viewer u/unpaginated-table}
[["sin" "cos" "tan"]
 [(u/fgraph m/sin) (u/fgraph m/cos) (u/fgraph m/tan)]
 ["cot" "sec" "csc"]
 [(u/fgraph m/cot) (u/fgraph m/sec) (u/fgraph m/csc)]
 ["asin" "acos" "atan"]
 [(u/fgraph m/asin [-1.1 1.1] nil) (u/fgraph m/acos [-1.1 1.1] [-0.5 nil]) (u/fgraph m/atan)]
 ["acot" "asec" "acsc"]
 [(u/fgraph m/acot) (u/fgraph m/asec) (u/fgraph m/acsc)]]

;; ### Hyperbolic

^{::clerk/visibility :hide}
(u/table
 [[sinh true (clerk/md "$\\frac{e^x-e^{-x}}{2}$")]
  [cosh true (clerk/md "$\\frac{e^x+e^{-x}}{2}$")]
  [tanh true (clerk/md "$\\frac{\\sinh(x)}{\\cosh(x)}$")]
  [coth false (clerk/md "$\\frac{1}{\\tanh(x)}$")]
  [sech false (clerk/md "$\\frac{1}{\\cosh(x)}$")]
  [csch false (clerk/md "$\\frac{1}{\\sinh(x)}$")]
  [asinh true (clerk/md "$\\operatorname{arsinh}(x)$")]
  [acosh true (clerk/md "$\\operatorname{arcosh}(x)$")]
  [atanh true (clerk/md "$\\operatorname{artanh}(x)$")]
  [acoth false (clerk/md "$\\operatorname{artanh}(\\frac{1}{x})$")]
  [asech false (clerk/md "$\\operatorname{arcosh}(\\frac{1}{x})$")]
  [acsch false (clerk/md "$\\operatorname{arsinh}(\\frac{1}{x})$")]])

^{::clerk/visibility :hide ::clerk/viewer u/unpaginated-table}
[["sinh" "cosh" "tanh"]
 [(u/fgraph m/sinh) (u/fgraph m/cosh) (u/fgraph m/tanh)]
 ["coth" "sech" "csch"]
 [(u/fgraph m/coth) (u/fgraph m/sech) (u/fgraph m/csch)]
 ["asinh" "acosh" "atanh"]
 [(u/fgraph m/asinh) (u/fgraph m/acosh) (u/fgraph m/atanh)]
 ["acoth" "asech" "acsch"]
 [(u/fgraph m/acoth) (u/fgraph m/asech) (u/fgraph m/acsch)]]

;; ### Historical

^{::clerk/visibility :hide}
(u/table
 [[crd false (clerk/md "$2\\sin(\\frac{x}{2})$")]
  [versin false (clerk/md "$1-\\cos(x)$")]
  [coversin false (clerk/md "$1-\\sin(x)$")]
  [vercos false (clerk/md "$1+\\cos(x)$")]
  [covercos false (clerk/md "$1+\\sin(x)$")]
  [haversin false (clerk/md "$\\frac{1-\\cos(x)}{2}$, see below")]
  [hacoversin false (clerk/md "$\\frac{1-\\sin(x)}{2}$")]
  [havercos false (clerk/md   "$\\frac{1+\\cos(x)}{2}$")]
  [hacovercos false (clerk/md "$\\frac{1+\\sin(x)}{2}$")]
  [exsec false (clerk/md "$\\sec(x)-1$")]
  [excsc false (clerk/md "$\\csc(x)-1$")]
  [acrd false (clerk/md "$2\\arcsin(\\frac{x}{2})$")]
  [aversin false (clerk/md "$\\arccos(1-x)$")]
  [acoversin false (clerk/md "$\\arcsin(1-x)$")]
  [avercos false (clerk/md "$\\arccos(x-1)$")]
  [acovercos false (clerk/md "$\\arcsin(x-1)$")]
  [ahaversin false (clerk/md   "$\\arccos(1-2x)$")]
  [ahacoversin false (clerk/md "$\\arcsin(1-2x)$")]
  [ahavercos false (clerk/md   "$\\arccos(2x-1)$")]
  [ahacovercos false (clerk/md "$\\arcsin(2x-1)$")]
  [aexsec false (clerk/md "$\\operatorname{arcsec}(x+1)$")]
  [aexcsc false (clerk/md "$\\operatorname{arccsc}(x+1)$")]])

;; Additionally there is a special case of `haversin` which accepts longitude and lattitude.

^{::clerk/visibility :hide}
(u/table
 [[haversine false "accepts lat/lon (in radians) pairs"]
  [haversine-dist false "caluclates distance between lat/lon (in radians) pairs"]])

^{::clerk/visibility :hide}
(clerk/example
 (m/haversine [0.3 0.3] [0.5 0.5])
 (m/haversine 0.3 0.3 0.5 0.5)
 (m/haversine-dist [0.3 0.3] [0.5 0.5])
 (m/haversine-dist 0.3 0.3 0.5 0.5))

;; ### Special

^{::clerk/visibility :hide}
(u/table
 [[sinc false (clerk/md "$\\frac{\\sin(x)}{x}$")]
  [Si false (clerk/md "$\\operatorname{Si}(x)=\\int_{0}^{x}\\frac{\\sin(x)}{x}dx$")]
  [Ci false (clerk/md "$\\operatorname{Ci}(x)=-\\int_{x}^{\\infty}\\frac{\\cos(x)}{x}dx$")]])

^{::clerk/visibility :hide ::clerk/viewer u/unpaginated-table}
[["sinc" "Si" "Ci"]
 [(u/fgraph m/sinc) (u/fgraph m/Si [-10 10] nil) (u/fgraph m/Ci [0 10] [-4.0 nil])]]

;; ## Power, roots and log

;; ### Power and roots

^{::clerk/visibility :hide}
(u/table
 [[pow true (clerk/md "$x^a$")]
  [qpow true "fast and less accurate pow"]
  [fpow true "fast pow for integer exponents"]
  [sq false (clerk/md "$x^2$")]
  [cb false (clerk/md "$x^3$")]
  [pow2 false "same as sq"]
  [pow3 false "same as cb"]
  [sqrt true (clerk/md "$\\sqrt{x}$")]
  [safe-sqrt false "returns 0 for negative values"]
  [qsqrt true "approximated sqrt with max relative error of about 3.41e-2"]
  [rqsqrt true "fast inverse square root, zero step Newton's method"]
  [cbrt true (clerk/md "$\\sqrt[3]{x}$")]])

^{::clerk/visibility :hide ::clerk/viewer u/unpaginated-table}
[["sqrt" "qsqrt"]
 [(u/fgraph m/sqrt [0 5] nil) (u/fgraph m/qsqrt [0 5] nil)]
 ["(pow x 0.5) " "(qpow x 0.5)"]
 [(u/fgraph #(m/pow % 0.5) [0 5] nil) (u/fgraph #(m/pow % 0.5) [0 5] nil)]
 ["1/sqrt" "rqsqrt"]
 [(u/fgraph #(/ 1.0 (m/sqrt %)) [0.09 5] [-0.5 nil]) (u/fgraph m/rqsqrt [0.09 5] [-0.5 nil])]] 

;; ### Exp and log

;; Basic functions

^{::clerk/visibility :hide}
(u/table
 [[exp true (clerk/md "$e^x$")]
  [log true (clerk/md "$\\ln(x)$")]
  [ln true (clerk/md "$\\ln(x)$")]
  [log2 false (clerk/md "$\\log_{2}(x)$")]
  [log10 true (clerk/md "$\\log_{10}(x)$")]
  [logb false (clerk/md "$\\log_b(x)$")]
  [qexp true "fast and less accurate exp"]
  [qlog true "fast and less accurate log"]])

^{::clerk/visibility :hide ::clerk/viewer u/unpaginated-table}
[["exp" "log"]
 [(u/fgraph m/exp [-5.0 1.0] nil) (u/fgraph m/log [0.01 5.0] nil)]
 ["qexp" "qlog"]
 [(u/fgraph m/qexp [-5.0 1.0] nil) (u/fgraph m/qlog [0.01 5.0] nil)]]

;; Various additional special functions based of log and exp. Some of them are optimized. Source [LogExpFunctions from Julia](https://juliastats.org/LogExpFunctions.jl/stable/)

^{::clerk/visibility :hide}
(u/table
 [[expm1 true (clerk/md "$e^x-1$")]
  [xexpx false (clerk/md "$xe^x$")]
  [xexpy false (clerk/md "$xe^y$")]
  [cexpexp false (clerk/md "$1-e^{-e^x}$")]
  [sigmoid false (clerk/md "$\\frac{1}{1+e^{-x}}$")]
  [logistic false "same as sigmoid"]
  [log1p true (clerk/md "$\\ln(1+x)$")]
  [log1pexp false (clerk/md "$\\ln(1+e^x)$")]
  [log1mexp false (clerk/md "$\\ln(1-e^x)$")]
  [log2mexp false (clerk/md "$\\ln(2-e^x)$")]
  [log1psq false (clerk/md "$\\ln(1+x^2)$")]
  [logexpm1 false (clerk/md "$\\ln(e^x-1)$")]
  [log1pmx false (clerk/md "$\\ln(1+x)-x$")]
  [logmxp1 false (clerk/md "$\\ln(x)-x+1$")]
  [logaddexp false (clerk/md "$\\ln(e^x+e^y)$")]
  [logsubexp false (clerk/md "$\\ln(e^x-e^y)$")]
  [logsumexp false (clerk/md "$\\ln(e^{x_1}+\\dots+e^{x_n})$")]
  [xlogx false (clerk/md "$x\\ln(x)$")]
  [xlogy false (clerk/md "$x\\ln(y)$")]
  [xlog1py false (clerk/md "$x\\ln(1+y)$")]
  [cloglog false (clerk/md "$\\ln(-\\ln(1-x))$")]
  [logit false (clerk/md "$\\ln(\\frac{x}{1-x})$")]
  [logcosh false (clerk/md "$\\ln(\\cosh(x))$")]])

^{::clerk/visibility :hide ::clerk/viewer u/unpaginated-table}
[["expm1" "xexpx" "cexpexp"]
 [(u/fgraph m/expm1) (u/fgraph m/xexpx) (u/fgraph m/cexpexp)]
 ["log1p" "log1pexp" "log1mexp"]
 [(u/fgraph m/log1p) (u/fgraph m/log1pexp) (u/fgraph m/log1mexp)]
 ["log2mexp" "log1psq" "logexpm1"]
 [(u/fgraph m/log2mexp) (u/fgraph m/log1psq) (u/fgraph m/logexpm1)]
 ["log1pmx" "xlogx" "cloglog"]
 [(u/fgraph m/log1pmx) (u/fgraph m/xlogx) (u/fgraph m/cloglog)]
 ["sigmoid" "logit" "logcosh"]
 [(u/fgraph m/sigmoid [-6 6] [-0.2 1.2]) (u/fgraph m/logit [-0.2 1.2] nil) (u/fgraph m/logcosh)]]

;; ## Distance

^{::clerk/visibility :hide}
(u/table
 [[dist false (clerk/md "$\\sqrt{(x_1-x_2)^2+(y_1-y_2)^2}$")]
  [qdist false (clerk/md "same as `dist` but `qsqrt` is used")]
  [hypot false (clerk/md "$\\sqrt{x^2+y^2}$ or $\\sqrt{x^2+y^2+z^2}$ without over/underflow")]
  [hypot-sqrt false (clerk/md "$\\sqrt{x^2+y^2}$ or $\\sqrt{x^2+y^2+z^2}$")]
  [haversine-dist false "spherical distance"]])

^{::clerk/visibility :hide}
(clerk/example
 (m/dist 1.1 2.2 3.3 4.4)
 (m/dist [1.1 2.2] [3.3 4.4])
 (m/qdist 1.1 2.2 3.3 4.4)
 (m/qdist [1.1 2.2] [3.3 4.4])
 (m/hypot 1.0 (m/pow 3 512))
 (m/hypot-sqrt 1.0 (m/pow 3 512)))

;; ## Special

^{::clerk/visibility :hide}
(u/table
 [[erf true "error function, for two arguments it's a difference between erf(x) and erf(y)"]
  [erfc true "complementary error function"]
  [inv-erf true "inverse error function"]
  [inv-erfc true "inverse complementary error function"]
  [gamma true "gamma(x) function"]
  [log-gamma true "log of gamma(x)"]
  [log-gamma-1p true "log of gamma(x+1) for -0.5<x<1.5"]
  [digamma true "derivative of log of gamma(x)"]
  [trigamma true "derivative of digamma"]
  [log-gamma-1p true "log of gamma(x+1) for -0.5<x<1.5"]
  [inv-gamma-1pm1 true (clerk/md "$\\frac{1}{\\Gamma(x+1)}-1$ for $x\\in(-0.5,1.5)$")]
  [regularized-gamma-p true "regularized gamma P"]
  [regularized-gamma-q true "regularized gamma Q"]
  [log-beta true "log of beta function"]
  [regularized-beta true "regularized beta function"]
  [bessel-j true "Bessel J function for given order and argument"]
  [jinc false "Bessel J of order 1 divided by x"]
  [I0 false "Modified Bessel function of the first kind and order 0"]
  [log-I0 false "log of I0"]
  [minkowski false "Minkowski's question mark function, ?(x)"]])

^{::clerk/visibility :hide ::clerk/viewer u/unpaginated-table}
[["erf" "inv-erf"]
 [(u/fgraph m/erf) (u/fgraph m/inv-erf [-1.1 1.1] nil)]
 ["gamma" "log-gamma"] 
 [(u/fgraph m/gamma [-5.0 3.0]) (u/fgraph m/log-gamma [-0.2 3.0])]
 ["digamma" "trigamma"]
 [(u/fgraph m/digamma [-0.2 3.0] [-5 1.1]) (u/fgraph m/trigamma [-0.2 3.0])]
 ["bessel-j 0" "bessel-j 1"]
 [(u/fgraph #(m/bessel-j 0 %) [0 5.0] nil) (u/fgraph #(m/bessel-j 1 %) [0 5.0] nil)]
 ["I0" "log-I0"]
 [(u/fgraph m/I0) (u/fgraph m/log-I0)]
 ["jinc" "minkowski"]
 [(u/fgraph m/jinc [0 5] nil) (u/fgraph m/minkowski [-0.1 1.1])]]

;; ## Interpolation

;; Various interpolations on $[x_1,x_2]$ interval, $t\in[0,1]$:

^{::clerk/visibility :hide}
(u/table
 [[lerp false (clerk/md "Linear, $\\operatorname{lerp}(t)=x_1+(x_2-x_1)t$")]
  [mlerp true "macro version of lerp"]
  [cos-interpolation false (clerk/md "$f(x)=x_1+(x_2-x_1)\\frac{1-\\cos(\\pi t)}{2}$")]
  [smooth-interpolation false (clerk/md "$f(x)=x_1+3(x_2-x_1)(t^2-2t^3)$")]
  [quad-interpolation false "quadratic interpolation"]])

^{::clerk/visibility :hide ::clerk/viewer u/unpaginated-table}
[["lerp" "cos" "smooth" "quad"]
 [(u/fgraph (partial m/lerp 0 1) [0 1] nil)
  (u/fgraph (partial m/cos-interpolation 0 1) [0 1] nil)
  (u/fgraph (partial m/smooth-interpolation 0 1) [0 1] nil)
  (u/fgraph (partial m/quad-interpolation 0 1) [0 1] nil)]]

^{::clerk/visibility :hide}
(clerk/example
 (map (partial m/lerp -1.0 1.0) (range 0.0 1.01 0.25))
 (map (partial m/cos-interpolation -1.0 1.0) (range 0.0 1.01 0.25))
 (map (partial m/smooth-interpolation -1.0 1.0) (range 0.0 1.01 0.25))
 (map (partial m/quad-interpolation -1.0 1.0) (range 0.0 1.01 0.25)))

;; ## Mapping

;; Mappings and conversions. Please note that:

;; * `smoothstep` expects value as the last argument (as in GLSL)
;; * `norm` with variants expect value as the first arguemnt (as `map` in Processing)

^{::clerk/visibility :hide}
(u/table
 [[constrain true "clamp value to a given range, (max (min value mx) mn)"]
  [norm false "map lineary value from given interval to a new interval (or [0,1] by default)"]
  [mnorm true "macro version of norm"]
  [cnorm false "norm which clapms result to a target interval"]
  [make-norm false "create mapping function"]
  [smoothstep false "maps a value from given interval to [0,1] using hermite mapping"]
  [radians false "convert degrees to radians"]
  [degrees false "convert radians to degrees"]])

^{::clerk/visibility :hide}
(clerk/example
 (m/constrain 3 -5 5)
 (m/constrain -10 -5 5)
 (m/constrain 10 -5 5)
 (m/norm 1.0 -10.0 10.0)
 (m/norm -100.0 -10.0 10.0)
 (m/norm 1.0 -10.0 10.0 100.0 1000.0)
 (m/norm 1.0 10.0 -10.0 100.0 1000.0)
 (m/mnorm 1.0 -10.0 10.0)
 (m/mnorm 1.0 -10.0 10.0 100.0 1000.0)
 (m/cnorm -100.0 -10.0 10.0 100.0 1000.0)
 (m/smoothstep -10.0 10.0 1.0)
 (m/radians 180.0)
 (m/degrees m/PI))

;; Let's create mapping function $[-10,10]$ -> $[100.0 1000.0]$

(def example-norm (m/make-norm -10.0 10.0 100.0 1000.0))

(map example-norm (range -20 20 3))

;; ## Slicing and sampling

^{::clerk/visibility :hide}
(u/table
 [[slice-range false "evenly distrubuted points from given interval (or [0,1] by default)"]
  [sample false "applies function on evenly distributed points, can return [x,y] pairs"]
  [cut false "cuts given range or data into even intervals"]
  [co-intervals false "returns overlapped intervals containing similar number of values"]
  [group-by-intervals false "create a map with intervals and seq of values"]])

;; Notes:

;; * `cut` - left endpoint of the first interval is slightly less than declared
;; * `co-intervals` - can produce less intervals than required, last (optional) argument controls overlap and defaults to 0.5.
;; * `group-by-intervals` - by default uses `co-intervals` to create initial intervals

(def data (repeatedly 200 (fn [] (m/sq (rand-int 10)))))

^{::clerk/visibility :hide}
(clerk/example
 (m/slice-range 5)
 (m/slice-range 10.0 100.0 6)
 (m/sample m/sec 5)
 (m/sample m/sec m/-PI m/PI 5)
 (m/sample m/sec 5 true)
 (m/sample m/sec m/-PI m/PI 5 true)
 (m/cut 10.0 100.0 6)
 (m/cut data 6)
 (m/co-intervals data 6)
 (m/co-intervals data 6 0.7)
 (m/group-by-intervals data)
 (m/group-by-intervals [[0 50] [50 100]] data))

;; ## Rank and order

^{::clerk/visibility :hide}
(u/table
 [[rank false "calculate rank, index is 0 based"]
  [rank1 false "calculate rank, index is 1 based"]
  [order false "ordering indexes, 0 based"]])

;; Rank can be calculated with different strategies for ties: `:average` (default), `:first`, `:last`, `:random`, `:min`, `:max` and `:dense`.

;; Both functions can work on descending order.

^{::clerk/visibility :hide}
(clerk/example
 (m/rank [1 2 2 2 7 7 1])
 (m/rank1 [1 2 2 2 7 7 1])
 (m/rank [1 2 2 2 7 7 1] :average)
 (m/rank [1 2 2 2 7 7 1] :average true)
 (m/rank [1 2 2 2 7 7 1] :first)
 (m/rank [1 2 2 2 7 7 1] :last)
 (m/rank [1 2 2 2 7 7 1] :random)
 (m/rank [1 2 2 2 7 7 1] :min)
 (m/rank [1 2 2 2 7 7 1] :max)
 (m/rank [1 2 2 2 7 7 1] :dense)
 (m/order [1 2 2 2 7 7 1])
 (m/order [1 2 2 2 7 7 1] true))

;; ## Double manipulations

^{::clerk/visibility :hide}
(u/table
 [[next-double false "next (or next nth) possible double value"]
  [prev-double false "previous (or previous nth) possible double value"]
  [double-high-bits false "high 32 bits of binary representation of double value"]
  [double-low-bits false "low 32 bits of binary representation of double value"]
  [double-bits false "64 bits (long) of binary representation of double value"]
  [bits->double false "convert 64 bits to double value"]])

^{::clerk/visibility :hide}
(clerk/example
 (m/next-double 1.0)
 (m/next-double 12345.999)
 (m/next-double 1.0 100)
 (m/prev-double 1.0)
 (m/prev-double 12345.999)
 (m/prev-double 1.0 100)
 (m/double-high-bits 1.1)
 (m/double-low-bits 1.1)
 (m/double-bits 1.1)
 (m/double-bits ##NaN)
 (m/bits->double 1234567890123456789)
 (m/bits->double (m/double-bits -0.991))
 (m/bits->double (m/double-bits ##NaN)))

;; ## Combinatorics

;; All functions (but `factorial20`) are able to work with big numbers as they use and return double.

^{::clerk/visibility :hide}
(u/table
 [[factorial20 false "factorial up to 20!, returns long"]
  [factorial false "factorial, accepts long and returns double"]
  [log-factorial false "log of factorial, accepts long and returns double"]
  [combinations false "n choose k, returns double"]
  [log-combinations false "log of binomial coefficient"]])

^{::clerk/visibility :hide}
(clerk/example
 (map m/factorial20 (range 21))
 (m/factorial 100)
 (m/log-factorial 9)
 (m/exp (m/log-factorial 9))
 (m/combinations 10 5)
 (m/combinations 1000 456)
 (m/log-combinations 1000 456)
 (m/exp (m/log-combinations 10 5)))

;; ## Seq <-> double array

;; Functions which help convert seq or seq of seqs into double-(double)-array and vice versa

^{::clerk/visibility :hide}
(u/table
 [[seq->double-array false "Converts sequence into double array"]
  [seq->double-double-array false "Converts seq of seqs into array of double arrays"]
  [double-array->seq false "Converts double array into sequence"]
  [double-double-array->seq false "Converts array of double arrays into seq of seqs"]])

^{::clerk/visibility :hide}
(clerk/example
 (m/seq->double-array [1 2 3 1/2])
 (m/double-array->seq (m/seq->double-array [1 2 3 1/2]))
 (m/seq->double-double-array [[1 2] [3 1/2]])
 (m/double-double-array->seq (m/seq->double-double-array [[1 2] [3 1/2]]))
 (= m/double-array-type (type (m/seq->double-array [1 2 3 1/2])))
 (= m/double-double-array-type (type (m/seq->double-double-array [1 2 3 1/2]))))

;; ## Constants

^{::clerk/visibility {:code :hide :result :hide}}
(def formulas '{-E "-e" E "e" M_E "e" EPSILON "\\epsilon"
              -HALF_PI "-\\frac{\\pi}{2}" -QUARTER_PI "-\\frac{\\pi}{4}" -THIRD_PI "-\\frac{\\pi}{3}"
              -PI "-\\pi" -TAU "-2\\pi" -TWO_PI "-2\\pi"
              HALF_PI "\\frac{\\pi}{2}" QUARTER_PI "\\frac{\\pi}{4}" THIRD_PI "\\frac{\\pi}{3}"
              PI "\\pi" TAU "2\\pi" TWO_PI "2\\pi"
              M_PI_2 "\\frac{\\pi}{2}" M_PI_4 "\\frac{\\pi}{4}"
              M_PI "\\pi" M_TWOPI "2\\pi"
              FOUR_INV_PI "\\frac{4}{\\pi}"
              GAMMA "\\text{Euler–Mascheroni constant}"
              CATALAN_G "\\text{Catalan G constant}"
              INV_LN2 "\\frac{1}{\\ln(2)}" INV_LOG_HALF "\\frac{1}{\\ln(0.5)}"
              INV_PI "\\frac{1}{\\pi}" INV_SQRT2PI "\\frac{1}{\\sqrt{2\\pi}}"
              INV_TWO_PI "\\frac{1}{2\\pi}" INV_SQRT_2 "\\frac{1}{\\sqrt2}"
              LANCZOS_G "\\text{g constant for Lanczos Gamma approx.}"
              LN10 "\\ln(10)" LN2 "\\ln(2)" LN2_2 "\\frac{\\ln(2)}{2}"
              LOG10E "\\log_{10}(e)" LOG2E "\\log_{2}(e)"
              LOG_HALF "\\ln(0.5)" LOG_PI "\\ln(\\pi)" LOG_TWO_PI "\\ln(2\\pi)"
              MACHINE-EPSILON "\\varepsilon>0, 1+\\varepsilon=1"
              M_1_PI "\\frac{1}{\\pi}" M_2_PI "\\frac{2}{\\pi}"
              M_2_SQRTPI "\\frac{2}{\\sqrt\\pi}"
              M_3PI_4 "\\frac{3\\pi}{4}"
              M_INVLN2 "\\frac{1}{\\ln(2)}" M_IVLN10 "\\frac{1}{\\ln(10)}"
              M_LN2 "\\ln(2)" M_LN10 "\\ln(10)"
              M_LOG10E "\\log_{10}(e)" M_LOG2E "\\log_{2}(e)" M_LOG2_E "\\ln(2)"
              M_SQRT1_2 "\\sqrt{0.5}" M_SQRT2 "\\sqrt2" M_SQRT3 "\\sqrt3" M_SQRT_PI "\\sqrt\\pi"
              PHI "\\phi=\\frac{1+\\sqrt5}{2}" SIXTH "\\frac{1}{6}" THIRD "\\frac{1}{6}"
              SILVER "\\delta_S=1+\\sqrt2" TWO_THIRD "\\frac{2}{3}"
              SQRT5 "\\sqrt5"
              SQRT2 "\\sqrt2" SQRT2_2 "\\frac{\\sqrt2}{2}"
              SQRT3 "\\sqrt3" SQRT3_2 "\\frac{\\sqrt3}{2}"
              SQRT3_3 "\\frac{\\sqrt3}{3}" SQRT3_4 "\\frac{\\sqrt3}{4}"
              SQRTPI "\\sqrt{\\pi}" SQRT2PI "\\sqrt{2\\pi}" SQRT_HALFPI "\\sqrt{\\frac{\\pi}{2}}"
              TWO_INV_PI "\\frac{2}{\\pi}"
              deg-in-rad "\\frac{\\pi}{180}" rad-in-deg "\\frac{180}{\\pi}"
              double-one-minus-epsilon "\\inf\\{x:x<1\\}"})

^{::clerk/visibility :hide ::clerk/viewer u/unpaginated-table}
{:head ["symbol" "formula" "value"]
 :rows (for [nm (->> (ns-publics 'fastmath.core)
                     (map (comp meta second))
                     (filter :const)
                     (map :name)
                     (sort))]
         [nm
          (clerk/md (str "$" (get formulas nm "-") "$"))
          (var-get (ns-resolve 'fastmath.core nm))])}

;; ## List of symbols

^{::clerk/visibility :hide}
(u/make-public-fns-table 'fastmath.core)
