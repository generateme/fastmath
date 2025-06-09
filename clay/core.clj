^:kindly/hide-code
(ns core
  (:require [fastmath.dev.codox :as codox]
            [fastmath.dev.ggplot :as gg]
            [fastmath.dev.clay :as utls]
            [fastmath.dev.helpers :as hp]
            [fastmath.core :as m]
            [scicloj.kindly.v4.kind :as kind]
            [fastmath.fields.n :as n]))

;; # Core {.unnumbered}

;; Collection of type hinted math macros and functions. Partially backed by Java static functions and exposed as macros. They are prepared to accept primitive long or double arguments and return long or double only.

;; There is a possibility to replace `clojure.core` functions with a selection of `fastmath.core` macros. Call:

;; * `(m/use-primitive-operators)` to replace functions with macros
;; * `(m/unuse-primitive-operators)` to revert replacement.

;; Be aware that there are some differences and `fastmath.core` versions shoudn't be treated as a drop-in replacement for `clojure.core` versions. Also, since Clojure 1.12, always call `unuse-primitive-operators` at the end of the namespace.

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

;; ## Basic operations

;; Basic math operations.

;; When used in an expression oprations are inlined and can accept mixture of `long` and `double` values. If all values are of `long` primitive type, `long` is returned, `double` otherwise.

;; When used in higher order function, `double` is returned always. To operate on `long` primitive type, reach for `long-` versions.

;; ::: {.callout-tip title="Defined functions"}
;; * `+`, `-`, `*`, `/`, `quot`
;; * `inc`, `dec`
;; * `min`, `max`, `smooth-max`, `constrain`
;; * `rem`, `mod`, `remainder`, `wrap`
;; * `abs`
;;
;; `long` versions
;;
;; * `long-add`, `long-sub`, `long-mult`, `long-div`, `long-quot`
;; * `long-inc`, `long-dec`
;; * `long-min`, `long-max`
;; * `long-rem`, `long-mod`
;; * `long-abs`
;; :::

;; ### Arithmetics

;; * addition, incrementation
;; * subtraction, decrementation
;; * multiplication
;; * division
;; * absolute value

;; ::: {.callout-warning title="Division differences"}
;; Please note that there some differences between division in `fastmath` and `clojure.core`
;;
;; * when called with one argument (`double` or `long`) `m//` always returns reciprocal (`clojure.core//` returns a ratio)
;; * when called on `long` arguments, `m//` is a long division (`clojure.core//` returns a ratio)
;; * `m//` for two `long` arguments is equivalent to `m/quot`
;; :::

(utls/callout "note" "Examples"
  (kind/md "Addition")
  (utls/examples
    (m/+)
    (m/+ 1 2 3 4)
    (m/+ 1.0 2.5 3 4)
    (reduce m/+ [1 2 3]))
  (utls/examples
    (m/long-add)
    (m/long-add 1 2 3 4)
    (m/long-add 1.0 2.5 3 4)
    (reduce m/long-add [1 2 3.5]))
  (kind/md "Subtraction")
  (utls/examples
    [(m/- 1) (m/- 1.0)]
    (m/- 1 2 3 4)
    (m/- 1.0 2.5 3 4)
    (reduce m/- [1 2 3]))
  (utls/examples
    (m/long-sub 1)
    (m/long-sub 1 2 3 4)
    (m/long-sub 1.0 2.5 3 4)
    (reduce m/long-sub [1 2 3.5]))
  (kind/md "Multiplication")
  (utls/examples
    (m/*)
    (m/* 1 2 3 4)
    (m/* 1.0 2.5 3 4)
    (reduce m/* [1 2 3]))
  (utls/examples
    (m/long-mult)
    (m/long-mult 1 2 3 4)
    (m/long-mult 1.0 2.5 3 4)
    (reduce m/long-mult [1 2 3.5]))
  (kind/md "Division")
  (utls/examples
    [(m// 2) (m// 2) (/ 2)]
    (m// 1 2 3 4)
    (m// 1.0 2.5 3 4)
    (reduce m// [1 2 3])
    (m/quot 10.5 -3))
  (utls/examples
    (m/long-div 2)
    (m/long-div 100 5 3)
    (m/long-div 100.5 2.5 3)
    (reduce m/long-div [100 2 3.5])
    (m/long-quot 10 -3))
  (kind/md "Increment and decrement")
  (utls/examples
    (m/inc 4)
    (m/inc 4.5)
    (m/dec 4)
    (m/dec 4.5)
    (map m/inc [1 2 3.5 4.5]))
  (utls/examples
    (m/long-inc 4)
    (m/long-inc 4.5)
    (m/long-dec 4)
    (m/long-dec 4.5)
    (map m/long-inc [1 2 3.5 4.5]))
  (kind/md "Absolute value")
  (utls/examples
    (m/abs -3)
    (m/long-abs -3)
    (m/abs -3.5)
    (m/long-abs -3.5)))

;; ### Remainders

;; * `rem` and `mod` are the same as in `clojure.core`,
;; * `remainder` returns $dividend - divisor * n$, where $n$ is the mathematical integer closest to $\frac{dividend}{divisor}$. Returned value is inside the $[\frac{-|divisor|}{2},\frac{|divisor|}{2}]$ range.
;; * `wrap` wraps the value to be within given interval (right open) $[a,b)$ `

^:kind/table
[[(gg/->image (gg/function #(m/mod % 2) {:x [-5 5] :ylim [-2 2] :title "mod(x,2)"}))
  (gg/->image (gg/function #(m/rem % 2) {:x [-5 5] :ylim [-2 2] :title "rem(x,2)"}))]
 [(gg/->image (gg/function #(m/remainder % 2) {:x [-5 5] :ylim [-2 2] :title "remainder(x,2)"}))
  (gg/->image (gg/function #(m/wrap -1 1 %) {:x [-5 5] :ylim [-2 2] :title "wrap(-1,1,x)"}))]]

(utls/examples-note
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
  (m/wrap -1.25 1.25 -1.25)
  (m/wrap -1.25 1.25 1.25)
  (m/wrap [-1.25 1.25] -1.35))

;; ### Min, max, constrain

;; Constrain is a macro which is equivalent to `(max (min value mx) mn)`

(utls/examples-note
  (m/min 1 2 -3)
  (m/min 1.0 2 -3)
  (m/max 1 2 -3)
  (m/max 1.0 2 -3)
  (m/constrain 10 -1 1)
  (m/constrain -10 -1 1)
  (m/constrain 0 -1 1))

;; #### Smooth maximum

;; Smooth maximum is a family of functions $\max_\alpha(xs)$ for which $\lim_{\alpha\to\infty}\max_\alpha(xs)=\max(xs)$. 

;; Five types of smooth maximum are defined (see [wikipedia](https://en.wikipedia.org/wiki/Smooth_maximum) for formulas):

;; * `:lse` - LogSumExp (default)
;; * `:boltzmann` - Boltzmann operator, works for small alpha values
;; * `:mellowmax`
;; * `:p-norm`
;; * `:smu` - smooth maximum unit, $\epsilon=\frac{1}{\alpha}$

;; `:lse`, `:boltzmann` and `:mellowmax` are also smooth minimum for negative $\alpha$ values.

;; The following plots show value of the smooth max for different $\alpha$ and set of the numbers equal to `[-3.5 -2 -1 0.1 3 4]`. Blue dashed horizontal lines are minimum (-3.5) and maximum values (4.0).

^:kind/table
[[(-> (gg/function #(m/smooth-max [-3.5 -2 -1 0.1 3 4] % :lse) {:x [-3.0 3.0]
                                                                :xlab "alpha"
                                                                :ylim [-15 20]
                                                                :title ":lse"})
      (gg/geom-hline :yintercept 4 :color gg/color-light :linetype 2)
      (gg/geom-hline :yintercept -3.5 :color gg/color-light :linetype 2)
      (gg/->image))

  (-> (gg/function #(m/smooth-max [-3.5 -2 -1 0.1 3 4] % :boltzmann) {:x [-10.0 10.0]
                                                                      :xlab "alpha"
                                                                      :title ":boltzmann"})
      (gg/geom-hline :yintercept 4 :color gg/color-light :linetype 2)
      (gg/geom-hline :yintercept -3.5 :color gg/color-light :linetype 2)
      (gg/->image))]

 [(-> (gg/function #(m/smooth-max [-3.5 -2 -1 0.1 3 4] % :mellowmax) {:x [-10.0 10.0]
                                                                      :xlab "alpha"
                                                                      :title ":mellowmax"})
      (gg/geom-hline :yintercept 4 :color gg/color-light :linetype 2)
      (gg/geom-hline :yintercept -3.5 :color gg/color-light :linetype 2)
      (gg/->image))]]

;; The following plots are defined only for positive $\alpha$.

^:kind/table
[[(-> (gg/function #(m/smooth-max [-3.5 -2 -1 0.1 3 4] % :p-norm) {:x [0.1 10.0]
                                                                   :ylim [nil 20]
                                                                   :xlab "alpha"
                                                                   :title ":p-norm"})
      (gg/geom-hline :yintercept 4 :color gg/color-light :linetype 2)
      (gg/->image))

  (-> (gg/function #(m/smooth-max [-3.5 -2 -1 0.1 3 4] % :smu) {:x [0.1 10.0]
                                                                :xlab "alpha"
                                                                :title ":smu"})
      (gg/geom-hline :yintercept 4 :color gg/color-light :linetype 2)
      (gg/->image))]]

(utls/examples-note
  (m/smooth-max [-3.5 -2 -1 0.1 3 4] 4.0 :lse)
  (m/smooth-max [-3.5 -2 -1 0.1 3 4] -4.0 :lse)
  (m/smooth-max [-3.5 -2 -1 0.1 3 4] 4.0 :boltzmann)
  (m/smooth-max [-3.5 -2 -1 0.1 3 4] -4.0 :boltzmann)
  (m/smooth-max [-3.5 -2 -1 0.1 3 4] 4.0 :mellowmax)
  (m/smooth-max [-3.5 -2 -1 0.1 3 4] -4.0 :mellowmax)
  (m/smooth-max [-3.5 -2 -1 0.1 3 4] 4.0 :p-norm)
  (m/smooth-max [-3.5 -2 -1 0.1 3 4] 4.0 :smu))

;; ### fma

;; Fused multiply-add $fma(a,b,c)=a+bc$ is the operation implemented with better accuracy in Java 9+ and as one instruction (see more [here](https://en.wikipedia.org/wiki/Multiply%E2%80%93accumulate_operation#Fused_multiply%E2%80%93add) and [here](https://docs.oracle.com/javase/9/docs/api/java/lang/Math.html#fma-double-double-double-)). When Java 8 is used `fma` is replaced with direct `a+bc` formula.

;; ::: {.callout-tip title="Defined functions"}
;; * `fma`, `muladd`, `negmuladd`
;; * `difference-of-products`, `sum-of-products`
;; :::

;; $$\operatorname{fma}(a,b,c)=\operatorname{muladd}(a,b,c)=a+bc$$
;; $$\operatorname{negmuladd}(a,b,c)=\operatorname{fma}(-a,b,c)$$

;; `difference-of-products` (dop) and `sum-of-products` (sop) are using Kahan's algorithm to avoid [catastrophic cancellation](https://en.wikipedia.org/wiki/Catastrophic_cancellation).

;; $$\operatorname{dop}(a,b,c,d)=ab-cd=\operatorname{fma}(a,b,-cd)+\operatorname{fma}(-c,d,cd)$$
;; $$\operatorname{sop}(a,b,c,d)=ab+cd=\operatorname{fma}(a,b,cd)+\operatorname{fma}(c,d,-cd)$$

;; The following example shows that $x^2-y^2$ differs from the best floating point approximation which is equal `1.8626451518330422e-9`.

(let [x (m/inc (m/pow 2 -29))
      y (m/inc (m/pow 2 -30))]
  {:proper-value (m/difference-of-products x x y y)
   :wrong-value (m/- (m/* x x) (m/* y y))})

(utls/examples-note
  (m/fma 3 4 5)
  (m/muladd 3 4 5)
  (m/negmuladd 3 4 5)
  (m/difference-of-products 3 3 4 4)
  (m/sum-of-products 3 3 4 4))

;; ## Rounding

;; Various rounding functions.

;; ::: {.callout-tip title="Defined functions"}
;; * `floor`, `ceil`
;; * `round`, `round-even`, `rint`, `approx`, `trunc`, `itrunc`
;; * `qfloor`, `qceil`, `qround`
;; * `frac`, `sfrac`
;; * `round-up-pow2`
;; :::

;; * `floor`, `ceil` and `rint` accept additional argument, `scale`, which allows to round to the nearest multiple of scale.
;; * `round` returns `long` while `rint` returns `double`
;; * `round-even` performs IEEE / IEC rounding (even-odd or bankers' rounding)
;; * `approx` rounds number to the given number of digits, uses `bigdec`
;; * `trunc` returns integer part of a number, `frac` returns fractional part
;; * `trunc` returns `double` while `itrunc` returns long
;; * `sfrac` keeps sign of the argument
;; * `qfloor`, `qceil` and `qround` are implemented using casting to `long`
;; * `round-up-pow2` rounds to the lowest power of 2 greater than an argument, $2^{\left\lceil{\log_2{x}}\right\rceil}$, returns `long`.

^:kind/table
[[(gg/->image (gg/function m/floor {:x [-2 2] :title "floor(x)"}))
  (gg/->image (gg/function m/ceil {:x [-2 2] :title "ceil(x)"}))]
 [(gg/->image (gg/function m/round {:x [-2 2] :title "round(x)"}))
  (gg/->image (gg/function m/trunc {:x [-2 2] :title "trunc(x)"}))]
 [(gg/->image (gg/function m/frac {:x [-2 2] :title "frac(x)"}))
  (gg/->image (gg/function m/sfrac {:x [-2 2] :title "sfrac(x)"}))]]

(utls/callout "note" "Examples"
  (utls/examples
    (map m/floor [-10.5 10.5])
    (m/floor 10.5 4.0)
    (map m/ceil [-10.5 10.5])
    (m/ceil 10.5 4.0)
    (map m/rint [-10.51 -10.5 -10.49 10.49 10.5 10.51])
    (m/rint 10.5 4.0)
    (m/rint 10.591 0.1)
    (map m/round [-10.51 -10.5 -10.49 10.49 10.5 10.51])
    (map m/round-even [-10.51 -10.5 -10.49 10.49 10.5 10.51])
    (map m/qfloor [-10.5 10.5])
    (map m/qceil [-10.5 10.5])
    (map m/qround [-10.51 -10.5 -10.49 10.49 10.5 10.51])
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
    (map m/round-up-pow2 (range 10)))
  
  (kind/md "The difference between `rint` and `round`. `round` is bounded by minimum and maximum `long` values.")
  (utls/examples
    (m/rint 1.234567890e30)
    (m/round 1.234567890e30)))

;; ## Sign

;; Sign of the number.

;; ::: {.callout-tip title="Defined functions"}
;; * `signum` and `sgn`
;; * `copy-sign`
;; :::

;; $$\operatorname{signum}(x)=\begin{cases}
;; -1 & x<0 \\
;; 1 & x>0 \\
;; 0 & x=0
;; \end{cases}$$

;; $$\operatorname{sgn}(x)=\begin{cases}
;; -1 & x<0 \\
;; 1 & x\geq 0
;; \end{cases}$$

;; `copy-sign` sets the sign of the second argument to the first. Please note that `-0.0` is negative and `0.0` is positive.

;; $$\operatorname{copy-sign}(x,y)=\begin{cases}
;; |x| & y>0 \lor y=0.0\\
;; -|x| & y<0 \lor y=-0.0
;; \end{cases}$$

^:kind/table
[[(gg/->image (gg/function m/signum {:x [-0.1 0.1] :title "signum(x)"}))
  (gg/->image (gg/function m/sgn {:x [-0.1 0.1] :title "sgn(x)"}))]]

(utls/examples-note
  (m/signum -2.5)
  (m/signum 2.5)
  (m/sgn -2.5)
  (m/sgn 2.5)
  (m/signum 0)
  (m/sgn 0)
  (m/copy-sign 123 -10)
  (m/copy-sign -123 10)
  (m/copy-sign 123 -0.0)
  (m/copy-sign -123 0.0))

;; ## Comparison and Predicates

;; Various predicates and comparison functions

;; ::: {.callout-tip title="Defined functions"}
;; * `==`, `eq`, `not==`, `<`, `>`, `<=`, `>=`
;; * `approx-eq`, `approx=`, `delta-eq`, `delta=`
;; * `zero?`, `negative-zero?`, `near-zero?`, `one?`
;; * `neg?`, `pos?`, `not-neg?`, `not-pos?`
;; * `even?`, `odd?`
;; * `integer?`
;; * `nan?`, `inf?`, `pos-inf?`, `neg-inf?`, `invalid-double?`, `valid-double?`
;; * `between?`, `between-?`
;; :::

;; ### Comparison

;; Comparison functions operate on primitive values and can handle multiple arguments, chaining the comparison (e.g., `(m/< 1 2 3)` is equivalent to `(and (m/< 1 2) (m/< 2 3))`).

;; Standard comparisons:

;; * `==`, `eq`: Primitive equality. Note that `(m/== 1.0 1)` is true. Multi-arity checks if the first argument is equal to all subsequent arguments.
;; * `not==`: Primitive inequality. Multi-arity checks if all arguments are pairwise unique.
;; * `<`, `>`, `<=`, `>=`: Standard primitive inequalities. Multi-arity checks if the values are monotonically increasing/decreasing.

;; Approximate equality:

;; * `approx-eq`, `approx=`: Checks for equality after rounding to a specified number of decimal digits. Can be inaccurate.
;; * `delta-eq`, `delta=`: Checks if the absolute difference between two numbers is within a given tolerance (absolute and/or relative). This is the recommended way to compare floating-point numbers for near equality
;;   * With absolute tolerance: $|a - b| < \text{abs-tol}$
;;   * With absolute and relative tolerance: $|a - b| < \max(\text{abs-tol}, \text{rel-tol} \cdot \max(|a|, |b|))$


;; Range checks:

;; * `between?`: Checks if a value is within a closed interval $[a, b]$, i.e., $a \le \text{value} \le b$.
;; * `between-?`: Checks if a value is within a half-open interval $(a, b]$, i.e., $a < \text{value} \le b$.


(utls/examples-note
  (m/== 1.0 1)
  (m/== 1.0 2)
  (m/== 1 1 2 3 4)
  (m/eq 1.0 1)
  (m/not== 1.0 1)
  (m/not== 1.0 2)
  (m/not== 1 2 3 4)
  (m/not== 1 2 3 1 4)
  (m/< 1.0 1)
  (m/< 1.0 2)
  (m/< 2 1.0)
  (m/< 1 2 3 4 10)
  (m/< 10 4 3 2 1)
  (m/<= 1.0 1)
  (m/<= 1.0 2)
  (m/<= 2 1.0)
  (m/> 1.0 1)
  (m/> 1.0 2)
  (m/> 2 1.0)
  (m/> 1 2 3 4 10)
  (m/> 10 4 3 2 1)
  (m/>= 1.0 1)
  (m/>= 1.0 2)
  (m/>= 2 1.0)
  (m/between? -1 1 -1)
  (m/between? -1 1 0)
  (m/between? -1 1 1)
  (m/between-? -1 1 -1)
  (m/between-? -1 1 0)
  (m/between-? -1 1 1))

(utls/examples-note
  (m/approx-eq 10 10.01)
  (m/approx-eq 10 10.001)
  (m/approx-eq 10 10.01 1)
  (m/approx-eq 10 10.001 5)
  (m/delta-eq 10 10.01)
  (m/delta-eq 10 10.01 0.1)
  (m/delta-eq 1.0e-6 1.01e-6)
  (m/delta-eq 1.0e-6 1.01e-6 1.0e-8)
  (m/delta-eq 1.0e-6 1.01e-6 1.0e-8 1.0e-2))

;; ### Predicates

;; Basic predicates for common number properties:

;; * `zero?`: Checks if value is 0 or 0.0.
;; * `negative-zero?`: Checks specifically for the floating point `-0.0`.
;; * `near-zero?`: Checks if absolute value is within absolute and/or relative tolerance of zero: $|x| < \text{abs-tol}$ or $|x| < \max(\text{abs-tol}, \text{rel-tol} \cdot |x|)$.
;; * `one?`: Checks if value is 1 or 1.0.
;; * `neg?`: Checks if value is $< 0$.
;; * `pos?`: Checks if value is $> 0$.
;; * `not-neg?`: Checks if value is $\ge 0$.
;; * `not-pos?`: Checks if value is $\le 0$.
;; * `even?`: Checks if a long is even.
;; * `odd?`: Checks if a long is odd.
;; * `integer?`: Checks if a number (long or double) has a zero fractional part.

(utls/examples-note
  (m/zero? 0)
  (m/zero? 0.0)
  (m/zero? -0.0)
  (m/zero? 1)
  (m/negative-zero? 0.0)
  (m/negative-zero? -0.0)
  (m/one? 1)
  (m/one? 1.0)
  (m/neg? -1)
  (m/neg? 0)
  (m/neg? 1.0)
  (m/pos? -1)
  (m/pos? 0)
  (m/pos? 1.0)
  (m/not-neg? -1)
  (m/not-neg? 0)
  (m/not-neg? 1.0)
  (m/not-pos? -1)
  (m/not-pos? 0)
  (m/not-pos? 1.0)
  (m/even? 0)
  (m/even? 1)
  (m/even? 2)
  (m/odd? 0)
  (m/odd? 1)
  (m/odd? 2)
  (m/integer? 1)
  (m/integer? 1.0)
  (m/integer? 1.1))

;; Predicates for floating point special values:

;; * `nan?`: Checks if value is Not-a-Number (NaN).
;; * `inf?`: Checks if value is positive or negative infinity (Inf or -Inf).
;; * `pos-inf?`: Checks if value is positive infinity (Inf).
;; * `neg-inf?`: Checks if value is negative infinity (-Inf).
;; * `invalid-double?`: Checks if value is not a finite double (NaN or ±Inf).
;; * `valid-double?`: Checks if value is a finite double (not NaN or ±Inf).

(utls/examples-note
  (m/nan? ##NaN)
  (m/nan? ##Inf)
  (m/nan? ##-Inf)
  (m/nan? 1)
  (m/inf? ##NaN)
  (m/inf? ##Inf)
  (m/inf? ##-Inf)
  (m/inf? 1)
  (m/pos-inf? ##NaN)
  (m/pos-inf? ##Inf)
  (m/pos-inf? ##-Inf)
  (m/pos-inf? 1)
  (m/neg-inf? ##NaN)
  (m/neg-inf? ##Inf)
  (m/neg-inf? ##-Inf)
  (m/neg-inf? 1)
  (m/valid-double? ##NaN)
  (m/valid-double? ##Inf)
  (m/valid-double? ##-Inf)
  (m/valid-double? 1)
  (m/invalid-double? ##NaN)
  (m/invalid-double? ##Inf)
  (m/invalid-double? ##-Inf)
  (m/invalid-double? 1))

;; ## Trigonometry

;; Trigonometric (with historical variants) and hyperbolic functions.

;; ::: {.callout-tip title="Defined functions"}
;; * `radians`, `degrees`
;; * `sin`, `cos`, `tan`, `cot`, `sec`, `csc`
;; * `qsin`, `qcos`
;; * `sinpi`, `cospi`, `tanpi`, `cotpi`, `secpi`, `cscpi`
;; * `asin`, `acos`, `atan`, `atan2`, `acot`, `asec`, `acsc`
;; * `sinh`, `cosh`, `tanh`, `coth`, `sech`, `scsh`
;; * `asinh`, `acosh`, `atanh`, `acoth`, `asech`, `ascsh`
;; * `crd`, `acrd`
;; * `versin`, `coversin`, `vercos`, `covercos`
;; * `aversin`, `acoversin`, `avercos`, `acovercos`
;; * `haversin`, `hacoversin`, `havercos`, `hacovercos`
;; * `ahaversin`, `ahacoversin`, `ahavercos`, `ahacovercos`
;; * `exsec`, `excsc`
;; * `aexsec`, `aexcsc`
;; * `sinc`
;; :::

;; ### Angle conversion

;; Convert between radians and degrees

;; * `radians(deg)`: Converts an angle from degrees to radians.
;; $$\text{radians} = \text{degrees} \cdot \frac{\pi}{180^\circ}$$
;; * `degrees(rad)`: Converts an angle from radians to degrees.
;; $$\text{degrees} = \text{radians} \cdot \frac{180^\circ}{\pi}$$

(utls/examples-note
  (m/radians 180)
  (m/degrees m/PI)
  (m/radians 90)
  (m/degrees m/HALF_PI))

;; ### Trigonometric

;; Standard trigonometric functions:

;; * `sin(x)`: Sine of `x`.
;; * `cos(x)`: Cosine of `x`.
;; * `tan(x)`: Tangent of `x`, $\tan(x) = \frac{\sin(x)}{\cos(x)}$.
;; * `cot(x)`: Cotangent of `x`, $\cot(x) = \frac{1}{\tan(x)}$.
;; * `sec(x)`: Secant of `x`, $\sec(x) = \frac{1}{\cos(x)}$.
;; * `csc(x)`: Cosecant of `x`, $\csc(x) = \frac{1}{\sin(x)}$.

(kind/table
 [[(gg/->image (gg/functions [["sin" m/sin] ["cos" m/cos]] {:x [-4 4] :title "sin(x), cos(x)"}))
   (gg/->image (gg/functions [["tan" m/tan] ["cot" m/cot]] {:x [-4 4] :ylim [-3 3] :title "tan(x), cot(x)"}))]
  [(gg/->image (gg/functions [["sec" m/sec] ["csc" m/csc]] {:x [-4 4] :ylim [-3 3] :title "sec(x), csc(x)"}))]])

(utls/examples-note
  (m/sin m/HALF_PI)
  (m/cos m/PI)
  (m/tan m/QUARTER_PI)
  (m/cot m/HALF_PI)
  (m/sec 0.0)
  (m/csc m/HALF_PI))

;; Quick trigonometric functions:

;; * `qsin(x)`: Fast, less accurate sine.
;; * `qcos(x)`: Fast, less accurate cosine.

^:kind/table
[[(gg/->image (gg/function m/qsin {:x [-4 4] :title "qsin(x)"}))
  (gg/->image (gg/function (hp/relative-error m/sin m/qsin) {:x [-4 4] :title "Relative error between sin(x) and qsin(x)."}))
  ]
 [(gg/->image (gg/function m/qcos {:x [-4 4] :title "qcos(x)"}))
  (gg/->image (gg/function (hp/relative-error m/cos m/qcos) {:x [-4 4] :title "Relative error between cos(x) and qcos(x)."}))]]

(utls/examples-note
  (m/qsin 0.1)
  (m/sin 0.1)
  (m/qcos 0.1)
  (m/cos 0.1))

;; Pi-scaled trigonometric functions:

;; * `sinpi(x)`: Sine of $\pi x$, $\sin(\pi x)$.
;; * `cospi(x)`: Cosine of $\pi x$, $\cos(\pi x)$.
;; * `tanpi(x)`: Tangent of $\pi x$, $\tan(\pi x)$.
;; * `cotpi(x)`: Cotangent of $\pi x$, $\cot(\pi x)$.
;; * `secpi(x)`: Secant of $\pi x$, $\sec(\pi x)$.
;; * `cscpi(x)`: Cosecant of $\pi x$, $\csc(\pi x)$.

(utls/examples-note
  (m/sinpi 0.5)   
  (m/cospi 1)     
  (m/tanpi 0.25)  
  (m/cotpi 0.5)   
  (m/secpi 0)     
  (m/cscpi 0.5))

;; ### Inverse trigonometric

;; Inverse trigonometric functions:

;; * `asin(x)`: Arcsine of `x`, $\arcsin(x)$.
;; * `acos(x)`: Arccosine of `x`, $\arccos(x)$.
;; * `atan(x)`: Arctangent of `x`, $\arctan(x)$.
;; * `atan2(y, x)`: Arctangent of $\frac{y}{x}$, returning the angle in the correct quadrant.
;; * `acot(x)`: Arccotangent of `x`, $\operatorname{arccot}(x) = \frac{\pi}{2} - \arctan(x)$.
;; * `asec(x)`: Arcsecant of `x`, $\operatorname{arcsec}(x) = \arccos(\frac{1}{x})$.
;; * `acsc(x)`: Arccosecant of `x`, $\operatorname{arccsc}(x) = \arcsin(\frac{1}{x})$.

(kind/table
 [[(gg/->image (gg/functions [["asin" m/asin] ["acos" m/acos]] {:x [-1.1 1.1] :title "asin(x), acos(x)"}))
   (gg/->image (gg/functions [["atan" m/atan] ["acot" m/acot]] {:x [-4 4] :title "atan(x), acot(x)"}))]
  [(gg/->image (gg/functions [["asec" m/asec] ["acsc" m/acsc]] {:x [-4 4] :title "asec(x), acsc(x)"}))]])

(utls/examples-note
  (m/asin 1)
  (m/acos -1)
  (m/atan 1)
  (m/acot 0)
  (m/asec 1)
  (m/acsc 1)
  (m/atan2 1 1))

;; ### Special

;; `sinc(x)`: Sinc function: $\operatorname{sinc}(x) = \frac{\sin(\pi x)}{\pi x}$ for $x \ne 0$, and $1$ for $x=0$.

(utls/examples-note
  (m/sinc 0.0)
  (m/sinc 1.0))

(gg/->image (gg/function m/sinc {:x [-4 4] :title "sinc(x)"}))

;; ### Hyperbolic

;; Hyperbolic functions

;; * `sinh(x)`: Hyperbolic sine, $\sinh(x) = \frac{e^x - e^{-x}}{2}$.
;; * `cosh(x)`: Hyperbolic cosine, $\cosh(x) = \frac{e^x + e^{-x}}{2}$.
;; * `tanh(x)`: Hyperbolic tangent, $\tanh(x) = \frac{\sinh(x)}{\cosh(x)}$.
;; * `coth(x)`: Hyperbolic cotangent, $\coth(x) = \frac{1}{\tanh(x)}$.
;; * `sech(x)`: Hyperbolic secant, $\operatorname{sech}(x) = \frac{1}{\cosh(x)}$.
;; * `csch(x)`: Hyperbolic cosecant, $\operatorname{csch}(x) = \frac{1}{\sinh(x)}$.

(kind/table
 [[(gg/->image (gg/functions [["sinh" m/sinh] ["cosh" m/cosh]] {:x [-4 4] :title "sinh(x), cosh(x)"}))
   (gg/->image (gg/functions [["tanh" m/tanh] ["coth" m/coth]] {:x [-4 4] :ylim [-3 3] :title "sinh(x), cosh(x)"}))]
  [(gg/->image (gg/functions [["sech" m/sech] ["csch" m/csch]] {:x [-4 4] :ylim [-3 3] :title "sech(x), csch(x)"}))]])

(utls/examples
  (m/sinh 0)
  (m/cosh 0)
  (m/tanh 0)
  (m/coth 1)
  (m/sech 0)
  (m/csch 1))

;; ### Inverse hyperbolic

;; Inverse hyperbolic functions:

;; * `asinh(x)`: Area hyperbolic sine, $\operatorname{arsinh}(x) = \ln(x + \sqrt{x^2 + 1})$.
;; * `acosh(x)`: Area hyperbolic cosine, $\operatorname{arcosh}(x) = \ln(x + \sqrt{x^2 - 1})$ for $x \ge 1$.
;; * `atanh(x)`: Area hyperbolic tangent, $\operatorname{artanh}(x) = \frac{1}{2} \ln(\frac{1 + x}{1 - x})$ for $-1 < x < 1$.
;; * `acoth(x)`: Area hyperbolic cotangent, $\operatorname{arcoth}(x) = \operatorname{artanh}(\frac{1}{x})$.
;; * `asech(x)`: Area hyperbolic secant, $\operatorname{arsech}(x) = \operatorname{arcosh}(\frac{1}{x})$.
;; * `acsch(x)`: Area hyperbolic cosecant, $\operatorname{arcsch}(x) = \operatorname{arsinh}(\frac{1}{x})$.

(utls/examples-note
  (m/asinh 0)
  (m/acosh 1)
  (m/atanh 0)
  (m/acoth 2)
  (m/asech 1)
  (m/acsch 1))

(kind/table
 [[(gg/->image (gg/functions [["asinh" m/asinh] ["acosh" m/acosh]] {:x [-4 4] :title "asinh(x), acosh(x)"}))
   (gg/->image (gg/functions [["atanh" m/atanh] ["acoth" m/acoth]] {:x [-4 4] :ylim [-3 3] :title "asinh(x), acosh(x)"}))]
  [(gg/->image (gg/functions [["asech" m/asech] ["acsch" m/acsch]] {:x [-4 4] :ylim [-3 3] :title "asech(x), acsch(x)"}))]])

;; ### Historical

;; Historical/Specific trigonometric functions:

;; * `crd(x)`: Chord, $\operatorname{crd}(x) = 2 \sin(\frac{x}{2})$.
;; * `acrd(x)`: Inverse chord, $\operatorname{acrd}(x) = 2 \arcsin(\frac{x}{2})$.
;; * `versin(x)`: Versine, $\operatorname{versin}(x) = 1 - \cos(x)$.
;; * `coversin(x)`: Coversine, $\operatorname{coversin}(x) = 1 - \sin(x)$.
;; * `vercos(x)`: Vercosine, $\operatorname{vercos}(x) = 1 + \cos(x)$.
;; * `covercos(x)`: Covercosine, $\operatorname{covercos}(x) = 1 + \sin(x)$.
;; * `haversin(x)` / `haversine(x)`: Haversine, $\operatorname{hav}(x) = \sin^2(\frac{x}{2}) = \frac{1 - \cos(x)}{2}$. Also computes the haversine value for pairs of geographic coordinates.
;; * `hacoversin(x)`: Hacoversine, $\operatorname{hacov}(x) = \frac{1 - \sin(x)}{2}$.
;; * `havercos(x)`: Havercosine, $\operatorname{haver}(x) = \frac{1 + \cos(x)}{2}$.
;; * `hacovercos(x)`: Hacovercosine, $\operatorname{hacover}(x) = \frac{1 + \sin(x)}{2}$.
;; * `exsec(x)`: Exsecant, $\operatorname{exsec}(x) = \sec(x) - 1$.
;; * `excsc(x)`: Excosecant, $\operatorname{excsc}(x) = \csc(x) - 1$.
;; * `aversin(x)`: Arc versine, $\operatorname{aversin}(x) = \arccos(1 - x)$.
;; * `acoversin(x)`: Arc coversine, $\operatorname{acoversin}(x) = \arcsin(1 - x)$.
;; * `avercos(x)`: Arc vercosine, $\operatorname{avercos}(x) = \arccos(x - 1)$.
;; * `acovercos(x)`: Arc covercosine, $\operatorname{acovercos}(x) = \arcsin(x - 1)$.
;; * `ahaversin(x)`: Arc haversine, $\operatorname{ahav}(x) = \arccos(1 - 2x)$.
;; * `ahacoversin(x)`: Arc hacoversine, $\operatorname{ahacov}(x) = \arcsin(1 - 2x)$.
;; * `ahavercos(x)`: Arc havercosine, $\operatorname{ahaver}(x) = \arccos(2x - 1)$.
;; * `ahacovercos(x)`: Arc hacovercosine, $\operatorname{ahacover}(x) = \arcsin(2x - 1)$.
;; * `aexsec(x)`: Arc exsecant, $\operatorname{aexsec}(x) = \operatorname{arcsec}(1 + x)$.
;; * `aexcsc(x)`: Arc excosecant, $\operatorname{aexcsc}(x) = \operatorname{arccsc}(1 + x)$.

(utls/examples-note
  (m/crd m/HALF_PI)
  (m/versin m/HALF_PI)
  (m/coversin m/PI)
  (m/vercos m/HALF_PI)
  (m/covercos m/PI)
  (m/haversin m/HALF_PI)
  (m/hacoversin m/PI)
  (m/havercos m/HALF_PI)
  (m/hacovercos m/PI)
  (m/exsec m/PI)
  (m/excsc m/HALF_PI)
  
  (m/acrd 2.0)
  (m/aversin 1.0)
  (m/acoversin 0.0)
  (m/avercos 1.0)
  (m/acovercos 0.0)
  (m/ahaversin 1.0)
  (m/ahacoversin 0.0)
  (m/ahavercos 0.0)
  (m/ahacovercos 1.0)
  (m/aexsec -2.0)
  (m/aexcsc 0.0))

(kind/table
 [[(gg/->image (gg/function m/crd
                            {:x [-7 7] :title "crd(x)"}))
   (gg/->image (gg/functions [["exsec" m/exsec]
                              ["excsc" m/excsc]]
                             {:x [-4 4] :ylim [-5 3]}))]
  [(gg/->image (gg/functions [["versin" m/versin]
                              ["coversin" m/coversin]
                              ["vercos" m/vercos]
                              ["covercos" m/covercos]]
                             {:x [-4 4]}))

   (gg/->image (gg/functions [["haversin" m/haversin]
                              ["hacoversin" m/hacoversin]
                              ["havercos" m/havercos]
                              ["hacovercos" m/hacovercos]]
                             {:x [-4 4] :ylim [0 2]}))]
  [(gg/->image (gg/function m/acrd {:x [-2.1 2.1] :title "acrd(x)"}))
   (gg/->image (gg/functions [["aexsec" m/aexsec]
                              ["aexcsc" m/aexcsc]] {:x [-0.1 2] :ylim [0 m/HALF_PI]}))]
  [(gg/->image (gg/functions [["aversin" m/aversin]
                              ["acoversin" m/acoversin]
                              ["avercos" m/avercos]
                              ["acovercos" m/acovercos]] {:x [-0.1 2]}))
   (gg/->image (gg/functions [["ahaversin" m/ahaversin]
                              ["ahacoversin" m/ahacoversin]
                              ["ahavercos" m/ahavercos]
                              ["ahacovercos" m/ahacovercos]] {:x [-0.1 2]}))]])

;; The haversine formula for calculating the square of half the chord length between two points $(\phi_1, \lambda_1)$ and $(\phi_2, \lambda_2)$ on a sphere (where $\phi$ is latitude and $\lambda$ is longitude, in radians) is:

;; $$ a = \sin^2\left(\frac{\phi_2 - \phi_1}{2}\right) + \cos(\phi_1) \cos(\phi_2) \sin^2\left(\frac{\lambda_2 - \lambda_1}{2}\right) $$

;; In this formula:
;; * $a$ is the square of half the chord length of the great-circle arc.
;; * $\phi_1, \phi_2$ are the latitudes of the two points.
;; * $\lambda_1, \lambda_2$ are the longitudes of the two points.

;; The actual great-circle distance $d$ between the points on a sphere of radius $R$ is then given by:
;; $$ d = R \cdot 2 \arcsin(\sqrt{a}) $$

;; The `fastmath.core/haversin` function with four arguments `[lat1 lon1 lat2 lon2]` calculates the value $a$. The `fastmath.core/haversine-dist` function then calculates the distance $d$ assuming $R=1$.

;; Let's calculate `haversin` value for two lat/lon points in degrees `[38.898N, 77.037E]` (White House) and `[48.858N, 2.294W]` (Eiffel Tower).

(m/haversin (m/radians 38.898) (m/radians 77.037) (m/radians 48.858) (m/radians -2.294))

;; Which is equivalent to a distance in km on Earth:

(* 6371.2
   (m/haversine-dist (m/radians 38.898) (m/radians 77.037) (m/radians 48.858) (m/radians -2.294)))


;; ## Power and logarithms

;; This section covers exponential, logarithmic, and power functions, including various specialized and numerically stable variants.

;; ::: {.callout-tip title="Defined functions"}
;; * `exp`, `exp2`, `exp10`, `qexp`
;; * `ln`, `log`, `logb`, `log2`, `log10`, `qlog`
;; * `expm1`, `exprel`, `xexpx`, `xexpy`, `cexpexp`, `expexp` 
;; * `log1p`, `log1pexp`, `log1mexp`, `log2mexp`, `log1psq`, `logexpm1`,`log1pmx`, `logmxp1`
;; * `xlogx`, `xlogy`, `xlog1py`, `cloglog`, `loglog`, `logcosh`
;; * `logaddexp`, `logsubexp`, `logsumexp`, `log2int`
;; * `sigmoid`, `logit`
;; * `sqrt`, `cbrt`, `sq`, `cb`
;; * `safe-sqrt`, `qsqrt`, `rqsqrt`
;; * `pow`, `spow`, `fpow`, `qpow`, `mpow`, `pow2`, `pow3`, `pow10`
;; * `low-2-exp`, `high-2-exp`, `low-exp`, `high-exp`
;; :::

;; ### Exponents

;; * `exp(x)`: Natural exponential function $e^x$.
;; * `exp2(x)`: Base-2 exponential function $2^x$.
;; * `exp10(x)`: Base-10 exponential function $10^x$.
;; * `qexp(x)`: Fast, less accurate version of `exp(x)$.

(kind/table
 [[(gg/->image (gg/functions [["exp" m/exp]
                              ["exp2" m/exp2]
                              ["exp10" m/exp10]] {:x [-3 3] :ylim [nil 5] :title "exp functions"}))
   (gg/->image (gg/function (hp/relative-error m/exp m/qexp) {:x [-10 10] :title "Relative error exp vs qexp"}))]])

(utls/examples-note
  (m/exp 1)
  (m/exp2 3)
  (m/exp10 2)
  (m/qexp 1.0))

;; ### Logarithms

;; * `ln(x)`: Natural logarithm $\ln(x)$. Alias for `log(x)` with one argument.
;; * `log(x)`: Natural logarithm $\ln(x)$. With two arguments, computes $\log_b(x)$.
;; * `logb(b, x)`: Logarithm of $x$ with base $b$, $\log_b(x)$.
;; * `log2(x)`: Base-2 logarithm $\log_2(x)$.
;; * `log2int(x)`: Integer base-2 logarithm, related to the exponent of the floating point representation. Returns `long`.
;; * `log10(x)`: Base-10 logarithm $\log_{10}(x)$.
;; * `qlog(x)`: Fast, less accurate version of `log(x)$.

(kind/table
 [[(gg/->image (gg/functions [["ln" m/ln]
                              ["log2" m/log2]
                              ["log10" m/log10]] {:x [0.1 3] :title "log functions"}))
   (gg/->image (gg/function (hp/relative-error m/log m/qlog) {:x [0.1 10] :title "Relative error log vs qlog"}))]])

(utls/examples-note
  (m/ln 10)
  (m/log 10)
  (m/log 2 8)
  (m/logb 2 8)
  (m/log2 8)
  (m/log2int 8)
  (m/log2int 7.1)
  (m/log10 100)
  (m/qlog 10.0))

;; ### Specialized Log/Exp functions

;; These functions provide numerically stable computations for expressions involving `exp` and `log` especially for small or large arguments.

;; * `expm1(x)`: $e^x - 1$, computed accurately for small $x$.
;; * `exprel(x)`: $(e^x - 1)/x$, computed accurately for small $x$. Returns 1 for $x=0$.
;; * `xexpx(x)`: $x e^x$.
;; * `xexpy(x,y)`: $x e^y$.
;; * `cexpexp(x)`: $1-e^{-e^x}$, inverse of `cloglog`.
;; * `expexp(x)`: $e^{-e^{-x}}$, inverse of `loglog`.
;; * `log1p(x)`: $\ln(1+x)$, computed accurately for small $x$.
;; * `log1pexp(x)`: $\ln(1+e^x)$.
;; * `log1mexp(x)`: $\ln(1-e^x)$, for $x < 0$.
;; * `log2mexp(x)`: $\ln(2-e^x)$.
;; * `log1psq(x)`: $\ln(1+x^2)$, computed accurately for small $x$.
;; * `logexpm1(x)`: $\ln(e^x - 1)$.
;; * `log1pmx(x)`: $\ln(1+x)-x$, computed accurately for small $x$.
;; * `logmxp1(x)`: $\ln(x)-x+1$, computed accurately for $x$ near 1.
;; * `xlogx(x)`: $x \ln(x)$. Returns 0 for $x=0$.
;; * `xlogy(x, y)`: $x \ln(y)$. Returns 0 for $x=0$.
;; * `xlog1py(x, y)`: $x \ln(1+y)$. Returns 0 for $x=0$.
;; * `cloglog(x)`: $\ln(-\ln(1-x))$. Used in complementary log-log models.
;; * `loglog(x)`: $-\ln(-\ln(x))$.
;; * `logcosh(x)`: $\ln(\cosh(x))$.

(kind/table
 [[(gg/->image (gg/function m/expm1 {:x [-2 2] :title "expm1"}))
   (gg/->image (gg/function (hp/relative-error m/expm1 #(dec (m/exp %))) {:x [-2 2] :title "Relative error between expm1(x) and exp(x)-1"}))]
  [(gg/->image (gg/function m/xexpx {:x [-3 2] :title "xexpx"}))
   (gg/->image (gg/functions [["-1" (partial m/xexpy -1)]
                              ["1" (partial m/xexpy 1)]
                              ["2" (partial m/xexpy 2)]] {:x [-3 2] :title "xexpy" :legend-name "x"
                                                          :xlab "y" :ylab "xexpy"}))]
  [(gg/->image (gg/function m/cexpexp {:x [-3 2] :title "cexpexp"}))
   (gg/->image (gg/function m/expexp {:x [-3 2] :title "expexp"}))]
  [(gg/->image (gg/function m/log1p {:x [-1 4] :title "log1p"}))
   (gg/->image (gg/function m/log1pexp {:x [-3 2] :title "log1pexp"}))]
  [(gg/->image (gg/function m/log1mexp {:x [-4 0] :title "log1mexp"}))
   (gg/->image (gg/function m/log2mexp {:x [-4 0] :title "log2mexp"}))]
  [(gg/->image (gg/function m/log1psq {:x [-4 4] :title "log1psq"}))
   (gg/->image (gg/function m/logexpm1 {:x [0 4] :title "logexpm1"}))]
  [(gg/->image (gg/function m/log1pmx {:x [-1 4] :title "log1pmx"}))
   (gg/->image (gg/functions [["-1" (partial m/xlog1py -1)]
                              ["1" (partial m/xlog1py 1)]
                              ["2" (partial m/xlog1py 2)]] {:x [-1 4] :title "xlog1py" :legend-name "x"
                                                            :xlab "y" :ylab "xlog1py"}))]
  [(gg/->image (gg/function m/xlogx {:x [0 4] :title "xlogx"}))
   (gg/->image (gg/functions [["-1" (partial m/xlogy -1)]
                              ["1" (partial m/xlogy 1)]
                              ["2" (partial m/xlogy 2)]] {:x [0 4] :title "xlogy" :legend-name "x"
                                                          :xlab "y" :ylab "xlogy"}))]
  [(gg/->image (gg/function m/cloglog {:x [0 1] :title "cloglog"}))
   (gg/->image (gg/function m/loglog {:x [0 1] :title "loglog"}))]
  [(gg/->image (gg/function m/logcosh {:x [-3 3] :title "logcosh"}))]])

(utls/examples-note
  (m/expm1 1e-9)
  (m/exprel 1e-9)
  (m/xexpx 0.5)
  (m/xexpy 0.5 -0.5)
  (m/cexpexp 1.0)
  (m/expexp 1.0)
  (m/log1p 1e-9)
  (m/log1pexp 0)
  (m/log1mexp -1)
  (m/log2mexp -1)
  (m/log1psq 1e-5)
  (m/logexpm1 1)
  (m/log1pmx 1)
  (m/xlogx 2)
  (m/xlogy 2 5)
  (m/xlog1py 0.5 -0.5)
  (m/cloglog 0.5)
  (m/loglog 0.5)
  (m/logcosh 1))

;; ### Log-sum-exp

;; These functions are used for numerically stable computation of sums and differences of exponents.

;; * `logaddexp(x, y)`: $\ln(e^x + e^y)$.
;; * `logsubexp(x, y)`: $\ln(|e^x - e^y|)$.
;; * `logsumexp(xs)`: $\ln(\sum_{i} e^{x_i})$.

(utls/examples-note
  (m/logaddexp 0 0) ; log(e^0 + e^0) = log(1+1) = log(2)
  (m/logaddexp 100 100) ; log(e^100 + e^100) = log(2 * e^100) = log(2) + 100
  (m/logsubexp 0 0) ; log(|e^0 - e^0|) = log(0) = -Inf
  (m/logsubexp 100 99) ; log(|e^100 - e^99|) = log(e^99(e - 1)) = 99 + log(e-1)
  (m/logsumexp [0 0 0]) ; log(e^0 + e^0 + e^0) = log(3)
  (m/logsumexp [100 100 100])) 

;; ### Sigmoid and Logit

;; Functions used in statistics and machine learning.

;; * `sigmoid(x)`: Sigmoid function $\sigma(x) = \frac{1}{1+e^{-x}}$. Also known as the logistic function.
;; * `logit(x)`: Logit function $\operatorname{logit}(x) = \ln(\frac{x}{1-x})$ for $0 < x < 1$.

(kind/table
 [[(gg/->image (gg/function m/sigmoid {:x [-5 5] :title "sigmoid(x)"}))
   (gg/->image (gg/function m/logit {:x [0.01 0.99] :title "logit(x)"}))]])

(utls/examples-note
  (m/sigmoid 0)
  (m/sigmoid 10)
  (m/sigmoid -10)
  (m/logit 0.5)
  (m/logit 0.9)
  (m/logit 0.1))

;; ### Roots and Powers

;; Basic root and power functions.

;; * `sqrt(x)`: Square root $\sqrt{x}$.
;; * `cbrt(x)`: Cubic root $\sqrt[3]{x}$.
;; * `sq(x)` or `pow2(x)`: Square $x^2$.
;; * `cb(x)` or `pow3(x)`: Cube $x^3$.
;; * `pow10(x)`: $x^{10}$.
;; * `pow(x, exponent)`: $x^{\text{exponent}}$.
;; * `spow(x, exponent)`: Symmetric power of $x$, keeping the sign: $\operatorname{sgn}(x) |x|^{\text{exponent}}$.
;; * `fpow(x, exponent)`: Fast integer power $x^n$, where $n$ is an integer.
;; * `mpow(x, exponent, modulus)`: Modular exponentiation $(x^e \\pmod m)$.
;; * `qsqrt(x)`: Fast, less accurate square root.
;; * `rqsqrt(x)`: Fast, less accurate reciprocal square root $1/\sqrt{x}$.
;; * `safe-sqrt(x)`: Square root, returning 0 for $x \le 0$.

(utls/examples-note
  (m/sqrt 9)
  (m/qsqrt 9)
  (m/rqsqrt 9)
  (m/cbrt 27)
  (m/sq 3)
  (m/cb 3)
  (m/pow 2 3)
  (m/pow 9 0.5)
  (m/pow 0.2 0.3)
  (m/qpow 0.2 0.3)
  (m/spow -8 1/3)
  (m/fpow 2 10)
  (m/mpow 30 112 74)
  (m/safe-sqrt -4))

(kind/table
 [[(gg/->image (gg/functions [["sqrt" m/sqrt]
                              ["qsqrt" m/qsqrt]] {:x [0 2] :title "sqrt"}))
   (gg/->image (gg/function (hp/relative-error m/sqrt m/qsqrt) {:x [0 2] :title "Relative error between sqrt(x) and qsqrt(x)"}))]
  [(gg/->image (gg/functions [["1/sqrt" (comp m// m/sqrt)]
                              ["rqsqrt" m/rqsqrt]] {:x [0.1 2] :title "1/sqrt"}))
   (gg/->image (gg/function (hp/relative-error (comp m// m/sqrt) m/rqsqrt) {:x [0.1 2] :title "Relative error between 1/sqrt(x) and rqsqrt(x)"}))]
  [(gg/->image (gg/functions [["x^(0.3)" #(m/pow % 0.3)]
                              ["cb" m/cb]] {:x [0.1 2] :title "pow"}))
   (gg/->image (gg/function (hp/relative-error #(m/pow % 0.3) #(m/qpow % 0.3)) {:x [0 2] :title "Relative error between pow and qpow for exponent=0.3" :steps 800}))]])

;; ### Exponent/Log Utilities

;; Functions related to powers of 2 and general bases.

;; * `low-2-exp(x)`: Finds the largest integer $n$ such that $2^n \le |x|$. Returns `long`.
;; * `high-2-exp(x)`: Finds the smallest integer $n$ such that $2^n \ge |x|$. Returns `long`.
;; * `low-exp(b, x)`: Finds the largest integer $n$ such that $b^n \le |x|$. Returns `long`.
;; * `high-exp(b, x)`: Finds the smallest integer $n$ such that $b^n \ge |x|$. Returns `long`.

(utls/examples-note
  (m/low-2-exp 10)
  (m/high-2-exp 10)
  (m/low-exp 10 150)
  (m/high-exp 10 150))

;; ## Bitwise operations

;; Functions for performing bitwise operations on `long` primitive types.

;; ::: {.callout-tip title="Defined functions"}
;; * `bit-and`, `bit-or`, `bit-xor`,
;; * `bit-not`, `bit-nand`, `bit-nor`, `bit-xnor`, `bit-and-not`
;; * `bit-set`, `bit-clear`, `bit-flip`, `bit-test`
;; * `<<`, `bit-shift-left`, `>>`, `bit-shift-right`, `>>>`, `unsigned-bit-shift-right`
;; :::

;; ### Logical Bitwise Operations

;; These functions perform standard logical operations on the individual bits of their `long` arguments. They are inlined and accept one or more arguments for multi-arity operations.

;; *   `bit-and`: Bitwise AND ($\land$).
;; *   `bit-or`: Bitwise OR ($\lor$).
;; *   `bit-xor`: Bitwise XOR ($\oplus$).
;; *   `bit-not`: Bitwise NOT ($\sim$).
;; *   `bit-nand`: Bitwise NAND ($\neg (x \land y)$).
;; *   `bit-nor`: Bitwise NOR ($\neg (x \lor y)$).
;; *   `bit-xnor`: Bitwise XNOR ($\neg (x \oplus y)$).
;; *   `bit-and-not`: Bitwise AND with complement of the second argument ($x \land \sim y$).

(utls/examples-note
  (m/bit-and 2r1100 2r1010)
  (m/bit-or 2r1100 2r1010)
  (m/bit-xor 2r1100 2r1010)
  (m/bit-not 2r1100)
  (m/bit-nand 2r1100 2r1010)
  (m/bit-nor 2r1100 2r1010)
  (m/bit-xnor 2r1100 2r1010)
  (m/bit-and-not 2r1100 2r1010)
  (kind/md "Multi-arity examples")
  (m/bit-and 2r1111 2r1100 2r1010)
  (m/bit-or 2r0000 2r1100 2r1010)
  (m/bit-xor 2r1111 2r1100 2r1010))

;; ### Bit Shift Operations

;; These functions shift the bits of a `long` value left or right.

;; *   `<<` / `bit-shift-left`: Signed left shift ($x \ll \text{shift}$). Bits shifted off the left are discarded, and zero bits are shifted in from the right.
;; *   `>>` / `bit-shift-right`: Signed right shift ($x \gg \text{shift}$). Bits shifted off the right are discarded. The sign bit (the leftmost bit) is extended to fill in from the left, preserving the number's sign. 
;; *   `>>>` / `unsigned-bit-shift-right`: Unsigned right shift ($x \ggg \text{shift}$). Bits shifted off the right are discarded. Zero bits are shifted in from the left, regardless of the number's sign.

(utls/examples-note
  (m/<< 2r1100 2)
  (m/bit-shift-left 2r1100 2)
  (m/>> 2r1100 2)
  (m/bit-shift-right 2r1100 2)
  (m/>> -2r1100 2)
  (m/bit-shift-right -2r1100 2)
  (m/>>> 2r1100 2)
  (m/unsigned-bit-shift-right 2r1100 2)
  (m/>>> -2r1100 2) 
  (m/unsigned-bit-shift-right -2r1100 2)) 

;; ### Bit Manipulation

;; Functions to manipulate individual bits within a `long` value.

;; *   `bit-set`: Sets a specific bit at the given index to `1`.
;; *   `bit-clear`: Clears a specific bit at the given index to `0`.
;; *   `bit-flip`: Flips the state of a specific bit at the given index (`0` becomes `1`, `1` becomes `0`).
;; *   `bit-test`: Tests the state of a specific bit at the given index. Returns `true` if the bit is `1`, `false` if it is `0`.

(utls/examples-note
  (m/bit-set 2r1010 1) 
  (m/bit-clear 2r1010 3)
  (m/bit-flip 2r1010 2) 
  (m/bit-test 2r1010 1) 
  (m/bit-test 2r1010 0))

;; ## Floating point

;; Functions for inspecting and manipulating the binary representation of floating-point numbers, and for finding adjacent floating-point values.

;; ::: {.callout-tip title="Defined functions"}
;; * `next-double`, `prev-double`, `ulp`
;; * `double-bits`, `double-high-bits`, `double-low-bits`, `bits->double`
;; * `double-exponent`, `double-significand`
;; :::

;; * `next-double(x)`: Returns the floating-point value adjacent to `x` in the direction of positive infinity. Can take an optional `delta` argument to find the value `delta` steps away.
;; * `prev-double(x)`: Returns the floating-point value adjacent to `x` in the direction of negative infinity. Can take an optional `delta` argument to find the value `delta` steps away.
;; * `ulp(x)`: Returns the size of an ulp of `x` - the distance between this floating-point value and the floating-point value adjacent to it. Formally, it's the spacing between floating-point numbers in the neighborhood of `x`.
;; * `double-bits(x)`: Returns the 64-bit `long` integer representation of the `double` value `x` according to the IEEE 754 floating-point "double format" bit layout.
;; * `double-high-bits(x)`: Returns the high 32 bits of the `double` value `x`'s IEEE 754 representation as a `long`.
;; * `double-low-bits(x)`: Returns the low 32 bits of the `double` value `x`'s IEEE 754 representation as a `long`.
;; * `bits->double(bits)`: Returns the `double` floating-point value corresponding to the given 64-bit `long` integer representation. This is the inverse of `double-bits`.
;; * `double-exponent(x)`: Returns the unbiased exponent of the `double` value `x`.
;; * `double-significand(x)`: Returns the significand (mantissa) of the `double` value `x` as a `long`. This is the 52 explicit bits of the significand for normalized numbers.
;; * `log2int(x)`: Returns an integer approximation of $\log_2(|x|)$. It's closely related to `double-exponent` but includes an adjustment based on the significand to provide a more precise floor-like integer log.


(utls/examples-note
  (m/next-double 0.0)
  (m/next-double -0.0)
  (m/next-double 1.0)
  (m/next-double 1.0 10)
  (m/next-double 1.0e20)
  (m/prev-double 0.0)
  (m/prev-double 1.0)
  (m/prev-double 1.0 10)
  (m/prev-double 1.0e20)
  (m/ulp 1.0)
  (m/ulp 2.0)
  (m/ulp 1.0e20)
  (m/log2int 8.0)
  (m/double-exponent 8.0)
  (m/log2int 7.1)
  (m/double-exponent 7.1))

;; Now let's convert `123.456` to internal representation.

(utls/zp
 (let [d 123.456]
   {:double d
    :bits (m/double-bits d)
    :high (m/double-high-bits d)
    :low (m/double-low-bits d)
    :exponent (m/double-exponent d)
    :significand (m/double-significand d)}))

;; Convert back to a double:

(m/bits->double 4638387860618067575)

;; ## Combinatorics

;; Functions for common combinatorial calculations, including factorials and binomial coefficients.

;; ::: {.callout-tip title="Defined functions"}
;; * `factorial20`, `factorial`, `inv-factorial`, `log-factorial`
;; * `falling-factorial`, `falling-factorial-int`, `rising-factorial`, `rising-factorial-int`
;; * `combinations`, `log-combinations`
;; :::

;; ### Factorials
;;
;; The factorial of a non-negative integer $n$, denoted by $n!$, is the product of all positive integers less than or equal to $n$. $0!$ is defined as $1$.
;;
;; * `factorial20(n)`: Computes $n!$ for $0 \le n \le 20$ using a precomputed table. Returns `long`.
;; * `factorial(n)`: Computes $n!$ for any non-negative integer $n$. For $n > 20$, it uses the Gamma function: $n! = \Gamma(n+1)$. Returns `double`.
;; * `inv-factorial(n)`: Computes the inverse factorial, $\frac{1}{n!}$. Returns `double`.
;; * `log-factorial(n)`: Computes the natural logarithm of the factorial, $\ln(n!) = \ln(\Gamma(n+1))$. Returns `double`.

(utls/examples-note
  (m/factorial 5)
  (m/factorial20 5)
  (m/factorial 21)
  (m/inv-factorial 5)
  (m/log-factorial 5))

;; ### Factorial-like products

;; These functions generalize the factorial to falling (descending) and rising (Pochhammer) products.
;;
;; * `falling-factorial-int(n, x)`: Computes the falling factorial $x^{\underline{n}} = x(x-1)\dots(x-n+1)$ for integer $n \ge 0$.
;; * `falling-factorial(n, x)`: Computes the falling factorial for real $n$, $x^{\underline{n}} = \frac{\Gamma(x+1)}{\Gamma(x-n+1)}$.
;; * `rising-factorial-int(n, x)`: Computes the rising factorial (Pochhammer symbol) $x^{\overline{n}} = x(x+1)\dots(x+n-1)$ for integer $n \ge 0$.
;; * `rising-factorial(n, x)`: Computes the rising factorial for real $n$, $x^{\overline{n}} = \frac{\Gamma(x+n)}{\Gamma(x)}$.

(utls/examples-note
  (m/falling-factorial-int 3 10)
  (m/falling-factorial 3.5 10.0) 
  (m/rising-factorial-int 3 10)
  (m/rising-factorial 3.5 10.0)) 

;; ### Combinations

;; The binomial coefficient $\binom{n}{k}$ represents the number of ways to choose $k$ elements from a set of $n$ distinct elements, without regard to the order of selection.
;;
;; * `combinations(n, k)`: Computes the binomial coefficient $\binom{n}{k} = \frac{n!}{k!(n-k)!}$. Returns `double`.
;; * `log-combinations(n, k)`: Computes the natural logarithm of the binomial coefficient, $\ln\binom{n}{k}$. Returns `double`.

(utls/examples-note
  (m/combinations 10 2)
  (m/combinations 10 8)
  (m/combinations 5 0)
  (m/combinations 5 6)
  (m/log-combinations 10 2)
  (m/log-combinations 1000 500))

;; ## Rank and order

;; Functions for determining the rank of elements within a collection and the order (indices) required to sort it.

;; ::: {.callout-tip title="Defined functions"}
;; * `rank`, `rank1`
;; * `order`
;; :::

;; `rank` computes the rank of each element in a collection. Ranks are 0-based indices indicating the position an element would have if the collection were sorted. It supports various tie-breaking strategies:

;; * `:average`: Assign the average rank to tied elements (default).
;; * `:first`: Assign ranks based on their original order for ties.
;; * `:last`: Assign ranks based on their original order (reverse of `:first`).
;; * `:random`: Assign random ranks for ties.
;; * `:min`: Assign the minimum rank to tied elements.
;; * `:max`: Assign the maximum rank to tied elements.
;; * `:dense`: Assign consecutive ranks without gaps for tied elements.

;; It also supports ascending (default) and descending order.

;; `rank1` is identical to `rank` but returns 1-based indices.

(def data [5 1 8 1 5 1 1 1])

(utls/callout "note" "Examples"
  (kind/md "Ascending")
  (utls/examples
    (m/rank data)
    (m/rank data :dense)
    (m/rank data :min)
    (m/rank data :max)
    (m/rank data :random)
    (m/rank data :random)
    (m/rank data :first)
    (m/rank data :last))
  (kind/md "Descending")
  (utls/examples
    (m/rank data :average true)
    (m/rank data :dense true)
    (m/rank data :min true)
    (m/rank data :max true)
    (m/rank data :random true)
    (m/rank data :random true)
    (m/rank data :first true)
    (m/rank data :last true))
  (kind/md "1-based indices")
  (utls/examples
    (m/rank1 data)))

;; `order` computes the indices that would sort the collection. Applying these indices to the original collection yields the sorted sequence.

(utls/examples-note
  (m/order data)
  (map data (m/order data))
  (m/order data true)
  (map data (m/order data true)))

;; ## Interpolation and mapping

;; Functions for mapping values between numerical ranges (normalization) and for various types of interpolation, smoothing the transition between values. These tools are useful for graphics, animation, signal processing, and data manipulation.

;; ::: {.callout-tip title="Defined functions"}
;; * `norm`, `mnorm`, `cnorm`, `make-norm`
;; * `lerp`, `mlerp`
;; * `smoothstep`
;; * `cos-interpolation`, `smooth-interpolation`, `quad-interpolation`
;; :::

;; ### Mapping and Normalization

;; These functions map a value from one numerical range to another.

;; * `norm(v, start, stop)`: Maps `v` from the range $[start, stop]$ to $[0, 1]$. The formula is $\frac{v - start}{stop - start}$.
;; * `norm(v, start1, stop1, start2, stop2)`: Maps `v` from the range $[start1, stop1]$ to $[start2, stop2]$. The formula is $start2 + (stop2 - start2) \frac{v - start1}{stop1 - start1}$.
;; * `mnorm`: Macro version of `norm` for inlining.
;; * `cnorm`: Constrained version of `norm`. The result is clamped to the target range $[0, 1]$ or $[start2, stop2]$.
;; * `make-norm(start, stop, [dstart, dstop])`: Creates a function that maps values from $[start, stop]$ to $[0, 1]$ or $[dstart, dstop]$.

(utls/examples-note
  (m/norm 5 0 10)     
  (m/norm 15 0 10)    
  (m/norm 5 0 10 100 200)
  (m/cnorm 5 0 10)     
  (m/cnorm 15 0 10)    
  (m/cnorm -5 0 10 100 200) 
  (let [map-fn (m/make-norm 0 10 100 200)] (map-fn 5))
  (let [map-fn (m/make-norm 0 10)] (map-fn 5 100 200)))

;; ### Interpolation

;; Interpolation functions find intermediate values between two points based on a blending factor, often called `t` (for time). The factor `t` typically ranges from 0 to 1, where `t=0` corresponds to the start value and `t=1` to the stop value.

;; * `lerp(start, stop, t)`: Linear interpolation. Blends `start` and `stop` linearly based on `t`. The formula is $start + t \cdot (stop - start)$.
;; * `mlerp`: Macro version of `lerp` for inlining.
;; * `smoothstep(edge0, edge1, x)`: Smooth interpolation between 0 and 1. If $x \le edge0$, returns 0. If $x \ge edge1$, returns 1. Otherwise, it performs a cubic Hermite interpolation for $x$ mapped from $[edge0, edge1]$ to $[0, 1]$. This provides a smooth transition with zero derivative at the edges. The formula for the interpolation step is $t^2 (3 - 2t)$, where $t = \operatorname{cnorm}(x, edge0, edge1)$.
;; * `cos-interpolation(start, stop, t)`: Cosine interpolation. Uses a cosine curve to smooth the transition between `start` and `stop`.
;; * `smooth-interpolation(start, stop, t)`: Uses the `smoothstep` interpolation curve (cubic $3t^2 - 2t^3$) to blend `start` and `stop`.
;; * `quad-interpolation(start, stop, t)`: Quadratic interpolation. Uses a parabolic curve, giving a faster initial and slower final rate of change compared to linear.

(kind/table
 [[(gg/->image (gg/function #(m/lerp 0 10 %)
                            {:x [0 1] :title "lerp(0,10,t)"}))
   (gg/->image (gg/function #(m/smoothstep 2 8 %)
                            {:x [0 10] :title "smoothstep(2,8,x)"}))]])

(gg/->image (gg/functions [["linear" #(m/lerp 0 10 %)]
                           ["cosine" #(m/cos-interpolation 0 10 %)]
                           ["smooth" #(m/smooth-interpolation 0 10 %)]
                           ["quadratic" #(m/quad-interpolation 0 10 %)]]
                          {:x [0 1] :title "Interpolation types (0 to 10)"}))

(utls/examples-note
  (m/lerp 0 10 0.5)
  (m/lerp 0 10 0.0)
  (m/lerp 0 10 1.0)
  (m/smoothstep 2 8 1)
  (m/smoothstep 2 8 5)
  (m/smoothstep 2 8 9)
  (m/cos-interpolation 0 10 0.5)
  (m/smooth-interpolation 0 10 0.5)
  (m/quad-interpolation 0 10 0.5))

;; ## Distance

;; Functions for calculating distances between points or the magnitude (Euclidean norm) of vectors, including specialized functions for geographic distances.


;; ::: {.callout-tip title="Defined functions"}
;; * `dist`, `qdist`, `hypot`, `hypot-sqrt`
;; * `haversine-dist`
;; :::

;; * `dist(x1, y1, x2, y2)`: Calculates the Euclidean distance between two 2D points $(x_1, y_1)$ and $(x_2, y_2)$. Also accepts pairs of coordinates `[x1 y1]` and `[x2 y2]`.
;; $$ \operatorname{dist}((x_1, y_1), (x_2, y_2)) = \sqrt{(x_2 - x_1)^2 + (y_2 - y_1)^2} $$
;; * `qdist`: A faster, less accurate version of `dist` using [[qsqrt]] instead of [[sqrt]].
;; * `hypot(x, y)` and `hypot(x, y, z)`: Calculates the Euclidean norm (distance from the origin) of a 2D or 3D vector $\sqrt{x^2 + y^2}$ or $\sqrt{x^2 + y^2 + z^2}$. This function uses a numerically stable algorithm to avoid potential overflow or underflow issues compared to a direct calculation.
;; $$ \operatorname{hypot}(x, y) = \sqrt{x^2 + y^2} $$
;; $$ \operatorname{hypot}(x, y, z) = \sqrt{x^2 + y^2 + z^2} $$
;; * `hypot-sqrt(x, y)` and `hypot-sqrt(x, y, z)`: Calculates the Euclidean norm using the direct formula $\sqrt{x^2+y^2}$ or $\sqrt{x^2+y^2+z^2}$. This may be less numerically stable than [[hypot]] for inputs with vastly different magnitudes.
;; * `haversine-dist(lat1, lon1, lat2, lon2)`: Calculates the great-circle distance between two points on a sphere given their latitude and longitude (in radians), assuming a sphere with radius $R=1$. This function uses the haversine formula component computed by [[haversin]] (described in the Trigonometry section) and the inverse haversine formula to find the angle, then scales by $R=1$. Also accepts coordinate pairs `[lat1 lon1]` and `[lat2 lon2]`. The distance is $d = 2 \arcsin(\sqrt{a})$, where $a$ is the value computed by `(haversin lat1 lon1 lat2 lon2)`.

(utls/examples-note
  (m/dist 0 0 3 4)
  (m/dist [0 0] [3 4])
  (m/qdist 0 0 3 4)
  (m/hypot 3 4)
  (m/hypot 1.0e150 1.0e250)
  (m/hypot-sqrt 1.0e150 1.0e250)
  (m/hypot 2 3 4)
  (m/haversine-dist (m/radians 38.898) (m/radians 77.037) (m/radians 48.858) (m/radians -2.294))
  (m/haversine-dist [(m/radians 38.898) (m/radians 77.037)] [(m/radians 48.858) (m/radians -2.294)]))

;; ## Intervals

;; Functions for partitioning numerical ranges and grouping data into intervals. These are useful for data analysis, histogram creation, and data binning.

;; ::: {.callout-tip title="Defined functions"}
;; * `slice-range`, `cut`
;; * `co-intervals`, `group-by-intervals`
;; * `sample`
;; :::

;; * `slice-range(cnt)` / `slice-range(start, end, cnt)`: Generates a sequence of `cnt` evenly spaced `double` values between `start` and `end` (inclusive). If only `cnt` is provided, the range `[0.0, 1.0]` is used. If `cnt` is 1, returns the midpoint.
;; * `cut(data, breaks)` / `cut(x1, x2, breaks)`: Divides a numerical range into `breaks` adjacent intervals. The range is either explicitly given by `x1, x2` or determined by the min/max finite values in `data`. Returns a sequence of 2-element vectors `[lower-bound upper-bound]`. The intervals are effectively $[min\_value, p_1], (p_1, p_2], \dots, (p_{breaks-1}, max\_value]$. The lower bound of the first interval is slightly adjusted downwards using [[prev-double]] to ensure the exact minimum is included.

(utls/examples-note
  (m/slice-range 5)
  (m/slice-range -10 10 5)
  (m/cut (range 10) 3)
  (m/cut 0 9 3))

;; * `co-intervals(data, [number, overlap])`: Generates `number` (default 6) overlapping intervals from the sorted finite values in `data`. Each interval aims to contain a similar number of data points. `overlap` (default 0.5) controls the proportion of overlap between consecutive intervals. Returns a sequence of 2-element vectors `[lower-bound upper-bound]`.
;; * `group-by-intervals(intervals, coll)` / `group-by-intervals(coll)`: Groups the values in `coll` into the provided `intervals`. Returns a map where keys are the interval vectors and values are sequences of numbers from `coll` falling into that interval (checked using [[between-?]] for $(lower, upper]$). If no intervals are provided, it first computes them using [[co-intervals]] on `coll`.

(def sample-data (repeatedly 20 #(+ (rand 10) (rand 10) (rand 10))))

(utls/examples-note
  (m/co-intervals sample-data 4)
  (m/co-intervals sample-data 4 0.2))

(utls/zp
 (into (sorted-map) (m/group-by-intervals (m/co-intervals sample-data 4) sample-data)))

;; ### Function sampling

;; `sample(f, number-of-values, [domain-min, domain-max, domain?])`: Samples a function `f` by evaluating it at `number-of-values` evenly spaced points within the specified `[domain-min, domain-max]` range. If `domain?` is true, returns `[x, (f x)]` pairs; otherwise, returns just `(f x)` values. Defaults to sampling 6 values in the `[0, 1]` range.

(utls/examples-note
  (m/sample m/sin 0 m/PI 5)
  (m/sample m/sin 0 m/PI 5 true))

(let [s (m/sample m/sin 0 m/TWO_PI 50 true)
      xs (map first s)
      ys (map second s)]
  (gg/->image (gg/scatter xs ys {:title "Sampling sin(x), 50 points"})))

;; ## Other

;; This section contains a collection of miscellaneous utility functions that don't fit neatly into the preceding categories, including integer number theory functions, boolean operations, data structure conversions, and error calculation utilities.

;; ::: {.callout-tip title="Defined functions"}
;; * `gcd`, `lcm`
;; * `bool-not`, `bool-xor`, `xor`
;; * `identity-long`, `identity-double`
;; * `relative-error`, `absolute-error`
;; * `seq->double-array`, `seq->double-double-array`, `double-array->seq`, `double-double-array->seq`
;; * `use-primitive-operators`, `unuse-primitive-operators`
;; :::

;; ### Number Theory

;; Functions for calculating the greatest common divisor (GCD) and least common multiple (LCM) of integers.

;; * `gcd(a, b)`: Computes the greatest common divisor of two `long` integers `a` and `b`. Uses the binary GCD algorithm (Stein's algorithm), which is typically faster than the Euclidean algorithm for integers.
;; $$ \operatorname{gcd}(a, b) = \max \{ d \in \mathbb{Z}^+ : d \mid a \land d \mid b \} $$
;; * `lcm(a, b)`: Computes the least common multiple of two `long` integers `a` and `b`.
;; $$ \operatorname{lcm}(a, b) = \frac{|a \cdot b|}{\operatorname{gcd}(a, b)} $$

(utls/examples-note
  (m/gcd 12 18)
  (m/gcd 35 49)
  (m/gcd 17 23)
  (m/lcm 12 18)
  (m/lcm 35 49)
  (m/lcm 17 23))

;; ### Boolean and Identity Utilities

;; Basic boolean operations and identity functions for primitive types.

;; * `bool-not(x)`: Primitive boolean NOT. Returns `true` if `x` is logically false, `false` otherwise.
;; * `bool-xor(x, y)` / `xor(x, y)`: Primitive boolean XOR (exclusive OR). Returns `true` if exactly one of `x` or `y` is logically true. Accepts multiple arguments, chaining the XOR operation.
;; * `identity-long(x)`: Returns its `long` argument unchanged. Useful for type hinting or guaranteeing a `long` primitive.
;; * `identity-double(x)`: Returns its `double` argument unchanged. Useful for type hinting or guaranteeing a `double` primitive.

(utls/examples-note
  (m/bool-not true)
  (m/bool-not false)
  (m/xor true false)
  (m/xor true true)
  (m/xor true false true)
  (m/identity-long 5)
  (m/identity-double 5.0))

;; ### Error Calculation

;; Functions to compute the difference between a value and its approximation.

;; * `absolute-error(v, v-approx)`: Computes the absolute difference between a value `v` and its approximation `v-approx`.
;; $$ \operatorname{absolute-error}(v, v_{\text{approx}}) = |v - v_{\text{approx}}| $$
;; * `relative-error(v, v-approx)`: Computes the relative difference between a value `v` and its approximation `v-approx`.
;; $$ \operatorname{relative-error}(v, v_{\text{approx}}) = \left| \frac{v - v_{\text{approx}}}{v} \right| $$

(utls/examples-note
  (m/absolute-error 10.0 10.01)
  (m/relative-error 10.0 10.01)
  (m/absolute-error 1.0e-6 1.01e-6)
  (m/relative-error 1.0e-6 1.01e-6))

;; ### Array and Sequence Conversion

;; Utilities for converting between Clojure sequences and primitive Java arrays (`double[]` and `double[][]`). These are often necessary for interoperation with Java libraries or for performance-critical operations.

;; * `seq->double-array(vs)`: Converts a sequence `vs` into a `double[]` array. If `vs` is already a `double[]`, it is returned directly. If `vs` is a single number, returns a `double[]` of size 1 containing that number.
;; * `double-array->seq(res)`: Converts a `double[]` array `res` into a sequence. This is an alias for Clojure's built-in `seq` function which works correctly for Java arrays.
;; * `seq->double-double-array(vss)`: Converts a sequence of sequences `vss` into a `double[][]` array. If `vss` is already a `double[][]`, it is returned directly.
;; * `double-double-array->seq(res)`: Converts a `double[][]` array `res` into a sequence of sequences.

(utls/examples-note
  (m/seq->double-array [1 2 3])
  (m/double-array->seq (double-array [1.0 2.0 3.0]))
  (m/seq->double-double-array [[1 2] [3 4]])
  (m/double-double-array->seq (into-array (map double-array [[1 2] [3 4]]))))

;; ### Primitive Operators Toggle

;; Macros to replace or restore core Clojure math functions with `fastmath.core` primitive-specialized versions. This can improve performance but should be used carefully as the `fastmath.core` versions have specific type and return value behaviors.

;; * `use-primitive-operators`: Replaces a select set of `clojure.core` functions (like `+`, `-`, `*`, `/`, comparison operators, bitwise operators, `inc`, `dec`, predicates) with their primitive-specialized `fastmath.core` macro equivalents within the current namespace.
;; * `unuse-primitive-operators`: Reverts the changes made by `use-primitive-operators`, restoring the original `clojure.core` functions in the current namespace. Recommended practice is to call this at the end of any namespace that calls `use-primitive-operators` (especially important in Clojure 1.12+).

;; Example usage (typically done at the top or bottom of a namespace):

;; ```clojure
;; (ns my-namespace (:require [fastmath.core :as m]))
;; (m/use-primitive-operators)
;; ... your code using primitive math ...
;; (m/unuse-primitive-operators) ;; or at the end of the file
;; ```


;; ## Constants

(utls/gen-constants 'fastmath.core #{'jvm-version})

;; ## Reference

(codox/make-public-fns-table-clay 'fastmath.core)

