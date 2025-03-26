^:kindly/hide-code
(ns core
  (:require [fastmath.dev.codox :as codox]
            [fastmath.dev.ggplot :as gg]
            [fastmath.dev.clay :as utls]
            [fastmath.dev.helpers :as hp]
            [fastmath.core :as m]
            [scicloj.kindly.v4.kind :as kind]))

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

;; ## Predicates

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

;; ## Trigonometry

;; ::: {.callout-tip title="Defined functions"}
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

^:kind/table
[[(gg/->image (gg/function m/sin {:x [-4 4] :title "sin(x)"}))
  (gg/->image (gg/function m/cos {:x [-4 4] :title "cos(x)"}))]
 [(gg/->image (gg/function m/tan {:x [-4 4] :ylim [-3 3] :title "tan(x)"}))
  (gg/->image (gg/function m/cot {:x [-4 4] :ylim [-3 3] :title "cot(x)"}))]
 [(gg/->image (gg/function m/sec {:x [-4 4] :ylim [-3 3] :title "sec(x)"}))
  (gg/->image (gg/function m/csc {:x [-4 4] :ylim [-3 3] :title "csc(x)"}))]]

^:kind/table
[[(gg/->image (gg/function m/qsin {:x [-4 4] :title "qsin(x)"}))
  (gg/->image (gg/function (hp/relative-error m/sin m/qsin) {:x [-4 4] :title "Relative error between sin(x) and qsin(x)."}))
  ]
 [(gg/->image (gg/function m/qcos {:x [-4 4] :title "qcos(x)"}))
  (gg/->image (gg/function (hp/relative-error m/cos m/qcos) {:x [-4 4] :title "Relative error between cos(x) and qcos(x)."}))]]

^:kind/table
[[(gg/->image (gg/function m/asin {:x [-2 2] :title "asin(x)"}))
  (gg/->image (gg/function m/acos {:x [-2 2] :title "acos(x)"}))]
 [(gg/->image (gg/function m/atan {:x [-4 4] :title "atan(x)"}))
  (gg/->image (gg/function m/acot {:x [-4 4] :title "acot(x)"}))]
 [(gg/->image (gg/function m/asec {:x [-4 4] :ylim [-3 3] :title "asec(x)"}))
  (gg/->image (gg/function m/acsc {:x [-4 4] :ylim [-3 3] :title "acsc(x)"}))]]


^:kind/table
[[(gg/->image (gg/function m/sinh {:x [-4 4] :title "sinh(x)"}))
  (gg/->image (gg/function m/cosh {:x [-4 4] :title "cosh(x)"}))]
 [(gg/->image (gg/function m/tanh {:x [-4 4] :title "tanh(x)"}))
  (gg/->image (gg/function m/coth {:x [-4 4] :ylim [-3 3] :title "coth(x)"}))]
 [(gg/->image (gg/function m/sech {:x [-4 4] :title "sech(x)"}))
  (gg/->image (gg/function m/csch {:x [-4 4] :ylim [-3 3] :title "csch(x)"}))]]

^:kind/table
[[(gg/->image (gg/function m/asinh {:x [-4 4] :title "asinh(x)"}))
  (gg/->image (gg/function m/acosh {:x [ 0 4] :title "acosh(x)"}))]
 [(gg/->image (gg/function m/atanh {:x [-2 2] :title "atanh(x)"}))
  (gg/->image (gg/function m/acoth {:x [-4 4] :ylim [-3 3] :title "acoth(x)"}))]
 [(gg/->image (gg/function m/asech {:x [ 0 2] :title "asech(x)"}))
  (gg/->image (gg/function m/acsch {:x [-4 4] :ylim [-3 3] :title "acsch(x)"}))]]

(gg/->image (gg/function m/sinc {:x [-4 4] :title "sinc(x)"}))


;; ## Power and logarithms

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
;; * `pow`, `spow`, `fpow`, `qpow`, `pow2`, `pow3`, `pow10`
;; * `low-2-exp`, `high-2-exp`, `low-exp`, `high-exp`
;; :::


;; ## Bitwise operations

;; ::: {.callout-tip title="Defined functions"}
;; * `bit-and`, `bit-or`, `bit-xor`,
;; * `bit-not`, `bit-nand`, `bit-nor`, `bit-xnor`, `bit-and-not`
;; * `bit-set`, `bit-clear`, `bit-flip`, `bit-test`
;; * `<<`, `bit-shift-left`, `>>`, `bit-shift-right`, `>>>`, `unsigned-bit-shift-right`
;; :::

;; ## Floating point

;; ::: {.callout-tip title="Defined functions"}
;; * `next-double`, `prev-double`, `ulp`
;; * `double-bits`, `double-high-bits`, `double-low-bits`, `bits->double`
;; * `double-exponent`, `double-significand`
;; :::

;; ## Combinatorics

;; ::: {.callout-tip title="Defined functions"}
;; * `factorial20`, `factorial`, `inv-factorial`, `log-factorial`
;; * `falling-factorial`, `falling-factorial-int`, `rising-factorial`, `rising-factorial-int`
;; * `combinations`, `log-combinations`
;; :::

;; ## Rank and order

;; ::: {.callout-tip title="Defined functions"}
;; * `rank`, `rank1`
;; * `order`
;; :::

;; ## Interpolation and mapping

;; ::: {.callout-tip title="Defined functions"}
;; * `norm`, `mnorm`, `cnorm`, `make-norm`
;; * `lerp`, `mlerp`
;; * `smoothstep`
;; * `cos-interpolation`, `smooth-interpolation`, `quad-interpolation`
;; :::

;; ## Distance

;; ::: {.callout-tip title="Defined functions"}
;; * `dist`, `qdist`, `hypot`, `hypot-sqrt`
;; * `haversine-dist`
;; :::

;; ## Intervals

;; ::: {.callout-tip title="Defined functions"}
;; * `slice-range`, `cut`
;; * `co-intervals`, `group-by-intervals`
;; :::

;; ### Other

;; ::: {.callout-tip title="Defined functions"}
;; * `gcd`, `lcm`
;; * `bool-not`, `bool-xor`, `xor`
;; * `identity-long`, `identity-double`
;; * `radians`, `degrees`
;; * `sample`
;; * `relative-error`, `absolute-error`
;; * `seq->double-array`, `seq->double-double-array`, `double-array->seq`, `double-double-array->seq`
;; * `use-primitive-operators`, `unuse-primitive-operators`
;; :::

;; ## Constants

(utls/gen-constants 'fastmath.core #{'jvm-version})

;; ## Reference

(codox/make-public-fns-table-clay 'fastmath.core)

