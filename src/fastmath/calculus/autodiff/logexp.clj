(ns ^:no-doc fastmath.calculus.autodiff.logexp
  "Internal: the exponential and logarithmic variants of `fastmath.core` for numbers and duals of [[fastmath.calculus.autodiff]] (base 2 and 10, logarithms with a base, `x log y` forms, the log-sum-exp family, the sigmoid and link functions, the smooth maximum), dispatched by `instance?`, and their derivative rules: the forms of DiffRules.jl and LogExpFunctions.jl where these define the function, otherwise derived from the definition (`logb` in the base, `exprel`, `expexp`, `cexpexp`, `cloglog`, `loglog`, the smooth maximum)."
  (:require [fastmath.core :as m]
            [fastmath.calculus.common :as cc]
            [fastmath.calculus.autodiff.dual :as d]
            [fastmath.calculus.autodiff.ops :as ops :refer [add sub mul div neg exp log]]
            [fastmath.calculus.autodiff.elementary :as el]
            [fastmath.calculus.autodiff.rules :as r])
  (:import [fastmath.calculus.autodiff.dual Dual]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(declare dual-exp2 dual-exp10 dual-logb dual-xlogx dual-xlogy dual-xlog1py dual-log1pexp dual-logexpm1 dual-log1mexp
         dual-log1pmx dual-log1psq dual-logcosh dual-exprel dual-expexp dual-cexpexp dual-cloglog dual-loglog
         dual-sigmoid dual-logit dual-logaddexp dual-logsubexp dual-logsumexp dual-smooth-max)

;; ## Operations on numbers and duals

(defn exp2 "`2^x`." [x] (r/dispatch dual-exp2 m/exp2 x))
(defn exp10 "`10^x`." [x] (r/dispatch dual-exp10 m/exp10 x))
(defn logb "Logarithm of `x` with the base `b`." [b x] (r/dispatch dual-logb m/logb b x))
(defn xlogx "`x log(x)`, `0` at `0`." [x] (r/dispatch dual-xlogx m/xlogx x))
(defn xlogy "`x log(y)`, `0` for `x = 0` and `y` not `NaN`." [x y] (r/dispatch dual-xlogy m/xlogy x y))
(defn xlog1py "`x log(1 + y)`, `0` for `x = 0` and `y` not `NaN`." [x y] (r/dispatch dual-xlog1py m/xlog1py x y))
(defn log1pexp "`log(1 + exp(x))`." [x] (r/dispatch dual-log1pexp m/log1pexp x))
(defn logexpm1 "`log(exp(x) - 1)`." [x] (r/dispatch dual-logexpm1 m/logexpm1 x))
(defn log1mexp "`log(1 - exp(x))`." [x] (r/dispatch dual-log1mexp m/log1mexp x))
(defn log1pmx "`log(1 + x) - x`." [x] (r/dispatch dual-log1pmx m/log1pmx x))
(defn log1psq "`log(1 + x^2)`." [x] (r/dispatch dual-log1psq m/log1psq x))
(defn logcosh "`log(cosh(x))`." [x] (r/dispatch dual-logcosh m/logcosh x))
(defn exprel "`(exp(x) - 1) / x`, `1` at `0`." [x] (r/dispatch dual-exprel m/exprel x))
(defn expexp "`exp(-exp(-x))`." [x] (r/dispatch dual-expexp m/expexp x))
(defn cexpexp "`1 - exp(-exp(x))`." [x] (r/dispatch dual-cexpexp m/cexpexp x))
(defn cloglog "`log(-log(1 - x))`." [x] (r/dispatch dual-cloglog m/cloglog x))
(defn loglog "`-log(-log(x))`." [x] (r/dispatch dual-loglog m/loglog x))
(defn sigmoid "`1 / (1 + exp(-x))`." [x] (r/dispatch dual-sigmoid m/sigmoid x))
(defn logit "`log(x / (1 - x))`." [x] (r/dispatch dual-logit m/logit x))
(defn logaddexp "`log(exp(x) + exp(y))`." [x y] (r/dispatch dual-logaddexp m/logaddexp x y))
(defn logsubexp "`log(|exp(x) - exp(y)|)`." [x y] (r/dispatch dual-logsubexp m/logsubexp x y))

(defn- non-empty-vec
  "`xs` as a vector; an empty sequence throws."
  [xs]
  (let [v (vec xs)]
    (when (empty? v) (cc/invalid-option! :xs xs "a non-empty sequence"))
    v))

(defn logsumexp
  "`log(sum(exp(xs)))` of a non-empty sequence."
  [xs]
  (let [xs (non-empty-vec xs)]
    (if (some (fn [x] (instance? Dual x)) xs)
      (dual-logsumexp xs)
      (m/logsumexp xs))))

(defn smooth-max
  "Smooth maximum of a non-empty sequence, `fastmath.core/smooth-max`."
  [xs alpha family]
  (let [xs (non-empty-vec xs)]
    (if (or (instance? Dual alpha) (some (fn [x] (instance? Dual x)) xs))
      (dual-smooth-max xs alpha family)
      (m/smooth-max xs (.doubleValue ^Number alpha) family))))

;; ## Derivative rules

(defn- dual-exp2 [x] (r/chain [v x, y (exp2 v)] :exp2 (mul y m/LN2)))
(defn- dual-exp10 [x] (r/chain [v x, y (exp10 v)] :exp10 (mul y m/LN10)))

(defn- dual-logb
  "`log(x) / log(b)`: `-log(x) / (b log(b)^2)` in `b` (as `-y / (b log(b))`), `1 / (x log(b))` in `x`."
  [b x]
  (r/chain [b b, x x, y (logb b x)] :logb (neg (div y (mul b (log b)))) (div 1.0 (mul x (log b)))))

(defn- dual-xlogx [x] (r/chain [v x, _y (xlogx v)] :xlogx (add 1.0 (log v))))

(defn- x-over
  "`x / y`, `0` when `x` is zero (with all its partials, as `iszero` of ForwardDiff.jl) and `y` is not `NaN`."
  [x y]
  (if (and (d/zero-element? x) (not (m/nan? (d/scalar y)))) 0.0 (div x y)))

(defn- dual-xlogy [x y] (r/chain [a x, b y, _z (xlogy a b)] :xlogy (log b) (x-over a b)))
(defn- dual-xlog1py [x y] (r/chain [a x, b y, _z (xlog1py a b)] :xlog1py (el/log1p b) (x-over a (add 1.0 b))))

(defn- dual-log1pexp [x] (r/chain [v x, _y (log1pexp v)] :log1pexp (sigmoid v)))
(defn- dual-logexpm1 [x] (r/chain [v x, y (logexpm1 v)] :logexpm1 (exp (sub v y))))
(defn- dual-log1mexp [x] (r/chain [v x, y (log1mexp v)] :log1mexp (neg (exp (sub v y)))))
(defn- dual-log1pmx [x] (r/chain [v x, _y (log1pmx v)] :log1pmx (neg (div v (add 1.0 v)))))
(defn- dual-log1psq [x] (r/chain [v x, _y (log1psq v)] :log1psq (div (mul 2.0 v) (add 1.0 (mul v v)))))
(defn- dual-logcosh [x] (r/chain [v x, _y (logcosh v)] :logcosh (el/tanh v)))

;; `exprel' = (exp(x) - exprel(x)) / x`, which cancels near `0` (0/0 at `0`, the limit is `1/2`): below the threshold
;; a Taylor series written with the operations on duals; above it, for a positive `x`, as `exprel(x) (1 - 1/x) + 1/x`
;; (more accurate on `[0.5, 1]` than the first form, `1.5e-16` against `1.1e-15`). Measured against mpmath
;; on `|x|` in `[1e-8, 10]`: at most `1.8e-16` (first derivative) and `2.4e-16` (second) below the threshold, `1.1e-15`
;; and `5.9e-15` above it.

(def ^{:private true :const true :tag 'double} exprel-series-threshold 0.5)

(def ^:private exprel-derivative-coefficients
  "`exprel'(x) = sum (j + 1) x^j / (j + 2)!`."
  (mapv (fn [^long j] (m// (double (m/long-inc j)) (m/factorial (double (m/long-add j 2))))) (range 20)))

(defn- exprel-derivative
  [v y]
  (let [s (d/scalar v)]
    (cond
      (m/< (m/abs s) exprel-series-threshold) (ops/horner exprel-derivative-coefficients v)
      (m/pos? s) (add (mul y (sub 1.0 (div 1.0 v))) (div 1.0 v))
      :else (div (sub (exp v) y) v))))

(defn- dual-exprel [x] (r/chain [v x, y (exprel v)] :exprel (exprel-derivative v y)))

(defn- dual-expexp [x] (r/chain [v x, _y (expexp v)] :expexp (exp (sub (neg v) (exp (neg v))))))
(defn- dual-cexpexp [x] (r/chain [v x, _y (cexpexp v)] :cexpexp (exp (sub v (exp v)))))
(defn- dual-cloglog [x] (r/chain [v x, _y (cloglog v)] :cloglog (neg (div 1.0 (mul (sub 1.0 v) (el/log1p (neg v)))))))
(defn- dual-loglog [x] (r/chain [v x, _y (loglog v)] :loglog (neg (div 1.0 (mul v (log v))))))

(defn- dual-sigmoid [x] (r/chain [v x, y (sigmoid v)] :sigmoid (mul y (sub 1.0 y))))
(defn- dual-logit [x] (r/chain [v x, _y (logit v)] :logit (div 1.0 (mul v (sub 1.0 v)))))

(defn- dual-logaddexp
  [x y]
  (r/chain [a x, b y, z (logaddexp a b)] :logaddexp (exp (sub a z)) (exp (sub b z))))

(defn- dual-logsubexp
  "The signs of the partials follow the larger argument, as in LogExpFunctions.jl."
  [x y]
  (r/chain [a x, b y, z (logsubexp a b)] :logsubexp
    (if (m/> (d/scalar a) (d/scalar b)) (exp (sub a z)) (neg (exp (sub a z))))
    (if (m/> (d/scalar a) (d/scalar b)) (neg (exp (sub b z))) (exp (sub b z)))))

(defn- softmax-derivative
  "The partial derivative of `logsumexp` in an element as a function of its value, for the values `vs` and the value
  `y`: `exp(v - y)`. For `y = +Inf` with exactly one element `+Inf` the limit, `1` in that element and `0` in the others
  (constants, as LogExpFunctions.jl; the formula is `exp(Inf - Inf) = NaN`)."
  [vs y]
  (if (and (m/== ##Inf (d/scalar y)) (m/== 1 (count (filter (fn [v] (m/== ##Inf (d/scalar v))) vs))))
    (fn [v] (if (m/== ##Inf (d/scalar v)) 1.0 0.0))
    (fn [v] (exp (sub v y)))))

(defn- dual-logsumexp
  "The partials are the softmax, `exp(x_i - logsumexp(xs))`."
  [xs]
  (r/chain-seq :logsumexp xs [] (fn [vs _] (logsumexp vs)) (fn [vs _ y] (softmax-derivative vs y))))

;; ## Smooth maximum
;;
;; The families with a smooth closed form derivative. The pairwise `:smu` is a fold of a binary rule, its value (`0.5
;; (a + b +/- sqrt((a - b)^2 + 1/|alpha|))`) is computed as in `fastmath.core/smooth-max`, so the folded value is the
;; same double.

(def ^:private differentiable-families #{:lse :boltzmann :mellowmax :p-norm :smu})

(defn- smu-pair-value
  "The value of a pair, `fastmath.core/smooth-max` of `[a b]` (infinite and far apart arguments included)."
  ^double [^double a ^double b ^double alpha]
  (m/smooth-max [a b] alpha :smu))

(declare dual-smu-pair)

(defn- smu-pair [a b alpha] (r/dispatch dual-smu-pair smu-pair-value a b alpha))

(defn- smu-slope
  "`q = (a - b) / sqrt((a - b)^2 + 1/|alpha|)`, the partials of a pair are `(1 +/- q) / 2`. Beyond `|a - b| = 1e150` (and
  for infinite arguments) as `sign(a - b) / sqrt(1 + (e/|d|)/|d|)`, where `(a - b)^2` would overflow."
  [a b alpha]
  (let [d (sub a b)
        epsilon (el/abs (div 1.0 alpha))
        sd (d/scalar d)]
    (if (m/> (m/abs sd) 1.0e150)
      (let [ad (el/abs d)]
        (div (if (m/neg? sd) -1.0 1.0) (el/sqrt (add 1.0 (div (div epsilon ad) ad)))))
      (div d (el/sqrt (add (el/sq d) epsilon))))))

(defn- dual-smu-pair
  [a b alpha]
  (r/chain [a a, b b, alpha alpha, _y (smu-pair a b alpha)] :smooth-max
    (let [q (smu-slope a b alpha)]
      (mul 0.5 (if (m/pos? (d/scalar alpha)) (add 1.0 q) (sub 1.0 q))))
    (let [q (smu-slope a b alpha)]
      (mul 0.5 (if (m/pos? (d/scalar alpha)) (sub 1.0 q) (add 1.0 q))))
    nil))

(defn- softmax-weight
  "The softmax weight `exp(alpha v) / sum(exp(alpha vs))` as a function of a value `v`, shifted by the hard extreme `m`
  of the scalar values (the maximum for positive `alpha`, the minimum otherwise; a constant shift, so derivatives of
  any order are exact) to avoid the overflow of `alpha v`. For an infinite `m` or `alpha` the limit when exactly one
  value equals `m`: `1` there and `0` elsewhere; `NaN` otherwise (no limit)."
  [vs alpha]
  (let [a (double (d/scalar alpha))
        ss (map (fn [v] (double (d/scalar v))) vs)
        m (double (reduce (fn [^double x ^double y] (if (m/pos? a) (m/max x y) (m/min x y))) ss))]
    (if (or (m/inf? m) (m/inf? a))
      (if (m/== 1 (count (filter (fn [^double s] (m/== s m)) ss)))
        (fn [v] (if (m/== (d/scalar v) m) 1.0 0.0))
        (fn [_] ##NaN))
      (let [den (reduce add (map (fn [v] (exp (mul alpha (sub v m)))) vs))]
        (fn [v] (div (exp (mul alpha (sub v m))) den))))))

(defn- smooth-max-derivative
  "The partial derivative in an element as a function of its value `v`, for the values `vs`, the sharpness `alpha` and
  the smooth maximum `y`. `:lse` and `:mellowmax` have the softmax as the gradient, `:boltzmann` the softmax times
  `1 + alpha (v - y)` (a zero weight gives 0, an infinite `y` or `alpha` the softmax limit). `:smu` comes here with one element
  only (its value is the element)."
  [family vs alpha y]
  (case family
    :smu (fn [_] 1.0)
    (:lse :mellowmax) (softmax-weight vs alpha)
    :boltzmann (let [w (softmax-weight vs alpha)]
                 (if (or (m/inf? (d/scalar y)) (m/inf? (d/scalar alpha)))
                   w
                   (fn [v] (let [wv (w v)]
                             (if (m/zero? (d/scalar wv)) 0.0 (mul wv (add 1.0 (mul alpha (sub v y)))))))))
    :p-norm (let [one-infinite? (and (m/inf? (d/scalar y))
                                     (m/== 1 (count (filter (fn [v] (m/inf? (d/scalar v))) vs))))]
              ;; a single infinite magnitude dominates: its partial is its sign, the others vanish
              (fn [v] (if (and one-infinite? (m/inf? (d/scalar v)))
                        (if (el/sign-bit? v) -1.0 1.0)
                        (let [r (ops/pow (div (el/abs v) y) (sub alpha 1.0))]
                          (if (el/sign-bit? v) (neg r) r)))))))

(defn- dual-smooth-max
  [xs alpha family]
  (when-not (contains? differentiable-families family)
    (throw (ex-info (str "The derivative of smooth-max is available for the families " differentiable-families
                         " only, not " family ".")
                    {:function :smooth-max :family family})))
  (if (and (= :smu family) (next xs))
    (reduce (fn [a b] (smu-pair a b alpha)) xs)
    (r/chain-seq :smooth-max xs [[:alpha alpha]]
                 (fn [vs [al]] (smooth-max vs al family))
                 (fn [vs [al] y] (smooth-max-derivative family vs al y)))))
