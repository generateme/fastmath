(ns ^:no-doc fastmath.calculus.autodiff.expint
  "Internal: the exponential, logarithmic and trigonometric integrals of `fastmath.special` for numbers and duals of [[fastmath.calculus.autodiff]], dispatched by `instance?`, and their derivative rules (the integrands, written without removable singularities). The order of [[En]] has no derivative."
  (:require [fastmath.core :as m]
            [fastmath.special :as s]
            [fastmath.calculus.autodiff.dual :as d]
            [fastmath.calculus.autodiff.ops :refer [add sub mul div neg exp log cos]]
            [fastmath.calculus.autodiff.elementary :as el]
            [fastmath.calculus.autodiff.logexp :as lx]
            [fastmath.calculus.autodiff.rules :as r]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(declare dual-E0 dual-E1 dual-En dual-Ei dual-Ein dual-Si dual-shifted-Si dual-Ci dual-Cin dual-li dual-offset-li)

;; ## Operations on numbers and duals

(defn E0 "`exp(-x)/x`." [x] (r/dispatch dual-E0 s/E0 x))
(defn E1 "Exponential integral `E1`." [x] (r/dispatch dual-E1 s/E1 x))
(defn En "Generalized exponential integral of order `n`." [n x] (r/dispatch dual-En s/En n x))
(defn Ei "Exponential integral `Ei`." [x] (r/dispatch dual-Ei s/Ei x))
(defn Ein "Entire exponential integral." [x] (r/dispatch dual-Ein s/Ein x))
(defn Si "Sine integral." [x] (r/dispatch dual-Si s/Si x))
(defn si "`Si(x) - pi/2`." [x] (r/dispatch dual-shifted-Si s/si x))
(defn Ci "Cosine integral, `x >= 0`." [x] (r/dispatch dual-Ci s/Ci x))
(defn Cin "Entire cosine integral." [x] (r/dispatch dual-Cin s/Cin x))
(defn li "Logarithmic integral." [x] (r/dispatch dual-li s/li x))
(defn Li "Offset logarithmic integral, `li(x) - li(2)`." [x] (r/dispatch dual-offset-li s/Li x))

;; ## Derivative rules

(defn- dual-E0 [x] (r/chain [v x, y (E0 v)] :E0 (neg (mul y (add 1.0 (div 1.0 v))))))
(defn- dual-E1 [x] (r/chain [v x, _y (E1 v)] :E1 (neg (E0 v))))
(defn- dual-En [n x] (r/chain [n n, x x, _y (En n x)] :En nil (neg (En (sub n 1.0) x))))
(defn- dual-Ei
  "`exp(x)/x`; above `709`, where `exp(x)` overflows before the quotient does (up to `x = 716.4`), as
  `exp(x/2) (exp(x/2)/x)`."
  [x]
  (r/chain [v x, _y (Ei v)] :Ei
    (if (m/> (d/scalar v) 709.0)
      (let [h (exp (mul 0.5 v))] (mul h (div h v)))
      (div (exp v) v))))

(defn- dual-Ein
  "`(1 - exp(-x))/x`, written as `exprel(-x)`: finite at `0` (`1`) and for large `x` (where `exp(-x) exprel(x)` is
  `0 * Inf`)."
  [x]
  (r/chain [v x, _y (Ein v)] :Ein (lx/exprel (neg v))))

;; The trigonometric integrands decay like `1/x`: at an infinite argument their limit `0` is used (`sin(Inf)` and
;; `cos(Inf)` are `NaN`, while the values of the integrals are finite, or infinite for `Cin`).

(defn- sine-integrand "`sin(x)/x` as `sinc(x/pi)`, `1` at `0`, `0` at an infinite `x`." [v] (if (d/scalar-infinite? v) 0.0 (el/sinc (div v m/PI))))

(defn- dual-Si [x] (r/chain [v x, _y (Si v)] :Si (sine-integrand v)))
(defn- dual-shifted-Si [x] (r/chain [v x, _y (si v)] :si (sine-integrand v)))
(defn- dual-Ci
  "`cos(x)/x`, `0` at `Inf`. The denominator is `x + 0.0`: `-0.0` becomes `0.0` (`Ci` is defined for `x >= 0` and its
  value at `-0.0` is the value at `0`, so the derivative is `+Inf`, not `-Inf`), every other value is unchanged."
  [x]
  (r/chain [v x, _y (Ci v)] :Ci (if (d/scalar-infinite? v) 0.0 (div (cos v) (add v 0.0)))))

(defn- dual-Cin
  "`(1 - cos(x))/x` as `(x/2) sinc(x/(2 pi))^2`, `0` at `0` and at an infinite `x`."
  [x]
  (r/chain [v x, _y (Cin v)] :Cin
    (if (d/scalar-infinite? v) 0.0 (let [s (el/sinc (div v m/TWO_PI))] (mul (mul 0.5 v) (mul s s))))))

(defn- dual-li [x] (r/chain [v x, _y (li v)] :li (div 1.0 (log v))))
(defn- dual-offset-li [x] (r/chain [v x, _y (Li v)] :Li (div 1.0 (log v))))
