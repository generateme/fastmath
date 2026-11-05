(ns ^:no-doc fastmath.calculus.autodiff.special
  "Internal: the gamma and beta functions, the error functions and their inverses, and the Lambert W functions of numbers and duals of [[fastmath.calculus.autodiff]], dispatched by `instance?`, and their derivative rules. The shape parameters of the incomplete functions and the order of [[polygamma]] have no derivative."
  (:require [fastmath.core :as m]
            [fastmath.special :as s]
            [fastmath.calculus.autodiff.dual :as d]
            [fastmath.calculus.autodiff.ops :as ops :refer [add sub mul div neg pow exp log]]
            [fastmath.calculus.autodiff.elementary :as el]
            [fastmath.calculus.autodiff.rules :as r])
  (:import [fastmath.calculus.autodiff.dual Dual]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(declare dual-gamma dual-log-gamma dual-log-gamma-1p dual-inv-gamma-1pm1 dual-digamma dual-trigamma dual-polygamma
         dual-beta dual-log-beta dual-lower-incomplete-gamma dual-upper-incomplete-gamma dual-regularized-gamma-p
         dual-regularized-gamma-q dual-incomplete-beta dual-regularized-beta
         dual-erf dual-erf2 dual-erfc dual-inv-erf dual-inv-erfc dual-lambert-W dual-lambert-W-1)

;; ## Operations on numbers and duals

(defn gamma "Gamma function." [x] (r/dispatch dual-gamma s/gamma x))
(defn log-gamma "Logarithm of the absolute value of the gamma function." [x] (r/dispatch dual-log-gamma s/log-gamma x))
(defn log-gamma-1p "`log-gamma(1 + x)`." [x] (r/dispatch dual-log-gamma-1p s/log-gamma-1p x))
(defn inv-gamma-1pm1 "`1 / gamma(1 + x) - 1`." [x] (r/dispatch dual-inv-gamma-1pm1 s/inv-gamma-1pm1 x))
(defn digamma "Digamma function." [x] (r/dispatch dual-digamma s/digamma x))
(defn trigamma "Trigamma function." [x] (r/dispatch dual-trigamma s/trigamma x))
(defn polygamma "Polygamma function of order `m`." [m x] (r/dispatch dual-polygamma s/polygamma m x))
(defn beta "Beta function." [p q] (r/dispatch dual-beta s/beta p q))
(defn log-beta "Logarithm of the beta function." [p q] (r/dispatch dual-log-beta s/log-beta p q))
(defn lower-incomplete-gamma "Lower incomplete gamma function." [s x] (r/dispatch dual-lower-incomplete-gamma s/lower-incomplete-gamma s x))
(defn upper-incomplete-gamma "Upper incomplete gamma function." [s x] (r/dispatch dual-upper-incomplete-gamma s/upper-incomplete-gamma s x))
(defn regularized-gamma-p "Regularized lower incomplete gamma function." [a x] (r/dispatch dual-regularized-gamma-p s/regularized-gamma-p a x))
(defn regularized-gamma-q "Regularized upper incomplete gamma function." [a x] (r/dispatch dual-regularized-gamma-q s/regularized-gamma-q a x))
(defn incomplete-beta "Incomplete beta function, `x` first." [x a b] (r/dispatch dual-incomplete-beta s/incomplete-beta x a b))
(defn regularized-beta "Regularized incomplete beta function, `x` first." [x a b] (r/dispatch dual-regularized-beta s/regularized-beta x a b))
(defn erf
  "Error function; with two arguments `erf(x2) - erf(x1)`."
  ([x] (r/dispatch dual-erf s/erf x))
  ([x1 x2] (r/dispatch dual-erf2 s/erf x1 x2)))
(defn erfc "Complementary error function." [x] (r/dispatch dual-erfc s/erfc x))
(defn inv-erf "Inverse error function." [x] (r/dispatch dual-inv-erf s/inv-erf x))
(defn inv-erfc "Inverse complementary error function." [x] (r/dispatch dual-inv-erfc s/inv-erfc x))
(defn lambert-W "Lambert W, the principal branch." [x] (r/dispatch dual-lambert-W s/lambert-W x))
(defn lambert-W-1 "Lambert W, the lower branch." [x] (r/dispatch dual-lambert-W-1 s/lambert-W-1 x))

;; ## Derivative rules: gamma, polygamma and beta (the forms of DiffRules.jl)

(defn- dual-gamma [x] (r/chain [v x, y (gamma v)] :gamma (mul (digamma v) y)))
(defn- dual-log-gamma [x] (r/chain [v x, _y (log-gamma v)] :log-gamma (digamma v)))
(defn- dual-log-gamma-1p [x] (r/chain [v x, _y (log-gamma-1p v)] :log-gamma-1p (digamma (add 1.0 v))))
(defn- dual-inv-gamma-1pm1
  "The derivative `-digamma(1 + x) / gamma(1 + x)`; not with `1 + inv-gamma-1pm1(x)`, which loses all digits when
  `1 / gamma(1 + x)` is below the rounding of `1` (`x > 18`)."
  [x]
  (r/chain [v x, _y (inv-gamma-1pm1 v)] :inv-gamma-1pm1 (neg (div (digamma (add 1.0 v)) (gamma (add 1.0 v))))))
(defn- dual-digamma [x] (r/chain [v x, _y (digamma v)] :digamma (trigamma v)))
(defn- dual-trigamma [x] (r/chain [v x, _y (trigamma v)] :trigamma (polygamma 2.0 v)))
(defn- dual-polygamma [m x] (r/chain [m m, x x, _y (polygamma m x)] :polygamma nil (polygamma (add m 1.0) x)))

(defn- dual-beta
  [p q]
  (r/chain [p p, q q, y (beta p q)] :beta
    (mul y (sub (digamma p) (digamma (add p q))))
    (mul y (sub (digamma q) (digamma (add p q))))))

(defn- dual-log-beta
  [p q]
  (r/chain [p p, q q, _y (log-beta p q)] :log-beta
    (sub (digamma p) (digamma (add p q)))
    (sub (digamma q) (digamma (add p q)))))

;; ## Derivative rules: incomplete functions, in `x` only
;;
;; The densities are computed in log space, stable for large shape parameters (`x^(a-1)` and `gamma(a)` overflow
;; above 171). At the ends of the domain (`x = 0`, and `x = 1` for beta) the log-space form is `0 * -Inf = NaN` for a
;; shape parameter equal to `1`, the power form is used there. The normalizations `log-gamma` and `log-beta` are the
;; logarithms of the absolute values, their signs (negative for some negative shape parameters) are applied
;; separately.

(defn- gamma-sign
  "Sign of `gamma(a)` for a non-pole `a`: `-1` for a negative `a` with an odd `floor(a)`, otherwise `1`. The parity is
  computed in doubles: `a` can be beyond the `long` range (every such double is an even integer)."
  ^double [a]
  (let [sa (d/scalar a)]
    (if (m/neg? sa)
      (let [f (m/floor sa)]
        (if (m/== 1.0 (m/abs (m/- f (m/* 2.0 (m/floor (m/* 0.5 f)))))) -1.0 1.0))
      1.0)))

(defn- beta-sign
  "Sign of `beta(a, b) = gamma(a) gamma(b) / gamma(a + b)`."
  ^double [a b]
  (m/* (gamma-sign a) (gamma-sign b) (gamma-sign (add a b))))

(defn- gamma-density
  "`x^(s-1) e^(-x) / exp(log-norm)`; `0` at `x = +Inf`, its limit (the log-space form is `Inf - Inf = NaN` there for
  `s > 1`)."
  [s x log-norm]
  (let [sx (d/scalar x)]
    (cond
      (m/zero? sx) (div (mul (pow x (sub s 1.0)) (exp (neg x))) (exp log-norm))
      (m/== sx ##Inf) 0.0
      :else (exp (sub (sub (mul (sub s 1.0) (log x)) x) log-norm)))))

(defn- scale-keeping-zeros
  "The number or dual `p` times the double `s`, every zero of the value and of the partials (recursively) kept zero,
  also for an infinite `s`. Used where the zeros are exact (positive powers of zero) and `s` is the reciprocal of an
  underflowed normalization: the true product of an exact zero and a finite, unrepresentable number is zero."
  [p ^double s]
  (if (instance? Dual p)
    (let [^Dual p p]
      (d/new-dual (.tag p) (scale-keeping-zeros (.value p) s) (d/p-map #(scale-keeping-zeros % s) (.partials p))))
    (let [v (.doubleValue ^Number p)]
      (if (m/zero? v) v (m/* v s)))))

(defn- beta-density
  "`x^(a-1) (1-x)^(b-1) / exp(log-norm)`.

  At the ends the power form is used. When the normalization underflows (`beta(a, b)` below the smallest double) its
  reciprocal is infinite and multiplies the power form keeping the zeros: they are exact zeros of positive powers
  (`a - 1` at `x = 0`, `b - 1` at `x = 1`), so `0 * Inf` is `0` and every order is exact (`0`, or `Inf` where the
  true value overflows)."
  [x a b log-norm]
  (let [sx (d/scalar x)
        ea (sub a 1.0)
        eb (sub b 1.0)]
    (if (or (m/zero? sx) (m/== 1.0 sx))
      (let [norm (exp log-norm)
            power-form (mul (pow x ea) (pow (sub 1.0 x) eb))]
        (if (m/zero? (d/scalar norm))
          (scale-keeping-zeros power-form (m/exp (m/- (d/scalar log-norm))))
          (div power-form norm)))
      (exp (sub (add (mul ea (log x)) (mul eb (el/log1p (neg x)))) log-norm)))))

(defn- dual-lower-incomplete-gamma
  [s x]
  (r/chain [s s, x x, _y (lower-incomplete-gamma s x)] :lower-incomplete-gamma nil (gamma-density s x 0.0)))

(defn- dual-upper-incomplete-gamma
  [s x]
  (r/chain [s s, x x, _y (upper-incomplete-gamma s x)] :upper-incomplete-gamma nil (neg (gamma-density s x 0.0))))

(defn- dual-regularized-gamma-p
  [a x]
  (r/chain [a a, x x, _y (regularized-gamma-p a x)] :regularized-gamma-p nil (mul (gamma-sign a) (gamma-density a x (log-gamma a)))))

(defn- dual-regularized-gamma-q
  [a x]
  (r/chain [a a, x x, _y (regularized-gamma-q a x)] :regularized-gamma-q nil (mul (m/- (gamma-sign a)) (gamma-density a x (log-gamma a)))))

(defn- dual-incomplete-beta
  [x a b]
  (r/chain [x x, a a, b b, _y (incomplete-beta x a b)] :incomplete-beta (beta-density x a b 0.0) nil nil))

(defn- dual-regularized-beta
  [x a b]
  (r/chain [x x, a a, b b, _y (regularized-beta x a b)] :regularized-beta (mul (beta-sign a b) (beta-density x a b (log-beta a b))) nil nil))

;; ## Derivative rules: the error functions and Lambert W

(def ^{:private true :const true :tag 'double} sqrt-pi-over-two (m// m/SQRTPI 2.0))

(defn- erf-derivative "`(2/sqrt(pi)) exp(-x^2)`." [x] (mul m/M_2_SQRTPI (exp (neg (mul x x)))))

(defn- dual-erf [x] (r/chain [v x, _y (erf v)] :erf (erf-derivative v)))
(defn- dual-erf2 [x1 x2] (r/chain [a x1, b x2, _y (erf a b)] :erf (neg (erf-derivative a)) (erf-derivative b)))
(defn- dual-erfc [x] (r/chain [v x, _y (erfc v)] :erfc (neg (erf-derivative v))))
(defn- dual-inv-erf [x] (r/chain [v x, y (inv-erf v)] :inv-erf (mul sqrt-pi-over-two (exp (mul y y)))))
(defn- dual-inv-erfc [x] (r/chain [v x, y (inv-erfc v)] :inv-erfc (neg (mul sqrt-pi-over-two (exp (mul y y))))))

(defn- lambert-W-derivative
  "The derivative of either branch at `x` with the value `W`: `1 / (exp(W) (1 + W))` for `|W| <= 1` (`1` at `0`, where
  the quotient form is `0/0`; infinite at the branch point `-1/e`), otherwise `W / (x (1 + W))` (`exp(W)` would amplify
  the rounding of `W` about `|W|` times: `3e-14` at `W = -697`). The limits at the infinite values: `-Inf` for
  `W = -Inf` (the lower branch at `0`), `0` for `W = Inf` (the quotient form is `Inf/Inf` there); the higher derivatives
  at these points are not limits."
  [x w]
  (let [sw (d/scalar w)]
    (cond
      (m/== ##-Inf sw) ##-Inf
      (m/== ##Inf sw) 0.0
      (m/<= (m/abs sw) 1.0) (div 1.0 (mul (exp w) (add 1.0 w)))
      :else (div w (mul x (add 1.0 w))))))

(defn- dual-lambert-W [x] (r/chain [v x, w (lambert-W v)] :lambert-W (lambert-W-derivative v w)))
(defn- dual-lambert-W-1 [x] (r/chain [v x, w (lambert-W-1 v)] :lambert-W-1 (lambert-W-derivative v w)))
