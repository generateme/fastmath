(ns ^:no-doc fastmath.calculus.autodiff.special
  "Internal: the gamma and beta functions of numbers and duals of [[fastmath.calculus.autodiff]], dispatched by `instance?`, and their derivative rules. The shape parameters of the incomplete functions and the order of [[polygamma]] have no derivative."
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
         dual-regularized-gamma-q dual-incomplete-beta dual-regularized-beta)

;; ## Operations on numbers and duals

(defmacro ^:private dispatch
  "A dual in any argument selects the dual rule, otherwise the arguments are coerced to doubles and passed to the
  `fastmath.special` function."
  [dual-rule special-fn & args]
  `(if (or ~@(map (fn [a] `(instance? Dual ~a)) args))
     (~dual-rule ~@args)
     (~special-fn ~@(map (fn [a] `(.doubleValue ~(with-meta a {:tag 'Number}))) args))))

(defn gamma "Gamma function." [x] (dispatch dual-gamma s/gamma x))
(defn log-gamma "Logarithm of the absolute value of the gamma function." [x] (dispatch dual-log-gamma s/log-gamma x))
(defn log-gamma-1p "`log-gamma(1 + x)`." [x] (dispatch dual-log-gamma-1p s/log-gamma-1p x))
(defn inv-gamma-1pm1 "`1 / gamma(1 + x) - 1`." [x] (dispatch dual-inv-gamma-1pm1 s/inv-gamma-1pm1 x))
(defn digamma "Digamma function." [x] (dispatch dual-digamma s/digamma x))
(defn trigamma "Trigamma function." [x] (dispatch dual-trigamma s/trigamma x))
(defn polygamma "Polygamma function of order `m`." [m x] (dispatch dual-polygamma s/polygamma m x))
(defn beta "Beta function." [p q] (dispatch dual-beta s/beta p q))
(defn log-beta "Logarithm of the beta function." [p q] (dispatch dual-log-beta s/log-beta p q))
(defn lower-incomplete-gamma "Lower incomplete gamma function." [s x] (dispatch dual-lower-incomplete-gamma s/lower-incomplete-gamma s x))
(defn upper-incomplete-gamma "Upper incomplete gamma function." [s x] (dispatch dual-upper-incomplete-gamma s/upper-incomplete-gamma s x))
(defn regularized-gamma-p "Regularized lower incomplete gamma function." [a x] (dispatch dual-regularized-gamma-p s/regularized-gamma-p a x))
(defn regularized-gamma-q "Regularized upper incomplete gamma function." [a x] (dispatch dual-regularized-gamma-q s/regularized-gamma-q a x))
(defn incomplete-beta "Incomplete beta function, `x` first." [x a b] (dispatch dual-incomplete-beta s/incomplete-beta x a b))
(defn regularized-beta "Regularized incomplete beta function, `x` first." [x a b] (dispatch dual-regularized-beta s/regularized-beta x a b))

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
