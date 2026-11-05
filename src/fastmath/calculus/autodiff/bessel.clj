(ns ^:no-doc fastmath.calculus.autodiff.bessel
  "Internal: the Bessel functions (first and second kind, modified, half odd order), the spherical Bessel functions and the Airy functions of `fastmath.special` for numbers and duals of [[fastmath.calculus.autodiff]], dispatched by `instance?`, and their derivative rules: the recurrences through the neighbouring orders (DLMF 10.6.2, 10.29.2, 10.51.2, 10.51.5), without a division by `x`, and `Ai'' = x Ai`, `Bi'' = x Bi`. The orders have no derivative."
  (:require [fastmath.core :as m]
            [fastmath.special :as s]
            [fastmath.calculus.autodiff.dual :as d]
            [fastmath.calculus.autodiff.ops :refer [add sub mul div neg]]
            [fastmath.calculus.autodiff.rules :as r])
  (:import [fastmath.calculus.autodiff.dual Dual]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(declare dual-bessel-J dual-bessel-Y dual-bessel-I dual-bessel-K dual-bessel-K-half-odd dual-bessel-K-half-odd-scaled
         dual-spherical-j dual-spherical-y dual-spherical-1-i dual-spherical-2-i dual-spherical-k
         dual-airy-Ai dual-airy-Ai' dual-airy-Bi dual-airy-Bi')

;; ## Operations on numbers and duals
;;
;; The functions of a fixed order are the functions of the general order with a constant order: the derivative
;; rules are those of the general order, the values those of the fixed-order functions of `fastmath.special`.

(defn bessel-J "Bessel function of the first kind of order `order`." [order x] (r/dispatch dual-bessel-J s/bessel-J order x))
(defn bessel-Y "Bessel function of the second kind of order `order`." [order x] (r/dispatch dual-bessel-Y s/bessel-Y order x))
(defn bessel-I "Modified Bessel function of the first kind of order `order`." [order x] (r/dispatch dual-bessel-I s/bessel-I order x))
(defn bessel-K "Modified Bessel function of the second kind of order `order`." [order x] (r/dispatch dual-bessel-K s/bessel-K order x))
(defn bessel-K-half-odd "`K` of the order `order/2`, `order` odd." [order x] (r/dispatch dual-bessel-K-half-odd s/bessel-K-half-odd order x))
(defn bessel-K-half-odd-scaled "`exp(x) K` of the order `order/2`, `order` odd." [order x] (r/dispatch dual-bessel-K-half-odd-scaled s/bessel-K-half-odd-scaled order x))
(defn spherical-bessel-j "Spherical Bessel function of the first kind." [order x] (r/dispatch dual-spherical-j s/spherical-bessel-j order x))
(defn spherical-bessel-y "Spherical Bessel function of the second kind." [order x] (r/dispatch dual-spherical-y s/spherical-bessel-y order x))
(defn spherical-bessel-1-i "First modified spherical Bessel function of the first kind." [order x] (r/dispatch dual-spherical-1-i s/spherical-bessel-1-i order x))
(defn spherical-bessel-2-i "Second modified spherical Bessel function of the first kind." [order x] (r/dispatch dual-spherical-2-i s/spherical-bessel-2-i order x))
(defn spherical-bessel-k "Modified spherical Bessel function of the second kind." [order x] (r/dispatch dual-spherical-k s/spherical-bessel-k order x))
(defn airy-Ai "Airy function `Ai`." [x] (r/dispatch dual-airy-Ai s/airy-Ai x))
(defn airy-Ai' "Derivative of the Airy function `Ai`." [x] (r/dispatch dual-airy-Ai' s/airy-Ai' x))
(defn airy-Bi "Airy function `Bi`." [x] (r/dispatch dual-airy-Bi s/airy-Bi x))
(defn airy-Bi' "Derivative of the Airy function `Bi`." [x] (r/dispatch dual-airy-Bi' s/airy-Bi' x))

;; ## Derivative rules: Bessel functions
;;
;; Each rule takes the value function `(value-fn n x)` as its last argument: the general-order function of this
;; namespace, or a fixed-order function (the order is then a constant equal to the fixed order).

(defn- bessel-J-rule
  [order x value-fn]
  (r/chain [order order, x x, _y (value-fn order x)] :bessel-J nil (mul 0.5 (sub (bessel-J (sub order 1.0) x) (bessel-J (add order 1.0) x)))))

(defn- rising-from-minus-infinity?
  "Is the value `y` at `x = 0` equal to `-Inf` (the functions of the second kind)? The derivative there is the limit
  `+Inf` (the recurrences are `-Inf - (-Inf)`); the higher derivatives at that point are not limits."
  [x y]
  (and (m/zero? (d/scalar x)) (m/== ##-Inf (d/scalar y))))

(defn- bessel-Y-rule
  [order x value-fn]
  (r/chain [order order, x x, y (value-fn order x)] :bessel-Y
    nil (if (rising-from-minus-infinity? x y)
          ##Inf
          (mul 0.5 (sub (bessel-Y (sub order 1.0) x) (bessel-Y (add order 1.0) x))))))

(defn- bessel-I-rule
  [order x value-fn]
  (r/chain [order order, x x, _y (value-fn order x)] :bessel-I nil (mul 0.5 (add (bessel-I (sub order 1.0) x) (bessel-I (add order 1.0) x)))))

(defn- bessel-K-rule
  [order x value-fn]
  (r/chain [order order, x x, _y (value-fn order x)] :bessel-K nil (neg (mul 0.5 (add (bessel-K (sub order 1.0) x) (bessel-K (add order 1.0) x))))))

(defn- dual-bessel-J [order x] (bessel-J-rule order x bessel-J))
(defn- dual-bessel-Y [order x] (bessel-Y-rule order x bessel-Y))
(defn- dual-bessel-I [order x] (bessel-I-rule order x bessel-I))
(defn- dual-bessel-K [order x] (bessel-K-rule order x bessel-K))

(defn- dual-bessel-K-half-odd
  "`K_a' = -(K_(a-1) + K_(a+1))/2` with `a = order/2`: the orders `order - 2` and `order + 2`, odd again."
  [order x]
  (r/chain [order order, x x, _y (bessel-K-half-odd order x)] :bessel-K-half-odd
    nil (neg (mul 0.5 (add (bessel-K-half-odd (sub order 2.0) x) (bessel-K-half-odd (add order 2.0) x))))))

(defn- dual-bessel-K-half-odd-scaled
  "`(exp(x) K_a)' = exp(x) K_a - (exp(x) K_(a-1) + exp(x) K_(a+1))/2`."
  [order x]
  (r/chain [order order, x x, y (bessel-K-half-odd-scaled order x)] :bessel-K-half-odd-scaled
    nil (sub y (mul 0.5 (add (bessel-K-half-odd-scaled (sub order 2.0) x) (bessel-K-half-odd-scaled (add order 2.0) x))))))

;; ## Derivative rules: spherical Bessel functions
;;
;; `(2n + 1) f_n' = n f_(n-1) -/+ (n + 1) f_(n+1)` (DLMF 10.51.2 for `j`, `y`; 10.51.5 for `i1`, `i2` and `(-1)^n k`),
;; finite at `x = 0` (`j_1'(0) = 1/3`). For the orders in `(-1, 0)` it cancels: `2n + 1` vanishes at `-1/2` and the
;; relative error grows like `eps / |2n + 1|` (3e-3 at `n = -1/2 + 1e-12`, 1e-13 at `-0.53`). There
;; `f_n' = f_(n-1) - ((n + 1)/x) f_n` (10.51.1, 10.51.4) is used, which is no less accurate on all of `(-1, 0)` (measured
;; against mpmath at `x` in `[0.3, 20]`: its error follows the error of the values); it divides by `x`, where the
;; functions of these orders are singular. The orders `n >= 0` keep the recurrence, which is finite at `x = 0`.

(defn- order-between-minus-one-and-zero?
  "Is the order `n` in `(-1, 0)`, where the recurrence through both neighbouring orders cancels?"
  [n]
  (m/< -1.0 (d/scalar n) 0.0))

(defn- lower-order-term
  "`n f_(n-1)(x)`, `0` for `n = 0` (where `f_(-1)` can be infinite: `0 * Inf` at `x = 0`)."
  [f n x]
  (if (m/zero? (d/scalar n)) 0.0 (mul n (f (sub n 1.0) x))))

(defn- upper-order-term
  "`(n + 1) f_(n+1)(x)`, `0` for `n = -1` (where `f_0` can be infinite: `0 * Inf` at `x = 0`, e.g. `y_(-1)'(0) = 0`)."
  [f n x]
  (if (m/== -1.0 (d/scalar n)) 0.0 (mul (add n 1.0) (f (add n 1.0) x))))

(defn- spherical-derivative
  "The derivative of `f` of the order `n` at `x` with the value `y`. `sign` is `-1` for `j`, `y` (`f_(n+1)` subtracted),
  `1` for `i1`, `i2`; `k` is handled by the caller through `(-1)^n k`."
  [f n x y sign]
  (if (order-between-minus-one-and-zero? n)
    (sub (f (sub n 1.0) x) (mul (div (add n 1.0) x) y))
    (div (add (lower-order-term f n x) (mul sign (upper-order-term f n x)))
         (add (mul 2.0 n) 1.0))))

(defn- spherical-j-rule [order x value-fn] (r/chain [order order, x x, y (value-fn order x)] :spherical-bessel-j nil (spherical-derivative spherical-bessel-j order x y -1.0)))
(defn- spherical-y-rule
  [order x value-fn]
  (r/chain [order order, x x, y (value-fn order x)] :spherical-bessel-y
    nil (if (rising-from-minus-infinity? x y) ##Inf (spherical-derivative spherical-bessel-y order x y -1.0))))
(defn- spherical-1-i-rule [order x value-fn] (r/chain [order order, x x, y (value-fn order x)] :spherical-bessel-1-i nil (spherical-derivative spherical-bessel-1-i order x y 1.0)))
(defn- spherical-2-i-rule [order x value-fn] (r/chain [order order, x x, y (value-fn order x)] :spherical-bessel-2-i nil (spherical-derivative spherical-bessel-2-i order x y 1.0)))

(defn- spherical-k-rule
  "`(2n + 1) k_n' = -n k_(n-1) - (n + 1) k_(n+1)`; for the orders in `(-1, 0)` `k_n' = -k_(n-1) - ((n + 1)/x) k_n`."
  [order x value-fn]
  (r/chain [order order, x x, y (value-fn order x)] :spherical-bessel-k
    nil (if (order-between-minus-one-and-zero? order)
          (sub (neg (spherical-bessel-k (sub order 1.0) x)) (mul (div (add order 1.0) x) y))
          (neg (div (add (lower-order-term spherical-bessel-k order x) (upper-order-term spherical-bessel-k order x))
                    (add (mul 2.0 order) 1.0))))))

(defn- dual-spherical-j [order x] (spherical-j-rule order x spherical-bessel-j))
(defn- dual-spherical-y [order x] (spherical-y-rule order x spherical-bessel-y))
(defn- dual-spherical-1-i [order x] (spherical-1-i-rule order x spherical-bessel-1-i))
(defn- dual-spherical-2-i [order x] (spherical-2-i-rule order x spherical-bessel-2-i))
(defn- dual-spherical-k [order x] (spherical-k-rule order x spherical-bessel-k))

;; ## Functions of a fixed order

(defn- fixed
  "The function of a fixed order `n`: `double-fn` (of `fastmath.special`) for a number, the rule of the general order
  for a dual, with this function as its value function (for the duals of the older layers in the value)."
  [rule ^double n double-fn]
  (letfn [(f [x] (if (instance? Dual x) (rule n x (fn [_ v] (f v))) (double-fn (.doubleValue ^Number x))))]
    f))

(def ^{:arglists '([x]) :doc "Bessel function of the first kind of order 0."} bessel-J0 (fixed bessel-J-rule 0.0 s/bessel-J0))
(def ^{:arglists '([x]) :doc "Bessel function of the first kind of order 1."} bessel-J1 (fixed bessel-J-rule 1.0 s/bessel-J1))
(def ^{:arglists '([x]) :doc "Bessel function of the second kind of order 0."} bessel-Y0 (fixed bessel-Y-rule 0.0 s/bessel-Y0))
(def ^{:arglists '([x]) :doc "Bessel function of the second kind of order 1."} bessel-Y1 (fixed bessel-Y-rule 1.0 s/bessel-Y1))
(def ^{:arglists '([x]) :doc "Modified Bessel function of the first kind of order 0."} bessel-I0 (fixed bessel-I-rule 0.0 s/bessel-I0))
(def ^{:arglists '([x]) :doc "Modified Bessel function of the first kind of order 1."} bessel-I1 (fixed bessel-I-rule 1.0 s/bessel-I1))
(def ^{:arglists '([x]) :doc "Modified Bessel function of the second kind of order 0."} bessel-K0 (fixed bessel-K-rule 0.0 s/bessel-K0))
(def ^{:arglists '([x]) :doc "Modified Bessel function of the second kind of order 1."} bessel-K1 (fixed bessel-K-rule 1.0 s/bessel-K1))
(def ^{:arglists '([x]) :doc "Spherical Bessel function of the first kind of order 0."} spherical-bessel-j0 (fixed spherical-j-rule 0.0 s/spherical-bessel-j0))
(def ^{:arglists '([x]) :doc "Spherical Bessel function of the first kind of order 1."} spherical-bessel-j1 (fixed spherical-j-rule 1.0 s/spherical-bessel-j1))
(def ^{:arglists '([x]) :doc "Spherical Bessel function of the first kind of order 2."} spherical-bessel-j2 (fixed spherical-j-rule 2.0 s/spherical-bessel-j2))
(def ^{:arglists '([x]) :doc "Spherical Bessel function of the second kind of order 0."} spherical-bessel-y0 (fixed spherical-y-rule 0.0 s/spherical-bessel-y0))
(def ^{:arglists '([x]) :doc "Spherical Bessel function of the second kind of order 1."} spherical-bessel-y1 (fixed spherical-y-rule 1.0 s/spherical-bessel-y1))
(def ^{:arglists '([x]) :doc "Spherical Bessel function of the second kind of order 2."} spherical-bessel-y2 (fixed spherical-y-rule 2.0 s/spherical-bessel-y2))
(def ^{:arglists '([x]) :doc "First modified spherical Bessel function of the first kind of order 0."} spherical-bessel-1-i0 (fixed spherical-1-i-rule 0.0 s/spherical-bessel-1-i0))
(def ^{:arglists '([x]) :doc "First modified spherical Bessel function of the first kind of order 1."} spherical-bessel-1-i1 (fixed spherical-1-i-rule 1.0 s/spherical-bessel-1-i1))
(def ^{:arglists '([x]) :doc "First modified spherical Bessel function of the first kind of order 2."} spherical-bessel-1-i2 (fixed spherical-1-i-rule 2.0 s/spherical-bessel-1-i2))
(def ^{:arglists '([x]) :doc "Second modified spherical Bessel function of the first kind of order 0."} spherical-bessel-2-i0 (fixed spherical-2-i-rule 0.0 s/spherical-bessel-2-i0))
(def ^{:arglists '([x]) :doc "Second modified spherical Bessel function of the first kind of order 1."} spherical-bessel-2-i1 (fixed spherical-2-i-rule 1.0 s/spherical-bessel-2-i1))
(def ^{:arglists '([x]) :doc "Second modified spherical Bessel function of the first kind of order 2."} spherical-bessel-2-i2 (fixed spherical-2-i-rule 2.0 s/spherical-bessel-2-i2))
(def ^{:arglists '([x]) :doc "Modified spherical Bessel function of the second kind of order 0."} spherical-bessel-k0 (fixed spherical-k-rule 0.0 s/spherical-bessel-k0))
(def ^{:arglists '([x]) :doc "Modified spherical Bessel function of the second kind of order 1."} spherical-bessel-k1 (fixed spherical-k-rule 1.0 s/spherical-bessel-k1))
(def ^{:arglists '([x]) :doc "Modified spherical Bessel function of the second kind of order 2."} spherical-bessel-k2 (fixed spherical-k-rule 2.0 s/spherical-bessel-k2))

;; ## Derivative rules: Airy functions

(defn- dual-airy-Ai [x] (r/chain [v x, _y (airy-Ai v)] :airy-Ai (airy-Ai' v)))
(defn- dual-airy-Ai' [x] (r/chain [v x, _y (airy-Ai' v)] :airy-Ai' (mul v (airy-Ai v))))
(defn- dual-airy-Bi [x] (r/chain [v x, _y (airy-Bi v)] :airy-Bi (airy-Bi' v)))
(defn- dual-airy-Bi' [x] (r/chain [v x, _y (airy-Bi' v)] :airy-Bi' (mul v (airy-Bi v))))
