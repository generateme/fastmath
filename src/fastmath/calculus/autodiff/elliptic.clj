(ns ^:no-doc fastmath.calculus.autodiff.elliptic
  "Internal: the Jacobi elliptic functions, their inverses, the amplitude and the incomplete elliptic integrals of `fastmath.special` for numbers and duals of [[fastmath.calculus.autodiff]], dispatched by `instance?`, and their derivative rules in `u` (Jacobi functions, DLMF 22.13), in `x` (inverses, DLMF 22.15 with the signs of the principal branches of `fastmath.special`) and in the amplitude `phi` (integrals, the integrands). The parameter `m` (`m = k^2`) and the characteristic `n` have no derivative."
  (:require [fastmath.core :as m]
            [fastmath.special :as s]
            [fastmath.calculus.autodiff.dual :as d]
            [fastmath.calculus.autodiff.ops :refer [add sub mul div neg sin]]
            [fastmath.calculus.autodiff.elementary :as el]
            [fastmath.calculus.autodiff.rules :as r])
  (:import [fastmath.calculus.autodiff.dual Dual]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(declare jacobi-sn jacobi-cn jacobi-dn jacobi-rule inverse-rule inverse-derivatives
         dual-jacobi-am dual-elliptic-F dual-elliptic-E dual-elliptic-D dual-elliptic-PI)

;; ## Jacobi elliptic functions
;;
;; The derivatives in `u` (DLMF 22.13.1-22.13.12) in terms of `sn`, `cn`, `dn` at `(u, m)`, written with the functions
;; of this namespace, so they nest; `m' = 1 - m`.

(defn- complement-m [m] (sub 1.0 m))

(def ^:private jacobi-derivatives
  "Function name -> the derivative in `u` as a function of `[sn cn dn m]`, `sn`, `cn`, `dn` delays (only the needed
  ones are evaluated)."
  {:sn (fn [_ cn dn _] (mul @cn @dn))
   :cn (fn [sn _ dn _] (neg (mul @sn @dn)))
   :dn (fn [sn cn _ m] (neg (mul m (mul @sn @cn))))
   :sc (fn [_ cn dn _] (div @dn (mul @cn @cn)))
   :sd (fn [_ cn dn _] (div @cn (mul @dn @dn)))
   :cs (fn [sn _ dn _] (neg (div @dn (mul @sn @sn))))
   :cd (fn [sn _ dn m] (neg (div (mul (complement-m m) @sn) (mul @dn @dn))))
   :ds (fn [sn cn _ _] (neg (div @cn (mul @sn @sn))))
   :dc (fn [sn cn _ m] (div (mul (complement-m m) @sn) (mul @cn @cn)))
   :ns (fn [sn cn dn _] (neg (div (mul @cn @dn) (mul @sn @sn))))
   :nc (fn [sn cn dn _] (div (mul @sn @dn) (mul @cn @cn)))
   :nd (fn [sn cn dn m] (div (mul m (mul @sn @cn)) (mul @dn @dn)))})

(defn- jacobi-dispatch
  "`(double-fn u m)` for numbers, the rule of the function `id` for a dual."
  [id f double-fn u m]
  (if (or (instance? Dual u) (instance? Dual m))
    (jacobi-rule id f u m)
    (double-fn (.doubleValue ^Number u) (.doubleValue ^Number m))))

(defn jacobi-sn "Jacobi elliptic function `sn(u|m)`." [u m] (jacobi-dispatch :sn jacobi-sn s/jacobi-sn u m))
(defn jacobi-cn "Jacobi elliptic function `cn(u|m)`." [u m] (jacobi-dispatch :cn jacobi-cn s/jacobi-cn u m))
(defn jacobi-dn "Jacobi elliptic function `dn(u|m)`." [u m] (jacobi-dispatch :dn jacobi-dn s/jacobi-dn u m))
(defn jacobi-sc "Jacobi elliptic function `sc(u|m)`." [u m] (jacobi-dispatch :sc jacobi-sc s/jacobi-sc u m))
(defn jacobi-sd "Jacobi elliptic function `sd(u|m)`." [u m] (jacobi-dispatch :sd jacobi-sd s/jacobi-sd u m))
(defn jacobi-cs "Jacobi elliptic function `cs(u|m)`." [u m] (jacobi-dispatch :cs jacobi-cs s/jacobi-cs u m))
(defn jacobi-cd "Jacobi elliptic function `cd(u|m)`." [u m] (jacobi-dispatch :cd jacobi-cd s/jacobi-cd u m))
(defn jacobi-ds "Jacobi elliptic function `ds(u|m)`." [u m] (jacobi-dispatch :ds jacobi-ds s/jacobi-ds u m))
(defn jacobi-dc "Jacobi elliptic function `dc(u|m)`." [u m] (jacobi-dispatch :dc jacobi-dc s/jacobi-dc u m))
(defn jacobi-ns "Jacobi elliptic function `ns(u|m)`." [u m] (jacobi-dispatch :ns jacobi-ns s/jacobi-ns u m))
(defn jacobi-nc "Jacobi elliptic function `nc(u|m)`." [u m] (jacobi-dispatch :nc jacobi-nc s/jacobi-nc u m))
(defn jacobi-nd "Jacobi elliptic function `nd(u|m)`." [u m] (jacobi-dispatch :nd jacobi-nd s/jacobi-nd u m))

(defn- jacobi-derivative
  "The derivative in `u` of the Jacobi function `id` at `(u, m)`."
  [id u m]
  ((jacobi-derivatives id) (delay (jacobi-sn u m)) (delay (jacobi-cn u m)) (delay (jacobi-dn u m)) m))

(defn- jacobi-rule
  [id f u m]
  (r/chain [u u, m m, _y (f u m)] (keyword (str "jacobi-" (name id))) (jacobi-derivative id u m) nil))

;; ## Amplitude

(defn jacobi-am "Jacobi amplitude `am(u|m)`." [u m] (r/dispatch dual-jacobi-am s/jacobi-am u m))
(defn- dual-jacobi-am [u m] (r/chain [u u, m m, _y (jacobi-am u m)] :jacobi-am (jacobi-dn u m) nil))

;; ## Inverses
;;
;; The derivatives in `x` (DLMF 22.15.12-22.15.23, the integrands of the integral representations of the principal
;; values), with `m' = 1 - m`. They are exact at the ends of the range (infinite) and `NaN` outside it (a negative
;; radicand); `1/f'(u)` at `u = f^-1(x)` has the same magnitude inside the range, but at its ends `f'(u)` is a rounding
;; error (`cn(K)` is about `1e-16`, with any sign). The magnitudes hold for every `m` (they follow from identities as
;; `dn^2 = 1 - m sn^2`); DLMF states the signs for `0 <= m <= 1`, and the principal branches of `fastmath.special` have
;; the opposite sign for `cd`, `dc` when `m > 1` and for `dn`, `nd` when `m < 0` (measured against `1/f'(u)` on
;; `x` in `[-4, 4]` and `m` in `{-3, -2, -0.5, 0.3, 0.8, 1.2, 1.5, 3}`: the sign is constant in each region).

(defn- x2 [x] (mul x x))

(defn- inverse-derivative
  "`sign / sqrt(a * b)`."
  [sign a b]
  (div sign (el/sqrt (mul a b))))

(def ^:private inverse-derivatives
  "Function name -> the derivative of its inverse in `x` as a function of `[x m]`."
  {:sn (fn [x m] (inverse-derivative 1.0 (sub 1.0 (x2 x)) (sub 1.0 (mul m (x2 x)))))
   :cn (fn [x m] (inverse-derivative -1.0 (sub 1.0 (x2 x)) (add (complement-m m) (mul m (x2 x)))))
   :dn (fn [x m] (inverse-derivative -1.0 (sub 1.0 (x2 x)) (sub (x2 x) (complement-m m))))
   :sc (fn [x m] (inverse-derivative 1.0 (add 1.0 (x2 x)) (add 1.0 (mul (complement-m m) (x2 x)))))
   :sd (fn [x m] (inverse-derivative 1.0 (sub 1.0 (mul (complement-m m) (x2 x))) (add 1.0 (mul m (x2 x)))))
   :cs (fn [x m] (inverse-derivative -1.0 (add (x2 x) 1.0) (add (x2 x) (complement-m m))))
   :cd (fn [x m] (inverse-derivative -1.0 (sub 1.0 (x2 x)) (sub 1.0 (mul m (x2 x)))))
   :ds (fn [x m] (inverse-derivative -1.0 (sub (x2 x) (complement-m m)) (add (x2 x) m)))
   :dc (fn [x m] (inverse-derivative 1.0 (sub (x2 x) 1.0) (sub (x2 x) m)))
   :ns (fn [x m] (inverse-derivative -1.0 (sub (x2 x) 1.0) (sub (x2 x) m)))
   :nc (fn [x m] (inverse-derivative 1.0 (sub (x2 x) 1.0) (add (mul (complement-m m) (x2 x)) m)))
   :nd (fn [x m] (inverse-derivative 1.0 (sub (x2 x) 1.0) (sub 1.0 (mul (complement-m m) (x2 x)))))})

(defn- inverse-dispatch
  "`(double-fn x m)` for numbers, the rule of the inverse of `id` for a dual."
  [id f double-fn x m]
  (if (or (instance? Dual x) (instance? Dual m))
    (inverse-rule id f x m)
    (double-fn (.doubleValue ^Number x) (.doubleValue ^Number m))))

(defn jacobi-asn "Inverse of the Jacobi elliptic function `sn`: `u` with `sn(u|m) = x`." [x m] (inverse-dispatch :sn jacobi-asn s/jacobi-asn x m))
(defn jacobi-acn "Inverse of the Jacobi elliptic function `cn`: `u` with `cn(u|m) = x`." [x m] (inverse-dispatch :cn jacobi-acn s/jacobi-acn x m))
(defn jacobi-adn "Inverse of the Jacobi elliptic function `dn`: `u` with `dn(u|m) = x`." [x m] (inverse-dispatch :dn jacobi-adn s/jacobi-adn x m))
(defn jacobi-asc "Inverse of the Jacobi elliptic function `sc`: `u` with `sc(u|m) = x`." [x m] (inverse-dispatch :sc jacobi-asc s/jacobi-asc x m))
(defn jacobi-asd "Inverse of the Jacobi elliptic function `sd`: `u` with `sd(u|m) = x`." [x m] (inverse-dispatch :sd jacobi-asd s/jacobi-asd x m))
(defn jacobi-acs "Inverse of the Jacobi elliptic function `cs`: `u` with `cs(u|m) = x`." [x m] (inverse-dispatch :cs jacobi-acs s/jacobi-acs x m))
(defn jacobi-acd "Inverse of the Jacobi elliptic function `cd`: `u` with `cd(u|m) = x`." [x m] (inverse-dispatch :cd jacobi-acd s/jacobi-acd x m))
(defn jacobi-ads "Inverse of the Jacobi elliptic function `ds`: `u` with `ds(u|m) = x`." [x m] (inverse-dispatch :ds jacobi-ads s/jacobi-ads x m))
(defn jacobi-adc "Inverse of the Jacobi elliptic function `dc`: `u` with `dc(u|m) = x`." [x m] (inverse-dispatch :dc jacobi-adc s/jacobi-adc x m))
(defn jacobi-ans "Inverse of the Jacobi elliptic function `ns`: `u` with `ns(u|m) = x`." [x m] (inverse-dispatch :ns jacobi-ans s/jacobi-ans x m))
(defn jacobi-anc "Inverse of the Jacobi elliptic function `nc`: `u` with `nc(u|m) = x`." [x m] (inverse-dispatch :nc jacobi-anc s/jacobi-anc x m))
(defn jacobi-and "Inverse of the Jacobi elliptic function `nd`: `u` with `nd(u|m) = x`." [x m] (inverse-dispatch :nd jacobi-and s/jacobi-and x m))

(defn- branch-sign
  "`-1` where the principal branch of `fastmath.special` has the opposite sign to the DLMF form: `cd`, `dc` for
  `m > 1`, `dn`, `nd` for `m < 0`; otherwise `1`."
  ^double [id m]
  (let [sm (d/scalar m)]
    (cond
      (and (#{:cd :dc} id) (m/> sm 1.0)) -1.0
      (and (#{:dn :nd} id) (m/neg? sm)) -1.0
      :else 1.0)))

(defn- inverse-rule
  [id f x m]
  (r/chain [x x, m m, _u (f x m)] (keyword (str "jacobi-a" (name id)))
    (mul (branch-sign id m) ((inverse-derivatives id) x m))
    nil))

;; ## Incomplete elliptic integrals: the derivatives in `phi` are the integrands, `Delta = sqrt(1 - m sin^2 phi)`

(defn elliptic-F "Incomplete elliptic integral of the first kind `F(phi|m)`." [phi m] (r/dispatch dual-elliptic-F s/elliptic-F phi m))
(defn elliptic-E "Incomplete elliptic integral of the second kind `E(phi|m)`." [phi m] (r/dispatch dual-elliptic-E s/elliptic-E phi m))
(defn elliptic-D "Incomplete Legendre integral `D(phi|m)`." [phi m] (r/dispatch dual-elliptic-D s/elliptic-D phi m))
(defn elliptic-PI "Incomplete elliptic integral of the third kind `PI(n; phi|m)`." [n phi m] (r/dispatch dual-elliptic-PI s/elliptic-PI n phi m))

(defn- sin-squared [phi] (let [s (sin phi)] (mul s s)))
(defn- delta [phi m] (el/sqrt (sub 1.0 (mul m (sin-squared phi)))))

(defn- dual-elliptic-F [phi m] (r/chain [phi phi, m m, _y (elliptic-F phi m)] :elliptic-F (div 1.0 (delta phi m)) nil))
(defn- dual-elliptic-E [phi m] (r/chain [phi phi, m m, _y (elliptic-E phi m)] :elliptic-E (delta phi m) nil))
(defn- dual-elliptic-D [phi m] (r/chain [phi phi, m m, _y (elliptic-D phi m)] :elliptic-D (div (sin-squared phi) (delta phi m)) nil))

(defn- dual-elliptic-PI
  [n phi m]
  (r/chain [n n, phi phi, m m, _y (elliptic-PI n phi m)] :elliptic-PI
    nil (div 1.0 (mul (sub 1.0 (mul n (sin-squared phi))) (delta phi m))) nil))
