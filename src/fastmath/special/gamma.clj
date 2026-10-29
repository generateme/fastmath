(ns fastmath.special.gamma
  "Gamma function over the whole double range (internal; see [[fastmath.special/gamma]])."
  (:require [fastmath.core :as m]
            [fastmath.polynomials :as poly])
  (:import [fastmath.java Array]
           [fastmath.vector Vec2]
           [org.apache.commons.math3.special Gamma]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; Commons Math computes `pow(y, x+0.5) * exp(-y)` which overflows at x ~ 141 although
;; Gamma(x) is finite up to x ~ 171.62. Above this limit the power is split in two halves.
(def ^{:private true :const true :tag 'double} COMMONS-LIMIT 140.0)

(defn sin-pi
  "`sin(pi x)` with the argument reduced exactly: zero at every integer."
  ^double [^double x]
  (if (m/zero? x)
    x
    (let [r (m/rint x)
          s (m/sin (m/* m/PI (m/- x r)))]
      (if (m/even? (long r)) s (m/- s)))))

(defn cos-pi
  "`cos(pi x)` with the argument reduced exactly: `+0.0` at every half-integer (an even function)."
  ^double [^double x]
  (let [ax (m/abs x)]
    (if (m/== (m/- ax (m/floor ax)) 0.5)
      0.0
      (sin-pi (m/+ ax 0.5)))))

(def ^{:private true :const true :tag 'double} GAMMA-OVERFLOW 171.7)

(defn- gamma-large-pos
  "Gamma for `x > COMMONS-LIMIT`; `Inf` above the overflow point (the split power itself gives NaN above x = 1490)."
  ^double [^double x]
  (if (m/> x GAMMA-OVERFLOW)
    ##Inf
    (let [y (m/+ x Gamma/LANCZOS_G 0.5)
        s (m/* (m/pow y (m/* 0.5 (m/+ x 0.5))) (m/exp (m/* -0.5 y)))
        t (m/* s (m/sqrt (m/* (m// m/SQRT2PI x) (Gamma/lanczos x))))]
      (m/* t t))))

(defn- gamma-neg
  "Reflection formula for non-integer `x < -20` (Commons Math reflects with an inexact `sin(pi x)` near the poles)."
  ^double [^double x]
  (let [z (m/- 1.0 x)
        sp (sin-pi x)
        g (if (m/> z COMMONS-LIMIT) (gamma-large-pos z) (Gamma/gamma z))]
    (if (m/inf? g)
      ;; Gamma(1-x) overflows but the result is a subnormal or zero
      (m/* (m/signum sp)
           (m/exp (m/- (m/log (m// m/PI (m/abs sp))) (Gamma/logGamma z))))
      (m// m/PI (m/* sp g)))))

(defn gamma
  "Gamma function for any double: `NaN` at the poles (zero and negative integers) and at `-Inf`, `Inf` at `Inf`, overflow to `Inf` above 171.62."
  ^double [^double x]
  (cond
    (m/nan? x) ##NaN
    (m/pos-inf? x) ##Inf
    (m/neg-inf? x) ##NaN
    (m/> x COMMONS-LIMIT) (gamma-large-pos x)
    (m/< x -20.0) (if (m/== x (m/rint x)) ##NaN (gamma-neg x))
    :else (Gamma/gamma x)))

(defn unsafe-gamma
  "`Gamma(x)` with `+Inf` (instead of `NaN`) at the poles, so that `1/Gamma` is `0` there (for the hypergeometric transformations)."
  ^double [^double x]
  (if (and (m/not-pos? x) (m/== x (m/rint x)))
    ##Inf
    (gamma x)))

(defn log-gamma
  "`log|Gamma(x)|` for any double: `Inf` at the poles and at `Inf`, `NaN` at `-Inf`."
  ^double [^double x]
  (cond
    (m/nan? x) ##NaN
    (m/pos-inf? x) ##Inf
    (m/neg-inf? x) ##NaN
    (m/pos? x) (Gamma/logGamma x)
    (m/== x (m/rint x)) ##Inf
    :else (m/- (m/log (m// m/PI (m/abs (sin-pi x)))) (Gamma/logGamma (m/- 1.0 x)))))

(defn gamma-sign
  "Sign of `Gamma(x)` (`+-1.0`) for a non-pole `x`."
  ^double [^double x]
  (if (m/pos? x) 1.0 (m/signum (sin-pi x))))

(defn gamma-ratio-fallback
  "Log-space value of the gamma ratio, used when the products overflow or underflow.

  `ns` and `ds` are the numerator and denominator arguments."
  ^double [^doubles ns ^doubles ds]
  (loop [i 0 sgn 1.0 lg 0.0]
    (if (m/< i (alength ns))
      (let [x (Array/aget ns i)]
        (recur (m/inc i) (m/* sgn (gamma-sign x)) (m/+ lg (log-gamma x))))
      (loop [j 0 sgn sgn lg lg]
        (if (m/< j (alength ds))
          (let [x (Array/aget ds j)]
            (recur (m/inc j) (m/* sgn (gamma-sign x)) (m/- lg (log-gamma x))))
          (m/* sgn (m/exp lg)))))))

(defn gamma-pole?
  "True for a pole of Gamma: zero and the negative integers."
  [^double x]
  (and (m/not-pos? x) (m/== x (m/rint x))))

(defmacro gamma-ratio
  "`Gamma(n1)...Gamma(nk) / (Gamma(d1)...Gamma(dl))`.

  Evaluated directly (a pole gives `Inf`, as [[unsafe-gamma]]); when the products overflow or underflow and no argument is a pole it is evaluated through `log|Gamma|`."
  [nums dens]
  (let [ns (repeatedly (count nums) #(gensym "n"))
        ds (repeatedly (count dens) #(gensym "d"))]
    `(let [~@(interleave ns nums)
           ~@(interleave ds dens)
           num# (m/* ~@(map (fn [s] `(unsafe-gamma ~s)) ns))
           den# (m/* ~@(map (fn [s] `(unsafe-gamma ~s)) ds))]
       (if (and (not (m/nan? num#)) (not (m/nan? den#))
                (or (m/inf? num#) (m/inf? den#) (m/zero? num#) (m/zero? den#))
                (not (or ~@(map (fn [x] `(gamma-pole? ~x)) (concat ns ds)))))
         (gamma-ratio-fallback (double-array [~@ns]) (double-array [~@ds]))
         (m// num# den#)))))

(defn digamma
  "Digamma function `psi(x)`: `NaN` at the negative integers, `-Inf` at `+0.0` and `Inf` at `-0.0` (one-sided limits)."
  ^double [^double x]
  (if (and (m/not-pos? x) (m/== x (m/rint x)))
    (if (m/zero? x)
      (m/- (m// 1.0 x))
      ##NaN)
    (let [^Vec2 psix (if (m/not-pos? x)
                       (Vec2. (m/* m/-PI (m// (cos-pi x) (sin-pi x))) (m/- 1.0 x))
                       (Vec2. 0.0 x))
          x (.y psix)
          ^Vec2 psix (if (m/< x 8.0)
                       (let [n (unchecked-long (m/- 8.0 (m/floor x)))]
                         (loop [v (long 1)
                                psi (.x psix)]
                           (if (m/== v n)
                             (Vec2. (m/- psi (m// x)) (m/+ x n))
                             (recur (m/inc v) (m/- psi (m// (m/+ x v)))))))
                       psix)
          t (m// (.y psix))
          psi (m/+ (.x psix) (m/- (m/log (.y psix)) (m/* 0.5 t)))
          t (m/* t t)]
      (m/- psi (m/* t (poly/mevalpoly t 0.08333333333333333,-0.008333333333333333,0.003968253968253968,-0.004166666666666667,0.007575757575757576,-0.021092796092796094,0.08333333333333333,-0.4432598039215686))))))
