(ns ^:no-doc fastmath.calculus.autodiff.piecewise
  "Internal: the piecewise functions of `fastmath.core` (rounding, fractional parts, `mod`, `rem`, sign functions, the symmetric power, `safe-sqrt`, `smoothstep`) for numbers and duals of [[fastmath.calculus.autodiff]], dispatched by `instance?`, and their derivative rules. The rounding functions have the derivative `0`; at the discontinuities of the fractional parts and of `mod` and `rem` the derivatives are `NaN`, as in DiffRules.jl."
  (:refer-clojure :exclude [mod rem])
  (:require [fastmath.core :as m]
            [fastmath.calculus.autodiff.dual :as d]
            [fastmath.calculus.autodiff.ops :as ops :refer [sub mul div neg log pow]]
            [fastmath.calculus.autodiff.elementary :as el]
            [fastmath.calculus.autodiff.rules :as r])
  (:import [fastmath.calculus.autodiff.dual Dual]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(declare dual-trunc dual-rint dual-round dual-sgn dual-frac dual-sfrac dual-mod dual-rem dual-copy-sign
         dual-spow dual-safe-sqrt dual-smoothstep)

;; ## Operations on numbers and duals

(defn- round-double ^double [^double x] (double (m/round x)))

(defn trunc "Rounding towards zero." [x] (r/dispatch dual-trunc m/trunc x))
(defn rint "Rounding to the nearest integer, half to even." [x] (r/dispatch dual-rint m/rint x))
(defn round "Rounding to the nearest integer, half up, as a double." [x] (r/dispatch dual-round round-double x))
(defn sgn "`-1` for a negative value, otherwise `1`." [x] (r/dispatch dual-sgn m/sgn x))
(defn frac "Unsigned fractional part." [x] (r/dispatch dual-frac m/frac x))
(defn sfrac "Signed fractional part." [x] (r/dispatch dual-sfrac m/sfrac x))
(defn mod "Floored modulus." [x y] (r/dispatch dual-mod m/mod x y))
(defn rem "Truncated remainder." [x y] (r/dispatch dual-rem m/rem x y))
(defn copy-sign "The magnitude of `magnitude` with the sign of `sign`." [magnitude sign] (r/dispatch dual-copy-sign m/copy-sign magnitude sign))
(defn spow "Symmetric power `sgn(x) |x|^e`." [x e] (r/dispatch dual-spow m/spow x e))
(defn safe-sqrt "Square root, `0` for a negative value." [x] (r/dispatch dual-safe-sqrt m/safe-sqrt x))
(defn smoothstep "Clamped Hermite step from `edge0` to `edge1`." [edge0 edge1 x] (r/dispatch dual-smoothstep m/smoothstep edge0 edge1 x))


;; ## Derivative rules: rounding and kinks

(defn- dual-trunc [x] (el/piecewise-constant trunc x))
(defn- dual-rint [x] (el/piecewise-constant rint x))
(defn- dual-round [x] (el/piecewise-constant round x))
(defn- dual-sgn [x] (el/piecewise-constant sgn x))

(defn- integer-valued?
  "Is the double `z` a finite integer (as `isinteger` of Julia)?"
  [^double z]
  (and (not (Double/isInfinite z)) (m/== z (m/floor z))))

(defn- nonzero-integer? [x] (let [s (d/scalar x)] (and (not (m/zero? s)) (integer-valued? s))))

(defn- dual-frac
  "The slope is the sign of the argument by its sign bit (`-1` at `-0.0`, as [[fastmath.calculus.autodiff/abs]]),
  `NaN` at the nonzero integers (the jumps)."
  [x]
  (r/chain [v x, _y (frac v)] :frac (cond (nonzero-integer? v) ##NaN (el/sign-bit? v) -1.0 :else 1.0)))

(defn- dual-sfrac [x] (r/chain [v x, _y (sfrac v)] :sfrac (if (nonzero-integer? v) ##NaN 1.0)))

(defn- quotient-slope
  "`(f (x/y))` of the scalar quotient, `NaN` at an integer quotient (the jumps)."
  ^double [f x y]
  (let [z (m// (d/scalar x) (d/scalar y))]
    (if (integer-valued? z) ##NaN (double (f z)))))

(defn- one [_] 1.0)
(defn- neg-floor [^double z] (m/- (m/floor z)))
(defn- neg-trunc [^double z] (m/- (m/trunc z)))

(defn- dual-mod
  "`(1, -floor(x/y))`, `NaN` in both at an integer quotient (the jumps), as DiffRules.jl."
  [x y]
  (r/chain [a x, b y, _z (mod a b)] :mod (quotient-slope one a b) (quotient-slope neg-floor a b)))

(defn- dual-rem
  "`(1, -trunc(x/y))`, `NaN` in both at an integer quotient (the jumps), as DiffRules.jl."
  [x y]
  (r/chain [a x, b y, _z (rem a b)] :rem (quotient-slope one a b) (quotient-slope neg-trunc a b)))

(defn- dual-copy-sign
  "`1` in the magnitude when the sign bits of the arguments agree, otherwise `-1`; `0` in the sign."
  [magnitude sign]
  (r/chain [a magnitude, b sign, _y (copy-sign a b)] :copy-sign
    (if (= (el/sign-bit? a) (el/sign-bit? b)) 1.0 -1.0)
    0.0))

;; ## Derivative rules: the symmetric power, safe-sqrt, smoothstep


(defn- dual-spow
  "`e |x|^(e-1)` in `x`; `spow(x, e) log|x|` in `e`, `0` at `x = 0` for `e > 0` (the limit, as [[pow]] with a
  constant base)."
  [x e]
  (r/chain [a x, b e, y (spow a b)] :spow
    (mul b (pow (el/abs a) (sub b 1.0)))
    (if (and (m/zero? (d/scalar a)) (m/pos? (d/scalar b))) 0.0 (mul y (log (el/abs a))))))

(defn- dual-safe-sqrt
  "Zero partials for a negative argument (the value is the constant `0`), otherwise the square root."
  [^Dual x]
  (if (m/neg? (d/scalar x))
    (el/piecewise-constant safe-sqrt x)
    (r/chain [v x, y (safe-sqrt v)] :safe-sqrt (div 1.0 (mul 2.0 y)))))

(defn- step-slope-and-position
  "`[s t]` with `t = (x - edge0) / w`, `s = 6 t (1 - t) / w`, `w = edge1 - edge0`, the slope of the Hermite step and
  the normalized position; `nil` where the position is clamped (outside the edges, at the edges, for equal edges)."
  [edge0 edge1 x]
  (when (m/< 0.0 (m/norm (d/scalar x) (d/scalar edge0) (d/scalar edge1)) 1.0)
    (let [w (sub edge1 edge0)
          t (div (sub x edge0) w)]
      [(div (mul 6.0 (mul t (sub 1.0 t))) w) t])))

(defn- dual-smoothstep
  "Inside the edges: `dx = s`, `dedge0 = -s (1 - t)`, `dedge1 = -s t`; `0` where clamped."
  [edge0 edge1 x]
  (r/chain [a edge0, b edge1, c x, _y (smoothstep a b c)] :smoothstep
    (if-let [[s t] (step-slope-and-position a b c)] (neg (mul s (sub 1.0 t))) 0.0)
    (if-let [[s t] (step-slope-and-position a b c)] (neg (mul s t)) 0.0)
    (if-let [[s] (step-slope-and-position a b c)] s 0.0)))
