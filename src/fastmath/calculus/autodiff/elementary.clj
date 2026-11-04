(ns ^:no-doc fastmath.calculus.autodiff.elementary
  "Internal: elementary functions of numbers and duals of [[fastmath.calculus.autodiff]] beyond the arithmetic, dispatched by `instance?`, and their derivative rules (the derivatives of DiffRules.jl, the kinks as in ForwardDiff.jl)."
  (:refer-clojure :exclude [abs min max])
  (:require [fastmath.core :as m]
            [fastmath.calculus.autodiff.dual :as d]
            [fastmath.calculus.autodiff.ops :as ops :refer [add sub mul div neg]]
            [fastmath.calculus.autodiff.rules :as r])
  (:import [fastmath.calculus.autodiff.dual Dual]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(declare dual-sqrt dual-cbrt dual-sq dual-tan dual-asin dual-acos dual-atan dual-atan2 dual-sinh dual-cosh dual-tanh
         dual-asinh dual-acosh dual-atanh dual-expm1 dual-log1p dual-log2 dual-log10 dual-hypot dual-hypot3
         dual-abs dual-signum dual-floor dual-ceil dual-min dual-max)

;; ## Operations on numbers and duals

(defmacro ^:private unary-op
  "A dual selects the dual rule, a number is coerced to a double and passed to the `fastmath.core` function."
  [x dual-rule core-fn]
  `(if (instance? Dual ~x) (~dual-rule ~x) (~core-fn (.doubleValue ~(with-meta x {:tag 'Number})))))

(defn sqrt "Square root." [x] (unary-op x dual-sqrt m/sqrt))
(defn cbrt "Cubic root." [x] (unary-op x dual-cbrt m/cbrt))
(defn sq "Square." [x] (unary-op x dual-sq m/sq))
(defn tan "Tangent." [x] (unary-op x dual-tan m/tan))
(defn asin "Arcus sine." [x] (unary-op x dual-asin m/asin))
(defn acos "Arcus cosine." [x] (unary-op x dual-acos m/acos))
(defn atan "Arcus tangent." [x] (unary-op x dual-atan m/atan))
(defn sinh "Hyperbolic sine." [x] (unary-op x dual-sinh m/sinh))
(defn cosh "Hyperbolic cosine." [x] (unary-op x dual-cosh m/cosh))
(defn tanh "Hyperbolic tangent." [x] (unary-op x dual-tanh m/tanh))
(defn asinh "Area hyperbolic sine." [x] (unary-op x dual-asinh m/asinh))
(defn acosh "Area hyperbolic cosine." [x] (unary-op x dual-acosh m/acosh))
(defn atanh "Area hyperbolic tangent." [x] (unary-op x dual-atanh m/atanh))
(defn expm1 "`exp(x) - 1`." [x] (unary-op x dual-expm1 m/expm1))
(defn log1p "`log(1 + x)`." [x] (unary-op x dual-log1p m/log1p))
(defn log2 "Logarithm base 2." [x] (unary-op x dual-log2 m/log2))
(defn log10 "Logarithm base 10." [x] (unary-op x dual-log10 m/log10))
(defn abs "Absolute value." [x] (unary-op x dual-abs m/abs))
(defn signum "Sign: `1`, `0` or `-1`." [x] (unary-op x dual-signum m/signum))
(defn floor "Rounding towards negative infinity." [x] (unary-op x dual-floor m/floor))
(defn ceil "Rounding towards positive infinity." [x] (unary-op x dual-ceil m/ceil))

(defn atan2
  "Angle of the point `(x, y)`, `Math/atan2(y, x)`: `y` first."
  [y x]
  (if (or (instance? Dual y) (instance? Dual x))
    (dual-atan2 y x)
    (m/atan2 (.doubleValue ^Number y) (.doubleValue ^Number x))))

(defn hypot
  "Euclidean norm of two or three arguments."
  ([x y]
   (if (or (instance? Dual x) (instance? Dual y))
     (dual-hypot x y)
     (m/hypot (.doubleValue ^Number x) (.doubleValue ^Number y))))
  ([x y z]
   (if (or (instance? Dual x) (instance? Dual y) (instance? Dual z))
     (dual-hypot3 x y z)
     (m/hypot (.doubleValue ^Number x) (.doubleValue ^Number y) (.doubleValue ^Number z)))))

(defn min
  "The smaller of two arguments."
  [x y]
  (if (or (instance? Dual x) (instance? Dual y))
    (dual-min x y)
    (m/min (.doubleValue ^Number x) (.doubleValue ^Number y))))

(defn max
  "The larger of two arguments."
  [x y]
  (if (or (instance? Dual x) (instance? Dual y))
    (dual-max x y)
    (m/max (.doubleValue ^Number x) (.doubleValue ^Number y))))

;; ## Derivative rules (the forms of DiffRules.jl, `x^2` as `x * x`)

(defn- dual-sqrt [x] (r/chain [v x, y (sqrt v)] :sqrt (div 1.0 (mul 2.0 y))))
(defn- dual-cbrt [x] (r/chain [v x, y (cbrt v)] :cbrt (div 1.0 (mul 3.0 (mul y y)))))
(defn- dual-sq [x] (r/chain [v x, _y (sq v)] :sq (add v v)))
(defn- dual-tan [x] (r/chain [v x, y (tan v)] :tan (add 1.0 (mul y y))))
(defn- dual-asin [x] (r/chain [v x, _y (asin v)] :asin (div 1.0 (sqrt (sub 1.0 (mul v v))))))
(defn- dual-acos [x] (r/chain [v x, _y (acos v)] :acos (neg (div 1.0 (sqrt (sub 1.0 (mul v v)))))))
(defn- dual-atan [x] (r/chain [v x, _y (atan v)] :atan (div 1.0 (add 1.0 (mul v v)))))
(defn- dual-sinh [x] (r/chain [v x, _y (sinh v)] :sinh (cosh v)))
(defn- dual-cosh [x] (r/chain [v x, _y (cosh v)] :cosh (sinh v)))
(defn- dual-tanh [x] (r/chain [v x, y (tanh v)] :tanh (sub 1.0 (mul y y))))
(defn- dual-asinh [x] (r/chain [v x, _y (asinh v)] :asinh (div 1.0 (sqrt (add (mul v v) 1.0)))))
(defn- dual-acosh [x] (r/chain [v x, _y (acosh v)] :acosh (div 1.0 (sqrt (sub (mul v v) 1.0)))))
(defn- dual-atanh [x] (r/chain [v x, _y (atanh v)] :atanh (div 1.0 (sub 1.0 (mul v v)))))
(defn- dual-expm1 [x] (r/chain [v x, _y (expm1 v)] :expm1 (ops/exp v)))
(defn- dual-log1p [x] (r/chain [v x, _y (log1p v)] :log1p (div 1.0 (add v 1.0))))
(defn- dual-log2 [x] (r/chain [v x, _y (log2 v)] :log2 (div (div 1.0 v) m/LN2)))
(defn- dual-log10 [x] (r/chain [v x, _y (log10 v)] :log10 (div (div 1.0 v) m/LN10)))

(defn- dual-atan2
  [y x]
  (r/chain [a y, b x, _z (atan2 a b)] :atan2
    (div b (add (mul a a) (mul b b)))
    (div (neg a) (add (mul a a) (mul b b)))))

(defn- dual-hypot
  [x y]
  (r/chain [a x, b y, h (hypot a b)] :hypot (div a h) (div b h)))

(defn- dual-hypot3
  [x y z]
  (r/chain [a x, b y, c z, h (hypot a b c)] :hypot (div a h) (div b h) (div c h)))

;; ## Kinks (as in ForwardDiff.jl)

(defn- sign-bit? [x] (m/neg? (m/copy-sign 1.0 (d/scalar x))))

(defn- dual-abs
  "`-x` when the sign bit of the value is set, otherwise `x`: the derivative is `1` at `0.0` and `-1` at `-0.0`."
  [^Dual x]
  (if (sign-bit? (.value x)) (neg x) x))

(defn- piecewise-constant
  "A dual with the value `(f value)` and zero partials."
  [f ^Dual x]
  (d/new-dual (.tag x) (f (.value x)) (d/p-zeros (.partials x))))

(defn- dual-signum [x] (piecewise-constant signum x))
(defn- dual-floor [x] (piecewise-constant floor x))
(defn- dual-ceil [x] (piecewise-constant ceil x))

(defn- dual-max
  "The value of `fastmath.core/max`, the partials of `x` when its value is larger, otherwise of `y` (at a tie and
  for `NaN`), as `max` of ForwardDiff.jl."
  [x y]
  (let [t (d/top-tag x y)
        vx (d/value-at x t)
        vy (d/value-at y t)
        px (d/partials-at x t)
        py (d/partials-at y t)
        p (if (m/> (d/scalar vx) (d/scalar vy)) px py)]
    (d/new-dual t (max vx vy) (or p (d/p-zeros (or px py))))))

(defn- dual-min
  "The value of `fastmath.core/min`, the partials of `y` when the value of `x` is larger, otherwise of `x` (at a tie
  and for `NaN`), as `min` of ForwardDiff.jl."
  [x y]
  (let [t (d/top-tag x y)
        vx (d/value-at x t)
        vy (d/value-at y t)
        px (d/partials-at x t)
        py (d/partials-at y t)
        p (if (m/> (d/scalar vx) (d/scalar vy)) py px)]
    (d/new-dual t (min vx vy) (or p (d/p-zeros (or px py))))))
