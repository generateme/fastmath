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
         dual-abs dual-signum dual-floor dual-ceil dual-min dual-max
         dual-sec dual-csc dual-cot dual-sech dual-csch dual-coth dual-asec dual-acsc dual-acot dual-asech dual-acsch
         dual-acoth dual-sinpi dual-cospi dual-tanpi dual-cotpi dual-secpi dual-cscpi dual-sinc)

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
(defn sec "Secant." [x] (unary-op x dual-sec m/sec))
(defn csc "Cosecant." [x] (unary-op x dual-csc m/csc))
(defn cot "Cotangent." [x] (unary-op x dual-cot m/cot))
(defn sech "Hyperbolic secant." [x] (unary-op x dual-sech m/sech))
(defn csch "Hyperbolic cosecant." [x] (unary-op x dual-csch m/csch))
(defn coth "Hyperbolic cotangent." [x] (unary-op x dual-coth m/coth))
(defn asec "Arcus secant." [x] (unary-op x dual-asec m/asec))
(defn acsc "Arcus cosecant." [x] (unary-op x dual-acsc m/acsc))
(defn acot "Arcus cotangent." [x] (unary-op x dual-acot m/acot))
(defn asech "Area hyperbolic secant." [x] (unary-op x dual-asech m/asech))
(defn acsch "Area hyperbolic cosecant." [x] (unary-op x dual-acsch m/acsch))
(defn acoth "Area hyperbolic cotangent." [x] (unary-op x dual-acoth m/acoth))
(defn sinpi "`sin(pi x)`." [x] (unary-op x dual-sinpi m/sinpi))
(defn cospi "`cos(pi x)`." [x] (unary-op x dual-cospi m/cospi))
(defn tanpi "`tan(pi x)`." [x] (unary-op x dual-tanpi m/tanpi))
(defn cotpi "`cot(pi x)`." [x] (unary-op x dual-cotpi m/cotpi))
(defn secpi "`sec(pi x)`." [x] (unary-op x dual-secpi m/secpi))
(defn cscpi "`csc(pi x)`." [x] (unary-op x dual-cscpi m/cscpi))
(defn sinc "Normalized sinc, `sin(pi x) / (pi x)`." [x] (unary-op x dual-sinc m/sinc))
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

(defn- dual-sec [x] (r/chain [v x, y (sec v)] :sec (mul y (tan v))))
(defn- dual-csc [x] (r/chain [v x, y (csc v)] :csc (neg (mul y (cot v)))))
(defn- dual-cot [x] (r/chain [v x, y (cot v)] :cot (neg (add 1.0 (mul y y)))))
(defn- dual-sech [x] (r/chain [v x, y (sech v)] :sech (neg (mul (tanh v) y))))
(defn- dual-csch [x] (r/chain [v x, y (csch v)] :csch (neg (mul (coth v) y))))
(defn- dual-coth [x] (r/chain [v x, _y (coth v)] :coth (neg (sq (csch v)))))
(defn- dual-asec [x] (r/chain [v x, _y (asec v)] :asec (div 1.0 (mul (abs v) (sqrt (sub (mul v v) 1.0))))))
(defn- dual-acsc [x] (r/chain [v x, _y (acsc v)] :acsc (neg (div 1.0 (mul (abs v) (sqrt (sub (mul v v) 1.0)))))))
(defn- dual-acot [x] (r/chain [v x, _y (acot v)] :acot (neg (div 1.0 (add 1.0 (mul v v))))))
(defn- dual-asech [x] (r/chain [v x, _y (asech v)] :asech (neg (div 1.0 (mul v (sqrt (sub 1.0 (mul v v))))))))
(defn- dual-acsch [x] (r/chain [v x, _y (acsch v)] :acsch (neg (div 1.0 (mul (abs v) (sqrt (add 1.0 (mul v v))))))))
(defn- dual-acoth [x] (r/chain [v x, _y (acoth v)] :acoth (div 1.0 (sub 1.0 (mul v v)))))

;; The pi-scaled functions: the derivatives of the plain functions times `pi`.

(defn- dual-sinpi [x] (r/chain [v x, _y (sinpi v)] :sinpi (mul m/PI (cospi v))))
(defn- dual-cospi [x] (r/chain [v x, _y (cospi v)] :cospi (neg (mul m/PI (sinpi v)))))
(defn- dual-tanpi [x] (r/chain [v x, y (tanpi v)] :tanpi (mul m/PI (add 1.0 (mul y y)))))
(defn- dual-cotpi [x] (r/chain [v x, y (cotpi v)] :cotpi (neg (mul m/PI (add 1.0 (mul y y))))))
(defn- dual-secpi [x] (r/chain [v x, y (secpi v)] :secpi (mul m/PI (mul y (tanpi v)))))
(defn- dual-cscpi [x] (r/chain [v x, y (cscpi v)] :cscpi (neg (mul m/PI (mul y (cotpi v))))))

;; `sinc' = cosc(x) = (cos(pi x) - sinc(x)) / x`, which cancels near `0` (0/0 at `0`): below the threshold a Taylor
;; series in `z = pi x` is used, written with the operations on duals, so the higher derivatives at `0` are those
;; of the series (`sinc''(0) = -pi^2/3`). Measured against mpmath on `|x|` in `[1e-8, 1]`: at most `4.7e-16` (first
;; derivative) and `3.4e-15` (second); the direct form loses `1.7e-13` below `0.1` and all digits below `0.01`.

(def ^{:private true :const true :tag 'double} cosc-series-threshold 0.5)

(def ^:private cosc-coefficients
  "`cosc(x) = pi^2 x sum (-1)^(j+1) (2j+2) (pi x)^(2j) / (2j+3)!`, the factor of `pi^2 x` and the powers of `(pi x)^2`."
  (mapv (fn [^long j] (m// (m/* (if (even? j) -1.0 1.0) (double (m/long-add (m/long-mult 2 j) 2)))
                           (m/factorial (double (m/long-add (m/long-mult 2 j) 3)))))
        (range 12)))

(defn- cosc
  "`(cos(pi x) - sinc(x)) / x`, the derivative of [[sinc]]."
  [v]
  (if (m/< (m/abs (d/scalar v)) cosc-series-threshold)
    (let [z (mul m/PI v)]
      (mul (mul m/PI z) (ops/horner cosc-coefficients (mul z z))))
    (div (sub (cospi v) (sinc v)) v)))

(defn- dual-sinc [x] (r/chain [v x, _y (sinc v)] :sinc (cosc v)))

;; The norm and the angle: the partials are built from the norm `hypot`, without the squares (which under- or
;; overflow at moderate magnitudes); at an infinite argument the limits are used (`Inf / Inf` would be `NaN`).

(defn- infinite-count
  "The number of arguments with an infinite value, `0` when an argument is `NaN` (no limit is taken then)."
  [& xs]
  (let [ss (map d/scalar xs)]
    (if (some m/nan? ss) 0 (count (filter (fn [^double v] (Double/isInfinite v)) ss)))))

(defn- norm-partial
  "`c / h`, the partial of the norm `h` in its argument `c`. At an infinite norm with infinite arguments the limit:
  `sign(c)` when `c` is the only infinite argument, `0` for a finite `c`, `NaN` when several arguments are infinite
  (the limit depends on the direction)."
  [c h ^long infinite-arguments]
  (if (and (m/pos? infinite-arguments) (d/scalar-infinite? h))
    (cond
      (not (d/scalar-infinite? c)) 0.0
      (m/== 1 infinite-arguments) (m/signum (d/scalar c))
      :else ##NaN)
    (div c h)))

(defn- angle-partial
  "`c / (y^2 + x^2)` as `(c / h) / h` with `h = hypot(y, x)`; `0` at an infinite argument (the limit `1/|(y, x)|`)."
  [c h ^long infinite-arguments]
  (if (m/pos? infinite-arguments)
    0.0
    (div (div c h) h)))

(defn- dual-atan2
  [y x]
  (r/chain [a y, b x, _z (atan2 a b)] :atan2
    (angle-partial b (hypot a b) (infinite-count a b))
    (angle-partial (neg a) (hypot a b) (infinite-count a b))))

(defn- dual-hypot
  [x y]
  (r/chain [a x, b y, h (hypot a b)] :hypot
    (norm-partial a h (infinite-count a b))
    (norm-partial b h (infinite-count a b))))

(defn- dual-hypot3
  [x y z]
  (r/chain [a x, b y, c z, h (hypot a b c)] :hypot
    (norm-partial a h (infinite-count a b c))
    (norm-partial b h (infinite-count a b c))
    (norm-partial c h (infinite-count a b c))))

;; ## Kinks (as in ForwardDiff.jl)

(defn sign-bit?
  "Is the sign bit of the value of `x` set (also for `-0.0`)?"
  [x] (m/neg? (m/copy-sign 1.0 (d/scalar x))))

(defn- dual-abs
  "`-x` when the sign bit of the value is set, otherwise `x`: the derivative is `1` at `0.0` and `-1` at `-0.0`."
  [^Dual x]
  (if (sign-bit? (.value x)) (neg x) x))

(defn piecewise-constant
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
