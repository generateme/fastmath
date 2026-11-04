(ns ^:no-doc fastmath.calculus.autodiff.ops
  "Internal: arithmetic and elementary functions on numbers and duals of [[fastmath.calculus.autodiff]], dispatched by `instance?`, and the derivative rules for duals."
  (:require [fastmath.core :as m]
            [fastmath.calculus.autodiff.dual :as d])
  (:import [fastmath.calculus.autodiff.dual Dual]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; `dual-cmp` is declared with its primitive signature, so `cmp` calls it without boxing the result.
(declare dual-add dual-sub dual-mul dual-div dual-pow
         dual-neg dual-sin dual-cos dual-exp dual-log)
(declare ^{:arglists '(^long [x y])} dual-cmp)

;; ## Operations on numbers and duals
;;
;; A dual in any argument selects the dual rule, otherwise the arguments are numbers, coerced to doubles. A protocol
;; is not used: its method cache keeps one receiver class, and the receivers alternate between `Dual` and `Double`
;; in nested computations, which makes the cache miss on almost every call.

(defn- cmp-doubles
  ^long [^double a ^double b]
  (cond (m/< a b) -1
        (m/> a b) 1
        (m/== a b) 0
        :else 2))

(defn add
  "Sum of `x` and `y`."
  [x y]
  (if (or (instance? Dual x) (instance? Dual y))
    (dual-add x y)
    (m/+ (.doubleValue ^Number x) (.doubleValue ^Number y))))

(defn sub
  "Difference of `x` and `y`."
  [x y]
  (if (or (instance? Dual x) (instance? Dual y))
    (dual-sub x y)
    (m/- (.doubleValue ^Number x) (.doubleValue ^Number y))))

(defn mul
  "Product of `x` and `y`."
  [x y]
  (if (or (instance? Dual x) (instance? Dual y))
    (dual-mul x y)
    (m/* (.doubleValue ^Number x) (.doubleValue ^Number y))))

(defn div
  "Quotient of `x` and `y`."
  [x y]
  (if (or (instance? Dual x) (instance? Dual y))
    (dual-div x y)
    (m// (.doubleValue ^Number x) (.doubleValue ^Number y))))

(defn pow
  "`x` raised to the power `y`."
  [x y]
  (if (or (instance? Dual x) (instance? Dual y))
    (dual-pow x y)
    (m/pow (.doubleValue ^Number x) (.doubleValue ^Number y))))

(defn cmp
  "Lexicographic comparison of (value, partials): `-1`, `0` or `1`, and `2` when unordered (NaN)."
  ^long [x y]
  (if (or (instance? Dual x) (instance? Dual y))
    (dual-cmp x y)
    (cmp-doubles (.doubleValue ^Number x) (.doubleValue ^Number y))))

(defn neg
  "Negation of `x`."
  [x]
  (if (instance? Dual x) (dual-neg x) (m/- (.doubleValue ^Number x))))

(defn sin
  "Sine of `x`."
  [x]
  (if (instance? Dual x) (dual-sin x) (m/sin (.doubleValue ^Number x))))

(defn cos
  "Cosine of `x`."
  [x]
  (if (instance? Dual x) (dual-cos x) (m/cos (.doubleValue ^Number x))))

(defn exp
  "Exponential of `x`."
  [x]
  (if (instance? Dual x) (dual-exp x) (m/exp (.doubleValue ^Number x))))

(defn log
  "Natural logarithm of `x`."
  [x]
  (if (instance? Dual x) (dual-log x) (m/log (.doubleValue ^Number x))))

(defn horner
  "The polynomial with the coefficients `cs` (a vector of doubles, the constant term first) at `u`, a number or a dual."
  [cs u]
  (reduce (fn [acc c] (add c (mul u acc))) 0.0 (rseq cs)))

;; ## Partials arithmetic (`nil` is a zero vector)

(defn p-add
  "Sum of partials, `nil` is a zero vector."
  [px py]
  (cond (nil? px) py
        (nil? py) px
        :else (do (d/check-same-length! px py) (d/p-map2 add px py))))

(defn- p-sub [px py]
  (cond (nil? py) px
        (nil? px) (d/p-map neg py)
        :else (do (d/check-same-length! px py) (d/p-map2 sub px py))))

(defn p-scale
  "Each partial times `s` (element on the left), without allocating a closure; `nil` for `nil`."
  [p s]
  (when p
    (let [^objects p p
          n (alength p)
          a (object-array n)]
      (loop [i 0]
        (when (m/< i n)
          (aset a i (mul (aget p i) s))
          (recur (m/long-inc i))))
      a)))

(defn- p-div
  "Each partial divided by `s`, without allocating a closure."
  [p s]
  (when p
    (let [^objects p p
          n (alength p)
          a (object-array n)]
      (loop [i 0]
        (when (m/< i n)
          (aset a i (div (aget p i) s))
          (recur (m/long-inc i))))
      a)))

(defn- p-cmp
  "Lexicographic comparison of partials, `nil` is a zero vector."
  ^long [px py]
  (when (and px py) (d/check-same-length! px py))
  (let [^objects p (or px py)
        n (alength p)]
    (loop [i 0]
      (if (m/< i n)
        (let [c (cmp (if px (aget ^objects px i) 0.0) (if py (aget ^objects py i) 0.0))]
          (if (m/zero? c) (recur (m/long-inc i)) c))
        0))))

;; ## Derivative rules (at least one argument is a dual)

(defn- dual-add [x y]
  (let [t (d/top-tag x y)]
    (d/new-dual t (add (d/value-at x t) (d/value-at y t)) (p-add (d/partials-at x t) (d/partials-at y t)))))

(defn- dual-sub [x y]
  (let [t (d/top-tag x y)]
    (d/new-dual t (sub (d/value-at x t) (d/value-at y t)) (p-sub (d/partials-at x t) (d/partials-at y t)))))

(defn- dual-mul [x y]
  (let [t (d/top-tag x y)
        vx (d/value-at x t)
        vy (d/value-at y t)]
    (d/new-dual t (mul vx vy) (p-add (p-scale (d/partials-at x t) vy) (p-scale (d/partials-at y t) vx)))))

(defn- dual-div [x y]
  (let [t (d/top-tag x y)
        vx (d/value-at x t)
        vy (d/value-at y t)
        q (div vx vy)]
    ;; d(x/y) = dx/y - dy * (x/y)/y
    (d/new-dual t q (p-sub (p-div (d/partials-at x t) vy) (p-scale (d/partials-at y t) (div q vy))))))

(defn- positive-element? [x] (m/== 1 (cmp x 0.0)))

(defn- pow-dual-const
  "`x` is tagged `t`, `y` is a constant in this layer."
  [x y ^long t]
  (let [vx (d/value-at x t)
        px (d/partials-at x t)]
    (d/new-dual t (pow vx y)
                (if (or (d/zero-element? y) (d/p-zero? px))
                  (d/p-zeros px)
                  (p-scale px (mul y (pow vx (sub y 1.0))))))))

(defn- pow-const-dual
  "`x` is a constant in this layer, `y` is tagged `t`."
  [x y ^long t]
  (let [vy (d/value-at y t)
        py (d/partials-at y t)
        e (pow x vy)]
    (d/new-dual t e
                (if (or (d/p-zero? py) (and (d/zero-element? x) (positive-element? vy)))
                  (d/p-zeros py)
                  (p-scale py (mul e (log x)))))))

(defn- pow-dual-dual
  "Both `x` and `y` are tagged `t`."
  [x y ^long t]
  (let [px (d/partials-at x t)
        py (d/partials-at y t)
        vx (d/value-at x t)
        vy (d/value-at y t)]
    (cond
      (d/p-zero? py) (pow-dual-const x vy t)
      (d/p-zero? px) (pow-const-dual vx y t)
      :else (let [e (pow vx vy)
                  dx (mul vy (pow vx (sub vy 1.0)))
                  dy (if (and (d/zero-element? vx) (positive-element? vy)) 0.0 (mul e (log vx)))]
              (d/new-dual t e (p-add (p-scale px dx) (p-scale py dy)))))))

(defn- dual-pow [x y]
  (let [t (d/top-tag x y)]
    (cond
      (nil? (d/partials-at y t)) (pow-dual-const x y t)
      (nil? (d/partials-at x t)) (pow-const-dual x y t)
      :else (pow-dual-dual x y t))))

(defn- dual-cmp
  ^long [x y]
  (let [t (d/top-tag x y)
        c (cmp (d/value-at x t) (d/value-at y t))]
    (if (m/zero? c) (p-cmp (d/partials-at x t) (d/partials-at y t)) c)))

(defn- dual-neg [^Dual x]
  (d/new-dual (.tag x) (neg (.value x)) (d/p-map neg (.partials x))))

(defn- dual-sin [^Dual x]
  (let [v (.value x)]
    (d/new-dual (.tag x) (sin v) (p-scale (.partials x) (cos v)))))

(defn- dual-cos [^Dual x]
  (let [v (.value x)]
    (d/new-dual (.tag x) (cos v) (p-scale (.partials x) (neg (sin v))))))

(defn- dual-exp [^Dual x]
  (let [e (exp (.value x))]
    (d/new-dual (.tag x) e (p-scale (.partials x) e))))

(defn- dual-log
  "Outside the domain (a negative argument, including `-Inf`) the value is `NaN` and so are the partials: `1/x` is
  scaled by `NaN`. Scaling the quotient, not the partials, keeps the perturbations of the outer layers in it, so the
  higher derivatives are `NaN` too."
  [^Dual x]
  (let [v (.value x)
        p (p-div (.partials x) v)]
    (d/new-dual (.tag x) (log v) (if (m/neg? (d/scalar v)) (p-scale p ##NaN) p))))
