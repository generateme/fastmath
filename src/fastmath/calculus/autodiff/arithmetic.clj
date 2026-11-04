(ns ^:no-doc fastmath.calculus.autodiff.arithmetic
  "Internal: the polynomial and linear functions of `fastmath.core` (powers, fused multiply-add, linear interpolation and range mapping) for numbers and duals of [[fastmath.calculus.autodiff]], dispatched by `instance?`, and their derivative rules."
  (:require [fastmath.core :as m]
            [fastmath.calculus.autodiff.dual :as d]
            [fastmath.calculus.autodiff.ops :refer [add sub mul div neg]]
            [fastmath.calculus.autodiff.rules :as r]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(declare dual-cb dual-pow2 dual-pow10 dual-fma dual-lerp dual-norm dual-norm5)

;; ## Operations on numbers and duals

(defn cb "`x^3`." [x] (r/dispatch dual-cb m/cb x))
(defn pow2 "`x^2`." [x] (r/dispatch dual-pow2 m/pow2 x))
(defn pow10 "`x^10`." [x] (r/dispatch dual-pow10 m/pow10 x))
(defn fma "`x y + z`, fused." [x y z] (r/dispatch dual-fma m/fma x y z))
(defn lerp "Linear interpolation `start + t (stop - start)`." [start stop t] (r/dispatch dual-lerp m/lerp start stop t))

(defn- norm5-double [v start1 stop1 start2 stop2] (m/norm v start1 stop1 start2 stop2))

(defn norm
  "Linear map of `v` from `[start, stop]` to `[0, 1]`, or from `[start1, stop1]` to `[start2, stop2]`."
  ([v start stop] (r/dispatch dual-norm m/norm v start stop))
  ([v start1 stop1 start2 stop2] (r/dispatch dual-norm5 norm5-double v start1 stop1 start2 stop2)))

;; ## Derivative rules

(defn- dual-cb [x] (r/chain [v x, _y (cb v)] :cb (mul 3.0 (mul v v))))
(defn- dual-pow2 [x] (r/chain [v x, _y (pow2 v)] :pow2 (add v v)))
(defn- dual-pow10 [x] (r/chain [v x, _y (pow10 v)] :pow10 (let [c (mul (mul v v) v)] (mul 10.0 (mul (mul c c) c)))))

(defn- dual-fma [x y z] (r/chain [a x, b y, c z, _w (fma a b c)] :fma b a 1.0))

(defn- dual-lerp
  [start stop t]
  (r/chain [a start, b stop, c t, _y (lerp a b c)] :lerp (sub 1.0 c) c (sub b a)))

(defn- degenerate-range? [start stop] (m/== (d/scalar start) (d/scalar stop)))

(defn- norm-position
  "`(v - start) / (stop - start)`; for equal ends the constant of `fastmath.core/norm` (`0` or `1`)."
  [v start stop]
  (if (degenerate-range? start stop)
    (m/norm (d/scalar v) (d/scalar start) (d/scalar stop))
    (div (sub v start) (sub stop start))))

(defn- dual-norm
  "`1/w` in `v`, `(t - 1)/w` in `start`, `-t/w` in `stop`; `0` for equal ends (the value is a constant)."
  [v start stop]
  (r/chain [a v, b start, c stop, _y (norm a b c)] :norm
    (if (degenerate-range? b c) 0.0 (div 1.0 (sub c b)))
    (if (degenerate-range? b c) 0.0 (div (sub (norm-position a b c) 1.0) (sub c b)))
    (if (degenerate-range? b c) 0.0 (neg (div (norm-position a b c) (sub c b))))))

(defn- dual-norm5
  "The partials of [[dual-norm]] times `stop2 - start2`, `1 - t` in `start2`, `t` in `stop2`."
  [v start1 stop1 start2 stop2]
  (r/chain [a v, b start1, c stop1, e start2, f stop2, _y (norm a b c e f)] :norm
    (if (degenerate-range? b c) 0.0 (div (sub f e) (sub c b)))
    (if (degenerate-range? b c) 0.0 (div (mul (sub f e) (sub (norm-position a b c) 1.0)) (sub c b)))
    (if (degenerate-range? b c) 0.0 (neg (div (mul (sub f e) (norm-position a b c)) (sub c b))))
    (sub 1.0 (norm-position a b c))
    (norm-position a b c)))
