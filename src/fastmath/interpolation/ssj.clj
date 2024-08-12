(ns fastmath.interpolation.ssj
  (:require [fastmath.core :as m])
  (:import [umontreal.ssj.functionfit BSpline PolInterp SmoothingCubicSpline]
           [umontreal.ssj.functions MathFunction]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn cubic-smoothing
  ([xs ys] (cubic-smoothing xs ys nil))
  ([xs ys {:keys [weights ^double rho]
           :or {rho 1.0}}]
   (let [xs (m/seq->double-array xs)
         ys (m/seq->double-array ys)
         ^MathFunction obj (if weights
                             (SmoothingCubicSpline. xs ys (m/seq->double-array weights) rho)
                             (SmoothingCubicSpline. xs ys rho))]
     (fn ^double [^double x] (.evaluate obj x)))))

(defn b-spline
  "B-spline interpolation or approximation.

  Parameters:
  * `:knots` or `:degree` - sequence of knots or degree
  * `:clamped?` - for clamped knots (calculated), default: `false`
  * `:hp1` - for clamped knots, if defined (degree < hp1 <= n) approximation function is created.

  For clamped knots, default degree is set to `3`, `N-1` otherwise."
  ([xs ys] (b-spline xs ys nil))
  ([xs ys {:keys [^int degree hp1 clamped? knots]
           :or {clamped? false}}]
   (let [xs (m/seq->double-array xs)
         ys (m/seq->double-array ys)
         ^MathFunction obj (cond
                             knots (BSpline. xs ys (m/seq->double-array knots))
                             clamped? (let [degree (unchecked-int (or degree 3))]
                                        (if hp1
                                          (BSpline/createApproxBSpline xs ys degree (unchecked-int hp1))
                                          (BSpline/createInterpBSpline xs ys degree)))
                             :else (BSpline. xs ys (unchecked-int (or degree (m/dec (count xs))))))]
     (fn ^double [^double x] (.evaluate obj x)))))

(defn polynomial
  [xs ys]
  (let [^MathFunction obj (PolInterp. (m/seq->double-array xs) (m/seq->double-array ys))]
    (fn ^double [^double x] (.evaluate obj x))))

