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

  If `hp1` is defined (degree < hp1 <= n) approximation function is created."
  ([xs ys] (b-spline xs ys nil))
  ([xs ys {:keys [^int degree hp1]
           :or {degree 3}}]
   (let [xs (m/seq->double-array xs)
         ys (m/seq->double-array ys)
         ^MathFunction obj (if hp1
                             (BSpline/createApproxBSpline xs ys degree (unchecked-int hp1))
                             (BSpline/createInterpBSpline xs ys degree))]
     (fn ^double [^double x] (.evaluate obj x)))))

(defn polynomial
  [xs ys]
  (let [^MathFunction obj (PolInterp. (m/seq->double-array xs) (m/seq->double-array ys))]
    (fn ^double [^double x] (.evaluate obj x))))
