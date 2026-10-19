(ns fastmath.interpolation.ssj
  (:require [fastmath.core :as m]
            [fastmath.interpolation.common :as ic])
  (:import [fastmath.java Array]
           [umontreal.ssj.functionfit SmoothingCubicSpline]
           [umontreal.ssj.functions MathFunction]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn cubic-smoothing
  "Creates a cubic smoothing spline (Reinsch), a natural cubic spline which trades the closeness to the data against the roughness.

  The spline minimizes `rho * sum of weighted squared residuals + (1 - rho) * integral of the squared second derivative`. It equals `scipy.interpolate.make_smoothing_spline(x, y, lam=lam)` with `lam = (1 - rho) / rho`.

  Parameters:

  - `xs` (sequence of numbers): nodes, at least `2`, strictly increasing and finite
  - `ys` (sequence of numbers): values, same length as `xs`
  - `opts` (map, optional):
      - `:rho` (double, default: `1.0`): smoothing parameter in `[0, 1]`
      - `:weights` (sequence of numbers, default: `nil` meaning all `1`): positive finite weights of the residuals, same length as `xs`; the weights of scipy are the square roots of these weights

  The limits of `rho`: `1.0` gives the natural interpolating cubic spline (the data are not smoothed, equal to [[fastmath.interpolation.cubic/cubic]]), `0.0` is accepted and gives the linear regression line (the least squares line, the roughness is not allowed); values between give smoother curves for smaller `rho`. The documented interval in earlier versions was open at zero, the value `0.0` is valid.

  Returns a function of one argument `x` which returns the value of the spline. Outside of the data range the spline continues linearly with the tangent of the end point (natural boundary), it does not continue the end cubic. `NaN` gives `NaN`.

  Throws `ex-info` (`:type :invalid-input`) when there are less than two points, `xs` and `ys` or `weights` have different lengths, `xs` is not finite or not strictly increasing (also for duplicates), `:rho` is outside of `[0, 1]` (also `NaN`) or a weight is not positive and finite.

  See also [[fastmath.interpolation.cubic/cubic]], [[fastmath.interpolation/interpolation]]."
  ([xs ys] (cubic-smoothing xs ys nil))
  ([xs ys {:keys [weights rho]}]
   (let [xs (m/seq->double-array xs)
         ys (m/seq->double-array ys)
         rho (double (or rho 1.0))]
     (ic/validate-xy! :cubic-smoothing xs ys {:min-points 2 :sorted? true})
     (when-not (and (m/>= rho 0.0) (m/<= rho 1.0))
       (ic/invalid-input! :cubic-smoothing :rho
                          (str ":rho has to be in [0, 1], got " rho)
                          {:option :rho :value rho :expected "[0, 1]"}))
     (let [weights (when weights (m/seq->double-array weights))]
       (when weights
         (when-not (m/== (alength ^doubles weights) (alength xs))
           (ic/invalid-input! :cubic-smoothing :length-mismatch
                              (str "xs and weights have different lengths: " (alength xs) " and " (alength ^doubles weights))
                              {:xs-count (alength xs) :weights-count (alength ^doubles weights)}))
         (dotimes [i (alength ^doubles weights)]
           (let [w (Array/aget ^doubles weights i)]
             (when-not (and (m/valid-double? w) (m/pos? w))
               (ic/invalid-input! :cubic-smoothing :weights
                                  (str "weights have to be positive and finite, got " w " at index " i)
                                  {:option :weights :index i :value w})))))
       (let [^MathFunction obj (if weights
                                 (SmoothingCubicSpline. xs ys ^doubles weights rho)
                                 (SmoothingCubicSpline. xs ys rho))]
         (fn ^double [^double x] (.evaluate obj x)))))))
