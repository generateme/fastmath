(ns fastmath.interpolation.acm
  (:require [fastmath.core :as m]
            [fastmath.interpolation.common :as ic])
  (:import [org.apache.commons.math3.analysis UnivariateFunction]
           [org.apache.commons.math3.analysis.interpolation DividedDifferenceInterpolator NevilleInterpolator]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn neville
  "Creates the interpolating polynomial through all points using Neville's algorithm.

  The nodes `xs` have to be finite and distinct, at least 2 points, the order is not required. Throws `ex-info` with `:property` `:length-mismatch`, `:too-few-points`, `:not-finite` or `:duplicate-points` otherwise.

  Returns a function of one argument `x` returning a `double`.

  See also [[divided-difference]]."
  [xs ys]
  (let [xs (m/seq->double-array xs)
        ys (m/seq->double-array ys)]
    (ic/validate-xy! :neville xs ys {:min-points 2 :sorted? false})
    (try
      (let [^UnivariateFunction obj (.interpolate (NevilleInterpolator.) xs ys)]
        (fn ^double [^double x] (.value obj x)))
      (catch org.apache.commons.math3.exception.NonMonotonicSequenceException e
        (let [v (.getPrevious e)]
          (ic/invalid-input! :neville :duplicate-points
                             (str "xs has to be distinct, " v " is repeated")
                             {:axis :xs :value v}))))))

(defn divided-difference
  "Creates the interpolating polynomial through all points in Newton's form (divided differences).

  The nodes `xs` have to be finite and strictly increasing, at least 2 points. Throws `ex-info` with `:property` `:length-mismatch`, `:too-few-points`, `:not-finite` or `:not-increasing` otherwise.

  Returns a function of one argument `x` returning a `double`.

  See also [[neville]]."
  [xs ys]
  (let [xs (m/seq->double-array xs)
        ys (m/seq->double-array ys)]
    (ic/validate-xy! :divided-difference xs ys {:min-points 2})
    (let [^UnivariateFunction obj (.interpolate (DividedDifferenceInterpolator.) xs ys)]
      (fn ^double [^double x] (.value obj x)))))
