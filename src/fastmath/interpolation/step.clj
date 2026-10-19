(ns fastmath.interpolation.step
  (:require [fastmath.core :as m]
            [fastmath.interpolation.common :as ic])
  (:import [org.apache.commons.math3.analysis.function StepFunction]
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn step-after
  "Creates a right-continuous piecewise constant (step) interpolator, the step occurs after each node.

  The value at `x` is the `y` of the nearest node on the left, that is the largest `xs[i]` that is less than or equal to `x`. Hence the function is continuous from the right: it jumps at each node and takes the value of that node there. It equals `scipy.interpolate.interp1d` with `kind='previous'` (and `steps-post` of matplotlib).

  Parameters:

  - `xs` (sequence of numbers): nodes, finite and strictly increasing; at least 1 point.
  - `ys` (sequence of numbers): values, the same length as `xs`.

  Returns a function `f` such that `(f x)` returns the `double` value at `x`. Left of the first node `f` returns the first value, right of the last node it returns the last value.

  Throws `ex-info` (with `:type :invalid-input`, `:method :step-after` and a `:property`) when the lengths differ, `xs` is empty, contains `NaN` or an infinity, or is not strictly increasing.

  See also [[step-before]], [[step]]."
  [xs ys]
  (let [xs (m/seq->double-array xs)
        ys (m/seq->double-array ys)
        _ (ic/validate-xy! :step-after xs ys {:min-points 1})
        ^StepFunction sf (StepFunction. xs ys)]
    (fn ^double [^double x] (.value sf x))))

(defn step-before
  "Creates a left-continuous piecewise constant (step) interpolator, the step occurs before each node.

  The value at `x` is the `y` of the nearest node on the right, that is the smallest `xs[i]` that is greater than or equal to `x`. Hence the function is continuous from the left: it jumps just after each node and takes the value of that node there. It equals `scipy.interpolate.interp1d` with `kind='next'` (and `steps-pre` of matplotlib).

  Parameters:

  - `xs` (sequence of numbers): nodes, finite and strictly increasing; at least 1 point.
  - `ys` (sequence of numbers): values, the same length as `xs`.

  Returns a function `f` such that `(f x)` returns the `double` value at `x`. Left of the first node `f` returns the first value, right of the last node it returns the last value.

  Throws `ex-info` (with `:type :invalid-input`, `:method :step-before` and a `:property`) when the lengths differ, `xs` is empty, contains `NaN` or an infinity, or is not strictly increasing.

  See also [[step-after]], [[step]]."
  [xs ys]
  (let [xs (m/seq->double-array xs)
        ys (m/seq->double-array ys)
        _ (ic/validate-xy! :step-before xs ys {:min-points 1})
        l (dec (alength xs))]
    (fn ^double [^double x]
      (let [b (java.util.Arrays/binarySearch xs x)
            ^long i (if (m/neg? b) (m/constrain (m/dec (m/- b)) 0 l) b)]
        (Array/aget ys i)))))

(defn step
  "Creates a piecewise constant interpolator with the step inside each interval.

  Within each interval between two neighbouring nodes `xs[i]` and `xs[i+1]` the function switches from `ys[i]` to `ys[i+1]` at the position `xs[i] + point * (xs[i+1] - xs[i])`. With the default `:point` 0.5 the value at `x` is the `y` of the nearest node and a tie (exactly the middle of the interval) goes to the right node, which equals `scipy.interpolate.interp1d` with `kind='nearest-up'`. The function takes the node value exactly at each node.

  Parameters:

  - `xs` (sequence of numbers): nodes, finite and strictly increasing; at least 2 points.
  - `ys` (sequence of numbers): values, the same length as `xs`.
  - `opts` (map, optional):
      - `:point` (double, default: 0.5): position of the step inside an interval as a fraction of its width. `0.0` gives the same function as [[step-before]] and `1.0` the same as [[step-after]].

  Returns a function `f` such that `(f x)` returns the `double` value at `x`. Left of the first node `f` returns the first value, right of the last node it returns the last value.

  Throws `ex-info` (with `:type :invalid-input`, `:method :step` and a `:property`) when the lengths differ, there are fewer than 2 points, `xs` contains `NaN` or an infinity, or is not strictly increasing.

  See also [[step-before]], [[step-after]]."
  ([xs ys] (step xs ys nil))
  ([xs ys {:keys [^double point]
           :or {point 0.5}}]
   (let [xs (m/seq->double-array xs)
         ys (m/seq->double-array ys)
         _ (ic/validate-xy! :step xs ys {:min-points 2})
         l (m/dec (alength xs))
         lst (m/dec (m/- l))]
     (fn ^double [^double x]
       (let [b (java.util.Arrays/binarySearch ^doubles xs x) 
             ^long i (cond
                       (m/== b -1) 0
                       (m/< b lst) l
                       (m/< b -1) (let [b-2 (m/- (m/- b) 2)
                                        b-1 (m/dec (m/- b))
                                        x1 (Array/aget xs b-2)
                                        x2 (Array/aget xs b-1)
                                        diff (m/+ x1 (m/* point (m/- x2 x1)))]
                                    (if (m/< x diff) b-2 b-1))
                       :else b)]
         (Array/aget ys i))))))

