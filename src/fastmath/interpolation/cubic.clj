;; Numerical recipes p.123
(ns fastmath.interpolation.cubic
  (:require [fastmath.core :as m]
            [fastmath.interpolation.common :as ic]
            [fastmath.interpolation.linear :as il])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(deftype Coeffs [^double x1 ^double x2
                 ^double y1 ^double y2
                 ^double y21 ^double y22
                 ^double h ^double h2])

;; not optimal, but calculated once
(defn- y''
  [^doubles xs ^doubles ys]
  (let [n (count xs)
        n- (m/dec n)
        u (double-array n-)
        y2 (double-array n)]
    (dotimes [iter (m/dec n-)]
      (let [i (inc iter)
            i- (m/dec i)
            i+ (m/inc i)
            x- (Array/aget xs i-)
            x (Array/aget xs i)
            x+ (Array/aget xs i+)
            sig (m// (m/- x x-) (m/- x+ x-))
            p (m/+ (m/* sig (Array/aget y2 i-)) 2.0)
            y2i (m// (m/dec sig) p)
            y (Array/aget ys i)
            ui (m/- (m// (m/- (Array/aget ys i+) y)
                         (m/- x+ x))
                    (m// (m/- y (Array/aget ys i-))
                         (m/- x x-)))
            ui (m// (m/- (m// (m/* 6.0 ui)
                              (m/- x+ x-))
                         (m/* sig (Array/aget u i-))) p)]
        (Array/aset u i ui)
        (Array/aset y2 i y2i)))
    (loop [k (m/dec n-)
           y2 y2]
      (if (m/neg? k)
        y2
        (recur (m/dec k)
               (Array/aset y2 k (m/+ (m/* (Array/aget y2 k)
                                          (Array/aget y2 (m/inc k ))) (Array/aget u k))))))))

(defn- coeffs-array
  [^doubles xs ^doubles ys]
  (let [n (m/dec (alength xs))
        ^doubles y2 (y'' xs ys)]
    (loop [i (long 0)
           arr []]
      (if (m/== i n)
        arr
        (let [i+ (m/inc i)
              x1 (Array/aget xs i)
              x2 (Array/aget xs i+)
              h (m/- x2 x1)]
          (recur (m/inc i)
                 (conj arr (Coeffs. x1 x2
                                    (Array/aget ys i) (Array/aget ys i+)
                                    (Array/aget y2 i) (Array/aget y2 i+)
                                    h (m// (* h h) 6.0)))))))))

(defn- natural-spline
  "Natural cubic spline function over double arrays, without any validation of the data."
  [^doubles xs ^doubles ys]
  (let [coeffs (coeffs-array xs ys)]
    (fn ^double [^double x]
      (let [id (ic/binary-search xs x)
            ^Coeffs segment (coeffs id)
            a (m// (m/- (.x2 segment) x) (.h segment))
            b (m// (m/- x (.x1 segment)) (.h segment))]
        (m/+ (m/* a (.y1 segment))
             (m/* b (.y2 segment))
             (m/* (m/+ (m/* (m/- (m/* a a a) a)
                            (.y21 segment))
                       (m/* (m/- (m/* b b b) b)
                            (.y22 segment))) (.h2 segment)))))))

(defn cubic
  "Creates a 1D natural cubic spline interpolator.

  The spline is piecewise cubic, passes through all data points and has continuous first and second derivatives. The natural boundary conditions set the second derivative to zero at both end points.

  Parameters:

  - `xs` (sequence of numbers): nodes, finite and strictly increasing; at least 2 points.
  - `ys` (sequence of numbers): values, the same length as `xs`.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`. It equals `scipy.interpolate.CubicSpline` with `bc_type='natural'`. Outside the domain the cubic of the first or the last interval is continued (extrapolation, like scipy; the second derivative is not zero there, so the extrapolation is cubic, not linear). With 2 points the spline is the straight line.

  Throws `ex-info` (with `:type :invalid-input`, `:method :cubic` and a `:property`) when the lengths differ, there are fewer than 2 points, `xs` contains `NaN` or an infinity, or `xs` is not strictly increasing (unsorted or duplicated).

  See also [[fastmath.interpolation.linear/linear]], [[fastmath.interpolation.monotone/monotone]], [[cubic-2d]]."
  [xs ys]
  (let [xs (m/seq->double-array xs)
        ys (m/seq->double-array ys)]
    (ic/validate-xy! :cubic xs ys {:min-points 2})
    (natural-spline xs ys)))

(defn cubic-2d
  "Creates a 2D natural cubic spline interpolator on a rectangular grid.

  The surface is a tensor product of natural cubic splines: for a query point the splines along `ys` through every row of `vss` are evaluated at `y`, and a natural cubic spline along `xs` through these values is evaluated at `x`. Note that the end conditions are natural (zero second derivative), not not-a-knot as in `scipy.interpolate.RectBivariateSpline`, so the two differ near the border of the grid. The result equals two passes of `scipy.interpolate.CubicSpline` with `bc_type='natural'`.

  Parameters:

  - `xs` (sequence of numbers): nodes along the first axis, finite and strictly increasing; at least 2 points.
  - `ys` (sequence of numbers): nodes along the second axis, finite and strictly increasing; at least 2 points.
  - `vss` (sequence of sequences): values, `(count xs)` rows with `(count ys)` values each; `vss[i][j]` is the value at `(xs[i], ys[j])`.

  Returns a function `f` such that `(f x y)` or `(f [x y])` returns the interpolated `double` value. Outside the grid the end cubics are continued (extrapolation). A spline along `xs` is built for every query, so evaluation costs O(number of rows) plus the spline set up.

  Throws `ex-info` (with `:type :invalid-input`, `:method :cubic-2d` and a `:property`) when an axis has fewer than 2 points, contains `NaN` or an infinity, is not strictly increasing, or when the shape of `vss` does not match the axes (`:shape-mismatch`).

  See also [[cubic]], [[fastmath.interpolation.linear/bilinear]]."
  [xs ys vss]
  (let [xs (m/seq->double-array xs)
        ys (m/seq->double-array ys)
        vss (m/seq->double-double-array vss)
        _ (il/check-grid! :cubic-2d xs ys vss)
        cbs (mapv (fn [^doubles row] (natural-spline ys row)) vss)
        n (alength xs)]
    (fn cubic-2d-interpolation
      (^double [[^double x ^double y]] (cubic-2d-interpolation x y))
      (^double [^double x ^double y]
       (let [nys (double-array n)]
         (dotimes [i n]
           (let [^clojure.lang.IFn$DD f (cbs i)]
             (Array/aset nys i (.invokePrim f y))))
         ((natural-spline xs nys) x))))))
