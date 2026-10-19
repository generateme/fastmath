(ns fastmath.interpolation.linear
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.interpolation.common :as ic])
  (:import [fastmath.vector Vec3]
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- coeffs-array
  [^doubles xs ^doubles ys]
  (let [n (m/dec (alength xs))]
    (loop [i (long 0)
           arr []]
      (if (m/== i n)
        arr
        (let [i+ (m/inc i)
              x1 (Array/aget xs i)
              y1 (Array/aget ys i)]
          (recur (m/inc i) (conj arr (Vec3. x1 y1 (m// (m/- (Array/aget ys i+) y1)
                                                       (m/- (Array/aget xs i+) x1))))))))))

(defn linear
  "Creates a 1D piecewise linear interpolator.

  Connects consecutive data points with straight line segments.

  Parameters:

  - `xs` (sequence of numbers): nodes, finite and strictly increasing; at least 2 points.
  - `ys` (sequence of numbers): values, the same length as `xs`.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`. It equals `numpy.interp` inside the domain. Outside the domain `f` extrapolates linearly with the slope of the first or the last segment (like `scipy.interpolate.interp1d` with `fill_value='extrapolate'`).

  Throws `ex-info` (with `:type :invalid-input`, `:method :linear` and a `:property`) when the lengths differ, there are fewer than 2 points, `xs` contains `NaN` or an infinity, or `xs` is not strictly increasing (unsorted or duplicated).

  See also [[fastmath.interpolation.cubic/cubic]], [[fastmath.interpolation.monotone/monotone]], [[bilinear]]."
  [xs ys]
  (let [xs (m/seq->double-array xs)
        ys (m/seq->double-array ys)
        _ (ic/validate-xy! :linear xs ys {:min-points 2})
        coeffs (coeffs-array xs ys)]
    (fn ^double [^double x]
      (let [id (ic/binary-search xs x)
            ^Vec3 segment (coeffs id)]
        (m/+ (.y segment)
             (m/* (m/- x (.x segment)) (.z segment)))))))

(defn check-grid!
  "Checks the axes and the shape of the values of a rectangular grid interpolator and throws `ex-info` when they are invalid.

  Parameters:

  - `method` (keyword): name of the interpolation method, used in the message
  - `xs`, `ys` (double arrays): the axes, at least 2 finite and strictly increasing coordinates each
  - `vss` (`double[][]`): values with `(alength xs)` rows of `(alength ys)` values each

  Returns `nil`. Throws `ex-info` with `:property` set to `:too-few-points`, `:not-finite`, `:not-increasing` or `:shape-mismatch`. Internal helper shared by `bilinear` and `cubic-2d`.

  See also [[bilinear]], [[fastmath.interpolation.cubic/cubic-2d]]."
  {:no-doc true}
  [method ^doubles xs ^doubles ys vss]
  (ic/check-axis! method :xs xs {:min-points 2})
  (ic/check-axis! method :ys ys {:min-points 2})
  (when-not (instance? m/double-double-array-type vss)
    (ic/invalid-input! method :shape-mismatch
                       "vss has to be a non-empty sequence of rows of numbers"
                       {:expected-rows (alength xs)}))
  (when-not (m/== (alength ^"[[D" vss) (alength xs))
    (ic/invalid-input! method :shape-mismatch
                       (str "vss has to have one row for every x: expected " (alength xs)
                            " rows, got " (alength ^"[[D" vss))
                       {:expected-rows (alength xs) :rows (alength ^"[[D" vss)}))
  (dotimes [i (alength ^"[[D" vss)]
    (let [^doubles row (aget ^"[[D" vss i)]
      (when-not (m/== (alength row) (alength ys))
        (ic/invalid-input! method :shape-mismatch
                           (str "every row of vss has to have one value for every y: expected "
                                (alength ys) " values, got " (alength row) " in row " i)
                           {:expected-columns (alength ys) :columns (alength row) :row i})))))

(defn bilinear
  "Creates a 2D bilinear interpolator on a rectangular grid.

  Within each cell of the grid the value is a bilinear form of the four corner values.

  Parameters:

  - `xs` (sequence of numbers): nodes along the first axis, finite and strictly increasing; at least 2 points.
  - `ys` (sequence of numbers): nodes along the second axis, finite and strictly increasing; at least 2 points.
  - `vss` (sequence of sequences): values, `(count xs)` rows with `(count ys)` values each; `vss[i][j]` is the value at `(xs[i], ys[j])`.

  Returns a function `f` such that `(f x y)` or `(f [x y])` returns the interpolated `double` value. It equals `scipy.interpolate.RegularGridInterpolator` with `method='linear'` and `fill_value=None`. Outside the grid `f` extrapolates linearly using the formula of the nearest (edge) cell.

  Throws `ex-info` (with `:type :invalid-input`, `:method :bilinear` and a `:property`) when an axis has fewer than 2 points, contains `NaN` or an infinity, is not strictly increasing, or when the shape of `vss` does not match the axes (`:shape-mismatch`).

  See also [[linear]], [[fastmath.interpolation.cubic/cubic-2d]]."
  [xs ys vss]
  (let [vss (m/seq->double-double-array vss)
        xs (m/seq->double-array xs)
        ys (m/seq->double-array ys)]
    (check-grid! :bilinear xs ys vss)
    (fn bilinear-interpolation
      (^double [[^double x ^double y]] (bilinear-interpolation x y))
      (^double [^double x ^double y]
       (let [i (ic/binary-search xs x)
             j (ic/binary-search ys y)
             i+ (m/inc i)
             j+ (m/inc j)
             xi (Array/aget xs i)
             yj (Array/aget ys j)
             t (m// (m/- x xi) (m/- (Array/aget xs i+) xi))
             u (m// (m/- y yj) (m/- (Array/aget ys j+) yj))
             t- (m/- 1.0 t)
             u- (m/- 1.0 u)]
         (m/+ (m/* t- u- (Array/aget2d vss i  j))
              (m/* t  u- (Array/aget2d vss i+ j))
              (m/* t- u  (Array/aget2d vss i  j+))
              (m/* t  u  (Array/aget2d vss i+ j+))))))))
