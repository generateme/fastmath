(ns fastmath.interpolation
  "Interpolation and extrapolation functions for 1D, 2D grid, and multivariate data.

  Every interpolator is constructed from a set of known data points and returns a Clojure
  function that evaluates the fitted surface at any query point. All constructors are available
  both as pure named functions and via the unified [[interpolation]] multimethod, which selects
  the implementation from a keyword.

  The data is validated at construction (equal lengths, minimal number of points, finite and, where
  the method needs it, strictly increasing `xs`); invalid data throws `ex-info` with `:type`
  `:invalid-input`, the `:method` and the violated `:property` in `ex-data`.

  1D interpolators (inputs `xs`, `ys`):

  - [[linear]] - piecewise linear; extrapolates beyond domain
  - [[cubic]] - natural cubic spline; C2-continuous
  - [[monotone]] - monotonicity-preserving cubic (Fritsch-Carlson)
  - [[pchip]] - shape-preserving cubic Hermite (Fritsch-Butland PCHIP)
  - [[akima]] - Akima and modified Akima (makima) spline; robust to outliers; extrapolates
  - [[sprague]] - fifth-order Sprague; accurate for smooth uniform data
  - [[neville]] - polynomial via Neville's algorithm
  - [[divided-difference]] - polynomial via Newton's divided differences
  - [[polynomial]] - polynomial in the barycentric form
  - [[barycentric]] - numerically stable barycentric rational interpolation (Floater-Hormann)
  - [[b-spline]] - interpolating B-spline of any degree
  - [[step-before]] - left-continuous piecewise constant
  - [[step-after]] - right-continuous piecewise constant
  - [[step]] - midpoint-blended piecewise constant

  1D smoothing interpolators (inputs `xs`, `ys`; do not pass through every point):

  - [[loess]] - locally weighted scatterplot smoothing (LOESS)
  - [[cubic-smoothing]] - roughness-penalised smoothing cubic spline

  1D isotonic regression interpolators (inputs `xs`, `ys`; enforce monotonicity):

  - [[pava]] - Pool-Adjacent-Violators Algorithm; fits `ys` only
  - [[cir]] - Centered Isotonic Regression; shrinks both `xs` and `ys`

  2D grid interpolators (inputs `xs`, `ys`, `vss`; return `(f x y)` or `(f [x y])`):

  - [[bilinear]] - bilinear interpolation on a rectangular grid
  - [[bicubic]] - bicubic Hermite interpolation on a rectangular grid
  - [[cubic-2d]] - 2D cubic spline on a rectangular grid

  Multivariate interpolators (inputs `xss`, `ys`; work in any number of dimensions):

  - [[microsphere-projection]] - microsphere projection
  - [[shepard]] - inverse distance weighting
  - [[rbf]] - radial basis function with optional regularisation and polynomial drift
  - [[kriging]] - Kriging with variogram; supports universal Kriging drift
  - [[gp]] - Gaussian process posterior mean

  Extrapolation:

  - [[extrapolation]] - wraps any 1D interpolator to handle values outside the domain"
  (:require [fastmath.interpolation.acm :as acm]
            [fastmath.interpolation.ssj :as ssj]
            [fastmath.interpolation.common :as ic]
            [fastmath.interpolation.linear :as linear]
            [fastmath.interpolation.cubic :as cubic]
            [fastmath.interpolation.akima :as akima]
            [fastmath.interpolation.bspline :as bspline]
            [fastmath.interpolation.bicubic :as bicubic]
            [fastmath.interpolation.loess :as loess]
            [fastmath.interpolation.microsphere :as microsphere]
            [fastmath.interpolation.barycentric :as bc]
            [fastmath.interpolation.shepard :as shepard]
            [fastmath.interpolation.rbf :as rbf]
            [fastmath.interpolation.kriging :as kriging]
            [fastmath.interpolation.gp :as gp]
            [fastmath.interpolation.step :as step]
            [fastmath.interpolation.monotone :as monotone]
            [fastmath.interpolation.sprague :as sprague]
            [fastmath.ml.regression :as reg]
            [fastmath.core :as m]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defmulti interpolation (fn [interpolation-name & _] interpolation-name))

;; 1d

(defmethod interpolation :linear [_ xs ys] (linear/linear xs ys))
(defmethod interpolation :cubic [_ xs ys] (cubic/cubic xs ys))
(defmethod interpolation :monotone [_ xs ys] (monotone/monotone xs ys))
(defmethod interpolation :pchip [_ xs ys] (monotone/pchip xs ys))
(defmethod interpolation :akima
  ([_ xs ys] (akima/akima xs ys))
  ([_ xs ys params] (akima/akima xs ys params)))
(defmethod interpolation :neville [_ xs ys] (acm/neville xs ys))
(defmethod interpolation :divided-difference [_ xs ys] (acm/divided-difference xs ys))
(defmethod interpolation :polynomial [_ xs ys] (bc/polynomial xs ys))
(defmethod interpolation :sprague [_ xs ys] (sprague/sprague xs ys))

(defmethod interpolation :step-before [_ xs ys] (step/step-before xs ys))
(defmethod interpolation :step-after [_ xs ys] (step/step-after xs ys))
(defmethod interpolation :step
  ([_ xs ys] (step/step xs ys))
  ([_ xs ys params] (step/step xs ys params)))

(defmethod interpolation :b-spline
  ([_ xs ys] (bspline/b-spline xs ys))
  ([_ xs ys params] (bspline/b-spline xs ys params)))

(defmethod interpolation :barycentric
  ([_ xs ys] (bc/barycentric xs ys))
  ([_ xs ys params] (bc/barycentric xs ys params)))

;; smoothing

(defmethod interpolation :loess
  ([_ xs ys] (loess/loess xs ys))
  ([_ xs ys params] (loess/loess xs ys params)))

(defmethod interpolation :cubic-smoothing
  ([_ xs ys] (ssj/cubic-smoothing xs ys))
  ([_ xs ys params] (ssj/cubic-smoothing xs ys params)))

;; isotonic regression

(defn- interpolation-method?
  "True when `method` is a function-like value, not a keyword naming a method of [[interpolation]]."
  [method]
  (and (ifn? method) (not (keyword? method))))

(defmethod interpolation :pava
  ([_ xs ys] (interpolation :pava xs ys nil))
  ([_ xs ys {:keys [order weights method]
             :or {order :asc method :linear}}]
   (let [nys (reg/pava ys weights order)]
     (if (interpolation-method? method)
       (method xs nys)
       (interpolation method xs nys)))))

(defmethod interpolation :cir
  ([_ xs ys] (interpolation :cir xs ys nil))
  ([_ xs ys {:keys [order weights method]
             :or {order :asc method :linear}}]
   (let [[nxs nys] (reg/cir xs ys weights order)]
     (if (interpolation-method? method)
       (method nxs nys)
       (interpolation method nxs nys)))))

;; 2d grid

(defmethod interpolation :bilinear [_ xs ys vss] (linear/bilinear xs ys vss))
(defmethod interpolation :bicubic [_ xs ys vss] (bicubic/bicubic xs ys vss))
(defmethod interpolation :cubic-2d [_ xs ys vss] (cubic/cubic-2d xs ys vss))

;; multidim and kernel

(defmethod interpolation :microsphere-projection
  ([_ xss ys] (microsphere/microsphere-projection xss ys))
  ([_ xss ys params] (microsphere/microsphere-projection xss ys params)))

(defmethod interpolation :shepard
  ([_ xss ys] (shepard/shepard xss ys))
  ([_ xss ys params] (shepard/shepard xss ys params)))

(defmethod interpolation :rbf
  ([_ xss ys] (rbf/rbf xss ys))
  ([_ xss ys kernel] (rbf/rbf xss ys kernel))
  ([_ xss ys kernel params] (rbf/rbf xss ys kernel params)))

(defmethod interpolation :kriging
  ([_ xss ys] (kriging/kriging xss ys))
  ([_ xss ys variogram] (kriging/kriging xss ys variogram))
  ([_ xss ys variogram params] (kriging/kriging xss ys variogram params)))

(defmethod interpolation :gp
  ([_ xss ys] (gp/gp xss ys))
  ([_ xss ys kernel] (gp/gp xss ys kernel))
  ([_ xss ys kernel params] (gp/gp xss ys kernel params)))

;; isotonic regression interpolation

(defn pava
  "Creates a 1D isotonic regression interpolator using the Pool-Adjacent-Violators Algorithm (PAVA).

  Fits a monotone step function to `ys` by minimising the (optionally weighted) L2 loss subject to
  the constraint that the result is monotone in the direction given by `:order`. The fitted values
  replace the original `ys`, and a 1D interpolator of type `:method` is then built over the same `xs`.

  Parameters:

  - `xs` - a sequence of x coordinates; the data order is the order of the fit, so `xs` should be sorted in ascending order.
  - `ys` - a sequence of response values corresponding to each x.
  - `opts` - (optional) a map of options:
    - `:order` - monotonicity direction; one of `:asc` / `:increasing` (default), `:desc` / `:decreasing`,
      `:non-decreasing`, or `:non-increasing`. `:desc` and `:decreasing` are the same. The fitted values
      do not depend on the weak or strict variant.
    - `:weights` - a sequence of non-negative observation weights; when `nil`, all points are weighted equally.
    - `:method` - function or keyword of the 1D interpolation method to apply to the isotonic-fitted values (default: `:linear`);
      any keyword accepted by [[interpolation]] that takes `[xs ys]` is valid, a function (or any other function-like value) is called with `xs` and the fitted values.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` with the accepted values for an unknown `:order`.

  See also [[cir]], [[interpolation]]."
  ([xs ys] (interpolation :pava xs ys nil))
  ([xs ys opts] (interpolation :pava xs ys opts)))

(defn cir
  "Creates a 1D isotonic regression interpolator using Centered Isotonic Regression (CIR).

  Unlike PAVA, CIR also adjusts the x coordinates alongside the y values, shrinking each
  constant block to a single representative point. This produces a sparser, centred knot
  set that can yield better-behaved interpolants when combined with smooth methods.
  A 1D interpolator of type `:method` is then built over the reduced `[xs ys]` pair.

  Parameters:

  - `xs` - a sequence of x coordinates; the data order is the order of the fit, so `xs` should be sorted in ascending order.
  - `ys` - a sequence of response values corresponding to each x.
  - `opts` - (optional) a map of options:
    - `:order` - monotonicity direction; one of `:asc` / `:increasing` (default), `:desc` / `:decreasing`,
      `:non-decreasing`, or `:non-increasing`. `:desc` and `:decreasing` are the same. `:asc` pools tied
      neighbouring values into one point, `:non-decreasing` keeps tied values as separate points.
    - `:weights` - a sequence of non-negative observation weights; when `nil`, all points are weighted equally.
    - `:method` - function or keyword of the 1D interpolation method to apply to the reduced knot set (default: `:linear`);
      any keyword accepted by [[interpolation]] that takes `[xs ys]` is valid, a function (or any other function-like value) is called with the reduced `xs` and `ys`.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` with the accepted values for an unknown `:order`.

  See also [[pava]], [[interpolation]]."
  ([xs ys] (interpolation :cir xs ys nil))
  ([xs ys opts] (interpolation :cir xs ys opts)))

;; pure wrapper functions

;; 1d

(defn linear
  "Creates a 1D piecewise linear interpolator.

  Connects consecutive data points with straight line segments. Equal to `numpy.interp` inside the
  domain. Extrapolates linearly beyond the domain boundaries using the slope of the nearest segment.

  Parameters:

  - `xs` - a sequence of x coordinates; at least 2 finite, strictly increasing values.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` for invalid data (different lengths, fewer than 2 points, `xs` not finite or not strictly increasing).

  See also [[cubic]], [[monotone]], [[interpolation]]."
  [xs ys]
  (linear/linear xs ys))

(defn cubic
  "Creates a 1D natural cubic spline interpolator.

  Fits a piecewise cubic polynomial through all data points such that the
  function, its first derivative, and its second derivative are continuous
  everywhere. The natural boundary conditions set the second derivative to
  zero at both endpoints. Equal to `scipy.interpolate.CubicSpline` with `bc_type` natural.

  Outside the domain the cubic of the nearest end segment is continued (as scipy does), it is
  not a linear extrapolation.

  Parameters:

  - `xs` - a sequence of x coordinates; at least 2 finite, strictly increasing values.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` for invalid data (different lengths, fewer than 2 points, `xs` not finite or not strictly increasing).

  See also [[linear]], [[monotone]], [[akima]], [[interpolation]]."
  [xs ys]
  (cubic/cubic xs ys))

(defn monotone
  "Creates a 1D monotone piecewise cubic (Fritsch-Carlson) interpolator.

  Produces a smooth curve that preserves the monotonicity of the data and never leaves the range of the
  neighbouring data values: tangents are the mean of the neighbouring secant slopes, set to zero where the
  secants change sign or vanish, and limited in every interval (including the last) to the circle of
  radius 3 of the Fritsch-Carlson monotonicity region. This avoids the spurious oscillations that can
  appear with standard cubic splines. The curve is C1.

  Outside the domain the cubic of the nearest end segment is continued. For the harmonic mean tangents of
  scipy see [[pchip]].

  Parameters:

  - `xs` - a sequence of x coordinates; at least 2 finite, strictly increasing values.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` for invalid data (different lengths, fewer than 2 points, `xs` not finite or not strictly increasing).

  See also [[pchip]], [[cubic]], [[akima]], [[interpolation]]."
  [xs ys]
  (monotone/monotone xs ys))

(defn pchip
  "Creates a 1D shape-preserving piecewise cubic Hermite interpolator (PCHIP, Fritsch-Butland).

  The tangent at an interior node is the weighted harmonic mean of the neighbouring secant slopes
  (zero where they change sign or vanish), the end tangents come from the non-centered three point formula
  clamped to preserve shape. Equal to `scipy.interpolate.PchipInterpolator`. The curve is C1 and monotone
  where the data is monotone.

  Outside the domain the cubic of the nearest end segment is continued.

  Parameters:

  - `xs` - a sequence of x coordinates; at least 2 finite, strictly increasing values.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` for invalid data (different lengths, fewer than 2 points, `xs` not finite or not strictly increasing).

  See also [[monotone]], [[akima]], [[cubic]], [[interpolation]]."
  [xs ys]
  (monotone/pchip xs ys))

(defn akima
  "Creates a 1D Akima spline interpolator.

  Fits a piecewise cubic Hermite polynomial that is less sensitive to outliers than the natural cubic
  spline. The slope at each knot is a weighted average of the adjacent secant slopes, with weights from the
  differences of the more distant secant slopes; the end slopes use two extra secant slopes obtained by
  linear extrapolation. The curve is C1. Equal to `scipy.interpolate.Akima1DInterpolator`; with
  `:method` `:makima` equal to its `method='makima'`.

  Outside the domain the cubic of the nearest end segment is continued.

  Parameters:

  - `xs` - a sequence of x coordinates; at least 5 finite, strictly increasing values.
  - `ys` - a sequence of y values corresponding to each x.
  - `opts` - (optional) a map of options:
    - `:method` - `:akima` (default) or `:makima` (modified Akima: adds half of the absolute sum of the two secant slopes to each weight).

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` for invalid data or an unknown `:method`.

  See also [[cubic]], [[monotone]], [[pchip]], [[interpolation]]."
  ([xs ys] (akima/akima xs ys))
  ([xs ys opts] (akima/akima xs ys opts)))

(defn neville
  "Creates a 1D polynomial interpolator using Neville's algorithm.

  Constructs the unique polynomial of degree `n-1` passing through all `n`
  data points and evaluates it via the recursive Neville scheme. Suitable for
  small datasets; for large datasets high-degree polynomials may oscillate
  (Runge's phenomenon).

  Parameters:

  - `xs` - a sequence of finite, distinct x coordinates, at least 2, in any order.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` for invalid data (different lengths, fewer than 2 points, `xs` not finite or duplicated).

  See also [[divided-difference]], [[polynomial]], [[barycentric]], [[interpolation]]."
  [xs ys]
  (acm/neville xs ys))

(defn divided-difference
  "Creates a 1D polynomial interpolator using Newton's divided differences.

  Builds the interpolating polynomial in Newton's form using a divided-difference
  table. Mathematically equivalent to [[neville]] but uses a different
  evaluation algorithm. Suitable for small to medium datasets; beware of oscillation
  for large `n`.

  Parameters:

  - `xs` - a sequence of x coordinates; at least 2 finite, strictly increasing values.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` for invalid data (different lengths, fewer than 2 points, `xs` not finite or not strictly increasing).

  See also [[neville]], [[polynomial]], [[interpolation]]."
  [xs ys]
  (acm/divided-difference xs ys))

(defn polynomial
  "Creates a 1D polynomial interpolator in the barycentric form.

  Constructs the unique interpolating polynomial of degree `n-1` through all data points
  (the same as `(barycentric xs ys {:order (dec (count xs))})`). Suitable for
  small datasets; high-degree polynomials may exhibit Runge's phenomenon for large `n`
  or unevenly spaced knots. Far outside the domain the barycentric evaluation is less accurate than
  [[neville]].

  Parameters:

  - `xs` - a sequence of finite, distinct x coordinates, at least 2, in any order.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` for invalid data (different lengths, fewer than 2 points, `xs` not finite or duplicated).

  See also [[neville]], [[divided-difference]], [[barycentric]], [[interpolation]]."
  [xs ys]
  (bc/polynomial xs ys))

(defn sprague
  "Creates a 1D Sprague fifth-order interpolator.

  Fits a fifth-degree polynomial within each interval using six surrounding
  data points (two on each side of the interval). Two synthetic boundary
  points are appended at each end to handle the edges (the CIE rule, equal to the `colour-science`
  `SpragueInterpolator`). Provides higher accuracy than cubic splines for smooth, uniformly sampled
  data: polynomials up to degree 4 are reproduced exactly in the interior intervals, in the boundary
  intervals the accuracy is only of second order (degree 1 is reproduced exactly).

  The spacing of `xs` is assumed to be uniform and it is not checked; non-uniform grids give silently
  wrong values. Outside the domain the polynomial of the nearest end interval is continued.

  Parameters:

  - `xs` - a sequence of x coordinates; at least 6 finite, strictly increasing, uniformly spaced values.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` for invalid data (different lengths, fewer than 6 points, `xs` not finite or not strictly increasing).

  See also [[cubic]], [[akima]], [[interpolation]]."
  [xs ys]
  (sprague/sprague xs ys))

;; step

(defn step-before
  "Creates a left-continuous (step-before) piecewise constant interpolator.

  On `(xs[i], xs[i+1]]` the value is `ys[i+1]`, the y of the nearest knot to the right
  (equal to `scipy.interpolate.interp1d` with kind `'next'` and to the `steps-pre` style of matplotlib).
  The value at a knot is its own y. Outside the domain, clamps to the first or last y value.

  Parameters:

  - `xs` - a sequence of x coordinates; at least 1 finite, strictly increasing value.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` for invalid data (different lengths, no points, `xs` not finite or not strictly increasing).

  See also [[step-after]], [[step]], [[interpolation]]."
  [xs ys]
  (step/step-before xs ys))

(defn step-after
  "Creates a right-continuous (step-after) piecewise constant interpolator.

  On `[xs[i], xs[i+1])` the value is `ys[i]`, the y of the nearest knot to the left
  (equal to `scipy.interpolate.interp1d` with kind `'previous'` and to the `steps-post` style of matplotlib).
  Outside the domain, clamps to the first or last y value.

  Parameters:

  - `xs` - a sequence of x coordinates; at least 1 finite, strictly increasing value.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` for invalid data (different lengths, no points, `xs` not finite or not strictly increasing).

  See also [[step-before]], [[step]], [[interpolation]]."
  [xs ys]
  (step/step-after xs ys))

(defn step
  "Creates a midpoint-blended piecewise constant interpolator.

  Within each interval `[xs[i], xs[i+1]]`, the step occurs at a fractional
  position controlled by `:point`. At the blend point the value switches from
  `ys[i]` to `ys[i+1]` (a tie goes to the right knot: with the default `:point`
  it is equal to `scipy.interpolate.interp1d` with kind `'nearest-up'`). Setting `:point` to `0.0` gives
  step-before behaviour and `1.0` gives step-after behaviour.

  Parameters:

  - `xs` - a sequence of x coordinates; at least 2 finite, strictly increasing values.
  - `ys` - a sequence of y values corresponding to each x.
  - `opts` - (optional) a map of options:
    - `:point` - fractional position within each interval where the step occurs (default: `0.5`).

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` for invalid data (different lengths, fewer than 2 points, `xs` not finite or not strictly increasing).

  See also [[step-before]], [[step-after]], [[interpolation]]."
  ([xs ys] (step/step xs ys))
  ([xs ys opts] (step/step xs ys opts)))

;; smoothing

(defn loess
  "Creates a 1D LOESS (locally weighted scatterplot smoothing) interpolator.

  Smooths the data using a locally weighted linear regression over the `k` nearest neighbours of every
  node with tricube weights, optionally repeated with bisquare robustness weights from residuals scaled by
  6 times their median, then fits a natural cubic spline through the smoothed points. The smoothed values equal
  `statsmodels` `lowess` with `frac` equal to `:bandwidth`, `it` equal to `:iters` and `delta` 0. Useful when the
  data contains noise and strict interpolation through every point is not desired.

  Parameters:

  - `xs` - a sequence of x coordinates; at least 2 finite, strictly increasing values.
  - `ys` - a sequence of y values corresponding to each x.
  - `opts` - (optional) a map of options:
    - `:bandwidth` - fraction of points used for each local fit, in `(0, 1]` (default: `0.3`); `k = floor(bandwidth * n)` has to be at least 2, so the default needs 7 points. `k = 2` returns the data unchanged.
    - `:iters` - number of robustness iterations (default: `2`).
    - `:accuracy` - the robustness iterations stop when the median absolute residual is below it (default: `1.0e-12`).
    - `:weights` - a sequence of non-negative observation weights, multiplied into the local weights; when `nil`, all points are weighted equally.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` for invalid data or options.

  See also [[cubic-smoothing]], [[cubic]], [[interpolation]]."
  ([xs ys] (loess/loess xs ys))
  ([xs ys opts] (loess/loess xs ys opts)))

(defn cubic-smoothing
  "Creates a 1D smoothing cubic spline (Reinsch).

  Fits a cubic spline that balances fidelity to the data against smoothness,
  controlled by the parameter `rho`. Equal to `scipy.interpolate.make_smoothing_spline` with
  `lam = (1 - rho) / rho`. When `rho` is `1.0` the spline is the natural interpolating cubic spline;
  smaller values produce smoother curves that do not pass through every point; `0.0` is accepted and gives
  the least squares straight line. Outside the domain the spline continues linearly with the end tangent.

  Parameters:

  - `xs` - a sequence of x coordinates; at least 2 finite, strictly increasing values.
  - `ys` - a sequence of y values corresponding to each x.
  - `opts` - (optional) a map of options:
    - `:rho` - smoothing parameter in `[0, 1]`; `1.0` gives exact interpolation (default: `1.0`).
    - `:weights` - a sequence of positive finite observation weights (scipy's weights are their square roots); when `nil`, all points are weighted equally.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` for invalid data or options.

  See also [[loess]], [[cubic]], [[interpolation]]."
  ([xs ys] (ssj/cubic-smoothing xs ys))
  ([xs ys opts] (ssj/cubic-smoothing xs ys opts)))

;; barycentric / b-spline

(defn barycentric
  "Creates a 1D barycentric rational interpolator (Floater-Hormann).

  Uses the barycentric form of the Floater-Hormann blend of local polynomial interpolants, which is
  numerically stable, has no poles on the real line and avoids solving a linear system. The `order`
  parameter is the degree of the blended local polynomials: `0` is Berrut's interpolant, `n-1` is the
  interpolating polynomial. The nodes are sorted internally, so the order of the data does not matter.

  Parameters:

  - `xs` - a sequence of finite, distinct x coordinates, in any order.
  - `ys` - a sequence of y values corresponding to each x.
  - `opts` - (optional) a map of options:
    - `:order` - integer in `[0, (count xs) - 1]` (default: `1`).

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` for invalid data, duplicated nodes or an invalid `:order`.

  See also [[polynomial]], [[neville]], [[b-spline]], [[interpolation]]."
  ([xs ys] (bc/barycentric xs ys))
  ([xs ys opts] (bc/barycentric xs ys opts)))

(defn b-spline
  "Creates a 1D interpolating B-spline.

  The spline of degree `:degree` passes through all data points. It is parametrised in x and equal to
  `scipy.interpolate.make_interp_spline`. The default knots are not-a-knot for an odd degree and the
  midpoints between data points for an even degree. Outside the data range the end polynomial is continued.

  Parameters:

  - `xs` - a sequence of x coordinates; at least `degree + 1` finite, strictly increasing values.
  - `ys` - a sequence of y values corresponding to each x.
  - `opts` - (optional) a map of options:
    - `:degree` - integer degree of the spline, at least `1` (default: `3`).
    - `:knots` - an optional full knot vector of `(count xs) + degree + 1` non-decreasing values satisfying the Schoenberg-Whitney conditions.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.
  Throws `ex-info` for invalid data or options.

  The former control polygon form and the options `:clamped?` and `:hp1` do not exist anymore.

  See also [[cubic]], [[barycentric]], [[polynomial]], [[interpolation]]."
  ([xs ys] (bspline/b-spline xs ys))
  ([xs ys opts] (bspline/b-spline xs ys opts)))

;; 2d grid

(defn bilinear
  "Creates a 2D bilinear grid interpolator.

  Linearly interpolates within each rectangular cell of a regular or irregular
  grid. Given a query point `(x, y)`, it locates the enclosing grid cell and
  blends the four corner values using bilinear weights. Equal to scipy `RegularGridInterpolator` with the
  linear method. Outside the grid the nearest border cell is continued linearly.

  Parameters:

  - `xs` - a sequence of x grid coordinates; at least 2 finite, strictly increasing values.
  - `ys` - a sequence of y grid coordinates; at least 2 finite, strictly increasing values.
  - `vss` - a 2D sequence of values; `(get-in vss [i j])` is the value at `(xs[i], ys[j])`; `(count xs)` rows of `(count ys)` values.

  Returns a function `f` accepting either `(f x y)` or `(f [x y])` that returns
  the interpolated `double` value at `(x, y)`.
  Throws `ex-info` for invalid axes or a wrong shape of `vss`.

  See also [[bicubic]], [[cubic-2d]], [[interpolation]]."
  [xs ys vss]
  (linear/bilinear xs ys vss))

(defn bicubic
  "Creates a 2D bicubic Hermite grid interpolator.

  Fits a bicubic Hermite patch within each rectangular cell of the grid; the surface is C1 and
  reproduces `a + b x + c y + d x y` exactly. The derivatives at the nodes are central differences between
  the neighbours at interior nodes (as in Apache Commons Math), second order one-sided three point
  differences at boundary nodes, the secant for an axis with two nodes (a 2x2 grid is bilinear).
  Outside the grid the nearest border patch is continued.

  Parameters:

  - `xs` - a sequence of x grid coordinates; at least 2 finite, strictly increasing values.
  - `ys` - a sequence of y grid coordinates; at least 2 finite, strictly increasing values.
  - `vss` - a 2D sequence of values; `(get-in vss [i j])` is the value at `(xs[i], ys[j])`; `(count xs)` rows of `(count ys)` values.

  Returns a function `f` accepting either `(f x y)` or `(f [x y])` that returns
  the interpolated `double` value at `(x, y)`.
  Throws `ex-info` for invalid axes or a wrong shape of `vss`.

  See also [[bilinear]], [[cubic-2d]], [[interpolation]]."
  [xs ys vss]
  (bicubic/bicubic xs ys vss))

(defn cubic-2d
  "Creates a 2D cubic spline grid interpolator.

  Applies 1D natural cubic splines first along each row (the `ys` axis) and then
  along the `xs` axis, giving a smooth (C2) surface that passes through all grid points. It is a
  natural tensor-product spline, not the not-a-knot spline of scipy `RectBivariateSpline`. The surface
  continues outside the grid with the end cubics.

  Parameters:

  - `xs` - a sequence of x grid coordinates; at least 2 finite, strictly increasing values.
  - `ys` - a sequence of y grid coordinates; at least 2 finite, strictly increasing values.
  - `vss` - a 2D sequence of values; `(get-in vss [i j])` is the value at `(xs[i], ys[j])`; `(count xs)` rows of `(count ys)` values.

  Returns a function `f` accepting either `(f x y)` or `(f [x y])` that returns
  the interpolated `double` value at `(x, y)`.
  Throws `ex-info` for invalid axes or a wrong shape of `vss`.

  See also [[bilinear]], [[bicubic]], [[interpolation]]."
  [xs ys vss]
  (cubic/cubic-2d xs ys vss))

;; multidim

(defn microsphere-projection
  "Creates a multivariate interpolator using the microsphere projection method (Apache Commons Math).

  A microsphere with facets is placed in the query point and lit by the nodes, with the brightness falling with the distance to the power of `:exponent`; the value is the illumination weighted mean of the values. It passes through the nodes, returns `:background` where the sphere is too dark, and is usually less accurate than [[rbf]] on smooth data.

  Parameters:

  - `xss` - a sequence of input points; numbers or one element sequences (1D) or sequences of coordinates.
  - `ys` - a sequence of output values corresponding to each point in `xss`.
  - `opts` - (optional) a map of options (a `nil` value is a missing key):
    - `:elements` - number of facets of the microsphere (default: `max(5, 4+dims, n/5)`).
    - `:exponent` - distance weighting exponent, not negative (default: `1.0`).
    - `:max-dark-friction` - maximum fraction of dark facets, above it `:background` is returned (default: `0.9`).
    - `:dark-threshold` - illumination below which a facet is dark (default: `0.01`).
    - `:background` - value returned where the sphere is dark (default: `0.0`).
    - `:no-interpolation-tolerance` - distance below which the value of the node is returned (default: `1.0e-6`).
    - `:rng` - random number generator of the facets (see `fastmath.random/rng`); seed it to get the same interpolator on every construction. In two dimensions the facets are deterministic and `:rng` is not used.

  Returns a function `f` of a query point (a number or a one element sequence in 1D, a sequence of coordinates otherwise; in 2D also `(f x y)`) returning a `double`. Throws `ex-info` for invalid input.

  See also [[shepard]], [[rbf]], [[interpolation]]."
  ([xss ys] (microsphere/microsphere-projection xss ys))
  ([xss ys opts] (microsphere/microsphere-projection xss ys opts)))

(defn shepard
  "Creates a multivariate Shepard (inverse distance weighting) interpolator.

  The value at a query point is the mean of the training values weighted by `1 / d^p`, where `d` is the distance to the node. The interpolator passes through the nodes and its result is always a weighted mean of `ys` (far from the data it tends to the mean of all the values). A query is a node when its distance to the node is exactly zero (a long hits a double node, a `Vec2` a vector node); for repeated nodes the last value wins. The weights are scaled by the nearest distance, so the result is never `Inf` or `NaN` for very small distances.

  Parameters:

  - `xss` - a sequence of input points (at least one); numbers (1D) or sequences of coordinates of an equal dimension.
  - `ys` - a sequence of output values corresponding to each point in `xss`.
  - `opts` - (optional) a map of options:
    - `:p` - the power of the distance, finite and not negative; larger values give a more local influence of the nearest nodes; `0` gives the mean of all values (default: `2.0`).
    - `:distance` - a two-argument distance function (default: the Euclidean distance).

  Returns a function `f` of a query point returning a `double`. Throws `ex-info` for different lengths of `xss` and `ys`, no points, an invalid `:p` or, with the default distance, a wrong dimension.

  See also [[microsphere-projection]], [[rbf]], [[interpolation]]."
  ([xss ys] (shepard/shepard xss ys))
  ([xss ys opts] (shepard/shepard xss ys opts)))

(defn rbf
  "Creates a multivariate Radial Basis Function (RBF) interpolator.

  Solves a linear system whose kernel matrix is built from a radial basis function of the pairwise distances between the nodes. A polynomial drift can be added via `:polynomial-terms`.

  Parameters:

  - `xss` - a sequence of input points (at least one); numbers (1D) or sequences of coordinates.
  - `ys` - a sequence of output values corresponding to each point in `xss`.
  - `kernel` - (optional) a radial kernel function `r -> double` (see `fastmath.kernel/rbf`), its `:scale` sets the length scale; defaults to the gaussian.
  - `opts` - (optional) a map of options:
    - `:kscale` - multiplier of all kernel evaluations; it cancels without regularization (default: `1.0`).
    - `:lambda` - regularization strength, finite and not negative; `0.0` is the exact interpolation (default: `0.0`).
    - `:regularization` - `:normal` (default, ridge on the normal equations `(A^T A + lambda I) w = A^T y`, squares the condition number, `lambda` acts as `lambda / kscale^2` without polynomial terms) or `:ridge` (kernel ridge `(K + lambda I) w = y`, the `smoothing` of the scipy `RBFInterpolator`, the polynomial block is not regularized, `lambda` acts as `lambda / kscale`).
    - `:distance` - a two-argument distance function (default: the Euclidean distance).
    - `:polynomial-terms` - a function `x -> seq` returning the polynomial basis at `x`; augments the system with a polynomial.

  Conditionally positive definite kernels (thin plate splines, cubic and other radial powers, linear, generalized multiquadratic) need `:polynomial-terms` of the matching degree for a guaranteed solution. The system is solved by LU; a (numerically) singular one, for example with repeated nodes, silently falls back to the minimum norm SVD solution. Smooth kernels on dense nodes are ill-conditioned: use a small `:lambda` with `:ridge`.

  Returns a function `f` of a query point returning a `double`. Throws `ex-info` for different lengths, no points, an invalid `:lambda` or `:regularization`, or a wrong query dimension.

  See also [[shepard]], [[kriging]], [[gp]], [[interpolation]]."
  ([xss ys] (rbf/rbf xss ys))
  ([xss ys kernel] (rbf/rbf xss ys kernel))
  ([xss ys kernel opts] (rbf/rbf xss ys kernel opts)))

(defn kriging
  "Creates a multivariate Kriging interpolator, ordinary by default and universal with `:polynomial-terms`.

  Solves the kriging system built from a semivariogram. The variogram has to be valid: zero at lag zero (a nugget is the limit at lag 0+) and conditionally negative definite. If none is supplied, the empirical variogram of the data is estimated and a gaussian model fitted to it by random search (use `:rng` to make it reproducible, `fastmath.kernel.variogram` to check or fit a model yourself).

  Parameters:

  - `xss` - a sequence of input points (at least two when the variogram is estimated); numbers (1D) or sequences of coordinates.
  - `ys` - a sequence of output values corresponding to each point in `xss`.
  - `variogram` - (optional) a semivariogram function `lag -> double`; estimated and fitted when `nil`.
  - `opts` - (optional) a map of options:
    - `:error` - standard deviation of the measurement error, a number or one per node (default: `0.0`); `error^2` is subtracted from the diagonal of the variogram matrix (the gstat `Err` is this variance), so the result smooths the data.
    - `:distance` - a two-argument distance function (default: the Euclidean distance); also used by the estimation of the default variogram.
    - `:polynomial-terms` - a function `x -> seq` returning the drift basis (default: `(constantly [1.0])`; `nil` removes the drift).
    - `:rng`, `:N` - random number generator and number of starting points of the default variogram fit.

  Repeated nodes give the mean of their values (silent SVD fallback). The kriging variance is not available.

  Returns a function `f` of a query point returning a `double`. Throws `ex-info` for different lengths, too few points, an invalid `:error` (also a sequence of a wrong length) or a wrong query dimension.

  See also [[rbf]], [[gp]], [[interpolation]]."
  ([xss ys] (kriging/kriging xss ys))
  ([xss ys variogram] (kriging/kriging xss ys variogram))
  ([xss ys variogram opts] (kriging/kriging xss ys variogram opts)))

(defn gp
  "Creates a multivariate Gaussian Process (GP) interpolator returning the posterior mean.

  For the standard deviation, the log marginal likelihood and samples use `fastmath.interpolation.gp/gaussian-process`, `predict`, `predict-all`, `L`, `prior-samples` and `posterior-samples`.

  Parameters:

  - `xss` - a sequence of input points (at least one); numbers (1D) or sequences of coordinates of an equal dimension.
  - `ys` - a sequence of output values corresponding to each point in `xss`.
  - `kernel` - (optional) a covariance kernel `(x1 x2) -> double` (see `fastmath.kernel/kernel`); defaults to `(fastmath.kernel/kernel :gaussian)` (sigma `1.0`).
  - `opts` - (optional) a map of options:
    - `:kscale` - positive multiplier of the kernel, the signal variance (default: `1.0`).
    - `:noise` - variance added to the diagonal of the covariance matrix, in the units of the normalised values (default: `1.0e-8`); a larger value smooths the data.
    - `:normalize?` - standardise `ys` to zero mean and unit sample standard deviation (n - 1) before the fit (default: `false`); a constant `ys` uses the standard deviation `1`.

  Returns a function `f` of a query point returning the posterior mean as a `double`. Throws `ex-info` for different lengths, no points, ragged points, an invalid `:kscale` or `:noise`, a covariance matrix which is not positive definite (increase `:noise`) or a query of a wrong dimension.

  See also [[kriging]], [[rbf]], [[interpolation]]."
  ([xss ys] (gp/gp xss ys))
  ([xss ys kernel] (gp/gp xss ys kernel))
  ([xss ys kernel opts] (gp/gp xss ys kernel opts)))

;;

(defn- extrapolation-sides
  "Creates the function returning the first of `sides` below `start`, the second above `end`, NaN for a NaN `x` and delegating to `interpolator` inside `[start, end]`."
  [interpolator ^double start ^double end sides]
  (let [left (double (first sides))
        right (double (second sides))]
    (fn ^double [^double x]
      (cond (m/nan? x) ##NaN
            (m/< x start) left
            (m/> x end) right
            :else (interpolator x)))))

(defn- extrapolation-side-value
  "Returns the value of one side of the extrapolation policy as a double, throws `ex-info` when it is not a number."
  ^double [method side value]
  (if (number? value)
    (double value)
    (ic/invalid-input! :extrapolation :method
                       (str "the " (name side) " value of the extrapolation method has to be a number, got " (pr-str value))
                       {:extrapolation method :side side :value value})))

(defn extrapolation
  "Wraps a 1D interpolator with a boundary policy that controls behaviour outside `[start, end]`.

  Inside the domain (the boundaries included) the wrapped function always delegates to `interpolator`.
  Outside it, the chosen `method` determines the returned value. When called with three arguments `method`
  defaults to `:skip`. A `NaN` argument gives `NaN` with every policy except `:error`, which throws.

  Parameters:

  - `interpolator` - a 1D function `(fn ^double [^double x] ...)` produced by any 1D interpolation constructor.
  - `start` - the lower boundary of the valid domain as a `double`.
  - `end` - the upper boundary of the valid domain as a `double`, greater than `start`.
  - `method` - the extrapolation policy; one of:
    - `:skip` (default) - returns `interpolator` itself: all `x` values are passed directly to the interpolator, disabling boundary enforcement.
    - `:constant` - returns `(interpolator start)` for `x < start` and `(interpolator end)` for `x > end`.
    - `:zero` - returns `0.0` for any `x` outside `[start, end]`.
    - `:error` - throws `IndexOutOfBoundsException` for any `x` outside `[start, end]`.
    - a number - returns that constant value for any `x` outside `[start, end]`.
    - a two-element sequential `[left right]` of numbers - returns `left` for `x < start` and `right` for `x > end`.
    - a map `{:left left-val :right right-val}` of numbers - same as the two-element vector form.

  Returns a new function `f` with the same signature as `interpolator`,
  i.e. `(fn ^double [^double x] ...)`.
  Throws `ex-info` (`:type` `:invalid-input`, `:method` `:extrapolation`) when `start` is not lower than `end`
  (also for `NaN`), for an unknown keyword, `nil` or any other `method`, and for a vector or a map without
  two numeric sides.

  See also [[linear]], [[interpolation]]."
  ([interpolator ^double start ^double end]
   (extrapolation interpolator :skip start end))
  ([interpolator method ^double start ^double end]
   (when-not (m/< start end)
     (ic/invalid-input! :extrapolation :domain
                        (str "start has to be lower than end, got " start " and " end)
                        {:start start :end end}))
   (cond
     (sequential? method) (if (= 2 (count method))
                            (let [[left right] method]
                              (extrapolation-sides interpolator start end
                                                   [(extrapolation-side-value method :left left)
                                                    (extrapolation-side-value method :right right)]))
                            (ic/invalid-input! :extrapolation :method
                                               (str "a sequential extrapolation method needs exactly two values [left right], got " (count method))
                                               {:extrapolation method}))
     (number? method) (let [v (double method)]
                        (extrapolation-sides interpolator start end [v v]))
     (map? method) (extrapolation-sides interpolator start end
                                        [(extrapolation-side-value method :left (:left method))
                                         (extrapolation-side-value method :right (:right method))])
     (keyword? method) (case method
                         :skip interpolator
                         :constant (extrapolation-sides interpolator start end
                                                        [(interpolator start) (interpolator end)])
                         :zero (extrapolation-sides interpolator start end [0.0 0.0])
                         :error (fn ^double [^double x]
                                  (if (m/between? start end x)
                                    (interpolator x)
                                    (throw (IndexOutOfBoundsException.
                                            (str "x=" x " not in [" start ", " end "] range")))))
                         (ic/invalid-input! :extrapolation :method
                                            (str "unknown extrapolation method " method)
                                            {:extrapolation method}))
     :else (ic/invalid-input! :extrapolation :method
                              (str "unknown extrapolation method " (pr-str method))
                              {:extrapolation method}))))
