(ns fastmath.interpolation
  "Interpolation and extrapolation functions for 1D, 2D grid, and multivariate data.

  Every interpolator is constructed from a set of known data points and returns a Clojure
  function that evaluates the fitted surface at any query point. All constructors are available
  both as pure named functions and via the unified [[interpolation]] multimethod, which selects
  the implementation from a keyword.

  1D interpolators (inputs `xs`, `ys`):

  - [[linear]] - piecewise linear; extrapolates beyond domain
  - [[cubic]] - natural cubic spline; C2-continuous
  - [[monotone]] - monotonicity-preserving cubic (Fritsch-Carlson)
  - [[akima]] - Akima spline; robust to outliers; extrapolates
  - [[sprague]] - fifth-order Sprague; accurate for smooth uniform data
  - [[neville]] - polynomial via Neville's algorithm
  - [[divided-difference]] - polynomial via Newton's divided differences
  - [[polynomial]] - polynomial via SSJ `PolInterp`
  - [[barycentric]] - numerically stable barycentric rational interpolation
  - [[b-spline]] - B-spline interpolation or least-squares approximation
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
  - [[bicubic]] - bicubic interpolation on a rectangular grid
  - [[cubic-2d]] - 2D cubic spline on a rectangular grid

  Multivariate interpolators (inputs `xss`, `ys`; work in any number of dimensions):

  - [[microsphere-projection]] - ACM microsphere projection
  - [[shepard]] - inverse distance weighting
  - [[rbf]] - radial basis function with optional regularisation and polynomial drift
  - [[kriging]] - Kriging with variogram; supports universal Kriging drift
  - [[gp]] - Gaussian process posterior mean

  Extrapolation:

  - [[extrapolation]] - wraps any 1D interpolator to handle values outside the domain"
  (:require [fastmath.interpolation.acm :as acm]
            [fastmath.interpolation.ssj :as ssj]
            [fastmath.interpolation.linear :as linear]
            [fastmath.interpolation.cubic :as cubic]
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
(defmethod interpolation :akima [_ xs ys] (acm/akima xs ys))
(defmethod interpolation :neville [_ xs ys] (acm/neville xs ys))
(defmethod interpolation :divided-difference [_ xs ys] (acm/divided-difference xs ys))
(defmethod interpolation :polynomial [_ xs ys] (ssj/polynomial xs ys))
(defmethod interpolation :sprague [_ xs ys] (sprague/sprague xs ys))

(defmethod interpolation :step-before [_ xs ys] (step/step-before xs ys))
(defmethod interpolation :step-after [_ xs ys] (step/step-after xs ys))
(defmethod interpolation :step
  ([_ xs ys] (step/step xs ys))
  ([_ xs ys params] (step/step xs ys params)))

(defmethod interpolation :b-spline
  ([_ xs ys] (ssj/b-spline xs ys))
  ([_ xs ys params] (ssj/b-spline xs ys params)))

(defmethod interpolation :barycentric
  ([_ xs ys] (bc/barycentric xs ys))
  ([_ xs ys params] (bc/barycentric xs ys params)))

;; smoothing

(defmethod interpolation :loess
  ([_ xs ys] (acm/loess xs ys))
  ([_ xs ys params] (acm/loess xs ys params)))

(defmethod interpolation :cubic-smoothing
  ([_ xs ys] (ssj/cubic-smoothing xs ys))
  ([_ xs ys params] (ssj/cubic-smoothing xs ys params)))

;; isotonic regression

(defmethod interpolation :pava
  ([_ xs ys] (interpolation :pava xs ys nil))
  ([_ xs ys {:keys [order weights method]
             :or {order :asc method :linear}}]
   (let [nys (reg/pava ys weights order)]
     (if (fn? method)
       (method xs nys)
       (interpolation method xs nys)))))

(defmethod interpolation :cir
  ([_ xs ys] (interpolation :cir xs ys nil))
  ([_ xs ys {:keys [order weights method]
             :or {order :asc method :linear}}]
   (let [[nxs nys] (reg/cir xs ys weights order)]
     (if (fn? method)
       (method nxs nys)
       (interpolation method nxs nys)))))

;; 2d grid

(defmethod interpolation :bilinear [_ xs ys vss] (linear/bilinear xs ys vss))
(defmethod interpolation :bicubic [_ xs ys vss] (acm/bicubic xs ys vss))
(defmethod interpolation :cubic-2d [_ xs ys vss] (cubic/cubic-2d xs ys vss))

;; multidim and kernel

(defmethod interpolation :microsphere-projection
  ([_ xss ys] (acm/microsphere-projection xss ys))
  ([_ xss ys params] (acm/microsphere-projection xss ys params)))

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

  - `xs` - a sequence of x coordinates; must be sorted in ascending order.
  - `ys` - a sequence of response values corresponding to each x.
  - `opts` - (optional) a map of options:
    - `:order` - monotonicity direction; one of `:asc` / `:increasing` (default), `:desc` / `:decreasing`,
      `:non-decreasing`, or `:non-increasing`.
    - `:weights` - a sequence of non-negative observation weights; when `nil`, all points are weighted equally.
    - `:method` - function or keyword of the 1D interpolation method to apply to the isotonic-fitted values (default: `:linear`);
      any keyword accepted by [[interpolation]] that takes `[xs ys]` is valid.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

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

  - `xs` - a sequence of x coordinates; must be sorted in ascending order.
  - `ys` - a sequence of response values corresponding to each x.
  - `opts` - (optional) a map of options:
    - `:order` - monotonicity direction; one of `:asc` / `:increasing` (default), `:desc` / `:decreasing`,
      `:non-decreasing`, or `:non-increasing`.
    - `:weights` - a sequence of non-negative observation weights; when `nil`, all points are weighted equally.
    - `:method` - function or keyword of the 1D interpolation method to apply to the reduced knot set (default: `:linear`);
      any keyword accepted by [[interpolation]] that takes `[xs ys]` is valid.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

  See also [[pava]], [[interpolation]]."
  ([xs ys] (interpolation :cir xs ys nil))
  ([xs ys opts] (interpolation :cir xs ys opts)))

;; pure wrapper functions

;; 1d

(defn linear
  "Creates a 1D piecewise linear interpolator.

  Connects consecutive data points with straight line segments. Extrapolates
  linearly beyond the domain boundaries using the slope of the nearest segment.

  Parameters:

  - `xs` - a sequence of x coordinates; must be sorted in ascending order.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

  See also [[cubic]], [[monotone]], [[interpolation]]."
  [xs ys]
  (linear/linear xs ys))

(defn cubic
  "Creates a 1D natural cubic spline interpolator.

  Fits a piecewise cubic polynomial through all data points such that the
  function, its first derivative, and its second derivative are continuous
  everywhere. The natural boundary conditions set the second derivative to
  zero at both endpoints.

  Parameters:

  - `xs` - a sequence of x coordinates; must be sorted in ascending order.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

  See also [[linear]], [[monotone]], [[akima]], [[interpolation]]."
  [xs ys]
  (cubic/cubic xs ys))

(defn monotone
  "Creates a 1D monotone piecewise cubic (Fritsch-Carlson) interpolator.

  Produces a smooth curve that preserves the monotonicity of the data:
  if the input values are monotonically increasing (or decreasing) over
  an interval, the interpolant will be too. This avoids the spurious
  oscillations that can appear with standard cubic splines.

  Parameters:

  - `xs` - a sequence of x coordinates; must be sorted in ascending order.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

  See also [[cubic]], [[akima]], [[interpolation]]."
  [xs ys]
  (monotone/monotone xs ys))

(defn akima
  "Creates a 1D Akima spline interpolator.

  Fits a piecewise cubic polynomial that is less sensitive to outliers than
  the natural cubic spline. The slope at each knot is computed from a weighted
  average of adjacent slopes, so a few extreme data points do not distort the
  whole curve. Extrapolation beyond the domain is supported.

  Parameters:

  - `xs` - a sequence of x coordinates; must be sorted in ascending order; requires at least 5 points.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

  See also [[cubic]], [[monotone]], [[interpolation]]."
  [xs ys]
  (acm/akima xs ys))

(defn neville
  "Creates a 1D polynomial interpolator using Neville's algorithm.

  Constructs the unique polynomial of degree `n-1` passing through all `n`
  data points and evaluates it via the recursive Neville scheme. Suitable for
  small datasets; for large datasets high-degree polynomials may oscillate
  (Runge's phenomenon).

  Parameters:

  - `xs` - a sequence of x coordinates; must be distinct.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

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

  - `xs` - a sequence of x coordinates; must be distinct.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

  See also [[neville]], [[polynomial]], [[interpolation]]."
  [xs ys]
  (acm/divided-difference xs ys))

(defn polynomial
  "Creates a 1D polynomial interpolator using the SSJ `PolInterp` implementation.

  Constructs the unique interpolating polynomial through all data points.
  Suitable for small datasets; high-degree polynomials may exhibit Runge's
  phenomenon for large `n` or unevenly spaced knots.

  Parameters:

  - `xs` - a sequence of x coordinates; must be distinct.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

  See also [[neville]], [[divided-difference]], [[b-spline]], [[interpolation]]."
  [xs ys]
  (ssj/polynomial xs ys))

(defn sprague
  "Creates a 1D Sprague fifth-order interpolator.

  Fits a fifth-degree polynomial within each interval using six surrounding
  data points (two on each side of the interval). Two synthetic boundary
  points are appended at each end to handle the edges. Provides higher
  accuracy than cubic splines for smooth, uniformly sampled data.

  Parameters:

  - `xs` - a sequence of x coordinates; must be sorted in ascending order; requires at least 6 points.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

  See also [[cubic]], [[akima]], [[interpolation]]."
  [xs ys]
  (sprague/sprague xs ys))

;; step

(defn step-before
  "Creates a left-continuous (step-before) piecewise constant interpolator.

  Returns the y value of the nearest knot to the left of `x` (i.e., the
  largest `xs[i]` that is less than or equal to `x`). Outside the domain,
  clamps to the first or last y value.

  Parameters:

  - `xs` - a sequence of x coordinates; must be sorted in ascending order.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

  See also [[step-after]], [[step]], [[interpolation]]."
  [xs ys]
  (step/step-before xs ys))

(defn step-after
  "Creates a right-continuous (step-after) piecewise constant interpolator.

  Returns the y value of the nearest knot to the right of `x` (i.e., the
  smallest `xs[i]` that is greater than or equal to `x`). Outside the domain,
  clamps to the first or last y value.

  Parameters:

  - `xs` - a sequence of x coordinates; must be sorted in ascending order.
  - `ys` - a sequence of y values corresponding to each x.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

  See also [[step-before]], [[step]], [[interpolation]]."
  [xs ys]
  (step/step-after xs ys))

(defn step
  "Creates a midpoint-blended piecewise constant interpolator.

  Within each interval `[xs[i], xs[i+1]]`, the step occurs at a fractional
  position controlled by `:point`. At the blend point the value switches from
  `ys[i]` to `ys[i+1]`. Setting `:point` to `0.0` gives step-before behaviour
  and `1.0` gives step-after behaviour.

  Parameters:

  - `xs` - a sequence of x coordinates; must be sorted in ascending order.
  - `ys` - a sequence of y values corresponding to each x.
  - `opts` - (optional) a map of options:
    - `:point` - fractional position within each interval where the step occurs (default: `0.5`).

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

  See also [[step-before]], [[step-after]], [[interpolation]]."
  ([xs ys] (step/step xs ys))
  ([xs ys opts] (step/step xs ys opts)))

;; smoothing

(defn loess
  "Creates a 1D LOESS (locally weighted scatterplot smoothing) interpolator.

  Smooths the data using a locally weighted regression, then fits a cubic spline
  through the smoothed points. Useful when the data contains noise and strict
  interpolation through every point is not desired.

  Parameters:

  - `xs` - a sequence of x coordinates; must be sorted in ascending order.
  - `ys` - a sequence of y values corresponding to each x.
  - `opts` - (optional) a map of options:
    - `:bandwidth` - fraction of points used for each local fit (default: ACM `LoessInterpolator` default, `0.3`).
    - `:iters` - number of robustness iterations (default: ACM default, `2`).
    - `:accuracy` - convergence criterion for robustness iterations (default: ACM default, `1.0e-12`).
    - `:weights` - a sequence of observation weights; when `nil`, all points are weighted equally.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

  See also [[cubic-smoothing]], [[cubic]], [[interpolation]]."
  ([xs ys] (acm/loess xs ys))
  ([xs ys opts] (acm/loess xs ys opts)))

(defn cubic-smoothing
  "Creates a 1D smoothing cubic spline interpolator using the SSJ library.

  Fits a cubic spline that balances fidelity to the data against smoothness,
  controlled by the roughness penalty parameter `rho`. When `rho` is `1.0` the
  spline interpolates exactly; smaller values produce smoother curves that do
  not pass through every point.

  Parameters:

  - `xs` - a sequence of x coordinates; must be sorted in ascending order.
  - `ys` - a sequence of y values corresponding to each x.
  - `opts` - (optional) a map of options:
    - `:rho` - roughness penalty in `(0, 1]`; `1.0` gives exact interpolation (default: `1.0`).
    - `:weights` - a sequence of positive observation weights; when `nil`, all points are weighted equally.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

  See also [[loess]], [[cubic]], [[interpolation]]."
  ([xs ys] (ssj/cubic-smoothing xs ys))
  ([xs ys opts] (ssj/cubic-smoothing xs ys opts)))

;; barycentric / b-spline

(defn barycentric
  "Creates a 1D barycentric rational interpolator.

  Uses the barycentric form of polynomial interpolation (Berrut and Trefethen),
  which is numerically stable and avoids solving a linear system. The `order`
  parameter controls the degree of the local polynomial used; higher values
  yield smoother curves but require more data points.

  Parameters:

  - `xs` - a sequence of x coordinates; must be sorted in ascending order.
  - `ys` - a sequence of y values corresponding to each x.
  - `opts` - (optional) a map of options:
    - `:order` - local polynomial degree; must satisfy `0 <= order < (count xs)` (default: `1`).

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

  See also [[neville]], [[b-spline]], [[interpolation]]."
  ([xs ys] (bc/barycentric xs ys))
  ([xs ys opts] (bc/barycentric xs ys opts)))

(defn b-spline
  "Creates a 1D B-spline interpolator or approximator using the SSJ library.

  Supports interpolation (exact fit) and approximation (least-squares fit)
  with configurable knots and degree. When `:clamped?` is `true`, clamped
  knots are computed automatically from the degree; otherwise the knots
  default to a uniform distribution of degree `N-1`.

  Parameters:

  - `xs` - a sequence of x coordinates.
  - `ys` - a sequence of y values corresponding to each x.
  - `opts` - (optional) a map of options:
    - `:knots` - an explicit sequence of knot positions; when provided, `:degree` and `:clamped?` are ignored.
    - `:degree` - polynomial degree of each spline segment; defaults to `3` when `:clamped?` is `true`, or `N-1` otherwise.
    - `:clamped?` - when `true`, uses clamped (pinned) knots computed from the degree (default: `false`).
    - `:hp1` - when set and `:clamped?` is `true`, creates an approximating B-spline with `hp1` basis functions; must satisfy `degree < hp1 <= N`.

  Returns a function `f` such that `(f x)` returns the interpolated `double` value at `x`.

  See also [[barycentric]], [[polynomial]], [[interpolation]]."
  ([xs ys] (ssj/b-spline xs ys))
  ([xs ys opts] (ssj/b-spline xs ys opts)))

;; 2d grid

(defn bilinear
  "Creates a 2D bilinear grid interpolator.

  Linearly interpolates within each rectangular cell of a regular or irregular
  grid. Given a query point `(x, y)`, it locates the enclosing grid cell and
  blends the four corner values using bilinear weights.

  Parameters:

  - `xs` - a sequence of x grid coordinates; must be sorted in ascending order.
  - `ys` - a sequence of y grid coordinates; must be sorted in ascending order.
  - `vss` - a 2D sequence of values; `(get-in vss [i j])` is the value at `(xs[i], ys[j])`.

  Returns a function `f` accepting either `(f x y)` or `(f [x y])` that returns
  the interpolated `double` value at `(x, y)`.

  See also [[bicubic]], [[cubic-2d]], [[interpolation]]."
  [xs ys vss]
  (linear/bilinear xs ys vss))

(defn bicubic
  "Creates a 2D bicubic grid interpolator using the ACM library.

  Fits a bicubic polynomial within each rectangular cell of the grid, providing
  smoother results than bilinear interpolation. Requires at least 3 grid points
  along each axis.

  Parameters:

  - `xs` - a sequence of x grid coordinates; must be sorted in ascending order; at least 3 points.
  - `ys` - a sequence of y grid coordinates; must be sorted in ascending order; at least 3 points.
  - `vss` - a 2D sequence of values; `(get-in vss [i j])` is the value at `(xs[i], ys[j])`.

  Returns a function `f` accepting either `(f x y)` or `(f [x y])` that returns
  the interpolated `double` value at `(x, y)`.

  See also [[bilinear]], [[cubic-2d]], [[interpolation]]."
  [xs ys vss]
  (acm/bicubic xs ys vss))

(defn cubic-2d
  "Creates a 2D cubic spline grid interpolator.

  Applies 1D natural cubic splines first along each row (the `ys` axis) and then
  along the `xs` axis, giving a smooth surface that passes through all grid points.

  Parameters:

  - `xs` - a sequence of x grid coordinates; must be sorted in ascending order.
  - `ys` - a sequence of y grid coordinates; must be sorted in ascending order.
  - `vss` - a 2D sequence of values; `(get-in vss [i j])` is the value at `(xs[i], ys[j])`.

  Returns a function `f` accepting either `(f x y)` or `(f [x y])` that returns
  the interpolated `double` value at `(x, y)`.

  See also [[bilinear]], [[bicubic]], [[interpolation]]."
  [xs ys vss]
  (cubic/cubic-2d xs ys vss))

;; multidim

(defn microsphere-projection
  "Creates a multivariate interpolator using the ACM microsphere projection method.

  Assigns each query point a synthetic microsphere and illuminates it from
  the direction of each training point, weighted by the inverse distance to
  the power of `:exponent`. The predicted value is the weighted average of
  the illuminated hemisphere. Works for any number of input dimensions; the
  dimensionality is inferred from the first element of `xss`.

  Parameters:

  - `xss` - a sequence of input points; each point is either a number (1D) or a sequence of coordinates.
  - `ys` - a sequence of output values corresponding to each point in `xss`.
  - `opts` - (optional) a map of options:
    - `:elements` - number of surface elements on the microsphere; defaults to `max(5, 4+dims, n/5)`.
    - `:exponent` - distance weighting exponent (default: `1.0`).
    - `:max-dark-friction` - maximum fraction of sphere that can remain dark (default: `0.9`).
    - `:dark-threshold` - minimum illumination level to consider a element lit (default: `0.01`).
    - `:background` - background illumination value (default: `0.0`).
    - `:no-interpolation-tolerance` - distance below which an exact match is returned (default: `1.0e-6`).

  Returns a function `f` that accepts a point in the same form as the training inputs and returns
  the interpolated `double` value. For 2D inputs `f` also accepts `(f x y)`.

  See also [[shepard]], [[rbf]], [[interpolation]]."
  ([xss ys] (acm/microsphere-projection xss ys))
  ([xss ys opts] (acm/microsphere-projection xss ys opts)))

(defn shepard
  "Creates a multivariate Shepard (inverse distance weighting) interpolator.

  Predicts the value at a query point as the weighted mean of all training values,
  where weights are the inverse `p`-th powers of distances. Exact data points are
  cached and returned directly without division by zero.

  Parameters:

  - `xss` - a sequence of input points; each point is either a number (1D) or a sequence of coordinates.
  - `ys` - a sequence of output values corresponding to each point in `xss`.
  - `opts` - (optional) a map of options:
    - `:p` - the distance power exponent; higher values give more local influence to nearby points (default: `2.0`).
    - `:distance` - a two-argument distance function; defaults to `euclidean-1d` for scalar inputs and `euclidean` for vector inputs.

  Returns a function `f` that accepts a point in the same form as the training inputs
  and returns the interpolated `double` value.

  See also [[microsphere-projection]], [[rbf]], [[interpolation]]."
  ([xss ys] (shepard/shepard xss ys))
  ([xss ys opts] (shepard/shepard xss ys opts)))

(defn rbf
  "Creates a multivariate Radial Basis Function (RBF) interpolator.

  Solves a linear system whose kernel matrix is built from a radial basis function
  evaluated on pairwise distances between training points. Optional Tikhonov
  regularisation (`:lambda`) stabilises the solution when the kernel matrix is
  ill-conditioned. A polynomial drift term can be added via `:polynomial-terms`.

  Parameters:

  - `xss` - a sequence of input points; each point is either a number (1D) or a sequence of coordinates.
  - `ys` - a sequence of output values corresponding to each point in `xss`.
  - `kernel` - (optional) a radial kernel function `r -> double`; defaults to a Gaussian RBF.
  - `opts` - (optional) a map of options:
    - `:kscale` - scalar multiplier applied to all kernel evaluations (default: `1.0`).
    - `:lambda` - Tikhonov regularisation parameter; `0.0` gives exact interpolation (default: `0.0`).
    - `:distance` - a two-argument distance function; defaults to `euclidean-1d` for scalar inputs and `euclidean` for vector inputs.
    - `:polynomial-terms` - a function `x -> seq` returning polynomial basis values at `x`; augments the kernel system with a polynomial drift.

  RBF kernels can be created with `fastmath.kernel/rbf` function.  

  Returns a function `f` that accepts a point in the same form as the training inputs
  and returns the interpolated `double` value.

  See also [[shepard]], [[kriging]], [[gp]], [[interpolation]]."
  ([xss ys] (rbf/rbf xss ys))
  ([xss ys kernel] (rbf/rbf xss ys kernel))
  ([xss ys kernel opts] (rbf/rbf xss ys kernel opts)))

(defn kriging
  "Creates a multivariate Kriging (Gaussian process regression with variogram) interpolator.

  Solves the Kriging system using a variogram model that describes the spatial
  covariance structure of the data. If no variogram is supplied, one is fitted
  automatically to the empirical variogram using a Gaussian model. A polynomial
  drift term (`:polynomial-terms`) enables universal Kriging.

  Parameters:

  - `xss` - a sequence of input points; each point is either a number (1D) or a sequence of coordinates.
  - `ys` - a sequence of output values corresponding to each point in `xss`.
  - `variogram` - (optional) a fitted variogram function `distance -> double`; auto-fitted when `nil`.
  - `opts` - (optional) a map of options:
    - `:error` - measurement error (nugget); can be a scalar or a per-point sequence (default: `0.0`).
    - `:distance` - a two-argument distance function; defaults to `euclidean-1d` for scalar inputs and `euclidean` for vector inputs.
    - `:polynomial-terms` - a function `x -> seq` returning polynomial basis values for universal Kriging drift; defaults to `(constantly [1.0])`.

  Returns a function `f` that accepts a point in the same form as the training inputs
  and returns the interpolated `double` value.

  `fastmath.kernel.variogram` contains various semi-variograms and fitting methods.  

  See also [[rbf]], [[gp]], [[interpolation]]."
  ([xss ys] (kriging/kriging xss ys))
  ([xss ys variogram] (kriging/kriging xss ys variogram))
  ([xss ys variogram opts] (kriging/kriging xss ys variogram opts)))

(defn gp
  "Creates a multivariate Gaussian Process (GP) interpolator.

  Builds a Gaussian Process conditioned on the training data using a covariance
  kernel, then returns a prediction function that queries the posterior mean.
  For full probabilistic output (mean and standard deviation) use
  `fastmath.interpolation.gp/gaussian-process` and `fastmath.interpolation.gp/predict` directly.

  Parameters:

  - `xss` - a sequence of input points; each point is either a number (1D) or a sequence of coordinates.
  - `ys` - a sequence of output values corresponding to each point in `xss`.
  - `kernel` - (optional) a covariance kernel function `(x1, x2) -> double`; defaults to a Gaussian kernel with scale `1.0`.
  - `opts` - (optional) a map of options:
    - `:kscale` - global scale factor applied to all kernel evaluations (default: `1.0`).
    - `:noise` - diagonal noise added to the covariance matrix for numerical stability (default: `1.0e-8`).
    - `:normalize?` - when `true`, standardises `ys` to zero mean and unit variance before fitting (default: `false`).

  Returns a function `f` that accepts a point in the same form as the training inputs
  and returns the predicted `double` mean value.

  `fastmath.kernel/kernel` can be used to create predefined kernels.  

  See also [[kriging]], [[rbf]], [[interpolation]]."
  ([xss ys] (gp/gp xss ys))
  ([xss ys kernel] (gp/gp xss ys kernel))
  ([xss ys kernel opts] (gp/gp xss ys kernel opts)))

;;

(defn extrapolation
  "Wraps a 1D interpolator with a boundary policy that controls behaviour outside `[start, end]`.

  Inside the domain the wrapped function always delegates to `interpolator`. Outside it, the
  chosen `method` determines the returned value. When called with three arguments `method`
  defaults to `:skip`.

  Parameters:

  - `interpolator` - a 1D function `(fn ^double [^double x] ...)` produced by any interpolation constructor.
  - `start` - the lower boundary of the valid domain as a `double`.
  - `end` - the upper boundary of the valid domain as a `double`.
  - `method` - the extrapolation policy; one of:
    - `:skip` (default) - passes all `x` values directly to `interpolator`, disabling boundary enforcement.
    - `:constant` - returns `(interpolator start)` for `x < start` and `(interpolator end)` for `x > end`.
    - `:zero` - returns `0.0` for any `x` outside `[start, end]`.
    - `:error` - throws `IndexOutOfBoundsException` for any `x` outside `[start, end]`.
    - a `double` number - returns that constant value for any `x` outside `[start, end]`.
    - a two-element sequential `[left right]` - returns `left` for `x < start` and `right` for `x > end`.
    - a map `{:left left-val :right right-val}` - same as the two-element vector form.

  Returns a new function `f` with the same signature as `interpolator`,
  i.e. `(fn ^double [^double x] ...)`.

  See also [[linear]], [[interpolation]]."
  ([interpolator ^double start ^double end]
   (extrapolation interpolator :skip start end))
  ([interpolator method ^double start ^double end]
   (cond
     (sequential? method) (let [[^double left ^double right] method]
                            (fn ^double [^double x]
                              (cond (m/< x start) left
                                    (m/> x end) right
                                    :else (interpolator x))))
     (number? method) (let [v (double method)]
                        (fn ^double [^double x]
                          (if (m/between? start end x) (interpolator x) v)))
     (map? method) (let [{:keys [left right]} method]
                     (extrapolation interpolator [left right] start end))
     (keyword? method) (case method
                         :constant (let [sv (interpolator start)
                                         ev (interpolator end)]
                                     (extrapolation interpolator [sv ev] start end))
                         :zero (extrapolation interpolator 0.0 start end)
                         :error (fn ^double [^double x]
                                  (if (m/between? start end x)
                                    (interpolator x)
                                    (throw (IndexOutOfBoundsException.
                                            (str "x=" x " not in [" start ", " end "] range")))))
                         interpolator)
     :else interpolator)))

