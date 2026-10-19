;; http://shihchinw.github.io/2018/10/data-interpolation-with-radial-basis-functions-rbfs.html
(ns fastmath.interpolation.rbf
  "Radial basis function interpolation of scattered data."
  (:require [fastmath.vector :as v]
            [fastmath.kernel :as k]
            [fastmath.core :as m]
            [fastmath.matrix :as mat]
            [fastmath.interpolation.common :as ic])
  (:import [fastmath.java Array]
           [org.apache.commons.math3.linear MatrixUtils
            RealMatrix RealVector]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- ->points
  "Converts the nodes to double arrays (a number is a one element point), checks their dimension."
  [xss]
  (let [pts (mapv m/seq->double-array xss)
        dim (alength ^doubles (first pts))]
    (doseq [^doubles pt pts]
      (when-not (m/== dim (alength pt))
        (ic/invalid-input! :rbf :dimension-mismatch
                           (str "all the nodes have to have the same dimension, got " dim " and " (alength pt))
                           {:expected dim :got (alength pt)})))
    pts))

(defn- query-converter
  "Returns a function converting a query point to a double array of the dimension `dim`."
  [^long dim]
  (fn ^doubles [x]
    (let [q (m/seq->double-array x)]
      (when-not (m/== dim (alength q))
        (ic/invalid-input! :rbf :dimension-mismatch
                           (str "the query has the dimension " (alength q) ", the nodes have " dim)
                           {:expected dim :got (alength q)}))
      q)))

(def ^:private regularizations #{:normal :ridge})

(defn- check-lambda!
  [lambda]
  (when-not (and (number? lambda) (m/valid-double? (double lambda)) (m/>= (double lambda) 0.0))
    (ic/invalid-input! :rbf :lambda
                       (str ":lambda has to be a finite number >= 0, got " (pr-str lambda))
                       {:option :lambda :value lambda :expected "a finite number >= 0"})))

(defn- check-regularization!
  [regularization]
  (when-not (regularizations regularization)
    (ic/invalid-input! :rbf :regularization
                       (str ":regularization has to be one of " (vec (sort regularizations))
                            ", got " (pr-str regularization))
                       {:option :regularization :value regularization :expected (vec (sort regularizations))})))

(defn- check-data!
  [xss ys]
  (let [nx (count xss)
        ny (count ys)]
    (when-not (m/== nx ny)
      (ic/invalid-input! :rbf :length-mismatch
                         (str "xss and ys have different lengths: " nx " and " ny)
                         {:xs-count nx :ys-count ny}))
    (when (m/zero? nx)
      (ic/invalid-input! :rbf :too-few-points "at least 1 point is required, got 0"
                         {:axis :xss :min-points 1 :count 0}))))

(defn- phi-matrix
  "Kernel matrix of the function `kfn` of two nodes, `diag` is added to the diagonal."
  ^RealMatrix [xss kfn ^double diag]
  (let [n (count xss)
        arr (make-array Double/TYPE n n)]
    (dotimes [i n]
      (let [^doubles row (aget ^"[[D" arr i)
            xi (xss i)]
        (dotimes [j n]
          (Array/aset row j (double (kfn xi (xss j)))))
        (Array/aset row i (m/+ (Array/aget row i) diag))))
    (MatrixUtils/createRealMatrix ^"[[D" arr)))

(defn- solve-normal
  "Solves the regularized normal equations `(A^T A + lambda I) x = A^T b`."
  [^RealMatrix A ^RealVector b ^double lambda]
  (let [At (mat/transpose A)]
    (ic/solve-LU-or-SVD (mat/add (mat/mulm At A)
                                 (MatrixUtils/createRealDiagonalMatrix
                                  (m/seq->double-array (repeat (.getDimension b) lambda))))
                        (mat/mulv At b))))

(defn rbf
  "Creates the radial basis function (RBF) interpolator of scattered data.

  The interpolant is `f(x) = sum_i w_i kscale kernel(d(x, xs_i))` (plus a polynomial `sum_j c_j p_j(x)` when `:polynomial-terms` are given), where the weights `w` (and the coefficients `c`) are found from the linear system built on the nodes. Without regularization the interpolant passes through all the nodes.

  Parameters:

  - `xss` (sequence): the nodes, at least one; numbers for one dimensional data or sequences of numbers
  - `ys` (sequence of numbers): the values at the nodes, as many as nodes
  - `kernel` (function, optional, default: `(fastmath.kernel/rbf :gaussian)`): a radial function of the distance (see [[fastmath.kernel/rbf]]); its own `:scale` sets the length scale of the interpolation
  - `opts` (map, optional; a `nil` value is the same as a missing key):
      - `:kscale` (double, default: `1.0`): a multiplier applied to all the kernel evaluations
      - `:lambda` (double, default: `0.0`): regularization strength, finite and not negative; `0.0` is the exact interpolation
      - `:regularization` (keyword, default: `:normal`): how `:lambda` enters the system
          - `:normal` - ridge on the normal equations: `(A^T A + lambda I) w = A^T y`, where `A` is the whole system matrix (the polynomial block included, which is regularized too); it squares the condition number. This is the original behaviour
          - `:ridge` - kernel ridge regression: `(K + lambda I) w = y`, `lambda` is added to the kernel matrix only (the polynomial block is not regularized). This is the `smoothing` of the scipy `RBFInterpolator`
      - `:distance` (function, default: Euclidean distance, `fastmath.distance/euclidean-1d` for numbers and `fastmath.distance/euclidean` for sequences): a binary distance function `(fn [x1 x2])`
      - `:polynomial-terms` (function, default: `nil`, no polynomial): a function of a point returning a sequence of the values of the polynomial basis (for example `(fn [[x y]] [1.0 x y])` for a linear polynomial in two dimensions); the system is augmented with the polynomial constraints

  Returns a function of a query point (a number for one dimensional data, a sequence of numbers otherwise) returning a double.

  The relation of `:kscale` and `:lambda`: without regularization `:kscale` cancels (the weights scale by its inverse) and has no effect on the values. With `:regularization :normal` and no polynomial terms `:lambda` acts as `lambda / kscale^2`, with `:ridge` (also with polynomial terms) as `lambda / kscale`, and with `:normal` and polynomial terms the effect is mixed (the polynomial block is not scaled by `:kscale`). Tune `:lambda` for a fixed `:kscale` (the default `1.0` is the simplest).

  Conditionally positive definite kernels (thin plate splines, cubic and other radial powers, linear, generalized multiquadratic of a positive order, shifted surface splines) need `:polynomial-terms` of the matching degree (at least a constant for the linear kernel and the multiquadratic, a linear polynomial for the thin plate splines and the cubic kernel) to have a guaranteed unique solution and to reproduce the polynomials. Without them the system is often still solvable but the result differs from the correct one. Positive definite kernels (gaussian, inverse multiquadratic, matern, compactly supported kernels) need no polynomial.

  The linear system is solved by the LU decomposition. When the matrix is (numerically) singular, for example for repeated nodes, the minimum norm least squares solution of the SVD is used silently, so repeated nodes with different values give the mean of the values. The system of a smooth kernel (gaussian, multiquadratic) is ill-conditioned for dense nodes and a large length scale: the result is accurate only up to the conditioning (often `1e-8` or worse). A small `:lambda` with `:ridge` stabilises it.

  Throws `ex-info` (with `:type :invalid-input`, `:method :rbf` and the key `:property`) when `xss` and `ys` have different lengths (`:length-mismatch`), there are no nodes (`:too-few-points`), `:lambda` is negative, infinite or `NaN` (`:lambda`) or `:regularization` is not one of `:normal`, `:ridge` (`:regularization`).

  See also [[fastmath.interpolation/rbf]], [[fastmath.kernel/rbf]], [[fastmath.interpolation.kriging/kriging]], [[fastmath.interpolation.shepard/shepard]]."
  ([xss ys] (rbf xss ys nil))
  ([xss ys kernel] (rbf xss ys kernel nil))
  ([xss ys kernel opts]
   (let [{:keys [^double kscale ^double lambda regularization distance polynomial-terms]
          :or {kscale 1.0 lambda 0.0 regularization :normal}}
         (into {} (remove (fn [[_ v]] (nil? v))) opts)
         _ (check-lambda! lambda)
         _ (check-regularization! regularization)
         _ (check-data! xss ys)
         kernel (or kernel (k/rbf :gaussian))
         xss (vec xss)
         ysize (count xss)
         ;; the default distance works on arrays: numbers, vectors, Vec2 and lists can be mixed
         pts (if distance xss (->points xss))
         dist (or distance v/dist)
         ->query (if distance identity (query-converter (alength ^doubles (first pts))))
         ridge? (and (= :ridge regularization) (m/pos? lambda))
         kfn (fn ^double [xs1 xs2] (m/* kscale (double (kernel (dist xs1 xs2)))))
         Phi (phi-matrix pts kfn (if ridge? lambda 0.0))
         [Phi ys] (if polynomial-terms (ic/pi-matrix-ys xss Phi ys polynomial-terms) [Phi ys])
         ysv (v/vec->RealVector ys)
         w (if (or ridge? (m/zero? lambda))
             (ic/solve-LU-or-SVD Phi ysv)
             (solve-normal Phi ysv lambda))
         ^doubles ws (m/seq->double-array (take ysize w))
         c (when polynomial-terms (vec (drop ysize w)))]
     (fn ^double [x]
       (let [q (->query x)
             s (double
                (loop [i (long 0)
                       sum 0.0]
                  (if (m/< i ysize)
                    (recur (m/inc i)
                           (m/+ sum (m/* (Array/aget ws i)
                                         (double (kfn (pts i) q)))))
                    sum)))]
         (if c
           (m/+ s (v/dot c (polynomial-terms x)))
           s))))))
