(ns fastmath.interpolation.gp
  "Gaussian process regression (interpolation) of scattered data.

  A zero mean Gaussian process with a covariance kernel, `kscale` times a kernel of `fastmath.kernel`, and an independent noise (a variance added to the diagonal). Training data can be normalised to the zero mean and the unit variance. The model is fixed (hyperparameters are not optimised).

  Functions:

  - [[gaussian-process]] creates the model, [[gp]] returns the posterior mean as a function.
  - [[predict]] and [[predict-all]] give the posterior mean and optionally the standard deviation at one or many points.
  - [[L]] is the log marginal likelihood of the model.
  - [[prior-samples]] and [[posterior-samples]] draw random functions at the given points.

  See more [here](https://nextjournal.com/generateme/gaussian-processes#gp%2B)"
  (:require [fastmath.core :as m]
            [fastmath.kernel :as k]
            [fastmath.vector :as v]
            [fastmath.matrix :as mat]
            [fastmath.random :as r]
            [fastmath.interpolation.common :as ic])
  (:import [clojure.lang IFn]
           [fastmath.java Array]
           [org.apache.commons.math3.exception MaxCountExceededException]
           [org.apache.commons.math3.random RandomGenerator]
           [org.apache.commons.math3.stat StatUtils]
           [org.apache.commons.math3.linear CholeskyDecomposition EigenDecomposition MatrixUtils
            NonPositiveDefiniteMatrixException RealMatrix RealVector SingularValueDecomposition]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

;; native gaussian processes
;; Based on: https://www.cs.ubc.ca/~nando/540-2013/lectures/gp.py

(defn- ->point
  "A point as a vector: a number is a one element vector."
  [x]
  (cond (number? x) [x]
        (instance? clojure.lang.PersistentVector x) x
        ;; Vec2 etc. are vectors, but the kernels need the same type of both points
        :else (into [] x)))

(defn- dimension-error!
  [expected got]
  (ic/invalid-input! :gp :dimension-mismatch
                     (str "the points have to have the dimension " expected ", got " got)
                     {:expected expected :got got}))

(defn- ->points
  "Converts the points to vectors, checks that all have the dimension `dim` (any when `dim` is negative)."
  [xss ^long dim]
  (let [pts (mapv ->point xss)
        dim (long (if (m/neg? dim) (count (first pts)) dim))]
    (doseq [pt pts]
      (when-not (m/== dim (count pt))
        (dimension-error! dim (count pt))))
    pts))

(defn- kernel-matrix
  "Matrix of `scale * kernel(x, x')` for `x` from `xss1` (rows) and `x'` from `xss2` (columns)."
  ^RealMatrix [kernel xss1 xss2 ^double scale]
  (let [n1 (count xss1)
        n2 (count xss2)
        arr (make-array Double/TYPE n1 n2)]
    (dotimes [i n1]
      (let [^doubles row (aget ^"[[D" arr i)
            xi (xss1 i)]
        (dotimes [j n2]
          (Array/aset row j (m/* scale (double (kernel (xss2 j) xi)))))))
    (MatrixUtils/createRealMatrix ^"[[D" arr)))

(defn- add-diagonal
  "Adds `diag` to the diagonal of the (mutable) square matrix, returns the matrix."
  ^RealMatrix [^RealMatrix a ^double diag]
  (dotimes [i (.getRowDimension a)]
    (.addToEntry a i i diag))
  a)

(declare predict)
(declare predict-all)

(defrecord GaussianProcess [^double kscale kernel ^double noise xss ys ^double ymean ^double ystddev
                            chol w ^long n ^double L ^long dim]
  IFn
  (invoke [this xs]
    (if (and (sequential? xs)
             (sequential? (first xs)))
      (predict-all this xs)
      (predict this xs)))
  (invoke [this xs stddev?]
    (if (and (sequential? xs)
             (sequential? (first xs)))
      (predict-all this xs stddev?)
      (predict this xs stddev?))))

(def ^{:const true :private true :tag 'double} LOG2PI (m/log m/TWO_PI))

(defn- log-det
  "Logarithm of the determinant of the matrix with the Cholesky factor `L`: twice the sum of the logarithms of its diagonal."
  ^double [^CholeskyDecomposition chol]
  (let [^RealMatrix l (.getL chol)
        n (.getRowDimension l)]
    (m/* 2.0 (double (loop [i (long 0)
                            sum 0.0]
                       (if (m/< i n)
                         (recur (m/inc i) (m/+ sum (m/log (.getEntry l i i))))
                         sum))))))

(defn L
  "Computes the log marginal likelihood of a Gaussian process.

  `L = -1/2 (y^T K^-1 y + log det K + n log 2 pi)`, where `K` is the covariance matrix of the training points (the kernel times `kscale` plus the noise on the diagonal) and `y` are the training values (normalised when the model was created with `:normalize? true`, so `L` is the likelihood of the normalised values). The determinant is computed stably as twice the sum of the logarithms of the diagonal of the Cholesky factor: it does not underflow or overflow for hundreds of points or an extreme `kscale`.

  Parameters:

  - `gp-object` (`GaussianProcess`): a model created by [[gaussian-process]]. The value computed at the creation is returned when it exists (the option `:L?` was true), otherwise it is computed
  - `ys`, `w`, `chol`, `n`: (long form) the training values (`RealVector`), the weights `K^-1 y` (`RealVector`), the Cholesky decomposition of `K` and the number of points

  Returns the log marginal likelihood as a double. It is `NaN` for a model whose `ys` are not finite.

  See also [[gaussian-process]], [[predict]]."
  (^double [ys w ^CholeskyDecomposition chol ^long n]
   (m/* -0.5 (m/+ (m/+ (v/dot ys w)
                       (log-det chol))
                  (m/* n LOG2PI))))
  (^double [^GaussianProcess gp-object]
   (if (m/invalid-double? (.L gp-object))
     (L (.ys gp-object)
        (.w gp-object)
        (.chol gp-object)
        (.n gp-object))
     (.L gp-object))))

(def ^:private ^:const absolute-positivity-threshold 0.0)

(defn- cholesky
  "Cholesky decomposition of the covariance matrix with a readable error."
  ^CholeskyDecomposition [^RealMatrix cov]
  (try
    (CholeskyDecomposition. cov CholeskyDecomposition/DEFAULT_RELATIVE_SYMMETRY_THRESHOLD
                            absolute-positivity-threshold)
    (catch NonPositiveDefiniteMatrixException e
      (throw (ex-info (str "gp: the covariance matrix of the training points is not positive definite "
                           "(numerically), increase :noise or check for repeated points")
                      {:type :invalid-input :method :gp :property :not-positive-definite}
                      e)))))

(defn- check-option!
  [option value ok? expected]
  (when-not (and (number? value) (ok? (double value)))
    (ic/invalid-input! :gp option
                       (str option " has to be " expected ", got " (pr-str value))
                       {:option option :value value :expected expected})))

(defn gaussian-process
  "Creates a Gaussian process model from the training data.

  The posterior is the one of a zero mean Gaussian process with the covariance `kscale * kernel(x, x')`, observed with an independent noise of the variance `noise`.

  Parameters:

  - `xss` (sequence): the training points, at least one; numbers for one dimensional data or sequences of numbers (vectors, `Vec2`, ...) of an equal dimension
  - `ys` (sequence of numbers): the values at the points, as many as points
  - `opts` (map, optional; a `nil` value is the same as a missing key):
      - `:kernel` (function, default: `(fastmath.kernel/kernel :gaussian)`, sigma `1.0`): a covariance kernel `(fn [x1 x2])` of two points given as vectors, see [[fastmath.kernel/kernel]]; its own parameters (for example `:sigma`) set the length scale
      - `:kscale` (double, default: `1.0`): a positive multiplier of the kernel, the signal variance
      - `:noise` (double, default: `1.0e-8`): a variance added to the diagonal of the covariance matrix of the training values, in the units of the normalised values (so relative to the signal variance `kscale` when `:normalize?` is true), finite and not negative. It is the `alpha` of the scikit-learn `GaussianProcessRegressor`. A small value (`1.0e-8`) keeps the matrix numerically positive definite; a larger one smooths the data (the model then does not pass through the points)
      - `:normalize?` (boolean, default: `false`): shift the values to the zero mean and scale them to the unit sample standard deviation (with `n - 1`, scikit-learn `normalize_y` uses `n`) before the fit, predictions are mapped back. A constant `ys` (or a single point) uses the standard deviation `1`: the predictions equal the constant, the standard deviation comes from the kernel only
      - `:L?` (boolean, default: `true`): compute the log marginal likelihood [[L]] now and store it in the model; `false` skips it (`NaN` is stored, [[L]] computes it on demand). [[gp]] uses `false`

  Returns a `GaussianProcess` record, which is also a function: `(gp-object x)` or `(gp-object x true)` calls [[predict]] for one point and [[predict-all]] for a sequence of points.

  Throws `ex-info` (with `:type :invalid-input`, `:method :gp` and the key `:property`) when `xss` and `ys` have different lengths (`:length-mismatch`), there are no points (`:too-few-points`), the training points have different dimensions (`:dimension-mismatch`), `:kscale` or `:noise` is out of range (`:kscale`, `:noise`), or the covariance matrix is not positive definite (`:not-positive-definite`: repeated points with a tiny `:noise`; increase it).

  See also [[gp]], [[predict]], [[predict-all]], [[L]], [[posterior-samples]], [[fastmath.interpolation/gp]]."
  ([xss ys] (gaussian-process xss ys nil))
  ([xss ys opts]
   (let [{:keys [kscale kernel noise normalize? L?]
          :or {kscale 1.0 normalize? false noise 1.0e-8 L? true}}
         (into {} (remove (fn [[_ v]] (nil? v))) opts)
         _ (check-option! :kscale kscale (fn [^double x] (and (m/valid-double? x) (m/pos? x))) "a finite number > 0")
         _ (check-option! :noise noise (fn [^double x] (and (m/valid-double? x) (m/>= x 0.0))) "a finite number >= 0")
         kscale (double kscale)
         noise (double noise)
         kernel (or kernel (k/kernel :gaussian))
         ys (if (nil? ys) (double-array 0) (m/seq->double-array ys))
         nx (count xss)]
     (when-not (m/== nx (alength ys))
       (ic/invalid-input! :gp :length-mismatch
                          (str "xss and ys have different lengths: " nx " and " (alength ys))
                          {:xs-count nx :ys-count (alength ys)}))
     (when (m/zero? nx)
       (ic/invalid-input! :gp :too-few-points "at least 1 point is required, got 0"
                          {:axis :xss :min-points 1 :count 0}))
     (let [xss (->points xss -1)
           dim (count (first xss))
           normalize? (if (m/== 1 (alength ys)) false normalize?)
           ymean (if normalize? (StatUtils/mean ys) 0.0)
           ystddev (if normalize?
                     (let [sd (m/sqrt (StatUtils/variance ys))]
                       ;; a constant ys has no scale
                       (if (or (m/zero? sd) (m/== (StatUtils/min ys) (StatUtils/max ys))) 1.0 sd))
                     1.0)
           ys (v/vec->RealVector (if normalize?
                                   (m/seq->double-array (map (fn [^double y] (m// (m/- y ymean) ystddev)) ys))
                                   ys))
           chol (cholesky (add-diagonal (kernel-matrix kernel xss xss kscale) noise))
           w (.solve (.getSolver chol) ^RealVector ys)
           n (count xss)
           l (if-not L? ##NaN (L ys w chol n))]
       (->GaussianProcess kscale kernel noise xss ys ymean ystddev chol w n l dim)))))

(defn- cross-cov
  "Matrix `m x n` of the covariances between `m` query points and the `n` training points."
  ^RealMatrix [^GaussianProcess gp-object xvals]
  (kernel-matrix (.kernel gp-object) xvals (.xss gp-object) (.kscale gp-object)))

(defn- query-vector
  "Vector of the covariances between the query point `xtest` and the training points."
  ^RealVector [^GaussianProcess gp-object xtest]
  (let [xss (.xss gp-object)
        kernel (.kernel gp-object)
        kscale (.kscale gp-object)
        n (.n gp-object)
        arr (double-array n)]
    (dotimes [i n]
      (Array/aset arr i (m/* kscale (double (kernel xtest (xss i))))))
    (MatrixUtils/createRealVector arr)))

(defn predict
  "Computes the posterior mean (and optionally the standard deviation) of a Gaussian process at one point.

  Parameters:

  - `gp-object` (`GaussianProcess`): a model created by [[gaussian-process]]
  - `xval` (number or sequence of numbers): the query point, a number for one dimensional data; its dimension has to be the dimension of the training points
  - `stddev?` (boolean, default: `false`): also return the standard deviation

  Returns the mean as a double, or, when `stddev?` is true, a `Vec2` of the mean and the standard deviation of the latent function (it does not include `:noise`; it is `0.0` up to the rounding at a training point for a small noise, the prior `sqrt(kscale)` times the standard deviation of `ys` far from the data). Both are in the original units of `ys` (also with `:normalize? true`). A tiny negative variance, a rounding artefact, gives the standard deviation `0.0`.

  Throws `ex-info` with `:property :dimension-mismatch` when the dimension of `xval` differs from the training dimension.

  See also [[predict-all]], [[gaussian-process]], [[posterior-samples]]."
  ([gp-object xval] (predict gp-object xval false))
  ([^GaussianProcess gp-object xval stddev?]
   (let [xtest (->point xval)
         _ (when-not (m/== (.dim gp-object) (count xtest))
             (dimension-error! (.dim gp-object) (count xtest)))
         cov-vector (query-vector gp-object xtest)
         mu (m/+ (.ymean gp-object)
                 (m/* (.ystddev gp-object)
                      (.dotProduct ^RealVector (.w gp-object) cov-vector)))]
     (if-not stddev?
       mu
       (v/vec2 mu (m/* (.ystddev gp-object)
                       (m/safe-sqrt
                        (m/- (m/* (.kscale gp-object) (double ((.kernel gp-object) xtest xtest)))
                             (.dotProduct cov-vector
                                          (.solve (.getSolver ^CholeskyDecomposition (.chol gp-object))
                                                  cov-vector))))))))))

(defn gp
  "Creates the Gaussian process interpolator, a function returning the posterior mean.

  Parameters:

  - `xss`, `ys`: the training points and values, see [[gaussian-process]]
  - `kernel` (function, optional, default: `(fastmath.kernel/kernel :gaussian)`): the covariance kernel
  - `opts` (map, optional): the options of [[gaussian-process]] (`:kscale`, `:noise`, `:normalize?`); `:kernel` is taken from the argument and `:L?` is always `false`

  Returns a function of a query point (a number for one dimensional data, a sequence of numbers otherwise) returning the posterior mean as a double. Throws the `ex-info` of [[gaussian-process]] and of [[predict]].

  See also [[gaussian-process]], [[predict]], [[fastmath.interpolation/gp]]."
  ([xss ys] (gp xss ys nil))
  ([xss ys kernel] (gp xss ys kernel nil))
  ([xss ys kernel opts]
   (let [gp-object (gaussian-process xss ys (assoc opts :kernel kernel :L? false))]
     (fn ^double [xs] (predict gp-object xs false)))))

(defn- post-variances
  "Posterior variances (normalised units) at the query points, `kv` is the cross covariance matrix."
  ^doubles [^GaussianProcess gp-object xvals ^RealMatrix kv]
  (let [m (count xvals)
        ^RealMatrix solved (.solve (.getSolver ^CholeskyDecomposition (.chol gp-object)) (.transpose kv))
        kernel (.kernel gp-object)
        kscale (.kscale gp-object)
        arr (double-array m)]
    (dotimes [i m]
      (Array/aset arr i (m/- (m/* kscale (double (kernel (xvals i) (xvals i))))
                             (.dotProduct (.getRowVector kv i) (.getColumnVector solved i)))))
    arr))

(defn- mean-values
  "Posterior means (original units) at the query points."
  ^doubles [^GaussianProcess gp-object ^RealMatrix kv]
  (let [^RealVector mus (.operate kv ^RealVector (.w gp-object))
        arr (.toArray mus)
        ystddev (.ystddev gp-object)
        ymean (.ymean gp-object)]
    (dotimes [i (alength arr)]
      (Array/aset arr i (m/+ ymean (m/* ystddev (Array/aget arr i)))))
    arr))

(defn- with-sd
  "Pairs `[value sd]` for the values and the posterior variances in the normalised units."
  [^double ystddev values ^doubles variances]
  (map (fn [mu ^double variance] [mu (m/* ystddev (m/safe-sqrt variance))])
       values variances))

(defn predict-all
  "Computes the posterior mean (and optionally the standard deviation) of a Gaussian process at many points.

  Parameters:

  - `gp-object` (`GaussianProcess`): a model created by [[gaussian-process]]
  - `xvals` (sequence): the query points (numbers for one dimensional data, sequences of numbers otherwise); their dimension has to be the dimension of the training points
  - `stddev?` (boolean, default: `false`): also return the standard deviations

  Returns a sequence of the means (doubles), or, when `stddev?` is true, a sequence of vectors `[mean standard-deviation]` in the original units of `ys`, the same values as [[predict]] gives for each point (a tiny negative variance, a rounding artefact, gives the standard deviation `0.0`, never `NaN`).

  Throws `ex-info` with `:property :dimension-mismatch` when a point has a different dimension than the training points.

  See also [[predict]], [[gaussian-process]], [[posterior-samples]]."
  ([gp-object xvals] (predict-all gp-object xvals false))
  ([^GaussianProcess gp-object xvals stddev?]
   (let [xvals (->points xvals (.dim gp-object))]
     (if (empty? xvals)
       ()
       (let [kv (cross-cov gp-object xvals)
             mus (mean-values gp-object kv)]
         (if-not stddev?
           (seq mus)
           (with-sd (.ystddev gp-object) mus (post-variances gp-object xvals kv))))))))

(defn- factor-eigen
  "Returns the matrix `B` with `B B^T` equal to the covariance, from the eigenvalues (negative ones are clipped to zero)."
  ^RealMatrix [^RealMatrix cov]
  (let [eig (EigenDecomposition. cov)
        lambdas (.getRealEigenvalues eig)
        ^RealMatrix vs (.getV eig)
        n (alength lambdas)
        b (.copy vs)]
    (dotimes [j n]
      (let [s (m/sqrt (m/max 0.0 (Array/aget lambdas j)))]
        (dotimes [i n]
          (.setEntry b i j (m/* s (.getEntry vs i j))))))
    b))

(defn- factor-svd
  "Returns the matrix `B` with `B B^T` equal to the covariance, from the singular value decomposition, for the matrices on which the eigenvalue iteration does not converge."
  ^RealMatrix [^RealMatrix cov]
  (let [svd (SingularValueDecomposition. cov)
        sigmas (.getSingularValues svd)
        ^RealMatrix u (.getU svd)
        n (alength sigmas)
        b (.copy u)]
    (dotimes [j n]
      (let [s (m/sqrt (Array/aget sigmas j))]
        (dotimes [i n]
          (.setEntry b i j (m/* s (.getEntry u i j))))))
    b))

(defn- sample-normal
  "Draws a sample from the normal distribution of the mean `mu` and the covariance `cov` (positive semidefinite up to rounding).

  The covariance is symmetrised and decomposed into the eigenvalues, negative eigenvalues (rounding errors of a rank deficient covariance) are clipped to zero, so a dense grid or a point at a training point is fine. The singular value decomposition is the fallback when the eigenvalue iteration of Apache Commons Math does not converge."
  ^doubles [^doubles mu ^RealMatrix cov ^RandomGenerator rng]
  (let [n (alength mu)
        sym (mat/symmetrize cov)
        ^RealMatrix b (try (factor-eigen sym)
                           (catch MaxCountExceededException _ (factor-svd sym)))
        z (double-array n)]
    (dotimes [i n]
      (Array/aset z i (r/grandom rng)))
    (let [res (.toArray (.operate b (MatrixUtils/createRealVector z)))]
      (dotimes [i n]
        (Array/aset res i (m/+ (Array/aget res i) (Array/aget mu i))))
      res)))

(defn- scale-matrix
  ^RealMatrix [^RealMatrix a ^double s]
  (.scalarMultiply a s))

(defn prior-samples
  "Draws a random sample of the prior of a Gaussian process at the given points.

  The sample is drawn from the normal distribution with the mean `ymean` and the covariance `ystddev^2 * (kscale * kernel + noise I)` at `xvals` (for a model created with `:normalize? false` this is the zero mean and the covariance `kscale * kernel + noise I`). The noise is part of the covariance, so the sample is a noisy observation of the prior. The training data do not take part (only the kernel, `kscale`, `noise` and the normalisation).

  Parameters:

  - `gp-object` (`GaussianProcess`): a model created by [[gaussian-process]]
  - `xvals` (sequence): the points (numbers for one dimensional data, sequences of numbers otherwise); their dimension has to be the dimension of the training points
  - `rng` (random number generator, optional, default: `nil`): the generator of the random numbers, see [[fastmath.random/ensure-rng]]; when missing or `nil` a new unseeded generator is created for every call (the shared `fastmath.random/default-rng` is never used), so give a seeded one to reproduce the samples, for example `(fastmath.random/rng :isaac 42)`

  Returns a sequence of values, one for each point. The covariance is decomposed into the eigenvalues with the negative ones clipped to zero, so a dense grid does not fail.

  Throws `ex-info` with `:property :dimension-mismatch` for a point of a wrong dimension.

  See also [[posterior-samples]], [[gaussian-process]]."
  ([gp-object xvals] (prior-samples gp-object xvals nil))
  ([^GaussianProcess gp-object xvals rng]
   (let [xvals (->points xvals (.dim gp-object))
         rng (r/ensure-rng rng)
         ystddev (.ystddev gp-object)]
     (if (empty? xvals)
       ()
       (let [cov (scale-matrix (add-diagonal (kernel-matrix (.kernel gp-object) xvals xvals (.kscale gp-object))
                                             (.noise gp-object))
                               (m/* ystddev ystddev))
             mu (double-array (count xvals) (.ymean gp-object))]
         (seq (sample-normal mu cov rng)))))))

(defn- posterior-cov
  "Posterior covariance matrix of the latent function at the query points (already converted), in the original units of `ys`."
  ^RealMatrix [^GaussianProcess gp-object xvals ^RealMatrix kv]
  (let [^RealMatrix solved (.solve (.getSolver ^CholeskyDecomposition (.chol gp-object)) (.transpose kv))
        ystddev (.ystddev gp-object)
        ^RealMatrix cov (.subtract ^RealMatrix (kernel-matrix (.kernel gp-object) xvals xvals (.kscale gp-object))
                                   (.multiply kv solved))]
    (scale-matrix cov (m/* ystddev ystddev))))

(defn posterior-samples
  "Draws a random sample of the posterior of a Gaussian process at the given points.

  The sample is drawn from the joint normal distribution of the posterior at `xvals`: the mean of [[predict-all]] and the full posterior covariance of the latent function (without `:noise`), in the original units of `ys` (the covariance is multiplied by `ystddev^2`, so the standard deviation of the samples at a point equals the one of [[predict]], also with `:normalize? true`).

  Parameters:

  - `gp-object` (`GaussianProcess`): a model created by [[gaussian-process]]
  - `xvals` (sequence): the points (numbers for one dimensional data, sequences of numbers otherwise); their dimension has to be the dimension of the training points
  - `stddev?` (boolean, default: `false`): return pairs `[value standard-deviation]` where the second element is the posterior standard deviation at the point
  - `rng` (random number generator, optional, default: `nil`): the generator of the random numbers, see [[fastmath.random/ensure-rng]]; when missing or `nil` a new unseeded generator is created for every call (the shared `fastmath.random/default-rng` is never used), so give a seeded one to reproduce the samples, for example `(fastmath.random/rng :isaac 42)`; to give it, give `stddev?` as well

  Returns a sequence of values, one for each point (or pairs). The covariance is symmetrised and decomposed into the eigenvalues, negative eigenvalues (rounding) are clipped to zero: a dense grid (hundreds of points) or points at the training points do not fail. The sample is one random function evaluated at the points, the values at different points are correlated.

  Throws `ex-info` with `:property :dimension-mismatch` for a point of a wrong dimension.

  See also [[prior-samples]], [[predict-all]], [[gaussian-process]]."
  ([gp-object xvals] (posterior-samples gp-object xvals false nil))
  ([gp-object xvals stddev?] (posterior-samples gp-object xvals stddev? nil))
  ([^GaussianProcess gp-object xvals stddev? rng]
   (let [xvals (->points xvals (.dim gp-object))
         rng (r/ensure-rng rng)]
     (if (empty? xvals)
       ()
       (let [kv (cross-cov gp-object xvals)
             mus (mean-values gp-object kv)
             cov (posterior-cov gp-object xvals kv)
             post (sample-normal mus cov rng)]
         (if-not stddev?
           (seq post)
           (map (fn [value i] [value (m/safe-sqrt (.getEntry cov (long i) (long i)))])
                post (range (alength post)))))))))
