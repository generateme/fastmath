;; Numerical Recipes
(ns fastmath.interpolation.kriging
  "Kriging interpolation of scattered data."
  (:require [fastmath.vector :as v]
            [fastmath.distance :as d]
            [fastmath.core :as m]
            [fastmath.kernel.variogram :as vgram]
            [fastmath.interpolation.common :as ic])
  (:import [fastmath.java Array]
           [org.apache.commons.math3.linear MatrixUtils]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- ->points
  "Converts the nodes to double arrays (a number is a one element point), checks their dimension."
  [xss]
  (let [pts (mapv m/seq->double-array xss)
        dim (alength ^doubles (first pts))]
    (doseq [^doubles pt pts]
      (when-not (m/== dim (alength pt))
        (ic/invalid-input! :kriging :dimension-mismatch
                           (str "all the nodes have to have the same dimension, got " dim " and " (alength pt))
                           {:expected dim :got (alength pt)})))
    pts))

(defn- query-converter
  "Returns a function converting a query point to a double array of the dimension `dim`."
  [^long dim]
  (fn ^doubles [x]
    (let [q (m/seq->double-array x)]
      (when-not (m/== dim (alength q))
        (ic/invalid-input! :kriging :dimension-mismatch
                           (str "the query has the dimension " (alength q) ", the nodes have " dim)
                           {:expected dim :got (alength q)}))
      q)))

(defn- check-data!
  [xss ys min-points]
  (let [nx (count xss)
        ny (count ys)]
    (when-not (m/== nx ny)
      (ic/invalid-input! :kriging :length-mismatch
                         (str "xss and ys have different lengths: " nx " and " ny)
                         {:xs-count nx :ys-count ny}))
    (when (m/< nx (long min-points))
      (ic/invalid-input! :kriging :too-few-points
                         (str "at least " min-points " points are required, got " nx)
                         {:axis :xss :min-points min-points :count nx}))))

(defn- error-variances
  "Returns the array of the squared errors of the `n` nodes or `nil` when there is no error."
  ^doubles [error ^long n]
  (cond
    (not error) nil
    (number? error) (let [e (double error)]
                      (when-not (m/valid-double? e)
                        (ic/invalid-input! :kriging :error (str ":error has to be finite, got " e)
                                           {:option :error :value error :expected "a finite number"}))
                      (when-not (m/zero? e)
                        (double-array n (m/sq e))))
    (or (sequential? error) (= (type error) m/double-array-type))
    (let [es (m/seq->double-array error)]
      (cond
        (m/zero? (alength es)) nil
        (not (m/== (alength es) n))
        (ic/invalid-input! :kriging :error
                           (str ":error as a sequence needs one value for each node: expected " n
                                ", got " (alength es))
                           {:option :error :value error :expected (str "a number or a sequence of " n " numbers")})
        (not (every? m/valid-double? es))
        (ic/invalid-input! :kriging :error ":error has to be finite"
                           {:option :error :value error :expected "finite numbers"})
        :else (let [arr (double-array n)]
                (dotimes [i n]
                  (Array/aset arr i (m/sq (Array/aget es i))))
                arr)))
    :else (ic/invalid-input! :kriging :error
                             (str ":error has to be a number or a sequence of numbers, got " (pr-str error))
                             {:option :error :value error :expected "a number or a sequence of numbers"})))

(defn- v-matrix
  "Variogram matrix of the nodes with the squared errors subtracted from the diagonal."
  [xss variogram distance error]
  (let [n (count xss)
        errs (error-variances error n)
        arr (make-array Double/TYPE n n)]
    (dotimes [i n]
      (let [^doubles row (aget ^"[[D" arr i)
            xi (xss i)]
        (dotimes [j n]
          (Array/aset row j (double (variogram (distance xi (xss j))))))
        (when errs
          (Array/aset row i (m/- (Array/aget row i) (Array/aget ^doubles errs i))))))
    (MatrixUtils/createRealMatrix ^"[[D" arr)))

(defn kriging
  "Creates the kriging interpolator of scattered data, ordinary or universal.

  The interpolant is `f(x) = sum_i w_i gamma(d(x, xs_i)) + sum_j c_j p_j(x)`, where `gamma` is the semivariogram and the weights `w` and the coefficients `c` solve the kriging system built on the nodes. The default polynomial is a constant (ordinary kriging).

  Parameters:

  - `xss` (sequence): the nodes, at least one (at least two when the variogram is estimated); numbers for one dimensional data or sequences of numbers
  - `ys` (sequence of numbers): the values at the nodes, as many as nodes
  - `variogram` (function, optional): the semivariogram `gamma(h)`, for example a model from `fastmath.kernel.variogram` (`(vgram/spherical {:nugget 0.0 :psill 1.0 :range 3.0})`). It has to be a valid semivariogram: `gamma(0) = 0` (a nugget is the limit at `0+`, not the value at `0`) and conditionally negative definite. When it is `nil` (the default), the empirical variogram of the data is estimated and a gaussian model is fitted to it by random search, see the options `:rng` and `:N`
  - `opts` (map, optional; a `nil` value is the same as a missing key, except for `:polynomial-terms`):
      - `:error` (double or sequence of doubles, default: `0.0`): the standard deviation of the measurement error of the values, the same for all the nodes or one per node. The variance `error^2` is subtracted from the diagonal of the variogram matrix, so the interpolator smooths the data instead of passing through the nodes. It is the `Err` component of the gstat variogram (its value is the variance there, so `Err = error^2`). An empty sequence means no error
      - `:distance` (function, default: Euclidean distance, `fastmath.distance/euclidean-1d` for numbers and `fastmath.distance/euclidean` for sequences): a binary distance function `(fn [x1 x2])`; it is used for the kriging system and also for the estimation of the empirical variogram of the default variogram
      - `:polynomial-terms` (function, default: a constant, `(constantly [1.0])`): a function of a point returning the sequence of the values of the drift basis, for example `(fn [[x y]] [1.0 x y])` for the linear drift of universal kriging; `nil` or `false` removes the drift terms and the unbiasedness constraint, the system is then the plain variogram matrix equation, which is rarely wanted
      - `:rng` (random number generator, default: a new unseeded generator for every call): the generator which draws the starting points of the variogram fit; it is used only when `variogram` is `nil`. Pass a seeded generator (see [[fastmath.random/rng]]) to get the same interpolator on every construction; the generator is consumed
      - `:N` (long, default: `100`): the number of the starting points scanned by the fit of the default variogram, increase it for a hard fit; used only when `variogram` is `nil`

  Returns a function of a query point (a number for one dimensional data, a sequence of numbers otherwise) returning a double. It does not return the kriging variance.

  Repeated nodes make the system singular: the minimum norm solution of the SVD is used silently, and the value at a repeated node is the mean of the values there (gstat gives `NA`). The default variogram is the best gaussian model for the binned estimate and is only a starting point: a visual check of the variogram or an explicit model fitted with `fastmath.kernel.variogram/fit` is recommended.

  Throws `ex-info` (with `:type :invalid-input`, `:method :kriging` and the key `:property`) when `xss` and `ys` have different lengths (`:length-mismatch`), there are too few nodes (`:too-few-points`) or `:error` is not a finite number or a sequence of finite numbers with one entry per node (`:error`). The estimation of the default variogram can throw its own `ex-info`, for example when all the nodes are at one location.

  See also [[fastmath.interpolation/kriging]], [[fastmath.kernel.variogram/fit]], [[fastmath.kernel.variogram/empirical]], [[fastmath.interpolation.rbf/rbf]], [[fastmath.interpolation.gp/gp]]."
  ([xss ys] (kriging xss ys nil))
  ([xss ys variogram] (kriging xss ys variogram nil))
  ([xss ys variogram opts]
   (let [{:keys [error distance rng N]} (into {} (remove (fn [[_ v]] (nil? v))) opts)
         polynomial-terms (if (contains? opts :polynomial-terms)
                            (:polynomial-terms opts)
                            (constantly [1.0]))
         _ (check-data! xss ys (if variogram 1 2))
         xss (vec xss)
         ysize (count xss)
         ;; the default distance works on arrays: numbers, vectors, Vec2 and lists can be mixed
         pts (if distance xss (->points xss))
         dist (or distance d/euclidean)
         ->query (if distance identity (query-converter (alength ^doubles (first pts))))
         variogram (or variogram
                       (vgram/fit (vgram/empirical (if distance xss (mapv vec pts)) ys {:distance (or distance d/euclidean)})
                                  :gaussian
                                  (cond-> {}
                                    rng (assoc :rng rng)
                                    N (assoc :N N))))
         V (v-matrix pts variogram dist error)
         [V ys] (if polynomial-terms (ic/pi-matrix-ys xss V ys polynomial-terms) [V ys])
         ysv (v/vec->RealVector ys)
         w (ic/solve-LU-or-SVD V ysv)
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
                                         (double (variogram (dist (pts i) q))))))
                    sum)))]
         (if c
           (m/+ s (v/dot c (polynomial-terms x)))
           s))))))
