(ns fastmath.ml.regression
  "Linear and generalized linear regression models with analysis, categorical encodings and isotonic regression.

  [[lm]] fits ordinary (OLS) and weighted (WLS) least squares models, [[glm]] fits generalized linear models with the IRLS method for the families and links stored in [[families]] and [[links]] (custom ones are created with [[->family]] and [[->link]]), [[glm-nbinomial]] additionally estimates theta of the negative binomial family. Both fit functions can use an offset, weights, a transformation of the terms (see [[polynomial-transformer]], [[trigonometric-transformer]] and [[b-spline-transformer]]) and regularization by data augmentation ([[add-penalty]]).

  A fitted model is a record which can be called as a function (see [[predict]]). [[analysis]] returns leverage, residuals and influence measures, [[quantile-residuals]] the randomized quantile residuals of a GLM and [[dose]] the lethal or effective dose of a binomial model. Categorical terms are encoded with `fastmath.ml.regression.contrast`.

  Isotonic regression is available as [[pava]] and [[cir]]."
  (:require [fastmath.core :as m]
            [fastmath.stats :as stats]
            [fastmath.vector :as v]
            [fastmath.protocols :as prot]
            [fastmath.random :as r]
            [fastmath.matrix :as mat]
            [fastmath.special :as special]
            [clojure.pprint :as pprint])
  (:import [org.apache.commons.math3.linear SingularValueDecomposition DiagonalMatrix
            QRDecomposition RealMatrix RealVector DecompositionSolver CholeskyDecomposition
            SingularMatrixException NonPositiveDefiniteMatrixException]
           [fastmath.java Array Monotone]
           [fastmath.vector Vec2]
           [clojure.lang IFn ArityException]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn with-offset
  "Splits a data point into an offset and terms.

  Parameters:

  - `xs` - a data point: a sequence of terms, `nil` or `[]` for a model without terms
  - `offset?` - if `true` the first element of `xs` is the offset

  Returns a pair `[offset terms]` where the terms are a vector. The offset is `0.0` when `offset?` is `false`.

  Throws `ex-info` with the `:offset` when `offset?` is `true` and the first element is missing or is not a number.

  See also [[predict]]."
  [xs offset?]
  (let [xs (v/vec->vector (or xs []))]
    (if offset?
      (let [offset (first xs)]
        (when-not (number? offset) (throw (ex-info "Offset should be a number." {:offset offset})))
        [(xs 0) (subvec xs 1)])
      [0.0 xs])))

(defn predict
  "Predict from the given model and data point.

  If `stderr?` is true, a map is returned with the keys `:fit`, `:stderr` (the standard error of the fit) and `:confidence-interval` (of the fit). For a GLM `:fit` is the mean, `:link` the linear predictor and the interval is calculated on the link scale and transformed back. For `lm` the map has also `:prediction-interval`, the interval of a new observation with the weight 1 (it is wider by the variance of the residual `sigma2`).
  If model is fitted with offset, first element of data point should contain provided offset.

  Expected data point:

  * `[x1,x2,...,xn]` - when model was trained without offset
  * `[offset,x1,x2,...,xn]` - when offset was used for training
  * `[]` or `nil` - when model was trained with intercept only
  * `[offset]` - when model was trained with intercept and offset"
  ([model xs] (prot/predict model xs))
  ([model xs stderr?] (prot/predict model xs stderr?)))

(declare ->string)

(defn- arity-error
  "A model is called with a data point and optionally `stderr?`."
  [model n]
  (throw (ArityException. (int n) (.getName (class model)))))

(defn- apply-model
  "Implementation of `applyTo` of a model: a data point and optionally `stderr?`."
  [^IFn model args]
  (case (count args)
    1 (.invoke model (first args))
    2 (.invoke model (first args) (second args))
    (arity-error model (count args))))

(defrecord LMData [model intercept? offset? transformer
                   ^RealMatrix xtxinv
                   ^double intercept beta coefficients
                   offset weights residuals fitted df ^long observations ^long nobs
                   names ^double cv
                   ^double r-squared ^double adjusted-r-squared ^double sigma2 ^double sigma
                   ^double tss ^double rss ^double regss ^double msreg
                   ^double qt
                   ^double f-statistic ^double p-value
                   ll analysis decomposition augmentation]
  Object
  (toString [this] (->string this))
  IFn
  (invoke [m xs] (prot/predict m xs false))
  (invoke [m xs stderr?] (prot/predict m xs stderr?))
  (invoke [m] (arity-error m 0))
  (invoke [m _a _b _c] (arity-error m 3))
  (invoke [m _a _b _c _d] (arity-error m 4))
  (invoke [m _a _b _c _d _e] (arity-error m 5))
  (invoke [m _a _b _c _d _e _f] (arity-error m 6))
  (invoke [m _a _b _c _d _e _f _g] (arity-error m 7))
  (invoke [m _a _b _c _d _e _f _g _h] (arity-error m 8))
  (invoke [m _a _b _c _d _e _f _g _h _i] (arity-error m 9))
  (invoke [m _a _b _c _d _e _f _g _h _i _j] (arity-error m 10))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k] (arity-error m 11))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l] (arity-error m 12))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l _m] (arity-error m 13))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l _m _n] (arity-error m 14))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l _m _n _o] (arity-error m 15))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l _m _n _o _p] (arity-error m 16))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l _m _n _o _p _q] (arity-error m 17))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l _m _n _o _p _q _r] (arity-error m 18))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l _m _n _o _p _q _r _s] (arity-error m 19))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l _m _n _o _p _q _r _s _t] (arity-error m 20))
  (applyTo [m args] (apply-model m args))
  prot/PredictProto
  (predict [m xs] (prot/predict m xs false))
  (predict [_ xs stderr?]
    (let [[^double off xs] (with-offset xs offset?)
          xs (if transformer (transformer xs) xs)]
      (when-not (m/== (count beta) (count xs))
        (throw (ex-info (str "Prediction expects " (count beta) " number of entries "
                             (when offset? "and offset ") "as input (after transformation).") {:offset? offset? :expected (count beta) :received (count xs)})))
      (if stderr?
        (let [arr (m/seq->double-array (if intercept? (cons 1.0 xs) xs))
              fit (m/+ off intercept (v/dot beta xs))
              stderr (m/sqrt (m/* sigma2 (v/dot arr (mat/mulv xtxinv arr))))
              scale (m/* stderr qt)
              pscale (m/* qt (m/sqrt (m/+ sigma2 (m/* stderr stderr))))]
          {:fit fit
           :stderr stderr
           :confidence-interval [(m/- fit scale) (m/+ fit scale)]
           :prediction-interval [(m/- fit pscale) (m/+ fit pscale)]})
        (m/+ off intercept (v/dot beta xs))))))

;;

(defn- add-diffs
  [xss {:keys [^double lambda ^long order]
        :or {lambda 1.0 order 2}}]
  (let [s (count (first xss))
        _ (when-not (m/< order s)
            (throw (ex-info "Order of differences should be smaller than the number of terms (with intercept)."
                            {:option :order :value order :expected (str "< " s)})))
        augm (-> (mat/eye s true)
                 (mat/differences order)
                 (mat/muls (m/sqrt lambda))
                 (mat/rows))]
    (concat xss (mapv v/vec->array augm))))

(defn- add-ridge
  [xss {:keys [^double lambda]
        :or {lambda 0.1}} intercept?]
  (let [s (count (first xss))
        zeros (repeat s 0.0)
        diag (vec (concat zeros (conj zeros (m/sqrt lambda))))
        icval (if intercept? -1 0)
        aug (map (fn [^long i]
                   (let [id (m/long-add icval (m/long-sub s i))]
                     (subvec diag id (m/long-add s id)))) (range (m/+ icval s)))]
    (vec (concat xss (if intercept? [zeros] '()) aug))))

(defn- add-one-penalty
  [xss penalty penalty-param intercept?]
  (case penalty
    :ridge (add-ridge xss penalty-param intercept?)
    :diffs (add-diffs xss penalty-param)
    xss))

(defn add-penalty
  "Adds rows with penalty data to a given seq of seqs (a row matrix).

  Penalties and parameters map:

  * `:ridge` - diagonal matrix with sqrt(`:lambda`) parameter, used in regularized regression.
  * `:diffs` - differences of given `:order` multiplied by sqrt(`:lambda`), used in penalized b-splines.

  Two add more than one penalty, create seq of penalty types and parameters, eg: `[:diffs :ridge]`"
  ([xss penalty penalty-param] (add-penalty xss penalty penalty-param false))
  ([xss penalty penalty-param intercept?]
   (if (and (sequential? penalty)
            (sequential? penalty-param))
     (reduce (fn [buff [p pp]]
               (add-one-penalty buff p pp intercept?)) xss (map vector penalty penalty-param))
     (add-one-penalty xss penalty penalty-param intercept?))))

(defn- invalid-option
  [option value expected]
  (throw (ex-info (str "Invalid " (name option) ", expected " expected ".")
                  {:option option :value value :expected expected})))

(defn- finite-number?
  [x]
  (and (number? x) (m/valid-double? (double x))))

(defn- first-invalid
  "Returns `[x]` for the first element `x` which does not satisfy `pred`, `nil` when all of them do (`x` can be `nil`)."
  [pred xs]
  (reduce (fn [_ x] (when-not (pred x) (reduced [x]))) nil xs))

(defn- check-terms
  "Converts rows of terms to vectors, applies the `transformer` (when given) and checks if all rows (after the
  transformation) have the same length and contain finite numbers. Returns a vector of vectors, or `:intercept`."
  [xss transformer]
  (if (= :intercept xss)
    xss
    (let [rows (vec (map-indexed
                     (fn [i row]
                       (let [terms (try (v/vec->seq row)
                                        (catch IllegalArgumentException e
                                          (throw (ex-info "A row of terms should be a number or a sequence of numbers."
                                                          {:option :xss :row i :value row} e))))]
                         (vec (if transformer (transformer terms) terms))))
                     xss))
          p (count (first rows))]
      (doseq [[i row] (map-indexed vector rows)]
        (when-not (m/== (count row) p)
          (throw (ex-info "All rows of terms should have the same length."
                          {:option :xss :row i :value (count row) :expected p})))
        (doseq [x row]
          (when-not (finite-number? x)
            (throw (ex-info "Terms should be finite numbers." {:option :xss :row i :value x :expected "a finite number"})))))
      rows)))

(defn- check-response
  "Checks if the response is a finite number, the family can require more (see [[check-family-response]])."
  [ys]
  (doseq [[i y] (map-indexed vector ys)]
    (when-not (finite-number? y)
      (throw (ex-info "Response should be finite numbers." {:option :ys :row i :value y :expected "a finite number"})))))

(defn- check-fit-options
  "Checks the options shared by `lm` and `glm`."
  [{:keys [tol alpha decomposition augmentation augmentation-params]}]
  (when-not (and (finite-number? tol) (m/>= (double tol) 0.0))
    (invalid-option :tol tol "a finite number >= 0"))
  (when-not (and (finite-number? alpha) (m/< 0.0 (double alpha) 1.0))
    (invalid-option :alpha alpha "a number between 0 and 1"))
  (when-not (#{:cholesky :qr :rrqr} decomposition)
    (invalid-option :decomposition decomposition ":cholesky, :qr or :rrqr"))
  (let [multi? (sequential? augmentation)
        penalties (cond (nil? augmentation) [] multi? augmentation :else [augmentation])
        params (if multi? augmentation-params [augmentation-params])]
    (when (and multi? (not (and (sequential? params) (m/== (count params) (count penalties)))))
      (invalid-option :augmentation-params augmentation-params "a sequence of the maps, one for each penalty"))
    (doseq [[penalty pp] (map vector penalties (if multi? params (repeat nil)))
            :let [pp (if multi? pp augmentation-params)]]
      (when-not (#{:ridge :diffs} penalty)
        (invalid-option :augmentation penalty "nil, :ridge, :diffs or a sequence of them"))
      (when-let [lambda (:lambda pp)]
        (when-not (and (finite-number? lambda) (m/>= (double lambda) 0.0))
          (invalid-option :lambda lambda "a finite number >= 0")))
      (when (and (= penalty :diffs) (contains? pp :order))
        (let [order (:order pp)]
          (when-not (and (finite-number? order) (m/== (double order) (double (long order))) (m/>= (double order) 1.0))
            (invalid-option :order order "an integer >= 1")))))))

(defn- check-data
  "Checks the lengths of the data and the values of `weights` and `offset`; nil `weights` and `offset` are missing options."
  [ys xss weights offset intercept?]
  (let [n (count ys)
        _ (when (m/zero? n)
            (throw (ex-info "Response should not be empty." {:option :ys :value n :expected "at least one observation"})))
        check-length (fn [option v]
                       (when (and (some? v) (not (m/== (count v) n)))
                         (throw (ex-info (str "Length of " (name option) " should be the same as the length of the response.")
                                         {:option option :value (count v) :expected n}))))]
    (when-not (= :intercept xss)
      (check-length :xss xss))
    (when (and (not intercept?)
               (or (= :intercept xss) (and (seq xss) (empty? (v/vec->seq (first xss))))))
      (throw (ex-info "Model has no terms." {:intercept? intercept?})))
    (check-length :weights weights)
    (check-length :offset offset)
    (when-let [[w] (first-invalid (fn [w] (and (finite-number? w) (m/>= (double w) 0.0))) weights)]
      (throw (ex-info "Weights should be finite and non-negative." {:option :weights :value w :expected "finite, >= 0"})))
    (when-let [[o] (first-invalid finite-number? offset)]
      (throw (ex-info "Offset should be finite." {:option :offset :value o :expected "finite"})))))

(defn- enhance-or-fill-seq
  [s ^long size ^double v]
  (take size (concat s (repeat v))))

(defn- augment-data
  "Add intercept and possible penalty"
  [ys xss weights offset intercept? augmentation augmentation-params]
  (let [xss (if (= xss :intercept) (repeat (count ys) '()) xss)
        xss (->  (if intercept? (map (fn [xs] (conj (seq xs) 1.0)) xss) xss)
                 (add-penalty augmentation augmentation-params intercept?))
        size (count xss)]
    [(enhance-or-fill-seq ys size 0.0)
     (mat/mat xss)
     (enhance-or-fill-seq weights size 1.0)
     (enhance-or-fill-seq offset size 0.0)]))

;;

(defn- small-singular-values-count
  ^long [^doubles sv ^double epsilon]
  (let [f (Array/aget sv 0)]
    (count
     (filter (fn [^double v] (m/< (m// v f) epsilon)) sv))))

(defn- svd
  [^RealMatrix xss tol]
  (let [^SingularValueDecomposition S (SingularValueDecomposition. xss)
        singular-values (.getSingularValues S)
        k (small-singular-values-count singular-values tol)]

    (if (m/pos? k)
      (throw (ex-info "Near rank-deficient model matrix"
                      {:rank (m/- (alength singular-values) k) :columns (mat/ncol xss)
                       :tol tol :singular-values singular-values}))
      [S (-> singular-values v/reciprocal v/vec->RealVector)])))

(defn- solver
  ^DecompositionSolver [^RealMatrix mat decomposition ^double tol]
  (:solver (case decomposition
             :qr (mat/qr-decomposition mat tol)
             :rrqr (mat/rrqr-decomposition mat tol)
             ;; the normal matrix is symmetric only up to rounding, the relative check fails for entries close to zero
             (mat/cholesky-decomposition (mat/symmetrize mat)
                                         CholeskyDecomposition/DEFAULT_RELATIVE_SYMMETRY_THRESHOLD
                                         1.0e-16))))

(defn- column-scales
  "Euclidean norms of the columns of the model matrix (`1.0` for a zero column). The columns are divided by them
  before the rank is checked and the system is solved, so the result does not depend on the units of the terms."
  ^doubles [^RealMatrix xss]
  (let [rows (.getRowDimension xss)
        cols (.getColumnDimension xss)
        scales (double-array cols)]
    (dotimes [j cols]
      (let [^double ss (loop [i (long 0) ss 0.0]
                 (if (m/< i rows)
                   (recur (m/inc i) (m/+ ss (m/sq (.getEntry xss i j))))
                   ss))]
        (Array/aset scales j (if (m/pos? ss) (m/sqrt ss) 1.0))))
    scales))

(defn- scale-columns
  "Divides the columns of the model matrix by `scales`."
  ^RealMatrix [^RealMatrix xss ^doubles scales]
  (let [rows (.getRowDimension xss)
        cols (.getColumnDimension xss)
        ^"[[D" arr (make-array Double/TYPE rows cols)]
    (dotimes [i rows]
      (let [^doubles row (aget arr i)]
        (dotimes [j cols]
          (Array/aset row j (m// (.getEntry xss i j) (Array/aget scales j))))))
    (mat/mat arr)))

(defn- unscale-coefficients
  "Coefficients of the original model matrix from the coefficients of the model matrix with scaled columns."
  ^doubles [coefficients ^doubles scales]
  (v/ediv (v/vec->array coefficients) scales))

(defn- mean-positive
  "Mean of the positive weights, `1.0` when there is none. Weights are divided by it before the system is solved,
  so the result does not depend on the scale of the weights."
  ^double [^doubles weights]
  (let [n (alength weights)]
    (loop [i (long 0) sum 0.0 cnt (long 0)]
      (if (m/< i n)
        (let [w (Array/aget weights i)]
          (if (m/pos? w)
            (recur (m/inc i) (m/+ sum w) (m/inc cnt))
            (recur (m/inc i) sum cnt)))
        (if (m/pos? cnt) (m// sum cnt) 1.0)))))

(defn- normalized-weights
  "Weights divided by their mean: `[diagonal-matrix scale]`."
  [^doubles weights]
  (let [scale (mean-positive weights)]
    [(DiagonalMatrix. ^doubles (v/mult weights (m// scale))) scale]))

(defn- xtxinv
  "Inverse of `X^T W X` of the original model matrix, from the matrix with scaled columns `xss` (see [[column-scales]])."
  ^RealMatrix [^RealMatrix xss ^doubles scales ^doubles weights ^double tol]
  (let [[^DiagonalMatrix dweights ^double wscale] (normalized-weights weights)
        ^RealMatrix inv (-> (mat/tmulm xss (mat/mulm dweights xss))
                            (QRDecomposition. tol)
                            (.getSolver)
                            (.getInverse))
        p (alength scales)
        ^"[[D" arr (make-array Double/TYPE p p)]
    (dotimes [i p]
      (let [^doubles row (aget arr i)]
        (dotimes [j p]
          (Array/aset row j (m// (.getEntry inv i j)
                                 (m/* (Array/aget scales i) (Array/aget scales j) wscale))))))
    (mat/mat arr)))

(defn- standard-errors
  [^RealMatrix xtxinv ^double dispersion]
  (-> xtxinv
      (mat/diag)
      (v/mult dispersion)
      (v/sqrt)
      (v/vec->seq)))

;; zero weights

(defn- count-positive
  "Number of positive weights, the number of observations used by a model."
  ^long [weights]
  (long (count (filter (fn [^double w] (m/pos? w)) weights))))

(defn- positive-only
  "Keeps the elements of `xs` with a positive weight."
  [xs weights]
  (->> (map vector xs weights)
       (filter (fn [[_ ^double w]] (m/pos? w)))
       (map first)))

(defn- mask-zero-weights
  "Replaces the entries of `xs` with `##NaN` where the weight is zero."
  [xs weights]
  (map (fn [x ^double w] (if (m/pos? w) x ##NaN)) xs weights))

(defn- mask-zero-weights-rows
  "As [[mask-zero-weights]], for a sequence of sequences (one per coefficient)."
  [xss weights]
  (map #(mask-zero-weights % weights) xss))

(defn- normality-analysis
  [{:keys [weighted raw]} weights]
  (let [weighted (positive-only weighted weights)
        raw (positive-only raw weights)
        skew (stats/skewness weighted :g1)
        kurt (stats/kurtosis weighted :kurt)
        normality (stats/normality-test weighted skew (- kurt 3.0) nil)]
    {:skewness skew
     :kurtosis kurt
     :durbin-watson {:raw (stats/durbin-watson raw)
                     :weighted (stats/durbin-watson weighted)}
     :jarque-berra (dissoc (stats/jarque-bera-test weighted skew kurt nil) :skewness :kurtosis)
     :omnibus (dissoc normality :skewness :kurtosis)}))

(defn- coefficients-analysis
  [coefficients stderrs distr ^double q]
  (let [label (if (= (r/distribution-id distr) :normal) :z-value :t-value)]
    (mapv (fn [^double b ^double err]
            (let [-value (m// b err)
                  scale (m/* q err)]
              {:estimate b :stderr err label -value
               :p-value (stats/p-value distr -value :both)
               :confidence-interval [(m/- b scale)
                                     (m/+ b scale)]})) coefficients stderrs)))

(defn- coefficients-names
  [names intercept? ^long term-count]
  (let [namesv (if-not names "X" names)
        term-count (if intercept? (m/inc term-count) term-count)
        step1 (cond-> []
                intercept? (conj "Intercept")
                (string? namesv) (into (take term-count (map #(str namesv "_" %) (range))))
                (sequential? namesv) (into namesv))
        curr-count (count step1)]
    (if (m/< curr-count term-count)
      (into step1 (take (m/long-sub term-count curr-count) (map #(str "X_" %) (range))))
      (subvec step1 0 term-count))))

;; new version

(defn- hat-matrix
  [^RealMatrix xtxinv ^RealMatrix xss weights ^long m]
  (let [sqrt-weights (m/seq->double-array (map m/sqrt weights))
        ^RealMatrix whx (mat/mulm (DiagonalMatrix. sqrt-weights) xss)]
    (->> (mat/rows whx)
         (map (fn [row] (let [h (v/dot row (mat/mulv xtxinv row))]
                         (m/min 1.0 h))))
         (take m))))

(defn- lm-transform-residuals
  [residuals hat sigmas]
  (map (fn [^double wr ^double h ^double sigma]
         (m// wr (m/* sigma (m/sqrt (m/- 1.0 h))))) residuals hat sigmas))

(defn- lm-cooks-distance
  [residuals hat ^double sigma ^long p]
  (map (fn [^double r ^double h]
         (m// (m/* h (m/sq (m// r (m/* (m/- 1.0 h) sigma)))) p)) residuals hat))

(defn- dffits
  [residuals hat sigmas]
  (map (fn [^double r ^double h ^double s]
         (if (m/< h 1.0)
           (m/* r (m// (m/sqrt h) (m/* s (m/- 1.0 h))))
           ##NaN)) residuals hat sigmas))

(defn- dfbetas
  [^RealMatrix leverage-coeffs ^RealMatrix xtx-1 sigmas]
  (let [d (v/vec->array (mat/diag xtx-1))]
    (map (fn [^RealVector v ^double s]
           (v/ediv (v/vec->seq v) (v/mult sigmas (m/sqrt s)))) (mat/rows leverage-coeffs) d)))

(defn- covratio
  [residuals hat sigmas p n]
  (let [p (long p)
        n-p (m/- (long n) p)
        n-p-1 (m/dec n-p)]
    (map (fn [^double r ^double h ^double s]
           (if (and (m/valid-double? s) (m/< h 1.0))
             (let [h- (m/- 1.0 h)]
               (m// (m/* h- (m/fpow (m// (m/+ n-p-1 (m/sq (m// r (m/* s (m/sqrt h-)))))
                                         n-p) p))))
             ##NaN))
         residuals hat sigmas)))

(defn- influential-rows
  [what how]
  (->> (map-indexed vector what)
       (filter (comp how second))
       (map first)))

(defn- measures-influential-rows
  [influence hat ^long n ^long p]
  (let [n-p (m/- n p)]
    (if (m/not-pos? n-p)
      {:cooks-distance '() :dffits '() :covratio '() :hat '()
       :dfbetas (repeat (count (:dfbetas influence)) '()) :combined {}}
      (let [pf (r/distribution :f {:numerator-degrees-of-freedom p
                                   :denominator-degrees-of-freedom n-p})
            dffits-v (m/* 3.0 (m/sqrt (m// (double p) n-p)))
            covratio-v (m// (m/* 3.0 p) n-p)
            hat-v (m// (* 3.0 p) n)
            res {:cooks-distance (influential-rows (:cooks-distance influence)
                                                   (fn [^double x] (m/> (r/cdf pf x) 0.5)))
                 :dffits (influential-rows (:dffits influence) (fn [^double x] (m/> (m/abs x) dffits-v)))
                 :covratio (influential-rows (:covratio influence) (fn [^double x] (m/> (m/abs (m/- 1.0 x)) covratio-v)))

                 :hat (influential-rows hat (fn [^double x] (m/> x hat-v)))
                 :dfbetas (map influential-rows (:dfbetas influence) (repeat (fn [^double x] (m/> (m/abs x) 1.0))))}]
        (assoc res :combined (frequencies (flatten (vals res))))))))

(defn- correlation
  [^RealMatrix xtxinv ^double sigma2]
  (let [r (mat/muls xtxinv sigma2)
        se (mat/rows->RealMatrix [(v/sqrt (mat/diag r))])
        outer (mat/tmulm se se)]
    (map (comp v/vec->seq v/ediv) (mat/rows r) (mat/rows outer))))

(defn- sigmas
  [residuals hat ^double rss ^long df-]
  (map (fn [^double wr ^double h]
         (if (m/< h 1.0)
           (m/sqrt (m// (m/- rss (m// (m/* wr wr)
                                      (m/- 1.0 h))) df-))
           (m/safe-sqrt (m// rss df-)))) residuals hat))

(defn- leverage-coeffs
  ^RealMatrix [xtxinv xss residuals weights hat p observations]
  (let [^RealMatrix bb (mat/mulm xtxinv (mat/transpose xss))]
    (doseq [^long row (range p)
            [^long col ^double f] (map (fn [^long id ^double r ^double h ^double w]
                                         [id (if (m/< h 1.0)
                                               (m// (m/* r w)
                                                    (m/- 1.0 h))
                                               0.0)])
                                       (range observations) residuals hat weights)]
      (.multiplyEntry bb row col f))
    bb))

(defn- lm-extra-analysis
  [{:keys [weights hat residuals ^double rss ^long observations ^long nobs
           ^double sigma ^double sigma2 ^RealMatrix xtxinv] :as model} ^RealMatrix xss]
  (let [rresiduals (:raw residuals)
        wresiduals (:weighted residuals)
        {^double df :residual} (:df model)
        p (mat/ncol xss)
        
        df- (m/dec df)

        sigmas (sigmas wresiduals hat rss (long df-))        
        
        ^RealMatrix leverage-coeffs (leverage-coeffs xtxinv xss rresiduals weights hat p observations)

        ;; rows with zero weight do not take part in the model, their diagnostics are NaN
        masked (fn [xs] (mask-zero-weights xs weights))
        influence {:cooks-distance (masked (lm-cooks-distance wresiduals hat sigma p))
                   :dffits (masked (dffits wresiduals hat sigmas))
                   :dfbetas (mask-zero-weights-rows (dfbetas leverage-coeffs xtxinv sigmas) weights)
                   :covratio (masked (covratio wresiduals hat sigmas p nobs))}]
    {:normality (normality-analysis residuals weights)
     :residuals {:standardized (masked (lm-transform-residuals wresiduals hat (repeat sigma)))
                 :studentized (masked (lm-transform-residuals wresiduals hat sigmas))}
     :leverage {:hat hat
                :sigmas (masked sigmas)
                :coefficients (mask-zero-weights-rows (map seq (mat/mat->array2d leverage-coeffs)) weights)}
     :influence influence
     :influential (measures-influential-rows influence hat nobs p)
     :correlation (correlation xtxinv sigma2)}))

(defn lm
  "Fit a linear model using ordinary (OLS) or weighted (WLS) least squares.

  Arguments:

  * `ys` - response vector
  * `xss` - terms of systematic component
  * optional parameters

  Parameters:

  * `:tol` - tolerance of the matrix decomposition: the relative threshold of small singular values (rank deficiency of the model matrix) and the threshold of the QR/RRQR decomposition, default: `1.0e-8`
  * `:weights` - optional weights for WLS, finite and non-negative; an observation with the weight 0 stays in `:fitted` and `:residuals` but it does not take part in the model (it is not counted in `:nobs` and degrees of freedom)
  * `:offset` - optional offset
  * `:alpha` - significance level, a number between 0 and 1, default: `0.05`
  * `:intercept?` - should intercept term be included, default: `true`
  * `:transformer` - an optional function which will be used to transform systematic component `xs` before fitting and prediction
  * `:names` - sequence or string, used as name for coefficient when pretty-printing model, default `'X'`
  * `:decomposition` - which matrix decomposition use to find solution, `:cholesky` (default), `:rrqr` (rank revealing) or `:qr`
  * `:augmentation` and `augmentation-params` - regularization by data augmentation
      - `:ridge` - adds ridge regresion penalty (intercept is not penalized), default parameters `{:lambda 0.1}` 
      - `:diffs` - adds differences penalty, use with `b-spline-transformation` for smoothing, default parameters `{:lambda 1.0 :order 2}` 

  Notes:

  * SVD decomposition is used instead of more common QR
  * intercept term is added implicitely if `intercept?` is set to `true` (by default)
  * `xss` equal to `:intercept` fits a model with the intercept only: `:f-statistic` and `:p-value` are `##NaN` and R-squared is `0.0`
  * Two variants of AIC/BIC are calculated, one based on log-likelihood, second on RSS/n

  The fit does not depend on the units of the terms nor on the scale of the weights: columns of the model matrix and weights are scaled before the rank is checked and the system is solved.

  Throws `ex-info` when an option is invalid (`:option`, `:value` and `:expected` in the data), when the response is empty or not finite, when the terms are not finite numbers or the rows have different lengths, when the lengths of `ys`, `xss`, `:weights` and `:offset` differ, when weights are negative or not finite, when the model has no terms, when the model matrix is rank deficient (`:rank`, `:columns`, `:tol` and `:singular-values` in the data) and when there are no residual degrees of freedom (`:observations`, `:parameters`).

  A constant response has no variance to explain: R-squared, adjusted R-squared, the F statistic and its p-value are `##NaN`.

  Returned record implementes `IFn` protocol and contains:

  * `:model` - `:ols` or `:wls`
  * `:intercept?` - whether intercept term is included or not
  * `:xtxinv` - (X^T X)^-1
  * `:intercept` - intercept term value
  * `:beta` - vector of model coefficients (without intercept)
  * `:coefficients` - coefficient analysis, a list of maps containing `:estimate`, `:stderr`, `:t-value`, `:p-value` and `:confidence-interval`
  * `:weights` - initial weights
  * `:residuals` - a map containing `:raw` and `:weighted` residuals, also `:loocv`
  * `:fitted` - fitted values for xss
  * `:df` - degrees of freedom: `:residual`, `:model` and `:intercept`
  * `:observations` - number of observations (rows)
  * `:nobs` - number of observations with a positive weight, used in degrees of freedom and information criteria
  * `:r-squared` and `:adjusted-r-squared`
  * `:sigma` and `:sigma2` - deviance and variance
  * `:msreg` - regression mean squared
  * `:rss`, `:regss`, `:tss` - residual, regression and total sum of squares
  * `:qt` - (1-alpha/2) quantile of T distribution for residual degrees of freedom
  * `:f-statistic` and `:p-value` - F statistic and respective p-value
  * `:ll` - a map containing log-likelihood and AIC/BIC in two variants: based on log-likelihood and RSS
  * `:analysis` - leverage, residual and influence analysis - a delay
  * `:decomposition` - decomposition used
  * `:augmentation` - augmentation used
  * `:cv` - cross validation statistic
  * `:effective-dimension` - effective dimension of the model (sum of hat matrix diagonal)
  
  Analysis, delay containing a map:

  * `:residuals` - `:standardized` and `:studentized` weighted residuals
  * `:leverage` - `:hat`, `:sigmas` and leveraged `:coefficients` (leave-one-out)
  * `:influence` - `:cooks-distance`, `:dffits`, `:dfbetas` and `:covratio`
  * `:influential` - list of influential observations (ids) for influence measures
  * `:correlation` - correlation matrix of estimated parameters
  * `:normality` - residuals normality tests: `:skewness`, `:kurtosis`, `:durbin-watson` (for raw and weighted), `:jarque-berra` and `:omnibus` (normality)"
  ([ys xss] (lm ys xss nil))
  ([ys xss {:keys [^double tol weights ^double alpha intercept? offset transformer names decomposition
                   augmentation augmentation-params]
            :or {tol 1.0e-8 alpha 0.05 intercept? true decomposition :cholesky}}]

   (check-fit-options {:tol tol :alpha alpha :decomposition decomposition
                       :augmentation augmentation :augmentation-params augmentation-params})
   (check-data ys xss weights offset intercept?)
   (check-response ys)
   (let [offset? (boolean offset)
         weights? (some? weights)
         xss (check-terms xss transformer)

         
         m (long (count ys))

         [ys ^RealMatrix xss weights offset] (augment-data ys xss weights offset intercept?
                                                           augmentation augmentation-params)


         n (mat/ncol xss)

         ;; observations with zero weight do not take part in the model
         nobs (count-positive (take m weights))
         _ (when (and (not augmentation) (m/<= nobs n))
             (throw (ex-info "No residual degrees of freedom." {:observations nobs :parameters n :df (m/- nobs n)})))

         ;; the rank and the system do not depend on the scale of terms and weights
         ^doubles scales (column-scales xss)
         ^RealMatrix xss-scaled (scale-columns xss scales)
         [^SingularValueDecomposition S singular-values] (svd xss-scaled tol)
         uts (.getUT S)
         ut (.getU S)

         ^doubles offset (m/seq->double-array offset)
         ^doubles ys-orig (m/seq->double-array ys)
         ^doubles ys (v/sub ys-orig offset)

         ^doubles daweights (m/seq->double-array weights)
         [^DiagonalMatrix dweights wscale] (normalized-weights daweights)
         ^RealVector rvwys (v/vec->RealVector (v/mult (v/emult daweights ys) (m// (double wscale))))


         ^RealVector result (let [new-t (->> (try (-> (mat/mulm uts (mat/mulm dweights ut))
                                                      ^DecompositionSolver (solver decomposition tol)
                                                      (.solve ^RealVector (mat/mulv uts rvwys)))
                                                  (catch SingularMatrixException e
                                                    (throw (ex-info "The normal equations are singular, try the :qr or :rrqr decomposition."
                                                                    {:decomposition decomposition} e)))
                                                  (catch NonPositiveDefiniteMatrixException e
                                                    (throw (ex-info "The normal equations are not positive definite, try the :qr or :rrqr decomposition."
                                                                    {:decomposition decomposition} e))))
                                             (mat/mulv ut))]

                              (->> singular-values
                                   (v/emult (mat/mulv uts new-t))
                                   (mat/mulv (.getV S))
                                   (#(unscale-coefficients % scales))
                                   (v/vec->RealVector)))

         ;; remove augmentation
         fitted (take m (v/vec->seq (mat/mulv xss result)))
         result-array (v/vec->array result)

         intercept (if intercept? (Array/aget result-array 0) 0.0)
         beta (vec (if intercept? (rest result-array) result-array))

         xtx-1 (xtxinv xss-scaled scales daweights tol)
         hat (hat-matrix xtx-1 xss weights m)

         ed (v/sum hat)

         ;; remove augmentation
         ys (take m ys)
         weights (take m weights)
         offset (take m offset)

         df (if augmentation (m/- nobs ed) (m/- nobs n))
         _ (when-not (m/pos? (double df))
             (throw (ex-info "No residual degrees of freedom." {:observations nobs :parameters n :df df})))
         model-df (count beta)
         intercept-df (long (if intercept? 1 0))

         raw-residuals (map m/- ys fitted)
         wresiduals (if weights?
                      (map (fn [^double r ^double w] (m/* r (m/sqrt w))) raw-residuals weights)
                      raw-residuals)
         
         loocv (map (fn [^double r ^double h] (m// r (m/- 1.0 h))) wresiduals hat)

         ;; loocv, cross validation 
         cv (m/sqrt (m// (v/sum (v/sq loocv)) nobs))

         rss (v/dot wresiduals wresiduals)
         sigma2 (m// rss (double df))
         sigma (m/sqrt sigma2)

         tss (if intercept?
               (if weights?
                 (let [mean-ys (stats/mean ys weights)]
                   (v/sum (map (fn [^double y ^double w] (m/* w (m/sq (m/- y mean-ys)))) ys weights)))
                 (let [mean-ys (stats/mean ys)]
                   (v/sum (map (fn [^double y] (m/sq (m/- y mean-ys))) ys))))
               (if weights?
                 (v/sum (map (fn [^double y ^double w] (m/* w y y)) ys weights))
                 (v/dot ys ys)))

         ;; a constant response has no variance to explain (the sum of squares is zero up to rounding)
         constant-response? (m/<= tss (m/* 1.0e-28 (if weights?
                                                     (v/sum (map (fn [^double y ^double w] (m/* w y y)) ys weights))
                                                     (v/dot ys ys))))
         ;; a model without terms (intercept only) explains nothing
         terms? (and (m/pos? model-df) (not constant-response?))
         regss (if terms? (m/- tss rss) 0.0)
         rsquared (double (cond constant-response? ##NaN
                                (m/pos? model-df) (m/- 1.0 (m// rss tss))
                                :else 0.0))
         adjusted-rsquared (m/- 1.0 (m/* (m// (m/- nobs intercept-df) (double df)) (m/- 1.0 rsquared)))
         ;; no F test for a model without terms (intercept only) and for a constant response
         f-statistic (if terms? (m// (m// regss model-df) sigma2) ##NaN)
         p-value (if terms?
                   (stats/p-value (r/distribution :f {:numerator-degrees-of-freedom model-df
                                                      :denominator-degrees-of-freedom df})
                                  f-statistic :one-sided-greater)
                   ##NaN)
         tdistr (r/distribution :t {:degrees-of-freedom df})
         qt (double (r/icdf tdistr (m/- 1.0 (m/* alpha 0.5))))

         ;; log-likelihood based measures
         ll (let [logm (m/log nobs)
                  nobs2 (m// nobs 2.0)
                  log-likelihood (m/+ (m/- (m/- (m/* nobs2 (m/log rss)))
                                           (m/* nobs2 (m/inc (m/log (m// m/PI nobs2)))))
                                      (m/* 0.5 (v/sum (v/log (filter m/pos? weights)))))
                  mlogsigma (m/* nobs (m/log (m// rss nobs)))]
              {:log-likelihood log-likelihood
               :aic (m/+ (m/* -2.0 log-likelihood) (* 2.0 (m/inc ed)))
               :bic (m/+ (m/* -2.0 log-likelihood) (* logm (inc ed)))
               :aic-rss (m/+ mlogsigma (m/* 2.0 ed))
               :bic-rss (m/+ mlogsigma (m/* logm ed))})

         stderrs (standard-errors xtx-1 sigma2)
         namev (coefficients-names names intercept? (count beta))

         model {:model (if weights? :wls :ols)
                :transformer transformer
                :intercept? intercept?
                :offset? offset?
                :xtxinv xtx-1
                :intercept intercept
                :beta beta
                :hat hat
                :offset offset
                :coefficients (coefficients-analysis result-array stderrs tdistr qt)
                :names namev
                :weights weights
                :residuals {:weighted wresiduals :raw raw-residuals :loocv loocv}
                :fitted (v/add fitted (seq offset))
                :df {:residual df :model model-df :intercept intercept-df}
                :observations m
                :nobs nobs
                :effective-dimension ed
                :cv cv
                :augmentation augmentation
                :r-squared rsquared
                :adjusted-r-squared adjusted-rsquared
                :sigma2 sigma2
                :sigma sigma
                :tss tss :rss rss :regss regss :msreg (m// regss model-df)
                :qt qt
                :f-statistic f-statistic :p-value p-value
                :ll ll
                :decomposition decomposition}

         analysis (delay (lm-extra-analysis model (if augmentation
                                                    (.getSubMatrix xss 0 (m/dec m) 0 (m/dec n))
                                                    xss)))]

     (map->LMData (assoc model :analysis analysis)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; GLM https://bwlewis.github.io/GLM/
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; link

;; .Machine$double.eps of R: the bounds used by the links to keep the mean in the domain of the family
(def ^{:private true :const true :tag 'double} link-eps 2.220446049250313E-16)
(def ^{:private true :const true :tag 'double} link-eps-upper 0.9999999999999998)

(defn- constantly-1 ^double [^double _] 1.0)

(defrecord Link [g mean derivative])

(defn ->link
  "Creates link record.

  Args:

  * `g` - link function
  * `mean` - mean, inverse link function
  * `mean-derivative` - derivative of mean"
  ([link-map] (map->Link link-map))
  ([g mean mean-derivative]
   (->Link g mean mean-derivative)))

(def ^:private log-link
  (->Link m/log
          (fn ^double [^double x] (m/max link-eps (m/exp x)))
          (fn ^double [^double x] (m/max link-eps (m/exp x)))))

(def ^{:doc "Link functions of GLM, a map from a keyword to a `Link` record (see [[->link]]) or to a function returning a `Link` from a map of parameters.

  Records: `:logit`, `:probit`, `:cauchit`, `:cloglog`, `:loglog`, `:identity`, `:log`, `:clog` (`log(1-mu)`), `:sqrt`, `:inversesq` (`1/mu^2`) and `:inverse`.

  Functions: `:power` (parameter `:power-exponent`, default `1.0`, the exponent `0.0` is the `:log` link), `:nbinomial` (parameter `:nbinomial-theta`, default `1.0`) and `:distribution` (parameter `:distribution`, the quantile function of a distribution as the link, default: standard normal).

  Use [[family-with-link]] to get a family merged with a link."}
  links {:logit (->Link m/logit
                         m/sigmoid
                         (fn ^double [^double x] (let [e (m/exp x)] (m// e (m/sq (m/inc e))))))
          :probit (->Link (fn ^double [^double x] (m/* m/SQRT2 (special/inv-erf (m/dec (m/* 2.0 x)))))
                          (fn ^double [^double x] (-> (m// x m/SQRT2)
                                                     (m/-)
                                                     (special/erfc)
                                                     (m/* 0.5)))
                          (fn ^double [^double x] (m/max link-eps (m/* m/INV_SQRT2PI
                                                                               (m/exp (m/* -0.5 x x))))))
          :cauchit (->Link (fn ^double [^double x] (m/tan (m/* m/PI (m/- x 0.5))))
                           (fn ^double [^double x] (-> x
                                                      (m/atan)
                                                      (m// m/PI)
                                                      (m/+ 0.5)))
                           (fn ^double [^double x] (m/* m/INV_PI (m// 1.0 (m/inc (m/* x x))))))
          :cloglog (->Link m/cloglog
                           (fn ^double [^double x] (m/constrain (m/cexpexp x)
                                                               link-eps
                                                               link-eps-upper))
                           (fn ^double [^double x] (m/max link-eps (m/exp (m/- x (m/exp x))))))
          :loglog (->Link m/loglog
                          m/expexp
                          (fn ^double [^double x] (m/max link-eps
                                                        (m/exp (m/- (m/- (m/exp (m/- x))) x)))))
          :identity (->Link m/identity-double
                            m/identity-double
                            constantly-1)
          :log log-link
          :clog (->Link (fn ^double [^double x] (m/log (m/- 1.0 x)))
                        (fn ^double [^double x] (m/- 1.0 (m/exp x)))
                        (fn ^double [^double x] (m/- (m/exp x))))
          :sqrt (->Link m/sqrt m/sq (fn ^double [^double x] (m/* 2.0 x)))
          :inversesq (->Link (fn ^double [^double x] (m// 1.0 (m/* x x)))
                             (fn ^double [^double x] (m// 1.0 (m/sqrt x)))
                             (fn ^double [^double x] (m// -1.0 (m/* 2.0 (m/pow x 1.5)))))
          :inverse (->Link m// m// (fn ^double [^double x] (m// -1.0 (m/sq x))))
          :nbinomial (fn [{:keys [^double nbinomial-theta]
                          :or {nbinomial-theta 1.0}}]
                       (->Link (fn ^double [^double x] (m/log (m// x (m/+ x nbinomial-theta))))
                               (fn ^double [^double x] (m// -1.0 (m// (m/- 1.0 (m/exp (m/- x)))
                                                                     nbinomial-theta)))
                               (fn ^double [^double x] (let [e (m/exp x)]
                                                        (m// e (m// (m/sq (m/- 1.0 e))
                                                                    nbinomial-theta))))))
          :power (fn [{:keys [^double power-exponent]
                      :or {power-exponent 1.0}}]
                   (if (m/zero? power-exponent)
                     ;; the limit of the power link for the exponent 0
                     log-link
                     (let [rexponent (m// power-exponent)]
                       (->Link (fn ^double [^double x] (m/pow x power-exponent))
                               (fn ^double [^double x] (m/pow x rexponent))
                               (fn ^double [^double x] (m/* rexponent (m/pow x (m/* rexponent
                                                                                   (m/- 1.0 power-exponent)))))))))
          :distribution (fn [{:keys [distribution]
                             :or {distribution (r/distribution :normal)}}]
                          (->Link (partial r/icdf distribution)
                                  (partial r/cdf distribution)
                                  (partial r/pdf distribution)))})

;; family

(defn- default-initialize [ys weights] [ys ys weights])

(defn- binomial-initialize
  [ys weights]
  (if (sequential? (first ys)) ;; check if we have pairs of successes and failures
    (let [;; rows added by data augmentation (regularization) are plain zeros: no success, one trial
          ys (map (fn [y] (if (sequential? y) y [0.0 1.0])) ys)
          n (map (fn [[^double s ^double f]] (m/+ s f)) ys)
          ;; no trials: the proportion is arbitrary, the weight is zero
          y (map (fn [[^double s _] ^double trials] (if (m/pos? trials) (m// s trials) 0.0)) ys n)]
      [y (map (fn [^double y ^double n]
                (m// (m/+ 0.5 (m/* n y))
                     (m/inc n))) y n)
       (v/emult weights n) n])
    [ys
     (map (fn [^double y ^double w]
            (m// (m/+ 0.5 (m/* w y))
                 (m/inc w))) ys weights)
     weights (repeat (count ys) 1.0)]))

(defn- binomial-log-pmf
  "Log of the binomial probability mass function; the number of trials can be outside of the int range of the distribution."
  ^double [^double k ^double n ^double p]
  (if (m/> n 2.0e9)
    (cond
      (or (m/< k 0.0) (m/> k n)) ##-Inf
      (m/<= p 0.0) (if (m/zero? k) 0.0 ##-Inf)
      (m/>= p 1.0) (if (m/== k n) 0.0 ##-Inf)
      :else (m/+ (m/- (special/log-gamma (m/inc n)) (special/log-gamma (m/inc k)) (special/log-gamma (m/inc (m/- n k))))
                 (m/* k (m/log p))
                 (m/* (m/- n k) (m/log1p (m/- p)))))
    (r/lpdf (r/distribution :binomial {:trials n :p p}) k)))

(defn- binomial-aic
  [ys fitted weights _deviance _observations rank ns]
  (let [ms (if (some (fn [^long n] (m/> n 1)) ns) ns weights)]
    (m/+ (m/* -2.0 (v/sum (map (fn [^double y ^double m ^double mu ^double wt]
                                 (if (m/not-pos? m)
                                   0.0
                                   ;; R rounds half to even
                                   (m/* (m// wt m)
                                        (binomial-log-pmf (m/rint (m/* m y)) (m/rint m) mu))))
                               ys ms fitted weights)))
         (m/* 2.0 (double rank)))))

(defn- binomial-quantile-residuals
  [{:keys [ys weights fitted]}]
  (map (fn [^double y ^double w ^double p]
         (let [distr (r/distribution :binomial {:trials w :p p})
               yw (m/* y w)
               a (r/cdf distr (m/dec yw))
               b (r/cdf distr yw)]
           (r/icdf r/default-normal (r/drand a b)))) ys (:initial weights) fitted))

(defn- ylogymu ^double [^double y ^double mu] (if (m/zero? y) 0.0 (m/* y (m/log (m// y mu)))))
(defn- binomial-residual-deviance
  [ys mus weights]
  (map (fn ^double [^double y ^double mu ^double w]
         (m/* 2.0 w (m/+ (ylogymu y mu)
                         (ylogymu (m/- 1.0 y)
                                  (m/- 1.0 mu))))) ys mus weights))

(defn- gaussian-residual-deviance
  [ys mus weights]
  (map (fn [^double y ^double mu ^double w]
         (m/* w (m/sq (m/- y mu))))
       ys mus weights))

(defn- gaussian-aic
  [_ys _fitted weights deviance observations rank _ns]
  (let [nobs (long observations)]
    (m/+ (m/* nobs (m/inc (m/log (m/* m/TWO_PI (m// (double deviance) nobs)))))
         2.0 (m/- (v/sum (v/log (filter m/pos? weights))))
         (m/* 2.0 (double rank)))))

(defn- gamma-residual-deviance
  [ys mus weights]
  (map (fn [^double y ^double mu ^double w]
         (m/* -2.0 w (m/- (m/log (if (m/zero? y) 1.0 (m// y mu)))
                          (m// (m/- y mu) mu))))
       ys mus weights))

(defn- gamma-aic
  [ys fitted weights deviance _observations rank _ns]
  (let [disp (m// (double deviance) (v/sum weights))
        rdisp (m// disp)]
    (m/+ (m/* -2.0 (v/sum (map (fn [^double y ^double mu ^double w]
                                 (m/* w (r/lpdf
                                         (r/distribution :gamma {:shape rdisp :scale (m/* mu disp)})
                                         y)))
                               ys fitted weights)))
         2.0 (m/* 2.0 (double rank)))))

(defn- gamma-quantile-residuals
  [{:keys [ys fitted ^double dispersion weights]}]
  (map (fn [^double y ^double mu ^double w]
         (->> (m// (m// (m/* w y) mu) dispersion)
              (r/cdf (r/distribution :gamma {:shape (m// w dispersion) :scale 1.0}))
              (r/icdf r/default-normal))) ys fitted (:initial weights)))

(defn- poisson-initialize
  [ys weights]
  [ys (v/shift ys 0.1) weights])

(defn- poisson-residual-deviance [ys mus weights]
  (map (fn [^double y ^double mu ^double w]
         (let [r (if (m/pos? y)
                   (m/-
                    (m/* y (m/log (m// y mu)))
                    (m/- y mu))
                   mu)]
           (m/* 2.0 w r)))
       ys mus weights))

(defn- poisson-aic
  [ys fitted weights _deviance _observations rank _ns]
  (m/+ (m/* -2.0 (v/sum (map (fn [^double y ^double mu ^double w]
                               (m/* w (r/lpdf (r/distribution :poisson {:p mu}) y))) ys fitted weights)))
       (m/* 2.0 (double rank))))

(defn- poisson-quantile-residuals
  [{:keys [ys fitted]}]
  (map (fn [^double y ^double mu]
         (let [distr (r/distribution :poisson {:p mu})
               a (r/cdf distr (m/dec y))
               b (r/cdf distr y)]
           (r/icdf r/default-normal (r/drand a b)))) ys fitted))

(defn- inverse-gaussian-residual-deviance
  [ys mus weights]
  (map (fn [^double y ^double mu ^double w]
         (m/* w (m// (m/sq (m/- y mu))
                     (m/* y mu mu))))
       ys mus weights))

(defn- inverse-gaussian-aic
  [ys _fitted weights deviance _observations rank _ns]
  (let [sumw (v/sum weights)
        d (double deviance)]
    (m/+ (m/* sumw (m/inc (m/log (m/* m/TWO_PI (m// d sumw)))))
         (m/* 3.0 (v/sum (map (fn [^double y ^double w] (m/* (m/log y) w)) ys weights)))
         2.0 (m/* 2.0 (double rank)))))

(defn- inverse-gaussian-quantile-residuals
  [{:keys [ys fitted dispersions]}]
  (let [rdispersion (m// (double (:mean-deviance dispersions)))]
    (map (fn [^double y ^double mu]
           (if (m/== y mu)
             y
             (let [distr (r/distribution :inverse-gaussian {:mu mu :lambda rdispersion})]
               (if (m/< y mu)
                 (r/icdf r/default-normal (r/cdf distr y))
                 (m/- (double (r/icdf r/default-normal (r/ccdf distr y)))))))) ys fitted)))

(defn- nbinomial-initialize
  [ys weights]
  [ys (map (fn ^double [^double y]
             (if (m/zero? y) m/SIXTH y)) ys) weights])

(defn- ->nbinomial-residual-deviance
  [^double theta]
  (fn [ys mus weights]
    (map (fn [^double y ^double mu ^double w]
           (let [y+theta (m/+ y theta)]
             (m/* 2.0 w (m/- (m/* y (m/log (m// (m/max y 1.0) mu)))
                             (m/* y+theta
                                  (m/log (m// y+theta (m/+ mu theta))))))))
         ys mus weights)))

(defn- ->nbinomial-quantile-residuals
  [^double theta]
  (fn [{:keys [ys fitted]}]
    (map (fn [^double y ^double mu]
           (let [p (m// theta (m/+ theta mu))
                 a (if (m/pos? y)
                     (r/cdf (r/distribution :beta {:alpha theta :beta (m/max y 1.0)}) p)
                     0.0)
                 b (r/cdf (r/distribution :beta {:alpha theta :beta (m/inc y)}) p)]
             (r/icdf r/default-normal (r/drand a b)))) ys fitted)))

(defn- ->nbinomial-aic
  [^double theta]
  (let [lgt-tlt (m/- (special/log-gamma theta) (m/* theta (m/log theta)))]
    (fn [ys fitted weights _deviance _observations rank _ns]
      (m/+ (m/* 2.0 (v/sum (map (fn [^double y ^double mu ^double w]
                                  (let [y+t (m/+ y theta)]
                                    (m/* w (m/- (m/+ (m/* y+t (m/log (m/+ mu theta)))
                                                     (special/log-gamma (m/inc y))
                                                     lgt-tlt)
                                                (m/* y (m/log mu))
                                                (special/log-gamma y+t))))) ys fitted weights)))
           (m/* 2.0 (m/inc (double rank)))))))

(defrecord Family [default-link variance initialize residual-deviance aic quantile-residuals-fun dispersion])

(defn ->family
  "Create `Family` record.

  Arguments:

  * `default-link` - canonical link function, default: `:identity`
  * `variance` - variance function in terms of mean
  * `initialize` - initialization of glm, default: the same as in `:gaussian`
  * `residual-deviance` - calculates residual deviance
  * `aic` - calculates AIC, default `(constantly ##NaN)`
  * `quantile-residuals-fun` - calculates quantile residuals, default as in `:gaussian`
  * `disperation` - value or `:estimate` (default), `:pearson` or `:mean-deviance`

  Initialization will be called with `ys` and `weights` and should return:

  * ys, possibly changed if any adjustment is necessary
  * init-mu, starting point
  * weights, possibly changes or orignal
  * (optional) any other data used to calculate AIC

  AIC function should accept: `ys`, `fitted`, `weights`, `deviance`, `observation`, `rank` (fitted parameters) and additional data created by initialization

  Minimum version should define `variance` and `residual-deviance`."
  ([family-map] (map->Family family-map))
  ([default-link variance initialize residual-deviance aic quantile-residuals-fun dispersion]
   (->Family default-link variance initialize residual-deviance aic quantile-residuals-fun dispersion))
  ([variance residual-deviance]
   (->Family :identity variance default-initialize residual-deviance (constantly ##NaN) nil :estimate)))

(def ^{:doc "Families of GLM, a map from a keyword to a `Family` record (see [[->family]]) or to a function returning a `Family` from a map of parameters.

  Records: `:binomial`, `:quasi-binomial`, `:gaussian`, `:gamma`, `:poisson`, `:quasi-poisson` and `:inverse-gaussian`. The `:quasi-` families estimate the dispersion, `:binomial` and `:poisson` fix it to `1.0`.

  Functions: `:nbinomial` (parameter `:nbinomial-theta`, default `1.0`).

  Each family has the canonical link as the default one: `:logit` for the binomial families, `:identity` for `:gaussian`, `:inverse` for `:gamma`, `:log` for the poisson families and `:nbinomial`, `:inversesq` for `:inverse-gaussian`. Use [[family-with-link]] to get a family merged with a link."}
  families {:binomial (->Family :logit (fn ^double [^double x] (m/* x (m/- 1.0 x)))
                                   binomial-initialize binomial-residual-deviance
                                   binomial-aic binomial-quantile-residuals 1.0)
               :quasi-binomial (->Family :logit (fn ^double [^double x] (m/* x (m/- 1.0 x)))
                                         binomial-initialize binomial-residual-deviance
                                         (constantly ##NaN) binomial-quantile-residuals :estimate)
               :gaussian (->Family :identity constantly-1
                                   default-initialize gaussian-residual-deviance gaussian-aic nil :estimate)
               :gamma (->Family :inverse m/sq
                                default-initialize gamma-residual-deviance gamma-aic
                                gamma-quantile-residuals :estimate)
               :poisson (->Family :log m/identity-double
                                  poisson-initialize poisson-residual-deviance poisson-aic
                                  poisson-quantile-residuals 1.0)
               :quasi-poisson (->Family :log m/identity-double
                                        poisson-initialize poisson-residual-deviance
                                        (constantly ##NaN) poisson-quantile-residuals :estimate)
               :inverse-gaussian (->Family :inversesq m/cb
                                           default-initialize inverse-gaussian-residual-deviance
                                           inverse-gaussian-aic inverse-gaussian-quantile-residuals
                                           :estimate)
               :nbinomial (fn [{:keys [^double nbinomial-theta]
                                :or {nbinomial-theta 1.0}}]
                            (->Family :log (fn ^double [^double x] (m/+ x (m// (m/* x x) nbinomial-theta)))
                                      nbinomial-initialize
                                      (->nbinomial-residual-deviance nbinomial-theta)
                                      (->nbinomial-aic nbinomial-theta)
                                      (->nbinomial-quantile-residuals nbinomial-theta)
                                      1.0))})

(defn family-with-link
  "Returns family with a link as single map."
  ([family] (family-with-link family nil))
  ([family params] (family-with-link family nil params))
  ([family link params]
   (let [fm (families family family)
         fm-data (if (fn? fm) (fm params) fm)
         link (or link (:default-link fm-data))
         lk (links link link)
         lk-data (if (fn? lk) (lk params) lk)]
     (assoc (merge fm-data lk-data) :family family :link link))))

(defn- estimate-dispersion
  ^double [residuals weights ^double df]
  (m// (v/sum (map (fn [^double r ^double w]
                     (m/* w r r))
                   residuals weights))
       df))

(defn- glm-cooks-distance
  [residuals hat ^double dp]
  (map (fn [^double r ^double h]
         (if (m/< h 1.0)
           (* (m/sq (m// r (m/- 1.0 h)))
              (m// h dp))
           ##NaN)) residuals hat))

(defn- glm-standardize-residuals
  [residuals hat ^double inv-sqrt-dispersion]
  (map (fn [^double r ^double h]
         (if (m/< h 1.0)
           (m/* r (m/sqrt (m// (m/- 1.0 h))) inv-sqrt-dispersion)
           ##NaN)) residuals hat))

(defn- glm-studentized-residuals
  [dresiduals presiduals sigmas hat family]
  (map (fn [^double dr ^double dp ^double s ^double h]
         (if (m/< h 1.0)
           (m// (m/* (m/signum dr)
                     (m/sqrt (m/+ (m/* dr dr)
                                  (m// (m/* h dp dp)
                                       (m/- 1.0 h))))) s)
           ##NaN)) dresiduals presiduals (if (#{:poisson :binomial} family)
                                           (repeat 1.0)
                                           sigmas) hat))

(defn- glm-extra-analysis
  [{:keys [weights hat residuals ^double dispersion
           ^long observations ^long nobs ^RealMatrix xtxinv family] :as model} ^RealMatrix xss]
  (let [initial-weights (:initial weights)
        weights (:weights weights)
        presiduals (:pearson residuals)
        dresiduals (:deviance residuals)
        rss (v/dot dresiduals dresiduals)
        inv-sqrt-dispersion (m// (m/sqrt dispersion))

        {^double df :residual} (:df model)
        p (mat/ncol xss)

        df- (m/dec df)
        sigmas (sigmas dresiduals hat rss (long df-))
        ^RealMatrix leverage-coeffs (leverage-coeffs xtxinv xss dresiduals
                                                     (map m/sqrt weights)
                                                     hat p observations)

        ;; rows with zero weight do not take part in the model, their diagnostics are NaN
        masked (fn [xs] (mask-zero-weights xs initial-weights))
        influence {:cooks-distance (masked (glm-cooks-distance presiduals hat (m/* dispersion p)))
                   :dffits (masked (dffits dresiduals hat sigmas))
                   :dfbetas (mask-zero-weights-rows (dfbetas leverage-coeffs xtxinv sigmas) initial-weights)
                   :covratio (masked (covratio dresiduals hat sigmas p nobs))}]
    {:residuals {:standardized {:pearson (masked (glm-standardize-residuals presiduals hat inv-sqrt-dispersion))
                                :deviance (masked (glm-standardize-residuals dresiduals hat inv-sqrt-dispersion))}
                 :studentized (masked (glm-studentized-residuals dresiduals presiduals sigmas hat family))}
     :leverage {:hat hat
                :sigmas (masked sigmas)
                :coefficients (mask-zero-weights-rows (map seq (mat/mat->array2d leverage-coeffs)) initial-weights)}
     :influence influence
     :influential (measures-influential-rows influence hat nobs p)
     :correlation (correlation xtxinv dispersion)}))

(defn- check-glm-options
  "Checks the options of `glm` which are not checked by [[check-fit-options]]; `m` is the number of observations."
  [{:keys [epsilon max-iters dispersion-estimator family link init-mu nbinomial-theta]
    :or {epsilon 1.0e-8 max-iters 25 family :gaussian}} m]
  (when-not (and (finite-number? epsilon) (m/>= (double epsilon) 0.0))
    (invalid-option :epsilon epsilon "a finite number >= 0"))
  (when-not (and (finite-number? max-iters) (m/>= (double max-iters) 1.0)
                 (m/== (double max-iters) (double (long max-iters))))
    (invalid-option :max-iters max-iters "an integer >= 1"))
  (when-not (or (nil? dispersion-estimator)
                (#{:pearson :mean-deviance} dispersion-estimator)
                (and (finite-number? dispersion-estimator) (m/pos? (double dispersion-estimator))))
    (invalid-option :dispersion-estimator dispersion-estimator ":pearson, :mean-deviance or a positive number"))
  (when-not (or (nil? nbinomial-theta) (and (finite-number? nbinomial-theta) (m/pos? (double nbinomial-theta))))
    (invalid-option :nbinomial-theta nbinomial-theta "a positive number"))
  (when-not (or (and (keyword? family) (contains? families family)) (instance? Family family) (fn? family))
    (invalid-option :family family (str "one of " (sort (keys families)) ", a Family record or a function")))
  (when-not (or (nil? link) (and (keyword? link) (contains? links link)) (instance? Link link) (fn? link))
    (invalid-option :link link (str "nil, one of " (sort (keys links)) ", a Link record or a function")))
  (when (some? init-mu)
    (when-not (and (or (sequential? init-mu) (.isArray (class init-mu)))
                   (m/== (count init-mu) (long m)) (every? finite-number? init-mu))
      (invalid-option :init-mu (if (sequential? init-mu) (count init-mu) init-mu) (str "a sequence of " m " finite numbers")))))

(defn- check-family-response
  "Checks if the response is in the domain of the family: proportions for binomial, non-negative for poisson
  and negative binomial, positive for gamma and inverse gaussian. A binomial response can be a pair
  `[successes failures]`. Custom families are not checked."
  [family ys]
  (let [domain (case family
                 (:binomial :quasi-binomial) [(fn [^double y] (m/<= 0.0 y 1.0)) "between 0 and 1"]
                 (:poisson :quasi-poisson :nbinomial) [(fn [^double y] (m/>= y 0.0)) ">= 0"]
                 (:gamma :inverse-gaussian) [(fn [^double y] (m/pos? y)) "> 0"]
                 :gaussian [(fn [^double _] true) "a finite number"]
                 nil)
        binomial? (#{:binomial :quasi-binomial} family)]
    (when domain
      (doseq [[i y] (map-indexed vector ys)]
        (if (sequential? y)
          (when-not (and binomial? (m/== (count y) 2) (every? #(and (finite-number? %) (m/>= (double %) 0.0)) y))
            (throw (ex-info "A pair of successes and failures should contain two non-negative numbers (binomial family only)."
                            {:option :ys :row i :value y :expected "[successes failures]"})))
          (when-not (and (finite-number? y) ((first domain) y))
            (throw (ex-info (str "Response is outside of the domain of the " (name family) " family.")
                            {:option :ys :row i :value y :expected (second domain)}))))))))

(defrecord GLMData [model transformer
                    ^RealMatrix xtxinv ys intercept? offset?
                    ^double intercept  beta coefficients ^long observations ^long nobs
                    residuals fitted weights offset ^double cv
                    names
                    deviance df ^double dispersion dispersions estimated-dispersion?
                    family link mean-fun link-fun iters quantile-residuals-fun
                    ^double q ^double chi2 ^double p-value
                    ll analysis decomposition]
  Object
  (toString [this] (->string this))
  IFn
  (invoke [m xs] (prot/predict m xs false))
  (invoke [m xs stderr?] (prot/predict m xs stderr?))
  (invoke [m] (arity-error m 0))
  (invoke [m _a _b _c] (arity-error m 3))
  (invoke [m _a _b _c _d] (arity-error m 4))
  (invoke [m _a _b _c _d _e] (arity-error m 5))
  (invoke [m _a _b _c _d _e _f] (arity-error m 6))
  (invoke [m _a _b _c _d _e _f _g] (arity-error m 7))
  (invoke [m _a _b _c _d _e _f _g _h] (arity-error m 8))
  (invoke [m _a _b _c _d _e _f _g _h _i] (arity-error m 9))
  (invoke [m _a _b _c _d _e _f _g _h _i _j] (arity-error m 10))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k] (arity-error m 11))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l] (arity-error m 12))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l _m] (arity-error m 13))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l _m _n] (arity-error m 14))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l _m _n _o] (arity-error m 15))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l _m _n _o _p] (arity-error m 16))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l _m _n _o _p _q] (arity-error m 17))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l _m _n _o _p _q _r] (arity-error m 18))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l _m _n _o _p _q _r _s] (arity-error m 19))
  (invoke [m _a _b _c _d _e _f _g _h _i _j _k _l _m _n _o _p _q _r _s _t] (arity-error m 20))
  (applyTo [m args] (apply-model m args))
  prot/PredictProto
  (predict [m xs] (prot/predict m xs false))
  (predict [_ xs stderr?]
    (let [[^double off xs] (with-offset xs offset?)
          xs (if transformer (transformer xs) xs)]
      (if stderr?
        (let [arr (double-array (if intercept? (cons 1.0 xs) xs))
              linear (m/+ off intercept (v/dot beta xs))
              fit (double (mean-fun linear))
              stderr (m/sqrt (m/* dispersion (v/dot arr (mat/mulv xtxinv arr))))
              scale (m/* stderr q)]
          {:fit fit
           :link linear
           :stderr stderr
           :confidence-interval [(mean-fun (m/- linear scale)) (mean-fun (m/+ linear scale))]})
        (mean-fun (m/+ off intercept (v/dot beta xs)))))))

;; https://bwlewis.github.io/GLM/
(defn glm
  "Fit a generalized linear model using IRLS method.

  Arguments:

  * `ys` - response vector
  * `xss` - terms of systematic component
  * optional parameters

  Parameters:

  * `:tol` - tolerance of the matrix decomposition: the relative threshold of small singular values (rank deficiency of the model matrix) and the threshold of the QR/RRQR decomposition, default: `1.0e-8`
  * `:epsilon` - tolerance for IRLS (stopping condition): relative change of the deviance, default: `1.0e-8`
  * `:max-iters` - maximum numbers of iterations, an integer >= 1, default: `25`
  * `:weights` - optional weights, finite and non-negative; an observation with the weight 0 stays in `:fitted` and `:residuals` but it does not take part in the model (it is not counted in `:nobs` and degrees of freedom)
  * `:offset` - optional offset
  * `:alpha` - significance level, default: `0.05`
  * `:intercept?` - should intercept term be included, default: `true`
  * `:init-mu` - initial response vector for IRLS (finite numbers, one for each observation)
  * `:simple?` - returns simplified result
  * `:dispersion-estimator` - `:pearson`, `:mean-deviance` or any positive number, replaces default one.
  * `:family` - family, default: `:gaussian`
  * `:link` - link
  * `:nbinomial-theta` - theta for `:nbinomial` family, default: `1.0`.
  * `:transformer` - an optional function which will be used to transform systematic component `xs` before fitting and prediction
  * `:names` - an optional vector of names to use when printing the model
  * `:decomposition` - which matrix decomposition use to find solution, `:cholesky` (default), `:rrqr` (rank revealing) or `:qr`
  * `:augmentation` and `augmentation-params` - regularization by data augmentation
      - `:ridge` - adds ridge regresion penalty (intercept is not penalized), default parameters `{:lambda 0.1}` 
      - `:diffs` - adds differences penalty, use with `b-spline-transformation` for smoothing, default parameters `{:lambda 1.0 :order 2}` 


  Family is one of the: `:gaussian` (default), `:binomial`, `:quasi-binomial`, `:poisson`, `:quasi-poisson`, `:gamma`, `:inverse-gaussian`, `:nbinomial`, custom `Family` record (see [[->family]]) or a function returning Family (accepting a map as an argument)

  Link is one of the: `:probit`, `:identity`, `:loglog`, `:sqrt`, `:inverse`, `:logit`, `:power` (parameter `:power-exponent`, default `1.0`; the exponent `0.0` is the `:log` link), `:nbinomial`, `:cauchit`, `:distribution` (parameter `:distribution`), `:cloglog`, `:inversesq`, `:log`, `:clog`, custom `Link` record (see [[->link]]) or a function returning Link (accepting a map as an argument)

  Notes:

  * SVD decomposition is used instead of more common QR
  * intercept term is added implicitely if `intercept?` is set to `true` (by default)
  * `:nbinomial` family requires `:nbinomial-theta` parameter
  * Each family has its own default (canonical) link.
  * For the `:binomial` families `ys` can be a sequence of pairs `[successes failures]`; the weights are multiplied by the number of trials.
  * Gaussian AIC uses the observations with a positive weight (R returns infinity when a weight is zero).

  The fit does not depend on the units of the terms nor on the scale of the weights (see [[lm]]).

  Throws `ex-info` when an option is invalid (`:option`, `:value` and `:expected` in the data), when the response is outside of the domain of the family (proportions for binomial, non-negative values for poisson and negative binomial, positive values for gamma and inverse gaussian, finite numbers otherwise), when the terms are not finite numbers or the rows have different lengths, when the lengths of `ys`, `xss`, `:weights`, `:offset` and `:init-mu` differ, when weights are negative or not finite, when the model has no terms, when the model matrix is rank deficient (`:rank`, `:columns`, `:tol` and `:singular-values` in the data), when the variance of the mean is invalid or the IRLS equations are singular, which usually means a diverging fit, eg. perfect separation (`:iteration` in the data) and when a family with an estimated dispersion has no residual degrees of freedom.

  Returned record implementes `IFn` protocol and contains:

  * `:model` - set to `:glm`
  * `:intercept?` - whether intercept term is included or not
  * `:xtxinv` - (X^T X)^-1
  * `:intercept` - intercept term value
  * `:beta` - vector of model coefficients (without intercept)
  * `:coefficients` - coefficient analysis, a list of maps containing `:estimate`, `:stderr`, `:t-value`, `:p-value` and `:confidence-interval`
  * `:weights` - weights, `:weights` (working) and `:initial`
  * `:residuals` - a map containing `:raw`, `:working`, `:pearsons` and `:deviance` residuals
  * `:fitted` - fitted values for xss
  * `:df` - degrees of freedom: `:residual`, `:null` and `:intercept`
  * `:observations` - number of observations (rows)
  * `:nobs` - number of observations with a positive weight, used in degrees of freedom and information criteria
  * `:deviance` - deviances: `:residual` and `:null`
  * `:dispersion` - default or calculated, used in a model
  * `:dispersions` - `:pearson` and `:mean-deviance`
  * `:family` - family used
  * `:link` - link used
  * `:link-fun` - link function, `g`
  * `:mean-fun` - mean function, `g^-1`
  * `:q` - (1-alpha/2) quantile of T or Normal distribution for residual degrees of freedom
  * `:chi2` and `:p-value` - Chi-squared statistic and respective p-value
  * `:ll` - a map containing log-likelihood and AIC/BIC (GLM and based on deviance, dev+2ED)
  * `:analysis` - leverage, residual and influence analysis - a delay
  * `:iters` and `:converged?` - number of iterations and convergence indicator (`true` when the change of the deviance is smaller than `:epsilon`, also at the last allowed iteration)
  * `:decomposition` - decomposition used

  Analysis, delay containing a map:

  * `:residuals` - `:standardized` and `:studentized` residuals (pearsons and deviance)
  * `:leverage` - `:hat`, `:sigmas` and leveraged `:coefficients` (leave-one-out)
  * `:influence` - `:cooks-distance`, `:dffits`, `:dfbetas` and `:covratio`
  * `:influential` - list of influential observations (ids) for influence measures
  * `:correlation` - correlation matrix of estimated parameters"
  ([ys xss] (glm ys xss nil))
  ([ys xss {:keys [^long max-iters ^double tol ^double epsilon family link weights ^double alpha offset
                   dispersion-estimator intercept? init-mu simple? transformer names decomposition
                   augmentation augmentation-params]
            :or {max-iters 25 tol 1.0e-8 epsilon 1.0e-8 family :gaussian alpha 0.05 intercept? true
                 simple? false decomposition :cholesky}
            :as params}]

   (check-fit-options {:tol tol :alpha alpha :decomposition decomposition
                       :augmentation augmentation :augmentation-params augmentation-params})
   (check-data ys xss weights offset intercept?)
   (check-glm-options params (count ys))
   (check-family-response family ys)
   (let [offset? (boolean offset)
         xss (check-terms xss transformer)

         m (long (count ys))
         
         [ys ^RealMatrix xss weights offset] (augment-data ys xss weights offset intercept?
                                                           augmentation augmentation-params)


         {link-fun :g link-mean :mean
          link-derivative :derivative link-variance :variance
          initialize :initialize
          residual-deviance :residual-deviance
          dispersion :dispersion
          link :link family :family
          aic-fun :aic
          :as family+link} (family-with-link family link params)
         
         new-m (mat/nrow xss)
         n (mat/ncol xss)

         [ys start-t weights optional-data] (initialize (seq ys) (seq weights))

         ;; observations with zero weight do not take part in the model
         nobs (count-positive (take m weights))
         estimated-dispersion? (boolean (or (= :estimate dispersion) dispersion-estimator))
         _ (when (or (m/zero? nobs)
                     (and estimated-dispersion? (not (number? dispersion-estimator))
                          (not augmentation) (m/<= nobs n)))
             (throw (ex-info "No residual degrees of freedom." {:observations nobs :parameters n :df (m/- nobs n)})))

         ;; the rank and the system do not depend on the scale of terms and weights
         ^doubles scales (column-scales xss)
         ^RealMatrix xss-scaled (scale-columns xss scales)
         [^SingularValueDecomposition S singular-values] (svd xss-scaled tol)
         uts (.getUT S)
         ut (.getU S)

         ^doubles ys (m/seq->double-array ys)
         ^doubles weights (m/seq->double-array weights)
         ^doubles offset (m/seq->double-array offset)
         ^RealVector rvoffset (v/vec->RealVector offset)

         init-t (double-array (map link-fun (or init-mu start-t)))

         ^doubles buff-g (double-array new-m)
         ^doubles buff-z (double-array new-m)
         ^doubles buff-W (double-array (repeat new-m 1.0))
         ^doubles buff-Wz (double-array new-m)         
         
         [^RealVector result t ^long iters converged?]
         (loop [iter (long 1)
                ^doubles t init-t
                dev ##Inf]

           ;; iterate over t
           (dotimes [idx m]
             (let [eta (Array/aget t idx)
                   g (double (link-mean eta))
                   v (double (link-variance g))
                   g' (double (link-derivative eta))]
               (when (or (m/zero? v) (m/nan? v) (m/nan? g))
                 (throw (ex-info "Invalid variance of mean, the fit diverged (perfect separation or an invalid starting point?)."
                                 {:eta eta :mean g :dmean g' :variance v :coeff idx :iteration iter})))

               (let [off (Array/aget offset idx)
                     z (m/+ (m/- eta off) (m// (m/- (Array/aget ys idx) g) g'))
                     w (m/* (Array/aget weights idx) (m// (m/* g' g') v))]
                 (Array/aset buff-g idx g)
                 (Array/aset buff-z idx z)
                 (Array/aset buff-W idx w)
                 (Array/aset buff-Wz idx (m/* w z)))))

           (let [wscale (mean-positive buff-W)
                 wfactor (m// wscale)
                 ^RealVector new-t (->> (try (-> (mat/mulm uts (mat/mulm (DiagonalMatrix. ^doubles (v/mult buff-W wfactor)) ut))
                                                 ^DecompositionSolver (solver decomposition tol)
                                                 (.solve (.operate uts (v/vec->RealVector (v/mult buff-Wz wfactor)))))
                                             (catch SingularMatrixException e
                                               (throw (ex-info "The IRLS equations are singular, check the response (eg. all values equal or all successes) and perfect separation."
                                                               {:iteration iter} e)))
                                             (catch NonPositiveDefiniteMatrixException e
                                               (throw (ex-info "The IRLS equations are not positive definite, check the response and perfect separation, or try the :qr or :rrqr decomposition."
                                                               {:iteration iter} e))))
                                        (mat/mulv ut))
                 new-dev (v/sum (residual-deviance ys buff-g weights))]
             
             (if (m/< (m// (m/abs (m/- new-dev dev))
                           (m/+ 0.1 (m/abs new-dev))) epsilon)
               [(->> singular-values
                     (v/emult (mat/mulv uts new-t))
                     (mat/mulv (.getV S)))
                (v/add (v/vec->array new-t) offset)
                iter
                true]
               (if (m/== iter max-iters)
                 [(->> singular-values
                       (v/emult (mat/mulv uts new-t))
                       (mat/mulv (.getV S)))
                  (v/add (v/vec->array new-t) offset)
                  iter
                  false]

                 (recur (m/inc iter)
                        (v/add (v/vec->array new-t) offset)
                        new-dev)))))

         result-array (unscale-coefficients result scales)
         result (v/vec->RealVector result-array)
         
         intercept (if intercept? (Array/aget result-array 0) 0.0)
         beta (vec (if intercept? (rest result-array) result-array))

         xtx-1 (xtxinv xss-scaled scales buff-W tol)
         hat (hat-matrix xtx-1 xss buff-W m)

         ed (v/sum hat)
         df (if augmentation (m/- nobs ed) (m/- nobs n))
         
         intercept-df (long (if intercept? 1 0))
         null-df (m/- nobs intercept-df)

         ;; remove augmentation
         fitted (map link-mean (take m (v/vec->seq (v/add rvoffset (mat/mulv xss result)))))
         
         ;; remove augmentation
         ys (take m ys)
         weights (take m weights)
         offset (take m offset)
         t (take m t)
         buff-W (take m buff-W)

         ;; fix deviance after removing augmentation
         dev (v/sum (residual-deviance ys fitted weights))
         
         raw-residuals (map m/- ys fitted)
         residuals (v/ediv raw-residuals (map link-derivative t))

         loocv (map (fn [^double r ^double h] (m// r (m/- 1.0 h))) residuals hat)

         ;; loocv, cross validation, zero weight observations are not a part of the model
         cv (m/sqrt (m// (v/sum (v/sq (positive-only loocv weights))) nobs))

         _ (when (and estimated-dispersion? (not (number? dispersion-estimator)) (not (m/pos? (double df))))
             (throw (ex-info "No residual degrees of freedom." {:observations nobs :parameters n :df df})))
         pearson-dispersion (estimate-dispersion residuals buff-W df)
         mean-deviance (m// dev (double df))

         dispersion-value (cond
                            (not estimated-dispersion?) dispersion
                            (number? dispersion-estimator) (double dispersion-estimator)
                            (= dispersion-estimator :mean-deviance) mean-deviance
                            :else pearson-dispersion)

         distr (if estimated-dispersion?
                 (r/distribution :t {:degrees-of-freedom df})
                 (r/distribution :normal))
         q (double (r/icdf distr (m/- 1.0 (m/* alpha 0.5))))

         stderrs (standard-errors xtx-1 dispersion-value)
         namev (coefficients-names names intercept? (count beta))
         
         model {:model :glm
                :transformer transformer
                :xtxinv xtx-1
                :intercept? intercept?
                :offset? offset?
                :estimated-dispersion? estimated-dispersion?
                :intercept intercept
                :beta beta
                :hat hat
                :observations m
                :nobs nobs
                :residuals {:working residuals
                            :raw raw-residuals
                            :loocv loocv}
                :deviance {:residual dev}
                :cv cv
                :df {:residual df :null null-df :intercept intercept-df}
                :effective-dimension ed
                :dispersion dispersion-value
                :family family
                :link link
                :names namev
                :mean-fun link-mean
                :link-fun link-fun
                :iters iters
                :converged? converged?
                :q q :chi2 ##NaN
                :p-value ##NaN
                :decomposition decomposition}]

     (if simple?
       (map->GLMData model)

       (let [pearson-residuals (map (fn [^double y ^double mu ^double w]
                                      (m/* (m/- y mu)
                                           (m/sqrt (m// w (double (link-variance mu))))))
                                    ys fitted weights)
             deviance-residuals (->> (residual-deviance ys fitted weights)
                                     (map m/safe-sqrt)
                                     (map (fn [^double y ^double mu ^double r]
                                            (if (m/<= y mu) (m/- r) r)) ys fitted))

             chi2 (v/sum (map m/sq pearson-residuals))

             null-deviance (if (and offset? intercept?)
                             ;; if offset if set and intercept is used
                             ;; null-deviance is a deviance from model fitted with intercept only
                             ;; ys and weights are already initialized (eg. successes divided by trials),
                             ;; the full model fit is the starting point, the same as in R
                             (get-in (glm ys :intercept
                                          (assoc params
                                                 :simple? true
                                                 :transformer nil
                                                 :weights weights
                                                 :augmentation nil
                                                 :augmentation-params nil
                                                 :init-mu fitted)) [:deviance :residual])
                             ;; offset is not set
                             (v/sum (residual-deviance ys (if intercept?
                                                            (repeat m (stats/mean ys weights))
                                                            (map link-mean offset)) weights)))

             nn (double (if augmentation ed n))
             aic (if aic-fun (double (aic-fun ys fitted weights dev nobs nn optional-data)) ##Inf)
             [log-likelihood bic] (let [p (if (#{:gaussian :inverse-gaussian
                                                 :gamma :nbinomial} family) (m/inc nn) nn)
                                        ll (m/- p (m// aic 2.0))]
                                    [ll (m/- (m/* p (m/log nobs)) (m/* 2.0 ll))])

             model (-> model
                       (assoc :ys (v/vec->seq ys)
                              :coefficients (coefficients-analysis result-array stderrs distr q)
                              :offset (v/vec->seq offset)
                              :weights {:weights (v/vec->seq buff-W)
                                        :initial (v/vec->seq weights)}
                              :fitted fitted
                              :dispersions {:pearson pearson-dispersion
                                            :mean-deviance mean-deviance}
                              :quantile-residuals-fun (:quantile-residuals-fun family+link)
                              :chi2 chi2
                              :p-value (if (m/zero? (double df))
                                         ##Inf
                                         (stats/p-value (r/distribution :chi-squared {:degrees-of-freedom df})
                                                        chi2 :right))
                              :effective-dimension ed
                              :ll {:log-likelihood log-likelihood
                                   :aic aic
                                   :bic bic
                                   ;; p-spline AIC, dev+2*ED
                                   :aic-dev (m/+ dev (m/* 2.0 ed))
                                   :bic-dev (m/+ dev (m/* (m/log nobs) ed))})
                       (assoc-in [:residuals :pearson] pearson-residuals)
                       (assoc-in [:residuals :deviance] deviance-residuals)
                       (assoc-in [:deviance :null] null-deviance))

             analysis (delay (glm-extra-analysis model (if augmentation
                                                         (.getSubMatrix xss 0 (m/dec m) 0 (m/dec n))
                                                         xss)))]
         
         (map->GLMData (assoc model :analysis analysis)))))))

;; fit nbinomial theta

(defn- nbinomial-theta-score [ys mus weights]
  (fn ^double [^double theta]
    (let [dg (special/digamma theta)
          lt (m/log theta)]
      (v/sum (map (fn [^double y ^double mu ^double w]
                    (let [y+th (m/+ y theta)
                          mu+th (m/+ mu theta)]
                      (m/* w (m/- (m/+ (special/digamma y+th) lt 1.0)
                                  dg (m/log mu+th) (m// y+th mu+th)))))
                  ys mus weights)))))

(defn- nbinomial-theta-score' [ys mus weights]
  (fn ^double [^double theta]
    (let [tg (special/trigamma theta)
          rt (m// theta)]
      (v/sum (map (fn [^double y ^double mu ^double w]
                    (let [y+th (m/+ y theta)
                          mu+th (m/+ mu theta)]
                      (m/* w (m/- (m/+ tg (m// 2.0 mu+th))
                                  (special/trigamma y+th) rt (m// y+th (m/sq mu+th))))))
                  ys mus weights)))))

(defn- nbinomial-theta-init
  ^double [ys mus weights ^long n]
  (m// n (v/sum (map (fn [^double y ^double mu ^double w]
                       (m/* w (m/sq (m/dec (m// y mu))))) ys mus weights))))

(defn- optimize-theta
  ^double [{:keys [ys fitted weights]} ^long max-iters ^double epsilon]
  (let [w (:initial weights)
        n (v/sum w)
        init (nbinomial-theta-init ys fitted w (long n))
        f    (nbinomial-theta-score ys fitted w)
        f'   (nbinomial-theta-score' ys fitted w)]
    (loop [i (long 0)
           theta init
           delta 1.0]
      (if (or (m/< delta epsilon) (m/== i max-iters))
        theta
        (let [ndelta (m// (double (f theta))
                          (double (f' theta)))]
          (recur (m/inc i) (m/+ theta ndelta) ndelta))))))

(defn- nbinomial-theta-stderr
  "Standard error of theta, from the observed information at the estimate (`nbinomial-theta-score'` is the information, the derivative of the score with the opposite sign)."
  ^double [{:keys [ys fitted weights]} ^double theta]
  (let [info (nbinomial-theta-score' ys fitted (:initial weights))]
    (m// (m/sqrt (double (info theta))))))

(defn glm-nbinomial
  "Fits theta for negative binomial glm in iterative process.

  Returns fitted model with the keys `:nbinomial-theta` and `:nbinomial-theta-stderr` (the standard error of theta from the observed information, `SE.theta` of `glm.nb` in R).

  Arguments and parameters are the same as for `glm`.

  Additional parameters:

  * `:nbinomial-theta` - initial theta used as a starting point for optimization."
  ([ys xss] (glm-nbinomial ys xss nil))
  ([ys xss {:keys [nbinomial-theta ^long max-iters ^double epsilon]
            :or {max-iters 25 epsilon 1.0e-8}
            :as params}]
   (when-not (or (nil? nbinomial-theta) (and (finite-number? nbinomial-theta) (m/pos? (double nbinomial-theta))))
     (invalid-option :nbinomial-theta nbinomial-theta "a positive number"))
   (let [fit0 (if nbinomial-theta
                (glm ys xss (assoc params :family :nbinomial))
                (glm ys xss (assoc params :family :poisson)))
         theta0 (optimize-theta fit0 max-iters epsilon)
         init-mu (:fitted fit0)]
     (loop [iter (long 0)
            model fit0
            theta theta0
            ll ##Inf
            check false]
       (if (or check (m/>= iter max-iters))
         (assoc model
                :nbinomial-theta theta
                :nbinomial-theta-stderr (nbinomial-theta-stderr model theta))
         (let [nmodel (glm ys xss (assoc params
                                         :family :nbinomial
                                         :nbinomial-theta theta
                                         :init-mu init-mu))
               ntheta (optimize-theta nmodel max-iters epsilon)
               nll (double (get-in nmodel [:ll :log-likelihood]))
               ncheck (m/+ (m/abs (m// (m/- nll ll)))
                           (m/abs (m/- ntheta theta)))]
           (recur (m/inc iter) nmodel ntheta nll (m/< ncheck epsilon))))))))

(defn quantile-residuals
  "Quantile residuals of a GLM.

  The residuals are standard normal when the model is correct. For discrete families (`:binomial`, `:poisson`, `:nbinomial`) they are randomized: a value is drawn in the interval of the cumulative distribution (with the global random generator, so the result differs between calls). For other families they are deterministic. A model without a family specific function gets standardized deviance residuals.

  Parameters:

  - `model` - a model returned by [[glm]]

  Returns a sequence of residuals.

  See also [[analysis]]."
  [{:keys [quantile-residuals-fun residuals dispersion] :as model}]
  (cond
    quantile-residuals-fun (quantile-residuals-fun model)
    (nil? (:deviance residuals)) (throw (ex-info "Residuals are not available for a model fitted with :simple? true."
                                                  {:model (:model model)}))
    :else (v/div (:deviance residuals) (m/sqrt dispersion))))

(defn analysis
  "Influence analysis, leverage, standardized and studentized residuals, correlation.

  Parameters:

  - `model` - a model returned by [[lm]] or [[glm]]

  Returns the map described in the documentation of [[lm]] and [[glm]] (it is calculated once, on the first call). Entries for observations with a zero weight are `##NaN` and these observations are never listed in `:influential`.

  See also [[quantile-residuals]]."
  [model]
  (if-let [a (:analysis model)]
    (deref a)
    (throw (ex-info "Analysis is not available for a model fitted with :simple? true." {:model (:model model)}))))

(defn dose
  "Predict Lethal/Effective dose for given `p` (default: p=0.5, median).

  * intercept-id - id of intercept, default: 0
  * coeff-id is the coefficient used for calculating dose, default: 1

  The model can be fitted with `:simple? true`. Throws `ex-info` when an index is not an index of a coefficient. For `p` of `0` or `1` the link is infinite: the dose and its `:stderr` are infinite.

  `:stderr` is scaled by the model's `:dispersion` (1.0 for fixed-dispersion
  families such as binomial/poisson, so a no-op there; the estimated value
  for gaussian/gamma/inverse-gaussian/quasi-* families, where the parameter
  covariance is `dispersion * xtxinv`, not `xtxinv` alone)."
  ([glm-model] (dose glm-model 0.5))
  ([glm-model ^double p] (dose glm-model p 1))
  ([glm-model ^double p ^long coeff-id] (dose glm-model p 0 coeff-id))
  ([{:keys [link-fun ^RealMatrix xtxinv intercept? ^double intercept beta ^double dispersion]}
    ^double p ^long intercept-id ^long coeff-id]
   (let [estimates (vec (if intercept? (cons intercept beta) beta))
         n (count estimates)]
     (doseq [[option ^long id] [[:intercept-id intercept-id] [:coeff-id coeff-id]]]
       (when (and (m/> n 1) (not (and (m/>= id 0) (m/< id n))))
         (invalid-option option id (str "an index of a coefficient, 0.." (dec n)))))
     (if (m/> n 1)
       (let [e (double (link-fun p))
             submatrix-idxs (int-array [intercept-id coeff-id])
             b0 (double (nth estimates intercept-id))
             b1 (double (nth estimates coeff-id))
             xp (m// (m/- e b0) b1)
             pd (mat/rows->RealMatrix [[(m// -1.0 b1) (m/- (m// xp b1))]])
             var-entry (mat/entry (->> (.getSubMatrix xtxinv submatrix-idxs submatrix-idxs)
                                       (mat/mulm pd)
                                       (mat/mulmt pd)) 0 0)
             ;; the link of p = 0 or p = 1 is infinite: so is the dose and its error
             err (if (m/valid-double? e) (m/sqrt (m/* dispersion var-entry)) ##Inf)]
         {:dose xp
          :p p
          :stderr err})
       {:dose ##Inf
        :p p
        :stderr ##Inf}))))

(defn- residuals-quantiles
  "Summary statistics used for printing residuals with model"
  [residuals]
  [(zipmap [:min :q1 :median :q3 :max]
           (map #(m/approx % 6) [(apply min residuals)
                                 (stats/quantile residuals 0.25)
                                 (stats/median residuals)
                                 (stats/quantile residuals 0.75)
                                 (apply max residuals)]))])

(defn- coefficients-table
  [nm coefficient]
  (let [tz (if (:t-value coefficient) :t-value :z-value)]
    (zipmap [:name :estimate :stderr tz :p-value :confidence-interval]
            [nm
             (m/approx (:estimate coefficient) 6)
             (m/approx (:stderr coefficient) 6)
             (m/approx (tz coefficient) 6)
             (m/approx (:p-value coefficient) 6)
             (mapv #(m/approx % 6) (:confidence-interval coefficient))])))

;; Model printing methods
(defmulti ^{:doc "Returns a text summary of the fitted model (`LMData` or `GLMData`): residuals, coefficient table and fit statistics. It is used by `print-method` and `toString` of the models."} ->string class)

(defmethod ->string LMData [{:keys [model names coefficients residuals
                                    r-squared adjusted-r-squared f-statistic
                                    p-value df sigma ll]}]
  (with-out-str
    (println (str (when (= model :wls) "Weighted ") "Residuals:")) ;; Print residuals
    (-> (residuals-quantiles (:weighted residuals))
        pprint/print-table)
    (println)
    (println "Coefficients:")
    (->> coefficients
         (map-indexed #(coefficients-table (nth names %1) %2))
         pprint/print-table)
    (println)
    (println (str "F-statistic: " f-statistic " on degrees of freedom: " df))
    (println (str "p-value: " p-value))
    (println)
    (println (str "R2: " r-squared))
    (println (str "Adjusted R2: " adjusted-r-squared))
    (println (str "Residual standard error: " sigma " on " (:residual df) " degrees of freedom"))
    (println (str "AIC: " (:aic ll)))))

(defmethod ->string GLMData [{:keys [residuals coefficients family link dispersion iters
                                     converged? ll deviance df names intercept? intercept beta]}]
  (with-out-str
    (when-let [dres (:deviance residuals)]
      (println "Deviance Residuals:")
      (-> (residuals-quantiles dres)
          pprint/print-table)
      (println))
    (println "Coefficients:")
    (if coefficients
      (->> coefficients
           (map-indexed #(coefficients-table (nth names %1) %2))
           pprint/print-table)
      ;; a model fitted with :simple? true has no coefficients analysis
      (->> (if intercept? (cons intercept beta) beta)
           (map (fn [nm ^double estimate] {:name nm :estimate (m/approx estimate 6)}) names)
           pprint/print-table))
    (println)
    (println "Dispersion: " dispersion)
    (println)
    (println "Null deviance: " (:null deviance) " on " (:null df) " degrees of freedom")
    (println "Residual deviance: " (:residual deviance) " on " (:residual df) " degrees of freedom")
    (println (str "AIC: " (:aic ll)))
    (println)
    (println (str "Family: " family ", Link: " link))
    (println (str "Fisher Scoring Iterations (converged? " (if converged? "yes" "no") "): " iters))))

(defmethod print-method LMData [data ^java.io.Writer w]
(.write w (str (->string data))))
(defmethod print-method GLMData [data  ^java.io.Writer w]
(.write w (str (->string data))))


;; common transformers

(defn polynomial-transformer
  "Creates a transformer which expands a number into orthogonal polynomial terms, as `poly(x, degree)` in R.

  Parameters:

  - `xs` - training values of the variable, the polynomials are orthogonal on them
  - `degree` - the degree of polynomial

  Returns a function which accepts a one element data point `[x]` (a model `:transformer`) and returns a vector of `degree` terms. The terms of the training values are orthonormal (columns of `poly(x, degree)`), the same function can be evaluated outside the training values.

  See also [[lm]], [[glm]], [[trigonometric-transformer]], [[b-spline-transformer]]."
  [xs ^long degree]
  (let [different (count (distinct xs))]
    (when-let [[x] (first-invalid finite-number? xs)]
      (invalid-option :xs x "finite numbers"))
    (when-not (m/>= different 2)
      (invalid-option :xs xs "at least two different values"))
    (when-not (and (m/>= degree 1) (m/< degree different))
      (invalid-option :degree degree (str "an integer between 1 and " (dec different) " (the number of different values of xs minus one)"))))
  (let [degree+ (m/inc degree)
        r (range degree+)
        m (stats/mean xs)
        xs (v/shift xs (m/- m))
        qrX (-> (mat/cols->RealMatrix (map (partial v/pow xs) r))
                (mat/qr-decomposition))
        R (mat/decomposition-component qrX :R)
        Q (mat/decomposition-component qrX :Q)
        Z (mat/cols->RealMatrix (map (fn [^long i] (v/mult (mat/col Q i) (mat/entry R i i))) r))
        Z2 (mat/emulm Z Z)
        norm2 (map v/sum (mat/cols Z2))
        alpha (butlast (map (fn [col ^double n2]
                              (m/+ m (m// (v/sum (v/emult (v/vec->seq col) xs)) n2))) (mat/cols Z2) norm2))
        alpha1 (double (first alpha))
        anormf (map (fn [^double a ^double n] (Vec2. a n))
                    (rest alpha)
                    (map (fn [[^double n1 ^double n2]] (m// n2 n1)) (partition 2 1 (butlast norm2))))
        normsqrt (mapv m/sqrt (rest norm2))]
    (fn [[^double x]]
      (let [f (m/- x alpha1)]
        (-> ((reduce (fn [[buff ^double zi- ^double zi] ^Vec2 an]
                       (let [zi+ (m/- (m/* (m/- x (.x an)) zi)
                                      (m/* (.y an) zi-))]
                         [(conj buff zi+) zi zi+])) [[f] 1.0 f] anormf) 0)
            (v/ediv normsqrt))))))

(defn trigonometric-transformer
  "Creates a transformer which expands a number into trigonometric (Fourier) terms.

  Parameters:

  - `period` - the period of the variable
  - `degree` - the number of harmonics

  Returns a function which accepts a one element data point `[x]` (a model `:transformer`) and returns a sequence of `2*degree` terms: `sin(2 pi k x / period)` and `cos(2 pi k x / period)` for `k=1..degree`.

  See also [[lm]], [[glm]], [[polynomial-transformer]], [[b-spline-transformer]]."
  [^double period ^double degree]
  (when-not (and (m/valid-double? period) (m/pos? period))
    (invalid-option :period period "a finite number > 0"))
  (when-not (and (m/valid-double? degree) (m/>= degree 1.0) (m/== degree (m/rint degree)))
    (invalid-option :degree degree "an integer >= 1"))
  (let [r (range 1 (m/inc degree))]
    (fn [[^double x]]
      (let [xx (m/* m/TWO_PI (m// x period))]
        (mapcat (fn [^long k]
                  (let [kxx (m/* k xx)]
                    [(m/sin kxx) (m/cos kxx)])) r)))))

(defn b-spline-transformer
  "Creates a transformer which expands a number into B-spline basis functions on equally spaced knots, as `splineDesign` in R.

  Parameters:

  - `xs` - training values of the variable, their minimum and maximum are `xl` and `xr`
  - `xl`, `xr` - the range of the variable (the arity without `xs`)
  - `nseg` - number of segments between `xl` and `xr`
  - `degree` - the degree of B-splines, at least `1` (`3` is cubic)

  Returns a function which accepts a one element data point `[x]` (a model `:transformer`) and returns a sequence of `nseg + degree` basis functions evaluated at `x`. The functions sum up to `1.0` for `x` in `[xl, xr]`.

  Use the `:diffs` penalty of the `:augmentation` option of [[lm]] and [[glm]] for smoothing (P-splines).

  See also [[polynomial-transformer]], [[trigonometric-transformer]], [[add-penalty]]."
  ([xs ^long nseg ^long degree]
   (let [[xl xr] (stats/extent xs)]
     (b-spline-transformer xl xr nseg degree)))
  ([^double xl ^double xr ^long nseg ^long degree]
   (when-not (and (m/valid-double? xl) (m/valid-double? xr) (m/< xl xr))
     (invalid-option :xl xl (str "finite numbers xl < xr, got xr " xr)))
   (when-not (m/>= nseg 1)
     (invalid-option :nseg nseg "an integer >= 1"))
   (when-not (m/>= degree 1)
     (invalid-option :degree degree "an integer >= 1"))
   (let [degree+ (m/long-inc degree)
         dx (m// (m/- xr xl) nseg)
         degdx (m/* dx degree)
         knots (double-array (range (m/- xl degdx) (m/+ xr degdx (m/* 0.1 dx)) dx))
         n (alength knots)
         b1 (if (m/even? degree+) 1 -1)
         denom (m// (m/* b1 (special/gamma degree+) (m/pow dx degree)))
         D (mapv v/vec->array (mat/rows (mat/muls (mat/differences (mat/eye n true) degree+) denom)))
         nb (m/long-sub n degree+)
         sk (double-array (take nb (drop degree+ knots)))]
     (fn [[^double x]]
       ;; unrolled tpow, 25% speedup when compared with simple `map`
       (if (m/nan? x)
         (repeat nb ##NaN)
         (let [p (double-array n)]
           (loop [i (long 0)]
             (when (m/< i n)
               (let [k (Array/aget knots i)
                     diff (m/- x k)]
                 (when (m/pos? diff)
                   (Array/aset p i (m/fpow diff degree))
                   (recur (m/inc i))))))         
           (map (fn [d ^double k]
                  (if (m/< x k)
                    (v/dot p d)
                    0.0)) D sk)))))))

;;;;;;;


;;;;;;


;; isotonic

;; https://www.jstatsoft.org/article/view/v102c01

(defn- ir-decode-order
  "Decodes the `order` of isotonic regression to the code of `Monotone`: 0 ascending, 1 descending, 2 and 3 the
  `:non-decreasing` and `:non-increasing` variants. `nil` is ascending."
  [order]
  (case order
    nil 0
    :asc 0
    :increasing 0
    :desc 1
    :decreasing 1
    :non-decreasing 2
    :non-increasing 3
    (throw (ex-info (str "Unknown order of isotonic regression: " (pr-str order)
                         ". Accepted values: nil, :asc, :increasing, :desc, :decreasing, :non-decreasing, :non-increasing.")
                    {:order order
                     :accepted [:asc :increasing :desc :decreasing :non-decreasing :non-increasing]}))))

(defn- check-isotonic
  "Checks the data of the isotonic regression: finite response, positive finite weights and finite, sorted `xs`
  (when given), all of the same length as `ys`."
  [ys ws xs]
  (let [n (count ys)]
    (when-let [[y] (first-invalid finite-number? ys)]
      (invalid-option :ys y "finite numbers"))
    (when (some? ws)
      (when-not (m/== (count ws) n)
        (throw (ex-info "Length of weights should be the same as the length of the response."
                        {:option :ws :value (count ws) :expected n})))
      (when-let [[w] (first-invalid #(and (finite-number? %) (m/pos? (double %))) ws)]
        (invalid-option :ws w "finite numbers > 0")))
    (when (some? xs)
      (when-not (m/== (count xs) n)
        (throw (ex-info "Length of xs should be the same as the length of the response."
                        {:option :xs :value (count xs) :expected n})))
      (when-let [[x] (first-invalid finite-number? xs)]
        (invalid-option :xs x "finite numbers"))
      (when-not (every? true? (map <= (map double xs) (map double (rest xs))))
        (invalid-option :xs "unsorted" "values in the ascending order")))))

(defn pava
  "Isotonic regression, pool-adjacent-violators algorithm with up-and-down-blocks variant.

  Isotonic regression minimizes the (weighted) L2 loss function with a constraint that result should be monotonic (ascending or descending).

  Arguments:
  
  - `ys` - response variable data
  - `ws` - weights (optional)
  - `order` - `:asc` or `:increasing` (default), `:desc` or `:decreasing`, `:non-decreasing` and `:non-increasing`

  Returns monotonic predicted values, an empty sequence for an empty `ys`.

  Throws `ex-info` for a not finite response, for weights which are not finite and positive, for the weights of a different length than `ys` and for an unknown order."
  ([ys] (pava ys :increasing))
  ([ys order] (pava ys nil order))
  ([ys ws order]
   (check-isotonic ys ws nil)
   (if (empty? ys)
     (do (ir-decode-order order) ())
     (let [^doubles y (double-array ys)
           n (alength y)
           w (double-array (or ws (repeat n 1.0)))
           r (int-array (m/inc n))]
       (->> (Monotone/pava_step1 y w r (ir-decode-order order))
            (Monotone/pava_step2 y r))
       (seq y)))))

;; https://arxiv.org/abs/1701.05964

(defn cir
  "Centered Isotonic Regression.

  Returns shrinked [`xs`,`ys`] pair.

  Arguments:
  - `xs` - regressor variable
  - `ys` - response variable
  - `ws` - weights (optional)
  - `order` - `:asc` or `:increasing` (default), `:desc` or `:decreasing`, `:non-decreasing` and `:non-increasing`.

  Returns a pair of empty sequences for an empty `ys`.

  Throws `ex-info` for a not finite response, for `xs` which are not finite or not sorted ascending, for weights which are not finite and positive, for the data of different lengths and for an unknown order."
  ([ys] (cir (range (count ys)) ys))
  ([xs ys] (cir xs ys :increasing))
  ([xs ys order] (cir xs ys nil order))
  ([xs ys ws order]
   (check-isotonic ys ws xs)
   (if (empty? ys)
     (do (ir-decode-order order) [() ()])
     (let [^doubles x (double-array xs)
           ^doubles y (double-array ys)
           n (alength y)
           ^doubles w (double-array (or ws (repeat n 1.0)))
           ^ints r (int-array (m/inc n))
           b (Monotone/pava_step1 x y w r (ir-decode-order order))]
       (loop [b (long b)
              f n
              buff '()]
         (if (m/neg? b)
           [(map first buff)
            (map second buff)]
           (let [t (aget r b)
                 yd (Array/aget y b)
                 nbuff (cond-> buff
                         (and (m/== f n) (m/> (m/- f t) 1)) (conj [(Array/aget x (m/dec f)) yd])
                         true (conj [(Array/aget x b) yd])
                         (and (m/zero? t) (m/> f 1)) (conj [(first xs) yd]))]
             (recur (m/dec b) t nbuff))))))))

