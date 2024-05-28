(ns fastmath.ml.regression
  "OLS, WLS and GLM regression models with analysis."
  (:require [fastmath.core :as m]
            [fastmath.stats :as stats]
            [fastmath.vector :as v]
            [fastmath.protocols :as prot]
            [fastmath.random :as r]
            [fastmath.matrix :as mat])
  (:import [org.apache.commons.math3.linear SingularValueDecomposition DiagonalMatrix
            CholeskyDecomposition RealMatrix RealVector]
           [fastmath.java Array]
           [clojure.lang IFn]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- with-offset
  [xs offset?]
  (if offset?
    (let [xs (v/vec->seq xs)]
      [(first xs) (rest xs)])
    [0.0 (v/vec->seq xs)]))

(defn predict
  "Predict from the given model and data point.

  If `stderr?` is true, standard error and confidence interval is added.
  If model is fitted with offset, first element of data point should contain provided offset.

  Expected data point:
  
  * `[x1,x2,...,xn]` - when model was trained without offset
  * `[offset,x1,x2,...,xn]` - when offset was used for training
  * `[]` or `nil` - when model was trained with intercept only
  * `[offset]` - when model was trained with intercept and offset"
  ([model xs] (prot/predict model xs))
  ([model xs stderr?] (prot/predict model xs stderr?)))

(defrecord LMData [model intercept? offset?
                   ^RealMatrix xtxinv
                   ^double intercept beta coefficients
                   offset weights residuals fitted df ^long observations
                   ^double r-squared ^double adjusted-r-squared ^double sigma2 ^double sigma
                   ^double tss ^double rss ^double regss ^double msreg
                   ^double qt
                   ^double f-statistic ^double p-value
                   ll analysis]
  IFn
  (invoke [m xs] (prot/predict m xs false))
  (invoke [m xs stderr?] (prot/predict m xs stderr?))
  prot/PredictProto
  (predict [m xs] (prot/predict m xs false))
  (predict [_ xs stderr?]
    (let [^double off xs] (with-offset xs offset?)
         (if stderr?
           (let [[^double off xs] (with-offset xs offset?)
                 arr (double-array (if intercept? (conj xs 1.0) xs))
                 fit (double (m/+ off intercept (v/dot beta xs)))
                 stderr (m/sqrt (m/* sigma2 (v/dot arr (mat/mulv xtxinv arr))))
                 scale (m/* stderr qt)]
             {:fit fit
              :stderr stderr
              :confidence-interval [(m/- fit scale) (m/+ fit scale)]})
           (m/+ off intercept (v/dot beta xs))))))

;;

(defn- small-singular-values-count
  ^long [^doubles sv ^double epsilon]
  (let [f (Array/aget sv 0)]
    (count
     (filter (fn [^double v] (m/< (m// v f) epsilon)) sv))))

(defn- svd
  [ys xss tol intercept?]
  (let [xss (if (= xss :intercept) (repeat (count ys) '()) xss)
        xss (if intercept?
              (map (fn [xs] (conj (v/vec->seq xs) 1.0)) xss)
              (map v/vec->seq xss))
        
        ^RealMatrix xss (mat/mat xss)
        ^SingularValueDecomposition S (SingularValueDecomposition. xss)
        singular-values (.getSingularValues S)
        k (small-singular-values-count singular-values tol)]

    (if (m/pos? k)
      (throw (ex-info "Near rank-deficient model matrix" {:data singular-values}))
      [xss S (-> singular-values v/reciprocal v/vec->RealVector)])))

(defn- xtxinv
  ^RealMatrix [^RealMatrix xss ^DiagonalMatrix dweights ^double tol]
  (-> (mat/tmulm xss (mat/mulm dweights xss))
      (CholeskyDecomposition. tol 1.0e-16)
      (.getSolver)
      (.getInverse)))

(defn- standard-errors
  [^RealMatrix xtxinv ^double dispersion]
  (-> xtxinv
      (mat/diag)
      (v/mult dispersion)
      (v/sqrt)
      (v/vec->seq)))

(defn- normality-analysis
  [{:keys [weighted raw]}]
  (let [skew (stats/skewness weighted :g1)
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

(def ^:private ^:const ME100 (m/- 1.0 (m/* 100.0 m/MACHINE-EPSILON)))

(defn- hat-matrix
  [^RealMatrix xtxinv ^RealMatrix xss weights]
  (let [sqrt-weights (m/seq->double-array (map m/sqrt weights))
        ^RealMatrix whx (mat/mulm (DiagonalMatrix. sqrt-weights) xss)]
    (->> (mat/mulmt xtxinv whx)
         (mat/mulm whx)
         (mat/diag)
         (v/vec->seq)
         (map (fn [^double h] (if (m/>= h ME100) 1.0 h))))))

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
  [^RealMatrix laverage-coeffs ^RealMatrix xtx-1 sigmas]
  (let [d (v/vec->array (mat/diag xtx-1))]
    (mapv (fn [^RealVector v ^double s]
            (v/ediv (v/vec->seq v) (v/mult sigmas (m/sqrt s)))) (mat/rows laverage-coeffs) d)))

(defn- covratio
  [residuals hat sigmas ^long p]
  (let [n (count residuals)
        n-p (m/- n p)
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
                 :dfbetas (mapv influential-rows (:dfbetas influence) (repeat (fn [^double x] (m/> (m/abs x) 1.0))))}]
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

(defn- laverage-coeffs
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
  [{:keys [weights residuals ^double rss ^long observations
           ^double sigma ^double sigma2 ^RealMatrix xtxinv] :as model} ^RealMatrix xss]
  (let [rresiduals (:raw residuals)
        wresiduals (:weighted residuals)
        {^double df :residual} (:df model)
        p (mat/ncol xss)
        hat (hat-matrix xtxinv xss weights)        
        
        df- (m/dec df)
        
        sigmas (sigmas wresiduals hat rss df-)
        ^RealMatrix laverage-coeffs (laverage-coeffs xtxinv xss rresiduals weights hat p observations)

        influence {:cooks-distance (lm-cooks-distance wresiduals hat sigma p)
                   :dffits (dffits wresiduals hat sigmas)
                   :dfbetas (dfbetas laverage-coeffs xtxinv sigmas)
                   :covratio (covratio wresiduals hat sigmas p)}]
    {:normality (normality-analysis residuals)
     :residuals {:standardized (lm-transform-residuals wresiduals hat (repeat sigma))
                 :studentized (lm-transform-residuals wresiduals hat sigmas)}
     :laverage {:hat hat
                :sigmas sigmas
                :coefficients (mapv seq (mat/mat->array2d laverage-coeffs))}
     :influence influence
     :influential (measures-influential-rows influence hat observations p)
     :correlation (correlation xtxinv sigma2)}))

(defn lm
  "Fit a linear model using ordinary (OLS) or weighted (WLS) least squares.

  Arguments:

  * `ys` - response vector
  * `xss` - terms of systematic component
  * optional parameters

  Parameters:

  * `:tol` - tolerance for matrix decomposition (SVD and Cholesky), default: `1.0e-8`
  * `:weights` - optional weights for WLS
  * `:offset` - optional offset
  * `:alpha` - significance level, default: `0.05`
  * `:intercept?` - should intercept term be included, default: `true`

  Notes:
  
  * SVD decomposition is used instead of more common QR
  * intercept term is added implicitely if `intercept?` is set to `true` (by default)
  * Two variants of AIC/BIC are calculated, one based on log-likelihood, second on RSS/n

  Returned record implementes `IFn` protocol and contains:
  
  * `:model` - `:ols` or `:wls`
  * `:intercept?` - whether intercept term is included or not
  * `:xtxinv` - (X^T X)^-1
  * `:intercept` - intercept term value
  * `:beta` - vector of model coefficients (without intercept)
  * `:coefficients` - coefficient analysis, a list of maps containing `:estimate`, `:stderr`, `:t-value`, `:p-value` and `:confidence-interval`
  * `:weights` - initial weights
  * `:residuals` - a map containing `:raw` and `:weighted` residuals
  * `:fitted` - fitted values for xss
  * `:df` - degrees of freedom: `:residual`, `:model` and `:intercept`
  * `:observations` - number of observations
  * `:r-squared` and `:adjusted-r-squared`
  * `:sigma` and `:sigma2` - deviance and variance
  * `:msreg` - regression mean squared 
  * `:rss`, `:regss`, `:tss` - residual, regression and total sum of squares
  * `:qt` - (1-alpha/2) quantile of T distribution for residual degrees of freedom
  * `:f-statistic` and `:p-value` - F statistic and respective p-value
  * `:ll` - a map containing log-likelihood and AIC/BIC in two variants: based on log-likelihood and RSS
  * `:analysis` - laverage, residual and influence analysis - a delay

  Analysis, delay containing a map:
  
  * `:residuals` - `:standardized` and `:studentized` weighted residuals
  * `:laverage` - `:hat`, `:sigmas` and laveraged `:coefficients` (leave-one-out)
  * `:influence` - `:cooks-distance`, `:dffits`, `:dfbetas` and `:covratio`
  * `:influential` - list of influential observations (ids) for influence measures
  * `:correlation` - correlation matrix of estimated parameters
  * `:normality` - residuals normality tests: `:skewness`, `:kurtosis`, `:durbin-watson` (for raw and weighted), `:jarque-berra` and `:omnibus` (normality)"
  ([ys xss] (lm ys xss nil))
  ([ys xss {:keys [^double tol weights ^double alpha intercept? offset]
            :or {tol 1.0e-8 alpha 0.05 intercept? true}}]

   
   (let [[^RealMatrix xss ^SingularValueDecomposition S singular-values] (svd ys xss tol intercept?)
         uts (.getUT S)
         ut (.getU S)

         m (mat/nrow xss)
         n (mat/ncol xss)

         offset? (boolean offset)
         ^doubles offset (m/seq->double-array (or offset (repeat m 0.0)))

         ^doubles ys-orig (m/seq->double-array ys)
         ^doubles ys (v/sub ys-orig offset)
         
         ;; maybe separate wls and ols...
         weights? (sequential? weights)           
         weights (vec (or weights (repeat m 1.0)))
         ^doubles daweights (m/seq->double-array weights)
         ^DiagonalMatrix dweights (DiagonalMatrix. daweights)
         ^RealVector rvwys (v/vec->RealVector (v/emult daweights ys)) 
         
         ^RealVector result (let [new-t (->> (-> (mat/mulm uts (mat/mulm dweights ut))
                                                 (CholeskyDecomposition. tol 1.0e-16)
                                                 (.getSolver)
                                                 (.solve ^RealVector (mat/mulv uts rvwys)))
                                             (mat/mulv ut))]

                              (->> singular-values
                                   (v/emult (mat/mulv uts new-t))
                                   (mat/mulv (.getV S))))

         fitted (v/vec->seq (mat/mulv xss result))
         result-array (v/vec->array result)
         
         intercept (if intercept? (Array/aget result-array 0) 0.0)
         beta (vec (if intercept? (rest result-array) result-array))
         
         df (m/- m n)
         model-df (count beta)
         intercept-df (long (if intercept? 1 0))

         raw-residuals (mapv m/- ys fitted)
         wresiduals (if weights?
                      (mapv (fn [^double r ^double w] (m/* r (m/sqrt w))) raw-residuals weights)
                      raw-residuals)
         
         rss (v/dot wresiduals wresiduals)
         sigma2 (m// rss df)
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
         
         regss (m/- tss rss)
         rsquared (m/- 1.0 (m// rss tss))
         adjusted-rsquared (m/- 1.0 (m/* (m// (m/- m intercept-df) (double df)) (m/- 1.0 rsquared)))
         f-statistic (m// (m// regss model-df) sigma2)
         p-value (stats/p-value (r/distribution :f {:numerator-degrees-of-freedom model-df
                                                    :denominator-degrees-of-freedom df})
                                f-statistic :one-sided-greater)
         tdistr (r/distribution :t {:degrees-of-freedom df})
         qt (double (r/icdf tdistr (m/- 1.0 (m/* alpha 0.5))))

         ;; log-likelihood based measures
         ll (let [nobs2 (m// m 2.0)
                  log-likelihood (m/+ (m/- (m/- (m/* nobs2 (m/log rss)))
                                           (m/* nobs2 (m/inc (m/log (m// m/PI nobs2)))))
                                      (m/* 0.5 (v/sum (v/log weights))))
                  mlogsigma (m/* m (m/log (m// rss m)))]
              {:log-likelihood log-likelihood
               :aic (m/+ (m/* -2.0 log-likelihood) (* 2.0 (m/inc n)))
               :bic (m/+ (m/* -2.0 log-likelihood) (* (m/log m) (inc n)))
               :aic-rss (m/+ mlogsigma (m/* 2.0 n))
               :bic-rss (m/+ mlogsigma (m/* (m/log m) n))})

         xtx-1 (xtxinv xss dweights tol)
         stderrs (standard-errors xtx-1 sigma2)
         
         model {:model (if weights? :wls :ols)
                :intercept? intercept?
                :offset? offset?
                :xtxinv xtx-1
                :intercept intercept
                :beta beta
                :offset offset
                :coefficients (coefficients-analysis result-array stderrs tdistr qt)
                :weights weights
                :residuals {:weighted wresiduals :raw raw-residuals}
                :fitted (v/add fitted (seq offset))
                :df {:residual df :model model-df :intercept intercept-df}
                :observations m
                :r-squared rsquared
                :adjusted-r-squared adjusted-rsquared
                :sigma2 sigma2
                :sigma sigma
                :tss tss :rss rss :regss regss :msreg (m// regss model-df)
                :qt qt
                :f-statistic f-statistic :p-value p-value
                :ll ll}

         analysis (delay (lm-extra-analysis model xss))]

     (map->LMData (assoc model :analysis analysis)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; GLM https://bwlewis.github.io/GLM/
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- constantly-1 ^double [^double _] 1.0)
(defn- default-initialize [ys weights ^long obs] [ys ys weights (repeat obs 1.0)])

(defn- binomial-initialize
  [ys weights ^long obs]
  (if (sequential? (first ys)) ;; check if we have pairs of successes and failures
    (let [n (map (fn [[^double s ^double f]] (m/+ s f)) ys)
          y (v/ediv (map first ys) n)]
      [y (map (fn [^double y ^double n]
                (m// (m/+ 0.5 (m/* n y))
                     (m/inc n))) y n)
       (v/emult weights n) n])
    [ys
     (map (fn [^double y ^double w]
            (m// (m/+ 0.5 (m/* w y))
                 (m/inc w))) ys weights)
     weights (repeat obs 1.0)]))

(defn- binomial-aic
  [ys ns fitted weights _deviance _observations rank]
  (let [ms (if (some (fn [^long n] (m/> n 1)) ns) ns weights)]
    (m/+ (m/* -2.0 (v/sum (map (fn [^double y ^double m ^double mu ^double wt]
                                 (if (m/not-pos? m)
                                   0.0
                                   (let [distr (r/distribution :binomial {:trials m
                                                                          :p mu})]
                                     (m/* (m// wt m)
                                          (r/lpdf distr (m/round (m/* m y)))))))
                               ys ms fitted weights)))
         (m/* 2.0 (int rank)))))

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
  [_ys _ns _fitted weights deviance observations rank]
  (let [nobs (long observations)]
    (m/+ (m/* nobs (m/inc (m/log (m/* m/TWO_PI (m// (double deviance) nobs)))))
         2.0 (m/- (v/sum (v/log weights)))
         (m/* 2.0 (int rank)))))

(defn- gamma-residual-deviance
  [ys mus weights]
  (map (fn [^double y ^double mu ^double w]
         (m/* -2.0 w (m/- (m/log (if (m/zero? y) 1.0 (m// y mu)))
                          (m// (m/- y mu) mu))))
       ys mus weights))

(defn- poisson-initalize
  [ys weights ^long obs]
  [ys (v/shift ys 0.1) weights (repeat obs 1.0)])

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
  [ys _ns fitted weights _deviance _observations rank]
  (m/+ (m/* -2.0 (v/sum (map (fn [^double y ^double mu ^double w]
                               (m/* w (r/lpdf (r/distribution :poisson {:p mu}) y))) ys fitted weights)))
       (m/* 2.0 (int rank))))

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

(defn- nbinomial-initialize
  [ys weights ^long obs]
  [ys (map (fn ^double [^double y]
             (if (m/zero? y) m/SIXTH y))) weights (repeat obs 1.0)])

(defn- nbinomial-residual-deviance
  [^double theta]
  (fn [ys mus weights]
    (map (fn [^double y ^double mu ^double w]
           (let [y+theta (m/+ y theta)]
             (m/* 2.0 w (m/- (m/* y (m/log (m// (m/max y 1.0) mu)))
                             (m/* y+theta
                                  (m/log (m// y+theta (m/+ mu theta))))))))
         ys mus weights)))

(defn- estimate-dispersion
  ^double [residuals weights ^long df]
  (m// (v/sum (map (fn [^double r ^double w]
                     (m/* w r r))
                   residuals weights))
       df))

(defrecord Link [g mean derivative])

(def links {:logit (->Link m/logit
                         m/sigmoid
                         (fn ^double [^double x] (let [e (m/exp x)] (m// e (m/sq (m/inc e))))))
          :probit (->Link (fn ^double [^double x] (m/* m/SQRT2 (m/inv-erf (m/dec (m/* 2.0 x)))))
                          (fn ^double [^double x] (-> (m// x m/SQRT2)
                                                     (m/-)
                                                     (m/erfc)
                                                     (m/* 0.5)))
                          (fn ^double [^double x] (m/max m/MACHINE-EPSILON (m/* m/INV_SQRT2PI
                                                                               (m/exp (m/* -0.5 x x))))))
          :cauchit (->Link (fn ^double [^double x] (m/tan (m/* m/PI (m/- x 0.5))))
                           (fn ^double [^double x] (-> x
                                                      (m/atan)
                                                      (m// m/PI)
                                                      (m/+ 0.5)))
                           (fn ^double [^double x] (m/* m/INV_PI (m// 1.0 (m/inc (m/* x x))))))
          :cloglog (->Link m/cloglog
                           (fn ^double [^double x] (m/constrain (m/cexpexp x)
                                                               m/MACHINE-EPSILON
                                                               0.9999999999999999))
                           (fn ^double [^double x] (m/max m/MACHINE-EPSILON (m/exp (m/- x (m/exp x))))))
          :loglog (->Link m/loglog
                          m/expexp
                          (fn ^double [^double x] (m/max m/MACHINE-EPSILON
                                                        (m/exp (m/- (m/- (m/exp (m/- x))) x)))))
          :identity (->Link m/identity-double
                            m/identity-double
                            constantly-1)
          :log (->Link m/log
                       (fn ^double [^double x] (m/max m/MACHINE-EPSILON (m/exp x)))
                       (fn ^double [^double x] (m/max m/MACHINE-EPSILON (m/exp x))))
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
          :power (fn [{:keys [^double exponent]
                      :or {exponent 1.0}}]
                   (let [rexponent (m// exponent)]
                     (->Link (fn ^double [^double x] (m/pow x exponent))
                             (fn ^double [^double x] (m/pow x rexponent))
                             (fn ^double [^double x] (m/* rexponent (m/pow x (m/* rexponent
                                                                                 (m/- 1.0 exponent))))))))
          :distribution (fn [{:keys [distribution]
                             :or {distribution (r/distribution :normal)}}]
                          (->Link (partial r/icdf distribution)
                                  (partial r/cdf distribution)
                                  (partial r/pdf distribution)))})

(defrecord Family [default-link variance initialize residual-deviance aic quantile-residuals dispersion])

(def families {:binomial (->Family :logit (fn ^double [^double x] (m/* x (m/- 1.0 x)))
                                 binomial-initialize binomial-residual-deviance
                                 binomial-aic binomial-quantile-residuals 1.0)
             :quasi-binomial (->Family :logit (fn ^double [^double x] (m/* x (m/- 1.0 x)))
                                       binomial-initialize binomial-residual-deviance
                                       (constantly ##NaN) nil :estimate)
             :gaussian (->Family :identity constantly-1
                                 default-initialize gaussian-residual-deviance gaussian-aic nil :estimate)
             :gamma (->Family :inverse m/sq
                              default-initialize gamma-residual-deviance nil nil :estimate)
             :poisson (->Family :log m/identity-double 
                                poisson-initalize poisson-residual-deviance poisson-aic
                                poisson-quantile-residuals 1.0)
             :quasi-poisson (->Family :log m/identity-double 
                                      poisson-initalize poisson-residual-deviance
                                      (constantly ##NaN) nil :estimate)
             :inverse-gaussian (->Family :inversesq m/cb
                                         default-initialize inverse-gaussian-residual-deviance nil nil
                                         :estimate)
             :nbinomial (fn [{:keys [^double nbinomial-theta]
                             :or {nbinomial-theta 1.0}}]
                          (->Family :log (fn ^double [^double x] (m/+ x (m// (m/* x x) nbinomial-theta)))
                                    nbinomial-initialize (nbinomial-residual-deviance nbinomial-theta) nil nil
                                    :estimate))})

(defn family-with-link
  "Returns family with optional link"
  ([family link params]
   (let [fm (families family family)
         fm-data (if (fn? fm) (fm params) fm)
         link (or link (:default-link fm))
         lk (links link link)
         lk-data (if (fn? lk) (lk params) lk)]
     (assoc (merge fm-data lk-data) :family family :link link))))

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
  [dresiduals presiduals sigmas hat dispersion]
  (let [d (double dispersion)]
    (map (fn [^double dr ^double dp ^double s ^double h]
           (if (m/< h 1.0)
             (m// (m/* (m/signum dr)
                       (m/sqrt (m/+ (m/* dr dr)
                                    (m// (m/* h dp dp)
                                         (m/- 1.0 h))))) s)
             ##NaN)) dresiduals presiduals (if (m/one? d)
                                             (repeat 1.0)
                                             sigmas) hat)))

(defn- glm-extra-analysis
  [{:keys [weights residuals ^double dispersion
           ^long observations ^RealMatrix xtxinv] :as model} ^RealMatrix xss]
  (let [weights (:weights weights)        
        presiduals (:pearson residuals)
        dresiduals (:deviance residuals)
        rss (v/dot dresiduals dresiduals)
        inv-sqrt-dispersion (m// (m/sqrt dispersion))
        
        {^double df :residual} (:df model)
        p (mat/ncol xss)
        hat (hat-matrix xtxinv xss weights)                
        
        df- (m/dec df)
        sigmas (sigmas dresiduals hat rss df-)
        ^RealMatrix laverage-coeffs (laverage-coeffs xtxinv xss dresiduals
                                                     (map m/sqrt weights)
                                                     hat p observations)

        influence {:cooks-distance (glm-cooks-distance presiduals hat (m/* dispersion p))
                   :dffits (dffits dresiduals hat sigmas)
                   :dfbetas (dfbetas laverage-coeffs xtxinv sigmas)
                   :covratio (covratio dresiduals hat sigmas p)}]
    {:residuals {:standardized {:pearson (glm-standardize-residuals presiduals hat inv-sqrt-dispersion)
                                :deviance (glm-standardize-residuals dresiduals hat inv-sqrt-dispersion)}
                 :studentized (glm-studentized-residuals dresiduals presiduals sigmas hat dispersion)}
     :laverage {:hat hat
                :sigmas sigmas
                :coefficients (mapv seq (mat/mat->array2d laverage-coeffs))}
     :influence influence
     :influential (measures-influential-rows influence hat observations p)
     :correlation (correlation xtxinv dispersion)}))

(defrecord GLMData [model
                    ^RealMatrix xtxinv ys intercept? offset?
                    ^double intercept  beta coefficients ^long observations
                    residuals fitted weights offset
                    deviance df ^double dispersion dispersions estimated-dispersion?
                    family link mean-fun link-fun iters
                    ^double q ^double chi2 ^double p-value
                    ll analysis]
  IFn
  (invoke [m xs] (prot/predict m xs false))
  (invoke [m xs stderr?] (prot/predict m xs stderr?))
  prot/PredictProto
  (predict [m xs] (prot/predict m xs false))
  (predict [_ xs stderr?]
    (let [[^double off xs] (with-offset xs offset?)]
      (if stderr?
        (let [[^double off xs] (with-offset xs offset?)
              arr (double-array (conj xs 1.0))
              linear (m/+ off intercept (v/dot beta xs))
              fit (double (mean-fun linear))
              stderr (m/sqrt (m/* dispersion (v/dot arr (mat/mulv xtxinv arr))))
              scale (m/* stderr q)]
          {:fit fit
           :link linear
           :stderr stderr
           :confidence-interval [(mean-fun (m/- linear scale)) (mean-fun (m/+ linear scale))]})
        (mean-fun (m/+ off intercept (v/dot beta xs)))))))

(defn glm
  ([ys xss] (glm ys xss nil))
  ([ys xss {:keys [^long max-iters ^double tol ^double epsilon family link weights ^double alpha offset
                   dispersion-estimator intercept? init-mu simple?]
            :or {max-iters 25 tol 1.0e-8 epsilon 1.0e-8 family :gaussian alpha 0.05 intercept? true
                 simple? false}
            :as params}]

   (let [[^RealMatrix xss ^SingularValueDecomposition S singular-values] (svd ys xss epsilon intercept?)
         uts (.getUT S)
         ut (.getU S)

         {link-fun :g link-mean :mean
          link-derivative :derivative link-variance :variance
          initialize :initialize
          residual-deviance :residual-deviance
          dispersion :dispersion
          link :link family :family
          aic-fun :aic} (family-with-link family link params)

         m (mat/nrow xss)
         n (mat/ncol xss)

         weights (or weights (repeat m 1.0))

         [ys start-t ^doubles weights ns] (initialize (seq ys) (seq weights) m)

         ^doubles ys (m/seq->double-array ys)
         ^doubles weights (m/seq->double-array weights)

         offset? (boolean offset)
         ^doubles offset (m/seq->double-array (or offset (repeat m 0.0)))
         ^RealVector rvoffset (v/vec->RealVector offset)

         init-t (double-array (map link-fun (or init-mu start-t)))
         
         ^doubles buff-g (double-array m)
         ^doubles buff-z (double-array m)
         ^doubles buff-W (double-array m)
         ^doubles buff-Wz (double-array m)

         [^RealVector result t ^double dev ^long iters]
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
                 (throw (ex-info "Invalid variance of mean." {:mean g :variance v :coeff idx})))
               
               (let [off (Array/aget offset idx)
                     z (m/+ (m/- eta off) (m// (m/- (Array/aget ys idx) g) g'))
                     w (m/* (Array/aget weights idx) (m// (m/* g' g') v))]
                 (Array/aset buff-g idx g)
                 (Array/aset buff-z idx z)
                 (Array/aset buff-W idx w)
                 (Array/aset buff-Wz idx (m/* w z)))))           

           (let [^RealVector new-t (->> (-> (mat/mulm uts (mat/mulm (DiagonalMatrix. buff-W) ut))
                                            #_(QRDecomposition. tol)
                                            (CholeskyDecomposition. tol 1.0e-16)
                                            (.getSolver)
                                            (.solve (.operate uts (v/vec->RealVector buff-Wz))))
                                        (mat/mulv ut))
                 new-dev (v/sum (residual-deviance ys buff-g weights))]
             
             (if (or (m/< (m// (m/abs (m/- new-dev dev))
                               (m/+ 0.1 (m/abs new-dev))) epsilon)
                     (m/== iter max-iters))
               [(->> singular-values
                     (v/emult (mat/mulv uts new-t))
                     (mat/mulv (.getV S)))
                (v/add (v/vec->array new-t) offset)
                new-dev
                iter]
               
               (recur (m/inc iter)
                      (v/add (v/vec->array new-t) offset)
                      new-dev))))

         result-array (v/vec->array result)

         intercept (if intercept? (Array/aget result-array 0) 0.0)
         beta (vec (if intercept? (rest result-array) result-array))

         df (m/- m n)
         intercept-df (long (if intercept? 1 0))
         null-df (m/- m intercept-df)

         fitted (map link-mean (v/vec->seq (v/add rvoffset (mat/mulv xss result))))

         raw-residuals (map m/- ys fitted)
         residuals (v/ediv raw-residuals (map link-derivative t))
         
         estimated-dispersion? (or (= :estimate dispersion) dispersion-estimator)
         pearson-dispersion (estimate-dispersion residuals buff-W df)
         mean-deviance (m// dev df)

         dispersion-value (cond
                            (not estimated-dispersion?) dispersion
                            (number? dispersion-estimator) (double dispersion-estimator)
                            (= dispersion-estimator :mean-deviance) mean-deviance
                            :else pearson-dispersion)
         
         distr (if estimated-dispersion?
                 (r/distribution :t {:degrees-of-freedom df})
                 (r/distribution :normal))
         q (double (r/icdf distr (m/- 1.0 (m/* (double alpha) 0.5))))
         
         xtx-1 (xtxinv xss (DiagonalMatrix. buff-W) tol)
         stderrs (standard-errors xtx-1 dispersion-value)

         model {:model :glm
                :xtxinv xtx-1
                :intercept? intercept?
                :offset? offset?
                :estimated-dispersion? estimated-dispersion?
                :intercept intercept
                :beta beta
                :observations m
                :residuals {:working residuals
                            :raw raw-residuals}
                :deviance {:residual dev}
                :df {:residual df :null null-df :intercept intercept-df}
                :dispersion dispersion-value
                :family family
                :link link
                :mean-fun link-mean
                :link-fun link-fun                
                :iters iters
                :converged? (m/< iters max-iters)
                :q q :chi2 ##NaN
                :p-value ##NaN}]

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
                             (get-in (glm ys :intercept
                                          (assoc params
                                                 :simple? true
                                                 :init-t fitted)) [:deviance :residual])
                             ;; offset is not set
                             (v/sum (residual-deviance ys (if intercept?
                                                            (repeat m (stats/mean ys weights))
                                                            (map link-mean offset)) weights)))

             aic (if aic-fun (double (aic-fun ys ns fitted weights dev m n)) ##Inf)
             [log-likelihood bic] (let [p (if (#{:gaussian :inverse-gaussian :gamma} family) (m/inc n) n)
                                        ll (m/- p (m// aic 2.0))]
                                    [ll (m/- (m/* p (m/log m)) (m/* 2.0 ll))])
             
             model (-> model
                       (assoc :ys (v/vec->seq ys)
                              :coefficients (coefficients-analysis result-array stderrs distr q)
                              :offset (v/vec->seq offset)
                              :weights {:weights (v/vec->seq buff-W)
                                        :initial (v/vec->seq weights)}
                              :fitted fitted
                              :dispersions {:pearson pearson-dispersion
                                            :mean-deviance mean-deviance}
                              :chi2 chi2
                              :p-value (if (m/zero? df)
                                         ##Inf
                                         (stats/p-value (r/distribution :chi-squared {:degrees-of-freedom df})
                                                        chi2 :right))
                              :ll {:log-likelihood log-likelihood
                                   :aic aic
                                   :bic bic})
                       (assoc-in [:residuals :pearson] pearson-residuals)
                       (assoc-in [:residuals :deviance] deviance-residuals)
                       (assoc-in [:deviance :null] null-deviance))
             
             analysis (delay (glm-extra-analysis model xss))]

         (map->GLMData (assoc model :analysis analysis)))))))

(defn quantile-residuals
  "Quantile residuals for a model, possibly randomized."
  [{:keys [family residuals ^double dispersion] :as model}]
  (if-let [qres (get-in families [family :quantile-residuals])]
    (qres model)
    (v/div (:deviance residuals) (m/sqrt dispersion))))

(defn dose
  "Predict Lethal/Effective dose for given `p` (default: p=0.5, median).
  
  * intercept-id - id of intercept, default: 0
  * coeff-id is the coefficient used for calculating dose, default: 1"
  ([glm-model] (dose glm-model 0.5))
  ([glm-model ^double p] (dose glm-model p 1))
  ([glm-model ^double p ^long coeff-id] (dose glm-model p 0 coeff-id))
  ([{:keys [link-fun ^RealMatrix xtxinv coefficients]} ^double p ^long intercept-id ^long coeff-id]
   (if (m/> (count coefficients) 1)
     (let [e (double (link-fun p))
           submatrix-idxs (int-array [intercept-id coeff-id])
           coeffs (map :estimate coefficients)
           b0 (double (nth coeffs intercept-id))
           b1 (double (nth coeffs coeff-id))
           xp (m// (m/- e b0) b1)
           pd (mat/rows->RealMatrix [[(m// -1.0 b1) (m/- (m// xp b1))]])
           err (m/sqrt (mat/entry (->> (.getSubMatrix xtxinv submatrix-idxs submatrix-idxs)
                                       (mat/mulm pd)
                                       (mat/mulmt pd)) 0 0))]
       {:dose xp
        :p p
        :stderr err})
     {:dose ##Inf
      :p p
      :stderr ##Inf})))
