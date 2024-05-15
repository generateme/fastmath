(ns fastmath.ml.regression
  (:require [fastmath.core :as m]
            [fastmath.stats :as stats]
            [fastmath.vector :as v]
            [fastmath.protocols :as prot]
            [fastmath.random :as r]
            [fastmath.matrix :as mat])
  (:import [org.apache.commons.math3.stat.regression OLSMultipleLinearRegression]
           [org.apache.commons.math3.linear SingularValueDecomposition DiagonalMatrix
            CholeskyDecomposition RealMatrix RealVector]
           [fastmath.java Array]
           [clojure.lang IFn]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn predict [model xs] (prot/predict model xs))

(defrecord LMData [^OLSMultipleLinearRegression model
                   ^double intercept beta coefficients
                   residual fitted df ^long observations
                   ^double r-squared ^double adjusted-r-squared
                   ^double stderr ^double tss ^double rss ^double ess ^double mse
                   ^double f-statistic ^double p-value
                   ^double log-likelihood ^double aic ^double bic]
  IFn
  (invoke [_ xs] (m/+ (v/dot beta xs) intercept))
  prot/PredictProto
  (predict [_ xs] (m/+ (v/dot beta xs) intercept)))

(defn- lm-residual-analysis
  [residuals ^double rss ^long df]
  (let [skew (stats/skewness residuals :g1)
        kurt (stats/kurtosis residuals :kurt)
        normality (stats/normality-test residuals skew (- kurt 3.0) nil)]
    {:ss rss :df df
     :residuals residuals
     :ms (m// rss df)
     :skewness skew
     :kurtosis kurt
     :durbin-watson (stats/durbin-watson residuals)
     :jarque-berra (dissoc (stats/jarque-bera-test residuals skew kurt nil) :skewness :kurtosis)
     :normality (dissoc normality :skewness :kurtosis)}))

(defn- lm-analysis
  [beta stderrs ^long df]
  (let [tdistr (r/distribution :t {:degrees-of-freedom df})]
    (map (fn [^double b ^double err]
           (let [t-value (m// b err)]
             {:estimate b :stderr err :t-value t-value
              :p-value (stats/p-value tdistr t-value :both)})) beta stderrs)))

(defn lm
  ([ys xss] (lm ys xss nil))
  ([ys xss {:keys [^double qr-threshold]
            :or {qr-threshold 0.0}}]
   (let [model (OLSMultipleLinearRegression. qr-threshold)
         xss (if (number? (first xss)) (map vector xss) xss)]
     (.newSampleData model (m/seq->double-array ys) (m/seq->double-double-array xss))
     (let [^doubles beta (.estimateRegressionParameters model)
           ^doubles stderrs (.estimateRegressionParametersStandardErrors model)
           n (count xss)           
           intercept (Array/aget beta 0)
           coefficients (vec (rest beta))
           residuals (.estimateResiduals model)
           fitted (map (fn [xs] (m/+ (v/dot xs coefficients) intercept)) xss)
           df (m/- (count xss) (alength beta))
           model-df (count coefficients)
           stderr (.estimateRegressionStandardError model)
           tss (.calculateTotalSumOfSquares model)
           rss (.calculateResidualSumOfSquares model)
           ess (m/- tss rss)
           r2 (.calculateRSquared model)
           f-statistic (m// (m// ess model-df) (m/sq stderr))
           log-likelihood (m/* (m// n -2) (m/+ m/LOG_TWO_PI (m/log (m// rss n)) 1.0))]
       (->LMData model intercept coefficients (lm-analysis beta stderrs df)
                 (lm-residual-analysis (seq residuals) rss df)
                 fitted {:residuals df :model model-df} n
                 r2 (.calculateAdjustedRSquared model)
                 stderr tss rss ess (m// ess model-df)
                 f-statistic (stats/p-value (r/distribution :f {:numerator-degrees-of-freedom model-df
                                                                :denominator-degrees-of-freedom df})
                                            f-statistic :one-sided-greater)
                 log-likelihood (m/* -2.0 (m/- log-likelihood (m/inc model-df)))
                 (m/+ (m/* -2.0 log-likelihood)
                      (m/* (m/log n) (m/inc model-df))))))))

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
                     (m/inc n))) ys n)
       (v/emult weights n) n])
    [ys
     (map (fn [^double y ^double w]
            (m// (m/+ 0.5 (m/* w y))
                 (m/inc w))) ys weights)
     weights (repeat obs 1.0)]))

(defn- ylogymu ^double [^double y ^double mu] (if (m/zero? y) 0.0 (m/* y (m/log (m// y mu)))))
(defn- binomial-dev-resids
  [ys mus weights]
  (map (fn ^double [^double y ^double mu ^double w]
         (m/* 2.0 w (m/+ (ylogymu y mu)
                         (ylogymu (m/- 1.0 y)
                                  (m/- 1.0 mu))))) ys mus weights))

(defn- gaussian-dev-resids
  [ys mus weights]
  (map (fn [^double y ^double mu ^double w]
         (m/* w (m/sq (m/- y mu))))
       ys mus weights))

(defn- gamma-dev-resids
  [ys mus weights]
  (map (fn [^double y ^double mu ^double w]
         (m/* -2.0 w (m/- (m/log (if (m/zero? y) 1.0 (m// y mu)))
                          (m// (m/- y mu) mu))))
       ys mus weights))

(defn- poisson-initalize
  [ys weights ^long obs]
  [ys (v/shift ys 0.1) weights (repeat obs 1.0)])

(defn- poisson-dev-resids [ys mus weights]
  (map (fn [^double y ^double mu ^double w]
         (let [r (if (m/pos? y)
                   (m/-
                    (m/* y (m/log (m// y mu)))
                    (m/- y mu))
                   mu)]
           (m/* 2.0 w r)))
       ys mus weights))

(defn- inverse-gaussian-dev-resids
  [ys mus weights]
  (map (fn [^double y ^double mu ^double w]
         (m/* w (m// (m/sq (m/- y mu))
                     (m/* y mu mu))))
       ys mus weights))

(defn- nbinomial-initialize
  [ys weights ^long obs]
  [ys (map (fn ^double [^double y]
             (if (m/zero? y) m/SIXTH y))) weights (repeat obs 1.0)])

(defn- nbinomial-dev-resids
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
                          (fn ^double [^double x] (-> (m/constrain x -8.209536151601387 8.209536151601387)
                                                     (m// m/SQRT2)
                                                     (m/-)
                                                     (m/erfc)
                                                     (m/* 0.5)))
                          (fn ^double [^double x] (m/max m/MACHINE-EPSILON (m/* m/INV_SQRT2PI
                                                                               (m/exp (m/* -0.5 x x))))))
          :cauchit (->Link (fn ^double [^double x] (m/tan (m/* m/PI (m/- x 0.5))))
                           (fn ^double [^double x] (-> (m/constrain x
                                                                   -1.978937966095219E15
                                                                   1.978937966095219E15)
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
                          (fn ^double [^double x] (m/exp (m/- (m/- (m/exp (m/- x))) x))))
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
          :nbinomial (fn [{:keys [^double alpha]
                          :or {alpha 1.0}}]
                       (let [ralpha (m// alpha)]
                         (->Link (fn ^double [^double x] (m/log (m// x (m/+ x ralpha))))
                                 (fn ^double [^double x] (m// -1.0 (m/* alpha (m/- 1.0 (m/exp (m/- x))))))
                                 (fn ^double [^double x] (let [e (m/exp x)]
                                                          (m// e (m/* alpha (m/sq (m/- 1.0 e)))))))))
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

(defrecord Family [default-link variance initialize dev-resids dispersion])

(def families {:binomial (->Family :logit (fn ^double [^double x] (m/* x (m/- 1.0 x)))
                                 binomial-initialize binomial-dev-resids 1.0)
             :gaussian (->Family :identity constantly-1
                                 default-initialize gaussian-dev-resids :estimate)
             :gamma (->Family :inverse m/sq
                              default-initialize gamma-dev-resids :estimate)
             :poisson (->Family :log m/identity-double 
                                poisson-initalize poisson-dev-resids 1.0)
             :inverse-gaussian (->Family :inversesq m/cb
                                         default-initialize inverse-gaussian-dev-resids :estimate)
             :nbinomial (fn [{:keys [^double alpha]
                             :or {alpha 1.0}}]
                          (->Family :log (fn ^double [^double x] (m/+ x (m/* alpha x x)))
                                    nbinomial-initialize (nbinomial-dev-resids (m// alpha))
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


(defrecord GLMData [^double intercept beta coefficients
                    residual fitted weights deviance df dispersion
                    chi2 family link mean iters]
  IFn
  (invoke [_ xs] (mean (m/+ (v/dot beta xs) intercept)))
  prot/PredictProto
  (predict [_ xs] (mean (m/+ (v/dot beta xs) intercept))))

(defn- small-singular-values-count
  ^long [^doubles sv ^double epsilon]
  (let [f (Array/aget sv 0)]
    (count
     (filter (fn [^double v] (m/< (m// v f) epsilon)) sv))))

(defn glm-residual-analysis
  [residuals deviance pearson]
  {:residuals residuals :deviance deviance
   :pearson pearson})

(defn- glm-analysis
  [beta stderrs ^long df estimated?]
  (let [tdistr (if estimated?
                 (r/distribution :t {:degrees-of-freedom df})
                 r/default-normal)]
    (map (fn [^double b ^double err]
           (let [t-value (m// b err)]
             {:estimate b :stderr err :t-value t-value
              :p-value (stats/p-value tdistr t-value :both)})) beta stderrs)))

(defn glm
  [ys xss {:keys [^long max-iters ^double tol family link weights]
           :or {max-iters 25 tol 1.0e-8 family :gaussian}
           :as params}]
  
  (let [^RealMatrix xss (->> xss
                             (map (fn [xs] (conj (v/vec->seq xs) 1.0)))
                             (mat/mat))

        ^SingularValueDecomposition S (SingularValueDecomposition. xss)
        singular-values (.getSingularValues S)
        k (small-singular-values-count singular-values tol)]

    (when (m/pos? k) (ex-info "Near rank-deficient model matrix" {:data singular-values}))
    
    (let [uts (.getUT S)
          ut (.getU S)

          {link-fun :g link-mean :mean
           link-derivative :derivative link-variance :variance
           initialize :initialize
           dev-resids :dev-resids
           dispersion :dispersion
           link :link family :family} (family-with-link family link params)

          m (mat/nrow xss)
          n (mat/ncol xss)

          weights (or weights (repeat m 1.0))
          
          [ys start-t weights ns] (initialize ys weights m)
          
          ^doubles ys (m/seq->double-array ys)
          ^doubles weights (m/seq->double-array weights)
          
          init-t (double-array (map link-fun start-t))
          
          ^doubles buff-g (double-array m)
          ^doubles buff-z (double-array m)
          ^doubles buff-W (double-array m)
          ^doubles buff-Wz (double-array m)
          
          [^RealVector result t ^double dev iters]
          (loop [iter (long 0)
                 ^doubles t init-t
                 dev ##Inf]

            ;; iterate over t
            (dotimes [idx m]
              (let [eta (Array/aget t idx)
                    g (double (link-mean eta))
                    v (double (link-variance g))
                    g' (double (link-derivative eta))]
                (when (or (m/zero? v) (m/nan? v) (m/nan? g))
                  (ex-info "Invalid variance of mean." {:mean g :variance v :coeff idx}))
                
                (let [z (m/+ eta (m// (m/- (Array/aget ys idx) g) g'))
                      w (m/* (Array/aget weights idx) (m// (m/* g' g') v))]
                  (Array/aset buff-g idx g)
                  (Array/aset buff-z idx z)
                  (Array/aset buff-W idx w)
                  (Array/aset buff-Wz idx (m/* w z)))))

            (let [new-t (->> (-> (.multiply uts (.multiply (DiagonalMatrix. buff-W) ut))
                                 (CholeskyDecomposition.  tol tol)
                                 (.getSolver)
                                 (.solve (.operate uts (v/vec->RealVector buff-Wz))))
                             (.operate ut))
                  new-dev (v/sum (dev-resids ys buff-g weights))]

              (if (or (m/< (m// (m/abs (m/- new-dev dev))
                                (m/+ 0.1 (m/abs new-dev))) tol)
                      (m/== iter max-iters))
                [(->> singular-values
                      (v/reciprocal)
                      (v/vec->RealVector)
                      (.ebeMultiply (.operate uts new-t))
                      (.operate (.getV S)))
                 (v/vec->array new-t)
                 new-dev
                 iter]
                
                (recur (m/inc iter)
                       (v/vec->array new-t)
                       new-dev))))

          fitted (map link-mean (v/vec->seq (.operate xss result)))
          residuals (v/ediv (map m/- ys fitted) (map link-derivative t))
          pearson (map (fn [^double y ^double mu ^double w]
                         (m/* (m/- y mu)
                              (m/sqrt (m// w (double (link-variance mu))))))
                       ys fitted weights)
          chi2 (v/sum (map m/sq pearson))
          
          result-array (v/vec->array result)
          
          intercept (Array/aget result-array 0)
          beta (vec (rest result-array))

          null-df (m/dec m)
          df (m/- m n)

          estimated-dispersion? (= :estimate dispersion)
          dispersion-value (double (if estimated-dispersion?
                                     (estimate-dispersion residuals buff-W df)
                                     dispersion))

          stderrs (-> xss
                      (.transpose)
                      (.multiply (.multiply (DiagonalMatrix. buff-W) xss))
                      (CholeskyDecomposition.)
                      (.getSolver)
                      (.getInverse)
                      (mat/diag)
                      (v/mult dispersion-value)
                      (v/sqrt)
                      (v/vec->seq))]
      
      (->GLMData intercept beta
                 (glm-analysis result-array stderrs df estimated-dispersion?)
                 (glm-residual-analysis
                  residuals
                  (->> (dev-resids ys fitted weights)
                       (map m/safe-sqrt)
                       (map (fn [^double y ^double mu ^double r]
                              (if (m/<= y mu) (m/- r) r)) ys fitted))
                  pearson)
                 fitted (seq buff-W)
                 {:residuals dev :null (v/sum (dev-resids ys (repeat m (stats/mean ys weights)) weights))}
                 {:residuals df :null null-df}
                 dispersion-value chi2 family link link-mean iters))))

(comment
  (def weight [ 0.520 0.700 1.000 1.170 1.198 1.480 1.617 1.693 1.720 2.340 2.516 2.796 2.804 3.108 3.204 3.353 3.478 3.587 3.612 3.390 3.740])
  (def age [22 23 25 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44])
  (def births [   1   1   1   1   6   1   3   6   7   7   7  29  43 114 222 353 401 247  53 9   1])

  (predict (glm weight age {:family :gaussian :weights births}) [30]))


(comment
  (def data-xss [[712.0    21.0  105.0  82.4  13566.0  12.3             14952.0]
               [643.0    26.5   97.0  80.2  13566.0  15.3             17039.5]
               [679.0    28.3  113.0  86.3   9611.0  13.9             19215.7]
               [801.0    27.1  109.0  80.4   9483.0  13.6             21707.1]
               [753.0    22.0  115.0  64.7   9265.0  14.6             16566.0]
               [714.0    24.3  107.0  79.0   9555.0  13.8             17350.2]
               [920.0    21.2  118.0  72.2   9611.0  13.3             19504.0]
               [779.0    20.5  114.0  75.2   9483.0  14.5             15969.5]
               [771.0    23.2  102.0  81.1   9483.0  14.2             17887.2]
               [724.0    20.5  112.0  80.3  12656.0  13.7             14842.0]
               [682.0    23.8   96.0  83.0   9483.0  14.6             16231.6]
               [837.0    22.1  111.0  74.5  12656.0  11.6             18497.7]
               [599.0    19.9  117.0  83.8   8298.0  15.1             11920.1]
               [680.0    21.5  121.0  77.6   9265.0  13.7             14620.0]
               [747.0    22.5  109.0  77.9   8314.0  14.4             16807.5]
               [982.0    19.4  137.0  65.3   9483.0  13.3             19050.8]
               [719.0    25.9  109.0  80.9   8298.0  14.9             18622.1]
               [831.0    18.5  138.0  80.2   9483.0  14.6             15373.5]
               [858.0    19.4  119.0  84.8  12656.0  14.3             16645.2]
               [652.0    27.2  108.0  86.4  13566.0  14.6             17734.4]
               [718.0    23.7  115.0  73.5   9483.0  15.0             17016.6]
               [787.0    20.8  126.0  74.7   9483.0  14.9             16369.6]
               [515.0    26.8  106.0  87.8   8298.0  15.3             13802.0]
               [732.0    23.0  103.0  86.6   9611.0  13.8             16836.0]
               [783.0    20.5  125.0  78.5   9483.0  14.1             16051.5]
               [612.0    23.7  100.0  80.6   9033.0  13.3             14504.4]
               [486.0    23.2  117.0  84.8   8298.0  15.9             11275.2]
               [765.0    23.6  105.0  79.2   9483.0  13.7             18054.0]
               [793.0    21.7  125.0  78.4   9483.0  14.5             17208.1]
               [776.0    23.0  110.0  77.2   9265.0  13.6             17848.0]
               [978.0    19.3  130.0  71.5   9483.0  15.3             18875.4]
               [792.0    21.2  126.0  82.2  12656.0  15.1             16790.4]])

  (def data-ys [60.3
              52.3
              53.4
              57.0
              68.7
              48.8
              65.5
              70.5
              59.1
              62.7
              51.6
              62.0
              68.4
              69.2
              64.7
              75.0
              62.1
              67.2
              67.7
              52.7
              65.7
              72.2
              47.4
              51.3
              63.6
              50.7
              51.6
              56.2
              67.6
              58.9
              74.7
              67.3])

  (def glm-example (glm data-ys data-xss {:family :gamma}))

  (predict glm-example [792 21 126 82 12656 15 16790])
  ;; => 71.48845692808746

  glm-example
  ;; => {:intercept -0.01776527027537515,
  ;;     :beta
  ;;     [4.96176829942184E-5
  ;;      0.002034422589585567
  ;;      -7.18142873678752E-5
  ;;      1.1185201293319783E-4
  ;;      -1.4675150420118665E-7
  ;;      -5.186831119354037E-4
  ;;      -2.4271749790784295E-6],
  ;;     :coefficients
  ;;     ({:estimate -0.01776527027537515,
  ;;       :stderr 0.011479217035549256,
  ;;       :t-value -1.5476029611043172,
  ;;       :p-value 0.13480435421482173}
  ;;      {:estimate 4.96176829942184E-5,
  ;;       :stderr 1.6215765104251235E-5,
  ;;       :t-value 3.0598422384158916,
  ;;       :p-value 0.005380919142093488}
  ;;      {:estimate 0.002034422589585567,
  ;;       :stderr 5.320801860007495E-4,
  ;;       :t-value 3.8235263088385354,
  ;;       :p-value 8.220108346288502E-4}
  ;;      {:estimate -7.18142873678752E-5,
  ;;       :stderr 2.7116639036702992E-5,
  ;;       :t-value -2.6483476536554886,
  ;;       :p-value 0.014073122618301244}
  ;;      {:estimate 1.1185201293319783E-4,
  ;;       :stderr 4.057690945483059E-5,
  ;;       :t-value 2.756543424227744,
  ;;       :p-value 0.01098067148177595}
  ;;      {:estimate -1.4675150420118665E-7,
  ;;       :stderr 1.2365685050481715E-7,
  ;;       :t-value -1.1867640458420847,
  ;;       :p-value 0.24693566328846595}
  ;;      {:estimate -5.186831119354037E-4,
  ;;       :stderr 2.402533747068508E-4,
  ;;       :t-value -2.1589004215582137,
  ;;       :p-value 0.04107317540274753}
  ;;      {:estimate -2.4271749790784295E-6,
  ;;       :stderr 7.460253329340098E-7,
  ;;       :t-value -3.253475280166025,
  ;;       :p-value 0.003373346804228008}),
  ;;     :residual
  ;;     {:residuals
  ;;      (-7.46911736792169E-4
  ;;       3.429625820340195E-4
  ;;       -0.0011094607670473703
  ;;       3.9091242724339814E-4
  ;;       3.5558471032353127E-4
  ;;       0.0024991959298770474
  ;;       2.953769069182333E-4
  ;;       -0.0010241699624186028
  ;;       -3.488351298312974E-4
  ;;       1.3309859496109607E-4
  ;;       8.203818425814939E-4
  ;;       -1.8902553083192738E-4
  ;;       -8.545972669001706E-4
  ;;       -0.0014207843004860694
  ;;       -0.0010865654964093644
  ;;       -2.944816406450419E-5
  ;;       -0.0017280691781362792
  ;;       9.289405904117575E-4
  ;;       -8.071950066136936E-4
  ;;       -2.495954668027426E-4
  ;;       -3.5142380463519264E-4
  ;;       -2.0727164934354302E-4
  ;;       -7.961111480553285E-4
  ;;       0.0012039765934349761
  ;;       7.53464030978887E-4
  ;;       4.894737413637737E-4
  ;;       0.0014522125574446458
  ;;       5.691341311035412E-4
  ;;       -4.857408216191105E-5
  ;;       4.3496186941919396E-4
  ;;       -1.6034147394135498E-4
  ;;       4.7771919901831445E-4),
  ;;      :deviance
  ;;      (0.04256856307258896
  ;;       -0.018383238680803323
  ;;       0.05508223793521002
  ;;       -0.022977588507490765
  ;;       -0.025268886860163545
  ;;       -0.14953290844286152
  ;;       -0.019868087616261642
  ;;       0.06616299065719772
  ;;       0.020073264893048
  ;;       -0.008439839773358937
  ;;       -0.044965039673061515
  ;;       0.011540922031546816
  ;;       0.054396012535400134
  ;;       0.08760809483907067
  ;;       0.0645547114120105
  ;;       0.0022021389967450647
  ;;       0.09474141566511653
  ;;       -0.06845344861830359
  ;;       0.051075224043423566
  ;;       0.012929286993512508
  ;;       0.022411014677977792
  ;;       0.014675644488433509
  ;;       0.03597724035997922
  ;;       -0.06765533826344133
  ;;       -0.05134214298909052
  ;;       -0.025684193710770496
  ;;       -0.08392236097387545
  ;;       -0.03345155572404423
  ;;       0.003269332290466995
  ;;       -0.026545925794190878
  ;;       0.011790996437022772
  ;;       -0.03363248313457027),
  ;;      :pearson
  ;;      (0.04317472117963417
  ;;       -0.018270763852120536
  ;;       0.05609819722808162
  ;;       -0.022801936668716875
  ;;       -0.025056497674960043
  ;;       -0.1421742901606863
  ;;       -0.01973672574632739
  ;;       0.06763014573796869
  ;;       0.0202078009536464
  ;;       -0.008416112859729973
  ;;       -0.04429362861290146
  ;;       0.011585362292139705
  ;;       0.05538675988537233
  ;;       0.09018494863358936
  ;;       0.06595122370896916
  ;;       0.0022037557653639996
  ;;       0.09775671971316503
  ;;       -0.06690048219337352
  ;;       0.05194845950301177
  ;;       0.012985069081454453
  ;;       0.02257874427164603
  ;;       0.014747523629403002
  ;;       0.036409981653792195
  ;;       -0.06613827000207885
  ;;       -0.050467256331928154
  ;;       -0.02546477336992414
  ;;       -0.08159130979883032
  ;;       -0.0330795979698636
  ;;       0.003272896105229182
  ;;       -0.026311551870628457
  ;;       0.011837384433197845
  ;;       -0.03325649666431903)},
  ;;     :fitted
  ;;     (57.80431482447356
  ;;      53.27334470063797
  ;;      50.56347993033019
  ;;      58.330037828448134
  ;;      70.46562168593833
  ;;      56.888012844871625
  ;;      66.81878401480326
  ;;      66.03410392769389
  ;;      57.929374726164475
  ;;      63.23216907126933
  ;;      53.991478496798656
  ;;      61.28993391078229
  ;;      64.81036393466701
  ;;      63.475468164125324
  ;;      60.69696113756204
  ;;      74.83508175712625
  ;;      56.5699110602814
  ;;      72.01804171752491
  ;;      64.3567651897932
  ;;      52.02445880844704
  ;;      64.24933079045786
  ;;      71.15070332151728
  ;;      45.7347968844956
  ;;      54.93318587979207
  ;;      66.98031260545203
  ;;      52.024799734864004
  ;;      56.18413735686392
  ;;      58.122674712420014
  ;;      67.37947398203181
  ;;      60.49162862428671
  ;;      73.82609216583235
  ;;      69.61515621029368),
  ;;     :weights
  ;;     (3341.3388123272575
  ;;      2838.0492555932015
  ;;      2556.6655026645876
  ;;      3402.393313068035
  ;;      4965.403839585406
  ;;      3236.2460054380454
  ;;      4464.749897216312
  ;;      4360.502881533235
  ;;      3355.812456164202
  ;;      3998.3072054583213
  ;;      2915.0797502701803
  ;;      3756.4559987881644
  ;;      4200.383273344021
  ;;      4029.1350586547396
  ;;      3684.1210913343866
  ;;      5600.289461594642
  ;;      3200.1548373675573
  ;;      5186.598332826947
  ;;      4141.793225694223
  ;;      2706.5443143118105
  ;;      4127.9765070213625
  ;;      5062.422583146225
  ;;      2091.671646066016
  ;;      3017.654910903618
  ;;      4486.362276723918
  ;;      2706.5797874526147
  ;;      3156.6572905348594
  ;;      3378.2453157254436
  ;;      4539.993514094771
  ;;      3659.237133618213
  ;;      5450.291884477022
  ;;      4846.26997418398),
  ;;     :deviance
  ;;     {:residuals 0.08738851641700845, :null 0.5360720799622841},
  ;;     :df {:residuals 24, :null 31},
  ;;     :dispersion 0.003584283173493368,
  ;;     :chi2 0.08602279616383574,
  ;;     :family :gamma,
  ;;     :link :inverse,
  ;;     :mean #function[fastmath.core//]}
  
  )


(comment
  (def dist [   2  10   4  22  16  10  18  26  34  17  28  14  20  24  28  26  34  34  46
           26  36  60  80  20  26  54  32  40  32  40  50  42  56  76  84  36  46  68
           32  48  52  56  64  66  54  70  92  93 120  85])

  (def speed [  4  4  7  7  8  9 10 10 10 11 11 12 12 12 12 13 13 13 13 14 14 14 14 15 15
            15 16 16 17 17 17 18 18 18 18 19 19 19 20 20 20 20 20 22 23 24 24 24 24 25])

  (def ctl [4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14])
  (def trt [4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69])
  (def weight (concat ctl trt))
  (def group (concat (repeat 10 0) (repeat 10 1)))

  (lm group weight)

  (dissoc (:residual (lm group weight)) :residuals)
  ;; => {:variance 0.45943421052631583,
  ;;     :skewness 0.49076129003476227,
  ;;     :kurtosis 2.545241591026403,
  ;;     :durbin-watson 2.7630301572300033,
  ;;     :jarque-berra {:p-value 0.6141108002291311, :Z 0.9751598214321114},
  ;;     :normality {:p-value 0.5591574636998055, :Z 1.162648314601708}}

  (dissoc (lm group weight) :residual :fitted)
  ;; => {:intercept 5.031999999999999,
  ;;     :beta [-0.37099999999999966],
  ;;     :adjusted-r-squared 0.021581910045406882,
  ;;     :r-squared 0.07307759899038546,
  ;;     :f-statistic 1.4191012973623132,
  ;;     :rss 8.729249999999999,
  ;;     :p-value 0.24902316597300622,
  ;;     :df {:residuals 18, :model 1},
  ;;     :tss 9.417454999999999,
  ;;     :ess 0.6882049999999983,
  ;;     :coefficients
  ;;     ({:estimate 5.031999999999999,
  ;;       :stderr 0.2202176953229084,
  ;;       :t-value 22.850116529561824,
  ;;       :p-value 9.547918011776346E-15}
  ;;      {:estimate -0.37099999999999966,
  ;;       :stderr 0.3114348514002032,
  ;;       :t-value -1.1912603818487013,
  ;;       :p-value 0.2490231659730061}),
  ;;     :stderr 0.6963894982933999,
  ;;     :model
  ;;     #object[org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression 0x4a1f3279 "org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression@4a1f3279"]}

  (prot/predict (lm group weight) [1])

  (dissoc (lm speed dist) :residual :fitted)
  ;; => {:intercept -17.57909489051092,
  ;;     :beta [3.932408759124086],
  ;;     :adjusted-r-squared 0.6438102011907144,
  ;;     :r-squared 0.6510793807582509,
  ;;     :df 48,
  ;;     :coefficients
  ;;     ({:estimate -17.57909489051092,
  ;;       :std-err 6.7584401693792335,
  ;;       :t-value -2.601058003022252,
  ;;       :p-value 0.012318816153808889}
  ;;      {:estimate 3.932408759124086,
  ;;       :std-err 0.4155127766571222,
  ;;       :t-value 9.463989990298368,
  ;;       :p-value 1.48991929904696E-12}),
  ;;     :model
  ;;     #object[org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression 0x3c3845ab "org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression@3c3845ab"]}


  (def xss [[ 0.00000000e+00,  0.00000000e+00,  2.50000000e+01],
          [ 4.08163265e-01,  3.96924149e-01,  2.10849646e+01],
          [ 8.16326531e-01,  7.28634783e-01,  1.75031237e+01],
          [ 1.22448980e+00,  9.40632785e-01,  1.42544773e+01],
          [ 1.63265306e+00,  9.98087482e-01,  1.13390254e+01],
          [ 2.04081633e+00,  8.91559230e-01,  8.75676801e+00],
          [ 2.44897959e+00,  6.38550320e-01,  6.50770512e+00],
          [ 2.85714286e+00,  2.80629400e-01,  4.59183673e+00],
          [ 3.26530612e+00, -1.23398137e-01,  3.00916285e+00],
          [ 3.67346939e+00, -5.07151709e-01,  1.75968347e+00],
          [ 4.08163265e+00, -8.07581691e-01,  8.43398584e-01],
          [ 4.48979592e+00, -9.75328286e-01,  2.60308205e-01],
          [ 4.89795918e+00, -9.82831204e-01,  1.04123282e-02],
          [ 5.30612245e+00, -8.28857736e-01,  9.37109538e-02],
          [ 5.71428571e+00, -5.38705288e-01,  5.10204082e-01],
          [ 6.12244898e+00, -1.60045086e-01,  1.25989171e+00],
          [ 6.53061224e+00,  2.44910071e-01,  2.34277384e+00],
          [ 6.93877551e+00,  6.09627196e-01,  3.75885048e+00],
          [ 7.34693878e+00,  8.74184299e-01,  5.50812162e+00],
          [ 7.75510204e+00,  9.95115395e-01,  7.59058726e+00],
          [ 8.16326531e+00,  9.52551848e-01,  1.00062474e+01],
          [ 8.57142857e+00,  7.53486727e-01,  1.27551020e+01],
          [ 8.97959184e+00,  4.30625870e-01,  1.58371512e+01],
          [ 9.38775510e+00,  3.70144015e-02,  1.92523948e+01],
          [ 9.79591837e+00, -3.62678429e-01,  2.30008330e+01],
          [ 1.02040816e+01, -7.02784220e-01,  2.70824656e+01],
          [ 1.06122449e+01, -9.27424552e-01,  3.14972928e+01],
          [ 1.10204082e+01, -9.99691655e-01,  3.62453145e+01],
          [ 1.14285714e+01, -9.07712248e-01,  4.13265306e+01],
          [ 1.18367347e+01, -6.66598288e-01,  4.67409413e+01],
          [ 1.22448980e+01, -3.15964115e-01,  5.24885464e+01],
          [ 1.26530612e+01,  8.65820672e-02,  5.85693461e+01],
          [ 1.30612245e+01,  4.74903061e-01,  6.49833403e+01],
          [ 1.34693878e+01,  7.85198826e-01,  7.17305289e+01],
          [ 1.38775510e+01,  9.66488646e-01,  7.88109121e+01],
          [ 1.42857143e+01,  9.88987117e-01,  8.62244898e+01],
          [ 1.46938776e+01,  8.48997803e-01,  9.39712620e+01],
          [ 1.51020408e+01,  5.69520553e-01,  1.02051229e+02],
          [ 1.55102041e+01,  1.96472687e-01,  1.10464390e+02],
          [ 1.59183673e+01, -2.08855085e-01,  1.19210746e+02],
          [ 1.63265306e+01, -5.79868557e-01,  1.28290296e+02],
          [ 1.67346939e+01, -8.55611267e-01,  1.37703040e+02],
          [ 1.71428571e+01, -9.90779466e-01,  1.47448980e+02],
          [ 1.75510204e+01, -9.63165404e-01,  1.57528113e+02],
          [ 1.79591837e+01, -7.77305991e-01,  1.67940441e+02],
          [ 1.83673469e+01, -4.63737404e-01,  1.78685964e+02],
          [ 1.87755102e+01, -7.39780734e-02,  1.89764681e+02],
          [ 1.91836735e+01,  3.27935645e-01,  2.01176593e+02],
          [ 1.95918367e+01,  6.75970465e-01,  2.12921699e+02],
          [ 2.00000000e+01,  9.12945251e-01,  2.25000000e+02]])

  (def ys [ 4.48933467,  5.34231784,  5.81681264,  6.69809894,  5.94032168,
         6.65120265,  6.59592036,  6.07958661,  6.36186999,  6.7775208 ,
         6.6836138 ,  6.56551043,  7.17869566,  7.25158876,  8.04903067,
         7.84930564,  8.39547421,  8.25715349,  9.42432217,  8.43407786,
         9.70661341,  8.17868934, 10.44401028,  8.72852257,  9.32853046,
         9.70458442,  8.58017153,  9.78906458, 10.20670793,  9.20752194,
         9.26055989, 10.66574522,  9.89388631, 11.87907875, 10.65688892,
         11.09346875, 10.77674283, 10.11844363, 10.75104318, 10.19883781,
         10.68264587, 10.73214468, 10.13324477,  9.6745733 , 10.05157595,
         10.74993283, 10.4821192 , 10.38354809, 11.52270194, 11.19546067])

  (:residual (lm ys xss))
  ;; => {:model
  ;;     #object[org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression 0x5a0f9d15 "org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression@5a0f9d15"],
  ;;     :intercept 5.205844957505818,
  ;;     :beta
  ;;     [0.4687244824101058 0.48360119133943785 -0.017404792025998428],
  ;;     :coefficients
  ;;     ({:estimate 5.205844957505818,
  ;;       :stderr 0.17121765098613817,
  ;;       :t-value 30.404838096553984,
  ;;       :p-value 0.0}
  ;;      {:estimate 0.4687244824101058,
  ;;       :stderr 0.02640602380582178,
  ;;       :t-value 17.75066499435501,
  ;;       :p-value 0.0}
  ;;      {:estimate 0.48360119133943785,
  ;;       :stderr 0.10380518021321035,
  ;;       :t-value 4.65873851715441,
  ;;       :p-value 2.7357988015230816E-5}
  ;;      {:estimate -0.017404792025998428,
  ;;       :stderr 0.0023184680741556492,
  ;;       :t-value -7.507022511982179,
  ;;       :p-value 1.5906164352043508E-9}),
  ;;     :residual
  ;;     {:skewness 0.20739432426143098,
  ;;      :durbin-watson 2.895610842471915,
  ;;      :df 46,
  ;;      :residuals
  ;;      (-0.281390486855857
  ;;       0.1201831997789844
  ;;       0.18060503046724197
  ;;       0.7115107121796109
  ;;       -0.3161106545012595
  ;;       0.21002773477713443
  ;;       0.04664027043646435
  ;;       -0.5212639045929288
  ;;       -0.2624545488744907
  ;;       0.12571689950623455
  ;;       -0.03016566394419673
  ;;       -0.26861126500748433
  ;;       0.15253688652965813
  ;;       -0.03889808824475338
  ;;       0.43415861572664927
  ;;       -0.12695489929314618
  ;;       0.050908099821604935
  ;;       -0.43045985438778267
  ;;       0.4478982783433221
  ;;       -0.7558996851690925
  ;;       0.38794759174036564
  ;;       -1.1871812187495188
  ;;       1.0966019242028828
  ;;       -0.5604093189118498
  ;;       0.10681516991152051
  ;;       0.5270685634295127
  ;;       -0.6031849766769817
  ;;       0.5319787291205937
  ;;       0.8022621531175034
  ;;       -0.41060627454152687
  ;;       -0.6184156884468752
  ;;       0.5066167977434599
  ;;       -0.5327165039046751
  ;;       1.22853381728563
  ;;       -0.14941147324256043
  ;;       0.2140037158843029
  ;;       -0.09150837760802233
  ;;       -0.6653379860026725
  ;;       0.10278114237650193
  ;;       -0.29249481169484604
  ;;       0.337447353275687
  ;;       0.49280638651751474
  ;;       -0.0624160472758728
  ;;       -0.5503326284207137
  ;;       -0.27330354015398584
  ;;       0.2691187059729234
  ;;       -0.1856763668723076
  ;;       -0.47130760645851844
  ;;       0.5126412325563283
  ;;       0.08970285913042986),
  ;;      :ss 11.607740226533197,
  ;;      :kurtosis 3.0259454940918746,
  ;;      :ms 0.2523421788376782,
  ;;      :normality {:p-value 0.7208635552500905, :Z 0.6546108067291813},
  ;;      :jarque-berra
  ;;      {:p-value 0.8353373910613211, :Z 0.35983914918144544}},
  ;;     :fitted
  ;;     (4.770725156855857
  ;;      5.2221346402210145
  ;;      5.636207609532758
  ;;      5.9865882278203895
  ;;      6.256432334501259
  ;;      6.441174915222865
  ;;      6.549280089563535
  ;;      6.600850514592928
  ;;      6.6243245388744905
  ;;      6.651803900493766
  ;;      6.7137794639441974
  ;;      6.834121695007484
  ;;      7.026158773470343
  ;;      7.290486848244754
  ;;      7.614872054273352
  ;;      7.976260539293145
  ;;      8.344566110178395
  ;;      8.687613344387781
  ;;      8.976423891656678
  ;;      9.189977545169095
  ;;      9.318665818259635
  ;;      9.365870558749517
  ;;      9.347408355797118
  ;;      9.28893188891185
  ;;      9.22171529008848
  ;;      9.177515856570487
  ;;      9.18335650667698
  ;;      9.257085850879404
  ;;      9.404445776882497
  ;;      9.618128214541526
  ;;      9.878975578446877
  ;;      10.15912842225654
  ;;      10.426602813904674
  ;;      10.650544932714368
  ;;      10.80630039324256
  ;;      10.879465034115698
  ;;      10.868251207608022
  ;;      10.783781616002674
  ;;      10.648262037623496
  ;;      10.491332621694847
  ;;      10.345198516724313
  ;;      10.239338293482485
  ;;      10.19566081727587
  ;;      10.224905928420714
  ;;      10.324879490153986
  ;;      10.480814124027077
  ;;      10.667795566872307
  ;;      10.854855696458518
  ;;      11.010060707443671
  ;;      11.10575781086957),
  ;;     :df {:residuals 46, :model 3},
  ;;     :observations 50,
  ;;     :r-squared 0.9325049619539212,
  ;;     :adjusted-r-squared 0.9281031116465681,
  ;;     :stderr 0.5023367185839377,
  ;;     :tss 171.97916413660798,
  ;;     :rss 11.607740226533197,
  ;;     :ess 160.37142391007478,
  ;;     :mse 53.45714130335826,
  ;;     :f-statistic 211.84386038667412,
  ;;     :p-value 0.0,
  ;;     :log-likelihood -34.438154937049326,
  ;;     :aic 76.87630987409865,
  ;;     :bic 84.52440189581124}

  (defn convert-x [x] [x (m/sin x) (m/sq (- x 5))])

  (predict (lm xss ys) (convert-x 2.123))
  ;; => 6.468609171503275
  )
