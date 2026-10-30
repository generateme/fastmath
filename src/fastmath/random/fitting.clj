(ns fastmath.random.fitting
  "Fitting of probability distributions to data.

  A distribution is fitted by its key (see [[fastmath.random/distribution]]) to a sample. What the engine needs to know about the distribution comes from [[fastmath.random/fit-spec]]: the fitted parameters with their domains and roles, the data validation and, optionally, closed form estimators and starting values. Distributions without a specification can't be fitted.

  The methods are the maximum likelihood estimate (`:mle`), the minimum goodness-of-fit estimates (`:ks`, `:cvm`, `:ad` and its variants), the maximum product of spacings (`:mps`), quantile matching (`:qme`), moment matching (`:mme`) and L-moment matching (`:lmme`). An estimate is a closed form one when the specification has an exact estimator for the method, otherwise it is found numerically: the box constrained quasi-Newton method L-BFGS-B (the simplex method for the minimum distance methods) started from the starting values of the specification or, when there are none, from the best points of a scan of the parameter space with a quantile matching objective.

  The result of [[fit]] is a map with the fitted distribution, its parameters and statistics of the fit (information criteria and, for continuous distributions, goodness-of-fit statistics)."
  (:refer-clojure :exclude [compare])
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.stats :as stats]
            [fastmath.optimization :as opt]
            [fastmath.calculus :as calc]
            [fastmath.polynomials :as poly]
            [fastmath.matrix :as mat]
            [fastmath.stats.bootstrap :as boot]
            [fastmath.vector :as v])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^{:private true :const true :tag 'double} lower-limit 1.0e-10)

(def ^{:private true :const true :tag 'double} l-moment-tail-limit 36.0)
(def ^{:private true :const true :tag 'double} l-moment-tail-limit-start 30.0)

(def ^:private seed-probabilities (vec (range 0.025 0.9751 0.05)))

(def ^:private methods-list [:mle :ks :cvm :ad :adr :adl :ad2r :ad2l :ad2 :mps :qme :mme :lmme])

;; objectives which are not smooth (:ks) or have plateaus and a huge gradient far from the optimum (the Anderson-Darling variants)
(def ^:private distance-methods #{:ks :cvm :ad :adr :adl :ad2r :ad2l :ad2})

(defn- invalid-input!
  "Throws `ex-info` for an invalid `property` of `fit`."
  [property message value expected]
  (throw (ex-info (str "fit, " (name property) ": " message)
                  {:type :invalid-input
                   :function 'fit
                   :property property
                   :value value
                   :expected expected})))

(defn- domain-bounds
  "Bounds `[lo hi]` of the domain of a parameter: the open ends are moved a little inside."
  [domain]
  (case domain
    :real [##-Inf ##Inf]
    :positive [lower-limit ##Inf]
    :non-negative [0.0 ##Inf]
    :unit [lower-limit (m/- 1.0 lower-limit)]
    domain))

(defn- in-domain?
  [domain ^double v]
  (let [[^double lo ^double hi] (domain-bounds domain)]
    (and (m/<= lo v) (m/<= v hi))))

(defn- robust-scale
  "Scale of the sorted data `xs`: IQR / 1.349, then the standard deviation, then `1.0` for constant data."
  ^double [xs]
  (let [iqr-scale (m// (stats/iqr xs) 1.349)]
    (cond
      (m/pos? iqr-scale) iqr-scale
      (m/pos? (stats/stddev xs)) (stats/stddev xs)
      :else 1.0)))

(defn- shape-range
  [domain]
  (case domain
    :real [-2.0 2.0]
    :positive [0.05 20.0]
    :non-negative [0.0 20.0]
    :unit [0.01 0.99]
    domain))

(defn- bound-role?
  "Is `role` the role of a parameter which is an end of the support of the distribution?"
  [role]
  (contains? #{:lower-bound :upper-bound} role))

(defn- codec
  "How a free parameter is seen by the optimizer: a parameter which is an end of the support is searched as the distance of the end from the sample minimum (the lower end) or maximum (the upper end), which can't be negative, the others as they are."
  [{:keys [role]} ^double xmin ^double xmax]
  (case role
    :lower-bound {:decode (fn [^double v] (m/- xmin v)) :encode (fn [^double a] (m/max 0.0 (m/- xmin a)))}
    :upper-bound {:decode (fn [^double v] (m/+ xmax v)) :encode (fn [^double a] (m/max 0.0 (m/- a xmax)))}
    {:decode identity :encode identity}))

(defn- search-bounds
  "Bounds of the free parameter as the optimizer sees it."
  [{:keys [role domain]}]
  (if (bound-role? role) [0.0 ##Inf] (domain-bounds domain)))

(defn- seed-range
  "Finite range of the scan for a parameter, from its role and the sorted data `xs`, inside its domain."
  [{:keys [domain role]} ^doubles xs ^double scale]
  (let [[^double dlo ^double dhi] (domain-bounds domain)
        mn (stats/minimum xs)
        mx (stats/maximum xs)
        [^double lo ^double hi] (case role
                                  :location (if (m/< mn mx) [mn mx] [(m/- mn 1.0) (m/+ mx 1.0)])
                                  :scale [(m/* 0.01 scale) (m/* 10.0 scale)]
                                  (:lower-bound :upper-bound) [0.0 scale]
                                  :probability [0.01 0.99]
                                  (shape-range domain))
        [^double dlo ^double dhi] (if (bound-role? role) [0.0 ##Inf] [dlo dhi])
        lo (m/max lo dlo)
        hi (m/min hi dhi)]
    (if (m/< lo hi) [lo hi] [(m/max dlo -1.0) (m/min dhi 1.0)])))

(defn- log-likelihood
  "Log likelihood of the sample, `##-Inf` outside of the support, `NaN` when a density is `NaN`."
  ^double [d ^doubles xs]
  (let [n (alength xs)]
    (loop [i (long 0)
           s 0.0]
      (if (m/< i n)
        (let [s (m/+ s (r/lpdf d (Array/aget xs i)))]
          (if (m/valid-double? s)
            (recur (m/inc i) s)
            s))
        s))))

(defn- objective-function
  "Function of a point (values of the free parameters as the optimizer sees them, `decode` turns them into the parameters) returning `(objective d params)` for the distribution `d` with the parameters `params`, `##Inf` when the value is not finite or the distribution can't be created."
  [dist-key fixed names decode objective]
  (fn ^double [point]
    (try
      (let [params (merge fixed (zipmap names (decode point)))
            value (double (objective (r/distribution dist-key params) params))]
        (if (m/valid-double? value) value ##Inf))
      (catch Exception _ ##Inf))))

(defn- quantile-matching
  "Objective: the mean squared distance between the quantiles of the sample (with the estimation `strategy`) and the quantiles of the distribution at the probabilities `probs`."
  [xs probs strategy]
  (let [data-quantiles (stats/quantiles xs probs strategy)]
    (fn [d _]
      (m// (v/sum (map (fn [^double p ^double q] (m/sq (m/- (double (r/icdf d p)) q)))
                       probs data-quantiles))
           (count probs)))))

(defn- spacings
  "Objective: minus the mean logarithm of the spacings of the cumulative distribution function at the sorted sample, with the ends `0` and `1`. A spacing between tied observations is zero: it is replaced by the density (mass) at the observation. A zero last spacing is left out."
  [^doubles xs]
  (let [n (alength xs)]
    (fn [d _]
      (loop [i (long 0)
             previous 0.0
             sum 0.0
             cnt (long 0)]
        (let [last? (m/== i n)
              current (if last? 1.0 (r/cdf d (Array/aget xs i)))
              diff (m/- current previous)
              term (double (cond
                             (m/pos? diff) (m/log diff)
                             last? ##NaN
                             :else (r/lpdf d (Array/aget xs i))))
              included? (not (m/nan? term))
              sum (if included? (m/+ sum term) sum)
              cnt (if included? (m/inc cnt) cnt)]
          (cond
            (not last?) (recur (m/inc i) current sum cnt)
            (m/zero? cnt) ##Inf
            :else (m/- (m// sum cnt))))))))

(defn- undefined-moments!
  [d what]
  (throw (ex-info (str "fit: the " what " of the distribution don't exist") {:type :moments-undefined :function 'fit :distribution d})))

(defn- raw-moments
  "First `k` raw moments of the distribution: from `mean` and `variance` up to the order 2, above by the integral of the quantile function raised to the order."
  [d ^long k]
  (let [mean (r/mean d)
        variance (if (m/> k 1) (r/variance d) 0.0)]
    (when-not (and (m/valid-double? mean) (m/valid-double? variance)) (undefined-moments! d "moments"))
    (let [res (mapv (fn [^long j]
                      (case j
                        1 mean
                        2 (m/+ variance (m/sq mean))
                        (calc/integrate (fn [^double u] (m/fpow (r/icdf d u) j)) 0.0 1.0)))
                    (range 1 (m/inc k)))]
      (when-not (every? m/valid-double? res) (undefined-moments! d "moments"))
      res)))

(defn- l-moment-integral
  "The L-moment of the order `r`: the integral of the quantile function times the shifted Legendre polynomial. Each half of the unit interval is integrated in the variable `v = -ln(distance from the end)`, `v` from `ln 2` to `tail-limit`, where the heavy tails of the quantile function become exponentially decaying integrands."
  ^double [d ^long r]
  (let [ln2 (m/log 2.0)
        part (fn ^double [upper?]
               (double (calc/integrate (fn ^double [^double v]
                                         (let [e (m/exp (m/- v))
                                               u (if upper? (m/- 1.0 e) e)
                                               x (m/* (double (r/icdf d u))
                                                      (poly/eval-legendre-P (m/long-dec r) (m/dec (m/* 2.0 u)))
                                                      e)]
                                           ;; the quantile function is infinite at the end, the weight underflows
                                           (if (and (m/invalid-double? x) (m/> v l-moment-tail-limit-start)) 0.0 x)))
                                       ln2 l-moment-tail-limit)))]
    (m/+ (double (part false)) (double (part true)))))

(defn- l-moments
  "First `k` L-moments of the distribution; they exist when the mean is finite."
  [d ^long k]
  (when-not (m/valid-double? (r/mean d)) (undefined-moments! d "L-moments"))
  (let [res (mapv #(l-moment-integral d %) (range 1 (m/inc k)))]
    (when-not (every? m/valid-double? res) (undefined-moments! d "L-moments"))
    res))

(defn- moment-matching
  "Objective: the sum of squared relative differences of the first `k` raw moments of the sample and of the distribution (`spec-moments` computes the latter for a map of parameters, or `nil` for the quantile function integrals)."
  [^doubles xs ^long k spec-moments]
  (let [orders (range 1 (m/inc k))
        sd (stats/stddev xs)
        sample (mapv (fn [^long j] (stats/mean (map (fn [^double x] (m/fpow x j)) xs))) orders)
        scales (mapv (fn [^double m ^long j] (m/+ (m/abs m) (m/fpow sd j))) sample orders)]
    (fn ^double [d params]
      (let [theoretical (if spec-moments (vec (take k (spec-moments params))) (raw-moments d k))]
        (v/sum (map (fn [^double t ^double s ^double c] (m/sq (m// (m/- t s) c)))
                    theoretical sample scales))))))

(defn- l-moment-matching
  "Objective: the sum of squared differences of the first `k` L-moments of the sample and of the distribution, in units of the second L-moment of the sample."
  [xs ^long k spec-l-moments]
  (let [sample (mapv (fn [^long r] (stats/l-moment xs r {:sorted? true})) (range 1 (m/inc k)))
        l2 (stats/l-moment xs 2 {:sorted? true})]
    (fn ^double [d params]
      (let [theoretical (if spec-l-moments (vec (take k (spec-l-moments params))) (l-moments d k))]
        (v/sum (map (fn [^double t ^double s] (m/sq (m// (m/- t s) l2))) theoretical sample))))))

(defn- method-objective
  "Objective of `method`, a function of a distribution and the map of its parameters returning the value to minimize."
  [method xs spec k {:keys [quantiles strategy] :or {quantiles 50 strategy :r7}}]
  (case method
    :mle (let [^doubles xs xs] (fn [d _] (m/- (log-likelihood d xs))))
    :ks (fn [d _] (stats/ks-statistic xs d {:sorted? true}))
    :cvm (fn [d _] (stats/cvm-statistic xs d {:sorted? true}))
    (:ad :adr :adl :ad2r :ad2l :ad2) (fn [d _] (stats/ad-statistic xs d {:variant method :sorted? true}))
    :mps (spacings xs)
    :qme (let [probs (if (number? quantiles) (m/slice-range 0.01 0.99 (long quantiles)) (vec quantiles))]
           (when-not (and (seq probs) (every? #(and (number? %) (m/< 0.0 (double %) 1.0)) probs))
             (invalid-input! :quantiles "probabilities have to be numbers inside (0, 1), or a positive count" quantiles "probabilities or a count"))
           (when (m/< (count probs) (long k))
             (invalid-input! :quantiles "not enough quantiles to fit the parameters" quantiles (str "at least " k)))
           (quantile-matching xs probs strategy))
    :mme (moment-matching xs k (:moments spec))
    :lmme (l-moment-matching xs k (:lmoments spec))))

(defn- clamp-to-bounds
  [bounds point]
  (mapv (fn [[^double lo ^double hi] ^double v] (m/constrain v lo hi)) bounds point))

(defn- scanned-starts
  "Starting points: the best results of a scan of the parameter space with the quantile matching objective."
  [dist-key fixed free decode xs scale restarts rng]
  (let [names (mapv :name free)
        bounds (mapv #(seed-range % xs scale) free)
        res (opt/scan-and-minimize :lbfgsb
                                   (objective-function dist-key fixed names decode (quantile-matching xs seed-probabilities :r7))
                                   {:bounds bounds :N (m/* 50 (count free)) :n restarts :take-last-n restarts
                                    :rng rng :stats? true})]
    (map :point (if (map? res) [res] res))))

(defn- simplex-bounds
  "Bounds around `start` which set the size of the initial simplex of `:nelder-mead` to a fifth of the value of each parameter."
  [start]
  (mapv (fn [^double x] (let [h (m/max (m/* 0.1 (m/abs x)) 1.0e-3)] [(m/- x h) (m/+ x h)])) start))

(defn- polish
  "Minimizes the objective from `start`, inside the domain `bounds`; `nil` when the optimization fails. A quasi-Newton step which ends above the starting value (a very steep or ill-conditioned objective, for example with the ends of the support as parameters) is repeated with the simplex method, and when nothing is better than the start the start is the result, as not converged."
  [optimizer objective bounds start]
  (let [start (clamp-to-bounds bounds start)
        start-value (double (objective start))
        run (fn [method]
              (try
                (let [res (if (= :lbfgsb method)
                            (opt/minimize method objective {:bounds bounds :initial start :stats? true})
                            (opt/minimize method objective {:bounds (simplex-bounds start) :length 1.0 :initial start :stats? true}))]
                  (when (m/valid-double? (:value res)) res))
                (catch Exception _ nil)))
        first-result (run optimizer)]
    (if (and first-result (m/<= (double (:value first-result)) start-value))
      first-result
      (let [fallback (when (= :lbfgsb optimizer) (run :nelder-mead))
            best (first (sort-by :value (filter #(and % (m/<= (double (:value %)) start-value)) [first-result fallback])))]
        (cond
          best best
          (m/valid-double? start-value) {:point start :value start-value :iterations 0 :status :max-iterations}
          :else nil)))))

(defn- information-criteria
  [^double ll ^long k ^long n]
  (let [aic (m/- (m/* 2.0 k) (m/* 2.0 ll))
        dof (m/- n k 1)]
    {:aic aic
     :aicc (if (m/pos? dof) (m/+ aic (m// (m/* 2.0 k (m/inc k)) dof)) ##Inf)
     :bic (m/- (m/* k (m/log n)) (m/* 2.0 ll))}))

(defn- automatic-breaks
  "Bounds of the cells of the chi-square statistic for the sorted sample `xs`, by the rule of `fitdistrplus::gofstat`: a bound is the observation at the position `round(n / (4 n)^(2/5))` of the observations which are left, the observations up to it (ties included) are taken off, until at most 1.5 times that number is left. Returns `nil` when the sample is too small."
  [^doubles xs]
  (let [n (alength xs)
        cell (m/round (m// n (m/pow (m/* 4.0 n) 0.4)))
        rest-limit (m/round (m/ceil (m/* 1.5 cell)))]
    (when (m/> n rest-limit)
      (loop [rest-xs (seq xs)
             breaks []]
        (if (m/> (count rest-xs) rest-limit)
          (let [limit (double (nth rest-xs (dec cell)))]
            (recur (drop-while (fn [^double v] (m/<= v limit)) rest-xs) (conj breaks limit)))
          breaks)))))

(defn- chi-square-statistics
  "`:chisq` and `:chisq-df` of the fitted discrete distribution `d` for the sorted sample `xs` and the cells bounded by `breaks` (or the automatic ones), as in `fitdistrplus::gofstat`: the cells are `(-Inf, b1]`, `(b1, b2]`, ..., and `(b_last, Inf)` unless the distribution and the sample end there. `k` is the number of fitted parameters. Empty when there are no cells."
  [d ^doubles xs ^long k breaks]
  (when-let [breaks (seq (or breaks (automatic-breaks xs)))]
    (let [n (alength xs)
          observed-cumulative (mapv (fn [^double b] (long (count (take-while (fn [^double v] (m/<= v b)) xs)))) breaks)
          theoretical-cumulative (mapv (fn [^double b] (r/cdf d b)) breaks)
          closed? (and (m/== 1.0 (double (peek theoretical-cumulative))) (m/== n (long (peek observed-cumulative))))
          differences (fn [cumulative end] (map m/- (if closed? cumulative (conj cumulative end)) (cons 0 cumulative)))
          observed (differences observed-cumulative n)
          expected (map #(m/* n ^double %) (differences theoretical-cumulative 1.0))
          cells (count observed)]
      {:chisq (v/sum (map (fn [^double o ^double e] (m// (m/sq (m/- o e)) e)) observed expected))
       :chisq-df (m/- cells 1 k)})))

(defn- fit-statistics
  [d ^doubles xs ^long k breaks]
  (let [n (alength xs)
        ll (log-likelihood d xs)]
    (merge {:k k :n n :loglik ll}
           (information-criteria ll k n)
           (if (r/continuous? d)
             {:ks (stats/ks-statistic xs d {:sorted? true})
              :cvm (stats/cvm-statistic xs d {:sorted? true})
              :ad (stats/ad-statistic xs d {:sorted? true})}
             (chi-square-statistics d xs k breaks)))))



(defn- build-result
  [dist-key spec fixed method params tier converged? iterations statistics]
  (let [d (try (r/distribution dist-key params)
               (catch Exception e
                 (throw (ex-info "fit: the fitted parameters don't make a distribution"
                                 {:type :not-converged :function 'fit :params params} e))))]
    {:distribution d
     :params (into {} (map (fn [{:keys [name]}] [name (double (get params name))])) (:params spec))
     :method method
     :tier tier
     :converged? converged?
     :iterations iterations
     :fixed fixed
     :stats (when statistics (statistics d))}))

(defn- estimate
  "Parameters from the closed form estimator `closed-form` (an entry of `:estimators`), or `nil`. An estimator declared for fixed parameters applies when exactly these are fixed, any other one when `fixed-only?` is false or no parameter is fixed."
  [closed-form xs fixed exact?]
  (when closed-form
    (if-let [required (:fixed closed-form)]
      (when (= required (set (keys fixed)))
        ((:estimate closed-form) xs fixed))
      (when (or (not exact?) (empty? fixed))
        ((:estimate closed-form) xs)))))

(defn- fit-fixed
  "Fit with the integer parameters fixed (`fixed` holds all the fixed parameters, `shown-fixed` the ones given by the user, `k` is the number of fitted parameters). Returns the map of `:result` and the `:value` of the objective of the method at it."
  [{:keys [dist-key spec xs method optimizer restarts closed-form? rng opts stats?]} fixed shown-fixed k]
  (let [statistics (when stats? #(fit-statistics % xs k (:breaks opts)))
        estimators (when closed-form? (:estimators spec))
        closed-form (get estimators method)
        free (vec (remove #(contains? fixed (:name %)) (:params spec)))
        free-names (mapv :name free)
        objective (method-objective method xs spec (count free) opts)
        value-of (fn [result] (double (objective (:distribution result) (:params result))))
        estimated (estimate closed-form xs fixed (:exact? closed-form))]
    (cond
      (and (:exact? closed-form) estimated)
      (let [result (build-result dist-key spec shown-fixed method (merge estimated fixed) :closed-form true 0 statistics)]
        {:result result :value (value-of result)})
      ;; every parameter is fixed (the searched integer was the last free one)
      (empty? free)
      (let [result (build-result dist-key spec shown-fixed method fixed :closed-form true 0 statistics)]
        {:result result :value (value-of result)})
      :else
      (let [bounds (mapv search-bounds free)
            scale (robust-scale xs)
            codecs (let [^doubles sorted xs]
                     (mapv #(codec % (Array/aget sorted 0) (Array/aget sorted (m/dec (alength sorted)))) free))
            decode (fn [point] (mapv (fn [c v] ((:decode c) v)) codecs point))
            objective-fn (objective-function dist-key fixed free-names decode objective)
            seed (or estimated
                     (when-let [initial (:initial spec)] (initial xs))
                     (estimate (:mle estimators) xs fixed false))
            starts (if (and seed (every? #(contains? seed %) free-names))
                     [(mapv (fn [c nm] ((:encode c) (seed nm))) codecs free-names)]
                     (scanned-starts dist-key fixed free decode xs scale restarts rng))
            best (->> starts
                      (keep #(polish optimizer objective-fn bounds %))
                      (filter #(m/valid-double? (double (:value %))))
                      (sort-by :value)
                      (first))]
        (when (or (nil? best) (m/inf? (:value best)))
          (if (#{:mme :lmme} method)
            (undefined-moments! dist-key "moments to match")
            (throw (ex-info "fit: no usable estimate was found" {:type :not-converged :function 'fit :distribution dist-key}))))
        (let [result (build-result dist-key spec shown-fixed method (merge fixed (zipmap free-names (decode (:point best))))
                                   :numeric (contains? #{nil :converged :stalled} (:status best)) (or (:iterations best) 0) statistics)]
          {:result result :value (:value best)})))))

(defn- standard-errors
  "Standard errors of the fitted parameters of `result` as a map to merge into its `:stats`: `:se` (a map from a parameter name to its standard error, `nil` when it can't be computed) and `:se-warning` (the reason, when `:se` is `nil`).

  They are the square roots of the diagonal of the inverse of the negative Hessian of the log likelihood at the estimate (central finite differences), the usual asymptotic result of the maximum likelihood method. Fixed parameters, integer parameters and the ends of the support (their log likelihood is not smooth) are held at their values and have no entry; the others are conditional on them."
  [dist-key spec xs {:keys [params fixed method]}]
  (let [names (into [] (comp (remove #(or (contains? fixed (:name %)) (= :integer (:role %)) (bound-role? (:role %))))
                             (map :name))
                    (:params spec))
        failure (fn [message] {:se nil :se-warning message})]
    (cond
      (not= :mle method) (failure "standard errors are available for the maximum likelihood method only")
      (empty? names) {:se {}}
      :else
      (let [loglik (fn ^double [point]
                     (try (log-likelihood (r/distribution dist-key (merge params (zipmap names point))) xs)
                          (catch Exception _ ##NaN)))
            hessian ((calc/hessian loglik) (mapv params names))]
        (if-not (every? (fn [row] (every? (fn [^double v] (m/valid-double? v)) row)) hessian)
          (failure "the log likelihood is not finite around the estimate")
          (let [information (into-array (map #(double-array (map m/- %)) hessian))
                variances (try (map-indexed (fn [^long i ^doubles row] (Array/aget row i))
                                 (mat/mat->array2d (mat/inverse (mat/cholesky-decomposition information))))
                               (catch Exception _ nil))]
            (if (and variances (every? (fn [^double v] (and (m/valid-double? v) (m/pos? v))) variances))
              {:se (zipmap names (map m/sqrt variances))}
              (failure "the Hessian of the log likelihood is singular or not positive definite at the estimate"))))))))

(defn- search-integer
  "Integer from `lo` to `hi` where `g` (a function of a long returning a number) is the smallest, found from `start`: the direction of descent is found, the step doubles while `g` decreases and the bracket is shrunk by the ternary search. `g` has to be unimodal (a plateau is allowed)."
  [g lo hi start]
  (let [lo (long lo)
        hi (long hi)
        g (memoize g)
        value (fn [i] (g i))
        start (long (m/constrain (long start) lo hi))
        argmin (fn [a b] (reduce (fn [best i] (if (m/< (double (value i)) (double (value best))) i best)) (range a (m/inc (long b)))))
        direction (long (cond
                    (and (m/> start lo) (m/< (double (value (m/dec start))) (double (value start)))) -1
                    (and (m/< start hi) (m/< (double (value (m/inc start))) (double (value start)))) 1
                    :else 0))]
    (if (or (m/== lo hi) (m/zero? direction))
      start
      (let [[a b] (loop [before start
                         previous start
                         k 1]
                    (let [next (long (m/constrain (m/+ start (m/* direction (m/dec (m/bit-shift-left 1 k)))) lo hi))]
                      (if (and (m/not== next previous) (m/< (double (value next)) (double (value previous))))
                        (recur previous next (m/inc k))
                        [(m/min before next) (m/max before next)])))]
        (loop [a (long a)
               b (long b)]
          (if (m/<= (m/- b a) 2)
            (argmin a b)
            (let [third (m/quot (m/- b a) 3)
                  m1 (m/+ a third)
                  m2 (m/- b third)]
              (if (m/< (double (value m1)) (double (value m2))) (recur a m2) (recur m1 b)))))))))

(defn fit
  "Fits a distribution to data.

  The distribution is given by its key (see [[fastmath.random/distribution]]), the fitted parameters are the ones listed by [[fastmath.random/fit-spec]]. When the specification of the distribution has an exact closed form estimator for the method (and the fixed parameters are the ones it is declared for, none for most estimators, and `:closed-form?` is not `false`) it is the result. Otherwise the objective of the method is minimized, from a closed form estimate that is not exact or the starting values of the specification or the maximum likelihood estimate of the specification, or from the best points of a scan of the parameter space (a quantile matching objective over ranges derived from the data and the roles of the parameters) when it has none.

  Methods (the `:method` option):

  - `:mle`: maximum likelihood, the default.
  - `:ks`, `:cvm`, `:ad`, `:adr`, `:adl`, `:ad2r`, `:ad2l`, `:ad2`: minimum distance, the minimum of the Kolmogorov-Smirnov, Cramer-von Mises or Anderson-Darling statistic or of one of its variants (see [[fastmath.stats/ks-statistic]], [[fastmath.stats/cvm-statistic]], [[fastmath.stats/ad-statistic]]).
  - `:mps`: maximum product of spacings, the maximum of the mean logarithm of the spacings of the cumulative distribution function at the sorted sample. A spacing between tied observations is zero and is replaced by the density (mass) at the observation, so the method can be used for discrete distributions.
  - `:qme`: quantile matching, the minimum of the mean squared distance of the sample quantiles and the distribution quantiles.
  - `:mme`: moment matching, the first `k` raw moments (`k` is the number of fitted parameters) of the distribution equal to the ones of the sample.
  - `:lmme`: L-moment matching, the first `k` L-moments of the distribution equal to the ones of the sample.

  Integer parameters (the `:integer` role, for example the trials of the binomial distribution) can't be found by the optimizer. When the specification gives a `:range` for the one which is not fixed, it is searched: for each integer the other parameters are fitted and the integer with the best value of the objective is chosen (the objective has to be unimodal in the integer). Otherwise an integer parameter has to be fixed.

  Parameters:

  - `dist-key` (keyword): A key of [[fastmath.random/distribution]] which has a specification.
  - `data` (sequence of numbers): The sample; it has to pass the validation of the specification (for example `:gamma` requires positive values and `:poisson` integers).
  - `opts` (map, optional):
    - `:method` (keyword, default `:mle`): The estimation method.
    - `:fixed` (map, default `nil`): Parameters which are not fitted, as a map from the name of a parameter to its value. At least one parameter has to stay free (or be searched).
    - `:optimizer` (keyword): The method of [[fastmath.optimization]] used for the numerical minimization, `:lbfgsb` or `:nelder-mead` (its initial simplex is a fifth of the starting values). The default is `:nelder-mead` for the minimum distance methods (`:ks`, `:cvm`, `:ad` and its variants), whose objectives are not smooth or have plateaus, and `:lbfgsb` otherwise.
    - `:stats?` (boolean, default `true`): Set to `false` to skip all the statistics (`:stats` is `nil`, `:se?` is ignored), as [[bootstrap]] does for its refits.
    - `:breaks` (sequence of numbers, default `nil`): For discrete distributions, the increasing bounds of the cells of the chi-square statistic (see `:chisq` below); the cells are `(-Inf, b1]`, `(b1, b2]`, ..., `(b_last, Inf)`. By default the bounds are chosen by the rule of `fitdistrplus::gofstat`: the observation at the position `round(n / (4 n)^(2/5))` among the observations which are left, until at most 1.5 times that number is left; no bounds (and no chi-square) for a small sample.
    - `:se?` (boolean, default `true`): Set to `false` to skip the standard errors (they take a few dozens of evaluations of the likelihood, slow for distributions with a numerically integrated density).
    - `:restarts` (long, default `3`): The number of starting points taken from the scan, when there are no other starting values.
    - `:closed-form?` (boolean, default `true`): Set to `false` to skip the closed form estimators and fit numerically.
    - `:rng`: The random number generator of the scan. The default is a new unseeded generator; pass a seeded one for a reproducible scan. The scan only provides starting points, so the fit is the same up to the tolerance of the optimizer. Fits which start from a closed form estimate or from the starting values of the specification don't use it.
    - `:quantiles` (number or sequence, default `50`) and `:strategy` (keyword, default `:r7`): For `:qme`. Either the probabilities of the matched quantiles (numbers inside `(0, 1)`, at least as many as the fitted parameters) or their count, spread evenly from `0.01` to `0.99`; and the quantile estimation strategy of [[fastmath.stats/quantiles]].

  Returns a map with the keys:

  - `:distribution`: The fitted distribution object.
  - `:params`: The map of all parameters (fitted, searched and fixed) to doubles.
  - `:method`: The estimation method.
  - `:tier`: `:closed-form` or `:numeric`. A fit in which the search of an integer parameter leaves no parameter to optimize (the chi distribution) is `:closed-form`.
  - `:converged?`: `true` when the optimizer reported convergence or a stalled objective (it stopped changing, the usual end at the optimum), `false` when it stopped at the iteration limit. A closed form result is always converged.
  - `:iterations`: The number of iterations of the optimizer, `0` for a closed form.
  - `:fixed`: The map of fixed parameters, as given.
  - `:stats`: A map with `:k` (the number of fitted parameters, a searched integer included), `:n`, `:loglik`, `:aic`, `:aicc`, `:bic`, and for continuous distributions the goodness-of-fit statistics `:ks`, `:cvm` and `:ad`, which are descriptive, not tests. They describe the fitted distribution, whatever the method. for discrete distributions `:chisq` (the Pearson chi-square statistic over the cells of `:breaks`, absent when the sample is too small for the automatic cells) and `:chisq-df` (cells - 1 - `k`, may be not positive, then the statistic is of no use), also descriptive, no p-value. `:aicc` is `##Inf` when `n` is not above `k + 1`. Unless `:se?` is `false` it holds also `:se`, the map from the name of a fitted parameter to its standard error, from the inverse of the negative Hessian of the log likelihood at the estimate, a maximum likelihood result (the other methods give `:se` equal to `nil`). Fixed parameters, a searched integer and the ends of the support (the `:lower-bound` and `:upper-bound` roles, whose log likelihood is not smooth) have no entry and the other standard errors are conditional on them. When the Hessian is not finite, singular or not positive definite (a flat ridge of the likelihood), `:se` is `nil` and `:se-warning` gives the reason; it never throws.

  Moment matching with the moments above the second order, and L-moment matching, integrate the quantile function numerically unless the specification has its own `:moments` or `:lmoments`.

  Throws `ex-info` with `:type :not-fittable` when the distribution has no specification, with `:type :invalid-input` for data which fail the validation, an unknown `:method`, wrong `:quantiles`, a wrong `:fixed` value, an integer parameter which has to be fixed or more than one integer parameter to search, with `:type :moments-undefined` when `:mme` or `:lmme` finds no parameters for which the moments exist (the Cauchy distribution has none), and with `:type :not-converged` when no usable estimate is found. An evaluation of the objective which throws or is not finite counts as an infinite objective value.

  See also [[fastmath.random/fit-spec]], [[fastmath.random/log-likelihood]]."
  ([dist-key data] (fit dist-key data nil))
  ([dist-key data {:keys [method fixed optimizer restarts closed-form? rng se? stats? breaks]
                   :or {method :mle restarts 3 closed-form? true se? true stats? true}
                   :as opts}]
   (when-not (some #{method} methods-list)
     (invalid-input! :method "unknown method" method (set methods-list)))
   (when-not (contains? #{nil :lbfgsb :nelder-mead} optimizer)
     (invalid-input! :optimizer "unsupported optimizer" optimizer #{:lbfgsb :nelder-mead}))
   (when-not (or (nil? breaks) (and (sequential? breaks) (every? number? breaks) (apply < breaks)))
     (invalid-input! :breaks "the bounds of the cells have to be increasing numbers" breaks "increasing numbers"))
   (let [spec (or (r/fit-spec dist-key)
                  (throw (ex-info (str "fit: the distribution " dist-key " can't be fitted")
                                  {:type :not-fittable :function 'fit :distribution dist-key})))
         optimizer (or optimizer (if (contains? distance-methods method) :nelder-mead :lbfgsb))
         names (set (map :name (:params spec)))
         fixed (or fixed {})
         _ (when-let [message ((:valid? spec) data)]
             (invalid-input! :data message data "data accepted by the specification"))
         _ (doseq [[nm v] fixed
                   :let [prm (first (filter #(= nm (:name %)) (:params spec)))]]
             (cond
               (not (names nm)) (invalid-input! :fixed (str nm " is not a parameter") nm names)
               (not (number? v)) (invalid-input! :fixed (str nm " is not a number") v "number")
               (not (in-domain? (:domain prm) v)) (invalid-input! :fixed (str nm " is outside of its domain") v (:domain prm))))
         xs (m/seq->double-array (sort data))
         context {:dist-key dist-key :spec spec :xs xs :method method :optimizer optimizer
                  :restarts restarts :closed-form? closed-form? :rng (or rng (r/rng :jdk)) :opts opts :stats? stats?}
         k (m/long-sub (count (:params spec)) (count fixed))
         closed-form (when closed-form? (get-in spec [:estimators method]))
         exact? (and (:exact? closed-form) (estimate closed-form xs fixed true))
         searched (remove #(contains? fixed (:name %)) (filter #(= :integer (:role %)) (:params spec)))
         result (cond
                 ;; the exact estimator also finds the integer parameters it is declared for
                 exact? (:result (fit-fixed context fixed fixed k))
                 (empty? searched) (do (when (m/zero? k)
                                         (invalid-input! :fixed "at least one parameter has to stay free" fixed "fewer fixed parameters"))
                                       (:result (fit-fixed context fixed fixed k)))
                 (m/> (count searched) 1) (invalid-input! :fixed "at most one integer parameter can be searched, fix the others" (mapv :name searched) "one integer parameter")
                 (not (:range (first searched))) (invalid-input! :fixed (str (:name (first searched)) " is an integer parameter and has to be fixed") (:name (first searched)) "a value in :fixed")
                 :else (let [{:keys [name range]} (first searched)
                             [lo hi start] (range xs)
                             fit-at (memoize (fn [i] (fit-fixed context (assoc fixed name i) fixed k)))
                             best (search-integer #(:value (fit-at %)) lo hi start)]
                         (:result (fit-at best))))]
     (if (and se? stats?)
       (update result :stats merge (standard-errors dist-key spec xs result))
       result))))

(def ^:private ci-functions
  {:percentile boot/ci-percentile
   :basic boot/ci-basic
   :normal boot/ci-normal
   :bc boot/ci-bc
   :mad-median (fn [boot-data _] (boot/ci-mad-median boot-data))})

(defn bootstrap
  "Bootstrap of the parameters of a fitted distribution: the nonparametric bootstrap with confidence intervals.

  The distribution is fitted to `data` with [[fit]], then to `:samples` resamples of `data` drawn with replacement (see [[fastmath.stats.bootstrap/bootstrap]]). The refits run in parallel (`pmap`), skip the statistics, and each gets its own generator seeded from `:rng`, so the result is the same for a seeded `:rng`. A refit which throws (a resample which the distribution can't be fitted to) is left out.

  Parameters:

  - `dist-key` (keyword): A key of [[fastmath.random/distribution]] which has a specification, see [[fit]].
  - `data` (sequence of numbers): The sample.
  - `opts` (map, optional): The options of [[fit]] (`:method`, `:fixed`, `:optimizer`, `:restarts`, `:closed-form?` and the options of the method) are used for every fit; and
    - `:samples` (long, default `200`): The number of resamples, at least `2`.
    - `:alpha` (double, default `0.05`): The significance level, inside `(0, 1)`; the intervals have the level `1 - alpha`.
    - `:ci-type` (keyword, default `:percentile`): The interval, `:percentile`, `:basic`, `:normal`, `:bc` or `:mad-median`, the median of the estimates -/+ their median absolute deviation (see [[fastmath.stats.bootstrap/ci-mad-median]]; the interval of `fitdistr`, it has no level and ignores `:alpha`); the others are described in [[fastmath.stats.bootstrap/ci-percentile]], [[fastmath.stats.bootstrap/ci-basic]], [[fastmath.stats.bootstrap/ci-normal]] and [[fastmath.stats.bootstrap/ci-bc]].
    - `:rng`: The random number generator of the resampling and of the seeds of the refits. The default is a new unseeded generator.

  Returns a map with the keys:

  - `:fit`: The result of [[fit]] for `data`.
  - `:replicates`: The number of refits which succeeded.
  - `:alpha` and `:ci-type`: As used.
  - `:params`: A map from the name of every parameter which was not fixed (a searched integer included) to a map with `:t0` (the estimate from `data`), `:ts` (the estimates of the resamples), `:ci` (the vector `[lower upper]` of the interval), `:bias` (the mean of `:ts` minus `:t0`) and `:stddev` (the standard deviation of `:ts`, the bootstrap standard error).

  Throws what [[fit]] throws for `data`, `ex-info` with `:type :invalid-input` for a wrong `:samples`, `:alpha` or `:ci-type`, and with `:type :not-converged` when fewer than two refits succeed.

  See also [[fit]], [[fastmath.stats.bootstrap/bootstrap]]."
  ([dist-key data] (bootstrap dist-key data nil))
  ([dist-key data {:keys [samples alpha ci-type rng]
                   :or {samples 200 alpha 0.05 ci-type :percentile}
                   :as opts}]
   (when-not (and (integer? samples) (m/>= (long samples) 2))
     (invalid-input! :samples "at least 2 resamples are needed" samples "an integer >= 2"))
   (when-not (and (number? alpha) (m/< 0.0 (double alpha) 1.0))
     (invalid-input! :alpha "the level has to be inside (0, 1)" alpha "number inside (0, 1)"))
   (when-not (contains? ci-functions ci-type)
     (invalid-input! :ci-type "unknown interval" ci-type (set (keys ci-functions))))
   (let [fit-opts (dissoc opts :samples :alpha :ci-type :rng)
         result (fit dist-key data fit-opts)
         rng (or rng (r/rng :jdk))
         resamples (mapv doall (:samples (boot/bootstrap data {:samples samples :rng rng})))
         seeds (doall (repeatedly samples #(r/irandom rng)))
         refit-opts (assoc fit-opts :stats? false)
         refits (->> (pmap (fn [resample seed]
                             (try (fit dist-key resample (assoc refit-opts :rng (r/rng :jdk seed)))
                                  (catch Exception _ nil)))
                           resamples seeds)
                     (remove nil?)
                     (doall))]
     (when (m/< (count refits) 2)
       (throw (ex-info "bootstrap: fewer than two refits succeeded"
                       {:type :not-converged :function 'bootstrap :distribution dist-key})))
     {:fit result
      :replicates (count refits)
      :alpha alpha
      :ci-type ci-type
      :params (into {}
                    (map (fn [nm]
                           (let [t0 (double (get-in result [:params nm]))
                                 ts (mapv #(double (get-in % [:params nm])) refits)
                                 ats (m/seq->double-array ts)]
                             [nm {:t0 t0
                                  :ts ts
                                  :ci (vec (take 2 ((ci-functions ci-type) {:t0 t0 :ts ts} alpha)))
                                  :bias (m/- (stats/mean ats) t0)
                                  :stddev (stats/stddev ats)}])))
                    (remove #(contains? (:fixed result) %) (keys (:params result))))})))

(defn- plotting-positions
  "Plotting positions of `n` ordered observations, `(i - a) / (n + 1 - 2 a)` with `a = 3/8` for `n <= 10` and `1/2` otherwise (`ppoints` of R)."
  [^long n]
  (let [a (if (m/<= n 10) 0.375 0.5)]
    (mapv (fn [^long i] (m// (m/- (m/inc i) a) (m/- (m/inc n) (m/* 2.0 a)))) (range n))))

(defn- fitted-distribution
  "The distribution of a result of [[fit]], or the distribution itself."
  [result-or-distribution]
  (if (r/distribution? result-or-distribution)
    result-or-distribution
    (or (:distribution result-or-distribution)
        (invalid-input! :distribution "a distribution or a result of fit is expected" result-or-distribution "distribution or fit result"))))

(defn qq
  "Data of the quantile-quantile plot of a sample against a distribution.

  The sample is sorted and its `i`-th value is paired with the quantile of the distribution at the `i`-th plotting position `(i - a) / (n + 1 - 2 a)` (`i` from 1, `a = 3/8` for `n <= 10`, `1/2` otherwise, the convention of `ppoints` and `qqnorm` of R).

  Parameters:

  - `result-or-distribution`: A result of [[fit]] or a distribution object.
  - `data` (sequence of numbers): The sample.

  Returns a vector of pairs `[theoretical-quantile sample-value]`, ordered by the sample. A sample from the distribution lies along the diagonal.

  See also [[pp]], [[fit]]."
  [result-or-distribution data]
  (let [d (fitted-distribution result-or-distribution)
        xs (sort data)]
    (mapv (fn [^double p ^double x] [(double (r/icdf d p)) x]) (plotting-positions (count xs)) xs)))

(defn pp
  "Data of the probability-probability plot of a sample against a distribution.

  The sample is sorted and the cumulative probability of the distribution at its `i`-th value is paired with the `i`-th plotting position `(i - a) / (n + 1 - 2 a)` (see [[qq]]).

  Parameters:

  - `result-or-distribution`: A result of [[fit]] or a distribution object.
  - `data` (sequence of numbers): The sample.

  Returns a vector of pairs `[theoretical-probability empirical-probability]`, ordered by the sample. A sample from the distribution lies along the diagonal.

  See also [[qq]], [[fit]]."
  [result-or-distribution data]
  (let [d (fitted-distribution result-or-distribution)
        xs (sort data)]
    (mapv (fn [^double p ^double x] [(r/cdf d x) p]) (plotting-positions (count xs)) xs)))

(def ^:private criteria #{:aic :aicc :bic})

(defn compare
  "Compares results of [[fit]] fitted to the same sample in a table.

  The rows are ordered from the best to the worst by the information criterion (the smallest value is the best).

  Parameters:

  - `results` (sequence of results of [[fit]]): At least one; all fitted to the same data, all continuous or all discrete (see `(fastmath.random/continuous? d)`; a density and a mass log likelihood can't be compared), with the statistics (not fitted with `:stats? false`).
  - `opts` (map, optional):
    - `:criterion` (keyword, default `:aic`): `:aic`, `:aicc` or `:bic`, the criterion of the order and of the `:rank`.

  Returns a vector of rows, maps with the keys: `:distribution` (the key of the distribution), `:method`, `:k`, `:n`, `:loglik`, `:aic`, `:aicc`, `:bic`, `:delta-aic` and `:delta-bic` (the difference from the smallest value), `:akaike-weight` (relative likelihood of the model, `exp(-delta-aic / 2)` normalized to sum `1`) and `:rank` (from `1`). Continuous results add `:ks`, `:cvm` and `:ad` and discrete results add `:chisq` and `:chisq-df`, a value which a result doesn't have is `nil`. There are no p-values.

  Throws `ex-info` with `:type :invalid-input` and `:property :kind` for results of both kinds, `:property :n` for results of different sample sizes, `:property :results` for no results or a result without the statistics and `:property :criterion` for an unknown criterion.

  See also [[fit]], [[qq]], [[pp]]."
  ([results] (compare results nil))
  ([results {:keys [criterion] :or {criterion :aic}}]
   (when-not (contains? criteria criterion)
     (invalid-input! :criterion "unknown criterion" criterion criteria))
   (when (or (empty? results) (not-every? :stats results))
     (invalid-input! :results "results of fit with the statistics are expected" (count results) "non empty results with :stats"))
   (when (not= 1 (count (distinct (map #(r/continuous? (:distribution %)) results))))
     (invalid-input! :kind "results of continuous and discrete distributions can't be compared" (mapv #(r/distribution-id (:distribution %)) results) "all continuous or all discrete"))
   (when (not= 1 (count (distinct (map #(get-in % [:stats :n]) results))))
     (invalid-input! :n "results have to be fitted to the same sample" (mapv #(get-in % [:stats :n]) results) "the same sample size"))
   (let [continuous? (r/continuous? (:distribution (first results)))
         min-aic (double (apply min (map #(get-in % [:stats :aic]) results)))
         min-bic (double (apply min (map #(get-in % [:stats :bic]) results)))
         weights (map #(m/exp (m/* -0.5 (m/- ^double (get-in % [:stats :aic]) min-aic))) results)
         total (v/sum weights)
         rows (map (fn [{:keys [distribution method stats]} ^double w]
                     (merge {:distribution (r/distribution-id distribution)
                             :method method}
                            (select-keys stats [:k :n :loglik :aic :aicc :bic])
                            {:delta-aic (m/- ^double (:aic stats) min-aic)
                             :delta-bic (m/- ^double (:bic stats) min-bic)
                             :akaike-weight (m// w total)}
                            (zipmap (if continuous? [:ks :cvm :ad] [:chisq :chisq-df])
                                    (map stats (if continuous? [:ks :cvm :ad] [:chisq :chisq-df])))))
                   results weights)]
     (vec (map-indexed (fn [^long i row] (assoc row :rank (m/inc i)))
                       (sort-by criterion rows))))))
