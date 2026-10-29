(ns fastmath.random.fitting
  "Fitting of probability distributions to data.

  A distribution is fitted by its key (see [[fastmath.random/distribution]]) to a sample. What the engine needs to know about the distribution comes from [[fastmath.random/fit-spec]]: the fitted parameters with their domains and roles, the data validation and, optionally, closed form estimators and starting values. Distributions without a specification can't be fitted.

  The methods are the maximum likelihood estimate (`:mle`), the minimum goodness-of-fit estimates (`:ks`, `:cvm`, `:ad` and its variants), the maximum product of spacings (`:mps`), quantile matching (`:qme`), moment matching (`:mme`) and L-moment matching (`:lmme`). An estimate is a closed form one when the specification has an exact estimator for the method, otherwise it is found numerically: the box constrained quasi-Newton method L-BFGS-B (the simplex method for the minimum distance methods) started from the starting values of the specification or, when there are none, from the best points of a scan of the parameter space with a quantile matching objective.

  The result of [[fit]] is a map with the fitted distribution, its parameters and statistics of the fit (information criteria and, for continuous distributions, goodness-of-fit statistics)."
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.stats :as stats]
            [fastmath.optimization :as opt]
            [fastmath.calculus :as calc]
            [fastmath.polynomials :as poly]
            [fastmath.vector :as v])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^{:private true :const true :tag 'double} lower-limit 1.0e-10)

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

(defn- seed-range
  "Finite range of the scan for a parameter, from its role and the sorted data `xs`, inside its domain."
  [{:keys [domain role]} ^doubles xs ^double scale]
  (let [[^double dlo ^double dhi] (domain-bounds domain)
        mn (stats/minimum xs)
        mx (stats/maximum xs)
        [^double lo ^double hi] (case role
                                  :location (if (m/< mn mx) [mn mx] [(m/- mn 1.0) (m/+ mx 1.0)])
                                  :scale [(m/* 0.01 scale) (m/* 10.0 scale)]
                                  :probability [0.01 0.99]
                                  (shape-range domain))
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
  "Function of a point (values of the free parameters) returning `(objective d params)` for the distribution `d` with the parameters `params`, `##Inf` when the value is not finite or the distribution can't be created."
  [dist-key fixed names objective]
  (fn ^double [point]
    (try
      (let [params (merge fixed (zipmap names point))
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

(defn- l-moments
  "First `k` L-moments of the distribution: the integrals of the quantile function times the shifted Legendre polynomials."
  [d ^long k]
  (when-not (m/valid-double? (r/mean d)) (undefined-moments! d "L-moments"))
  (let [res (mapv (fn [^long r]
                    (calc/integrate (fn [^double u] (m/* (double (r/icdf d u)) (poly/eval-legendre-P (m/long-dec r) (m/dec (m/* 2.0 u))))) 0.0 1.0))
                  (range 1 (m/inc k)))]
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
  [dist-key fixed free xs scale restarts rng]
  (let [names (mapv :name free)
        bounds (mapv #(seed-range % xs scale) free)
        res (opt/scan-and-minimize :lbfgsb
                                   (objective-function dist-key fixed names (quantile-matching xs seed-probabilities :r7))
                                   {:bounds bounds :N (m/* 50 (count free)) :n restarts :take-last-n restarts
                                    :rng rng :stats? true})]
    (map :point (if (map? res) [res] res))))

(defn- simplex-bounds
  "Bounds around `start` which set the size of the initial simplex of `:nelder-mead` to a fifth of the value of each parameter."
  [start]
  (mapv (fn [^double x] (let [h (m/max (m/* 0.1 (m/abs x)) 1.0e-3)] [(m/- x h) (m/+ x h)])) start))

(defn- polish
  "Minimizes the objective from `start`, inside the domain `bounds`; `nil` when the optimization fails."
  [optimizer objective bounds start]
  (try
    (let [start (clamp-to-bounds bounds start)
          res (if (= :lbfgsb optimizer)
                (opt/minimize optimizer objective {:bounds bounds :initial start :stats? true})
                (opt/minimize optimizer objective {:bounds (simplex-bounds start) :length 1.0 :initial start :stats? true}))]
      (when (m/valid-double? (:value res)) res))
    (catch Exception _ nil)))

(defn- information-criteria
  [^double ll ^long k ^long n]
  (let [aic (m/- (m/* 2.0 k) (m/* 2.0 ll))
        dof (m/- n k 1)]
    {:aic aic
     :aicc (if (m/pos? dof) (m/+ aic (m// (m/* 2.0 k (m/inc k)) dof)) ##Inf)
     :bic (m/- (m/* k (m/log n)) (m/* 2.0 ll))}))

(defn- fit-statistics
  [d ^doubles xs ^long k]
  (let [n (alength xs)
        ll (log-likelihood d xs)]
    (merge {:k k :n n :loglik ll}
           (information-criteria ll k n)
           (when (r/continuous? d)
             {:ks (stats/ks-statistic xs d {:sorted? true})
              :cvm (stats/cvm-statistic xs d {:sorted? true})
              :ad (stats/ad-statistic xs d {:sorted? true})}))))



(defn- build-result
  [dist-key spec fixed method params tier converged? iterations xs k]
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
     :stats (fit-statistics d xs k)}))

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
  [{:keys [dist-key spec xs method optimizer restarts closed-form? rng opts]} fixed shown-fixed k]
  (let [estimators (when closed-form? (:estimators spec))
        closed-form (get estimators method)
        free (vec (remove #(contains? fixed (:name %)) (:params spec)))
        free-names (mapv :name free)
        objective (method-objective method xs spec (count free) opts)
        value-of (fn [result] (double (objective (:distribution result) (:params result))))
        estimated (estimate closed-form xs fixed (:exact? closed-form))]
    (cond
      (and (:exact? closed-form) estimated)
      (let [result (build-result dist-key spec shown-fixed method (merge estimated fixed) :closed-form true 0 xs k)]
        {:result result :value (value-of result)})
      ;; every parameter is fixed (the searched integer was the last free one)
      (empty? free)
      (let [result (build-result dist-key spec shown-fixed method fixed :closed-form true 0 xs k)]
        {:result result :value (value-of result)})
      :else
      (let [bounds (mapv (comp domain-bounds :domain) free)
            scale (robust-scale xs)
            objective-fn (objective-function dist-key fixed free-names objective)
            seed (or estimated
                     (when-let [initial (:initial spec)] (initial xs))
                     (estimate (:mle estimators) xs fixed false))
            starts (if (and seed (every? #(contains? seed %) free-names))
                     [(mapv seed free-names)]
                     (scanned-starts dist-key fixed free xs scale restarts rng))
            best (->> starts
                      (keep #(polish optimizer objective-fn bounds %))
                      (filter #(m/valid-double? (double (:value %))))
                      (sort-by :value)
                      (first))]
        (when (or (nil? best) (m/inf? (:value best)))
          (if (#{:mme :lmme} method)
            (undefined-moments! dist-key "moments to match")
            (throw (ex-info "fit: no usable estimate was found" {:type :not-converged :function 'fit :distribution dist-key}))))
        (let [result (build-result dist-key spec shown-fixed method (merge fixed (zipmap free-names (:point best)))
                                   :numeric (contains? #{nil :converged :stalled} (:status best)) (or (:iterations best) 0) xs k)]
          {:result result :value (:value best)})))))

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
    - `:restarts` (long, default `3`): The number of starting points taken from the scan, when there are no other starting values.
    - `:closed-form?` (boolean, default `true`): Set to `false` to skip the closed form estimators and fit numerically.
    - `:rng`: The random number generator of the scan. The default is a new unseeded generator; pass a seeded one for a reproducible scan. The scan only provides starting points, so the fit is the same up to the tolerance of the optimizer. Fits which start from a closed form estimate or from the starting values of the specification don't use it.
    - `:quantiles` (number or sequence, default `50`) and `:strategy` (keyword, default `:r7`): For `:qme`. Either the probabilities of the matched quantiles (numbers inside `(0, 1)`, at least as many as the fitted parameters) or their count, spread evenly from `0.01` to `0.99`; and the quantile estimation strategy of [[fastmath.stats/quantiles]].

  Returns a map with the keys:

  - `:distribution`: The fitted distribution object.
  - `:params`: The map of all parameters (fitted, searched and fixed) to doubles.
  - `:method`: The estimation method.
  - `:tier`: `:closed-form` or `:numeric`. A fit in which the search of an integer parameter leaves no parameter to optimize (the chi distribution) is `:closed-form`. A fit in which the search of an integer parameter leaves no parameter to optimize (the chi distribution) is `:closed-form`.
  - `:converged?`: `true` when the optimizer reported convergence or a stalled objective (it stopped changing, the usual end at the optimum), `false` when it stopped at the iteration limit. A closed form result is always converged.
  - `:iterations`: The number of iterations of the optimizer, `0` for a closed form.
  - `:fixed`: The map of fixed parameters, as given.
  - `:stats`: A map with `:k` (the number of fitted parameters, a searched integer included), `:n`, `:loglik`, `:aic`, `:aicc`, `:bic`, and for continuous distributions the goodness-of-fit statistics `:ks`, `:cvm` and `:ad`, which are descriptive, not tests. They describe the fitted distribution, whatever the method. `:aicc` is `##Inf` when `n` is not above `k + 1`.

  Moment matching with the moments above the second order, and L-moment matching, integrate the quantile function numerically unless the specification has its own `:moments` or `:lmoments`.

  Throws `ex-info` with `:type :not-fittable` when the distribution has no specification, with `:type :invalid-input` for data which fail the validation, an unknown `:method`, wrong `:quantiles`, a wrong `:fixed` value, an integer parameter which has to be fixed or more than one integer parameter to search, with `:type :moments-undefined` when `:mme` or `:lmme` finds no parameters for which the moments exist (the Cauchy distribution has none), and with `:type :not-converged` when no usable estimate is found. An evaluation of the objective which throws or is not finite counts as an infinite objective value.

  See also [[fastmath.random/fit-spec]], [[fastmath.random/log-likelihood]]."
  ([dist-key data] (fit dist-key data nil))
  ([dist-key data {:keys [method fixed optimizer restarts closed-form? rng]
                   :or {method :mle restarts 3 closed-form? true}
                   :as opts}]
   (when-not (some #{method} methods-list)
     (invalid-input! :method "unknown method" method (set methods-list)))
   (when-not (contains? #{nil :lbfgsb :nelder-mead} optimizer)
     (invalid-input! :optimizer "unsupported optimizer" optimizer #{:lbfgsb :nelder-mead}))
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
                  :restarts restarts :closed-form? closed-form? :rng (or rng (r/rng :jdk)) :opts opts}
         k (m/long-sub (count (:params spec)) (count fixed))
         closed-form (when closed-form? (get-in spec [:estimators method]))
         exact? (and (:exact? closed-form) (estimate closed-form xs fixed true))
         searched (remove #(contains? fixed (:name %)) (filter #(= :integer (:role %)) (:params spec)))]
     (cond
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
               (:result (fit-at best)))))))
