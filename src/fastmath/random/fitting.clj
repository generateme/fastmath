(ns fastmath.random.fitting
  "Fitting of probability distributions to data.

  A distribution is fitted by its key (see [[fastmath.random/distribution]]) to a sample. What the engine needs to know about the distribution comes from [[fastmath.random/fit-spec]]: the fitted parameters with their domains and roles, the data validation and, optionally, closed form estimators and starting values. Distributions without a specification can't be fitted.

  The maximum likelihood estimate (`:mle`) is found with a closed form estimator when the specification has an exact one, otherwise numerically: the box constrained quasi-Newton method L-BFGS-B started from the starting values of the specification or, when there are none, from the best points of a scan of the parameter space with a quantile matching objective.

  The result of [[fit]] is a map with the fitted distribution, its parameters and statistics of the fit (information criteria and, for continuous distributions, goodness-of-fit statistics)."
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.stats :as stats]
            [fastmath.optimization :as opt])
  (:import [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^{:private true :const true :tag 'double} lower-limit 1.0e-10)

(def ^:private seed-probabilities (vec (range 0.025 0.9751 0.05)))

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
  (let [[^double q1 ^double q3] (stats/quantiles xs [0.25 0.75])
        iqr-scale (m// (m/- q3 q1) 1.349)]
    (cond
      (m/pos? iqr-scale) iqr-scale
      (m/pos? (double (stats/stddev xs))) (double (stats/stddev xs))
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

(defn- negative-log-likelihood
  "Function of a point (values of the free parameters) returning minus the log likelihood, `##Inf` when it is not finite or the distribution can't be created."
  [dist-key fixed names ^doubles xs]
  (fn ^double [point]
    (try
      (let [ll (log-likelihood (r/distribution dist-key (merge fixed (zipmap names point))) xs)]
        (if (m/valid-double? ll) (m/- ll) ##Inf))
      (catch Exception _ ##Inf))))

(defn- quantile-objective
  "Function of a point returning the squared distance between the quantiles of the data and the quantiles of the distribution, in units of `scale`."
  [dist-key fixed names scale data-quantiles]
  (fn ^double [point]
    (try
      (let [^double scale scale
            d (r/distribution dist-key (merge fixed (zipmap names point)))
            sum (double (reduce m/+ 0.0 (map (fn [^double p ^double q] (m/sq (m/- (double (r/icdf d p)) q)))
                                             seed-probabilities data-quantiles)))
            res (m// sum (m/* scale scale (count seed-probabilities)))]
        (if (m/valid-double? res) res ##Inf))
      (catch Exception _ ##Inf))))

(defn- clamp-to-bounds
  [bounds point]
  (mapv (fn [[^double lo ^double hi] ^double v] (m/constrain v lo hi)) bounds point))

(defn- scanned-starts
  "Starting points: the best results of a scan of the parameter space with the quantile matching objective."
  [dist-key fixed free xs scale restarts rng]
  (let [names (mapv :name free)
        data-quantiles (stats/quantiles xs seed-probabilities)
        bounds (mapv #(seed-range % xs scale) free)
        res (opt/scan-and-minimize :lbfgsb (quantile-objective dist-key fixed names scale data-quantiles)
                                   {:bounds bounds :N (m/* 50 (count free)) :n restarts :take-last-n restarts
                                    :rng rng :stats? true})]
    (map :point (if (map? res) [res] res))))

(defn- polish
  "Minimizes minus the log likelihood from `start`; `nil` when the optimization fails."
  [optimizer nll bounds start]
  (try
    (let [res (opt/minimize optimizer nll {:bounds bounds :initial (clamp-to-bounds bounds start) :stats? true})]
      (when (m/valid-double? (double (:value res))) res))
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
  [dist-key spec fixed params tier converged? iterations xs]
  (let [d (try (r/distribution dist-key params)
               (catch Exception e
                 (throw (ex-info "fit: the fitted parameters don't make a distribution"
                                 {:type :not-converged :function 'fit :params params} e))))]
    {:distribution d
     :params (into {} (map (fn [{:keys [name]}] [name (double (get params name))])) (:params spec))
     :method :mle
     :tier tier
     :converged? converged?
     :iterations iterations
     :fixed fixed
     :stats (fit-statistics d xs (m/long-sub (count (:params spec)) (count fixed)))}))

(defn fit
  "Fits a distribution to data with the maximum likelihood method.

  The distribution is given by its key (see [[fastmath.random/distribution]]), the fitted parameters are the ones listed by [[fastmath.random/fit-spec]]. When the specification of the distribution has an exact closed form estimator (and `:fixed` is empty and `:closed-form?` is not `false`) it is the result. Otherwise minus the log likelihood is minimized, from the starting values of the specification, or from the best points of a scan of the parameter space (a quantile matching objective over ranges derived from the data and the roles of the parameters) when it has none.

  Parameters:

  - `dist-key` (keyword): A key of [[fastmath.random/distribution]] which has a specification.
  - `data` (sequence of numbers): The sample; it has to pass the validation of the specification (for example `:gamma` requires positive values).
  - `opts` (map, optional):
    - `:method` (keyword, default `:mle`): The estimation method, `:mle` is the only one for now.
    - `:fixed` (map, default `nil`): Parameters which are not fitted, as a map from the name of a parameter to its value. At least one parameter has to stay free. A parameter with the `:integer` role has to be given here.
    - `:optimizer` (keyword, default `:lbfgsb`): The method of [[fastmath.optimization]] used for the numerical maximization. It has to accept infinite bounds, which only `:lbfgsb` does.
    - `:restarts` (long, default `3`): The number of starting points taken from the scan, when the specification has no starting values.
    - `:closed-form?` (boolean, default `true`): Set to `false` to skip the closed form estimator and fit numerically.
    - `:rng`: The random number generator of the scan. The default is a generator with a fixed seed, so the fit is deterministic.

  Returns a map with the keys:

  - `:distribution`: The fitted distribution object.
  - `:params`: The map of all parameters (fitted and fixed) to doubles.
  - `:method`: The estimation method.
  - `:tier`: `:closed-form` or `:numeric`.
  - `:converged?`: `true` when the optimizer reported convergence or a stalled objective (it stopped changing, the usual end at the optimum), `false` when it stopped at the iteration limit. A closed form result is always converged.
  - `:iterations`: The number of iterations of the optimizer, `0` for a closed form.
  - `:fixed`: The map of fixed parameters.
  - `:stats`: A map with `:k` (the number of fitted parameters), `:n`, `:loglik`, `:aic`, `:aicc`, `:bic`, and for continuous distributions the goodness-of-fit statistics `:ks`, `:cvm` and `:ad` (see [[fastmath.stats/ks-statistic]], [[fastmath.stats/cvm-statistic]], [[fastmath.stats/ad-statistic]]), which are descriptive, not tests. `:aicc` is `##Inf` when `n` is not above `k + 1`.

  Throws `ex-info` with `:type :not-fittable` when the distribution has no specification, with `:type :invalid-input` for data which fail the validation, an unknown `:method`, a wrong `:fixed` value or a missing `:integer` parameter, and with `:type :not-converged` when no usable estimate is found. An evaluation of the likelihood which throws or is not finite counts as an infinite objective value.

  See also [[fastmath.random/fit-spec]], [[fastmath.random/log-likelihood]]."
  ([dist-key data] (fit dist-key data nil))
  ([dist-key data {:keys [method fixed optimizer restarts closed-form? rng]
                   :or {method :mle optimizer :lbfgsb restarts 3 closed-form? true}}]
   (when-not (= :mle method)
     (invalid-input! :method "unknown method" method #{:mle}))
   (let [spec (or (r/fit-spec dist-key)
                  (throw (ex-info (str "fit: the distribution " dist-key " can't be fitted")
                                  {:type :not-fittable :function 'fit :distribution dist-key})))
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
         _ (doseq [{:keys [name role]} (:params spec)]
             (when (and (= :integer role) (not (contains? fixed name)))
               (invalid-input! :fixed (str name " is an integer parameter and has to be fixed") name "a value in :fixed")))
         free (vec (remove #(contains? fixed (:name %)) (:params spec)))
         _ (when (empty? free)
             (invalid-input! :fixed "at least one parameter has to stay free" fixed "fewer fixed parameters"))
         xs (m/seq->double-array (sort data))
         closed-form (when (and closed-form? (empty? fixed)) (get-in spec [:estimators method]))]
     (if (:exact? closed-form)
       (build-result dist-key spec fixed ((:estimate closed-form) xs) :closed-form true 0 xs)
       (let [free-names (mapv :name free)
             bounds (mapv (comp domain-bounds :domain) free)
             scale (robust-scale xs)
             rng (or rng (r/rng :jdk 1))
             nll (negative-log-likelihood dist-key fixed free-names xs)
             seed (cond
                    closed-form ((:estimate closed-form) xs)
                    (:initial spec) ((:initial spec) xs))
             starts (if seed
                      [(mapv seed free-names)]
                      (scanned-starts dist-key fixed free xs scale restarts rng))
             best (->> starts
                       (keep #(polish optimizer nll bounds %))
                       (sort-by :value)
                       (first))]
         (when-not best
           (throw (ex-info "fit: no usable estimate was found" {:type :not-converged :function 'fit :distribution dist-key})))
         (build-result dist-key spec fixed (merge fixed (zipmap free-names (:point best)))
                       :numeric (contains? #{nil :converged :stalled} (:status best)) (or (:iterations best) 0) xs))))))
