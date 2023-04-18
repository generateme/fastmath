(ns fastmath.stats.bootstrap
  "Bootstrap methods and confidence intervals"
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.vector :as v]
            [fastmath.stats :as stats]))

(defn- bootstrap-antithetic
  [rng model method samples size]
  (->> (repeatedly #(let [s (r/->seq rng size method)]
                      [s (map (fn [^double v] (- 1.0 v)) s)]))
       (mapcat identity)
       (take samples)
       (map (fn [xs] (map (partial r/icdf model) xs)))))

(defn- build-model
  [data {:keys [rng smoothing kernel bandwidth distribution]
         :or {distribution :real-discrete-distribution}}]
  (if (= :kde smoothing)
    (r/distribution :continuous-distribution
                    {:data data
                     :kde kernel
                     :bandwidth bandwidth
                     :rng rng})
    (r/distribution distribution {:data data
                                  :rng rng})))

(defn bootstrap-stats
  "Calculate bootstrap analysis.

  Arguments:

  * a map containing `:data` and `:samples`
  * `statistic` - bootstrap statistic

  Returns a map containing:

  * `:data`, `:samples`, `model` and `:statistic` - from input
  * `:t0` - statistic from data (single value)
  * `:ts` - statistic from bootstrap samples (sequence)
  * `:bias` - difference between mean of ts and t0
  * `:mean`, `:median`, `:variance`, `:stddev`, `:sem` - statistics from ts"
  [{:keys [data samples] :as input} statistic]
  (let [t0 (statistic data)
        ts (map statistic samples)
        ats (double-array ts)
        m (stats/mean ats)
        variance (stats/variance ats)
        stddev (m/sqrt variance)
        bias (- m t0)]
    (assoc input 
           :statistic statistic
           :t0 t0
           :ts ts
           :bias bias
           :mean m
           :median (stats/median ats)
           :sem (/ stddev (m/sqrt (count data)))
           :stddev stddev
           :variance variance)))

(defn jackknife
  "Generates set of samples using jackknife leave-one-out method"
  [vs]
  (map (fn [id]
         (let [[a b] (split-at id vs)]
           (concat a (rest b)))) (range (count vs))))

(defn jackknife+
  "Generates set of samples using jackknife positive method"
  [vs]
  (map (fn [id]
         (conj vs (nth vs id))) (range (count vs))))


(defn bootstrap
  "Create set of samples from given data (nonparametric) or model (parametric).

  Input:
  * sequence of values (any type) or sequence of sequences for multidimensional data
  * a map containing:
      * `:data` - sequence
      * `:model` - model for parametric bootstrap (optional)
  * `statistic` function which returns statistic value (optional)

  Parameters:
  * `:samples` - number of bootstrapped samples (default: 500)
  * `:size` - forced size of individual sample (default: same as source)
  * `:method`
      *  `nil` (default) - random
      * `:jackknife` for leave-one-out jackknife
      * `:jackknife+` for positive jackknife
      * any accepted by `r/->seq` method
  * `:rng` - random number generator (see: `r/rng`)
  * `:smoothing` - smoothing bootstrap:
      * `:kde` - kernel density estimation, additional options are: `:kernel` (default) and `:bandwidth` (auto)
      * `:gaussian` - add random value from N(0,standard error)
  * `:distribution` - distribution used to auto generate model (distribution) from data:
      * `:real-discrete-distribution` - default
      * `:integer-discrete-distribution` - for integer values
      * `:categorical-distribution` - for any other type
  * `:dimensions` - if set to `:multi` - multidimensional data and models are created

  As a model it can be:
  * any distribution object
  * any 0-arity function which returns random sample

  When model is ommited, function creates discrete distribution.
  When multidimensional data are provided, models should be created for every dimension (as a sequence)."
  ([input] (bootstrap input nil))
  ([input statistic] (bootstrap input statistic {}))
  ([input statistic {:keys [rng ^long samples ^long size method antithetic?
                            smoothing dimensions]
                     :or {samples 500}
                     :as params}]
   (let [{:keys [data model] :as input} (if (map? input) input {:data input})]
     (if-not (or model (#{:jackknife :jackknife+} method))
       (bootstrap (assoc input :model (if (= :multi dimensions)
                                        (map #(build-model % params) data)
                                        (build-model data params))) statistic params)
       (let [xsss (case method
                    :jackknife (jackknife data)
                    :jackknife+ (jackknife+ data)
                    (let [rng (or rng (r/rng :jvm))
                          data0 data
                          data (if (not= :multi dimensions) [data] data)
                          models (if (not= :multi dimensions) [model] model)
                          sizes (if size (repeat (count data) size) (map count data))
                          gen (fn [model size]
                                (if (r/distribution? model)
                                  (if antithetic?
                                    (bootstrap-antithetic rng model method samples size)
                                    (repeatedly samples #(r/->seq model size method)))
                                  (repeatedly samples #(repeatedly size model))))
                          xsss (map gen models sizes)
                          xsss (if (= :gaussian smoothing)
                                 (let [sds (map #(m/sqrt (/ (stats/variance %1) %2)) data sizes)]
                                   (map (fn [xss sd size]
                                          (map (fn [xs]                                                 
                                                 (v/add xs (repeatedly size #(r/grandom rng sd)))) xss))
                                        xsss sds sizes))
                                 xsss)
                          xsss (if (not= :multi dimensions)
                                 (first xsss)
                                 (apply map vector xsss))]
                      (conj xsss data0)))
             res (assoc input :samples xsss)]
         (if statistic
           (bootstrap-stats res statistic)
           res))))))

(defn ci-normal
  "Normal (gaussian) bias-corrected confidence interval

  `:t0` and `:ts` are obligatory"
  ([boot-data] (ci-normal boot-data 0.05))
  ([{:keys [^double t0 ts stddev bias]} ^double alpha]
   (let [a (- 1.0 (/ alpha 2.0))
         ats (delay (m/seq->double-array ts))
         ^double stddev (or stddev (stats/stddev @ats))
         ^double bias (or bias (- (stats/mean @ats) t0))
         tb (- t0 bias)
         merr (* stddev (r/icdf r/default-normal a))]
     [(- tb merr) (+ tb merr) t0])))

(defn ci-basic
  "Basic percentile confidence interval

  `:t0` and `:ts` are obligatory"
  ([boot-data] (ci-basic boot-data 0.05))
  ([boot-data ^double alpha] (ci-basic boot-data alpha :legacy))
  ([{:keys [^double t0 ts]} ^double alpha estimation-strategy]
   (let [a (/ alpha 2.0)
         t2 (* 2.0 t0)
         [^double q1 ^double q2] (stats/quantiles ts [(- 1.0 a) a] estimation-strategy)]
     [(- t2 q1) (- t2 q2) t0])))

(defn ci-percentile
  "Percentile confidence interval

  `:t0` and `:ts` are obligatory"
  ([boot-data] (ci-percentile boot-data 0.05))
  ([boot-data ^double alpha] (ci-percentile boot-data alpha :legacy))
  ([{:keys [^double t0 ts]} ^double alpha estimation-strategy]
   (let [a (/ alpha 2.0)]
     (conj (stats/quantiles ts [a (- 1.0 a)] estimation-strategy) t0))))

(defn- cdf-accelerated-quantile
  ^double [^double z0 ^double z ^double a]
  (let [num (+ z0 z)
        denom (if (zero? a) 1.0 (- 1.0 (* a num)))]
    (->> (+ z0 (/ num denom))
         (r/cdf r/default-normal))))

(defn- empirical-cdf
  ^double [vs ^double value]
  (/ (double (count (filter (fn [^double v] (< v value)) vs))) (count vs)))

(defn- bca-common
  ([ats q1 q2 t0 estimation-strategy]
   (bca-common ats q1 q2 t0 0.0 estimation-strategy))
  ([ats q1 q2 t0 accel estimation-strategy]
   (let [e (empirical-cdf ats t0)
         ^double z0 (r/icdf r/default-normal e)
         ^double z1 (r/icdf r/default-normal q1)
         ^double z2 (r/icdf r/default-normal q2)
         q1 (cdf-accelerated-quantile z0 z1 accel)
         q2 (cdf-accelerated-quantile z0 z2 accel)]
     (conj (stats/quantiles ats [q1 q2] estimation-strategy) t0))))

(defn ci-bc
  "Bias-corrected confidence interval"
  ([boot-data] (ci-bc boot-data 0.05))
  ([boot-data ^double alpha] (ci-bc boot-data alpha :legacy))
  ([{:keys [^double t0 ts]} ^double alpha estimation-strategy]
   (let [ats (m/seq->double-array ts)
         a (/ alpha 2.0)]
     (bca-common ats a (- 1.0 a) t0 estimation-strategy))))

(defn- acceleration [ts] (/ (stats/skewness ts :skew) -6.0))

(defn ci-bca
  "Bias-corrected and accelerated confidence interval.

  There are two ways to calculate acceleration:
  * jackknife method (when boot-data contains `:data` and `:statistic`)
  * empirical from bootstrap estimations `ts` otherwise"
  ([boot-data] (ci-bca boot-data 0.05))
  ([boot-data ^double alpha] (ci-bca boot-data alpha :legacy))
  ([{:keys [^double t0 ts data statistic]} ^double alpha estimation-strategy]
   (let [ats (m/seq->double-array ts)
         a (/ alpha 2.0)
         acc (if-not (and data statistic)
               (acceleration ats)
               (acceleration (->> data
                                  (jackknife)
                                  (map statistic))))]
     (bca-common ats a (- 1.0 a) t0 acc estimation-strategy))))

(defn ci-studentized
  "Confidence interval from studentized data.

  `:t0`, `:ts`, `:data` and `:samples` are obligatory in `boot-data`"
  ([boot-data] (ci-studentized boot-data 0.05))
  ([boot-data ^double alpha] (ci-studentized boot-data alpha :legacy))
  ([{:keys [^double t0 ts data samples]} ^double alpha estimation-strategy]
   (assert (and (seq data)
                (seq samples)) "Bootstrap samples can't be empty.")
   (let [a (/ alpha 2.0)
         z (map (fn [^double t s]
                  (/ (- t t0) (stats/stddev s))) ts samples)
         stddev (stats/stddev data)
         [^double q1 ^double q2] (stats/quantiles z [(- 1.0 a) a] estimation-strategy)]
     [(- t0 (* stddev q1))
      (- t0 (* stddev q2))
      t0])))

(defn ci-t
  "Student's T confidence interval.

  `:t0` and `:ts` are obligatory"
  ([boot-data] (ci-t boot-data 0.05))
  ([{:keys [^double t0 ts stddev]} ^double alpha]
   (let [a (/ alpha 2.0)
         ^double stddev (or stddev (stats/stddev ts))
         t (* stddev (r/icdf (r/distribution :t {:degrees-of-freedom (dec (count ts))}) (- 1.0 a)))]
     [(- t0 t) (+ t0 t) t0])))
