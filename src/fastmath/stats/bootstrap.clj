(ns fastmath.stats.bootstrap
  "Bootstrap methods and confidence intervals"
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.stats :as stats]
            [fastmath.protocols :as prot]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

(defn- bootstrap-antithetic
  [rng model method samples size]
  (->> (repeatedly #(let [s (r/->seq rng size method)]
                      [s (map (fn [^double v] (- 1.0 v)) s)]))
       (mapcat identity)
       (take samples)
       (map (fn [xs] (map (partial r/icdf model) xs)))))

(defn- build-model
  [data {:keys [rng smoothing kernel bandwidth distribution size]
         :or {distribution :real-discrete-distribution
              kernel :gaussian}}]
  (case smoothing
    :kde (r/distribution :continuous-distribution
                         {:data data
                          :kde kernel
                          :bandwidth bandwidth
                          :rng rng})
    :gaussian (let [sd (m/sqrt (/ (stats/variance data) (double (or size (count data)))))
                    d (r/distribution distribution {:data data :rng rng})]
                (reify
                  prot/DistributionIdProto
                  (distribution? [_] true)
                  prot/DistributionProto
                  (icdf [_ p] (m/+ (r/grandom rng sd) (double (prot/icdf d p))))
                  prot/RNGProto
                  (->seq [_] (repeatedly #(m/+ (r/grandom rng sd) (double (prot/sample d)))))
                  (->seq [_  n] (repeatedly n #(m/+ (r/grandom rng sd) (double (prot/sample d)))))))
    (if (every? number? data)
      (r/distribution distribution {:data data :rng rng})
      (r/distribution :categorical-distribution {:data data :rng rng}))))

(defn bootstrap-stats
  "Calculates summary statistics for bootstrap results.

  Takes bootstrap output (typically from [[bootstrap]]) and a statistic function,
  computes the statistic on the original data (`t0`) and on each bootstrap sample (`ts`),
  and derives various descriptive statistics from the distribution of `ts`.

  Parameters:

  * `boot-data` (map): A map containing:
      * `:data`: The original dataset.
      * `:samples`: A collection of bootstrap samples (e.g., from [[bootstrap]]).
      * (optional) other keys like `:model` from bootstrap generation.
  * `statistic` (function): A function that accepts a sequence of data and returns
    a single numerical statistic (e.g., `fastmath.stats/mean`, `fastmath.stats/median`).

  Returns a map which is the input `boot-data` augmented with bootstrap analysis results:

  * `:statistic`: The statistic function applied.
  * `:t0`: The statistic calculated on the original `:data`.
  * `:ts`: A sequence of the statistic calculated on each bootstrap sample in `:samples`.
  * `:bias`: The estimated bias of the statistic: `mean(:ts) - :t0`.
  * `:mean`, `:median`, `:variance`, `:stddev`, `:sem`: Descriptive statistics
    (mean, median, variance, standard deviation, standard error of the mean)
    calculated from the distribution of `:ts`.

  This function prepares the results for calculating various bootstrap
  confidence intervals (e.g., [[ci-normal]], [[ci-percentile]], etc.)."
  [{:keys [data samples] :as input} statistic]
  (let [^double t0 (statistic data)
        ts (map statistic samples)
        ats (double-array (filter m/valid-double? ts))
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
  "Generates a set of samples from a given sequence using the jackknife leave-one-out method.

  For an input sequence `vs` of size `n`, this method creates `n` samples. Each sample is formed by removing a single observation from the original sequence.

  Parameters:

  * `vs` (sequence): The input data sequence.

  Returns a sequence of sequences. The i-th inner sequence is `vs` with the i-th element removed.

  These samples are commonly used for estimating the bias and standard error of a statistic (e.g., via [[bootstrap-stats]])."
  [vs]
  (map (fn [id]
         (let [[a b] (split-at id vs)]
           (concat a (rest b)))) (range (count vs))))

(defn jackknife+
  "Generates a set of samples from a sequence using the 'jackknife positive' method.

  For an input sequence `vs` of size `n`, this method creates `n` samples. Each sample is formed by duplicating a single observation from the original sequence and adding it back to the original sequence. Thus, each sample has size `n+1`.

  Parameters:

  * `vs` (sequence): The input data sequence.

  Returns a sequence of sequences. The i-th inner sequence is `vs` with an additional copy of the i-th element of `vs`.

  This method is used in specific resampling techniques for estimating bias and variance of a statistic."
  [vs]
  (map (fn [id]
         (conj vs (nth vs id))) (range (count vs))))

(defn bootstrap
  "Generates bootstrap samples from a given dataset or probabilistic model for resampling purposes.

  This function supports both **nonparametric bootstrap** (resampling directly from the data)
  and **parametric bootstrap** (resampling from a statistical model estimated from or
  provided for the data). It can optionally apply a statistic function to the original
  data and each sample, returning summary statistics for the bootstrap distribution.

  The primary input can be:

  *   A sequence of data values (for nonparametric bootstrap).
  *   A map containing:
      *   `:data`: The sequence of data values.
      *   `:model`: An optional model for parametric bootstrap. If not provided,
          a default discrete distribution is built from the data (see `:distribution`
          and `:smoothing`).

  The function offers various parameters to control the sampling process and
  model generation.

  Parameters:

  *   `input` (sequence or map): The data source. Can be a sequence of numbers
      or a map containing `:data` and optionally `:model`. Can be sequence of
      sequences for multidimensional data (when `:dimensions` is `:multi`).
  *   `statistic` (function, optional): A function that takes a sequence of data
      and returns a single numerical value (e.g., `fastmath.stats/mean`,
      `fastmath.stats/median`). If provided, `bootstrap-stats` is called on the
      results.
  *   `params` (map, optional): An options map to configure the bootstrap process.
      Keys include:
      *   `:samples` (long, default: 500): The number of bootstrap samples to generate.
      *   `:size` (long, optional): The size of each individual bootstrap sample.
          Defaults to the size of the original data.
      *   `:method` (keyword, optional): Specifies the sampling method.
          *   `nil` (default): Standard random sampling with replacement.
          *   `:jackknife`: Performs leave-one-out jackknife resampling (ignores
              `:samples` and `:size`).
          *   `:jackknife+`: Performs positive jackknife resampling (duplicates
              each observation once; ignores `:samples`).
          *   Other keywords are passed to `fastmath.random/->seq` for sampling
              from a distribution (only relevant if a `:model` is used or built).
      *   `:rng` (random number generator, optional): An instance of a random number
          generator (see `fastmath.random/rng`). A default JVM RNG is used if not provided.
      *   `:smoothing` (keyword, optional): Applies smoothing to the bootstrap process.
          *   `:kde`: Uses Kernel Density Estimation to smooth the empirical distribution
              before sampling. Requires specifying `:kernel` (default `:gaussian`)
              and optionally `:bandwidth` (auto-estimated by default).
          *   `:gaussian`: Adds random noise drawn from N(0, standard error) to each
              resampled value.
      *   `:distribution` (keyword, default: `:real-discrete-distribution`): The type
          of discrete distribution to build automatically from the data if no explicit
          `:model` is provided. Other options include `:integer-discrete-distribution`
          (for integer data) and `:categorical-distribution` (for any data type).
      *   `:dimensions` (keyword, optional): If set to `:multi`, treats the input
          `:data` as a sequence of sequences (multidimensional data). Models are
          built or used separately for each dimension, and samples are generated
          as sequences of vectors.
      *   `:antithetic?` (boolean, default: `false`): If `true`, uses antithetic sampling
          for variance reduction (paired samples are generated as `x` and `1-x` from a uniform
          distribution, then transformed by the inverse CDF of the model). Requires sampling
          from a distribution model.
      *   `:include?` (boolean, default: `false`): If `true`, the original dataset
          is included as one of the samples in the output collection.

  Model for parametric bootstrap:

  The `:model` parameter in the input map can be:

  *   Any `fastmath.random` distribution object (e.g., `(r/distribution :normal {:mu 0 :sd 1})`).
  *   Any 0-arity function that returns a random sample when called.

  If `:model` is omitted from the input map, a default discrete distribution
  (`:real-discrete-distribution` by default, see `:distribution` param) is built
  from the `:data`. Smoothing options (`:smoothing`) apply to this automatically
  built model.

  When `:dimensions` is `:multi`, `:model` should be a sequence of models, one for
  each dimension.

  Returns:

  *   If `statistic` is provided: A map containing the original input map augmented
      with analysis results from `bootstrap-stats` (e.g., `:t0`, `:ts`, `:bias`,
      `:mean`, `:stddev`).
  *   If `statistic` is `nil`: A map containing the original input map augmented
      with the generated bootstrap samples in the `:samples` key. The `:samples`
      value is a collection of sequences, where each inner sequence is one
      bootstrap sample. If `:dimensions` is `:multi`, samples are sequences of vectors.

  See also [[jackknife]], [[jackknife+]], [[bootstrap-stats]],
  [[ci-normal]], [[ci-basic]], [[ci-percentile]], [[ci-bc]], [[ci-bca]],
  [[ci-studentized]], [[ci-t]]."
  ([input] (bootstrap input nil))
  ([input params-or-statistic]
   (if (fn? params-or-statistic)
     (bootstrap input params-or-statistic nil)
     (bootstrap input nil params-or-statistic)))
  ([input statistic {:keys [rng ^long samples size method antithetic? dimensions include? multi?]
                     :or {samples 500}
                     :as params}]
   (let [{:keys [data model] :as input} (if (map? input) input {:data input})
         params (assoc params :rng (or rng (r/rng :jvm)))
         multi? (or multi? (= :multi dimensions))] ;; ensure we have rng before model creation
     (if-not (or model (#{:jackknife :jackknife+} method))
       (bootstrap (assoc input :model (if multi?
                                        (map #(build-model % params) data)
                                        (build-model data params))) statistic params)
       (let [xsss (case method
                    :jackknife (jackknife data)
                    :jackknife+ (jackknife+ data)
                    (let [data0 data
                          data (if multi? data [data])
                          models (if multi? model [model])
                          sizes (if size (repeat (count data) size) (map count data))
                          gen (fn [model size]
                                (if (r/distribution? model)
                                  (if antithetic?
                                    (bootstrap-antithetic rng model method samples size)
                                    (repeatedly samples #(r/->seq model size method)))
                                  (repeatedly samples #(repeatedly size model))))]
                      (as-> (map gen models sizes) xsss                        
                        (if multi?
                          (apply map vector xsss)
                          (first xsss))
                        (if include?
                          (conj xsss data0)
                          xsss))))
             res (assoc input :samples xsss)]
         (if statistic
           (bootstrap-stats res statistic)
           res))))))

(defn ci-normal
  "Calculates a Normal (Gaussian) approximation bias-corrected confidence interval.

  This method assumes the distribution of the bootstrap replicates of the statistic (`:ts`)
  is approximately normal. It computes a confidence interval centered around the
  mean of the bootstrap statistics, adjusted by the estimated bias (`mean(:ts) - :t0`),
  and uses the standard error of the bootstrap statistics for scaling.

  Parameters:

  * `boot-data` (map): A map containing bootstrap results. Typically produced by [[bootstrap-stats]].
    Requires keys:
      * `:t0` (double): The statistic calculated on the original data.
      * `:ts` (sequence of numbers): The statistic calculated on each bootstrap sample.
    May optionally include pre-calculated `:stddev` (standard deviation of `:ts`)
    and `:bias` for efficiency.
  * `alpha` (double, optional): The significance level for the interval.
    Defaults to `0.05` (for a 95% CI). The interval is based on the `alpha/2`
    and `1 - alpha/2` quantiles of the standard normal distribution.

  Returns a vector `[lower-bound, upper-bound, t0]`.

  * `lower-bound` (double): The lower limit of the confidence interval.
  * `upper-bound` (double): The upper limit of the confidence interval.
  * `t0` (double): The statistic calculated on the original data (from `boot-data`).

  See also [[bootstrap-stats]] for input preparation and other confidence interval methods:
  [[ci-basic]], [[ci-percentile]], [[ci-bc]], [[ci-bca]], [[ci-studentized]], [[ci-t]]."
  ([boot-data] (ci-normal boot-data 0.05))
  ([{:keys [^double t0 ts stddev bias]} ^double alpha]
   (let [a (- 1.0 (/ alpha 2.0))
         ats (delay (m/seq->double-array ts))
         ^double stddev (or stddev (stats/stddev @ats))
         ^double bias (or bias (- (stats/mean @ats) t0))
         tb (- t0 bias)
         merr (* stddev ^double (r/icdf r/default-normal a))]
     [(- tb merr) (+ tb merr) t0])))

(defn ci-basic
  "Calculates the Basic (or Percentile-t) bootstrap confidence interval.

  This method is based on the assumption that the distribution of the bootstrap
  replicates (`:ts`) centered around the true statistic (`t`) is approximately the
  same as the distribution of the original statistic (`:t0`) centered around the mean
  of the bootstrap replicates (`mean(:ts)`).

  The interval is constructed using the quantiles of the bootstrap replicates (`:ts`)
  relative to the original statistic (`:t0`). Specifically, the lower bound is
  `2 * :t0 - q_upper` and the upper bound is `2 * :t0 - q_lower`, where `q_lower` and
  `q_upper` are the `alpha/2` and `1 - alpha/2` quantiles of `:ts`, respectively.

  Parameters:

  * `boot-data` (map): A map containing bootstrap results, typically from [[bootstrap-stats]].
    Requires keys:
      * `:t0` (double): The statistic calculated on the original data.
      * `:ts` (sequence of numbers): The statistic calculated on each bootstrap sample.
  * `alpha` (double, optional): The significance level for the interval.
    Defaults to `0.05` (for a 95% CI). The interval is based on the `alpha/2`
    and `1 - alpha/2` quantiles of the `:ts` distribution.
  * `estimation-strategy` (keyword, optional): Specifies the quantile estimation strategy
    used to calculate the quantiles of `:ts`. Defaults to `:legacy`. See [[quantiles]]
    for available options (e.g., `:r1` through `:r9`).

  Returns a vector `[lower-bound, upper-bound, t0]`.

  * `lower-bound` (double): The lower limit of the confidence interval.
  * `upper-bound` (double): The upper limit of the confidence interval.
  * `t0` (double): The statistic calculated on the original data (from `boot-data`).

  See also [[bootstrap-stats]] for input preparation and other confidence interval methods:
  [[ci-normal]], [[ci-percentile]], [[ci-bc]], [[ci-bca]], [[ci-studentized]], [[ci-t]], [[quantiles]]."
  ([boot-data] (ci-basic boot-data 0.05))
  ([boot-data ^double alpha] (ci-basic boot-data alpha :legacy))
  ([{:keys [^double t0 ts]} ^double alpha estimation-strategy]
   (let [a (/ alpha 2.0)
         t2 (* 2.0 t0)
         [^double q1 ^double q2] (stats/quantiles ts [(- 1.0 a) a] estimation-strategy)]
     [(- t2 q1) (- t2 q2) t0])))

(defn ci-percentile
  "Calculates the Percentile bootstrap confidence interval.

  This is the simplest bootstrap confidence interval method. It directly uses
  the quantiles of the bootstrap replicates of the statistic (`:ts`) as the
  confidence interval bounds.

  For a confidence level of `1 - alpha`, the interval is formed by taking the
  `alpha/2` and `1 - alpha/2` quantiles of the distribution of bootstrap
  replicates (`:ts`).

  Parameters:

  * `boot-data` (map): A map containing bootstrap results, typically from [[bootstrap-stats]].
    Requires keys:
      * `:t0` (double): The statistic calculated on the original data.
      * `:ts` (sequence of numbers): The statistic calculated on each bootstrap sample.
  * `alpha` (double, optional): The significance level for the interval.
    Defaults to `0.05` (for a 95% CI). The interval is based on the `alpha/2`
    and `1 - alpha/2` quantiles of the `:ts` distribution.
  * `estimation-strategy` (keyword, optional): Specifies the quantile estimation strategy
    used to calculate the quantiles of `:ts`. Defaults to `:legacy`. See [[quantiles]]
    for available options (e.g., `:r1` through `:r9`).

  Returns a vector `[lower-bound, upper-bound, t0]`.

  * `lower-bound` (double): The `alpha/2` quantile of `:ts`.
  * `upper-bound` (double): The `1 - alpha/2` quantile of `:ts`.
  * `t0` (double): The statistic calculated on the original data (from `boot-data`).

  See also [[bootstrap-stats]] for input preparation and other confidence interval methods:
  [[ci-normal]], [[ci-basic]], [[ci-bc]], [[ci-bca]], [[ci-studentized]], [[ci-t]], [[quantiles]]."
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
  "Calculates the Bias-Corrected (BC) bootstrap confidence interval.

  This method adjusts the standard Percentile bootstrap interval ([[ci-percentile]])
  to account for potential bias in the statistic's distribution. The correction
  is based on the proportion of bootstrap replicates of the statistic (`:ts`)
  that are less than the statistic calculated on the original data (`:t0`).

  The procedure involves:
  1.  Calculating a bias correction factor ($z_0$) based on the empirical cumulative
      distribution function (CDF) of the bootstrap replicates at the point of the
      original statistic ($z_0 = \\Phi^{-1}(\\text{Proportion of } t^* < t_0)$,
      where $\\Phi^{-1}$ is the inverse standard normal CDF).
  2.  Shifting the standard normal quantiles corresponding to the desired confidence
      level ($\\alpha/2$ and $1-\\alpha/2$) by $z_0$.
  3.  Finding the corresponding quantiles in the distribution of bootstrap
      replicates (`:ts`) based on these shifted probabilities.

  Parameters:

  * `boot-data` (map): A map containing bootstrap results, typically from [[bootstrap-stats]].
    Requires keys:
      * `:t0` (double): The statistic calculated on the original data.
      * `:ts` (sequence of numbers): The statistic calculated on each bootstrap sample.
  * `alpha` (double, optional): The significance level for the interval.
    Defaults to `0.05` (for a 95% CI). The interval is based on quantiles of the
    `:ts` distribution, adjusted by the bias correction factor.
  * `estimation-strategy` (keyword, optional): Specifies the quantile estimation strategy
    used to calculate the final interval bounds from `:ts` after applying corrections.
    Defaults to `:legacy`. See [[quantiles]] for available options (e.g., `:r1` through `:r9`).

  Returns a vector `[lower-bound, upper-bound, t0]`.

  * `lower-bound` (double): The lower limit of the confidence interval.
  * `upper-bound` (double): The upper limit of the confidence interval.
  * `t0` (double): The statistic calculated on the original data (from `boot-data`).

  See also [[bootstrap-stats]] for input preparation and other confidence interval methods:
  [[ci-normal]], [[ci-basic]], [[ci-percentile]], [[ci-bca]], [[ci-studentized]], [[ci-t]], [[quantiles]]."
  ([boot-data] (ci-bc boot-data 0.05))
  ([boot-data ^double alpha] (ci-bc boot-data alpha :legacy))
  ([{:keys [^double t0 ts]} ^double alpha estimation-strategy]
   (let [ats (m/seq->double-array ts)
         a (/ alpha 2.0)]
     (bca-common ats a (- 1.0 a) t0 estimation-strategy))))

(defn- acceleration [ts] (/ (stats/skewness ts :skew) -6.0))

(defn ci-bca
  "Calculates the Bias-Corrected and Accelerated (BCa) bootstrap confidence interval.

  The BCa interval is a sophisticated method that corrects for both bias and
  skewness in the distribution of the bootstrap statistic replicates. It is
  considered a more accurate interval, particularly when the bootstrap
  distribution is skewed.

  The calculation requires two components:
  1.  A **bias correction factor** ($z_0$) based on the proportion of bootstrap
      replicates less than the original statistic ($t_0$).
  2.  An **acceleration factor** ($a$) which quantifies the rate of change of the
      standard error of the statistic with respect to the true parameter value.

  The function uses one of two methods to calculate the acceleration factor:

  *   **Jackknife method**: If the input `boot-data` map contains the original
      `:data` and the `:statistic` function used to compute `:t0` and `:ts`,
      the acceleration factor is estimated using the jackknife method (by computing
      the statistic on leave-one-out jackknife samples).
  *   **Empirical method**: If `:data` or `:statistic` are missing from `boot-data`,
      the acceleration factor is estimated empirically from the distribution of
      the bootstrap replicates (`:ts`) using its skewness.

  Parameters:

  * `boot-data` (map): A map containing bootstrap results, typically from [[bootstrap-stats]].
    Requires keys:
      * `:t0` (double): The statistic calculated on the original data.
      * `:ts` (sequence of numbers): The statistic calculated on each bootstrap sample.
    May optionally include:
      * `:data` (sequence): The original dataset (required for jackknife acceleration).
      * `:statistic` (function): The function used to calculate the statistic (required for jackknife acceleration).
  * `alpha` (double, optional): The significance level for the interval.
    Defaults to `0.05` (for a 95% CI). The BCa method uses quantiles of the
    normal distribution and the bootstrap replicates, adjusted by the bias
    and acceleration factors.
  * `estimation-strategy` (keyword, optional): Specifies the quantile estimation strategy
    used to calculate the quantiles of the bootstrap replicates (`:ts`) for the
    final interval bounds after applying corrections. Defaults to `:legacy`.
    See [[quantiles]] for available options (e.g., `:r1` through `:r9`).

  Returns a vector `[lower-bound, upper-bound, t0]`.

  * `lower-bound` (double): The lower limit of the confidence interval.
  * `upper-bound` (double): The upper limit of the confidence interval.
  * `t0` (double): The statistic calculated on the original data (from `boot-data`).

  See also [[bootstrap-stats]] for input preparation and other confidence interval methods:
  [[ci-normal]], [[ci-basic]], [[ci-percentile]], [[ci-bc]], [[ci-studentized]], [[ci-t]], [[jackknife]], [[quantiles]]."
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
  "Calculates the Studentized (or Bootstrap-t) confidence interval.

  This method is based on the distribution of the studentized pivotal quantity
  ` (statistic(sample) - statistic(data)) / standard_error(statistic(sample)) `.
  It estimates the quantiles of this distribution using bootstrap replicates
  and then uses them to construct a confidence interval around the statistic calculated
  on the original data (`:t0`), scaled by the standard error of the statistic calculated
  on the original data (`stddev(:data)`).

  Parameters:

  * `boot-data` (map): A map containing bootstrap results and necessary inputs.
    This map typically comes from [[bootstrap-stats]] and augmented with `:data` and `:samples`
    from the original [[bootstrap]] call if not already present.
    Requires the following keys:
      * `:t0` (double): The statistic calculated on the original data.
      * `:ts` (sequence of numbers): The statistic calculated on each bootstrap sample.
      * `:data` (sequence): The original dataset used for bootstrapping. Needed to estimate the standard error of the statistic for scaling the interval.
      * `:samples` (collection of sequences): The collection of bootstrap samples. Needed to calculate the standard error of the statistic for each bootstrap sample.
  * `alpha` (double, optional): The significance level for the interval.
    Defaults to `0.05` (for a 95% CI). The interval is based on the `alpha/2`
    and `1 - alpha/2` quantiles of the studentized bootstrap replicates.
  * `estimation-strategy` (keyword, optional): Specifies the quantile estimation strategy
    used to calculate the quantiles of the studentized replicates. Defaults to `:legacy`.
    See [[quantiles]] for available options (e.g., `:r1` through `:r9`).

  Returns a vector `[lower-bound, upper-bound, t0]`.

  * `lower-bound` (double): The lower limit of the confidence interval.
  * `upper-bound` (double): The upper limit of the confidence interval.
  * `t0` (double): The statistic calculated on the original data (from `boot-data`).

  See also [[bootstrap-stats]] for input preparation and other confidence interval methods:
  [[ci-normal]], [[ci-basic]], [[ci-percentile]], [[ci-bc]], [[ci-bca]], [[ci-t]], [[stats/stddev]], [[stats/quantiles]]."
  ([boot-data] (ci-studentized boot-data 0.05))
  ([boot-data ^double alpha] (ci-studentized boot-data alpha :legacy))
  ([{:keys [^double t0 ts data samples]} ^double alpha estimation-strategy]
   (assert (and (seq data) (seq samples)) "Bootstrap samples can't be empty.")
   (let [a (/ alpha 2.0)
         z (map (fn [^double t s]
                  (/ (- t t0) (stats/stddev s))) ts samples)
         stddev (stats/stddev data)
         [^double q1 ^double q2] (stats/quantiles z [(- 1.0 a) a] estimation-strategy)]
     [(- t0 (* stddev q1))
      (- t0 (* stddev q2))
      t0])))

(defn ci-t
  "Calculates a confidence interval based on Student's t-distribution, centered at the original statistic value.

  This method constructs a confidence interval centered at the statistic calculated on the original data (`:t0`). The width of the interval is determined by the standard deviation of the bootstrap replicates (`:ts`), scaled by a critical value from a Student's t-distribution. The degrees of freedom for the t-distribution are based on the number of bootstrap replicates (`count(:ts) - 1`).

  This interval does not explicitly use the Studentized bootstrap pivotal quantity. Instead, it applies a standard t-interval structure using components derived from the bootstrap results and the original data.

  Parameters:

  * `boot-data` (map): A map containing bootstrap results, typically from [[bootstrap-stats]]. Requires keys:
      * `:t0` (double): The statistic calculated on the original data.
      * `:ts` (sequence of numbers): The statistic calculated on each bootstrap sample.
      May optionally include pre-calculated `:stddev` (standard deviation of `:ts`) for efficiency.
  * `alpha` (double, optional): The significance level for the interval. Defaults to `0.05` (for a 95% CI). The interval is based on the `alpha/2` and `1 - alpha/2` quantiles of the Student's t-distribution with `count(:ts) - 1` degrees of freedom.

  Returns a vector `[lower-bound, upper-bound, t0]`.

  * `lower-bound` (double): The lower limit of the confidence interval.
  * `upper-bound` (double): The upper limit of the confidence interval.
  * `t0` (double): The statistic calculated on the original data (from `boot-data`).

  See also [[bootstrap-stats]] for input preparation and other confidence interval methods:
  [[ci-normal]], [[ci-basic]], [[ci-percentile]], [[ci-bc]], [[ci-bca]], [[ci-studentized]]."
  ([boot-data] (ci-t boot-data 0.05))
  ([{:keys [^double t0 ts stddev]} ^double alpha]
   (let [a (/ alpha 2.0)
         ^double stddev (or stddev (stats/stddev ts))
         t (* stddev ^double (r/icdf (r/distribution :t {:degrees-of-freedom (dec (count ts))}) (- 1.0 a)))]
     [(- t0 t) (+ t0 t) t0])))

(m/unuse-primitive-operators)
