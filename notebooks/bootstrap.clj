^{:nextjournal.clerk/visibility :hide-ns
  :nextjournal.clerk/toc true} 
(ns bootstrap
  {:clj-kondo/config '{:config-in-call {utils/table {:ignore [:unresolved-symbol]}
                                        utils/table2 {:ignore [:unresolved-symbol]}
                                        utils/table3 {:ignore [:unresolved-symbol]}}}}
  (:require [nextjournal.clerk :as clerk]
            [utils :as u]
            [fastmath.stats :as stats]
            [fastmath.stats.bootstrap :as boot]
            [fastmath.random :as r]
            [fastmath.core :as m]
            [fastmath.kernel :as k]))

;; # fastmath.stats.bootstrap

;; Bootstrap methods and confidence interval for bootstrapped data. Algorithms are mostly based on `boot` and `bootstrap` R packages.

#_:clj-kondo/ignore
(require '[fastmath.stats.bootstrap :as boot])

;; Datasets used:

;; Petal width for Setosa from `iris` dataset

(def iris-data (first (u/by u/iris :species :petal-width)))

;; Times between failures of the air-conditioning equipment in Boeing 720 aircraft, `aircondit` dataset from `boot` R package.

(def aircondit [3 5 7 18 43 85 91 98 100 130 230 487])

;; 1920 and 1930 population from 10 US cities, `city` dataset from `boot` R package

(def city [[138 143] [93 104] [61 69] [179 260] [48 75]
         [37 63] [29 50] [23 48] [30 111] [2 50]])

;; Multidimensional data, `:mpg` from `mtcars` grouped by `:am`

(def mtcars-pair (u/by u/mtcars :am :mpg))

;; ## Generating samples

^{::clerk/visibility :hide}
(u/table2
 [[bootstrap "generate n random samples and apply statistic"]
  [bootstrap-stats "analyse bootstrap samples"]
  [jackknife "generate leave-one-out samples"]
  [jackknife+ "generate positive jackknife (add-one) samples"]])

;; ### Bootstrap 

;; `bootstrap` function is the main entry point for getting samples and statistics. It can be used for both nonparametric and parametric bootstrap.

;; #### Input

;; Two variants for samples are possible:
;; * sequence of values or sequence of tuples (for multidimensional sampling)
;; * map with two keys:
;;     * `:data` - sequence of values or sequence of tuples
;;     * `:model` (optional) - distribution object (see `fastmath.random/distribution`) or 0-arity function which returns a value when called. For multidimensional sampling, sequence of models should be provided.

;; `statistic` if provided will be called on bootstrapped samples, if `statistic` is nil, only samples and model (if available) are returned.

;; Parameters:
;; * `:samples` - number of bootstrapped samples (default: 500)
;; * `:size` - forced size of individual sample (default: same as source, unless jackknife used)
;; * `:method`
;;    *  `nil` (default) - uniform random
;;    * `:jackknife` for leave-one-out jackknife
;;    * `:jackknife+` for positive jackknife
;;    * any method accepted by `fastmath.random/->seq`
;; * `:rng` - random number generator (see: `fastmath.random/rng`)
;; * `:smoothing` - smoothing bootstrap:
;;    * `:kde` - kernel density estimation, additional options are: `:kernel` (default) and `:bandwidth` (auto)
;;    * `:gaussian` - add random value from N(0,standard error)
;; * `:distribution` - distribution used to auto generate model (distribution) from data:
;;    * `:real-discrete-distribution` - default
;;    * `:integer-discrete-distribution` - for integer values
;;    * `:categorical-distribution` - for any other type
;; * `:antithetic?` - antithetic sampling (default: false)
;; * `:dimensions` - if set to `:multi` - multidimensional data and models are created
;; * `:include?` - if set to `true` (default: `false`) original dataset is included in samples

;; #### Output

;; Function returns a map
;; * `:data` - orignal dataset
;; * `:samples` - sequence of samples
;; * `:model` - model used

;; #### Sampling

;; By default uniform data sampling is done. For jackknife, set `:jackknife` or `:jackknife+` (see below).

;; If `:antithetic?` is set to true, samples are drawn in pairs with complementary probabilities from model iCDF.

;; If `:dimensions` is set to `:multi`, each dimension is sampled independently.

(boot/bootstrap mtcars-pair nil {:dimensions :multi :samples 2})

;; `:size` can be used to force different than default sample size.

(boot/bootstrap aircondit nil {:size 50 :samples 2})

;; #### Nonparametric

;; For nonparametric case the input is the sequence of values or a map with `:data` key containing sequence. Unless jackknife sampling method is selected, model will be created automatically as `:real-discrete-distribution` assuming that data values are real numbers. 

(boot/bootstrap aircondit)

;; * For integers use `:integer-discrete-distribution`
;; * For any other object use `:categorical-distribution`, for example `city` dataset

(boot/bootstrap aircondit nil {:distribution :integer-discrete-distribution :samples 2})
(boot/bootstrap city nil {:distribution :categorical-distribution :samples 2})

;; If `:kde` is set for `:smoothing`: `:continuous-distribution` will be created as a model.

(boot/bootstrap aircondit nil {:smoothing :kde :samples 2})

;; Eventually `:jackknife` and `:jackknife+` can be used as a sampling method. In this case `:model` is not created:

(boot/bootstrap aircondit nil {:method :jackknife})
(boot/bootstrap aircondit nil {:method :jackknife+})

;; #### Parametric

;; For parametric bootstrap model should be provided as a `:model` key in input map. A model can be any distribution object and 0-arity function which returns random values.

(def aircondit-model (r/distribution :exponential {:mean (stats/mean aircondit)}))

(boot/bootstrap {:data aircondit :model aircondit-model} nil {:samples 2})

;; #### Smoothing

;; By setting `:smoothing` parameter, smoothed bootstrap is selected, two options are possible:
;; * `:kde` - kernel density estimation is applied (continuous distribution is created), this method can be used when model is not provided
;; * `:gaussian` - after sampling, gaussian noise $N(0,se)$ is added to samples, where `se` is standard error of mean

(boot/bootstrap iris-data nil {:smoothing :gaussian :samples 2})
(boot/bootstrap {:data aircondit :model aircondit-model} nil {:smoothing :gaussian :samples 2})

;; #### Statistics

;; When statistic function is provided, additionally following information is returned:
;; * `:t0` - statistic of original data
;; * `:ts` - statistic of bootstrapped data
;; * `:mean` - mean of `ts`
;; * `:bias` - difference between mean of `ts` and `t0`
;; * `:variance` - variance of `ts`
;; * `:stddev` - standard deviation of `ts` (`boot` R package names this as standard error)
;; * `:sem` - standard error of mean of `ts`
;; * `:median` - median of `ts`

(dissoc (boot/bootstrap aircondit stats/mean {:samples 10}) :samples :data)

;; `bootstrap-stats` function also can be used, a map with `:data` and `:samples` keys and statistic function is necessary. 

(let [data [1 2 3 4 5 10 100]]
  (boot/bootstrap-stats {:data data
                         :samples (boot/jackknife+ data)} stats/stddev))

;; ### Jackknife

;; `jackknife` and `jackknife+` just generate seq of samples. They can be called indirectly by `bootstrap`. Number of generated samples is the same as data size. Sample size is one less than original for `jackknife` and one more for `jackknife+`.

^{::clerk/visibility :hide}
(clerk/example
 (boot/jackknife [1 2 3 4])
 (boot/jackknife+ [1 2 3 4]))

;; ## Confidence intervals

^{::clerk/visibility :hide}
(u/table2
 [[ci-normal "normal bias-corrected"]
  [ci-basic "basic percentile"]
  [ci-percentile "percentile"]
  [ci-bc "bias-corrected"]
  [ci-bca "bias-corrected and accelerated"]
  [ci-studentized "studentized data"]
  [ci-t "Student's T"]])

;; All functions accept result returned by `bootstrap` function. Most of the functions accept a map containing `:t0` and `:ts` keys only.

;; Some notes:
;; * `ci-normal` uses `bias` to correct intervals
;; * `ci-bca` acceleration term is calculated from jackknife, if input doesn't contain original data or statistic function, acceleration is calculated from `ts`
;; * `ci-studentized` calculates variance for each bootstrapped sample

^{::clerk/visibility :hide}
(u/bgraph-int (u/h->bars iris-data) [0 1] [0 30])

(def bootstrap-result (boot/bootstrap iris-data stats/skewness {:rng (r/rng :jdk 7)
                                                              :samples 200}))

;; Skewness of the data (`t0`)

(:t0 bootstrap-result)

;; Mean of ts

(:mean bootstrap-result)

;; Histogram and density of `ts`

^{::clerk/visibility :hide}
(clerk/table
 [[(u/bgraph-int (u/h->bars (:ts bootstrap-result) 20) [-0.2 3] [0 33])
   (u/fgraph-int (k/kernel-density :gaussian (:ts bootstrap-result)) [-0.2 3] [0 nil])]])

^{::clerk/visibility :hide}
(clerk/table
 {:head ["method" "lower" "upper"]
  :rows (for [[m f] [["normal" boot/ci-normal]
                     ["basic" boot/ci-basic]
                     ["percentile" boot/ci-percentile]
                     ["bc" boot/ci-bc]
                     ["bca (acc from jackknife)" boot/ci-bca]
                     ["bca (acc from ts)" (fn [in] (boot/ci-bca (dissoc in :data)))]
                     ["studentized" boot/ci-studentized]
                     ["t" boot/ci-t]]
              :let [[a b] (f bootstrap-result)]]
          [m a b])})

;; Confidence intervals illustrated:
;; * `t0` - green
;; * `ts` - blue, as density
;; * confidence interval - orange

^{::clerk/visibility {:code :hide :result :hide}}
(defn ci-chart
  [bootstrap-result ci v [dx dy :as dom] r]
  (let [d (k/kernel-density :gaussian (:ts bootstrap-result))
        t0 (:t0 bootstrap-result)
        mx (d t0)
        diff (/ (- dy dx) 100)
        diff2 (/ diff 10)
        mkf (fn [b] (let [[a b] (b bootstrap-result)]
                     (fn [x] (cond
                              (< a x b) v
                              (or (m/delta-eq a x diff)
                                  (m/delta-eq b x diff)) 0.0
                              :else ##NaN))))]
    (u/fgraphs-int [d (mkf ci) (fn [v] (cond
                                        (m/delta-eq v t0 diff2) mx
                                        (m/delta-eq v t0 diff) 0.0
                                        :else ##NaN))]
                   dom r)))

^{::clerk/visibility :hide}
(let [g (fn [b] (ci-chart bootstrap-result b 0.1 [-0.2 3] [0 nil]))]
  (clerk/table
   [["normal" "basic" "percentile"]
    [(g boot/ci-basic) (g boot/ci-basic) (g boot/ci-percentile)]
    ["bc" "bca (jackkinfe)" "bca (ts)"]
    [(g boot/ci-bc) (g boot/ci-bca) (g (fn [in] (boot/ci-bca (dissoc in :data))))]
    ["studentized" "t"]
    [(g boot/ci-studentized) (g boot/ci-t)]]))

;; # Examples

;; ## aircondit

aircondit

;; Let's estimate mean with nonparametric and parametric bootstrap.

(def aircondit-means-np (boot/bootstrap aircondit stats/mean {:rng (r/rng :jdk 10)}))
(def aircondit-means-s (boot/bootstrap aircondit stats/mean {:rng (r/rng :jdk 10)
                                                           :smoothing :gaussian}))
(def aircondit-means-kde (boot/bootstrap aircondit stats/mean {:rng (r/rng :jdk 10)
                                                             :smoothing :kde}))
(def aircondit-means-p (boot/bootstrap {:data aircondit
                                      :model (r/distribution :exponential {:mean (stats/mean aircondit)})}
                                     stats/mean
                                     {:rng (r/rng :jdk 10)}))

^{::clerk/visibility :hide}
(clerk/table
 {:head ["key" "nonparametric" "parametric" "smoothed" "kde"]
  :rows (for [k [:t0 :mean :bias :median :stddev]]
          [k (aircondit-means-np k) (aircondit-means-p k) (aircondit-means-s k) (aircondit-means-kde k)])})

;; Density of `ts` with BCa confidence intervals

^{::clerk/visibility :hide}
(let [g (fn [d] (ci-chart d boot/ci-bca 0.001 [0 300] [0 0.015]))]
  (clerk/table
   [["nonparametric" "parametric"]
    [(g aircondit-means-np) (g aircondit-means-p)]
    ["gaussian smoothing" "kde smoothing"]
    [(g aircondit-means-s) (g aircondit-means-kde)]]))

;; ## city

city

;; City data set contains pairs of data, we will use reatio of means statistic

(defn rom
  [data]
  (/ (stats/mean (map second data))
     (stats/mean (map first data))))

;; Our `t0` is:

(rom city)

;; We need to use categorical distribution to deal with pairs of numbers.

(def city-bootstrap (boot/bootstrap city rom {:distribution :categorical-distribution
                                            :rng (r/rng :jdk 3)}))

^{::clerk/visibility :hide}
(clerk/table
 {:head ["key" "nonparametric"]
  :rows (for [k [:t0 :mean :bias :median :stddev]]
          [k (city-bootstrap k)])})

(ci-chart city-bootstrap boot/ci-bca 0.15 [0.5 2.8] [0 nil])

;; ## mtcars

mtcars-pair

;; `mtcars-pair` is 2-dimensinal data (two series) which can be sampled separately and used to evaluate `cohens-d-corrected` effect size

(def mtcars-bootstrap (boot/bootstrap mtcars-pair stats/cohens-d-corrected {:dimensions :multi
                                                                          :rng (r/rng :jdk 2)}))
(def mtcars-bootstrap-s (boot/bootstrap mtcars-pair stats/cohens-d-corrected {:dimensions :multi
                                                                            :smoothing :gaussian
                                                                            :rng (r/rng :jdk 2)}))
(def mtcars-bootstrap-a (boot/bootstrap mtcars-pair stats/cohens-d-corrected {:dimensions :multi
                                                                            :antithetic? true
                                                                            :rng (r/rng :jdk 2)}))


^{::clerk/visibility :hide}
(clerk/table
 {:head ["key" "nonparametric" "smoothed" "antithetic"]
  :rows (for [k [:t0 :mean :bias :median :stddev]]
          [k (mtcars-bootstrap k) (mtcars-bootstrap-s k) (mtcars-bootstrap-a k)])})

;; Density of `ts` with BCa confidence intervals

^{::clerk/visibility :hide}
(let [g (fn [d] (ci-chart (dissoc d :data) boot/ci-bca 0.1 [-3 0] [0 1]))]
  (clerk/table
   [["nonparametric" "smoothed" "antithetic"]
    [(g mtcars-bootstrap) (g mtcars-bootstrap-s) (g mtcars-bootstrap-a)]]))

;; # List of symbols

^{::clerk/visibility :hide}
(u/make-public-fns-table 'fastmath.stats.bootstrap)
