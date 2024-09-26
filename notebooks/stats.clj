^{:nextjournal.clerk/visibility :hide-ns
  :nextjournal.clerk/toc true} 
(ns stats
  {:clj-kondo/config '{:config-in-call {utils/table {:ignore [:unresolved-symbol]}
                                        utils/table2 {:ignore [:unresolved-symbol]}
                                        utils/table3 {:ignore [:unresolved-symbol]}}}}
  (:require [nextjournal.clerk :as clerk]
            [utils :as u]
            [fastmath.stats :as stats]
            [fastmath.random :as r]
            [fastmath.core :as m]
            [fastmath.kernel :as k]))

;; # fastmath.stats

;; Statistical functions

#_:clj-kondo/ignore
(require '[fastmath.stats :as stats])

;; ## Datasets

;; Two datasets are used in following chapters. They are taken from R:
;; * `mtcars`
;; * `iris`

^{::clerk/viewer clerk/table}
(u/mtcars)

^{::clerk/viewer clerk/table}
(u/iris)

;; Both `u/mtcars` and `u/iris` are helper functions which can return column and filter a dataset. Additionally `u/by` can group datasets by some variable(s).

;; * Select `:mpg` column:

(u/mtcars :mpg)

;; * Select `:mpg` column where `:am` is `0`:

(u/mtcars :mpg (comp zero? :am))

;; * Group by `:am` and select `:mpg`:

(u/by u/mtcars :am :mpg)

;; ## Descriptive statistics

;; Collection of basic statistical functions and their variants:
;; * `basic` - sum, max, min
;; * `means` - regular, geometric, harmonic, weighted and general power means
;; * `moments` - raw, central, absolute moments, skewness and kurtosis
;; * `deviations` - variance, standard deviation, MAD
;; * `quantiles` - quantile/percentile, median and their weighted variants
;; * `mode` - mode and modes for multimodal sample

;; ### Basic

^{::clerk/visibility :hide}
(u/table2
 [[minimum "minimum value"]
  [maximum "maximum value"]
  [sum "sum of values"]])

^{::clerk/visibility :hide}
(clerk/example
 (stats/minimum (u/mtcars :mpg))
 (stats/maximum (u/mtcars :mpg))
 (stats/sum (u/mtcars :mpg)))

;; ### Means

^{::clerk/visibility :hide}
(u/table2
 [[mean "mean of sample"]
  [geomean "geometric mean, for positive values only"]
  [harmean "harmonic mean"]
  [powmean "power mean"]
  [wmean "weighted mean"]])

^{::clerk/visibility :hide}
(clerk/example
 (stats/mean (u/mtcars :carb))
 (stats/geomean (u/mtcars :carb))
 (stats/harmean (u/mtcars :carb))
 (stats/wmean (u/mtcars :carb) (reverse (u/mtcars :carb))))

;; #### Power mean

;; Power mean generalizes mean calculation, exponent:
;; * `-1.0` - harmonic mean
;; * `0.0` - geometric mean
;; * `1.0` - regular mean

^{::clerk/visibility :hide}
(clerk/example
 (stats/powmean (u/mtcars :carb) -10)
 (stats/powmean (u/mtcars :carb) -1)
 (stats/powmean (u/mtcars :carb) 0)
 (stats/powmean (u/mtcars :carb) 0.1)
 (stats/powmean (u/mtcars :carb) 0.5)
 (stats/powmean (u/mtcars :carb) 1)
 (stats/powmean (u/mtcars :carb) 2)
 (stats/powmean (u/mtcars :carb) 10))

;; Plot of `powmean` for exponents in `[-5,10]` range.

^{::clerk/visibility :hide}
(u/fgraph (partial stats/powmean (u/mtcars :carb)) [-5 10])

;; ### Deviations

^{::clerk/visibility :hide}
(u/table2
 [[population-variance "population variance"]
  [variance "sample variance"]
  [population-stddev "population standard deviation"]
  [stddev "standard deviation"]
  [variation "coefficient of variation (CV)"]
  [median-absolute-deviation "MAD around median or any center"]
  [mean-absolute-deviation "average absolute deviation around the mean or any center"]
  [sem "standard error of mean"]])

;; Desity estimation of `:drat` column from `mtcars`

^{::clerk/visibility :hide}
(u/fgraph (k/kernel-density :gaussian (u/mtcars :drat)) [0 6] [0 1.1])

^{::clerk/visibility :hide}
(clerk/example
 (stats/population-variance (u/mtcars :drat))
 (stats/variance (u/mtcars :drat))
 (stats/population-stddev (u/mtcars :drat))
 (stats/stddev (u/mtcars :drat))
 (stats/variation (u/mtcars :drat))
 (stats/median-absolute-deviation (u/mtcars :drat))
 (stats/median-absolute-deviation (u/mtcars :drat) (stats/mean (u/mtcars :drat)))
 (stats/median-absolute-deviation (u/mtcars :drat) 0.0)
 (stats/mean-absolute-deviation (u/mtcars :drat))
 (stats/mean-absolute-deviation (u/mtcars :drat) (stats/median (u/mtcars :drat)))
 (stats/mean-absolute-deviation (u/mtcars :drat) 0.0)
 (stats/sem (u/mtcars :drat)))

;; #### Pooled variance

^{::clerk/visibility :hide}
(u/table2
 [[pooled-variance "pooled variance for samples"]
  [pooled-stddev "pooled standard deviation"]])

;; Pooled variance and standard deviations are calculated for many groups. They combine individual variances. There are three methods used:

;; * `:unbiased` (default) - $\sigma_p^2=\frac{\sum_{i=1}^m (n_i-1)\sigma_i^2}{\sum_{i=1}^m (n_i-1)}$
;; * `:biased` - $\sigma_p^2=\frac{\sum_{i=1}^m (n_i-1)\sigma_i^2}{\sum_{i=1}^m n_i}$
;; * `:avg` - $\sigma_p^2=\frac{\sum_{i=1}^m\sigma_i^2}{m}$, used for the same number of observation in each group.

;; Standard deviation is just square root of variance, $\sigma_p=\sqrt{\sigma_p^2}$.

;; More info on [wiki](https://en.wikipedia.org/wiki/Pooled_variance)

(def iris-petal-lengths (u/by u/iris :species :petal-length))

^{::clerk/visibility :hide}
(clerk/example
 (stats/pooled-variance iris-petal-lengths :unbiased)
 (stats/pooled-variance iris-petal-lengths :biased)
 (stats/pooled-variance iris-petal-lengths :avg)
 (stats/pooled-stddev iris-petal-lengths :unbiased)
 (stats/pooled-stddev iris-petal-lengths :biased)
 (stats/pooled-stddev iris-petal-lengths :avg))

;; ### Moments

^{::clerk/visibility :hide}
(u/table2
 [[moment "Calculate central or raw moment, can be absolute and normalized"]
  [skewness "Skewness, various algorithms"]
  [kurtosis "Kurtosis, various algorithms"]])

;; For the data we sample some skewed disctribution

(def moment-data (r/->seq (r/distribution :triangular {:c 0.9
                                                     :rng (r/rng :mersenne 1)}) 200))

^{::clerk/visibility :hide}
(clerk/example
 (stats/mean moment-data)
 (stats/median moment-data)
 (stats/mode moment-data :histogram))

;; Desity estimation

^{::clerk/visibility :hide}
(u/fgraph (k/kernel-density :gaussian moment-data 0.05) [-1.3 1.3])

;; #### moment

;; Calculate raw, central or standardized moment of any order with following variants:
;; * `:order` - order / degree of moment, default `2`
;; * `:center` - center value, set to `0` for raw moment, default mean (central)
;; * `:absolute?` - if true, uses absolute values of differences, default `false`
;; * `:mean?` - if true, returns mean, otherwise just sum of differences, default `true`
;; * `:normalize?` - divide result by standard deviation raised to order power, returns standardized moment, default `false`

^{::clerk/visibility :hide}
(clerk/table
 {:head [:absolute? :normalize? "order=1, center=0" "order=2" "order=3" "order=4" "order=5"]
  :rows [[false false
          (m/approx (stats/moment moment-data 1 {:center 0.0}) 4)
          (m/approx (stats/moment moment-data 2) 4)
          (m/approx (stats/moment moment-data 3) 4)
          (m/approx (stats/moment moment-data 4) 4)
          (m/approx (stats/moment moment-data 5) 4)]
         [true false
          (m/approx (stats/moment moment-data 1 {:center 0.0 :absolute? true}) 4)
          (m/approx (stats/moment moment-data 2 {:absolute? true}) 4)
          (m/approx (stats/moment moment-data 3 {:absolute? true}) 4)
          (m/approx (stats/moment moment-data 4 {:absolute? true}) 4)
          (m/approx (stats/moment moment-data 5 {:absolute? true}) 4)]
         [false true
          (m/approx (stats/moment moment-data 1 {:center 0.0 :normalize? true}) 4)
          (m/approx (stats/moment moment-data 2 {:normalize? true}) 4)
          (m/approx (stats/moment moment-data 3 {:normalize? true}) 4)
          (m/approx (stats/moment moment-data 4 {:normalize? true}) 4)
          (m/approx (stats/moment moment-data 5 {:normalize? true}) 4)]
         [true true
          (m/approx (stats/moment moment-data 1 {:normalize? true :center 0.0 :absolute? true}) 4)
          (m/approx (stats/moment moment-data 2 {:normalize? true :absolute? true}) 4)
          (m/approx (stats/moment moment-data 3 {:normalize? true :absolute? true}) 4)
          (m/approx (stats/moment moment-data 4 {:normalize? true :absolute? true}) 4)
          (m/approx (stats/moment moment-data 5 {:normalize? true :absolute? true}) 4)]]})

;; #### skewness

;; A measure of a skewness

;; Following symbols are used:
;; * $m_3$ - third central moment
;; * $m'_2$, $m'_3$ - second and third raw moment
;; * $s$ - standard deviation
;; * $n$ - sample size
;; * $\mu$ - mean
;; * $\nu$ - median
;; * $Q$ - quantile

;; Methods for mode from continuous data are: `:kde`, `:histogram` and `:pearson`. Refer [mode](#continuous) for more info.

^{::clerk/visibility :hide}
(clerk/table
 {:head ["Method" "Description" "Formula"]
  :rows [[:G1 "adjusted Fisher-Pearson, default" (clerk/md "$\\frac{n^2}{(n-1)(n-2)}\\frac{m_3}{s^3}$")]
         [:g1 "method of moments estimator" (clerk/md "$\\frac{n-2}{\\sqrt{n(n-1)}}G_1$")]
         [:b1 "natural estimator" (clerk/md "$\\frac{m_3}{s^3}$")]
         [:B1  "Yule skewness" (clerk/md "$\\frac{Q(0.75)+Q(0.25)-2Q(0.5)}{Q(0.75)-Q(0.25)}$")]
         [[:B1 'u] "generalized Yule formula" (clerk/md "$\\frac{Q(1-u)+Q(u)-2Q(0.5)}{Q(1-u)-Q(u)}$")]
         [:B3 "Groeneveld and Meeden's coefficient" (clerk/md "$\\frac{\\mu-\\nu}{E(|X-\\nu|)}$")]
         [:skew "BCa based skew" (clerk/md "$\\frac{1}{\\sqrt{n}}\\frac{m'_3}{m'_2}$")]
         [:mode "Pearson's mode skewness" (clerk/md "$\\frac{\\mu-\\operatorname{mode}}{s}$")]
         [[:mode 'method 'optional-args] "mode calculated from continuous data" (clerk/md "more info: [mode](#continuous)")]
         [:median "Pearson's median skewness" (clerk/md "$\\frac{3(\\mu-\\nu)}{s}$")]]})

;; Values of skewness for `moment-data`.

^{::clerk/visibility :hide}
(clerk/table
 {:head ["Method" "Value"]
  :rows (for [m [:G1 :g1 :b1 :B1 [:B1 0.1] [:B1 0.2] :B3 :skew
                 [:mode :histogram] [:mode :kde] [:mode :kde {:bandwidth 0.2}] [:mode :pearson] :median]]
          [m (str (stats/skewness moment-data m))])})

;; #### kurtosis

;; A measure of kurtosis

;; Following symbols are used:
;; * $m_2$, $m_4$ - second and fourth central moment
;; * $\mu$ - mean
;; * $s$ - standard deviation
;; * $n$ - sample size

^{::clerk/visibility :hide}
(clerk/table
 {:head ["Method" "Description" "Formula"]
  :rows [[:G2 "unbiased estimator, default" (clerk/md "$\\frac{n^2(n+1)}{(n-1)(n-2)(n-3)}\\frac{m_4}{s^4}-\\frac{3(n-1)^2}{(n-2)(n-3)}$")]
         [:g2 "biased estimator, excess kurtosis" (clerk/md "$\\frac{m_4}{m_2^2}-3$")]
         [:kurt "fourth standardized moment" (clerk/md "$\\frac{m_4}{m_2^2}$")]
         [:geary "Geary's kurtosis" (clerk/md "$\\frac{E(|X-\\mu|)}{\\sqrt{m_2}}$")]]})

;; Values of kurtosis for `moment-data`.

^{::clerk/visibility :hide}
(clerk/table
 {:head ["Method" "Value"]
  :rows (for [m [:G2 :g2 :kurt :geary]]
          [m (str (stats/kurtosis moment-data m))])})

;; ### Quantiles

;; Returns quantile(s) or percentile(s) for given sample. Different variants of estimation strategies are possible. Also weighted quantile algorithm is implemented, based on [spatstat.geom R package](https://search.r-project.org/CRAN/refmans/spatstat.geom/html/weighted.median.html).

;; List of the estimation strategies:
;; * `:legacy` (default) - described [here](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.html)
;; * `:r1`, `:r2`, `:r3`, `:r4`, `:r5`, `:r6`, `:r7`, `:r8`, `:r9` - described [here](https://en.wikipedia.org/wiki/Quantile#Estimating_quantiles_from_a_sample)

;; List of weighted quantile types:
;; * `:linear` (default) - linear interpolation
;; * `:step` - next largest value (step interpolation)
;; * `:average` - average of two nearest values

;; Default quantiles are `[0.25, 0.5, 0.75, 1.0]` or `[25 50 75 100] for percentiles.

^{::clerk/visibility :hide}
(u/table2
 [[quantile "quantile, input from 0 to 1"]
  [quantiles "list of quantiles"]
  [percentile "same as quantile, input from 0 to 100"]
  [percentiles "list of percentiles"]
  [median "median, quantile for q=0.5"]
  [median3 "median of three values"]
  [wquantile "weighted quantile"]
  [wquantiles "list of weighted quantiles"]
  [wmedian "weighted median"]])

(def data [1 3 6 6 6 50 6 7 7 12 12 13 17])

^{::clerk/visibility :hide}
(clerk/example
 (stats/quantile data 0.6)
 (map (fn [es] (m/approx (stats/quantile data 0.6 es)))
      [:legacy :r1 :r2 :r3 :r4 :r5 :r6 :r7 :r8 :r9])
 (stats/quantiles data)
 (stats/quantiles data [0.05 0.1 0.3 0.7 0.9 0.95])
 (stats/quantiles data [0.05 0.1 0.3 0.7 0.9 0.95] :r1)
 (stats/percentile data 69)
 (stats/percentiles data [25 50 75 100])
 (stats/median data)
 (stats/median data :r3))

(def weights [1 1 2 1 1 2 3 4 3 4 1 1 9])

^{::clerk/visibility :hide}
(clerk/example
 (stats/wquantile data weights 0.75)
 (stats/wquantile data weights 0.75 :step)
 (stats/wquantile data weights 0.75 :average)
 (stats/wquantiles data weights)
 (stats/wquantiles data weights [0.25 0.5 0.75 1] :step)
 (stats/wquantiles data weights [0.25 0.5 0.75 1] :average)
 (stats/wmedian data weights)
 (stats/wmedian data weights :step)
 (stats/wmedian data weights :average))

;; For even weights, `:step` method is equivalent to `:r1`

(def weights1 (repeat (count data) 1.0))
(def deciles [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0])

^{::clerk/visibility :hide}
(clerk/example
 (stats/wquantiles data weights1 deciles)
 (stats/wquantiles data weights1 deciles :step)
 (stats/wquantiles data weights1 deciles :average)
 (stats/quantiles data deciles)
 (stats/quantiles data deciles :r1)
 (stats/quantiles data deciles :r6))

;; ### Mode

;; Value or values that occurs most often in the sample.

^{::clerk/visibility :hide}
(u/table2
 [[mode "only first occurence when multimodal"]
  [modes "all modes as a sequence"]])

^{::clerk/visibility :hide}
(clerk/example
 (stats/mode [1 3 6 6 6 6 7 7 12 12 17])
 (stats/mode [1 1 2 4 4])
 (stats/modes [1 1 2 4 4]))

;; #### Continuous

;; There are three options to calculuate modes for continuous sample:
;; * `:histogram` - based on histogram, optional parameters are:
;;    * `:bins` - number of bins or method for calculation number of beans, default: `:rice`
;; * `:kde` - based on finite difference method on KDE function
;;    * `:kde` - kernel, default: `:gaussian`
;;    * `:bandwidth` - kernel width, default: auto
;; * `pearson` - $\operatorname{mode}=3\operatorname{median}-2\operatorname{mean}$
;;    * `:estimation-strategy` - estimation strategy for percentiles

;; All estimates can be highly innacurate and depends on samples and parameter tuning. `mode-data` should have mode=`0.75`.

;; Mode from histogram is calculated from the following formula

;; $$\operatorname{mode}=L+h(\frac{f_m-f_1}{2f_m-f_1-f_2})$$

;; Where:
;;  * $L$ - left bound of bin with highest frequency
;;  * $h$ - bin width
;;  * $f_m$ - highest frequency
;;  * $f_1$ - frequency in preceding bin
;;  * $f_2$ - frequency in succeeding bin

;; Generated, skewed data:

(def mode-data (r/->seq (r/distribution :triangular {:c 0.75 :rng (r/rng :mersenne 3)}) 200))

;; Density from data

^{::clerk/visibility :hide}
(clerk/table
 {:head ["gaussian" "epanechnikov, bandwidth=0.1" "epanechnikov, bandwidth=0.05"]
  :rows [[(u/fgraph (k/kernel-density :gaussian mode-data) [-1.5 1.5])
          (u/fgraph (k/kernel-density :epanechnikov mode-data 0.02) [-1.5 1.5])
          (u/fgraph (k/kernel-density :epanechnikov mode-data 0.1) [-1.5 1.5])]]})

^{::clerk/visibility :hide}
(clerk/example
 (stats/mode mode-data :pearson)
 (stats/mode mode-data :kde)
 (stats/mode mode-data :kde {:kde :epanechnikov :bandwidth 0.1})
 (stats/mode mode-data :kde {:kde :epanechnikov :bandwidth 0.05})
 (stats/mode mode-data :histogram)
 (stats/mode mode-data :histogram {:bins 30})
 (stats/mode mode-data :histogram {:bins :freedman-diaconis}))

;; ## Outliers

^{::clerk/visibility :hide}
(u/table2
 [[outliers "return lazy sequence of outliers, values outside inner fence"]])

;; Outliers are defined as samples which are outside inner fence. Let `Q1` is 25% percentile, `Q3` is 75% percentile. `IQR=Q3-Q1`. Then lower inner fence `LIF=Q1-1.5*IQR` and upper inner fence `UIF=Q3+1.5*IQR`.

(stats/outliers [0 0 0 1 2 3 4 5 1 2 3 4 1 2 3 1 2 1 10 20 30 40 50 -50])

;; ## Bootstrap

;; All bootstrap functions are moved to the `fastmath.stats.bootstrap` namespace.

;; ## Intervals

^{::clerk/visibility :hide}
(u/table2
 [[extent "min, max, mean"]
  [stddev-extent "mean-stddev, mean+stddev, mean"]
  [mad-extent "median-mad, median+mad, median"]
  [sem-extent "mean-sem, mean+sem, mean"]
  [quantile-extent "q1, q2, median"]
  [percentile-extent "p1, p2, median"]
  [pi-extent "sample size based quantile extent"]
  [pi "sample size based quantile extent as a map"]
  [hpdi-extent "higher posterior density interval"]
  [adjacent-values "LAV, UAV, median"]
  [inner-fence-extent "LIF, UIF, median"]
  [outer-fence-extent "LOF, UOF, median"]
  [ci "Student's t confidence interval for mean"]
  [bootstrap-ci "Basic percentile confidence interval from bootstrapped samples, deprecated"]
  [percentile-bca-extent "bias-corrected and accelerated"]
  [percentile-bc-extent "bias-corrected extent"]])

;; Some notes:

;; * `pi`, `pi-extent` and `hpdi-extent` are based on R `rethinking` package
;; * `pi` is actually `quantile-extent`
;; * `LAV` - lower adjacent value is the lowest data value which is greater than LIF
;; * `UAV` - upper adjacent value is the biggest data value which is lower than UIF
;; * `LIF`, `UIF` - lower/upper inner fence is $Q_1-1.5(Q_3-Q_1)$ and $Q_3+1.5(Q_3-Q_1)$
;; * `LOF`, `UOF` - lower/upper outer fence is $Q_1-3(Q_3-Q_1)$ and $Q_3+3(Q_3-Q_1)$

(def intervals-data [0 0 0 1 2 3 4 5 1 2 3 4 1 2 3 1 2 1 10 20 30 40 50 -50])

^{::clerk/visibility :hide}
(clerk/example
 (stats/extent intervals-data)
 (stats/stddev-extent intervals-data)
 (stats/mad-extent intervals-data)
 (stats/sem-extent intervals-data)
 (stats/quantile-extent intervals-data)
 (stats/quantile-extent intervals-data 0.2)
 (stats/quantile-extent intervals-data 0.3 0.9)
 (stats/percentile-extent intervals-data)
 (stats/percentile-extent intervals-data 20)
 (stats/percentile-extent intervals-data 30 90)
 (stats/percentile-bca-extent intervals-data)
 (stats/percentile-bca-extent intervals-data 5)
 (stats/percentile-bc-extent intervals-data)
 (stats/percentile-bc-extent intervals-data 5)
 (stats/pi-extent intervals-data)
 (stats/pi-extent intervals-data 0.9)
 (stats/pi intervals-data)
 (stats/pi intervals-data 0.9)
 (stats/hpdi-extent intervals-data)
 (stats/hpdi-extent intervals-data 0.5)
 (stats/adjacent-values intervals-data)
 (stats/inner-fence-extent intervals-data)
 (stats/outer-fence-extent intervals-data)
 (stats/ci intervals-data)
 (stats/bootstrap-ci intervals-data)
 (stats/bootstrap-ci intervals-data 0.9)
 (stats/bootstrap-ci intervals-data 0.9 50000)
 (stats/bootstrap-ci intervals-data 0.9 1000 stats/median))

;; ### Binomial

^{::clerk/visibility :hide}
(u/table2
 [[binamial-ci "Confidence intervals for binomial data"]])

;; Parameters:
;; * `number-of-successes`
;; * `number-of-trials`
;; * `method` - default: `:asymptotic`
;; * `alpha` - default: `0.05`

;; Full list of methods (use `:all` to get everything in a map)

stats/binomial-ci-methods

;; For `number-of successes=16` and `number-of-trials=100`, `p=0.16` here are the list of lower and upper interval values

^{::clerk/visibility :hide}
(clerk/table
 {:head ["method" "lower" "upper"]
  :rows (for [m stats/binomial-ci-methods
              :let [[a b] (stats/binomial-ci 16 100 m)]]
          [m a b])})

;; To get all values use `:all` as a method (the same `p=0.16` but higher number of trials)

(stats/binomial-ci 1600 10000 :all)

;; ### Range

^{::clerk/visibility :hide}
(u/table2
 [[iqr "interquartile range, difference between Q3 (75% percentile) and Q1 (25% percentile)"]
  [span "width of the samples, difference between maximum and minimum values"]])

^{::clerk/visibility :hide}
(clerk/example
 (stats/iqr intervals-data)
 (stats/span intervals-data))

;; ## Metrics

;; For given two samples calculate distance or value of given measure. Second argument can be either sequence or a number.

^{::clerk/visibility :hide}
(u/table2
 [[me "mean error"]
  [mae "mean absolute error"]
  [mape "mean absolute percentage error"]
  [rss "residual sum of squares"]
  [r2 "R2"]
  [mse "mean squared error"]
  [rmse "root mean squared error"]
  [psnr "peak signal to noise"]
  [count= "count equals, L0"]
  [L0 "count equals, count="]
  [L1 "manhattan distance"]
  [L2 "euclidean distance"]
  [L2sq "euclidean distance squared"]
  [LInf "chebyshev distance"]])

;; As a data we'll use `mpg` column and the same column with distortion

(def mpg (u/mtcars :mpg))
(def mpg+ (map (fn [v] (+ v (r/grand))) mpg))

^{::clerk/visibility :hide}
(clerk/example
  (stats/me mpg mpg+)
  (stats/me mpg 20.09)
  (stats/mae mpg mpg+)
  (stats/mae mpg 20.09)
  (stats/mape mpg mpg+)
  (stats/mape mpg 20.09)
  (stats/rss mpg mpg+)
  (stats/rss mpg 20.09)
  (stats/r2 mpg mpg+)
  (stats/r2 mpg 20.09)
  (stats/mse mpg mpg+)
  (stats/mse mpg 20.09)
  (stats/rmse mpg mpg+)
  (stats/rmse mpg 20.09)
  (stats/count= mpg mpg+)
  (stats/count= mpg 21)
  (stats/L0 mpg mpg+)
  (stats/L0 mpg 21)
  (stats/L1 mpg mpg+)
  (stats/L1 mpg 20.09)
  (stats/L2 mpg mpg+)
  (stats/L2 mpg 20.09)
  (stats/L2sq mpg mpg+)
  (stats/L2sq mpg 20.09)
  (stats/LInf mpg mpg+)
  (stats/LInf mpg 20.09))

;; ## Dissimilarity / similarity

;; There are two functions to compare PDFs of histograms against other histogram or data against PDF of any distribution.

^{::clerk/visibility :hide}
(u/table2
 [[dissimilarity "Distance between two histograms (PDFs)"]
  [similarity "Similarity of two histograms (PDFs)"]] )

;; As the input you can provide one of:
;;
;; * counts, which will be converted to probabilities (PDFs)
;; * data which will be compared to a distribution (histogram will be calculated for data)

;; Possible methods are desribed in the the [Comprehensive Survey on Distance/Similarity Measures between Probability Density Functions](https://www.naun.org/main/NAUN/ijmmas/mmmas-49.pdf) paper by Sung-Hyuk Cha or in the [philentropy R](https://search.r-project.org/CRAN/refmans/philentropy/html/distance.html) package.

(def normal-sample (repeatedly 1000 r/grand))
(def P-observed-sample [0 2 10 55 2 0 1])
(def Q-expected-sample [2 0 20 44 1 0 4])

(stats/dissimilarity :jensen-shannon P-observed-sample Q-expected-sample)
(stats/dissimilarity :jensen-shannon normal-sample (r/distribution :normal))

;; We can select also number of bins for second case (set to 5 bins only)

(stats/dissimilarity :jensen-shannon normal-sample (r/distribution :normal) {:bins 5})

^{::clerk/visibility :hide ::clerk/viewer u/unpaginated-table}
{:head ["method" "distance between P and Q" "distance between normal sample and distribution"]
 :rows (vec (for [ds [:euclidean, :city-block, :manhattan, :chebyshev, :minkowski, :sorensen, :gower, :soergel, :kulczynski, :canberra, :lorentzian, :non-intersection, :wave-hedges, :czekanowski, :motyka, :tanimoto, :jaccard, :dice, :bhattacharyya, :hellinger, :matusita, :squared-chord, :euclidean-sq, :squared-euclidean, :pearson-chisq, :chisq, :neyman-chisq, :squared-chisq, :symmetric-chisq, :divergence, :clark, :additive-symmetric-chisq, :kullback-leibler, :jeffreys, :k-divergence, :topsoe, :jensen-shannon, :jensen-difference, :taneja, :kumar-johnson, :avg]]
              [ds (stats/dissimilarity ds P-observed-sample Q-expected-sample)
               (stats/dissimilarity ds normal-sample r/default-normal)]))}

;; And similarity

(stats/similarity :ruzicka P-observed-sample Q-expected-sample)
(stats/similarity :ruzicka normal-sample (r/distribution :normal))

;; With changed number of bins to 5

(stats/similarity :ruzicka normal-sample (r/distribution :normal) {:bins 5})

^{::clerk/visibility :hide ::clerk/viewer u/unpaginated-table}
{:head ["method" "similarity of P and Q" "similarity of normal sample and distribution"]
 :rows (vec (for [ds [:intersection, :czekanowski, :motyka, :kulczynski, :ruzicka, :inner-product, :harmonic-mean, :cosine, :jaccard, :dice, :fidelity, :squared-chord]]
              [ds (stats/similarity ds P-observed-sample Q-expected-sample)
               (stats/similarity ds normal-sample r/default-normal)]))}


;; ## Correlation

^{::clerk/visibility :hide}
(u/table2
 [[covariance "covariance of two samples"]
  [correlation "correlation of two samples (Pearson's)"]
  [spearman-correlation "Spearman's correlation"]
  [pearson-correlation "Pearson's correlation"]
  [kendall-correlation "Kendall's correlation"]
  [kullback-leibler-divergence "Kullback-Leibler divergence, deprecated, use: dissimilarity function"]
  [jensen-shannon-divergence "Jensen-Shannon divergence, deprecated, use: dissimilarity function"]
  [correlation-matrix "correlation matrix"]
  [covariance-matrix "covariance martrix"]
  [coefficient-matrix "create any coefficient matrix for seq of seqs"]])

^{::clerk/visibility :hide}
(clerk/example
  (stats/covariance (u/mtcars :mpg) (u/mtcars :cyl))
  (stats/correlation (u/mtcars :mpg) (u/mtcars :cyl))
  (stats/pearson-correlation (u/mtcars :mpg) (u/mtcars :cyl))
  (stats/spearman-correlation (u/mtcars :mpg) (u/mtcars :cyl))
  (stats/kendall-correlation (u/mtcars :mpg) (u/mtcars :cyl))
  (stats/kullback-leibler-divergence (u/mtcars :mpg) (u/mtcars :cyl))
  (stats/jensen-shannon-divergence (u/mtcars :mpg) (u/mtcars :cyl))) 

;; Matrices are created by applying given correlation method (on any function) for every pair of given sequences. For:
;; * `covariance-matrix` - covariance
;; * `correlation-matrix` - correlation (or divergence), Pearson's by default, other options:
;;    * `:pearsons` - default
;;    * `:spearman`
;;    * `:kendall`
;; * `coefficient-matrix` - apply any function (`pearson-correlation` by default) to every pair of sequences. Last argument indicates if function is commutative or not (to speedup calculations)

(def three-irises (u/by u/iris :species :sepal-width))

^{::clerk/visibility :hide}
(clerk/example
 (stats/covariance-matrix three-irises)
 (stats/correlation-matrix three-irises)
 (stats/correlation-matrix three-irises :spearman)
 (stats/correlation-matrix three-irises :kendall)
 (stats/coefficient-matrix three-irises stats/L1))

;; ### ACF/PACF

^{::clerk/visibility :hide}
(u/table2
 [[acf "autocorrelation function"]
  [acf-ci "autocorrelation function with confidence intervals"]
  [pacf "partial autocorrelation function"]
  [pacf-ci "partial autocorrelation function with confidence intervals"]])

;; Please note: `acf` can return selected lags

^{::clerk/visibility :hide}
(clerk/example
 (stats/acf (u/mtcars :mpg))
 (stats/acf (u/mtcars :mpg) 5)
 (stats/acf (u/mtcars :mpg) [0 1 2 5 10])
 (stats/acf-ci (u/mtcars :mpg) 5)
 (stats/pacf (u/mtcars :mpg))
 (stats/pacf (u/mtcars :mpg) 5)
 (stats/pacf-ci (u/mtcars :mpg) 5))

^{::clerk/visibility :hide}
(clerk/table
 [["ACF (mpg)" "PACF (mpg)"]
  [(u/bgraph-int (zipmap (range) (stats/acf (u/mtcars :mpg))) [0 32] [-1.1 1.1])
   (u/bgraph-int (zipmap (range) (stats/pacf (u/mtcars :mpg))) [0 32] [-1.1 1.1])]
  ["ACF (uniform random)" "PACF (uniform random)"]
  [(u/bgraph-int (zipmap (range) (stats/acf (take-nth 10 (repeatedly 400 r/drand)))) [0 32] [-1.1 1.1])
   (u/bgraph-int (zipmap (range) (stats/pacf (take-nth 10 (repeatedly 400 r/drand)))) [0 32] [-1.1 1.1])]])

;; ## Contingency table

;; Contingency table is a map with frequencies for given tuples.

^{::clerk/visibility :hide}
(u/table2
 [[contingency-table "create table from variables"]
  [rows->contingency-table "create table from rows"]
  [contingency-table->marginals "calculate marginals"]])

;; For example:

^{::clerk/visibility :hide}
(clerk/table
 [["-" :right :left "TOTAL"]
  [:top 10 20 30]
  [:bottom 30 40 70]
  ["TOTAL" 40 60 100]])

;; is represented by

(def ct-example {[:top :right] 10 [:top :left] 20 [:bottom :right] 30 [:bottom :left] 40})

;; To calculate marginals call `contingency-table->marginals`. Diagonals are calculated when col/row names are the same.

(stats/contingency-table->marginals ct-example)

;; Contingency table can be build from:
;; * observations - the input contains values (as seqs) for each variable
;; * actual contingency table rows - seq of seqs

(stats/contingency-table [:a :b :a :a :a :b] [:c :a :a :a :c :b])
(stats/rows->contingency-table [[1 2 3] [2 3 4] [4 3 2]])

;; ### 2x2

;; #### Contingency

;; 2x2 contingency table is used to compare two groups and two outcomes. For example:

^{::clerk/visibility :hide}
(clerk/table
 [["" (clerk/md "`outcome 1`") (clerk/md "`outcome 2`")]
  [(clerk/md "`group 1`") "A" "B"]
  [(clerk/md "`group 2`") "C" "D"]])

;; In case of risk analysis in epidemiology or clinical trials, the table is structured as follows:

^{::clerk/visibility :hide}
(clerk/table
 [["" (clerk/md "`events`") (clerk/md "`non-events`")]
  [(clerk/md "`experimental group`") "A (EE)" "B (EN)"]
  [(clerk/md "`control group`") "C (CE)" "D (CN)"]])

;; 2x2 contingency table can be build from:
;; * 4 individual values for a,b,c and d
;; * sequence of 4 values (a,b,c,d)
;; * two pairs for group1 and group2
;; * a map with keys `:a`, `:b`, `:c` and `:d`

;; There are two functions:

^{::clerk/visibility :hide}
(u/table2
 [[contingency-2x2-measures-all "Returns a full collection of measures for 2x2 contingency table"]
  [contingency-2x2-measures "Returns a `:measures` key with odd-ratio added"]])

^{::clerk/opts {:auto-expand-results? true}}
(def ctable (stats/contingency-2x2-measures-all 41 64 216 180))

;; Functions return a map with following list of measures:

^{::clerk/visibility :hide}
(clerk/table
 {:head ["field name" "description" "value"]
  :rows [[:n "total number of observations" (:n ctable)]
         [:table "table itself" "see below"]
         [:marginals "table marginals" "see below"]
         [:proportions "proportions for cells, rows, cols and marginals" "see below"]
         [:expected "expected frequencies" "see below"]
         [:p-values "p-values for chi2 statistics" "see below"]
         [:OR "odds ratio" (:OR ctable)]
         [:lOR "log of odds ratio" (:lOR ctable)]
         [:SE "standard error for the log of odds ratio" (:SE ctable)]
         [:risk "epidemiology or medical treatment related measures, see below" "see below"]
         [:measures "collection of measures" "see below"]]})

^{::clerk/visibility :hide
  ::clerk/opts {:auto-expand-results? true}} 
(clerk/example
 (:table ctable)
 (:marginals ctable)
 (:proportions ctable)
 (:expected ctable)
 (:p-values ctable))

;; `:risk` key contains the following structure and values for example table

^{::clerk/visibility {:result :hide :code :hide}}
(def ct-risk (:risk ctable))

^{::clerk/visibility :hide}
(clerk/table
{:head ["field name" "description" "value"]
 :rows [[:RR "relative risk" (:RR ct-risk)]
        [:RD "risk difference" (:RD ct-risk)]
        [:ES "total subjects in experimental group" (:ES ct-risk)]
        [:CS "total subjects in control group" (:CS ct-risk)]
        [:EER "experimental event rate" (:EER ct-risk)]
        [:CER "control event rate" (:CER ct-risk)]
        [:ARR "absolute risk reduction" (:ARR ct-risk)]
        [:ARI "absolute risk increase" (:ARI ct-risk)]
        [:RRR "relative risk reduction" (:RRR ct-risk)]
        [:RRI "relative risk increase" (:RRI ct-risk)]
        [:NNT "number needed to treat" (:NNT ct-risk)]
        [:NNH "number needed to harm" (:NNH ct-risk)]
        [:AFe "attributable risk percent among the exposed" (:AFe ct-risk)]
        [:PFu "preventable fraction among the unexposed" (:PFu ct-risk)]]})

;; `:measures` kesy contains the following measures and values:

^{::clerk/visibility {:result :hide :code :hide}}
(def ct-m (:measures ctable))

^{::clerk/visibility :hide}
(clerk/table
 {:head ["measure" "description" "value"]
  :rows [[:chi2 "Pearson's chi squared" (:chi2 ct-m)]
         [:mcnemars-chi2 "McNemar’s chi-squared" (:mcnemars-chi2 ct-m)]
         [:yates "Yates corrected chi squared" (:yates ct-m)]
         [:cochran-mantel-haenszel "Cochran-Mantel-Haenszel" (:cochran-mantel-haenszel ct-m)]
         [:cohens-kappa "Cohen's kappa" (:cohens-kappa ct-m)]
         [:cohens-h "Cohen's H" (:cohens-h ct-m)]
         [:huberts-gamma "Hubert’s Γ" (:huberts-gamma ct-m)]
         [:yules-q "Yule's Q" (:yules-q ct-m)]
         [:yules-y "Yule's Y" (:yules-y ct-m)]
         [:holley-guilfords-g "Holley and Guilford’s G, SAC" (:holley-guilfords-g ct-m)]
         [:youdens-j "Youden's J" (:youdens-j ct-m)]
         [:cramers-v "Cramer's V" (:cramers-v ct-m)]
         [:bangdiwalas-b "Shankar and Bangdiwala’s B" (:bangdiwalas-b ct-m)]
         [:phi "Pearson's phi, Yule's phi, correlation" (:phi ct-m)]
         [:scotts-pi "Scott’s pi" (:scotts-pi ct-m)]
         [:F1 "F-score" (:F1 ct-m)]
         [:gwets-ac1 "Gwet’s AC1" (:gwets-ac1 ct-m)]
         [:PCC "Pearson’s Contingent Coefficient" (:PCC ct-m)]
         [:PCC-adjusted "adjusted Pearson’s Contingent Coefficient (PCC * sqrt(2))" (:PCC-adjusted ct-m)]
         [:TCC "Tetrachoric Correlation Coefficient" (:TCC ct-m)]]})

;; read more: [arxiv](https://arxiv.org/abs/2203.09628), [wiki OR](https://en.wikipedia.org/wiki/Odds_ratio), [stat pages](https://statpages.info/ctab2x2.html), [wiki tables](https://en.wikipedia.org/wiki/Contingency_table)

;; #### Confusion

;; Confusion matrix (error matrix, matching matrix) contains information about correctness of classification.
;; Matrix layout:

^{::clerk/visibility :hide}
(clerk/table
 [["" (clerk/md "`PP - predicted positive`") (clerk/md "`PN - predicted negative`")]
  [(clerk/md "`P - actual positive`") "TP - true positive" "FN - false negative"]
  [(clerk/md "`N - actual negative`") "FP - false positive" "TN - true negative"]])

^{::clerk/visibility :hide}
(u/table2
 [[->confusion-matrix "create confusion matrix from the input"]
  [binary-measures-all "full set confusion matrix measures"]
  [binary-measyres "selection of confusion matrix measures"]])

;; Confusion matrix can be created from:

;; * `tp`, `fn` `fp` and `tn` values
;; * a sequence of `tp`, `fn`, `fp` and `tn` values
;; * pair of sequences with `[tp fn]` and `[fp tn]`
;; * a map with `:tp`, `:fn`, `:fp` and `:tn` keys
;; * two sequences of actual and predicted values (as `true` and `false`)

;; In case if actual/predicted sequences do not contain true/false values, you can decode it with third argument:
;; * any sequence of true values
;; * map or function which will be applied to input

;; For example, if actual and predicted list contain `0` as `false` and `1` as `true`, you should provide one of the following as the third argument, to encode `1` as true:
;; * `#{1}`
;; * `(fn [v] (not (zero? v)))`
;; * `{1 true}`

^{::clerk/visibility :hide}
(clerk/example
 (stats/->confusion-matrix 10 20 30 40)
 (stats/->confusion-matrix [10 20 30 40])
 (stats/->confusion-matrix [[10 20] [30 40]])
 (stats/->confusion-matrix [true false true false true]
                           [false true true false true])
 (stats/->confusion-matrix [1 0 1 0 1]
                           [0 1 1 0 1] #{1})
 (stats/->confusion-matrix [1 0 1 0 1]
                           [0 1 1 0 1]
                           (fn [v] (not (zero? v))))
 (stats/->confusion-matrix [1 0 1 0 1]
                           [0 1 1 0 1] {1 true})
 (stats/->confusion-matrix [:a :b nil :c :d]
                           [nil nil :c :d nil]
                           identity))

;; List of all measures for given table.

^{::clerk/visibility :hide}
(clerk/table
 [["" "PP - predicted positive" "PN - predicted negative" "total"]
  ["P - actual positive" 50 20 70]
  ["N - actual negative" 10 30 40]
  ["total" 60 50 110]])

(def cm2x2 (stats/binary-measures-all {:tp 50 :fn 20 :fp 10 :tn 30}))

^{::clerk/visibility :hide ::clerk/viewer u/unpaginated-table}
{:head ["field name" "description" "value"]
 :rows [[:tp "true positive" (:tp cm2x2)]
        [:fn "false negative" (:fn cm2x2)]
        [:fp "false positive" (:fp cm2x2)]
        [:tn "true negative" (:tn cm2x2)]
        [:p "actual positive" (:p cm2x2)]
        [:n "actual negative" (:n cm2x2)]
        [:pp "predicted positive" (:pp cm2x2)]
        [:pn "predicted negative" (:pn cm2x2)]
        [:total "data count" (:total cm2x2)]

        [:tpr "true positive rate (recall, sensitivity, hit-rate)" (:tpr cm2x2)]
        [:recall "recall (tpr, sensitivity, hit-rate)" (:recall cm2x2)]
        [:sensitivity "sensitivity (tpr, recall, hit-rate)" (:sensitivity cm2x2)]
        [:hit-rate "hit rate (tpr, recall, sensitivity)" (:hit-rate cm2x2)]

        [:fnr "false negative rate (miss-rate)" (:fnr cm2x2)]
        [:miss-rate "miss rate (fnr)" (:miss-rate cm2x2)]

        [:fpr "false positive rate (fall-out)" (:fpr cm2x2)]
        [:fall-out "fall-out (fpr)" (:fall-out cm2x2)]

        [:tnr "true negative rate (specificity, selectivity)" (:tnr cm2x2)]
        [:specificity "specificity (tnr, selectivity)" (:specificity cm2x2)]
        [:selectivity "selectivity (tnr, specificity)" (:selectivity cm2x2)]
        
        [:ppv "positive predictive value" (:ppv cm2x2)]
        [:for "false omission rate" (:for cm2x2)]
        [:fdr "false discovery rate" (:fdr cm2x2)]
        [:npv "negative predictive value" (:npv cm2x2)]
        
        [:accuracy "accuracy" (:accuracy cm2x2)]
        [:ba "balanced accuracy" (:ba cm2x2)]
        [:prevalence "prevalence" (:prevalence cm2x2)]
        [:f1-score "F1 score (F measure)" (:f1-score cm2x2)]
        [:f-measure "F measure (F1 score)" (:f-measure cm2x2)]
        [:fm "Fowlkes-Mallows index" (:fm cm2x2)]
        [:mcc "Matthews correlation coefficient, phi" (:mcc cm2x2)]
        [:phi "Matthews correlation coefficient, phi" (:phi cm2x2)]
        [:kappa "Cohen's kappa" (:kappa cm2x2)]

        [:mk "markedness" (:mk cm2x2)]
        [:bm "bookmarker informedness" (:bm cm2x2)]
        [:pt "prevalence threshold" (:pt cm2x2)]
        [:ts "treat score, critical succes score, Jaccard index" (:ts cm2x2)]
        [:jaccard "treat score, critical succes score, Jaccard index" (:jaccard cm2x2)]
        
        [:lr+ "positive likelihood ratio" (:lr+ cm2x2)]
        [:lr- "negative likelihood ratio" (:lr- cm2x2)]
        [:dor "diagnostic odds ratio" (:dor cm2x2)]

        [:f-beta "F_beta score, a function" "-"]]}

;; `f-beta` is generalized F-score

;; $$F_{\beta}=(1+\beta^2)\frac{\operatorname{precision}\cdot\operatorname{recall}}{(\beta^2\cdot\operatorname{precision})+\operatorname{recall}}$$

^{::clerk/visibility :hide}
(clerk/example
 ((:f-beta cm2x2) 0.0)
 ((:f-beta cm2x2) 0.5)
 ((:f-beta cm2x2) 1.0)
 ((:f-beta cm2x2) 2.0))

;; more info [wiki](https://en.wikipedia.org/wiki/Confusion_matrix)

;; ## Effect Size

;; Relationship measures between two variables.

;; Data used in examples:

(def mpgs (u/by u/mtcars :am :mpg))

;; ### Differences

^{::clerk/visibility :hide}
(u/table2
[[cohens-d "Cohen's d"]
 [cohens-d-corrected "Cohen's d with bias correction"]
 [hedges-g "Hedges' g"]
 [hedges-g-corrected "Hedges' g with bias correction"]
 [hedges-g* "Hedges' g with bias correction using exact J term"]
 [glass-delta "Glass' delta"]
 [means-ratio "RoM, ratio of means"]
 [means-ratio-corrected "RoM, ratio of means with bias correction"]
 [cliffs-delta "Cliff's delta"]])

;; Please note:
;; * *Cohen's d* can be calculated using different pooled standard deviation methods.
;; * *Cehen's d* by default uses `:unbiased` pooled standard deviation, which is the same as in *Hedges' g*.
;; * Correction term `-corrected` versions is $1-\frac{3}{4k-1}$ where `k` are degrees of freedom ($n$ for biased and $n-2$ for unbiased pooled standard deviation)  
;; * Correction term for *Hedges' g** equals: $J(x)=\left(\frac{\Gamma(\frac{n-2}{2})}{\sqrt{\frac{n-2}{2}}\Gamma(\frac{n-3}{2})}\right)$.
;; * *Glass' $\Delta$* is not commutative, standarization is done by standard deviation of second group only.

;; More info:
;; * *Cohen's d* - [wiki](https://en.wikipedia.org/wiki/Effect_size#Cohen's_d)
;; * *Hedges' g* - [wiki](https://en.wikipedia.org/wiki/Effect_size#Hedges'_g)
;; * *Glass' $\Delta$* - [wiki](https://en.wikipedia.org/wiki/Effect_size#Glass'_%CE%94)
;; * *RoM* - [wiki](https://en.wikipedia.org/wiki/Ratio_estimator), [R effectsize](https://easystats.github.io/effectsize/reference/means_ratio.html)
;; * *Cliff's $\delta$* - [wiki](https://en.wikipedia.org/wiki/Effect_size#Effect_size_for_ordinal_data), [R effectsize](https://easystats.github.io/effectsize/reference/rank_biserial.html) 

^{::clerk/visibility :hide}
(clerk/example
(stats/cohens-d (mpgs 0) (mpgs 1))
(stats/cohens-d (mpgs 0) (mpgs 1) :biased)
(stats/cohens-d (mpgs 0) (mpgs 1) :avg)
(stats/cohens-d-corrected (mpgs 0) (mpgs 1))
(stats/cohens-d-corrected (mpgs 0) (mpgs 1) :biased)
(stats/hedges-g (mpgs 0) (mpgs 1))
(stats/hedges-g-corrected (mpgs 0) (mpgs 1))
(stats/hedges-g* (mpgs 0) (mpgs 1))
(stats/glass-delta (mpgs 0) (mpgs 1))
(stats/glass-delta (mpgs 1) (mpgs 0))
(stats/means-ratio (mpgs 0) (mpgs 1))
(stats/means-ratio-corrected (mpgs 0) (mpgs 1))
(stats/cliffs-delta (mpgs 0) (mpgs 1)))

;; ### Common Language

^{::clerk/visibility :hide}
(u/table2
[[ameasure "Vargha-Delaney A measure"]
 [wmw-odds "Wilcoxon-Mann-Whitney odds"]
 [p-overlap "overlapping index (PDF overlapping), kernel density estimation"]
 [cohens-u2 "proportion of one group exceeding the same proportion of the other group"]
 [cohens-u3 "proportion of group2 which is lower than median of group1"]])

;; * `p-overlap` uses kernel density estimation of groups and finds aproximate common area under PDFs.  Optional parameters are ([more info](https://www.frontiersin.org/articles/10.3389/fpsyg.2019.01089/full)):
;;     * `:kde` - density estinamation kernel, default: `:gaussian`
;;     * `:bandwidth` - kde bandwidth, default: `nil` (estimated)
;;     * `:min-iterations` - integration iterations, default: `3`
;;     * `:steps` - integration steps, default: `500`
;; * `cohens-u2` and `cohens-u3` are beased on quantiles and accept optional `estimation-strategy`
;; * `cohens-u2` uses `:brent` optimization to find a value

^{::clerk/visibility :hide}
(clerk/example
(stats/ameasure (mpgs 0) (mpgs 1))
(stats/wmw-odds (mpgs 0) (mpgs 1))
(stats/p-overlap (mpgs 0) (mpgs 1))
(stats/p-overlap (mpgs 0) (mpgs 1) {:kde :epanechnikov})
(stats/p-overlap (mpgs 0) (mpgs 0))
(stats/p-overlap (mpgs 0) (mpgs 0) {:kde :epanechnikov})
(stats/cohens-u2 (mpgs 0) (mpgs 1))
(stats/cohens-u2 (mpgs 0) (mpgs 1) :r7)
(stats/cohens-u3 (mpgs 0) (mpgs 1)))

;; ### Correlation

^{::clerk/visibility :hide}
(u/table2
[[pearson-r "R, Pearson's correlation"]
 [r2-determination "R^2, coefficient of determination"]
 [eta-sq "R^2"]
 [omega-sq "adjusted R^2"]
 [epsilon-sq "less biased R^2"]
 [cohens-f2 "f^2, based on eta, omega or epsilon"]
 [cohens-f "sqrt of f^2"]
 [cohens-q "comparison for two correlations"]])

;; Cohen's $f^2=\frac{R^2}{1-R^2}$

^{::clerk/visibility :hide}
(clerk/example
(stats/pearson-r (u/mtcars :cyl) (u/mtcars :mpg))
(stats/r2-determination (u/mtcars :cyl) (u/mtcars :mpg))
(stats/eta-sq (u/mtcars :cyl) (u/mtcars :mpg))
(stats/omega-sq (u/mtcars :cyl) (u/mtcars :mpg))
(stats/epsilon-sq (u/mtcars :cyl) (u/mtcars :mpg))
(stats/cohens-f2 (u/mtcars :cyl) (u/mtcars :mpg))
(stats/cohens-f2 (u/mtcars :cyl) (u/mtcars :mpg) :omega)
(stats/cohens-f2 (u/mtcars :cyl) (u/mtcars :mpg) :epsilon)
(stats/cohens-f (u/mtcars :cyl) (u/mtcars :mpg))
(stats/cohens-f (u/mtcars :cyl) (u/mtcars :mpg) :omega)
(stats/cohens-f (u/mtcars :cyl) (u/mtcars :mpg) :epsilon))

;; Cohen's q compares two pearson correlation values and equals $q=\operatorname{atanh}(r_1)-\operatorname{atanh}(r_2)$, possible arguments
;; * 2-arity - two correlation values
;; * 3-arity - correlation of group1 and group2 vs correlation of group1 and group3
;; * 4-arity - correlation of two first groups vs last two groups

^{::clerk/visibility :hide}
(clerk/example
(stats/cohens-q (stats/pearson-correlation
                 (u/iris :petal-width) (u/iris :petal-length))
                (stats/pearson-correlation
                 (u/iris :sepal-width) (u/iris :sepal-length)))
(stats/cohens-q (u/iris :petal-width) (u/iris :petal-length)
                (u/iris :sepal-width) (u/iris :sepal-length))
(stats/cohens-q (u/mtcars :mpg)
                (u/mtcars :cyl) (u/mtcars :am)))

;; ### Kruskal

;; Rank based effect size based on Kruskal-Wallis test

^{::clerk/visibility :hide}
(u/table2
[[rank-epsilon-sq "Epsilon squared R"]
 [rank-eta-sq "Eta squared H"]])

;; For:
;; * H - Kruskal-Wallis test H statistic
;; * n - data size
;; * k - number of groups

;; $$\Epsilon^2_R=\frac{H(n+1)}{n^2-1}$$
;; $$\eta^2_H=\frac{H-k+1}{n-k}$$

(def mpgs-by-cyl (u/by u/mtcars :cyl :mpg))

^{::clerk/visibility :hide}
(clerk/example
 (stats/rank-epsilon-sq mpgs-by-cyl)
 (stats/rank-eta-sq mpgs-by-cyl))

;; ### Contingency

;; Effect size or measures of contingency tables.

^{::clerk/visibility :hide}
(u/table2
 [[mcc "Mathew correlation coefficient, phi"]
  [cramers-c "Cramer's C"]
  [cramers-v "Cramer's V"]
  [cramers-v-corrected "Cramer's V corrected"]
  [cohens-w "Cohen's W"]
  [cohens-kappa "Cohen's kappa"]
  [weighted-kappa "Cohen's weighted kappa"]
  [tschuprows-t "Tschuprow's T"]])

^{::clerk/visibility :hide}
(def ct-table
  (stats/rows->contingency-table [[150 100 165 130]
                                  [50 65 35 10]
                                  [2 55 40 25]]))

^{::clerk/visibility :hide}
(clerk/example
 (stats/cramers-c ct-table)
 (stats/cramers-v ct-table)
 (stats/cramers-v-corrected ct-table)
 (stats/cohens-w ct-table)
 (stats/cohens-kappa ct-table)
 (stats/tschuprows-t ct-table))

;; Weighted kappa assumes that contingency table is indexed, rows and columns should be integers. There are two predefined methods for weights:
;; * `:equal-spacing` - linear weights, default (also known as Cicchetti-Allison)
;; * `:fleiss-cohen` - quadratic weights

;; You can also provide own weights as map or custom 3-arity function: `R` (maximum index), `row-id` and `col-id`.

(def wk-table
  (stats/rows->contingency-table [[11 3 1 0]
                                  [1 9 0 1]
                                  [0 1 10 0]
                                  [1 2 0 10]]))

^{::clerk/visibility :hide}
(clerk/example
 (stats/cohens-kappa wk-table)
 (stats/weighted-kappa wk-table)
 (stats/weighted-kappa wk-table :fleiss-cohen)
 (stats/weighted-kappa wk-table {[0 0] 1 [1 1] 1 [2 2] 1 [3 3] 1})
 (stats/weighted-kappa wk-table (fn [R row col]
                                  (- 1.0 (m/sqrt (/ (m/abs (- row col)) R))))))

;; ## Scaling and trimming

^{::clerk/visibility :hide}
(u/table2
 [[demean "subtract mean from the data"]
  [standardize "normalize samples to make mean=0 and stddev=1"]
  [robust-standardize "normalize samples to make median=0 and MAD=1"]
  [winsor "constrain data to given quantiles: [q, 1-q]"]
  [trim "remove data which are outside given quantiles: [q, 1-q]"]
  [rescale "lineary scale data to given range"]])

(def some-data [1 2 3 6 21])

^{::clerk/visibility :hide}
(clerk/example
 (stats/mean some-data)
 (stats/stddev some-data)
 (stats/median some-data)
 (stats/median-absolute-deviation some-data))

^{::clerk/visibility :hide}
(clerk/example
 (stats/demean some-data)
 (stats/mean (stats/demean some-data))
 (stats/stddev (stats/demean some-data))
 (stats/standardize some-data)
 (stats/mean (stats/standardize some-data))
 (stats/stddev (stats/standardize some-data))
 (stats/robust-standardize some-data)
 (stats/median (stats/robust-standardize some-data))
 (stats/median-absolute-deviation (stats/robust-standardize some-data))
 (stats/winsor some-data)
 (stats/winsor some-data 0.4)
 (stats/trim some-data)
 (stats/trim some-data 0.4)
 (stats/rescale some-data)
 (stats/rescale some-data 10 100))

;; `trim` and `winsor` also replace NaN value with median. 4-arity version accept data, low value, high value and NaN replacement.

^{::clerk/visibility :hide}
(clerk/example
 (stats/winsor [1 2 ##NaN 4 5 6 7 10])
 (stats/winsor [1 2 ##NaN 4 5 6 7 10] 2 6 nil)
 (stats/trim [1 2 ##NaN 4 5 6 7 10])
 (stats/trim [1 2 ##NaN 4 5 6 7 10] 2 6 nil))

;; ## Histogram

^{::clerk/visibility :hide}
(u/table2
 [[histogram "build histogram data"]
  [estimate-bins "estimate number of bins using different algorithms"]])

;; `histogram` for given data and (optional) number of bins or bins' estimation algorithm returns a map with following content:
;; * `:size` - actual number of bins
;; * `:step` - bin width
;; * `:min` - minimal value
;; * `:max` - maximum value
;; * `:samples` - total data size
;; * `:bins` - seq of pairs containing lower value of the bin and count

;; If the last argument is pair `[min value, max value]` input data is filtered to be inside the given range

;; Here is the list of bin estimation method (`:freedman-diaconis` is default)

^{::clerk/visibility :hide}
(clerk/table
 {:head ["method" "bin count for sepal-length (iris)"]
  :rows (for [m [:freedman-diaconis :sqrt :sturges :rice :doane :scott]]
          [m (stats/estimate-bins (u/iris :sepal-length) m)])})

^{::clerk/visibility :hide}
(clerk/example
 (stats/histogram (u/iris :sepal-length))
 (stats/histogram (u/iris :sepal-length) :sqrt)
 (stats/histogram (u/iris :sepal-length) 20)
 (stats/histogram (u/iris :sepal-length) :doane [4.8 6.2]))

^{::clerk/visibility :hide}
(clerk/table
 [["default" ":sqrt"]
  [(u/bgraph-int (u/h->bars (u/iris :sepal-length)) [4 8] [0 40])
   (u/bgraph-int (u/h->bars (u/iris :sepal-length) :sqrt) [4 8] [0 40])]
  ["20 bins" ":doane, 4.8<=v<=6.2"]
  [(u/bgraph-int (u/h->bars (u/iris :sepal-length) 20) [4 8] [0 40])
   (u/bgraph-int (u/h->bars (u/iris :sepal-length) :doane [4.8 6.2]) [4 8] [0 40])]])

;; ## Tests

;; Collection of various statistical tests

^{::clerk/visibility :hide}
(u/table2
 [[p-value "calculate p-value"]])

;; `p-value` calculates probability of any statistic for given distribution, N(0,1) by default. `p-value` can be calculated as:
;; * `:two-sided` or `:both` - using both sides of distribution (default)
;; * `:one-sided-greater` or `:right` - using right side of the distribution
;; * `:one-sided`, `:one-sided-lower` or `:left` - using left side of the distribution

^{::clerk/visibility :hide}
(clerk/example
 (stats/p-value 2.0)
 (stats/p-value (r/distribution :t) 2.0 :right)
 (stats/p-value (r/distribution :t) 2.0 :left)
 (stats/p-value (r/distribution :t) 2.0 :both))

;; All tests accept a map:
;; * `:sides` to indicate which distribution side to use to calculate `p-value`
;; * `:alpha` significance, if confidence interval is calculated

;; All tests return a map with calculated test information and input configuration, common keys are:
;; * `:stat` - test statistic
;; * `:p-value` - probability value of test statistic
;; * `:sides` - sides
;; * `:n` - data size (or sizes if more than one sample is used)
;; * `:confidence-intarval` - confidence interval (if applicable)

;; Data used:

mpg
(def sepal-widths (vec (u/by u/iris :species :sepal-width)))


;; ### T and Z

;; One and two samples (paired or equal/unequal variance) `t-test` and `z-test`. 

;; #### One sample

^{::clerk/visibility :hide}
(u/table2
 [[t-test-one-sample "one sample Student's t-test"]
  [z-test-one-sample "one sample z-test"]])

;; Additional functions arguments:
;; * `:mu` - mean (default: 0.0)

;; Functions return (additionally):
;; * `:df` - degrees of freedom (for `t-test`)
;; * `:estimate` - mean estimation
;; * `:t` or `:z` - test statistic
;; * `:stderr` - standard error

^{::clerk/visibility :hide}
(clerk/example
 (stats/t-test-one-sample mpg)
 (stats/z-test-one-sample mpg)
 (stats/t-test-one-sample mpg)
 (stats/z-test-one-sample mpg))

;; #### Two samples

^{::clerk/visibility :hide}
(u/table2
 [[t-test-two-samples "two samples Student's t-test"]
  [z-test-two-samples "two samples z-test"]])

;; Addtitional functions arguments:
;; * `:mu` - mean (default: 0.0)
;; * `:paired?` - paired test? (default: `false`)
;; * `:equal-variances?` - equal or unequal variances (default: `false`)

;; Functions return (additionally):
;; * `:nx`, `:ny` - samples sizes
;; * `:df` - degrees of freedom (for `t-test`)
;; * `:estimate` - difference in means
;; * `:estimated-mu` - estimated means of `xs` and `ys`
;; * `:t` :or `:z` - test statistic
;; * `:stderr` - standard error

^{::clerk/visibility :hide}
(clerk/example
 (stats/t-test-two-samples (sepal-widths 0) (sepal-widths 1))
 (stats/t-test-two-samples (sepal-widths 0) (sepal-widths 1) {:paired? true})
 (stats/t-test-two-samples (sepal-widths 0) (sepal-widths 1) {:equal-variances? true})
 (stats/z-test-two-samples (sepal-widths 0) (sepal-widths 1))
 (stats/z-test-two-samples (sepal-widths 0) (sepal-widths 1) {:paired? true})
 (stats/z-test-two-samples (sepal-widths 0) (sepal-widths 1) {:equal-variances? true}))

;; ### Binomial

^{::clerk/visibility :hide}
(u/table2
 [[binomial-test "binomial test"]])

;; Test if data are from binomial distribution with `p` (default: 0.5) probability.
;; As a data function accepts one of below:
;; * `xs` - sequence of true/false values
;; * `number-of-successes` and `number-of-trials`

;; Additional parameters include:
;; * `:p` - probability (default: 0.5)
;; * `:ci-method` - one of the binomial confidence interval methods (default: `:asymptotic`).

;; List of confidence interval methods:

stats/binomial-ci-methods

;; Data:

(def binomial-data [true false true true true false true true true])

^{::clerk/visibility :hide}
(clerk/example
 (stats/binomial-test binomial-data)
 (stats/binomial-test binomial-data {:p 0.2})
 (stats/binomial-test 100 201)
 (stats/binomial-test 100 201 {:ci-method :wilson}))

;; ### Variances

;; #### F

;; F-test analyses variances ratio

^{::clerk/visibility :hide}
(u/table2
 [[f-test "two samples f-test"]])

;; Function returns (additionally):
;; * `:F` - test statistic
;; * `:df` - degreeses of freedom

^{::clerk/visibility :hide}
(clerk/example
 (stats/f-test (sepal-widths 0) (sepal-widths 1))
 (stats/f-test (sepal-widths 1) (sepal-widths 0)))

;; #### ANOVA

;; Variance tests for multiple samples

^{::clerk/visibility :hide}
(u/table2
 [[one-way-anova-test "One way ANOVA"]
  [levene-test "Levene test"]
  [brown-forsythe-test "Brown Forsythe test"]
  [fligner-killeen-test "Flinger Killen test"]])

;; Functions return (additionally)
;; * `:SSt` - sum of squares (treatments)
;; * `:SSe` - sum of squares (error)
;; * `:DFt` - degrees of freedom (treatments)
;; * `:DFe` - degrees of freedom (error)
;; * `:MSt` - mean square (treatments)
;; * `:MSe` - mean square (error)
;; * `:F`, `:W` or `:chi2` - statistic (F - anova, W - levene or brown-forsythe, chi2 - flinger-killeen)

;; Please note:
;; * `levene-test` can accept `statistic` parameter (statistical function for spread, default `mean`) and `scorediff` function applied to score (spread) (default: `abs`)
;; * `brown-forsythe-test` uses `median` as the statistic

^{::clerk/visibility :hide}
(clerk/example
 (stats/one-way-anova-test sepal-widths)
 (stats/levene-test sepal-widths)
 (stats/brown-forsythe-test sepal-widths)
 (stats/fligner-killeen-test sepal-widths))

;; ### Power divergence

;; Family of power divergence tests

^{::clerk/visibility :hide}
(u/table3
 "lambda"
 [[power-divergence-test "general, power divergence test" "$\\lambda$"]
  [chisq-test "Chi2 test" 1]
  [multinomial-likelihood-ratio-test "multinomial test, log likelihood test, G-test" 0]
  [minimum-discrimination-information-test "modified log likelihood test" -1]
  [neyman-modified-chisq-test "Neyman chi2 test, lambda=-2" -2]
  [freeman-tukey-test "Freeman-Tukey test" -1/2]
  [cressie-read-test "Cressie-Read test" 2/3]])

;; Divergence tests work accept contingency table or sequence of numbers (in this case goodness-of-fit is performed for given `p` probabilities). Additional options:

;; * `:lambda` - test parameter (for `power-divergence-test` only)
;; * `:p` - probabilites/weights for goodness-of-fit test (default: $p_i=\frac{1}{n}$)
;; * `:ci-sides` - sides for confidence intervals (the same values as for `:sides`) (default: `:both`)
;; * `:ddof` - degrees of freedom correction (default: 0.0)
;; * `:bootstrap-samples` - bootstrap samples for confidence intervals

;; ### Contingency

(def ctab (stats/rows->contingency-table [[1 2 10]
                                        [3 20 50]
                                        [1 4 30]]))

^{::clerk/visibility :hide}
(clerk/example
  (stats/power-divergence-test ctab {:lambda 2})
  (stats/chisq-test ctab)
  (stats/multinomial-likelihood-ratio-test ctab)
  (stats/minimum-discrimination-information-test ctab)
  (stats/neyman-modified-chisq-test ctab)
  (stats/freeman-tukey-test ctab)
  (stats/cressie-read-test ctab))

;; Chart of `chi2` statistic for different `lambda`
^{::clerk/visibility :hide}
(u/fgraph (fn [l] (:chi2 (stats/power-divergence-test ctab {:lambda l :bootstrap-samples 1}))) [-4 4] [0 nil])

;; ### Goodness of fit

;; When goodness of fit is performed, `p` should contain probabilities or weights for given counts. There is a possibility to use power divergence tests for continuous data, see the next section.

(def counts [1 1 2 10])

^{::clerk/visibility :hide}
(clerk/example
  (stats/chisq-test counts)
  (stats/chisq-test counts {:p [1 10 2 50]})
  (stats/chisq-test counts {:p [1 1 2 11]}))

;; Chart of `chi2` statistic for different `lambda`
^{::clerk/visibility :hide}
(u/fgraph (fn [l] (:chi2 (stats/power-divergence-test counts {:lambda l :bootstrap-samples 1}))) [-4 4] [0 nil])

;; ### Distribution

^{::clerk/visibility :hide}
(u/table2
 [[ad-test-one-sample "Anderson-Darling test for one sample"]
  [ks-test-one-sample "Kolmogorov-Smirnov test for one sample"]
  [ks-test-two-samples "Kolmogorov-Smirnov test for two samples"]])

;; For one sample test, verify if a sample is from given distribution. Second parameter can be actual distribution (default: N(0,1)) or data. In case of data `continuous-distribution` will be created. Parameters:
;; * `:kernel` - kde kernel (see `fastmath.random` for kernels) or `:enumerated`
;; * `:bandwidth` - bandwidth for kde
;; * `:distinct?` - make `xs` distinct for Kolmogorov-Smirnov

(def normal-data1 (repeatedly 50 r/grand))
(def normal-data2 (repeatedly 100 r/grand))

^{::clerk/visibility :hide}
(clerk/example
  (stats/ad-test-one-sample normal-data1)
  (stats/ad-test-one-sample normal-data1 normal-data1)
  (stats/ad-test-one-sample normal-data1 normal-data2)
  (stats/ad-test-one-sample normal-data1 normal-data2 {:kernel :epanechnikov})
  (stats/ks-test-one-sample normal-data1)
  (stats/ks-test-one-sample normal-data1 normal-data1)
  (stats/ks-test-one-sample normal-data1 normal-data2)
  (stats/ks-test-one-sample normal-data1 normal-data2 {:kernel :epanechnikov})
  (stats/ks-test-two-samples normal-data1 normal-data1)
  (stats/ks-test-two-samples normal-data1 normal-data2)) 

;; #### Goodness of fit

;; Power divergence tests can be used on continuous data. In such case, histogram will be created for input and `p` option has to be a distribution object.

(stats/chisq-test normal-data2 {:p (r/distribution :normal)})

;; ### Rank

^{::clerk/visibility :hide}
(u/table2
 [[kruskal-test "Kruskal-Wallis rank sum test"]])

^{::clerk/visibility :hide}
(clerk/example
 (stats/kruskal-test sepal-widths))

;; ## List of symbols

^{::clerk/visibility :hide}
(u/make-public-fns-table 'fastmath.stats)
