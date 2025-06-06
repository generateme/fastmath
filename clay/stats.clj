^:kindly/hide-code
(ns stats
  (:require [fastmath.stats :as stats]
            [fastmath.stats.bootstrap :as boot]
            [fastmath.dev.ggplot :as gg]
            [fastmath.dev.codox :as codox]
            [fastmath.dev.clay :as utls]
            [fastmath.dev.dataset :as ds]
            [fastmath.core :as m]
            [scicloj.kindly.v4.kind :as kind]
            [fastmath.random :as r]
            [fastmath.kernel :as kernel]
            [fastmath.vector :as v]))

;; # Statistics {.unnumbered}

;; This notebook provides a comprehensive overview and examples of the statistical
;; functions available in the `fastmath.stats` namespace. It covers a wide array
;; of descriptive statistics, measures of spread, quantiles, moments, correlation,
;; distance metrics, contingency table analysis, binary classification metrics,
;; effect sizes, various statistical tests (normality, binomial, t/z, variance,
;; goodness-of-fit, ANOVA, autocorrelation), time series analysis tools (ACF, PACF),
;; data transformations, and histogram generation. Each section introduces the
;; relevant concepts and demonstrates the use of `fastmath.stats` functions
;; with illustrative datasets (`mtcars`, `iris`, `winequality`).

;; ::: {.callout-warning}
;; This notebook is written with the support of Gemini LLM models:
;; 
;; * `gemini-2.5-pro-exp-03-25`
;; * `gemini-2.5-flash-preview-04-17`
;;
;; I did my best to verify the output of LLMs however I don't guarantee absence of the model hallucinations or incorrectnesses.
;; :::

;; ## Datasets

;; To illustrate statistical functions three dataset will be used:

;; * `mtcars`
;; * `iris`
;; * `winequality`

;; A `fastmath.dev.dataset` as `ds` will be used to access data.

;; Select the `:mpg` column from `mtcars`

(ds/mtcars :mpg)

;; Select the `:mpg` column with a row predicate

(ds/mtcars :mpg (fn [row] (and (< (row :mpg) 20.0) (zero? (row :am)))))

;; Group by the `:am` column and select `:mpg`

(ds/by ds/mtcars :am :mpg)

;; ### `mtcars`

;; 11 attributes of 32 different cars comparison

(kind/table (ds/mtcars) {:use-datatables true})

;; ### `iris`

;; Sepal and Petal iris species comparison 

(kind/table (ds/iris) {:use-datatables true})

;; ### `winequality`

;; White [Portugese Vinho Verde](https://archive.ics.uci.edu/dataset/186/wine+quality) wine consisting set of quality parameters (physicochemical and sensory).

;; ### Data

;; Let's bind some common columns to a global vars

;; #### `mtcars`

;; Miles per US gallon

(def mpg (ds/mtcars :mpg))

;; Horsepower

(def hp (ds/mtcars :hp))

;; Weight of the car

(def wt (ds/mtcars :wt))

(kind/table
 [[(gg/->image (gg/density mpg {:title "MPG - Miles Per Gallon"}))
   (gg/->image (gg/density hp {:title "HP - Horsepower"}))
   (gg/->image (gg/density wt {:title "WT - Weight"}))]])

;; #### `iris`

;; Sepal lengths

(def sepal-lengths (ds/by ds/iris :species :sepal-length))
(def setosa-sepal-length (sepal-lengths :setosa))
(def virginica-sepal-length (sepal-lengths :virginica))

(kind/table
 [[(gg/->image (gg/density setosa-sepal-length {:title "Setosa sepal length"}))
   (gg/->image (gg/density virginica-sepal-length {:title "Virginica sepal length"}))]])

;; Sepal widths

(def sepal-widths (ds/by ds/iris :species :sepal-width))
(def setosa-sepal-width (sepal-widths :setosa))
(def virginica-sepal-width (sepal-widths :virginica))

(kind/table
 [[(gg/->image (gg/density setosa-sepal-width {:title "Setosa sepal width"}))
   (gg/->image (gg/density virginica-sepal-width {:title "Virginica sepal width"}))]])


;; #### `winequality`

(def residual-sugar (ds/winequality "residual sugar"))
(def alcohol (ds/winequality "alcohol"))

(kind/table
 [[(gg/->image (gg/density residual-sugar {:title "Residual sugar"}))
   (gg/->image (gg/density alcohol {:title "Alcohol"}))]])


;; ## Basic Descriptive Statistics

;; General summary statistics for describing the center and range of data.

;; ::: {.callout-tip title="Defined functions"}
;; * `minimum`, `maximum`
;; * `sum`
;; * `mean`
;; * `geomean`, `harmean`, `powmean`
;; * `mode`, `modes`
;; * `wmode`, `wmodes`
;; * `stats-map`
;; :::

;; ### Basic

;; This section covers fundamental descriptive statistics including finding the smallest (`minimum`) and largest (`maximum`) values, and calculating the total (`sum`).

(utls/examples-note
  (stats/minimum mpg)
  (stats/maximum mpg)
  (stats/sum mpg))

;; Compensated summation can be used to reduce numerical error. There are three algorithms implemented:

;; * `:kahan`: The classic algorithm using a single correction variable to reduce numerical error.
;; * `:neumayer`: An improvement on Kahan, also using one correction variable but often providing better accuracy.
;; * `:klein`: A higher-order method using two correction variables, typically offering the highest accuracy at a slight computational cost.

;; As you can see below, all compensated summation give accurate result for `mpg` data.

(utls/examples-note
  (stats/sum mpg)
  (stats/sum mpg :kahan)
  (stats/sum mpg :neumayer)
  (stats/sum mpg :klein))

;; But here is the example in which normal summation and `:kahan` fails.

(utls/examples-note
  (stats/sum [1.0 1.0e100 1.0 -1.0e100])
  (stats/sum [1.0 1.0e100 1.0 -1.0e100] :kahan)
  (stats/sum [1.0 1.0e100 1.0 -1.0e100] :neumayer)
  (stats/sum [1.0 1.0e100 1.0 -1.0e100] :klein))

;; ### Mean

;; The mean is a measure of central tendency. `fastmath.stats` provides several types of means:
;; 
;; * **Arithmetic Mean** (`mean`): The sum of values divided by their count. It's the most common type of average.

;; $$\mu = \frac{1}{n} \sum_{i=1}^{n} x_i$$
;; 
;; * **Geometric Mean** (`geomean`): The n-th root of the product of n numbers. Suitable for averaging ratios, growth rates, or values that are multiplicative in nature. Requires all values to be positive. It's less affected by extreme large values than the arithmetic mean, but more affected by extreme small values.

;; $$G = \left(\prod_{i=1}^{n} x_i\right)^{1/n} = \exp\left(\frac{1}{n}\sum_{i=1}^{n} \ln(x_i)\right)$$

;; * **Harmonic Mean** (`harmean`): The reciprocal of the arithmetic mean of the reciprocals of the observations. Appropriate for averaging rates (e.g., speeds). It is sensitive to small values and requires all values to be positive and non-zero.

;; $$H = \frac{n}{\sum_{i=1}^{n} \frac{1}{x_i}}$$

;; * **Power Mean** (`powmean`): Also known as the generalized mean or Hölder mean. It generalizes the arithmetic, geometric, and harmonic means. Defined by a power parameter $p$.

;; $$M_p = \left(\frac{1}{n} \sum_{i=1}^{n} x_i^p\right)^{1/p} \text{ for } (p \neq 0)$$

;; Special cases:

;; *   $p \to 0$: Geometric Mean
;; *   $p = 1$: Arithmetic Mean
;; *   $p = -1$: Harmonic Mean
;; *   $p = 2$: Root Mean Square (RMS)
;; *   $p \to \infty$: Maximum value
;; *   $p \to -\infty$: Minimum value

;; The behavior depends on $p$: higher $p$ gives more weight to larger values, lower $p$ gives more weight to smaller values.

(utls/examples-note
  (stats/mean residual-sugar)
  (stats/geomean residual-sugar)
  (stats/harmean residual-sugar)
  (stats/powmean residual-sugar ##-Inf)
  (stats/powmean residual-sugar -4.5)
  (stats/powmean residual-sugar -1.0)
  (stats/powmean residual-sugar 0.0)
  (stats/powmean residual-sugar 1.0)
  (stats/powmean residual-sugar 4.5)
  (stats/powmean residual-sugar ##Inf))

;; All values of power mean for `residual-sugar` data and range of the power from `-5` to `5`.

(gg/->image (gg/function (partial stats/powmean residual-sugar)
                         {:x [-5 5]
                          :title "Power mean for residual-sugar"
                          :xlab "power"
                          :ylab "powmean"}))

;; #### Weighted

;; Every mean function accepts optional weights vector. Formulas for weighted means are as follows.

;; * **Arithmetic Mean** (`mean`):

;; $$\mu_w = \frac{\sum_{i=1}^{n} w_i x_i}{\sum_{i=1}^{n} w_i}$$

;; * **Geometric Mean** (`geomean`):

;; $$G_w = \left(\prod_{i=1}^{n} x_i^{w_i}\right)^{1/\sum w_i} = \exp\left(\frac{\sum_{i=1}^{n} w_i \ln(x_i)}{\sum_{i=1}^{n} w_i}\right)$$

;; * **Harmonic Mean** (`harmean`):

;; $$H_w = \frac{\sum_{i=1}^{n} w_i}{\sum_{i=1}^{n} \frac{w_i}{x_i}}$$

;; * **Power Mean** (`powmean`):

;; $$M_{w,p} = \left(\frac{\sum_{i=1}^{n} w_i x_i^p}{\sum_{i=1}^{n} w_i}\right)^{1/p} \text{ for } (p \neq 0)$$

;; Let's calculate mean of `hp` (horsepower) weighted by `wt` (car weight).

(utls/examples-note
  (stats/mean hp)
  (stats/mean hp wt)
  (stats/geomean hp)
  (stats/geomean hp wt)
  (stats/harmean hp)
  (stats/harmean hp wt)
  (stats/powmean hp -4.5)
  (stats/powmean hp wt -4.5)
  (stats/powmean hp 4.5)
  (stats/powmean hp wt 4.5))

;; #### Expectile

;; The `expectile` is a measure of location, related to both the [[mean]] and [[quantile]]. For a given level `τ` (tau, a value between 0 and 1), the `τ`-th expectile is the value `t` that minimizes an asymmetrically weighted sum of *squared* differences from `t`. This is distinct from quantiles, which minimize an asymmetrically weighted sum of *absolute* differences.

;; A key property of expectiles is that the 0.5-expectile is identical to the arithmetic [[mean]]. As `τ` varies from 0 to 1, expectiles span a range of values, typically from the minimum (`τ=0`) to the maximum (`τ=1`) of the dataset. Like the mean, expectiles are influenced by the magnitude of all data points, making them generally more sensitive to outliers than corresponding quantiles (e.g., the median).

(utls/examples-note
  (stats/expectile residual-sugar 0.1)
  (stats/expectile residual-sugar 0.5)
  (stats/mean residual-sugar)
  (stats/expectile residual-sugar 0.9)
  (stats/expectile hp wt 0.25))

;; Plotting expectiles for `residual-sugar` data across the range of `τ` from 0.0 to 1.0.

(gg/->image (gg/function (partial stats/expectile residual-sugar)
                         {:title "Expectiles for residual-sugar"
                          :xlab "τ (tau)"
                          :ylab "Expectile value"}))

;; ### Mode

;; The mode is the value that appears most frequently in a dataset.

;; For numeric data, `mode` returns the first mode encountered, while `modes` returns a sequence of all modes (in increasing order for the default method).

;; Let's see that mode returns elements with the highest frequency of `hp` (showing only first 10 values)

^:kind/hidden
(->> (frequencies hp)
     (sort-by second >)
     (take 10))

(kind/table
 {:column-names ["value" "frequency"]
  :row-vectors (->> (frequencies hp)
                    (sort-by second >)
                    (take 10))})

(utls/examples-note
 (stats/mode hp)
 (stats/modes hp))

;; When dealing with data potentially from a continuous distribution, these functions can estimate the mode using different methods:
;; 
;; *   `:histogram`: Mode(s) based on the peak(s) of a histogram.
;; *   `:pearson`: Mode estimated using Pearson's second skewness coefficient (`mode ≈ 3 * median - 2 * mean`).
;; *   `:kde`: Mode(s) based on Kernel Density Estimation, finding original data points with the highest estimated density.
;; *   The default method finds the exact most frequent value(s), suitable for discrete data.

(utls/examples-note
  (stats/mode residual-sugar)
  (stats/mode residual-sugar :histogram)
  (stats/mode residual-sugar :pearson)
  (stats/mode residual-sugar :kde))

;; For weighted data, or data of any type (not just numeric), use `wmode` and `wmodes`.
;; `wmode` returns the first weighted mode (the one with the highest total weight encountered first), and `wmodes` returns all values that share the highest total weight. If weights are omitted, they default to 1.0 for each value, effectively calculating unweighted modes for any data type.

(utls/examples-note
  (stats/wmode [:a :b :c :d] [1 2.5 1 2.5])
  (stats/wmodes [:a :b :c :d] [1 2.5 1 2.5]))

;; ### Stats

;; The `stats-map` function provides a comprehensive summary of descriptive statistics for a given dataset. It returns a map where keys are statistic names (as keywords) and values are their calculated measures. This function is a convenient way to get a quick overview of the data's characteristics.

;; The resulting map contains the following key-value pairs:

;; *   `:Size` - The number of data points in the sequence.
;; *   `:Min` - The smallest value in the sequence (see [[minimum]]).
;; *   `:Max` - The largest value in the sequence (see [[maximum]]).
;; *   `:Range` - The difference between the maximum and minimum values.
;; *   `:Mean` - The arithmetic average of the values (see [[mean]]).
;; *   `:Median` - The middle value of the sorted sequence (see [[median]]).
;; *   `:Mode` - The most frequently occurring value(s) (see [[mode]]).
;; *   `:Q1` - The first quartile (25th percentile) of the data (see [[percentile]]).
;; *   `:Q3` - The third quartile (75th percentile) of the data (see [[percentile]]).
;; *   `:Total` - The sum of all values in the sequence (see [[sum]]).
;; *   `:SD` - The sample standard deviation, a measure of data dispersion around the mean.
;; *   `:Variance` - The sample variance, the square of the standard deviation.
;; *   `:MAD` - The Median Absolute Deviation, a robust measure of variability (see [[median-absolute-deviation]]).
;; *   `:SEM` - The Standard Error of the Mean, an estimate of the standard deviation of the sample mean.
;; *   `:LAV` - The Lower Adjacent Value, the smallest observation that is not an outlier (see [[adjacent-values]]).
;; *   `:UAV` - The Upper Adjacent Value, the largest observation that is not an outlier (see [[adjacent-values]]).
;; *   `:IQR` - The Interquartile Range, the difference between Q3 and Q1.
;; *   `:LOF` - The Lower Outer Fence, a threshold for identifying extreme low outliers (`Q1 - 3.0 * IQR`).
;; *   `:UOF` - The Upper Outer Fence, a threshold for identifying extreme high outliers (`Q3 + 3.0 * IQR`).
;; *   `:LIF` - The Lower Inner Fence, a threshold for identifying mild low outliers (`Q1 - 1.5 * IQR`).
;; *   `:UIF` - The Upper Inner Fence, a threshold for identifying mild high outliers (`Q3 + 1.5 * IQR`).
;; *   `:Outliers` - A list of data points considered outliers (values outside the inner fences, see [[outliers]]).
;; *   `:Kurtosis` - A measure of the "tailedness" or "peakedness" of the distribution (see [[kurtosis]]).
;; *   `:Skewness` - A measure of the asymmetry of the distribution (see [[skewness]]).

(utls/zp
 (stats/stats-map residual-sugar))

;; ## Quantiles and Percentiles

;; Statistics related to points dividing the distribution of data.

;; ::: {.callout-tip title="Defined functions"}
;; * `percentile`, `percentiles`
;; * `quantile`, `quantiles`
;; * `wquantile`, `wquantiles`
;; * `median`, `median-3`
;; * `wmedian`
;; :::

;; Quantiles and percentiles are statistics that divide the range of a probability distribution into continuous intervals with equal probabilities, or divide the observations in a sample in the same way.

;; `fastmath.stats` provides several functions to calculate these measures:

;; *   `percentile`: Calculates the p-th percentile (a value from 0 to 100) of a sequence.
;; *   `percentiles`: Calculates multiple percentiles for a sequence.
;; *   `quantile`: Calculates the q-th quantile (a value from 0.0 to 1.0) of a sequence. This is equivalent to `(percentile vs (* q 100.0))`.
;; *   `quantiles`: Calculates multiple quantiles for a sequence.
;; *   `median`: Calculates the median (0.5 quantile or 50th percentile) of a sequence.
;; *   `median-3`: A specialized function that calculates the median of exactly three numbers.
;;
;; All `percentile`, `quantiles`, `quantile`, `quantiles`, and `median` functions accept an optional `estimation-strategy` keyword. This parameter determines how the quantile is estimated, particularly how interpolation is handled when the desired quantile falls between data points in the sorted sequence.

(utls/examples-note
  (stats/quantile mpg 0.25)
  (stats/quantiles mpg [0.1 0.25 0.5 0.75 0.9])
  (stats/percentile mpg 50.0)
  (stats/percentiles residual-sugar [10 25 50 75 90])
  (stats/median mpg)
  (stats/median-3 15 5 10))

;; The available `estimation-strategy` values are `:legacy` (default), `:r1`, `:r2`, `:r3`, `:r4`, `:r5`, `:r6`, `:r7`, `:r8` and `:r9`. Formulas for all of them can be found on this [Wikipedia article](https://en.wikipedia.org/wiki/Quantile#Estimating_quantiles_from_a_sample). `:legacy` uses estimate of the form: $Q_p = x_{\lceil p(N+1) - 1/2 \rceil}$

;; The plot below illustrates the differences between these estimation strategies for the sample `vs = [1 10 10 30]`.

(def vs [1 10 10 30])

(gg/->image (gg/functions [[":legacy" #(stats/quantile vs %)]
                           [":r1" #(stats/quantile vs % :r1)]
                           [":r2" #(stats/quantile vs % :r2)]
                           [":r3" #(stats/quantile vs % :r3)]
                           [":r4" #(stats/quantile vs % :r4)]
                           [":r5" #(stats/quantile vs % :r5)]
                           [":r6" #(stats/quantile vs % :r6)]
                           [":r7" #(stats/quantile vs % :r7)]
                           [":r8" #(stats/quantile vs % :r8)]
                           [":r9" #(stats/quantile vs % :r9)]]
                          {:x [m/MACHINE-EPSILON 1.0]}))

(utls/examples-note
  (stats/quantiles vs [0.1 0.25 0.5 0.75 0.9] :legacy)
  (stats/quantiles vs [0.1 0.25 0.5 0.75 0.9] :r1)
  (stats/quantiles vs [0.1 0.25 0.5 0.75 0.9] :r2)
  (stats/quantiles vs [0.1 0.25 0.5 0.75 0.9] :r3)
  (stats/quantiles vs [0.1 0.25 0.5 0.75 0.9] :r4)
  (stats/quantiles vs [0.1 0.25 0.5 0.75 0.9] :r5)
  (stats/quantiles vs [0.1 0.25 0.5 0.75 0.9] :r6)
  (stats/quantiles vs [0.1 0.25 0.5 0.75 0.9] :r7)
  (stats/quantiles vs [0.1 0.25 0.5 0.75 0.9] :r8)
  (stats/quantiles vs [0.1 0.25 0.5 0.75 0.9] :r9))

;; ### Weighted

;; There are also functions to calculate weighted quantiles and medians. These are useful when individual data points have different levels of importance or contribution.

;; *   `wquantile`: Calculates the q-th weighted quantile for a sequence `vs` with corresponding `weights`.
;; *   `wquantiles`: Calculates multiple weighted quantiles for a sequence `vs` with `weights`.
;; *   `wmedian`: Calculates the weighted median (0.5 weighted quantile) for `vs` with `weights`.

;; All these functions accept an optional `method` keyword argument that specifies the interpolation strategy when a quantile falls between points in the weighted empirical cumulative distribution function (ECDF). The available methods are:

;; *   `:linear` (Default): Performs linear interpolation between the data values corresponding to the cumulative weights surrounding the target quantile.
;; *   `:step`: Uses a step function (specifically, step-before interpolation) based on the weighted ECDF. The result is the data value whose cumulative weight range includes the target quantile.
;; *   `:average`: Computes the average of the step-before and step-after interpolation methods. This can be useful when a quantile corresponds exactly to a cumulative weight boundary.

;; Let's define a sample dataset and weights:

(def sample-data [10 15 30 50 100])
(def sample-weights [1 2 5 1 1])

(utls/examples-note
  ;; Calculate the 0.25, 0.5 (median), and 0.75 weighted quantiles
  (stats/wquantile sample-data sample-weights 0.25)
  (stats/wquantile sample-data sample-weights 0.5) ;; same as wmedian
  (stats/wmedian sample-data sample-weights)
  (stats/wquantile sample-data sample-weights 0.75)

  ;; Using different interpolation methods for the weighted median
  (stats/wmedian sample-data sample-weights :linear)
  (stats/wmedian sample-data sample-weights :step)
  (stats/wmedian sample-data sample-weights :average)

  ;; Calculate multiple weighted quantiles at once
  (stats/wquantiles sample-data sample-weights [0.2 0.5 0.8])
  (stats/wquantiles sample-data sample-weights [0.2 0.5 0.8] :step)
  (stats/wquantiles sample-data sample-weights [0.2 0.5 0.8] :average))

;; Using `mpg` data and `wt` (car weight) as weights:

(utls/examples-note
  (stats/wmedian mpg wt)
  (stats/wquantile mpg wt 0.25)
  (stats/wquantiles mpg wt [0.1 0.25 0.5 0.75 0.9] :average))

;; When weights are equal to 1.0, then:

;; * `:linear` method is the same as `:r4` estimation strategy in `quantiles`
;; * `:step` is the same as `:r1`
;; * `:average` has no corresponding strategy 

(utls/examples-note
  (stats/quantiles mpg [0.1 0.25 0.5 0.75 0.9] :r4)
  (stats/wquantiles mpg (repeat (count mpg) 1.0) [0.1 0.25 0.5 0.75 0.9])
  (stats/quantiles mpg [0.1 0.25 0.5 0.75 0.9] :r1)
  (stats/wquantiles mpg (repeat (count mpg) 1.0) [0.1 0.25 0.5 0.75 0.9] :step)
  (stats/wquantiles mpg (repeat (count mpg) 1.0) [0.1 0.25 0.5 0.75 0.9] :average))

;; ## Measures of Dispersion/Deviation

;; Statistics describing the spread or variability of data.

;; ::: {.callout-tip title="Defined functions"}
;; * `variance`, `population-variance`
;; * `stddev`, `population-stddev`
;; * `wvariance`, `population-wvariance`
;; * `wstddev`, `population-wstddev`
;; * `pooled-variance`, `pooled-stddev`
;; * `variation`, `l-variation`
;; * `mean-absolute-deviation`
;; * `median-absolute-deviation`, `mad`
;; * `pooled-mad`
;; * `sem`
;; :::

;; ### Variance and standard deviation

;; Variance and standard deviation are fundamental measures of the dispersion or spread of a dataset around its mean.

;; * **Variance** quantifies the average squared difference of each data point from the mean. A higher variance indicates that the data points are more spread out, while a lower variance indicates they are clustered closer to the mean.

;; * **Standard Deviation** is the square root of the variance. It is expressed in the same units as the data, making it more interpretable than variance as a measure of spread.


;; * **Sample Variance** (`variance`) and **Sample Standard Deviation** (`stddev`): These are estimates of the population variance and standard deviation, calculated from a sample of data. They use a denominator of $N-1$ (Bessel's correction) to provide an unbiased estimate of the population variance.

;; $$s^2 = \frac{\sum_{i=1}^{n} (x_i - \bar{x})^2}{n-1}$$

;; $$s = \sqrt{s^2}$$

;; Both functions can optionally accept a pre-calculated mean (`mu`) as a second argument.

;; * **Population Variance** (`population-variance`) and **Population Standard Deviation** (`population-stddev`): These are used when the data represents the entire population of interest, or when a biased estimate (maximum likelihood estimate) from a sample is desired. They use a denominator of $N$.

;; $$\sigma^2 = \frac{\sum_{i=1}^{N} (x_i - \mu)^2}{N}$$

;; $$\sigma = \sqrt{\sigma^2}$$

;; These also accept an optional pre-calculated population mean (`mu`).

;; * **Weighted Variance** (`wvariance`, `population-wvariance`) and **Weighted Standard Deviation** (`wstddev`, `population-wstddev`): These calculate variance and standard deviation when each data point has an associated weight.
;;     For weighted sample variance (unbiased form, where $w_i$ are weights):

;; $$\bar{x}_w = \frac{\sum w_i x_i}{\sum w_i}$$

;; $$s_w^2 = \frac{\sum w_i (x_i - \bar{x}_w)^2}{(\sum w_i) - 1}$$

;; For weighted population variance:

;; $$\sigma_w^2 = \frac{\sum w_i (x_i - \bar{x}_w)^2}{\sum w_i}$$

;; Weighted standard deviations are the square roots of their respective variances.

;; * **Pooled Variance** (`pooled-variance`) and **Pooled Standard Deviation** (`pooled-stddev`): These are used to estimate a common variance when data comes from several groups that are assumed to have the same population variance. Following methods can be used (where $k$ is the number of groups, each with $n_i$ number of observations and sample variance $s_i^2$)
;;     - `:unbiased` (default)
;; $$s_p^2 = \frac{\sum_{i=1}^{k} (n_i-1)s_i^2}{\sum_{i=1}^{k} n_i - k}$$
;;     - `:biased`
;; $$s_p^2 = \frac{\sum_{i=1}^{k} (n_i-1)s_i^2}{\sum_{i=1}^{k} n_i}$$
;;     - `:avg` - simple average of group variances.
;; $$s_p^2 = \frac{\sum_{i=1}^{k} s_i^2}{k}$$

;; Pooled standard deviation is $\sqrt{s_p^2}$.

(utls/examples-note
  (stats/variance mpg)
  (stats/stddev mpg)

  (stats/population-variance mpg)
  (stats/population-stddev mpg)

  (stats/variance mpg (stats/mean mpg))
  (stats/population-variance mpg (stats/mean mpg)))

;; Weighted variance and standard deviation

(utls/examples-note
  (stats/wvariance hp wt)
  (stats/wstddev hp wt)
  (stats/population-wvariance hp wt)
  (stats/population-wstddev hp wt))

;; Pooled variance and standard deviation

(utls/examples-note
  (stats/pooled-variance [setosa-sepal-length virginica-sepal-length])
  (stats/pooled-stddev [setosa-sepal-length virginica-sepal-length])
  (stats/pooled-variance [setosa-sepal-length virginica-sepal-length] :biased)
  (stats/pooled-variance [setosa-sepal-length virginica-sepal-length] :avg))

;; Beyond variance and standard deviation, we have three additional functions:
;;
;; * **Coefficient of Variation** (`variation`): This is a standardized measure of dispersion, calculated as the ratio of the standard deviation $s$ to the mean $\bar{x}$.

;; $$CV = \frac{s}{\bar{x}}$$

;; The CV is unitless, making it useful for comparing the variability of datasets with different means or units. It's most meaningful for data measured on a ratio scale (i.e., with a true zero point) and where all values are positive.
;;
;; * **Standard Error of the Mean** (`sem`): The SEM estimates the standard deviation of the sample mean if you were to draw multiple samples from the same population. It indicates how precisely the sample mean estimates the true population mean.

;; $$SEM = \frac{s}{\sqrt{n}}$$

;; where $s$ is the sample standard deviation and $n$ is the sample size. A smaller SEM suggests a more precise estimate of the population mean.

;; * **L-variation** (`l-variation`): Calculates the coefficient of L-variation. This is a dimensionless measure of dispersion, analogous to the coefficient of variation.

;; $$\tau_2 = \lambda_2 / \lambda_1$$

(utls/examples-note
  (stats/variation mpg)
  (stats/variation residual-sugar)
  (stats/l-variation mpg)
  (stats/l-variation residual-sugar)
  (stats/sem mpg)
  (stats/sem residual-sugar))

;; ### MAD

;; MAD typically refers to Median Absolute Deviation, a robust measure of statistical dispersion. `fastmath.stats` also provides the Mean Absolute Deviation.
;; 
;; * **Median Absolute Deviation** (`median-absolute-deviation` or `mad`): This is a robust measure of the variability of a univariate sample. It is defined as the median of the absolute deviations from the data's median.
;; 
;; $$MAD = \text{median}(|X_i - \text{median}(X)|)$$
;; 
;; If a specific center $c$ is provided, it's $MAD_c = \text{median}(|X_i - c|)$. Also, different estimation strategies can be used, see [[median]]
;; MAD is less sensitive to outliers than the standard deviation.
;; 
;; * **Mean Absolute Deviation** (`mean-absolute-deviation`): This measures variability as the average of the absolute deviations from a central point, typically the data's mean.
;; 
;; $$MeanAD = \frac{1}{n} \sum_{i=1}^{n} |X_i - \text{mean}(X)|$$
;; 
;; If a specific center $c$ is provided, it's $MeanAD_c = \frac{1}{n} \sum_{i=1}^{n} |X_i - c|$.
;; MeanAD is more sensitive to outliers than MAD but less sensitive than the standard deviation.
;; 
;; * **Pooled MAD** (`pooled-mad`): This function calculates a pooled estimate of the Median Absolute Deviation when data comes from several groups. For each group $i$, absolute deviations from its median $M_i$ are calculated: $Y_{ij} = |X_{ij} - M_i|$. The pooled MAD is then the median of all such $Y_{ij}$ values, scaled by a constant `const` (which defaults to approximately 1.4826, to make it comparable to the standard deviation for normal data).
;; 
;; $$PooledMAD = \text{const} \cdot \text{median}(\{Y_{ij} \mid \text{for all groups } i \text{ and observations } j \text{ in group } i\})$$

(utls/examples-note
  (stats/mad mpg)
  (stats/median-absolute-deviation mpg)
  (stats/median-absolute-deviation mpg (stats/median mpg) :r3)
  (stats/median-absolute-deviation mpg (stats/mean mpg))

  (stats/mean-absolute-deviation mpg)
  (stats/mean-absolute-deviation mpg (stats/median mpg))

  (stats/pooled-mad [setosa-sepal-length virginica-sepal-length])
  (stats/pooled-mad [setosa-sepal-length virginica-sepal-length] 1.0))

;; ## Moments and Shape

;; Moments and shape statistics describe the form of a dataset's distribution, particularly its symmetry and peakedness.
;;
;; ::: {.callout-tip title="Defined functions"}
;; * `moment`
;; * `skewness`
;; * `kurtosis`
;; * `l-moment`
;; :::

;; ### Conventional Moments (`moment`)

;; The `moment` function calculates statistical moments of a dataset. Moments can be central (around the mean), raw (around zero), or around a specified center. They can also be absolute and/or normalized.

;; * **k-th Central Moment**: $\mu_k = E[(X - \mu)^k] \approx \frac{1}{n} \sum_{i=1}^{n} (x_i - \bar{x})^k$. Calculated when `center` is `nil` (default) and `:mean?` is `true` (default).

;; * **k-th Raw Moment** (about origin): $\mu'_k = E[X^k] \approx \frac{1}{n} \sum_{i=1}^{n} x_i^k$. Calculated if `center` is `0.0`.

;; * **k-th Absolute Central Moment**: $E[|X - \mu|^k] \approx \frac{1}{n} \sum_{i=1}^{n} |x_i - \bar{x}|^k$. Calculated if `:absolute?` is `true`.

;; * **Normalization**: If `:normalize?` is `true`, the moment is divided by $\sigma^k$ (where $\sigma$ is the standard deviation), yielding a scale-invariant measure. For example, the 3rd normalized central moment is related to skewness, and the 4th to kurtosis.

;; * **Power of sum of differences**:  If `:mean?` is `false`, the function returns the sum $\sum (x_i - c)^k$ (or sum of absolute values) instead of the mean.

;; The `order` parameter specifies $k$. For example, the 2nd central moment for $k=2$ is the variance.

(utls/examples-note
  (stats/moment mpg 2) ;; Variance (almost, as it's mean of squared devs, not sample variance)
  (stats/variance mpg)
  (stats/moment mpg 3) ;; Related to skewness
  (stats/moment mpg 4 {:normalize? true}) ;; Normalized 4th central moment, related to kurtosis
  (stats/moment mpg 1 {:absolute? true :center (stats/median mpg)}) ;; Mean absolute deviation from median
  )

;; ### Skewness

;; Skewness measures the asymmetry of a probability distribution about its mean. `fastmath.stats/skewness` offers several types:

;; **Moment-based (sensitive to outliers):**

;; * `:G1` (Default): Sample skewness based on the 3rd standardized moment, adjusted for sample bias (via Apache Commons Math).

;; $$G_1 = \frac{n}{(n-1)(n-2)} \sum_{i=1}^n \left(\frac{x_i - \bar{x}}{s}\right)^3$$

;; * `:g1` or `:pearson`: Pearson's moment coefficient of skewness, another bias-adjusted version of the 3rd standardized central moment $m_3$.

;; $$g_1 = \frac{m_3}{m_2^{3/2}}$$

;; * `:b1`: Sample skewness coefficient, related to $g_1$.

;; $$b_1 = \frac{m_3}{s^3}$$

;; * `:skew`: Skewness used in BCa

;; $$SKEW = \frac{\sum_{i=1}^n (x_i - \bar{x})^3}{(\sum_{i=1}^n (x_i - \bar{x})^2)^{3/2}}  = \frac{g_1}{\sqrt{n}}$$

;; **Robust (less sensitive to outliers):**

;; * `:median`: Median Skewness / Pearson's first skewness coefficient.

;; $$S_P = 3 \frac{\text{mean} - \text{median}}{\text{stddev}}$$

;; * `:mode`: Pearson's second skewness coefficient. Mode estimation method can be specified.

;; $$S_K = \frac{\text{mean} - \text{mode}}{\text{stddev}}$$

;; * `:bowley` or `:yule` (with $u=0.25$): Based on quartiles $Q_1, Q_2, Q_3$.

;; $$S_B = \frac{(Q_3 - Q_2) - (Q_2 - Q_1)}{Q_3 - Q_1} = \frac{Q_3 + Q_1 - 2Q_2}{Q_3 - Q_1}$$

;; * `:yule` or `:B1` (Yule's coefficient): Generalization of Bowley's, using quantiles $Q_u, Q_{0.5}, Q_{1-u}$.

;; $$B_1 = S_Y(u) = \frac{(Q_{1-u} - Q_{0.5}) - (Q_{0.5} - Q_u)}{Q_{1-u} - Q_u}$$

;; * `:B3`: Robust measure by Groeneveld and Meeden

;; $$B_3 = \frac{\text{mean} - \text{median}}{\text{mean}(|X_i - \text{median}|)}$$

;; * `:hogg`: Based on comparing trimmed means ($U_{0.05}$: mean of top 5%, $L_{0.05}$: mean of bottom 5%, $M_{0.25}$: 25% trimmed mean).

;; $$S_H = \frac{U_{0.05} - M_{0.25}}{M_{0.25} - L_{0.05}}$$

;; * `:l-skewness`: L-moments based skewness.

;; $$\tau_3 = \lambda_3 / \lambda_2$$

;; Positive skewness indicates a tail on the right side of the distribution; negative skewness indicates a tail on the left. Zero indicates symmetry.

(utls/examples-note
  (stats/skewness residual-sugar)
  (stats/skewness residual-sugar :G1)
  (stats/skewness residual-sugar :g1)
  (stats/skewness residual-sugar :pearson)
  (stats/skewness residual-sugar :b1)
  (stats/skewness residual-sugar :skew)

  (stats/skewness residual-sugar :median)
  (stats/skewness residual-sugar :mode)
  (stats/skewness residual-sugar [:mode :histogram])
  (stats/skewness residual-sugar [:mode :kde])
  (stats/skewness residual-sugar :bowley)
  (stats/skewness residual-sugar :yule)
  (stats/skewness residual-sugar [:yule 0.1])
  (stats/skewness residual-sugar :B3)
  (stats/skewness residual-sugar :hogg)
  (stats/skewness residual-sugar :l-skewness))

;; Effect of an outlier is visible for moment based skewness, while has no effect when robust method is used.

(utls/examples-note
  (stats/skewness (conj residual-sugar -1000))
  (stats/skewness (conj residual-sugar -1000) :l-skewness))

;; ### Kurtosis

;; Kurtosis measures the "tailedness" or "peakedness" of a distribution. High kurtosis means heavy tails (more outliers) and a sharp peak (leptokurtic); low kurtosis means light tails and a flatter peak (platykurtic). `fastmath.stats/kurtosis` offers several types:

;; **Moment-based (sensitive to outliers):**

;; * `:G2` (Default): Sample kurtosis (Fisher's definition, not excess), adjusted for sample bias (via Apache Commons Math). For a normal distribution, this is approximately 3.

;; $$G_2 = \frac{(n+1)n}{(n-1)(n-2)(n-3)} \sum_{i=1}^n \left(\frac{x_i - \bar{x}}{s}\right)^4 - 3\frac{(n-1)^2}{(n-2)(n-3)}$$

;; * `:g2` or `:excess`: Sample excess kurtosis. For a normal distribution, this is approximately 0.

;; $$g_2 = \frac{m_4}{m_2^2}-3$$

;; *   `:kurt`: Kurtosis defined as $g_2 + 3$.

;; $$g_{kurt} = \frac{m_4}{m_2^2} = g_2 + 3$$

;; *   `:b2`: Sample kurtosis

;; $$b_2 = \frac{m_4}{s^4}-3$$

;; **Robust (less sensitive to outliers):**

;; * `:geary`: Geary's 'g' measure of kurtosis. Normal $\approx \sqrt{2/\pi} \approx 0.798$.

;; $$g = \frac{MeanAD}{\sigma^2}$$

;; * `:moors`: Based on octiles $E_i$ (quantiles $i/8$) and centered by subtracting $1.233$ (Moors' constant for normality).

;; $$M_0 = \frac{(E_7-E_5) + (E_3-E_1)}{E_6-E_2}-1.233$$

;; * `:crow` (Crow-Siddiqui): Based on quantiles $Q_\alpha, Q_{1-\alpha}, Q_\beta, Q_{1-\beta}$ and centered for normality. By default $\alpha=0.025$ and $\beta=0.25$.

;; $$CS(\alpha, \beta) = \frac{Q_{1-\alpha} - Q_{\alpha}}{Q_{1-\beta} - Q_{\beta}}-2.906$$

;; * `:hogg`: Based on trimmed means $U_p$ (mean of top $p\%$) and $L_p$ (mean of bottom $p\%$) and centered by subtracting $2.585$. By default $\alpha=0.005$ and $\beta=0.5$.

;; $$K_H(\alpha, \beta) = \frac{U_{\alpha} - L_{\alpha}}{U_{\beta} - L_{\beta}}-2.585$$

;; * `:l-kurtosis`: L-moments based kurtosis.

;; $$\tau_4 = \lambda_4 / \lambda_2$$

(utls/examples-note
  (stats/kurtosis residual-sugar)
  (stats/kurtosis residual-sugar :G2)
  (stats/kurtosis residual-sugar :g1)
  (stats/kurtosis residual-sugar :excess)
  (stats/kurtosis residual-sugar :kurt)
  (stats/kurtosis residual-sugar :b2)
  
  (stats/kurtosis residual-sugar :geary)
  (stats/kurtosis residual-sugar :moors)
  (stats/kurtosis residual-sugar :crow)
  (stats/kurtosis residual-sugar [:crow 0.05 0.25])
  (stats/kurtosis residual-sugar :hogg)
  (stats/kurtosis residual-sugar [:hogg 0.025 0.45])
  (stats/kurtosis residual-sugar :l-kurtosis))

;; Effect of an outlier is visible for moment based kurtosis, while has no effect when robust method is used.

(utls/examples-note
  (stats/kurtosis (conj residual-sugar -1000 1000))
  (stats/kurtosis (conj residual-sugar -1000 1000) :l-kurtosis))

;; ### L-moment

;; L-moments are summary statistics analogous to conventional moments but are computed from linear combinations of order statistics (sorted data). They are more robust to outliers and provide better estimates for small samples compared to conventional moments.

;; *   `l-moment vs order`: Calculates the L-moment of a specific `order`.
;;     *   $\lambda_1$: L-location (identical to the mean).
;;     *   $\lambda_2$: L-scale (a measure of dispersion).
;;     *   Higher orders relate to shape.
;; *   Trimmed L-moments (TL-moments) can be calculated by specifying `:s` (left trim) and `:t` (right trim) as number of trimmed samples
;; *   L-moment Ratios: If `:ratio? true`, normalized L-moments are returned.
;;     *   $\tau_3 = \lambda_3 / \lambda_2$: Coefficient of L-skewness (same as `(stats/skewness vs :l-skewness)`).
;;     *   $\tau_4 = \lambda_4 / \lambda_2$: Coefficient of L-kurtosis (same as `(stats/kurtosis vs :l-kurtosis)`).

;; L-moments often provide more reliable inferences about the underlying distribution shape, especially when data may contain outliers or come from heavy-tailed distributions.

(utls/examples-note
  (stats/l-moment mpg 1)
  (stats/mean mpg)
  (stats/l-moment mpg 2)
  (stats/l-moment mpg 3)

  (stats/l-moment residual-sugar 3)
  (stats/l-moment residual-sugar 3 {:s 10})
  (stats/l-moment residual-sugar 3 {:t 10})
  (stats/l-moment residual-sugar 3 {:s 10 :t 10}))

;; Relation to skewness and kurtosis

(utls/examples-note
  (stats/l-moment residual-sugar 3 {:ratio? true})
  (stats/skewness residual-sugar :l-skewness)
  (stats/l-moment residual-sugar 4 {:ratio? true})
  (stats/kurtosis residual-sugar :l-kurtosis))

;; ## Intervals and Extents

;; This section provides functions to describe the spread or define specific ranges and intervals within a dataset.

;; ::: {.callout-tip title="Defined functions"}
;; * `span`, `iqr`
;; * `extent`, 
;; * `stddev-extent`, `mad-extent`, `sem-extent`
;; * `percentile-extent`, `quantile-extent`
;; * `pi`, `pi-extent`
;; * `hpdi-extent`
;; * `adjacent-values`
;; * `inner-fence-extent`, `outer-fence-extent`
;; * `percentile-bc-extent`, `percentile-bca-extent`
;; * `ci`
;; :::

;; * **Basic Range:** Functions like `span` ($max - min$) and `extent` (providing $[min, max]$ and optionally the mean) offer simple measures of the total spread of the data.
;;
;; * **Interquartile Range:** `iqr` ($Q_3 - Q_1$) specifically measures the spread of the middle 50% of the data, providing a robust alternative to the total range.
;;
;; * **Symmetric Spread Intervals:** Functions ending in `-extent` such as `stddev-extent`, `mad-extent`, and `sem-extent` define intervals typically centered around the mean or median. They represent a range defined by adding/subtracting a multiple (usually 1) of a measure of dispersion (Standard Deviation, Median Absolute Deviation, or Standard Error of the Mean) from the central point.
;;
;; * **Quantile-Based Intervals:** `percentile-extent`, `quantile-extent`, `pi`, `pi-extent`, and `hpdi-extent` define intervals based on quantiles or percentiles of the data. These functions capture specific ranges containing a certain percentage of the data points (e.g., the middle 95% defined by quantiles 0.025 and 0.975). `hpdi-extent` calculates the shortest interval containing a given proportion of data, based on empirical density.
;;
;; * **Box Plot Boundaries:** `adjacent-values` (LAV, UAV) and fence functions (`inner-fence-extent`, `outer-fence-extent`) calculate specific bounds based on quartiles and multiples of the IQR. These are primarily used in box plot visualization and as a conventional method for identifying potential outliers.
;;
;; * **Confidence and Prediction Intervals:** `ci`, `percentile-bc-extent`, and `percentile-bca-extent` provide inferential intervals. `ci` estimates a confidence interval for the population mean using the t-distribution. `percentile-bc-extent` and `percentile-bca-extent` (Bias-Corrected and Bias-Corrected Accelerated) are advanced bootstrap methods for estimating confidence intervals for statistics, offering robustness against non-normality and bias.

;; Note that:

;; * $IQR = Q_3-Q_1$
;; * $LIF=Q_1-1.5 \times IQR$
;; * $UIF=Q_3+1.5 \times IQR$
;; * $LOF=Q_1-3\times IQR$
;; * $UOF=Q_3+3\times IQR$
;; * $CI=\bar{x} \pm t_{\alpha/2, n-1} \frac{s}{\sqrt{n}}$

(kind/table
 {:column-names ["Function" "Returned value"]
  :row-vectors (utls/wrap-into-markdown
                [["`span`" "$max-min$"]
                 ["`iqr`" "$Q_3-Q_1$"]
                 ["`extent`" "`[min, max, mean]` or `[min, max]` (when `:mean?` is `false`)"]
                 ["`stddev-extent`" "`[mean - stddev, mean + stddev, mean]`"]
                 ["`mad-extent`" "`[median - mad, median + mad, median]`"]
                 ["`sem-extent`" "`[mean - sem, mean + sem, mean]`"]
                 ["`percentile-extent`" "`[p1-val, p2-val, median]` with default `p1=25` and `p2=75`"]
                 ["`quantile-extent`" "`[q1-val, q2-val, median]` with default `q1=0.25` and `q2=0.75`"]
                 ["`pi`" "`{p1 p1-val p2 p2-val}` defined by `size=p2-p1`"]
                 ["`pi-extent`" "`[p1-val, p2-val, median]` defined by `size=p2-p1`"]
                 ["`hdpi-extent`" "`[p1-val, p2-val, median]` defined by `size=p2-p1`"]
                 ["`adjacent-values`" "`[LAV, UAV, median]`"]
                 ["`inner-fence-extent`" "`[LIF, UIF, median]`"]
                 ["`outer-fence-extent`" "`[LOF, UOF, median]`"]
                 ["`ci`" "`[lower upper mean]`"]
                 ["`percentile-bc-extent`" "`[lower upper mean]`"]
                 ["`percentile-bca-extent`" "`[lower upper mean]`"]])})

(utls/examples-note
  (stats/span mpg)
  (stats/iqr mpg)
  (stats/extent mpg)
  (stats/extent mpg false)
  (stats/stddev-extent mpg)
  (stats/mad-extent mpg)
  (stats/sem-extent mpg)
  (stats/percentile-extent mpg)
  (stats/percentile-extent mpg 2.5 97.5)
  (stats/percentile-extent mpg 2.5 97.5 :r9)
  (stats/quantile-extent mpg)
  (stats/quantile-extent mpg 0.025 0.975)
  (stats/quantile-extent mpg 0.025 0.975 :r9)
  (stats/pi mpg 0.95)
  (stats/pi-extent mpg 0.95)
  (stats/hpdi-extent mpg 0.95)
  (stats/adjacent-values mpg)
  (stats/inner-fence-extent mpg)
  (stats/outer-fence-extent mpg)
  (stats/ci mpg)
  (stats/ci mpg 0.1)
  (stats/percentile-bc-extent mpg)
  (stats/percentile-bc-extent mpg 10.0)
  (stats/percentile-bca-extent mpg)
  (stats/percentile-bca-extent mpg 10.0))

(gg/->image (gg/errorbars [[:extent (stats/extent mpg)]
                           [:stddev (stats/stddev-extent mpg)]
                           [:mad (stats/mad-extent mpg)]
                           [:sem (stats/sem-extent mpg)]
                           [:quantile (stats/quantile-extent mpg)]
                           [:pi-0.5 (stats/pi-extent mpg)]
                           [:pi-0.95 (stats/pi-extent mpg 0.95)]
                           [:hpdi-0.5 (stats/hpdi-extent mpg 0.5)]
                           [:hpdi-0.95 (stats/hpdi-extent mpg)]
                           [:adjacent-vals (stats/adjacent-values mpg)]
                           [:inner-fence (stats/inner-fence-extent mpg)]
                           [:outer-fence (stats/outer-fence-extent mpg)]
                           [:ci-0.05 (stats/ci mpg)]
                           [:ci-0.1 (stats/ci mpg 0.1)]
                           [:bc (stats/percentile-bc-extent mpg)]
                           [:bca (stats/percentile-bca-extent mpg)]]
                          {:ylab "Extents"
                           :title "Various extents for MPG"}))

;; ## Outlier Detection

;; Outlier detection involves identifying data points that are significantly different from other observations. Outliers can distort statistical analyses and require careful handling. `fastmath.stats` provides functions to find and optionally remove such values based on the Interquartile Range (IQR) method.

;; ::: {.callout-tip title="Defined functions"}
;; * `outliers`
;; :::

;; `outliers` function use the inner fence rule based on the IQR and returns a sequence containing only the data points identified as outliers.

;; *   **Lower Inner Fence (LIF):** $Q_1 - 1.5 \times IQR$
;; *   **Upper Inner Fence (UIF):** $Q_3 + 1.5 \times IQR$

;; Where $Q_1$ is the first quartile (25th percentile) and $Q_3$ is the third quartile (75th percentile). Points falling below the LIF or above the UIF are considered outliers.

;; Function accepts an optional `estimation-strategy` keyword (see [[quantile]]) to control how quartiles are calculated, which affects the fence boundaries.

;; Let's find the outliers in the `residual-sugar` data.

(stats/outliers residual-sugar)

;; ## Data Transformation

;; Functions to modify data (scaling, normalizing, transforming).

;; ::: {.callout-tip title="Defined functions"}
;; * `standardize`, `robust-standardize`, `demean`
;; * `rescale`
;; * `remove-outliers`
;; * `trim`, `trim-lower`, `trim-upper`, `winsor`
;; * `box-cox-infer-lambda`, `box-cox-transformation`
;; * `yeo-johnson-infer-lambda`, `yeo-johnson-transformation`
;; :::

;; Data transformations are often necessary preprocessing steps in statistical analysis and machine learning. They can help meet the assumptions of certain models (e.g., normality, constant variance), improve interpretability, or reduce the influence of outliers. `fastmath.stats` offers several functions for these purposes, broadly categorized into linear scaling/centering, outlier handling, and power transformations for normality.

;; Let's demonstrate some of these transformations using the `residual-sugar` data from the wine quality dataset.

(utls/examples-note
  (stats/mean residual-sugar)
  (stats/stddev residual-sugar)
  (stats/median residual-sugar)
  (stats/mad residual-sugar)
  (stats/extent residual-sugar false)
  (count residual-sugar))

(gg/->image (gg/density residual-sugar
                        {:title "Residual Sugar dataset"}))

;; **Linear Transformations:** `standardize`, `robust-standardize`, `demean`, and `rescale` linearly transform data, preserving its shape but changing its location and/or scale.

;; *   `demean` centers the data by subtracting the mean, resulting in a dataset with a mean of zero.
;; *   `standardize` scales the demeaned data by dividing by the standard deviation, resulting in data with mean zero and standard deviation one (z-score normalization). This makes the scale of different features comparable.
;; *   `robust-standardize` provides a version less sensitive to outliers by centering around the median and scaling by the Median Absolute Deviation (MAD) or a quantile range (like the IQR).
;; *   `rescale` linearly maps the data to a specific target range (e.g., [0, 1]), useful for algorithms sensitive to input scale.

(def residual-sugar-demeaned (-> residual-sugar stats/demean))
(def residual-sugar-standardized (-> residual-sugar stats/standardize))
(def residual-sugar-robust-standardized (-> residual-sugar stats/robust-standardize))
(def residual-sugar-rescaled (-> residual-sugar stats/rescale))

(utls/examples-note
  (stats/mean residual-sugar-demeaned)
  (stats/stddev residual-sugar-demeaned)
  (stats/median residual-sugar-demeaned)
  (stats/mad residual-sugar-demeaned)
  (stats/extent residual-sugar-demeaned false)
  
  (stats/mean residual-sugar-standardized)
  (stats/stddev residual-sugar-standardized)
  (stats/median residual-sugar-standardized)
  (stats/mad residual-sugar-standardized)
  (stats/extent residual-sugar-standardized false)
  
  (stats/mean residual-sugar-robust-standardized)
  (stats/stddev residual-sugar-robust-standardized)
  (stats/median residual-sugar-robust-standardized)
  (stats/mad residual-sugar-robust-standardized)
  (stats/extent residual-sugar-robust-standardized false)
  
  (stats/mean residual-sugar-rescaled)
  (stats/stddev residual-sugar-rescaled)
  (stats/median residual-sugar-rescaled)
  (stats/mad residual-sugar-rescaled)
  (stats/extent residual-sugar-rescaled false))

;; **Outlier Handling:** `remove-outliers`, `trim`, `trim-lower`, `trim-upper`, and `winsor` address outliers based on quantile fences.

;; * `remove-outliers` returns a sequence containing the data points from the original sequence *excluding* those identified as outliers.
;; * `trim` removes values outside a specified quantile range (defaulting to 0.2 quantile, removing the bottom and top 20%). `trim-lower` and `trim-upper` remove only below or above a single quantile.
;; * `winsor` caps values outside a quantile range to the boundary values instead of removing them. This retains the sample size but reduces the influence of extreme values.

(def residual-sugar-no-outliers (stats/remove-outliers residual-sugar))
(def residual-sugar-trimmed (stats/trim residual-sugar))
(def residual-sugar-winsorized (stats/winsor residual-sugar))

(utls/examples-note
  (stats/mean residual-sugar-no-outliers)
  (stats/stddev residual-sugar-no-outliers)
  (stats/median residual-sugar-no-outliers)
  (stats/mad residual-sugar-no-outliers)
  (stats/extent residual-sugar-no-outliers false)
  (count residual-sugar-no-outliers)

  (stats/mean residual-sugar-trimmed)
  (stats/stddev residual-sugar-trimmed)
  (stats/median residual-sugar-trimmed)
  (stats/mad residual-sugar-trimmed)
  (stats/extent residual-sugar-trimmed false)
  (count residual-sugar-trimmed)

  (stats/mean residual-sugar-winsorized)
  (stats/stddev residual-sugar-winsorized)
  (stats/median residual-sugar-winsorized)
  (stats/mad residual-sugar-winsorized)
  (stats/quantiles residual-sugar [0.2 0.8])
  (stats/extent residual-sugar-winsorized false)
  (count residual-sugar-winsorized))

(kind/table
 [[(gg/->image (gg/density residual-sugar-no-outliers
                           {:title "Residual Sugar without outliers"}))
   (gg/->image (gg/density residual-sugar-trimmed
                           {:title "Residual Sugar trimmed"}))
   (gg/->image (gg/density residual-sugar-winsorized
                           {:title "Residual Sugar winsorized"}))]])

;; **Power Transformations:** `box-cox-transformation` and `yeo-johnson-transformation` (and their `infer-lambda` counterparts) are non-linear transformations that can change the shape of the distribution to be more symmetric or normally distributed. They are particularly useful for data that is skewed or violates assumptions of linear models. Both are invertable.

;; *   `box-cox-transformation` works in general for strictly positive data. It includes the log transformation as a special case (when lambda is $0.0$) and generalizes square root, reciprocal, and other power transformations. `box-cox-infer-lambda` helps find the optimal lambda parameter. Optional parameters:
;;     * `:negative?` (default: `false`), when set to `true` specific transformation is performed to keep information about sign.
;;     * `:scaled?` (default: `false`), when set to `true`, scale data by geometric mean, when is a number, this number is used as a scale.

;; $$y_{BC}^{(\lambda)}=\begin{cases}
;; \frac{y^\lambda-1}{\lambda} & \lambda\neq 0 \\
;; \log(y) & \lambda = 0
;; \end{cases}$$

;; Scaled version, with default scale set to geometric mean (GM):

;; $$y_{BC}^{(\lambda, s)}=\begin{cases}
;; \frac{y^\lambda-1}{\lambda s^{\lambda - 1}} & \lambda\neq 0 \\
;; s\log(y) & \lambda = 0
;; \end{cases}$$

;; When `:negative?` is set to true, formula takes the following form:

;; $$y_{BCneg}^{(\lambda)}=\begin{cases}
;; \frac{\operatorname{sgn}(y)|y|^\lambda-1}{\lambda} & \lambda\neq 0 \\
;; \operatorname{sgn}(y)\log(|y|+1) & \lambda = 0
;; \end{cases}$$

;; *   `yeo-johnson-transformation` extends Box-Cox to handle zero and negative values. `yeo-johnson-infer-lambda` finds the optimal lambda for this transformation.

;; $$y_{YJ}^{(\lambda)}=\begin{cases}
;; \frac{(y+1)^\lambda - 1}{\lambda} & \lambda \neq 0, y\geq 0 \\
;; \log(y+1) & \lambda = 0, y\geq 0 \\
;; \frac{(1-y)^{2-\lambda} - 1}{\lambda - 2} & \lambda \neq 2, y\geq 0 \\
;; -\log(1-y) & \lambda = 2, y\geq 0  
;; \end{cases}$$

;; Both fuctions accept additional parameters:

;; * `:alpha` (dafault: `0.0`): perform dataset shift by value of the `:alpha` before transformation.
;; * `:inversed?` (default: `false`): perform inverse transformation for given `lambda`. 

;; When `lambda` is set to `nil` optimal lambda will be calculated (only when `:inversed?` is `false`).

(utls/examples-note
  (stats/box-cox-transformation [0 1 10] 0.0)
  (stats/box-cox-transformation [0 1 10] 2.0)
  (stats/box-cox-transformation [0 1 10] -2.0 {:alpha 2})
  (stats/box-cox-transformation [0.375 0.444 0.497] -2.0 {:alpha 2 :inverse? true})
  (stats/box-cox-transformation [0 1 10] nil {:alpha 1})
  (stats/box-cox-transformation [0 1 10] nil {:scaled? true :alpha 1})
  (stats/box-cox-transformation [0 1 10] nil {:alpha -5 :negative? true})
  (stats/box-cox-transformation [0 1 10] 2.0 {:alpha -5 :negative? true :scaled? 2})
  (stats/box-cox-transformation [-6.5 -4.25 6.0] 2.0 {:alpha -5 :negative? true :scaled? 2 :inverse? true})
  (stats/yeo-johnson-transformation [0 1 10])
  (stats/yeo-johnson-transformation [0 1 10] 0.0)
  (stats/yeo-johnson-transformation [0 1 10] 2.0 {:alpha -5})
  (stats/yeo-johnson-transformation [-1.79 -1.61 17.5] 2.0 {:alpha -5 :inverse? true}))

;; Let's illustrate how real data look after transformation. We'll start with finding an optimal lambda parameter for both transformations.

(stats/box-cox-infer-lambda residual-sugar)
(stats/yeo-johnson-infer-lambda residual-sugar)

(def residual-sugar-box-cox (stats/box-cox-transformation residual-sugar nil))
(def residual-sugar-yeo-johnson (stats/yeo-johnson-transformation residual-sugar nil))

(utls/examples-note
  (stats/mean residual-sugar-box-cox)
  (stats/stddev residual-sugar-box-cox)
  (stats/median residual-sugar-box-cox)
  (stats/mad residual-sugar-box-cox)
  (stats/extent residual-sugar-box-cox false)

  (stats/mean residual-sugar-yeo-johnson)
  (stats/stddev residual-sugar-yeo-johnson)
  (stats/median residual-sugar-yeo-johnson)
  (stats/mad residual-sugar-yeo-johnson)
  (stats/extent residual-sugar-yeo-johnson false))

(kind/table
 [[(gg/->image (gg/density residual-sugar-box-cox
                           {:title "Residual Sugar Box-Cox, optimal lambda"}))
   (gg/->image (gg/density residual-sugar-yeo-johnson
                           {:title "Residual Sugar Yeo-Johnson, optimal lambda"}))]
  [(gg/->image (gg/density (stats/box-cox-transformation residual-sugar -1.0)
                           {:title "Residual Sugar Box-Cox, lambda=-1.0"}))
   (gg/->image (gg/density (stats/yeo-johnson-transformation residual-sugar -1.0)
                           {:title "Residual Sugar Yeo-Johnson, lambda=-1.0"}))]
  [(gg/->image (gg/density (stats/box-cox-transformation residual-sugar 0.0)
                           {:title "Residual Sugar Box-Cox, lambda=0.0"}))
   (gg/->image (gg/density (stats/yeo-johnson-transformation residual-sugar 0.0)
                           {:title "Residual Sugar Yeo-Johnson, lambda=0.0"}))]])

;; As you can see, the Yeo-Johnson transformation with the inferred lambda has made the `residual-sugar` distribution appear more symmetric and perhaps closer to a normal distribution shape.

;; Both power transformation can work on negative data as well. When Box-Cox is used, `:negative?` option should be set to `true`.

(stats/box-cox-infer-lambda residual-sugar
                            nil {:alpha (- (stats/mean residual-sugar)) :negative? true})

(stats/yeo-johnson-infer-lambda residual-sugar nil {:alpha (- (stats/mean residual-sugar))})

(def residual-sugar-box-cox-demeaned
  (stats/box-cox-transformation
   residual-sugar nil {:alpha (- (stats/mean residual-sugar)) :negative? true}))

(def residual-sugar-yeo-johnson-demeaned (stats/yeo-johnson-transformation
                                        residual-sugar nil
                                        {:alpha (- (stats/mean residual-sugar))}))

(utls/examples-note
  (stats/mean residual-sugar-box-cox-demeaned)
  (stats/stddev residual-sugar-box-cox-demeaned)
  (stats/median residual-sugar-box-cox-demeaned)
  (stats/mad residual-sugar-box-cox-demeaned)
  (stats/extent residual-sugar-box-cox-demeaned false)

  (stats/mean residual-sugar-yeo-johnson-demeaned)
  (stats/stddev residual-sugar-yeo-johnson-demeaned)
  (stats/median residual-sugar-yeo-johnson-demeaned)
  (stats/mad residual-sugar-yeo-johnson-demeaned)
  (stats/extent residual-sugar-yeo-johnson-demeaned false))

(kind/table
 [[(gg/->image (gg/density residual-sugar-box-cox-demeaned
                           {:title "Residual Sugar demeaned Box-Cox, optimal lambda"}))
   (gg/->image (gg/density residual-sugar-yeo-johnson-demeaned
                           {:title "Residual Sugar demeaned Yeo-Johnson, optimal lambda"}))]])

;; ## Correlation and Covariance

;; Measures of the relationship between two or more variables.

;; ::: {.callout-tip title="Defined functions"}
;; * `covariance`, `correlation`
;; * `pearson-correlation`, `spearman-correlation`, `kendall-correlation`
;; * `coefficient-matrix`, `correlation-matrix`, `covariance-matrix`
;; :::

;; **Covariance vs. Correlation:**

;; * `covariance` measures the extent to which two variables change together. A positive covariance means they tend to increase or decrease simultaneously. A negative covariance means one tends to increase when the other decreases. A covariance near zero suggests no linear relationship. The magnitude of covariance depends on the scales of the variables, making it difficult to compare covariances between different pairs of variables.
;; The sample covariance between two sequences $X = \{x_1, \dots, x_n\}$ and $Y = \{y_1, \dots, y_n\}$ is calculated as:

;; $$ \text{Cov}(X, Y) = \frac{1}{n-1} \sum_{i=1}^n (x_i - \bar{x})(y_i - \bar{y}) $$

;; where $\bar{x}$ and $\bar{y}$ are the sample means.

;; * `correlation` standardizes the covariance, resulting in a unitless measure that ranges from -1 to +1. It indicates both the direction and strength of a relationship. A correlation of +1 indicates a perfect positive relationship, -1 a perfect negative relationship, and 0 no linear relationship. The `correlation` function in `fastmath.stats` defaults to computing the Pearson correlation coefficient.

;; **Types of Correlation:**

;; *   `pearson-correlation`: The most common correlation coefficient, also known as the Pearson product-moment correlation coefficient ($r$). It measures the strength and direction of a *linear* relationship between two continuous variables. It assumes the variables are approximately normally distributed and that the relationship is linear. It is sensitive to outliers.
;;         The formula for the sample Pearson correlation coefficient is:

;; $$r = \frac{\text{Cov}(X, Y)}{s_x s_y} = \frac{\sum_{i=1}^n (x_i - \bar{x})(y_i - \bar{y})}{\sqrt{\sum_{i=1}^n (x_i - \bar{x})^2 \sum_{i=1}^n (y_i - \bar{y})^2}}$$

;; where $s_x$ and $s_y$ are the sample standard deviations.

;; *   `spearman-correlation`: Spearman's rank correlation coefficient ($\rho$) is a non-parametric measure of the strength and direction of a *monotonic* relationship between two variables. A monotonic relationship is one that is either consistently increasing or consistently decreasing, but not necessarily linear. Spearman's correlation is calculated by applying the Pearson formula to the *ranks* of the data values rather than the values themselves. This makes it less sensitive to outliers than Pearson correlation and suitable for ordinal data or when the relationship is monotonic but non-linear.
;; *   `kendall-correlation`: Kendall's Tau rank correlation coefficient ($\tau$) is another non-parametric measure of the strength and direction of a *monotonic* relationship. It is based on the number of concordant and discordant pairs of observations. A pair of data points is concordant if their values move in the same direction (both increase or both decrease) relative to each other, and discordant if they move in opposite directions. Kendall's Tau is generally preferred over Spearman's Rho for smaller sample sizes or when there are many tied ranks.
;; One common formulation, Kendall's Tau-a, is:

;; $$\tau_A = \frac{N_c - N_d}{n(n-1)/2}$$

;; where $N_c$ is the number of concordant pairs and $N_d$ is the number of discordant pairs.

;; **Comparison of Correlation Methods:**

;; *   Use **Pearson** for measuring *linear* relationships between continuous, normally distributed variables.
;; *   Use **Spearman** or **Kendall** for measuring *monotonic* relationships (linear or non-linear) between variables, especially when data is not normally distributed, contains outliers, or is ordinal. Kendall is often more robust with ties and smaller samples.
;; 
;; **Matrix Functions for Multiple Variables:**

;; *   `coefficient-matrix`: A generic function that computes a specified pairwise measure (defined by a function passed as an argument) between all pairs of sequences in a collection. Useful for generating matrices of custom similarity, distance, or correlation measures.
;; *   `covariance-matrix`: A specialization that computes the pairwise `covariance` for all sequences in a collection. The output is a symmetric matrix where the element at row `i`, column `j` is the covariance between sequence `i` and sequence `j`. The diagonal elements are the variances of the individual sequences.
;; *   `correlation-matrix`: A specialization that computes the pairwise `correlation` (Pearson by default, or specified via keyword like `:spearman` or `:kendall`) for all sequences in a collection. The output is a symmetric matrix where the element at row `i`, column `j` is the correlation between sequence `i` and sequence `j`. The diagonal elements are always 1.0 (a variable is perfectly correlated with itself).

;; Let's examine the correlations between the numerical features in the `iris` dataset.

(utls/examples-note
  (stats/covariance virginica-sepal-length setosa-sepal-length)
  (stats/correlation virginica-sepal-length setosa-sepal-length)
  (stats/pearson-correlation virginica-sepal-length setosa-sepal-length)
  (stats/spearman-correlation virginica-sepal-length setosa-sepal-length)
  (stats/kendall-correlation virginica-sepal-length setosa-sepal-length))

;; To generate matrices we'll use three sepal lengths samples. The last two examples use custom measure function: Euclidean distance between samples and Glass' delta.

(utls/examples-note
  (stats/covariance-matrix (vals sepal-lengths))
  (stats/correlation-matrix (vals sepal-lengths))
  (stats/correlation-matrix (vals sepal-lengths) :kendall)
  (stats/correlation-matrix (vals sepal-lengths) :spearman)
  (stats/coefficient-matrix (vals sepal-lengths) stats/L2 true)
  (stats/coefficient-matrix (vals sepal-lengths) stats/glass-delta))

;; ## Distance and Similarity Metrics

;; Measures of distance, error, or similarity between sequences or distributions.

;; ::: {.callout-tip title="Defined functions"}
;; * `me`, `mae`, `mape` 
;; * `rss`, `mse`, `rmse`
;; * `r2`
;; * `count=`, `L0`, `L1`, `L2sq`, `L2`, `LInf`
;; * `psnr`
;; * `dissimilarity`, `similarity`
;; :::

;; Distance metrics quantify how *far apart* or *different* two data sequences or probability distributions are. Similarity metrics, conversely, measure how *close* or *alike* they are, often being the inverse or a transformation of a distance. `fastmath.stats` provides a range of these measures suitable for comparing numerical sequences, observed counts (histograms), or theoretical probability distributions.

;; ### Error Metrics
;;
;; These functions typically quantify the difference between an observed sequence and a predicted or reference sequence, focusing on the magnitude of errors. All can accept a constant as a second argument.
;;
;; *   `me` (Mean Error): The average of the differences between corresponding elements.
;;     $$ ME = \frac{1}{n} \sum_{i=1}^n (x_i - y_i) $$
;; *   `mae` (Mean Absolute Error): The average of the absolute differences. More robust to outliers than squared error.
;;     $$ MAE = \frac{1}{n} \sum_{i=1}^n |x_i - y_i| $$
;; *   `mape` (Mean Absolute Percentage Error): The average of the absolute percentage errors. Useful for relative error assessment, but undefined if the reference value $x_i$ is zero.
;;     $$ MAPE = \frac{1}{n} \sum_{i=1}^n \left| \frac{x_i - y_i}{x_i} \right| \times 100\% $$
;; *   `rss` (Residual Sum of Squares): The sum of the squared differences. Used in least squares regression.
;;     $$ RSS = \sum_{i=1}^n (x_i - y_i)^2 $$
;; *   `mse` (Mean Squared Error): The average of the squared differences. Penalizes larger errors more heavily.
;;     $$ MSE = \frac{1}{n} \sum_{i=1}^n (x_i - y_i)^2 $$
;; *   `rmse` (Root Mean Squared Error): The square root of the MSE. Has the same units as the original data.
;;     $$ RMSE = \sqrt{\frac{1}{n} \sum_{i=1}^n (x_i - y_i)^2} $$
;; *   `r2` (Coefficient of Determination): Measures the proportion of the variance in the dependent variable that is predictable from the independent variable(s). Calculated as $1 - (RSS / TSS)$, where TSS is the Total Sum of Squares. It ranges from 0 to 1 for linear regression.
;;     $$ R^2 = 1 - \frac{\sum (x_i - y_i)^2}{\sum (x_i - \bar{x})^2} $$
;; * adjusted `r2`: A modified version of $R^2$ that has been adjusted for the number of predictors in the model. It increases only if the new term improves the model more than would be expected by chance.
;; $$ R^2_{adj} = 1 - (1 - R^2) \frac{n-1}{n-p-1} $$

;; Let's use `setosa-sepal-length` as observed and `virginica-sepal-length` as predicted (though they are independent samples, not predictions) to illustrate error measures.

(utls/examples-note
  (stats/me setosa-sepal-length virginica-sepal-length)
  (stats/mae setosa-sepal-length virginica-sepal-length)
  (stats/mape setosa-sepal-length virginica-sepal-length)
  (stats/rss setosa-sepal-length virginica-sepal-length)
  (stats/mse setosa-sepal-length virginica-sepal-length)
  (stats/rmse setosa-sepal-length virginica-sepal-length)
  (stats/r2 setosa-sepal-length virginica-sepal-length)
  (stats/r2 setosa-sepal-length virginica-sepal-length 2)
  (stats/r2 setosa-sepal-length virginica-sepal-length 5))

;; Also we can compare an observed sequence to a constant value. For example to a mean of the virginica sepal length.

(def vsl-mean (stats/mean virginica-sepal-length))

(utls/examples-note
  (stats/me setosa-sepal-length vsl-mean)
  (stats/mae setosa-sepal-length vsl-mean)
  (stats/mape setosa-sepal-length vsl-mean)
  (stats/rss setosa-sepal-length vsl-mean)
  (stats/mse setosa-sepal-length vsl-mean)
  (stats/rmse setosa-sepal-length vsl-mean)
  (stats/r2 setosa-sepal-length vsl-mean)
  (stats/r2 setosa-sepal-length vsl-mean 2)
  (stats/r2 setosa-sepal-length vsl-mean 5))

;; ### Distance Metrics (L-p Norms and others)
;;
;; These functions represent common distance measures, often related to L-p norms between vectors (sequences).
;;
;; *   `count=`, `L0`: Counts the number of elements that are *equal* in both sequences. While related to the L0 "norm" (which counts non-zero elements), this implementation counts *equal* elements after subtraction.
;;     $$ Count= = \sum_{i=1}^n \mathbb{I}(x_i = y_i) $$
;; *   `L1` (Manhattan/City Block Distance): The sum of the absolute differences.
;;     $$ L_1 = \sum_{i=1}^n |x_i - y_i| $$
;; *   `L2sq` (Squared Euclidean Distance): The sum of the squared differences. Equivalent to `rss`.
;;     $$ L_2^2 = \sum_{i=1}^n (x_i - y_i)^2 $$
;; *   `L2` (Euclidean Distance): The square root of the sum of the squared differences. The most common distance metric.
;;     $$ L_2 = \sqrt{\sum_{i=1}^n (x_i - y_i)^2} $$
;; *   `LInf` (Chebyshev Distance): The maximum absolute difference between corresponding elements.
;;     $$ L_\infty = \max_{i} |x_i - y_i| $$
;; *   `psnr` (Peak Signal-to-Noise Ratio): A measure of signal quality often used in image processing, derived from the MSE. Higher PSNR indicates better quality (less distortion). Calculated based on the maximum possible value of the data and the MSE.
;;     $$ PSNR = 10 \cdot \log_{10} \left( \frac{MAX^2}{MSE} \right) $$
;;
;; Using the sepal length samples again:
;;
(utls/examples-note
  (stats/count= setosa-sepal-length virginica-sepal-length)
  (stats/L0 setosa-sepal-length virginica-sepal-length)
  (stats/L1 setosa-sepal-length virginica-sepal-length)
  (stats/L2sq setosa-sepal-length virginica-sepal-length)
  (stats/L2 setosa-sepal-length virginica-sepal-length)
  (stats/LInf setosa-sepal-length virginica-sepal-length)
  (stats/psnr setosa-sepal-length virginica-sepal-length))
;;
;; ### Dissimilarity and Similarity

;; `dissimilarity` and `similarity` functions provide measures for comparing probability distributions or frequency counts (like histograms). They quantify how 'far apart' or 'alike' two data sequences, interpreted as distributions, are. They take a method keyword specifying the desired measure. Many methods exist, each with different properties and interpretations. They can accept raw data sequences, automatically creating histograms for comparison (controlled by `:bins`), or they can take pre-calculated frequency sequences or a data sequence and a `fastmath.random` distribution object.

;; Parameters:

;; *   `method` - The specific distance or similarity method to use.
;; *   `P-observed` - Frequencies, probabilities, or raw data (when `Q-expected` is a distribution or `:bins` is set).
;; *   `Q-expected` - Frequencies, probabilities, or a distribution object (when `P-observed` is raw data or `:bins` is set).
;; *   `opts` (map, optional) - Configuration options, including:
;;     *   `:probabilities?` (boolean, default: `true`): If `true`, input sequences are normalized to sum to 1.0 before calculating the measure, treating them as probability distributions.
;;     *   `:epsilon` (double, default: `1.0e-6`): A small number used to replace zero values in denominators or logarithms to avoid division-by-zero or log-of-zero errors.
;;     *   `:log-base` (double, default: `m/E`): The base for logarithms in information-theoretic measures.
;;     *   `:power` (double, default: `2.0`): The exponent for the `:minkowski` distance.
;;     *   `:remove-zeros?` (boolean, default: `false`): Removes pairs where both `P` and `Q` are zero before calculation.
;;     *   `:bins` (number, keyword, or seq): Used for comparisons involving raw data or distributions. Specifies the number of histogram bins, an estimation method (see [[histogram]]), or explicit bin edges for histogram creation.

;; #### Dissimilarity Methods

;; Higher values generally indicate greater difference.

;; *   **L-p Norms and Related:**
;;     *   `:euclidean`: Euclidean distance ($L_2$ norm) between the frequency/probability vectors.
;;         $$ D(P, Q) = \sqrt{\sum_i (P_i - Q_i)^2} $$
;;     *   `:city-block` / `:manhattan`: Manhattan distance ($L_1$ norm).
;;         $$ D(P, Q) = \sum_i |P_i - Q_i| $$
;;     *   `:chebyshev`: Chebyshev distance ($L_\infty$ norm), the maximum absolute difference.
;;         $$ D(P, Q) = \max_i |P_i - Q_i| $$
;;     *   `:minkowski`: Minkowski distance (generalized $L_p$ norm, controlled by `:power`).
;;         $$ D(P, Q) = \left(\sum_i |P_i - Q_i|^p\right)^{1/p} $$
;;     *   `:euclidean-sq` / `:squared-euclidean`: Squared Euclidean distance.
;;         $$ D(P, Q) = \sum_i (P_i - Q_i)^2 $$
;;     *   `:squared-chord`: Squared chord distance, related to Hellinger distance.
;;         $$ D(P, Q) = \sum_i (\sqrt{P_i} - \sqrt{Q_i})^2 $$

;; *   **Set-based/Overlap:** Measures derived from the concept of set overlap applied to frequencies/probabilities.
;;     *   `:sorensen`: Sorensen-Dice dissimilarity (1 - Dice similarity).
;;         $$ D(P, Q) = \frac{\sum_i |P_i - Q_i|}{\sum_i (P_i + Q_i)} $$
;;     *   `:gower`: Gower distance (average Manhattan distance).
;;         $$ D(P, Q) = \frac{1}{N} \sum_i |P_i - Q_i| $$
;;     *   `:soergel`: Soergel distance (1 - Jaccard similarity).
;;         $$ D(P, Q) = \frac{\sum_i |P_i - Q_i|}{\sum_i \max(P_i, Q_i)} $$
;;     *   `:kulczynski`: Kulczynski dissimilarity (1 - Kulczynski similarity, can be > 1).
;;         $$ D(P, Q) = \frac{\sum_i |P_i - Q_i|}{\sum_i \min(P_i, Q_i)} $$
;;     *   `:canberra`: Canberra distance, sensitive to small values.
;;         $$ D(P, Q) = \sum_i \frac{|P_i - Q_i|}{P_i + Q_i} $$
;;     *   `:lorentzian`: Lorentzian distance.
;;         $$ D(P, Q) = \sum_i \ln(1 + |P_i - Q_i|) $$
;;     *   `:non-intersection`: Non-intersection measure.
;;         $$ D(P, Q) = \frac{1}{2} \sum_i |P_i - Q_i| $$
;;     *   `:wave-hedges`: Wave Hedges distance.
;;         $$ D(P, Q) = \sum_i \frac{|P_i - Q_i|}{\max(P_i, Q_i)} $$
;;     *   `:czekanowski`: Czekanowski dissimilarity (same as Sorensen).
;;     *   `:motyka`: Motyka dissimilarity.
;;         $$ D(P, Q) = 1 - \frac{\sum_i \min(P_i, Q_i)}{\sum_i (P_i + Q_i)} $$
;;     *   `:tanimoto`: Tanimoto dissimilarity (extended Jaccard or Dice).
;;         $$ D(P, Q) = \frac{\sum_i (\max(P_i, Q_i) - \min(P_i, Q_i))}{\sum_i \max(P_i, Q_i)} $$
;;     *   `:jaccard`: Jaccard dissimilarity (1 - Jaccard similarity).
;;         $$ D(P, Q) = \frac{\sum_i (P_i - Q_i)^2}{\sum_i P_i^2 + \sum_i Q_i^2 - \sum_i P_i Q_i} $$
;;     *   `:dice`: Dice dissimilarity (1 - Dice similarity).
;;         $$ D(P, Q) = \frac{\sum_i (P_i - Q_i)^2}{\sum_i P_i^2 + \sum_i Q_i^2} $$
;;     *   `:bhattacharyya`: Bhattacharyya distance.
;;         $$ D(P, Q) = -\ln \left( \sum_i \sqrt{P_i Q_i} \right) $$
;;     *   `:hellinger`: Hellinger distance, derived from Bhattacharyya coefficient.
;;         $$ D(P, Q) = \sqrt{2 \sum_i (\sqrt{P_i} - \sqrt{Q_i})^2} $$
;;     *   `:matusita`: Matusita distance.
;;         $$ D(P, Q) = \sqrt{\sum_i (\sqrt{P_i} - \sqrt{Q_i})^2} $$

;; *   **Chi-squared based:**
;;     *   `:pearson-chisq` / `:chisq`: Pearson's Chi-squared statistic.
;;         $$ D(P, Q) = \sum_i \frac{(P_i - Q_i)^2}{Q_i} $$
;;     *   `:neyman-chisq`: Neyman's Chi-squared statistic.
;;         $$ D(P, Q) = \sum_i \frac{(P_i - Q_i)^2}{P_i} $$
;;     *   `:squared-chisq`: Squared Chi-squared distance.
;;         $$ D(P, Q) = \sum_i \frac{(P_i - Q_i)^2}{P_i + Q_i} $$
;;     *   `:symmetric-chisq`: Symmetric Chi-squared distance.
;;         $$ D(P, Q) = 2 \sum_i \frac{(P_i - Q_i)^2}{P_i + Q_i} $$
;;     *   `:divergence`: Divergence statistic.
;;         $$ D(P, Q) = 2 \sum_i \frac{(P_i - Q_i)^2}{(P_i + Q_i)^2} $$
;;     *   `:clark`: Clark distance.
;;         $$ D(P, Q) = \sqrt{\sum_i \left(\frac{P_i - Q_i}{P_i + Q_i}\right)^2} $$
;;     *   `:additive-symmetric-chisq`: Additive Symmetric Chi-squared distance.
;;         $$ D(P, Q) = \sum_i \frac{(P_i - Q_i)^2 (P_i + Q_i)}{P_i Q_i} $$

;; *   **Information Theory based (Divergences):** Measure the difference in information content.
;;     *   `:kullback-leibler`: Kullback-Leibler divergence (not symmetric, $KL(P||Q)$).
;;         $$ D(P, Q) = \sum_i P_i \ln\left(\frac{P_i}{Q_i}\right) $$
;;     *   `:jeffreys`: Jeffreys divergence (symmetric KL).
;;         $$ D(P, Q) = \sum_i (P_i - Q_i) \ln\left(\frac{P_i}{Q_i}\right) $$
;;     *   `:k-divergence`: K divergence (related to KL).
;;         $$ D(P, Q) = \sum_i P_i \ln\left(\frac{2 P_i}{P_i + Q_i}\right) $$
;;     *   `:topsoe`: Topsoe divergence.
;;         $$ D(P, Q) = \sum_i \left( P_i \ln\left(\frac{2 P_i}{P_i + Q_i}\right) + Q_i \ln\left(\frac{2 Q_i}{P_i + Q_i}\right) \right) $$
;;     *   `:jensen-shannon`: Jensen-Shannon divergence (symmetric, finite, based on KL).
;;         $$ D(P, Q) = \frac{1}{2} \left( KL(P || M) + KL(Q || M) \right), \text{ where } M = \frac{P+Q}{2} $$
;;     *   `:jensen-difference`: Jensen difference divergence.
;;         $$ D(P, Q) = \sum_i \left( \frac{P_i \ln P_i + Q_i \ln Q_i}{2} - \frac{(P_i+Q_i)}{2} \ln\left(\frac{P_i+Q_i}{2}\right) \right) $$
;;     *   `:taneja`: Taneja divergence.
;;         $$ D(P, Q) = \sum_i \frac{P_i + Q_i}{2} \ln\left(\frac{(P_i+Q_i)/2}{\sqrt{P_i Q_i}}\right) $$
;;     *   `:kumar-johnson`: Kumar-Johnson divergence.
;;         $$ D(P, Q) = \sum_i \frac{(P_i^2 - Q_i^2)^2}{2 (P_i Q_i)^{3/2}} $$

;; *   **Other:**
;;     *   `:avg`: Average of Manhattan and Chebyshev distances.
;;         $$ D(P, Q) = \frac{1}{2} \left( \sum_i |P_i - Q_i| + \max_i |P_i - Q_i| \right) $$

;; Let's use the sepal length samples from the iris dataset.

(utls/examples-note
  (stats/dissimilarity :euclidean setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :manhattan setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :chebyshev setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :minkowski setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :minkowski setosa-sepal-length virginica-sepal-length {:power 0.5})
  (stats/dissimilarity :euclidean-sq setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :squared-chord setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :sorensen setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :gower setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :kulczynski setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :canberra setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :lorentzian setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :non-intersection setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :wave-hedges setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :czekanowski setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :motyka setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :tanimoto setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :jaccard setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :dice setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :bhattacharyya setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :hellinger setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :matusita setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :pearson-chisq setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :neyman-chisq setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :squared-chisq setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :symmetric-chisq setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :divergence setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :clark setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :additive-symmetric-chisq setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :kullback-leibler setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :jeffreys setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :k-divergence setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :topsoe setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :jensen-shannon setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :jensen-difference setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :taneja setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :kumar-johnson setosa-sepal-length virginica-sepal-length)
  (stats/dissimilarity :avg setosa-sepal-length virginica-sepal-length))

;; We can compare our data to a distribution. The method used here is based on building histogram for P and quantize distribution for Q.

(utls/examples-note
  (stats/dissimilarity :chisq (stats/standardize setosa-sepal-length) r/default-normal)
  (stats/dissimilarity :chisq setosa-sepal-length (r/distribution :normal {:mu (stats/mean setosa-sepal-length) :sd (stats/stddev setosa-sepal-length)}))
  (stats/dissimilarity :chisq (repeatedly 1000 r/grand) r/default-normal))

;; In case when counts of samples are not equal we can use histograms. Also we can bin our data before comparison.

(utls/examples-note
  (stats/dissimilarity :gower (repeatedly 1000 r/grand) (repeatedly 800 r/grand) {:bins :auto})
  (stats/dissimilarity :gower setosa-sepal-length virginica-sepal-length {:bins :auto})
  (stats/dissimilarity :gower setosa-sepal-length virginica-sepal-length {:bins 10}))

;; #### Similarity Methods

;; Higher values generally indicate greater similarity.

;; *   **Overlap/Set-based:**
;;     *   `:intersection`: Intersection measure (sum of element-wise minimums).
;;         $$ S(P, Q) = \sum_i \min(P_i, Q_i) $$
;;     *   `:czekanowski`: Czekanowski similarity (same as Sorensen-Dice).
;;         $$ S(P, Q) = \frac{2 \sum_i \min(P_i, Q_i)}{\sum_i (P_i + Q_i)} $$
;;     *   `:motyka`: Motyka similarity.
;;         $$ S(P, Q) = \frac{\sum_i \min(P_i, Q_i)}{\sum_i (P_i + Q_i)} $$
;;     *   `:kulczynski`: Kulczynski similarity (can be > 1).
;;         $$ S(P, Q) = \frac{\sum_i \min(P_i, Q_i)}{\sum_i |P_i - Q_i|} $$
;;     *   `:ruzicka`: Ruzicka similarity.
;;         $$ S(P, Q) = \frac{\sum_i \min(P_i, Q_i)}{\sum_i \max(P_i, Q_i)} $$
;;     *   `:fidelity`: Probability fidelity (Bhattacharyya coefficient).
;;         $$ S(P, Q) = \sum_i \sqrt{P_i Q_i} $$
;;     *   `:squared-chord`: Squared chord similarity (1 - Squared Chord dissimilarity).
;;         $$ S(P, Q) = 2 \sum_i \sqrt{P_i Q_i} - 1 $$

;; *   **Inner Product / Angle:**
;;     *   `:inner-product`: Inner product of the vectors.
;;         $$ S(P, Q) = \sum_i P_i Q_i $$
;;     *   `:cosine`: Cosine similarity.
;;         $$ S(P, Q) = \frac{\sum_i P_i Q_i}{\sqrt{\sum_i P_i^2} \sqrt{\sum_i Q_i^2}} $$

;; *   **Set-based (adapted):**
;;     *   `:jaccard`: Jaccard similarity (generalized to distributions).
;;         $$ S(P, Q) = \frac{\sum_i P_i Q_i}{\sum_i P_i^2 + \sum_i Q_i^2 - \sum_i P_i Q_i} $$
;;     *   `:dice`: Dice similarity (generalized to distributions).
;;         $$ S(P, Q) = \frac{2 \sum_i P_i Q_i}{\sum_i P_i^2 + \sum_i Q_i^2} $$

;; *   **Harmonic Mean:**
;;     *   `:harmonic-mean`: Harmonic mean similarity.
;;         $$ S(P, Q) = 2 \sum_i \frac{P_i Q_i}{P_i + Q_i} $$

(utls/examples-note
  (stats/similarity :intersection setosa-sepal-length virginica-sepal-length)
  (stats/similarity :czekanowski setosa-sepal-length virginica-sepal-length)
  (stats/similarity :motyka setosa-sepal-length virginica-sepal-length)
  (stats/similarity :kulczynski setosa-sepal-length virginica-sepal-length)
  (stats/similarity :ruzicka setosa-sepal-length virginica-sepal-length)
  (stats/similarity :fidelity setosa-sepal-length virginica-sepal-length)
  (stats/similarity :squared-chord setosa-sepal-length virginica-sepal-length)
  (stats/similarity :inner-product setosa-sepal-length virginica-sepal-length)
  (stats/similarity :cosine setosa-sepal-length virginica-sepal-length)
  (stats/similarity :jaccard setosa-sepal-length virginica-sepal-length)
  (stats/similarity :dice setosa-sepal-length virginica-sepal-length)
  (stats/similarity :harmonic-mean setosa-sepal-length virginica-sepal-length))

;; As for `dissimilarity`, we can compare our data to a distribution. The method used here is based on building histogram for P and quantize distribution for Q.

(utls/examples-note
  (stats/similarity :dice (stats/standardize setosa-sepal-length) r/default-normal)
  (stats/similarity :dice setosa-sepal-length (r/distribution :normal {:mu (stats/mean setosa-sepal-length) :sd (stats/stddev setosa-sepal-length)}))
  (stats/similarity :dice (repeatedly 10000 r/grand) r/default-normal))

;; In case when counts of samples are not equal we can use histograms. Also we can bin our data before comparison.

(utls/examples-note
  (stats/similarity :ruzicka (repeatedly 1000 r/grand) (repeatedly 800 r/grand) {:bins :auto})
  (stats/similarity :ruzicka setosa-sepal-length virginica-sepal-length {:bins :auto})
  (stats/similarity :ruzicka setosa-sepal-length virginica-sepal-length {:bins 10}))

;; ## Contingency Tables

;; Functions for creating and analyzing contingency tables.

;; ::: {.callout-tip title="Defined functions"}
;; * `contingency-table`, `rows->contingency-table`, `contingency-table->marginals`
;; * `contingency-2x2-measures-all`, `contingency-2x2-measures`
;; * `mcc`
;; * `cramers-c`, `cramers-v`, `cramers-v-corrected`
;; * `cohens-w`, `tschuprows-t`
;; * `cohens-kappa`, `weighted-kappa`
;; :::

;; Contingency tables, also known as cross-tabulations or cross-tabs, are a fundamental tool in statistics for displaying the frequency distribution of two or more categorical variables. They help visualize and analyze the relationship or association between these variables.
;; 
;; For two variables, a contingency table typically looks like this:
;; 
;; |          | Category 1 (Col) | Category 2 (Col) | ... | Column Totals (Marginals) |
;; |:---------|:-----------------|:-----------------|:----|:--------------------------|
;; | **Cat A (Row)** | $n_{11}$         | $n_{12}$         | ... | $R_1 = \sum_j n_{1j}$     |
;; | **Cat B (Row)** | $n_{21}$         | $n_{22}$         | ... | $R_2 = \sum_j n_{2j}$     |
;; | **...**  | ...              | ...              | ... | ...                       |
;; | **Row Totals (Marginals)** | $C_1 = \sum_i n_{i1}$ | $C_2 = \sum_i n_{i2}$ | ... | $N = \sum_i \sum_j n_{ij}$ |
;; 
;; Where $n_{ij}$ is the count of observations falling into row category $i$ and column category $j$. $R_i$ are row marginal totals, $C_j$ are column marginal totals, and $N$ is the grand total number of observations.

;; Let's use data from [wikipedia article](https://en.wikipedia.org/wiki/Contingency_table#Example) as the 2x2 example.

(kind/table
 [[" " "Right-handed" "Left-handed" "Total"]
  ["Male" 43 9 52]
  ["Female" 44 4 48]
  ["Total" 87 13 100]])

(def ct-data-1 [[43 9] [44 4]])

;; Another example will be from [openstax book](https://openstax.org/books/introductory-business-statistics-2e/pages/3-4-contingency-tables-and-probability-trees#M05_ch03-tbl005)

(kind/table
 [[" " "Lake Path" "Hilly Path" "Wooded Path" "Total"]
  ["Younger" 45 38 27 110]
  ["Older" 26 52 12 90]
  ["Total" 71 90 39 200]])

(def ct-data-2 [[45 38 27] [26 52 12]])

;; The last example will consist two sequences for gender and exam outcome

(def gender [:male :female :male :female :male :male :female :female :female :male :male])
(def outcome [:pass :fail :pass :pass :fail :pass :fail :pass :pass :pass :pass])


;; `fastmath.stats` provides functions for creating and analyzing these tables:
;; 
;; *   `contingency-table`: Creates a frequency map from one or more sequences. If given two sequences of equal length, say `vars1` and `vars2`, it produces a map where keys are pairs `[value-from-vars1, value-from-vars2]` and values are the counts of these co-occurrences.
;; *   `rows->contingency-table`: Takes a sequence of sequences, interpreted as rows of counts in a grid, and converts it into a map format where keys are `[row-index, column-index]` and values are the non-zero counts. This is useful for inputting tables structured as lists of lists.
;; *   `contingency-table->marginals`: Calculates the row totals (`:rows`), column totals (`:cols`), grand total (`:n`), and diagonal elements (`:diag`) from a contingency table (either in map format or sequence of sequences format).
;; 
;; Let's create a simple contingency table from above data.


(def ct-gender-outcome (stats/contingency-table gender outcome))
(def contingency-table-1 (stats/rows->contingency-table ct-data-1))
(def contingency-table-2 (stats/rows->contingency-table ct-data-2))

(utls/examples-note
  ct-gender-outcome
  contingency-table-1
  contingency-table-2
  (stats/contingency-table->marginals ct-gender-outcome)
  (stats/contingency-table->marginals contingency-table-1)
  (stats/contingency-table->marginals contingency-table-2))

;; **Measures of Association:** These statistics quantify the strength and nature of the relationship between the variables in the table. They are often derived from the Pearson's Chi-squared statistic ($\chi^2$, obtainable via [[chisq-test]]), which tests for independence.
;; 
;; *   `cramers-c`: Cramer's C is a measure of association for any $R \times K$ contingency table. It ranges from 0 to 1, where 0 indicates no association and 1 indicates perfect association.
;;     $$ C = \sqrt{\frac{\chi^2}{N + \chi^2}} $$
;; *   `cramers-v`: Cramer's V is another measure of association, also ranging from 0 to 1. It is widely used and is often preferred over Tschuprow's T because it can attain the value 1 even for non-square tables.
;;     $$ V = \sqrt{\frac{\chi^2/N}{\min(R-1, K-1)}} $$
;; *   `cramers-v-corrected`: Corrected Cramer's V ($V^*$) applies a bias correction, which is particularly important for small sample sizes or tables with many cells having low expected counts.
;; *   `cohens-w`: Cohen's W is a measure of effect size for Chi-squared tests. It quantifies the magnitude of the difference between observed and expected frequencies. It ranges from 0 upwards, with 0 indicating no difference (independence).
;;     $$ W = \sqrt{\frac{\chi^2}{N}} $$
;; *   `tschuprows-t`: Tschuprow's T is a measure of association ranging from 0 to 1. However, it can only reach 1 in square tables ($R=K$).
;;     $$ T = \sqrt{\frac{\chi^2/N}{\sqrt{(R-1)(K-1)}}} $$
;; 
;; `fastmath.stats` allows you to calculate these measures by providing either the raw sequences or a pre-calculated contingency table:
;; 

(utls/examples-note
  (stats/cramers-c ct-gender-outcome)
  (stats/cramers-v ct-gender-outcome)
  (stats/cramers-v-corrected ct-gender-outcome)
  (stats/cohens-w ct-gender-outcome)
  (stats/tschuprows-t ct-gender-outcome)
  
  (stats/cramers-c contingency-table-1)
  (stats/cramers-v contingency-table-1)
  (stats/cramers-v-corrected contingency-table-1)
  (stats/cohens-w contingency-table-1)
  (stats/tschuprows-t contingency-table-1)

  (stats/cramers-c contingency-table-2)
  (stats/cramers-v contingency-table-2)
  (stats/cramers-v-corrected contingency-table-2)
  (stats/cohens-w contingency-table-2)
  (stats/tschuprows-t contingency-table-2))

;; **Measures of Agreement:** These statistics are typically used for square tables ($R=K$) to assess the consistency of agreement between two independent raters or methods that categorize items.
;; 
;; *   `cohens-kappa`: Cohen's Kappa ($\kappa$) measures the agreement between two raters for nominal (or ordinal) categories, correcting for the agreement that would be expected by chance. It ranges from -1 (perfect disagreement) to +1 (perfect agreement), with 0 indicating agreement equivalent to chance.
;;     $$ \kappa = \frac{p_0 - p_e}{1 - p_e} $$
;;     where $p_0$ is the observed proportional agreement (sum of diagonal cells divided by N) and $p_e$ is the proportional agreement expected by chance (based on marginal probabilities).
;; *   `weighted-kappa`: Weighted Kappa ($\kappa_w$) is an extension of Cohen's Kappa that allows for different levels of disagreement penalties, suitable for ordinal categories. Disagreements between categories that are closer together (e.g., "Good" vs "Very Good") are penalized less than disagreements between categories that are further apart (e.g., "Poor" vs "Excellent"). It requires specifying a weighting scheme (e.g., `:equal-spacing`, `:fleiss-cohen`).
;; 
;; Examples for agreement measures (using the [[rows->contingency-table]] format for clarity on cell positions):

(utls/examples-note
  (stats/cohens-kappa ct-gender-outcome)
  (stats/cohens-kappa contingency-table-1)

  (stats/weighted-kappa contingency-table-1 :equal-spacing)
  (stats/weighted-kappa contingency-table-1 :fleiss-cohen)
  (stats/weighted-kappa contingency-table-1 {[0 0] 0.1 [1 0] 0.2 [0 1] 0.4 [1 1] 0.3})
  (stats/weighted-kappa contingency-table-1 (fn [dim row-id col-id] (/ (max row-id col-id) dim))))

;; **2x2 Specific Measures:** These functions are tailored for 2x2 tables.

;; *   `mcc`: Calculates the [[mcc]], which is equivalent to the Phi coefficient for a 2x2 table. It is a single, balanced measure of classification performance, ranging from -1 to +1.
;; *   `contingency-2x2-measures`: A convenience function that returns a map containing a selection of commonly used statistics specifically for 2x2 tables (Chi-squared, Kappa, Phi, Yule's Q, Odds Ratio, etc.).
;; *   `contingency-2x2-measures-all`: Provides a very comprehensive map of measures for a 2x2 table, including various Chi-squared statistics and p-values, measures of association, agreement, and risk/effect size measures like Odds Ratio (OR), Relative Risk (RR), and Number Needed to Treat (NNT). It's the most detailed summary for a 2x2 table.

(utls/examples-note
  (stats/mcc contingency-table-1))

;; The `contingency-2x2-measures-all` function accepts the counts `a, b, c, d` directly, a sequence `[a b c d]`, a matrix `[[a b] [c d]]`, or a map `{:a a :b b :c c :d d}`.

;; The map returned by `contingency-2x2-measures-all` is extensive. Here's a description of its main keys and their contents:
;; 
;; *   `:n`: The grand total number of observations in the table ($N$).
;; *   `:table`: A map representation of the input counts, typically `{:a a, :b b, :c c, :d d}` corresponding to the top-left, top-right, bottom-left, bottom-right cells.
;; *   `:expected`: A map of the expected counts `{:a exp_a, :b exp_b, :c exp_c, :d exp_d}` for each cell under the assumption of independence.
;; *   `:marginals`: A map containing the row totals (`:row1`, `:row2`), column totals (`:col1`, `:col2`), and the grand total (`:total`).
;; *   `:proportions`: A nested map providing proportions:
;;     *   `:table`: Cell counts divided by the grand total (`N`).
;;     *   `:rows`: Cell counts divided by their respective row totals (conditional probabilities of columns given rows).
;;     *   `:cols`: Cell counts divided by their respective column totals (conditional probabilities of rows given columns).
;;     *   `:marginals`: Row and column totals divided by the grand total (`N`).
;; *   `:p-values`: A map containing p-values for various Chi-squared tests:
;;     *   `:chi2`: Pearson's Chi-squared test p-value.
;;     *   `:yates`: Yates' continuity corrected Chi-squared test p-value.
;;     *   `:cochran-mantel-haenszel`: Cochran-Mantel-Haenszel statistic p-value (useful for stratified data, but calculated for the single table here).
;; *   `:OR`: The Odds Ratio ($OR$) for the 2x2 table, typically $(a \times d) / (b \times c)$. Quantifies the strength of association, representing the odds of an outcome occurring in one group compared to the odds in another group.
;; *   `:lOR`: The natural logarithm of the Odds Ratio ($\ln(OR)$). Useful for constructing confidence intervals for the OR.
;; *   `:RR`: The Relative Risk ($RR$) for the 2x2 table, typically $(a/(a+b)) / (c/(c+d))$. Compares the probability of an outcome in one group to the probability in another group. Requires row marginals to represent the exposed/unexposed or intervention/control groups.
;; *   `:risk`: A map containing various risk-related measures, especially relevant in epidemiology or clinical trials, derived from the counts assuming a specific row structure (e.g., exposure/outcome):
;;     *   `:RR`: Relative Risk (same as top-level `:RR`).
;;     *   `:RRR`: Relative Risk Reduction ($1 - RR$).
;;     *   `:RD`: Risk Difference ($a/(a+b) - c/(c+d)$).
;;     *   `:ES`: Exposure Sample size (row 1 total).
;;     *   `:CS`: Control Sample size (row 2 total).
;;     *   `:EER`: Experimental Event Rate ($a/(a+b)$).
;;     *   `:CER`: Control Event Rate ($c/(c+d)$).
;;     *   `:ARR`: Absolute Risk Reduction ($CER - EER$).
;;     *   `:NNT`: Number Needed to Treat ($1 / ARR$). The average number of patients who need to be treated to prevent one additional bad outcome.
;;     *   `:ARI`: Absolute Risk Increase ($-ARR$).
;;     *   `:NNH`: Number Needed to Harm ($-NNT$).
;;     *   `:RRI`: Relative Risk Increase ($RR - 1$).
;;     *   `:AFe`: Attributable Fraction among the exposed ($(EER - CER)/EER$).
;;     *   `:PFu`: Prevented Fraction among the unexposed ($1 - RR$).
;; *   `:SE`: Standard Error related values used in some calculations.
;; *   `:measures`: A map containing a wide array of statistical measures:
;;     *   `:chi2`: Pearson's Chi-squared statistic (same as `:p-values :chi2` but the statistic value).
;;     *   `:yates`: Yates' continuity corrected Chi-squared statistic.
;;     *   `:cochran-mantel-haenszel`: Cochran-Mantel-Haenszel statistic.
;;     *   `:cohens-kappa`: Cohen's Kappa coefficient.
;;     *   `:yules-q`: Yule's Q, a measure of association based on the Odds Ratio. Ranges -1 to +1.
;;     *   `:holley-guilfords-g`: Holley-Guilford's G.
;;     *   `:huberts-gamma`: Hubert's Gamma.
;;     *   `:yules-y`: Yule's Y, based on square roots of cell proportions. Ranges -1 to +1.
;;     *   `:cramers-v`: Cramer's V measure of association.
;;     *   `:phi`: Phi coefficient (equivalent to MCC for 2x2 tables), calculated as $\frac{ad-bc}{\sqrt{(a+b)(c+d)(a+c)(b+d)}}$. Ranges -1 to +1.
;;     *   `:scotts-pi`: Scott's Pi, another measure of inter-rater reliability, similar to Kappa.
;;     *   `:cohens-h`: Cohen's H, a measure of effect size for comparing two proportions.
;;     *   `:PCC`: Pearson's Contingency Coefficient.
;;     *   `:PCC-adjusted`: Adjusted Pearson's Contingency Coefficient.
;;     *   `:TCC`: Tschuprow's Contingency Coefficient.
;;     *   `:F1`: F1 Score, common in binary classification, the harmonic mean of precision and recall. Calculated assuming `a=TP, b=FP, c=FN, d=TN`.
;;     *   `:bangdiwalas-b`: Bangdiwala's B, a measure of agreement/reproducibility for clustered data, but calculated for the 2x2 table.
;;     *   `:mcnemars-chi2`: McNemar's Chi-squared test statistic, specifically for paired nominal data (e.g., before/after) to test for changes in proportions. Calculated from off-diagonal elements ($b$ and $c$).
;;     *   `:gwets-ac1`: Gwet's AC1, another measure of inter-rater reliability claimed to be more robust to prevalence issues than Kappa.

(utls/zp
 (stats/contingency-2x2-measures-all (flatten ct-data-1)))

;; ### Binary Classification Metrics

;; Metrics derived from a 2x2 confusion matrix are essential for evaluating the performance of algorithms in binary classification tasks (problems where the outcome belongs to one of two classes, e.g., "positive" or "negative", "success" or "failure"). These metrics quantify how well a classifier distinguishes between the two classes.

;; ::: {.callout-tip title="Defined functions"}
;; * `confusion-matrix`
;; * `binary-measures-all`, `binary-measures`
;; :::

;; The foundation for these metrics is the **confusion matrix**, a 2x2 table summarizing the results by comparing the actual class of each instance to the predicted class:
;; 
;; |                | Predicted Positive (P') | Predicted Negative (N') | Total
;; |:---------------|:------------------------|:------------------------|:----
;; | **Actual Positive (P)**  | True Positives (TP)     | False Negatives (FN)    | P
;; | **Actual Negative (N)** | False Positives (FP)    | True Negatives (TN)     | N
;; | **Total**      | P'                      | N'                      | Total
;; 
;; *   **True Positives (TP):** Instances correctly predicted as positive.
;; *   **False Negatives (FN):** Instances incorrectly predicted as negative (Type II error). These are positive instances missed by the classifier.
;; *   **False Positives (FP):** Instances incorrectly predicted as positive (Type I error). These are negative instances misclassified as positive.
;; *   **True Negatives (TN):** Instances correctly predicted as negative.
;; 
;; The counts from these four cells form the basis for almost all binary classification metrics.
;; 
;; `fastmath.stats` provides functions to generate and analyze these matrices:
;; 
;; *   `confusion-matrix`: This function constructs the 2x2 confusion matrix counts (`:tp`, `:fn`, `:fp`, `:tn`). It can take the four counts directly, a structured map/sequence representation, or two sequences of actual and predicted outcomes. It can also handle different data types for outcomes using an `encode-true` parameter.

(utls/examples-note
  (stats/confusion-matrix 10 2 5 80)
  (stats/confusion-matrix [10 2 5 80])
  (stats/confusion-matrix [[10 5] [2 80]])
  (stats/confusion-matrix {:tp 10 :fn 2 :fp 5 :tn 80})
  (stats/confusion-matrix [:pos :neg :pos :pos :neg :pos :neg :pos :neg :pos :pos]
                          [:pos :neg :pos :pos :neg :pos :neg :neg :pos :pos :pos]
                          #{:pos})
  (stats/confusion-matrix [1 0 1 1 0 1 0 1 0 1 1]
                          [1 0 1 1 0 1 0 0 1 1 1])
  (stats/confusion-matrix [true false true true false true false true false true true]
                          [true false true true false true false false true true true]))

;; *   `binary-measures-all`: This is the primary function for calculating a wide range of metrics from a 2x2 confusion matrix. It accepts the same input formats as `confusion-matrix`. It returns a map containing a comprehensive set of derived statistics.
;; 
;; *   `binary-measures`: A convenience function that calls `binary-measures-all` but returns a smaller, more commonly used subset of the metrics. It accepts the same input formats.
;;

;; The map returned by `binary-measures-all` contains numerous metrics. Key metrics include:
;; 
;; *   `:tp`, `:fn`, `:fp`, `:tn`: The raw counts from the confusion matrix.
;; *   `:total`: Total number of instances ($TP + FN + FP + TN$).
;; *   `:cp`, `:cn`, `:pcp`, `:pcn`: Marginal totals (e.g., `:cp` is actual positives $TP + FN$).
;; *   `:accuracy`: Overall proportion of correct predictions:
;;     $$ Accuracy = \frac{TP + TN}{TP + FN + FP + TN} $$
;; *   `:sensitivity`, `:recall`, `:tpr`: True Positive Rate (proportion of actual positives correctly identified):
;;     $$ Sensitivity = \frac{TP}{TP + FN} $$
;; *   `:specificity`, `:tnr`: True Negative Rate (proportion of actual negatives correctly identified):
;;     $$ Specificity = \frac{TN}{FP + TN} $$
;; *   `:precision`, `:ppv`: Positive Predictive Value (proportion of positive predictions that were actually positive):
;;     $$ Precision = \frac{TP}{TP + FP} $$
;; *   `:fdr`: False Discovery Rate ($1 - Precision$).
;; *   `:npv`: Negative Predictive Value (proportion of negative predictions that were actually negative):
;;     $$ NPV = \frac{TN}{FN + TN} $$
;; *   `:for`: False Omission Rate ($1 - NPV$).
;; *   `:fpr`: False Positive Rate (proportion of actual negatives incorrectly identified as positive):
;;     $$ FPR = \frac{FP}{FP + TN} $$
;; *   `:fnr`: False Negative Rate (proportion of actual positives incorrectly identified as negative):
;;     $$ FNR = \frac{FN}{TP + FN} $$
;; *   `:f1-score`, `:f-measure`: The harmonic mean of Precision and Recall. It provides a balance between the two:
;;     $$ F1 = 2 \cdot \frac{Precision \cdot Recall}{Precision + Recall} $$
;; *   `:mcc`, `:phi`: Matthews Correlation Coefficient / Phi coefficient. A balanced measure ranging from -1 to +1, considered robust for imbalanced datasets:
;;     $$ MCC = \frac{TP \cdot TN - FP \cdot FN}{\sqrt{(TP+FP)(TP+FN)(TN+FP)(TN+FN)}} $$
;; *   `:prevalence`: Proportion of the positive class in the actual data:
;;     $$ Prevalence = \frac{TP + FN}{TP + FN + FP + TN} $$
;; *   Other metrics like `:ba` (Balanced Accuracy), `:lr+`, `:lr-`, `:dor`, `:fm`, `:pt`, `:bm`, `:kappa` (Cohen's Kappa for 2x2), `:mk`, and `:f-beta` (a function to calculate F-beta score for any beta value).
;;
;; `binary-measures` returns a map containing `:tp`, `:tn`, `:fp`, `:fn`, `:accuracy`, `:fdr`, `:f-measure`, `:fall-out` (FPR), `:precision`, `:recall` (Sensitivity/TPR), `:sensitivity`, `:specificity` (TNR), and `:prevalence`. This provides a useful quick summary. For a deeper dive or specific niche metrics, use `binary-measures-all`.

;; Let's calculate the metrics for a sample confusion matrix `{:tp 10 :fn 2 :fp 5 :tn 80}`.

(utls/zp
 (stats/binary-measures-all {:tp 10 :fn 2 :fp 5 :tn 80}))

(utls/zp
 (stats/binary-measures {:tp 10 :fn 2 :fp 5 :tn 80}))

;; ## Effect Size

;; Measures quantifying the magnitude of a phenomenon or relationship.

;; ### Difference Family

;; Measures quantifying the magnitude of the difference between two group means, standardized by a measure of variability. These are widely used effect size statistics, particularly in comparing experimental and control groups or two distinct conditions.

;; ::: {.callout-tip title="Defined functions"}
;; * `cohens-d`, `cohens-d-corrected`
;; * `hedges-g`, `hedges-g-corrected`, `hedges-g*`
;; * `glass-delta`
;; :::

;; These functions primarily quantify the difference between the means of two groups (usually group 1 and group 2: $\bar{x}_1 - \bar{x}_2$), scaled by some estimate of the population standard deviation.
;; 
;; * **Cohen's d (`cohens-d`)**: Measures the difference between two means divided by the pooled standard deviation. It assumes equal variances between the groups. The formula is:
;;     $$ d = \frac{\bar{x}_1 - \bar{x}_2}{s_p} $$
;;     where $s_p = \sqrt{\frac{(n_1-1)s_1^2 + (n_2-1)s_2^2}{n_1+n_2-2}}$ is the *unbiased* pooled standard deviation. An optional `method` argument can be passed to select the pooled standard deviation calculation method (see [[pooled-stddev]]).
;; 
;; * **Hedges' g (`hedges-g`)**: In `fastmath.stats`, this function is equivalent to `cohens-d` using the default (unbiased) pooled standard deviation. Conceptually, Hedges' g also divides the mean difference by a pooled standard deviation, but often refers to a bias-corrected version for small samples (see below).
;;     $$ g = \frac{\bar{x}_1 - \bar{x}_2}{s_p} $$
;; 
;; * **Bias-Corrected Cohen's d (`cohens-d-corrected`) / Bias-Corrected Hedges' g (`hedges-g-corrected`)**: These functions apply a small-sample bias correction factor to Cohen's d or Hedges' g. This correction factor (an approximation) is multiplied by the calculated d or g value to provide a less biased estimate of the population effect size, especially for small sample sizes ($n < 20$). The correction factor is approximately $1 - \frac{3}{4\nu - 1}$, where $\nu = n_1+n_2-2$ are the degrees of freedom.
;; 
;; * **Exact Bias-Corrected Hedges' g (`hedges-g*`)**: Applies the precise bias correction factor $J(\nu)$ (based on the Gamma function) to Hedges' g (or Cohen's d with unbiased pooling). This provides the most theoretically accurate bias correction for small samples.
;;     $$ g^* = g \cdot J(\nu) $$
;;     where $\nu = n_1+n_2-2$ and $J(\nu) = \frac{\Gamma(\nu/2)}{\sqrt{\nu/2} \Gamma((\nu-1)/2)}$.
;; 
;; * **Glass's Delta (`glass-delta`)**: Measures the difference between two means divided by the standard deviation of the *control* group (conventionally group 2).
;;     $$ \Delta = \frac{\bar{x}_1 - \bar{x}_2}{s_2} $$
;;     This is useful when the control group's variance is considered a better estimate of the population variance than a pooled variance, or when the intervention is expected to affect variance.
;;

;; The **returned value** from these functions represents the difference between the means of Group 1 and Group 2, expressed in units of the relevant standard deviation (pooled for Cohen's/Hedges', control group's for Glass's).
;; 
;; *   A positive value indicates that the mean of Group 1 is greater than the mean of Group 2.
;; *   A negative value indicates that the mean of Group 1 is less than the mean of Group 2.
;; *   The magnitude of the value indicates the size of the difference relative to the spread of the data. Cohen's informal guidelines for interpreting the magnitude of *d* (and often *g*) are:
;;     *   $|d| \approx 0.2$: small effect
;;     *   $|d| \approx 0.5$: medium effect
;;     *   $|d| \approx 0.8$: large effect

;; **Comparison:**
;; 
;; *   `cohens-d` and `hedges-g` (as implemented here) are standard measures for the mean difference scaled by pooled standard deviation, assuming equal variances. Cohen's d is more commonly cited, while Hedges' g is often used when sample size is small.
;; *   `cohens-d-corrected`, `hedges-g-corrected`, and `hedges-g*` are preferred over the standard `d`/`g` when sample sizes are small, as they provide less biased estimates of the population effect size. `hedges-g*` uses the most accurate correction formula.
;; *   `glass-delta` is distinct in using only the control group's standard deviation for scaling. Use this when the variance of the control group is more appropriate as a baseline (e.g., comparing an experimental group to a control, or if the intervention might affect variance).
;; 
;; Let's illustrate these with the sepal lengths from the iris dataset, using 'virginica' as group 1 and 'setosa' as group 2.

(utls/examples-note
  (stats/cohens-d virginica-sepal-length setosa-sepal-length)
  (stats/cohens-d-corrected virginica-sepal-length setosa-sepal-length)
  (stats/cohens-d virginica-sepal-length setosa-sepal-length :biased)
  (stats/cohens-d-corrected virginica-sepal-length setosa-sepal-length :biased)
  (stats/cohens-d virginica-sepal-length setosa-sepal-length :avg)
  (stats/cohens-d-corrected virginica-sepal-length setosa-sepal-length :avg)
  (stats/hedges-g virginica-sepal-length setosa-sepal-length)
  (stats/hedges-g-corrected virginica-sepal-length setosa-sepal-length)
  (stats/hedges-g* virginica-sepal-length setosa-sepal-length)
  (stats/glass-delta virginica-sepal-length setosa-sepal-length))


;; ### Ratio Family

;; Effect size measures in the ratio family quantify the difference between two group means as a multiplicative factor, rather than a difference. This is useful when comparing values that are inherently ratios or when expressing the effect in terms of relative change.

;; ::: {.callout-tip title="Defined functions"}
;; * `means-ratio`, `means-ratio-corrected`
;; :::

;; * **Ratio of Means** (`means-ratio`): Calculates the ratio of the mean of the first group ($\bar{x}_1$) to the mean of the second group ($\bar{x}_2$).
;;     $$ Ratio = \frac{\bar{x}_1}{\bar{x}_2} $$
;;     A value greater than 1 means the first group's mean is larger; less than 1 means it's smaller. A value of 1 indicates equal means.
;; 
;; * **Corrected Ratio of Means** (`means-ratio-corrected`): Applies a small-sample bias correction to the simple ratio of means. This provides a less biased estimate of the population ratio, especially for small sample sizes. The correction is based on incorporating the variances of the two groups. This function is equivalent to calling `(means-ratio group1 group2 true)`.
;; 
;; Use `means-ratio` for a direct, uncorrected ratio. Use `means-ratio-corrected` for a less biased estimate of the population ratio, particularly advisable with small sample sizes.
;; 
;; Let's calculate these for the virginica and setosa sepal lengths from the iris dataset.

(utls/examples-note
  (stats/means-ratio virginica-sepal-length setosa-sepal-length)
  (stats/means-ratio-corrected virginica-sepal-length setosa-sepal-length))

;; ### Ordinal / Non-parametric Family

;; Effect size measures in the ordinal or non-parametric family are suitable for data that are not normally distributed, have outliers, or are measured on an ordinal scale. They are often based on comparing pairs of observations between two groups.

;; ::: {.callout-tip title="Defined functions"}
;; * `cliffs-delta`
;; * `ameasure`, `wmw-odds`
;; :::

;; *   **Cliff's Delta** (`cliffs-delta`): A non-parametric measure of the amount of separation or overlap between two distributions. It is calculated as the probability that a randomly selected observation from one group is greater than a randomly selected observation from the other group, minus the probability of the reverse:
;;
;;     $$\delta = P(X > Y) - P(Y > X)$$
;;
;;     It ranges from -1 (complete separation, with values in the second group always greater than values in the first) to +1 (complete separation, with values in the first group always greater than values in the second). A value of 0 indicates complete overlap (stochastic equality). It is a robust alternative to Cohen's d when assumptions for parametric tests are not met.
;;
;; *   **Vargha-Delaney A** (`ameasure`): A non-parametric measure of stochastic superiority. It quantifies the probability that a randomly chosen value from the first sample (`group1`) is greater than a randomly chosen value from the second sample (`group2`).
;;
;;     $$A = P(X > Y)$$
;;
;;     It ranges from 0 to 1. A value of 0.5 indicates stochastic equality (distributions are overlapping). Values > 0.5 mean `group1` tends to have larger values; values < 0.5 mean `group2` tends to have larger values. It is directly related to Cliff's Delta: $A = (\delta + 1)/2$.
;;
;; *   **Wilcoxon-Mann-Whitney Odds** (`wmw-odds`): This non-parametric measure quantifies the odds that a randomly chosen observation from the first group (`group1`) is greater than a randomly chosen observation from the second group (`group2`).
;;
;;     $$\psi = \frac{P(X > Y)}{P(Y > X)} = \frac{1 + \delta}{1 - \delta}$$
;;
;;     It ranges from 0 to infinity. A value of 1 indicates stochastic equality. Values > 1 mean `group1` values tend to be larger; values < 1 mean `group1` values tend to be smaller. The natural logarithm of ψ is the log-odds that a random observation from group 1 is greater than a random observation from group 2.
;;
;; These three measures (`cliffs-delta`, `ameasure`, `wmw-odds`) are non-parametric and robust, making them suitable for ordinal data or when assumptions of normality and equal variances are violated. They are inter-related and can be transformed into one another. Cliff's Delta is centered around 0, A is centered around 0.5, and WMW Odds is centered around 1 (or 0 on a log scale), offering different interpretations of the effect magnitude and direction.
;;
;; Let's illustrate these with the sepal lengths from the iris dataset, using 'virginica' as group 1 and 'setosa' as group 2.

(utls/examples-note
  (stats/cliffs-delta virginica-sepal-length setosa-sepal-length)
  (stats/ameasure virginica-sepal-length setosa-sepal-length)
  (stats/wmw-odds virginica-sepal-length setosa-sepal-length))

;; ### Overlap Family

;; Measures quantifying the degree to which two distributions share common values. These metrics assess how much the data from one group tends to blend or overlap with the data from another group. They provide an alternative perspective on effect size compared to mean differences or ratios, which is particularly useful when the focus is on the proportion of scores that are shared or distinct between groups.

;; ::: {.callout-tip title="Defined functions"}
;; * `p-overlap`
;; * `cohens-u1-normal`, `cohens-u2-normal`, `cohens-u3-normal`
;; * `cohens-u1`, `cohens-u2`, `cohens-u3`
;; :::

;; These functions provide different ways to quantify overlap:
;; 
;; *   **`p-overlap`**: Estimates the proportion of overlap between two distributions based on Kernel Density Estimation (KDE). It calculates the area under the minimum of the two estimated density functions.
;; 
;;     $$ \text{Overlap} = \int_{-\infty}^{\infty} \min(f_1(x), f_2(x)) \, dx $$
;; 
;;     This measure is non-parametric, meaning it doesn't assume the data comes from a specific distribution (like normal), and it is symmetric (`p-overlap(A, B) == p-overlap(B, A)`). The result is a proportion between 0 (no overlap) and 1 (complete overlap, identical distributions).

(utls/examples-note
  (stats/p-overlap virginica-sepal-length setosa-sepal-length)
  (stats/p-overlap virginica-sepal-length setosa-sepal-length {:kde :epanechnikov :bandwidth 0.7}))

(kind/table
 [[(gg/->image (gg/functions [["virginica" (fastmath.kernel/kernel-density :gaussian virginica-sepal-length)]
                              ["setosa" (fastmath.kernel/kernel-density :gaussian setosa-sepal-length)]]
                             {:x [2 10]
                              :title "KDE with gaussian kernel"}))

   (gg/->image (gg/functions [["virginica" (fastmath.kernel/kernel-density :epanechnikov virginica-sepal-length 0.7)]
                              ["setosa" (fastmath.kernel/kernel-density :epanechnikov setosa-sepal-length 0.7)]]
                             {:x [2 10]
                              :title "KDE with epanechnikov(0.7) kernel"}))]])

;; *   **Cohen's U measures (`cohens-u1-normal`, `cohens-u2-normal`, `cohens-u3-normal`)**: These are measures of non-overlap or overlap derived from Cohen's d, assuming **normal distributions** with **equal variances**. They quantify the proportion of one group's scores that fall beyond a certain point in relation to the other group.
;;     *   **`cohens-u1-normal`**: Quantifies the proportion of scores in the *lower*-scoring group that overlap with the scores in the *higher*-scoring group. Calculated from Cohen's d ($d$) using the standard normal CDF ($\Phi$) as $(2\Phi(|d|/2) - 1) / \Phi(|d|/2)$. Range is [0, 1]. 0 means no overlap, 1 means complete overlap.
;;     *   **`cohens-u2-normal`**: Quantifies the proportion of scores in the *lower*-scoring group that are below the point located halfway between the means of the two groups. Calculated as $\Phi(|d|/2)$. Range is [0, 1]. 0.5 means perfect overlap (medians are at the same point), 0 or 1 means no overlap (distributions are far apart).
;;     *   **`cohens-u3-normal`**: Quantifies the proportion of scores in the *lower*-scoring group that fall below the mean of the *higher*-scoring group. Calculated as $\Phi(d)$. Range is [0, 1]. This measure is **asymmetric**: `U3(A, B)` is the proportion of B below A's mean; `U3(B, A)` is the proportion of A below B's mean.

(utls/examples-note
  (stats/cohens-u1-normal virginica-sepal-length setosa-sepal-length)
  (stats/cohens-u2-normal virginica-sepal-length setosa-sepal-length)
  (stats/cohens-u3-normal virginica-sepal-length setosa-sepal-length))

(gg/->image (gg/functions [["virginica" (partial r/pdf
                                                 (r/distribution
                                                  :normal
                                                  {:mu (stats/mean virginica-sepal-length)
                                                   :sd (stats/stddev virginica-sepal-length)}))]
                           ["setosa" (partial r/pdf
                                              (r/distribution
                                               :normal
                                               {:mu (stats/mean setosa-sepal-length)
                                                :sd (stats/stddev setosa-sepal-length)}))]]
                          {:x [2 10]
                           :title "Sepal Length distributions assuming normality."}))

;; *   **Non-parametric Cohen's U measures (`cohens-u1`, `cohens-u2`, `cohens-u3`)**: These are analogous measures that do **not assume normality** or equal variances. They are based on comparing percentiles or medians of the empirical distributions.
;;     *   **`cohens-u1`**: Non-parametric measure of difference/separation, related to [[cohens-u2]]. Its value indicates greater separation as it increases towards +1, and more overlap as it approaches -1. Symmetric.
;;     *   **`cohens-u2`**: Quantifies the minimum overlap between corresponding quantiles. It finds the smallest distance between quantiles $q$ and $1-q$ of the two distributions across $q \in [0.5, 1.0]$. Range is typically [0, max_diff]. It is symmetric.
;;     *   **`cohens-u3`**: Quantifies the proportion of values in the second group (`group2`) that are less than the median of the first group (`group1`). This is a non-parametric version of the concept behind `cohens-u3-normal`, but is calculated directly from the empirical distributions' medians and CDFs. Range is [0, 1]. This measure is also **asymmetric**.
;;

(utls/examples-note
  (stats/cohens-u1 virginica-sepal-length setosa-sepal-length)
  (stats/cohens-u2 virginica-sepal-length setosa-sepal-length)
  (stats/cohens-u3 virginica-sepal-length setosa-sepal-length))

;; **Comparison:**

;; - Use `p-overlap` for a robust, non-parametric, symmetric measure of direct area overlap based on estimated densities.
;; - Use `cohens-u*-normal` functions when you are comfortable assuming normality and equal variances, and want measures directly derived from Cohen's d. `u1`, `u2`, and `u3` offer slightly different views of overlap/separation relative to means and standard deviations.
;; - Use non-parametric `cohens-u*` functions when assumptions of normality are not met or when working with ordinal data. `u2` is a symmetric overlap measure based on quantiles, while `u3` is an asymmetric measure based on comparing to the median.

;; ### Correlation / Association Family

;; Effect size measures in the correlation and association family quantify the strength
;; and magnitude of the relationship between variables. Unlike difference-based measures
;; (like Cohen's d), these focus on how well one variable predicts or is associated
;; with another, often in terms of shared or explained variance. They are widely used
;; in regression, ANOVA, and general correlation analysis.

;; ::: {.callout-tip title="Defined functions"}
;; * `pearson-r`, `r2-determination`
;; * `eta-sq`, `omega-sq`, `epsilon-sq`
;; * `cohens-f2`, `cohens-f`, `cohens-q`
;; * `rank-eta-sq`, `rank-epsilon-sq`
;; :::

;; *   **Pearson r** (`pearson-r`): The Pearson product-moment correlation coefficient.
;;     An alias for [[pearson-correlation]]. Measures the strength and direction of a
;;     *linear* relationship between two continuous variables ($r$). Ranges from -1 to +1.
;;     $$ r = \frac{\sum(x_i - \bar{x})(y_i - \bar{y})}{\sqrt{\sum(x_i - \bar{x})^2 \sum(y_i - \bar{y})^2}} $$
;; *   **R² Determination** (`r2-determination`): The coefficient of determination.
;;     An alias for the standard [[r2]] calculation between two sequences. For a simple
;;     linear relationship, this is the square of the Pearson correlation ($R^2 = r^2$).
;;     It quantifies the proportion of the variance in one variable that is
;;     linearly predictable from the other. Ranges from 0 to 1.
;;     $$ R^2 = 1 - \frac{RSS}{TSS} = \frac{SS_{regression}}{SS_{total}} $$
;; *   **Eta-squared** (`eta-sq`): A measure of the proportion of the variance in
;;     a dependent variable that is associated with an independent variable. In the
;;     context of two numerical sequences, this function calculates $R^2$ from the
;;     simple linear regression of `group1` on `group2`.
;;     $$ \eta^2 = \frac{SS_{regression}}{SS_{total}} $$
;; *   **Omega-squared** (`omega-sq`): A less biased estimate of the population
;;     Eta-squared ($\omega^2$). It estimates the proportion of variance in the
;;     dependent variable accounted for by the independent variable *in the population*.
;;     Often preferred over Eta-squared as a population effect size estimate, especially for small sample sizes.
;;     $$ \omega^2 = \frac{SS_{regression} - p \cdot MSE}{SS_{total} + MSE} $$
;;     where $p$ is the number of predictors (1 for simple linear regression) and $MSE$ is the Mean Squared Error of the residuals.
;; *   **Epsilon-squared** (`epsilon-sq`): Another less biased estimate of population
;;     Eta-squared ($\varepsilon^2$), similar to adjusted $R^2$. Also aims to estimate
;;     the proportion of variance explained in the population.
;;     $$ \varepsilon^2 = \frac{SS_{regression} - MSE}{SS_{total}} $$
;; *   **Cohen's f²** (`cohens-f2`): Measures the ratio of the variance explained by
;;     the effect to the unexplained variance. It can be calculated based on Eta-squared,
;;     Omega-squared, or Epsilon-squared (controlled by the `:type` parameter). Ranges from 0 upwards.
;;     $$ f^2 = \frac{\text{Proportion of Variance Explained}}{\text{Proportion of Variance Unexplained}} = \frac{\eta^2}{1-\eta^2} $$
;; *   **Cohen's f** (`cohens-f`): The square root of Cohen's f² ($f = \sqrt{f^2}$).
;;     It quantifies the magnitude of the effect size in standardized units. Ranges from 0 upwards.

;; Let's look at the relationship between `mpg` (Miles Per Gallon) and `hp` (Horsepower)
;; from the `mtcars` dataset. We expect a negative relationship (higher HP means lower MPG).

;; The Pearson r value shows a strong negative linear correlation (-0.77). Squaring
;; this gives the R² (`r2-determination`, `eta-sq`), indicating that about 59% of the
;; variance in MPG can be explained by the linear relationship with HP in this sample.
;; `omega-sq` and `epsilon-sq` provide adjusted estimates for the population, which
;; are slightly lower, as expected. Cohen's f² and f quantify the magnitude of this
;; effect, suggesting a large effect size according to conventional guidelines.

(utls/examples-note
  (stats/pearson-r mpg hp)
  (stats/r2-determination mpg hp)
  (stats/eta-sq mpg hp)
  (stats/omega-sq mpg hp)
  (stats/epsilon-sq mpg hp)
  (stats/cohens-f2 mpg hp)
  (stats/cohens-f2 mpg hp :omega)
  (stats/cohens-f2 mpg hp :epsilon-sq)
  (stats/cohens-f mpg hp)
  (stats/cohens-f mpg hp :omega)
  (stats/cohens-f mpg hp :epsilon-sq))

;; *   **Cohen's q** (`cohens-q`): Measures the difference between two correlation
;;     coefficients after applying the Fisher z-transformation (`atanh`). Useful for
;;     comparing the strength of correlation between different pairs of variables or in different samples.
;;     $$ q = \text{atanh}(r_1) - \text{atanh}(r_2) $$

;; Let's compare correlations using `cohens-q`. We'll compare the correlation
;; between `mpg` and `hp` against the correlation between `mpg` and `wt` (weight).
;; Second test is to compare correlation between setosa sepal widths and lengths against correlation between virginical sepal widths and lengths.

(def r-mpg-hp (stats/pearson-correlation mpg hp))
(def r-mpg-wt (stats/pearson-correlation mpg wt))

(utls/examples-note
  r-mpg-hp
  r-mpg-wt
  (stats/cohens-q r-mpg-hp r-mpg-wt)
  (stats/cohens-q mpg hp wt)
  (stats/cohens-q setosa-sepal-width setosa-sepal-length
                  virginica-sepal-width virginica-sepal-length))

;; The `cohens-q` value quantifies the difference in the strength of these two
;; (dependent) correlations. A larger absolute value of `q` suggests a more
;; substantial difference between the correlation coefficients.

;; *   **Rank Eta-squared** (`rank-eta-sq`): Effect size measure for the Kruskal-Wallis
;;     test. It represents the proportion of variation in the dependent variable accounted
;;     for by group membership, based on ranks. Ranges from 0 to 1. Calculated based on
;;     the Kruskal-Wallis H statistic ($H$), number of groups ($k$), and total sample size ($n$).
;;     $$ \eta^2_H = \frac{\max(0, H - (k-1))}{n - k} $$
;; *   **Rank Epsilon-squared** (`rank-epsilon-sq`): Another effect size measure
;;     for the Kruskal-Wallis test, also based on ranks and providing a less biased
;;     estimate than Rank Eta-squared. Ranges from 0 to 1. Calculated based on $H$ and $n$.
;;     $$ \varepsilon^2_H = \frac{H}{n-1} $$

;; We can treat the different species (setosa, virginica and versicolor) as groups and compare
;; their sepal lengths using a non-parametric approach conceptually related to
;; Kruskal-Wallis (though the function itself is a generic effect size calculation).

(utls/examples-note
  (stats/rank-eta-sq (vals sepal-lengths))
  (stats/rank-eta-sq (vals sepal-widths))
  (stats/rank-epsilon-sq (vals sepal-lengths))
  (stats/rank-epsilon-sq (vals sepal-widths)))

;; **Comparison:**

;; *   `pearson-r` measures linear relationship direction and strength. `r2-determination`
;;     quantifies the proportion of shared variance ($r^2$) for a simple linear relationship.
;; *   `eta-sq`, `omega-sq`, and `epsilon-sq` all quantify the proportion of variance
;;     explained. `eta-sq` is a sample-based measure (equal to $R^2$ here), while
;;     `omega-sq` and `epsilon-sq` are less biased population estimates, preferred when
;;     generalizing beyond the sample.
;; *   `cohens-f2` and `cohens-f` measure effect magnitude relative to unexplained
;;     variance. They are useful for power analysis and interpreting the practical
;;     significance of a regression or ANOVA effect. `cohens-f` is often more
;;     interpretable as it's in the same units as standard deviations (though it's
;;     not directly scaled by a single standard deviation like Cohen's d).
;; *   `cohens-q` is specifically for comparing correlation coefficients, assessing
;;     if the strength of relationship between two variables differs significantly
;;     from the strength between another pair or in another context.
;; *   `rank-eta-sq` and `rank-epsilon-sq` are specific to rank-based tests like Kruskal-Wallis,
;;     providing effect sizes analogous to Eta-squared/Epsilon-squared but suitable for
;;     non-parametric data.

;; **Interpretation Guidelines:**

;; *   **Correlation (r):**
;;     *   $|r| < 0.3$: weak/small linear relationship
;;     *   $0.3 \le |r| < 0.7$: moderate/medium linear relationship
;;     *   $|r| \ge 0.7$: strong/large linear relationship
;; *   **Proportion of Variance Explained ($R^2, \eta^2, \omega^2, \varepsilon^2$):**
;;     *   $0.01$: small effect (1% of variance explained)
;;     *   $0.06$: medium effect (6% of variance explained)
;;     *   $0.14$: large effect (14% of variance explained)
;; *   **Cohen's f (and f²):**
;;     *   $f = 0.1$ ($f^2 = 0.01$): small effect
;;     *   $f = 0.25$ ($f^2 \approx 0.06$): medium effect
;;     *   $f = 0.4$ ($f^2 = 0.16$): large effect

;; ## Statistical Tests

;; Functions for hypothesis testing.

;; ### Normality and Shape Tests

;; These tests assess whether a dataset deviates significantly from a normal distribution, or exhibits specific characteristics related to skewness (asymmetry) and kurtosis (tailedness). Deviations from normality are important to check as many statistical methods assume normal data.

;; ::: {.callout-tip title="Defined functions"}
;; * `skewness-test`
;; * `kurtosis-test`
;; * `normality-test`
;; * `jarque-bera-test`
;; * `bonett-seier-test`
;; :::

;; `fastmath.stats` provides several functions for these tests:
;;
;; *   `skewness-test`: Tests if the sample skewness significantly differs from the zero skewness expected of a normal distribution. It calculates a standardized test statistic (approximately Z) and its p-value. By default, it uses the `:g1` type of skewness (Pearson's moment coefficient).
;; *   `kurtosis-test`: Tests if the sample kurtosis significantly differs from the value expected of a normal distribution (3 for raw kurtosis, 0 for excess kurtosis). It also calculates a standardized test statistic (approximately Z) and its p-value. By default, it uses the `:kurt` type of kurtosis (Excess Kurtosis + 3).
;; *   `normality-test`: The D'Agostino-Pearson K² test. This is an omnibus test that combines the skewness and kurtosis tests into a single statistic ($K^2 = Z_{skewness}^2 + Z_{kurtosis}^2$). $K^2$ follows approximately a Chi-squared distribution with 2 degrees of freedom under the null hypothesis of normality. It tests for overall departure from normality. It internally uses the default `:g1` skewness and default `:kurt` kurtosis types from the individual tests.
;; *   `jarque-bera-test`: Another omnibus test for normality based on sample skewness ($S$, specifically type `:g1`) and excess kurtosis ($K$, specifically type `:g2`). The test statistic is $JB = \frac{n}{6}(S^2 + \frac{1}{4}K^2)$, which also follows approximately a Chi-squared distribution with 2 degrees of freedom under the null hypothesis. Similar to the K² test, it assesses whether the combined deviation in skewness and kurtosis is significant. Note that the specific `:g1` and `:g2` types are used *for the JB statistic formula* regardless of any `:type` option passed to the underlying `skewness` or `kurtosis` functions.
;; *   `bonett-seier-test`: A test for normality based on Geary's 'g' measure of kurtosis, which is more robust to outliers than standard moment-based kurtosis. It calculates a Z-statistic comparing the sample Geary's 'g' to its expected value for a normal distribution ($\sqrt{2/\pi}$). This test specifically uses the `:geary` kurtosis type.
;; 
;; **Interpretation:** The output of these functions is typically a map containing:
;; 
;; *   `:stat` (or an alias like `:Z`, `:K2`, `:JB`, `:F`, `:chi2`): The calculated value of the test statistic.
;; *   `:p-value`: The probability of observing a test statistic as extreme as, or more extreme than, the one calculated, assuming the null hypothesis (that the data comes from a normal distribution, or has the expected shape property) is true.
;; *   Other keys might include `:df` (degrees of freedom), `:skewness` (the value *of the type used in the test*), `:kurtosis` (the value *of the type used in the test*), `:n` (sample size), and `:sides` (the alternative hypothesis used).
;; 
;; A small p-value (typically less than the chosen significance level, e.g., 0.05) suggests that the observed data is unlikely to have come from a normal distribution (or satisfy the specific shape property being tested), leading to rejection of the null hypothesis. A large p-value indicates that there is not enough evidence to reject the null hypothesis.
;; 
;; Let's apply some of these tests to the `residual-sugar` data, which we observed earlier to be right-skewed. We'll show the specific skewness/kurtosis values used by each test for clarity.

(utls/examples-note
  (stats/skewness residual-sugar :g1) ;; Skewness used by skewness-test, normality-test, jarque-bera-test
  (stats/kurtosis residual-sugar :kurt) ;; Kurtosis used by kurtosis-test, normality-test
  (stats/kurtosis residual-sugar :g2) ;; Excess Kurtosis used by jarque-bera-test
  (stats/kurtosis residual-sugar :geary) ;; Kurtosis used by bonett-seier-test
  
  (stats/skewness-test residual-sugar)
  (stats/kurtosis-test residual-sugar)
  (stats/normality-test residual-sugar)
  (stats/jarque-bera-test residual-sugar)
  (stats/bonett-seier-test residual-sugar))

;; As expected for the right-skewed `residual-sugar` data, the `skewness-test`, `normality-test` (K²), and `jarque-bera-test` yield very small p-values, strongly suggesting that the data is not normally distributed. The `kurtosis-test` and `bonett-seier-test` examine tailedness; their p-values indicate whether the deviation from normal kurtosis is significant.

;; ### Binomial Tests

;; These tests are designed for data representing counts of successes in a fixed number of trials, where each trial has only two possible outcomes (success or failure). They are used to make inferences about the true underlying proportion of successes in the population.

;; ::: {.callout-tip title="Defined functions"}
;; * `binomial-test`
;; * `binomial-ci`
;; :::

;; **`binomial-test`**: Performs an exact hypothesis test on a binomial proportion. This test is used to determine if the observed number of successes in a given number of trials is significantly different from what would be expected under a specific hypothesized probability of success ($p_0$). The null hypothesis is typically $H_0: p = p_0$.
;;
;;     The test calculates a p-value based on the binomial probability distribution, assessing the likelihood of observing the sample results (or more extreme) if $p_0$ were the true population proportion.

;; Let's illustrate with some examples. Imagine we check the `am` column of the `mtcars` dataset, where 0 represents automatic and 1 represents manual transmission. We want to test if the proportion of cars with manual transmission is significantly different from 0.5 (an equal split). We also want to estimate a confidence interval for this proportion.

;; Let's get the counts:

(def mtcars-am (ds/mtcars :am))
(def manual-count (count (filter m/one? mtcars-am)))
(def total-count (count mtcars-am))

(utls/examples-note
  manual-count
  total-count
  (stats/mean mtcars-am))

;; Now, let's perform the binomial test:

(utls/zp (stats/binomial-test manual-count total-count))
(utls/zp (stats/binomial-test manual-count total-count {:p 0.5 :sides :one-sided-greater}))
(utls/zp (stats/binomial-test manual-count total-count {:p 0.8 :sides :one-sided-less :alpha 0.01}))

;; The output map from `binomial-test` contains:

;; *   `:p-value`: The probability of observing the sample result or more extreme outcomes if the true population proportion were equal to `:p`. A small p-value (typically < 0.05) suggests evidence against the null hypothesis.
;; *   `:stat`: The observed number of successes.
;; *   `:estimate`: The observed proportion of successes (`:successes / :trials`).
;; *   `:confidence-interval`: A confidence interval for the true population proportion, based on the observed data and using the specified `:ci-method` and `:alpha`.

;; #### Binomial Confidence Intervals

;; **`binomial-ci`**: Calculates a confidence interval for a binomial proportion. Given the observed number of successes and trials, this function estimates a range of values that is likely to contain the true population probability of success ($p$) with a specified level of confidence.
;;
;; The `binomial-ci` function offers various methods for calculating the confidence interval for a binomial proportion, controlled by the optional `method` keyword. These methods differ in their underlying assumptions and formulas, leading to intervals that can vary in width and coverage properties, particularly with small sample sizes or proportions close to 0 or 1.

;; Available `method` values:

;; *   `:asymptotic`: **Normal Approximation (Wald) Interval**. Based on the Central Limit Theorem, using the sample proportion and its estimated standard error. Simple to calculate but can have poor coverage (too narrow) for small sample sizes or proportions near 0 or 1.
;; *   `:agresti-coull`: **Agresti-Coull Interval**. An adjusted Wald interval that adds 'pseudo-counts' (typically 2 successes and 2 failures) to the observed counts. This adjustment improves coverage and performance, especially for small samples.
;; *   `:clopper-pearson`: **Clopper-Pearson Interval**. An 'exact' method based on inverting binomial tests. It provides guaranteed minimum coverage probability (i.e., the true proportion is included in the interval at least 100 * (1-alpha)% of the time). However, it is often wider than necessary and can be overly conservative.
;; *   `:wilson`: **Wilson Score Interval**. Derived from the score test, which is based on the null hypothesis standard error rather than the observed standard error. It performs well across different sample sizes and probabilities and is generally recommended over the Wald interval.
;; *   `:prop.test`: **Continuity-Corrected Wald Interval**. Applies a continuity correction to the Wald interval, similar to what is often used with the Chi-squared test for 2x2 tables.
;; *   `:cloglog`: **Complementary Log-log Transformation Interval**. Calculates the interval on the clog-log scale and then transforms it back to the probability scale. Can be useful when dealing with skewed data or probabilities close to 0.
;; *   `:logit`: **Logit Transformation Interval**. Calculates the interval on the logit scale (`log(p/(1-p))`) and transforms it back. Also useful for handling probabilities near boundaries.
;; *   `:probit`: **Probit Transformation Interval**. Calculates the interval on the probit scale (inverse of the standard normal CDF) and transforms it back.
;; *   `:arcsine`: **Arcsine Transformation Interval**. Calculates the interval using the arcsine square root transformation (`asin(sqrt(p))`) and transforms it back.
;; *   `:all`: Calculates and returns a map of intervals from all the above methods.

(utls/examples-note
  (stats/binomial-ci manual-count total-count)
  (stats/binomial-ci manual-count total-count :asymptotic 0.01)
  (stats/binomial-ci manual-count total-count :asymptotic 0.1)
  (stats/binomial-ci manual-count total-count :agresti-coull)
  (stats/binomial-ci manual-count total-count :clopper-pearson)
  (stats/binomial-ci manual-count total-count :wilson)
  (stats/binomial-ci manual-count total-count :prop.test)
  (stats/binomial-ci manual-count total-count :cloglog)
  (stats/binomial-ci manual-count total-count :logit)
  (stats/binomial-ci manual-count total-count :probit)
  (stats/binomial-ci manual-count total-count :arcsine))

;; For methods other than `:all`, the function returns a vector `[lower-bound, upper-bound, estimated-p]`.
;; *   `lower-bound`, `upper-bound`: The ends of the calculated confidence interval. We are 100 * (1 - alpha)% confident that the true population proportion lies within this range.
;; *   `estimated-p`: The observed proportion of successes in the sample.
;;
;; If the method is `:all`, a map is returned where keys are the method keywords and values are the `[lower, upper, estimate]` vectors for each method. Note that different methods can produce different interval widths and positions.

(utls/zp
 (stats/binomial-ci manual-count total-count :all))

;; ### Location Tests (T/Z Tests)

;; Location tests, specifically t-tests and z-tests, are fundamental statistical tools used to compare means. They help determine if the mean of a sample is significantly different from a known or hypothesized value (one-sample tests), or if the means of two different samples are significantly different from each other (two-sample tests). The choice between a t-test and a z-test typically depends on whether the population standard deviation is known and the sample size.

;; ::: {.callout-tip title="Defined functions"}
;; * `t-test-one-sample`, `z-test-one-sample`
;; * `t-test-two-samples`, `z-test-two-samples`
;; * `p-value`
;; :::

;; `fastmath.stats` provides the following functions for location tests:
;;
;; *   `t-test-one-sample`: Performs a one-sample Student's t-test. This tests the null hypothesis that the true population mean ($\mu$) is equal to a hypothesized value ($\mu_0$). It is typically used when the population standard deviation is unknown and estimated from the sample. The test statistic follows a t-distribution with $n-1$ degrees of freedom, where $n$ is the sample size.
;;
;;     $$ t = \frac{\bar{x} - \mu_0}{s/\sqrt{n}} $$
;;
;; *   `z-test-one-sample`: Performs a one-sample Z-test. This tests the null hypothesis that the true population mean ($\mu$) is equal to a hypothesized value ($\mu_0$). It is typically used when the population standard deviation is known or when the sample size is large (generally $n > 30$), allowing the sample standard deviation to be used as a reliable estimate for the population standard deviation, and the test statistic approximates a standard normal distribution.
;;
;;     $$ z = \frac{\bar{x} - \mu_0}{\sigma/\sqrt{n}} \quad \text{or} \quad z \approx \frac{\bar{x} - \mu_0}{s/\sqrt{n}} $$
;;
;; *   `t-test-two-samples`: Performs a two-sample t-test to compare the means of two samples. It can perform:
;;     *   **Unpaired tests**: For independent samples. By default, it performs Welch's t-test (`:equal-variances? false`), which does not assume equal population variances. It can also perform Student's t-test (`:equal-variances? true`), assuming equal variances and using a pooled standard deviation.
;;     *   **Paired tests**: For dependent samples (`:paired? true`), essentially performing a one-sample t-test on the differences between pairs.
;;     The null hypothesis is typically that the true difference between population means is zero, or equal to a hypothesized value ($\mu_0$). The test statistic follows a t-distribution with degrees of freedom calculated appropriately for the specific variant (pooled for Student's, Satterthwaite approximation for Welch's, $n-1$ for paired).
;;
;; *   `z-test-two-samples`: Performs a two-sample Z-test to compare the means of two independent or paired samples. Similar to the one-sample Z-test, it is used when population variances are known or samples are large. It can handle independent samples (with or without assuming equal variances, affecting the standard error calculation) or paired samples (by performing a one-sample Z-test on differences). The test statistic approximates a standard normal distribution.
;;
;; **Comparison:**

;; *   **T vs Z**: Choose Z-tests primarily when population standard deviations are known or with large samples (often $n > 30$ for each group in two-sample tests), as the sampling distribution of the mean is well-approximated by the normal distribution. Use T-tests when population standard deviations are unknown and estimated from samples, especially with small to moderate sample sizes, as the sampling distribution is better described by the t-distribution.
;; *   **One-Sample vs Two-Sample**: Use one-sample tests to compare a single sample mean against a known or hypothesized constant. Use two-sample tests to compare the means of two distinct samples.
;; *   **Paired vs Unpaired (Two-Sample)**: Use paired tests when the two samples consist of paired observations (e.g., measurements before and after an intervention on the same subjects). Use unpaired tests for independent samples (e.g., comparing a treatment group to a control group with different subjects in each).
;; *   **Welch's vs Student's (Unpaired Two-Sample T-test)**: Welch's test is generally recommended as it does not assume equal variances, a common violation in practice. Student's t-test requires the assumption of equal population variances.
;;
;; All these functions return a map containing key results, typically including the test statistic (`:t` or `:z`), `:p-value`, `:confidence-interval` for the mean (or mean difference), `:estimate` (sample mean or mean difference), sample size(s) (`:n`, `:nx`, `:ny`), `:mu` (hypothesized value), `:stderr` (standard error of the estimate), `:alpha` (significance level), and `:sides` (alternative hypothesis type). Two-sample tests also report `:paired?` and, if unpaired, `:equal-variances?`.
;;
;; **Arguments:**
;;
;; *   `xs` (sequence of numbers): The sample data for one-sample tests or the first sample for two-sample tests.
;; *   `ys` (sequence of numbers): The second sample for two-sample tests.
;; *   `params` (map, optional): A map of options. Common keys:
;;     *   `:alpha` (double, default `0.05`): Significance level. Confidence level is $1 - \alpha$.
;;     *   `:sides` (keyword, default `:two-sided`): Alternative hypothesis. Can be `:two-sided`, `:one-sided-greater`, or `:one-sided-less`.
;;     *   `:mu` (double, default `0.0`): Hypothesized population mean (one-sample) or hypothesized difference in population means (two-sample).
;;     *   `:paired?` (boolean, default `false`): For two-sample tests, specifies if samples are paired.
;;     *   `:equal-variances?` (boolean, default `false`): For unpaired two-sample tests, specifies if equal population variances are assumed (Student's vs Welch's).
;;
;; **Return Value (Map):**
;;
;; *   `:stat` (double): The calculated test statistic (alias `:t` or `:z`).
;; *   `:p-value` (double): The probability of observing the test statistic or more extreme values under the null hypothesis.
;; *   `:confidence-interval` (vector of doubles): `[lower-bound, upper-bound]`. A range likely containing the true population mean (one-sample) or mean difference (two-sample) with $1-\alpha$ confidence. Includes the estimate as a third value (e.g., `[lower upper estimate]`).
;; *   `:estimate` (double): The sample mean (one-sample) or the difference between sample means (two-sample).
;; *   `:n` (long or vector of longs): Sample size (one-sample) or sample sizes `[nx ny]` (two-sample).
;; *   `:nx`, `:ny` (long): Sample sizes of `xs` and `ys` respectively (two-sample unpaired).
;; *   `:estimated-mu` (vector of doubles): Sample means `[mean xs, mean ys]` (two-sample unpaired).
;; *   `:mu` (double): The hypothesized value used in the null hypothesis.
;; *   `:stderr` (double): The standard error of the sample mean or mean difference.
;; *   `:alpha` (double): The significance level used.
;; *   `:sides` / `:test-type` (keyword): The alternative hypothesis side.
;; *   `:df` (long or vector of longs): Degrees of freedom (t-tests only).
;; *   `:paired?` (boolean): Indicates if a paired test was performed (two-sample only).
;; *   `:equal-variances?` (boolean): Indicates if equal variances were assumed (two-sample unpaired only).
;;
;; **Helper function `p-value`:**
;;
;; The `p-value` function is a general utility used internally by many statistical tests (including t-tests and z-tests) to calculate the p-value from a given test statistic, its null distribution, and the specified alternative hypothesis (sides). It determines the probability of observing a statistic value as extreme as or more extreme than the one obtained from the sample, assuming the null hypothesis is true. The interpretation of 'extreme' depends on whether a two-sided, one-sided-greater, or one-sided-less test is performed. For discrete distributions, a continuity correction is applied.
;;
;; Let's apply these tests. First, one-sample tests on `mpg` data against a hypothesized mean of 20.0:

(utls/zp (stats/t-test-one-sample mpg {:mu 20.0}))
(utls/zp (stats/z-test-one-sample mpg {:mu 20.0 :sides :one-sided-greater}))

;; Now, two-sample tests comparing `setosa-sepal-length` and `virginica-sepal-length`. We'll use unpaired tests first, comparing Welch's and Student's:

(utls/zp (stats/t-test-two-samples setosa-sepal-length virginica-sepal-length))
(utls/zp (stats/t-test-two-samples setosa-sepal-length virginica-sepal-length {:equal-variances? true}))

;; And the corresponding Z-tests (note the Z-test doesn't need degrees of freedom):

(utls/zp (stats/z-test-two-samples setosa-sepal-length virginica-sepal-length))
(utls/zp (stats/z-test-two-samples setosa-sepal-length virginica-sepal-length {:equal-variances? true}))

;; If we had paired data (e.g., `hp` measurements before and after a modification for the same cars, stored in `hp-before` and `hp-after`), we would use the paired option:

;; Assume `hp-before` and `hp-after` are sequences of the same length

(def hp-before [100 120 150])
(def hp-after [110 125 165])
(utls/zp (stats/t-test-two-samples hp-before hp-after {:paired? true}))
(utls/zp (stats/z-test-two-samples hp-before hp-after {:paired? true}))

;; ### Variance Tests

;; Variance tests are used to assess whether the variances of two or more independent samples are statistically different. The null hypothesis for these tests is typically that the population variances are equal (homogeneity of variances, or homoscedasticity). This assumption is important for many statistical procedures, such as ANOVA and Student's t-test.

;; ::: {.callout-tip title="Defined functions"}
;; * `f-test`
;; * `levene-test`
;; * `brown-forsythe-test`
;; * `fligner-killeen-test`
;; :::

;; *   `f-test`: Performs an F-test for the equality of variances of **two** independent samples.
;;
;;     It calculates the ratio of the sample variances ($s_1^2 / s_2^2$) which follows an F-distribution under the null hypothesis of equal population variances, assuming the data in both groups are normally distributed.
;;
;;     $$ F = \frac{s_1^2}{s_2^2} $$
;;
;;     Parameters:
;;     - `xs` (seq of numbers): The first sample.
;;     - `ys` (seq of numbers): The second sample.
;;     - `params` (map, optional): Options map with `:sides` (default `:two-sided`) and `:alpha` (default `0.05`).
;;
;;     Returns a map with keys like `:F` (the F-statistic), `:p-value`, `:df` ([numerator-df, denominator-df]), `:estimate` (the variance ratio), `:n` ([nx, ny]), and `:sides`.
;;
;;     Assumptions: Independent samples, normality within each sample. Sensitive to departures from normality.
;;
;; *   `levene-test`: Performs Levene's test for homogeneity of variances across **two or more** independent groups.
;;
;;     This test is more robust to departures from normality than the F-test. It performs a one-way ANOVA on the absolute deviations of each observation from its group's **mean**.
;;
;;     $$ W = \frac{N-k}{k-1} \frac{\sum_{i=1}^k n_i (\bar{Z}_{i\cdot} - \bar{Z}_{\cdot\cdot})^2}{\sum_{i=1}^k \sum_{j=1}^{n_i} (Z_{ij} - \bar{Z}_{i\cdot})^2} $$
;;
;;     where $Z_{ij} = |X_{ij} - \bar{X}_{i\cdot}|$ are the absolute deviations from the group means, $N$ is total sample size, $k$ is number of groups, $n_i$ is size of group $i$, $\bar{Z}_{i\cdot}$ is mean of absolute deviations in group $i$, and $\bar{Z}_{\cdot\cdot}$ is the grand mean of absolute deviations.
;;
;;     Parameters:
;;     - `xss` (sequence of sequences): Collection of groups.
;;     - `params` (map, optional): Options map with `:sides` (default `:one-sided-greater`), `:statistic` (default [[mean]]), and `:scorediff` (default [[abs]]).
;;
;;     Returns a map with keys like `:W` (the F-statistic for the ANOVA on deviations), `:p-value`, `:df` ([DFt, DFe]), `:n` (group sizes), and standard ANOVA output keys.
;;
;;     Assumptions: Independent samples. Less sensitive to non-normality than F-test.
;;
;; *   `brown-forsythe-test`: Performs the Brown-Forsythe test, a modification of Levene's test using the **median** as the center.
;;
;;     This version is even more robust to non-normality than the standard Levene's test. It performs an ANOVA on the absolute deviations from group **medians**.
;;
;;     $$ F^* = \frac{N-k}{k-1} \frac{\sum_{i=1}^k n_i (\tilde{Z}_{i\cdot} - \tilde{Z}_{\cdot\cdot})^2}{\sum_{i=1}^k \sum_{j=1}^{n_i} (Z_{ij} - \tilde{Z}_{i\cdot})^2} $$
;;
;;     where $Z_{ij} = |X_{ij} - \tilde{X}_{i\cdot}|$ are the absolute deviations from the group medians, $\tilde{Z}_{i\cdot}$ is the median of absolute deviations in group $i$, and $\tilde{Z}_{\cdot\cdot}$ is the grand median of absolute deviations.
;;
;;     Parameters:
;;     - `xss` (sequence of sequences): Collection of groups.
;;     - `params` (map, optional): Options map with `:sides` (default `:one-sided-greater`) and `:scorediff` (default [[abs]]). Internally sets `:statistic` to [[median]] in [[levene-test]].
;;
;;     Returns a map similar to `levene-test`, with an F-statistic derived from the ANOVA on median deviations.
;;
;;     Assumptions: Independent samples. Most robust to non-normality among parametric-based tests.
;;
;; *   `fligner-killeen-test`: Performs the Fligner-Killeen test for homogeneity of variances across **two or more** independent groups.
;;
;;     This is a non-parametric test based on ranks of the absolute deviations from group medians. It is generally considered one of the most robust tests for homogeneity of variances when assumptions about the underlying distribution are questionable.
;;
;;     The test statistic is based on the ranks of $|X_{ij} - \tilde{X}_{i\cdot}|$, where $\tilde{X}_{i\cdot}$ is the median of group $i$.
;;
;;     Parameters:
;;     - `xss` (sequence of sequences): Collection of groups.
;;     - `params` (map, optional): Options map with `:sides` (default `:one-sided-greater`).
;;
;;     Returns a map with keys like `:chi2` (the Chi-squared statistic), `:p-value`, `:df` (number of groups - 1), `:n` (group sizes), and standard ANOVA output keys (calculated on transformed ranks).
;;
;;     Assumptions: Independent samples. Non-parametric, robust to non-normality. Requires distributions to have similar shape for valid inference.
;;
;; Comparison:
;; - The F-test is simple but highly sensitive to non-normality. Use cautiously if data isn't clearly normal.
;; - Levene's and Brown-Forsythe tests are more robust to non-normality. Brown-Forsythe (median-based) is generally preferred over standard Levene's (mean-based) when distributions are potentially heavy-tailed or skewed.
;; - Fligner-Killeen is a non-parametric alternative and the most robust to non-normality, based on ranks.
;;
;; The choice depends on sample size, suspected distribution shape, and the need for robustness.
;;

;; **Function-specific keys:**
;;
;; *   `f-test` additionally returns:
;;     *   `:F`: The calculated F-statistic (ratio of sample variances, Var(xs) / Var(ys)). Alias for `:stat`.
;;     *   `:nx`, `:ny`: Sample sizes of the first (`xs`) and second (`ys`) groups.
;;     *   `:estimate`: The ratio of the sample variances (same value as `:F`).
;;     *   `:confidence-interval`: A confidence interval for the true ratio of population variances (Var(xs) / Var(ys)).
;; *   `levene-test` and `brown-forsythe-test` additionally return:
;;     *   `:W`: The calculated test statistic (which is an F-statistic derived from the ANOVA on deviations). Alias for `:stat`.
;;     *   `:F`: Alias for `:W` (as it's an F-statistic).
;;     *   `:SSt`: Sum of Squares Treatment (between groups of deviations).
;;     *   `:SSe`: Sum of Squares Error (within groups of deviations).
;;     *   `:DFt`: Degrees of Freedom Treatment (number of groups - 1).
;;     *   `:DFe`: Degrees of Freedom Error (total sample size - number of groups).
;;     *   `:MSt`: Mean Square Treatment (SSt / DFt).
;;     *   `:MSe`: Mean Square Error (SSe / DFe).
;; *   `fligner-killeen-test` additionally returns:
;;     *   `:chi2`: The calculated Chi-squared statistic. Alias for `:stat`.
;;     *   `:SSt`: Sum of Squares Treatment, calculated from transformed ranks.
;;     *   `:SSe`: Sum of Squares Error, calculated from transformed ranks.
;;     *   `:DFt`: Degrees of Freedom Treatment (number of groups - 1).
;;     *   `:DFe`: Degrees of Freedom Error (total sample size - number of groups).
;;     *   `:MSt`: Mean Square Treatment (SSt / DFt).
;;     *   `:MSe`: Mean Square Error (SSe / DFe).


;; Let's apply these tests to compare variances of sepal lengths for 'setosa' and 'virginica' species, and then across all three species (including 'versicolor' from the original `sepal-lengths` map).

(utls/zp (stats/f-test setosa-sepal-length virginica-sepal-length))

(utls/zp (stats/levene-test (vals sepal-lengths)))

(utls/zp (stats/brown-forsythe-test (vals sepal-lengths)))

(utls/zp (stats/fligner-killeen-test (vals sepal-lengths)))

;; We can also compare sepal widths across species using these tests.

(utls/zp (stats/f-test (sepal-widths :setosa) (sepal-widths :virginica)))
(utls/zp (stats/levene-test (vals sepal-widths)))
(utls/zp (stats/brown-forsythe-test (vals sepal-widths)))
(utls/zp (stats/fligner-killeen-test (vals sepal-widths)))

;; ### Goodness-of-Fit and Independence Tests

;; These tests are used to evaluate how well an observed distribution of data matches a hypothesized theoretical distribution (Goodness-of-Fit) or whether there is a statistically significant association between two or more categorical variables (Independence). They are commonly applied to categorical data or data that has been grouped into categories (e.g., histograms).

;; ::: {.callout-tip title="Defined functions"}
;; * `power-divergence-test`
;; * `chisq-test`
;; * `multinomial-likelihood-ratio-test`
;; * `minimum-discrimination-information-test`
;; * `neyman-modified-chisq-test`
;; * `freeman-tukey-test`
;; * `cressie-read-test`
;; * `ad-test-one-sample`
;; * `ks-test-one-sample`, `ks-test-two-samples`
;; :::

;; **Power Divergence Tests:**
;; 
;; The `power-divergence-test` is a generalized framework for several statistical tests that compare observed frequencies to expected frequencies. The specific test performed is determined by the `lambda` parameter.
;; 
;; The general test statistic is:
;; 
;; $$ 2 \sum_i O_i \left( \frac{\left(\frac{O_i}{E_i}\right)^\lambda - 1}{\lambda(\lambda+1)} \right) $$
;; 
;; where $O_i$ are the observed counts, $E_i$ are the expected counts, and $\lambda$ is the power parameter.
;; 
;; This test can be used for two main purposes:
;; 
;; *   **Goodness-of-Fit:** Testing if a sequence of observed counts matches a set of expected counts (proportions/weights) or if a dataset matches a theoretical distribution.
;; *   **Independence:** Testing if there is an association between categorical variables in a contingency table.
;; 
;; The test statistic approximately follows a Chi-squared distribution with degrees of freedom determined by the specific application (e.g., number of categories - 1 - number of estimated parameters for GOF, $(rows-1)(cols-1)$ for independence).
;; 
;; `fastmath.stats` provides the general `power-divergence-test` and aliases for common `lambda` values:
;; 
;; *   `chisq-test`: Pearson's Chi-squared test ($\lambda=1$). The most common test in this family.
;;     $$ \chi^2 = \sum_i \frac{(O_i - E_i)^2}{E_i} $$
;; *   `multinomial-likelihood-ratio-test`: G-test ($\lambda=0$). Based on the ratio of likelihoods under the null and alternative hypotheses.
;;     $$ G = 2 \sum_i O_i \ln\left(\frac{O_i}{E_i}\right) $$
;; *   `minimum-discrimination-information-test`: Minimum Discrimination Information test ($\lambda=-1$). Also known as the G-test on expected counts vs observed.
;;     $$ I = 2 \sum_i E_i \ln\left(\frac{E_i}{O_i}\right) $$
;; *   `neyman-modified-chisq-test`: Neyman Modified Chi-squared test ($\lambda=-2$).
;;     $$ NM = \sum_i \frac{(O_i - E_i)^2}{O_i} $$
;; *   `freeman-tukey-test`: Freeman-Tukey test ($\lambda=-0.5$).
;;     $$ FT = \sum_i (\sqrt{O_i} - \sqrt{E_i})^2 $$
;; *   `cressie-read-test`: Cressie-Read test ($\lambda=2/3$, default for `power-divergence-test`). A compromise test.
;; 
;; **Arguments:**
;; 
;; *   `contingency-table-or-xs`: For independence tests, a contingency table (sequence of sequences or map). For goodness-of-fit with counts, a sequence of observed counts. For goodness-of-fit with data, a sequence of raw data.
;; *   `params` (map, optional): Options including:
;;     *   `:lambda` (double): Power parameter (default 2/3).
;;     *   `:p`: For GOF, expected probabilities/weights (seq) or a `fastmath.random` distribution. Ignored for independence tests.
;;     *   `:sides` (keyword, default `:one-sided-greater`): For p-value calculation.
;;     *   `:alpha` (double, default 0.05): For confidence intervals.
;;     *   `:ci-sides` (keyword, default `:two-sided`): For confidence intervals.
;;     *   `:bootstrap-samples` (long, default 1000): For bootstrap CI.
;;     *   `:ddof` (long, default 0): Delta degrees of freedom subtracted from the calculated DF.
;;     *   `:bins`: For GOF with data vs distribution, histogram bins (see [[histogram]]).
;; 
;; **Returns:** A map with keys:
;; 
;; *   `:stat` (or `:chi2`): The calculated test statistic.
;; *   `:df`: Degrees of freedom.
;; *   `:p-value`: P-value.
;; *   `:n`: Total number of observations.
;; *   `:estimate`: Observed proportions.
;; *   `:expected`: Expected counts or proportions.
;; *   `:confidence-interval`: Bootstrap confidence intervals for observed proportions.
;; *   `:lambda`, `:alpha`, `:sides`, `:ci-sides`: Input parameters.
;; 
;; **Examples:**
;; 
;; Goodness-of-Fit test comparing observed counts to expected proportions (e.g., testing if a six-sided die is fair based on 60 rolls).
;; 

(def observed-rolls [10 12 8 11 9 10])
(def expected-proportions [1 1 1 1 1 1]) ;; Using weights proportional to expected, they will be normalized

(utls/zp
 (stats/chisq-test observed-rolls {:p expected-proportions}))

;; Goodness-of-Fit test comparing sample data (`mpg`) to a theoretical distribution (Normal).
;; This implicitly creates a histogram from `mpg` and compares its counts to the counts expected from a Normal distribution within those bins.
;; 

(utls/zp
 (stats/chisq-test mpg {:p (r/distribution :normal {:mu (stats/mean mpg) :sd (stats/stddev mpg)})}))

;; Independence test using a contingency table (from the previous section).
;; 

(utls/zp
 (stats/chisq-test ct-data-1))

(utls/zp
 (stats/multinomial-likelihood-ratio-test ct-data-2))

;; **Distribution Comparison Tests (AD/KS):**
;; 
;; These tests directly compare empirical cumulative distribution functions (ECDFs) or density estimates, providing non-parametric ways to assess goodness-of-fit or compare two samples without strong assumptions about the underlying distribution shape.
;; 
;; *   `ad-test-one-sample`: Anderson-Darling test. Primarily a Goodness-of-Fit test. It is particularly sensitive to differences in the tails of the distributions being compared. It tests if a sample `xs` comes from a specified distribution or an empirical distribution estimated from `ys`.
;;     $$ A^2 = -n - \sum_{i=1}^n \frac{2i-1}{n} [\ln(F(X_i)) + \ln(1-F(X_{n-i+1}))] $$
;;     where $X_i$ are the ordered data, $n$ is the sample size, and $F$ is the CDF of the hypothesized distribution.
;;     
;;     **Arguments:** `xs` (sample), `distribution-or-ys` (reference distribution or sample), `opts` (map: `:sides`, `:kernel`, `:bandwidth`).
;;
;;     **Returns:** Map with `:A2`, `:stat`, `:p-value`, `:n`, `:mean`, `:stddev`, `:sides`.
;; 
;; *   `ks-test-one-sample`: One-sample Kolmogorov-Smirnov test. Compares the ECDF of a sample `xs` to a theoretical CDF or the ECDF of another sample `ys`. It is sensitive to the largest vertical difference between the two CDFs.
;;     $$ D = \max_i(|F_n(X_i) - F(X_i)|) $$
;;     where $F_n$ is the ECDF of the sample, $F$ is the reference CDF, and $X_i$ are the ordered data.
;;     
;;     **Arguments:** `xs` (sample), `distribution-or-ys` (reference distribution or sample), `opts` (map: `:sides`, `:kernel`, `:bandwidth`, `:distinct?`).
;;
;;     **Returns:** Map with `:n`, `:dp`, `:dn`, `:d`, `:stat`, `:KS`, `:p-value`, `:sides`.
;; 
;; *   `ks-test-two-samples`: Two-sample Kolmogorov-Smirnov test. Compares the ECDFs of two independent samples `xs` and `ys`. It is sensitive to the largest vertical difference between the two ECDFs.
;;     $$ D_{n,m} = \max_i(|F_n(X_i) - G_m(X_i)|) $$
;;     where $F_n$ and $G_m$ are the ECDFs of the two samples.
;;     
;;     **Arguments:** `xs` (sample 1), `ys` (sample 2), `opts` (map: `:method` (exact/approx), `:sides`, `:distinct?`, `:correct?`).
;;
;;     **Returns:** Map with `:nx`, `:ny`, `:n`, `:dp`, `:dn`, `:d`, `:stat`, `:KS`, `:p-value`, `:sides`, `:method`.
;; 
;; **Notes:**
;; 
;; *   The AD test is generally more sensitive to differences in the tails of the distributions than the KS test.
;; *   The KS test (both one-sample and two-sample) is based on the maximum difference between CDFs, making it sensitive to any point where the distributions diverge, but perhaps less so specifically in the extreme tails compared to AD.
;; *   The KS test has issues with discrete data and ties; the `distinct?` option attempts to mitigate this.
;; *   The Power Divergence tests are typically applied to counts or binned data, while AD/KS can be applied directly to continuous data (or compared against continuous distributions).
;;
;; Handling Ties in KS Tests:
;;

;; The Kolmogorov-Smirnov (KS) test is theoretically defined for continuous distributions, where the probability of observing duplicate values (ties) is zero. In practice, real-world data often contains ties, which can affect the accuracy of the KS test, particularly for the exact p-value calculation methods. The `ks-test-one-sample` and `ks-test-two-samples` functions in `fastmath.stats` provide the `:distinct?` option to address this.
;; 
;; The `:distinct?` option controls how ties are handled:
;; 
;; *   `:ties` (default): This is the default behavior. It keeps all observed data points, including ties. If the `:method` is `:exact`, information about the ties is passed to the underlying exact calculation algorithm to attempt a correction. The accuracy depends on the specific algorithm's tie-handling capabilities.
;; *   `true`: This applies the Clojure `distinct` function *separately* to each input sequence (`xs` and `ys` for the two-sample test) before combining and processing. This removes duplicate values *within* each sample but does not guarantee that ties between the samples are resolved or correctly handled for the exact method.
;; *   `:jitter`: This adds a small amount of random noise to each data point to break all ties. This is a common pragmatic approach for dealing with ties in continuous data tests when an exact tie correction is unavailable or complex, but it slightly alters the original data.
;; *   `false`: The data is used exactly as provided, without any specific handling or correction for ties.
;; 
;; The presence and handling of ties can significantly influence the calculated p-value, especially when using exact calculation methods. Choosing the appropriate method depends on the nature of the data and the required accuracy.

;; **Examples:**
;; 
;; AD test comparing `setosa-sepal-length` data to a normal distribution.
;; 

(utls/zp
 (stats/ad-test-one-sample setosa-sepal-length
                           (r/distribution :normal {:mu (stats/mean setosa-sepal-length)
                                                    :sd (stats/stddev setosa-sepal-length)})))

;; KS test comparing `setosa-sepal-length` data to an estimated empirical distribution build from normal samples (KDE based).

(utls/zp
 (let [d (r/distribution :normal {:mu (stats/mean setosa-sepal-length)
                                  :sd (stats/stddev setosa-sepal-length)})]
   (stats/ks-test-one-sample setosa-sepal-length (r/->seq d 100) {:kernel :epanechnikov})))

;; Two-sample KS test comparing `setosa-sepal-length` and `virginica-sepal-length`.


(utls/zp
 (stats/ks-test-two-samples setosa-sepal-length virginica-sepal-length))

;; ### ANOVA and Rank Sum Tests

;; ;; Analysis of Variance (ANOVA) and rank sum tests are used to determine if there are statistically significant differences between the means (ANOVA) or distributions (rank sum tests) of two or more independent groups. One-way tests are used when comparing groups based on a single categorical independent variable.

;; ::: {.callout-tip title="Defined functions"}
;; * `one-way-anova-test`
;; * `kruskal-test`
;; :::

;; * **One-way ANOVA Test** (`one-way-anova-test`): A parametric test used to compare the means of two or more independent groups. It assesses whether the variation among group means is larger than the variation within groups.
;;     -   **Null Hypothesis ($H_0$)**: The means of all groups are equal.
;;     -   **Alternative Hypothesis ($H_1$)**: At least one group mean is different from the others.
;;     -   **Assumptions**: Independence of observations, normality within each group, and homogeneity of variances (equal variances across groups).
;;     -   **Test Statistic (F-statistic)**: Calculated as the ratio of the variance between groups (Mean Square Treatment, $MSt$) to the variance within groups (Mean Square Error, $MSe$).
;;         $$ F = \frac{MSt}{MSe} $$
;;         Where $MSt = SSt / DFt$ and $MSe = SSe / DFe$. $SSt$ is the Sum of Squares Treatment (between groups), $SSe$ is the Sum of Squares Error (within groups), $DFt$ is the degrees of freedom for the treatment (number of groups - 1), and $DFe$ is the degrees of freedom for the error (total observations - number of groups).
;;
;; * **Kruskal-Wallis H-Test** (`kruskal-test`): A non-parametric alternative to the one-way ANOVA. It is used to compare the distributions of two or more independent groups when the assumptions of ANOVA (especially normality) are not met. It tests whether the groups come from the same distribution.
;;     -   **Null Hypothesis ($H_0$)**: The distributions of all groups are identical (specifically, that the median ranks are equal).
;;     -   **Alternative Hypothesis ($H_1$)**: At least one group's distribution is different from the others.
;;     -   **Assumptions**: Independence of observations. While it doesn't assume normality, it assumes that the distributions have similar shapes if the intent is to compare medians.
;;     -   **Test Statistic (H)**: Calculated based on the ranks of the combined data from all groups. The sum of ranks is computed for each group, and H is derived from deviations of these rank sums from what would be expected under the null hypothesis.
;;
;; **Comparison:**
;;
;; -   **`one-way-anova-test`** is a parametric test that compares *means*. It assumes normality and equal variances.
;; -   **`kruskal-test`** is a non-parametric test that compares *distributions* (or median ranks). It does not assume normality or equal variances, making it suitable for ordinal data or data violating ANOVA assumptions.
;;
;; Both functions accept a sequence of sequences (`xss`) where each inner sequence represents a group. Both return a map containing the test statistic, degrees of freedom, and p-value.
;;
;; **Arguments:**
;;
;; *   `xss` (sequence of sequences of numbers): The groups of data to compare.
;; *   `params` (map, optional): An options map.
;;     *   `:sides` (keyword, default `:one-sided-greater` for both): Specifies the alternative hypothesis side for the test statistic's distribution.
;;
;; **Return Value (Map):**
;;
;; *   `:stat` (double): The calculated test statistic (F for ANOVA, H for Kruskal-Wallis).
;; *   `:p-value` (double): The probability of observing the test statistic or more extreme values under the null hypothesis.
;; *   `:df` (long or vector): Degrees of freedom. Vector `[DFt, DFe]` for ANOVA, scalar `DFt` for Kruskal-Wallis.
;; *   `:n` (sequence of longs): Sample sizes of each group.
;; *   `:k` (long): Number of groups (Kruskal-Wallis only).
;; *   `:sides` (keyword): The alternative hypothesis side used.
;; *   **ANOVA only**: `:F` (alias for `:stat`), `:SSt`, `:SSe`, `:DFt`, `:DFe`, `:MSt`, `:MSe`.
;; *   **Kruskal-Wallis only**: `:H` (alias for `:stat`).
;;
;; Let's compare the sepal lengths across the three iris species using both tests.
;;

(utls/zp (stats/one-way-anova-test (vals sepal-lengths)))

(utls/zp (stats/kruskal-test (vals sepal-lengths)))

;; ### Autocorrelation Tests

;; Autocorrelation tests examine whether values in a sequence are correlated with
;; past values in the same sequence. This is particularly important when analyzing
;; time series data or the residuals from a regression analysis, as autocorrelated
;; residuals violate the assumption of independent errors, which can invalidate
;; statistical inferences. `fastmath.stats` provides the Durbin-Watson test,
;; a common method for detecting first-order (lag-1) autocorrelation in regression residuals.

;; ::: {.callout-tip title="Defined functions"}
;; * `durbin-watson`
;; :::

;; *   **Durbin-Watson Test** (`durbin-watson`): Calculates the Durbin-Watson statistic (d),
;;     which is used to test for the presence of serial correlation, especially
;;     first-order (lag-1) autocorrelation, in the residuals of a regression analysis.
;;
;;     The Durbin-Watson statistic is calculated as:
;;     $$ d = \frac{\sum_{t=2}^T (e_t - e_{t-1})^2}{\sum_{t=1}^T e_t^2} $$
;;     where $e_t$ are the residuals at time $t$, and $T$ is the total number of observations.
;;
;;     The value of the statistic $d$ ranges from 0 to 4.
;;
;;     *   Values near 2 suggest no first-order autocorrelation.
;;     *   Values less than 2 suggest positive autocorrelation (residuals tend to be followed by residuals of the same sign).
;;     *   Values greater than 2 suggest negative autocorrelation (residuals tend to be followed by residuals of the opposite sign).
;;
;;     Testing for significance requires comparing the calculated statistic to lower ($d_L$)
;;     and upper ($d_U$) critical values from Durbin-Watson tables, which depend on the
;;     sample size and number of predictors in the regression model.
;;
;;     Parameters:
;;
;;     *   `rs` (sequence of numbers): The sequence of residuals from a regression model.
;;         The sequence should represent observations ordered by time or sequence index.
;;
;;     Returns the calculated Durbin-Watson statistic as a double.
;;
;;     Note: This function only calculates the statistic. Determining statistical
;;     significance typically requires consulting Durbin-Watson critical value tables
;;     based on the sample size and the number of independent variables in the model.
;;
;;     Let's calculate the Durbin-Watson statistic for a few example sequences representing different scenarios: no autocorrelation, positive autocorrelation, and negative autocorrelation.
;;

(def residuals-no-autocorrelation (repeatedly 100 r/grand))
(def residuals-positive-autocorrelation [2.0 1.8 1.5 1.0 -0.3 -0.8 -1.2 -1.5])
(def residuals-negative-autocorrelation [1.0 -1.0 0.5 -0.5 0.2 -0.2 0.1 -0.1])

(utls/examples-note
  (stats/durbin-watson residuals-no-autocorrelation)
  (stats/durbin-watson residuals-positive-autocorrelation)
  (stats/durbin-watson residuals-negative-autocorrelation))

;; ## Time Series Analysis

;; Functions specifically for analyzing sequential data, such as time series.
;; These functions help identify patterns like trends, seasonality, and autocorrelation,
;; which are crucial for understanding the underlying process generating the data
;; and for building forecasting models (like ARIMA).

;; ::: {.callout-tip title="Defined functions"}
;; * `acf`
;; * `pacf`
;; * `acf-ci`, `pacf-ci`
;; :::

;; Autocorrelation and Partial Autocorrelation functions are fundamental tools for
;; analyzing the dependence structure of a time series. They measure the correlation
;; between a series and its lagged values.
;;
;; *   **Autocorrelation Function (ACF)** (`acf`): Measures the linear dependence
;;     between a time series and its lagged values. Specifically, the ACF at lag $k$
;;     is the correlation between the series and itself shifted by $k$ time units.
;;     It captures both direct and indirect dependencies.
;;     For a stationary series, the sample ACF at lag $k$ is estimated as:
;;     $$ \rho_k = \frac{\sum_{t=k+1}^n (y_t - \bar{y})(y_{t-k} - \bar{y})}{\sum_{t=1}^n (y_t - \bar{y})^2} $$
;;     where $y_t$ is the observation at time $t$, $\bar{y}$ is the sample mean, and $n$ is the series length.
;;
;; *   **Partial Autocorrelation Function (PACF)** (`pacf`): Measures the linear
;;     dependence between a time series and its lagged values *after removing the
;;     linear dependence from the intermediate lags*. The PACF at lag $k$ is the
;;     correlation between $y_t$ and $y_{t-k}$, conditional on the intermediate
;;     observations $y_{t-1}, y_{t-2}, \dots, y_{t-k+1}$. It isolates the direct
;;     relationship at each lag.
;;     The PACF values ($\phi_{kk}$) are the last coefficients in a sequence of
;;     autoregressive models of increasing order. For example, $\phi_{11}$ is the
;;     correlation at lag 1, $\phi_{22}$ is the correlation at lag 2 after accounting
;;     for lag 1, etc.
;;
;; **Comparison:**
;;
;; *   The **ACF** shows correlations at various lags, including those that are
;;     simply due to preceding correlations. It tends to decay gradually for
;;     autoregressive (AR) processes and cut off sharply for moving average (MA) processes.
;; *   The **PACF** shows only the *direct* correlation at each lag, after removing
;;     the effects of shorter lags. It tends to cut off sharply for AR processes
;;     and decay gradually for MA processes.
;;
;; These patterns in ACF and PACF plots are diagnostic tools for identifying the
;; order of AR or MA components in time series models (e.g., ARIMA).
;;
;; **Confidence Intervals:**
;;
;; *   `acf-ci`, `pacf-ci`: Calculate ACF and PACF values, respectively, and
;;     provide approximate confidence intervals around these estimates. These
;;     intervals help determine whether the autocorrelation/partial autocorrelation
;;     at a given lag is statistically significant (i.e., unlikely to be zero
;;     in the population). For a stationary series, the standard error of the
;;     sample ACF at lag $k$ is approximately $1/\sqrt{n}$ for $k>0$. The PACF
;;     has a similar standard error for lags $k>0$. The confidence interval
;;     at level $\alpha$ is typically $\pm z_{\alpha/2} / \sqrt{n}$.
;;     `acf-ci` also provides cumulative confidence intervals (`:cis`) based on
;;     the variance of the sum of squared autocorrelations up to each lag.
;;
;; **Arguments:**
;;
;; *   `data` (seq of numbers): The time series data.
;; *   `lags` (long or seq of longs, optional): The maximum lag to calculate ACF/PACF for (if a number), or a specific list of lags (for `acf` only). Defaults to `(dec (count data))`.
;; *   `alpha` (double, optional, for `*-ci` functions): Significance level for the confidence intervals (default 0.05 for 95% CI).
;;
;; **Returns:**
;;
;; *   `acf`: A sequence of doubles representing the autocorrelation coefficients
;;     at lags 0 (always 1.0) up to `lags` (or specified lags).
;; *   `pacf`: A sequence of doubles representing the partial autocorrelation
;;     coefficients at lags 0 (always 0.0) up to `lags`.
;; *   `acf-ci`, `pacf-ci`: A map containing:
;;     *   `:ci` (double): The value of the standard confidence interval (e.g., $1.96/\sqrt{n}$ for 95% CI).
;;     *   `:acf` or `:pacf` (seq of doubles): The calculated ACF or PACF values.
;;     *   `:cis` (seq of doubles, only for `acf-ci`): The cumulative confidence intervals for ACF.

;; Let's illustrate ACF and PACF with the white noise data ($sigma=0.5$) and AR, MA and ARFIMA processes.

(def white-noise (take 1000 (r/white-noise (r/rng :mersenne 1) 0.5)))

(utls/examples-note
  (stats/acf white-noise 10)
  (stats/pacf white-noise 10)
  (stats/acf-ci white-noise 10)
  (stats/pacf-ci white-noise 10))

(gg/->image (gg/line (range 200) white-noise {:title "white noise"}) {:width 800 :height 200})

(kind/table
 [[(gg/->image (gg/acf white-noise "white noise"))
   (gg/->image (gg/pacf white-noise "white noise"))]])

;; And now for AR(30) time series.

(def ar20 (take 10000 (drop 100 (r/ar (v/sq (m/slice-range 0.35 0.001 20)) white-noise))))

(utls/examples-note
  (stats/acf ar20 10)
  (stats/pacf ar20 10)
  (stats/acf-ci ar20 10)
  (stats/pacf-ci ar20 10))

(gg/->image (gg/line (range 200) ar20 {:title "AR(20)"}) {:width 800 :height 200})

(kind/table
 [[(gg/->image (gg/acf ar20 "AR(20)"))
   (gg/->image (gg/pacf ar20 "AR(20)"))]])

;; MA(10) time series.

(def ma10 (take 1000 (drop 100 (r/ma [0.1 0.1 0.1 2 1 0.1 0.1 0.1 -1 -2]))))

(utls/examples-note
  (stats/acf ma10 10)
  (stats/pacf ma10 10)
  (stats/acf-ci ma10 10)
  (stats/pacf-ci ma10 10))

(gg/->image (gg/line (range 200) ma10 {:title "MA(10)"}) {:width 800 :height 200})

(kind/table
 [[(gg/->image (gg/acf ma10 "MA(10)"))
   (gg/->image (gg/pacf ma10 "MA(10)"))]])

;; ARFIMA(3,0.1,3)

(def arfima (take 1000 (drop 100 (r/arfima [0.1 0.1 -0.2] 0.1 [0.9 0.1 0.01]))))

(utls/examples-note
  (stats/acf arfima 10)
  (stats/pacf arfima 10)
  (stats/acf-ci arfima 10)
  (stats/pacf-ci arfima 10))

(gg/->image (gg/line (range 200) arfima {:title "ARFIMA"}) {:width 800 :height 200})

(kind/table
 [[(gg/->image (gg/acf arfima "ARFIMA"))
   (gg/->image (gg/pacf arfima "ARFIMA"))]])

;; ## Histograms

;; Histograms are fundamental graphical and statistical tools used to represent the distribution of numerical data. They group data into bins along the x-axis and display the frequency (count or proportion) of data points falling into each bin as bars on the y-axis. Histograms provide a visual summary of the shape, center, and spread of the data, highlighting features like modality, symmetry, and outliers.

;; ::: {.callout-tip title="Defined functions"}
;; * `histogram`
;; * `estimate-bins`
;; :::

;; 
;; `fastmath.stats` provides functions to construct histograms and assist in choosing appropriate binning strategies:
;; 
;; *   `histogram`: Computes and returns the data structure representing a histogram. It takes the input data and parameters defining the bins (number, estimation method, or explicit intervals). It can also process collections of data sequences (for grouped histograms).
;; *   `estimate-bins`: A utility function to recommend the number of bins for a given dataset based on various commonly used heuristic rules.
;; 
;; The `histogram` function calculates the counts of data points falling into predefined intervals (bins).
;; 
;; **Parameters:**
;; 
;; *   `vs` (sequence of numbers or sequence of sequences): The input data. Can be a single sequence for a simple histogram or a collection of sequences for grouped histograms.
;; *   `bins-or-estimate-method` (number, keyword, or sequence, optional): Defines the histogram bins.
;;     *   A number: The desired number of bins. The function calculates equally spaced intervals between the minimum and maximum values.
;;     *   A keyword (default: `:freedman-diaconis`): Uses a specific heuristic to estimate the number of bins (see [[estimate-bins]]).
;;     *   A sequence of numbers: Explicitly defines the bin edges (intervals). Data points are counted if they fall into `[edge_i, edge_{i+1})`.
;;     *   If omitted, uses the default `:freedman-diaconis` estimation method.
;; *   `mn`, `mx` (doubles, optional): Explicit minimum and maximum values to consider for binning. Data outside `[mn, mx]` are excluded. If omitted, the minimum and maximum of the data are used. When explicit `bins-or-estimate-method` (sequence of edges) is provided, `mn` and `mx` are inferred from the provided edges and `vs` data is filtered to fit.
;; 
;; **Return Value (Map):**
;; 
;; Returns a map describing the histogram structure and counts. Key elements include:
;; 
;; *   `:size`: The number of bins.
;; *   `:step`: The average width of the bins (only applicable for equally spaced bins).
;; *   `:samples`: The total number of data points included in the histogram (may be less than input if `mn`/`mx` are specified).
;; *   `:min`, `:max`: The minimum and maximum data values used for binning.
;; *   `:intervals`: The sequence of numbers defining the bin edges.
;; *   `:bins`: A sequence of `[lower-edge, count]` pairs for each bin.
;; *   `:frequencies`: A map where keys are the average value of each bin and values are the counts.
;; *   `:bins-maps`: A sequence of detailed maps for each bin, including `:min`, `:max`, `:step` (bin width), `:count`, `:avg` (mean value within the bin), and `:probability` (count / total samples).
;; 
;; If the input `vs` is a sequence of sequences, the function returns a sequence of such maps, one for each inner sequence.

(utls/zp (stats/histogram mpg))

(gg/->image (gg/from-histogram mpg {:title "MPG histogram (bins = :freedman-diaconis, default)"}))

;; Histogram with 3 bins

(utls/zp (stats/histogram mpg 3))

(gg/->image (gg/from-histogram mpg {:title "MPG histogram (bins = 3)"
                                    :bins 3}))

;; Histogram with irregular bins

(utls/zp (stats/histogram mpg [10 20 25 27 35]))

(gg/->image (gg/from-histogram mpg {:title "MPG histogram (custom bins)"
                                    :bins [10 20 25 27 35]}))

;; Histogram with given min and max range

(gg/->image (gg/from-histogram mpg {:title "MPG histogram (bins = 30, range = [0,40])"
                                    :bins 30
                                    :span [0 50]}))

;; The `estimate-bins` function provides recommendations for the number of bins in a histogram based on common rules.
;; 
;; **Parameters:**
;; 
;; *   `vs` (sequence of numbers): The input data.
;; *   `bins-or-estimate-method` (keyword, optional): The rule to use for estimation.
;;     *   A keyword (default: `:freedman-diaconis`): Specifies the rule.
;;         *   `:sqrt`: Square root rule (simple). $k = \lceil\sqrt{n}\rceil$
;;         *   `:sturges`: Sturges' rule (assumes approximately normal data). $k = \lceil\log_2(n) + 1\rceil$
;;         *   `:rice`: Rice rule. $k = \lceil 2 n^{1/3} \rceil$
;;         *   `:doane`: Doane's rule (modification of Sturges' for non-normal data). $k = \lceil \log_2(n) + 1 + \log_2(1 + \frac{|g_1|}{\sigma_{g_1}}) \rceil$, where $g_1$ is sample skewness and $\sigma_{g_1}$ is its standard error.
;;         *   `:scott`: Scott's normal reference rule (bin width $h = 3.5 \hat{\sigma} n^{-1/3}$). $k = \lceil (max - min)/h \rceil$
;;         *   `:freedman-diaconis` (default): Freedman-Diaconis rule (robust to outliers, bin width $h = 2 \cdot IQR \cdot n^{-1/3}$). $k = \lceil (max - min)/h \rceil$
;;     *   A number: Returns the number itself (useful for passing a fixed value through).
;;     *   If omitted, uses the default `:freedman-diaconis` rule.
;; 
;; **Return Value (Long):**
;; 
;; Returns the estimated number of bins as a long integer. The returned value is constrained to be no greater than the number of samples.

;; ---

;; Estimate bins for `alcohol` data using different methods

(utls/examples-note
  (stats/estimate-bins alcohol)
  (stats/estimate-bins alcohol :sqrt)
  (stats/estimate-bins alcohol :sturges)
  (stats/estimate-bins alcohol :rice)
  (stats/estimate-bins alcohol :doane)
  (stats/estimate-bins alcohol :scott)
  (stats/estimate-bins alcohol 3))

(kind/table
 [[(gg/->image (gg/from-histogram alcohol {:title "Alcohol histogram (bins = :freedman-diaconis, default)"}))
   (gg/->image (gg/from-histogram alcohol {:title "Alcohol histogram (bins = :sqrt)"
                                           :bins :sqrt}))]
  [(gg/->image (gg/from-histogram alcohol {:title "Alcohol histogram (bins = :sturges)"
                                           :bins :sturges}))
   (gg/->image (gg/from-histogram alcohol {:title "Alcohol histogram (bins = :rice)"
                                           :bins :rice}))]
  [(gg/->image (gg/from-histogram alcohol {:title "Alcohol histogram (bins = :doane)"
                                           :bins :doane}))
   (gg/->image (gg/from-histogram alcohol {:title "Alcohol histogram (bins = :scott)"
                                           :bins :scott}))]])

;; ## Bootstrap

;; Bootstrap is a widely used resampling technique in statistics. It involves repeatedly drawing samples with replacement from an original dataset (or simulating data from a model) to create many "bootstrap samples". By analyzing the distribution of a statistic computed from each of these bootstrap samples, one can estimate the sampling distribution of the statistic, its standard error, bias, and construct confidence intervals without relying on strong parametric assumptions.

;; ::: {.callout-tip title="Defined functions"}
;; * `bootstrap`
;; * `jackknife`, `jackknife+`
;; * `bootstrap-stats`
;; * `ci-normal`, `ci-basic`, `ci-percentile`, `ci-bc`, `ci-bca`, `ci-studentized`, `ci-t`
;; :::

;; The core of the bootstrap process is generating the resamples. `fastmath.stats` provides the main `bootstrap` function for general resampling, and dedicated functions for specific resampling techniques like Jackknife.
;;
;; *   `bootstrap`: Generates bootstrap samples from data or a probabilistic model. It supports various sampling methods (standard resampling with replacement, Jackknife variants, or sampling from a distribution) and options like smoothing, antithetic sampling, and handling multidimensional data.
;; *   `jackknife`: Generates samples using the leave-one-out Jackknife method. For a dataset of size $n$, it produces $n$ samples, each by removing one observation from the original dataset. This is a less computationally intensive resampling method than standard bootstrap, primarily used for bias and variance estimation.
;; *   `jackknife+`: Generates samples using the positive Jackknife method. For a dataset of size $n$, it produces $n$ samples, each by adding one extra copy of an observation to the original dataset, resulting in samples of size $n+1$.

;; The primary function for generating bootstrap samples. It's flexible, supporting both nonparametric (resampling from data) and parametric (sampling from a model) approaches, and various options.
;;
;; **Parameters:**
;;
;; *   `input` (sequence or map): The data source.
;;     *   If a sequence of data values (e.g., `[1.2 3.4 5.0]`), it's treated as nonparametric input.
;;     *   If a map `{:data data-sequence :model model-object (optional)}`, it supports parametric bootstrap. If `:model` is omitted, a distribution is automatically built from `:data`.
;;     *   Can be a sequence of sequences (e.g., `[[1 2] [3 4]]`) for multidimensional data when `:dimensions` is `:multi`.
;; *   `statistic` (function, optional): A function `(fn [sample-sequence])` that calculates a statistic (e.g., `fastmath.stats/mean`). If provided, `bootstrap-stats` is automatically called on the results. If `nil`, the raw bootstrap samples are returned.
;; *   `params` (map, optional): Configuration options:
;;     *   `:samples` (long, default: 500): Number of bootstrap samples. Ignored for `:jackknife`/`:jackknife+`.
;;     *   `:size` (long, optional): Size of each sample. Defaults to original data size. Ignored for `:jackknife`/`:jackknife+`.
;;     *   `:method` (keyword, optional): Sampling method.
;;         *   `nil` (default): Standard resampling with replacement from data or model.
;;         *   `:jackknife`: Uses [[jackknife]].
;;         *   `:jackknife+`: Uses [[jackknife+]].
;;         *   Other keywords are passed to `fastmath.random/->seq` for sampling from a distribution model.
;;     *   `:rng` (random number generator, optional): `fastmath.random` RNG. Defaults to JVM RNG.
;;     *   `:smoothing` (keyword, optional): Applies smoothing. `:kde` for Kernel Density Estimation on the model; `:gaussian` adds noise to resampled values.
;;     *   `:distribution` (keyword, default `:real-discrete-distribution`): Type of distribution to build if `:model` is missing.
;;     *   `:dimensions` (keyword, optional): `:multi` for multidimensional data (sequence of sequences).
;;     *   `:antithetic?` (boolean, default `false`): Uses antithetic sampling (requires distribution model).
;;     *   `:include?` (boolean, default `false`): Includes original dataset as one sample.
;;
;; **Returns:**
;;
;; *   If `statistic` is provided: A map including original input + analysis results from `bootstrap-stats` (e.g., `:t0`, `:ts`, `:bias`, `:mean`, `:stddev`).
;; *   If `statistic` is `nil`: A map including original input + `:samples` (a collection of bootstrap sample sequences).
;;

;; Let's demonstrate various sampling methods for `mpg` data and see the outcome.

;; First we'll generate two samples using default (resampling with replacement) method without calculating statistics.

(def rng (r/rng :mersenne 12345))

(def boot-mpg (boot/bootstrap mpg {:samples 2 :rng rng}))

boot-mpg

(kind/table
 [[(gg/->image (gg/density (first (:samples boot-mpg))
                           {:title "Density of the first sample"}))
   (gg/->image (gg/density (second (:samples boot-mpg))
                           {:title "Density of the second sample"}))]])

;; Now we'll create set of samples using various bootstrapping methods: jackknife, jackknife+, gaussian and kde (epanechnikov kernel) smoothing, antithetic, systematic, stratified.

;; **jackknife**

;; Removes one selected data point from each sample.

(def boot-mpg-jacknife (boot/bootstrap mpg {:method :jackknife}))

(count (:samples boot-mpg-jacknife))

(kind/table
 [[(gg/->image (gg/density (nth (:samples boot-mpg-jacknife) 10)
                           {:title "Density of the 10th sample for jackknife method"}))
   (gg/->image (gg/density (nth (:samples boot-mpg-jacknife) 30)
                           {:title "Density of the 30th sample for jackknife method"}))]])

;; **jackknife+**

;; Duplicates one selected data point to each sample.

(def boot-mpg-jacknife+ (boot/bootstrap mpg {:method :jackknife+}))

(count (:samples boot-mpg-jacknife))

(kind/table
 [[(gg/->image (gg/density (nth (:samples boot-mpg-jacknife+) 10)
                           {:title "Density of the 10th sample for jackknife+ method"}))
   (gg/->image (gg/density (nth (:samples boot-mpg-jacknife+) 30)
                           {:title "Density of the 30th sample for jackknife+ method"}))]])

;; **gaussian smoothing**

;; Adds random gaussian noise to bootstrapped samples.

(def boot-mpg-gaussian (boot/bootstrap mpg {:samples 2 :rng rng :smoothing :gaussian}))

boot-mpg-gaussian

(kind/table
 [[(gg/->image (gg/density (first (:samples boot-mpg-gaussian))
                           {:title "Density of the first sample, gaussian smoothing"}))
   (gg/->image (gg/density (second (:samples boot-mpg-gaussian))
                           {:title "Density of the second sample, gaussian smoothing"}))]])


;; **KDE (epanechnikov) smoothing**

;; Builds KDE based continuous distributions and samples from it.

(def boot-mpg-kde (boot/bootstrap mpg {:samples 2 :rng rng :smoothing :kde
                                     :kernel :epanechnikov}))

boot-mpg-kde

(kind/table
 [[(gg/->image (gg/density (first (:samples boot-mpg-kde))
                           {:title "Density of the first sample, kde smoothing"}))
   (gg/->image (gg/density (second (:samples boot-mpg-kde))
                           {:title "Density of the second sample, kde smoothing"}))]])

;; **antithetic**

;; Samples in pairs.

(def boot-mpg-antithetic (boot/bootstrap mpg {:samples 2 :rng rng :antithetic? true}))

boot-mpg-antithetic

(kind/table
 [[(gg/->image (gg/density (first (:samples boot-mpg-antithetic))
                           {:title "Density of the first sample, antithetic"}))
   (gg/->image (gg/density (second (:samples boot-mpg-antithetic))
                           {:title "Density of the second sample, antithetic"}))]])

;; **systematic**

;; When `:systematic` method is used, the result will be visible when model is continuous.

(def boot-mpg-systematic (boot/bootstrap mpg {:samples 2 :rng rng :method :systematic
                                            :smoothing :kde}))

boot-mpg-systematic

(kind/table
 [[(gg/->image (gg/density (first (:samples boot-mpg-systematic))
                           {:title "Density of the first sample, systematic"}))
   (gg/->image (gg/density (second (:samples boot-mpg-systematic))
                           {:title "Density of the second sample, systematic"}))]])

;; **stratified**

;; Similarly to a `:systematic`, we should sample from conitinuous distribution.

(def boot-mpg-stratified (boot/bootstrap mpg {:samples 2 :rng rng :method :stratified
                                            :smoothing :kde}))

boot-mpg-stratified

(kind/table
 [[(gg/->image (gg/density (first (:samples boot-mpg-stratified))
                           {:title "Density of the first sample, stratified"}))
   (gg/->image (gg/density (second (:samples boot-mpg-systematic))
                           {:title "Density of the second sample, stratified"}))]])

;; **Model**

;; Automatically created model is in most cases a distribution object and it can be used to generate more data.

(r/->seq (:model boot-mpg) 10)

(r/->seq (:model boot-mpg) 10 :stratified)

;; Let's extract median

(r/icdf (:model boot-mpg) 0.5)

;; **Model based sampling**

;; Now let's simulate our data when model is given. Let's assume that `mpg` follows normal distribution. So instead of resampling we will draw samples from fitted normal distribution. 

(def boot-mpg-model (boot/bootstrap {:data mpg
                                   :model (r/distribution :normal {:mu (stats/mean mpg)
                                                                   :sd (stats/stddev mpg)})}
                                  stats/mean
                                  {:samples 10000 :rng rng}))

(kind/table
 [[(gg/->image (gg/density (first (:samples boot-mpg-model))
                           {:title "Density of the first sample drawn from normal distribution."}))
   (gg/->image (gg/density (:ts boot-mpg-model)
                           {:vline {:xintercept (:t0 boot-mpg-model) :color "#0000FF"}
                            :title (str "Bootstraped mean, bias="
                                        (m/approx (:bias boot-mpg-model) 6)
                                        ", method=default with provided model")}))]])

;; **Non-numerical data**

;; Bootstrap can be also called on categorical data.

(utls/zp (boot/bootstrap [:cat :dog :cat :dog :dog] {:samples 2 :size 20}))

;; ### Statistics

;; Once bootstrap samples are generated, one typically applies the statistic of interest to each sample to get the bootstrap distribution of the statistic.
;;
;; *   `bootstrap-stats`: Computes a specified statistic on the original data (`t0`) and on each bootstrap sample (`ts`), and calculates descriptive statistics (mean, median, variance, stddev, SEM, bias) for the distribution of `ts`. This function is called automatically by `bootstrap` when a statistic function is provided.
;;
;; **Parameters:**
;;
;; *   `boot-data` (map): A map containing `:data` (original data) and `:samples` (collection of bootstrap samples).
;; *   `statistic` (function): The function `(fn [sample-sequence])` to apply to the data and samples.
;;
;; **Returns:** The input `boot-data` map augmented with results:

;; * `:statistic` - statistic calculation function
;; * `:t0` - statistic for original data
;; * `:ts` - statistics for all bootstrapped samples
;; * `:bias` - difference between `t0` and mean of `ts`
;; * `:mean` - mean of `ts`
;; * `:median` - median of `ts`
;; * `:variance` - variance of `ts`
;; * `:stddev` - standard deviation of `ts`
;; * `:sem` - standard error of mean of `ts`
;;
;; Example (usually called internally by `bootstrap`):
;;
(def raw-boot-samples (boot/bootstrap mpg {:samples 1000}))
(def analyzed-boot-data (boot/bootstrap-stats raw-boot-samples stats/mean))
(utls/zp (dissoc analyzed-boot-data :samples :ts))

;; Now we'll see bootstrapped `ts` for different types of sampling methods and sample size = 10000. Vertical line marks `t0`.

(kind/table
 (let [t0 (stats/mean mpg)]
   [[(let [b (boot/bootstrap mpg stats/mean {:samples 10000 :rng rng})]
       (gg/->image (gg/density (:ts b)
                               {:vline {:xintercept t0 :color "#0000FF"}
                                :title (str "Bootstraped mean, bias="
                                            (m/approx (:bias b) 6)
                                            ", method=default")})))
     (let [b (boot/bootstrap mpg stats/mean {:samples 10000 :rng rng :antithetic? true})]
       (gg/->image (gg/density (:ts b)
                               {:vline {:xintercept t0 :color "#0000FF"}
                                :title (str "Bootstraped mean, bias="
                                            (m/approx (:bias b) 6)
                                            ", method=default, antithetic")})))]
    [(let [b (boot/bootstrap mpg stats/mean {:samples 10000 :rng rng :smoothing :gaussian})]
       (gg/->image (gg/density (:ts b)
                               {:vline {:xintercept t0 :color "#0000FF"}
                                :title (str "Bootstraped mean, bias="
                                            (m/approx (:bias b) 6)
                                            ", method=smoothing: gaussian")})))
     (let [b (boot/bootstrap mpg stats/mean {:samples 10000 :rng rng :smoothing :kde :kernel :epanechnikov})]
       (gg/->image (gg/density (:ts b)
                               {:vline {:xintercept t0 :color "#0000FF"}
                                :title (str "Bootstraped mean, bias="
                                            (m/approx (:bias b) 6)
                                            ", method=smoothing: KDE, epanechnikov")})))]
    [(let [b (boot/bootstrap mpg stats/mean {:samples 10000 :rng rng :method :systematic :smoothing :kde})]
       (gg/->image (gg/density (:ts b)
                               {:vline {:xintercept t0 :color "#0000FF"}
                                :title (str "Bootstraped mean, bias="
                                            (m/approx (:bias b) 6)
                                            ", method=systematic, smoothing: KDE, gaussian")})))
     (let [b (boot/bootstrap mpg stats/mean {:samples 10000 :rng rng :method :stratified :smoothing :kde})]
       (gg/->image (gg/density (:ts b)
                               {:vline {:xintercept t0 :color "#0000FF"}
                                :title (str "Bootstraped mean, bias="
                                            (m/approx (:bias b) 6)
                                            ", method=stratified, smoothing: KDE, gaussian")})))]
    [(let [b (boot/bootstrap mpg stats/mean {:rng rng :method :jackknife})]
       (gg/->image (gg/density (:ts b)
                               {:vline {:xintercept t0 :color "#0000FF"}
                                :title (str "MPG mean, bias="
                                            (m/approx (:bias b) 6)
                                            ", method=jacknife")})))
     (let [b (boot/bootstrap mpg stats/mean {:rng rng :method :jackknife+})]
       (gg/->image (gg/density (:ts b)
                               {:vline {:xintercept t0 :color "#0000FF"}
                                :title (str "MPG mean, bias="
                                            (m/approx (:bias b) 6)
                                            ", method=jacknife+")})))]]))


;; **Different statistic**

;; We can also estimate median instead of mean from smoothed samples.

(def boot-mpg-median (boot/bootstrap mpg stats/median {:samples 10000 :smoothing :gaussian :rng rng}))

(gg/->image (gg/density (:ts boot-mpg-median)
                        {:vline {:xintercept (:t0 boot-mpg-median) :color "#0000FF"}
                         :title (str "Bootstraped median, bias="
                                     (m/approx (:bias boot-mpg-median) 6)
                                     ", method=smoothing, gaussian")}))

;; ### CI

;; Bootstrap confidence intervals are constructed from the distribution of the statistic calculated on the bootstrap samples (`:ts`). `fastmath.stats.bootstrap` provides several methods, varying in complexity and robustness. All CI functions take the results map (typically from `bootstrap-stats`) and an optional significance level `alpha` (default 0.05 for 95% CI). They return a vector `[lower-bound, upper-bound, t0]`.
;;
;; *   `ci-normal`: Normal approximation interval. Assumes the distribution of `ts` is normal.
;; *   `ci-basic`: Basic (Percentile-t) interval.
;; *   `ci-percentile`: Percentile interval. The simplest method, directly using the quantiles of the bootstrap distribution of the statistic.
;; *   `ci-bc`: Bias-Corrected (BC) interval. Adjusts the Percentile interval by correcting for median bias in the bootstrap distribution.
;; *   `ci-bca`: Bias-Corrected and Accelerated (BCa) interval. Corrects for both bias ($z_0$) and skewness/non-constant variance (acceleration factor $a$). This is often considered the most accurate method but is more complex. Requires original `:data` and `:statistic` in the input map to calculate $a$ via Jackknife, or it estimates $a$ from `ts` skewness if data/statistic are missing.
;; *   `ci-studentized`: Studentized (Bootstrap-t) interval.
;; *   `ci-t`: T-distribution based interval.
;;
;; **Parameters (common to most CI functions):**
;;
;; *   `boot-data` (map): Bootstrap analysis results, typically from [[bootstrap-stats]]. Must contain `:t0` and `:ts`. Some methods (`ci-bca`, `ci-studentized`) require additional keys like `:data`, `:samples`, `:statistic`.
;; *   `alpha` (double, optional): Significance level. Default is 0.05 (for 95% CI).
;; *   `estimation-strategy` (keyword, optional, for `ci-basic`, `ci-percentile`, `ci-bc`, `ci-bca`): Quantile estimation strategy for `:s`. Defaults to `:legacy`. See [[quantiles]].
;;
;; **Returns (common):** A vector `[lower-bound, upper-bound, t0]`.
;;
;; Let's calculate various confidence intervals for the mean of `mpg` using the `boot-results` map obtained earlier (which includes `:t0`, `:ts`, etc.):

(def boot-results (boot/bootstrap mpg stats/mean {:samples 500 :rng rng}))

(utls/examples-note
  (:t0 boot-results)
  (boot/ci-normal boot-results)
  (boot/ci-basic boot-results)
  (boot/ci-percentile boot-results)
  (boot/ci-bc boot-results)
  (boot/ci-bca boot-results)
  (boot/ci-studentized boot-results)
  (boot/ci-t boot-results))

(kind/table
 [[(gg/->image (gg/density (:ts boot-results)
                           {:vline {:xintercept (:t0 boot-results) :color "#0000FF"}
                            :title "Density of bootstrapped mean (ts)"}))
   (gg/->image (gg/errorbars [[:normal (boot/ci-normal boot-results)]
                              [:basic (boot/ci-basic boot-results)]
                              [:percentile (boot/ci-percentile boot-results)]
                              [:bc (boot/ci-bc boot-results)]
                              [:bca (boot/ci-bca boot-results)]
                              [:studentized (boot/ci-studentized boot-results)]
                              [:t (boot/ci-t boot-results)]]
                             {:ylab "Bootstrap CI"
                              :title "Confidence intervals for MPG boostrap mean"
                              :vline {:xintercept (:t0 boot-results)}}))]])

;; ## Reference


(codox/make-public-fns-table-clay 'fastmath.stats)
(codox/make-public-fns-table-clay 'fastmath.stats.bootstrap)
