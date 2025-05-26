(ns stats
  (:require [fastmath.stats :as stats]
            [fastmath.stats.bootstrap :as boot]
            [fastmath.dev.ggplot :as gg]
            [fastmath.dev.codox :as codox]
            [fastmath.dev.clay :as utls]
            [fastmath.dev.dataset :as ds]
            [fastmath.core :as m]
            [scicloj.kindly.v4.kind :as kind]
            [fastmath.random :as r]))

;; # Statistics {.unnumbered}

;; Statistical functions

;; ::: {.callout-note}
;; This notebook is written with the support of Gemini LLM models:
;; 
;; * `gemini-2.5-pro-exp-03-25`
;; * `gemini-2.5-flash-preview-04-17`
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
;; **Arithmetic Mean** (`mean`): The sum of values divided by their count. It's the most common type of average.

;; $$\mu = \frac{1}{n} \sum_{i=1}^{n} x_i$$
;; 
;; **Geometric Mean** (`geomean`): The n-th root of the product of n numbers. Suitable for averaging ratios, growth rates, or values that are multiplicative in nature. Requires all values to be positive. It's less affected by extreme large values than the arithmetic mean, but more affected by extreme small values.

;; $$G = \left(\prod_{i=1}^{n} x_i\right)^{1/n} = \exp\left(\frac{1}{n}\sum_{i=1}^{n} \ln(x_i)\right)$$

;; **Harmonic Mean** (`harmean`): The reciprocal of the arithmetic mean of the reciprocals of the observations. Appropriate for averaging rates (e.g., speeds). It is sensitive to small values and requires all values to be positive and non-zero.

;; $$H = \frac{n}{\sum_{i=1}^{n} \frac{1}{x_i}}$$

;; **Power Mean** (`powmean`): Also known as the generalized mean or Hölder mean. It generalizes the arithmetic, geometric, and harmonic means. Defined by a power parameter $p$.

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

;; **Arithmetic Mean** (`mean`):

;; $$\mu_w = \frac{\sum_{i=1}^{n} w_i x_i}{\sum_{i=1}^{n} w_i}$$

;; **Geometric Mean** (`geomean`):

;; $$G_w = \left(\prod_{i=1}^{n} x_i^{w_i}\right)^{1/\sum w_i} = \exp\left(\frac{\sum_{i=1}^{n} w_i \ln(x_i)}{\sum_{i=1}^{n} w_i}\right)$$

;; **Harmonic Mean** (`harmean`):

;; $$H_w = \frac{\sum_{i=1}^{n} w_i}{\sum_{i=1}^{n} \frac{w_i}{x_i}}$$

;; **Power Mean** (`powmean`):

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

(stats/stats-map residual-sugar)

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

;; **Variance** quantifies the average squared difference of each data point from the mean. A higher variance indicates that the data points are more spread out, while a lower variance indicates they are clustered closer to the mean.

;; **Standard Deviation** is the square root of the variance. It is expressed in the same units as the data, making it more interpretable than variance as a measure of spread.


;; **Sample Variance** (`variance`) and **Sample Standard Deviation** (`stddev`): These are estimates of the population variance and standard deviation, calculated from a sample of data. They use a denominator of $N-1$ (Bessel's correction) to provide an unbiased estimate of the population variance.

;; $$s^2 = \frac{\sum_{i=1}^{n} (x_i - \bar{x})^2}{n-1}$$

;; $$s = \sqrt{s^2}$$

;; Both functions can optionally accept a pre-calculated mean (`mu`) as a second argument.

;; **Population Variance** (`population-variance`) and **Population Standard Deviation** (`population-stddev`): These are used when the data represents the entire population of interest, or when a biased estimate (maximum likelihood estimate) from a sample is desired. They use a denominator of $N$.

;; $$\sigma^2 = \frac{\sum_{i=1}^{N} (x_i - \mu)^2}{N}$$

;; $$\sigma = \sqrt{\sigma^2}$$

;; These also accept an optional pre-calculated population mean (`mu`).

;; **Weighted Variance** (`wvariance`, `population-wvariance`) and **Weighted Standard Deviation** (`wstddev`, `population-wstddev`): These calculate variance and standard deviation when each data point has an associated weight.
;;     For weighted sample variance (unbiased form, where $w_i$ are weights):

;; $$\bar{x}_w = \frac{\sum w_i x_i}{\sum w_i}$$

;; $$s_w^2 = \frac{\sum w_i (x_i - \bar{x}_w)^2}{(\sum w_i) - 1}$$

;; For weighted population variance:

;; $$\sigma_w^2 = \frac{\sum w_i (x_i - \bar{x}_w)^2}{\sum w_i}$$

;; Weighted standard deviations are the square roots of their respective variances.

;; **Pooled Variance** (`pooled-variance`) and **Pooled Standard Deviation** (`pooled-stddev`): These are used to estimate a common variance when data comes from several groups that are assumed to have the same population variance. Following methods can be used (where $k$ is the number of groups, each with $n_i$ number of observations and sample variance $s_i^2$)

;; * `:unbiased` (default)

;; $$s_p^2 = \frac{\sum_{i=1}^{k} (n_i-1)s_i^2}{\sum_{i=1}^{k} n_i - k}$$

;; * `:biased`

;; $$s_p^2 = \frac{\sum_{i=1}^{k} (n_i-1)s_i^2}{\sum_{i=1}^{k} n_i}$$

;; * `:avg` - simple average of group variances.

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
;; **Coefficient of Variation** (`variation`): This is a standardized measure of dispersion, calculated as the ratio of the standard deviation $s$ to the mean $\bar{x}$.

;; $$CV = \frac{s}{\bar{x}}$$

;; The CV is unitless, making it useful for comparing the variability of datasets with different means or units. It's most meaningful for data measured on a ratio scale (i.e., with a true zero point) and where all values are positive.
;;
;; **Standard Error of the Mean** (`sem`): The SEM estimates the standard deviation of the sample mean if you were to draw multiple samples from the same population. It indicates how precisely the sample mean estimates the true population mean.

;; $$SEM = \frac{s}{\sqrt{n}}$$

;; where $s$ is the sample standard deviation and $n$ is the sample size. A smaller SEM suggests a more precise estimate of the population mean.

;; **L-variation** (`l-variation`): Calculates the coefficient of L-variation. This is a dimensionless measure of dispersion, analogous to the coefficient of variation.

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
;; **Median Absolute Deviation** (`median-absolute-deviation` or `mad`): This is a robust measure of the variability of a univariate sample. It is defined as the median of the absolute deviations from the data's median.
;; 
;; $$MAD = \text{median}(|X_i - \text{median}(X)|)$$
;; 
;; If a specific center $c$ is provided, it's $MAD_c = \text{median}(|X_i - c|)$. Also, different estimation strategies can be used, see [[median]]
;; MAD is less sensitive to outliers than the standard deviation.
;; 
;; **Mean Absolute Deviation** (`mean-absolute-deviation`): This measures variability as the average of the absolute deviations from a central point, typically the data's mean.
;; 
;; $$MeanAD = \frac{1}{n} \sum_{i=1}^{n} |X_i - \text{mean}(X)|$$
;; 
;; If a specific center $c$ is provided, it's $MeanAD_c = \frac{1}{n} \sum_{i=1}^{n} |X_i - c|$.
;; MeanAD is more sensitive to outliers than MAD but less sensitive than the standard deviation.
;; 
;; **Pooled MAD** (`pooled-mad`): This function calculates a pooled estimate of the Median Absolute Deviation when data comes from several groups. For each group $i$, absolute deviations from its median $M_i$ are calculated: $Y_{ij} = |X_{ij} - M_i|$. The pooled MAD is then the median of all such $Y_{ij}$ values, scaled by a constant `const` (which defaults to approximately 1.4826, to make it comparable to the standard deviation for normal data).
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

;; **k-th Central Moment**: $\mu_k = E[(X - \mu)^k] \approx \frac{1}{n} \sum_{i=1}^{n} (x_i - \bar{x})^k$. Calculated when `center` is `nil` (default) and `:mean?` is `true` (default).

;; **k-th Raw Moment** (about origin): $\mu'_k = E[X^k] \approx \frac{1}{n} \sum_{i=1}^{n} x_i^k$. Calculated if `center` is `0.0`.

;; **k-th Absolute Central Moment**: $E[|X - \mu|^k] \approx \frac{1}{n} \sum_{i=1}^{n} |x_i - \bar{x}|^k$. Calculated if `:absolute?` is `true`.

;; **Normalization**: If `:normalize?` is `true`, the moment is divided by $\sigma^k$ (where $\sigma$ is the standard deviation), yielding a scale-invariant measure. For example, the 3rd normalized central moment is related to skewness, and the 4th to kurtosis.

;; **Power of sum of differences**:  If `:mean?` is `false`, the function returns the sum $\sum (x_i - c)^k$ (or sum of absolute values) instead of the mean.

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

;; * `:crow` (Crow-Siddiqui): Based on quantiles $Q_\alpha, Q_{1-\alpha}, Q_\beta, Q_{1-\beta}$ and centered by subtracting $2.906$. By default $\alpha=0.025$ and $\beta=0.25$.

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

;; **Basic Range:** Functions like `span` ($max - min$) and `extent` (providing $[min, max]$ and optionally the mean) offer simple measures of the total spread of the data.
;;
;; **Interquartile Range:** `iqr` ($Q_3 - Q_1$) specifically measures the spread of the middle 50% of the data, providing a robust alternative to the total range.
;;
;; **Symmetric Spread Intervals:** Functions ending in `-extent` such as `stddev-extent`, `mad-extent`, and `sem-extent` define intervals typically centered around the mean or median. They represent a range defined by adding/subtracting a multiple (usually 1) of a measure of dispersion (Standard Deviation, Median Absolute Deviation, or Standard Error of the Mean) from the central point.
;;
;; **Quantile-Based Intervals:** `percentile-extent`, `quantile-extent`, `pi`, `pi-extent`, and `hpdi-extent` define intervals based on quantiles or percentiles of the data. These functions capture specific ranges containing a certain percentage of the data points (e.g., the middle 95% defined by quantiles 0.025 and 0.975). `hpdi-extent` calculates the shortest interval containing a given proportion of data, based on empirical density.
;;
;; **Box Plot Boundaries:** `adjacent-values` (LAV, UAV) and fence functions (`inner-fence-extent`, `outer-fence-extent`) calculate specific bounds based on quartiles and multiples of the IQR. These are primarily used in box plot visualization and as a conventional method for identifying potential outliers.
;;
;; **Confidence and Prediction Intervals:** `ci`, `percentile-bc-extent`, and `percentile-bca-extent` provide inferential intervals. `ci` estimates a confidence interval for the population mean using the t-distribution. `percentile-bc-extent` and `percentile-bca-extent` (Bias-Corrected and Bias-Corrected Accelerated) are advanced bootstrap methods for estimating confidence intervals for statistics, offering robustness against non-normality and bias.

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
;; * `mcc` ;; Matthews correlation coefficient / Phi
;; * `cramers-c`, `cramers-v`, `cramers-v-corrected`
;; * `cohens-w`, `tschuprows-t`
;; * `cohens-kappa`, `weighted-kappa`
;; :::

;; ## Binary Classification Metrics

;; Metrics derived from a 2x2 confusion matrix for binary classification tasks.

;; ::: {.callout-tip title="Defined functions"}
;; * `->confusion-matrix`
;; * `binary-measures-all`, `binary-measures`
;; :::

;; ## Effect Size

;; Measures quantifying the magnitude of a phenomenon or relationship.

;; ### Difference Family

;; Based on differences between means.

;; ::: {.callout-tip title="Defined functions"}
;; * `cohens-d`, `cohens-d-corrected`
;; * `hedges-g`, `hedges-g-corrected`, `hedges-g*`
;; * `glass-delta`
;; :::

;; ### Ratio Family

;; Based on the ratio of means.

;; ::: {.callout-tip title="Defined functions"}
;; * `means-ratio`, `means-ratio-corrected`
;; :::

;; ### Ordinal / Non-parametric Family

;; Based on ranks or ordinal comparisons.

;; ::: {.callout-tip title="Defined functions"}
;; * `cliffs-delta`
;; * `ameasure`, `wmw-odds`
;; :::

;; ### Overlap Family

;; Based on distribution overlap.

;; ::: {.callout-tip title="Defined functions"}
;; * `p-overlap`
;; * `cohens-u1-normal`, `cohens-u2-normal`, `cohens-u3-normal`
;; * `cohens-u1`, `cohens-u2`, `cohens-u3`
;; :::

;; ### Correlation / Association Family

;; Based on correlation coefficients or explained variance.

;; ::: {.callout-tip title="Defined functions"}
;; * `pearson-r`, `r2-determination`
;; * `eta-sq`, `omega-sq`, `epsilon-sq`
;; * `cohens-f2`, `cohens-f`, `cohens-q`
;; :::

;; ### Rank-based Association Family

;; Effect sizes for rank-based tests like Kruskal-Wallis.

;; ::: {.callout-tip title="Defined functions"}
;; * `rank-eta-sq`, `rank-epsilon-sq`
;; :::

;; ## Statistical Tests

;; Functions for hypothesis testing.

;; ### Normality and Shape Tests

;; Tests checking if data follows a normal distribution or has specific shape characteristics.

;; ::: {.callout-tip title="Defined functions"}
;; * `skewness-test`
;; * `kurtosis-test`
;; * `normality-test` ;; D'Agostino-Pearson K² test
;; * `jarque-bera-test`
;; * `bonett-seier-test` ;; Test based on Geary's kurtosis
;; :::

;; ### Binomial Tests

;; Tests related to binomial distributions (success/failure).

;; ::: {.callout-tip title="Defined functions"}
;; * `binomial-test`
;; * `binomial-ci` ;; Confidence intervals for binomial proportions
;; :::

;; ### Location Tests (T/Z Tests)

;; Tests comparing means of one or two groups.

;; ::: {.callout-tip title="Defined functions"}
;; * `t-test-one-sample`, `z-test-one-sample`
;; * `t-test-two-samples`, `z-test-two-samples` ;; (Includes paired and unpaired, Welch's)
;; * `p-value` ;; Helper function
;; :::

;; ### Variance Tests

;; Tests comparing variances of two or more groups.

;; ::: {.callout-tip title="Defined functions"}
;; * `f-test` ;; Test for equality of two variances
;; * `levene-test` ;; Test for homogeneity of variances (mean-based)
;; * `brown-forsythe-test` ;; Test for homogeneity of variances (median-based)
;; * `fligner-killeen-test` ;; Non-parametric test for homogeneity of variances
;; :::

;; ### Goodness-of-Fit and Independence Tests

;; Tests comparing observed data to expected distributions or testing independence.

;; ::: {.callout-tip title="Defined functions"}
;; * `power-divergence-test` ;; General framework
;; * `chisq-test` ;; Chi-squared test (λ=1)
;; * `multinomial-likelihood-ratio-test` ;; G-test (λ=0)
;; * `minimum-discrimination-information-test` ;; (λ=-1)
;; * `neyman-modified-chisq-test` ;; (λ=-2)
;; * `freeman-tukey-test` ;; (λ=-0.5)
;; * `cressie-read-test` ;; (λ=2/3)
;; * `ad-test-one-sample` ;; Anderson-Darling test
;; * `ks-test-one-sample`, `ks-test-two-samples` ;; Kolmogorov-Smirnov test
;; :::

;; ### ANOVA and Rank Sum Tests

;; Tests comparing means/distributions across multiple groups.

;; ::: {.callout-tip title="Defined functions"}
;; * `one-way-anova-test`
;; * `kruskal-test` ;; Kruskal-Wallis rank sum test
;; :::

;; ### Autocorrelation Tests

;; Tests for autocorrelation in sequences (e.g., residuals).

;; ::: {.callout-tip title="Defined functions"}
;; * `durbin-watson`
;; :::

;; ## Time Series Analysis

;; Functions specifically for analyzing sequential data.

;; ::: {.callout-tip title="Defined functions"}
;; * `acf` ;; Autocorrelation function
;; * `pacf` ;; Partial autocorrelation function
;; * `acf-ci`, `pacf-ci` ;; ACF/PACF with confidence intervals
;; :::

;; ## Histograms

;; Functions for creating histograms from data.

;; ::: {.callout-tip title="Defined functions"}
;; * `histogram`
;; * `estimate-bins`
;; :::


;; ## Bootstrap

;; Resampling methods for estimating statistics and confidence intervals.
;; Functions are primarily defined in the `fastmath.stats.bootstrap` namespace.

;; ::: {.callout-tip title="Defined functions"}
;; * `bootstrap` ;; Core function for generating bootstrap samples
;; * `jackknife`, `jackknife+` ;; Jackknife resampling methods
;; * `bootstrap-stats` ;; Helper to calculate stats on bootstrap samples
;; * `ci-normal` ;; Normal approximation CI
;; * `ci-basic` ;; Basic (reverse percentile) CI
;; * `ci-percentile` ;; Percentile CI
;; * `ci-bc` ;; Bias-corrected CI
;; * `ci-bca` ;; Bias-corrected and accelerated CI
;; * `ci-studentized` ;; Studentized CI
;; * `ci-t` ;; Student's T CI
;; :::

;; ## Reference

(codox/make-public-fns-table-clay 'fastmath.stats)
(codox/make-public-fns-table-clay 'fastmath.stats.bootstrap)
