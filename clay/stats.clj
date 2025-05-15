(ns stats
  (:require [fastmath.stats :as stats]
            [fastmath.stats.bootstrap :as boot]
            [fastmath.dev.ggplot :as gg]
            [fastmath.dev.codox :as codox]
            [fastmath.dev.clay :as utls]
            [fastmath.dev.dataset :as ds]
            [fastmath.core :as m]
            [scicloj.kindly.v4.kind :as kind]
            [fastmath.fields.n :as n]
            [fastmath.stats :as stat]))

;; # Statistics {.unnumbered}

;; Statistical functions

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

;; *   **k-th Central Moment**: $\mu_k = E[(X - \mu)^k] \approx \frac{1}{n} \sum_{i=1}^{n} (x_i - \bar{x})^k$. Calculated when `center` is `nil` (default) and `:mean?` is `true` (default).
;; *   **k-th Raw Moment** (about origin): $\mu'_k = E[X^k] \approx \frac{1}{n} \sum_{i=1}^{n} x_i^k$. Calculated if `center` is `0.0`.
;; *   **k-th Absolute Central Moment**: $E[|X - \mu|^k] \approx \frac{1}{n} \sum_{i=1}^{n} |x_i - \bar{x}|^k$. Calculated if `:absolute?` is `true`.
;; *   **Normalization**: If `:normalize?` is `true`, the moment is divided by $\sigma^k$ (where $\sigma$ is the standard deviation), yielding a scale-invariant measure. For example, the 3rd normalized central moment is related to skewness, and the 4th to kurtosis.
;; *   If `:mean?` is `false`, the function returns the sum $\sum (x_i - c)^k$ (or sum of absolute values) instead of the mean.

;; The `order` parameter specifies $k$. For example, the 2nd central moment is the variance.

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

;; *   **Robust (less sensitive to outliers):**

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

;; ### L-moments

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

;; *   **Basic Range:** Functions like `span` ($max - min$) and `extent` (providing $[min, max]$ and optionally the mean) offer simple measures of the total spread of the data.
;;
;; *   **Interquartile Range:** `iqr` ($Q_3 - Q_1$) specifically measures the spread of the middle 50% of the data, providing a robust alternative to the total range.
;;
;; *   **Symmetric Spread Intervals:** Functions ending in `-extent` such as `stddev-extent`, `mad-extent`, and `sem-extent` define intervals typically centered around the mean or median. They represent a range defined by adding/subtracting a multiple (usually 1) of a measure of dispersion (Standard Deviation, Median Absolute Deviation, or Standard Error of the Mean) from the central point.
;;
;; *   **Quantile-Based Intervals:** `percentile-extent`, `quantile-extent`, `pi`, `pi-extent`, and `hpdi-extent` define intervals based on quantiles or percentiles of the data. These functions capture specific ranges containing a certain percentage of the data points (e.g., the middle 95% defined by quantiles 0.025 and 0.975). `hpdi-extent` calculates the shortest interval containing a given proportion of data, based on empirical density.
;;
;; *   **Box Plot Boundaries:** `adjacent-values` (LAV, UAV) and fence functions (`inner-fence-extent`, `outer-fence-extent`) calculate specific bounds based on quartiles and multiples of the IQR. These are primarily used in box plot visualization and as a conventional method for identifying potential outliers.
;;
;; *   **Confidence and Prediction Intervals:** `ci`, `percentile-bc-extent`, and `percentile-bca-extent` provide inferential intervals. `ci` estimates a confidence interval for the population mean using the t-distribution. `percentile-bc-extent` and `percentile-bca-extent` (Bias-Corrected and Bias-Corrected Accelerated) are advanced bootstrap methods for estimating confidence intervals for statistics, offering robustness against non-normality and bias.

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
                           :title "Various extents"}))

;; ## Outlier Detection
;;
;; Outlier detection involves identifying data points that are significantly different from other observations. Outliers can distort statistical analyses and require careful handling. `fastmath.stats` provides functions to find and optionally remove such values based on the Interquartile Range (IQR) method.
;;
;; ::: {.callout-tip title="Defined functions"}
;; * `outliers`
;; * `remove-outliers`
;; :::
;;
;; Both functions use the inner fence rule based on the IQR:
;;
;; *   **Lower Inner Fence (LIF):** $Q_1 - 1.5 \times IQR$
;; *   **Upper Inner Fence (UIF):** $Q_3 + 1.5 \times IQR$
;;
;; Where $Q_1$ is the first quartile (25th percentile) and $Q_3$ is the third quartile (75th percentile). Points falling below the LIF or above the UIF are considered outliers.
;;
;; *   `outliers`: Returns a sequence containing only the data points identified as outliers.
;; *   `remove-outliers`: Returns a sequence containing the data points from the original sequence *excluding* those identified as outliers.
;;
;; Both functions accept an optional `estimation-strategy` keyword (see [[quantile]]) to control how quartiles are calculated, which affects the fence boundaries.
;;
;; Let's find the outliers in the `residual-sugar` data and then see how the data looks after removing them.
;;
(utls/examples-note
  (stats/outliers residual-sugar)
  (count residual-sugar)
  (count (stats/remove-outliers residual-sugar)))

;; We can visualize the effect of removing outliers on the distribution shape.
(kind/table
 [[(gg/->image (gg/density residual-sugar {:title "Original Residual Sugar Distribution"}))
   (gg/->image (gg/density (stats/remove-outliers residual-sugar) {:title "Residual Sugar Distribution (Outliers Removed)"}))]])

;; ## Data Transformation

;; Functions to modify data (scaling, normalizing, transforming).

;; ::: {.callout-tip title="Defined functions"}
;; * `standardize`, `robust-standardize`, `demean`
;; * `rescale`
;; * `trim`, `trim-lower`, `trim-upper`, `winsor`
;; * `box-cox-infer-lambda`, `box-cox-transformation`
;; * `yeo-johnson-infer-lambda`, `yeo-johnson-transformation`
;; :::

;; ## Correlation and Covariance

;; Measures of the relationship between two or more variables.

;; ::: {.callout-tip title="Defined functions"}
;; * `covariance`, `correlation`
;; * `pearson-correlation`, `spearman-correlation`, `kendall-correlation`
;; * `coefficient-matrix`, `correlation-matrix`, `covariance-matrix`
;; :::

;; ## Distance and Similarity Metrics

;; Measures of distance, error, or similarity between sequences or distributions.

;; ::: {.callout-tip title="Defined functions"}
;; * `dissimilarity` ;; (Covers many specific distances like Euclidean, KL, JSD etc.)
;; * `similarity` ;; (Covers many specific similarities like Cosine, Jaccard etc.)
;; * `me`, `mae`, `mape` ;; Error metrics
;; * `rss`, `mse`, `rmse` ;; Squared error metrics
;; * `r2` ;; Coefficient of determination
;; * `count=`, `L0`, `L1`, `L2sq`, `L2`, `LInf` ;; L-norms / distances
;; * `psnr` ;; Peak signal-to-noise ratio
;; :::

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
