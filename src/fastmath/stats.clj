(ns fastmath.stats
  "Namespace provides a comprehensive collection of functions for
  performing statistical analysis in Clojure. It focuses on providing efficient
  implementations for common statistical tasks, leveraging fastmath's underlying
  numerical capabilities.

  This namespace covers a wide range of statistical methods, including:

  *   **Descriptive Statistics**: Measures of central tendency (mean, median, mode, expectile),
      dispersion (variance, standard deviation, MAD, SEM), and shape (skewness, kurtosis, L-moments).
  *   **Quantiles and Percentiles**: Functions for calculating percentiles, quantiles, and the median,
      including weighted versions and various estimation strategies.
  *   **Intervals and Extents**: Methods for defining ranges within data, such as span, IQR,
      standard deviation/MAD/SEM extents, percentile/quantile intervals, prediction intervals (PI, HPDI),
      and fence boundaries for outlier detection.
  *   **Outlier Detection**: Functions for identifying data points outside conventional fence boundaries.
  *   **Data Transformation**: Utilities for scaling, centering, trimming, winsorizing,
      and applying power transformations (Box-Cox, Yeo-Johnson) to data.
  *   **Correlation and Covariance**: Measures of the linear and monotonic relationship
      between two or more variables (Pearson, Spearman, Kendall), and functions for
      generating covariance and correlation matrices.
  *   **Distance and Similarity Metrics**: Functions for quantifying differences or
      likeness between data sequences or distributions, including error metrics (MAE, MSE, RMSE),
      L-p norms, and various distribution dissimilarity/similarity measures.
  *   **Contingency Tables**: Functions for creating, analyzing, and deriving measures
      of association and agreement (Cramer's V, Cohen's Kappa) from contingency tables,
      including specialized functions for 2x2 tables.
  *   **Binary Classification Metrics**: Functions for generating confusion matrices
      and calculating a wide array of performance metrics (Accuracy, Precision, Recall, F1, MCC, etc.).
  *   **Effect Size**: Measures quantifying the magnitude of statistical effects,
      including difference-based (Cohen's d, Hedges' g, Glass's delta), ratio-based,
      ordinal/non-parametric (Cliff's Delta, Vargha-Delaney A), and overlap-based (Cohen's U, p-overlap),
      as well as measures related to explained variance (Eta-squared, Omega-squared, Cohen's f²).
  *   **Statistical Tests**: Functions for performing hypothesis tests, including:
      -   Normality and Shape tests (Skewness, Kurtosis, D'Agostino-Pearson K², Jarque-Bera, Bonett-Seier).
      -   Binomial tests and confidence intervals.
      -   Location tests (one-sample and two-sample T/Z tests, paired/unpaired).
      -   Variance tests (F-test, Levene's, Brown-Forsythe, Fligner-Killeen).
      -   Goodness-of-Fit and Independence tests (Power Divergence family including Chi-squared, G-test; AD/KS tests).
      -   ANOVA and Rank Sum tests (One-way ANOVA, Kruskal-Wallis).
      -   Autocorrelation tests (Durbin-Watson).
  *   **Time Series Analysis**: Functions for analyzing the dependence structure of
      time series data, such as Autocorrelation (ACF) and Partial Autocorrelation (PACF).
  *   **Histograms**: Functions for computing histograms and estimating optimal binning strategies.

  This namespace aims to provide a robust set of statistical tools for data analysis
  and modeling within the Clojure ecosystem."
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.distance :as d]
            [fastmath.vector :as v]
            [fastmath.interpolation.step :as step-interp]
            [fastmath.interpolation.linear :as linear-interp]
            [fastmath.optimization.lbfgsb :as lbfgsb]
            [fastmath.optimization :as opt]
            [fastmath.kernel.density :as kd]
            [fastmath.special :as special]
            [fastmath.solver :as solver]

            [fastmath.stats.bins :as bins]
            [fastmath.stats.binary :as binary]
            [fastmath.stats.logmean :as logmean])
  (:import [org.apache.commons.math3.stat StatUtils]
           [org.apache.commons.math3.stat.descriptive.rank Percentile Percentile$EstimationType]
           [org.apache.commons.math3.stat.descriptive.moment Kurtosis Skewness]
           [org.apache.commons.math3.stat.correlation KendallsCorrelation SpearmansCorrelation PearsonsCorrelation]
           [org.apache.commons.math3.stat.regression SimpleRegression]
           [org.apache.commons.math3.analysis.integration RombergIntegrator]
           [org.apache.commons.math3.analysis UnivariateFunction]
           [fastmath.java Array]
           [fastmath.vector Vec2 Vec3]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn minimum
  "Finds the smallest value in a sequence of numbers.

  Accepts either a `double-array` (handled with an optimized, array-specific path) or any other sequence of numbers.

  Parameters:

  - `vs` (sequence of numbers or `double-array`): Input data.

  Returns the minimum value as a double. Throws an exception when `vs` is empty.

  See also [[maximum]], [[extent]]."
  ^double [vs]
  (if (= (type vs) m/double-array-type)
    (Array/min ^doubles vs)
    (reduce m/min vs)))

(defn maximum
  "Finds the largest value in a sequence of numbers.

  Accepts either a `double-array` (handled with an optimized, array-specific path) or any other sequence of numbers.

  Parameters:

  - `vs` (sequence of numbers or `double-array`): Input data.

  Returns the maximum value as a double. Throws an exception when `vs` is empty.

  See also [[minimum]], [[extent]]."
  ^double [vs]
  (if (= (type vs) m/double-array-type)
    (Array/max ^doubles vs)
    (reduce m/max vs)))

(def ^{:doc "Map of available estimation strategies usable by [[percentile]] and [[quantile]].

  Each strategy defines how a percentile position is computed from a sorted sample and how the value at a (possibly fractional) position is estimated, following Hyndman and Fan's classification of sample quantile methods. Pass the corresponding keyword as the `estimation-strategy` argument.

  Available strategies:

  - `:legacy` - Apache Commons Math's original method, close to `:r6` but with min/max returned for out-of-range positions. Default when no strategy is specified.
  - `:r1` - Inverse of the empirical CDF, no interpolation (a step function).
  - `:r2` - Like `:r1`, but averages the two candidate values at discontinuities.
  - `:r3` - Nearest, rounding to even (SAS default).
  - `:r4` - Linear interpolation of the empirical CDF.
  - `:r5` - Hazen's method, linear interpolation through the midpoints of the order statistics.
  - `:r6` - Weibull's method, linear interpolation of the expectations of the order statistics (used by SPSS and Minitab).
  - `:r7` - Linear interpolation of the modes of the order statistics (default in R and NumPy).
  - `:r8` - Linear interpolation giving an estimate that is approximately median-unbiased regardless of the underlying distribution.
  - `:r9` - Linear interpolation giving an estimate that is approximately unbiased for a normal distribution.

  See also [[percentile]], [[quantile]]."}
  estimation-strategies-list {:legacy Percentile$EstimationType/LEGACY
                              :r1 Percentile$EstimationType/R_1
                              :r2 Percentile$EstimationType/R_2
                              :r3 Percentile$EstimationType/R_3
                              :r4 Percentile$EstimationType/R_4
                              :r5 Percentile$EstimationType/R_5
                              :r6 Percentile$EstimationType/R_6
                              :r7 Percentile$EstimationType/R_7
                              :r8 Percentile$EstimationType/R_8
                              :r9 Percentile$EstimationType/R_9})

(defn- kahan-step
  ^Vec2 [^Vec2 b ^double v]
  (let [av (m/- v (.y b))
        nsum (m/+ av (.x b))]
    (Vec2. nsum (m/- nsum (.x b) av))))

(defn- neumayer-step
  ^Vec2 [^Vec2 b ^double v]
  (let [t (m/+ (.x b) v)]
    (Vec2. t (m/+ (.y b) (if (m/>= (m/abs (.x b)) (m/abs v))
                           (m/+ (m/- (.x b) t) v)
                           (m/+ (m/- v t) (.x b)))))))

(defn- klein-step
  ^Vec3 [^Vec3 b ^double v]
  (let [t (m/+ (.x b) v)
        c (if (m/>= (m/abs (.x b)) (m/abs v))
            (m/+ (m/- (.x b) t) v)
            (m/+ (m/- v t) (.x b)))
        nt (m/+ (.y b) c)
        cc (if (m/>= (m/abs (.y b)) (m/abs c))
             (m/+ (m/- (.y b) nt) c)
             (m/+ (m/- c nt) (.y b)))]
    (Vec3. t nt (m/+ (.z b) cc))))

(defn sum
  "Calculates the sum of all values in `vs`.

  By default, plain summation is used (an optimized array-specific path for `double-array`, or a simple reduction otherwise), which is fast but can accumulate floating-point rounding error for long sequences or values of widely differing magnitude. An optional compensated summation algorithm can be selected instead to improve numerical accuracy at the cost of extra computation.

  Parameters:

  - `vs` (sequence of numbers or `double-array`): Values to sum.
  - `compensation-method` (keyword, optional): Compensated summation algorithm. One of `:kahan`, `:neumayer` (Neumaier's improved Kahan algorithm) or `:klein` (second-order compensated summation). Any other value falls back to plain summation.

  Returns the sum as a double."
  (^double [vs]
   (if (= (type vs) m/double-array-type)
     (Array/sum ^doubles vs)
     (reduce m/+ vs)))
  (^double [vs compensation-method]
   (case compensation-method
     :kahan (let [^Vec2 r (reduce kahan-step (Vec2. 0.0 0.0) vs)] (.x r))
     :neumayer (v/sum (reduce neumayer-step (Vec2. 0.0 0.0) vs))
     :klein (v/sum (reduce klein-step (Vec3. 0.0 0.0 0.0) vs)))))

;; https://www.amherst.edu/media/view/129116/original/Sample+Quantiles.pdf

(defn percentile
  "Calculates the p-th percentile of a sequence `vs`.

  The percentile `p` is a value between 0 and 100, inclusive.

  An optional `estimation-strategy` keyword can be provided to specify the
  method used for estimating the percentile, particularly how interpolation is
  handled when the desired percentile falls between data points in the sorted
  sequence.

  Available `estimation-strategy` values:

  - `:legacy` (Default): The original method used in Apache Commons Math.
  - `:r1` through `:r9`: Correspond to the nine quantile estimation algorithms recommended by Hyndman and Fan (1996). Each strategy differs slightly in how it calculates the index (e.g., using `np` or `(n+1)p`) and how it interpolates between points.

  For detailed mathematical descriptions of each estimation strategy, refer to
  the [Apache Commons Math Percentile documentation](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.EstimationType.html).

  See also [[quantile]] (which uses a 0.0-1.0 range) and [[percentiles]]."
  (^double [vs ^double p]
   (if (m/zero? p)
     (minimum vs)
     (StatUtils/percentile (m/seq->double-array vs) p)))
  (^double [vs ^double p estimation-strategy]
   (if (m/zero? p)
     (minimum vs)
     (let [^Percentile perc (.withEstimationType (Percentile.) (get estimation-strategies-list estimation-strategy Percentile$EstimationType/LEGACY))]
       (.evaluate perc (m/seq->double-array vs) p)))))

(defn percentiles
  "Calculates the sequence of p-th percentiles of a sequence `vs`.

  Percentiles `ps` is sequence of values between 0 and 100, inclusive.

  An optional `estimation-strategy` keyword can be provided to specify the
  method used for estimating the percentile, particularly how interpolation is
  handled when the desired percentile falls between data points in the sorted
  sequence.

  Available `estimation-strategy` values:

  - `:legacy` (Default): The original method used in Apache Commons Math.
  - `:r1` through `:r9`: Correspond to the nine quantile estimation algorithms recommended by Hyndman and Fan (1996). Each strategy differs slightly in how it calculates the index (e.g., using `np` or `(n+1)p`) and how it interpolates between points.

  For detailed mathematical descriptions of each estimation strategy, refer to
  the [Apache Commons Math Percentile documentation](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.EstimationType.html).

  See also [[quantiles]] (which uses a 0.0-1.0 range) and [[percentile]]."
  ([vs] (percentiles vs [25 50 75 100]))
  ([vs ps] (percentiles vs ps nil))
  ([vs ps estimation-strategy]
   (let [^Percentile perc (.withEstimationType (Percentile.) (or (estimation-strategies-list estimation-strategy) Percentile$EstimationType/LEGACY))
         d (m/seq->double-array vs)]
     (.setData perc d)
     (mapv (fn [^double p] (if (m/zero? p) (minimum d) (.evaluate perc p))) ps))))

(defn quantile
  "Calculates the q-th quantile of a sequence `vs`.

  The quantile `q` is a value between 0.0 and 1.0, inclusive.

  An optional `estimation-strategy` keyword can be provided to specify the
  method used for estimating the quantile, particularly how interpolation is
  handled when the desired quantile falls between data points in the sorted
  sequence.

  Available `estimation-strategy` values:

  - `:legacy` (Default): The original method used in Apache Commons Math.
  - `:r1` through `:r9`: Correspond to the nine quantile estimation algorithms recommended by Hyndman and Fan (1996). Each strategy differs slightly in how it calculates the index (e.g., using `np` or `(n+1)p`) and how it interpolates between points.

  For detailed mathematical descriptions of each estimation strategy, refer to
  the [Apache Commons Math Percentile documentation](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.EstimationType.html).

  See also [[percentile]] (which uses a 0-100 range) and [[quantiles]]."
  (^double [vs ^double q]
   (percentile vs (m/constrain (m/* q 100.0) 0.0 100.0)))
  (^double [vs ^double q estimation-strategy]
   (percentile vs (m/constrain (m/* q 100.0) 0.0 100.0) estimation-strategy)))

(defn quantiles
  "Calculates the sequence of q-th quantiles of a sequence `vs`.

  Quantiles `q` is a sequence of values between 0.0 and 1.0, inclusive.

  An optional `estimation-strategy` keyword can be provided to specify the
  method used for estimating the quantile, particularly how interpolation is
  handled when the desired quantile falls between data points in the sorted
  sequence.

  Available `estimation-strategy` values:

  - `:legacy` (Default): The original method used in Apache Commons Math.
  - `:r1` through `:r9`: Correspond to the nine quantile estimation algorithms recommended by Hyndman and Fan (1996). Each strategy differs slightly in how it calculates the index (e.g., using `np` or `(n+1)p`) and how it interpolates between points.

  For detailed mathematical descriptions of each estimation strategy, refer to
  the [Apache Commons Math Percentile documentation](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.EstimationType.html).

  See also [[percentiles]] (which uses a 0-100 range) and [[quantile]]."
  ([vs] (quantiles vs [0.25 0.5 0.75 1.0]))
  ([vs qs]
   (percentiles vs (map (fn [^double v] (m/constrain (m/* v 100.0) 0.0 100.0)) qs)))
  ([vs qs estimation-strategy]
   (percentiles vs (map (fn [^double v] (m/constrain (m/* v 100.0) 0.0 100.0)) qs) estimation-strategy)))

(defn- wquantile-interpolator
  [vs ws method]
  (let [sorted (sort-by first (map vector vs ws))
        probabilities (map second sorted)
        data (map first sorted)
        data (conj data (first data))
        wsum (sum probabilities)
        weights (conj (reductions m/+ (map (fn [^double p] (m// p wsum)) probabilities)) 0.0)]
    (case method
      :linear (linear-interp/linear weights data)
      :average (let [interp1 (step-interp/step-before weights data)
                     interp2 (step-interp/step-after weights data)]
                 (fn [^double x] (m/* 0.5 (m/+ (double (interp1 x))
                                              (double (interp2 x))))))
      :step (step-interp/step-before weights data))))

;; based on spatstat.geom::weighted.quantile

(defn wquantile
  "Calculates the q-th weighted quantile of a sequence `vs` with corresponding weights `ws`.

  The quantile `q` is a value between 0.0 and 1.0, inclusive.

  The calculation involves constructing a weighted empirical cumulative distribution
  function (ECDF) and interpolating to find the value at quantile `q`.

  Parameters:

  - `vs`: Sequence of data values.
  - `ws`: Sequence of corresponding non-negative weights. Must have the same count as `vs`.
  - `q`: The quantile level (0.0 < q <= 1.0).
  - `method` (optional keyword): Specifies the interpolation method used when `q` falls
    between points in the weighted ECDF. Defaults to `:linear`.
      - `:linear`: Performs linear interpolation between the data values corresponding
        to the cumulative weights surrounding `q`.
      - `:step`: Uses a step function (specifically, step-before) based on the
        weighted ECDF. The result is the data value whose cumulative weight range
        includes `q`.
      - `:average`: Computes the average of the step-before and step-after
        interpolation methods. Useful when `q` corresponds exactly to a cumulative
        weight boundary.

  See also: [[wmedian]], [[wquantiles]], [[quantile]]."
  (^double [vs ws ^double q] (wquantile vs ws q :linear))
  (^double [vs ws ^double q method]
   (let [interp (wquantile-interpolator vs ws method)]
     (interp q))))

(defn wquantiles
  "Calculates the sequence of q-th weighted quantiles of a sequence `vs` with corresponding weights `ws`.

  Quantiles `qs` is a sequence of values between 0.0 and 1.0, inclusive.

  The calculation involves constructing a weighted empirical cumulative distribution
  function (ECDF) and interpolating to find the value at quantiles `qs`.

  Parameters:

  - `vs`: Sequence of data values.
  - `ws`: Sequence of corresponding non-negative weights. Must have the same count as `vs`.
  - `qs`: Sequence of quantiles level (0.0 < q <= 1.0).
  - `method` (optional keyword): Specifies the interpolation method used when `qs` falls
    between points in the weighted ECDF. Defaults to `:linear`.
      - `:linear`: Performs linear interpolation between the data values corresponding
        to the cumulative weights surrounding `q`.
      - `:step`: Uses a step function (specifically, step-before) based on the
        weighted ECDF. The result is the data value whose cumulative weight range
        includes `q`.
      - `:average`: Computes the average of the step-before and step-after
        interpolation methods. Useful when `q` corresponds exactly to a cumulative
        weight boundary.

  See also: [[wquantile]], [[quantiles]]."
  ([vs ws] (wquantiles vs ws [0.25 0.5 0.75 1.0]))
  ([vs ws qs] (wquantiles vs ws qs :linear))
  ([vs ws qs method]
   (let [interp (wquantile-interpolator vs ws method)]
     (mapv interp qs))))

(defn wmedian
  "Calculates median of a sequence `vs` with corresponding weights `ws`.

  Parameters:

  - `vs`: Sequence of data values.
  - `ws`: Sequence of corresponding non-negative weights. Must have the same count as `vs`.
  - `method` (optional keyword): Specifies the interpolation method used when `qs` falls
    between points in the weighted ECDF. Defaults to `:linear`.
      - `:linear`: Performs linear interpolation between the data values corresponding
        to the cumulative weights surrounding `q=0.5`.
      - `:step`: Uses a step function (specifically, step-before) based on the
        weighted ECDF. The result is the data value whose cumulative weight range
        includes `q=0.5`.
      - `:average`: Computes the average of the step-before and step-after
        interpolation methods.

  See also: [[wquantile]], [[quantile]]."
  (^double [vs ws] (wquantile vs ws 0.5))
  (^double [vs ws method] (wquantile vs ws 0.5 method)))

(defn median
  "Calculates median of a sequence `vs`.

  An optional `estimation-strategy` keyword can be provided to specify the
  method used for estimating the quantile, particularly how interpolation is
  handled when the desired quantile falls between data points in the sorted
  sequence.

  Available `estimation-strategy` values:

  - `:legacy` (Default): The original method used in Apache Commons Math.
  - `:r1` through `:r9`: Correspond to the nine quantile estimation algorithms
      recommended by Hyndman and Fan (1996). Each strategy differs slightly in how it calculates the index (e.g., using `np` or `(n+1)p`) and how it interpolates between points.

  For detailed mathematical descriptions of each estimation strategy, refer to
  the [Apache Commons Math Percentile documentation](http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/stat/descriptive/rank/Percentile.EstimationType.html).

  See also [[quantile]], [[median-3]]"
  (^double [vs estimation-strategy]
   (percentile vs 50.0 estimation-strategy))
  (^double [vs]
   (percentile vs 50.0)))

(defn median-3
  "Median of three values. See [[median]]."
  ^double [^double a ^double b ^double c]
  (m/max (m/min a b) (m/min (m/max a b) c)))

(defn mean
  "Calculates the arithmetic mean (average) of a sequence `vs`.

  If `weights` are provided, calculates the weighted arithmetic mean.

  Parameters:

  - `vs`: Sequence of numbers.
  - `weights` (optional): Sequence of non-negative weights corresponding to `vs`.
    Must have the same count as `vs`.

  Returns the calculated mean as a double.

  See also [[geomean]], [[harmean]], [[powmean]], [[median]]."
  (^double [vs] (StatUtils/mean (m/seq->double-array vs)))
  (^double [vs weights] (m// (v/dot vs weights) (sum weights))))

(defn geomean
  "Calculates the geometric mean of a sequence `vs`.

  The geometric mean is suitable for averaging ratios or rates of change and requires
  all values in the sequence to be positive. It is calculated as the n-th root
  of the product of n numbers.

  Parameters:

  - `vs`: Sequence of numbers. Non-positive values will result in `NaN` or `0.0` due
          to the internal use of `log`.
  - `weights` (optional): Sequence of non-negative weights corresponding to `vs`.
    Must have the same count as `vs`.

  Returns the calculated geometric mean as a double.

  See also [[mean]], [[harmean]], [[powmean]]."
  (^double [vs] (m/exp (mean (map m/log vs))))
  (^double [vs weights] (m/exp (mean (map m/log vs) weights))))

(defn harmean
  "Calculates the harmonic mean of a sequence `vs`.

  The harmonic mean is the reciprocal of the arithmetic mean of the reciprocals
  of the observations.

  Parameters:

  - `vs`: Sequence of numbers. Values must be non-zero.
  - `weights` (optional): Sequence of non-negative weights corresponding to `vs`.
    Must have the same count as `vs`.

  Returns the calculated harmonic mean as a double.

  See also [[mean]], [[geomean]], [[powmean]]."
  (^double [vs] (m// (mean (map m// vs))))
  (^double [vs weights] (m// (mean (map m// vs) weights))))

(defn powmean
  "Calculates the generalized power mean (also known as the Hölder mean) of a sequence `vs`.

  The power mean is a generalization of the Pythagorean means (arithmetic, geometric, harmonic)
  and other means like the quadratic mean (RMS). It is defined for a non-zero real number `power`.

  Parameters:

  - `vs`: Sequence of numbers. Constraints depend on the `power`:
    - For `power > 0`, values should be non-negative.
    - For `power = 0`, values must be positive (reduces to geometric mean).
    - For `power < 0`, values must be positive and non-zero.
  - `weights` (optional): Sequence of non-negative weights corresponding to `vs`.
    Must have the same count as `vs`.
  - `power` (double): The exponent defining the mean.

  Special Cases:

  - `power = 0`: Returns the [[geomean]].
  - `power = 1`: Returns the arithmetic [[mean]].
  - `power = -1`: Equivalent to the [[harmean]]. (Handled by the general formula)
  - `power = 2`: Returns the Root Mean Square (RMS) or quadratic mean.
  - `power = inf`: Returns maximum.
  - `power = -inf`: Returms minimum.
  - The implementation includes optimized paths for `power` values 1/3, 0.5, 2, and 3.

  Returns the calculated power mean as a double.

  See also [[mean]], [[geomean]], [[harmean]]."
  (^double [vs ^double power]
   (cond
     (m/zero? power) (geomean vs)
     (m/one? power) (mean vs)
     (m/== power m/THIRD) (m/cb (mean (map m/cbrt vs)))
     (m/== power 0.5) (m/sq (mean (map m/sqrt vs)))
     (m/== power 2.0) (m/sqrt (mean (map m/sq vs)))
     (m/== power 3.0) (m/cbrt (mean (map m/cb vs)))
     (m/pos-inf? power) (maximum vs)
     (m/neg-inf? power) (minimum vs)
     :else (m/pow (mean (map (fn [^double v] (m/pow v power)) vs)) (m// power))))
  (^double [vs weights ^double power]
   (cond
     (m/zero? power) (geomean vs weights)
     (m/one? power) (mean vs weights)
     (m/== power m/THIRD) (m/cb (mean (map m/cbrt vs) weights))
     (m/== power 0.5) (m/sq (mean (map m/sqrt vs) weights))
     (m/== power 2.0) (m/sqrt (mean (map m/sq vs) weights))
     (m/== power 3.0) (m/cbrt (mean (map m/cb vs) weights))
     (m/pos-inf? power) (maximum vs)
     (m/neg-inf? power) (minimum vs)
     :else (m/pow (mean (map (fn [^double v] (m/pow v power)) vs) weights) (m// power)))))

;; https://www.survo.fi/papers/logmean.pdf

(defn logmean
  "Calculates the generalized logarithmic mean of a sequence of positive numbers.

  The logarithmic mean generalizes the two-argument logarithmic mean `L(x,y) = (x-y)/(ln(x)-ln(y))`
  to `n` arguments. For `n=1` the single value is returned as-is; for `n=2` the classical
  two-argument formula is used; for `n>=3` one of two methods is applied.

  Parameters:

  - `xs` - a sequence of positive numbers (length `n >= 1`)
  - `opts` - optional map of configuration keys:
    - `:method` - algorithm to use, one of `:integral` (default) or `:mean-value`
      (aliases `:divided`, `:divided-differences`); `:integral` uses the integral
      representation of the generalized log-mean; `:mean-value` uses divided differences
    - `:tol` - convergence tolerance for the `:integral` series expansion (default `1.0e-15`)
    - `:max-iters` - maximum number of series iterations for the `:integral` method (default `1000`)

  Returns the generalized logarithmic mean as a `double`. Returns `##NaN` for an empty sequence.
  All values in `xs` must be strictly positive.

  Note 1: `:mean-value` can be unstable for large number of entries (`n >= 100`) or when differences are very small.
  Note 2: `:mean-value` and `:integral` are two different interpretations (definitions) and produce different results.

  See also [[mean]], [[geometric-mean]], [[harmonic-mean]].

  Generalized integral method is made with Opus 4.8"
  (^double [xs] (logmean xs nil))
  (^double [xs {:keys [^double tol ^long max-iters method]
                :or {tol 1.0e-15 max-iters 1000 method :integral}}]
   (let [xs (vec xs)
         n (count xs)]
     (case n
       0 ##NaN
       1 (xs 0)
       2 (logmean/logmean2 xs)
       (case method
         (:mean-value :divided :divided-differences) (if (m/== n 3)
                                                       (logmean/logmean3-mean-value xs)
                                                       (logmean/logmean-mean-value xs n))
         :integral (if (m/== n 3)
                     (logmean/logmean3-integral xs)
                     (logmean/logmean-integral xs n max-iters tol)))))))

(defn wmean
  "Weighted mean"
  {:deprecated "Use `mean`"}
  (^double [vs] (mean vs))
  (^double [vs weights]
   (m// (v/dot vs weights) (sum weights))))

(defn population-variance
  "Calculates the population (biased) variance of `vs`.

  The mean of the squared deviations from the mean is divided by the number of observations `n`, unlike the sample variance [[variance]], which divides by `n-1` (Bessel's correction). Use this version when `vs` represents the entire population rather than a sample drawn from it.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `mu` (double, optional): Precomputed mean of `vs`. When omitted, the mean is computed from `vs`.

  Returns the population variance as a double.

  See also [[variance]], [[population-stddev]], [[population-wvariance]]."
  (^double [vs]
   (StatUtils/populationVariance (m/seq->double-array vs)))
  (^double [vs ^double mu]
   (StatUtils/populationVariance (m/seq->double-array vs) mu)))

(defn population-wvariance
  "Calculates the weighted population (biased) variance of `vs`.

  Each value is weighted by the corresponding entry in `freqs`, and the weighted mean of the squared deviations from the weighted mean is divided by the sum of the weights. This is the weighted analogue of [[population-variance]]; for the unbiased version, dividing by `(sum freqs) - 1` instead, see [[wvariance]].

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `freqs` (sequence of numbers): Weights (e.g. frequencies) corresponding to each value in `vs`, same length as `vs`.

  Returns the weighted population variance as a double.

  See also [[wvariance]], [[population-variance]], [[population-wstddev]]."
  ^double [vs freqs]
  (let [sw (sum freqs)
        mu (m// (v/dot vs freqs) sw)
        v (sum (map (fn [^double x ^double w]
                      (m/* w (m/sq (m/- x mu)))) vs freqs))]
    (m// v sw)))

(defn variance
  "Calculates the sample (unbiased) variance of `vs`.

  The sum of squared deviations from the mean is divided by `n-1` (Bessel's correction), where `n` is the number of observations. Use this version when `vs` is a sample drawn from a larger population; for the biased version dividing by `n`, see [[population-variance]].

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `mu` (double, optional): Precomputed mean of `vs`. When omitted, the mean is computed from `vs`.

  Returns the sample variance as a double.

  See also [[population-variance]], [[stddev]], [[wvariance]]."
  (^double [vs]
   (StatUtils/variance (m/seq->double-array vs)))
  (^double [vs ^double mu]
   (StatUtils/variance (m/seq->double-array vs) mu)))

(defn wvariance
  "Calculates the weighted sample (unbiased) variance of `vs`.

  Each value is weighted by the corresponding entry in `freqs`, and the weighted sum of squared deviations from the weighted mean is divided by `(sum freqs) - 1`. This is the weighted analogue of [[variance]]; for the biased version dividing by `(sum freqs)`, see [[population-wvariance]].

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `freqs` (sequence of numbers): Weights (e.g. frequencies) corresponding to each value in `vs`, same length as `vs`.

  Returns the weighted sample variance as a double.

  See also [[population-wvariance]], [[variance]], [[wstddev]]."
  ^double [vs freqs]
  (let [sw (sum freqs)
        mu (m// (v/dot vs freqs) sw)
        v (sum (map (fn [^double x ^double w]
                      (m/* w (m/sq (m/- x mu)))) vs freqs))]
    (m// v (m/dec sw))))

(defn population-stddev
  "Calculates the population (biased) standard deviation of `vs`.

  Computed as the square root of [[population-variance]], i.e. it divides the sum of squared deviations by `n` rather than `n-1`. Use this version when `vs` represents the entire population rather than a sample drawn from it.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `mu` (double, optional): Precomputed mean of `vs`. When omitted, the mean is computed from `vs`.

  Returns the population standard deviation as a double.

  See also [[stddev]], [[population-variance]], [[population-wstddev]]."
  (^double [vs]
   (m/sqrt (population-variance vs)))
  (^double [vs ^double mu]
   (m/sqrt (population-variance vs mu))))

(defn population-wstddev
  "Calculates the weighted population (biased) standard deviation of `vs`.

  Computed as the square root of [[population-wvariance]], i.e. each value is weighted by the corresponding entry in `freqs` and the sum of squared deviations is divided by the sum of the weights.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `weights` (sequence of numbers): Weights (e.g. frequencies) corresponding to each value in `vs`, same length as `vs`.

  Returns the weighted population standard deviation as a double.

  See also [[wstddev]], [[population-wvariance]], [[population-stddev]]."
  ^doubles [vs weights]
  (m/sqrt (population-wvariance vs weights)))

(defn stddev
  "Calculates the sample (unbiased) standard deviation of `vs`.

  Computed as the square root of [[variance]], i.e. it divides the sum of squared deviations by `n-1` (Bessel's correction). Use this version when `vs` is a sample drawn from a larger population.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `mu` (double, optional): Precomputed mean of `vs`. When omitted, the mean is computed from `vs`.

  Returns the sample standard deviation as a double.

  See also [[population-stddev]], [[variance]], [[wstddev]]."
  (^double [vs]
   (m/sqrt (variance vs)))
  (^double [vs ^double mu]
   (m/sqrt (variance vs mu))))

(defn wstddev
  "Calculates the weighted sample (unbiased) standard deviation of `vs`.

  Computed as the square root of [[wvariance]], i.e. each value is weighted by the corresponding entry in `freqs` and the sum of squared deviations is divided by `(sum freqs) - 1`.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `freqs` (sequence of numbers): Weights (e.g. frequencies) corresponding to each value in `vs`, same length as `vs`.

  Returns the weighted sample standard deviation as a double.

  See also [[population-wstddev]], [[wvariance]], [[stddev]]."
  ^doubles [vs freqs]
  (m/sqrt (wvariance vs freqs)))

(defn variation
  "Calculates the coefficient of variation (CV) for a sequence `vs`.

  The CV is a standardized measure of dispersion of a probability distribution or
  frequency distribution. It is defined as the ratio of the standard deviation
  to the mean:

  `CV = stddev(vs) / mean(vs)`

  This measure is unitless and allows for comparison of variability between
  datasets with different means or different units.

  Parameters:

  - `vs`: Sequence of numbers.

  Returns the calculated coefficient of variation as a double.

  Note: The CV is undefined if the mean is zero, and may be misleading if the
  mean is close to zero or if the data can take both positive and negative values.
  All values in `vs` should ideally be positive.

  See also [[stddev]], [[mean]]."
  ^double [vs]
  (let [vs (m/seq->double-array vs)]
    (m// (stddev vs)
         (mean vs))))

(defn median-absolute-deviation
  "Calculates the Median Absolute Deviation (MAD) of a sequence `vs`.

  MAD is a robust measure of the variability of a univariate sample of quantitative
  data. It is defined as the median of the absolute deviations from the data's median
  (or a specified center).

  `MAD = median(|X_i - median(X)|)`

  Parameters:

  - `vs`: Sequence of numbers.
  - `center-or-estimation-strategy` (optional): The central point from which to calculate deviations or estimation strategy.
    If `nil` or not provided, the [[median]] of `vs` is used as the center. If keyword, it's treated as estimation strategy for median.
  - `estimation-strategy` (optional, keyword): The estimation strategy to use for
    calculating the median(s). This applies to the calculation of the central
    value (if `center` is not provided) and to the final median of the absolute
    deviations. See [[median]] or [[quantile]] for available strategies (e.g.,
    `:legacy`, `:r1` through `:r9`).

  Returns the calculated MAD as a double.

  MAD is less sensitive to outliers than the standard deviation.

  See also [[mean-absolute-deviation]], [[stddev]], [[median]], [[quantile]]."
  (^double [vs] (median-absolute-deviation vs nil))
  (^double [vs center-or-estimation-strategy]
   (let [[es center] (if (keyword? center-or-estimation-strategy)
                       [center-or-estimation-strategy nil]
                       [nil center-or-estimation-strategy])
         m (double (or center (median vs es)))]
     (median (map (fn [^double x] (m/abs (m/- x m))) vs) es)))
  (^double [vs center estimation-strategy]
   (let [m (double (or center (median vs estimation-strategy)))]
     (median (map (fn [^double x] (m/abs (m/- x m))) vs) estimation-strategy))))

(def ^{:doc "Alias for [[median-absolute-deviation]]"}
  mad median-absolute-deviation)

(defn mean-absolute-deviation
  "Calculates the Mean Absolute Deviation of a sequence `vs`.

  MeanAD is a measure of the variability of a univariate sample of quantitative data.
  It is defined as the mean of the absolute deviations from a central point,
  typically the data's mean.

  `MeanAD = mean(|X_i - center|)`

  Parameters:

  - `vs`: Sequence of numbers.
  - `center` (optional, double): The central point from which to calculate deviations.
    If `nil` or not provided, the arithmetic [[mean]] of `vs` is used as the center.

  Returns the calculated Mean Absolute Deviation as a double.

  Unlike [[median-absolute-deviation]], which uses the median of absolute deviations
  from the median, the Mean Absolute Deviation uses the mean of absolute deviations
  from the mean (or specified center). This makes it more sensitive to outliers
  than [[median-absolute-deviation]] but less sensitive than the standard deviation.

  See also [[median-absolute-deviation]], [[stddev]], [[mean]]."
  (^double [vs] (mean-absolute-deviation vs nil))
  (^double [vs center]
   (let [m (double (or center (mean vs)))]
     (mean (map (fn [^double x] (m/abs (m/- x m))) vs)))))

(defn sem
  "Calculates the Standard Error of the Mean (SEM) for a sequence `vs`.

  The SEM estimates the standard deviation of the sample mean, providing an
  indication of how accurately the sample mean represents the population mean.
  It is calculated as:

  `SEM = stddev(vs) / sqrt(count(vs))`

  where `stddev(vs)` is the sample standard deviation and `count(vs)` is the
  sample size.

  Parameters:

  - `vs`: Sequence of numbers.

  Returns the calculated SEM as a double.

  A smaller SEM indicates that the sample mean is likely to be a more precise
  estimate of the population mean.

  See also [[stddev]], [[mean]]."
  ^double [vs]
  (m// (stddev vs)
       (m/sqrt (count vs))))

(defn stddev-extent
  "Calculates the mean of `vs` together with its `-/+` standard deviation extent.

  Parameters:

  - `vs` (sequence of numbers): Input data.

  Returns a 3-element vector `[(m/- mean stddev), (m/+ mean stddev), mean]`.

  See also [[mean]], [[stddev]], [[mad-extent]], [[sem-extent]]."
  [vs]
  (let [vs (m/seq->double-array vs)
        m (mean vs)
        s (stddev vs)]
    [(m/- m s) (m/+ m s) m]))

(defn mad-extent
  "Calculates the median of `vs` together with its `-/+` median absolute deviation (MAD) extent.

  Parameters:

  - `vs` (sequence of numbers): Input data.

  Returns a 3-element vector `[(m/- median mad), (m/+ median mad), median]`.

  See also [[median]], [[median-absolute-deviation]], [[stddev-extent]], [[sem-extent]]."
  [vs]
  (let [vs (m/seq->double-array vs)
        m (median vs)
        s (median-absolute-deviation vs)]
    [(m/- m s) (m/+ m s) m]))

(defn sem-extent
  "Calculates the mean of `vs` together with its `-/+` standard error of the mean (SEM) extent.

  Parameters:

  - `vs` (sequence of numbers): Input data.

  Returns a 3-element vector `[(m/- mean sem), (m/+ mean sem), mean]`.

  See also [[mean]], [[sem]], [[stddev-extent]], [[mad-extent]]."
  [vs]
  (let [vs (m/seq->double-array vs)
        m (mean vs)
        s (sem vs)]
    [(m/- m s) (m/+ m s) m]))

(defn percentile-extent
  "Calculates a pair of percentiles of `vs` together with its median.

  By default, computes the `p`-th and `(100-p)`-th percentiles, forming a symmetric interval around the median; two independent percentiles `p1` and `p2` can be given instead.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `p` (double): Lower percentile, in range `[0,100]`. The upper percentile is `100-p`. Defaults to `25.0`.
  - `p1`, `p2` (doubles): Two independent percentiles, in range `[0,100]`, used instead of `p`/`100-p`.
  - `estimation-strategy` (keyword): Percentile estimation method, one of the keys of [[estimation-strategies-list]]. Defaults to `:legacy`.

  Returns a 3-element vector `[(percentile vs p1), (percentile vs p2), (median vs)]`.

  See also [[percentile]], [[median]], [[quantile-extent]], [[pi-extent]]."
  ([vs] (percentile-extent vs 25.0))
  ([vs ^double p] (percentile-extent vs p (m/- 100.0 p)))
  ([vs p1 p2] (percentile-extent vs p1 p2 :legacy))
  ([vs ^double p1 ^double p2 estimation-strategy]
   (let [avs (m/seq->double-array vs)]
     [(percentile avs p1 estimation-strategy)
      (percentile avs p2 estimation-strategy)
      (median avs)])))

(defn quantile-extent
  "Calculates a pair of quantiles of `vs` together with its median.

  By default, computes the `q`-th and `(1.0-q)`-th quantiles, forming a symmetric interval around the median; two independent quantiles `q1` and `q2` can be given instead.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `q` (double): Lower quantile, in range `[0,1]`. The upper quantile is `1.0-q`. Defaults to `0.25`.
  - `q1`, `q2` (doubles): Two independent quantiles, in range `[0,1]`, used instead of `q`/`1.0-q`.
  - `estimation-strategy` (keyword): Quantile estimation method, one of the keys of [[estimation-strategies-list]]. Defaults to `:legacy`.

  Returns a 3-element vector `[(quantile vs q1), (quantile vs q2), (median vs)]`.

  See also [[quantile]], [[median]], [[percentile-extent]], [[pi-extent]]."
  ([vs] (quantile-extent vs 0.25))
  ([vs ^double q] (quantile-extent vs q (m/- 1.0 q)))
  ([vs q1 q2] (quantile-extent vs q1 q2 :legacy))
  ([vs ^double q1 ^double q2 estimation-strategy]
   (let [avs (m/seq->double-array vs)]
     [(quantile avs q1 estimation-strategy)
      (quantile avs q2 estimation-strategy)
      (median avs)])))

(defn pi
  "Calculates the Percentile Interval (PI), a symmetric quantile-based credible interval of `vs`.

  Given a target probability mass `size`, the interval spans the quantiles `(1-size)/2` and `1-(1-size)/2`, so that `size` proportion of the data lies within it. Unlike [[hpdi-extent]], this interval is not necessarily the narrowest one covering `size`, but it is symmetric with respect to probability mass on each tail.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `size` (double): Target probability content of the interval, in range `[0,1]`. Defaults to `0.5`.
  - `estimation-strategy` (keyword): Quantile estimation method, one of the keys of [[estimation-strategies-list]]. Defaults to `:legacy`.

  Returns a map of two entries, keyed by the lower and upper bound percentage mapping to the corresponding quantile values.

  See also [[pi-extent]] (vector form including the median), [[hpdi-extent]], [[quantile]]."
  ([vs] (pi vs 0.5))
  ([vs ^double size] (pi vs size :legacy))
  ([vs ^double size estimation-strategy]
   (let [a (m/* 0.5 (m/- 1.0 size))
         avs (m/seq->double-array vs)
         q1 (quantile avs a estimation-strategy)
         q2 (quantile avs (m/- 1.0 a) estimation-strategy)]
     {(m/* a 100.0) q1
      (m/* (m/- 1.0 a) 100.0) q2})))

(defn pi-extent
  "Calculates the Percentile Interval (PI) of `vs` together with its median.

  Given a target probability mass `size`, the interval spans the quantiles `(1-size)/2` and `1-(1-size)/2`. This is the vector-returning counterpart of [[pi]], following the same convention as [[percentile-extent]] and [[quantile-extent]].

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `size` (double): Target probability content of the interval, in range `[0,1]`. Defaults to `0.5`.
  - `estimation-strategy` (keyword): Quantile estimation method, one of the keys of [[estimation-strategies-list]]. Defaults to `:legacy`.

  Returns a 3-element vector `[lower-bound upper-bound median]`.

  See also [[pi]], [[hpdi-extent]], [[quantile-extent]]."
  ([vs] (pi-extent vs 0.5))
  ([vs ^double size] (pi-extent vs size :legacy))
  ([vs ^double size estimation-strategy]
   (let [a (m/* 0.5 (m/- 1.0 size))]
     (quantile-extent vs a (m/- 1.0 a) estimation-strategy))))

(defn hpdi-extent
  "Calculates the Highest Posterior Density Interval (HPDI) of `vs` together with its median.

  Unlike the symmetric [[pi-extent]], which fixes the probability mass on each tail, the HPDI is the narrowest interval that contains the target probability content `size` of the (sorted) samples. It is found by sliding a fixed-size window across the sorted data and picking the position with the smallest width.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `size` (double): Target probability content of the interval, in range `(0,1)`. Defaults to `0.95`.

  Returns a 3-element vector `[lower-bound upper-bound median]`.

  See also [[pi-extent]], [[quantile-extent]]."
  ([vs] (hpdi-extent vs 0.95))
  ([vs ^double size]
   (let [avs (m/seq->double-array vs)
         nsamp (alength avs)
         gap (m/constrain (m/round (m/* nsamp size)) 1 (m/dec nsamp))
         max-idx (m/- nsamp gap)]
     (java.util.Arrays/sort avs)
     (loop [idx (long 0)
            min-idx (long 0)
            mn Double/MAX_VALUE]
       (if (m/< idx max-idx)
         (let [diff (m/- (aget avs (m/+ idx gap))
                         (aget avs idx))]
           (if (m/< diff mn)
             (recur (m/inc idx) idx diff)
             (recur (m/inc idx) min-idx mn)))
         [(aget avs min-idx) (aget avs (m/+ min-idx gap)) (median avs)])))))

(defn iqr
  "Calculates the interquartile range (IQR) of `vs`.

  The IQR is the difference between the third quartile (75th percentile) and the first quartile (25th percentile), and is a robust measure of statistical dispersion that is insensitive to outliers.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `estimation-strategy` (keyword): Percentile estimation method, one of the keys of [[estimation-strategies-list]]. Defaults to `:legacy`.

  Returns the IQR as a double.

  See also [[percentile-extent]], [[inner-fence-extent]], [[outer-fence-extent]]."
  (^double [vs] (iqr vs :legacy))
  (^double [vs estimation-strategy]
   (let [[^double q1 ^double q3] (percentile-extent vs 25.0 75.0 estimation-strategy)]
     (m/- q3 q1))))

(defn adjacent-values
  "Calculates the lower and upper adjacent values (LAV and UAV) of `vs`, together with its median.

  LAV and UAV are the actual data values used as whisker endpoints in a Tukey box-and-whisker plot: rather than the inner-fence thresholds themselves (see [[inner-fence-extent]]), they are the most extreme observed values that still fall within those thresholds. Let `Q1` and `Q3` be the 25th and 75th percentiles and `IQR = Q3 - Q1`:

  - LAV is the smallest value in `vs` that is greater than or equal to the lower inner fence, `Q1 - 1.5*IQR`.
  - UAV is the largest value in `vs` that is less than or equal to the upper inner fence, `Q3 + 1.5*IQR`.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `estimation-strategy` (keyword): Percentile estimation method, one of the keys of [[estimation-strategies-list]]. Defaults to `:legacy`.
  - `q1`, `q3`, `m` (doubles): Precomputed first quartile, third quartile, and median, used instead of computing them from `vs`.

  Returns a 3-element vector `[LAV UAV median]`.

  See also [[inner-fence-extent]], [[outer-fence-extent]], [[iqr]], [[percentile-extent]]."
  ([vs]
   (adjacent-values vs :legacy))
  ([vs estimation-strategy]
   (let [avs (m/seq->double-array vs)
         [q1 m q3] (percentiles avs [25.0 50.0 75.0] estimation-strategy)]
     (adjacent-values avs q1 q3 m)))
  ([vs ^double q1 ^double q3 ^double m]
   (let [avs (m/seq->double-array vs)
         iqr (m/* 1.5 (m/- q3 q1))
         lav-thr (m/- q1 iqr)
         uav-thr (m/+ q3 iqr)]
     (java.util.Arrays/sort avs)
     [(first (filter (fn [^double v] (m/>= v lav-thr)) avs))
      (last (filter (fn [^double v] (m/<= v uav-thr)) avs))
      m])))

(defn inner-fence-extent
  "Calculates the lower and upper inner fence thresholds (LIF and UIF) of `vs`, together with its median.

  The inner fences are the classic Tukey box-and-whisker plot outlier thresholds, computed from `Q1` and `Q3`, the 25th and 75th percentiles, and `IQR = Q3 - Q1`:

  - LIF is the lower inner fence, `Q1 - 1.5*IQR`.
  - UIF is the upper inner fence, `Q3 + 1.5*IQR`.

  Unlike [[adjacent-values]], these are threshold values and not necessarily values present in `vs`. Values below LIF or above UIF are considered (mild) outliers, see [[outliers]] and [[remove-outliers]].

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `estimation-strategy` (keyword): Percentile estimation method, one of the keys of [[estimation-strategies-list]]. Defaults to `:legacy`.

  Returns a 3-element vector `[LIF UIF median]`.

  See also [[outer-fence-extent]], [[adjacent-values]], [[iqr]], [[percentile-extent]]."
  ([vs] (inner-fence-extent vs :legacy))
  ([vs estimation-strategy]
   (let [[^double q1 ^double m ^double q3] (percentiles vs [25.0 50.0 75.0] estimation-strategy)
         iqr+ (m/* 1.5 (m/- q3 q1))]
     [(m/- q1 iqr+) (m/+ q3 iqr+) m])))

(defn outer-fence-extent
  "Calculates the lower and upper outer fence thresholds (LOF and UOF) of `vs`, together with its median.

  The outer fences are wider outlier thresholds than the inner fences (see [[inner-fence-extent]]), computed from `Q1` and `Q3`, the 25th and 75th percentiles, and `IQR = Q3 - Q1`:

  - LOF is the lower outer fence, `Q1 - 3*IQR`.
  - UOF is the upper outer fence, `Q3 + 3*IQR`.

  Values beyond the outer fences are considered extreme (far) outliers, in contrast to the milder outliers flagged by the inner fences.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `estimation-strategy` (keyword): Percentile estimation method, one of the keys of [[estimation-strategies-list]]. Defaults to `:legacy`.

  Returns a 3-element vector `[LOF UOF median]`.

  See also [[inner-fence-extent]], [[adjacent-values]], [[iqr]], [[percentile-extent]]."
  ([vs] (outer-fence-extent vs :legacy))
  ([vs estimation-strategy]
   (let [[^double q1 ^double m ^double q3] (percentiles vs [25.0 50.0 75.0] estimation-strategy)
         iqr+ (m/* 3.0 (m/- q3 q1))]
     [(m/- q1 iqr+) (m/+ q3 iqr+) m])))

(defn span
  "Calculates the span (range width) of `vs`.

  The span is the difference between the maximum and the minimum values of the data, and gives the total width covered by the sample.

  Parameters:

  - `vs` (sequence of numbers): Input data.

  Returns the span as a double. Returns `0.0` when all values in `vs` are equal.

  See also [[extent]], [[maximum]], [[minimum]]."
  ^double [vs]
  (let [avs (m/seq->double-array vs)]
    (m/- (maximum avs) (minimum avs))))

(defn extent
  "Calculates the minimum and maximum values of `vs`, optionally together with the mean.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `mean?` (boolean): Whether to append the mean value to the result. Defaults to `true`.

  Returns a 2-element vector `[min max]` when `mean?` is `false`, or a 3-element vector `[min max mean]` when `mean?` is `true`.

  See also [[span]], [[maximum]], [[minimum]], [[mean]]."
  ([vs] (extent vs true))
  ([vs mean?]
   (let [fv (double (first vs))
         mm (reduce (fn [^Vec2 curr ^double v]
                      (Vec2. (m/min (.x curr) v) (m/max (.y curr) v))) (Vec2. fv fv) (rest vs))]
     (if mean? (conj mm (mean vs)) mm))))

(defn- process-outliers
  ([f vs] (process-outliers f vs :legacy))
  ([f vs estimation-strategy]
   (let [avs (m/seq->double-array vs)
         q1 (percentile avs 25.0 estimation-strategy)
         q3 (percentile avs 75.0 estimation-strategy)]
     (process-outliers f avs q1 q3)))
  ([f vs ^double q1 ^double q3]
   (let [iqr (m/* 1.5 (m/- q3 q1))
         lif-thr (m/- q1 iqr)
         uif-thr (m/+ q3 iqr)]
     (f (fn [^double v]
          (or (m/< v lif-thr)
              (m/> v uif-thr))) vs))))

(defn outliers
  "Finds outliers in `vs`, defined as values falling outside the inner fences.

  Let `Q1` and `Q3` be the 25th and 75th percentiles and `IQR = Q3 - Q1`. A value is considered an outlier when it is below the lower inner fence, `Q1 - 1.5*IQR`, or above the upper inner fence, `Q3 + 1.5*IQR`.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `estimation-strategy` (keyword): Percentile estimation method, one of the keys of [[estimation-strategies-list]]. Defaults to `:legacy`.
  - `q1`, `q3` (doubles): Precomputed first and third quartiles, used instead of computing them from `vs`.

  Returns a lazy sequence of the values from `vs` that lie outside the inner fences, in their original order.

  See also [[remove-outliers]], [[inner-fence-extent]], [[adjacent-values]], [[iqr]]."
  ([vs] (process-outliers filter vs))
  ([vs estimation-strategy] (process-outliers filter vs estimation-strategy))
  ([vs ^double q1 ^double q3] (process-outliers filter vs q1 q3)))

(defn remove-outliers
  "Removes outliers from `vs`, defined as values falling outside the inner fences.

  Let `Q1` and `Q3` be the 25th and 75th percentiles and `IQR = Q3 - Q1`. A value is considered an outlier when it is below the lower inner fence, `Q1 - 1.5*IQR`, or above the upper inner fence, `Q3 + 1.5*IQR`.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `estimation-strategy` (keyword): Percentile estimation method, one of the keys of [[estimation-strategies-list]]. Defaults to `:legacy`.
  - `q1`, `q3` (doubles): Precomputed first and third quartiles, used instead of computing them from `vs`.

  Returns a lazy sequence of the values from `vs` that lie within the inner fences, in their original order.

  See also [[outliers]], [[inner-fence-extent]], [[adjacent-values]], [[iqr]]."
  ([vs] (process-outliers remove vs))
  ([vs estimation-strategy] (process-outliers remove vs estimation-strategy))
  ([vs ^double q1 ^double q3] (process-outliers remove vs q1 q3)))

(defn wmodes
  "Returns the weighted mode(s) of a sequence `vs`.

  The mode is the value that appears most often in a dataset. This function
  generalizes the mode concept by considering weights associated with each value.
  A value's contribution to the mode calculation is proportional to its weight.

  Parameters:

  - `vs`: Sequence of data values. Can contain any data type (numbers, keywords, etc.).
  - `weights` (optional): Sequence of non-negative weights corresponding to `vs`.
    Must have the same count as `vs`. Defaults to a sequence of 1.0s if omitted,
    effectively calculating the unweighted modes.

  Returns a sequence containing all values that have the highest total weight.
  If there are ties (multiple values share the same maximum total weight), all
  tied values are included in the returned sequence. The order of modes in the
  returned sequence is not guaranteed.
  
  See also [[wmode]] (returns only one mode in case of ties) and [[modes]] (for unweighted numeric data)."
  ([vs] (wmodes vs (repeat 1.0)))
  ([vs weights]
   (let [all (->> (map vector vs weights)
                  (reduce (fn [b [v w]]
                            (update b v (fnil m/+ 0.0) w)) {})
                  (sort-by second m/>))
         best-score (-> all first second double)]
     (->> all
          (take-while (fn [p] (m/== best-score (double (second p)))))
          (map first)))))

(defn wmode
  "Returns the primary weighted mode of a sequence `vs`.

  The mode is the value that appears most often in a dataset. This function
  generalizes the mode concept by considering weights associated with each value.
  A value's contribution to the mode calculation is proportional to its weight.

  If multiple values share the same highest total weight (i.e., there are ties
  for the mode), this function returns only the first one encountered during
  processing. The specific mode returned in case of a tie is not guaranteed
  to be stable across different runs or environments. Use [[wmodes]] if you
  need all tied modes.

  Parameters:

  - `vs`: Sequence of data values. Can contain any data type (numbers, keywords, etc.).
  - `weights` (optional): Sequence of non-negative weights corresponding to `vs`.
    Must have the same count as `vs`. Defaults to a sequence of 1.0s if omitted,
    effectively calculating the unweighted mode.

  Returns a single value representing the mode (or one of the modes if ties exist).

  See also [[wmodes]] (returns all modes) and [[mode]] (for unweighted numeric data)."
  ([vs] (wmode vs (repeat 1.0)))
  ([vs weights]
   (first (wmodes vs weights))))

(declare histogram)

(defn modes
  "Find the values that appear most often in a dataset `vs`.

  Returns sequence with all most appearing values. For the default method
  (discrete data), modes are sorted in increasing order.

  For samples potentially drawn from a continuous distribution, simply finding the
  most frequent exact value might not be meaningful. Several estimation methods
  are provided via the `method` argument:

  * `:histogram`: Calculates the mode(s) based on the peak(s) of a histogram constructed
    from `vs`. Uses interpolation within the bin(s) with the highest frequency.
    Accepts options via `opts`, primarily `:bins` to control histogram
    construction (see [[histogram]]).
  * `:pearson`: Estimates the mode using Pearson's second skewness coefficient
    formula: `mode ≈ 3 * median - 2 * mean`. Accepts `:estimation-strategy`
    in `opts` for median calculation (see [[median]]). Returns a single estimated mode.
  * `:kde`: Estimates the mode(s) by finding the original data points in `vs`
    with the highest estimated probability density, based on Kernel Density
    Estimation (KDE). Accepts KDE options in `opts` like `:kernel`, `:bandwidth`,
    etc. (passed to `fastmath.kernel.density/kernel-density`).
  * `:default` (or when `method` is omitted): Finds the exact value(s) that occur
    most frequently in `vs`. Suitable for discrete data.

  The optional `opts` map provides method-specific configuration.

  See also [[mode]] (returns only the first mode) and [[wmodes]] (for weighted data)."
  ([vs method] (modes vs method {}))
  ([vs method opts]
   (let [avs (m/seq->double-array vs)]
     (case method
       :histogram (let [{:keys [bins ^double step]} (histogram avs (get opts :bins :rice))
                        ibins (vec (map-indexed #(conj %2 %1) bins))]
                    (->> ibins
                         (map (fn [[^double L ^long fm ^long id]]
                                (let [f1 (long (second (get ibins (m/dec id) [0 0])))
                                      f2 (long (second (get ibins (m/inc id) [0 0])))]
                                  [(m/+ L (m/* step (m// (m/- fm f1)
                                                         (m/- (m/* 2.0 fm) f1 f2)))) (m/- fm)])))
                         (sort-by second)
                         (map first)))
       :pearson (let [mu (mean avs)
                      m (median avs (:estimation-strategy opts))]
                  [(m/- (m/* 3.0 m) (m/* 2.0 mu))])
       :kde (let [kde (kd/kernel-density (get opts :kernel :gaussian) avs opts)]
              (->> (map (fn [^double v] [v (m/- (double (kde v)))]) vs)
                   (sort-by second)
                   (map first)))
       (seq (StatUtils/mode avs)))))
  ([vs] (seq (StatUtils/mode (m/seq->double-array vs)))))

(defn mode
  "Find the value that appears most often in a dataset `vs`.

  If multiple values share the same highest frequency (or estimated density/histogram peak),
  this function returns only the *first* one encountered during processing. The specific
  mode returned in case of a tie is not guaranteed to be stable. Use [[modes]] if
  you need all tied modes.

  For samples potentially drawn from a continuous distribution, several estimation
  methods are provided via the `method` argument:

  * `:histogram`: Calculates the mode based on the peak of a histogram constructed from `vs`.
    Uses interpolation within the bin with the highest frequency.
    Accepts options via `opts`, primarily `:bins` to control histogram
    construction (see [[histogram]]).
  * `:pearson`: Estimates the mode using Pearson's second skewness coefficient
    formula: `mode ≈ 3 * median - 2 * mean`. Accepts `:estimation-strategy`
    in `opts` for median calculation (see [[median]]).
  * `:kde`: Estimates the mode by finding the original data point in `vs`
    with the highest estimated probability density, based on Kernel Density
    Estimation (KDE). Accepts KDE options in `opts` like `:kernel`, `:bandwidth`,
    etc. (passed to `fastmath.kernel.density/kernel-density`).
  * `:default` (or when `method` is omitted): Finds the exact value that occurs
    most frequently in `vs`. Suitable for discrete data.

  The optional `opts` map provides method-specific configuration.

  See also [[modes]] (returns all modes) and [[wmode]] (for weighted data)."
  (^double [vs method] (mode vs method {}))
  (^double [vs method opts]
   (first (modes vs method opts)))
  (^double [vs]
   (let [m (StatUtils/mode (m/seq->double-array vs))]
     (aget ^doubles m 0))))

(defn moment
  "Calculates the statistical moment of `vs` of the given `order` (default: `2.0`).

  This is a general building block used to compute variance-like, skewness-like, and kurtosis-like statistics by combining a few switches: whether deviations are taken as absolute values, around which center, whether the result is averaged, and whether it is normalized by the standard deviation. For example, `order=2.0` with default options gives the population variance, `order=3.0` with `:normalize?` set to `true` gives a skewness-like statistic, and `order=4.0` with `:normalize?` set to `true` gives a kurtosis-like statistic.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `order` (double): Power to raise the (optionally absolute) deviations to. Defaults to `2.0`.
  - `opts` (map, optional):
    - `:absolute?` (boolean) - when `true`, deviations are taken as absolute values before raising to `order`. Defaults to `false`.
    - `:mean?` (boolean) - when `true`, averages the deviations (a proper moment); when `false`, returns their plain sum. Defaults to `true`.
    - `:center` (double) - the center to compute deviations from. Defaults to `nil`, meaning the mean of `vs` is used.
    - `:normalize?` (boolean) - when `true`, divides the result by the standard deviation of `vs` raised to `order`. Defaults to `false`.

  Returns the computed moment as a double.

  See also [[variance]], [[skewness]], [[kurtosis]]."
  (^double [vs] (moment vs 2.0 nil))
  (^double [vs ^double order] (moment vs order nil))
  (^double [vs ^double order {:keys [absolute? center mean? normalize?]
                              :or {mean? true}}]
   (let [in (m/seq->double-array vs)
         cin (alength in)
         out (double-array cin)
         nf (if normalize? (m/pow (variance in) (m/* 0.5 order)) 1.0)
         center (double (or center (mean in)))
         f (cond
             (m/one? order) m/identity-double
             (m/== order 2.0) m/sq
             (m/== order 3.0) m/cb
             (m/== order 4.0) (fn ^double [^double diff] (m/sq (m/sq diff)))
             :else (fn ^double [^double diff] (m/pow diff order)))
         a (if absolute? m/abs m/identity-double)]
     (loop [idx (long 0)]
       (when (m/< idx cin)
         (aset out idx (double (f (a (m/- (aget in idx) center)))))
         (recur (m/inc idx))))
     (m// (if mean? (mean out) (sum out)) nf))))

(def ^{:deprecated "Use [[moment]] function"} second-moment moment)

;;

(defn l-moment
  "Calculates the L-moment, TL-moment (trimmed L-moment), or (T)L-moment ratio of `vs` for a given `order`.

  L-moments are summary statistics for describing the shape of a probability distribution, computed as linear combinations of the expected order statistics of the sample. Unlike conventional moments, they exist whenever the mean of the distribution exists, are more robust to outliers, and suffer less from bias in small samples. Trimming (`:s`, `:t`) further excludes extreme order statistics, yielding TL-moments, which remain well-defined even for distributions without a finite mean.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `order` (long): The order of the L-moment. Order `0` always returns `1.0`.
  - `opts` (map, optional):
    - `:s` (long) - number of smallest (left) values trimmed from the sample. Defaults to `0`.
    - `:t` (long) - number of largest (right) values trimmed from the sample. Defaults to `0`.
    - `:sorted?` (boolean) - set to `true` when `vs` is already sorted ascending, to skip re-sorting. Defaults to `false`.
    - `:ratio?` (boolean) - when `true`, returns the (T)L-moment ratio, i.e. the L-moment of `order` normalized by the second L-moment, instead of the raw L-moment. Defaults to `false`.

  Returns the computed (T)L-moment or (T)L-moment ratio as a double.

  See also [[moment]], [[variance]], [[skewness]]."
  (^double [vs ^long order] (l-moment vs order nil))
  (^double [vs ^long order {:keys [^long s ^long t sorted? ratio?]
                            :or {s 0 t 0}
                            :as opts}]
   (if (m/zero? order)
     1.0
     (let [^doubles svs (m/seq->double-array (if sorted? vs (sort vs)))]
       (if ratio?
         (let [nopts (assoc opts :sorted? true :ratio? false)
               l2 (l-moment svs 2 nopts)]
           (m// (l-moment svs order nopts) l2))
         (let [r- (m/long-dec order)
               s+ (m/inc s)
               n (alength svs)
               n-t+ (m/- n t -1)]
           (m// (double (reduce (fn [^double b1 ^long k]
                                  (let [c1 (m/long-sub (m/long-add r- s) k)
                                        c2 (m/long-add t k)]
                                    (-> (if (m/even? k) (m/combinations r- k) (m/- (m/combinations r- k)))
                                        (m/* (double (reduce (fn [^double b2 ^long j]
                                                               (let [j- (m/long-dec j)]
                                                                 (-> (m/* (m/combinations j- c1)
                                                                          (m/combinations (m/long-sub n j) c2)
                                                                          (Array/aget svs j-))
                                                                     (m/+ b2)))) 0.0 (range s+ n-t+))))
                                        (m/+ b1))))
                                0.0 (range order)))
                (m/* order (m/combinations n (m/long-add order s t))))))))))

(defn l-variation
  "Calculates the coefficient of L-variation (L-CV) of `vs`.

  This is a robust analogue of the ordinary coefficient of variation [[variation]]: instead of the ratio of the standard deviation to the mean, it is the ratio of the second L-moment (a measure of scale, see [[l-moment]]) to the mean.

  Parameters:

  - `vs` (sequence of numbers): Input data.

  Returns the L-CV as a double.

  See also [[l-moment]], [[variation]], [[mean]]."
  ^double [vs]
  (let [svs (m/seq->double-array vs)]
    (java.util.Arrays/sort svs)
    (m// (l-moment svs 2 {:sorted? true})
         (mean svs))))

;;; expectile

(defn- expectile-map
  [vs ^double tau ^double tau- ^double t]
  (map (fn [^double v]
         (m/* (if (m/<= v t) tau- tau) (m/- t v))) vs))

(defn- expectile-target
  ([vs weights ^double tau ^double tau-]
   (fn ^double [^double t]
     (let [m (expectile-map vs tau tau- t)]
       (if weights (mean m weights) (mean m)))))  )

(defn expectile
  "Calculate the tau-th expectile of a sequence `vs`.

  Expectiles are related to quantiles but are determined by minimizing an
  asymmetrically weighted sum of squared differences, rather than absolute
  differences. The `tau` parameter controls the asymmetry.

  A key property is that the expectile for `tau = 0.5` is equal to the [[mean]].

  The calculation involves finding the value `t` such that the weighted sum
  of `w_i * (v_i - t)` is zero, where the effective weights depend on `tau` and whether
  `v_i` is above or below `t`.

  Parameters:

  - `vs`: Sequence of data values.
  - `weights` (optional): Sequence of corresponding non-negative weights.
    Must have the same count as `vs`. If omitted, calculates the unweighted expectile.
  - `tau`: The expectile level, a value between 0.0 and 1.0 (inclusive).

  Returns the calculated expectile as a double.

  See also [[quantile]], [[mean]], [[median]]."
  (^double [vs ^double tau] (expectile vs nil tau))
  (^double [vs weights ^double tau]
   (let [avg (if weights (mean vs weights) (mean vs))]
     (if (m/== tau 0.5)
       avg
       (let [[^double x0 ^double x1] (if (m/> tau 0.5)
                                       [avg (maximum vs)]
                                       [(minimum vs) avg])]
         (if (m/== x0 x1)
           x0
           (solver/find-root (expectile-target vs weights tau (m/- 1.0 tau)) x0 x1)))))))

(defn winsor
  "Winsorizes `vs`, clamping extreme values to a pair of quantile-based bounds.

  Values below the lower bound are replaced by the lower bound, and values above the upper bound are replaced by the upper bound, so the size of the sequence is unchanged. `NaN` entries are always replaced by a given replacement value (the median, unless given explicitly).

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `quantile` (double): Trimming proportion applied to each tail, e.g. `0.2` clamps the bottom 20% and top 20% of the data to the 20th and 80th percentile values respectively. Defaults to `0.2`.
  - `estimation-strategy` (keyword): Quantile estimation method, one of the keys of [[estimation-strategies-list]]. Defaults to `:legacy`.
  - `low`, `high` (doubles): Explicit lower and upper bounds, used instead of computing them from `quantile`. Order does not matter, the smaller value is always used as the lower bound.
  - `nan` (double): Replacement value for `NaN` entries in `vs`.

  Returns a lazy sequence of the same length as `vs`, with out-of-bound values clamped and `NaN` values replaced.

  See also [[trim]], [[remove-outliers]], [[quantile]]."
  ([vs] (winsor vs 0.2))
  ([vs quantile] (winsor vs quantile :legacy))
  ([vs ^double quantile estimation-strategy]
   (let [[qlow qmid qhigh] (quantiles (remove m/nan? vs)
                                      [quantile 0.5 (m/- 1.0 quantile)] estimation-strategy)]
     (winsor vs qlow qhigh qmid)))
  ([vs ^double low ^double high nan]
   (let [[^double low ^double high] (if (m/< low high) [low high] [high low])]
     (map (fn [^double v]
            (if (m/nan? v)
              nan
              (m/constrain v low high))) vs))))

(defn trim
  "Trims `vs`, discarding values falling outside a pair of quantile-based bounds.

  Unlike [[winsor]], which clamps out-of-bound values, this function removes them entirely, so the resulting sequence may be shorter than `vs`. `NaN` entries are always kept but replaced by a given replacement value (the median, unless given explicitly).

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `quantile` (double): Trimming proportion applied to each tail, e.g. `0.2` discards the bottom 20% and top 20% of the data, below the 20th and above the 80th percentile. Defaults to `0.2`.
  - `estimation-strategy` (keyword): Quantile estimation method, one of the keys of [[estimation-strategies-list]]. Defaults to `:legacy`.
  - `low`, `high` (doubles): Explicit lower and upper bounds, used instead of computing them from `quantile`. Order does not matter, the smaller value is always used as the lower bound.
  - `nan` (double): Replacement value for `NaN` entries in `vs`.

  Returns a lazy sequence containing only the values from `vs` within `[low,high]`, with `NaN` values replaced.

  See also [[trim-lower]], [[trim-upper]], [[winsor]], [[remove-outliers]], [[quantile]]."
  ([vs] (trim vs 0.2))
  ([vs quantile] (trim vs quantile :legacy))
  ([vs ^double quantile estimation-strategy]
   (let [[qlow qmid qhigh] (quantiles (remove m/nan? vs)
                                      [quantile 0.5 (m/- 1.0 quantile)] estimation-strategy)]
     (trim vs qlow qhigh qmid)))
  ([vs ^double low ^double high nan]
   (let [[^double low ^double high] (if (m/< low high) [low high] [high low])]
     (->> vs
          (filter (fn [^double v]
                    (or (m/nan? v)
                        (m/<= low v high))))
          (map (fn [^double v] (if (m/nan? v) nan v)))))))

(defn trim-lower
  "Trims `vs`, discarding values below a given quantile.

  This is a one-sided version of [[trim]]: only values below the cutoff are discarded, values above it are all kept.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `quantile` (double): Trimming proportion, e.g. `0.2` discards values below the 20th percentile. Defaults to `0.2`.
  - `estimation-strategy` (keyword): Quantile estimation method, one of the keys of [[estimation-strategies-list]]. Defaults to `:legacy`.

  `NaN` entries are always kept but replaced by the median of `vs`.

  Returns a lazy sequence containing only the values from `vs` at or above the cutoff quantile, with `NaN` values replaced.

  See also [[trim-upper]], [[trim]], [[winsor]]."
  ([vs] (trim-lower vs 0.2))
  ([vs quantile] (trim-lower vs quantile :legacy))
  ([vs ^double quantile estimation-strategy]
   (let [[q qmid] (quantiles (remove m/nan? vs) [quantile 0.5] estimation-strategy)]
     (trim vs q ##Inf qmid))))

(defn trim-upper
  "Trims `vs`, discarding values above a given quantile.

  This is a one-sided version of [[trim]]: only values above the cutoff are discarded, values below it are all kept.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `quantile` (double): Trimming proportion, e.g. `0.2` discards values above the 80th percentile (`1.0-quantile`). Defaults to `0.2`.
  - `estimation-strategy` (keyword): Quantile estimation method, one of the keys of [[estimation-strategies-list]]. Defaults to `:legacy`.

  `NaN` entries are always kept but replaced by the median of `vs`.

  Returns a lazy sequence containing only the values from `vs` at or below the cutoff quantile, with `NaN` values replaced.

  See also [[trim-lower]], [[trim]], [[winsor]]."
  ([vs] (trim-upper vs 0.2))
  ([vs quantile] (trim-upper vs quantile :legacy))
  ([vs ^double quantile estimation-strategy]
   (let [[q qmid] (quantiles (remove m/nan? vs) [quantile 0.5] estimation-strategy)]
     (trim vs ##-Inf q qmid))))

;; On More Robust Estimation of Skewness and Kurtosis: Simulation and Application to the S&P500 Index

(defn- yule-skewness
  ^double [vs ^double u]
  (let [[^double q1 ^double q2 ^double q3] (quantiles vs [u 0.5 (m/- 1.0 u)])]
    (m// (m/+ q3 (m/* -2.0 q2) q1)
         (m/- q3 q1))))

(defn- bowley-skewness
  ^double [vs]
  (yule-skewness vs 0.25))

(defn- hogg-skewness
  ^double [vs]
  (let [m25 (mean (trim vs 0.25))
        u005 (mean (trim-lower vs 0.95))
        l005 (mean (trim-upper vs 0.05))]
    (m// (m/- u005 m25) (m/- m25 l005))))

(defn skewness
  "Calculates the skewness of `vs`, a measure of the asymmetry of a probability distribution about its mean.

  Several skewness definitions are supported through the `typ` argument, ranging from moment-based estimators to robust, quantile- or trimmed-mean-based estimators.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `typ` (keyword or sequence): Skewness measure to calculate. Defaults to `:G1`. One of:
    - `:G1` (default): sample skewness based on the third standardized moment, adjusted for sample size bias, as implemented by Apache Commons Math `Skewness`.
    - `:g1` or `:pearson`: Pearson's moment coefficient of skewness, a bias-adjusted version of the third standardized moment. Expected value `0` for symmetric distributions.
    - `:b1`: sample skewness coefficient (b1), related to `:g1`.
    - `:B1` or `:yule`: Yule's coefficient (robust), based on quantiles. The sequence form `[:B1 u]` or `[:yule u]` sets the quantile `u` used, default `0.25`.
    - `:B3`: robust measure comparing the mean and median relative to the mean absolute deviation around the median.
    - `:skew`: an adjusted skewness definition sometimes used in bootstrap (BCa) calculations.
    - `:mode`: Pearson's second skewness coefficient, `(mean - mode) / stddev`. The sequence form `[:mode method opts]` sets the mode estimation `method` and its `opts`, see [[mode]].
    - `:median`: robust measure, `3 * (mean - median) / stddev`.
    - `:bowley`: Bowley's coefficient (robust), also known as the Yule-Bowley coefficient, based on quartiles `Q1`, `Q2`, `Q3`: `(Q3 + Q1 - 2*Q2) / (Q3 - Q1)`.
    - `:hogg`: Hogg's robust measure, based on the ratio of differences between trimmed means.
    - `:l-skewness`: L-skewness (τ₃), the ratio of the third L-moment (λ₃) to the second L-moment (λ₂, L-scale), a robust measure of asymmetry. Calculated directly using [[l-moment]] with `:ratio?` set to `true`. Expected value `0` for symmetric distributions.

  Positive values generally indicate a distribution skewed to the right (a longer tail on the right), negative values indicate a distribution skewed to the left (a longer tail on the left), and values near `0` suggest relative symmetry.

  Returns the calculated skewness as a double.

  See also [[skewness-test]], [[kurtosis]], [[normality-test]], [[jarque-bera-test]], [[l-moment]], [[moment]]."
  (^double [vs] (skewness vs :G1))
  (^double [vs typ]
   (let [vs (m/seq->double-array vs)]
     (if (sequential? typ)
       (let [[typ a b] typ]
         (case typ
           :mode (m// (m/- (mean vs) (mode vs a b)) (stddev vs))
           (:B1 :yule) (yule-skewness vs a)))
       (case typ
         :mode (m// (m/- (mean vs) (mode vs)) (stddev vs))
         :median (m// (m/* 3.0 (m/- (mean vs) (median vs))) (stddev vs))
         :bowley (bowley-skewness vs)
         :hogg (hogg-skewness vs)
         (:B1 :yule) (yule-skewness vs 0.25)
         :B3 (let [v (median vs)]
               (m// (m/- (mean vs) v)
                    (moment vs 1.0 {:absolute? true :center v})))
         :l-skewness (l-moment vs 3 {:ratio? true})
         (let [^Skewness k (Skewness.)
               n (alength vs)
               v (.evaluate k vs)]
           (case typ
             :b1 (m/* v (m// (m/* (m/- n 2.0) (m/dec n)) (m/* n n)))
             (:pearson :g1) (m/* v (m// (m/- n 2.0) (m/sqrt (m/* n (m/dec n)))))
             :skew (m/* v (m// (m/- n 2.0) (m/* n (m/sqrt (m/dec n))))) ;; BCa based, g1/sqrt(n)
             v)))))))

;; centered
(defn- moors-kurtosis
  ^double [vs]
  (let [[^double e1 ^double e2 ^double e3
         ^double e5 ^double e6 ^double e7] (quantiles vs [0.125 0.25 0.375 0.625 0.75 0.875])]
    (m/- (m// (m/+ (m/- e7 e5) (m/- e3 e1))
            (m/- e6 e2)) 1.233)))

;; centered
(defn- crow-kurtosis
  (^double [vs] (crow-kurtosis vs 0.025 0.25))
  (^double [vs ^double alpha ^double beta]
   (let [[^double a1 ^double a2 ^double b1 ^double b2] (quantiles vs [alpha (m/- 1.0 alpha)
                                                                      beta (m/- 1.0 beta)])
         c (m// (double (r/icdf r/default-normal alpha))
                (double (r/icdf r/default-normal beta)))]
     (m/- (m// (m/- a2 a1) (m/- b2 b1)) c))))

;; centered
(defn- hogg-kurtosis
  (^double [vs] (hogg-kurtosis vs 0.05 0.5))
  (^double [vs ^double alpha ^double beta]
   (let [ua (mean (trim-lower vs (m/- 1.0 alpha)))
         ub (mean (trim-lower vs (m/- 1.0 beta)))
         la (mean (trim-upper vs alpha))
         lb (mean (trim-upper vs beta))]
     (m/- (m// (m/- ua la) (m/- ub lb)) 2.585))))


;; https://aakinshin.net/posts/misleading-kurtosis/
(defn kurtosis
  "Calculates the kurtosis of `vs`, a measure of the tailedness or peakedness of a distribution compared to a normal distribution.

  Several kurtosis definitions are supported through the `typ` argument, ranging from moment-based estimators to robust, quantile- or trimmed-mean-based estimators. Different types use different algorithms and have different expected values under normality (`0` or `3`).

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `typ` (keyword or sequence): Kurtosis measure to calculate. Defaults to `:G2`. One of:
    - `:G2` (default): sample excess kurtosis based on the fourth standardized moment, as implemented by Apache Commons Math `Kurtosis`. Expected value `0` for a normal distribution.
    - `:g2` or `:excess`: an alternative sample excess kurtosis, obtained from `:G2` by inverting its finite-sample bias correction. Expected value `0` for a normal distribution.
    - `:kurt`: the classical kurtosis definition where normal is `3`, calculated as `:g2 + 3`.
    - `:b2`: another excess kurtosis variant, obtained by applying an additional finite-sample correction factor to `:kurt`. Expected value `0` for a normal distribution.
    - `:geary`: Geary's g, a robust measure calculated as `mean-absolute-deviation / population-stddev`. Expected value for normal is `sqrt(2/pi) ≈ 0.798`, lower values indicate leptokurtosis.
    - `:moors`: Moors' robust kurtosis measure based on octiles, centered so that the expected value for normal is `0`.
    - `:crow`: Crow-Siddiqui robust kurtosis measure based on quantiles, centered so that the expected value for normal is `0`. The sequence form `[:crow alpha beta]` sets the quantile parameters `alpha` and `beta`.
    - `:hogg`: Hogg's robust kurtosis measure based on trimmed means, centered so that the expected value for normal is `0`. The sequence form `[:hogg alpha beta]` sets the trimming parameters `alpha` and `beta`.
    - `:l-kurtosis`: L-kurtosis (τ₄), the ratio of the fourth L-moment (λ₄) to the second L-moment (λ₂, L-scale), a robust measure. Calculated directly using [[l-moment]] with `:ratio?` set to `true`. Expected value for a normal distribution is approximately `0.1226`.

  For the excess-kurtosis variants (`:G2`, `:g2`, `:excess`, `:b2` and the centered robust measures), positive values indicate a leptokurtic distribution (heavier tails, more peaked than normal), negative values indicate a platykurtic distribution (lighter tails, flatter than normal), and values near `0` suggest kurtosis similar to a normal distribution.

  Returns the calculated kurtosis as a double.

  See also [[kurtosis-test]], [[bonett-seier-test]], [[skewness]], [[normality-test]], [[jarque-bera-test]], [[l-moment]], [[moment]]."
  (^double [vs] (kurtosis vs :G2))      ; Default to :G2 as per code
  (^double [vs typ]
   (let [vs (m/seq->double-array vs)
         n (alength vs)]
     (if (sequential? typ)
       (let [[typ ^double a ^double b] typ]
         (case typ
           :crow (crow-kurtosis vs a b)
           :hogg (hogg-kurtosis vs a b)))
       (case typ
         :geary (m// (mean-absolute-deviation vs)
                     (population-stddev vs))
         :moors (moors-kurtosis vs)
         :crow (crow-kurtosis vs)
         :hogg (hogg-kurtosis vs)
         ;; Use l-moment directly with :ratio? true for L-kurtosis (λ₄ / λ₂)
         :l-kurtosis (l-moment vs 4 {:ratio? true})
         ;; Default case for moment-based kurtosis
         (let [^Kurtosis k (Kurtosis.)
               v (.evaluate k vs)]
           (case typ
             (:excess :g2) (m// (m/- (m// (m/* v (m/- n 2) (m/- n 3)) (m/dec n)) 6.0) 
                              (m/inc n))
             :kurt (m/+ 3.0 (m// (m/- (m// (m/* v (m/- n 2) (m/- n 3)) (m/dec n)) 6.0) ; g2 + 3
                             (m/inc n)))
             :b2 (m/- (m/* (m/+ 3.0 (m// (m/- (m// (m/* v (m/- n 2) (m/- n 3)) (m/dec n)) 6.0) ; g2 + 3
                                 (m/inc n)))
                       (m/sq (m/- 1.0 (m// 1.0 n)))) 3.0)
             v    ; Default is :G2 (sample kurtosis from Commons Math)
             )))))))

(defn ci
  "Calculates the Student's t-distribution based confidence interval for the mean of `vs`, together with the mean itself.

  The interval is built as `mean ± critical-value * stddev / sqrt(n)`, where the critical value is the `1-alpha/2` quantile of the t-distribution with `n-1` degrees of freedom (`n` being the sample size). This is the standard parametric confidence interval for the mean, assuming approximately normal data.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `alpha` (double): Significance level, in range `(0,1)`. Defaults to `0.05`, giving a 95% confidence interval.

  Returns a 3-element vector `[lower-bound upper-bound mean]`.

  See also [[bootstrap-ci]], [[mean]], [[stddev]], [[sem-extent]]."
  ([vs] (ci vs 0.05))
  ([vs ^double alpha]
   (let [vsa (m/seq->double-array vs)
         cnt (count vs)
         dist (r/distribution :t {:degrees-of-freedom (m/dec cnt)})
         crit-val (double (r/icdf dist (m/- 1.0 (m/* 0.5 alpha))))
         mean-ci (m// (m/* crit-val (stddev vsa)) (m/sqrt cnt))
         mn (mean vsa)]
     [(m/- mn mean-ci) (m/+ mn mean-ci) mn])))

;; https://ocw.mit.edu/courses/mathematics/18-05-introduction-to-probability-and-statistics-spring-2014/readings/MIT18_05S14_Reading24.pdf
(defn bootstrap-ci
  "Bootstrap method to calculate confidence interval.

  Alpha defaults to 0.98, samples to 1000.
  Last parameter is statistical function used to measure, default: [[mean]].

  Returns ci and statistical function value."
  {:deprecated "Please use fastmath.stats.boostrap/ci-basic instead"}
  ([vs] (bootstrap-ci vs 0.98))
  ([vs alpha] (bootstrap-ci vs alpha 1000))
  ([vs alpha samples] (bootstrap-ci vs alpha samples mean))
  ([vs ^double alpha ^long samples stat-fn]
   (let [vsa (m/seq->double-array vs)
         cnt (count vs)
         dist (r/distribution :enumerated-real {:data vsa})
         m (double (stat-fn vsa))
         deltas (m/seq->double-array (repeatedly samples #(m/- (mean (r/->seq dist cnt)) m)))
         q1 (quantile deltas alpha)
         q2 (quantile deltas (m/- 1.0 alpha))]
     [(m/- m q1) (m/- m q2) m])))

(defn bootstrap
  {:doc "Generate set of samples of given size from provided data.

  Default `samples` is 200, number of `size` defaults to sample size."
   :deprecated "Please use fastmath.stats.bootstrap/bootstrap instead"}
  ([vs] (bootstrap vs 200))
  ([vs samples] (bootstrap vs samples (count vs)))
  ([vs samples size]
   (let [dist (r/distribution :enumerated-real {:data vs})]
     (repeatedly samples #(r/->seq dist size)))))

(defn stats-map
  "Calculates a comprehensive set of descriptive statistics for `vs` and returns them as a map, giving a quick overview of the data's central tendency, dispersion, shape, and potential outliers.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `estimation-strategy` (keyword): Quantile estimation method used for the median, quartiles, and fence-related statistics, one of the keys of [[estimation-strategies-list]]. Defaults to `:legacy`.

  Returns a map with the following keys:

  - `:Size`: number of data points, `count`.
  - `:Min`: minimum value, see [[minimum]].
  - `:Max`: maximum value, see [[maximum]].
  - `:Range`: `Max - Min`.
  - `:Mean`: arithmetic average, see [[mean]].
  - `:Median`: middle value, see [[median]].
  - `:Mode`: most frequent value, see [[mode]] (default method).
  - `:Q1`: first quartile (25th percentile), see [[percentile]].
  - `:Q3`: third quartile (75th percentile), see [[percentile]].
  - `:Total`: sum of all values, see [[sum]].
  - `:SD`: sample standard deviation, see [[stddev]].
  - `:Variance`: sample variance, `SD^2`, see [[variance]].
  - `:MAD`: median absolute deviation, see [[median-absolute-deviation]].
  - `:SEM`: standard error of the mean, see [[sem]].
  - `:LAV`: lower adjacent value, the smallest value within the lower inner fence, see [[adjacent-values]].
  - `:UAV`: upper adjacent value, the largest value within the upper inner fence, see [[adjacent-values]].
  - `:IQR`: interquartile range, `Q3 - Q1`.
  - `:LOF`: lower outer fence, `Q1 - 3*IQR`, see [[outer-fence-extent]].
  - `:UOF`: upper outer fence, `Q3 + 3*IQR`, see [[outer-fence-extent]].
  - `:LIF`: lower inner fence, `Q1 - 1.5*IQR`, see [[inner-fence-extent]].
  - `:UIF`: upper inner fence, `Q3 + 1.5*IQR`, see [[inner-fence-extent]].
  - `:Outliers`: sequence of data points falling outside the inner fences, see [[outliers]].
  - `:Kurtosis`: measure of tailedness/peakedness, see [[kurtosis]] (default `:G2` type).
  - `:Skewness`: measure of asymmetry, see [[skewness]] (default `:G1` type).

  See also [[extent]], [[percentile-extent]], [[adjacent-values]], [[inner-fence-extent]], [[outer-fence-extent]], [[outliers]]."
  ([vs] (stats-map vs :legacy))
  ([vs estimation-strategy]
   (let [avs (m/seq->double-array vs)
         sz (alength avs)
         mn (Array/min avs)
         mx (Array/max avs)
         sm (Array/sum avs)
         u (m// sm sz)
         mdn (median avs)
         q1 (percentile avs 25.0 estimation-strategy)
         q3 (percentile avs 75.0 estimation-strategy)
         iqr (m/- q3 q1)
         sd (stddev avs)
         mad (median-absolute-deviation avs)
         [lav uav] (adjacent-values avs q1 q3 mdn)]
     {:Size sz
      :Min mn
      :Max mx
      :Range (m/- mx mn)
      :Mean u
      :Median mdn
      :Mode (mode avs)
      :Q1 q1
      :Q3 q3
      :Total sm
      :SD sd
      :Variance (m/* sd sd)
      :MAD mad
      :SEM (m// sd (m/sqrt sz))
      :LAV lav
      :UAV uav
      :IQR iqr
      :LOF (m/- q1 (m/* 3.0 iqr))
      :UOF (m/+ q3 (m/* 3.0 iqr))
      :LIF (m/- q1 (m/* 1.5 iqr))
      :UIF (m/+ q3 (m/* 1.5 iqr))
      :Outliers (outliers avs q1 q3)
      :Kurtosis (kurtosis avs)
      :Skewness (skewness avs)})))

;;

(defn standardize
  "Standardizes `vs` (z-scores) so that the result has a mean of `0` and a standard deviation of `1`.

  Each value `x` in `vs` is transformed to `(x - mean) / stddev`, using the sample mean and sample standard deviation of `vs`.

  Parameters:

  - `vs` (sequence of numbers): Input data.

  Returns a lazy sequence of standardized values, of the same length as `vs`.

  Not robust to outliers; use [[robust-standardize]] for a median/MAD-based alternative.

  See also [[robust-standardize]], [[mean]], [[stddev]]."
  [vs]
  (seq (StatUtils/normalize (m/seq->double-array vs))))

(defn robust-standardize
  "Standardizes `vs` using robust, outlier-resistant statistics instead of mean and standard deviation.

  Each value `x` in `vs` is transformed to `(x - median) / scale`. By default `scale` is the median absolute deviation (MAD, see [[median-absolute-deviation]]), giving a result centered at `0` with `MAD = 1`. This is a robust counterpart of [[standardize]], less sensitive to outliers.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `q` (double, optional): When provided, `scale` is instead the quantile difference `Q(1-q) - Q(q)`, in range `(0,0.5)`. For example, `0.25` uses the interquartile range (IQR) as the scale.

  Returns a lazy sequence of standardized values, of the same length as `vs`.

  See also [[standardize]], [[median]], [[median-absolute-deviation]], [[iqr]], [[quantiles]]."
  ([vs]
   (let [avs (m/seq->double-array vs)
         mad (median-absolute-deviation avs)
         md (median avs)]
     (map (fn [^double x] (m// (m/- x md) mad)) vs)))
  ([vs ^double q]
   (let [[^double q1 ^double md ^double q2] (quantiles vs [q 0.5 (m/- 1.0 q)])
         diff (m/abs (m/- q2 q1))]
     (map (fn [^double x] (m// (m/- x md) diff)) vs))))

(defn demean
  "Subtracts the mean of `vs` from each of its values, centering the data around `0`.

  Parameters:

  - `vs` (sequence of numbers): Input data.

  Returns a lazy sequence of the same length as `vs`, with mean `0`.

  See also [[standardize]], [[mean]]."
  [vs]
  (let [m (mean vs)]
    (map (fn [^double v]
           (m/- v m)) vs)))

(defn rescale
  "Linearly rescales `vs` from its observed range to a target range, `[0,1]` by default.

  Each value is mapped from `[min(vs), max(vs)]` to `[low, high]` using linear interpolation, see [[fastmath.core/make-norm]].

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `low`, `high` (doubles): Target range bounds. Default to `0.0` and `1.0`.

  Returns a lazy sequence of the same length as `vs`, rescaled to `[low, high]`.

  Values are not clamped beyond `[low, high]`, since `min(vs)` and `max(vs)` always map exactly to `low` and `high`. Constant input (`min(vs) = max(vs)`) leads to division by zero.

  See also [[standardize]], [[demean]], [[extent]]."
  ([vs] (rescale vs 0.0 1.0))
  ([vs ^double low ^double high]
   (let [avs (m/seq->double-array vs)
         mn (minimum avs)
         mx (maximum avs)
         n (m/make-norm mn mx low high)]
     (map n vs))))

(defn covariance
  "Covariance of two sequences.

  This function calculates the *sample* covariance.

  Parameters:

  - `[vs1 vs2]` (sequence of two sequences): A sequence containing the two sequences of numbers.
  - `vs1`, `vs2` (sequences): The two sequences of numbers directly as arguments.

  Both sequences must have the same length.

  Returns the calculated sample covariance as a double.

  See also [[correlation]], [[covariance-matrix]]."
  (^double [[vs1 vs2]] (covariance vs1 vs2))
  (^double [vs1 vs2]
   (let [avs1 (m/seq->double-array vs1)
         avs2 (m/seq->double-array vs2)]
     (m// (v/dot (v/shift avs1 (m/- (mean avs1)))
               (v/shift avs2 (m/- (mean avs2))))
        (m/dec (alength avs1))))))

(defn correlation
  "Calculates the correlation coefficient between two sequences.

  By default, this function calculates the Pearson product-moment correlation
  coefficient, which measures the linear relationship between two datasets.

  This function handles the standard deviation normalization based on whether
  the inputs `vs1` and `vs2` are treated as samples or populations (it uses
  sample standard deviation derived from [[variance]]).

  Parameters:

  - `[vs1 vs2]` (sequence of two sequences): A sequence containing the two sequences of numbers.
  - `vs1`, `vs2` (sequences): The two sequences of numbers directly as arguments.

  Both sequences must have the same length.

  Returns the calculated correlation coefficient (a value between -1.0 and 1.0) as a double.
  Returns `NaN` if one or both sequences have zero variance (are constant).

  See also [[covariance]], [[pearson-correlation]], [[spearman-correlation]], [[kendall-correlation]], [[correlation-matrix]]."
  (^double [[vs1 vs2]] (correlation vs1 vs2))
  (^double [vs1 vs2]
   (let [avs1 (m/seq->double-array vs1)
         avs2 (m/seq->double-array vs2)
         cov (covariance avs1 avs2)
         v1 (variance avs1)
         v2 (variance avs2)]
     (if (or (m/zero? v1) (m/zero? v2))
       ##NaN
       (m// cov (m/sqrt (m/* v1 v2)))))))

(defn spearman-correlation
  "Calculates Spearman's rank correlation coefficient between two sequences.

  Spearman's rank correlation is a non-parametric measure of the monotonic
  relationship between two datasets. It assesses how well the relationship between
  two variables can be described using a monotonic function. It does not require
  the data to be linearly related or follow a specific distribution. The
  coefficient is calculated on the ranks of the data rather than the raw values.

  Parameters:

  - `[vs1 vs2]` (sequence of two sequences): A sequence containing the two sequences of numbers.
  - `vs1`, `vs2` (sequences): The two sequences of numbers directly as arguments.

  Both sequences must have the same length.

  Returns the calculated Spearman rank correlation coefficient (a value between -1.0 and 1.0) as a double.
  A value of 1 indicates a perfect monotonic increasing relationship, -1 a perfect
  monotonic decreasing relationship, and 0 no monotonic relationship.

  See also [[pearson-correlation]], [[kendall-correlation]], [[correlation]]."
  (^double [[vs1 vs2]] (spearman-correlation vs1 vs2))
  (^double [vs1 vs2]
   (.correlation (SpearmansCorrelation.) (m/seq->double-array vs1) (m/seq->double-array vs2))))

(defn pearson-correlation
  "Calculates the Pearson product-moment correlation coefficient between two sequences.

  This function measures the linear relationship between two datasets. The coefficient
  value ranges from -1.0 (perfect negative linear correlation) to 1.0 (perfect
  positive linear correlation), with 0.0 indicating no linear correlation.

  Parameters:

  - `[vs1 vs2]` (sequence of two sequences): A sequence containing the two sequences of numbers.
  - `vs1`, `vs2` (sequences): The two sequences of numbers directly as arguments.

  Both input sequences must contain only numbers and must have the same length.

  Returns the calculated Pearson correlation coefficient as a double. Returns `NaN` if
  either sequence has zero variance (i.e., all elements are the same).

  See also [[correlation]] (general correlation, defaults to Pearson), [[spearman-correlation]],
  [[kendall-correlation]], [[correlation-matrix]]."
  (^double [[vs1 vs2]] (pearson-correlation vs1 vs2))
  (^double [vs1 vs2]
   (.correlation (PearsonsCorrelation.) (m/seq->double-array vs1) (m/seq->double-array vs2))))

(defn kendall-correlation
  "Calculates Kendall's rank correlation coefficient (Kendall's Tau) between two sequences.

  Kendall's Tau is a non-parametric statistic used to measure the ordinal association
  between two measured quantities. It assesses the degree of similarity between the
  orderings of data when ranked by each of the quantities.

  The coefficient value ranges from -1.0 (perfect disagreement in ranking) to 1.0
  (perfect agreement in ranking), with 0.0 indicating no monotonic relationship.
  Unlike Pearson correlation, it does not require the relationship to be linear.

  Parameters:

  - `[vs1 vs2]` (sequence of two sequences): A sequence containing the two sequences of numbers.
  - `vs1`, `vs2` (sequences): The two sequences of numbers directly as arguments.

  Both input sequences must contain only numbers and must have the same length.

  Returns the calculated Kendall's Tau coefficient as a double.

  See also [[pearson-correlation]], [[spearman-correlation]], [[correlation]]."
  (^double [[vs1 vs2]] (kendall-correlation vs1 vs2))
  (^double [vs1 vs2]
   (.correlation (KendallsCorrelation.) (m/seq->double-array vs1) (m/seq->double-array vs2))))

(defn ^{:deprecated "Use [[dissimilarity]]."} kullback-leibler-divergence
  "Kullback-Leibler divergence of two sequences."
  (^double [[vs1 vs2]] (kullback-leibler-divergence vs1 vs2))
  (^double [vs1 vs2]
   (let [res (->> (map vector vs1 vs2)
                  (remove #(m/zero? (v/prod %)))
                  (map (fn [[^double p ^double q]] (m/* p (m/log (m// p q))))))]
     (if (seq res) (sum res) ##Inf))))

(defn ^{:deprecated "Use [[dissimilarity]]."} jensen-shannon-divergence
  "Jensen-Shannon divergence of two sequences."
  (^double [[vs1 vs2]] (jensen-shannon-divergence vs1 vs2))
  (^double [vs1 vs2]
   (let [m (v/mult (mapv + vs1 vs2) 0.5)]
     (m/* 0.5 (m/+ (kullback-leibler-divergence vs1 m)
               (kullback-leibler-divergence vs2 m))))))

(defn coefficient-matrix
  "Generates a matrix of pairwise coefficients from a sequence of sequences.

  This function calculates a matrix where the element at row `i` and column `j`
  is the result of applying the provided `measure-fn` to the `i`-th sequence
  and the `j`-th sequence from the input `vss`.

  Parameters:

  - `vss` (sequence of sequences of numbers): The collection of data sequences. Each
    inner sequence is treated as a variable or set of observations. All inner
    sequences should ideally have the same length if the `measure-fn` expects it.
  - `measure-fn` (function, optional): A function of two arguments (sequences)
    that returns a double representing the coefficient or measure between them.
    Defaults to [[pearson-correlation]].
  - `symmetric?` (boolean, optional): If `true`, the function assumes that
    `measure-fn(a, b)` is equal to `measure-fn(b, a)`. It calculates the upper
    (or lower) triangle of the matrix and mirrors the values to the other side.
    This is an optimization for symmetric measures like correlation and covariance.
    If `false` (default), all pairwise combinations `(i, j)` are calculated independently.

  Returns a sequence of sequences (a matrix) of doubles.

  Note: While this function's `symmetric?` parameter defaults to `false`,
  convenience functions like [[correlation-matrix]] and [[covariance-matrix]]
  wrap this function and explicitly set `symmetric?` to `true` as their
  respective measures are symmetric.

  See also [[correlation-matrix]], [[covariance-matrix]]."
  ([vss] (coefficient-matrix vss pearson-correlation))
  ([vss measure-fn] (coefficient-matrix vss measure-fn false))
  ([vss measure-fn symmetric?]
   (if symmetric?
     (let [avss (map-indexed (fn [id v] [id (m/seq->double-array v)]) vss)
           cache (atom {})]
       (for [[^long id1 ^doubles a] avss]
         (mapv (fn [[^long id2 ^doubles b]]
                 (let [key (if (m/< id1 id2) [id1 id2] [id2 id1])]
                   (if (contains? @cache key)
                     (@cache key)
                     (let [cov (measure-fn a b)]
                       (swap! cache assoc key cov)
                       cov)))) avss)))
     (let [avss (map m/seq->double-array vss)]
       (for [^doubles a avss]
         (mapv #(measure-fn a ^doubles %) avss))))))

(defn correlation-matrix
  "Generates a matrix of pairwise correlation coefficients from a sequence of sequences.

  Given a collection of data sequences `vss`, where each inner sequence represents
  a variable, this function calculates a square matrix where the element at row `i`
  and column `j` is the correlation coefficient between the `i`-th and `j`-th
  sequences in `vss`.

  Parameters:

  - `vss` (sequence of sequences of numbers): The collection of data sequences.
    Each inner sequence is treated as a variable. All inner sequences must have the same length.
  - `measure` (keyword, optional): Specifies the type of correlation coefficient to calculate.
    Defaults to `:pearson`.
    - `:pearson` (default): Calculates the Pearson product-moment correlation coefficient.
    - `:kendall`: Calculates Kendall's Tau rank correlation coefficient.
    - `:spearman`: Calculates Spearman's rank correlation coefficient.

  Returns a sequence of sequences (a matrix) of doubles representing the correlation matrix.
  The matrix is symmetric, as correlation is a symmetric measure.

  See also [[pearson-correlation]], [[spearman-correlation]], [[kendall-correlation]],
  [[covariance-matrix]], [[coefficient-matrix]]."
  ([vss] (correlation-matrix vss :pearson))
  ([vss measure]
   (let [measure (get {:pearson pearson-correlation
                       :kendall kendall-correlation
                       :spearman spearman-correlation} measure pearson-correlation)]
     (coefficient-matrix vss measure true))))

(defn covariance-matrix
  "Generates a matrix of pairwise covariance coefficients from a sequence of sequences.

  Given a collection of data sequences `vss`, where each inner sequence represents
  a variable, this function calculates a square matrix where the element at row `i`
  and column `j` is the sample covariance between the `i`-th and `j`-th sequences
  in `vss`.

  Parameters:

  - `vss` (sequence of sequences of numbers): The collection of data sequences.
    Each inner sequence is treated as a variable. All inner sequences must have the same length.

  Returns a sequence of sequences (a matrix) of doubles representing the covariance matrix.
  The matrix is symmetric, as covariance is a symmetric measure ($Cov(X,Y) = Cov(Y,X)$).

  Internally uses [[coefficient-matrix]] with the [[covariance]] function and `symmetric?` set to `true`.

  See also [[covariance]], [[correlation-matrix]], [[coefficient-matrix]]."
  [vss] (coefficient-matrix vss covariance true))

(defn- maybe-number->seq [vs1 vs2] (if (number? vs2) [(seq vs1) (repeat (count vs1) vs2)] [vs1 vs2]))

(defn me
  "Calculates the Mean Error (ME) between two sequences or a sequence and constant value.

  Parameters:

  - `vs1` (sequence of numbers): The first sequence.
  - `vs2-or-val` (sequence of numbers or single number): The second sequence of
    numbers, or a single number to compare against each element of `vs1`.

  Both sequences (`vs1` and `vs2`) must have the same length if both are sequences.
  If `vs2-or-val` is a single number, it is compared element-wise to `vs1`.

  Returns the calculated Mean Error as a double.

  Note: Positive ME indicates that `vs1` values tend to be greater than `vs2` values
  on average, while negative ME indicates `vs1` values tend to be smaller. ME can be
  influenced by the magnitude of errors and their signs. It does not directly measure
  the magnitude of the typical error due to potential cancellation of positive and
  negative differences.

  See also [[mae]] (Mean Absolute Error), [[mse]] (Mean Squared Error), [[rmse]] (Root Mean Squared Error)."
  (^double [[vs1 vs2-or-val]] (me vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (mean (map - v1 v2))))
  (^double [vs1 vs2-or-val weights]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (m// (sum (map (fn [^double a ^double b ^double w] (m/* w (m/- a b))) v1 v2 weights))
          (sum weights)))))

(defn mae
  "Calculates the Mean Absolute Error (MAE) between two sequences or a sequence and constant value.

  MAE is a measure of the difference between two sequences of values. It quantifies
  the average magnitude of the errors, without considering their direction.

  Parameters:

  - `vs1` (sequence of numbers): The first sequence (often the observed or true values).
  - `vs2-or-val` (sequence of numbers or single number): The second sequence
    (often the predicted or reference values), or a single number to compare
    against each element of `vs1`.

  If both inputs are sequences, they must have the same length. If `vs2-or-val`
  is a single number, it is effectively treated as a sequence of that number
  repeated `count(vs1)` times.

  Returns the calculated Mean Absolute Error as a double.

  Note: MAE is less sensitive to large outliers than metrics like Mean Squared Error (MSE)
  because it uses the absolute value of differences rather than the squared difference.

  See also [[me]] (Mean Error), [[mse]] (Mean Squared Error), [[rmse]] (Root Mean Squared Error)."
  (^double [[vs1 vs2-or-val]] (mae vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (mean (map (comp m/abs -) v1 v2))))
  (^double [vs1 vs2-or-val weights]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (m// (sum (map (fn [^double a ^double b ^double w] (m/* w (m/abs (m/- a b)))) v1 v2 weights))
          (sum weights)))))

(defn mape
  "Calculates the Mean Absolute Percentage Error (MAPE) between two sequences
  or a sequence and a constant value.

  MAPE is a measure of prediction accuracy of a forecasting method, for example
  in time series analysis. It is calculated as the average of the absolute
  percentage errors.

  Parameters:

  - `vs1` (sequence of numbers): The first sequence (conventionally, the actual or true values).
  - `vs2-or-val` (sequence of numbers or single number): The second sequence
    (conventionally, the predicted or reference values), or a single number to
    compare against each element of `vs1`.

  If both inputs are sequences, they must have the same length. If `vs2-or-val`
  is a single number, it is effectively treated as a sequence of that number
  repeated `count(vs1)` times.

  Returns the calculated Mean Absolute Percentage Error as a double.

  Note: MAPE is scale-independent and useful for comparing performance across
  different datasets. However, it is undefined if any of the actual values (`x_i`)
  are zero, and can be skewed by small actual values.

  See also [[me]] (Mean Error), [[mae]] (Mean Absolute Error), [[mse]] (Mean Squared Error), [[rmse]] (Root Mean Squared Error)."
  (^double [[vs1 vs2-or-val]] (mape vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (mean (map (fn [^double a ^double b]
                  (m/abs (m// (m/- a b) a))) v1 v2))))
  (^double [vs1 vs2-or-val weights]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (m// (sum (map (fn [^double a ^double b ^double w]
                      (m/abs (m/* w (m/- a b)))) v1 v2 weights))
          (v/dot v1 weights)))))

(defn rss
  "Calculates the Residual Sum of Squares (RSS) between two sequences or a sequence and a constant value.

  RSS is a measure of the discrepancy between data and a model, often used in
  regression analysis to quantify the total squared difference between observed
  values and predicted (or reference) values.

  Parameters:

  - `vs1` (sequence of numbers): The first sequence (often observed values).
  - `vs2-or-val` (sequence of numbers or single number): The second sequence
    (often predicted or reference values), or a single number to compare
    against each element of `vs1`.

  If both sequences (`vs1` and `vs2`) are provided, they must have the same length.
  If `vs2-or-val` is a single number, it is effectively treated as a sequence
  of that number repeated `count(vs1)` times.

  Returns the calculated Residual Sum of Squares as a double.

  See also [[mse]] (Mean Squared Error), [[rmse]] (Root Mean Squared Error), [[r2]] (Coefficient of Determination)."
  (^double [[vs1 vs2-or-val]] (rss vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (sum (map (comp m/sq -) v1 v2)))))

(defn r2
  "Calculates the Coefficient of Determination ($R^2$) or adjusted version between two sequences or a sequence and a constant value.

  $R^2$ is a statistical measure that represents the proportion of the variance in the
  dependent variable that is predictable from the independent variable(s) in a
  statistical model. It indicates how well the model fits the observed data.

  The standard $R^2$ is calculated as $1 - (RSS / TSS)$, where:
  - $RSS$ (Residual Sum of Squares) is the sum of the squared differences
    between the observed values (`vs1`) and the predicted/reference values (`vs2` or `vs2-or-val`).
    See [[rss]].
  - $TSS$ (Total Sum of Squares) is the sum of the squared differences
    between the observed values (`vs1`) and their mean. This is calculated
    using [[moment]] of order 2 with `:mean?` set to `false`.

  This function has two arities:

  1.  `(r2 vs1 vs2-or-val)`: Calculates the standard $R^2$.

      - `vs1` (seq of numbers): The sequence of observed or actual values.
      - `vs2-or-val` (seq of numbers or single number): The sequence of predicted or
        reference values, or a single constant value to compare against.

      Returns the calculated standard $R^2$ as a double. For simple linear regression,
      this is equal to the square of the Pearson correlation coefficient ([[r2-determination]]).
      $R^2$ typically ranges from 0 to 1 in this context, but can be negative
      if the chosen model fits the data worse than a horizontal line through the mean
      of the observed data.

  2.  `(r2 vs1 vs2-or-val no-of-variables)`: Calculates the **Adjusted $R^2$**.
      The adjusted $R^2$ is a modified version of $R^2$ that has been adjusted
      for the number of predictors in the model. It increases only if the new term
      improves the model more than would be expected by chance.
      The formula for adjusted $R^2$ is:
      $$ R^2_{adj} = 1 - (1 - R^2) \\frac{n-1}{n-p-1} $$
      where $n$ is the number of observations (length of `vs1`) and $p$ is the
      number of independent variables (`no-of-variables`).

      - `vs1` (seq of numbers): The sequence of observed or actual values.
      - `vs2-or-val` (seq of numbers or single number): The sequence of predicted or
        reference values, or a single constant value to compare against.
      - `no-of-variables` (double): The number of independent variables ($p$) used in the model
        that produced the `vs2-or-val` predictions.

      Returns the calculated adjusted $R^2$ as a double.

  Both `vs1` and `vs2` (if `vs2-or-val` is a sequence) must have the same length.

  See also [[rss]], [[mse]], [[rmse]], [[pearson-correlation]], [[r2-determination]]."
  (^double [[vs1 vs2-or-val]] (r2 vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (m/- 1.0 (m// (rss vs1 vs2-or-val)
             (moment vs1 2 {:mean? false}))))
  (^double [vs1 vs2-or-val ^double no-of-variables]
   (let [rr (r2 vs1 vs2-or-val)
         n (count vs1)]
     (m/- 1.0 (m/* (m/- 1.0 rr)
               (m// (m/dec n)
                  (m/- n no-of-variables 1.0)))))))

(defn mse
  "Calculates the Mean Squared Error (MSE) between two sequences or a sequence and a constant value.

  MSE is a measure of the quality of an estimator or predictor. It quantifies
  the average of the squared differences between corresponding elements of the
  input sequences.

  Parameters:

  - `vs1` (sequence of numbers): The first sequence (often the observed or true values).
  - `vs2-or-val` (sequence of numbers or single number): The second sequence
    (often the predicted or reference values), or a single number to compare
    against each element of `vs1`.

  If both inputs are sequences, they must have the same length. If `vs2-or-val`
  is a single number, it is effectively treated as a sequence of that number
  repeated `count(vs1)` times.

  Returns the calculated Mean Squared Error as a double.

  Note: MSE penalizes larger errors more heavily than smaller errors because the
  errors are squared. This makes it sensitive to outliers. It is the average
  of the [[rss]] (Residual Sum of Squares). Its square root is the [[rmse]].

  See also [[rss]] (Residual Sum of Squares), [[rmse]] (Root Mean Squared Error),
  [[me]] (Mean Error), [[mae]] (Mean Absolute Error), [[r2]] (Coefficient of Determination)."
  (^double [[vs1 vs2-or-val]] (mse vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (mean (map (comp m/sq -) v1 v2))))
  (^double [vs1 vs2-or-val weights]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (m// (sum (map (fn [^double a ^double b ^double w] (m/* w (m/sq (m/- a b)))) v1 v2 weights))
          (sum (weights))))))

(defn rmse
  "Calculates the Root Mean Squared Error (RMSE) between two sequences or a sequence and a constant value.

  RMSE is the square root of the [[mse]] (Mean Squared Error). It represents the
  standard deviation of the residuals (prediction errors) and has the same units
  as the original data, making it more interpretable than MSE. It measures the
  average magnitude of the errors, penalizing larger errors more than smaller ones
  due to the squaring involved.

  Parameters:

  - `vs1` (sequence of numbers): The first sequence (often the observed or true values).
  - `vs2-or-val` (sequence of numbers or single number): The second sequence
    (often the predicted or reference values), or a single number to compare
    against each element of `vs1`.

  If both inputs are sequences, they must have the same length. If `vs2-or-val`
  is a single number, it is effectively treated as a sequence of that number
  repeated `count(vs1)` times.

  Returns the calculated Root Mean Squared Error as a double.

  See also [[mse]] (Mean Squared Error), [[rss]] (Residual Sum of Squares),
  [[me]] (Mean Error), [[mae]] (Mean Absolute Error), [[r2]] (Coefficient of Determination)."
  (^double [[vs1 vs2-or-val]] (rmse vs1 vs2-or-val))
  (^double [vs1 vs2-or-val] (m/sqrt (mse vs1 vs2-or-val)))
  (^double [vs1 vs2-or-val weights] (m/sqrt (mse vs1 vs2-or-val weights))))

(defn count=
  "Count equal values in both seqs. Same as [[L0]]

  Calculates the number of pairs of corresponding elements that are equal between
  two sequences, or between a sequence and a single scalar value.

  Parameters:

  - `vs1` (sequence of numbers): The first sequence.
  - `vs2-or-val` (sequence of numbers or single number): The second sequence of
    numbers, or a single number to compare against each element of `vs1`.

  If both inputs are sequences, they must have the same length. If `vs2-or-val`
  is a single number, it is effectively treated as a sequence of that number
  repeated `count(vs1)` times.

  Returns the count of equal elements as a long integer."
  (^long [[vs1 vs2-or-val]] (count= vs1 vs2-or-val))
  (^long [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (count (filter (fn [^double v] (m/zero? v)) (map m/- v1 v2))))))

(def ^{:doc "Count equal values in both seqs. Alias for [[count==]]"} L0 count=)

(defn L1
  "Calculates the L1 distance (Manhattan or City Block distance) between two sequences or a sequence and a constant value.

  The L1 distance is the sum of the absolute differences between corresponding elements.

  Parameters:

  - `vs1` (sequence of numbers): The first sequence.
  - `vs2-or-val` (sequence of numbers or single number): The second sequence of numbers, or a single number to compare against each element of `vs1`.

  If both inputs are sequences, they must have the same length. If `vs2-or-val`
  is a single number, it is effectively treated as a sequence of that number
  repeated `count(vs1)` times.

  Returns the calculated L1 distance as a double.

  See also [[L2]], [[L2sq]], [[LInf]], [[mae]] (Mean Absolute Error)."
  (^double [[vs1 vs2-or-val]] (L1 vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (d/manhattan v1 v2))))

(defn L2sq
  "Calculates the Squared Euclidean distance between two sequences or a sequence and a constant value.

  This is the sum of the squared differences between corresponding elements.
  It is equivalent to the [[rss]] (Residual Sum of Squares).

  Parameters:

  - `vs1` (sequence of numbers): The first sequence.
  - `vs2-or-val` (sequence of numbers or single number): The second sequence of numbers, or a single number to compare against each element of `vs1`.

  If both inputs are sequences, they must have the same length. If `vs2-or-val`
  is a single number, it is effectively treated as a sequence of that number
  repeated `count(vs1)` times.

  Returns the calculated Squared Euclidean distance as a double.

  See also [[L1]], [[L2]], [[LInf]], [[rss]] (Residual Sum of Squares), [[mse]] (Mean Squared Error)."
  (^double [[vs1 vs2-or-val]] (L2sq vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (d/euclidean-sq v1 v2))))

(defn L2
  "Calculates the L2 distance (Euclidean distance) between two sequences or a sequence and a constant value.

  This is the standard straight-line distance between two points (vectors) in Euclidean space.
  It is the square root of the [[L2sq]] distance.

  Parameters:

  - `vs1` (sequence of numbers): The first sequence.
  - `vs2-or-val` (sequence of numbers or single number): The second sequence of numbers, or a single number to compare against each element of `vs1`.

  If both inputs are sequences, they must have the same length. If `vs2-or-val`
  is a single number, it is effectively treated as a sequence of that number
  repeated `count(vs1)` times.

  Returns the calculated L2 distance as a double.

  See also [[L1]], [[L2sq]], [[LInf]], [[rmse]] (Root Mean Squared Error)."
  (^double [[vs1 vs2-or-val]] (L2 vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (d/euclidean v1 v2))))

(defn LInf
  "Calculates the L-infinity distance (Chebyshev distance) between two sequences or a sequence and a constant value.

  The Chebyshev distance is the maximum absolute difference between corresponding elements.

  Parameters:

  - `vs1` (sequence of numbers): The first sequence.
  - `vs2-or-val` (sequence of numbers or single number): The second sequence of numbers, or a single number to compare against each element of `vs1`.

  If both inputs are sequences, they must have the same length. If `vs2-or-val`
  is a single number, it is effectively treated as a sequence of that number
  repeated `count(vs1)` times.

  Returns the calculated L-infinity distance as a double.

  See also [[L1]], [[L2]], [[L2sq]]."
  (^double [[vs1 vs2-or-val]] (LInf vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [[v1 v2] (maybe-number->seq vs1 vs2-or-val)]
     (d/chebyshev v1 v2))))

(defn psnr
  "Peak Signal-to-Noise Ratio (PSNR).

  PSNR is a measure used to quantify the quality of reconstruction of lossy
  compression codecs (e.g., for images or video). It is calculated using the
  Mean Squared Error (MSE) between the original and compressed images/signals.
  A higher PSNR generally indicates a higher quality signal reconstruction
  (i.e., less distortion).

  Parameters:

  - `vs1` (sequence of numbers): The first sequence (conventionally, the original or reference signal/data).
  - `vs2-or-val` (sequence of numbers or single number): The second sequence
    (conventionally, the reconstructed or noisy signal/data), or a single number
    to compare against each element of `vs1`.
  - `max-value` (optional, double): The maximum possible value of a sample in the data.
    If not provided, the function automatically determines the maximum value present
    across both input sequences (`vs1` and `vs2` if a sequence, or `vs1` and the scalar value
    if `vs2-or-val` is a number). Providing an explicit `max-value` is often more
    appropriate based on the data type's theoretical maximum range (e.g., 255 for 8-bit).

  If `vs2-or-val` is a sequence, both `vs1` and `vs2` must have the same length.

  Returns the calculated Peak Signal-to-Noise Ratio as a double. Returns `-Double/Infinity`
  if the MSE is zero (perfect match). Returns `NaN` if MSE is non-positive.

  See also [[mse]], [[rmse]]."
  (^double [[vs1 vs2-or-val]] (psnr vs1 vs2-or-val))
  (^double [vs1 vs2-or-val]
   (let [mx1 (maximum vs1)
         mx2 (double (if (number? vs2-or-val) vs2-or-val (maximum vs2-or-val)))]
     (psnr vs1 vs2-or-val (m/max mx1 mx2))))
  (^double [vs1 vs2-or-val ^double max-value]
   (m/- (m/* 20.0 (m/log10 max-value))
        (m/* 10.0 (m/log10 (mse vs1 vs2-or-val))))))

;;

(defn estimate-bins
  "Estimates a suitable number of bins for a histogram of `vs`.

  Several standard rules of thumb are available, balancing resolution against noise: some rely only on sample size, while others also take the spread or shape of the data into account.

  Parameters:

  - `vs` (sequence of numbers): Input data.
  - `bins-or-estimate-method` (keyword or long): Estimation method to use, or an explicit bin count. Defaults to `:freedman-diaconis`. One of:
    - `:sqrt`: square root of the sample size, `sqrt(n)`.
    - `:sturges`: `1 + log2(n)`, assumes roughly normal data, tends to undersmooth for large or non-normal samples.
    - `:rice`: `2 * cbrt(n)`, depends only on sample size.
    - `:doane`: a refinement of Sturges' rule that accounts for the sample skewness (see [[skewness]]), better suited for non-normal data.
    - `:scott`: bin width `3.5 * stddev / cbrt(n)`, assumes roughly normal data.
    - `:freedman-diaconis` (default): bin width `2 * IQR / cbrt(n)` (see [[iqr]]), robust to outliers.
    - a `long`: used directly as the number of bins, bypassing estimation.

  Returns the estimated number of bins as a long, never higher than the number of samples in `vs`.

  See also [[histogram]]."
  (^long [vs] (estimate-bins vs :freedman-diaconis))
  (^long [vs bins-or-estimate-method]
   (if-not (keyword? bins-or-estimate-method)
     (or bins-or-estimate-method (estimate-bins vs))
     (let [n (count vs)]
       (m/min n (int (case bins-or-estimate-method
                       :sqrt (m/max 1 (m/sqrt n))
                       :sturges (bins/sturges n)
                       :rice (bins/rice n)
                       :doane (bins/doane (m/seq->double-array vs) n)
                       :scott (bins/scott (m/seq->double-array vs) n)
                       :freedman-diaconis (bins/freedman-diaconis (m/seq->double-array vs) n))))))))

(defn- constrain-data
  [vs ^double mn ^double mx]
  (filter (fn [^double v] (m/<= mn v mx)) vs))

(defn- process-vs-and-bins
  [vs bins-or-estimate-method ^double mn ^double mx]
  (if (sequential? bins-or-estimate-method)
    (let [[^double nmn ^double nmx] (extent bins-or-estimate-method)]
      [(constrain-data vs nmn nmx)
       nmn nmx (sort bins-or-estimate-method)])
    (let [nvs (constrain-data vs mn mx)
          bins (if (m/== mn mx) 1 (estimate-bins nvs bins-or-estimate-method))]
      [nvs mn mx (m/slice-range mn mx (m/long-inc bins))])))

(defn- histogram-internal
  [[vs ^double mn ^double mx intervals]]
  (let [diff (m/- mx mn)
        bins (if (m/zero? diff) 1 (m/dec (count intervals)))
        bins- (m/dec bins)
        step (m// diff bins)
        search-array (double-array intervals)
        buff (long-array bins)
        sum (double-array bins)
        size (count vs)
        dsize (double size)]

    (doseq [^double v vs]
      (let [b (java.util.Arrays/binarySearch ^doubles search-array v)
            pos (unchecked-int (m/min bins- (long (if (m/neg? b) (m/abs (m/+ b 2)) b))))]
        (fastmath.java.Array/inc ^longs buff pos)
        (fastmath.java.Array/add ^doubles sum pos v)))

    (let [bins-map (map (fn [[^double mn ^double mx] ^long cnt ^double s]
                          (let [step (m/- mx mn)]
                            {:min mn :max mx :count cnt
                             :step step
                             :mid (m/+ mn (m/* 0.5 step))
                             :avg (m// s cnt)
                             :probability (m// cnt dsize)})) (partition 2 1 search-array) buff sum)]
      {:size bins
       :step step
       :samples size
       :min mn
       :max mx
       :bins (map (juxt :min :count) bins-map)
       :bins-maps bins-map
       :intervals intervals
       :frequencies (into {} (map (juxt :avg :count) bins-map))})))

(defn histogram
  "Builds a histogram of `vs`, partitioning the data into equal-width bins and counting the samples in each.

  `vs` can also be a sequence of sequences, in which case a shared set of bin boundaries is derived across all of them (from their combined extent and the largest estimated bin count), and a histogram is computed for each inner sequence using those boundaries, letting the results be compared directly.

  Parameters:

  - `vs` (sequence of numbers, or sequence of sequences of numbers): Input data.
  - `bins-or-estimate-method` (keyword, long or sequence): Number of bins, or bin-count estimation method, see [[estimate-bins]] (default `:freedman-diaconis`). Alternatively, a sequence of numbers to use directly as explicit bin boundaries (intervals).
  - `mn`, `mx` (doubles, or a 2-element vector `[mn mx]`): Optional bounds. When provided, data is filtered to keep only values in `[mn,mx]` before binning, and this range is used instead of the data extent. Ignored when `bins-or-estimate-method` is a sequence of boundaries, since its own extent is used instead.

  Returns a map (or, when `vs` is a sequence of sequences, a sequence of such maps, one per inner sequence, sharing the same `:intervals`) with the following keys:

  - `:size`: number of bins.
  - `:step`: average distance between bins.
  - `:samples`: number of samples actually used (after filtering to `[mn,mx]`).
  - `:min`, `:max`: lower and upper bounds used for binning.
  - `:bins`: sequence of pairs `[bin-lower-bound bin-count]`.
  - `:bins-maps`: sequence of maps, one per bin, with keys `:min` (lower bound), `:max` (upper bound), `:mid` (middle value), `:step` (actual bin width), `:count` (number of elements), `:avg` (average value of the elements in the bin), and `:probability` (fraction of samples in the bin).
  - `:intervals`: the bin boundaries used to build the histogram.
  - `:frequencies`: a map from each bin's `:avg` to its `:count`.

  If the difference between `mn` and `mx` is `0`, the number of bins is set to `1`.

  See also [[estimate-bins]]."
  ([vs] (histogram vs :freedman-diaconis))
  ([vs bins-or-estimate-method] (histogram vs bins-or-estimate-method (extent
                                                                       (if (sequential? (first vs))
                                                                         (flatten vs) vs))))
  ([vs bins-or-estimate-method [^double mn ^double mx]] (histogram vs bins-or-estimate-method mn mx))
  ([vs bins-or-estimate-method ^double mn ^double mx]
   (if (sequential? (first vs))
     (let [sbins? (sequential? bins-or-estimate-method)
           [^double nmn ^double nmx] (if sbins? (extent bins-or-estimate-method) [mn mx])
           nvs (map #(constrain-data % nmn nmx) vs)
           intervals (if sbins?
                       (sort bins-or-estimate-method)
                       (if (m/== nmn nmx)
                         [nmn nmx]
                         (->> nvs
                              (map #(estimate-bins % bins-or-estimate-method))
                              (reduce m/max 1.0)
                              (long)
                              (m/long-inc)
                              (m/slice-range nmn nmx))))]
       (map (fn [vs] (histogram-internal [vs nmn nmx intervals])) nvs))
     (histogram-internal (process-vs-and-bins vs bins-or-estimate-method mn mx)))))

;; distances

(defn- quantize-distribution
  "Find probabilities from distribution for given intervals from histogram."
  [xs distr bins]
  (let [{:keys [^double step bins]} (if (map? xs) xs (histogram xs bins))
        last-idx (m/dec (count bins))
        counts (map second bins)]
    [counts (map-indexed (fn [^long id [^double s]]
                           (condp = id
                             0 (r/cdf distr (m/+ s step))
                             last-idx (m/- 1.0 (r/cdf distr s))
                             (r/cdf distr s (m/+ s step)))) bins)]))
(defn- pq-from-histograms
  [P Q bins]
  (let [[Ph Qh] (histogram [P Q] bins)]
    [(map second (:bins Ph))
     (map second (:bins Qh))]))

(defn- normalize-PQ
  [P-observed Q-expected bins probabilities?]
  (let [[P Q] (cond
                (r/distribution? Q-expected) (quantize-distribution P-observed Q-expected bins)
                bins (pq-from-histograms P-observed Q-expected bins) 
                :else [P-observed Q-expected])]
    (cond
      (r/distribution? Q-expected) [(v/div P (v/sum P)) Q]
      probabilities? [(v/div P (v/sum P)) (v/div Q (v/sum Q))]
      :else [P Q])))

(defn- remove-zeros-pairwise
  [[P Q]]
  (let [pairs (->> (map v/vec2 P Q)
                   (remove (fn [^Vec2 v] (or (m/zero? (.x v))
                                            (m/zero? (.y v))))))]
    [(map first pairs) (map second pairs)]))

(defn- safe-div
  ^double [^double n ^double d ^double e]
  (if (m/zero? d) (m// n e) (m// n d)))

(defn- make-safe-log
  [^double ex ^double e]
  (cond
    (m/== ex m/E) (fn ^double [^double v] (m/log (if (m/zero? v) e v)))
    (m/== ex 2.0) (fn ^double [^double v] (m/log2 (if (m/zero? v) e v)))
    (m/== ex 10.0) (fn ^double [^double v] (m/log10 (if (m/zero? v) e v)))
    :else (fn ^double [^double v] (m/logb ex (if (m/zero? v) e v)))))

(defn dissimilarity
  "Calculates a dissimilarity (distance) measure between two probability density functions, `P-observed` and `Q-expected`, given as histograms, frequencies, probabilities, or raw data.

  If `Q-expected` is a distribution object, a histogram is built from `P-observed` and compared against it directly. If `:bins` is set (and `Q-expected` is not a distribution), both `P-observed` and `Q-expected` are treated as raw data and turned into a matching pair of histograms.

  Parameters:

  - `method` (keyword): Dissimilarity method to use, see below.
  - `P-observed` (sequence of numbers): Frequencies, probabilities, or raw data.
  - `Q-expected` (sequence of numbers, or distribution object): Frequencies, probabilities, raw data, or a distribution object to compare `P-observed` against.
  - `opts` (map, optional):
    - `:probabilities?` (boolean): Whether `P-observed`/`Q-expected` are normalized to probabilities (summing to `1`) before computing the distance. Defaults to `true`.
    - `:bins` (long or keyword): Number of bins, or a bin-count estimation method, used to turn raw `P-observed`/`Q-expected` data into histograms, see [[histogram]] and [[estimate-bins]].
    - `:remove-zeros?` (boolean): Whether to drop paired entries where either `P-observed` or `Q-expected` is zero before computing the distance. Defaults to `false`.
    - `:epsilon` (double): Small number substituted for `0.0` wherever a division or logarithm would otherwise be undefined. Defaults to `1.0e-6`.
    - `:log-base` (double): Base of the logarithm used by log-based methods. Defaults to `e`.
    - `:power` (double): Exponent used by the `:minkowski` method. Defaults to `2.0`.

  `method` is one of: `:euclidean`, `:city-block`, `:manhattan`, `:chebyshev`, `:minkowski`, `:sorensen`, `:gower`, `:soergel`, `:kulczynski`, `:canberra`, `:lorentzian`, `:non-intersection`, `:wave-hedges`, `:czekanowski`, `:motyka`, `:tanimoto`, `:jaccard`, `:dice`, `:bhattacharyya`, `:hellinger`, `:matusita`, `:squared-chord`, `:euclidean-sq`, `:squared-euclidean`, `:pearson-chisq`, `:chisq`, `:neyman-chisq`, `:squared-chisq`, `:symmetric-chisq`, `:divergence`, `:clark`, `:additive-symmetric-chisq`, `:kullback-leibler`, `:jeffreys`, `:k-divergence`, `:topsoe`, `:jensen-shannon`, `:jensen-difference`, `:taneja`, `:kumar-johnson`, `:avg` (the mean of the `:city-block` and `:chebyshev` distances).

  Returns the calculated dissimilarity as a double.

  Definitions for the individual methods are given in the Comprehensive Survey on Distance/Similarity Measures between Probability Density Functions by Sung-Hyuk Cha.

  See also [[similarity]], [[histogram]], [[estimate-bins]]."
  (^double [method P-observed Q-expected] (dissimilarity method P-observed Q-expected nil))
  (^double [method P-observed Q-expected {:keys [bins probabilities? ^double epsilon ^double log-base ^double power remove-zeros?]
                                          :or {probabilities? true epsilon 1.0e-6 log-base m/E power 2.0}}]
   (let [pq (normalize-PQ P-observed Q-expected bins probabilities?)
         [P Q] (if remove-zeros? (remove-zeros-pairwise pq) pq)
         log-fn (make-safe-log log-base epsilon)]
     (case method
       :euclidean (L2 P Q)
       :city-block (L1 P Q)
       :manhattan (L1 P Q)
       :chebyshev (LInf P Q)
       :minkowski (m/pow (v/sum (map (fn [^double p ^double q] (m/pow (m/abs (m/- p q)) power)) P Q))
                         (m// power))
       :sorensen (safe-div (L1 P Q) (m/+ (v/sum P) (v/sum Q)) epsilon)
       :gower (safe-div (L1 P Q) (count P) epsilon)
       :soergel (safe-div (L1 P Q) (v/sum (v/emx P Q)) epsilon)
       :kulczynski (safe-div (L1 P Q) (v/sum (v/emn P Q)) epsilon)
       :canberra (v/sum (map (fn [^double p ^double q]
                               (safe-div (m/abs (m/- p q)) (m/+ p q) epsilon)) P Q))
       :lorentzian (v/sum (map (fn [^double p ^double q] (log-fn (m/inc (m/abs (m/- p q))))) P Q))
       :non-intersection (m/* 0.5 (L1 P Q))
       :wave-hedges (v/sum (map (fn [^double p ^double q]
                                  (safe-div (m/abs (m/- p q)) (m/max p q) epsilon)) P Q))
       :czekanowski (safe-div (L1 P Q) (m/+ (v/sum P) (v/sum Q)) epsilon)
       :motyka (safe-div (v/sum (v/emx P Q)) (m/+ (v/sum P) (v/sum Q)) epsilon)
       :tanimoto (let [mx (v/emx P Q)]
                   (safe-div (v/sum (v/sub mx (v/emn P Q))) (v/sum mx) epsilon))
       :jaccard (safe-div (L2sq P Q) (m/- (m/+ (v/sum (v/sq P))
                                               (v/sum (v/sq Q)))
                                          (v/sum (v/emult P Q))) epsilon)
       :dice (safe-div (L2sq P Q) (m/+ (v/sum (v/sq P))
                                       (v/sum (v/sq Q))) epsilon)
       :bhattacharyya (m/- (double (log-fn (v/sum (v/sqrt (v/emult P Q))))))
       :hellinger (m/* 2.0 (m/sqrt (m/- 1.0 (v/sum (v/sqrt (v/emult P Q))))))
       :matusita (m/sqrt (m/- 2.0 (m/* 2.0 (v/sum (v/sqrt (v/emult P Q))))))
       :squared-chord (L2sq (v/sqrt P) (v/sqrt Q))
       :euclidean-sq (L2sq P Q)
       :squared-euclidean (L2sq P Q)
       :pearson-chisq (v/sum (map (fn [^double p ^double q]
                                    (safe-div (m/sq (m/- p q)) q epsilon)) P Q))
       :chisq (v/sum (map (fn [^double p ^double q]
                            (safe-div (m/sq (m/- p q)) q epsilon)) P Q))
       :neyman-chisq (v/sum (map (fn [^double p ^double q]
                                   (safe-div (m/sq (m/- p q)) p epsilon)) P Q))
       :squared-chisq (v/sum (map (fn [^double p ^double q]
                                    (safe-div (m/sq (m/- p q)) (m/+ p q) epsilon)) P Q))
       :symmetric-chisq (m/* 2.0 (v/sum (map (fn [^double p ^double q]
                                               (safe-div (m/sq (m/- p q)) (m/+ p q) epsilon)) P Q)))
       :divergence (m/* 2.0 (v/sum (map (fn [^double p ^double q]
                                          (safe-div (m/sq (m/- p q)) (m/sq (m/+ p q)) epsilon)) P Q)))
       :clark (m/sqrt (v/sum (map (fn [^double p ^double q]
                                    (m/sq (safe-div (m/abs (m/- p q)) (m/+ p q) epsilon))) P Q)))
       :additive-symmetric-chisq (v/sum (map (fn [^double p ^double q]
                                               (safe-div (m/* (m/sq (m/- p q)) (m/+ p q))
                                                         (m/* p q)
                                                         epsilon)) P Q))
       :kullback-leibler (v/sum (map (fn [^double p ^double q]
                                       (m/* p (double (log-fn (safe-div p q epsilon))))) P Q))
       :jeffreys (v/sum (map (fn [^double p ^double q]
                               (m/* (m/- p q) (double (log-fn (safe-div p q epsilon))))) P Q))
       :k-divergence (v/sum (map (fn [^double p ^double q]
                                   (m/* p (double (log-fn (safe-div (m/* 2.0 p) (m/+ p q) epsilon))))) P Q))
       :topsoe (v/sum (map (fn [^double p ^double q]
                             (m/+ (m/* p (double (log-fn (safe-div (m/* 2.0 p) (m/+ p q) epsilon))))
                                  (m/* q (double (log-fn (safe-div (m/* 2.0 q) (m/+ p q) epsilon)))))) P Q))
       :jensen-shannon (m/* 0.5 (m/+ (v/sum (map (fn [^double p ^double q]
                                                   (m/* p (double (log-fn (safe-div (m/* 2.0 p) (m/+ p q) epsilon))))) P Q))
                                     (v/sum (map (fn [^double p ^double q]
                                                   (m/* q (double (log-fn (safe-div (m/* 2.0 q) (m/+ p q) epsilon))))) P Q))))
       :jensen-difference (v/sum (map (fn [^double p ^double q]
                                        (let [pq2 (m/* 0.5 (m/+ p q))]
                                          (m/- (m/* 0.5 (m/+ (m/* p (double (log-fn p)))
                                                             (m/* q (double (log-fn q))))) (m/* pq2 (double (log-fn pq2)))))) P Q))
       :taneja (v/sum (map (fn [^double p ^double q]
                             (let [pq2 (m/* 0.5 (m/+ p q))]
                               (m/* pq2 (double (log-fn (safe-div pq2 (m/sqrt (m/* p q)) epsilon)))))) P Q))
       :kumar-johnson (v/sum (map (fn [^double p ^double q]
                                    (safe-div (m/sq (m/- (m/sq p) (m/sq q)))
                                              (m/* 2.0 (m/sqrt (m/cb (m/* p q)))) epsilon)) P Q))
       :avg (m/* 0.5 (m/+ (L1 P Q) (LInf P Q)))))))

(defn similarity
  "Calculates a similarity measure between two probability density functions, `P-observed` and `Q-expected`, given as histograms, frequencies, probabilities, or raw data.

  If `Q-expected` is a distribution object, a histogram is built from `P-observed` and compared against it directly. If `:bins` is set (and `Q-expected` is not a distribution), both `P-observed` and `Q-expected` are treated as raw data and turned into a matching pair of histograms.

  Parameters:

  - `method` (keyword): Similarity method to use, see below.
  - `P-observed` (sequence of numbers): Frequencies, probabilities, or raw data.
  - `Q-expected` (sequence of numbers, or distribution object): Frequencies, probabilities, raw data, or a distribution object to compare `P-observed` against.
  - `opts` (map, optional):
    - `:probabilities?` (boolean): Whether `P-observed`/`Q-expected` are normalized to probabilities (summing to `1`) before computing the similarity. Defaults to `true`.
    - `:bins` (long or keyword): Number of bins, or a bin-count estimation method, used to turn raw `P-observed`/`Q-expected` data into histograms, see [[histogram]] and [[estimate-bins]].
    - `:epsilon` (double): Small number substituted for `0.0` wherever a division would otherwise be undefined. Defaults to `1.0e-6`.

  `method` is one of: `:intersection`, `:czekanowski`, `:motyka`, `:kulczynski`, `:ruzicka`, `:inner-product`, `:harmonic-mean`, `:cosine`, `:jaccard`, `:dice`, `:fidelity`, `:squared-chord`.

  Returns the calculated similarity as a double. Higher values generally indicate more similar distributions, though the range and interpretation depend on the chosen `method`.

  Definitions for the individual methods are given in the Comprehensive Survey on Distance/Similarity Measures between Probability Density Functions by Sung-Hyuk Cha.

  See also [[dissimilarity]], [[histogram]], [[estimate-bins]]."
  (^double [method P-observed Q-expected] (similarity method P-observed Q-expected nil))
  (^double [method P-observed Q-expected {:keys [bins probabilities? ^double epsilon]
                                          :or {probabilities? true epsilon 1.0e-6}}]
   (let [[P Q] (normalize-PQ P-observed Q-expected bins probabilities?)]
     (case method
       :intersection (v/sum (v/emn P Q))
       :czekanowski (safe-div (m/* 2.0 (v/sum (v/emn P Q))) (m/+ (v/sum P) (v/sum Q)) epsilon)
       :motyka (safe-div (v/sum (v/emn P Q)) (m/+ (v/sum P) (v/sum Q)) epsilon)
       :kulczynski (safe-div (v/sum (v/emn P Q)) (L1 P Q) epsilon)
       :ruzicka (safe-div (v/sum (v/emn P Q)) (v/sum (v/emx P Q)) epsilon)
       :inner-product (v/dot P Q)
       :harmonic-mean (m/* 2.0 (v/sum (map (fn [^double p ^double q]
                                           (safe-div (m/* p q) (m/+ p q) epsilon)) P Q)))
       :cosine (safe-div (v/sum (v/emult P Q)) (m/* (v/mag P) (v/mag Q)) epsilon)
       :jaccard (let [pq (v/sum (v/emult P Q))]
                  (safe-div pq (m/- (m/+ (v/magsq P) (v/magsq Q)) pq) epsilon))
       :dice (let [pq (v/sum (v/emult P Q))]
               (safe-div (m/* 2.0 pq) (m/+ (v/magsq P) (v/magsq Q)) epsilon))
       :fidelity (v/sum (v/sqrt (v/emult P Q)))
       :squared-chord (m/dec (m/* 2.0 (v/sum (v/sqrt (v/emult P Q)))))))))

;;

(defn- weighted-variance-average
  ^double [groups]
  (reduce (fn [^double sum xs]
            (m/+ sum (m/* (m/dec (count xs))
                      (variance xs)))) 0.0 groups))

(defn pooled-variance
  "Calculates the pooled variance of several samples, a weighted combination of their individual variances, used when the samples are assumed to share a common (unknown) variance.

  Parameters:

  - `groups` (sequence of sequences of numbers): The independent samples to pool.
  - `method` (keyword): Pooling method. Defaults to `:unbiased`. One of:
    - `:unbiased` (default): each group's variance is weighted by its degrees of freedom `n_i - 1`, then the weighted sum is divided by the total degrees of freedom `sum(n_i) - k`, where `k` is the number of groups. This is the classic unbiased pooled-variance estimator.
    - `:biased`: same weighted sum of variances as `:unbiased`, but divided by the total sample size `sum(n_i)` instead of the total degrees of freedom.
    - `:avg`: plain, unweighted average of the individual group variances, each group contributing equally regardless of its size.

  Returns the pooled variance as a double.

  See also [[pooled-stddev]], [[variance]], [[cohens-d]]."
  (^double [groups] (pooled-variance groups :unbiased))
  (^double [groups method]
   (let [agroups (map m/seq->double-array groups)]
     (case method
       :biased (m// (weighted-variance-average agroups)
                  (sum (map alength agroups)))
       :avg (m// (sum (map variance agroups))
               (count groups))
       (m// (weighted-variance-average agroups)
          (m/- (sum (map alength agroups)) (count groups)))))))

(defn pooled-stddev
  "Calculates the pooled standard deviation of several samples, the square root of [[pooled-variance]].

  Parameters:

  - `groups` (sequence of sequences of numbers): The independent samples to pool.
  - `method` (keyword): Pooling method, see [[pooled-variance]] for details. Defaults to `:unbiased`.

  Returns the pooled standard deviation as a double.

  See also [[pooled-variance]], [[stddev]], [[cohens-d]]."
  (^double [groups] (m/sqrt (pooled-variance groups)))
  (^double [groups method] (m/sqrt (pooled-variance groups method))))

(defn pooled-mad
  "Calculates a pooled, robust measure of scale across several samples, based on the median absolute deviation (MAD) of their pooled residuals.

  Each group is first centered by subtracting its own median, and all the centered groups are then concatenated into a single pooled sample of residuals, whose MAD is computed and scaled by `const`.

  Parameters:

  - `groups` (sequence of sequences of numbers): The independent samples to pool.
  - `const` (double): Scaling constant applied to the MAD of the pooled residuals. Defaults to `1.4826022185056023`, the constant that makes MAD a consistent estimator of the standard deviation for normally distributed data, see [[median-absolute-deviation]].

  Returns the pooled MAD as a double.

  See also [[median-absolute-deviation]], [[pooled-variance]], [[pooled-stddev]]."
  (^double [groups] (pooled-mad groups 1.4826022185056023))
  (^double [groups ^double const]
   (let [Y (mapcat (fn [g] (let [md (median g)]
                             (v/shift g (m/- md)))) groups)]
     (m/* const (median-absolute-deviation Y)))))

;; effect size

(defn cohens-d
  "Calculate Cohen's d effect size between two groups.

  Cohen's d is a standardized measure used to quantify the magnitude of the
  difference between the means of two independent groups. It expresses the mean
  difference in terms of standard deviation units.

  The most common formula for Cohen's d is:

      d = (mean(group1) - mean(group2)) / pooled_stddev

  where `pooled_stddev` is the pooled standard deviation of the two groups,
  calculated under the assumption of equal variances.

  Parameters:

  - `group1` (seq of numbers): The first independent sample.
  - `group2` (seq of numbers): The second independent sample.
  - `method` (optional keyword): Specifies the method for calculating the pooled standard deviation,
    affecting the denominator of the formula. Possible values are `:unbiased` (default),
    `:biased`, or `:avg`. See [[pooled-stddev]] for details on these methods.

  Returns the calculated Cohen's d effect size as a double.

  Interpretation guidelines (approximate for normal distributions):
  - |d| = 0.2: small effect
  - |d| = 0.5: medium effect
  - |d| = 0.8: large effect

  Assumptions:
  - The two samples are independent.
  - Data within each group are approximately normally distributed.
  - The choice of `:method` implies assumptions about equal variances (default `:unbiased` and `:biased` assume equal variances, while `:avg` does not but might be less standard).

  See also [[hedges-g]] (a version bias-corrected for small sample sizes),
  [[glass-delta]] (an alternative effect size measure using the control group standard deviation),
  [[pooled-stddev]]."
  (^double [[group1 group2]] (cohens-d group1 group2))
  (^double [group1 group2] (cohens-d group1 group2 :unbiased))
  (^double [group1 group2 method]
   (let [group1 (m/seq->double-array group1)
         group2 (m/seq->double-array group2)
         diff (m/- (mean group1) (mean group2))]
     (m// diff (pooled-stddev [group1 group2] method)))))

(defn- effect-size-correction
  ^double [^long df]
  (m/- 1.0 (m// 3.0 (m/dec (m/* 4.0 df)))))

(defn cohens-d-corrected
  "Calculates Cohen's d effect size corrected for bias in small sample sizes.

  This function applies a correction factor (derived from the gamma function) to
  Cohen's d ([[cohens-d]]) to provide a less biased estimate of the population
  effect size when sample sizes are small. This corrected measure is sometimes
  referred to as Hedges' g, though this function specifically implements the
  correction applied to Cohen's d.

  The correction factor is `(1 - 3 / (4 * df - 1))` where `df` is the degrees of
  freedom used in the standard Cohen's d calculation.

  Parameters:

  - `group1` (seq of numbers): The first independent sample.
  - `group2` (seq of numbers): The second independent sample.
  - `method` (optional keyword): Specifies the method for calculating the pooled
    standard deviation, affecting the denominator of the formula (passed to
    [[cohens-d]]). Possible values are `:unbiased` (default), `:biased`, or `:avg`.
    See [[pooled-stddev]] for details on these methods.

  Returns the calculated bias-corrected Cohen's d effect size as a double.

  Note: While this function is named `cohens-d-corrected`, Hedges' g (calculated
  by [[hedges-g-corrected]]) also applies a similar small-sample bias correction.
  Differences might exist based on the specific correction formula or degree of
  freedom definition used. This function uses `(count group1) + (count group2) - 2`
  as the degrees of freedom for the correction by default (when `:unbiased` method
  is used for `cohens-d`).

  See also [[cohens-d]], [[hedges-g]], [[hedges-g-corrected]]."
  (^double [[group1 group2]] (cohens-d-corrected group1 group2))
  (^double [group1 group2] (cohens-d-corrected group1 group2 :unbiased))
  (^double [group1 group2 method]
   (m/* (effect-size-correction (long (if (= method :biased)
                                        (m/+ (count group1)
                                             (count group2))
                                        (m/+ (count group1)
                                             (count group2) -2))))
        (cohens-d group1 group2 method))))

(defn hedges-g
  "Calculates Hedges's g effect size for comparing the means of two independent groups.

  Hedges's g is a standardized measure quantifying the magnitude of the difference
  between the means of two independent groups. It is similar to Cohen's d but
  uses the *unbiased* pooled standard deviation in the denominator.

  This implementation calculates g using the unbiased pooled standard deviation as the denominator.

  Parameters:

  - `group1`, `group2` (sequences): The two independent samples directly as arguments.

  Returns the calculated Hedges's g effect size as a double.

  Note: This specific function uses the unbiased pooled standard deviation but does
  *not* apply the small-sample bias correction factor (often denoted as J)
  sometimes associated with Hedges's g. For a bias-corrected version, see [[hedges-g-corrected]].
  This function is equivalent to calling `(cohens-d group1 group2 :unbiased)`.

  See also [[cohens-d]], [[hedges-g-corrected]], [[glass-delta]], [[pooled-stddev]]."
  (^double [[group1 group2]] (hedges-g group1 group2))
  (^double [group1 group2]
   (cohens-d group1 group2 :unbiased)))

(defn hedges-g-corrected
  "Calculates a small-sample bias-corrected effect size for comparing the means
  of two independent groups, often referred to as a form of Hedges's g.

  This function calculates Cohen's d ([[cohens-d]]) using the *unbiased*
  pooled standard deviation (equivalent to [[hedges-g]]), and then applies
  a specific correction factor designed to reduce the bias in the effect size
  estimate for small sample sizes.

  The correction factor applied is `(1 - 3 / (4 * df - 1))`, where `df` is the
  degrees of freedom for the unbiased pooled variance calculation (`n1 + n2 - 2`).
  This corresponds to calling [[cohens-d-corrected]] with the `:unbiased` method
  for pooled standard deviation.

  Parameters:

  - `group1` (seq of numbers): The first independent sample.
  - `group2` (seq of numbers): The second independent sample.

  Returns the calculated bias-corrected effect size as a double.

  Note: This function applies *a* correction factor. For the more
  standard Hedges's g bias correction using the exact gamma function
  based correction factor, see [[hedges-g*]].

  See also [[cohens-d]], [[cohens-d-corrected]], [[hedges-g]], [[hedges-g*]],
  [[pooled-stddev]]."
  (^double [[group1 group2]] (hedges-g-corrected group1 group2))
  (^double [group1 group2]
   (cohens-d-corrected group1 group2 :unbiased)))

(defn hedges-g*
  "Calculates a less biased estimate of Hedges's g effect size for comparing the means of two independent groups, using the exact J bias correction.

  Hedges's g is a standardized measure of the difference between two means. For small sample sizes, the standard Hedges's g (and Cohen's d) can overestimate the true population effect size. This function applies a specific correction factor, often denoted as J, to mitigate this bias.

  The calculation involves:
  1. Calculating the standard Hedges's g (equivalent to [[hedges-g]], which uses the unbiased pooled standard deviation).
  2. Calculating the J correction factor based on the degrees of freedom (`n1 + n2 - 2`) using the gamma function.
  3. Multiplying the standard Hedges's g by the J factor.

  The J factor is calculated as `(Gamma(df/2) / (sqrt(df/2) * Gamma((df-1)/2)))`.

  Parameters:

  - `group1` (seq of numbers): The first independent sample.
  - `group2` (seq of numbers): The second independent sample.

  Returns the calculated bias-corrected Hedges's g effect size as a double.

  This version of Hedges's g is generally preferred over the standard version or Cohen's d when working with small sample sizes, as it provides a more accurate estimate of the population effect size.

  Assumptions:
  - The two samples are independent.
  - Data within each group are approximately normally distributed.
  - Equal variances are assumed for calculating the pooled standard deviation.

  See also [[cohens-d]], [[hedges-g]] (uncorrected), [[hedges-g-corrected]] (another correction method)."
  (^double [[group1 group2]] (hedges-g* group1 group2))
  (^double [group1 group2]
   (let [df (m/+ (count group1) (count group2) -2)
         df2 (m/* 0.5 df)
         j (m/exp (m/- (special/log-gamma df2)
                     (m/log (m/sqrt df2))
                     (special/log-gamma (m/* 0.5 (m/dec df)))))]
     (m/* j (hedges-g group1 group2)))))

(defn glass-delta
  "Calculates Glass's delta (Δ), an effect size measure for the difference
  between two group means, using the standard deviation of the control group.

  Glass's delta is used to quantify the magnitude of the difference between an
  experimental group and a control group, specifically when the control group's
  standard deviation is considered a better estimate of the population
  standard deviation than a pooled variance.

  Parameters:

  - `group1` (seq of numbers): The experimental group.
  - `group2` (seq of numbers): The control group.

  Returns the calculated Glass's delta as a double.

  This measure is less common than [[cohens-d]] or [[hedges-g]] but is preferred
  when the intervention is expected to affect the variance or when group2 (the control)
  is clearly the baseline against which variability should be assessed.

  See also [[cohens-d]], [[hedges-g]]."
  (^double [[group1 group2]] (glass-delta group1 group2))
  (^double [group1 group2]
   (let [group2 (m/seq->double-array group2)]
     (m// (m/- (mean group1) (mean group2)) (stddev group2)))))

(defn means-ratio
  "Calculates the ratio of the mean of `group1` to the mean of `group2`.

  This is a measure of effect size in the 'Ratio Family', comparing the central tendency
  of two groups multiplicatively.

  Parameters:

  - `group1` (seq of numbers): The first independent sample. The mean of this group is the numerator.
  - `group2` (seq of numbers): The second independent sample. The mean of this group is the denominator.
  - `adjusted?` (boolean, optional): If `true`, applies a small-sample bias correction to the ratio.
    Defaults to `false`.

  Returns the calculated ratio of means as a double.

  A value greater than 1 indicates that `group1` has a larger mean than `group2`.
  A value less than 1 indicates `group1` has a smaller mean.
  A value close to 1 indicates similar means.

  The `adjusted?` version attempts to provide a less biased estimate of the population
  mean ratio, particularly for small sample sizes, by incorporating variances into the calculation
  (based on Bickel and Doksum, see also [[means-ratio-corrected]]).

  See also [[means-ratio-corrected]] (which is equivalent to calling this with `adjusted?` set to `true`)."
  (^double [[group1 group2]] (means-ratio group1 group2))
  (^double [group1 group2] (means-ratio group1 group2 false))
  (^double [group1 group2 adjusted?]
   (let [ag1 (m/seq->double-array group1)
         ag2 (m/seq->double-array group2)
         m1 (mean ag1)
         m2 (mean ag2)]
     (if-not adjusted?
       (m// m1 m2)
       (let [v1 (variance ag1)
             v2 (variance ag2)
             n1 (alength ag1)
             n2 (alength ag2)
             J (m/* 0.5 (m/- (m// v1 (m/* n1 m1 m1))
                         (m// v2 (m/* n2 m2 m2))))]
         (m/exp (m/+ (m/- (m/log m1) (m/log m2)) J)))))))

(defn means-ratio-corrected
  "Calculates a bias-corrected ratio of the mean of `group1` to the mean of `group2`.

  This function applies a correction (based on Bickel and Doksum) to the simple
  ratio `mean(group1) / mean(group2)` to reduce bias, particularly for small
  sample sizes.

  It is equivalent to calling `(means-ratio group1 group2 true)`.

  Parameters:

  - `group1` (seq of numbers): The first independent sample. The mean of this group
    is the numerator.
  - `group2` (seq of numbers): The second independent sample. The mean of this group
    is the denominator.

  Returns the calculated bias-corrected ratio of means as a double.

  See also [[means-ratio]] (for the simple, uncorrected ratio)."
  (^double [[group1 group2]] (means-ratio-corrected group1 group2))
  (^double [group1 group2]
   (means-ratio group1 group2 true)))

(defn cliffs-delta
  "Calculates Cliff's Delta (δ), a non-parametric effect size measure for assessing the difference between two groups of ordinal or continuous data.

  Cliff's Delta quantifies the degree of overlap between two distributions. It represents the probability that a randomly chosen value from the first group is greater than a randomly chosen value from the second group, minus the reverse probability.

  Parameters:

  - `group1` (seq of numbers): The first sample.
  - `group2` (seq of numbers): The second sample.

  Returns the calculated Cliff's Delta value as a double.

  Interpretation:

  - A value of +1 indicates complete separation where every value in `group1` is greater than every value in `group2`.
  - A value of -1 indicates complete separation where every value in `group2` is greater than every value in `group1`.
  - A value of 0 indicates complete overlap between the distributions.
  - Values between -1 and 1 indicate varying degrees of overlap. Cohen (1988) suggested guidelines for effect size: |δ| < 0.147 (negligible), 0.147 ≤ |δ| < 0.33 (small), 0.33 ≤ |δ| < 0.474 (medium), |δ| ≥ 0.474 (large).

  Cliff's Delta is a robust measure, suitable for ordinal data or when assumptions of parametric tests (like normality or equal variances) are violated. It is closely related to the [[wmw-odds]] (Wilcoxon-Mann-Whitney odds) and the [[ameasure]] (Vargha-Delaney A).

  See also [[wmw-odds]], [[ameasure]], [[cohens-d]], [[glass-delta]]."
  (^double [[group1 group2]] (cliffs-delta group1 group2))
  (^double [group1 group2]
   (m// (sum (for [a group1
                 b group2]
             (m/signum (compare a b))))
      (m/* (count group1) (count group2)))))

;;

(defn ameasure
  "Calculates the Vargha-Delaney A measure for two independent samples.

  A non-parametric effect size measure quantifying the probability that a randomly chosen value from the first sample (`group1`) is greater than a randomly chosen value from the second sample (`group2`).

  Parameters:

  - `group1`: The first independent sample.
  - `group2`: The second independent sample.

  Returns the calculated A measure (a double) in the range [0, 1].
  A value of 0.5 indicates stochastic equality (distributions are overlapping). Values > 0.5 mean `group1` tends to be larger; values < 0.5 mean `group2` tends to be larger.

  Related to [[cliffs-delta]] and the Wilcoxon-Mann-Whitney U test statistic.

  See also [[cliffs-delta]], [[wmw-odds]]."
  (^double [[group1 group2]] (ameasure group1 group2))
  (^double [group1 group2]
   (let [m (count group1)
         n (count group2)
         r1 (sum (take m (m/rank1 (concat group1 group2))))]
     (m// (m/- (m/+ r1 r1) (m/* m (m/inc m)))
        (m/* 2.0 m n)))))

(defn wmw-odds
  "Calculates the Wilcoxon-Mann-Whitney odds (often denoted as ψ) for two independent samples.

  This non-parametric effect size measure quantifies the odds that a randomly chosen
  observation from the first group (`group1`) is greater than a randomly chosen
  observation from the second group (`group2`).

  The statistic is directly related to [[cliffs-delta]] (δ): ψ = (1 + δ) / (1 - δ).

  Parameters:

  - `group1` (seq of numbers): The first independent sample.
  - `group2` (seq of numbers): The second independent sample.

  Returns the calculated WMW odds as a double.

  Interpretation:

  - A value greater than 1 indicates that values from `group1` tend to be larger than values from `group2`.
  - A value less than 1 indicates that values from `group1` tend to be smaller than values from `group2`.
  - A value of 1 indicates stochastic equality between the distributions (50/50 odds).

  This measure is robust to violations of normality and is suitable for ordinal data.
  It is closely related to Cliff's Delta (δ) and the Mann-Whitney U test statistic.

  See also [[cliffs-delta]], [[ameasure]]."
  (^double [[group1 group2]] (wmw-odds group1 group2))
  (^double [group1 group2]
   (m/exp (m/logit (m/* 0.5 (m/inc (cliffs-delta group1 group2)))))))

(defn- integrate-kde
  [iterations kde ranges]
  (let [^RombergIntegrator integrator
        (RombergIntegrator. iterations RombergIntegrator/ROMBERG_MAX_ITERATIONS_COUNT)
        
        f (reify UnivariateFunction (value [_ x] (kde x)))]
    (map (fn [[^double mn ^double mx]]
           (.integrate integrator Integer/MAX_VALUE f mn mx)) ranges)))

(defn p-overlap
  "Calculates the overlapping index between the estimated distributions of two samples using Kernel Density Estimation (KDE).

  This function estimates the probability density function (PDF) for `group1` and `group2` using KDE and then calculates the area of overlap between the two estimated PDFs. The area of overlap is the integral of the minimum of the two density functions.

  Parameters:

  - `group1` (seq of numbers): The first sample.
  - `group2` (seq of numbers): The second sample.
  - `opts` (map, optional): Options map for KDE and integration:
    - `:kde` (keyword, default `:gaussian`): The kernel function to use for KDE. See `fastmath.kernel.density/kernel-density+` for options.
    - `:bandwidth` (double, optional): The bandwidth for KDE. If omitted, it is automatically estimated.
    - `:min-iterations` (long, default 3): Minimum number of iterations for Romberg integration.
    - `:steps` (long, default 500): Number of steps (subintervals) for numerical integration over the relevant range.

  Returns the calculated overlapping index as a double, representing the area of overlap between the two estimated distributions. A value closer to 1 indicates greater overlap, while a value closer to 0 indicates less overlap.

  This measure quantifies the degree to which two distributions share common values and can be seen as a measure of similarity."
  (^double [[group1 group2]] (p-overlap group1 group2))
  (^double [group1 group2] (p-overlap group1 group2 {}))
  (^double [group1 group2 {:keys [kde bandwidth ^long min-iterations ^long steps]
                           :or {kde :gaussian min-iterations 3 steps 500}}]
   (let [{kde1 :kde ^double h1 :h
          ^double mn1 :mn ^double mx1 :mx} (kd/kernel-density+ kde group1 {:bandwidth bandwidth})
         {kde2 :kde ^double h2 :h
          ^double mn2 :mn ^double mx2 :mx} (kd/kernel-density+ kde group2 {:bandwidth bandwidth})
         h (m/* 2.0 (m/+ h1 h2))
         mn (m/- (m/min mn1 mn2) h)
         mx (m/+ (m/max mx1 mx2) h)
         iters (m/max 2 min-iterations)
         ranges (partition 2 1 (m/slice-range mn mx steps))
         i1 (integrate-kde iters kde1 ranges)
         i2 (integrate-kde iters kde2 ranges)]
     (sum (map m/min i1 i2)))))

;; https://www.psy.gla.ac.uk/~steve/best/effectsize.ppt.pdf

(defn cohens-u1-normal
  "Calculates Cohen's U1, a measure of non-overlap between two distributions assumed to be normal with equal variances.

  Cohen's U1 quantifies the proportion of scores in the lower-scoring group that overlap
  with the scores in the higher-scoring group. A U1 of 0 means no overlap, while
  a U1 of 1 means complete overlap (distributions are identical).

  This measure is calculated directly from Cohen's d statistic ([[cohens-d]]) assuming
  normal distributions and equal variances.

  Parameters:

  - `group1` (seq of numbers): The first sample.
  - `group2` (seq of numbers): The second sample.
  - `method` (optional keyword): Specifies the method for calculating the pooled standard deviation
    used in the underlying [[cohens-d]] calculation. Possible values are `:unbiased` (default),
    `:biased`, or `:avg`. See [[pooled-stddev]] for details.
  - `d` (double): A pre-calculated Cohen's d value. If provided, `group1`, `group2`, and `method` are ignored.

  Returns the calculated Cohen's U1 as a double [0, 1].

  Assumptions:
  - Both samples are drawn from normally distributed populations.
  - The populations have equal variances (homoscedasticity).

  See also [[cohens-d]], [[cohens-u2-normal]], [[cohens-u3-normal]], [[p-overlap]] (a non-parametric overlap measure)."
  (^double [group1 group2] (cohens-u1-normal group1 group2 :unbiased))
  (^double [group1 group2 method] (cohens-u1-normal (cohens-d group1 group2 method)))
  (^double [^double d]
   (let [p (r/cdf r/default-normal (m/* 0.5 (m/abs d)))]
     (m// (m/dec (m/* 2.0 p)) p))))

(defn cohens-u2-normal
  "Calculates Cohen's U2, a measure of overlap between two distributions assumed to be normal with equal variances.

  Cohen's U2 quantifies the proportion of scores in the lower-scoring group that are below the point located halfway between the means of the two groups (or equivalently, the proportion of scores in the higher-scoring group that are above this halfway point). This measure is calculated from Cohen's d statistic ([[cohens-d]]) using the standard normal cumulative distribution function ($\\Phi$): $\\Phi(0.5 |d|)$.

  Parameters:

  - `group1` (seq of numbers): The first sample.
  - `group2` (seq of numbers): The second sample.
  - `method` (optional keyword): Specifies the method for calculating the pooled standard deviation used in the underlying [[cohens-d]] calculation. Possible values are `:unbiased` (default), `:biased`, or `:avg`. See [[pooled-stddev]] for details.
  - `d` (double): A pre-calculated Cohen's d value. If provided, `group1`, `group2`, and `method` are ignored.

  Returns the calculated Cohen's U2 as a double [0.0, 1.0]. A value closer to 0.5 indicates greater overlap between the distributions; values closer to 0 or 1 indicate less overlap.

  Assumptions:
  - Both samples are drawn from normally distributed populations.
  - The populations have equal variances (homoscedasticity).

  See also [[cohens-d]], [[cohens-u1-normal]], [[cohens-u3-normal]], [[p-overlap]] (a non-parametric overlap measure)."
  (^double [group1 group2] (cohens-u2-normal group1 group2 :unbiased))
  (^double [group1 group2 method] (cohens-u2-normal (cohens-d group1 group2 method)))
  (^double [^double d] (r/cdf r/default-normal (m/* 0.5 (m/abs d)))))

(defn cohens-u3-normal
  "Calculates Cohen's U3, a measure of overlap between two distributions assumed to be normal with equal variances.

  Cohen's U3 quantifies the proportion of scores in the lower-scoring group that fall
  below the mean of the higher-scoring group. It is calculated from Cohen's d statistic
  ([[cohens-d]]) using the standard normal cumulative distribution function ($\\Phi$):
  `U3 = Φ(d)`.

  The measure is asymmetric: `U3(group1, group2)` is not necessarily equal to
  `U3(group2, group1)`. The interpretation depends on which group is considered
  the 'higher-scoring' one based on the sign of d. By convention, the result
  often represents the proportion of the *first* group (`group1`) that is below the
  mean of the *second* group (`group2`) if d is negative, or the proportion of the
  *second* group (`group2`) that is below the mean of the *first* group (`group1`) if d is positive.

  Parameters:

  - `group1` (seq of numbers): The first sample.
  - `group2` (seq of numbers): The second sample.
  - `method` (optional keyword): Specifies the method for calculating the pooled standard deviation
    used in the underlying [[cohens-d]] calculation. Possible values are `:unbiased` (default),
    `:biased`, or `:avg`. See [[pooled-stddev]] for details.
  - `d` (double): A pre-calculated Cohen's d value. If provided, `group1`, `group2`, and `method` are ignored.

  Returns the calculated Cohen's U3 as a double [0.0, 1.0].
  A value close to 0.5 suggests significant overlap. Values closer to 0 or 1 suggest
  less overlap (greater separation between the means).

  Assumptions:
  - Both samples are drawn from normally distributed populations.
  - The populations have equal variances (homoscedasticity).

  See also [[cohens-d]], [[cohens-u1-normal]], [[cohens-u2-normal]], [[p-overlap]] (a non-parametric overlap measure)."
  (^double [group1 group2] (cohens-u3-normal group1 group2 :unbiased))
  (^double [group1 group2 method] (cohens-u3-normal (cohens-d group1 group2 method)))
  (^double [^double d] (r/cdf r/default-normal d)))

(defn cohens-u2
  "Calculates a measure of overlap between two samples, referred to as Cohen's U2.

  This function quantifies the degree to which the distributions of `group1` and `group2` overlap. It is related to comparing values at corresponding percentile levels across the two groups or the proportion of values in one group that are below the median of the other. A value of 0 indicates no overlap, while a value of 1 indicates complete overlap (distributions are identical).

  The measure is symmetric, meaning `(cohens-u2 group1 group2)` is equal to `(cohens-u2 group2 group1)`.

  This is a non-parametric measure, suitable for any data samples, and does not assume normality, unlike [[cohens-u2-normal]].

  Parameters:

  - `group1`, `group2` (sequences): The two samples directly as arguments.

  Returns the calculated Cohen's U2 value as a double. The value typically ranges from 0 to 1. A value closer to 0.5 indicates substantial overlap between the distributions (e.g., the median of one group is near the median of the other); values closer to 0 or 1 indicate less overlap (greater separation between the distributions)."
  (^double [[group1 group2]] (cohens-u2 group1 group2))
  (^double [group1 group2]
   (let [g1 (r/distribution :real-discrete-distribution {:data group1})
         g2 (r/distribution :real-discrete-distribution {:data group2})
         target-fn (fn [^double p]
                     (let [p (m/constrain p 4.9E-324 0.9999999999999999)
                           p- (m/- 1.0 p)
                           q1 (Vec2. (r/icdf g1 p) (r/icdf g1 p-))
                           q2 (Vec2. (r/icdf g2 p-) (r/icdf g2 p))]
                       (-> (v/sub q1 q2) v/abs v/mn)))]
     (-> (opt/minimize :brent target-fn {:bounds [[0.5 1.0]]}) ffirst double))))

(defn cohens-u1
  "Calculates a non-parametric measure of difference or separation between two samples.

  This function computes a value derived from [[cohens-u2]], which internally
  quantifies a minimal difference between corresponding quantiles of the two
  empirical distributions.

  Parameters:

  - `group1` (seq of numbers): The first sample.
  - `group2` (seq of numbers): The second sample.

  Returns the calculated measure as a double.

  Interpretation:

  - Values close to -1 indicate high similarity or maximum overlap between the
    distributions (as the minimal difference between quantiles approaches zero).
  - Increasing values indicate greater difference or separation between the
    distributions (as the minimal difference between quantiles is larger).

  This measure is symmetric, meaning the order of `group1` and `group2` does not
  affect the result. It is a non-parametric measure applicable to any data samples.

  See also [[cohens-u2]] (the measure this calculation is based on),
  [[cohens-u3]] (related non-parametric measure), [[cohens-u1-normal]]
  (the version applicable to normal data)."
  (^double [[group1 group2]] (cohens-u1 group1 group2))
  (^double [group1 group2]
   (let [u2 (cohens-u2 group1 group2)]
     (m// (m/dec (m/* 2.0 u2)) u2))))

(defn cohens-u3
  "Calculates Cohen's U3 for two samples.

  In this implementation, Cohen's U3 is defined as the proportion of values
  in the second sample (`group2`) that are less than the median of the first
  sample (`group1`).

  Parameters:

  - `group1` (seq of numbers): The first sample. The median of this sample is used as the threshold.
  - `group2` (seq of numbers): The second sample. Values from this sample are counted if they fall below the median of `group1`.
  - `estimation-strategy` (optional keyword): The strategy used to estimate the median of `group1`.
    Defaults to `:legacy`. See [[median]] or [[quantile]] for available strategies
    (e.g., `:r1` through `:r9`).

  Returns the calculated proportion as a double between 0.0 and 1.0.

  Interpretation:

  - A value close to 0 means most values in `group2` are greater than or equal to the median of `group1`.
  - A value close to 0.5 means approximately half the values in `group2` are below the median of `group1`.
  - A value close to 1 means most values in `group2` are less than the median of `group1`.

  Note: This measure is **not symmetric**. `(cohens-u3 group1 group2)` is generally
  not equal to `(cohens-u3 group2 group1)`.

  This is a non-parametric measure, suitable for any data samples, and does not
  assume normality, unlike [[cohens-u3-normal]].

  See also [[cohens-u3-normal]] (the version applicable to normal data), [[cohens-u2]]
  (a related symmetric non-parametric measure), [[median]], [[quantile]]."
  (^double [[group1 group2]] (cohens-u3 group1 group2))
  (^double [group1 group2] (cohens-u3 group1 group2 :legacy))
  (^double [group1 group2 estimation-strategy]
   (let [m (median group1 estimation-strategy)]
     (m// (count (filter (fn [^double v] (m/< v m)) group2)) (double (count group2))))))

;;

(defn pearson-r
  "Calculates the Pearson `r` correlation coefficient between two sequences.

  This function is an alias for [[pearson-correlation]].

  See [[pearson-correlation]] for detailed documentation, parameters, and usage examples."
  (^double [[group1 group2]] (pearson-r group1 group2))
  (^double [group1 group2]
   (pearson-correlation group1 group2)))

(defn r2-determination
  "Calculates the Coefficient of Determination ($R^2$) between two sequences.

  This function computes the square of the Pearson product-moment correlation
  coefficient ([[pearson-correlation]]) between `group1` and `group2`.

  $R^2$ measures the proportion of the variance in one variable that is predictable
  from the other variable in a linear relationship. For a simple linear regression
  with one independent variable, this value is equivalent to the $R^2$ calculated
  from the Residual Sum of Squares (RSS) and Total Sum of Squares (TSS).

  Parameters:

  - `group1` (seq of numbers): The first sequence.
  - `group2` (seq of numbers): The second sequence.

  Both sequences must have the same length.

  Returns the calculated $R^2$ value (a double between 0.0 and 1.0) as a double.
  Returns `NaN` if the Pearson correlation cannot be calculated (e.g., one sequence is constant).

  See also [[r2]] (for general $R^2$ and adjusted $R^2$), [[pearson-correlation]]."
  (^double [[group1 group2]] (r2-determination group1 group2))
  (^double [group1 group2]
   (m/sq (pearson-correlation group1 group2))))

(defn- local-linear-regression
  ^SimpleRegression [group1 group2]
  (let [lm (SimpleRegression. true)]
    (.addData lm (m/seq->double-double-array (map vector group1 group2)))
    (.regress lm)
    lm))

(defn eta-sq
  "Calculates a measure of association between two sequences, named `eta-sq` (Eta-squared).

  *Note*: The current implementation calculates the R-squared coefficient of determination from a simple linear regression where the first input sequence (`group1`) is treated as the dependent variable and the second (`group2`) as the independent variable. In this context, it quantifies the proportion of the variance in `group1` that is linearly predictable from `group2`.

  Parameters:

  - `group1` (seq of numbers): The first sequence (treated as dependent variable).
  - `group2` (seq of numbers): The second sequence (treated as independent variable).

  Returns the calculated R-squared value as a double [0.0, 1.0].

  Interpretation:

  - 0.0 indicates that `group2` explains none of the variance in `group1` linearly.
  - 1.0 indicates that `group2` linearly explains all the variance in `group1`.

  While Eta-squared ($\\eta^2$) is commonly used in ANOVA to quantify the proportion of variance in a dependent variable explained by group membership, this function's calculation method differs from the standard ANOVA $\\eta^2$ unless `group2` explicitly represents numeric codes for two groups.

  See also [[r2-determination]] (which is equivalent to this function), [[pearson-correlation]], [[omega-sq]], [[epsilon-sq]], [[one-way-anova-test]]."
  (^double [[group1 group2]] (eta-sq group1 group2))
  (^double [group1 group2]
   (.getRSquare (local-linear-regression group1 group2))))

(defn omega-sq
  "Calculates Omega squared (ω²), an effect size measure for the simple linear regression of `group1` on `group2`.

  Omega squared estimates the proportion of variance in the dependent variable (`group1`) that is accounted for by the independent variable (`group2`) in the population. It is considered a less biased alternative to [[r2-determination]].

  Parameters:

  - `group1` (seq of numbers): The dependent variable.
  - `group2` (seq of numbers): The independent variable. Must have the same length as `group1`.
  - `degrees-of-freedom` (double, optional): The degrees of freedom for the regression model. Defaults to 1.0, which is standard for simple linear regression and used in the 2-arity version. Providing a different value allows calculating ω² for cases with multiple predictors if the sums of squares are computed for the overall model.

  Returns the calculated Omega squared value as a double. The value typically ranges from 0.0 to 1.0.

  Interpretation:

  - 0.0 indicates that `group2` explains none of the variance in `group1` in the population.
  - 1.0 indicates that `group2` perfectly explains the variance in `group1` in the population.

  Note: While often presented in the context of ANOVA, this implementation applies the formula to the sums of squares obtained from a simple linear regression between the two sequences. The 3-arity version allows specifying a custom degrees of freedom for regression, which might be relevant for calculating overall $\\omega^2$ in multiple regression contexts (where `degrees-of-freedom` would be the number of predictors).

  See also [[eta-sq]] (Eta-squared, often based on $R^2$), [[epsilon-sq]] (another adjusted R²-like measure), [[r2-determination]] (R-squared)."
  (^double [[group1 group2]] (omega-sq group1 group2))
  (^double [group1 group2]
   (let [lm (local-linear-regression group1 group2)
         mse (.getMeanSquareError lm)]
     (m// (m/- (.getRegressionSumSquares lm) mse)
        (m/+ (.getTotalSumSquares lm) mse))))
  (^double [group1 group2 ^double degrees-of-freedom]
   (let [lm (local-linear-regression group1 group2)
         mse (.getMeanSquareError lm)]
     (m// (m/- (.getRegressionSumSquares lm) (m/* degrees-of-freedom mse))
        (m/+ (.getTotalSumSquares lm) mse)))))

(defn epsilon-sq
  "Calculates Epsilon squared (ε²), an effect size measure for the simple linear regression of `group1` on `group2`.

  Epsilon squared estimates the proportion of variance in the dependent variable (`group1`)
  that is accounted for by the independent variable (`group2`) in the population. It is
  considered a less biased alternative to the sample R-squared ([[r2-determination]]).

  The calculation is based on the sums of squares from the simple linear regression of
  `group1` on `group2`.

  Parameters:

  - `group1` (seq of numbers): The dependent variable.
  - `group2` (seq of numbers): The independent variable. Must have the same length as `group1`.

  Returns the calculated Epsilon squared value as a double. The value typically ranges
  from 0.0 to 1.0.

  Interpretation:

  - 0.0 indicates that `group2` explains none of the variance in `group1` in the population.
  - 1.0 indicates that `group2` perfectly explains the variance in `group1` in the population.

  Note: While often presented in the context of ANOVA, this implementation applies the
  formula to the sums of squares obtained from a simple linear regression between the
  two sequences.

  See also [[eta-sq]] (Eta-squared, often based on $R^2$), [[omega-sq]] (another adjusted
  R²-like measure), [[r2-determination]] (R-squared)."
  (^double [[group1 group2]] (epsilon-sq group1 group2))
  (^double [group1 group2]
   (let [lm (local-linear-regression group1 group2)]
     (m// (m/- (.getRegressionSumSquares lm) (.getMeanSquareError lm))
        (.getTotalSumSquares lm)))))

(defn cohens-f2
  "Calculates Cohen's f², a measure of effect size often used in ANOVA or regression.

  Cohen's f² quantifies the magnitude of the effect of an independent variable or set
  of predictors on a dependent variable, expressed as the ratio of the variance
  explained by the effect to the unexplained variance.

  This function allows calculating f² using different measures for the 'Proportion of Variance Explained',
  specified by the `type` parameter:

  - `:eta` (default): Uses [[eta-sq]] (Eta-squared), which in this implementation is
    equivalent to the sample $R^2$ from a linear regression of `group1` on `group2`.
    This is a measure of the proportion of variance explained in the sample.
  - `:omega`: Uses [[omega-sq]] (Omega-squared), a less biased estimate of the
    proportion of variance explained in the population.
  - `:epsilon`: Uses [[epsilon-sq]] (Epsilon-squared), another less biased estimate
    of the proportion of variance explained in the population, similar to adjusted $R^2$.
  - Any function: A function accepting `group1` and `group2` and returning a double representing the proportion of variance explained.

  Parameters:

  - `group1` (seq of numbers): The dependent variable.
  - `group2` (seq of numbers): The independent variable (or predictor). Must have the same length as `group1`.
  - `type` (keyword, optional): Specifies the measure of 'Proportion of Variance Explained' to use (`:eta`, `:omega`, `:epsilon` or any function). Defaults to `:eta`.

  Returns the calculated Cohen's f² effect size as a double. Values range from 0 upwards.

  Interpretation Guidelines (approximate, often used for F-tests in ANOVA/regression):
  - $f^2 = 0.02$: small effect
  - $f^2 = 0.15$: medium effect
  - $f^2 = 0.35$: large effect

  See also [[cohens-f]], [[eta-sq]], [[omega-sq]], [[epsilon-sq]]."
  (^double [[group1 group2]] (cohens-f2 group1 group2))
  (^double [group1 group2] (cohens-f2 group1 group2 :eta))
  (^double [group1 group2 type]
   (let [f (if (keyword? type) (case type
                                 :omega omega-sq
                                 :epsilon epsilon-sq
                                 eta-sq)
               type)
         v (double (f group1 group2))]
     (m// v (m/- 1.0 v)))))

(defn cohens-f
  "Calculates Cohen's f, a measure of effect size derived as the square root of Cohen's f² ([[cohens-f2]]).

  Cohen's f is a standardized measure quantifying the magnitude of an effect,
  often used in the context of ANOVA or regression. It is the square root of
  the ratio of the variance explained by the effect to the unexplained variance.

  Parameters:

  - `group1` (seq of numbers): The dependent variable.
  - `group2` (seq of numbers): The independent variable (or predictor). Must have the same length as `group1`.
  - `type` (keyword, optional): Specifies the measure of 'Proportion of Variance Explained'
    used in the underlying [[cohens-f2]] calculation. Defaults to `:eta`.
    - `:eta` (default): Uses Eta-squared (sample R²), a measure of variance explained in the sample.
    - `:omega`: Uses Omega-squared, a less biased estimate of variance explained in the population.
    - `:epsilon`: Uses Epsilon-squared, another less biased estimate of variance explained in the population.
    - Any function: A function accepting `group1` and `group2` and returning a double representing the proportion of variance explained.

  Returns the calculated Cohen's f effect size as a double. Values range from 0 upwards.

  Interpretation:

  - Values are positive. Larger values indicate a stronger effect (more variance in `group1` explained by `group2`).
  - Cohen's guidelines for interpreting the magnitude of f² (and by extension, f) are:
    - $f = 0.10$ (approx. $f^2 = 0.01$): small effect
    - $f = 0.25$ (approx. $f^2 = 0.0625$): medium effect
    - $f = 0.40$ (approx. $f^2 = 0.16$): large effect
    (Note: Guidelines are often quoted for f², interpret f as $\\sqrt{f^2}$)

  See also [[cohens-f2]], [[eta-sq]], [[omega-sq]], [[epsilon-sq]]."
  (^double [[group1 group2]] (cohens-f group1 group2))
  (^double [group1 group2] (cohens-f group1 group2 :eta))
  (^double [group1 group2 type] (m/sqrt (cohens-f2 group1 group2 type))))

(defn cohens-q
  "Compares two correlation coefficients by calculating the difference between their Fisher z-transformations.

  The Fisher z-transformation (`atanh`) of a correlation coefficient `r` helps normalize the sampling distribution of correlation coefficients. The difference between two z'-transformed correlations is often used as a test statistic.

  The function supports comparing correlations in different scenarios via its arities:

  - `(cohens-q r1 r2)`: Calculates the difference between the Fisher z-transformations of two correlation values `r1` and `r2` provided directly. This is typically used when comparing two *independent* correlation coefficients (e.g., correlations from two separate studies). Returns `atanh(r1) - atanh(r2)`.
    - `r1`, `r2` (double): Correlation coefficient values (-1.0 to 1.0).

  - `(cohens-q group1 group2a group2b)`: Calculates the difference between the correlation of `group1` with `group2a` and the correlation of `group1` with `group2b`. This is commonly used for comparing *dependent* correlations (where `group1` is a common variable). Calculates `atanh(pearson-correlation(group1, group2a)) - atanh(pearson-correlation(group1, group2b))`.
    - `group1`, `group2a`, `group2b` (sequences): Data sequences from which Pearson correlations are computed.

  - `(cohens-q group1a group2a group1b group2b)`: Calculates the difference between the correlation of `group1a` with `group2a` and the correlation of `group1b` with `group2b`. This is typically used for comparing two *independent* correlations obtained from two distinct pairs of variables (all four sequences are independent). Calculates `atanh(pearson-correlation(group1a, group2a)) - atanh(pearson-correlation(group1b, group2b))`.
    - `group1a`, `group2a`, `group1b`, `group2b` (sequences): Data sequences from which Pearson correlations are computed.

  Returns the difference between the Fisher z-transformed correlation values as a double.

  Note: For comparing dependent correlations (3-arity case), standard statistical tests (e.g., Steiger's test) are more complex than a simple difference of z-transforms and involve the correlation between `group2a` and `group2b`. This function provides the basic difference value."
  (^double [^double r1 ^double r2]
   (m/- (m/atanh r1) (m/atanh r2)))
  (^double [group1 group2a group2b]
   (cohens-q (pearson-correlation group1 group2a)
             (pearson-correlation group1 group2b)))
  (^double [group1a group2a group1b group2b]
   (cohens-q (pearson-correlation group1a group2a)
             (pearson-correlation group1b group2b))))

(declare kruskal-test)

(defn rank-eta-sq
  "Calculates the Rank Eta-squared (η²), an effect size measure for the Kruskal-Wallis H-test.

  Rank Eta-squared is a non-parametric measure quantifying the proportion of the
  total variability (based on ranks) in the dependent variable that is associated
  with group membership (the independent variable). It is analogous to Eta-squared
  in one-way ANOVA but used for the rank-based Kruskal-Wallis test.

  The statistic is calculated based on the Kruskal-Wallis H statistic, the number
  of groups (`k`), and the total number of observations (`n`).

  Parameters:

  - `xss` (sequence of sequences): A collection where each element is a sequence
    representing a group of observations, as used in [[kruskal-test]].

  Returns the calculated Rank Eta-squared value as a double, ranging from 0 to 1.

  Interpretation:
  - A value of 0 indicates no difference in the distributions across groups (all variability is within groups).
  - A value closer to 1 indicates that a large proportion of the variability is
    due to differences between group ranks.

  Rank Eta-squared is a useful supplement to the Kruskal-Wallis test, providing a
  measure of the magnitude of the group effect that is not sensitive to assumptions
  about the data distribution shape (beyond having similar shapes for valid
  interpretation of the Kruskal-Wallis test itself).

  See also [[kruskal-test]], [[rank-epsilon-sq]] (another rank-based effect size)."
  ^double [xss]
  (let [{:keys [^double stat ^long k ^long n]} (kruskal-test xss)]
    (m/max 0.0 (m// (m/inc (m/- stat k))
                    (m/- n k)))))

(defn rank-epsilon-sq
  "Calculates Rank Epsilon-squared (ε²), a measure of effect size for the Kruskal-Wallis H-test.

  Rank Epsilon-squared is a non-parametric measure quantifying the proportion of the
  total variability (based on ranks) in the dependent variable that is associated
  with group membership (the independent variable). It is analogous to Eta-squared
  or Epsilon-squared in one-way ANOVA but used for the rank-based Kruskal-Wallis test.

  This function calculates Epsilon-squared based on the Kruskal-Wallis H statistic (`H`)
  and the total number of observations (`n`) across all groups.

  Parameters:

  - `xss` (sequence of sequences): A collection where each element is a sequence
    representing a group of observations, as used in [[kruskal-test]].

  Returns the calculated Rank Epsilon-squared value as a double, ranging from 0 to 1.

  Interpretation:

  - A value of 0 indicates no difference in the distributions across groups.
  - A value closer to 1 indicates that a large proportion of the variability is
    due to differences between group ranks.

  Rank Epsilon-squared is a useful supplement to the Kruskal-Wallis test, providing a
  measure of the magnitude of the group effect that is not sensitive to assumptions
  about the data distribution shape (beyond having similar shapes for valid
  interpretation of the Kruskal-Wallis test itself).

  See also [[kruskal-test]], [[rank-eta-sq]] (another rank-based effect size)."
  ^double [xss]
  (let [{:keys [^double stat ^long n]} (kruskal-test xss)]
    (m// stat (m// (m/dec (m/* n n)) (m/inc n)))))

;;

(defn contingency-table
  "Creates a frequency map (contingency table) from one or more sequences.

  If one sequence `xs` is provided, it returns a simple frequency map of the values
  in `xs`.

  If multiple sequences `s1, s2, ..., sn` are provided, it creates a contingency table
  of the tuples formed by the corresponding elements `[s1_i, s2_i, ..., sn_i]` at
  each index `i`. The returned map keys are these tuples, and values are their
  frequencies.

  Parameters:

  - `seqs` (one or more sequences): The input sequences. All sequences should ideally
    have the same length, as elements are paired by index.

  Returns a map where keys represent unique combinations of values (or single values
  if only one sequence is input) and values are the counts of these combinations.

  See also [[rows->contingency-table]], [[contingency-table->marginals]]."
  [& seqs]
  (if (= 1 (count seqs))
    (frequencies (first seqs))
    (frequencies (apply map vector seqs))))

(defn rows->contingency-table
  "Converts a sequence of sequences (representing rows of counts) into a map-based contingency table.

  This function takes a collection where each inner sequence is treated as a row
  of counts in a grid or matrix. It transforms this matrix representation into a
  map where keys are `[row-index, column-index]` tuples and values are the
  non-zero counts at that intersection.

  This is particularly useful for converting structured count data, like the
  output of some grouping or tabulation processes, into a format suitable for
  functions expecting a contingency table map (like `contingency-table->marginals`
  or chi-squared tests).

  Parameters:

  - `xss` (sequence of sequences of numbers): A collection where each inner
    sequence `xs_i` contains counts for row `i`. Values within `xs_i` are
    interpreted as counts for columns `0, 1, ...`.

  Returns a map where keys are `[row-index, column-index]` vectors and values
  are the corresponding non-zero counts from the input matrix. Zero counts are
  omitted from the output map.

  See also [[contingency-table]] (for building tables from raw data), [[contingency-table->marginals]]."
  [xss]
  (->> (for [[row-id row] (map-indexed vector xss)
             [col-id ^long val] (map-indexed vector row)
             :when (not (m/zero? val))]
         [[row-id col-id] val])
       (into {})))

(defn- ct-marginals-sum
  [m]
  (->> (map (fn [[k v]]
              [k (sum (map second v))]) m)
       (sort-by first)))

(defn- infer-ct
  [ct]
  (if (map? ct) ct (rows->contingency-table ct)))

(defn contingency-table->marginals
  "Calculates marginal sums (row and column totals) and the grand total from a contingency table.

  A contingency table represents the frequency distribution of observations for two or
  more categorical variables. This function summarizes these frequencies along the
  rows and columns.

  The function accepts two main input formats for the contingency table:

  1.  A map where keys are `[row-index, column-index]` tuples and values are counts (e.g., `{[0 0] 10, [0 1] 5, [1 0] 3, [1 1] 12}`). This format is produced by [[contingency-table]] when given multiple sequences or by [[rows->contingency-table]].
  2.  A sequence of sequences representing the rows of the table, where each inner sequence contains counts for the columns in that row (e.g., `[[10 5] [3 12]]`). The function internally converts this format to the map format.

  Parameters:

  - `ct` (map or sequence of sequences): The contingency table input.

  Returns a map containing:

  - `:rows`: A sequence of `[row-index, row-total]` pairs.
  - `:cols`: A sequence of `[column-index, column-total]` pairs.
  - `:n`: The grand total of all counts in the table.
  - `:diag`: A sequence of `[[index, index], count]` pairs for cells on the diagonal
    (where row index equals column index). This is useful for square tables like
    confusion matrices.

  See also [[contingency-table]], [[rows->contingency-table]]."
  [ct]
  (let [ct (infer-ct ct)
        rows (ct-marginals-sum (group-by ffirst ct))
        cols (ct-marginals-sum (group-by (comp second first) ct))
        n (sum (vals ct))
        diag (filter (fn [[[a b]]] (= a b)) ct)]
    {:rows rows :cols cols :n n :diag diag}))

(defn mcc
  "Calculates the Matthews Correlation Coefficient (MCC), also known as the Phi coefficient,
  for a 2x2 contingency table or binary classification outcomes.

  MCC is a measure of the quality of binary classifications. It is a balanced
  measure which can be used even if the classes are of very different sizes.
  Its value ranges from -1 to +1.

  - A coefficient of +1 represents a perfect prediction.
  - 0 represents a prediction no better than random.
  - -1 represents a perfect inverse prediction.

  The function can be called in two ways:

  1.  With two sequences `group1` and `group2`:
      The function will automatically construct a 2x2 contingency table from
      the unique values in the sequences (assuming they represent two binary
      variables). The mapping of values to table cells (e.g., what corresponds
      to TP, TN, FP, FN) depends on how `contingency-table` orders the unique values.
      For direct control over which cell is which, use the contingency table input.

  2.  With a contingency table:
      The contingency table can be provided as:
      - A map where keys are `[row-index, column-index]` tuples and values are counts
        (e.g., `{[0 0] TP, [0 1] FP, [1 0] FN, [1 1] TN}`). This is the output format
        of [[contingency-table]] with two inputs.
      - A sequence of sequences representing the rows of the table
        (e.g., `[[TP FP] [FN TN]]`). This is equivalent to `rows->contingency-table`.

  Parameters:

  - `group1` (sequence): The first sequence of binary outcomes/categories.
  - `group2` (sequence): The second sequence of binary outcomes/categories.
    Must have the same length as `group1`.
  - `contingency-table` (map or sequence of sequences): A pre-computed 2x2 contingency table.

  Returns the calculated Matthews Correlation Coefficient as a double.

  Note: The implementation uses marginal sums from the contingency table, which
  is mathematically equivalent to the standard formula but avoids potential
  division by zero in the denominator product if any marginal sum is zero.

  See also [[contingency-table]], [[contingency-2x2-measures]], [[binary-measures-all]]."
  ([group1 group2] (mcc (contingency-table group1 group2)))
  ([ct]
   (let [{:keys [diag rows cols ^long n]} (contingency-table->marginals (infer-ct ct))
         t (map second rows)
         p (map second cols)
         d (map second diag)
         s2 (m/* n n)]
     (m// (m/- (m/* n (sum d)) (v/dot t p))
        (m/* (m/sqrt (m/- s2 (v/dot p p)))
           (m/sqrt (m/- s2 (v/dot t t))))))))

(declare chisq-test)

(defn cramers-c
  "Calculates Cramer's C, a measure of association (effect size) between two
  nominal variables represented in a contingency table.

  Its value ranges from 0 to 1, where 0 indicates no association and 1 indicates
  a perfect association. It is particularly useful for tables larger than 2x2.

  The function can be called in two ways:

  1.  With two sequences `group1` and `group2`:
      The function will automatically construct a contingency table from
      the unique values in the sequences.
  2.  With a contingency table:
      The contingency table can be provided as:
      - A map where keys are `[row-index, column-index]` tuples and values are counts
        (e.g., `{[0 0] 10, [0 1] 5, [1 0] 3, [1 1] 12}`). This is the output format
        of [[contingency-table]] with two inputs.
      - A sequence of sequences representing the rows of the table
        (e.g., `[[10 5] [3 12]]`). This is equivalent to `rows->contingency-table`.

  Parameters:

  - `group1` (sequence): The first sequence of categorical data.
  - `group2` (sequence): The second sequence of categorical data. Must have the same length as `group1`.
  - `contingency-table` (map or sequence of sequences): A pre-computed contingency table.

  Returns the calculated Cramer's C coefficient as a double.

  See also [[chisq-test]], [[cramers-v]], [[cohens-w]], [[tschuprows-t]], [[contingency-table]]."
  (^double [group1 group2] (cramers-c (contingency-table group1 group2)))
  (^double [contingency-table]
   (let [{:keys [^double chi2 ^long n]} (chisq-test (infer-ct contingency-table))]
     (m/sqrt (m// chi2 (m/+ n chi2))))))

(defn cramers-v
  "Calculates Cramer's V, a measure of association (effect size) between two
  nominal variables represented in a contingency table.

  Its value ranges from 0 to 1, where 0 indicates no association and 1 indicates
  a perfect association. It is related to the Pearson's Chi-squared statistic
  and is useful for tables of any size.

  The function can be called in two ways:

  1.  With two sequences `group1` and `group2`:
      The function will automatically construct a contingency table from
      the unique values in the sequences.
  2.  With a contingency table:
      The contingency table can be provided as:
      - A map where keys are `[row-index, column-index]` tuples and values are counts
        (e.g., `{[0 0] 10, [0 1] 5, [1 0] 3, [1 1] 12}`). This is the output format
        of [[contingency-table]] with two inputs.
      - A sequence of sequences representing the rows of the table
        (e.g., `[[10 5] [3 12]]`). This is equivalent to `rows->contingency-table`.

  Parameters:

  - `group1` (sequence): The first sequence of categorical data.
  - `group2` (sequence): The second sequence of categorical data. Must have the same length as `group1`.
  - `contingency-table` (map or sequence of sequences): A pre-computed contingency table.

  Returns the calculated Cramer's V coefficient as a double.

  See also [[chisq-test]], [[cramers-c]], [[cohens-w]], [[tschuprows-t]], [[contingency-table]]."
  (^double [group1 group2] (cramers-v (contingency-table group1 group2)))
  (^double [contingency-table]
   (let [{:keys [^double chi2 ^long k ^long r ^long n]} (chisq-test (infer-ct contingency-table))]
     (m/sqrt (m// (m// chi2 n)
                  (m/min (m/dec k) (m/dec r)))))))

(defn cramers-v-corrected
  "Calculates the **corrected Cramer's V**, a measure of association (effect size)
  between two nominal variables represented in a contingency table, with a correction
  to reduce bias, particularly for small sample sizes or tables with many cells
  having small expected counts.

  Like the uncorrected Cramer's V ([[cramers-v]]), its value ranges from 0 to 1,
  where 0 indicates no association and 1 indicates a perfect association. The
  correction tends to yield a value closer to the true population value in
  biased situations.

  The function can be called in two ways:

  1.  With two sequences `group1` and `group2`:
      The function will automatically construct a contingency table from
      the unique values in the sequences.
  2.  With a contingency table:
      The contingency table can be provided as:
      - A map where keys are `[row-index, column-index]` tuples and values are counts
        (e.g., `{[0 0] 10, [0 1] 5, [1 0] 3, [1 1] 12}`). This is the output format
        of [[contingency-table]] with two inputs.
      - A sequence of sequences representing the rows of the table
        (e.g., `[[10 5] [3 12]]`). This is equivalent to [[rows->contingency-table]].

  Parameters:

  - `group1` (sequence): The first sequence of categorical data.
  - `group2` (sequence): The second sequence of categorical data. Must have the same length as `group1`.
  - `contingency-table` (map or sequence of sequences): A pre-computed contingency table.

  Returns the calculated corrected Cramer's V coefficient as a double.

  See also [[chisq-test]], [[cramers-v]] (uncorrected), [[cramers-c]], [[cohens-w]],
  [[tschuprows-t]], [[contingency-table]]."
  (^double [group1 group2] (cramers-v-corrected (contingency-table group1 group2)))
  (^double [contingency-table]
   (let [{:keys [^double chi2 ^long k ^long r ^long n]} (chisq-test (infer-ct contingency-table))
         k1 (m/dec k)
         r1 (m/dec r)
         n1 (m/dec n)
         phi2_ (m/max 0.0 (m/- (m// chi2 n) (m// (m/* k1 r1) n1)))
         k_ (m/- k (m// (m/* k1 k1) n1))
         r_ (m/- r (m// (m/* r1 r1) n1))]
     (m/sqrt (m// phi2_
                  (m/min (m/dec k_) (m/dec r_)))))))

(defn cohens-w
  "Calculates Cohen's W effect size for the association between two nominal
  variables represented in a contingency table.

  Cohen's W is a measure of association derived from the Pearson's Chi-squared
  statistic. It quantifies the magnitude of the difference between the observed
  frequencies and the frequencies expected under the assumption of independence
  between the variables.

  Its value ranges from 0 upwards:
  - A value of 0 indicates no association between the variables.
  - Larger values indicate a stronger association.

  The function can be called in two ways:

  1.  With two sequences `group1` and `group2`:
      The function will automatically construct a contingency table from
      the unique values in the sequences.
  2.  With a contingency table:
      The contingency table can be provided as:
      - A map where keys are `[row-index, column-index]` tuples and values are counts
        (e.g., `{[0 0] 10, [0 1] 5, [1 0] 3, [1 1] 12}`). This is the output format
        of [[contingency-table]] with two inputs.
      - A sequence of sequences representing the rows of the table
        (e.g., `[[10 5] [3 12]]`). This is equivalent to [[rows->contingency-table]].

  Parameters:

  - `group1` (sequence): The first sequence of categorical data.
  - `group2` (sequence): The second sequence of categorical data. Must have the same length as `group1`.
  - `contingency-table` (map or sequence of sequences): A pre-computed contingency table.

  Returns the calculated Cohen's W coefficient as a double.

  See also [[chisq-test]], [[cramers-v]], [[cramers-c]], [[tschuprows-t]], [[contingency-table]]."
  (^double [group1 group2] (cohens-w (contingency-table group1 group2)))
  (^double [contingency-table]
   (let [{:keys [^double chi2 ^long n]} (chisq-test (infer-ct contingency-table))]
     (m/sqrt (m// chi2 n)))))

(defn tschuprows-t
  "Calculates Tschuprow's T, a measure of association between two nominal variables
  represented in a contingency table.

  Tschuprow's T is derived from the Pearson's Chi-squared statistic and measures
  the strength of the association. Its value ranges from 0 to 1.

  - A value of 0 indicates no association between the variables.
  - A value of 1 indicates perfect association, but only when the number of rows
    (`r`) equals the number of columns (`k`) in the contingency table. If `r != k`,
    Tschuprow's T cannot reach 1, making Cramer's V ([[cramers-v]]) often preferred
    as it can reach 1 for any table size.

  The function can be called in two ways:

  1.  With two sequences `group1` and `group2`:
      The function will automatically construct a contingency table from
      the unique values in the sequences.
  2.  With a contingency table:
      The contingency table can be provided as:
      - A map where keys are `[row-index, column-index]` tuples and values are counts
        (e.g., `{[0 0] 10, [0 1] 5, [1 0] 3, [1 1] 12}`). This is the output format
        of [[contingency-table]] with two inputs.
      - A sequence of sequences representing the rows of the table
        (e.g., `[[10 5] [3 12]]`). This is equivalent to [[rows->contingency-table]].

  Parameters:

  - `group1` (sequence): The first sequence of categorical data.
  - `group2` (sequence): The second sequence of categorical data. Must have the same length as `group1`.
  - `contingency-table` (map or sequence of sequences): A pre-computed contingency table.

  Returns the calculated Tschuprow's T coefficient as a double.

  See also [[chisq-test]], [[cramers-c]], [[cramers-v]], [[cohens-w]], [[contingency-table]]."
  (^double [group1 group2] (tschuprows-t (contingency-table group1 group2)))
  (^double [contingency-table]
   (let [{:keys [^double chi2 ^long k ^long r ^long n]} (chisq-test (infer-ct contingency-table))]
     (m/sqrt (m// (m// chi2 n)
                (m/sqrt (m/* (m/dec k) (m/dec r))))))))

(defn cohens-kappa
  "Calculates Cohen's Kappa coefficient (κ), a statistic that measures inter-rater
  agreement for categorical items, while correcting for chance agreement.

  It is often used to assess the consistency of agreement between two raters or
  methods. Its value typically ranges from -1 to +1:

  - `κ = 1`: Perfect agreement.
  - `κ = 0`: Agreement is no better than chance.
  - `κ < 0`: Agreement is worse than chance.

  The function can be called in two ways:

  1.  With two sequences `group1` and `group2`:
      The function will automatically construct a 2x2 contingency table from
      the unique values in the sequences (assuming they represent two binary
      variables). The mapping of values to table cells (e.g., what corresponds
      to TP, TN, FP, FN) depends on how `contingency-table` orders the unique values.
      For direct control over which cell is which, use the contingency table input.

  2.  With a contingency table:
      The contingency table can be provided as:
      - A map where keys are `[row-index, column-index]` tuples and values are counts
        (e.g., `{[0 0] TP, [0 1] FP, [1 0] FN, [1 1] TN}`). This is the output format
        of [[contingency-table]] with two inputs. The mapping of indices to TP/TN/FP/FN
        depends on the order of unique values in the original data if generated by
        [[contingency-table]], or the explicit structure if created manually or via
        [[rows->contingency-table]]. Standard convention maps `[0 0]` to TP, `[0 1]` to FP,
        `[1 0]` to FN, and `[1 1]` to TN for binary outcomes.
      - A sequence of sequences representing the rows of the table
        (e.g., `[[TP FP] [FN TN]]`). This is equivalent to [[rows->contingency-table]].

  Parameters:

  - `group1` (sequence): The first sequence of binary outcomes/categories.
  - `group2` (sequence): The second sequence of binary outcomes/categories.
    Must have the same length as `group1`.
  - `contingency-table` (map or sequence of sequences): A pre-computed 2x2 contingency table.
    The cell values should represent counts (e.g., TP, FN, FP, TN).

  Returns the calculated Cohen's Kappa coefficient as a double.

  See also [[weighted-kappa]] (for ordinal data with partial agreement), [[contingency-table]], [[contingency-2x2-measures]], [[binary-measures-all]]."
  (^double [group1 group2] (cohens-kappa (contingency-table group1 group2)))
  (^double [contingency-table]
   (let [{:keys [^long n rows cols diag]} (contingency-table->marginals (infer-ct contingency-table))
         n2 (m/* n n)
         pe (m// (sum (map (fn [r c]
                             (m/* (double (second r)) (double (second c)))) rows cols)) n2)
         p0 (m// (sum (map second diag)) n)]
     (m// (m/- p0 pe) (m/- 1.0 pe)))))

(defn- weighted-kappa-equal-spacing
  ^double [^long R ^long id1 ^long id2]
  (m/- 1.0 (m// (m/abs (m/- id1 id2)) R)))

(defn- weighted-kappa-fleiss-cohen
  ^double [^long R ^long id1 ^long id2]
  (m/- 1.0 (m// (m/sq (m/- id1 id2)) (m/sq R))))

(defn weighted-kappa
  "Calculates Cohen's weighted Kappa coefficient (κ) for a contingency table,
  allowing for partial agreement between categories, typically used for ordinal data.

  Weighted Kappa measures inter-rater agreement, similar to [[cohens-kappa]],
  but assigns different penalties to disagreements based on their magnitude.
  Disagreements between closely related categories are penalized less than
  disagreements between distantly related categories.

  The function can be called in two ways:

  1.  With two sequences `group1` and `group2`:
      The function will automatically construct a contingency table from
      the unique values in the sequences. These values are assumed to be ordinal
      and their position in the sorted unique value list determines their index.
      The mapping of values to table indices might need verification.

  2.  With a contingency table:
      The contingency table can be provided as:
      - A map where keys are `[row-index, column-index]` tuples and values are counts
        (e.g., `{[0 0] 10, [0 1] 5, [1 0] 3, [1 1] 12}`). This is the output format
        of [[contingency-table]] with two inputs. Indices are assumed to represent
        the ordered categories.
      - A sequence of sequences representing the rows of the table
        (e.g., `[[10 5] [3 12]]`). This is equivalent to [[rows->contingency-table]].
      The row and column indices are assumed to correspond to the ordered categories.

  Parameters:

  - `group1` (sequence): The first sequence of ordinal outcomes/categories.
  - `group2` (sequence): The second sequence of ordinal outcomes/categories.
    Must have the same length as `group1`.
  - `contingency-table` (map or sequence of sequences): A pre-computed contingency table
    where row and column indices correspond to ordered categories.
  - `weights` (keyword, function, or map, optional): Specifies the weighting scheme
    to quantify the difference between categories. Defaults to `:equal-spacing`.
    - `:equal-spacing` (default, linear weights): Penalizes disagreements linearly
      with the distance between categories. Weight is `1 - |i-j|/R`, where `i` is
      row index, `j` is column index, and `R` is the maximum dimension of the table (max(max_row_index, max_col_index)).
    - `:fleiss-cohen` (quadratic weights): Penalizes disagreements quadratically
      with the distance. Weight is `1 - (|i-j|/R)^2`.
    - (function `(fn [R id1 id2])`): A custom function that takes the maximum
      dimension `R`, row index `id1`, and column index `id2` and returns the weight
      (typically between 0 and 1, where 1 is perfect agreement).
    - (map `{[id1 id2] weight}`): A custom map providing weights for specific
      `[row-index, column-index]` pairs. Missing pairs default to a weight of 0.0.

  Returns the calculated weighted Cohen's Kappa coefficient as a double.

  Interpretation:

  - `κ_w = 1`: Perfect agreement.
  - `κ_w = 0`: Agreement is no better than chance.
  - `κ_w < 0`: Agreement is worse than chance.

  See also [[cohens-kappa]] (unweighted Kappa)."
  (^double [contingency-table] (weighted-kappa contingency-table :equal-spacing))
  (^double [contingency-table weights]
   (let [ct (infer-ct contingency-table)
         R (long (reduce (fn [^long c [^long a ^long b]]
                           (m/max c a b)) 0 (keys ct)))
         {:keys [^long n rows cols]} (contingency-table->marginals ct)
         wfn (cond
               (= weights :equal-spacing) weighted-kappa-equal-spacing
               (= weights :fleiss-cohen) weighted-kappa-fleiss-cohen
               (fn? weights) weights
               :else (fn [_ id1 id2] (get weights [id1 id2] 0.0)))
         n2 (m/* n n)
         pe (m// (sum (for [[^long idr ^long r] rows
                            [^long idc ^long c] cols]
                        (m/* (double (wfn R idr idc)) r c))) n2)
         p0 (m// (sum (map (fn [[[^long idr ^long idc] ^long v]]
                             (m/* (double (wfn R idr idc)) v)) ct)) n)]
     (m// (m/- p0 pe) (m/- 1.0 pe)))))

(defn durbin-watson
  "Calculates the Durbin-Watson statistic (d) for a sequence of residuals.

  This statistic is used to test for the presence of serial correlation,
  especially first-order (lag-1) autocorrelation, in the residuals from a
  regression analysis. Autocorrelation violates the assumption of independent errors.

  Parameters:

  - `rs` (sequence of numbers): The sequence of residuals from a regression model.
    The sequence should represent observations ordered by time or sequence index.

  Returns the calculated Durbin-Watson statistic as a double. The value ranges from 0 to 4.

  Interpretation:

  - Values near 2 suggest no first-order autocorrelation.
  - Values less than 2 suggest positive autocorrelation (residuals tend to be followed by residuals of the same sign).
  - Values greater than 2 suggest negative autocorrelation (residuals tend to be followed by residuals of the opposite sign)."
  [rs]
  (let [es (map (fn [[^double x1 ^double x2]] (m/- x1 x2)) (partition 2 1 rs))]
    (m// (v/dot es es) (v/dot rs rs))))

;; binary classification statistics

(defn confusion-matrix
  "Creates a 2x2 confusion matrix for binary classification.

  A confusion matrix summarizes the results of a binary classification problem, showing
  the counts of True Positives (TP), False Positives (FP), False Negatives (FN),
  and True Negatives (TN).

  TP: Actual is True, Predicted is True
  FP: Actual is False, Predicted is True (Type I error)
  FN: Actual is True, Predicted is False (Type II error)
  TN: Actual is False, Predicted is False

  The function supports several input formats:

  1.  `(confusion-matrix tp fn fp tn)`: Direct input of the four counts.
      - `tp` (long): True Positive count.
      - `fn` (long): False Negative count.
      - `fp` (long): False Positive count.
      - `tn` (long): True Negative count.

  2.  `(confusion-matrix confusion-matrix-representation)`: Input as a structured representation.
      - `confusion-matrix-representation`: Can be:
        - A map with keys like `:tp`, `:fn`, `:fp`, `:tn` (e.g., `{:tp 10 :fn 2 :fp 5 :tn 80}`).
        - A sequence of sequences representing rows `[[TP FP] [FN TN]]` (e.g., `[[10 5] [2 80]]`).
        - A flat sequence `[TP FN FP TN]` (e.g., `[10 2 5 80]`).

  3.  `(confusion-matrix actual prediction)`: Input as two sequences of outcomes.
      - `actual` (sequence): Sequence of true outcomes.
      - `prediction` (sequence): Sequence of predicted outcomes. Must have the same length as `actual`.
      Values in `actual` and `prediction` are compared element-wise. By default,
      any non-`nil` or non-zero value is treated as `true`, and `nil` or `0.0` is
      treated as `false`.

  4.  `(confusion-matrix actual prediction encode-true)`: Input as two sequences with a specified encoding for `true`.
      - `actual`, `prediction`: Sequences as in the previous arity.
      - `encode-true`: Specifies how values in `actual` and `prediction` are converted to boolean `true` or `false`.
        - `nil` (default): Non-`nil`/non-zero is true.
        - Any sequence/set: Values found in this collection are true.
        - A map: Values are mapped according to the map; if a key is not found or maps to `false`, the value is false.
        - A predicate function: Returns `true` if the value satisfies the predicate.

  Returns a map with keys `:tp`, `:fn`, `:fp`, and `:tn` representing the counts.

  This function is commonly used to prepare input for binary classification
  metrics like those provided by [[binary-measures-all]] and [[binary-measures]]."
  ([tp fn fp tn] {:tp tp :fn fn :fp fp :tn tn})
  ([confusion-mat] (binary/infer-confusion-matrix confusion-mat))
  ([actual prediction] (confusion-matrix actual prediction nil))
  ([actual prediction encode-true]
   (let [truth (binary/binary-process-list actual encode-true)
         prediction (binary/binary-process-list prediction encode-true)]
     (merge {:tp 0 :fn 0 :fp 0 :tn 0}
            (frequencies (map binary/binary-confusion truth prediction))))))

(def ^{:deprecated "Use `confusion-matrix`"} ->confusion-matrix confusion-matrix)

(defn binary-measures-all
  "Calculates a comprehensive set of evaluation metrics for binary classification results.

  This function computes various statistics derived from a 2x2 confusion matrix,
  summarizing the performance of a binary classifier.

  The 2x2 confusion matrix is based on True Positives (TP), False Positives (FP),
  False Negatives (FN), and True Negatives (TN):

  |                | Predicted True | Predicted False |
  |:---------------|:---------------|:----------------|
  | **Actual True**  | TP             | FN              |
  | **Actual False** | FP             | TN              |

  The function supports several input formats:

  1.  `(binary-measures-all tp fn fp tn)`: Direct input of the four counts as arguments.
      - `tp` (long): True Positive count.
      - `fn` (long): False Negative count.
      - `fp` (long): False Positive count.
      - `tn` (long): True Negative count.

  2.  `(binary-measures-all confusion-matrix)`: Input as a structured representation of the confusion matrix.
      - `confusion-matrix`: Can be:
        - A map with keys like `:tp`, `:fn`, `:fp`, `:tn` (e.g., `{:tp 10 :fn 2 :fp 5 :tn 80}`).
        - A sequence of sequences representing rows `[[TP FP] [FN TN]]` (e.g., `[[10 5] [2 80]]`).
        - A flat sequence `[TP FN FP TN]` (e.g., `[10 2 5 80]`).

  3.  `(binary-measures-all actual prediction)`: Input as two sequences of outcomes.
      - `actual` (sequence): Sequence of true outcomes.
      - `prediction` (sequence): Sequence of predicted outcomes. Must have the same length as `actual`.
      Values in `actual` and `prediction` are converted to boolean `true`/`false`. By default,
      any non-`nil` or non-zero numeric value is treated as `true`, and `nil` or `0.0` is
      treated as `false`.

  4.  `(binary-measures-all actual prediction true-value)`: Input as two sequences with a specified encoding for `true`.
      - `actual`, `prediction`: Sequences as in the previous arity.
      - `true-value` (optional): Specifies how values in `actual` and `prediction` are converted to boolean `true` (success) or `false` (failure).
        - `nil` (default): Non-`nil`/non-zero (for numbers) is true.
        - Any sequence/set: Values found in this collection are true.
        - A map: Values are mapped according to the map; if a key is not found or maps to `false`, the value is false.
        - A predicate function: Returns `true` if the value satisfies the predicate.

  Returns a map containing a wide array of calculated metrics. This includes, but is not limited to:

  - Basic Counts: `:tp`, `:fn`, `:fp`, `:tn`
  - Totals: `:cp` (Actual Positives), `:cn` (Actual Negatives), `:pcp` (Predicted Positives), `:pcn` (Predicted Negatives), `:total` (Grand Total)
  - Rates (often ratios of counts):
    - `:tpr` (True Positive Rate, Recall, Sensitivity, Hit Rate)
    - `:fnr` (False Negative Rate, Miss Rate)
    - `:fpr` (False Positive Rate, Fall-out)
    - `:tnr` (True Negative Rate, Specificity, Selectivity)
    - `:ppv` (Positive Predictive Value, Precision)
    - `:fdr` (False Discovery Rate, `1 - ppv`)
    - `:npv` (Negative Predictive Value)
    - `:for` (False Omission Rate, `1 - npv`)
  - Ratios/Odds:
    - `:lr+` (Positive Likelihood Ratio)
    - `:lr-` (Negative Likelihood Ratio)
    - `:dor` (Diagnostic Odds Ratio)
  - Combined Scores:
    - `:accuracy`
    - `:ba` (Balanced Accuracy)
    - `:fm` (Fowlkes–Mallows index)
    - `:pt` (Prevalence Threshold)
    - `:ts` (Threat Score, Jaccard index)
    - `:f-measure` / `:f1-score` (F1 Score, special case of F-beta score)
    - `:f-beta` (Function to calculate F-beta for any beta)
    - `:mcc` / `:phi` (Matthews Correlation Coefficient, Phi coefficient)
    - `:bm` (Bookmaker Informedness)
    - `:kappa` (Cohen's Kappa, for 2x2 table)
    - `:mk` (Markedness)

  Metrics are generally calculated using standard formulas based on the TP, FN, FP, TN counts.
  For more details on specific metrics, refer to standard classification literature or
  the Wikipedia page on [Precision and recall](https://en.wikipedia.org/wiki/Precision_and_recall),
  which covers many of these concepts.

  See also [[confusion-matrix]], [[binary-measures]] (for a selected subset of metrics),
  [[mcc]], [[contingency-2x2-measures-all]] (for a broader set of 2x2 table measures).
  "
  ([tp fn fp tn] (binary/binary-measures-all-calc (binary-measures-all {:tp tp :fn fn :fp fp :tn tn})))
  ([confusion-matrix] (binary/binary-measures-all-calc (binary/infer-confusion-matrix confusion-matrix)))
  ([actual prediction] (binary-measures-all actual prediction nil))
  ([actual prediction true-value]
   (let [truth (binary/binary-process-list actual true-value)
         prediction (binary/binary-process-list prediction true-value)]
     (binary/binary-measures-all-calc (frequencies (map binary/binary-confusion truth prediction))))))

(defn- cm-select-keys
  [cm]
  (select-keys cm [:tp :tn :fp :fn :accuracy :fdr :f-measure :fall-out :precision :recall :sensitivity :specificity :prevalence]))

(defn binary-measures
  "Calculates a selected subset of common evaluation metrics for binary classification results.

  This function is a convenience wrapper around [[binary-measures-all]], providing
  a map containing the most frequently used metrics derived from a 2x2 confusion matrix.

  The 2x2 confusion matrix is based on True Positives (TP), False Positives (FP),
  False Negatives (FN), and True Negatives (TN):

  |                | Predicted True | Predicted False |
  |:---------------|:---------------|:----------------|
  | **Actual True**  | TP             | FN              |
  | **Actual False** | FP             | TN              |

  The function accepts the same input formats as [[binary-measures-all]]:

  1.  `(binary-measures tp fn fp tn)`: Direct input of the four counts.
  2.  `(binary-measures confusion-matrix)`: Input as a structured representation
      (map with keys like `:tp`, `:fn`, `:fp`, `:tn`; sequence of sequences
      `[ [TP FP] [FN TN] ]`; or flat sequence `[TP FN FP TN]`).
  3.  `(binary-measures actual prediction)`: Input as two sequences of outcomes.
  4.  `(binary-measures actual prediction true-value)`: Input as two sequences with
      a specified encoding for `true` (success).

  Parameters:

  - `tp, fn, fp, tn` (long): Counts from the confusion matrix.
  - `confusion-matrix` (map or sequence): Representation of the confusion matrix.
  - `actual`, `prediction` (sequences): Sequences of true and predicted outcomes.
  - `true-value` (optional): Specifies how outcomes are converted to boolean `true`/`false`.

  Returns a map containing the following selected metrics:

  - `:tp` (True Positives)
  - `:tn` (True Negatives)
  - `:fp` (False Positives)
  - `:fn` (False Negatives)
  - `:accuracy`
  - `:fdr` (False Discovery Rate, 1 - Precision)
  - `:f-measure` (F1 Score, harmonic mean of Precision and Recall)
  - `:fall-out` (False Positive Rate)
  - `:precision` (Positive Predictive Value)
  - `:recall` (True Positive Rate / Sensitivity)
  - `:sensitivity` (Alias for Recall/TPR)
  - `:specificity` (True Negative Rate)
  - `:prevalence` (Proportion of positive cases)

  See also [[confusion-matrix]], [[binary-measures-all]], [[mcc]], [[contingency-2x2-measures-all]]."
  ([tp fn fp tn] (cm-select-keys (binary-measures-all tp fn fp tn)))
  ([confusion-matrix] (cm-select-keys (binary-measures-all confusion-matrix)))
  ([actual prediction] (binary-measures actual prediction nil))
  ([actual prediction true-value] (cm-select-keys (binary-measures-all actual prediction true-value))))

;; https://www.evidentlyai.com/classification-metrics/multi-class-metrics

(defn- ->metric-fn
  "Converts metric keyword to a function, processing `:f-beta` or `:f-inv-beta` case."
  [f ^double beta]
  (cond
    (= f :f-beta) (binary/->f-beta beta)
    (= f :f-inv-beta) (binary/->f-inv-beta beta)
    (keyword? f) (binary/measures f)
    :else f))

(defn multiclass-measure
  "Calculates an average of selected metric for multiclass binary classification results using one vs the rest strategy. Possible variants are `macro`, `weighted` and `micro`. 

  Possible metrics are the same as in [[binary-measures-all]].  

  Options:

  - `:average` - average function, `nil` or `:micro` (default: `mean`)
  - `:metric` - measure to calculate (default: `:f1-score`).
  - `:weighted?` - `weighted` average variant (default: `false`)
  - `:beta` - beta for `:f-beta` (default: `0.5`)"
  ([actual prediction] (multiclass-measure actual prediction nil))
  ([actual prediction {:keys [average metric weighted? ^double beta]
                       :or {average mean metric :f1-score weighted? false beta 1.0}}]
   (let [metric-fn (->metric-fn metric beta)
         all-classes (distinct (concat actual prediction))
         cm (partial confusion-matrix actual prediction)]
     (when-not metric-fn (throw (ex-info "Unknown metric!" {:metric metric})))
     (if (= average :micro)
       (->> (map cm all-classes)
            (apply merge-with m/+)
            (metric-fn))
       (let [measures (->> all-classes
                           (map (comp metric-fn cm))
                           (map (fn [^double v] (if (m/nan? v) 0.0 v))))]
         (cond
           (not average) (zipmap all-classes measures)
           weighted? (let [weights (map (frequencies actual) all-classes)]
                       (average measures weights))
           :else (average measures)))))))

(defn binary-measures-thr
  "Calculate binary metrics at given thresholds (a performance).
  
  Each measure is a sequence with value for given threshold. Ties are interpolated lineary.

  A `true-value` determines which label(s) is true.  

  Curves:

  * ROC (Receiver Operating Characteristic) - `:fpr` and `:tpr`
  * PR or PrecRec (Precision-Recall) - `:recall` and `:precision`
  * DET (Detection Error Tradeoff) - `:fpr` and `:fnr`

  Returned measures contain:

  * `:thr` - thresholds  
  * `:p`, `:n`, `:total`
  * `:tp`, `:fp`, `:fn`, `:tn`
  * `:tpr`/`:recall`/`:sensitivity`, `:fnr`/`:miss-rate`, `:fpr`/`:fallout`, `:tnr`/`:specificity`
  * `:ppv`/`:precision`, `:fdr`, `:for`, `:npv`
  * `:mcc`
  * `:f1-score`
  * `:kappa`
  * `:fm`
  * `:ts`/`:jaccard`"
  ([labels scores] (binary-measures-thr labels scores nil))
  ([labels scores true-value] (binary/binary-measures-thr labels scores true-value)))

(defn- auc-axis
  [performance axis]
  (cond
    (keyword? axis) (performance axis)
    (fn? axis) (apply map axis ((juxt :tp :fp :fn :tn) performance))
    :else axis))

(defn auc
  "Calculates the Area Under the Curve (AUC) using trapezoidal integration.

  The curve is defined by an x-axis and a y-axis, each resolved from a performance map,
  a plain sequence of values, or a metric function. Common use cases are:

  - ROC curve (default): `x-axis` = `:fpr`, `y-axis` = `:tpr`
  - Precision-Recall curve: `x-axis` = `:recall`, `y-axis` = `:precision`
  - DET curve: `x-axis` = `:fpr`, `y-axis` = `:fnr`

  Parameters:

  - `performance` - a performance map returned by [[binary-measures-thr]]; when called with
    one argument, defaults to the ROC curve (`:fpr` vs `:tpr`).
  - `x-axis`, `y-axis` - each can be:
    - a keyword looked up in `performance` (e.g. `:fpr`, `:tpr`, `:recall`, `:precision`, `:fnr`)
    - a sequence of numeric values used directly as axis coordinates
    - a function of `tp`, `fp`, `fn`, `tn` (called with the four columns from `performance`) that returns a sequence of values
  - `[x-axis y-axis]` - when called with two arguments and no `performance` map, treats both as
    plain sequences of coordinates and integrates directly.

  Returns the AUC as a double. The value is not constrained; it depends on the direction and
  range of the input axes. For a standard ROC curve the result lies in `[0.0, 1.0]`.
  Returns `0.0` for an empty or single-point input.

  See also [[auc-roc]], [[binary-measures-thr]], [[multiclass-auc]]."
  ([performance] (auc performance :fpr :tpr))
  ([performance x-axis y-axis]
   (auc (auc-axis performance x-axis)
        (auc-axis performance y-axis)))
  ([x-axis y-axis]
   (->> (map vector x-axis y-axis)
        (partition 2 1)
        (reduce (fn [^double curr [[^double x1 ^double y1] [^double x2 ^double y2]]]
                  (m/+ curr (m/* 0.5 (m/- x2 x1)
                                 (m/+ y1 y2)))) 0.0))))

(defn auc-roc
  "Calculates the ROC AUC (Area Under the Receiver Operating Characteristic Curve) using the U-statistic (Wilcoxon-Mann-Whitney statistic).

  This method computes AUC directly from labels and continuous scores without constructing the ROC curve explicitly.
  It is equivalent to the probability that a randomly chosen positive instance is ranked higher than a randomly chosen negative instance.
  The result is constrained to `[0.0, 1.0]`.

  Parameters:

  - `labels` - a sequence of class labels (binary).
  - `scores` - a sequence of continuous scores or predicted probabilities, one per label.
  - `true-value` - (optional) the label value to treat as positive; when `nil`, the positive class is inferred via [[binary-process-list]].

  Returns the AUC score as a double in the range `[0.0, 1.0]`. A value of `1.0` means perfect ranking,
  `0.5` corresponds to a random classifier, and values below `0.5` indicate worse-than-random ranking.
  Returns `0.0` when there are no positive or no negative examples.

  See also [[auc]], [[binary-measures-thr]], [[multiclass-auc]]."
  ([labels scores] (auc-roc labels scores nil))
  ([labels scores true-value]
   (let [labels (binary/binary-process-list labels true-value)
         {^long t true ^long f false :or {t 0 f 0}} (frequencies labels)
         rank (m/rank1 scores)
         sum (v/sum (map (fn [l r] (if l r 0.0)) labels rank))]
     (m/constrain (m// (m/- sum (m/* t (m/inc t) 0.5))
                       (m/* t f)) 0.0 1.0))))

(defn multiclass-auc
  "Calculates AUC for a multiclass problem using a one-vs-rest strategy.

  Each class is treated as the positive class in turn, and its AUC is computed against all remaining classes.
  The final result is aggregated across classes according to the chosen averaging strategy.

  Parameters:

  - `classes` - a sequence of true class labels (any comparable values).
  - `scores` - either a sequence of score vectors (one per instance, shared across all classes)
    or a map from class label to a sequence of per-instance scores for that class.
  - `opts` - (optional) a map of options:
    - `:metric` - the curve to integrate over (default: `:roc`); supported values:
      - `:roc` - ROC curve (`:fpr` vs `:tpr`)
      - `:pr` - Precision-Recall curve (`:recall` vs `:precision`)
      - `:det` - Detection Error Tradeoff curve (`:fpr` vs `:fnr`)
      - a two-element vector `[x-axis y-axis]` of keywords or metric functions accepted by [[auc]]
    - `:average` - how to aggregate per-class AUC values (default: `mean`); supported values:
      - a function (e.g. `mean`) - macro average, applies the function to per-class AUC scores
      - `:micro` - micro average, concatenates all classes before computing a single AUC
      - `nil` - returns a map of `{class auc}` with no aggregation
    - `:weighted?` - when `true`, passes per-class instance counts as weights to `:average` function (default: `false`); ignored when `:average` is `:micro` or `nil`

  Returns a single aggregated AUC double when `:average` is a function or `:micro`,
  or a map of `{class auc}` when `:average` is `nil`. `NaN` per-class AUC values are replaced with `0.0`.

  See also [[auc]], [[auc-roc]], [[binary-measures-thr]]."
  ([classes scores] (multiclass-auc classes scores nil))
  ([classes scores {:keys [average metric weighted?]
                    :or {average mean metric :roc weighted? false}}]
   (let [[x-axis y-axis] (case metric
                           :roc [:fpr :tpr]
                           :pr [:recall :precision]
                           :det [:fpr :fnr]
                           metric)
         all-classes (distinct classes)
         scores-map (if (map? scores) scores (zipmap all-classes (repeat scores)))]
     (if (= average :micro)
       (let [combined-classes (mapcat (fn [cl] (binary/binary-process-list classes #{cl})) all-classes)
             combined-scores (mapcat scores-map all-classes)]
         (auc (binary-measures-thr combined-classes combined-scores) x-axis y-axis))
       (let [measures (->> all-classes
                           (map (fn [cl] (-> (binary-measures-thr classes (scores-map cl) #{cl})
                                             (auc x-axis y-axis))))
                           (map (fn [^double v] (if (m/nan? v) 0.0 v))))]
         (cond
           (not average) (zipmap all-classes measures)
           weighted? (let [weights (map (frequencies classes) all-classes)]
                       (average measures weights))
           :else (average measures)))))))

;; contingency 2x2

(defn- contingency-2x2-measures-calc
  [^long a ^long b ^long c ^long d]
  (let [fields [:a :b :c :d]
        a (or a 0) b (or b 0) c (or c 0) d (or d 0)
        table [a b c d]
        r1 (m/long-add a b) r2 (m/long-add c d)
        c1 (m/long-add a c) c2 (m/long-add b d)
        n (m/long-add a b c d)
        dr1 (double r1) dr2 (double r2)
        dc1 (double c1) dc2 (double c2)
        dn (double n)
        expected [(m// (m/* r1 c1) dn) (m// (m/* r1 c2) dn) (m// (m/* r2 c1) dn) (m// (m/* r2 c2) dn)]
        proportions (v/div table n)
        chi2distr (r/distribution :chi-squared {:degrees-of-freedom 1})
        chi2 (let [diff (v/sub table expected)]
               (v/sum (v/ediv (v/emult diff diff) expected)))
        yates (let [diff (v/shift (v/abs (v/sub table expected)) -0.5)]
                (v/sum (v/ediv (v/emult diff diff) expected)))
        cmh (let [nmt (m// (m/* r1 c1) dn)]
              (m// (m/sq (m/- a nmt))
                   (m/* nmt (m// (m/* r2 c2) (m/* dn (m/dec dn))))))
        phi (m// (m/- (m/* a d) (m/* b c))
                 (m/sqrt (m/* r1 r2 c1 c2)))
        kappa (m// (m/* 2.0 (m/- (m/* a d) (m/* b c)))
                   (m/+ (m/* c1 r2) (m/* r1 c2)))
        G (m// (m/- (m/+ a d) (m/+ b c)) dn)
        OR (m// (m/* a d) (double (m/long-mult b c)))
        RD (m/- (m// a dr1) (m// c dr2))
        EER (m// a dr1) CER (m// c dr2) ARR (m/- CER EER)
        NNT (if (m/zero? ARR) ##Inf (m// ARR))
        RR (m// (m// a dr1) (m// c dr2))
        ad2 (m/+ (m/sq a) (m/sq d))
        hbc (m/* 0.5 (m/+ b c))
        pcc (m/sqrt (m// chi2 (m/+ chi2 dn)))]
    {:n n
     :table (zipmap fields table)
     :expected (zipmap fields expected)
     :marginals {:row1 (m/+ a b) :row2 (m/+ c d)
                 :col1 (m/+ a c) :col2 (m/+ b d)
                 :total n}
     :proportions {:table (zipmap fields proportions)
                   :rows (zipmap fields [(m// a dr1) (m// b dr1) (m// c dr2) (m// d dr2)])
                   :cols (zipmap fields [(m// a dc1) (m// b dc2) (m// c dc1) (m// d dc2)])
                   :marginals (merge (zipmap [:row1 :row2] (v/div [r1 r2] n))
                                     (zipmap [:col1 :col2] (v/div [c1 c2] n)))}
     :p-values {:chi2 (r/ccdf chi2distr chi2)
                :yates (r/ccdf chi2distr yates)
                :cochran-mantel-haenszel (r/ccdf chi2distr cmh)}
     :OR OR :lOR (m/log OR) :RR RR
     :risk {:RR RR :RRR (m/- 1.0 RR)
            :RD RD :ES r1 :CS r2
            :EER EER :CER CER :ARR ARR :NNT NNT
            :ARI (m/- ARR) :NNH (m/- NNT) :RRI (m/dec RR)
            :AFe (m// (m/dec RR) RR) :PFu (m/- 1.0 RR)}
     :SE (m/sqrt (v/sum (v/reciprocal table)))
     :measures {:chi2 chi2
                :yates yates
                :cochran-mantel-haenszel cmh
                :cohens-kappa kappa
                :yules-q (m// (m/dec OR) (m/inc OR))
                :holley-guilfords-g G
                :huberts-gamma (m/sq G)
                :youdens-j (m// (m/- (m/* a d) (m/* b c))
                                (m/* dr1 dr2))
                :yules-y (let [sad (m/sqrt (m/* a d))
                               sbc (m/sqrt (m/* b c))]
                           (m// (m/- sad sbc) (m/+ sad sbc)))
                :cramers-v (m/abs phi)
                :phi phi
                :scotts-pi (m// (m/- (m/* a d) (m/* hbc hbc))
                                (m/* (m/+ a hbc) (m/+ d hbc)))
                :cohens-h (m/* 2.0 (m/- (m/asin (m/sqrt (m// a dr1)))
                                        (m/asin (m/sqrt (m// c dr2)))))
                :PCC pcc
                :PCC-adjusted (m/* m/SQRT2 pcc)
                :TCC (m/cos (m// m/PI (m/inc (m/sqrt OR))))                
                :F1 (m// (m/+ a a)
                         (double (m/long-add a a b c)))
                :bangdiwalas-b (m// ad2 (m/+ (m/* dr1 dc1) (m/* dr2 dc2)))
                :mcnemars-chi2 (m// (m/sq (m/- b c)) (m/+ b c))
                :gwets-ac1 (let [bc (m/+ b c)
                                 hbc2 (m/* bc hbc)]
                             (m// (m/- ad2 hbc2)
                                  (m/+ ad2 hbc2 (m/* (m/+ a d) (m/+ b c)))))}}))

(defn contingency-2x2-measures-all
  "Calculates a comprehensive set of statistics and measures for a 2x2 contingency table.

  A 2x2 contingency table cross-tabulates two categorical variables, each with two levels.
  The table counts are typically represented as:

  +---+---+
  | a | b |
  +---+---+
  | c | d |
  +---+---+

  Where `a, b, c, d` are the counts in the respective cells. 

  This function calculates numerous measures, including:

  *   Chi-squared statistics (Pearson, Yates' corrected, CMH) and their p-values.
  *   Measures of association (Phi, Yule's Q, Holley-Guilford's G, Hubert's Gamma, Yule's Y, Cramer's V, Scott's Pi, Cohen's H, Pearson/Tschuprow's CC).
  *   Measures of agreement (Cohen's Kappa).
  *   Risk and effect size measures (Odds Ratio (OR), Relative Risk (RR), Risk Difference (RD), NNT, etc.).
  *   Table marginals and proportions.

  The function can be called with the four counts directly or with a representation
  of the contingency table:

  1.  `(contingency-2x2-measures-all a b c d)`: Takes the four counts as arguments.
  2.  `(contingency-2x2-measures-all [a b c d])`: Takes a sequence of the four counts.
  3.  `(contingency-2x2-measures-all [[a b] [c d]])`: Takes a sequence of sequences representing the rows.
  4.  `(contingency-2x2-measures-all {:a a :b b :c c :d d})`: Takes a map of counts (accepts `:a/:b/:c/:d` keys).

  Parameters:

  - `a` (long): Count in the top-left cell.
  - `b` (long): Count in the top-right cell.
  - `c` (long): Count in the bottom-left cell.
  - `d` (long): Count in the bottom-right cell.
  - `map-or-seq` (map or sequence): A representation of the 2x2 table as described above.

  Returns a map containing a wide range of calculated statistics. Keys include:
  `:n`, `:table`, `:expected`, `:marginals`, `:proportions`, `:p-values` (map), `:OR`, `:lOR`, `:RR`, `:risk` (map), `:SE`, `:measures` (map).

  See also [[contingency-2x2-measures]] for a selected subset of these measures,
  [[mcc]] for the Matthews Correlation Coefficient (Phi), and [[binary-measures-all]]
  for metrics derived from a confusion matrix (often a 2x2 table in binary classification)."
  ([^long a ^long b ^long c ^long d] (contingency-2x2-measures-calc a b c d))
  ([map-or-seq]
   (if (map? map-or-seq)
     (let [{:keys [a b c d] :or {a 0 b 0 c 0 d 0}} map-or-seq]
       (contingency-2x2-measures-all a b c d))
     (apply contingency-2x2-measures-calc map-or-seq)))
  ([[^long a ^long b] [^long c ^long d]] (contingency-2x2-measures-all a b c d)))

(defn contingency-2x2-measures
  "Calculates a subset of common statistics and measures for a 2x2 contingency table.

  This function provides a selection of the most frequently used measures from the
  more comprehensive [[contingency-2x2-measures-all]].

  The function accepts the same input formats as [[contingency-2x2-measures-all]]:

  1.  `(contingency-2x2-measures a b c d)`: Takes the four counts as arguments.
  2.  `(contingency-2x2-measures [a b c d])`: Takes a sequence of the four counts.
  3.  `(contingency-2x2-measures [[a b] [c d]])`: Takes a sequence of sequences representing the rows.
  4.  `(contingency-2x2-measures {:a a :b b :c c :d d})`: Takes a map of counts (accepts `:a/:b/:c/:d` keys).

  Parameters:

  - `a, b, c, d` (long): Counts in the 2x2 table cells.
  - `map-or-seq` (map or sequence): A representation of the 2x2 table.

  Returns a map containing a selection of measures:

  - `:OR`: Odds Ratio (Odds Ratio)
  - `:chi2`: Pearson's Chi-squared statistic
  - `:yates`: Yates' continuity corrected Chi-squared statistic
  - `:cochran-mantel-haenszel`: Cochran-Mantel-Haenszel statistic
  - `:cohens-kappa`: Cohen's Kappa coefficient
  - `:yules-q`: Yule's Q measure of association
  - `:holley-guilfords-g`: Holley-Guilford's G measure
  - `:huberts-gamma`: Hubert's Gamma measure
  - `:yules-y`: Yule's Y measure of association
  - `:cramers-v`: Cramer's V measure of association
  - `:phi`: Phi coefficient (Matthews Correlation Coefficient)
  - `:scotts-pi`: Scott's Pi measure of agreement
  - `:cohens-h`: Cohen's H measure
  - `:PCC`: Pearson's Contingency Coefficient
  - `:PCC-adjusted`: Adjusted Pearson's Contingency Coefficient
  - `:TCC`: Tschuprow's Contingency Coefficient
  - `:F1`: F1 Score
  - `:bangdiwalas-b`: Bangdiwala's B statistic
  - `:mcnemars-chi2`: McNemar's Chi-squared test statistic
  - `:gwets-ac1`: Gwet's AC1 measure

  For a more comprehensive set of 2x2 measures and their detailed descriptions, see [[contingency-2x2-measures-all]]."
  [& args]
  (let [m (apply contingency-2x2-measures-all args)]
    (-> (:measures m)
        (assoc :OR (:OR m)))))

;; acf/pacf

(defn- cov-for-acf
  ^double [xs1 xs2]
  (reduce + 0.0 (map * xs1 xs2)))

;; http://feldman.faculty.pstat.ucsb.edu/174-03/lectures/l12
(defn acf
  "Calculates the Autocorrelation Function (ACF) for a given time series `data`.

  The ACF measures the linear dependence between a time series and its lagged values.
  It helps identify patterns (like seasonality or trend) and inform the selection of
  models for time series analysis (e.g., in ARIMA modeling).

  Parameters:

  * `data` (seq of numbers): The time series data.
  * `lags` (long or seq of longs, optional):
    * If a number, calculates ACF for lags from 0 up to this maximum lag.
    * If a sequence of numbers, calculates ACF for each lag specified in the sequence.
    * If omitted (1-arity call), calculates ACF for lags from 0 up to `(m/dec (count data))`.

  Returns a sequence of doubles: the autocorrelation coefficients for the specified lags.
  The value at lag 0 is always 1.0.

  See also [[acf-ci]] (Calculates ACF with confidence intervals), [[pacf]], [[pacf-ci]]."
  ([data] (acf data (m/dec (count data))))
  ([data lags]
   (let [vdata (vec (demean data))
         rcov0 (m// (cov-for-acf vdata vdata))]
     (map (fn [^long lag]
            (if (m/zero? lag)
              1.0
              (let [v2 (subvec vdata lag)
                    v1 (subvec vdata 0 (count v2))]
                (m/* rcov0 (cov-for-acf v1 v2))))) (if (number? lags)
                                                     (range (m/inc (long lags)))
                                                     (seq lags))))))

;; http://feldman.faculty.pstat.ucsb.edu/174-03/lectures/l13
(defn pacf
  "Calculates the Partial Autocorrelation Function (PACF) for a given time series `data`.

  The PACF measures the linear dependence between a time series and its lagged values *after removing* the effects of the intermediate lags. It helps identify the direct relationship at each lag and is used to determine the order of autoregressive (AR) components in time series models (e.g., ARIMA).

  Parameters:

  * `data` (seq of numbers): The time series data.
  * `lags` (long, optional): The maximum lag for which to calculate the PACF. If omitted, calculates PACF for lags from 0 up to `(dec (count data))`.

  Returns a sequence of doubles representing the partial autocorrelation coefficients for the specified lags. The value at lag 0 is always 0.0.

  See also [[acf]], [[acf-ci]], [[pacf-ci]]."
  ([data] (pacf data (m/long-dec (count data))))
  ([data ^long lags]
   (let [acfs (vec (acf data lags))
         phis (reductions (fn [curr ^long id]
                            (let [phi (m// (m/- (double (acfs id))
                                                (sum
                                                 (map-indexed (fn [^long idx ^double c]
                                                                (m/* c (double (acfs (m/dec (m/- id idx)))))) curr)))
                                           (m/- 1.0
                                                (sum (map-indexed (fn [^long id ^double c]
                                                                    (m/* c (double (acfs (m/inc id))))) curr))))]

                              (conj (mapv (fn [^double p1 ^double p2]
                                            (m/- p1 (m/* phi p2))) curr (reverse curr)) phi))) [(acfs 1)] (range 2 (m/inc lags)))]
     (conj (map last phis) 0.0))))

(defn- p-acf-ci-value
  ^double [data ^double alpha]
  (m/* (m// (m/sqrt (count data)))
       (double (r/icdf r/default-normal (m/* 0.5 (m/inc (m/- 1.0 alpha)))))))

(defn pacf-ci
  "Calculates the Partial Autocorrelation Function (PACF) for a time series and provides approximate confidence intervals.

  This function computes the PACF of the input time series `data` for specified lags
  (see [[pacf]]) and includes approximate confidence intervals around the PACF
  estimates. These intervals help determine whether the partial autocorrelation at
  a specific lag is statistically significant (i.e., likely non-zero in the population).

  Parameters:

  * `data` (seq of numbers): The time series data.
  * `lags` (long, optional): The maximum lag for which to calculate the PACF and CI.
    If omitted, calculates for lags up to `(dec (count data))`.
  * `alpha` (double, optional): The significance level for the confidence intervals.
    Defaults to `0.05` (for a 95% CI).

  Returns a map containing:

  * `:ci` (double): The value of the approximate standard confidence interval bound
    for lags > 0. If the absolute value of a PACF
    coefficient at lag `k > 0` exceeds this value, it is considered statistically significant.
  * `:pacf` (seq of doubles): The sequence of partial autocorrelation coefficients
    at lags from 0 up to `lags` (calculated using [[pacf]]).

  See also [[pacf]], [[acf]], [[acf-ci]]."
  ([data] (pacf-ci data (m/dec (count data))))
  ([data lags] (pacf-ci data lags 0.05))
  ([data ^long lags ^double alpha]
   (let [pacf-data (pacf data lags)
         ci (p-acf-ci-value data alpha)]
     {:ci ci
      :pacf pacf-data})))

(defn acf-ci
  "Calculates the Autocorrelation Function (ACF) for a time series and provides approximate confidence intervals.

  This function computes the ACF of the input time series `data` for specified lags
  (see [[acf]]) and includes approximate confidence intervals around the ACF
  estimates. These intervals help determine whether the autocorrelation at a
  specific lag is statistically significant (i.e., likely non-zero in the population).

  Parameters:

  * `data` (seq of numbers): The time series data.
  * `lags` (long or seq of longs, optional):
    * If a number, calculates ACF for lags from 0 up to this maximum lag.
    * If a sequence of numbers, calculates ACF for each lag specified in the sequence.
    * If omitted (1-arity call), calculates ACF for lags from 0 up to `(dec (count data))`.
  * `alpha` (double, optional): The significance level for the confidence intervals.
    Defaults to `0.05` (for a 95% CI).

  Returns a map containing:

  * `:ci` (double): The value of the approximate standard confidence interval bound
    for lags > 0. If the absolute value of an ACF
    coefficient at lag `k > 0` exceeds this value, it is considered statistically significant.
  * `:acf` (seq of doubles): The sequence of autocorrelation coefficients
    at lags from 0 up to `lags` (or specified lags if `lags` is a sequence), calculated
    using [[acf]].
  * `:cis` (seq of doubles): Cumulative confidence intervals for ACF. These are based on the
    variance of the sum of squared sample autocorrelations up to each lag.

  See also [[acf]], [[pacf]], [[pacf-ci]]."
  ([data] (acf-ci data (m/dec (count data))))
  ([data lags] (acf-ci data lags 0.05))
  ([data ^long lags ^double alpha]
   (let [acf-data (acf data lags)
         ci (p-acf-ci-value data alpha)]
     {:ci ci
      :acf acf-data
      :cis (map (fn [^double r]
                  (m/* ci (m/sqrt (m/dec (m/+ r r))))) (reductions (fn [^double acc ^double s]
                                                               (m/+ acc (m/* s s))) acf-data))})))

;;

(defn- estimate-acceleration
  "Estimates acceleration for BCA bootstrap confidence interval computation"
  ^double [avs]
  (m// (skewness avs :skew) -6.0))

(defn- cdf-accelerated-quantile
  ^double [^double z0 ^double z ^double a]
  (let [num (m/+ z0 z)
        denom (m/- 1.0 (m/* a num))]
    (->> (m/+ z0 (m// num denom))
         (r/cdf r/default-normal))))

(defn- empirical-cdf
  ^double [vs ^double value]
  (m// (double (count (filter (fn [^double v] (m/< v value)) vs))) (count vs)))

(defn- percentile-bca-common
  [avs p1 p2 m accel estimation-strategy]
  (let [z0 (double (r/icdf r/default-normal (empirical-cdf avs m)))
        z1 (double (r/icdf r/default-normal (m// (double p1) 100.0)))
        z2 (double (r/icdf r/default-normal (m// (double p2) 100.0)))
        q1 (cdf-accelerated-quantile z0 z1 accel)
        q2 (cdf-accelerated-quantile z0 z2 accel)]
    [(quantile avs q1 estimation-strategy)
     (quantile avs q2 estimation-strategy)
     m]))

(defn percentile-bca-extent
  "Return bias corrected percentile range and mean for bootstrap samples. Also accounts for variance
   variations throught the accelaration parameter.
  See https://projecteuclid.org/euclid.ss/1032280214

  `p` - calculates extent of bias corrected `p` and `100-p` (default: `p=2.5`)

  Set `estimation-strategy` to `:r7` to get the same result as in R `coxed::bca`."
  ([vs] (percentile-bca-extent vs 2.5))
  ([vs ^double p] (percentile-bca-extent vs p (m/- 100.0 p)))
  ([vs p1 p2] (percentile-bca-extent vs p1 p2 :legacy))
  ([vs p1 p2 estimation-strategy]
   (let [avs (m/seq->double-array vs)
         accel (estimate-acceleration avs)]
     (percentile-bca-common avs p1 p2 (mean avs) accel estimation-strategy)))
  ([vs p1 p2 accel estimation-strategy]
   (let [avs (m/seq->double-array vs)]
     (percentile-bca-common avs p1 p2 (mean avs) accel estimation-strategy))))

(defn percentile-bc-extent
  "Return bias corrected percentile range and mean for bootstrap samples.
  See https://projecteuclid.org/euclid.ss/1032280214

  `p` - calculates extent of bias corrected `p` and `100-p` (default: `p=2.5`)

  Set `estimation-strategy` to `:r7` to get the same result as in R `coxed::bca`."
  ([vs] (percentile-bc-extent vs 2.5))
  ([vs ^double p] (percentile-bc-extent vs p (m/- 100.0 p)))
  ([vs p1 p2] (percentile-bc-extent vs p1 p2 :legacy))
  ([vs p1 p2 estimation-strategy]
   (percentile-bca-extent vs p1 p2 0.0 estimation-strategy)))

;;

(def binomial-ci-methods (sort [:asymptotic :agresti-coull :clopper-pearson :wilson :prop.test
                                :cloglog :logit :probit :arcsine]))

(defn binomial-ci
  "Calculates a confidence interval for a binomial proportion.

  Given the number of observed `successes` in a fixed number of `trials`, this function
  estimates a confidence interval for the true underlying probability of success (`p`).

  Different statistical methods are available for calculating the interval, as the
  accuracy and behavior of the interval can vary, especially for small sample sizes
  or probabilities close to 0 or 1.

  Parameters:

  - `number-of-successes` (long): The count of successful outcomes.
  - `number-of-trials` (long): The total number of independent trials.
  - `method` (keyword, optional): The method used to calculate the confidence interval.
    Defaults to `:asymptotic`.
  - `alpha` (double, optional): The significance level (alpha) for the interval.
    The confidence level is `1 - alpha`. Defaults to `0.05` (yielding a 95% CI).

  Available `method` values:

  - `:asymptotic`: Normal approximation interval (Wald interval), based on the Central Limit Theorem. Simple but can be inaccurate for small samples or probabilities near 0 or 1.
  - `:agresti-coull`: An adjustment to the asymptotic interval, adding 'pseudo-counts' to improve performance for small samples.
  - `:clopper-pearson`: An exact method based on inverting binomial tests. Provides guaranteed coverage but can be overly conservative (wider than necessary).
  - `:wilson`: Score interval, derived from the score test. Generally recommended as a good balance of accuracy and coverage for various sample sizes.
  - `:prop.test`: Interval typically used with `prop.test` in R, applies a continuity correction.
  - `:cloglog`: Confidence interval based on the complementary log-log transformation.
  - `:logit`: Confidence interval based on the logit transformation.
  - `:probit`: Confidence interval based on the probit transformation (inverse of standard normal CDF).
  - `:arcsine`: Confidence interval based on the arcsine transformation.
  - `:all`: Applies all available methods and returns a map where keys are method keywords and values are their respective confidence intervals (as triplets).

  Returns:

  - A vector `[lower-bound, upper-bound, estimated-p]`.
    - `lower-bound` (double): The lower limit of the confidence interval.
    - `upper-bound` (double): The upper limit of the confidence interval.
    - `estimated-p` (double): The observed proportion of successes (`number-of-successes / number-of-trials`).

  If `method` is `:all`, returns a map of results from each method.

  See also [[binomial-test]] for performing a hypothesis test on a binomial proportion."
  ([^long number-of-successes ^long number-of-trials]
   (binomial-ci number-of-successes number-of-trials :asymptotic))
  ([^long number-of-successes ^long number-of-trials method]
   (binomial-ci number-of-successes number-of-trials method 0.05))
  ([^long number-of-successes ^long number-of-trials method ^double alpha]
   (let [p (m// (double number-of-successes) number-of-trials)
         alpha2 (m/* 0.5 alpha)
         z (double (r/icdf r/default-normal (m/- 1.0 alpha2)))
         z2 (m/* z z)
         x0? (m/zero? number-of-successes)
         xn? (m/== number-of-trials number-of-successes)]
     (case method
       :all (into {} (map #(vector %1 (binomial-ci number-of-successes number-of-trials % alpha))
                          binomial-ci-methods))
       :cloglog (let [logp (m/log p)
                      mu (m/log (m/- logp))
                      sd (m/* z (m/sqrt (-> (m/- 1.0 p) (m// number-of-trials) (m// p) (m// (m/* logp logp)))))
                      lcl (cond
                            x0? 0.0
                            xn? (m/pow alpha2 (m// 1.0 number-of-trials))
                            :else (m/exp (m/- (m/exp (m/+ mu sd)))))
                      ucl (cond
                            x0? (m/- 1.0 (m/pow alpha2 (m// 1.0 number-of-trials)))
                            xn? 1.0
                            :else (m/exp (m/- (m/exp (m/- mu sd)))))]
                  [lcl ucl p])
       :logit (let [logitp (m/- (m/log p) (m/log1p (m/- p)))
                    sd (m/* z (m/sqrt (-> (m// 1.0 number-of-trials) (m// p) (m// (m/- 1.0 p)))))
                    lcl (cond
                          x0? 0.0
                          xn? (m/pow alpha2 (m// 1.0 number-of-trials))
                          :else (let [lcl (m/exp (m/- logitp sd))]
                                  (m// lcl (m/inc lcl))))
                    ucl (cond
                          x0? (m/- 1.0 (m/pow alpha2 (m// 1.0 number-of-trials)))
                          xn? 1.0
                          :else (let [ucl (m/exp (m/+ logitp sd))]
                                  (m// ucl (m/inc ucl))))]
                [lcl ucl p])
       :probit (let [probitp (double (r/icdf r/default-normal p))
                     sd (m/* z (m/sqrt (-> (m/* p (m/- 1.0 p))
                                           (m// number-of-trials)
                                           (m// (m/sq (r/pdf r/default-normal probitp))))))
                     lcl (cond
                           x0? 0.0
                           xn? (m/pow alpha2 (m// 1.0 number-of-trials))
                           :else (r/cdf r/default-normal (m/- probitp sd)))
                     ucl (cond
                           x0? (m/- 1.0 (m/pow alpha2 (m// 1.0 number-of-trials)))
                           xn? 1.0
                           :else (r/cdf r/default-normal (m/+ probitp sd)))]
                 [lcl ucl p])
       :prop.test (let [yatesn (m// (m/min 0.5 (m/abs (m/- number-of-successes (m/* number-of-trials 0.5))))
                                    number-of-trials)
                        nn (m/* 2.0 number-of-trials)
                        z22n (m// z2 nn)
                        z22n+ (m/inc (m/* 2.0 z22n))
                        z22n2n (m// z22n nn)
                        pc (m/- p yatesn)
                        pl (if-not (m/pos? pc) 0.0
                                   (m// (m/- (m/+ pc z22n) (m/* z (m/sqrt (m/+ (m/* pc (m// (m/- 1.0 pc) number-of-trials))
                                                                               z22n2n))))
                                        z22n+))
                        pc (m/+ p yatesn)
                        pu (if (m/>= pc 1.0) 1.0
                               (m// (m/+ pc z22n (m/* z (m/sqrt (m/+ (m/* pc (m// (m/- 1.0 pc) number-of-trials))
                                                                     z22n2n))))
                                    z22n+))]
                    [pl pu p])
       :wilson (let [z2n (m// z2 number-of-trials)
                     p1 (m/+ p (m/* 0.5 z2n))
                     p2 (m/* z (m/sqrt (m// (m/+ (m/* p (m/- 1.0 p))
                                                 (m/* 0.25 z2n)) number-of-trials)))
                     p3 (m/inc z2n)]
                 [(m// (m/- p1 p2) p3) (m// (m/+ p1 p2) p3) p])
       :clopper-pearson (let [diff (m/- number-of-trials
                                        number-of-successes)
                              lclbeta (if x0? 1.0 (double (r/icdf
                                                           (r/distribution :beta
                                                                           {:alpha (m/inc diff)
                                                                            :beta number-of-successes})
                                                           (m/- 1.0 alpha2))))
                              uclbeta (if xn? 0.0 (double (r/icdf
                                                           (r/distribution :beta
                                                                           {:alpha diff
                                                                            :beta (m/inc number-of-successes)})
                                                           alpha2)))]
                          [(m/- 1.0 lclbeta) (m/- 1.0 uclbeta) p])
       :agresti-coull (let [x (m/+ number-of-successes (m/* 0.5 z2))
                            n (m/+ number-of-trials z2)
                            p' (m// x n)
                            zse (m/* z (m/sqrt (m/* p' (m// (m/- 1.0 p') n))))]
                        [(m/- p' zse) (m/+ p' zse) p])
       :arcsine (let [ap (m/asin (m/sqrt p))
                      zn (m// z (m/* 2.0 (m/sqrt number-of-trials)))]
                  [(m/sq (m/sin (m/max 0.0 (m/- ap zn))))
                   (m/sq (m/sin (m/min m/HALF_PI (m/+ ap zn)))) p])
       (let [zse (m/* z (m/sqrt (m/* p (m// (m/- 1.0 p) number-of-trials))))]
         [(m/- p zse) (m/+ p zse) p])))))

;; tests

;; t-test, reimplementation of R version

(defmacro ^:private sides-case
  [sides both right left]
  `(case ~sides
     (:two-sided :both) ~both
     (:one-sided-greater :right) ~right
     ~left))

(defn p-value
  "Calculates the p-value for a given test statistic based on a reference probability distribution.

  The p-value represents the probability of observing a test statistic as extreme as,
  or more extreme than, the provided `stat`, assuming the null hypothesis is true
  (where the null hypothesis implies `stat` follows the given `distribution`).

  Parameters:

  - `distribution` (distribution object, optional): The probability distribution object
    (from `fastmath.random`) that the test statistic follows under the null
    hypothesis. Defaults to the standard normal distribution (`fastmath.random/default-normal`)
    if omitted.
  - `stat` (double): The observed value of the test statistic.
  - `sides` (keyword, optional): Specifies the type of alternative hypothesis and
    how 'extremeness' is defined. Defaults to `:two-sided`.
    - `:two-sided` or `:both`: Alternative hypothesis is that the true parameter is
      different from the null value (tests for extremeness in either tail).
      Calculates `2 * min(CDF(stat), CCDF(stat))` (adjusted for discrete).
    - `:one-sided-greater` or `:right`: Alternative hypothesis is that the true
      parameter is greater than the null value (tests for extremeness in the right tail).
      Calculates `CCDF(stat)` (adjusted for discrete).
    - `:one-sided-less`, `:left`, or `:one-sided`: Alternative hypothesis is that the true
      parameter is less than the null value (tests for extremeness in the left tail).
      Calculates `CDF(stat)`.

  Note: For discrete distributions, a continuity correction (`stat - 1` for CCDF calculations)
  is applied when calculating right-tail or two-tail probabilities involving the
  upper tail. This ensures the probability mass *at* the statistic value is correctly
  accounted for.

  Returns the calculated p-value (a double between 0.0 and 1.0)."
  ([^double stat] (p-value r/default-normal stat))
  ([distribution ^double stat] (p-value distribution stat :two-sided))
  ([distribution ^double stat sides]
   (let [stat2 (if (r/continuous? distribution) stat (m/dec stat))]
     (sides-case sides
                 (m/min 1.0 (m/* 2.0 (m/min (r/cdf distribution stat)
                                            (r/ccdf distribution stat2))))
                 (r/ccdf distribution stat2)
                 (r/cdf distribution stat)))))

(defn skewness-test
  "Performs the D'Agostino test for normality based on sample skewness.

  This test assesses the null hypothesis that the data comes from a normally
  distributed population by checking if the sample skewness significantly deviates
  from the zero skewness expected under normality.

  The test works by:

  1. Calculating the sample skewness (type configurable via `:type`, default `:g1`).
  2. Standardizing the sample skewness relative to its expected value (0) and
     standard error under the null hypothesis.
  3. Applying a further transformation (inverse hyperbolic sine based) to this
     standardized score to yield a final test statistic `Z` that more closely
     follows a standard normal distribution under the null hypothesis.

  Parameters:

  - `xs` (seq of numbers): The sample data.
  - `skew` (double, optional): A pre-calculated skewness value. If omitted, it's calculated from `xs`.
  - `params` (map, optional): Options map:
    - `:sides` (keyword, default `:two-sided`): Specifies the alternative hypothesis.
      - `:two-sided` (default): The population skewness is different from 0.
      - `:one-sided-greater`: The population skewness is greater than 0 (right-skewed).
      - `:one-sided-less`: The population skewness is less than 0 (left-skewed).
    - `:type` (keyword, default `:g1`): The type of skewness to calculate if `skew` is not provided. Note that the internal normalization constants are derived based on `:g1`. See [[skewness]] for options.

  Returns a map containing:

  - `:Z`: The final test statistic, approximately standard normal under H0.
  - `:stat`: Alias for `:Z`.
  - `:p-value`: The p-value associated with `Z` and the specified `:sides`.
  - `:skewness`: The sample skewness value used in the test (either provided or calculated).

  See also [[kurtosis-test]], [[normality-test]], [[jarque-bera-test]]."
  ([xs] (skewness-test xs nil))
  ([xs params] (skewness-test xs nil params))
  ([xs skew {:keys [sides type]
             :or {sides :two-sided type :g1}}]
   (let [skew (double (or skew (skewness xs type)))
         n (count xs)
         y (m/* skew (m/sqrt (m// (m/* (m/inc n) (m/+ n 3))
                                  (m/* 6.0 (m/- n 2)))))
         beta2- (m/dec (m// (m/* 3.0 (m/+ (m/* n n) (m/* 27 n) -70) (m/+ n 1) (m/+ n 3))
                            (m/* (m/- n 2) (m/+ n 5) (m/+ n 7) (m/+ n 9))))
         w2 (m/dec (m/sqrt (m/* 2.0 beta2-)))
         delta (m// 1.0 (m/sqrt (m/* 0.5 (m/log w2))))
         alpha (m/sqrt (m// 2.0 (m/dec w2)))
         ya (double (if (m/zero? y) (m// 1.0 alpha) (m// y alpha)))
         Z (m/* delta (m/log (m/+ ya (m/sqrt (m/inc (m/* ya ya))))))]
     {:p-value (p-value r/default-normal Z sides)
      :Z Z
      :skewness skew})))

(defn kurtosis-test
  "Performs a test for normality based on sample kurtosis.

  This test assesses the null hypothesis that the data comes from a normally
  distributed population by checking if the sample kurtosis significantly deviates
  from the kurtosis expected under normality (approximately 3).

  The test works by:

  1. Calculating the sample kurtosis (type configurable via `:type`, default `:kurt`).
  2. Standardizing the difference between the sample kurtosis and the expected
     kurtosis under normality using the theoretical standard error.
  3. Applying a further transformation (e.g., Anscombe-Glynn/D'Agostino) to this standardized
     score to yield a final test statistic `Z` that more closely follows a
     standard normal distribution under the null hypothesis, especially for
     smaller sample sizes.

  Parameters:

  - `xs` (seq of numbers): The sample data.
  - `kurt` (double, optional): A pre-calculated kurtosis value. If omitted, it's calculated from `xs`.
  - `params` (map, optional): Options map:
    - `:sides` (keyword, default `:two-sided`): Specifies the alternative hypothesis.
      - `:two-sided` (default): The population kurtosis is different from normal.
      - `:one-sided-greater`: The population kurtosis is greater than normal (leptokurtic).
      - `:one-sided-less`: The population kurtosis is less than normal (platykurtic).
    - `:type` (keyword, default `:kurt`): The type of kurtosis to calculate if `kurt` is not provided. See [[kurtosis]] for options (e.g., `:kurt`, `:G2`, `:g2`).

  Returns a map containing:

  - `:Z`: The final test statistic, approximately standard normal under H0.
  - `:stat`: Alias for `:Z`.
  - `:p-value`: The p-value associated with `Z` and the specified `:sides`.
  - `:kurtosis`: The sample kurtosis value used in the test (either provided or calculated).

  See also [[skewness-test]], [[normality-test]], [[jarque-bera-test]], [[bonett-seier-test]]."
  ([xs] (kurtosis-test xs nil))
  ([xs params] (kurtosis-test xs nil params))
  ([xs kurt {:keys [sides type]
             :or {sides :two-sided type :kurt}}]
   (let [kurt (double (or kurt (kurtosis xs type)))
         n (count xs)
         e (m// (m/* 3.0 (m/dec n)) (m/inc n))
         varb2 (m// (m/* 24.0 n (m/- n 2) (m/- n 3))
                  (m/* (m/sq (m/inc n)) (m/+ n 3) (m/+ n 5)))
         x (m// (m/- kurt e) (m/sqrt varb2))
         sqrtbeta1 (m/* (m// (m/* 6.0 (m/+ (m/* n n) (m/* -5 n) 2))
                         (m/* (m/+ n 7) (m/+ n 9)))
                      (m/sqrt (m// (m/* 6.0 (m/+ n 3) (m/+ n 5))
                                 (m/* n (m/- n 2) (m/- n 3)))))
         a (m/+ 6.0 (m/* (m// 8.0 sqrtbeta1) (m/+ (m// 2.0 sqrtbeta1)
                                          (m/sqrt (m/inc (m// 4.0 (m/* sqrtbeta1 sqrtbeta1)))))))
         term1 (m/- 1.0 (m// 2.0 (m/* 9.0 a)))
         denom (m/inc (m/* x (m/sqrt (m// 2.0 (m/- a 4.0)))))
         term2 (m/* (m/signum denom) (m/cbrt (m// (m/- 1.0 (m// 2.0 a))
                                              (m/abs denom))))
         Z (m// (m/- term1 term2)
              (m/sqrt (m// 2.0 (m/* 9.0 a))))]
     {:p-value (p-value r/default-normal Z sides)
      :Z Z
      :kurtosis kurt})))

(defn normality-test
  "Performs the D'Agostino-Pearson K² omnibus test for normality.

  This test combines the results of the skewness and kurtosis tests to provide
  an overall assessment of whether the sample data deviates from a normal distribution
  in terms of either asymmetry or peakedness/tailedness.

  The test works by:
  1. Calculating a normalized test statistic (Z₁) for skewness using [[skewness-test]].
  2. Calculating a normalized test statistic (Z₂) for kurtosis using [[kurtosis-test]].
  3. Combining these into an omnibus statistic: K² = Z₁² + Z₂².
  4. Under the null hypothesis that the data comes from a normal distribution,
     K² approximately follows a Chi-squared distribution with 2 degrees of freedom.

  Parameters:

  - `xs` (seq of numbers): The sample data.
  - `skew` (double, optional): A pre-calculated skewness value (type `:g1` used by default in underlying test).
  - `kurt` (double, optional): A pre-calculated kurtosis value (type `:kurt` used by default in underlying test).
  - `params` (map, optional): Options map:
    - `:sides` (keyword, default `:one-sided-greater`): Specifies the side(s) of the
      Chi-squared(2) distribution used for p-value calculation.
      - `:one-sided-greater` (default and standard): Tests if K² is significantly large,
        indicating departure from normality in skewness, kurtosis, or both.
      - `:one-sided-less`: Tests if the K² statistic is significantly small.
      - `:two-sided`: Tests if the K² statistic is extreme in either tail.

  Returns a map containing:

  - `:Z`: The calculated K² omnibus test statistic (labeled `:Z` for consistency,
           though it follows Chi-squared(2)).
  - `:stat`: Alias for `:Z`.
  - `:p-value`: The p-value associated with the K² statistic and `:sides`.
  - `:skewness`: The sample skewness value used (either provided or calculated).
  - `:kurtosis`: The sample kurtosis value used (either provided or calculated).

  See also [[skewness-test]], [[kurtosis-test]], [[jarque-bera-test]]."
  ([xs] (normality-test xs nil))
  ([xs params] (normality-test xs nil nil params))
  ([xs skew kurt {:keys [sides]
                  :or {sides :one-sided-greater}}]
   (let [{^double skew-Z :Z skew :skewness} (skewness-test xs skew nil)
         {^double kurt-Z :Z kurt :kurtosis} (kurtosis-test xs kurt nil)
         Z (m/+ (m/* skew-Z skew-Z)
              (m/* kurt-Z kurt-Z))]
     {:p-value (p-value (r/distribution :chi-squared {:degrees-of-freedom 2}) Z sides)
      :Z Z
      :skewness skew
      :kurtosis kurt})))

(defn jarque-bera-test
  "Performs the Jarque-Bera goodness-of-fit test to determine if sample data
  exhibits skewness and kurtosis consistent with a normal distribution.

  The test assesses the null hypothesis that the data comes from a normally
  distributed population (i.e., population skewness is 0 and population excess
  kurtosis is 0).

  The test statistic is calculated as:
  `JB = (n/6) * (S^2 + (1/4)*K^2)`
  where `n` is the sample size, `S` is the sample skewness (using `:g1` type),
  and `K` is the excess kurtosis `:g2`.
  Under the null hypothesis, the JB statistic asymptotically follows a Chi-squared
  distribution with 2 degrees of freedom.

  Parameters:

  - `xs` (seq of numbers): The sample data.
  - `skew` (double, optional): A pre-calculated sample skewness value (type `:g1`).
    If omitted, it's calculated from `xs`.
  - `kurt` (double, optional): A pre-calculated sample *excess* kurtosis value (type `:g2`).
    If omitted, it's calculated from `xs`.
  - `params` (map, optional): Options map:
    - `:sides` (keyword, default `:one-sided-greater`): Specifies the side(s) of the
      Chi-squared(2) distribution used for p-value calculation.
      - `:one-sided-greater` (default and standard for JB): Tests if the JB statistic is
        significantly large, indicating departure from normality.
      - `:one-sided-less`: Tests if the statistic is significantly small.
      - `:two-sided`: Tests if the statistic is extreme in either tail.

  Returns a map containing:

  - `:Z`: The calculated Jarque-Bera test statistic (labeled `:Z` for consistency,
           though it follows Chi-squared(2)).
  - `:stat`: Alias for `:Z`.
  - `:p-value`: The p-value associated with the test statistic and `:sides`, derived
                 from the Chi-squared(2) distribution.
  - `:skewness`: The sample skewness (type `:g1`) used in the calculation.
  - `:kurtosis`: The sample kurtosis (type `:g2`) used in the calculation.

  See also [[skewness-test]], [[kurtosis-test]], [[normality-test]], [[bonett-seier-test]]."
  ([xs] (jarque-bera-test xs nil))
  ([xs params] (jarque-bera-test xs nil nil params))
  ([xs skew kurt {:keys [sides]
                  :or {sides :one-sided-greater}}]
   (let [skew (double (or skew (skewness xs :g1)))
         kurt (double (or kurt (kurtosis xs :g2)))
         n (count xs)
         Z (m/* n m/SIXTH (m/+ (m/* skew skew) (m/* 0.25 kurt kurt)))]
     {:p-value (p-value (r/distribution :chi-squared {:degrees-of-freedom 2}) Z sides)
      :Z Z
      :skewness skew
      :kurtosis kurt})))

(defn bonett-seier-test
  "Performs the Bonett-Seier test for normality based on Geary's 'g' kurtosis measure.

  This test assesses the null hypothesis that the data comes from a normally
  distributed population by checking if the sample Geary's 'g' statistic
  significantly deviates from the value expected under normality (`sqrt(2/pi)`).

  Parameters:

  - `xs` (seq of numbers): The sample data. Requires `(count xs) > 3` for variance calculation.
  - `geary-kurtosis` (double, optional): A pre-calculated Geary's 'g' kurtosis value.
    If omitted, it's calculated from `xs`.
  - `params` (map, optional): Options map:
    - `:sides` (keyword, default `:two-sided`): Specifies the alternative hypothesis
      regarding the deviation from normal kurtosis.
      - `:two-sided` (default): The population kurtosis (measured by 'g') is different from normal.
      - `:one-sided-greater`: Population is leptokurtic ('g' < sqrt(2/pi)). Note Geary's 'g' decreases with peakedness.
      - `:one-sided-less`: Population is platykurtic ('g' > sqrt(2/pi)). Note Geary's 'g' increases with flatness.

  Returns a map containing:

  - `:Z`: The final test statistic (approximately standard normal under H0).
  - `:stat`: Alias for `:Z`.
  - `:p-value`: The p-value associated with `Z` and the specified `:sides`.
  - `:kurtosis`: The Geary's 'g' kurtosis value used in the test.
  - `:n`: The sample size.
  - `:sides`: The alternative hypothesis side used.

  References:
  - Bonett, D. G., & Seier, E. (2002). A test of normality with high uniform power.
    Computational Statistics & Data Analysis, 40(3), 435-445. (Provides theoretical basis)

  See also [[kurtosis]], [[kurtosis-test]], [[normality-test]], [[jarque-bera-test]]."
  ([xs] (bonett-seier-test xs nil))
  ([xs params] (bonett-seier-test xs nil params))
  ([xs geary-kurtosis {:keys [sides] :or {sides :two-sided}}]
   (let [n (count xs)]
     (when (m/<= n 3) (throw (ex-info "Test requires sample size > 3 for variance calculation." {:n n})))
     (let [g (double (or geary-kurtosis (kurtosis xs :geary)))
           omega (m/* -13.29 (m/log g))
           Z (m// (m/* (m/sqrt (m/+ n 2))
                       (m/- omega 3.0)) 3.54)]
       {:p-value (p-value r/default-normal Z sides)
        :stat Z :Z Z
        :kurtosis g
        :n n
        :sides sides}))))

(defn binomial-test
  "Performs an exact test of a simple null hypothesis about the probability of success
  in a Bernoulli experiment, based on the binomial distribution.

  This test assesses the null hypothesis that the true probability of success (`p`)
  in the underlying population is equal to a specified value (default 0.5).

  The function can be called in two ways:

  1. With counts: `(binomial-test number-of-successes number-of-trials params)`
  2. With data: `(binomial-test xs params)`, where `xs` is a sequence of outcomes.
     In this case, the outcomes in `xs` are converted to true/false based on the
     `:true-false-conv` parameter (if provided, otherwise numeric 1s are true),
     and the number of successes and total trials are derived from `xs`.

  Parameters:

  - `number-of-successes` (long): Observed number of successful outcomes.
  - `number-of-trials` (long): Total number of trials.
  - `xs` (sequence): Sample data (used in the alternative call signature).
  - `params` (map, optional): Options map:
    - `:p` (double, default `0.5`): The hypothesized probability of success under the null hypothesis.
    - `:alpha` (double, default `0.05`): Significance level for confidence interval calculation.
    - `:sides` (keyword, default `:two-sided`): Specifies the alternative hypothesis.
      - `:two-sided` (default): True probability `p` is not equal to the hypothesized `p`.
      - `:one-sided-greater`: True probability `p` is greater than the hypothesized `p`.
      - `:one-sided-less`: True probability `p` is less than the hypothesized `p`.
    - `:ci-method` (keyword, default `:asymptotic`): Method used to calculate the confidence interval for the probability of success. See [[binomial-ci]] and [[binomial-ci-methods]] for available options (e.g., `:wilson`, `:clopper-pearson`).
    - `:true-false-conv` (optional, used only with `xs`): A function, set, or map to convert elements of `xs` into boolean `true` (success) or `false` (failure). See [[binary-measures-all]] documentation for details. If `nil` and `xs` contains numbers, `1.0` is treated as success.

  Returns a map containing:

  - `:p-value`: The probability of observing a result as extreme as, or more extreme than, the observed number of successes, assuming the null hypothesis is true. Calculated using the binomial distribution.
  - `:p`: The hypothesized probability of success used in the test.
  - `:successes`: The observed number of successes.
  - `:trials`: The total number of trials.
  - `:alpha`: Significance level used for the confidence interval.
  - `:level`: Confidence level (`1 - alpha`).
  - `:sides` / `:test-type`: Alternative hypothesis side used.
  - `:stat`: The test statistic (the observed number of successes).
  - `:estimate`: The observed proportion of successes (`successes / trials`).
  - `:ci-method`: Confidence interval method used.
  - `:confidence-interval`: A confidence interval for the true probability of success, calculated using the specified `:ci-method` and adjusted for the `:sides` parameter."
  ([xs] (binomial-test xs {}))
  ([xs maybe-params]
   (if (map? maybe-params)
     (let [{:keys [true-false-conv] :as params} maybe-params
           xxs (binary/binary-process-list xs true-false-conv)
           nos (count (filter identity xxs))
           not (count xxs)]
       (binomial-test nos not params))
     (binomial-test xs maybe-params {})))
  ([^long number-of-successes ^long number-of-trials {:keys [^double alpha ^double p ci-method sides]
                                                      :or {alpha 0.05 p 0.5
                                                           ci-method :asymptotic sides :two-sided}}]
   (let [distr (r/distribution :binomial {:trials number-of-trials :p p})]
     {:p-value (p-value distr (double number-of-successes) sides)
      :p p
      :successes number-of-successes
      :trials number-of-trials
      :alpha alpha
      :level (m/- 1.0 alpha)
      :test-type sides
      :stat number-of-successes
      :estimate (m// (double number-of-successes) number-of-trials)
      :ci-method ci-method
      :confidence-interval (let [bci (partial binomial-ci number-of-successes number-of-trials ci-method)]
                             (sides-case sides
                                         (vec (butlast (bci (m/- 1.0 alpha))))
                                         [(first (bci (m/- 1.0 (m/* alpha 2.0)))) 1.0]
                                         [0.0 (second (bci (m/- 1.0 (m/* alpha 2.0))))]))})))

;; t/z

(defn- test-update-ci
  [^double mu ^double stderr [^double l ^double r]]
  [(m/+ mu (m/* l stderr)) (m/+ mu (m/* r stderr))])

(defn- test-pvalue-ci
  [d sides ^double stat ^double alpha]
  {:confidence-interval (sides-case sides
                                    (let [cint (double (r/icdf d (m/- 1.0 (m/* 0.5 alpha))))]
                                      [(m/- stat cint) (m/+ stat cint)])
                                    [(m/- stat (double (r/icdf d (m/- 1.0 alpha)))) ##Inf]
                                    [##-Inf (m/+ stat (double (r/icdf d (m/- 1.0 alpha))))])
   :p-value (p-value d stat sides)})

(defn- test-one-sample
  [xs {:keys [^double alpha sides ^double mu]
       :or {alpha 0.05 sides :two-sided mu 0.0}}]
  (let [axs (m/seq->double-array xs)
        n (alength axs)
        m (mean axs)
        v (variance axs)
        stderr (m/sqrt (m// v n))]
    (when (m/< stderr (m/* 10.0 m/MACHINE-EPSILON (m/abs m))) (throw (ex-info "Constant data, can't perform test." {:stderr stderr :mean m})))
    {:n n
     :estimate m
     :mu mu
     :stat (m// (m/- m mu) stderr)
     :test-type sides
     :stderr stderr
     :alpha alpha
     :level (m/- 1.0 alpha)}))

(defn t-test-one-sample
  "Performs a one-sample Student's t-test to compare the sample mean against a hypothesized population mean.

  This test assesses the null hypothesis that the true population mean is equal to `mu`.
  It is suitable when the population standard deviation is unknown and is estimated
  from the sample.

  Parameters:

  - `xs` (seq of numbers): The sample data.
  - `params` (map, optional): Options map:
    - `:alpha` (double, default `0.05`): Significance level for the confidence interval.
    - `:sides` (keyword, default `:two-sided`): Specifies the alternative hypothesis.
      - `:two-sided` (default): The true mean is not equal to `mu`.
      - `:one-sided-greater`: The true mean is greater than `mu`.
      - `:one-sided-less`: The true mean is less than `mu`.
    - `:mu` (double, default `0.0`): The hypothesized population mean under the null hypothesis.

  Returns a map containing:

  - `:t`: The calculated t-statistic.
  - `:stat`: Alias for `:t`.
  - `:df`: Degrees of freedom (`n-1`).
  - `:p-value`: The p-value associated with the t-statistic and `:sides`.
  - `:confidence-interval`: Confidence interval for the true population mean.
  - `:estimate`: The calculated sample mean.
  - `:n`: The sample size.
  - `:mu`: The hypothesized population mean used in the test.
  - `:stderr`: The standard error of the mean (calculated from the sample).
  - `:alpha`: Significance level used.
  - `:sides`: Alternative hypothesis side used.
  - `:test-type`: Alias for `:sides`.

  Assumptions:

  - The data are independent observations.
  - The data are drawn from a population that is approximately normally distributed.
    (The t-test is relatively robust to moderate violations, especially with larger sample sizes).

  See also [[z-test-one-sample]] for large samples or known population standard deviation."
  ([xs] (t-test-one-sample xs {}))
  ([xs m]
   (let [{:keys [^long n ^double stat test-type ^double alpha ^double mu ^double stderr]
          :as res} (test-one-sample xs m)
         df (m/dec n)
         pvals (-> (test-pvalue-ci (r/distribution :t {:degrees-of-freedom df}) test-type stat alpha)
                   (update :confidence-interval (partial test-update-ci mu stderr)))]
     (assoc (merge pvals res) :df df :t stat))))

(def ^{:deprecated "Use [[t-test-one-sample]]"} ttest-one-sample t-test-one-sample)

(defn z-test-one-sample
  "Performs a one-sample Z-test to compare the sample mean against a hypothesized population mean.

  This test assesses the null hypothesis that the true population mean is equal to `mu`.
  It typically assumes either a known population standard deviation or relies on a
  large sample size (e.g., n > 30) where the sample standard deviation provides a
  reliable estimate. This implementation uses the sample standard deviation to calculate
  the standard error.

  Parameters:

  - `xs` (seq of numbers): The sample data.
  - `params` (map, optional): Options map:
    - `:alpha` (double, default `0.05`): Significance level for the confidence interval.
    - `:sides` (keyword, default `:two-sided`): Specifies the alternative hypothesis.
      - `:two-sided` (default): The true mean is not equal to `mu`.
      - `:one-sided-greater`: The true mean is greater than `mu`.
      - `:one-sided-less`: The true mean is less than `mu`.
    - `:mu` (double, default `0.0`): The hypothesized population mean under the null hypothesis.

  Returns a map containing:

  - `:z`: The calculated Z-statistic.
  - `:stat`: Alias for `:z`.
  - `:p-value`: The p-value associated with the Z-statistic and the specified `:sides`.
  - `:confidence-interval`: Confidence interval for the true population mean.
  - `:estimate`: The calculated sample mean.
  - `:n`: The sample size.
  - `:mu`: The hypothesized population mean used in the test.
  - `:stderr`: The standard error of the mean (calculated using sample standard deviation).
  - `:alpha`: Significance level used.
  - `:sides`: Alternative hypothesis side used.
  - `:test-type`: Alias for `:sides`.

  See also [[t-test-one-sample]] for smaller samples or when the population standard deviation is unknown."
  ([xs] (z-test-one-sample xs {}))
  ([xs m]
   (let [{:keys [^double stat test-type ^double alpha ^double mu ^double stderr]
          :as res} (test-one-sample xs m)
         pvals (-> (test-pvalue-ci r/default-normal test-type stat alpha)
                   (update :confidence-interval (partial test-update-ci mu stderr)))]
     (assoc (merge pvals res) :z stat))))

(defn- test-equal-variances
  [^double nx ^double ny ^double vx ^double vy]
  (let [df (m/- (m/+ nx ny) 2.0)
        v (m// (m/+ (m/* vx (m/dec nx))
                (m/* vy (m/dec ny))) df)]
    [df (m/sqrt (m/* v (m/+ (m// 1.0 nx)
                        (m// 1.0 ny))))]))

(defn- test-not-equal-variances
  [^double nx ^double ny ^double vx ^double vy]
  (let [stderrx (m/sqrt (m// vx nx))
        stderry (m/sqrt (m// vy ny))
        stderr (m/hypot-sqrt stderrx stderry)
        df (m// (m/sq (m/sq stderr))
              (m/+ (m// (m/sq (m/sq stderrx)) (m/dec nx))
                 (m// (m/sq (m/sq stderry)) (m/dec ny))))]
    [df stderr]))

(defn- test-two-samples-not-paired
  [xs ys {:keys [^double alpha sides ^double mu equal-variances?]
          :or {alpha 0.05 sides :two-sided mu 0.0 equal-variances? false}}]
  (let [axs (m/seq->double-array xs)
        ays (m/seq->double-array ys)
        nx (alength axs)
        ny (alength ays)
        mx (mean axs)
        my (mean ays)
        vx (variance axs)
        vy (variance ays)
        [df ^double stderr] (if equal-variances?
                              (test-equal-variances nx ny vx vy)
                              (test-not-equal-variances nx ny vx vy))]
    {:n [nx ny] :nx nx :ny ny
     :estimated-mu [mx my]
     :mu mu
     :estimate (m/- mx my mu)
     :stat (m// (m/- mx my mu) stderr)
     :sides sides
     :test-type sides
     :stderr stderr
     :alpha alpha
     :level (m/- 1.0 alpha)
     :df df
     :paired? false
     :equal-variances? equal-variances?}))

(defn t-test-two-samples
  "Performs a two-sample Student's t-test to compare the means of two samples.

  This function can perform:

  - An **unpaired t-test** (assuming independent samples) using either:
    - **Welch's t-test** (default: `:equal-variances? false`): Does not assume equal population variances. Uses the Satterthwaite approximation for degrees of freedom. Recommended unless variances are known to be equal.
    - **Student's t-test** (`:equal-variances? true`): Assumes equal population variances and uses a pooled variance estimate.
  - A **paired t-test** (`:paired? true`): Assumes observations in `xs` and `ys` are paired (e.g., before/after measurements on the same subjects). This performs a one-sample t-test on the differences between paired observations.

  The test assesses the null hypothesis that the true difference between the population
  means (or the mean of the differences for paired test) is equal to `mu`.

  Parameters:

  - `xs` (seq of numbers): The first sample.
  - `ys` (seq of numbers): The second sample.
  - `params` (map, optional): Options map:
    - `:alpha` (double, default `0.05`): Significance level for the confidence interval.
    - `:sides` (keyword, default `:two-sided`): Specifies the alternative hypothesis.
      - `:two-sided` (default): The true difference in means is not equal to `mu`.
      - `:one-sided-greater`: The true difference (`mean(xs) - mean(ys)` or `mean(diff)`) is greater than `mu`.
      - `:one-sided-less`: The true difference (`mean(xs) - mean(ys)` or `mean(diff)`) is less than `mu`.
    - `:mu` (double, default `0.0`): The hypothesized difference in means under the null hypothesis.
    - `:paired?` (boolean, default `false`): If `true`, performs a paired t-test (requires `xs` and `ys` to have the same length). If `false`, performs an unpaired test.
    - `:equal-variances?` (boolean, default `false`): Used only when `paired?` is `false`. If `true`, assumes equal population variances (Student's). If `false`, does not assume equal variances (Welch's).

  Returns a map containing:

  - `:t`: The calculated t-statistic.
  - `:stat`: Alias for `:t`.
  - `:df`: Degrees of freedom used for the t-distribution.
  - `:p-value`: The p-value associated with the t-statistic and `:sides`.
  - `:confidence-interval`: Confidence interval for the true difference in means.
  - `:estimate`: The observed difference between sample means (`mean(xs) - mean(ys)` or `mean(differences)`).
  - `:n`: Sample sizes as `[count xs, count ys]` (or `count diffs` if paired).
  - `:nx`: Sample size of `xs` (if unpaired).
  - `:ny`: Sample size of `ys` (if unpaired).
  - `:estimated-mu`: Observed sample means as `[mean xs, mean ys]` (if unpaired).
  - `:mu`: The hypothesized difference under the null hypothesis.
  - `:stderr`: The standard error of the difference between the means (or of the mean difference if paired).
  - `:alpha`: Significance level used.
  - `:sides`: Alternative hypothesis side used.
  - `:test-type`: Alias for `:sides`.
  - `:paired?`: Boolean indicating if a paired test was performed.
  - `:equal-variances?`: Boolean indicating the variance assumption used (if unpaired).

  Assumptions:
  - Independence of observations (within and between groups for unpaired).
  - Normality of the underlying populations (or of the differences for paired). The t-test is relatively robust to violations of normality, especially with larger sample sizes.
  - Equal variances (only if `:equal-variances? true`)."
  ([xs ys] (t-test-two-samples xs ys {}))
  ([xs ys {:keys [paired? equal-variances?]
           :or {paired? false equal-variances? false}
           :as params}]
   (let [nx (count xs)
         ny (count ys)]
     (when-not (or (and equal-variances? (m/< 2 (m/+ nx ny)) (m/pos? nx) (m/pos? ny))
                   (and (not equal-variances?)
                        (m/> nx 1) (m/> ny 1))) (throw (ex-info "Not enough observations." {:nx nx :ny ny :equal-variances? equal-variances?})))
     (when (and paired? (m/not== nx ny)) (throw (ex-info "Lengths of xs and ys should be equal." {:nx nx :ny ny})))
     (if paired?
       (-> (t-test-one-sample (map - xs ys) params)
           (assoc :paired? true))
       (let [{:keys [test-type ^double stat ^double alpha ^double df ^double mu ^double stderr]
              :as res} (test-two-samples-not-paired xs ys params)
             pvals (-> (test-pvalue-ci (r/distribution :t {:degrees-of-freedom df}) test-type stat alpha)
                       (update :confidence-interval (partial test-update-ci mu stderr)))]
         (assoc (merge pvals res) :t stat))))))

(def ^{:deprecated "Use [[t-test-two-samples]]"} ttest-two-samples t-test-two-samples)

(defn z-test-two-samples
  "Performs a two-sample Z-test to compare the means of two independent or paired samples.

  This test assesses the null hypothesis that the difference between the population
  means is equal to `mu` (default 0). It typically assumes known population variances
  or relies on large sample sizes where sample variances provide good estimates.
  This implementation calculates the standard error using the provided sample variances.

  Parameters:

  - `xs` (seq of numbers): The first sample.
  - `ys` (seq of numbers): The second sample.
  - `params` (map, optional): Options map:
    - `:alpha` (double, default `0.05`): Significance level for the confidence interval.
    - `:sides` (keyword, default `:two-sided`): Specifies the alternative hypothesis.
      - `:two-sided` (default): The true difference in means is not equal to `mu`.
      - `:one-sided-greater`: The true difference in means (`mean(xs) - mean(ys)`) is greater than `mu`.
      - `:one-sided-less`: The true difference in means (`mean(xs) - mean(ys)`) is less than `mu`.
    - `:mu` (double, default `0.0`): The hypothesized difference in means under the null hypothesis.
    - `:paired?` (boolean, default `false`): If `true`, performs a paired Z-test by applying [[z-test-one-sample]] to the differences between paired observations in `xs` and `ys` (requires `xs` and `ys` to have the same length). If `false`, performs a two-sample test assuming independence.
    - `:equal-variances?` (boolean, default `false`): Used only when `paired?` is `false`. If `true`, assumes population variances are equal and calculates a pooled standard error. If `false`, calculates the standard error without assuming equal variances (Welch's approach adapted for Z-test). This affects the standard error calculation but the standard normal distribution is still used for inference.

  Returns a map containing:

  - `:z`: The calculated Z-statistic.
  - `:stat`: Alias for `:z`.
  - `:p-value`: The p-value associated with the Z-statistic and the specified `:sides`.
  - `:confidence-interval`: Confidence interval for the true difference in means.
  - `:estimate`: The observed difference between sample means (`mean(xs) - mean(ys)`).
  - `:n`: Sample sizes as `[count xs, count ys]`.
  - `:nx`: Sample size of `xs`.
  - `:ny`: Sample size of `ys`.
  - `:estimated-mu`: The observed sample means as `[mean xs, mean ys]`.
  - `:mu`: The hypothesized difference under the null hypothesis.
  - `:stderr`: The standard error of the difference between the means.
  - `:alpha`: Significance level used.
  - `:sides`: Alternative hypothesis side used.
  - `:test-type`: Alias for `:sides`.
  - `:paired?`: Boolean indicating if a paired test was performed.
  - `:equal-variances?`: Boolean indicating the assumption used for standard error calculation (if unpaired).

  See also [[t-test-two-samples]] for smaller samples or when population variances are unknown."
  ([xs ys] (z-test-two-samples xs ys {}))
  ([xs ys {:keys [paired? equal-variances?]
           :or {paired? false equal-variances? false}
           :as params}]
   (let [nx (count xs)
         ny (count ys)]
     (when-not (or (and equal-variances? (m/< 2 (m/+ nx ny)) (m/pos? nx) (m/pos? ny))
                   (and (not equal-variances?)
                        (m/> nx 1) (m/> ny 1))) (throw (ex-info "Not enough observations." {:nx nx :ny ny :equal-variances? equal-variances?})))
     (when (and paired? (m/not== nx ny)) (throw (ex-info "Lengths of xs and ys should be equal." {:nx nx :ny ny})))
     (if paired?
       (-> (z-test-one-sample (map - xs ys) params)
           (assoc :paired? true))
       (let [{:keys [test-type ^double stat ^double alpha ^double ^double mu ^double stderr]
              :as res} (test-two-samples-not-paired xs ys params)
             pvals (-> (test-pvalue-ci r/default-normal test-type stat alpha)
                       (update :confidence-interval (partial test-update-ci mu stderr)))]
         (-> (merge res pvals)
             (assoc :z stat)
             (dissoc :df)))))))

(defn f-test
  "Performs an F-test to compare the variances of two independent samples.

  The test assesses the null hypothesis that the variances of the populations
  from which `xs` and `ys` are drawn are equal.

  Assumes independence of samples. The test is sensitive to departures from
  the assumption that both populations are normally distributed.

  Parameters:

  - `xs` (seq of numbers): The first sample.
  - `ys` (seq of numbers): The second sample.
  - `params` (map, optional): Options map:
    - `:sides` (keyword, default `:two-sided`): Specifies the alternative hypothesis
      regarding the ratio of variances (Var(xs) / Var(ys)).
      - `:two-sided` (default): Variances are not equal (ratio != 1).
      - `:one-sided-greater`: Variance of `xs` is greater than variance of `ys` (ratio > 1).
      - `:one-sided-less`: Variance of `xs` is less than variance of `ys` (ratio < 1).
    - `:alpha` (double, default `0.05`): Significance level for the confidence interval.

  Returns a map containing:

  - `:F`: The calculated F-statistic (ratio of sample variances: Var(xs) / Var(ys)).
  - `:stat`: Alias for `:F`.
  - `:estimate`: Alias for `:F`, representing the estimated ratio of variances.
  - `:df`: Degrees of freedom as `[numerator-df, denominator-df]`, corresponding to `[(count xs)-1, (count ys)-1]`.
  - `:n`: Sample sizes as `[count xs, count ys]`.
  - `:nx`: Sample size of `xs`.
  - `:ny`: Sample size of `ys`.
  - `:sides`: The alternative hypothesis side used (`:two-sided`, `:one-sided-greater`, or `:one-sided-less`).
  - `:test-type`: Alias for `:sides`.
  - `:p-value`: The p-value associated with the F-statistic and the specified `:sides`.
  - `:confidence-interval`: A confidence interval for the true ratio of the population variances (Var(xs) / Var(ys))."
  ([xs ys] (f-test xs ys {}))
  ([xs ys {:keys [sides ^double alpha]
           :or {sides :two-sided alpha 0.05}}]
   (let [nx (count xs)
         ny (count ys)
         dfx (m/dec nx)
         dfy (m/dec ny)
         F (m// (variance xs) (variance ys))
         distr (r/distribution :f {:denominator-degrees-of-freedom dfy
                                   :numerator-degrees-of-freedom dfx})]
     {:F F
      :stat F
      :estimate F
      :df [dfx dfy]
      :n [nx ny] :nx nx :ny ny
      :sides sides
      :test-type sides
      :p-value (p-value distr F sides)
      :confidence-interval (sides-case sides
                                       [(m// F (double (r/icdf distr (m/- 1.0 (m/* alpha 0.5)))))
                                        (m// F (double (r/icdf distr (m/* alpha 0.5))))]
                                       [(m// F (double (r/icdf distr (m/- 1.0 alpha)))) ##Inf]
                                       [0.0 (m// F (double (r/icdf distr alpha)))])})))

(defn- pdt-gof
  "Goodness of fit"
  [xs p ^double lambda]
  (let [cnt (count xs)
        n (sum xs)
        df (m/dec cnt)
        p (or p (repeat cnt 1.0))
        psum (sum p)
        p (map (fn [^double p] (m// p psum)) p)
        xhat (map (fn [^double p] (m/* n p)) p)
        stat (condp = lambda
               0.0 (m/* 2.0 (sum (map (fn [^long a ^double b]
                                      (m/* a (m/- (m/log a) (m/log b)))) xs xhat)))
               -1.0 (m/* 2.0 (sum (map (fn [^double a ^long b]
                                       (m/* a (m/- (m/log a) (m/log b)))) xhat xs)))
               (m/* (m// 2.0 (m/* lambda (m/inc lambda)))
                  (sum (map (fn [^long a ^double b]
                              (m/* a (m/dec (m/pow (m// a b) lambda)))) xs xhat))))]
    {:stat stat :df df :n n :expected xhat :p p
     :estimate (map (fn [^double v] (m// v n)) xs)}))

(defn- pdt-distribution
  [xs distr bins ^double lambda]
  (let [[counts probabilities] (quantize-distribution xs distr bins)]
    (pdt-gof counts probabilities lambda)))

(defn- pdt-bootstrap-ci
  [{:keys [estimate ^long n ^double alpha ci-sides]} samples]
  (let [alpha (if (#{:both :two-sided} ci-sides) alpha (m/* alpha 2.0))
        d (r/distribution :multinomial {:trials n :ps (if (map? estimate) (vals estimate) estimate)})
        rands (apply map (comp m/seq->double-array vector) (r/->seq d samples))
        vs (sides-case ci-sides
                       (let [qs [(m// alpha 2.0) (m/- 1.0 (m// alpha 2.0))]]
                         (map #(v/div (quantiles % qs) (double n)) rands))
                       (let [q (m// alpha 2.0)]
                         (map #(vector (m// (quantile % q) n) 1.0) rands))
                       (let [q (m/- 1.0 (m// alpha 2.0))]
                         (map #(vector 0.0 (m// (quantile % q) n)) rands)))]
    (if (map? estimate) (zipmap (keys estimate) vs) vs)))

(defn- pdt-multi
  [ct ^double lambda]
  (let [xs (infer-ct ct)
        {:keys [^long n cols rows]} (contingency-table->marginals xs)
        n (double n)
        xhat (->> (for [[k1 ^long v1] rows
                        [k2 ^long v2] cols]
                    [[k1 k2] (m// (m/* v1 v2) n)])
                  (into {}))
        n1 (count (map first rows))
        n2 (count (map second cols))
        df (m/* (m/dec n1) (m/dec n2))
        stat (condp = lambda
               0.0 (m/* 2.0 (double (reduce (fn [^double sum [k ^long cnt]]
                                              (m/+ sum (m/* cnt (m/- (m/log cnt) (m/log (xhat k)))))) 0.0 xs)))
               -1.0 (m/* 2.0 (double (reduce (fn [^double sum [k ^double xhv]]
                                               (let [^double cnt (get xs k 0.0)]
                                                 (m/+ sum (m/* xhv (m/- (m/log xhv) (m/log cnt) ))))) 0.0 xhat)))
               (m/* (m// 2.0 (m/* lambda (m/inc lambda)))
                    (double (reduce (fn [^double sum [k ^long cnt]]
                                      (m/+ sum (m/* cnt (m/dec (m/pow (m// cnt ^double (xhat k)) lambda))))) 0.0 xs))))]
    {:stat stat :df df :n n :k n1 :r n2
     :expected xhat
     :estimate (into {} (map (fn [[k ^long v]] [k (m// v n)]) xs))}))

(defn power-divergence-test
  "Performs a power divergence test, which encompasses several common statistical tests
  like Chi-squared, G-test (likelihood ratio), etc., based on the lambda parameter.
  This function can perform either a goodness-of-fit test or a test for independence
  in a contingency table.

  Usage:

  1.  **Goodness-of-Fit (GOF):**
      - Input: `observed-counts` (sequence of numbers) and `:p` (expected probabilities/weights).
      - Input: `data` (sequence of numbers) and `:p` (a distribution object).
        In this case, a histogram of `data` is created (controlled by `:bins`) and
        compared against the probability mass/density of the distribution in those bins.

  2.  **Test for Independence:**
      - Input: `contingency-table` (2D sequence or map format). The `:p` option is ignored.

  Options map:

  * `:lambda` (double, default: `2/3`): Determines the specific test statistic. Common values:
      * `1.0`: Pearson Chi-squared test ([[chisq-test]]).
      * `0.0`: G-test / Multinomial Likelihood Ratio test ([[multinomial-likelihood-ratio-test]]).
      * `-0.5`: Freeman-Tukey test ([[freeman-tukey-test]]).
      * `-1.0`: Minimum Discrimination Information test ([[minimum-discrimination-information-test]]).
      * `-2.0`: Neyman Modified Chi-squared test ([[neyman-modified-chisq-test]]).
      * `2/3`: Cressie-Read test (default, [[cressie-read-test]]).
  * `:p` (seq of numbers or distribution): Expected probabilities/weights (for GOF with counts)
    or a `fastmath.random` distribution object (for GOF with data). Ignored for independence tests.
  * `:alpha` (double, default: `0.05`): Significance level for confidence intervals.
  * `:ci-sides` (keyword, default: `:two-sided`): Sides for bootstrap confidence intervals
    (`:two-sided`, `:one-sided-greater`, `:one-sided-less`).
  * `:sides` (keyword, default: `:one-sided-greater`): Alternative hypothesis side for the p-value calculation
    against the Chi-squared distribution (`:one-sided-greater`, `:one-sided-less`, `:two-sided`).
  * `:bootstrap-samples` (long, default: `1000`): Number of bootstrap samples for confidence interval estimation.
  * `:ddof` (long, default: `0`): Delta degrees of freedom. Adjustment subtracted from the calculated degrees of freedom.
  * `:bins` (number, keyword, or seq): Used only for GOF test against a distribution.
    Specifies the number of bins, an estimation method (see [[histogram]]), or explicit bin edges for histogram creation.

  Returns a map containing:

  - `:stat`: The calculated power divergence test statistic.
  - `:chi2`: Alias for `:stat`.
  - `:df`: Degrees of freedom for the test.
  - `:p-value`: The p-value associated with the test statistic.
  - `:n`: Total number of observations.
  - `:estimate`: Observed proportions.
  - `:expected`: Expected counts or proportions under the null hypothesis.
  - `:confidence-interval`: Bootstrap confidence intervals for the observed proportions.
  - `:lambda`, `:alpha`, `:sides`, `:ci-sides`: Input options used."
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {}))
  ([contingency-table-or-xs {:keys [^double lambda ci-sides sides p ^double alpha ^long bootstrap-samples
                                    ^long ddof bins]
                             :or {lambda m/TWO_THIRD sides :one-sided-greater ci-sides :two-sided
                                  alpha 0.05 bootstrap-samples 1000 ddof 0}}]
   (let [{:keys [df stat] :as res} (-> (cond
                                         (and p (r/distribution? p))
                                         (pdt-distribution contingency-table-or-xs p bins lambda)

                                         (and (sequential? contingency-table-or-xs)
                                              (every? number? contingency-table-or-xs))
                                         (pdt-gof contingency-table-or-xs p lambda)

                                         :else (pdt-multi contingency-table-or-xs lambda))
                                       (update :df (fn [^long df] (m/- df ddof))))
         distr (r/distribution :chi-squared {:degrees-of-freedom df})
         res (assoc res :lambda lambda :sides sides :test-type sides :ci-sides ci-sides :chi2 stat :alpha alpha :level (m/- 1.0 alpha)
                    :p-value (p-value distr stat sides))]
     (assoc res :confidence-interval (pdt-bootstrap-ci res bootstrap-samples)))))

(defn chisq-test
  "Chi square test, a power divergence test for `lambda` 1.0"
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda 1.0}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda 1.0))))

(defn multinomial-likelihood-ratio-test
  "Multinomial likelihood ratio test, a power divergence test for `lambda` 0.0"
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda 0.0}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda 0.0))))

(defn minimum-discrimination-information-test
  "Minimum discrimination information test, a power divergence test for `lambda` -1.0"
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda -1.0}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda -1.0))))

(defn neyman-modified-chisq-test
  "Neyman modifield chi square test, a power divergence test for `lambda` -2.0"
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda -2.0}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda -2.0))))

(defn freeman-tukey-test
  "Freeman-Tukey test, a power divergence test for `lambda` -0.5"
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda -0.5}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda -0.5))))

(defn cressie-read-test
  "Cressie-Read test, a power divergence test for `lambda` 2/3"
  ([contingency-table-or-xs] (power-divergence-test contingency-table-or-xs {:lambda m/TWO_THIRD}))
  ([contingency-table-or-xs params]
   (power-divergence-test contingency-table-or-xs (assoc params :lambda m/TWO_THIRD))))

;; copy docs
(doseq [v [#'chisq-test #'multinomial-likelihood-ratio-test #'minimum-discrimination-information-test
           #'neyman-modified-chisq-test #'freeman-tukey-test #'cressie-read-test]]
  (alter-meta! v update :doc str "\n\n" (:doc (meta #'power-divergence-test))))

;;

(defn- anova
  [xss]
  (let [Ni (map count xss)
        Zi (map mean xss)
        Z (m// (sum (mapcat identity xss)) (sum Ni))
        SSt (sum (map (fn [^double n ^double z]
                        (m/* n (m/sq (m/- z Z)))) Ni Zi))
        SSe (sum (map (fn [xs ^double zi]
                        (v/magsq (map (fn [^double v]
                                        (m/- v zi)) xs))) xss Zi))
        k (count Ni)
        DFt (m/dec k)
        DFe (m/- (sum Ni) k)
        MSe (m// SSe DFe)]
    {:n Ni :SSt SSt :SSe SSe :DFt DFt :DFe (int DFe) :MSt (m// SSt DFt) :MSe MSe}))

(defn- update-f-p-value
  [{:keys [DFt DFe ^double MSt ^double MSe] :as aov} sides]
  (let [F (m// MSt MSe)
        distr (r/distribution :f {:numerator-degrees-of-freedom DFt
                                  :denominator-degrees-of-freedom DFe})]
    (assoc aov
           :F F :df [DFt DFe] :stat F
           :p-value (p-value distr F sides))))

(defn one-way-anova-test
  "Performs a one-way analysis of variance (ANOVA) test.

  ANOVA tests the null hypothesis that the means of two or more independent groups
  are equal. It assumes that the data within each group are normally distributed
  and have equal variances.

  Parameters:

  - `xss` (sequence of sequences): A collection where each element is a sequence
    representing a group of observations.
  - `params` (map, optional): Options map with the following key:
    - `:sides` (keyword, default `:one-sided-greater`): Alternative hypothesis side for the F-test.
      Possible values: `:one-sided-greater`, `:one-sided-less`, `:two-sided`.

  Returns a map containing:

  - `:F`: The F-statistic for the test.
  - `:stat`: Alias for `:F`.
  - `:p-value`: The p-value for the test.
  - `:df`: Degrees of freedom for the F-statistic ([DFt, DFe]).
  - `:n`: Sequence of sample sizes for each group.
  - `:SSt`: Sum of squares between groups (treatment).
  - `:SSe`: Sum of squares within groups (error).
  - `:DFt`: Degrees of freedom between groups.
  - `:DFe`: Degrees of freedom within groups.
  - `:MSt`: Mean square between groups.
  - `:MSe`: Mean square within groups.
  - `:sides`: Test side used."
  ([xss] (one-way-anova-test xss {}))
  ([xss {:keys [sides]
         :or {sides :one-sided-greater}}]
   (update-f-p-value (anova xss) sides)))

(defn levene-test
  "Performs Levene's test for homogeneity of variances across two or more groups.

  Levene's test assesses the null hypothesis that the variances of the groups are equal.
  It calculates an ANOVA on the absolute deviations of the data points from their group
  center (mean by default).

  Parameters:

  - `xss` (sequence of sequences): A collection where each element is a sequence representing a group of observations.
  - `params` (map, optional): Options map with the following keys:
    - `:sides` (keyword, default `:one-sided-greater`): Alternative hypothesis side for the F-test.
      Possible values: `:one-sided-greater`, `:one-sided-less`, `:two-sided`.
    - `:statistic` (fn, default [[mean]]): Function to calculate the center of each group (e.g., [[mean]], [[median]]). Using [[median]] results in the Brown-Forsythe test.
    - `:scorediff` (fn, default [[abs]]): Function applied to the difference between each data point and its group center (e.g., [[abs]], [[sq]]).

  Returns a map containing:

  - `:W`: The Levene test statistic (which is an F-statistic).
  - `:stat`: Alias for `:W`.
  - `:p-value`: The p-value for the test.
  - `:df`: Degrees of freedom for the F-statistic ([DFt, DFe]).
  - `:n`: Sequence of sample sizes for each group.
  - `:SSt`: Sum of squares between groups (treatment).
  - `:SSe`: Sum of squares within groups (error).
  - `:DFt`: Degrees of freedom between groups.
  - `:DFe`: Degrees of freedom within groups.
  - `:MSt`: Mean square between groups.
  - `:MSe`: Mean square within groups.
  - `:sides`: Test side used.

  See also [[brown-forsythe-test]]."
  ([xss] (levene-test xss {}))
  ([xss {:keys [sides statistic scorediff]
         :or {sides :one-sided-greater statistic mean scorediff abs}}]
   (let [res (update-f-p-value (anova (map (fn [xs]
                                             (let [s (double (statistic xs))]
                                               (map (fn [^double v]
                                                      (scorediff (m/- v s))) xs))) xss)) sides)]
     (-> (assoc res :W (:F res))
         (dissoc :F)))))

(defn brown-forsythe-test
  "Brown-Forsythe test for homogeneity of variances.

  This test is a modification of Levene's test, using the median instead of the mean
  for calculating the spread within each group. This makes the test more robust
  against non-normally distributed data.

  Calls [[levene-test]] with `:statistic` set to [[median]]. Accepts the same parameters
  as [[levene-test]], except for `:statistic`.

  Parameters:
  - `xss` (sequence of sequences): A collection of data groups.
  - `params` (map, optional): Options map (see [[levene-test]])."
  ([xss] (levene-test xss {:statistic median}))
  ([xss params] (levene-test xss (assoc params :statistic median))))

(defn fligner-killeen-test
  "Performs the Fligner-Killeen test for homogeneity of variances across two or more groups.

  The Fligner-Killeen test is a non-parametric test that assesses the null hypothesis
  that the variances of the groups are equal. It is robust against departures from normality.
  The test is based on ranks of the absolute deviations from the group medians.

  Parameters:
  
  - `xss` (sequence of sequences): A collection where each element is a sequence representing a group of observations.
  - `params` (map, optional): Options map with the following key:
    - `:sides` (keyword, default `:one-sided-greater`): Alternative hypothesis side for the Chi-squared test.
      Possible values: `:one-sided-greater`, `:one-sided-less`, `:two-sided`.

  Returns a map containing:

  - `:chi2`: The Fligner-Killeen test statistic (Chi-squared value).
  - `:stat`: Alias for `:chi2`.
  - `:p-value`: The p-value for the test.
  - `:df`: Degrees of freedom for the test (number of groups - 1).
  - `:n`: Sequence of sample sizes for each group.
  - `:SSt`: Sum of squares between groups (treatment) based on transformed ranks.
  - `:SSe`: Sum of squares within groups (error) based on transformed ranks.
  - `:DFt`: Degrees of freedom between groups.
  - `:DFe`: Degrees of freedom within groups.
  - `:MSt`: Mean square between groups.
  - `:MSe`: Mean square within groups.
  - `:sides`: Test side used."
  ([xss] (fligner-killeen-test xss {}))
  ([xss {:keys [sides]
         :or {sides :one-sided-greater}}]
   (let [Z (mapcat (fn [xs]
                     (let [s (median xs)]
                       (map (fn [^double v] (abs (m/- v s))) xs))) xss)
         ranks (m/rank Z)
         rden (m// (m/* 2.0 (m/inc (count ranks))))
         qij (map (fn [^double r]
                    (r/icdf r/default-normal (m/+ 0.5 (m/* rden (m/inc r))))) ranks)
         {:keys [^double SSt ^double SSe ^int DFt ^int DFe]
          :as res} (->> (map count xss)
                        (reductions m/+ 0)
                        (partition 2 1)
                        (map (fn [[d t]] (drop d (take t qij))))
                        (anova))
         y (m// SSt SSe)
         chi2 (m// (m/* y (m/+ DFt DFe)) (m/inc y))
         distr (r/distribution :chi-squared {:degrees-of-freedom DFt})]
     (assoc res
            :chi2 chi2 :df DFt :stat chi2
            :p-value (p-value distr chi2 sides)))))

;; ad/ks tests

(defn- a2-stat
  ^double [^doubles xs d]
  (let [n (alength xs)]
    (reduce - (m/- n) (map (fn [^long idx]
                             (m/* (m// (m/+ idx idx 1.0) n)
                                  (m/+ (m/log (r/cdf d (aget xs idx)))
                                       (m/log (r/ccdf d (aget xs (m/- n idx 1))))))) (range n)))))

(defn ad-test-one-sample
  "Performs the Anderson-Darling (AD) test for goodness-of-fit.

  This test assesses the null hypothesis that a sample `xs` comes from a
  specified theoretical distribution or another empirical distribution. It is
  sensitive to differences in the tails of the distributions.

  Parameters:

  - `xs` (seq of numbers): The sample data to be tested.
  - `distribution-or-ys` (optional):
    - A `fastmath.random` distribution object to test against. If omitted, defaults
      to the standard normal distribution (`fastmath.random/default-normal`).
    - A sequence of numbers (`ys`). In this case, an empirical distribution is
      estimated from `ys` using Kernel Density Estimation (KDE) or an enumerated
      distribution (see `:kernel` option).
  - `opts` (map, optional): Options map:
    - `:sides` (keyword, default `:right`): Specifies the side(s) of the
      A^2 statistic's distribution used for p-value calculation.
      - `:right` (default): Tests if the observed A^2 statistic is significantly
        large (standard approach for AD test, indicating poor fit).
      - `:left`: Tests if the observed A^2 statistic is significantly small.
      - `:two-sided`: Tests if the observed A^2 statistic is extreme in either tail.
    - `:kernel` (keyword, default `:gaussian`): Used only when `distribution-or-ys`
      is a sequence. Specifies the method to estimate the empirical distribution:
        - `:gaussian` (or other KDE kernels): Uses Kernel Density Estimation.
        - `:enumerated`: Creates a discrete empirical distribution from `ys`.
    - `:bandwidth` (double, optional): Bandwidth for KDE (if applicable).

  Returns a map containing:

  - `:A2`: The Anderson-Darling test statistic (A^2).
  - `:stat`: Alias for `:A2`.
  - `:p-value`: The p-value associated with the test statistic and the specified `:sides`.
  - `:n`: Sample size of `xs`.
  - `:mean`: Mean of the sample `xs` (for context).
  - `:stddev`: Standard deviation of the sample `xs` (for context).
  - `:sides`: The alternative hypothesis side used for p-value calculation."
  ([xs] (ad-test-one-sample xs r/default-normal))
  ([xs distribution-or-ys] (ad-test-one-sample xs distribution-or-ys {}))
  ([xs distribution-or-ys {:keys [sides kernel bandwidth] :or {sides :right kernel :gaussian}}]
   (let [d (cond
             (r/distribution? distribution-or-ys) distribution-or-ys
             (= kernel :enumerated) (r/distribution :enumerated-real {:data distribution-or-ys})
             :else (r/distribution :continuous-distribution {:data distribution-or-ys :kde kernel
                                                             :bandwidth bandwidth}))
         axs (m/seq->double-array (sort xs))
         stat (a2-stat axs d)
         n (alength axs)
         distr (r/distribution :anderson-darling {:n n})]
     {:stat stat :A2 stat :mean (mean axs) :stddev (stddev axs) :n n :sides sides
      :p-value (p-value distr stat sides)})))

(defn- ks-jitter-range
  ^double [xs]
  (m/* 0.1 (double (->> (sort xs)
                        (partition 2 1)
                        (map (fn [[^double x ^double y]] (m/- y x)))
                        (filter m/pos?)
                        (reduce m/min)))))

(defn- ks-jitter-seq
  [xs ^double mdiff]
  (let [mdiff- (m/- mdiff)]
    (map m/+ xs (repeatedly #(r/randval (r/drand mdiff- -1.0e-15) (r/drand 1.0e-15 mdiff))))))

(defn- ks-jitter
  ([xs] (ks-jitter-seq xs (ks-jitter-range xs)))
  ([xs ys]
   (let [mdiff (ks-jitter-range (concat xs ys))
         jxs (ks-jitter-seq xs mdiff)
         jys (ks-jitter-seq ys mdiff)]
     [jxs jys])))

(defn ks-test-one-sample
  "Performs the one-sample Kolmogorov-Smirnov (KS) test.

  This test compares the empirical cumulative distribution function (ECDF) of a
  sample `xs` against a specified theoretical distribution or the ECDF of
  another empirical sample. It assesses the null hypothesis that `xs` is drawn
  from the reference distribution.

  Parameters:

  - `xs` (seq of numbers): The sample data to be tested.
  - `distribution-or-ys` (optional):
    - A `fastmath.random` distribution object to test against. If omitted, defaults
      to the standard normal distribution (`fastmath.random/default-normal`).
    - A sequence of numbers (`ys`). In this case, an empirical distribution is
      estimated from `ys` using Kernel Density Estimation (KDE) or an enumerated
      distribution (see `:kernel` option).
  - `opts` (map, optional): Options map:
    - `:sides` (keyword, default `:two-sided`): Specifies the alternative hypothesis
      regarding the difference between the ECDF of `xs` and the reference CDF.
      - `:two-sided` (default): Tests if the ECDF of `xs` is different from the reference CDF.
      - `:right`: Tests if the ECDF of `xs` is significantly *below* the reference CDF (i.e., `xs` tends to have larger values, stochastically greater).
      - `:left`: Tests if the ECDF of `xs` is significantly *above* the reference CDF (i.e., `xs` tends to have smaller values, stochastically smaller).
    - `:kernel` (keyword, default `:gaussian`): Used only when `distribution-or-ys`
      is a sequence. Specifies the method to estimate the empirical distribution:
      - `:gaussian` (or other KDE kernels): Uses Kernel Density Estimation.
      - `:enumerated`: Creates a discrete empirical distribution from `ys`.
    - `:bandwidth` (double, optional): Bandwidth for KDE (if applicable).
    - `:distinct?` (boolean or keyword, default `true`): How to handle duplicate values in `xs`.
      - `true` (default): Removes duplicate values from `xs` before computation.
      - `false`: Uses all values in `xs`, including duplicates.
      - `:jitter`: Adds a small amount of random noise to each value in `xs` to break ties.

  Returns a map containing:

  - `:n`: Sample size of `xs` (after applying `:distinct?`).
  - `:dp`: Maximum positive difference (ECDF(xs) - CDF(ref)).
  - `:dn`: Maximum positive difference (CDF(ref) - ECDF(xs)).
  - `:d`: The KS test statistic (max absolute difference: `max(dp, dn)`).
  - `:stat`: The specific statistic used for p-value calculation, depending on `:sides` (`d`, `dp`, or `dn`).
  - `:p-value`: The p-value associated with the test statistic and the specified `:sides`.
  - `:sides`: The alternative hypothesis side used."
  ([xs] (ks-test-one-sample xs r/default-normal))
  ([xs distribution-or-ys] (ks-test-one-sample xs distribution-or-ys {}))
  ([xs distribution-or-ys {:keys [sides kernel bandwidth distinct?]
                           :or {sides :two-sided kernel :gaussian distinct? true}}]
   (let [d (cond
             (r/distribution? distribution-or-ys) distribution-or-ys
             (= kernel :enumerated) (r/distribution :enumerated-real {:data distribution-or-ys})
             :else (r/distribution :continuous-distribution {:data distribution-or-ys :kde kernel
                                                             :bandwidth bandwidth}))
         xs (cond
              (= :jitter distinct?) (ks-jitter xs)
              distinct? (distinct xs)
              :else xs)
         n (count xs)
         dn (m// (double n))
         idxs (map (fn [^long i] (m/* i dn)) (range (m/inc n)))
         cdfs (map (partial r/cdf d) (sort xs))
         dp (double (reduce m/max (map m/- (rest idxs) cdfs)))
         dn (m/- (double (reduce m/min (map m/- (butlast idxs) cdfs))))
         d (m/max dp dn)]
     {:n n :dp dp :dn dn :d d :sides sides
      :stat (sides-case sides d dp dn) 
      :p-value (sides-case sides
                           (p-value (r/distribution :kolmogorov-smirnov {:n n}) d :right)
                           (p-value (r/distribution :kolmogorov-smirnov+ {:n n}) dp :right)
                           (p-value (r/distribution :kolmogorov-smirnov+ {:n n}) dn :right))})))

(defn- process-ks-diffs
  [vs ^long nx ^long ny]
  (let [vs (vec vs)           ;; concatenated xs and ys
        os (vec (m/order vs)) ;; order of it
        cnt- (m/dec (count vs))
        dx (m// 1.0 nx)
        dy (m// -1.0 ny)]
    (loop [i (long 0)
           d 0.0
           dn 0.0
           dp 0.0]
      (let [id (long (os i))
            nd (m/+ d (if (m/< id nx) dx dy))]
        (if (m/== i cnt-)
          [(m/min nd dn) (m/max nd dp)]
          (let [v1 (double (vs id))
                v2 (double (vs (os (m/inc i))))]
            (if (m/not== v1 v2)
              (recur (m/inc i) nd (m/min nd dn) (m/max nd dp))
              (recur (m/inc i) nd dn dp))))))))

(defn- ks-c-test
  [^double v ^double x ^double y abs?]
  (if abs?
    (m/>= (m/abs (m/- x y)) v)
    (m/>= (m/- x y) v)))

(defn- ks-exact
  ^double [^double d ^long m ^long n {:keys [abs? ties]}]
  (let [not-ties? (not ties)
        lag (double-array (m/inc n))
        md (double m)
        nd (double n)]
    (dotimes [j- n]
      (let [j (m/inc j-)]
        (if (and (ks-c-test d 0.0 (m// j nd) abs?) (or not-ties? (ties j)))
          (Array/aset lag j 1.0)
          (Array/aset lag j (Array/aget lag j-)))))
    (dotimes [i- m]
      (let [i (m/inc i-)
            idmd (m// i md)]
        (when (and (ks-c-test d idmd 0.0 abs?) (or not-ties? (ties i)))
          (Array/aset lag 0 1.0))
        (dotimes [j- n]
          (let [j (m/inc j-)
                i+j (m/+ i j)]
            (if (and (ks-c-test d idmd (m// j nd) abs?) (or not-ties? (ties i+j)))
              (Array/aset lag j 1.0)
              (let [di+j (double i+j)
                    v (m// i di+j)
                    w (m// j di+j)]
                (Array/aset lag j (m/+ (m/* v (Array/aget lag j))
                                       (m/* w (Array/aget lag j-))))))))))
    (Array/aget lag n)))

(defn- ks-correction
  ^double [^double d ^long m ^long n]
  (let [mn (m/* m n)]
    (m// (m/+ 0.5 (m/floor (m/- (m/* d mn) 1.0e-7))) mn)))

(defn- ks-distinct
  [xs ys distinct?]
  (cond
    (= :jitter distinct?) (ks-jitter xs ys)
    (and (not (keyword? distinct?)) distinct?) [(distinct xs) (distinct ys)]
    :else [xs ys]))

(defn- ks-find-ties
  [vs]
  (let [svs (sort vs)
        ties (vec (conj (map (comp m/not-zero? m/-) (rest svs) svs) false))]
    (when (some identity ties) (conj ties true))))

(defn ks-test-two-samples
  "Performs the two-sample Kolmogorov-Smirnov (KS) test.

  This test compares the empirical cumulative distribution functions (ECDFs) of two
  independent samples, `xs` and `ys`, to assess the null hypothesis that they
  are drawn from the same continuous distribution.

  Parameters:

  - `xs` (seq of numbers): The first sample.
  - `ys` (seq of numbers): The second sample.
  - `opts` (map, optional): Options map:
    - `:method` (keyword, optional): Specifies the calculation method for the p-value.
        - `:exact`: Attempts an exact calculation (suitable for small samples, sensitive to ties). Default if `nx * ny < 10000`.
        - `:approximate`: Uses the asymptotic Kolmogorov distribution (suitable for larger samples). Default otherwise.
    - `:sides` (keyword, default `:two-sided`): Specifies the alternative hypothesis.
        - `:two-sided` (default): Tests if the distributions differ (ECDFs are different).
        - `:right`: Tests if `xs` is stochastically greater than `ys` (ECDF(xs) is below ECDF(ys)).
        - `:left`: Tests if `xs` is stochastically smaller than `ys` (ECDF(xs) is above ECDF(ys)).
    - `:distinct?` (keyword or boolean, default `:ties`): How to handle duplicate values (ties).
        - `:ties` (default): Includes all points. Passes information about ties to the `:exact` calculation method. Accuracy depends on the exact method's tie handling.
        - `:jitter`: Adds a small amount of random noise to break ties before comparison. A practical approach if exact tie handling is complex or not required.
        - `true`: Applies `distinct` to `xs` and `ys` separately before combining. May not resolve all ties between the combined samples.
        - `false`: Uses the data as-is, without attempting to handle ties explicitly (may lead to less accurate p-values, especially with the exact method).
    - `:correct?` (boolean, default `true`): Apply continuity correction when using the `:exact` calculation method for a more accurate p-value especially for smaller sample sizes.

  Returns a map containing:

  - `:nx`: Number of observations in `xs` (after `:distinct?` processing if applicable).
  - `:ny`: Number of observations in `ys` (after `:distinct?` processing if applicable).
  - `:n`: Effective sample size used for asymptotic calculation (`nx*ny / (nx+ny)`).
  - `:dp`: Maximum positive difference (ECDF(xs) - ECDF(ys)).
  - `:dn`: Maximum positive difference (ECDF(ys) - ECDF(xs)).
  - `:d`: The KS test statistic (max absolute difference: `max(dp, dn)`).
  - `:stat`: The specific statistic used for p-value calculation (`d`, `dp`, or `dn` for exact; scaled version for approximate).
  - `:KS`: Alias for `:stat`.
  - `:p-value`: The p-value associated with the test statistic and `:sides`.
  - `:sides`: The alternative hypothesis side used.
  - `:method`: The calculation method used (`:exact` or `:approximate`).

  Note on Ties: The KS test is strictly defined for continuous distributions where ties have zero probability.
  The presence of ties in sample data affects the p-value calculation. The `:distinct?` option provides ways to manage this, with `:jitter` being a common pragmatic choice."
  ([xs ys] (ks-test-two-samples xs ys {}))
  ([xs ys {:keys [method sides distinct? correct?]
           :or {sides :two-sided distinct? :ties  correct? true}}]
   (let [[xs ys] (ks-distinct xs ys distinct?)
         nx (count xs)
         ny (count ys)
         method (or method (if (m/< (m/* nx ny) 10000) :exact :approximate))
         vs (concat xs ys)
         ties (when (= distinct? :ties) (ks-find-ties vs))
         [dn dp] (process-ks-diffs vs nx ny)
         dn (m/- (double dn))
         dp (double dp)
         d (m/max dn dp)
         res {:nx nx :ny ny :dp dp :dn dn :d d
              :method method
              :sides sides}]
     (if (= method :exact)
       (let [n (m/+ nx ny)
             stat (sides-case sides d dp dn)
             corrected-stat (if correct? (ks-correction stat nx ny) stat)]
         (assoc res :n n :stat stat :KS stat
                :p-value (sides-case sides
                                     (ks-exact corrected-stat nx ny {:abs? true :ties ties})
                                     (ks-exact corrected-stat nx ny {:abs? false :ties ties})
                                     (ks-exact corrected-stat nx ny {:abs? false :ties ties}))))
       (let [n (m// (m/* nx ny) (double (m/long-add nx ny)))
             stat (m/* (m/sqrt n) (sides-case sides d dp dn))]
         (assoc res :n n :stat stat :KS stat
                :p-value (sides-case sides
                                     (p-value (r/distribution :kolmogorov) stat :right)
                                     (m/exp (m/* -2.0 stat stat))
                                     (m/exp (m/* -2.0 stat stat)))))))))

(defn kruskal-test
  "Performs the Kruskal-Wallis H-test (rank sum test) for independent samples.

  The Kruskal-Wallis test is a non-parametric alternative to one-way ANOVA.
  It determines whether there is a statistically significant difference between the distributions of two or more independent groups. It does not assume normality but requires that distributions have a similar shape for the test to be valid.

  Parameters:

  - `data-groups` (vector of sequences): A collection where each element is a sequence 
    representing a group of observations.
  - a map containing `:sides` key with values of: `:right` (default), `:left` or `:both`

  Returns a map containing:

  - `:stat`: The Kruskal-Wallis H statistic.
  - `:n`: Total number of observations across all groups.
  - `:df`: Degrees of freedom (number of groups - 1).
  - `:k`: Number of groups.
  - `:sides`: Test side
  - `:p-value`: The p-value for the test (null hypothesis: all groups have the same distribution)."
  ([xss] (kruskal-test xss {}))
  ([xss {:keys [sides] :or {sides :right}}] ;; as in R
   (let [k (count xss)
         df (m/dec k)
         xs (flatten xss)
         groups (mapcat (fn [[xs id]] (repeat (count xs) id)) (map vector xss (range)))
         n (count xs)
         r (m/rank1 xs)
         ties (vals (frequencies xs))
         stat (->> (map vector r groups)
                   (group-by second)
                   (vals)
                   (map #(let [ranks (map first %)]
                           (m// (m/sq (sum ranks)) (count ranks))))
                   (sum))
         stat (m// (m/- (m// (m/* 12.0 stat)
                             (m/* n (m/inc n)))
                        (m/* 3.0 (m/inc n)))
                   (m/- 1.0 (m// (sum (map (fn [^long t] (m/- (m/cb t) t)) ties))
                                 (m/- (m/cb n) n))))]
     {:stat stat :n n :df df :k k :sides sides
      :p-value (p-value (r/distribution :chi-squared {:degrees-of-freedom df}) stat sides)})))

;; transformations

(defn- box-cox-scaled
  [nxs sgn ^double lambda gm]
  (let [gm (if (number? gm) (double gm) (geomean nxs))]
    (if (seq sgn)
      (if (m/zero? lambda)
        (map (fn [^double s ^double x] (m/* gm s (m/log (m/inc x)))) sgn nxs)
        (let [fact (m/* lambda (m/pow gm (m/dec lambda)))]
          (map (fn [^double s ^double x] (m// (m/dec (m/* s (m/pow x lambda))) fact)) sgn nxs)))
      (if (m/zero? lambda)
        (map (fn [^double x] (m/* gm (m/log x))) nxs)
        (let [fact (m/* lambda (m/pow gm (m/dec lambda)))]
          (map (fn [^double x] (m// (m/dec (m/pow x lambda)) fact)) nxs))))))

(defn- box-cox-scaled-inv
  [xs ^double lambda {:keys [^double alpha negative? scaled?] :or {alpha 0.0}}]
  (-> (let [gm (double scaled?)]
        (if negative?
          (if (m/zero? lambda)
            (map (fn [^double x] (m/* (m/sgn x) (m/dec (m/exp (m// (m/abs x) gm))))) xs)
            (let [rl (m// lambda)
                  fact (m/* lambda (m/pow gm (m/dec lambda)))]
              (map (fn [^double x]
                     (let [y (m/inc (m/* fact x))]
                       (m/* (m/sgn y) (m/pow (m/abs y) rl)))) xs))            )
          (if (m/zero? lambda)
            (map (fn [^double x] (m/exp (m// x gm))) xs)
            (let [fact (m/* lambda (m/pow gm (m/dec lambda)))
                  rl (m// lambda)]
              (map (fn [^double x] (m/pow (m/inc (m/* fact x)) rl)) xs)))))
      (v/shift (m/- alpha))))

(defn- box-cox-not-scaled
  [nxs sgn ^double lambda]
  (if (seq sgn)
    (if (m/zero? lambda)
      (map (fn [^double s ^double x] (m/* s (m/log (m/inc x)))) sgn nxs)
      (map (fn [^double s ^double x] (m// (m/dec (m/* s (m/pow x lambda))) lambda)) sgn nxs))
    (if (m/zero? lambda)
      (map m/log nxs)
      (map (fn [^double x] (m// (m/dec (m/pow x lambda)) lambda)) nxs))))

(defn- box-cox-not-scaled-inv
  [xs ^double lambda {:keys [^double alpha negative?] :or {alpha 0.0}}]
  (-> (if negative?
        (if (m/zero? lambda)
          (map (fn [^double x] (m/* (m/sgn x) (m/dec (m/exp (m/abs x))))) xs)
          (let [rl (m// lambda)]
            (map (fn [^double x]
                   (let [y (m/inc (m/* lambda x))]
                     (m/* (m/sgn y) (m/pow (m/abs y) rl)))) xs)))
        (if (m/zero? lambda)
          (map m/exp xs)
          (let [rl (m// lambda)]
            (map (fn [^double x] (m/pow (m/inc (m/* lambda x)) rl)) xs))))
      (v/shift (m/- alpha))))

(defn- box-cox-prepare-data
  [xs {:keys [^double alpha negative?]
       :or {alpha 0.0}}]
  (let [d (if (m/zero? alpha) xs (v/shift xs alpha))]
    [(v/abs d)
     (when negative? (map m/sgn d))]))

(defn- box-cox-maximize-ll
  ^double [nxs sgn lambda-range]
  (let [[^double lambda-min ^double lambda-max :as lr] (or lambda-range [-3.0 3.0])
        lxs (v/sum (map m/log nxs))
        n- (m/- (m/* 0.5 (count nxs)))
        target (fn ^double [^double l]
                 (let [res (box-cox-not-scaled nxs sgn l)
                       v (variance res)]
                   (m/+ (m/* n- (m/log v))
                        (m/* (m/dec l) lxs))))]
    (-> (lbfgsb/maximize target {:bounds [lr]
                                 :initial [(m/lerp lambda-min lambda-max 0.51)]})
        (ffirst))))

(defn box-cox-infer-lambda
  "Finds the optimal lambda (λ) parameter for the Box-Cox transformation of a dataset using the Maximum Likelihood Estimation (MLE) method.

  The Box-Cox transformation is a family of power transformations often applied to positive data to make it more closely resemble a normal distribution and stabilize variance. This function estimates the lambda value that maximizes the log-likelihood function of the transformed data, assuming the transformed data is normally distributed.

  Parameters:

  - `xs` (sequence of numbers): The input numerical data sequence.
  - `lambda-range` (vector of two numbers, optional): A sequence `[min-lambda, max-lambda]` defining the closed interval within which the optimal lambda is searched. Defaults to `[-3.0, 3.0]`.
  - `opts` (map, optional): Additional options affecting the data used for the likelihood calculation. These options are passed to the internal data preparation step. Key options include:
    - `:alpha` (double, default 0.0): A constant value added to `xs` before estimating lambda. This is often used when `xs` contains zero or negative values and the standard Box-Cox (which requires positive input) is desired, or to explore transformations around a shifted location.
    - `:negative?` (boolean, default `false`): If `true`, indicates that the likelihood is estimated based on the modified Box-Cox transformation (Bickel and Doksum approach) suitable for negative values. The estimation process will work with the absolute values of the data shifted by `:alpha`.

  Returns the estimated optimal lambda value as a double.

  The inferred lambda value can then be used as the `lambda` parameter for the [[box-cox-transformation]] function to apply the actual transformation to the dataset.

  See also [[box-cox-transformation]], [[yeo-johnson-infer-lambda]], [[yeo-johnson-transformation]]."
  (^double [xs] (box-cox-infer-lambda xs nil))
  (^double [xs lambda-range] (box-cox-infer-lambda xs lambda-range nil))
  (^double [xs lambda-range opts]
   (let [[nxs sgn] (box-cox-prepare-data xs opts)]
     (box-cox-maximize-ll nxs sgn lambda-range))))

(defn box-cox-transformation
  "Applies Box-Cox transformation to a data.

   The Box-Cox transformation is a family of power transformations used to stabilize variance and make data more normally distributed.

  Parameters:

  - `xs` (seq of numbers): The input data.
  - `lambda` (default `0.0`): The power parameter. If `nil` or `[lambda-min, lambda-max]`, `lambda` is inferred using maximum log likelihood.
  - Options map:
    - `alpha` (optional): A shift parameter applied before transformation.
    - `scaled?` (default `false`): Scale by geometric mean or any other number
    - `negative?` (default `false`): Allow negative values
    - `inverse?` (default: `false`): Perform inverse operation, `lambda` can't be inferred.

  Returns transformed data.

  Related: `yeo-johnson-transformation`"
  ([xs] (box-cox-transformation xs nil))
  ([xs lambda] (box-cox-transformation xs lambda nil))
  ([xs lambda {:keys [scaled? inverse?] :as opts}]
   (if inverse?
     (if scaled?
       (box-cox-scaled-inv xs lambda opts)
       (box-cox-not-scaled-inv xs lambda opts))
     (let [[nxs sgn] (box-cox-prepare-data xs opts)
           lambda (if (number? lambda) lambda (box-cox-maximize-ll nxs sgn lambda))]
       (if scaled?
         (box-cox-scaled nxs sgn lambda scaled?)
         (box-cox-not-scaled nxs sgn lambda))))))

(defn- yeo-johnson
  [nxs ^double lambda]
  (let [l2 (m/- 2.0 lambda)]
    (map (fn [^double x]
           (if (m/neg? x)
             (if (m/== lambda 2.0)
               (m/- (m/log (m/- 1.0 x)))
               (m/- (m// (m/dec (m/pow (m/- 1.0 x) l2)) l2)))
             (if (m/zero? lambda)
               (m/log (m/inc x))
               (m// (m/dec (m/pow (m/inc x) lambda)) lambda)))) nxs)))

(defn- yeo-johnson-inv
  [nxs ^double lambda]
  (let [l2 (m/- 2.0 lambda)]
    (map (fn [^double x]
           (if (m/neg? x)
             (if (m/== lambda 2.0)
               (m/- 1.0 (m/exp (m/- x)))
               (m/- 1.0 (m/pow (m/inc (m/* l2 (m/- x))) (m// l2))))
             (if (m/zero? lambda)
               (m/dec (m/exp x))
               (m/dec (m/pow (m/inc (m/* lambda x)) (m// lambda)))))) nxs)))

(defn- yeo-johnson-maximize-ll
  ^double [nxs lambda-range]
  (let [[^double lambda-min ^double lambda-max :as lr] (or lambda-range [-3.0 3.0])
        lxs (v/sum (map (fn [^double x] (m/* (m/signum x) (m/log (m/inc (m/abs x))))) nxs))
        n- (m/- (m/* 0.5 (count nxs)))
        target (fn ^double [^double l]
                 (let [res (yeo-johnson nxs l)
                       v (variance res)]
                   (m/+ (m/* n- (m/log v))
                        (m/* (m/dec l) lxs))))]
    (-> (lbfgsb/maximize target {:bounds [lr]
                                 :initial [(m/lerp lambda-min lambda-max 0.51)]})
        (ffirst))))

(defn yeo-johnson-infer-lambda
  "Find optimal `lambda` parameter for Yeo-Johnson tranformation using maximum log likelihood method."
  (^double [xs] (yeo-johnson-infer-lambda xs nil))
  (^double [xs lambda-range] (yeo-johnson-infer-lambda xs lambda-range nil))
  (^double [xs lambda-range {:keys [^double alpha] :or {alpha 0.0}}]
   (let [nxs (if (m/zero? alpha) xs (v/shift xs alpha))]
     (yeo-johnson-maximize-ll nxs lambda-range))))

(defn yeo-johnson-transformation
  "Applies the Yeo-Johnson transformation to a dataset.

  This transformation is used to stabilize variance and make data more normally distributed. It extends the Box-Cox transformation to allow for zero and negative values.

  Parameters:

  - `xs`: The input dataset.
  - `lambda` (default: 0.0): The power parameter controlling the transformation. If `lambda` is `nil` or a range `[lambda-min, lambda-max]` it will be inferred using maximum log-likelihood method.
  - Options map:
    - `:alpha` (optional): A shift parameter applied before transformation.
    - `:inverse?` (optional): Perform inverse operation, `lambda` should be provided (can't be inferred). 

  Returns:

  - A transformed sequence of numbers.

  Related: `box-cox-tranformation`"
  ([xs] (yeo-johnson-transformation xs nil))
  ([xs lambda] (yeo-johnson-transformation xs lambda nil))
  ([xs lambda {:keys [^double alpha inverse?] :or {alpha 0.0}}]
   (if inverse?
     (let [nxs (yeo-johnson-inv xs lambda)]
       (if (m/zero? alpha) nxs (v/shift nxs (m/- alpha))))
     (let [nxs (if (m/zero? alpha) xs (v/shift xs alpha))
           lambda (if (number? lambda) lambda (yeo-johnson-maximize-ll nxs lambda))]
       (yeo-johnson nxs lambda)))))

;;

(defn power-transformation
  "Applies a power transformation to a data."
  {:deprecated "Use `(box-cox-transformation xs lambda {:scaled true})"}
  ([xs] (power-transformation xs 0.0))
  ([xs ^double lambda]
   (box-cox-transformation xs lambda {:scaled? true}))
  ([xs ^double lambda ^double alpha]
   (box-cox-transformation xs lambda {:scaled? true :alpha alpha})))

(defn modified-power-transformation
  "Applies a modified power transformation (Bickel and Doksum) to a data."
  {:deprecated "Use `(box-cox-transformation xs lambda {:negative? true})"}
  ([xs] (modified-power-transformation xs 0.0))
  ([xs ^double lambda]
   (box-cox-transformation xs lambda {:negative? true}))
  ([xs ^double lambda ^double alpha]
   (box-cox-transformation xs lambda {:negative? true :alpha alpha})))

