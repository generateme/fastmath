(ns stats
  (:require [fastmath.stats :as stats]
            [fastmath.stats.bootstrap :as boot]
            [fastmath.dev.ggplot :as gg]
            [fastmath.dev.codox :as codox]
            [fastmath.dev.clay :as utls]
            [fastmath.dev.dataset :as ds]
            [fastmath.random :as r]
            [fastmath.core :as m]
            [scicloj.kindly.v4.kind :as kind]
            [fastmath.vector :as v]))

;; # Statistics {.unnumbered}

;; Statistical functions

;; ## Datasets

;; To illustrate statistical functions two dataset will be used: `mtcars` and `iris`. ;; A `fastmath.dev.dataset` as `ds` will be used to access data.

;; Select the `:mpg` column

(ds/mtcars :mpg)

;; Select the `:mpg` column with a row predicate

(ds/mtcars :mpg (fn [row] (and (< (row :mpg) 20.0) (zero? (row :am)))))

;; Group by the `:am` column and select `:mpg`

(ds/by ds/mtcars :am :mpg)

;; ### mtcars

;; 11 attributes of 32 different cars comparison

(kind/table (ds/mtcars) {:use-datatables true})

;; ### iris

;; Sepal and Petal iris species comparison 

(kind/table (ds/iris) {:use-datatables true})

;; ## Basic Descriptive Statistics

;; General summary statistics for describing the center and range of data.

;; ::: {.callout-tip title="Defined functions"}
;; * `minimum`, `maximum`
;; * `sum`
;; * `mean`
;; * `median`, `median-3`
;; * `mode`, `modes`
;; * `wmode`, `wmodes`
;; * `geomean`, `harmean`, `powmean`
;; * `stats-map`
;; :::

;; ## Measures of Dispersion/Deviation

;; Statistics describing the spread or variability of data.

;; ::: {.callout-tip title="Defined functions"}
;; * `variance`, `population-variance`
;; * `stddev`, `population-stddev`
;; * `wvariance`, `population-wvariance`
;; * `wstddev`, `population-wstddev`
;; * `pooled-variance`, `pooled-stddev`
;; * `variation`
;; * `mean-absolute-deviation`
;; * `median-absolute-deviation`, `mad`
;; * `pooled-mad`
;; * `sem` ;; Standard Error of Mean
;; * `iqr` ;; Interquartile Range
;; :::

;; ## Quantiles and Percentiles

;; Statistics related to points dividing the distribution of data.

;; ::: {.callout-tip title="Defined functions"}

;; * `percentile`, `percentiles`
;; * `quantile`, `quantiles`
;; * `wquantile`, `wquantiles`
;; * `wmedian`
;; * `expectile`
;; :::

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

;; ## Moments and Shape

;; Statistics describing the shape of the distribution (symmetry, peakedness).

;; ::: {.callout-tip title="Defined functions"}
;; * `moment`
;; * `skewness`
;; * `kurtosis`
;; * `l-moment`, `l-variation`
;; :::

;; ## Intervals and Extents

;; Functions defining ranges or intervals based on the data.

;; ::: {.callout-tip title="Defined functions"}
;; * `span`, `extent`
;; * `stddev-extent`, `mad-extent`, `sem-extent`
;; * `percentile-extent`, `quantile-extent`
;; * `pi`, `pi-extent` ;; Percentile Interval
;; * `hpdi-extent` ;; Highest Posterior Density Interval
;; * `adjacent-values`
;; * `inner-fence-extent`, `outer-fence-extent`
;; * `ci` ;; Confidence Interval (Student's t)
;; :::

;; ## Outlier Detection

;; Functions to identify or remove potential outliers.

;; ::: {.callout-tip title="Defined functions"}
;; * `outliers`
;; * `remove-outliers`
;; :::

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
