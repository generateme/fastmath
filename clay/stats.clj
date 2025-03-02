^:kindly/hide-code
(ns stats
  (:require [fastmath.core :as m]
            [fastmath.stats :as stats]
            [fastmath.vector :as v]
            [fastmath.random :as r]
            [clojure.string :as str]
            [clojure.java.io :as io]

            [scicloj.kindly.v4.kind :as kind]
            [fastmath.dev.ggplot :as gg]
            [fastmath.dev.clay :as utls]
            [fastmath.dev.codox :as codox]))

;; # Statistical Functions {.unnumbered}

;; The `fastmath.stats` namespace provides a comprehensive set of statistical functions for data analysis. This tutorial will cover the main components, demonstrate real-world applications, and explore most of the API functions.

;; First, let's require the necessary namespaces:

(require '[fastmath.stats :as stats]
         '[fastmath.core :as m]
         '[fastmath.random :as r]
         '[clojure.string :as str]
         '[clojure.java.io :as io])

;; ## Loading Real-World Datasets

;; For this tutorial, we'll load datasets from multiple sources. We'll use:
;; 
;; 1. **Iris Dataset**: Comes with the fastmath library (local resource)
;; 2. **Motor Trend Car Road Tests (mtcars)**: Comes with the fastmath library (local resource)
;; 3. **Rock Dataset**: Loaded from a public URL (measurements of 48 rock samples)
;; 4. **Diamonds Dataset**: Loaded from a public URL (information about diamond prices and attributes)

;; ### CSV Processing Functions

;; First, let's add the Charred library dependency for proper CSV parsing
(require '[charred.api :as csv])

(defn parse-number 
  "Parse a string to a number, if possible"
  [s]
  (if (re-matches #"^-?\d+(\.\d+)?$" s)
    (Double/parseDouble s)
    s))

(defn process-csv-data 
  "Process CSV data with header into maps using Charred for proper CSV parsing"
  [content]
  (let [csv-data (cond
                   (string? content) (csv/read-csv content)
                   (sequential? content) (csv/read-csv (str/join "\n" content))
                   :else (throw (Exception. "Invalid CSV content format")))
        headers (map keyword (first csv-data))
        rows (rest csv-data)]
    (map (fn [row]
           (zipmap headers 
                   (map parse-number row)))
         rows)))

(defn load-csv-from-url
  "Load CSV data from a URL as a string (to be parsed with Charred)"
  [url]
  (let [connection (-> url java.net.URL. .openConnection)
        _ (.setRequestProperty connection "User-Agent" "Mozilla/5.0")
        response (slurp (.getInputStream connection))]
    response))

;; ### Local Datasets

;; Load and process the Iris dataset
(def iris-raw 
  (with-open [reader (io/reader (io/resource "iris.csv"))]
    (slurp reader)))

(def iris-data (process-csv-data iris-raw))

;; Extract just the numerical data from iris for some analyses
(def iris-measurements 
  (map (fn [row] 
         [(row :sepal_length) 
          (row :sepal_width) 
          (row :petal_length) 
          (row :petal_width)])
       iris-data))

;; Create a vector for each measurement type
(def sepal-length (map :sepal_length iris-data))
(def sepal-width (map :sepal_width iris-data))
(def petal-length (map :petal_length iris-data))
(def petal-width (map :petal_width iris-data))

;; Group by species
(def iris-by-species (group-by :species iris-data))

;; Extract data for each species
(def setosa-data (map #(vector (:sepal_length %) (:sepal_width %) (:petal_length %) (:petal_width %)) 
                     (iris-by-species "setosa")))
(def versicolor-data (map #(vector (:sepal_length %) (:sepal_width %) (:petal_length %) (:petal_width %)) 
                         (iris-by-species "versicolor")))
(def virginica-data (map #(vector (:sepal_length %) (:sepal_width %) (:petal_length %) (:petal_width %)) 
                        (iris-by-species "virginica")))

;; Load and process the mtcars dataset
(def mtcars-raw 
  (with-open [reader (io/reader (io/resource "mtcars.csv"))]
    (slurp reader)))

;; Process mtcars data which has quoted column names (using Charred)
(def mtcars-parsed (csv/read-csv mtcars-raw))
(def mtcars-headers (map keyword (first mtcars-parsed)))

(def mtcars-data 
  (map (fn [row]
         (let [car-name (first row)
               values (map parse-number (rest row))]
           (assoc (zipmap (rest mtcars-headers) values) :name car-name)))
       (rest mtcars-parsed)))

;; Extract some useful columns from mtcars
(def mpg (map :mpg mtcars-data))
(def hp (map :hp mtcars-data))
(def wt (map :wt mtcars-data))
(def qsec (map :qsec mtcars-data))

;; Extract MPG data by transmission type for later analysis
(def auto-trans-mpg (map :mpg (filter #(zero? (:am %)) mtcars-data)))
(def manual-trans-mpg (map :mpg (filter #(= 1.0 (:am %)) mtcars-data)))

;; ### Remote Datasets (loaded at runtime)

;; Load and process the Rock dataset from a public URL
(def rock-raw 
  (try
    (load-csv-from-url "https://vincentarelbundock.github.io/Rdatasets/csv/datasets/rock.csv")
    (catch Exception e
      (println "Failed to load Rock dataset, using placeholder data")
      "X,area,peri,shape,perm\n1,4990,2791.90,0.0903296,6.3\n2,7002,3892.60,0.1486220,6.3")))

;; Process the Rock dataset into maps (using Charred)
(def rock-data (process-csv-data rock-raw))

;; Extract columns for analysis
(def rock-area (map :area rock-data))
(def rock-peri (map :peri rock-data))
(def rock-shape (map :shape rock-data))
(def rock-perm (map :perm rock-data))

;; Load and process the Diamonds dataset from a public URL (using just a sample of 200 rows)
(def diamonds-url "https://raw.githubusercontent.com/tidyverse/ggplot2/master/data-raw/diamonds.csv")
(def diamonds-raw 
  (try
    (load-csv-from-url diamonds-url)
    (catch Exception e
      (println "Failed to load Diamonds dataset, using placeholder data")
      "carat,cut,color,clarity,depth,table,price,x,y,z\n0.23,\"Ideal\",\"E\",\"SI2\",61.5,55,326,3.95,3.98,2.43\n0.21,\"Premium\",\"E\",\"SI1\",59.8,61,326,3.89,3.84,2.31")))

;; Process the Diamonds dataset into maps
(def diamonds-data (process-csv-data diamonds-raw))

;; Extract numerical columns for analysis
(def diamond-carat (map :carat diamonds-data))
(def diamond-depth (map :depth diamonds-data))
(def diamond-table (map :table diamonds-data))
(def diamond-price (map :price diamonds-data))
(def diamond-x (map :x diamonds-data)) ; length in mm
(def diamond-y (map :y diamonds-data)) ; width in mm
(def diamond-z (map :z diamonds-data)) ; depth in mm

;; ## Comprehensive Descriptive Statistics

;; ### Basic Summary Statistics

;; Let's start by exploring some key statistics from the Iris dataset.
;; We'll examine sepal length to understand these fundamental measures:

(utls/examples-note
 ;; Count - Number of observations
 (count sepal-length)

 ;; Sum - Total of all values
 (stats/sum sepal-length)
 
 ;; Minimum - Smallest value
 (stats/minimum sepal-length)
 
 ;; Maximum - Largest value
 (stats/maximum sepal-length)
 
 ;; Range - Difference between max and min
 (- (stats/maximum sepal-length) (stats/minimum sepal-length)))

;; The [minimum](https://en.wikipedia.org/wiki/Sample_maximum_and_minimum) and [maximum](https://en.wikipedia.org/wiki/Sample_maximum_and_minimum) values define the range of the data, which is useful for understanding the spread.

;; ### Central Tendency Measures

;; Central tendency measures help identify the "typical" value in a dataset:

(utls/examples-note
 ;; Mean - arithmetic average
 (stats/mean sepal-length)
 
 ;; Median - middle value when sorted
 (stats/median sepal-length)
 
 ;; Mode - most frequent value(s)
 (stats/mode sepal-length)
 
 ;; Modes - all most frequent values
 (stats/modes sepal-length))

;; The [mean](https://en.wikipedia.org/wiki/Arithmetic_mean) is the average value but can be influenced by outliers. The [median](https://en.wikipedia.org/wiki/Median) is more robust to outliers. The [mode](https://en.wikipedia.org/wiki/Mode_(statistics)) identifies the most common value(s).

;; There are also different types of means for special cases:

(utls/examples-note
 ;; Geometric mean - for data representing multiplicative processes
 ;; Note: All values must be positive
 (stats/geomean (map (partial max 0.1) petal-width))
 
 ;; Harmonic mean - often used for rates and ratios
 (stats/harmean (map (partial max 0.1) petal-width))
 
 ;; Power mean with various powers
 (stats/powmean petal-width 2)  ;; Quadratic mean (or root-mean-square)
 (stats/powmean petal-width 0)  ;; Geometric mean
 (stats/powmean petal-width -1)) ;; Harmonic mean

;; The [geometric mean](https://en.wikipedia.org/wiki/Geometric_mean) is the nth root of the product of n values, useful for growth rates. The [harmonic mean](https://en.wikipedia.org/wiki/Harmonic_mean) is the reciprocal of the arithmetic mean of the reciprocals, useful for averaging rates.

;; ### Weighted Statistics

;; Sometimes observations have different importance or frequency. 
;; Weighted statistics account for these differences:

;; Let's create a weighted dataset to demonstrate
(def weights [5 2 1 3 4 1 3 2 5 4])
(def values [10 15 20 25 30 35 40 45 50 55])

(utls/examples-note
 ;; Weighted mean
 (stats/mean values weights)
 
 ;; Standard mean (for comparison)
 (stats/mean values)
 
 ;; Weighted median
 (stats/wmedian values weights)
 
 ;; Weighted quantiles (25%, 50%, 75%)
 (stats/wquantiles values weights [0.25 0.5 0.75])
 
 ;; Weighted variance
 (stats/wvariance values weights))

;; Weighted statistics are useful when some data points should have more influence than others, such as when working with measurements of different reliability or importance.

;; ### Dispersion and Variability Measures

;; Dispersion measures quantify how spread out the data is:

(utls/examples-note
 ;; Variance - average squared deviation from the mean
 (stats/variance sepal-length)
 
 ;; Standard deviation - square root of variance
 (stats/stddev sepal-length)
 
 ;; Population variance and standard deviation 
 ;; (when data represents the entire population rather than a sample)
 (stats/population-variance sepal-length)
 (stats/population-stddev sepal-length)
 
 ;; Coefficient of variation - standardized measure of dispersion
 (stats/variation sepal-length)
 
 ;; Median Absolute Deviation (MAD) - robust measure of variability
 (stats/median-absolute-deviation sepal-length)
 
 ;; Mean Absolute Deviation - another measure of variability
 (stats/mean-absolute-deviation sepal-length)
 
 ;; Interquartile Range (IQR) - difference between 75th and 25th percentiles
 (stats/iqr sepal-length))

;; [Standard deviation](https://en.wikipedia.org/wiki/Standard_deviation) and [variance](https://en.wikipedia.org/wiki/Variance) quantify the average deviation from the mean. [Median Absolute Deviation](https://en.wikipedia.org/wiki/Median_absolute_deviation) (MAD) is more robust to outliers. The [Interquartile Range](https://en.wikipedia.org/wiki/Interquartile_range) covers the middle 50% of the data.

;; ### Shape Measures

;; Shape measures describe the distribution's form:

(utls/examples-note
 ;; Skewness - asymmetry of the distribution
 (stats/skewness sepal-length)
 (stats/skewness sepal-length :g1)  ;; Pearson's moment coefficient of skewness
 (stats/skewness sepal-length :b1)  ;; Adjusted Fisher-Pearson standardized moment coefficient
 
 ;; Kurtosis - "tailedness" of the distribution
 (stats/kurtosis sepal-length)
 (stats/kurtosis sepal-length :g2)  ;; Excess kurtosis
 (stats/kurtosis sepal-length :moors)) ;; Moors' kurtosis (robust measure)

;; [Skewness](https://en.wikipedia.org/wiki/Skewness) measures asymmetry. Positive skewness indicates a right tail (higher values stretch out), while negative skewness indicates a left tail. A normal distribution has a skewness of 0.

;; [Kurtosis](https://en.wikipedia.org/wiki/Kurtosis) measures the "heaviness" of the tails. Higher kurtosis indicates more outliers. A normal distribution has a kurtosis of 3 (or excess kurtosis of 0).

;; ### Percentiles and Quantiles

;; Percentiles and quantiles show the distribution by dividing data into equal parts:

(utls/examples-note
 ;; Single percentile - value below which X% of observations fall
 (stats/percentile sepal-length 25)  ;; 25th percentile (Q1)
 
 ;; Multiple percentiles
 (stats/percentiles sepal-length [25 50 75 90 95 99])
 
 ;; Quantiles - same as percentiles but in range [0,1]
 (stats/quantile sepal-length 0.25)  ;; 0.25 quantile (Q1)
 (stats/quantile sepal-length 0.5)   ;; 0.5 quantile (median)
 
 ;; Multiple quantiles
 (stats/quantiles sepal-length [0.25 0.5 0.75 0.9 0.95 0.99])
 
 ;; Different estimation strategies affect how quantiles are calculated
 ;; This matters most when there aren't exact data points at the quantile
 (stats/quantile sepal-length 0.5 :legacy)
 (stats/quantile sepal-length 0.5 :r7))  ;; R-7 method (commonly used in R)

;; [Percentiles](https://en.wikipedia.org/wiki/Percentile) and [quantiles](https://en.wikipedia.org/wiki/Quantile) divide the data into equal portions. They help understand the distribution without assuming a specific shape.

;; ### Extents and Ranges

;; Extents provide ranges around central values:

(utls/examples-note
 ;; Standard deviation extent (mean ± stddev) - contains ~68% of data if normal
 (stats/stddev-extent sepal-length)
 
 ;; MAD extent (median ± MAD) - robust alternative to stddev extent
 (stats/mad-extent sepal-length)
 
 ;; Standard error of mean extent - useful for confidence intervals
 (stats/sem-extent sepal-length)
 
 ;; Percentile extent - based on percentiles (default: 25%, 75%, and median)
 (stats/percentile-extent sepal-length)
 
 ;; Inner fence extent - based on IQR, used for outlier detection
 (stats/inner-fence-extent sepal-length)
 
 ;; Adjacent values - non-outlier extremes
 (stats/adjacent-values sepal-length))

;; Extents provide a quick way to understand the typical range of values in a dataset. Inner fences are used in box plots to identify potential outliers.

;; ### Complete Statistical Summary

;; The `stats-map` function provides a comprehensive statistical summary:

(def summary (stats/stats-map sepal-length))

;; Let's look at some key parts of the summary
(select-keys summary [:Mean :Median :Mode :SD :Skewness :Kurtosis :IQR])

;; The full summary includes over 20 different statistics, making it a powerful way to understand your data at a glance:
summary

;; ## Data Transformation and Preprocessing

;; ### Standardization and Normalization

;; Standardization and normalization are essential for many machine learning algorithms:

(utls/examples-note
 ;; Standardize - transform to mean=0, stddev=1 (z-scores)
 (take 5 (stats/standardize sepal-length))
 
 ;; Robust standardization - using median and MAD instead of mean and stddev
 (take 5 (stats/robust-standardize sepal-length))
 
 ;; Rescale to a specific range (default: [0,1])
 (take 5 (stats/rescale sepal-length))
 
 ;; Rescale to [-1,1] range
 (take 5 (stats/rescale sepal-length -1 1)))

;; [Standardization](https://en.wikipedia.org/wiki/Standard_score) makes data comparable across features with different scales. [Rescaling](https://en.wikipedia.org/wiki/Feature_scaling) is useful for algorithms that expect inputs in a specific range.

;; ### Data Cleaning and Processing

;; Functions for preparing data for analysis:

(utls/examples-note
 ;; Remove the mean
 (take 5 (stats/demean sepal-length))
 
 ;; Trim outliers (default: remove values outside middle 60%)
 (stats/trim sepal-length)
 
 ;; Winsorize data (replace outliers with boundary values)
 (stats/winsor sepal-length)
 
 ;; Find outliers
 (stats/outliers (concat sepal-length [10 10.5])) ;; add outliers for demonstration
 
 ;; Remove outliers
 (count (stats/remove-outliers (concat sepal-length [10 10.5]))))

;; [Trimming](https://en.wikipedia.org/wiki/Trimmed_estimator) and [winsorizing](https://en.wikipedia.org/wiki/Winsorizing) help deal with outliers. The former removes them, while the latter replaces them with boundary values.

;; ### Statistical Transformations

;; Transformations can make non-normal data more normal and stabilize variance:

(utls/examples-note
 ;; Power transformation (Box-Cox transformation)
 ;; Lambda parameter controls the transformation
 ;; Lambda = 0 is log transform, Lambda = 1 is no transformation
 ;; Note: Data must be positive
 (take 5 (stats/power-transformation petal-width 0.5))  ;; Square root-like transformation
 
 ;; Modified power transformation - allows negative values
 (take 5 (stats/modified-power-transformation petal-width 0.5))
 
 ;; Yeo-Johnson transformation - an alternative that handles negative values
 (take 5 (stats/yeo-johnson-transformation petal-length 0.5)))

;; [Power transformations](https://en.wikipedia.org/wiki/Power_transform) like Box-Cox help normalize skewed data. [Yeo-Johnson](https://en.wikipedia.org/wiki/Power_transform#Yeo%E2%80%93Johnson_transformation) is a generalization that handles negative values.

;; ## Correlation and Covariance

;; ### Understanding Relationships Between Variables

;; Correlation and covariance measure relationships between variables:

(utls/examples-note
 ;; Covariance - measures how two variables change together
 (stats/covariance sepal-length petal-length)
 
 ;; Pearson correlation - measures linear correlation
 (stats/pearson-correlation sepal-length petal-length)
 (stats/correlation sepal-length petal-length) ;; Same as Pearson
 
 ;; Spearman rank correlation - measures monotonic relationships
 (stats/spearman-correlation sepal-length petal-length)
 
 ;; Kendall rank correlation - another rank-based measure
 (stats/kendall-correlation sepal-length petal-length))

;; [Pearson correlation](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) measures linear relationships, while [Spearman](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient) and [Kendall](https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient) measure rank-based relationships, making them more robust to outliers and non-linear patterns.

;; ### Correlation and Covariance Matrices

;; For multivariate data, we can compute full correlation and covariance matrices:

;; Let's use the iris measurements data
(def iris-measurements-transposed (apply map vector iris-measurements))

;; Correlation matrix shows correlations between all pairs of variables
(stats/correlation-matrix iris-measurements-transposed)

;; Covariance matrix shows covariances between all pairs of variables
(stats/covariance-matrix iris-measurements-transposed)

;; We can specify which correlation method to use
(stats/correlation-matrix iris-measurements-transposed :spearman)

;; [Correlation matrices](https://en.wikipedia.org/wiki/Correlation_and_dependence#Correlation_matrices) show all pairwise correlations at once, useful for multivariate analysis. [Covariance matrices](https://en.wikipedia.org/wiki/Covariance_matrix) are fundamental in many statistical methods, including principal component analysis.

;; ## Regression Error Metrics

;; ### Evaluating Prediction Accuracy

;; These metrics quantify the error between predicted and actual values:

;; Let's create some predicted values for demonstration
(def actual mpg)
(def predicted (map #(+ % (- (rand 3) 1.5)) actual))  ;; actual with some noise

(utls/examples-note
 ;; Mean Error - average of errors (can have cancellation)
 (stats/me [actual predicted])
 
 ;; Mean Absolute Error - average of absolute errors
 (stats/mae [actual predicted])
 
 ;; Mean Absolute Percentage Error - percentage errors
 (stats/mape [actual predicted])
 
 ;; Residual Sum of Squares - sum of squared errors
 (stats/rss [actual predicted])
 
 ;; R-squared - coefficient of determination
 (stats/r2 [actual predicted])
 
 ;; Mean Squared Error - average of squared errors
 (stats/mse [actual predicted])
 
 ;; Root Mean Squared Error - square root of MSE
 (stats/rmse [actual predicted]))

;; These metrics help evaluate the performance of regression models. [RMSE](https://en.wikipedia.org/wiki/Root-mean-square_deviation) penalizes large errors more than MAE. [R²](https://en.wikipedia.org/wiki/Coefficient_of_determination) indicates the proportion of variance explained.

;; ### Distance Metrics

;; Distance metrics quantify the difference between vectors:

(utls/examples-note
 ;; L0 - count of different elements
 (stats/L0 [1 2 3 4] [1 2 4 4])
 
 ;; L1 - Manhattan distance (sum of absolute differences)
 (stats/L1 [1 2 3 4] [2 3 4 5])
 
 ;; L2 - Euclidean distance (square root of sum of squared differences)
 (stats/L2 [1 2 3 4] [2 3 4 5])
 
 ;; L2 squared - Sum of squared differences
 (stats/L2sq [1 2 3 4] [2 3 4 5])
 
 ;; L-infinity - Maximum absolute difference
 (stats/LInf [1 2 3 4] [2 3 4 8])
 
 ;; Peak Signal-to-Noise Ratio - used in image quality assessment
 (stats/psnr [1 2 3 4] [2 3 4 5]))

;; [Distance metrics](https://en.wikipedia.org/wiki/Metric_(mathematics)) quantify the difference between vectors in different ways. Each has specific properties making it suitable for different applications.

;; ## Histograms and Distributions

;; ### Creating and Analyzing Histograms

;; Histograms divide data into bins to show the distribution:

;; Create a histogram of sepal lengths
(def sepal-length-hist (stats/histogram sepal-length))

;; Basic information about the histogram
(select-keys sepal-length-hist [:size :min :max :step :samples])

;; Bin information - range and count for each bin
(:bins sepal-length-hist)

;; We can control the number of bins or the binning method
(def sepal-length-hist-10bins (stats/histogram sepal-length 10))
(:bins sepal-length-hist-10bins)

(def sepal-length-hist-sturges (stats/histogram sepal-length :sturges))
(:bins sepal-length-hist-sturges)

;; Different bin estimation methods for histograms
(stats/estimate-bins sepal-length :sqrt)      ;; Square root rule
(stats/estimate-bins sepal-length :sturges)   ;; Sturges' formula
(stats/estimate-bins sepal-length :rice)      ;; Rice rule
(stats/estimate-bins sepal-length :doane)     ;; Doane's formula
(stats/estimate-bins sepal-length :scott)     ;; Scott's normal reference rule

;; [Histograms](https://en.wikipedia.org/wiki/Histogram) visualize data distribution by counting values in each bin. The choice of bin count affects how the distribution appears.

;; ### Comparing Distributions

;; There are several ways to compare distributions:

;; Calculate the similarity between two distributions
(def dist1 petal-length)
(def dist2 petal-width)

(utls/examples-note
 ;; Dissimilarity metrics
 ;; Kullback-Leibler divergence
 (stats/dissimilarity :kullback-leibler dist1 dist2)
 
 ;; Jensen-Shannon divergence
 (stats/dissimilarity :jensen-shannon dist1 dist2)
 
 ;; Hellinger distance
 (stats/dissimilarity :hellinger dist1 dist2)
 
 ;; Euclidean-based distance
 (stats/dissimilarity :euclidean dist1 dist2)
 
 ;; Similarity metrics (higher means more similar)
 ;; Intersection similarity
 (stats/similarity :intersection dist1 dist2)
 
 ;; Cosine similarity
 (stats/similarity :cosine dist1 dist2))

;; [Distribution comparisons](https://en.wikipedia.org/wiki/Statistical_distance) are useful in many domains, including machine learning, information retrieval, and statistical testing.

;; ## Bootstrap Methods

;; ### Resampling for Robust Statistics

;; Bootstrapping creates resamples from data to estimate statistics and their uncertainty:

;; Generate bootstrap resamples of sepal length
(def bootstrap-samples 
  (repeatedly 100 
              #(stats/mean (r/->seq (r/distribution :enumerated-real 
                                                    {:data sepal-length}) 
                                    (count sepal-length)))))

;; Calculate bootstrap confidence intervals (BCa - bias-corrected and accelerated)
(stats/percentile-bca-extent bootstrap-samples)

;; Or use simpler percentile-based intervals
(stats/percentile-extent bootstrap-samples 2.5)

;; ### Advanced Bootstrap Examples with Real-World Data

;; Let's apply bootstrapping to our public datasets to solve practical problems.

;; #### Example 1: Confidence Interval for Diamond Price/Carat Ratio

;; The price-per-carat ratio is an important metric in diamond valuation.
;; Let's use bootstrapping to estimate the confidence interval for the mean price/carat.

;; Calculate price per carat for each diamond
(def price-per-carat (map #(/ %1 %2) diamond-price diamond-carat))

;; Direct calculation of the mean
(stats/mean price-per-carat)

;; Bootstrap to estimate confidence interval for the mean price/carat
(def ppc-bootstrap-samples
  (repeatedly 100
              #(stats/mean (r/->seq (r/distribution :enumerated-real
                                                    {:data price-per-carat})
                                    (count price-per-carat)))))

;; Calculate 95% confidence interval
(stats/percentile-extent ppc-bootstrap-samples 2.5)

;; #### Example 2: Robust Correlation Estimation for Rock Properties

;; We can use bootstrapping to get a more robust estimate of the correlation
;; between rock area and permeability, along with confidence intervals.

;; Function to calculate correlation from resampled data
(defn bootstrap-correlation [xs ys]
  (let [n (count xs)
        indices (repeatedly n #(r/irand n))
        x-sample (map #(nth xs %) indices)
        y-sample (map #(nth ys %) indices)]
    (stats/correlation x-sample y-sample)))

;; Generate bootstrap samples of correlation
(def area-perm-corr-samples
  (repeatedly 100 #(bootstrap-correlation rock-area rock-perm)))

;; Calculate the mean and confidence interval
(stats/mean area-perm-corr-samples)
(stats/percentile-extent area-perm-corr-samples 2.5)

;; #### Example 3: Comparing Medians with Bootstrapped Hypothesis Testing

;; We'll compare the median MPG of automatic vs manual transmission cars
;; using a bootstrapped hypothesis test approach.

;; Calculate the observed difference in medians
(def observed-diff (- (stats/median manual-trans-mpg) (stats/median auto-trans-mpg)))

;; Function to simulate under the null hypothesis (no difference)
(defn simulate-null-hypothesis [xs ys]
  (let [pooled (concat xs ys)
        n1 (count xs)
        n2 (count ys)
        shuffled (shuffle pooled)
        sim-xs (take n1 shuffled)
        sim-ys (drop n1 shuffled)]
    (- (stats/median sim-ys) (stats/median sim-xs))))


;; Generate distribution under null hypothesis
(def null-distribution
  (repeatedly 100 #(simulate-null-hypothesis auto-trans-mpg manual-trans-mpg)))

;; Calculate the p-value (two-tailed test)
(def p-value
  (/ (count (filter #(>= (m/abs %) (m/abs observed-diff)) null-distribution))
     100))

;; Print results
observed-diff
p-value

;; [Bootstrapping](https://en.wikipedia.org/wiki/Bootstrapping_(statistics)) is a powerful technique for estimating the sampling distribution of a statistic, especially when the underlying distribution is unknown or complex. It allows us to make robust inferences without relying on distributional assumptions.

;; ## Hypothesis Testing

;; ### Statistical Tests

;; Hypothesis tests determine if observed data supports a statistical hypothesis:

;; t-tests for comparing means
(utls/examples-note
 ;; One-sample t-test - compare a sample mean to a hypothesized value
 (stats/t-test-one-sample sepal-length {:mu 5.5})
 
 ;; Two-sample t-test - compare means of two samples
 (stats/t-test-two-samples 
  (map first setosa-data)
  (map first versicolor-data))
 
 ;; Paired t-test - compare matched pairs
 (stats/t-test-two-samples 
  (map first setosa-data) 
  (map second setosa-data) 
  {:paired? true}))

;; z-tests (for when population standard deviation is known)
(utls/examples-note
 ;; One-sample z-test
 (stats/z-test-one-sample sepal-length {:mu 5.8})
 
 ;; Two-sample z-test
 (stats/z-test-two-samples
  (map first setosa-data)
  (map first versicolor-data)))

;; F-test for comparing variances
(stats/f-test (map first setosa-data)
             (map first versicolor-data))

;; ANOVA for comparing multiple groups
(stats/one-way-anova-test 
 [(map first setosa-data)
  (map first versicolor-data)
  (map first virginica-data)])

;; Normality tests
(utls/examples-note
 ;; Test based on skewness and kurtosis
 (stats/normality-test sepal-length)
 
 ;; Jarque-Bera test
 (stats/jarque-bera-test sepal-length))

;; Binomial test (for binary data)
(stats/binomial-test [true false true true false true true])

;; [Hypothesis testing](https://en.wikipedia.org/wiki/Statistical_hypothesis_testing) is foundational in statistical inference. The p-value from these tests indicates the probability of observing the data (or more extreme) if the null hypothesis is true.

;; ### Effect Size Measures

;; Effect size quantifies the magnitude of an effect or relationship:

(utls/examples-note
 ;; Cohen's d - standardized mean difference
 (stats/cohens-d 
  (map first setosa-data)
  (map first versicolor-data))
 
 ;; Hedges' g - less biased alternative to Cohen's d
 (stats/hedges-g
  (map first setosa-data)
  (map first versicolor-data))
 
 ;; Glass's delta - uses only the control group's standard deviation
 (stats/glass-delta
  (map first setosa-data)
  (map first versicolor-data))
 
 ;; Cliff's delta - for ordinal data
 (stats/cliffs-delta
  (map first setosa-data)
  (map first versicolor-data)))

;; Effect sizes for correlations and regression
(utls/examples-note
 ;; R-squared - proportion of variance explained
 (stats/r2-determination [sepal-length petal-length])
 
 ;; Cohen's f-squared - for multiple regression
 (stats/cohens-f2 [sepal-length petal-length])
 
 ;; Cohen's q - comparison of correlations
 (stats/cohens-q
  (stats/correlation sepal-length petal-length)
  (stats/correlation sepal-width petal-width)))

;; [Effect size](https://en.wikipedia.org/wiki/Effect_size) measures complement p-values by quantifying the magnitude of effects, which is often more informative than statistical significance alone.

;; ## Time Series Analysis

;; ### Autocorrelation Analysis

;; Autocorrelation measures correlation between a time series and its lags:

;; Create a simple time series with a pattern
(def time-series (concat (repeat 10 1) (repeat 10 2) (repeat 10 1) (repeat 10 2)))

;; Calculate autocorrelation function (ACF)
(stats/acf time-series 10)

;; Get ACF with confidence intervals
(stats/acf-ci time-series 10)

;; Calculate partial autocorrelation function (PACF)
(stats/pacf time-series 10)

;; Get PACF with confidence intervals
(stats/pacf-ci time-series 10)

;; Calculate Durbin-Watson statistic (test for autocorrelation in residuals)
(stats/durbin-watson time-series)

;; [Autocorrelation](https://en.wikipedia.org/wiki/Autocorrelation) and [partial autocorrelation](https://en.wikipedia.org/wiki/Partial_autocorrelation_function) help identify patterns in time series data, which is crucial for forecasting and model selection.

;; ## Binary Classification Metrics

;; ### Evaluating Classification Performance

;; These metrics assess the performance of binary classification models:

;; Create a confusion matrix from actual and predicted classes
(def confusion-matrix 
  (stats/->confusion-matrix 
   [true true false true false false true true]  ;; actual
   [true false false true false true true false])) ;; predicted

;; Calculate classification metrics
(def metrics (stats/binary-measures-all confusion-matrix))

;; Key performance metrics
(select-keys metrics [:accuracy :precision :recall :f1-score])

;; The full set of metrics
metrics

;; Alternatively, specify true positives, false negatives, etc. directly
(stats/binary-measures 3 2 1 2)

;; Contingency table measures
(stats/contingency-2x2-measures 3 2 1 2)

;; [Classification metrics](https://en.wikipedia.org/wiki/Confusion_matrix) provide a comprehensive evaluation of a classifier's performance, emphasizing different aspects (precision, recall, etc.).

;; ## Statistical Distributions

;; The power-divergence test and its variants help analyze contingency tables and compare distributions:

;; Create a contingency table
(def cont-table {[0 0] 30 [0 1] 10 [1 0] 15 [1 1] 25})

(utls/examples-note
 ;; Power-divergence test
 (stats/power-divergence-test cont-table)
 
 ;; Chi-squared test
 (stats/chisq-test cont-table)
 
 ;; Multinomial likelihood ratio test (G-test)
 (stats/multinomial-likelihood-ratio-test cont-table)
 
 ;; Cressie-Read test
 (stats/cressie-read-test cont-table))

;; KS and AD tests compare samples to distributions
(utls/examples-note
 ;; Kolmogorov-Smirnov test - one sample vs distribution
 (stats/ks-test-one-sample sepal-length r/default-normal)
 
 ;; Kolmogorov-Smirnov test - two samples
 (stats/ks-test-two-samples (map first setosa-data) (map first versicolor-data))
 
 ;; Anderson-Darling test
 (stats/ad-test-one-sample sepal-length r/default-normal))

;; [Goodness-of-fit tests](https://en.wikipedia.org/wiki/Goodness_of_fit) like chi-squared and Kolmogorov-Smirnov determine if data follows a particular distribution or if two samples come from the same distribution.

;; ## Advanced Statistics

;; ### Pooled Statistics

;; Pooled statistics combine data from multiple groups:

(utls/examples-note
 ;; Pooled variance
 (stats/pooled-variance 
  [(map first setosa-data) 
   (map first versicolor-data)])
 
 ;; Pooled standard deviation
 (stats/pooled-stddev
  [(map first setosa-data) 
   (map first versicolor-data)]))

;; [Pooled variance](https://en.wikipedia.org/wiki/Pooled_variance) combines variances from multiple groups, assuming equal variances, useful in t-tests and ANOVA.

;; ### Confidence Intervals

;; Confidence intervals quantify the uncertainty in estimates:

(utls/examples-note
 ;; t-distribution based CI for mean
 (stats/ci sepal-length)
 
 ;; Binomial proportion confidence interval
 (stats/binomial-ci 40 100)
 
 ;; Different methods for binomial CIs
 (stats/binomial-ci 40 100 :wilson)
 (stats/binomial-ci 40 100 :clopper-pearson)
 (stats/binomial-ci 40 100 :agresti-coull))

;; [Confidence intervals](https://en.wikipedia.org/wiki/Confidence_interval) provide a range of plausible values for a parameter, along with a confidence level (typically 95%).

;; ## Putting It All Together: Real-World Examples

;; ### Example 1: Car Performance Analysis

;; Let's analyze the relationship between car weight and miles per gallon from the mtcars dataset:

;; First, let's look at the summary statistics
(stats/stats-map mpg)
(stats/stats-map wt)

;; Calculate the correlation between weight and mpg
(stats/correlation wt mpg)

;; Test if the correlation is statistically significant
;; We can do this by examining the coefficient of determination (r-squared)
(stats/r2-determination [wt mpg])

;; Is there a significant difference in MPG between cars with automatic vs manual transmission?
;; We defined auto-trans-mpg and manual-trans-mpg earlier

;; Run a t-test to compare the means
(stats/t-test-two-samples auto-trans-mpg manual-trans-mpg)

;; Calculate the effect size (Cohen's d)
(stats/cohens-d auto-trans-mpg manual-trans-mpg)

;; Visualize the distribution with a histogram
(def mpg-hist (stats/histogram mpg))
(:bins mpg-hist)

;; Check if MPG follows a normal distribution
(stats/normality-test mpg)

;; ### Example 2: Rock Properties Analysis

;; Now, let's analyze the rock dataset we loaded from a public URL.
;; This dataset contains measurements of 48 rock samples from a petroleum reservoir.

;; First, let's get a statistical summary of the rock permeability
(stats/stats-map rock-perm)

;; And look at the rock shape factor
(stats/stats-map rock-shape)

;; Let's examine the relationship between rock shape and permeability
(stats/correlation rock-shape rock-perm)

;; Is there a significant correlation?
(stats/r2-determination [rock-shape rock-perm])

;; We can divide rocks into groups by permeability to compare their characteristics
(def high-perm-rocks (filter #(> (:perm %) 100) rock-data))
(def low-perm-rocks (filter #(<= (:perm %) 100) rock-data))

;; Compare the area of high vs low permeability rocks
(def high-perm-area (map :area high-perm-rocks))
(def low-perm-area (map :area low-perm-rocks))

;; Run a t-test to see if there's a significant difference
(when (and (seq high-perm-area) (seq low-perm-area))
  (stats/t-test-two-samples high-perm-area low-perm-area))

;; Create a histogram of rock permeability
(def perm-hist (stats/histogram rock-perm))
(:bins perm-hist)

;; ### Example 3: Diamond Price Analysis

;; Finally, let's explore the diamonds dataset to understand price determinants.

;; Get summary statistics for diamond prices
(stats/stats-map diamond-price)

;; And for diamond carat (weight)
(stats/stats-map diamond-carat)

;; Examine the correlation between carat and price
(stats/correlation diamond-carat diamond-price)

;; Calculate the R² to understand how much of price variation is explained by carat
(stats/r2-determination [diamond-carat diamond-price])

;; Let's examine how diamond proportions affect price
;; Calculate the volume (approximated as x * y * z) and analyze its relationship with price
(def diamond-volume (map #(* (:x %) (:y %) (:z %)) diamonds-data))

;; Correlation between volume and price
(stats/correlation diamond-volume diamond-price)

;; Let's check if there's a relationship between table percentage and price
(stats/correlation diamond-table diamond-price)

;; Standardize the price and carat to compare their distributions
(def std-price (stats/standardize diamond-price))
(def std-carat (stats/standardize diamond-carat))

;; How similar are these distributions?
(stats/dissimilarity :jensen-shannon std-price std-carat)

;; Create a histogram of diamond prices
(def price-hist (stats/histogram diamond-price))
(:bins price-hist)

;; Test if diamond prices follow a normal distribution
(stats/normality-test diamond-price)

;; ### Example 4: Cross-Dataset Comparison

;; Let's compare the distribution shapes across our different datasets

;; Calculate skewness for different variables across datasets
(utls/examples-note
 ;; Iris petal width skewness
 (stats/skewness petal-width)
 
 ;; Car mpg skewness
 (stats/skewness mpg)
 
 ;; Rock permeability skewness
 (stats/skewness rock-perm)
 
 ;; Diamond price skewness
 (stats/skewness diamond-price))

;; Calculate kurtosis for different variables across datasets
(utls/examples-note
 ;; Iris petal width kurtosis
 (stats/kurtosis petal-width)
 
 ;; Car mpg kurtosis
 (stats/kurtosis mpg)
 
 ;; Rock permeability kurtosis
 (stats/kurtosis rock-perm)
 
 ;; Diamond price kurtosis
 (stats/kurtosis diamond-price))

;; Compare the coefficient of variation across datasets
(utls/examples-note
 ;; Iris petal width variation
 (stats/variation petal-width)
 
 ;; Car mpg variation
 (stats/variation mpg)
 
 ;; Rock permeability variation
 (stats/variation rock-perm)
 
 ;; Diamond price variation
 (stats/variation diamond-price))

;; ## Additional Statistical Techniques

;; ### Compensation Summation Methods

;; For high precision summation of floating-point numbers, FastMath provides compensation methods:

(utls/examples-note
 ;; Regular sum
 (stats/sum (repeat 1000 0.1))
 
 ;; Kahan summation algorithm - reduces numerical error
 (stats/sum (repeat 1000 0.1) :kahan)
 
 ;; Neumaier algorithm - improved Kahan algorithm
 (stats/sum (repeat 1000 0.1) :neumaier)
 
 ;; Klein algorithm - further improved summation
 (stats/sum (repeat 1000 0.1) :klein))

;; [Kahan summation](https://en.wikipedia.org/wiki/Kahan_summation_algorithm) and related algorithms help reduce numerical errors in floating-point summation by tracking the lost low-order bits.

;; ### Advanced Statistical Tests

;; FastMath includes several tests for checking equality of variances across groups:

;; For these tests, let's define three samples with different dispersions
(def normal1 (repeatedly 50 #(r/grand 100 10))) ;; Normal with SD=10
(def normal2 (repeatedly 50 #(r/grand 100 20))) ;; Normal with SD=20
(def skewed (map #(m/sq %) (repeatedly 50 #(r/grand 0 15)))) ;; Skewed distribution

(utls/examples-note
 ;; Levene's test - tests equality of variances across groups
 ;; Less sensitive to departures from normality
 (stats/levene-test [normal1 normal2 skewed])
 
 ;; Brown-Forsythe test - modification of Levene's test
 ;; Uses median instead of mean, more robust
 (stats/brown-forsythe-test [normal1 normal2 skewed])
 
 ;; Fligner-Killeen test - non-parametric test for homogeneity of variances
 ;; Very robust to departures from normality
 (stats/fligner-killeen-test [normal1 normal2 skewed]))

;; These [tests for homogeneity of variance](https://en.wikipedia.org/wiki/Homogeneity_of_variance) are important for validating assumptions in ANOVA and other statistical procedures.

;; ### Advanced Contingency Table Analysis

;; We can analyze relationships in categorical data through contingency tables:

;; Create a 3x3 contingency table
(def contingency-3x3 {[0 0] 30 [0 1] 10 [0 2] 5
                      [1 0] 15 [1 1] 45 [1 2] 10
                      [2 0] 5  [2 1] 15 [2 2] 25})

;; Extract marginal totals (row and column sums)
(stats/contingency-table->marginals contingency-3x3)

;; Calculate Cramer's V (measure of association)
(utls/examples-note
 ;; Cramer's V - effect size measure for contingency tables
 (stats/cramers-v contingency-3x3)
 
 ;; With bias correction for small samples
 (stats/cramers-v-corrected contingency-3x3)
 
 ;; Tschuprow's T - alternative to Cramer's V
 (stats/tschuprows-t contingency-3x3)
 
 ;; Cohen's kappa - agreement measure for categorical data
 (stats/cohens-kappa {[0 0] 40 [0 1] 10 
                       [1 0] 5 [1 1] 45}))

;; [Contingency table analysis](https://en.wikipedia.org/wiki/Contingency_table) helps understand relationships between categorical variables. Effect size measures like [Cramer's V](https://en.wikipedia.org/wiki/Cram%C3%A9r%27s_V) quantify the strength of these relationships.

;; ### Histogram Features and Bin Selection

;; Creating histograms involves choosing an appropriate number of bins:

(utls/examples-note
 ;; Different bin estimation methods
 (stats/estimate-bins sepal-length :sturges)   ;; Sturges' formula
 (stats/estimate-bins sepal-length :rice)      ;; Rice rule
 (stats/estimate-bins sepal-length :scott)     ;; Scott's rule
 (stats/estimate-bins sepal-length :fd)        ;; Freedman-Diaconis rule
 (stats/estimate-bins sepal-length :sqrt))     ;; Square root rule

;; The choice of binning method affects how the histogram represents the data:
;; - **Sturges' rule** works well for normally distributed data with ~30-200 observations
;; - **Scott's rule** adapts bin width based on the data's standard deviation
;; - **Freedman-Diaconis rule** is more robust against outliers, using the IQR
;; - **Square root rule** is simple, using √n bins for n data points

;; Creating histograms with different bin strategies:
(def hist-sturges (stats/histogram sepal-length :sturges))
(def hist-scott (stats/histogram sepal-length :scott))
(def hist-fd (stats/histogram sepal-length :fd))

;; Compare the bin counts
(map #(count (:bins %)) [hist-sturges hist-scott hist-fd])

;; ### Effect Size Interpretation Guidelines

;; Effect sizes have conventional interpretations for practical significance:

(utls/examples-note
 ;; For Cohen's d and Hedges' g:
 ;;   ~0.2 = small effect
 ;;   ~0.5 = medium effect
 ;;   ~0.8 = large effect
 (let [d (stats/cohens-d sepal-length petal-length)]
   (cond
     (< (m/abs d) 0.2) "Negligible effect"
     (< (m/abs d) 0.5) "Small effect"
     (< (m/abs d) 0.8) "Medium effect"
     :else "Large effect"))
 
 ;; For correlation coefficients (r):
 ;;   ~0.1 = small effect
 ;;   ~0.3 = medium effect
 ;;   ~0.5 = large effect
 (let [r (stats/correlation sepal-length petal-length)]
   (cond
     (< (m/abs r) 0.1) "Negligible correlation"
     (< (m/abs r) 0.3) "Small correlation"
     (< (m/abs r) 0.5) "Medium correlation"
     :else "Large correlation"))
 
 ;; For Cohen's f (ANOVA):
 ;;   ~0.1 = small effect
 ;;   ~0.25 = medium effect
 ;;   ~0.4 = large effect
 (let [f (stats/cohens-f [sepal-length petal-length petal-width])]
   (cond
     (< f 0.1) "Negligible effect"
     (< f 0.25) "Small effect"
     (< f 0.4) "Medium effect"
     :else "Large effect")))

;; ### Specialized Overlap Measures

;; These measures help quantify the overlap or separation between two distributions:

;; Create two groups with partial overlap
(def group-a (repeatedly 50 #(r/grand 100 15)))
(def group-b (repeatedly 50 #(r/grand 130 20)))

(utls/examples-note
 ;; Probability of superiority (a.k.a. Common Language Effect Size)
 ;; Probability that a random value from group B exceeds a random value from group A
 (stats/wmw-odds group-a group-b)
 
 ;; Cohen's U3 - percentage of group A less than the median of group B
 (stats/cohens-u3 group-a group-b)
 
 ;; Probability of overlap - proportion of shared area between distributions
 (stats/p-overlap group-a group-b))

;; These measures provide intuitive interpretations of effect size in terms of distribution overlap, which can be easier to understand than abstract statistical measures.

;; ## Reference

;; For a complete list of available functions in the `fastmath.stats` namespace:

(codox/make-public-fns-table-clay 'fastmath.stats)
