(ns fastmath.stats-examples
  (:require [metadoc.examples :refer :all]
            [fastmath.stats :refer :all]
            [fastmath.random :as r]
            [fastmath.core :as m]))

(add-examples mode
  (example "Example" (mode [1 2 3 -1 -1 2 -1 11 111]))
  (example "Returns lowest value when every element appears equally." (mode [5 1 2 3 4])))

(add-examples modes
  (example "Example" (modes [1 2 3 -1 -1 2 -1 11 111]))
  (example "Returns lowest value when every element appears equally." (modes [5 5 1 1 2 3 4 4])))

(add-examples estimation-strategies-list
  (example "List of estimation strategies for percentile" (sort (keys estimation-strategies-list))))

(add-examples percentile
  (example "Percentile 25%" (percentile [1 2 3 -1 -1 2 -1 11 111] 25.0))
  (example "Percentile 50% (median)" (percentile [1 2 3 -1 -1 2 -1 11 111] 50.0))
  (example "Percentile 75%" (percentile [1 2 3 -1 -1 2 -1 11 111] 75.0))
  (example "Percentile 90%" (percentile [1 2 3 -1 -1 2 -1 11 111] 90.0))
  (example-session "Various estimation strategies."
    (percentile [1 2 3 -1 -1 2 -1 11 111] 85.0 :legacy)
    (percentile [1 2 3 -1 -1 2 -1 11 111] 85.0 :r1)
    (percentile [1 2 3 -1 -1 2 -1 11 111] 85.0 :r2)
    (percentile [1 2 3 -1 -1 2 -1 11 111] 85.0 :r3)
    (percentile [1 2 3 -1 -1 2 -1 11 111] 85.0 :r4)
    (percentile [1 2 3 -1 -1 2 -1 11 111] 85.0 :r5)
    (percentile [1 2 3 -1 -1 2 -1 11 111] 85.0 :r6)
    (percentile [1 2 3 -1 -1 2 -1 11 111] 85.0 :r7)
    (percentile [1 2 3 -1 -1 2 -1 11 111] 85.0 :r8)
    (percentile [1 2 3 -1 -1 2 -1 11 111] 85.0 :r9)))

(add-examples quantile
  (example "Quantile 0.25" (quantile [1 2 3 -1 -1 2 -1 11 111] 0.25))
  (example "Quantile 0.5 (median)" (quantile [1 2 3 -1 -1 2 -1 11 111] 0.5))
  (example "Quantile 0.75" (quantile [1 2 3 -1 -1 2 -1 11 111] 0.75))
  (example "Quantile 0.9" (quantile [1 2 3 -1 -1 2 -1 11 111] 0.9))
  (example-session "Various estimation strategies."
    (quantile [1 11 111 1111] 0.7 :legacy)
    (quantile [1 11 111 1111] 0.7 :r1)
    (quantile [1 11 111 1111] 0.7 :r2)
    (quantile [1 11 111 1111] 0.7 :r3)
    (quantile [1 11 111 1111] 0.7 :r4)
    (quantile [1 11 111 1111] 0.7 :r5)
    (quantile [1 11 111 1111] 0.7 :r6)
    (quantile [1 11 111 1111] 0.7 :r7)
    (quantile [1 11 111 1111] 0.7 :r8)
    (quantile [1 11 111 1111] 0.7 :r9)))

(add-examples median
  (example "Median (percentile 50%)." (median [1 2 3 -1 -1 2 -1 11 111]))
  (example "For three elements use faster [[median-3]]." (median [7 1 4])))

(add-examples median-3 (example "Median of [7 1 4]" (median-3 7 1 4)))
(add-examples mean (example "Mean (average value)" (mean [1 2 3 -1 -1 2 -1 11 111])))
(add-examples population-variance (example "Population variance" (population-variance [1 2 3 -1 -1 2 -1 11 111])))
(add-examples variance (example "Variance." (variance [1 2 3 -1 -1 2 -1 11 111])))
(add-examples population-stddev (example "Population standard deviation." (population-stddev [1 2 3 -1 -1 2 -1 11 111])))
(add-examples stddev (example "Standard deviation." (stddev [1 2 3 -1 -1 2 -1 11 111])))
(add-examples median-absolute-deviation (example "MAD" (median-absolute-deviation [1 2 3 -1 -1 2 -1 11 111])))

(add-examples adjacent-values
  (example "[LAV, UAV]" (adjacent-values [1 2 3 -1 -1 2 -1 11 111]))
  (example "Gaussian distribution [LAV, UAV]" (adjacent-values (repeatedly 1000000 r/grand))))

(add-examples iqr
  (example "IQR" (iqr (repeatedly 100000 r/grand))))

(add-examples outliers
              (example "Outliers" (outliers [1 2 3 -1 -1 2 -1 11 111]))
              (example "Gaussian distribution outliers" (count (outliers (repeatedly 3000000 r/grand)))))

(add-examples minimum (example "Minimum value" (minimum [1 2 3 -1 -1 2 -1 11 111])))
(add-examples maximum (example "Maximum value" (maximum [1 2 3 -1 -1 2 -1 11 111])))
(add-examples extent (example "min/max from gaussian distribution" (extent (repeatedly 100000 r/grand))))
(add-examples sum (example "Sum" (sum [1 2 3 -1 -1 2 -1 11 111])))

(add-examples kurtosis (example "Kurtosis" (kurtosis [1 2 3 -1 -1 2 -1 11 111])))
(add-examples second-moment (example "Second Moment" (second-moment [1 2 3 -1 -1 2 -1 11 111])))
(add-examples skewness (example "Skewness" (skewness [1 2 3 -1 -1 2 -1 11 111])))

(add-examples stats-map (example "Stats" (stats-map [1 2 3 -1 -1 2 -1 11 111])))

(add-examples standardize (example "Standardize" (standardize [1 2 3 -1 -1 2 -1 11 111])))

(add-examples covariance (example "Covariance of uniform and gaussian distribution samples." (covariance (repeatedly 100000 (partial r/grand 1.0 10.0)) (repeatedly 100000 (partial r/drand -10.0 -5.0)))))
(add-examples correlation (example "Correlation of uniform and gaussian distribution samples." (correlation (repeatedly 100000 (partial r/grand 1.0 10.0)) (repeatedly 100000 (partial r/drand -10.0 -5.0)))))
(add-examples spearman-correlation (example "Spearsman's correlation of uniform and gaussian distribution samples." (spearman-correlation (repeatedly 100000 (partial r/grand 1.0 10.0)) (repeatedly 100000 (partial r/drand -10.0 -5.0)))))
(add-examples pearson-correlation (example "Pearson's correlation of uniform and gaussian distribution samples." (pearson-correlation (repeatedly 100000 (partial r/grand 1.0 10.0)) (repeatedly 100000 (partial r/drand -10.0 -5.0)))))
(add-examples kendall-correlation (example "Kendall's correlation of uniform and gaussian distribution samples." (kendall-correlation (repeatedly 100000 (partial r/grand 1.0 10.0)) (repeatedly 100000 (partial r/drand -10.0 -5.0)))))

(add-examples kullback-leibler-divergence (example "Kullback-Leibler divergence." (kullback-leibler-divergence (repeatedly 100 #(r/irand 100)) (repeatedly 100 #(r/irand 100)))))
(add-examples jensen-shannon-divergence (example "Jensen-Shannon divergence" (jensen-shannon-divergence (repeatedly 100 #(r/irand 100)) (repeatedly 100 #(r/irand 100)))))

(add-examples estimate-bins
  (let [d (r/distribution :log-normal)
        vs (repeatedly 1000 #(r/drandom d))]
    (example-session "Estimate number of bins for various methods. `vs` contains 1000 random samples from Log-Normal distribution."
      (estimate-bins vs :sqrt)
      (estimate-bins vs :sturges)
      (estimate-bins vs :rice)
      (estimate-bins vs :doane)
      (estimate-bins vs :scott)
      (estimate-bins vs :freedman-diaconis))))

(add-examples histogram
  (example "3 bins from uniform distribution." (histogram (repeatedly 1000 rand) 3))
  (example "3 bins from uniform distribution for given range." (histogram (repeatedly 10000 rand) 3 [0.1 0.5]))
  (example "5 bins from normal distribution." (histogram (repeatedly 10000 r/grand) 5))
  (example "Estimate number of bins" (:size (histogram (repeatedly 10000 r/grand))))
  (example "Estimate number of bins, Rice rule" (:size (histogram (repeatedly 10000 r/grand) :rice))))

(add-examples kernel-density
  (example (let [kd (kernel-density [0 10 10 10 10 10 10 10 10 0 0 0 0 1 1 1] 1)]
             (map (comp m/approx kd) (range -5 15)))))
