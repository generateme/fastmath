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

(add-examples percentiles
  (example (percentiles [1 2 3 -1 -1 2 -1 11 111] [25 50 75 90])))

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

(add-examples quantiles
  (example (quantiles [1 2 3 -1 -1 2 -1 11 111] [0.25 0.5 0.75 0.9])))

(add-examples median
  (example "Median (percentile 50%)." (median [1 2 3 -1 -1 2 -1 11 111]))
  (example "For three elements use faster [[median-3]]." (median [7 1 4])))

(add-examples median-3 (example "Median of [7 1 4]" (median-3 7 1 4)))

(add-examples mean (example "Mean (average value)" (mean [1 2 3 -1 -1 2 -1 11 111])))
(add-examples geomean (example "Geometric mean" (geomean [1 2 3 1 1 2 1 11 111])))
(add-examples harmean (example "Harmonic mean" (harmean [1 2 3 -1 -1 2 -1 11 111])))
(add-examples powmean (example-session "Power mean"
                        (powmean [1 2 3 1 1 2 1 11 111] 0.0)
                        (powmean [1 2 3 1 1 2 1 11 111] 0.1)
                        (powmean [1 2 3 1 1 2 1 11 111] 0.5)
                        (powmean [1 2 3 1 1 2 1 11 111] 1.0)
                        (powmean [1 2 3 1 1 2 1 11 111] 2.0)
                        (powmean [1 2 3 1 1 2 1 11 111] 3.0)
                        (powmean [1 2 3 1 1 2 1 11 111] 5.5)))


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
(add-examples extent (example "min/max and mean of gaussian distribution" (extent (repeatedly 100000 r/grand))))
(add-examples mad-extent (example "median absolute deviation from median for gaussian distribution" (mad-extent (repeatedly 100000 r/grand))))
(add-examples sem-extent (example "standard error of mean and mean for gaussian distribution" (sem-extent (repeatedly 100000 r/grand))))
(add-examples stddev-extent (example "standard deviation from mean and mean for gaussian distribution" (stddev-extent (repeatedly 100000 r/grand))))
(add-examples percentile-extent (example-session "for samples from gaussian distribution"
                                  (percentile-extent (repeatedly 100000 r/grand))
                                  (percentile-extent (repeatedly 100000 r/grand) 10)
                                  (percentile-extent (repeatedly 100000 r/grand) 30 70)))
(add-examples percentile-bc-extent (example-session "for samples from gaussian distribution"
                                     (percentile-bc-extent (repeatedly 100000 r/grand))
                                     (percentile-bc-extent (repeatedly 100000 r/grand) 10)
                                     (percentile-bc-extent (repeatedly 100000 r/grand) 30 70)))
(add-examples percentile-bca-extent (example-session "for samples from gaussian distribution"
                                      (percentile-bca-extent (repeatedly 100000 r/grand))
                                      (percentile-bca-extent (repeatedly 100000 r/grand) 10)
                                      (percentile-bca-extent (repeatedly 100000 r/grand) 30 70)))

(add-examples sum (example "Sum" (sum [1 2 3 -1 -1 2 -1 11 111])))

(add-examples kurtosis
  (example-session "Kurtosis"
    (kurtosis [1 2 3 -1 -1 2 -1 11 111])
    (kurtosis [1 2 3 -1 -1 2 -1 11 111] :G2)
    (kurtosis [1 2 3 -1 -1 2 -1 11 111] :g2)
    (kurtosis [1 2 3 -1 -1 2 -1 11 111] :excess)
    (kurtosis [1 2 3 -1 -1 2 -1 11 111] :kurt)))

(add-examples moment
  (example-session "Usage"
    (moment [3 7 5 9 -8])
    (moment [3 7 5 9 -8] 1.0)
    (moment [3 7 5 9 -8] 4.0)
    (moment [3 7 5 9 -8] 3.0)
    (moment [3 7 5 9 -8] 3.0 {:center 0.0})
    (moment [3 7 5 9 -8] 3.0 {:mean? false})
    (moment [3 7 5 9 -8] 3.0 {:absolute? true})
    (moment [3 7 5 9 -8] 3.0 {:center -3.0})
    (moment [3 7 5 9 -8] 0.5 {:absolute? true})))

(add-examples skewness
  (example-session "Skewness"
    (skewness [1 2 3 -1 -1 2 -1 11 111])
    (skewness [1 2 3 -1 -1 2 -1 11 111] :G1)
    (skewness [1 2 3 -1 -1 2 -1 11 111] :g1)
    (skewness [1 2 3 -1 -1 2 -1 11 111] :pearson)
    (skewness [1 2 3 -1 -1 2 -1 11 111] :b1)
    (skewness [1 2 3 -1 -1 2 -1 11 111] :B1)
    (skewness [1 2 3 -1 -1 2 -1 11 111] :yule)
    (skewness [1 2 3 -1 -1 2 -1 11 111] :B3)
    (skewness [1 2 3 -1 -1 2 -1 11 111] :mode)
    (skewness [1 2 3 -1 -1 2 -1 11 111] :median)
    (skewness [1 2 3 -1 -1 2 -1 11 111] :skew)))

(add-examples sem (example "SEM" (sem [1 2 3 -1 -1 2 -1 11 111])))

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

(add-examples acf
  (example-session "Usage"
    (acf (repeatedly 1000 r/grand) 5)
    (acf (repeatedly 1000 r/grand) [10 20 100 500])
    (acf [1 2 3 4 5 4 3 2 1])))

(add-examples acf-ci
  (example-session "Usage"
    (acf-ci (repeatedly 1000 r/grand) 3)
    (acf-ci [1 2 3 4 5 4 3 2 1] 3)))

(add-examples pacf
  (example-session "Usage"
    (pacf (repeatedly 1000 r/grand) 10)
    (pacf [1 2 3 4 5 4 3 2 1])))

(add-examples pacf-ci
  (example-session "Usage"
    (pacf-ci (repeatedly 1000 r/grand) 3)
    (pacf-ci [1 2 3 4 5 4 3 2 1] 3)))

(add-examples ameasure
  (example (let [t [10,10,20,20,20,30,30,30,40,50]
                 c [10,20,30,40,40,50]]
             (ameasure t c))))

(add-examples cohens-d
  (example (let [t [10,10,20,20,20,30,30,30,40,50]
                 c [10,20,30,40,40,50]]
             (cohens-d t c))))

(add-examples cohens-d-corrected
  (example (let [t [10,10,20,20,20,30,30,30,40,50]
                 c [10,20,30,40,40,50]]
             (cohens-d-corrected t c))))

(add-examples glass-delta
  (example (let [t [10,10,20,20,20,30,30,30,40,50]
                 c [10,20,30,40,40,50]]
             (glass-delta t c))))

(add-examples cliffs-delta
  (example-session "Usage"
    (let [t [10,10,20,20,20,30,30,30,40,50]
          c [10,20,30,40,40,50]]
      (cliffs-delta t c))
    (let [t [:a :b :c :D :e :f]
          c [:a :z :X :y :x :y]]
      (cliffs-delta t c))))

(add-examples hedges-g
  (example (let [t [10,10,20,20,20,30,30,30,40,50]
                 c [10,20,30,40,40,50]]
             (hedges-g t c))))

(add-examples hedges-g-corrected
  (example (let [t [10,10,20,20,20,30,30,30,40,50]
                 c [10,20,30,40,40,50]]
             (hedges-g-corrected t c))))

(add-examples hedges-g*
  (example (let [t [10,10,20,20,20,30,30,30,40,50]
                 c [10,20,30,40,40,50]]
             (hedges-g* t c))))

(add-examples pearson-r
  (example (let [t [10,10,20,20,20,30,30,30,40,50]
                 c [-50,20,30,40,40,50,10,20,30,10]]
             (pearson-r t c))))

(add-examples r2-determination
  (example (let [t [10,10,20,20,20,30,30,30,40,50]
                 c [-50,20,30,40,40,50,10,20,30,10]]
             (r2-determination t c))))

(add-examples eta-sq
  (example (let [t [10,10,20,20,20,30,30,30,40,50]
                 c [-50,20,30,40,40,50,10,20,30,10]]
             (eta-sq t c))))

(add-examples omega-sq
  (example (let [t [10,10,20,20,20,30,30,30,40,50]
                 c [-50,20,30,40,40,50,10,20,30,10]]
             (omega-sq t c))))

(add-examples epsilon-sq
  (example (let [t [10,10,20,20,20,30,30,30,40,50]
                 c [-50,20,30,40,40,50,10,20,30,10]]
             (epsilon-sq t c))))

(add-examples cohens-f2
  (example (let [t [10,10,20,20,20,30,30,30,40,50]
                 c [-50,20,30,40,40,50,10,20,30,10]]
             {:default (cohens-f2 t c)
              :eta (cohens-f2 :eta t c)
              :omega (cohens-f2 :omega t c)
              :epsilon (cohens-f2 :epsilon t c)})))

(add-examples cohens-q
  (example (let [t [10,10,20,20,20,30,30,30,40,50]
                 c [-50,20,30,40,40,50,10,20,30,10]
                 d [5 2 3 4 4 5 4 2 3 -1]
                 e (range 10)]
             {:arity-2 (cohens-q 0.5 -0.25)
              :arity-3 (cohens-q t c d)
              :arity-4 (cohens-q t c d e)})))

(add-examples cramers-v
  (example (let [a [:a :a :b :b :f :a :a :b :b :c :a :a :b :b :c :a :a :b :b :c]
                 b [:b :f :a :a :b :b :y :z :c :b :b :c :a :a :b :b :c :a :a :b]]
             (cramers-v a b))))

(add-examples cramers-v-corrected
  (example (let [a [:a :a :b :b :f :a :a :b :b :c :a :a :b :b :c :a :a :b :b :c]
                 b [:b :f :a :a :b :b :y :z :c :b :b :c :a :a :b :b :c :a :a :b]]
             (cramers-v-corrected a b))))

(add-examples cohens-w
  (example (let [a [:a :a :b :b :f :a :a :b :b :c :a :a :b :b :c :a :a :b :b :c]
                 b [:b :f :a :a :b :b :y :z :c :b :b :c :a :a :b :b :c :a :a :b]]
             (cohens-w a b))))

(add-examples tschuprows-t
  (example (let [a [:a :a :b :b :f :a :a :b :b :c :a :a :b :b :c :a :a :b :b :c]
                 b [:b :f :a :a :b :b :y :z :c :b :b :c :a :a :b :b :c :a :a :b]]
             (tschuprows-t a b))))

(add-examples binary-measures
  (example (binary-measures [true false true false true false true false]
                            [true false false true false false false true]))
  (example "Treat `1` as `true` value." (binary-measures [1 0 1 0 1 0 1 0]
                                                         [1 0 0 1 0 0 0 1] [1]))
  (example "Treat `:a` and `:b` as `true` value." (binary-measures [:a :b :c :d :e :f :a :b]
                                                                   [:a :b :a :b :a :f :d :b] {:a true
                                                                                              :b true
                                                                                              :c false})))

(add-examples binary-measures-all
  (example (binary-measures-all [true false true false true false true false]
                                [true false false true false false false true]))
  (example "Treat `1` as `true` value." (binary-measures-all [1 0 1 0 1 0 1 0]
                                                             [1 0 0 1 0 0 0 1] [1]))
  (example "Treat `:a` and `:b` as `true` value." (binary-measures-all [:a :b :c :d :e :f :a :b]
                                                                       [:a :b :a :b :a :f :d :b] {:a true
                                                                                                  :b true
                                                                                                  :c false}))
  (example "F-beta is a function. When `beta` is equal `1.0`, you get `f1-score`."
    (let [fbeta (:f-beta (binary-measures-all [true false true false true false true false]
                                              [true false false true false false false true]))]
      [(fbeta 1.0)
       (fbeta 2.0)
       (fbeta 0.5)])))

(add-examples bootstrap
  (example-session "Usage"
    (bootstrap [1 2 3 4 1 2 3 1 2 1] 2 20)
    (let [data [1 2 3 4 1 2 3 1 2 1]
          fdata (frequencies data)
          bdata (bootstrap data 5 1000)]
      {:source fdata
       :bootstrapped (map frequencies bdata)})))

(add-examples bootstrap-ci
  (example-session "Usage"
    (bootstrap-ci [-5 1 1 1 1 2 2 5 11 71])
    (bootstrap-ci [-5 1 1 1 1 2 2 5 11 71] 0.8)
    (bootstrap-ci [-5 1 1 1 1 2 2 5 11 71] 0.8 100000)
    (bootstrap-ci [-5 1 1 1 1 2 2 5 11 71] 0.98 1000 median)))

(add-examples ci
  (example-session "Usage"
    (ci [-5 1 1 1 1 2 2 5 11 71])
    (ci [-5 1 1 1 1 2 2 5 11 71] 0.8)))

(add-examples covariance-matrix
  (example (covariance-matrix [[1 2 3 4 5 11] [3 2 3 2 3 4]])))

(add-examples demean
  (example (demean [-5 1 1 1 1 2 2 5 11 71])))

(add-examples ttest-one-sample
  (example-session "Usage"
    (ttest-one-sample [1 2 3 4 5 6 7 8 9 10])
    (ttest-one-sample [1 2 3 4 5 6 7 8 9 10] {:alpha 0.2})
    (ttest-one-sample [1 2 3 4 5 6 7 8 9 10] {:sides :one-sided})
    (ttest-one-sample [1 2 3 4 5 6 7 8 9 10] {:mu 5.0})))

(add-examples ttest-two-samples
  (example-session "Usage"
    (ttest-two-samples [1 2 3 4 5 6 7 8 9 10] [7 8 9 10 11 12 13 14 15 16 17 18 19 20])
    (ttest-two-samples [1 2 3 4 5 6 7 8 9 10] [7 8 9 10 11 12 13 14 15 16 17 18 19 20 200])
    (ttest-two-samples [1 2 3 4 5 6 7 8 9 10] [7 8 9 10 11 12 13 14 15 16 17 18 19 20] {:equal-variances? true})
    (ttest-two-samples [1 2 3 4 5 6 7 8 9 10] [200 11 200 11 200 11 200 11 200 11] {:paired? true})))
