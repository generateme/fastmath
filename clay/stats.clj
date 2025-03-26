(ns stats
  (:require [fastmath.stats :as stats]
            [fastmath.stats.bootstrap :as boot]
            [fastmath.dev.ggplot :as gg]
            [fastmath.dev.codox :as codox]
            [fastmath.dev.clay :as utls]
            [fastmath.dev.dataset :as ds]
            [fastmath.random :as r]
            [fastmath.core :as m]
            [clojure.pprint :as pp]
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

;; ## Basics

;; ::: {.callout-tip title="Defined functions"}
;; * `sum`, `maximum`, `minimum`
;; * `percentile`, `percentiles`, `quantile`, `quantiles`, `wquantile`, `wquantiles`, `pi`
;; * `expectile`
;; * `median`, `median-3`, `wmedian`
;; * `mean`, `geomean`, `harmean`, `powmean`
;; * `variance`, `population-variance`, `wvariance`, `wpopulation-variance`
;; * `stddev`, `population-stddev`, `wstddev`, `wpopultion-stddev`
;; * `pooled-variance`, `pooled-stddev`
;; * `variation`, `l-variation`
;; * `median-absolute-deviation`, `mad`, `mean-absolute-deviation`, `sem`
;; * `mode`, `modes`, `wmode`, `wmodes`
;; * `moment`, `l-moment`
;; * `skewness`, `kurtosis`
;; * `outliers`, `remove-outliers`
;; * `stats-map`
;; :::

;; ### Mean

;; ### Deviation

;; ### Quantiles

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

(def g1 (double-array [1.492282 -0.868731 1.424883 1.113509 0.3186633 0.4629741 1.291994 -0.4128517 0.724347 0.5558941 1.320532 1.13727 1.106706 -0.7290745 1.337817 0.3213244 0.5938525 -0.1433358 -0.7021112 1.203495 0.3286861 0.8253163 0.8289214 0.2008333 0.3994702 1.491169 -0.2580156 -1.227364 0.1812415 0.1493682 2.076148 0.6176723 0.8686895 -0.4396303 -0.6171961 1.596263 -0.7315334 0.8303997 0.6807137 -0.38303 1.374767 1.968743 -1.335768 -1.028934 -0.4564536 0.2601586 -0.1765424 -0.2541709 0.04285849 0.1677205 0.3909273 0.685582 1.107281 -0.1948993 1.614563 2.045843 0.5028711 0.4768879 -1.059179 0.7950796 -1.29715 -1.116765 -0.5888742 1.461546 0.1650975 0.7369283 -0.3403211 0.4315137 0.2894102 0.39839 -1.51201 -2.014084 0.914831 -0.008679345 -0.5604248 0.2005498 -1.039049 0.933337 -0.4382221 1.0336 -0.4258954 1.831763 -1.276483 0.3397851 -1.354471 0.6484439 -0.5189286 0.8063258 1.519506 1.1948 -1.190062 0.2310371 0.3753072 0.6445605 -0.8111008 0.6360471 -0.3907083 -1.983018 -1.441653 -0.2097716]))
(def g2 (double-array [0.8053097 23.47817 24.49539 12.47836 13.3117 12.07374 19.40675 15.15028 11.45262 7.659371 -3.167234 18.82419 3.520448 8.522478 6.691219 7.525913 5.73451 3.936463 22.55598 15.20478 8.030727 1.448121 24.426 15.27894 19.08787 2.361847 17.09298 3.694654 1.914526 17.18825 -11.71995 14.26035 -1.320718 18.60305 8.500502 13.13979 20.71251 28.16612 15.10173 2.473348 -9.78168 31.48709 12.16099 17.09689 11.77993 -8.886036 21.81709 26.34823 35.31201 11.91587 3.962687 14.21603 -6.621846 3.130057 -6.571614 5.247542 24.21332 6.030832 13.60529 11.30866 2.127471 1.524765 10.34462 -3.611007 2.697738 15.49237 -1.229746 8.202995 11.24158 18.87931 7.428771 2.710082 1.721434 20.79805 14.03751 32.51092 7.506515 -2.640752 14.34556 8.090143 4.693302 14.27261 9.123078 9.137117 11.79224 21.37185 13.72709 6.689889 18.63855 4.821303 25.1361 4.713299 15.86299 16.73875 20.54414 21.65894 8.027159 -10.78236 6.016482 17.96367]))

(let [g1 (r/distribution :real-discrete-distribution {:data (repeatedly 1000 r/grand)})
      g2 (r/distribution :real-discrete-distribution {:data (repeatedly 1000 #(r/grand 2 1))})]
  (gg/->file (gg/function (fn [^double p]
                            (let [p (m/constrain p 4.9E-324 0.9999999999999999)
                                  p- (m/- 1.0 p)
                                  q1 (v/vec2 (r/icdf g1 p) (r/icdf g1 p-))
                                  q2 (v/vec2 (r/icdf g2 p-) (r/icdf g2 p))]
                              (-> (v/sub q1 q2) v/abs v/mn))))))

;; ### Mode

;; ### Moment

;; ### Outliers

;; ### Other

;; ## Intervals and extents

;; ::: {.callout-tip title="Defined functions"}
;; * `stddev-extent`, `mad-extent`, `sem-extent`
;; * `percentile-extent`, `quantile-extent`
;; * `pi-extent`, `hpdi-extent`
;; * `iqr`, `adjacent-values`
;; * `inner-fence-extent`, `outer-fence-extent`
;; * `span`, `extent`
;; * `ci`
;; :::


;; ## Metrics

;; ::: {.callout-tip title="Defined functions"}
;; * `me`, `mae`, `mape`, `mse`, `rss`, `r2`, `mse`, `rmse`
;; * `count=`, `L0`, `L1`, `L1sq`, `L2`, `LInf`
;; * `psnr`
;; :::

;; ## Dissimilarity and similarity

;; ::: {.callout-tip title="Defined functions"}
;; * `dissimilarity`
;; * `similarity`
;; :::

;; ## Correlation

;; ::: {.callout-tip title="Defined functions"}
;; * `correlation`, `spearman-correlation`, `pearson-correlation`, `kendall-correlation`, `covariance`
;; * `correlation-matrix`, `coefficient-matrix`, `covariance-matrix`
;; :::

;; ## Contingency table

;; ## Effect size

;; Effect size estimates a strength of the relationship between samples.

;; ::: {.callout-tip title="Defined functions"}
;; * `cohens-d`, `cohens-d-corrected`
;; * `hedges-g`, `hedges-g-corrected`, `hedges-g*`
;; * `glass-delta`
;; * `means-ratio`, `means-ratio-corrected`
;; * `cliffs-delta`
;; * `ameasure`, `wmw-odds`
;; * `p-overlap`
;; * `cohens-u2`, `cohens-u3`
;; * `pearson-r`, `r2-determination`, `eta-sq`, `omega-sq`, `epsilon-sq`
;; * `cohens-f2`, `cohens-f`, `cohens-q`
;; * `rank-eta-sq`, `rank-epsilon-sq`
;; :::

;; ### Difference family

;; Effect sizes based on a difference between means.

(utls/examples-note
  )

;; ## Transformation

;; ::: {.callout-tip title="Defined functions"}
;; * `trim`, `trim-lower`, `trim-upper`, `winsor`
;; * `standardize`, `robust-standardize`, `demean`
;; * `rescale`
;; :::


;; ## Histogram

;; ::: {.callout-tip title="Defined functions"}
;; * `histogram`
;; * `estimate-bins`
;; :::


;; ## Tests

;; ## Bootstrap

;; ## Reference

(codox/make-public-fns-table-clay 'fastmath.stats)
(codox/make-public-fns-table-clay 'fastmath.stats.bootstrap)
