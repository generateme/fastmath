^{:nextjournal.clerk/visibility :hide-ns
  :nextjournal.clerk/toc true}
(ns stats
  {:clj-kondo/config '{:config-in-call {utils/table {:ignore [:unresolved-symbol]}
                                        utils/table2 {:ignore [:unresolved-symbol]}}}}
  (:require [nextjournal.clerk :as clerk]
            [utils :as u]
            [fastmath.stats :as stats]))

;; # fastmath.stats

;; Statistical functions

#_:clj-kondo/ignore
(require '[fastmath.stats :as stats])

;; ## Descriptive statistics

;; ### mean

;; ### measure

;; ### deviation

;; ### quantile

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
 (stats/quantile data 0.69)
 (map (fn [es] (stats/quantile data 0.69 es)) [:legacy :r1 :r2 :r3 :r4 :r5 :r6 :r7 :r8 :r9])
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

;; ### mode

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

;; ### Intervals

;; #### Binomial

;; ## Correlation

;; ### ACF

;; ### PACF

;; ## Metrics

;; ## Effect Size

;; ## Data manipulation

;; ### Bootstrap

;; ## Confusion matrix

;; ## Histogram

;; ## Tests

;; ## List of symbols

^{::clerk/visibility :hide}
(u/make-public-fns-table 'fastmath.stats)
