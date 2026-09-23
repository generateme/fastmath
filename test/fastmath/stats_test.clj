(ns fastmath.stats-test
  (:require [fastmath.stats :as sut]
            [fastmath.stats.binary :as bm]
            [clojure.test :as t]
            [fastmath.core :as m]
            [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [fastmath.vector :as v]
            [fastmath.random :as r]))

(defn transform [spec data]
  (if (sequential? spec)
    (let [[nm f] spec] [nm (f data)])
    [spec (read-string data)]))
(defn zip [header row] (into {} (map transform header row)))
(defn parse-rows  [header rows]  (map (partial zip header) rows))
(defn read-csv
  [fname header]
  (->> (io/resource fname)
       (slurp)
       (csv/read-csv)
       (rest)
       (parse-rows header)))

(defn data->fn
  [data]
  (memoize (fn
             ([] data)
             ([selector] (map selector data))
             ([selector filter-pred] (map selector (filter filter-pred data))))))

(defn by
  ([data f]
   (mapv second (sort-by first (group-by f (data)))))
  ([data f selector]
   (vec (for [group (by data f)]
          (map selector group)))))

(def iris (data->fn (read-csv "iris.csv"
                            [:sepal-length :sepal-width :petal-length :petal-width [:species keyword]])))
(def mtcars (data->fn (read-csv "mtcars.csv"
                              [[:name identity] :mpg :cyl :disp :hp :drat :wt :qsec :vs :am :gear :carb])))


;; basic
;; reference values from R: min(mtcars$mpg), max(mtcars$mpg), sum(mtcars$mpg)

(t/deftest basic-test
  (t/are [f res] (m/delta= res (f (mtcars :mpg)))
    sut/minimum 10.4
    sut/maximum 33.9
    sut/sum 642.9)
  (t/testing "sum compensation methods"
    (t/testing "cancellation-then-regrowth: classic Kahan does not recover, Neumaier/Klein do"
      ;; true sum = 2.0, confirmed via Python's math.fsum([1e16, 1.0, -1e16, 1.0])
      (let [xs [1e16 1.0 -1e16 1.0]]
        (t/is (m/delta-eq 1.0 (sut/sum xs)) "plain summation loses both +1.0 terms")
        (t/is (m/delta-eq 1.0 (sut/sum xs :kahan))
              "classic Kahan's compensation term is itself swallowed by the -1e16 cancellation step (textbook limitation motivating Neumaier's variant), not a fastmath bug")
        (t/is (m/delta-eq 2.0 (sut/sum xs :neumayer)))
        (t/is (m/delta-eq 2.0 (sut/sum xs :klein)))))
    (t/testing "many sub-ulp increments: all three compensated methods recover, plain does not"
      ;; true sum = 1.000000000001, confirmed via Python's math.fsum([1.0] + [1e-16]*10000)
      (let [ys (into [1.0] (repeat 10000 1e-16))]
        (t/is (m/delta-eq 1.0 (sut/sum ys)) "plain summation: each 1e-16 increment is below 1.0's ulp")
        (t/is (m/delta-eq 1.000000000001 (sut/sum ys :kahan)))
        (t/is (m/delta-eq 1.000000000001 (sut/sum ys :neumayer)))
        (t/is (m/delta-eq 1.000000000001 (sut/sum ys :klein)))))))

;; percentiles / quantiles
;; reference values from R: quantile(mtcars$mpg, probs, type=1:9), one per Hyndman-Fan
;; estimation strategy (fastmath's :r1..:r9); weighted-quantile reference values from R:
;; spatstat.geom::weighted.quantile(vs, w=ws, probs, type=4) for :linear, type=1 for :step

(def qprobs [10.0 25.0 50.0 75.0 90.0])

(t/deftest percentiles-quantiles-test
  (t/testing "percentile matches R's quantile() for all 9 Hyndman-Fan types, :legacy == :r6"
    (t/are [strategy res] (v/delta-eq res (mapv #(sut/percentile (mtcars :mpg) % strategy) qprobs))
      :r1     [14.3 15.2 19.2 22.8 30.4]
      :r2     [14.3 15.35 19.2 22.8 30.4]
      :r3     [13.3 15.2 19.2 22.8 30.4]
      :r4     [13.5 15.2 19.2 22.8 29.78]
      :r5     [14.0 15.35 19.2 22.8 30.4]
      :r6     [13.6 15.275 19.2 22.8 30.4]
      :r7     [14.34 15.425 19.2 22.8 30.09]
      :r8     [13.866666666666667 15.325 19.2 22.8 30.4]
      :r9     [13.9 15.331249999999999 19.2 22.8 30.4]
      :legacy [13.6 15.275 19.2 22.8 30.4]))
  (t/testing "boundary: p=0/p=100 equal minimum/maximum regardless of strategy"
    (t/is (m/delta= 10.4 (sut/percentile (mtcars :mpg) 0.0)))
    (t/is (m/delta= 33.9 (sut/percentile (mtcars :mpg) 100.0)))
    (t/is (m/delta= 33.9 (sut/percentile (mtcars :mpg) 100.0 :r7))))
  (t/testing "quantile (0.0-1.0 scale) agrees with percentile (0-100 scale)"
    (t/are [q] (m/delta= (sut/percentile (mtcars :mpg) (* q 100.0) :r7)
                         (sut/quantile (mtcars :mpg) q :r7))
      0.1 0.25 0.5 0.75 0.9))
  (t/testing "percentiles/quantiles (vector forms) agree with repeated single-value calls"
    (t/is (= (mapv #(sut/percentile (mtcars :mpg) % :r7) qprobs)
             (sut/percentiles (mtcars :mpg) qprobs :r7)))
    (t/is (= (mapv #(sut/quantile (mtcars :mpg) % :r7) [0.1 0.25 0.5 0.75 0.9])
             (sut/quantiles (mtcars :mpg) [0.1 0.25 0.5 0.75 0.9] :r7)))
    (t/is (= (sut/percentiles (mtcars :mpg)) (sut/percentiles (mtcars :mpg) [25 50 75 100]))
          "default probs")
    (t/is (= (sut/quantiles (mtcars :mpg)) (sut/quantiles (mtcars :mpg) [0.25 0.5 0.75 1.0]))
          "default quantiles"))
  (t/testing "median"
    (t/is (m/delta= 19.2 (sut/median (mtcars :mpg))) "R: median(mtcars$mpg)")
    (t/is (m/delta= 19.2 (sut/median (mtcars :mpg) :r7))
          "default (:legacy) agrees with R's type=7 default here")
    (t/is (m/delta= (sut/percentile (mtcars :mpg) 50.0) (sut/median (mtcars :mpg)))))
  (t/testing "median-3"
    (t/are [a b c res] (m/delta= res (sut/median-3 a b c))
      1.0 2.0 3.0 2.0
      3.0 1.0 2.0 2.0
      2.0 3.0 1.0 2.0
      5.0 5.0 1.0 5.0
      -1.0 -5.0 -3.0 -3.0))
  (t/testing "weighted quantiles: :linear and :step methods"
    (let [vs [2.0 5.0 3.0 8.0 1.0 9.0 4.0]
          ws [1.0 2.0 0.5 1.5 3.0 1.0 2.0]]
      (t/are [method res] (v/delta-eq res (sut/wquantiles vs ws [0.1 0.25 0.5 0.75 0.9] method))
        :linear [1.0 1.0 3.5 4.875 7.8]
        :step   [1.0 1.0 4.0 5.0 8.0])
      (t/is (m/delta= 3.5 (sut/wmedian vs ws)) "0.5 quantile, :linear")
      (t/is (= (sut/wquantiles vs ws) (sut/wquantiles vs ws [0.25 0.5 0.75 1.0]))
            "default probs")
      (t/is (= (sut/wquantiles vs ws [0.1 0.25 0.5 0.75 0.9])
               (mapv #(sut/wquantile vs ws %) [0.1 0.25 0.5 0.75 0.9]))
            "vector form agrees with repeated single-value calls")
      ;; :average ("midpoint of step-before/step-after") intentionally left untested here:
      ;; no equivalent found in R's spatstat.geom::weighted.quantile type options (1-4 tried,
      ;; none matched); see the audit note for this group.
      )))

;; means
;; reference values from R: mean(carb), exp(mean(log(carb))), 1/mean(1/carb), mean(carb^p)^(1/p)

(def wts [1.0 3 0.5 1.5 1 2 2.5 1 1 3 4 2 1 1 2 3 1 1 2 1 1 2 1 1 1 2 1 1 1 1 1 2])

(t/deftest means-test
  (t/are [f res] (m/delta= res (f (mtcars :carb)))
    sut/mean 2.8125
    sut/geomean 2.395987
    sut/harmean 2.026385)
  (t/are [power res] (m/delta= res (sut/powmean (mtcars :carb) power))
    -10.0 1.163977
    -1.0 2.026385
    0.0 2.395987
    0.1 2.43644
    0.5 2.60179
    1.0 2.8125
    2.0 3.230712
    10.0 5.69326)
  (t/testing "weighted mean and power mean"
    ;; reference: R weighted.mean(mtcars$mpg, wts), weighted.mean(mtcars$mpg^p, wts)^(1/p)
    (t/is (m/delta= 19.3484848485 (sut/mean (mtcars :mpg) wts)))
    (t/is (m/delta= 20.1893136049 (sut/powmean (mtcars :mpg) wts 2.0)))
    (t/is (m/delta= 17.7101894299 (sut/powmean (mtcars :mpg) wts -1.0))))
  (t/testing "wmean (deprecated, delegates to mean)"
    (t/is (= (sut/mean (mtcars :mpg)) (sut/wmean (mtcars :mpg))))
    (t/is (= (sut/mean (mtcars :mpg) wts) (sut/wmean (mtcars :mpg) wts))))
  (t/testing "logmean"
    ;; 2-value case is the classical closed form L(x,y)=(x-y)/(ln x - ln y);
    ;; reference: R (30-20)/(log(30)-log(20))
    (t/is (m/delta-eq 24.663034623 (sut/logmean [30 20])))
    ;; paper (n=200000 case; algorithm-internal-consistency coverage lives in
    ;; test/fastmath/stats/logmean_test.clj, out of this audit's scope)
    (t/is (m/delta-eq 73578.65538616560 (sut/logmean (range 1 200001))))
    (t/is (m/nan? (sut/logmean nil)))
    (t/is (m/== 1.0123 (sut/logmean [1.0123])))
    (t/is (m/== (sut/logmean [4 2]) (m// 2.0 m/M_LN2)))))

;; deviations
;; reference values from R: var/sd (+ n/(n-1) rescale for population variant),
;; sd(x)/mean(x), mad(x, constant=1), mean(abs(x-center)), sd(x)/sqrt(n)

(t/deftest deviations-test
  (t/are [f res] (m/delta= res (f (mtcars :drat)))
    sut/population-variance 0.2769476
    sut/variance 0.2858814
    sut/population-stddev 0.5262581
    sut/stddev 0.5346787
    sut/variation 0.1486638
    sut/mean-absolute-deviation 0.4532422
    sut/median-absolute-deviation 0.475
    sut/sem 0.09451874)
  (let [d (mtcars :drat)]
    (t/is (m/delta= 0.4459375 (sut/mean-absolute-deviation d (sut/median d))))
    (t/is (m/delta= 3.5965625 (sut/mean-absolute-deviation d 0.0)))
    (t/is (m/delta= 0.465 (sut/median-absolute-deviation d (sut/mean d))))
    (t/is (m/delta= 0.46 (sut/median-absolute-deviation d nil :r1)))
    (t/is (m/delta= 3.695 (sut/median-absolute-deviation d 0))))
  (t/testing "weighted variance/stddev, population (biased) and sample (unbiased)"
    ;; reference values from R Hmisc::wtd.var(x, w, method="ML") (population) and
    ;; method="unbiased" (sample), cross-checked against the manual formula
    ;; sum(w*(x-mu)^2)/sum(w) resp. /(sum(w)-1)
    (let [wx [21.0 22.8 18.7 24.4 14.3 19.2 17.8]
          ww [2.0 1 3 0.5 2.5 1 1]]
      (t/is (m/delta= 8.462995867768592 (sut/population-wvariance wx ww)))
      (t/is (m/delta= 9.309295454545452 (sut/wvariance wx ww)))
      (t/is (m/delta= 2.9091228691426205 (sut/population-wstddev wx ww)))
      (t/is (m/delta= 3.0511138055709184 (sut/wstddev wx ww))))))

;; extent family (mean/median +/- dispersion, percentile/quantile pairs, Tukey fences, HPDI)
;; reference values from R on mtcars$mpg: mean(x)+/-sd(x), median(x)+/-mad(x,constant=1),
;; mean(x)+/-sd(x)/sqrt(n), quantile(x, type=6) (== :legacy, see Percentiles & Quantiles group),
;; Tukey fences/adjacent-values from Q1/Q3 (type=6) +/- 1.5*IQR / 3*IQR, and a from-scratch R
;; reimplementation of the HPDI sliding-window-of-width-`gap`-order-statistics algorithm

(t/deftest extent-family-test
  (let [mpg (mtcars :mpg)]
    (t/is (v/delta-eq [14.0636769479 26.1175730521 20.090625] (sut/stddev-extent mpg)))
    (t/is (v/delta-eq [15.55 22.85 19.2] (sut/mad-extent mpg)))
    (t/is (v/delta-eq [19.0252010406 21.1560489594 20.090625] (sut/sem-extent mpg)))
    (t/is (v/delta-eq [15.275 22.8 19.2] (sut/percentile-extent mpg)))
    (t/is (v/delta-eq [15.275 22.8 19.2] (sut/quantile-extent mpg)))
    (t/is (= {25.0 15.274999999999999 75.0 22.8} (sut/pi mpg)))
    (t/is (v/delta-eq [15.275 22.8 19.2] (sut/pi-extent mpg)))
    (t/is (v/delta-eq [10.4 32.4 19.2] (sut/hpdi-extent mpg)) "size=0.95, gap=30 order statistics")
    (t/is (m/delta= 7.525 (sut/iqr mpg)))
    (t/is (v/delta-eq [3.9875 34.0875 19.2] (sut/inner-fence-extent mpg)))
    (t/is (v/delta-eq [-7.3 45.375 19.2] (sut/outer-fence-extent mpg)))
    (t/is (m/delta= 23.5 (sut/span mpg)))
    (t/is (v/delta-eq [10.4 33.9 20.090625] (sut/extent mpg)))
    (t/is (= [10.4 33.9] (sut/extent mpg false)) "mean? = false")
    (t/testing "adjacent-values clips to nearest in-fence data value, distinct from the fence thresholds"
      ;; mtcars$mpg has no outliers beyond the inner fence, so adjacent-values == min/max there
      ;; (already checked above); this dataset has a real outlier (50) to exercise the clip
      (t/is (v/delta-eq [10.4 33.9 19.2] (sut/adjacent-values mpg))
            "no outliers here: LAV/UAV == min/max")
      (t/is (= [2.0 8.0 5.0] (sut/adjacent-values [2.0 3 4 4 5 5 5 6 6 7 8 50]))
            "outlier (50) correctly excluded from UAV"))))

;; outliers
;; same inner-fence definition and reference values as adjacent-values, above (Extent Range group)

(t/deftest outliers-test
  (let [xout [2.0 3 4 4 5 5 5 6 6 7 8 50]]
    (t/is (= '(50.0) (sut/outliers xout)) "fence [-0.125, 10.875] (type=6/:legacy Q1,Q3)")
    (t/is (= '(2.0 3.0 4.0 4.0 5.0 5.0 5.0 6.0 6.0 7.0 8.0) (sut/remove-outliers xout)))
    (t/is (= (count xout) (+ (count (sut/outliers xout)) (count (sut/remove-outliers xout))))
          "outliers and remove-outliers partition vs"))
  (t/testing "mtcars$mpg has no outliers beyond its inner fence (see Extent Range group)"
    (t/is (empty? (sut/outliers (mtcars :mpg))))
    (t/is (= 32 (count (sut/remove-outliers (mtcars :mpg)))))))

;; ci, stats-map, bootstrap, bootstrap-ci
;; reference values from R on mtcars$mpg: ci <-> t.test(mpg)$conf.int (mean +/- qt(1-alpha/2,
;; n-1)*sd/sqrt(n)); stats-map's 24 keys cross-checked against the already-R-verified
;; primitives from earlier groups (Core Extremes & Sum, Percentiles & Quantiles, Variance &
;; Dispersion, Extent Range Family, Outliers, Mode), plus fresh e1071::skewness/kurtosis
;; (type=2) and DescTools::Mode values; bootstrap (deprecated, stochastic) checked
;; structurally, not by value; bootstrap-ci (deprecated, stochastic) checked for statistical
;; consistency against an independent R replication of the same basic/pivotal-bootstrap
;; procedure at N=20000 resamples -- exact agreement isn't expected (two independently-seeded
;; Monte Carlo runs), a tolerance of ~1% of the mean is used, and the deterministic `mean`
;; component (the third vector element) is checked exactly

(t/deftest ci-bootstrap-stats-map-test
  (t/testing "ci: Student's t-distribution CI for the mean"
    (t/is (v/delta-eq [17.917678508746246 22.263571491253753 20.090625] (sut/ci (mtcars :mpg))))
    (t/is (= (sut/ci (mtcars :mpg)) (sut/ci (mtcars :mpg) 0.05)) "default alpha == 0.05"))
  (t/testing "stats-map: all 24 keys, cross-checked against already-verified primitives + R"
    (let [m (sut/stats-map (mtcars :mpg))]
      (t/are [k v] (m/delta= v (k m))
        :Size 32
        :Min 10.4
        :Max 33.9
        :Range 23.5
        :Mean 20.090625
        :Median 19.2
        :Mode 10.4
        :Q1 15.275
        :Q3 22.8
        :Total 642.9
        :SD 6.026948052089105
        :Variance 36.32410282258065
        :MAD 3.65
        :SEM 1.0654239593728148
        :LAV 10.4
        :UAV 33.9
        :IQR 7.525
        :LOF -7.3
        :UOF 45.375
        :LIF 3.9875
        :UIF 34.0875
        :Kurtosis -0.0220062914
        :Skewness 0.6723771376)
      (t/is (empty? (:Outliers m)))))
  (t/testing "bootstrap (deprecated): structural properties, not exact values (stochastic)"
    (let [mpg-set-d (set (map double (mtcars :mpg)))
          samples (sut/bootstrap (mtcars :mpg) 50)]
      (t/is (= 50 (count samples)))
      (t/is (every? #(= 32 (count %)) samples) "default size == input count")
      (t/is (every? (fn [s] (every? mpg-set-d s)) samples)
            "every resampled value is drawn from the original data")
      (t/is (every? #(= 5 (count %)) (sut/bootstrap (mtcars :mpg) 10 5)) "explicit size")))
  (t/testing "bootstrap-ci (deprecated): statistical consistency vs. an independent R
              replication of the same basic-bootstrap procedure at N=20000, not exact match"
    (let [[lo hi mn] (sut/bootstrap-ci (mtcars :mpg) 0.98 20000)]
      (t/is (m/delta= 20.090625 mn) "the deterministic mean component matches exactly")
      (t/is (m/delta-eq 17.94375 lo 0.2) "R replication: 17.94375, generous MC tolerance")
      (t/is (m/delta-eq 22.1875 hi 0.2) "R replication: 22.1875, generous MC tolerance"))))

;; standardize/robust-standardize/demean/rescale
;; reference values from R on mtcars$mpg: standardize <-> scale(x); robust-standardize <->
;; (x-median)/mad(x,constant=1), and with q: (x-median)/abs(quantile(1-q)-quantile(q)),
;; type=6; demean <-> x-mean(x); rescale <-> linear map from [min,max] to [low,high]

(t/deftest scaling-standardization-test
  (let [mpg (mtcars :mpg)
        first5 (fn [s] (take 5 s))]
    (t/is (v/delta-eq [0.15088482 0.15088482 0.44954345 0.21725341 -0.23073453] (first5 (sut/standardize mpg))))
    (t/is (m/delta-eq 0.0 (sut/sum (sut/standardize mpg)) 1.0e-10))
    (t/is (m/delta-eq 1.0 (sut/stddev (sut/standardize mpg))))
    (t/is (v/delta-eq [0.49315068 0.49315068 0.98630137 0.60273973 -0.13698630] (first5 (sut/robust-standardize mpg))))
    (t/is (v/delta-eq [0.23920266 0.23920266 0.47840532 0.29235880 -0.06644518] (first5 (sut/robust-standardize mpg 0.25))))
    (t/is (v/delta-eq [0.90937500 0.90937500 2.70937500 1.30937500 -1.39062500] (first5 (sut/demean mpg))))
    (t/is (m/delta-eq 0.0 (sut/sum (sut/demean mpg)) 1.0e-10))
    (t/is (v/delta-eq [0.45106383 0.45106383 0.52765957 0.46808511 0.35319149] (first5 (sut/rescale mpg))))
    (t/is (v/delta-eq [3.53191489 3.53191489 5.82978723 4.04255319 0.59574468] (first5 (sut/rescale mpg -10.0 20.0))))))

;; covariance/correlation family, incl. matrix forms
;; reference values from R on mtcars$mpg vs $cyl (pairwise) and mpg/cyl/hp (matrix forms):
;; cov(), cor(method="pearson"/"spearman"/"kendall"); kendall confirmed to match R's default
;; tau-b convention (ties present in cyl)

(t/deftest correlation-covariance-test
  (let [mpg (mtcars :mpg) cyl (mtcars :cyl) hp (mtcars :hp)]
    (t/is (m/delta= -9.1723790323 (sut/covariance mpg cyl)))
    (t/is (= (sut/covariance mpg cyl) (sut/covariance [mpg cyl])) "2-seq-wrapped arity")
    (t/is (m/delta= -0.8521619594 (sut/correlation mpg cyl)))
    (t/is (m/delta= -0.8521619594 (sut/pearson-correlation mpg cyl)))
    (t/is (m/delta= -0.9108013109 (sut/spearman-correlation mpg cyl)))
    (t/is (m/delta= -0.7953134086 (sut/kendall-correlation mpg cyl)) "R's default tau-b, ties present in cyl")
    (t/testing "matrix forms (mpg, cyl, hp)"
      (let [vss [mpg cyl hp]]
        (t/is (v/delta-eq [36.32410282 -9.17237903 -320.73205645] (first (sut/covariance-matrix vss))))
        (t/is (v/delta-eq [-9.17237903 3.18951613 101.93145161] (second (sut/covariance-matrix vss))))
        (t/is (v/delta-eq [1.0 -0.85216196 -0.77616837] (first (sut/correlation-matrix vss))))
        (t/is (v/delta-eq [-0.85216196 1.0 0.83244745] (second (sut/correlation-matrix vss))))
        (t/is (v/delta-eq [1.0 -0.79531341 -0.74281251] (first (sut/correlation-matrix vss :kendall))))
        (t/is (every? true? (map v/delta-eq (sut/correlation-matrix vss) (sut/coefficient-matrix vss)))
              "coefficient-matrix defaults to pearson-correlation (symmetric?=true caches/mirrors one
               computed value per pair, symmetric?=false computes each direction independently --
               both correct, differing only in the last ULP for a non-perfectly-symmetric measure-fn)")))))

;; kullback-leibler-divergence / jensen-shannon-divergence (deprecated)
;; reference values: manual R computation of sum(p*log(p/q)) (p=0 terms skipped by convention)

(t/deftest divergence-deprecated-test
  (let [p1 [0.2 0.3 0.5] q1 [0.1 0.4 0.5]
        p2 [0.0 0.5 0.5] q2 [0.2 0.3 0.5]]
    (t/is (m/delta-eq 0.0523248144 (sut/kullback-leibler-divergence p1 q1)))
    (t/is (= (sut/kullback-leibler-divergence p1 q1) (sut/kullback-leibler-divergence [p1 q1]))
          "2-seq-wrapped arity")
    (t/is (m/delta-eq 0.2554128119 (sut/kullback-leibler-divergence p2 q2)) "p=0 term skipped by convention")
    (t/is (m/delta-eq 0.0120786284 (sut/jensen-shannon-divergence p1 q1)))
    ;; bug fixed this session: when q=0 but p>0, the term was silently dropped instead of
    ;; making the whole divergence +Inf, producing a mathematically impossible negative KL
    ;; divergence (KL >= 0 always, by Gibbs' inequality) -- confirmed pre-fix: -0.3466 for
    ;; p=[0.5 0.5] q=[0.0 1.0]; post-fix: correctly +Inf
    (t/is (m/pos-inf? (sut/kullback-leibler-divergence [0.5 0.5] [0.0 1.0]))
          "q=0 with p>0 forces +Inf, per the standard KL-divergence definition")
    (t/is (m/delta= 0.0 (sut/kullback-leibler-divergence [0.0 0.0] [0.5 0.5]))
          "degenerate all-p-zero case is 0.0 (was the incorrect +Inf fallback pre-fix)")))

;; regression error metrics: me/mae/mape/rss/r2/mse/rmse/count=/L1/L2sq/L2/LInf/psnr
;; reference values: manual R computation of each metric's standard formula

(t/deftest regression-error-metrics-test
  (let [v1 [3.0 5 2 8 10 1 7]
        v2 [2.5 5.5 1.0 9.0 8.0 1.5 7.5]
        w [1.0 2 0.5 1.5 1 3 2]]
    (t/is (m/delta= 0.1428571429 (sut/me v1 v2)))
    (t/is (= (sut/me v1 v2) (sut/me [v1 v2])) "2-seq-wrapped arity")
    (t/is (m/delta= 0.8571428571 (sut/mae v1 v2)))
    (t/is (m/delta= 0.2375850340 (sut/mape v1 v2)))
    (t/is (m/delta= 7.0 (sut/rss v1 v2)))
    (t/is (m/delta= 1.0 (sut/mse v1 v2)))
    (t/is (m/delta= 1.0 (sut/rmse v1 v2)))
    (t/is (m/delta= 0.8952991453 (sut/r2 v1 v2)))
    (t/is (m/delta= 0.8429487179 (sut/r2 v1 v2 2.0)) "adjusted R^2, p=2")
    (t/is (= 0 (sut/count= v1 v2)))
    (t/is (= 2 (sut/count= [1.0 2.0 3.0] [1.0 5.0 3.0])))
    (t/is (m/delta= 6.0 (sut/L1 v1 v2)))
    (t/is (m/delta= 7.0 (sut/L2sq v1 v2)))
    (t/is (m/delta= 2.6457513111 (sut/L2 v1 v2)))
    (t/is (m/delta= 2.0 (sut/LInf v1 v2)))
    (t/testing "weighted arities"
      ;; bug fixed this session: mse's weighted arity called (sum (weights)) -- invoking
      ;; `weights` (a vector) as a 0-argument function -- instead of (sum weights), throwing
      ;; an ArityException on every call; also affects rmse's weighted arity, built on mse
      (t/is (m/delta= -0.1818181818 (sut/me v1 v2 w)))
      (t/is (m/delta= 0.7272727273 (sut/mae v1 v2 w)))
      (t/is (m/delta= 0.7272727273 (sut/mse v1 v2 w)) "previously threw ArityException")
      (t/is (m/delta= 0.8528028654 (sut/rmse v1 v2 w)) "previously threw ArityException (via mse)")
      (t/is (m/delta= 0.1509433962 (sut/mape v1 v2 w))))
    (t/testing "scalar (vs2-or-val) broadcast arity"
      (t/is (m/delta= 0.1428571429 (sut/me v1 5.0)))
      (t/is (m/delta= 2.7142857143 (sut/mae v1 5.0))))
    (t/testing "psnr"
      (t/is (m/delta= 20.0 (sut/psnr v1 v2)) "auto max-value == max(v1,v2) == 10.0")
      (t/is (m/delta= 48.1308036087 (sut/psnr v1 v2 255.0))))))

;; effect size
;; reference values: R, group1=mtcars$mpg[am==0] (n=19), group2=mtcars$mpg[am==1] (n=13).
;; cohens-d/-corrected: manual pooled-sd formulas, `:unbiased` cross-checked vs the
;; pooled-variance group's ANOVA-MSE identity. hedges-g/-corrected/*: same formulas plus
;; the exact J = exp(lgamma(df/2) - log(sqrt(df/2)) - lgamma((df-1)/2)) correction.
;; glass-delta: (mean1-mean2)/sd(group2). means-ratio(-corrected): mean1/mean2, Bickel-Doksum
;; J-adjusted variant. cliffs-delta/ameasure/wmw-odds: manual sign-sum / rank-sum formulas.
;; cohens-u3: proportion of group2 below median(group1).
;; cohens-u2/-u1: see the dedicated cohens-u2-global-minimum-test below (bug fixed this
;; session: Brent's method converged to a non-global local minimum of this non-unimodal,
;; piecewise-constant objective; replaced with an exact finite-breakpoint search).

(def mpgs (by mtcars :am :mpg))

(t/deftest effect-size-test
  (t/are [res f attr] (m/delta= res (apply f (mpgs 0) (mpgs 1) attr))
    -1.477947 sut/cohens-d nil
    -1.477947 sut/cohens-d [:unbiased]
    -1.526417 sut/cohens-d [:biased]
    -1.411046 sut/cohens-d [:avg]
    -1.440688 sut/cohens-d-corrected nil
    -1.440688 sut/cohens-d-corrected [:unbiased]
    -1.490360 sut/cohens-d-corrected [:biased]
    -1.4779470958015888 sut/hedges-g nil
    -1.4406879253191958 sut/hedges-g-corrected nil
    -1.4406354024018648 sut/hedges-g* nil
    -1.174886 sut/glass-delta nil
    0.7029826 sut/means-ratio nil
    0.7021799 sut/means-ratio-corrected nil
    -0.659919 sut/cliffs-delta nil
    0.1700405 sut/ameasure nil
    0.204878 sut/wmw-odds nil
    0.1538462 sut/cohens-u3 nil))

;; cohens-u1-normal/u2-normal/u3-normal: R pnorm() applied to fastmath's own cohens-d
;; (-1.4779470958015888, full double precision, :unbiased); d=-1.4779470958015888 =>
;; p=pnorm(0.5*abs(d))=0.77003847002150, u1=(2p-1)/p, u2=p, u3=pnorm(d).
(t/deftest cohens-u-normal-test
  (let [g1 (mpgs 0) g2 (mpgs 1)]
    (t/is (m/delta-eq 0.70136358255960 (sut/cohens-u1-normal g1 g2)))
    (t/is (m/delta-eq 0.77003847002150 (sut/cohens-u2-normal g1 g2)))
    (t/is (m/delta-eq 0.06971096959831 (sut/cohens-u3-normal g1 g2)))))

;; cohens-u2/cohens-u1: bug fixed this session (blame: implementation bug). cohens-u2's
;; objective -- min(|icdf1(p)-icdf2(1-p)|, |icdf1(1-p)-icdf2(p)|) over p in [0.5,1], where
;; icdf1/icdf2 are empirical (real-discrete-distribution) step-function quantiles -- is
;; piecewise-constant and generally non-unimodal, so the old `opt/minimize :brent` search
;; could (and did, on the mtcars case below) converge to a non-global local minimum.
;; Fixed by replacing it with an exact search over the finite set of breakpoints where the
;; objective can change value (each sample's empirical cumulative-probability breakpoints
;; and their 1-p reflections): each resulting shelf is sampled at its *midpoint* (never
;; exactly on a jump boundary, avoiding floating-point ambiguity about which side of a jump
;; a boundary value falls on) to find the minimal shelf, but its *left edge* is reported
;; (so two identical/fully-overlapping samples, whose minimal shelf always starts exactly at
;; p=0.5, canonically report u2=0.5). Every edge is *also* evaluated directly: two of the
;; objective's four icdf lookups are taken at `1-p`, which is non-increasing in p, so whenever
;; a direct breakpoint of one sample coincides with 1 minus a breakpoint of either sample, the
;; objective has an isolated single point (invisible to shelf-midpoint sampling) that can be
;; the true global minimum -- see the "isolated single-point minimum" case below. Ties are
;; broken toward the smallest p. All reference values below are cross-checked against a
;; brute-force grid search (100k-5M points) over [0.5,1] confirming the reported p attains the
;; true global minimum of the objective -- not against R, whose default (continuous,
;; interpolated) quantile convention differs from fastmath's step-function empirical icdf, so
;; an exact numeric match to R is not expected (this was the pre-existing disabled assertion's
;; root cause, now correctly characterized rather than just flagged).
(t/deftest cohens-u2-global-minimum-test
  (t/testing "mtcars case: old Brent-based value (0.7977457514062631) was NOT the true minimum"
    (let [g1 (mpgs 0) g2 (mpgs 1)]
      (t/is (m/delta-eq 0.7692307692307692 (sut/cohens-u2 g1 g2)))
      (t/is (m/delta-eq 0.7692307692307692 (sut/cohens-u2 g2 g1)) "symmetric")
      (t/is (m/delta-eq 0.6999999999999998 (sut/cohens-u1 g1 g2)))))
  (t/testing "isolated single-point minimum: group1's breakpoint 11/16=0.6875 exactly equals
1 minus group2's breakpoint 5/16=0.3125, so p=0.6875 alone attains distance 0 while both
0.68749999 and 0.68750001 attain distance 1 -- no shelf-midpoint sees this point"
    (let [g1 [9.0 7.0 5.0 5.0 6.0 5.0 6.0 6.0 10.0 5.0 6.0 9.0 9.0 1.0 2.0 6.0]
          g2 [8.0 7.0 6.0 9.0 7.0 9.0 8.0 7.0 7.0 9.0 2.0 10.0 1.0 3.0 9.0 1.0]]
      (t/is (m/delta-eq 0.6875 (sut/cohens-u2 g1 g2)))
      (t/is (m/delta-eq 0.6875 (sut/cohens-u2 g2 g1)) "symmetric")))
  (t/testing "identical groups: complete overlap, minimal distance is 0 on a shelf starting at 0.5"
    (t/is (m/delta-eq 0.5 (sut/cohens-u2 [1.0 2.0 3.0] [1.0 2.0 3.0])))
    (t/is (m/delta-eq 0.5 (sut/cohens-u2 [1.0 2.0 3.0 4.0 5.0] [1.0 2.0 3.0 4.0 5.0])))
    (t/is (m/delta-eq 0.0 (sut/cohens-u1 [1.0 2.0 3.0] [1.0 2.0 3.0]))))
  (t/testing "fully separated groups: minimal shelf is not the first one tried"
    (t/is (m/delta-eq 0.6666666666666666 (sut/cohens-u2 [1.0 2.0 3.0] [10.0 11.0 12.0])))
    (t/is (m/delta-eq 0.4999999999999999 (sut/cohens-u1 [1.0 2.0 3.0] [10.0 11.0 12.0])))))

;; p-overlap: group1=[0..4], group2=[2..6], explicit gaussian bandwidth h=1.0 (fixed, to make
;; the reference reproducible independent of fastmath's own bandwidth-estimation rule).
;; Reference: independent R re-implementation of the exact same Gaussian-KDE formula
;; (f(x) = mean(dnorm((x-data)/h))/h), pointwise min of f1,f2, Simpson-integrated over the
;; matching widened domain -> 0.60964001. fastmath's own integration approximates
;; integral(min(f1,f2)) by summing, over `:steps` partition pieces, min(integral(f1 over
;; piece), integral(f2 over piece)) -- exact only in the limit of fine partitioning (an
;; approximation of the true value, not a bug); verified by observing monotonic convergence
;; toward the R reference as `:steps` increases (500 -> 0.60965834, 20000 -> 0.60964002).
(t/deftest p-overlap-test
  (t/is (m/delta-eq 0.60964001 (sut/p-overlap [0.0 1.0 2.0 3.0 4.0] [2.0 3.0 4.0 5.0 6.0]
                                               {:bandwidth 1.0 :steps 20000}))))

;; effect size - correlations

(t/deftest effect-size-correlation-test
  (t/are [res f] (m/delta= res (f (mtcars :cyl) (mtcars :mpg)))
    -0.852162 sut/pearson-r
    0.72618 sut/r2-determination
    0.72618 sut/eta-sq
    0.7105671 sut/omega-sq
    0.7170527 sut/epsilon-sq)
  (t/are [method res] (m/delta= res (sut/cohens-f2 (mtcars :cyl) (mtcars :mpg) method))
    :eta 2.652034
    :omega 2.455032
    :epsilon 2.534227)
  (t/are [method res] (m/delta= res (sut/cohens-f (mtcars :cyl) (mtcars :mpg) method))
    :eta 1.628507
    :omega 1.566854
    :epsilon 1.5919255)
  (t/is (m/delta= 2.0921305 (sut/cohens-q (sut/pearson-correlation
                                           (iris :petal-width) (iris :petal-length))
                                          (sut/pearson-correlation
                                           (iris :sepal-width) (iris :sepal-length)))))
  (t/is (m/delta= 2.0921305 (sut/cohens-q (iris :petal-width) (iris :petal-length)
                                          (iris :sepal-width) (iris :sepal-length))))
  (t/is (m/delta= -1.9568811 (sut/cohens-q (mtcars :mpg)
                                           (mtcars :cyl) (mtcars :am)))))

;; contingency-table builders and marginals
;; reference values are manually-counted expected frequencies/sums (pure data-munging functions,
;; independently hand-counted from the fixtures rather than via an external statistics package)

(t/deftest contingency-table-builders-test
  (t/is (= {:a 3 :b 2 :c 1} (sut/contingency-table [:a :b :a :c :a :b])))
  (t/is (= {[:x 1] 3, [:y 2] 2, [:y 1] 1, [:x 2] 1}
           (sut/contingency-table [:x :x :y :y :x :y :x] [1 1 2 1 2 2 1])))
  (let [rows [[10 5 0] [3 12 7]]]
    (t/is (= {[0 0] 10 [0 1] 5 [1 0] 3 [1 1] 12 [1 2] 7}
             (sut/rows->contingency-table rows)))
    (let [{:keys [rows cols n diag]} (sut/contingency-table->marginals rows)]
      (t/is (= [[0 15.0] [1 22.0]] rows))
      (t/is (= [[0 13.0] [1 17.0] [2 7.0]] cols))
      (t/is (m/delta= 37.0 n))
      (t/is (= [[[0 0] 10] [[1 1] 12]] diag)))))

;; mcc (Matthews Correlation Coefficient)
;; reference: 2x2 classic formula (a*d-b*c)/sqrt((a+b)(a+c)(b+d)(c+d)); RxC multiclass generalization
;; cross-checked against Python sklearn.metrics.matthews_corrcoef on a 3-class confusion matrix,
;; verified both via a pre-built confusion matrix and via the raw-label 2-sequence arity

(t/deftest mcc-test
  (t/is (m/delta= 0.6847367880 (sut/mcc [[30 10] [5 55]])))
  (let [cm [[3 2 0] [0 2 2] [1 0 4]]
        y-true [0 0 0 1 1 1 2 2 2 2 0 1 2 0]
        y-pred [0 0 1 1 1 2 2 2 0 2 0 2 2 1]]
    (t/is (m/delta= 0.4651302547 (sut/mcc cm)))
    (t/is (m/delta= 0.4651302547 (sut/mcc y-true y-pred)))))

;; nominal association (chi-squared based)
;; reference values from R: DescTools::{CramerV,ContCoef,TschuprowT}(mtcars$cyl, mtcars$am),
;; cross-checked against Python's dython.nominal.cramers_v(bias_correction=False) (dython's default
;; bias_correction=True additionally applies scipy's Yates continuity correction for 2x2 tables,
;; which is not applicable here since cyl has 3 levels)

(t/deftest nominal-association-test
  (t/are [f res] (m/delta= res (f (mtcars :cyl) (mtcars :am)))
    sut/cramers-v           0.5226355372
    sut/cramers-v-corrected 0.4643125760
    sut/cramers-c           0.4631903544
    sut/tschuprows-t        0.4394823497
    sut/cohens-w            0.5226355372))

;; power divergence tests, incl. Yates' continuity correction for 2x2 tables
;; reference values from R: chisq.test(matrix(c(70,4,2,40), nrow=2, byrow=TRUE), correct=TRUE/FALSE)
;; and Python: scipy.stats.chi2_contingency(tab, correction=TRUE/FALSE, lambda_='log-likelihood')

(def t2x2 [[70 2] [4 40]])

(t/deftest power-divergence-test-yates
  (let [res (sut/chisq-test t2x2)]
    (t/is (m/delta= 91.8380458380 (:stat res)))
    (t/is (= 1 (:df res)))
    (t/is (m/delta= 88.0620571871 (:yates res)) "Yates-corrected statistic, matches R's chisq.test(correct=TRUE)")
    (t/is (contains? res :yates-p-value)))
  (let [res (sut/multinomial-likelihood-ratio-test t2x2)]
    (t/is (m/delta= 106.7810675108 (:stat res)))
    (t/is (m/delta= 101.1087556019 (:yates res)) "Yates correction generalizes to other lambda values (here: G-test)"))
  (t/testing "Yates' correction only applies to 2x2 tables (df = 1)"
    (let [res (sut/chisq-test [[3 8] [4 3] [12 2]])]
      (t/is (= 2 (:df res)))
      (t/is (not (contains? res :yates)))
      (t/is (not (contains? res :yates-p-value))))))

;; remaining power-divergence lambda variants (lambda=-2,-1,-0.5,2/3), incl. the base
;; power-divergence-test fn called directly (goodness-of-fit and independence dispatch branches)
;; reference: independent Cressie-Read formula reimplementation in R:
;;   lambda=0:  2*sum(obs*log(obs/exp))
;;   lambda=-1: 2*sum(exp*log(exp/obs))
;;   else:      (2/(lambda*(lambda+1))) * sum(obs*((obs/exp)^lambda - 1))
;; applied manually to both a goodness-of-fit fixture (observed counts vs custom :p probabilities)
;; and an independence (3x3 contingency table) fixture; also spot-checks Yates' correction
;; generalizes correctly to lambda=-1 on a 2x2 table.

(t/deftest power-divergence-remaining-lambdas-test
  (let [obs [10 20 30 25 15]
        p [0.15 0.25 0.30 0.20 0.10]]
    (t/are [f stat] (let [res (f obs {:p p})]
                      (and (m/delta= stat (:stat res)) (= 4 (:df res))))
      sut/minimum-discrimination-information-test 6.2860865942
      sut/neyman-modified-chisq-test              6.4166666667
      sut/freeman-tukey-test                      6.2699441774
      sut/cressie-read-test                       6.3583136471)
    (t/is (m/delta= 6.2860865942 (:stat (sut/power-divergence-test obs {:p p :lambda -1.0}))))
    (t/is (= 4 (:df (sut/power-divergence-test obs {:p p :lambda -1.0})))
          "power-divergence-test's own goodness-of-fit dispatch branch"))
  (let [tab [[20 15 30] [10 25 5] [8 12 22]]]
    (t/are [f stat] (let [res (f tab)]
                      (and (m/delta= stat (:stat res)) (= 4 (:df res)) (m/delta= 147.0 (:n res))))
      sut/minimum-discrimination-information-test 27.8910368298
      sut/neyman-modified-chisq-test              35.5357903623
      sut/freeman-tukey-test                      25.7249037347
      sut/cressie-read-test                       23.2659726491)
    (t/is (m/delta= 27.8910368298 (:stat (sut/power-divergence-test tab {:lambda -1.0})))
          "power-divergence-test's own independence dispatch branch"))
  (t/testing "Yates' correction generalizes to lambda=-1 (minimum-discrimination-information-test)"
    (let [res (sut/minimum-discrimination-information-test [[18 7] [5 20]])]
      (t/is (m/delta= 15.9732388129 (:stat res)))
      (t/is (m/delta= 13.2498081331 (:yates res))))))

;; power-divergence-test's third dispatch branch: raw data + a distribution object + :bins
;; (goodness-of-fit against a distribution's binned probability mass, via the private
;; quantize-distribution helper). No independent R replication of the binning/RNG state is
;; feasible here, so this is a composition check: the same expected counts are recomputed from
;; scratch using fastmath's own already-independently-verified `histogram` (group 16) and
;; `fastmath.random/cdf` (distribution namespace's own test suite), fed through the same
;; Cressie-Read formula verified above, and compared against power-divergence-test's own result.

(t/deftest power-divergence-distribution-gof-test
  (let [distr (r/distribution :normal {:mu 0.0 :sd 1.0})
        data (vec (r/->seq distr 200))
        res (sut/power-divergence-test data {:p distr :bins 5 :lambda -1.0})
        {:keys [step bins]} (sut/histogram data 5)
        last-idx (dec (count bins))
        counts (map second bins)
        n (double (reduce + counts))
        probs (map-indexed (fn [id [^double s]]
                             (cond
                               (= id 0) (r/cdf distr (+ s step))
                               (= id last-idx) (- 1.0 (r/cdf distr s))
                               :else (- (r/cdf distr (+ s step)) (r/cdf distr s))))
                           bins)
        expected (map #(* n %) probs)
        manual-stat (* 2.0 (reduce + (map (fn [^double e ^double o] (* e (- (Math/log e) (Math/log o))))
                                          expected counts)))]
    (t/is (m/delta= manual-stat (:stat res)))
    (t/is (= 4 (:df res)))
    (t/is (m/delta= 200.0 (:n res)))))

;; pairwise-regression correlation-based effect sizes
;; reference values from R: lm(g1 ~ g2); manual SSreg/SStot/MSE formulas for omega2/epsilon2;
;; cor(g1,g2) for pearson-r/r2-determination

(t/deftest effect-size-pairwise-regression-test
  (let [g1 [2.3 4.5 3.1 5.6 6.2 4.8 3.9 5.1 4.4 6.0]
        g2 [1.1 2.5 1.8 3.2 3.9 2.7 2.1 3.0 2.4 3.5]]
    (t/is (m/delta= 0.9914320117 (sut/pearson-r g1 g2)))
    (t/is (m/delta= 0.9829374338 (sut/r2-determination g1 g2)))
    (t/is (m/delta= 0.9829374338 (sut/eta-sq g1 g2)) "eta-sq equals r2-determination by construction")
    (t/is (m/delta= 0.9787171847 (sut/omega-sq g1 g2)))
    (t/is (m/delta= 0.9808046130 (sut/epsilon-sq g1 g2)))
    (t/are [type f2 f] (and (m/delta= f2 (sut/cohens-f2 g1 g2 type))
                           (m/delta= f (sut/cohens-f g1 g2 type)))
      :eta     57.6078312108 7.5899822932
      :omega   45.9862649686 6.7813173476
      :epsilon 51.0958499652 7.1481361183)))

;; cohens-q: difference of Fisher z-transformed correlations
;; reference: R atanh(r1)-atanh(r2); 3-/4-arity checked structurally against fastmath's own
;; already-verified pearson-correlation composed with atanh (R's RNG can't be reproduced in Clojure,
;; so the composition itself, not a specific dataset, is the thing under test)

(t/deftest cohens-q-test
  (t/is (m/delta= 0.2397865401 (sut/cohens-q 0.5 0.3)))
  (let [gg1 [1.2 3.4 2.1 5.6 4.4 3.3 6.1 2.8 4.9 3.7 5.2 1.9 4.1 3.0 5.8]
        gg2a [2.1 4.0 1.8 5.2 4.9 2.7 6.5 3.1 4.2 3.9 5.5 2.2 3.8 2.9 6.0]
        gg2b [6.0 1.1 5.5 2.0 1.9 4.8 0.9 5.1 2.5 3.6 1.2 5.9 3.0 4.5 0.8]
        gg1b [3.1 2.2 4.5 1.8 5.0 2.9 3.6 4.1 2.0 5.3 1.5 3.8 2.6 4.9 3.3]
        gg2bb [1.9 4.2 2.1 5.5 1.2 4.8 2.7 3.0 5.1 1.6 4.4 2.3 5.0 1.8 3.5]
        ra (sut/pearson-correlation gg1 gg2a)
        rb (sut/pearson-correlation gg1 gg2b)
        rd (sut/pearson-correlation gg1b gg2bb)]
    (t/is (m/delta= (m/- (m/atanh ra) (m/atanh rb)) (sut/cohens-q gg1 gg2a gg2b)))
    (t/is (m/delta= (m/- (m/atanh ra) (m/atanh rd)) (sut/cohens-q gg1 gg2a gg1b gg2bb)))))

;; kruskal effect size

(t/deftest effect-size-kruskal
  (t/is (m/delta= 0.8305211 (sut/rank-epsilon-sq (by mtcars :cyl :mpg))))
  (t/is (m/delta= 0.818833 (sut/rank-eta-sq (by mtcars :cyl :mpg))))
  (t/testing "manual small fixture cross-checked against R kruskal.test(H) formulas"
    (let [xss [[1 2 3 4 5] [6 7 8 9] [2.5 10 11 12]]]
      (t/is (m/delta= 0.4398901099 (sut/rank-eta-sq xss)))
      (t/is (m/delta= 0.5332417582 (sut/rank-epsilon-sq xss))))))

;; one-way ANOVA (correlation ratio) effect size
;; reference values from R: aov(mpg ~ factor(cyl), data = mtcars); effectsize::{eta,omega,epsilon}_squared, cohens_f

(t/deftest effect-size-anova
  (let [xss (by mtcars :cyl :mpg)]
    (t/is (m/delta= 0.7324600596 (sut/anova-eta-sq xss)))
    (t/is (m/delta= 0.7074821420 (sut/anova-omega-sq xss)))
    (t/is (m/delta= 0.7140090293 (sut/anova-epsilon-sq xss)))
    (t/are [type f2 f] (and (m/delta= f2 (sut/anova-cohens-f2 xss type))
                           (m/delta= f (sut/anova-cohens-f xss type)))
      :eta     2.7377596728 1.6546176818
      :omega   2.4185947035 1.5551831736
      :epsilon 2.4966138875 1.5800676845))
  (let [xss [[45 70 29 15 21] [40 20 30 42] [65 95 80 70 85 73]]]
    (t/is (m/delta= 0.7033195021 (sut/anova-eta-sq xss)))
    (t/is (m/delta= 0.6380968449 (sut/anova-omega-sq xss)))
    (t/is (m/delta= 0.6538727524 (sut/anova-epsilon-sq xss)))
    (t/is (m/delta= 2.3706293706 (sut/anova-cohens-f2 xss)))
    (t/is (m/delta= 1.5396848283 (sut/anova-cohens-f xss))))
  (t/testing "omega/epsilon are clamped at zero for a near-null effect"
    (let [xss [[5 5 5 6] [5 4 5 5] [6 5 5 4]]]
      (t/is (m/delta= 0.125 (sut/anova-eta-sq xss)))
      (t/is (m/delta= 0.0 (sut/anova-omega-sq xss)))
      (t/is (m/delta= 0.0 (sut/anova-epsilon-sq xss))))))

;; entropy, mutual information and Theil's U (uncertainty coefficient)
;; reference values from R: DescTools::UncertCoef(mtcars$cyl, mtcars$am) (fully populated 3x2 table);
;; for the sparse cyl x gear table (which has empty cells) reference values were computed by summing
;; -sum(p*log(p)) over the observed (non-zero) cells only, matching this implementation, since
;; DescTools's own zero-cell correction perturbs the result by a tiny, cell-count-dependent amount.

(t/deftest entropy-and-theils-u-test
  (t/is (m/delta= 1.0612039760 (sut/entropy (mtcars :cyl))))
  (t/is (m/delta= 0.6754645825 (sut/entropy (mtcars :am))))
  (t/is (m/delta= 1.5914372257 (sut/joint-entropy (mtcars :cyl) (mtcars :am))))
  (t/is (m/delta= 0.1452313328 (sut/mutual-information (mtcars :cyl) (mtcars :am))))
  (t/is (m/delta= 1.5309937135 (sut/entropy (mtcars :cyl) 2.0)) "entropy in bits")
  (t/is (m/delta= 0.6931471806 (sut/entropy {:a 0.5 :b 0.5})) "accepts a probability/frequency map")
  (t/is (m/delta= 0.1672527922 (sut/theils-u (mtcars :cyl) (mtcars :am))) "defaults to :symmetric")
  (t/are [dir res] (m/delta= res (sut/theils-u (mtcars :cyl) (mtcars :am) dir))
    :symmetric 0.1672527922
    :group1    0.1368552475
    :group2    0.2150095453)
  (t/testing "sparse contingency table (some cyl x gear combinations do not occur)"
    (t/are [dir res] (m/delta= res (sut/theils-u (mtcars :cyl) (mtcars :gear) dir))
      :symmetric 0.3504371541
      :group1    0.3424817995
      :group2    0.3587708805))
  (t/testing "a constant sequence has zero entropy, and normalizing by it gives NaN"
    (t/is (m/delta= 0.0 (sut/entropy (repeat 10 :a))))
    (t/is (m/nan? (sut/theils-u (repeat 10 :a) (mtcars :am) :group1)))))

;; 2x2 contingency

(def c2x2 (sut/contingency-2x2-measures-all 70 2 4 40))
(def c2x2p (sut/contingency-2x2-measures-all 70 40 40 40))
(def c2x2rr (sut/contingency-2x2-measures-all 15 135 100 150))
(def c2x2ri (sut/contingency-2x2-measures-all 75 75 100 150))

(defn seq-delta-eq
  ([a b] (seq-delta-eq a b 1.0e-6))
  ([a b ^double acc] (every? identity (map #(m/delta-eq %1 %2 acc) a b))))

(defn map-delta-eq
  ([keys a b] (map-delta-eq keys a b 1.0e-6))
  ([keys a b acc]
   (let [f (apply juxt keys)]
     (seq-delta-eq (f a) (f b) acc))))

;; source: wiki, https://statpages.info/ctab2x2.html, https://arxiv.org/pdf/2203.09628.pdf --
;; re-verified 2026-09-21 by independently reproducing every :measures/:p-values formula from
;; standard textbook definitions in R (not by reading fastmath's implementation) and comparing;
;; :p-values tightened from the pre-existing 1e-3 literals to exact R pchisq() values.
(t/deftest contingency-2x2-measures-test
  (t/is (map-delta-eq [:chi2 :yates :cochran-mantel-haenszel]
                      (:p-values c2x2p)
                      {:chi2 0.06015674717006436 :yates 0.08348067296735406
                       :cochran-mantel-haenszel 0.06083537556771568}))
  (t/are [ks vs] (seq-delta-eq ((juxt :a :b :c :d) (get-in c2x2 ks)) vs 1.0e-3)
    [:expected] [45.931 26.069 28.069 15.931]
    [:proportions :table] [0.603 0.017 0.034 0.345]
    [:proportions :rows]  [0.972 0.028 0.091 0.909]
    [:proportions :cols]  [0.946 0.048 0.054 0.952])
  (t/is (map-delta-eq [:row1 :row2 :col1 :col2]
                      (:marginals (:proportions c2x2))
                      {:row1 0.621 :row2 0.379 :col1 0.638 :col2 0.362} 1.0e-3))
  (t/is (= (:table c2x2) {:a 70 :b 2 :c 4 :d 40}))
  (t/is (= (:marginals c2x2) {:row1 72 :row2 44 :col1 74 :col2 42 :total 116}))
  (t/are [v k] (m/delta= v (c2x2 k))
    350.0 :OR
    (m/log 350.0) :lOR
    116 :n
    10.6944444 :RR
    (m/sqrt (+ (/ 70.0) (/ 2.0) (/ 4.0) (/ 40.0))) :SE)
  (t/are [v k] (m/delta= v (get-in c2x2 [:measures k]))
    0.889172 :cohens-kappa
    0.994302 :yules-q
    0.8897794 :phi
    0.8891367 :scotts-pi
    0.9589041 :F1
    0.9030371 :gwets-ac1
    0.8965517 :holley-guilfords-g
    0.6666667 :mcnemars-chi2
    0.987322 :TCC
    0.664735 :PCC
    0.940078 :PCC-adjusted
    0.9589041 :F1
    0.8985198 :yules-y
    0.9030371 :gwets-ac1
    0.8813131 :youdens-j
    0.803805 :huberts-gamma
    0.9057971 :bangdiwalas-b
    91.838046 :chi2
    88.062057 :yates
    91.046338 :cochran-mantel-haenszel
    (m/sqrt (/ 91.838046 116.0)) :cramers-v
    2.1941417568110757 :cohens-h)
  (t/are [v k] (m/delta= v (get-in c2x2rr [:risk k]))
    150 :ES
    250 :CS
    0.1 :EER
    0.4 :CER
    0.3 :ARR
    3.333333 :NNT
    0.25 :RR
    0.75 :RRR
    0.75 :PFu
    (- (/ 15.0 150.0)
       (/ 100.0 250.0)) :RD)
  (t/are [v k] (m/delta= v (get-in c2x2ri [:risk k]))
    150 :ES
    250 :CS
    0.5 :EER
    0.4 :CER
    0.1 :ARI
    10.0 :NNH
    1.25 :RR
    0.25 :RRI
    0.2 :AFe
    (- (/ 75.0 150.0)
       (/ 100.0 250.0)) :RD))

;; contingency-2x2-measures-all / contingency-2x2-measures: all 4 documented input formats
;; (4 args; flat 4-seq; nested rows [[a b] [c d]]; map) must agree.
(t/deftest contingency-2x2-input-formats-test
  (t/testing "contingency-2x2-measures-all: all input formats agree"
    (t/is (= c2x2 (sut/contingency-2x2-measures-all [70 2 4 40])))
    (t/is (= c2x2 (sut/contingency-2x2-measures-all {:a 70 :b 2 :c 4 :d 40})))
    ;; bug fixed this session (blame: implementation bug): the single-argument nested-rows
    ;; form documented in the docstring -- (contingency-2x2-measures-all [[a b] [c d]]) --
    ;; threw an ArityException; only the *separate*-arguments 2-arity form
    ;; (contingency-2x2-measures-all [a b] [c d]) actually worked. Fixed by disambiguating
    ;; the single-arg case (mirroring fastmath.stats.binary/infer-confusion-matrix's existing
    ;; "pair of pairs" vs "flat sequence" pattern) so both forms now agree.
    (t/is (= c2x2 (sut/contingency-2x2-measures-all [[70 2] [4 40]])) "previously threw ArityException")
    (t/is (= c2x2 (sut/contingency-2x2-measures-all [70 2] [4 40]))))
  (t/testing "contingency-2x2-measures: subset matches contingency-2x2-measures-all, all input formats agree"
    (let [expected (assoc (:measures c2x2) :OR (:OR c2x2))]
      (t/is (= expected (sut/contingency-2x2-measures 70 2 4 40)))
      (t/is (= expected (sut/contingency-2x2-measures [70 2 4 40])))
      (t/is (= expected (sut/contingency-2x2-measures {:a 70 :b 2 :c 4 :d 40})))
      (t/is (= expected (sut/contingency-2x2-measures [[70 2] [4 40]])) "previously threw ArityException")
      (t/is (= expected (sut/contingency-2x2-measures [70 2] [4 40]))))))

;; acf / pacf / acf-ci / pacf-ci: reference R base acf()/pacf() (type="correlation"), and a
;; manual R reproduction of the CI formulas -- flat white-noise band qnorm(0.975)/sqrt(n), and
;; acf-ci's :cis cumulative Bartlett band qnorm(0.975)/sqrt(n) * sqrt(2*cumsum(r^2)-1).
(def ts-data [2.0 4.0 3.0 6.0 5.0 8.0 7.0 9.0 10.0 8.0 7.0 9.0 11.0 10.0 12.0 13.0 11.0 14.0 13.0 15.0])

(t/deftest acf-pacf-test
  (t/testing "acf matches R acf(x, type=\"correlation\")"
    (t/is (seq-delta-eq (sut/acf ts-data 10)
                         [1.0 0.70425843 0.61510427 0.40121809 0.33545118 0.16385695
                          0.09243812 0.02667121 -0.06423699 -0.06101150 -0.09052816]
                         1.0e-6)))
  (t/testing "acf: sequence-of-lags arity picks the same values out of the full sequence"
    (t/is (= [1.0 0.6151042681738454 0.16385694796336 -0.09052816215162739]
             (sut/acf ts-data [0 2 5 10]))))
  (t/testing "acf: 1-arity form defaults to lags 0..(dec (count data))"
    (t/is (= 20 (count (sut/acf ts-data)))))
  (t/testing "pacf matches R pacf(x) (lag 0 is always 0.0, R's pacf omits it)"
    (t/is (seq-delta-eq (sut/pacf ts-data 10)
                         [0.0 0.704258429 0.236348394 -0.201822500 0.066873071 -0.143192059
                          -0.041601440 0.060101533 -0.164169657 0.104372101 -0.002817981]
                         1.0e-6)))
  (t/testing "acf-ci / pacf-ci: :ci matches the flat white-noise band qnorm(0.975)/sqrt(n)"
    (t/is (m/delta-eq 0.43826127028829076 (:ci (sut/acf-ci ts-data 10))))
    (t/is (m/delta-eq 0.43826127028829076 (:ci (sut/pacf-ci ts-data 10)))))
  (t/testing "acf-ci: :cis matches the cumulative Bartlett formula"
    (t/is (seq-delta-eq (:cis (sut/acf-ci ts-data 10))
                         [0.4382613 0.6185480 0.7265979 0.7679731 0.7956190 0.8020746
                          0.8041182 0.8042881 0.8052729 0.8061603 0.8081105]
                         1.0e-6))))

;; binomial-ci: reference R binom::binom.confint(x, n, methods=c(8 shared methods)) for 8 of the
;; 9 methods; :arcsine (not in that package) verified via its standard closed-form manual R
;; formula (sin(asin(sqrt(p)) +- z/(2*sqrt(n)))^2). x=23, n=50.
(t/deftest binomial-ci-test
  (let [res (sut/binomial-ci 23 50 :all)]
    (t/are [k lo hi] (and (m/delta-eq lo (first (res k)))
                           (m/delta-eq hi (second (res k))))
      :asymptotic     0.3218538187 0.5981461813
      :agresti-coull  0.3296681238 0.5960396842
      :clopper-pearson 0.3181491786 0.6067580239
      :wilson         0.3296965219 0.5960112860
      :prop.test      0.3206340998 0.6054718510
      :cloglog        0.3188005138 0.5900967075
      :logit          0.3281671048 0.5976784840
      :probit         0.3269085628 0.5977773762
      :arcsine        0.3251427495 0.5979107922)
    (t/is (every? #(m/delta-eq 0.46 (nth (res %) 2)) (keys res)) "estimated p is x/n for every method")))

;; percentile-bca-extent / percentile-bc-extent: reference is an independent R reimplementation
;; of the documented BCa/BC algorithm (Efron & Tibshirani), NOT a third-party R package (coxed,
;; the package the docstring cites for :r7 parity, is unavailable for this R version). BC (no
;; acceleration, accel=0.0) needs no external input and matches exactly. BCa additionally needs
;; an acceleration value derived from skewness; fastmath's own `:skew` estimator (a distinct,
;; already-verified BCa-specific skewness variant, see [[Skewness & Kurtosis]]) is taken as a
;; given input to isolate what this group is responsible for: the BCa mechanics themselves
;; (bias-correction z0, the accelerated-quantile formula, and the final quantile lookup).
(def bca-data [3.928902 -1.728267 -2.623711 -2.700109 -1.711872 2.502012 3.924331 -0.041672
               4.873639 2.502897 -11.008574 -7.343458 1.626829 9.413745 5.920651 -10.092515
               2.872223 -4.960665 -5.12469 14.326941 -6.075085 -6.681371 7.410476 12.560565
               3.556289 3.668704 -6.397729 1.74468 -9.703922 9.338795 -5.756261 -7.942128
               7.399011 5.800307 2.796319 -10.523582 3.610312 -2.200884 -2.671818 0.027861
               -7.317924 3.338304 -7.83735 -3.56692 6.369877 5.29273 -7.177113 7.996789
               4.410622 8.149898 11.373292 -7.1747 3.417676 12.670954 -8.840221 -3.873605
               3.541493 -0.302805 -0.282437 -8.060906])

(t/deftest percentile-bca-bc-extent-test
  (t/testing "percentile-bc-extent (accel=0.0, no skewness dependency): default (:legacy/type-6)
and :r7 (type-7) estimation strategies"
    (t/is (seq-delta-eq (sut/percentile-bc-extent bca-data)
                         [-10.753952938566199 13.45754787543103 0.2774138335000996] 1.0e-6))
    (t/is (seq-delta-eq (sut/percentile-bc-extent bca-data 2.5 97.5 :r7)
                         [-10.31882528877211 12.618518997218853 0.2774138335000996] 1.0e-6)))
  (t/testing "percentile-bca-extent: with fastmath's own :skew estimator as the (separately
verified) acceleration input"
    (t/is (seq-delta-eq (sut/percentile-bca-extent bca-data)
                         [-10.768583157145699 13.407193188901768 0.2774138335000996] 1.0e-6))
    (t/is (seq-delta-eq (sut/percentile-bca-extent bca-data 2.5 97.5 :r7)
                         [-10.331402455977047 12.615272367305515 0.2774138335000996] 1.0e-6))))

;; p-value: reference is manual R pnorm() (continuous, standard normal) and manual R ppois()
;; (discrete, Poisson lambda=5, exercising the k-1 continuity correction for CCDF-based sides).
(t/deftest p-value-test
  (t/testing "continuous (default normal distribution): two-sided/greater/less and their aliases"
    (t/are [res sides] (m/delta-eq res (sut/p-value r/default-normal 1.5 sides))
      0.1336144025 :two-sided
      0.1336144025 :both
      0.0668072013 :one-sided-greater
      0.0668072013 :right
      0.9331927987 :one-sided-less
      0.9331927987 :left
      0.9331927987 :one-sided)
    (t/is (m/delta-eq 0.1336144025 (sut/p-value 1.5)) "1-arity defaults to :two-sided vs default-normal")
    (t/is (m/delta-eq 0.1336144025 (sut/p-value r/default-normal 1.5)) "2-arity defaults to :two-sided"))
  (t/testing "discrete (Poisson, lambda=5): CCDF-based sides apply the k-1 continuity correction"
    (let [pois (r/distribution :poisson {:p 5.0})]
      (t/are [res sides] (m/delta-eq res (sut/p-value pois 8.0 sides))
        0.2667433481 :two-sided
        0.1333716741 :one-sided-greater
        0.9319063653 :one-sided-less))))

;; skewness-test / kurtosis-test / normality-test / jarque-bera-test / bonett-seier-test /
;; binomial-test. Reference: R packages moments (agostino.test, anscombe.test, jarque.test) and
;; tseries (jarque.bera.test) for the first/second/fourth; normality-test (K^2 omnibus)
;; cross-checked by composing the already-verified skew-Z/kurt-Z into K^2=Z1^2+Z2^2 and an
;; independent R pchisq(K2,2) call; bonett-seier-test's omega/Z/p-value transform verified via
;; a from-scratch R reimplementation of the formula, using fastmath's own (separately verified
;; in Skewness & Kurtosis) `:geary` kurtosis estimator as the given input; binomial-test vs R
;; base binom.test. n=50 dataset.
(def norm-xs [3.9 1.7 2.6 2.7 1.7 2.5 3.9 0.0 4.9 2.5 11.0 7.3 1.6 9.4 5.9 10.1 2.9 5.0 5.1 14.3
              6.1 6.7 7.4 12.6 3.6 3.7 6.4 1.7 9.7 9.3 5.8 7.9 7.4 5.8 2.8 10.5 3.6 2.2 2.7 0.0
              7.3 3.3 7.8 3.6 6.4 5.3 7.2 8.0 4.4 8.1])

(t/deftest skewness-kurtosis-normality-tests-test
  (t/testing "skewness-test vs R moments::agostino.test"
    (let [{:keys [Z p-value]} (sut/skewness-test norm-xs)]
      (t/is (m/delta-eq 1.6558087948 Z))
      (t/is (m/delta-eq 0.0977605462 p-value))))
  (t/testing "kurtosis-test vs R moments::anscombe.test"
    (let [{:keys [Z p-value]} (sut/kurtosis-test norm-xs)]
      (t/is (m/delta-eq 0.1687582795 Z))
      (t/is (m/delta-eq 0.8659867757 p-value))))
  (t/testing "normality-test (K^2 omnibus): composed from the verified skew-Z/kurt-Z above,
cross-checked against an independent R pchisq(Z1^2+Z2^2, df=2) call"
    (let [{:keys [Z p-value]} (sut/normality-test norm-xs)]
      (t/is (m/delta-eq 2.7701821217 Z))
      (t/is (m/delta-eq 0.2503010061 p-value))))
  (t/testing "jarque-bera-test vs R moments::jarque.test / tseries::jarque.bera.test"
    (let [{:keys [Z p-value]} (sut/jarque-bera-test norm-xs)]
      (t/is (m/delta-eq 2.4522418247 Z))
      (t/is (m/delta-eq 0.2934286082 p-value))))
  (t/testing "bonett-seier-test: formula verified via independent R reimplementation, given
fastmath's own :geary kurtosis estimator as input"
    (t/are [res sides] (m/delta-eq res (:p-value (sut/bonett-seier-test norm-xs {:sides sides})))
      0.3594890055 :two-sided
      0.8202554972 :one-sided-greater
      0.1797445028 :one-sided-less)
    (t/is (m/delta-eq -0.9163392158 (:Z (sut/bonett-seier-test norm-xs))))))

;; binomial-test: reference R base binom.test(37, 50, p=0.5). p-value is independent of
;; :ci-method; :confidence-interval is checked with :ci-method :clopper-pearson (R's own
;; default) and with fastmath's default :asymptotic (verified against a manual Wald-interval
;; R computation).
(t/deftest binomial-test-test
  (t/testing "p-value (independent of :ci-method) vs R binom.test"
    (t/are [res sides] (m/delta-eq res (:p-value (sut/binomial-test 37 50 {:sides sides})))
      0.0009362229 :two-sided
      0.0004681115 :one-sided-greater
      0.9998470680 :one-sided-less))
  (t/testing "confidence-interval, :clopper-pearson (matches R binom.test's default exact CI)"
    ;; bug fixed this session (blame: implementation bug): the confidence-interval computation
    ;; passed (1.0 - alpha) -- a confidence *level* -- into binomial-ci's `alpha` parameter,
    ;; which expects a significance level directly (binomial-ci's own docstring: "confidence
    ;; level is 1 - alpha"). This inverted alpha=0.05 into effectively alpha=0.95, producing a
    ;; drastically too-narrow interval (e.g. two-sided [0.7245,0.7521] instead of the correct
    ;; [0.5966,0.8537]). Fixed by passing alpha (and 2*alpha for the one-sided cases) directly.
    (t/is (seq-delta-eq (:confidence-interval (sut/binomial-test 37 50 {:ci-method :clopper-pearson}))
                         [0.5965523213 0.8536994156] 1.0e-6)
          "previously [0.7245230438471877 0.7521353831707929]")
    (t/is (seq-delta-eq (:confidence-interval (sut/binomial-test 37 50 {:ci-method :clopper-pearson
                                                                         :sides :one-sided-greater}))
                         [0.6187364440 1.0] 1.0e-6))
    (t/is (seq-delta-eq (:confidence-interval (sut/binomial-test 37 50 {:ci-method :clopper-pearson
                                                                         :sides :one-sided-less}))
                         [0.0 0.8388254035] 1.0e-6)))
  (t/testing "confidence-interval, default :asymptotic (Wald interval, manual R computation)"
    (t/is (seq-delta-eq (:confidence-interval (sut/binomial-test 37 50))
                         [0.6184190248 0.8615809752] 1.0e-6))))

;; t-test-one-sample / z-test-one-sample: reference R t.test(x, mu=5) (one-sample) for the
;; t-test, and a manual R z-test reimplementation (same stderr formula, qnorm/pnorm instead of
;; qt/pt) for the z-test, since R has no built-in one/two-sample z-test.
(t/deftest one-sample-t-z-test-test
  (t/testing "t-test-one-sample vs R t.test(x, mu=5)"
    (t/are [res k sides] (m/delta-eq res (get (sut/t-test-one-sample norm-xs {:mu 5.0 :sides sides}) k))
      1.1551049817 :t :two-sided
      0.2536493753 :p-value :two-sided
      0.1268246876 :p-value :one-sided-greater
      0.8731753124 :p-value :one-sided-less)
    (t/is (= 49 (:df (sut/t-test-one-sample norm-xs {:mu 5.0}))))
    (t/is (seq-delta-eq (:confidence-interval (sut/t-test-one-sample norm-xs {:mu 5.0}))
                         [4.6108999948 6.4411000052] 1.0e-6)))
  (t/testing "z-test-one-sample vs manual R (same stderr, qnorm/pnorm)"
    (let [{:keys [z p-value confidence-interval]} (sut/z-test-one-sample norm-xs {:mu 5.0})]
      (t/is (m/delta-eq 1.1551049817 z))
      (t/is (m/delta-eq 0.2480474384 p-value))
      (t/is (seq-delta-eq confidence-interval [4.6334914642 6.4185085358] 1.0e-6)))))

;; t-test-two-samples / z-test-two-samples / f-test: unequal-n groups (nx=40, ny=25), to
;; distinguish the pooled (:equal-variances? true) vs Welch/Satterthwaite (default) standard
;; error formulas (which happen to coincide algebraically when nx=ny, so an nx=ny fixture
;; wouldn't actually exercise the difference). Reference: R t.test(x,y)/t.test(x,y,var.equal=T)/
;; t.test(x,y,paired=T)/var.test(x,y) for the t/F tests; manual R z-test reimplementation
;; (same stderr formulas as the t-tests, qnorm/pnorm) for the z-tests.
(def two-sample-x [3.9 1.7 2.6 2.7 1.7 2.5 3.9 0.0 4.9 2.5 11.0 7.3 1.6 9.4 5.9 10.1 2.9 5.0 5.1 14.3
                   6.1 6.7 7.4 12.6 3.6 3.7 6.4 1.7 9.7 9.3 5.8 7.9 7.4 5.8 2.8 10.5 3.6 2.2 2.7 0.0])
(def two-sample-y [5.1 3.2 4.0 6.1 2.9 3.5 7.2 1.1 6.8 4.4 10.3 8.1 2.0 8.9 6.4 9.2 3.5 5.8 6.0 13.1
                   6.9 7.5 8.0 11.8 4.1])

;; n=50, paired 1:1 with norm-xs, for the paired-test / f-test cases below (which need equal-n
;; samples).
(def ph-y [5.1 3.2 4.0 6.1 2.9 3.5 7.2 1.1 6.8 4.4 10.3 8.1 2.0 8.9 6.4 9.2 3.5 5.8 6.0 13.1
           6.9 7.5 8.0 11.8 4.1 4.5 6.9 2.4 8.9 9.9 6.2 8.5 7.9 6.3 3.4 9.8 4.2 2.9 3.3 1.4
           7.8 3.9 8.3 4.2 7.1 5.9 7.7 8.6 5.0 8.7])

(t/deftest two-sample-t-z-f-test-test
  (t/testing "t-test-two-samples: Welch (default) vs R t.test(x,y)"
    (let [{:keys [t df p-value confidence-interval]} (sut/t-test-two-samples two-sample-x two-sample-y)]
      (t/is (m/delta-eq -1.0625991202 t))
      (t/is (m/delta-eq 56.8024865621 df))
      (t/is (m/delta-eq 0.2924606052 p-value))
      (t/is (seq-delta-eq confidence-interval [-2.4908862895 0.7638862895] 1.0e-6))))
  (t/testing "t-test-two-samples: Student (:equal-variances? true) vs R t.test(x,y,var.equal=TRUE)"
    (let [{:keys [t df p-value confidence-interval]} (sut/t-test-two-samples two-sample-x two-sample-y
                                                                              {:equal-variances? true})]
      (t/is (m/delta-eq -1.0256630948 t))
      (t/is (m/delta-eq 63.0 df))
      (t/is (m/delta-eq 0.3089725590 p-value))
      (t/is (seq-delta-eq confidence-interval [-2.5458916812 0.8188916812] 1.0e-6))))
  (t/testing "t-test-two-samples: paired vs R t.test(x,y,paired=TRUE) (nx=ny=50 fixture)"
    (let [{:keys [t df p-value confidence-interval]} (sut/t-test-two-samples norm-xs ph-y {:paired? true})]
      (t/is (m/delta-eq -5.4948149937 t))
      (t/is (= 49 df))
      (t/is (m/delta-eq 1.3903E-6 p-value 1.0e-9))
      (t/is (seq-delta-eq confidence-interval [-0.9095709136 -0.4224290864] 1.0e-6))))
  (t/testing "z-test-two-samples: Welch-stderr vs manual R z-test"
    (let [{:keys [z p-value confidence-interval]} (sut/z-test-two-samples two-sample-x two-sample-y)]
      (t/is (m/delta-eq -1.0625991202 z))
      (t/is (m/delta-eq 0.2879637862 p-value))
      (t/is (seq-delta-eq confidence-interval [-2.4562256747 0.7292256747] 1.0e-6))))
  (t/testing "z-test-two-samples: pooled-stderr (:equal-variances? true) vs manual R z-test"
    (let [{:keys [z p-value confidence-interval]} (sut/z-test-two-samples two-sample-x two-sample-y
                                                                           {:equal-variances? true})]
      (t/is (m/delta-eq -1.0256630948 z))
      (t/is (m/delta-eq 0.3050504157 p-value))
      (t/is (seq-delta-eq confidence-interval [-2.5135826725 0.7865826725] 1.0e-6))))
  (t/testing "f-test vs R var.test(x,y) (nx=ny=50 fixture)"
    (t/are [res k sides] (m/delta-eq res (get (sut/f-test norm-xs ph-y {:sides sides}) k))
      1.4052907085 :F :two-sided
      0.2372174573 :p-value :two-sided
      0.1186087286 :p-value :one-sided-greater)
    (t/is (seq-delta-eq (:confidence-interval (sut/f-test norm-xs ph-y))
                         [0.7974689881 2.4763871751] 1.0e-6))
    (t/is (seq-delta-eq (:confidence-interval (sut/f-test norm-xs ph-y {:sides :one-sided-greater}))
                         [0.8743233506 ##Inf] 1.0e-6))))

;; one-way-anova-test / levene-test / brown-forsythe-test / fligner-killeen-test: reference R
;; oneway.test(y~g, var.equal=TRUE), car::leveneTest(y, g, center=mean/median), and base R
;; fligner.test(y, g). 3 unequal-size groups.
(def anova-g1 [3.9 1.7 2.6 2.7 1.7 2.5 3.9 0.0 4.9 2.5 11.0 7.3 1.6 9.4 5.9])
(def anova-g2 [5.1 3.2 4.0 6.1 2.9 3.5 7.2 1.1 6.8 4.4 10.3 8.1 2.0 8.9 6.4 9.2 3.5 5.8])
(def anova-g3 [6.1 6.7 7.4 12.6 3.6 3.7 6.4 1.7 9.7 9.3 5.8 7.9])
(def anova-groups [anova-g1 anova-g2 anova-g3])

(t/deftest anova-variance-homogeneity-test
  (t/testing "one-way-anova-test vs R oneway.test(y~g, var.equal=TRUE)"
    (let [{:keys [F df p-value]} (sut/one-way-anova-test anova-groups)]
      (t/is (m/delta-eq 2.8234530556 F))
      (t/is (= [2 42] df))
      (t/is (m/delta-eq 0.0707141429 p-value)))
    (t/is (m/delta-eq 0.1414282858 (:p-value (sut/one-way-anova-test anova-groups {:sides :two-sided}))))
    (t/is (m/delta-eq 0.9292858571 (:p-value (sut/one-way-anova-test anova-groups {:sides :one-sided-less})))))
  (t/testing "levene-test vs R car::leveneTest(y, g, center=mean)"
    (let [{:keys [W p-value]} (sut/levene-test anova-groups)]
      (t/is (m/delta-eq 0.0812630690 W))
      (t/is (m/delta-eq 0.9220957202 p-value))))
  (t/testing "brown-forsythe-test vs R car::leveneTest(y, g) (default center=median)"
    (let [{:keys [W p-value]} (sut/brown-forsythe-test anova-groups)]
      (t/is (m/delta-eq 0.0062233907 W))
      (t/is (m/delta-eq 0.9937968507 p-value))))
  (t/testing "fligner-killeen-test vs R base fligner.test(y, g)"
    (let [{:keys [chi2 df p-value]} (sut/fligner-killeen-test anova-groups)]
      (t/is (m/delta-eq 0.0454724155 chi2))
      (t/is (= 2 df))
      (t/is (m/delta-eq 0.9775203121 p-value)))))

;; ad-test-one-sample / ks-test-one-sample: reference R goftest::ad.test(x, "pnorm", 0, 1)
;; (fixed, not estimated, parameters -- matches fastmath's convention of testing against the
;; distribution object directly) and base R ks.test(x, "pnorm", 0, 1). Tie-free n=30 sample
;; (KS test warns/is inaccurate with ties, and fastmath's default :distinct? true silently
;; dedupes, so a tie-free fixture avoids an apples-to-oranges comparison).
(def rank-x [1.9452 -0.3776 0.7358 1.0594 0.7851 0.1727 2.1138 0.1864 2.7221 0.2247 1.8658 3.044
             -1.3666 -0.0345 0.14 1.0631 -0.0411 -2.8877 -2.6286 1.8841 -0.068 -1.8376 0.0937
             1.7576 2.5742 -0.2166 -0.0087 -1.8158 0.8521 -0.468])

(t/deftest ad-ks-one-sample-test
  (t/testing "ad-test-one-sample vs R goftest::ad.test (fixed N(0,1) parameters)"
    (let [{:keys [A2 p-value]} (sut/ad-test-one-sample rank-x)]
      (t/is (m/delta-eq 4.6172523880 A2))
      (t/is (m/delta-eq 0.0044597241 p-value))))
  (t/testing "ks-test-one-sample vs R ks.test(x, \"pnorm\", 0, 1)"
    (t/is (m/delta-eq 0.2272588700 (:d (sut/ks-test-one-sample rank-x))))
    (t/is (m/delta-eq 0.0763805103 (:p-value (sut/ks-test-one-sample rank-x))))
    (let [greater (sut/ks-test-one-sample rank-x r/default-normal {:sides :right})
          less (sut/ks-test-one-sample rank-x r/default-normal {:sides :left})]
      (t/is (m/delta-eq 0.0986328047 (:stat greater)))
      (t/is (m/delta-eq 0.5239632376 (:p-value greater)))
      (t/is (m/delta-eq 0.2272588700 (:stat less)))
      (t/is (m/delta-eq 0.0381915796 (:p-value less))))))

;; ks-test-two-samples / kruskal-test: reference base R ks.test(x,y) (:exact, nx*ny=600<10000
;; matches fastmath's default) and ks.test(x,y,exact=FALSE) (:approximate), and
;; kruskal.test(list(g1,g2,g3)) (reuses group 33's anova-groups fixture). kruskal-test's
;; pre-existing reference values had an LLM-sourced provenance flag (elevated re-verification
;; priority); both the H statistic and R's default one-tailed p-value convention (fastmath's
;; :sides :right default, "as in R" per the source comment) are confirmed exact here.
(def rank-y [0.4 1.1 -0.9 0.3 2.2 -1.5 0.8 1.7 -0.2 0.05 1.3 -0.6 0.95 2.6 -1.1 0.15 1.9 -2.0 0.65 3.1])

(t/deftest ks-two-sample-kruskal-test-test
  (t/testing "ks-test-two-samples, :exact method (default, nx*ny=600<10000) vs R ks.test(x,y)"
    (let [two-sided (sut/ks-test-two-samples rank-x rank-y)
          greater (sut/ks-test-two-samples rank-x rank-y {:sides :right})
          less (sut/ks-test-two-samples rank-x rank-y {:sides :left})]
      (t/is (= :exact (:method two-sided)))
      (t/is (m/delta-eq 0.1666666667 (:d two-sided)))
      (t/is (m/delta-eq 0.8602250677 (:p-value two-sided)))
      (t/is (m/delta-eq 0.4801172711 (:p-value greater)))
      (t/is (m/delta-eq 0.0833333333 (:stat less)))
      (t/is (m/delta-eq 0.8097102703 (:p-value less)))))
  (t/testing "ks-test-two-samples, :approximate method vs R ks.test(x,y,exact=FALSE)"
    (t/is (m/delta-eq 0.8927783373 (:p-value (sut/ks-test-two-samples rank-x rank-y {:method :approximate})))))
  (t/testing "kruskal-test vs R kruskal.test(list(g1,g2,g3)) -- previously LLM-sourced reference
values, now confirmed exact via genuine R computation"
    (let [{:keys [stat df p-value]} (sut/kruskal-test anova-groups)]
      (t/is (m/delta-eq 5.6919380356 stat))
      (t/is (= 2 df))
      (t/is (m/delta-eq 0.0580779609 p-value)))))

;; box-cox-transformation / yeo-johnson-transformation / power-transformation /
;; modified-power-transformation. Reference: R car::bcPower/yjPower for the not-scaled forms.
;; :scaled? verified via a manual R reproduction of the documented Jacobian-rescaling formula
;; (Ecfun::BoxCox's own reference source). :negative? verified against the *original*
;; Bickel & Doksum (1981) formula sign(x)*|x|^lambda-1)/lambda (Pengfei Li's "Box-Cox
;; Transformation: An Overview", which explicitly attributes this exact formula to Bickel &
;; Doksum 1981) -- an initial manual R attempt used a different, superficially similar formula
;; (sign(x)*(|x|^lambda-1)/lambda, the convention some forecasting packages also label
;; "Bickel-Doksum"), which is a *different* named transformation, not what fastmath implements
;; or documents; confirmed fastmath matches the original 1981 paper's own formula, not a bug.
(def pt-x [0.5 1.2 2.3 3.1 4.8 5.5 6.2 7.9 8.1 9.3])
(def pt-xn [-3.2 -1.5 0.0 0.5 1.2 2.3 3.1 4.8 5.5 6.2])
(def pt-xn2 [-3.2 -1.5 0.5 1.2 2.3 3.1])

(t/deftest power-transformation-test
  (t/testing "box-cox-transformation, not-scaled (default) vs R car::bcPower"
    (t/is (seq-delta-eq (sut/box-cox-transformation pt-x 0.5)
                         [-0.5857864376 0.1908902300 1.0331501776 1.5213633723 2.3817804600
                          2.6904157598 2.9799598392 3.6213877290 3.6920997883 4.0991802728] 1.0e-6))
    (t/is (seq-delta-eq (sut/box-cox-transformation pt-x 0.0)
                         [-0.6931471806 0.1823215568 0.8329091229 1.1314021115 1.5686159179
                          1.7047480922 1.8245492921 2.0668627595 2.0918640617 2.2300144002] 1.0e-6)))
  (t/testing "box-cox-transformation, :scaled? true, vs manual R (Jacobian-rescaled) formula"
    (t/is (seq-delta-eq (sut/box-cox-transformation pt-x 0.5 {:scaled? true})
                         [-1.1187444279 0.3645652536 1.9731269456 2.9055244135 4.5487629058
                          5.1381996010 5.6911755742 6.9161849489 7.0512319852 7.8286808889] 1.0e-6))
    (t/is (seq-delta-eq (sut/box-cox-transformation pt-x 0.0 {:scaled? true})
                         [-2.5281835813 0.6649992661 3.0379509986 4.1266737027 5.7213664288
                          6.2178946376 6.6548572840 7.5386709199 7.6298606175 8.1337498740] 1.0e-6)))
  (t/testing "box-cox-transformation, :negative? true, vs the original Bickel & Doksum (1981) formula"
    (t/is (seq-delta-eq (sut/box-cox-transformation pt-xn2 0.5 {:negative? true})
                         [-5.5777087640 -4.4494897428 -0.5857864376 0.1908902300 1.0331501776 1.5213633723] 1.0e-6))
    (t/is (seq-delta-eq (sut/box-cox-transformation pt-xn2 0.0 {:negative? true})
                         [-1.4350845253 -0.9162907319 0.4054651081 0.7884573604 1.1939224685 1.4109869737] 1.0e-6)))
  (t/testing "yeo-johnson-transformation vs R car::yjPower"
    (t/is (seq-delta-eq (sut/yeo-johnson-transformation pt-xn 0.5)
                         [-5.0716257623 -1.9685647168 0.0 0.4494897428 0.9664793948 1.6331804249
                          2.0496913463 2.8166378315 3.0990195136 3.3665631460] 1.0e-6))
    (t/is (seq-delta-eq (sut/yeo-johnson-transformation pt-xn 2.0)
                         [-1.4350845253 -0.9162907319 0.0 0.6250000000 1.9200000000 4.9450000000
                          7.9050000000 16.3200000000 20.6250000000 25.4200000000] 1.0e-6))
    (t/is (seq-delta-eq (sut/yeo-johnson-transformation pt-xn 0.0)
                         [-8.3200000000 -2.6250000000 0.0 0.4054651081 0.7884573604 1.1939224685
                          1.4109869737 1.7578579176 1.8718021769 1.9740810260] 1.0e-6)))
  (t/testing "inverse round-trips (forward then :inverse? true recovers the original data)"
    (t/is (v/delta-eq pt-x (sut/box-cox-transformation (sut/box-cox-transformation pt-x 0.5) 0.5 {:inverse? true})))
    (let [gm (sut/geomean pt-x)]
      (t/is (v/delta-eq pt-x (sut/box-cox-transformation
                              (sut/box-cox-transformation pt-x 0.5 {:scaled? true})
                              0.5 {:scaled? gm :inverse? true}))))
    (t/is (v/delta-eq pt-xn2 (sut/box-cox-transformation
                              (sut/box-cox-transformation pt-xn2 0.5 {:negative? true})
                              0.5 {:negative? true :inverse? true})))
    (t/is (v/delta-eq pt-xn (sut/yeo-johnson-transformation
                             (sut/yeo-johnson-transformation pt-xn 0.5) 0.5 {:inverse? true}))))
  (t/testing "scaled inverse validation"
    ;; bug fixed this session (blame: implementation bug): inverting a :scaled? true (the
    ;; auto-compute-geometric-mean sentinel used by the forward transform) transformation threw
    ;; a raw ClassCastException ("Boolean cannot be cast to Number") instead of a clear error --
    ;; the geometric mean genuinely can't be recovered from already-transformed data, so the
    ;; caller must supply the actual numeric gm used forward; now throws an informative ex-info.
    (t/is (thrown-with-msg? clojure.lang.ExceptionInfo #"requires the actual numeric geometric mean"
                             (sut/box-cox-transformation [1.0] 0.5 {:scaled? true :inverse? true}))))
  (t/testing "power-transformation / modified-power-transformation (deprecated) delegate correctly"
    (t/is (= (sut/box-cox-transformation pt-x 0.5 {:scaled? true}) (sut/power-transformation pt-x 0.5)))
    (t/is (= (sut/box-cox-transformation pt-xn2 0.5 {:negative? true}) (sut/modified-power-transformation pt-xn2 0.5)))))

;; cohens-kappa / weighted-kappa: general (3x3, ordinal) contingency table, not just the 2x2
;; special case exercised above. Reference: R irr::kappa2(df, weight=...) -- "unweighted",
;; "equal" (fastmath's :equal-spacing), "squared" (fastmath's :fleiss-cohen).
(def agreement-r1 [1 1 1 1 2 2 2 2 2 3 3 3 1 2 3 1 2 3 2 1])
(def agreement-r2 [1 1 2 1 2 2 1 3 2 3 3 2 1 2 3 2 2 3 3 1])

(t/deftest agreement-test
  (let [ct (sut/contingency-table agreement-r1 agreement-r2)]
    (t/testing "cohens-kappa"
      (t/is (m/delta-eq 0.5454545454545454 (sut/cohens-kappa agreement-r1 agreement-r2)))
      (t/is (m/delta-eq 0.5454545454545454 (sut/cohens-kappa ct))))
    (t/testing "weighted-kappa"
      ;; bug fixed this session (blame: implementation bug): weighted-kappa-equal-spacing's
      ;; `1 - |id1-id2|/R` used integer division (all its arguments are longs), truncating
      ;; every fractional weight to 0 and leaving every off-diagonal weight at 1.0 -- silently
      ;; degenerating :equal-spacing (the *default* weighting scheme) toward unweighted
      ;; behavior, and producing ##NaN outright whenever that made pe reach exactly 1.0
      ;; (as happened on this fixture, category indices 1..3 with R=3: |id1-id2| never
      ;; reaches R, so integer division truncated *every* off-diagonal weight to 0).
      (t/is (m/delta-eq 0.6428571428571433 (sut/weighted-kappa ct)) "default is :equal-spacing")
      (t/is (m/delta-eq 0.6428571428571433 (sut/weighted-kappa ct :equal-spacing)))
      (t/is (m/delta-eq 0.7500000000000009 (sut/weighted-kappa ct :fleiss-cohen)))
      (t/testing "custom weight sources agree with the equivalent built-in keyword"
        (t/is (m/delta-eq 0.7500000000000009 (sut/weighted-kappa ct (var-get #'sut/weighted-kappa-fleiss-cohen))))
        (let [equal-spacing-map (into {} (for [i [1 2 3] j [1 2 3]]
                                            [[i j] (- 1.0 (/ (double (Math/abs (- i j))) 3.0))]))]
          (t/is (m/delta-eq 0.6428571428571433 (sut/weighted-kappa ct equal-spacing-map))))))))

;; durbin-watson: reference computed two independent ways in R -- the manual formula
;; sum(diff(rs)^2)/sum(rs^2), and the car::durbinWatsonTest(rs) package function.
(t/deftest durbin-watson-test
  (t/is (m/delta-eq 2.666486486486486
                     (sut/durbin-watson [2.1 -0.5 1.3 -1.8 0.2 0.9 -1.1 2.4 -0.3 1.0]))))

;; classification: confusion-matrix / binary-measures(-all) / multiclass-measure / auc /
;; auc-roc / binary-measures-thr / multiclass-auc. Reference: Python sklearn.metrics
;; (confusion_matrix, accuracy_score/precision_score/recall_score/f1_score,
;; matthews_corrcoef, balanced_accuracy_score, roc_auc_score, roc_curve, precision_recall_curve).
;; The individual metric formulas inside binary-measures-all live in fastmath.stats.binary
;; (out of this audit's scope, has its own test suite); what's verified here is that these
;; fastmath.stats wrapper functions dispatch correctly across their input-format arities and
;; produce values matching an independent reference implementation.
(def cls-actual [1 1 1 0 0 1 0 0 1 0 1 0 0 1 0])
(def cls-pred   [1 0 1 0 1 1 0 0 1 0 0 0 1 1 0])

(t/deftest confusion-matrix-and-binary-measures-test
  (t/testing "confusion-matrix: all input formats agree, and match sklearn.metrics.confusion_matrix"
    (let [cm (sut/confusion-matrix cls-actual cls-pred)]
      (t/is (= {:tp 5 :fn 2 :fp 2 :tn 6} cm))
      (t/is (= cm (sut/confusion-matrix (:tp cm) (:fn cm) (:fp cm) (:tn cm))))
      (t/is (= cm (sut/confusion-matrix [(:tp cm) (:fn cm) (:fp cm) (:tn cm)])))
      (t/is (= cm (sut/confusion-matrix [[(:tp cm) (:fp cm)] [(:fn cm) (:tn cm)]])))))
  (t/testing "binary-measures / binary-measures-all: headline metrics vs sklearn"
    (let [bm (sut/binary-measures cls-actual cls-pred)
          bma (sut/binary-measures-all cls-actual cls-pred)]
      (t/is (m/delta-eq 0.7333333333333333 (:accuracy bm)))
      (t/is (m/delta-eq 0.7142857142857143 (:precision bm)))
      (t/is (m/delta-eq 0.7142857142857143 (:recall bm)))
      (t/is (m/delta-eq 0.7142857142857143 (:f-measure bm)))
      (t/is (m/delta-eq 0.4642857142857143 (:mcc bma)) "matthews_corrcoef")
      (t/is (m/delta-eq 0.7321428571428572 (:ba bma)) "balanced_accuracy_score"))))

;; multiclass-measure: 3-class macro/micro/weighted averages of f1/precision/recall vs
;; sklearn.metrics.{f1,precision,recall}_score(average=...).
(def mc-actual ["a" "a" "a" "b" "b" "b" "c" "c" "c" "c" "a" "b" "c" "a" "b"])
(def mc-pred   ["a" "b" "a" "b" "b" "c" "c" "a" "c" "c" "a" "b" "b" "a" "b"])

(t/deftest multiclass-measure-test
  (t/are [res metric average weighted?] (m/delta-eq res (sut/multiclass-measure
                                                           mc-actual mc-pred
                                                           {:metric metric :average average :weighted? weighted?}))
    0.7313131313131312 :f1-score sut/mean false
    0.7333333333333333 :f1-score :micro false
    0.7313131313131312 :f1-score sut/mean true
    0.7388888888888889 :precision sut/mean false
    0.7333333333333333 :precision :micro false
    0.7388888888888888 :precision sut/mean true
    0.7333333333333334 :recall sut/mean false
    0.7333333333333333 :recall :micro false
    0.7333333333333333 :recall sut/mean true))

;; binary-measures-thr / auc / auc-roc: reference via sklearn.metrics.roc_curve/
;; precision_recall_curve/roc_auc_score. Uses a tie-free score sequence: fastmath's
;; documented linear tie-interpolation in binary-measures-thr produces a (correctly)
;; different, finer polyline than sklearn's one-point-per-tied-group convention when scores
;; are tied, verified separately below by comparing AUC convergence as ties are removed.
(def roc-labels [1 0 1 1 0 0 1 0 1 0 0 1 0 1 0 1 0 0 1 0])
(def roc-scores [0.91 0.41 0.31 0.65 0.81 0.51 0.71 0.61 0.21 0.45 0.55 0.35 0.15 0.75 0.25 0.85 0.11 0.95 0.605 0.305])

(t/deftest binary-measures-thr-and-auc-test
  (t/testing "binary-measures-thr: fpr/tpr at every threshold match sklearn.metrics.roc_curve
(sklearn additionally drops non-corner intermediate points; every threshold sklearn does keep
matches exactly)"
    (let [perf (sut/binary-measures-thr roc-labels roc-scores)
          by-thr (zipmap (:thr perf) (map vector (:fpr perf) (:tpr perf)))]
      (t/are [thr fpr tpr] (let [[afpr atpr] (by-thr thr)]
                             (and (m/delta-eq fpr afpr) (m/delta-eq tpr atpr)))
        0.95  0.09090909090909091 0.0
        0.85  0.09090909090909091 0.2222222222222222
        0.81  0.18181818181818182 0.2222222222222222
        0.65  0.18181818181818182 0.5555555555555556
        0.61  0.2727272727272727  0.5555555555555556
        0.605 0.2727272727272727  0.6666666666666666
        0.41  0.6363636363636364  0.6666666666666666
        0.31  0.6363636363636364  0.8888888888888888
        0.25  0.8181818181818182  0.8888888888888888
        0.21  0.8181818181818182  1.0
        0.11  1.0                 1.0)))
  (t/testing "auc-roc (direct U-statistic) matches sklearn.metrics.roc_auc_score"
    (t/is (m/delta-eq 0.6565656565656566 (sut/auc-roc roc-labels roc-scores))))
  (t/testing "auc (trapezoidal, via binary-measures-thr curve) matches roc_auc_score and
matches auc-roc"
    (let [perf (sut/binary-measures-thr roc-labels roc-scores)]
      (t/is (m/delta-eq 0.6565656565656566 (sut/auc perf)) "default is ROC (:fpr :tpr)")
      (t/is (m/delta-eq 0.5373006895065718 (sut/auc perf :recall :precision))
            "matches sklearn.metrics.auc(*precision_recall_curve(...))"))))

;; multiclass-auc: 3-class one-vs-rest AUC, macro/weighted/micro, vs
;; sklearn.metrics.roc_auc_score(y_true_bin, scores, multi_class="ovr", average=...).
(def auc-classes ["a" "a" "a" "b" "b" "b" "c" "c" "c" "c" "a" "b" "c" "a" "b"])
(def auc-scores {"a" [0.7 0.3 0.6 0.2 0.3 0.1 0.2 0.5 0.1 0.2 0.6 0.2 0.3 0.5 0.2]
                  "b" [0.2 0.5 0.3 0.6 0.4 0.2 0.2 0.3 0.1 0.1 0.3 0.7 0.3 0.4 0.6]
                  "c" [0.1 0.2 0.1 0.2 0.3 0.7 0.6 0.2 0.8 0.7 0.1 0.1 0.4 0.1 0.2]})

(t/deftest multiclass-auc-test
  (t/is (m/delta-eq 0.8866666666666667 (sut/multiclass-auc auc-classes auc-scores {:average sut/mean})))
  (t/is (m/delta-eq 0.8866666666666667 (sut/multiclass-auc auc-classes auc-scores {:average sut/mean :weighted? true})))
  (t/is (m/delta-eq 0.8777777777777779 (sut/multiclass-auc auc-classes auc-scores {:average :micro}))))

;; histogram

(def one {:size 1, :step 0.0, :samples 1, :min 1.0, :max 1.0, :bins '([1.0 1])
        :bins-maps '({:min 1.0, :max 1.0, :mid 1.0, :count 1, :step 0.0, :avg 1.0 :probability 1.0})
        :frequencies {1.0 1} :intervals '(1.0 1.0)})
(def one2 {:size 1, :step 1.0, :samples 2, :min 1.0, :max 2.0, :bins '([1.0 2])
         :bins-maps '({:min 1.0, :mid 1.5, :max 2.0, :count 2, :step 1.0, :avg 1.5 :probability 1.0})
         :frequencies {1.5 2} :intervals '(1.0 2.0)})
(def two {:size 2, :step 0.5, :samples 2, :min 1.0, :max 2.0, :bins '([1.0 1] [1.5 1])
        :bins-maps '({:min 1.0, :mid 1.25, :max 1.5, :count 1, :step 0.5, :avg 1.0 :probability 0.5}
                     {:min 1.5, :mid 1.75, :max 2.0, :count 1, :step 0.5, :avg 2.0 :probability 0.5})
        :frequencies {1.0 1, 2.0 1} :intervals '(1.0 1.5 2.0)})

(t/deftest histogram-tests
  (t/are [in method res] (= res (sut/histogram in method))
    [1 1 1 1] :sqrt  {:size 1, :step 0.0, :samples 4, :min 1.0, :max 1.0, :bins '([1.0 4])
                      :bins-maps '({:min 1.0, :max 1.0, :mid 1.0, :count 4, :step 0.0, :avg 1.0
                                    :probability 1.0}),
                      :frequencies {1.0 4} :intervals '(1.0 1.0)}
    [1 1] :sturges {:size 1, :step 0.0, :samples 2, :min 1.0, :max 1.0, :bins '([1.0 2])
                    :bins-maps '({:min 1.0, :max 1.0, :mid 1.0, :count 2, :step 0.0, :avg 1.0 :probability 1.0}),
                    :frequencies {1.0 2} :intervals '(1.0 1.0)}
    [1] :rice one
    [1] :doane one
    [1 2] :sqrt one2
    [1 2] :sturges two
    [1 2] :rice two
    [1 2] :doane one2
    [1 2] :scott one2
    [1 2] :freedman-diaconis one2))

;; estimate-bins (all 6 methods) and a full histogram on real data (mtcars$mpg)
;; reference values from R: nclass.Sturges/nclass.scott (sqrt/sturges/scott match R's exact
;; convention); rice/freedman-diaconis instead verified against fastmath's own documented
;; formula (ceil(2*cbrt(n)) resp. ceil(range/(2*IQR(type=6)/cbrt(n)))) computed in R -- these
;; deliberately differ from R's own nclass.Sturges-family rounding/IQR-type conventions (the
;; same ceiling-rounding and :legacy/type-6-percentile patterns already established elsewhere
;; in this file, e.g. Percentiles & Quantiles, Mode groups), not bugs; doane verified against
;; R's e1071::skewness(type=1)-based formula: raw value ceil(7.390837)=8 (Fastmath Stats Bins
;; Audit, 2026-09-23 -- corrected from a previously-recorded 7, which reflected a real
;; fastmath.stats.bins/doane bug (truncating instead of ceiling the raw formula value,
;; unlike its sibling bin-count estimators), independently confirmed via numpy's
;; histogram_bin_edges(bins='doane'); now fixed in fastmath.stats.bins/doane, see CHANGELOG)

(t/deftest estimate-bins-and-full-histogram-test
  (let [mpg (mtcars :mpg)]
    (t/are [method n] (= n (sut/estimate-bins mpg method))
      :sqrt 5
      :sturges 6
      :rice 7
      :scott 4
      :freedman-diaconis 5
      :doane 8)
    (t/is (= 5 (sut/estimate-bins mpg 5)) "explicit long bypasses estimation")
    (t/is (= (sut/estimate-bins mpg) (sut/estimate-bins mpg :freedman-diaconis)) "default method")
    (t/testing "full histogram (:sturges, 6 bins) against manual R binning"
      (let [h (sut/histogram mpg :sturges)]
        (t/is (= 6 (:size h)))
        (t/is (m/delta-eq 3.9166666667 (:step h)))
        (t/is (= 32 (:samples h)))
        (t/is (= [4 10 9 4 1 4] (mapv second (:bins h))))
        (t/is (v/delta-eq [12.1 16.1 20.344444444 24.0 27.3 31.775] (mapv :avg (:bins-maps h))))
        (t/is (= 32 (reduce + (map second (:bins h)))) "counts partition all samples")))))

;; dissimilarity (40 methods) / similarity (12 methods)
;; reference values: independent from-scratch Python re-implementation of every formula listed
;; in dissimilarity's/similarity's own docstrings (Cha's survey); P=[0.1 0.2 0.3 0.4],
;; Q=[0.4 0.3 0.2 0.1], :probabilities? false (already sum to 1). Edge cases (epsilon
;; substitution, :remove-zeros?, :bins raw-data-to-histogram path, distribution-object Q)
;; verified separately below by hand/R computation.

(def diss-P [0.1 0.2 0.3 0.4])
(def diss-Q [0.4 0.3 0.2 0.1])

(t/deftest dissimilarity-similarity-test
  (t/testing "dissimilarity, all 40 methods"
    (t/are [method res] (m/delta-eq res (sut/dissimilarity method diss-P diss-Q {:probabilities? false}))
      :euclidean 0.447213595500
      :city-block 0.800000000000
      :manhattan 0.800000000000
      :chebyshev 0.300000000000
      :minkowski 0.447213595500
      :sorensen 0.400000000000
      :gower 0.200000000000
      :soergel 0.571428571429
      :kulczynski 1.333333333333
      :canberra 1.600000000000
      :lorentzian 0.715348888544
      :non-intersection 0.400000000000
      :wave-hedges 2.166666666667
      :czekanowski 0.400000000000
      :motyka 0.700000000000
      :tanimoto 0.571428571429
      :jaccard 0.500000000000
      :dice 0.333333333333
      :bhattacharyya 0.116648487374
      :hellinger 0.663632583418
      :matusita 0.469259099951
      :squared-chord 0.220204102887
      :euclidean-sq 0.200000000000
      :squared-euclidean 0.200000000000
      :pearson-chisq 1.208333333333
      :chisq 1.208333333333
      :neyman-chisq 1.208333333333
      :squared-chisq 0.400000000000
      :symmetric-chisq 0.800000000000
      :divergence 1.600000000000
      :clark 0.894427191000
      :additive-symmetric-chisq 2.416666666667
      :kullback-leibler 0.456434819147
      :jeffreys 0.912869638294
      :k-divergence 0.106440135286
      :topsoe 0.212880270572
      :jensen-shannon 0.106440135286
      :jensen-difference 0.106440135286
      :taneja 0.121777274287
      :kumar-johnson 2.982603454360
      :avg 0.550000000000))
  (t/testing "similarity, all 12 methods"
    (t/are [method res] (m/delta-eq res (sut/similarity method diss-P diss-Q {:probabilities? false}))
      :intersection 0.600000000000
      :czekanowski 0.600000000000
      :motyka 0.300000000000
      :kulczynski 0.750000000000
      :ruzicka 0.428571428571
      :inner-product 0.200000000000
      :harmonic-mean 0.800000000000
      :cosine 0.666666666667
      :jaccard 0.500000000000
      :dice 0.666666666667
      :fidelity 0.889897948557
      :squared-chord 0.779795897113))
  (t/testing "edge cases"
    (let [Pz [0.2 0.3 0.5] Qz [0.2 0.0 0.8]]
      (t/is (m/delta= 90000.1125 (sut/dissimilarity :chisq Pz Qz {:probabilities? false}))
            "epsilon (1e-6) substituted for the zero denominator")
      (t/is (m/delta= 0.1125 (sut/dissimilarity :chisq Pz Qz {:probabilities? false :remove-zeros? true}))
            "the q=0 pair is dropped entirely instead"))
    (t/testing ":bins raw-data-to-histogram path"
      ;; hand-verified: both sequences binned into 4 shared-edge bins -> counts [2 2 3 1] and
      ;; [1 2 3 2], normalized to probabilities, then euclidean distance = sqrt(2*0.125^2)
      (t/is (m/delta-eq 0.1767766953
                        (sut/dissimilarity :euclidean [1.0 1 2 2 3 3 3 4] [1.0 2 2 3 3 3 4 4] {:bins 4}))))
    (t/testing "distribution-object Q-expected path"
      ;; reference: R pnorm() quantization of Normal(20,6) over the mtcars$mpg :sturges bin
      ;; edges (already R-verified in the Histogram group), euclidean vs. the mpg histogram
      (t/is (m/delta-eq 0.1757709379
                        (sut/dissimilarity :euclidean (mtcars :mpg)
                                           (r/distribution :normal {:mu 20 :sd 6}) {:bins 6}))))))

;; chatgpt :)
(t/deftest kruskal-wallis-test
  (t/testing "Basic Kruskal-Wallis test with distinct groups"
    (let [result (sut/kruskal-test [[1 2 3 4] [2 3 4 5] [6 7 8 9]])]
      (t/is (map? result))
      (t/is (contains? result :stat))
      (t/is (contains? result :p-value))
      (t/is (contains? result :df))
      (t/is (= (:df result) 2))
      (t/is (= (:k result) 3))
      (t/is (= (:sides result) :right))
      (t/is (< (:p-value result) 0.05) "Should detect significant difference")))

  (t/testing "Kruskal-Wallis test with identical groups (no difference expected)"
    (let [result (sut/kruskal-test [[5 5 5 5] [5 5 5 5] [5 5 5 1]])]
      (t/is (map? result))
      (t/is (contains? result :stat))
      (t/is (contains? result :p-value))
      (t/is (contains? result :df))
      (t/is (= (:df result) 2))
      (t/is (= (:k result) 3))
      (t/is (> (:p-value result) 0.05) "Should not detect significant difference"))))

(t/deftest transformations
  (let [ys [0.24 0.61 1 1.88 11.86 29.46 84.07 164.82 247.68]
        ys- [-1 0.24 0.61 1 1.88 11.86 29.46 84.07 164.82 247.68]]
    (t/is (v/delta-eq (sut/box-cox-transformation ys 0.65)
                      [-0.9300129 -0.4227523  0.0000000  0.7804770  6.1394250 12.3319046 25.8838035
                       40.9374869 53.8112999]))
    (t/is (v/delta-eq (sut/box-cox-transformation ys- 0.65 {:alpha 2.0})
                      [0 1.06018885 1.331677747 1.603605588 2.175413313 6.95789528 12.93691302
                       26.30609562 41.27180377 54.10140563]))
    (t/is (v/delta-eq (sut/box-cox-transformation ys 0.65 {:scaled? true})
                      [-1.9941886  -0.9064904   0.0000000   1.6735450  13.1645181  26.4427989
                       55.5015817  87.7805795 115.3853705]))
    (t/is (v/delta-eq (sut/box-cox-transformation ys- 0.65 {:scaled? true :alpha 2.0})
                      [0 2.578237269 3.238461899 3.899753984 5.290313778 16.9206693 31.46083955
                       63.97290082 100.3674984 131.5673715]))
    (t/is (v/delta-eq (sut/box-cox-transformation ys- 0.65 {:negative? true})
                      [-3.0769230769230766 -0.9300128623708235 -0.4227522822458623 0.0 0.7804770384569027
                       6.139425022491503 12.33190461964338 25.883803510261142 40.937486879512214
                       53.8112999422251]))
    (t/is (v/delta-eq (sut/yeo-johnson-transformation ys- 0.65)
                      [-1.147497226 0.2308761917 0.5581765097 0.8756433781 1.521329444 6.554235326
                       12.63614713 26.09538393 41.10482174 53.95645486]))
    ;; same as Julia BoxCoxTrans
    (t/is (m/delta-eq (sut/box-cox-infer-lambda ys) 0.01235497))
    ;; same as Python Scipy
    (t/is (m/delta-eq (sut/yeo-johnson-infer-lambda ys-) -0.004272723))))

;; reference values (all 46 assertions below): untrimmed L-moments/ratios/L-variation from R
;; lmom::samlmu(ys, nmom=8, ratios=FALSE/TRUE); trimmed (TL-moment) variants from R
;; lmomco::TLmoms(ys, nmom=8, trim=1 / leftrim=1,rightrim=0 / leftrim=0,rightrim=1)

(t/deftest l-moments
  (let [ys (sort [-1.7728323 , -1.34577139, -0.06550838, -3.06360177,  0.28895007,
                  1.36727242, -2.39761596, -0.63039948,  0.86878218,  0.68354586,
                  2.82116707,  1.04778304, -0.17578999, -2.32261317,  2.91884183,
                  0.59641717, -0.53668897, -1.09227543,  0.51320672, -4.122064  ,
                  0.20198867,  0.70337326, -0.09336329,  0.21434579,  0.06893888,
                  1.69565577,  2.27596985,  0.60637383,  0.27447569,  0.18262615,
                  -0.24239664, -0.47757852, -4.22660875,  0.92345263,  1.58138485,
                  -1.68156114,  0.0093263 , -0.5041279 ,  3.81543149, -0.82702232,
                  -0.43378152, -2.02045542, -1.07154893, -2.3172839 ,  0.96072838,
                  -0.18502813,  0.526105  ,  0.404828  , -0.18582972, -2.65091393,
                  -0.31282541, -0.19424187, -0.19954772, -0.16113933,  0.42566295,
                  1.44705338, -0.09691397,  0.28710832,  1.28102335,  0.09519895,
                  -1.61865923,  1.12958884,  0.06267352,  0.07283967, -0.66745461,
                  -0.48761821,  2.66665185, -0.27096777,  0.28232669,  0.72006836,
                  -1.39375443, -4.38246225,  7.93752809,  2.78756317, -2.2537608 ,
                  -4.89626616,  0.64561325, -2.4040938 ,  0.02983087, -1.28748567,
                  -0.74376858, -0.16892732, -0.03936628, -0.80297845,  0.69381382,
                  -0.22758131,  1.96957425,  0.29550671,  0.30425838,  3.40986819,
                  0.43560846,  1.83062494,  0.47820415, -2.35097606, -0.06590768,
                  -1.17337123,  0.39161468,  0.48678179,  3.5270126 ])]
    (t/testing "First 8 L-moments"
      (t/are [order r] (m/delta-eq r (sut/l-moment ys order {:sorted? true}))
        0 1.0
        1 -0.01412282
        2 0.94063132
        3 -0.00167452
        4 0.27196273
        5 0.04125875
        6 0.0800777
        7 0.02874381
        8 0.04678928))
    (t/testing "L-moment ratios"
      (t/are [order r] (m/delta-eq r (sut/l-moment ys order {:sorted? true :ratio? true}))
        3 -0.0017802
        4 0.28912787
        5 0.04386283
        6 0.08513187
        7 0.030558
        8 0.04974242))
    (t/testing "First 8 L-moments, trimmed=1"
      (t/are [order r] (m/delta-eq r (sut/l-moment ys order {:sorted? true :s 1 :t 1}) 1.0e-8)
        0 1.00000000e+00
        1 -1.24483014e-02
        2 4.01201153e-01
        3 -2.04444131e-02,
        4 7.99520964e-02
        5 4.77843216e-03
        6 1.19496896e-02
        7 4.17554371e-03,
        8 -2.80236176e-04))
    (t/testing "First 8 L-moments, left trimmed=1"
      (t/are [order r] (m/delta-eq r (sut/l-moment ys order {:sorted? true :s 1 :t 0}))
        0 1.0
        1 0.9265085
        2 0.7042176
        3 0.18019214
        4 0.19576343,
        5 0.07280187
        6 0.06347921
        7 0.04316176
        8 0.03563686))
    (t/testing "First 8 L-moments, right trimmed=1"
      (t/are [order r] (m/delta-eq r (sut/l-moment ys order {:sorted? true :s 0 :t 1}))
        0 1.0
        1 -0.95475414
        2 0.70672938
        3 -0.18242483
        4 0.14418999,
        5 -0.02329137
        6 0.02994477
        7 -0.0103117
        8 0.01700108))
    (t/testing "L-variation"
      (t/is (m/delta-eq -66.603658638321 (sut/l-variation ys) 1.0e-5)))))

;; reference values from R: uniroot() applied to the expectile-defining equation
;; tau*sum(max(x-t,0)) == (1-tau)*sum(max(t-x,0))

(t/deftest expectile
  (let [ys [-1 0.24 0.61 1 1.88 11.86 29.46 84.07 164.82 247.68]]
    (t/is (m/delta-eq (sut/expectile ys 0.5) (sut/mean ys)))
    (t/is (m/delta-eq 25.9 (sut/expectile ys 0.25)))
    (t/is (m/delta-eq 97.54428571428572 (sut/expectile ys 0.75)))
    (t/is (m/delta-eq 11.2492 (sut/expectile ys 0.1)))
    (t/is (m/delta-eq 147.71615385 (sut/expectile ys 0.9) 1.0e-5)
          "R uniroot() tol vs. fastmath's own solver differ past ~1e-6, both converging to the same root")))

;; moment
;; reference values from R on mtcars$drat: order=2 default == population variance (already
;; verified, see Variance & Dispersion group); order=2 normalized == 1.0 (self-consistency:
;; a moment normalized by its own scale); order=3/4 normalized == the classical "g1"/"g2"
;; population-moment skewness/kurtosis coefficients mean((x-mu)^3)/popvar^1.5,
;; mean((x-mu)^4)/popvar^2; :absolute?/:center/:mean?=false == sum(abs(x-center)^order)

(t/deftest moment-test
  (let [d (mtcars :drat)]
    (t/is (m/delta= (sut/population-variance d) (sut/moment d)) "order=2.0 default")
    (t/is (m/delta= (sut/population-variance d) (sut/moment d 2.0)))
    ;; bug fixed this session: :normalize? was dividing by the sample (n-1) variance instead
    ;; of the population variance used everywhere else in this function, silently biasing
    ;; every :normalize? result by a hidden (n-1)/n factor
    (t/is (m/delta-eq 1.0 (sut/moment d 2.0 {:normalize? true}))
          "a moment normalized by its own scale is ~1.0")
    (t/is (m/delta-eq 0.2788734320 (sut/moment d 3.0 {:normalize? true})) "g1, population-moment skewness")
    (t/is (m/delta-eq 2.4351160973 (sut/moment d 4.0 {:normalize? true})) "g2, population-moment kurtosis (non-excess)")
    (t/is (m/delta-eq 220.0094882255 (sut/moment d 1.5 {:absolute? true :center 0.0 :mean? false})))))

;; winsor/trim/trim-lower/trim-upper
;; reference values from R on mtcars$mpg: quantile(mpg, c(0.2,0.5,0.8), type=6) = q20/q50/q80,
;; then pmin/pmax (winsor) or direct filtering (trim*) against those bounds; NaN case verified
;; via R's ifelse(is.na(x), median, x) + filtering on the non-NaN quantile bounds

(t/deftest trimming-winsorizing-test
  (let [mpg (mtcars :mpg)]
    (t/is (m/delta= 625.36 (sut/sum (sut/winsor mpg))))
    (t/is (= 32 (count (sut/winsor mpg))) "winsor never changes length")
    (t/is (m/delta= 384.4 (sut/sum (sut/trim mpg))))
    (t/is (= 20 (count (sut/trim mpg))))
    (t/is (m/delta= 564.8 (sut/sum (sut/trim-lower mpg))))
    (t/is (= 26 (count (sut/trim-lower mpg))))
    ;; bug fixed this session: trim-upper used the `quantile`-th percentile as its cutoff
    ;; instead of the documented `(1.0-quantile)`-th, making it discard the bottom N% and
    ;; keep the top (1-N)%, the exact opposite of "discard values above the 80th percentile"
    (t/is (m/delta= 462.5 (sut/sum (sut/trim-upper mpg))))
    (t/is (= 26 (count (sut/trim-upper mpg)))))
  (t/testing "NaN handling: kept, replaced by median of the non-NaN values"
    (let [x [1.0 2 3 ##NaN 5 100]]
      (t/is (v/delta-eq [1.2 2.0 3.0 3.0 5.0 81.0] (sut/winsor x)))
      (t/is (= [2.0 3.0 3.0 5.0] (sut/trim x))))))

;; skewness/kurtosis, all estimator variants
;; reference values from R on mtcars$drat:
;; - :G1/:g1/:pearson/:b1  <-> e1071::skewness(x, type=2/1/1/3)
;; - :G2/:g2/:excess/:kurt/:b2 <-> e1071::kurtosis(x, type=2/1/1) with :kurt = type1+3
;; - :skew (BCa) <-> e1071::skewness(type=2) * (n-2)/(n*sqrt(n-1))
;; - :bowley/:yule/[:yule u] <-> (q(1-u)-2*q(0.5)+q(u))/(q(1-u)-q(u)), quantile(type=6)
;; - :B3 <-> (mean-median)/mean(abs(x-median))
;; - :mode <-> (mean-DescTools::Mode(x))/sd(x); :median <-> 3*(mean-median)/sd(x)
;; - :moors <-> classical octile formula, quantile(type=6)
;; - :crow <-> Crow-Siddiqui formula, quantile(type=6) + qnorm(alpha)/qnorm(beta)
;; - :geary <-> mean(abs(x-mean(x)))/sqrt(mean((x-mean(x))^2)) (mad/popstddev, both
;;   already-verified primitives, see Variance & Dispersion group)
;; - :l-skewness/:l-kurtosis <-> l-moment(x,3/4,{:ratio? true}) (already R-verified, see
;;   Moments & L-moments group) -- checked here for correct delegation only
;; - :hogg (both fns) <-> depends on trim-lower/trim-upper (see Trimming & Winsorizing group);
;;   **bug found and fixed this session**: hogg-skewness/hogg-kurtosis were calling
;;   trim-upper with the wrong quantile parameter (alpha/beta instead of 1.0-alpha/1.0-beta) --
;;   this was previously masked by trim-upper's own now-fixed bug (the two bugs cancelled
;;   out); fixing trim-upper alone (Trimming & Winsorizing group) broke :hogg here as a
;;   direct side effect, caught and fixed in this same pass

(t/deftest skewness-kurtosis-test
  (let [d (mtcars :drat)]
    (t/testing "skewness"
      (t/are [typ res] (m/delta-eq res (sut/skewness d typ))
        :G1 0.2927802132
        :g1 0.2788734320
        :pearson 0.2788734320
        :b1 0.2659039046
        :skew 0.0492983237
        :bowley -0.4642857143
        :B3 -0.2207428171
        :hogg 1.6607629428
        :mode 0.9848203500
        :median -0.5523176444
        :l-skewness 0.0452120451)
      (t/is (m/delta-eq 0.2927802132 (sut/skewness d)) "default == :G1")
      (t/is (m/delta-eq -0.4642857143 (sut/skewness d [:yule 0.25])))
      (t/is (m/delta-eq -0.1725768322 (sut/skewness d [:yule 0.1])))
      (t/is (m/delta-eq 0.0452120451 (sut/skewness d :l-skewness) 1.0e-8)
            "delegates to l-moment order 3 ratio")
      (t/is (= (sut/l-moment d 3 {:ratio? true}) (sut/skewness d :l-skewness))))
    (t/testing "kurtosis"
      (t/are [typ res] (m/delta-eq res (sut/kurtosis d typ))
        :G2 -0.4504324511
        :g2 -0.5648839027
        :excess -0.5648839027
        :b2 -0.7147006157
        :kurt 2.4351160973
        :geary 0.8612546038
        :moors -0.5975833333
        :crow -0.3225136183
        :hogg -0.1519236160
        :l-kurtosis 0.0352669499)
      (t/is (m/delta-eq -0.4504324511 (sut/kurtosis d)) "default == :G2")
      (t/is (m/delta-eq -0.3225136183 (sut/kurtosis d [:crow 0.025 0.25])))
      (t/is (= (sut/l-moment d 4 {:ratio? true}) (sut/kurtosis d :l-kurtosis))))))

;; wmode(s)/mode(s): default (:default) method
;; reference values from R: DescTools::Mode(c(4,1,4,2,4,2)) == 4; weighted totals via
;; aggregate(ws, by=list(ys), FUN=sum) -> value 2 has the highest total weight (3.5)

(t/deftest wmode
  (let [ys [4 1 4 2 4 2]
        ys2 [:a :b :a :c :a :c]
        ws [1, 3, 0.5, 1.5, 1, 2]]
    (t/is (== 4 (sut/wmode ys)))
    (t/is (== 4 (sut/mode ys)))
    (t/is (== 2 (sut/wmode ys ws)))
    (t/is (= :a (sut/wmode ys2)))
    (t/is (= :c (sut/wmode ys2 ws)))
    (t/is (= '(4) (sut/wmodes ys)))
    (t/is (= '(2) (sut/wmodes ys ws)))
    (t/is (= [1.0 2.0] (sut/modes [1 1 2 2 3])) "tie: R DescTools::Mode(c(1,1,2,2,3)) == c(1,2)")))

;; mode(s): :pearson, :histogram, :kde estimation methods
;; reference values from R on y = c(2,3,3,4,4,4,5,5,6,20):
;; :pearson  -> 3*median(y)-2*mean(y)
;; :histogram -> rice-rule bins (k=ceiling(2*n^(1/3))=5) + the classical grouped-data mode
;;   formula L+h*(fm-f1)/(2fm-f1-f2), matching fastmath's own histogram bin edges/counts
;; :kde -> fastmath's :nrd bandwidth (Silverman's rule, scale 1.06) uses Apache Commons'
;;   :legacy/type-6 IQR internally, NOT R's own bw.nrd() (which uses IQR type=7) -- confirmed
;;   by hand-deriving R's quantile(y, type=6) reproduces fastmath's bandwidth exactly
;;   (1.1230099004 both ways); this is a genuine, correctly-implemented convention difference,
;;   not a bug. The KDE argmax-among-data-points logic itself is then verified independently
;;   via R's density(y, bw=<that confirmed bandwidth>, kernel="gaussian").

(t/deftest mode-estimation-methods-test
  (let [y [2.0 3 3 4 4 4 5 5 6 20]]
    (t/is (m/delta-eq 0.8 (sut/mode y :pearson)))
    (t/is (= [0.8000000000000007] (sut/modes y :pearson)))
    (t/is (m/delta-eq 3.92 (sut/mode y :histogram)))
    (t/is (m/delta-eq 4.0 (sut/mode y :kde)))))

;; comparison with R test
(t/deftest ks-two-samples
  (t/testing "distinct values"
    (let [x (range -100 101)
          y [21.125 22.925 21.525 18.825 18.225 14.425 24.525 19.325 17.925 16.525 17.425 15.325 10.525 14.825 32.525 30.525 34.025 21.625 15.625 13.425 27.425 26.125 15.925 19.825 15.125]]
      (t/testing "exact method"
        (let [res (sut/ks-test-two-samples x y)
              res+ (sut/ks-test-two-samples x y {:sides :right})
              res- (sut/ks-test-two-samples x y {:sides :left})]
          (t/is (m/delta-eq (:p-value res) 6.675949e-07 1.0e-12 1.0e-12))
          (t/is (m/delta-eq (:stat res) 0.5522388))
          (t/is (m/delta-eq (:p-value res+) 3.337974e-07 1.0e-12 1.0e-12))
          (t/is (m/delta-eq (:stat res+) 0.5522388))
          (t/is (m/delta-eq (:p-value res-) 0.005897501))
          (t/is (m/delta-eq (:stat res-) 0.3283582))))
      (t/testing "approximate method"
        (let [res (sut/ks-test-two-samples x y {:method :approximate})
              res+ (sut/ks-test-two-samples x y {:sides :right :method :approximate})
              res- (sut/ks-test-two-samples x y {:sides :left :method :approximate})]
          (t/is (m/delta-eq (:p-value res) 2.57807e-06 1.0e-12 1.0e-12))
          (t/is (m/delta-eq (:d res) 0.5522388))
          (t/is (m/delta-eq (:p-value res+)  1.289035e-06 1.0e-12 1.0e-12))
          (t/is (m/delta-eq (:dp res+) 0.5522388))
          (t/is (m/delta-eq (:p-value res-) 0.008274217))
          (t/is (m/delta-eq (:dn res-) 0.3283582))))))
  (t/testing "ties"
    (let [{:keys [p-value stat]} (sut/ks-test-two-samples [1 1 1 1 1 20] [1 1 1 -1])]
      (t/is (m/delta-eq 0.6666667 p-value))
      (t/is (m/delta-eq 0.25 stat)))
    (let [{:keys [p-value stat]} (sut/ks-test-two-samples [1 1 1 -1] [1 1 1 1 1 20])]
      (t/is (m/delta-eq 0.6666667 p-value))
      (t/is (m/delta-eq 0.25 stat)))))

(t/deftest ks-one-sample
  (t/testing "normality"
    (t/are [side distr s pv]
        (let [{:keys [stat p-value]} (sut/ks-test-one-sample (iris :petal-width)
                                                             (r/distribution distr)
                                                             {:sides side})]
          (and (m/delta-eq stat s)
               (m/delta-eq p-value pv 1.0e-6 1.0e-14)))
      :both :normal 0.5686175 2.981394e-07
      :left :normal 0.5686175 1.490697e-07
      :right :normal 0.006209665 0.9929283))
  (t/testing "cauchy"
    (t/are [side distr s pv]
        (let [{:keys [stat p-value]} (sut/ks-test-one-sample (iris :petal-width)
                                                             (r/distribution distr)
                                                             {:sides side})]
          (and (m/delta-eq stat s)
               (m/delta-eq p-value pv 1.0e-6 1.0e-14)))
      :both :cauchy 0.5317255 2.469265e-06
      :left :cauchy 0.5317255  1.234632e-06
      :right :cauchy 0.1211189 0.4858579)))

(t/deftest auc
  (let [labels [0 1 0 1 1]
        scores [0 1 1 2 3]
        m (sut/binary-measures-thr labels scores)]
    (t/is (m/delta-eq 0.9166666667 (sut/auc m) 1.0e-10))
    (t/is (m/delta-eq 0.9513888889 (sut/auc m :recall :precision) 1.0e-10))
    (t/is (m/delta-eq 0.3585624955 (sut/auc m (bm/->f-beta 2) :recall) 1.0e-10))
    (t/is (m/delta-eq 2.25 (sut/auc m :tp [1 1 0.5 0.5 0.5]) 1.0e-10))
    (t/testing "U based auc-roc"
      (t/is (m/delta-eq 0.9166666667 (sut/auc-roc labels scores) 1.0e-10)))))

(t/deftest multiclass-measures
  (t/are [res m] (= res (sut/multiclass-measure [0, 1, 2, 1, 1, 2] [0, 1, 1, 0, 0, 2] m))
    {0 0.5, 1 0.4, 2 0.6666666666666666} {:average nil}
    {0 0.6666666666666666, 1 0.5, 2 0.8333333333333334} {:metric :accuracy :average nil}
    {0 2.6632016632016633, 1 2.1683991683991684, 2 1.1683991683991684} {:metric :f-inv-beta :beta 0.45 :average nil})
  (t/are [res m] (m/delta-eq res (sut/multiclass-measure [0, 1, 2, 1, 1, 2] [0, 1, 1, 0, 0, 2] m) 1.0e-10)
    0.5222222222222223 {}
    0.5 {:average sut/harmean}
    0.5108729549290354 {:average sut/geomean}
    0.5055555555555555 {:weighted? true}
    0.6111111111111112 {:metric :recall}
    0.3777777777777777 {:metric :mk}
    0.25 {:metric :mk :average :micro}
    0.5641764963265006 {:metric (bm/->f-beta 0.45)}
    0.5641764963265006 {:metric :f-beta :beta 0.45}
    2.0 {:metric :f-inv-beta :beta 0.45}
    0.5 {:average :micro}
    1.0 {:metric :tp}))

;; python

(t/deftest multiclass-auc
  (let [cls [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
             2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
             2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        probs (zipmap [0 1 2] (apply map vector [[9.81571915e-01, 1.84280700e-02, 1.45178280e-08],
                                                 [9.71303935e-01, 2.86960349e-02, 3.02480322e-08],
                                                 [9.85256375e-01, 1.47436122e-02, 1.23637595e-08],
                                                 [9.76035540e-01, 2.39644199e-02, 3.97833127e-08],
                                                 [9.85224691e-01, 1.47752968e-02, 1.20202360e-08],
                                                 [9.70237188e-01, 2.97627377e-02, 7.39922282e-08],
                                                 [9.86761074e-01, 1.32389059e-02, 2.00173574e-08],
                                                 [9.76132059e-01, 2.38679136e-02, 2.77705876e-08],
                                                 [9.79594414e-01, 2.04055553e-02, 3.06806396e-08],
                                                 [9.68730674e-01, 3.12692937e-02, 3.17879211e-08],
                                                 [9.76228194e-01, 2.37717864e-02, 1.93851435e-08],
                                                 [9.75192202e-01, 2.48077541e-02, 4.40211987e-08],
                                                 [9.74194503e-01, 2.58054754e-02, 2.15639157e-08],
                                                 [9.91855741e-01, 8.14425510e-03, 3.90224543e-09],
                                                 [9.88016494e-01, 1.19835030e-02, 2.84388421e-09],
                                                 [9.86664404e-01, 1.33355827e-02, 1.29289141e-08],
                                                 [9.87967825e-01, 1.20321655e-02, 9.26622416e-09],
                                                 [9.81326796e-01, 1.86731838e-02, 1.97939406e-08],
                                                 [9.56133678e-01, 4.38662530e-02, 6.89813101e-08],
                                                 [9.83979353e-01, 1.60206265e-02, 2.06537553e-08],
                                                 [9.46162653e-01, 5.38372604e-02, 8.69361054e-08],
                                                 [9.81572481e-01, 1.84274864e-02, 3.30183499e-08],
                                                 [9.95955798e-01, 4.04420107e-03, 1.31942828e-09],
                                                 [9.51852640e-01, 4.81471223e-02, 2.37432063e-07],
                                                 [9.51599615e-01, 4.84001781e-02, 2.06841197e-07],
                                                 [9.50935312e-01, 4.90646012e-02, 8.71140079e-08],
                                                 [9.69326068e-01, 3.06738459e-02, 8.66087146e-08],
                                                 [9.74634012e-01, 2.53659631e-02, 2.51118716e-08],
                                                 [9.77037142e-01, 2.29628404e-02, 1.75183184e-08],
                                                 [9.70949888e-01, 2.90500534e-02, 5.86690765e-08],
                                                 [9.63897392e-01, 3.61025375e-02, 7.06065611e-08],
                                                 [9.64477230e-01, 3.55227119e-02, 5.78003442e-08],
                                                 [9.88293432e-01, 1.17065609e-02, 7.09086352e-09],
                                                 [9.88961024e-01, 1.10389705e-02, 5.35397576e-09],
                                                 [9.68320259e-01, 3.16796972e-02, 4.33328420e-08],
                                                 [9.84418938e-01, 1.55810539e-02, 8.03169211e-09],
                                                 [9.78616140e-01, 2.13838501e-02, 9.70819839e-09],
                                                 [9.86732169e-01, 1.32678227e-02, 8.55691000e-09],
                                                 [9.85673700e-01, 1.43262842e-02, 1.55572690e-08],
                                                 [9.73798331e-01, 2.62016402e-02, 2.85802792e-08],
                                                 [9.86459575e-01, 1.35404140e-02, 1.14219407e-08],
                                                 [9.61558525e-01, 3.84414086e-02, 6.60696542e-08],
                                                 [9.88894343e-01, 1.11056459e-02, 1.13028599e-08],
                                                 [9.72218536e-01, 2.77813260e-02, 1.37484231e-07],
                                                 [9.60092902e-01, 3.99068746e-02, 2.23467087e-07],
                                                 [9.73508856e-01, 2.64911043e-02, 4.00774441e-08],
                                                 [9.80164438e-01, 1.98355369e-02, 2.54747511e-08],
                                                 [9.83156866e-01, 1.68431138e-02, 2.01948304e-08],
                                                 [9.78350311e-01, 2.16496704e-02, 1.88317688e-08],
                                                 [9.78396413e-01, 2.16035675e-02, 1.93702474e-08],
                                                 [2.13987432e-03, 8.73951959e-01, 1.23908167e-01],
                                                 [5.85216912e-03, 8.59695656e-01, 1.34452175e-01],
                                                 [1.06758790e-03, 7.25097759e-01, 2.73834653e-01],
                                                 [1.54552953e-02, 9.39649064e-01, 4.48956405e-02],
                                                 [2.38989912e-03, 8.15215488e-01, 1.82394613e-01],
                                                 [7.01568413e-03, 8.60085940e-01, 1.32898376e-01],
                                                 [3.79562357e-03, 7.16470249e-01, 2.79734128e-01],
                                                 [1.47908548e-01, 8.49036974e-01, 3.05447806e-03],
                                                 [2.79169950e-03, 8.96639557e-01, 1.00568743e-01],
                                                 [4.16030884e-02, 9.11777354e-01, 4.66195574e-02],
                                                 [5.58421379e-02, 9.37689054e-01, 6.46880852e-03],
                                                 [1.52500349e-02, 8.98756062e-01, 8.59939028e-02],
                                                 [9.11131071e-03, 9.76561002e-01, 1.43276872e-02],
                                                 [3.06586048e-03, 7.79254339e-01, 2.17679800e-01],
                                                 [7.46682507e-02, 9.14880929e-01, 1.04508203e-02],
                                                 [5.31148888e-03, 9.26323353e-01, 6.83651586e-02],
                                                 [8.75944810e-03, 7.74608924e-01, 2.16631628e-01],
                                                 [1.64929247e-02, 9.65142381e-01, 1.83646947e-02],
                                                 [1.81849807e-03, 8.01061537e-01, 1.97119965e-01],
                                                 [2.40187127e-02, 9.59384631e-01, 1.65966560e-02],
                                                 [2.32188990e-03, 4.40086504e-01, 5.57591606e-01],
                                                 [1.68885300e-02, 9.56721430e-01, 2.63900402e-02],
                                                 [7.19118786e-04, 5.96336663e-01, 4.02944218e-01],
                                                 [3.05498792e-03, 8.59895356e-01, 1.37049656e-01],
                                                 [7.10773506e-03, 9.42926499e-01, 4.99657660e-02],
                                                 [5.10304023e-03, 9.20078724e-01, 7.48182358e-02],
                                                 [1.12624118e-03, 8.01558033e-01, 1.97315726e-01],
                                                 [5.82649343e-04, 4.81159124e-01, 5.18258227e-01],
                                                 [5.51186421e-03, 8.13058474e-01, 1.81429662e-01],
                                                 [6.17716388e-02, 9.34849170e-01, 3.37919145e-03],
                                                 [2.92013412e-02, 9.57187455e-01, 1.36112039e-02],
                                                 [3.72110271e-02, 9.55257338e-01, 7.53163510e-03],
                                                 [2.52496079e-02, 9.56412697e-01, 1.83376948e-02],
                                                 [4.52127608e-04, 3.49726620e-01, 6.49821253e-01],
                                                 [1.02808960e-02, 7.50809181e-01, 2.38909923e-01],
                                                 [1.00508595e-02, 7.88600958e-01, 2.01348183e-01],
                                                 [2.27821616e-03, 8.05162838e-01, 1.92558946e-01],
                                                 [2.77362712e-03, 9.13185469e-01, 8.40409035e-02],
                                                 [2.71657101e-02, 9.28408022e-01, 4.44262680e-02],
                                                 [1.99679531e-02, 9.38027271e-01, 4.20047761e-02],
                                                 [8.74954409e-03, 8.97865098e-01, 9.33853579e-02],
                                                 [4.67010088e-03, 8.28230808e-01, 1.67099091e-01],
                                                 [1.76280668e-02, 9.56967313e-01, 2.54046202e-02],
                                                 [1.21786612e-01, 8.75164475e-01, 3.04891349e-03],
                                                 [1.45196867e-02, 9.20414079e-01, 6.50662344e-02],
                                                 [2.00671541e-02, 9.38012634e-01, 4.19202121e-02],
                                                 [1.71686071e-02, 9.25354578e-01, 5.74768153e-02],
                                                 [8.53517335e-03, 9.35086554e-01, 5.63782730e-02],
                                                 [2.43745722e-01, 7.54963149e-01, 1.29112920e-03],
                                                 [1.92301930e-02, 9.35969076e-01, 4.48007311e-02],
                                                 [9.12432594e-07, 3.91398404e-03, 9.96085104e-01],
                                                 [2.44409256e-04, 1.62562032e-01, 8.37193559e-01],
                                                 [2.50456633e-06, 2.55801238e-02, 9.74417372e-01],
                                                 [3.14684752e-05, 8.16928482e-02, 9.18275683e-01],
                                                 [3.76935037e-06, 1.74447068e-02, 9.82551524e-01],
                                                 [5.59246256e-08, 4.64120216e-03, 9.95358742e-01],
                                                 [5.79752413e-03, 5.13669097e-01, 4.80533378e-01],
                                                 [6.27420734e-07, 2.13546452e-02, 9.78644727e-01],
                                                 [5.26779981e-06, 5.33042449e-02, 9.46690487e-01],
                                                 [6.61332577e-07, 5.74260730e-03, 9.94256731e-01],
                                                 [3.04457419e-04, 2.10442697e-01, 7.89252845e-01],
                                                 [7.32666069e-05, 1.37317523e-01, 8.62609210e-01],
                                                 [2.14474813e-05, 6.52667101e-02, 9.34711842e-01],
                                                 [2.30768730e-04, 1.45237851e-01, 8.54531380e-01],
                                                 [6.95555419e-05, 4.34975598e-02, 9.56432885e-01],
                                                 [5.20470570e-05, 5.40241368e-02, 9.45923816e-01],
                                                 [5.60392603e-05, 1.22915359e-01, 8.77028602e-01],
                                                 [8.56887541e-08, 3.56143008e-03, 9.96438484e-01],
                                                 [3.18606600e-09, 1.00100928e-03, 9.98998988e-01],
                                                 [3.91566006e-04, 4.52005093e-01, 5.47603341e-01],
                                                 [5.65329877e-06, 2.38481898e-02, 9.76146157e-01],
                                                 [6.17556478e-04, 1.90422288e-01, 8.08960156e-01],
                                                 [3.17040922e-08, 4.65751504e-03, 9.95342453e-01],
                                                 [5.89152668e-04, 3.93054898e-01, 6.06355950e-01],
                                                 [1.29504371e-05, 3.86041120e-02, 9.61382938e-01],
                                                 [4.90165820e-06, 5.15012509e-02, 9.48493847e-01],
                                                 [1.07899404e-03, 4.56445696e-01, 5.42475310e-01],
                                                 [1.02969477e-03, 3.85315989e-01, 6.13654316e-01],
                                                 [1.07399358e-05, 3.63419125e-02, 9.63647348e-01],
                                                 [1.70281075e-05, 1.42022336e-01, 8.57960636e-01],
                                                 [1.07347467e-06, 2.92100364e-02, 9.70788890e-01],
                                                 [7.11029023e-07, 1.74223056e-02, 9.82576983e-01],
                                                 [7.94942894e-06, 2.72639644e-02, 9.72728086e-01],
                                                 [5.32052441e-04, 4.75589879e-01, 5.23878068e-01],
                                                 [6.30714650e-05, 1.88649594e-01, 8.11287335e-01],
                                                 [3.97353603e-07, 1.17478837e-02, 9.88251719e-01],
                                                 [1.17294484e-05, 1.73241394e-02, 9.82664131e-01],
                                                 [6.81993651e-05, 1.19491575e-01, 8.80440225e-01],
                                                 [1.63122392e-03, 4.40321653e-01, 5.58047123e-01],
                                                 [4.00086159e-05, 9.34804514e-02, 9.06479540e-01],
                                                 [6.36270339e-06, 2.02880236e-02, 9.79705614e-01],
                                                 [1.00474356e-04, 1.20620465e-01, 8.79279060e-01],
                                                 [2.44409256e-04, 1.62562032e-01, 8.37193559e-01],
                                                 [2.06886328e-06, 1.25871870e-02, 9.87410744e-01],
                                                 [3.84717412e-06, 1.21026066e-02, 9.87893546e-01],
                                                 [5.63724539e-05, 8.00990354e-02, 9.19844592e-01],
                                                 [2.28496543e-04, 2.51916686e-01, 7.47854818e-01],
                                                 [1.39414414e-04, 1.57115825e-01, 8.42744760e-01],
                                                 [4.60438469e-05, 3.84191847e-02, 9.61534771e-01],
                                                 [4.78871489e-04, 2.34864482e-01, 7.64656646e-01]]))]
    (t/are [res m] (m/delta-eq res (sut/multiclass-auc cls probs m) 1.0e-10)
      0.9991777777777777 {:average :micro}
      0.9983333333333334 {}
      0.0016666666666666 {:metric :det}
      0.9967842720318594 {:metric :pr}
      0.1778127938399307 {:metric [(bm/->f-beta 2) :recall]})))

;; pooled-variance / pooled-stddev / pooled-mad
;; reference: 3 unequal-size groups (n=5,7,4), R manual formulas; `:unbiased` cross-checked
;; against `summary(aov(y~g))$"Mean Sq"[2]` (ANOVA within-group mean square == classic
;; pooled variance). `pooled-mad`'s default constant (1.4826022185056023 = 1/qnorm(0.75))
;; cross-checked against R `mad(resid, constant=1/qnorm(0.75))`; R's own `mad()` default
;; constant 1.4826 is a rounded approximation, so `mad(resid)` alone differs in the 5th decimal
;; (not a bug).
(def pooled-g1 [2.0 4.0 6.0 8.0 10.0])
(def pooled-g2 [1.0 3.0 3.0 5.0 7.0 9.0 11.0])
(def pooled-g3 [10.0 20.0 15.0 25.0])
(def pooled-groups [pooled-g1 pooled-g2 pooled-g3])

(t/deftest pooled-statistics-test
  (t/testing "pooled-variance"
    (t/is (m/delta-eq 18.67032967032967 (sut/pooled-variance pooled-groups)))
    (t/is (m/delta-eq 18.67032967032967 (sut/pooled-variance pooled-groups :unbiased)))
    (t/is (m/delta-eq 15.169642857142858 (sut/pooled-variance pooled-groups :biased)))
    (t/is (m/delta-eq 21.53968253968254 (sut/pooled-variance pooled-groups :avg))))
  (t/testing "pooled-stddev"
    (t/is (m/delta-eq 4.32091768844648 (sut/pooled-stddev pooled-groups)))
    (t/is (m/delta-eq 3.8948225706882793 (sut/pooled-stddev pooled-groups :biased)))
    (t/is (m/delta-eq 4.641086353396426 (sut/pooled-stddev pooled-groups :avg))))
  (t/testing "pooled-mad"
    (t/is (m/delta-eq 3.7065055462640055 (sut/pooled-mad pooled-groups)))
    (t/is (m/delta-eq 2.5 (sut/pooled-mad pooled-groups 1.0)))))
