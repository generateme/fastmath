(ns fastmath.stats.bins-test
  (:require [fastmath.stats.bins :as sut]
            [fastmath.core :as m]
            [fastmath.grid :as grid]
            [clojure.test :as t]
            [clojure.edn :as edn]
            [clojure.java.io :as io]))

;; Group 1: counting estimators (sturges, rice, doane)
;; reference: Python numpy.histogram_bin_edges(x, bins=...) and R e1071::skewness(type=1)
;; for doane's raw formula (both independently re-derived, not trusted from
;; fastmath's own output). Two datasets: n100 = 1..100 (synthetic), mpg = mtcars$mpg
;; (n=32, matching the fixture already used elsewhere in this repo's test suite).
;; numpy: sturges(n100)=8, rice(n100)=10, doane(n100)=8
;;        sturges(mpg)=6,  rice(mpg)=7,   doane(mpg)=8

(def n100 (double-array (range 1 101)))
(def mpg (double-array [21.0 21.0 22.8 21.4 18.7 18.1 14.3 24.4 22.8 19.2 17.8 16.4 17.3 15.2
                         10.4 10.4 14.7 32.4 30.4 33.9 21.5 15.5 15.2 13.3 19.2 27.3 26.0 30.4
                         15.8 19.7 15.0 21.4]))

(t/deftest sturges-test
  (t/are [n bins] (= bins (sut/sturges n))
    100 8
    32  6
    1   1
    2   2)
  (t/testing "n < 1: log2 is not finite, guarded to 1 (confirmed bug: used to throw)"
    (t/is (= 1 (sut/sturges 0)))
    (t/is (= 1 (sut/sturges -5)))))

(t/deftest rice-test
  (t/are [n bins] (= bins (sut/rice n))
    100 10
    32  7
    1   2
    0   1))

;; sqrt / terrell-scott: added post-audit at user request (both known, simple,
;; data-independent bin-count rules, per the same numpy.histogram_bin_edges oracle
;; used throughout this audit for sturges/rice; terrell-scott has no numpy/R builtin,
;; verified against its well-known closed form ceil(cbrt(2n)) directly).
;; numpy: sqrt(n100)=10, sqrt(mpg)=6

(t/deftest sqrt-test
  (t/are [n bins] (= bins (sut/sqrt n))
    100 10
    32  6
    1   1
    0   1
    -5  1))

(t/deftest terrell-scott-test
  (t/are [n bins] (= bins (sut/terrell-scott n))
    100 6
    32  4
    1   2
    0   1
    -5  1))

(t/deftest doane-test
  (t/testing "n < 3: skewness undefined, always 1"
    (t/is (= 1 (sut/doane (double-array [1.0]) 1)))
    (t/is (= 1 (sut/doane (double-array [1.0 2.0]) 2))))
  (t/testing "ceiling, not truncation, of the raw formula value -- the confirmed bug"
    (t/is (= 8 (sut/doane n100 100)) "raw value 7.6438 -> ceil 8, not trunc 7")
    (t/is (= 8 (sut/doane mpg 32)) "raw value ~7.39-7.43 -> ceil 8, not trunc 7 (R e1071::skewness(type=1) confirms)")))

;; Group 2: width-based estimators (scott, freedman-diaconis, scott-fd-helper)
;; reference: numpy.histogram_bin_edges(x, bins='scott'/'fd')
;; numpy: scott(n100)=5, fd(n100)=5; scott(mpg)=4, fd(mpg)=5

(t/deftest scott-test
  (t/are [avs n bins] (= bins (sut/scott avs n))
    n100 100 5
    mpg  32  4))

;; scott-2d: multivariate normal-reference rule (Scott 1992, p.82), b_x=c*stddev(xs)*n^-1/4,
;; b_y=c*stddev(ys)*n^-1/4, c=(48*pi)^(1/4) (the exact 2d constant), area=b_x*b_y. Hand-computed independently against the same
;; formula for a fixed synthetic dataset (not against fastmath's own output).
(t/deftest scott-2d-test
  (let [xs (double-array [1.0 2.0 3.0 4.0 5.0 1.5 2.5 3.5])
        ys (double-array [2.0 3.0 1.0 5.0 2.0 4.0 0.5 3.5])
        n (alength xs)
        sx (Math/sqrt (org.apache.commons.math3.stat.StatUtils/variance xs))
        sy (Math/sqrt (org.apache.commons.math3.stat.StatUtils/variance ys))
        factor (Math/pow n -0.25)
        c (Math/pow (* 48.0 Math/PI) 0.25)
        expected (* c sx factor c sy factor)]
    (t/is (m/delta-eq expected (sut/scott-2d xs ys))))
  (t/testing "degenerate input: no crash, propagates non-positive/non-finite gracefully"
    (t/is (not (pos? (sut/scott-2d (double-array []) (double-array [])))))
    (t/is (zero? (sut/scott-2d (double-array [1.0]) (double-array [2.0]))))
    (t/is (zero? (sut/scott-2d (double-array [1.0 1.0 1.0]) (double-array [2.0 2.0 2.0]))) "two constant axes")))

(t/deftest freedman-diaconis-test
  (t/are [avs n bins] (= bins (sut/freedman-diaconis avs n))
    n100 100 5
    mpg  32  5))

(t/deftest scott-fd-degenerate-fallback-test
  (t/testing "confirmed bug: heavily-tied-but-not-fully-constant data (IQR=0, MAD=0,
              range<>0) used to crash with IllegalArgumentException: Value out of
              range for long: Infinity -- now falls back to 1 bin instead of crashing"
    (let [degen (double-array [1 1 1 1 1 1 1 1 1 2])]
      (t/is (pos? (sut/scott degen 10)))
      (t/is (pos? (sut/freedman-diaconis degen 10)))))
  (t/testing "fully-constant data: still 1, now via an explicit guard rather than an
              accidental NaN->0 long-cast"
    (let [const-data (double-array (repeat 10 5.0))]
      (t/is (= 1 (sut/scott const-data 10)))
      (t/is (= 1 (sut/freedman-diaconis const-data 10)))))
  (t/testing "n < 1 / empty data (Unknown C, confirmed for scott/freedman-diaconis):
              used to throw NullPointerException-style errors on (first empty-array);
              now guarded to 1, matching doane's own n<3 early-return pattern"
    (t/is (= 1 (sut/scott (double-array []) 0)))
    (t/is (= 1 (sut/freedman-diaconis (double-array []) 0)))))

;; ## 2d normal-reference rules and the data-driven estimators
;; Reference: utils/fastmath/dev/generate_bins_reference.py (numpy 2.5.2, percentiles 'weibull' =
;; Commons Math StatUtils/percentile), an independent implementation of the same specification,
;; plus numpy `stone` (leave-one-out cross-validation) for the 1d bin counts.

(def ^:private reference
  (delay (-> "stats/bins_reference.edn" io/resource slurp edn/read-string)))

(defn- arr [xs] (double-array xs))

(t/deftest normal-reference-2d-test
  (doseq [[set-name {:keys [xs ys scott scott-correlated freedman-diaconis]}] (:two-d @reference)
          :let [xs (arr xs) ys (arr ys)]]
    (t/is (m/delta-eq scott (sut/scott-2d xs ys)) (str set-name))
    (t/is (m/delta-eq scott-correlated (sut/scott-correlated xs ys)) (str set-name))
    (t/is (m/delta-eq freedman-diaconis (sut/freedman-diaconis-2d xs ys)) (str set-name)))
  (t/testing "scott-correlated equals scott-2d for uncorrelated axes with equal spreads"
    (let [xs (arr [1.0 -1.0 0.0 0.0 2.0 -2.0 0.0 0.0])
          ys (arr [0.0 0.0 1.0 -1.0 0.0 0.0 2.0 -2.0])]
      (t/is (m/delta-eq (sut/scott-2d xs ys) (sut/scott-correlated xs ys)))))
  (t/testing "correlation gives smaller cells"
    (let [{:keys [xs ys]} (:correlated (:two-d @reference))
          xs (arr xs) ys (arr ys)]
      (t/is (< (sut/scott-correlated xs ys) (* 0.5 (sut/scott-2d xs ys))))))
  (t/testing "degenerate input: not finite or not positive, no exception"
    (doseq [f [sut/scott-2d sut/scott-correlated sut/freedman-diaconis-2d]
            [xs ys] [[[] []] [[1.0] [2.0]] [[1.0 1.0 1.0] [2.0 2.0 2.0]]]]
      (t/is (not (pos? (f (arr xs) (arr ys)))) (str xs ys))))
  (t/testing "nearly correlated is not a line"
    (let [xs (arr [1.0 2.0 3.0 4.0]) ys (arr [2.0 4.0 6.0 8.001])]
      (t/is (< (sut/scott-correlated xs ys) (* 1.0e-3 (#'sut/scott-line-area 7.0 4)))))))

(defn- line-width
  "1d bin width `c * spread * n^(-1/3)` (closed form for points on a line)."
  ^double [c spread n]
  (* c spread (Math/pow n (/ -1.0 3.0))))

(defn- sample-variance [xs] (let [m (/ (reduce + xs) (count xs))] (/ (reduce + (map #(let [d (- % m)] (* d d)) xs)) (dec (count xs)))))

(t/deftest line-branch-test
  (let [t (map double (range 1000))
        n 1000]
    (t/testing "scott-correlated: any line, 1d Scott width along it"
      (doseq [[label xs ys] [["diagonal" t (map #(* 2.0 %) t)]
                             ["anti-diagonal" t (map #(- 3.0 (* 0.1 %)) t)]
                             ["axis" t (repeat n 5.0)]
                             ["small scale" (map #(* 1.0e-6 %) t) (map #(* 2.0e-6 %) t)]]]
        (let [w (line-width 3.5 (Math/sqrt (+ (sample-variance xs) (sample-variance ys))) n)]
          (t/is (m/delta-eq (* w w) (sut/scott-correlated (arr xs) (arr ys)) 0.0 1.0e-9) label))))
    (t/testing "scott-2d and freedman-diaconis-2d: a line parallel to an axis, 1d widths of the other axis"
      (let [ys (repeat n 5.0)
            iqr (- (org.apache.commons.math3.stat.StatUtils/percentile (arr t) 75.0)
                   (org.apache.commons.math3.stat.StatUtils/percentile (arr t) 25.0))
            ws (line-width 3.5 (Math/sqrt (sample-variance t)) n)
            wf (line-width 2.0 iqr n)]
        (t/is (m/delta-eq (* ws ws) (sut/scott-2d (arr t) (arr ys)) 0.0 1.0e-9))
        (t/is (m/delta-eq (* ws ws) (sut/scott-2d (arr ys) (arr t)) 0.0 1.0e-9) "constant x")
        (t/is (m/delta-eq (* wf wf) (sut/freedman-diaconis-2d (arr t) (arr ys)) 0.0 1.0e-9))
        (t/is (m/delta-eq (* wf wf) (sut/freedman-diaconis-2d (arr ys) (arr t)) 0.0 1.0e-9) "constant x")))
    (t/testing "data-driven: a line parallel to an axis gives the line rule, any type"
      (doseq [f [sut/cross-validation-2d sut/shimazaki-shinomoto-2d]
              type grid/cell-names]
        (t/is (m/delta-eq (sut/scott-correlated (arr t) (arr (repeat n 5.0)))
                          (f (arr t) (arr (repeat n 5.0)) type)) (str type))))
    (t/testing "data-driven: a diagonal line is searched (finite, positive, not the line rule)"
      (let [xs (arr t) ys (arr (map #(* 2.0 %) t))]
        (doseq [f [sut/cross-validation-2d sut/shimazaki-shinomoto-2d]]
          (let [a (f xs ys :square)]
            (t/is (and (m/valid-double? a) (pos? a)))
            (t/is (not (m/delta-eq a (sut/scott-correlated xs ys))))))))))

(t/deftest options-validation-test
  (let [{:keys [xs ys]} (:mixture (:two-d @reference))
        xs (arr xs) ys (arr ys)
        option-of (fn [f] (try (f) ::no-exception (catch clojure.lang.ExceptionInfo e (:option (ex-data e)))))]
    (doseq [[k v] [[:shifts 0] [:shifts -2] [:shifts 1.5] [:shifts "2"] [:shifts 2N]
                   [:candidates 0] [:candidates -5] [:candidates 2.0]]
            f [#(sut/cross-validation-bins xs {k v}) #(sut/shimazaki-shinomoto-bins xs {k v})
               #(sut/cross-validation-2d xs ys :square {k v}) #(sut/shimazaki-shinomoto-2d xs ys :square {k v})]]
      (t/is (= k (option-of f)) (str k " " (pr-str v))))
    (t/testing "boundary valid values and nil"
      (doseq [o [{:shifts 1} {:candidates 1} {:shifts nil :candidates nil} {:shifts (int 2)}]]
        (t/is (pos? (sut/cross-validation-bins xs o)) (pr-str o))
        (t/is (pos? (sut/shimazaki-shinomoto-2d xs ys :square o)) (pr-str o))))))

(t/deftest data-driven-1d-reference-test
  (doseq [[set-name {:keys [data candidates cv-costs ss-costs cv ss cv-shifts-5 ss-shifts-5 stone stone-in-range]}]
          (:one-d @reference)
          :let [avs (arr data)
                label (str set-name)]]
    (t/is (= candidates (#'sut/bins-candidates avs (alength avs) (- (apply max data) (apply min data)) 200)) label)
    (let [sample (fastmath.stats.bins.Sample1d. avs (apply min data) (apply max data) 1)]
      (t/is (every? true? (map #(m/delta-eq %1 (#'sut/cost-1d :cross-validation sample %2) 1.0e-9 1.0e-12) cv-costs candidates)) label)
      (t/is (every? true? (map #(m/delta-eq %1 (#'sut/cost-1d :shimazaki-shinomoto sample %2) 1.0e-9 1.0e-12) ss-costs candidates)) label))
    (t/is (= cv (sut/cross-validation-bins avs)) label)
    (t/is (= ss (sut/shimazaki-shinomoto-bins avs)) label)
    (t/is (= cv-shifts-5 (sut/cross-validation-bins avs {:shifts 5})) label)
    (t/is (= ss-shifts-5 (sut/shimazaki-shinomoto-bins avs {:shifts 5})) label)
    (t/is (= cv (sut/cross-validation-bins avs {:shifts 1 :candidates 200})) (str label " explicit defaults"))
    (when stone-in-range
      (t/is (= stone (sut/cross-validation-bins avs)) (str label " numpy stone")))))

(t/deftest data-driven-2d-reference-test
  ;; :square cells only, the other cell types are checked in fastmath.stats-test against histogram2d counts
  (doseq [[set-name {:keys [xs ys] :as ref}] (:two-d @reference)
          [f k opts] [[sut/cross-validation-2d :cv nil]
                      [sut/cross-validation-2d :cv-sv {:sv [0.3 0.7]}]
                      [sut/cross-validation-2d :cv-shifts-3 {:shifts 3}]
                      [sut/shimazaki-shinomoto-2d :ss nil]
                      [sut/shimazaki-shinomoto-2d :ss-sv {:sv [0.3 0.7]}]
                      [sut/shimazaki-shinomoto-2d :ss-shifts-3 {:shifts 3}]]]
    (t/is (m/delta-eq (ref k) (f (arr xs) (arr ys) :square opts) 1.0e-9 1.0e-12) (str set-name " " k))))

(t/deftest data-driven-candidates-test
  (t/testing "1d: log-spaced counts between a quarter and twenty times Freedman-Diaconis, distinct, ascending"
    (let [avs (arr (:data (:bimodal (:one-d @reference))))
          fd (sut/freedman-diaconis avs (alength avs))
          cs (#'sut/bins-candidates avs (alength avs) (- (apply max avs) (apply min avs)) 200)]
      (t/is (= (max 1 (quot fd 4)) (first cs)))
      (t/is (= (min (* 20 fd) (alength avs)) (peek cs)))
      (t/is (apply < cs))
      (t/is (<= (count cs) 200))
      (t/is (= [(first cs)] (#'sut/bins-candidates avs (alength avs) (- (apply max avs) (apply min avs)) 1)) "one candidate")))
  (t/testing "2d: sizes between a twentieth and four times the Freedman-Diaconis size, descending"
    (let [{:keys [xs ys]} (:mixture (:two-d @reference))
          xs (arr xs) ys (arr ys)
          base (grid/area->size :flat-hex (sut/freedman-diaconis-2d xs ys))
          cs (#'sut/grid-candidates xs ys :flat-hex 50)]
      (t/is (= 50 (count cs)))
      (t/is (m/delta-eq (* 4.0 base) (first cs)))
      (t/is (m/delta-eq (/ base 20.0) (last cs)))
      (t/is (apply > cs)))))

(t/deftest data-driven-discretized-test
  (t/testing "1d: at most range / (2 * smallest gap) bins"
    (let [{:keys [data]} (:integers (:one-d @reference))
          limit (Math/floor (/ (- (apply max data) (apply min data)) 2.0))]
      (doseq [f [sut/cross-validation-bins sut/shimazaki-shinomoto-bins]]
        (t/is (<= (f (arr data)) limit)))))
  (t/testing "2d: area at least 4 * gap_x * gap_y"
    (let [{:keys [xs ys]} (:integers (:two-d @reference))]
      (doseq [f [sut/cross-validation-2d sut/shimazaki-shinomoto-2d]
              type grid/cell-names]
        (t/is (>= (* (f (arr xs) (arr ys) type) (+ 1.0 1.0e-12)) 4.0) (str type))))))

(t/deftest data-driven-ties-test
  ;; equally spaced values give equal costs for several bin counts: the chosen count is the
  ;; smallest one with the lowest cost
  (let [avs (arr [0.0 1.0 2.0 3.0])
        sample (fastmath.stats.bins.Sample1d. avs 0.0 3.0 1)
        cands (#'sut/bins-candidates avs 4 3.0 200)
        costs (map #(#'sut/cost-1d :cross-validation sample %) cands)
        best (apply min costs)]
    (t/is (= (first (keep (fn [[c cost]] (when (== cost best) c)) (map vector cands costs)))
             (sut/cross-validation-bins avs)))))

(t/deftest data-driven-properties-test
  (let [{:keys [data]} (:bimodal (:one-d @reference))
        {:keys [xs ys]} (:mixture (:two-d @reference))]
    (t/testing "1d: scaling and translation do not change the bin count"
      (doseq [f [sut/cross-validation-bins sut/shimazaki-shinomoto-bins]]
        (t/is (= (f (arr data)) (f (arr (map #(* 4.0 %) data))) (f (arr (map #(+ 1024.0 %) data)))))))
    (t/testing "2d: scaling by 4 multiplies the area by 16, translation with :sv keeps it"
      (doseq [f [sut/cross-validation-2d sut/shimazaki-shinomoto-2d]
              type grid/cell-names]
        (let [a (f (arr xs) (arr ys) type)]
          (t/is (m/delta-eq (* 16.0 a) (f (arr (map #(* 4.0 %) xs)) (arr (map #(* 4.0 %) ys)) type) 1.0e-9 1.0e-9) (str type))
          (t/is (m/delta-eq a (f (arr (map #(+ 8.0 %) xs)) (arr (map #(+ 8.0 %) ys)) type {:sv [8.0 8.0]}) 1.0e-9 1.0e-9) (str type)))))
    (t/testing "closed forms: scaling by 4 multiplies the area by 16, translation keeps it"
      (doseq [f [sut/scott-2d sut/scott-correlated sut/freedman-diaconis-2d]]
        (t/is (m/delta-eq (* 16.0 (f (arr xs) (arr ys))) (f (arr (map #(* 4.0 %) xs)) (arr (map #(* 4.0 %) ys)))))
        (t/is (m/delta-eq (f (arr xs) (arr ys)) (f (arr (map #(+ 8.0 %) xs)) (arr (map #(+ 8.0 %) ys)))))))))

(t/deftest data-driven-degenerate-test
  (doseq [f [sut/cross-validation-bins sut/shimazaki-shinomoto-bins]
          avs [[] [1.0] [2.0 2.0 2.0]]]
    (t/is (= 1 (f (arr avs))) (str avs)))
  (doseq [f [sut/cross-validation-2d sut/shimazaki-shinomoto-2d]
          [xs ys] [[[] []] [[1.0] [2.0]] [[1.0 1.0 1.0] [2.0 2.0 2.0]]]]
    (t/is (m/nan? (f (arr xs) (arr ys) :square)) (str xs ys))))

