(ns fastmath.stats.bins-test
  (:require [fastmath.stats.bins :as sut]
            [fastmath.core :as m]
            [clojure.test :as t]))

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

;; scott-2d: multivariate normal-reference rule (Scott 1992, p.82), b_x=3.5*stddev(xs)*n^-1/4,
;; b_y=3.5*stddev(ys)*n^-1/4, area=b_x*b_y. Hand-computed independently against the same
;; formula for a fixed synthetic dataset (not against fastmath's own output).
(t/deftest scott-2d-test
  (let [xs (double-array [1.0 2.0 3.0 4.0 5.0 1.5 2.5 3.5])
        ys (double-array [2.0 3.0 1.0 5.0 2.0 4.0 0.5 3.5])
        n (alength xs)
        sx (Math/sqrt (org.apache.commons.math3.stat.StatUtils/variance xs))
        sy (Math/sqrt (org.apache.commons.math3.stat.StatUtils/variance ys))
        factor (Math/pow n -0.25)
        expected (* 3.5 sx factor 3.5 sy factor)]
    (t/is (m/delta-eq expected (sut/scott-2d xs ys))))
  (t/testing "degenerate input: no crash, propagates non-positive/non-finite gracefully"
    (t/is (not (pos? (sut/scott-2d (double-array []) (double-array [])))))
    (t/is (zero? (sut/scott-2d (double-array [1.0]) (double-array [2.0]))))
    (t/is (zero? (sut/scott-2d (double-array [1.0 1.0 1.0]) (double-array [1.0 2.0 3.0]))))))

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
