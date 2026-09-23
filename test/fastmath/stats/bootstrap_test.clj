(ns fastmath.stats.bootstrap-test
  (:require [fastmath.stats.bootstrap :as sut]
            [fastmath.stats :as stats]
            [fastmath.random :as r]
            [fastmath.core :as m]
            [clojure.test :as t]))

;; Group 1: jackknife family

(t/deftest jackknife-test
  (t/testing "leave-one-out: n samples of size n-1, i-th sample excludes i-th element"
    (let [vs [10 20 30 40 50]
          out (sut/jackknife vs)]
      (t/is (= 5 (count out)))
      (t/is (every? #(= 4 (count %)) out))
      (t/is (= [[20 30 40 50] [10 30 40 50] [10 20 40 50]
                 [10 20 30 50] [10 20 30 40]] (map vec out))))))

(t/deftest jackknife-plus-test
  (t/testing "add-one: n samples of size n+1, i-th sample duplicates i-th element at the end"
    (let [vs [10 20 30 40 50]
          out (sut/jackknife+ vs)]
      (t/is (= 5 (count out)))
      (t/is (every? #(= 6 (count %)) out))
      (t/is (= [[10 20 30 40 50 10] [10 20 30 40 50 20] [10 20 30 40 50 30]
                 [10 20 30 40 50 40] [10 20 30 40 50 50]] (map vec out))))))

;; Group 2: generation & aggregation
;; `aircondit` dataset from R `boot` package (also used in notebooks/bootstrap.clj)

(def aircondit [3 5 7 18 43 85 91 98 100 130 230 487])
(def aircondit-set (set (map double aircondit)))

(t/deftest bootstrap-nonparametric-test
  (t/testing "default: :samples 500, size == input count, every value drawn from original data"
    (let [{:keys [samples data model]} (sut/bootstrap aircondit nil {:samples 50 :rng (r/rng :jvm 42)})]
      (t/is (= 50 (count samples)))
      (t/is (every? #(= 12 (count %)) samples))
      (t/is (every? (fn [s] (every? aircondit-set s)) samples))
      (t/is (= aircondit data))
      (t/is (r/distribution? model))))
  (t/testing ":size forces a different sample size"
    (let [{:keys [samples]} (sut/bootstrap aircondit nil {:samples 10 :size 5 :rng (r/rng :jvm 42)})]
      (t/is (every? #(= 5 (count %)) samples)))))

(t/deftest bootstrap-jackknife-method-test
  (t/testing ":method :jackknife delegates to jackknife, no :model built"
    (let [{:keys [samples] :as res} (sut/bootstrap aircondit nil {:method :jackknife})]
      (t/is (not (contains? res :model)))
      (t/is (= (sut/jackknife aircondit) samples))))
  (t/testing ":method :jackknife+ delegates to jackknife+, no :model built"
    (let [{:keys [samples] :as res} (sut/bootstrap aircondit nil {:method :jackknife+})]
      (t/is (not (contains? res :model)))
      (t/is (= (sut/jackknife+ aircondit) samples)))))

(t/deftest bootstrap-smoothing-test
  (t/testing ":gaussian smoothing adds noise -- resampled values need not be in the original data"
    (let [{:keys [samples]} (sut/bootstrap aircondit nil {:samples 20 :smoothing :gaussian :rng (r/rng :jvm 1)})]
      (t/is (every? #(= 12 (count %)) samples))
      (t/is (not (every? (fn [s] (every? aircondit-set s)) samples)))))
  (t/testing ":kde smoothing builds a continuous model -- resampled values need not be in the original data"
    (let [{:keys [samples]} (sut/bootstrap aircondit nil {:samples 20 :smoothing :kde :rng (r/rng :jvm 1)})]
      (t/is (every? #(= 12 (count %)) samples))
      (t/is (not (every? (fn [s] (every? aircondit-set s)) samples))))))

(t/deftest bootstrap-distribution-test
  (t/testing ":integer-discrete-distribution for integer data"
    (let [idata [1 2 3 4 5 6 7 8]
          {:keys [samples]} (sut/bootstrap idata nil {:samples 10 :distribution :integer-discrete-distribution :rng (r/rng :jvm 1)})]
      (t/is (every? #(= 8 (count %)) samples))
      (t/is (every? (fn [s] (every? (set idata) s)) samples))))
  (t/testing ":categorical-distribution for non-numeric data, also the auto-selected default"
    (let [cdata [:a :b :c :a :b]
          explicit (sut/bootstrap cdata nil {:samples 10 :distribution :categorical-distribution :rng (r/rng :jvm 1)})
          auto (sut/bootstrap cdata nil {:samples 10 :rng (r/rng :jvm 1)})]
      (t/is (= (:samples explicit) (:samples auto)) "same seed, same result -> confirms auto-dispatch to categorical")
      (t/is (every? (fn [s] (every? (set cdata) s)) (:samples explicit))))))

(t/deftest bootstrap-parametric-test
  (t/testing "explicit :model, any fastmath.random distribution"
    (let [model (r/distribution :exponential {:mean (stats/mean aircondit)})
          {:keys [samples data model]} (sut/bootstrap {:data aircondit :model model} nil {:samples 10 :rng (r/rng :jvm 3)})]
      (t/is (= 10 (count samples)))
      (t/is (every? #(= 12 (count %)) samples))
      (t/is (= aircondit data))
      (t/is (r/distribution? model)))))

(t/deftest bootstrap-multi-test
  (t/testing ":dimensions :multi samples each dimension independently, default per-dimension size"
    (let [g0 [21 21 22.8 21.4 18.7]
          g1 [26 30.4 15.8 19.7 15 21.4]
          g0-set (set (map double g0))
          g1-set (set (map double g1))
          {:keys [samples]} (sut/bootstrap [g0 g1] nil {:dimensions :multi :samples 8 :rng (r/rng :jvm 1)})]
      (t/is (= 8 (count samples)))
      (t/is (every? (fn [s] (= [5 6] (map count s))) samples))
      (t/is (every? (fn [s] (every? g0-set (first s))) samples))
      (t/is (every? (fn [s] (every? g1-set (second s))) samples)))))

(t/deftest bootstrap-antithetic-test
  (t/testing ":antithetic? true still yields the requested number/size of samples (parametric only)"
    (let [model (r/distribution :exponential {:mean (stats/mean aircondit)})
          {:keys [samples]} (sut/bootstrap {:data aircondit :model model} nil
                                            {:samples 7 :antithetic? true :rng (r/rng :jvm 5)})]
      (t/is (= 7 (count samples)))
      (t/is (every? #(= 12 (count %)) samples)))))

(t/deftest bootstrap-include-test
  (t/testing ":include? true adds the original dataset as one of the returned samples"
    (let [{:keys [samples]} (sut/bootstrap aircondit nil {:samples 3 :include? true :rng (r/rng :jvm 1)})]
      (t/is (= 4 (count samples)))
      (t/is (some #(= aircondit (vec %)) samples)))))

(t/deftest bootstrap-statistic-test
  (t/testing "when a statistic is given, bootstrap delegates to bootstrap-stats"
    (let [{:keys [t0 ts mean median variance stddev bias sem]}
          (sut/bootstrap aircondit stats/mean {:samples 200 :rng (r/rng :jvm 7)})
          ats (double-array ts)]
      (t/is (m/delta-eq t0 (stats/mean aircondit)))
      (t/is (= 200 (count ts)))
      (t/is (m/delta-eq mean (stats/mean ats)))
      (t/is (m/delta-eq median (stats/median ats)))
      (t/is (m/delta-eq variance (stats/variance ats)))
      (t/is (m/delta-eq stddev (m/sqrt (stats/variance ats))))
      (t/is (m/delta-eq bias (- (stats/mean ats) t0)))
      (t/is (m/delta-eq sem (/ (m/sqrt (stats/variance ats)) (m/sqrt (count aircondit))))))))

(t/deftest bootstrap-stats-test
  (t/testing "bootstrap-stats applied directly to a :data/:samples map (e.g. from jackknife+)"
    (let [data [1 2 3 4 5 10 100]
          {:keys [t0 ts mean variance]}
          (sut/bootstrap-stats {:data data :samples (sut/jackknife+ data)} stats/stddev)
          ats (double-array ts)]
      (t/is (m/delta-eq t0 (stats/stddev data)))
      (t/is (= 7 (count ts)))
      (t/is (m/delta-eq mean (stats/mean ats)))
      (t/is (m/delta-eq variance (stats/variance ats))))))

;; Group 3: simple CIs
;; reference: R `boot` package, same `aircondit` dataset, `mean` statistic, R=20000
;; replicates, independently seeded from fastmath's run -- exact agreement isn't
;; expected (two independent Monte Carlo runs), tolerance ~2% of the mean (~2.2)
;; #+begin_src R
;; b <- boot(aircondit$hours, function(d,i) mean(d[i]), R=20000)  # set.seed(123)
;; boot.ci(b, type=c("norm","basic","perc"))
;; #+end_src
;; Normal: (34.2, 182.0); Basic: (25.6, 170.1); Percentile: (46.1, 190.6)

(def ^:private aircondit-boot20k
  (delay (sut/bootstrap aircondit stats/mean {:samples 20000 :rng (r/rng :jvm 123)})))

(t/deftest ci-normal-test
  (let [[lo hi t0] (sut/ci-normal @aircondit-boot20k)]
    (t/is (m/delta-eq 108.08333333333333 t0) "t0 is deterministic, exact")
    (t/is (m/delta-eq 34.2 lo 2.0) "R: 34.2")
    (t/is (m/delta-eq 182.0 hi 2.0) "R: 182.0")))

(t/deftest ci-basic-test
  (let [[lo hi t0] (sut/ci-basic @aircondit-boot20k)]
    (t/is (m/delta-eq 108.08333333333333 t0))
    (t/is (m/delta-eq 25.6 lo 2.0) "R: 25.6")
    (t/is (m/delta-eq 170.1 hi 2.0) "R: 170.1")))

(t/deftest ci-percentile-test
  (let [[lo hi t0] (sut/ci-percentile @aircondit-boot20k)]
    (t/is (m/delta-eq 108.08333333333333 t0))
    (t/is (m/delta-eq 46.1 lo 2.0) "R: 46.1")
    (t/is (m/delta-eq 190.6 hi 2.0) "R: 190.6")))

;; Group 4: t-based CIs

(t/deftest ci-t-test
  (t/testing "deterministic given :ts -- independent re-derivation of the qt() multiplier"
    (let [[lo hi t0] (sut/ci-t @aircondit-boot20k)
          ^double stddev (:stddev @aircondit-boot20k)
          ;; R: qt(0.975, 19999) = 1.960083 (df = count(ts)-1 = 19999)
          merr (* stddev 1.960083)]
      (t/is (m/delta-eq 108.08333333333333 t0))
      (t/is (m/delta-eq (- t0 merr) lo 1.0e-3))
      (t/is (m/delta-eq (+ t0 merr) hi 1.0e-3)))))

;; reference: R `boot` package, `aircondit$hours`, statistic returning [mean, var/n]
;; (as `boot.ci(type="stud")` requires), R=20000, independently seeded from fastmath.
;; Studentized/bootstrap-t intervals are known to have higher Monte Carlo variance than
;; the other CI types (division by small per-resample stddev estimates on skewed data);
;; empirically confirmed here by re-running fastmath's own ci-studentized across 3 seeds,
;; observing a ~4-5 unit spread at the upper bound alone -- tolerance widened to 6.0
;; accordingly (vs. 2.0 for Group 3's CI types).
;; #+begin_src R
;; b <- boot(x, function(d,i) c(mean(d[i]), var(d[i])/length(i)), R=20000)
;; boot.ci(b, type="stud")
;; #+end_src
;; -> (46.8, 290.0)

(t/deftest ci-studentized-test
  (let [[lo hi t0] (sut/ci-studentized @aircondit-boot20k)]
    (t/is (m/delta-eq 108.08333333333333 t0))
    (t/is (m/delta-eq 46.8 lo 6.0) "R: 46.8")
    (t/is (m/delta-eq 290.0 hi 6.0) "R: 290.0")))

;; Group 5: BC/BCa CIs
;; `ci-bc` and `ci-bca`'s ts-skewness path are deterministic given `:ts`/`:t0`/`acc` --
;; verified exact against a from-scratch R reimplementation of the Efron & Tibshirani
;; BC/BCa formula (same approach the Stats Audit used for the sibling
;; percentile-bc-extent/percentile-bca-extent, R's coxed::bca being unavailable),
;; :legacy <-> R quantile type=6 (established correspondence, Stats Audit).
;; #+begin_src R
;; bca_ci <- function(t0, ts, alpha, acc=0, qtype=6) {
;;   z0 <- qnorm(mean(ts < t0)); p <- c(alpha/2, 1-alpha/2); zp <- qnorm(p)
;;   num <- z0 + zp; denom <- if (acc==0) 1.0 else (1 - acc*num)
;;   adj_p <- pnorm(z0 + num/denom)
;;   as.numeric(quantile(ts, probs=adj_p, type=qtype))
;; }
;; #+end_src
;; on fastmath's own 500-replicate `ts` (t0=108.08333333333333):
;;   BC (acc=0): (49.25683, 198.47975)
;;   BCa (acc=-0.003847763633842666, fastmath's own skewness(ts,:skew)/-6): (48.54447, 198.25057)

(def ^:private aircondit-boot500
  (delay (sut/bootstrap aircondit stats/mean {:samples 500 :rng (r/rng :jvm 321)})))

(t/deftest ci-bc-test
  (let [[lo hi t0] (sut/ci-bc @aircondit-boot500)]
    (t/is (m/delta-eq 108.08333333333333 t0))
    (t/is (m/delta-eq 49.25683362322621 lo 1.0e-6))
    (t/is (m/delta-eq 198.47975272000716 hi 1.0e-6))))

(t/deftest ci-bca-ts-skewness-path-test
  (t/testing "acceleration from ts skewness when :data/:statistic absent"
    (let [[lo hi t0] (sut/ci-bca (dissoc @aircondit-boot500 :data))]
      (t/is (m/delta-eq 108.08333333333333 t0))
      (t/is (m/delta-eq 48.544468288281884 lo 1.0e-6))
      (t/is (m/delta-eq 198.25057130568084 hi 1.0e-6)))))

(t/deftest ci-bca-jackknife-path-test
  (t/testing "default: acceleration from jackknife replicates of :data/:statistic"
    ;; reference: R boot::boot.ci(type="bca") on aircondit/mean, R=20000,
    ;; independently seeded (set.seed(789)) -> (57.5, 223.8)
    (let [[lo hi t0] (sut/ci-bca @aircondit-boot20k)]
      (t/is (m/delta-eq 108.08333333333333 t0))
      (t/is (m/delta-eq 57.5 lo 3.0) "R: 57.5")
      (t/is (m/delta-eq 223.8 hi 3.0) "R: 223.8"))))

(t/deftest ci-bca-multi-dim-fallback-test
  (t/testing "multi-dimensional :data no longer crashes ci-bca's jackknife-acceleration
              path (regression for the confirmed :dimensions :multi bug); falls back to,
              and matches exactly, the ts-skewness path"
    (let [g0 [21 21 22.8 21.4 18.7 18.1 14.3 24.4 22.8 19.2 17.8]
          g1 [26 30.4 15.8 19.7 15 21.4 15.5 30.4 21.5 32.4 27.3]
          rom (fn [[a b]] (- (stats/mean b) (stats/mean a)))
          bs (sut/bootstrap [g0 g1] rom {:dimensions :multi :rng (r/rng :jvm 1)})]
      (t/is (contains? bs :data) "sanity: :data is present, this is the path that used to crash")
      (t/is (= (sut/ci-bca bs) (sut/ci-bca (dissoc bs :data)))))))

