(ns fastmath.random.truncated-mixture-test
  "Tests of `ccdf`, `truncated` and `mixture` of `fastmath.random.distributions` (Random Audit, work packages 5.4 and 6.1-6.10).

  The oracle is `test/resources/random/truncated_mixture_reference.edn` (utils/fastmath/dev/generate_random_truncated_mixture_reference.py):
  scipy `sf` for `ccdf`, mpmath (60 digits) for the truncated distributions and the mixture quantiles. The seeded outputs of the
  code before the audit are `test/resources/random/baseline_reference.edn`."
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.random.distributions :as distr]
            [fastmath.stats :as stats]))

(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "random/truncated_mixture_reference.edn")))))

(def ^:private baseline
  (delay (edn/read-string (slurp (io/resource "random/baseline_reference.edn")))))

(defn- rel-eq?
  "Relative equality `abs(a - b) <= tol * abs(b)`, for `b = 0` the absolute one."
  [a b tol]
  (m/delta-eq (double a) (double b) 0.0 tol))

(defn- thrown-input?
  "True when `f` throws `ex-info` of the type `:invalid-input` with the given `:property`."
  [property f]
  (try (f) false
       (catch clojure.lang.ExceptionInfo e
         (let [data (ex-data e)]
           (and (= :invalid-input (:type data)) (= property (:property data)))))
       (catch Throwable _ false)))

;; ---------------------------------------------------------------- ccdf

(def ^:private ccdf-closed-form-keys
  [:normal :t :logistic :laplace :cauchy :exponential :weibull :pareto :log-normal :gamma :chi-squared])

(def ^:private ccdf-fallback-keys
  "Keys without the survival function: `ccdf` is `1 - cdf`."
  [:beta :f :gumbel :levy :nakagami :triangular :uniform-real :folded-normal :inverse-gaussian
   :poisson :binomial :geometric :negative-binomial :logarithmic :truncated :mixture])

(defn- accurate-levels
  "The levels at which `ccdf` is expected to be accurate, the SSJ `barF` of the Pearson VI loses digits below 1e-12."
  [k]
  (case k
    :pearson-6 [1.0e-6 1.0e-12]
    [1.0e-6 1.0e-12 1.0e-30 1.0e-100]))

(deftest ccdf-accuracy
  (doseq [{:keys [key params median cdf-at-median points]} (:ccdf @reference)
          :when (not= key :inverse-gaussian)] ;; the inverse Gaussian uses the fallback
    (let [d (r/distribution key params)]
      (testing (str key " " params)
        ;; the parameters of the reference are the parameters of the constructor
        (is (m/delta-eq (r/cdf d median) cdf-at-median 1.0e-12) "parameter mapping")
        (doseq [[level x sf] points
                :when (some #{level} (accurate-levels key))]
          ;; relative error observed: 7e-14 at most (normal, sf 1e-100), the limit is the tolerance of the task
          (is (m/delta-eq (distr/ccdf d x) sf 0.0 1.0e-9) (str key " sf level " level)))))))

(deftest ccdf-closed-forms-cover-the-listed-keys
  (is (= (set ccdf-closed-form-keys)
         (set (filter (set ccdf-closed-form-keys) (map :key (:ccdf @reference)))))
      "every closed form has reference rows")
  (doseq [k ccdf-closed-form-keys
          :let [d (r/distribution k)]]
    (testing (str k)
      (is (m/nan? (distr/ccdf d ##NaN)))
      (is (m/delta-eq (distr/ccdf d ##Inf) 0.0 1.0e-300))
      (is (m/delta-eq (distr/ccdf d ##-Inf) 1.0 1.0e-15))
      (is (m/delta-eq (+ (distr/ccdf d 0.7) (r/cdf d 0.7)) 1.0 1.0e-14)))))

(deftest ccdf-outside-the-support
  (is (== 1.0 (distr/ccdf (r/exponential) -3.0)))
  (is (== 1.0 (distr/ccdf (r/gamma) -3.0)))
  (is (== 1.0 (distr/ccdf (r/chi-squared) 0.0)))
  (is (== 1.0 (distr/ccdf (r/weibull) -1.0)))
  (is (== 1.0 (distr/ccdf (r/log-normal) 0.0)))
  (is (== 1.0 (distr/ccdf (r/pareto {:scale 2.0 :shape 3.0}) 1.5)))
  (is (== 1.0 (distr/ccdf (r/pareto {:scale 2.0 :shape 3.0}) 2.0))))

(deftest ccdf-fallback
  (doseq [k ccdf-fallback-keys
          :let [d (r/distribution k)
                xs (keep (fn [p] (try (double (r/icdf d p)) (catch Throwable _ nil))) [0.1 0.5 0.9])]]
    (testing (str k)
      (is (seq xs))
      (doseq [x xs]
        (is (== (- 1.0 (double (r/cdf d x))) (distr/ccdf d x)) "bitwise 1 - cdf")))))

(deftest ccdf-folded-normal-and-inverse-gaussian
  ;; SSJ `barF` of the folded normal is wrong (negative, -0.575 at 0.7979), the fallback is correct
  (let [d (r/folded-normal)]
    (is (m/delta-eq (distr/ccdf d 0.7979) 0.4249 1.0e-3))
    (is (== (- 1.0 (double (r/cdf d 0.7979))) (distr/ccdf d 0.7979))))
  ;; barF of the inverse Gaussian is `1 - cdf` of the SSJ cdf: not worth more than the fallback
  (let [d (r/inverse-gaussian {:mu 2.0 :lambda 5.0})]
    (is (== (- 1.0 (double (r/cdf d 5.0))) (distr/ccdf d 5.0)))))

(deftest ccdf-ssj-keys-use-bar-f
  ;; the far tail where `1 - cdf` is exactly 0.0
  (let [d (r/erlang {:k 3 :lambda 2.0})
        x 40.0]
    (is (zero? (- 1.0 (double (r/cdf d x)))))
    (is (pos? (distr/ccdf d x)))
    (is (== (.barF ^umontreal.ssj.probdist.ContinuousDistribution (r/source-object d) x) (distr/ccdf d x)))))

(deftest ccdf-normal-tail
  ;; scipy.stats.norm.sf(8.5) = 9.479534822203318e-18, sf(15) = 3.670966298e-51; old `1 - cdf` was 0.0
  (is (rel-eq? (distr/ccdf (r/normal) 8.5) 9.479534822203318e-18 1.0e-12))
  (is (m/delta-eq (distr/ccdf (r/normal) 0.0) 0.5 1.0e-15)))

;; ---------------------------------------------------------------- truncated

(def ^:private truncated-bases
  {:n-8-9 (r/normal) :n-5-6 (r/normal) :n-10-11 (r/normal) :n-m9-m8 (r/normal)
   :n-inf-m3 (r/normal) :n-2-inf (r/normal) :n-m1-2 (r/normal) :n-3-4 (r/normal)
   :exp-30-31 (r/exponential {:mean 1.0}) :cauchy-m5-5 (r/cauchy)
   :beta-2-5-01-04 (r/beta {:alpha 2.0 :beta 5.0}) :gamma-3-2-1-10 (r/gamma {:shape 3.0 :scale 2.0})})

(defn- finite-or-nil [x] (when (Double/isFinite (double x)) x))

(deftest truncated-continuous-accuracy
  (doseq [[k base] truncated-bases
          :let [{:keys [left right xs pdf cdf mean variance ps icdf]} (get-in @reference [:truncated k])
                d (r/truncated {:distr base :left (finite-or-nil left) :right (finite-or-nil right)})]]
    (testing (str k)
      (doseq [[x p c] (map vector xs pdf cdf)]
        ;; observed: pdf rel 1.6e-14, cdf abs 1.9e-15 (N[8,9], N[-9,-8]); the limits are the tolerances of the task
        (is (m/delta-eq (r/pdf d x) p 0.0 1.0e-9) (str "pdf " x))
        (is (m/delta-eq (r/cdf d x) c 1.0e-9) (str "cdf " x)))
      ;; observed mean rel 1.5e-15 at most
      (is (m/delta-eq (r/mean d) mean 1.0e-12 1.0e-9) "mean") ;; the absolute part covers the zero mean of the Cauchy
      ;; observed variance rel 1.6e-14 at most (N[8,9]); the limit of the task is 1e-6
      (is (m/delta-eq (r/variance d) variance 0.0 1.0e-6) "variance")
      (doseq [[p q] (map vector ps icdf)]
        ;; the Cauchy median is 0 up to 1e-41 in the reference: the absolute tolerance covers it
        (is (m/delta-eq (r/icdf d p) q 1.0e-12 1.0e-9) (str "icdf " p))
        (is (<= (m/abs (- (r/cdf d (r/icdf d p)) p)) 1.0e-10) (str "round trip " p))))))

(deftest truncated-normal-upper-tail-anchors
  (let [d (r/truncated {:distr (r/normal) :left 8.0 :right 9.0})]
    (is (m/delta-eq (r/mean d) 8.121189 1.0e-6))
    (is (< 8.0 (r/mean d) 9.0))
    (is (< (r/variance d) 0.0142))
    (is (> (r/variance d) 0.0141))
    (is (m/delta-eq (r/icdf d 1.0e-12) 8.0 1.0e-9))
    (is (<= 8.0 (r/icdf d 1.0e-12)))
    (is (< (r/icdf d (- 1.0 1.0e-12)) 9.0)))
  (testing "the interval [10, 11] has a finite density"
    (let [d (r/truncated {:distr (r/normal) :left 10.0 :right 11.0})]
      (is (m/valid-double? (r/pdf d 10.5)))
      (is (< 10.0 (r/mean d) 10.2)))))

(deftest truncated-icdf-speed
  (let [d (r/truncated {:distr (r/normal) :left 8.0 :right 9.0})
        ps [1.0e-12 1.0e-3 0.5 0.999 (- 1.0 1.0e-12)]
        _ (dotimes [_ 200] (doseq [p ps] (r/icdf d p))) ;; warm up
        times (vec (sort (for [_ (range 101)]
                           (let [t0 (System/nanoTime)]
                             (r/icdf d 0.3)
                             (- (System/nanoTime) t0)))))
        median-ms (/ (double (nth times 50)) 1.0e6)]
    ;; observed 0.03-0.25 ms
    (is (<= median-ms 1.0) (str "median icdf time " median-ms " ms"))))

(deftest truncated-sampler-tail
  (doseq [seed [1 2 3]]
    (let [d (r/truncated {:distr (r/normal {:rng (r/rng :mersenne seed)}) :left 8.0 :right 9.0})
          t0 (System/nanoTime)
          xs (vec (r/->seq d 1000))
          seconds (/ (- (System/nanoTime) t0) 1.0e9)]
      (is (<= seconds 5.0) (str "1000 samples in " seconds " s"))
      (is (every? #(<= 8.0 % 9.0) xs))
      (is (> (:p-value (stats/ks-test-one-sample xs d)) 0.01) (str "KS seed " seed))))
  (testing "reproducible through the seed of the base distribution"
    (let [mk #(r/truncated {:distr (r/normal {:rng (r/rng :mersenne 5)}) :left 8.0 :right 9.0})]
      (is (= (vec (r/->seq (mk) 20)) (vec (r/->seq (mk) 20)))))))

(deftest truncated-sampler-methods
  (testing "rejection sampler for Z >= 1e-3 is unchanged: seeded streams of the code before the audit"
    (let [b (:dist-samples @baseline)
          s (fn [d] (vec (r/->seq d 10)))]
      (is (= (:truncated-n-3-4 b)
             (s (r/truncated {:distr (r/normal {:rng (r/rng :mersenne 3)}) :left 3.0 :right 4.0}))))
      (is (= (:truncated-exp b)
             (s (r/truncated {:distr (r/exponential {:mean 2.0 :rng (r/rng :mersenne 3)}) :left 0.5 :right 3.0}))))
      (is (= (:truncated-n-moderate b)
             (s (r/truncated {:distr (r/normal {:rng (r/rng :mersenne 3)}) :left -1.0 :right 2.0}))))))
  (testing "inverse cdf sampler follows the truncated distribution: Exp[30, 31], N(-inf, -9]"
    (doseq [[d lo hi] [[(r/truncated {:distr (r/exponential {:rng (r/rng :mersenne 8)}) :left 30.0 :right 31.0}) 30.0 31.0]
                       [(r/truncated {:distr (r/normal {:rng (r/rng :mersenne 9)}) :right -9.0}) ##-Inf -9.0]]
          :let [xs (vec (r/->seq d 500))]]
      (is (every? #(<= lo % hi) xs))
      (is (> (:p-value (stats/ks-test-one-sample xs d)) 0.01)))))

(deftest truncated-matches-baseline-for-moderate-intervals
  (let [b (get-in @baseline [:dist-values :truncated-n-moderate])
        d (r/truncated {:distr (r/normal) :left -1.0 :right 2.0})]
    (is (= (:pdf b) (mapv #(r/pdf d %) (:xs b))))
    (is (= (:cdf b) (mapv #(r/cdf d %) (:xs b))))
    (is (== (:mean b) (r/mean d)))
    (is (== (:variance b) (r/variance d)))
    ;; icdf is a root of the accurate cdf now: observed difference 1.1e-15 at most
    (is (every? true? (map #(m/delta-eq %1 (r/icdf d %2) 1.0e-12) (:icdf b) [0.01 0.25 0.5 0.75 0.99])))))

(deftest truncated-validation
  (let [n (r/normal)
        poisson (r/poisson {:p 4.0})]
    (testing "bounds"
      (is (thrown-input? :interval #(distr/truncated n 2 -2)))
      (is (thrown-input? :interval #(distr/truncated n 0 0)) "left = right, continuous")
      (is (thrown-input? :left #(distr/truncated n ##NaN 1)))
      (is (thrown-input? :right #(distr/truncated n 0 ##NaN)))
      (is (thrown-input? :left #(distr/truncated n "a" 1)))
      (is (thrown-input? :right #(distr/truncated n 0 :b))))
    (testing "interval outside of the support"
      (is (thrown-input? :interval #(distr/truncated (r/beta) 2 3)))
      (is (thrown-input? :interval #(distr/truncated (r/beta) 1 3)))
      (is (thrown-input? :interval #(distr/truncated (r/exponential) -3 -1)))
      (is (thrown-input? :interval #(distr/truncated (r/uniform-real {:lower 0.0 :upper 1.0}) 1.0 2.0))))
    (testing "zero mass"
      (is (thrown-input? :interval #(distr/truncated n 40 41)) "Z underflows to zero")
      (is (thrown-input? :interval #(distr/truncated poisson 2.2 2.8)) "no integer inside")
      (is (thrown-input? :interval #(distr/truncated poisson -5 -1))))
    (testing "the base"
      (is (thrown-input? :distr #(distr/truncated nil 0 1)))
      (is (thrown-input? :distr #(distr/truncated 5 0 1)))
      (is (thrown-input? :distr #(distr/truncated (r/multi-normal) 0 1)))
      (is (thrown-input? :distr #(distr/truncated (r/enumerated-real {:data [1.0 2.0]}) 0 5)) "discrete on reals"))
    (testing "valid"
      (is (some? (distr/truncated n 0 1)))
      (is (some? (distr/truncated n nil nil)))
      (is (some? (distr/truncated n ##-Inf ##Inf)))
      (is (some? (distr/truncated poisson 3 3)) "discrete single point")
      (is (some? (distr/truncated n 0 1.0e-9)) "tiny interval")
      (is (some? (distr/truncated n 1/2 3))))))

(deftest truncated-bounds-and-identity
  (testing "bounds are the interval intersected with the support"
    (let [d (r/truncated {:distr (r/beta) :left -1 :right 0.5})]
      (is (== 0.0 (r/lower-bound d)))
      (is (== 0.5 (r/upper-bound d)))
      (is (m/delta-eq (r/cdf d 0.25) 0.3125 1.0e-12)))
    (let [d (r/truncated {:distr (r/exponential) :left -5 :right 2})]
      (is (== 0.0 (r/lower-bound d)))
      (is (== 2.0 (r/upper-bound d))))
    (let [d (r/truncated {:distr (r/normal) :left ##-Inf :right ##Inf})]
      (is (== ##-Inf (r/lower-bound d)))
      (is (== ##Inf (r/upper-bound d))))
    (let [d (r/truncated {:distr (r/normal)})]
      (is (== ##-Inf (r/lower-bound d)))
      (is (== ##Inf (r/upper-bound d)))
      (is (m/delta-eq (r/mean d) 0.0 1.0e-12))
      (is (m/delta-eq (r/variance d) 1.0 1.0e-9))
      (is (m/delta-eq (r/icdf d 0.975) 1.959963984540054 1.0e-9))))
  (testing "id and parameters"
    (let [d (r/truncated {:distr (r/normal) :left 0 :right 1})]
      (is (= :truncated (r/distribution-id d)))
      (is (= [:distr :left :right :rng] (r/distribution-parameters d)))
      (is (true? (r/continuous? d)))))
  (testing "outside of the interval"
    (let [d (r/truncated {:distr (r/normal) :left 8.0 :right 9.0})]
      (is (== 0.0 (r/pdf d 7.9)))
      (is (== 0.0 (r/pdf d 9.1)))
      (is (== 0.0 (r/cdf d 7.0)))
      (is (== 1.0 (r/cdf d 10.0)))
      (is (== 8.0 (r/icdf d 0.0)))
      (is (== 9.0 (r/icdf d 1.0)))))
  (testing "truncated of a truncated distribution"
    (let [d (r/truncated {:distr (r/truncated {:distr (r/normal) :left -3 :right 3}) :left 0 :right 1})
          n (r/normal)
          z (- (r/cdf n 1.0) (r/cdf n 0.0))]
      (is (m/delta-eq (r/cdf d 0.5) (/ (- (r/cdf n 0.5) 0.5) z) 1.0e-12)))))

(deftest truncated-set-seed
  (let [mk #(r/truncated {:distr (r/normal {:rng (r/rng :mersenne 1)}) :left 0.0 :right 1.0})
        a (r/set-seed! (mk) 77)
        b (r/set-seed! (mk) 77)]
    (is (= (vec (r/->seq a 30)) (vec (r/->seq b 30))) "set-seed! reseeds the base")))

(deftest truncated-discrete
  (let [bases {:poisson-4-2-6 (r/poisson {:p 4.0}) :poisson-4-0-6 (r/poisson {:p 4.0}) :poisson-4-3-3 (r/poisson {:p 4.0})
               :poisson-4-2-inf (r/poisson {:p 4.0}) :binomial-10-04-3-7 (r/binomial {:trials 10 :p 0.4})
               :geometric-03-2-inf (r/geometric {:p 0.3})}]
    (doseq [[k base] bases
            :let [{:keys [lo hi ks pmf cdf mean variance ps icdf]} (get-in @reference [:discrete k])
                  d (r/truncated {:distr base :left lo :right (finite-or-nil hi)})]]
      (testing (str k)
        (is (false? (r/continuous? d)))
        (is (== lo (r/lower-bound d)))
        (doseq [[x p c] (map vector ks pmf cdf)]
          ;; observed pmf rel 3e-15, cdf abs 3.9e-14 (unbounded Poisson)
          (is (m/delta-eq (r/pdf d x) p 0.0 1.0e-12) (str "pmf " x))
          (is (m/delta-eq (r/pdf d (+ x 0.5)) 0.0 1.0e-300) (str "pmf at a non integer " x))
          (is (m/delta-eq (r/cdf d x) c 1.0e-12) (str "cdf " x)))
        ;; observed mean rel 5e-15, variance rel 3e-13 (the summation of the infinite support stops at 1e-16)
        (is (m/delta-eq (r/mean d) mean 0.0 1.0e-12))
        (is (m/delta-eq (r/variance d) variance 0.0 1.0e-11))
        (is (= icdf (mapv #(long (r/icdf d %)) ps)) "icdf grid")))
    (testing "Poisson(4) on [2, 6]"
      (let [d (r/truncated {:distr (bases :poisson-4-2-6) :left 2 :right 6})]
        (is (m/delta-eq (reduce + (map #(r/pdf d %) (range 2 7))) 1.0 1.0e-15))
        (is (m/delta-eq (r/mean d) 3.8449 1.0e-4))
        (is (m/delta-eq (r/variance d) 1.6739 1.0e-4))
        (is (== 0.0 (r/cdf d 1)))
        (is (== 0.0 (r/cdf d 1.9)))
        (is (== 1.0 (r/cdf d 6)))
        (is (== 2 (r/icdf d 0.0)))
        (is (== 6 (r/icdf d 1.0)))
        (is (== 0.0 (r/pdf d 1)))
        (is (== ##-Inf (r/lpdf d 7)))))
    (testing "binomial(10, 0.4) on [3, 7]"
      (let [d (r/truncated {:distr (bases :binomial-10-04-3-7) :left 3 :right 7})]
        (is (m/delta-eq (r/mean d) 4.410 1.0e-3))
        (is (m/delta-eq (r/variance d) 1.348 1.0e-3))))
    (testing "single point"
      (let [d (r/truncated {:distr (bases :poisson-4-3-3) :left 3 :right 3})]
        (is (== 1.0 (r/pdf d 3)))
        (is (== 0.0 (r/cdf d 2.5)))
        (is (== 1.0 (r/cdf d 3)))
        (is (== 3.0 (r/mean d)))
        (is (== 0.0 (r/variance d)))
        (is (every? #(== 3 %) (r/->seq d 20)))))
    (testing "unbounded right side: the bounds"
      (let [d (r/truncated {:distr (bases :geometric-03-2-inf) :left 2})]
        (is (== 2.0 (r/lower-bound d)))
        (is (> (r/upper-bound d) 1.0e9))
        (is (m/delta-eq (r/pdf d 2) 0.3 1.0e-12))))
    (testing "samples follow the pmf (seeded)"
      (let [d (r/truncated {:distr (r/poisson {:p 4.0 :rng (r/rng :mersenne 3)}) :left 2 :right 6})
            xs (vec (r/->seq d 20000))
            freq (frequencies xs)]
        (is (every? #(<= 2 % 6) xs))
        (is (every? integer? xs))
        ;; the standard error of a frequency is below 0.0035 here
        (is (every? (fn [k] (< (m/abs (- (/ (double (get freq k 0)) 20000.0) (r/pdf d k))) 0.015)) (range 2 7)))))
    (testing "rejection limit: an interval with a negligible probability"
      (let [d (r/truncated {:distr (r/poisson {:p 1.0 :rng (r/rng :mersenne 1)}) :left 30 :right 40})
            ex (try (r/sample d) nil (catch clojure.lang.ExceptionInfo e e))]
        (is (some? ex))
        (is (= :sampling-failure (:type (ex-data ex))))))))

;; ---------------------------------------------------------------- mixture

(deftest mixture-icdf
  (testing "reference quantiles (mpmath)"
    (doseq [[k comps] {:n0-n5 [(r/normal) (r/normal {:mu 5.0})]
                       :n0-n0-1e6 [(r/normal) (r/normal {:sd 1.0e6})]
                       :n0-n100 [(r/normal) (r/normal {:mu 100.0})]}
            :let [{:keys [weights ps icdf]} (get-in @reference [:mixture k])
                  d (r/mixture {:distrs comps :weights weights})]]
      (testing (str k)
        (doseq [[p q] (map vector ps icdf)]
          ;; observed difference 5e-13 (relative to max(1, |q|)); the stop tolerance of the bisection is 1e-12
          (is (m/delta-eq (r/icdf d p) q 1.0e-9 1.0e-9) (str "p " p))))))
  (testing "round trip on the components' tails"
    (let [d (r/mixture {:distrs [(r/normal) (r/normal {:mu 5.0})] :weights [0.3 0.7]})]
      (doseq [p [1.0e-12 0.001 0.25 0.5 0.75 0.999 (- 1.0 1.0e-12)]]
        ;; observed 6.3e-13 at most
        (is (<= (m/abs (- (r/cdf d (r/icdf d p)) p)) 1.0e-10) (str "p " p)))))
  (testing "formerly failing brackets"
    (is (m/delta-eq (r/icdf (r/mixture) 0.5) 0.0 1.0e-9))
    (is (m/delta-eq (r/icdf (r/mixture {:distrs [(r/normal {:mu 3.0})]}) 0.3)
                    (r/icdf (r/normal {:mu 3.0}) 0.3) 1.0e-9) "a single component")
    (is (m/delta-eq (r/icdf (r/mixture {:distrs [(r/normal) (r/normal)] :weights [1 1]}) 0.5) 0.0 1.0e-9) "identical components")
    (is (m/delta-eq (r/icdf (r/mixture {:distrs [(r/normal) (r/normal {:sd 2.0})] :weights [1 1]}) 0.5) 0.0 1.0e-9) "symmetric")
    (is (m/delta-eq (r/icdf (r/mixture {:distrs [(r/normal) (r/normal {:sd 1.0e6})] :weights [1 1]}) 0.5) 0.0 1.0e-9))
    (is (== ##-Inf (r/icdf (r/mixture) 0.0)) "p = 0 is the lower bound")
    (is (== ##Inf (r/icdf (r/mixture) 1.0)) "p = 1 is the upper bound")
    (is (thrown? clojure.lang.ExceptionInfo (r/icdf (r/mixture) ##NaN)) "NaN p is rejected by the icdf proxy (WP 5.3)"))
  (testing "plateau: the smallest x with cdf(x) >= p"
    (let [d (r/mixture {:distrs [(r/normal) (r/normal {:mu 100.0})] :weights [1 1]})
          x (r/icdf d 0.5)]
      (is (< x 10.0))
      (is (m/delta-eq x 8.2924 1.0e-3))
      (is (>= (r/cdf d x) 0.5))
      (is (apply <= (map #(r/icdf d %) (range 0.01 1.0 0.01))) "nondecreasing")))
  (testing "bounded components"
    (let [d (r/mixture {:distrs [(r/uniform-real {:lower 0.0 :upper 1.0}) (r/uniform-real {:lower 5.0 :upper 6.0})] :weights [1 1]})]
      (is (== 0.0 (r/icdf d 0.0)))
      (is (== 6.0 (r/icdf d 1.0)))
      (is (m/delta-eq (r/icdf d 0.25) 0.5 1.0e-9))
      (is (m/delta-eq (r/icdf d 0.75) 5.5 1.0e-9)))))

(deftest mixture-discrete
  (let [{:keys [weights ps icdf]} (get-in @reference [:mixture :pois3-pois10])
        d (r/mixture {:distrs [(r/poisson {:p 3.0}) (r/poisson {:p 10.0})] :weights weights})
        qs (mapv #(r/icdf d %) ps)]
    (is (= icdf qs) "the smallest integer with cdf >= p")
    (is (every? integer? qs))
    (is (integer? (r/sample d)))
    (is (every? integer? (r/->seq d 50)))
    (is (false? (r/continuous? d)))
    (is (true? (r/discrete? d))))
  (testing "mixed discrete and continuous"
    (let [d (r/mixture {:distrs [(r/normal) (r/poisson {:p 3.0})] :weights [1 1]})]
      (is (true? (r/continuous? d)))
      ;; the cdf jumps at the integers: the quantile is the smallest x with cdf(x) >= p
      (let [x (r/icdf d 0.6)]
        (is (>= (r/cdf d x) 0.6))
        (is (< (r/cdf d (- x 1.0e-6)) 0.6))))))

(deftest mixture-validation
  (let [n (r/normal)
        mix (fn [ds ws] #(distr/mixture ds ws nil))]
    (is (thrown-input? :weights (mix [n n] [1])) "fewer weights")
    (is (thrown-input? :weights (mix [n] [1 1])) "more weights")
    (is (thrown-input? :weights (mix [n n] [1 -1])) "negative")
    (is (thrown-input? :weights (mix [n n] [1 ##NaN])) "NaN")
    (is (thrown-input? :weights (mix [n n] [1 ##Inf])) "infinite")
    (is (thrown-input? :weights (mix [n n] [0 0])) "all zero")
    (is (thrown-input? :weights (mix [n] ["a"])) "not a number")
    (is (thrown-input? :weights (mix [n] {:a 1})) "weights as a map")
    (is (thrown-input? :distrs (mix [] [])) "empty")
    (is (thrown-input? :distrs (mix {:a n} [1])) "distributions as a map")
    (is (thrown-input? :distrs (mix [n 5] [1 1])) "not a distribution")
    (is (thrown-input? :distrs (mix [n nil] [1 1])) "nil")
    (is (thrown-input? :distrs (mix [(r/multi-normal)] [1])) "multivariate")
    (is (thrown-input? :distrs (mix nil [1])))
    (testing "valid"
      (is (== ##-Inf (r/lower-bound (distr/mixture (map #(r/normal {:mu %}) [0 1 2]) (map inc [0 1 2]) nil))) "lazy sequences")
      (is (some? (distr/mixture [n] [3] nil)))
      (is (some? (distr/mixture (list n n) (double-array [1.0 2.0]) nil)) "an array of weights"))))

(deftest mixture-zero-weight-and-continuity
  (let [d (r/mixture {:distrs [(r/uniform-real {:lower 0.0 :upper 1.0}) (r/uniform-real {:lower 5.0 :upper 6.0})]
                      :weights [1 0]})]
    (is (== 0.0 (r/lower-bound d)))
    (is (== 1.0 (r/upper-bound d)) "the component with the zero weight is excluded from bounds")
    (is (m/delta-eq (r/mean d) 0.5 1.0e-12))
    (is (every? #(<= 0.0 % 1.0) (r/->seq d 500))))
  (testing "continuous? is a boolean"
    (is (true? (r/continuous? (r/mixture))))
    (is (false? (r/continuous? (r/mixture {:distrs [(r/poisson) (r/binomial)]}))))
    (is (true? (r/discrete? (r/mixture {:distrs [(r/poisson) (r/binomial)]}))))
    (is (true? (r/continuous? (r/mixture {:distrs [(r/poisson) (r/normal)]})))))
  (testing "id and parameters"
    (let [d (r/mixture)]
      (is (= :mixture (r/distribution-id d)))
      (is (= [:distrs :weights :rng] (r/distribution-parameters d))))))

(deftest mixture-set-seed
  (let [mk (fn [] (r/mixture {:distrs [(r/normal) (r/normal {:mu 5.0})] :weights [0.5 0.5] :rng (r/rng :mersenne 11)}))
        a (r/set-seed! (mk) 7)
        b (r/set-seed! (mk) 7)]
    (is (= (vec (r/->seq a 500)) (vec (r/->seq b 500))) "equal after set-seed!")
    (is (not= (vec (r/->seq (r/set-seed! (mk) 8) 20)) (vec (r/->seq (r/set-seed! (mk) 7) 20)))))
  (testing "the components get different streams"
    (let [d (r/mixture {:distrs [(r/normal) (r/normal {:mu 5.0})] :weights [0.5 0.5]})
          _ (r/set-seed! d 7)
          xs (vec (r/->seq d 4000))
          from-1 (filterv #(< % 2.5) xs)
          from-2 (set (map #(Math/round (* 1.0e9 (- % 5.0))) (filter #(>= % 2.5) xs)))
          matches (count (filter #(contains? from-2 (Math/round (* 1.0e9 (double %)))) from-1))
          rate (/ matches (double (count from-1)))]
      (is (> (count from-1) 1500))
      ;; with the same seed in both components 97-99% of the draws repeat the other component
      (is (<= rate 0.01) (str "matching rate " rate)))))

(deftest mixture-matches-baseline
  (let [b @baseline]
    (testing "seeded construction"
      (is (= (get-in b [:dist-samples :mixture-2])
             (vec (r/->seq (r/mixture {:distrs [(r/normal {:rng (r/rng :mersenne 1)}) (r/normal {:mu 5.0 :rng (r/rng :mersenne 2)})]
                                       :weights [0.3 0.7] :rng (r/rng :mersenne 3)}) 10)))))
    (testing "pdf, cdf, mean, variance"
      (let [v (get-in b [:dist-values :mixture])
            d (r/mixture {:distrs [(r/normal) (r/normal {:mu 5.0})] :weights [0.3 0.7]})]
        (is (= (:pdf v) (mapv #(r/pdf d %) (:xs v))))
        (is (= (:cdf v) (mapv #(r/cdf d %) (:xs v))))
        (is (== (:mean v) (r/mean d)))
        (is (== (:variance v) (r/variance d)))
        ;; the old solver was accurate to 1e-8 absolute; the observed difference is 2.7e-8 at most
        (is (every? true? (map #(m/delta-eq %1 (r/icdf d %2) 1.0e-7) (:icdf v) [0.1 0.25 0.9])))
        (doseq [p [0.1 0.25 0.9]]
          (is (<= (m/abs (- (r/cdf d (r/icdf d p)) p)) 1.0e-10)))))))
