(ns fastmath.random.distribution-overrides-test
  "Tests of the Java subclasses of Commons Math distributions in `fastmath.java.dist` (Random Audit, work package 5.1)
  and of `fastmath.java.SyncedRandomGenerator` (work package 1.8), used directly without `fastmath.random`.

  The five overridden behaviours are checked against scipy (`test/resources/random/dist_override_reference.edn`,
  `utils/fastmath/dev/generate_random_override_reference.py`); everything that is not overridden is checked to be
  bitwise equal to the Commons Math parent class (pdf, log-pdf, cdf, icdf, moments, seeded samples)."
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m])
  (:import [fastmath.java SyncedRandomGenerator]
           [fastmath.java.dist FFixed GumbelFixed LevyFixed LogisticFixed NakagamiFixed]
           [java.lang.reflect Constructor]
           [java.util Arrays]
           [org.apache.commons.math3.distribution AbstractRealDistribution FDistribution GumbelDistribution LevyDistribution
            LogisticDistribution NakagamiDistribution]
           [org.apache.commons.math3.random MersenneTwister RandomGenerator SynchronizedRandomGenerator]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "random/dist_override_reference.edn")))))

(defn- mt ^MersenneTwister [seed] (MersenneTwister. (long seed)))

(def ^:private families
  "Per family: the override class, the Commons parent class, constructors taking a generator and the parameter map of the
  reference file."
  {:logistic {:fixed LogisticFixed :parent LogisticDistribution
              :make-fixed (fn ^AbstractRealDistribution [^RandomGenerator r {:keys [mu s]}] (LogisticFixed. r (double mu) (double s)))
              :make-parent (fn ^AbstractRealDistribution [^RandomGenerator r {:keys [mu s]}] (LogisticDistribution. r (double mu) (double s)))
              :unchanged [#(.getNumericalMean ^AbstractRealDistribution %)]}
   :gumbel {:fixed GumbelFixed :parent GumbelDistribution
            :make-fixed (fn ^AbstractRealDistribution [^RandomGenerator r {:keys [mu beta]}] (GumbelFixed. r (double mu) (double beta)))
            :make-parent (fn ^AbstractRealDistribution [^RandomGenerator r {:keys [mu beta]}] (GumbelDistribution. r (double mu) (double beta)))
            :unchanged [#(.getNumericalVariance ^AbstractRealDistribution %)]}
   :nakagami {:fixed NakagamiFixed :parent NakagamiDistribution
              :make-fixed (fn ^AbstractRealDistribution [^RandomGenerator r {:keys [mu omega]}]
                            (NakagamiFixed. r (double mu) (double omega) NakagamiDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY))
              :make-parent (fn ^AbstractRealDistribution [^RandomGenerator r {:keys [mu omega]}]
                             (NakagamiDistribution. r (double mu) (double omega) NakagamiDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY))
              :unchanged [#(.getNumericalMean ^AbstractRealDistribution %)
                          #(.getNumericalVariance ^AbstractRealDistribution %)]}
   :f {:fixed FFixed :parent FDistribution
       :make-fixed (fn ^AbstractRealDistribution [^RandomGenerator r {:keys [d1 d2]}]
                     (FFixed. r (double d1) (double d2) FDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY))
       :make-parent (fn ^AbstractRealDistribution [^RandomGenerator r {:keys [d1 d2]}]
                      (FDistribution. r (double d1) (double d2) FDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY))
       :unchanged [#(.getNumericalMean ^AbstractRealDistribution %)
                   #(.getNumericalVariance ^AbstractRealDistribution %)]}
   :levy {:fixed LevyFixed :parent LevyDistribution
          :make-fixed (fn ^AbstractRealDistribution [^RandomGenerator r {:keys [mu c]}] (LevyFixed. r (double mu) (double c)))
          :make-parent (fn ^AbstractRealDistribution [^RandomGenerator r {:keys [mu c]}] (LevyDistribution. r (double mu) (double c)))
          :unchanged [#(.getNumericalMean ^AbstractRealDistribution %)
                      #(.getNumericalVariance ^AbstractRealDistribution %)]}})

(defn- bits ^long [^double x] (Double/doubleToLongBits x))

(defn- same-bits? [a b] (= (bits (double a)) (bits (double b))))

(defn- rel-close?
  "`abs(a - b) <= tol * abs(b)`."
  [a b tol]
  (let [a (double a)
        b (double b)]
    (<= (m/abs (- a b)) (* (double tol) (m/abs b)))))

(defn- same-value?
  "`==` for numbers (so `0.0` equals `-0.0`), `##Inf` and `##-Inf` must match, `##NaN` is never the same."
  [a b]
  (== (double a) (double b)))

(defn- fixed [family params] ((get-in families [family :make-fixed]) (mt 1) params))
(defn- parent [family params] ((get-in families [family :make-parent]) (mt 1) params))

(defn- icdf-or-class
  "Inverse cdf, the class of the exception when it throws."
  [^AbstractRealDistribution d ^double p]
  (try (.inverseCumulativeProbability d p)
       (catch Exception e (class e))))

(defn- same-result? [a b]
  (if (number? a)
    (and (number? b) (same-bits? a b))
    (= a b)))

(defn- interior-points
  "About 500 points inside the support: 250 quantiles of the parent in `(0, 1)` and 250 equally spaced points between
  the 0.001 and 0.999 quantiles (the Commons solver returns 0.0 for the tiniest quantiles: points at the support bound
  are dropped, they are the overridden cases)."
  [^AbstractRealDistribution parent]
  (let [lb (.getSupportLowerBound parent)
        qs (->> (range 250)
                (map #(icdf-or-class parent (/ (+ 0.5 (double %)) 250.0)))
                (filter #(and (number? %) (Double/isFinite (double %))))
                (mapv double))
        lo (double (icdf-or-class parent 0.001))
        hi (double (icdf-or-class parent 0.999))
        grid (mapv #(+ lo (* (/ (double %) 249.0) (- hi lo))) (range 250))]
    (filterv #(> (double %) lb) (into qs grid))))

(def ^:private probabilities
  (vec (concat [0.0 1.0e-12 1.0e-6] (map #(/ (+ 0.5 (double %)) 250.0) (range 250)) [0.999999 (- 1.0 1.0e-12) 1.0])))

;; ---------------------------------------------------------------------------------------------
;; structure

(deftest classes-and-constructors-test
  (doseq [[family {:keys [fixed parent]}] families
          :let [^Class fixed fixed
                ^Class parent parent
                ctor-types (fn [^Class c] (set (map (fn [^Constructor k] (vec (.getParameterTypes k))) (.getConstructors c))))]]
    (testing (name family)
      (is (= parent (.getSuperclass fixed)))
      (is (= (ctor-types parent) (ctor-types fixed)) "same public constructors as the parent")))
  (testing "constructors without a generator"
    (is (instance? LogisticDistribution (LogisticFixed. 0.0 1.0)))
    (is (instance? GumbelDistribution (GumbelFixed. 1.0 2.0)))
    (is (instance? NakagamiDistribution (NakagamiFixed. 1.0 1.0)))
    (is (instance? NakagamiDistribution (NakagamiFixed. 1.0 1.0 1.0e-6)))
    (is (instance? FDistribution (FFixed. 1.0 1.0)))
    (is (instance? LevyDistribution (LevyFixed. 0.0 1.0)))))

;; ---------------------------------------------------------------------------------------------
;; overridden behaviours against scipy

(deftest logistic-variance-test
  (doseq [{:keys [params variance mean]} (:logistic @reference)
          :let [^AbstractRealDistribution d (fixed :logistic params)
                ^AbstractRealDistribution p (parent :logistic params)]]
    (testing (str params)
      (is (rel-close? (.getNumericalVariance d) variance 1.0e-12))
      (is (same-value? mean (.getNumericalMean d)))
      (when-not (== 1.0 (double (:s params)))
        (is (not (rel-close? (.getNumericalVariance p) variance 1.0e-3)) "Commons Math 3.6.1 returns pi^2 / (3 s^2)")))))

(deftest gumbel-mean-test
  (doseq [{:keys [params mean variance]} (:gumbel @reference)
          :let [^AbstractRealDistribution d (fixed :gumbel params)
                ^AbstractRealDistribution p (parent :gumbel params)
                beta (double (:beta params))]]
    (testing (str params)
      (is (rel-close? (.getNumericalMean d) mean 1.0e-12))
      (is (rel-close? (.getNumericalVariance d) variance 1.0e-12) "variance is unchanged")
      (is (> (m/abs (- (.getNumericalMean p) (double mean))) (* 6.0e-4 beta)) "Commons Math 3.6.1 uses gamma = 0.5778636..."))))

(deftest nakagami-cdf-test
  (doseq [{:keys [params xs cdf]} (:nakagami @reference)
          :let [^AbstractRealDistribution d (fixed :nakagami params)
                ^AbstractRealDistribution p (parent :nakagami params)]]
    (testing (str params)
      (is (every? true? (map (fn [x e] (same-value? e (.cumulativeProbability d (double x)))) xs cdf)))
      (is (every? zero? (map (fn [x] (.cumulativeProbability d (double x))) xs)))
      (is (every? zero? (map (fn [x] (.density d (double x))) xs)) "density is 0 for x <= 0")
      (is (every? #(= Double/NEGATIVE_INFINITY %) (map (fn [x] (.logDensity d (double x))) xs)) "log-density is -Inf for x <= 0")
      (is (zero? (.cumulativeProbability d Double/NEGATIVE_INFINITY)))
      (is (m/nan? (.cumulativeProbability d Double/NaN)))
      (is (> (.cumulativeProbability p -1.0e6) 0.5) "Commons Math 3.6.1 returns the cdf of |x|"))))

(deftest f-density-test
  (doseq [{:keys [params xs pdf logpdf cdf]} (:f @reference)
          :let [^AbstractRealDistribution d (fixed :f params)
                ^AbstractRealDistribution p (parent :f params)]]
    (testing (str params)
      (is (every? true? (map (fn [x e] (same-value? e (.density d (double x)))) xs pdf)) "density")
      (is (every? true? (map (fn [x e] (same-value? e (.logDensity d (double x)))) xs logpdf)) "log-density")
      (is (every? true? (map (fn [x e] (same-value? e (.cumulativeProbability d (double x)))) xs cdf)) "cdf")
      (is (m/nan? (.density p -1.0)) "Commons Math 3.6.1: NaN")
      (is (m/nan? (.density d Double/NaN)))
      (is (m/nan? (.logDensity d Double/NaN)))
      (is (zero? (.density d Double/NEGATIVE_INFINITY))))))

(deftest f-limit-at-zero-test
  (testing "density at 0 by the numerator degrees of freedom"
    (is (= Double/POSITIVE_INFINITY (.density ^AbstractRealDistribution (FFixed. (mt 1) 1.0 5.0) 0.0)))
    (is (= Double/POSITIVE_INFINITY (.density ^AbstractRealDistribution (FFixed. (mt 1) 1.999999 5.0) 0.0)))
    (is (= 1.0 (.density ^AbstractRealDistribution (FFixed. (mt 1) 2.0 5.0) 0.0)))
    (is (= 0.0 (.density ^AbstractRealDistribution (FFixed. (mt 1) 2.000001 5.0) 0.0)))
    (is (= 0.0 (.density ^AbstractRealDistribution (FFixed. (mt 1) 10.0 5.0) 0.0)))
    (is (= 0.0 (.density ^AbstractRealDistribution (FFixed. (mt 1) 10.0 5.0) -0.0)))
    (is (= 0.0 (.logDensity ^AbstractRealDistribution (FFixed. (mt 1) 2.0 5.0) 0.0)))
    (is (= Double/NEGATIVE_INFINITY (.logDensity ^AbstractRealDistribution (FFixed. (mt 1) 10.0 5.0) 0.0))))
  (testing "constructors without a generator, with a generator and accuracy"
    (is (= 1.0 (.density ^AbstractRealDistribution (FFixed. 2.0 3.0) 0.0)))
    (is (= 1.0 (.density ^AbstractRealDistribution (FFixed. 2.0 3.0 1.0e-6) 0.0)))
    (is (= 1.0 (.density ^AbstractRealDistribution (FFixed. (mt 1) 2.0 3.0) 0.0)))))

(deftest levy-support-test
  (doseq [{:keys [params xs pdf logpdf cdf]} (:levy @reference)
          :let [^AbstractRealDistribution d (fixed :levy params)
                ^AbstractRealDistribution p (parent :levy params)
                mu (double (:mu params))]]
    (testing (str params)
      (is (every? #(<= (double %) mu) xs) "the reference points are outside the open support")
      (is (every? true? (map (fn [x e] (same-value? e (.density d (double x)))) xs pdf)) "density")
      (is (every? true? (map (fn [x e] (same-value? e (.logDensity d (double x)))) xs logpdf)) "log-density")
      (is (every? true? (map (fn [x e] (same-value? e (.cumulativeProbability d (double x)))) xs cdf)) "cdf")
      (is (= 0.0 (.density d mu)) "density at x = mu")
      (is (= Double/NEGATIVE_INFINITY (.logDensity d mu)))
      (is (= 0.0 (.cumulativeProbability d mu)))
      (is (= 0.0 (.cumulativeProbability d Double/NEGATIVE_INFINITY)))
      (is (m/nan? (.density p (- mu 1.0))) "Commons Math 3.6.1: NaN")
      (is (m/nan? (.density d Double/NaN)))
      (is (m/nan? (.cumulativeProbability d Double/NaN)))
      (is (== 1.0 (.cumulativeProbability d Double/POSITIVE_INFINITY))))))

;; ---------------------------------------------------------------------------------------------
;; unchanged behaviours

(deftest interior-bitwise-equal-test
  (doseq [[family {:keys [unchanged]}] families
          {:keys [params]} (get @reference family)
          :let [^AbstractRealDistribution d (fixed family params)
                ^AbstractRealDistribution p (parent family params)
                xs (interior-points p)]]
    (testing (str (name family) " " params)
      (is (>= (count xs) 450))
      (is (every? true? (map (fn [x] (same-bits? (.density d (double x)) (.density p (double x)))) xs)) "density")
      (is (every? true? (map (fn [x] (same-bits? (.logDensity d (double x)) (.logDensity p (double x)))) xs)) "log-density")
      (is (every? true? (map (fn [x] (same-bits? (.cumulativeProbability d (double x)) (.cumulativeProbability p (double x)))) xs))
          "cdf")
      (is (every? true? (map (fn [q] (same-result? (icdf-or-class d (double q)) (icdf-or-class p (double q)))) probabilities))
          "icdf")
      (is (every? true? (map (fn [f] (same-bits? (f d) (f p))) unchanged)) "moments")
      (is (same-bits? (.getSupportLowerBound d) (.getSupportLowerBound p)))
      (is (same-bits? (.getSupportUpperBound d) (.getSupportUpperBound p)))
      (is (= (.isSupportConnected d) (.isSupportConnected p))))))

(deftest seeded-samples-equal-test
  (doseq [[family {:keys [make-fixed make-parent]}] families
          :let [entries (get @reference family)]
          {:keys [params]} [(first entries) (last entries)]]
    (testing (str (name family) " " params)
      (let [^doubles a (.sample ^AbstractRealDistribution (make-fixed (mt 20260507) params) 200000)
            ^doubles b (.sample ^AbstractRealDistribution (make-parent (mt 20260507) params) 200000)]
        (is (= 200000 (alength a)))
        (is (Arrays/equals a b) "identical samples for the same seed")))))

;; ---------------------------------------------------------------------------------------------
;; SyncedRandomGenerator

(defn- ints-of
  "`n` consecutive `nextInt` results of `rng`."
  ^"[I" [^RandomGenerator rng ^long n]
  (let [a (int-array n)]
    (dotimes [i n] (aset a i (.nextInt rng)))
    a))

(defn- ints=? [a b] (Arrays/equals ^ints a ^ints b))

(deftest synced-random-generator-test
  (testing "delegate accessor and class hierarchy"
    (let [base (mt 42)
          s (SyncedRandomGenerator. base)]
      (is (instance? SynchronizedRandomGenerator s))
      (is (instance? RandomGenerator s))
      (is (identical? base (.getDelegate s)))
      (is (instance? MersenneTwister (.getDelegate s)))))
  (testing "same stream as the plain generator, seeding goes to the delegate"
    (let [s (SyncedRandomGenerator. (mt 42))]
      (is (ints=? (ints-of s 100) (ints-of (mt 42) 100)))
      (.setSeed s 7)
      (is (ints=? (ints-of s 100) (ints-of (mt 7) 100)))
      (.setSeed s (int-array [1 2 3]))
      (is (ints=? (ints-of s 100) (ints-of (doto (mt 0) (.setSeed (int-array [1 2 3]))) 100)))
      (.setSeed s 7)
      (is (== (.nextDouble s) (.nextDouble (mt 7))))
      (is (== (.nextLong s) (.nextLong (doto (mt 7) (.nextDouble)))))))
  (testing "the delegate is shared, a second wrapper on it continues the same stream"
    (let [base (mt 5)
          s (SyncedRandomGenerator. base)
          expected (ints-of (mt 5) 20)
          a (ints-of s 10)
          b (ints-of (SyncedRandomGenerator. (.getDelegate s)) 10)]
      (is (ints=? expected (int-array (concat a b))))))
  (testing "8 threads x 10^5 nextInt: the multiset of results equals the first 8*10^5 results of a single thread"
    (let [per-thread 100000
          threads 8
          s (SyncedRandomGenerator. (mt 42))
          futures (mapv (fn [_] (future (ints-of s per-thread))) (range threads))
          results (mapv deref futures)
          ^ints merged (int-array (mapcat vec results))
          ^ints single (ints-of (mt 42) (* per-thread threads))]
      (is (every? #(= per-thread (alength ^ints %)) results))
      (is (= (* per-thread threads) (alength merged)))
      (Arrays/sort merged)
      (Arrays/sort single)
      (is (Arrays/equals merged single)))))
