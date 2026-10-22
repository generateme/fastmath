(ns fastmath.random.utilities-test
  "Tests of the utility constructors of `fastmath.random` (Random Audit, work packages 7.3, 7.8, 7.9, 7.10, 7.11 and 7.12).

  - 7.3  `nbi` / `nbii`: the pmf, cdf and quantiles against R `gamlss.dist` and the Poisson fallback region against the
         exact negative binomial law (mpmath); the validation of the parameters.
  - 7.8, 7.9  `enumerated-real`, `enumerated-int` (sorted pairs) and the validation of the five finite-support families
         (`enumerated-real`, `enumerated-int`, `integer-discrete`, `real-discrete`, `categorical`).
  - 7.10 int typed constructor parameters (2.5, out of the int range, negative, `nil`, a string).
  - 7.11 validation of the parameters of `reciprocal`, `logarithmic`, `negative-binomial`, `power`, `hypoexponential-equal`,
         `hypoexponential` and `beta-symmetrical`.
  - 7.12 `empirical` validation and the documented behaviour of its kernel smoothing, `multinomial` `:n`, the unsupported
         operations of `dirichlet` and `multinomial` and the key phrases of the docstrings.

  The oracle is `test/resources/random/utilities_reference.edn`
  (utils/fastmath/dev/generate_random_utilities_reference.py and .R: R gamlss.dist, mpmath). The enumerated
  distributions are compared with the definition (the smallest `x` with `cdf(x) >= p`) evaluated by hand in exact rational
  arithmetic in the test. The seeded outputs of the code before the audit are `test/resources/random/baseline_reference.edn`."
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m]
            [fastmath.protocols :as prot]
            [fastmath.random :as r]
            [fastmath.special :as spec])
  (:import [umontreal.ssj.probdist AndersonDarlingDist AndersonDarlingDistQuick BetaSymmetricalDist ChiDist CramerVonMisesDist
            HypoExponentialDistEqual KolmogorovSmirnovDist KolmogorovSmirnovDistQuick KolmogorovSmirnovPlusDist
            WatsonGDist WatsonUDist]))

(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "random/utilities_reference.edn")))))

(def ^:private baseline
  (delay (edn/read-string (slurp (io/resource "random/baseline_reference.edn")))))

(defn- rel-err
  "Relative error `abs(a - b) / abs(b)`, the absolute one for `b = 0`."
  ^double [a b]
  (let [a (double a) b (double b)]
    (if (m/zero? b) (m/abs a) (m/abs (/ (- a b) b)))))

(defn- rel-eq? [a b tol] (<= (rel-err a b) tol))

(defn- max-of ^double [xs] (reduce max 0.0 xs))

(defn- bits=
  "Bitwise equality of two doubles (or of two sequences of doubles)."
  [a b]
  (if (and (number? a) (number? b))
    (== (Double/doubleToLongBits (double a)) (Double/doubleToLongBits (double b)))
    (and (= (count a) (count b)) (every? true? (map bits= a b)))))

(defn- nan? [x] (Double/isNaN (double x)))

(defn- thrown-input?
  "True when `f` throws `ex-info` of the type `:invalid-input`."
  [f]
  (try (f) false
       (catch clojure.lang.ExceptionInfo e (= :invalid-input (:type (ex-data e))))
       (catch Throwable _ false)))

(defn- throws?
  "True when `f` throws anything (for the backends which validate on their own)."
  [f]
  (try (f) false (catch Throwable _ true)))

(defn- outcome
  "A value comparable between calls: the cdf at a few points of the distribution built by `f`, `[:ex-info]` for an invalid input
  exception and `[:other class]` for any other exception."
  [f]
  (try (let [d (f)] (mapv #(double (prot/cdf d %)) [0.3 0.7 1.0 2.0 3.0 5.0]))
       (catch clojure.lang.ExceptionInfo _ [:ex-info])
       (catch Throwable t [:other (class t)])))

(defn- ->long [x] (long x))

;; ---------------------------------------------------------------- 7.3 nbi, nbii

(def ^:private nb-constructors {:nbi r/nbi :nbii r/nbii})

(deftest nbi-nbii-against-gamlss
  (doseq [[k ctor] nb-constructors
          {:keys [mu sigma ks pmf lpmf cdf sf ps icdf var]} (get (:gamlss @reference) k)
          :let [d (ctor {:mu mu :sigma sigma})
                ks (map ->long ks)
                label (str k " mu=" mu " sigma=" sigma)]]
    (testing label
      ;; observed maxima over the 12 parameter sets: pmf 3.4e-11 (relative), lpdf 2.3e-10 (abs), cdf and ccdf 2.1e-11 (abs),
      ;; variance 1.1e-13; the tolerances are 30-50 times above, the gamlss values are lgamma differences
      (let [e (max-of (map #(rel-err (prot/pdf d %1) %2) ks pmf))]
        (is (<= e 1.0e-9) (str "pmf rel error " e)))
      (let [e (max-of (map #(m/abs (- (double (prot/lpdf d %1)) %2)) ks lpmf))]
        (is (<= e 1.0e-8) (str "lpdf abs error " e)))
      (let [e (max-of (map #(m/abs (- (double (prot/cdf d %1)) %2)) ks cdf))]
        (is (<= e 1.0e-9) (str "cdf abs error " e)))
      (let [e (max-of (map #(m/abs (- (r/ccdf d %1) %2)) ks sf))]
        (is (<= e 1.0e-9) (str "ccdf abs error " e)))
      (let [bad (remove nil? (map (fn [p q] (when (not= (long q) (long (prot/icdf d p))) [p q (prot/icdf d p)])) ps icdf))]
        (is (empty? bad) (str "icdf mismatches " (vec bad))))
      (is (rel-eq? (prot/mean d) mu 1.0e-9))
      (is (rel-eq? (prot/variance d) var 1.0e-9) (str "variance " (prot/variance d) " expected " var)))))

(deftest nbi-nbii-poisson-fallback-region-against-exact-law
  (doseq [[k ctor] nb-constructors
          {:keys [mu sigma ks pmf cdf lpmf var]} (get (:fallback @reference) k)
          :let [d (ctor {:mu mu :sigma sigma})
                ks (map ->long ks)]]
    (testing (str k " mu=" mu " sigma=" sigma)
      ;; the requirement of the specification: max relative pmf error <= 1e-6 over sigma in [1e-9, 1e-3], mu in {1, 100, 1e4}.
      ;; Observed maxima over the 30 cases: pmf 5.0e-8 (nbii mu=100 sigma=1e-9, the Poisson cells), cdf abs 1.2e-8, lpdf abs
      ;; 5.0e-8, variance relative 1.0e-9; the other tolerances are 10 times above
      (let [e (max-of (map #(rel-err (prot/pdf d %1) %2) ks pmf))]
        (is (<= e 1.0e-6) (str "pmf rel error " e)))
      (let [e (max-of (map #(m/abs (- (double (prot/cdf d %1)) %2)) ks cdf))]
        (is (<= e 1.0e-7) (str "cdf abs error " e)))
      (let [e (max-of (map #(m/abs (- (double (prot/lpdf d %1)) %2)) ks lpmf))]
        (is (<= e 1.0e-6) (str "lpdf abs error " e)))
      (is (rel-eq? (prot/mean d) mu 1.0e-9))
      (is (<= (rel-err (prot/variance d) var) 1.0e-8) (str "variance " (prot/variance d) " expected " var)))))

(deftest fallback-grid-maximum
  (let [errs (for [[k ctor] nb-constructors
                   {:keys [mu sigma ks pmf]} (get (:fallback @reference) k)
                   :let [d (ctor {:mu mu :sigma sigma})]]
               (max-of (map #(rel-err (prot/pdf d (->long %1)) %2) ks pmf)))]
    (is (= 30 (count errs)))
    ;; observed maximum 5.0e-8 (the requirement of the specification is 1e-6)
    (is (<= (max-of errs) 1.0e-6) (str "max pmf rel error over the 30 (family, sigma, mu) cases: " (max-of errs)))))

(deftest nbi-large-mu-small-sigma
  ;; before the audit: the Poisson fallback gave the variance 1e4 (true value 15000) and pdf(1e4) = 0.003989 (true 0.003257)
  (let [d (r/nbi {:mu 1.0e4 :sigma 5.0e-5})
        {:keys [ks pmf]} (first (filter #(and (== 1.0e4 (:mu %)) (== 5.0e-5 (:sigma %))) (get-in @reference [:fallback :nbi])))
        exact (get (zipmap (map long ks) pmf) 10000)]
    (is (rel-eq? (prot/variance d) 15000.0 1.0e-6))
    (is (rel-eq? (prot/mean d) 1.0e4 1.0e-9))
    (is (rel-eq? (prot/pdf d 10000) exact 1.0e-6))
    (is (< (prot/pdf d 10000) 0.0035) "the Poisson value is 0.003989")
    (is (not= :poisson (r/distribution-id d)))
    (is (rel-eq? (prot/variance (r/nbii {:mu 1.0e4 :sigma 5.0e-5})) 10000.5 1.0e-6))))

(deftest nbi-nbii-small-sigma-is-valid
  (testing "tiny positive parameters are valid (the Poisson limit)"
    (doseq [[k ctor] nb-constructors]
      (let [d (ctor {:mu 3.0 :sigma 1.0e-300})]
        (is (rel-eq? (prot/pdf d 3) (/ (* 27.0 (m/exp -3.0)) 6.0) 1.0e-9) (str k)))
      (let [d (ctor {:mu 1.0e-10 :sigma 1.0})]
        (is (rel-eq? (prot/pdf d 0) 1.0 1.0e-9) (str k))))))

(def ^:private invalid-numbers
  [["zero" 0.0] ["negative zero" -0.0] ["negative" -1.0] ["tiny negative" -1.0e-9] ["NaN" ##NaN] ["+Inf" ##Inf] ["-Inf" ##-Inf]])

(deftest nbi-nbii-validation
  (doseq [[k ctor] nb-constructors
          [what v] invalid-numbers]
    (testing (str k " sigma " what)
      (is (thrown-input? #(ctor {:mu 4.0 :sigma v}))))
    (testing (str k " mu " what)
      (is (thrown-input? #(ctor {:mu v :sigma 0.5})))))
  (testing "an explicit nil mu or sigma is an invalid input (documented in the docstring of nbi), the missing key is the default"
    (doseq [[_ ctor] nb-constructors
            opts [{:mu nil} {:sigma nil} {:mu nil :sigma nil}]]
      (is (thrown-input? #(ctor opts)) (str opts)))
    (doseq [[_ ctor] nb-constructors]
      (is (rel-eq? (prot/mean (ctor {:sigma 0.5})) 1.0 1.0e-12))
      (is (rel-eq? (prot/mean (ctor {})) 1.0 1.0e-12)))))

(deftest nbi-nbii-baseline
  (let [{:keys [pdf cdf]} (get-in @baseline [:dist-values :nbi])
        samples (:dist-samples @baseline)
        d (r/nbi {:mu 4.0 :sigma 0.5})]
    (testing "ordinary parameters give the values from before the audit"
      (is (bits= pdf (mapv #(prot/pdf d %) [0 1 5 20])))
      (is (bits= cdf (mapv #(prot/cdf d %) [0 1 5 20])))
      (is (= (:nbi samples) (mapv ->long (r/->seq (r/nbi {:mu 4.0 :sigma 0.5 :rng (r/rng :mersenne 9)}) 20))))
      (is (= (:nbii samples) (mapv ->long (r/->seq (r/nbii {:mu 4.0 :sigma 0.5 :rng (r/rng :mersenne 9)}) 20)))))))

;; ---------------------------------------------------------------- 7.8, 7.9 finite-support families

(defn- hand-support
  "Distinct values of `data` with the cumulative probability, ascending, in exact rational arithmetic (`probs` `nil` = uniform)."
  [data probs]
  (let [w (map rationalize (or probs (repeat (count data) 1)))
        total (reduce + w)
        merged (->> (map vector data w)
                    (group-by first)
                    (map (fn [[x ws]] [x (reduce + (map second ws))]))
                    (sort-by first))]
    (loop [acc [] cum 0N [[x p] & more] merged]
      (if x
        (let [cum (+ cum (/ p total))] (recur (conj acc {:x x :p (/ p total) :cum cum}) cum more))
        acc))))

(defn- hand-icdf
  "The definition: the smallest `x` of the support with `cdf(x) >= p`."
  [support p]
  (some (fn [{:keys [x cum]}] (when (>= cum p) x)) support))

(defn- hand-cdf [support x] (double (reduce + 0 (map :p (filter #(<= (:x %) x) support)))))

(defn- p-grid
  "Quantile levels as rationals: k/100, the ends 0 and 1 (when `ends?`) and a step on each side of every cdf step. The levels
  within 1e-9 of a step of the cdf are dropped (the doubles can not decide there)."
  [support ends?]
  (let [cums (map :cum support)
        eps (/ 1 1000000)
        near? (fn [p] (some #(< (m/abs (double (- p %))) 1.0e-9) cums))
        grid (concat (map #(/ % 100) (range 1 100))
                     (mapcat (fn [c] [(+ c eps) (- c eps)]) cums))
        ends (if (and ends? (every? pos? (map :p support))) [0 1] [])]
    (->> grid (filter #(<= 0 % 1)) (remove near?) (concat ends) distinct sort)))

(def ^:private enumerated-cases
  "[name data probabilities]: unsorted, sorted, reversed, duplicates, uniform, negative and fractional values."
  [["unsorted with a repeated value" [3 1 2 2 5] [0.1 0.2 0.3 0.1 0.3]]
   ["sorted" [1 2 3 5] [0.1 0.4 0.3 0.2]]
   ["reversed" [5 3 1] [0.2 0.3 0.5]]
   ["duplicates" [2 2 2 1 1] [0.1 0.2 0.3 0.15 0.25]]
   ["uniform" [4 2 3 1] nil]
   ["negative values" [-2 10 -7 0 3] [1 2 3 4 5]]
   ["probabilities not summing to 1" [9 7 8] [2 3 5]]
   ["a zero probability" [3 1 2] [0.5 0.0 0.5]]
   ["a single value" [42] nil]])

(deftest enumerated-icdf-is-the-smallest-x-with-cdf-geq-p
  (doseq [[what data probs] enumerated-cases
          [k ctor] {:enumerated-real r/enumerated-real :enumerated-int r/enumerated-int}
          :let [data (if (= k :enumerated-real) (mapv double data) data)
                d (ctor (cond-> {:data data} probs (assoc :probabilities probs)))
                support (hand-support data probs)
                grid (p-grid support true)]]
    (testing (str k " " what)
      (let [bad (remove nil? (map (fn [p] (let [e (hand-icdf support p) g (prot/icdf d (double p))]
                                            (when-not (== e g) [(double p) e g])))
                                  grid))]
        (is (empty? bad) (str (count bad) " of " (count grid) " icdf mismatches, first " (vec (take 3 bad)))))
      (doseq [{:keys [x]} support]
        (is (m/delta-eq (hand-cdf support x) (prot/cdf d x) 1.0e-15) (str "cdf " x))
        (is (m/delta-eq (double (:p (first (filter #(= x (:x %)) support)))) (prot/pdf d x) 1.0e-15) (str "pdf " x))))))

(deftest enumerated-icdf-specification-example
  (let [ps [0.0 0.1 0.2 0.21 0.6 0.61 0.8 0.81 1.0]
        expected [1 1 1 2 2 3 5 5 5]]
    (doseq [[k ctor] {:enumerated-real r/enumerated-real :enumerated-int r/enumerated-int}
            :let [d (ctor {:data [3 1 2 2 5] :probabilities [0.1 0.2 0.3 0.1 0.3]})]]
      ;; before the audit `enumerated-real` returned 3 3 1 1 2 2 5 5 5 (the order of the input)
      (is (= expected (mapv #(long (prot/icdf d %)) ps)) (str k)))
    (let [d (r/enumerated-real {:data [5.0 3.0 1.0] :probabilities [0.2 0.3 0.5]})]
      (is (= [1 1 1 3 3 5 5] (mapv #(long (prot/icdf d %)) [0.0 0.2 0.5 0.51 0.8 0.81 1.0]))
          "reversed data"))))

(deftest discrete-families-icdf-is-the-smallest-x-with-cdf-geq-p
  (doseq [[what data probs] (remove (comp nil? #(nth % 2)) enumerated-cases)
          [k ctor] {:integer-discrete r/integer-discrete :real-discrete r/real-discrete}
          :let [data (if (= k :real-discrete) (mapv double data) data)
                d (ctor {:data data :probabilities probs})
                support (hand-support data probs)]]
    (testing (str k " " what)
      (let [bad (remove nil? (map (fn [p] (let [e (hand-icdf support p) g (prot/icdf d (double p))]
                                            (when-not (== e g) [(double p) e g])))
                                  (p-grid support false)))]
        (is (empty? bad) (str (count bad) " icdf mismatches, first " (vec (take 3 bad)))))
      (doseq [{:keys [x cum p]} support]
        (is (m/delta-eq (double p) (prot/pdf d x) 1.0e-15) (str "pdf " x))
        (is (m/delta-eq (double cum) (prot/cdf d x) 1.0e-15) (str "cdf " x))))))

(deftest enumerated-real-streams
  (let [samples (:dist-samples @baseline)
        s (fn [data probs] (vec (r/->seq (r/enumerated-real {:data data :probabilities probs :rng (r/rng :mersenne 9)}) 20)))]
    (testing "sorted data: the stream of the code before the audit"
      (is (= (:enumerated-real-sorted samples) (s [1.0 2.0 3.0 5.0] [0.1 0.4 0.3 0.2]))))
    (testing "unsorted data: the pairs are sorted, the stream is the one of the sorted pairs"
      (is (= (s [1.0 2.0 3.0 5.0] [0.4 0.3 0.1 0.2]) (s [3.0 1.0 2.0 5.0] [0.1 0.4 0.3 0.2])))
      (is (= (s [1.0 2.0 3.0 5.0] [0.4 0.3 0.1 0.2]) (s [5.0 2.0 1.0 3.0] [0.2 0.3 0.4 0.1]))))))

(deftest enumerated-valid-inputs
  (testing "enumerated-int at the edges of the int range"
    (let [d (r/enumerated-int {:data [2147483647 -2147483648] :probabilities [0.25 0.75]})]
      (is (== 0.25 (prot/pdf d 2147483647)))
      (is (== 0.75 (prot/cdf d -2147483648)))))
  (testing "integer-discrete keeps the longs outside of the int range"
    (let [d (r/integer-discrete {:data [1 2 1099511627776] :probabilities [0.5 0.3 0.2]})]
      (is (m/delta-eq 0.2 (prot/pdf d 1099511627776) 1.0e-15))
      (is (rel-eq? (prot/mean d) (+ 0.5 0.6 (* 0.2 1099511627776.0)) 1.0e-12))))
  (testing "mean and variance of the numeric families"
    (doseq [ctor [r/enumerated-real r/enumerated-int r/integer-discrete r/real-discrete]
            :let [d (ctor {:data [3 1 2 2 5] :probabilities [0.1 0.2 0.3 0.1 0.3]})
                  mean (+ 0.2 (* 2 0.4) 0.3 1.5)
                  var (- (+ 0.2 (* 4 0.4) (* 9 0.1) (* 25 0.3)) (* mean mean))]]
      (is (rel-eq? (prot/mean d) mean 1.0e-12))
      (is (rel-eq? (prot/variance d) var 1.0e-12)))
    (testing "huge range of the weights"
      (let [d (r/integer-discrete {:data [1 2] :probabilities [1.0e-300 1.0]})]
        (is (== 1.0e-300 (prot/pdf d 1)))))))

(def ^:private finite-support-constructors
  {:enumerated-real r/enumerated-real
   :enumerated-int r/enumerated-int
   :integer-discrete r/integer-discrete
   :integer-discrete-distribution r/integer-discrete-distribution
   :real-discrete r/real-discrete
   :real-discrete-distribution r/real-discrete-distribution
   :categorical r/categorical
   :categorical-distribution r/categorical-distribution})

(def ^:private invalid-finite-support
  [["empty data" {:data []}]
   ["empty data and probabilities" {:data [] :probabilities []}]
   ["fewer probabilities" {:data [1 2 3] :probabilities [0.5 0.5]}]
   ["more probabilities" {:data [1 2 3] :probabilities [0.2 0.3 0.5 0.5]}]
   ["one probability for three values" {:data [1 2 3] :probabilities [1.0]}]
   ["negative probability" {:data [1 2 3] :probabilities [0.5 -0.1 0.6]}]
   ["NaN probability" {:data [1 2 3] :probabilities [0.5 ##NaN 0.5]}]
   ["+Inf probability" {:data [1 2 3] :probabilities [0.5 ##Inf 0.5]}]
   ["-Inf probability" {:data [1 2 3] :probabilities [0.5 ##-Inf 0.5]}]
   ["all zero probabilities" {:data [1 2 3] :probabilities [0 0 0]}]
   ["all zero double probabilities" {:data [1 2 3] :probabilities [0.0 0.0 0.0]}]
   ["negative zero sum" {:data [1 2] :probabilities [0.0 -0.0]}]])

(deftest finite-support-validation-table
  (doseq [[k ctor] finite-support-constructors
          [what opts] invalid-finite-support]
    (testing (str k " " what)
      (is (thrown-input? #(ctor opts))))))

(deftest finite-support-type-specific-validation
  (testing "integer data"
    (doseq [[k ctor] {:enumerated-int r/enumerated-int :integer-discrete r/integer-discrete}
            [what data] [["fractional" [1.5 2.5]] ["one fractional" [1 2.5 3]] ["a string" [1 "a"]]]]
      (is (thrown-input? #(ctor {:data data :probabilities (vec (repeat (count data) 1.0))})) (str k " " what))))
  (testing "enumerated-int data must fit into an int"
    (doseq [data [[1 2 2147483648] [-2147483649 1] [1 2 1099511627776] [1 2 3000000000]]]
      (is (thrown-input? #(r/enumerated-int {:data data})) (str data))
      (is (thrown-input? #(r/enumerated-int {:data data :probabilities (repeat (count data) 1.0)})) (str data))))
  (testing "integer-discrete data must fit into a long"
    (is (thrown-input? #(r/integer-discrete {:data [1 (biginteger "9223372036854775808")]})))
    (is (thrown-input? #(r/integer-discrete {:data [-1 (biginteger "-9223372036854775809")]})))
    (is (some? (r/integer-discrete {:data [4611686018427387904 -4611686018427387904]})))))

(deftest integer-discrete-long-range-edges
  ;; the specification (WP 7.9): the data are integers within the long range, so the ends of the range are valid
  (doseq [x [Long/MAX_VALUE Long/MIN_VALUE (dec Long/MAX_VALUE) (inc Long/MIN_VALUE)]]
    (let [d (try (r/integer-discrete {:data [x 0]}) (catch Throwable t t))]
      (is (not (instance? Throwable d)) (str x " -> " (when (instance? Throwable d) (ex-message d))))
      ;; the pdf/cdf arguments are doubles, and every long within 1024 of the ends rounds to +-2^63, which is outside the
      ;; exactly representable long range, so the mass at the ends is observable through `sample` only
      (when-not (instance? Throwable d)
        ;; the data are held as doubles: beyond 2^53 a value is returned as the nearest double (saturating at the ends)
        (is (every? #{(long (double x)) 0} (take 50 (prot/->seq d))) (str x))))))

(deftest finite-support-valid-inputs-baseline
  (testing "the values of the code before the audit"
    (let [{:keys [pdf cdf]} (get-in @baseline [:dist-values :integer-discrete])
          d (r/integer-discrete {:data [1 2 4] :probabilities [0.2 0.3 0.5]})]
      (is (bits= pdf (mapv #(prot/pdf d %) [0 1 2 3 4 5])))
      (is (bits= cdf (mapv #(prot/cdf d %) [0 1 2 3 4 5])))))
  (testing "no probabilities = uniform, a categorical with any values"
    (doseq [[k ctor] finite-support-constructors
            :let [d (ctor {:data [1 2 3 4]})]]
      (is (m/delta-eq 0.25 (prot/pdf d 3) 1.0e-15) (str k))
      (is (m/delta-eq 0.5 (prot/cdf d 2) 1.0e-15) (str k)))
    (let [d (r/categorical {:data [:a :b :a] :probabilities [1 1 2]})]
      (is (m/delta-eq 0.75 (prot/pdf d :a) 1.0e-15))
      (is (m/delta-eq 0.25 (prot/pdf d :b) 1.0e-15)))))

;; ---------------------------------------------------------------- 7.10 int typed parameters

(def ^:private cdf-points [0.3 0.7 1.0 2.0 3.0 5.0])

(defn- choose [n k] (reduce * 1 (map #(/ (- n %) (inc %)) (range k))))

(def ^:private int-constructors
  "One entry per constructor: the keys typed `int`, a valid parameter set, `:negative` the keys which are invalid at -1
  and `:expected` the function building the same law from the closed form or from the backend class directly."
  [{:name :binomial :ctor r/binomial :base {:trials 10 :p 0.3} :keys [:trials] :negative [:trials]
    :expected (fn [x] (double (reduce + (for [k (range (inc (long (m/floor x))))] (* (choose 10 k) (m/pow 0.3 k) (m/pow 0.7 (- 10 k)))))))}
   {:name :pascal :ctor r/pascal :base {:r 4 :p 0.3} :keys [:r] :negative [:r]
    :expected (fn [x] (double (reduce + (for [k (range (inc (long (m/floor x))))] (* (choose (+ k 3) k) (m/pow 0.3 4) (m/pow 0.7 k))))))}
   {:name :zipf :ctor r/zipf :base {:number-of-elements 20 :exponent 1.5} :keys [:number-of-elements] :negative [:number-of-elements]
    :expected (fn [x] (let [h (reduce + (map #(m/pow % -1.5) (range 1 21)))]
                        (double (/ (reduce + 0.0 (map #(m/pow % -1.5) (range 1 (inc (min 20 (long (m/floor x))))))) h))))}
   {:name :hypergeometric :ctor r/hypergeometric :base {:population-size 30 :number-of-successes 12 :sample-size 8}
    :keys [:population-size :number-of-successes :sample-size] :negative [:population-size :number-of-successes :sample-size]
    :expected (fn [x] (double (reduce + (for [k (range (inc (long (m/floor x))))]
                                          (/ (* (choose 12 k) (choose 18 (- 8 k))) (choose 30 8))))))}
   {:name :uniform-int :ctor r/uniform-int :base {:lower 2 :upper 9} :keys [:lower :upper] :negative []
    :expected (fn [x] (double (/ (m/constrain (- (m/floor x) 1.0) 0.0 8.0) 8.0)))}
   {:name :chi :ctor r/chi :base {:nu 3} :keys [:nu] :negative [:nu] :expected #(.cdf (ChiDist. 3) (double %))}
   {:name :erlang :ctor r/erlang :base {:k 3 :lambda 2.0} :keys [:k] :negative [:k]
    :expected (fn [x] (let [lx (* 2.0 (double x))] (- 1.0 (* (m/exp (- lx)) (+ 1.0 lx (/ (* lx lx) 2.0))))))}
   {:name :kolmogorov-smirnov :ctor r/kolmogorov-smirnov :base {:n 10} :keys [:n] :negative [:n]
    :expected #(.cdf (KolmogorovSmirnovDist. 10) (double %))}
   {:name :kolmogorov-smirnov+ :ctor r/kolmogorov-smirnov+ :base {:n 10} :keys [:n] :negative [:n]
    :expected #(.cdf (KolmogorovSmirnovPlusDist. 10) (double %))}
   {:name :kolmogorov-smirnov-quick :ctor r/kolmogorov-smirnov-quick :base {:n 10} :keys [:n] :negative [:n]
    :expected #(.cdf (KolmogorovSmirnovDistQuick. 10) (double %))}
   {:name :anderson-darling :ctor r/anderson-darling :base {:n 10} :keys [:n] :negative [:n]
    :expected #(.cdf (AndersonDarlingDist. 10) (double %))}
   {:name :anderson-darling-quick :ctor r/anderson-darling-quick :base {:n 10} :keys [:n] :negative [:n]
    :expected #(.cdf (AndersonDarlingDistQuick. 10) (double %))}
   {:name :cramer-von-mises :ctor r/cramer-von-mises :base {:n 10} :keys [:n] :negative [:n]
    :expected #(.cdf (CramerVonMisesDist. 10) (double %))}
   {:name :watson-g :ctor r/watson-g :base {:n 10} :keys [:n] :negative [:n]
    :expected #(.cdf (WatsonGDist. 10) (double %))}
   {:name :watson-u :ctor r/watson-u :base {:n 10} :keys [:n] :negative [:n]
    :expected #(.cdf (WatsonUDist. 10) (double %))}
   {:name :beta-symmetrical :ctor r/beta-symmetrical :base {:alpha 2.5 :d 10} :keys [:d] :negative [:d]
    :expected #(.cdf (BetaSymmetricalDist. 2.5 10) (double %))}
   {:name :hypoexponential-equal :ctor r/hypoexponential-equal :base {:n 5 :k 3 :h 0.5} :keys [:n :k] :negative [:n :k]
    :expected #(.cdf (HypoExponentialDistEqual. 5 3 0.5) (double %))}])

(deftest int-parameters-valid-values-unchanged
  (doseq [{:keys [name ctor base expected]} int-constructors
          :let [d (ctor base)]]
    (testing name
      (doseq [x (conj cdf-points 0.1 0.9 4.0)]
        (is (m/delta-eq (double (expected x)) (double (prot/cdf d x)) 1.0e-13) (str "cdf " x))))))

(deftest int-parameters-ssj-pass-through-is-bitwise
  (doseq [{:keys [name ctor base expected]} int-constructors
          :when (#{:chi :kolmogorov-smirnov :kolmogorov-smirnov+ :kolmogorov-smirnov-quick :anderson-darling :anderson-darling-quick
                   :cramer-von-mises :watson-g :watson-u :beta-symmetrical :hypoexponential-equal} name)
          :let [d (ctor base)]]
    (testing name
      (is (bits= (mapv #(double (expected %)) cdf-points) (mapv #(prot/cdf d %) cdf-points))))))

(deftest int-parameters-invalid-values
  (doseq [{:keys [name ctor base keys negative]} int-constructors
          k keys]
    (testing (str name " " k)
      (testing "fractional"
        (is (thrown-input? #(ctor (assoc base k 2.5))) "2.5")
        (is (thrown-input? #(ctor (assoc base k 0.5))) "0.5"))
      (testing "out of the int range"
        (is (thrown-input? #(ctor (assoc base k 2147483648))) "2^31")
        (is (thrown-input? #(ctor (assoc base k 3000000000))) "3000000000")
        (is (thrown-input? #(ctor (assoc base k -2147483649))) "-2^31-1")
        (is (thrown-input? #(ctor (assoc base k Long/MAX_VALUE))) "Long/MAX_VALUE"))
      (testing "not a number"
        (is (thrown-input? #(ctor (assoc base k "a"))) "\"a\"")
        (is (thrown-input? #(ctor (assoc base k :a))) ":a"))
      (testing "negative, where the parameter is invalid (the backend may throw its own exception)"
        (when (some #{k} negative)
          (is (throws? #(ctor (assoc base k -1))) "-1")))
      (testing "nil is the missing option"
        (is (= (outcome #(ctor (dissoc base k))) (outcome #(ctor (assoc base k nil)))))))))

(deftest int-parameters-edges
  (testing "the int range itself is valid"
    (is (some? (r/binomial {:trials Integer/MAX_VALUE})))
    (is (== 0.0 (prot/cdf (r/binomial {:trials 0 :p 0.3}) -1)))
    (is (== 1.0 (prot/pdf (r/binomial {:trials 0 :p 0.3}) 0)))
    (let [d (r/uniform-int {:lower Integer/MIN_VALUE :upper -2147483641})]
      (is (m/delta-eq 0.125 (prot/pdf d -2147483645) 1.0e-15)))
    (is (thrown-input? #(r/uniform-int {:lower 0 :upper -2147483649})))
    (is (thrown-input? #(r/uniform-int {:lower 2147483648 :upper 2147483649}))))
  (testing "uniform-int default: nil lower and upper"
    (is (= (outcome #(r/uniform-int)) (outcome #(r/uniform-int {:lower nil :upper nil}))))))

;; ---------------------------------------------------------------- 7.11 constructor validation

(deftest reciprocal-validation
  (doseq [[what opts] [["a = 0" {:a 0.0 :b 1.0}] ["a negative" {:a -1.0 :b 1.0}] ["a = -0.0" {:a -0.0 :b 1.0}]
                       ["a = NaN" {:a ##NaN :b 1.0}] ["b = NaN" {:a 1.0 :b ##NaN}] ["equal bounds" {:a 2.0 :b 2.0}]
                       ["reversed bounds" {:a 5.0 :b 2.0}] ["both negative" {:a -5.0 :b -2.0}]]]
    (testing what (is (thrown-input? #(r/reciprocal opts)))))
  (testing "the first valid values"
    (let [d (r/reciprocal {:a 1.0e-300 :b 1.0})]
      (is (rel-eq? (prot/cdf d 1.0e-150) 0.5 1.0e-9)))
    (let [b (m/next-double 1.0)
          d (r/reciprocal {:a 1.0 :b b})]
      (is (== 1.0 (prot/lower-bound d)))
      (is (== b (prot/upper-bound d)))
      (is (== 1.0 (prot/cdf d b))))
    (let [d (r/reciprocal {:a 1.0 :b 10.0})]
      (is (rel-eq? (prot/cdf d 3.0) (/ (m/log 3.0) (m/log 10.0)) 1.0e-14)))))

(deftest logarithmic-validation
  (doseq [[what p] [["p = 0" 0.0] ["p = -0.0" -0.0] ["p = 1" 1.0] ["negative" -0.1] ["above 1" 1.1] ["NaN" ##NaN] ["+Inf" ##Inf]
                    ["-Inf" ##-Inf]]]
    (testing what (is (thrown-input? #(r/logarithmic {:p p})))))
  (testing "the first valid values: pmf(k) = -p^k / (k ln(1-p))"
    ;; observed: relative error 8.3e-8 at p = 1e-10 (cancellation in log(1 - p)), below 1e-12 for the other p
    (doseq [p [1.0e-10 (m/prev-double 1.0) 0.5]
            :let [d (r/logarithmic {:p p})
                  c (/ -1.0 (m/log1p (- p)))]]
      (is (rel-eq? (prot/pdf d 1) (* c p) 1.0e-6) (str p))
      (is (rel-eq? (prot/pdf d 3) (/ (* c p p p) 3.0) 1.0e-6) (str p)))))

(deftest negative-binomial-validation
  (doseq [[what opts] [["r = 0" {:r 0.0 :p 0.5}] ["r negative" {:r -1.0 :p 0.5}] ["r NaN" {:r ##NaN :p 0.5}] ["p = 0" {:r 5.0 :p 0.0}] ["p negative" {:r 5.0 :p -0.1}]
                       ["p above 1" {:r 5.0 :p 1.5}] ["p just above 1" {:r 5.0 :p (m/next-double 1.0)}] ["p NaN" {:r 5.0 :p ##NaN}]]]
    (testing what (is (thrown-input? #(r/negative-binomial opts)))))
  (testing "the first valid values"
    (let [d (r/negative-binomial {:r 5.0 :p 1.0})]
      (is (== 1.0 (prot/pdf d 0)))
      (is (== 0.0 (prot/pdf d 1)))
      (is (== 0.0 (prot/mean d))))
    (let [d (r/negative-binomial {:r 1.0e-10 :p 0.5})]
      (is (rel-eq? (prot/pdf d 0) (m/pow 0.5 1.0e-10) 1.0e-12)))
    (let [d (r/negative-binomial {:r 5.0 :p (m/prev-double 1.0)})]
      (is (> (prot/pdf d 0) 0.999999)))
    (testing "r is a real number: pmf(k) = Gamma(k+r) / (k! Gamma(r)) p^r (1-p)^k"
      (let [d (r/negative-binomial {:r 2.5 :p 0.4})]
        (is (rel-eq? (prot/mean d) (/ (* 2.5 0.6) 0.4) 1.0e-12))
        (doseq [k [0 1 4]]
          (is (rel-eq? (prot/pdf d k)
                       (m/exp (+ (- (spec/log-gamma (+ k 2.5)) (spec/log-gamma 2.5) (spec/log-gamma (inc k)))
                                 (* 2.5 (m/log 0.4)) (* k (m/log 0.6))))
                       1.0e-10) (str k)))))))

(deftest power-validation
  (doseq [[what opts] [["a = b" {:a 1.0 :b 1.0 :c 2.0}] ["a > b" {:a 1.0 :b 0.0 :c 2.0}] ["c = 0" {:a 0.0 :b 1.0 :c 0.0}]
                       ["c negative" {:a 0.0 :b 1.0 :c -1.0}] ["c NaN" {:a 0.0 :b 1.0 :c ##NaN}] ["a NaN" {:a ##NaN :b 1.0 :c 2.0}]
                       ["b NaN" {:a 0.0 :b ##NaN :c 2.0}] ["c = -0.0" {:a 0.0 :b 1.0 :c -0.0}]]]
    (testing what (is (thrown-input? #(r/power opts)))))
  (testing "the first valid values: cdf((x-a)/(b-a))^c"
    (let [d (r/power {:a 1.0 :b (m/next-double 1.0) :c 2.0})]
      (is (== 1.0 (prot/cdf d (m/next-double 1.0)))))
    (let [d (r/power {:a -1.0 :b 3.0 :c 1.0e-3})]
      (is (rel-eq? (prot/cdf d 0.0) (m/pow 0.25 1.0e-3) 1.0e-12)))
    (let [d (r/power {:a 1.0 :b 3.0 :c 3.0})]
      (is (rel-eq? (prot/cdf d 2.0) 0.125 1.0e-12))
      (is (rel-eq? (prot/mean d) (+ 1.0 (* 2.0 0.75)) 1.0e-12)))))

(deftest hypoexponential-equal-validation
  (doseq [[what opts] [["k = 0" {:n 3 :k 0 :h 1.0}] ["k negative" {:n 3 :k -1 :h 1.0}] ["n < k" {:n 2 :k 3 :h 1.0}]
                       ["n = 0" {:n 0 :k 0 :h 1.0}] ["n negative" {:n -2 :k 1 :h 1.0}] ["h = 0" {:n 3 :k 2 :h 0.0}]
                       ["h negative" {:n 3 :k 2 :h -1.0}] ["h NaN" {:n 3 :k 2 :h ##NaN}]
                       ["h = -0.0" {:n 3 :k 2 :h -0.0}]]]
    (testing what (is (thrown-input? #(r/hypoexponential-equal opts)))))
  (testing "the first valid values: n = k = 1 and n = k"
    (is (rel-eq? (prot/mean (r/hypoexponential-equal {:n 1 :k 1 :h 2.0})) 0.5 1.0e-12))
    (is (rel-eq? (prot/mean (r/hypoexponential-equal {:n 4 :k 4 :h 1.0})) (+ 1.0 0.5 (/ 1.0 3.0) 0.25) 1.0e-12)))
  (testing "the rates are (n+1-i) h, i = 1..k"
    (doseq [[n k h] [[5 3 0.5] [4 1 2.0] [6 6 0.25] [10 2 1.5]]
            :let [rates (map #(* (- (inc n) %) h) (range 1 (inc k)))
                  d (r/hypoexponential-equal {:n n :k k :h h})]]
      (is (rel-eq? (prot/mean d) (reduce + (map #(/ 1.0 %) rates)) 1.0e-12) (str [n k h]))
      (is (rel-eq? (prot/variance d) (reduce + (map #(/ 1.0 (* % %)) rates)) 1.0e-12) (str [n k h])))))

(deftest hypoexponential-validation
  (doseq [[what opts] [["empty" {:lambdas []}] ["zero rate" {:lambdas [1.0 0.0]}] ["negative rate" {:lambdas [-1.0 2.0]}]
                       ["NaN rate" {:lambdas [1.0 ##NaN]}] ["a number" {:lambdas 2.0}]]]
    (testing what (is (thrown-input? #(r/hypoexponential opts)))))
  (testing "valid"
    (is (rel-eq? (prot/mean (r/hypoexponential {:lambdas [1.0]})) 1.0 1.0e-12))
    (is (rel-eq? (prot/mean (r/hypoexponential {:lambdas [1.0 2.0 4.0]})) 1.75 1.0e-12))
    (is (rel-eq? (prot/variance (r/hypoexponential {:lambdas [1.0 2.0 4.0]})) (+ 1.0 0.25 0.0625) 1.0e-12))))

(deftest beta-symmetrical-validation
  (doseq [[what opts] [["alpha = 0" {:alpha 0.0 :d 10}] ["alpha negative" {:alpha -1.0 :d 10}] ["alpha NaN" {:alpha ##NaN :d 10}] ["d negative" {:alpha 2.0 :d -1}] ["d just above the range (36)" {:alpha 2.0 :d 36}]
                       ["d far above" {:alpha 2.0 :d 100}]]]
    (testing what (is (thrown-input? #(r/beta-symmetrical opts)))))
  (testing "d is valid from 0 to 35 (the size of the table of SSJ)"
    (doseq [d [0 35]]
      (is (m/delta-eq 0.5 (prot/cdf (r/beta-symmetrical {:alpha 2.5 :d d}) 0.5) 1.0e-9) (str "d=" d))))
  (testing "valid: the cdf of the symmetric beta at 1/2 is 1/2"
    (doseq [d [1 5 10 14]]
      (is (m/delta-eq 0.5 (prot/cdf (r/beta-symmetrical {:alpha 2.5 :d d}) 0.5) 1.0e-9) (str "d=" d)))))

;; ---------------------------------------------------------------- 7.12 empirical, multinomial, dirichlet

(deftest empirical-validation
  (doseq [[what opts] [["empty data" {:data []}] ["empty lazy data" {:data (range 0)}]
                       ["bin-count 0" {:data [1 2 3 4 5] :bin-count 0}] ["bin-count negative" {:data [1 2 3 4 5] :bin-count -3}]
                       ["bin-count fractional" {:data [1 2 3 4 5] :bin-count 4.7}] ["bin-count 0.5" {:data [1 2 3 4 5] :bin-count 0.5}]
                       ["bin-count a string" {:data [1 2 3 4 5] :bin-count "a"}] ["bin-count NaN" {:data [1 2 3 4 5] :bin-count ##NaN}]]]
    (testing what (is (thrown-input? #(r/empirical opts)))))
  (testing "nil is the missing option, the data of the default distribution or an invalid input"
    (let [o (outcome #(r/empirical {:data nil}))]
      (is (or (= [:ex-info] o) (= (outcome #(r/empirical)) o)) (str "nil data: " o)))
    (is (= (outcome #(r/empirical {:data [1 2 3 4 5]})) (outcome #(r/empirical {:data [1 2 3 4 5] :bin-count nil})))))
  (testing "valid"
    (let [d (r/empirical {:data [1 2 3 4 5] :bin-count 1})]
      (is (rel-eq? (prot/mean d) 3.0 1.0e-12)))
    (let [d (r/empirical {:data [1 2 3 4 5 6 7 8] :bin-count 4})]
      (is (== 1.0 (prot/lower-bound d)))
      (is (== 8.0 (prot/upper-bound d))))))

(deftest empirical-documented-behaviour
  ;; Characterisation tests: they pin the behaviour which the docstring of `empirical` documents (Commons
  ;; EmpiricalDistribution smooths with a Gaussian kernel per bin, it is not a histogram). The sampling statistics are
  ;; deliberately not tested (kernel smoothing is not uniform within a bin).
  (testing "pdf is NaN in a bin with fewer than 2 observations or no spread: [1 1 1 1 9], 4 bins"
    (let [d (r/empirical {:data [1 1 1 1 9] :bin-count 4})]
      (is (nan? (prot/pdf d 2.0)))
      (is (nan? (prot/pdf d 3.0)))
      (is (nan? (prot/pdf d 5.0)))
      (is (nan? (prot/pdf d 6.0)))
      (is (== 0.0 (prot/pdf d 8.0)))))
  (testing "cdf at the lower bound is not 0"
    (let [d (r/empirical {:data [1 1 1 1 9] :bin-count 4})]
      (is (not (m/zero? (double (prot/cdf d (prot/lower-bound d)))))))
    (let [g (r/rng :mersenne 3)
          data (vec (repeatedly 5000 #(r/grandom g)))
          d (r/empirical {:data data :rng (r/rng :mersenne 4)})]
      (is (pos? (prot/cdf d (prot/lower-bound d))))
      (is (< (prot/cdf d (prot/lower-bound d)) 1.0e-3)))))

(deftest multinomial-n
  (testing "invalid n"
    (doseq [n [2.5 -1 "a" ##NaN 2147483648 3000000000]]
      (is (thrown-input? #(r/multinomial {:n n})) (str n))))
  (testing "n = 0 is the distribution of the empty trial"
    (let [d (r/multinomial {:n 0 :ps [1 1]})]
      (is (== 1.0 (prot/pdf d [0 0])))
      (is (== 0.0 (prot/pdf d [1 0])))))
  (testing "nil n is the default 10"
    (is (rel-eq? (prot/pdf (r/multinomial {:n nil}) [5 5]) 0.24609375 1.0e-12))
    (is (rel-eq? (prot/pdf (r/multinomial {:n nil :ps [1 1 2]}) [3 3 4])
                 (/ (* 4200.0 (m/pow 0.25 6) (m/pow 0.5 4)) 1.0) 1.0e-12))
    (is (rel-eq? (prot/pdf (r/multinomial) [5 5]) 0.24609375 1.0e-12))))

(deftest unsupported-operations-are-documented
  ;; The multivariate distributions expose no icdf (dirichlet no cdf either): the call fails with the class which the
  ;; docstring names. The class itself is pinned (a change of the backend would change both).
  (doseq [[v f what] [[#'r/dirichlet #(prot/cdf (r/dirichlet) [0.5 0.5]) "dirichlet cdf"]
                      [#'r/dirichlet #(prot/icdf (r/dirichlet) 0.5) "dirichlet icdf"]
                      [#'r/multinomial #(prot/icdf (r/multinomial) 0.5) "multinomial icdf"]]]
    (let [t (try (f) nil (catch Throwable t t))
          doc (str/lower-case (:doc (meta v)))]
      (testing what
        (is (some? t) "throws")
        (is (= "AbstractMethodError" (some-> t class .getSimpleName)) "the class of the exception")
        (is (str/includes? doc (str/lower-case (some-> t class .getSimpleName))) "the docstring names the class")))))

(def ^:private docstring-phrases
  "var -> case-insensitive substrings of its docstring (a short keyword per documented fact, agreed after reading the
  finished docstrings)."
  {#'r/multinomial ["AbstractMethodError" "IllegalArgumentException" "n"]
   #'r/dirichlet ["AbstractMethodError" "simplex"]
   #'r/uniform-int ["Integer/MAX_VALUE" "overflows"]
   #'r/chi-squared-noncentral ["lambda = 0.0" "not accepted"]
   #'r/log-logistic ["shape parameter" "scale parameter"]
   #'r/nbi ["quadratic" "poisson"]
   #'r/folded-normal ["non-negative"]
   #'r/hypoexponential-equal ["(n + 1 - i) * h"]
   #'r/empirical ["gaussian kernel" "not uniformly" "fewer than two observations"]
   #'r/kde ["alias"]
   #'r/continuous-distribution ["not the moments" "steps"]
   #'r/half-cauchy ["lower-bound"]
   #'r/logarithmic ["smallest"]
   #'r/enumerated-real ["sorted ascending" "integer-discrete"]
   #'r/enumerated-int ["int range"]
   #'r/integer-discrete-distribution ["long range"]
   #'r/real-discrete-distribution ["sorted"]
   #'r/categorical-distribution ["first appearance"]})

(deftest docstrings-contain-the-key-phrases
  (doseq [[v phrases] docstring-phrases
          :let [doc (str/lower-case (or (:doc (meta v)) ""))]
          phrase phrases]
    (is (str/includes? doc (str/lower-case phrase)) (str (:name (meta v)) " docstring contains \"" phrase "\""))))
