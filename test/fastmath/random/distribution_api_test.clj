(ns fastmath.random.distribution-api-test
  "Tests of the protocol layer of `fastmath.random` (Random Audit, work packages 5.1-5.8): the user facing proxies `cdf`,
  `ccdf`, `pdf`, `lpdf`, `probability`, `icdf`, the wiring of the Commons overrides (`fastmath.java.dist`), of
  `inverse-gaussian` and `multi-normal` through the public constructors, `distribution-id` and `distribution-parameters`,
  and the documentation of the proxies.

  Oracles (all outside the implementation):
  - `test/resources/random/dist_override_reference.edn`: scipy, the five Commons overrides (5.1);
  - `test/resources/random/native_distributions_reference.edn`: mpmath, `inverse-gaussian` (5.2), scipy `multivariate_normal` (5.5);
  - `test/resources/random/truncated_mixture_reference.edn`: scipy `sf`, `ccdf` of the closed-form and SSJ families (5.4);
  - `test/resources/random/dist_api_reference.edn` (utils/fastmath/dev/generate_random_dist_api_reference.py): scipy
    `multivariate_normal` with the default zero means, `norm.sf` in the tail, three-argument `cdf`;
  - `test/resources/random/baseline_reference.edn`: seeded outputs of the code before the audit (unchanged values);
  - the Commons Math parent classes and the SSJ static functions, constructed directly, for the values which have to be unchanged.

  The loops over all keys of the `distribution` multimethod call every function under a `future` with a time limit, so a
  regression of the former hangs (`chi`, `erlang`, `chi-squared-noncentral` with NaN) fails fast instead of blocking the run."
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m]
            [fastmath.protocols :as prot]
            [fastmath.random :as r]
            [fastmath.random.native-distributions-test :as native])
  (:import [fastmath.java.dist FFixed GumbelFixed LevyFixed LogisticFixed NakagamiFixed]
           [org.apache.commons.math3.distribution AbstractRealDistribution FDistribution GumbelDistribution LevyDistribution
            LogisticDistribution NakagamiDistribution]
           [org.apache.commons.math3.random RandomGenerator]
           [umontreal.ssj.probdist InverseGaussianDist]))

(set! *warn-on-reflection* true)

;; ---------------------------------------------------------------- reference data and helpers

(defn- read-resource [name] (edn/read-string (slurp (io/resource (str "random/" name)))))

(def ^:private reference (delay (read-resource "dist_api_reference.edn")))
(def ^:private override-reference (delay (read-resource "dist_override_reference.edn")))
(def ^:private native-reference (delay (read-resource "native_distributions_reference.edn")))
(def ^:private ccdf-reference (delay (read-resource "truncated_mixture_reference.edn")))
(def ^:private baseline (delay (read-resource "baseline_reference.edn")))

(defn- nan? [x] (Double/isNaN (double x)))

(defn- rel-eq?
  "Relative equality `abs(a - b) <= tol * abs(b)`, for `b = 0` the absolute one."
  [a b tol]
  (m/delta-eq (double a) (double b) 0.0 tol))

(defn- bits=
  "Bitwise equality of two doubles (or of two sequences of doubles)."
  [a b]
  (if (and (number? a) (number? b))
    (== (Double/doubleToLongBits (double a)) (Double/doubleToLongBits (double b)))
    (and (= (count a) (count b)) (every? true? (map bits= a b)))))

(defn- same-value?
  "`==` for numbers (`0.0` equals `-0.0`, the infinities must match), two NaN are the same value."
  [a b]
  (or (and (nan? a) (nan? b)) (== (double a) (double b))))

(defn- attempt
  "`[:ok value]` or `[:throw class message]`."
  [f]
  (try [:ok (f)] (catch Throwable t [:throw (class t) (ex-message t)])))

(defn- same-outcome?
  "Two results of `attempt`: the same value (NaN equals NaN) or the same exception class and message."
  [a b]
  (and (= (first a) (first b))
       (if (= :ok (first a))
         (let [x (second a) y (second b)]
           (if (and (number? x) (number? y)) (same-value? x y) (= x y)))
         (= (rest a) (rest b)))))

(defn- timed
  "Calls `f` under a `future`: `[:ok value]`, `[:throw class message data]` or `[:timeout]` after `timeout-ms`."
  [timeout-ms f]
  (let [fu (future (try [:ok (f)] (catch Throwable t [:throw (class t) (ex-message t) (ex-data t)])))
        v (deref fu timeout-ms [:timeout])]
    (when (= [:timeout] v) (future-cancel fu))
    v))

(defn- thrown-input?
  "True when `f` throws `ex-info` of the type `:invalid-input` with the given `:property` (any when `property` is `nil`)."
  [property f]
  (try (f) false
       (catch clojure.lang.ExceptionInfo e
         (let [data (ex-data e)]
           (and (= :invalid-input (:type data)) (or (nil? property) (= property (:property data))))))
       (catch Throwable _ false)))

(def ^:private all-keys (vec (sort (keys (methods r/distribution)))))

(def ^:private multivariate-keys #{:multi-normal :dirichlet :multinomial})

(def ^:private non-numeric-keys
  "Keys whose argument is not a number (the outcomes of the categorical distribution are arbitrary objects)."
  #{:categorical :categorical-distribution})

(def ^:private univariate-keys (vec (remove (into multivariate-keys non-numeric-keys) all-keys)))

(defn- mt ^RandomGenerator [seed] (r/rng :mersenne seed))

;; ---------------------------------------------------------------- 5.1 wiring of the Commons overrides

(def ^:private override-families
  "Per family: constructor, override class, Commons parent class, translation of the reference parameters to the options of
  the constructor and the Commons parent built directly from the reference parameters."
  {:logistic {:ctor r/logistic :fixed LogisticFixed :parent LogisticDistribution
              :opts identity
              :make-parent (fn ^AbstractRealDistribution [^RandomGenerator g {:keys [mu s]}]
                             (LogisticDistribution. g (double mu) (double s)))
              :moments [#(.getNumericalMean ^AbstractRealDistribution %)]}
   :gumbel {:ctor r/gumbel :fixed GumbelFixed :parent GumbelDistribution
            :opts identity
            :make-parent (fn ^AbstractRealDistribution [^RandomGenerator g {:keys [mu beta]}]
                           (GumbelDistribution. g (double mu) (double beta)))
            :moments [#(.getNumericalVariance ^AbstractRealDistribution %)]}
   :nakagami {:ctor r/nakagami :fixed NakagamiFixed :parent NakagamiDistribution
              :opts identity
              :make-parent (fn ^AbstractRealDistribution [^RandomGenerator g {:keys [mu omega]}]
                             (NakagamiDistribution. g (double mu) (double omega)
                                                    NakagamiDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY))
              :moments [#(.getNumericalMean ^AbstractRealDistribution %)
                        #(.getNumericalVariance ^AbstractRealDistribution %)]}
   :f {:ctor r/f :fixed FFixed :parent FDistribution
       :opts (fn [{:keys [d1 d2]}] {:numerator-degrees-of-freedom d1 :denominator-degrees-of-freedom d2})
       :make-parent (fn ^AbstractRealDistribution [^RandomGenerator g {:keys [d1 d2]}]
                      (FDistribution. g (double d1) (double d2) FDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY))
       :moments [#(.getNumericalMean ^AbstractRealDistribution %)
                 #(.getNumericalVariance ^AbstractRealDistribution %)]}
   :levy {:ctor r/levy :fixed LevyFixed :parent LevyDistribution
          :opts identity
          :make-parent (fn ^AbstractRealDistribution [^RandomGenerator g {:keys [mu c]}]
                         (LevyDistribution. g (double mu) (double c)))
          :moments [#(.getNumericalMean ^AbstractRealDistribution %)
                    #(.getNumericalVariance ^AbstractRealDistribution %)]}})

(defn- make-override
  "The object from the public constructor of `family` for the reference parameters `params`."
  [family params rng]
  (let [{:keys [ctor opts]} (get override-families family)]
    (ctor (assoc (opts params) :rng rng))))

(defn- make-parent
  ^AbstractRealDistribution [family params rng]
  ((get-in override-families [family :make-parent]) rng params))

(deftest override-constructors-return-fixed-classes
  (doseq [[family {:keys [ctor ^Class fixed ^Class parent]}] override-families]
    (testing (name family)
      (is (instance? fixed (ctor)) "default parameters")
      (is (instance? fixed (ctor nil)))
      (is (instance? parent (ctor {})) "a subclass of the Commons class")
      (is (= fixed (class (ctor {}))))
      (doseq [{:keys [params]} (get @override-reference family)]
        (is (= fixed (class (make-override family params nil))) (str params)))
      (let [d (ctor)]
        (is (identical? d (r/source-object d)) "source-object of a Commons key is the object itself")
        (is (= family (r/distribution-id d)))))))

(deftest override-values-against-scipy
  (testing "variance of logistic and mean of gumbel (relative 1e-12), the other moment is unchanged"
    (doseq [{:keys [params variance mean]} (:logistic @override-reference)
            :let [d (make-override :logistic params nil)]]
      (is (rel-eq? (r/variance d) variance 1.0e-12) (str "logistic variance " params))
      (is (same-value? mean (r/mean d)) (str "logistic mean " params)))
    (doseq [{:keys [params variance mean]} (:gumbel @override-reference)
            :let [d (make-override :gumbel params nil)]]
      (is (rel-eq? (r/mean d) mean 1.0e-12) (str "gumbel mean " params))
      (is (rel-eq? (r/variance d) variance 1.0e-12) (str "gumbel variance " params))))
  (testing "values outside of the support (relative 1e-9; the reference values are 0, -Inf, exactly)"
    (doseq [{:keys [params xs cdf]} (:nakagami @override-reference)
            :let [d (make-override :nakagami params nil)]]
      (is (every? true? (map (fn [x e] (rel-eq? (r/cdf d x) e 1.0e-9)) xs cdf)) (str "nakagami cdf " params))
      (is (every? #(m/zero? (r/pdf d %)) xs) (str "nakagami pdf " params))
      (is (every? #(m/neg-inf? (r/lpdf d %)) xs) (str "nakagami lpdf " params)))
    (doseq [family [:f :levy]
            {:keys [params xs pdf logpdf cdf]} (get @override-reference family)
            :let [d (make-override family params nil)]]
      (is (every? true? (map (fn [x e] (same-value? e (r/pdf d x))) xs pdf)) (str family " pdf " params))
      (is (every? true? (map (fn [x e] (same-value? e (r/lpdf d x))) xs logpdf)) (str family " lpdf " params))
      (is (every? true? (map (fn [x e] (same-value? e (r/cdf d x))) xs cdf)) (str family " cdf " params))
      (is (every? #(not (nan? (r/probability d %))) xs) (str family " probability " params)))))

(defn- interior-points
  "Quantiles of the Commons parent in the open support (the points at the lower bound are the overridden cases)."
  [^AbstractRealDistribution parent probabilities]
  (let [lb (.getSupportLowerBound parent)]
    (vec (for [p probabilities
               :let [x (try (.inverseCumulativeProbability parent (double p)) (catch Exception _ ##NaN))]
               :when (and (Double/isFinite x) (> x lb))]
           x))))

(def ^:private interior-probabilities (mapv #(/ (+ 0.5 (double %)) 200.0) (range 200)))

(deftest override-interior-equals-commons-parent
  (doseq [[family {:keys [moments]}] override-families
          {:keys [params]} (get @override-reference family)
          :let [d (make-override family params (mt 1))
                parent (make-parent family params (mt 1))
                xs (interior-points parent interior-probabilities)]]
    (testing (str (name family) " " params)
      (is (>= (count xs) 150))
      (is (every? #(bits= (.density parent (double %)) (r/pdf d %)) xs) "pdf")
      (is (every? #(bits= (.logDensity parent (double %)) (r/lpdf d %)) xs) "lpdf")
      (is (every? #(bits= (.cumulativeProbability parent (double %)) (r/cdf d %)) xs) "cdf")
      (is (every? #(bits= (.cumulativeProbability parent (double %) (+ (double %) 1.0))
                          (r/cdf d % (+ (double %) 1.0))) xs) "cdf, three arguments")
      (is (every? (fn [p] (same-outcome? (attempt #(.inverseCumulativeProbability parent (double p)))
                                         (attempt #(r/icdf d p))))
                  interior-probabilities) "icdf")
      (is (every? #(bits= (% parent) (% d)) moments) "the moment which is not overridden")
      (is (bits= (.getSupportLowerBound parent) (r/lower-bound d)))
      (is (bits= (.getSupportUpperBound parent) (r/upper-bound d))))))

(deftest override-seeded-samples-equal-commons-parent
  (doseq [[family _] override-families
          :let [entries (get @override-reference family)]
          {:keys [params]} [(first entries) (last entries)]
          :let [d (make-override family params (mt 20260507))
                ^doubles expected (.sample (make-parent family params (mt 20260507)) 2000)
                ^doubles actual (double-array (repeatedly 2000 #(r/sample d)))]]
    (testing (str (name family) " " params)
      (is (java.util.Arrays/equals expected actual)))))

(deftest override-baseline-values-unchanged
  (let [values (:dist-values @baseline)
        check (fn [k d]
                (let [{:keys [xs pdf cdf]} (get values k)]
                  (is (bits= pdf (mapv #(r/pdf d %) xs)) (str k " pdf"))
                  (is (bits= cdf (mapv #(r/cdf d %) xs)) (str k " cdf"))))]
    (check :logistic (r/logistic {:mu 0.5 :s 2.0}))
    (check :logistic-default (r/logistic))
    (check :gumbel (r/gumbel {:mu 1.0 :beta 2.0}))
    (check :nakagami (r/nakagami {:mu 1.5 :omega 2.0}))
    (check :f (r/f {:numerator-degrees-of-freedom 5.0 :denominator-degrees-of-freedom 7.0}))
    (check :levy (r/levy {:mu 0.1 :c 1.5}))
    (is (bits= (get-in values [:logistic :mean]) (r/mean (r/logistic {:mu 0.5 :s 2.0}))) "mean is unchanged")
    (is (bits= (get-in values [:gumbel :variance]) (r/variance (r/gumbel {:mu 1.0 :beta 2.0}))) "variance is unchanged")
    (is (rel-eq? (r/variance (r/logistic)) (get-in values [:logistic-default :variance]) 1.0e-15)
        "s = 1: the corrected variance pi^2/3 equals the old value")))

;; ---------------------------------------------------------------- 5.2 inverse-gaussian

(def ^:private inverse-gaussian-ratios [1.0 100.0 354.0 355.0 400.0 1000.0 10000.0])

(deftest inverse-gaussian-cdf-against-mpmath
  (let [cases (:inverse-gaussian @native-reference)]
    (is (every? (fn [ratio] (some #(== ratio (/ (double (:lambda %)) (double (:mu %)))) cases)) inverse-gaussian-ratios)
        "all required ratios lambda/mu are in the reference")
    (doseq [{:keys [mu lambda xs cdf]} cases
            :let [d (r/inverse-gaussian {:mu mu :lambda lambda})]]
      (testing (str "mu=" mu " lambda=" lambda)
        (doseq [[x c] (map vector xs cdf)]
          (is (rel-eq? (r/cdf d x) c 1.0e-9) (str "cdf(" x ")")))))))

(deftest inverse-gaussian-round-trip-and-samples
  (doseq [ratio inverse-gaussian-ratios
          :let [d (r/inverse-gaussian {:mu 1.0 :lambda ratio :rng (mt 5)})
                ps [1.0e-9 1.0e-6 0.001 0.1 0.5 0.9 0.999 (- 1.0 1.0e-9)]
                xs (mapv #(r/icdf d %) ps)]]
    (testing (str "lambda/mu=" ratio)
      (is (apply < xs) "increasing quantiles, not degenerate")
      (is (every? true? (map (fn [p x] (<= (m/abs (- p (r/cdf d x))) 1.0e-9)) ps xs))
          (str "max |cdf(icdf p) - p| " (apply max (map (fn [p x] (m/abs (- p (r/cdf d x)))) ps xs))))
      (let [n 20000
            samples (vec (r/->seq d n))
            mean-v (/ (reduce + samples) (double n))]
        (is (> (count (distinct samples)) (* 0.99 n)) "not degenerate")
        (is (every? #(and (pos? %) (m/valid-double? %)) samples))
        (is (< (m/abs (- mean-v 1.0)) (* 5.0 (/ (m/sqrt (r/variance d)) (m/sqrt n)))) "sample mean within 5 standard errors")))))

(deftest inverse-gaussian-small-ratio-equals-ssj
  (testing "lambda/mu <= 100: the SSJ cdf and inverse cdf, bitwise"
    (doseq [[mu lambda] [[1.0 1.0] [1.0 50.0] [1.0 100.0] [2.0 200.0] [100.0 1.0] [3.0 30.0]]
            :let [d (r/inverse-gaussian {:mu mu :lambda lambda})
                  xs (mapv #(* mu %) [0.05 0.2 0.5 0.8 1.0 1.2 2.0 5.0])]]
      (testing (str "mu=" mu " lambda=" lambda)
        (is (every? #(bits= (InverseGaussianDist/cdf mu lambda (double %)) (r/cdf d %)) xs) "cdf")
        (is (every? #(bits= (InverseGaussianDist/inverseF mu lambda (double %)) (r/icdf d %)) [0.001 0.1 0.5 0.9 0.999]) "icdf"))))
  (testing "seeded samples and interior values equal the code before the audit"
    (let [samples (:dist-samples @baseline)
          v (get-in @baseline [:dist-values :inverse-gaussian-1-50])
          d (r/inverse-gaussian {:mu 1.0 :lambda 50.0})]
      (is (= (:inverse-gaussian-1-50 samples)
             (vec (r/->seq (r/inverse-gaussian {:mu 1.0 :lambda 50.0 :rng (mt 9)}) 10))))
      (is (= (:inverse-gaussian-2-100 samples)
             (vec (r/->seq (r/inverse-gaussian {:mu 2.0 :lambda 100.0 :rng (mt 9)}) 10))))
      (is (bits= (:pdf v) (mapv #(r/pdf d %) (:xs v))))
      (is (bits= (:cdf v) (mapv #(r/cdf d %) (:xs v)))))))

;; ---------------------------------------------------------------- 5.3 hardening of the proxies

(def ^:private limit-cases
  "[name proxy [[argument expected]]] for the univariate distributions."
  [["cdf" r/cdf [[##NaN ##NaN] [##Inf 1.0] [##-Inf 0.0]]]
   ["pdf" r/pdf [[##NaN ##NaN] [##Inf 0.0] [##-Inf 0.0]]]
   ["lpdf" r/lpdf [[##NaN ##NaN] [##Inf ##-Inf] [##-Inf ##-Inf]]]
   ["probability" r/probability [[##NaN ##NaN] [##Inf 0.0] [##-Inf 0.0]]]
   ["ccdf" r/ccdf [[##NaN ##NaN] [##Inf 0.0] [##-Inf 1.0]]]])

(deftest nan-and-infinity-all-univariate-keys
  (is (> (count univariate-keys) 120))
  (let [start (System/nanoTime)
        slow (atom [])]
    (doseq [k univariate-keys
            :let [d (r/distribution k)]
            [fname f expected] limit-cases]
      (let [t0 (System/nanoTime)
            failures (vec (for [[x e] expected
                                :let [o (timed 5000 #(f d x))]
                                :when (not (and (= :ok (first o)) (same-value? e (second o))))]
                            [x o]))
            ms (/ (double (- (System/nanoTime) t0)) 1.0e6)]
        (when (> ms 1000.0) (swap! slow conj [k fname ms]))
        (is (empty? failures) (str k " " fname ": " (pr-str failures)))))
    (let [seconds (/ (double (- (System/nanoTime) start)) 1.0e9)]
      (is (< seconds 60.0) (str "the loop over " (count univariate-keys) " keys took " seconds " s"))
      (is (empty? @slow) (str "calls slower than 1 s: " @slow)))))

(deftest former-hangs-return-quickly
  (doseq [[k f x expected] [[:chi r/cdf ##NaN ##NaN]
                            [:chi r/ccdf ##NaN ##NaN]
                            [:erlang r/cdf ##NaN ##NaN]
                            [:erlang r/ccdf ##NaN ##NaN]
                            [:chi-squared-noncentral r/cdf ##NaN ##NaN]
                            [:logarithmic r/cdf ##Inf 1.0]
                            [:logarithmic r/cdf ##-Inf 0.0]
                            [:nakagami r/cdf ##Inf 1.0]
                            [:gamma r/cdf ##Inf 1.0]
                            [:gamma r/ccdf ##Inf 0.0]
                            [:chi-squared r/cdf ##Inf 1.0]
                            [:zaga r/cdf ##Inf 1.0]
                            [:pascal r/pdf ##Inf 0.0]
                            [:watson-g r/cdf ##NaN ##NaN]
                            [:f r/cdf ##Inf 1.0]
                            [:inverse-gaussian r/cdf ##Inf 1.0]
                            [:pearson-6 r/cdf ##Inf 1.0]]
          :let [o (timed 1000 #(f (r/distribution k) x))]]
    (is (= :ok (first o)) (str k " " x " " (pr-str o)))
    (is (and (= :ok (first o)) (same-value? expected (second o))) (str k " " x " expected " expected))))

(deftest nan-and-infinity-discrete-and-integer-arguments
  (testing "integers (Long) are numbers too"
    (let [d (r/poisson {:p 4.0})]
      (is (= 1.0 (r/cdf d Long/MAX_VALUE)))
      (is (= 0.0 (r/cdf d Long/MIN_VALUE)))
      (is (= 0.0 (r/pdf d Long/MAX_VALUE)))
      (is (= ##-Inf (r/lpdf d Long/MIN_VALUE)))
      (is (= 0.0 (r/ccdf d Long/MAX_VALUE)))
      (is (= 1.0 (r/ccdf d Long/MIN_VALUE)))))
  (testing "the pmf of a discrete distribution at NaN"
    (is (nan? (r/pdf (r/poisson) ##NaN)))
    (is (nan? (r/probability (r/binomial) ##NaN)))
    (is (nan? (r/lpdf (r/geometric) ##NaN)))))

(defn- sample-point
  "A valid argument for the multivariate distribution."
  [d]
  (vec (r/sample d)))

(deftest vector-arguments-pass-through
  (doseq [k (sort multivariate-keys)
          :let [d (r/distribution k)
                ok (sample-point d)
                points [ok [] nil [1.0] [1.0 2.0 3.0] (vec (repeat 5 0.5)) [##NaN ##NaN] [##Inf ##-Inf]]]
          [fname f p] [["cdf" r/cdf prot/cdf] ["pdf" r/pdf prot/pdf] ["lpdf" r/lpdf prot/lpdf]
                       ["probability" r/probability prot/probability]]
          v points]
    (is (same-outcome? (attempt #(p d v)) (attempt #(f d v))) (str k " " fname " " v)))
  (testing "ccdf is 1 - cdf for vectors"
    (doseq [k (sort multivariate-keys)
            :let [d (r/distribution k)
                  ok (sample-point d)]
            v [ok [] nil [1.0]]
            :let [expected (attempt #(- 1.0 (double (prot/cdf d v))))
                  actual (attempt #(r/ccdf d v))]]
      (is (= (first expected) (first actual)) (str k " ccdf " v))
      (when (= :ok (first expected)) (is (bits= (second expected) (second actual)) (str k " ccdf " v)))
      (when (= :throw (first expected)) (is (= (rest expected) (rest actual)) (str k " ccdf " v)))))
  (testing "the outcomes of the categorical distribution are not numbers"
    (let [d (r/categorical {:data [:a :b :c] :probabilities [0.2 0.3 0.5]})]
      (doseq [[fname f p] [["cdf" r/cdf prot/cdf] ["pdf" r/pdf prot/pdf] ["lpdf" r/lpdf prot/lpdf]
                           ["probability" r/probability prot/probability]]
              v [:a :b :c :z nil "a" [:a]]]
        (is (same-outcome? (attempt #(p d v)) (attempt #(f d v))) (str fname " " v)))
      (is (= 0.3 (r/pdf d :b)))
      (is (= ##-Inf (r/lpdf d :z))))))

(deftest icdf-out-of-range
  (let [bad [-0.1 1.1 ##NaN -5.0e-324 1.0000000000000002 ##Inf ##-Inf -1.0 2.0 Double/MAX_VALUE]]
    (doseq [k all-keys
            :let [d (r/distribution k)]]
      (testing k
        (doseq [p bad
                :let [o (timed 5000 #(r/icdf d p))]]
          (is (= :throw (first o)) (str "icdf " p " " (pr-str o)))
          (is (= clojure.lang.ExceptionInfo (second o)) (str "icdf " p " " (pr-str o)))
          (let [data (nth o 3 nil)]
            (is (= :invalid-input (:type data)) (str "icdf " p " " data))
            (is (= :p (:property data)) (str "icdf " p " " data))))))))

(deftest icdf-end-points-pass-through
  (doseq [k [:normal :gamma :beta :poisson :binomial :exponential :t :weibull :mixture :integer-discrete :truncated :uniform-real]
          :let [d (r/distribution k)]
          p [0.0 -0.0 1.0 0.5 1.0e-300 0.9999999999999999]]
    (is (same-outcome? (attempt #(prot/icdf d p)) (attempt #(r/icdf d p))) (str k " icdf " p)))
  (testing "values"
    (is (= ##-Inf (r/icdf (r/normal) 0.0)))
    (is (= ##Inf (r/icdf (r/normal) 1.0)))
    (is (= 0.0 (r/icdf (r/normal) 0.5)))
    (is (= 0 (r/icdf (r/poisson) 0.0))))
  (testing "integers and ratios are numbers"
    (is (= ##-Inf (r/icdf (r/normal) 0)))
    (is (= ##Inf (r/icdf (r/normal) 1)))
    (is (= 0.0 (apply r/icdf [(r/normal) 1/2]))))
  (testing "nil is not a probability"
    (doseq [k [:normal :poisson :mixture]]
      (is (= :throw (first (attempt #(apply r/icdf [(r/distribution k) nil])))) (str k)))))

(def ^:private three-arity-keys [:normal :gamma :beta :poisson :binomial :exponential :laplace :weibull :logistic :rayleigh])

(deftest cdf-three-arguments
  (testing "P(v1 < X <= v2) = cdf(v2) - cdf(v1), unchanged by the hardening"
    (doseq [k three-arity-keys
            :let [d (r/distribution k)
                  v1 (double (r/icdf d 0.2))
                  v2 (double (r/icdf d 0.8))]]
      (testing k
        (is (< v1 v2))
        (is (m/delta-eq (r/cdf d v1 v2) (- (r/cdf d v2) (r/cdf d v1)) 1.0e-15))
        (is (m/delta-eq 0.0 (r/cdf d v1 v1) 1.0e-15)))))
  (testing "scipy"
    (doseq [{:keys [key params ranges]} (:cdf-range @reference)
            :let [d (r/distribution key params)]
            [v1 v2 e] ranges]
      ;; Commons Math cdf of the discrete families is accurate to ~1e-14 only (observed: poisson(4) 3.5e-14)
      (is (m/delta-eq (r/cdf d v1 v2) e 1.0e-12) (str key " " v1 " " v2)))))

(deftest boundary-numeric-arguments
  (testing "zero, one, minus one, the smallest and the largest numbers, Long limits: always a number in range"
    (let [xs [0.0 1.0 -1.0 4.9e-324 Double/MAX_VALUE (- Double/MAX_VALUE) Long/MAX_VALUE Long/MIN_VALUE 0 1 -1 1/2 -0.0]
          ;; [key function (double argument)] of the defects of the backend (Commons Math 3.6.1) at finite extreme arguments,
          ;; not covered by the specification (5.3 hardens NaN and the infinities only, 5.9 WONT):
          ;; logistic pdf: exp(-x) overflows; levy pdf: 0 * Inf at the smallest denormal; nakagami: the incomplete gamma
          ;; function does not converge at the largest double (cdf throws, pdf is NaN); weibull log-density: log(0) - Inf at 0
          backend-defects {[:logistic :pdf] #(< % -700.0)
                           [:logistic :lpdf] #(< % -700.0)
                           [:levy :pdf] #(== % 4.9e-324)
                           [:levy :lpdf] #(== % 4.9e-324)
                           [:nakagami :pdf] #(== % Double/MAX_VALUE)
                           [:nakagami :lpdf] #(== % Double/MAX_VALUE)
                           [:nakagami :cdf] #(== % Double/MAX_VALUE)
                           [:nakagami :ccdf] #(== % Double/MAX_VALUE)
                           [:weibull :lpdf] #(== % 0.0)}]
      (doseq [k [:normal :gamma :beta :poisson :binomial :exponential :t :weibull :mixture :truncated :integer-discrete
                 :logistic :f :zip :negative-binomial :logarithmic :cauchy :levy :nakagami :uniform-real]
              :let [d (r/distribution k)]
              x xs
              :let [dx (double x)
                    defect? (fn [fname] ((get backend-defects [k fname] (constantly false)) dx))
                    c (timed 5000 #(r/cdf d x))
                    cc (timed 5000 #(r/ccdf d x))
                    p (timed 5000 #(r/pdf d x))
                    lp (timed 5000 #(r/lpdf d x))]]
        (when-not (defect? :cdf)
          (is (and (= :ok (first c)) (<= 0.0 (second c) 1.0)) (str k " cdf " x " " (pr-str c))))
        (when-not (defect? :ccdf)
          (is (and (= :ok (first cc)) (<= 0.0 (second cc) 1.0)) (str k " ccdf " x " " (pr-str cc))))
        (when-not (defect? :pdf)
          (is (and (= :ok (first p)) (<= 0.0 (second p))) (str k " pdf " x " " (pr-str p))))
        (when-not (defect? :lpdf)
          ;; the density of the f distribution at 0 is +Inf (d1 = 1), so is its log
          (is (and (= :ok (first lp)) (not (nan? (second lp)))) (str k " lpdf " x " " (pr-str lp))))))))

(deftest boundary-non-numeric-arguments
  (testing "nil, a string and vectors are not arguments of a univariate distribution: an exception, no value, no hang"
    (doseq [k [:normal :gamma :beta :poisson :binomial :mixture :truncated :integer-discrete :logarithmic :f :zip]
            :let [d (r/distribution k)]
            [fname f] [["cdf" r/cdf] ["pdf" r/pdf] ["lpdf" r/lpdf] ["probability" r/probability] ["ccdf" r/ccdf]]
            v [nil "a" [] [1.0] [1.0 2.0] :a]
            :let [o (timed 5000 #(f d v))]]
      (is (= :throw (first o)) (str k " " fname " " (pr-str v) " " (pr-str o)))))
  (testing "the functions without an argument on nil, a number and an empty vector"
    (doseq [[fname f] [["sample" r/sample] ["dimensions" r/dimensions] ["source-object" r/source-object]
                       ["continuous?" r/continuous?] ["discrete?" r/discrete?] ["mean" r/mean] ["variance" r/variance]
                       ["lower-bound" r/lower-bound] ["upper-bound" r/upper-bound] ["distribution-id" r/distribution-id]]
            v [nil 0 [] 1.5]
            :let [o (timed 5000 #(f v))]]
      (is (= :throw (first o)) (str fname " " (pr-str v) " " (pr-str o)))))
  (testing "distribution? is false for everything which is not a distribution"
    (doseq [v [nil 0 1.5 [] {} "normal" :normal (mt 1) ##NaN]]
      (is (false? (r/distribution? v)) (pr-str v)))
    (is (true? (r/distribution? (r/normal))))
    (is (true? (r/distribution? (r/multi-normal))))
    (is (true? (r/distribution? (r/categorical)))))
  (testing "means and covariance are for the multivariate distributions only"
    (doseq [k [:normal :poisson :gamma :mixture]]
      (is (= :throw (first (attempt #(r/means (r/distribution k))))) (str k))
      (is (= :throw (first (attempt #(r/covariance (r/distribution k))))) (str k)))
    (is (= [0.0 0.0] (vec (r/means (r/multi-normal)))))))

(deftest likelihood-boundary
  (let [d (r/normal)]
    (is (= 0.0 (r/log-likelihood d [])))
    (is (= 0.0 (r/log-likelihood d nil)))
    (is (= 1.0 (r/likelihood d [])))
    (is (= 0.0 (r/observe d [])))
    (is (nan? (r/log-likelihood d [0.0 ##NaN 1.0])))
    (is (= ##-Inf (r/log-likelihood d [0.0 ##Inf])))
    (is (= ##-Inf (r/observe d [##-Inf 0.0])))
    (is (= 0.0 (r/likelihood d [##Inf])))
    (is (nan? (r/likelihood d [##NaN])))
    (is (bits= (r/lpdf d 0.3) (r/observe1 d 0.3)))
    (is (nan? (r/observe1 d ##NaN)))
    (is (= ##-Inf (r/observe1 d ##Inf)))
    (is (m/delta-eq (+ (r/lpdf d 0.1) (r/lpdf d 0.2)) (r/log-likelihood d [0.1 0.2]) 1.0e-15))
    (is (m/delta-eq (* (r/pdf d 0.1) (r/pdf d 0.2)) (r/likelihood d [0.1 0.2]) 1.0e-15)))
  (testing "multivariate observations are vectors, a wrong length is an error"
    (let [d (r/multi-normal)]
      (is (m/delta-eq (+ (r/lpdf d [0.0 0.0]) (r/lpdf d [1.0 1.0])) (r/log-likelihood d [[0.0 0.0] [1.0 1.0]]) 1.0e-15))
      (is (= :throw (first (attempt #(r/log-likelihood d [[0.0]]))))))))

;; ---------------------------------------------------------------- 5.4 ccdf

(def ^:private ccdf-closed-form-keys
  [:normal :t :logistic :laplace :cauchy :exponential :weibull :pareto :log-normal :gamma :chi-squared])

(def ^:private ccdf-fallback-keys
  "Keys without a survival function: `ccdf` is `1 - cdf`."
  [:beta :f :gumbel :levy :nakagami :triangular :uniform-real :folded-normal :inverse-gaussian
   :poisson :binomial :geometric :negative-binomial :logarithmic])

(defn- accurate-levels
  "The levels at which `ccdf` is expected to be accurate: the SSJ `barF` of the Pearson VI loses digits below 1e-12
  (reported by the truncated/mixture tests)."
  [k]
  (case k
    :pearson-6 [1.0e-6 1.0e-12]
    [1.0e-6 1.0e-12 1.0e-30 1.0e-100]))

(deftest ccdf-against-scipy
  (let [rows (remove #(= :inverse-gaussian (:key %)) (:ccdf @ccdf-reference))]
    (is (every? (set (map :key rows)) ccdf-closed-form-keys) "every closed form has reference rows")
    (is (>= (count rows) 25))
    (doseq [{:keys [key params points]} rows
            :let [d (r/distribution key params)]]
      (testing (str key " " params)
        (doseq [[level x sf] points
                :when (some #{level} (accurate-levels key))]
          ;; the observed relative error is below 1e-13 (normal at 1e-100); the tolerance is the one of the specification
          (is (rel-eq? (r/ccdf d x) sf 1.0e-9) (str key " sf level " level " x " x)))))))

(deftest ccdf-normal-tail
  (let [{:keys [xs sf]} (:normal-sf @reference)
        d (r/normal)]
    (doseq [[x e] (map vector xs sf)]
      (is (rel-eq? (r/ccdf d x) e 1.0e-9) (str "ccdf(normal, " x ")")))
    (is (rel-eq? (r/ccdf d 8.3) 5.2055697448902465e-17 1.0e-9) "the old 1 - cdf was exactly 0.0")
    (is (= 1.0 (r/ccdf d -8.3)))
    (is (pos? (r/ccdf d 37.0)))
    (is (= 0.0 (r/ccdf d ##Inf)))
    (is (= 1.0 (r/ccdf d ##-Inf)))))

(deftest ccdf-fallback-is-one-minus-cdf
  (doseq [k ccdf-fallback-keys
          :let [d (r/distribution k)
                xs (concat (keep (fn [p] (try (double (r/icdf d p)) (catch Throwable _ nil))) [0.001 0.1 0.5 0.9 0.999])
                           [0.3 2.5 7.0])]]
    (testing k
      (is (seq xs))
      (doseq [x xs]
        (is (bits= (- 1.0 (r/cdf d x)) (r/ccdf d x)) (str "x = " x))))))

(deftest ccdf-symmetry-and-complement
  (doseq [k ccdf-closed-form-keys
          :let [d (r/distribution k)]
          x [-2.0 -0.5 0.5 1.5 3.0]]
    (is (m/delta-eq (+ (r/cdf d x) (r/ccdf d x)) 1.0 1.0e-14) (str k " " x)))
  (is (m/delta-eq (r/ccdf (r/normal) 1.7) (r/cdf (r/normal) -1.7) 1.0e-15) "symmetry of the normal distribution"))

;; ---------------------------------------------------------------- 5.5 multi-normal

(deftest multi-normal-default-means
  (testing "covariances only: the means are zeros of the dimension of the covariances"
    (let [d (r/multi-normal {:covariances [[1 0 0] [0 1 0] [0 0 1]]})]
      (is (= [0.0 0.0 0.0] (vec (r/means d))))
      (is (= 3 (r/dimensions d))))
    (is (= 1 (r/dimensions (r/multi-normal {:covariances [[2.0]]}))))
    (is (= 2 (r/dimensions (r/multi-normal))))
    (is (= [0.0 0.0] (vec (r/means (r/multi-normal nil)))))
    (is (= [0.0 0.0] (vec (r/means (r/multi-normal {})))))
    (is (= [1.0 2.0 3.0] (vec (r/means (r/multi-normal {:means [1 2 3]})))) "means only: identity covariance")
    (is (= 3 (r/dimensions (r/multi-normal {:means [1 2 3]}))))))

(deftest multi-normal-zero-mean-against-scipy
  (doseq [{:keys [covariances points pdf lpdf]} (:multi-normal-zero-mean @reference)
          :let [d (r/multi-normal {:covariances covariances})]]
    (testing (str (count covariances) " dimensions")
      (doseq [[x p lp] (map vector points pdf lpdf)]
        (is (rel-eq? (r/pdf d x) p 1.0e-12) (str "pdf " x))
        (is (rel-eq? (r/lpdf d x) lp 1.0e-12) (str "lpdf " x))))))

(deftest multi-normal-with-means-against-scipy
  (doseq [{:keys [means covariances points pdf lpdf]} (get-in @native-reference [:multi-normal :interior])
          :let [d (r/multi-normal {:means means :covariances covariances})]]
    (testing (str (count means) " dimensions")
      (doseq [[x p lp] (map vector points pdf lpdf)]
        (is (rel-eq? (r/pdf d x) p 1.0e-12) (str "pdf " x))
        (is (rel-eq? (r/lpdf d x) lp 1.0e-12) (str "lpdf " x))
        (is (bits= (m/log (r/pdf d x)) (r/lpdf d x)) "lpdf is log(pdf) bitwise while the density does not underflow")))))

(deftest multi-normal-lpdf-far-tail
  (let [d (r/multi-normal {:means [0.0 0.0] :covariances [[1.0 0.5] [0.5 2.0]]})]
    (is (= 0.0 (r/pdf d [40.0 -40.0])) "the density underflows")
    (is (rel-eq? (r/lpdf d [40.0 -40.0]) -1830.6891135318056 1.0e-12))
    (is (rel-eq? (r/lpdf d [40 -40]) -1830.6891135318056 1.0e-12) "integer coordinates")
    (is (= 1.0 (r/ccdf (r/normal) -40.0)) "sanity of the helper"))
  (doseq [{:keys [means covariances point lpdf]} (get-in @native-reference [:multi-normal :far])
          :let [d (r/multi-normal {:means means :covariances covariances})]]
    (is (rel-eq? (r/lpdf d point) lpdf 1.0e-12) (str point))
    (is (= 0.0 (r/pdf d point)) (str point))
    (is (bits= (r/lpdf d point) (r/log-likelihood d [point])) "log-likelihood does not underflow")
    (is (m/valid-double? (r/log-likelihood d [point])))))

(deftest multi-normal-baseline-unchanged
  (let [d (r/multi-normal {:means [0.0 1.0] :covariances [[1.0 0.5] [0.5 2.0]]})
        v (get-in @baseline [:dist-values :multi-normal])]
    (is (bits= (:pdf v) (mapv #(r/pdf d %) (:points v))))
    (is (bits= (:lpdf v) (mapv #(r/lpdf d %) (:points v))))))

(deftest multi-normal-validation
  (testing "invalid covariances: ex-info, not Commons exceptions"
    (doseq [[label cov] {"asymmetric" [[1 0.9] [0.1 1]]
                         "ragged" [[1 0] [0]]
                         "not square" [[1 0 0] [0 1 0]]
                         "indefinite" [[1 2] [2 1]]
                         "singular" [[1 1] [1 1]]
                         "zero matrix" [[0 0] [0 0]]
                         "negative variance" [[-1 0] [0 1]]
                         "NaN entry" [[1 ##NaN] [##NaN 1]]
                         "empty" []
                         "a vector" [1 2 3]}]
      (is (thrown-input? nil #(r/multi-normal {:covariances cov})) label)
      (is (thrown-input? nil #(r/multi-normal {:means (vec (repeat (count cov) 0)) :covariances cov})) label)))
  (testing "sizes of means and covariances"
    (is (thrown-input? nil #(r/multi-normal {:means [0 0 0] :covariances [[1 0] [0 1]]})))
    (is (thrown-input? nil #(r/multi-normal {:means [0] :covariances [[1 0] [0 1]]}))))
  (testing "the point"
    (let [d (r/multi-normal {:covariances [[1 0.5] [0.5 2]]})]
      (doseq [x [[1] [1 2 3] nil [] 3.0 "ab"]]
        (is (thrown-input? nil #(r/pdf d x)) (str "pdf " (pr-str x)))
        (is (thrown-input? nil #(r/lpdf d x)) (str "lpdf " (pr-str x))))
      (is (nan? (r/pdf d [##NaN 0.0])))
      (is (nan? (r/lpdf d [##NaN 0.0]))))))

;; ---------------------------------------------------------------- 5.6 id and parameters

(deftest id-and-parameters-all-keys
  (let [mismatches @#'native/known-mismatches
        aliases @#'native/constructor-aliases
        ids @#'native/documented-ids]
    (is (> (count all-keys) 120))
    (is (every? (set all-keys) (keys mismatches)) "known mismatches are keys of the multimethod")
    (is (every? #(and (string? %) (> (count %) 20)) (vals mismatches)) "every exception has a reason")
    (doseq [k all-keys
            :let [v (#'native/constructor-var k)
                  d (r/distribution k)
                  params (r/distribution-parameters d)
                  ctor-keys (some-> v (#'native/constructor-keys))
                  canonical (keyword (name (or (aliases k) (symbol (name k)))))]]
      (testing k
        (is (some? v) "the constructor is found")
        (is (contains? (conj (get ids k #{}) k canonical) (r/distribution-id d)) (str "id " (r/distribution-id d)))
        (is (every? keyword? params))
        (is (= (count params) (count (set params))) "no repeated parameters")
        (cond
          (nil? ctor-keys) (is (#{:continuous-distribution :kde} k) "constructors without destructuring")
          (contains? mismatches k) (is (not= ctor-keys (set params)) (str "fixed? remove " k " from known-mismatches"))
          :else (is (= ctor-keys (set params)) (str "parameters " (sort params) " constructor " (sort ctor-keys))))))))

(deftest id-and-parameters-fixed-entries
  (is (= :half-cauchy (r/distribution-id (r/distribution :half-cauchy))))
  (is (= [:distr :left :right :rng] (r/distribution-parameters (r/truncated))))
  (is (= :truncated (r/distribution-id (r/truncated))))
  (is (= [:distrs :weights :rng] (r/distribution-parameters (r/mixture))))
  (is (= :mixture (r/distribution-id (r/mixture))))
  (doseq [k [:fishers-noncentral-hypergeometric :wallenius-noncentral-hypergeometric]]
    (is (some #{:nf} (r/distribution-parameters (r/distribution k))) (str k))
    (is (not-any? #{:fn} (r/distribution-parameters (r/distribution k))) (str k)))
  (testing "a keyword is accepted"
    (is (= (r/distribution-parameters (r/normal)) (r/distribution-parameters :normal)))
    (is (= (r/distribution-parameters (r/distribution :gev)) (r/distribution-parameters :gev)))))

;; ---------------------------------------------------------------- 5.7 moments at the proxy level

(deftest mean-and-variance-are-doubles
  (doseq [k [:negative-binomial :logarithmic :nbi :beta-binomial :zip :generalized-extreme-value :generalized-pareto
             :f-noncentral :t-noncentral :half-cauchy :inverse-gaussian :integer-discrete]
          :let [d (r/distribution k)]]
    (testing k
      (is (double? (r/mean d)))
      (is (double? (r/variance d)))
      (is (not= ##-Inf (r/variance d)))))
  (testing "undefined moments through the proxies"
    (is (nan? (r/mean (r/distribution :t-noncentral {:df 1.0 :ncp 1.5}))))
    (is (nan? (r/variance (r/distribution :t-noncentral {:df 1.0 :ncp 1.5}))))
    (is (= ##Inf (r/mean (r/distribution :generalized-extreme-value {:mu 0.0 :sigma 1.0 :xi 1.5}))))
    (is (nan? (r/mean (r/distribution :half-cauchy))))
    (is (= ##Inf (r/variance (r/distribution :f-noncentral {:df1 5.0 :df2 4.0 :ncp 1.0}))))))

;; ---------------------------------------------------------------- 5.8 documentation

(def ^:private protocol-proxies
  "The 28 vars of the protocol layer of `fastmath.random`: the 24 functions of the distribution protocols, `set-seed!` and
  the three vars `default-normal`, `distributions-list`, `distribution`."
  '[distribution? cdf ccdf pdf lpdf icdf probability sample dimensions source-object continuous? discrete? observe1
    log-likelihood observe likelihood mean means variance covariance lower-bound upper-bound distribution-id
    distribution-parameters set-seed! default-normal distributions-list distribution])

(defn- doc-of [sym] (:doc (meta (ns-resolve 'fastmath.random sym))))

(def ^:private non-function-proxies
  "The vars among the 28 which are a value or a multimethod: their docstring is a description, it has no \"Returns\" sentence."
  '#{default-normal distributions-list distribution})

(deftest proxy-docstrings
  (is (= 28 (count (distinct protocol-proxies))))
  (doseq [sym protocol-proxies
          :let [doc (doc-of sym)]]
    (testing sym
      (is (string? doc))
      (is (> (count doc) 80) (str "one-liner: " (pr-str doc)))
      (when-not (non-function-proxies sym)
        (is (str/includes? (str doc) "Returns") "the docstring describes the return value"))))
  (testing "specific content of the key proxies"
    (let [doc #(str (doc-of %))]
      (is (str/includes? (doc 'cdf) "v1") "three arguments")
      (is (str/includes? (doc 'cdf) "NaN"))
      (is (str/includes? (doc 'icdf) ":invalid-input"))
      (is (str/includes? (doc 'ccdf) "survival function"))
      (is (str/includes? (doc 'ccdf) "closed forms") "the list of the accurate keys")
      (is (str/includes? (doc 'ccdf) "SSJ"))
      (is (str/includes? (doc 'distribution) ":inverse-abs-accuracy"))
      (is (str/includes? (doc 'distribution) "ignored") "unknown option names are ignored")
      (is (str/includes? (doc 'distribution) "IllegalArgumentException") "unknown key"))))

(defn- ns-docs
  "[var-or-namespace-name docstring] of the namespace and of all its interned vars."
  [ns-sym]
  (cons [(str ns-sym) (:doc (meta (the-ns ns-sym)))]
        (for [[sym v] (ns-interns (the-ns ns-sym))] [(str sym) (:doc (meta v))])))

(deftest no-broken-links
  (doseq [[nm doc] (ns-docs 'fastmath.random)
          broken ["[[RNGProto]]" "[[DistributionProto]]" "[[integrate-pdf]]"]]
    (is (not (str/includes? (str doc) broken)) (str nm " contains " broken))))

(def ^:private external-link-allow-list
  "Targets of `[[...]]` which are not vars of `fastmath.random` or `fastmath.random.distributions` and are valid external
  references (the Commons Math and SSJ classes named in the docstrings of the backed constructors, Clojure core)."
  #{})

(defn- link-targets [doc] (map second (re-seq #"\[\[([^\[\]]+)\]\]" (str doc))))

(defn- resolves?
  "A target is valid when it is a public var of `fastmath.random` or `fastmath.random.distributions`, a `fastmath.*`
  qualified name or an alias-qualified name of an existing public var."
  [target]
  (let [target (str/trim target)]
    (cond
      (contains? external-link-allow-list target) true
      (str/starts-with? target "fastmath.") true
      (str/includes? target "/")
      (let [[a n] (str/split target #"/" 2)
            the-ns' (or (get (ns-aliases (the-ns 'fastmath.random)) (symbol a)) (find-ns (symbol a)))]
        (boolean (and the-ns' (get (ns-publics the-ns') (symbol n)))))
      :else (boolean (or (get (ns-publics (the-ns 'fastmath.random)) (symbol target))
                         (get (ns-publics (the-ns 'fastmath.random.distributions)) (symbol target)))))))

(defn- public-docs
  "[name docstring] of the namespace and of its public vars."
  [ns-sym]
  (cons [(str ns-sym) (:doc (meta (the-ns ns-sym)))]
        (for [[sym v] (ns-publics (the-ns ns-sym))] [(str sym) (:doc (meta v))])))

(defn- unresolved-links [ns-sym]
  (vec (distinct (for [[nm doc] (public-docs ns-sym)
                       target (link-targets doc)
                       :when (not (resolves? target))]
                   [nm target]))))

(deftest links-resolve
  (let [unresolved (unresolved-links 'fastmath.random)]
    (is (empty? unresolved) (str (count unresolved) " unresolved links: " (pr-str (take 40 unresolved)))))
  (is (> (count (mapcat (comp link-targets second) (public-docs 'fastmath.random))) 1000) "the resolver sees the links"))

(def ^:private known-unresolved-distributions-links
  "Links of the public docstrings of `fastmath.random.distributions` to private helper functions (not documented by
  Codox: dead links in the generated docs), found by the resolver. Remove an entry when the docstring is fixed (this
  test fails for a fixed one)."
  #{["inverse-gaussian" "inverse-gaussian-cdf"]
    ["normal-inverse-gaussian" "generalized-hyperbolic-core"]})

(deftest links-resolve-in-distributions
  (is (= known-unresolved-distributions-links (set (unresolved-links 'fastmath.random.distributions)))))

(deftest distributions-namespace-documented
  (is (string? (:doc (meta (the-ns 'fastmath.random.distributions)))))
  (is (> (count (:doc (meta (the-ns 'fastmath.random.distributions)))) 200))
  (let [undocumented (vec (sort (for [[sym v] (ns-publics (the-ns 'fastmath.random.distributions))
                                      :when (str/blank? (:doc (meta v)))]
                                  sym)))]
    (is (empty? undocumented) (str (count undocumented) " public vars without a docstring: " undocumented))))

(deftest list-vars-documented-as-snapshots
  (doseq [sym '[distributions-list]]
    (is (re-find #"(?i)snapshot|load time|loaded" (doc-of sym)) (str sym))))

(deftest distributions-list-is-the-multimethod-keys
  (is (= (set (keys (methods r/distribution))) (set r/distributions-list)))
  (is (= (sort r/distributions-list) (seq r/distributions-list)) "sorted")
  (is (> (count r/distributions-list) 120)))
