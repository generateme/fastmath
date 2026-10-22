(ns fastmath.random.native-distributions-test
  "Tests of the native distributions of `fastmath.random.distributions` (Random Audit, work packages 5.2, 5.5, 5.6, 5.7,
  7.1, 7.2, 7.4, 7.5, 7.6, 7.7 and 7.9).

  The oracle is `test/resources/random/native_distributions_reference.edn`
  (utils/fastmath/dev/generate_random_native_distributions_reference.py: mpmath for `logarithmic` and `inverse-gaussian`, scipy for
  `dirichlet`, `multinomial` and `multivariate_normal`). The seeded outputs of the code before the audit are
  `test/resources/random/baseline_reference.edn`.

  The functions of the protocols (`fastmath.protocols`) are called directly, the NaN and infinity handling of the user facing
  functions of `fastmath.random` is not a part of these tests."
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.set :as set]
            [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m]
            [fastmath.protocols :as prot]
            [fastmath.random :as r]
            [fastmath.random.distributions :as distr]))

(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "random/native_distributions_reference.edn")))))

(def ^:private baseline
  (delay (edn/read-string (slurp (io/resource "random/baseline_reference.edn")))))

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

(defn- thrown-input?
  "True when `f` throws `ex-info` of the type `:invalid-input` with the given `:property` (when `property` is not `nil`)."
  [property f]
  (try (f) false
       (catch clojure.lang.ExceptionInfo e
         (let [data (ex-data e)]
           (and (= :invalid-input (:type data)) (or (nil? property) (= property (:property data))))))
       (catch Throwable _ false)))

(defn- nan? [x] (Double/isNaN (double x)))

;; ---------------------------------------------------------------- 7.1 native discrete distributions

(def ^:private discrete-constructors
  "The keys with the native discrete distributions of `fastmath.random.distributions` with the constructor parameters."
  {:negative-binomial {:r 5.0 :p 0.4}
   :nbi {:mu 4.0 :sigma 0.5}
   :nbii {:mu 4.0 :sigma 0.5}
   :logarithmic {:p 0.5}
   :beta-binomial {:alpha 2.0 :beta 3.0 :n 10}
   :bb nil
   :zibb nil
   :zabb nil
   :zibi nil
   :zabi nil
   :zinbi nil
   :zanbi nil
   :zip nil
   :zip2 nil
   :zap nil
   :integer-discrete {:data [1 2 4] :probabilities [0.2 0.3 0.5]}
   :fishers-noncentral-hypergeometric nil
   :wallenius-noncentral-hypergeometric nil})

(defn- discrete [k] (r/distribution k (get discrete-constructors k)))

(deftest discrete-non-integer-arguments
  (doseq [k (keys discrete-constructors)
          :let [d (discrete k)]]
    (testing k
      (testing "pmf is 0, lpdf -Inf at non-integer arguments"
        (doseq [x [2.5 -0.5 0.5 1e-300 -1e-300 1.0000000000000002 3.9999999999999996]]
          (is (m/zero? (prot/pdf d x)) (str "pdf " x))
          (is (m/zero? (prot/probability d x)) (str "probability " x))
          (is (m/neg-inf? (prot/lpdf d x)) (str "lpdf " x))))
      (testing "cdf floors and is 0 below the support"
        (doseq [x [0.5 1.5 2.5 3.999 10.75 -0.5 -1e-300 -3.2]]
          (is (bits= (prot/cdf d (m/floor x)) (prot/cdf d x)) (str "cdf " x)))
        (is (m/zero? (prot/cdf d -0.5)))
        (is (m/zero? (prot/cdf d -1)))
        (is (m/zero? (prot/cdf d -1e300))))
      (testing "range cdf is the probability of (v1, v2]"
        (is (m/delta-eq (- (prot/cdf d 4) (prot/cdf d 1)) (prot/cdf d 1.5 4.9) 1.0e-15))))))

(deftest discrete-boundary-arguments
  (doseq [k (keys discrete-constructors)
          :let [d (discrete k)]]
    (testing k
      (testing "zero, one, minus one, -0.0"
        (is (bits= (prot/pdf d 0) (prot/pdf d -0.0)))
        (is (bits= (prot/lpdf d 0) (prot/lpdf d -0.0)))
        (is (bits= (prot/cdf d 0) (prot/cdf d -0.0)))
        (is (bits= (prot/pdf d 1) (prot/pdf d 1.0)))
        (is (m/zero? (prot/pdf d -1)))
        (is (m/neg-inf? (prot/lpdf d -1))))
      (testing "large arguments"
        (doseq [x [1.0e9 2147483648.0 4.0e18]]
          (is (m/not-neg? (prot/pdf d x)) (str "pdf " x))
          (is (<= (prot/pdf d x) 1.0) (str "pdf " x))
          (is (m/delta-eq 1.0 (prot/cdf d x) 1.0e-12) (str "cdf " x))))
      (testing "arguments out of the long range and infinities"
        (doseq [x [9.3e18 -9.3e18 1e300 -1e300 9.223372036854775807E18 ##Inf ##-Inf]]
          (is (m/zero? (prot/pdf d x)) (str "pdf " x))
          (is (m/neg-inf? (prot/lpdf d x)) (str "lpdf " x))
          (is (m/delta-eq (if (pos? x) 1.0 0.0) (prot/cdf d x) 1.0e-15) (str "cdf " x))))
      (testing "NaN"
        (is (nan? (prot/pdf d ##NaN)))
        (is (nan? (prot/lpdf d ##NaN)))
        (is (nan? (prot/probability d ##NaN)))
        (is (nan? (prot/cdf d ##NaN)))))))

(deftest discrete-integer-arguments-equal-baseline
  (testing "pmf and cdf at the integers are the values before the audit (bitwise)"
    (let [values (:dist-values @baseline)
          check (fn [baseline-key k ks pdf-key]
                  (let [d (discrete k)
                        v (get values baseline-key)]
                    (is (bits= (get v pdf-key) (mapv #(prot/pdf d %) ks)) (str baseline-key " pdf"))
                    (is (bits= (:cdf v) (mapv #(prot/cdf d %) ks)) (str baseline-key " cdf"))))]
      (check :neg-binomial :negative-binomial [0 1 5 20] :pdf)
      (check :nbi :nbi [0 1 5 20] :pdf)
      (check :beta-binomial :beta-binomial [0 1 5 10] :pdf)
      (check :integer-discrete :integer-discrete [0 1 2 3 4 5] :pdf)
      (let [d (r/zero-inflated-poisson)
            v (:zip values)]
        (is (bits= (:pdf v) (mapv #(prot/pdf d %) (:ks v))))
        (is (bits= (:cdf v) (mapv #(prot/cdf d %) (:ks v)))))
      (let [d (r/logarithmic {:p 0.5})
            v (:logarithmic values)]
        (is (bits= (:pdf v) (mapv #(prot/pdf d %) (:ks v))))))))

(deftest discrete-seeded-streams-equal-baseline
  (let [samples (:dist-samples @baseline)
        s (fn [d n] (mapv long (r/->seq d n)))]
    (is (= (:logarithmic-5 samples) (s (r/logarithmic {:p 0.5 :rng (r/rng :mersenne 9)}) 20)))
    (is (= (:zip-default samples) (s (r/zero-inflated-poisson {:rng (r/rng :mersenne 9)}) 20)))
    (is (= (:nbi samples) (s (r/nbi {:mu 4.0 :sigma 0.5 :rng (r/rng :mersenne 9)}) 20)))
    (is (= (:nbii samples) (s (r/nbii {:mu 4.0 :sigma 0.5 :rng (r/rng :mersenne 9)}) 20)))))

;; ---------------------------------------------------------------- 7.2 logarithmic

(deftest logarithmic-cdf-icdf-against-mpmath
  (doseq [{:keys [p ks cdf icdf-p icdf]} (:logarithmic @reference)
          :let [d (r/logarithmic {:p p})]]
    (testing (str "p=" p)
      (let [errs (mapv (fn [k c] (m/abs (- (prot/cdf d k) c))) ks cdf)]
        ;; observed: below 2e-15 for every k, the tolerance is the one of the specification
        (is (every? #(<= % 1.0e-12) errs) (str "cdf, max abs error " (reduce max errs))))
      (let [bad (remove nil? (map (fn [q k] (when (not= k (prot/icdf d q)) [q k (prot/icdf d q)])) icdf-p icdf))]
        (if (<= p 0.99)
          (is (empty? bad) (str "icdf mismatches " (count bad) " of " (count icdf-p) ": " (vec (take 3 bad))))
          ;; p = .999 is outside of the requirement: a quantile level within 1e-16 of a cdf step (observed: one
          ;; of 756, 1 - 1e-12 at k = 22530) can fall on the other side of the double precision cdf
          (is (and (<= (count bad) 2) (every? (fn [[_ k k']] (<= (m/abs (- k k')) 1)) bad)) (str "icdf mismatches " (vec bad))))))))

(deftest logarithmic-icdf-is-the-smallest-x-with-cdf-geq-p
  (doseq [p [0.3 0.9 0.99]
          :let [d (r/logarithmic {:p p})]]
    (testing (str "p=" p)
      ;; only the steps of the cdf which are visible in double precision (below 1.0 and above the previous value)
      (doseq [k (range 2 400)
              :let [c (prot/cdf d k)
                    c- (prot/cdf d (dec k))]
              :when (and (< c 1.0) (> c c-))]
        (is (= k (prot/icdf d c)) (str "icdf(cdf(" k "))"))
        (is (= k (prot/icdf d (m/next-double c-))) (str "icdf(next(cdf(" (dec k) ")))"))))))

(deftest logarithmic-limits
  (let [d (r/logarithmic {:p 0.99})]
    (testing "pmf outside of the support"
      (doseq [x [0 -1 -3 0.5 1.5 -0.5]]
        (is (m/zero? (prot/pdf d x)) (str x))
        (is (m/zero? (prot/probability d x)) (str x))
        (is (m/neg-inf? (prot/lpdf d x)) (str x))))
    (testing "cdf at the limits returns"
      (is (m/zero? (prot/cdf d ##-Inf)))
      (is (= 1.0 (prot/cdf d ##Inf)))
      (is (= 1.0 (prot/cdf d 1e300)))
      (is (= 1.0 (prot/cdf d 1e9)))
      (is (m/zero? (prot/cdf d 0)))
      (is (m/zero? (prot/cdf d -1)))
      (is (m/zero? (prot/cdf d -1e300))))
    (testing "lpdf in the far tail where the pmf underflows"
      (let [k 200000
            lp (prot/lpdf d k)]
        (is (m/zero? (prot/pdf d 1e8)))
        (is (m/delta-eq (- (+ (m/log (/ -1.0 (m/log 0.01))) (* 1e8 (m/log 0.99))) (m/log 1e8)) (prot/lpdf d 1e8) 1.0e-6))
        (is (m/delta-eq (- (+ (m/log (/ -1.0 (m/log 0.01))) (* k (m/log 0.99))) (m/log k)) lp 1.0e-9))))
    (testing "upper bound unchanged, icdf at the ends"
      (is (== Integer/MAX_VALUE (prot/upper-bound d)))
      (is (== 1 (prot/lower-bound d)))
      (is (= 1 (prot/icdf d 0.0)))
      (is (= ##Inf (prot/icdf d 1.0))))
    (testing "invalid p"
      (doseq [p [0.0 1.0 -0.5 1.5 ##NaN]]
        (is (thrown-input? :p #(r/logarithmic {:p p})) (str p))))))

(deftest logarithmic-sampler-speed
  (doseq [p [0.99 0.999]]
    (let [d (r/logarithmic {:p p :rng (r/rng :mersenne 1)})
          n 20000
          start (System/nanoTime)
          xs (vec (r/->seq d n))
          seconds (/ (- (System/nanoTime) start) 1.0e9)]
      ;; observed: ~2e6 samples/s for p=.99, the requirement is 1e4 samples/s
      (is (< seconds (/ n 1.0e4)) (str "p=" p ", " seconds " s for " n " samples"))
      (is (every? #(>= (long %) 1) xs)))))

(deftest logarithmic-sampler-distribution
  (let [p 0.9
        d (r/logarithmic {:p p :rng (r/rng :mersenne 7)})
        n 200000
        xs (vec (r/->seq d n))
        mean-v (/ (reduce + xs) (double n))
        sd (m/sqrt (r/variance d))]
    (is (< (m/abs (- mean-v (r/mean d))) (* 4.0 (/ sd (m/sqrt n)))) "mean within 4 standard errors")
    (let [p1 (prot/pdf d 1)]
      (is (< (m/abs (- (/ (count (filter #(= 1 %) xs)) (double n)) p1)) (* 4.0 (m/sqrt (/ (* p1 (- 1.0 p1)) n)))) "P(X=1)"))))

(deftest logarithmic-baseline
  (let [d (r/logarithmic {:p 0.5})
        v (get-in @baseline [:dist-values :logarithmic])]
    (testing "pdf is bitwise equal, cdf equals within the rounding of the new series (observed difference 2.1e-15)"
      (is (bits= (:pdf v) (mapv #(prot/pdf d %) (:ks v))))
      (is (every? true? (map #(m/delta-eq %1 (prot/cdf d %2) 1.0e-13) (:cdf v) (:ks v)))))
    (testing "icdf at the baseline points"
      (is (= (:icdf v) (mapv #(prot/icdf d %) [0.1 0.5 0.9 0.99]))))))

;; ---------------------------------------------------------------- 7.4 zero-inflated-poisson

(defn- smallest-with-cdf [d p]
  (first (filter #(>= (prot/cdf d %) p) (range 0 1000))))

(deftest zero-inflated-poisson-icdf
  (doseq [[k d] [[:zip (r/distribution :zip)]
                 [:zip2 (r/distribution :zip2)]
                 [:zip-sigma (r/distribution :zip {:mu 3.5 :sigma 0.25})]
                 [:zip-sigma-high (r/distribution :zip {:mu 0.5 :sigma 0.9})]]]
    (testing k
      (let [cdfs (mapv #(prot/cdf d %) (range 31))
            ;; cdf = 1.0 (the upper end of the support) is not a quantile level of the search
            ps (remove #(>= % 1.0)
                       (concat (map #(/ % 20.0) (range 1 20))
                               (map #(/ % 60.0) (range 1 60))
                               cdfs
                               (map #(- % 1.0e-12) cdfs)
                               [1.0e-9 0.999999 (- 1.0 1.0e-12)]))
            bad (remove #(= (prot/icdf d %) (smallest-with-cdf d %)) ps)]
        (is (empty? bad) (str "mismatches: " (vec (take 5 bad))))
        (is (every? true? (keep-indexed (fn [k c] (when (and (< c 1.0) (or (zero? k) (> c (nth cdfs (dec k)))))
                                                    (= k (prot/icdf d c))))
                                        cdfs))
            "icdf(cdf(k)) = k where the cdf has a step")))))

(deftest zero-inflated-poisson-baseline
  (let [d (r/zero-inflated-poisson)
        v (get-in @baseline [:dist-values :zip])]
    (is (= (:icdf v) (mapv #(prot/icdf d %) (map #(/ % 20.0) (range 1 20)))))))

;; ---------------------------------------------------------------- 7.5 dirichlet

(deftest dirichlet-scipy
  (doseq [{:keys [alpha points pdf lpdf]} (:dirichlet @reference)
          :let [d (r/dirichlet {:alpha alpha})]]
    (testing alpha
      (doseq [[x p lp] (map vector points pdf lpdf)]
        (is (rel-eq? (prot/pdf d x) p 1.0e-12) (str "pdf " x))
        (is (m/delta-eq (prot/lpdf d x) lp 1.0e-12) (str "lpdf " x))))))

(deftest dirichlet-boundary
  (testing "alpha = 1: finite limit"
    (let [d (r/dirichlet {:alpha [1 1]})]
      (is (m/delta-eq 1.0 (prot/pdf d [1.0 0.0])))
      (is (m/delta-eq 1.0 (prot/pdf d [0.0 1.0])))
      (is (m/delta-eq 0.0 (prot/lpdf d [1.0 0.0])))
      (is (m/delta-eq 1.0 (prot/pdf d [0.3 0.7])))))
  (testing "alpha < 1: +Inf at the boundary, the limit grows without bound"
    (let [d (r/dirichlet {:alpha [0.5 0.5]})]
      (is (= ##Inf (prot/pdf d [1.0 0.0])))
      (is (= ##Inf (prot/lpdf d [0.0 1.0])))
      (is (>= (prot/pdf d [(- 1.0 1.0e-12) 1.0e-12]) 3.0e5))
      ;; true value 1/(pi sqrt(x (1-x))) = 318309.886...
      (is (rel-eq? (prot/pdf d [(- 1.0 1.0e-12) 1.0e-12]) 318309.8861839 1.0e-8))))
  (testing "alpha > 1: zero at the boundary"
    (let [d (r/dirichlet {:alpha [2 2]})]
      (is (m/zero? (prot/pdf d [1.0 0.0])))
      (is (m/neg-inf? (prot/lpdf d [1.0 0.0])))))
  (testing "mixed alphas: the zero wins"
    (let [d (r/dirichlet {:alpha [0.5 2.0]})]
      (is (= ##Inf (prot/pdf d [0.0 1.0])))
      (is (m/zero? (prot/pdf d [1.0 0.0])))))
  (testing "three dimensions, one zero coordinate"
    (let [d (r/dirichlet {:alpha [1.0 1.0 1.0]})]
      (is (m/delta-eq 2.0 (prot/pdf d [0.5 0.5 0.0])))
      (is (m/delta-eq 2.0 (prot/pdf d [1.0 0.0 0.0])))))
  (testing "NaN"
    (is (nan? (prot/pdf (r/dirichlet {:alpha [2 3]}) [##NaN 0.5])))))

(deftest dirichlet-invalid-points
  (let [d (r/dirichlet {:alpha [2 3]})]
    (doseq [[label x] {"sum above 1" [0.9 0.9]
                       "sum below 1" [0.3 0.3]
                       "negative entry" [-0.1 1.1]
                       "too short" [0.5]
                       "too long" [0.2 0.3 0.5]
                       "empty" []
                       "nil" nil
                       "number" 0.5}]
      (is (thrown-input? :x #(prot/pdf d x)) (str "pdf " label))
      (is (thrown-input? :x #(prot/lpdf d x)) (str "lpdf " label))
      (is (thrown-input? :x #(prot/probability d x)) (str "probability " label)))
    (testing "sum within 1e-8 is accepted"
      (is (m/pos? (prot/pdf d [0.4 (+ 0.6 5.0e-9)]))))))

(deftest dirichlet-sampler
  (testing "unchanged stream for the ordinary concentration (bitwise equal to the baseline)"
    (is (= (:dirichlet-234 (:dist-samples @baseline))
           (vec (r/->seq (r/dirichlet {:alpha [2.0 3.0 4.0] :rng (r/rng :mersenne 9)}) 10)))))
  (testing "tiny concentration: the mean of x1 is alpha1 / sum(alpha) = 0.25 (the unit vector fallback gave 0.484)"
    (doseq [seed [11 12 13]]
      (let [d (r/dirichlet {:alpha [0.001 0.003] :rng (r/rng :mersenne seed)})
            n 100000
            xs (mapv first (r/->seq d n))
            mean-v (/ (reduce + xs) (double n))
            se (/ (m/sqrt (* mean-v (- 1.0 mean-v))) (m/sqrt n))]
        (is (< (m/abs (- mean-v 0.25)) (* 4.0 se)) (str "seed " seed ", mean " mean-v)))))
  (testing "samples sum to 1 and stay in the simplex"
    (doseq [alpha [[2.0 3.0 4.0] [0.001 0.003 0.01] [0.05 0.15] [0.5 0.5 0.5 0.5]]]
      (let [xs (r/->seq (r/dirichlet {:alpha alpha :rng (r/rng :mersenne 3)}) 5000)]
        (is (every? #(m/delta-eq 1.0 (reduce + %) 1.0e-12) xs) (str alpha))
        (is (every? #(every? (fn [x] (<= 0.0 x 1.0)) %) xs) (str alpha))
        (is (every? #(= (count alpha) (count %)) xs))))))

(deftest dirichlet-baseline
  (let [d (r/dirichlet {:alpha [2.0 3.0 4.0]})
        v (get-in @baseline [:dist-values :dirichlet])
        pts [[0.2 0.3 0.5] [0.1 0.1 0.8]]]
    (is (bits= (:pdf v) (mapv #(prot/pdf d %) pts)))
    (is (bits= (:lpdf v) (mapv #(prot/lpdf d %) pts)))))

;; ---------------------------------------------------------------- 7.6 multinomial

(deftest multinomial-scipy
  (doseq [{:keys [n ps points pmf cdf-points cdf]} (:multinomial @reference)
          :let [d (r/multinomial {:n n :ps ps})]]
    (testing [n ps]
      (doseq [[x p] (map vector points pmf)]
        (is (rel-eq? (prot/pdf d x) p 1.0e-12) (str "pdf " x))
        (is (rel-eq? (prot/probability d x) p 1.0e-12) (str "probability " x))
        (is (rel-eq? (m/exp (prot/lpdf d x)) p 1.0e-12) (str "lpdf " x)))
      (doseq [[x c] (map vector cdf-points cdf)]
        (is (rel-eq? (prot/cdf d x) c 1.0e-12) (str "cdf " x))))))

(deftest multinomial-zero-probability-category
  (let [d (r/multinomial {:n 10 :ps [1 0 1]})]
    (is (m/delta-eq 0.24609375 (prot/pdf d [5 0 5]) 1.0e-12))
    (is (m/delta-eq (m/log 0.24609375) (prot/lpdf d [5 0 5]) 1.0e-12))
    (is (m/delta-eq 0.24609375 (prot/cdf d [5 0 5]) 1.0e-12))
    (is (m/delta-eq 1.0 (prot/cdf d [10 10 10]) 1.0e-12))
    (testing "positive count in the zero probability category"
      (is (m/zero? (prot/pdf d [4 1 5])))
      (is (m/neg-inf? (prot/lpdf d [4 1 5]))))
    (testing "samples keep the zero count"
      (is (every? #(and (zero? (second %)) (= 10 (reduce + %))) (r/->seq d 200))))))

(deftest multinomial-invalid-outcomes
  (let [d (r/multinomial {:n 10 :ps [1 2 3]})]
    (testing "negative count: pmf 0, cdf 0"
      (is (m/zero? (prot/pdf d [-1 5 6])))
      (is (m/neg-inf? (prot/lpdf d [-1 5 6])))
      (is (m/zero? (prot/cdf d [-1 5 6]))))
    (testing "non-integer counts and a wrong sum have the probability 0"
      (is (m/zero? (prot/pdf d [2.5 2.5 5])))
      (is (m/zero? (prot/pdf d [1 1 1])))
      (is (m/zero? (prot/pdf d [11 0 0])))
      (is (m/zero? (prot/pdf d [1e300 0 0]))))
    (testing "NaN"
      (is (nan? (prot/pdf d [##NaN 5 5]))))
    (testing "wrong length and nil are errors"
      (is (thrown-input? :x #(prot/pdf d [5 5])))
      (is (thrown-input? :x #(prot/pdf d nil)))
      (is (thrown-input? :x #(prot/lpdf d nil)))
      (is (thrown-input? :x #(prot/cdf d [1 2 3 4]))))
    (testing "lpdf where the pmf underflows"
      (let [d (r/multinomial {:n 1000 :ps [1 1]})]
        (is (m/zero? (prot/pdf (r/multinomial {:n 5000 :ps [1 1]}) [0 5000])))
        ;; log(0.5^1000) = -1000 log 2, the multinomial coefficient is 1
        (is (m/delta-eq (* 1000.0 (- (m/log 2.0))) (prot/lpdf d [0 1000]) 1.0e-9))))))

(deftest multinomial-zero-trials
  (let [d (r/multinomial {:n 0 :ps [1 1]})]
    (is (= 1.0 (prot/pdf d [0 0])))
    (is (= 0.0 (prot/lpdf d [0 0])))
    (is (m/zero? (prot/pdf d [1 0])))
    (is (= 1.0 (prot/cdf d [0 0])))
    (is (= [0 0] (vec (r/sample d))))
    (is (= [0.0 0.0] (r/means d))))
  (testing "n has to be an integer >= 0"
    (doseq [n [-1 2.5 ##NaN ##Inf 3.0e10 nil "5"]]
      (is (thrown-input? :n #(distr/multinomial n [1 1] r/binomial nil)) (str n)))))

(deftest multinomial-invalid-probabilities
  (doseq [ps [[] nil [1 -1] [0 0] [##NaN 1] [##Inf 1] 5 [1 :a]]]
    (is (thrown-input? :ps #(distr/multinomial 5 ps r/binomial nil)) (str ps))))

(deftest multinomial-baseline
  (let [d (r/multinomial {:n 10 :ps [1 2 3]})
        v (get-in @baseline [:dist-values :multinomial])]
    (is (bits= (:pdf v) (mapv #(prot/pdf d %) [[1 3 6] [2 2 6] [0 5 5]])))
    (is (bits= (:cdf v) (mapv #(prot/cdf d %) [[1 3 6] [5 5 5]])))))

;; ---------------------------------------------------------------- 7.7, 5.5 multi-normal

(deftest multi-normal-scipy
  (doseq [{:keys [means covariances points pdf lpdf]} (get-in @reference [:multi-normal :interior])
          :let [d (distr/multi-normal means covariances nil)]]
    (testing (str (count means) " dimensions")
      (doseq [[x p lp] (map vector points pdf lpdf)]
        (is (rel-eq? (prot/pdf d x) p 1.0e-12) (str "pdf " x))
        (is (rel-eq? (prot/lpdf d x) lp 1.0e-12) (str "lpdf " x)))))
  (testing "the log density where the density underflows"
    (doseq [{:keys [means covariances point lpdf]} (get-in @reference [:multi-normal :far])
            :let [d (distr/multi-normal means covariances nil)]]
      (is (m/zero? (prot/pdf d point)) "the density underflows")
      (is (rel-eq? (prot/lpdf d point) lpdf 1.0e-12) (str point)))
    (let [d (distr/multi-normal [0.0 0.0] [[1.0 0.5] [0.5 2.0]] nil)]
      (is (rel-eq? (prot/lpdf d [40.0 -40.0]) -1830.6891135318056 1.0e-12)))))

(deftest multi-normal-baseline
  (let [d (r/multi-normal {:means [0.0 1.0] :covariances [[1.0 0.5] [0.5 2.0]]})
        v (get-in @baseline [:dist-values :multi-normal])]
    (is (bits= (:pdf v) (mapv #(prot/pdf d %) (:points v))))
    (is (bits= (:lpdf v) (mapv #(prot/lpdf d %) (:points v))))))

(deftest multi-normal-validation
  (testing "invalid covariances"
    (doseq [[label cov] {"asymmetric" [[1 0.9] [0.1 1]]
                         "ragged" [[1 0] [0]]
                         "not square" [[1 0 0] [0 1 0]]
                         "not positive definite" [[1 2] [2 1]]
                         "singular" [[1 1] [1 1]]
                         "zero matrix" [[0 0] [0 0]]
                         "negative variance" [[-1 0] [0 1]]
                         "NaN entry" [[1 ##NaN] [##NaN 1]]
                         "infinite entry" [[##Inf 0] [0 1]]
                         "empty" []
                         "not a matrix" [1 2 3]
                         "number" 5}]
      (is (thrown-input? :covariances #(distr/multi-normal nil cov nil)) label)))
  (testing "invalid means"
    (doseq [[label means] {"not finite" [##NaN 0]
                           "empty" []
                           "number" 3}]
      (is (thrown-input? :means #(distr/multi-normal means [[1 0] [0 1]] nil)) label))
    (is (thrown-input? :means #(distr/multi-normal [0 0 0] [[1 0] [0 1]] nil)) "sizes do not match"))
  (testing "valid inputs and defaults"
    (let [d (distr/multi-normal nil [[1 0 0] [0 1 0] [0 0 1]] nil)]
      (is (= [0.0 0.0 0.0] (r/means d)) "zero means of the dimension of covariances"))
    (let [d (distr/multi-normal [1 2 3] nil nil)]
      (is (= [[1.0 0.0 0.0] [0.0 1.0 0.0] [0.0 0.0 1.0]] (mapv vec (r/covariance d))) "identity covariance"))
    (let [d (distr/multi-normal nil nil nil)]
      (is (= 2 (r/dimensions d))))
    (is (some? (distr/multi-normal [0 0] [[1 0.5] [(+ 0.5 1.0e-11) 2]] nil)) "symmetric within the tolerance")
    (is (some? (distr/multi-normal [0] [[2.0]] nil)) "one dimension")
    (is (= 3 (r/dimensions (distr/multi-normal [0 0 0] (mapv vec (partition 3 [2 0.5 0 0.5 2 0 0 0 2])) (r/rng :mersenne 1))))))
  (testing "the point of pdf and lpdf"
    (let [d (distr/multi-normal [0 0] [[1 0.5] [0.5 2]] nil)]
      (doseq [x [[1] [1 2 3] nil [] 3.0]]
        (is (thrown-input? :x #(prot/pdf d x)) (str "pdf " x))
        (is (thrown-input? :x #(prot/lpdf d x)) (str "lpdf " x)))
      (is (nan? (prot/pdf d [##NaN 0.0]))))))

;; ---------------------------------------------------------------- 7.9 finite support distributions

(def ^:private invalid-discrete-inputs
  "[label data probabilities] of the inputs rejected by `integer-discrete`, `real-discrete` and `categorical`."
  [["empty data" [] nil]
   ["nil data" nil nil]
   ["empty data and probabilities" [] []]
   ["data is a number" 5 nil]
   ["data is a map" {1 2} nil]
   ["fewer probabilities" [1 2 3] [0.5 0.5]]
   ["more probabilities" [1 2 3] [0.2 0.3 0.5 0.5]]
   ["negative probability" [1 2 3] [0.6 -0.1 0.5]]
   ["NaN probability" [1 2 3] [0.5 ##NaN 0.5]]
   ["infinite probability" [1 2 3] [0.5 ##Inf 0.5]]
   ["negative infinite probability" [1 2 3] [0.5 ##-Inf 0.5]]
   ["zero sum" [1 2 3] [0 0 0]]
   ["non-numeric probability" [1 2 3] [0.5 :a 0.5]]
   ["probabilities is a number" [1 2 3] 0.5]])

(deftest finite-support-invalid-inputs
  (doseq [[label data probabilities] invalid-discrete-inputs]
    (is (thrown-input? nil #(r/integer-discrete {:data data :probabilities probabilities})) (str "integer-discrete: " label))
    (is (thrown-input? nil #(r/real-discrete {:data data :probabilities probabilities})) (str "real-discrete: " label))
    (is (thrown-input? nil #(r/categorical {:data data :probabilities probabilities})) (str "categorical: " label)))
  (testing "the properties"
    (is (thrown-input? :data #(r/integer-discrete {:data []})))
    (is (thrown-input? :probabilities #(r/integer-discrete {:data [1 2 3] :probabilities [0.5 0.5]})))
    (is (thrown-input? :probabilities #(r/real-discrete {:data [1.0 2.0] :probabilities [0.5 -0.5]})))
    (is (thrown-input? :probabilities #(r/categorical {:data [:a :b] :probabilities [0 0]}))))
  (testing "data of the integer distribution have to be integers"
    (doseq [data [[1 2.5 3] [1.0 ##NaN] [1 ##Inf] [1 :a] [1 "2"] [1 9.3e18]]]
      (is (thrown-input? :data #(r/integer-discrete {:data data})) (str data))))
  (testing "data of the real distribution cannot be NaN"
    (is (thrown-input? :data #(r/real-discrete {:data [1.0 ##NaN] :probabilities [0.5 0.5]})))
    (is (thrown-input? :data #(r/real-discrete {:data [1.0 :a]})))))

(deftest finite-support-valid-inputs
  (testing "integer-valued doubles and longs outside the int range are accepted"
    (let [d (r/integer-discrete {:data [1.0 2.0 4.0] :probabilities [0.2 0.3 0.5]})]
      (is (m/delta-eq 0.3 (prot/pdf d 2) 1.0e-15)))
    (let [d (r/integer-discrete {:data [1 2 1099511627776] :probabilities [0.2 0.3 0.5]})]
      (is (m/delta-eq 0.5 (prot/pdf d 1099511627776) 1.0e-15))))
  (testing "missing probabilities are equal"
    (is (m/delta-eq 0.25 (prot/pdf (r/integer-discrete {:data [1 2 3 4]}) 3) 1.0e-15))
    (is (m/delta-eq 0.25 (prot/pdf (r/real-discrete {:data [1.5 2 3 4]}) 1.5) 1.0e-15))
    (is (m/delta-eq 0.25 (prot/pdf (r/categorical {:data [:a :b :c :d]}) :c) 1.0e-15)))
  (testing "zero probabilities are accepted"
    (is (m/zero? (prot/pdf (r/integer-discrete {:data [1 2 3] :probabilities [0.5 0 0.5]}) 2)))
    (is (m/zero? (prot/pdf (r/categorical {:data [:a :b :c] :probabilities [1 0 1]}) :b)))))

(deftest finite-support-baseline
  (let [d (r/integer-discrete {:data [1 2 4] :probabilities [0.2 0.3 0.5]})
        v (get-in @baseline [:dist-values :integer-discrete])]
    (is (bits= (:pdf v) (mapv #(prot/pdf d %) [0 1 2 3 4 5])))
    (is (bits= (:cdf v) (mapv #(prot/cdf d %) [0 1 2 3 4 5])))))

;; ---------------------------------------------------------------- 5.2 inverse-gaussian

(deftest inverse-gaussian-cdf-against-mpmath
  (doseq [{:keys [mu lambda xs cdf]} (:inverse-gaussian @reference)
          :let [d (distr/inverse-gaussian mu lambda nil)]]
    (testing (str "mu=" mu " lambda=" lambda " (lambda/mu = " (/ lambda mu) ")")
      (doseq [[x c] (map vector xs cdf)]
        (is (rel-eq? (prot/cdf d x) c 1.0e-9) (str "cdf(" x ")"))))))

(deftest inverse-gaussian-icdf-round-trip
  (doseq [[mu lambda] [[1.0 1.0] [1.0 100.0] [1.0 354.0] [1.0 355.0] [1.0 400.0] [1.0 1000.0] [1.0 10000.0]
                       [0.01 100.0] [5.0 1775.0] [3.0 3.0e5]]
          :let [d (distr/inverse-gaussian mu lambda nil)
                ps [1.0e-9 1.0e-6 0.001 0.1 0.5 0.9 0.999 (- 1.0 1.0e-9)]
                xs (mapv #(prot/icdf d %) ps)]]
    (testing (str "mu=" mu " lambda=" lambda)
      (is (every? #(<= 0.0 % ##Inf) xs))
      (is (apply < xs) "increasing, not degenerate")
      (is (every? true? (map #(<= (m/abs (- %1 (prot/cdf d %2))) 1.0e-9) ps xs))
          (str "|cdf(icdf p) - p| " (apply max (map #(m/abs (- %1 (prot/cdf d %2))) ps xs))))))
  (testing "ends"
    (let [d (distr/inverse-gaussian 1.0 1000.0 nil)]
      (is (= 0.0 (prot/icdf d 0.0)))
      (is (= ##Inf (prot/icdf d 1.0)))
      (is (= 0.0 (prot/cdf d 0.0)))
      (is (= 0.0 (prot/cdf d -1.0)))
      (is (= 1.0 (prot/cdf d ##Inf)))
      (is (= 1.0 (prot/cdf d 1e300)))
      (is (= 0.0 (prot/cdf d 1e-300)))
      (is (m/delta-eq 0.5 (prot/cdf d 1.0) 0.02) "median is close to the mean for a concentrated distribution")
      (is (m/delta-eq (- (prot/cdf d 1.1) (prot/cdf d 0.9)) (prot/cdf d 0.9 1.1) 1.0e-15)))))

(defn- ks-p-value
  "Two-sided one-sample Kolmogorov-Smirnov p-value of `xs` against `cdf` (asymptotic Kolmogorov distribution)."
  [xs cdf]
  (let [xs (vec (sort xs))
        n (count xs)
        dn (reduce max (map-indexed (fn [i x]
                                      (let [c (cdf x)]
                                        (max (- c (/ i (double n))) (- (/ (inc i) (double n)) c))))
                                    xs))
        sn (m/sqrt n)
        z (* (+ sn 0.12 (/ 0.11 sn)) dn)]
    (- 1.0 (prot/cdf (distr/kolmogorov nil) z))))

(deftest inverse-gaussian-sampler
  (doseq [[mu lambda] [[1.0 354.0] [1.0 355.0] [1.0 1000.0] [1.0 10000.0] [0.01 100.0]]
          seed [1 2 3]]
    (let [d (distr/inverse-gaussian mu lambda (r/rng :mersenne seed))
          xs (vec (r/->seq d 200000))]
      (testing (str "mu=" mu " lambda=" lambda " seed=" seed)
        (is (> (count (distinct xs)) 199000) "not degenerate")
        (is (> (ks-p-value xs #(prot/cdf d %)) 0.01) "Kolmogorov-Smirnov")
        (is (< (m/abs (- (/ (reduce + xs) 200000.0) mu)) (* 4.0 (/ (m/sqrt (r/variance d)) (m/sqrt 200000.0)))) "mean")))))

(deftest inverse-gaussian-small-ratio-equals-ssj
  (testing "lambda/mu <= 100: SSJ values, bitwise equal seeded samples"
    (let [samples (:dist-samples @baseline)]
      (is (= (:inverse-gaussian-1-50 samples)
             (vec (r/->seq (distr/inverse-gaussian 1.0 50.0 (r/rng :mersenne 9)) 10))))
      (is (= (:inverse-gaussian-2-100 samples)
             (vec (r/->seq (distr/inverse-gaussian 2.0 100.0 (r/rng :mersenne 9)) 10)))))
    (let [d (distr/inverse-gaussian 1.0 50.0 nil)
          v (get-in @baseline [:dist-values :inverse-gaussian-1-50])]
      (is (bits= (:pdf v) (mapv #(prot/pdf d %) (:xs v))))
      (is (bits= (:cdf v) (mapv #(prot/cdf d %) (:xs v))))))
  (testing "the protocols of the large ratio object"
    (let [d (distr/inverse-gaussian 1.0 1000.0 (r/rng :mersenne 1))]
      (is (= :inverse-gaussian (prot/distribution-id d)))
      (is (= [:mu :lambda :rng] (prot/distribution-parameters d)))
      (is (prot/continuous? d))
      (is (= 1 (prot/dimensions d)))
      (is (m/delta-eq 1.0 (prot/mean d) 1.0e-12))
      (is (m/delta-eq 0.001 (prot/variance d) 1.0e-12))
      (is (= 0.0 (prot/lower-bound d)))
      (is (= ##Inf (prot/upper-bound d)))
      (let [x1 (do (prot/set-seed! d 5) (prot/sample d))
            x2 (do (prot/set-seed! d 5) (prot/sample d))]
        (is (= x1 x2) "reseeding repeats the sample"))
      (is (pos? (prot/pdf d 1.0)))
      (is (m/delta-eq (prot/pdf d 1.0) (prot/probability d 1.0) 0.0)))))

;; ---------------------------------------------------------------- 5.6 distribution-id / distribution-parameters

(def ^:private constructor-aliases
  "Keys of `distribution` whose constructor in `fastmath.random` has a different name."
  {:zibi 'zero-inflated-binomial :zabi 'zero-adjusted-binomial :zibb 'zero-inflated-beta-binomial
   :zabb 'zero-adjusted-beta-binomial :zinbi 'zero-inflated-negative-binomial :zanbi 'zero-adjusted-negative-binomial
   :zip 'zero-inflated-poisson :zip2 'zero-inflated-poisson2 :zap 'zero-adjusted-poisson
   :zaga 'zero-adjusted-gamma :zaig 'zero-adjusted-inverse-gaussian :gev 'generalized-extreme-value
   :gpd 'generalized-pareto :ge 'generalized-exponential :gg 'generalized-gamma :gnd 'generalized-normal
   :gig 'generalized-inverse-gaussian :gh 'generalized-hyperbolic :ghl 'generalized-half-logistic})

(defn- constructor-var [k]
  (ns-resolve 'fastmath.random (or (constructor-aliases k) (symbol (name k)))))

(defn- constructor-keys
  "The keys destructured (`:keys`) by the constructor from its `:arglists`, `nil` when there is no destructuring."
  [v]
  (let [ks (for [arglist (:arglists (meta v))
                 form arglist
                 :when (and (map? form) (:keys form))]
             (set (map keyword (:keys form))))]
    (when (seq ks) (apply set/union ks))))

(def ^:private documented-ids
  "Keys whose `distribution-id` is another (canonical) key: aliases of the constructors with another name, and keys sharing the
  object of another key (`distribution-id` of an object does not know the key it was created with)."
  {:bb #{:beta-binomial}
   :nbi #{:negative-binomial}
   :nbii #{:negative-binomial}
   :exgaus #{:ex-gaussian}
   :bernoulli #{:binomial}
   :zip2 #{:zero-inflated-poisson}
   :zero-inflated-poisson2 #{:zero-inflated-poisson}
   :categorical #{:categorical-distribution}
   :integer-discrete #{:integer-discrete-distribution}
   :real-discrete #{:real-discrete-distribution}
   :kde #{:continuous-distribution}})

(def ^:private known-mismatches
  "Keys whose `distribution-parameters` differs from the keys destructured by the constructor in `fastmath.random`, the value is
  the reason. Remove the key when the constructor is fixed (this test fails for a fixed key to keep the map up to date).
  `bb`, `nbi`, `nbii` and `exgaus` were fixed by passing the parameter list of the alias to the shared constructor (WP 5.6)."
  {:bernoulli "the Commons-derived object reports the binomial parameters [:trials :p :rng], the constructor takes [:p :rng] (documented, the id is :binomial)"})

(deftest distribution-id-and-parameters
  (let [ks (sort (keys (methods r/distribution)))]
    (is (> (count ks) 120))
    (doseq [k ks
            :let [v (constructor-var k)
                  d (r/distribution k)
                  id (r/distribution-id d)
                  params (set (r/distribution-parameters d))
                  ctor-keys (some-> v constructor-keys)]]
      (testing k
        (is (some? v) "the constructor is found")
        (is (contains? (conj (get documented-ids k #{}) k (keyword (name (or (constructor-aliases k) (symbol (name k)))))) id)
            (str "id " id))
        (is (= (count (r/distribution-parameters d)) (count params)) "no repeated parameters")
        (cond
          (nil? ctor-keys) (is (#{:continuous-distribution :kde} k) "constructors without destructuring")
          (contains? known-mismatches k) (is (not= ctor-keys params) (str "fixed? remove " k " from known-mismatches"))
          :else (is (= ctor-keys params) (str "parameters " (sort params) " constructor " (sort ctor-keys))))))))

(deftest fixed-ids-and-parameters
  (is (= :half-cauchy (r/distribution-id (r/distribution :half-cauchy))))
  (doseq [k [:fishers-noncentral-hypergeometric :wallenius-noncentral-hypergeometric]]
    (is (= [:ns :nf :n :omega :rng] (r/distribution-parameters (r/distribution k))) (str k))
    (is (= k (r/distribution-id (r/distribution k)))))
  (testing "the optional parameters of the shared constructors"
    (is (= [:mu :sigma :rng] (prot/distribution-parameters (distr/negative-binomial 5.0 0.4 nil [:mu :sigma :rng]))))
    (is (= [:r :p :rng] (prot/distribution-parameters (distr/negative-binomial 5.0 0.4 nil))))
    (is (= [:mu :sigma :bd :rng] (prot/distribution-parameters (distr/beta-binomial 2.0 3.0 10 nil [:mu :sigma :bd :rng]))))
    (is (= [:alpha :beta :n :rng] (prot/distribution-parameters (distr/beta-binomial 2.0 3.0 10 nil))))))

;; ---------------------------------------------------------------- 5.7 undefined moments

(def ^:private undefined-moments
  "[key parameters mean-expectation variance-expectation]: `:nan`, `:inf` (+Infinity) or `:finite`."
  [[:generalized-extreme-value {:mu 0.0 :sigma 1.0 :xi 1.0} :inf :inf]
   [:generalized-extreme-value {:mu 0.0 :sigma 1.0 :xi 1.5} :inf :inf]
   [:generalized-extreme-value {:mu 0.0 :sigma 1.0 :xi 0.5} :finite :inf]
   [:generalized-extreme-value {:mu 0.0 :sigma 1.0 :xi 0.7} :finite :inf]
   [:generalized-extreme-value {:mu 0.0 :sigma 1.0 :xi 0.3} :finite :finite]
   [:generalized-pareto {:mu 0.0 :sigma 1.0 :xi 1.0} :inf :inf]
   [:generalized-pareto {:mu 0.0 :sigma 1.0 :xi 2.0} :inf :inf]
   [:generalized-pareto {:mu 0.0 :sigma 1.0 :xi 0.5} :finite :inf]
   [:generalized-pareto {:mu 0.0 :sigma 1.0 :xi 0.8} :finite :inf]
   [:generalized-gamma {:mu 1.0 :sigma 1.0 :nu -1.5} :inf :inf]
   [:generalized-gamma {:mu 1.0 :sigma 1.0 :nu -0.8} :finite :inf]
   [:f-noncentral {:df1 5.0 :df2 2.0 :ncp 1.0} :inf :inf]
   [:f-noncentral {:df1 5.0 :df2 1.0 :ncp 1.0} :inf :inf]
   [:f-noncentral {:df1 5.0 :df2 3.0 :ncp 1.0} :finite :inf]
   [:f-noncentral {:df1 5.0 :df2 4.0 :ncp 1.0} :finite :inf]
   [:f-noncentral {:df1 5.0 :df2 5.0 :ncp 1.0} :finite :finite]
   [:t-noncentral {:df 1.0 :ncp 1.5} :nan :nan]
   [:t-noncentral {:df 0.5 :ncp 1.5} :nan :nan]
   [:t-noncentral {:df 1.5 :ncp 1.5} :finite :inf]
   [:t-noncentral {:df 2.0 :ncp 1.5} :finite :inf]
   [:t-noncentral {:df 3.0 :ncp 1.5} :finite :finite]
   [:half-cauchy {} :nan :nan]
   [:half-cauchy {:mu 1.0 :scale 3.0} :nan :nan]])

(deftest undefined-moments-of-native-keys
  (doseq [[k params mean-kind variance-kind] undefined-moments
          :let [d (r/distribution k params)
                check (fn [what v kind]
                        (case kind
                          :nan (is (nan? v) (str k " " params " " what " = " v))
                          :inf (is (= ##Inf v) (str k " " params " " what " = " v))
                          :finite (is (and (m/valid-double? v) (or (= what "mean") (m/not-neg? v))) (str k " " params " " what " = " v))))]]
    (testing (str k " " params)
      (let [mean-v (try (prot/mean d) (catch Throwable e e))
            variance-v (try (prot/variance d) (catch Throwable e e))]
        (is (number? mean-v) (str "mean does not throw: " mean-v))
        (is (number? variance-v) (str "variance does not throw: " variance-v))
        (when (number? mean-v) (check "mean" mean-v mean-kind))
        (when (number? variance-v)
          (check "variance" variance-v variance-kind)
          (is (not= ##-Inf variance-v) "never a negative infinite variance"))))))
