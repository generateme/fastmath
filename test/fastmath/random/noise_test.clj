(ns fastmath.random.noise-test
  "Tests of the noise functions of `fastmath.random` (Random Audit, work packages 3.2-3.8, 3.11, 9.4): the presets `noise`, `vnoise`,
  `simplex`, `single-noise`, `fbm-noise`, `billow-noise`, `ridgedmulti-noise`, `warp-noise(-fn)`, `random-noise-cfg`, `random-noise(-fn)`
  and `discrete-noise`. The Java classes alone are covered by `noise-java-test`.

  Oracles
  - `test/resources/random/noise_reference.edn` (`utils/fastmath/dev/generate_random_noise_reference.{clj,py}`): 216 cases x 12 points
    from an independent Python re-implementation of the published algorithms (Gustavson simplex, fade/lerp gradient and value noise,
    libnoise style fbm/billow/ridged multifractal with the NEW ridged normaliser `1 / sum lacunarity^-k`). Only the lookup tables of
    the seed (`perm`, `perm12`, `valueLUT`) come from the library.
  - the seeded outputs of the code before the audit in `test/resources/random/baseline_reference.edn` (`:noise`, 96 cases x 12 points):
    bitwise equal, except ridged noise with lacunarity 2.1 / gain 0.45 where new = old * (sum gain^k / sum lacunarity^-k).
  - structural facts of the noise definitions (lattice zeros of gradient noise, continuity classes of the interpolations, ranges),
    the preset goldens pin the library output (properties independent of the golden: preset = `fbm-noise` with the seed 1337)."
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.string :as str]
            [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m]
            [fastmath.random :as r])
  (:import [java.util Arrays SplittableRandom]))

(set! *warn-on-reflection* true)

(def ^:private reference
  (delay (edn/read-string (slurp (io/resource "random/noise_reference.edn")))))

(def ^:private baseline
  (delay (:noise (edn/read-string (slurp (io/resource "random/baseline_reference.edn"))))))

(def ^:private noise-points
  "The points of the baseline by dimension (copy of `noise-points` in generate_random_baseline.clj)."
  {1 [[0.3] [-4.7] [123.456] [-0.0001]]
   2 [[0.3 0.4] [-4.7 2.2] [123.456 -78.9] [10.0 10.0]]
   3 [[0.3 0.4 0.5] [-4.7 2.2 9.9] [123.4 -78.9 1.1] [0.5 0.5 0.5]]})

(def ^:private generator-fns
  {:single r/single-noise
   :fbm r/fbm-noise
   :billow r/billow-noise
   :ridgemulti r/ridgedmulti-noise})

(def ^:private generators (keys generator-fns))

;; deprecated vars are resolved at run time (no linter warning, they stay tested)
(def ^:private random-noise-fn* (requiring-resolve 'fastmath.random/random-noise-fn))
(def ^:private warp-noise-fn* (requiring-resolve 'fastmath.random/warp-noise-fn))

(defn- noise-fn
  "Noise function of the generator `gen` for the configuration map."
  [gen cfg]
  ((generator-fns gen) cfg))

(defn- value
  "Value of the noise function `f` at the point `pt` (1 to 3 coordinates), as a double."
  ^double [f pt]
  (double (apply f pt)))

(defn- bits ^long [^double x] (Double/doubleToLongBits x))

(defn- same-bits? [a b] (= (bits (double a)) (bits (double b))))

(defn- ex-data-of
  "`ex-data` of the `ex-info` thrown by `f`, `::other` for another exception, `nil` when nothing is thrown."
  [f]
  (try (f) nil
       (catch clojure.lang.ExceptionInfo e (or (ex-data e) {}))
       (catch Throwable _ ::other)))

(defn- invalid-input?
  "True when `f` throws `ex-info` of the type `:invalid-input`; with `property` it must be the `:property` of the error
  (a non-keyword `property` only requires some property)."
  ([f] (let [data (ex-data-of f)] (and (map? data) (= :invalid-input (:type data)) (some? (:property data)))))
  ([property f] (let [data (ex-data-of f)] (and (map? data) (= :invalid-input (:type data)) (= property (:property data))))))

(defn- ex-info?
  [f]
  (map? (ex-data-of f)))

(defn- z-score
  "Standardised distance of the observed share `k / n` from the probability `p`."
  ^double [^long k ^long n ^double p]
  (/ (- (/ (double k) n) p) (m/sqrt (/ (* p (- 1.0 p)) n))))

(defn- points
  "`n` seeded random points of dimension `dim`, coordinates uniform in `[-lim, lim]`."
  [^long dim ^long n ^double lim ^long seed]
  (let [rnd (SplittableRandom. seed)]
    (vec (repeatedly n (fn [] (vec (repeatedly dim #(* lim (- (* 2.0 (.nextDouble rnd)) 1.0)))))))))

;; ---------------------------------------------------------------------------------------------
;; 3.11: the 216-case reference (library vs independent Python implementation)

(deftest reference-216-cases
  (let [{:keys [cases configs points]} @reference
        seen (atom {})]
    (is (= 216 (count cases)))
    (is (= #{1 2 3} (set (keys points))))
    (is (every? #(= 12 (count (points %))) [1 2 3]))
    (doseq [{:keys [generator noise-type interpolation dim normalize? values]} cases
            :let [cfg (assoc (configs normalize?) :noise-type noise-type :interpolation interpolation)
                  f (noise-fn generator cfg)
                  actual (mapv #(value f %) (points dim))
                  diff (double (apply max (map (fn [a e] (m/abs (- (double a) (double e)))) actual values)))]]
      (swap! seen update generator (fnil max 0.0) diff)
      ;; observed max difference 1.3e-13 (ridged simplex 3-D): different order of floating point operations in Python and Java
      (is (< diff 1.0e-12) (str generator " " noise-type " " interpolation " dim " dim " normalize? " normalize? " diff " diff)))
    (is (every? #(< (double %) 1.0e-12) (vals @seen)))
    (testing "the cases cover every combination"
      (is (= 216 (count (set (map (juxt :generator :noise-type :interpolation :dim :normalize?) cases)))))
      (is (= #{:single :fbm :billow :ridgemulti} (set (map :generator cases))))
      (is (= #{:none :linear :hermite :quintic} (set (map :interpolation cases)))))))

(deftest reference-values-are-not-trivial
  (testing "the reference distinguishes the cases (a constant function would pass an equality check against equal values)"
    (let [{:keys [cases]} @reference]
      (is (> (count (distinct (map (comp first :values) cases))) 100))
      (is (every? (fn [c] (> (count (distinct (:values c))) 4)) cases)))))

;; ---------------------------------------------------------------------------------------------
;; 3.2 / 3.3: the baseline of the code before the audit

(deftest baseline-unchanged-values
  (let [cases @baseline
        ridged-new-normaliser? (fn [{:keys [generator cfg]}]
                                 (and (= :ridgemulti generator)
                                      (not (and (== 2.0 (:lacunarity cfg)) (== 0.5 (:gain cfg))))))]
    (is (= 96 (count cases)))
    (testing "bitwise equal to the pre-change values (every case except ridged with lacunarity 2.1 / gain 0.45)"
      (doseq [{:keys [generator cfg values] :as c} cases
              :when (not (ridged-new-normaliser? c))]
        (let [f (noise-fn generator cfg)]
          (doseq [[d expected] values]
            (is (every? true? (map (fn [pt e] (same-bits? e (value f pt))) (noise-points d) expected))
                (str generator " " cfg " dim " d))))))
    (testing "ridged, lacunarity 2.1 / gain 0.45: new = old * (sum gain^k / sum lacunarity^-k)"
      (let [cs (filter ridged-new-normaliser? cases)]
        (is (= 12 (count cs)))
        (doseq [{:keys [cfg values]} cs
                :let [{:keys [octaves lacunarity gain normalize?]} cfg
                      f (noise-fn :ridgemulti cfg)
                      k (/ (double (reduce + (map #(m/pow gain %) (range octaves))))
                           (double (reduce + (map #(m/pow lacunarity (- (double %))) (range octaves)))))]
                [d expected] values
                [pt e] (map vector (noise-points d) expected)
                :let [old (double e)
                      ;; normalised: s * bounding; raw: 2 s * bounding - 1
                      expected-new (if normalize? (* k old) (dec (* k (inc old))))]]
          (is (m/delta-eq expected-new (value f pt) 1.0e-14) (str cfg " " pt)))
        (is (every? (fn [{:keys [cfg values]}]
                      (let [f (noise-fn :ridgemulti cfg)]
                        (some (fn [[d es]] (some true? (map (fn [pt e] (> (m/abs (- (value f pt) (double e))) 1.0e-3)) (noise-points d) es)))
                              values)))
                    cs)
            "the new normaliser changes the values visibly (k = 0.958)")))))

;; ---------------------------------------------------------------------------------------------
;; 3.5: validation of the configuration

(def ^:private invalid-configs
  "[label constructor-argument expected-property]; a property of `nil` only requires some property."
  [["noise-type :perlin" {:noise-type :perlin} :noise-type]
   ["noise-type :cellular" {:noise-type :cellular} :noise-type]
   ["noise-type string" {:noise-type "gradient"} :noise-type]
   ["noise-type number" {:noise-type 1} :noise-type]
   ["noise-type capitalised" {:noise-type :Gradient} :noise-type]
   ["interpolation :cubic" {:interpolation :cubic} :interpolation]
   ["interpolation :smooth" {:interpolation :smooth} :interpolation]
   ["interpolation string" {:interpolation "hermite"} :interpolation]
   ["interpolation number" {:interpolation 2} :interpolation]
   ["octaves 0" {:octaves 0} :octaves]
   ["octaves -1" {:octaves -1} :octaves]
   ["octaves 65" {:octaves 65} :octaves]
   ["octaves 2.5" {:octaves 2.5} :octaves]
   ["octaves string" {:octaves "a"} :octaves]
   ["octaves NaN" {:octaves ##NaN} :octaves]
   ["octaves 1e10" {:octaves 1.0e10} :octaves]
   ["octaves Integer/MAX_VALUE" {:octaves Integer/MAX_VALUE} :octaves]
   ["lacunarity 0" {:lacunarity 0.0} :lacunarity]
   ["lacunarity negative" {:lacunarity -2.0} :lacunarity]
   ["lacunarity NaN" {:lacunarity ##NaN} :lacunarity]
   ["lacunarity Inf" {:lacunarity ##Inf} :lacunarity]
   ["lacunarity string" {:lacunarity "2"} :lacunarity]
   ["gain negative" {:gain -0.1} :gain]
   ["gain NaN" {:gain ##NaN} :gain]
   ["gain Inf" {:gain ##Inf} :gain]
   ["gain -Inf" {:gain ##-Inf} :gain]
   ["gain string" {:gain "0.5"} :gain]
   ["normalize? 0" {:normalize? 0} :normalize?]
   ["normalize? 1" {:normalize? 1} :normalize?]
   ["normalize? keyword" {:normalize? :yes} :normalize?]
   ["normalize? string" {:normalize? "true"} :normalize?]
   ["seed double" {:seed 1.5} :seed]
   ["seed string" {:seed "a"} :seed]
   ["seed NaN" {:seed ##NaN} :seed]
   ["seed keyword" {:seed :a} :seed]
   ["seed ratio" {:seed 1/2} :seed]
   ["config a vector" [1 2] nil]
   ["config a number" 5 nil]
   ["config a string" "x" nil]
   ["config a keyword" :fbm nil]
   ["config a sequence" '(:seed 1) nil]])

(deftest config-validation-table
  (is (>= (count invalid-configs) 25))
  (doseq [gen generators
          [label cfg property] invalid-configs]
    (testing (str gen " " label)
      (if property
        (is (invalid-input? property #(noise-fn gen cfg)) (str gen " " (pr-str cfg)))
        (is (invalid-input? #(noise-fn gen cfg)) (str gen " " (pr-str cfg))))))
  (testing "errors are raised at construction (not at the call): the constructor throws, not the evaluation"
    (doseq [[_ cfg _] invalid-configs]
      (is (ex-info? #(r/fbm-noise cfg)))))
  (testing "ridged noise requires gain > 0 (the weight feedback divides by it), the other generators accept 0"
    (is (invalid-input? :gain #(r/ridgedmulti-noise {:gain 0.0})))
    (is (invalid-input? :gain #(r/ridgedmulti-noise {:gain 0})))
    (is (Double/isFinite (value (r/fbm-noise {:gain 0.0 :seed 1}) [0.3 0.4])))
    (is (Double/isFinite (value (r/billow-noise {:gain 0 :seed 1}) [0.3 0.4])))
    (is (Double/isFinite (value (r/single-noise {:gain 0.0 :seed 1}) [0.3 0.4])))))

(deftest config-valid-inputs
  (testing "boundary values are accepted"
    (doseq [cfg [{:octaves 1} {:octaves 64 :lacunarity 1.0} {:octaves 3N} {:octaves (int 4)} {:octaves (long 4)}
                 {:lacunarity 0.5} {:lacunarity 1e-3} {:lacunarity 2} {:lacunarity 1/2 :octaves 3}
                 {:gain 1e-9} {:gain 3.0} {:gain 1}
                 {:seed 0} {:seed -1} {:seed Long/MAX_VALUE} {:seed Long/MIN_VALUE} {:seed 3N} {:seed (int 7)}
                 {:normalize? false} {:normalize? true}
                 {:foo 1} {}]
            gen generators]
      (is (Double/isFinite (value (noise-fn gen cfg) [0.3 0.4])) (str gen " " (pr-str cfg)))))
  (testing "64 octaves with lacunarity 2: the limit is 2^-33 (maxCoord), the function is defined near the origin"
    (doseq [gen generators
            :let [f (noise-fn gen {:octaves 64 :seed 3})]]
      (is (Double/isFinite (value f [1.0e-11 1.0e-11])) (str gen))
      (is (m/nan? (value f [0.3 0.4])) (str gen))))
  (testing "nil configuration and zero arity"
    (doseq [gen generators]
      (is (<= -0.0065 (value (noise-fn gen nil) [0.3 0.4]) 1.0065))
      (is (<= -0.0065 (value ((generator-fns gen)) [0.3 0.4]) 1.0065))))
  (testing "a nil option behaves as a missing one"
    (doseq [gen generators
            :let [base {:seed 99 :noise-type :value}
                  v (fn [cfg] (value (noise-fn gen cfg) [0.3 0.4 0.5]))]]
      (is (same-bits? (v base) (v (assoc base :octaves nil))) (str gen " octaves"))
      (is (same-bits? (v base) (v (assoc base :lacunarity nil))) (str gen " lacunarity"))
      (is (same-bits? (v base) (v (assoc base :gain nil))) (str gen " gain"))
      (is (same-bits? (v base) (v (assoc base :normalize? nil))) (str gen " normalize?"))
      (is (same-bits? (v base) (v (assoc base :interpolation nil))) (str gen " interpolation"))
      (is (same-bits? (v (dissoc base :noise-type)) (v (assoc base :noise-type nil))) (str gen " noise-type"))
      (is (same-bits? (v (assoc base :octaves 6 :lacunarity 2.0 :gain 0.5 :normalize? true :interpolation :hermite))
                      (v (assoc base :octaves nil :lacunarity nil :gain nil :normalize? nil :interpolation nil)))
          (str gen " all defaults")))
    (testing "a nil seed means a random seed"
      (let [fs (repeatedly 5 #(r/fbm-noise {:seed nil}))
            vs (map #(value % [0.3 0.4 0.5]) fs)]
        (is (every? #(<= -0.0065 (double %) 1.0065) vs))
        (is (> (count (distinct vs)) 1))))))

(deftest seed-is-32-bit
  (testing "seeds outside the int range are truncated to 32 bits (documented), the result is deterministic"
    (let [v (fn [seed] (mapv #(value (r/fbm-noise {:seed seed}) %) [[0.3] [0.3 0.4] [0.3 0.4 0.5]]))]
      (is (= (v 1) (v 4294967297)))
      (is (= (v 1) (v (+ 1 (* 7 4294967296)))))
      (is (= (v -1) (v Long/MAX_VALUE)))
      (is (= (v 0) (v Long/MIN_VALUE)))
      (is (= (v Integer/MIN_VALUE) (v 2147483648)))
      (is (not= (v 1) (v 2)))
      (is (not= (v 1) (v -1)))
      (is (not= (v 2147483648) (v 2147483649)))
      (is (= (v 5) (v 5N)) "BigInt")))
  (testing "docstring says 32 bits"
    (is (re-find #"(?i)32" (:doc (meta #'r/fbm-noise))))))

;; ---------------------------------------------------------------------------------------------
;; 3.6: presets

(def ^:private preset-points
  [[0.3] [0.3 0.4] [0.3 0.4 0.5] [-2.7 5.1] [10.5 -3.3 7.7] [100.25]])

(def ^:private preset-goldens
  "Output of the presets for `preset-points` (computed from the library after the audit; seed 1337)."
  {:noise [0.23192126984126982 0.5103539094349206 0.5567364845307936 0.4413407644850797 0.5489575183847618 0.5262896825396826]
   :vnoise [0.3306977386555427 0.71398856148373 0.5328319719246631 0.3448364715293522 0.5444553636560396 0.5561008117942512]
   :simplex [0.24500770455619048 0.6084496923106212 0.3494353439153439 0.3516570697490189 0.7189892711017042 0.5809174601236979]})

(def ^:private presets
  {:noise [r/noise {:noise-type :gradient :interpolation :quintic}]
   :vnoise [r/vnoise {:noise-type :value :interpolation :hermite}]
   :simplex [r/simplex {:noise-type :simplex :interpolation :hermite}]})

(deftest presets-golden-values
  (doseq [[k [f _]] presets]
    (testing (name k)
      ;; 1e-12 tolerates last-ulp differences of Math.pow between JVMs in the value table
      (is (every? true? (map (fn [pt e] (m/delta-eq (value f pt) e 1.0e-12)) preset-points (preset-goldens k)))
          (str k " " (mapv #(value f %) preset-points))))))

(deftest presets-equal-fbm-with-the-fixed-seed
  (doseq [[k [f cfg]] presets
          :let [g (r/fbm-noise (assoc cfg :seed 1337))]]
    (testing (name k)
      (doseq [d [1 2 3]
              pt (points d 500 50.0 (+ d 17))]
        (is (same-bits? (value f pt) (value g pt)) (str k " " pt)))
      (doseq [pt preset-points]
        (is (same-bits? (value f pt) (value g pt)))))
    (testing "the seed is not 1336 or 1338"
      (is (not= (value f [0.3 0.4]) (value (r/fbm-noise (assoc cfg :seed 1338)) [0.3 0.4]))))))

(deftest presets-are-documented-and-fixed
  (doseq [v [#'r/noise #'r/vnoise #'r/simplex]]
    (is (str/includes? (:doc (meta v)) "1337") (str (:name (meta v)) " documents the seed"))))

;; ---------------------------------------------------------------------------------------------
;; 3.7: random-noise-cfg / random-noise / random-noise-fn with :rng

(defn- default-rng-untouched?
  "True when calling `f` leaves the shared `default-rng` stream untouched (compared with a control run)."
  [f]
  (r/set-seed! 4321)
  (f)
  (let [after (r/drand)]
    (r/set-seed! 4321)
    (= after (r/drand))))

(deftest random-noise-cfg-with-rng
  (testing "equal seeds give equal configurations, the cfg has no :rng"
    (doseq [seed [1 2 3 99 12345]]
      (let [a (r/random-noise-cfg {:rng (r/rng :mersenne seed)})
            b (r/random-noise-cfg {:rng (r/rng :mersenne seed)})]
        (is (= a b) (str seed))
        (is (not (contains? a :rng)))
        (is (= #{:seed :generator :noise-type :interpolation :octaves :lacunarity :gain :warp-scale :warp-depth :normalize?}
               (set (keys a)))))))
  (testing "different seeds differ, consecutive draws from one rng differ"
    (is (not= (r/random-noise-cfg {:rng (r/rng :mersenne 1)}) (r/random-noise-cfg {:rng (r/rng :mersenne 2)})))
    (let [g (r/rng :mersenne 1)]
      (is (not= (r/random-noise-cfg {:rng g}) (r/random-noise-cfg {:rng g})))))
  (testing "other generators as :rng, pre-config keys override"
    (doseq [k [:jdk :isaac :well512a]]
      (is (= (r/random-noise-cfg {:rng (r/rng k 5)}) (r/random-noise-cfg {:rng (r/rng k 5)}))))
    (let [c (r/random-noise-cfg {:rng (r/rng :mersenne 1) :seed 5 :octaves 3 :generator :billow})]
      (is (= 5 (:seed c)))
      (is (= 3 (:octaves c)))
      (is (= :billow (:generator c)))
      (is (not (contains? c :rng))))
    (let [c (r/random-noise-cfg {:rng (r/rng :mersenne 1) :custom 42})]
      (is (= 42 (:custom c)) "unknown keys are kept")))
  (testing "no :rng: a fresh generator, results differ, default-rng is untouched"
    (is (not= (r/random-noise-cfg) (r/random-noise-cfg)))
    (is (not= (r/random-noise-cfg nil) (r/random-noise-cfg {})))
    (is (not= (r/random-noise-cfg {:rng nil}) (r/random-noise-cfg {:rng nil})) "nil rng: a fresh generator")
    (is (default-rng-untouched? #(dotimes [_ 20] (r/random-noise-cfg))))
    (is (default-rng-untouched? #(dotimes [_ 20] (r/random-noise-cfg {:octaves 4}))))
    (is (default-rng-untouched? #(dotimes [_ 20] (r/random-noise-cfg {:rng (r/rng :mersenne 1)})))))
  (testing "reproducible with :rng whatever the state of the default rng and of clojure.core/rand (no rand-nth)"
    (let [c (fn [] (r/set-seed! 7) (r/random-noise-cfg {:rng (r/rng :mersenne 3)}))]
      (is (= (c) (c))))))

(deftest random-noise-cfg-frequencies
  ;; 10^5 configurations from one generator; shares compared with the docstring probabilities, |z| < 4
  (let [n 100000
        g (r/rng :mersenne 20241021)
        cfgs (vec (repeatedly n #(r/random-noise-cfg {:rng g})))
        share (fn [pred] (count (filter pred cfgs)))]
    (testing "uniform choices"
      (doseq [k [:single :fbm :billow :ridgemulti]]
        (is (< (m/abs (z-score (share #(= k (:generator %))) n 0.25)) 4.0) (str k)))
      (doseq [k [:value :gradient :simplex]]
        (is (< (m/abs (z-score (share #(= k (:noise-type %))) n (/ 1.0 3.0))) 4.0) (str k)))
      (doseq [k [:none :linear :hermite :quintic]]
        (is (< (m/abs (z-score (share #(= k (:interpolation %))) n 0.25)) 4.0) (str k)))
      (doseq [o (range 1 10)]
        (is (< (m/abs (z-score (share #(= o (:octaves %))) n (/ 1.0 9.0))) 4.0) (str "octaves " o)))
      (is (= #{1 2 3 4 5 6 7 8 9} (set (map :octaves cfgs)))))
    (testing "warp scale: 0.0 w.p. 0.8, 4.0 w.p. 0.1, otherwise uniform in [0.1, 10)"
      (let [zero (share #(== 0.0 (double (:warp-scale %))))
            four (share #(== 4.0 (double (:warp-scale %))))
            other (filter #(let [s (double (:warp-scale %))] (and (not= s 0.0) (not= s 4.0))) cfgs)]
        (is (< (m/abs (- (/ zero (double n)) 0.8)) 0.01))
        (is (< (m/abs (z-score zero n 0.8)) 4.0))
        (is (< (m/abs (z-score four n 0.1)) 4.0))
        (is (< (m/abs (z-score (count other) n 0.1)) 4.0))
        (is (every? #(<= 0.1 (double (:warp-scale %)) 10.0) other))
        (is (< (m/abs (- (/ (reduce + (map #(double (:warp-scale %)) other)) (count other)) 5.05)) 0.15) "mean of U(0.1, 10)")))
    (testing "warp depth: 1 w.p. 0.8 + 0.2/3, 2 and 3 w.p. 0.2/3 each"
      (is (< (m/abs (z-score (share #(= 1 (:warp-depth %))) n (+ 0.8 (/ 0.2 3.0)))) 4.0))
      (is (< (m/abs (z-score (share #(= 2 (:warp-depth %))) n (/ 0.2 3.0))) 4.0))
      (is (< (m/abs (z-score (share #(= 3 (:warp-depth %))) n (/ 0.2 3.0))) 4.0))
      (is (= #{1 2 3} (set (map :warp-depth cfgs)))))
    (testing "uniform doubles: lacunarity [1.5, 2.5), gain [0.2, 0.8)"
      (let [lac (map #(double (:lacunarity %)) cfgs)
            gain (map #(double (:gain %)) cfgs)]
        (is (every? #(and (<= 1.5 (double %)) (< (double %) 2.5)) lac))
        (is (every? #(and (<= 0.2 (double %)) (< (double %) 0.8)) gain))
        ;; mean of U(a, b): sd = (b - a) / sqrt(12 n)
        (is (< (m/abs (- (/ (reduce + lac) n) 2.0)) (* 4.0 (/ 1.0 (m/sqrt (* 12.0 n))))))
        (is (< (m/abs (- (/ (reduce + gain) n) 0.5)) (* 4.0 (/ 0.6 (m/sqrt (* 12.0 n))))))))
    (testing "always normalized, seeds are integers of the int range"
      (is (every? #(true? (:normalize? %)) cfgs))
      (is (every? #(and (integer? (:seed %)) (<= Integer/MIN_VALUE (long (:seed %)) Integer/MAX_VALUE)) cfgs))
      (is (> (count (distinct (map :seed cfgs))) (* 0.99 n))))))

(def ^:private fixed-cfg
  {:seed 17 :noise-type :gradient :interpolation :quintic :octaves 4 :lacunarity 2.0 :gain 0.5 :warp-scale 0.0 :warp-depth 1})

(deftest random-noise-with-rng
  (doseq [[nm mk] [["random-noise" r/random-noise] ["random-noise-fn" random-noise-fn*]]]
    (testing nm
      (doseq [seed (range 1 31)
              :let [f1 (mk {:rng (r/rng :mersenne seed)})
                    f2 (mk {:rng (r/rng :mersenne seed)})]]
        (is (every? #(same-bits? (value f1 %) (value f2 %)) (concat (points 1 5 5.0 seed) (points 2 5 5.0 seed) (points 3 5 5.0 seed)))
            (str "seed " seed)))
      (testing "different seeds give different functions"
        (let [vs (for [seed (range 1 31)] (value (mk {:rng (r/rng :mersenne seed)}) [0.3 0.4 0.5]))]
          (is (> (count (distinct vs)) 25))))
      (testing "default-rng untouched without :rng and with :rng"
        (is (default-rng-untouched? #(dotimes [_ 10] (mk))))
        (is (default-rng-untouched? #(dotimes [_ 10] (mk nil))))
        (is (default-rng-untouched? #(dotimes [_ 10] (mk {:rng (r/rng :mersenne 3)})))))
      (testing "a function of 1, 2 and 3 coordinates"
        (let [f (mk {:rng (r/rng :mersenne 5)})]
          (is (number? (f 0.3)))
          (is (number? (f 0.3 0.4)))
          (is (number? (f 0.3 0.4 0.5)))))))
  (testing "the returned function is the generator of the cfg (:rng is not part of the configuration)"
    (doseq [gen generators
            :let [cfg (assoc fixed-cfg :generator gen)
                  f (r/random-noise (assoc cfg :rng (r/rng :mersenne 1)))
                  g (noise-fn gen cfg)]]
      (is (every? #(same-bits? (value f %) (value g %)) (points 2 100 40.0 3)) (str gen)))
    (testing "with warp"
      (let [cfg (assoc fixed-cfg :generator :fbm :warp-scale 4.0 :warp-depth 2)
            f (r/random-noise (assoc cfg :rng (r/rng :mersenne 1)))
            g (r/warp-noise (r/fbm-noise cfg) 4.0 2)]
        (is (every? #(same-bits? (value f %) (value g %)) (points 3 100 40.0 4)))))))

(deftest random-noise-generator-option
  (testing ":ridgedmulti is an alias of :ridgemulti"
    (doseq [nm [r/random-noise random-noise-fn*]
            :let [a (nm (assoc fixed-cfg :generator :ridgedmulti))
                  b (nm (assoc fixed-cfg :generator :ridgemulti))
                  c (r/ridgedmulti-noise fixed-cfg)]]
      (is (every? #(same-bits? (value a %) (value c %)) (points 2 100 40.0 5)))
      (is (every? #(same-bits? (value b %) (value c %)) (points 2 100 40.0 5)))))
  (testing "the generator keys of noise-generators"
    (is (= #{:fbm :single :billow :ridgemulti} (set (keys r/noise-generators)))))
  (testing "unknown generator: ex-info, in every form of the option"
    (doseq [g [:zzz :perlin "fbm" 1 :FBM :ridgedmulti-noise]
            nm [r/random-noise random-noise-fn*]]
      (is (ex-info? #(nm {:generator g})) (str (pr-str g)))
      (is (ex-info? #(nm {:generator g :rng (r/rng :mersenne 1)})) (str (pr-str g)))))
  (testing "nil generator means a random one"
    (doseq [nm [r/random-noise random-noise-fn*]]
      (is (number? ((nm {:generator nil}) 0.3 0.4)))
      (is (number? ((nm {:generator nil :rng (r/rng :mersenne 4)}) 0.3 0.4)))
      (let [vs (for [seed (range 40)]
                 (value (nm (assoc fixed-cfg :generator nil :rng (r/rng :mersenne seed))) [0.3 0.4]))]
        (is (> (count (distinct vs)) 3) "random generators give different values"))
      (is (every? (fn [seed]
                    (let [a (nm {:generator nil :rng (r/rng :mersenne seed)})
                          b (nm {:generator nil :rng (r/rng :mersenne seed)})]
                      (same-bits? (value a [0.3 0.4]) (value b [0.3 0.4]))))
                  (range 20))))
    (let [gens (set (for [seed (range 200)] (:generator (r/random-noise-cfg {:generator nil :rng (r/rng :mersenne seed)}))))]
      (is (= #{:single :fbm :billow :ridgemulti} gens)))))

;; ---------------------------------------------------------------------------------------------
;; 3.8: warp-noise

(defn- counting-noise
  "A noise function of 1 to 3 arguments counting its calls."
  [counter]
  (fn
    ([x] (swap! counter inc) (* 0.5 (double x)))
    ([x y] (swap! counter inc) (* 0.5 (+ (double x) (double y))))
    ([x y z] (swap! counter inc) (* 0.5 (+ (double x) (double y) (double z))))))

(deftest warp-noise-validation
  (let [n (r/fbm-noise {:seed 1})]
    (doseq [[nm warp] [["warp-noise" r/warp-noise] ["warp-noise-fn" warp-noise-fn*]]]
      (testing nm
        (testing "depth"
          (doseq [depth [-1 -100 2.5 ##NaN ##Inf "a" :a nil 1/2]]
            (is (ex-info? #(warp n 4.0 depth)) (str "depth " (pr-str depth)))))
        (testing "scale"
          (doseq [scale [##NaN ##Inf ##-Inf "a" :a nil]]
            (is (ex-info? #(warp n scale 1)) (str "scale " (pr-str scale)))
            (is (ex-info? #(warp n scale)) (str "scale " (pr-str scale) " 2 args"))))
        (testing "noise"
          (doseq [bad [nil 5 "f" :k {:a 1} [1 2]]]
            (is (ex-info? #(warp bad 4.0 1)) (str "noise " (pr-str bad)))
            (is (ex-info? #(warp bad)) (str "noise " (pr-str bad) " 1 arg"))))
        (testing "valid arguments are accepted: 0 depth, negative and zero scale, integer-valued numbers"
          (is (fn? (warp n 4.0 0)))
          (is (fn? (warp n 0.0 3)))
          (is (fn? (warp n -2.0 1)))
          (is (fn? (warp n 2 2)))
          (is (fn? (warp n 4.0 (int 1))))
          (is (fn? (warp n)))
          (is (fn? (warp n 3.0)))
          (is (fn? (warp))))))
    (testing "the error is raised at construction, with :invalid-input"
      (is (invalid-input? #(r/warp-noise n 4.0 -1)))
      (is (invalid-input? #(r/warp-noise n ##NaN 1))))))

(deftest warp-noise-depth-zero-and-scale-zero
  (doseq [[nm warp] [["warp-noise" r/warp-noise] ["warp-noise-fn" warp-noise-fn*]]
          [gen cfg] [[:fbm {:seed 4}] [:billow {:seed 5 :noise-type :value}] [:ridgemulti {:seed 6 :noise-type :simplex}]]
          :let [n (noise-fn gen cfg)]]
    (testing (str nm " " gen)
      (doseq [d [1 2 3]
              pt (points d 300 30.0 d)]
        (is (same-bits? (value n pt) (value (warp n 4.0 0) pt)) "depth 0 is the identity")
        (is (same-bits? (value n pt) (value (warp n 0.0 2) pt)) "scale 0 is the unwarped noise")
        (is (same-bits? (value n pt) (value (warp n -0.0 3) pt)) "scale -0.0"))
      (testing "a real warp changes the values"
        (is (some #(not (same-bits? (value n %) (value (warp n 4.0 1) %))) (points 2 50 30.0 2)))))))

(deftest warp-noise-evaluation-count
  (doseq [[nm warp] [["warp-noise" r/warp-noise] ["warp-noise-fn" warp-noise-fn*]]
          depth (range 5)]
    (testing (str nm " depth " depth)
      (let [counter (atom 0)
            w (warp (counting-noise counter) 4.0 depth)]
        (w 0.3)
        (is (= (inc depth) @counter) "1-D: depth + 1")
        (reset! counter 0)
        (w 0.3 0.4)
        (is (= (dec (long (m/pow 2 (inc depth)))) @counter) "2-D: 2^(depth+1) - 1")
        (reset! counter 0)
        (w 0.3 0.4 0.5)
        (is (= (quot (dec (long (m/pow 3 (inc depth)))) 2) @counter) "3-D: (3^(depth+1) - 1) / 2")))))

(deftest warp-noise-keeps-the-range
  (let [n (r/fbm-noise {:seed 3})
        w (r/warp-noise n 4.0 2)]
    (is (every? #(<= -0.0065 (value w %) 1.0065) (concat (points 1 2000 100.0 1) (points 2 2000 100.0 2) (points 3 2000 100.0 3))))))

;; ---------------------------------------------------------------------------------------------
;; 3.11: structure of the noise definitions

(deftest gradient-noise-is-zero-on-the-lattice
  (doseq [it [:linear :hermite :quintic]
          seed [1 2 3 4 5]
          :let [f (r/single-noise {:seed seed :noise-type :gradient :interpolation it :normalize? false})
                g (r/single-noise {:seed seed :noise-type :gradient :interpolation it})]]
    (testing (str it " seed " seed)
      (doseq [x (range -6 7)]
        (is (< (m/abs (double (f x))) 1.0e-12) (str "1-D " x))
        (is (m/delta-eq 0.5 (double (g x)) 1.0e-12))
        (doseq [y (range -4 5)]
          (is (< (m/abs (double (f x y))) 1.0e-12) (str "2-D " x " " y))
          (doseq [z (range -3 4)]
            (is (< (m/abs (double (f x y z))) 1.0e-12) (str "3-D " x " " y " " z)))))))
  (testing "not zero between the lattice points, and not for value noise"
    (let [f (r/single-noise {:seed 1 :noise-type :gradient :normalize? false})]
      (is (some #(> (m/abs (double (f % (- 1.0 %)))) 1.0e-3) [0.3 0.4 0.6 0.7]))
      (is (some #(> (m/abs (double (f %))) 0.1) [0.3 0.5 0.7])))
    (let [f (r/single-noise {:seed 1 :noise-type :value :normalize? false})]
      (is (some #(> (m/abs (double (f % 0))) 0.1) (range -6 7))))))

(defn- sample-at
  "Value of `f` at the first coordinate `x` and the other coordinates `others`."
  ^double [f others ^double x]
  (double (apply f x others)))

(defn- crossing-jumps
  "Jumps of the value (J0) and of the first derivative (J1) of `f` along the first coordinate across the integer `n`
  (the other coordinates are `others`): one-sided values at `n -+ 1e-7`, one-sided second order derivatives with step 1e-3."
  [f ^long n others]
  (let [h 1.0e-3
        delta 1.0e-7
        l0 (sample-at f others (- n delta))
        l1 (sample-at f others (- n delta h))
        l2 (sample-at f others (- n delta (* 2.0 h)))
        r0 (sample-at f others (+ n delta))
        r1 (sample-at f others (+ n delta h))
        r2 (sample-at f others (+ n delta (* 2.0 h)))
        dl (/ (+ (* 3.0 l0) (* -4.0 l1) l2) (* 2.0 h))
        dr (/ (+ (* -3.0 r0) (* 4.0 r1) (- r2)) (* 2.0 h))]
    [(m/abs (- l0 r0)) (m/abs (- dl dr))]))

(defn- median ^double [xs] (let [s (vec (sort xs))] (double (nth s (quot (count s) 2)))))

(deftest continuity-classes
  ;; the class of the interpolation: none = discontinuous, linear = C0, hermite and quintic = C1 (finite differences along x, crossings
  ;; at the integers -9..9 including 0 and the negative ones, 8 random values of the remaining coordinates)
  (doseq [nt [:gradient :value]
          dim [1 2 3]
          :let [rnd (SplittableRandom. 99)
                others (fn [] (vec (repeatedly (dec dim) #(+ 0.1 (* 0.8 (.nextDouble rnd))))))
                crossings (fn [it]
                            (let [f (r/single-noise {:seed 7 :noise-type nt :interpolation it :normalize? false})]
                              (vec (for [n (range -9 10) _ (range (if (= dim 1) 1 8))]
                                     (crossing-jumps f n (others))))))]]
    (testing (str nt " " dim "-D")
      (let [none (crossings :none)
            linear (crossings :linear)
            hermite (crossings :hermite)
            quintic (crossings :quintic)]
        (is (> (median (map first none)) 0.05) (str "none: value jumps, median J0 " (median (map first none))))
        (doseq [[it cs] [[:linear linear] [:hermite hermite] [:quintic quintic]]]
          ;; one-sided values at n -+ 1e-7: the jump of a continuous function is up to 2e-7 * slope; observed max 8e-7
          (is (< (apply max (map first cs)) 1.0e-5) (str it ": continuous, max J0 " (apply max (map first cs)))))
        ;; 1-D: the jump is g0 - g1' of the two gradients +-2 around the lattice point, 0 or 4 with probability 1/2
        (is (> (/ (count (filter #(> (double %) 0.3) (map second linear))) (count linear)) 0.2)
            "linear: derivative jumps at a large share of the crossings")
        (is (> (apply max (map second linear)) 1.0))
        (doseq [[it cs] [[:hermite hermite] [:quintic quintic]]]
          ;; the finite difference error is of the order h^2 f''' ~ 1e-5; observed max 7.9e-5 (gradient quintic 1-D), linear jumps are >= 2
          (is (< (apply max (map second cs)) 1.0e-3) (str it ": C1, max J1 " (apply max (map second cs)))))))))

(deftest simplex-is-continuous
  (let [f (r/single-noise {:seed 7 :noise-type :simplex :normalize? false})
        rnd (SplittableRandom. 5)]
    (doseq [n (range -9 10)
            _ (range 4)
            :let [y (+ 0.1 (* 0.8 (.nextDouble rnd)))]]
      (is (< (m/abs (double (- (double (f (- n 1.0e-7) y)) (double (f (+ n 1.0e-7) y))))) 1.0e-5)))))

(deftest nan-above-max-coordinate
  (let [lim 33554432.0 ;; 2^30 / 2^5 for 6 octaves, lacunarity 2
        types [:value :gradient :simplex]]
    (doseq [t types
            gen generators
            :let [f (noise-fn gen {:seed 1 :noise-type t})]
            d [1 2 3]
            k (range d)
            :let [at (fn [v] (vec (assoc (vec (repeat d 0.3)) k v)))]]
      (testing (str t " " gen " dim " d " coordinate " k)
        (is (m/nan? (value f (at 6.8e7))))
        (is (m/nan? (value f (at -6.8e7))))
        (is (m/nan? (value f (at 1.0e300))))
        (is (m/nan? (value f (at ##Inf))))
        (is (m/nan? (value f (at ##-Inf))))
        (is (m/nan? (value f (at ##NaN))))
        (is (Double/isFinite (value f (at 3.3e7))))
        (is (Double/isFinite (value f (at -3.3e7))))
        (is (Double/isFinite (value f (at lim))))
        (is (m/nan? (value f (at (Math/nextUp lim)))))))
    (testing "values below the limit are the same as before (bitwise, see baseline-unchanged-values) and in the range"
      (let [f (r/fbm-noise {:seed 1})]
        (is (every? #(<= -0.0065 (value f %) 1.0065) (points 3 2000 3.3e7 1)))))
    (testing "the limit depends on lacunarity and octaves"
      (let [f (r/fbm-noise {:seed 1 :octaves 5 :lacunarity 2.1})] ;; 2^30 / 2.1^4 = 5.52e7
        (is (Double/isFinite (value f [5.4e7 0.3])))
        (is (m/nan? (value f [5.6e7 0.3]))))
      (let [f (r/fbm-noise {:seed 1 :octaves 1})]
        (is (Double/isFinite (value f [1.0e9 -1.0e9])))
        (is (m/nan? (value f [1.1e9 0.0])))))
    (testing "the presets"
      (doseq [[nm f] [["noise" r/noise] ["vnoise" r/vnoise] ["simplex" r/simplex]]
              big [6.8e7 -6.8e7 1.0e300 ##Inf ##-Inf ##NaN]]
        (is (m/nan? (double (f big))) (str nm " 1-D " big))
        (is (m/nan? (double (f 0.3 big))) (str nm " 2-D " big))
        (is (m/nan? (double (f 0.3 0.4 big))) (str nm " 3-D " big))
        (is (m/nan? (double (f big 0.3 0.4))) (str nm " 3-D first " big))
        (is (<= -0.0065 (double (f 3.3e7)) 1.0065))
        (is (<= -0.0065 (double (f 3.3e7 -3.3e7 3.3e7)) 1.0065))))
    (testing "single-noise has the same limit (it is evaluated with the same configuration)"
      (let [f (r/single-noise {:seed 1})]
        (is (m/nan? (value f [6.8e7 0.3])))
        (is (m/nan? (value f [##NaN])))
        (is (Double/isFinite (value f [3.3e7 0.3])))))))

(defn- range-stats
  "[min max non-finite-count] of `f` over `n` seeded random points of dimension `dim`, coordinates in `[-lim, lim]`."
  [f dim n lim seed]
  (let [rnd (SplittableRandom. (long seed))
        dim (long dim)
        n (long n)
        lim (double lim)]
    (loop [i 0 mn Double/POSITIVE_INFINITY mx Double/NEGATIVE_INFINITY bad 0]
      (if (< i n)
        (let [x (* lim (- (* 2.0 (.nextDouble rnd)) 1.0))
              y (* lim (- (* 2.0 (.nextDouble rnd)) 1.0))
              z (* lim (- (* 2.0 (.nextDouble rnd)) 1.0))
              v (double (case dim 1 (f x) 2 (f x y) (f x y z)))]
          (if (Double/isFinite v)
            (recur (inc i) (Math/min mn v) (Math/max mx v) bad)
            (recur (inc i) mn mx (inc bad))))
        [mn mx bad]))))

(deftest ranges
  ;; 2 * 10^5 points per combination; gradient noise in 3-D reaches 1.0061 (linear), the bound 1.0065 is the one of the Spec
  (let [combos (for [gen generators t [:value :gradient :simplex] d [1 2 3] [norm it] [[true :linear] [false :quintic]]]
                 [gen t d norm it])
        results (doall (pmap (fn [[gen t d norm it :as combo]]
                               (let [f (noise-fn gen {:seed 11 :noise-type t :interpolation it :normalize? norm})]
                                 [combo (range-stats f d 200000 1000.0 (+ 3 d))]))
                             combos))]
    (is (= 72 (count results)))
    (doseq [[[gen t d norm it] [mn mx bad]] results]
      (let [mn (double mn) mx (double mx)]
        (testing (str gen " " t " dim " d " normalize? " norm " " it " min " mn " max " mx)
          (is (zero? bad))
          (if norm
            (do (is (<= -0.0065 mn)) (is (<= mx 1.0065)))
            (do (is (<= -1.0065 mn)) (is (<= mx 1.0065))))
          (when (= gen :ridgemulti)
            (if norm
              (is (<= -1.0e-12 mn) "normalised ridged noise is not negative")
              (is (<= -1.0 (+ mn 1.0e-12)) "raw ridged noise is 2 s - 1 >= -1"))))))))

(deftest ridged-normaliser
  (testing "gain 0.2 (lacunarity 2): the lattice origin has all ridge signals 1, the value must not exceed 1 (old normaliser: 1.575)"
    (doseq [[lac gain] [[2.0 0.2] [1.5 0.2] [3.0 0.8] [2.0 0.5] [2.4 0.25] [1.0 0.5]]
            :let [f (r/ridgedmulti-noise {:seed 7 :noise-type :gradient :interpolation :hermite :lacunarity lac :gain gain})
                  v (value f [0.0 0.0])]]
      (is (<= v (+ 1.0 1.0e-12)) (str lac " " gain " " v))
      (is (m/delta-eq 1.0 v 1.0e-12) (str lac " " gain " " v))
      (is (m/delta-eq 1.0 (value f [0.0 0.0 0.0]) 1.0e-12))
      (is (m/delta-eq 1.0 (value f [0.0]) 1.0e-12))
      (is (m/delta-eq 1.0 (value (r/ridgedmulti-noise {:seed 7 :noise-type :gradient :interpolation :hermite :lacunarity lac :gain gain
                                                       :normalize? false})
                                 [0.0 0.0])
                      1.0e-12) "raw: 2 * 1 - 1")))
  (testing "the maximum over random points is close to 1 and not above 1.0065"
    (doseq [[lac gain] [[1.5 0.2] [3.0 0.8] [2.0 0.5]]
            :let [f (r/ridgedmulti-noise {:seed 3 :lacunarity lac :gain gain})
                  [_ mx bad] (range-stats f 2 200000 1000.0 11)]]
      (is (zero? bad))
      (is (<= 0.8 (double mx) 1.0065) (str lac " " gain " max " mx)))))

(deftest thread-safety
  (let [fns {:fbm-gradient (r/fbm-noise {:seed 5})
             :billow-value (r/billow-noise {:seed 6 :noise-type :value})
             :ridged-simplex (r/ridgedmulti-noise {:seed 7 :noise-type :simplex})
             :single (r/single-noise {:seed 8})
             :warp (r/warp-noise (r/fbm-noise {:seed 9}) 3.0 1)
             :noise r/noise
             :vnoise r/vnoise
             :simplex r/simplex}
        run (fn [f ^long seed]
              (let [rnd (SplittableRandom. seed)
                    out (double-array 100000)]
                (dotimes [i 100000]
                  (aset out i (double (f (* 100.0 (.nextDouble rnd)) (* 100.0 (.nextDouble rnd)) (* 100.0 (.nextDouble rnd))))))
                out))]
    (doseq [[k f] fns]
      (testing (name k)
        (let [expected (mapv #(run f %) (range 8))
              actual (mapv deref (mapv #(future (run f %)) (range 8)))]
          (is (every? true? (map (fn [^doubles a ^doubles e] (Arrays/equals a e)) actual expected))
              "8 threads give the single-thread values bitwise"))))
    (testing "all 8 threads on the same points"
      (let [f (fns :fbm-gradient)
            expected (run f 1)
            actual (mapv deref (mapv (fn [_] (future (run f 1))) (range 8)))]
        (is (every? #(Arrays/equals ^doubles % ^doubles expected) actual))))
    (testing "construction in several threads"
      (let [vs (mapv deref (mapv (fn [i] (future (value (r/fbm-noise {:seed (mod i 3)}) [0.3 0.4]))) (range 24)))]
        (is (= 3 (count (distinct vs))))
        (is (every? (fn [i] (= (nth vs i) (nth vs (mod i 3)))) (range 24)))))))

;; ---------------------------------------------------------------------------------------------
;; 3.9: discrete-noise

(deftest discrete-noise-range-and-determinism
  (testing "range [0, 1], determinism"
    (let [vs (for [x (range -200 200) y (range -50 50)] (r/discrete-noise x y))]
      (is (every? #(<= 0.0 (double %) 1.0) vs))
      (is (= vs (for [x (range -200 200) y (range -50 50)] (r/discrete-noise x y)))))
    (is (every? #(<= 0.0 (r/discrete-noise %) 1.0) (range -100000 100000))))
  (testing "1-arg = 2-arg with 0"
    (is (every? #(same-bits? (r/discrete-noise %) (r/discrete-noise % 0)) (range -1000 1000)))
    (is (same-bits? (r/discrete-noise 123456789) (r/discrete-noise 123456789 0))))
  (testing "pinned value at the origin (unchanged hash, Hugo Elias integer noise)"
    (is (m/delta-eq 0.6408954922300277 (r/discrete-noise 0 0) 1.0e-15))
    (is (m/delta-eq 0.6408954922300277 (r/discrete-noise 0) 1.0e-15)))
  (testing "values are spread over [0, 1]: mean 1/2 within 4 sigma over 10^5 values"
    (let [n 100000
          mean (/ (reduce + (map #(r/discrete-noise %) (range n))) n)]
      (is (< (m/abs (- mean 0.5)) (* 4.0 (/ (m/sqrt (/ 1.0 12.0)) (m/sqrt n))))))))

(deftest discrete-noise-int-truncation
  (testing "the arguments are 32-bit ints: longs are truncated (documented)"
    (is (same-bits? (r/discrete-noise 1) (r/discrete-noise (+ 4294967296 1))))
    (is (same-bits? (r/discrete-noise 1 2) (r/discrete-noise (+ 4294967296 1) (+ 4294967296 2))))
    (is (same-bits? (r/discrete-noise -1) (r/discrete-noise Long/MAX_VALUE)))
    (is (same-bits? (r/discrete-noise 0) (r/discrete-noise Long/MIN_VALUE)))
    (is (same-bits? (r/discrete-noise Integer/MIN_VALUE) (r/discrete-noise 2147483648)))
    (is (same-bits? (r/discrete-noise 5 7) (r/discrete-noise (+ 5 (* 3 4294967296)) (- 7 4294967296))))
    (is (not (same-bits? (r/discrete-noise 1) (r/discrete-noise 2))))))

(deftest discrete-noise-aliasing
  (testing "the 2-D hash is n = X + 57 Y: (x, y) and (x + 57, y - 1) collide (documented)"
    (doseq [x (range -20 20) y (range -20 20)]
      (is (same-bits? (r/discrete-noise x y) (r/discrete-noise (+ x 57) (dec y))) (str x " " y))
      (is (same-bits? (r/discrete-noise x y) (r/discrete-noise (- x 114) (+ y 2))))
      (is (same-bits? (r/discrete-noise (+ x (* 57 y)) 0) (r/discrete-noise x y)) "depends on X + 57 Y only")))
  (testing "no aliasing without the lattice relation"
    (is (not (same-bits? (r/discrete-noise 0 0) (r/discrete-noise 1 0))))
    (is (not (same-bits? (r/discrete-noise 0 1) (r/discrete-noise 0 0))))
    (is (not (same-bits? (r/discrete-noise 10 3) (r/discrete-noise 11 2)))))
  (testing "a 300 x 300 grid has far fewer distinct values than points (17343 of 90000 measured before the audit)"
    (let [distinct-count (count (distinct (for [x (range 300) y (range 300)] (r/discrete-noise x y))))]
      (is (= 17343 distinct-count))))
  (testing "documented in the docstring"
    (is (str/includes? (:doc (meta #'r/discrete-noise)) "57"))))
