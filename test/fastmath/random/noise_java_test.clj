(ns fastmath.random.noise-java-test
  "Tests of the Java noise classes `fastmath.java.noise.*` (Random Audit, work packages 3.1-3.4), used directly
  without `fastmath.random`: `NoiseConfig` (clamped octave count, `maxCoord`, `ridgedBounding`, removed cellular
  leftovers), `Noise`/`FBM`/`Billow`/`RidgedMulti` (NaN above `maxCoord`, ridged normaliser).

  The oracle for unchanged values is the pre-audit build, recorded in `test/resources/random/baseline_reference.edn`
  (`utils/fastmath/dev/generate_random_baseline.clj`, key `:noise`, 96 cases x 12 points)."
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.test :refer [deftest is testing]]
            [fastmath.core :as m])
  (:import [fastmath.java.noise Billow FBM Noise NoiseConfig RidgedMulti]
           [java.util SplittableRandom]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(def ^:private baseline
  (delay (:noise (edn/read-string (slurp (io/resource "random/baseline_reference.edn"))))))

(def ^:private noise-points
  "The points of the baseline, by dimension (copy of `noise-points` in generate_random_baseline.clj)."
  {1 [[0.3] [-4.7] [123.456] [-0.0001]]
   2 [[0.3 0.4] [-4.7 2.2] [123.456 -78.9] [10.0 10.0]]
   3 [[0.3 0.4 0.5] [-4.7 2.2 9.9] [123.4 -78.9 1.1] [0.5 0.5 0.5]]})

(def ^:private noise-types
  {:value NoiseConfig/NOISE_VALUE
   :gradient NoiseConfig/NOISE_GRADIENT
   :simplex NoiseConfig/NOISE_SIMPLEX})

(def ^:private interpolations
  {:none NoiseConfig/INTERPOLATE_NONE
   :linear NoiseConfig/INTERPOLATE_LINEAR
   :hermite NoiseConfig/INTERPOLATE_HERMITE
   :quintic NoiseConfig/INTERPOLATE_QUINTIC})

(defn- make-cfg
  "`NoiseConfig` from the map used in the baseline file (the mapping of `fastmath.random/noise-config-obj`)."
  ^NoiseConfig [{:keys [seed noise-type interpolation octaves lacunarity gain normalize?]}]
  (NoiseConfig. (int seed) (int (noise-types noise-type)) (int (interpolations interpolation))
                (int octaves) (double lacunarity) (double gain) (boolean normalize?)))

(defn- cfg
  "Default fractal configuration with overrides: seed 1, 6 octaves, lacunarity 2, gain 0.5, normalized."
  (^NoiseConfig [] (cfg {}))
  (^NoiseConfig [opts]
   (make-cfg (merge {:seed 1 :noise-type :gradient :interpolation :hermite :octaves 6
                     :lacunarity 2.0 :gain 0.5 :normalize? true}
                    opts))))

(def ^:private generators
  "Generator -> functions of 1, 2 and 3 coordinates (config first)."
  {:single [(fn ^double [^NoiseConfig c ^double x] (Noise/noise c x))
            (fn ^double [^NoiseConfig c ^double x ^double y] (Noise/noise c x y))
            (fn ^double [^NoiseConfig c ^double x ^double y ^double z] (Noise/noise c x y z))]
   :fbm [(fn ^double [^NoiseConfig c ^double x] (FBM/noise c x))
         (fn ^double [^NoiseConfig c ^double x ^double y] (FBM/noise c x y))
         (fn ^double [^NoiseConfig c ^double x ^double y ^double z] (FBM/noise c x y z))]
   :billow [(fn ^double [^NoiseConfig c ^double x] (Billow/noise c x))
            (fn ^double [^NoiseConfig c ^double x ^double y] (Billow/noise c x y))
            (fn ^double [^NoiseConfig c ^double x ^double y ^double z] (Billow/noise c x y z))]
   :ridgemulti [(fn ^double [^NoiseConfig c ^double x] (RidgedMulti/noise c x))
                (fn ^double [^NoiseConfig c ^double x ^double y] (RidgedMulti/noise c x y))
                (fn ^double [^NoiseConfig c ^double x ^double y ^double z] (RidgedMulti/noise c x y z))]})

(defn- value
  "Value of the generator `gen` at the point `pt` (1 to 3 coordinates)."
  ^double [gen c pt]
  (double (apply (get-in generators [gen (dec (count pt))]) c pt)))

(defn- bits ^long [^double x] (Double/doubleToLongBits x))

(defn- same-bits? [a b] (= (bits (double a)) (bits (double b))))

(defn- abs-diff ^double [a b] (m/abs (- (double a) (double b))))

(defn- sample-stats
  "Minimum, maximum and number of non-finite values of `n` seeded random points, each coordinate
  uniform in `[-lim, lim]` (default `maxCoord`), for `FBM` in `dim` dimensions."
  [^NoiseConfig c ^long dim ^long n ^double lim]
  (let [r (SplittableRandom. 12345)]
    (loop [i 0 mn Double/POSITIVE_INFINITY mx Double/NEGATIVE_INFINITY bad 0]
      (if (< i n)
        (let [x (* lim (- (* 2.0 (.nextDouble r)) 1.0))
              y (* lim (- (* 2.0 (.nextDouble r)) 1.0))
              z (* lim (- (* 2.0 (.nextDouble r)) 1.0))
              v (case dim
                  1 (FBM/noise c x)
                  2 (FBM/noise c x y)
                  (FBM/noise c x y z))]
          (if (Double/isFinite v)
            (recur (inc i) (Math/min mn v) (Math/max mx v) bad)
            (recur (inc i) mn mx (inc bad))))
        [mn mx bad]))))

(defn- ridged-max
  "Maximum of `RidgedMulti` over `n` seeded random 2-D points in `[-1000, 1000]^2`."
  ^double [^NoiseConfig c ^long n]
  (let [r (SplittableRandom. 777)]
    (loop [i 0 mx Double/NEGATIVE_INFINITY]
      (if (< i n)
        (let [x (* 1000.0 (- (* 2.0 (.nextDouble r)) 1.0))
              y (* 1000.0 (- (* 2.0 (.nextDouble r)) 1.0))]
          (recur (inc i) (Math/max mx (RidgedMulti/noise c x y))))
        mx))))

;; ---------------------------------------------------------------------------------------------
;; 3.2 / 3.3: unchanged values

(deftest baseline-bitwise-test
  (let [cases @baseline]
    (is (= 96 (count cases)))
    (doseq [{:keys [generator values] cfg-map :cfg} cases
            :let [c (make-cfg cfg-map)
                  ridged-changed? (and (= :ridgemulti generator)
                                       (not (and (= 2.0 (:lacunarity cfg-map)) (= 0.5 (:gain cfg-map)))))]]
      (when-not ridged-changed?
        (testing (pr-str generator cfg-map)
          (doseq [[d expected] values]
            (is (every? true? (map (fn [pt e] (same-bits? e (value generator c pt))) (noise-points d) expected))
                (str generator " " cfg-map " dim " d))))))))

(deftest baseline-ridged-new-normaliser-test
  (testing "ridged noise with lacunarity 2.1 / gain 0.45: new value = old value * (ridgedBounding / fractalBounding)"
    (let [cases (filter #(and (= :ridgemulti (:generator %)) (= 2.1 (get-in % [:cfg :lacunarity]))) @baseline)
          diffs (for [{:keys [values] cfg-map :cfg} cases
                      :let [c (make-cfg cfg-map)
                            k (/ (.ridgedBounding c) (.fractalBounding c))]
                      [d expected] values
                      [pt e] (map vector (noise-points d) expected)]
                  [(value :ridgemulti c pt) (* k (double e)) e])]
      (is (= 12 (count cases)))
      (is (every? (fn [[a b _]] (m/delta-eq a b 1.0e-14)) diffs))
      (let [largest (double (apply max (map (fn [[a _ old]] (abs-diff a old)) diffs)))]
        (is (< 0.01 largest 1.0) "the new normaliser is not the old one here: |value - old value| must be visible"))))
  (testing "gain 0.5, lacunarity 2: the old formula (fractalBounding) gives bitwise the same ridged values"
    (doseq [norm [true false]
            :let [new-cfg (cfg {:seed 31 :octaves 6 :normalize? norm})
                  old-cfg (cfg {:seed 31 :octaves 6 :normalize? norm})
                  _ (set! (.-ridgedBounding old-cfg) (.fractalBounding old-cfg))
                  r (SplittableRandom. 4)
                  pts (repeatedly 2000 #(vector (* 100.0 (.nextDouble r)) (* 100.0 (.nextDouble r)) (* 100.0 (.nextDouble r))))]]
      (is (= (bits (.ridgedBounding new-cfg)) (bits (.fractalBounding new-cfg))))
      (is (every? (fn [pt] (and (same-bits? (value :ridgemulti new-cfg pt) (value :ridgemulti old-cfg pt))
                                (same-bits? (value :ridgemulti new-cfg (take 2 pt)) (value :ridgemulti old-cfg (take 2 pt)))
                                (same-bits? (value :ridgemulti new-cfg (take 1 pt)) (value :ridgemulti old-cfg (take 1 pt)))))
                  pts)
          (str "normalize? " norm)))))

;; ---------------------------------------------------------------------------------------------
;; 3.1: NoiseConfig

(deftest noiseconfig-fields-test
  (testing "maxCoord = 2^30 / max(1, lacunarity)^(octaves - 1)"
    (is (= 33554432.0 (.maxCoord (cfg))) "6 octaves, lacunarity 2: 2^25")
    (is (= (/ 1073741824.0 (Math/pow 2.1 4.0)) (.maxCoord (cfg {:octaves 5 :lacunarity 2.1}))))
    (is (= (/ 1073741824.0 (Math/pow 3.0 2.0)) (.maxCoord (cfg {:octaves 3 :lacunarity 3.0}))))
    (is (= 1073741824.0 (.maxCoord (cfg {:octaves 1})) (.maxCoord (cfg {:octaves 1 :lacunarity 7.0}))))
    (is (= 1073741824.0 (.maxCoord (cfg {:lacunarity 1.0}))) "lacunarity 1")
    (is (= 1073741824.0 (.maxCoord (cfg {:lacunarity 0.5}))) "lacunarity 0.5 does not enlarge the limit")
    (is (= 1073741824.0 (.maxCoord (cfg {:lacunarity 0.5 :octaves 64}))))
    (is (= (Math/scalb 1.0 (int -33)) (.maxCoord (cfg {:octaves 64 :lacunarity 2.0}))) "64 octaves: 2^30 / 2^63"))
  (testing "octaves are clamped to at least 1 and the arrays are sized by the clamped count"
    (doseq [o [0 -1 -100 Integer/MIN_VALUE]
            :let [c (cfg {:octaves o})]]
      (is (= 1 (.octaves c)) (str "octaves " o))
      (is (= 1 (alength ^doubles (.spectralWeights c))))
      (is (= 1073741824.0 (.maxCoord c)))
      (is (= 1.0 (.fractalBounding c)))
      (is (= 1.0 (.ridgedBounding c)))
      (is (Double/isFinite (value :ridgemulti c [0.3 0.4])))
      (is (Double/isFinite (value :fbm c [0.3 0.4])))
      (is (same-bits? (value :fbm c [0.3 0.4]) (value :fbm (cfg {:octaves 1}) [0.3 0.4])))))
  (testing "64 octaves can be constructed and evaluated"
    (let [c (cfg {:octaves 64 :lacunarity 1.0})]
      (is (= 64 (.octaves c)))
      (is (= 64 (alength ^doubles (.spectralWeights c))))
      (is (every? #(Double/isFinite (value % c [0.3 0.4 0.5])) [:single :fbm :billow :ridgemulti])))
    (let [c (cfg {:octaves 64})]
      (is (= 64 (alength ^doubles (.spectralWeights c))))
      (is (Double/isFinite (value :fbm c [0.0 0.0 0.0])))
      (is (Double/isFinite (value :fbm c [1.0e-10])))
      (is (m/nan? (value :fbm c [1.0e-9])) "above 2^-33 = 1.16e-10")))
  (testing "spectral weights and ridgedBounding"
    (let [c (cfg {:octaves 6 :lacunarity 2.0})]
      (is (= [1.0 0.5 0.25 0.125 0.0625 0.03125] (vec (.spectralWeights c))))
      (is (= (/ 32.0 63.0) (.ridgedBounding c))))
    (let [c (cfg {:octaves 5 :lacunarity 1.5 :gain 0.2})
          w (map #(Math/pow 1.5 (- (double %))) (range 5))]
      (is (every? true? (map #(m/delta-eq %1 %2 1.0e-15) w (.spectralWeights c))))
      (is (m/delta-eq (/ 1.0 (double (reduce + w))) (.ridgedBounding c) 1.0e-15))
      (is (not (m/delta-eq (.fractalBounding c) (.ridgedBounding c) 1.0e-3))
          "gain-based and lacunarity-based normalisers differ"))))

(deftest cellular-removed-test
  (testing "removed members"
    (doseq [f ["cell2dX" "cell2dY" "rev_octaves"]]
      (is (thrown? NoSuchFieldException (.getField NoiseConfig ^String f)) f))
    (is (thrown? NoSuchFieldException (.getField NoiseConfig "NOISE_CELLULAR")))
    (is (thrown? ClassNotFoundException (Class/forName "fastmath.java.noise.Cellular"))))
  (testing "perm and valueLUT are not changed by removing the cell tables (baseline value noise is bitwise equal)"
    (let [c (cfg {:seed 12345 :noise-type :value})]
      (is (= 512 (alength ^ints (.perm c))))
      (is (= 256 (alength ^doubles (.valueLUT c)))))))

;; ---------------------------------------------------------------------------------------------
;; 3.2: NaN above maxCoord

(deftest nan-above-maxcoord-test
  (let [lim (.maxCoord (cfg))
        types [:value :gradient :simplex]]
    (is (= 33554432.0 lim))
    (doseq [t types
            gen [:single :fbm :billow :ridgemulti]
            :let [c (cfg {:noise-type t})]
            d [1 2 3]
            k (range d)
            :let [at (fn [v] (vec (assoc (vec (repeat d 0.3)) k v)))]]
      (testing (str t " " gen " dim " d " coordinate " k)
        (is (m/nan? (value gen c (at 6.8e7))))
        (is (m/nan? (value gen c (at -6.8e7))))
        (is (m/nan? (value gen c (at (Math/nextUp lim)))))
        (is (m/nan? (value gen c (at (- (Math/nextUp lim))))))
        (is (m/nan? (value gen c (at ##Inf))))
        (is (m/nan? (value gen c (at ##-Inf))))
        (is (m/nan? (value gen c (at ##NaN))))
        (is (m/nan? (value gen c (at 1.0e300))))
        (is (Double/isFinite (value gen c (at lim))) "exactly maxCoord is accepted")
        (is (Double/isFinite (value gen c (at (- lim)))))
        (is (Double/isFinite (value gen c (at (Math/nextDown lim)))))
        (is (Double/isFinite (value gen c (at (- 12345.678)))) "negative coordinates")))
    (testing "all coordinates at +-maxCoord (corners)"
      (doseq [t types gen [:single :fbm :billow :ridgemulti] s [1.0 -1.0] d [1 2 3]
              :let [c (cfg {:noise-type t})
                    v (value gen c (vec (repeat d (* (double s) lim))))]]
        (is (Double/isFinite v) (str t " " gen " " s " " d))
        (is (<= -1.0065 v 1.0065) (str t " " gen " " s " " d " " v))))))

(deftest nan-other-configs-test
  (testing "the limit follows the lacunarity and octaves"
    (let [c (cfg {:octaves 5 :lacunarity 2.1})
          lim (.maxCoord c)]
      (is (< 5.5e7 lim 5.6e7))
      (is (Double/isFinite (value :fbm c [lim 0.0 lim])))
      (is (m/nan? (value :fbm c [(Math/nextUp lim) 0.0 0.0])))
      (is (m/nan? (value :billow c [0.0 0.0 (- (Math/nextUp lim))]))))
    (let [c (cfg {:octaves 1})]
      (is (Double/isFinite (value :fbm c [1.0e9 -1.0e9])))
      (is (m/nan? (value :fbm c [(* 2.0 1073741824.0) 0.0]))))
    (let [c (cfg {:lacunarity 1.0})]
      (is (Double/isFinite (value :ridgemulti c [1073741824.0 0.0 0.0])))
      (is (m/nan? (value :ridgemulti c [1073741825.0 0.0 0.0]))))
    (let [c (cfg {:lacunarity 0.5})]
      (is (Double/isFinite (value :billow c [-1073741824.0])))
      (is (m/nan? (value :billow c [-1073741825.0])))))
  (testing "-0.0 is an ordinary coordinate"
    (doseq [t [:value :gradient :simplex] gen [:single :fbm :billow :ridgemulti] d [1 2 3]
            :let [c (cfg {:noise-type t})]]
      (is (== (value gen c (vec (repeat d 0.0))) (value gen c (vec (repeat d -0.0))))
          (str t " " gen " " d)))))

(deftest range-below-maxcoord-test
  (doseq [t [:value :gradient :simplex]
          d [1 2 3]
          :let [c (cfg {:noise-type t})
                [mn mx bad] (sample-stats c d 1000000 (.maxCoord c))
                mn (double mn)
                mx (double mx)]]
    (testing (str t " dim " d " (10^6 points in [-maxCoord, maxCoord]) min " mn " max " mx)
      (is (zero? bad))
      (is (<= -0.0065 mn))
      (is (<= mx 1.0065)))))

;; ---------------------------------------------------------------------------------------------
;; 3.3: ridged normaliser

(deftest ridged-bounding-test
  (testing "gain 0.2: the lattice origin (all ridged signals 1) is not above 1"
    (let [c (cfg {:seed 7 :noise-type :gradient :interpolation :hermite :gain 0.2})
          raw (cfg {:seed 7 :noise-type :gradient :interpolation :hermite :gain 0.2 :normalize? false})
          old (cfg {:seed 7 :noise-type :gradient :interpolation :hermite :gain 0.2})]
      (set! (.-ridgedBounding old) (.fractalBounding old))
      (is (<= (value :ridgemulti c [0.0 0.0]) (+ 1.0 1.0e-12)))
      (is (m/delta-eq 1.0 (value :ridgemulti c [0.0 0.0]) 1.0e-12))
      (is (m/delta-eq 1.0 (value :ridgemulti raw [0.0 0.0]) 1.0e-12))
      (is (> (value :ridgemulti old [0.0 0.0]) 1.5) "the old normaliser gave 2.24 here")))
  (testing "other lacunarity/gain: the maximum over 2*10^6 points is close to 1, never above 1.0065"
    (doseq [[lac gain] [[1.5 0.2] [3.0 0.8] [2.0 0.5]]
            :let [mx (ridged-max (cfg {:seed 3 :lacunarity lac :gain gain}) 2000000)]]
      (testing (str "lacunarity " lac " gain " gain " max " mx)
        (is (<= 0.95 mx 1.0065))))))

(deftest ridged-1d-3d-test
  (testing "ridged noise in 1-D and 3-D stays in [0, 1.0065] for the normalised configuration"
    (doseq [[lac gain] [[1.5 0.2] [3.0 0.8]]
            :let [c (cfg {:seed 5 :lacunarity lac :gain gain})
                  r (SplittableRandom. 8)
                  vs (for [_ (range 20000)]
                       [(RidgedMulti/noise c (* 500.0 (.nextDouble r)))
                        (RidgedMulti/noise c (* 500.0 (.nextDouble r)) (* 500.0 (.nextDouble r)) (* 500.0 (.nextDouble r)))])]]
      (is (every? (fn [[a b]] (let [a (double a) b (double b)] (and (<= -1.0e-12 a 1.0065) (<= -1.0e-12 b 1.0065)))) vs)
          (str "lacunarity " lac " gain " gain)))))
