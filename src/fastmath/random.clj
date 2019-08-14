(ns fastmath.random
  "Various random and noise functions.

  Namespace defines various random number generators (RNGs), different types of random functions, sequence generators and noise functions.

  ### RNGs

  You can use a selection of various RNGs defined in [Apache Commons Math](http://commons.apache.org/proper/commons-math/apidocs/org/apache/commons/math3/random/package-summary.html) library.

  Currently supported RNGs:

  * `:jdk` - default java.util.Random
  * `:mersenne` - MersenneTwister
  * `:isaac` - ISAAC
  * `:well512a`, `:well1024a`, `:well19937a`, `:well19937c`, `:well44497a`, `:well44497b` - several WELL variants

  To create your RNG use [[rng]] multimethod. Pass RNG name and (optional) seed. Returned RNG is equipped with [[RNGProto]] protocol with methods: [[irandom]], [[lrandom]], [[frandom]] [[drandom]], [[grandom]], [[brandom]] which return random primitive value with given RNG.

  ```
  (let [rng (rng :isaac 1337)]
    (irandom rng))
  ```

  For conveniency default RNG (`:jdk`) with following functions are created: [[irand]], [[lrand]], [[frand]], [[drand]], [[grand]], [[brand]].

  Each prefix denotes returned type:

  * i - int
  * l - long
  * f - float
  * d - double
  * g - gaussian (double)
  * b - boolean

  Check individual function for parameters description.

  ### Random Vector Sequences

  Couple of functions to generate sequences of numbers or vectors.

  To create generator call [[sequence-generator]] with generator name and vector size [1,4].
  Following generators are available:

  * `:halton` - Halton low-discrepancy sequence; range [0,1]
  * `:sobol` - Sobol low-discrepancy sequence; range [0,1]
  * `:r2` - R2 low-discrepancy sequence; range [0,1], [more...](http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/)
  * `:sphere` - uniformly random distributed on unit sphere
  * `:gaussian` - gaussian distributed (mean=0, stddev=1)
  * `:default` - uniformly random; range:[0,1]

  `:halton`, `:sobol` and `:r2` can be also randomly jittered according to this [article](http://extremelearning.com.au/a-simple-method-to-construct-isotropic-quasirandom-blue-noise-point-sequences/). Call [[jittered-sequence-generator]].
  
  After creation you get lazy sequence

  ### Noise

  List of continuous noise functions (1d, 2d and 3d):

  * `:value` - value noise
  * `:gradient` - gradient noise (improved Ken Perlin version)
  * `:simplex` - simplex noise

  First two (`:value` and `:gradient`) can use 4 different interpolation types: `:none`, `:linear`, `:hermite` (cubic) and `:quintic`.
  
  All can be used as into:

  * Noise - pure noise value, create with [[single-noise]]
  * FBM - fractal brownian motion, create with [[fbm-noise]]
  * Billow - billow noise, [[billow-noise]]
  * RidgedMulti - ridged multi, [[ridgedmulti-noise]]

  Noise creation requires detailed configuration which is simple map of following keys:

  * `:seed` - seed as integer
  * `:noise-type` - type of noise: `:value`, `:gradient` (default), `:simplex`
  * `:interpolation` - type of interpolation (for value and gradient): `:none`, `:linear`, `:hermite` (default) or `:quintic`
  * `:octaves` - number of octaves for combined noise (like FBM), default: 6
  * `:lacunarity` - scaling factor for combined noise, default: 2.00
  * `:gain` - amplitude scaling factor for combined noise, default: 0.5
  * `:normalize?` - should be normalized to `[0,1]` range (true, default) or to `[-1,1]` range (false)

  For usage convenience 3 ready to use functions are prepared. Return is normalized to `[0,1]` range:

  * [[noise]] - Perlin Noise (gradient noise, 6 octaves, quintic interpolation)
  * [[vnoise]] - Value Noise (as in Processing, 6 octaves, hermite interpolation)
  * [[simplex]] - Simpled Noise (6 octaves)

  #### Discrete Noise

  [[discrete-noise]] is a 1d or 2d hash function for given integers. Returns double from `[0,1]` range.

  ### Distribution

  Various real and integer distributions. See [[DistributionProto]] and [[RNGProto]] for functions.

  To create distribution call [[distribution]] multimethod with name as a keyword and map as parameters."
  {:metadoc/categories {:rand "Random number generation"
                        :noise "Noise functions"
                        :gen "Random sequence generation"
                        :dist "Distributions"}}
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.kernel :as k])
  (:import [org.apache.commons.math3.random RandomGenerator ISAACRandom JDKRandomGenerator MersenneTwister
            Well512a Well1024a Well19937a Well19937c Well44497a Well44497b
            RandomVectorGenerator HaltonSequenceGenerator SobolSequenceGenerator UnitSphereRandomVectorGenerator
            EmpiricalDistribution SynchronizedRandomGenerator]
           [fastmath.java R2]
           [umontreal.ssj.probdist ContinuousDistribution DiscreteDistributionInt InverseGammaDist AndersonDarlingDistQuick ChiDist ChiSquareNoncentralDist CramerVonMisesDist ErlangDist FatigueLifeDist FoldedNormalDist FrechetDist HyperbolicSecantDist InverseGaussianDist HypoExponentialDist HypoExponentialDistEqual JohnsonSBDist JohnsonSLDist JohnsonSUDist KolmogorovSmirnovDistQuick KolmogorovSmirnovPlusDist LogarithmicDist LoglogisticDist NormalInverseGaussianDist Pearson6Dist PowerDist RayleighDist WatsonGDist WatsonUDist]
           [umontreal.ssj.probdistmulti DirichletDist]
           [fastmath.java.noise Billow RidgedMulti FBM NoiseConfig Noise Discrete]
           [smile.stat.distribution Distribution DiscreteDistribution NegativeBinomialDistribution]
           [org.apache.commons.math3.distribution AbstractRealDistribution RealDistribution BetaDistribution CauchyDistribution ChiSquaredDistribution EnumeratedRealDistribution ExponentialDistribution FDistribution GammaDistribution, GumbelDistribution, LaplaceDistribution, LevyDistribution, LogisticDistribution, LogNormalDistribution, NakagamiDistribution, NormalDistribution, ParetoDistribution, TDistribution, TriangularDistribution, UniformRealDistribution WeibullDistribution MultivariateNormalDistribution]
           [org.apache.commons.math3.distribution IntegerDistribution AbstractIntegerDistribution BinomialDistribution EnumeratedIntegerDistribution, GeometricDistribution, HypergeometricDistribution, PascalDistribution, PoissonDistribution, UniformIntegerDistribution, ZipfDistribution]))

(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;; Type hinted functions generating random value
(defn- next-random-value-long
  "Generate next long.

  * arity 0 - from 0 to maximum long value
  * arity 1 - from 0 to provided integer (excluded)
  * arity 2 - from the provided range (included, excluded)"
  (^long [^RandomGenerator r] (.nextLong r))
  (^long [^RandomGenerator r ^long mx] (mod (.nextLong r) mx))
  (^long [r ^long mn ^long mx]
   (let [diff (- mx mn)]
     (if (zero? diff) mn
         (+ mn (next-random-value-long r diff))))))

(defn- next-random-value-double
  "Generate next double.

  * arity 0 - from 0 to 1 (exluded)
  * arity 1 - from 0 to provided double (excluded)
  * arity 2 - from the provided range (included, excluded)"
  (^double [^RandomGenerator r] (.nextDouble r))
  (^double [^RandomGenerator r ^double mx] (* (.nextDouble r) mx))
  (^double [r ^double mn ^double mx]
   (let [diff (- mx mn)]
     (if (zero? diff) mn
         (+ mn (next-random-value-double r diff))))))

(defn- next-random-value-gaussian
  "Generate next random value from normal distribution.

  * arity 0 - N(0,1)
  * arity 1 - N(0,par)
  * arity 2 - N(par1,par2)"
  (^double [^RandomGenerator r] (.nextGaussian r))
  (^double [^RandomGenerator r ^double mx] (* (.nextGaussian r) mx))
  (^double [r ^double mn ^double mx]
   (let [diff (- mx mn)]
     (if (zero? diff) mn
         (+ mn (next-random-value-gaussian r diff))))))

(defprotocol RNGProto
  "Defines set of random functions for different RNGs or distributions returning primitive values."
  (^{:metadoc/categories #{:rand}}
   irandom [t] [t mx] [t mn mx]
   "Random integer.

For RNGs:
As default returns random integer from full integer range. 
When `mx` is passed, range is set to `[0, mx)`. When `mn` is passed, range is set to `[mn, mx)`.

See [[irand]].

For distributions, just returns random integer (call without parameters).")
  (^{:metadoc/categories #{:rand}}
   drandom [t] [t mx] [t mn mx]
   "Random double.

For RNGs:
As default returns random double from `[0,1)` range.
When `mx` is passed, range is set to `[0, mx)`. When `mn` is passed, range is set to `[mn, mx)`.

See [[drand]].

For distributions, just returns random double (call without parameters).")
  (^{:metadoc/categories #{:rand}} lrandom [t] [t mx] [t mn mx]
   "Random long.

For RNGs:
As default returns random long from full long range. 
When `mx` is passed, range is set to `[0, mx)`. When `mn` is passed, range is set to `[mn, mx)`.

See [[lrand]].

For distributions, just returns random long (call without parameters).")
  (^{:metadoc/categories #{:rand}} frandom [t] [t mx] [t mn mx]
   "Random float.

For RNGs:
As default returns random float from `[0,1)` range.
When `mx` is passed, range is set to `[0, mx)`. When `mn` is passed, range is set to `[mn, mx)`.

See [[frand]].

For distributions, just returns random float (call without parameters).")
  (^{:metadoc/categories #{:rand}} grandom [t] [t std] [t mean std]
   "Random double from gaussian distribution.
As default returns random double from `N(0,1)`. 
When `std` is passed, `N(0,std)` is used. When `mean` is passed, distribution is set to `N(mean, std)`.

See [[grand]].")
  (^{:metadoc/categories #{:rand}} brandom [t] [t thr]
   "Boolean random.
Returns true or false with equal probability. You can set probability for `true` setting `thr` (from `[0-1]` range).

See [[brand]].")
  (^{:metadoc/categories #{:rand}} set-seed! [t v] "Sets seed. Returns RNG or distribution itself.")
  (^{:metadoc/categories #{:rand}} ->seq [t] [t n] "Returns sequence of random samples limited to optional `n` values."))

;; Extend RandomGenerator interface with functions created by macro `next-random-value-fn`. This way all RNG classes are enriched with new, more convenient functions.
;;
;; Note that `grandom` is under special care due to different [mn mx] range meaning.

(extend RandomGenerator 
  RNGProto
  {:irandom (comp unchecked-int next-random-value-long)
   :lrandom next-random-value-long
   :frandom (comp float next-random-value-double)
   :drandom next-random-value-double
   :grandom (fn
              ([t] (next-random-value-gaussian t))
              ([t std] (next-random-value-gaussian t std))
              ([t ^double mean ^double std] (next-random-value-gaussian t mean (+ mean std))))
   :brandom (fn
              ([^RandomGenerator t] (.nextBoolean t))
              ([t ^double thr] (< (next-random-value-double t) thr)))
   :set-seed! #(do
                 (.setSeed ^RandomGenerator %1 (long %2))
                 %1)
   :->seq (fn
            ([^RandomGenerator t] (repeatedly #(next-random-value-double t)))
            ([^RandomGenerator t n] (repeatedly n #(next-random-value-double t))))})

;; Helper macro which creates RNG object of given class and/or seed.
(defmacro ^:private create-object-with-seed
  "Create object of the class with (or not) given seed. Used to create RNG."
  [cl seed]
  `(if-let [arg# ~seed]
     (new ~cl (int arg#))
     (new ~cl)))

(defmulti rng
  "Create RNG for given name (as keyword) and optional seed. Return object enhanced with [[RNGProto]]. See: [[rngs-list]] for names."
  {:metadoc/categories #{:rand}}
  (fn [m & _] m))

(defmethod rng :mersenne [_ & [seed]]
  (create-object-with-seed MersenneTwister seed))
(defmethod rng :isaac [_ & [seed]]
  (create-object-with-seed ISAACRandom seed))
(defmethod rng :well512a [_ & [seed]]
  (create-object-with-seed Well512a seed))
(defmethod rng :well1024a [_ & [seed]]
  (create-object-with-seed Well1024a seed))
(defmethod rng :well19937a [_ & [seed]]
  (create-object-with-seed Well19937a seed))
(defmethod rng :well19937c [_ & [seed]]
  (create-object-with-seed Well19937c seed))
(defmethod rng :well44497a [_ & [seed]]
  (create-object-with-seed Well44497a seed))
(defmethod rng :well44497b [_ & [seed]]
  (create-object-with-seed Well44497b seed))
(defmethod rng :jdk [_ & [seed]]
  (create-object-with-seed JDKRandomGenerator seed))
(defmethod rng :default [_ & [seed]]
  (rng :jdk seed))

(defn synced-rng
  "Create synchronized RNG for given name and optional seed. Wraps [[rng]] method."
  ([m] (SynchronizedRandomGenerator. (rng m)))
  ([m seed] (SynchronizedRandomGenerator. (rng m seed))))

;; List of randomizers
(def ^{:metadoc/categories #{:rand}
       :doc "List of all possible RNGs."}
  rngs-list (remove #{:default} (keys (methods rng))))

;; ### Default RNG

(defonce ^{:doc "Default RNG - JDK"
           :metadoc/categories #{:rand}}
  default-rng (rng :jdk))

(def ^{:doc "Random float number with Mersenne Twister RNG."
       :metadoc/categories #{:rand}}
  frand (partial frandom default-rng))

(def ^{:doc "Random boolean with Mersenne Twister RNG."
       :metadoc/categories #{:rand}} 
  brand (partial brandom default-rng))

(defn drand
  "Random double number with Mersenne Twister RNG."
  {:metadoc/categories #{:rand}}
  (^double [] (drandom default-rng))
  (^double [mx] (drandom default-rng mx))
  (^double [mn mx] (drandom default-rng mn mx)))

(defn grand
  "Random gaussian double number with Mersenne Twister RNG."
  {:metadoc/categories #{:rand}}
  (^double [] (grandom default-rng))
  (^double [stddev] (grandom default-rng stddev))
  (^double [mean stddev] (grandom default-rng mean stddev)))

(defn irand
  "Random integer number with Mersenne Twister RNG."
  {:metadoc/categories #{:rand}}
  (^long [] (irandom default-rng))
  (^long [mx] (irandom default-rng mx))
  (^long [mn mx] (irandom default-rng mn mx)))

(defn lrand
  "Random long number with Mersenne Twister RNG."
  {:metadoc/categories #{:rand}}
  (^long [] (lrandom default-rng))
  (^long [mx] (lrandom default-rng mx))
  (^long [mn mx] (lrandom default-rng mn mx)))

(defmacro randval
  "Retrun value with given probability (default 0.5)"
  {:metadoc/categories #{:rand}}
  ([v1 v2]
   `(if (brandom default-rng) ~v1 ~v2))
  ([prob v1 v2]
   `(if (brandom default-rng ~prob) ~v1 ~v2))
  ([prob]
   `(brandom default-rng ~prob))
  ([]
   `(brandom default-rng)))

(defn flip
  "Returns 1 with given probability, 0 otherwise"
  (^long [p]
   (randval p 1 0))
  (^long []
   (randval 0.5 1 0)))

(defn flipb
  "Returns true with given probability, false otherwise"
  ([p] (randval p))
  ([] (randval)))

;; generators

(defn- rv-generators
  "Generators from commons math and custom classes."
  [seq-generator ^long dimensions]
  (let [s (case seq-generator
            :halton (m/constrain dimensions 1 40)
            :sobol (m/constrain dimensions 1 1000)
            :r2 (m/constrain dimensions 1 4)
            dimensions)
        ^RandomVectorGenerator g (case seq-generator
                                   :halton (HaltonSequenceGenerator. s)
                                   :sobol (SobolSequenceGenerator. s)
                                   :sphere (UnitSphereRandomVectorGenerator. s)
                                   :r2 (R2. s))]
    (repeatedly (case s
                  1 #(aget (.nextVector g) 0)
                  2 #(v/array->vec2 (.nextVector g))
                  3 #(v/array->vec3 (.nextVector g))
                  4 #(v/array->vec4 (.nextVector g))
                  #(vec (.nextVector g))))))

;; R2
;; http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/

(defn- random-generators
  "Random generators"
  [seq-generator ^long dimensions]
  (let [g (case seq-generator
            :default drand
            :gaussian grand)]
    (repeatedly (case dimensions
                  1 g
                  2 (partial v/generate-vec2 g)
                  3 (partial v/generate-vec3 g)
                  4 (partial v/generate-vec4 g)
                  #(vec (repeatedly dimensions g))))))

;; jittering
;; http://extremelearning.com.au/a-simple-method-to-construct-isotropic-quasirandom-blue-noise-point-sequences/

(defn- jitter-generator
  "Generate random jitter"
  [seq-generator ^long dimensions ^double jitter]
  (let [[^double d0 ^double i0 ^double f ^double p] (case seq-generator
                                                      :r2 [0.76 0.7 0.25 -0.5]
                                                      :halton [0.9 0.7 0.25 -0.5]
                                                      :sobol [0.16 0.58 0.4 -0.2]
                                                      [0.5 0.5 0.25 -0.5])
        c (* jitter m/SQRTPI d0 f)
        g (random-generators :default dimensions)
        mod-fn (if (== dimensions 1)
                 (fn [^long i ^double v] (* v c (m/pow (- (inc i) i0) p)))
                 (fn [^long i v] (v/mult v (* c (m/pow (- (inc i) i0) p)))))]
    (map-indexed mod-fn g)))

;; Sequence creators

(defmulti
  ^{:doc "Create Sequence generator. See [[sequence-generators-list]] for names.

Values:

* `:r2`, `:halton`, `:sobol`, `:default` - range `[0-1] for each dimension`
* `:gaussian` - from `N(0,1)` distribution
* `:sphere` -  from surface of unit sphere (ie. euclidean distance from origin equals 1.0)

Possible dimensions:

* `:r2` - 1-4
* `:halton` - 1-40
* `:sobol` - 1-1000
* the rest - 1+

See also [[jittered-sequence-generator]]."
    :metadoc/categories #{:gen}}
  sequence-generator (fn [seq-generator _] seq-generator))
(defmethod sequence-generator :halton [seq-generator dimensions] (rv-generators seq-generator dimensions))
(defmethod sequence-generator :sobol [seq-generator dimensions] (rv-generators seq-generator dimensions))
(defmethod sequence-generator :r2 [seq-generator dimensions] (rv-generators seq-generator dimensions))
(defmethod sequence-generator :sphere [seq-generator dimensions] (rv-generators seq-generator dimensions))
(defmethod sequence-generator :gaussian [seq-generator dimensions] (random-generators seq-generator dimensions))
(defmethod sequence-generator :default [seq-generator dimensions] (random-generators seq-generator dimensions))

(defn jittered-sequence-generator
  "Create jittered sequence generator.

  Suitable for `:r2`, `:sobol` and `:halton` sequences.

  `jitter` parameter range is from `0` (no jitter) to `1` (full jitter).

  See also [[sequence-generator]]."
  [seq-generator ^long dimensions jitter]
  (let [s (sequence-generator seq-generator dimensions)
        j (jitter-generator seq-generator dimensions jitter)
        mod-fn (if (== dimensions 1)
                 (fn [^double v ^double vj] (m/frac (+ v vj)))
                 (fn [v vj] (v/fmap (v/add v vj) m/frac)))]
    (map mod-fn s j)))

(def ^{:doc "List of random sequence generator. See [[sequence-generator]]."
       :metadoc/categories #{:gen}}
  sequence-generators-list (keys (methods sequence-generator)))

;; ## Noise

(def ^{:doc "List of possible noise interpolations as a map of names and values."
       :metadoc/categories #{:noise}}
  interpolations {:none NoiseConfig/INTERPOLATE_NONE
                  :linear NoiseConfig/INTERPOLATE_LINEAR
                  :hermite NoiseConfig/INTERPOLATE_HERMITE
                  :quintic NoiseConfig/INTERPOLATE_QUINTIC})

(def ^{:doc "List of possible noise types as a map of names and values."
       :metadoc/categories #{:noise}}
  noise-types {:value NoiseConfig/NOISE_VALUE
               :gradient NoiseConfig/NOISE_GRADIENT
               :simplex NoiseConfig/NOISE_SIMPLEX})

(defn- noise-config-obj
  "Create noise configuration object based on map."
  [{:keys [seed noise-type interpolation octaves lacunarity gain normalize?]}]
  (NoiseConfig. seed
                (or (noise-types noise-type) NoiseConfig/NOISE_GRADIENT)
                (or (interpolations interpolation) NoiseConfig/INTERPOLATE_HERMITE)
                octaves lacunarity gain normalize?))

(defn- noise-config
  "Create FBM noise function for given configuration."
  ([] (noise-config {}))
  ([cfg]
   (noise-config-obj (merge {:seed (irand)
                             :noise-type :gradient
                             :interpolation :hermite
                             :octaves 6
                             :lacunarity 2.00
                             :gain 0.5
                             :normalize? true} cfg))))

(def ^:private perlin-noise-config (noise-config {:interpolation :quintic}))
(def ^:private simplex-noise-config (noise-config {:noise-type :simplex}))
(def ^:private value-noise-config (noise-config {:noise-type :value}))

(defn vnoise
  "Value Noise.

  6 octaves, Hermite interpolation (cubic, h01)."
  {:metadoc/categories #{:noise}}
  (^double [^double x] (FBM/noise value-noise-config x))
  (^double [^double x ^double y] (FBM/noise value-noise-config x y))
  (^double [^double x ^double y ^double z] (FBM/noise value-noise-config x y z)))

(defn noise
  "Create improved Perlin Noise.

  6 octaves, quintic interpolation."
  {:metadoc/categories #{:noise}}
  (^double [^double x] (FBM/noise perlin-noise-config x))
  (^double [^double x ^double y] (FBM/noise perlin-noise-config x y))
  (^double [^double x ^double y ^double z] (FBM/noise perlin-noise-config x y z)))

(defn simplex
  "Create Simplex noise. 6 octaves."
  {:metadoc/categories #{:noise}}
  (^double [^double x] (FBM/noise simplex-noise-config x))
  (^double [^double x ^double y] (FBM/noise simplex-noise-config x y))
  (^double [^double x ^double y ^double z] (FBM/noise simplex-noise-config x y z)))

(defmacro ^:private gen-noise-function
  "Generate various noise for static function"
  [method]
  `(fn gen-noise-name#
     ([] (gen-noise-name# {}))
     ([cfg#]
      (let [ncfg# (noise-config cfg#)]
        (fn
          (^double [x#] (~method ncfg# x#))
          (^double [x# y#] (~method ncfg# x# y#))
          (^double [x# y# z#] (~method ncfg# x# y# z#)))))))

(def ^{:metadoc/categories #{:noise}} single-noise (gen-noise-function Noise/noise))
(def ^{:metadoc/categories #{:noise}} fbm-noise (gen-noise-function FBM/noise))
(def ^{:metadoc/categories #{:noise}} billow-noise (gen-noise-function Billow/noise))
(def ^{:metadoc/categories #{:noise}} ridgedmulti-noise (gen-noise-function RidgedMulti/noise))

(defn random-noise-cfg
  "Create random noise configuration."
  {:metadoc/categories #{:noise}}
  []
  {:seed (irand)
   :noise-type (rand-nth (keys noise-types))
   :interpolation (rand-nth (keys interpolations))
   :octaves (irand 1 10)
   :lacunarity (drand 1.5 2.5)
   :gain (drand 0.2 0.8)
   :normalize? true})

(defn random-noise-fn
  "Create random noise function from all possible options.

  Optionally provide own configuration `cfg`. In this case one of 4 different blending methods will be selected."
  {:metadoc/categories #{:noise}}
  ([cfg]
   (let [cfg (merge (random-noise-cfg) cfg)]
     (rand-nth [(single-noise cfg)
                (fbm-noise cfg)
                (billow-noise cfg)
                (ridgedmulti-noise cfg)])))
  ([] (random-noise-fn (random-noise-cfg))))


;; ### Discrete noise

(defmacro discrete-noise
  "Discrete noise. Parameters:

  * X (long)
  * Y (long, optional)

  Returns double value from [0,1] range"
  {:metadoc/categories #{:noise}}
  ([X Y]
   `(Discrete/value ~X ~Y))
  ([X]
   `(discrete-noise ~X 0)))

;; Distribution

(defprotocol DistributionProto
  "Get information from distributions."
  (^{:metadoc/categories #{:dist}} cdf [d v] [d v1 v2] "Cumulative probability.")
  (^{:metadoc/categories #{:dist}} pdf [d v] "Density")
  (^{:metadoc/categories #{:dist}} lpdf [d v] "Log density")
  (^{:metadoc/categories #{:dist}} icdf [d p] "Inversed cumulative probability")
  (^{:metadoc/categories #{:dist}} probability [d v] "Probability (PMF)")
  (^{:metadoc/categories #{:dist}} sample [d] "Returns random sample.")
  (^{:metadoc/categories #{:dist}} dimensions [d] "Returns dimensions")
  (^{:metadoc/categories #{:dist}} source-object [d] "Returns Java object from backend library")
  (^{:metadoc/categories #{:dist}} continuous? [d] "true if distributions supports real values"))

(defprotocol UnivariateDistributionProto
  (^{:metadoc/categories #{:dist}} mean [d] "Mean")
  (^{:metadoc/categories #{:dist}} variance [d] "Variance")
  (^{:metadoc/categories #{:dist}} lower-bound [d] "Lower value")
  (^{:metadoc/categories #{:dist}} upper-bound [d] "Higher value"))

(defprotocol MultivariateDistributionProto
  "Get information from distributions."
  (^{:metadoc/categories #{:dist}} means [d] "Mean")
  (^{:metadoc/categories #{:dist}} covariance [d] "Variance"))

;; apache commons math
(extend RealDistribution
  DistributionProto
  {:cdf (fn
          ([^RealDistribution d ^double v] (.cumulativeProbability d v))
          ([^RealDistribution d ^double v1 ^double v2] (.cumulativeProbability d v1 v2)))
   :pdf (fn [^RealDistribution d ^double v] (.density d v))
   :lpdf (fn [^AbstractRealDistribution d ^double v] (.logDensity d v))
   :icdf (fn [^RealDistribution d ^double p] (.inverseCumulativeProbability d p))
   :probability (fn [^RealDistribution d ^double p] (.density d p))
   :sample (fn [^RealDistribution d] (.sample d))
   :dimensions (constantly 1)
   :source-object identity
   :continuous? (constantly true)} 
  UnivariateDistributionProto
  {:mean (fn [^RealDistribution d] (.getNumericalMean d))
   :variance (fn [^RealDistribution d] (.getNumericalVariance d))
   :lower-bound (fn [^RealDistribution d] (.getSupportLowerBound d))
   :upper-bound (fn [^RealDistribution d] (.getSupportUpperBound d))}
  RNGProto
  {:drandom (fn [^RealDistribution d] (.sample d))
   :frandom (fn [^RealDistribution d] (unchecked-float (.sample d)))
   :lrandom (fn [^RealDistribution d] (unchecked-long (.sample d)))
   :irandom (fn [^RealDistribution d] (unchecked-int (.sample d)))
   :->seq (fn
            ([^RealDistribution d] (repeatedly #(.sample d)))
            ([^RealDistribution d n] (repeatedly n #(.sample d))))
   :set-seed! (fn [^RealDistribution d ^double seed] (.reseedRandomGenerator d seed) d)})

;; ssj

(defn- reify-continuous-ssj
  [^ContinuousDistribution d ^RandomGenerator rng]
  (reify
    DistributionProto
    (pdf [_ v] (.density d v))
    (lpdf [_ v] (m/log (.density d v)))
    (cdf [_ v] (.cdf d v))
    (cdf [_ v1 v2] (- (.cdf d v2) (.cdf d v1)))
    (icdf [_ v] (.inverseF d v))
    (probability [_ v] (.density d v))
    (sample [_] (.inverseF d (drandom rng)))
    (dimensions [_] 1)
    (source-object [_] d)
    (continuous? [_] true)
    UnivariateDistributionProto
    (mean [_] (.getMean d))
    (variance [_] (.getVariance d))
    (lower-bound [_] (.getXinf d))
    (upper-bound [_] (.getXsup d))
    RNGProto
    (drandom [_] (.inverseF d (drandom rng)))
    (frandom [_] (unchecked-float (.inverseF d (drandom rng))))
    (lrandom [_] (unchecked-long (.inverseF d (drandom rng))))
    (irandom [_] (unchecked-int (.inverseF d (drandom rng))))
    (->seq [_] (repeatedly #(.inverseF d (drandom rng))))
    (->seq [_ n] (repeatedly n #(.inverseF d (drandom rng))))
    (set-seed! [d seed] (set-seed! rng seed) d)))

(defn- reify-integer-ssj
  [^DiscreteDistributionInt d ^RandomGenerator rng]
  (reify
    DistributionProto
    (pdf [_ v] (.prob d (m/floor v)))
    (lpdf [_ v] (m/log (.prob d (m/floor v))))
    (cdf [_ v] (.cdf d (m/floor v)))
    (cdf [_ v1 v2] (- (.cdf d (m/floor v2)) (.cdf d (m/floor v1))))
    (icdf [_ v] (.inverseF d v))
    (probability [_ v] (.prob d (m/floor v)))
    (sample [_] (.inverseF d (drandom rng)))
    (dimensions [_] 1)
    (source-object [_] d)
    (continuous? [_] false)
    UnivariateDistributionProto
    (mean [_] (.getMean d))
    (variance [_] (.getVariance d))
    (lower-bound [_] (.getXinf d))
    (upper-bound [_] (.getXsup d))
    RNGProto
    (drandom [_] (.inverseF d (drandom rng)))
    (frandom [_] (unchecked-float (.inverseF d (drandom rng))))
    (lrandom [_] (unchecked-long (.inverseF d (drandom rng))))
    (irandom [_] (unchecked-int (.inverseF d (drandom rng))))
    (->seq [_] (repeatedly #(.inverseF d (drandom rng))))
    (->seq [_ n] (repeatedly n #(.inverseF d (drandom rng))))
    (set-seed! [d seed] (set-seed! rng seed) d)))

;; smile
(extend DiscreteDistribution
  DistributionProto
  {:cdf (fn
          ([^Distribution d ^double v] (.cdf d (m/floor v)))
          ([^Distribution d ^double v1 ^double v2] (- (.cdf d (m/floor v2)) (.cdf d (m/floor v1)))))
   :pdf (fn [^Distribution d ^double v] (.p d (m/floor v)))
   :lpdf (fn [^Distribution d ^double v] (.logp d (m/floor v)))
   :icdf (fn [^Distribution d ^double p] (.quantile d p))
   :probability (fn [^Distribution d ^double v] (.p d (m/floor v)))
   :sample (fn [^Distribution d] (.rand d))
   :dimensions (constantly 1)
   :source-object identity
   :continuous? (constantly false)}  
  UnivariateDistributionProto
  {:mean (fn [^Distribution d] (.mean d))
   :variance (fn [^Distribution d] (.var d))}
  RNGProto
  {:drandom (fn [^Distribution d] (.rand d))
   :frandom (fn [^Distribution d] (unchecked-float (.rand d)))
   :lrandom (fn [^Distribution d] (unchecked-long (.rand d)))
   :irandom (fn [^Distribution d] (unchecked-int (.rand d)))
   :->seq (fn
            ([^Distribution d] (repeatedly #(.rand d)))
            ([^Distribution d n] (repeatedly n #(.rand d))))})

(extend IntegerDistribution
  DistributionProto
  {:cdf (fn
          ([^IntegerDistribution d ^double v] (.cumulativeProbability d (m/floor v)))
          ([^IntegerDistribution d ^double v1 ^double v2] (.cumulativeProbability d (m/floor v1) (m/floor v2))))
   :icdf (fn [^IntegerDistribution d ^double p] (.inverseCumulativeProbability d p))
   :pdf (fn [^IntegerDistribution d ^double p] (.probability d (m/floor p)))
   :lpdf (fn [^AbstractIntegerDistribution d ^double p] (.logProbability d (m/floor p)))
   :probability (fn [^IntegerDistribution d ^double p] (.probability d (m/floor p)))
   :sample (fn [^IntegerDistribution d] (.sample d))
   :dimensions (constantly 1)
   :source-object identity
   :continuous? (constantly false)}
  UnivariateDistributionProto
  {:mean (fn [^IntegerDistribution d] (.getNumericalMean d))
   :variance (fn [^IntegerDistribution d] (.getNumericalVariance d))
   :lower-bound (fn [^IntegerDistribution d] (.getSupportLowerBound d))
   :upper-bound (fn [^IntegerDistribution d] (.getSupportUpperBound d))}
  RNGProto
  {:drandom (fn [^IntegerDistribution d] (unchecked-double (.sample d)))
   :frandom (fn [^IntegerDistribution d] (unchecked-float (.sample d)))
   :lrandom (fn [^IntegerDistribution d] (unchecked-long (.sample d)))
   :irandom (fn [^IntegerDistribution d] (.sample d))
   :->seq (fn
            ([^IntegerDistribution d] (repeatedly #(.sample d)))
            ([^IntegerDistribution d n] (repeatedly n #(.sample d))))
   :set-seed! (fn [^IntegerDistribution d ^double seed] (.reseedRandomGenerator d seed) d)})

(extend MultivariateNormalDistribution
  DistributionProto
  {:pdf (fn [^MultivariateNormalDistribution d v] (.density d (m/seq->double-array v)))
   :lpdf (fn [^MultivariateNormalDistribution d v] (m/log (.density d (m/seq->double-array v))))
   :sample (fn [^MultivariateNormalDistribution d] (vec (.sample d)))
   :dimensions (fn [^MultivariateNormalDistribution d] (.getDimension d))
   :source-object identity
   :continuous? (constantly true)}
  MultivariateDistributionProto
  {:means (fn [^MultivariateNormalDistribution d] (vec (.getMeans d)))
   :covariance (fn [^MultivariateNormalDistribution d]
                 (let [^org.apache.commons.math3.linear.Array2DRowRealMatrix cv (.getCovariances d)]
                   (m/double-double-array->seq (.getDataRef cv))))}
  RNGProto
  {:drandom (fn [^MultivariateNormalDistribution d] (vec (.sample d)))
   :frandom (fn [^MultivariateNormalDistribution d] (mapv unchecked-float (.sample d)))
   :lrandom (fn [^MultivariateNormalDistribution d] (mapv unchecked-long (.sample d)))
   :irandom (fn [^MultivariateNormalDistribution d] (mapv unchecked-int (.sample d)))
   :->seq (fn
            ([^MultivariateNormalDistribution d] (repeatedly #(vec (.sample d))))
            ([^MultivariateNormalDistribution d n] (repeatedly n #(vec (.sample d)))))
   :set-seed! (fn [^MultivariateNormalDistribution d ^double seed] (.reseedRandomGenerator d seed) d)})

(defn observe
  "Log of probability/density of the value. Alias for [[lpdf]]."
  ^double [d ^double v]
  (lpdf d v))

(defn log-likelihood
  "Log likelihood of samples"
  ^double [d vs] 
  (reduce #(if (m/inf? ^double %1)
             (reduced %1)
             (+ ^double %1 ^double %2)) 0.0 (map #(lpdf d %) vs)))

(defn likelihood
  "Likelihood of samples"
  ^double [d vs]
  (m/exp (log-likelihood d vs)))

(defmulti
  ^{:doc "Create distribution object.

First parameter is distribution as a `:key`.
Second parameter is a map with configuration.
All distributions accept `rng` under `:rng` key (default: [[default-rng]]) and some of them accept `inverse-cumm-accuracy` (default set to `1e-9`).

Distributions should be called using [[DistributionProto]] and [[RNGProto]].

The rest parameters goes as follows:

#### Real distributions

* `:beta` - `:alpha` (default: 2.0) and `:beta` (default: 5.0)
* `:cauchy` - `:median` (default: 0.0) and `:scale` (default: 1.0)
* `:chi-squared` - `:degrees-of-freedom` (default: 1.0)
* `:empirical` - `:bean-count` (default: 1000) and `:data` as a sequence
* `:enumerated-real` - `:data` as a sequence and `:probabilities` as a optional sequence
* `:exponential` - `:mean` (default: 1.0)
* `:f` - `:numerator-degrees-of-freedom` (default: 1.0) and `:denominator-degrees-of-freedom` (default: 1.0)
* `:gamma` - `:shape` (default: 2.0) and `:scale` (default: 2.0)
* `:gumbel` - `:mu` (default: 1.0) and `:beta` (default: 2.0)
* `:laplace` - `:mu` (default: 1.0) and `:beta` (default: 2.0)
* `:levy` - `:mu` (default: 0.0) and `:c` (default: 1.0)
* `:logistic` - `:mu` (default: 1.0) and `:s` (default: 2.0)
* `:log-normal` - `:scale` (default: 1.0) and `:shape` (default: 1.0)
* `:nakagami` - `:mu` (default: 1.0) and `:omega` (default: 1.0)
* `:normal` - `:mu` (default: 0.0) and `:sd` (default: 1.0)
* `:pareto` - `:scale` (default: 1.0) and `:shape` (default: 1.0)
* `:t` - `:degrees-of-freedom` (default: 1.0)
* `:triangular` - `:a` (default: -1.0), `:c` (default: 0.0) and `:b` (default: 1.0)
* `:uniform-real` - `:lower` (default: 0.0) and `:upper` (default: 1.0)
* `:weibull` - `:alpha` (default: 2.0) and `:beta` (default: 1.0)

#### Integer distributions

* `:binomial` - `:trials` (default: 20) and `:p` (default: 0.5)
* `:bernoulli` - `:p` (default: 0.5) 
* `:enumerated-int` - `:data` and `:probabilities` as a sequences
* `:geometric` - `:p` (default: 0.5)
* `:hypergeometric` - `:population-size` (default: 100), `:number-of-successes` (default: 50) and `:sample-size` (default: 25)
* `:pascal` - `:r` (default: 5) and `:p` (default: 0.5)
* `:poisson` - `:p` (default: 0.5), `:epsilon` (default: 1.0e-12), `:max-iterations` (default: 10000000)
* `:uniform-int` - `:lower` (default: 0) and `:upper` (default: `Integer/MAX_VALUE`)
* `:zipf` - `:number-of-elements` (default: 100) and `:exponent` (default: 3.0)
"
    :metadoc/categories #{:dist}}
  distribution (fn ([k _] k) ([k] k)))

(defmacro ^:private make-acm-distr
  [nm obj ks vs]
  (let [or-map (zipmap ks vs)] 
    `(defmethod distribution ~nm
       ([n# {:keys [~@ks]
             :or ~or-map
             :as all#}]
        (let [^RandomGenerator r# (or (:rng all#) (rng :jvm))]
          (new ~obj r# ~@ks)))
       ([n#] (distribution ~nm {})))))

(make-acm-distr :beta BetaDistribution
                [alpha beta inverse-cumm-accuracy]
                [2.0 5.0 BetaDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY])

(make-acm-distr :cauchy CauchyDistribution
                [median scale inverse-cumm-accuracy]
                [0.0 1.0 CauchyDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY])

(make-acm-distr :chi-squared ChiSquaredDistribution
                [degrees-of-freedom inverse-cumm-accuracy]
                [1.0 ChiSquaredDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY])

(make-acm-distr :exponential ExponentialDistribution
                [mean inverse-cumm-accuracy]
                [1.0 ExponentialDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY])

(make-acm-distr :f FDistribution
                [numerator-degrees-of-freedom denominator-degrees-of-freedom inverser-cumm-accuracy]
                [1.0 1.0 FDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY])

(make-acm-distr :gamma GammaDistribution
                [shape scale inverse-cumm-accuracy]
                [2.0 2.0 GammaDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY])

(make-acm-distr :gumbel GumbelDistribution [mu beta] [1.0 2.0])
(make-acm-distr :laplace LaplaceDistribution [mu beta] [0.0 1.0])
(make-acm-distr :levy LevyDistribution [mu c] [0.0 1.0])
(make-acm-distr :logistic LogisticDistribution [mu s] [0.0 1.0])

(make-acm-distr :log-normal LogNormalDistribution
                [scale shape inverse-cumm-accuracy]
                [1.0 1.0 LogNormalDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY])

(make-acm-distr :nakagami NakagamiDistribution
                [mu omega inverse-cumm-accuracy]
                [1.0 1.0 NakagamiDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY])

(make-acm-distr :normal NormalDistribution
                [mu sd inverse-cumm-accuracy]
                [0.0 1.0 NormalDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY])

(make-acm-distr :pareto ParetoDistribution
                [scale shape inverse-cumm-accuracy]
                [1.0 1.0 ParetoDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY])

(make-acm-distr :t TDistribution
                [degrees-of-freedom inverser-cumm-accuracy]
                [1.0 TDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY])

(make-acm-distr :triangular TriangularDistribution [a c b] [-1.0 0.0 1.0])
(make-acm-distr :uniform-real UniformRealDistribution [^double lower ^double upper] [0.0 1.0])

(make-acm-distr :weibull WeibullDistribution
                [alpha beta inverse-cumm-accuracy]
                [2.0 1.0 WeibullDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY])

(defmethod distribution :empirical
  ([_ {:keys [^long bin-count data]
       :or {bin-count EmpiricalDistribution/DEFAULT_BIN_COUNT}
       :as all}]
   (let [^RandomGenerator r (or (:rng all) (rng :jvm))
         ^EmpiricalDistribution d (EmpiricalDistribution. bin-count r)]
     (.load d ^doubles (m/seq->double-array data))
     d))
  ([_] (distribution :empirical {})))

(defmethod distribution :enumerated-real
  ([_ {:keys [data probabilities]
       :as all}]
   (let [^RandomGenerator r (or (:rng all) (rng :jvm))]
     (if probabilities
       (EnumeratedRealDistribution. r (m/seq->double-array data) (m/seq->double-array probabilities))
       (EnumeratedRealDistribution. r ^doubles (m/seq->double-array data)))))
  ([_] (distribution :enumerated-real {})))

;; integer

(defmethod distribution :negative-binomial
  ([_ {:keys [^double r ^double p]
       :or {r 20.0 p 0.5}}]
   (NegativeBinomialDistribution. r p))
  ([_] (distribution :negative-binomial {})))

(defmethod distribution :bernoulli
  ([_ {:keys [^double p]
       :or {p 0.5}
       :as all}]
   (BinomialDistribution. (or (:rng all) (rng :jvm)) 1 p))
  ([_] (distribution :bernoulli {})))

(defmethod distribution :enumerated-int
  ([_ {:keys [data probabilities]
       :as all}]
   (let [^RandomGenerator r (or (:rng all) (rng :jvm))]
     (if probabilities
       (EnumeratedIntegerDistribution. r (int-array data) (m/seq->double-array probabilities))
       (EnumeratedIntegerDistribution. r (int-array data)))))
  ([_] (distribution :enumerated-int {})))

(make-acm-distr :binomial BinomialDistribution [trials p] [20 0.5])
(make-acm-distr :geometric GeometricDistribution [p] [0.5])
(make-acm-distr :hypergeometric HypergeometricDistribution
                [population-size number-of-successes sample-size] [100 50 25])
(make-acm-distr :pascal PascalDistribution [r p] [5 0.5])
(make-acm-distr :poisson PoissonDistribution
                [p epsilon max-iterations]
                [0.5 PoissonDistribution/DEFAULT_EPSILON PoissonDistribution/DEFAULT_MAX_ITERATIONS])
(make-acm-distr :uniform-int UniformIntegerDistribution [lower upper] [0 Integer/MAX_VALUE])
(make-acm-distr :zipf ZipfDistribution [number-of-elements exponent] [100 3.0])

;; ssj

(defmacro ^:private make-ssj-distr
  [rf nm obj ks vs]
  (let [or-map (zipmap ks vs)]
    `(defmethod distribution ~nm
       ([n# {:keys [~@ks]
             :or ~or-map
             :as all#}]
        (let [^RandomGenerator r# (or (:rng all#) (rng :jvm))]
          (~rf (new ~obj ~@ks) r#)))
       ([n#] (distribution ~nm {})))))

(defmacro ^:private make-ssjc-distr
  [nm obj ks vs] `(make-ssj-distr reify-continuous-ssj ~nm ~obj ~ks ~vs))
(defmacro ^:private make-ssji-distr
  [nm obj ks vs] `(make-ssj-distr reify-integer-ssj ~nm ~obj ~ks ~vs))


(make-ssjc-distr :anderson-darling AndersonDarlingDistQuick [n] [1.0])
(make-ssjc-distr :inverse-gamma InverseGammaDist [alpha beta] [2.0 1.0])
(make-ssjc-distr :chi ChiDist [nu] [1.0])
(make-ssjc-distr :chi-squared-noncentral ChiSquareNoncentralDist [nu lambda] [1.0 1.0])
(make-ssjc-distr :cramer-von-mises CramerVonMisesDist [n] [1.0])
(make-ssjc-distr :erlang ErlangDist [k lambda] [1 1])
(make-ssjc-distr :fatigue-life FatigueLifeDist [mu beta gamma] [0.0 1.0 1.0])
(make-ssjc-distr :folded-normal FoldedNormalDist [mu sigma] [0.0 1.0])
(make-ssjc-distr :frechet FrechetDist [alpha beta delta] [1.0 1.0 0.0])
(make-ssjc-distr :hyperbolic-secant HyperbolicSecantDist [mu sigma] [0.0 1.0])
(make-ssjc-distr :inverse-gaussian InverseGaussianDist [mu lambda] [1.0 1.0])
(make-ssjc-distr :hypoexponential-equal HypoExponentialDistEqual [n k h] [1.0 1.0 1.0])
(make-ssjc-distr :johnson-sb JohnsonSBDist [gamma delta xi lambda] [0.0 1.0 0.0 1.0])
(make-ssjc-distr :johnson-sl JohnsonSLDist [gamma delta xi lambda] [0.0 1.0 0.0 1.0])
(make-ssjc-distr :johnson-su JohnsonSUDist [gamma delta xi lambda] [0.0 1.0 0.0 1.0])
(make-ssjc-distr :kolmogorov-smirnov KolmogorovSmirnovDistQuick [n] [1.0])
(make-ssjc-distr :kolmogorov-smirnov+ KolmogorovSmirnovPlusDist [n] [1.0])
(make-ssji-distr :logarithmic LogarithmicDist [theta] [0.5])
(make-ssjc-distr :log-logistic LoglogisticDist [alpha beta] [3.0 1.0])
(make-ssjc-distr :normal-inverse-gaussian NormalInverseGaussianDist [alpha beta mu delta] [1.0 0.0 0.0 1.0])
(make-ssjc-distr :pearson-6 Pearson6Dist [alpha1 alpha2 beta] [1.0 1.0 1.0])
(make-ssjc-distr :power PowerDist [a b c] [0.0 1.0 2.0])
(make-ssjc-distr :rayleigh RayleighDist [a beta] [0.0 1.0])
(make-ssjc-distr :watson-g WatsonGDist [n] [2.0])
(make-ssjc-distr :watson-u WatsonUDist [n] [2.0])

(defmethod distribution :hypoexponential
  ([_ {:keys [lambdas]
       :or {lambdas [1.0]}
       :as all}]
   (reify-continuous-ssj (HypoExponentialDist. (m/seq->double-array lambdas)) (or (:rng all) (rng :jvm))))
  ([_] (distribution :hypoexponential {})))

(defmethod distribution :reciprocal-sqrt
  ([_ {:keys [^double a]
       :or {a 0.5}
       :as all}]
   (let [f (* 2.0 (m/sqrt a))
         icdf-fn (fn [^double x]
                   (cond
                     (zero? x) a
                     :else (m/sq (* 0.5 (+ x f)))))
         ^double b (icdf-fn 1.0)
         m (* (/ 2.0 3.0) (- (m/pow b 1.5) (m/pow a 1.5)))
         m1 (* 15.0 m m)
         m2 (* -10.0 m)
         v (* (/ 2.0 15.0)
              (- (* (m/sqrt b)
                    (+ m1 (* m2 b) (* 3.0 b b)))
                 (* (m/sqrt a)
                    (+ m1 (* m2 a) (* 3.0 a a)))))
         r (or (:rng all) (rng :jvm))]
     (reify
       DistributionProto
       (pdf [_ v] (if (<= a ^double v b) (/ (m/sqrt v)) 0.0))
       (lpdf [d v] (m/log (pdf d v)))
       (cdf [_ v] (cond
                    (< ^double v a) 0.0
                    (> ^double v b) 1.0
                    :else (- (* 2.0 (m/sqrt ^double v)) f)))
       (cdf [d v1 v2] (- ^double (cdf d v2) ^double (cdf d v1)))
       (icdf [_ v] (icdf-fn v))
       (probability [d v] (pdf d v))
       (sample [d] (icdf-fn (drandom r)))
       (dimensions [_] 1)
       (source-object [d] d)
       (continuous? [_] true)
       UnivariateDistributionProto
       (mean [_] m)
       (variance [_] v)
       (lower-bound [_] a)
       (upper-bound [_] b)
       RNGProto
       (drandom [_] (icdf-fn (drandom r)))
       (frandom [_] (unchecked-float (icdf-fn (drandom r))))
       (lrandom [_] (unchecked-long (icdf-fn (drandom r))))
       (irandom [_] (unchecked-int (icdf-fn (drandom r))))
       (->seq [_] (repeatedly #(icdf-fn (drandom r))))
       (->seq [_ n] (repeatedly n #(icdf-fn (drandom r))))
       (set-seed! [d seed] (set-seed! r seed) d)))))

;;

(defmethod distribution :multi-normal
  ([_ {:keys [means covariances] :as all}]
   (let [covariances (if-not covariances
                       [[1.0 0.0] [0.0 1.0]]
                       covariances)
         means (if-not means (repeat (count (first covariances)) 0.0) means)]
     (assert (= (count means) (count (first covariances)))
             "Means and covariances sizes do not match.")
     (MultivariateNormalDistribution. (or (:rng all) (rng :jvm)) (m/seq->double-array means) (m/seq->double-double-array covariances))))
  ([_] (distribution :multi-normal {})))

(defonce ^:const ^:private ^double zero+epsilon (m/next-double 0.0))
(defonce ^:const ^:private ^double one-epsilon (m/prev-double 1.0))

(defn- dirichlet-rev-log-beta
  ^double [alpha]
  (let [d (m/log-gamma (reduce m/fast+ alpha))
        ^double n (reduce m/fast+ (map #(m/log-gamma %) alpha))]
    (- d n)))

(defn- dirichlet-lpdf
  ^double [alpha- values ^double lbeta]
  (let [^double p (reduce m/fast+ (mapv (fn [^double ai ^double x]
                                          (* ai (m/log x))) alpha- values))]
    (+ lbeta p)))

(defmethod distribution :dirichlet
  ([_ {:keys [alpha] :as all}]
   (let [alpha (cond
                 (nil? alpha) (double-array [1])
                 (seqable? alpha) (m/seq->double-array alpha)
                 (integer? alpha) (double-array alpha 1.0)
                 :else (double-array [1]))
         r (or (:rng all) (rng :jvm))
         sampler (mapv #(distribution :gamma {:shape % :scale 1.0 :rng r}) alpha)

         lbeta (dirichlet-rev-log-beta alpha)
         alpha- (map clojure.core/dec alpha)
         
         m (delay (seq (DirichletDist/getMean alpha)))
         cv (delay (mapv vec (DirichletDist/getCovariance alpha)))
         dim (count alpha)]
     (reify
       DistributionProto
       (pdf [_ v] (m/exp (dirichlet-lpdf alpha- v lbeta)))
       (lpdf [_ v] (dirichlet-lpdf alpha- v lbeta))
       (probability [d v] (m/exp (dirichlet-lpdf alpha- v lbeta)))
       (sample [_] (let [samples (mapv #(sample %) sampler)
                         ^double s (v/sum samples)]
                     (mapv (fn [^double s] (cond
                                            (zero? s) zero+epsilon
                                            (== s 1.0) one-epsilon
                                            :else s))
                           (if (> s 1.0e-30)
                             (v/div samples (v/sum samples))
                             (let [a (int-array dim)]
                               (aset ^ints a (irand dim) 1)
                               a)))))
       (dimensions [_] dim)
       (source-object [this] this)
       (continuous? [_] true)
       MultivariateDistributionProto
       (means [_] @m)
       (covariance [_] @cv)
       RNGProto
       RNGProto
       (drandom [d] (sample d))
       (frandom [d] (mapv unchecked-float (sample d)))
       (lrandom [d] (mapv unchecked-long (sample d)))
       (irandom [d] (mapv unchecked-int (sample d)))
       (->seq [d] (repeatedly #(sample d)))
       (->seq [d n] (repeatedly n #(sample d)))
       (set-seed! [d seed] (set-seed! r seed) d)))) 
  ([_] (distribution :dirichlet {})))

;; 

(defmethod distribution :continuous-distribution
  ([_ {:keys [data kernel h bin-count probabilities]
       :or {kernel :smile}
       :as all}]
   (let [d (m/seq->double-array data)]
     (java.util.Arrays/sort d)
     (let [r (or (:rng all) (rng :jvm)) 
           kde (if h
                 (k/kernel-density kernel d h)
                 (k/kernel-density kernel d))
           ^RealDistribution enumerated (distribution :enumerated-real {:data d :probabilities probabilities :rng r})
           ^RealDistribution empirical (distribution :empirical (if-not bin-count
                                                                  {:data d}
                                                                  {:data d :bin-count bin-count}))]
       (reify
         DistributionProto
         (pdf [_ v] (kde v))
         (lpdf [_ v] (m/log (kde v)))
         (cdf [_ v] (.cumulativeProbability enumerated ^double v))
         (cdf [_ v1 v2] (.cumulativeProbability enumerated ^double v1 ^double v2))
         (icdf [_ v] (.inverseCumulativeProbability empirical v))
         (probability [_ v] (kde v))
         (sample [_] (.sample enumerated))
         (dimensions [_] 1)
         (source-object [d] {:enumerated enumerated
                             :empirical empirical})
         (continuous? [_] true)
         UnivariateDistributionProto
         (mean [_] (mean enumerated))
         (variance [_] (variance enumerated))
         (lower-bound [_] (lower-bound enumerated))
         (upper-bound [_] (upper-bound enumerated))
         RNGProto
         (drandom [_] (drandom enumerated))
         (frandom [_] (frandom enumerated))
         (lrandom [_] (lrandom enumerated))
         (irandom [_] (irandom enumerated))
         (->seq [_] (->seq enumerated))
         (->seq [_ n] (->seq enumerated n))
         (set-seed! [d seed] (set-seed! r seed) d))))))

(defmethod distribution :integer-discrete-distribution [_ d]
  (distribution :enumerated-int (update d :data #(map (fn [^double v]
                                                        (int (m/floor v))) %))))

(defmethod distribution :real-discrete-distribution [_ d]
  (distribution :enumerated-real d))

(defmethod distribution :categorical-distribution
  ([_ {:keys [data probabilities]
       :as all}]
   (let [r (or (:rng all) (rng :jvm))
         
         ^clojure.lang.ILookup unique (vec (distinct data))
         ^clojure.lang.ILookup dict (zipmap unique (range (count unique)))
         
         ^AbstractIntegerDistribution enumerated (distribution :enumerated-int {:data (map dict data) :probabilities probabilities :rng r})]
     (reify
       DistributionProto
       (pdf [_ v] (.probability enumerated (.valAt dict v -1)))
       (lpdf [_ v] (.logProbability enumerated (.valAt dict v -1)))
       (cdf [_ v] (.cumulativeProbability enumerated (.valAt dict v -1)))
       (icdf [_ v] (.valAt unique (.inverseCumulativeProbability enumerated ^double v)))
       (probability [_ v] (.probability enumerated (.valAt dict v -1)))
       (sample [_] (.valAt unique (.sample enumerated)))
       (dimensions [_] 1)
       (source-object [d] enumerated)
       (continuous? [_] false)
       UnivariateDistributionProto
       (mean [_] (mean enumerated))
       (variance [_] (variance enumerated))
       RNGProto
       (->seq [_] (map #(.valAt unique %) (->seq enumerated)))
       (->seq [_ n] (map #(.valAt unique %) (->seq enumerated n)))
       (set-seed! [d seed] (set-seed! r seed) d)))))

;;

(def ^{:doc "List of distributions."
       :metadoc/categories #{:dist}}
  distributions-list
  (into (sorted-set) (keys (methods distribution))))

(defonce ^{:doc "Default normal distribution (u=0.0, sigma=1.0)."
           :metadoc/categories #{:dist}}
  default-normal (distribution :normal))

(type (zipmap [:a] [:b]))
