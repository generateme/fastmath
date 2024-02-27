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

  To create generator call [[sequence-generator]] with generator name and vector size.
  Following generators are available:

  * `:halton` - Halton low-discrepancy sequence; range [0,1]
  * `:sobol` - Sobol low-discrepancy sequence; range [0,1]
  * `:r2` - R2 low-discrepancy sequence; range [0,1], [more...](http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/)
  * `:sphere` - uniformly random distributed on unit sphere
  * `:ball` - uniformly random distributed from unit ball
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
  
  All can be combined in following variants:

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

  For usage convenience 3 ready to use functions are prepared. Returning value from `[0,1]` range:

  * [[noise]] - Perlin Noise (gradient noise, 6 octaves, quintic interpolation)
  * [[vnoise]] - Value Noise (as in Processing, 6 octaves, hermite interpolation)
  * [[simplex]] - Simplex Noise (6 octaves)

  For random noise generation you can use [[random-noise-cfg]] and [[random-noise-fn]]. Both can be feed with configuration. Additional configuration:

  * `:generator` can be set to one of the noise variants, defaults to `:fbm`
  * `:warp-scale` - 0.0 - do not warp, >0.0 warp
  * `:warp-depth` - depth for warp (default 1.0, if warp-scale is positive)
  
  #### Discrete Noise

  [[discrete-noise]] is a 1d or 2d hash function for given integers. Returns double from `[0,1]` range.

  ### Distribution

  Various real and integer distributions. See [[DistributionProto]] and [[RNGProto]] for functions.

  To create distribution call [[distribution]] multimethod with name as a keyword and map as parameters."  
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.kernel :as k]
            [fastmath.protocols :as prot]
            [fastmath.interpolation :as i]
            [fastmath.solver :as solver]
            [clojure.data.int-map :as im]            )
  (:import [org.apache.commons.math3.random RandomGenerator ISAACRandom JDKRandomGenerator MersenneTwister
            Well512a Well1024a Well19937a Well19937c Well44497a Well44497b
            RandomVectorGenerator HaltonSequenceGenerator SobolSequenceGenerator UnitSphereRandomVectorGenerator
            EmpiricalDistribution SynchronizedRandomGenerator]
           [fastmath.java R2]
           [umontreal.ssj.probdist ContinuousDistribution DiscreteDistributionInt InverseGammaDist AndersonDarlingDistQuick ChiDist ChiSquareNoncentralDist CramerVonMisesDist ErlangDist FatigueLifeDist FoldedNormalDist FrechetDist HyperbolicSecantDist InverseGaussianDist HypoExponentialDist HypoExponentialDistEqual JohnsonSBDist JohnsonSLDist JohnsonSUDist KolmogorovSmirnovDistQuick KolmogorovSmirnovPlusDist LogarithmicDist LoglogisticDist NormalInverseGaussianDist Pearson6Dist PowerDist RayleighDist WatsonGDist WatsonUDist]
           [umontreal.ssj.probdistmulti DirichletDist MultinomialDist]
           [fastmath.java.noise Billow RidgedMulti FBM NoiseConfig Noise Discrete]
           [org.apache.commons.math3.stat StatUtils]
           [org.apache.commons.math3.distribution AbstractRealDistribution RealDistribution BetaDistribution CauchyDistribution ChiSquaredDistribution EnumeratedRealDistribution ExponentialDistribution FDistribution GammaDistribution, GumbelDistribution, LaplaceDistribution, LevyDistribution, LogisticDistribution, LogNormalDistribution, NakagamiDistribution, NormalDistribution, ParetoDistribution, TDistribution, TriangularDistribution, UniformRealDistribution WeibullDistribution MultivariateNormalDistribution]
           [org.apache.commons.math3.distribution IntegerDistribution AbstractIntegerDistribution BinomialDistribution EnumeratedIntegerDistribution, GeometricDistribution, HypergeometricDistribution, PascalDistribution, PoissonDistribution, UniformIntegerDistribution, ZipfDistribution]
           [org.apache.commons.math3.analysis UnivariateFunction]
           [org.apache.commons.math3.analysis.integration RombergIntegrator]
           [smile.math MathEx]))

(set! *unchecked-math* :warn-on-boxed)
(m/use-primitive-operators)

;; Helper macro which creates RNG object of given class and/or seed.
(defmacro ^:private create-object-with-seed
  "Create object of the class with (or not) given seed. Used to create RNG."
  [cl seed]
  `(if-let [arg# ~seed]
     (prot/set-seed! (new ~cl) (long arg#)) ;; seeding via protocols
     (new ~cl)))

(defmulti rng
  "Create RNG for given name (as keyword) and optional seed. Return object enhanced with [[RNGProto]]. See: [[rngs-list]] for names."
  (fn [m & _] m))

(def ^:private rng-class->keyword {MersenneTwister :mersenne
                                 ISAACRandom :isaac
                                 Well512a :well512a
                                 Well1024a :well1024a
                                 Well19937a :well19937a
                                 Well19937c :well19937c
                                 Well44497a :well44497a
                                 Well44497b :well44497b
                                 JDKRandomGenerator :jdk})

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
(defonce ^{:doc "List of all possible RNGs."}
  rngs-list (remove #{:default} (keys (methods rng))))

;; protocol proxies
(defn frandom
  "Random double number with provided RNG"
  ([rng] (prot/frandom rng))
  ([rng mx] (prot/frandom rng mx))
  ([rng mn mx] (prot/frandom rng mn mx)))

(defn drandom
  "Random double number with provided RNG"
  (^double [rng] (prot/drandom rng))
  (^double [rng mx] (prot/drandom rng mx))
  (^double [rng mn mx] (prot/drandom rng mn mx)))

(defn grandom
  "Random gaussian double number with provided RNG"
  (^double [rng] (prot/grandom rng))
  (^double [rng stddev] (prot/grandom rng stddev))
  (^double [rng mean stddev] (prot/grandom rng mean stddev)))

(defn irandom
  "Random integer number with provided RNG"
  (^long [rng] (prot/irandom rng))
  (^long [rng mx] (prot/irandom rng mx))
  (^long [rng mn ^long mx] (prot/irandom rng mn mx)))

(defn lrandom
  "Random long number with provided RNG"
  (^long [rng] (prot/lrandom rng))
  (^long [rng mx] (prot/lrandom rng mx))
  (^long [rng mn mx] (prot/lrandom rng mn mx)))

(defn brandom
  "Random boolean with provided RNG"
  ([rng] (prot/brandom rng))
  ([rng p] (prot/brandom rng p)))

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

(defn- next-random-value-int
  "Generate next int.

  * arity 0 - from 0 to maximum int value
  * arity 1 - from 0 to provided integer (excluded)
  * arity 2 - from the provided range (included, excluded)"
  (^long [^RandomGenerator r] (.nextInt r))
  (^long [^RandomGenerator r ^long mx] (.nextInt r mx))
  (^long [r ^long mn ^long mx]
   (let [diff (- mx mn)]
     (if (zero? diff) mn
         (+ mn (next-random-value-int r diff))))))

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

(defn- next-random-value-float
  "Generate next float.

  * arity 0 - from 0 to 1 (exluded)
  * arity 1 - from 0 to provided float (excluded)
  * arity 2 - from the provided range (included, excluded)"
  ([^RandomGenerator r] (.nextFloat r))
  ([^RandomGenerator r ^double mx] (unchecked-float (* (.nextFloat r) mx)))
  ([r ^double mn ^double mx]
   (let [diff (- mx mn)]
     (unchecked-float (if (zero? diff) mn
                          (+ mn ^float (next-random-value-float r diff)))))))

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

;; Extend RandomGenerator interface with functions created by macro `next-random-value-fn`. This way all RNG classes are enriched with new, more convenient functions.
;;
;; Note that `grandom` is under special care due to different [mn mx] range meaning.

(extend RandomGenerator 
  prot/RNGProto
  {:irandom next-random-value-int
   :lrandom next-random-value-long
   :frandom next-random-value-float
   :drandom next-random-value-double
   :grandom (fn
              ([t] (next-random-value-gaussian t))
              ([t std] (next-random-value-gaussian t std))
              ([t ^double mean ^double std] (next-random-value-gaussian t mean (+ mean std))))
   :brandom (fn
              ([^RandomGenerator t] (.nextBoolean t))
              ([t ^double thr] (< (next-random-value-double t) thr)))
   :set-seed! (fn [^RandomGenerator t ^long seed]
                (.setSeed t seed)
                t)
   :set-seed #(let [rng-name (rng-class->keyword (class %1))]
                (rng rng-name (long %2)))
   :->seq (fn
            ([^RandomGenerator t] (repeatedly #(next-random-value-double t)))
            ([^RandomGenerator t n] (repeatedly n #(next-random-value-double t))))})

;; ### Default RNG

(defonce ^{:doc "Default RNG - JDK"} default-rng (rng :jdk))

(def ^{:doc "Random boolean with default RNG.

Returns true or false with equal probability. You can set `p` probability for `true`"} 
  brand (partial prot/brandom default-rng))

(defn frand
  "Random double number with default RNG.

  As default returns random float from `[0,1)` range.
  When `mx` is passed, range is set to `[0, mx)`. When `mn` is passed, range is set to `[mn, mx)`."
  ([] (prot/frandom default-rng))
  ([mx] (prot/frandom default-rng mx))
  ([mn mx] (prot/frandom default-rng mn mx)))

(defn drand
  "Random double number with default RNG.

  As default returns random double from `[0,1)` range.
  When `mx` is passed, range is set to `[0, mx)`. When `mn` is passed, range is set to `[mn, mx)`."
  (^double [] (prot/drandom default-rng))
  (^double [^double mx] (prot/drandom default-rng mx))
  (^double [^double mn ^double mx] (prot/drandom default-rng mn mx)))

(defn grand
  "Random gaussian double number with default RNG.

  As default returns random double from `N(0,1)`.
  When `std` is passed, `N(0,std)` is used. When `mean` is passed, distribution is set to `N(mean, std)`."
  (^double [] (prot/grandom default-rng))
  (^double [^double stddev] (prot/grandom default-rng stddev))
  (^double [^double mean ^double stddev] (prot/grandom default-rng mean stddev)))

(defn irand
  "Random integer number with default RNG.

  As default returns random integer from full integer range. 
  When `mx` is passed, range is set to `[0, mx)`. When `mn` is passed, range is set to `[mn, mx)`."
  (^long [] (prot/irandom default-rng))
  (^long [mx] (prot/irandom default-rng mx))
  (^long [mn mx] (prot/irandom default-rng mn mx)))

(defn lrand
  "Random long number with default RNG.

  As default returns random long from full integer range. 
  When `mx` is passed, range is set to `[0, mx)`. When `mn` is passed, range is set to `[mn, mx)`."
  (^long [] (prot/lrandom default-rng))
  (^long [^long mx] (prot/lrandom default-rng mx))
  (^long [^long mn ^long mx] (prot/lrandom default-rng mn mx)))

(defmacro randval
  "Return value with given probability (default 0.5)"
  ([v1 v2]
   `(if (prot/brandom default-rng) ~v1 ~v2))
  ([prob v1 v2]
   `(if (prot/brandom default-rng ~prob) ~v1 ~v2))
  ([prob]
   `(prot/brandom default-rng ~prob))
  ([]
   `(prot/brandom default-rng)))

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

(defn roll-a-dice
  "Roll a dice with given sides"
  (^long [sides]
   (inc (irand sides)))
  (^long [dices sides]
   (reduce clojure.core/+ (repeatedly dices #(inc (irand sides)))))  )

;; rng versions

(defmacro randval-rng
  "Return value with given probability (default 0.5), for given rng"
  ([rng v1 v2]
   `(if (prot/brandom ~rng) ~v1 ~v2))
  ([rng prob v1 v2]
   `(if (prot/brandom ~rng ~prob) ~v1 ~v2))
  ([rng prob]
   `(prot/brandom ~rng ~prob))
  ([rng]
   `(prot/brandom ~rng)))

(defn flip-rng
  "Returns 1 with given probability, 0 otherwise, for given rng"
  (^long [rng p]
   (randval-rng rng p 1 0))
  (^long [rng]
   (randval-rng rng 0.5 1 0)))

(defn flipb-rng
  "Returns true with given probability, false otherwise, for given rng"
  ([rng p] (randval-rng rng p))
  ([rng] (randval-rng rng)))

(defn roll-a-dice-rng
  "Roll a dice with given sides and given rng"
  (^long [rng sides]
   (inc (irandom rng sides)))
  (^long [rng dices sides]
   (reduce clojure.core/+ (repeatedly dices #(inc (irandom rng sides))))))

;; generators

;; http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/#more-2165
(defn ball-random
  "Return random vector from a ball"
  ([^long dims] (ball-random default-rng dims))
  ([rng ^long dims]
   (let [u (double-array (repeatedly (+ dims 2) #(grandom rng)))
         ^doubles n (v/div u (v/mag u))]
     (case dims
       1 (aget n 0)
       2 (v/array->vec2 n)
       3 (v/array->vec3 n)
       4 (v/array->vec4 n)
       (vec (take dims n))))))

(defn- rv-generators
  "Generators from commons math and custom classes."
  [seq-generator ^long dimensions]
  (assert (case seq-generator
            :halton (m/<= 1 dimensions 40)
            :sobol (m/<= 1 dimensions 1000)
            :r2 (m/<= 1 dimensions 15)
            true) (str "Number of dimensions for " seq-generator " should be less or equal than "
                       ({:halton 40 :sobol 1000 :r2 15} seq-generator)))
  (let [^RandomVectorGenerator g (case seq-generator
                                   :halton (HaltonSequenceGenerator. dimensions)
                                   :sobol (SobolSequenceGenerator. dimensions)
                                   :sphere (UnitSphereRandomVectorGenerator. dimensions)
                                   :r2 (R2. dimensions))]
    (repeatedly (case dimensions
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
  (let [g (if (= seq-generator :gaussian)
            grand
            drand)]
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
        g (random-generators :default dimensions)]
    (map-indexed (fn [^long i v] (v/mult v (* c (m/pow (- (inc i) i0) p)))) g)))


;; Sequence creators

(defmulti
  ^{:doc "Create Sequence generator. See [[sequence-generators-list]] for names.

Values:

* `:r2`, `:halton`, `:sobol`, `:default`/`:uniform` - range `[0-1] for each dimension`
* `:gaussian` - from `N(0,1)` distribution
* `:sphere` -  from surface of unit sphere (ie. euclidean distance from origin equals 1.0)
* `:ball` - from an unit ball

Possible dimensions:

* `:r2` - 1-15
* `:halton` - 1-40
* `:sobol` - 1-1000
* the rest - 1+

See also [[jittered-sequence-generator]]."}
  sequence-generator (fn [seq-generator _] seq-generator))
(defmethod sequence-generator :halton [seq-generator dimensions] (rv-generators seq-generator dimensions))
(defmethod sequence-generator :sobol [seq-generator dimensions] (rv-generators seq-generator dimensions))
(defmethod sequence-generator :r2 [seq-generator dimensions] (rv-generators seq-generator dimensions))
(defmethod sequence-generator :sphere [seq-generator dimensions] (rv-generators seq-generator dimensions))
(defmethod sequence-generator :gaussian [seq-generator dimensions] (random-generators seq-generator dimensions))

(defmethod sequence-generator :default [seq-generator dimensions] (random-generators seq-generator dimensions))
(defmethod sequence-generator :ball [_ dimensions] (repeatedly (partial ball-random dimensions)))

(defn jittered-sequence-generator
  "Create jittered sequence generator.

  Suitable for `:r2`, `:sobol` and `:halton` sequences.

  `jitter` parameter range is from `0` (no jitter) to `1` (full jitter). Default: 0.25.

  See also [[sequence-generator]]."
  ([seq-generator ^long dimensions] (jittered-sequence-generator seq-generator dimensions 0.25))
  ([seq-generator ^long dimensions ^double jitter]
   (let [s (sequence-generator seq-generator dimensions) 
         [j mod-fn] (if (#{:sphere :gaussian} seq-generator)
                      (let [j (sequence-generator :gaussian dimensions)
                            jitter-low (* m/SQRTPI 0.5 0.25 jitter)]
                        [j (if (m/one? dimensions)
                             (fn [^double v ^double vj] (+ v (* jitter-low vj)))
                             (fn [v vj] (v/add v (v/mult vj jitter-low))))])
                      (let [j (jitter-generator seq-generator dimensions jitter)]
                        [j (if (m/one? dimensions)
                             (fn [^double v ^double vj] (m/frac (+ v vj)))
                             (fn [v vj] (v/fmap (v/add v vj) m/frac)))]))]
     (map mod-fn s j))))

(def ^{:doc "List of random sequence generator. See [[sequence-generator]]."}
  sequence-generators-list (keys (methods sequence-generator)))

;; ## Noise

(def ^{:doc "List of possible noise interpolations as a map of names and values."}
  noise-interpolations {:none NoiseConfig/INTERPOLATE_NONE
                        :linear NoiseConfig/INTERPOLATE_LINEAR
                        :hermite NoiseConfig/INTERPOLATE_HERMITE
                        :quintic NoiseConfig/INTERPOLATE_QUINTIC})

(def ^{:doc "List of possible noise types as a map of names and values."}
  noise-types {:value NoiseConfig/NOISE_VALUE
               :gradient NoiseConfig/NOISE_GRADIENT
               :simplex NoiseConfig/NOISE_SIMPLEX})

(defn- noise-config-obj
  "Create noise configuration object based on map."
  [{:keys [seed noise-type interpolation octaves lacunarity gain normalize?]}]
  (NoiseConfig. seed
                (or (noise-types noise-type) NoiseConfig/NOISE_GRADIENT)
                (or (noise-interpolations interpolation) NoiseConfig/INTERPOLATE_HERMITE)
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

(defonce ^:private perlin-noise-config (noise-config {:interpolation :quintic}))
(defonce ^:private simplex-noise-config (noise-config {:noise-type :simplex}))
(defonce ^:private value-noise-config (noise-config {:noise-type :value}))

(defn vnoise
  "Value Noise.

  6 octaves, Hermite interpolation (cubic, h01)."
  (^double [^double x] (FBM/noise value-noise-config x))
  (^double [^double x ^double y] (FBM/noise value-noise-config x y))
  (^double [^double x ^double y ^double z] (FBM/noise value-noise-config x y z)))

(defn noise
  "Improved Perlin Noise.

  6 octaves, quintic interpolation."
  (^double [^double x] (FBM/noise perlin-noise-config x))
  (^double [^double x ^double y] (FBM/noise perlin-noise-config x y))
  (^double [^double x ^double y ^double z] (FBM/noise perlin-noise-config x y z)))

(defn simplex
  "Simplex noise. 6 octaves."
  (^double [^double x] (FBM/noise simplex-noise-config x))
  (^double [^double x ^double y] (FBM/noise simplex-noise-config x y))
  (^double [^double x ^double y ^double z] (FBM/noise simplex-noise-config x y z)))

(defmacro ^:private gen-noise-function
  "Generate various noise as static function"
  [noise-type method]
  `(defn ~noise-type
     ~(str "Create " noise-type " function with optional configuration.")
     ([] (~noise-type {}))
     ([cfg#]
      (let [ncfg# (noise-config cfg#)]
        (fn
          ([x#] (~method ncfg# x#))
          ([x# y#] (~method ncfg# x# y#))
          ([x# y# z#] (~method ncfg# x# y# z#)))))))

(gen-noise-function single-noise Noise/noise)
(gen-noise-function fbm-noise FBM/noise)
(gen-noise-function billow-noise Billow/noise)
(gen-noise-function ridgedmulti-noise RidgedMulti/noise)

(defn- make-warp-1d
  [n ^double scale ^long depth]
  (let [warp-noise-1d-proto (fn warp-noise-1d
                              (^double [^double x ^long depth]
                               (if (zero? depth)
                                 (n x)
                                 (let [q1 (* scale ^double (warp-noise-1d (+ x depth 0.321) (dec depth)))]
                                   (n (+ x q1))))))]
    (fn [^double x] (warp-noise-1d-proto x depth))))

(defn- make-warp-2d
  [n ^double scale ^long depth]
  (let [warp-noise-2d-proto (fn warp-noise-2d
                              (^double [^double x ^double y ^long depth]
                               (if (zero? depth)
                                 (n x y)
                                 (let [q1 (* scale ^double (warp-noise-2d (+ x depth 0.321) (+ y depth 4.987) (dec depth)))
                                       q2 (* scale ^double (warp-noise-2d (+ x depth 3.591) (+ y depth -2.711) (dec depth)))]
                                   (n (+ x q1) (+ y q2))))))]
    (fn [^double x ^double y] (warp-noise-2d-proto x y depth))))

(defn- make-warp-3d
  [n ^double scale ^long depth]
  (let [warp-noise-3d-proto (fn warp-noise-3d
                              (^double [^double x ^double y ^double z ^long depth]
                               (if (zero? depth)
                                 (n x y z)
                                 (let [q1 (* scale ^double (warp-noise-3d (+ x depth 0.321) (+ y depth 4.987) (+ z depth 2.12) (dec depth)))
                                       q2 (* scale ^double (warp-noise-3d (+ x depth 3.591) (+ y depth -2.711) (+ z depth -5.4321) (dec depth)))
                                       q3 (* scale ^double (warp-noise-3d (+ x depth -1.591) (+ y depth 12.1711) (+ z depth 3.1) (dec depth)))]
                                   (n (+ x q1) (+ y q2) (+ z q3))))))]
    (fn [^double x ^double y ^double z] (warp-noise-3d-proto x y z depth))))

(defn warp-noise-fn
  "Create warp noise (see [Inigo Quilez article](http://www.iquilezles.org/www/articles/warp/warp.htm)).

  Parameters:

  * noise function, default: vnoise
  * scale factor, default: 4.0
  * depth (1 or 2), default 1

  Normalization of warp noise depends on normalization of noise function."
  ([noise ^double scale ^long depth]
   (let [n1 (make-warp-1d noise scale depth)
         n2 (make-warp-2d noise scale depth)
         n3 (make-warp-3d noise scale depth)]
     (fn
       (^double [^double x] (n1 x))
       (^double [^double x ^double y] (n2 x y))
       (^double [^double x ^double y ^double z] (n3 x y z)))))
  ([noise ^double scale] (warp-noise-fn noise scale 1))
  ([noise] (warp-noise-fn noise 4.0 1))
  ([] (warp-noise-fn vnoise 4.0 1)))

(defonce ^{:doc "List of possible noise generators as a map of names and functions."}
  noise-generators
  {:fbm fbm-noise
   :single single-noise
   :billow billow-noise
   :ridgemulti ridgedmulti-noise})

(defn random-noise-cfg
  "Create random noise configuration.

  Optional map with fixed values."
  ([pre-config]
   (merge {:seed (irand)
           :generator (rand-nth [:single :fbm :billow :ridgemulti])
           :noise-type (rand-nth (keys noise-types))
           :interpolation (rand-nth (keys noise-interpolations))
           :octaves (irand 1 10)
           :lacunarity (drand 1.5 2.5)
           :gain (drand 0.2 0.8)
           :warp-scale (randval 0.8 0.0 (randval 0.5 4.0 (drand 0.1 10.0)))
           :warp-depth (randval 0.8 1 (irand 1 4))
           :normalize? true} pre-config))
  ([] (random-noise-cfg nil)))

(defn random-noise-fn
  "Create random noise function from all possible options.

  Optionally provide own configuration `cfg`. In this case one of 4 different blending methods will be selected."
  ([cfg]
   (let [cfg (random-noise-cfg cfg)
         gen-fn (noise-generators (get cfg :generator :fbm))
         noise (gen-fn cfg)]
     (if (pos? ^double (:warp-scale cfg))
       (warp-noise-fn noise (:warp-scale cfg) (:warp-depth cfg))
       noise)))
  ([] (random-noise-fn nil)))

;; ### Discrete noise

(defn discrete-noise
  "Discrete noise. Parameters:

  * X (long)
  * Y (long, optional)

  Returns double value from [0,1] range"
  (^double [^long X ^long Y] (Discrete/value X Y))
  (^double [^long X] (Discrete/value X 0)))

;; Distribution

(defmulti
  ^{:doc "Create distribution object.

* First parameter is distribution as a `:key`.
* Second parameter is a map with configuration.

All distributions accept `rng` under `:rng` key (default: [[default-rng]]) and some of them accept `inverse-cumm-accuracy` (default set to `1e-9`)."}
  distribution (fn ([k _] k) ([k] k)))

(extend Object
  prot/DistributionIdProto
  {:distribution? (constantly false)})

(defn distribution?
  "Checks if `distr` is a distribution object."
  [distr]
  (prot/distribution? distr))

;; protocol proxies
(defn cdf
  "Cumulative probability."
  (^double [d v] (prot/cdf d v))
  (^double [d v1 v2] (prot/cdf d v1 v2)))

(defn ccdf
  "Complementary cumulative probability."
  ^double [d v] (- 1.0 (cdf d v)))

(defn pdf
  "Density"
  ^double [d v] (prot/pdf d v))

(defn lpdf
  "Log density"
  ^double [d v] (prot/lpdf d v))

(defn icdf
  "Inverse cumulative probability"
  [d ^double v] (prot/icdf d v))

(defn probability
  "Probability (PMF)"
  ^double [d v] (prot/probability d v))

(defn sample
  "Random sample"
  [d] (prot/sample d))

(defn dimensions
  "Distribution dimensionality"
  ^long [d] (prot/dimensions d))

(defn source-object
  "Returns Java or proxy object from backend library (if available)"
  [d] (prot/source-object d))

(defn continuous?
  "Does distribution support continuous domain?"
  [d] (prot/continuous? d))

(defn observe1
  "Log of probability/density of the value. Alias for [[lpdf]]."
  ^double [d v]
  (prot/lpdf d v))

(defn log-likelihood
  "Log likelihood of samples"
  ^double [d vs] 
  (reduce (fn [^double s ^double v] (if (m/invalid-double? s)
                                     (reduced s)
                                     (+ s v))) 0.0 (map #(prot/lpdf d %) vs)))

(defmacro observe
  "Log likelihood of samples. Alias for [[log-likelihood]]."
  [d vs]
  `(log-likelihood ~d ~vs))

(defn likelihood
  "Likelihood of samples"
  ^double [d vs]
  (m/exp (log-likelihood d vs)))

(defn mean
  "Distribution mean"
  ^double [d] (prot/mean d))

(defn means
  "Distribution means (for multivariate distributions)"
  [d] (prot/means d))

(defn variance
  "Distribution variance"
  ^double [d] (prot/variance d))

(defn covariance
  "Distribution covariance matrix (for multivariate distributions)"
  [d] (prot/covariance d))

(defn lower-bound
  "Distribution lowest supported value"
  ^double [d] (prot/lower-bound d))

(defn upper-bound
  "Distribution highest supported value"
  ^double [d] (prot/upper-bound d))

(defn distribution-id
  "Distribution identifier as keyword."
  [d] (prot/distribution-id d))

(defn distribution-parameters
  "Distribution highest supported value.

  When `all?` is true, technical parameters are included, ie: `:rng` and `:inverser-cumm-accuracy`."
  ([d] (distribution-parameters d false))
  ([d all?]
   (let [d' (if (keyword? d) (distribution d) d)]
     (if-not all?
       (-> (prot/distribution-parameters d')
           (set)
           (disj :rng :inverse-cumm-accuracy :epsilon :max-iterations)
           (vec))
       (prot/distribution-parameters d')))))

(defn integrate-pdf
  "Integrate PDF function, returns CDF and iCDF

  Parameters:
  * `pdf-func` - univariate function
  * `mn` - lower bound for integration, value of pdf-func should be 0.0 at this point
  * `mx` - upper bound for integration
  * `steps` - how much subintervals to integrate (default 1000)
  * `min-iterations` - minimum iterations for RombergIntegrator (default 3)
  * `interpolator` - interpolation method between integrated points (default :spline)

  Possible interpolation methods: `:linear` (default), `:spline`, `:monotone` or any function from `fastmath.interpolation`"
  ([pdf-func mn mx steps]
   (integrate-pdf pdf-func {:mn mn :mx mx :steps steps}))
  ([pdf-func {:keys [^double mn ^double mx ^long steps
                     interpolator ^long min-iterations]
              :or {mn 0.0 mx 1.0 steps 1000
                   min-iterations 3 interpolator :linear}}]
   (let [step (/ (- mx mn) steps)
         u-pdf-func (reify UnivariateFunction
                      (value [_ x] (pdf-func x)))
         ^RombergIntegrator romberg-integrator (RombergIntegrator. (max 2 min-iterations) RombergIntegrator/ROMBERG_MAX_ITERATIONS_COUNT)
         ;; go through the intervals and integrate them, assuming that kde of `mn` is 0.0
         points (second (reduce (fn [[^double curr lst] [^double x1 ^double x2]]
                                  (let [i (.integrate romberg-integrator Integer/MAX_VALUE u-pdf-func x1 x2) ;; integration can be very slow on very narrow spikes
                                        curr-new (m/constrain (+ i curr) 0.0 1.0)
                                        res (if (> curr-new curr)
                                              [curr-new (conj lst [x2 curr-new])]
                                              [curr-new lst])]
                                    (if (== curr-new 1.0) ;; avoid overflow of integration (usually it's underestimated)
                                      (reduced res)
                                      res))) [0.0 [[mn 0.0]]]
                                (partition 2 1 (m/slice-range mn mx steps))))
         [^double lx ^double ly] (last points)
         points (if (< ly 1.0)
                  (if (< lx mx)
                    (conj points [mx 1.0])
                    (conj points [(+ mx step) 1.0]))
                  (if (< lx mx)
                    (conj points [mx (m/next-double 1.0)])
                    (conj points [(+ mx step) (m/next-double 1.0)]))) ;; fix upper endpoint
         xs (m/seq->double-array (map first points))
         ys (m/seq->double-array (map second points))
         intpol (case interpolator
                  :linear i/linear-smile
                  :spline i/cubic-spline
                  :monotone i/monotone
                  (if (fn? interpolator) interpolator i/linear))]
     ;; interpolate points lineary, return cdf and icdf
     [(intpol xs ys)
      (intpol ys xs)])))

;; apache commons math
(extend RealDistribution
  prot/DistributionProto
  {:cdf (fn
          (^double [^RealDistribution d ^double v] (.cumulativeProbability d v))
          (^double [^RealDistribution d ^double v1 ^double v2] (.cumulativeProbability d v1 v2)))
   :pdf (fn ^double [^RealDistribution d ^double v] (.density d v))
   :lpdf (fn ^double [^AbstractRealDistribution d ^double v] (.logDensity d v))
   :icdf (fn ^double [^RealDistribution d ^double p] (.inverseCumulativeProbability d p))
   :probability (fn ^double [^RealDistribution d ^double p] (.density d p))
   :sample (fn ^double [^RealDistribution d] (.sample d))
   :dimensions (constantly 1)
   :source-object identity
   :continuous? (constantly true)} 
  prot/UnivariateDistributionProto
  {:mean (fn ^double [^RealDistribution d] (.getNumericalMean d))
   :variance (fn ^double [^RealDistribution d] (.getNumericalVariance d))
   :lower-bound (fn ^double [^RealDistribution d] (.getSupportLowerBound d))
   :upper-bound (fn ^double [^RealDistribution d] (.getSupportUpperBound d))}
  prot/RNGProto
  {:drandom (fn ^double [^RealDistribution d] (.sample d))
   :frandom (fn [^RealDistribution d] (unchecked-float (.sample d)))
   :lrandom (fn ^long [^RealDistribution d] (m/round-even (.sample d)))
   :irandom (fn ^long [^RealDistribution d] (unchecked-int (m/round-even (.sample d))))
   :->seq (fn
            ([^RealDistribution d] (repeatedly #(.sample d)))
            ([^RealDistribution d n] (repeatedly n #(.sample d))))
   :set-seed! (fn [^RealDistribution d ^long seed] (.reseedRandomGenerator d seed) d)})

;; ssj

(defn- reify-continuous-ssj
  [^ContinuousDistribution d ^RandomGenerator rng nm & ks]
  (let [kss (vec (conj ks :rng))]
    (reify
      prot/DistributionProto
      (pdf [_ v] (.density d v))
      (lpdf [_ v] (m/log (.density d v)))
      (probability [_ v] (.density d v))
      (cdf [_ v] (.cdf d v))
      (cdf [_ v1 v2] (- (.cdf d v2) (.cdf d v1)))
      (icdf [_ v] (.inverseF d v))
      (sample [_] (.inverseF d (prot/drandom rng)))
      (dimensions [_] 1)
      (source-object [_] d)
      (continuous? [_] true)
      prot/DistributionIdProto
      (distribution? [_] true)
      (distribution-id [_] nm)
      (distribution-parameters [_] kss)
      prot/UnivariateDistributionProto
      (mean [_] (.getMean d))
      (variance [_] (.getVariance d))
      (lower-bound [_] (.getXinf d))
      (upper-bound [_] (.getXsup d))
      prot/RNGProto
      (drandom [_] (.inverseF d (prot/drandom rng)))
      (frandom [_] (unchecked-float (.inverseF d (prot/drandom rng))))
      (lrandom [_] (m/round-even (.inverseF d (prot/drandom rng))))
      (irandom [_] (unchecked-int (m/round-even (.inverseF d (prot/drandom rng)))))
      (->seq [_] (repeatedly #(.inverseF d (prot/drandom rng))))
      (->seq [_ n] (repeatedly n #(.inverseF d (prot/drandom rng))))
      (set-seed! [d seed] (prot/set-seed! rng seed) d))))

(defn- reify-continuous-ssj-no-pdf
  [^ContinuousDistribution d ^RandomGenerator rng nm & ks]
  (let [kss (vec (conj ks :rng))]
    (reify
      prot/DistributionProto
      (pdf [_ v] (/ (- (.cdf d (+ ^double v 0.5e-6))
                       (.cdf d (- ^double v 0.5e-6)))
                    1.0e-6))
      (lpdf [rd v] (m/log (prot/pdf rd v)))
      (probability [rd v] (prot/pdf rd v))
      (cdf [_ v] (.cdf d v))
      (cdf [_ v1 v2] (- (.cdf d v2) (.cdf d v1)))
      (icdf [_ v] (.inverseF d v))
      (sample [_] (.inverseF d (prot/drandom rng)))
      (dimensions [_] 1)
      (source-object [_] d)
      (continuous? [_] true)
      prot/DistributionIdProto
      (distribution? [_] true)
      (distribution-id [_] nm)
      (distribution-parameters [_] kss)
      prot/UnivariateDistributionProto
      (mean [_] (.getMean d))
      (variance [_] (.getVariance d))
      (lower-bound [_] (.getXinf d))
      (upper-bound [_] (.getXsup d))
      prot/RNGProto
      (drandom [_] (.inverseF d (prot/drandom rng)))
      (frandom [_] (unchecked-float (.inverseF d (prot/drandom rng))))
      (lrandom [_] (m/round-even (.inverseF d (prot/drandom rng))))
      (irandom [_] (unchecked-int (m/round-even (.inverseF d (prot/drandom rng)))))
      (->seq [_] (repeatedly #(.inverseF d (prot/drandom rng))))
      (->seq [_ n] (repeatedly n #(.inverseF d (prot/drandom rng))))
      (set-seed! [d seed] (prot/set-seed! rng seed) d))))

(defn- reify-integer-ssj
  [^DiscreteDistributionInt d ^RandomGenerator rng nm m]
  (let [kss (vec (conj (keys m) :rng))
        icdf-fn (case nm
                  :logarithmic (let [upper (/ 25.0 (- 1.0 ^double (get m :theta 0.5)))
                                     r (range 0 (inc upper))]                                 
                                 (i/step-before (rest (reductions
                                                       (fn [^double s ^double v]
                                                         (+ s (.prob d v))) 0.0 r)) r))
                  (fn [^double v] (.inverseF d v)))]
    (reify
      prot/DistributionProto
      (pdf [_ v] (.prob d (m/floor v)))
      (lpdf [_ v] (m/log (.prob d (m/floor v))))
      (cdf [_ v] (.cdf d (m/floor v)))
      (cdf [_ v1 v2] (- (.cdf d (m/floor v2)) (.cdf d (m/floor v1))))
      (icdf [_ v] (unchecked-long (icdf-fn v)))
      (probability [_ v] (.prob d (m/floor v)))
      (sample [_] (unchecked-long (icdf-fn (prot/drandom rng))))
      (dimensions [_] 1)
      (source-object [_] d)
      (continuous? [_] false)
      prot/DistributionIdProto
      (distribution? [_] true)
      (distribution-id [_] nm)
      (distribution-parameters [_] kss)
      prot/UnivariateDistributionProto
      (mean [_] (.getMean d))
      (variance [_] (.getVariance d))
      (lower-bound [_] (.getXinf d))
      (upper-bound [_] (.getXsup d))
      prot/RNGProto
      (drandom [_] (.inverseF d (prot/drandom rng)))
      (frandom [_] (unchecked-float (icdf-fn (prot/drandom rng))))
      (lrandom [_] (unchecked-long (icdf-fn (prot/drandom rng))))
      (irandom [_] (unchecked-int (icdf-fn (prot/drandom rng))))
      (->seq [_] (repeatedly #(unchecked-long (icdf-fn (prot/drandom rng)))))
      (->seq [_ n] (repeatedly n #(unchecked-long (icdf-fn (prot/drandom rng)))))
      (set-seed! [d seed] (prot/set-seed! rng seed) d))))

(extend IntegerDistribution
  prot/DistributionProto
  {:cdf (fn
          (^double [^IntegerDistribution d ^double v] (.cumulativeProbability d (m/floor v)))
          (^double [^IntegerDistribution d ^double v1 ^double v2] (.cumulativeProbability d (m/floor v1) (m/floor v2))))
   :icdf (fn ^long [^IntegerDistribution d ^double p] (.inverseCumulativeProbability d p))
   :pdf (fn ^double [^IntegerDistribution d ^double p] (.probability d (m/floor p)))
   :lpdf (fn ^double [^AbstractIntegerDistribution d ^double p] (.logProbability d (m/floor p)))
   :probability (fn ^double [^IntegerDistribution d ^double p] (.probability d (m/floor p)))
   :sample (fn ^long [^IntegerDistribution d] (.sample d))
   :dimensions (constantly 1)
   :source-object identity
   :continuous? (constantly false)}
  prot/UnivariateDistributionProto
  {:mean (fn ^double [^IntegerDistribution d] (.getNumericalMean d))
   :variance (fn ^double [^IntegerDistribution d] (.getNumericalVariance d))
   :lower-bound (fn ^long [^IntegerDistribution d] (.getSupportLowerBound d))
   :upper-bound (fn ^long [^IntegerDistribution d] (.getSupportUpperBound d))}
  prot/RNGProto
  {:drandom (fn ^double [^IntegerDistribution d] (unchecked-double (.sample d)))
   :frandom (fn [^IntegerDistribution d] (unchecked-float (.sample d)))
   :lrandom (fn ^long [^IntegerDistribution d] (unchecked-long (.sample d)))
   :irandom (fn ^long [^IntegerDistribution d] (.sample d))
   :->seq (fn
            ([^IntegerDistribution d] (repeatedly #(.sample d)))
            ([^IntegerDistribution d n] (repeatedly n #(.sample d))))
   :set-seed! (fn [^IntegerDistribution d ^long seed] (.reseedRandomGenerator d seed) d)})

(extend MultivariateNormalDistribution
  prot/DistributionProto
  {:pdf (fn ^double [^MultivariateNormalDistribution d v] (.density d (m/seq->double-array v)))
   :lpdf (fn ^double [^MultivariateNormalDistribution d v] (m/log (.density d (m/seq->double-array v))))
   :sample (fn [^MultivariateNormalDistribution d] (vec (.sample d)))
   :dimensions (fn ^long [^MultivariateNormalDistribution d] (.getDimension d))
   :source-object identity
   :continuous? (constantly true)}
  prot/MultivariateDistributionProto
  {:means (fn [^MultivariateNormalDistribution d] (vec (.getMeans d)))
   :covariance (fn [^MultivariateNormalDistribution d]
                 (let [^org.apache.commons.math3.linear.Array2DRowRealMatrix cv (.getCovariances d)]
                   (m/double-double-array->seq (.getDataRef cv))))}
  prot/RNGProto
  {:->seq (fn
            ([^MultivariateNormalDistribution d] (repeatedly #(vec (.sample d))))
            ([^MultivariateNormalDistribution d n] (repeatedly n #(vec (.sample d)))))
   :set-seed! (fn [^MultivariateNormalDistribution d ^long seed] (.reseedRandomGenerator d seed) d)})

(defmacro ^:private make-acm-distr
  [nm obj ks vs]
  (let [or-map (zipmap ks vs)] 
    `(do
       (extend ~obj
         prot/DistributionIdProto
         {:distribution? (constantly true)
          :distribution-id (fn [d#] ~nm)
          :distribution-parameters (fn [d#] [~@(conj (map keyword ks) :rng)])})
       (defmethod distribution ~nm
         ([n# {:keys [~@ks]
               :or ~or-map
               :as all#}]
          (let [^RandomGenerator r# (or (:rng all#) (rng :jvm))]
            (new ~obj r# ~@ks)))
         ([n#] (distribution ~nm {}))))))

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
                [numerator-degrees-of-freedom denominator-degrees-of-freedom inverse-cumm-accuracy]
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
                [degrees-of-freedom inverse-cumm-accuracy]
                [1.0 TDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY])

(make-acm-distr :triangular TriangularDistribution [a c b] [-1.0 0.0 1.0])
(make-acm-distr :uniform-real UniformRealDistribution [^double lower ^double upper] [0.0 1.0])

(make-acm-distr :weibull WeibullDistribution
                [alpha beta inverse-cumm-accuracy]
                [2.0 1.0 WeibullDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY])

(extend EmpiricalDistribution
  prot/DistributionIdProto
  {:distribution? (constantly true)
   :distribution-id (fn [_] :empirical)
   :distribution-parameters (fn [_] [:rng :bin-count :data])})

(defmethod distribution :empirical
  ([_ {:keys [^long bin-count data]
       :or {data [1.0]}
       :as all}]
   (let [bin-count (or bin-count (unchecked-long (max (* 0.1 (count data)) 1.0)))
         ^RandomGenerator r (or (:rng all) (rng :jvm))
         ^EmpiricalDistribution d (EmpiricalDistribution. bin-count r)]
     (.load d ^doubles (m/seq->double-array data))
     d))
  ([_] (distribution :empirical {})))

(extend EnumeratedRealDistribution
  prot/DistributionProto
  {:cdf (fn
          (^double [^EnumeratedRealDistribution d ^double v] (.cumulativeProbability d v))
          (^double [^EnumeratedRealDistribution d ^double v1 ^double v2] (.probability d v1 v2)))
   :icdf (fn ^double [^EnumeratedRealDistribution d ^double p] (.inverseCumulativeProbability d p))
   :pdf (fn ^double [^EnumeratedRealDistribution d ^double p] (.probability d p))
   :lpdf (fn ^double [^EnumeratedRealDistribution d ^double p] (.logDensity d p))
   :probability (fn ^double [^EnumeratedRealDistribution d ^double p] (.probability d p))
   :sample (fn ^double [^EnumeratedRealDistribution d] (.sample d))
   :dimensions (constantly 1)
   :source-object identity
   :continuous? (constantly false)}
  prot/DistributionIdProto
  {:distribution? (constantly true)
   :distribution-id (fn [_] :enumerated-real)
   :distribution-parameters (fn [_] [:rng :data :probabilities])})

(defmethod distribution :enumerated-real
  ([_ {:keys [data probabilities]
       :or {data [1.0]}
       :as all}]
   (let [^RandomGenerator r (or (:rng all) (rng :jvm))]
     (if probabilities
       (EnumeratedRealDistribution. r (m/seq->double-array data) (m/seq->double-array probabilities))
       (EnumeratedRealDistribution. r ^doubles (m/seq->double-array data)))))
  ([_] (distribution :enumerated-real {})))

;; integer

(defmethod distribution :bernoulli
  ([_ {:keys [^double p]
       :or {p 0.5}
       :as all}]
   (BinomialDistribution. (or (:rng all) (rng :jvm)) 1 p))
  ([_] (distribution :bernoulli {})))

(extend EnumeratedIntegerDistribution
  prot/DistributionIdProto
  {:distribution? (constantly true)
   :distribution-id (fn [_] :enumerated-int)
   :distribution-parameters (fn [_] [:data :probabilities :rng])})

(defmethod distribution :enumerated-int
  ([_ {:keys [data probabilities]
       :or {data [1]}
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
(make-acm-distr :pascal PascalDistribution [r p] [20 0.5])
(make-acm-distr :poisson PoissonDistribution
                [p epsilon max-iterations]
                [0.5 PoissonDistribution/DEFAULT_EPSILON PoissonDistribution/DEFAULT_MAX_ITERATIONS])
(make-acm-distr :uniform-int UniformIntegerDistribution [lower upper] [0 Integer/MAX_VALUE])
(make-acm-distr :zipf ZipfDistribution [number-of-elements exponent] [100 3.0])

;; ssj

(defmacro ^:private make-ssj-distr
  [rf nm obj ks vs]
  (let [or-map (zipmap ks vs)
        k-map (zipmap (map keyword ks) vs)]
    `(defmethod distribution ~nm
       ([n# {:keys [~@ks]
             :or ~or-map
             :as all#}]
        (let [^RandomGenerator r# (or (:rng all#) (rng :jvm))]
          (~rf (new ~obj ~@ks) r# ~nm (merge ~k-map all#))))
       ([n#] (distribution ~nm {})))))

(defmacro ^:private make-ssjc-distr
  [nm obj ks vs] `(make-ssj-distr reify-continuous-ssj ~nm ~obj ~ks ~vs))
(defmacro ^:private make-ssjc-distr-no-pdf
  [nm obj ks vs] `(make-ssj-distr reify-continuous-ssj-no-pdf ~nm ~obj ~ks ~vs))
(defmacro ^:private make-ssji-distr
  [nm obj ks vs] `(make-ssj-distr reify-integer-ssj ~nm ~obj ~ks ~vs))

(make-ssjc-distr :anderson-darling AndersonDarlingDistQuick [n] [1.0])
(make-ssjc-distr :inverse-gamma InverseGammaDist [alpha beta] [2.0 1.0])
(make-ssjc-distr :chi ChiDist [nu] [1.0])
(make-ssjc-distr :chi-squared-noncentral ChiSquareNoncentralDist [nu lambda] [1.0 1.0])
(make-ssjc-distr-no-pdf :cramer-von-mises CramerVonMisesDist [n] [1.0])
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
  ([k {:keys [lambdas]
       :or {lambdas [1.0]}
       :as all}]
   (reify-continuous-ssj (HypoExponentialDist. (m/seq->double-array lambdas)) (or (:rng all) (rng :jvm)) k :lambdas))
  ([_] (distribution :hypoexponential {})))

(defmethod distribution :reciprocal-sqrt
  ([_ {:keys [^double a]
       :or {a 0.5}
       :as all}]
   (let [f (* 2.0 (m/sqrt a))
         icdf-fn (fn [^double x]
                   (cond
                     (m/not-pos? x) a
                     :else (m/sq (* 0.5 (+ (min x 1.0) f)))))
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
       prot/DistributionProto
       (pdf [_ v] (if (<= a ^double v b) (/ (m/sqrt v)) 0.0))
       (lpdf [d v] (m/log (prot/pdf d v)))
       (cdf [_ v] (cond
                    (< ^double v a) 0.0
                    (> ^double v b) 1.0
                    :else (- (* 2.0 (m/sqrt ^double v)) f)))
       (cdf [d v1 v2] (- ^double (prot/cdf d v2) ^double (prot/cdf d v1)))
       (icdf [_ v] (icdf-fn v))
       (probability [d v] (prot/pdf d v))
       (sample [_] (icdf-fn (prot/drandom r)))
       (dimensions [_] 1)
       (source-object [d] d)
       (continuous? [_] true)
       prot/DistributionIdProto
       (distribution? [_] true)
       (distribution-id [_] :reciprocal-sqrt)
       (distribution-parameters [_] [:a :rng])
       prot/UnivariateDistributionProto
       (mean [_] m)
       (variance [_] v)
       (lower-bound [_] a)
       (upper-bound [_] b)
       prot/RNGProto
       (drandom [_] (icdf-fn (prot/drandom r)))
       (frandom [_] (unchecked-float (icdf-fn (prot/drandom r))))
       (lrandom [_] (m/round-even (icdf-fn (prot/drandom r))))
       (irandom [_] (unchecked-int (m/round-even (icdf-fn (prot/drandom r)))))
       (->seq [_] (repeatedly #(icdf-fn (prot/drandom r))))
       (->seq [_ n] (repeatedly n #(icdf-fn (prot/drandom r))))
       (set-seed! [d seed] (prot/set-seed! r seed) d))))
  ([_] (distribution :reciprocal-sqrt {})))

;;

(extend MultivariateNormalDistribution
  prot/DistributionIdProto
  {:distribution? (constantly true)
   :distribution-id (fn [_] :multi-normal)
   :distribution-parameters (fn [_] [:means :covariances :rng])})

(defmethod distribution :multi-normal
  ([_ {:keys [means covariances] :as all}]
   (let [covariances (cond
                       (and means (not covariances)) (for [id (range (count means))
                                                           :let [a (double-array (count means))]]
                                                       (do (aset a id 1.0)
                                                           a))
                       (not covariances) [[1.0 0.0] [0.0 1.0]]
                       :else covariances)
         means (if-not means (repeat (count (first covariances)) 0.0) means)]
     (assert (= (count means) (count (first covariances)))
             "Means and covariances sizes do not match.")
     (MultivariateNormalDistribution. (or (:rng all) (rng :jvm)) (m/seq->double-array means) (m/seq->double-double-array covariances))))
  ([_] (distribution :multi-normal {})))

(def ^{:const true :private true :tag 'double} zero+epsilon (m/next-double 0.0))
(def ^{:const true :private true :tag 'double} one-epsilon (m/prev-double 1.0))

(defn- dirichlet-rev-log-beta
  ^double [alpha]
  (let [d (m/log-gamma (reduce m/fast+ alpha))
        ^double n (reduce m/fast+ (map #(m/log-gamma %) alpha))]
    (- d n)))

#_(defn- dirichlet-lpdf
    ^double [alpha- values ^double lbeta]
    (if (every? #(< 0.0 ^double % 1.0) values)
      (let [^double p (reduce m/fast+ (mapv (fn [^double ai ^double x] 
                                              (* ai (m/log x))) alpha- values))]
        (+ lbeta p))
      ##-Inf))

(defn- dirichlet-lpdf
  ^double [alpha- values ^double lbeta]
  (let [v (reduce m/fast+ lbeta (map (fn [^double ai ^double x] 
                                       (* ai (m/log x))) alpha- values))]
    (if (m/invalid-double? v) ##-Inf v)))

(defmethod distribution :dirichlet
  ([_ {:keys [alpha] :as all}]
   (let [alpha (cond
                 (nil? alpha) (double-array [1 1])
                 (seqable? alpha) (do (assert (> (count alpha) 1))
                                      (m/seq->double-array alpha))
                 (integer? alpha) (do (assert (> (int alpha) 1))
                                      (double-array alpha 1.0))
                 :else (double-array [1 1]))
         r (or (:rng all) (rng :jvm))
         sampler (mapv #(distribution :gamma {:shape % :scale 1.0 :rng r}) alpha)

         lbeta (dirichlet-rev-log-beta alpha)
         alpha- (map clojure.core/dec alpha)
         
         m (delay (seq (DirichletDist/getMean alpha)))
         cv (delay (mapv vec (DirichletDist/getCovariance alpha)))
         dim (count alpha)]
     (reify
       prot/DistributionProto
       (pdf [_ v] (m/exp (dirichlet-lpdf alpha- v lbeta)))
       (lpdf [_ v] (dirichlet-lpdf alpha- v lbeta))
       (probability [_ v] (m/exp (dirichlet-lpdf alpha- v lbeta)))
       (sample [_] (let [samples (map #(prot/sample %) sampler)
                         s (v/sum samples)]
                     (mapv (fn [^double v] (cond
                                            (zero? v) zero+epsilon
                                            (== v 1.0) one-epsilon
                                            :else v))
                           (if (> s 1.0e-6)
                             (v/div samples s)
                             (let [a (int-array dim)]
                               (aset ^ints a (irand dim) 1)
                               a)))))
       (dimensions [_] dim)
       (source-object [this] this)
       (continuous? [_] true)
       prot/DistributionIdProto
       (distribution? [_] true)
       (distribution-id [_] :dirichlet)
       (distribution-parameters [_] [:alpha :rng])
       prot/MultivariateDistributionProto
       (means [_] @m)
       (covariance [_] @cv)
       prot/RNGProto
       (->seq [d] (repeatedly #(prot/sample d)))
       (->seq [d n] (repeatedly n #(prot/sample d)))
       (set-seed! [d seed] (prot/set-seed! r seed) d)))) 
  ([_] (distribution :dirichlet {})))

(defmethod distribution :multinomial
  ([_ {:keys [^int trials ps]
       :or {trials 20 ps [0.5 0.5]}:as all}]
   (let [ps (m/seq->double-array (v/div ps (v/sum ps)))
         r (or (:rng all) (rng :jvm))
         
         m (delay (seq (MultinomialDist/getMean trials ps)))
         cv (delay (mapv vec (MultinomialDist/getCovariance trials ps)))
         dim (count ps)
         binom-probs (mapv (fn [^double prob ^double sum]
                             (/ prob (- 1.0 sum))) ps (reductions m/fast+ 0 ps))]
     (reify
       prot/DistributionProto
       (pdf [_ v] (MultinomialDist/prob trials ps (int-array v)))
       (lpdf [_ v] (m/log (MultinomialDist/prob trials ps (int-array v))))
       (probability [_ v] (MultinomialDist/prob trials ps (int-array v)))
       (cdf [_ v] (MultinomialDist/cdf trials ps (int-array v)))
       (sample [_] (first (reduce (fn [[buf ^int curr] ^double prob]
                                    (let [res (int (prot/sample (distribution :binomial
                                                                              {:trials curr
                                                                               :p (m/constrain prob 0.0 1.0)})))]
                                      [(conj buf res) (- curr res)])) [[] trials] binom-probs)))
       (dimensions [_] dim)
       (source-object [this] this)
       (continuous? [_] false)
       prot/DistributionIdProto
       (distribution? [_] true)
       (distribution-id [_] :multinomial)
       (distribution-parameters [_] [:n :ps :rng])
       prot/MultivariateDistributionProto
       (means [_] @m)
       (covariance [_] @cv)
       prot/RNGProto
       (->seq [d] (repeatedly #(prot/sample d)))
       (->seq [d n] (repeatedly n #(prot/sample d)))
       (set-seed! [d seed] (prot/set-seed! r seed) d)))) 
  ([_] (distribution :multinomial {})))

;; 

(defn- find-first-non-zero
  ^double [f xs]
  (or (->> xs
           (map #(vector % (f %)))
           (filter #(pos? ^double (second %)))
           (ffirst))
      (first xs)))

(defn- narrow-range
  [kd [^double mn ^double mx ^double step] ^long steps]
  [(- (find-first-non-zero kd (m/slice-range mn mx steps)) step)
   (+ (find-first-non-zero kd (m/slice-range mx mn steps)) step)])

(defmethod distribution :continuous-distribution
  ([_ {:keys [data ^long steps kde bandwidth]
       :or {data [-1 0 1] steps 5000 kde :epanechnikov}
       :as all}]
   (let [[kd _ _ ^double mn ^double mx] (k/kernel-density kde data bandwidth true)
         step (/ (- mx mn) steps)
         [^double mn ^double mx] (narrow-range kd [mn mx step] (* 4 steps))
         [cdf-fn icdf-fn] (integrate-pdf kd (merge all {:mn mn :mx mx :steps steps}))
         r (or (:rng all) (rng :jvm))
         m (delay (StatUtils/mean (m/seq->double-array data)))
         v (delay (StatUtils/variance (m/seq->double-array data) ^double @m))]
     (reify
       prot/DistributionProto
       (pdf [_ v] (kd v))
       (lpdf [_ v] (m/log (kd v)))
       (cdf [_ v] (m/constrain ^double (cdf-fn v) 0.0 1.0))
       (cdf [d v1 v2] (- ^double (prot/cdf d v2) ^double (prot/cdf d v1)))
       (icdf [_ v] (icdf-fn (m/constrain ^double v 0.0 1.0)))
       (probability [_ v] (kd v))
       (sample [_] (icdf-fn (prot/drandom r)))
       (dimensions [_] 1)
       (source-object [d] d)
       (continuous? [_] true)
       prot/DistributionIdProto
       (distribution? [_] true)
       (distribution-id [_] :kde)
       (distribution-parameters [_] [:data :steps :kde :bandwidth :rng])
       prot/UnivariateDistributionProto
       (mean [_] @m)
       (variance [_] @v)
       (lower-bound [_] (- mn step))
       (upper-bound [_] (+ mx step))
       prot/RNGProto
       (drandom [_] (icdf-fn (prot/drandom r)))
       (frandom [_] (unchecked-float (icdf-fn (prot/drandom r))))
       (lrandom [_] (m/round-even (icdf-fn (prot/drandom r))))
       (irandom [_] (unchecked-int (m/round-even (icdf-fn (prot/drandom r)))))
       (->seq [_] (repeatedly #(icdf-fn (prot/drandom r))))
       (->seq [_ n] (repeatedly n #(icdf-fn (prot/drandom r))))
       (set-seed! [d seed] (prot/set-seed! r seed) d))))
  ([_] (distribution :continuous-distribution {})))

(defmethod distribution :kde
  [_ & r] (apply distribution :continuous-distribution r))

;;

(defn- diff-cdf
  ^double [cdf-fn ^double v1 ^double v2]
  (- ^double (cdf-fn v2) ^double (cdf-fn v1)))

(defmacro ^:private distribution-template
  [d-name {:keys [pdf? dimensions continuous? distribution-parameters mean variance lower-bound upper-bound]
           :or {pdf? false dimensions 1 continuous? true mean 'mean variance 'variance}} & let-body]
  `(defmethod distribution ~d-name
     ([_# ~'args]
      (let [~'r (or (:rng ~'args) (rng :jvm))
            ~@let-body]
        (reify
          prot/DistributionProto
          (pdf [_# ~'v] ~(if pdf? `(~'pdf-fn ~'v) `(m/exp (~'lpdf-fn ~'v))))
          (lpdf [_# ~'v] ~(if pdf? `(m/log (~'pdf-fn ~'v)) `(~'lpdf-fn ~'v)))
          (cdf [_# v#] (~'cdf-fn v#))
          (cdf [_# v1# v2#] (diff-cdf ~'cdf-fn v1# v2#))
          (icdf [_# v#] (~'icdf-fn v#))
          (probability [_# ~'v] ~(if pdf? `(~'pdf-fn ~'v) `(m/exp (~'lpdf-fn ~'v))))
          (sample [_#] (~'icdf-fn (prot/drandom ~'r)))
          (dimensions [_#] ~dimensions)
          (source-object [d#] d#)
          (continuous? [_#] ~continuous?)
          prot/DistributionIdProto
          (distribution? [_#] true)
          (distribution-id [_#] ~d-name)
          (distribution-parameters [_#] ~distribution-parameters)
          prot/UnivariateDistributionProto
          (mean [_#] ~mean)
          (variance [_#] ~variance)
          (lower-bound [_#] ~lower-bound)
          (upper-bound [_#] ~upper-bound)
          prot/RNGProto
          (drandom [_#] (~'icdf-fn (prot/drandom ~'r)))
          (frandom [_#] (unchecked-float (~'icdf-fn (prot/drandom ~'r))))
          (lrandom [_#] (m/round-even (~'icdf-fn (prot/drandom ~'r))))
          (irandom [_#] (unchecked-int (m/round-even (~'icdf-fn (prot/drandom ~'r)))))
          (->seq [_#] (repeatedly #(~'icdf-fn (prot/drandom ~'r))))
          (->seq [_# n#] (repeatedly n# #(~'icdf-fn (prot/drandom ~'r))))
          (set-seed! [d# seed#] (prot/set-seed! ~'r seed#) d#))))
     ([_#] (distribution ~d-name {}))))

(defn- discrete-binary-search
  ([cdf-fn ^double p [mid step]] (discrete-binary-search cdf-fn step p [0 mid]))
  ([cdf-fn ^long step ^double p [^long mn ^long mx]]
   (cond
     (> p ^double (cdf-fn mx)) (recur cdf-fn (* 2 step) p [mx (+ mx step)])
     (m/one? (- mx mn)) (if (>= ^double (cdf-fn mn) p) mn mx)
     :else (let [mid (/ (+ mn mx) 2)]
             (if (> ^double (cdf-fn mid) p)
               (recur cdf-fn step p [mn mid])
               (recur cdf-fn step p [mid mx]))))))

(distribution-template :negative-binomial
    {:continuous? false :lower-bound 0 :upper-bound Integer/MAX_VALUE :mean mmean
     :distribution-parameters [:r :p :rng]}
  {:keys [^double r ^double p] :or {r 20.0 p 0.5}} args
  p- (- 1.0 p)
  mmean (/ (* r p-) p)
  variance (/ mmean p)
  lgr (m/log-gamma r)
  lp- (m/log (- 1.0 p))
  lpr (* r (m/log p))
  lpdf-fn (fn [^long k]
            (if (neg? k)
              ##-Inf
              (+ (- (m/log-gamma (+ r k))
                    (+ (m/log-factorial k) lgr))
                 (* k lp-) lpr)))
  cdf-fn (fn [^double k]
           (if (neg? k)
             0.0
             (m/regularized-beta p r (inc (m/rint k)))))
  icdf-fn (fn [^double p]
            (cond
              (m/not-pos? p) 0
              (>= p 1.0) ##Inf
              :else (discrete-binary-search cdf-fn p [(long mmean) (long (m/sqrt variance))]))))

(defn- build-discrete
  [kind data probabilities]
  (let [cnt (count data)
        probabilities (or probabilities (repeat cnt 1))
        ^double sum (reduce m/fast+ probabilities)
        [emptymap upd corr] (if (= :int kind)
                              [(im/int-map) im/update unchecked-long]
                              [(sorted-map) update unchecked-double])
        pmf (reduce (fn [m [v ^double p]]
                      (upd m v (fnil m/fast+ 0.0) (/ p sum))) emptymap (map vector data probabilities))
        cumsum (reductions m/fast+ (vals pmf))
        ks (keys pmf)
        ^double mnk (first ks)
        icdf (let [stepf (i/step-before cumsum ks)]
               (fn [^double x]
                 (corr (stepf x))))
        cdf (let [stepf (i/step-after ks cumsum)]
              (fn [^double x]
                (if (< x mnk) 0.0 (stepf x))))]
    [pmf cdf icdf mnk (last ks)]))

(distribution-template :integer-discrete-distribution
                       {:mean mmean :variance @variance :distribution-parameters [:data :probabilities :rng]
                        :lower-bound lower-bound :upper-bound upper-bound :continuous? false :pdf? true}
                       {:keys [data probabilities]
                        :or {data [0]}} args
                       [pmf cdf-fn icdf-fn
                        ^long lower-bound ^long upper-bound] (build-discrete :int data probabilities)
                       pdf-fn (fn [^long k] (get pmf k 0.0))
                       ^double mmean (reduce-kv (fn [^double s ^long k ^double v]
                                                  (+ s (* k v))) 0.0 pmf)
                       variance (delay (- ^double (reduce-kv (fn [^double s ^long k ^double v]
                                                               (+ s (* k k v))) 0.0 pmf) (* mmean mmean))))

(distribution-template :real-discrete-distribution
                       {:mean mmean :variance @variance :distribution-parameters [:data :probabilities :rng]
                        :lower-bound lower-bound :upper-bound upper-bound :continuous? false :pdf? true}
                       {:keys [data probabilities]
                        :or {data [0]}} args
                       [pmf cdf-fn icdf-fn
                        ^double lower-bound ^double upper-bound] (build-discrete :double data probabilities)
                       pdf-fn (fn [^double k] (get pmf k 0.0))
                       ^double mmean (reduce-kv (fn [^double s ^double k ^double v]
                                                  (+ s (* k v))) 0.0 pmf)
                       variance (delay (- ^double (reduce-kv (fn [^double s ^double k ^double v]
                                                               (+ s (* k k v))) 0.0 pmf) (* mmean mmean))))

(defmethod distribution :categorical-distribution
  ([_ {:keys [data probabilities]
       :or {data [0]}
       :as all}]
   (let [r (or (:rng all) (rng :jvm))
         
         ^clojure.lang.ILookup unique (vec (distinct data))
         ^clojure.lang.ILookup dict (zipmap unique (range (count unique)))

         enumerated (distribution :integer-discrete-distribution
                                  {:data (map dict data) :probabilities probabilities :rng r})]
     (reify
       prot/DistributionProto
       (pdf [_ v] (prot/pdf enumerated (.valAt dict v -1)))
       (lpdf [_ v] (prot/lpdf enumerated (.valAt dict v -1)))
       (cdf [_ v] (prot/cdf enumerated (.valAt dict v -1)))
       (icdf [_ v] (.valAt unique (prot/icdf enumerated v)))
       (probability [_ v] (prot/probability enumerated (.valAt dict v -1)))
       (sample [_] (.valAt unique (prot/sample enumerated)))
       (dimensions [_] 1)
       (source-object [_] enumerated)
       (continuous? [_] false)
       prot/DistributionIdProto
       (distribution? [_] true)
       (distribution-id [_] :categorical-distribution)
       (distribution-parameters [_] [:data :probabilities :rng])
       prot/UnivariateDistributionProto
       (mean [_] ##NaN)
       (variance [_] ##NaN)
       prot/RNGProto
       (->seq [_] (map #(.valAt unique %) (prot/->seq enumerated)))
       (->seq [_ n] (map #(.valAt unique %) (prot/->seq enumerated n)))
       (set-seed! [d seed] (prot/set-seed! r seed) d))))
  ([_] (distribution :categorical-distribution {})))

(def ^{:const true :private true :tag 'double} LOG_M_2_PI (m/log m/M_2_PI))

(distribution-template :half-cauchy
                       {:mean ##NaN :variance ##NaN :distribution-parameters [:scale :rng]
                        :lower-bound 0.0 :upper-bound ##Inf}
                       {:keys [^double scale]
                        :or {scale 1.0}} args
                       ls (m/log scale)
                       lpdf-fn (fn [^double x]
                                 (if (neg? x)
                                   ##-Inf
                                   (- LOG_M_2_PI ls (m/log1p (m/sq (/ x scale))))))
                       icdf-fn (fn [^double p]
                                 (cond
                                   (m/not-pos? p) 0.0
                                   (>= p 1.0) 1.0
                                   :else (* scale (m/tan (* m/HALF_PI p)))))
                       cdf-fn (fn [^double v]
                                (* m/M_2_PI (m/atan (/ v scale)))))

(distribution-template :half-normal
                       {:mean mmean :variance variance
                        :distribution-parameters [:sigma :rng]
                        :lower-bound 0.0 :uppor-bound ##Inf}
                       {:keys [^double sigma] :or {sigma 1.0}} args
                       mmean (* sigma (/ m/SQRT2 m/M_SQRT_PI))
                       variance (* sigma sigma (- 1.0 (/ 2.0 m/PI)))
                       dist (distribution :normal {:mu 0.0 :sd sigma})
                       lpdf-fn (fn ^double [^double x]
                                 (if (neg? x)
                                   ##-Inf
                                   (+ m/M_LN2 (lpdf dist x))))
                       cdf-fn (fn ^double [^double x]
                                (if (neg? x)
                                  0.0
                                  (dec (* 2.0 (cdf dist x)))))
                       icdf-fn (fn ^double [^double x]
                                 (icdf dist (* 0.5 (inc x)))))

;; source: https://github.com/cran/gamlss.dist

(distribution-template :zaga
                       {:mean mmean :distribution-parameters [:mu :sigma :nu :lower-tail? :rng]
                        :lower-bound 0.0 :upper-bound ##Inf}
                       {:keys [^double mu ^double sigma ^double nu lower-tail?]
                        :or {mu 1.0 sigma 1.0 nu 0.1 lower-tail? true}} args
                       mmean (* (- 1.0 nu) mu)
                       s2 (* sigma sigma)
                       rs2 (/ s2)
                       lgrs2 (m/log-gamma rs2)
                       mus2 (* s2 mu)
                       rmus2 (/ mus2)
                       variance (* mmean mu (+ s2 nu))
                       lnu (m/log nu)
                       -nu (- 1.0 nu)
                       l1nu (m/log -nu)
                       gamma-dist (distribution :gamma (assoc args :rng r :shape rs2 :scale mus2))
                       lpdf-fn (fn [^double x]
                                 (if (zero? x)
                                   lnu
                                   (let [xx (* x rmus2)]
                                     (- (+ l1nu (* rs2 (m/log xx))) xx (m/log x) lgrs2))))
                       cdf-fn (fn [^double x]
                                (let [cdf (if (zero? x)
                                            nu
                                            (+ nu (* -nu (cdf gamma-dist x))))]
                                  (if lower-tail? cdf (- 1.0 cdf))))
                       icdf-fn (fn [^double x]
                                 (let [p (if lower-tail? x (- 1.0 x))
                                       p (if (<= p nu) nu p)]
                                   (prot/icdf gamma-dist (/ (- p nu) -nu)))))

(distribution-template :nbi
    {:mean mu
     :distribution-parameters [:mu :sigma :rng]
     :continuous? false
     :lower-bound 0.0 :upper-bound ##Inf}
  {:keys [^double mu ^double sigma]
   :or {mu 1.0 sigma 1.0}} args
  variance (+ mu (* sigma mu mu))
  distr (if (< sigma 0.0001)
          (distribution :poisson {:p mu :rng r})
          (let [nbinom-r (/ sigma)]
            (distribution :negative-binomial {:r nbinom-r
                                              :p (/ nbinom-r (+ nbinom-r mu))
                                              :rng r})))
  lpdf-fn (fn ^double [^double v] (prot/lpdf distr v))
  cdf-fn (fn ^double [^double v] (prot/cdf distr v))
  icdf-fn (fn ^double [^double v] (prot/icdf distr v)))

(distribution-template :zinbi
    {:mean mmean :distribution-parameters [:mu :sigma :nu :rng]
     :continuous? false :lower-bound 0.0 :upper-bound ##Inf}
  {:keys [^double mu ^double sigma ^double nu]
   :or {mu 1.0 sigma 1.0 nu 0.3}} args
  nu- (- 1.0 nu)
  lnu- (m/log nu-)
  mmean (* nu- mu)
  variance (+ mmean (* mmean mu (+ sigma nu)))
  distr (distribution :nbi {:mu mu :sigma sigma :rng r})
  lpdf-fn (fn ^double [^double x]
            (let [fy (lpdf distr x)]
              (if (zero? x)
                (m/log (+ nu (* nu- (m/exp fy))))
                (+ lnu- fy))))
  cdf-fn (fn ^double [^double x]
           (+ nu (* nu- (cdf distr x))))
  icdf-fn (fn ^double [^double x]
            (let [pnew (max 0.0 (- (/ (- x nu) nu-) 1.0e-7))]
              (prot/icdf distr pnew))))

(distribution-template :zanbi
    {:mean mmean :distribution-parameters [:mu :sigma :nu :rng]
     :continuous? false :lower-bound 0.0 :upper-bound ##Inf}
  {:keys [^double mu ^double sigma ^double nu]
   :or {mu 1.0 sigma 1.0 nu 0.3}} args
  lnu (m/log nu)
  nu- (- 1.0 nu)
  lnu- (m/log nu-)
  c (/ nu- (- 1.0 (m/pow (inc (* mu sigma)) (- (/ sigma)))))
  mmean (* mu c)
  variance (+ mmean (* mmean mu (inc (- sigma c))))
  distr (distribution :nbi {:mu mu :sigma sigma :rng r})
  lfy0 (- (m/log (- 1.0 (m/exp (lpdf distr 0.0)))))
  lpdf-fn (fn ^double [^double x]
            (let [fy (lpdf distr x)]                                   
              (if (zero? x) lnu (+ lnu- fy lfy0))))
  cdf0 (cdf distr 0.0)
  rcdf0- (/ (- 1.0 cdf0))
  cdf-fn (fn ^double [^double x]
           (if (zero? x)
             nu
             (+ nu (* nu- (- (cdf distr x) cdf0) rcdf0-))))
  icdf-fn (fn ^double [^double x]
            (let [pnew (- (/ (- x nu) nu-) 1.0e-10)
                  pnew2 (+ (* cdf0 (- 1.0 pnew)) pnew)]
              (prot/icdf distr (max 0.0 pnew2)))))

(distribution-template :zip
                       {:mean mmean :distribution-parameters [:mu :sigma :rng]
                        :continuous? false :lower-bound 0.0 :upper-bound ##Inf}
                       {:keys [^double mu ^double sigma]
                        :or {mu 5 sigma 0.1}} args
                       sigma- (- 1.0 sigma)
                       mmean (* sigma- mu)
                       variance (* mmean (inc (* mu sigma)))
                       lsigma-mu (- (m/log sigma-) mu)
                       lmu (m/log mu)
                       lpdf0 (m/log (+ sigma (* sigma- (m/exp (- mu)))))
                       lpdf-fn (fn ^double [^long x]
                                 (if (zero? x) lpdf0 (- (+ lsigma-mu (* x lmu))
                                                        (m/log-gamma (inc x)))))
                       dist (distribution :poisson {:p mu :rng r})
                       cdf-fn (fn ^double [^long x]
                                (+ sigma (* sigma- (cdf dist x))))
                       icdf-fn (fn ^double [^double x]
                                 (let [pnew (- (/ (- x sigma) sigma-) 1.0e-7)]
                                   (prot/icdf dist (max 0.0 pnew)))))

(distribution-template :zip2
                       {:mean mu
                        :distribution-parameters [:mu :sigma :rng]
                        :continuous? false :lower-bound 0.0 :upper-bound ##Inf}
                       {:keys [^double mu ^double sigma]
                        :or {mu 5.0 sigma 0.1}} args
                       sigma- (- 1.0 sigma)
                       lsigma- (m/log sigma-)
                       variance (* mu (inc (/ (* mu sigma) sigma-)))
                       mus (/ mu sigma-)
                       lmu (m/log mu)
                       lpdf0 (m/log (+ sigma (* sigma- (m/exp (- mus)))))
                       lpdf-fn (fn ^double [^long x]
                                 (if (zero? x) lpdf0 (+ (- (* (- 1.0 x) lsigma-) mus
                                                           (m/log-gamma (inc x)))
                                                        (* x lmu))))
                       dist (distribution :poisson {:p mus :rng r})
                       cdf-fn (fn ^double [^long x]
                                (+ sigma (* sigma- (cdf dist x))))
                       icdf-fn (fn ^double [^double x]
                                 (let [pnew (- (/ (- x sigma) sigma-) 1.0e-7)]
                                   (if (pos? pnew) (prot/icdf dist pnew) 0.0))))

(distribution-template :exgaus
                       {:distribution-parameters [:mu :sigma :nu :rng]
                        :continuous? true :lower-bound ##-Inf :upper-bound ##Inf}
                       {:keys [^double mu ^double sigma ^double nu]
                        :or {mu 0.0 sigma 1.0 nu 1.0}} args
                       mean (+ mu nu)
                       sigma2 (* sigma sigma)
                       variance (+ sigma2 (* nu nu))
                       -lnu (- (m/log nu))
                       sigma2nu (/ sigma2 nu)
                       dist (distribution :normal {:mu mu :sd sigma :rng r})
                       ndist (distribution :normal {:rng r})
                       lpdf-fn (if (> nu (* sigma 0.05))
                                 (fn ^double [^double x]
                                   (let [z (- x mu sigma2nu)]                                   
                                     (+ (- -lnu (/ (+ z (* 0.5 sigma2nu)) nu))
                                        (m/log (cdf ndist (/ z sigma))))))
                                 (fn ^double [^double x] (prot/lpdf dist x)))
                       cdf-fn (if (> nu (* sigma 0.05))
                                (let [exppart (- (m/sq (+ mu sigma2nu)) (* mu mu))]
                                  (fn ^double [^double q]
                                    (let [z (- q mu sigma2nu)
                                          pnorm1 (cdf ndist (/ (- q mu) sigma))
                                          pnorm2 (cdf ndist (/ z sigma))]
                                      (- pnorm1 (* pnorm2 (m/exp (/ (- exppart (* 2.0 q sigma2nu))
                                                                    (* 2.0 sigma2))))))))
                                (fn ^double [^double q]
                                  (prot/cdf dist q)))
                       ^double hmu (cdf-fn mu)
                       icdf-fn (fn ^double [^double p]
                                 (let [h1 (fn ^double [^double q] (- ^double (cdf-fn q) p))]
                                   (if (< hmu p)
                                     (loop [interval (+ mu sigma)
                                            j 2]
                                       (if (< ^double (cdf-fn interval) p)
                                         (recur (+ mu (* j sigma)) (inc j))
                                         (solver/find-root h1 mu interval)))
                                     (loop [interval (- mu sigma)
                                            j 2]
                                       (if (> ^double (cdf-fn interval) p)
                                         (recur (- mu (* j sigma)) (inc j))
                                         (solver/find-root h1 interval mu)))))))

(distribution-template :bb
                       {:mean mmean :distribution-parameters [:mu :sigma :bd :rng]
                        :continuous? false :lower-bound 0.0 :upper-bound bd}
                       {:keys [^double mu ^double sigma ^long bd]
                        :or {mu 0.5 sigma 1.0 bd 10}} args
                       mmean (* bd mu)
                       variance (* mmean (- 1.0 mu) (inc (/ (* sigma (dec bd))
                                                            (inc sigma))))
                       dist (distribution :binomial {:p mu :trials bd :rng r})
                       lpdf-fn (if (< sigma 0.00001)
                                 (fn ^double [^long x] (prot/lpdf dist x))
                                 (let [rsigma (/ sigma)
                                       mursigma (* mu rsigma)
                                       mu-rsigma (* (- 1.0 mu) rsigma)
                                       lgamma-part (- (+ (m/log-gamma (inc bd))
                                                         (m/log-gamma rsigma))
                                                      (m/log-gamma mursigma)
                                                      (m/log-gamma mu-rsigma)
                                                      (m/log-gamma (+ bd rsigma)))]
                                   (fn ^double [^long x]
                                     (+ (- lgamma-part
                                           (m/log-gamma (inc x))
                                           (m/log-gamma (inc (- bd x))))
                                        (m/log-gamma (+ x mursigma))
                                        (m/log-gamma (- (+ bd mu-rsigma) x))))))
                       cdf-fn (if (< sigma 0.00001)
                                (fn ^double [^long x] (prot/cdf dist x))
                                (memoize (fn ^double [^long q]
                                           (reduce m/fast+ (map #(m/exp (lpdf-fn %)) (range (inc (long q))))))))
                       icdf-fn (if (< sigma 0.00001)
                                 (fn ^double [^double p] (prot/icdf dist p))
                                 (let [r (range 0 (inc bd))]
                                   (i/step-before (rest (reductions
                                                         (fn [^double s ^double v]
                                                           (+ s (m/exp (lpdf-fn v)))) 0.0 r)) r))))


(distribution-template :zabi
                       {:mean mmean
                        :distribution-parameters [:mu :sigma :bd :rng]
                        :continuous? false :lower-bound 0.0 :upper-bound bd}
                       {:keys [^double mu ^double sigma ^long bd]
                        :or {mu 0.5 sigma 0.1 bd 1}} args
                       sigma- (- 1.0 sigma)
                       mmean (/ (* sigma- bd mu)
                                (- 1.0 (m/pow (- 1.0 mu) bd)))
                       variance (- (* mmean (+ (- 1.0 mu) (* bd mu))) (* mmean mmean))
                       lsigma (m/log sigma)
                       lsigma- (m/log sigma-)
                       dist (distribution :binomial {:trials bd :p mu :rng r})
                       lpdf0- (m/log (- 1.0 (pdf dist 0.0)))
                       lpdf-fn (fn ^double [^double x]
                                 (if (zero? x)
                                   lsigma
                                   (- (+ lsigma- (lpdf dist x)) lpdf0-)))
                       cdf2 (cdf dist 0.0)
                       rcdf2- (/ (- 1.0 cdf2))
                       cdf-fn (fn ^double [^double q]
                                (if (zero? q)
                                  sigma
                                  (let [cdf1 (cdf dist q)]
                                    (+ sigma (* sigma- (- cdf1 cdf2) rcdf2-)))))
                       icdf-fn (fn ^double [^double p]
                                 (let [pnew (- (/ (- p sigma) sigma-) 1.0e-10)]
                                   (if (pos? pnew)
                                     (prot/icdf dist (+ (* cdf2 (- 1.0 pnew)) pnew))
                                     0.0))))

(distribution-template :zibi
                       {:mean mmean
                        :distribution-parameters [:mu :sigma :bd :rng]
                        :continuous? false :lower-bound 0.0 :upper-bound bd}
                       {:keys [^double mu ^double sigma ^long bd]
                        :or {mu 0.5 sigma 0.1 bd 1}} args
                       sigma- (- 1.0 sigma)
                       lsigma- (m/log sigma-)
                       mmean (* sigma- bd mu)
                       variance (* mmean (+ (- 1.0 mu) (* sigma bd mu)))
                       dist (distribution :binomial {:trials bd :p mu :rng r})
                       pdf0 (m/log (+ sigma (* sigma- (pdf dist 0.0))))
                       lpdf-fn (fn ^double [^double x]
                                 (if (zero? x)
                                   pdf0
                                   (+ lsigma- (lpdf dist x))))
                       cdf-fn (fn ^double [^double q]
                                (+ sigma (* sigma- (cdf dist q))))
                       icdf-fn (fn ^double [^double p]
                                 (let [pnew (- (/ (- p sigma) sigma-) 1.0e-10)]
                                   (if (pos? pnew) (prot/icdf dist pnew) 0.0))))

;; mean and variance from the paper: https://www.gamlss.com/wp-content/uploads/2018/01/DistributionsForModellingLocationScaleandShape.pdf

(distribution-template :zabb
                       {:mean mmean :distribution-parameters [:mu :sigma :bd :nu :rng]
                        :continuous? false :lower-bound 0.0 :upper-bound bd}
                       {:keys [^double mu ^double sigma ^double bd ^double nu]
                        :or {mu 0.5 sigma 0.1 nu 0.1 bd 1.0}} args
                       lnu (m/log nu)
                       nu- (- 1.0 nu)
                       lnu- (m/log nu-)
                       dist (distribution :bb {:mu mu :sigma sigma :bd bd :rng r})
                       pdf0- (- 1.0 (pdf dist 0.0))
                       mmean (/ (* nu- mu bd) pdf0-)
                       variance (- (/ (* nu- (+ (* bd mu (- 1.0 mu) (inc (/ (* sigma (dec bd)) (inc sigma))))
                                                (* bd bd mu mu))) pdf0-)
                                   (* mmean mmean))
                       -lpdf0- (- (m/log pdf0-))
                       cdf2 (cdf dist 0.0)
                       rcdf2- (/ (- 1.0 cdf2))
                       lpdf-fn (fn ^double [^double x]
                                 (if (zero? x)
                                   lnu
                                   (+ lnu- (lpdf dist x) -lpdf0-)))
                       cdf-fn (fn ^double [^double q]
                                (if (zero? q)
                                  nu
                                  (min 1.0 (+ nu (* nu- (- (cdf dist q) cdf2) rcdf2-)))))
                       icdf-fn (fn ^double [^double p]
                                 (let [pnew (max 0.0 (- (/ (- p nu) nu-)1.0e-7))]
                                   (if (pos? pnew)
                                     (prot/icdf dist (+ (* cdf2 (- 1.0 pnew)) pnew)) 0.0))))

(distribution-template :zibb
                       {:mean mmean
                        :distribution-parameters [:mu :sigma :bd :nu :rng]
                        :continuous? false :lower-bound 0.0 :upper-bound bd}
                       {:keys [^double mu ^double sigma ^double bd ^double nu]
                        :or {mu 0.5 sigma 0.5 nu 0.1 bd 1.0}} args
                       lnu (m/log nu)
                       nu- (- 1.0 nu)
                       lnu- (m/log nu-)
                       mmean (* nu- bd mu)
                       variance (+ (* mmean (- 1.0 mu) (inc (/ (* sigma (dec bd)) (inc sigma))))
                                   (* nu nu- bd bd mu mu))
                       dist (distribution :bb {:mu mu :sigma sigma :bd bd :rng r})
                       lpdf0- (m/log (+ nu (* nu- (pdf dist 0.0))))
                       lpdf-fn (fn ^double [^double x]
                                 (if (zero? x) lpdf0- (+ lnu- (lpdf dist x))))
                       cdf-fn (fn ^double [^double q]
                                (min 1.0 (+ nu (* nu- (cdf dist q)))))
                       icdf-fn (fn ^double [^double p]
                                 (let [pnew (max 0.0 (- (/ (- p nu) nu-)1.0e-7))]
                                   (if (pos? pnew)
                                     (prot/icdf dist pnew) 0.0))))

(defonce ^{:doc "Default normal distribution (u=0.0, sigma=1.0)."} default-normal (distribution :normal))

(distribution-template :truncated
                       {:mean ##NaN :variance ##NaN
                        :distribution-parameters [:distr :left :right :rng]
                        :continuous? (continuous? distribution)
                        :lower-bound left-bound :upper-bound right-bound}
                       {:keys [distr left right]
                        :or {distr default-normal}} args
                       left-cdf (if left (cdf distr left) 0.0)
                       right-cdf (if right (cdf distr right) 1.0)
                       cdf-diff (- right-cdf left-cdf)
                       lcdf-diff (m/log cdf-diff)
                       left-bound (or left (lower-bound distr))
                       right-bound (or right (upper-bound distr))
                       ^double left (or left ##-Inf)
                       ^double right (or right ##Inf)
                       lpdf-fn (fn ^double [^double x]
                                 (if (and (<= left x) (<= x right))
                                   (- (lpdf distr x) lcdf-diff)
                                   ##-Inf))
                       cdf-fn (fn ^double [^double x]
                                (m/constrain (/ (- (cdf distr x) left-cdf) cdf-diff) 0.0 1.0))
                       icdf-fn (fn ^double [^double x]
                                 (icdf distr (+ left-cdf (* x cdf-diff)))))

(distribution-template :mixture
                       {:pdf? true :distribution-parameters [:distrs :weights :rng]
                        :mean mean-val
                        :continuous? continuous? :lower-bound lower-bound :upper-bound upper-bound}
                       {:keys [distrs weights]
                        :or {distrs [default-normal]}} args
                       cnt (count distrs)
                       weights (vec (or weights (repeat cnt 1.0)))
                       weights (v/div weights (v/sum weights))
                       continuous? (continuous? (first distrs))
                       lower-bound (reduce m/fast-min (map lower-bound distrs))
                       upper-bound (reduce m/fast-max (map upper-bound distrs))
                       mean-val (reduce m/fast+ (map (fn ^double [^double w d]
                                                       (* w (mean d))) weights distrs))
                       variance (- ^double (reduce m/fast+ (map (fn ^double [^double w d]
                                                                  (* w (+ (m/sq (mean d))
                                                                          (variance d)))) weights distrs))
                                   (m/sq mean-val))
                       pdf-fn (fn ^double [^double x]
                                (reduce m/fast+ (map (fn ^double [^double w d]
                                                       (* w (pdf d x))) weights distrs)))
                       cdf-fn (fn ^double [^double x]
                                (reduce m/fast+ (map (fn ^double [^double w d]
                                                       (* w (cdf d x))) weights distrs)))
                       icdf-fn (fn [^double x]
                                 (let [icdfs (map #(icdf % x) distrs)
                                       mn (reduce m/fast-min icdfs)
                                       mx (reduce m/fast-max icdfs)
                                       target-fn (fn ^double [^double v] (- ^double (cdf-fn v) x))]
                                   (solver/find-root target-fn mn mx))))

;; from Julia
(distribution-template :kolmogorov
                       {:pdf? true :distribution-parameters [:rng]
                        :mean 0.8687311606361591 :variance 0.0677732039638651
                        :continuous? true :lower-bound 0.0 :upper-bound ##Inf}
                       cdf-raw (fn ^double [^double x]
                                 (let [a (- (m/sq (/ m/PI x)))
                                       f (m/exp a)
                                       f2 (* f f)
                                       u (inc (* f (inc f2)))]
                                   (/ (* m/SQRT2PI (m/exp (* 0.125 a)) u) x)))
                       ccdf-raw (fn [^double ^double x]
                                  (let [f (m/exp (* -2.0 x x))
                                        f2 (* f f)
                                        f3 (* f f2)
                                        f5 (* f2 f3)
                                        f7 (* f2 f5)
                                        u (- 1.0 (* f3 (- 1.0 (* f5 (- 1.0 f7)))))]
                                    (* 2.0 f u)))
                       cdf-fn (fn ^double [^double x]
                                (cond
                                  (not (pos? x)) 0.0
                                  (<= x 1.0) (cdf-raw x)
                                  :else (- 1.0 ^double (ccdf-raw x))))
                       pdf-fn (fn ^double [^double x]
                                (cond
                                  (not (pos? x)) 0.0
                                  (<= x 1.0) (let [c (/ m/PI (* 2.0 x))
                                                   ks (map (fn [^long i]
                                                             (let [k (m/sq (* i c))]
                                                               (* (dec k) (m/exp (* -0.5 k)))))
                                                           (range 1 40 2))
                                                   ^double s (reduce m/fast+ ks)]
                                               (/ (* m/SQRT2PI s) (* x x)))
                                  :else (let [ks (map (fn [^double a ^long i]
                                                        (* a i i (m/exp (* -2.0 (m/sq (* i x))))))
                                                      (cycle [1.0 -1.0]) (range 1 21))
                                              ^double s (reduce m/fast+ ks)]
                                          (* 8.0 x s))))
                       icdf-fn (fn ^double [^double p]
                                 (let [h1 (fn ^double [^double q] (- ^double (cdf-fn q) p))]
                                   (if (< 0.5626816957524641 p) ;; cdf(mean)
                                     (loop [interval 1.1290640320985719 ;; mean + sigma
                                            j 2]
                                       (if (< ^double (cdf-fn interval) p)
                                         (recur (+ 0.8687311606361591 (* j 0.26033287146241274)) (inc j))
                                         (solver/find-root h1 0.8687311606361591 interval)))
                                     (loop [interval 0.6083982891737463 ;; mean - sigma
                                            j 2]
                                       (if (> ^double (cdf-fn interval) p)
                                         (recur (- 0.8687311606361591 (* j 0.26033287146241274)) (inc j))
                                         (solver/find-root h1 interval 0.8687311606361591)))))))

;; from Julia
(distribution-template :fishers-noncentral-hypergeometric
                       {:pdf? true :distribution-parameters [:rng :ns :nf :n :omega]
                        :continuous? false :lower-bound mlower-bound :upper-bound mupper-bound
                        :mean mmean}
                       {:keys [^long ns ^long nf ^long n ^double omega]
                        :or {ns 10 nf 10 n 5 omega 1.0}} args
                       mlower-bound (max 0 (- n nf))
                       mupper-bound (min ns n)
                       mode (let [A (dec omega)
                                  B (- n nf (* (+ ns n 2) omega))
                                  C (* (inc ns) (* (inc n)) omega)]
                              (long (m/floor (/ (* -2.0 C)
                                                (- B (m/sqrt (- (* B B) (* 4.0 A C))))))))
                       fri (fn [^long i]
                             (* (/ (* (inc (- ns i)) omega)
                                   (* i (+ (- nf n) i)))
                                (inc (- n i))))
                       fri+ (fn [^long i]
                              (* (/ (* (- ns i) omega)
                                    (* (inc i)
                                       (inc (+ (- nf n) i))))
                                 (- n i)))
                       expectation (fn [f]
                                     (let [[^double s ^double m] (loop [m (double (f mode))
                                                                        fi 1.0
                                                                        s 1.0
                                                                        i (inc mode)]
                                                                   (if (> i mupper-bound)
                                                                     [s m]
                                                                     (let [^double ri (fri i)
                                                                           nfi (* fi ri)
                                                                           sfi (+ s nfi)]
                                                                       (if (== sfi s)
                                                                         [s m]
                                                                         (recur (+ m (* ^double (f i) nfi))
                                                                                nfi sfi (inc i))))))]
                                       (loop [m m
                                              fi 1.0
                                              s s
                                              i (dec mode)]
                                         (if (< i mlower-bound)
                                           (/ m s)
                                           (let [^double ri (fri+ i) 
                                                 nfi (/ fi ri)
                                                 sfi (+ s nfi)]
                                             (if (== sfi s)
                                               (/ m s)
                                               (recur (+ m (* ^double (f i) nfi))
                                                      nfi sfi (dec i))))))))
                       ^double mmean (expectation identity)
                       variance (expectation (fn [^double t] (m/sq (- t mmean))))
                       pdf-fn (fn [^double k]
                                (if-not (m/between? mlower-bound mupper-bound k)
                                  0.0
                                  (let [k (int k)
                                        [^double s ^double fk] (loop [fk 1.0
                                                                      fi 1.0
                                                                      s 1.0
                                                                      i (inc mode)]
                                                                 (if (> i mupper-bound)
                                                                   [s fk]
                                                                   (let [^double ri (fri i)
                                                                         nfi (* fi ri)
                                                                         sfi (+ s nfi)]
                                                                     (if (and (== sfi s)
                                                                              (> i k))
                                                                       [s fk]
                                                                       (recur (if (== k i) nfi fk)
                                                                              nfi sfi (inc i))))))]
                                    (loop [fk fk
                                           fi 1.0
                                           s s
                                           i (dec mode)]
                                      (if (< i mlower-bound)
                                        (/ fk s)
                                        (let [^double ri (fri+ i) 
                                              nfi (/ fi ri)
                                              sfi (+ s nfi)]
                                          (if (and (== sfi s)
                                                   (< i k))
                                            (/ fk s)
                                            (recur (if (== k i) nfi fk)
                                                   nfi sfi (dec i)))))))))
                       cdf-fn (fn [^double k]                                
                                (cond
                                  (< k mlower-bound) 0.0
                                  (>= k mupper-bound) 1.0
                                  :else (let [k (int k)
                                              [^double s ^double fk] (loop [fk (if (>= k mode) 1.0 0.0)
                                                                            fi 1.0
                                                                            s 1.0
                                                                            i (inc mode)]
                                                                       (if (> i mupper-bound)
                                                                         [s fk]
                                                                         (let [^double ri (fri i)
                                                                               nfi (* fi ri)
                                                                               sfi (+ s nfi)]
                                                                           (if (and (== sfi s)
                                                                                    (> i k))
                                                                             [s fk]
                                                                             (recur (if (<= i k)
                                                                                      (+ fk nfi) fk)
                                                                                    nfi sfi (inc i))))))]
                                          (loop [fk fk
                                                 fi 1.0
                                                 s s
                                                 i (dec mode)]
                                            (if (< i mlower-bound)
                                              (/ fk s)
                                              (let [^double ri (fri+ i) 
                                                    nfi (/ fi ri)
                                                    sfi (+ s nfi)]
                                                (if (and (== sfi s)
                                                         (< i k))
                                                  (/ fk s)
                                                  (recur (if (<= i k) (+ fk nfi) fk)
                                                         nfi sfi (dec i)))))))))
                       icdf-fn (fn [^double q]
                                 (if-not (<= 0.0 q 1.0)
                                   ##NaN
                                   (loop [i mlower-bound]
                                     (cond
                                       (> i mupper-bound) mupper-bound
                                       (> q ^double (cdf-fn i)) (recur (inc i)) 
                                       :else i)))))

(defonce ^{:doc "List of distributions."}
  distributions-list
  (into (sorted-set) (keys (methods distribution))))
;;

(defn set-seed
  "Create and return new RNG"
  ([]
   (prot/set-seed default-rng (lrand)))
  ([^long v]
   (prot/set-seed default-rng v))
  ([rng ^long v]
   (prot/set-seed rng v)))

(defn set-seed!
  "Sets seed.

  If `rng` is `:smile` calls `smile.math.MathEx/setSeed()`.

  Without `rng` sets both `:smile` and `default-rng`"
  ([]
   (MathEx/setSeed (lrand))
   (prot/set-seed! default-rng (lrand)))
  ([^long v]
   (MathEx/setSeed v)
   (prot/set-seed! default-rng v))
  ([rng ^long v]
   (if (= rng :smile)
     (MathEx/setSeed v)
     (prot/set-seed! rng v))))

;;

(defn- uniform-spacings
  ([^long n] (uniform-spacings default-rng n))
  ([rng ^long n]
   (let [xs (reductions m/fast+ (repeatedly (inc n) #(- (m/log (drandom rng)))))
         l (/ ^double (last xs))]
     (map (fn [^double x] (* x l)) (butlast xs)))))

(defn- systematic-spacings
  ([^long n] (systematic-spacings default-rng n))
  ([rng ^long n]
   (let [l (/ 1.0 n)
         d (drandom rng)]
     (map (fn [^long x] (* (+ x d) l)) (range n)))))

(defn- stratified-spacings
  ([^long n] (systematic-spacings default-rng n))
  ([rng ^long n]
   (let [l (/ 1.0 n)]
     (map (fn [^long x] (* (+ x (drandom rng)) l)) (range n)))))

(defn- antithetic-sampling
  ([^long n] (antithetic-sampling default-rng n))
  ([rng ^long n]
   (->> (repeatedly (fn [] (let [r1 (drandom rng)]
                            [r1 (- 1.0 r1)])))
        (mapcat identity)
        (take n))))

(defn- jittered-sequence-sampling
  ([kind ^long n] (jittered-sequence-sampling kind nil n))
  ([kind _ ^long n]
   (take n (jittered-sequence-generator kind 1))))


(def ^:private spacings
  {:uniform uniform-spacings
   :systematic systematic-spacings
   :stratified stratified-spacings
   :antithetic antithetic-sampling
   :r2 (partial jittered-sequence-sampling :r2)
   :sobol (partial jittered-sequence-sampling :sobol)
   :halton (partial jittered-sequence-sampling :halton)})

(defn ->seq
  "Returns lazy sequence of random samples (can be limited to optional `n` values).

  Additionally one of the sampling methods can be provided, ie: `:uniform`, `:antithetic`, `:systematic` and `:stratified`."
  ([] (prot/->seq default-rng))
  ([rng] (prot/->seq rng))
  ([rng n] (prot/->seq rng n))
  ([rng n sampling-method]
   (if-not sampling-method
     (->seq rng n)
     (if (distribution? rng)
       (map (partial icdf rng) ((spacings sampling-method uniform-spacings) n))
       ((spacings sampling-method uniform-spacings) rng n)))))

(m/unuse-primitive-operators)
