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

  For random noise generation you can use [[random-noise-cfg]] and [[random-noise]]. Both can be feed with configuration. Additional configuration:

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
            [fastmath.matrix :as mat]
            [fastmath.protocols :as prot]
            [fastmath.special :as special]
            [fastmath.random.distributions :as distr]
            [fastmath.solver :as solver]
            [fastmath.interpolation.step :as step-interp]
            [fastmath.stats.bins :as bins])
  (:import [org.apache.commons.math3.random RandomGenerator ISAACRandom JDKRandomGenerator MersenneTwister
            Well512a Well1024a Well19937a Well19937c Well44497a Well44497b
            RandomVectorGenerator HaltonSequenceGenerator SobolSequenceGenerator UnitSphereRandomVectorGenerator
            EmpiricalDistribution SynchronizedRandomGenerator]
           [fastmath.java R2]
           [umontreal.ssj.probdist AndersonDarlingDist AndersonDarlingDistQuick BetaSymmetricalDist
            InverseGammaDist  ChiDist ChiSquareNoncentralDist CramerVonMisesDist ErlangDist FatigueLifeDist FoldedNormalDist FrechetDist HalfNormalDist HyperbolicSecantDist InverseGaussianDist HypoExponentialDist HypoExponentialDistEqual JohnsonSBDist JohnsonSLDist JohnsonSUDist KolmogorovSmirnovDist KolmogorovSmirnovDistQuick KolmogorovSmirnovPlusDist LoglogisticDist Pearson6Dist PowerDist RayleighDist WatsonGDist WatsonUDist]
           [fastmath.java.noise Billow RidgedMulti FBM NoiseConfig Noise Discrete]
           [org.apache.commons.math3.distribution BetaDistribution CauchyDistribution ChiSquaredDistribution ConstantRealDistribution EnumeratedRealDistribution ExponentialDistribution FDistribution GammaDistribution, GumbelDistribution, LaplaceDistribution, LevyDistribution, LogisticDistribution, LogNormalDistribution, NakagamiDistribution, NormalDistribution, ParetoDistribution, TDistribution, TriangularDistribution, UniformRealDistribution WeibullDistribution MultivariateNormalDistribution]
           [org.apache.commons.math3.distribution BinomialDistribution EnumeratedIntegerDistribution, GeometricDistribution, HypergeometricDistribution, PascalDistribution, PoissonDistribution, UniformIntegerDistribution, ZipfDistribution]))


(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)
(m/use-primitive-operators)

;; Helper macro which creates RNG object of given class and/or seed.
(defmacro ^:private create-object-with-seed
  "Create object of the class with (or not) given seed. Used to create RNG."
  [cl seed]
  `(if-let [arg# ~seed]
     (prot/set-seed! (new ~cl) (long arg#)) ;; seeding via protocols
     (new ~cl)))

(defmulti rng
  "Creates a random number generator, from Apache Commons Math, for the given algorithm.

  Dispatches on `rng-name` to construct the underlying generator object; see [[rngs-list]] for the complete list of registered names. Use [[synced-rng]] to get a thread-safe wrapper around one of these generators.

  Parameters:

  - `rng-name` (keyword): algorithm to use:
      - `:jdk` - `java.util.Random` (`JDKRandomGenerator`).
      - `:mersenne` - Mersenne Twister.
      - `:isaac` - ISAAC.
      - `:well512a`, `:well1024a`, `:well19937a`, `:well19937c`, `:well44497a`, `:well44497b` - WELL generator variants.
      - `:default` - alias for `:jdk`.
  - `seed` (optional, long): initial seed. When given, the generator produces a fully reproducible sequence of values for a given `rng-name`; when omitted, the generator seeds itself (implementation-dependent, typically from system entropy or clock).

  Returns a new, mutable generator object implementing [[RNGProto]] (`irandom`, `lrandom`, `frandom`, `drandom`, `grandom`, `brandom`, `set-seed`, `set-seed!`, `->seq`), usable wherever an `rng` parameter is expected throughout the namespace, for example `(irandom (rng :isaac 1337))`.

  See also [[rngs-list]], [[synced-rng]], [[default-rng]]."
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
   (reduce m/+ (repeatedly dices #(inc (irand sides)))))  )

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
   (reduce m/+ (repeatedly dices #(inc (irandom rng sides))))))

;; generators

;; http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/#more-2165
(defn ball-random
  "Draws a uniformly random point from the interior of a `dims`-dimensional unit ball.

  Uses the method described in this [article](http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/#more-2165): draws `dims + 2` independent `N(0,1)` samples, normalizes them to a unit vector (a uniformly random point on the `(dims + 2)`-dimensional unit sphere), and keeps only its first `dims` coordinates. Projecting a uniform sample from a sphere two dimensions higher down onto `dims` coordinates yields a point uniformly distributed by volume within the `dims`-ball, unlike naively rejecting or normalizing `dims`-dimensional samples, which biases the result.

  Parameters:

  - `rng` (optional): random number generator to draw the underlying gaussian samples from. Default: [[default-rng]].
  - `dims` (long): dimensionality of the ball.

  Returns a double when `dims` is `1`, a `Vec2`, `Vec3` or `Vec4` when `dims` is `2`, `3` or `4`, or a plain vector of doubles otherwise; magnitude of the result is always less than or equal to `1.0`.

  Used internally by the `:ball` method of [[sequence-generator]]."
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
  ^{:doc "Creates a lazy, infinite sequence of random or quasi-random points.

  Dispatches on `seq-generator` to pick the sampling method; every point in the resulting sequence has `dimensions` components. See [[sequence-generators-list]] for the full list of registered `seq-generator` keys, and [[jittered-sequence-generator]] to add blue-noise jitter to `:r2`, `:halton` and `:sobol` sequences.

  Parameters:

  - `seq-generator`: keyword selecting the generator:
      - `:r2`, `:halton`, `:sobol` - low-discrepancy (quasi-random) sequences filling `[0,1]` for each dimension more evenly than pseudo-random sampling.
      - `:default`, `:uniform` or any other unregistered keyword - independent uniform pseudo-random points from `[0,1]` for each dimension.
      - `:gaussian` - independent pseudo-random points, each component drawn from `N(0,1)`.
      - `:sphere` - pseudo-random points on the surface of a unit sphere (euclidean distance from origin equals `1.0`).
      - `:ball` - pseudo-random points uniformly distributed within a unit ball.
  - `dimensions` (long): number of components per point. Limited to `1-15` for `:r2`, `1-40` for `:halton` and `1-1000` for `:sobol`; unrestricted (`1+`) for the other generators.

  Returns a lazy, infinite sequence of points: a double when `dimensions` is `1`, a `Vec2`, `Vec3` or `Vec4` when `dimensions` is `2`, `3` or `4`, or a plain vector of doubles otherwise.

  Throws an assertion error when `dimensions` exceeds the allowed range for `:r2`, `:halton` or `:sobol`.

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
  "Creates a lazy, infinite sequence of jittered [[sequence-generator]] points.

  Perturbs each point of the underlying `seq-generator` sequence by a small amount of noise, according to this [article](http://extremelearning.com.au/a-simple-method-to-construct-isotropic-quasirandom-blue-noise-point-sequences/), breaking up the perfectly regular structure of low-discrepancy sequences while keeping their approximately uniform coverage, closer to a blue-noise point distribution. Two jittering strategies are used depending on `seq-generator`:

  - `:sphere` and `:gaussian` - each coordinate is offset by independent `N(0,1)` noise scaled by `jitter`; the result is not renormalized back onto the sphere or wrapped, so points may drift away from it.
  - every other generator (including `:r2`, `:halton`, `:sobol`, `:default`/`:uniform` and `:ball`) - each coordinate is offset by scaled quasi-random noise (tuned per `seq-generator`, falling back to generic constants for unlisted ones) and wrapped back into `[0,1]` with `frac`.

  Parameters:

  - `seq-generator`: keyword selecting the base generator, same as for [[sequence-generator]]. Intended for `:r2`, `:sobol` and `:halton`, whose evenly spaced points benefit most from jittering, but works with any registered generator.
  - `dimensions` (long): number of components per point, passed through to [[sequence-generator]].
  - `jitter` (double, optional): jitter amount, from `0.0` (no jitter, identical to the unjittered [[sequence-generator]] sequence) to `1.0` (full jitter). Default: `0.25`.

  Returns a lazy, infinite sequence of points, same shape as [[sequence-generator]] for the given `dimensions` (a double, `Vec2`, `Vec3`, `Vec4` or a plain vector).

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
  "Value noise, 6-octave FBM with Hermite interpolation.

  A ready-to-use, zero-argument preset of [[fbm-noise]] with `:noise-type` `:value`, `:interpolation` `:hermite`, `:octaves` `6`, `:lacunarity` `2.0` and `:gain` `0.5`. The underlying noise interpolates randomly assigned values at integer lattice points, giving a blockier, less directional look than gradient-based [[noise]]. The `:seed` is fixed once when the namespace loads, so repeated calls with the same arguments always return the same value.

  Accepts 1, 2 or 3 double arguments (`x`, `x y` or `x y z`) and returns a double from the `[0,1]` range.

  See also [[noise]], [[simplex]], [[single-noise]], [[fbm-noise]], [[billow-noise]], [[ridgedmulti-noise]]."
  (^double [^double x] (FBM/noise value-noise-config x))
  (^double [^double x ^double y] (FBM/noise value-noise-config x y))
  (^double [^double x ^double y ^double z] (FBM/noise value-noise-config x y z)))

(defn noise
  "Improved Perlin noise, 6-octave FBM with quintic interpolation.

  A ready-to-use, zero-argument preset of [[fbm-noise]] with `:noise-type` `:gradient`, `:interpolation` `:quintic`, `:octaves` `6`, `:lacunarity` `2.0` and `:gain` `0.5`. Gradient noise interpolates dot products of pseudo-random gradient vectors at lattice points, and the quintic (`6t^5 - 15t^4 + 10t^3`) interpolation removes second-derivative discontinuities at cell boundaries, giving the smoother look of Ken Perlin's improved noise. The `:seed` is fixed once when the namespace loads, so repeated calls with the same arguments always return the same value.

  Accepts 1, 2 or 3 double arguments (`x`, `x y` or `x y z`) and returns a double from the `[0,1]` range.

  See also [[vnoise]], [[simplex]], [[single-noise]], [[fbm-noise]], [[billow-noise]], [[ridgedmulti-noise]]."
  (^double [^double x] (FBM/noise perlin-noise-config x))
  (^double [^double x ^double y] (FBM/noise perlin-noise-config x y))
  (^double [^double x ^double y ^double z] (FBM/noise perlin-noise-config x y z)))

(defn simplex
  "Simplex noise, 6-octave FBM.

  A ready-to-use, zero-argument preset of [[fbm-noise]] with `:noise-type` `:simplex`, `:octaves` `6`, `:lacunarity` `2.0` and `:gain` `0.5`. Simplex noise evaluates gradients on a simplectic (triangular/tetrahedral) lattice rather than a square/cubic grid, which reduces directional artifacts and scales better to higher dimensions than gradient noise. The `:seed` is fixed once when the namespace loads, so repeated calls with the same arguments always return the same value.

  Accepts 1, 2 or 3 double arguments (`x`, `x y` or `x y z`) and returns a double from the `[0,1]` range.

  See also [[noise]], [[vnoise]], [[single-noise]], [[fbm-noise]], [[billow-noise]], [[ridgedmulti-noise]]."
  (^double [^double x] (FBM/noise simplex-noise-config x))
  (^double [^double x ^double y] (FBM/noise simplex-noise-config x y))
  (^double [^double x ^double y ^double z] (FBM/noise simplex-noise-config x y z)))

(defn single-noise
  "Creates a single-octave noise function.

  Produces raw, unblended noise: one evaluation of the underlying noise type (value, gradient or simplex, selected by `:noise-type`) per point, without combining multiple octaves. This is the base building block used by [[fbm-noise]], [[billow-noise]] and [[ridgedmulti-noise]] to build multi-octave noise.

  Parameters:

  - `cfg` (optional, map): noise configuration, all keys optional:
      - `:seed` - long seed for the noise's internal RNG, default: random.
      - `:noise-type` - `:value`, `:gradient` or `:simplex`, default: `:gradient`.
      - `:interpolation` - `:none`, `:linear`, `:hermite` or `:quintic`, used only by `:value` and `:gradient` noise types, default: `:hermite`.
      - `:normalize?` - normalize result to `[0,1]` range (`true`, default) or leave it in `[-1,1]` (`false`).
      - `:octaves`, `:lacunarity`, `:gain` - accepted for consistency with [[fbm-noise]], [[billow-noise]] and [[ridgedmulti-noise]] but have no effect here, since only a single octave is evaluated.

  Returns a function of one, two or three double arguments (`x`, `[x y]` or `[x y z]`), returning a double from the `[0,1]` range when `:normalize?` is true (the default), or `[-1,1]` otherwise.

  See also [[fbm-noise]], [[billow-noise]], [[ridgedmulti-noise]], [[noise]], [[vnoise]], [[simplex]], [[random-noise-cfg]], [[random-noise]]."
  ([] (single-noise nil))
  ([cfg]
   (let [ncfg (noise-config cfg)]
     (fn
       (^double [^double x] (Noise/noise ncfg x))
       (^double [^double x ^double y] (Noise/noise ncfg x y))
       (^double [^double x ^double y ^double z] (Noise/noise ncfg x y z))))))

(defn fbm-noise
  "Creates a Fractal Brownian Motion (FBM) noise function.

  Sums several octaves of the underlying noise type, each with frequency scaled by `:lacunarity` and amplitude scaled by `:gain` relative to the previous octave, producing natural, self-similar terrain-like noise. [[noise]], [[vnoise]] and [[simplex]] are ready-to-use FBM noise functions with fixed presets.

  Parameters:

  - `cfg` (optional, map): noise configuration, all keys optional:
      - `:seed` - long seed for the noise's internal RNG, default: random.
      - `:noise-type` - `:value`, `:gradient` or `:simplex`, default: `:gradient`.
      - `:interpolation` - `:none`, `:linear`, `:hermite` or `:quintic`, used only by `:value` and `:gradient` noise types, default: `:hermite`.
      - `:octaves` - number of octaves summed together, default: `6`.
      - `:lacunarity` - frequency multiplier applied to each successive octave, default: `2.0`.
      - `:gain` - amplitude multiplier applied to each successive octave, default: `0.5`.
      - `:normalize?` - normalize result to `[0,1]` range (`true`, default) or leave it in `[-1,1]` (`false`).

  Returns a function of one, two or three double arguments (`x`, `[x y]` or `[x y z]`), returning a double from the `[0,1]` range when `:normalize?` is true (the default), or `[-1,1]` otherwise.

  See also [[single-noise]], [[billow-noise]], [[ridgedmulti-noise]], [[noise]], [[vnoise]], [[simplex]], [[random-noise-cfg]], [[random-noise]]."
  ([] (fbm-noise nil))
  ([cfg]
   (let [ncfg (noise-config cfg)]
     (fn
       (^double [^double x] (FBM/noise ncfg x))
       (^double [^double x ^double y] (FBM/noise ncfg x y))
       (^double [^double x ^double y ^double z] (FBM/noise ncfg x y z))))))

(defn billow-noise
  "Creates a billow noise function.

  Similar to [[fbm-noise]], sums several octaves of the underlying noise type scaled by `:lacunarity` and `:gain`, but folds each octave's value through `abs(v) * 2 - 1` before summing. This produces puffy, billowy cloud-like patterns instead of the smoother look of plain FBM noise.

  Parameters:

  - `cfg` (optional, map): noise configuration, all keys optional:
      - `:seed` - long seed for the noise's internal RNG, default: random.
      - `:noise-type` - `:value`, `:gradient` or `:simplex`, default: `:gradient`.
      - `:interpolation` - `:none`, `:linear`, `:hermite` or `:quintic`, used only by `:value` and `:gradient` noise types, default: `:hermite`.
      - `:octaves` - number of octaves summed together, default: `6`.
      - `:lacunarity` - frequency multiplier applied to each successive octave, default: `2.0`.
      - `:gain` - amplitude multiplier applied to each successive octave, default: `0.5`.
      - `:normalize?` - normalize result to `[0,1]` range (`true`, default) or leave it in `[-1,1]` (`false`).

  Returns a function of one, two or three double arguments (`x`, `[x y]` or `[x y z]`), returning a double from the `[0,1]` range when `:normalize?` is true (the default), or `[-1,1]` otherwise.

  See also [[single-noise]], [[fbm-noise]], [[ridgedmulti-noise]], [[random-noise-cfg]], [[random-noise]]."
  ([] (billow-noise nil))
  ([cfg]
   (let [ncfg (noise-config cfg)]
     (fn
       (^double [^double x] (Billow/noise ncfg x))
       (^double [^double x ^double y] (Billow/noise ncfg x y))
       (^double [^double x ^double y ^double z] (Billow/noise ncfg x y z))))))

(defn ridgedmulti-noise
  "Creates a ridged multifractal noise function.

  Combines several octaves of the underlying noise type into sharp, ridge-like features: each octave's value is folded through `(1 - abs(v))^2`, weighted by the strength of the previous octave's signal (scaled by `:gain` and clamped to `[0,1]`), with frequency scaled by `:lacunarity` as usual. This feedback between octaves produces the jagged mountain-range look typical of ridged multifractal noise.

  Parameters:

  - `cfg` (optional, map): noise configuration, all keys optional:
      - `:seed` - long seed for the noise's internal RNG, default: random.
      - `:noise-type` - `:value`, `:gradient` or `:simplex`, default: `:gradient`.
      - `:interpolation` - `:none`, `:linear`, `:hermite` or `:quintic`, used only by `:value` and `:gradient` noise types, default: `:hermite`.
      - `:octaves` - number of octaves combined together, default: `6`.
      - `:lacunarity` - frequency multiplier applied to each successive octave, default: `2.0`.
      - `:gain` - factor scaling the inter-octave weight feedback, default: `0.5`.
      - `:normalize?` - normalize result to `[0,1]` range (`true`, default) or leave it in `[-1,1]` (`false`).

  Returns a function of one, two or three double arguments (`x`, `[x y]` or `[x y z]`), returning a double from the `[0,1]` range when `:normalize?` is true (the default), or `[-1,1]` otherwise.

  See also [[single-noise]], [[fbm-noise]], [[billow-noise]], [[random-noise-cfg]], [[random-noise]]."
  ([] (ridgedmulti-noise nil))
  ([cfg]
   (let [ncfg (noise-config cfg)]
     (fn
       (^double [^double x] (RidgedMulti/noise ncfg x))
       (^double [^double x ^double y] (RidgedMulti/noise ncfg x y))
       (^double [^double x ^double y ^double z] (RidgedMulti/noise ncfg x y z))))))

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
  "Creates a domain-warped noise function.

  Applies domain warping as described in [Inigo Quilez's article](http://www.iquilezles.org/www/articles/warp/warp.htm): at each of `depth` warp levels, the input coordinates are recursively perturbed by evaluating the warp itself (with fixed constant offsets) at the previous level, scaled by `scale` and added to the coordinates, before evaluating `noise` at the warped position. A `depth` of `0` performs no warping and simply evaluates `noise` directly.

  Parameters:

  - `noise` (optional, function): the noise function to warp; must accept 1, 2 or 3 double arguments, same shape as [[vnoise]], [[noise]] or [[simplex]]. Defaults to [[vnoise]].
  - `scale` (optional, double): strength of the coordinate perturbation added at each warp level. Defaults to `4.0`.
  - `depth` (optional, long): number of recursive warp levels; each additional level multiplies the number of `noise` evaluations, so larger values are increasingly expensive. Defaults to `1`.

  Returns a function of one, two or three double arguments (`x`, `[x y]` or `[x y z]`), returning a double. The normalization (value range) of the returned function follows the normalization of `noise`.

  See also [[vnoise]], [[noise]], [[simplex]], [[random-noise]]."
  {:deprecated "Use `warp-noise` instead."}
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

(defn warp-noise
  "Creates a domain-warped noise function.

  Applies domain warping as described in [Inigo Quilez's article](http://www.iquilezles.org/www/articles/warp/warp.htm): at each of `depth` warp levels, the input coordinates are recursively perturbed by evaluating the warp itself (with fixed constant offsets) at the previous level, scaled by `scale` and added to the coordinates, before evaluating `noise` at the warped position. A `depth` of `0` performs no warping and simply evaluates `noise` directly.

  Parameters:

  - `noise` (optional, function): the noise function to warp; must accept 1, 2 or 3 double arguments, same shape as [[vnoise]], [[noise]] or [[simplex]]. Defaults to [[vnoise]].
  - `scale` (optional, double): strength of the coordinate perturbation added at each warp level. Defaults to `4.0`.
  - `depth` (optional, long): number of recursive warp levels; each additional level multiplies the number of `noise` evaluations, so larger values are increasingly expensive. Defaults to `1`.

  Returns a function of one, two or three double arguments (`x`, `[x y]` or `[x y z]`), returning a double. The normalization (value range) of the returned function follows the normalization of `noise`.

  See also [[vnoise]], [[noise]], [[simplex]], [[random-noise]]."
  ([noise ^double scale ^long depth]
   (let [n1 (make-warp-1d noise scale depth)
         n2 (make-warp-2d noise scale depth)
         n3 (make-warp-3d noise scale depth)]
     (fn
       (^double [^double x] (n1 x))
       (^double [^double x ^double y] (n2 x y))
       (^double [^double x ^double y ^double z] (n3 x y z)))))
  ([noise ^double scale] (warp-noise noise scale 1))
  ([noise] (warp-noise noise 4.0 1))
  ([] (warp-noise vnoise 4.0 1)))

(defonce ^{:doc "List of possible noise generators as a map of names and functions."}
  noise-generators
  {:fbm fbm-noise
   :single single-noise
   :billow billow-noise
   :ridgemulti ridgedmulti-noise})

(defn random-noise-cfg
  "Generates a randomized noise configuration map.

  Produces a configuration suitable for [[random-noise]] (or, for its `:seed`, `:noise-type`, `:interpolation`, `:octaves`, `:lacunarity`, `:gain` and `:normalize?` keys, for [[fbm-noise]], [[single-noise]], [[billow-noise]] and [[ridgedmulti-noise]] directly), with every key assigned a random but sensible value, so calling it repeatedly yields varied noise behavior without manual tuning.

  Parameters:

  - `pre-config` (optional, map): fixed values to use instead of randomizing; any keys present here override the corresponding randomly generated value. Defaults to an empty map, meaning every key is randomized.

  The randomized keys are:

  - `:seed` - a random integer seed.
  - `:generator` - one of `:single`, `:fbm`, `:billow`, `:ridgemulti`, used by [[random-noise]] to pick the blending method.
  - `:noise-type` - one of `:value`, `:gradient`, `:simplex`.
  - `:interpolation` - one of `:none`, `:linear`, `:hermite`, `:quintic`.
  - `:octaves` - integer between `1` and `9`.
  - `:lacunarity` - double between `1.5` and `2.5`.
  - `:gain` - double between `0.2` and `0.8`.
  - `:warp-scale` - `0.0` (no warp) with 80% probability, otherwise either `4.0` or a random double between `0.1` and `10.0`.
  - `:warp-depth` - `1` with 80% probability, otherwise a random integer between `1` and `3`.
  - `:normalize?` - always `true`.

  Returns a configuration map.

  See also [[random-noise]], [[fbm-noise]], [[single-noise]], [[billow-noise]], [[ridgedmulti-noise]], [[warp-noise]]."
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
  "Generates a fully random noise function.

  Combines [[random-noise-cfg]] with one of the noise blending methods (`:single`, `:fbm`, `:billow`, `:ridgemulti`, see `noise-generators`) and, when the resulting configuration requests warping (`:warp-scale` greater than `0.0`), wraps the noise with [[warp-noise]]. The result is a ready-to-use noise function in the same shape as [[noise]], [[vnoise]] and [[simplex]].

  Parameters:

  - `cfg` (optional, map): configuration overrides passed to [[random-noise-cfg]]; any key not provided is filled in randomly. Defaults to `nil`, meaning a fully random configuration.

  Returns a function of one, two or three double arguments (`x`, `[x y]` or `[x y z]`) returning a double.

  See also [[random-noise-cfg]], [[warp-noise]], [[fbm-noise]], [[single-noise]], [[billow-noise]], [[ridgedmulti-noise]]."
  {:deprecated "Use `random-noise` instead."}
  ([cfg]
   (let [cfg (random-noise-cfg cfg)
         gen-fn (noise-generators (get cfg :generator :fbm))
         noise (gen-fn cfg)]
     (if (pos? ^double (:warp-scale cfg))
       (warp-noise noise (:warp-scale cfg) (:warp-depth cfg))
       noise)))
  ([] (random-noise-fn nil)))

(defn random-noise
  "Generates a fully random noise function.

  Combines [[random-noise-cfg]] with one of the noise blending methods (`:single`, `:fbm`, `:billow`, `:ridgemulti`, see `noise-generators`) and, when the resulting configuration requests warping (`:warp-scale` greater than `0.0`), wraps the noise with [[warp-noise-fn]]. The result is a ready-to-use noise function in the same shape as [[noise]], [[vnoise]] and [[simplex]].

  Parameters:

  - `cfg` (optional, map): configuration overrides passed to [[random-noise-cfg]]; any key not provided is filled in randomly. Defaults to `nil`, meaning a fully random configuration.

  Returns a function of one, two or three double arguments (`x`, `[x y]` or `[x y z]`) returning a double.

  See also [[random-noise-cfg]], [[warp-noise]], [[fbm-noise]], [[single-noise]], [[billow-noise]], [[ridgedmulti-noise]]."
  ([cfg]
   (let [cfg (random-noise-cfg cfg)
         gen-fn (noise-generators (get cfg :generator :fbm))
         noise (gen-fn cfg)]
     (if (pos? ^double (:warp-scale cfg))
       (warp-noise noise (:warp-scale cfg) (:warp-depth cfg))
       noise)))
  ([] (random-noise nil)))


;; ### Discrete noise

(defn discrete-noise
  "Computes a deterministic hash-based noise value for one or two integer coordinates.

  Unlike stochastic random number generators, `discrete-noise` is a pure function: calling it repeatedly with the same `X` (and `Y`) always returns the same value, which makes it useful as a coordinate-based building block for procedural noise and hashing.

  Parameters:

  - `X` (long): first coordinate.
  - `Y` (long, optional): second coordinate. Defaults to `0`, producing a 1d hash of `X` alone.

  Returns a double in the `[0,1]` range.

  See also [[vnoise]], [[simplex]], [[random-noise]]."
  (^double [^long X ^long Y] (Discrete/value X Y))
  (^double [^long X] (Discrete/value X 0)))

;; Distributions

(defmulti
  ^{:doc "Create distribution object.

* First parameter is distribution as a `:key`.
* Second parameter is a map with configuration.

All distributions accept `rng` under `:rng` key (default: [[default-rng]]) and some of them accept `inverse-cumm-accuracy` (default set to `1e-9`).

Below is the full list of supported `:key`s, grouped by kind. For each: its accepted configuration-map parameters (`rng` included; every entry also has its own dedicated function of the same name, e.g. [[beta]], with a complete docstring covering formulas, defaults and cross-references) and a one-line description of what it models.

**Common continuous**

* `:beta` ([[beta]]) - `alpha`, `beta`, `inverse-abs-accuracy`, `rng` - proportions/probabilities on `[0,1]`, shaped by two parameters.
* `:cauchy` ([[cauchy]]) - `median`, `scale`, `inverse-abs-accuracy`, `rng` - symmetric, heavy-tailed, undefined mean/variance.
* `:chi-squared` ([[chi-squared]]) - `degrees-of-freedom`, `inverse-abs-accuracy`, `rng` - sum of squares of independent standard normals.
* `:exponential` ([[exponential]]) - `mean`, `inverse-abs-accuracy`, `rng` - memoryless waiting time between events.
* `:f` ([[f]]) - `numerator-degrees-of-freedom`, `denominator-degrees-of-freedom`, `inverse-abs-accuracy`, `rng` - ratio of two independent chi-squared variables.
* `:gamma` ([[gamma]]) - `shape`, `scale`, `inverse-abs-accuracy`, `rng` - waiting times, sums of exponentials.
* `:gumbel` ([[gumbel]]) - `mu`, `beta`, `rng` - type-I extreme value distribution for maxima/minima.
* `:laplace` ([[laplace]]) - `mu`, `beta`, `rng` - double exponential, sharper peak/heavier tails than normal.
* `:levy` ([[levy]]) - `mu`, `c`, `rng` - heavy-tailed, closed-form pdf/cdf but undefined moments.
* `:logistic` ([[logistic]]) - `mu`, `s`, `rng` - symmetric, sigmoid cdf, underlies logistic regression.
* `:log-normal` ([[log-normal]]) - `scale`, `shape`, `inverse-abs-accuracy`, `rng` - variable whose logarithm is normal.
* `:nakagami` ([[nakagami]]) - `mu`, `omega`, `inverse-abs-accuracy`, `rng` - amplitude of fading wireless signals.
* `:normal` ([[normal]]) - `mu`, `sd`, `inverse-abs-accuracy`, `rng` - the Gaussian bell curve.
* `:pareto` ([[pareto]]) - `scale`, `shape`, `inverse-abs-accuracy`, `rng` - heavy-tailed wealth/file-size type quantities.
* `:t` ([[t]]) - `degrees-of-freedom`, `inverse-abs-accuracy`, `rng` - Student's t, heavier tails than normal.
* `:triangular` ([[triangular]]) - `a`, `c`, `b`, `rng` - triangular density from lower bound, mode, upper bound.
* `:uniform-real` ([[uniform-real]]) - `lower`, `upper`, `rng` - continuous uniform over an interval.
* `:weibull` ([[weibull]]) - `alpha`, `beta`, `inverse-abs-accuracy`, `rng` - time-to-failure, reliability/survival analysis.
* `:constant` ([[constant]]) - `value` - degenerate (Dirac) distribution always returning the same value.

**Empirical / enumerated**

* `:empirical` ([[empirical]]) - `data`, `bin-count`, `rng` - histogram-estimated distribution from a data sample.
* `:enumerated-real` ([[enumerated-real]]) - `data`, `probabilities`, `rng` - explicit finite set of real values with probabilities.
* `:enumerated-int` ([[enumerated-int]]) - `data`, `probabilities`, `rng` - explicit finite set of integer values with probabilities.

**Common discrete**

* `:bernoulli` ([[bernoulli]]) - `p`, `rng` - single yes/no trial, `1` with probability `p`.
* `:binomial` ([[binomial]]) - `trials`, `p`, `rng` - number of successes across independent trials.
* `:geometric` ([[geometric]]) - `p`, `rng` - number of failures before the first success.
* `:hypergeometric` ([[hypergeometric]]) - `population-size`, `number-of-successes`, `sample-size`, `rng` - successes drawn without replacement.
* `:pascal` ([[pascal]]) - `r`, `p`, `rng` - failures observed before accumulating `r` successes.
* `:poisson` ([[poisson]]) - `p`, `epsilon`, `max-iterations`, `rng` - event count at a constant average rate (`p` is the rate/mean, traditionally lambda).
* `:uniform-int` ([[uniform-int]]) - `lower`, `upper`, `rng` - discrete uniform over an integer range.
* `:zipf` ([[zipf]]) - `number-of-elements`, `exponent`, `rng` - power-law rank/frequency distribution.

**Multivariate**

* `:multi-normal` ([[multi-normal]]) - `means`, `covariances`, `rng` - multivariate Gaussian over correlated vector components.

**Goodness-of-fit test statistics**

* `:anderson-darling` ([[anderson-darling]]) - `n`, `rng` - sampling distribution of the Anderson-Darling statistic.
* `:anderson-darling-quick` ([[anderson-darling-quick]]) - `n`, `rng` - same statistic, faster algorithm.
* `:cramer-von-mises` ([[cramer-von-mises]]) - `n`, `rng` - sampling distribution of the Cramer-von Mises statistic.
* `:kolmogorov-smirnov` ([[kolmogorov-smirnov]]) - `n`, `rng` - sampling distribution of the two-sided KS statistic.
* `:kolmogorov-smirnov+` ([[kolmogorov-smirnov+]]) - `n`, `rng` - one-sided (D+) KS statistic.
* `:kolmogorov-smirnov-quick` ([[kolmogorov-smirnov-quick]]) - `n`, `rng` - two-sided KS statistic, faster algorithm.
* `:kolmogorov` ([[kolmogorov]]) - `rng` - parameter-free limiting distribution of the scaled KS statistic.
* `:watson-g` ([[watson-g]]) - `n`, `rng` - Watson G statistic for circular (directional) goodness-of-fit.
* `:watson-u` ([[watson-u]]) - `n`, `rng` - Watson U-squared statistic for circular goodness-of-fit.

**Less-common continuous**

* `:beta-symmetrical` ([[beta-symmetrical]]) - `alpha`, `d`, `rng` - symmetric special case of beta (equal shape parameters).
* `:chi` ([[chi]]) - `nu`, `rng` - square root of a sum of squares of `nu` standard normals.
* `:erlang` ([[erlang]]) - `k`, `lambda`, `rng` - sum of `k` independent, identically-rated exponential stages.
* `:fatigue-life` ([[fatigue-life]]) - `alpha`, `beta`, `gamma`, `rng` - Birnbaum-Saunders model, time to failure under cyclic stress.
* `:folded-normal` ([[folded-normal]]) - `mu`, `sigma`, `rng` - absolute value of a normal random variable.
* `:frechet` ([[frechet]]) - `alpha`, `beta`, `delta`, `rng` - heavy-tailed type-II extreme value distribution.
* `:half-normal` ([[half-normal]]) - `mu`, `sigma`, `rng` - one-sided fold of a normal distribution.
* `:hyperbolic-secant` ([[hyperbolic-secant]]) - `mu`, `sigma`, `rng` - symmetric shape between normal and Cauchy.
* `:hypoexponential-equal` ([[hypoexponential-equal]]) - `n`, `k`, `h`, `rng` - sum of `k` phases with equally-spaced exponential rates.
* `:hypoexponential` ([[hypoexponential]]) - `lambdas`, `rng` - sum of independent exponentials with arbitrary, distinct rates.
* `:inverse-gamma` ([[inverse-gamma]]) - `alpha`, `beta`, `rng` - reciprocal of a gamma-distributed variable.
* `:inverse-gaussian` ([[inverse-gaussian]]) - `mu`, `lambda`, `rng` - Wald distribution, first-passage time of a drifting Brownian motion.
* `:johnson-sb` ([[johnson-sb]]) - `gamma`, `delta`, `xi`, `lambda`, `rng` - bounded member of the Johnson system.
* `:johnson-sl` ([[johnson-sl]]) - `gamma`, `delta`, `xi`, `lambda`, `rng` - semi-bounded, log-normal-like member of the Johnson system.
* `:johnson-su` ([[johnson-su]]) - `gamma`, `delta`, `xi`, `lambda`, `rng` - unbounded member of the Johnson system.
* `:log-logistic` ([[log-logistic]]) - `alpha`, `beta`, `rng` - Fisk distribution, logarithm of the variable is logistic.
* `:normal-inverse-gaussian` ([[normal-inverse-gaussian]]) - `alpha`, `beta`, `mu`, `delta`, `rng` - heavy-tailed normal-variance mixture, asset returns.
* `:pearson-6` ([[pearson-6]]) - `alpha1`, `alpha2`, `beta`, `rng` - scaled beta distribution of the second kind.
* `:power` ([[power]]) - `a`, `b`, `c`, `rng` - power-function distribution generalizing the uniform on `[a,b]`.
* `:rayleigh` ([[rayleigh]]) - `a`, `beta`, `rng` - magnitude of a 2D vector with independent zero-mean normal components.
* `:half-cauchy` ([[half-cauchy]]) - `mu`, `scale`, `rng` - right half of a Cauchy distribution folded at its center.
* `:reciprocal` ([[reciprocal]]) - `a`, `b`, `rng` - log-uniform distribution for scale-invariant quantities.
* `:ex-gaussian` ([[ex-gaussian]]) - `mu`, `sigma`, `tau`, `rng` - sum of a normal and an exponential variable, models reaction times.
* `:exgaus` ([[exgaus]]) - `mu`, `sigma`, `nu`, `rng` - `ex-gaussian` with gamlss-style parameter naming (`nu` for `tau`).
* `:von-mises` ([[von-mises]]) - `mu`, `kappa`, `rng` - circular analogue of the normal distribution, for angular data.

**Noncentral family**

* `:chi-squared-noncentral` ([[chi-squared-noncentral]]) - `nu`, `lambda`, `rng` - chi-squared generalized to noncentered normal components.
* `:f-noncentral` ([[f-noncentral]]) - `df1`, `df2`, `ncp`, `rng` - F-distribution generalized to a noncentral numerator.
* `:t-noncentral` ([[t-noncentral]]) - `df`, `ncp`, `rng` - Student's t generalized to a noncentral numerator.
* `:beta-noncentral` ([[beta-noncentral]]) - `alpha`, `beta`, `ncp`, `rng` - beta distribution generalized via a noncentral chi-squared numerator.
* `:fishers-noncentral-hypergeometric` ([[fishers-noncentral-hypergeometric]]) - `ns`, `nf`, `n`, `omega`, `rng` - biased hypergeometric, conditioned on a 2x2-table odds ratio.
* `:wallenius-noncentral-hypergeometric` ([[wallenius-noncentral-hypergeometric]]) - `ns`, `nf`, `n`, `omega`, `rng` - biased sequential (urn) hypergeometric sampling.

**Multivariate discrete**

* `:multinomial` ([[multinomial]]) - `n`, `ps`, `rng` - category counts from `n` trials with per-category probabilities.
* `:dirichlet` ([[dirichlet]]) - `alpha`, `rng` - distribution over the probability simplex, conjugate prior for multinomial.
* `:categorical-distribution`/`:categorical` ([[categorical-distribution]]/[[categorical]]) - `data`, `probabilities`, `rng` - discrete distribution over an arbitrary, non-numeric set of values.

**Data-driven**

* `:continuous-distribution`/`:kde` ([[continuous-distribution]]/[[kde]]) - `data`, `kde`, `bandwidth`, `steps`, `interpolator`, `rng` - nonparametric continuous distribution via kernel density estimation.

**Other discrete**

* `:negative-binomial` ([[negative-binomial]]) - `r`, `p`, `rng` - generalized (Polya) negative binomial, failures before `r` (possibly non-integer) successes.
* `:nbi` ([[nbi]]) - `mu`, `sigma`, `rng` - negative binomial (type I), gamlss mean/dispersion reparametrization.
* `:nbii` ([[nbii]]) - `mu`, `sigma`, `rng` - negative binomial (type II), gamlss mean/dispersion reparametrization.
* `:logarithmic` ([[logarithmic]]) - `p`, `rng` - log-series distribution over positive integers, species-abundance data.
* `:integer-discrete-distribution`/`:integer-discrete` ([[integer-discrete-distribution]]/[[integer-discrete]]) - `data`, `probabilities`, `rng` - arbitrary discrete distribution over a finite integer support.
* `:real-discrete-distribution`/`:real-discrete` ([[real-discrete-distribution]]/[[real-discrete]]) - `data`, `probabilities`, `rng` - arbitrary discrete distribution over a finite real-valued support.
* `:beta-binomial` ([[beta-binomial]]) - `alpha`, `beta`, `n`, `rng` - binomial with a beta-distributed success probability, overdispersed counts.
* `:bb` ([[bb]]) - `mu`, `sigma`, `bd`, `rng` - beta-binomial, gamlss mean/dispersion reparametrization.

**Zero-inflated / zero-adjusted (gamlss) family**

* `:zero-inflated-binomial`/`:zibi` ([[zero-inflated-binomial]]) - `mu`, `sigma`, `bd`, `rng` - binomial with an extra point mass at zero.
* `:zero-adjusted-binomial`/`:zabi` ([[zero-adjusted-binomial]]) - `mu`, `sigma`, `bd`, `rng` - hurdle model: exact zero w.p. `sigma`, else zero-truncated binomial.
* `:zero-inflated-beta-binomial`/`:zibb` ([[zero-inflated-beta-binomial]]) - `mu`, `sigma`, `bd`, `nu`, `rng` - beta-binomial with an extra point mass at zero.
* `:zero-adjusted-beta-binomial`/`:zabb` ([[zero-adjusted-beta-binomial]]) - `mu`, `sigma`, `bd`, `nu`, `rng` - hurdle model over beta-binomial.
* `:zero-inflated-negative-binomial`/`:zinbi` ([[zero-inflated-negative-binomial]]) - `mu`, `sigma`, `nu`, `rng` - negative binomial (`nbi`) with an extra point mass at zero.
* `:zero-adjusted-negative-binomial`/`:zanbi` ([[zero-adjusted-negative-binomial]]) - `mu`, `sigma`, `nu`, `rng` - hurdle model over `nbi`.
* `:zero-inflated-poisson`/`:zip` ([[zero-inflated-poisson]]) - `mu`, `sigma`, `rng` - Poisson with an extra point mass at zero.
* `:zero-inflated-poisson2`/`:zip2` ([[zero-inflated-poisson2]]) - `mu`, `sigma`, `rng` - mean-parameterized reparametrization of zero-inflated Poisson.
* `:zero-adjusted-poisson`/`:zap` ([[zero-adjusted-poisson]]) - `mu`, `sigma`, `rng` - hurdle model over Poisson.
* `:zero-adjusted-gamma`/`:zaga` ([[zero-adjusted-gamma]]) - `mu`, `sigma`, `nu`, `rng` - hurdle model mixing a point mass at zero with a gamma distribution.
* `:zero-adjusted-inverse-gaussian`/`:zaig` ([[zero-adjusted-inverse-gaussian]]) - `mu`, `sigma`, `nu`, `rng` - hurdle model mixing a point mass at zero with an inverse Gaussian.

**Generalized family**

* `:generalized-extreme-value`/`:gev` ([[generalized-extreme-value]]) - `mu`, `sigma`, `xi`, `rng` - unifies Gumbel/Frechet/Weibull as limits of normalized maxima.
* `:generalized-logistic` ([[generalized-logistic]]) - `mu`, `sigma`, `alpha`, `rng` - skewed generalization of the logistic distribution.
* `:generalized-pareto`/`:gpd` ([[generalized-pareto]]) - `mu`, `sigma`, `xi`, `rng` - limiting distribution of excesses over a threshold.
* `:generalized-exponential`/`:ge` ([[generalized-exponential]]) - `alpha`, `lambda`, `rng` - exponentiated exponential distribution.
* `:generalized-gamma`/`:gg` ([[generalized-gamma]]) - `mu`, `sigma`, `nu`, `rng` - Stacy distribution generalizing gamma, Weibull, log-normal.
* `:generalized-normal`/`:gnd` ([[generalized-normal]]) - `mu`, `alpha`, `beta`, `rng` - exponential power (Subbotin) distribution.
* `:generalized-inverse-gaussian`/`:gig` ([[generalized-inverse-gaussian]]) - `chi`, `psi`, `lambda`, `rng` - generalizes gamma, inverse-gamma, inverse Gaussian.
* `:generalized-hyperbolic`/`:gh` ([[generalized-hyperbolic]]) - `mu`, `delta`, `alpha`, `beta`, `lambda`, `rng` - flexible skew/kurtosis distribution, financial returns.
* `:half-logistic` ([[half-logistic]]) - `scale`, `rng` - absolute value of a logistic random variable.
* `:generalized-half-logistic`/`:ghl` ([[generalized-half-logistic]]) - `alpha`, `lambda`, `rng` - exponentiated half-logistic distribution.

**Meta / combinator distributions**

* `:truncated` ([[truncated]]) - `distr`, `left`, `right`, `rng` - restricts an existing distribution's support to `[left, right]`.
* `:mixture` ([[mixture]]) - `distrs`, `weights`, `rng` - finite weighted mixture combining several component distributions."}
  distribution (fn ([k _] k) ([k] k)))

(defmacro ^:private add-distr-method
  ([d]
   `(add-distr-method ~d ~(keyword d)))
  ([d kd]
   `(defmethod distribution ~kd
      ([_#] (distribution ~kd nil))
      ([_# opts#] (~d opts#)))))

;;

(defn beta
  "Creates a beta distribution object.

  The beta distribution is a continuous distribution supported on the interval `[0, 1]`, shaped by two positive parameters `alpha` and `beta` that control its skewness. It is commonly used to model random proportions and probabilities, and as a conjugate prior for the Bernoulli/binomial parameter in Bayesian statistics.

  Parameters (single, optional map):

  - `alpha` (double): first shape parameter. Default: `2.0`.
  - `beta` (double): second shape parameter. Default: `2.0`.
  - `inverse-abs-accuracy` (double): accuracy used when inverting the cumulative distribution function (`icdf`). Default: the implementation's default accuracy.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `BetaDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[gamma]], [[uniform-real]]."
  (^BetaDistribution [] (beta nil))
  (^BetaDistribution [{:keys [^double alpha ^double beta ^double inverse-abs-accuracy rng]
                       :or {alpha 2.0 beta 2.0 inverse-abs-accuracy BetaDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY}}]
   (BetaDistribution. (or rng (JDKRandomGenerator.)) alpha beta inverse-abs-accuracy)))

(add-distr-method beta)

(defn beta-noncentral
  "Creates a noncentral beta distribution object.

  The noncentral beta distribution is a continuous distribution over `[0, 1]`, generalizing the [[beta]] distribution to the case where the underlying chi-squared random variable in the numerator has a nonzero, shared mean (encoded through the noncentrality parameter `ncp`): if `X1` follows a [[chi-squared-noncentral]] distribution with `2*alpha` degrees of freedom and noncentrality `ncp`, and `X2` follows an independent central [[chi-squared]] distribution with `2*beta` degrees of freedom, then `X1/(X1+X2)` follows this distribution. It is used, among others, in power calculations for tests on proportions and correlation coefficients under a non-null alternative hypothesis.

  Internally it is computed as a Poisson(`ncp/2`)-weighted mixture of central [[beta]] distributions with first shape parameter `alpha + k`, truncated once the cumulative Poisson mass is within `1e-15` of 1 (the same mixture representation used by [[f-noncentral]] and [[chi-squared-noncentral]]).

  Parameters (single, optional map):

  - `alpha` (double): first shape parameter, strictly positive. Default: `2.0`.
  - `beta` (double): second shape parameter, strictly positive. Default: `2.0`.
  - `ncp` (double): noncentrality parameter, non-negative; `ncp = 0.0` reduces the distribution exactly to a central [[beta]] distribution with the same `alpha`/`beta`. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Note: `pdf(0)` follows a boundary convention based on `alpha` (`0` for `alpha > 1`, the finite value of the `k=0` mixture term for `alpha = 1`, `##Inf` for `alpha < 1`), and `pdf(1)` follows the mirrored convention based on `beta`. These are computed directly from the known analytic limits rather than the underlying Apache Commons `BetaDistribution` implementation, which throws an exception at `x=0`/`x=1` when the corresponding shape parameter is below `1`, and silently returns the wrong value (`0.0` instead of the true finite limit) when it is exactly `1`.

  `mean`/`variance` are exact weighted sums of the mixture's own central-beta component means/variances (not an approximation, and not restricted to any parameter range, unlike [[f-noncentral]]'s). There is no closed-form `cdf`/`icdf`; `cdf` sums the Poisson-weighted mixture terms directly (each a call into the well-tested Apache Commons `BetaDistribution` implementation away from the `0`/`1` boundaries), and `icdf` root-finds on that `cdf` over its known `[0, 1]` domain.

  Matches the `(shape1, shape2, ncp)` parameterization used by R's base `stats` package (`beta`/`pbeta`/`qbeta`).

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[beta]], [[f-noncentral]], [[chi-squared-noncentral]]."
  ([] (beta-noncentral nil))
  ([{:keys [^double alpha ^double beta ^double ncp rng]
     :or {alpha 2.0 beta 2.0 ncp 1.0}}]
   (distr/beta-noncentral alpha beta ncp rng)))

(add-distr-method beta-noncentral)

(defn cauchy
  "Creates a Cauchy distribution object.

  The Cauchy distribution is a continuous, symmetric, heavy-tailed distribution centered on `median` with spread controlled by `scale`. Its mean and variance are undefined due to the heaviness of its tails, which makes it a useful example and stress-test case for statistical methods that assume finite moments.

  Parameters (single, optional map):

  - `median` (double): location parameter, the center and median of the distribution. Default: `0.0`.
  - `scale` (double): scale parameter controlling the spread of the distribution. Default: `1.0`.
  - `inverse-abs-accuracy` (double): accuracy used when inverting the cumulative distribution function (`icdf`). Default: the implementation's default accuracy.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `CauchyDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[levy]], [[t]]."
  (^CauchyDistribution [] (cauchy nil))
  (^CauchyDistribution [{:keys [^double median ^double scale ^double inverse-abs-accuracy rng]
                         :or {median 0.0 scale 1.0 inverse-abs-accuracy CauchyDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY}}]
   (CauchyDistribution. (or rng (JDKRandomGenerator.)) median scale inverse-abs-accuracy)))

(add-distr-method cauchy)

(defn chi-squared
  "Creates a chi-squared distribution object.

  The chi-squared distribution is a continuous distribution over non-negative reals, arising as the distribution of a sum of squares of `degrees-of-freedom` independent standard normal random variables. It is widely used in hypothesis testing (chi-squared tests, goodness-of-fit) and interval estimation.

  Parameters (single, optional map):

  - `degrees-of-freedom` (double): number of degrees of freedom, shaping the distribution. Default: `1.0`.
  - `inverse-abs-accuracy` (double): accuracy used when inverting the cumulative distribution function (`icdf`). Default: the implementation's default accuracy.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `ChiSquaredDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[f]], [[t]], [[gamma]]."
  (^ChiSquaredDistribution [] (chi-squared nil))
  (^ChiSquaredDistribution [{:keys [^double degrees-of-freedom ^double inverse-abs-accuracy rng]
                             :or {degrees-of-freedom 1.0 inverse-abs-accuracy ChiSquaredDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY}}]
   (ChiSquaredDistribution. (or rng (JDKRandomGenerator.)) degrees-of-freedom inverse-abs-accuracy)))

(add-distr-method chi-squared)

(defn constant
  "Creates a constant, Dirac (degenerate) distribution object.

  The constant distribution always returns the same `value` with probability `1.0`; it has zero variance and represents a deterministic random variable. It is mostly useful as a placeholder or edge case where a distribution object is expected but the value should not actually vary.

  Parameters (single, optional map):

  - `value` (double): the single value always returned by the distribution. Default: `0.0`.

  Called with no arguments or with `nil`, creates the distribution with default parameter values. Unlike other continuous distributions here, this constructor does not accept an `rng`, since no randomness is involved.

  Returns a `ConstantRealDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[uniform-real]]."
  (^ConstantRealDistribution [] (constant nil))
  (^ConstantRealDistribution [{:keys [^double value]
                               :or {value 0.0}}]
   (ConstantRealDistribution. value)))

(add-distr-method constant)

(defn exponential
  "Creates an exponential distribution object.

  The exponential distribution is a continuous distribution over non-negative reals, modelling the waiting time between independent events occurring at a constant average rate. It is the continuous analogue of the [[geometric]] distribution and is memoryless.

  Parameters (single, optional map):

  - `mean` (double): mean of the distribution, ie. the average waiting time (the reciprocal of the rate). Default: `1.0`.
  - `inverse-abs-accuracy` (double): accuracy used when inverting the cumulative distribution function (`icdf`). Default: the implementation's default accuracy.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns an `ExponentialDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[geometric]], [[gamma]], [[weibull]]."
  (^ExponentialDistribution [] (exponential nil))
  (^ExponentialDistribution [{:keys [^double mean ^double inverse-abs-accuracy rng]
                              :or {mean 1.0 inverse-abs-accuracy ExponentialDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY}}]
   (ExponentialDistribution. (or rng (JDKRandomGenerator.)) mean inverse-abs-accuracy)))

(add-distr-method exponential)

(defn f
  "Creates an F-distribution (Fisher-Snedecor) object.

  The F-distribution is a continuous distribution over non-negative reals, arising as the ratio of two independent chi-squared random variables, each divided by their own degrees of freedom. It is widely used for comparing variances and in analysis of variance (ANOVA) and regression F-tests.

  Parameters (single, optional map):

  - `numerator-degrees-of-freedom` (double): degrees of freedom of the numerator chi-squared variable. Default: `1.0`.
  - `denominator-degrees-of-freedom` (double): degrees of freedom of the denominator chi-squared variable. Default: `1.0`.
  - `inverse-abs-accuracy` (double): accuracy used when inverting the cumulative distribution function (`icdf`). Default: the implementation's default accuracy.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns an `FDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[chi-squared]], [[t]]."
  (^FDistribution [] (f nil))
  (^FDistribution [{:keys [^double numerator-degrees-of-freedom ^double denominator-degrees-of-freedom ^double inverse-abs-accuracy rng]
                    :or {numerator-degrees-of-freedom 1.0 denominator-degrees-of-freedom 1.0
                         inverse-abs-accuracy FDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY}}]
   (FDistribution. (or rng (JDKRandomGenerator.)) numerator-degrees-of-freedom denominator-degrees-of-freedom inverse-abs-accuracy)))

(add-distr-method f)

(defn gamma
  "Creates a gamma distribution object.

  The gamma distribution is a continuous distribution over positive reals, shaped by a `shape` parameter and a `scale` parameter. It generalizes the [[exponential]] and [[chi-squared]] distributions and is commonly used to model waiting times, sums of exponential variables, and as a conjugate prior in Bayesian statistics.

  Parameters (single, optional map):

  - `shape` (double): shape parameter, controlling the form of the distribution. Default: `2.0`.
  - `scale` (double): scale parameter, stretching or shrinking the distribution. Default: `2.0`.
  - `inverse-abs-accuracy` (double): accuracy used when inverting the cumulative distribution function (`icdf`). Default: the implementation's default accuracy.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `GammaDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[exponential]], [[chi-squared]], [[beta]]."
  (^GammaDistribution [] (gamma nil))
  (^GammaDistribution [{:keys [^double shape ^double scale ^double inverse-abs-accuracy rng]
                        :or {shape 2.0 scale 2.0 inverse-abs-accuracy GammaDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY}}]
   (GammaDistribution. (or rng (JDKRandomGenerator.)) shape scale inverse-abs-accuracy)))

(add-distr-method gamma)

(defn gumbel
  "Creates a Gumbel distribution object.

  The Gumbel distribution is a continuous, right-skewed distribution used to model the maximum (or minimum) of a number of samples of other distributions; it is a type-I extreme value distribution. It is commonly used in flood, wind and other extreme-event analysis.

  Parameters (single, optional map):

  - `mu` (double): location parameter, shifting the mode of the distribution. Default: `1.0`.
  - `beta` (double): scale parameter, controlling the spread of the distribution. Default: `2.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `GumbelDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[levy]], [[laplace]]."
  (^GumbelDistribution [] (gumbel nil))
  (^GumbelDistribution [{:keys [^double mu ^double beta rng]
                         :or {mu 1.0 beta 2.0}}]
   (GumbelDistribution. (or rng (JDKRandomGenerator.)) mu beta)))

(add-distr-method gumbel)

(defn laplace
  "Creates a Laplace distribution object.

  The Laplace distribution, also known as the double exponential distribution, is a continuous, symmetric distribution formed from two exponential distributions back to back around a location parameter. Compared to the [[normal]] distribution it has a sharper peak and heavier tails, making it useful for modelling data with more extreme outliers.

  Parameters (single, optional map):

  - `mu` (double): location parameter, the center and mean of the distribution. Default: `0.0`.
  - `beta` (double): scale parameter, controlling the spread of the distribution. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `LaplaceDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[normal]], [[logistic]]."
  (^LaplaceDistribution [] (laplace nil))
  (^LaplaceDistribution [{:keys [^double mu ^double beta rng]
                          :or {mu 0.0 beta 1.0}}]
   (LaplaceDistribution. (or rng (JDKRandomGenerator.)) mu beta)))

(add-distr-method laplace)

(defn levy
  "Creates a Levy distribution object.

  The Levy distribution is a continuous, heavy-tailed distribution supported on `[mu, Infinity)`, notable for having a closed-form probability density and cumulative distribution function despite an undefined mean and variance. It is a special case of the inverse-gamma distribution family and appears in random walk and stable-distribution theory.

  Parameters (single, optional map):

  - `mu` (double): location parameter, the lower bound of the support. Default: `0.0`.
  - `c` (double): scale parameter, controlling the spread of the distribution. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `LevyDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[cauchy]], [[gumbel]]."
  (^LevyDistribution [] (levy nil))
  (^LevyDistribution [{:keys [^double mu ^double c rng]
                       :or {mu 0.0 c 1.0}}]
   (LevyDistribution. (or rng (JDKRandomGenerator.)) mu c)))

(add-distr-method levy)

(defn logistic
  "Creates a logistic distribution object.

  The logistic distribution is a continuous, symmetric distribution whose cumulative distribution function is the logistic (sigmoid) function. It resembles the [[normal]] distribution in shape but has heavier tails, and underlies logistic regression as well as some growth models.

  Parameters (single, optional map):

  - `mu` (double): location parameter, the center and mean of the distribution. Default: `0.0`.
  - `s` (double): scale parameter, controlling the spread of the distribution. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `LogisticDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[normal]], [[laplace]]."
  (^LogisticDistribution [] (logistic nil))
  (^LogisticDistribution [{:keys [^double mu ^double s rng]
                           :or {mu 0.0 s 1.0}}]
   (LogisticDistribution. (or rng (JDKRandomGenerator.)) mu s)))

(add-distr-method logistic)

(defn log-normal
  "Creates a log-normal distribution object.

  The log-normal distribution is a continuous distribution over positive reals whose logarithm follows a [[normal]] distribution. It is commonly used to model quantities that result from the multiplicative combination of many independent positive factors, such as incomes, stock prices or particle sizes.

  Parameters (single, optional map):

  - `scale` (double): location parameter of the underlying normal distribution of the logarithm of the variable (often denoted `mu`). Default: `1.0`.
  - `shape` (double): scale parameter of the underlying normal distribution of the logarithm of the variable (often denoted `sigma`); controls the spread. Default: `1.0`.
  - `inverse-abs-accuracy` (double): accuracy used when inverting the cumulative distribution function (`icdf`). Default: the implementation's default accuracy.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `LogNormalDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[normal]], [[pareto]]."
  (^LogNormalDistribution [] (log-normal nil))
  (^LogNormalDistribution [{:keys [^double scale ^double shape ^double inverse-abs-accuracy rng]
                            :or {scale 1.0 shape 1.0 inverse-abs-accuracy LogNormalDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY}}]
   (LogNormalDistribution. (or rng (JDKRandomGenerator.)) scale shape inverse-abs-accuracy)))

(add-distr-method log-normal)

(defn nakagami
  "Creates a Nakagami distribution object.

  The Nakagami distribution is a continuous distribution over positive reals, often used to model the amplitude of fading wireless communication signals. Its shape parameter `mu` controls the fading severity, while `omega` controls the average signal power (spread).

  Parameters (single, optional map):

  - `mu` (double): shape parameter controlling the severity of fading; must be at least `0.5`. Default: `1.0`.
  - `omega` (double): spread parameter, the average of the squared random variable. Default: `1.0`.
  - `inverse-abs-accuracy` (double): accuracy used when inverting the cumulative distribution function (`icdf`). Default: the implementation's default accuracy.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `NakagamiDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[weibull]], [[gamma]]."
  (^NakagamiDistribution [] (nakagami nil))
  (^NakagamiDistribution [{:keys [^double mu ^double omega ^double inverse-abs-accuracy rng]
                           :or {mu 1.0 omega 1.0 inverse-abs-accuracy NakagamiDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY}}]
   (NakagamiDistribution. (or rng (JDKRandomGenerator.)) mu omega inverse-abs-accuracy)))

(add-distr-method nakagami)

(defn normal
  "Creates a normal (Gaussian) distribution object.

  The normal distribution is a continuous, symmetric, bell-shaped distribution fully described by its mean `mu` and standard deviation `sd`. It is the most widely used distribution in statistics, arising naturally as the limit of sums of many independent random variables via the central limit theorem.

  Parameters (single, optional map):

  - `mu` (double): mean of the distribution, its center of symmetry. Default: `0.0`.
  - `sd` (double): standard deviation, controlling the spread of the distribution. Default: `1.0`.
  - `inverse-abs-accuracy` (double): accuracy used when inverting the cumulative distribution function (`icdf`). Default: the implementation's default accuracy.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `NormalDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[log-normal]], [[laplace]], [[logistic]]."
  (^NormalDistribution [] (normal nil))
  (^NormalDistribution [{:keys [^double mu ^double sd ^double inverse-abs-accuracy rng]
                         :or {mu 0.0 sd 1.0 inverse-abs-accuracy NormalDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY}}]
   (NormalDistribution. (or rng (JDKRandomGenerator.)) mu sd inverse-abs-accuracy)))

(add-distr-method normal)

(defn pareto
  "Creates a Pareto (Type I) distribution object.

  The Pareto distribution is a continuous, heavy-tailed distribution supported on `[scale, Infinity)`, modelling quantities where a small share of values accounts for a large share of the total, such as wealth distributions, file sizes or the classic 80/20 rule. Its `shape` parameter controls how quickly probability decays for larger values.

  Parameters (single, optional map):

  - `scale` (double): scale parameter, the minimum possible value of the distribution. Default: `1.0`.
  - `shape` (double): shape parameter (also known as the tail index or `alpha`); smaller values give heavier tails. Default: `1.0`.
  - `inverse-abs-accuracy` (double): accuracy used when inverting the cumulative distribution function (`icdf`). Default: the implementation's default accuracy.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `ParetoDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[zipf]], [[log-normal]]."
  (^ParetoDistribution [] (pareto nil))
  (^ParetoDistribution [{:keys [^double scale ^double shape ^double inverse-abs-accuracy rng]
                         :or {scale 1.0 shape 1.0 inverse-abs-accuracy ParetoDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY}}]
   (ParetoDistribution. (or rng (JDKRandomGenerator.)) scale shape inverse-abs-accuracy)))

(add-distr-method pareto)

(defn t
  "Creates a Student's t-distribution object.

  Student's t-distribution is a continuous, symmetric, bell-shaped distribution similar to the [[normal]] distribution but with heavier tails, controlled by its `degrees-of-freedom`. As degrees of freedom grow, the distribution approaches the standard normal distribution; it is widely used for inference about a mean when the sample size is small or the variance is unknown.

  Parameters (single, optional map):

  - `degrees-of-freedom` (double): number of degrees of freedom, shaping the heaviness of the tails. Default: `1.0`.
  - `inverse-abs-accuracy` (double): accuracy used when inverting the cumulative distribution function (`icdf`). Default: the implementation's default accuracy.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `TDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[normal]], [[cauchy]], [[f]]."
  (^TDistribution [] (t nil))
  (^TDistribution [{:keys [^double degrees-of-freedom ^double inverse-abs-accuracy rng]
                    :or {degrees-of-freedom 1.0 inverse-abs-accuracy TDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY}}]
   (TDistribution. (or rng (JDKRandomGenerator.)) degrees-of-freedom inverse-abs-accuracy)))

(add-distr-method t)

(defn triangular
  "Creates a triangular distribution object.

  The triangular distribution is a continuous distribution supported on `[a, b]`, with probability density rising linearly from `a` to the mode `c` and then falling linearly to `b`, forming a triangle shape. It is often used as a simple model when only a minimum, maximum and most likely value are known, such as in project estimation.

  Parameters (single, optional map):

  - `a` (double): lower limit of the support. Default: `-1.0`.
  - `c` (double): mode, the most likely value; must lie between `a` and `b`. Default: `0.0`.
  - `b` (double): upper limit of the support. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `TriangularDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[uniform-real]]."
  (^TriangularDistribution [] (triangular nil))
  (^TriangularDistribution [{:keys [^double a ^double c ^double b ^double rng]
                             :or {a -1.0 c 0.0 b 1.0}}]
   (TriangularDistribution. (or rng (JDKRandomGenerator.)) a c b)))

(add-distr-method triangular)

(defn uniform-real
  "Creates a continuous uniform distribution object over the interval `[lower, upper]`.

  Every value within `[lower, upper]` is equally likely; the density is constant across the support and zero outside of it. It is the continuous analogue of [[uniform-int]].

  Parameters (single, optional map):

  - `lower` (double): lower bound of the support, inclusive. Default: `0.0`.
  - `upper` (double): upper bound of the support, inclusive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `UniformRealDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[uniform-int]], [[triangular]]."
  (^UniformRealDistribution [] (uniform-real nil))
  (^UniformRealDistribution [{:keys [^double lower ^double upper rng]
                              :or {lower 0.0 upper 1.0}}]
   (let [^RandomGenerator rng (or rng (JDKRandomGenerator.))]
     (UniformRealDistribution. rng lower upper))))

(add-distr-method uniform-real)

(defn weibull
  "Creates a Weibull distribution object.

  The Weibull distribution is a continuous distribution over non-negative reals, widely used in reliability engineering and survival analysis to model time-to-failure data. Its shape parameter `alpha` determines whether the failure rate increases, decreases or stays constant over time (with `alpha` equal to `1.0` reducing to the [[exponential]] distribution).

  Parameters (single, optional map):

  - `alpha` (double): shape parameter, controlling the failure-rate behaviour over time. Default: `1.0`.
  - `beta` (double): scale parameter, stretching or shrinking the distribution. Default: `1.0`.
  - `inverse-abs-accuracy` (double): accuracy used when inverting the cumulative distribution function (`icdf`). Default: the implementation's default accuracy.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `WeibullDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[exponential]], [[nakagami]], [[gumbel]]."
  (^WeibullDistribution [] (weibull nil))
  (^WeibullDistribution [{:keys [^double alpha ^double beta ^double inverse-abs-accuracy rng]
                          :or {alpha 1.0 beta 1.0 inverse-abs-accuracy WeibullDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY}}]
   (WeibullDistribution. (or rng (JDKRandomGenerator.)) alpha beta inverse-abs-accuracy)))

(add-distr-method weibull)

(defn empirical
  "Creates an empirical distribution object estimated from a sample of `data`.

  The empirical distribution builds a histogram out of the provided `data` and treats it as a piecewise, continuous distribution: sampling picks a bin according to its observed frequency and then draws uniformly within that bin. It is useful for approximating an unknown continuous distribution directly from observations, without assuming any particular parametric form.

  Parameters (single, optional map):

  - `data` (sequence of doubles): sample used to build the histogram. Default: `[1.0]`.
  - `bin-count` (long): number of histogram bins. Default: estimated automatically from `data` using the Freedman-Diaconis rule.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns an `EmpiricalDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[enumerated-real]]."
  (^EmpiricalDistribution [] (empirical nil))
  (^EmpiricalDistribution [{:keys [bin-count data rng]
                            :or {data [1.0]}}]
   (let [^doubles data (m/seq->double-array data)
         bin-count (int (if (number? bin-count) bin-count (bins/freedman-diaconis data (alength data))))
         ^RandomGenerator rng (or rng (JDKRandomGenerator.))
         ^EmpiricalDistribution d (EmpiricalDistribution. bin-count rng)]
     (.load d data)
     d)))

(add-distr-method empirical)

(defn enumerated-real
  "Creates an enumerated real distribution object over an explicit, finite set of double values.

  The distribution puts all of its probability mass on the values listed in `data`. Each value's probability is given by the corresponding entry in `probabilities`, or, when `probabilities` is not supplied, all values are treated as equally likely. Repeated values in `data` accumulate probability mass.

  Parameters (single, optional map):

  - `data` (sequence of doubles): the finite set of values the distribution can take. Default: `[1.0]`.
  - `probabilities` (sequence of doubles): probability associated with each corresponding value in `data`; does not need to be normalized. Default: `nil`, meaning uniform probabilities.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns an `EnumeratedRealDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[enumerated-int]], [[empirical]]."
  (^EnumeratedRealDistribution [] (enumerated-real nil))
  (^EnumeratedRealDistribution [{:keys [data probabilities rng]
                                 :or {data [1.0]}}]
   (let [^RandomGenerator r (or rng (JDKRandomGenerator.))]
     (if probabilities
       (EnumeratedRealDistribution. r (m/seq->double-array data) (m/seq->double-array probabilities))
       (EnumeratedRealDistribution. r ^doubles (m/seq->double-array data))))))

(add-distr-method enumerated-real)

(defn enumerated-int
  "Creates an enumerated integer distribution object over an explicit, finite set of integer values.

  The distribution puts all of its probability mass on the values listed in `data`. Each value's probability is given by the corresponding entry in `probabilities`, or, when `probabilities` is not supplied, all values are treated as equally likely. Repeated values in `data` accumulate probability mass.

  Parameters (single, optional map):

  - `data` (sequence of longs): the finite set of values the distribution can take. Default: `[1]`.
  - `probabilities` (sequence of doubles): probability associated with each corresponding value in `data`; does not need to be normalized. Default: `nil`, meaning uniform probabilities.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns an `EnumeratedIntegerDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[enumerated-real]], [[uniform-int]]."
  (^EnumeratedIntegerDistribution [] (enumerated-int nil))
  (^EnumeratedIntegerDistribution [{:keys [data probabilities rng]
                                    :or {data [1]}}]
   (let [^RandomGenerator r (or rng (JDKRandomGenerator.))]
     (if probabilities
       (EnumeratedIntegerDistribution. r (int-array data) (m/seq->double-array probabilities))
       (EnumeratedIntegerDistribution. r (int-array data))))))

(add-distr-method enumerated-int)

(defn bernoulli
  "Creates a Bernoulli distribution object.

  The Bernoulli distribution is a discrete distribution over the two outcomes `0` and `1`, taking the value `1` with probability `p` and `0` with probability `(- 1.0 p)`. It is the special case of the [[binomial]] distribution with a single trial and is commonly used to model a single yes/no or success/failure event.

  Parameters (single, optional map):

  - `p` (double): probability of success (drawing `1`). Default: `0.5`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `BinomialDistribution` object (with a single trial) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[binomial]], [[geometric]]."
  (^BinomialDistribution [] (bernoulli nil))
  (^BinomialDistribution [{:keys [^double p rng]
                           :or {p 0.5}}]
   (BinomialDistribution. (or rng (JDKRandomGenerator.)) 1 p)))

(add-distr-method bernoulli)

(defn binomial
  "Creates a binomial distribution object.

  The binomial distribution is a discrete distribution over the integers `0` to `trials`, giving the probability of observing a given number of successes out of `trials` independent yes/no experiments, each succeeding with probability `p`. It generalizes the [[bernoulli]] distribution to repeated trials.

  Parameters (single, optional map):

  - `trials` (long): number of independent trials. Default: `20`.
  - `p` (double): probability of success on each trial. Default: `0.5`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `BinomialDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[bernoulli]], [[pascal]], [[hypergeometric]]."
  (^BinomialDistribution [] (binomial nil))
  (^BinomialDistribution [{:keys [^long trials ^double p rng]
                           :or {trials 20 p 0.5}}]
   (BinomialDistribution. (or rng (JDKRandomGenerator.)) trials p)))

(add-distr-method binomial)

(defn geometric
  "Creates a geometric distribution object.

  The geometric distribution is a discrete distribution over the non-negative integers, modelling the number of failures before the first success in a sequence of independent yes/no trials, each succeeding with probability `p`. It is the discrete analogue of the [[exponential]] distribution.

  Parameters (single, optional map):

  - `p` (double): probability of success on each trial. Default: `0.5`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `GeometricDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[bernoulli]], [[pascal]], [[exponential]]."
  (^GeometricDistribution [] (geometric nil))
  (^GeometricDistribution [{:keys [^double p rng]
                            :or {p 0.5}}]
   (GeometricDistribution. (or rng (JDKRandomGenerator.)) p)))

(add-distr-method geometric)

(defn hypergeometric
  "Creates a hypergeometric distribution object.

  The hypergeometric distribution is a discrete distribution modelling the number of successes obtained when drawing `sample-size` elements without replacement from a finite population of `population-size` elements, of which `number-of-successes` are considered successes. Unlike [[binomial]], draws are not independent since sampling is done without replacement.

  Parameters (single, optional map):

  - `population-size` (long): total size of the population being sampled from. Default: `100`.
  - `number-of-successes` (long): number of success elements present in the population. Default: `50`.
  - `sample-size` (long): number of elements drawn from the population without replacement. Default: `25`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `HypergeometricDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[binomial]], [[pascal]]."
  (^HypergeometricDistribution [] (hypergeometric nil))
  (^HypergeometricDistribution [{:keys [^long population-size ^long number-of-successes ^long sample-size rng]
                                 :or {population-size 100 number-of-successes 50 sample-size 25}}]
   (HypergeometricDistribution. (or rng (JDKRandomGenerator.)) population-size number-of-successes sample-size)))

(add-distr-method hypergeometric)

(defn pascal
  "Creates a Pascal distribution object.

  The Pascal distribution is a discrete distribution (a form of the negative binomial distribution) over the non-negative integers, modelling the number of failures observed before accumulating `r` successes in a sequence of independent yes/no trials, each succeeding with probability `p`. It generalizes the [[geometric]] distribution to more than one required success.

  Parameters (single, optional map):

  - `r` (long): number of successes to accumulate. Default: `20`.
  - `p` (double): probability of success on each trial. Default: `0.5`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `PascalDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[geometric]], [[binomial]], [[hypergeometric]]."
  (^PascalDistribution [] (pascal nil))
  (^PascalDistribution [{:keys [^long r ^double p rng]
                         :or {r 20 p 0.5}}]
   (PascalDistribution. (or rng (JDKRandomGenerator.)) r p)))

(add-distr-method pascal)

(defn poisson
  "Creates a Poisson distribution object.

  The Poisson distribution is a discrete distribution over the non-negative integers, modelling the number of events occurring in a fixed interval when events happen independently at a constant average rate. It is commonly used for count data, such as arrivals, defects or occurrences per unit of time or space.

  Parameters (single, optional map):

  - `p` (double): the distribution's rate parameter (mean number of events, traditionally denoted lambda). Default: `0.5`.
  - `epsilon` (double): convergence criterion used internally when computing probabilities. Default: the implementation's default epsilon.
  - `max-iterations` (long): maximum number of iterations used internally when computing probabilities. Default: the implementation's default maximum.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `PoissonDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[binomial]], [[exponential]]."
  (^PoissonDistribution [] (poisson nil))
  (^PoissonDistribution [{:keys [^double p ^double epsilon ^long max-iterations rng]
                          :or {p 0.5 epsilon PoissonDistribution/DEFAULT_EPSILON max-iterations PoissonDistribution/DEFAULT_MAX_ITERATIONS}}]
   (PoissonDistribution. (or rng (JDKRandomGenerator.)) p epsilon max-iterations)))

(add-distr-method poisson)

(defn uniform-int
  "Creates a discrete uniform distribution object over the integers from `lower` to `upper`, inclusive.

  Every integer in the `[lower, upper]` range is equally likely to be drawn. It is the discrete analogue of [[uniform-real]].

  Parameters (single, optional map):

  - `lower` (long): lowest value in the support, inclusive. Default: `0`.
  - `upper` (long): highest value in the support, inclusive. Default: `Integer/MAX_VALUE`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `UniformIntegerDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[uniform-real]], [[enumerated-int]]."
  (^UniformIntegerDistribution [] (uniform-int nil))
  (^UniformIntegerDistribution [{:keys [^long lower ^long upper rng]
                                 :or {lower 0 upper Integer/MAX_VALUE}}]
   (UniformIntegerDistribution. (or rng (JDKRandomGenerator.)) lower upper)))

(add-distr-method uniform-int)

(defn zipf
  "Creates a Zipf distribution object.

  The Zipf distribution is a discrete distribution over the ranks `1` to `number-of-elements`, where the probability of rank `k` is proportional to `k` raised to the negative `exponent`. It is commonly used to model frequency data in which a small number of elements occur very often and the rest occur rarely, such as word frequencies in natural language or city population sizes.

  Parameters (single, optional map):

  - `number-of-elements` (long): number of elements `N`, ie. the highest rank in the distribution support. Default: `100`.
  - `exponent` (double): exponent characterizing the distribution (also known as `s`); higher values make the distribution decay faster for higher ranks. Default: `3.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `ZipfDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[pascal]], [[geometric]], [[hypergeometric]]."
  (^ZipfDistribution [] (zipf nil))
  (^ZipfDistribution [{:keys [^long number-of-elements ^double exponent rng]
                       :or {number-of-elements 100 exponent 3.0}}]
   (ZipfDistribution. (or rng (JDKRandomGenerator.)) number-of-elements exponent)))

(add-distr-method zipf)

(defn multi-normal
  "Creates a multivariate normal (Gaussian) distribution object.

  The multivariate normal distribution is a continuous distribution over `n`-dimensional real vectors, generalizing the [[normal]] distribution to multiple, possibly correlated dimensions. It is fully described by a mean vector and a covariance matrix, and is central to multivariate statistics, sampling correlated random vectors, and many machine learning models.

  Parameters (single, optional map):

  - `means` (sequence of doubles): mean vector, one value per dimension. Default: `[0.0 0.0]`; when only `covariances` is supplied, defaults instead to a zero vector matching its dimension.
  - `covariances` (sequence of sequences of doubles): covariance matrix; must be square, symmetric and positive semi-definite, with a dimension matching `means`. Default: the identity matrix sized to match `means` (a `2x2` identity matrix when neither `means` nor `covariances` is given).
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates a 2-dimensional standard normal distribution with zero means and identity covariance.

  Throws an exception (`ex-info`) if the dimensions of `means` and `covariances` do not match each other.

  Returns a `MultivariateNormalDistribution` object which can be used with [[pdf]], [[lpdf]], [[sample]], [[means]], [[covariance]] and other distribution protocol functions. Unlike univariate distributions, it exposes no `cdf` or `icdf`, and `[[sample]]` returns a vector rather than a scalar.

  See also [[distribution]], [[normal]], [[means]], [[covariance]]."
  ([] (multi-normal nil))
  ([{:keys [means covariances rng]
     :or {means [0.0 0.0]}}]
   (let [covariances (cond
                       (and means (not covariances)) (-> (mat/eye (count means) true)
                                                         (mat/mat->array2d))
                       (not covariances) [[1.0 0.0] [0.0 1.0]]
                       :else covariances)
         means (if-not means (repeat (count (first covariances)) 0.0) means)]
     (when (not= (count means) (count (first covariances)) (count covariances)) (throw (ex-info "Means and covariances sizes do not match."
                                                                                                {:means means :covariances covariances})))
     (MultivariateNormalDistribution. (or rng (JDKRandomGenerator.)) (m/seq->double-array means) (m/seq->double-double-array covariances)))))

(add-distr-method multi-normal)

;; SSJ

(defn anderson-darling
  "Creates a distribution object for the Anderson-Darling goodness-of-fit test statistic.

  Given a sample of `n` independent uniform(0,1) random variables, the Anderson-Darling statistic compares the sorted sample values against the ideal uniform distribution, weighting discrepancies in the tails more heavily than the Kolmogorov-Smirnov statistic does. This distribution describes the sampling distribution of that statistic for a sample of size `n`, and is used to obtain critical values or p-values when testing whether a sample follows a fully specified distribution.

  Parameters (single, optional map):

  - `n` (long): sample size the statistic is computed for. Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `AndersonDarlingDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[anderson-darling-quick]]."
  ([] (anderson-darling nil))
  ([{:keys [^long n rng]
     :or {n 1}}]
   (distr/ssj-continuous :anderson-darling (AndersonDarlingDist. n) rng [:n :rng])))

(add-distr-method anderson-darling)

(defn anderson-darling-quick
  "Creates a distribution object for the Anderson-Darling goodness-of-fit test statistic, using a faster computation algorithm.

  This is a variant of [[anderson-darling]] describing the same underlying statistic for a sample of size `n`, but relying on an alternative, quicker algorithm to evaluate the distribution's functions. It is preferable when the distribution has to be evaluated many times and computation speed matters more than using the original, reference algorithm.

  Parameters (single, optional map):

  - `n` (long): sample size the statistic is computed for. Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `AndersonDarlingDistQuick` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[anderson-darling]]."
  ([] (anderson-darling-quick nil))
  ([{:keys [^long n rng]
     :or {n 1}}]
   (distr/ssj-continuous :anderson-darling-quick (AndersonDarlingDistQuick. n) rng [:n :rng])))

(add-distr-method anderson-darling-quick)

(defn beta-symmetrical
  "Creates a symmetrical beta distribution object.

  The symmetrical beta distribution is a continuous distribution supported on the interval `[0, 1]`, the special case of the [[beta]] distribution where both shape parameters are equal to `alpha`. Being symmetrical, its density is centered and mirrored around `0.5`.

  Parameters (single, optional map):

  - `alpha` (double): shared shape parameter, `alpha = beta` in the general beta distribution. Default: `2.0`.
  - `d` (long): approximate number of decimal digits of precision used internally when computing the distribution, complementary distribution and inverse functions. Default: `14`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `BetaSymmetricalDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[beta]]."
  ([] (beta-symmetrical nil))
  ([{:keys [^double alpha ^long d rng]
     :or {alpha 2.0 d 14}}]
   (distr/ssj-continuous :beta-symmetrical (BetaSymmetricalDist. alpha d) rng [:alpha :d :rng])))

(add-distr-method beta-symmetrical)

(defn chi
  "Creates a chi distribution object.

  The chi distribution is a continuous distribution over non-negative reals, describing the distribution of the square root of a sum of squares of `nu` independent standard normal random variables. It generalizes the [[normal]] distribution's absolute value (for `nu` equal to `1`) and the Rayleigh distribution (for `nu` equal to `2`), and is related to the [[chi-squared]] distribution.

  Parameters (single, optional map):

  - `nu` (long): number of degrees of freedom, ie. the number of underlying standard normal variables. Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `ChiDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[chi-squared]], [[normal]]."
  ([] (chi nil))
  ([{:keys [^long nu rng]
     :or {nu 1}}]
   (distr/ssj-continuous :chi (ChiDist. nu) rng [:nu :rng])))

(add-distr-method chi)

(defn chi-squared-noncentral
  "Creates a noncentral chi-squared distribution object.

  The noncentral chi-squared distribution is a continuous distribution over non-negative reals, generalizing the [[chi-squared]] distribution to the case where the underlying normal random variables have a nonzero, shared mean (encoded through the noncentrality parameter `lambda`). It is used, among others, in power calculations for chi-squared tests and in analyses involving sums of squares of non-centered normal variables.

  Parameters (single, optional map):

  - `nu` (double): number of degrees of freedom. Default: `1.0`.
  - `lambda` (double): noncentrality parameter; `lambda` equal to `0.0` reduces the distribution to a central [[chi-squared]] distribution. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `ChiSquareNoncentralDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[chi-squared]], [[chi]]."
  ([] (chi-squared-noncentral nil))
  ([{:keys [^double nu ^double lambda rng]
     :or {nu 1.0 lambda 1.0}}]
   (distr/ssj-continuous :chi-squared-noncentral (ChiSquareNoncentralDist. nu lambda) rng [:nu :lambda :rng])))

(add-distr-method chi-squared-noncentral)

(defn f-noncentral
  "Creates a noncentral F-distribution object.

  The noncentral F-distribution is a continuous distribution over non-negative reals, generalizing the [[f]] (Fisher-Snedecor) distribution to the case where the numerator's underlying normal random variables have a nonzero, shared mean (encoded through the noncentrality parameter `ncp`): if `X1` follows a [[chi-squared-noncentral]] distribution with `df1` degrees of freedom and noncentrality `ncp`, and `X2` follows an independent central [[chi-squared]] distribution with `df2` degrees of freedom, then `(X1/df1)/(X2/df2)` follows this distribution. It is used, among others, in power calculations for F-tests (ANOVA, regression) under a non-null alternative hypothesis.

  Internally it is computed as a Poisson(`ncp/2`)-weighted mixture of central F-distributions with numerator degrees of freedom `df1 + 2*k` (the same mixture representation that makes [[chi-squared-noncentral]] a Poisson mixture of central [[chi-squared]] distributions), truncated once the cumulative Poisson mass is within `1e-15` of 1.

  Parameters (single, optional map):

  - `df1` (double): numerator degrees of freedom, strictly positive. Default: `1.0`.
  - `df2` (double): denominator degrees of freedom, strictly positive. Default: `1.0`.
  - `ncp` (double): noncentrality parameter, non-negative; `ncp = 0.0` reduces the distribution exactly to a central [[f]] distribution with the same `df1`/`df2`. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Note: `pdf(0)` follows the same boundary convention as the central F-distribution based on `df1`: it is `0` for `df1 > 2`, `e^(-ncp/2)` for `df1 = 2`, and `##Inf` for `df1 < 2`.

  `mean` is finite only for `df2 > 2` (`mean = df2*(df1+ncp) / (df1*(df2-2))`), and `variance` only for `df2 > 4`; outside those ranges they are `##Inf`. There is no closed-form `cdf`/`icdf`; `cdf` sums the Poisson-weighted mixture terms directly (each a call into the well-tested Apache Commons `FDistribution` implementation), and `icdf` root-finds on that `cdf`.

  Matches the `(df1, df2, ncp)` parameterization used by R's base `stats` package (`df`/`pf`/`qf`).

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[f]], [[chi-squared-noncentral]], [[chi-squared]]."
  ([] (f-noncentral nil))
  ([{:keys [^double df1 ^double df2 ^double ncp rng]
     :or {df1 1.0 df2 1.0 ncp 1.0}}]
   (distr/f-noncentral df1 df2 ncp rng)))

(add-distr-method f-noncentral)

(defn t-noncentral
  "Creates a noncentral Student's t-distribution object.

  The noncentral t-distribution is a continuous distribution over the whole real line, generalizing the central [[t]] distribution to the case where the underlying normal numerator has a nonzero mean: if `Z` is a standard normal random variable, `V` an independent chi-squared random variable with `df` degrees of freedom, and `ncp` the noncentrality parameter, then `(Z + ncp) / sqrt(V/df)` follows this distribution. It arises, among others, in power calculations for one- and two-sample t-tests under a non-null alternative hypothesis.

  Internally it is computed by rewriting the defining ratio in terms of `W = sqrt(V)`, which follows a chi distribution with `df` degrees of freedom, and numerically integrating (Gauss-Kronrod quadrature) over `W`: `pdf(x) = integral of chi-pdf(w) * (w/sqrt(df)) * phi(x*w/sqrt(df) - ncp) dw` and `cdf(x) = integral of chi-pdf(w) * Phi(x*w/sqrt(df) - ncp) dw` (`phi`/`Phi` the standard normal density/cdf), both over `w` in `[0, Infinity)`.

  Parameters (single, optional map):

  - `df` (double): degrees of freedom, strictly positive. Default: `1.0`.
  - `ncp` (double): noncentrality parameter; `ncp = 0.0` reduces the distribution exactly to a central [[t]] distribution with the same `df`. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Note: `mean` is finite only for `df > 1` (`mean = ncp*sqrt(df/2)*Gamma((df-1)/2)/Gamma(df/2)`), and `variance` only for `df > 2`; outside those ranges they are `##NaN` (as with the central [[t]] distribution, the moments genuinely do not exist there rather than diverging to infinity). There is no closed-form `cdf`/`icdf`; `cdf` is computed by the quadrature above, and `icdf` root-finds on that `cdf`.

  Matches the `(df, ncp)` parameterization used by R's base `stats` package (`dt`/`pt`/`qt`).

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[t]], [[f-noncentral]], [[chi-squared-noncentral]]."
  ([] (t-noncentral nil))
  ([{:keys [^double df ^double ncp rng]
     :or {df 1.0 ncp 1.0}}]
   (distr/t-noncentral df ncp rng)))

(add-distr-method t-noncentral)

(defn cramer-von-mises
  "Creates a distribution object for the Cramer-von Mises goodness-of-fit test statistic.

  Given a sample of `n` independent uniform(0,1) random variables, the Cramer-von Mises statistic measures the squared distance between the empirical distribution function of the sorted sample and the ideal uniform distribution. This distribution describes the sampling distribution of that statistic for a sample of size `n`, and is used to obtain critical values or p-values when testing whether a sample follows a fully specified distribution. Unlike most other SSJ-backed distributions here, it exposes no closed-form density, so its `pdf` is approximated numerically from the cumulative distribution function.

  Parameters (single, optional map):

  - `n` (long): sample size the statistic is computed for. Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `CramerVonMisesDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[anderson-darling]], [[anderson-darling-quick]]."
  ([] (cramer-von-mises nil))
  ([{:keys [^long n rng]
     :or {n 1}}]
   (distr/ssj-continuous-no-pdf :cramer-von-mises (CramerVonMisesDist. n) rng [:n :rng])))

(add-distr-method cramer-von-mises)

(defn erlang
  "Creates an Erlang distribution object.

  The Erlang distribution is a continuous distribution over non-negative reals, describing the sum of `k` independent [[exponential]] random variables each with rate `lambda`. It is the special case of the [[gamma]] distribution restricted to an integer shape parameter, and commonly models waiting times for multiple sequential events, such as queueing systems.

  Parameters (single, optional map):

  - `k` (long): number of exponential stages summed together (shape parameter); must be a positive integer. Default: `2`.
  - `lambda` (double): rate parameter of each underlying exponential stage. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `ErlangDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[exponential]], [[gamma]]."
  ([] (erlang nil))
  ([{:keys [^long k ^double lambda rng]
     :or {k 2 lambda 1.0}}]
   (distr/ssj-continuous :erlang (ErlangDist. k lambda) rng [:k :lambda :rng])))

(add-distr-method erlang)

(defn fatigue-life
  "Creates a fatigue life (Birnbaum-Saunders) distribution object.

  The fatigue life distribution is a continuous distribution over reals greater than `alpha`, originally derived to model the time to failure of materials under cyclic stress caused by crack growth. It is constructed from a transformation of a standard normal random variable and is widely used in reliability and survival analysis.

  Parameters (single, optional map):

  - `alpha` (double): location parameter, the lower bound of the support. Default: `0.0`.
  - `beta` (double): scale parameter; must be positive. Default: `1.0`.
  - `gamma` (double): shape parameter; must be positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `FatigueLifeDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[weibull]], [[normal]]."
  ([] (fatigue-life nil))
  ([{:keys [^double alpha ^double beta ^double gamma rng]
     :or {alpha 0.0 beta 1.0 gamma 1.0}}]
   (distr/ssj-continuous :fatigue-life (FatigueLifeDist. alpha beta gamma) rng [:alpha :beta :gamma :rng])))

(add-distr-method fatigue-life)

(defn folded-normal
  "Creates a folded normal distribution object.

  The folded normal distribution is a continuous distribution over non-negative reals, describing the distribution of the absolute value of a [[normal]] random variable with mean `mu` and standard deviation `sigma`. It is used, among others, to model measurement magnitudes or absolute errors when only the size, not the sign, of an underlying normally distributed quantity can be observed.

  Parameters (single, optional map):

  - `mu` (double): mean of the underlying normal distribution, before folding. Default: `0.0`.
  - `sigma` (double): standard deviation of the underlying normal distribution; must be positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `FoldedNormalDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[normal]], [[half-normal]]."
  ([] (folded-normal nil))
  ([{:keys [^double mu ^double sigma rng]
     :or {mu 0.0 sigma 1.0}}]
   (distr/ssj-continuous :folded-normal (FoldedNormalDist. mu sigma) rng [:mu :sigma :rng])))

(add-distr-method folded-normal)

(defn frechet
  "Creates a Frechet distribution object.

  The Frechet distribution is a continuous, heavy-tailed distribution supported on `(delta, Infinity)`, one of the three families of extreme value distributions (the type II extreme value distribution). It commonly models the maximum of a number of samples, such as extreme rainfall, flood levels or maximum returns in finance.

  Parameters (single, optional map):

  - `alpha` (double): shape parameter; must be positive. Default: `1.0`.
  - `beta` (double): scale parameter; must be positive. Default: `1.0`.
  - `delta` (double): location parameter, the lower bound of the support. Default: `0.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `FrechetDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[gumbel]], [[weibull]], [[pareto]]."
  ([] (frechet nil))
  ([{:keys [^double alpha ^double beta ^double delta rng]
     :or {alpha 1.0 beta 1.0 delta 0.0}}]
   (distr/ssj-continuous :frechet (FrechetDist. alpha beta delta) rng [:alpha :beta :delta :rng])))

(add-distr-method frechet)

(defn half-normal
  "Creates a half-normal distribution object.

  The half-normal distribution is a continuous distribution over reals not less than `mu`, obtained by folding a [[normal]] distribution at its location `mu` and keeping only the non-negative half. It is the one-sided special case of the [[folded-normal]] distribution and is often used to model magnitudes or absolute errors that cannot be negative.

  Parameters (single, optional map):

  - `mu` (double): location parameter, the lower bound of the support. Default: `0.0`.
  - `sigma` (double): scale parameter, related to the spread of the distribution; must be positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `HalfNormalDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[folded-normal]], [[normal]], [[chi]]."
  ([] (half-normal nil))
  ([{:keys [^double mu ^double sigma rng]
     :or {mu 0.0 sigma 1.0}}]
   (distr/ssj-continuous :half-normal (HalfNormalDist. mu sigma) rng [:mu :sigma :rng])))

(add-distr-method half-normal)

(defn hyperbolic-secant
  "Creates a hyperbolic secant distribution object.

  The hyperbolic secant distribution is a continuous, symmetric, bell-shaped distribution whose density is proportional to the hyperbolic secant of the standardized variable. Its shape lies between the [[normal]] and [[cauchy]] distributions: it is more peaked and has slightly heavier tails than the normal distribution, while remaining much lighter-tailed than the Cauchy distribution.

  Parameters (single, optional map):

  - `mu` (double): location parameter, the center of the distribution. Default: `0.0`.
  - `sigma` (double): scale parameter, controlling the spread of the distribution; must be positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `HyperbolicSecantDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[normal]], [[logistic]], [[cauchy]]."
  ([] (hyperbolic-secant nil))
  ([{:keys [^double mu ^double sigma rng]
     :or {mu 0.0 sigma 1.0}}]
   (distr/ssj-continuous :hyperbolic-secant (HyperbolicSecantDist. mu sigma) rng [:mu :sigma :rng])))

(add-distr-method hyperbolic-secant)

(defn hypoexponential-equal
  "Creates a hypoexponential distribution object with equally spaced rates.

  The hypoexponential distribution describes the sum of independent [[exponential]] random variables (phases), possibly with different rates. This variant models the sum of `k` phases whose rates are equally spaced with common difference `h`, out of `n` available equally spaced rates (with `n` at least `k`); it is a compact way to specify such a distribution without listing every individual rate, useful for approximating other positive-valued distributions in queueing and reliability models.

  Parameters (single, optional map):

  - `n` (long): total number of equally spaced rates the distribution is built from; must be at least `k`. Default: `1`.
  - `k` (long): number of phases, out of the `n` available, that are summed together. Default: `1`.
  - `h` (double): common spacing between consecutive rates. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `HypoExponentialDistEqual` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[hypoexponential]], [[erlang]], [[exponential]]."
  ([] (hypoexponential-equal nil))
  ([{:keys [^long n ^long k ^double h rng]
     :or {n 1 k 1 h 1.0}}]
   (distr/ssj-continuous :hypoexponential-equal (HypoExponentialDistEqual. n k h) rng [:n :k :h :rng])))

(add-distr-method hypoexponential-equal)

(defn hypoexponential
  "Creates a hypoexponential distribution object with explicit rates.

  The hypoexponential distribution is a continuous distribution over non-negative reals, describing the sum of independent [[exponential]] random variables, one for each rate listed in `lambdas`. Unlike the [[erlang]] distribution, the rates need not be equal, which lets the distribution capture more general shapes for multi-stage processes such as sequential tasks or phase-type service times.

  Parameters (single, optional map):

  - `lambdas` (sequence of doubles): rate of each independent exponential phase summed together. Default: `[1.0]`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `HypoExponentialDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[hypoexponential-equal]], [[erlang]], [[exponential]]."
  ([] (hypoexponential nil))
  ([{:keys [lambdas rng]
     :or {lambdas [1.0]}}]
   (distr/ssj-continuous :hypoexponential (HypoExponentialDist. (m/seq->double-array lambdas))
                         rng [:lambdas :rng])))

(add-distr-method hypoexponential)

(defn inverse-gamma
  "Creates an inverse gamma distribution object.

  The inverse gamma distribution is a continuous distribution over positive reals, describing the distribution of the reciprocal of a [[gamma]]-distributed random variable. It is commonly used as a conjugate prior for the variance parameter of a normal distribution in Bayesian statistics.

  Parameters (single, optional map):

  - `alpha` (double): shape parameter; must be positive. Default: `1.0`.
  - `beta` (double): scale parameter; must be positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `InverseGammaDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[gamma]]."
  ([] (inverse-gamma nil))
  ([{:keys [^double alpha ^double beta rng]
     :or {alpha 1.0 beta 1.0}}]
   (distr/ssj-continuous :inverse-gamma (InverseGammaDist. alpha beta) rng [:alpha :beta :rng])))

(add-distr-method inverse-gamma)

(defn inverse-gaussian
  "Creates an inverse Gaussian (Wald) distribution object.

  The inverse Gaussian distribution is a continuous distribution over positive reals, describing, among other interpretations, the first passage time of a Brownian motion with positive drift towards a fixed positive threshold. Despite its name, it is not obtained by inverting the [[normal]] distribution's cdf, but rather because its cumulant generating function is the functional inverse of that of the normal distribution.

  Parameters (single, optional map):

  - `mu` (double): mean of the distribution; must be positive. Default: `1.0`.
  - `lambda` (double): shape parameter, controlling how concentrated the distribution is around its mean; must be positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `InverseGaussianDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[gamma]], [[normal]]."
  ([] (inverse-gaussian nil))
  ([{:keys [^double mu ^double lambda rng]
     :or {mu 1.0 lambda 1.0}}]
   (distr/ssj-continuous :inverse-gaussian (InverseGaussianDist. mu lambda) rng [:mu :lambda :rng])))

(add-distr-method inverse-gaussian)

(defn johnson-sb
  "Creates a Johnson SB (bounded) distribution object.

  The Johnson SB distribution is a continuous, four-parameter distribution supported on the bounded interval `(xi, xi + lambda)`, obtained by applying a logit-like transformation to a standard normal random variable. It belongs to the flexible Johnson system of distributions, alongside [[johnson-sl]] and [[johnson-su]], and can approximate a wide range of bounded, possibly skewed shapes by fitting its shape parameters to observed data.

  Parameters (single, optional map):

  - `gamma` (double): first shape parameter, controlling skewness. Default: `0.0`.
  - `delta` (double): second shape parameter; must be positive. Default: `1.0`.
  - `xi` (double): location parameter, the lower bound of the support. Default: `0.0`.
  - `lambda` (double): scale parameter, the width of the support; must be positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `JohnsonSBDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[johnson-sl]], [[johnson-su]], [[beta]]."
  ([] (johnson-sb nil))
  ([{:keys [^double gamma ^double delta ^double xi ^double lambda rng]
     :or {gamma 0.0 delta 1.0 xi 0.0 lambda 1.0}}]
   (distr/ssj-continuous :johnson-sb (JohnsonSBDist. gamma delta xi lambda) rng [:gamma :delta :xi :lambda :rng])))

(add-distr-method johnson-sb)

(defn johnson-sl
  "Creates a Johnson SL (log-normal, semi-bounded) distribution object.

  The Johnson SL distribution is a continuous, four-parameter distribution supported on the semi-bounded interval `(xi, Infinity)`, obtained by applying a logarithmic transformation to a standard normal random variable; it is equivalent to a shifted and rescaled [[log-normal]] distribution. It belongs to the Johnson system of distributions, alongside [[johnson-sb]] and [[johnson-su]].

  Parameters (single, optional map):

  - `gamma` (double): first shape parameter, controlling skewness. Default: `0.0`.
  - `delta` (double): second shape parameter; must be positive. Default: `1.0`.
  - `xi` (double): location parameter, the lower bound of the support. Default: `0.0`.
  - `lambda` (double): scale parameter; must be positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `JohnsonSLDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[johnson-sb]], [[johnson-su]], [[log-normal]]."
  ([] (johnson-sl nil))
  ([{:keys [^double gamma ^double delta ^double xi ^double lambda rng]
     :or {gamma 0.0 delta 1.0 xi 0.0 lambda 1.0}}]
   (distr/ssj-continuous :johnson-sl (JohnsonSLDist. gamma delta xi lambda) rng [:gamma :delta :xi :lambda :rng])))

(add-distr-method johnson-sl)

(defn johnson-su
  "Creates a Johnson SU (unbounded) distribution object.

  The Johnson SU distribution is a continuous, four-parameter distribution supported on the whole real line, obtained by applying an inverse hyperbolic sine transformation to a standard normal random variable. It belongs to the Johnson system of distributions, alongside [[johnson-sb]] and [[johnson-sl]], and its flexible shape makes it useful for fitting skewed or heavy-tailed data that a plain [[normal]] distribution cannot capture.

  Parameters (single, optional map):

  - `gamma` (double): first shape parameter, controlling skewness. Default: `0.0`.
  - `delta` (double): second shape parameter; must be positive. Default: `1.0`.
  - `xi` (double): location parameter. Default: `0.0`.
  - `lambda` (double): scale parameter; must be positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `JohnsonSUDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[johnson-sb]], [[johnson-sl]], [[normal]]."
  ([] (johnson-su nil))
  ([{:keys [^double gamma ^double delta ^double xi ^double lambda rng]
     :or {gamma 0.0 delta 1.0 xi 0.0 lambda 1.0}}]
   (distr/ssj-continuous :johnson-su (JohnsonSUDist. gamma delta xi lambda) rng [:gamma :delta :xi :lambda :rng])))

(add-distr-method johnson-su)

(defn kolmogorov-smirnov
  "Creates a distribution object for the (two-sided) Kolmogorov-Smirnov goodness-of-fit test statistic.

  Given a sample of `n` independent uniform(0,1) random variables, the Kolmogorov-Smirnov statistic measures the largest absolute distance between the empirical distribution function of the sorted sample and the ideal uniform distribution. This distribution describes the sampling distribution of that statistic for a sample of size `n`, and is used to obtain critical values or p-values when testing whether a sample follows a fully specified distribution.

  Parameters (single, optional map):

  - `n` (long): sample size the statistic is computed for. Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `KolmogorovSmirnovDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[kolmogorov-smirnov+]], [[kolmogorov-smirnov-quick]], [[anderson-darling]], [[cramer-von-mises]]."
  ([] (kolmogorov-smirnov nil))
  ([{:keys [^long n rng]
     :or {n 1}}]
   (distr/ssj-continuous :kolmogorov-smirnov (KolmogorovSmirnovDist. n) rng [:n :rng])))

(add-distr-method kolmogorov-smirnov)

(defn kolmogorov-smirnov+
  "Creates a distribution object for the one-sided Kolmogorov-Smirnov (D+) goodness-of-fit test statistic.

  This is the one-sided variant of [[kolmogorov-smirnov]], describing the sampling distribution of the largest positive (rather than absolute) deviation between the empirical distribution function of a sample of `n` independent uniform(0,1) variables and the ideal uniform distribution. It is sensitive to departures where the empirical distribution lies above the hypothesized one.

  Parameters (single, optional map):

  - `n` (long): sample size the statistic is computed for. Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `KolmogorovSmirnovPlusDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[kolmogorov-smirnov]], [[kolmogorov-smirnov-quick]]."
  ([] (kolmogorov-smirnov+ nil))
  ([{:keys [^long n rng]
     :or {n 1}}]
   (distr/ssj-continuous :kolmogorov-smirnov+ (KolmogorovSmirnovPlusDist. n) rng [:n :rng])))

(add-distr-method kolmogorov-smirnov+)

(defn kolmogorov-smirnov-quick
  "Creates a distribution object for the (two-sided) Kolmogorov-Smirnov goodness-of-fit test statistic, using a faster computation algorithm.

  This is a variant of [[kolmogorov-smirnov]] describing the same underlying statistic for a sample of size `n`, but relying on an alternative, quicker algorithm to evaluate the distribution's functions, similar in spirit to [[anderson-darling-quick]].

  Parameters (single, optional map):

  - `n` (long): sample size the statistic is computed for. Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `KolmogorovSmirnovDistQuick` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[kolmogorov-smirnov]], [[kolmogorov-smirnov+]]."
  ([] (kolmogorov-smirnov-quick nil))
  ([{:keys [^long n rng]
     :or {n 1}}]
   (distr/ssj-continuous :kolmogorov-smirnov-quick (KolmogorovSmirnovDistQuick. n) rng [:n :rng])))

(add-distr-method kolmogorov-smirnov-quick)

(defn kolmogorov
  "Creates a Kolmogorov distribution object.

  The Kolmogorov distribution is a continuous, parameter-free distribution over positive reals equal to the limiting distribution of `(* (m/sqrt n) (kolmogorov-smirnov n))` as `n` grows to infinity, that is, of the scaled two-sided Kolmogorov-Smirnov statistic; equivalently, it is the distribution of the supremum of the absolute value of a standard Brownian bridge on `[0, 1]`. It underlies the asymptotic (large-sample) Kolmogorov-Smirnov goodness-of-fit test, used when the exact, finite-sample distribution given by [[kolmogorov-smirnov]] is impractical to compute.

  Parameters (single, optional map):

  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[kolmogorov-smirnov]], [[kolmogorov-smirnov-quick]], [[kolmogorov-smirnov+]]."
  ([] (kolmogorov nil))
  ([{:keys [rng]}]
   (distr/kolmogorov rng)))

(add-distr-method kolmogorov)

(defn log-logistic
  "Creates a log-logistic (Fisk) distribution object.

  The log-logistic distribution is a continuous distribution over positive reals whose logarithm follows a [[logistic]] distribution. It resembles the [[log-normal]] and [[weibull]] distributions in shape and is commonly used in survival analysis and to model income and other positively-skewed data, since, unlike the Weibull distribution, its hazard function can be non-monotonic.

  Parameters (single, optional map):

  - `alpha` (double): scale parameter; must be positive. Default: `1.0`.
  - `beta` (double): shape parameter; must be positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `LoglogisticDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[logistic]], [[log-normal]], [[weibull]]."
  ([] (log-logistic nil))
  ([{:keys [^double alpha ^double beta rng]
     :or {alpha 1.0 beta 1.0}}]
   (distr/ssj-continuous :log-logistic (LoglogisticDist. alpha beta) rng [:alpha :beta :rng])))

(add-distr-method log-logistic)

(defn normal-inverse-gaussian
  "Creates a normal-inverse Gaussian distribution object.

  The normal-inverse Gaussian distribution is a continuous, heavy-tailed distribution over the whole real line, constructed as a normal-variance mixture where the mixing variance follows an [[inverse-gaussian]] distribution. It can flexibly model both skewness and excess kurtosis and is widely used in finance to model asset returns.

  Parameters (single, optional map):

  - `alpha` (double): tail heaviness parameter; must be positive. Default: `1.0`.
  - `beta` (double): asymmetry parameter; must satisfy `(< (m/abs beta) alpha)`. Default: `0.0`.
  - `mu` (double): location parameter. Default: `0.0`.
  - `delta` (double): scale parameter; must be positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  `pdf` is closed-form; `cdf`/`icdf` are obtained by numerically integrating `pdf` (this is exactly the `lambda = -0.5` special case of [[generalized-hyperbolic]], and reuses its machinery internally), and `mean`/`variance` have closed forms.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[inverse-gaussian]], [[normal]], [[generalized-hyperbolic]]."
  ([] (normal-inverse-gaussian nil))
  ([{:keys [^double alpha ^double beta ^double mu ^double delta rng]
     :or {alpha 1.0 beta 0.0 mu 0.0 delta 1.0}}]
   (distr/normal-inverse-gaussian {:alpha alpha :beta beta :mu mu :delta delta :rng rng})))

(add-distr-method normal-inverse-gaussian)

(defn pearson-6
  "Creates a Pearson type VI distribution object.

  The Pearson type VI distribution (also known as the scaled beta distribution of the second kind) is a continuous distribution over positive reals, shaped by two shape parameters `alpha1` and `alpha2` and a scale parameter `beta`. It is related to the [[f]] and [[beta]] distributions through a change of variables, and is used, among others, to model positive, right-skewed quantities.

  Parameters (single, optional map):

  - `alpha1` (double): first shape parameter; must be positive. Default: `1.0`.
  - `alpha2` (double): second shape parameter; must be positive. Default: `1.0`.
  - `beta` (double): scale parameter; must be positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `Pearson6Dist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[beta]], [[f]], [[gamma]]."
  ([] (pearson-6 nil))
  ([{:keys [^double alpha1 ^double alpha2 ^double beta rng]
     :or {alpha1 1.0 alpha2 1.0 beta 1.0}}]
   (distr/ssj-continuous :pearson-6 (Pearson6Dist. alpha1 alpha2 beta) rng
                         [:alpha1 :alpha2 :beta :rng])))

(add-distr-method pearson-6)

(defn power
  "Creates a power-function distribution object.

  The power-function distribution is a continuous distribution over the bounded interval `[a, b]`, whose cumulative distribution function grows as the `c`-th power of the normalized position within the interval. It generalizes the [[uniform-real]] distribution, recovered when `c` is equal to `1.0`, and provides a simple way to model quantities concentrated towards one end of a bounded range.

  Parameters (single, optional map):

  - `a` (double): lower bound of the support. Default: `0.0`.
  - `b` (double): upper bound of the support. Default: `1.0`.
  - `c` (double): shape parameter, controlling how the density concentrates near `b` (for values greater than `1.0`) or near `a` (for values less than `1.0`); must be positive. Default: `2.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `PowerDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[uniform-real]], [[triangular]]."
  ([] (power nil))
  ([{:keys [^double a ^double b ^double c rng]
     :or {a 0.0 b 1.0 c 2.0}}]
   (distr/ssj-continuous :power (PowerDist. a b c) rng
                         [:a :b :c :rng])))

(add-distr-method power)

(defn rayleigh
  "Creates a Rayleigh distribution object.

  The Rayleigh distribution is a continuous distribution over reals not less than `a`, arising as the magnitude of a two-dimensional vector whose components are independent, identically-distributed, zero-mean normal random variables. It is commonly used to model wind speeds, wave heights and the magnitude of radio signals affected by multipath fading.

  Parameters (single, optional map):

  - `a` (double): location parameter, the lower bound of the support. Default: `0.0`.
  - `beta` (double): scale parameter; must be positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `RayleighDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[chi]], [[weibull]], [[normal]]."
  ([] (rayleigh nil))
  ([{:keys [^double a ^double beta rng]
     :or {a 0.0 beta 1.0}}]
   (distr/ssj-continuous :rayleigh (RayleighDist. a beta) rng
                         [:a :beta :rng])))

(add-distr-method rayleigh)

(defn watson-g
  "Creates a distribution object for the Watson G goodness-of-fit test statistic.

  The Watson G statistic is a variant of the [[kolmogorov-smirnov]] statistic adjusted to be invariant to the choice of origin, making it especially suited for testing goodness-of-fit on circular (directional) data, such as angles or times of day. This distribution describes the sampling distribution of that statistic for a sample of size `n`.

  Parameters (single, optional map):

  - `n` (long): sample size the statistic is computed for; must be at least `2` (the backing SSJ class requires it). Default: `2`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `WatsonGDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[watson-u]], [[kolmogorov-smirnov]]."
  ([] (watson-g nil))
  ([{:keys [^long n rng]
     :or {n 2}}]
   (distr/ssj-continuous :watson-g (WatsonGDist. n) rng [:n :rng])))

(add-distr-method watson-g)

(defn watson-u
  "Creates a distribution object for the Watson U-squared goodness-of-fit test statistic.

  The Watson U-squared statistic is a variant of the [[cramer-von-mises]] statistic adjusted to be invariant to the choice of origin, making it especially suited for testing goodness-of-fit on circular (directional) data, such as angles or times of day. This distribution describes the sampling distribution of that statistic for a sample of size `n`.

  Parameters (single, optional map):

  - `n` (long): sample size the statistic is computed for; must be at least `2` (the backing SSJ class requires it). Default: `2`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `WatsonUDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[watson-g]], [[cramer-von-mises]]."
  ([] (watson-u nil))
  ([{:keys [^long n rng]
     :or {n 2}}]
   (distr/ssj-continuous :watson-u (WatsonUDist. n) rng [:n :rng])))

(add-distr-method watson-u)

(defn von-mises
  "Creates a von Mises distribution object.

  The von Mises distribution is the circular (directional) analogue of the [[normal]] distribution: a continuous distribution over angles, symmetric and unimodal around a mean direction `mu`, with concentration controlled by `kappa` (the higher `kappa`, the tighter the distribution clusters around `mu`; `kappa = 0` gives the uniform distribution on the circle, and for large `kappa` it approaches a normal distribution with standard deviation `1/sqrt(kappa)`). It is widely used to model angular data such as wind directions, compass bearings, or times of day.

  Its density is `f(x) = e^(kappa*cos(x-mu)) / (2*pi*I_0(kappa))`, where `I_0` is the modified Bessel function of the first kind of order 0.

  Parameters (single, optional map):

  - `mu` (double): mean direction, in radians. Default: `0.0`.
  - `kappa` (double): concentration parameter, non-negative. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Note: since the underlying quantity is angular (periodic with period `2*pi`), the distribution is represented here on the principal branch `[mu-pi, mu+pi]`: `pdf`/`cdf` are `0`/`0`/`1` outside that range (the usual convention for every other bounded distribution in this library), rather than wrapping `x` back into range as e.g. R's `circular` package's `dvonmises`/`pvonmises` do. `mean` is exactly `mu` (by symmetry); `variance` is the ordinary linear `E[(X-mu)^2]` restricted to `[mu-pi, mu+pi]` (computed numerically, no closed form) - this is *not* the same as the circular-statistics notion of circular variance (`1 - I_1(kappa)/I_0(kappa)`), which is a different, bounded-in-`[0,1]` quantity.

  For large `kappa` the density concentrates into an increasingly narrow peak around `mu` (characteristic width `~1/sqrt(kappa)`); `pdf`'s normalizing constant, and the internal numerics behind `cdf`/`icdf`/`variance`, are all scaled with `kappa` so that they stay accurate (verified against direct numerical integration) even for very large `kappa`, rather than silently losing precision (or overflowing) once the peak becomes narrower than a fixed-resolution scheme would resolve.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[normal]], [[cramer-von-mises]], [[watson-g]], [[watson-u]]."
  ([] (von-mises nil))
  ([{:keys [^double mu ^double kappa rng]
     :or {mu 0.0 kappa 1.0}}]
   (distr/von-mises mu kappa rng)))

(add-distr-method von-mises)

(defn multinomial
  "Creates a multinomial distribution object.

  The multinomial distribution is a discrete, multivariate distribution generalizing the [[binomial]] distribution to more than two possible outcomes: it describes the counts of each of several categories obtained from `n` independent trials, where each trial falls into category `i` with probability `ps[i]`. It is commonly used to model counts of outcomes across multiple categories, such as votes for several candidates or colors drawn from a bag.

  Parameters (single, optional map):

  - `n` (long): number of independent trials distributed among the categories. Default: `10`.
  - `ps` (sequence of doubles): relative probability of each category; does not need to sum to `1.0`, since values are normalized internally. Default: `[1 1]`, two equally likely categories (paired with [[dirichlet]]'s own default).
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `MultinomialDist` class), with as many dimensions as there are entries in `ps`, which can be used with [[pdf]], [[cdf]], [[sample]], [[means]], [[covariance]] and other distribution protocol functions. Each value returned by `sample` (and accepted by `pdf`/`cdf`) is a vector of non-negative integer counts summing to `n`, one per category. Unlike univariate distributions, it exposes no `icdf`, `mean` or `variance`.

  See also [[distribution]], [[binomial]], [[multi-normal]]."
  ([] (multinomial nil))
  ([{:keys [^long n ps rng]
     :or {n 10 ps [1 1]}}]
   (distr/multinomial n ps binomial rng)))

(add-distr-method multinomial)

;; custom

(defn dirichlet
  "Creates a Dirichlet distribution object.

  The Dirichlet distribution is a continuous, multivariate distribution over the open probability simplex, ie. over vectors of positive values summing to `1.0`, generalizing the [[beta]] distribution to more than two categories. It is commonly used as a prior over probability vectors in Bayesian statistics, for example over the category probabilities of a [[multinomial]] distribution.

  Parameters (single, optional map):

  - `alpha` (sequence of doubles, or long): concentration parameter for each dimension of the simplex. Larger values concentrate samples closer to the uniform point, smaller values push mass towards the corners of the simplex. When given as a plain integer `n`, it is treated as a symmetric Dirichlet distribution with `n` dimensions, each with concentration `1.0` (ie. uniform over the simplex). Default: `[1 1]`, ie. uniform over the 2-dimensional simplex.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object, with as many dimensions as entries in `alpha`, which can be used with [[pdf]], [[lpdf]], [[sample]], [[means]], [[covariance]] and other distribution protocol functions. Sampled vectors, and the values accepted by `pdf`/`lpdf`, are vectors of positive doubles summing (approximately) to `1.0`. Unlike univariate distributions, it exposes no `cdf`, `icdf`, `mean` or `variance`.

  See also [[distribution]], [[beta]], [[multinomial]], [[gamma]]."
  ([] (dirichlet nil))
  ([{:keys [alpha rng]
     :or {alpha [1 1]}}]
   (distr/dirichlet alpha gamma rng)))

(add-distr-method dirichlet)

(defn continuous-distribution
  "Creates a continuous distribution object estimated from a sample of `data` using kernel density estimation (KDE).

  Rather than assuming a parametric family, this builds a smooth, nonparametric estimate of the underlying density from `data` using a chosen kernel and bandwidth (see [[fastmath.kernel.density]]), then numerically integrates it to obtain matching `cdf` and `icdf` functions. It is useful for approximating an unknown continuous distribution directly from observations, offering a smoother alternative to [[empirical]].

  Parameters (single, optional map):

  - `data` (sequence of doubles): sample used to estimate the density. Default: `[-1 0 1]`.
  - `kde` (keyword or function): kernel used for the density estimation, either a keyword naming a kernel (eg. `:epanechnikov`, `:gaussian`) or a custom kernel function; see [[fastmath.kernel.density/kernel-density]]. Default: `:epanechnikov`.
  - `bandwidth` (double or keyword): bandwidth used by the kernel density estimator; either a fixed positive number or one of `:nrd`, `:nrd0`, `:nrd-adjust`, `:rlcv`, `:lcv`, `:lscv` to infer it from `data`. Default: inferred automatically.
  - `steps` (long): number of subintervals used when numerically integrating the density into a `cdf`/`icdf` pair, see [[integrate-pdf]]. Default: `5000`.
  - `interpolator` (keyword or function): interpolation method used between the integrated `cdf`/`icdf` points; one of `:linear`, `:cubic`, `:monotone`, or a custom function. Default: `:linear`.
  - `rng`: random number generator used for sampling.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[kde]], [[empirical]], [[enumerated-real]]."
  ([] (continuous-distribution nil))
  ([options]
   (distr/continuous-distribution (merge {:data [-1 0 1] :steps 5000 :kde :epanechnikov :rng (:rng options)} options))))

(add-distr-method continuous-distribution)

(defn kde
  "Creates a continuous distribution object estimated from a sample of data using kernel density estimation (KDE).

  Alias for [[continuous-distribution]]; see its docstring for the accepted parameters (`:data`, `:kde`, `:bandwidth`, `:steps`, `:interpolator`, `:rng`) and their defaults.

  See also [[distribution]], [[continuous-distribution]], [[empirical]]."
  ([] (kde nil))
  ([options] (continuous-distribution options)))

(add-distr-method kde)

(defn negative-binomial
  "Creates a negative binomial distribution object.

  The negative binomial distribution is a discrete distribution over the non-negative integers, modelling the number of failures observed before accumulating `r` successes in a sequence of independent yes/no trials, each succeeding with probability `p`. Unlike [[pascal]], which is restricted to an integer number of successes, `r` here may be any positive real number (a generalized, or Polya, negative binomial distribution).

  Parameters (single, optional map):

  - `r` (double): number of successes to accumulate; need not be an integer. Default: `20`.
  - `p` (double): probability of success on each trial. Default: `0.5`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[pascal]], [[geometric]], [[binomial]]."
  ([] (negative-binomial nil))
  ([{:keys [^long r ^double p rng]
     :or {r 20 p 0.5}}]
   (distr/negative-binomial r p rng)))

(add-distr-method negative-binomial)

(defn nbi
  "Creates a negative binomial (type I) distribution object, using gamlss-style parameter names.

  Alias for [[negative-binomial]], reparameterized in terms of the mean `mu` and a dispersion parameter `sigma`, matching the naming and parametrization used by R's `gamlss.dist` package (its `NBI` family, also known as the `NB1` parametrization). The underlying [[negative-binomial]] parameters are recovered as `r = 1 / sigma` and `p = r / (r + mu)`, giving mean `mu` and variance `(+ mu (* sigma mu mu))`, i.e. a variance that is linear in `mu`.

  Parameters (single, optional map):

  - `mu` (double): mean of the distribution, strictly positive. Default: `1.0`.
  - `sigma` (double): dispersion parameter, strictly positive; larger values give more overdispersion relative to a [[poisson]] distribution with the same mean, which is approached as `sigma` tends to `0`. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  When `sigma` is below `1.0e-4`, `r = 1 / sigma` grows so large that evaluating [[negative-binomial]] directly loses numerical precision; in that regime this falls back to an exact [[poisson]] distribution with rate `mu`, which is the limit of `nbi` as `sigma` tends to `0`.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[negative-binomial]], [[nbii]], [[poisson]]."
  ([] (nbi nil))
  ([{:keys [^double mu ^double sigma rng]
     :or {mu 1.0 sigma 1.0}}]
   (if (m/< sigma 1.0e-4)
     (poisson {:p mu :rng rng})
     (let [r (m// 1.0 sigma)
           p (m// r (m/+ r mu))]
       (distr/negative-binomial r p rng)))))

(add-distr-method nbi)

(defn nbii
  "Creates a negative binomial (type II) distribution object, using gamlss-style parameter names.

  Alias for [[negative-binomial]], reparameterized in terms of the mean `mu` and a dispersion parameter `sigma`, matching the naming and parametrization used by R's `gamlss.dist` package (its `NBII` family). The underlying [[negative-binomial]] parameters are recovered as `r = mu / sigma` and `p = r / (r + mu)`, giving mean `mu` and variance `(* mu (+ 1.0 sigma))`, i.e. a variance that is linear in the dispersion parameter, unlike [[nbi]] whose variance is quadratic in `mu`.

  Parameters (single, optional map):

  - `mu` (double): mean of the distribution, strictly positive. Default: `1.0`.
  - `sigma` (double): dispersion parameter, strictly positive; larger values give more overdispersion relative to a [[poisson]] distribution with the same mean, which is approached as `sigma` tends to `0`. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  When `sigma` is below `1.0e-4`, `r = mu / sigma` grows so large that evaluating [[negative-binomial]] directly loses numerical precision; in that regime this falls back to an exact [[poisson]] distribution with rate `mu`, which is the limit of `nbii` as `sigma` tends to `0`.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[negative-binomial]], [[nbi]], [[poisson]]."
  ([] (nbii nil))
  ([{:keys [^double mu ^double sigma rng]
     :or {mu 1.0 sigma 1.0}}]
   (if (m/< sigma 1.0e-4)
     (poisson {:p mu :rng rng})
     (let [r (m// mu sigma)
           p (m// r (m/+ r mu))]
       (distr/negative-binomial r p rng)))))

(add-distr-method nbii)

(defn logarithmic
  "Creates a logarithmic (log-series) distribution object.

  The logarithmic distribution is a discrete distribution over the positive integers `1, 2, 3, ...`, whose probability mass decays proportionally to `(/ (m/pow p k) k)`. It was introduced by Fisher to model species-abundance data (the number of species represented by a given number of individuals) and is also used to build the negative binomial distribution as a Poisson mixture.

  Parameters (single, optional map):

  - `p` (double): shape parameter, strictly between `0.0` and `1.0`; larger values give a heavier tail. Default: `0.5`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[geometric]], [[negative-binomial]]."
  ([] (logarithmic nil))
  ([{:keys [^double p rng]
     :or {p 0.5}}]
   (distr/logarithmic p rng)))

(add-distr-method logarithmic)

;;

(defn half-cauchy
  "Creates a half-Cauchy distribution object.

  The half-Cauchy distribution is the continuous distribution of the absolute deviation `(- x mu)` of a Cauchy-distributed variable restricted to values not below `mu`, i.e. the right half of a Cauchy distribution centered at `mu` and folded onto `[mu, ##Inf]`. Like the Cauchy distribution, it is heavy-tailed, and it is commonly used as a weakly informative prior for scale parameters in Bayesian modelling.

  Parameters (single, optional map):

  - `mu` (double): location parameter, the lower bound of the support. Default: `0.0`.
  - `scale` (double): scale parameter controlling the spread of the distribution. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Density is `0.0` for `x` below `mu`. Mean and variance are undefined (`##NaN`), since they are undefined for the parent Cauchy distribution.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[cauchy]], [[half-normal]], [[levy]]."
  ([] (half-cauchy nil))
  ([{:keys [^double mu ^double scale rng]
     :or {mu 0.0 scale 1.0}}]
   (distr/half-cauchy mu scale rng)))

(add-distr-method half-cauchy)

;;

(defn integer-discrete-distribution
  "Creates a discrete distribution object over a fixed, finite set of integer values.

  Given a collection of `data` values and matching `probabilities`, the resulting probability mass function assigns to each distinct value the sum of the weights of its occurrences, normalized to sum to `1.0`. It is a general way to build an arbitrary discrete distribution from a user-supplied support and mass, as opposed to distributions parameterized analytically, such as [[poisson]] or [[binomial]].

  Parameters (single, optional map):

  - `data` (sequence of longs): finite integer support of the distribution; repeated values accumulate probability mass. Default: `[1]`.
  - `probabilities` (sequence of doubles): weight associated with each element of `data`, in the same order; weights need not sum to `1.0`, as they are normalized internally. Default: equal weight for every element of `data`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values, a degenerate distribution concentrated on `1`.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[integer-discrete]], [[real-discrete-distribution]], [[categorical-distribution]], [[enumerated-int]]."
  ([] (integer-discrete-distribution nil))
  ([{:keys [data probabilities rng]
     :or {data [1]}}]
   (distr/integer-discrete-distribution data probabilities rng)))

(add-distr-method integer-discrete-distribution)

(defn integer-discrete
  "Creates a discrete distribution object over a fixed, finite set of integer values.

  Alias for [[integer-discrete-distribution]]; see its docstring for the accepted parameters (`:data`, `:probabilities`, `:rng`) and their defaults.

  See also [[distribution]], [[integer-discrete-distribution]]."
  ([] (integer-discrete nil))
  ([{:keys [data probabilities rng]
     :or {data [1]}}]
   (distr/integer-discrete-distribution data probabilities rng)))

(add-distr-method integer-discrete)

(defn real-discrete-distribution
  "Creates a discrete distribution object over a fixed, finite set of real (double) values.

  Given a collection of `data` values and matching `probabilities`, the resulting probability mass function assigns to each distinct value the sum of the weights of its occurrences, normalized to sum to `1.0`. It is the real-valued counterpart of [[integer-discrete-distribution]], useful for building an arbitrary discrete distribution over a fixed, non-integer support.

  Parameters (single, optional map):

  - `data` (sequence of doubles): finite support of the distribution; repeated values accumulate probability mass. Default: `[1.0]`.
  - `probabilities` (sequence of doubles): weight associated with each element of `data`, in the same order; weights need not sum to `1.0`, as they are normalized internally. Default: equal weight for every element of `data`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values, a degenerate distribution concentrated on `1.0`.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[real-discrete]], [[integer-discrete-distribution]], [[categorical-distribution]], [[enumerated-real]]."
  ([] (real-discrete-distribution nil))
  ([{:keys [data probabilities rng]
     :or {data [1.0]}}]
   (distr/real-discrete-distribution data probabilities rng)))

(add-distr-method real-discrete-distribution)

(defn real-discrete
  "Creates a discrete distribution object over a fixed, finite set of real (double) values.

  Alias for [[real-discrete-distribution]]; see its docstring for the accepted parameters (`:data`, `:probabilities`, `:rng`) and their defaults.

  See also [[distribution]], [[real-discrete-distribution]]."
  ([] (real-discrete nil))
  ([{:keys [data probabilities rng]
     :or {data [1.0]}}]
   (distr/real-discrete-distribution data probabilities rng)))

(add-distr-method real-discrete)

(defn categorical-distribution
  "Creates a discrete distribution object over an arbitrary, finite set of category values.

  Given a collection of `data` values of any type, such as keywords, strings or numbers, and matching `probabilities`, the resulting probability mass function assigns to each distinct value the sum of the weights of its occurrences, normalized to sum to `1.0`. 

  Parameters (single, optional map):

  - `data` (sequence of any values): finite support of the distribution, not necessarily numeric; repeated values accumulate probability mass. Default: `[1]`.
  - `probabilities` (sequence of doubles): weight associated with each element of `data`, in the same order; weights need not sum to `1.0`, as they are normalized internally. Default: equal weight for every element of `data`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values, a degenerate distribution concentrated on `1`.

  Mean and variance are undefined (`##NaN`), since the support is not necessarily numeric.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]] and [[sample]].

  See also [[distribution]], [[categorical]], [[integer-discrete-distribution]], [[real-discrete-distribution]]."
  ([] (categorical-distribution nil))
  ([{:keys [data probabilities rng]
     :or {data [1]}}]
   (distr/categorical-distribution data probabilities rng)))

(add-distr-method categorical-distribution)

(defn categorical
  "Creates a discrete distribution object over an arbitrary, finite set of category values.

  Alias for [[categorical-distribution]]; see its docstring for the accepted parameters (`:data`, `:probabilities`, `:rng`) and their defaults.

  See also [[distribution]], [[categorical-distribution]]."
  ([] (categorical nil))
  ([{:keys [data probabilities rng]
     :or {data [1]}}]
   (distr/categorical-distribution data probabilities rng)))

(add-distr-method categorical)

(defn fishers-noncentral-hypergeometric
  "Creates Fisher's noncentral hypergeometric distribution object.

  Fisher's noncentral hypergeometric distribution is a discrete distribution over the number of successes obtained when drawing `n` elements without replacement from a finite population of `ns` success and `nf` failure elements, in which each success element is `omega` times as likely as a failure element to be included among the draws. Equivalently, it is the conditional distribution of one cell of a 2x2 contingency table given both of its margins and an odds ratio `omega`. Setting `omega` to `1.0` recovers the ordinary [[hypergeometric]] distribution. It arises in exact tests and confidence intervals for the odds ratio of a 2x2 table, such as Fisher's exact test.

  Parameters (single, optional map):

  - `ns` (long): number of success elements in the population. Default: `5`.
  - `nf` (long): number of failure elements in the population. Default: `5`.
  - `n` (long): number of elements drawn without replacement. Default: `5`.
  - `omega` (double): odds ratio, the relative likelihood of a success versus a failure element being among the draws; `1.0` corresponds to the ordinary hypergeometric distribution. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  The support is the integer range from `(max 0 (- n nf))` to `(min ns n)`.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[wallenius-noncentral-hypergeometric]], [[hypergeometric]]."
  ([] (fishers-noncentral-hypergeometric nil))
  ([{:keys [^long ns ^long nf ^long n ^double omega rng]
     :or {ns 5 nf 5 n 5 omega 1.0}}]
   (distr/fishers-noncentral-hypergeometric {:ns ns :nf nf :n n :omega omega :rng rng})))

(add-distr-method fishers-noncentral-hypergeometric)

(defn wallenius-noncentral-hypergeometric
  "Creates Wallenius' noncentral hypergeometric distribution object.

  Wallenius' noncentral hypergeometric distribution is a discrete distribution over the number of successes obtained by sequentially drawing `n` elements without replacement from a finite population of `ns` success and `nf` failure elements, where at each individual draw a remaining success element is `omega` times as likely to be picked as a remaining failure element. Unlike [[fishers-noncentral-hypergeometric]], which conditions a pair of counts on their sum, this distribution arises from an explicit biased sequential sampling (urn) process, making it the appropriate model for selection sampling where items are removed one at a time under a persistent bias. Setting `omega` to `1.0` recovers the ordinary [[hypergeometric]] distribution.

  Parameters (single, optional map):

  - `ns` (long): number of success elements in the population. Default: `5`.
  - `nf` (long): number of failure elements in the population. Default: `5`.
  - `n` (long): number of elements drawn without replacement. Default: `5`.
  - `omega` (double): odds ratio, the relative likelihood of a remaining success versus a remaining failure element being picked at each draw; `1.0` corresponds to the ordinary hypergeometric distribution. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  The support is the integer range from `(max 0 (- n nf))` to `(min ns n)`. Unlike [[fishers-noncentral-hypergeometric]], the probability mass function is evaluated via numerical integration, so construction may be relatively slow.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[fishers-noncentral-hypergeometric]], [[hypergeometric]]."
  ([] (wallenius-noncentral-hypergeometric nil))
  ([{:keys [^long ns ^long nf ^long n ^double omega rng]
     :or {ns 5 nf 5 n 5 omega 1.0}}]
   (distr/wallenius-noncentral-hypergeometric {:ns ns :nf nf :n n :omega omega :rng rng})))

(add-distr-method wallenius-noncentral-hypergeometric)

;;

(defn reciprocal
  "Creates a reciprocal (log-uniform) distribution object.

  The reciprocal distribution is a continuous distribution over the positive interval `[a, b]` whose density is proportional to `(/ 1.0 x)`, equivalently, its logarithm, `(m/log x)`, is uniformly distributed over `[(m/log a), (m/log b)]`. It is commonly used to sample scale-invariant quantities spanning several orders of magnitude, such as hyperparameters in a search space or physical quantities of unknown scale, since it assigns equal probability mass to intervals of equal relative (rather than absolute) width.

  Parameters (single, optional map):

  - `a` (double): lower bound of the support, must be strictly positive. Default: `1.0`.
  - `b` (double): upper bound of the support, must be greater than `a`. Default: `10.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[uniform-real]], [[log-normal]], [[pareto]]."
  ([] (reciprocal nil))
  ([{:keys [^double a ^double b rng]
     :or {a 1 b 10}}]
   (distr/reciprocal a b rng)))

(add-distr-method reciprocal)

(defn ex-gaussian
  "Creates an ex-Gaussian (exponentially modified Gaussian) distribution object.

  The ex-Gaussian distribution is the continuous distribution of the sum of an independent [[normal]] variable with mean `mu` and standard deviation `sigma`, and an [[exponential]] variable with mean `tau`. It is right-skewed, with the normal component shaping its bell-like core and the exponential component stretching its right tail, and is widely used to model human reaction-time data as well as chromatography peak shapes.

  Parameters (single, optional map):

  - `mu` (double): mean of the normal component. Default: `0.0`.
  - `sigma` (double): standard deviation of the normal component. Default: `1.0`.
  - `tau` (double): mean of the exponential component, controlling the length of the right tail. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[normal]], [[exponential]], [[exgaus]]."
  ([] (ex-gaussian nil))
  ([{:keys [^double mu ^double sigma ^double tau rng]
     :or {mu 0.0 sigma 1.0 tau 1.0}}]
   (distr/ex-gaussian {:mu mu :sigma sigma :tau tau :rng rng :normal normal :exponential exponential})))

(add-distr-method ex-gaussian)

(defn exgaus
  "Creates an ex-Gaussian (exponentially modified Gaussian) distribution object, using gamlss-style parameter names.

  Alias for [[ex-gaussian]], renaming its `tau` parameter (mean of the exponential component) to `nu`, matching the naming used by R's `gamlss.dist` package. Registered under the `:exgaus` key, distinct from `ex-gaussian`'s `:ex-gaussian` key, so that `(distribution :exgaus {:nu ...})` and `(distribution :ex-gaussian {:tau ...})` are equivalent ways of building the same distribution.

  Parameters (single, optional map):

  - `mu` (double): mean of the normal component. Default: `0.0`.
  - `sigma` (double): standard deviation of the normal component. Default: `1.0`.
  - `nu` (double): mean of the exponential component, controlling the length of the right tail. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[ex-gaussian]]."
  ([] (exgaus nil))
  ([{:keys [^double mu ^double sigma ^double nu rng]
     :or {mu 0.0 sigma 1.0 nu 1.0}}]
   (distr/ex-gaussian {:mu mu :sigma sigma :tau nu :rng rng :normal normal :exponential exponential})))

(add-distr-method exgaus)

(defn beta-binomial
  "Creates a beta-binomial distribution object.

  The beta-binomial distribution is a discrete distribution over the integers `0` to `n`, obtained as a compound of the [[binomial]] distribution whose success probability `p` is itself random, drawn from a [[beta]] distribution with shape parameters `alpha` and `beta`. Compared to the binomial distribution, it is overdispersed, making it a common model for count data (e.g. successes out of `n` trials) exhibiting more variability than the binomial would predict.

  Parameters (single, optional map):

  - `alpha` (double): first shape parameter of the underlying beta distribution. Default: `0.5`.
  - `beta` (double): second shape parameter of the underlying beta distribution. Default: `0.5`.
  - `n` (long): number of trials, the upper bound of the support. Default: `10`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[bb]], [[binomial]], [[beta]], [[hypergeometric]]."
  ([] (beta-binomial nil))
  ([{:keys [^double alpha ^double beta ^long n rng]
     :or {alpha 0.5 beta 0.5 n 10}}]
   (distr/beta-binomial alpha beta n rng)))

(add-distr-method beta-binomial)

(defn bb
  "Creates a beta-binomial distribution object, using gamlss-style parameter names.

  Alias for [[beta-binomial]], reparameterized in terms of the mean success probability `mu` and a dispersion parameter `sigma`, matching the naming used by R's `gamlss.dist` package (its `BB` family). The underlying [[beta-binomial]] shape parameters are recovered as `alpha = mu / sigma` and `beta = (1 - mu) / sigma`.

  Parameters (single, optional map):

  - `mu` (double): mean success probability, in `(0, 1)`. Default: `0.5`.
  - `sigma` (double): dispersion parameter; larger values give more overdispersion relative to the [[binomial]] distribution. Default: `1.0`.
  - `bd` (long): binomial denominator, i.e. number of trials, the upper bound of the support. Default: `10`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[beta-binomial]]."
  ([] (bb nil))
  ([{:keys [^double mu ^double sigma ^long bd rng]
     :or {mu 0.5 sigma 1.0 bd 10}}]
   (distr/beta-binomial (m// mu sigma) (m// (m/- 1.0 mu) sigma) bd rng)))

(add-distr-method bb)

(defn zero-inflated-binomial
  "Creates a zero-inflated binomial distribution object, using gamlss-style parameter names.

  The zero-inflated binomial distribution augments the [[binomial]] distribution with an extra point mass at `0`, on top of whatever probability the binomial itself already places there. It is a common model for count data (successes out of `bd` trials) exhibiting more zeros than a plain [[binomial]] would predict. With probability `sigma` the outcome is forced to `0`; with probability `(1 - sigma)` it is drawn from a binomial distribution with success probability `mu` and `bd` trials. Matches the naming and parameterization used by R's `gamlss.dist` package (its `ZIBI` family).

  Parameters (single, optional map):

  - `mu` (double): success probability of the underlying binomial component, in `(0, 1)`. Default: `0.5`.
  - `sigma` (double): zero-inflation probability, in `[0, 1)`, i.e. the extra probability of observing `0` beyond the binomial's own mass at `0`. Default: `0.1`.
  - `bd` (long): binomial denominator, i.e. number of trials, the upper bound of the support. Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[binomial]], [[zero-inflated-beta-binomial]]."
  ([] (zero-inflated-binomial nil))
  ([{:keys [^double mu ^double sigma ^long bd rng]
     :or {mu 0.5 sigma 0.1 bd 1}}]
   (distr/zero-inflated-binomial {:mu mu :sigma sigma :bd bd :rng rng :binomial binomial})))

(add-distr-method zero-inflated-binomial)
(add-distr-method zero-inflated-binomial :zibi)

(defn zero-adjusted-binomial
  "Creates a zero-adjusted binomial distribution object, using gamlss-style parameter names.

  The zero-adjusted binomial distribution is a hurdle model over the integers `0` to `bd`: with probability `sigma` the outcome is exactly `0`, and with probability `(1 - sigma)` it is drawn from a zero-truncated binomial distribution with success probability `mu` and `bd` trials, i.e. the binomial's own positive-value probabilities rescaled to sum to `(1 - sigma)`. Unlike [[zero-inflated-binomial]], which adds an extra point mass at `0` on top of the binomial's own mass there, this distribution replaces that mass entirely, so `sigma` is exactly the probability of observing `0`, not merely an addition to it. Matches the naming and parameterization used by R's `gamlss.dist` package (its `ZABI` family).

  Parameters (single, optional map):

  - `mu` (double): success probability of the underlying binomial component, in `(0, 1)`. Default: `0.5`.
  - `sigma` (double): probability of observing `0`, in `[0, 1)`. Default: `0.1`.
  - `bd` (long): binomial denominator, i.e. number of trials, the upper bound of the support. Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[binomial]], [[zero-inflated-binomial]]."
  ([] (zero-adjusted-binomial nil))
  ([{:keys [^double mu ^double sigma ^long bd rng]
     :or {mu 0.5 sigma 0.1 bd 1}}]
   (distr/zero-adjusted-binomial {:mu mu :sigma sigma :bd bd :rng rng :binomial binomial})))

(add-distr-method zero-adjusted-binomial)
(add-distr-method zero-adjusted-binomial :zabi)


(defn zero-inflated-beta-binomial
  "Creates a zero-inflated beta-binomial distribution object, using gamlss-style parameter names.

  The zero-inflated beta-binomial distribution augments the [[bb]]/[[beta-binomial]] distribution with an extra point mass at `0`, on top of whatever probability the beta-binomial itself already places there. It is a common model for count data (successes out of `bd` trials) exhibiting more zeros than a plain beta-binomial or [[binomial]] would predict. With probability `nu` the outcome is forced to `0`; with probability `(1 - nu)` it is drawn from a beta-binomial distribution parameterized by mean success probability `mu`, dispersion `sigma` and `bd` trials (as in [[bb]]). Matches the naming and parameterization used by R's `gamlss.dist` package (its `ZIBB` family).

  Parameters (single, optional map):

  - `mu` (double): mean success probability of the underlying beta-binomial component, in `(0, 1)`. Default: `0.5`.
  - `sigma` (double): dispersion parameter of the underlying beta-binomial component; larger values give more overdispersion relative to the [[binomial]] distribution. Default: `0.1`.
  - `nu` (double): zero-inflation probability, in `[0, 1)`, i.e. the extra probability of observing `0` beyond the beta-binomial's own mass at `0`. Default: `0.1`.
  - `bd` (long): binomial denominator, i.e. number of trials, the upper bound of the support. Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[bb]], [[beta-binomial]]."
  ([] (zero-inflated-beta-binomial nil))
  ([{:keys [^double mu ^double sigma ^long bd ^double nu rng]
     :or {mu 0.5 sigma 0.1 nu 0.1 bd 1}}]
   (distr/zero-inflated-beta-binomial {:mu mu :sigma sigma :nu nu :bd bd :rng rng})))

(add-distr-method zero-inflated-beta-binomial)
(add-distr-method zero-inflated-beta-binomial :zibb)

(defn zero-adjusted-beta-binomial
  "Creates a zero-adjusted beta-binomial distribution object, using gamlss-style parameter names.

  The zero-adjusted beta-binomial distribution is a hurdle model over the integers `0` to `bd`: with probability `nu` the outcome is exactly `0`, and with probability `(1 - nu)` it is drawn from a zero-truncated beta-binomial distribution parameterized by mean success probability `mu`, dispersion `sigma` and `bd` trials (as in [[bb]]), i.e. the beta-binomial's own positive-value probabilities rescaled to sum to `(1 - nu)`. Unlike [[zero-inflated-beta-binomial]], which adds an extra point mass at `0` on top of the beta-binomial's own mass there, this distribution replaces that mass entirely, so `nu` is exactly the probability of observing `0`, not merely an addition to it. Matches the naming and parameterization used by R's `gamlss.dist` package (its `ZABB` family).

  Parameters (single, optional map):

  - `mu` (double): mean success probability of the underlying beta-binomial component, in `(0, 1)`. Default: `0.5`.
  - `sigma` (double): dispersion parameter of the underlying beta-binomial component; larger values give more overdispersion relative to the [[binomial]] distribution. Default: `0.1`.
  - `nu` (double): probability of observing `0`, in `[0, 1)`. Default: `0.1`.
  - `bd` (long): binomial denominator, i.e. number of trials, the upper bound of the support. Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[bb]], [[beta-binomial]], [[zero-inflated-beta-binomial]]."
  ([] (zero-adjusted-beta-binomial nil))
  ([{:keys [^double mu ^double sigma ^long bd ^double nu rng]
     :or {mu 0.5 sigma 0.1 nu 0.1 bd 1}}]
   (distr/zero-adjusted-beta-binomial {:mu mu :sigma sigma :nu nu :bd bd :rng rng})))

(add-distr-method zero-adjusted-beta-binomial)
(add-distr-method zero-adjusted-beta-binomial :zabb)

(defn zero-inflated-negative-binomial
  "Creates a zero-inflated negative binomial distribution object, using gamlss-style parameter names.

  The zero-inflated negative binomial distribution augments the [[nbi]]/[[negative-binomial]] distribution with an extra point mass at `0`, on top of whatever probability the negative binomial itself already places there. It is a common model for unbounded count data exhibiting more zeros than a plain [[nbi]] or [[poisson]] would predict. With probability `nu` the outcome is forced to `0`; with probability `(1 - nu)` it is drawn from a negative binomial distribution parameterized by mean `mu` and dispersion `sigma` (as in [[nbi]]). Matches the naming and parameterization used by R's `gamlss.dist` package (its `ZINBI` family).

  Parameters (single, optional map):

  - `mu` (double): mean of the underlying negative binomial component, strictly positive. Default: `1.0`.
  - `sigma` (double): dispersion parameter of the underlying negative binomial component; larger values give more overdispersion relative to a [[poisson]] distribution with the same mean. Default: `1.0`.
  - `nu` (double): zero-inflation probability, in `[0, 1)`, i.e. the extra probability of observing `0` beyond the negative binomial's own mass at `0`. Default: `0.3`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[nbi]], [[negative-binomial]], [[zero-inflated-binomial]]."
  ([] (zero-inflated-negative-binomial nil))
  ([{:keys [^double mu ^double sigma ^double nu rng]
     :or {mu 1.0 sigma 1.0 nu 0.3}}]
   (distr/zero-inflated-negative-binomial {:mu mu :sigma sigma :nu nu :rng rng :nbi nbi})))

(add-distr-method zero-inflated-negative-binomial)
(add-distr-method zero-inflated-negative-binomial :zinbi)

(defn zero-adjusted-negative-binomial
  "Creates a zero-adjusted negative binomial distribution object, using gamlss-style parameter names.

  The zero-adjusted negative binomial distribution is a hurdle model over the non-negative integers: with probability `nu` the outcome is exactly `0`, and with probability `(1 - nu)` it is drawn from a zero-truncated negative binomial distribution parameterized by mean `mu` and dispersion `sigma` (as in [[nbi]]), i.e. the negative binomial's own positive-value probabilities rescaled to sum to `(1 - nu)`. Unlike [[zero-inflated-negative-binomial]], which adds an extra point mass at `0` on top of the negative binomial's own mass there, this distribution replaces that mass entirely, so `nu` is exactly the probability of observing `0`, not merely an addition to it. Matches the naming and parameterization used by R's `gamlss.dist` package (its `ZANBI` family).

  Parameters (single, optional map):

  - `mu` (double): mean of the underlying negative binomial component, strictly positive. Default: `1.0`.
  - `sigma` (double): dispersion parameter of the underlying negative binomial component; larger values give more overdispersion relative to a [[poisson]] distribution with the same mean. Default: `1.0`.
  - `nu` (double): probability of observing `0`, in `[0, 1)`. Default: `0.3`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[nbi]], [[negative-binomial]], [[zero-inflated-negative-binomial]]."
  ([] (zero-adjusted-negative-binomial nil))
  ([{:keys [^double mu ^double sigma ^double nu rng]
     :or {mu 1.0 sigma 1.0 nu 0.3}}]
   (distr/zero-adjusted-negative-binomial {:mu mu :sigma sigma :nu nu :rng rng :nbi nbi})))

(add-distr-method zero-adjusted-negative-binomial)
(add-distr-method zero-adjusted-negative-binomial :zanbi)

(defn zero-inflated-poisson
  "Creates a zero-inflated Poisson distribution object, using gamlss-style parameter names.

  The zero-inflated Poisson distribution augments the [[poisson]] distribution with an extra point mass at `0`, on top of whatever probability the Poisson itself already places there. It is a common model for count data exhibiting more zeros than a plain [[poisson]] would predict. With probability `sigma` the outcome is forced to `0`; with probability `(1 - sigma)` it is drawn from a [[poisson]] distribution with rate `mu`. Note that the mean of the resulting distribution is `(1 - sigma) * mu`, not `mu` itself; see [[zero-inflated-poisson2]] for a reparameterization where `mu` is the mean directly. Matches the naming and parameterization used by R's `gamlss.dist` package (its `ZIP` family).

  Parameters (single, optional map):

  - `mu` (double): rate of the underlying [[poisson]] component, strictly positive. Default: `5.0`.
  - `sigma` (double): zero-inflation probability, in `[0, 1)`, i.e. the extra probability of observing `0` beyond the Poisson's own mass at `0`. Default: `0.1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[poisson]], [[zero-inflated-poisson2]], [[zero-inflated-negative-binomial]]."
  ([] (zero-inflated-poisson nil))
  ([{:keys [^double mu ^double sigma rng]
     :or {mu 5.0 sigma 0.1}}]
   (distr/zero-inflated-poisson {:mu mu :sigma sigma :rng rng :poisson poisson})))

(add-distr-method zero-inflated-poisson)
(add-distr-method zero-inflated-poisson :zip)

(defn zero-inflated-poisson2
  "Creates a zero-inflated Poisson distribution object, mean-parameterized, using gamlss-style parameter names.

  Alias for [[zero-inflated-poisson]], reparameterized so that `mu` is the mean of the resulting distribution directly (rather than the rate of the underlying [[poisson]] component). The underlying [[poisson]] rate is recovered as `mu / (1 - sigma)`, so that the mean works out to `mu` exactly, matching the naming and parameterization used by R's `gamlss.dist` package (its `ZIP2` family).

  Parameters (single, optional map):

  - `mu` (double): mean of the resulting distribution, strictly positive. Default: `5.0`.
  - `sigma` (double): zero-inflation probability, in `[0, 1)`. Default: `0.1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[poisson]], [[zero-inflated-poisson]]."
  ([] (zero-inflated-poisson2 nil))
  ([{:keys [^double mu ^double sigma rng]
     :or {mu 5.0 sigma 0.1}}]
   (let [nmu (m// mu (m/- 1.0 sigma))]
     (distr/zero-inflated-poisson {:mu nmu :sigma sigma :rng rng :poisson poisson}))))

(add-distr-method zero-inflated-poisson2)
(add-distr-method zero-inflated-poisson2 :zip2)

(defn zero-adjusted-poisson
  "Creates a zero-adjusted Poisson distribution object, using gamlss-style parameter names.

  The zero-adjusted Poisson distribution is a hurdle model over the non-negative integers: with probability `sigma` the outcome is exactly `0`, and with probability `(1 - sigma)` it is drawn from a zero-truncated [[poisson]] distribution with rate `mu`, i.e. the Poisson's own positive-value probabilities rescaled to sum to `(1 - sigma)`. Unlike [[zero-inflated-poisson]], which adds an extra point mass at `0` on top of the Poisson's own mass there, this distribution replaces that mass entirely, so `sigma` is exactly the probability of observing `0`, not merely an addition to it. Matches the naming and parameterization used by R's `gamlss.dist` package (its `ZAP` family).

  Parameters (single, optional map):

  - `mu` (double): rate of the underlying [[poisson]] component, strictly positive. Default: `5.0`.
  - `sigma` (double): probability of observing `0`, in `[0, 1)`. Default: `0.1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[poisson]], [[zero-inflated-poisson]]."
  ([] (zero-adjusted-poisson nil))
  ([{:keys [^double mu ^double sigma rng]
     :or {mu 5.0 sigma 0.1}}]
   (distr/zero-adjusted-poisson {:mu mu :sigma sigma :rng rng :poisson poisson})))

(add-distr-method zero-adjusted-poisson)
(add-distr-method zero-adjusted-poisson :zap)

(defn zero-adjusted-gamma
  "Creates a zero-adjusted gamma distribution object, using gamlss-style parameter names.

  The zero-adjusted gamma distribution is a hurdle model mixing a discrete point mass at `0` with an otherwise continuous [[gamma]] distribution: with probability `nu` the outcome is exactly `0`, and with probability `(1 - nu)` it is drawn from a [[gamma]] distribution parameterized by mean `mu` and coefficient of variation `sigma` (shape `1/sigma^2`, scale `sigma^2 * mu`). Since the underlying [[gamma]] distribution places no mass at `0` itself, `nu` is exactly `P(X = 0)`, with no rescaling needed for the positive part (unlike the discrete hurdle distributions such as [[zero-adjusted-negative-binomial]]). Matches the naming and parameterization used by R's `gamlss.dist` package (its `ZAGA` family).

  Parameters (single, optional map):

  - `mu` (double): mean of the underlying gamma component, strictly positive. Default: `1.0`.
  - `sigma` (double): coefficient of variation of the underlying gamma component, strictly positive. Default: `1.0`.
  - `nu` (double): probability of observing `0`, in `[0, 1)`. Default: `0.1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[gamma]], [[zero-adjusted-negative-binomial]]."
  ([] (zero-adjusted-gamma nil))
  ([{:keys [^double mu ^double sigma ^double nu rng]
     :or {mu 1.0 sigma 1.0 nu 0.1}}]
   (distr/zero-adjusted-gamma {:mu mu :sigma sigma :nu nu :rng rng :gamma gamma})))

(add-distr-method zero-adjusted-gamma)
(add-distr-method zero-adjusted-gamma :zaga)

(defn zero-adjusted-inverse-gaussian
  "Creates a zero-adjusted inverse Gaussian distribution object, using gamlss-style parameter names.

  The zero-adjusted inverse Gaussian distribution is a hurdle model mixing a discrete point mass at `0` with an otherwise continuous [[inverse-gaussian]] distribution: with probability `nu` the outcome is exactly `0`, and with probability `(1 - nu)` it is drawn from an [[inverse-gaussian]] distribution parameterized by mean `mu` and dispersion `sigma` (`lambda = 1/sigma^2`). Since the underlying [[inverse-gaussian]] distribution places no mass at `0` itself, `nu` is exactly `P(X = 0)`, with no rescaling needed for the positive part (unlike the discrete hurdle distributions such as [[zero-adjusted-negative-binomial]]). Matches the naming and parameterization used by R's `gamlss.dist` package (its `ZAIG` family).

  Parameters (single, optional map):

  - `mu` (double): mean of the underlying inverse Gaussian component, strictly positive. Default: `1.0`.
  - `sigma` (double): dispersion of the underlying inverse Gaussian component, strictly positive. Default: `1.0`.
  - `nu` (double): probability of observing `0`, in `[0, 1)`. Default: `0.1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[inverse-gaussian]], [[zero-adjusted-gamma]]."
  ([] (zero-adjusted-inverse-gaussian nil))
  ([{:keys [^double mu ^double sigma ^double nu rng]
     :or {mu 1.0 sigma 1.0 nu 0.1}}]
   (distr/zero-adjusted-inverse-gaussian {:mu mu :sigma sigma :nu nu :rng rng :inverse-gaussian inverse-gaussian})))

(add-distr-method zero-adjusted-inverse-gaussian)
(add-distr-method zero-adjusted-inverse-gaussian :zaig)

(defn generalized-extreme-value
  "Creates a generalized extreme value (GEV) distribution object.

  The GEV distribution is the limiting distribution of normalized maxima (or, via negation, minima) of independent, identically distributed samples, and unifies the three classical extreme-value families: Gumbel (`xi = 0`), Fréchet (`xi > 0`, heavy right tail), and (reversed) Weibull (`xi < 0`, bounded above). It is standard in extreme-value analysis of e.g. flood, temperature, or wind-speed maxima.

  For shape `xi != 0` and `z = 1 + xi * (x - mu) / sigma`, the CDF is `exp(-z^(-1/xi))` on the support where `z > 0`; for `xi = 0` it degenerates to the Gumbel CDF `exp(-exp(-(x - mu) / sigma))`.

  Parameters (single, optional map):

  - `mu` (double): location parameter. Default: `0.0`.
  - `sigma` (double): scale parameter, strictly positive. Default: `1.0`.
  - `xi` (double): shape parameter. `xi = 0` gives the Gumbel (type I) distribution, unbounded on both sides; `xi > 0` gives the Fréchet (type II) distribution, bounded below at `mu - sigma/xi`, with a heavy right tail; `xi < 0` gives the (reversed) Weibull (type III) distribution, bounded above at `mu - sigma/xi`. Default: `0.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Note: `mean` and `variance` are finite only for `xi < 1` and `xi < 0.5` respectively; outside those ranges they are `##Inf`.

  Matches the standard `(mu, sigma, xi)` parameterization used in extreme-value statistics; note that R's `EnvStats` package (`GEVD` family) uses a shape parameter `kappa` such that `xi = -kappa`.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]]."
  ([] (generalized-extreme-value nil))
  ([{:keys [^double mu ^double sigma ^double xi rng]
     :or {mu 0.0 sigma 1.0 xi 0.0}}]
   (distr/generalized-extreme-value mu sigma xi rng)))

(add-distr-method generalized-extreme-value)
(add-distr-method generalized-extreme-value :gev)

(defn generalized-logistic
  "Creates a generalized logistic distribution object (type I / skew-logistic).

  The type I generalized logistic distribution is a skewed generalization of the standard logistic distribution, with CDF equal to the standard logistic CDF raised to the power `alpha`: `F(x) = (1 + exp(-z))^(-alpha)`, where `z = (x - mu) / sigma`. Equivalently, it is the distribution of the maximum of `alpha` i.i.d. standard logistic random variables (with `alpha` generalized to any positive real). `alpha = 1` recovers the standard (symmetric) logistic distribution; `alpha != 1` introduces skew, right-skewed for `alpha > 1` and left-skewed for `alpha < 1`.

  Parameters (single, optional map):

  - `mu` (double): location parameter. Default: `0.0`.
  - `sigma` (double): scale parameter, strictly positive. Default: `1.0`.
  - `alpha` (double): shape parameter, strictly positive. `alpha = 1` gives the standard logistic distribution. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Support is the whole real line for any `alpha > 0`; mean and variance are always finite, given by `mu + sigma * (digamma(alpha) + Euler-Mascheroni constant)` and `sigma^2 * (trigamma(alpha) + pi^2/6)` respectively.

  Matches the `(location, scale, shape)` parameterization used by R's `glogis` package.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[logistic]]."
  ([] (generalized-logistic nil))
  ([{:keys [^double mu ^double sigma ^double alpha rng]
     :or {mu 0.0 sigma 1.0 alpha 1.0}}]
   (distr/generalized-logistic mu sigma alpha rng)))

(add-distr-method generalized-logistic)

(defn generalized-pareto
  "Creates a generalized Pareto distribution (GPD) object.

  The GPD is the limiting distribution of excesses over a high threshold (i.e. of `X - u` given `X > u`, as `u` increases), and is the standard model for peaks-over-threshold extreme value analysis, complementing the [[generalized-extreme-value]] distribution's block-maxima approach. It unifies the exponential (`xi = 0`), Pareto type II / Lomax (`xi > 0`, heavy right tail), and bounded (`xi < 0`) families.

  For shape `xi != 0` and `t = 1 + xi * (x - mu) / sigma`, the CDF is `1 - t^(-1/xi)` on the support `x >= mu` where `t > 0`; for `xi = 0` it degenerates to the shifted exponential CDF `1 - exp(-(x - mu) / sigma)`.

  Parameters (single, optional map):

  - `mu` (double): location parameter; also the (always finite) lower bound of the support. Default: `0.0`.
  - `sigma` (double): scale parameter, strictly positive. Default: `1.0`.
  - `xi` (double): shape parameter. `xi = 0` gives the shifted exponential distribution, unbounded above; `xi > 0` gives a Pareto type II (Lomax) distribution, unbounded above with a heavy right tail; `xi < 0` gives a distribution bounded above at `mu - sigma/xi`. Default: `0.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Note: `mean` and `variance` are finite only for `xi < 1` and `xi < 0.5` respectively; outside those ranges they are `##Inf`.

  Matches the standard `(mu, sigma, xi)` / `(loc, scale, shape)` parameterization used in extreme-value statistics, e.g. R's `evd` package (`gpd` family).

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[generalized-extreme-value]], [[pareto]]."
  ([] (generalized-pareto nil))
  ([{:keys [^double mu ^double sigma ^double xi rng]
     :or {mu 0.0 sigma 1.0 xi 0.0}}]
   (distr/generalized-pareto mu sigma xi rng)))

(add-distr-method generalized-pareto)
(add-distr-method generalized-pareto :gpd)

(defn generalized-exponential
  "Creates a generalized exponential distribution object (Gupta-Kundu exponentiated exponential distribution).

  The generalized exponential distribution raises the standard exponential CDF to a shape power `alpha`: `F(x) = (1 - exp(-lambda * x))^alpha`, for `x >= 0`. It was introduced by Gupta and Kundu (1999) as a flexible alternative to the gamma and Weibull distributions, sharing the same two-parameter (shape, scale/rate) structure. `alpha = 1` recovers the standard exponential distribution with rate `lambda`.

  Parameters (single, optional map):

  - `alpha` (double): shape parameter, strictly positive. For `alpha < 1` the pdf is strictly decreasing (with an unbounded value at `x = 0`); for `alpha > 1` it is unimodal, rising from `0` at `x = 0` to a peak and then decaying; `alpha = 1` gives the (monotonically decreasing) exponential density. Default: `1.0`.
  - `lambda` (double): rate (scale) parameter, strictly positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Mean and variance are always finite, given by `(digamma(alpha + 1) + Euler-Mascheroni constant) / lambda` and `(pi^2/6 - trigamma(alpha + 1)) / lambda^2` respectively (Gupta & Kundu 1999).

  Matches the `(alpha, lambda)` (shape, scale) parameterization used by R's `reliaR` package (`gen.exp` family).

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[exponential]], [[weibull]], [[gamma]]."
  ([] (generalized-exponential nil))
  ([{:keys [^double alpha ^double lambda rng]
     :or {alpha 1.0 lambda 1.0}}]
   (distr/generalized-exponential alpha lambda rng)))

(add-distr-method generalized-exponential)
(add-distr-method generalized-exponential :ge)

(defn generalized-gamma
  "Creates a generalized gamma distribution object (Stacy distribution), using gamlss-style parameter names.

  The generalized gamma distribution is a flexible three-parameter continuous distribution over positive reals that includes the [[gamma]] (`nu = 1`), Weibull (`sigma = 1`), and log-normal (`nu = 0`, in the limit) distributions as special cases. Following the parameterization used by R's `gamlss.dist` package (its `GG` family, Lopatatzidis & Green 2000), its density is `f(y) = theta^theta * z^theta * |nu| * exp(-theta * z) / (Gamma(theta) * y)`, where `z = (y/mu)^nu` and `theta = 1 / (sigma^2 * nu^2)`, for `y > 0`.

  Equivalently, `Y = mu * (U / theta)^(1/nu)` where `U` follows a [[gamma]] distribution with shape `theta` and scale `1`; this relationship is used internally for `cdf`/`icdf`/sampling. For `nu = 0`, the distribution degenerates to [[log-normal]] with `scale = log(mu)`, `shape = sigma`.

  Parameters (single, optional map):

  - `mu` (double): scale parameter, strictly positive. For `nu = 1` (i.e. the ordinary [[gamma]] distribution) this is exactly the mean. Default: `1.0`.
  - `sigma` (double): dispersion parameter, strictly positive; smaller values concentrate the distribution more tightly. Default: `0.5`.
  - `nu` (double): shape parameter, any real number; controls the skewness/tail behavior and which special case the distribution reduces to (`nu = 1`: gamma; `nu = 0`: log-normal; general `nu`: full Stacy generalized gamma family). Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Note: `mean` is finite only for `nu > -1/sigma^2`, and `variance` only for `nu > -1/(2*sigma^2)`; outside those ranges they are `##Inf`.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[gamma]], [[weibull]], [[log-normal]]."
  ([] (generalized-gamma nil))
  ([{:keys [^double mu ^double sigma ^double nu rng]
     :or {mu 1.0 sigma 0.5 nu 1.0}}]
   (distr/generalized-gamma {:mu mu :sigma sigma :nu nu :rng rng :gamma gamma :log-normal log-normal})))

(add-distr-method generalized-gamma)
(add-distr-method generalized-gamma :gg)

(defn generalized-normal
  "Creates a generalized normal distribution object (exponential power / Subbotin distribution).

  The generalized normal distribution is a symmetric, unimodal continuous distribution over the whole real line with density `f(x) = beta / (2 * alpha * Gamma(1/beta)) * exp(-(|x - mu| / alpha)^beta)`. It generalizes several common distributions via its shape parameter `beta`: `beta = 2` gives the [[normal]] distribution (with standard deviation `alpha / sqrt(2)`), `beta = 1` gives the [[laplace]] distribution (with scale `alpha`), and `beta -> Infinity` approaches a uniform distribution on `[mu - alpha, mu + alpha]`. Larger `beta` produces flatter, more platykurtic shapes; smaller `beta` produces more peaked, heavier-tailed shapes.

  Parameters (single, optional map):

  - `mu` (double): location parameter (and, by symmetry, the mean). Default: `0.0`.
  - `alpha` (double): scale parameter, strictly positive. Default: `1.0`.
  - `beta` (double): shape parameter, strictly positive. Default: `1.0` (i.e. [[laplace]] by default).
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  `variance` is `alpha^2 * Gamma(3/beta) / Gamma(1/beta)`, always finite for `beta > 0`.

  Matches the `(mu, alpha, beta)` parameterization used by R's `gnorm` package.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[normal]], [[laplace]]."
  ([] (generalized-normal nil))
  ([{:keys [^double mu ^double alpha ^double beta rng]
     :or {mu 0.0 alpha 1.0 beta 1.0}}]
   (distr/generalized-normal {:mu mu :alpha alpha :beta beta :rng rng :gamma gamma})))

(add-distr-method generalized-normal)
(add-distr-method generalized-normal :gnd)

(defn generalized-inverse-gaussian
  "Creates a generalized inverse Gaussian distribution object (GIG, Sichel's distribution).

  The generalized inverse Gaussian distribution is a flexible three-parameter continuous distribution over positive reals whose density involves the modified Bessel function of the second kind, `K_lambda`: `f(x) = (psi/chi)^(lambda/2) / (2 * K_lambda(sqrt(chi*psi))) * x^(lambda - 1) * exp(-(chi/x + psi*x) / 2)`, for `x > 0`. It includes the [[gamma]] distribution as a limiting case (`chi -> 0`, `lambda > 0`), the inverse gamma distribution as a limiting case (`psi -> 0`, `lambda < 0`), and the ordinary inverse Gaussian distribution as a special case (`lambda = -1/2`). It is widely used in finance (as the mixing distribution of the generalized hyperbolic distribution family) and in actuarial/ecological modelling (Sichel's distribution).

  Parameters (single, optional map):

  - `chi` (double): strictly positive; controls the behavior near `0`.
  - `psi` (double): strictly positive; controls the tail decay.
  - `lambda` (double): shape parameter, any real number.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Note: this implementation requires `chi > 0` and `psi > 0` strictly (the general, non-degenerate GIG); the gamma/inverse-gamma limiting cases (`chi = 0` or `psi = 0`) are not supported directly — use [[gamma]] or an appropriately-parameterized [[gamma]]-based distribution instead.

  `mean` and `variance` are always finite, computed from ratios of modified Bessel functions of the second kind: `mean = sqrt(chi/psi) * K_(lambda+1)(omega) / K_lambda(omega)` where `omega = sqrt(chi*psi)`, and similarly for the second moment via `K_(lambda+2)`. There is no closed-form `cdf`/`icdf`, so at construction time the density is numerically integrated once over a generous range (25 standard deviations past the mean) into a monotone-interpolated table (following the same approach as [[continuous-distribution]]), which `cdf` and `icdf` then look up; this makes construction itself relatively more expensive (single-digit milliseconds), but individual `cdf`/`icdf` calls fast, at the cost of a small interpolation error (empirically below `1e-5` in absolute terms across a range of tested parameters).

  Matches the `(chi, psi, lambda)` parameterization used by R's `GeneralizedHyperbolic` package (`dgig`/`pgig`/`qgig`).

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[gamma]], [[inverse-gaussian]]."
  ([] (generalized-inverse-gaussian nil))
  ([{:keys [^double chi ^double psi ^double lambda rng]
     :or {chi 1.0 psi 1.0 lambda 1.0}}]
   (distr/generalized-inverse-gaussian {:chi chi :psi psi :lambda lambda :rng rng})))

(add-distr-method generalized-inverse-gaussian)
(add-distr-method generalized-inverse-gaussian :gig)

(defn generalized-hyperbolic
  "Creates a generalized hyperbolic distribution object (GH).

  The generalized hyperbolic distribution is a flexible five-parameter continuous distribution over the whole real line that can model both skewness and (leptokurtic or platykurtic) excess kurtosis. It arises as a normal variance-mean mixture: `X = mu + beta*W + sqrt(W)*Z`, where `Z` is standard normal and `W` follows a [[generalized-inverse-gaussian]] distribution with `chi = delta^2`, `psi = alpha^2 - beta^2`, and the same `lambda`. Its density is `f(x) = c * K_(lambda-1/2)(alpha*s) * (s/alpha)^(lambda-1/2) * exp(beta*(x-mu))`, where `s = sqrt(delta^2 + (x-mu)^2)`, `K_v` is the modified Bessel function of the second kind, and `c` is a normalizing constant depending on `alpha`, `beta`, `delta` and `lambda`.

  It includes several well-known distributions as special or limiting cases: the [[normal-inverse-gaussian]] distribution (`lambda = -1/2`), the hyperbolic distribution (`lambda = 1`), the variance-gamma distribution (`delta = 0`, `lambda > 0`, not supported directly by this implementation — see the note below), and approaches the [[normal]] distribution as `delta -> Infinity` with `delta/alpha^2` held fixed. It is widely used in finance to model asset returns.

  Parameters (single, optional map):

  - `mu` (double): location parameter. Default: `0.0`.
  - `delta` (double): scale parameter, strictly positive. Default: `1.0`.
  - `alpha` (double): tail heaviness parameter, strictly positive. Default: `1.0`.
  - `beta` (double): asymmetry/skewness parameter; must satisfy `(< (m/abs beta) alpha)`. Default: `0.0`.
  - `lambda` (double): shape parameter, any real number; controls which named sub-family the distribution resembles. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Note: this implementation requires `delta > 0` and `(< (m/abs beta) alpha)` strictly (the general, non-degenerate GH, mirroring the same restriction as its [[generalized-inverse-gaussian]] mixing distribution); boundary/degenerate cases such as the variance-gamma distribution (`delta = 0`) or the skewed Student's t-like boundary case (`(== (m/abs beta) alpha)`, `lambda < 0`) are not supported directly.

  `mean` and `variance` are always finite, computed from the mixture representation as `mean = mu + beta * E[W]` and `variance = E[W] + beta^2 * Var[W]`, where `E[W]`/`Var[W]` are the mean/variance of the mixing [[generalized-inverse-gaussian]] distribution (themselves ratios of modified Bessel functions of the second kind). There is no closed-form `cdf`/`icdf`, so — following the same approach as [[generalized-inverse-gaussian]] — the density is numerically integrated once at construction time (over a generous range of 25 standard deviations either side of the mean) into a monotone-interpolated table which `cdf` and `icdf` then look up. `sample` bypasses this table and instead draws directly from the mixture representation (sampling `W` from the mixing [[generalized-inverse-gaussian]] distribution and `Z` from a standard normal), which is both exact and independent of the table's interpolation error.

  Matches the `(mu, delta, alpha, beta, lambda)` parameterization used by R's `GeneralizedHyperbolic` package (`dghyp`/`pghyp`/`qghyp`).

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[generalized-inverse-gaussian]], [[normal-inverse-gaussian]], [[normal]]."
  ([] (generalized-hyperbolic nil))
  ([{:keys [^double mu ^double delta ^double alpha ^double beta ^double lambda rng]
     :or {mu 0.0 delta 1.0 alpha 1.0 beta 0.0 lambda 1.0}}]
   (distr/generalized-hyperbolic {:mu mu :delta delta :alpha alpha :beta beta :lambda lambda :rng rng})))

(add-distr-method generalized-hyperbolic)
(add-distr-method generalized-hyperbolic :gh)

(defn half-logistic
  "Creates a half-logistic distribution object.

  The half-logistic distribution is the distribution of `|X|` for `X` a [[logistic]] random variable scaled by `scale`: a continuous distribution over non-negative reals with density `f(x) = (2/scale) * e^(-x/scale) / (1+e^(-x/scale))^2`, cdf `F(x) = (1-e^(-x/scale))/(1+e^(-x/scale))`, and quantile function `icdf(p) = 2*scale*atanh(p)`, all closed-form, for `x >= 0`.

  It is exactly the `alpha = 1` special case of [[generalized-half-logistic]] (with `lambda = 1/scale`); unlike the generalized version, `mean` and `variance` here have simple closed forms: `mean = 2*ln(2)*scale` and `variance = (pi^2/3 - 4*ln(2)^2) * scale^2`.

  Parameters (single, optional map):

  - `scale` (double): scale parameter, strictly positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Matches the `scale` parameterization used by R's `bayesmeta` package (`dhalflogistic`/`phalflogistic`/`qhalflogistic`).

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[generalized-half-logistic]], [[logistic]], [[half-normal]], [[half-cauchy]]."
  ([] (half-logistic nil))
  ([{:keys [^double scale rng]
     :or {scale 1.0}}]
   (distr/half-logistic scale rng)))

(add-distr-method half-logistic)

(defn generalized-half-logistic
  "Creates a generalized half-logistic distribution object (exponentiated half-logistic distribution).

  The generalized half-logistic distribution is a two-parameter continuous distribution over non-negative reals, obtained by exponentiating the cdf of the half-logistic distribution (the distribution of `|Z|` for `Z` standard logistic): `F(x) = ((1 - e^(-lambda*x)) / (1 + e^(-lambda*x)))^alpha`, for `x >= 0`. It is the exact structural analogue, within the half-logistic family, of [[generalized-exponential]] within the exponential family (Gupta-Kundu exponentiated exponential): both exponentiate a base cdf by a shape parameter `alpha` while `lambda` controls the underlying rate.

  Its density is `f(x) = alpha * lambda * ((1-e^(-lambda*x))/(1+e^(-lambda*x)))^(alpha-1) * 2*e^(-lambda*x) / (1+e^(-lambda*x))^2`, and its quantile function is `icdf(p) = (2/lambda) * atanh(p^(1/alpha))`, both closed-form.

  When `alpha = 1`, it reduces exactly to the ordinary [[half-logistic]] distribution with rate `lambda` (i.e. `scale = 1/lambda`).

  Parameters (single, optional map):

  - `alpha` (double): shape parameter, strictly positive. Default: `1.0`.
  - `lambda` (double): rate parameter, strictly positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Note: `pdf(0)` follows the usual `0^(alpha-1)` convention for exponentiated families: it is `0` for `alpha > 1`, `lambda/2` for `alpha = 1` (matching the ordinary half-logistic), and `##Inf` for `alpha < 1`.

  There is no elementary closed form for `mean`/`variance` (unlike [[generalized-exponential]]'s digamma-based moments); they are computed once, lazily, via numerical integration of `x * pdf(x)` (and `x^2 * pdf(x)`) over `x` in `[0, Infinity)`.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[half-logistic]], [[generalized-exponential]], [[generalized-logistic]]."
  ([] (generalized-half-logistic nil))
  ([{:keys [^double alpha ^double lambda rng]
     :or {alpha 1.0 lambda 1.0}}]
   (distr/generalized-half-logistic alpha lambda rng)))

(add-distr-method generalized-half-logistic)
(add-distr-method generalized-half-logistic :ghl)

;;

(defn truncated
  "Creates a truncated version of an existing distribution object, restricted to a `[left, right]` interval.

  Given a base distribution `distr`, this conditions it on lying within `[left, right]`: the resulting density is the density of `distr` renormalized by the probability mass `distr` originally assigned to that interval, and is zero outside of it. Sampling is done by rejection: `distr` is repeatedly sampled until a value falling inside `[left, right]` is obtained. It is useful for restricting the support of an otherwise unbounded or wide distribution, for example to build a bounded prior in Bayesian modelling.

  Parameters (single, optional map):

  - `distr`: the distribution to truncate. Default: a standard [[normal]] distribution.
  - `left` (double): lower truncation bound. Default: `nil`, meaning `distr`'s own [[lower-bound]] is used, i.e. no truncation on the left.
  - `right` (double): upper truncation bound. Default: `nil`, meaning `distr`'s own [[upper-bound]] is used, i.e. no truncation on the right.
  - `rng`: random number generator, used only to build the default `distr` when it is not supplied; ignored when `distr` is given explicitly (in that case `distr` itself is reused as the source of randomness).

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  `mean` and `variance` are computed by numerical (Gauss-Kronrod) integration over `[left, right]`, so construction may be relatively slow and their accuracy depends on the smoothness of `distr`'s density on that interval.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[normal]]."
  ([] (truncated nil))
  ([{:keys [distr left right rng]}]
   (let [distr (or distr (normal {:rng rng}))]
     (distr/truncated distr left right))))

(add-distr-method truncated)

(defn mixture
  "Creates a finite mixture distribution object combining several component distributions.

  A mixture distribution samples one of several component distributions, `distrs`, each chosen with its own probability (proportional to the corresponding entry in `weights`), and then samples from the selected component. Its density, cumulative distribution function and mean are therefore the `weights`-weighted average of those of the individual components, while its variance additionally accounts for the variability between component means, via the law of total variance. Components need not share the same type, domain (discrete or continuous) or support, making mixtures a flexible way to model multimodal or heterogeneous data.

  Parameters (single, optional map):

  - `distrs` (sequence of distributions): the component distributions to mix. Default: a single standard [[normal]] distribution.
  - `weights` (sequence of doubles): relative weight of each entry in `distrs`, in the same order; weights need not sum to `1.0`, as they are normalized internally. Default: equal weight for every component.
  - `rng`: random number generator used to select components and for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values, equivalent to a plain standard normal distribution.

  `continuous?` is `true` as soon as at least one component is continuous. `icdf` is obtained by numerically inverting `cdf` via root finding, so it may be relatively slow to evaluate and only approximately accurate.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[categorical]], [[truncated]]."
  ([] (mixture nil))
  ([{:keys [distrs weights rng]}]
   (let [distrs (or distrs [(normal {:rng rng})])
         weights (or weights (repeat (count distrs) 1.0))]
     (distr/mixture distrs weights rng))))

(add-distr-method mixture)

;;

(defonce ^{:doc "Default normal distribution (u=0.0, sigma=1.0)."} default-normal (normal))


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

(defn discrete?
  "Does distribution support discrete domain?"
  [d] (not (prot/continuous? d)))

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

(defn observe
  "Log likelihood of samples. Alias for [[log-likelihood]]."
  ^double [d vs]
  (log-likelihood d vs))

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
  "Distribution parameters"
  [d]
  (let [d' (if (keyword? d) (distribution d) d)]
    (prot/distribution-parameters d')))


;;

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
  "Sets seed."
  ([] (prot/set-seed! default-rng (lrand)))
  ([^long v] (prot/set-seed! default-rng v))
  ([rng ^long v] (prot/set-seed! rng v)))

;;

(defn- uniform-spacings
  ([^long n] (uniform-spacings default-rng n))
  ([rng ^long n]
   (let [xs (reductions m/+ (repeatedly (inc n) #(- (m/log (drandom rng)))))
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
  "Returns a lazy sequence of random samples.

  Works with either a random number generator or a distribution as `rng`, producing raw uniform draws or values following the given distribution, respectively.

  Parameters:

  - `rng` (optional): a random number generator or a distribution to sample from. Defaults to the default random number generator.
  - `n` (optional, long): limits the returned sequence to `n` values. When omitted, an infinite lazy sequence is returned.
  - `sampling-method` (optional, keyword): selects a sampling scheme applied to the underlying uniform draws, instead of plain independent sampling. Requires `n` to be given; `nil` falls back to plain sampling. One of:
    - `:uniform` - order statistics of `n` uniform draws, simulating a sorted uniform sample.
    - `:systematic` - low-variance systematic sampling using a single shared random offset.
    - `:stratified` - stratified sampling, one random draw per equal-width stratum.
    - `:antithetic` - antithetic sampling, pairing each draw `r` with its complement `1 - r`.
    - `:r2`, `:sobol`, `:halton` - jittered low-discrepancy sequences.
    When `rng` is a distribution, the resulting uniform values are transformed through its inverse cumulative distribution function; otherwise they are returned as-is.

  Returns a lazy sequence of samples.

  See also [[white-noise]], [[distribution?]], [[icdf]]."
  ([] (prot/->seq default-rng))
  ([rng] (prot/->seq rng))
  ([rng n] (prot/->seq rng n))
  ([rng n sampling-method]
   (if-not sampling-method
     (->seq rng n)
     (if (spacings sampling-method)
       (if (distribution? rng)
         (map (partial icdf rng) ((spacings sampling-method) n))
         ((spacings sampling-method) rng n))
       (throw (ex-info "Wrong sampling method" {:sampling-method sampling-method}))))))

(defn white-noise
  "Generates Gaussian white noise.

  The generated sequence consists of independent random numbers drawn from a normal
  distribution with mean 0.0 and a specified standard deviation.

  Parameters:

  - `sigma` (double, optional): The standard deviation of the generated noise.
  Defaults to 1.0 (standard normal distribution).
  - `rng` (fastmath.random.IRandomGenerator, optional): The random number generator to use.
  Defaults to a new JDK generator (`fastmath.random/rng :jdk`).

  Returns a lazy sequence of random numbers."
  ([] (white-noise 1.0))
  ([^double sigma] (white-noise (rng :jdk) sigma))
  ([rng ^double sigma] (repeatedly #(grandom rng sigma))))

;; arfima

(defn- ma-or-ar
  [ma? coeffs signal]
  (let [coeffs (if (number? coeffs) [coeffs] (vec (reverse coeffs)))]
    (if (every? m/zero? coeffs)
      signal
      (let [cnt (count coeffs)
            [es rsignal] (split-at cnt signal)]
        (letfn [(step [sig history]
                  (when (seq sig)
                    (let [e (double (first sig))
                          a (m/+ (v/dot history coeffs) e)]
                      (cons a (lazy-seq (step (rest sig) (subvec (conj history (if ma? e a)) 1)))))))]
          (drop cnt (lazy-seq (step rsignal (vec es)))))))))

(defn ma
  "Generates a Moving Average (MA) stochastic process.

  An MA(q) process expresses the current value as a linear combination of the current and past white noise error terms `et`. A single `theta` coefficient produces an MA(1) process `Yt = et + theta1 et-1`; a sequence of q coefficients produces an MA(q) process `Yt = et + theta1 et-1 + theta2 et-2 + ... + thetaq et-q`.

  Parameters:

  - `theta` (optional, number or sequence of numbers): the MA coefficient(s). A single number produces an MA(1) process; a sequence `[theta1 theta2 ... thetaq]` produces an MA(q) process, where `theta[i]` is the coefficient for `et-(i+1)`. Defaults to `0.0` (pure white noise, `signal` is returned unchanged).
  - `signal` (optional, sequence of numbers): the underlying white noise series (`et`). Defaults to a sequence of standard normal random numbers, see [[white-noise]].

  Returns a lazy sequence representing the generated MA process. Leading elements of `signal` are consumed to seed the initial history and to warm up the recursion, so the returned sequence is shorter than `signal`.

  See also [[ar]], [[arma]], [[fi]], [[arfima]], [[white-noise]]."
  ([] (ma 0.0))
  ([theta] (ma theta (white-noise)))
  ([theta signal] (ma-or-ar true theta signal)))

(defn ar
  "Generates an Autoregressive (AR) stochastic process.

  An AR(p) process expresses the current value as a linear combination of the current white noise error term and past values of the process itself. A single `phi` coefficient produces an AR(1) process `Yt = phi1 Yt-1 + et`; a sequence of p coefficients produces an AR(p) process `Yt = phi1 Yt-1 + phi2 Yt-2 + ... + phip Yt-p + et`.

  Parameters:

  - `phi` (optional, number or sequence of numbers): the AR coefficient(s). A single number produces an AR(1) process; a sequence `[phi1 phi2 ... phip]` produces an AR(p) process, where `phi[i]` is the coefficient for `Yt-(i+1)`. Defaults to `0.0` (pure white noise, `signal` is returned unchanged).
  - `signal` (optional, sequence of numbers): the driving white noise series (`et`). Defaults to a sequence of standard normal random numbers, see [[white-noise]].

  Returns a lazy sequence representing the generated AR process. Leading elements of `signal` are consumed to seed the initial history and to warm up the recursion, so the returned sequence is shorter than `signal`.

  See also [[ma]], [[arma]], [[fi]], [[arfima]], [[white-noise]]."
  ([] (ar 0.0))
  ([phi] (ar phi (white-noise)))
  ([phi signal] (ma-or-ar false phi signal)))

(defn arma
  "Generates an Autoregressive Moving Average (ARMA) stochastic process.

  An ARMA(p,q) process combines an MA(q) filter and an AR(p) filter: the driving `signal` is first passed through an MA(q) filter using `theta` (see [[ma]]), and the resulting series is then passed through an AR(p) filter using `phi` (see [[ar]]), giving `Yt = phi1 Yt-1 + ... + phip Yt-p + et + theta1 et-1 + ... + thetaq et-q`.

  Parameters:

  - `phi` (optional, number or sequence of numbers): the AR coefficient(s), see [[ar]]. Defaults to `0.0`.
  - `theta` (optional, number or sequence of numbers): the MA coefficient(s), see [[ma]]. Defaults to `0.0`.
  - `signal` (optional, sequence of numbers): the underlying white noise series. Defaults to a sequence of standard normal random numbers, see [[white-noise]].

  Returns a lazy sequence representing the generated ARMA process. As with [[ma]] and [[ar]], leading elements of `signal` are consumed while seeding both filters, so the returned sequence is noticeably shorter than `signal`.

  See also [[ma]], [[ar]], [[fi]], [[arfima]], [[white-noise]]."
  ([] (arma 0.0 0.0))
  ([phi theta] (arma phi theta (white-noise)))
  ([phi theta signal] (->> signal (ma theta) (ar phi))))

(defn- i-diffs
  [^long d]
  (->> (range d)
       (map (fn [^long k]
              (m/* (if (m/even? k) 1.0 -1.0)
                   (m/combinations d (m/long-inc k)))))))

(defn- fi-diffs
  [^double d ^long limit]
  (-> (reduce (fn [buff ^long k]
                (conj buff (m/* (m/- (double (if (m/zero? k) -1.0 (buff (m/dec k)))))
                                (m// (m/- d k) (m/inc k))))) [] (range limit))
      (vec)))

(defn fi
  "Generates a Fractionally Integrated (FI) process, or applies fractional integration to a series.

  Filters `signal` through the fractional differencing operator `(1-B)^-d`, where `B` is the backshift operator. Positive `d` adds long-range dependence (long memory) to the series, negative `d` produces anti-persistent, over-differenced behavior, and an integer `d` reduces to plain cumulative summation of order `d`.

  Parameters:

  - `d` (double, optional): the fractional differencing order. Defaults to `0.0` (identity, `signal` is returned unchanged).
  - `dlimit` (long, optional): for non-integer `d`, truncates the fractional differencing filter to this many lag coefficients (computed with Hosking's recursive formula). Ignored when `d` is an integer. When called with two arguments they are interpreted as `d signal` and `dlimit` defaults to `0`; when called with three arguments they are interpreted as `d dlimit signal`.
  - `signal` (optional, sequence of numbers): the input series to integrate. Defaults to a sequence of standard normal random numbers, see [[white-noise]].

  Returns a lazy sequence representing the fractionally integrated series. Integer `d` is handled by an exact finite expansion using binomial coefficients, applied with [[ar]]. Non-integer `d` with `dlimit` greater than `0` uses a truncated approximation, also applied with [[ar]]. Non-integer `d` with `dlimit` equal to `0` (the default) uses an exact recursive expansion whose coefficient count keeps growing as more of `signal` is consumed, trading speed for accuracy. In every non-identity case, leading elements of `signal` are consumed to seed the initial history, so the returned sequence is shorter than `signal`.

  See also [[ar]], [[ma]], [[arma]], [[arfima]]."
  ([] (fi 0.0))
  ([^double d] (fi d (white-noise)))
  ([^double d signal] (fi d 0 signal))
  ([^double d ^long dlimit signal]
   (cond
     (m/zero? d) signal
     (m/integer? d) (ar (i-diffs (long d)) signal)
     (m/pos? dlimit) (ar (fi-diffs d dlimit) signal)
     :else (letfn [(step [sig history coeffs ^long k]
                     (when (seq sig)
                       (let [phi1 (double (if (m/zero? k) -1.0 (coeffs (m/dec k))))
                             phi2 (m/* (m/- phi1)
                                       (m// (m/- d k) (m/inc k)))
                             ncoeffs (conj coeffs phi2)
                             a (m/+ (v/dot history ncoeffs) (double (first sig)))
                             nhistory (conj history a)]
                         (cons a (lazy-seq (step (rest sig) nhistory ncoeffs (m/inc k)))))))]
             (lazy-seq (step (rest signal) (list (first signal)) [] 0))))))

(defn arfima
  "Generates an Autoregressive Fractionally Integrated Moving Average (ARFIMA) stochastic process.

  An ARFIMA(p,d,q) process combines an MA(q) filter, fractional integration of order `d` and an AR(p) filter: the driving `signal` is first passed through an MA(q) filter using `theta` (see [[ma]]), then through the fractional integration filter using `d` and `dlimit` (see [[fi]]), and finally through an AR(p) filter using `phi` (see [[ar]]).

  Parameters:

  - `phi` (optional, number or sequence of numbers): the AR coefficient(s), see [[ar]]. Defaults to `0`.
  - `d` (double, optional): the fractional differencing order, see [[fi]]. Defaults to `0`.
  - `dlimit` (long, optional): truncation limit for the fractional integration step when `d` is not an integer, see [[fi]]. Defaults to `0` (no truncation). When called with four arguments they are interpreted as `phi d theta signal`; when called with five arguments they are interpreted as `phi d dlimit theta signal`.
  - `theta` (optional, number or sequence of numbers): the MA coefficient(s), see [[ma]]. Defaults to `0`.
  - `signal` (optional, sequence of numbers): the underlying white noise series. Defaults to a sequence of standard normal random numbers, see [[white-noise]].

  Returns a lazy sequence representing the generated ARFIMA process. As with [[ma]], [[fi]] and [[ar]], leading elements of `signal` are consumed while seeding each of the three filters in turn, so the returned sequence is noticeably shorter than `signal`.

  See also [[ma]], [[ar]], [[fi]], [[arma]], [[white-noise]]."
  ([] (arfima 0 0 0))
  ([phi d theta] (arfima phi d theta (white-noise)))
  ([phi d theta signal] (arfima phi d 0 theta signal))
  ([phi d dlimit theta signal] (->> signal (ma theta) (fi d dlimit) (ar phi))))

(m/unuse-primitive-operators)
