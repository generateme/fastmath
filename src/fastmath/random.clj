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

  To create your RNG use [[rng]] multimethod. Pass RNG name and (optional) integer seed; an unknown RNG name throws. Returned RNG is equipped with the `fastmath.protocols/RNGProto` protocol with methods: [[irandom]], [[lrandom]], [[frandom]] [[drandom]], [[grandom]], [[brandom]] which return random primitive value with given RNG. [[synced-rng]] creates a thread-safe generator, a plain generator must not be shared between threads.

  ```
  (let [rng (rng :isaac 1337)]
    (irandom rng))
  ```

  For conveniency default RNG (`:jdk`, [[default-rng]], one shared mutable object) with following functions are created: [[irand]], [[lrand]], [[frand]], [[drand]], [[grand]], [[brand]].

  Each prefix denotes returned type:

  * i - int (a long in Clojure)
  * l - long
  * f - float (a `Float` object)
  * d - double
  * g - gaussian (double)
  * b - boolean

  The draw functions are on the critical path and do not validate their arguments: a range or probability outside the documented one gives an unspecified result.

  Check individual function for parameters description.

  ### Random Vector Sequences

  Couple of functions to generate sequences of numbers or vectors.

  To create generator call [[sequence-generator]] with generator name and vector size.
  Following generators are available:

  * `:halton` - Halton low-discrepancy sequence (the weighted Kocis-Whiten variant of Apache Commons Math, equal to the standard one only in 1d); range [0,1), starts at the origin
  * `:sobol` - Sobol low-discrepancy sequence; range [0,1), starts at the origin
  * `:r2` - R2 low-discrepancy sequence; range [0,1), [more...](http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/)
  * `:sphere` - uniformly random distributed on unit sphere
  * `:ball` - uniformly random distributed from unit ball
  * `:gaussian` - gaussian distributed (mean=0, stddev=1)
  * `:default` - uniformly random; range:[0,1)

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

  Noise creation requires detailed configuration which is simple map of following keys, all optional (a `nil` value is a missing key, unknown keys are ignored, an invalid value throws `ex-info`):

  * `:seed` - seed as integer; the noise uses its low 32 bits (an integer outside of the 32-bit range is truncated); when missing, a random integer is drawn from the shared [[default-rng]]
  * `:noise-type` - type of noise, a key of [[noise-types]]: `:value`, `:gradient` (default), `:simplex`
  * `:interpolation` - type of interpolation (for value and gradient), a key of [[noise-interpolations]]: `:none`, `:linear`, `:hermite` (default) or `:quintic`
  * `:octaves` - number of octaves for combined noise (like FBM), an integer from 1 to 64, default: 6
  * `:lacunarity` - frequency scaling factor for combined noise, a finite number > 0, default: 2.00
  * `:gain` - amplitude scaling factor for combined noise, a finite number >= 0, default: 0.5 (for [[ridgedmulti-noise]] it divides the feedback and must be > 0)
  * `:normalize?` - should be normalized to `[0,1]` range (true, default) or to `[-1,1]` range (false)

  The values are in the given range but the realised range is narrower (the default [[noise]] in 2-D covers about `[0.16, 0.84]`). The noise of 1, 2 and 3 dimensions are unrelated functions, the 1-D noise has the period 256, and a coordinate which is `NaN`, infinite or has an absolute value above `2^30 / lacunarity^(octaves-1)` (`3.36e7` for the defaults) gives `NaN`.

  For usage convenience 3 ready to use functions are prepared, they use the fixed seed 1337, so the values are the same in every session. Returning value from `[0,1]` range:

  * [[noise]] - gradient (Perlin-like) noise (6 octaves, quintic interpolation)
  * [[vnoise]] - Value Noise (as in Processing, 6 octaves, hermite interpolation)
  * [[simplex]] - Simplex Noise (6 octaves)

  For random noise generation you can use [[random-noise-cfg]] and [[random-noise]]. Both can be feed with configuration. Additional configuration:

  * `:generator` can be set to one of the noise variants (`:single`, `:fbm`, `:billow`, `:ridgemulti`), always random when absent
  * `:warp-scale` - 0.0 - do not warp, >0.0 warp
  * `:warp-depth` - depth for warp, an integer >= 0 (random when absent, 1 most of the time), used if warp-scale is positive
  * `:rng` - random number generator which makes all the random choices, a fresh generator when missing; the shared [[default-rng]] is not used
  
  #### Discrete Noise

  [[discrete-noise]] is a 1d or 2d hash function for given integers (32-bit, the 2-d hash is of `X + 57Y` and collides). Returns double from `[0,1]` range.

  ### Distribution

  Various real and integer distributions (129 keys, see [[distributions-list]]), backed by Apache Commons Math, SSJ or implemented in this library. The functions working on a distribution are the proxies of `fastmath.protocols/DistributionProto` and `fastmath.protocols/RNGProto` defined in this namespace: [[pdf]], [[lpdf]], [[cdf]], [[ccdf]], [[icdf]], [[probability]], [[sample]], [[mean]], [[variance]], [[lower-bound]], [[upper-bound]], [[log-likelihood]] and others (the multivariate distributions have [[means]] and [[covariance]]).

  To create distribution call [[distribution]] multimethod with name as a keyword and map as parameters, or call the function of the same name as the key. The random number generator of a new distribution is a fresh unseeded one unless the `:rng` key is given; [[set-seed!]] reseeds a distribution.

  The numeric arguments of [[cdf]], [[ccdf]], [[pdf]], [[lpdf]], [[probability]] are checked for `NaN` and infinities before the call, [[icdf]] throws for a probability outside of `[0, 1]`."  
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.protocols :as prot]
            [fastmath.random.distributions :as distr]
            [fastmath.special :as special]
            [fastmath.stats.bins :as bins])
  (:import [org.apache.commons.math3.random RandomGenerator ISAACRandom JDKRandomGenerator MersenneTwister
            Well512a Well1024a Well19937a Well19937c Well44497a Well44497b
            RandomVectorGenerator HaltonSequenceGenerator SobolSequenceGenerator UnitSphereRandomVectorGenerator
            EmpiricalDistribution]
           [fastmath.java Array R2 SyncedRandomGenerator]
           [fastmath.java.dist FFixed GumbelFixed LevyFixed LogisticFixed NakagamiFixed]
           [umontreal.ssj.probdist AndersonDarlingDist AndersonDarlingDistQuick BetaSymmetricalDist
            InverseGammaDist  ChiDist ChiSquareNoncentralDist CramerVonMisesDist ErlangDist FatigueLifeDist FoldedNormalDist FrechetDist HalfNormalDist HyperbolicSecantDist HypoExponentialDist HypoExponentialDistEqual JohnsonSBDist JohnsonSLDist JohnsonSUDist KolmogorovSmirnovDist KolmogorovSmirnovDistQuick KolmogorovSmirnovPlusDist LoglogisticDist Pearson6Dist PowerDist RayleighDist WatsonGDist WatsonUDist]
           [fastmath.java.noise Billow RidgedMulti FBM NoiseConfig Noise Discrete]
           [org.apache.commons.math3.distribution BetaDistribution CauchyDistribution ChiSquaredDistribution ConstantRealDistribution EnumeratedRealDistribution ExponentialDistribution FDistribution GammaDistribution, GumbelDistribution, LaplaceDistribution, LevyDistribution, LogisticDistribution, LogNormalDistribution, NakagamiDistribution, NormalDistribution, ParetoDistribution, TDistribution, TriangularDistribution, UniformRealDistribution WeibullDistribution]
           [org.apache.commons.math3.distribution BinomialDistribution EnumeratedIntegerDistribution, GeometricDistribution, HypergeometricDistribution, PascalDistribution, PoissonDistribution, UniformIntegerDistribution, ZipfDistribution]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)
(m/use-primitive-operators)

(defn- invalid-input
  "Returns `ex-info` with `:type :invalid-input`, the name of the failing function, the property, the offending value and the description of expected values."
  [function property value expected message]
  (ex-info message {:type :invalid-input
                    :function function
                    :property property
                    :value value
                    :expected expected}))

(defn- long-integer?
  "Is `x` an integer which fits into a long?"
  [x]
  (or (instance? Long x)
      (instance? Integer x)
      (instance? Short x)
      (instance? Byte x)
      (and (integer? x) (< (.bitLength (biginteger x)) 64))))

(defn- check-seed
  "Returns the seed as a long or throws `ex-info` when `seed` is not an integer within the long range. `function` names the public function for the error."
  ^long [function seed]
  (if (long-integer? seed)
    (.longValue ^Number seed)
    (throw (invalid-input function :seed seed "an integer within the long range"
                          (str "Seed must be an integer within the long range, got " (pr-str seed))))))

(defn- check-parameter!
  "Throws `ex-info` (`:type :invalid-input`) naming the `function` and the `property` when `ok?` is false. `expected` describes the valid values."
  [function property value ok? expected]
  (when-not ok?
    (throw (invalid-input function property value expected
                          (str (name function) ", " (name property) ": expected " expected ", got " (pr-str value))))))

(defn- finite-number?
  "Is `x` a number with a finite double value?"
  [x]
  (and (number? x) (Double/isFinite (double x))))

(defn- int-parameter
  "Returns the integer parameter `value` of the constructor `function` as a long. `nil` gives `default`.

  The parameter has to be an integer value (an integral double like `20.0` is accepted) within the `int` range, because the backing classes narrow it to `int`. Throws `ex-info` (`:type :invalid-input`, `:property` the parameter) otherwise: a non-integer, `NaN`, infinity, a value out of the `int` range or a non-number."
  ^long [function property value ^long default]
  (if (nil? value)
    default
    (do (check-parameter! function property value
                          (and (finite-number? value)
                               (m/== (double value) (m/floor (double value)))
                               (m/<= (double Integer/MIN_VALUE) (double value) (double Integer/MAX_VALUE)))
                          "an integer within the int range [-2^31, 2^31-1]")
        (long value))))

(defn- long-parameter
  "Returns the integer parameter `value` of the constructor `function` as a long. `nil` gives `default`.

  The parameter has to be an integer value (an integral double like `20.0` is accepted) within the long range. Throws `ex-info` (`:type :invalid-input`, `:property` the parameter) otherwise: a non-integer, `NaN`, infinity, a value out of the long range or a non-number."
  ^long [function property value ^long default]
  (if (nil? value)
    default
    (do (check-parameter! function property value
                          (and (finite-number? value)
                               (m/== (double value) (m/floor (double value)))
                               (m/< (m/abs (double value)) 9.223372036854775807E18))
                          "an integer within the long range")
        (long value))))

;; Helper macro which creates RNG object of given class and/or seed.
(defmacro ^:private create-object-with-seed
  "Create object of the class with (or not) given seed. `nil` seed means unseeded generator, any other seed has to be an integer. Used to create RNG."
  [cl seed]
  `(let [arg# ~seed]
     (if (nil? arg#)
       (new ~cl)
       (prot/set-seed! (new ~cl) (check-seed ~''rng arg#))))) ;; seeding via protocols

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
  - `seed` (optional, integer): initial seed, any integer within the long range (`0` and negative values included). When given, the generator produces a fully reproducible sequence of values for a given `rng-name`; when omitted or `nil`, the generator seeds itself (implementation-dependent, typically from system entropy or clock). `:jdk` uses only the low 48 bits of the seed (seeds `s` and `s + 2^48` give the same stream), the other algorithms use all 64 bits.

  Returns a new, mutable generator object implementing the `fastmath.protocols/RNGProto` protocol (`irandom`, `lrandom`, `frandom`, `drandom`, `grandom`, `brandom`, `set-seed`, `set-seed!`, `->seq`), usable wherever an `rng` parameter is expected throughout the namespace, for example `(irandom (rng :isaac 1337))`.

  Thread safety: a generator is not safe to be shared between threads, only a single `:jdk` call is atomic and multi-step draws (`lrandom`, `drandom`) of any algorithm can interleave. Create one generator per thread, see [[child-rngs]], or use [[synced-rng]].

  Throws `ex-info` (`:type :invalid-input`, with `:expected` set to [[rngs-list]]) when `rng-name` is not a registered name, e.g. `nil`, a string or a misspelled keyword, and `ex-info` (`:property :seed`) when `seed` is neither `nil` nor an integer within the long range: `false`, a double, a ratio, a string, `NaN` and a big integer beyond the long range are rejected.

  See also [[rngs-list]], [[synced-rng]], [[default-rng]], [[child-rngs]]."
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
(declare rngs-list)

(defmethod rng :default [rng-name & [seed]]
  (if (= :default rng-name)
    (rng :jdk seed)
    (throw (invalid-input 'rng :rng-name rng-name rngs-list
                          (str "Unknown random number generator name " (pr-str rng-name)
                               ", expected one of " (pr-str rngs-list))))))

(defn synced-rng
  "Creates a thread-safe random number generator for the given algorithm and optional seed.

  Wraps the generator created by [[rng]] in a synchronized `fastmath.java.SyncedRandomGenerator` (a Commons Math `SynchronizedRandomGenerator` which also exposes the wrapped generator). Every call draws under a lock, so a generator can be shared between threads: the values drawn by several threads are the values of the single-thread stream of the same seed, only their distribution among the threads is not defined.

  Parameters:

  - `m` (keyword): algorithm name, see [[rng]] and [[rngs-list]].
  - `seed` (optional, integer): initial seed, as for [[rng]].

  Returns the synchronized generator, usable wherever an `rng` is expected. [[set-seed]] and [[child-rngs]] keep the algorithm and the synchronization: the new generators are synchronized generators of the same algorithm. [[set-seed!]] reseeds the generator in place.

  Throws `ex-info` for an unknown algorithm name or an invalid seed, see [[rng]].

  See also [[rng]], [[child-rngs]], [[ensure-rng]]."
  ([m] (SyncedRandomGenerator. (rng m)))
  ([m seed] (SyncedRandomGenerator. (rng m seed))))

;; List of randomizers
(defonce ^{:doc "Sequence of the names of all random number generators registered by [[rng]] (`:default` excluded).

The sequence is a snapshot taken when the namespace is loaded, in no particular order: algorithms added later with `defmethod` on [[rng]] do not appear in it. Each name can be passed to [[rng]] and [[synced-rng]]."}
  rngs-list (remove #{:default} (keys (methods rng))))

;; protocol proxies
(defn frandom
  "Draws a random float with the given generator.

  Parameters:

  - `rng`: a random number generator, see [[rng]].
  - `mx` (optional, number): upper bound, the range is `[0, mx)`.
  - `mn` (optional, number): lower bound, the range is `[mn, mx)`. Give `mn` and `mx` together.

  Without bounds the range is `[0, 1)`. Returns a `java.lang.Float` (not a double), the resolution is 23 bits for the Commons Math generators and 24 bits for `:jdk`; use [[drandom]] for a double.

  The function is on the critical path and does not validate its arguments: `mx <= mn`, infinite or `NaN` bounds give an unspecified result, and for a span of a few float ulps the rounded result can reach `mx`.

  See also [[drandom]], [[frand]]."
  ([rng] (prot/frandom rng))
  ([rng mx] (prot/frandom rng mx))
  ([rng mn mx] (prot/frandom rng mn mx)))

(defn drandom
  "Draws a random double with the given generator.

  Parameters:

  - `rng`: a random number generator, see [[rng]]. A distribution is accepted in the 1-argument form only and then draws a variate of the distribution.
  - `mx` (optional, number): upper bound, the range is `[0, mx)`.
  - `mn` (optional, number): lower bound, the range is `[mn, mx)`. Give `mn` and `mx` together.

  Without bounds the range is `[0, 1)`; `0.0` is possible, with probability 2^-52 for Commons Math generators and 2^-53 for `:jdk`. Returns a double.

  The function is on the critical path and does not validate its arguments: `mx < mn` gives values in `(mx, mn]`, a span which overflows or an infinite or `NaN` bound gives infinity or `NaN`, and for a span of a few ulps the rounded result can reach `mx`.

  See also [[frandom]], [[grandom]], [[drand]]."
  (^double [rng] (prot/drandom rng))
  (^double [rng mx] (prot/drandom rng mx))
  (^double [rng mn mx] (prot/drandom rng mn mx)))

(defn grandom
  "Draws a random number from the normal distribution with the given generator.

  Parameters:

  - `rng`: a random number generator, see [[rng]].
  - `stddev` (optional, number): standard deviation, the result is drawn from `N(0, stddev)`.
  - `mean`, `stddev` (optional, numbers): the result is `mean + stddev * g` with `g` drawn from `N(0, 1)`, so `stddev` is the standard deviation and not a bound, also for `mean` much larger than `stddev`.

  Without arguments the result is drawn from `N(0, 1)`. Returns a double. A negative `stddev` mirrors the distribution and `stddev = 0` returns `mean`, one number is drawn in every case.

  The function is on the critical path and does not validate its arguments: `NaN` or infinite parameters give `NaN` or infinity.

  See also [[drandom]], [[grand]]."
  (^double [rng] (prot/grandom rng))
  (^double [rng stddev] (prot/grandom rng stddev))
  (^double [rng mean stddev] (prot/grandom rng mean stddev)))

(defn irandom
  "Draws a random integer with the given generator.

  Parameters:

  - `rng`: a random number generator, see [[rng]].
  - `mx` (optional, integer): exclusive upper bound, the range is `[0, mx)`; must be positive.
  - `mn` (optional, integer): inclusive lower bound, the range is `[mn, mx)`. Give `mn` and `mx` together; `mx > mn`, and `mx = mn` returns `mn`.

  Without bounds the result is any 32-bit integer. The result is a long in Clojure. A range (the bound `mx`, or the span `mx - mn`) up to `Integer/MAX_VALUE` is drawn from the 32-bit generator, a wider range is drawn without bias from the 64-bit one as in [[lrandom]], so the full `[mn, mx)` range is covered also for bounds beyond 32 bits. Non-integer bounds are truncated.

  The function is on the critical path and does not validate its arguments: `mx <= 0` or `mx < mn` is invalid and gives an exception (the type depends on the generator) or an unspecified value, never a documented result.

  See also [[lrandom]], [[irand]]."
  (^long [rng] (prot/irandom rng))
  (^long [rng mx] (prot/irandom rng mx))
  (^long [rng mn ^long mx] (prot/irandom rng mn mx)))

(defn lrandom
  "Draws a random long with the given generator.

  Parameters:

  - `rng`: a random number generator, see [[rng]].
  - `mx` (optional, long): exclusive upper bound, the range is `[0, mx)`; must be positive. The bound is applied without modulo bias (draws from the incomplete tail of the 64-bit range are rejected).
  - `mn` (optional, long): inclusive lower bound, the range is `[mn, mx)`. Give `mn` and `mx` together; `mx > mn`, and `mx = mn` returns `mn`. Spans of 2^63 and more, up to the full `[Long/MIN_VALUE, Long/MAX_VALUE)`, are supported.

  Without bounds the result is any long (negative values included). Non-integer bounds are truncated.

  The function is on the critical path and does not validate its arguments: `mx = 0` throws `ArithmeticException`, `mx < 0` or `mx < mn` is invalid and gives an unspecified value in `(mx, 0]` or `(mx, mn]`.

  See also [[irandom]], [[lrand]]."
  (^long [rng] (prot/lrandom rng))
  (^long [rng mx] (prot/lrandom rng mx))
  (^long [rng mn mx] (prot/lrandom rng mn mx)))

(defn brandom
  "Draws a random boolean with the given generator.

  Parameters:

  - `rng`: a random number generator, see [[rng]].
  - `p` (optional, number): probability of `true`, drawn as `u < p` for `u` from `[0, 1)`.

  Without `p` the result is `true` or `false` with equal probability (one bit of the generator, so the stream differs from the one of `p = 0.5`). Returns a boolean.

  The function is on the critical path and does not validate `p`: `p <= 0` and `NaN` always give `false`, `p >= 1` always gives `true`. A probability below the resolution of the generator (2^-52) is not honoured.

  See also [[brand]], [[randval-rng]], [[flip-rng]]."
  ([rng] (prot/brandom rng))
  ([rng p] (prot/brandom rng p)))

(defn ensure-rng
  "Validates a random number generator given as an argument and replaces a missing one with a new generator.

  Functions which accept an optional `rng` use it to avoid the shared [[default-rng]], so a result can be reproduced with a seeded generator and a call never disturbs the stream of the shared one.

  Parameters:

  - `rng`: a random number generator created by [[rng]] or [[synced-rng]] (any Apache Commons Math `RandomGenerator`), or `nil`.

  Returns `rng` itself when it is a generator, or a new unseeded `JDKRandomGenerator` when `rng` is `nil`. A distribution is not accepted.

  Throws `ex-info` with `{:rng value}` for any other value.

  See also [[child-rngs]], [[rng]], [[synced-rng]]."
  ^RandomGenerator [rng]
  (cond
    (nil? rng) (JDKRandomGenerator.)
    (instance? RandomGenerator rng) rng
    :else (throw (ex-info (str "Expected a random number generator (RandomGenerator) or nil, got " (pr-str rng))
                          {:rng rng}))))

(defn- child-rng
  "Creates a new generator seeded with a number drawn from `rng`."
  [rng]
  (prot/set-seed rng (lrandom rng)))

(defn child-rngs
  "Creates new random number generators seeded from a given generator, for independent runs which must not share one mutable generator.

  Every child draws one seed from `rng` with [[lrandom]] and is created by `set-seed` of the protocol, so children are distinct objects and the same seeded `rng` gives the same children. The seeds are drawn when the children are created, in order, from the thread which creates them. Create the children before handing them to other threads.

  Parameters:

  - `rng`: a random number generator created by [[rng]] or [[synced-rng]] (any Apache Commons Math `RandomGenerator`), checked by [[ensure-rng]]. `nil` stands for a new generator.
  - `n` (optional, long): the number of children.

  Returns a lazy infinite sequence of children when `n` is not given, where realizing an element draws one seed. With `n` returns a vector of `n` children and draws exactly `n` seeds, none for `n` equal to `0`.

  A child has the class of `rng` for the generators created by [[rng]]. A generator created by [[synced-rng]] gives synchronized children (see [[synced-rng]]) of the same algorithm. A generator of another kind (a custom `RandomGenerator`) gives children of the default type (`JDKRandomGenerator`). `n` is truncated to a long.

  Throws `ex-info` when `rng` is neither a generator nor `nil` (a distribution included), or when `n` is negative.

  See also [[ensure-rng]], [[set-seed]], [[rng]]."
  ([rng]
   (let [rng (ensure-rng rng)]
     (repeatedly #(child-rng rng))))
  ([rng ^long n]
   (when (neg? n) (throw (ex-info "The number of children must not be negative" {:n n})))
   (let [rng (ensure-rng rng)]
     (vec (repeatedly n #(child-rng rng))))))

;; Type hinted functions generating random value
(defn- next-random-value-long
  "Generate next long.

  * arity 0 - any long
  * arity 1 - from 0 to provided integer (excluded), without modulo bias for positive bound
  * arity 2 - from the provided range (included, excluded), spans up to the full long range

  For a positive bound the draws `x` above `Long/MAX_VALUE - (2^64 mod mx)` are rejected, what leaves a block of `2^64 - (2^64 mod mx)` consecutive values, a whole number of residue cycles. The value is `(mod x mx)`, as without rejection. The exact test needs a division, the cheap sufficient one `x <= Long/MAX_VALUE - mx` does not."
  (^long [^RandomGenerator r] (.nextLong r))
  (^long [^RandomGenerator r ^long mx]
   (if (< 0 mx)
     (let [fast-limit (- Long/MAX_VALUE mx)]
       (loop [x (.nextLong r)]
         (if (or (<= x fast-limit)
                 (<= x (- Long/MAX_VALUE (Long/remainderUnsigned (- mx) mx))))
           (mod x mx)
           (recur (.nextLong r)))))
     (mod (.nextLong r) mx)))
  (^long [^RandomGenerator r ^long mn ^long mx]
   (let [diff (- mx mn)]
     (cond
       (zero? diff) mn
       (< 0 diff) (+ mn (next-random-value-long r diff))
       ;; span of 2^63 or more: diff wrapped to negative, it is the unsigned span
       (< mn mx) (loop [x (.nextLong r)]
                   (if (< (Long/compareUnsigned x diff) 0)
                     (+ mn x)
                     (recur (.nextLong r))))
       :else (+ mn (next-random-value-long r diff))))))

(defn- next-random-value-int
  "Generate next int.

  * arity 0 - from 0 to maximum int value
  * arity 1 - from 0 to provided integer (excluded), bounds above `Integer/MAX_VALUE` use the long generator
  * arity 2 - from the provided range (included, excluded), spans above `Integer/MAX_VALUE` use the long generator"
  (^long [^RandomGenerator r] (.nextInt r))
  (^long [^RandomGenerator r ^long mx]
   (if (<= mx Integer/MAX_VALUE)
     (.nextInt r mx)
     (next-random-value-long r mx)))
  (^long [^RandomGenerator r ^long mn ^long mx]
   (let [diff (- mx mn)]
     (cond
       (zero? diff) mn
       (and (<= 0 diff) (<= diff Integer/MAX_VALUE)) (+ mn (.nextInt r diff))
       (< mn mx) (next-random-value-long r mn mx)
       :else (+ mn (.nextInt r diff))))))

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
  * arity 2 - N(mean,std), computed as `mean + std * g`"
  (^double [^RandomGenerator r] (.nextGaussian r))
  (^double [^RandomGenerator r ^double mx] (* (.nextGaussian r) mx))
  (^double [^RandomGenerator r ^double mean ^double std]
   (+ mean (* std (.nextGaussian r)))))

;; Extend RandomGenerator interface with functions created by macro `next-random-value-fn`. This way all RNG classes are enriched with new, more convenient functions.
;;
;; Note that `grandom` is under special care due to different [mn mx] range meaning.

(extend RandomGenerator 
  prot/RNGProto
  {:irandom next-random-value-int
   :lrandom next-random-value-long
   :frandom next-random-value-float
   :drandom next-random-value-double
   :grandom next-random-value-gaussian
   :brandom (fn
              ([^RandomGenerator t] (.nextBoolean t))
              ([t ^double thr] (< (next-random-value-double t) thr)))
   :set-seed! (fn [^RandomGenerator t ^long seed]
                (.setSeed t seed)
                t)
   ;; synchronized generator gives synchronized generator of the same algorithm, unknown classes give the default `:jdk`
   :set-seed (fn [^RandomGenerator t ^long seed]
               (if (instance? SyncedRandomGenerator t)
                 (SyncedRandomGenerator. (prot/set-seed (.getDelegate ^SyncedRandomGenerator t) seed))
                 (rng (get rng-class->keyword (class t) :jdk) seed)))
   :->seq (fn
            ([^RandomGenerator t] (repeatedly #(next-random-value-double t)))
            ([^RandomGenerator t n] (repeatedly n #(next-random-value-double t))))})

;; ### Default RNG

(defonce ^{:doc "The default random number generator, a `:jdk` generator (`JDKRandomGenerator`) created when the namespace is loaded.

It is one shared, mutable object: [[irand]], [[lrand]], [[frand]], [[drand]], [[grand]], [[brand]], [[randval]], [[flip]], [[roll-a-dice]], [[set-seed!]] and [[set-seed]] without a generator use it, so seeding it (`(set-seed! 42)`) makes them reproducible, and the calls of all threads draw from the same stream. The functions of this namespace which take an optional `rng` never use it for a missing `rng`, see [[ensure-rng]].

The functions read the var at call time, so `with-redefs` replaces the generator."}
  default-rng (rng :jdk))

(defn brand
  "Draws a random boolean with the default generator [[default-rng]].

  Parameters:

  - `p` (optional, number): probability of `true`, drawn as `u < p` for `u` from `[0, 1)`.

  Without `p` the result is `true` or `false` with equal probability. Returns a boolean. The generator is read at call time, so `with-redefs` of [[default-rng]] changes the result.

  The function is on the critical path and does not validate `p`: `p <= 0` and `NaN` always give `false`, `p >= 1` always gives `true`.

  See also [[brandom]], [[randval]], [[flip]]."
  ([] (prot/brandom default-rng))
  ([p] (prot/brandom default-rng p)))

(defn frand
  "Draws a random float with the default generator [[default-rng]].

  Parameters:

  - `mx` (optional, number): upper bound, the range is `[0, mx)`.
  - `mn` (optional, number): lower bound, the range is `[mn, mx)`. Give `mn` and `mx` together.

  Without bounds the range is `[0, 1)`. Returns a `java.lang.Float` (not a double), use [[drand]] for a double.

  The function does not validate its arguments, see [[frandom]].

  See also [[frandom]], [[drand]]."
  ([] (prot/frandom default-rng))
  ([mx] (prot/frandom default-rng mx))
  ([mn mx] (prot/frandom default-rng mn mx)))

(defn drand
  "Draws a random double with the default generator [[default-rng]].

  Parameters:

  - `mx` (optional, double): upper bound, the range is `[0, mx)`.
  - `mn` (optional, double): lower bound, the range is `[mn, mx)`. Give `mn` and `mx` together.

  Without bounds the range is `[0, 1)`. Returns a double.

  The function does not validate its arguments, see [[drandom]].

  See also [[drandom]], [[frand]], [[grand]]."
  (^double [] (prot/drandom default-rng))
  (^double [^double mx] (prot/drandom default-rng mx))
  (^double [^double mn ^double mx] (prot/drandom default-rng mn mx)))

(defn grand
  "Draws a random number from the normal distribution with the default generator [[default-rng]].

  Parameters:

  - `stddev` (optional, double): standard deviation, the result is drawn from `N(0, stddev)`.
  - `mean`, `stddev` (optional, doubles): the result is drawn from `N(mean, stddev)`, computed as `mean + stddev * g`.

  Without arguments the result is drawn from `N(0, 1)`. Returns a double.

  The function does not validate its arguments, see [[grandom]].

  See also [[grandom]], [[drand]]."
  (^double [] (prot/grandom default-rng))
  (^double [^double stddev] (prot/grandom default-rng stddev))
  (^double [^double mean ^double stddev] (prot/grandom default-rng mean stddev)))

(defn irand
  "Draws a random integer with the default generator [[default-rng]].

  Parameters:

  - `mx` (optional, integer): exclusive upper bound, the range is `[0, mx)`; must be positive.
  - `mn` (optional, integer): inclusive lower bound, the range is `[mn, mx)`. Give `mn` and `mx` together; `mx > mn`.

  Without bounds the result is any 32-bit integer. The result is a long in Clojure; ranges wider than `Integer/MAX_VALUE` are supported, see [[irandom]].

  The function does not validate its arguments: an invalid range throws or gives an unspecified value.

  See also [[irandom]], [[lrand]]."
  (^long [] (prot/irandom default-rng))
  (^long [mx] (prot/irandom default-rng mx))
  (^long [mn mx] (prot/irandom default-rng mn mx)))

(defn lrand
  "Draws a random long with the default generator [[default-rng]].

  Parameters:

  - `mx` (optional, long): exclusive upper bound, the range is `[0, mx)`; must be positive.
  - `mn` (optional, long): inclusive lower bound, the range is `[mn, mx)`. Give `mn` and `mx` together; `mx > mn`, spans up to the full long range are supported.

  Without bounds the result is any long.

  The function does not validate its arguments: `mx = 0` throws `ArithmeticException`, an invalid range gives an unspecified value, see [[lrandom]].

  See also [[lrandom]], [[irand]]."
  (^long [] (prot/lrandom default-rng))
  (^long [^long mx] (prot/lrandom default-rng mx))
  (^long [^long mn ^long mx] (prot/lrandom default-rng mn mx)))

(defmacro randval
  "Chooses between two expressions, or draws a boolean, with the default generator [[default-rng]].

  The meaning of the arguments depends on their number:

  - `(randval)` - returns a boolean, `true` with probability `0.5`.
  - `(randval prob)` - returns a boolean (not a value), `true` with probability `prob`.
  - `(randval v1 v2)` - evaluates and returns `v1` or `v2` with equal probability; only the chosen expression is evaluated. Note that `(randval 0.3 :a)` is this two-argument form and returns `0.3` or `:a`, not a value with probability `0.3`.
  - `(randval prob v1 v2)` - evaluates and returns `v1` with probability `prob`, otherwise `v2`.

  The probability is not validated, see [[brand]]: `prob <= 0` and `NaN` always choose `v2`, `prob >= 1` always chooses `v1`.

  See also [[randval-rng]], [[flip]], [[flipb]], [[brand]]."
  ([v1 v2]
   `(if (prot/brandom default-rng) ~v1 ~v2))
  ([prob v1 v2]
   `(if (prot/brandom default-rng ~prob) ~v1 ~v2))
  ([prob]
   `(prot/brandom default-rng ~prob))
  ([]
   `(prot/brandom default-rng)))

(defn flip
  "Returns `1` with given probability and `0` otherwise, drawn with the default generator [[default-rng]].

  Parameters:

  - `p` (optional, number): probability of `1`, default `0.5`. It is not validated: `p <= 0` and `NaN` always give `0`, `p >= 1` always gives `1`.

  Returns a long.

  See also [[flipb]], [[flip-rng]], [[brand]]."
  (^long [p]
   (randval p 1 0))
  (^long []
   (randval 0.5 1 0)))

(defn flipb
  "Returns `true` with given probability and `false` otherwise, drawn with the default generator [[default-rng]].

  Parameters:

  - `p` (optional, number): probability of `true`, default `0.5`. It is not validated, see [[flip]].

  Returns a boolean.

  See also [[flip]], [[flipb-rng]], [[brand]]."
  ([p] (randval p))
  ([] (randval)))

(defn roll-a-dice
  "Rolls dice with the default generator [[default-rng]] and returns the sum of the dots.

  Parameters:

  - `sides` (integer): number of sides of a die, a roll is a number from `1` to `sides`; must be positive.
  - `dices` (optional, integer): number of dice, default `1`. The argument order is `(roll-a-dice dices sides)`. `dices <= 0` returns `0`.

  Returns the sum as a long. The arguments are not validated: a non-positive `sides` throws, a non-integer `dices` is rounded up.

  See also [[roll-a-dice-rng]], [[irand]]."
  (^long [sides]
   (inc (irand sides)))
  (^long [dices sides]
   (reduce m/+ (repeatedly dices #(inc (irand sides)))))  )

;; rng versions

(defmacro randval-rng
  "Chooses between two expressions, or draws a boolean, with the given generator.

  The meaning of the arguments depends on their number (`rng` is evaluated at every use):

  - `(randval-rng rng)` - returns a boolean, `true` with probability `0.5`.
  - `(randval-rng rng prob)` - returns a boolean (not a value), `true` with probability `prob`.
  - `(randval-rng rng v1 v2)` - evaluates and returns `v1` or `v2` with equal probability; only the chosen expression is evaluated.
  - `(randval-rng rng prob v1 v2)` - evaluates and returns `v1` with probability `prob`, otherwise `v2`.

  The probability is not validated, see [[brandom]].

  See also [[randval]], [[flip-rng]], [[flipb-rng]], [[brandom]]."
  ([rng v1 v2]
   `(if (prot/brandom ~rng) ~v1 ~v2))
  ([rng prob v1 v2]
   `(if (prot/brandom ~rng ~prob) ~v1 ~v2))
  ([rng prob]
   `(prot/brandom ~rng ~prob))
  ([rng]
   `(prot/brandom ~rng)))

(defn flip-rng
  "Returns `1` with given probability and `0` otherwise, drawn with the given generator.

  Parameters:

  - `rng`: a random number generator, see [[rng]].
  - `p` (optional, number): probability of `1`, default `0.5`. It is not validated, see [[flip]].

  Returns a long.

  See also [[flip]], [[flipb-rng]], [[brandom]]."
  (^long [rng p]
   (randval-rng rng p 1 0))
  (^long [rng]
   (randval-rng rng 0.5 1 0)))

(defn flipb-rng
  "Returns `true` with given probability and `false` otherwise, drawn with the given generator.

  Parameters:

  - `rng`: a random number generator, see [[rng]].
  - `p` (optional, number): probability of `true`, default `0.5`. It is not validated, see [[flip]].

  Returns a boolean.

  See also [[flipb]], [[flip-rng]], [[brandom]]."
  ([rng p] (randval-rng rng p))
  ([rng] (randval-rng rng)))

(defn roll-a-dice-rng
  "Rolls dice with the given generator and returns the sum of the dots.

  Parameters:

  - `rng`: a random number generator, see [[rng]].
  - `sides` (integer): number of sides of a die, a roll is a number from `1` to `sides`; must be positive.
  - `dices` (optional, integer): number of dice, default `1`. The argument order is `(roll-a-dice-rng rng dices sides)`. `dices <= 0` returns `0`.

  Returns the sum as a long. The arguments are not validated, see [[roll-a-dice]].

  See also [[roll-a-dice]], [[irandom]]."
  (^long [rng sides]
   (inc (irandom rng sides)))
  (^long [rng dices sides]
   (reduce m/+ (repeatedly dices #(inc (irandom rng sides))))))

;; generators

(defn- check-dimensions
  "Returns the number of dimensions as a long or throws `ex-info` (`:type :invalid-input`, `:property :dimensions`) when it is not an integer from 1 to the limit of the generator: 40 for `:halton`, 1000 for `:sobol`, 15 for `:r2`, no limit for the others. `function` names the public function for the error."
  ^long [function seq-generator dimensions]
  (let [limit (case seq-generator :halton 40 :sobol 1000 :r2 15 Long/MAX_VALUE)]
    (if (and (long-integer? dimensions)
             (<= 1 (.longValue ^Number dimensions) limit))
      (.longValue ^Number dimensions)
      (let [expected (if (== limit Long/MAX_VALUE)
                       "an integer >= 1"
                       (str "an integer from 1 to " limit))]
        (throw (invalid-input function :dimensions dimensions expected
                              (str "Number of dimensions for " (pr-str seq-generator) " must be " expected
                                   ", got " (pr-str dimensions))))))))

(defn- check-jitter
  "Returns the jitter as a double or throws `ex-info` (`:type :invalid-input`, `:property :jitter`) when it is not a finite number `>= 0`."
  ^double [jitter]
  (if (and (number? jitter)
           (Double/isFinite (double jitter))
           (>= (double jitter) 0.0))
    (double jitter)
    (throw (invalid-input 'jittered-sequence-generator :jitter jitter "a finite number >= 0"
                          (str "Jitter must be a finite number >= 0, got " (pr-str jitter))))))

(defn- ball-point
  "Draws a point from the interior of the unit ball of `dims` dimensions (already validated) using `rng` (already checked)."
  [rng ^long dims]
  (let [u (double-array (repeatedly (+ dims 2) #(grandom rng)))
        ^doubles n (v/div u (v/mag u))]
    (case dims
      1 (aget n 0)
      2 (v/array->vec2 n)
      3 (v/array->vec3 n)
      4 (v/array->vec4 n)
      (vec (take dims n)))))

;; http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/#more-2165
(defn ball-random
  "Draws a uniformly random point from the interior of a `dims`-dimensional unit ball.

  Uses the method described in this [article](http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/#more-2165): draws `dims + 2` independent `N(0,1)` samples, normalizes them to a unit vector (a uniformly random point on the `(dims + 2)`-dimensional unit sphere), and keeps only its first `dims` coordinates. Projecting a uniform sample from a sphere two dimensions higher down onto `dims` coordinates yields a point uniformly distributed by volume within the `dims`-ball, unlike naively rejecting or normalizing `dims`-dimensional samples, which biases the result.

  Parameters:

  - `rng` (optional): random number generator to draw the underlying gaussian samples from, checked by [[ensure-rng]]: `nil` stands for a new unseeded generator. When the argument is missing (the one-argument form) the shared [[default-rng]] is used.
  - `dims` (integer): dimensionality of the ball, an integer `>= 1`.

  Returns a double when `dims` is `1`, a `Vec2`, `Vec3` or `Vec4` when `dims` is `2`, `3` or `4`, or a plain vector of doubles otherwise; magnitude of the result is always less than or equal to `1.0`, every coordinate is in `[-1, 1]`.

  Throws `ex-info` (`:property :dimensions`) when `dims` is not an integer `>= 1` (also with assertions disabled), and `ex-info` with `{:rng value}` when `rng` is neither a generator nor `nil` (a distribution included).

  Used internally by the `:ball` method of [[sequence-generator]].

  See also [[sequence-generator]], [[ensure-rng]]."
  ([dims] (ball-random default-rng dims))
  ([rng dims]
   (let [dims (check-dimensions 'ball-random :ball dims)]
     (ball-point (ensure-rng rng) dims))))

(defn- rv-generators
  "Generators from commons math and custom classes. `rng` is used only by `:sphere`."
  [seq-generator dimensions ^RandomGenerator rng]
  (let [dimensions (check-dimensions 'sequence-generator seq-generator dimensions)
        ^RandomVectorGenerator g (case seq-generator
                                   :halton (HaltonSequenceGenerator. dimensions)
                                   :sobol (SobolSequenceGenerator. dimensions)
                                   :sphere (UnitSphereRandomVectorGenerator. dimensions rng)
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
  "Random generators drawing from `rng`."
  [seq-generator dimensions rng]
  (let [dimensions (check-dimensions 'sequence-generator seq-generator dimensions)
        g (if (= seq-generator :gaussian)
            #(grandom rng)
            #(drandom rng))]
    (repeatedly (case dimensions
                  1 g
                  2 (partial v/generate-vec2 g)
                  3 (partial v/generate-vec3 g)
                  4 (partial v/generate-vec4 g)
                  #(vec (repeatedly dimensions g))))))

;; jittering
;; http://extremelearning.com.au/a-simple-method-to-construct-isotropic-quasirandom-blue-noise-point-sequences/

;; Points of 1 to 4 dimensions are a number, `Vec2`, `Vec3` and `Vec4`. Points of more dimensions are plain vectors,
;; where the generic vector operations are slow, so the jitter of these points is an array and is added in one loop.
(def ^:private ^:const plain-vector-dimensions 5)

(defn- jitter-array
  "Returns `dimensions` numbers drawn from `rng`, each multiplied by `scale`, as a `double` array."
  ^doubles [^long dimensions ^double scale rng]
  (let [a (double-array dimensions)]
    (dotimes [d dimensions]
      (Array/set a d (m/* (drandom rng) scale)))
    a))

(defn- wrap-ball
  "Wraps `x` periodically into `[-1,1)`: `((x + 1) mod 2) - 1`."
  ^double [^double x]
  (m/- (mod (m/+ x 1.0) 2.0) 1.0))

(defn- add-jitter
  "Adds the `jitter` array to the point `v` (a vector) coordinate by coordinate, wraps the sums into [0,1) with `frac` (into [-1,1) when `ball?`) and returns a plain vector."
  [v ^doubles jitter ball?]
  (let [n (alength jitter)
        x (double-array n)]
    (if ball?
      (dotimes [d n]
        (Array/set x d (wrap-ball (m/+ (double (nth v d)) (Array/get jitter d)))))
      (dotimes [d n]
        (Array/set x d (m/frac (m/+ (double (nth v d)) (Array/get jitter d))))))
    (vec x)))

(defn- jitter-generator
  "Generate random jitter drawn from `rng`: a lazy sequence of numbers, vectors or, for `plain-vector-dimensions` and more dimensions, `double` arrays. Every element draws its `dimensions` numbers when it is realized, one element at a time."
  [seq-generator ^long dimensions ^double jitter rng]
  (let [[^double d0 ^double i0 ^double f ^double p] (case seq-generator
                                                      :r2 [0.76 0.7 0.25 -0.5]
                                                      :halton [0.9 0.7 0.25 -0.5]
                                                      :sobol [0.16 0.58 0.4 -0.2]
                                                      [0.5 0.5 0.25 -0.5])
        c (* jitter m/SQRTPI d0 f)
        scale (fn ^double [^long i] (* c (m/pow (- (inc i) i0) p)))]
    (if (m/>= dimensions plain-vector-dimensions)
      (map-indexed (fn [^long i _] (jitter-array dimensions (scale i) rng)) (repeat nil))
      (map-indexed (fn [^long i v] (v/mult v (scale i))) (random-generators :default dimensions rng)))))


;; Sequence creators

(defmulti
  ^{:doc "Creates a lazy, infinite sequence of random or quasi-random points.

  Dispatches on `seq-generator` to pick the sampling method; every point in the resulting sequence has `dimensions` components. See [[sequence-generators-list]] for the full list of registered `seq-generator` keys, and [[jittered-sequence-generator]] to add blue-noise jitter to `:r2`, `:halton` and `:sobol` sequences.

  Parameters:

  - `seq-generator`: keyword selecting the generator:
      - `:r2` - R2 low-discrepancy sequence (additive recurrence with the generalized golden ratio), points in `[0,1)` for each dimension; the first point is `0.5 + alpha`. The constants have a double rounding error which makes the points drift by about `3.6e-10` after `10^7` points (linear in the index).
      - `:halton` - Halton low-discrepancy sequence in the weighted variant of Kocis and Whiten as implemented in Apache Commons Math (the digits of the radical inverses are multiplied by weights). It is equal to the standard Halton sequence (for example the one of scipy) only for `dimensions = 1`; the points for more dimensions differ. Points in `[0,1)`.
      - `:sobol` - Sobol low-discrepancy sequence (Joe-Kuo direction numbers, equal to the unscrambled `scipy.stats.qmc.Sobol`), points in `[0,1)`.
      - `:default`, `:uniform` or any other unregistered keyword (a typo included) - independent uniform pseudo-random points from `[0,1)` for each dimension.
      - `:gaussian` - independent pseudo-random points, each component drawn from `N(0,1)`.
      - `:sphere` - pseudo-random points on the surface of a unit sphere (euclidean distance from origin equals `1.0` up to one ulp).
      - `:ball` - pseudo-random points uniformly distributed within a unit ball, every coordinate in `[-1,1]`.
  - `dimensions` (integer): number of components per point. An integer from `1` to `15` for `:r2`, `1` to `40` for `:halton` and `1` to `1000` for `:sobol`; an integer `>= 1` for the other generators.
  - `rng` (optional): random number generator (see [[ensure-rng]]) which the pseudo-random generators `:default`, `:gaussian`, `:sphere` and `:ball` draw from. The deterministic `:r2`, `:halton` and `:sobol` ignore it. When it is missing or `nil`, a new `JDKRandomGenerator` is created for every call, the shared [[default-rng]] is never used. The returned sequence holds the generator, so it must be consumed from one thread.

  Returns a lazy, infinite sequence of points: a double when `dimensions` is `1`, a `Vec2`, `Vec3` or `Vec4` when `dimensions` is `2`, `3` or `4`, or a plain vector of doubles otherwise.

  The `:halton` and `:sobol` sequences start with the origin (all coordinates `0.0`), so as uniforms for an inverse cdf the first point gives `-Inf` for an unbounded distribution; `:r2` does not start there. Every call starts the deterministic sequences again from the first point. Apache Commons Math indexes `:halton` and `:sobol` with an `int`, the points after `2^31 - 1` are not valid; this is out of reach in practice.

  Throws `ex-info` (`:type :invalid-input`, `:property :dimensions`, `:expected` the allowed range) when `dimensions` is outside the range above, is not an integer (`nil`, `2.5`) - also when assertions are disabled - and `ex-info` with `{:rng value}` when `rng` is neither a generator nor `nil`. The dimensions are checked when the sequence is created.

  See also [[jittered-sequence-generator]], [[ensure-rng]], [[ball-random]]."}
  sequence-generator (fn [seq-generator _ & _] seq-generator))
(defmethod sequence-generator :halton [seq-generator dimensions & [rng]]
  (rv-generators seq-generator dimensions (ensure-rng rng)))
(defmethod sequence-generator :sobol [seq-generator dimensions & [rng]]
  (rv-generators seq-generator dimensions (ensure-rng rng)))
(defmethod sequence-generator :r2 [seq-generator dimensions & [rng]]
  (rv-generators seq-generator dimensions (ensure-rng rng)))
(defmethod sequence-generator :sphere [seq-generator dimensions & [rng]]
  (rv-generators seq-generator dimensions (ensure-rng rng)))
(defmethod sequence-generator :gaussian [seq-generator dimensions & [rng]]
  (random-generators seq-generator dimensions (ensure-rng rng)))

(defmethod sequence-generator :default [seq-generator dimensions & [rng]]
  (random-generators seq-generator dimensions (ensure-rng rng)))
(defmethod sequence-generator :ball [seq-generator dimensions & [rng]]
  (let [dimensions (check-dimensions 'sequence-generator seq-generator dimensions)
        rng (ensure-rng rng)]
    (repeatedly #(ball-point rng dimensions))))

(defn jittered-sequence-generator
  "Creates a lazy, infinite sequence of jittered [[sequence-generator]] points.

  Perturbs each point of the underlying `seq-generator` sequence by a small amount of noise, according to this [article](http://extremelearning.com.au/a-simple-method-to-construct-isotropic-quasirandom-blue-noise-point-sequences/), breaking up the perfectly regular structure of low-discrepancy sequences while keeping their approximately uniform coverage, closer to a blue-noise point distribution. Three jittering strategies are used depending on `seq-generator`:

  - `:sphere` and `:gaussian` - each coordinate is offset by independent `N(0,1)` noise scaled by `jitter`; the result is not renormalized back onto the sphere or wrapped, so points may drift away from it.
  - `:ball` - each coordinate is offset by scaled quasi-random noise (generic constants) and wrapped periodically into `[-1,1)` with `((x + 1) mod 2) - 1`, so the points keep both signs. The wrapping is per coordinate, so a jittered point is in the cube `[-1,1)` but not necessarily in the unit ball.
  - every other generator (including `:r2`, `:halton`, `:sobol` and `:default`/`:uniform`) - each coordinate is offset by scaled quasi-random noise (tuned per `seq-generator`, falling back to generic constants for unlisted ones) and wrapped back into `[0,1)` with `frac`.

  The constants of the offset are calibrated in the article for the two-dimensional `:r2`, `:halton` and `:sobol` and are used unchanged for every number of dimensions and for the other generators. The offset of a coordinate is one-sided, uniform in `[0, scale)`, with a `scale` which decreases with the index of the point.

  Parameters:

  - `seq-generator`: keyword selecting the base generator, same as for [[sequence-generator]]. Intended for `:r2`, `:sobol` and `:halton`, whose evenly spaced points benefit most from jittering, but works with any registered generator.
  - `dimensions` (integer): number of components per point, passed through to [[sequence-generator]]; the same limits apply.
  - `jitter` (optional, finite number `>= 0`): jitter amount. `0.0` gives exactly the unjittered [[sequence-generator]] sequence for `:r2`, `:halton` and `:sobol` only: for the pseudo-random generators the jitter numbers are still drawn from `rng` between the base points, so the base stream is shifted, and the wrapped `:ball` is not the identity either. `1.0` is full jitter in the sense of the article, larger values are allowed (the points stay in the wrapped range). Default: `0.25`.
  - `rng` (optional): random number generator (see [[ensure-rng]]) which draws the jitter noise and, for the pseudo-random generators, the base points. The deterministic `:r2`, `:halton` and `:sobol` base points do not depend on it. When it is missing or `nil`, new `JDKRandomGenerator`s are created for every call, the shared [[default-rng]] is never used, so seed the generator to reproduce a sequence. The returned sequence holds the generator, so it must be consumed from one thread.

  Returns a lazy, infinite sequence of points, same shape as [[sequence-generator]] for the given `dimensions` (a double, `Vec2`, `Vec3`, `Vec4` or a plain vector).

  Throws `ex-info` (`:type :invalid-input`) when `dimensions` is invalid (`:property :dimensions`), when `jitter` is not a finite number `>= 0` (`NaN`, infinity, a negative number, `nil`, a non-number, for example a generator passed in its place - `:property :jitter`), and `ex-info` with `{:rng value}` when `rng` is neither a generator nor `nil`.

  See also [[sequence-generator]], [[ensure-rng]]."
  ([seq-generator dimensions] (jittered-sequence-generator seq-generator dimensions 0.25 nil))
  ([seq-generator dimensions jitter] (jittered-sequence-generator seq-generator dimensions jitter nil))
  ([seq-generator dimensions jitter rng]
   (let [dimensions (check-dimensions 'jittered-sequence-generator seq-generator dimensions)
         jitter (check-jitter jitter)
         ball? (= :ball seq-generator)
         rng-given? (some? rng)
         rng (ensure-rng rng)
         ;; without an rng the sequence generator creates its own, which keeps methods defined elsewhere with two arguments working
         s (if rng-given?
             (sequence-generator seq-generator dimensions rng)
             (sequence-generator seq-generator dimensions))
         [j mod-fn] (if (#{:sphere :gaussian} seq-generator)
                      (let [j (sequence-generator :gaussian dimensions rng)
                            jitter-low (* m/SQRTPI 0.5 0.25 jitter)]
                        [j (if (m/one? dimensions)
                             (fn [^double v ^double vj] (+ v (* jitter-low vj)))
                             (fn [v vj] (v/add v (v/mult vj jitter-low))))])
                      (let [j (jitter-generator seq-generator dimensions jitter rng)]
                        [j (cond
                             (m/one? dimensions) (if ball?
                                                   (fn [^double v ^double vj] (wrap-ball (+ v vj)))
                                                   (fn [^double v ^double vj] (m/frac (+ v vj))))
                             (m/>= dimensions plain-vector-dimensions) (fn [v vj] (add-jitter v vj ball?))
                             ball? (fn [v vj] (v/fmap (v/add v vj) wrap-ball))
                             :else (fn [v vj] (v/fmap (v/add v vj) m/frac)))]))]
     (map mod-fn s j))))

(def ^{:doc "Sequence of the keys of the generators registered by [[sequence-generator]]: `:halton`, `:sobol`, `:r2`, `:sphere`, `:gaussian`, `:ball` and `:default`, in no particular order.

The sequence is a snapshot taken when the namespace is loaded: generators added later with `defmethod` on [[sequence-generator]] do not appear in it. It contains `:default`, the generator of uniform points, which any unregistered key (`:uniform` included) selects as well, so `:uniform` is not in the list."}
  sequence-generators-list (keys (methods sequence-generator)))

;; ## Noise

(def ^{:doc "The interpolations of the `:interpolation` option of the noise configuration, as a map of the keyword to the integer constant of `fastmath.java.noise.NoiseConfig`.

  Keys: `:none` (the nearest lattice value, discontinuous), `:linear` (C0), `:hermite` (cubic `3t^2 - 2t^3`, C1, the default) and `:quintic` (`6t^5 - 15t^4 + 10t^3`, C2). Only `:value` and `:gradient` noise (see [[noise-types]]) use the interpolation, `:simplex` ignores it. Any other key in the `:interpolation` option throws, see [[fbm-noise]].

  The map is a constant: the keys are the values accepted by the noise functions, the integers are an implementation detail."}
  noise-interpolations {:none NoiseConfig/INTERPOLATE_NONE
                        :linear NoiseConfig/INTERPOLATE_LINEAR
                        :hermite NoiseConfig/INTERPOLATE_HERMITE
                        :quintic NoiseConfig/INTERPOLATE_QUINTIC})

(def ^{:doc "The noise types of the `:noise-type` option of the noise configuration, as a map of the keyword to the integer constant of `fastmath.java.noise.NoiseConfig`.

  Keys: `:value` (interpolated random lattice values), `:gradient` (interpolated dot products of random lattice gradients, the default; not Ken Perlin's improved noise, it has 4 diagonal gradients in 2-D and a different hash) and `:simplex` (Gustavson's simplex noise, ignores the interpolation). Any other key in the `:noise-type` option throws, see [[fbm-noise]].

  The map is a constant: the keys are the values accepted by the noise functions, the integers are an implementation detail."}
  noise-types {:value NoiseConfig/NOISE_VALUE
               :gradient NoiseConfig/NOISE_GRADIENT
               :simplex NoiseConfig/NOISE_SIMPLEX})

(defn- noise-option
  "Returns the option `property` of the configuration map `cfg` when `valid?` accepts it, `default` when the option is missing or `nil`, otherwise throws `ex-info` (`:type :invalid-input`) naming `function`."
  [function cfg property default valid? expected]
  (let [v (get cfg property)]
    (cond
      (nil? v) default
      (valid? v) v
      :else (throw (invalid-input function property v expected
                                  (str "Noise option " property " must be " expected ", got " (pr-str v)))))))

(defn- noise-config
  "Creates the `NoiseConfig` object of the configuration map `cfg` (`nil` = empty), the options are validated, a missing or `nil` option takes its default and unknown keys are ignored. `function` names the public function for the errors. A missing `:seed` is `(irand)` of the shared default generator."
  ^NoiseConfig [function cfg]
  (when-not (or (nil? cfg) (map? cfg))
    (throw (invalid-input function :cfg cfg "a map or nil"
                          (str "The noise configuration must be a map or nil, got " (pr-str cfg)))))
  (let [seed (noise-option function cfg :seed nil integer? "an integer")
        noise-type (noise-option function cfg :noise-type :gradient #(contains? noise-types %)
                                 (str "one of " (pr-str (keys noise-types))))
        interpolation (noise-option function cfg :interpolation :hermite #(contains? noise-interpolations %)
                                    (str "one of " (pr-str (keys noise-interpolations))))
        octaves (noise-option function cfg :octaves 6 #(and (long-integer? %) (<= 1 (long %) 64))
                              "an integer from 1 to 64")
        lacunarity (noise-option function cfg :lacunarity 2.0 #(and (finite-number? %) (> (double %) 0.0))
                                 "a finite number > 0")
        gain (noise-option function cfg :gain 0.5 #(and (finite-number? %) (>= (double %) 0.0))
                           "a finite number >= 0")
        normalize? (noise-option function cfg :normalize? true #(instance? Boolean %) "a boolean")]
    (NoiseConfig. (unchecked-int (long (if (nil? seed) (irand) (.longValue ^Number seed))))
                  (int (noise-types noise-type))
                  (int (noise-interpolations interpolation))
                  (.intValue ^Number octaves)
                  (.doubleValue ^Number lacunarity)
                  (.doubleValue ^Number gain)
                  (boolean normalize?))))

(def ^:private ^:const preset-noise-seed 1337)

(defonce ^:private perlin-noise-config (noise-config 'noise {:interpolation :quintic :seed preset-noise-seed}))
(defonce ^:private simplex-noise-config (noise-config 'simplex {:noise-type :simplex :seed preset-noise-seed}))
(defonce ^:private value-noise-config (noise-config 'vnoise {:noise-type :value :seed preset-noise-seed}))

(defn vnoise
  "Value noise, 6-octave FBM with Hermite interpolation.

  A ready-to-use preset equal to `(fbm-noise {:seed 1337 :noise-type :value})`: `:interpolation` `:hermite`, `:octaves` `6`, `:lacunarity` `2.0`, `:gain` `0.5`, `:normalize?` `true`. The underlying noise interpolates randomly assigned values at integer lattice points, giving a blockier, less directional look than gradient-based [[noise]]. The `:seed` is fixed to `1337` (the same seed as [[noise]] and [[simplex]]), so the values are deterministic: the same arguments give the same value in every call, session and JVM. The seed cannot be changed, build your own function with [[fbm-noise]] for another seed.

  Parameters:

  - `x`, `y`, `z` (double): coordinates, one, two or three of them. The noise of 1, 2 and 3 dimensions are unrelated functions (`(vnoise x)`, `(vnoise x 0.0)` and `(vnoise x 0.0 0.0)` differ), and the 1-D noise repeats with the period 256.

  Returns a double. The result lies in `[0,1]` but the realised range is narrower (the octaves rarely reach their extremes together). It is `NaN` when a coordinate is `NaN`, infinite or has an absolute value above `2^30 / 2^5 = 3.36e7` (the limit of 6 octaves with lacunarity 2, see [[fbm-noise]]).

  See also [[noise]], [[simplex]], [[single-noise]], [[fbm-noise]], [[billow-noise]], [[ridgedmulti-noise]]."
  (^double [^double x] (FBM/noise value-noise-config x))
  (^double [^double x ^double y] (FBM/noise value-noise-config x y))
  (^double [^double x ^double y ^double z] (FBM/noise value-noise-config x y z)))

(defn noise
  "Gradient (Perlin-like) noise, 6-octave FBM with quintic interpolation.

  A ready-to-use preset equal to `(fbm-noise {:seed 1337 :interpolation :quintic})`: `:noise-type` `:gradient`, `:octaves` `6`, `:lacunarity` `2.0`, `:gain` `0.5`, `:normalize?` `true`. Gradient noise interpolates dot products of pseudo-random gradient vectors at lattice points, and the quintic (`6t^5 - 15t^4 + 10t^3`) interpolation removes second-derivative discontinuities at cell boundaries, giving the smooth look of Ken Perlin's improved noise (it is not that algorithm: the hash and the gradient set differ, 2-D uses 4 diagonal gradients). The `:seed` is fixed to `1337` (the same seed as [[vnoise]] and [[simplex]]), so the values are deterministic: the same arguments give the same value in every call, session and JVM. The seed cannot be changed, build your own function with [[fbm-noise]] for another seed.

  Parameters:

  - `x`, `y`, `z` (double): coordinates, one, two or three of them. The noise of 1, 2 and 3 dimensions are unrelated functions (`(noise x)`, `(noise x 0.0)` and `(noise x 0.0 0.0)` differ), and the 1-D noise repeats with the period 256.

  Returns a double. The result lies in `[0,1]` but the realised range is narrower, the 2-D noise is about `[0.16, 0.84]` (the octaves rarely reach their extremes together, and gradient noise is `0.5` at the lattice points). It is `NaN` when a coordinate is `NaN`, infinite or has an absolute value above `2^30 / 2^5 = 3.36e7` (the limit of 6 octaves with lacunarity 2, see [[fbm-noise]]).

  See also [[vnoise]], [[simplex]], [[single-noise]], [[fbm-noise]], [[billow-noise]], [[ridgedmulti-noise]]."
  (^double [^double x] (FBM/noise perlin-noise-config x))
  (^double [^double x ^double y] (FBM/noise perlin-noise-config x y))
  (^double [^double x ^double y ^double z] (FBM/noise perlin-noise-config x y z)))

(defn simplex
  "Simplex noise, 6-octave FBM.

  A ready-to-use preset equal to `(fbm-noise {:seed 1337 :noise-type :simplex})`: `:octaves` `6`, `:lacunarity` `2.0`, `:gain` `0.5`, `:normalize?` `true`. Simplex noise evaluates gradients on a simplectic (triangular/tetrahedral) lattice rather than a square/cubic grid, which reduces directional artifacts and scales better to higher dimensions than gradient noise. The `:seed` is fixed to `1337` (the same seed as [[noise]] and [[vnoise]]), so the values are deterministic: the same arguments give the same value in every call, session and JVM. The seed cannot be changed, build your own function with [[fbm-noise]] for another seed.

  Parameters:

  - `x`, `y`, `z` (double): coordinates, one, two or three of them. The noise of 1, 2 and 3 dimensions are unrelated functions (`(simplex x)`, `(simplex x 0.0)` and `(simplex x 0.0 0.0)` differ), and the 1-D noise repeats with the period 256.

  Returns a double. The result lies in `[0,1]` but the realised range is narrower (the octaves rarely reach their extremes together). It is `NaN` when a coordinate is `NaN`, infinite or has an absolute value above `2^30 / 2^5 = 3.36e7` (the limit of 6 octaves with lacunarity 2, see [[fbm-noise]]).

  See also [[noise]], [[vnoise]], [[single-noise]], [[fbm-noise]], [[billow-noise]], [[ridgedmulti-noise]]."
  (^double [^double x] (FBM/noise simplex-noise-config x))
  (^double [^double x ^double y] (FBM/noise simplex-noise-config x y))
  (^double [^double x ^double y ^double z] (FBM/noise simplex-noise-config x y z)))

(defn single-noise
  "Creates a single-octave noise function.

  Produces raw, unblended noise: one evaluation of the underlying noise type (value, gradient or simplex, selected by `:noise-type`) per point, without combining multiple octaves. It is the building block of [[fbm-noise]], [[billow-noise]] and [[ridgedmulti-noise]], but it is not equal to `(fbm-noise (assoc cfg :octaves 1))`: the first octave of the multi-octave noises uses a different lattice offset derived from the seed.

  Parameters:

  - `cfg` (optional, map or `nil`): noise configuration, all keys optional, a `nil` value stands for a missing key (the default is used), unknown keys are ignored:
      - `:seed` - integer seed, default: a random integer drawn from the shared [[default-rng]] (so `(set-seed! n)` makes it reproducible). The noise uses the low 32 bits of the seed: an integer outside of the 32-bit range is accepted and truncated, `1` and `4294967297` give the same noise.
      - `:noise-type` - a key of [[noise-types]]: `:value`, `:gradient` or `:simplex`, default: `:gradient`.
      - `:interpolation` - a key of [[noise-interpolations]]: `:none`, `:linear`, `:hermite` or `:quintic`, used only by `:value` and `:gradient` noise types, default: `:hermite`.
      - `:normalize?` - a boolean, normalize result to `[0,1]` range (`true`, default) or leave it in `[-1,1]` (`false`).
      - `:octaves`, `:lacunarity`, `:gain` - validated as in [[fbm-noise]] but have no effect on the values, since only a single octave is evaluated. `:octaves` and `:lacunarity` still set the limit of the coordinates below (as for [[fbm-noise]]), with the default `6` and `2.0` it is `3.36e7`.

  Returns a function of one, two or three double arguments (`x`, `[x y]` or `[x y z]`), returning a double. The result lies in `[0,1]` when `:normalize?` is true (the default), or `[-1,1]` otherwise (gradient noise of 3 dimensions can exceed the range by about 0.6%, 1.0061 was observed). The noise of 1, 2 and 3 dimensions are unrelated functions, and the 1-D noise repeats with the period 256. Evaluation never throws: a coordinate which is `NaN`, infinite or has an absolute value above the limit `2^30 / lacunarity^(octaves-1)` (for `:lacunarity` below 1 the limit is `2^30`) gives `NaN`.

  Throws `ex-info` (`:type :invalid-input`, `:property` the option) for an invalid option, see [[fbm-noise]].

  See also [[fbm-noise]], [[billow-noise]], [[ridgedmulti-noise]], [[noise]], [[vnoise]], [[simplex]], [[random-noise-cfg]], [[random-noise]]."
  ([] (single-noise nil))
  ([cfg]
   (let [ncfg (noise-config 'single-noise cfg)]
     (fn
       (^double [^double x] (Noise/noise ncfg x))
       (^double [^double x ^double y] (Noise/noise ncfg x y))
       (^double [^double x ^double y ^double z] (Noise/noise ncfg x y z))))))

(defn fbm-noise
  "Creates a Fractal Brownian Motion (FBM) noise function.

  Sums several octaves of the underlying noise type, each with frequency scaled by `:lacunarity` and amplitude scaled by `:gain` relative to the previous octave, producing natural, self-similar terrain-like noise. The sum is divided by the sum of the amplitudes `1 + gain + ... + gain^(octaves-1)`. [[noise]], [[vnoise]] and [[simplex]] are ready-to-use FBM noise functions with fixed presets.

  Parameters:

  - `cfg` (optional, map or `nil`): noise configuration, all keys optional, a `nil` value stands for a missing key (the default is used), unknown keys are ignored:
      - `:seed` - integer seed, default: a random integer drawn from the shared [[default-rng]] (so `(set-seed! n)` makes it reproducible). The noise uses the low 32 bits of the seed: an integer outside of the 32-bit range is accepted and truncated, `1` and `4294967297` give the same noise.
      - `:noise-type` - a key of [[noise-types]]: `:value`, `:gradient` or `:simplex`, default: `:gradient`.
      - `:interpolation` - a key of [[noise-interpolations]]: `:none`, `:linear`, `:hermite` or `:quintic`, used only by `:value` and `:gradient` noise types, default: `:hermite`.
      - `:octaves` - number of octaves summed together, an integer from `1` to `64`, default: `6`.
      - `:lacunarity` - frequency multiplier applied to each successive octave, a finite number `> 0`, default: `2.0`.
      - `:gain` - amplitude multiplier applied to each successive octave, a finite number `>= 0`, default: `0.5`.
      - `:normalize?` - a boolean, normalize result to `[0,1]` range (`true`, default) or leave it in `[-1,1]` (`false`).

  Returns a function of one, two or three double arguments (`x`, `[x y]` or `[x y z]`), returning a double. The result lies in `[0,1]` when `:normalize?` is true (the default), or `[-1,1]` otherwise, but the realised range is narrower, because the octaves rarely reach their extremes together (the default gradient noise in 2-D, see [[noise]], covers about `[0.16, 0.84]`). The noise of 1, 2 and 3 dimensions are unrelated functions (`(f x)`, `(f x 0.0)` and `(f x 0.0 0.0)` differ), and the 1-D noise repeats with the period 256.

  Evaluation never throws: a coordinate which is `NaN`, infinite or has an absolute value above the limit `2^30 / lacunarity^(octaves-1)` (for `:lacunarity` below 1 the limit is `2^30`) gives `NaN`, because the lattice index of the last octave would overflow. For 6 octaves and lacunarity 2 the limit is `3.36e7`.

  Throws `ex-info` (`:type :invalid-input`, `:property` the option, `:value`, `:expected`) when `cfg` is not a map or `nil` or an option is invalid: an unknown `:noise-type` or `:interpolation`, `:octaves` outside `1..64` or not an integer, `:lacunarity` not a finite number `> 0`, `:gain` not a finite number `>= 0`, `:normalize?` not a boolean, `:seed` not an integer.

  See also [[single-noise]], [[billow-noise]], [[ridgedmulti-noise]], [[noise]], [[vnoise]], [[simplex]], [[random-noise-cfg]], [[random-noise]]."
  ([] (fbm-noise nil))
  ([cfg]
   (let [ncfg (noise-config 'fbm-noise cfg)]
     (fn
       (^double [^double x] (FBM/noise ncfg x))
       (^double [^double x ^double y] (FBM/noise ncfg x y))
       (^double [^double x ^double y ^double z] (FBM/noise ncfg x y z))))))

(defn billow-noise
  "Creates a billow noise function.

  Similar to [[fbm-noise]], sums several octaves of the underlying noise type scaled by `:lacunarity` and `:gain`, but folds each octave's value through `abs(v) * 2 - 1` before summing. This produces puffy, billowy cloud-like patterns instead of the smoother look of plain FBM noise.

  Parameters:

  - `cfg` (optional, map or `nil`): noise configuration with the same keys, defaults and validation as [[fbm-noise]] (`:seed`, `:noise-type`, `:interpolation`, `:octaves` `1..64`, `:lacunarity` `> 0`, `:gain` `>= 0`, `:normalize?`). A `nil` value stands for a missing key, unknown keys are ignored. A missing `:seed` is a random integer drawn from the shared [[default-rng]], the noise uses its low 32 bits.

  Returns a function of one, two or three double arguments (`x`, `[x y]` or `[x y z]`), returning a double. The result lies in `[0,1]` when `:normalize?` is true (the default), or `[-1,1]` otherwise; the realised range is narrower and not centred (the default gradient billow noise has a mean of about `0.2` in 3-D). As for [[fbm-noise]], the noise of different dimensions are unrelated, the 1-D noise repeats with the period 256 and a coordinate which is `NaN`, infinite or above the limit `2^30 / lacunarity^(octaves-1)` in absolute value gives `NaN`.

  Throws `ex-info` (`:type :invalid-input`, `:property` the option) for an invalid option, see [[fbm-noise]].

  See also [[single-noise]], [[fbm-noise]], [[ridgedmulti-noise]], [[random-noise-cfg]], [[random-noise]]."
  ([] (billow-noise nil))
  ([cfg]
   (let [ncfg (noise-config 'billow-noise cfg)]
     (fn
       (^double [^double x] (Billow/noise ncfg x))
       (^double [^double x ^double y] (Billow/noise ncfg x y))
       (^double [^double x ^double y ^double z] (Billow/noise ncfg x y z))))))

(defn ridgedmulti-noise
  "Creates a ridged multifractal noise function.

  Combines several octaves of the underlying noise type into sharp, ridge-like features: each octave's value is folded through `(1 - abs(v))^2`, weighted by the strength of the previous octave's signal divided by `:gain` (and clamped to `[0,1]`), with frequency scaled by `:lacunarity` as usual. This feedback between octaves produces the jagged mountain-range look typical of ridged multifractal noise.

  Note that `:gain` here divides the signal (a larger gain gives a weaker feedback; the default `0.5` corresponds to the libnoise gain `2.0`), and it is not the amplitude ratio of the octaves as in [[fbm-noise]]: the amplitude of octave `i` is `lacunarity^-i`. The sum is divided by `1 + 1/lacunarity + ... + lacunarity^-(octaves-1)`, so the normalised result never exceeds `1`; the maximum `1` is reached at the lattice points of gradient noise when `:gain <= 1`, for a larger gain the weights shrink and the realised maximum is smaller.

  Parameters:

  - `cfg` (optional, map or `nil`): noise configuration with the same keys, defaults and validation as [[fbm-noise]] (`:seed`, `:noise-type`, `:interpolation`, `:octaves` `1..64`, `:lacunarity` `> 0`, `:normalize?`), with `:gain` a finite number `> 0` (default `0.5`; `0` is rejected because the signal is divided by it). A `nil` value stands for a missing key, unknown keys are ignored. A missing `:seed` is a random integer drawn from the shared [[default-rng]], the noise uses its low 32 bits.

  Returns a function of one, two or three double arguments (`x`, `[x y]` or `[x y z]`), returning a double. The result lies in `[0,1]` when `:normalize?` is true (the default), or `[-1,1]` otherwise. As for [[fbm-noise]], the noise of different dimensions are unrelated, the 1-D noise repeats with the period 256 and a coordinate which is `NaN`, infinite or above the limit `2^30 / lacunarity^(octaves-1)` in absolute value gives `NaN`.

  Throws `ex-info` (`:type :invalid-input`, `:property` the option) for an invalid option, see [[fbm-noise]].

  See also [[single-noise]], [[fbm-noise]], [[billow-noise]], [[random-noise-cfg]], [[random-noise]]."
  ([] (ridgedmulti-noise nil))
  ([cfg]
   (let [ncfg (noise-config 'ridgedmulti-noise cfg)]
     (when-not (> (.gain ncfg) 0.0)
       (throw (invalid-input 'ridgedmulti-noise :gain (.gain ncfg) "a finite number > 0"
                             (str "The gain of ridged multifractal noise must be a finite number > 0, got " (.gain ncfg)))))
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

(defn- noise-function?
  "Is `x` a function, or a var holding one?"
  [x]
  (or (fn? x) (and (var? x) (fn? @x))))

(defn- warp-parameters
  "Validates the arguments of [[warp-noise]] and [[warp-noise-fn]], returns the vector `[noise scale depth]` with `scale` a double and `depth` a long, or throws `ex-info` (`:type :invalid-input`) naming `function`."
  [function noise scale depth]
  (when-not (noise-function? noise)
    (throw (invalid-input function :noise noise "a function of 1, 2 or 3 doubles"
                          (str "The noise to warp must be a function, got " (pr-str noise)))))
  (when-not (finite-number? scale)
    (throw (invalid-input function :scale scale "a finite number"
                          (str "The warp scale must be a finite number, got " (pr-str scale)))))
  (when-not (and (long-integer? depth) (>= (long depth) 0))
    (throw (invalid-input function :depth depth "an integer >= 0"
                          (str "The warp depth must be an integer >= 0, got " (pr-str depth)))))
  [noise (double scale) (long depth)])

(defn- make-warp
  [function noise scale depth]
  (let [[noise scale depth] (warp-parameters function noise scale depth)
        n1 (make-warp-1d noise scale depth)
        n2 (make-warp-2d noise scale depth)
        n3 (make-warp-3d noise scale depth)]
    (fn
      (^double [^double x] (n1 x))
      (^double [^double x ^double y] (n2 x y))
      (^double [^double x ^double y ^double z] (n3 x y z)))))

(defn warp-noise-fn
  "Deprecated, use [[warp-noise]].

  Creates a domain-warped noise function; it is the same function as [[warp-noise]] (including the validation of the arguments) and will be removed in a future version.

  Applies domain warping as described in [Inigo Quilez's article](http://www.iquilezles.org/www/articles/warp/warp.htm): at each of `depth` warp levels, the input coordinates are recursively perturbed by evaluating the warp itself (with fixed constant offsets) at the previous level, scaled by `scale` and added to the coordinates, before evaluating `noise` at the warped position. A `depth` of `0` performs no warping and simply evaluates `noise` directly.

  Parameters:

  - `noise` (optional, function): the noise function to warp; must accept 1, 2 or 3 double arguments, same shape as [[vnoise]], [[noise]] or [[simplex]]. Defaults to [[vnoise]].
  - `scale` (optional, finite number): strength of the coordinate perturbation added at each warp level. Defaults to `4.0`.
  - `depth` (optional, integer `>= 0`): number of recursive warp levels. Defaults to `1`.

  Returns a function of one, two or three double arguments (`x`, `[x y]` or `[x y z]`), returning a double. The normalization (value range) of the returned function follows the normalization of `noise`.

  Throws `ex-info` (`:type :invalid-input`, `:property` `:noise`, `:scale` or `:depth`) at construction when `noise` is not a function, `scale` is not a finite number or `depth` is not an integer `>= 0`.

  See also [[warp-noise]], [[vnoise]], [[noise]], [[simplex]], [[random-noise]]."
  {:deprecated "Use `warp-noise` instead."}
  ([noise scale depth] (make-warp 'warp-noise-fn noise scale depth))
  ([noise scale] (make-warp 'warp-noise-fn noise scale 1))
  ([noise] (make-warp 'warp-noise-fn noise 4.0 1))
  ([] (make-warp 'warp-noise-fn vnoise 4.0 1)))

(defn warp-noise
  "Creates a domain-warped noise function.

  Applies domain warping as described in [Inigo Quilez's article](http://www.iquilezles.org/www/articles/warp/warp.htm): at each of `depth` warp levels, the input coordinates are recursively perturbed by evaluating the warp itself (with fixed constant offsets) at the previous level, scaled by `scale` and added to the coordinates, before evaluating `noise` at the warped position. A `depth` of `0` performs no warping and simply evaluates `noise` directly. At `depth` 1 this is Quilez's `f(p + scale * q)`, deeper levels evaluate the warp vector itself with the warped noise at offset points.

  Parameters:

  - `noise` (optional, function): the noise function to warp; must accept 1, 2 or 3 double arguments, same shape as [[vnoise]], [[noise]] or [[simplex]] (a function which accepts fewer arities fails when called with the missing one, not at construction). Defaults to [[vnoise]].
  - `scale` (optional, finite number): strength of the coordinate perturbation added at each warp level. Defaults to `4.0`.
  - `depth` (optional, integer `>= 0`): number of recursive warp levels; each additional level multiplies the number of `noise` evaluations (`depth + 1` in 1-D, `2^(depth+1) - 1` in 2-D, `(3^(depth+1) - 1) / 2` in 3-D), so larger values are increasingly expensive. Defaults to `1`.

  Returns a function of one, two or three double arguments (`x`, `[x y]` or `[x y z]`), returning a double. The normalization (value range) of the returned function follows the normalization of `noise`; the warped coordinates are passed to `noise`, so the limit on the coordinates of `noise` (`NaN` above it) applies to them.

  Throws `ex-info` (`:type :invalid-input`, `:property` `:noise`, `:scale` or `:depth`) at construction when `noise` is not a function, `scale` is not a finite number or `depth` is not an integer `>= 0` (a negative depth used to overflow the stack at the first call).

  See also [[vnoise]], [[noise]], [[simplex]], [[random-noise]]."
  ([noise scale depth] (make-warp 'warp-noise noise scale depth))
  ([noise scale] (make-warp 'warp-noise noise scale 1))
  ([noise] (make-warp 'warp-noise noise 4.0 1))
  ([] (make-warp 'warp-noise vnoise 4.0 1)))

(defonce ^{:doc "The noise generators of the `:generator` option of [[random-noise-cfg]] and [[random-noise]], as a map of the key to the function which creates the noise from a configuration map.

  Keys: `:fbm` ([[fbm-noise]]), `:single` ([[single-noise]]), `:billow` ([[billow-noise]]) and `:ridgemulti` ([[ridgedmulti-noise]]; note that the key is `:ridgemulti`, not the name of the function; [[random-noise-cfg]] and [[random-noise]] also accept `:ridgedmulti` as its alias, which is not a key of this map)."}
  noise-generators
  {:fbm fbm-noise
   :single single-noise
   :billow billow-noise
   :ridgemulti ridgedmulti-noise})

(defn- noise-generator-fn
  "Returns the noise creating function for the `:generator` option (a key of [[noise-generators]] or `:ridgedmulti`) or throws `ex-info` (`:type :invalid-input`, `:property :generator`) naming `function`."
  [function generator]
  (or (get noise-generators generator)
      (when (= generator :ridgedmulti) ridgedmulti-noise)
      (throw (invalid-input function :generator generator
                            "one of :single, :fbm, :billow, :ridgemulti (alias :ridgedmulti) or nil"
                            (str "Unknown noise generator " (pr-str generator)
                                 ", expected one of :single, :fbm, :billow, :ridgemulti (alias :ridgedmulti) or nil")))))

(defn random-noise-cfg
  "Generates a randomized noise configuration map.

  Produces a configuration suitable for [[random-noise]] (or, for its `:seed`, `:noise-type`, `:interpolation`, `:octaves`, `:lacunarity`, `:gain` and `:normalize?` keys, for [[fbm-noise]], [[single-noise]], [[billow-noise]] and [[ridgedmulti-noise]] directly), with every key assigned a random but sensible value, so calling it repeatedly yields varied noise behavior without manual tuning.

  All the randomness comes from one generator, given by the `:rng` option: the configuration is reproducible with a seeded generator (equal generators give equal maps) and the shared [[default-rng]] is never used or advanced. The draws are always made for every key, in a fixed order, also for the keys fixed by `pre-config`, so a fixed key does not shift the other keys.

  Parameters:

  - `pre-config` (optional, map or `nil`): fixed values to use instead of randomizing; any key present here with a non-`nil` value overrides the corresponding randomly generated value (`nil` stands for a missing key), other keys are copied to the result. Defaults to an empty map, meaning every key is randomized. The key `:rng` is the random number generator (see [[rng]], checked by [[ensure-rng]]); a missing or `nil` `:rng` stands for a new unseeded generator.

  The randomized keys are:

  - `:seed` - a random integer from the range of `int` (the noise uses 32 bits of the seed).
  - `:generator` - one of `:single`, `:fbm`, `:billow`, `:ridgemulti` (each with probability 1/4), used by [[random-noise]] to pick the blending method.
  - `:noise-type` - one of `:value`, `:gradient`, `:simplex` (each with probability 1/3).
  - `:interpolation` - one of `:none`, `:linear`, `:hermite`, `:quintic` (each with probability 1/4).
  - `:octaves` - integer between `1` and `9`, uniform.
  - `:lacunarity` - double between `1.5` and `2.5`, uniform.
  - `:gain` - double between `0.2` and `0.8`, uniform.
  - `:warp-scale` - `0.0` (no warp) with 80% probability, otherwise either `4.0` (50%) or a random double between `0.1` and `10.0`.
  - `:warp-depth` - `1` with 80% probability, otherwise a random integer between `1` and `3` (so `1` has the total probability of about 87%). It is used only when `:warp-scale` is positive.
  - `:normalize?` - always `true`.

  Returns a configuration map. It does not contain `:rng`.

  Throws `ex-info` (`:type :invalid-input`) when `pre-config` is not a map or `nil` (`:property :pre-config`) or its `:generator` is not one of the four above, `:ridgedmulti` (an alias of `:ridgemulti`) or `nil` (`:property :generator`), and `ex-info` with `{:rng value}` when `:rng` is neither a generator nor `nil`. The other keys are validated when the noise is created.

  See also [[random-noise]], [[fbm-noise]], [[single-noise]], [[billow-noise]], [[ridgedmulti-noise]], [[warp-noise]], [[ensure-rng]]."
  ([pre-config]
   (when-not (or (nil? pre-config) (map? pre-config))
     (throw (invalid-input 'random-noise-cfg :pre-config pre-config "a map or nil"
                           (str "The pre-configuration must be a map or nil, got " (pr-str pre-config)))))
   (let [g (ensure-rng (get pre-config :rng))
         seed (irandom g)
         generator (nth [:single :fbm :billow :ridgemulti] (irandom g 4))
         noise-type (nth [:value :gradient :simplex] (irandom g 3))
         interpolation (nth [:none :linear :hermite :quintic] (irandom g 4))
         octaves (irandom g 1 10)
         lacunarity (drandom g 1.5 2.5)
         gain (drandom g 0.2 0.8)
         warp-scale (if (brandom g 0.8) 0.0 (if (brandom g 0.5) 4.0 (drandom g 0.1 10.0)))
         warp-depth (if (brandom g 0.8) 1 (irandom g 1 4))
         given (into {} (remove (fn [e] (nil? (val e)))) (dissoc pre-config :rng))]
     (when (contains? given :generator) (noise-generator-fn 'random-noise-cfg (:generator given)))
     (merge {:seed seed
             :generator generator
             :noise-type noise-type
             :interpolation interpolation
             :octaves octaves
             :lacunarity lacunarity
             :gain gain
             :warp-scale warp-scale
             :warp-depth warp-depth
             :normalize? true}
            given)))
  ([] (random-noise-cfg nil)))

(defn- make-random-noise
  [pre-config]
  (let [cfg (random-noise-cfg pre-config)
        gen-fn (noise-generator-fn 'random-noise (:generator cfg))
        noise (gen-fn cfg)
        warp-scale (:warp-scale cfg)]
    (when-not (finite-number? warp-scale)
      (throw (invalid-input 'random-noise :warp-scale warp-scale "a finite number"
                            (str "The warp scale must be a finite number, got " (pr-str warp-scale)))))
    (if (pos? (double warp-scale))
      (warp-noise noise warp-scale (:warp-depth cfg))
      noise)))

(defn random-noise-fn
  "Deprecated, use [[random-noise]].

  Generates a fully random noise function; it is the same function as [[random-noise]] and will be removed in a future version.

  Combines [[random-noise-cfg]] with one of the noise blending methods (`:single`, `:fbm`, `:billow`, `:ridgemulti`, see [[noise-generators]]) and, when the resulting configuration requests warping (`:warp-scale` greater than `0.0`), wraps the noise with [[warp-noise]]. The result is a ready-to-use noise function in the same shape as [[noise]], [[vnoise]] and [[simplex]].

  Parameters:

  - `cfg` (optional, map or `nil`): configuration overrides passed to [[random-noise-cfg]], including the `:rng` option; any key not provided is filled in randomly. Defaults to `nil`, meaning a fully random configuration.

  Returns a function of one, two or three double arguments (`x`, `[x y]` or `[x y z]`) returning a double.

  Throws `ex-info` for an invalid configuration, see [[random-noise]].

  See also [[random-noise]], [[random-noise-cfg]], [[warp-noise]], [[fbm-noise]]."
  {:deprecated "Use `random-noise` instead."}
  ([cfg] (make-random-noise cfg))
  ([] (make-random-noise nil)))

(defn random-noise
  "Generates a fully random noise function.

  Combines [[random-noise-cfg]] with one of the noise blending methods (`:single`, `:fbm`, `:billow`, `:ridgemulti`, see [[noise-generators]]) and, when the resulting configuration requests warping (`:warp-scale` greater than `0.0`), wraps the noise with [[warp-noise]]. The result is a ready-to-use noise function in the same shape as [[noise]], [[vnoise]] and [[simplex]].

  Parameters:

  - `cfg` (optional, map or `nil`): configuration overrides passed to [[random-noise-cfg]]; any key not provided (or `nil`) is filled in randomly, so by default a random blending method is used (`:generator` is randomised when absent: one of `:single`, `:fbm`, `:billow`, `:ridgemulti`, which is not the `:fbm` default of the other functions) and the warp is random too (`:warp-scale` is `0.0` with 80% probability, `:warp-depth` is used only for a positive `:warp-scale` and is always randomised when absent). The `:generator` may be `:ridgedmulti`, an alias of `:ridgemulti`. The option `:rng` is the random number generator for all the random choices (see [[random-noise-cfg]]): with a seeded generator the same noise function is created every time, without it a new unseeded generator is used and the shared [[default-rng]] is not touched. Defaults to `nil`, meaning a fully random configuration.

  Returns a function of one, two or three double arguments (`x`, `[x y]` or `[x y z]`) returning a double. The value range is `[0,1]` (the realised range is narrower, see [[fbm-noise]]), `NaN` above the limit of the coordinates described in [[fbm-noise]].

  Throws `ex-info` (`:type :invalid-input`, `:property` the option) for an unknown `:generator`, an invalid noise option (see [[fbm-noise]]; a `:gain` of `0` with `:ridgemulti`), a `:warp-scale` which is not a finite number or an invalid `:warp-depth` (see [[warp-noise]]), and `ex-info` with `{:rng value}` when `:rng` is not a generator.

  See also [[random-noise-cfg]], [[warp-noise]], [[fbm-noise]], [[single-noise]], [[billow-noise]], [[ridgedmulti-noise]]."
  ([cfg] (make-random-noise cfg))
  ([] (make-random-noise nil)))

;; ### Discrete noise

(defn discrete-noise
  "Computes a deterministic hash-based noise value for one or two integer coordinates.

  Unlike stochastic random number generators, `discrete-noise` is a pure function: calling it repeatedly with the same `X` (and `Y`) always returns the same value, which makes it useful as a coordinate-based building block for procedural noise and hashing. It is not related to the lattice noise of [[fbm-noise]] and has no seed.

  Parameters:

  - `X` (integer): first coordinate. The hash works on a Java `int`: a `long` argument is truncated to its low 32 bits (`2^32 + 1` gives the same value as `1`, `2^31` the same as `-2^31`), a double is truncated toward zero.
  - `Y` (integer, optional): second coordinate, converted as `X`. Defaults to `0`, producing a 1d hash of `X` alone (`(discrete-noise x)` equals `(discrete-noise x 0)`).

  Returns a double in the `[0,1]` range (closed in theory, the extreme values `0.0` and `1.0` have not been observed). The 2-D hash is computed from `X + 57 * Y`, so different points collide: `(x, y)` and `(x + 57, y - 1)` give the same value, and a 1000 x 1000 grid has only about 58 000 distinct values. Use it for a cheap, deterministic hash of 1-D integer indices; it is not a good 2-D hash.

  See also [[vnoise]], [[simplex]], [[random-noise]]."
  (^double [^long X ^long Y] (Discrete/value X Y))
  (^double [^long X] (Discrete/value X 0)))

;; Distributions

(defmulti
  ^{:doc "Create distribution object.

* First parameter is distribution as a `:key`, one of [[distributions-list]].
* Second parameter is a map with configuration, optional: `nil` or a missing map gives the default parameters. Unknown keys of the map are ignored.

All distributions except `:constant` accept the random number generator under the `:rng` key. The default is a freshly created, unseeded generator for every object: the constructors never use or advance [[default-rng]]. `:constant` draws no random numbers and ignores `:rng`. The Commons Math backed distributions (`:beta`, `:cauchy`, `:chi-squared`, `:exponential`, `:f`, `:gamma`, `:log-normal`, `:nakagami`, `:normal`, `:pareto`, `:t`, `:weibull`) accept the accuracy of the numerical inversion of the cdf under `:inverse-abs-accuracy` (`1.0e-9` by default); a key with another name, like `inverse-cumm-accuracy`, is silently ignored. The accuracy is ABSOLUTE, in the units of the variable, not relative: a quantile is found only up to about `1e-9`, so for a distribution with a scale (or mass) below that the quantiles of the keys without a closed form quantile are quantised and the samplers which invert the cdf give wrong draws. Examples: [[gamma]] `{:shape 2.0 :scale 1e-12}` has `icdf(0.3)` `1.46e-12` instead of `1.10e-12`, [[log-normal]] `{:scale -30.0}` gives `0.0` instead of `5.5e-14`, and the median of 4001 draws of [[chi-squared]] with `0.01` degrees of freedom is `1.1e-9` where it is `1.2e-20` with `:inverse-abs-accuracy` `1e-20` (the same for [[f]] with `0.01` numerator degrees of freedom): the mass below `1e-9` is sampled wrongly. Pass a smaller `:inverse-abs-accuracy` or rescale the variable; `:normal`, `:cauchy`, `:exponential` and `:weibull` have a closed form quantile and ignore it.

An unknown key (also `nil` or a string) throws `IllegalArgumentException` (`No method in multimethod 'distribution' for dispatch value`). Every key has also its own function of the same name which has a complete docstring; `(distribution :normal {:mu 1.0})` is `(normal {:mu 1.0})`.

Not every key supports every protocol function; see the individual functions ([[cdf]], [[ccdf]], [[icdf]], [[mean]], [[variance]], [[lower-bound]], [[means]] ...) for the details. In short: the multivariate keys (`:multi-normal`, `:dirichlet`, `:multinomial`) have no `icdf`, no `mean`, `variance` and bounds (they throw) and no `cdf` and `ccdf`, except for `:multinomial`, whose `cdf` is the joint cumulative probability `P(X_1 <= x_1, ..., X_k <= x_k)` (`ccdf` is `1 - cdf`), but they have `means` and `covariance`; the univariate keys have no `means` and `covariance`; `:categorical` has no bounds and gives `NaN` as `mean` and `variance`; `:anderson-darling`, `:anderson-darling-quick`, `:kolmogorov-smirnov`, `:kolmogorov-smirnov+`, `:kolmogorov-smirnov-quick` and `:watson-g` have no `mean` and `variance`; `set-seed` is not supported by any distribution (use [[set-seed!]]).

Discrete distributions on the integers and a non-integer argument: the Apache Commons Math keys (`:poisson`, `:binomial`, `:bernoulli`, `:geometric`, `:pascal`, `:hypergeometric`, `:zipf`, `:uniform-int`) take the integer part of the argument (`pdf` at `2.5` is `pdf` at `2`); the keys implemented in this library (`:negative-binomial`, `:nbi`, `:nbii`, `:logarithmic`, `:beta-binomial`, `:bb`, `:zero-inflated-*`, `:zero-adjusted-*`, `:integer-discrete`, the noncentral hypergeometric ones) return `0.0` for `pdf` at a non-integer; `cdf` floors the argument for both kinds (`cdf` at `1.999` is `cdf` at `1`), `NaN` gives `NaN`. [[uniform-int]] with the default `upper` has a negative `pdf` (see there).

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
* `:mixture` ([[mixture]]) - `distrs`, `weights`, `rng` - finite weighted mixture combining several component distributions."
   :arglists '([key] [key options])}
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
  - `inverse-abs-accuracy` (double): ABSOLUTE accuracy (in the units of the variable, not relative) of the numerical inversion of the cumulative distribution function done by `icdf` and by the samplers based on it; see [[distribution]] for the consequences. Default: `1.0e-9`.
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
  - `inverse-abs-accuracy` (double): ABSOLUTE accuracy (in the units of the variable, not relative) of the numerical inversion of the cumulative distribution function done by `icdf` and by the samplers based on it; see [[distribution]] for the consequences. Default: `1.0e-9`.
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
  - `inverse-abs-accuracy` (double): ABSOLUTE accuracy (in the units of the variable, not relative) of the numerical inversion of the cumulative distribution function done by `icdf` and by the samplers based on it; see [[distribution]] for the consequences. Default: `1.0e-9`.
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
  - `inverse-abs-accuracy` (double): ABSOLUTE accuracy (in the units of the variable, not relative) of the numerical inversion of the cumulative distribution function done by `icdf` and by the samplers based on it; see [[distribution]] for the consequences. Default: `1.0e-9`.
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
  - `inverse-abs-accuracy` (double): ABSOLUTE accuracy (in the units of the variable, not relative) of the numerical inversion of the cumulative distribution function done by `icdf` and by the samplers based on it; see [[distribution]] for the consequences. Default: `1.0e-9`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  `pdf` and `lpdf` are `0.0` and `##-Inf` for `x < 0`, and at `0` they take the limit (`##Inf` for numerator degrees of freedom below 2, `1.0` for 2, `0.0` above).

  Returns an `FDistribution` object (the subclass `fastmath.java.dist.FFixed`, Apache Commons Math 3.6.1 returns `NaN` for `x <= 0`) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[chi-squared]], [[t]]."
  (^FDistribution [] (f nil))
  (^FDistribution [{:keys [^double numerator-degrees-of-freedom ^double denominator-degrees-of-freedom ^double inverse-abs-accuracy rng]
                    :or {numerator-degrees-of-freedom 1.0 denominator-degrees-of-freedom 1.0
                         inverse-abs-accuracy FDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY}}]
   (FFixed. (or rng (JDKRandomGenerator.)) numerator-degrees-of-freedom denominator-degrees-of-freedom inverse-abs-accuracy)))

(add-distr-method f)

(defn gamma
  "Creates a gamma distribution object.

  The gamma distribution is a continuous distribution over positive reals, shaped by a `shape` parameter and a `scale` parameter. It generalizes the [[exponential]] and [[chi-squared]] distributions and is commonly used to model waiting times, sums of exponential variables, and as a conjugate prior in Bayesian statistics.

  Parameters (single, optional map):

  - `shape` (double): shape parameter, controlling the form of the distribution. Default: `2.0`.
  - `scale` (double): scale parameter, stretching or shrinking the distribution. Default: `2.0`.
  - `inverse-abs-accuracy` (double): ABSOLUTE accuracy (in the units of the variable, not relative) of the numerical inversion of the cumulative distribution function done by `icdf` and by the samplers based on it; see [[distribution]] for the consequences. Default: `1.0e-9`.
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

  Returns a `GumbelDistribution` object (the subclass `fastmath.java.dist.GumbelFixed`, whose `mean` is the correct `mu + beta * 0.5772156649015329`, Apache Commons Math 3.6.1 gives a wrong value) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[levy]], [[laplace]]."
  (^GumbelDistribution [] (gumbel nil))
  (^GumbelDistribution [{:keys [^double mu ^double beta rng]
                         :or {mu 1.0 beta 2.0}}]
   (GumbelFixed. (or rng (JDKRandomGenerator.)) mu beta)))

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

  `pdf`, `lpdf` and `cdf` are `0.0`, `##-Inf` and `0.0` for `x <= mu`; `mean` and `variance` are `##Inf`.

  Returns a `LevyDistribution` object (the subclass `fastmath.java.dist.LevyFixed`, Apache Commons Math 3.6.1 returns `NaN` below `mu`) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[cauchy]], [[gumbel]]."
  (^LevyDistribution [] (levy nil))
  (^LevyDistribution [{:keys [^double mu ^double c rng]
                       :or {mu 0.0 c 1.0}}]
   (LevyFixed. (or rng (JDKRandomGenerator.)) mu c)))

(add-distr-method levy)

(defn logistic
  "Creates a logistic distribution object.

  The logistic distribution is a continuous, symmetric distribution whose cumulative distribution function is the logistic (sigmoid) function. It resembles the [[normal]] distribution in shape but has heavier tails, and underlies logistic regression as well as some growth models.

  Parameters (single, optional map):

  - `mu` (double): location parameter, the center and mean of the distribution. Default: `0.0`.
  - `s` (double): scale parameter, controlling the spread of the distribution. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `LogisticDistribution` object (the subclass `fastmath.java.dist.LogisticFixed`, whose `variance` is the correct `pi^2 s^2 / 3`, Apache Commons Math 3.6.1 gives a wrong value) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[normal]], [[laplace]]."
  (^LogisticDistribution [] (logistic nil))
  (^LogisticDistribution [{:keys [^double mu ^double s rng]
                           :or {mu 0.0 s 1.0}}]
   (LogisticFixed. (or rng (JDKRandomGenerator.)) mu s)))

(add-distr-method logistic)

(defn log-normal
  "Creates a log-normal distribution object.

  The log-normal distribution is a continuous distribution over positive reals whose logarithm follows a [[normal]] distribution. It is commonly used to model quantities that result from the multiplicative combination of many independent positive factors, such as incomes, stock prices or particle sizes.

  Parameters (single, optional map):

  - `scale` (double): location parameter of the underlying normal distribution of the logarithm of the variable (often denoted `mu`). Default: `1.0`.
  - `shape` (double): scale parameter of the underlying normal distribution of the logarithm of the variable (often denoted `sigma`); controls the spread. Default: `1.0`.
  - `inverse-abs-accuracy` (double): ABSOLUTE accuracy (in the units of the variable, not relative) of the numerical inversion of the cumulative distribution function done by `icdf` and by the samplers based on it; see [[distribution]] for the consequences. Default: `1.0e-9`.
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
  - `inverse-abs-accuracy` (double): ABSOLUTE accuracy (in the units of the variable, not relative) of the numerical inversion of the cumulative distribution function done by `icdf` and by the samplers based on it; see [[distribution]] for the consequences. Default: `1.0e-9`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  `cdf` is `0.0` for `x <= 0`. Calling the `cumulativeProbability` method of the object directly at `##Inf` throws a Commons `ConvergenceException`, [[cdf]] returns `1.0` there.

  Returns a `NakagamiDistribution` object (the subclass `fastmath.java.dist.NakagamiFixed`, Apache Commons Math 3.6.1 returns a positive `cdf` for `x < 0`) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[weibull]], [[gamma]]."
  (^NakagamiDistribution [] (nakagami nil))
  (^NakagamiDistribution [{:keys [^double mu ^double omega ^double inverse-abs-accuracy rng]
                           :or {mu 1.0 omega 1.0 inverse-abs-accuracy NakagamiDistribution/DEFAULT_INVERSE_ABSOLUTE_ACCURACY}}]
   (NakagamiFixed. (or rng (JDKRandomGenerator.)) mu omega inverse-abs-accuracy)))

(add-distr-method nakagami)

(defn normal
  "Creates a normal (Gaussian) distribution object.

  The normal distribution is a continuous, symmetric, bell-shaped distribution fully described by its mean `mu` and standard deviation `sd`. It is the most widely used distribution in statistics, arising naturally as the limit of sums of many independent random variables via the central limit theorem.

  Parameters (single, optional map):

  - `mu` (double): mean of the distribution, its center of symmetry. Default: `0.0`.
  - `sd` (double): standard deviation, controlling the spread of the distribution. Default: `1.0`.
  - `inverse-abs-accuracy` (double): ABSOLUTE accuracy (in the units of the variable, not relative) of the numerical inversion of the cumulative distribution function done by `icdf` and by the samplers based on it; see [[distribution]] for the consequences. Default: `1.0e-9`.
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
  - `inverse-abs-accuracy` (double): ABSOLUTE accuracy (in the units of the variable, not relative) of the numerical inversion of the cumulative distribution function done by `icdf` and by the samplers based on it; see [[distribution]] for the consequences. Default: `1.0e-9`.
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
  - `inverse-abs-accuracy` (double): ABSOLUTE accuracy (in the units of the variable, not relative) of the numerical inversion of the cumulative distribution function done by `icdf` and by the samplers based on it; see [[distribution]] for the consequences. Default: `1.0e-9`.
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
  - `inverse-abs-accuracy` (double): ABSOLUTE accuracy (in the units of the variable, not relative) of the numerical inversion of the cumulative distribution function done by `icdf` and by the samplers based on it; see [[distribution]] for the consequences. Default: `1.0e-9`.
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

  The empirical distribution is the Apache Commons Math `EmpiricalDistribution`: the range of the `data` is divided into `bin-count` equal bins and the observations of every bin are summarised by their mean and standard deviation. Sampling picks a bin according to its observed frequency and then draws from a Gaussian kernel with the mean and the standard deviation of the observations in that bin (not uniformly within the bin; a bin with a single observation gives that observation). The samples lie within [[lower-bound]] and [[upper-bound]] (the minimum and the maximum of `data`). It is a rough way of approximating an unknown continuous distribution directly from observations; for a smooth density use [[kde]] or [[continuous-distribution]].

  Parameters (single, optional map):

  - `data` (sequence of numbers): sample used to build the histogram, a non-empty sequence of finite numbers. Default: `[1.0]`.
  - `bin-count` (integer): number of histogram bins, an integer from `1` to `2^31 - 1` (an integral double like `50.0` is accepted, `4.7` is not). Default: `nil`, estimated automatically from `data` using the Freedman-Diaconis rule (50 bins for 5000 standard normal observations).
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values (the single observation `1.0`: `pdf` at `1.0` is `##Inf`, `variance` is `0.0`).

  Throws `ex-info` (`:type :invalid-input`, `:property :data` or `:bin-count`) for `data` which is empty (`[]`, `nil`), is not a sequence (a number, a string, a map) or has a non-finite or non-number element, and for a `bin-count` which is not an integer `>= 1` within the `int` range (`0`, `-1`, `4.7`, `NaN`, `3000000000`, a keyword).

  Returns an `EmpiricalDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions, with the following properties of the Commons implementation:

  - `pdf` is not a histogram density: it is the density of the Gaussian kernel of the bin multiplied by the bin frequency, `NaN` for a point in a bin with fewer than two observations (8% of the points of the range for 5000 standard normal observations in 50 bins; with `data` `[1 1 1 1 9]` and 4 bins `pdf` at 2, 3 and 5 is `NaN`) and `##Inf` at the minimum of the data; hence it does not integrate to 1.
  - `cdf` is not `0` at [[lower-bound]] (`2e-4`, one observation of 5000, for the example above) and differs from the empirical cdf of the data by up to `7e-3` there. [[icdf]] at `0` and `1` are the lower and the upper bound.
  - [[mean]] and [[variance]] are the sample mean and the sample variance (with the denominator `n - 1`) of `data`.

  See also [[distribution]], [[kde]], [[continuous-distribution]], [[enumerated-real]]."
  (^EmpiricalDistribution [] (empirical nil))
  (^EmpiricalDistribution [{:keys [bin-count data rng]
                            :or {data [1.0]}}]
   (check-parameter! 'empirical :data data
                     (and (some? data) (not (number? data)) (not (map? data)) (not (string? data)) (seqable? data)
                          (seq data)
                          (every? finite-number? data))
                     "a non-empty sequence of finite numbers")
   (let [^doubles data (m/seq->double-array data)
         bin-count (if (nil? bin-count)
                     (bins/freedman-diaconis data (alength data))
                     (int-parameter 'empirical :bin-count bin-count 1))
         _ (check-parameter! 'empirical :bin-count bin-count (m/>= bin-count 1) "an integer >= 1")
         ^RandomGenerator rng (or rng (JDKRandomGenerator.))
         ^EmpiricalDistribution d (EmpiricalDistribution. (int bin-count) rng)]
     (.load d data)
     d)))

(add-distr-method empirical)

(defn- enumerated-input
  "Validates the `data` and `probabilities` of the finite support distribution `function`
  and returns `[data probabilities]` as vectors (`probabilities` is `nil` when it is missing).

  `data` has to be a non-empty sequence of numbers, other than `NaN` for `:real` or integers within the `int` range for `:integer`
  (integral doubles are accepted). `probabilities`, when given, has to be a sequence with the count of `data` of finite numbers
  `>= 0` with a finite positive sum. Throws `ex-info` (`:type :invalid-input`) otherwise."
  [function data probabilities data-type]
  (let [sequence-like? (fn [x] (and (some? x) (not (number? x)) (not (map? x)) (not (string? x)) (seqable? x)))
        data' (when (sequence-like? data) (vec data))]
    (check-parameter! function :data data (seq data') "a non-empty sequence")
    (check-parameter! function :data data
                      (every? (fn [x] (and (number? x)
                                           (case data-type
                                             :integer (and (finite-number? x)
                                                           (m/== (double x) (m/floor (double x)))
                                                           (m/<= (double Integer/MIN_VALUE) (double x) (double Integer/MAX_VALUE)))
                                             (not (Double/isNaN (double x))))))
                              data')
                      (if (= data-type :integer) "integers within the int range [-2^31, 2^31-1]" "numbers other than NaN"))
    (let [probs' (when (some? probabilities)
                   (check-parameter! function :probabilities probabilities (sequence-like? probabilities) "a sequence of numbers or nil")
                   (vec probabilities))]
      (when probs'
        (check-parameter! function :probabilities probabilities (== (count probs') (count data'))
                          (str "a sequence with the same count as data (" (count data') ")"))
        (check-parameter! function :probabilities probabilities (every? (fn [p] (and (finite-number? p) (m/not-neg? (double p)))) probs')
                          "finite numbers >= 0")
        (check-parameter! function :probabilities probabilities (let [s (double (reduce m/+ 0.0 (map double probs')))]
                                                                  (and (m/pos? s) (Double/isFinite s)))
                          "numbers with a finite, positive sum"))
      [data' probs'])))

(defn enumerated-real
  "Creates an enumerated real distribution object over an explicit, finite set of double values.

  The distribution puts all of its probability mass on the values listed in `data`. Each value's probability is given by the corresponding entry in `probabilities`, or, when `probabilities` is not supplied, all values are treated as equally likely. Repeated values in `data` accumulate probability mass.

  Parameters (single, optional map):

  - `data` (sequence of numbers): the finite set of values the distribution can take, in any order, `NaN` is not allowed. Default: `[1.0]`.
  - `probabilities` (sequence of doubles): probability associated with each corresponding value in `data`, finite and not negative, with a positive sum; does not need to be normalized. Default: `nil`, meaning uniform probabilities.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  The (value, probability) pairs are sorted ascending by value (a stable sort, a probability stays with its value, repeated values are not merged) before the Apache Commons object is created, because its `icdf` assumes the sorted data and returns wrong quantiles otherwise. `data` `[3 1 2 2 5]` with `probabilities` `[0.1 0.2 0.3 0.1 0.3]` has `icdf` `1 1 1 2 2 3 5 5 5` at `p` 0, 0.1, 0.2, 0.21, 0.6, 0.61, 0.8, 0.81, 1. `pdf`, `mean` and `variance` do not depend on the order (`cdf` may differ in the last digit, the probabilities are summed in the sorted order), the stream of the samples does for unsorted `data` (the sampler is the inversion of a uniform number over the sorted values); the samples of sorted `data` with `probabilities` are the same as before. Without `probabilities` the Commons object itself keeps the values in the order of a hash map, so its `icdf` and samples were wrong even for sorted `data`; the equal probabilities are now passed explicitly.

  Throws `ex-info` (`:type :invalid-input`, `:property` `:data` or `:probabilities`) for empty `data` (also `nil`), `data` which is not a sequence of numbers (or has `NaN`), `probabilities` of a different count than `data`, with a negative, `NaN`, infinite or non-number entry and with a zero (or infinite) sum.

  Returns an `EnumeratedRealDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  Choosing the finite support family:

  - `enumerated-real`, [[enumerated-int]]: Apache Commons objects, `data` are doubles (`int` values, an `int` range), [[distribution-id]] `:enumerated-real`/`:enumerated-int`; the other three are implemented in this library.
  - [[real-discrete-distribution]]/[[real-discrete]]: the same law as `enumerated-real`, the pdf is evaluated at an exact value of `data` (any other number gives `0.0`), `cdf` is a step function on the reals.
  - [[integer-discrete-distribution]]/[[integer-discrete]]: integers of the long range, the pmf accepts any number (a non-integer gives `0.0`); `icdf` returns longs.
  - [[categorical-distribution]]/[[categorical]]: any values (keywords, strings), `mean` and `variance` are `NaN`; `cdf` and `icdf` follow the order of the first appearance of the distinct values in `data`.

  See also [[distribution]], [[enumerated-int]], [[real-discrete-distribution]], [[empirical]]."
  (^EnumeratedRealDistribution [] (enumerated-real nil))
  (^EnumeratedRealDistribution [{:keys [data probabilities rng]
                                 :or {data [1.0]}}]
   (let [[data probabilities] (enumerated-input 'enumerated-real data probabilities :real)
         ^RandomGenerator r (or rng (JDKRandomGenerator.))
         ;; the constructor without probabilities counts the values in a hash map and the support is not sorted,
         ;; the same counts (the probabilities `count / n`) are passed in the sorted order
         counts (when-not probabilities
                  (let [^java.util.TreeMap tm (java.util.TreeMap.)]
                    (run! (fn [x] (let [k (double x)]
                                    (.put tm k (m/inc (long (or (.get tm k) 0))))))
                          data)
                    (vec tm)))
         pairs (if probabilities
                 (sort-by first (fn [^double a ^double b] (Double/compare a b))
                          (map vector (map double data) probabilities))
                 (map (fn [[x c]] [x (m// (double c) (count data))]) counts))
         ^doubles xs (m/seq->double-array (map first pairs))]
     (EnumeratedRealDistribution. r xs ^doubles (m/seq->double-array (map second pairs))))))

(add-distr-method enumerated-real)

(defn enumerated-int
  "Creates an enumerated integer distribution object over an explicit, finite set of integer values.

  The distribution puts all of its probability mass on the values listed in `data`. Each value's probability is given by the corresponding entry in `probabilities`, or, when `probabilities` is not supplied, all values are treated as equally likely. Repeated values in `data` accumulate probability mass.

  Parameters (single, optional map):

  - `data` (sequence of integers): the finite set of values the distribution can take, integers within the int range `[-2^31, 2^31-1]` (an integral double like `2.0` is accepted); a value outside of the range used to wrap silently. Default: `[1]`.
  - `probabilities` (sequence of doubles): probability associated with each corresponding value in `data`, finite and not negative, with a positive sum; does not need to be normalized. Default: `nil`, meaning uniform probabilities.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Throws `ex-info` (`:type :invalid-input`, `:property` `:data` or `:probabilities`) for empty `data` (also `nil`), `data` which is not a sequence of integers within the `int` range, `probabilities` of a different count than `data`, with a negative, `NaN`, infinite or non-number entry and with a zero (or infinite) sum. Unlike [[enumerated-real]] the order of `data` does not matter for the Commons object, which sorts the support itself.

  Returns an `EnumeratedIntegerDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions. For the differences between the finite support families see [[enumerated-real]].

  See also [[distribution]], [[enumerated-real]], [[integer-discrete-distribution]], [[uniform-int]]."
  (^EnumeratedIntegerDistribution [] (enumerated-int nil))
  (^EnumeratedIntegerDistribution [{:keys [data probabilities rng]
                                    :or {data [1]}}]
   (let [[data probabilities] (enumerated-input 'enumerated-int data probabilities :integer)
         ^RandomGenerator r (or rng (JDKRandomGenerator.))
         ^ints xs (int-array (map long data))]
     (if probabilities
       (EnumeratedIntegerDistribution. r xs (m/seq->double-array probabilities))
       (EnumeratedIntegerDistribution. r xs)))))

(add-distr-method enumerated-int)

(defn bernoulli
  "Creates a Bernoulli distribution object.

  The Bernoulli distribution is a discrete distribution over the two outcomes `0` and `1`, taking the value `1` with probability `p` and `0` with probability `(- 1.0 p)`. It is the special case of the [[binomial]] distribution with a single trial and is commonly used to model a single yes/no or success/failure event.

  Parameters (single, optional map):

  - `p` (double): probability of success (drawing `1`). Default: `0.5`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `BinomialDistribution` object (with a single trial) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions. Because the object is a Commons `BinomialDistribution`, [[distribution-id]] is `:binomial` and [[distribution-parameters]] is `[:trials :p :rng]`, although only `p` and `rng` are accepted.

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

  - `trials` (integer): number of independent trials, an integer within the `int` range (an integral double like `20.0` is accepted, `nil` is the default). Default: `20`.
  - `p` (double): probability of success on each trial. Default: `0.5`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Throws `ex-info` (`:type :invalid-input`, `:property :trials`) when `trials` is not an integer within the `int` range: `20.7`, `3000000000`, `NaN`, infinity, a string (the value used to be truncated or wrapped silently). Negative `trials` or `p` outside of `[0, 1]` are rejected by Apache Commons Math.

  Returns a `BinomialDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[bernoulli]], [[pascal]], [[hypergeometric]]."
  (^BinomialDistribution [] (binomial nil))
  (^BinomialDistribution [{:keys [trials ^double p rng]
                           :or {p 0.5}}]
   (BinomialDistribution. (or rng (JDKRandomGenerator.)) (int-parameter 'binomial :trials trials 20) p)))

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

  - `population-size` (integer): total size of the population being sampled from. Default: `100`.
  - `number-of-successes` (integer): number of success elements present in the population. Default: `50`.
  - `sample-size` (integer): number of elements drawn from the population without replacement. Default: `25`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  The three sizes have to be integers within the `int` range (an integral double like `20.0` is accepted, `nil` is the default), otherwise `ex-info` (`:type :invalid-input`, `:property` the size) is thrown; a non-integer or an out-of-range value used to be truncated or wrapped silently. Inconsistent sizes (a negative size, more successes or draws than the population) are rejected by Apache Commons Math.

  Returns a `HypergeometricDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[binomial]], [[pascal]]."
  (^HypergeometricDistribution [] (hypergeometric nil))
  (^HypergeometricDistribution [{:keys [population-size number-of-successes sample-size rng]}]
   (HypergeometricDistribution. (or rng (JDKRandomGenerator.))
                                (int-parameter 'hypergeometric :population-size population-size 100)
                                (int-parameter 'hypergeometric :number-of-successes number-of-successes 50)
                                (int-parameter 'hypergeometric :sample-size sample-size 25))))

(add-distr-method hypergeometric)

(defn pascal
  "Creates a Pascal distribution object.

  The Pascal distribution is a discrete distribution (a form of the negative binomial distribution) over the non-negative integers, modelling the number of failures observed before accumulating `r` successes in a sequence of independent yes/no trials, each succeeding with probability `p`. It generalizes the [[geometric]] distribution to more than one required success.

  Parameters (single, optional map):

  - `r` (integer): number of successes to accumulate, an integer within the `int` range (an integral double like `20.0` is accepted, `nil` is the default). Default: `20`.
  - `p` (double): probability of success on each trial. Default: `0.5`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Throws `ex-info` (`:type :invalid-input`, `:property :r`) when `r` is not an integer within the `int` range (it used to be truncated or wrapped silently); use [[negative-binomial]] for a non-integer `r`. A non-positive `r` or `p` outside of `(0, 1]` are rejected by Apache Commons Math.

  Returns a `PascalDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[geometric]], [[binomial]], [[hypergeometric]]."
  (^PascalDistribution [] (pascal nil))
  (^PascalDistribution [{:keys [r ^double p rng]
                         :or {p 0.5}}]
   (PascalDistribution. (or rng (JDKRandomGenerator.)) (int-parameter 'pascal :r r 20) p)))

(add-distr-method pascal)

(defn poisson
  "Creates a Poisson distribution object.

  The Poisson distribution is a discrete distribution over the non-negative integers, modelling the number of events occurring in a fixed interval when events happen independently at a constant average rate. It is commonly used for count data, such as arrivals, defects or occurrences per unit of time or space.

  Parameters (single, optional map):

  - `p` (double): the distribution's rate parameter (mean number of events, traditionally denoted lambda). Default: `0.5`.
  - `epsilon` (double): convergence criterion used internally when computing probabilities. Default: the implementation's default epsilon.
  - `max-iterations` (integer): maximum number of iterations used internally when computing probabilities, an integer within the `int` range (`ex-info` otherwise, `nil` is the default). Default: the implementation's default maximum.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a `PoissonDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[binomial]], [[exponential]]."
  (^PoissonDistribution [] (poisson nil))
  (^PoissonDistribution [{:keys [^double p ^double epsilon max-iterations rng]
                          :or {p 0.5 epsilon PoissonDistribution/DEFAULT_EPSILON}}]
   (PoissonDistribution. (or rng (JDKRandomGenerator.)) p epsilon
                         (int-parameter 'poisson :max-iterations max-iterations PoissonDistribution/DEFAULT_MAX_ITERATIONS))))

(add-distr-method poisson)

(defn uniform-int
  "Creates a discrete uniform distribution object over the integers from `lower` to `upper`, inclusive.

  Every integer in the `[lower, upper]` range is equally likely to be drawn. It is the discrete analogue of [[uniform-real]].

  Parameters (single, optional map):

  - `lower` (integer): lowest value in the support, inclusive, an integer within the `int` range. Default: `0`.
  - `upper` (integer): highest value in the support, inclusive, an integer within the `int` range. Default: `Integer/MAX_VALUE`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Throws `ex-info` (`:type :invalid-input`, `:property` `:lower` or `:upper`) when a bound is not an integer within the `int` range (it used to be truncated or wrapped silently; `nil` is the default). `lower > upper` is rejected by Apache Commons Math.

  The default `upper` makes the support `0 .. 2147483647` (2^31 values), which overflows the `int` arithmetic of Apache Commons Math in the pmf: `pdf` is then negative (`-4.66e-10` instead of `4.66e-10`) and `lpdf` is `NaN`. `cdf` (`1.397e-9` at 2), `mean` and sampling are not affected. Pass an `upper` below `Integer/MAX_VALUE` (`upper - lower + 1` has to be below 2^31) when the pmf is needed (`:upper 2147483646` gives `4.66e-10`).

  Returns a `UniformIntegerDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[uniform-real]], [[enumerated-int]]."
  (^UniformIntegerDistribution [] (uniform-int nil))
  (^UniformIntegerDistribution [{:keys [lower upper rng]}]
   (UniformIntegerDistribution. (or rng (JDKRandomGenerator.))
                                (int-parameter 'uniform-int :lower lower 0)
                                (int-parameter 'uniform-int :upper upper Integer/MAX_VALUE))))

(add-distr-method uniform-int)

(defn zipf
  "Creates a Zipf distribution object.

  The Zipf distribution is a discrete distribution over the ranks `1` to `number-of-elements`, where the probability of rank `k` is proportional to `k` raised to the negative `exponent`. It is commonly used to model frequency data in which a small number of elements occur very often and the rest occur rarely, such as word frequencies in natural language or city population sizes.

  Parameters (single, optional map):

  - `number-of-elements` (integer): number of elements `N`, ie. the highest rank in the distribution support, an integer within the `int` range (`ex-info` otherwise, `nil` is the default). Default: `100`.
  - `exponent` (double): exponent characterizing the distribution (also known as `s`); higher values make the distribution decay faster for higher ranks. Default: `3.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Throws `ex-info` (`:type :invalid-input`, `:property :number-of-elements`) when `number-of-elements` is not an integer within the `int` range (it used to be truncated or wrapped silently). A non-positive value is rejected by Apache Commons Math.

  Returns a `ZipfDistribution` object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[pascal]], [[geometric]], [[hypergeometric]]."
  (^ZipfDistribution [] (zipf nil))
  (^ZipfDistribution [{:keys [number-of-elements ^double exponent rng]
                       :or {exponent 3.0}}]
   (ZipfDistribution. (or rng (JDKRandomGenerator.)) (int-parameter 'zipf :number-of-elements number-of-elements 100) exponent)))

(add-distr-method zipf)

(defn multi-normal
  "Creates a multivariate normal (Gaussian) distribution object.

  The multivariate normal distribution is a continuous distribution over `n`-dimensional real vectors, generalizing the [[normal]] distribution to multiple, possibly correlated dimensions. It is fully described by a mean vector and a covariance matrix, and is central to multivariate statistics, sampling correlated random vectors, and many machine learning models.

  Parameters (single, optional map):

  - `means` (sequence of finite doubles): mean vector, one value per dimension. Default: a zero vector of the dimension of `covariances` (`nil` is the same as missing).
  - `covariances` (sequence of sequences of finite doubles): covariance matrix; must be square, symmetric (within `1e-10` of the largest entry) and positive definite, with the dimension of `means`. A singular (positive semi-definite only) matrix is rejected. Default: the identity matrix of the dimension of `means`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates a 2-dimensional standard normal distribution with zero means and identity covariance. When only `covariances` is given the dimension is its size, when both are given they must have the same size.

  Throws `ex-info` (`:type :invalid-input`, `:property` `:means` or `:covariances`) for a covariance matrix which is not square, ragged, empty, non-finite, asymmetric, not positive definite (or numerically singular), for invalid `means` and for a size mismatch of `means` and `covariances`.

  Returns a `MultivariateNormalDistribution` object which can be used with [[pdf]], [[lpdf]], [[sample]], [[means]], [[covariance]], [[dimensions]] and other distribution protocol functions. A point given to [[pdf]] or [[lpdf]] is a sequence with one number per dimension, any other point (a wrong length, a number, `nil`) throws `ex-info`. [[lpdf]] is computed in log space when the density underflows. Unlike univariate distributions it exposes no `cdf`, `icdf`, [[mean]], [[variance]] or bounds (these throw), and [[sample]] returns a vector rather than a scalar.

  See also [[distribution]], [[normal]], [[means]], [[covariance]]."
  ([] (multi-normal nil))
  ([{:keys [means covariances rng]}]
   (distr/multi-normal means covariances rng)))

(add-distr-method multi-normal)

;; SSJ

(defn anderson-darling
  "Creates a distribution object for the Anderson-Darling goodness-of-fit test statistic.

  Given a sample of `n` independent uniform(0,1) random variables, the Anderson-Darling statistic compares the sorted sample values against the ideal uniform distribution, weighting discrepancies in the tails more heavily than the Kolmogorov-Smirnov statistic does. This distribution describes the sampling distribution of that statistic for a sample of size `n`, and is used to obtain critical values or p-values when testing whether a sample follows a fully specified distribution.

  Parameters (single, optional map):

  - `n` (integer): sample size the statistic is computed for, an integer within the `int` range (`ex-info` `:property :n` otherwise - a non-integer, `NaN`, infinity, an out-of-range value used to be truncated or wrapped silently; `nil` is the default). Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  The `cdf`, `icdf` and the sampler are approximations for a small `n` (the absolute error of the `cdf` against a Monte Carlo simulation is about `1.2e-2` for `n` 2; SSJ uses the Marsaglia-Marsaglia algorithm, equal to R's `goftest::pAD` to `1e-5` for `n` of 5 and more). [[mean]] and [[variance]] are not implemented by SSJ and throw `UnsupportedOperationException` (the mean of the statistic is `1` for every `n`).

  Returns a distribution object (backed by the SSJ `AndersonDarlingDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[anderson-darling-quick]]."
  ([] (anderson-darling nil))
  ([{:keys [n rng]}]
   (distr/ssj-continuous :anderson-darling (AndersonDarlingDist. (int-parameter 'anderson-darling :n n 1)) rng [:n :rng])))

(add-distr-method anderson-darling)

(defn anderson-darling-quick
  "Creates a distribution object for the Anderson-Darling goodness-of-fit test statistic, using a faster computation algorithm.

  This is a variant of [[anderson-darling]] describing the same underlying statistic for a sample of size `n`, but relying on an alternative, quicker algorithm to evaluate the distribution's functions. It is preferable when the distribution has to be evaluated many times and computation speed matters more than using the original, reference algorithm.

  Parameters (single, optional map):

  - `n` (integer): sample size the statistic is computed for, an integer within the `int` range (`ex-info` `:property :n` otherwise - a non-integer, `NaN`, infinity, an out-of-range value used to be truncated or wrapped silently; `nil` is the default). Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  The `cdf`, `icdf` and the sampler are approximations for a small `n`, and the lower tail is not more precise than that of [[anderson-darling]] for `n` up to 10, whatever the name suggests (the absolute error of the `cdf` against a Monte Carlo simulation is `2.8e-2` at the probability `0.001` for `n` 2 and `2.1e-3` for `n` 5). [[mean]] and [[variance]] are not implemented by SSJ and throw `UnsupportedOperationException`.

  Returns a distribution object (backed by the SSJ `AndersonDarlingDistQuick` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[anderson-darling]]."
  ([] (anderson-darling-quick nil))
  ([{:keys [n rng]}]
   (distr/ssj-continuous :anderson-darling-quick (AndersonDarlingDistQuick. (int-parameter 'anderson-darling-quick :n n 1)) rng [:n :rng])))

(add-distr-method anderson-darling-quick)

(defn beta-symmetrical
  "Creates a symmetrical beta distribution object.

  The symmetrical beta distribution is a continuous distribution supported on the interval `[0, 1]`, the special case of the [[beta]] distribution where both shape parameters are equal to `alpha`. Being symmetrical, its density is centered and mirrored around `0.5`.

  Parameters (single, optional map):

  - `alpha` (double): shared shape parameter, `alpha = beta` in the general beta distribution; a finite number greater than `0`. Default: `2.0`.
  - `d` (integer): approximate number of decimal digits of precision used internally when computing the distribution, complementary distribution and inverse functions; an integer from `0` to `35` (the results for `d` of 15 and more are equal). Default: `14`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Throws `ex-info` (`:type :invalid-input`, `:property` `:alpha` or `:d`) when `alpha` is not a finite number greater than `0` or `d` is not an integer from `0` to `35` (a non-integer or an out-of-range `d` used to be truncated, wrapped or to fail with an `ArrayIndexOutOfBoundsException` on the first use; `nil` is the default for `d`).

  Returns a distribution object (backed by the SSJ `BetaSymmetricalDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[beta]]."
  ([] (beta-symmetrical nil))
  ([{:keys [^double alpha d rng]
     :or {alpha 2.0}}]
   (check-parameter! 'beta-symmetrical :alpha alpha (and (finite-number? alpha) (m/pos? alpha)) "a finite number > 0")
   (let [d (int-parameter 'beta-symmetrical :d d 14)]
     (check-parameter! 'beta-symmetrical :d d (m/<= 0 d 35) "an integer from 0 to 35")
     (distr/ssj-continuous :beta-symmetrical (BetaSymmetricalDist. alpha d) rng [:alpha :d :rng]))))

(add-distr-method beta-symmetrical)

(defn chi
  "Creates a chi distribution object.

  The chi distribution is a continuous distribution over non-negative reals, describing the distribution of the square root of a sum of squares of `nu` independent standard normal random variables. It generalizes the [[normal]] distribution's absolute value (for `nu` equal to `1`) and the Rayleigh distribution (for `nu` equal to `2`), and is related to the [[chi-squared]] distribution.

  Parameters (single, optional map):

  - `nu` (integer): number of degrees of freedom, ie. the number of underlying standard normal variables; an integer from `1` to `2^31 - 1`. Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Throws `ex-info` (`:type :invalid-input`, `:property :nu`) when `nu` is not an integer within the `int` range (`2.5`, `3000000000`, `NaN`; `(chi {:nu 2.5})` used to be the chi distribution with `nu` 2); `nu` below `1` is rejected by SSJ (`IllegalArgumentException`).

  Returns a distribution object (backed by the SSJ `ChiDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[chi-squared]], [[normal]]."
  ([] (chi nil))
  ([{:keys [nu rng]}]
   (distr/ssj-continuous :chi (ChiDist. (int-parameter 'chi :nu nu 1)) rng [:nu :rng])))

(add-distr-method chi)

(defn chi-squared-noncentral
  "Creates a noncentral chi-squared distribution object.

  The noncentral chi-squared distribution is a continuous distribution over non-negative reals, generalizing the [[chi-squared]] distribution to the case where the underlying normal random variables have a nonzero, shared mean (encoded through the noncentrality parameter `lambda`). It is used, among others, in power calculations for chi-squared tests and in analyses involving sums of squares of non-centered normal variables.

  Parameters (single, optional map):

  - `nu` (double): number of degrees of freedom. Default: `1.0`.
  - `lambda` (double): noncentrality parameter, has to be positive; `lambda = 0.0` is NOT accepted (SSJ throws `IllegalArgumentException` \"lambda = 0\"), use [[chi-squared]] for the central distribution. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Accuracy (measured against R `pchisq`/`qchisq` with `ncp` and scipy `ncx2`): the absolute error of `cdf` is about `4e-13` for `lambda` up to 500, jumps to `1.9e-6` at `lambda = 1000` (`nu` 10, near the median) and decreases slowly afterwards (`6.7e-7` at 2000, `1.7e-7` at 5000, `2.1e-8` at 2e4); `pdf` is accurate to `1e-12` (relative) for `lambda` up to 1000 but `NaN` for `x` below about `1e-16` when `nu = 1`; the `icdf` at small probabilities is inexact (relative error `1.8e-5` at `p = 1e-3`, `1.5e-11` instead of `4.3e-20` at `1e-10`) and `lpdf` is `-Inf` once the log-density is below `-745`. SSJ throws `IllegalArgumentException` (not `ex-info`) for a non-positive `nu` or `lambda`.

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

  - `n` (integer): sample size the statistic is computed for, an integer within the `int` range (`ex-info` `:property :n` otherwise - a non-integer, `NaN`, infinity, an out-of-range value used to be truncated or wrapped silently; `nil` is the default). Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  The `cdf`, `icdf` and the sampler are approximations for a small `n` (the absolute error of the `cdf` against a Monte Carlo simulation is `9.1e-3` for `n` 2 and `2.6e-3` for `n` 5). [[mean]] and [[variance]] are exact: `1/6` and `(4n - 3) / (180 n)`.

  Returns a distribution object (backed by the SSJ `CramerVonMisesDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[anderson-darling]], [[anderson-darling-quick]]."
  ([] (cramer-von-mises nil))
  ([{:keys [n rng]}]
   (distr/ssj-continuous-no-pdf :cramer-von-mises (CramerVonMisesDist. (int-parameter 'cramer-von-mises :n n 1)) rng [:n :rng])))

(add-distr-method cramer-von-mises)

(defn erlang
  "Creates an Erlang distribution object.

  The Erlang distribution is a continuous distribution over non-negative reals, describing the sum of `k` independent [[exponential]] random variables each with rate `lambda`. It is the special case of the [[gamma]] distribution restricted to an integer shape parameter, and commonly models waiting times for multiple sequential events, such as queueing systems.

  Parameters (single, optional map):

  - `k` (integer): number of exponential stages summed together (shape parameter); an integer from `1` to `2^31 - 1` (`ex-info` `:property :k` for a non-integer such as `2.5`, which used to be truncated to 2, or an out-of-range value; `k < 1` is rejected by SSJ). Default: `2`.
  - `lambda` (double): rate parameter of each underlying exponential stage. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object (backed by the SSJ `ErlangDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[exponential]], [[gamma]]."
  ([] (erlang nil))
  ([{:keys [k ^double lambda rng]
     :or {lambda 1.0}}]
   (distr/ssj-continuous :erlang (ErlangDist. (int-parameter 'erlang :k k 2) lambda) rng [:k :lambda :rng])))

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

  - `mu` (double): mean of the underlying normal distribution, before folding; must be non-negative (SSJ throws `IllegalArgumentException` \"mu < 0\" for a negative `mu`, although the folded law is symmetric in `mu`: pass `(abs mu)`). Default: `0.0`.
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

  - `n` (integer): total number of equally spaced rates the distribution is built from; an integer within the `int` range, at least `k`. Default: `1`.
  - `k` (integer): number of phases, out of the `n` available, that are summed together; an integer from `1` to `n`. Default: `1`.
  - `h` (double): common spacing of the rates, a finite number greater than `0`. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  The phases have the rates `(n + 1 - i) * h` for `i` = 1 .. `k`, i.e. the `k` largest of the `n` rates `h, 2h, ..., n*h`: `n*h, (n-1)*h, ..., (n-k+1)*h` (for example `n` 3, `k` 2, `h` 1 is the sum of exponentials with the rates 3 and 2, mean `1/3 + 1/2`). `k = n` gives all the rates `h .. n*h`; the rates are distinct, unlike for [[erlang]].

  Throws `ex-info` (`:type :invalid-input`, `:property` `:n`, `:k` or `:h`) when `n` or `k` is not an integer within the `int` range, when `k < 1` or `n < k`, or when `h` is not a finite number greater than `0` (SSJ reports these only at the first use, `n < k` and `k = 0` with an `IllegalArgumentException` and `h = 0` or `h < 0` with a sample of `Infinity` or a negative number).

  Returns a distribution object (backed by the SSJ `HypoExponentialDistEqual` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[hypoexponential]], [[erlang]], [[exponential]]."
  ([] (hypoexponential-equal nil))
  ([{:keys [n k ^double h rng]
     :or {h 1.0}}]
   (let [n (int-parameter 'hypoexponential-equal :n n 1)
         k (int-parameter 'hypoexponential-equal :k k 1)]
     (check-parameter! 'hypoexponential-equal :k k (m/>= k 1) "an integer >= 1")
     (check-parameter! 'hypoexponential-equal :n n (m/>= n k) "an integer >= k")
     (check-parameter! 'hypoexponential-equal :h h (and (finite-number? h) (m/pos? h)) "a finite number > 0")
     (distr/ssj-continuous :hypoexponential-equal (HypoExponentialDistEqual. n k h) rng [:n :k :h :rng]))))

(add-distr-method hypoexponential-equal)

(defn hypoexponential
  "Creates a hypoexponential distribution object with explicit rates.

  The hypoexponential distribution is a continuous distribution over non-negative reals, describing the sum of independent [[exponential]] random variables, one for each rate listed in `lambdas`. Unlike the [[erlang]] distribution, the rates need not be equal, which lets the distribution capture more general shapes for multi-stage processes such as sequential tasks or phase-type service times.

  Parameters (single, optional map):

  - `lambdas` (sequence of doubles): rate of each independent exponential phase summed together, a non-empty sequence of finite numbers greater than `0`. Repeated rates are accepted (`[1.0 1.0]` is the [[erlang]] distribution with `k` 2). Default: `[1.0]`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Throws `ex-info` (`:type :invalid-input`, `:property :lambdas`) when `lambdas` is not a non-empty sequence of finite numbers greater than `0` (an empty `lambdas` used to fail with an `ArrayIndexOutOfBoundsException` at the first use).

  Returns a distribution object (backed by the SSJ `HypoExponentialDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[hypoexponential-equal]], [[erlang]], [[exponential]]."
  ([] (hypoexponential nil))
  ([{:keys [lambdas rng]
     :or {lambdas [1.0]}}]
   (check-parameter! 'hypoexponential :lambdas lambdas
                     (and (some? lambdas) (not (number? lambdas)) (not (map? lambdas)) (not (string? lambdas)) (seqable? lambdas)
                          (seq lambdas)
                          (every? (fn [l] (and (finite-number? l) (m/pos? (double l)))) lambdas))
                     "a non-empty sequence of finite numbers > 0")
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

  For `lambda / mu <= 100` the distribution is the SSJ `InverseGaussianDist` object. For a larger ratio (the SSJ `cdf` overflows from about 355 and its `icdf` and samples degenerate) `cdf` and `icdf` are computed in log space and the samples are obtained by inverting `cdf` of uniform numbers; `pdf`, `mean`, `variance` and the bounds still come from SSJ. [[ccdf]] of this distribution is `1 - cdf`.

  Returns a distribution object (backed by the SSJ `InverseGaussianDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[gamma]], [[normal]]."
  ([] (inverse-gaussian nil))
  ([{:keys [^double mu ^double lambda rng]
     :or {mu 1.0 lambda 1.0}}]
   (distr/inverse-gaussian mu lambda rng)))

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

  - `n` (integer): sample size the statistic is computed for, an integer within the `int` range (`ex-info` `:property :n` otherwise - a non-integer, `NaN`, infinity, an out-of-range value used to be truncated or wrapped silently; `nil` is the default). Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  The computation is exact but slow for a large `n` (about 1.4 ms per `icdf` call for `n` 100). [[mean]] and [[variance]] are not implemented by SSJ and throw `UnsupportedOperationException`.

  Returns a distribution object (backed by the SSJ `KolmogorovSmirnovDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[kolmogorov-smirnov+]], [[kolmogorov-smirnov-quick]], [[anderson-darling]], [[cramer-von-mises]]."
  ([] (kolmogorov-smirnov nil))
  ([{:keys [n rng]}]
   (distr/ssj-continuous :kolmogorov-smirnov (KolmogorovSmirnovDist. (int-parameter 'kolmogorov-smirnov :n n 1)) rng [:n :rng])))

(add-distr-method kolmogorov-smirnov)

(defn kolmogorov-smirnov+
  "Creates a distribution object for the one-sided Kolmogorov-Smirnov (D+) goodness-of-fit test statistic.

  This is the one-sided variant of [[kolmogorov-smirnov]], describing the sampling distribution of the largest positive (rather than absolute) deviation between the empirical distribution function of a sample of `n` independent uniform(0,1) variables and the ideal uniform distribution. It is sensitive to departures where the empirical distribution lies above the hypothesized one.

  Parameters (single, optional map):

  - `n` (integer): sample size the statistic is computed for, an integer within the `int` range (`ex-info` `:property :n` otherwise - a non-integer, `NaN`, infinity, an out-of-range value used to be truncated or wrapped silently; `nil` is the default). Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  [[mean]] and [[variance]] are not implemented by SSJ and throw `UnsupportedOperationException`. The `icdf` has an absolute tolerance of about `1e-8` (so the very small quantiles are quantised, `icdf(1e-10)` is about `1e-8`) and the samples inherit it.

  Returns a distribution object (backed by the SSJ `KolmogorovSmirnovPlusDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[kolmogorov-smirnov]], [[kolmogorov-smirnov-quick]]."
  ([] (kolmogorov-smirnov+ nil))
  ([{:keys [n rng]}]
   (distr/ssj-continuous :kolmogorov-smirnov+ (KolmogorovSmirnovPlusDist. (int-parameter 'kolmogorov-smirnov+ :n n 1)) rng [:n :rng])))

(add-distr-method kolmogorov-smirnov+)

(defn kolmogorov-smirnov-quick
  "Creates a distribution object for the (two-sided) Kolmogorov-Smirnov goodness-of-fit test statistic, using a faster computation algorithm.

  This is a variant of [[kolmogorov-smirnov]] describing the same underlying statistic for a sample of size `n`, but relying on an alternative, quicker algorithm to evaluate the distribution's functions, similar in spirit to [[anderson-darling-quick]].

  Parameters (single, optional map):

  - `n` (integer): sample size the statistic is computed for, an integer within the `int` range (`ex-info` `:property :n` otherwise - a non-integer, `NaN`, infinity, an out-of-range value used to be truncated or wrapped silently; `nil` is the default). Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  The `pdf` is a finite difference of the `cdf` with the step `1/64` and the `icdf` has a relative error of `4e-5` to `2.5e-4`. [[mean]] and [[variance]] are not implemented by SSJ and throw `UnsupportedOperationException`.

  Returns a distribution object (backed by the SSJ `KolmogorovSmirnovDistQuick` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[kolmogorov-smirnov]], [[kolmogorov-smirnov+]]."
  ([] (kolmogorov-smirnov-quick nil))
  ([{:keys [n rng]}]
   (distr/ssj-continuous :kolmogorov-smirnov-quick (KolmogorovSmirnovDistQuick. (int-parameter 'kolmogorov-smirnov-quick :n n 1)) rng [:n :rng])))

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

  - `alpha` (double): SHAPE parameter; must be positive. Default: `1.0`.
  - `beta` (double): SCALE parameter (the median); must be positive. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  The cdf is `1 / (1 + (x / beta)^(-alpha))`, equal to scipy's `fisk(c=alpha, scale=beta)`; `(cdf (log-logistic {:alpha 3.0 :beta 2.0}) 2.0)` is `0.5`. (The documentation of the previous versions had the roles of `alpha` and `beta` swapped.) The mean exists for `alpha > 1` and the variance for `alpha > 2` (SSJ throws `IllegalArgumentException` otherwise).

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

  - `a` (double): lower bound of the support, a finite number. Default: `0.0`.
  - `b` (double): upper bound of the support, a finite number greater than `a`. Default: `1.0`.
  - `c` (double): shape parameter, controlling how the density concentrates near `b` (for values greater than `1.0`) or near `a` (for values less than `1.0`); a finite number greater than `0`. Default: `2.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Throws `ex-info` (`:type :invalid-input`, `:property` `:a`, `:b` or `:c`) unless `a < b` and `c > 0` with finite parameters (`a >= b` used to give a wrong mean, `c <= 0` an `IllegalArgumentException` at the first use).

  Returns a distribution object (backed by the SSJ `PowerDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[uniform-real]], [[triangular]]."
  ([] (power nil))
  ([{:keys [^double a ^double b ^double c rng]
     :or {a 0.0 b 1.0 c 2.0}}]
   (check-parameter! 'power :a a (finite-number? a) "a finite number")
   (check-parameter! 'power :b b (and (finite-number? b) (m/> b a)) "a finite number > a")
   (check-parameter! 'power :c c (and (finite-number? c) (m/pos? c)) "a finite number > 0")
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

  The `cdf`, `icdf` and the sampler are approximations for a small `n` (the absolute error of the `cdf` against a Monte Carlo simulation is about `9.4e-2` for `n` 2, `3.2e-2` for `n` 5, `7.7e-3` for `n` 20 and `1.2e-3` for `n` 100). [[mean]] and [[variance]] are not implemented by SSJ and throw `UnsupportedOperationException`.

  Returns a distribution object (backed by the SSJ `WatsonGDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[watson-u]], [[kolmogorov-smirnov]]."
  ([] (watson-g nil))
  ([{:keys [n rng]}]
   (distr/ssj-continuous :watson-g (WatsonGDist. (int-parameter 'watson-g :n n 2)) rng [:n :rng])))

(add-distr-method watson-g)

(defn watson-u
  "Creates a distribution object for the Watson U-squared goodness-of-fit test statistic.

  The Watson U-squared statistic is a variant of the [[cramer-von-mises]] statistic adjusted to be invariant to the choice of origin, making it especially suited for testing goodness-of-fit on circular (directional) data, such as angles or times of day. This distribution describes the sampling distribution of that statistic for a sample of size `n`.

  Parameters (single, optional map):

  - `n` (long): sample size the statistic is computed for; must be at least `2` (the backing SSJ class requires it). Default: `2`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  The `cdf`, `icdf` and the sampler are approximations for a small `n` (the absolute error of the `cdf` against a Monte Carlo simulation is `3.6e-3` for `n` 5). [[mean]] and [[variance]] are exact: `1/12` and `(n - 1) / (360 n)`.

  Returns a distribution object (backed by the SSJ `WatsonUDist` class) which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[watson-g]], [[cramer-von-mises]]."
  ([] (watson-u nil))
  ([{:keys [n rng]}]
   (distr/ssj-continuous :watson-u (WatsonUDist. (int-parameter 'watson-u :n n 2)) rng [:n :rng])))

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

  - `n` (integer): number of independent trials distributed among the categories, an integer from `0` to `2^31 - 1`; it is passed to the distribution as it is, so a non-integer (`2.5`), a negative number and an out-of-range value are rejected (an integral double like `5.0` is accepted); `nil` is the default. Default: `10`.
  - `ps` (sequence of doubles): relative probability of each category, finite and not negative with a positive sum; does not need to sum to `1.0`, since values are normalized internally. Default: `[1 1]`, two equally likely categories (paired with [[dirichlet]]'s own default).
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Throws `ex-info` (`:type :invalid-input`, `:property :n` or `:ps`) for an invalid `n` or `ps` (empty, a negative or non-finite entry, a zero sum).

  Returns a distribution object (the pmf is computed by the SSJ `MultinomialDist` class), with as many dimensions as there are entries in `ps`, which can be used with [[pdf]], [[lpdf]], [[cdf]], [[sample]], [[means]], [[covariance]] and other multivariate distribution functions. Each value returned by `sample` (and accepted by `pdf`/`cdf`) is a vector of non-negative integer counts summing to `n`, one per category. `pdf` is `0.0` for a vector whose sum is not `n` or which has a negative or non-integer count. `cdf` is the joint cumulative probability `P(X_1 <= x_1, ..., X_k <= x_k)` (the coordinates are floored, a negative one gives `0.0`), it is not the pmf: for `n` 5 and `ps` `[1 2 1]` `cdf` at `[3 3 3]` is `0.78125` while `pdf` there is `0`. The quantile function and the moments are not supported: [[icdf]] throws `AbstractMethodError` and [[mean]], [[variance]], [[lower-bound]], [[upper-bound]] throw `IllegalArgumentException` (use [[means]] and [[covariance]]).

  See also [[distribution]], [[binomial]], [[multi-normal]]."
  ([] (multinomial nil))
  ([{:keys [n ps rng]
     :or {ps [1 1]}}]
   ;; `n` is not coerced, a non-integer is rejected by the check of `distr/multinomial`
   (distr/multinomial (if (nil? n) 10 n) ps binomial rng)))

(add-distr-method multinomial)

;; custom

(defn dirichlet
  "Creates a Dirichlet distribution object.

  The Dirichlet distribution is a continuous, multivariate distribution over the open probability simplex, ie. over vectors of positive values summing to `1.0`, generalizing the [[beta]] distribution to more than two categories. It is commonly used as a prior over probability vectors in Bayesian statistics, for example over the category probabilities of a [[multinomial]] distribution.

  Parameters (single, optional map):

  - `alpha` (sequence of doubles, or long): concentration parameter for each dimension of the simplex. Larger values concentrate samples closer to the uniform point, smaller values push mass towards the corners of the simplex. When given as a plain integer `n`, it is treated as a symmetric Dirichlet distribution with `n` dimensions, each with concentration `1.0` (ie. uniform over the simplex). Default: `[1 1]`, ie. uniform over the 2-dimensional simplex.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object, with as many dimensions as entries in `alpha`, which can be used with [[pdf]], [[lpdf]], [[probability]], [[sample]], [[means]], [[covariance]] and other multivariate distribution functions. Sampled vectors are vectors of positive doubles (a zero or one coordinate is moved to the neighbouring double) summing to `1.0` up to rounding.

  `pdf`, `lpdf` and `probability` take a point on the simplex: a sequence with one non-negative coordinate per dimension whose sum is within `1e-8` of `1`; any other point (a wrong length, a negative coordinate, a sum off by more than `1e-8`, `nil`) throws `ex-info` (`:type :invalid-input`, `:property :x`), and a `NaN` coordinate gives `NaN`. The boundary of the simplex is evaluated by the limits: a coordinate `0` with `alpha < 1` gives `##Inf` (`(pdf (dirichlet {:alpha [0.5 0.5]}) [0 1])`), with `alpha = 1` the finite density (`1.0` for `[1 1]`), with `alpha > 1` `0.0`.

  Unsupported: `cdf` and `icdf` throw `AbstractMethodError` and [[mean]], [[variance]], [[lower-bound]], [[upper-bound]] throw `IllegalArgumentException` (use [[means]] and [[covariance]]).

  See also [[distribution]], [[beta]], [[multinomial]], [[gamma]]."
  ([] (dirichlet nil))
  ([{:keys [alpha rng]
     :or {alpha [1 1]}}]
   (distr/dirichlet alpha gamma rng)))

(add-distr-method dirichlet)

(defn continuous-distribution
  "Creates a continuous distribution object estimated from a sample of `data` using kernel density estimation (KDE).

  Rather than assuming a parametric family, this builds a smooth, nonparametric estimate of the underlying density from `data` using a chosen kernel and bandwidth (see `fastmath.kernel.density`), then numerically integrates it to obtain matching `cdf` and `icdf` functions. It is useful for approximating an unknown continuous distribution directly from observations, offering a smoother alternative to [[empirical]].

  Parameters (single, optional map):

  - `data` (sequence of doubles): sample used to estimate the density. Default: `[-1 0 1]`.
  - `kde` (keyword or function): kernel used for the density estimation, either a keyword naming a kernel (eg. `:epanechnikov`, `:gaussian`) or a custom kernel function; see [[fastmath.kernel.density/kernel-density]]. Default: `:epanechnikov`.
  - `bandwidth` (double or keyword): bandwidth used by the kernel density estimator; either a fixed positive number or one of `:nrd`, `:nrd0`, `:nrd-adjust`, `:sj`, `:sj-dpi`, `:rlcv`, `:lcv`, `:lscv` to infer it from `data` (see [[fastmath.kernel.density/bandwidth]]; the cross validation selectors are slow for a large `data`). Default: `:nrd`.
  - `steps` (long): number of subintervals used when numerically integrating the density into a `cdf`/`icdf` pair, see `fastmath.random.distributions/integrate-pdf`. Default: `5000`.
  - `interpolator` (keyword or function): interpolation method used between the integrated `cdf`/`icdf` points; one of `:linear`, `:cubic`, `:monotone`, or a custom function. Default: `:monotone`. The `:cubic` interpolator needs enough `steps`: with 100 steps its `icdf` is garbage (`icdf(0.001)` is `5.1e8` and `icdf(0.999)` is `1.3e9` for 1000 standard normal observations), with 1000 steps it is fine for that sample but not for a tiny one (`icdf(0.001)` is `-59.9` for the 5 observations `[1 2 3 4 8]`), the default 5000 is fine for both; `:monotone` and `:linear` stay sensible with 100 steps.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Throws `ex-info` for `data` with fewer than 2 observations (a bandwidth rule of thumb needs at least 2 points) and for `data` with a non-finite value.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions. `pdf` is the kernel density estimate, `cdf` and `icdf` are interpolated from its numerical integral (the estimate integrates to `0.99999`), so they are accurate to about `1e-5` only (for 1000 standard normal observations `cdf` differs from the integral of `pdf` by up to `2.3e-5`). [[mean]] and [[variance]] are the mean and the variance (with the denominator `n - 1`) of `data`, not the moments of the estimated density: the variance of the density is larger by the variance of the kernel multiplied by `bandwidth^2` (`0.9185` for the sample and `0.9296` for the density of 1000 standard normal observations with the Epanechnikov kernel and the default bandwidth, in which the kernel adds `h^2/5`). [[icdf]] at `0` and `1` are the ends of the integration range, which lie slightly beyond [[lower-bound]] and [[upper-bound]] (by `0.006` in the example above).

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

(defn- negative-binomial-point-mass
  "The negative binomial distribution with `p = 1`: the point mass at `0` (the log-pmf of `fastmath.random.distributions/negative-binomial` is `NaN` there, `0 * log(0)`)."
  [rng]
  (distr/->distribution {:lpdf (fn ^double [^long k] (if (m/zero? k) 0.0 ##-Inf))
                         :cdf (fn ^double [^double k] (if (m/neg? k) 0.0 1.0))
                         :icdf (fn [^double _] 0)
                         :sampler :long
                         :rng rng
                         :dimensions 1
                         :continuous? false
                         :name :negative-binomial
                         :parameters [:r :p :rng]
                         :mean 0.0
                         :variance 0.0
                         :lower-bound 0
                         :upper-bound 0}))

(defn negative-binomial
  "Creates a negative binomial distribution object.

  The negative binomial distribution is a discrete distribution over the non-negative integers, modelling the number of failures observed before accumulating `r` successes in a sequence of independent yes/no trials, each succeeding with probability `p`. Unlike [[pascal]], which is restricted to an integer number of successes, `r` here may be any positive real number (a generalized, or Polya, negative binomial distribution).

  Parameters (single, optional map):

  - `r` (double): number of successes to accumulate; need not be an integer (it used to be truncated to an integer silently), a finite number greater than `0`. Default: `20`.
  - `p` (double): probability of success on each trial, a number from `(0, 1]`; `1.0` gives the point mass at `0`. Default: `0.5`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Throws `ex-info` (`:type :invalid-input`, `:property :r` or `:p`) when `r` is not a finite number greater than `0` or `p` is not in `(0, 1]` (including `NaN`).

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[pascal]], [[geometric]], [[binomial]]."
  ([] (negative-binomial nil))
  ([{:keys [^double r ^double p rng]
     :or {r 20.0 p 0.5}}]
   (if (m/== p 1.0)
     (do (check-parameter! 'negative-binomial :r r (and (finite-number? r) (m/pos? r)) "a finite number > 0")
         (negative-binomial-point-mass rng))
     (distr/negative-binomial r p rng))))

(add-distr-method negative-binomial)

;; Poisson is used by `nbi` (`sigma * mu`) and `nbii` (`sigma`) below this relative excess of the variance
(def ^{:private true :const true :tag 'double} negative-binomial-poisson-limit 1.0e-8)

;; the size `r` from which the negative binomial pmf is evaluated by the Poisson log-pmf and a correction
(def ^{:private true :const true :tag 'double} negative-binomial-large-size 1.0e5)

(defn- stirling-correction
  "Error of the Stirling formula of log-gamma, `log(gamma(z)) - ((z - 1/2) log(z) - z + log(2 pi) / 2)`, for `z >= 10`."
  ^double [^double z]
  (let [z2 (m/* z z)]
    (m// (m/- (m// 1.0 12.0)
              (m// (m/- (m// 1.0 360.0)
                        (m// (m/- (m// 1.0 1260.0) (m// (m// 1.0 1680.0) z2)) z2))
                   z2))
         z)))

(defn- negative-binomial-large-size-distribution
  "Negative binomial distribution with the mean `mu` and the size `r >= 1e5` (`p = r / (r + mu)`).

  The pmf is `log P(k) = k log(mu) - mu - log(k!) + d` with the correction to the Poisson pmf
  `d = ((r - 1/2) log1p(k/r) - k) + (mu - r log1p(mu/r)) + k log1p((k - mu) / (r + mu)) + s(r + k) - s(r)`, where `s` is the error of the Stirling formula.
  Each term is small, the log-gamma difference `log(gamma(r + k)) - log(gamma(r))` (of the order of `r log r`, which is
  not accurate enough for a large `r`) is not evaluated. `cdf` and the quantile are those of `fastmath.random.distributions/negative-binomial`."
  [mu r variance rng parameters]
  (let [mu (double mu)
        r (double r)
        variance (double variance)
        p (m// r (m/+ r mu))
        r+mu (m/+ r mu)
        log-mu (m/log mu)
        s-r (stirling-correction r)
        mu-term (m/- mu (m/* r (m/log1p (m// mu r))))
        sd (m/sqrt variance)
        cdf (fn ^double [^double k]
              (if (m/neg? k)
                0.0
                (special/regularized-beta p r (m/inc (m/rint k)))))]
    (distr/->distribution {:lpdf (fn ^double [^long k]
                                   (if (m/neg? k)
                                     ##-Inf
                                     (let [kd (double k)]
                                       (m/+ (m/- (m/* kd log-mu) mu (m/log-factorial k))
                                            mu-term
                                            (m/- (m/* (m/- r 0.5) (m/log1p (m// kd r))) kd)
                                            (m/* kd (m/log1p (m// (m/- kd mu) r+mu)))
                                            (m/- (stirling-correction (m/+ r kd)) s-r)))))
                           :cdf cdf
                           :icdf (fn [^double q]
                                   (cond
                                     (m/not-pos? q) 0
                                     (m/>= q 1.0) ##Inf
                                     :else (distr/discrete-binary-search cdf q [(long mu) (long sd)])))
                           :sampler :long
                           :rng rng
                           :dimensions 1
                           :continuous? false
                           :name :negative-binomial
                           :parameters parameters
                           :mean mu
                           :variance variance
                           :lower-bound 0
                           :upper-bound Integer/MAX_VALUE})))

(defn- negative-binomial-mean
  "Negative binomial distribution with the mean `mu`, the size `r` and the given `variance`; reports `parameters`."
  [mu r variance rng parameters]
  (let [mu (double mu)
        r (double r)]
    (if (m/< r negative-binomial-large-size)
      (distr/negative-binomial r (m// r (m/+ r mu)) rng parameters)
      (negative-binomial-large-size-distribution mu r variance rng parameters))))

(defn nbi
  "Creates a negative binomial (type I) distribution object, using gamlss-style parameter names.

  Alias for [[negative-binomial]], reparameterized in terms of the mean `mu` and a dispersion parameter `sigma`, matching the naming and parametrization used by R's `gamlss.dist` package (its `NBI` family; in the usual NB1/NB2 nomenclature this is the quadratic, NB2, variance form, [[nbii]] is the linear one). The underlying [[negative-binomial]] parameters are recovered as `r = 1 / sigma` and `p = r / (r + mu)`, giving mean `mu` and variance `(+ mu (* sigma mu mu))`, i.e. a variance that is quadratic in `mu` (for example `mu` 100 and `sigma` 0.01 give the variance 200).

  Parameters (single, optional map):

  - `mu` (double): mean of the distribution, a finite number greater than `0`. Default: `1.0`.
  - `sigma` (double): dispersion parameter, a finite number greater than `0`; larger values give more overdispersion relative to a [[poisson]] distribution with the same mean, which is approached when `sigma * mu` tends to `0`. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Throws `ex-info` (`:type :invalid-input`, `:property` `:mu` or `:sigma`) when `mu` or `sigma` is not a finite number greater than `0` (`0`, a negative number, `NaN`, infinity, `nil`; `sigma = 0` used to give a Poisson distribution silently).

  The variance of `nbi` exceeds the variance of the Poisson distribution by the factor `1 + sigma * mu`, so the Poisson limit needs `sigma * mu` small, not `sigma` (`mu` 1e4 and `sigma` 5e-5 have the variance 15000, not 10000). When `sigma * mu` is below `1.0e-8` the object is the [[poisson]] distribution with rate `mu` (its `distribution-id` is `:poisson`), otherwise a negative binomial distribution (`distribution-id` `:negative-binomial`). The threshold is chosen so that the relative error of the pmf against the exact negative binomial pmf (computed with 50 digits) stays below `1e-6` for `sigma` from `1e-9` to `1e-3` and `mu` of 1, 100 and 1e4, at the values within six standard deviations from `mu`: the largest error is `5.3e-8` for the Poisson objects and `2.7e-10` for the negative binomial ones. The log-gamma differences used by the negative binomial pmf lose digits for a large `r = 1 / sigma` (the error is `4.6e-6` at `r = 1e9`), so for `r` of `1e5` or more the pmf is evaluated as the Poisson log-pmf plus a correction in which no large terms are subtracted. The `cdf` is the regularized incomplete beta function of the double `p`, its absolute error is below `1.2e-6` on the same grid.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[negative-binomial]], [[nbii]], [[poisson]]."
  ([] (nbi nil))
  ([{:keys [mu sigma rng]
     :or {mu 1.0 sigma 1.0}}]
   (check-parameter! 'nbi :mu mu (and (finite-number? mu) (m/pos? (double mu))) "a finite number > 0")
   (check-parameter! 'nbi :sigma sigma (and (finite-number? sigma) (m/pos? (double sigma))) "a finite number > 0")
   (let [mu (double mu)
         sigma (double sigma)]
     (if (m/< (m/* sigma mu) negative-binomial-poisson-limit)
       (poisson {:p mu :rng rng})
       (negative-binomial-mean mu (m// 1.0 sigma) (m/+ mu (m/* sigma mu mu)) rng [:mu :sigma :rng])))))

(add-distr-method nbi)

(defn nbii
  "Creates a negative binomial (type II) distribution object, using gamlss-style parameter names.

  Alias for [[negative-binomial]], reparameterized in terms of the mean `mu` and a dispersion parameter `sigma`, matching the naming and parametrization used by R's `gamlss.dist` package (its `NBII` family). The underlying [[negative-binomial]] parameters are recovered as `r = mu / sigma` and `p = r / (r + mu)`, giving mean `mu` and variance `(* mu (+ 1.0 sigma))`, i.e. a variance that is linear in the dispersion parameter, unlike [[nbi]] whose variance is quadratic in `mu`.

  Parameters (single, optional map):

  - `mu` (double): mean of the distribution, a finite number greater than `0`. Default: `1.0`.
  - `sigma` (double): dispersion parameter, a finite number greater than `0`; larger values give more overdispersion relative to a [[poisson]] distribution with the same mean, which is approached as `sigma` tends to `0`. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Throws `ex-info` (`:type :invalid-input`, `:property` `:mu` or `:sigma`) when `mu` or `sigma` is not a finite number greater than `0` (`0`, a negative number, `NaN`, infinity, `nil`).

  The variance is `mu * (1 + sigma)`, so the Poisson limit needs `sigma` small. When `sigma` is below `1.0e-8` the object is the [[poisson]] distribution with rate `mu` (its `distribution-id` is `:poisson`), otherwise a negative binomial distribution (`distribution-id` `:negative-binomial`) with the accuracy described in [[nbi]]: for `r = mu / sigma` of `1e5` or more the pmf is evaluated in a way which does not lose digits to the difference of the log-gamma values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[negative-binomial]], [[nbi]], [[poisson]]."
  ([] (nbii nil))
  ([{:keys [mu sigma rng]
     :or {mu 1.0 sigma 1.0}}]
   (check-parameter! 'nbii :mu mu (and (finite-number? mu) (m/pos? (double mu))) "a finite number > 0")
   (check-parameter! 'nbii :sigma sigma (and (finite-number? sigma) (m/pos? (double sigma))) "a finite number > 0")
   (let [mu (double mu)
         sigma (double sigma)]
     (if (m/< sigma negative-binomial-poisson-limit)
       (poisson {:p mu :rng rng})
       (negative-binomial-mean mu (m// mu sigma) (m/* mu (m/inc sigma)) rng [:mu :sigma :rng])))))

(add-distr-method nbii)

(defn logarithmic
  "Creates a logarithmic (log-series) distribution object.

  The logarithmic distribution is a discrete distribution over the positive integers `1, 2, 3, ...`, whose probability mass decays proportionally to `(/ (m/pow p k) k)`. It was introduced by Fisher to model species-abundance data (the number of species represented by a given number of individuals) and is also used to build the negative binomial distribution as a Poisson mixture.

  Parameters (single, optional map):

  - `p` (double): shape parameter, strictly between `0.0` and `1.0`; larger values give a heavier tail. Default: `0.5`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Throws `ex-info` (`:type :invalid-input`, `:property :p`) when `p` is not strictly between `0.0` and `1.0` (`NaN` included).

  `pmf` is `0.0` at a non-integer or non-positive argument (`pdf` at `2.5` and at `-1` is `0.0`) and `cdf` floors its argument. `cdf` is the exact running sum of the terms `p^j / j` divided by `-log(1 - p)` (no cancellation, it is `1.0` once the remaining tail is below `1e-17`), and `icdf` and the sampler scan the same running sum, which gives the exact smallest `k` with `cdf(k) >= q`: the cost of one `icdf` call or one draw is proportional to the quantile, around the mean `-p / ((1 - p) log(1 - p))` (so draws for `p` close to `1` are slow), the cost of `cdf` is proportional to `min(k, 40 / (1 - p))`.

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

  Limitation: [[lower-bound]] is `0.0` and [[upper-bound]] is `##Inf` whatever `mu` is (`(lower-bound (half-cauchy {:mu 3.0}))` is `0.0` although the support starts at `3.0`; `icdf(0.0)` is `mu` and `cdf` and `pdf` are `0.0` below `mu`), so the bounds are not the support when `mu` is not `0.0`.

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

  Throws `ex-info` (`:type :invalid-input`, `:property :data` or `:probabilities`) for empty `data` (also `nil`), `data` which are not integers within the long range, `probabilities` of a different count than `data`, with a negative, `NaN`, infinite or non-number entry and with a zero sum.

  `pdf` accepts any number (`0.0` at a non-integer), `cdf` is a step function on the reals and `icdf` is the smallest value of `data` with `cdf >= p` (a long). The support is a long, the Apache Commons object [[enumerated-int]] is limited to the `int` range; for the differences of the finite support families see [[enumerated-real]].

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

  Throws `ex-info` (`:type :invalid-input`, `:property :data` or `:probabilities`) for empty `data` (also `nil`), `data` with a `NaN` or a non-number, `probabilities` of a different count than `data`, with a negative, `NaN`, infinite or non-number entry and with a zero sum.

  `pdf` needs an exact value of `data` (any other number gives `0.0`), `cdf` is a step function on the reals and `icdf` is the smallest value of `data` with `cdf >= p`. The data are sorted by the distribution itself, so, unlike [[enumerated-real]], the order of `data` does not matter; for the differences of the finite support families see [[enumerated-real]].

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

  Throws `ex-info` (`:type :invalid-input`, `:property :data` or `:probabilities`) for empty `data` (also `nil`), `probabilities` of a different count than `data`, with a negative, `NaN`, infinite or non-number entry and with a zero sum.

  The values are compared by Clojure equality (`(pdf d 1)` and `(pdf d 1.0)` differ: `1.0` is not found in the data `[1]` and gets `0.0`; a value outside of `data` has `pdf` `0.0`, `lpdf` `-Inf` and `cdf` `0.0`). There is no order of the values: `cdf` is the cumulative probability in the order of the first appearance of the distinct values in `data` (data `[:a :b :c]` with equal weights give `(cdf d :b)` = 2/3), and `icdf` returns the value at which that cumulative probability reaches `p`. There are no bounds ([[lower-bound]] and [[upper-bound]] throw `AbstractMethodError`). For the differences of the finite support families see [[enumerated-real]].

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

  - `ns` (integer): number of success elements in the population, an integer within the long range (`ex-info` `:property :ns` for a non-integer, `NaN` or out-of-range value, `nil` is the default). Default: `5`.
  - `nf` (integer): number of failure elements in the population, an integer within the long range (`ex-info` `:property :nf` otherwise). Default: `5`.
  - `n` (integer): number of elements drawn without replacement, an integer within the long range (`ex-info` `:property :n` otherwise). Default: `5`.
  - `omega` (double): odds ratio, the relative likelihood of a success versus a failure element being among the draws; `1.0` corresponds to the ordinary hypergeometric distribution. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  The support is the integer range from `(max 0 (- n nf))` to `(min ns n)`.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[wallenius-noncentral-hypergeometric]], [[hypergeometric]]."
  ([] (fishers-noncentral-hypergeometric nil))
  ([{:keys [ns nf n ^double omega rng]
     :or {omega 1.0}}]
   (distr/fishers-noncentral-hypergeometric {:ns (long-parameter 'fishers-noncentral-hypergeometric :ns ns 5)
                                             :nf (long-parameter 'fishers-noncentral-hypergeometric :nf nf 5)
                                             :n (long-parameter 'fishers-noncentral-hypergeometric :n n 5)
                                             :omega omega :rng rng})))

(add-distr-method fishers-noncentral-hypergeometric)

(defn wallenius-noncentral-hypergeometric
  "Creates Wallenius' noncentral hypergeometric distribution object.

  Wallenius' noncentral hypergeometric distribution is a discrete distribution over the number of successes obtained by sequentially drawing `n` elements without replacement from a finite population of `ns` success and `nf` failure elements, where at each individual draw a remaining success element is `omega` times as likely to be picked as a remaining failure element. Unlike [[fishers-noncentral-hypergeometric]], which conditions a pair of counts on their sum, this distribution arises from an explicit biased sequential sampling (urn) process, making it the appropriate model for selection sampling where items are removed one at a time under a persistent bias. Setting `omega` to `1.0` recovers the ordinary [[hypergeometric]] distribution.

  Parameters (single, optional map):

  - `ns` (integer): number of success elements in the population, an integer within the long range (`ex-info` `:property :ns` for a non-integer, `NaN` or out-of-range value, `nil` is the default). Default: `5`.
  - `nf` (integer): number of failure elements in the population, an integer within the long range (`ex-info` `:property :nf` otherwise). Default: `5`.
  - `n` (integer): number of elements drawn without replacement, an integer within the long range (`ex-info` `:property :n` otherwise). Default: `5`.
  - `omega` (double): odds ratio, the relative likelihood of a remaining success versus a remaining failure element being picked at each draw; `1.0` corresponds to the ordinary hypergeometric distribution. Default: `1.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  The support is the integer range from `(max 0 (- n nf))` to `(min ns n)`. Unlike [[fishers-noncentral-hypergeometric]], the probability mass function is evaluated via numerical integration, so construction may be relatively slow.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[fishers-noncentral-hypergeometric]], [[hypergeometric]]."
  ([] (wallenius-noncentral-hypergeometric nil))
  ([{:keys [ns nf n ^double omega rng]
     :or {omega 1.0}}]
   (distr/wallenius-noncentral-hypergeometric {:ns (long-parameter 'wallenius-noncentral-hypergeometric :ns ns 5)
                                               :nf (long-parameter 'wallenius-noncentral-hypergeometric :nf nf 5)
                                               :n (long-parameter 'wallenius-noncentral-hypergeometric :n n 5)
                                               :omega omega :rng rng})))

(add-distr-method wallenius-noncentral-hypergeometric)

;;

(defn reciprocal
  "Creates a reciprocal (log-uniform) distribution object.

  The reciprocal distribution is a continuous distribution over the positive interval `[a, b]` whose density is proportional to `(/ 1.0 x)`, equivalently, its logarithm, `(m/log x)`, is uniformly distributed over `[(m/log a), (m/log b)]`. It is commonly used to sample scale-invariant quantities spanning several orders of magnitude, such as hyperparameters in a search space or physical quantities of unknown scale, since it assigns equal probability mass to intervals of equal relative (rather than absolute) width.

  Parameters (single, optional map):

  - `a` (double): lower bound of the support, a finite number greater than `0`. Default: `1.0`.
  - `b` (double): upper bound of the support, a finite number greater than `a`. Default: `10.0`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Throws `ex-info` (`:type :invalid-input`, `:property :a` or `:b`) unless `0 < a < b` with finite `a` and `b` (`a <= 0` used to give a `NaN` cdf, `a > b` a wrong mean, `a = b` an infinite density).

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[uniform-real]], [[log-normal]], [[pareto]]."
  ([] (reciprocal nil))
  ([{:keys [^double a ^double b rng]
     :or {a 1.0 b 10.0}}]
   (check-parameter! 'reciprocal :a a (and (finite-number? a) (m/pos? a)) "a finite number > 0")
   (check-parameter! 'reciprocal :b b (and (finite-number? b) (m/> b a)) "a finite number > a")
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
   (distr/ex-gaussian {:mu mu :sigma sigma :tau nu :rng rng :normal normal :exponential exponential
                       :parameters [:mu :sigma :nu :rng]})))

(add-distr-method exgaus)

(defn beta-binomial
  "Creates a beta-binomial distribution object.

  The beta-binomial distribution is a discrete distribution over the integers `0` to `n`, obtained as a compound of the [[binomial]] distribution whose success probability `p` is itself random, drawn from a [[beta]] distribution with shape parameters `alpha` and `beta`. Compared to the binomial distribution, it is overdispersed, making it a common model for count data (e.g. successes out of `n` trials) exhibiting more variability than the binomial would predict.

  Parameters (single, optional map):

  - `alpha` (double): first shape parameter of the underlying beta distribution. Default: `0.5`.
  - `beta` (double): second shape parameter of the underlying beta distribution. Default: `0.5`.
  - `n` (integer): number of trials, the upper bound of the support; an integer within the long range (`ex-info` `:property :n` for a non-integer, `NaN` or out-of-range value, `nil` is the default). Default: `10`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[bb]], [[binomial]], [[beta]], [[hypergeometric]]."
  ([] (beta-binomial nil))
  ([{:keys [^double alpha ^double beta n rng]
     :or {alpha 0.5 beta 0.5}}]
   (distr/beta-binomial alpha beta (long-parameter 'beta-binomial :n n 10) rng)))

(add-distr-method beta-binomial)

(defn bb
  "Creates a beta-binomial distribution object, using gamlss-style parameter names.

  Alias for [[beta-binomial]], reparameterized in terms of the mean success probability `mu` and a dispersion parameter `sigma`, matching the naming used by R's `gamlss.dist` package (its `BB` family). The underlying [[beta-binomial]] shape parameters are recovered as `alpha = mu / sigma` and `beta = (1 - mu) / sigma`.

  Parameters (single, optional map):

  - `mu` (double): mean success probability, in `(0, 1)`. Default: `0.5`.
  - `sigma` (double): dispersion parameter; larger values give more overdispersion relative to the [[binomial]] distribution. Default: `1.0`.
  - `bd` (integer): binomial denominator, i.e. number of trials, the upper bound of the support; an integer within the long range (`ex-info` `:property :bd` otherwise, `nil` is the default). Default: `10`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[beta-binomial]]."
  ([] (bb nil))
  ([{:keys [^double mu ^double sigma bd rng]
     :or {mu 0.5 sigma 1.0}}]
   (distr/beta-binomial (m// mu sigma) (m// (m/- 1.0 mu) sigma) (long-parameter 'bb :bd bd 10) rng [:mu :sigma :bd :rng])))

(add-distr-method bb)

(defn zero-inflated-binomial
  "Creates a zero-inflated binomial distribution object, using gamlss-style parameter names.

  The zero-inflated binomial distribution augments the [[binomial]] distribution with an extra point mass at `0`, on top of whatever probability the binomial itself already places there. It is a common model for count data (successes out of `bd` trials) exhibiting more zeros than a plain [[binomial]] would predict. With probability `sigma` the outcome is forced to `0`; with probability `(1 - sigma)` it is drawn from a binomial distribution with success probability `mu` and `bd` trials. Matches the naming and parameterization used by R's `gamlss.dist` package (its `ZIBI` family).

  Parameters (single, optional map):

  - `mu` (double): success probability of the underlying binomial component, in `(0, 1)`. Default: `0.5`.
  - `sigma` (double): zero-inflation probability, in `[0, 1)`, i.e. the extra probability of observing `0` beyond the binomial's own mass at `0`. Default: `0.1`.
  - `bd` (integer): binomial denominator, i.e. number of trials, the upper bound of the support; an integer (`ex-info` `:property :bd` for a non-integer, `NaN` or out-of-range value, `nil` is the default; the `int` range for the binomial based families, the long range for the beta-binomial based ones). Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[binomial]], [[zero-inflated-beta-binomial]]."
  ([] (zero-inflated-binomial nil))
  ([{:keys [^double mu ^double sigma bd rng]
     :or {mu 0.5 sigma 0.1}}]
   (distr/zero-inflated-binomial {:mu mu :sigma sigma :bd (int-parameter 'zero-inflated-binomial :bd bd 1) :rng rng :binomial binomial})))

(add-distr-method zero-inflated-binomial)
(add-distr-method zero-inflated-binomial :zibi)

(defn zero-adjusted-binomial
  "Creates a zero-adjusted binomial distribution object, using gamlss-style parameter names.

  The zero-adjusted binomial distribution is a hurdle model over the integers `0` to `bd`: with probability `sigma` the outcome is exactly `0`, and with probability `(1 - sigma)` it is drawn from a zero-truncated binomial distribution with success probability `mu` and `bd` trials, i.e. the binomial's own positive-value probabilities rescaled to sum to `(1 - sigma)`. Unlike [[zero-inflated-binomial]], which adds an extra point mass at `0` on top of the binomial's own mass there, this distribution replaces that mass entirely, so `sigma` is exactly the probability of observing `0`, not merely an addition to it. Matches the naming and parameterization used by R's `gamlss.dist` package (its `ZABI` family).

  Parameters (single, optional map):

  - `mu` (double): success probability of the underlying binomial component, in `(0, 1)`. Default: `0.5`.
  - `sigma` (double): probability of observing `0`, in `[0, 1)`. Default: `0.1`.
  - `bd` (integer): binomial denominator, i.e. number of trials, the upper bound of the support; an integer (`ex-info` `:property :bd` for a non-integer, `NaN` or out-of-range value, `nil` is the default; the `int` range for the binomial based families, the long range for the beta-binomial based ones). Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[binomial]], [[zero-inflated-binomial]]."
  ([] (zero-adjusted-binomial nil))
  ([{:keys [^double mu ^double sigma bd rng]
     :or {mu 0.5 sigma 0.1}}]
   (distr/zero-adjusted-binomial {:mu mu :sigma sigma :bd (int-parameter 'zero-adjusted-binomial :bd bd 1) :rng rng :binomial binomial})))

(add-distr-method zero-adjusted-binomial)
(add-distr-method zero-adjusted-binomial :zabi)


(defn zero-inflated-beta-binomial
  "Creates a zero-inflated beta-binomial distribution object, using gamlss-style parameter names.

  The zero-inflated beta-binomial distribution augments the [[bb]]/[[beta-binomial]] distribution with an extra point mass at `0`, on top of whatever probability the beta-binomial itself already places there. It is a common model for count data (successes out of `bd` trials) exhibiting more zeros than a plain beta-binomial or [[binomial]] would predict. With probability `nu` the outcome is forced to `0`; with probability `(1 - nu)` it is drawn from a beta-binomial distribution parameterized by mean success probability `mu`, dispersion `sigma` and `bd` trials (as in [[bb]]). Matches the naming and parameterization used by R's `gamlss.dist` package (its `ZIBB` family).

  Parameters (single, optional map):

  - `mu` (double): mean success probability of the underlying beta-binomial component, in `(0, 1)`. Default: `0.5`.
  - `sigma` (double): dispersion parameter of the underlying beta-binomial component; larger values give more overdispersion relative to the [[binomial]] distribution. Default: `0.1`.
  - `nu` (double): zero-inflation probability, in `[0, 1)`, i.e. the extra probability of observing `0` beyond the beta-binomial's own mass at `0`. Default: `0.1`.
  - `bd` (integer): binomial denominator, i.e. number of trials, the upper bound of the support; an integer (`ex-info` `:property :bd` for a non-integer, `NaN` or out-of-range value, `nil` is the default; the `int` range for the binomial based families, the long range for the beta-binomial based ones). Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[bb]], [[beta-binomial]]."
  ([] (zero-inflated-beta-binomial nil))
  ([{:keys [^double mu ^double sigma bd ^double nu rng]
     :or {mu 0.5 sigma 0.1 nu 0.1}}]
   (distr/zero-inflated-beta-binomial {:mu mu :sigma sigma :nu nu :bd (long-parameter 'zero-inflated-beta-binomial :bd bd 1) :rng rng})))

(add-distr-method zero-inflated-beta-binomial)
(add-distr-method zero-inflated-beta-binomial :zibb)

(defn zero-adjusted-beta-binomial
  "Creates a zero-adjusted beta-binomial distribution object, using gamlss-style parameter names.

  The zero-adjusted beta-binomial distribution is a hurdle model over the integers `0` to `bd`: with probability `nu` the outcome is exactly `0`, and with probability `(1 - nu)` it is drawn from a zero-truncated beta-binomial distribution parameterized by mean success probability `mu`, dispersion `sigma` and `bd` trials (as in [[bb]]), i.e. the beta-binomial's own positive-value probabilities rescaled to sum to `(1 - nu)`. Unlike [[zero-inflated-beta-binomial]], which adds an extra point mass at `0` on top of the beta-binomial's own mass there, this distribution replaces that mass entirely, so `nu` is exactly the probability of observing `0`, not merely an addition to it. Matches the naming and parameterization used by R's `gamlss.dist` package (its `ZABB` family).

  Parameters (single, optional map):

  - `mu` (double): mean success probability of the underlying beta-binomial component, in `(0, 1)`. Default: `0.5`.
  - `sigma` (double): dispersion parameter of the underlying beta-binomial component; larger values give more overdispersion relative to the [[binomial]] distribution. Default: `0.1`.
  - `nu` (double): probability of observing `0`, in `[0, 1)`. Default: `0.1`.
  - `bd` (integer): binomial denominator, i.e. number of trials, the upper bound of the support; an integer (`ex-info` `:property :bd` for a non-integer, `NaN` or out-of-range value, `nil` is the default; the `int` range for the binomial based families, the long range for the beta-binomial based ones). Default: `1`.
  - `rng`: random number generator used for sampling. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[bb]], [[beta-binomial]], [[zero-inflated-beta-binomial]]."
  ([] (zero-adjusted-beta-binomial nil))
  ([{:keys [^double mu ^double sigma bd ^double nu rng]
     :or {mu 0.5 sigma 0.1 nu 0.1}}]
   (distr/zero-adjusted-beta-binomial {:mu mu :sigma sigma :nu nu :bd (long-parameter 'zero-adjusted-beta-binomial :bd bd 1) :rng rng})))

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
     (distr/zero-inflated-poisson {:mu nmu :sigma sigma :rng rng :poisson poisson
                                   :parameters [:mu :sigma :rng]}))))

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

  Given a base distribution `distr`, this conditions it on lying within `[left, right]`: the density (probability mass) is the one of `distr` divided by the probability that `distr` assigns to the interval, and is zero outside of it. It is useful for restricting the support of an otherwise unbounded or wide distribution, for example to build a bounded prior in Bayesian modelling.

  Parameters (single, optional map):

  - `distr`: the univariate distribution to truncate. Default: a standard [[normal]] distribution. A continuous distribution, or a discrete one with integer support (the Commons discrete keys, [[negative-binomial]], [[integer-discrete]], [[enumerated-int]] and the other native integer keys). [[enumerated-real]], [[real-discrete]] and [[categorical]] are not supported.
  - `left` (number): lower truncation bound. Default: `nil`, meaning the [[lower-bound]] of `distr`, i.e. no truncation on the left.
  - `right` (number): upper truncation bound. Default: `nil`, meaning the [[upper-bound]] of `distr`, i.e. no truncation on the right.
  - `rng`: random number generator, used only to build the default `distr`; it is ignored when `distr` is given (the generator of `distr` is the source of randomness of the truncated distribution).

  Called with no arguments or with `nil`, creates the distribution with default parameter values.

  The normalising probability of a continuous base is computed from the survival function when the interval lies in the upper half of the distribution, so the tails are accurate: the truncation of a standard normal to `[8, 9]` has the correct density, [[cdf]] and mean `8.1212`. For a discrete base the interval is `[ceil(left), floor(right)]` and its probability is summed from the cdf; the bounds of a discrete base have to be finite on the left.

  `lower-bound` and `upper-bound` are the intersection of `[left, right]` with the support of `distr` (for a discrete base the integers inside), [[distribution-id]] is `:truncated` and [[distribution-parameters]] is `[:distr :left :right :rng]`. [[icdf]] is obtained by a root finding on the truncated cdf (continuous) or by the smallest integer with the cdf not below `p` (discrete). [[mean]] and [[variance]] are computed by numerical (Gauss-Kronrod) integration of the density for a continuous base (the normalised truncated density when the probability of the interval is below `0.1`) and by summation of the probability mass for a discrete one, so the first call may be relatively slow (the result is computed once, on the first request).

  Sampling takes values of `distr` and rejects those outside of the interval while the probability of the interval is at least `1e-3`. For a continuous base with a smaller probability the draw of `distr` is mapped by its own cdf and then by the inverse cdf of the truncated distribution (probability integral transform), so sampling stays fast in the far tails. A discrete base always uses rejection, which gives up after `10^7` rejected draws.

  [[set-seed!]] of the truncated distribution reseeds `distr`, which is shared with the original object: sampling from the truncated distribution advances the generator of `distr`.

  Throws `ex-info` (`:type :invalid-input`) when `distr` is not a univariate distribution or is a discrete distribution without integer support, when `left` or `right` is not a number or `nil` or is `NaN`, when `left > right`, when the interval is empty or does not intersect the support (a continuous base needs `left < right`, a discrete one an integer inside), or when its probability is zero (for example the standard normal truncated to `[40, 41]`). The sampler throws `ex-info` with `:type :sampling-failure` when the rejection limit of a discrete base is reached.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[normal]], [[mixture]], [[ccdf]]."
  ([] (truncated nil))
  ([{:keys [distr left right rng]}]
   (let [distr (or distr (normal {:rng rng}))]
     (distr/truncated distr left right))))

(add-distr-method truncated)

(defn mixture
  "Creates a finite mixture distribution object combining several component distributions.

  A mixture distribution samples one of several component distributions, `distrs`, each chosen with its own probability (proportional to the corresponding entry in `weights`), and then samples from the selected component. Its density, cumulative distribution function and mean are therefore the `weights`-weighted average of those of the individual components, while its variance additionally accounts for the variability between component means, via the law of total variance. Components need not share the same type, domain (discrete or continuous) or support, making mixtures a flexible way to model multimodal or heterogeneous data.

  Parameters (single, optional map):

  - `distrs` (sequence of univariate distributions): the component distributions to mix; a non-empty vector, list or lazy sequence. Default: a single standard [[normal]] distribution.
  - `weights` (sequence of numbers): relative weight of each entry in `distrs`, in the same order; finite and non-negative, with a positive sum and the same count as `distrs`. They need not sum to `1.0`, as they are normalized internally. A component with the weight zero is never sampled and is left out of the bounds. Default: equal weight for every component.
  - `rng`: random number generator used only to select the component. Default: a freshly created generator.

  Called with no arguments or with `nil`, creates the distribution with default parameter values, equivalent to a plain standard normal distribution.

  The components keep their own generators and draw the values with them: `rng` does not make the samples reproducible unless the components are seeded too. To reproduce a sample either create the components with seeded `:rng` and pass a seeded `rng`, or call [[set-seed!]] on the mixture. [[set-seed!]] reseeds the selector with `seed` and gives the component `i` its own seed (the `i`-th `long` of a Mersenne Twister generator started with `seed`), so the components never share a stream. It reseeds the component objects in place.

  [[pdf]] is the weighted average of the densities; for a mixture of discrete and continuous components it adds probability masses and densities, which is not a density of any kind. [[icdf]] is the smallest `x` with `cdf(x) >= p` found by bisection with the tolerance `1e-12` (relative or absolute), so on a plateau of the cdf (a gap between components) it returns the left end of it, to the resolution of the double cdf (`(icdf (mixture {:distrs [(normal) (normal {:mu 100.0})] :weights [1 1]}) 0.5)` is `8.29`, where the first component reaches `1.0`), and it may be relatively slow to evaluate. When all the components are integer-valued discrete distributions, the mixture is discrete, [[sample]] and [[icdf]] return integers (`Long`) and [[icdf]] is exact. [[continuous?]] is `true` when at least one component is continuous and `false` otherwise. [[lower-bound]] and [[upper-bound]] are the smallest lower and the largest upper bound of the components with a positive weight. [[distribution-parameters]] is `[:distrs :weights :rng]`.

  Throws `ex-info` (`:type :invalid-input`) when `distrs` is empty or not a sequence (a map is not accepted), contains a value which is not a univariate distribution, when `weights` is not a sequence, has a different count, contains a negative, infinite, `NaN` or non-number value, or sums to zero.

  Returns a distribution object which can be used with [[pdf]], [[cdf]], [[icdf]], [[sample]], [[mean]], [[variance]] and other distribution protocol functions.

  See also [[distribution]], [[categorical]], [[truncated]]."
  ([] (mixture nil))
  ([{:keys [distrs weights rng]}]
   (let [distrs (or distrs [(normal {:rng rng})])
         weights (or weights (repeat (count distrs) 1.0))]
     (distr/mixture distrs weights rng))))

(add-distr-method mixture)

;;

(defonce ^{:doc "Default normal distribution (`mu` = 0.0, `sd` = 1.0).

It is one shared mutable object created with an unseeded generator: sampling from it advances the stream seen by every user of the namespace, and [[set-seed!]] reseeds it for all of them. Create your own [[normal]] when reproducibility matters."} default-normal (normal))

(defn distribution?
  "Checks if `distr` is a distribution object.

  Parameters:

  - `distr`: any value.

  Returns `true` for the objects created by [[distribution]] and by the distribution constructors of this namespace (including [[truncated]] and [[mixture]]), `false` for everything else: `nil`, numbers, keywords, random number generators. Never throws.

  See also [[distribution]], [[distribution-id]], [[distribution-parameters]]."
  [distr]
  (prot/distribution? distr))

;; protocol proxies
(defmacro ^:private with-numeric-argument
  "Expands to the dispatch `call` hardened for a numeric argument `v`: `nan` for `NaN`, `neg-inf` and `pos-inf` for the infinities. A non-numeric `v` goes to `call` untouched."
  [v nan neg-inf pos-inf call]
  `(if (number? ~v)
     (let [x# (double ~v)]
       (cond
         (m/nan? x#) ~nan
         (m/pos-inf? x#) ~pos-inf
         (m/neg-inf? x#) ~neg-inf
         :else ~call))
     ~call))

(defn cdf
  "Cumulative distribution function: the probability that the random variable is less than or equal to `v`, or (3 arguments) lies in the interval `(v1, v2]`.

  Parameters:

  - `d`: a univariate distribution.
  - `v` (number): the point. For an integer-valued discrete distribution the step function is evaluated at the integer part of `v` (`2.7` gives the same as `2`). For [[categorical]] it is one of the values of the data and the result is the cumulative probability up to it, in the order of the data.
  - `v1`, `v2` (numbers): the ends of the interval, `P(v1 < X <= v2)`, that is `cdf(v2) - cdf(v1)`.

  The 2-argument form never hangs or fails on non-finite numbers: `NaN` gives `NaN`, `##Inf` gives `1.0` and `##-Inf` gives `0.0`, for every distribution (some backends loop forever or throw there). A non-numeric `v` is passed to the distribution untouched. The 3-argument form is dispatched as it is, without these checks.

  Returns a double from `[0, 1]`. For `v1 > v2` the keys backed by Apache Commons Math (for example [[normal]], [[gamma]], [[poisson]], [[binomial]], [[constant]], [[enumerated-real]], [[empirical]]) throw `NumberIsTooLargeException`, any other key returns the negative number `cdf(v2) - cdf(v1)`.

  Most multivariate distributions have no cdf: the call throws for [[multi-normal]] (`IllegalArgumentException`) and [[dirichlet]] (`AbstractMethodError`), [[multinomial]] returns the joint cumulative probability `P(X_1 <= x_1, ..., X_k <= x_k)` (not the probability mass: for `n` 5 and `ps` `[1 2 1]` it is `0.78125` at `[3 3 3]`, where the mass is `0`). For the upper tail use [[ccdf]], `1 - cdf` loses precision there.

  See also [[ccdf]], [[icdf]], [[pdf]], [[probability]]."
  (^double [d v]
   (with-numeric-argument v ##NaN 0.0 1.0 (prot/cdf d v)))
  (^double [d v1 v2] (prot/cdf d v1 v2)))

(defn ccdf
  "Complementary cumulative probability (survival function): the probability that the random variable is greater than `v`.

  Parameters:

  - `d`: a univariate distribution.
  - `v` (number): the point, as for [[cdf]]. `NaN` gives `NaN`, `##Inf` gives `0.0` and `##-Inf` gives `1.0` for every distribution.

  Returns a double from `[0, 1]`. The upper tail is computed without cancellation, so the result keeps its relative accuracy for tiny probabilities (`(ccdf (normal) 8.3)` is `5.2056e-17`), for these keys:

  - the keys backed by SSJ, which use the survival function of the SSJ distribution (except [[folded-normal]] and [[inverse-gaussian]]),
  - [[normal]], [[t]], [[logistic]], [[laplace]], [[cauchy]], [[exponential]], [[weibull]], [[pareto]], [[log-normal]], [[gamma]] and [[chi-squared]], which use closed forms.

  Every other key (all the discrete ones, [[truncated]], [[mixture]], [[folded-normal]], [[inverse-gaussian]] and the remaining native ones) returns `1 - cdf(v)`, which loses all digits for probabilities below about `1e-16`: the result is `0.0` as soon as `cdf` rounds to `1.0`.

  For a non-numeric `v` (a value of a [[categorical]] distribution) the result is `1 - cdf(v)`. A multivariate distribution is not supported, as for [[cdf]].

  See also [[cdf]], [[icdf]]."
  ^double [d v]
  (if (number? v)
    (let [x (double v)]
      (cond
        (m/nan? x) ##NaN
        (m/pos-inf? x) 0.0
        (m/neg-inf? x) 1.0
        :else (distr/ccdf d x)))
    (- 1.0 (double (prot/cdf d v)))))

(defn pdf
  "Probability density function; for a discrete distribution the probability mass function, the probability of the value `v`.

  Parameters:

  - `d`: a distribution.
  - `v`: the point. A number for a univariate distribution. For a multivariate one a sequence with one number per dimension: the point for [[multi-normal]], a point of the simplex (non-negative coordinates summing to 1) for [[dirichlet]], the counts per category for [[multinomial]]. For [[categorical]] one of its values.

  `NaN` gives `NaN`, `##Inf` and `##-Inf` give `0.0` for a numeric `v`. A non-numeric `v` is passed to the distribution untouched.

  Returns a double. It is `0.0` outside of the support for the keys which check it, but some backends return `NaN` or a negative number there. A non-integer `v` given to an integer-valued discrete distribution depends on the backend: the Apache Commons Math keys ([[poisson]], [[binomial]], [[geometric]] ...) use its integer part (`pdf` at `2.5` is `pdf` at `2`), the keys implemented in this library ([[negative-binomial]], [[logarithmic]], [[nbi]], [[integer-discrete]] ...) return `0.0`; [[real-discrete]], [[enumerated-real]] and [[categorical]] need an exact value of the data. Pass integers.

  Throws `ex-info` (`:type :invalid-input`) for a point of a multivariate distribution which has a wrong length or lies outside of the allowed set, a `ClassCastException` or `NullPointerException` for a vector or `nil` given to a univariate distribution.

  See also [[lpdf]], [[probability]], [[cdf]]."
  ^double [d v]
  (with-numeric-argument v ##NaN 0.0 0.0 (prot/pdf d v)))

(defn lpdf
  "Natural logarithm of the density (the probability mass for a discrete distribution) at `v`.

  Parameters:

  - `d`: a distribution.
  - `v`: the point, as for [[pdf]].

  `NaN` gives `NaN`, `##Inf` and `##-Inf` give `##-Inf` for a numeric `v`.

  Returns a double, `##-Inf` where the density is zero. The keys which compute `log` of the density (for example [[laplace]]) underflow to `##-Inf` in the far tails where the log-density is still representable, [[multi-normal]] switches to a log-space formula when its density underflows. The errors are those of [[pdf]].

  See also [[pdf]], [[observe1]], [[log-likelihood]]."
  ^double [d v]
  (with-numeric-argument v ##NaN ##-Inf ##-Inf (prot/lpdf d v)))

(defn icdf
  "Inverse of the cumulative distribution function (quantile function).

  Parameters:

  - `d`: a univariate distribution.
  - `p` (number): the probability, from `[0, 1]`.

  Returns the quantile. For a continuous distribution the `x` such that `cdf(x) = p`, with the accuracy of the backend (the Commons keys which invert the cdf numerically use the `inverse-abs-accuracy` of the constructor, `1e-9` by default). For a discrete distribution the smallest `x` such that `cdf(x) >= p`, an integer for integer-valued ones and the value of the data for [[categorical]]. `p = 0` and `p = 1` give the ends of the support of the backend: `##-Inf` and `##Inf` for an unbounded continuous distribution, for an unbounded discrete one `2147483647` (Commons) or `##Inf` ([[negative-binomial]], [[logarithmic]]).

  Throws `ex-info` (`:type :invalid-input`, `:function icdf`, `:property :p`) when `p` is `NaN` or outside of `[0, 1]`, the same for every distribution (the backends throw different exceptions or return garbage there). `nil` gives `NullPointerException`. Multivariate distributions have no quantile function, the call throws.

  See also [[cdf]], [[ccdf]], [[sample]]."
  [d ^double p]
  (if (and (m/<= 0.0 p) (m/<= p 1.0))
    (prot/icdf d p)
    (throw (invalid-input 'icdf :p p "a number from [0, 1]"
                          (str "icdf, p: has to be a number from [0, 1], got " p)))))

(defn probability
  "Probability of the value `v`: for a discrete distribution the probability mass function, the same as [[pdf]].

  For a continuous distribution this is the density, not a probability (the probability of a single point is zero); use [[cdf]] with 3 arguments for the probability of an interval.

  Parameters:

  - `d`: a distribution.
  - `v`: the point, as for [[pdf]]. `NaN` gives `NaN`, `##Inf` and `##-Inf` give `0.0` for a numeric `v`.

  Returns a double. [[multi-normal]] does not support it (throws).

  See also [[pdf]], [[cdf]]."
  ^double [d v]
  (with-numeric-argument v ##NaN 0.0 0.0 (prot/probability d v)))

(defn sample
  "Draws one random value from the distribution, with the random number generator held by the distribution.

  Parameters:

  - `d`: a distribution.

  Returns a `Double` for a continuous univariate distribution and for the real-valued discrete ones, a `Long` for an integer-valued discrete one, a vector for a multivariate one ([[multi-normal]], [[dirichlet]], [[multinomial]]) and one of the values of the data for [[categorical]]. The call advances the generator of the distribution; the stream is reproducible when the distribution was created with a seeded `:rng` or is reseeded with [[set-seed!]]. A distribution is not thread-safe unless its generator is (see [[synced-rng]]).

  See also [[set-seed!]], [[irandom]], [[drandom]]."
  [d] (prot/sample d))

(defn dimensions
  "Dimensionality of the distribution.

  Parameters:

  - `d`: a distribution.

  Returns a long: `1` for a univariate distribution, the length of the point for [[multi-normal]], [[dirichlet]] and [[multinomial]].

  See also [[sample]], [[means]]."
  ^long [d] (prot/dimensions d))

(defn source-object
  "Returns the object which implements the distribution.

  Parameters:

  - `d`: a distribution.

  Returns the Apache Commons Math distribution (a subclass from `fastmath.java.dist` for [[logistic]], [[gumbel]], [[nakagami]], [[f]] and [[levy]]), the SSJ `ContinuousDistribution` for the keys backed by SSJ, or the distribution itself for the keys implemented in this library (for [[categorical]] the inner [[integer-discrete]] distribution). Use it to reach the methods of the backend.

  See also [[distribution]], [[distribution-id]]."
  [d] (prot/source-object d))

(defn continuous?
  "Does the distribution have a continuous domain?

  Parameters:

  - `d`: a distribution.

  Returns `true` for the distributions on the real line with a density, including [[empirical]], [[kde]] and the multivariate [[multi-normal]] and [[dirichlet]], `false` for the discrete ones ([[enumerated-real]], [[real-discrete]], [[categorical]] and [[multinomial]] included). [[truncated]] follows its base distribution, [[mixture]] is `true` as soon as one component is continuous.

  See also [[discrete?]], [[pdf]]."
  [d]
  (prot/continuous? d))

(defn discrete?
  "Does the distribution have a discrete domain? The negation of [[continuous?]].

  Parameters:

  - `d`: a distribution.

  Returns `true` or `false`.

  See also [[continuous?]]."
  [d]
  (not (prot/continuous? d)))

(defn observe1
  "Log of the probability (density) of one observation. Alias for [[lpdf]].

  Parameters:

  - `d`: a distribution.
  - `v`: the observation.

  Returns a double, see [[lpdf]] for the special values and errors.

  See also [[observe]], [[log-likelihood]]."
  ^double [d v]
  (lpdf d v))

(defn log-likelihood
  "Log likelihood of a sample: the sum of the log densities (log masses) of the observations.

  Parameters:

  - `d`: a distribution.
  - `vs`: a sequence of observations; for a multivariate distribution each observation is a vector (a vector of vectors). `nil` and an empty sequence give `0.0`.

  Returns a double. An observation outside of the support gives `##-Inf`, a `NaN` observation gives `NaN`. The summation stops as soon as the running sum is not finite, but the observations are evaluated in chunks of up to 32, so an invalid element which follows closely after such an observation can still throw.

  Throws `NullPointerException` for a `nil` observation, `ClassCastException` for a vector observation of a univariate distribution and `IllegalArgumentException` when `vs` is not a sequence.

  See also [[likelihood]], [[observe]], [[lpdf]]."
  ^double [d vs]
  (reduce (fn [^double s ^double v] (if (m/invalid-double? s)
                                     (reduced s)
                                     (+ s v))) 0.0 (map #(prot/lpdf d %) vs)))

(defn observe
  "Log likelihood of a sample. Alias for [[log-likelihood]].

  Parameters:

  - `d`: a distribution.
  - `vs`: a sequence of observations, see [[log-likelihood]].

  Returns a double.

  See also [[observe1]], [[likelihood]]."
  ^double [d vs]
  (log-likelihood d vs))

(defn likelihood
  "Likelihood of a sample: the product of the densities (masses) of the observations, `exp` of [[log-likelihood]].

  Parameters:

  - `d`: a distribution.
  - `vs`: a sequence of observations, see [[log-likelihood]].

  Returns a double. It underflows to `0.0` as soon as the log likelihood is below about `-745`: for the standard [[normal]] this happens for roughly 800 observations at the mode and for about 500 typical ones. Compare samples with [[log-likelihood]].

  See also [[log-likelihood]]."
  ^double [d vs]
  (m/exp (log-likelihood d vs)))

(defn mean
  "Mean of the distribution.

  Parameters:

  - `d`: a univariate distribution (use [[means]] for a multivariate one).

  Returns a double. The keys implemented in this library return `##NaN` when the mean does not exist and `##Inf` or `##-Inf` when it diverges (for example [[half-cauchy]] gives `NaN`, [[generalized-pareto]] with `xi >= 1` gives `##Inf`); they never throw. The keys backed by Apache Commons Math or SSJ give what the backend gives: `##NaN` ([[cauchy]], [[t]] with one degree of freedom, [[f]] with small denominator degrees of freedom), `##Inf` ([[levy]], [[pareto]] with `shape <= 1`), an `IllegalArgumentException` ([[frechet]], [[log-logistic]] and [[pearson-6]] where the mean does not exist) or an `UnsupportedOperationException` ([[anderson-darling]], [[anderson-darling-quick]], [[kolmogorov-smirnov]], [[kolmogorov-smirnov+]], [[kolmogorov-smirnov-quick]], [[watson-g]]). [[categorical]] gives `NaN`. Multivariate distributions throw `IllegalArgumentException`.

  See also [[variance]], [[means]]."
  ^double [d] (prot/mean d))

(defn means
  "Mean vector of a multivariate distribution.

  Parameters:

  - `d`: a multivariate distribution ([[multi-normal]], [[dirichlet]], [[multinomial]]).

  Returns a vector with one mean per dimension. A univariate distribution throws `IllegalArgumentException`.

  See also [[covariance]], [[mean]], [[dimensions]]."
  [d] (prot/means d))

(defn variance
  "Variance of the distribution.

  Parameters:

  - `d`: a univariate distribution (use [[covariance]] for a multivariate one).

  Returns a double. The keys implemented in this library return `##NaN` when the variance does not exist and `##Inf` when it diverges, they never throw and never return a negative infinity. The keys backed by Apache Commons Math or SSJ give what the backend gives, which can be `##NaN`, `##Inf`, a negative infinity ([[inverse-gamma]] with `alpha = 1`, [[pearson-6]] with `alpha2 = 1`) or an exception, as listed for [[mean]]. [[categorical]] gives `NaN`. Multivariate distributions throw `IllegalArgumentException`.

  See also [[mean]], [[covariance]]."
  ^double [d] (prot/variance d))

(defn covariance
  "Covariance matrix of a multivariate distribution.

  Parameters:

  - `d`: a multivariate distribution ([[multi-normal]], [[dirichlet]], [[multinomial]]).

  Returns the matrix as a sequence of rows (one row per dimension). A univariate distribution throws `IllegalArgumentException`.

  See also [[means]], [[variance]]."
  [d] (prot/covariance d))

(defn lower-bound
  "Lowest value of the support.

  Parameters:

  - `d`: a univariate distribution.

  Returns a double, `##-Inf` when the support is unbounded below. [[truncated]] reports the intersection of its interval with the support of the base distribution. The bounds of some SSJ keys are not the support ([[power]], [[half-normal]], [[folded-normal]] and [[frechet]] report `##-Inf`). [[categorical]] and the multivariate distributions throw.

  See also [[upper-bound]]."
  ^double [d] (prot/lower-bound d))

(defn upper-bound
  "Highest value of the support.

  Parameters:

  - `d`: a univariate distribution.

  Returns a double, `##Inf` when the support is unbounded above (for the unbounded integer-valued Commons keys such as [[poisson]] it is `2.147483647E9`, the largest `int`). The remarks of [[lower-bound]] apply.

  See also [[lower-bound]]."
  ^double [d] (prot/upper-bound d))

(defn distribution-id
  "Identifier of the distribution as a keyword.

  Parameters:

  - `d`: a distribution.

  Returns the key of [[distribution]] for most distributions. An alias gives its canonical key (`:generalized-extreme-value` for `:gev`, `:integer-discrete-distribution` for `:integer-discrete`, `:zero-inflated-poisson` for `:zip`). The reparametrised and special keys give the id of the object they are built from: `:bernoulli` is `:binomial` (a [[binomial]] with one trial), `:bb` is `:beta-binomial`, `:nbi` and `:nbii` are `:negative-binomial` (`:poisson` when `sigma * mu` (`:nbi`) or `sigma` (`:nbii`) is below `1.0e-8`), `:exgaus` is `:ex-gaussian`, `:zip2` and `:zero-inflated-poisson2` are `:zero-inflated-poisson`. The other aliases of the list of [[distributions-list]] with a different id: `:categorical` (`:categorical-distribution`), `:real-discrete` (`:real-discrete-distribution`), `:kde` (`:continuous-distribution`), `:ge`, `:gg`, `:gh`, `:ghl`, `:gig`, `:gnd`, `:gpd` (the long names of the generalized family), `:zabb`, `:zabi`, `:zaga`, `:zaig`, `:zanbi`, `:zap`, `:zibb`, `:zibi`, `:zinbi` (the long names of the zero-adjusted and zero-inflated families). [[distribution-parameters]] of `:bernoulli` is `[:trials :p :rng]`, of `:nbi` and `:nbii` `[:mu :sigma :rng]` (of the `:poisson` fallback `[:p :epsilon :max-iterations :rng]`) and of `:bb` `[:mu :sigma :bd :rng]`.

  See also [[distribution-parameters]], [[source-object]], [[distribution]]."
  [d] (prot/distribution-id d))

(defn distribution-parameters
  "Keys of the configuration map of the distribution, as a vector of keywords.

  Parameters:

  - `d`: a distribution or a key of [[distribution]] (a keyword; then the distribution with the default parameters is created).

  Returns the keys the constructor accepts, usually with `:rng` last. For the reparametrised keys these are the keys of the key itself (`[:mu :sigma :bd :rng]` for `:bb`, `[:mu :sigma :rng]` for `:nbi` and `:nbii`, `[:mu :sigma :nu :rng]` for `:exgaus`). The exception is `:bernoulli`, which reports the keys of the [[binomial]] object it is built from (`[:trials :p :rng]`) although it accepts `:p` and `:rng` only. For [[truncated]] it is `[:distr :left :right :rng]` and for [[mixture]] `[:distrs :weights :rng]`.

  See also [[distribution-id]], [[distribution]]."
  [d]
  (let [d' (if (keyword? d) (distribution d) d)]
    (prot/distribution-parameters d')))

;;

(defonce ^{:doc "Sorted set of all the keys accepted by [[distribution]] (aliases included).

It is a snapshot taken when the namespace is loaded: a method added to [[distribution]] later does not appear here. Use `(keys (methods distribution))` for the live list."}
  distributions-list
  (into (sorted-set) (keys (methods distribution))))
;;

(defn set-seed
  "Creates a new random number generator of the same algorithm as a given one, seeded with a given seed.

  The given generator is not changed. Use [[set-seed!]] to reseed a generator in place.

  Parameters:

  - `rng` (optional): a generator created by [[rng]] or [[synced-rng]]. Default: [[default-rng]], which is left untouched. A distribution is not supported (it has [[set-seed!]] only).
  - `seed` (optional, integer): any integer within the long range, as for [[rng]]. Without a seed one is drawn from the default generator with [[lrand]], which advances [[default-rng]] by one `long`, and the new generator is effectively unseeded.

  Returns a new generator, reproducible from the seed: `(set-seed (rng :mersenne 1) 42)` equals `(rng :mersenne 42)`. The generator of the default is a new `:jdk` generator. The class of a generator from [[rng]] is kept, a generator from [[synced-rng]] gives a new synchronized generator of the same algorithm, and any other kind of `RandomGenerator` gives a `:jdk` one.

  Throws `ex-info` (`:type :invalid-input`, `:property :seed`) when `seed` is not an integer within the long range: `nil`, `false`, a double, a ratio, a string, `NaN` and a big integer beyond the long range.

  See also [[set-seed!]], [[rng]], [[child-rngs]]."
  ([]
   (prot/set-seed default-rng (lrand)))
  ([seed]
   (prot/set-seed default-rng (check-seed 'set-seed seed)))
  ([rng seed]
   (prot/set-seed rng (check-seed 'set-seed seed))))

(defn set-seed!
  "Reseeds a random number generator in place.

  Parameters:

  - `rng` (optional): a generator created by [[rng]] or [[synced-rng]], or a distribution. Default: [[default-rng]], the shared generator of the whole namespace.
  - `seed` (optional, integer): any integer within the long range, as for [[rng]]. Without a seed one is drawn from the default generator with [[lrand]].

  Returns `rng` itself (a distribution included), which from now on repeats the stream of a fresh generator of the same algorithm created with the same seed. The cached second value of a gaussian draw is discarded. A generator from [[synced-rng]] keeps its class.

  Throws `ex-info` (`:type :invalid-input`, `:property :seed`) when `seed` is not an integer within the long range: `nil`, `false`, a double, a ratio, a string, `NaN` and a big integer beyond the long range.

  See also [[set-seed]], [[rng]], [[default-rng]]."
  ([] (prot/set-seed! default-rng (lrand)))
  ([seed] (prot/set-seed! default-rng (check-seed 'set-seed! seed)))
  ([rng seed] (prot/set-seed! rng (check-seed 'set-seed! seed))))

;;

;; The sampling helpers take an explicit generator, they never use `default-rng`.

(defn- uniform-spacings
  [rng ^long n]
  (let [xs (reductions m/+ (repeatedly (inc n) #(- (m/log (drandom rng)))))
        l (/ ^double (last xs))]
    (map (fn [^double x] (* x l)) (butlast xs))))

(defn- systematic-spacings
  [rng ^long n]
  (let [l (/ 1.0 n)
        d (drandom rng)]
    (map (fn [^long x] (* (+ x d) l)) (range n))))

(defn- stratified-spacings
  [rng ^long n]
  (let [l (/ 1.0 n)]
    (map (fn [^long x] (* (+ x (drandom rng)) l)) (range n))))

(defn- antithetic-sampling
  [rng ^long n]
  (->> (repeatedly (fn [] (let [r1 (drandom rng)]
                           [r1 (- 1.0 r1)])))
       (mapcat identity)
       (take n)))

(defn- jittered-sequence-sampling
  [kind rng ^long n]
  (take n (jittered-sequence-generator kind 1 0.25 rng)))

(def ^:private spacings
  {:uniform uniform-spacings
   :systematic systematic-spacings
   :stratified stratified-spacings
   :antithetic antithetic-sampling
   :r2 (partial jittered-sequence-sampling :r2)
   :sobol (partial jittered-sequence-sampling :sobol)
   :halton (partial jittered-sequence-sampling :halton)})

(def ^:private sampling-methods [:uniform :systematic :stratified :antithetic :r2 :sobol :halton])

(defn- check-n
  "Returns the number of samples as a long or throws `ex-info` when it is not an integer within the long range."
  ^long [n]
  (if (long-integer? n)
    (.longValue ^Number n)
    (throw (invalid-input '->seq :n n "an integer >= 0"
                          (str "The number of samples must be an integer, got " (pr-str n))))))

(defn ->seq
  "Returns a lazy sequence of random samples.

  Works with either a random number generator or a distribution as `source`. A random number generator produces uniform draws, a distribution produces values following it. An optional sampling method replaces the independent draws by a more regular scheme.

  Parameters:

  - `source` (optional): a random number generator (see [[rng]]) or a distribution to sample from. Defaults to [[default-rng]] in the zero-argument form.
  - `n` (optional, integer): the number of samples. When omitted or `nil`, an infinite lazy sequence is returned. A negative `n` gives an empty sequence.
  - `sampling-method` (optional, keyword): selects a sampling scheme applied to the underlying uniform draws, instead of plain independent sampling. Requires `n`; `nil` or `false` falls back to plain sampling. One of:
    - `:uniform` - order statistics of `n` uniform draws, simulating a sorted uniform sample.
    - `:systematic` - low-variance systematic sampling using a single shared random offset.
    - `:stratified` - stratified sampling, one random draw per equal-width stratum.
    - `:antithetic` - antithetic sampling, pairing each draw `r` with its complement `1 - r` (for an odd `n` the last draw is unpaired).
    - `:r2`, `:sobol`, `:halton` - jittered one-dimensional low-discrepancy sequences (see [[jittered-sequence-generator]]), the first point of `:sobol` and `:halton` is the origin before the jitter is added.
  - `rng` (optional, 4-argument form only): the generator which draws the uniforms when `source` is a distribution, checked by [[ensure-rng]]. When it is missing or `nil` a new unseeded generator is used - never [[default-rng]] and never the generator of the distribution - so pass a seeded generator to reproduce the result. It must be `nil` (or missing) when `source` is a random number generator, because the source is already the generator.

  Without a sampling method a random number generator gives doubles from `[0,1)` (as repeated [[drandom]]) and a distribution gives its own variates, using the generator of the distribution. With a sampling method the uniforms are `u` from `[0,1)` for a random number generator, and `(icdf source u)` for a distribution. The uniform values can be exactly `0.0` (for the origin of `:sobol` and `:halton` or a zero draw), where an unbounded distribution has an infinite quantile.

  Returns a lazy sequence of samples.

  Throws `ex-info` (`:type :invalid-input`) when `n` is not an integer within the long range or `nil` (`:property :n`), when `sampling-method` is given without `n` (`:property :n`), when `sampling-method` is not one of the methods above (`:property :sampling-method`, `:expected` the list of the 7 methods), when a non-`nil` `rng` is given for a random number generator source (`:property :rng`), and `ex-info` with `{:rng value}` when `rng` is neither a generator nor `nil`.

  See also [[white-noise]], [[distribution?]], [[icdf]], [[ensure-rng]]."
  ([] (prot/->seq default-rng))
  ([source] (prot/->seq source))
  ([source n]
   (if (nil? n)
     (prot/->seq source)
     (let [n (check-n n)]
       (if (neg? n)
         ()
         (prot/->seq source n)))))
  ([source n sampling-method] (->seq source n sampling-method nil))
  ([source n sampling-method rng]
   (let [distribution-source? (distribution? source)]
     (when (and (not distribution-source?) (some? rng))
       (throw (invalid-input '->seq :rng rng "nil for a random number generator source"
                             "A random number generator source is already the generator, the rng argument must be nil")))
     (cond
       (not sampling-method) (->seq source n)
       (nil? n) (throw (invalid-input '->seq :n n "an integer >= 0 (a sampling method requires n)"
                                      (str "The sampling method " (pr-str sampling-method) " requires the number of samples n")))
       :else (if-let [sampling-fn (get spacings sampling-method)]
               (let [n (check-n n)]
                 (cond
                   (neg? n) ()
                   distribution-source? (map (partial icdf source) (sampling-fn (ensure-rng rng) n))
                   :else (sampling-fn source n)))
               (throw (ex-info (str "Wrong sampling method " (pr-str sampling-method) ", expected one of " (pr-str sampling-methods))
                               {:type :invalid-input
                                :function '->seq
                                :property :sampling-method
                                :value sampling-method
                                :expected sampling-methods
                                :sampling-method sampling-method})))))))

(defn white-noise
  "Generates Gaussian white noise.

  The generated sequence consists of independent random numbers drawn from a normal distribution with mean 0.0 and a specified standard deviation `sigma`: the element is `sigma` times a standard normal draw of [[grandom]], so the stream for `sigma` is the stream for `1.0` multiplied by `sigma`.

  Parameters (the arity decides the meaning of the single argument):

  - `()` - standard normal noise, drawn from a new unseeded `:jdk` generator.
  - `(white-noise x)` - when `x` is a number it is `sigma` and a new unseeded generator is used; otherwise `x` is the generator (`sigma` is `1.0`).
  - `(white-noise rng sigma)` - both given.
  - `rng` (random number generator or `nil`): the generator which draws the numbers, see [[rng]] and [[ensure-rng]]; `nil` stands for a new unseeded `:jdk` generator. The shared [[default-rng]] is never used. Pass a seeded generator for reproducible noise.
  - `sigma` (finite number `>= 0`): the standard deviation of the noise. Defaults to `1.0`. `0` gives zeros.

  Returns a lazy infinite sequence of doubles. The sequence draws from the generator while it is consumed, so a generator should not be shared with other consumers if the stream matters.

  Throws `ex-info` (`:type :invalid-input`, `:property :sigma`) when `sigma` is not a finite number `>= 0`, and `ex-info` with `{:rng value}` when `rng` is neither a generator nor `nil` (a distribution is not accepted).

  See also [[ma]], [[ar]], [[arma]], [[fi]], [[arfima]], [[grandom]], [[->seq]]."
  ([] (white-noise nil 1.0))
  ([x] (if (number? x) (white-noise nil x) (white-noise x 1.0)))
  ([rng sigma]
   (let [g (ensure-rng rng)]
     (when-not (and (finite-number? sigma) (>= (double sigma) 0.0))
       (throw (invalid-input 'white-noise :sigma sigma "a finite number >= 0"
                             (str "The standard deviation must be a finite number >= 0, got " (pr-str sigma)))))
     (let [sigma (double sigma)]
       (repeatedly #(grandom g sigma))))))

;; arfima

(defn- series-options
  "Returns the options map of a time series function (`nil` = empty) or throws `ex-info` naming `function`."
  [function opts]
  (cond
    (nil? opts) {}
    (map? opts) opts
    :else (throw (invalid-input function :opts opts "a map or nil"
                                (str "The options must be a map, got " (pr-str opts))))))

(defn- burn-in-option
  "Returns the `:burn-in` option as a long, or `nil` when it is missing (the stage default is used). Throws `ex-info` when it is not an integer `>= 0`."
  [function opts]
  (let [k (get opts :burn-in)]
    (cond
      (nil? k) nil
      (and (long-integer? k) (>= (long k) 0)) (long k)
      :else (throw (invalid-input function :burn-in k "an integer >= 0"
                                  (str "The burn-in must be an integer >= 0, got " (pr-str k)))))))

(defn- array?
  [x]
  (and (some? x) (.isArray (class x))))

(defn- check-signal
  "Returns the signal or throws `ex-info` (`:property :signal`) naming `function` when it is not `nil`, a sequential collection or an array. The elements are not checked."
  [function signal]
  (if (or (nil? signal) (sequential? signal) (array? signal))
    signal
    (throw (invalid-input function :signal signal "a sequence of numbers"
                          (str "The signal must be a sequence of numbers, got " (pr-str signal))))))

(defn- default-signal
  "Standard normal white noise (or with `:sigma`) drawn from the generator `:rng` of the options (a new generator when missing)."
  [opts]
  (let [sigma (get opts :sigma)]
    (white-noise (ensure-rng (get opts :rng)) (if (nil? sigma) 1.0 sigma))))

(defn- series-coefficients
  "Returns the coefficients as a double array in the user order, `nil` is no coefficients. Throws `ex-info` (`:property` is `property`) naming `function` when they are not a number or a sequential collection or an array of finite numbers."
  ^doubles [function property coeffs]
  (let [expected "a finite number or a sequence of finite numbers"
        fail (fn [v]
               (invalid-input function property v expected
                              (str "The coefficients must be " expected ", got " (pr-str v))))
        xs (cond
             (nil? coeffs) ()
             (number? coeffs) (list coeffs)
             (or (sequential? coeffs) (array? coeffs)) (seq coeffs)
             :else (throw (fail coeffs)))]
    (doseq [x xs]
      (when-not (finite-number? x) (throw (fail x))))
    (double-array (map double xs))))

(defn- history-dot
  "Dot product of the history and the coefficients, summed from the first (oldest) element like `fastmath.vector/dot` of vectors does: `((h0 c0 + h1 c1) + h2 c2) + ...`. Both arrays have the same positive length."
  ^double [^doubles history ^doubles coeffs]
  (let [n (alength coeffs)]
    (loop [i 1
           acc (* (Array/aget history 0) (Array/aget coeffs 0))]
      (if (< i n)
        (recur (inc i) (+ acc (* (Array/aget history i) (Array/aget coeffs i))))
        acc))))

(defn- push-history
  "New history array: `history` without its oldest element and with `x` as the newest one."
  ^doubles [^doubles history ^double x]
  (let [n (alength history)
        nh (double-array n)]
    (System/arraycopy history 1 nh 0 (dec n))
    (Array/aset nh (dec n) x)
    nh))

(defn- ma-or-ar
  "Moving average (`ma?` true) or autoregressive filter of the signal. The first `cnt` (number of coefficients) elements of the signal seed the history, then `burn-in` outputs (`nil` is `cnt`) are dropped. The history has a fixed size and each element is computed from the previous history array (the memory does not grow). All-zero or no coefficients is the identity: the same elements are dropped."
  [function ma? coeffs signal burn-in]
  (let [^doubles user (series-coefficients function (if ma? :theta :phi) coeffs)
        cnt (long (alength user))
        burn (if (nil? burn-in) cnt (long burn-in))
        signal (check-signal function signal)]
    (if (every? zero? user)
      (drop (+ cnt burn) signal)
      (let [^doubles cs (double-array (reverse user))
            ;; cs[j] multiplies the history element j, the newest (index cnt-1) is multiplied by the first coefficient
            step (fn step [sig ^doubles history]
                   (when-let [s (seq sig)]
                     (let [e (double (first s))
                           a (+ (history-dot history cs) e)]
                       (cons a (lazy-seq (step (rest s) (push-history history (if ma? e a))))))))]
        (drop burn
              (lazy-seq
               (let [history (double-array cnt)]
                 (loop [i 0
                        s (seq signal)]
                   (if (< i cnt)
                     (when s
                       (Array/aset history i (double (first s)))
                       (recur (inc i) (next s)))
                     (step s history))))))))))

(defn ma
  "Generates a Moving Average (MA) stochastic process.

  An MA(q) process expresses the current value as a linear combination of the current and past white noise error terms `et`. A single `theta` coefficient produces an MA(1) process `Yt = et + theta1 et-1`; a sequence of q coefficients produces an MA(q) process `Yt = et + theta1 et-1 + theta2 et-2 + ... + thetaq et-q`. The sign convention is the one of R `arima.sim` (`+theta`); `fracdiff::fracdiff.sim` uses `-theta`.

  Parameters (`signal` and the options map can be given as the last arguments, an options map is told from a signal by `map?`):

  - `theta` (optional, number or sequence/array of numbers): the MA coefficient(s). A single number produces an MA(1) process; a sequence `[theta1 theta2 ... thetaq]` produces an MA(q) process, where `theta[i]` is the coefficient for `et-(i+1)`. Defaults to `0.0`. `nil` or `[]` is no coefficients, as all zeros it gives the white noise itself. The coefficients must be finite; there is no invertibility check.
  - `signal` (optional, sequence of numbers): the underlying white noise series (`et`). Defaults to standard normal random numbers, see [[white-noise]]. A number, a map or a string is not a signal; the elements are not checked before they are consumed.
  - `opts` (optional, map or `nil`, the last argument): `:rng` (random number generator or `nil`) and `:sigma` (standard deviation, default `1.0`) define the default signal `(white-noise (ensure-rng rng) sigma)` when no `signal` is given (they are ignored when it is), so `{:rng (rng :mersenne 1)}` makes the process reproducible and never uses [[default-rng]]. `:burn-in` (integer `>= 0`) is the number of outputs dropped after the history is seeded, the default is the number of coefficients `q` (the length of the output is then `N - 2q`). Unknown keys are ignored.

  Calls: `(ma)`, `(ma theta)`, `(ma theta signal)`, `(ma theta opts)`, `(ma theta signal opts)`.

  Returns a lazy sequence of the MA process. The first `q` elements of `signal` seed the history of past errors, the following `burn-in` outputs are dropped. For a signal of `N` elements the default length is `N - 2q` (empty when `N < 2q`; `N` elements less `q` seeds less `burn-in`), and the first output is `e[2q] + theta1 e[2q-1] + ... + thetaq e[q]`. Each output is the exact MA value because the history holds the true errors (the start is stationary; the `q` extra dropped outputs are not needed, `{:burn-in 0}` returns them). With all coefficients zero (or none) `(ma 0.0 s)` is `(drop 2 s)`: the identity, which still drops `2q` elements for `q` coefficients. The summation is oldest to newest, `theta` is applied as `theta1` to the newest past error. The memory used is constant for a lazily consumed signal.

  Throws `ex-info` (`:type :invalid-input`, `:property` `:theta`, `:signal`, `:opts` or `:burn-in`) for invalid coefficients (not a number or a collection of finite numbers), a signal which is not a sequence, options which are not a map, or a negative or non-integer `:burn-in`.

  See also [[ar]], [[arma]], [[fi]], [[arfima]], [[white-noise]]."
  ([] (ma 0.0))
  ([theta] (ma theta {}))
  ([theta signal-or-opts]
   (if (map? signal-or-opts)
     (ma theta (default-signal signal-or-opts) signal-or-opts)
     (ma theta signal-or-opts nil)))
  ([theta signal opts]
   (ma-or-ar 'ma true theta signal (burn-in-option 'ma (series-options 'ma opts)))))

(defn ar
  "Generates an Autoregressive (AR) stochastic process.

  An AR(p) process expresses the current value as a linear combination of the current white noise error term and past values of the process itself. A single `phi` coefficient produces an AR(1) process `Yt = phi1 Yt-1 + et`; a sequence of p coefficients produces an AR(p) process `Yt = phi1 Yt-1 + phi2 Yt-2 + ... + phip Yt-p + et`. The sign convention is the one of R `arima.sim` (`+phi`); `statsmodels` `ArmaProcess` uses the polynomial `[1, -phi1, ...]`.

  Parameters (`signal` and the options map can be given as the last arguments, an options map is told from a signal by `map?`):

  - `phi` (optional, number or sequence/array of numbers): the AR coefficient(s). A single number produces an AR(1) process; a sequence `[phi1 phi2 ... phip]` produces an AR(p) process, where `phi[i]` is the coefficient for `Yt-(i+1)`. Defaults to `0.0`. `nil` or `[]` is no coefficients, as all zeros it gives the white noise itself. The coefficients must be finite; there is no stationarity check, `(ar 1.2 signal)` runs and explodes.
  - `signal` (optional, sequence of numbers): the driving white noise series (`et`). Defaults to standard normal random numbers, see [[white-noise]]. A number, a map or a string is not a signal; the elements are not checked before they are consumed.
  - `opts` (optional, map or `nil`, the last argument): `:rng` (random number generator or `nil`) and `:sigma` (standard deviation, default `1.0`) define the default signal `(white-noise (ensure-rng rng) sigma)` when no `signal` is given (they are ignored when it is), so `{:rng (rng :mersenne 1)}` makes the process reproducible and never uses [[default-rng]]. `:burn-in` (integer `>= 0`) is the number of outputs dropped after the history is seeded, the default is the number of coefficients `p` (the length of the output is then `N - 2p`). Unknown keys are ignored.

  Calls: `(ar)`, `(ar phi)`, `(ar phi signal)`, `(ar phi opts)`, `(ar phi signal opts)`.

  Returns a lazy sequence of the AR process. The first `p` elements of `signal` seed the history as if they were `Y[0] ... Y[p-1]` (the raw errors, zero initial conditions for AR(1)), the following `burn-in` outputs are dropped. For a signal of `N` elements the default length is `N - 2p` (empty when `N < 2p`), and the first output is `Y[2p] = phi1 Y[2p-1] + ... + phip Y[p] + e[2p]`. `{:burn-in 0}` returns `N - p` elements starting with `Y[p]`. With all coefficients zero (or none) `(ar 0.0 s)` is `(drop 2 s)`, `(ar [0.0 0.0] s)` is `(drop 4 s)` and `(ar [] s)` is `s`.

  The start is not stationary: the first outputs have a smaller variance than the stationary one, for AR(1) the output `j` of the default sequence has the variance `sigma^2 (1 - phi^(2(j+3))) / (1 - phi^2)`, and for p >= 2 the seeding differs from the process with zero initial conditions by a transient decaying like the largest root modulus of the process. Drop the first elements of the default output, `(drop n (ar phi signal))` (or add `n` to `:burn-in`), to get close to the stationary process. The number of leading outputs to drop for AR(1) to reach 95% / 99% of the stationary variance is:

  | phi | 0.9 | 0.95 | 0.99 | 0.999 |
  |---|---|---|---|---|
  | 95% | 12 | 27 | 147 | 1495 |
  | 99% | 19 | 42 | 227 | 2299 |

  The summation is oldest to newest and the memory used is constant for a lazily consumed signal.

  Throws `ex-info` (`:type :invalid-input`, `:property` `:phi`, `:signal`, `:opts` or `:burn-in`) for invalid coefficients (not a number or a collection of finite numbers), a signal which is not a sequence, options which are not a map, or a negative or non-integer `:burn-in`.

  See also [[ma]], [[arma]], [[fi]], [[arfima]], [[white-noise]]."
  ([] (ar 0.0))
  ([phi] (ar phi {}))
  ([phi signal-or-opts]
   (if (map? signal-or-opts)
     (ar phi (default-signal signal-or-opts) signal-or-opts)
     (ar phi signal-or-opts nil)))
  ([phi signal opts]
   (ma-or-ar 'ar false phi signal (burn-in-option 'ar (series-options 'ar opts)))))

(defn arma
  "Generates an Autoregressive Moving Average (ARMA) stochastic process.

  An ARMA(p,q) process combines an MA(q) filter and an AR(p) filter: the driving `signal` is first passed through an MA(q) filter using `theta` (see [[ma]]), and the resulting series is then passed through an AR(p) filter using `phi` (see [[ar]]), giving `Yt = phi1 Yt-1 + ... + phip Yt-p + et + theta1 et-1 + ... + thetaq et-q`. The sign conventions are the ones of R `arima.sim` (`+phi`, `+theta`).

  Parameters (`signal` and the options map can be given as the last arguments, an options map is told from a signal by `map?`):

  - `phi` (optional, number or sequence/array of numbers): the AR coefficient(s), see [[ar]]. Defaults to `0.0`.
  - `theta` (optional, number or sequence/array of numbers): the MA coefficient(s), see [[ma]]. Defaults to `0.0`.
  - `signal` (optional, sequence of numbers): the underlying white noise series. Defaults to standard normal random numbers, see [[white-noise]].
  - `opts` (optional, map or `nil`, the last argument): `:rng` and `:sigma` define the default signal when no `signal` is given (see [[ma]]), `:burn-in` (integer `>= 0`) is passed to both stages: the number of outputs dropped after seeding each stage, the default is the number of coefficients of the stage. Unknown keys are ignored.

  Calls: `(arma)`, `(arma phi theta)`, `(arma phi theta signal)`, `(arma phi theta opts)`, `(arma phi theta signal opts)`.

  Returns a lazy sequence of the ARMA process. For a signal of `N` elements the default length is `N - 2q - 2p` (the sum of the lengths dropped by the two stages; `N - p - q - 2 burn-in` for an explicit `:burn-in`). The start is not stationary (the AR stage is seeded with the raw MA output, see [[ar]] for the burn-in needed), the difference to the zero-initial-condition process decays geometrically with the largest root modulus of the AR polynomial.

  Throws `ex-info` (`:type :invalid-input`) as [[ma]] and [[ar]], with the name `arma`.

  See also [[ma]], [[ar]], [[fi]], [[arfima]], [[white-noise]]."
  ([] (arma 0.0 0.0))
  ([phi theta] (arma phi theta {}))
  ([phi theta signal-or-opts]
   (if (map? signal-or-opts)
     (arma phi theta (default-signal signal-or-opts) signal-or-opts)
     (arma phi theta signal-or-opts nil)))
  ([phi theta signal opts]
   (let [burn-in (burn-in-option 'arma (series-options 'arma opts))]
     (ma-or-ar 'arma false phi (ma-or-ar 'arma true theta signal burn-in) burn-in))))

(defn- fi-diffs
  [^double d ^long limit]
  (-> (reduce (fn [buff ^long k]
                (conj buff (m/* (m/- (double (if (m/zero? k) -1.0 (buff (m/dec k)))))
                                (m// (m/- d k) (m/inc k))))) [] (range limit))
      (vec)))

(defn- cumulative-sum
  "Cumulative sums of the signal with the zero initial condition, as a lazy sequence of doubles."
  [signal]
  (rest (reductions (fn [acc e] (+ (double acc) (double e))) 0.0 signal)))

(defn- differences
  "First differences `e[t] - e[t-1]` of the signal with `e[-1] = 0` (the first element is `e[0]`), as a lazy sequence of doubles."
  [signal]
  (map (fn [a b] (- (double a) (double b))) signal (cons 0.0 signal)))

(defn- check-fi-order
  "Returns `d` as a double or throws `ex-info` (`:property :d`) naming `function` when it is not a finite number."
  ^double [function d]
  (if (finite-number? d)
    (double d)
    (throw (invalid-input function :d d "a finite number"
                          (str "The fractional differencing order d must be a finite number, got " (pr-str d))))))

(defn- check-dlimit
  "Returns `dlimit` as a long or throws `ex-info` (`:property :dlimit`) naming `function` when it is not an integer `>= 0`."
  ^long [function dlimit]
  (if (and (long-integer? dlimit) (>= (long dlimit) 0))
    (long dlimit)
    (throw (invalid-input function :dlimit dlimit "an integer >= 0"
                          (str "The dlimit must be an integer >= 0, got " (pr-str dlimit))))))

(defn- fi-stage
  "The fractional integration of the signal (`fi` without the default signal and the options): `burn-in` is `nil` (default) or the validated number of dropped outputs."
  [function d dlimit signal burn-in]
  (let [d (check-fi-order function d)
        dlimit (check-dlimit function dlimit)
        signal (check-signal function signal)]
    (cond
      (m/zero? d) signal
      (m/integer? d) (do
                       (when (> (Math/abs d) 100.0)
                         (throw (invalid-input function :d d "an integer with an absolute value <= 100 or a non-integer"
                                               (str "An integer fractional differencing order must have an absolute value <= 100, got " (pr-str d)))))
                       (let [k (Math/abs (long d))
                             f (if (pos? d) cumulative-sum differences)]
                         (drop (+ k (if (nil? burn-in) k (long burn-in)))
                               (nth (iterate f signal) k))))
      (pos? dlimit) (ma-or-ar function false (fi-diffs d dlimit) signal burn-in)
      (some? burn-in) (throw (invalid-input function :burn-in burn-in "nil (the exact fractional integration, non-integer d with dlimit 0, has no burn-in)"
                                            "The :burn-in option does not apply to the exact fractional integration (non-integer d, dlimit 0); use a positive dlimit"))
      :else (letfn [(step [sig history coeffs ^long k]
                      (when (seq sig)
                        (let [phi1 (double (if (m/zero? k) -1.0 (coeffs (m/dec k))))
                              phi2 (m/* (m/- phi1)
                                        (m// (m/- d k) (m/inc k)))
                              ncoeffs (conj coeffs phi2)
                              a (m/+ (v/dot history ncoeffs) (double (first sig)))
                              nhistory (conj history a)]
                          (cons a (lazy-seq (step (rest sig) nhistory ncoeffs (m/inc k)))))))]
              (lazy-seq
               (when-let [s (seq signal)]
                 (step (rest s) (list (first s)) [] 0)))))))

(defn fi
  "Generates a Fractionally Integrated (FI) process, or applies fractional integration to a series.

  Filters `signal` through the fractional differencing operator `(1-B)^-d`, where `B` is the backshift operator. Positive `d` adds long-range dependence (long memory) to the series, negative `d` produces anti-persistent, over-differenced behavior. Three cases are distinguished:

  - integer `d` (`1.0` is an integer too): `d >= 1` is the `d`-fold cumulative sum, `d <= -1` the `|d|`-fold first difference, both with zero initial conditions; `d = 0` is the identity;
  - non-integer `d` and `dlimit = 0` (the default): the exact fractional filter with all the weights `psi[k] = Gamma(k+d) / (Gamma(d) k!)` up to the current time, computed by Hosking's recursion;
  - non-integer `d` and `dlimit = L > 0`: the filter truncated to `L` lags, an AR(L) process with the coefficients of the series expansion of `(1-B)^d` (the same recursion).

  Parameters (`signal` and the options map can be given as the last arguments, an options map is told from a signal by `map?`):

  - `d` (optional, finite number): the fractional differencing order. Defaults to `0.0` (identity, `signal` is returned unchanged). For an integer `d` the absolute value must not exceed `100`. Values outside of `(-0.5, 0.5)` are accepted (non-stationary), there is no check.
  - `dlimit` (optional, integer `>= 0`): for non-integer `d`, truncates the fractional filter to this many lag coefficients. Ignored for an integer `d`. A call with two arguments is `d signal` (or `d opts`) with `dlimit` equal to `0`; with three arguments it is `d dlimit signal` or, when the third is a map, `d signal opts`.
  - `signal` (optional, sequence of numbers): the input series to integrate. Defaults to standard normal random numbers, see [[white-noise]].
  - `opts` (optional, map or `nil`, the last argument): `:rng` and `:sigma` define the default signal when no `signal` is given (see [[ma]]). `:burn-in` (integer `>= 0`) is the number of outputs dropped after seeding the history, the default is the number of coefficients (`L` for the truncated filter, `|d|` for an integer `d`: `|d|` seed elements and `|d|` outputs are dropped). It does not apply to the exact filter (non-integer `d`, `dlimit` 0): a given `:burn-in` throws there, use a positive `dlimit`. It is ignored for `d = 0`. Unknown keys are ignored.

  Calls: `(fi)`, `(fi d)`, `(fi d signal)`, `(fi d opts)`, `(fi d dlimit signal)`, `(fi d signal opts)`, `(fi d dlimit signal opts)`.

  Returns a lazy sequence. For a signal of `N` elements the length is:

  - exact filter: `N - 1`. The first element is only the seed and is not emitted, output `j` is `y[j+1]` with `y[t] = sum over k <= t of psi[k] e[t-k]` (zero initial conditions, so the variance converges to the stationary one slowly: for `d = 0.45` the sum of `psi^2` is 55% of its limit at `t = 100` and 72% at `t = 10^4`);
  - truncated filter: `N - 2L` (the seeding is raw as for [[ar]], with the same non-stationary start);
  - integer `d`: `N - 2|d|`, the elements `|d|` ... `2|d|-1` of the cumulative sum (or difference) are dropped; `(fi 1 s)` is `(drop 2 (reductions + s))`, `(fi -1 s)` is `(drop 2 (map - s (cons 0 s)))`;
  - `d = 0`: `N`, the signal itself.

  The exact filter takes `O(n^2)` time for `n` consumed elements (about 36 s for 30000 elements) and `O(n)` memory: pass a `dlimit` (a few hundred at most; the cost is `O(n L)`) for long series. The truncated process is not ARFIMA(0,d,0) for `d` near `0.5`: for `d = 0.45` and `L = 50` its variance is only 54% of the ARFIMA one (`L = 500`: 64%), for `d = 0.35` 91% (95%), for `d = 0.2` and `L = 500` 99.8%.

  Throws `ex-info` (`:type :invalid-input`, `:property` `:d`, `:dlimit`, `:signal`, `:opts` or `:burn-in`) when `d` is not a finite number, an integer `d` has an absolute value above `100`, `dlimit` is not an integer `>= 0`, the signal is not a sequence, the options are not a map, or `:burn-in` is invalid or given for the exact filter.

  See also [[ar]], [[ma]], [[arma]], [[arfima]]."
  ([] (fi 0.0))
  ([d] (fi d {}))
  ([d signal-or-opts]
   (if (map? signal-or-opts)
     (fi d 0 (default-signal signal-or-opts) signal-or-opts)
     (fi d 0 signal-or-opts nil)))
  ([d dlimit-or-signal signal-or-opts]
   (if (map? signal-or-opts)
     (fi d 0 dlimit-or-signal signal-or-opts)
     (fi d dlimit-or-signal signal-or-opts nil)))
  ([d dlimit signal opts]
   (fi-stage 'fi d dlimit signal (burn-in-option 'fi (series-options 'fi opts)))))

(defn arfima
  "Generates an Autoregressive Fractionally Integrated Moving Average (ARFIMA) stochastic process.

  An ARFIMA(p,d,q) process combines an MA(q) filter, fractional integration of order `d` and an AR(p) filter: the driving `signal` is first passed through an MA(q) filter using `theta` (see [[ma]]), then through the fractional integration filter using `d` and `dlimit` (see [[fi]]), and finally through an AR(p) filter using `phi` (see [[ar]]). The sign conventions are the ones of R `arima.sim` (`+phi`, `+theta`); `fracdiff::fracdiff.sim` has the opposite sign of `theta`.

  Parameters:

  - `phi` (optional, number or sequence/array of numbers): the AR coefficient(s), see [[ar]]. Defaults to `0`.
  - `d` (optional, finite number): the fractional differencing order, see [[fi]]. Defaults to `0`.
  - `dlimit` (optional, integer `>= 0`): truncation limit for the fractional integration step when `d` is not an integer, see [[fi]]. Defaults to `0` (the exact filter).
  - `theta` (optional, number or sequence/array of numbers): the MA coefficient(s), see [[ma]]. Defaults to `0`.
  - `signal` (optional, sequence of numbers): the underlying white noise series. Defaults to standard normal random numbers, see [[white-noise]].
  - `opts` (optional, map or `nil`, the last argument): `:rng` and `:sigma` define the default signal when no `signal` is given (see [[ma]]); `:burn-in` (integer `>= 0`) is passed to every stage (the number of outputs dropped after seeding each stage, the default is the number of coefficients of the stage); it throws for the exact `fi` stage (non-integer `d`, `dlimit` 0). Unknown keys are ignored.

  Calls (the forms differ by the number of arguments, and by the last argument being a map):

  - `(arfima)`, `(arfima phi d theta)` - the default signal;
  - `(arfima phi d theta opts)` - the default signal with options (a map as the 4th argument);
  - `(arfima phi d theta signal)` - the 4-argument form is `phi d theta signal`, without `dlimit` (so `(arfima 0.5 0.2 100 0.3)` is not a call with `dlimit`: use the 5-argument form);
  - `(arfima phi d dlimit theta signal)` - the 5-argument form, `dlimit` is the 3rd argument;
  - `(arfima phi d theta signal opts)` - the 5-argument form with a map as the last argument;
  - `(arfima phi d dlimit theta signal opts)` - all.

  Returns a lazy sequence of the ARFIMA process. For a signal of `N` elements the default length is `N - 2q - 1 - 2p` for the exact filter and `N - 2q - 2L - 2p` for a truncated one (the sum of the lengths dropped by the three stages, see [[ma]], [[fi]] and [[ar]]; for an integer `d` `N - 2q - 2|d| - 2p`). The start is not stationary (see [[ar]] and [[fi]]); the exact `fi` stage costs `O(n^2)`, use a `dlimit` for long series.

  Throws `ex-info` (`:type :invalid-input`) as [[ma]], [[fi]] and [[ar]], with the name `arfima`.

  See also [[ma]], [[ar]], [[fi]], [[arma]], [[white-noise]]."
  ([] (arfima 0 0 0))
  ([phi d theta] (arfima phi d theta {}))
  ([phi d theta signal-or-opts]
   (if (map? signal-or-opts)
     (arfima phi d 0 theta (default-signal signal-or-opts) signal-or-opts)
     (arfima phi d 0 theta signal-or-opts nil)))
  ([phi d x y z]
   (if (map? z)
     (arfima phi d 0 x y z)
     (arfima phi d x y z nil)))
  ([phi d dlimit theta signal opts]
   (let [burn-in (burn-in-option 'arfima (series-options 'arfima opts))]
     (ma-or-ar 'arfima false phi
               (fi-stage 'arfima d dlimit (ma-or-ar 'arfima true theta signal burn-in) burn-in)
               burn-in))))

(m/unuse-primitive-operators)
