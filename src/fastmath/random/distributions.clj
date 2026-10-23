(ns fastmath.random.distributions
  "Distribution objects of `fastmath.random`: the protocol implementations and the constructors behind the `fastmath.random/distribution` keys.

  The namespace extends the Apache Commons Math `RealDistribution` and `IntegerDistribution` and the multivariate normal distribution with the `fastmath.protocols` distribution protocols, and defines the distributions which have no backend or need a fix: the wrappers of SSJ continuous distributions (`ssj-continuous`), the distributions implemented here (zero-inflated and zero-adjusted families, generalized families, noncentral families, `von-mises`, `kolmogorov`, finite support discrete distributions, `multinomial`, `dirichlet`, `multi-normal`, `negative-binomial`, `logarithmic`, `beta-binomial` and others) and the combinators `truncated` and `mixture`. `->distribution` builds a distribution object from a map of functions, `integrate-pdf` and `icdf-solver` build a cdf and an inverse cdf from a density or a cdf, `ccdf` computes the upper tail probability.

  The functions are the building blocks for `fastmath.random`, which wraps them with the user-facing constructors (defaults, docstrings, the `distribution` multimethod) and the protocol proxies (`cdf`, `pdf`, `icdf` ... with the checks of the arguments). Positional arguments of the constructors here are not optional, `rng` can be `nil` (a fresh generator is created). Invalid arguments throw `ex-info` with `:type :invalid-input`."
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.special :as special]
            [fastmath.polynomials :as poly]
            [fastmath.kernel.density :as k]
            [fastmath.solver :as solver]
            
            [fastmath.protocols :as prot]

            [fastmath.interpolation.linear :as linear-interp]
            [fastmath.interpolation.cubic :as cubic-interp]
            [fastmath.interpolation.monotone :as monotone-interp]
            [fastmath.interpolation.step :as step-interp]
            
            [fastmath.calculus.quadrature :as quad]

            [clojure.data.int-map :as im])
  (:import [org.apache.commons.math3.random RandomGenerator EmpiricalDistribution JDKRandomGenerator MersenneTwister]
           [org.apache.commons.math3.distribution AbstractRealDistribution RealDistribution BetaDistribution CauchyDistribution ChiSquaredDistribution ConstantRealDistribution EnumeratedRealDistribution ExponentialDistribution FDistribution GammaDistribution, GumbelDistribution, LaplaceDistribution, LevyDistribution, LogisticDistribution, LogNormalDistribution, NakagamiDistribution, NormalDistribution, ParetoDistribution, TDistribution, TriangularDistribution, UniformRealDistribution WeibullDistribution MultivariateNormalDistribution]
           [org.apache.commons.math3.distribution IntegerDistribution AbstractIntegerDistribution BinomialDistribution EnumeratedIntegerDistribution, GeometricDistribution, HypergeometricDistribution, PascalDistribution, PoissonDistribution, UniformIntegerDistribution, ZipfDistribution]
           [org.apache.commons.math3.stat StatUtils]
           [org.apache.commons.math3.special Erf Gamma]

           [umontreal.ssj.probdist ContinuousDistribution FoldedNormalDist InverseGaussianDist]
           [umontreal.ssj.probdistmulti DirichletDist MultinomialDist]
           
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- invalid-input!
  "Throws `ex-info` for an invalid `property` of `function`."
  [function property message value expected]
  (throw (ex-info (str (name function) ", " (name property) ": " message)
                  {:type :invalid-input
                   :function function
                   :property property
                   :value value
                   :expected expected})))

(extend Object
  prot/DistributionIdProto
  {:distribution? (constantly false)})

(extend nil
  prot/DistributionIdProto
  {:distribution? (constantly false)})

(extend RealDistribution
  prot/DistributionProto
  {:cdf (fn (^double [^RealDistribution d ^double v] (.cumulativeProbability d v))
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

(defn- integer-argument?
  "True when `x` is an integer within the `int` range (the argument of the Commons Math integer distributions)."
  [^double x]
  (and (m/== x (m/floor x))
       (m/<= (m/- 2147483648.0) x 2147483647.0)))

(extend IntegerDistribution
  prot/DistributionProto
  {:cdf (fn
          (^double [^IntegerDistribution d ^double v] (.cumulativeProbability d (m/floor v)))
          (^double [^IntegerDistribution d ^double v1 ^double v2] (.cumulativeProbability d (m/floor v1) (m/floor v2))))
   :icdf (fn ^long [^IntegerDistribution d ^double p] (.inverseCumulativeProbability d p))
   :pdf (fn ^double [^IntegerDistribution d ^double p]
          (cond (m/nan? p) ##NaN
                (integer-argument? p) (.probability d (unchecked-int p))
                :else 0.0))
   :lpdf (fn ^double [^AbstractIntegerDistribution d ^double p]
           (cond (m/nan? p) ##NaN
                 (integer-argument? p) (.logProbability d (unchecked-int p))
                 :else ##-Inf))
   :probability (fn ^double [^IntegerDistribution d ^double p]
                  (cond (m/nan? p) ##NaN
                        (integer-argument? p) (.probability d (unchecked-int p))
                        :else 0.0))
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
   :continuous? (constantly false)})

;; multivariate normal: validation, point check and log-space density

(def ^:private ^:const multi-normal-symmetry-tolerance 1.0e-10)
(def ^:private ^:const multi-normal-positivity-tolerance 1.0e-13)

(defn- cholesky-factor
  "Lower triangular Cholesky factor `L` (`double[][]`, `L L^T = cov`) of the square matrix `cov` (`double[][]`).
  Returns `nil` when the matrix is not positive definite (a pivot is not above `1e-13` of its diagonal entry)."
  [^"[[D" cov]
  (let [n (alength cov)
        ^"[[D" l (into-array ^Class (Class/forName "[D") (repeatedly n #(double-array n)))]
    (loop [j 0]
      (if (m/== j n)
        l
        (let [^doubles lj (aget l j)
              ^doubles cj (aget cov j)
              cjj (aget cj j)
              s (double (loop [k 0
                               s cjj]
                          (if (m/== k j)
                            s
                            (recur (m/inc k) (m/- s (m/sq (aget lj k)))))))]
          (when (and (m/pos? cjj) (m/> s (m/* multi-normal-positivity-tolerance cjj)))
            (let [ljj (m/sqrt s)]
              (aset lj j ljj)
              (dotimes [i (m/- n j 1)]
                (let [i (m/+ j 1 i)
                      ^doubles li (aget l i)
                      t (double (loop [k 0
                                       t (aget ^doubles (aget cov i) j)]
                                  (if (m/== k j)
                                    t
                                    (recur (m/inc k) (m/- t (m/* (aget li k) (aget lj k)))))))]
                  (aset li j (m// t ljj))))
              (recur (m/inc j)))))))))

(defn- multi-normal-covariance!
  "Validates `covariances` (a non-empty square matrix of finite numbers, symmetric within `1e-10` relative to the largest entry, positive definite) and returns it as `double[][]`; throws `ex-info` otherwise."
  [covariances]
  (let [rows (when (and (some? covariances) (not (number? covariances)) (seqable? covariances))
               (try (mapv m/seq->double-array covariances) (catch Exception _ nil)))
        n (count rows)]
    (when (or (m/zero? n) (not-every? (fn [^doubles r] (m/== (alength r) n)) rows))
      (invalid-input! 'multi-normal :covariances "has to be a non-empty square matrix" covariances "n x n matrix"))
    (let [^"[[D" cov (into-array ^Class (Class/forName "[D") rows)]
      (when-not (every? (fn [^doubles r] (every? #(Double/isFinite (double %)) r)) rows)
        (invalid-input! 'multi-normal :covariances "entries have to be finite" covariances "finite numbers"))
      (let [mx (double (reduce (fn [^double a ^doubles r] (m/max a (double (reduce (fn [^double b ^double x] (m/max b (m/abs x))) 0.0 r)))) 0.0 rows))]
        (dotimes [i n]
          (dotimes [j i]
            (when (m/> (m/abs (m/- (aget ^doubles (aget cov i) j) (aget ^doubles (aget cov j) i)))
                       (m/* multi-normal-symmetry-tolerance mx))
              (invalid-input! 'multi-normal :covariances "has to be symmetric" covariances "symmetric matrix")))))
      (when-not (cholesky-factor cov)
        (invalid-input! 'multi-normal :covariances "has to be positive definite" covariances "positive definite matrix"))
      cov)))

(def ^:private ^java.util.Map multi-normal-cholesky-cache (java.util.Collections/synchronizedMap (java.util.WeakHashMap.)))

(defn- multi-normal-cholesky
  "Cholesky factor of the covariance of `d`, computed on the first request and kept (weakly) for the next ones."
  [^MultivariateNormalDistribution d]
  (or (.get multi-normal-cholesky-cache d)
      (when-let [l (cholesky-factor (.getData (.getCovariances d)))]
        (.put multi-normal-cholesky-cache d l)
        l)))

(defn- multi-normal-point
  "Converts the point `v` to `double[]`; throws `ex-info` for `nil`, non-sequences or a wrong length."
  ^doubles [^MultivariateNormalDistribution d v]
  (let [a (when (and (some? v) (not (number? v)) (seqable? v))
            (try (m/seq->double-array v) (catch Exception _ nil)))
        dim (.getDimension d)]
    (when (or (nil? a) (m/not== (alength ^doubles a) dim))
      (invalid-input! 'multi-normal :x "point has to be a sequence with one value per dimension" v dim))
    a))

(defn- multi-normal-lpdf
  "Log density of `d` at `x`: `log` of the Commons density when it is positive, the log-space formula (Cholesky based) when the density underflows to zero."
  ^double [^MultivariateNormalDistribution d ^doubles x]
  (let [density (.density d x)]
    (if (m/zero? density)
      (if-let [^"[[D" l (multi-normal-cholesky d)]
        (let [n (alength x)
              mu (.getMeans d)
              z (double-array n)]
          (dotimes [i n]
            (let [^doubles li (aget l i)
                  t (double (loop [k 0
                                   t (m/- (aget x i) (aget mu i))]
                              (if (m/== k i)
                                t
                                (recur (m/inc k) (m/- t (m/* (aget li k) (aget z k)))))))]
              (aset z i (m// t (aget li i)))))
          (let [maha (double (loop [i 0 s 0.0] (if (m/== i n) s (recur (m/inc i) (m/+ s (m/sq (aget z i)))))))
                logdet (double (loop [i 0 s 0.0] (if (m/== i n) s (recur (m/inc i) (m/+ s (m/* 2.0 (m/log (aget ^doubles (aget l i) i))))))))]
            (m/* -0.5 (m/+ (m/* n m/LOG_TWO_PI) logdet maha))))
        ##-Inf)
      (m/log density))))

(extend MultivariateNormalDistribution
  prot/DistributionProto
  {:pdf (fn ^double [^MultivariateNormalDistribution d v] (.density d (multi-normal-point d v)))
   :lpdf (fn ^double [^MultivariateNormalDistribution d v] (multi-normal-lpdf d (multi-normal-point d v)))
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

(defn multi-normal
  "Creates `MultivariateNormalDistribution` (Commons) after validation.

  * `means` - sequence of finite numbers, `nil`: zeros of the dimension of `covariances` (when both are `nil`: 2 dimensions),
  * `covariances` - square, symmetric (within `1e-10` relative to the largest entry), positive definite matrix of finite numbers,
    `nil`: identity matrix of the dimension of `means`,
  * `rng` - random generator, `nil`: a fresh one.

  Throws `ex-info` (`:type :invalid-input`) for invalid arguments and for the size mismatch of `means` and `covariances`."
  [means covariances rng]
  (let [means (when (some? means)
                (let [a (when (and (not (number? means)) (seqable? means))
                          (try (m/seq->double-array means) (catch Exception _ nil)))]
                  (when (or (nil? a) (m/zero? (alength ^doubles a)) (not (every? #(Double/isFinite (double %)) a)))
                    (invalid-input! 'multi-normal :means "has to be a non-empty sequence of finite numbers" means "sequence of numbers"))
                  a))
        ^"[[D" cov (cond
                     (some? covariances) (multi-normal-covariance! covariances)
                     (some? means) (let [n (alength ^doubles means)
                                         c (repeatedly n #(double-array n))]
                                     (dotimes [i n] (aset ^doubles (nth c i) i 1.0))
                                     (into-array ^Class (Class/forName "[D") c))
                     :else (multi-normal-covariance! [[1.0 0.0] [0.0 1.0]]))
        n (alength cov)
        ^doubles means (or means (double-array n))]
    (when (m/not== (alength means) n)
      (invalid-input! 'multi-normal :means "size does not match the dimension of covariances" (vec means) n))
    (try
      (MultivariateNormalDistribution. (or rng (JDKRandomGenerator.)) means cov)
      (catch RuntimeException e
        (throw (ex-info "multi-normal, covariances: the matrix is numerically singular or not positive definite"
                        {:type :invalid-input :function 'multi-normal :property :covariances
                         :value (mapv vec cov) :expected "positive definite matrix"}
                        e))))))

(extend-protocol prot/DistributionIdProto
  BetaDistribution
  (distribution? [_] true)
  (distribution-id [_] :beta)
  (distribution-parameters [_] [:alpha :beta :inverse-abs-accuracy :rng])
  CauchyDistribution
  (distribution? [_] true)
  (distribution-id [_] :cauchy)
  (distribution-parameters [_] [:median :scale :inverse-abs-accuracy :rng])
  ChiSquaredDistribution
  (distribution? [_] true)
  (distribution-id [_] :chi-squared)
  (distribution-parameters [_] [:degrees-of-freedom :inverse-abs-accuracy :rng])
  ConstantRealDistribution
  (distribution? [_] true)
  (distribution-id [_] :constant)
  (distribution-parameters [_] [:value])
  ExponentialDistribution
  (distribution? [_] true)
  (distribution-id [_] :exponential)
  (distribution-parameters [_] [:mean :inverse-abs-accuracy :rng])
  FDistribution
  (distribution? [_] true)
  (distribution-id [_] :f)
  (distribution-parameters [_] [:numerator-degrees-of-freedom :denominator-degrees-of-freedom :inverse-abs-accuracy :rng])
  GammaDistribution
  (distribution? [_] true)
  (distribution-id [_] :gamma)
  (distribution-parameters [_] [:shape :scale :inverse-abs-accuracy :rng])
  GumbelDistribution
  (distribution? [_] true)
  (distribution-id [_] :gumbel)
  (distribution-parameters [_] [:mu :beta :rng])
  LaplaceDistribution
  (distribution? [_] true)
  (distribution-id [_] :laplace)
  (distribution-parameters [_] [:mu :beta :rng])
  LevyDistribution
  (distribution? [_] true)
  (distribution-id [_] :levy)
  (distribution-parameters [_] [:mu :c :rng])
  LogisticDistribution
  (distribution? [_] true)
  (distribution-id [_] :logistic)
  (distribution-parameters [_] [:mu :s :rng])
  LogNormalDistribution
  (distribution? [_] true)
  (distribution-id [_] :log-normal)
  (distribution-parameters [_] [:scale :shape :inverse-abs-accuracy :rng])
  NakagamiDistribution
  (distribution? [_] true)
  (distribution-id [_] :nakagami)
  (distribution-parameters [_] [:mu :omega :inverse-abs-accuracy :rng])
  NormalDistribution
  (distribution? [_] true)
  (distribution-id [_] :normal)
  (distribution-parameters [_] [:mu :sd :inverse-abs-accuracy :rng])
  ParetoDistribution
  (distribution? [_] true)
  (distribution-id [_] :pareto)
  (distribution-parameters [_] [:scale :shape :inverse-abs-accuracy :rng])
  TDistribution
  (distribution? [_] true)
  (distribution-id [_] :t)
  (distribution-parameters [_] [:degrees-of-freedom :inverse-abs-accuracy :rng])
  TriangularDistribution
  (distribution? [_] true)
  (distribution-id [_] :triangular)
  (distribution-parameters [_] [:a :c :b :rng])
  UniformRealDistribution
  (distribution? [_] true)
  (distribution-id [_] :uniform-real)
  (distribution-parameters [_] [:lower :upper :rng])
  WeibullDistribution
  (distribution? [_] true)
  (distribution-id [_] :weibull)
  (distribution-parameters [_] [:alpha :beta :inverse-abs-accuracy :rng])

  EmpiricalDistribution
  (distribution? [_] true)
  (distribution-id [_] :empirical)
  (distribution-parameters [_] [:bin-count :data :rng])

  EnumeratedRealDistribution
  (distribution? [_] true)
  (distribution-id [_] :enumerated-real)
  (distribution-parameters [_] [:data :probabilities :rng])
  EnumeratedIntegerDistribution
  (distribution? [_] true)
  (distribution-id [_] :enumerated-int)
  (distribution-parameters [_] [:data :probabilities :rng])

  BinomialDistribution
  (distribution? [_] true)
  (distribution-id [_] :binomial)
  (distribution-parameters [_] [:trials :p :rng])
  GeometricDistribution
  (distribution? [_] true)
  (distribution-id [_] :geometric)
  (distribution-parameters [_] [:p :rng])
  HypergeometricDistribution
  (distribution? [_] true)
  (distribution-id [_] :hypergeometric)
  (distribution-parameters [_] [:population-size :number-of-successes :sample-size :rng])
  PascalDistribution
  (distribution? [_] true)
  (distribution-id [_] :pascal)
  (distribution-parameters [_] [:r :p :rng])
  PoissonDistribution
  (distribution? [_] true)
  (distribution-id [_] :poisson)
  (distribution-parameters [_] [:p :epsilon :max-iterations :rng])
  UniformIntegerDistribution
  (distribution? [_] true)
  (distribution-id [_] :uniform-int)
  (distribution-parameters [_] [:lower :upper :rng])
  ZipfDistribution
  (distribution? [_] true)
  (distribution-id [_] :zipf)
  (distribution-parameters [_] [:number-of-elements :exponent :rng])

  MultivariateNormalDistribution
  (distribution? [_] true)
  (distribution-id [_] :multi-normal)
  (distribution-parameters [_] [:means :covariances :rng]))

;;

(defn ssj-continuous
  "Wraps the SSJ continuous distribution `d` into a distribution object.

  Parameters:

  - `nm`: keyword, the value of `distribution-id`
  - `d`: an SSJ `ContinuousDistribution`, it is also the `source-object`
  - `rng`: random generator (`nil`: a fresh one)
  - `params`: vector of keywords reported by `distribution-parameters`

  `pdf`, `cdf`, `icdf`, `mean`, `variance` and the bounds are the ones of SSJ, `lpdf` is the logarithm of the density. A sample is the SSJ quantile of a uniform number from `rng`. `grandom`, `brandom` and `set-seed` throw `UnsupportedOperationException`, `set-seed!` reseeds `rng`. Returns a univariate continuous distribution."
  [nm ^ContinuousDistribution d ^RandomGenerator rng params]
  (let [rng (or rng (JDKRandomGenerator.))]
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
      (distribution-parameters [_] params)
      prot/UnivariateDistributionProto
      (mean [_] (.getMean d))
      (variance [_] (.getVariance d))
      (lower-bound [_] (.getXinf d))
      (upper-bound [_] (.getXsup d))
      prot/RNGProto
      (grandom [_] (throw (java.lang.UnsupportedOperationException. "Gaussian random is not supported.")))
      (brandom [_] (throw (java.lang.UnsupportedOperationException. "Boolean random is not supported.")))
      (set-seed [_ _] (throw (java.lang.UnsupportedOperationException. "Immutable seeding is not supported, use `set-seed!` instead.")))
      (drandom [_] (.inverseF d (prot/drandom rng)))
      (frandom [_] (unchecked-float (.inverseF d (prot/drandom rng))))
      (lrandom [_] (m/round-even (.inverseF d (prot/drandom rng))))
      (irandom [_] (unchecked-int (m/round-even (.inverseF d (prot/drandom rng)))))
      (->seq [_] (repeatedly #(.inverseF d (prot/drandom rng))))
      (->seq [_ n] (repeatedly n #(.inverseF d (prot/drandom rng))))
      (set-seed! [d seed] (prot/set-seed! rng seed) d))))

(defn ssj-continuous-no-pdf
  "Wraps the SSJ continuous distribution `d` which has no density into a distribution object.

  It is `ssj-continuous` with the `pdf` replaced by the central difference of the cdf, `(cdf(x + 5e-7) - cdf(x - 5e-7)) / 1e-6`, so the density is approximate. Parameters and returned value as for `ssj-continuous`."
  [nm ^ContinuousDistribution d ^RandomGenerator rng params]
  (let [rng (or rng (JDKRandomGenerator.))]
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
      (distribution-parameters [_] params)
      prot/UnivariateDistributionProto
      (mean [_] (.getMean d))
      (variance [_] (.getVariance d))
      (lower-bound [_] (.getXinf d))
      (upper-bound [_] (.getXsup d))
      prot/RNGProto
      (grandom [_] (throw (java.lang.UnsupportedOperationException. "Gaussian random is not supported.")))
      (brandom [_] (throw (java.lang.UnsupportedOperationException. "Boolean random is not supported.")))
      (set-seed [_ _] (throw (java.lang.UnsupportedOperationException. "Immutable seeding is not supported, use `set-seed!` instead.")))
      (drandom [_] (.inverseF d (prot/drandom rng)))
      (frandom [_] (unchecked-float (.inverseF d (prot/drandom rng))))
      (lrandom [_] (m/round-even (.inverseF d (prot/drandom rng))))
      (irandom [_] (unchecked-int (m/round-even (.inverseF d (prot/drandom rng)))))
      (->seq [_] (repeatedly #(.inverseF d (prot/drandom rng))))
      (->seq [_ n] (repeatedly n #(.inverseF d (prot/drandom rng))))
      (set-seed! [d seed] (prot/set-seed! rng seed) d))))

(defn- log-erfcx
  "`log(erfcx(b)) = b^2 + log(erfc(b))` for `b >= 0`: direct for `b < 5`, the continued fraction
  `erfcx(b) = 1/(sqrt(pi) (b + (1/2)/(b + 1/(b + (3/2)/(b + ...)))))` (60 levels, error below 1e-16 there) for larger `b`
  (which neither overflows nor underflows)."
  ^double [^double b]
  (if (m/< b 5.0)
    (m/+ (m/* b b) (m/log (special/erfc b)))
    (let [t (loop [n 60
                   t b]
              (if (m/zero? n)
                t
                (recur (m/dec n) (m/+ b (m// (m/* 0.5 n) t)))))]
      (m/- (m/- (m/* 0.5 m/LOG_PI)) (m/log t)))))

(defn- inverse-gaussian-cdf
  "cdf of the inverse Gaussian distribution, `Phi(z) + exp(2 lambda/mu) Phi(-w)` with the second term computed in log-space:
  `exp(2 lambda/mu) Phi(-w) = exp(-(lambda/(2 mu)) (x/mu - 2 + mu/x)) erfcx(w/sqrt(2)) / 2`. Does not overflow for any `lambda/mu`."
  ^double [^double mu ^double lambda ^double x]
  (cond
    (m/not-pos? x) 0.0
    (m/pos-inf? x) 1.0
    :else (let [xt (m// x mu)
                s (m/sqrt (m// lambda x))
                t1 (m/* 0.5 (special/erfc (m// (m/- (m/* s (m/dec xt))) m/SQRT2)))
                b (m// (m/* s (m/inc xt)) m/SQRT2)
                lt2 (m/- (m/+ (m/- (m/* 0.5 (m// lambda mu) (m/+ (m/- xt 2.0) (m// 1.0 xt)))) (log-erfcx b)) m/LN2)]
            (m/min 1.0 (m/+ t1 (m/exp lt2))))))

(defn- inverse-gaussian-cdf-log
  "`cdf(exp(t)) - p` of the inverse Gaussian distribution, the function solved by [[inverse-gaussian-icdf]]."
  ^double [^double mu ^double lambda ^double p ^double t]
  (m/- (inverse-gaussian-cdf mu lambda (m/exp t)) p))

(defn- inverse-gaussian-icdf
  "Quantile of the inverse Gaussian distribution by the bracketed Brent root finding of the cdf in `log x`."
  ^double [^double mu ^double lambda ^double p]
  (cond
    (m/not-pos? p) 0.0
    (m/>= p 1.0) ##Inf
    :else (let [h (fn [^double t] (inverse-gaussian-cdf-log mu lambda p t))
                t0 (m/log mu)
                step (m// 1.0 (m/sqrt (m// lambda mu)))
                h0 (inverse-gaussian-cdf-log mu lambda p t0)]
            (cond
              (m/zero? h0) mu
              (m/neg? h0) (let [hi (loop [t (m/+ t0 step) k 2.0]
                                     (if (or (m/>= (inverse-gaussian-cdf-log mu lambda p t) 0.0) (m/> t 709.0))
                                       t
                                       (recur (m/+ t0 (m/* k step)) (m/* 2.0 k))))]
                            (m/exp (solver/find-root h t0 hi {:absolute-accuracy 1.0e-14 :max-iters 500})))
              :else (let [lo (loop [t (m/- t0 step) k 2.0]
                               (if (or (m/<= (inverse-gaussian-cdf-log mu lambda p t) 0.0) (m/< t -744.0))
                                 t
                                 (recur (m/- t0 (m/* k step)) (m/* 2.0 k))))]
                      (if (m/> (inverse-gaussian-cdf-log mu lambda p lo) 0.0)
                        0.0
                        (m/exp (solver/find-root h lo t0 {:absolute-accuracy 1.0e-14 :max-iters 500}))))))))

(defn inverse-gaussian
  "Inverse Gaussian (Wald) distribution `(mu, lambda)`.

  For `lambda/mu <= 100` the SSJ `InverseGaussianDist` is used for everything (same as [[ssj-continuous]]). For larger `lambda/mu`
  (the SSJ cdf overflows from about 355 and `icdf`/samples degenerate) `cdf` and `icdf` are computed by [[inverse-gaussian-cdf]]
  (log-space second term) and a bracketed root finding in `log x`; the samples are `icdf` of uniform numbers.
  `pdf`, moments and bounds come from SSJ."
  [^double mu ^double lambda rng]
  (let [rng (or rng (JDKRandomGenerator.))
        d (InverseGaussianDist. mu lambda)
        params [:mu :lambda :rng]]
    (if (m/<= (m// lambda mu) 100.0)
      (ssj-continuous :inverse-gaussian d rng params)
      (reify
        prot/DistributionProto
        (pdf [_ v] (.density d v))
        (lpdf [_ v] (m/log (.density d v)))
        (probability [_ v] (.density d v))
        (cdf [_ v] (inverse-gaussian-cdf mu lambda v))
        (cdf [_ v1 v2] (m/- (inverse-gaussian-cdf mu lambda v2) (inverse-gaussian-cdf mu lambda v1)))
        (icdf [_ v] (inverse-gaussian-icdf mu lambda v))
        (sample [_] (inverse-gaussian-icdf mu lambda (prot/drandom rng)))
        (dimensions [_] 1)
        (source-object [_] d)
        (continuous? [_] true)
        prot/DistributionIdProto
        (distribution? [_] true)
        (distribution-id [_] :inverse-gaussian)
        (distribution-parameters [_] params)
        prot/UnivariateDistributionProto
        (mean [_] (.getMean d))
        (variance [_] (.getVariance d))
        (lower-bound [_] (.getXinf d))
        (upper-bound [_] (.getXsup d))
        prot/RNGProto
        (grandom [_] (throw (java.lang.UnsupportedOperationException. "Gaussian random is not supported.")))
        (brandom [_] (throw (java.lang.UnsupportedOperationException. "Boolean random is not supported.")))
        (set-seed [_ _] (throw (java.lang.UnsupportedOperationException. "Immutable seeding is not supported, use `set-seed!` instead.")))
        (drandom [_] (inverse-gaussian-icdf mu lambda (prot/drandom rng)))
        (frandom [_] (unchecked-float (inverse-gaussian-icdf mu lambda (prot/drandom rng))))
        (lrandom [_] (m/round-even (inverse-gaussian-icdf mu lambda (prot/drandom rng))))
        (irandom [_] (unchecked-int (m/round-even (inverse-gaussian-icdf mu lambda (prot/drandom rng)))))
        (->seq [_] (repeatedly #(inverse-gaussian-icdf mu lambda (prot/drandom rng))))
        (->seq [_ n] (repeatedly n #(inverse-gaussian-icdf mu lambda (prot/drandom rng))))
        (set-seed! [d seed] (prot/set-seed! rng seed) d)))))

(defn- multinomial-point
  "Converts the point `v` to `double[]` of the length `dim`; throws `ex-info` for `nil`, non-sequences or a wrong length."
  ^doubles [v ^long dim]
  (let [a (when (and (some? v) (not (number? v)) (seqable? v))
            (try (m/seq->double-array v) (catch Exception _ nil)))]
    (when (or (nil? a) (m/not== (alength ^doubles a) dim))
      (invalid-input! 'multinomial :x "point has to be a sequence with one count per category" v dim))
    a))

(defn multinomial
  "Multinomial distribution of `n` (integer, >= 0) trials over the categories with the (relative) probabilities `ps`.

  The probability mass function uses `0 * log 0 = 0` (a category with the probability zero and the count zero does not contribute),
  returns 0 for negative or non-integer counts, counts with a wrong sum and a positive count in a category with the probability zero."
  [n ps binomial rng]
  (let [rng (or rng (JDKRandomGenerator.))
        _ (when-not (and (number? n) (m/not-neg? (double n)) (m/== (double n) (m/floor (double n))) (m/<= (double n) Integer/MAX_VALUE))
            (invalid-input! 'multinomial :n "has to be an integer from the range [0, 2^31)" n "integer >= 0"))
        n (long n)
        _ (when-not (and (some? ps) (not (number? ps)) (seqable? ps) (seq ps)
                         (every? (fn [p] (and (number? p) (Double/isFinite (double p)) (m/not-neg? (double p)))) ps)
                         (m/pos? (double (reduce m/+ 0.0 (map double ps)))))
            (invalid-input! 'multinomial :ps "has to be a non-empty sequence of finite, non-negative numbers with a positive sum" ps "probabilities"))
        ps (m/seq->double-array (v/div ps (v/sum ps)))
        dim (alength ps)
        ;; categories with the probability zero never occur and are removed for SSJ (0*log 0 = NaN there)
        positive (filterv (fn [^long i] (m/pos? (Array/aget ps i))) (range dim))
        ps+ (m/seq->double-array (map (fn [^long i] (Array/aget ps i)) positive))
        n-int (int n)
        ;; returns `nil` when the point is not a possible outcome, NaN when it contains NaN,
        ;; otherwise the counts of the categories with a positive probability (int[])
        outcome (fn [v]
                  (let [^doubles x (multinomial-point v dim)]
                    (cond
                      (some m/nan? x) ##NaN
                      (not-every? (fn [^double c] (and (m/>= c 0.0) (m/<= c n) (m/== c (m/floor c)))) x) nil
                      (m/not== (long (reduce m/+ 0.0 x)) n) nil
                      (not-every? (fn [^long i] (or (m/pos? (Array/aget ps i)) (m/zero? (Array/aget x i)))) (range dim)) nil
                      :else (int-array (map (fn [^long i] (Array/aget x i)) positive)))))
        pmf (fn ^double [v]
              (let [o (outcome v)]
                (cond
                  (nil? o) 0.0
                  (number? o) ##NaN
                  (m/zero? n) 1.0
                  :else (MultinomialDist/prob n-int ps+ ^ints o))))
        lpmf (fn ^double [v]
               (let [o (outcome v)]
                 (cond
                   (nil? o) ##-Inf
                   (number? o) ##NaN
                   (m/zero? n) 0.0
                   :else (let [p (MultinomialDist/prob n-int ps+ ^ints o)]
                           (if (m/pos? p)
                             (m/log p)
                             ;; underflow: log-space formula
                             (reduce m/+ (m/log-factorial n)
                                     (map (fn [^long c ^double pi] (if (m/zero? c) 0.0 (m/- (m/* c (m/log pi)) (m/log-factorial c))))
                                          o ps+)))))))
        cdf (fn ^double [v]
              (let [^doubles x (multinomial-point v dim)]
                (cond
                  (some m/nan? x) ##NaN
                  (some m/neg? x) 0.0
                  (m/zero? n) 1.0
                  :else (MultinomialDist/cdf n-int ps+ (int-array (map (fn [^long i] (m/min (Array/aget x i) (double n))) positive))))))
        m (delay (if (m/zero? n) (vec (repeat dim 0.0)) (vec (MultinomialDist/getMean n-int ps))))
        cv (delay (if (m/zero? n) (vec (repeat dim (vec (repeat dim 0.0)))) (mapv vec (MultinomialDist/getCovariance n-int ps))))
        binom-probs (mapv (fn [^double prob ^double sum]
                            (let [q (m// prob (m/- 1.0 sum))]
                              (if (m/nan? q) 0.0 q))) ps (reductions m/+ 0.0 ps))]
    (reify
      prot/DistributionProto
      (pdf [_ v] (pmf v))
      (lpdf [_ v] (lpmf v))
      (probability [_ v] (pmf v))
      (cdf [_ v] (cdf v))
      (sample [_] (first (reduce (fn [[buf ^int curr] ^double prob]
                                   (let [res (long (prot/sample (binomial {:trials curr :p (m/constrain prob 0.0 1.0) :rng rng})))]
                                     [(conj buf res) (m/- curr res)])) [[] n] binom-probs)))
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
      (set-seed! [d seed] (prot/set-seed! rng seed) d))))

;;

(def ^{:const true :private true :tag 'double} zero+epsilon (m/next-double 0.0))
(def ^{:const true :private true :tag 'double} one-epsilon (m/prev-double 1.0))

(defn- dirichlet-rev-log-beta
  ^double [alpha]
  (let [d (special/log-gamma (reduce m/+ alpha))
        ^double n (reduce m/+ (map (fn [^double v] (special/log-gamma v)) alpha))]
    (m/- d n)))

(def ^:private ^:const dirichlet-simplex-tolerance 1.0e-8)

(defn- dirichlet-lpdf
  "Log density of the Dirichlet distribution at the point `values`.

  The boundary (`x_i = 0`) is handled by the limits: `alpha_i < 1` gives `+Inf`, `alpha_i = 1` a finite value (the term is dropped),
  `alpha_i > 1` gives `-Inf` (which wins when both kinds of terms are present). A point with NaN gives NaN. Throws `ex-info` for a point
  of a wrong length, `nil`, with a negative coordinate or with the sum farther than `1e-8` from 1."
  ^double [^long dim alpha- values ^double lbeta]
  (let [x (when (and (some? values) (not (number? values)) (seqable? values))
            (try (m/seq->double-array values) (catch Exception _ nil)))]
    (when (or (nil? x) (m/not== (alength ^doubles x) dim))
      (invalid-input! 'dirichlet :x "point has to be a sequence with one coordinate per dimension" values dim))
    (if (some m/nan? x)
      ##NaN
      (do
        (when (some m/neg? x)
          (invalid-input! 'dirichlet :x "coordinates have to be non-negative" values "x_i >= 0"))
        (when (m/> (m/abs (m/- (double (reduce m/+ 0.0 x)) 1.0)) dirichlet-simplex-tolerance)
          (invalid-input! 'dirichlet :x "coordinates have to sum to 1 (tolerance 1e-8)" values "sum(x) = 1"))
        (let [v (reduce m/+ lbeta (map (fn [^double ai ^double xi]
                                         (if (m/zero? xi)
                                           (cond (m/zero? ai) 0.0
                                                 (m/neg? ai) ##Inf
                                                 :else ##-Inf)
                                           (m/* ai (m/log xi)))) alpha- x))]
          (if (m/nan? v) ##-Inf v))))))

(defn- dirichlet-log-space-draw
  "Dirichlet draw in log-space: `G_a = G_{a+1} * U^(1/a)` for the gamma variates, normalised with log-sum-exp.
  Used when the plain gamma draws are too small to be normalised."
  [alpha gamma+ rng]
  (let [lg (mapv (fn [^double a g]
                   (let [u (m/max Double/MIN_VALUE (double (prot/drandom rng)))]
                     (m/+ (m/log (double (prot/sample g))) (m// (m/log u) a))))
                 alpha gamma+)
        lse (m/logsumexp lg)]
    (mapv (fn [^double l] (m/exp (m/- l lse))) lg)))

(defn dirichlet
  "Dirichlet distribution on the simplex.

  Parameters:

  - `alpha`: sequence of positive concentration parameters, an integer `n` means `n` parameters equal to `1.0`; the dimension is their count
  - `gamma`: the constructor of the gamma distribution (a function of a map with `:shape`, `:scale` and `:rng`), used for sampling
  - `rng`: random generator (`nil`: a fresh one)

  A sample is a vector of normalised gamma variates (computed in log space when their sum is below `1e-8`, so small `alpha` do not degenerate); zero and one coordinates are moved to the neighbouring double. `pdf`, `lpdf` and `probability` take a sequence with one coordinate per dimension, non-negative and summing to 1 within `1e-8`, and use the limits at the boundary (a coordinate zero with `alpha < 1` gives `+Inf`); other points throw `ex-info`. `means` and `covariance` are available, `cdf`, `icdf`, `mean`, `variance` and bounds are not. Returns a multivariate continuous distribution."
  [alpha gamma rng]
  (let [rng (or rng (JDKRandomGenerator.))
        alpha (if (integer? alpha)
                (double-array alpha 1.0)
                (m/seq->double-array alpha))
        sampler (mapv (fn [^double shape] (gamma {:shape shape :scale 1.0 :rng rng})) alpha)
        sampler+ (delay (mapv (fn [^double shape] (gamma {:shape (m/inc shape) :scale 1.0 :rng rng})) alpha))

        lbeta (dirichlet-rev-log-beta alpha)
        alpha- (mapv m/dec alpha)
        
        m (delay (vec (DirichletDist/getMean alpha)))
        cv (delay (mapv vec (DirichletDist/getCovariance alpha)))
        dim (count alpha)]
    (reify
      prot/DistributionProto
      (pdf [_ v] (m/exp (dirichlet-lpdf dim alpha- v lbeta)))
      (lpdf [_ v] (dirichlet-lpdf dim alpha- v lbeta))
      (probability [_ v] (m/exp (dirichlet-lpdf dim alpha- v lbeta)))
      (sample [_] (let [samples (map prot/sample sampler)
                        s (v/sum samples)]
                    (mapv (fn [^double v] (cond
                                           (m/zero? v) zero+epsilon
                                           (m/== v 1.0) one-epsilon
                                           :else v))
                          (if (m/> s 1.0e-8)
                            (v/div samples s)
                            (dirichlet-log-space-draw alpha @sampler+ rng)))))
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
      (set-seed! [d seed] (prot/set-seed! rng seed) d))))

;;

(defn- diff-cdf
  ^double [cdf ^double v1 ^double v2]
  (m/- (double (cdf v2)) (double (cdf v1))))

(def ^{:const true :private true :tag 'double} long-range-limit 9.223372036854775807E18)

(defn- integer-valued-pmf
  "Wraps the function `f` (called with a `long`) of a discrete pmf (`value-outside` = 0.0 for `pdf`, -Inf for `lpdf`).
  Returns a function of a double argument: NaN gives NaN, a non-integer or a value out of the `long` range gives `value-outside`,
  an integer value is passed to `f` unchanged (the results at integers are those of `f`)."
  [f ^double value-outside]
  (fn ^double [^double x]
    (cond
      (m/nan? x) ##NaN
      (m/not== x (m/floor x)) value-outside
      (m/>= (m/abs x) long-range-limit) value-outside
      :else (double (f (unchecked-long x))))))

(defn- integer-valued-cdf
  "Wraps the function `f` (called with a `long`) of a discrete cdf: the argument is floored, the result is 0.0 below `lower-bound`
  (or when the floor is below the `long` range), 1.0 for the values out of the `long` range above, NaN for NaN."
  [f ^double lower-bound]
  (fn ^double [^double x]
    (cond
      (m/nan? x) ##NaN
      (m/< x lower-bound) 0.0
      (m/<= x (m/- long-range-limit)) 0.0
      (m/>= x long-range-limit) 1.0
      :else (double (f (unchecked-long (m/floor x)))))))

(defn ->distribution
  "Creates a distribution object from the given functions (`:pdf`, `:lpdf`, `:cdf`, `:icdf`, `:sampler`, ...).

  Parameters (a map):

  - `pdf`, `lpdf`: density (mass) and its logarithm, functions of a double (of a long for discrete distributions on the integers); one of them is enough, the other is derived
  - `cdf`: cumulative distribution function
  - `icdf`: quantile function; it is also the default sampler (applied to a uniform number from `rng`)
  - `sampler`: function without arguments returning a sample, or `:long` for a quantile based sampler returning longs
  - `seeder`: function of a seed called by `set-seed!`, default reseeds `rng`
  - `rng`: random generator (`nil`: a fresh one)
  - `dimensions` (long, default 1), `continuous?`, `integer-support?` (discrete distributions: `false` for the real valued ones)
  - `name`, `parameters`: values of `distribution-id` and `distribution-parameters`
  - `mean`, `variance`: numbers or delays (evaluated on first request)
  - `lower-bound`, `upper-bound`: bounds of the support

  The result is a univariate distribution (`dimensions` is only reported); `source-object` is the object itself, `grandom`, `brandom` and `set-seed` throw `UnsupportedOperationException`.

  For discrete distributions on the integers (`:continuous?` false, unless `:integer-support?` is false) the `:pdf`, `:lpdf` and `:cdf`
  functions are called with integer (`long`) arguments only: `pdf` and `probability` return 0.0 and `lpdf` -Inf for a non-integer
  argument, `cdf` floors its argument and is 0.0 below `:lower-bound`, NaN is passed through as NaN. Discrete distributions on
  arbitrary real numbers set `:integer-support?` to false and receive the argument unchanged."
  [{:keys [pdf lpdf cdf icdf rng ^long dimensions continuous? integer-support? name parameters mean variance lower-bound upper-bound sampler seeder]
    :or {dimensions 1}}]
  (let [rng (or rng (JDKRandomGenerator.))
        pdf (or pdf (fn ^double [^double v] (m/exp (double (lpdf v)))))
        lpdf (or lpdf (fn ^double [^double v] (m/log (double (pdf v)))))
        [pdf lpdf cdf] (if (if (nil? integer-support?) continuous? (not integer-support?))
                         [pdf lpdf cdf]
                         [(integer-valued-pmf pdf 0.0)
                          (integer-valued-pmf lpdf ##-Inf)
                          (integer-valued-cdf cdf (double (or lower-bound ##-Inf)))])
        sampler (if (= sampler :long)
                  (fn ^long [] (icdf (prot/drandom rng)))
                  (or sampler (fn ^double [] (icdf (prot/drandom rng)))))
        seeder (if seeder
                 (fn [^long seed] (seeder seed))
                 (fn [^long seed] (prot/set-seed! rng seed)))]
    (reify
      prot/DistributionProto
      (pdf [_ v] (pdf v))
      (lpdf [_ v] (lpdf v))
      (cdf [_ v] (cdf v))
      (cdf [_ v1 v2] (diff-cdf cdf v1 v2))
      (icdf [_ v] (icdf v))
      (probability [_ v] (pdf v))
      (sample [_] (sampler))
      (dimensions [_] dimensions)
      (source-object [d] d)
      (continuous? [_] continuous?)
      prot/DistributionIdProto
      (distribution? [_] true)
      (distribution-id [_] name)
      (distribution-parameters [_] parameters)
      prot/UnivariateDistributionProto
      (mean [_] (if (delay? mean) @mean mean))
      (variance [_] (if (delay? variance) @variance variance))
      (lower-bound [_] lower-bound)
      (upper-bound [_] upper-bound)
      prot/RNGProto
      (drandom [_] (sampler))
      (frandom [_] (unchecked-float (sampler)))
      (lrandom [_] (m/round-even (sampler)))
      (irandom [_] (unchecked-int (m/round-even (sampler))))
      (->seq [_] (repeatedly sampler))
      (->seq [_ n] (repeatedly n sampler))
      (set-seed! [d seed] (seeder seed) d)
      (grandom [_] (throw (java.lang.UnsupportedOperationException. "Gaussian random is not supported.")))
      (brandom [_] (throw (java.lang.UnsupportedOperationException. "Boolean random is not supported.")))
      (set-seed [_ _] (throw (java.lang.UnsupportedOperationException. "Immutable seeding is not supported, use `set-seed!` instead."))))))

;;

(defn- strictly-increasing!
  "Returns the cumulative values as a `double[]` in which every value is below its successor.

  Increments below half an ulp of the values near 1.0 (the zero tail of a pdf with a bounded support) are lost in the cumulative
  sum, which leaves equal values at the end and makes the inverse interpolation `NaN`. Going from the end, every value which is
  not below its successor is replaced by the double just below it, so the last value (1.0) is kept."
  ^doubles [ys]
  (let [^doubles ys (double-array ys)]
    (loop [i (m/long-sub (alength ys) 2)]
      (when (m/>= i 0)
        (let [next (Array/aget ys (m/long-inc i))]
          (when (m/>= (Array/aget ys i) next)
            (Array/aset ys i (m/prev-double next))))
        (recur (m/long-dec i))))
    ys))

(defn integrate-pdf
  "Integrates the density `pdf-func` numerically and returns the cdf and the inverse cdf as functions.

  Parameters:

  - `pdf-func`: density, a function of a double
  - `mn`: lower bound of the integration (default `0.0`), the density should be `0.0` there
  - `mx`: upper bound of the integration (default `1.0`)
  - `steps`: number of subintervals (default `1000`); each is integrated with the Gauss-Kronrod quadrature
  - `interpolator`: interpolation of the integrated values, `:monotone` (default), `:linear`, `:cubic` or an interpolation function of `fastmath.interpolation` (a function of `xs` and `ys`)

  The second arity takes `pdf-func` and a map with the keys `:mn`, `:mx`, `:steps` and `:interpolator`; other keys are passed to the quadrature. The range is extended by five subintervals on both sides; the values are normalised so the cdf is `1.0` at the end of the extended range; it is `0.0` below and `1.0` above that range.

  Returns a vector `[cdf icdf]` of two functions of a double. `icdf` expects a probability from `[0, 1]`, both are interpolations, so they are approximate."
  ([pdf-func mn mx steps]
   (integrate-pdf pdf-func {:mn mn :mx mx :steps steps}))
  ([pdf-func {:keys [^double mn ^double mx ^long steps interpolator]
              :or {mn 0.0 mx 1.0 steps 1000 interpolator :monotone}
              :as options}]
   (let [diff5 (* 5.0 (m// (m/- mx mn) steps))
         mn (m/- mn diff5)
         mx (m/+ mx diff5)
         xs (m/slice-range mn mx steps)
         f (fn [^double x] (if (m/<= mn x mx) (pdf-func x) 0.0))
         int-options (assoc options :info? false)
         ys (->> (partition 2 1 xs)
                 (map (fn [[^double x1 ^double x2]]
                        (m/max m/MACHINE-EPSILON ;; in case if integration is zero
                               ^double (quad/gk-quadrature f x1 x2 int-options))))
                 (reductions m/+ 0.0)
                 (m/seq->double-array))
         ys (strictly-increasing! (v/div ys (Array/aget ys (dec steps)))) ;; normalize to ensure 1 at the endpoint
         intpol (case interpolator
                  :monotone monotone-interp/monotone
                  :linear linear-interp/linear
                  :cubic cubic-interp/cubic
                  (if (fn? interpolator) interpolator linear-interp/linear))]
     [(let [i (intpol xs ys)]
        ;; outside of the integrated range the interpolator extrapolates, the cdf is 0 below and 1 above
        (fn [^double x] (cond
                          (m/<= x mn) 0.0
                          (m/>= x mx) 1.0
                          :else (m/constrain (double (i x)) 0.0 1.0))))
      (intpol ys xs)])))

(defn- find-first-non-zero
  ^double [f xs]
  (or (->> xs (filter (comp m/pos? f)) (first))
      (first xs)))

(defn- narrow-range
  [kd [^double mn ^double mx ^double step] ^long steps]
  [(m/- (find-first-non-zero kd (m/slice-range mn mx steps)) step)
   (m/+ (find-first-non-zero kd (m/slice-range mx mn steps)) step)])

(defn continuous-distribution
  "Nonparametric continuous distribution built from the data by a kernel density estimate.

  Parameters (a map):

  - `data`: sequence of numbers
  - `kde`: kernel density estimation method, see `fastmath.kernel.density/kernel-density+`
  - `bandwidth`: bandwidth of the estimate, `nil` selects the default one
  - `steps` (long): number of subintervals of the numerical integration of the density, from which `cdf` and `icdf` are interpolated
  - `interpolator`: interpolation of the integrated density, see `integrate-pdf`
  - `rng`: random generator (`nil`: a fresh one)

  `pdf` is the estimate, `cdf` and `icdf` come from `integrate-pdf` on the range of the data extended where the estimate is positive, `mean` and `variance` are the sample mean and the sample variance of `data`, not the moments of the estimate. Returns a univariate continuous distribution."
  [{:keys [data ^long steps kde bandwidth rng]
    :as all}]
  (let [{:keys [kde ^double mn ^double mx]} (k/kernel-density+ kde data {:bandwidth bandwidth})
        step (m// (m/- mx mn) steps)
        [^double mn ^double mx] (narrow-range kde [mn mx step] (m/long-mult 4 steps))
        [cdf-fn icdf-fn] (integrate-pdf kde (merge all {:mn mn :mx mx :steps steps}))
        m (delay (StatUtils/mean (m/seq->double-array data)))
        v (delay (StatUtils/variance (m/seq->double-array data) (double @m)))]
    (->distribution {:pdf kde
                     :cdf (fn ^double [^double v] (m/constrain (double (cdf-fn v)) 0.0 1.0))
                     :icdf (fn ^double [^double v] (icdf-fn (m/constrain v 0.0 1.0)))
                     :rng rng
                     :dimensions 1
                     :continuous? true
                     :name :continuous-distribution
                     :parameters [:data :steps :kde :bandwidth :interpolator :rng]
                     :mean m
                     :variance v
                     :lower-bound (m/- mn step)
                     :upper-bound (m/+ mx step)})))

;;

(defn discrete-binary-search
  "Smallest integer `k` with `cdf-fn(k) >= p`, found by an exponential and a binary search.

  Parameters:

  - `cdf-fn`: cdf of an integer-valued distribution, a function of a long
  - `p`: probability
  - `[mid step]` (3 arguments): the start, the search begins on `[0, mid]` and extends the right end with the initial `step` (at least 1) while `cdf-fn` is below `p`
  - `step`, `[mn mx]` (4 arguments): the current step and bracket

  Returns a long."
  ([cdf-fn ^double p [^long mid ^long step]] (discrete-binary-search cdf-fn (m/long-max 1 step) p [0 mid]))
  ([cdf-fn ^long step ^double p [^long mn ^long mx]]
   (cond
     (m/< (double (cdf-fn mx)) p) (recur cdf-fn (m/* 2 step) p [mx (m/+ mx step)])
     (m/<= (m/- mx mn) 1) (if (m/>= (double (cdf-fn mn)) p) mn mx)
     :else (let [mid (m// (m/+ mn mx) 2)]
             (if (m/> (double (cdf-fn mid)) p)
               (recur cdf-fn step p [mn mid])
               (recur cdf-fn step p [mid mx]))))))

(defn- smallest-with-cdf
  "Moves the candidate `k` (an approximate quantile) to the smallest integer `x >= lower` with `(cdf x) >= p` by checking the `cdf` itself.
  `p` has to be below 1. The search is limited to 1000 steps in each direction."
  ^long [cdf ^double p ^long lower ^long k]
  (loop [k (m/max lower k)
         steps 0]
    (cond
      (m/>= steps 1000) k
      (m/< (double (cdf k)) p) (recur (m/inc k) (m/inc steps))
      (and (m/> k lower) (m/>= (double (cdf (m/dec k))) p)) (recur (m/dec k) (m/inc steps))
      :else k)))

(declare negative-binomial-core)

(defn negative-binomial
  "Negative binomial distribution `(r, p)`: the number of failures before `r` successes with the success probability `p`
  (`r > 0`, `0 < p <= 1`). The optional `parameters` (default `[:r :p :rng]`) is what `distribution-parameters` reports,
  for the constructors accepting other parameters (`nbi`, `nbii`)."
  ([^double r ^double p rng] (negative-binomial r p rng [:r :p :rng]))
  ([^double r ^double p rng parameters]
   (when-not (and (m/pos? r) (m/< r ##Inf))
     (invalid-input! 'negative-binomial :r "has to be positive and finite" r "r > 0"))
   (when-not (and (m/pos? p) (m/<= p 1.0))
     (invalid-input! 'negative-binomial :p "has to be from the range (0, 1]" p "0 < p <= 1"))
   (negative-binomial-core r p rng parameters)))

(defn- negative-binomial-core
  [^double r ^double p rng parameters]
  (let [p- (m/- 1.0 p)
        mean (m// (m/* r p-) p)
        variance (m// mean p)
        lgr (special/log-gamma r)
        lp- (m/log (m/- 1.0 p))
        lpr (m/* r (m/log p))
        cdf (fn ^double [^double k]
              (if (m/neg? k)
                0.0
                (special/regularized-beta p r (m/inc (m/rint k)))))]
    (->distribution {:lpdf (fn [^long k]
                             (if (m/neg? k)
                               ##-Inf
                               (m/+ (m/- (special/log-gamma (m/+ r k))
                                         (m/+ (m/log-factorial k) lgr))
                                    (m/* k lp-) lpr)))
                     :cdf cdf
                     :icdf (fn [^double p]
                             (cond
                               (m/not-pos? p) 0
                               (m/>= p 1.0) ##Inf
                               :else (discrete-binary-search cdf p [(long mean) (long (m/sqrt variance))])))
                     :sampler :long
                     :rng rng
                     :dimensions 1
                     :continuous? false
                     :name :negative-binomial
                     :parameters parameters
                     :mean mean
                     :variance variance
                     :lower-bound 0
                     :upper-bound Integer/MAX_VALUE})))

(def ^{:const true :private true :tag 'double} logarithmic-plateau 1.0e-17)

(defn logarithmic
  "Logarithmic (log-series) distribution, `P(k) = -p^k / (k log(1-p))`, `k >= 1`, `0 < p < 1`.

  The cdf is the running sum of the pmf terms `p^j/j` (direct, no cancellation) normalised by `-log(1-p)`; it is 1.0 as soon as the
  remaining tail is below `1e-17`. `icdf` scans the same running sum, which makes it the exact smallest `k` with `cdf(k) >= p`
  and gives a fast inversion sampler (the cost is proportional to the quantile, i.e. around the mean `-p/((1-p)log(1-p))` per sample).
  The cost of `cdf` is proportional to `min(k, 40/(1-p))`."
  [^double p rng]
  (when-not (and (m/< 0.0 p) (m/< p 1.0))
    (invalid-input! 'logarithmic :p "has to be from the range (0, 1)" p "0 < p < 1"))
  (let [p- (m/- 1.0 p)
        logp- (m/log p-)
        r (m// logp-)
        nr (m/- r)
        logp (m/log p)
        mean (m/* -1.0 r (m// p p-))
        variance (m/- (m// (m/+ (m/* p p) (m/* p logp-))
                           (m/sq (m/* p- logp-))))
        ;; cumulative probability after the term `j` (`pw` = p^j, `s` = sum of the terms up to j)
        cum-at (fn [^double s ^double pw ^long j]
                 (if (m/< (m/* nr (m// (m/* pw p) (m/* (m/inc j) p-))) logarithmic-plateau)
                   1.0
                   (m/min 1.0 (m/* s nr))))
        ;; p^j: multiplication, refreshed every 32 steps to bound the rounding error
        next-pw (fn [^double pw ^long j]
                  (if (m/zero? (m/bit-and j 31)) (m/fpow p j) (m/* pw p)))
        cdf (fn ^double [^long k]
              (if (m/< k 1)
                0.0
                (loop [j 1
                       pw p
                       s p]
                  (let [c (double (cum-at s pw j))]
                    (if (or (m/== j k) (m/== c 1.0))
                      c
                      (let [j1 (m/inc j)
                            pw1 (double (next-pw pw j1))]
                        (recur j1 pw1 (m/+ s (m// pw1 j1)))))))))
        ;; fpow handles exponents within the int range only
        pdf (fn [^long k]
              (if (m/< k 1)
                0.0
                (m/* -1.0 r (m// (if (m/<= k Integer/MAX_VALUE) (m/fpow p k) (m/exp (m/* k logp))) k))))]
    (->distribution {:pdf pdf
                     :lpdf (fn ^double [^long k]
                             (if (m/< k 1)
                               ##-Inf
                               (let [pd (double (pdf k))]
                                 (if (m/pos? pd)
                                   (m/log pd)
                                   (m/- (m/+ (m/log nr) (m/* k logp)) (m/log k))))))
                     :cdf cdf
                     :icdf (fn [^double q]
                             (cond
                               (m/not-pos? q) 1
                               (m/>= q 1.0) ##Inf
                               :else (loop [j 1
                                            pw p
                                            s p]
                                       (if (m/>= (double (cum-at s pw j)) q)
                                         j
                                         (let [j1 (m/inc j)
                                               pw1 (double (next-pw pw j1))]
                                           (recur j1 pw1 (m/+ s (m// pw1 j1))))))))
                     :sampler :long
                     :rng rng
                     :dimensions 1
                     :continuous? false
                     :name :logarithmic
                     :parameters [:p :rng]
                     :mean mean
                     :variance variance
                     :lower-bound 1
                     :upper-bound Integer/MAX_VALUE})))


(def ^{:const true :private true :tag 'double} LOG_M_2_PI (m/log m/M_2_PI))

(defn half-cauchy
  "Half-Cauchy distribution with the location `mu` and the `scale`: the Cauchy density folded at `mu`, `2 / (pi scale (1 + ((x - mu) / scale)^2))` for `x >= mu`.

  `rng` is the random generator (`nil`: a fresh one). `mean` and `variance` are `NaN`. The reported bounds are `0.0` and `+Inf` whatever `mu` is. Returns a univariate continuous distribution."
  [^double mu ^double scale rng]
  (let [ls (m/log scale)]
    (->distribution {:lpdf (fn [^double x]
                             (if (m/< x mu )
                               ##-Inf
                               (m/- LOG_M_2_PI ls (m/log1p (m/sq (m// (m/- x mu) scale))))))
                     :cdf (fn [^double v] (if (m/< v mu)
                                           0.0
                                           (m/* m/M_2_PI (m/atan (m// (m/- v mu) scale)))))
                     :icdf (fn [^double p]
                             (cond
                               (m/not-pos? p) mu
                               (m/>= p 1.0) ##Inf
                               :else (m/+ mu (m/* scale (m/tan (m/* m/HALF_PI p))))))
                     :rng rng
                     :dimensions 1
                     :continuous? true
                     :name :half-cauchy
                     :parameters [:mu :scale :rng]
                     :mean ##NaN
                     :variance ##NaN
                     :lower-bound 0.0
                     :upper-bound ##Inf})))

;;

(defn- inverse-cdf-nodes
  "Returns `[cumulative-probabilities keys]` without the repeated cumulative probabilities (zero probability keys, or probabilities lost in the rounding of the cumulative sum): the first key of every run of equal values stays, which is the smallest key with that cumulative probability. The step interpolator needs strictly increasing nodes."
  [cumsum ks]
  (let [pairs (->> (map vector cumsum ks)
                   (partition-by first)
                   (map first))]
    [(map first pairs) (map second pairs)]))

(defn- discrete-input!
  "Validates the `data` and `probabilities` of a finite support distribution of `function`.
  Returns `[data probabilities]` (vectors, probabilities as doubles; missing `probabilities` are equal).

  Throws `ex-info` for: `data` not being a non-empty sequence, `probabilities` of a different count than `data`, a negative,
  NaN or infinite probability and for the zero sum of the probabilities. With `data-type` `:integer` the data have to be
  finite integers within the `long` range (returned as `long`s), with `:real` numbers other than NaN."
  [function data probabilities data-type]
  (let [sequence-like? (fn [x] (and (some? x) (not (number? x)) (not (map? x)) (seqable? x)))
        data' (when (sequence-like? data) (vec data))]
    (when (empty? data')
      (invalid-input! function :data "has to be a non-empty sequence" data "non-empty sequence"))
    (case data-type
      :integer (when-not (every? (fn [x] (and (number? x)
                                              (if (integer? x)
                                                (<= Long/MIN_VALUE x Long/MAX_VALUE)
                                                (and (Double/isFinite (double x))
                                                     (m/== (double x) (m/floor (double x)))
                                                     (m/< (m/abs (double x)) long-range-limit)))))
                                 data')
                 (invalid-input! function :data "has to contain integers only" data "integers within the long range"))
      :real (when-not (every? (fn [x] (and (number? x) (not (m/nan? (double x))))) data')
              (invalid-input! function :data "has to contain numbers other than NaN" data "numbers"))
      nil)
    (let [probs (if (nil? probabilities)
                  (vec (repeat (count data') 1.0))
                  (when (sequence-like? probabilities) (vec probabilities)))]
      (when (nil? probs)
        (invalid-input! function :probabilities "has to be a sequence of numbers or nil" probabilities "sequence of numbers"))
      (when (m/not== (count probs) (count data'))
        (invalid-input! function :probabilities (str "has to have the same count as data (" (count data') ")") probabilities (count data')))
      (when-not (every? (fn [p] (and (number? p) (Double/isFinite (double p)) (m/not-neg? (double p)))) probs)
        (invalid-input! function :probabilities "has to contain finite, non-negative numbers only" probabilities "finite numbers >= 0"))
      (let [probs (mapv double probs)]
        (when-not (m/pos? (double (reduce m/+ 0.0 probs)))
          (invalid-input! function :probabilities "has to have a positive sum" probabilities "sum > 0"))
        [(if (= data-type :integer) (mapv long data') data')
         probs]))))

(defn integer-discrete-distribution
  "Discrete distribution on a finite set of integers.

  Parameters:

  - `data`: non-empty sequence of integers within the long range; repeated values are summed
  - `probabilities`: sequence of non-negative finite weights of the same count, normalised internally; `nil`: equal weights
  - `rng`: random generator (`nil`: a fresh one)

  `icdf` is the smallest value with `cdf >= p`, `cdf` and the pmf accept any number (a non-integer `pdf` argument gives `0.0`), samples are `Long`s. Values of `data` with an absolute value above `2^53` are held as the nearest double (`pdf` there is `0.0`, `icdf` returns the rounded value). Bounds are the smallest and the largest value. Throws `ex-info` for invalid `data` or `probabilities` (empty data, a different count, negative, `NaN` or infinite weights, a zero sum). Returns a univariate discrete distribution."
  [data probabilities rng]
  (let [[data probabilities] (discrete-input! 'integer-discrete-distribution data probabilities :integer)
        sum (v/sum probabilities)
        pmf (->> (map vector data probabilities)
                 (reduce (fn [m [v ^double p]]
                           (im/update m v (fnil m/+ 0.0) (m// p sum))) (im/int-map)))
        probs (vals pmf)
        cumsum (reductions m/+ probs)
        ks (keys pmf)
        mnk (double (first ks))
        [icdf-cumsum icdf-ks] (inverse-cdf-nodes cumsum ks)
        step-before (step-interp/step-before icdf-cumsum icdf-ks)
        step-after (step-interp/step-after ks cumsum)
        mean (delay (v/dot ks probs))
        icdf (fn ^long [^double x] (step-before x))
        rng (or rng (JDKRandomGenerator.))]
    (->distribution {:pdf (fn ^double [^long k] (get pmf k 0.0))
                     :icdf icdf
                     :cdf (fn ^double [^double x] (if (m/< x mnk) 0.0 (step-after x)))
                     :sampler :long
                     :rng rng
                     :dimensions 1
                     :continuous? false
                     :name :integer-discrete-distribution
                     :parameters [:data :probabilities :rng]
                     :mean mean
                     :variance (delay (m/- (v/dot (v/sq ks) probs)
                                           (m/sq @mean)))
                     :lower-bound mnk
                     :upper-bound (last ks)})))


(defn real-discrete-distribution
  "Discrete distribution on a finite set of real numbers.

  Parameters, validation and the quantile as for `integer-discrete-distribution`, but `data` are any numbers other than `NaN`. `pdf` needs an exact value of the data (anything else is `0.0`), `cdf` is the step function on the reals. Samples are doubles. Returns a univariate discrete distribution."
  [data probabilities rng]
  (let [[data probabilities] (discrete-input! 'real-discrete-distribution data probabilities :real)
        sum (v/sum probabilities)
        pmf (->> (map vector data probabilities)
                 (reduce (fn [m [v ^double p]]
                           (update m v (fnil m/+ 0.0) (m// p sum))) (sorted-map)))
        probs (vals pmf)
        cumsum (reductions m/+ probs)
        ks (keys pmf)
        mnk (double (first ks))
        [icdf-cumsum icdf-ks] (inverse-cdf-nodes cumsum ks)
        step-before (step-interp/step-before icdf-cumsum icdf-ks)
        step-after (step-interp/step-after ks cumsum)
        mean (delay (v/dot ks probs))]
    (->distribution {:pdf (fn ^double [^double k] (get pmf k 0.0))
                     :icdf (fn ^double [^double x] (step-before x))
                     :cdf (fn ^double [^double x] (if (m/< x mnk) 0.0 (step-after x)))
                     :rng rng
                     :dimensions 1
                     :continuous? false
                     :integer-support? false
                     :name :real-discrete-distribution
                     :parameters [:data :probabilities :rng]
                     :mean mean
                     :variance (delay (m/- (v/dot (v/sq ks) probs)
                                           (m/sq @mean)))
                     :lower-bound mnk
                     :upper-bound (last ks)})))

(defn categorical-distribution
  "Discrete distribution on an arbitrary finite set of values (keywords, strings, any objects compared by equality).

  Parameters:

  - `data`: non-empty sequence of values; repeated values are summed
  - `probabilities`: non-negative finite weights of the same count, normalised internally; `nil`: equal weights
  - `rng`: random generator (`nil`: a fresh one)

  `pdf`, `cdf` and `probability` take a value of the data (a value outside of it gives `0.0`), `cdf` is the cumulative probability in the order in which the distinct values first appear in `data`, `icdf` returns the value at which that cumulative probability reaches `p` and `sample` a value. The `source-object` is the integer-discrete distribution over the positions. `mean` and `variance` are `NaN`, the bounds are not supported. Throws `ex-info` as `integer-discrete-distribution`."
  [data probabilities rng]
  (let [rng (or rng (JDKRandomGenerator.))
        [data probabilities] (discrete-input! 'categorical-distribution data probabilities nil)
        ^clojure.lang.ILookup unique (vec (distinct data))
        ^clojure.lang.ILookup dict (zipmap unique (range (count unique)))

        enumerated (integer-discrete-distribution (map dict data) probabilities rng)]
    (reify
      prot/DistributionProto
      (pdf [_ v] (prot/pdf enumerated (.valAt dict v -1)))
      (lpdf [_ v] (prot/lpdf enumerated (.valAt dict v -1)))
      (cdf [_ v] (prot/cdf enumerated (.valAt dict v -1)))
      (cdf [d v1 v2] (m/- (double (prot/cdf d v2))
                          (double (prot/cdf d v1))))
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
      (set-seed! [d seed] (prot/set-seed! rng seed) d))))

;;

(defn icdf-solver
  "Solves `cdf(x) = p` for `x`, for a cdf which is not bounded on at least one side.

  Parameters:

  - `cdf`: cdf function of a double
  - `p`: target probability
  - `init`: the starting point
  - `step`: the first step of the search for a bracket of the root

  The bracket is found by moving away from `init` with steps growing geometrically, then the root is found by a bracketing solver with the absolute accuracy `1e-10`. Returns the root as a double, `##Inf` or `##-Inf` when the bracket leaves the range of doubles."
  ^double [cdf ^double p ^double init ^double step]
  (let [h1 (fn ^double [^double q] (m/- (double (cdf q)) p))]
    (if (m/< (double (cdf init)) p)
      ;; bracket offset grows geometrically (step, 2*step, 4*step, ...)
      ;; rather than linearly, so heavy-tailed cdfs whose true root is many
      ;; orders of magnitude past `init` (e.g. f-noncentral with small df2)
      ;; are bracketed in a few dozen iterations instead of needing a
      ;; number of iterations proportional to the root's own magnitude.
      (loop [interval (m/+ init step)
             mult 2.0]
        (cond
          (m/pos-inf? interval) ##Inf ;; bracket grew past double range (or cdf
          ;; can never numerically reach p below Infinity, e.g. a quadrature-based
          ;; cdf that plateaus short of 1.0 for all finite x): find-root can't take
          ;; an infinite bound, and the true root is unresolvable anyway, so return
          ;; the boundary directly instead of throwing.
          (m/< (double (cdf interval)) p) (recur (m/+ init (m/* mult step)) (m/* mult 2.0))
          :else (solver/find-root h1 init interval {:absolute-accuracy 1.0e-10})))
      (loop [interval (m/- init step)
             mult 2.0]
        (cond
          (m/neg-inf? interval) ##-Inf
          (m/> (double (cdf interval)) p) (recur (m/- init (m/* mult step)) (m/* mult 2.0))
          :else (solver/find-root h1 interval init {:absolute-accuracy 1.0e-10}))))))

;;

(defn- kolmogorov-pdf
  ^double [^double x]
  (cond
    (m/not-pos? x) 0.0
    (m/<= x 1.0) (let [c (m// m/PI (m/* 2.0 x))
                       ks (map (fn [^long i]
                                 (let [k (m/sq (m/* i c))]
                                   (m/* (m/dec k) (m/exp (m/* -0.5 k)))))
                               (range 1 40 2))]
                   (m// (m/* m/SQRT2PI (v/sum ks)) (m/* x x)))
    :else (let [ks (map (fn [^double a ^long i]
                          (m/* a i i (m/exp (m/* -2.0 (m/sq (m/* i x))))))
                        (cycle [1.0 -1.0]) (range 1 21))]
            (* 8.0 x (v/sum ks)))))

(defn- kolmogorov-cdf-raw
  ^double [^double x]
  (let [a (m/- (m/sq (m// m/PI x)))
        f (m/exp a)
        f2 (m/* f f)
        u (m/inc (m/* f (m/inc f2)))]
    (m// (m/* m/SQRT2PI (m/exp (m/* 0.125 a)) u) x)))

(defn- kolmogorov-ccdf-raw
  ^double [^double x]
  (let [f (m/exp (m/* -2.0 x x))
        f2 (m/* f f)
        f3 (m/* f f2)
        f5 (m/* f2 f3)
        f7 (m/* f2 f5)
        u (m/- 1.0 (m/* f3 (m/- 1.0 (m/* f5 (m/- 1.0 f7)))))]
    (m/* 2.0 f u)))

(defn- kolmogorov-cdf
  ^double [^double x]
  (cond
    (m/not-pos? x) 0.0
    (m/<= x 1.0) (kolmogorov-cdf-raw x)
    :else (m/- 1.0 (kolmogorov-ccdf-raw x))))

(defn- kolmogorov-icdf
  ^double [^double p]
  (icdf-solver kolmogorov-cdf p 0.8687311606361591 0.26033287146241274))

(defn kolmogorov
  "Kolmogorov distribution, the limit of the scaled Kolmogorov-Smirnov statistic `sqrt(n) D_n` for large `n`.

  `rng` is the random generator (`nil`: a fresh one). Continuous on `[0, Inf)`, `mean` is `0.8687311606361591` and `variance` is `0.0677732039638651`. Returns a univariate continuous distribution."
  [rng]
  (->distribution {:pdf kolmogorov-pdf
                   :icdf kolmogorov-icdf
                   :cdf kolmogorov-cdf
                   :rng rng
                   :dimensions 1
                   :continuous? true
                   :name :kolmogorov
                   :parameters [:rng]
                   :mean 0.8687311606361591
                   :variance 0.0677732039638651
                   :lower-bound 0.0
                   :upper-bound ##Inf}))
;;

(defn- fnh-mode
  ^long [^double omega ^long n ^long ns ^long nf]
  (let [A (m/dec omega)
        B (m/- n nf (m/* (m/+ ns n 2) omega))
        C (m/* (inc ns) (inc n) omega)]
    (long (m/floor (m// (m/* -2.0 C)
                        (m/- B (m/sqrt (m/- (m/* B B) (m/* 4.0 A C)))))))))

(defn fishers-noncentral-hypergeometric
  "Fisher's noncentral hypergeometric distribution.

  The number of items of the first kind in a sample of size `n` from an urn with `ns` items of the first and `nf` of the second kind, when the odds of choosing the first kind are `omega` times the odds of the second (`omega = 1` is the ordinary hypergeometric distribution).

  Parameters (a map): `ns`, `nf`, `n` (longs), `omega` (positive double), `rng` (`nil`: a fresh generator). The support is `[max(0, n - nf), min(ns, n)]`; the probabilities are computed for the whole support at construction, so it should be of moderate size. Returns a univariate discrete distribution."
  [{:keys [^long ns ^long nf ^long n ^double omega rng]}]
  (let [rng (or rng (JDKRandomGenerator.))
        lower-bound (m/max 0 (m/- n nf))
        upper-bound (m/min ns n)
        lower-bound (m/min lower-bound upper-bound)
        upper-bound (m/max lower-bound upper-bound)
        mode (fnh-mode omega n ns nf)
        fri (fn ^double [^long i]
              (m/* (m// (m/* (m/inc (m/- ns i)) omega)
                        (m/* i (m/+ (m/- nf n) i)))
                   (m/inc (m/- n i))))
        fri+ (fn ^double [^long i]
               (m/* (m// (m/* (m/- ns i) omega)
                         (m/* (m/inc i)
                              (m/inc (m/+ (m/- nf n) i))))
                    (m/- n i)))
        pmf (fn ^double [^long k]
              (let [[^double s ^double fk] (loop [fk 1.0
                                                  fi 1.0
                                                  s 1.0
                                                  i (m/inc mode)]
                                             (if (m/> i upper-bound)
                                               [s fk]
                                               (let [ri (double (fri i))
                                                     nfi (m/* fi ri)
                                                     sfi (m/+ s nfi)]
                                                 (if (and (m/== sfi s)
                                                          (m/> i k))
                                                   [s fk]
                                                   (recur (if (m/== k i) nfi fk)
                                                          nfi sfi (m/inc i))))))]
                (loop [fk fk
                       fi 1.0
                       s s
                       i (m/dec mode)]
                  (if (m/< i lower-bound)
                    (m// fk s)
                    (let [ri (double (fri+ i)) 
                          nfi (m// fi ri)
                          sfi (m/+ s nfi)]
                      (if (and (m/== sfi s)
                               (m/< i k))
                        (m// fk s)
                        (recur (if (m/== k i) nfi fk)
                               nfi sfi (m/dec i))))))))
        xs (range lower-bound (m/inc upper-bound))
        probs (map pmf xs)
        d (integer-discrete-distribution xs probs rng)]
    (->distribution {:pdf (fn ^double [^long x] (prot/pdf d x))
                     :cdf (fn ^double [^long x] (prot/cdf d x))
                     :icdf (fn ^long [^double p] (prot/icdf d p))
                     :sampler (fn ^long [] (prot/sample d))
                     :rng rng
                     :dimensions 1
                     :continuous? false
                     :name :fishers-noncentral-hypergeometric
                     :parameters [:ns :nf :n :omega :rng]
                     :mean (delay (prot/mean d))
                     :variance (delay (prot/variance d))
                     :lower-bound lower-bound
                     :upper-bound upper-bound})))

;; integration in log space and center a quadrature at peak + substitution t = u^D to spread peaks
(defn- wnh-lpmf
  [^long ns ^long nf ^long n ^double omega]
  (fn [^long x]
    (let [n-x (m/long-sub n x)
          l1 (m/log-combinations ns x)
          l2 (m/log-combinations nf n-x)
          
          D (m/+ (m/* omega (m/- ns x)) (m/- nf n-x))
          D- (m/dec D)

          logg (fn ^double [^double u]
                 (m/+ (m/* D- (m/log u))
                      (m/* x (m/log1p (m/- (m/pow u omega))))
                      (m/* n-x (m/log1p (m/- u)))))

          dlogg (fn ^double [^double u]
                  (m/- (m// D- u)
                       (m// (m/* x omega (m/pow u (m/dec omega)))
                            (m/- (m/expm1 (m/* omega (m/log u)) )))
                       (m// n-x (m/- 1.0 u))))

          d1 (double (dlogg zero+epsilon))
          d2 (double (dlogg one-epsilon))
          
          M (double (cond (or (and (m/not-pos? d1) (m/not-neg? d2))
                              (and (m/not-neg? d1) (m/not-pos? d2))) (logg (solver/find-root dlogg zero+epsilon one-epsilon))
                          (m/< (m/abs d1) (m/abs d2)) (logg zero+epsilon)
                          :else (logg one-epsilon)))

          integral (quad/gk-quadrature (fn ^double [^double u] (m/exp (m/- (double (logg u)) M))) 0.0 1.0
                                       {:abs 1.0e-16 :max-iters 500 :initdiv 2})]

      (m/+ l1 l2 (m/log D) M (m/log integral)))))

(defn wallenius-noncentral-hypergeometric
  "Wallenius' noncentral hypergeometric distribution.

  The number of items of the first kind in a sample of size `n` drawn one by one without replacement from an urn with `ns` items of the first and `nf` of the second kind, where the probability of taking an item is proportional to its weight, `omega` for the first kind and `1` for the second.

  Parameters (a map): `ns`, `nf`, `n` (longs), `omega` (positive double), `rng` (`nil`: a fresh generator). The support is `[max(0, n - nf), min(ns, n)]`; the probabilities are computed (by numerical integration) for the whole support at construction. Returns a univariate discrete distribution."
  [{:keys [^long ns ^long nf ^long n ^double omega rng]}]
  (let [rng (or rng (JDKRandomGenerator.))
        lower-bound (m/max 0 (m/- n nf))
        upper-bound (m/min ns n)
        lower-bound (m/min lower-bound upper-bound)
        upper-bound (m/max lower-bound upper-bound)
        lpmf (wnh-lpmf ns nf n omega)
        xs (range lower-bound (m/inc upper-bound))
        probs (map (comp m/exp lpmf) xs)
        d (integer-discrete-distribution xs probs rng)]
    (->distribution {:pdf (fn ^double [^long x] (prot/pdf d x))
                     :cdf (fn ^double [^long x] (prot/cdf d x))
                     :icdf (fn ^long [^double p] (prot/icdf d p))
                     :sampler (fn ^long [] (prot/sample d))
                     :rng rng
                     :dimensions 1
                     :continuous? false
                     :name :wallenius-noncentral-hypergeometric
                     :parameters [:ns :nf :n :omega :rng]
                     :mean (delay (prot/mean d))
                     :variance (delay (prot/variance d))
                     :lower-bound lower-bound
                     :upper-bound upper-bound})))


;;

(defn reciprocal
  "Reciprocal (log-uniform) distribution on `[a, b]`, `0 < a < b`, with the density `1 / (x log(b / a))`.

  `rng` is the random generator (`nil`: a fresh one). Returns a univariate continuous distribution."
  [^double a ^double b rng]
  (let [loga (m/log a)
        logb (m/log b)
        ldiff (m/- logb loga)
        rldiff (m// ldiff)
        mean (m/* (m/- b a) rldiff)]
    (->distribution {:pdf (fn ^double [^double x] (if (m/<= a x b) (m// rldiff x) 0.0))
                     :cdf (fn ^double [^double x] (cond
                                                   (m/< x a) 0.0
                                                   (m/> x b) 1.0
                                                   :else (m/* (m/- (m/log x) loga) rldiff)))
                     :icdf (fn ^double [^double p] (m/exp (m/+ (m/* p ldiff) loga)))
                     :rng rng
                     :dimensions 1
                     :continuous? true
                     :name :reciprocal
                     :parameters [:a :b :rng]
                     :mean mean
                     :variance (m/- (m/* 0.5 (m/- (m/sq b) (m/sq a)) rldiff) (m/sq mean))
                     :lower-bound a
                     :upper-bound b})))

(defn- log-erfcx-asymptotic
  "log(erfcx(b)) = b^2 + log(erfc(b)) for large positive b, via the standard
  asymptotic expansion of erfc, avoiding the overflow (`exp(b^2)`) / underflow
  (`erfc(b) -> 0`) that plain `b^2 + log(erfc(b))` would hit directly. Valid
  (and needed) for `b >= 6.0` or so; relative error there is already < 3e-8
  and keeps shrinking as `b` grows."
  ^double [^double b]
  (let [u (m// 1.0 (m/* 2.0 b b))
        S (m/+ 1.0 (m/* u (m/+ -1.0 (m/* u (m/+ 3.0 (m/* u (m/+ -15.0 (m/* u (m/+ 105.0 (m/* u -945.0))))))))))]
    (m/- (m/log S) (m/log b) (m/* 0.5 m/LOG_PI))))

(defn ex-gaussian
  "Ex-Gaussian distribution, the sum of an independent normal variable `(mu, sigma)` and an exponential variable with the mean `tau`.

  Parameters (a map): `mu`, `sigma`, `tau` (doubles), `rng` (`nil`: a fresh generator), `normal` and `exponential` (the constructors of these distributions, functions of a map, passed by `fastmath.random`), `parameters` (vector reported by `distribution-parameters`, default `[:mu :sigma :tau :rng]`).

  The density and the cdf are closed forms with `erfc`, evaluated in log space in the far tails; `icdf` is found by root finding. `mean` is `mu + tau`, `variance` is `sigma^2 + tau^2`. Returns a univariate continuous distribution."
  [{:keys [^double mu ^double sigma ^double tau rng normal exponential parameters]
    :or {parameters [:mu :sigma :tau :rng]}}]
  (let [rng (or rng (JDKRandomGenerator.))
        N (normal {:mu mu :sd sigma :rng rng})
        E (exponential {:mean tau :rng rng})
        sigma2 (m/* sigma sigma)
        sigma2tau (m// sigma2 tau)
        s1 (m// (m/+ mu (m/* 0.5 sigma2tau)) tau)
        denom (m/* m/SQRT2 sigma)
        s2 (m// (m/+ mu sigma2tau) denom)
        log-half (m/log 0.5)
        ;; f(x) = 0.5 * exp(a) * erfc(b), a = s1 - x/tau, b = s2 - x/denom;
        ;; computed via log-space when b is large so that erfc(b)'s underflow
        ;; towards 0 doesn't collide with exp(a)'s overflow towards Infinity
        ;; and produce Infinity * 0 -> NaN (hit for x deep in the left tail,
        ;; or any x when mu is large relative to sigma/tau).
        f (fn ^double [^double x]
            (let [a (m/- s1 (m// x tau))
                  b (m/- s2 (m// x denom))
                  log-f (if (m/>= b 6.0)
                          (m/+ log-half a (m/- (m/* b b)) (log-erfcx-asymptotic b))
                          (m/+ log-half a (m/log (special/erfc b))))]
              (m/exp log-f)))
        cdf (fn ^double [^double x]
              (cond
                (m/neg-inf? x) 0.0
                (m/pos-inf? x) 1.0
                :else (m/constrain (m/- (double (prot/cdf N x)) (double (f x))) 0.0 1.0)))]
    (->distribution {:pdf (fn ^double [^double x]
                            (if (or (m/neg-inf? x) (m/pos-inf? x))
                              0.0
                              (m// (double (f x)) tau)))
                     :cdf cdf
                     :icdf (fn ^double [^double p] (icdf-solver cdf p mu sigma))
                     :sampler (fn ^double [] (m/+ (double (prot/sample N))
                                                 (double (prot/sample E))))
                     :rng rng
                     :mean (m/+ mu tau)
                     :variance (m/+ sigma2 (m/* tau tau))
                     :dimensions 1
                     :continuous? true
                     :name :ex-gaussian
                     :parameters parameters
                     :lower-bound ##-Inf
                     :upper-bound ##Inf})))

(declare beta-binomial-core)

(defn beta-binomial
  "Beta-binomial distribution `(alpha, beta, n)`. The optional `parameters` (default `[:alpha :beta :n :rng]`) is what
  `distribution-parameters` reports, for the constructors accepting other parameters (`bb`)."
  ([^double alpha ^double beta ^long n rng] (beta-binomial alpha beta n rng [:alpha :beta :n :rng]))
  ([alpha beta n rng parameters] (beta-binomial-core alpha beta n rng parameters)))

(defn- beta-binomial-core
  [alpha beta n rng parameters]
  (let [alpha (double alpha)
        beta (double beta)
        n (long n)
        rng (or rng (JDKRandomGenerator.))
        a+b (m/+ alpha beta)
        n+ (m/inc n)
        lpmf-const (m/- (m/+ (special/log-gamma n+)
                             (special/log-gamma a+b))
                        (special/log-gamma (m/+ n a+b))
                        (special/log-gamma alpha)
                        (special/log-gamma beta))
        lpmf (fn ^double [^double x]
               (m/- (m/+ lpmf-const
                         (special/log-gamma (m/+ alpha x))
                         (special/log-gamma (m/- (m/+ n beta) x)))
                    (special/log-gamma (m/inc x))
                    (special/log-gamma (m/- n+ x))))
        xs (range (m/inc n))
        d (integer-discrete-distribution xs (v/exp (map lpmf xs)) rng)
        mean (m// (m/* n alpha) a+b)]
    (->distribution {:pdf (fn ^double [^long x] (prot/pdf d x))
                     :cdf (fn ^double [^long x] (prot/cdf d x))
                     :icdf (fn ^long [^double p] (prot/icdf d p))
                     :sampler (fn ^long [] (prot/sample d))
                     :rng rng
                     :mean mean
                     :variance (m/* mean (m// (m/* beta (m/+ a+b n))
                                              (m/* a+b (m/inc a+b))))
                     :dimensions 1
                     :continuous? false
                     :name :beta-binomial
                     :parameters parameters
                     :lower-bound 0
                     :upper-bound n})))

(defn zero-inflated-beta-binomial
  "Zero-inflated beta-binomial distribution (gamlss `ZIBB`): the probability `nu` of an extra zero, otherwise a beta-binomial with the mean `mu`, the dispersion `sigma` and `bd` trials.

  Parameters (a map): `mu`, `sigma`, `nu` (doubles), `bd` (long), `rng` (`nil`: a fresh generator). Returns a univariate discrete distribution on `[0, bd]`."
  [{:keys [^double mu ^double sigma ^long bd ^double nu rng]}]
  (let [rng (or rng (JDKRandomGenerator.))
        alpha (m// mu sigma)
        beta (m// (m/- 1.0 mu) sigma)
        dist (beta-binomial alpha beta bd rng)
        p0 (double (prot/pdf dist 0))
        nu- (m/- 1.0 nu)
        pdf0 (m/+ nu (m/* nu- p0))
        mean (m/* nu- bd mu)]
    (->distribution {:pdf (fn ^double [^long x] (if (m/zero? x)
                                                 pdf0
                                                 (m/* nu- (double (prot/pdf dist x)))))
                     :cdf (fn ^double [^long x] (if (m/neg? x)
                                                 0.0
                                                 (m/+ nu (m/* nu- (double (prot/cdf dist x))))))
                     :icdf (fn ^long [^double p]
                             (if (m/<= p pdf0)
                               0.0
                               (prot/icdf dist (m// (m/- p nu) nu-))))
                     :sampler (fn ^long []
                                (let [v (double (prot/drandom rng))]
                                  (if (m/< v nu)
                                    0.0
                                    (prot/sample dist))))
                     :rng rng
                     :mean mean
                     :variance (delay (m/+ (m/* nu- (double (prot/variance dist)))
                                           (m/* mean nu bd mu)))
                     :dimensions 1
                     :continuous? false
                     :name :zero-inflated-beta-binomial
                     :parameters [:mu :sigma :bd :nu :rng]
                     :lower-bound 0
                     :upper-bound bd})))

(defn zero-adjusted-beta-binomial
  "Zero-adjusted (hurdle) beta-binomial distribution (gamlss `ZABB`): zero with the probability `nu`, otherwise a zero-truncated beta-binomial with the mean `mu`, the dispersion `sigma` and `bd` trials.

  Parameters (a map): `mu`, `sigma`, `nu` (doubles), `bd` (long), `rng` (`nil`: a fresh generator). Returns a univariate discrete distribution on `[0, bd]`."
  [{:keys [^double mu ^double sigma ^long bd ^double nu rng]}]
  (let [rng (or rng (JDKRandomGenerator.))
        alpha (m// mu sigma)
        beta (m// (m/- 1.0 mu) sigma)
        dist (beta-binomial alpha beta bd rng)
        p0 (double (prot/pdf dist 0))
        p0- (m/- 1.0 p0)
        nu- (m/- 1.0 nu)
        mean (m// (m/* nu- bd mu) p0-)]
    (->distribution {:pdf (fn ^double [^long x] (if (m/zero? x)
                                                 nu
                                                 (m// (m/* nu- (double (prot/pdf dist x))) p0-)))
                     :cdf (fn ^double [^long x] (cond
                                                 (m/neg? x) 0.0
                                                 (m/zero? x) nu
                                                 :else (m/+ nu (m// (m/* nu- (m/- (double (prot/cdf dist x)) p0 )) p0-))))
                     :icdf (fn ^long [^double p]
                             (cond
                               (m/<= p nu) 0
                               (m/>= p 1.0) bd
                               :else (let [np (m/+ p0 (m// (m/* p0- (m/- p nu)) nu-))]
                                       (prot/icdf dist np))))
                     :rng rng
                     :mean mean
                     :variance (delay (m/- (m// (m/* nu- (m/+ (m/* bd mu (m/- 1.0 mu) (m/inc (m// (m/* sigma (m/dec bd)) (m/inc sigma))))
                                                              (m/* bd bd mu mu))) p0-)
                                           (m/* mean mean)))
                     :dimensions 1
                     :continuous? false
                     :name :zero-adjusted-beta-binomial
                     :parameters [:mu :sigma :bd :nu :rng]
                     :lower-bound 0
                     :upper-bound bd})))

(defn zero-inflated-binomial
  "Zero-inflated binomial distribution (gamlss `ZIBI`): an extra zero with the probability `sigma`, otherwise a binomial with `bd` trials and the success probability `mu`.

  Parameters (a map): `mu`, `sigma` (doubles), `bd` (long), `rng` (`nil`: a fresh generator), `binomial` (the constructor of the binomial distribution, a function of a map). Returns a univariate discrete distribution on `[0, bd]`."
  [{:keys [^double mu ^double sigma ^long bd rng binomial]}]
  (let [rng (or rng (JDKRandomGenerator.))
        dist (binomial {:trials bd :p mu :rng rng})
        p0 (double (prot/pdf dist 0))
        sigma- (m/- 1.0 sigma)
        pdf0 (m/+ sigma (m/* sigma- p0))
        mean (m/* sigma- bd mu)]
    (->distribution {:pdf (fn ^double [^long x] (if (m/zero? x)
                                                 pdf0
                                                 (m/* sigma- (double (prot/pdf dist x)))))
                     :cdf (fn ^double [^long x] (if (m/neg? x)
                                                 0.0
                                                 (m/+ sigma (m/* sigma- (double (prot/cdf dist x))))))
                     :icdf (fn ^long [^double p]
                             (if (m/<= p pdf0)
                               0.0
                               (prot/icdf dist (m// (m/- p sigma) sigma-))))
                     :sampler (fn ^long []
                                (let [v (double (prot/drandom rng))]
                                  (if (m/< v sigma)
                                    0.0
                                    (prot/sample dist))))
                     :rng rng
                     :mean mean
                     :variance (delay (m/+ (m/* sigma- (double (prot/variance dist)))
                                           (m/* mean sigma bd mu)))
                     :dimensions 1
                     :continuous? false
                     :name :zero-inflated-binomial
                     :parameters [:mu :sigma :bd :rng]
                     :lower-bound 0
                     :upper-bound bd})))

(defn zero-adjusted-binomial
  "Zero-adjusted (hurdle) binomial distribution (gamlss `ZABI`): zero with the probability `sigma`, otherwise a zero-truncated binomial with `bd` trials and the success probability `mu`.

  Parameters (a map): `mu`, `sigma` (doubles), `bd` (long), `rng` (`nil`: a fresh generator), `binomial` (the constructor of the binomial distribution, a function of a map). Returns a univariate discrete distribution on `[0, bd]`."
  [{:keys [^double mu ^double sigma ^long bd rng binomial]}]
  (let [rng (or rng (JDKRandomGenerator.))
        dist (binomial {:trials bd :p mu :rng rng})
        p0 (double (prot/pdf dist 0))
        p0- (m/- 1.0 p0)
        sigma- (m/- 1.0 sigma)
        mean (m// (m/* sigma- bd mu) p0-)]
    (->distribution {:pdf (fn ^double [^long x] (if (m/zero? x)
                                                 sigma
                                                 (m// (m/* sigma- (double (prot/pdf dist x))) p0-)))
                     :cdf (fn ^double [^long x] (cond
                                                 (m/neg? x) 0.0
                                                 (m/zero? x) sigma
                                                 :else (m/+ sigma (m// (m/* sigma- (m/- (double (prot/cdf dist x)) p0)) p0-))))
                     :icdf (fn ^long [^double p]
                             (cond
                               (m/<= p sigma) 0
                               (m/>= p 1.0) bd
                               :else (let [np (m/+ p0 (m// (m/* p0- (m/- p sigma)) sigma-))]
                                       (prot/icdf dist np))))
                     :rng rng
                     :mean mean
                     :variance (m/* mean (m/- (m/inc (m/* bd mu)) mu mean))
                     :dimensions 1
                     :continuous? false
                     :name :zero-adjusted-binomial
                     :parameters [:mu :sigma :bd :rng]
                     :lower-bound 0
                     :upper-bound bd})))

(defn zero-inflated-negative-binomial
  "Zero-inflated negative binomial distribution (gamlss `ZINBI`): an extra zero with the probability `nu`, otherwise a type I negative binomial (`nbi`) with the mean `mu` and the dispersion `sigma`.

  Parameters (a map): `mu`, `sigma`, `nu` (doubles), `rng` (`nil`: a fresh generator), `nbi` (the constructor of the `nbi` distribution, a function of a map). Returns a univariate discrete distribution on the non-negative integers."
  [{:keys [^double mu ^double sigma ^double nu rng nbi]}]
  (let [rng (or rng (JDKRandomGenerator.))
        dist (nbi {:mu mu :sigma sigma :rng rng})
        p0 (double (prot/pdf dist 0))
        nu- (m/- 1.0 nu)
        pdf0 (m/+ nu (m/* nu- p0))
        mean (m/* nu- mu)]
    (->distribution {:pdf (fn ^double [^long x] (if (m/zero? x)
                                                 pdf0
                                                 (m/* nu- (double (prot/pdf dist x)))))
                     :cdf (fn ^double [^long x] (if (m/neg? x)
                                                 0.0
                                                 (m/+ nu (m/* nu- (double (prot/cdf dist x))))))
                     :icdf (fn ^long [^double p]
                             (if (m/<= p pdf0)
                               0.0
                               (prot/icdf dist (m// (m/- p nu) nu-))))
                     :sampler (fn ^long []
                                (let [v (double (prot/drandom rng))]
                                  (if (m/< v nu)
                                    0.0
                                    (prot/sample dist))))
                     :rng rng
                     :mean mean
                     :variance (delay (m/+ (m/* nu- (double (prot/variance dist)))
                                           (m/* mean nu mu)))
                     :dimensions 1
                     :continuous? false
                     :name :zero-inflated-negative-binomial
                     :parameters [:mu :sigma :nu :rng]
                     :lower-bound 0
                     :upper-bound Integer/MAX_VALUE})))

(defn zero-adjusted-negative-binomial
  "Zero-adjusted (hurdle) negative binomial distribution (gamlss `ZANBI`): zero with the probability `nu`, otherwise a zero-truncated type I negative binomial (`nbi`) with the mean `mu` and the dispersion `sigma`.

  Parameters (a map): `mu`, `sigma`, `nu` (doubles), `rng` (`nil`: a fresh generator), `nbi` (the constructor of the `nbi` distribution, a function of a map). Returns a univariate discrete distribution on the non-negative integers."
  [{:keys [^double mu ^double sigma ^double nu rng nbi]}]
  (let [rng (or rng (JDKRandomGenerator.))
        dist (nbi {:mu mu :sigma sigma :rng rng})
        p0 (double (prot/pdf dist 0))
        p0- (m/- 1.0 p0)
        nu- (m/- 1.0 nu)
        mean (m// (m/* nu- mu) p0-)]
    (->distribution {:pdf (fn ^double [^long x] (if (m/zero? x)
                                                 nu
                                                 (m// (m/* nu- (double (prot/pdf dist x))) p0-)))
                     :cdf (fn ^double [^long x] (cond
                                                 (m/neg? x) 0.0
                                                 (m/zero? x) nu
                                                 :else (m/+ nu (m// (m/* nu- (m/- (double (prot/cdf dist x)) p0)) p0-))))
                     :icdf (fn ^long [^double p]
                             (if (m/<= p nu)
                               0
                               (let [np (m/+ p0 (m// (m/* p0- (m/- p nu)) nu-))]
                                 (prot/icdf dist np))))
                     :rng rng
                     :mean mean
                     :variance (delay (m/- (m// (m/* nu- (m/+ (double (prot/variance dist)) (m/* mu mu))) p0-) (m/* mean mean)))
                     :dimensions 1
                     :continuous? false
                     :name :zero-adjusted-negative-binomial
                     :parameters [:mu :sigma :nu :rng]
                     :lower-bound 0
                     :upper-bound Integer/MAX_VALUE})))

(defn zero-inflated-poisson
  "Zero-inflated Poisson distribution (gamlss `ZIP`): an extra zero with the probability `sigma`, otherwise a Poisson distribution with the rate `mu`. The mean is `(1 - sigma) mu`.

  Parameters (a map): `mu`, `sigma` (doubles), `rng` (`nil`: a fresh generator), `poisson` (the constructor of the Poisson distribution, a function of a map), `parameters` (vector reported by `distribution-parameters`, default `[:mu :sigma :rng]`). `icdf` is the smallest integer with `cdf >= p`. Returns a univariate discrete distribution on the non-negative integers."
  [{:keys [^double mu ^double sigma rng poisson parameters]
    :or {parameters [:mu :sigma :rng]}}]
  (let [rng (or rng (JDKRandomGenerator.))
        dist (poisson {:p mu :rng rng})
        p0 (double (prot/pdf dist 0))
        sigma- (m/- 1.0 sigma)
        pdf0 (m/+ sigma (m/* sigma- p0))
        mean (m/* sigma- mu)
        cdf (fn [^long x] (if (m/neg? x)
                            0.0
                            (m/+ sigma (m/* sigma- (double (prot/cdf dist x))))))]
    (->distribution {:pdf (fn ^double [^long x] (if (m/zero? x)
                                                 pdf0
                                                 (m/* sigma- (double (prot/pdf dist x)))))
                     :cdf cdf
                     ;; the quantile of the Poisson part is only a start, the final value is the smallest x with cdf(x) >= p
                     :icdf (fn ^long [^double p]
                             (cond
                               (m/<= p (double (cdf 0))) 0
                               (m/>= p 1.0) (prot/icdf dist (m// (m/- p sigma) sigma-))
                               :else (smallest-with-cdf cdf p 0 (long (prot/icdf dist (m// (m/- p sigma) sigma-))))))
                     :sampler (fn ^long []
                                (let [v (double (prot/drandom rng))]
                                  (if (m/< v sigma)
                                    0.0
                                    (prot/sample dist))))
                     :rng rng
                     :mean mean
                     :variance (m/* mean (m/inc (m/* mu sigma)))
                     :dimensions 1
                     :continuous? false
                     :name :zero-inflated-poisson
                     :parameters parameters
                     :lower-bound 0
                     :upper-bound Integer/MAX_VALUE})))

(defn zero-adjusted-poisson
  "Zero-adjusted (hurdle) Poisson distribution (gamlss `ZAP`): zero with the probability `sigma`, otherwise a zero-truncated Poisson with the rate `mu`.

  Parameters (a map): `mu`, `sigma` (doubles), `rng` (`nil`: a fresh generator), `poisson` (the constructor of the Poisson distribution, a function of a map). Returns a univariate discrete distribution on the non-negative integers."
  [{:keys [^double mu ^double sigma rng poisson]}]
  (let [rng (or rng (JDKRandomGenerator.))
        dist (poisson {:p mu :rng rng})
        p0 (double (prot/pdf dist 0))
        p0- (m/- 1.0 p0)
        sigma- (m/- 1.0 sigma)
        mean (m// (m/* sigma- mu) p0-)]
    (->distribution {:pdf (fn ^double [^long x] (if (m/zero? x)
                                                 sigma
                                                 (m// (m/* sigma- (double (prot/pdf dist x))) p0-)))
                     :cdf (fn ^double [^long x] (cond
                                                 (m/neg? x) 0.0
                                                 (m/zero? x) sigma
                                                 :else (m/+ sigma (m// (m/* sigma- (m/- (double (prot/cdf dist x)) p0)) p0-))))
                     :icdf (fn ^long [^double p]
                             (if (m/<= p sigma)
                               0
                               (let [np (m/+ p0 (m// (m/* p0- (m/- p sigma)) sigma-))]
                                 (prot/icdf dist np))))
                     :rng rng
                     :mean mean
                     :variance (delay (m/- (m// (m/* sigma- (m/+ (double (prot/variance dist)) (m/* mu mu))) p0-) (m/* mean mean)))
                     :dimensions 1
                     :continuous? false
                     :name :zero-adjusted-poisson
                     :parameters [:mu :sigma :rng]
                     :lower-bound 0
                     :upper-bound Integer/MAX_VALUE})))

;;

(defn zero-adjusted-gamma
  "Zero-adjusted gamma distribution (gamlss `ZAGA`): the point mass `nu` at zero, otherwise a gamma distribution with the mean `mu` and the coefficient of variation `sigma` (shape `1 / sigma^2`, scale `sigma^2 mu`).

  Parameters (a map): `mu`, `sigma`, `nu` (doubles), `rng` (`nil`: a fresh generator), `gamma` (the constructor of the gamma distribution, a function of a map). Returns a univariate continuous (mixed) distribution on `[0, Inf)`."
  [{:keys [^double mu ^double sigma ^double nu rng gamma]}]
  (let [rng (or rng (JDKRandomGenerator.))
        sigma2 (m/sq sigma)
        shape (m// sigma2)
        scale (m/* sigma2 mu)
        dist (gamma {:shape shape :scale scale :rng rng})
        nu- (m/- 1.0 nu)
        mean (m/* nu- mu)]
    (->distribution {:pdf (fn ^double [^double x] (if (m/zero? x)
                                                   nu
                                                   (m/* nu- (double (prot/pdf dist x)))))
                     :cdf (fn ^double [^double x] (cond
                                                   (m/neg? x) 0.0
                                                   (m/zero? x) nu
                                                   :else (m/+ nu (m/* nu- (double (prot/cdf dist x))))))
                     :icdf (fn ^double [^double p]
                             (if (m/<= p nu)
                               0.0
                               (let [np (m// (m/- p nu) nu-)]
                                 (double (prot/icdf dist np)))))
                     :sampler (fn ^double []
                                (let [u (double (prot/drandom rng))]
                                  (if (m/< u nu)
                                    0.0
                                    (prot/sample dist))))
                     :rng rng
                     :mean mean
                     :variance (m/* mean mu (m/+ sigma2 nu))
                     :dimensions 1
                     :continuous? true
                     :name :zero-adjusted-gamma
                     :parameters [:mu :sigma :nu :rng]
                     :lower-bound 0
                     :upper-bound ##Inf})))

(defn zero-adjusted-inverse-gaussian
  "Zero-adjusted inverse Gaussian distribution (gamlss `ZAIG`): the point mass `nu` at zero, otherwise an inverse Gaussian distribution with the mean `mu` and `lambda = 1 / sigma^2`.

  Parameters (a map): `mu`, `sigma`, `nu` (doubles), `rng` (`nil`: a fresh generator), `inverse-gaussian` (the constructor of the inverse Gaussian distribution, a function of a map). Returns a univariate continuous (mixed) distribution on `[0, Inf)`."
  [{:keys [^double mu ^double sigma ^double nu rng inverse-gaussian]}]
  (let [rng (or rng (JDKRandomGenerator.))
        sigma2 (m/sq sigma)
        lambda (m// sigma2)
        dist (inverse-gaussian {:mu mu :lambda lambda :rng rng})
        nu- (m/- 1.0 nu)
        mean (m/* nu- mu)]
    (->distribution {:pdf (fn ^double [^double x] (if (m/zero? x)
                                                   nu
                                                   (m/* nu- (double (prot/pdf dist x)))))
                     :cdf (fn ^double [^double x] (cond
                                                   (m/neg? x) 0.0
                                                   (m/zero? x) nu
                                                   :else (m/+ nu (m/* nu- (double (prot/cdf dist x))))))
                     :icdf (fn ^double [^double p]
                             (if (m/<= p nu)
                               0.0
                               (let [np (m// (m/- p nu) nu-)]
                                 (double (prot/icdf dist np)))))
                     :sampler (fn ^double []
                                (let [u (double (prot/drandom rng))]
                                  (if (m/< u nu)
                                    0.0
                                    (prot/sample dist))))
                     :rng rng
                     :mean mean
                     :variance (m/* mean mu (m/+ (m/* sigma2 mu) nu))
                     :dimensions 1
                     :continuous? true
                     :name :zero-adjusted-inverse-gaussian
                     :parameters [:mu :sigma :nu :rng]
                     :lower-bound 0
                     :upper-bound ##Inf})))

(defn generalized-extreme-value
  "Generalized extreme value distribution with the location `mu`, the scale `sigma` and the shape `xi`, `cdf(x) = exp(-(1 + xi (x - mu) / sigma)^(-1 / xi))`; `xi = 0` is the Gumbel distribution.

  `rng` is the random generator (`nil`: a fresh one). `mean` is `##Inf` for `xi >= 1` and `variance` is `##Inf` for `xi >= 0.5`. The support is bounded on the left for `xi > 0` and on the right for `xi < 0`. Returns a univariate continuous distribution."
  [^double mu ^double sigma ^double xi rng]
  (let [boundary (m/- mu (m// sigma xi))
        pdf (if (m/zero? xi)
              (fn ^double [^double x]
                (let [s (m// (m/- x mu) sigma)
                      es (m/exp (m/- s))]
                  (if (m/pos-inf? es)
                    0.0
                    (m// (m/* es (m/exp (m/- es))) sigma))))
              (let [-rxi (m// -1.0 xi)
                    -rxi- (m/dec -rxi)                    ]
                (fn ^double [^double x]
                  (let [s (m// (m/- x mu) sigma)
                        xis (m/* xi s)]
                    (if (m/<= xis -1.0)
                      0.0
                      (let [t (m/inc xis)]
                        (if (m/pos-inf? t)
                          0.0
                          (m// (m/* (m/pow t -rxi-) (m/exp (m/- (m/pow t -rxi)))) sigma))))))))
        cdf (if (m/zero? xi)
              (fn ^double [^double x]
                (let [-s (m// (m/- mu x) sigma)
                      es (m/exp -s)]
                  (m/exp (m/- es))))
              (let [-rxi (m// -1.0 xi)]
                (fn ^double [^double x]
                  (let [s (m// (m/- x mu) sigma)
                        xis (m/* xi s)]
                    (cond
                      (m/> xis -1.0) (m/exp (m/- (m/pow (m/inc xis) -rxi)))
                      (and (m/pos? xi) (m/<= s (m// -1.0 xi))) 0.0
                      :else 1.0)))))
        icdf (if (m/zero? xi)
               (fn ^double [^double p]
                 (m/- mu (m/* sigma (m/log (m/- (m/log p))))))
               (fn ^double [^double p]
                 (m/+ mu (m/* (m// sigma xi) (m/dec (m/pow (m/- (m/log p)) (m/- xi)))))))]
    (->distribution {:pdf pdf
                     :cdf cdf
                     :icdf icdf
                     :rng rng
                     :mean (delay (cond
                                    (m/zero? xi) (m/+ mu (m/* sigma m/GAMMA))
                                    (m/< xi 1.0) (m/+ mu (m/* (m// sigma xi) (m/dec (special/gamma (m/- 1.0 xi)))))
                                    :else ##Inf))
                     :variance (delay (cond
                                        (m/zero? xi) (m/* m/PI2 m/SIXTH sigma sigma)
                                        (m/< xi 0.5) (m/* sigma sigma (m// (m/- (special/gamma (m/- 1.0 xi xi))
                                                                                (m/sq (special/gamma (m/- 1.0 xi))))
                                                                           (m/* xi xi)))
                                        :else ##Inf))
                     :dimensions 1
                     :continuous? true
                     :name :generalized-extreme-value
                     :parameters [:mu :sigma :xi :rng]
                     :lower-bound (if (m/not-pos? xi) ##-Inf boundary)
                     :upper-bound (if (m/not-neg? xi) ##Inf  boundary)})))

(defn generalized-logistic
  "Generalized logistic distribution (type I) with the location `mu`, the scale `sigma` and the shape `alpha`, `cdf(x) = (1 + exp(-(x - mu) / sigma))^(-alpha)`.

  `rng` is the random generator (`nil`: a fresh one). Returns a univariate continuous distribution on the reals."
  [^double mu ^double sigma ^double alpha rng]
  (let [diff (m/- (m/log alpha) (m/log sigma))
        -ralpha (m// -1.0 alpha)
        -alpha (m/- alpha)]
    (->distribution {:lpdf (fn ^double [^double x] (let [z (m// (m/- x mu) sigma)]
                                                    (m/- diff z (m/* (m/inc alpha) (m/logaddexp 0.0 (m/- z)))))) 
                     :cdf (fn ^double [^double x] (let [z (m// (m/- x mu) sigma)]
                                                   (m/pow (m/inc (m/exp (m/- z))) -alpha)))
                     :icdf (fn ^double [^double p] (m/- mu (m/* sigma (m/log (m/dec (m/pow p -ralpha))))))
                     :rng rng
                     :mean (delay (m/+ mu (m/* sigma (m/+ (special/digamma alpha) m/GAMMA))))
                     :variance (delay (m/* sigma sigma (m/+ (special/trigamma alpha) (m/* m/PI2 m/SIXTH))))
                     :dimensions 1
                     :continuous? true
                     :name :generalized-logistic
                     :parameters [:mu :sigma :alpha :rng]
                     :lower-bound ##-Inf
                     :upper-bound ##Inf})))

(defn generalized-pareto
  "Generalized Pareto distribution with the location `mu`, the scale `sigma` and the shape `xi`; `xi = 0` is the exponential distribution shifted by `mu`.

  `rng` is the random generator (`nil`: a fresh one). `mean` is `##Inf` for `xi >= 1` and `variance` is `##Inf` for `xi >= 0.5`. The support is `[mu, Inf)`, and `[mu, mu - sigma / xi]` for `xi < 0`. Returns a univariate continuous distribution."
  [^double mu ^double sigma ^double xi rng]
  (let [boundary (m/- mu (m// sigma xi))
        pdf (if (m/zero? xi)
              (fn ^double [^double x]
                (if (m/< x mu)
                  0.0
                  (m// (m/exp (m/- (m// (m/- x mu) sigma))) sigma)))
              (let [-rxi (m// -1.0 xi)
                    -rxi- (m/dec -rxi)]
                (fn ^double [^double x]
                  (if (m/< x mu)
                    0.0
                    (let [z (m// (m/- x mu) sigma)
                          t (m/inc (m/* xi z))]
                      (if (m/<= t 0.0)
                        0.0
                        (m// (m/pow t -rxi-) sigma)))))))
        cdf (if (m/zero? xi)
              (fn ^double [^double x]
                (if (m/< x mu)
                  0.0
                  (m/- 1.0 (m/exp (m/- (m// (m/- x mu) sigma))))))
              (let [-rxi (m// -1.0 xi)]
                (fn ^double [^double x]
                  (if (m/< x mu)
                    0.0
                    (let [z (m// (m/- x mu) sigma)
                          t (m/inc (m/* xi z))]
                      (if (m/<= t 0.0)
                        1.0
                        (m/- 1.0 (m/pow t -rxi))))))))
        icdf (if (m/zero? xi)
               (fn ^double [^double p]
                 (m/- mu (m/* sigma (m/log (m/- 1.0 p)))))
               (fn ^double [^double p]
                 (m/+ mu (m/* (m// sigma xi) (m/dec (m/pow (m/- 1.0 p) (m/- xi)))))))]
    (->distribution {:pdf pdf
                     :cdf cdf
                     :icdf icdf
                     :rng rng
                     :mean (delay (cond
                                    (m/zero? xi) (m/+ mu sigma)
                                    (m/< xi 1.0) (m/+ mu (m// sigma (m/- 1.0 xi)))
                                    :else ##Inf))
                     :variance (delay (cond
                                        (m/zero? xi) (m/sq sigma)
                                        (m/< xi 0.5) (m// (m/sq sigma) (m/* (m/sq (m/- 1.0 xi)) (m/- 1.0 (m/* 2.0 xi))))
                                        :else ##Inf))
                     :dimensions 1
                     :continuous? true
                     :name :generalized-pareto
                     :parameters [:mu :sigma :xi :rng]
                     :lower-bound mu
                     :upper-bound (if (m/neg? xi) boundary ##Inf)})))

(defn generalized-exponential
  "Generalized (exponentiated) exponential distribution with the shape `alpha` and the rate `lambda`, `cdf(x) = (1 - exp(-lambda x))^alpha` for `x >= 0`.

  `rng` is the random generator (`nil`: a fresh one). Returns a univariate continuous distribution on `[0, Inf)`."
  [^double alpha ^double lambda rng]
  (let [alpha-1 (m/dec alpha)
        log-alpha-lambda (m/+ (m/log alpha) (m/log lambda))]
    (->distribution {:lpdf (fn ^double [^double x]
                             (cond
                               (m/neg? x) ##-Inf
                               (m/zero? alpha-1) (m/+ log-alpha-lambda (m/- (m/* lambda x)))
                               :else (let [z (m/* lambda x)
                                           em1z (m/- 1.0 (m/exp (m/- z)))]
                                       (if (m/zero? em1z)
                                         (if (m/pos? alpha-1) ##-Inf ##Inf)
                                         (m/+ log-alpha-lambda (m/- z) (m/* alpha-1 (m/log em1z)))))))
                     :cdf (fn ^double [^double x]
                            (if (m/neg? x)
                              0.0
                              (m/pow (m/- 1.0 (m/exp (m/- (m/* lambda x)))) alpha)))
                     :icdf (fn ^double [^double p]
                             (m// (m/- (m/log (m/- 1.0 (m/pow p (m// 1.0 alpha))))) lambda))
                     :rng rng
                     :mean (delay (m// (m/+ (special/digamma (m/inc alpha)) m/GAMMA) lambda))
                     :variance (delay (m// (m/- (m/* m/PI2 m/SIXTH) (special/trigamma (m/inc alpha))) (m/* lambda lambda)))
                     :dimensions 1
                     :continuous? true
                     :name :generalized-exponential
                     :parameters [:alpha :lambda :rng]
                     :lower-bound 0
                     :upper-bound ##Inf})))

(defn generalized-gamma
  "Generalized gamma distribution in the gamlss `GG` parameterisation: `u = theta (y / mu)^nu` is gamma distributed with the shape `theta = 1 / (sigma^2 nu^2)`.

  Parameters (a map): `mu` (scale), `sigma`, `nu` (doubles; `nu = 0` is the log-normal distribution with `scale = log(mu)` and `shape = sigma`), `rng` (`nil`: a fresh generator), `gamma` and `log-normal` (the constructors of these distributions, functions of a map). `mean` and `variance` are `##Inf` when they diverge. Returns a univariate continuous distribution on `[0, Inf)`."
  [{:keys [^double mu ^double sigma ^double nu rng gamma log-normal]}]
  (let [rng (or rng (JDKRandomGenerator.))]
    (if (m/zero? nu)
      (log-normal {:scale (m/log mu) :shape sigma :rng rng})
      (let [theta (m// 1.0 (m/* sigma sigma nu nu))
            log-theta (m/log theta)
            log-mu (m/log mu)
            log-abs-nu (m/log (m/abs nu))
            log-gamma-theta (special/log-gamma theta)
            pos-nu? (m/pos? nu)
            r-nu (m// 1.0 nu)
            dist (gamma {:shape theta :scale 1.0 :rng rng})
            mean-fn (fn ^double []
                      (m/* mu (m/pow theta (m/- r-nu))
                           (m/exp (m/- (special/log-gamma (m/+ theta r-nu)) log-gamma-theta))))]
        (->distribution {:lpdf (fn ^double [^double y]
                                 (if (m/pos? y)
                                   (let [log-y (m/log y)
                                         log-u (m/+ log-theta (m/* nu (m/- log-y log-mu)))]
                                     (if (m/pos-inf? log-u)
                                       ##-Inf
                                       (m/- (m/+ log-abs-nu (m/* theta log-u)) (m/exp log-u) log-gamma-theta log-y)))
                                   ##-Inf))
                         :cdf (fn ^double [^double y]
                                (if (m/not-pos? y)
                                  0.0
                                  (let [u (m/* theta (m/pow (m// y mu) nu))]
                                    (if (m/pos-inf? u)
                                      (if pos-nu? 1.0 0.0)
                                      (let [c (double (prot/cdf dist u))]
                                        (if pos-nu? c (m/- 1.0 c)))))))
                         :icdf (fn ^double [^double p]
                                 (let [u (double (prot/icdf dist (if pos-nu? p (m/- 1.0 p))))]
                                   (m/* mu (m/pow (m// u theta) r-nu))))
                         :rng rng
                         :mean (delay (if (m/> nu (m/- (m// 1.0 (m/* sigma sigma))))
                                        (double (mean-fn))
                                        ##Inf))
                         :variance (delay (if (m/> nu (m/- (m// 1.0 (m/* 2.0 sigma sigma))))
                                            (let [m1 (double (mean-fn))
                                                  m2 (m/* mu mu (m/pow theta (m/* -2.0 r-nu))
                                                          (m/exp (m/- (special/log-gamma (m/+ theta (m/* 2.0 r-nu)))
                                                                      (special/log-gamma theta))))]
                                              (m/- m2 (m/* m1 m1)))
                                            ##Inf))
                         :dimensions 1
                         :continuous? true
                         :name :generalized-gamma
                         :parameters [:mu :sigma :nu :rng]
                         :lower-bound 0
                         :upper-bound ##Inf})))))

(defn generalized-normal
  "Generalized normal (exponential power, Subbotin) distribution with the location `mu`, the scale `alpha` and the shape `beta`; `beta = 2` is a normal distribution and `beta = 1` the Laplace distribution.

  Parameters (a map): `mu`, `alpha`, `beta` (doubles), `rng` (`nil`: a fresh generator), `gamma` (the constructor of the gamma distribution, a function of a map). `mean` is `mu` and `variance` is `alpha^2 Gamma(3 / beta) / Gamma(1 / beta)`. Returns a univariate continuous distribution on the reals."
  [{:keys [^double mu ^double alpha ^double beta rng gamma]}]
  (let [rng (or rng (JDKRandomGenerator.))
        inv-beta (m// 1.0 beta)
        log-gamma-inv-beta (special/log-gamma inv-beta)
        log-const (m/- (m/log beta) (m/log 2.0) (m/log alpha) log-gamma-inv-beta)
        variance (m/* alpha alpha (m/exp (m/- (special/log-gamma (m/* 3.0 inv-beta)) log-gamma-inv-beta)))
        dist (gamma {:shape inv-beta :scale 1.0 :rng rng})]
    (->distribution
     {:lpdf (fn ^double [^double x]
              (let [z (m/pow (m// (m/abs (m/- x mu)) alpha) beta)]
                (m/- log-const z)))
      :cdf (fn ^double [^double x]
             (let [d (m/- x mu)
                   z (m/pow (m// (m/abs d) alpha) beta)
                   s (m/signum d)]
               (if (m/pos-inf? z)
                 (m/* 0.5 (m/inc s))
                 (m/* 0.5 (m/inc (m/* s (double (prot/cdf dist z))))))))
      :icdf (fn ^double [^double p]
              (if (m/>= p 0.5)
                (let [q (m/dec (m/* 2.0 p))
                      z (double (prot/icdf dist q))]
                  (m/+ mu (m/* alpha (m/pow z inv-beta))))
                (let [q (m/- 1.0 (m/* 2.0 p))
                      z (double (prot/icdf dist q))]
                  (m/- mu (m/* alpha (m/pow z inv-beta))))))
      :rng rng
      :mean mu
      :variance variance
      :dimensions 1
      :continuous? true
      :name :generalized-normal
      :parameters [:mu :alpha :beta :rng]
      :lower-bound ##-Inf
      :upper-bound ##Inf})))

(defn generalized-inverse-gaussian
  "Generalized inverse Gaussian distribution with the density proportional to `x^(lambda - 1) exp(-(chi / x + psi x) / 2)` on `x > 0`.

  Parameters (a map): `chi`, `psi` (positive doubles), `lambda` (double), `rng` (`nil`: a fresh generator). `mean` and `variance` are closed forms with the modified Bessel function of the second kind; `cdf` and `icdf` are obtained by numerical integration of the density on `[0, mean + 25 sd]`, so they are approximate and `cdf` is `1.0` beyond that. Returns a univariate continuous distribution."
  [{:keys [^double chi ^double psi ^double lambda rng]}]
  (let [omega (m/sqrt (m/* chi psi))
        log-kv-lambda (m/log (special/bessel-K lambda omega))
        log-const (m/- (m/* 0.5 lambda (m/- (m/log psi) (m/log chi))) (m/log 2.0) log-kv-lambda)
        scale-est (m/sqrt (m// chi psi))
        lpdf (fn ^double [^double x]
               (cond
                 (m/not-pos? x) ##-Inf
                 (m/pos-inf? x) ##-Inf
                 :else (let [chi-term (m// chi x)]
                         (if (m/pos-inf? chi-term)
                           ##-Inf
                           (m/- (m/+ log-const (m/* (m/dec lambda) (m/log x)))
                                (m/* 0.5 (m/+ chi-term (m/* psi x))))))))
        pdf (fn ^double [^double x] (m/exp (lpdf x)))
        mean-v (double (m/* scale-est (m// (special/bessel-K (m/inc lambda) omega)
                                           (special/bessel-K lambda omega))))
        variance-v (m/- (m/* scale-est scale-est (m// (special/bessel-K (m/+ lambda 2.0) omega)
                                                      (special/bessel-K lambda omega)))
                        (m/* mean-v mean-v))
        mx (m/+ mean-v (m/* 25.0 (m/sqrt variance-v)))
        [cdf-fn icdf-fn] (integrate-pdf pdf {:mn 0.0 :mx mx :steps 2000}) ;; monotone is default
        cdf (fn ^double [^double x]
              (cond
                (m/not-pos? x) 0.0
                (m/>= x mx) 1.0
                :else (m/constrain (double (cdf-fn x)) 0.0 1.0)))
        icdf (fn ^double [^double p]
               (cond
                 (m/not-pos? p) 0.0
                 (m/>= p 1.0) ##Inf
                 :else (double (icdf-fn (m/constrain p 0.0 1.0)))))]
    (->distribution
     {:lpdf lpdf
      :cdf cdf
      :icdf icdf
      :rng rng
      :mean mean-v
      :variance variance-v
      :dimensions 1
      :continuous? true
      :name :generalized-inverse-gaussian
      :parameters [:chi :psi :lambda :rng]
      :lower-bound 0
      :upper-bound ##Inf})))

(defn- generalized-hyperbolic-core
  "Shared machinery for the generalized hyperbolic family: pdf, cdf/icdf (via
  numerical integration of the pdf with `integrate-pdf`), closed-form
  mean/variance, and a sampler (normal-variance mixture: `W` a
  generalized-inverse-gaussian mixing variable, `Z` standard normal). Returns
  a plain map meant to be merged into a `->distribution` options map by the
  caller, which supplies its own `:name`/`:parameters` (and, implicitly,
  `lambda`) - used both by [[generalized-hyperbolic]] itself and by
  [[normal-inverse-gaussian]], its `lambda = -0.5` special case.

  Takes a single map argument (rather than positional primitive args) since
  Clojure functions taking primitives support at most 4 such arguments, and
  this one logically has 5 (`mu`, `delta`, `alpha`, `beta`, `lambda`) plus
  `rng`."
  [{:keys [^double mu ^double delta ^double alpha ^double beta ^double lambda rng]}]
  (let [rng (or rng (JDKRandomGenerator.))
        gam (m/sqrt (m/- (m/* alpha alpha) (m/* beta beta)))
        chi (m/* delta delta)
        psi (m/* gam gam)
        omega (m/* delta gam)
        log-alpha (m/log alpha)
        log-const (m/- (m/* lambda (m/- (m/log gam) (m/log delta)))
                       (m/* 0.5 m/LOG_TWO_PI)
                       (m/log (special/bessel-K lambda omega)))
        nu (m/- lambda 0.5)
        W (generalized-inverse-gaussian {:chi chi :psi psi :lambda lambda :rng rng})
        scale-est (m// delta gam)
        kv0 (special/bessel-K lambda omega)
        kv1 (special/bessel-K (m/inc lambda) omega)
        kv2 (special/bessel-K (m/+ lambda 2.0) omega)
        Ew (m/* scale-est (m// kv1 kv0))
        Varw (m/- (m/* scale-est scale-est (m// kv2 kv0)) (m/* Ew Ew))
        mean-v (m/+ mu (m/* beta Ew))
        variance-v (m/+ Ew (m/* beta beta Varw))
        lpdf (fn ^double [^double x]
               (let [d (m/- x mu)
                     s2 (m/+ chi (m/* d d))
                     s (m/sqrt s2)]
                 (if (m/pos-inf? s) ;; x = +-Infinity (or squaring d overflowed)
                   ##-Inf
                   (m/+ log-const (m/log (special/bessel-K nu (m/* alpha s)))
                        (m/* nu (m/- (m/log s) log-alpha))
                        (m/* beta d)))))
        pdf (fn ^double [^double x] (m/exp (lpdf x)))
        sd (m/sqrt variance-v)
        mn-bound (m/- mean-v (m/* 25.0 sd))
        mx-bound (m/+ mean-v (m/* 25.0 sd))
        [cdf-fn icdf-fn] (integrate-pdf pdf {:mn mn-bound :mx mx-bound :steps 2000}) ;; monotone is default
        cdf (fn ^double [^double x]
              (cond
                (m/<= x mn-bound) 0.0
                (m/>= x mx-bound) 1.0
                :else (m/constrain (double (cdf-fn x)) 0.0 1.0)))
        icdf (fn ^double [^double p]
               (cond
                 (m/<= p 0.0) ##-Inf
                 (m/>= p 1.0) ##Inf
                 :else (double (icdf-fn (m/constrain p 0.0 1.0)))))]
    {:lpdf lpdf
     :cdf cdf
     :icdf icdf
     :sampler (fn ^double [] (let [w (double (prot/sample W))
                                  z (double (prot/grandom rng))]
                              (m/+ mu (m/* beta w) (m/* (m/sqrt w) z))))
     :rng rng
     :mean mean-v
     :variance variance-v
     :dimensions 1
     :continuous? true
     :lower-bound ##-Inf
     :upper-bound ##Inf}))

(defn generalized-hyperbolic
  "Generalized hyperbolic distribution.

  Parameters (a map): `mu` (location), `delta` (scale), `alpha` (tail), `beta` (skewness, `|beta| < alpha`), `lambda` (shape), `rng` (`nil`: a fresh generator). A sample is a normal variance-mean mixture with a generalized inverse Gaussian mixing variable. `mean` and `variance` are closed forms; `cdf` and `icdf` come from the numerical integration of the density on the mean plus or minus 25 standard deviations, so they are approximate. Returns a univariate continuous distribution on the reals."
  [opts]
  (->distribution (assoc (generalized-hyperbolic-core opts)
                         :name :generalized-hyperbolic
                         :parameters [:mu :delta :alpha :beta :lambda :rng])))

(defn normal-inverse-gaussian
  "Normal-inverse Gaussian distribution, in its own `(alpha, beta, mu, delta)`
  parameterization: exactly the `lambda = -0.5` special case of
  [[generalized-hyperbolic]], reusing its fully numerically-integrated
  pdf/cdf/icdf/mean/variance/sampler machinery ([[generalized-hyperbolic-core]])
  rather than the backing SSJ `NormalInverseGaussianDist` class, whose `cdf`
  is not implemented (and whose `sample`, going through `inverseF`/`icdf`,
  therefore always threw regardless of RNG/seeding)."
  [{:keys [alpha beta mu delta rng]}]
  (->distribution (assoc (generalized-hyperbolic-core {:mu mu :delta delta :alpha alpha :beta beta :lambda -0.5 :rng rng})
                         :name :normal-inverse-gaussian
                         :parameters [:alpha :beta :mu :delta :rng])))

(defn half-logistic
  "Half-logistic distribution with the `scale`: the absolute value of a logistic variable, `cdf(x) = (1 - exp(-x / scale)) / (1 + exp(-x / scale))` for `x >= 0`.

  `rng` is the random generator (`nil`: a fresh one). `mean` is `2 log(2) scale`. Returns a univariate continuous distribution on `[0, Inf)`."
  [^double scale rng]
  (let [log-2-scale (m/- (m/log 2.0) (m/log scale))
        r-scale (m// 1.0 scale)
        mean-v (m/* 2.0 m/LN2 scale)
        variance-v (m/* (m/- (m/* m/PI2 m/THIRD) (m/* 4.0 m/LN2 m/LN2)) scale scale)]
    (->distribution {:lpdf (fn ^double [^double x]
                             (if (m/neg? x)
                               ##-Inf
                               (let [e (m/exp (m/- (m/* x r-scale)))]
                                 (m/- log-2-scale (m/* x r-scale) (m/* 2.0 (m/log1p e))))))
                     :cdf (fn ^double [^double x]
                            (if (m/neg? x)
                              0.0
                              (let [e (m/exp (m/- (m/* x r-scale)))]
                                (m// (m/- 1.0 e) (m/+ 1.0 e)))))
                     :icdf (fn ^double [^double p] (m/* 2.0 scale (m/atanh p)))
                     :rng rng
                     :mean mean-v
                     :variance variance-v
                     :dimensions 1
                     :continuous? true
                     :name :half-logistic
                     :parameters [:scale :rng]
                     :lower-bound 0
                     :upper-bound ##Inf})))

(defn generalized-half-logistic
  "Generalized (exponentiated) half-logistic distribution with the shape `alpha` and the rate `lambda`: the half-logistic cdf raised to the power `alpha`.

  `rng` is the random generator (`nil`: a fresh one). `mean` and `variance` are computed by numerical integration on first request. Returns a univariate continuous distribution on `[0, Inf)`."
  [^double alpha ^double lambda rng]
  (let [alpha-1 (m/dec alpha)
        r-alpha (m// 1.0 alpha)
        log-lambda2 (m/+ (m/log lambda) (m/log 2.0))
        log-alpha-lambda2 (m/+ (m/log alpha) log-lambda2)
        lpdf (if (m/one? alpha)
               (fn ^double [^double x]
                 (if (m/neg? x)
                   ##-Inf
                   (let [e (m/exp (m/- (m/* lambda x)))]
                     (m/- log-lambda2 (m/* lambda x) (m/* 2.0 (m/log1p e))))))
               (fn ^double [^double x]
                 (if (m/neg? x)
                   ##-Inf
                   (let [e (m/exp (m/- (m/* lambda x)))
                         log-u (m/- (m/log1p (m/- e)) (m/log1p e))]
                     (m/- (m/+ log-alpha-lambda2 (m/* alpha-1 log-u))
                          (m/* lambda x) (m/* 2.0 (m/log1p e)))))))
        cdf (fn ^double [^double x]
              (if (m/neg? x)
                0.0
                (let [e (m/exp (m/- (m/* lambda x)))
                      u (m// (m/- 1.0 e) (m/+ 1.0 e))]
                  (m/pow u alpha))))
        icdf (fn ^double [^double p] (m// (m/* 2.0 (m/atanh (m/pow p r-alpha))) lambda))
        pdf (fn ^double [^double x] (m/exp (lpdf x)))
        mean (delay (double (quad/gk-quadrature (fn ^double [^double x] (m/* x (double (pdf x)))) 0.0 ##Inf)))
        variance (delay (m/- (double (quad/gk-quadrature (fn ^double [^double x] (m/* x x (double (pdf x)))) 0.0 ##Inf))
                             (m/sq (double @mean))))]
    (->distribution {:lpdf lpdf
                     :cdf cdf
                     :icdf icdf
                     :rng rng
                     :mean mean
                     :variance variance
                     :dimensions 1
                     :continuous? true
                     :name :generalized-half-logistic
                     :parameters [:alpha :lambda :rng]
                     :lower-bound 0
                     :upper-bound ##Inf})))

(defn f-noncentral
  "Noncentral F distribution with `df1` and `df2` degrees of freedom and the noncentrality `ncp`.

  `rng` is the random generator (`nil`: a fresh one). `pdf`, `cdf` and `icdf` are sums over a Poisson mixture of scaled central F distributions (terms are added until their weights sum to `1 - 1e-15`, at most 100000). `mean` is `##Inf` for `df2 <= 2` and `variance` for `df2 <= 4`. Returns a univariate continuous distribution on `[0, Inf)`."
  [^double df1 ^double df2 ^double ncp rng]
  (let [half-ncp (m/* 0.5 ncp)
        eps 1.0e-15
        max-terms 100000
        [ks probs] (loop [k (long 0) pk (m/exp (m/- half-ncp)) cum pk ks [0] probs [pk]]
                     (if (or (m/>= cum (m/- 1.0 eps)) (m/>= k max-terms))
                       [ks probs]
                       (let [k' (m/inc k)
                             pk' (m/* pk (m// half-ncp (double k')))
                             cum' (m/+ cum pk')]
                         (recur k' pk' cum' (conj ks k') (conj probs pk')))))
        ;; [probability, scale=df1/(df1+2k), central-F(df1+2k, df2) dist]
        terms (mapv (fn [^long k ^double p]
                      (let [df1k (m/+ df1 (m/* 2.0 k))]
                        [p (m// df1 df1k) (FDistribution. df1k df2)]))
                    ks probs)
        pdf-at-0 (cond (m/< df1 2.0) ##Inf
                       (m/== df1 2.0) (m/exp (m/- half-ncp))
                       :else 0.0)
        pdf (fn ^double [^double x]
              (cond
                (m/neg? x) 0.0
                (m/zero? x) pdf-at-0
                (m/pos-inf? x) 0.0
                :else (double (reduce (fn [^double acc [^double p ^double scale ^FDistribution d]]
                                        (m/+ acc (m/* p scale (.density d (m/* x scale)))))
                                      0.0 terms))))
        cdf (fn ^double [^double x]
              (cond
                (m/not-pos? x) 0.0
                (m/pos-inf? x) 1.0
                :else (double (reduce (fn [^double acc [^double p ^double scale ^FDistribution d]]
                                        (m/+ acc (m/* p (.cumulativeProbability d (m/* x scale)))))
                                      0.0 terms))))
        init-guess (m/max 1.0 (m/+ 1.0 (m// ncp df1)))
        icdf (fn ^double [^double p]
               (cond
                 (m/not-pos? p) 0.0
                 (m/>= p 1.0) ##Inf
                 :else (icdf-solver cdf p init-guess init-guess)))]
    (->distribution
     {:pdf pdf
      :cdf cdf
      :icdf icdf
      :rng rng
      :mean (delay (if (m/> df2 2.0)
                     (m// (m/* df2 (m/+ df1 ncp)) (m/* df1 (m/- df2 2.0)))
                     ##Inf))
      :variance (delay (if (m/> df2 4.0)
                         (let [d1nc (m/+ df1 ncp)
                               ratio (m// df2 df1)]
                           (m// (m/* 2.0 (m/+ (m/* d1nc d1nc) (m/* (m/+ df1 (m/* 2.0 ncp)) (m/- df2 2.0)))
                                     (m/* ratio ratio))
                                (m/* (m/* (m/- df2 2.0) (m/- df2 2.0)) (m/- df2 4.0))))
                         ##Inf))
      :dimensions 1
      :continuous? true
      :name :f-noncentral
      :parameters [:df1 :df2 :ncp :rng]
      :lower-bound 0
      :upper-bound ##Inf})))

(defn t-noncentral
  "Noncentral Student t distribution with `df` degrees of freedom and the noncentrality `ncp`.

  `rng` is the random generator (`nil`: a fresh one). `pdf` and `cdf` are computed by numerical integration over the chi mixing variable, so they are relatively slow and `icdf` is found by root finding. `mean` is `NaN` for `df <= 1`, `variance` is `NaN` for `df <= 1` and `##Inf` for `1 < df <= 2`. Returns a univariate continuous distribution on the reals."
  [^double df ^double ncp rng]
  (let [delta ncp
        half-df (m/* 0.5 df)
        df-1 (m/dec df)
        scale (m// 1.0 (m/sqrt df))
        log-c (m/- (m/* (m/- 1.0 half-df) m/LN2) (special/log-gamma half-df))
        ;; log-density of W = sqrt(chi-squared(df)), i.e. the chi(df) distribution,
        ;; used as the mixing variable: T = (Z+delta) / (W/sqrt(df)), Z ~ N(0,1)
        chi-lpdf (fn ^double [^double w]
                   (cond
                     (m/neg? w) ##-Inf
                     (m/zero? w) (cond (m/> df 1.0) ##-Inf
                                       (m/== df 1.0) log-c ;; avoids 0 * (-Infinity) -> NaN
                                       :else ##Inf)
                     :else (m/+ (m/* df-1 (m/log w)) (m/* -0.5 w w) log-c)))
        chi-pdf (fn ^double [^double w] (m/exp (chi-lpdf w)))
        normal-pdf (fn ^double [^double x] (m/* m/INV_SQRT2PI (m/exp (m/* -0.5 x x))))
        normal-cdf (fn ^double [^double x] (m/* 0.5 (m/inc (special/erf (m/* x m/INV_SQRT_2)))))
        pdf (fn ^double [^double x]
              (if (or (m/neg-inf? x) (m/pos-inf? x))
                0.0
                (double (quad/gk-quadrature
                         (fn ^double [^double w]
                           (let [z (m/- (m/* x w scale) delta)]
                             (m/* (double (chi-pdf w)) w scale (double (normal-pdf z)))))
                         0.0 ##Inf))))
        cdf (fn ^double [^double x]
              (cond
                (m/neg-inf? x) 0.0
                (m/pos-inf? x) 1.0
                :else (m/constrain
                       (double (quad/gk-quadrature
                                (fn ^double [^double w]
                                  (let [z (m/- (m/* x w scale) delta)]
                                    (m/* (double (chi-pdf w)) (double (normal-cdf z)))))
                                0.0 ##Inf))
                       0.0 1.0)))
        init-guess (m/max 1.0 (m/abs delta))
        icdf (fn ^double [^double p]
               (cond
                 (m/not-pos? p) ##-Inf
                 (m/>= p 1.0) ##Inf
                 :else (icdf-solver cdf p delta init-guess)))
        gamma-ratio (delay (m/exp (m/- (special/log-gamma (m/* 0.5 df-1)) (special/log-gamma half-df))))]
    (->distribution
     {:pdf pdf
      :cdf cdf
      :icdf icdf
      :rng rng
      :mean (delay (if (m/> df 1.0)
                     (m/* delta (m/sqrt (m/* 0.5 df)) (double @gamma-ratio))
                     ##NaN))
      ;; df <= 1: the mean does not exist (NaN), 1 < df <= 2: the second moment diverges (+Inf)
      :variance (delay (cond
                         (m/> df 2.0) (m/- (m// (m/* df (m/inc (m/* delta delta))) (m/- df 2.0))
                                           (m/* delta delta half-df (m/sq (double @gamma-ratio))))
                         (m/> df 1.0) ##Inf
                         :else ##NaN))
      :dimensions 1
      :continuous? true
      :name :t-noncentral
      :parameters [:df :ncp :rng]
      :lower-bound ##-Inf
      :upper-bound ##Inf})))

(defn beta-noncentral
  "Noncentral beta distribution with the shapes `alpha`, `beta` and the noncentrality `ncp`.

  `rng` is the random generator (`nil`: a fresh one). The distribution is a Poisson mixture of central beta distributions with the shape `alpha + k` (terms are added until their weights sum to `1 - 1e-15`, at most 100000), `mean` and `variance` are the exact moments of that mixture and `icdf` is found by root finding. Returns a univariate continuous distribution on `[0, 1]`."
  [^double alpha ^double beta ^double ncp rng]
  (let [half-ncp (m/* 0.5 ncp)
        eps 1.0e-15
        max-terms 100000
        [ks probs] (loop [k (long 0) pk (m/exp (m/- half-ncp)) cum pk ks [0] probs [pk]]
                     (if (or (m/>= cum (m/- 1.0 eps)) (m/>= k max-terms))
                       [ks probs]
                       (let [k' (m/inc k)
                             pk' (m/* pk (m// half-ncp (double k')))
                             cum' (m/+ cum pk')]
                         (recur k' pk' cum' (conj ks k') (conj probs pk')))))
        ;; [probability, shifted alpha1=alpha+k, central Beta(alpha+k, beta) dist]
        terms (mapv (fn [^long k ^double p]
                      (let [alphak (m/+ alpha (double k))]
                        [p alphak (BetaDistribution. alphak beta)]))
                    ks probs)
        ;; Apache Commons' BetaDistribution.density throws at x=0 when its
        ;; alpha<1 (and at x=1 when its beta<1), and silently returns the
        ;; wrong value (0.0 instead of the true finite limit) at x=0/x=1
        ;; when alpha/beta is exactly 1 - so both boundaries are handled
        ;; here directly from the known analytic limits instead of ever
        ;; calling .density there. Only the k=0 term can have alpha<=1 at
        ;; x=0 (every other term's shape is alpha+k>1 since alpha>0); every
        ;; term shares the same `beta` at x=1, so all of them matter there.
        pdf-at-0 (cond (m/> alpha 1.0) 0.0
                       (m/== alpha 1.0) (m/* (double (first probs)) beta) ;; B(1,beta)=1/beta
                       :else ##Inf)
        pdf-at-1 (cond (m/> beta 1.0) 0.0
                       (m/== beta 1.0) (double (reduce (fn [^double acc [^double p ^double alphak _]]
                                                         (m/+ acc (m/* p alphak))) ;; B(alphak,1)=1/alphak
                                                       0.0 terms))
                       :else ##Inf)
        pdf (fn ^double [^double x]
              (cond
                (m/neg? x) 0.0
                (m/> x 1.0) 0.0
                (m/zero? x) pdf-at-0
                (m/== x 1.0) pdf-at-1
                :else (double (reduce (fn [^double acc [^double p _ ^BetaDistribution d]]
                                        (m/+ acc (m/* p (.density d x))))
                                      0.0 terms))))
        cdf (fn ^double [^double x]
              (cond
                (m/not-pos? x) 0.0
                (m/>= x 1.0) 1.0
                :else (double (reduce (fn [^double acc [^double p _ ^BetaDistribution d]]
                                        (m/+ acc (m/* p (.cumulativeProbability d x))))
                                      0.0 terms))))
        icdf (fn ^double [^double p]
               (cond
                 (m/not-pos? p) 0.0
                 (m/>= p 1.0) 1.0
                 ;; the domain is already known to be exactly [0,1] (cdf(0)=0, cdf(1)=1),
                 ;; so root-finding can bracket directly on it, unlike icdf-solver's
                 ;; geometric bracket search (built for unbounded domains).
                 :else (m/constrain (solver/find-root (fn ^double [^double q] (m/- (double (cdf q)) p))
                                                       0.0 1.0 {:absolute-accuracy 1.0e-10})
                                    0.0 1.0)))
        ;; mean/variance are exact weighted sums of the mixture's own central-Beta
        ;; component moments (mean_k=alphak/(alphak+beta), var_k=alphak*beta/((alphak+beta)^2*(alphak+beta+1))),
        ;; not an approximation - just as exact as the truncated pdf/cdf mixture itself.
        mean (delay (double (reduce (fn [^double acc [^double p ^double alphak _]]
                                      (m/+ acc (m/* p (m// alphak (m/+ alphak beta)))))
                                    0.0 terms)))
        ex2 (delay (double (reduce (fn [^double acc [^double p ^double alphak _]]
                                     (let [s (m/+ alphak beta)
                                           mk (m// alphak s)
                                           vk (m// (m/* alphak beta) (m/* s s (m/inc s)))]
                                       (m/+ acc (m/* p (m/+ vk (m/* mk mk))))))
                                   0.0 terms)))]
    (->distribution
     {:pdf pdf
      :cdf cdf
      :icdf icdf
      :rng rng
      :mean mean
      :variance (delay (m/- (double @ex2) (m/sq (double @mean))))
      :dimensions 1
      :continuous? true
      :name :beta-noncentral
      :parameters [:alpha :beta :ncp :rng]
      :lower-bound 0
      :upper-bound 1})))

(defn- von-mises-log-I0
  "log(I_0(kappa)), safe for any kappa: `special/bessel-I0` itself computes
  `exp(kappa)` directly for large kappa and overflows to `##Inf` above ~709
  (I_0(kappa) is then genuinely too large to represent as a plain double, but
  its log is not). Above that threshold, use the standard large-argument
  asymptotic expansion of I_0 directly in log-space instead."
  ^double [^double kappa]
  (if (m/< kappa 700.0)
    (m/log (special/bessel-I0 kappa))
    (let [u (m// 1.0 (m/* 8.0 kappa))
          S (poly/mevalpoly u 1.0 1.0 4.5 37.5 459.375)]
      (m/- (m/+ kappa (m/log S)) (m/* 0.5 (m/log (m/* m/TWO_PI kappa)))))))

(defn von-mises
  "Von Mises distribution on the interval `[mu - pi, mu + pi]` with the mean direction `mu` and the concentration `kappa`, `pdf(x) = exp(kappa cos(x - mu)) / (2 pi I0(kappa))`.

  `rng` is the random generator (`nil`: a fresh one). `cdf` and `icdf` are obtained by numerical integration of the density. `mean` is `mu` and `variance` is the ordinary (linear) variance on the interval, not the circular variance. Returns a univariate continuous distribution."
  [^double mu ^double kappa rng]
  (let [mn (m/- mu m/PI)
        mx (m/+ mu m/PI)
        log-const (m/- 0.0 (von-mises-log-I0 kappa) m/LOG_TWO_PI)
        lpdf (fn ^double [^double x]
               (if (m/<= mn x mx)
                 (m/+ (m/* kappa (m/cos (m/- x mu))) log-const)
                 ##-Inf))
        pdf (fn ^double [^double x] (m/exp (lpdf x)))
        ;; the density concentrates into an ever-narrower peak (width ~1/sqrt(kappa))
        ;; around `mu` as `kappa` grows; a fixed step count eventually becomes wider
        ;; than the peak itself, silently missing it entirely, so scale resolution
        ;; with kappa (steps -> 100*sqrt(kappa+1), clamped to a sane range).
        steps (long (m/constrain (m/ceil (m/* 100.0 (m/sqrt (m/inc kappa)))) 2000 100000))
        [cdf-fn icdf-fn] (integrate-pdf pdf {:mn mn :mx mx :steps steps}) ;; monotone is default
        cdf (fn ^double [^double x]
              (cond
                (m/<= x mn) 0.0
                (m/>= x mx) 1.0
                :else (m/constrain (double (cdf-fn x)) 0.0 1.0)))
        icdf (fn ^double [^double p]
               (cond
                 (m/<= p 0.0) mn
                 (m/>= p 1.0) mx
                 :else (double (icdf-fn (m/constrain p 0.0 1.0)))))
        ;; variance: don't blindly integrate the whole [mu-pi,mu+pi] domain - once the
        ;; peak is much narrower than that range, the adaptive quadrature's own coarse
        ;; initial sampling can miss it entirely and confidently report a near-zero
        ;; (wrong) result instead of erroring. Integrate only a symmetric window sized
        ;; to the peak's own characteristic scale (~40 standard deviations, negligible
        ;; truncation error) and exploit the pdf's exact symmetry around `mu` to halve
        ;; the work.
        half-width (m/min m/PI (m// 40.0 (m/sqrt (m/inc kappa))))]
    (->distribution {:lpdf lpdf
                     :pdf pdf
                     :cdf cdf
                     :icdf icdf
                     :rng rng
                     :mean mu ;; exact, by symmetry of pdf around mu
                     ;; no simple closed form for the (linear) variance restricted to
                     ;; [mu-pi,mu+pi]; the usual circular-statistics analogue is the
                     ;; circular variance 1 - I1(kappa)/I0(kappa), which is NOT what's
                     ;; returned here (this follows the rest of the library's ordinary
                     ;; E[(X-mean)^2] convention for `variance`).
                     :variance (delay (m/* 2.0 (double (quad/gk-quadrature
                                                        (fn ^double [^double t] (m/* t t (double (pdf (m/+ mu t)))))
                                                        0.0 half-width))))
                     :dimensions 1
                     :continuous? true
                     :name :von-mises
                     :parameters [:mu :kappa :rng]
                     :lower-bound mn
                     :upper-bound mx})))

;; ---- ccdf

(defn- normal-ccdf
  ^double [^NormalDistribution d ^double x]
  (m/* 0.5 (Erf/erfc (m// (m/- x (.getMean d)) (m/* (.getStandardDeviation d) m/SQRT2)))))

(defn- laplace-ccdf
  ^double [^LaplaceDistribution d ^double x]
  (let [z (m// (m/- x (.getLocation d)) (.getScale d))]
    (if (m/>= z 0.0)
      (m/* 0.5 (m/exp (m/- z)))
      (m/- 1.0 (m/* 0.5 (m/exp z))))))

(defn- logistic-ccdf
  ^double [^LogisticDistribution d ^double x]
  (let [z (m// (m/- x (.getLocation d)) (.getScale d))]
    (if (m/pos? z)
      (let [e (m/exp (m/- z))]
        (m// e (m/inc e)))
      (m// 1.0 (m/inc (m/exp z))))))

(defn- cauchy-ccdf
  ^double [^CauchyDistribution d ^double x]
  (let [z (m// (m/- x (.getMedian d)) (.getScale d))]
    (if (m/pos? z)
      (m// (m/atan (m// 1.0 z)) m/PI)
      (m/- 0.5 (m// (m/atan z) m/PI)))))

(defn ccdf
  "Returns the complementary cumulative probability `P(X > v)` (the upper tail) of the univariate distribution `d`.

  Parameters:

  - `d`: univariate distribution object
  - `v` (double): the point

  The value is computed with the most accurate method available:

  - distributions backed by an SSJ continuous distribution use its survival function `barF` (the folded normal, whose `barF` is wrong, and the inverse Gaussian, whose `barF` is just `1 - cdf` of an inaccurate SSJ cdf, use `1 - cdf` of the distribution instead),
  - `normal`, `t`, `logistic`, `laplace`, `cauchy`, `exponential`, `weibull`, `pareto`, `log-normal`, `gamma` and `chi-squared` use closed forms (the regularized incomplete gamma function for `gamma` and `chi-squared`), accurate in the far upper tail where `1 - cdf` is exactly `0.0`,
  - any other distribution, discrete ones included, uses `1 - cdf`, which loses all digits once the tail probability is below about `1e-16`.

  Returns `##NaN` for a `NaN` argument. No range check is made: outside the support the result is `1.0` (left of it) or `0.0` (right of it) for the closed forms.

  See also [[truncated]], [[mixture]]."
  ^double [d ^double v]
  (if (m/nan? v)
    ##NaN
    (let [o (prot/source-object d)]
      (cond
        (instance? ContinuousDistribution o) (if (or (instance? FoldedNormalDist o)
                                                     (instance? InverseGaussianDist o))
                                               (m/- 1.0 (double (prot/cdf d v)))
                                               (.barF ^ContinuousDistribution o v))
        (instance? NormalDistribution o) (normal-ccdf o v)
        (instance? TDistribution o) (.cumulativeProbability ^TDistribution o (m/- v))
        (instance? LogisticDistribution o) (logistic-ccdf o v)
        (instance? LaplaceDistribution o) (laplace-ccdf o v)
        (instance? CauchyDistribution o) (cauchy-ccdf o v)
        (instance? ExponentialDistribution o) (if (m/neg? v)
                                                1.0
                                                (m/exp (m/- (m// v (.getMean ^ExponentialDistribution o)))))
        (instance? WeibullDistribution o) (if (m/not-pos? v)
                                            1.0
                                            (let [^WeibullDistribution w o]
                                              (m/exp (m/- (m/pow (m// v (.getScale w)) (.getShape w))))))
        (instance? ParetoDistribution o) (let [^ParetoDistribution p o
                                               scale (.getScale p)]
                                           (if (m/<= v scale)
                                             1.0
                                             (m/pow (m// scale v) (.getShape p))))
        (instance? LogNormalDistribution o) (if (m/not-pos? v)
                                              1.0
                                              (let [^LogNormalDistribution l o]
                                                (m/* 0.5 (Erf/erfc (m// (m/- (m/log v) (.getScale l))
                                                                        (m/* (.getShape l) m/SQRT2))))))
        (instance? GammaDistribution o) (if (m/not-pos? v)
                                          1.0
                                          (let [^GammaDistribution g o]
                                            (if (m/pos-inf? v)
                                              0.0
                                              (Gamma/regularizedGammaQ (.getShape g) (m// v (.getScale g))))))
        (instance? ChiSquaredDistribution o) (if (m/not-pos? v)
                                               1.0
                                               (if (m/pos-inf? v)
                                                 0.0
                                                 (Gamma/regularizedGammaQ (m/* 0.5 (.getDegreesOfFreedom ^ChiSquaredDistribution o))
                                                                          (m/* 0.5 v))))
        :else (m/- 1.0 (double (prot/cdf d v)))))))

;; ---- truncated and mixture

(def ^{:private true :const true :tag 'double} truncated-rejection-threshold 1.0e-3)
(def ^{:private true :const true :tag 'double} truncated-moments-threshold 0.1)
(def ^{:private true :const true :tag 'long} truncated-max-rejections 10000000)

;; discrete distributions without integer support, `truncated` and `mixture` treat integers only
(def ^:private non-integer-discrete-ids
  #{:enumerated-real :real-discrete-distribution :categorical-distribution})

(defn- cdf-at
  "cdf of `distr` which is also defined at infinities."
  ^double [distr ^double x]
  (cond
    (m/neg-inf? x) 0.0
    (m/pos-inf? x) 1.0
    (m/>= x (double (prot/upper-bound distr))) 1.0
    :else (double (prot/cdf distr x))))

(defn- ccdf-at
  "ccdf of `distr` which is also defined at infinities."
  ^double [distr ^double x]
  (cond
    (m/neg-inf? x) 1.0
    (m/pos-inf? x) 0.0
    (m/>= x (double (prot/upper-bound distr))) 0.0
    :else (ccdf distr x)))

(defn- interval-mass
  "Probability `P(a < X <= b)` of `distr` computed from the half of the distribution which is accurate: from cdf if `cdf(a) <= 0.5` or from ccdf otherwise."
  ^double [distr ^double a ^double b]
  (let [ca (cdf-at distr a)]
    (if (m/<= ca 0.5)
      (m/- (cdf-at distr b) ca)
      (m/- (ccdf-at distr a) (ccdf-at distr b)))))

(defn- truncation-bound
  "Checks a user supplied bound, returns double or nil."
  [property v]
  (cond
    (nil? v) nil
    (not (number? v)) (invalid-input! 'truncated property "has to be a number or nil." v "number or nil")
    (m/nan? (double v)) (invalid-input! 'truncated property "can not be NaN." v "number or nil")
    :else (double v)))

(defn- expanded-bracket
  "Finite bracket `[a b]` of the increasing function `f` with `f(a) <= 0 <= f(b)` (if exists) for the interval `[lb, ub]` with possibly infinite bounds."
  [f ^double lb ^double ub]
  (let [f (fn [^double x] (double (f x)))
        a (if (m/valid-double? lb)
            lb
            (let [start (if (m/valid-double? ub) ub 0.0)]
              (loop [step 1.0]
                (let [x (m/- start step)]
                  (if (or (m/<= (double (f x)) 0.0) (m/> step 1.0e300)) x (recur (m/* step 2.0)))))))
        b (if (m/valid-double? ub)
            ub
            (let [start (if (m/valid-double? lb) lb 0.0)]
              (loop [step 1.0]
                (let [x (m/+ start step)]
                  (if (or (m/>= (double (f x)) 0.0) (m/> step 1.0e300)) x (recur (m/* step 2.0)))))))]
    [a b]))

(defn- increasing-root
  "Root of the increasing function `f` on `[lb, ub]` (infinite bounds allowed). Returns the end of the bracket if `f` doesn't change sign."
  ^double [f ^double lb ^double ub]
  (let [[^double a ^double b] (expanded-bracket f lb ub)
        fa (double (f a))
        fb (double (f b))]
    (cond
      (m/>= fa 0.0) a
      (m/<= fb 0.0) b
      :else (solver/find-root f a b {:absolute-accuracy (m/* 1.0e-15 (m/max 1.0 (m/abs a) (m/abs b)))
                                     :relative-accuracy 1.0e-14
                                     :max-iters 500}))))

(defn- truncated-continuous-moments
  "Returns `[mean variance]` delays of the continuous distribution truncated to `[lb, ub]` with the mass `z`."
  [distr lb ub z log-z]
  (let [lb (double lb)
        ub (double ub)
        z (double z)
        log-z (double log-z)]
    (if (m/>= z truncated-moments-threshold)
      (let [mean (delay (m// (double (quad/gk-quadrature (fn ^double [^double x] (m/* x (double (prot/pdf distr x))))
                                                         lb ub))
                             z))]
        [mean
         (delay (m/- (m// (double (quad/gk-quadrature (fn ^double [^double x] (m/* x x (double (prot/pdf distr x))))
                                                      lb ub))
                          z)
                     (m/sq (double @mean))))])
      ;; small mass: the (normalized) truncated density is integrated, an integral of the density of the
      ;; base distribution would be dominated by the absolute tolerance of the quadrature
      (let [opts {:rel 1.0e-12 :abs 1.0e-14}
            tpdf (fn [^double x] (m/exp (m/- (double (prot/lpdf distr x)) log-z)))
            shift (double (cond (m/valid-double? lb) lb (m/valid-double? ub) ub :else 0.0))
            mean (delay (m/+ shift (double (quad/gk-quadrature (fn ^double [^double x] (m/* (m/- x shift) (double (tpdf x))))
                                                               lb ub opts))))]
        [mean
         (delay (let [mu (double @mean)]
                  (double (quad/gk-quadrature (fn ^double [^double x] (m/* (m/sq (m/- x mu)) (double (tpdf x))))
                                              lb ub opts))))]))))

(defn- truncated-continuous
  [distr l r]
  (let [base-lb (double (prot/lower-bound distr))
        base-ub (double (prot/upper-bound distr))
        lower-bound (if l (m/max (double l) base-lb) base-lb)
        upper-bound (if r (m/min (double r) base-ub) base-ub)
        _ (when-not (m/< lower-bound upper-bound)
            (invalid-input! 'truncated :interval
                            "the interval has to be non-empty and has to intersect the support of the distribution."
                            [(or l base-lb) (or r base-ub)] "left < right inside the support"))
        left-cdf (cdf-at distr lower-bound)
        right-cdf (cdf-at distr upper-bound)
        left-ccdf (ccdf-at distr lower-bound)
        right-ccdf (ccdf-at distr upper-bound)
        upper? (m/> left-cdf 0.5) ;; the interval is in the upper half, use the survival function
        z (if upper? (m/- left-ccdf right-ccdf) (m/- right-cdf left-cdf))
        _ (when-not (m/> z 0.0)
            (invalid-input! 'truncated :interval
                            "the interval has zero probability under the distribution."
                            [lower-bound upper-bound] "interval with positive probability"))
        log-z (m/log z)
        [mean variance] (truncated-continuous-moments distr lower-bound upper-bound z log-z)
        tcdf (fn [^double x]
               (cond
                 (m/nan? x) ##NaN
                 (m/<= x lower-bound) 0.0
                 (m/>= x upper-bound) 1.0
                 upper? (m/constrain (m// (m/- left-ccdf (ccdf-at distr x)) z) 0.0 1.0)
                 :else (m/constrain (m// (m/- (cdf-at distr x) left-cdf) z) 0.0 1.0)))
        tsurv (fn [^double x]
                (cond
                  (m/nan? x) ##NaN
                  (m/<= x lower-bound) 1.0
                  (m/>= x upper-bound) 0.0
                  :else (m/constrain (m// (interval-mass distr x upper-bound) z) 0.0 1.0)))
        icdf (fn ^double [^double p]
               (cond
                 (m/nan? p) ##NaN
                 (m/<= p 0.0) lower-bound
                 (m/>= p 1.0) upper-bound
                 (m/<= p 0.5) (increasing-root (fn ^double [^double x] (m/- (double (tcdf x)) p)) lower-bound upper-bound)
                 :else (let [q (m/- 1.0 p)]
                         (increasing-root (fn ^double [^double x] (m/- q (double (tsurv x)))) lower-bound upper-bound))))]
    (->distribution {:lpdf (fn ^double [^double x] (if (m/<= lower-bound x upper-bound)
                                                    (m/- (double (prot/lpdf distr x)) log-z)
                                                    ##-Inf))
                     :cdf tcdf
                     :icdf icdf
                     :sampler (if (m/>= z truncated-rejection-threshold)
                                (fn ^double [] (let [v (double (prot/sample distr))]
                                                 (if (m/<= lower-bound v upper-bound) v (recur))))
                                ;; probability integral transform: the base sample, mapped by its own cdf,
                                ;; is a uniform number to feed the inverse cdf of the truncated distribution
                                (fn ^double [] (icdf (cdf-at distr (double (prot/sample distr))))))
                     :rng distr
                     :dimensions 1
                     :continuous? true
                     :name :truncated
                     :parameters [:distr :left :right :rng]
                     :mean mean
                     :variance variance
                     :lower-bound lower-bound
                     :upper-bound upper-bound})))

(defn- discrete-moments
  "Mean and variance of the (unnormalized) integer pmf of `distr` on `[lo, hi]`, summed until the terms fall below `1e-16` of the sum on the decreasing side."
  [distr ^double lo ^double hi]
  (loop [k lo
         m0 0.0
         m1 0.0
         m2 0.0
         prev 0.0]
    (let [p (double (prot/pdf distr k))
          d (m/- k lo)
          m0 (m/+ m0 p)
          m1 (m/+ m1 (m/* d p))
          m2 (m/+ m2 (m/* d d p))]
      (if (or (m/>= k hi)
              (m/> d 1.0e9)
              (and (m/> m0 0.0) (m/<= p prev) (m/<= p (m/* 1.0e-16 m0))))
        (let [mu (m// m1 m0)]
          [(m/+ lo mu) (m/max 0.0 (m/- (m// m2 m0) (m/* mu mu)))])
        (recur (m/inc k) m0 m1 m2 p)))))

(defn- discrete-icdf
  "Smallest integer `k` from `[lo, hi]` for which `ok?` is true, `ok?` is monotone and true at `hi`."
  [ok? ^double lo ^double hi]
  (if (ok? lo)
    (long lo)
    ;; exponential search for the end of the bracket, then bisection
    (loop [a lo
           step 1.0]
      (let [b (m/min hi (m/+ a step))]
        (cond
          (ok? b) (loop [a a
                         b b]
                    (if (m/<= (m/- b a) 1.0)
                      (long b)
                      (let [mid (m/floor (m/* 0.5 (m/+ a b)))]
                        (if (ok? mid) (recur a mid) (recur mid b)))))
          (m/>= b hi) (if (m/valid-double? hi) (long hi) hi)
          :else (recur b (m/* step 2.0)))))))

(def ^{:private true :const true :tag 'long} truncated-table-limit 20000)

(defn- discrete-sampler
  "Rejection sampler of the integers from `[lo, hi]` drawn from `distr`, limited by `truncated-max-rejections` iterations."
  [distr lo hi z]
  (let [lo (double lo)
        hi (double hi)]
    (fn []
      (loop [i 0]
        (let [v (prot/sample distr)
              x (double v)]
          (cond
            (m/<= lo x hi) v
            (m/>= i truncated-max-rejections)
            (throw (ex-info "truncated: the rejection sampler reached the limit of iterations, the interval has too small probability."
                            {:type :sampling-failure
                             :function 'truncated
                             :iterations truncated-max-rejections
                             :interval [lo hi]
                             :probability z}))
            :else (recur (m/inc i))))))))

(defn- truncated-discrete-table
  "Truncated integer distribution with the finite support `[lo, hi]` (not too long): the pmf is tabulated and normalised by its own sum, so pmf, cdf and icdf are consistent exactly."
  [distr ^double lo ^double hi]
  (let [n (int (m/inc (m/- hi lo)))
        raw (double-array n)
        _ (dotimes [i n] (aset raw i (double (prot/pdf distr (m/+ lo i)))))
        z (double (loop [i 0 s 0.0]
                    (if (m/< i n) (recur (m/inc i) (m/+ s (aget raw i))) s)))
        _ (when-not (m/> z 0.0)
            (invalid-input! 'truncated :interval "the interval has zero probability under the distribution."
                            [lo hi] "interval with positive probability"))
        pmf (double-array n)
        cum (double-array n)
        sv (double-array n)
        _ (loop [i 0 s 0.0]
            (when (m/< i n)
              (let [r (aget raw i)
                    s (m/+ s r)]
                (aset pmf i (m// r z))
                (aset cum i (m// s z))
                (recur (m/inc i) s))))
        _ (aset cum (m/dec n) 1.0)
        _ (loop [i (m/dec n) s 0.0]
            (when (m/>= i 0)
              (aset sv i (m// s z))
              (recur (m/dec i) (m/+ s (aget raw i)))))
        mean-offset (double (loop [i 0 s 0.0]
                              (if (m/< i n) (recur (m/inc i) (m/+ s (m/* i (aget pmf i)))) s)))
        variance (double (loop [i 0 s 0.0]
                           (if (m/< i n) (recur (m/inc i) (m/+ s (m/* (m/sq (m/- i mean-offset)) (aget pmf i)))) s)))
        log-z (m/log z)
        index (fn [^double x] (if (and (m/<= lo x hi) (m/== x (m/floor x))) (unchecked-long (m/- x lo)) -1))
        first-index (fn [pred] ;; smallest i with (pred i); pred is monotone and true at n - 1
                      (loop [a 0 b (m/dec n)]
                        (if (m/>= a b)
                          b
                          (let [mid (m/+ a (bit-shift-right (m/- b a) 1))]
                            (if (pred mid) (recur a mid) (recur (m/inc mid) b))))))]
    {:pdf (fn ^double [^double x] (let [i (long (index x))] (if (m/neg? i) 0.0 (aget pmf i))))
     :lpdf (fn ^double [^double x] (let [i (long (index x))] (if (m/neg? i) ##-Inf (m/log (aget pmf i)))))
     :cdf (fn [^double x]
            (cond
              (m/nan? x) ##NaN
              (m/< x lo) 0.0
              (m/>= x hi) 1.0
              :else (aget cum (unchecked-int (m/- (m/floor x) lo)))))
     :icdf (fn [^double p]
             (cond
               (m/nan? p) ##NaN
               (m/<= p 0.0) (long lo)
               (m/>= p 1.0) (long hi)
               (m/<= p 0.5) (m/+ (long lo) (long (first-index (fn [^long i] (m/>= (aget cum i) p)))))
               :else (let [q (m/- 1.0 p)]
                       (m/+ (long lo) (long (first-index (fn [^long i] (m/<= (aget sv i) q))))))))
     :mean (m/+ lo mean-offset)
     :variance variance
     :z z
     :log-z log-z}))

(defn- truncated-discrete-sums
  "Truncated integer distribution with a long or infinite support `[lo, hi]`: cdf based normalisation, moments by summation."
  [distr ^double lo ^double hi]
  (let [before (m/dec lo)
        before-cdf (cdf-at distr before)
        before-ccdf (ccdf-at distr before)
        upper? (m/> before-cdf 0.5)
        z (if upper?
            (m/- before-ccdf (ccdf-at distr hi))
            (m/- (cdf-at distr hi) before-cdf))
        _ (when-not (m/> z 0.0)
            (invalid-input! 'truncated :interval "the interval has zero probability under the distribution."
                            [lo hi] "interval with positive probability"))
        log-z (m/log z)
        moments (delay (discrete-moments distr lo hi))
        tcdf (fn [^double x]
               (cond
                 (m/nan? x) ##NaN
                 (m/< x lo) 0.0
                 (m/>= x hi) 1.0
                 :else (let [k (m/floor x)]
                         (m/constrain (m// (if upper?
                                             (m/- before-ccdf (ccdf-at distr k))
                                             (m/- (cdf-at distr k) before-cdf))
                                           z)
                                      0.0 1.0))))
        tsurv (fn [^double k]
                (m/constrain (m// (interval-mass distr k hi) z) 0.0 1.0))
        in-range? (fn [^double x] (and (m/<= lo x hi) (m/== x (m/floor x))))]
    {:pdf (fn ^double [^double x] (if (in-range? x) (m// (double (prot/pdf distr x)) z) 0.0))
     :lpdf (fn ^double [^double x] (if (in-range? x) (m/- (double (prot/lpdf distr x)) log-z) ##-Inf))
     :cdf tcdf
     :icdf (fn [^double p]
             (cond
               (m/nan? p) ##NaN
               (m/<= p 0.0) (long lo)
               (m/>= p 1.0) (if (m/valid-double? hi) (long hi) hi)
               (m/<= p 0.5) (discrete-icdf (fn [^double k] (m/>= (double (tcdf k)) p)) lo hi)
               :else (let [q (m/- 1.0 p)]
                       (discrete-icdf (fn [^double k] (m/<= (double (tsurv k)) q)) lo hi))))
     :mean (delay (first @moments))
     :variance (delay (second @moments))
     :z z
     :log-z log-z}))

(defn- truncated-discrete
  [distr l r]
  (when (contains? non-integer-discrete-ids (prot/distribution-id distr))
    (invalid-input! 'truncated :distr "only distributions on integers can be truncated, this one is discrete on arbitrary reals."
                    (prot/distribution-id distr) "continuous or integer valued distribution"))
  (let [base-lb (double (prot/lower-bound distr))
        base-ub (double (prot/upper-bound distr))
        lo (m/ceil (if l (m/max (double l) base-lb) base-lb))
        hi (m/floor (if r (m/min (double r) base-ub) base-ub))
        _ (when-not (m/valid-double? lo)
            (invalid-input! 'truncated :left "the lower bound of a discrete distribution has to be finite." l "finite number"))
        _ (when-not (m/<= lo hi)
            (invalid-input! 'truncated :interval
                            "the interval has to contain an integer from the support of the distribution."
                            [(or l base-lb) (or r base-ub)] "interval with an integer inside the support"))
        {:keys [pdf lpdf cdf icdf mean variance z]} (if (m/< (m/- hi lo) truncated-table-limit)
                                                      (truncated-discrete-table distr lo hi)
                                                      (truncated-discrete-sums distr lo hi))]
    (->distribution {:pdf pdf
                     :lpdf lpdf
                     :cdf cdf
                     :icdf icdf
                     :sampler (discrete-sampler distr lo hi z)
                     :rng distr
                     :dimensions 1
                     :continuous? false
                     :name :truncated
                     :parameters [:distr :left :right :rng]
                     :mean mean
                     :variance variance
                     :lower-bound lo
                     :upper-bound hi})))

(defn truncated
  "Truncates the distribution `distr` to `[left, right]`, see `fastmath.random/truncated` for the behaviour.

  Parameters:

  - `distr`: a univariate distribution, continuous or discrete on integers
  - `left`, `right`: numbers or `nil` (the bound of `distr`)

  Returns the truncated distribution, throws `ex-info` (`:type :invalid-input`) for invalid arguments or an interval with zero probability."
  [distr left right]
  (when-not (prot/distribution? distr)
    (invalid-input! 'truncated :distr "has to be a distribution." distr "distribution"))
  (when-not (== 1 (long (prot/dimensions distr)))
    (invalid-input! 'truncated :distr "has to be a univariate distribution." distr "univariate distribution"))
  (let [l (truncation-bound :left left)
        r (truncation-bound :right right)]
    (when (and l r (m/> (double l) (double r)))
      (invalid-input! 'truncated :interval "left can not be greater than right." [left right] "left <= right"))
    (if (prot/continuous? distr)
      (truncated-continuous distr l r)
      (truncated-discrete distr l r))))

(defn- sequence-input?
  [x]
  (or (sequential? x) (and (some? x) (.isArray (class x)))))

(defn- bisect-smallest
  "Smallest `x` in `[lb, ub]` (up to the tolerance `1e-12`, relative or absolute) with `(ok? x)` for a monotone predicate, starting from the bracket `[a0, b0]` which is expanded (doubling the width) if needed."
  [ok? lb ub a0 b0]
  (let [lb (double lb)
        ub (double ub)
        clip (fn [^double x] (m/max -1.0e308 (m/min 1.0e308 (m/max lb (m/min ub x)))))
        a1 (double (clip (double a0)))
        b1 (double (clip (double b0)))
        w (m/max 1.0 (m/abs a1))
        a (if (m/< a1 b1) a1 (double (clip (m/- a1 w))))
        b (if (m/< a1 b1) b1 (double (clip (m/+ b1 w))))]
    (loop [phase 0
           a a
           b b
           i 0]
      (case phase
        ;; 0: `ok?` has to be false at `a`
        0 (cond
            (not (ok? a)) (recur 1 a b 0)
            (or (m/<= a lb) (m/<= a -1.0e308)) a
            :else (recur 0 (double (clip (m/- a (m/* 2.0 (m/- b a))))) a 0))
        ;; 1: `ok?` has to be true at `b`
        1 (cond
            (ok? b) (recur 2 a b 0)
            (or (m/>= b ub) (m/>= b 1.0e308)) b
            :else (recur 1 b (double (clip (m/+ b (m/* 2.0 (m/- b a))))) 0))
        ;; 2: bisection
        (if (or (m/<= (m/- b a) (m/* 1.0e-12 (m/max 1.0 (m/abs a) (m/abs b)))) (m/>= i 2000))
          b
          (let [mid (m/+ a (m/* 0.5 (m/- b a)))]
            (if (ok? mid)
              (recur 2 a mid (m/inc i))
              (recur 2 mid b (m/inc i)))))))))

(defn mixture
  "Finite mixture of the univariate distributions `distrs` with the `weights`, see `fastmath.random/mixture` for the behaviour.

  Parameters:

  - `distrs`: non-empty sequence of univariate distributions
  - `weights`: sequence of finite non-negative numbers of the same count, with a positive sum
  - `rng`: random generator which selects the component (`nil`: a fresh one)

  Returns the mixture, throws `ex-info` (`:type :invalid-input`) for invalid arguments."
  [distrs weights rng]
  (when-not (sequence-input? distrs)
    (invalid-input! 'mixture :distrs "has to be a sequence of distributions." distrs "sequence of distributions"))
  (when-not (sequence-input? weights)
    (invalid-input! 'mixture :weights "has to be a sequence of numbers." weights "sequence of numbers"))
  (let [all-ds (vec distrs)
        all-ws (vec weights)
        n (count all-ds)]
    (when (zero? n)
      (invalid-input! 'mixture :distrs "can not be empty." distrs "non-empty sequence of distributions"))
    (doseq [d all-ds]
      (when-not (and (prot/distribution? d) (== 1 (long (prot/dimensions d))))
        (invalid-input! 'mixture :distrs "all elements have to be univariate distributions." d "univariate distribution")))
    (when-not (== n (count all-ws))
      (invalid-input! 'mixture :weights "the number of weights has to be equal to the number of distributions."
                      (count all-ws) n))
    (doseq [w all-ws]
      (when-not (and (number? w) (m/valid-double? (double w)) (m/not-neg? (double w)))
        (invalid-input! 'mixture :weights "all weights have to be finite and non-negative numbers." w "finite number >= 0")))
    (when-not (let [s (double (reduce m/+ 0.0 (map double all-ws)))]
                (and (m/valid-double? s) (m/pos? s)))
      (invalid-input! 'mixture :weights "the sum of weights has to be finite and positive." all-ws "positive finite sum"))
    (let [used (filterv (fn [i] (m/pos? (double (nth all-ws i)))) (range n))
          distrs (mapv all-ds used)
          weights (mapv (fn [i] (double (nth all-ws i))) used)
          rng (or rng (JDKRandomGenerator.))
          probs (v/normalize-L1 weights)
          enum (categorical-distribution distrs probs rng)
          all-integer? (every? (fn [d] (and (not (prot/continuous? d))
                                        (not (contains? non-integer-discrete-ids (prot/distribution-id d)))))
                           distrs)
          lower-bound (double (reduce m/min (map (fn [d] (double (prot/lower-bound d))) distrs)))
          upper-bound (double (reduce m/max (map (fn [d] (double (prot/upper-bound d))) distrs)))
          cdf (fn [^double x] (v/dot (map (fn ^double [d] (prot/cdf d x)) distrs) probs))
          sv (fn [^double x] (v/dot (map (fn ^double [d] (ccdf d x)) distrs) probs))
          mean (delay (v/dot (map prot/mean distrs) probs))
          variance (delay (m/- (v/dot (v/add (v/sq (map prot/mean distrs))
                                             (map prot/variance distrs)) probs) (m/sq @mean)))
          icdf-fn (fn [^double p]
                    (cond
                      (m/nan? p) ##NaN
                      (m/<= p 0.0) (if all-integer? (long lower-bound) lower-bound)
                      (m/>= p 1.0) (if (and all-integer? (m/valid-double? upper-bound)) (long upper-bound) upper-bound)
                      (and all-integer? (m/valid-double? lower-bound))
                      (discrete-icdf (fn [^double k] (m/>= (double (cdf k)) p)) lower-bound upper-bound)
                      :else
                      (let [qs (filter (fn [^double q] (m/valid-double? q))
                                       (map (fn [d] (try (double (prot/icdf d p)) (catch Exception _ ##NaN))) distrs))
                            a (if (seq qs) (double (reduce m/min qs)) 0.0)
                            b (if (seq qs) (double (reduce m/max qs)) 0.0)]
                        (bisect-smallest (if (m/<= p 0.5)
                                           (fn [^double x] (m/>= (double (cdf x)) p))
                                           ;; the upper half is searched on the survival function (accurate tails)
                                           (let [q (m/- 1.0 p)]
                                             (fn [^double x] (m/<= (double (sv x)) q))))
                                         lower-bound upper-bound a b))))
          select (fn [] (prot/sample (prot/sample enum)))]
      (->distribution {:pdf (fn ^double [^double x] (v/dot (map (fn ^double [d] (prot/pdf d x)) distrs) probs))
                       :cdf cdf
                       :icdf icdf-fn
                       :sampler (if all-integer?
                                  select
                                  (fn ^double [] (select)))
                       :seeder (fn [^long seed]
                                 ;; the selector gets the seed, the i-th component the i-th long of a generator seeded with it
                                 (prot/set-seed! rng seed)
                                 (let [^MersenneTwister g (MersenneTwister. seed)]
                                   (doseq [d distrs]
                                     (prot/set-seed! d (.nextLong g)))))
                       :dimensions 1
                       :rng rng
                       :continuous? (boolean (some prot/continuous? distrs))
                       :name :mixture
                       :parameters [:distrs :weights :rng]
                       :mean mean
                       :variance variance
                       :lower-bound lower-bound
                       :upper-bound upper-bound}))))
