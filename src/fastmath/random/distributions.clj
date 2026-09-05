(ns fastmath.random.distributions
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
  (:import [org.apache.commons.math3.random RandomGenerator EmpiricalDistribution JDKRandomGenerator]
           [org.apache.commons.math3.distribution AbstractRealDistribution RealDistribution BetaDistribution CauchyDistribution ChiSquaredDistribution ConstantRealDistribution EnumeratedRealDistribution ExponentialDistribution FDistribution GammaDistribution, GumbelDistribution, LaplaceDistribution, LevyDistribution, LogisticDistribution, LogNormalDistribution, NakagamiDistribution, NormalDistribution, ParetoDistribution, TDistribution, TriangularDistribution, UniformRealDistribution WeibullDistribution MultivariateNormalDistribution]
           [org.apache.commons.math3.distribution IntegerDistribution AbstractIntegerDistribution BinomialDistribution EnumeratedIntegerDistribution, GeometricDistribution, HypergeometricDistribution, PascalDistribution, PoissonDistribution, UniformIntegerDistribution, ZipfDistribution]
           [org.apache.commons.math3.stat StatUtils]

           [umontreal.ssj.probdist ContinuousDistribution]
           [umontreal.ssj.probdistmulti DirichletDist MultinomialDist]
           
           [fastmath.java Array]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

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

(defn multinomial
  [^long n ps binomial rng]
  (let [rng (or rng (JDKRandomGenerator.))
        ps (m/seq->double-array (v/div ps (v/sum ps)))
        
        m (delay (vec (MultinomialDist/getMean n ps)))
        cv (delay (mapv vec (MultinomialDist/getCovariance n ps)))
        dim (count ps)
        binom-probs (mapv (fn [^double prob ^double sum]
                            (m// prob (m/- 1.0 sum))) ps (reductions m/+ 0.0 ps))]
    (reify
      prot/DistributionProto
      (pdf [_ v] (MultinomialDist/prob n ps (int-array v)))
      (lpdf [_ v] (m/log (MultinomialDist/prob n ps (int-array v))))
      (probability [_ v] (MultinomialDist/prob n ps (int-array v)))
      (cdf [_ v] (MultinomialDist/cdf n ps (int-array v)))
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

(defn- dirichlet-lpdf
  ^double [^long dim alpha- values ^double lbeta]
  (when (m/not== dim (count values)) (throw (ex-info "Invalid input size." {:is (count values) :expected dim})))
  (let [v (reduce m/+ lbeta (map (fn [^double ai ^double x] 
                                   (m/* ai (m/log x))) alpha- values))]
    (if (m/invalid-double? v) ##-Inf v)))

(defn dirichlet
  [alpha gamma rng]
  (let [rng (or rng (JDKRandomGenerator.))
        alpha (if (integer? alpha)
                (double-array alpha 1.0)
                (m/seq->double-array alpha))
        sampler (mapv (fn [^double shape] (gamma {:shape shape :scale 1.0 :rng rng})) alpha)

        lbeta (dirichlet-rev-log-beta alpha)
        alpha- (map m/dec alpha)
        
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
                            (let [a (double-array dim)] ;; any position
                              (aset ^doubles a (int (prot/irandom rng dim)) 1.0)
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
      (set-seed! [d seed] (prot/set-seed! rng seed) d))))

;;

(defn- diff-cdf
  ^double [cdf ^double v1 ^double v2]
  (m/- (double (cdf v2)) (double (cdf v1))))

(defn ->distribution
  [{:keys [pdf lpdf cdf icdf rng ^long dimensions continuous? name parameters mean variance lower-bound upper-bound sampler seeder]
    :or {dimensions 1}}]
  (let [rng (or rng (JDKRandomGenerator.))
        pdf (or pdf (fn ^double [^double v] (m/exp (double (lpdf v)))))
        lpdf (or lpdf (fn ^double [^double v] (m/log (double (pdf v)))))
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

(defn integrate-pdf
  "Integrate PDF function, returns CDF and iCDF

  Parameters:
  * `pdf-func` - univariate function
  * `mn` - lower bound for integration, value of pdf-func should be 0.0 at this point
  * `mx` - upper bound for integration
  * `steps` - how much subintervals to integrate (default 1000)
  * `interpolator` - interpolation method between integrated points (default :linear)

  Also other integration related parameters are accepted (`:gauss-kronrod` integration is used).

  Possible interpolation methods: `:linear`, `:cubic`, `:monotone` (default) or any function from `fastmath.interpolation`"
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
         ys (v/div ys (Array/aget ys (dec steps))) ;; normalize to ensure 1 at the endpoint
         intpol (case interpolator
                  :monotone monotone-interp/monotone
                  :linear linear-interp/linear
                  :cubic cubic-interp/cubic
                  (if (fn? interpolator) interpolator linear-interp/linear))]
     [(let [i (intpol xs ys)] (fn [^double x] (m/constrain (double (i x)) 0.0 1.0)))
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
                     :variance (delay (StatUtils/variance (m/seq->double-array data) (double @m)))
                     :lower-bound (m/- mn step)
                     :upper-bound (m/+ mx step)})))

;;

(defn discrete-binary-search
  ([cdf-fn ^double p [^long mid ^long step]] (discrete-binary-search cdf-fn (m/long-max 1 step) p [0 mid]))
  ([cdf-fn ^long step ^double p [^long mn ^long mx]]
   (cond
     (m/< (double (cdf-fn mx)) p) (recur cdf-fn (m/* 2 step) p [mx (m/+ mx step)])
     (m/<= (m/- mx mn) 1) (if (m/>= (double (cdf-fn mn)) p) mn mx)
     :else (let [mid (m// (m/+ mn mx) 2)]
             (if (m/> (double (cdf-fn mid)) p)
               (recur cdf-fn step p [mn mid])
               (recur cdf-fn step p [mid mx]))))))

(defn negative-binomial
  [^double r ^double p rng]
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
                     :parameters [:r :p :rng]
                     :mean mean
                     :variance variance
                     :lower-bound 0
                     :upper-bound Integer/MAX_VALUE})))

(defn logarithmic
  [^double p rng]
  (let [p- (m/- 1.0 p)
        logp- (m/log p-)
        r (m// logp-)
        mean (m/* -1.0 r (m// p p-))
        variance (m/- (m// (m/+ (m/* p p) (m/* p logp-))
                           (m/sq (m/* p- logp-))))
        cdf (fn [^long k]
              (m/inc (m/* r (special/incomplete-beta p (m/inc k) 0.0))))]
    (->distribution {:pdf (fn [^long k]
                            (if (m/< k 1) ##NaN
                                (m/* -1.0 r (m// (m/fpow p k) k))))
                     :cdf cdf
                     :icdf (fn [^double p]
                             (cond
                               (m/not-pos? p) 1
                               (m/>= p 1.0) ##Inf
                               :else (discrete-binary-search cdf p [(long mean) (long (m/sqrt variance))])))
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
                     :name :half-cachy
                     :parameters [:mu :scale :rng]
                     :mean ##NaN
                     :variance ##NaN
                     :lower-bound 0.0
                     :upper-bound ##Inf})))

;;

(defn integer-discrete-distribution
  [data probabilities rng]
  (let [cnt (count data)
        probabilities (or probabilities (repeat cnt 1))
        sum (v/sum probabilities)
        pmf (->> (map vector data probabilities)
                 (reduce (fn [m [v ^double p]]
                           (im/update m v (fnil m/+ 0.0) (m// p sum))) (im/int-map)))
        probs (vals pmf)
        cumsum (reductions m/+ probs)
        ks (keys pmf)
        mnk (double (first ks))
        step-before (step-interp/step-before cumsum ks)
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
  [data probabilities rng]
  (let [cnt (count data)
        probabilities (or probabilities (repeat cnt 1))
        sum (v/sum probabilities)
        pmf (->> (map vector data probabilities)
                 (reduce (fn [m [v ^double p]]
                           (update m v (fnil m/+ 0.0) (m// p sum))) (sorted-map)))
        probs (vals pmf)
        cumsum (reductions m/+ probs)
        ks (keys pmf)
        mnk (double (first ks))
        step-before (step-interp/step-before cumsum ks)
        step-after (step-interp/step-after ks cumsum)
        mean (delay (v/dot ks probs))]
    (->distribution {:pdf (fn ^double [^double k] (get pmf k 0.0))
                     :icdf (fn ^double [^double x] (step-before x))
                     :cdf (fn ^double [^double x] (if (m/< x mnk) 0.0 (step-after x)))                     
                     :rng rng
                     :dimensions 1
                     :continuous? false
                     :name :real-discrete-distribution
                     :parameters [:data :probabilities :rng]
                     :mean mean
                     :variance (delay (m/- (v/dot (v/sq ks) probs)
                                           (m/sq @mean)))
                     :lower-bound mnk
                     :upper-bound (last ks)})))

(defn categorical-distribution
  [data probabilities rng]
  (let [rng (or rng (JDKRandomGenerator.))
        
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
  "Solves for x, where cdf(x)=p"
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
                     :parameters [:ns :fn :n :omega :rng]
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
                     :parameters [:ns :fn :n :omega :rng]
                     :mean (delay (prot/mean d))
                     :variance (delay (prot/variance d))
                     :lower-bound lower-bound
                     :upper-bound upper-bound})))


;;

(defn reciprocal
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
  [{:keys [^double mu ^double sigma ^double tau rng normal exponential]}]
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
                     :parameters [:mu :sigma :tau :rng]
                     :lower-bound ##-Inf
                     :upper-bound ##Inf})))

(defn beta-binomial
  [^double alpha ^double beta ^long n rng]
  (let [rng (or rng (JDKRandomGenerator.))
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
                     :parameters [:alpha :beta :n :rng]
                     :lower-bound 0
                     :upper-bound n})))

(defn zero-inflated-beta-binomial
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
  [{:keys [^double mu ^double sigma rng poisson]}]
  (let [rng (or rng (JDKRandomGenerator.))
        dist (poisson {:p mu :rng rng})
        p0 (double (prot/pdf dist 0))
        sigma- (m/- 1.0 sigma)
        pdf0 (m/+ sigma (m/* sigma- p0))
        mean (m/* sigma- mu)]
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
                     :variance (m/* mean (m/inc (m/* mu sigma)))
                     :dimensions 1
                     :continuous? false
                     :name :zero-inflated-poisson
                     :parameters [:mu :sigma :rng]
                     :lower-bound 0
                     :upper-bound Integer/MAX_VALUE})))

(defn zero-adjusted-poisson
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
     :variance variance-v}))

(defn generalized-hyperbolic
  [opts]
  (->distribution (assoc (generalized-hyperbolic-core opts)
                         :dimensions 1
                         :continuous? true
                         :name :generalized-hyperbolic
                         :parameters [:mu :delta :alpha :beta :lambda :rng]
                         :lower-bound ##-Inf
                         :upper-bound ##Inf)))

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
                         :dimensions 1
                         :continuous? true
                         :name :normal-inverse-gaussian
                         :parameters [:alpha :beta :mu :delta :rng]
                         :lower-bound ##-Inf
                         :upper-bound ##Inf)))

(defn half-logistic
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
      :variance (delay (if (m/> df 2.0)
                          (m/- (m// (m/* df (m/inc (m/* delta delta))) (m/- df 2.0))
                               (m/* delta delta half-df (m/sq (double @gamma-ratio))))
                          ##NaN))
      :dimensions 1
      :continuous? true
      :name :t-noncentral
      :parameters [:df :ncp :rng]
      :lower-bound ##-Inf
      :upper-bound ##Inf})))

(defn beta-noncentral
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

;; ---- truncated and mixture

(defn truncated
  [distr left right]
  (let [nname (keyword (str "truncated-" (name (prot/distribution-id distr))))
        lower-bound (double (or left (prot/lower-bound distr)))
        upper-bound (double (or right (prot/upper-bound distr)))
        left-cdf (double (prot/cdf distr lower-bound))
        right-cdf (double (prot/cdf distr upper-bound))
        diff-cdf (m/- right-cdf left-cdf)
        ldiff-cdf (m/log diff-cdf)
        mean (delay (m// (double (quad/gk-quadrature (fn ^double [^double x] (m/* x (double (prot/pdf distr x))))
                                                     lower-bound upper-bound))
                         diff-cdf))
        variance (delay (m/- (m// (double (quad/gk-quadrature (fn ^double [^double x] (m/* x x (double (prot/pdf distr x))))
                                                              lower-bound upper-bound))
                                  diff-cdf)
                             (m/sq @mean)))]
    (->distribution {:lpdf (fn ^double [^double x] (if (m/<= lower-bound x upper-bound)
                                                    (m/- (double (prot/lpdf distr x)) ldiff-cdf)
                                                    ##-Inf))
                     :cdf (fn ^double [^double x] (m/constrain (m// (m/- (double (prot/cdf distr x)) left-cdf) diff-cdf) 0.0 1.0))
                     :icdf (fn ^double [^double p] (prot/icdf distr (m/+ left-cdf (m/* p diff-cdf))))
                     :sampler (fn ^double [] (let [v (double (prot/sample distr))]
                                              (if (m/<= lower-bound v upper-bound) v (recur))))
                     :rng distr
                     :dimensions 1
                     :continuous? (prot/continuous? distr)
                     :name nname
                     :parameters (prot/distribution-parameters distr)
                     :mean mean
                     :variance variance
                     :lower-bound lower-bound
                     :upper-bound upper-bound})))

(defn mixture
  [distrs weights rng]
  (let [rng (or rng (JDKRandomGenerator.))
        probs (v/normalize-L1 weights)
        enum (categorical-distribution distrs probs rng)
        cdf (fn ^double [^double x] (v/dot (map (fn ^double [d] (prot/cdf d x)) distrs) probs))
        mean (delay (v/dot (map prot/mean distrs) probs))
        variance (delay (m/- (v/dot (v/add (v/sq (map prot/mean distrs))
                                           (map prot/variance distrs)) probs) (m/sq @mean)))]
    (->distribution {:pdf (fn ^double [^double x] (v/dot (map (fn ^double [d] (prot/pdf d x)) distrs) probs))
                     :cdf cdf
                     :icdf (fn ^double [^double p]
                             (let [icdfs (map (fn ^double [d] (prot/icdf d p)) distrs)
                                   mn (v/mn icdfs)
                                   mx (v/mx icdfs)
                                   target-fn (fn ^double [^double v] (m/- (double (cdf v)) p))]
                               (solver/find-root target-fn mn mx)))
                     :sampler (fn ^double [] (prot/sample (prot/sample enum)))
                     :seeder (fn [^long seed]
                               (prot/set-seed! rng seed)
                               (doseq [d distrs]
                                 (prot/set-seed! d seed)))
                     :dimensions 1
                     :rng rng
                     :continuous? (some identity (map prot/continuous? distrs))
                     :name :mixture
                     :parameters nil
                     :mean mean
                     :variance variance
                     :lower-bound (v/mn (map prot/lower-bound distrs))
                     :upper-bound (v/mx (map prot/upper-bound distrs))})))
