(ns fastmath.random.distributions
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.special :as special]
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
  [{:keys [pdf lpdf cdf icdf rng ^long dimensions continuous? name parameters mean variance lower-bound upper-bound sampler]
    :or {dimensions 1}}]
  (let [rng (or rng (JDKRandomGenerator.))
        pdf (or pdf (fn ^double [^double v] (m/exp (double (lpdf v)))))
        lpdf (or lpdf (fn ^double [^double v] (m/log (double (pdf v)))))
        sampler (or sampler (fn ^double [] (icdf (prot/drandom rng))))]
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
      (set-seed! [d seed] (prot/set-seed! rng seed) d)
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

  Possible interpolation methods: `:linear` (default), `:spline`, `:monotone` or any function from `fastmath.interpolation`"
  ([pdf-func mn mx steps]
   (integrate-pdf pdf-func {:mn mn :mx mx :steps steps}))
  ([pdf-func {:keys [^double mn ^double mx ^long steps interpolator]
              :or {mn 0.0 mx 1.0 steps 1000 interpolator :linear}
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
                  :linear linear-interp/linear
                  :cubic cubic-interp/cubic
                  :monotone monotone-interp/monotone
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
                     :upper-bound (m/+ mn step)})))

;;

(defn discrete-binary-search
  ([cdf-fn ^double p [mid step]] (discrete-binary-search cdf-fn step p [0 mid]))
  ([cdf-fn ^long step ^double p [^long mn ^long mx]]
   (cond
     (m/> p (double (cdf-fn mx))) (recur cdf-fn (m/* 2 step) p [mx (m/+ mx step)])
     (m/one? (m/- mx mn)) (if (m/>= (double (cdf-fn mn)) p) mn mx)
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
                     :sampler (fn [] (unchecked-long (icdf (prot/drandom rng))))
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
      (loop [interval (m/+ init step) 
             j (long 2)]
        (if (m/< (double (cdf interval)) p)
          (recur (m/+ init (m/* j step)) (m/inc j))
          (solver/find-root h1 init interval)))
      (loop [interval (m/- init step)
             j 2]
        (if (m/> (double (cdf interval)) p)
          (recur (m/- init (m/* j step)) (m/inc j))
          (solver/find-root h1 interval init))))))

;;

(defn- kolmogorov-pdf
  [^double [^double x]]
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
                     :continuous? true
                     :name nname
                     :parameters (prot/distribution-parameters distr)
                     :mean mean
                     :variance variance
                     :lower-bound lower-bound
                     :upper-bound upper-bound})))
