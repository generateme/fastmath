^:kindly/hide-code
(ns random
  {:clj-kondo/config '{:config-in-call {fastmath.dev.clay/symbol-info-table {:ignore [:unresolved-symbol]}}}}
  (:require [fastmath.core :as m]
            [fastmath.random :as r]
            [fastmath.kernel :as k]
            [scicloj.kindly.v4.kind :as kind]
            [fastmath.dev.ggplot :as gg]
            [fastmath.dev.clay :as utls]
            [fastmath.dev.codox :as codox]))

;; # Random {.unnumbered}

;; Collection of functions which deal with randomness. There are four groups:

;; * random number generation
;; * probability distributions
;; * sequences generators
;; * noise (value, gradient, simplex)

^:kindly/hide-code
(kind/code "(require '[fastmath.random :as r])")

;; Set initial seed
(r/set-seed! 31337)

;; ## Common functions

;; List of the functions which work with PRNG and distribution objects:

;; ### Random number generation

;; * First argument should be PRNG or distribution object
;; * With no additional arguments, functions return:
;;     * full range for integers and longs
;;     * $[0,1)$ for floats and doubles
;; * Only for PRNGs
;;     * With one additional argument, number from $[0,max)$ range is returned
;;     * With two additional arguments, number form $[min,max)$ range is returned

(utls/symbol-info-table
 [[irandom "random integer, returns long"]
  [lrandom "random long"]
  [frandom "random float, returns boxed float"]
  [drandom "random double"]])

;; ### Sampling

;; To generate infinite lazy sequence of random values, call `->seq` on given PRNG or distribution.
;; Second (optional) argument can limit number of returned numbers.
;; Third (optional) argument controls sampling method (for given `n` number of samples), there are the following options:

;; * `:antithetic` - random values in pairs `[r1,1-r1,r2,1-r2,...]`
;; * `:uniform` - spacings between numbers follow uniform distribution
;; * `:systematic` - the same spacing with random starting point
;; * `:stratified` - divide $[0,1]$ into `n` intervals and get random value from each subinterval

(utls/symbol-info-table
 [[->seq "sequence of random doubles, sampling"]])

;; ### Seed

(utls/symbol-info-table
 [[set-seed! "sets seed, possibly mutating, returns PRNG or distribution object"]
  [set-seed "creates new instance with given seed"]])

;; * `set-seed!` mutates object if it's supported by underlying Java implementation, can also return newly created object.
;; * `set-seed` is implemented only for PRNGs currently

;; ## PRNG

;; Random number generation is based on PRNG (Pseudorandom Numeber Generator) objects which are responsible for keeping the state. All PRNGs are based on Apache Commons Math 3.6.1 algorithms.

(kind/table
 {:column-names ["Algorithm" "Description"]
  :row-vectors [[(kind/code ":mersenne") "Mersenne Twister"]
                [(kind/code ":isaac") "ISAAC, cryptographically secure PRNG"]
                [(kind/code ":well512a") (kind/md "WELL $2^{512}-1$ period, variant a")]
                [(kind/code ":well1024a") (kind/md "WELL $2^{1024}-1$ period, variant a")]
                [(kind/code ":well19937a") (kind/md "WELL $2^{19937}-1$ period, variant a")]
                [(kind/code ":well19937c") (kind/md "WELL $2^{19937}-1$ period, variant c")]
                [(kind/code ":well44497a") (kind/md "WELL $2^{44497}-1$ period, variant a")]
                [(kind/code ":well44497b") (kind/md "WELL $2^{44497}-1$ period, variant b")]
                [(kind/code ":jdk") "java.util.Random instance, thread safe"]]})

;; To create PRNG, call `rng` with algorithm name (as a keyword) and optional seed parameter.

(utls/symbol-info-table
 [[rng "Create PRNG with optional seed"]
  [synced-rng "Wrap PRNG to ensure thread safety"]])

(utls/examples-note
 (r/rng :isaac)
 (r/rng :isaac 5336)
 (r/synced-rng (r/rng :isaac)))

;; ### Functions

;; Two additional functions are supported by PRNGs

(utls/symbol-info-table
 [[grandom "random gaussian"]
  [brandom "random boolean"]])

;; * `grandom`
;;    * 1-arity - returns random number from Normal (Gaussian) distribution, N(0,1)
;;    * 2-arity - N(0,stddev)
;;    * 3-arity - N(mean, stddev)
;; * `brandom`
;;    * 1-arity - returns true/false with probability=0.5
;;    * 2-arity - returns true with given probability

;; All examples will use Mersenne Twister PRNG with random seed

(def my-prng (r/rng :mersenne))

(utls/examples-note
 (r/irandom my-prng)
 (r/irandom my-prng 5)
 (r/irandom my-prng 10 20)
 (r/lrandom my-prng)
 (r/lrandom my-prng 5)
 (r/lrandom my-prng 10 20)
 (r/frandom my-prng)
 (r/frandom my-prng 2.5)
 (r/frandom my-prng 10.1 10.11)
 (r/drandom my-prng)
 (r/drandom my-prng 2.5)
 (r/drandom my-prng 10.1 10.11)
 (r/grandom my-prng)
 (r/grandom my-prng 20.0)
 (r/grandom my-prng 10.0 20.0)
 (r/brandom my-prng)
 (r/brandom my-prng 0.01)
 (r/brandom my-prng 0.99))

;; ### Sampling

(utls/examples-note
 (take 2 (r/->seq my-prng))
 (r/->seq my-prng 2)
 (r/->seq my-prng 5 :antithetic)
 (r/->seq my-prng 5 :uniform)
 (r/->seq my-prng 5 :systematic)
 (r/->seq my-prng 5 :stratified))

;; ### Seed 

;; Let's define two copies of the same PRNG with the same seed. Second one is obtained by setting a seed new seed.

(def isaac-prng (r/rng :isaac 9988))

(def isaac-prng2 ;; new instance
  (r/set-seed isaac-prng 12345))

(utls/examples-note
 (r/->seq isaac-prng 3)
 (r/->seq isaac-prng2 3)
 (r/->seq isaac-prng 3)
 (r/->seq isaac-prng2 3))

;; Let's now reseed both PRNGs

(r/set-seed! isaac-prng 9988)

(r/set-seed! isaac-prng2 9988)

(utls/examples-note
 (r/->seq isaac-prng 3)
 (r/->seq isaac-prng2 3))

;; ### Default PRNG

;; There is defined one global variable `default-rng` which is synchonized `:jvm` PRNG. Following set of functions work on this particular PRNG. The are the same as `xrandom`, ie `(irand)` is the same as `(irandom default-rng)`.

(utls/symbol-info-table
 [[irand "random integer, as long"]
  [lrand "random long"]
  [frand "random float as boxed Float"]
  [drand "random double"]
  [grand "random gaussian"]
  [brand "random boolean"]
  [->seq "returns infinite lazy sequence"]
  [set-seed "seeds default-rng, returns new instance"]
  [set-seed! "seeds default-rng and Smile's RNG"]])

(utls/examples-note
 (r/irand)
 (r/irand 5)
 (r/irand 10 20)
 (r/lrand)
 (r/lrand 5)
 (r/lrand 10 20)
 (r/frand)
 (r/frand 2.5)
 (r/frand 10.1 10.11)
 (r/drand)
 (r/drand 2.5)
 (r/drand 10.1 10.11)
 (r/grand)
 (r/grand 20.0)
 (r/grand 10.0 20.0)
 (r/brand)
 (r/brand 0.01)
 (r/brand 0.99)
 (take 3 (r/->seq))
 r/default-rng
 (r/set-seed)
 (r/set-seed 9988)
 (r/set-seed!)
 (r/set-seed! 9988))

;; Additionally there are some helpers:

(utls/symbol-info-table
 [[randval "A macro, returns value with given probability (default true/false with prob=0.5)"]
  [flip "Returns 1 with given probability (or 0)"]
  [flipb "Returns true with given probability"]
  [roll-a-dice "Returns a result of rolling a n-sides dice(s)"]])

(utls/examples-note
 (r/randval)
 (r/randval 0.99)
 (r/randval 0.01)
 (r/randval :value1 :value2)
 (r/randval 0.99 :highly-probable-value :less-probably-value)
 (r/randval 0.01 :less-probably-value :highly-probable-value)
 (repeatedly 10 r/flip)
 (repeatedly 10 #(r/flip 0.8))
 (repeatedly 10 #(r/flip 0.2))
 (repeatedly 10 r/flipb)
 (repeatedly 10 #(r/flipb 0.8))
 (repeatedly 10 #(r/flipb 0.2))
 (r/roll-a-dice 6)
 (r/roll-a-dice 100)
 (r/roll-a-dice 10 6))

;; Above helpers can accept custom PRNG

(def isaac3-prng (r/rng :isaac 1234))

(utls/symbol-info-table
 [[randval-rng "A macro, returns value with given probability (default true/false with prob=0.5)"]
  [flip-rng "Returns 1 with given probability (or 0)"]
  [flipb-rng "Returns true with given probability"]
  [roll-a-dice-rng "Returns a result of rolling a n-sides dice(s)"]])

(utls/examples-note
 (r/randval-rng isaac3-prng)
 (r/randval-rng isaac3-prng 0.99)
 (r/randval-rng isaac3-prng 0.01)
 (r/randval-rng isaac3-prng :value1 :value2)
 (r/randval-rng isaac3-prng 0.99 :highly-probable-value :less-probably-value)
 (r/randval-rng isaac3-prng 0.01 :less-probably-value :highly-probable-value)
 (repeatedly 10 #(r/flip-rng isaac3-prng))
 (repeatedly 10 #(r/flip-rng isaac3-prng 0.8))
 (repeatedly 10 #(r/flip-rng isaac3-prng 0.2))
 (repeatedly 10 #(r/flipb-rng isaac3-prng))
 (repeatedly 10 #(r/flipb-rng isaac3-prng 0.8))
 (repeatedly 10 #(r/flipb-rng isaac3-prng 0.2))
 (r/roll-a-dice-rng isaac3-prng 6)
 (r/roll-a-dice-rng isaac3-prng 100)
 (r/roll-a-dice-rng isaac3-prng 10 6))

;; ## Distributions

;; Collection of probability distributions. 

(utls/symbol-info-table
 [[distribution "Distribution creator, a multimethod"]])

(utls/examples-note
 (r/distribution :cauchy)
 (r/distribution :cauchy {:median 2.0 :scale 2.0}))

;; ### Functions

(utls/symbol-info-table
 [[distribution? "checks if given object is a distribution"]
  [pdf "probability density (for continuous) or probability mass (for discrete)"]
  [lpdf "log of pdf"]
  [observe1 "log of pdf, alias of lpdf"]
  [cdf "cumulative density or probability, distribution"]
  [ccdf "complementary cdf, 1-cdf"]
  [icdf "inverse cdf, quantiles"]
  [sample "returns random sample"]
  [dimensions "number of dimensions"]
  [continuous? "is distribution continuous (true) or discrete (false)?"]
  [log-likelihood "sum of lpdfs for given seq of samples"]
  [observe "macro version of log-likelihood"]
  [likelihood "exp of log-likelihood"]
  [mean "distribution mean"]
  [means "distribution means for multivariate distributions"]
  [variance "distribution variance"]
  [covariance "covariance for multivariate distributions"]
  [lower-bound "support lower bound"]
  [upper-bound "support upper bound"]
  [distribution-id "name of distribution"]
  [distribution-parameters "list of parameters"]
  [integrate-pdf "construct cdf and icdf out of pdf function"]])

;; * `distribution-parameters` by default returns only obligatory parameters, when last argument is `true`, returns also optional parameters, for example `:rng`
;; * `drandom`, `frandom`, `lrandom` and `irandom` work only on univariate distributions
;; * `lrandom` and `irandom` use `round-even` to convert double to integer.
;; * `cdf`, `ccdf` and `icdf` are defined only for univariate distributions

;; #### Examples

;; Let's define Beta, Dirichlet and Bernoulli distributions as examples of univariate continuous, multivariate continuous and univariate discrete.

;; ##### Log-logistic, univariate continuous

(def log-logistic (r/distribution :log-logistic {:alpha 3 :beta 7}))

(utls/examples-note
 (r/drandom log-logistic)
 (r/frandom log-logistic)
 (r/lrandom log-logistic)
 (r/irandom log-logistic)
 (r/sample log-logistic)
 (r/distribution? log-logistic)
 (r/distribution? 2.3)
 (r/pdf log-logistic 1)
 (r/lpdf log-logistic 1)
 (r/observe1 log-logistic 1)
 (r/cdf log-logistic 1)
 (r/ccdf log-logistic 1)
 (r/icdf log-logistic 0.002907)
 (r/dimensions log-logistic)
 (r/continuous? log-logistic)
 (r/log-likelihood log-logistic (range 1 10))
 (r/observe log-logistic (range 1 10))
 (r/likelihood log-logistic (range 1 10))
 (r/mean log-logistic)
 (r/variance log-logistic)
 (r/lower-bound log-logistic)
 (r/upper-bound log-logistic)
 (r/distribution-id log-logistic)
 (r/distribution-parameters log-logistic)
 (r/distribution-parameters log-logistic true))

;; ##### Dirichlet, multivariate continuous

(def dirichlet (r/distribution :dirichlet {:alpha [1 2 3]}))

(utls/examples-note
 (r/sample dirichlet)
 (reduce + (r/sample dirichlet))
 (r/distribution? dirichlet)
 (r/pdf dirichlet [0.5 0.1 0.4])
 (r/lpdf dirichlet [0.5 0.1 0.4])
 (r/observe1 dirichlet [0.5 0.1 0.4])
 (r/dimensions dirichlet)
 (r/continuous? dirichlet)
 (r/log-likelihood dirichlet (repeatedly 10 #(r/sample dirichlet)))
 (r/observe dirichlet (repeatedly 10 #(r/sample dirichlet)))
 (r/likelihood dirichlet (repeatedly 10 #(r/sample dirichlet)))
 (r/means dirichlet)
 (r/covariance dirichlet)
 (r/distribution-id dirichlet)
 (r/distribution-parameters dirichlet)
 (r/distribution-parameters dirichlet true))

;; * multivariate distributions don't implement some of the functions, like `drandom` or `cdf`, `icdf`, etc.

;; ##### Poisson, univariate discrete

(def poisson (r/distribution :poisson {:p 10}))

(utls/examples-note
 (r/drandom poisson)
 (r/frandom poisson)
 (r/lrandom poisson)
 (r/irandom poisson)
 (r/sample poisson)
 (r/distribution? poisson)
 (r/pdf poisson 5)
 (r/lpdf poisson 5)
 (r/observe1 poisson 5)
 (r/cdf poisson 5)
 (r/ccdf poisson 5)
 (r/icdf poisson 0.999)
 (r/dimensions poisson)
 (r/continuous? poisson)
 (r/log-likelihood poisson (range 1 10))
 (r/observe poisson (range 1 10))
 (r/likelihood poisson (range 1 10))
 (r/mean poisson)
 (r/variance poisson)
 (r/lower-bound poisson)
 (r/upper-bound poisson)
 (r/distribution-id poisson)
 (r/distribution-parameters poisson)
 (r/distribution-parameters poisson true))

;; #### PDF integration

;; `integrate-pdf` returns a pair of CDF and iCDF functions using Romberg integration and interpolation. Given interval is divided into `steps` number of subinterval. Each subinteval is integrated and added to cumulative sum. All points a later interpolated to build CDF and iCDF.

;; Parameters list:

(utls/symbol-info-table
 [[pdf-func "PDF, univariate double->double function"]
  [:mn "lower bound, value of PDF(mn) should be 0.0"]
  [:mx "upper bound"]
  [:steps "number of subintervals"]
  [:min-iterations "minimum number of Romberg integrator iterations (default: 3, minimum 2)"]
  [:interpolator "one of: `:linear` (default), `:spline`, `:monotone` or any `fastmath.interpolator` method"]])

;; Let's compare cdf and icdf to integrated pdf of beta distribution

(def beta (r/distribution :beta))

(def integrated (r/integrate-pdf (partial r/pdf beta)
                                 {:mn -0.001 :mx 1.001 :steps 10000
                                  :interpolator :spline}))

(def beta-cdf (first integrated))

(def beta-icdf (second integrated))

(utls/examples-note
 (r/cdf beta 0.5)
 (beta-cdf 0.5)
 (r/icdf beta 0.5)
 (beta-icdf 0.5))

;; ### Sampling

(def weibull (r/distribution :weibull))

(utls/examples-note
 (r/sample weibull)
 (r/->seq weibull 3)
 (r/->seq weibull 5 :antithetic)
 (r/->seq weibull 5 :uniform)
 (r/->seq weibull 5 :systematic)
 (r/->seq weibull 5 :stratified))

;; ### PRNG and seed

;; All distributions accept `:rng` optional parameter which is used internally for random number generation. By default, creator constructs custom PRNG instance.

(def cauchy-mersenne-2288 (r/distribution :cauchy {:rng (r/rng :mersenne 2288)}))

(utls/examples-note
 (vec (r/->seq cauchy-mersenne-2288 3))
 (vec (r/->seq cauchy-mersenne-2288 3))
 (r/set-seed! cauchy-mersenne-2288 2288)
 (vec (r/->seq cauchy-mersenne-2288 3)))

;; ### Default Normal

(utls/symbol-info-table
 [[default-normal "public var which is a normal distribution N(0,1), thread-safe"]])

(utls/examples-note
 (r/->seq r/default-normal 3)
 (r/pdf r/default-normal 0.0)
 (r/cdf r/default-normal 0.0)
 (r/icdf r/default-normal 0.5)
 (r/mean r/default-normal)
 (r/variance r/default-normal))

^:kindly/hide-code
(defn build-distribution-list
  [multivariate? continuous?]
  (for [nm (->> r/distribution methods keys sort (remove #{:truncated :mixture :continuous-distribution
                                                           :empirical :enumerated-real :kde
                                                           :real-discrete-distribution
                                                           :integer-discrete-distribution
                                                           :categorical-distribution :enumerated-int}))
        :let [di (r/distribution nm)]
        :when (and (= multivariate? (not (m/one? (r/dimensions di))))
                   (= continuous? (r/continuous? di)))]
    [nm (r/distribution-parameters di true)]))

;; ### Univariate, cont.

(kind/table
 {:column-names ["name" "parameters"]
  :row-vectors (mapv (fn [[k v]]
                       [(kind/code (str k)) v]) (build-distribution-list false true))})

;; #### Anderson-Darling

;; Distribution of Anderson-Darling statistic $A^2$ on $n$ independent uniforms $U[0,1]$.

;; * Name: `:anderson-darling`
;; * Default parameters:
;;    * `:n`: $1$
;; * [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1AndersonDarlingDist.html)

(kind/table
 [[{:n 1} {:n 5}]
  [(gg/dgraph (r/distribution :anderson-darling {:n 1}) {:pdf [0 3]})
   (gg/dgraph (r/distribution :anderson-darling {:n 5}) {:pdf [0 3]})]])

;; #### Beta

;; * Name: `:beta`
;; * Default parameters:
;;    * `:alpha`: $2.0$
;;    * `:beta`: $5.0$
;; * [wiki](https://en.wikipedia.org/wiki/Beta_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/BetaDistribution.html)

(kind/table
 [[{:alpha 1 :beta 1} {:alpha 0.5 :beta 0.5}]
  [(gg/dgraph (r/distribution :beta {:alpha 1 :beta 1}) {:pdf [0.01 0.99]})
   (gg/dgraph (r/distribution :beta {:alpha 0.5 :beta 0.5}) {:pdf [0.01 0.99]})]
  [{:alpha 3 :beta 3} {:alpha 5 :beta 2}]
  [(gg/dgraph (r/distribution :beta {:alpha 3 :beta 3}) {:pdf [0 1]})
   (gg/dgraph (r/distribution :beta {:alpha 5 :beta 2}) {:pdf [0 1]})]])

;; #### Cauchy

;; * Name `:cauchy`
;; * Default parameters:
;;    * `:median`, location: $0.0$
;;    * `:scale`: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Cauchy_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/CauchyDistribution.html)

(kind/table
 [[{:median 0 :scale 1} {:median 1 :scale 0.5}]
  [(gg/dgraph (r/distribution :cauchy {:median 0 :scale 1}) {:pdf [-4 4] :icdf [0.02 0.95]})
   (gg/dgraph (r/distribution :cauchy {:median 1 :scale 0.5}) {:pdf [-3 4] :icdf [0.02 0.95]})]])

;; #### Chi

;; * Name: `:chi`
;; * Default parameters:
;;    * `:nu`, degrees of freedom: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Chi_distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1ChiDist.html)

(kind/table
 [[{:nu 1} {:nu 3}]
  [(gg/dgraph (r/distribution :chi {:nu 1}) {:pdf [0.01 5]})
   (gg/dgraph (r/distribution :chi {:nu 3}) {:pdf [0 5]})]])

;; #### Chi-squared

;; * Name: `:chi-squered`
;; * Default parameters:
;;    * `:degrees-of-freedom`: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Chi-squared_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/ChiSquaredDistribution.html)

(kind/table
 [[{:degrees-of-freedom 1} {:degrees-of-freedom 3}]
  [(gg/dgraph (r/distribution :chi-squared) {:pdf [0 5]})
   (gg/dgraph (r/distribution :chi-squared {:degrees-of-freedom 3}) {:pdf [0 6]})]])

;; #### Chi-squared noncentral

;; * Name: `:chi-squared-noncentral`
;; * Default parameters:
;;    * `:nu`, degrees-of-freedom: $1.0$
;;    * `:lambda`, noncentrality: $1.0$
;; * [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1ChiSquareNoncentralDist.html)

(kind/table
 [[{:nu 1 :lambda 1} {:nu 3 :lambda 2}]
  [(gg/dgraph (r/distribution :chi-squared-noncentral) {:pdf [0 5]})
   (gg/dgraph (r/distribution :chi-squared-noncentral {:nu 3 :lambda 2}) {:pdf [0 6]})]])

;; #### Cramer-von Mises

;; Distribution of Cramer-von Mises statistic $W^2$ on $n$ independent uniforms $U[0,1]$.

;; * Name: `:cramer-von-mises`
;; * Default parameters
;;    * `:n`: $1$
;; * [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1CramerVonMisesDist.html)

;; Note: PDF is calculated using finite difference method from CDF.

(kind/table
 [[{:n 1} {:n 5}]
  [(gg/dgraph (r/distribution :cramer-von-mises {:n 1}) {:pdf [0 0.5]})
   (gg/dgraph (r/distribution :cramer-von-mises {:n 5}) {:pdf [0 1]})]])

;; #### Erlang

;; * Name: `:erlang`
;; * Default parameters
;;    * `:k`, shape: $1.0$
;;    * `:lambda`, scale: $1.0$ 
;; * [wiki](https://en.wikipedia.org/wiki/Erlang_distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1ErlangDist.html)

(kind/table
 [[{:k 1 :lambda 1} {:k 7 :lambda 2.0}]
  [(gg/dgraph (r/distribution :erlang) {:pdf [0 5]})
   (gg/dgraph (r/distribution :erlang {:k 7 :lambda 2.0}) {:pdf [0 8]})]])

;; #### ex-Gaussian

;; * Name: `:exgaus`
;; * Default parameters
;;    * `:mu`, mean: $5.0$
;;    * `:sigma`, standard deviation: $1.0$
;;    * `:nu`, mean of exponential variable: $1.0$ 
;; * [wiki](https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution), [source](https://search.r-project.org/CRAN/refmans/gamlss.dist/html/exGAUS.html)

(kind/table
 [[{:mu 0 :sigma 1 :nu 1} {:mu -2 :sigma 0.5 :nu 4}]
  [(gg/dgraph (r/distribution :exgaus) {:pdf [-5 5] :icdf [0.001 0.999]})
   (gg/dgraph (r/distribution :exgaus {:mu -2 :sigma 0.5 :nu 4}) {:pdf [-6 10] :icdf [0.001 0.999]})]])

;; #### Exponential

;; * Name: `:exponential`
;; * Default parameters
;;    * `:mean`: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Exponential_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/ExponentialDistribution.html)

(kind/table
 [[{:mean 1} {:mean 3}]
  [(gg/dgraph (r/distribution :exponential) {:pdf [0 5]})
   (gg/dgraph (r/distribution :exponential {:mean 3}) {:pdf [0 10]})]])

;; #### F

;; * Name: `:f`
;; * Default parameters
;;    * `:numerator-degrees-of-freedom`: $1.0$
;;    * `:denominator-degrees-of-freedom`: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/F-distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/FDistribution.html)

(kind/table
 [[{:numerator-degrees-of-freedom 1 :denominator-degrees-of-freedom 1}]
  [(gg/dgraph (r/distribution :f) {:pdf [0 10] :icdf [0.0 0.9]})]
  [{:numerator-degrees-of-freedom 10 :denominator-degrees-of-freedom 15}]
  [(gg/dgraph (r/distribution :f {:numerator-degrees-of-freedom 10 :denominator-degrees-of-freedom 15}) {:pdf [0 10]})]])

;; #### Fatigue life

;; * Name: `:fatigue-life`
;; * Default parameters
;;    * `:mu`, location: $0.0$
;;    * `:beta`, scale: $1.0$
;;    * `:gamma`, shape: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1FatigueLifeDist.html)

(kind/table
 [[{:mu 0 :beta 1 :gamma 1} {:mu -1 :beta 3 :gamma 0.5}]
  [(gg/dgraph (r/distribution :fatigue-life) {:pdf [-0.5 5]})
   (gg/dgraph (r/distribution :fatigue-life {:mu -1 :beta 3 :gamma 0.5}) {:pdf [-2 7]})]])

;; #### Folded Normal

;; * Name: `:folded-normal`
;; * Default parameters
;;    * `:mu`: $0.0$
;;    * `:sigma`: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Folded_normal_distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1FoldedNormalDist.html)

(kind/table
 [[{:mu 0 :sigma 1} {:mu 2 :sigma 1}]
  [(gg/dgraph (r/distribution :folded-normal) {:pdf [-0.5 3]})
   (gg/dgraph (r/distribution :folded-normal {:mu 2 :sigma 1}) {:pdf [-0.5 5]})]])

;; #### Frechet

;; * Name: `:frechet`
;; * Default parameters
;;    * `:delta`, location: $0.0$
;;    * `:alpha`, shape: $1.0$
;;    * `:beta`, scale: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Fr%C3%A9chet_distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1FrechetDist.html)

(kind/table
 [[{:delta 0 :alpha 1 :beta 1} {:delta 1 :alpha 3 :beta 0.5}]
  [(gg/dgraph (r/distribution :frechet) {:pdf [-0.5 5] :icdf [0 0.9]})
   (gg/dgraph (r/distribution :frechet {:delta 1 :alpha 3 :beta 2}) {:pdf [-0.1 7]})]])

;; #### Gamma

;; * Name: `:gamma`
;; * Default parameters
;;    * `:shape`: $2.0$
;;    * `:scale`: $2.0$
;; * [wiki](https://en.wikipedia.org/wiki/Gamma_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/GammaDistribution.html)

(kind/table
 [[{:shape 2 :scale 2} {:shape 5 :scale 0.5}]
  [(gg/dgraph (r/distribution :gamma) {:pdf [-0.5 15]})
   (gg/dgraph (r/distribution :gamma {:shape 5 :scale 0.5}) {:pdf [-0.1 7]})]])

;; #### Gumbel

;; * Name: `:gumbel`
;; * Default parameters
;;    * `:mu`, location: $1.0$
;;    * `:beta`, scale: $2.0$
;; * [wiki](https://en.wikipedia.org/wiki/Gumbel_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/GumbelDistribution.html)

(kind/table
 [[{:mu 1 :beta 2} {:mu 1 :beta 0.5}]
  [(gg/dgraph (r/distribution :gumbel) {:pdf [-5 10]})
   (gg/dgraph (r/distribution :gumbel {:mu 1 :beta 0.5}) {:pdf [-1 5]})]])

;; #### Half Cauchy

;; * Name: `:half-cauchy`
;; * Default parameters
;;    * `:scale`: $1.0$
;; * [info](https://distribution-explorer.github.io/continuous/halfcauchy.html)

(kind/table
 [[{:scale 1} {:scale 2}]
  [(gg/dgraph (r/distribution :half-cauchy) {:pdf [0 8] :icdf [0 0.95]})
   (gg/dgraph (r/distribution :half-cauchy {:scale 2}) {:pdf [0 8] :icdf [0 0.95]})]])

;; #### Half Normal

;; * Name: `:half-normal`
;; * Default parameters
;;    * `:sigma`: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Half-normal_distribution)

(kind/table
 [[{:sigma 1} {:sigma 2}]
  [(gg/dgraph (r/distribution :half-normal) {:pdf [0 5]})
   (gg/dgraph (r/distribution :half-normal {:sigma 2}) {:pdf [0 7]})]])

;; #### Hyperbolic secant

;; * Name: `:hyperbolic-secant`
;; * Default parameters
;;    * `:mu`: $0.0$
;;    * `:sigma`, scale: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1HyperbolicSecantDist.html)

(kind/table
 [[{:mu 0 :sigma 1} {:mu 1 :sigma 2}]
  [(gg/dgraph (r/distribution :hyperbolic-secant) {:pdf [-5 5]})
   (gg/dgraph (r/distribution :hyperbolic-secant {:mu 1 :sigma 2}) {:pdf [-7 7]})]])

;; #### Hypoexponential

;; * Name: `:hypoexponential`
;; * Default parameters
;;    * `:lambdas`, list of rates: `[1.0]`
;; * [wiki](https://en.wikipedia.org/wiki/Hypoexponential_distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1HypoExponentialDist.html)

(kind/table
 [[{:lambdas [1]} {:lambdas [1 2 3 4 1]}]
  [(gg/dgraph (r/distribution :hypoexponential) {:pdf [0 5]})
   (gg/dgraph (r/distribution :hypoexponential {:lambdas [1 2 3 4 1]}) {:pdf [0 10]})]])

;; #### Hypoexponential equal

;; Hypoexponential distribution, where $\lambda_i=(n+1-i)h,\text{ for } i=1\dots k$

;; * Name: `:hypoexponential-equal`
;; * Default parameters
;;    * `:k`, number of rates: $1$
;;    * `:h`, difference between rates: $1$
;;    * `:n` $=\frac{\lambda_1}{h}$: $1$
;; * [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1HypoExponentialDistEqual.html)

(kind/table
 [[{:n 1 :k 1 :h 1} {:n 5 :h 0.5 :k 6}]
  [(gg/dgraph (r/distribution :hypoexponential-equal) {:pdf [0 5]})
   (gg/dgraph (r/distribution :hypoexponential-equal {:n 8 :h 0.5 :k 6}) {:pdf [0 10]})]])

;; #### Inverse Gamma

;; * Name: `:inverse-gamma`
;; * Default parameters
;;    * `:alpha`, shape: $2.0$
;;    * `:beta`, scale: $1.0$
;; *[wiki](https://en.wikipedia.org/wiki/Inverse-gamma_distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1InverseGammaDist.html)

(kind/table
 [[{:alpha 2 :beta 1} {:alpha 1.5 :beta 2}]
  [(gg/dgraph (r/distribution :inverse-gamma) {:pdf [0 3]})
   (gg/dgraph (r/distribution :inverse-gamma {:alpha 1.5 :beta 2}) {:pdf [0 7] :icdf [0.01 0.95]})]])

;; #### Inverse Gaussian

;; * Name: `:inverse-gaussian`
;; * Default parameters
;;    * `:mu`, location: $1.0$
;;    * `:lambda`, scale: $1.0$
;; *[wiki](https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1InverseGaussianDist.html)

(kind/table
 [[{:mu 1 :lambda 1} {:mu 2 :lambda 0.5}]
  [(gg/dgraph (r/distribution :inverse-gaussian) {:pdf [0 3]})
   (gg/dgraph (r/distribution :inverse-gaussian {:mu 2 :lambda 0.5}) {:pdf [0 3]})]])

;; #### Johnson Sb

;; * Name: `:johnson-sb`
;; * Default parameters
;;    * `:gamma`, shape: $0.0$
;;    * `:delta`, shape: $1.0$
;;    * `:xi`, location: $0.0$
;;    * `:lambda`, scale: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Johnson%27s_SU-distribution#Johnson's_SB-distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1JohnsonSBDist.html)

(kind/table
 [[{:gamma 0 :delta 1 :xi 0 :lambda 1} {:gamma 0.5 :delta 2 :xi 0 :lambda 2}]
  [(gg/dgraph (r/distribution :johnson-sb) {:pdf [0 1]})
   (gg/dgraph (r/distribution :johnson-sb {:gamma 0.5 :delta 2 :xi 0 :lambda 2}) {:pdf [0 3]})]])

;; #### Johnson Sl

;; * Name: `:johnson-sl`
;; * Default parameters
;;    * `:gamma`, shape: $0.0$
;;    * `:delta`, shape: $1.0$
;;    * `:xi`, location: $0.0$
;;    * `:lambda`, scale: $1.0$
;; * [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1JohnsonSLDist.html)

(kind/table
 [[{:gamma 0 :delta 1 :xi 0 :lambda 1} {:gamma 0.5 :delta 2 :xi 0 :lambda 2}]
  [(gg/dgraph (r/distribution :johnson-sl) {:pdf [0 4]})
   (gg/dgraph (r/distribution :johnson-sl {:gamma 0.5 :delta 2 :xi 1 :lambda 2}) {:pdf [0 6]})]])

;; #### Johnson Su

;; * Name: `:johnson-su`
;; * Default parameters
;;    * `:gamma`, shape: $0.0$
;;    * `:delta`, shape: $1.0$
;;    * `:xi`, location: $0.0$
;;    * `:lambda`, scale: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Johnson%27s_SU-distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1JohnsonSUDist.html)

(kind/table
 [[{:gamma 0 :delta 1 :xi 0 :lambda 1} {:gamma 0.5 :delta 2 :xi 0 :lambda 2}]
  [(gg/dgraph (r/distribution :johnson-su) {:pdf [-3 3]})
   (gg/dgraph (r/distribution :johnson-su {:gamma 0.5 :delta 2 :xi 1 :lambda 2}) {:pdf [-3 4]})]])

;; #### Kolmogorov

;; * Name: `:kolmogorov`
;; * [wiki](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test#Kolmogorov_distribution), [info](https://www.math.ucla.edu/~tom/distributions/Kolmogorov.html)

(gg/dgraph (r/distribution :kolmogorov) {:pdf [0 2]})

;; #### Kolmogorov-Smirnov

;; * Name: `:kolmogorov-smirnov`
;; * Default parameters
;;    * `:n`, sample size: 1
;; * [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1KolmogorovSmirnovDist.html)

(kind/table
 [[{:n 1} {:n 10}]
  [(gg/dgraph (r/distribution :kolmogorov-smirnov) {:pdf [0 1]})
   (gg/dgraph (r/distribution :kolmogorov-smirnov {:n 10}) {:pdf [0 1]})]])

;; #### Kolmogorov-Smirnov+

;; * Name: `:kolmogorov-smirnov+`
;; * Default parameters
;;    * `:n`, sample size: 1
;; * [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1KolmogorovSmirnovPlusDist.html)

(kind/table
 [[{:n 1} {:n 10}]
  [(gg/dgraph (r/distribution :kolmogorov-smirnov+) {:pdf [0 1]})
   (gg/dgraph (r/distribution :kolmogorov-smirnov+ {:n 10}) {:pdf [0 1]})]])

;; #### Laplace

;; * Name: `:laplace`
;; * Default parameters
;;    * `:mu`: $0.0$
;;    * `:beta`, scale: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Laplace_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/LaplaceDistribution.html)

(kind/table
 [[{:mu 0 :beta 1} {:mu 1 :beta 2}]
  [(gg/dgraph (r/distribution :laplace) {:pdf [-4 4]})
   (gg/dgraph (r/distribution :laplace {:mu 1 :beta 2}) {:pdf [-5 7]})]])

;; #### Levy

;; * Name: `:levy`
;; * Default parameters
;;    * `:mu`: $0.0$
;;    * `:c`, scale: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/L%C3%A9vy_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/LevyDistribution.html)

(kind/table
 [[{:mu 0 :c 1} {:mu 1 :c 2}]
  [(gg/dgraph (r/distribution :levy) {:pdf [0 7] :icdf [0.0 0.85]})
   (gg/dgraph (r/distribution :levy {:mu 1 :c 2}) {:pdf [0 10] :icdf [0.0 0.85]})]])

;; #### Log Logistic

;; * Name: `:log-logistic`
;; * Default parameters
;;    * `:alpha`, shape: $3.0$
;;    * `:beta`, scale: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Log-logistic_distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1LoglogisticDist.html)

(kind/table
 [[{:alpha 3 :beta 1} {:alpha 5 :beta 2}]
  [(gg/dgraph (r/distribution :log-logistic) {:pdf [0 3]})
   (gg/dgraph (r/distribution :log-logistic {:alpha 5 :beta 2}) {:pdf [0 5]})]])

;; #### Log Normal

;; * Name: `:log-normal`
;; * Default parameters
;;    * `:scale`: $1.0$
;;    * `:shape`: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Log-normal_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/LogNormalDistribution.html)

(kind/table
 [[{:scale 1 :shape 1} {:scale 2 :shape 2}]
  [(gg/dgraph (r/distribution :log-normal) {:pdf [0 8]})
   (gg/dgraph (r/distribution :log-normal {:scale 0.5 :shape 2}) {:pdf [0 8] :icdf [0.0 0.9]})]])

;; #### Logistic

;; * Name: `:logistic`
;; * Default parameters
;;    * `:mu`, location: $0.0$
;;    * `:s`, scale: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Logistic_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/LogisticDistribution.html)

(kind/table
 [[{:mu 0 :scale 1} {:mu 1 :scale 0.5}]
  [(gg/dgraph (r/distribution :logistic) {:pdf [-6 6]})
   (gg/dgraph (r/distribution :logistic {:mu 1 :scale 0.5}) {:pdf [-6 8]})]])

;; #### Nakagami

;; * Name: `:nakagami`
;; * Default parameters
;;    * `:mu`, shape: $1.0$
;;    * `:omega`, spread: $1.0$ 
;; * [wiki](https://en.wikipedia.org/wiki/Nakagami_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/NakagamiDistribution.html)

(kind/table
 [[{:mu 1 :omega 1} {:mu 0.5 :omega 0.5}]
  [(gg/dgraph (r/distribution :nakagami) {:pdf [0 3]})
   (gg/dgraph (r/distribution :nakagami {:mu 0.5 :omega 0.5}) {:pdf [0 3]})]])

;; #### Normal

;; * Name: `:normal`
;; * Default parameters
;;    * `:mu`, mean: $0.0$
;;    * `:sd`, standard deviation: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Normal_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/NormalDistribution.html)

(kind/table
 [[{:mu 0 :sd 1} {:mu 0.5 :sd 0.5}]
  [(gg/dgraph (r/distribution :normal) {:pdf [-3 3]})
   (gg/dgraph (r/distribution :normal {:mu 0.5 :sd 0.5}) {:pdf [-3 3]})]])

;; #### Normal-Inverse Gaussian

;; * Name: `:normal-inverse-gaussian`
;; * Default parameters
;;    * `:alpha`, tail heavyness: $1.0$
;;    * `:beta`, assymetry: $0.0$
;;    * `:mu`, location: $0.0$
;;    * `:delta`, scale: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Normal-inverse_Gaussian_distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1NormalInverseGaussianDist.html)

;; Only PDF is supported, you may call `integrate-pdf` to get CDF and iCDF pair.

(kind/table
 [[{:alpha 1 :beta 0 :mu 0 :delta 1} {:alpha 2 :beta 1 :mu 0 :delta 0.5}]
  [(gg/fgraph (partial r/pdf (r/distribution :normal-inverse-gaussian)) [-3 3])
   (gg/fgraph (partial r/pdf (r/distribution :normal-inverse-gaussian {:alpha 5 :beta 4 :mu 0 :delta 0.5})) [-2 3])]])

(let [[cdf icdf] (r/integrate-pdf
                  (partial r/pdf (r/distribution :normal-inverse-gaussian))
                  {:mn -800.0 :mx 800.0 :steps 5000
                   :interpolator :monotone})]
  [(cdf 0.0) (icdf 0.5)])

(let [[cdf icdf] (r/integrate-pdf
                  (partial r/pdf (r/distribution :normal-inverse-gaussian))
                  {:mn -800.0 :mx 800.0 :steps 5000
                   :interpolator :monotone})]
  (kind/table
   [["CDF" "iCDF"]
    [(gg/fgraph cdf [-3 3])
     (gg/fgraph icdf [0.01 0.99])]]))

;; #### Pareto

;; * Name: `:pareto`
;; * Default parameters
;;    * `:shape`: $1.0$
;;    * `:scale`: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Pareto_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/ParetoDistribution.html)

(kind/table
 [[{:scale 1 :shape 1} {:scale 2 :shape 2}]
  [(gg/dgraph (r/distribution :pareto) {:pdf [0 8] :icdf [0.01 0.95]})
   (gg/dgraph (r/distribution :pareto {:scale 2 :shape 2}) {:pdf [0 8]})]])

;; #### Pearson VI

;; * Name: `:pearson-6`
;; * Default parameters:
;;    * `:alpha1`: $1.0$
;;    * `:alpha2`: $1.0$
;;    * `:beta`: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Pearson_distribution#The_Pearson_type_VI_distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1Pearson6Dist.html)

(kind/table
 [[{:alpha1 1 :alpha2 1 :beta 1} {:scale 2 :shape 2}]
  [(gg/dgraph (r/distribution :pearson-6) {:pdf [0 8] :icdf [0.01 0.95]})
   (gg/dgraph (r/distribution :pearson-6 {:alpha1 2 :alpha2 2 :beta 2}) {:pdf [0 8] :icdf [0.0 0.98]})]])

;; #### Power

;; * Name: `:power`
;; * Default parameters
;;    * `:a`: $0.0$
;;    * `:b`: $1.0$
;;    * `:c`: $2.0$
;; * [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1PowerDist.html)

(kind/table
 [[{:a 0 :b 1 :c 2} {:a 1 :b 2 :c 1.25}]
  [(gg/dgraph (r/distribution :power) {:pdf [0 2]})
   (gg/dgraph (r/distribution :power {:a 1 :b 2 :c 1.25}) {:pdf [0 2]})]])

;; #### Rayleigh

;; * Name: `:rayleigh`
;; * Default parameters
;;    * `:a`, location: $0.0$
;;    * `:beta`, scale: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Rayleigh_distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1RayleighDist.html)

(kind/table
 [[{:a 0 :b 1 :c 2} {:a 1 :b 2 :c 1.25}]
  [(gg/dgraph (r/distribution :rayleigh) {:pdf [0 4]})
   (gg/dgraph (r/distribution :rayleigh {:a 1 :beta 0.5}) {:pdf [0 3]})]])

;; #### Reciprocal Sqrt

;; $\operatorname{PDF}(x)=\frac{1}{\sqrt{x}}, x\in(a,(\frac{1}{2}(1+2\sqrt{a}))^2)$

;; * Name: `:reciprocal-sqrt`
;; * Default parameters
;;    * `:a`, location, lower limit: $0.5$

(kind/table
 [[{:a 0.5} {:a 2}]
  [(gg/dgraph (r/distribution :reciprocal-sqrt) {:pdf [0 2]})
   (gg/dgraph (r/distribution :reciprocal-sqrt {:a 2}) {:pdf [0 4]})]])

;; #### Student's t

;; * Name: `:t`
;; * Default parameters
;;    * `:degrees-of-freedom`: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Student%27s_t-distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/TDistribution.html)

(kind/table
 [[{:degrees-of-freedom 1} {:degrees-of-freedom 50}]
  [(gg/dgraph (r/distribution :t) {:pdf [-5 5] :icdf [0.04 0.96]})
   (gg/dgraph (r/distribution :t {:degrees-of-freedom 50}) {:pdf [-5 5]})]])

;; #### Triangular

;; * Name: `:triangular`
;; * Default parameters
;;    * `:a`, lower limit: $-1.0$
;;    * `:b`, mode: $0.0$
;;    * `:c`, upper limit: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Triangular_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/TriangularDistribution.html)

(kind/table
 [[{:a -1 :b 1 :c 0} {:a -0.5 :b 1 :c 0.5}]
  [(gg/dgraph (r/distribution :triangular) {:pdf [-1.5 1.5]})
   (gg/dgraph (r/distribution :triangular {:a -0.5 :c 0.5}) {:pdf [-1.5 1.5]})]])

;; #### Uniform

;; * Name: `:uniform-real`
;; * Default parameters
;;    * `:lower`, lower limit: $0.0$
;;    * `:upper`, upper limit: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Continuous_uniform_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/UniformRealDistribution.html)

(kind/table
 [[{:lower 0 :upper 1} {:lower -1 :upper 0.5}]
  [(gg/dgraph (r/distribution :uniform-real) {:pdf [-1.1 1.1]})
   (gg/dgraph (r/distribution :uniform-real {:lower -1 :upper 0.5}) {:pdf [-1.1 1.1]})]])

;; #### Watson G

;; * Name: `:watson-g`
;; * Default parameters
;;    * `:n`: $2$
;; * [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1WatsonGDist.html)

(kind/table
 [[{:n 2} {:n 10}]
  [(gg/dgraph (r/distribution :watson-g) {:pdf [0 1.5]})
   (gg/dgraph (r/distribution :watson-g {:n 10}) {:pdf [0 1.5]})]])

;; #### Watson U

;; * Name: `:watson-u`
;; * Default parameters
;;    * `:n`: $2$
;; * [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1WatsonUDist.html)

(kind/table
 [[{:n 2} {:n 10}]
  [(gg/dgraph (r/distribution :watson-u) {:pdf [0 0.5]})
   (gg/dgraph (r/distribution :watson-u {:n 10}) {:pdf [0 0.5]})]])

;; #### Weibull

;; * Name: `:weibull`
;; * Default parameters
;;    * `:alpha`, shape: $2.0$
;;    * `:beta`, scale: $1.0$
;; * [wiki](https://en.wikipedia.org/wiki/Weibull_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/WeibullDistribution.html)

(kind/table
 [[{:alpha 2 :beta 1} {:alpha 1.2 :beta 0.8}]
  [(gg/dgraph (r/distribution :weibull) {:pdf [0 3]})
   (gg/dgraph (r/distribution :weibull {:alpha 1.2 :beta 0.8}) {:pdf [0 3]})]])

;; #### Zero adjusted Gamma (zaga)

;; * Name: `:zaga`
;; * Default parameters:
;;    * `:mu`, location: $0.0$
;;    * `:sigma`, scale: $1.0$
;;    * `:nu`, density at 0.0: $0.1$
;;    * `:lower-tail?` - true
;; * [source](https://search.r-project.org/CRAN/refmans/gamlss.dist/html/ZAGA.html), [book](https://www.gamlss.com/wp-content/uploads/2018/01/DistributionsForModellingLocationScaleandShape.pdf)

(kind/table
 [[{:mu 0 :sigma 1 :nu 0.1} {:mu 1 :sigma 1 :nu 0.5}]
  [(gg/dgraph (r/distribution :zaga) {:pdf [0 3]})
   (gg/dgraph (r/distribution :zaga {:mu 1 :sigma 0.5 :nu 0.7}) {:pdf [0 3]})]])

;; ### Univariate, discr.

(kind/table
 {:column-names ["name" "parameters"]
  :row-vectors (map (fn [[k v]]
                      [(kind/code (str k)) v])
                    (build-distribution-list false false))})

;; #### Beta Binomial (bb)

;; * Name: `:bb`
;; * Default parameters
;;    * `:mu`, probability: $0.5$
;;    * `:sigma`, dispersion: $1.0$
;;    * `:bd`, binomial denominator: $10$
;; * [wiki](https://en.wikipedia.org/wiki/Beta-binomial_distribution), [source](https://search.r-project.org/CRAN/refmans/gamlss.dist/html/BB.html)

;; Parameters $\mu,\sigma$ in terms of $\alpha, \beta$ (Wikipedia definition)

;; * probability: $\mu=\frac{\alpha}{\alpha+\beta}$
;; * dispersion: $\sigma=\frac{1}{\alpha+\beta}$

(kind/table
 [[{:mu 0.5 :sigma 1 :bd 10} {:mu 0.65 :sigma 0.3 :bd 20}]
  [(gg/dgraphi (r/distribution :bb) {:pdf [0 11]})
   (gg/dgraphi (r/distribution :bb {:mu 0.65 :sigma 0.3 :bd 20}) {:pdf [0 22]})]])

;; #### Bernoulli

;; The same as Binomial with trials=1.

;; * Name: `:bernoulli`
;; * Default parameters
;;    * `:p`, probability, $0.5$ 
;; * [wiki](https://en.wikipedia.org/wiki/Bernoulli_distribution)

(kind/table
 [[{:p 0.5} {:p 0.25}]
  [(gg/dgraphi (r/distribution :bernoulli) {:pdf [0 2]})
   (gg/dgraphi (r/distribution :bernoulli {:p 0.25}) {:pdf [0 2]})]])

;; #### Binomial

;; * Name: `:binomial`
;; * Default parameters
;;    * `:p`, probability: $0.5$
;;    * `:trials`: $20$
;; * [wiki](https://en.wikipedia.org/wiki/Binomial_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/BinomialDistribution.html)

(kind/table
 [[{:trials 20 :p 0.5} {:trials 50 :p 0.25}]
  [(gg/dgraphi (r/distribution :binomial) {:pdf [0 20]})
   (gg/dgraphi (r/distribution :binomial {:trials 50 :p 0.25}) {:pdf [0 30]})]])

;; #### Fisher's noncentral hypergeometric

;; * Name: `:fishers-noncentral-hypergeometric`
;; * Default parameters
;;    * `:ns`, number of sucesses: $10$
;;    * `:nf`, number of failures: $10$
;;    * `:n`, sample size, ($n<ns+nf$): $5$
;;    * `:omega`, odds ratio: $1$
;; * [wiki](https://en.wikipedia.org/wiki/Fisher%27s_noncentral_hypergeometric_distribution), [source](https://github.com/JuliaStats/Distributions.jl/blob/master/src/univariate/discrete/noncentralhypergeometric.jl)

(kind/table
 [[{:ns 10 :nf 10 :n 5 :omega 1} {:ns 30 :nf 60 :n 20 :omega 0.75}]
  [(gg/dgraphi (r/distribution :fishers-noncentral-hypergeometric) {:pdf [0 6]})
   (gg/dgraphi (r/distribution :fishers-noncentral-hypergeometric
                               {:ns 30 :nf 60 :n 20 :omega 0.75}) {:pdf [0 20]})]])

;; #### Geometric

;; * Name: `:geometric`
;; * Default parameters
;;    * `:p`, probability: $0.5$
;; * [wiki](https://en.wikipedia.org/wiki/Geometric_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/GeometricDistribution.html)

(kind/table
 [[{:p 0.5} {:p 0.15}]
  [(gg/dgraphi (r/distribution :geometric) {:pdf [0 10]})
   (gg/dgraphi (r/distribution :geometric {:p 0.15}) {:pdf [0 20]})]])

;; #### Hypergeometric

;; * Name: `:hypergeometric`
;; * Default parameters
;;    * `:population-size`: $100$
;;    * `:number-of-successes`: $50%
;;    * `:sample-size`: $25%
;; * [wiki](https://en.wikipedia.org/wiki/Hypergeometric_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/HypergeometricDistribution.html)

(kind/table
 [[{:population-size 100
    :number-of-successes 50
    :sample-size 25}]
  [(gg/dgraphi (r/distribution :hypergeometric) {:pdf [0 26]})]
  [{:population-size 2000
    :number-of-successes 20
    :sample-size 200}]
  [(gg/dgraphi (r/distribution :hypergeometric {:population-size 2000
                                                :number-of-successes 20
                                                :sample-size 200}) {:pdf [0 20]})]])

;; #### Logarithmic

;; * Name: `:logarithmic`
;; * Default parameters
;;    * `:theta`, shape: $0.5$
;; * [wiki](https://en.wikipedia.org/wiki/Logarithmic_distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdist_1_1LogarithmicDist.html)

(kind/table
 [[{:theta 0.5} {:theta 0.99}]
  [(gg/dgraphi (r/distribution :logarithmic) {:pdf [0 10]})
   (gg/dgraphi (r/distribution :logarithmic {:theta 0.9}) {:pdf [0 20]})]])

;; #### Negative binomial

;; * Name: `:negative-binomial`
;; * Default parameters
;;    * `:r`, number of successes: $20$, can be a real number.
;;    * `:p`, probability of success: $0.5$
;; * [wiki](https://en.wikipedia.org/wiki/Negative_binomial_distribution)

(kind/table
 [[{:r 20 :p 0.5} {:r 100 :p 0.95} {:r 21.2345 :p 0.7}]
  [(gg/dgraphi (r/distribution :negative-binomial) {:pdf [0 40]})
   (gg/dgraphi (r/distribution :negative-binomial {:r 100 :p 0.95}) {:pdf [0 20]})
   (gg/dgraphi (r/distribution :negative-binomial {:r 21.2345 :p 0.7}) {:pdf [0 30]})]])

;; #### Pascal

;; The same as `:negative-binomial` but `r` is strictly integer

;; * Name: `:pascal`
;; * Default parameters
;;    * `:r`, number of successes: $20$
;;    * `:p`, probability of success: $0.5$
;; * [wiki](https://en.wikipedia.org/wiki/Negative_binomial_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/PascalDistribution.html)

(kind/table
 [[{:r 20 :p 0.5} {:r 100 :p 0.95}]
  [(gg/dgraphi (r/distribution :pascal {:r 20 :p 0.5}) {:pdf [0 40]})
   (gg/dgraphi (r/distribution :pascal {:r 100 :p 0.95}) {:pdf [0 20]})]])

;; #### Poisson

;; * Name: `:poisson`
;; * Default parameters
;;    * `:p`, lambda, mean: $0.5$
;; * [wiki](https://en.wikipedia.org/wiki/Poisson_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/PoissonDistribution.html)

(kind/table
 [[{:p 0.5} {:p 4}]
  [(gg/dgraphi (r/distribution :poisson) {:pdf [0 5]})
   (gg/dgraphi (r/distribution :poisson {:p 4}) {:pdf [0 11]})]])

;; #### Uniform

;; * Name: `:uniform-int`
;; * Default parameters
;;    * `:lower`, lower bound: $0$
;;    * `:upper`, upper bound: $2147483647$
;; * [wiki](https://en.wikipedia.org/wiki/Discrete_uniform_distribution),[source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/UniformIntegerDistribution.html) 

(kind/table
 [[{:lower 0 :upper 20} {:lower -5 :upper 5}]
  [(gg/dgraphi (r/distribution :uniform-int {:upper 20}) {:pdf [-1 22]})
   (gg/dgraphi (r/distribution :uniform-int {:lower -5 :upper 5}) {:pdf [-6 7]})]])

;; #### Zero Adjusted Beta Binomial (zabb)

;; * Name: `:zabb`
;; * Default parameters
;;    * `:mu`, probability: $0.5$
;;    * `:sigma`, dispersion: $0.1$
;;    * `:nu`, probability at 0.0: $0.1$
;;    * `:bd`, binomial denominator: $1$
;; * [source](https://search.r-project.org/CRAN/refmans/gamlss.dist/html/ZABB.html), [book](https://www.gamlss.com/wp-content/uploads/2018/01/DistributionsForModellingLocationScaleandShape.pdf)

(kind/table
 [[{:mu 0.5 :sigma 0.1 :bd 10 :nu 0.1} {:nu 0.1 :mu 0.65 :sigma 0.3 :bd 20}]
  [(gg/dgraphi (r/distribution :zabb {:bd 10}) {:pdf [0 12]})
   (gg/dgraphi (r/distribution :zabb {:nu 0.1 :mu 0.65 :sigma 0.3 :bd 20}) {:pdf [0 22]})]])

;; #### Zero Adjusted Binomial (zabi)

;; * Name: `:zabi`
;; * Default parameters
;;    * `:mu`, probability: $0.5$
;;    * `:sigma`, probability at 0.0: $0.1$
;;    * `:bd`, binomial denominator: $1$
;; * [source](https://search.r-project.org/CRAN/refmans/gamlss.dist/html/ZABI.html), [book](https://www.gamlss.com/wp-content/uploads/2018/01/DistributionsForModellingLocationScaleandShape.pdf)

(kind/table
 [[{:mu 0.5 :sigma 0.1 :bd 10} {:mu 0.65 :sigma 0.3 :bd 20}]
  [(gg/dgraphi (r/distribution :zabi {:bd 10}) {:pdf [0 11]})
   (gg/dgraphi (r/distribution :zabi {:mu 0.65 :sigma 0.3 :bd 20}) {:pdf [0 21]})]])

;; #### Zero Adjusted Negative Binomial (zanbi)

;; * Name: `:zanbi`
;; * Default parameters
;;    * `:mu`, mean: $1.0$
;;    * `:sigma`, dispersion: $1.0$
;;    * `:nu`, probability at 0.0: $0.3$
;; * [source](https://search.r-project.org/CRAN/refmans/gamlss.dist/html/ZANBI.html), [book](https://www.gamlss.com/wp-content/uploads/2018/01/DistributionsForModellingLocationScaleandShape.pdf)

(kind/table
 [[{:mu 1 :sigma 1 :nu 0.3} {:nu 0.1 :mu 2 :sigma 0.5}]
  [(gg/dgraphi (r/distribution :zanbi) {:pdf [0 9]})
   (gg/dgraphi (r/distribution :zanbi {:nu 0.1 :mu 2 :sigma 0.5}) {:pdf [0 13]})]])

;; #### Zero Inflated Beta Binomial (zibb)

;; * Name: `:zibb`
;; * Default parameters
;;    * `:mu`, probability: $0.5$
;;    * `:sigma`, dispersion: $0.1$
;;    * `:nu`, probability factor at 0.0: $0.1$
;;    * `:bd`, binomial denominator: $1$
;; * [source](https://search.r-project.org/CRAN/refmans/gamlss.dist/html/ZABB.html), [book](https://www.gamlss.com/wp-content/uploads/2018/01/DistributionsForModellingLocationScaleandShape.pdf)

(kind/table
 [[{:mu 0.5 :sigma 0.5 :bd 10 :nu 0.1} {:nu 0.1 :mu 0.65 :sigma 1 :bd 20}]
  [(gg/dgraphi (r/distribution :zibb {:bd 10}) {:pdf [0 15]})
   (gg/dgraphi (r/distribution :zibb {:nu 0.1 :mu 0.65 :sigma 1 :bd 20}) {:pdf [0 21]})]])

;; #### Zero Inflated Binomial (zibi)

;; * Name: `:zibi`
;; * Default parameters
;;    * `:mu`, probability: $0.5$
;;    * `:sigma`, probability factor at 0.0: $0.1$
;;    * `:bd`, binomial denominator: $1$
;; * [source](https://search.r-project.org/CRAN/refmans/gamlss.dist/html/ZABI.html), [book](https://www.gamlss.com/wp-content/uploads/2018/01/DistributionsForModellingLocationScaleandShape.pdf)

(kind/table
 [[{:mu 0.5 :sigma 0.1 :bd 10} {:mu 0.65 :sigma 0.3 :bd 20}]
  [(gg/dgraphi (r/distribution :zibi {:bd 10}) {:pdf [0 11]})
   (gg/dgraphi (r/distribution :zibi {:mu 0.65 :sigma 0.3 :bd 20}) {:pdf [0 21]})]])

;; #### Zero Inflated Negative Binomial (zinbi)

;; * Name: `:zinbi`
;; * Default parameters
;;    * `:mu`, mean: $1.0$
;;    * `:sigma`, dispersion: $1.0$
;;    * `:nu`, probability factor at 0.0: $0.3$
;; * [source](https://search.r-project.org/CRAN/refmans/gamlss.dist/html/ZANBI.html), [book](https://www.gamlss.com/wp-content/uploads/2018/01/DistributionsForModellingLocationScaleandShape.pdf)

(kind/table
 [[{:mu 1 :sigma 1 :nu 0.3} {:nu 0.1 :mu 2 :sigma 0.5}]
  [(gg/dgraphi (r/distribution :zinbi) {:pdf [0 9]})
   (gg/dgraphi (r/distribution :zinbi {:nu 0.1 :mu 2 :sigma 0.5}) {:pdf [0 13]})]])

;; #### Zero Inflated Poisson (zip)

;; * Name: `:zip`
;; * Default parameters
;;    * `:mu`, mean: $5$
;;    * `:sigma`, probability at 0.0: $0.1$
;; * [source](https://search.r-project.org/CRAN/refmans/gamlss.dist/html/ZIP.html), [book](https://www.gamlss.com/wp-content/uploads/2018/01/DistributionsForModellingLocationScaleandShape.pdf)

(kind/table
 [[{:mu 5 :sigma 0.1} {:mu 2 :sigma 0.5}]
  [(gg/dgraphi (r/distribution :zip) {:pdf [0 14]})
   (gg/dgraphi (r/distribution :zip {:mu 2 :sigma 0.5}) {:pdf [0 7]})]])

;; #### Zero Inflated Poisson, type 2 (zip2)

;; * Name: `:zip2`
;; * Default parameters
;;    * `:mu`, mean: $5$
;;    * `:sigma`, probability at 0.0: $0.1$
;; * [source](https://search.r-project.org/CRAN/refmans/gamlss.dist/html/ZIP2.html), [book](https://www.gamlss.com/wp-content/uploads/2018/01/DistributionsForModellingLocationScaleandShape.pdf)

(kind/table
 [[{:mu 5 :sigma 0.1} {:mu 2 :sigma 0.5}]
  [(gg/dgraphi (r/distribution :zip2) {:pdf [0 14]})
   (gg/dgraphi (r/distribution :zip2 {:mu 2 :sigma 0.5}) {:pdf [0 9]})]])

;; #### Zipf

;; * Name: `:zipf`
;; * Default parameters
;;    * `:number-of-elements`: $100$
;;    * `:exponent`: $3.0$
;; * [wiki](https://en.wikipedia.org/wiki/Zipf%27s_law), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/ZipfDistribution.html)

(kind/table
 [[{:number-of-elements 100 :exponent 3} {:number-of-elements 20 :exponent 0.5}]
  [(gg/dgraphi (r/distribution :zipf) {:pdf [0 8]})
   (gg/dgraphi (r/distribution :zipf {:number-of-elements 20 :exponent 0.5}) {:pdf [0 21]})]])

;; ### Multivariate

(kind/table
 {:column-names ["name" "parameters" "continuous?"]
  :row-vectors (concat (map #(conj % true) (mapv (fn [[k v]] [(kind/code (str k)) v])
                                                 (build-distribution-list true true)))
                       (map #(conj % false) (mapv (fn [[k v]] [(kind/code (str k)) v])
                                                  (build-distribution-list true false))))})

;; #### Dirichlet

;; * Name: `:dirichlet`
;; * Default parameters
;;    * `:alpha`, concentration, vector: `[1 1]`
;; * [wiki](https://en.wikipedia.org/wiki/Dirichlet_distribution)

;; Please note, `PDF` doesn't validate input.

^:kindly/hide-code
(defn dirichlet-1d [alpha]
  (let [d (r/distribution :dirichlet {:alpha alpha})]
    (gg/fgraph (fn [^double x]
                 (let [v [x (- 1.0 x)]]
                   (r/pdf d v))) [0.0 1.0])))

^:kindly/hide-code
(defn dirichlet-2d [alpha]
  (let [d (r/distribution :dirichlet {:alpha alpha})]
    (gg/graph2d (fn [^double x ^double y]
                  (let [lv (- 1.0 x y)
                        v [x y lv]]
                    (if (neg? lv)
                      ##NaN
                      (r/pdf d v))))
                nil nil {:varg? false})))

;; Projections of the 2d and 3d Dirichlet distributions.

;; * 2d case - all vectors $[x,1-x]$
;; * 3d case - all (supported) vectors $[x,y,1-x-y]$

(kind/table
 [[{:alpha [0.6 0.6]} {:alpha [3 3]} {:alpha [0.5 2]}]
  [(dirichlet-1d [0.6 0.6]) (dirichlet-1d [3 3]) (dirichlet-1d [0.5 2])]
  [{:alpha [3 1 3]} {:alpha [3 3 3]} {:alpha [1 3 1]}]
  [(dirichlet-2d [3 1 3]) (dirichlet-2d [3 3 3]) (dirichlet-2d [1 3 1])]])

(def dirichlet3 (r/distribution :dirichlet {:alpha [3 1 3]}))

(utls/examples-note
 (r/sample dirichlet3)
 (r/sample dirichlet3)
 (r/pdf dirichlet3 [0.2 0.3 0.5])
 (r/pdf dirichlet3 [0.3 0.5 0.2])
 (r/pdf dirichlet3 [0.5 0.2 0.3])
 (r/means dirichlet3)
 (r/covariance dirichlet3))

;; #### Multi normal

;; * Name: `:multi-normal`
;; * Default parameters
;;    * `:means`, vector: `[0 0]`
;;    * `:covariances`, vector of vectors (row-wise matrix): `[[1 0] [0 1]]`
;; * [wiki](https://en.wikipedia.org/wiki/Multivariate_normal_distribution), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/MultivariateNormalDistribution.html)

^:kindly/hide-code
(defn multi-normal-2d [means covariances]
  (let [d (r/distribution :multi-normal {:means means :covariances covariances})]
    (gg/graph2d (fn [^double x ^double y]
                  (r/pdf d [x y])) [-3 3] [-3 3] {:varg? false})))

(kind/table
 [[{:means [0 0] :convariances [[1 0] [0 1]]}]
  [(multi-normal-2d nil nil)]
  [{:means [0.5 0] :convariances [[0.5 -0.1] [0.1 0.1]]}]
  [(multi-normal-2d [0.5 0] [[0.5 -0.1] [0.1 0.1]])]
  [{:means [0 -0.5] :convariances [[1 0.2] [0.3 1]]}]
  [(multi-normal-2d [0 -0.5] [[1 0.2] [0.3 1]])]])

;; #### Multinomial

;; * Name: `:multinomial`
;; * Default parameters
;;    * `:ps`, probabilities or weights, vector: `[0.5 0.5]`
;;    * `:trials`: $20$
;; * [wiki](https://en.wikipedia.org/wiki/Multinomial_distribution), [source](http://umontreal-simul.github.io/ssj/docs/master/classumontreal_1_1ssj_1_1probdistmulti_1_1MultinomialDist.html)

(def multinomial (r/distribution :multinomial {:trials 150 :ps [1 2 3 4 5]}))

(utls/examples-note
 (r/sample multinomial)
 (r/sample multinomial)
 (r/pdf multinomial [10 10 10 10 110])
 (r/pdf multinomial [10 20 30 40 50])
 (r/pdf multinomial [110 10 10 10 10])
 (r/means multinomial)
 (r/covariance multinomial))

;; ### Mixture

(kind/table
 {:column-names ["name" "parameters"]
  :row-vectors [[(kind/code ":mixture") (r/distribution-parameters (r/distribution :mixture) true)]]})

;; Mixture distribution  

;; Creates a distribution from other distributions and weights.

;; * Name: `:mixture`
;; * Default parameters:
;;    * `:distrs`, list of distributions: `[default-normal]`
;;    * `:weights`, list of weights: `[1.0]`
;; * [wiki](https://en.wikipedia.org/wiki/Mixture_distribution)

;; Please note: `set-seed!` doesn't affect distributions which are part of the mixture 

(def three-normals
  (r/distribution :mixture {:distrs [(r/distribution :normal {:mu -2 :sd 2})
                                     (r/distribution :normal)
                                     (r/distribution :normal {:mu 2 :sd 0.5})]
                            :weights [2 1 3]}))

(def mixture-of-three
  (r/distribution :mixture {:distrs [(r/distribution :gamma)
                                     (r/distribution :laplace)
                                     (r/distribution :log-logistic)]
                            :weights [2 1 3]}))

(kind/table
 [["three normals" "gamma, laplace and log-logistic"]
  [(gg/dgraph three-normals {:pdf [-5 5]})
   (gg/dgraph mixture-of-three {:pdf [-5 5]})]])

(utls/examples-note
 (r/sample mixture-of-three)
 (r/sample mixture-of-three)
 (r/pdf mixture-of-three 0)
 (r/pdf mixture-of-three 1)
 (r/pdf mixture-of-three 2)
 (r/cdf mixture-of-three 0)
 (r/cdf mixture-of-three 1)
 (r/cdf mixture-of-three 2)
 (r/icdf mixture-of-three 0.01)
 (r/icdf mixture-of-three 0.5)
 (r/icdf mixture-of-three 0.99)
 (r/mean mixture-of-three)
 (r/variance mixture-of-three)
 (r/lower-bound mixture-of-three)
 (r/upper-bound mixture-of-three))

;; ### Truncated

(kind/table
 {:column-names ["name" "parameters"]
  :row-vectors [[:truncated (r/distribution-parameters (r/distribution :truncated) true)]]})

;; * Name: `:truncated`
;; * Default parameters
;;    * `:distr`, distribution to truncate: `default-normal`
;;    * `:left`, lower boundary
;;    * `:right`, upper boundary
;; * [wiki](https://en.wikipedia.org/wiki/Truncated_distribution)

;; By default boundaries are the same as boundaries from a distributions. This way you can make one side truncation.

;; Please note: derived `mean` or `variance` is not calculated. Also, `set-seed!` doesn't affect original distribution.

(def truncated-normal (r/distribution :truncated {:distr r/default-normal
                                                  :left -2 :right 2}))

(def left-truncated-laplace (r/distribution :truncated {:distr (r/distribution :laplace)
                                                        :left -0.5}))

(kind/table
 [["truncated normal" "trucated levy (left side)"]
  [(gg/dgraph truncated-normal {:pdf [-2.5 2.5]})
   (gg/dgraph left-truncated-laplace {:pdf [-1 6]})]])

(utls/examples-note
 (r/sample truncated-normal)
 (r/sample truncated-normal)
 (r/pdf truncated-normal 2.0000001)
 (r/pdf truncated-normal 2.0)
 (r/cdf truncated-normal 2.0)
 (r/cdf truncated-normal 0.0)
 (map (partial r/icdf truncated-normal) [0.0001 0.5 0.9999])
 (r/lower-bound truncated-normal)
 (r/upper-bound truncated-normal))


;; ### From data

;; All below distributions can be constructed from datasets or list of values with probabilities.

(kind/table
 {:column-names ["name" "parameters" "continuous?"]
  :row-vectors (mapv (fn [[k v]]
                       [(kind/code (str k)) v])
                     [[:continuous-distribution (r/distribution-parameters (r/distribution :continuous-distribution) true) true]
                      [:kde (r/distribution-parameters (r/distribution :kde) true) true]
                      [:empirical (r/distribution-parameters (r/distribution :empirical) true) true]
                      [:real-discrete-distribution (r/distribution-parameters (r/distribution :real-discrete-distribution) true) false]
                      [:integer-discrete-distribution (r/distribution-parameters (r/distribution :integer-discrete-distribution) true) false]
                      [:categorical-distribution (r/distribution-parameters (r/distribution :categorical-distribution) true) false]
                      [:enumerated-real (r/distribution-parameters (r/distribution :enumerated-real) true) false]
                      [:enumerated-int (r/distribution-parameters (r/distribution :enumerated-int) true) false]])})


;; #### Continuous

;; Continous distribution build from data is based on KDE (Kernel Density Estimation) and PDF integration for CDF and iCDF. Mean and variance are calculated from samples.

;; `:continous-distribution` and `:kde` are two names for the same distribution

;; * Name: `:continuous-distribution` or `:kde`
;; * Default parameters:
;;    * `:data`, samples, sequence of numbers
;;    * `:kde`, density estimation kernel: `:epenechnikov`
;;    * `:bandwidth`, KDE bandwidth, smoothing parameter: `nil` (auto)
;;    * `:steps`, number of steps for PDF integration: 5000
;;    * `:min-iterations`, number of PDF integrator iterations: `nil` (default)
;;    * `:interpolator`, CDF/iCDF interpolator for PDF integration: `nil` (default, `:linear`)
;; * [wiki](https://en.wikipedia.org/wiki/Kernel_density_estimation), [pdf integration](#pdf-integration)

(def random-data (repeatedly 1000 (fn [] (+ (* (r/drand -2 2) (r/drand -2 2))
                                            (m/sqrt (* (r/drand) (r/drand)))))))

(def kde-distr (r/distribution :continuous-distribution {:data random-data}))

(let [build-distr (fn [pars]
                    (r/distribution :continuous-distribution (assoc pars :data random-data)))]
  (kind/table
   [["default", "bandwidth=1.0"]
    [(gg/dgraph kde-distr {:pdf [-5 5]})
     (gg/dgraph (build-distr {:bandwidth 1.0}) {:pdf [-5 5]})]
    ["gaussian kernel" "triangular kernel, bandwidth=0.1"]
    [(gg/dgraph (build-distr {:kde :gaussian}) {:pdf [-5 5]})
     (gg/dgraph (build-distr {:kde :triangular :bandwidth 0.1}) {:pdf [-5 5]})]]))

(utls/examples-note
 (r/sample kde-distr)
 (r/sample kde-distr)
 (r/pdf kde-distr 0.0)
 (r/cdf kde-distr 0.0)
 (map (partial r/icdf kde-distr) [0.0001 0.5 0.9999])
 (r/lower-bound kde-distr)
 (r/upper-bound kde-distr)
 (r/mean kde-distr)
 (r/variance kde-distr))

;; ##### Kernels

;; Distributions for various kernels:
;; * `:data`: `[-2 -2 -2 -1 0 1 2 -1 0 1 2 0 1 2 1 2 2]`
;; * `:steps`: `100`
;; * `:bandwidth`: auto

;; `fastmath.kernel/kernel-list` contains three types of kernels: RBF, KDE and what we can call "vector kernels" which includes Marcer, positive definite, and similar.

(kind/table
 {:column-names ["Type" "Kernels"]
  :row-vectors (->> k/kernel-list
                    (map (fn [[k v]]
                           [(kind/code (str k)) v]))
                    (into []))})

;; KDE kernels

^:kindly/hide-code
(def kde-kernels (:kde k/kernel-list))

(kind/table
 (for [k kde-kernels]
   [(kind/code (str k)) (gg/dgraph (r/distribution :continuous-distribution {:data [-2 -2 -2 -1 0 1 2 -1 0 1 2 0 1 2 1 2 2]
                                                                             :kde k :steps 100}) {:pdf [-3 3]})]))

;; #### Empirical

;; Empirical distribution calculates PDF, CDF and iCDF from a histogram.

;; * Name: `:empirical`
;; * Default parameters
;;    * `:data`, samples, sequence of numbers
;;    * `:bin-count`, number of bins for histogram: 10% of the size of the data
;; * [wiki](https://en.wikipedia.org/wiki/Empirical_distribution_function), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/random/EmpiricalDistribution.html)

(def empirical-distr (r/distribution :empirical {:data random-data}))

(kind/table
 [["default", "bin-count=10"]
  [(gg/dgraph empirical-distr {:pdf [-5 5]})
   (gg/dgraph (r/distribution :empirical {:data random-data :bin-count 10}) {:pdf [-5 5]})]])

(utls/examples-note
 (r/sample empirical-distr)
 (r/sample empirical-distr)
 (r/pdf empirical-distr 0.0)
 (r/cdf empirical-distr 0.0)
 (map (partial r/icdf empirical-distr) [0.0001 0.5 0.9999])
 (r/lower-bound empirical-distr)
 (r/upper-bound empirical-distr)
 (r/mean empirical-distr)
 (r/variance empirical-distr))

;; #### Discrete

;; * Default parameters
;;    * `:data`, sequence of numbers (integers/longs or doubles)
;;    * `:probabilities`, optional, probabilities or weights
;; * [wiki](https://en.wikipedia.org/wiki/Probability_distribution#Discrete_probability_distribution), [source1](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/EnumeratedRealDistribution.html), [source2](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/distribution/EnumeratedIntegerDistribution.html)

;; Please note: data can contain duplicates.

;; There are four discrete distributions:
;; * `:enumerated-int` for integers, backed by Apache Commons Math
;; * `:enumerated-real` for doubles, backed by Apache Commons Math
;; * `:integer-discrete-distribution` - for longs, custom implementation
;; * `:real-discrete-distribution` - for doubles, custom implementation

;; Please note:

;; * Apache Commons Math implementation have some issues with iCDF.
;; * `:integer-discrete-distribution` is backed by `clojure.data.int-map`

;; ##### Doubles

(def data-doubles (repeatedly 100 #(m/sq (m/approx (r/drand 2.0) 1))))

(kind/table
 [[(kind/code ":real-discrete-distribution") (kind/code ":enumerated-real")]
  (let [enumerated-real (r/distribution :enumerated-real {:data data-doubles})
        real-discrete (r/distribution :real-discrete-distribution {:data data-doubles})
        dist (distinct data-doubles)]
    [(gg/dgraphi real-discrete {:pdf [0 5] :data (into {} (map #(vector % (r/pdf real-discrete %))
                                                               dist))})
     (gg/dgraphi enumerated-real {:pdf [0 5] :data (into {} (map #(vector % (r/pdf enumerated-real %))
                                                                 dist))})])])

(def real-distr (r/distribution :real-discrete-distribution {:data [0.1 0.2 0.3 0.4 0.3 0.2 0.1]
                                                             :probabilities [5 4 3 2 1 5 4]}))

(utls/examples-note
 (r/sample real-distr)
 (r/sample real-distr)
 (r/pdf real-distr 0.1)
 (r/pdf real-distr 0.15)
 (r/cdf real-distr 0.2)
 (map (partial r/icdf real-distr) [0.0001 0.5 0.9999])
 (r/lower-bound real-distr)
 (r/upper-bound real-distr)
 (r/mean real-distr)
 (r/variance real-distr))

;; ##### Integers / Longs

(def data-ints (repeatedly 500 #(int (m/sqrt (r/drand 100.0)))))

(kind/table
 [[(kind/code ":integer-discrete-distribution") (kind/code ":enumerated-int")]
  (let [enumerated-int (r/distribution :enumerated-int {:data data-ints})
        int-discrete (r/distribution :integer-discrete-distribution {:data data-ints})
        dist (distinct data-ints)]
    [(gg/dgraphi int-discrete {:pdf [0 10] :data (into {} (map #(vector % (r/pdf int-discrete %))
                                                               dist))})
     (gg/dgraphi enumerated-int {:pdf [0 10] :data (into {} (map #(vector % (r/pdf enumerated-int %))
                                                                 dist))})])])

(def int-distr (r/distribution :integer-discrete-distribution {:data [10 20 30 40 30 20 10]
                                                               :probabilities [5 4 3 2 1 5 4]}))

(utls/examples-note
 (r/sample int-distr)
 (r/sample int-distr)
 (r/pdf int-distr 10)
 (r/pdf int-distr 15)
 (r/cdf int-distr 20)
 (map (partial r/icdf int-distr) [0.0001 0.5 0.9999])
 (r/lower-bound int-distr)
 (r/upper-bound int-distr)
 (r/mean int-distr)
 (r/variance int-distr))

;; #### Categorical

;; Categorical distribution is a discrete distribution which accepts any data.

;; * Name: `:categorical-distribution`
;; * Default parameters:
;;    * `:data`, sequence of any values
;;    * `:probabilities`, optional, probabilities or weights

;; Order for CDF/iCDF is created by calling `(distinct data)`. If sorted data is needed, external sort is necessary. `lower-bound` and `upper-bound` are not defined though.

(def cat-distr (r/distribution :categorical-distribution {:data (repeatedly 100 #(rand-nth [:a :b nil "s"]))}))

(utls/examples-note
 (r/sample cat-distr)
 (r/sample cat-distr)
 (r/pdf cat-distr nil)
 (r/pdf cat-distr "ss")
 (r/cdf cat-distr :b)
 (map (partial r/icdf cat-distr) [0.0001 0.5 0.9999]))

;; ## Sequences

;; Sequence generators can create random or quasi random vector with certain property like low discrepancy or from ball/sphere.

;; There are two multimethods:

(utls/symbol-info-table
 [[sequence-generator "Lazy sequence of generated vectors (for dim>1) or primitives (for dim=1)"]
  [jittered-sequence-generator "Adds jittering, works only for low discrepancy sequences"]])


^:kindly/hide-code
(kind/md "
Parameters:

   * `seq-generator` - generator name
   * `dimensions` - vector dimensionality, 1 for primitive
   * `jitter` - only for jittered sequences, from 0.0 to 1.0, default 0.25

For given dimensionality, returns sequence of:

   * 1 - doubles
   * 2 - `Vec2` type
   * 3 - `Vec3` type
   * 4 - `Vec4` type
   * n>4 - Clojure vector

`Vec2`, `Vec3` and `Vec4` are fixed size vectors optimized for speed. They act exactly like 2,3 and 4 elements Clojure vectors
")


;; ### Low discrepancy

;; There are 3 types of sequences:

;; * `:sobol` - up to 1000 dimensions, [wiki](https://en.wikipedia.org/wiki/Sobol_sequence), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/random/SobolSequenceGenerator.html)
;; * `:halton` - up to 40 dimensions, [wiki](https://en.wikipedia.org/wiki/Halton_sequence), [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/random/HaltonSequenceGenerator.html)
;; * `:r2` - up to 15 dimensions, [info](http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/)

;; 1000 samples from each of the sequence type without and with jittering

(kind/table
 [[":sobol" ":halton" ":r2"]
  [(gg/graph-scatter (take 1000 (r/sequence-generator :sobol 2)))
   (gg/graph-scatter (take 1000 (r/sequence-generator :halton 2)))
   (gg/graph-scatter (take 1000 (r/sequence-generator :r2 2)))]
  [":sobol (jittered)" ":halton (jittered)" ":r2 (jittered)"]
  [(gg/graph-scatter (take 1000 (r/jittered-sequence-generator :sobol 2)))
   (gg/graph-scatter (take 1000 (r/jittered-sequence-generator :halton 2)))
   (gg/graph-scatter (take 1000 (r/jittered-sequence-generator :r2 2)))]])

(utls/examples-note
 (first (r/sequence-generator :sobol 4))
 (first (r/jittered-sequence-generator :sobol 4))
 (first (r/sequence-generator :halton 3))
 (first (r/jittered-sequence-generator :halton 3))
 (first (r/sequence-generator :r2 2))
 (first (r/jittered-sequence-generator :r2 2)))

;; 15 dimensional sequence
(take 2 (r/sequence-generator :r2 15))

;; One dimensional sequence is just a sequence of numbers

(take 20 (r/sequence-generator :sobol 1))

;; ### Sphere and ball

;; Unit sphere or unit ball sequences can generate any dimension.

;; * `:sphere` - [source](https://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/org/apache/commons/math3/random/UnitSphereRandomVectorGenerator.html)
;; * `:ball` - dropped coordinates, [info]( http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls)

;; 500 samples

(kind/table
 [[(kind/code ":sphere") (kind/code ":ball")]
  [(gg/graph-scatter (take 1000 (r/sequence-generator :sphere 2)) {:aspect-ratio 1})
   (gg/graph-scatter (take 500 (r/sequence-generator :ball 2)) {:aspect-ratio 1})]])


(utls/examples-note
 (first (r/sequence-generator :sphere 4))
 (first (r/sequence-generator :ball 3)))

;; 20 dimensional sequence

(take 2 (r/sequence-generator :sphere 20))

;; ### Uniform and Gaussian

;; Additionally uniform and gaussian N(0,1) sequences can generate any number of dimensions. They rely on `default-rng`

;; * `:default` - uniform distribution U(0,1)
;; * `:gaussian` - gaussian, normal distribution, N(0,1) for each dimension

;; 1000 samples

(kind/table
 [[":default" ":gaussian"]
  [(gg/graph-scatter (take 1000 (r/sequence-generator :default 2)))
   (gg/graph-scatter (take 1000 (r/sequence-generator :gaussian 2)))]])

(utls/examples-note
 (first (r/sequence-generator :default 4))
 (first (r/sequence-generator :gaussian 3)))

;; ## Noise

;; Value, gradient and simplex noises plus various combinations. 1d ,2d and 3d versions are prepared.

;; ### Generation

;; There are four main methods of noise creation:

(utls/symbol-info-table
 [[single-noise "single frequency (octave) noise"]
  [fbm-noise "multi frequency (octaves), fractal brownian motion"]
  [billow-noise "multi frequency, \"billowy\" noise"]
  [ridgemulti-noise "multi frequency, ridged multi-fractal"]])

;; Each noise can be configured in, here is the list of options:

;; * `:seed` - seed for noise randomness
;; * `:noise-type`
;;    * `:value` - value noise
;;    * `:gradient` (default) - gradient noise (Perlin)
;;    * `:simplex` - OpenSimplex noise
;; * `:interpolation` - interpolation between knots, only for value and gradient noise
;;    * `:none`
;;    * `:linear`
;;    * `:hermite` (default)
;;    * `:quintic`
;; * `:octaves` (default: 6) - number of frequencies/octaves for multi frequency creators
;; * `:lacunarity` (default: 2) - noise length (1/frequency) for each octave
;; * `:gain` (default: 0.5) - amplitude factor for each octave
;; * `:normalize?` (default: true) - if true, range is `[0,1]`, `[-1,1]` otherwise.

;; [more info](http://www.campi3d.com/External/MariExtensionPack/userGuide5R4v1/Understandingsomebasicnoiseterms.html) about octaves, gain and lacunarity.

;; #### Single

(def single-g-noise (r/single-noise {:noise-type :gradient :seed 1}))

(utls/examples-note
 (single-g-noise 0.2)
 (single-g-noise 0.2 0.3)
 (single-g-noise 0.2 0.3 0.4))

;; Single octave of simplex noise:

(gg/graph2d (r/single-noise {:seed 1 :noise-type :simplex}) [-2 2] [-2 2] {:varg? false})

;; Value and gradient single noise for different interpolations
(kind/table
 [["" (kind/code ":none") (kind/code ":linear") (kind/code ":hermite") (kind/code ":quintic")]
  [:value
   (gg/graph2d (r/single-noise {:seed 1 :noise-type :value :interpolation :none}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/single-noise {:seed 1 :noise-type :value :interpolation :linear}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/single-noise {:seed 1 :noise-type :value :interpolation :hermite}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/single-noise {:seed 1 :noise-type :value :interpolation :quintic}) [-2 2] [-2 2] {:varg? false})]
  [:gradient
   (gg/graph2d (r/single-noise {:seed 1 :noise-type :gradient :interpolation :none}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/single-noise {:seed 1 :noise-type :gradient :interpolation :linear}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/single-noise {:seed 1 :noise-type :gradient :interpolation :hermite}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/single-noise {:seed 1 :noise-type :gradient :interpolation :quintic}) [-2 2] [-2 2] {:varg? false})]])

;; #### FBM

(def fbm-noise (r/fbm-noise {:noise-type :gradient :octaves 3 :seed 1}))

(utls/examples-note
 (fbm-noise 0.2)
 (fbm-noise 0.2 0.3)
 (fbm-noise 0.2 0.3 0.4))

;; 6 octave of simplex noise:

(gg/graph2d (r/fbm-noise {:seed 1 :noise-type :simplex}) [-2 2] [-2 2] {:varg? false})

;; Value and gradient FBM noise for different interpolations

(kind/table
 [(mapv (comp kind/code str) ["" :none :linear :hermite :quintic])
  [:value
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :value :interpolation :none}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :value :interpolation :linear}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :value :interpolation :hermite}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :value :interpolation :quintic}) [-2 2] [-2 2] {:varg? false})]
  [:gradient
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :interpolation :none}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :interpolation :linear}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :interpolation :hermite}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :interpolation :quintic}) [-2 2] [-2 2] {:varg? false})]])

;; Different number of octaves for FBM gradient noise

(kind/table
 [["octaves=2" "octaves=4" "octaves=6" "octaves=8"]
  [(gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :octaves 2}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :octaves 4}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :octaves 6}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :octaves 8}) [-2 2] [-2 2] {:varg? false})]])

;; Different gains and lacunarities for FBM gradient noise

(kind/table
 [["" "lacunarity=0.5" "lacunarity=2" "lacunarity=5" "lacunarity=8"]
  ["gain=0.25"
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :gain 0.25 :lacunarity 0.5}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :gain 0.25 :lacunarity 2}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :gain 0.25 :lacunarity 5}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :gain 0.25 :lacunarity 8}) [-2 2] [-2 2] {:varg? false})]
  ["gain=0.5"
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :lacunarity 0.5}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :lacunarity 2}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :lacunarity 5}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :lacunarity 8}) [-2 2] [-2 2] {:varg? false})]
  ["gain=0.75"
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :gain 0.75 :lacunarity 0.5}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :gain 0.75 :lacunarity 2}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :gain 0.75 :lacunarity 5}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/fbm-noise {:seed 1 :noise-type :gradient :gain 0.75 :lacunarity 8}) [-2 2] [-2 2] {:varg? false})]])

;; #### Billow

(def billow-noise (r/billow-noise {:seed 1}))

(gg/graph2d billow-noise [-2 2] [-2 2] {:varg? false})

(kind/table
 [["simplex noise" "value noise" "gradient noise, 1 octave"]
  [(gg/graph2d (r/billow-noise {:seed 1 :noise-type :simplex}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/billow-noise {:seed 1 :noise-type :value}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/billow-noise {:seed 1 :noise-type :gradient :octaves 1}) [-2 2] [-2 2] {:varg? false})]])

(utls/examples-note
 (billow-noise 0.2)
 (billow-noise 0.2 0.3)
 (billow-noise 0.2 0.3 0.4))

;; #### Ridged Multi

(def ridgedmulti-noise (r/ridgedmulti-noise {:seed 1}))

(gg/graph2d ridgedmulti-noise [-2 2] [-2 2] {:varg? false})

(kind/table
 [["simplex noise" "value noise" "gradient noise, 1 octave"]
  [(gg/graph2d (r/ridgedmulti-noise {:seed 1 :noise-type :simplex}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/ridgedmulti-noise {:seed 1 :noise-type :value}) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/ridgedmulti-noise {:seed 1 :noise-type :gradient :octaves 1}) [-2 2] [-2 2] {:varg? false})]])

(utls/examples-note
 (ridgedmulti-noise 0.2)
 (ridgedmulti-noise 0.2 0.3)
 (ridgedmulti-noise 0.2 0.3 0.4))

;; ### Predefined

;; There are three ready to use preconfigured noises:

(utls/symbol-info-table
 [[vnoise "FBM value noise, 6 octaves, hermite interpolation"]
  [noise "Perlin noise, FBM gradient noise, 6 octaves, quintic interpolation"]
  [simplex "FBM simplex noise"]])

(kind/table
 [["vnoise" "noise" "simplex"]
  [(gg/graph2d r/vnoise [-2 2] [-2 2] {:varg? false})
   (gg/graph2d r/noise [-2 2] [-2 2] {:varg? false})
   (gg/graph2d r/simplex [-2 2] [-2 2] {:varg? false})]])

(utls/examples-note
 (r/vnoise 0.2)
 (r/vnoise 0.2 0.3)
 (r/vnoise 0.2 0.3 0.4)
 (r/noise 0.2)
 (r/noise 0.2 0.3)
 (r/noise 0.2 0.3 0.4)
 (r/simplex 0.2)
 (r/simplex 0.2 0.3)
 (r/simplex 0.2 0.3 0.4))

;; ### Warping

;; Warp noise [info](http://www.iquilezles.org/www/articles/warp/warp.htm)

(utls/symbol-info-table
 [[warp-noise-fn "Create warp noise"]])

;; Default parameters:
;;  * `noise` - any noise function: `vnoise`
;;  * `scale`: $4.0$
;;  * `depth`: $1$

(kind/table
 [["" "vnoise" "noise" "simplex"]
  ["scale=2"
   (gg/graph2d (r/warp-noise-fn r/vnoise 2.0 1) [-2 2] [-2 2]  {:varg? false})
   (gg/graph2d (r/warp-noise-fn r/noise 2.0 1) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/warp-noise-fn r/simplex 2.0 1) [-2 2] [-2 2] {:varg? false})]
  ["scale=4"
   (gg/graph2d (r/warp-noise-fn r/vnoise 4.0 1) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/warp-noise-fn r/noise 4.0 1) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/warp-noise-fn r/simplex 4.0 1) [-2 2] [-2 2] {:varg? false})]])

;; ### Random configuration

;; For generative art purposes it's good to generate random configuration and noise based on it.

(utls/symbol-info-table
 [[random-noise-cfg "Create random configuration"]
  [random-noise-fn "Create random noise from random configuration"]])

;; Optional parameter is a map with values user wants to fix.

(r/random-noise-cfg)
(r/random-noise-cfg {:seed 1})

(def some-random-noise (r/random-noise-fn {:seed 1}))

(gg/graph2d some-random-noise [-2 2] [-2 2] {:varg? false})

(utls/examples-note
 (some-random-noise 0.2)
 (some-random-noise 0.2 0.3)
 (some-random-noise 0.2 0.3 0.4))

(kind/table
 [[(gg/graph2d (r/random-noise-fn) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/random-noise-fn) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/random-noise-fn) [-2 2] [-2 2] {:varg? false})]
  [(gg/graph2d (r/random-noise-fn) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/random-noise-fn) [-2 2] [-2 2] {:varg? false})
   (gg/graph2d (r/random-noise-fn) [-2 2] [-2 2] {:varg? false})]])

;; ### Discrete noise

;; Discrete noise is a function which hashes long or two longs and converts it to a double from `[0,1]` range.

(utls/symbol-info-table
 [[discrete-noise "1d or 2d hashing function"]])

(utls/examples-note
 (r/discrete-noise 100)
 (r/discrete-noise 101)
 (r/discrete-noise 200 100)
 (r/discrete-noise 200 101))

;; ## Reference

(codox/make-public-fns-table-clay 'fastmath.random)
