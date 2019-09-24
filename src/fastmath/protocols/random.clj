(ns fastmath.protocols.random
  "Set of protocols for fastmath.random namespace.

  Includes:

  * random generator protocol
  * distribution protocols")

(defprotocol RNGProto
  "Defines set of random functions for different RNGs or distributions returning primitive values."
  (irandom [rng] [rng mx] [rng mn mx]
    "Random integer.

As default returns random integer from full integer range. 
When `mx` is passed, range is set to `[0, mx)`. When `mn` is passed, range is set to `[mn, mx)`.")
  (drandom [rng] [rng mx] [rng mn mx]
    "Random double.

As default returns random double from `[0,1)` range.
When `mx` is passed, range is set to `[0, mx)`. When `mn` is passed, range is set to `[mn, mx)`.")
  (lrandom [rng] [rng mx] [rng mn mx]
    "Random long.

As default returns random long from full long range. 
When `mx` is passed, range is set to `[0, mx)`. When `mn` is passed, range is set to `[mn, mx)`.")
  (frandom [rng] [rng mx] [rng mn mx]
    "Random float.

As default returns random float from `[0,1)` range.
When `mx` is passed, range is set to `[0, mx)`. When `mn` is passed, range is set to `[mn, mx)`.")
  (grandom [rng] [rng std] [rng mean std]
    "Random double from gaussian distribution.

As default returns random double from `N(0,1)`. 
When `std` is passed, `N(0,std)` is used. When `mean` is passed, distribution is set to `N(mean, std)`.")
  (brandom [rng] [rng p]
    "Boolean random.

Returns true or false with equal probability. You can set `p` probability for `true`")
  (set-seed! [rng v] "Sets seed. Returns `rng`")
  (->seq [rng] [rng n] "Returns lazy sequence of random samples (can be limited to optional `n` values)."))

(defprotocol DistributionProto
  "Get information from distributions."
  (cdf [d v] [d v1 v2] "Cumulative probability.")
  (pdf [d v] "Density")
  (lpdf [d v] "Log density")
  (icdf [d p] "Inverse cumulative probability")
  (probability [d v] "Probability (PMF)")
  (sample [d] "Returns random sample.")
  (dimensions [d] "Returns dimensions")
  (source-object [d] "Returns Java object from backend library")
  (continuous? [d] "Does distribution support continuous range?"))

(defprotocol UnivariateDistributionProto
  (mean [d] "Mean")
  (variance [d] "Variance")
  (lower-bound [d] "Lower value")
  (upper-bound [d] "Higher value"))

(defprotocol MultivariateDistributionProto
  "Get information from distributions."
  (means [d] "Mean")
  (covariance [d] "Variance"))

(defprotocol DistributionIdProto
  "Get name and parameter names from distribution"
  (distribution-id [d] "Distribution id as keyword")
  (distribution-parameters [d] "List of distribution parameter names"))
