(ns fastmath.protocols
  "Set of protocols for fastmath.

  Includes:

  * random generator protocol
  * distribution protocols
  * vector protocol"
  (:refer-clojure :exclude [abs]))

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
  (set-seed [rng v] "Sets seed. Returns new `rng` object")
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
  (distribution-parameters [d] "List of distribution parameter names")
  (distribution? [d] "Returns true if the object is distribution"))

;; vector

(defprotocol VectorProto
  "Vector operations"
  (to-vec [v] "Convert to Clojure primitive vector `Vec`.")
  (to-acm-vec [v] "Convert to Apache Commons Math ArrayRealVector")
  (as-vec [v] [v xs] "Create vector from sequence as given type.")
  (fmap [v f] "Apply function to all vector values (like map but returns the same type).")
  (approx [v] [v d] "Round to 2 (or `d`) decimal places")
  (magsq [v1] "Length of the vector squared.")
  (mag [v1] "length of the vector.")
  (dot [v1 v2] "Dot product of two vectors.")
  (add [v1] [v1 v2] "Sum of two vectors.")
  (sub [v1] [v1 v2] "Subtraction of two vectors.")
  (shift [v1 v] "Add `v` value to every vector element.")
  (mult [v1 v] "Multiply vector by number `v`.")
  (emult [v1 v2] "Element-wise vector multiplication (Hadamard product).")
  (abs [v1] "Absolute value of vector elements")
  (mx [v1] "Maximum value of vector elements")
  (mn [v1] "Minimum value of vector elements")
  (emx [v1 v2] "Element-wise max from two vectors.")
  (emn [v1 v2] "Element-wise min from two vectors.")
  (maxdim [v] "Index of maximum value.")
  (mindim [v] "Index of minimum value.")
  (base-from [v] "List of perpendicular vectors (basis)")
  (sum [v1] "Sum of elements")
  (prod [v1] "Product of elements")
  (permute [v idxs] "Permute vector elements with given indices.")
  (reciprocal [v] "Reciprocal of elements.")
  (interpolate [v1 v2 t f] "Interpolate vectors, optionally set interpolation fn")
  (einterpolate [v1 v2 v f] "Interpolate vectors element-wise, optionally set interpolation fn")
  (econstrain [v val1 val2] "Element-wise constrain")
  (is-zero? [v1] "Is vector zero?")
  (is-near-zero? [v1] [v1 tol] "Is vector almost zero? (all absolute values of elements are less than `tol` tolerance or `1.0e-6`)")
  (heading [v1] "Angle between vector and unit vector `[1,0,...]`")
  (cross [v1 v2] "Cross product")
  (rotate [v1 angle] [v1 anglex angley anglez] "Rotate vector")
  (perpendicular [v1] [v1 v2] "Perpendicular vector (only 2d).")
  (axis-rotate [v1 angle axis] [v1 angle axis pivot] "Rotate around axis, 3d only")
  (transform [v1 o vx vy] [v1 o vx vy vz] "Transform vector; map point to coordinate system defined by origin, vx and vy (as bases), d and 3d only.")
  (to-polar [v1] "To polar coordinates (2d, 3d only), first element is length, the rest angle.")
  (from-polar [v1] "From polar coordinates (2d, 3d only)"))

(defprotocol TransformProto
  "Transformer functions."
  (forward-1d [t xs] "Forward transform of sequence or array.")
  (reverse-1d [t xs] "Reverse transform of sequence or array.")
  (forward-2d [t xss] "Forward transform of sequence of sequences.")
  (reverse-2d [t xss] "Reverse transform of sequence of sequences."))

(defprotocol GridProto
  "Common grid conversion functions."
  (coords->cell [g coords] [g x y] "Converts 2d space coordinates to cell coordinates.")
  (cell->anchor [g cell] [g q r] "Converts cell coordinates to anchor coordinates.")
  (coords->mid [g coords] [g x y] "Converts 2d space into cell midpoint.")
  (cell->mid [g cell] [g q r] "Converts cell coordinates to cell midpoint")
  (grid-type [g] "Returns type of the cell.")
  (corners [g coords] [g coords scale] [g x y scale] "Returns list of cell vertices for given 2d space coordinates."))
