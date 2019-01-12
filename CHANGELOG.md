# Change Log
All notable changes to this project will be documented in this file. This change log follows the conventions of [keepachangelog.com](http://keepachangelog.com/).

## [1.3.0]

### Added

* `make-vector` returns vector for given number of dimensions and optional sequence. 
* vectors implements now `IPersistentVector` to work with `vector?`
* core.matrix protocols for vectors
* math functions can operate on vectors now (like [[sin]] etc.)
* new functions for vectors `clamp`, `zero-count`, `nonzero-count`, `as-vec`
* classification bindings for SMILE, XGBoost and LIBLINEAR

### Changed

* use `fmap` instead of `applyf` (now deprecated)
* rbf function already implements Smile interface

### Removed

* rbf-obj converter

## [1.2.0]

### Added

* Various grid operations
* hashCode for vectors

### Changed

* Breaking: sequence generator creatators return lazy sequence now (instead of function returning sequence)

## [1.1.4]

### Added

* `lcm` - least common multiplier
* Vec types implement `Reversible` `Indexed` `ILookup`
* negative values for some constants (`PI`, `E` etc.)

### Fixed

* `GAMMA` constant name clash

### Changed

* `stats-map` returns list of outliers rather than number of outliers

### Removed

* `:bessel` vector field removed (hard to limit input range)

## [1.1.3]

### Added

* `sample` can now return pairs of `[x,(f x)]`.

## [1.1.2]

### Added

* estimators for number of bins in histogram
* kernel-density function

## [1.1.1]

### Changed

* histogram data enhanced

## [1.1.0]

### Added

* `fastmath.clustering` namespace

### Changed

* removed `k-means` from `fastmath.stats`

## [1.0.3]

### Fixed

* `norm` when domain is a point, returns range start when value is less or equal domain, range end otherwise
* small fix in one field

## [1.0.2]

### Added

* curl for 2d vector fields

## [1.0.1]

### Added

* :default (:jdk) rng is added

### Fixed

* MersenneTwister is not synchronized, default RNG is :jdk now

## [1.0.0]

### Added

* easings namespace
* `sample` function
* extent stat ([min,max] pair)
* histgram
* vector field functions moved from clojure2d, api changed

### Changed

* interpolator have shorter names now (without `-interpolator` suffix)
* rbf functions are moved to separate namespace
* changed parameters order in interpolators (xs and ys are last)
* next-float-up -> next-double
* next-float-down -> prev-double

## [0.1.1]

### Added

* BesselJ
* all Gamma function variants (regularized, etc)
* regularized Beta
* log1p
* `estimation-strategy` can be passed to any function calculating quantiles/percentiles
* ->seq can accept number of samples now
* low-exp, high-exp to find lower/greater exponent for given base and number

### Changed

* `stats-map` doesn't contain select keys option anymore

## [0.1.0]

Initial version created from `Clojure2d` library.
