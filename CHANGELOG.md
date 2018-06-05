# Change Log
All notable changes to this project will be documented in this file. This change log follows the conventions of [keepachangelog.com](http://keepachangelog.com/).

## [0.1.2-SNAPSHOT]

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
