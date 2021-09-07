# Change Log
All notable changes to this project will be documented in this file. This change log follows the conventions of [keepachangelog.com](http://keepachangelog.com/).

## [2.1.6]

### Added

* `jinc` - besselj1(x)/x
* `muladd` macro
* `evalpoly` and `mevalpoly` - evaluate polynomial (from Julia)
* `makepoly` create polynomial function
* `Si` and `Ci` - https://dlmf.nist.gov/6.2#ii
* much more types of `skewness` and `kurtosis`
* much more effect size functions: `r2-determination`, `eta-sq`, `omega-sq`, `epsilon-sq`, `cohens-f2`, `cohens-q`, `cramers-v`, `cohens-w`, `tschuprows-t`
* `trim` and `winsor` data based on quantiles
* distances:
  - `dist-ang` - angular distance
  - `sim-cos` - cosine similarity
  - `angular` - refers to angular distance
* `->seq` in `fastmath.random` can accept sampling scheme (`:uniform`, `:stratified` and `:systematic`)
* new vector functions: `softmax`, `logsoftmax`, `logsumexp`, `logmeanexp`, `shift` (adds a value to all elements), `average` (mean / weighted average)

### Fixed

* `moment` didn't work properly for certain cases
* `effect-size` methods fixed
* [breaking] fixes around cosine distance and similarity

### Changed

* `fastmath.distance/cosine` - refers to cosine similarity now
* changes around data based distributions:
  - `continuous-distribution` is the same as `kde`
  - `integer-discrete-distribution` and `real-discrete-distribution` are sorted version of `enumerated-int` and `enumerated-real`

### Removed

* [breaking] `fastmath.vector/dist-cos` - wrong implementation

## [2.1.5]

### Fixed

* grid: `rhombus` is now a real rhombus, `triangle` is equilateral now (cont. #8)

## [2.1.4]

### Fixed

* `cell->mid` for `:triangle` grid returned wrong mid point for `down` triangles (#8)

### Changed

* [breaking] `:triangle` grid anchor is 3d now, last coordinate indicates triangle position (up/down)

## [2.1.3]

### Fixed

* better histogram split

## [2.1.2]

### Added

* bias corrected and accelerated percentile method for boostrap ([PR7](https://github.com/generateme/fastmath/pull/7))

## [2.1.1]

### Added

* bias corrected percentile method for boostrap ([PR6](https://github.com/generateme/fastmath/pull/6))

## [2.1.0]

### Added 

* `l-bfgs-b` optimizer as `:bfgs`

### Changed

* SMILE version bump to 2.6.0
* `fastmath.gp` relies now on SMILE backend (openblas/mkl)

## [2.0.5]

### Added

* `set-seed!` accept `:smile` RNG which calls `smile.math.MathEx/setSeed` function. Also, without RNG function will set seed to both: Smile and `default-rng`.

## [2.0.4]

### Added

* `:kde` distribution - create distribution based on kernel density estimation

## [2.0.3]

### Changed

* SMILE version bumped to 2.5.0
* netlib replaced by mkl (smile)

## [2.0.2]

### Added

* `moving-average-filter` - moving average smoothing
* `kernel-smoothing-filter` - signal smoothing filter using kernel methods

## [2.0.1]

### Added

* `savgol-filter` - Savitzky-Golay smoothing filter

## [2.0.0-alpha1]

Breaking change: due to significant change of SMILE API I decided to remove three ML namespaces: classification and regression. They probably be back as a separated bindings or incorporeted into SMILE directly.

### Added

* `lloyd` k-means variants
* `spectral` clustering
* `gp` namespace (GaussianProcesses) - moved from regression namespace

### Removed

* `neural-gas` clustering

### Changed

* [breaking] removed namespaces: classification and regression
* `:outliers?` key is removed from regrouped data for clusters, now `:outliers` key is created when outliers are present

## [1.5.3]

### Added

* `round-even` - even (or IEEE/IEC) rounding
* `cut` - cut range into even intervals
* `slice-range` - cut range into even steps
* `:dense` rank ties method

### Changed

* R2 sequence is able to generate up to dimensions=15
* [breaking] `rank` and `order` chagned to be zero based.
* Licence changed to MIT

### Fixed

* `co-intervals` works the same as in R now
* `group-by-intervals` checks all intervals instead the first one

## [1.5.2]

### Fixed

* `moment` when center is nil should use mean

## [1.5.1]

### Added

* `moment` function

### Removed

* `second-moment` (invalid implementation)

## [1.5.0]

### Changed

* removed slf4j deps 

## [1.5.0-alpha4]

### Changed

* some defaults for scan optimization

## [1.5.0-alpha3]

### Added

* B-Spline and polynomial interpolation
* `cb` as cube of the number

## [1.5.0-alpha2]

### Added

* `mnorm` - macro version of `norm`

### Changed

* protocol's wrapping functions shouldn't have primitive math hints (to avoid too many primitive->boxed conversions)

## [1.5.0-alpha1]

Cleaned documentation with more usage examples.

### Added

* `demean` function
* `seq->vec2`, `seq->vec3` and `seq->vec4`
* `warp-noise`
* fastmath.signal - audio signal processing with effects + oscillators

### Changed

* [possibly breaking] all protocols put in separate namespace, wrapping functions introduced (with type hints where possible)
* binary measure parameters

### Fixed

* acf calculation
* MAD SMILE bindings replaced (SMILE mutates an array)

### Removed

* `TOLERANCE` constant in `fastmath.vector`

## [1.4.0 SNAPSHOT]

### Added

* `optimization` package with various optimization methods
* BayesianOptimization
* NegativeBinomial distribution
* bunch of SSJ based distributions
* bootstrap datasets
* t-test

### Changed

* [breaking] `classification` and `regression` refactored
* [breaking] `kernel-density` moved to `kernel` namespace

### Fixed

* m/seq->double-array didn't recognize array type properly

## [1.3.0 SNAPSHOT]

### Added

* `make-vector` returns vector for given number of dimensions and optional sequence. 
* vectors implements now `IPersistentVector` to work with `vector?`
* core.matrix protocols for vectors
* math functions can operate on vectors now (like [[sin]] etc.)
* new functions for vectors `clamp`, `zero-count`, `nonzero-count`, `as-vec`
* classification bindings for SMILE, ~~XGBoost~~ and LIBLINEAR
* `expm1` function
* monotone interpolation
* various extent stats funcions
* `haversine` and `haversine-dist` (distance)
* predicates `nan?`, `inf?`, `valid-double?`, `invalid-double?`, `between?`
* calculate intervals for set of values `co-intervals` (same as R's function). Also `group-by-intervals`.
* various kernel-density methods and `kernel-density-ci` for kernel density with confidence interval
* various effect size functions
* binary measures like fn, tp, etc... (around 30 statistics)
* multivariate normal distribution added
* gaussian processes reimplemented
* kernels consolidated in one namespace

### Changed

* use `fmap` instead of `applyf` (now deprecated)
* outliers are samples which are outside inner fence instead of outer fence
* Vectors are associative now
* [breaking] `kernel-density` is now multimethod
* [breaking] `histogram` `:bins` contain starting point and sample counts only

### Removed

* rbf-obj converter
* rbf and mercer namespaces (use fastmath.kernel instead)

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
