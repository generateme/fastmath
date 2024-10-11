# Change Log
All notable changes to this project will be documented in this file. This change log follows the conventions of [keepachangelog.com](http://keepachangelog.com/).

## [3.0.0 alpha3 (unreleased)]

### Added

	* Vectors implement full IPersistantVector interface (to allow `conj`)

## [3.0.0 alpha2]

First release without SMILE dependency.
Also a lot of refactoring [breaking] was done.

### Added

* kringing interpolation is enhanced to work with various variograms
* semi-variograms for kringing
* rbf interpolation accepts polynomial term
* barycentric interpolation
* `:frequencies` and `:bins-maps` data in histogram result
* JTransform fft,dct,dst,dht transformations
* Data transformations (Power/Box-Cox/Yeo-Johnson)
* all *mean functions accept optional weights
* `skewness-test`, `kurtosis-test`, `normality-test`
* `durbin-watson` stat
* shifting and scaling matrix columns/rows
* `lm` and `glm`  in `fastmath.ml.regression` with full analysis
* `fastmath.ml.clustering` - (fuzzy) kmeans++ and dbscan clustering algorithms
* `dissimilarity` and `similarity` can accept data for histograms
* `bessel-j0`, `bessel-j1`, `bessel-j`, `bessel-i0`, `bessel-i1`, `bessel-i`, `bessel-y0`, `bessel-y1`, `bessel-y`, `bessel-k0`, `bessel-k1`, `bessel-k` add to `fastmath.special`
* `ince-C` and `ince-S` (with radial versions) polynomials

### Changed

Removed all functions based on SMILE library, most of them are replaced/rewritten.

* [breaking] `fastmath.clustering` namespace is moved to the separate [fastmath-clustering](https://github.com/generateme/fastmath-clustering) library.
* [breaking] `fastmath.gp` moved to the `fastmath.interpolation.gp`
* [breaking] `fastmath.kernel` is rewritten with changed API for `rbf` and `kernel`
* [breaking] `fastmath.interpolation` is rewritten with changed API, part of the methods are implemented directly in Clojure, names are unified
* [breaking] `fastmath.distance` rewritten
* `integrate-pdf` rewritten to use `:gauss-kronrod` integration (safer and faster)
* `fastmath.core` functionality based on functions with inline code not macros
* [breaking] polynomial funcitons moved to `fastmath.polynomials`
* [breaking] special funcitons (gamma, beta, erf, kummers, bessel, Si/Ci etc) moved to `fastmath.special`

### Fixed

* `binary-measures` accepts numerical sequence where `1` is treated as `true` and other values as `false` (#33)

## [2.4.0]

### Added

* Power divergence tests can work on histograms / distribution
* `dissimilarity` and `similarity` distance function of PDF (stats)

### Changed

* `SMILE` downgraded to 2.6.0 due to licence issue (3.0.2 is GPL, 2.6.0 is LGPL)
* [breaking] `correlation-matrix` can use only correlation methods (`:pearson`, `:kendall` and `:spearman`)
* some SMILE based methods are rewritten to Java or Clojure

### Deprecated

* `kullback-leibler-divergence` and `jensen-shannon-divergence`, use `dissimilarity` function instead.

### Fixed

* stats distances type inconsistency on some cases

## [2.3.0]

### Added

* Apache Commons Math `RealMatrix` and double[][] matrix interfaces
* `vec->array` converts any vector to double array
* Rotation matrix, quaternions, vectors consistency

### Changed

* `SMILE` deps updated

## [2.2.2-SNAPSHOT]

### Added

* `mat->float-array2d`, `mat->array`, `mat->float-array` - row-wise matrix representations as arrays
* `debug?` option for `:lbfgsb` optimizer
* `rint` accepts scale argument (same as `floor` and `ceil`)
* `roll-a-dice` to simulate dice rolling
* `double-exponent`, `double-significand` and `ulp` functions
* `fastmath.solver/quadratic` solver of quadratic equation ax^2+bx+c=0
* `fastmath.quaternion` namespace
* Cholesky decomposition for fixed size matrices
* Eigenvectors and eigenvalues for matrices
* `negative-zero?` tests for `-0.0` value
* `linear-optimizer` - linear optimization
* `mat`, `rows->mat` and `cols->mat` - additional matrix creators
* `norm` - various matrix norms (0, 1, Inf, max, generalized L_pq, Shatten p-norm, Frobenius)
* `condition` - condition number of the matrix
* `smooth-maximum` - boltzmann operator, lse, smu, p-norm, mellowmax
* `fastmath.calculus` - integration (VEGAS+ MC, h-adaptive Cubature and Gauss-Kronrod, ACM wrappers) and differentiations: derivative, gradient and hessian + Richardson extrapolation

### Changed

* `fastmath.complex` reimplemented

### Fixed

* `matrix` object name in printing, indentation
* speed of element matrix accessing is improved
* eastwood boxed math and reflection warnings (removed problematic macros, shadowing global vars, etc.)
* ensuring working with Clojure 1.10.3
* problem with field methods propagation in uberjar [#25](https://github.com/generateme/fastmath/pull/25), thanks @pfeodrippe

## [2.2.1]

### Added

* `:lbfgsb` java implementation of L-BFGS-B optimization algorithm, use it instead of `:bfgs` (buggy SMILE implementation)

## [2.2.0]

`2.1.9-SNAPSHOT` renamed to be `2.2.0` - see changes below

## [2.1.9-SNAPSHOT]

### Added

* `fma` - fused multiply and add `(+ (* x y) z)` (for java 9+)
* `fastmath.efloat` - floating point operations with error bounds
* `project` - vector projection
* Apache Commons Math `ArrayRealVector` extended with `VectorProto`
* `fastmath.matrix` - matrix 2x2, 3x3 and 4x4 basic operations (unrolled)
* fields refactored
* `binomial-ci` - confidence interval methods for binomial distribution
* new distributions (based on R package `gamlss.dist`):
    - `:zaga` - zero adjusted gamma distribution
	- `:nbi` - negative binomial type I
	- `:zinbi` - zero inflated negative binomial type I
	- `:zanbi` - zero altered negative binomial type I
	- `:zip` - poisson zero inflated
	- `:zip2` - zero inflated poisson 2
	- `:exgaus` - ex-Gaussian, exponentially modified gaussian
	- `:zabi` - zero adjusted binomial
	- `:zibi` - binomial zero inflated
	- `:bb` - beta binomial
	- `:zabb` - zero adjusted beta binomial
	- `:zibb` - zero inflated beta binomial
* other distributions
    - `:kolmogorov` - Kolmogorov distribution
	- `:half-normal` - half-normal distribution
	- `:truncated` and `:mixture` distributions
	- `:multinomial` distribution
* `hpdi-extent` and `pi-extent` - new extents based on HPDI and PI from `rethinking` R package
* `distribution-parameters` accept a keyword additionally
* `fastmath.solver/find-root` - root solver for univariate functions
* `delta-eq` and `delta=` for checking eqality with given accuracy.
* `coefficient-matrix` and `correlation-matrix`
* grid functions accept separate coordinates
* many additional `log-exp` functions
* weighted mean (`wmean`), median (`wmedian`) and quantiles (`wquantile`,`wquantiles`)
* Vec2/3/4 implement IReduce and IReduceInit
* historical trigonometric functions (versine with variants, chord, exsec/excsc) with their inverses
* Clerk notebook for: fastmath.core
* Added PDF to `:cramer-von-misses` distribution (based on finite difference)
* `integrate-pdf` converted from private to public
* `contingency-2x2-measures` - selection of 2x2 contingency measures
* clj-kondo hooks
* `fastmath.stats.bootstrap` introduced as the replacement of `fastmath.stats/bootstrap` function (now deprecated)

#### Tests

* `p-value` - helper for calculating p-value
* `contingency-table` - calculate frequencies
* `binomial-test`
* `t-test-one-sample`, `t-test-two-samples` (more information returned)
* `z-test-one-sample`, `z-test-two-samples`
* `f-test`
* `power-divergence-test` + special cases
    - `chisq-test`
	- `multinomial-likelihood-ratio-test`
	- `minimum-discrimination-information-test`
	- `neyman-modified-chisq-test`
	- `freeman-tukey-test`
	- `cressie-read-test`
* distribution tests
    - `ad-test-one-sample` - Anderson Darling test
	- `ks-test-one-sample`, `ks-test-two-samples` - Kolmogorov Smirnov test
* ANOVA tests
	- `one-way-anova-test`
	- `levene-test`
	- `brown-forsythe-test`
	- `fligner-killeen-test`
* `kruskal-test` - Kruskal-Wallis rank sum test

### Changed

* `(set! *warn-on-reflection* true)` removed ([motivation](https://github.com/clojure-emacs/refactor-nrepl/issues/347))
* names changed: `ttest-` -> `t-test-`
* RidgedMulti.java refactored
* effect size functions revisited

### Fixed

* a nasty bug with primitive macro generation, how could I missed that?
* `:histogram` issues with low number of samples
* `fast-max` and `fast-min` had wrong inline operation
* wrong primitive hinting and protocol extensions
* `set-seed` calls protocol for PRNGs
* `integrate-pdf` accuracy increased
* `nrd` fixed when bandwidth was 0
* kde estimator revisited

## [2.1.8]

### Added

* `use-primitive-operators` accepts optional set of symbols which shouldn't be imported

### Fixed

* `abs` redefinition in `vector` was not properly implemented (problem was visible when `vector` was precompiled with javac for Clojure2d)

## [2.1.7 - do not use]

### Added

* `ball-random` - unit ball random vector
* `sequence-generator` `:ball` added - unit ball random sampling
* some more fields
* `ceil` and `floor` can snap to the nearest multiply of scale parameter (optional)

### Fixed

* Vec3 rotation
* [breaking] ensure proper behaviour for Clojure interfaces for fastmath.vector custom types (see: https://github.com/nextjournal/clerk/issues/64)
* `next-double` and `prev-double` to properly cross 0.0
* Clojure 1.11.0 fixes (`abs`)

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
