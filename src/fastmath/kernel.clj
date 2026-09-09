(ns fastmath.kernel
  "Various kernel functions.

  * RBF (double -> double functions)
  * vector kernels (vector x vector -> double function; may be positive definite, conditional positive definite, positive semi-definite, mercer)
  * density estimation
  * some kernel operations"
  (:require [fastmath.core :as m]
            [fastmath.kernel.rbf :as rbf]
            [fastmath.kernel.vector :as vk]
            [fastmath.kernel.density :as dens]
            [fastmath.kernel.window :as win]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)
(m/use-primitive-operators)

(defmacro ^:private emit
  ([kind nm f] `(emit ~kind ~nm ~f nil))
  ([kind nm f params]
   (if params
     `(defmethod ~kind ~nm
        ([_#] (~f ~params))
        ([_# params#] (~f (merge params# ~params))))
     `(defmethod ~kind ~nm
        ([_#] (~f))
        ([_# params#] (~f params#))))))


(defmulti rbf
  "Creates a radial basis function (RBF) kernel, a `double->double` function of the radius (distance from a center).

  RBFs are commonly used as the basis for scattered-data interpolation and approximation (see `fastmath.kernel.rbf` for the underlying implementations and formulas).

  Parameters:

  - `k` - a keyword selecting the kernel family or one of its parameter presets (dispatch value).
  - `params` (optional map) - kernel-specific parameters, always including `:scale` (default `1.0`), which rescales the input radius before evaluating the kernel (i.e. the kernel receives `x/scale`). When called with no `params` (or with `nil`), default parameter values are used, where defined.

  Available `k` values, grouped by family (each preset fixes some of the family's parameters, but `:scale` still applies):

  - `:linear`, `:gaussian` - no extra parameters besides `:scale`.
  - `:truncated-power` - needs `:k` exponent; presets `:truncated-power-1`, `:truncated-power-2`, `:truncated-power-3`, `:truncated-power-half`, `:truncated-power-third` fix `:k` to `1`, `2`, `3`, `0.5` and `1/3` respectively.
  - `:gaussians-laguerre` - needs `:dimension` and `:degree`; presets `:gaussians-laguerre-11`, `-12`, `-21`, `-22` fix `dimension`/`degree` to `1`/`1`, `1`/`2`, `2`/`1`, `2`/`2`.
  - `:poisson` - needs `:d`; presets `:poisson-2`, `:poisson-3`, `:poisson-4` fix `:d` to `2`, `3`, `4`.
  - `:matern` - needs odd `:order` (i.e. `order/2` such as `1/2`, `3/2`, `5/2`); presets `:matern-c0`, `:matern-c2`, `:matern-c4` fix `:order` to `1`, `3`, `5`.
  - `:generalized-multiquadratic` - needs `:beta` exponent and optional `:negate?` (`true`, `false` or `:auto`, default `:auto`); presets `:multiquadratic` and `:inverse-multiquadratic` fix `:beta` to `0.5` and `-0.5` with `:negate?` `false`.
  - `:radial-powers` - needs `:beta` and optional `:negate?`; preset `:radial-powers-3` fixes `:beta` to `3.0` with `:negate?` `false`.
  - `:thin-plate-splines` - needs `:beta` and optional `:negate?`; preset `:thin-plate` fixes `:beta` to `1.0` with `:negate?` `false`.
  - `:shifted-surface-splines` - needs `:s` (long) and `:beta` (double, greater than half of `s`).
  - `:wendland` - needs `:s` (double) and `:k` (long, only values less than `5` are supported).
  - `:gneiting` - needs `:s` and `:l`.
  - `:wu` - needs `:l` (double) and `:k` (long); presets `:wu-10` through `:wu-33` fix `:l` to `1`, `2` or `3` and `:k` from `0` up to `l`. Combinations outside the presets fall back to numerical integration/differentiation and may be less stable.
  - `:whittaker` - needs `:alpha`, `:k` and `:beta`; closed-form expressions exist for `alpha` `0` or `1` combined with `k` `2` or `3`, other combinations use numerical integration.

  Returns a function of a single `double` (the radius) returning a `double`. Throws when `k` is not a registered kernel name.

  See also [[kernel-list]] for the full, up-to-date list of registered RBF names, [[kernel]] for vector (non-radial) kernels."
  (fn [k & _] k))

(defmacro ^:private emit-rbf
  ([nm f] `(emit ~'rbf ~nm ~f))
  ([nm f params] `(emit ~'rbf ~nm ~f ~params)))

(emit-rbf :linear rbf/linear)
(emit-rbf :gaussian rbf/gaussian)

(emit-rbf :truncated-power rbf/truncated-power)
(emit-rbf :truncated-power-1 rbf/truncated-power {:k 1.0})
(emit-rbf :truncated-power-2 rbf/truncated-power {:k 2.0})
(emit-rbf :truncated-power-3 rbf/truncated-power {:k [3.0]})
(emit-rbf :truncated-power-half rbf/truncated-power {:k 0.5})
(emit-rbf :truncated-power-third rbf/truncated-power {:k m/THIRD})

(emit-rbf :gaussians-laguerre rbf/gaussians-laguerre)
(emit-rbf :gaussians-laguerre-11 rbf/gaussians-laguerre {:dimension 1.0 :degree 1.0})
(emit-rbf :gaussians-laguerre-12 rbf/gaussians-laguerre {:dimension 1.0 :degree 2.0})
(emit-rbf :gaussians-laguerre-21 rbf/gaussians-laguerre {:dimension 2.0 :degree 1.0})
(emit-rbf :gaussians-laguerre-22 rbf/gaussians-laguerre {:dimension 2.0 :degree 2.0})

(emit-rbf :poisson rbf/poisson)
(emit-rbf :poisson-2 rbf/poisson {:d 2.0})
(emit-rbf :poisson-3 rbf/poisson {:d 3.0})
(emit-rbf :poisson-4 rbf/poisson {:d 4.0})

(emit-rbf :matern rbf/matern)
(emit-rbf :matern-c0 rbf/matern {:order 1.0})
(emit-rbf :matern-c2 rbf/matern {:order 3.0})
(emit-rbf :matern-c4 rbf/matern {:order 5.0})

(emit-rbf :generalized-multiquadratic rbf/generalized-multiquadratic)
(emit-rbf :multiquadratic rbf/generalized-multiquadratic {:beta 0.5 :negate? false})
(emit-rbf :inverse-multiquadratic rbf/generalized-multiquadratic {:beta -0.5 :negate? false})

(emit-rbf :radial-powers rbf/radial-powers)
(emit-rbf :radial-powers-3 rbf/radial-powers {:beta 3.0 :negate? false})

(emit-rbf :thin-plate-splines rbf/thin-plate-splines)
(emit-rbf :thin-plate rbf/thin-plate-splines {:beta 1.0 :negate? false})

(emit-rbf :shifted-surface-splines rbf/shifted-surface-splines)

(emit-rbf :wendland rbf/wendland)
(emit-rbf :gneiting rbf/gneiting)

(emit-rbf :wu rbf/wu)
(emit-rbf :wu-10 rbf/wu {:l 1.0 :k 0.0})
(emit-rbf :wu-11 rbf/wu {:l 1.0 :k 1.0})
(emit-rbf :wu-20 rbf/wu {:l 2.0 :k 0.0})
(emit-rbf :wu-21 rbf/wu {:l 2.0 :k 1.0})
(emit-rbf :wu-22 rbf/wu {:l 2.0 :k 2.0})
(emit-rbf :wu-30 rbf/wu {:l 3.0 :k 0.0})
(emit-rbf :wu-31 rbf/wu {:l 3.0 :k 1.0})
(emit-rbf :wu-32 rbf/wu {:l 3.0 :k 2.0})
(emit-rbf :wu-33 rbf/wu {:l 3.0 :k 3.0})

(emit-rbf :whittaker rbf/whittaker)

;;

(defmulti kernel
  "Creates a vector kernel, a `(x y) -> double` similarity/distance-based function for two vectors (or two numbers).

  Vector kernels are used in kernel methods (SVMs, Gaussian processes, kriging, kernel regression, etc.) to measure similarity between data points; depending on the family they may be Mercer, positive definite, conditional positive definite, positive semi-definite kernels or none of the above.

  Parameters:

  - `k` - a keyword selecting the kernel family or one of its parameter presets (dispatch value).
  - `params` (optional map) - kernel-specific parameters (see below); most distance-based kernels accept a `:distance` key (a `(x y) -> double` function, default euclidean, see `fastmath.vector`) used to measure the separation between `x` and `y`. When called with no `params` (or with `nil`), default parameter values are used.

  Available `k` values, grouped by family:

  - `:linear` - plain dot product, no parameters.
  - `:polynomial` - power of the (shifted, scaled) dot product; params `:alpha` (default `1.0`), `:c` (default `0.0`), `:p` (default `2.0`).
  - `:hyperbolic-tangent` - hyperbolic tangent of the (shifted, scaled) dot product; params `:alpha` (default `1.0`), `:c` (default `0.0`).
  - `:gaussian`, `:exponential`, `:laplacian`, `:wave`, `:cauchy` - exponential-family kernels of the distance; param `:sigma` (default `1.0`) plus `:distance`.
  - `:periodic` - exponential of a squared sine of the distance; params `:sigma`, `:periodicity` (both default `1.0`) plus `:distance`.
  - `:power`, `:log` - negative power (resp. negative log of one plus power) of the distance; param `:p` (default `2.0`) plus `:distance`.
  - `:rational-quadratic`, `:multiquadratic`, `:inverse-multiquadratic` - functions of the squared distance shifted by `:c` (default `1.0`) plus `:distance`.
  - `:generalized-t-student` - reciprocal of one plus a power of the distance; params `:p`, `:distance`.
  - `:pearson` - Pearson VII kernel; params `:sigma`, `:omega` (both default `1.0`) plus `:distance`.
  - `:hyperbolic-secant` - hyperbolic secant of a scaled distance; param `:a` (default `1.0`) plus `:distance`.
  - `:bessel` - Bessel function of the first kind of the (scaled) distance; params `:sigma`, `:n`, `:v` (defaults `1.0`, `2.0`, `-1.0`) plus `:distance`.
  - `:bessel2` - alternative (R kernlab) Bessel kernel; params `:sigma`, `:degree`, `:order` (defaults `1.0`, `1.0`, `0.0`) plus `:distance`.
  - `:matern` - Matern kernel; params `:order` (odd, default `1`, use `5` for Matern 5/2 etc.), `:theta` (default `1.0`) plus `:distance`; presets `:matern-12`, `:matern-32`, `:matern-52` fix `:order` to `1`, `3`, `5`.
  - `:geometric` - compactly supported kernel; params `:n` (dimension), `:r` (default `1.0`) plus `:distance`; presets `:triangular`, `:circular`, `:spherical` fix `:n` to `1`, `2`, `3`.
  - `:anova` - sum, over elementwise powers of the vectors, of a Gaussian-like term; params `:sigma`, `:k`, `:d` (all default `1.0`).
  - `:spline` - product, over vector elements, of a fixed cubic-spline expression, no parameters.
  - `:b-spline` - product, over vector elements, of a B-spline basis of degree `:n` (default `2.0`).
  - `:chi-square`, `:chi-square2`, `:histogram` - elementwise, sum-based kernels for non-negative vectors (e.g. histograms), no parameters.
  - `:generalized-histogram` - elementwise minimum of powered absolute values; param `:p` (default `2.0`).
  - `:dirichlet` - Dirichlet kernel; param `:n` dimensionality (default `1.0`).
  - `:hellinger` - Hellinger kernel (sum of the square roots of elementwise products), no parameters.

  Returns a function of two arguments `x` and `y` (vectors or numbers, matching `fastmath.vector` conventions) returning a `double`. Throws when `k` is not a registered kernel name.

  See also [[kernel-list]] for the full, up-to-date list of registered kernel names, [[rbf]] for radial (single-argument) kernels, [[exp]], [[scale]], [[mult]], [[wmean]] for combining kernels."
  (fn [k & _] k))

(defmacro ^:private emit-kernel
  ([nm f] `(emit ~'kernel ~nm ~f))
  ([nm f params] `(emit ~'kernel ~nm ~f ~params)))

(emit-kernel :linear vk/linear)
(emit-kernel :polynomial vk/polynomial)
(emit-kernel :gaussian vk/gaussian)
(emit-kernel :exponential vk/exponential)
(emit-kernel :laplacian vk/laplacian)
(emit-kernel :anova vk/anova)

(emit-kernel :hyperbolic-tangent vk/hyperbolic-tangent)
(emit-kernel :hyperbolic-secant vk/hyperbolic-secant)

(emit-kernel :rational-quadratic vk/rational-quadratic)
(emit-kernel :multiquadratic vk/multiquadratic)
(emit-kernel :inverse-multiquadratic vk/inverse-multiquadratic)

(emit-kernel :triangular vk/geometric {:n 1})
(emit-kernel :circular vk/geometric {:n 2})
(emit-kernel :spherical vk/geometric {:n 3})
(emit-kernel :geometric vk/geometric)

(emit-kernel :wave vk/wave)
(emit-kernel :periodic vk/periodic)

(emit-kernel :power vk/power)
(emit-kernel :log vk/log)

(emit-kernel :spline vk/spline)
(emit-kernel :b-spline vk/b-spline)

(emit-kernel :bessel vk/bessel)
(emit-kernel :bessel2 vk/bessel2)

(emit-kernel :cauchy vk/cauchy)

(emit-kernel :chi-square vk/chi-square)
(emit-kernel :chi-square2 vk/chi-square2)

(emit-kernel :histogram vk/histogram)
(emit-kernel :generalized-histogram vk/generalized-histogram)
(emit-kernel :generalized-t-student vk/generalized-t-student)

(emit-kernel :dirichlet vk/dirichlet)
(emit-kernel :hellinger vk/hellinger)
(emit-kernel :pearson vk/pearson)

(emit-kernel :matern vk/matern)
(emit-kernel :matern-12 vk/matern {:order 1})
(emit-kernel :matern-32 vk/matern {:order 3})
(emit-kernel :matern-52 vk/matern {:order 5})

(defn exp
  "Kernel wraper. exp of kernel `k` with optional scaling value `t`."
  ([k] (exp k 1.0))
  ([k ^double t]
   (fn [x y] (m/exp (* t ^double (k x y))))))

(defn scale
  "Kernel wrapper. Scale kernel result."
  [k ^double scale]
  (fn [x y] (* scale ^double (k x y))))

(defn mult
  "Kernel wrapper. Multiply two or more kernels."
  ([k1] k1)
  ([k1 k2] (fn [x y] (* ^double (k1 x y) ^double (k2 x y))))
  ([k1 k2 k3] (fn [x y] (* ^double (k1 x y) ^double (k2 x y) ^double (k3 x y))))
  ([k1 k2 k3 & r]
   (let [k (mult k1 k2 k3)]
     (if-not (seq r) k
             (apply mult k r)))))

(defn wmean
  "Kernel wrapper. (Weighted) mean of kernel results."
  ([kernels] (wmean kernels (repeat (count kernels) 1.0)))
  ([kernels weights]
   (fn [x y] (wmean (map (fn [k] (k x y)) kernels) weights))))

;; kernel density estimation

(defn bandwidth
  "Estimates the 1d kernel density bandwidth (`h`) for `data`.

  The bandwidth controls how much each data point is spread out by the kernel; smaller values follow the data more closely (risking overfitting/noise), larger values produce a smoother, more biased estimate.

  Parameters:

  - `kernel` - a keyword naming a kernel registered in [[kernel-list]] under `:kde`, used by the `:nrd-adjust` method and by the cross-validation targets (`:rlcv`, `:lcv`, `:lscv`).
  - `data` - a sequence of numbers to estimate the bandwidth from.
  - `h` - selects the estimation method, a keyword, one of:
      - `:nrd` - rule-of-thumb, scaled by `1.06` (Silverman's rule).
      - `:nrd0` - rule-of-thumb, scaled by `0.9`.
      - `:nrd-adjust` - `:nrd` further adjusted by a kernel-specific canonical bandwidth factor; not supported for `silverman` and `cauchy` kernels.
      - `:rlcv` - robust likelihood cross-validation.
      - `:lcv` - likelihood cross-validation.
      - `:lscv` - least squares cross-validation.

  Returns the estimated bandwidth as a `double`.

  See also [[kernel-density]], [[kernel-density-ci]] which accept the same set of keywords (or a plain number) as a bandwidth."
  [kernel data h] (dens/bandwidth kernel data h))

(defn kernel-density
  "Returns a 1d kernel density estimation (KDE) for `data`.

  A kernel density estimate is a smoothed histogram: it places a scaled copy of `kernel` at every data point and sums the contributions, weighted by a bandwidth that controls how much each point spreads out.

  Parameters:

  - `kernel` - a keyword naming a kernel registered in [[kernel-list]] under `:kde`, or a custom kernel function accepting and returning a `double`.
  - `data` - a sequence of numbers to estimate the density from.
  - `bandwidth` (optional) - either a plain number (used directly as `h`), or a map accepted as `params`:
      - `:bandwidth` - the bandwidth `h`, a number or one of the keywords accepted by [[bandwidth]] (`:nrd`, `:nrd0`, `:nrd-adjust`, `:rlcv`, `:lcv`, `:lscv`). Default: `:nrd`.
      - `:binned?` - whether `data` should be pre-binned before evaluation to speed up estimation on large datasets. When `true`, bin width is `h` divided by `5`; when a number, that number is used as the divisor instead. Default: `false`.
  - `info?` (optional, default `false`) - when `true`, returns the estimator together with its fitting details instead of the bare function.

  Returns a function of a single `double` returning the estimated density at that point as a `double`. When `info?` is `true`, returns a map instead:

  - `:kde` - the density function described above.
  - `:factor` - the normalizing factor, equal to `1/(n*h)`.
  - `:h` - the bandwidth actually used.
  - `:mn` and `:mx` - the inferred support extent of the estimator (data range expanded by the kernel radius).

  See also [[bandwidth]] (bandwidth estimation alone), [[kernel-density-ci]] (density with confidence intervals)."
  ([kernel data] (kernel-density kernel data {:bandwidth :nrd}))
  ([kernel data bandwidth] (kernel-density kernel data bandwidth false))
  ([kernel data bandwidth info?]
   (let [h (if (number? bandwidth) {:bandwidth bandwidth} bandwidth)
         f (if info? dens/kernel-density+ dens/kernel-density)]
     (f kernel data h))))

(defn kernel-density-ci
  "Returns a 1d KDE for `data` that also reports pointwise asymptotic confidence intervals.

  Uses the standard asymptotic normal approximation of the kde variance (see section 6.1.5 of http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/tutorials/xlghtmlnode33.html): for a density estimate `f(x)` obtained with normalizing factor `1/(n*h)`, the variance is approximated as `k2 * f(x) / (n*h)`, where `k2` is the kernel roughness (integral of the squared kernel), and confidence bounds are derived from the normal quantile at the requested `alpha`.

  Parameters:

  - `kernel` - a keyword naming a kernel registered in [[kernel-list]] under `:kde` (a custom kernel function is not accepted, since its roughness constant `k2` must be looked up).
  - `data` - a sequence of numbers to estimate the density from.
  - `bandwidth-or-params` (optional) - either a plain number (used directly as bandwidth `h`), or a map with the same keys as [[kernel-density]] plus:
      - `:alpha` - the confidence level parameter, the resulting interval has coverage `1 - alpha`. Default: `0.05`.

  Returns a function of a single `double` `x` returning a 3-element vector `[fx lower upper]`, where `fx` is the estimated density at `x`, and `lower`/`upper` are the corresponding confidence bounds.

  Throws an assertion error when `kernel` is not one of the kernels known to [[kernel-list]] `:kde`.

  See also [[kernel-density]], [[bandwidth]]."
  ([kernel data] (kernel-density-ci kernel data {:bandwidth :nrd}))
  ([kernel data bandwidth-or-params]
   (let [p (if (number? bandwidth-or-params) {:bandwidth bandwidth-or-params} bandwidth-or-params)]
     (dens/kernel-density-ci kernel data p))))

;; window

(def ^:private windows
  {:rectangular05 win/rectangular05
   :rectangular win/rectangular
   :triangular win/triangular
   :parzen win/parzen
   :b-spline win/b-spline
   :welch win/welch
   :connes win/connes
   :parzen-algebraic win/parzen-algebraic
   :singla-singh win/singla-singh
   :sinc win/sinc
   :fejer win/fejer
   :de-la-vallee-poussin win/de-la-vallee-poussin
   :lanczos win/lanczos
   :hamming win/hamming
   :hamming-exact win/hamming-exact
   :hann win/hann
   :raised-cosine win/raised-cosine
   :webster-hamming win/webster-hamming
   :power-of-cosine win/power-of-cosine
   :raised-power-of-cosine win/raised-power-of-cosine
   :parzen-cosine win/parzen-cosine
   :bohman win/bohman
   :trapezoid win/trapezoid
   :tukey win/tukey
   :bartlett-hann win/bartlett-hann
   :blackman-harris-family win/blackman-harris-family
   :blackman win/blackman
   :blackman-exact win/blackman-exact
   :blackman-harris win/blackman-harris
   :blackman-harris-61db win/blackman-harris-61db
   :blackman-harris-67db win/blackman-harris-67db
   :blackman-harris-74db win/blackman-harris-74db
   :blackman-harris-92db win/blackman-harris-92db
   :nutall-3-1st win/nutall-3-1st
   :nutall-3-3rd win/nutall-3-3rd
   :blackman-nutall win/blackman-nutall
   :nutall-1st win/nutall-1st
   :nutall-3rd win/nutall-3rd
   :nutall-5th win/nutall-5th
   :mottaghi-kashtiban-shayesteh win/mottaghi-kashtiban-shayesteh
   :low-sidelobe win/low-sidelobe
   :exponential win/exponential
   :hanning-poisson win/hanning-poisson
   :gaussian win/gaussian
   :parzen-exponential win/parzen-exponential
   :dolph-chebyshev win/dolph-chebyshev
   :taylor win/taylor
   :cauchy win/cauchy
   :parzen-geometric win/parzen-geometric
   :kaiser-bessel win/kaiser-bessel
   :cosh win/cosh
   :avci-nacaroglu win/avci-nacaroglu
   :knab win/knab
   :ultraspherical win/ultraspherical
   :saramaki win/saramaki
   :legendre win/legendre
   :bessel-I1 win/bessel-I1
   :shayesteh-kashtiban win/shayesteh-kashtiban
   :kaiser-bessel-derived win/kaiser-bessel-derived
   :vorbis win/vorbis
   :flat-top win/flat-top
   :flat-top-3 win/flat-top-3})

(defn window
  "Returns tapering window coefficients, `N` numbers shaped to fade toward the edges.

  Windows (also called taper functions) are used to reduce spectral leakage before a Fourier transform, to design FIR filters, or to smooth signal edges.

  Parameters:

  - `window-name` - a keyword naming one of the windows registered in [[kernel-list]] under `:window`, or a custom function `(fn [N options] ...)` returning a sequence of `N` coefficients.
  - `N` - the number of samples to generate (default: `256`).
  - `options` (optional map), common to every window:
      - `:symmetric?` - `true` (default) for a symmetric window (suitable for FIR filter design), `false` for a periodic window (suitable for spectral analysis / FFT-based methods). A periodic window of size `N` is obtained by sampling a symmetric window of size `N+1` and dropping the last coefficient.
      - `:normalize?` - controls how a continuous window is rescaled after sampling. `true` (the default for most windows) rescales it to a maximum value of `1.0`; `false` leaves the raw sampled values (whose integral over the window is `1.0`); other accepted values are `:L1`, `:L2`, `:LInf`, `:N`, or a plain number to divide by (see `fastmath.kernel.window/normalize-coefficients`). Exceptions to the `true` default are `:rectangular05` and `:kaiser-bessel-derived`, which default to `false`.

  Some windows accept additional `options` keys beyond `:symmetric?`/`:normalize?`:

  - `:triangular` - `:shift` (long, default `2`).
  - `:parzen` - `:discrete?` (default `true`).
  - `:b-spline` - `:order` (long, default `3`).
  - `:connes`, `:exponential`, `:hanning-poisson`, `:gaussian`, `:cauchy`, `:kaiser-bessel`, `:kaiser-bessel-derived` - `:alpha` shape parameter (default `1.0`).
  - `:parzen-algebraic` - `:gamma` (default `1.0`), `:u` (default `3.0`).
  - `:singla-singh` - `:order` (double, default `1.0`).
  - `:lanczos` - `:L` (double, default `3.0`).
  - `:raised-cosine` - `:alpha` (default `0.5`).
  - `:webster-hamming` - `:v` (default `1.0`).
  - `:power-of-cosine` - `:m` (default `1.0`).
  - `:raised-power-of-cosine` - `:alpha` (default `0.05`), `:m` (default `1.0`).
  - `:parzen-cosine` - `:gamma` (default `1.0`), `:m` (default `2.0`).
  - `:trapezoid` - `:alpha` (default `0.25`).
  - `:tukey` - `:alpha` (default `0.5`).
  - `:blackman-harris-family` - `:coeffs` vector of coefficients (default `[0.5 0.5]`).
  - `:parzen-exponential`, `:parzen-geometric` - `:alpha` (default `1.0`), `:r` (default `1.0`).
  - `:dolph-chebyshev`, `:saramaki`, `:legendre` - `:level` target sidelobe attenuation in dB (default `-50.0`).
  - `:taylor` - `:level` (default `-50.0`), `:n` number of nearly constant-level sidelobes, a long (default `4`).
  - `:cosh`, `:avci-nacaroglu`, `:knab`, `:bessel-I1` - `:alpha` shape parameter (default `2.0`).
  - `:ultraspherical` - `:level` (default `-50.0`), `:alpha` (default `2.0`).

  The remaining registered windows (e.g. `:rectangular`, `:welch`, `:sinc`, `:fejer`, `:hamming`, `:hann`, `:bohman`, `:blackman` and its variants, the `:nutall-*`/`:blackman-*` family, `:vorbis`, `:flat-top`, `:flat-top-3`, and others) only accept the common `:symmetric?`/`:normalize?` options.

  Returns a sequence of `N` `double` window coefficients.

  Throws `ex-info` when `window-name` is neither a registered keyword nor a function.

  See also [[kernel-list]] for the full, up-to-date list of registered window names."
  ([window-name] (window window-name 256))
  ([window-name ^long N] (window window-name N nil))
  ([window-name ^long N {:keys [symmetric?]
                         :or {symmetric? true}
                         :as options}]
   (if-let [wind (or (windows window-name)
                     (when (fn? window-name) window-name))]
     (if symmetric? (wind N options) (butlast (wind (m/inc N) options)))
     (throw (ex-info "Unknown window." {:window-name window-name})))))

;;

(def kernel-list ^{:doc "List of available kernels: vector, rbf, kde, windows (tapering)"}
  {:vector (sort (keys (methods kernel)))
   :rbf (sort (keys (methods rbf)))
   :kde (sort (keys dens/kde-data))
   :window (sort (keys windows))})

;; variograms



(m/unuse-primitive-operators)
