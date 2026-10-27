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

(defn- same-value?
  "True when a preset value and a user value are the same: numbers are compared by value (`3` equals `3.0`), the rest with `=`."
  [a b]
  (if (and (number? a) (number? b))
    (m/== (double a) (double b))
    (= a b)))

(defn- resolve-params
  "Validates the parameters `params` given to the registered kernel `nm` of the registry `kind` (`:rbf` or `:kernel`) and merges them with the `preset`.

  Returns the preset (possibly `nil`) for `nil` params, otherwise the merged map. Throws `ex-info` (`:type :invalid-input`) when `params` is neither `nil` nor a map, or when it holds a value which differs from the value fixed by the preset."
  [kind nm params preset]
  (cond
    (nil? params) preset
    (map? params) (do (doseq [[pk pv] preset
                              :let [uv (get params pk)]
                              :when (some? uv)]
                        (when-not (same-value? pv uv)
                          (throw (ex-info (str "(" (name kind) " " nm ") fixes " pk " to " (pr-str pv)
                                               ", conflicting value given: " (pr-str uv))
                                          {:type :invalid-input
                                           :kind kind
                                           :name nm
                                           :option pk
                                           :value uv
                                           :expected pv}))))
                      (merge params preset))
    :else (throw (ex-info (str "(" (name kind) " " nm ") expects a map of parameters or nil, got: " (pr-str params))
                          {:type :invalid-input
                           :kind kind
                           :name nm
                           :value params
                           :expected "a map or nil"}))))

(defmacro ^:private emit
  ([kind nm f] `(emit ~kind ~nm ~f nil))
  ([kind nm f preset]
   (let [kw (keyword (name kind))]
     (if preset
       `(defmethod ~kind ~nm
          ([_#] (~f ~preset))
          ([_# params#] (~f (resolve-params ~kw ~nm params# ~preset))))
       `(defmethod ~kind ~nm
          ([_#] (~f))
          ([_# params#]
           (let [p# (resolve-params ~kw ~nm params# nil)]
             (if (nil? p#) (~f) (~f p#)))))))))


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
  - `:poisson` - needs `:d` (at least `1`), the kernel is normalised to `1` at `0`; presets `:poisson-2`, `:poisson-3`, `:poisson-4` fix `:d` to `2`, `3`, `4`.
  - `:matern` - needs odd `:order` (i.e. `order/2` such as `1/2`, `3/2`, `5/2`); presets `:matern-c0`, `:matern-c2`, `:matern-c4` fix `:order` to `1`, `3`, `5`.
  - `:generalized-multiquadratic` - needs `:beta` exponent and optional `:negate?` (`true`, `false` or `:auto`, default `:auto`: the sign is `(-1)^ceil(beta)` for positive `beta`, `+1` otherwise); presets `:multiquadratic` and `:inverse-multiquadratic` fix `:beta` to `0.5` and `-0.5` with `:negate?` `false`.
  - `:radial-powers` - needs `:beta` and optional `:negate?`; preset `:radial-powers-3` fixes `:beta` to `3.0` with `:negate?` `false`.
  - `:thin-plate-splines` - needs `:beta` (greater than `0`) and optional `:negate?` (`:auto` gives `(-1)^(ceil(beta) + 1)`, valid for integer `beta`); preset `:thin-plate` fixes `:beta` to `1.0` with `:negate?` `false`.
  - `:shifted-surface-splines` - needs `:s` (positive integer) and `:beta` (greater than half of `s`, an integer for even `s`).
  - `:wendland` - needs `:s` (dimension, at least `1`) and `:k` (integer from `0` to `4`); the degree is `floor(s/2) + k + 1`.
  - `:gneiting` - needs `:s` and `:l`.
  - `:wu` - needs `:l` (at least `0`) and `:k` (integer, at least `0`); presets `:wu-10` through `:wu-33` fix `:l` to `1`, `2` or `3` and `:k` from `0` up to `l`. Combinations outside the presets fall back to numerical integration/differentiation and may be less stable.
  - `:whittaker` - needs `:alpha`, `:k` and `:beta`; closed-form expressions exist for `alpha` `0` or `1` combined with `k` `2` or `3`, other combinations use numerical integration.

  Returns a function of a single `double` (the radius) returning a `double`.

  Corner cases and exceptions:

  - `params` must be a map or `nil`. A number (for example `(rbf :gaussian 2.0)`) or any other value throws `ex-info` with `:type :invalid-input`, the message names the keyword and the received value. Use `{:scale 2.0}` instead.
  - A preset keyword (for example `:matern-c2`) fixes some parameters. A `params` map holding a conflicting value for a fixed key (`(rbf :matern-c2 {:order 5})`) throws `ex-info` (`:type :invalid-input`, with `:option`, `:value` and `:expected`); the same value (`{:order 3}`, `{:order 3.0}`) or other keys are accepted.
  - A missing required parameter throws `ex-info` from the kernel constructor in `fastmath.kernel.rbf`, naming the key.
  - Throws `IllegalArgumentException` when `k` is not a registered kernel name.

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
(emit-rbf :truncated-power-3 rbf/truncated-power {:k 3.0})
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
  - `:power`, `:log` - negative power (resp. negative log of one plus power) of the distance; param `:p` (default `2.0`) plus `:distance`. Conditionally positive definite for `0 < p <= 2`.
  - `:rational-quadratic`, `:multiquadratic`, `:inverse-multiquadratic` - functions of the squared distance shifted by `:c` (default `1.0`) plus `:distance`.
  - `:generalized-t-student` - reciprocal of one plus a power of the distance; params `:p` (default `2.0`), `:distance`.
  - `:pearson` - Pearson VII kernel; params `:sigma`, `:omega` (both default `1.0`) plus `:distance`.
  - `:hyperbolic-secant` - hyperbolic secant of a scaled distance; param `:a` (default `1.0`) plus `:distance`.
  - `:bessel` - Bessel function of the first kind of the (scaled) distance; params `:sigma`, `:n`, `:v` (defaults `1.0`, `2.0`, `-1.0`) plus `:distance`.
  - `:bessel2` - alternative (R kernlab) Bessel kernel; params `:sigma`, `:degree`, `:order` (defaults `1.0`, `1.0`, `0.0`) plus `:distance`.
  - `:matern` - Matern kernel; params `:order` (odd positive integer, default `1`, use `5` for Matern 5/2 etc.), `:theta` (default `1.0`, equal to `sqrt(2)` times the length scale of the scikit-learn Matern kernel) plus `:distance`; presets `:matern-12`, `:matern-32`, `:matern-52` fix `:order` to `1`, `3`, `5`.
  - `:geometric` - compactly supported kernel (support radius `2 r`); params `:n` (dimension, required, minimal map `{:n 4}`), `:r` (default `1.0`) plus `:distance`; presets `:triangular`, `:circular`, `:spherical` fix `:n` to `1`, `2`, `3`.
  - `:anova` - sum, over elementwise powers of the vectors, of a Gaussian-like term; params `:sigma`, `:k`, `:d` (all default `1.0`).
  - `:spline` - product, over vector elements, of a fixed cubic-spline expression, no parameters.
  - `:b-spline` - product, over vector elements, of a B-spline basis of degree `:n` (default `2.0`).
  - `:chi-square`, `:chi-square2`, `:histogram` - elementwise, sum-based kernels for non-negative vectors (e.g. histograms), no parameters.
  - `:generalized-histogram` - elementwise minimum of powered absolute values; param `:p` (default `2.0`).
  - `:dirichlet` - Dirichlet kernel; param `:n` dimensionality (default `1.0`).
  - `:hellinger` - Hellinger kernel (sum of the square roots of elementwise products), no parameters.

  Returns a function of two arguments `x` and `y` (vectors or numbers, matching `fastmath.vector` conventions) returning a `double`.

  Corner cases and exceptions:

  - `params` must be a map or `nil`. A number (for example `(kernel :gaussian 1.0)`, which was silently ignored before) or any other value throws `ex-info` with `:type :invalid-input`, the message names the keyword and the received value. Use `{:sigma 1.0}` instead.
  - A preset keyword (for example `:matern-32`) fixes some parameters. A `params` map holding a conflicting value for a fixed key (`(kernel :matern-32 {:order 5})`) throws `ex-info` (`:type :invalid-input`, with `:option`, `:value` and `:expected`); the same value (`{:order 3}`, `{:order 3.0}`) or other keys are accepted.
  - A missing required parameter throws `ex-info` naming the key (`:geometric` needs `:n`).
  - Throws `IllegalArgumentException` when `k` is not a registered kernel name.

  See also [[kernel-list]] for the full, up-to-date list of registered kernel names, [[rbf]] for radial (single-argument) kernels, [[exp]], [[scale]], [[sum]], [[product]], [[mult]], [[wmean]] for combining kernels, [[rbf->kernel]]."
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

(defn- check-kernels
  "Throws `ex-info` (`:type :invalid-input`) when `kernels` is empty or holds a value which is not a function."
  [who kernels]
  (when (empty? kernels)
    (throw (ex-info (str who " needs at least one kernel")
                    {:type :invalid-input :function who :value kernels :expected "at least one kernel"})))
  (doseq [k kernels]
    (when-not (ifn? k)
      (throw (ex-info (str who " expects kernel functions, got: " (pr-str k))
                      {:type :invalid-input :function who :value k :expected "a function of two arguments"})))))

(defn exp
  "Kernel wrapper. Creates the kernel `exp(t k(x, y))` of the kernel `k`.

  Parameters:

  - `k` - a kernel, a function of two arguments returning a double
  - `t` (optional) - multiplier of the kernel value inside the exponential (default: `1.0`)

  Returns a function `(fn [x y])` returning a double. The exponential of a positive definite kernel with `t >= 0` is positive definite.

  See also [[scale]], [[sum]], [[product]]."
  ([k] (exp k 1.0))
  ([k ^double t]
   (fn ^double [x y] (m/exp (m/* t (double (k x y)))))))

(defn scale
  "Kernel wrapper. Creates the kernel `factor k(x, y)` of the kernel `k`.

  Parameters:

  - `k` - a kernel, a function of two arguments returning a double
  - `factor` - multiplier of the kernel value

  Returns a function `(fn [x y])` returning a double. A positive definite kernel stays positive definite for `factor >= 0`.

  See also [[exp]], [[sum]], [[product]]."
  [k ^double factor]
  (fn ^double [x y] (m/* factor (double (k x y)))))

(defn sum
  "Kernel wrapper. Creates the kernel which is the sum of the values of the given kernels.

  Parameters:

  - `k1`, `ks` - one or more kernels, functions of two arguments returning a double

  Returns a function `(fn [x y])` returning `k1(x, y) + k2(x, y) + ...` as a double, summed from left to right. With one kernel the kernel itself is returned. The sum of positive definite kernels is positive definite.

  Throws `ex-info` when an argument is not a function.

  See also [[product]], [[wmean]], [[scale]]."
  ([k1] (check-kernels "sum" [k1]) k1)
  ([k1 k2]
   (check-kernels "sum" [k1 k2])
   (fn ^double [x y] (m/+ (double (k1 x y)) (double (k2 x y)))))
  ([k1 k2 & ks]
   (let [kernels (vec (list* k1 k2 ks))]
     (check-kernels "sum" kernels)
     (fn ^double [x y] (double (reduce (fn [^double acc k] (m/+ acc ^double (k x y))) 0.0 kernels))))))

(defn product
  "Kernel wrapper. Creates the kernel which is the product of the values of the given kernels.

  Parameters:

  - `k1`, `ks` - one or more kernels, functions of two arguments returning a double

  Returns a function `(fn [x y])` returning `k1(x, y) * k2(x, y) * ...` as a double, multiplied from left to right. With one kernel the kernel itself is returned. The product of positive definite kernels is positive definite.

  Throws `ex-info` when an argument is not a function.

  See also [[mult]] (the same function), [[sum]], [[wmean]], [[scale]]."
  ([k1] (check-kernels "product" [k1]) k1)
  ([k1 k2]
   (check-kernels "product" [k1 k2])
   (fn ^double [x y] (m/* (double (k1 x y)) (double (k2 x y)))))
  ([k1 k2 & ks]
   (let [kernels (vec (list* k1 k2 ks))]
     (check-kernels "product" kernels)
     (fn ^double [x y] (double (reduce (fn [^double acc k] (m/* acc (double (k x y)))) 1.0 kernels))))))

(defn mult
  "Kernel wrapper. Multiplies one or more kernels, the alias of [[product]] kept for compatibility.

  Parameters:

  - `k1`, `ks` - one or more kernels, functions of two arguments returning a double

  Returns a function `(fn [x y])` returning the product of the values of the kernels as a double. With one kernel the kernel itself is returned.

  See also [[product]], [[sum]]."
  [k1 & ks]
  (apply product k1 ks))

(defn wmean
  "Kernel wrapper. Creates the weighted mean of the values of the given kernels: `sum(w_i k_i(x, y)) / sum(w_i)`.

  Parameters:

  - `kernels` - a sequence of one or more kernels, functions of two arguments returning a double
  - `weights` (optional) - a sequence of numbers of the same length as `kernels` (default: all `1.0`, the arithmetic mean); the weights are normalised by their sum, so only their ratios matter

  Returns a function `(fn [x y])` returning the weighted mean as a double. For non-negative weights with a positive sum the mean of positive definite kernels is positive definite.

  Throws `ex-info` when `kernels` is empty, holds a value which is not a function, the number of weights differs from the number of kernels, a weight is not finite or the weights sum to zero.

  See also [[sum]], [[product]], [[scale]]."
  ([kernels] (wmean kernels (repeat (count kernels) 1.0)))
  ([kernels weights]
   (let [ks (vec kernels)
         ws (mapv double weights)
         total (double (reduce m/+ 0.0 ws))]
     (check-kernels "wmean" ks)
     (when (not= (count ks) (count ws))
       (throw (ex-info (str "wmean needs one weight per kernel, got " (count ws) " weights for " (count ks) " kernels")
                       {:type :invalid-input :function "wmean" :value weights :expected (count ks)})))
     (when (or (not (every? m/valid-double? ws)) (m/zero? total))
       (throw (ex-info (str "wmean needs finite weights with a non-zero sum, got: " (pr-str weights))
                       {:type :invalid-input :function "wmean" :value weights :expected "finite weights with a non-zero sum"})))
     (fn ^double [x y]
       (m// (double (reduce (fn [^double acc ^long i] (m/+ acc (m/* (double (ws i)) (double ((ks i) x y)))))
                            0.0 (range (count ks))))
            total)))))

;; kernel density estimation

(defn bandwidth
  "Estimates the 1d kernel density bandwidth (`h`) for `data`.

  The bandwidth controls how much each data point is spread out by the kernel; smaller values follow the data more closely (risking overfitting/noise), larger values produce a smoother, more biased estimate.

  Parameters:

  - `kernel` - a keyword naming a kernel registered in [[kernel-list]] under `:kde`, or a kernel function (the selectors `:nrd-adjust`, `:sj` and `:sj-dpi` need a registered kernel).
  - `data` - a sequence of numbers to estimate the bandwidth from.
  - `h` - a number (returned as it is) or the keyword of a selector (`nil` means `:nrd`):
      - `:nrd` - the normal reference rule `1.06 * s * n^-0.2`, `s = min(sd, IQR / 1.34)` with the type 7 quantiles of R.
      - `:nrd0` - Silverman's rule of thumb `0.9 * s * n^-0.2`.
      - `:nrd-adjust` - `:nrd` multiplied by the canonical bandwidth factor of the kernel; not supported for `:silverman`, `:cauchy` and custom kernels.
      - `:sj`, `:sj-dpi` - the Sheather-Jones solve-the-equation and direct plug-in selectors (`bw.SJ` of R) scaled to the kernel.
      - `:lscv` - least squares cross-validation.
      - `:lcv` - likelihood cross-validation.
      - `:rlcv` - robust likelihood cross-validation (Wu, 2019).
  - `opts` (optional map): `:info?` (return a map with `:h`, `:value`, `:at-bound?` and `:status`), `:bounds`, `:grid`, `:bins`, `:rlcv-form`, `:radius`, `:extent`; see [[fastmath.kernel.density/bandwidth]] for the algorithms, the search interval of the cross validation and the defaults.

  Returns the estimated bandwidth as a `double` (a map with `:info?`).

  Throws `ex-info` (`:type :invalid-input`) for invalid data (empty, not finite, a single point for the selectors), an unknown kernel or selector and invalid options.

  See also [[kernel-density]], [[kernel-density-ci]]."
  ([kernel data h] (dens/bandwidth kernel data h))
  ([kernel data h opts] (dens/bandwidth kernel data h opts)))

(defn- density-params
  "The map of parameters of the density functions: a number or a keyword is a bandwidth."
  [bandwidth]
  (if (or (number? bandwidth) (keyword? bandwidth))
    {:bandwidth bandwidth}
    bandwidth))

(defn kernel-density
  "Returns a 1d kernel density estimation (KDE) for `data`.

  A kernel density estimate is a smoothed histogram: it places a scaled copy of `kernel` at every data point and sums the contributions, weighted by a bandwidth that controls how much each point spreads out. The evaluation is exact: every data point within the radius of the kernel contributes (all points for `:cauchy` and for a custom kernel without `:radius`).

  Parameters:

  - `kernel` - a keyword naming a kernel registered in [[kernel-list]] under `:kde`, or a custom kernel function accepting and returning a `double`.
  - `data` - a sequence of numbers to estimate the density from (finite, not empty).
  - `bandwidth` (optional) - a number or a keyword of a selector (used as the bandwidth `h`), or a map `params`:
      - `:bandwidth` - the bandwidth `h`, a number or one of the keywords accepted by [[bandwidth]] (`:nrd`, `:nrd0`, `:nrd-adjust`, `:sj`, `:sj-dpi`, `:rlcv`, `:lcv`, `:lscv`). Default: `:nrd`.
      - `:bandwidth-options` - a map of the options of [[bandwidth]] (`:bounds`, `:grid`, `:bins`, `:rlcv-form`).
      - `:binned?` - `false` (default) for the exact sum. `true` or a number: the data are replaced by a linear binning on a regular grid with the spacing `h / divisor` (`5` for `true`, the number otherwise); a finite number greater than 0.
      - `:radius`, `:extent` - for a custom kernel: the distance beyond which the kernel is ignored (default `##Inf`) and beyond which the density is zero for the support `:mn`, `:mx` (default `5.0`).
  - `info?` (optional, default `false`) - when `true`, returns the estimator together with its fitting details instead of the bare function.

  Returns a function of a single `double` returning the density estimate at that point. When `info?` is `true`, returns a map instead:

  - `:kde` - the density function described above.
  - `:factor` - the normalizing factor, equal to `1/(n*h)`.
  - `:h` - the bandwidth actually used.
  - `:mn` and `:mx` - the support of the estimator (data range expanded by `1.01 * h * extent` of the kernel).

  Throws `ex-info` (`:type :invalid-input`) for invalid input.

  See also [[bandwidth]] (bandwidth estimation alone), [[kernel-density-ci]] (density with confidence intervals)."
  ([kernel data] (kernel-density kernel data {:bandwidth :nrd}))
  ([kernel data bandwidth] (kernel-density kernel data bandwidth false))
  ([kernel data bandwidth info?]
   (let [params (density-params bandwidth)
         f (if info? dens/kernel-density+ dens/kernel-density)]
     (f kernel data params))))

(defn kernel-density-ci
  "Returns a 1d KDE for `data` that also reports pointwise asymptotic confidence intervals.

  Uses the asymptotic normal approximation of the kde variance (see section 6.1.5 of http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/tutorials/xlghtmlnode33.html): for a density estimate `f(x)` obtained with normalizing factor `1/(n*h)`, the variance is approximated as `k2 * f(x) / (n*h)`, where `k2` is the kernel roughness (integral of the squared kernel), and confidence bounds are derived from the normal quantile at the requested `alpha`. The bias of the estimate is ignored and the lower bound is not clipped at zero.

  Parameters:

  - `kernel` - a keyword naming a kernel registered in [[kernel-list]] under `:kde` (a custom kernel function is not accepted, since its roughness constant `k2` must be looked up).
  - `data` - a sequence of numbers to estimate the density from.
  - `bandwidth-or-params` (optional) - a number or a keyword of a selector (the bandwidth), or a map with the same keys as [[kernel-density]] plus:
      - `:alpha` - the confidence level parameter, a number in `(0, 1)`, the resulting interval has coverage `1 - alpha`. Default: `0.05`.

  Returns a function of a single `double` `x` returning a 3-element vector `[fx lower upper]`, where `fx` is the estimated density at `x`, and `lower`/`upper` are the corresponding confidence bounds.

  Throws `ex-info` (`:type :invalid-input`) when `kernel` is not one of the kernels registered under `:kde`, for `alpha` outside of `(0, 1)` and for invalid input of [[kernel-density]].

  See also [[kernel-density]], [[bandwidth]]."
  ([kernel data] (kernel-density-ci kernel data {:bandwidth :nrd}))
  ([kernel data bandwidth-or-params]
   (dens/kernel-density-ci kernel data (density-params bandwidth-or-params))))

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
   :cosine win/cosine
   :kaiser win/kaiser
   :general-gaussian win/general-gaussian
   :flat-top win/flat-top
   :flat-top-3 win/flat-top-3
   :flat-top-hp win/flat-top-hp
   :flat-top-srs win/flat-top-srs
   :flat-top-bk win/flat-top-bk
   :hft70 win/hft70
   :hft90d win/hft90d
   :hft95 win/hft95
   :hft116d win/hft116d
   :hft144d win/hft144d
   :hft169d win/hft169d
   :hft196d win/hft196d
   :hft223d win/hft223d
   :hft248d win/hft248d})

(defn- window-length!
  "Returns `N` as a long or throws `ex-info` unless it is an integer valued number `>= 0`."
  ^long [N]
  (if (and (number? N)
           (let [d (double N)]
             (and (m/valid-double? d) (m/== d (Math/floor d)) (m/<= 0.0 d 2.147483647E9))))
    (long (double N))
    (throw (ex-info (str "window, N: must be an integer >= 0, got: " (pr-str N))
                    {:type :invalid-input :function :window :property :N :value N :expected "integer >= 0"}))))

(defn- window-result!
  "Returns the coefficients `res` of a custom window function or throws `ex-info` when their number is not `expected`."
  [res expected]
  (let [res (if (nil? res) () res)
        len (long (cond
                    (sequential? res) (count res)
                    (.isArray ^Class (class res)) (java.lang.reflect.Array/getLength res)
                    :else -1))]
    (when-not (m/== len (long expected))
      (throw (ex-info (str "window, custom function: expected a sequence of " expected " coefficients, got "
                           (if (m/neg? len) (pr-str res) (str len " coefficients")))
                      {:type :invalid-input :function :window :property :window-name
                       :value res :expected (str "a sequence of " expected " numbers")})))
    res))

(defn window
  "Returns tapering window coefficients, `N` numbers shaped to fade toward the edges.

  Windows (also called taper functions) are used to reduce spectral leakage before a Fourier transform, to design FIR filters, or to smooth signal edges. The windows, their formulas, the scale and the equivalents in scipy and Octave are documented in [[fastmath.kernel.window]].

  Parameters:

  - `window-name` - a keyword naming one of the windows registered in [[kernel-list]] under `:window`, or a custom function `(fn [N options] ...)` returning a sequence of exactly `N` coefficients.
  - `N` - the number of samples, an integer `>= 0` (default: `256`). `N = 0` gives `()` and `N = 1` gives `(1.0)` for the default `:normalize?` (the symmetric and the periodic window).
  - `options` (optional map), common to every window:
      - `:symmetric?` - `true` (default) for a symmetric window (suitable for FIR filter design), `false` for a periodic window (suitable for spectral analysis / FFT-based methods). A periodic window of size `N` is obtained by sampling a symmetric window of size `N+1` and dropping the last coefficient. `nil` is the default.
      - `:normalize?` - controls the scale of the window. `true` (the default for most windows) divides by the value of the window function at the middle (continuous model windows: for an even `N` the maximum is below `1.0`, as in scipy) or by the maximum (windows defined on the index), so that the maximum is `1.0`; `false` leaves the raw values (the integral of the continuous function over the window is `1.0` and the mean of the Doerry discrete windows is `1.0`); other values are `:LInf` (alias `:Linf`), `:L1`, `:L2`, `:N` and a finite non-zero number to divide by, see [[fastmath.kernel.window/normalize-coefficients]]. `nil` is the default of the window, any other value throws. The default is `false` for `:rectangular05`, `:kaiser-bessel-derived`, `:vorbis` and `:cosine`.

  Some windows accept additional `options` keys; a value out of the documented range (given in [[fastmath.kernel.window]]) or a `NaN` throws `ex-info`, `nil` is the default:

  - `:triangular` - `:shift` (long, default `2`).
  - `:parzen` - `:discrete?` (default `true`).
  - `:b-spline` - `:order` (long, default `3`).
  - `:connes`, `:exponential`, `:hanning-poisson`, `:gaussian`, `:cauchy`, `:kaiser-bessel`, `:kaiser-bessel-derived` - `:alpha` shape parameter (default `1.0`).
  - `:kaiser` - `:beta` shape parameter of the Kaiser window (default `pi`, which is `:kaiser-bessel` with `:alpha 1.0`).
  - `:general-gaussian` - `:p` shape (default `1.5`), `:sig` standard deviation relative to the window length (default `0.25`).
  - `:parzen-algebraic` - `:gamma` (default `1.0`), `:u` (default `3.0`).
  - `:singla-singh` - `:order` (default `1.0`).
  - `:lanczos` - `:L` (double, default `3.0`; `1.0` is the window of scipy).
  - `:raised-cosine` - `:alpha` (default `0.5`).
  - `:webster-hamming` - `:v` (default `1.0`).
  - `:power-of-cosine` - `:m` (default `1.0`).
  - `:raised-power-of-cosine` - `:alpha` (default `0.05`), `:m` (default `1.0`).
  - `:parzen-cosine` - `:gamma` (default `1.0`), `:m` (default `2.0`).
  - `:trapezoid` - `:alpha` (default `0.25`).
  - `:tukey` - `:alpha` (default `0.5`).
  - `:blackman-harris-family` - `:coeffs` vector of the amplitudes of the cosine terms (default `[0.5 0.5]`); the window is not clipped at zero, so a negative lobe is kept.
  - `:parzen-exponential`, `:parzen-geometric` - `:alpha` (default `1.0`), `:r` (default `1.0`).
  - `:dolph-chebyshev`, `:saramaki`, `:legendre` - `:level` target sidelobe attenuation in dB (default `-50.0`).
  - `:taylor` - `:level` (default `-50.0`), `:n` number of nearly constant-level sidelobes, a long (default `4`).
  - `:cosh`, `:avci-nacaroglu`, `:knab`, `:bessel-I1` - `:alpha` shape parameter (default `2.0`).
  - `:ultraspherical` - `:level` (default `-50.0`), `:alpha` (default `2.0`).

  The remaining registered windows only accept the common options: `:rectangular`, `:welch`, `:sinc`, `:fejer`, `:de-la-vallee-poussin`, `:hamming`, `:hamming-exact`, `:hann`, `:bohman`, `:bartlett-hann`, `:cosine`, `:vorbis`, `:blackman` and `:blackman-exact`, the `:blackman-harris*`, `:nutall-*`, `:blackman-nutall` and `:mottaghi-kashtiban-shayesteh`, `:shayesteh-kashtiban`, and the flat-top windows `:flat-top`, `:flat-top-3`, `:flat-top-hp`, `:flat-top-srs`, `:flat-top-bk` (Octave flattopwin), `:hft70`, `:hft90d`, `:hft95`, `:hft116d`, `:hft144d`, `:hft169d`, `:hft196d`, `:hft223d`, `:hft248d`.

  Returns a sequence of `N` `double` window coefficients.

  Throws `ex-info` (`:type :invalid-input`) when `window-name` is neither a registered keyword (the message lists the valid keys) nor a function, when `N` is not an integer `>= 0`, when a custom window function does not return `N` coefficients, and for an invalid option value.

  See also [[kernel-list]] for the full, up-to-date list of registered window names."
  ([window-name] (window window-name 256))
  ([window-name N] (window window-name N nil))
  ([window-name N options]
   (let [wind (or (windows window-name)
                  (when (fn? window-name) window-name))
         _ (when-not wind
             (throw (ex-info (str "Unknown window " (pr-str window-name) ", expected a function or one of " (vec (sort (keys windows))))
                             {:type :invalid-input :function :window :property :window-name
                              :value window-name :expected "a window function or one of the keys of kernel-list :window"})))
         N (window-length! N)
         symmetric? (if (nil? (:symmetric? options)) true (boolean (:symmetric? options)))]
     (if (or symmetric? (m/<= N 1))
       (window-result! (wind N options) N)
       (take N (window-result! (wind (m/inc N) options) (m/inc N)))))))

;;

(def ^{:doc "Map of the names of the available kernels, each value is a sorted sequence of keywords.

  Keys:

  - `:vector` - names accepted by [[kernel]] (vector kernels and their presets)
  - `:rbf` - names accepted by [[rbf]] (radial basis functions and their presets)
  - `:kde` - names accepted by [[kernel-density]], [[kernel-density-ci]] and [[bandwidth]]
  - `:window` - names accepted by [[window]] (tapering windows)

  The lists are built when the namespace is loaded; a kernel registered later with `defmethod` is not included."}
  kernel-list
  {:vector (sort (keys (methods kernel)))
   :rbf (sort (keys (methods rbf)))
   :kde (sort (keys dens/kde-data))
   :window (sort (keys windows))})

(def ^{:doc "Converts a radial basis function into a vector kernel using a distance function; the function of `fastmath.kernel.vector`.

  Parameters:

  - `rbf-kernel` - a function of one double (the distance), for example created by [[rbf]]
  - `distance` (optional) - distance function of two points (default: euclidean)

  Returns a function `(fn [x y])` which applies `rbf-kernel` to the distance of `x` and `y`.

  See also [[rbf]], [[kernel]]."
       :arglists '([rbf-kernel] [rbf-kernel distance])}
  rbf->kernel vk/rbf->kernel)

