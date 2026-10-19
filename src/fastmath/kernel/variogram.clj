(ns fastmath.kernel.variogram
  "Parametric semivariogram models, empirical estimators, and model fitting for geostatistical analysis.

  A semivariogram describes how spatial dependence between observations decays with lag
  distance `h`. It is the primary structure function in kriging interpolation: given a set
  of spatially distributed observations, the empirical semivariogram is first estimated
  from the data, then a parametric model is fitted to it, and the fitted model is used by
  a kriging algorithm to derive optimal interpolation weights. Each parametric model is
  governed by up to four parameters: `nugget` (discontinuity at the origin representing
  micro-scale variation or measurement error), `psill` (partial sill; the structured
  spatial variance), `range` (the characteristic lag scale), and in some models `beta`
  (a shape or power exponent).

  Parametric semivariogram models — each takes a parameter map and returns a
  `(fn ^double [^double x])` computing semivariance at lag `x`:

  - Bounded models with a hard sill at `range`: `spherical`, `circular`, `linear`,
    `cubic`, `pentaspherical`, `tpower`
  - Bounded models that asymptotically approach the sill: `exponential`, `gaussian`,
    `expower`, `cauchy`, `rational`, `hole`
  - Unbounded models (no sill): `power` (valid for the exponent `beta` in `(0, 2]`)
  - Generalised families with a configurable shape order, each available as a
    two-step creator `->name` / convenience wrapper `name`: `->bessel` / `bessel`,
    `->matern` / `matern`, `->hyperspherical` / `hyperspherical`,
    `->superspherical` / `superspherical`, `->tplstable` / `tplstable`
  - RBF-based model: `rbf->variogram` wraps any suitable RBF kernel from
    `fastmath.kernel.rbf` into a semivariogram constructor

  Empirical semivariogram:

  - `empirical` — bins all observation pairs by lag distance and applies a chosen
    estimator to compute `gamma` per bin; returns a sequence of `{:h :gamma :n}` maps
    (`:n` is the number of pairs in the bin)

  Estimators used by `empirical` to compute `gamma` per bin from pairwise differences:

  - `matheron-estimator` — classical (Matheron) estimator
  - `cressie-estimator` — robust Cressie estimator
  - `highly-robust-estimator` — Genton highly robust estimator
  - `dowd-estimator` — Dowd estimator
  - `->quantile-estimator` — creates a quantile (Armstrong and Delfiner) estimator for a given level
  - `robust-m-estimator` — Gunst and Hartfield robust M-estimator

  Model fitting (L-BFGS-B optimisation against the empirical semivariogram):

  - `fit` — fits a model and returns the semivariogram function
  - `fit-params` — fits a model and returns a detailed result map with the fitted
    parameters, weights, and loss value

  Utilities:

  - `variogram->covariance` — converts a semivariogram function and the total sill to the covariance function
  - `bounding-box-diagonal` — computes the diagonal of the bounding box of spatial points;
    used to derive the default cutoff in `empirical`
  - `remove-outliers-fence` — removes outliers using Tukey's fences criterion
  - `remove-outliers-mad` — removes outliers using median absolute deviation
  - `maybe-remove-outliers` — applies one of the two filters (or none) to the points before `empirical` forms the pairs

  Conventions: a model is built from a parameter map (`:nugget`, `:psill`, `:range`, and `:beta` or
  an order where needed), the parameters are checked once at the construction and `ex-info` with
  `:type :invalid-input` is thrown for a missing or invalid one. The returned semivariogram
  function is exactly `0.0` at lag zero (the nugget acts just after zero), the bounded models are constant
  beyond `range`, the lag is expected to be non-negative. An RBF kernel given to `fit` or `rbf->variogram`
  is a function of the scaled lag; a model constructor such as `gaussian` is not a kernel."
  (:require [fastmath.core :as m]
            [fastmath.vector :as v]
            [fastmath.distance :as dist]
            [fastmath.stats :as stats]
            [fastmath.optimization :as optim]
            [fastmath.random :as r]
            [fastmath.special :as special]
            [fastmath.special.poly :as spoly])
  (:import [fastmath.vector Vec2]))

(set! *unchecked-math* :warn-on-boxed)
(set! *warn-on-reflection* true)

(defn- invalid-parameter
  [model property value expected]
  (throw (ex-info (str "Invalid parameter " property " of the " model " variogram model: expected "
                       expected ", got " (pr-str value))
                  {:type :invalid-input :model model :property property :value value :expected expected})))

(defn- check-finite
  ^double [model property value]
  (if (and (number? value) (Double/isFinite (double value)))
    (double value)
    (double (invalid-parameter model property value "a finite number"))))

(defn- check-nonneg!
  [model property value]
  (when (m/neg? (check-finite model property value))
    (invalid-parameter model property value "a number >= 0")))

(defn- check-pos!
  [model property value]
  (when-not (m/pos? (check-finite model property value))
    (invalid-parameter model property value "a number > 0")))

(defn- check-sill-range!
  [model nugget psill range]
  (check-nonneg! model :nugget nugget)
  (check-nonneg! model :psill psill)
  (check-pos! model :range range))

(defn- check-beta!
  [model beta max-beta]
  (let [b (check-finite model :beta beta)
        mx (double max-beta)]
    (when-not (and (m/pos? b) (m/<= b mx))
      (invalid-parameter model :beta beta (if (m/pos-inf? mx) "a number > 0" (str "a number in (0, " mx "]"))))))

(defn- check-hurst!
  [model h]
  (let [v (check-finite model :H h)]
    (when-not (and (m/pos? v) (m/< v 1.0))
      (invalid-parameter model :H h "a number in (0, 1)"))))

(defn power
  "Power semivariogram model.

  The formula is: `nugget + psill * (x / range)^beta`

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the scaling coefficient for the power term;
        must be non-negative
      - `:range` (double): the range parameter controlling the horizontal scale of the
        power curve; must be positive
      - `:beta` (double): the power exponent controlling the rate of growth; must be in
        `(0.0, 2.0]`, the only interval where the model is conditionally negative definite (a valid
        semivariogram); `2.0` is the parabola, `1.0` the linear model without a sill

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero. The model has no sill: it grows without bound (`Infinity` for an infinite lag).

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range (including `:beta` outside `(0.0, 2.0]`), throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[exponential]], [[gaussian]], [[spherical]], [[expower]], [[tpower]]."
  [{:keys [^double nugget ^double psill ^double range ^double beta]}]
  (check-sill-range! :power nugget psill range)
  (check-beta! :power beta 2.0)
  (fn ^double [^double x]
    (if (m/zero? x) 0.0
        (->> (m/pow (m// x range) beta)
             (m/* psill)
             (m/+ nugget)))))

(defn exponential
  "Exponential semivariogram model.

  The formula is: `nugget + psill * (1 - exp(-x / range))`

  `range` is the scale of the exponential (the gstat `Exp` range), not the practical range: the structured part reaches 95% of `psill` at about `3 * range` (exactly `1 - exp(-3)` = 0.9502 of it at `3 * range`) and never reaches the sill.

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the maximum contribution of the structured
        spatial variance above the nugget; must be non-negative
      - `:range` (double): the scale parameter controlling how fast the semivariance
        rises towards the sill (the sill is approached asymptotically); must be positive

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[gaussian]], [[spherical]], [[power]], [[expower]], [[hole]]."
  [{:keys [^double nugget ^double psill ^double range]}]
  (check-sill-range! :exponential nugget psill range)
  (fn ^double [^double x]
    (if (m/zero? x) 0.0
        (->> (m// x range)
             (m/-)
             (m/exp)
             (m/- 1.0)
             (m/* psill)
             (m/+ nugget)))))

(defn expower
  "Exponential-power (stable) semivariogram model.

  The formula is: `nugget + psill * (1 - exp(-(x / range)^beta))`

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the maximum contribution of the structured
        spatial variance above the nugget; must be non-negative
      - `:range` (double): the range parameter controlling the horizontal scale; must be positive
      - `:beta` (double): the power exponent controlling the curvature near the origin;
        must be in `(0.0, 2.0]` (the valid range of the model); `1.0` is [[exponential]], `2.0` is [[gaussian]]

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[exponential]], [[gaussian]], [[power]], [[cauchy]]."
  [{:keys [^double nugget ^double psill ^double range ^double beta]}]
  (check-sill-range! :expower nugget psill range)
  (check-beta! :expower beta 2.0)
  (fn ^double [^double x]
    (if (m/zero? x) 0.0
        (->> (m/pow (m// x range) beta)
             (m/-)
             (m/exp)
             (m/- 1.0)
             (m/* psill)
             (m/+ nugget)))))

(defn gaussian
  "Gaussian semivariogram model.

  The formula is: `nugget + psill * (1 - exp(-(x / range)^2))`

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the maximum contribution of the structured
        spatial variance above the nugget; must be non-negative
      - `:range` (double): the range parameter controlling the horizontal scale; must be positive

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[exponential]], [[expower]], [[spherical]], [[cubic]]."
  [{:keys [^double nugget ^double psill ^double range]}]
  (check-sill-range! :gaussian nugget psill range)
  (fn ^double [^double x]
    (if (m/zero? x) 0.0
        (->> (m// x range)
             (m/sq)
             (m/-)
             (m/exp)
             (m/- 1.0)
             (m/* psill)
             (m/+ nugget)))))

(defn rational
  "Rational quadratic semivariogram model.

  The formula is: `nugget + psill * (1 - (1 + (x / range)^2 / beta)^(-beta))`

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the maximum contribution of the structured
        spatial variance above the nugget; must be non-negative
      - `:range` (double): the range parameter controlling the horizontal scale; must be positive
      - `:beta` (double): the shape exponent controlling curvature and convergence rate;
        must be positive; values below `EPSILON` are treated as `EPSILON`

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[cauchy]], [[gaussian]], [[expower]], [[power]]."
  [{:keys [^double nugget ^double psill ^double range ^double beta]}]
  (check-sill-range! :rational nugget psill range)
  (check-beta! :rational beta ##Inf)
  (let [b (m/max beta m/EPSILON)
        -b (m/- b)]
    (fn ^double [^double x]
      (if (m/zero? x) 0.0
          (->> (m/pow (m/inc (m// (m/sq (m// x range)) b)) -b)
               (m/- 1.0)
               (m/* psill)
               (m/+ nugget))))))

(defn spherical
  "Spherical semivariogram model.

  The formula is:
  - `x = 0`: `0.0`
  - `x <= range`: `nugget + psill * (1.5 * (x / range) - 0.5 * (x / range)^3)`
  - `x > range`: `nugget + psill`

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the structured spatial variance above the nugget;
        must be non-negative
      - `:range` (double): the range at which the sill is reached exactly; must be positive

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[linear]], [[cubic]], [[pentaspherical]], [[circular]], [[gaussian]]."
  [{:keys [^double nugget ^double psill ^double range]}]
  (check-sill-range! :spherical nugget psill range)
  (fn ^double [^double x]
    (cond
      (m/zero? x) 0.0
      (m/< range x) (m/+ nugget psill)
      :else (let [xr (m// x range)]
              (->> (m/cb xr)
                   (m/* 0.5)
                   (m/- (m/* 1.5 xr))
                   (m/* psill)
                   (m/+ nugget))))))

(defn cubic
  "Cubic semivariogram model.

  The formula is:
  - `x = 0`: `0.0`
  - `x <= range`: `nugget + psill * (7*(x/range)^2 - 8.75*(x/range)^3 + 3.5*(x/range)^5 - 0.75*(x/range)^7)`
  - `x > range`: `nugget + psill`

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the structured spatial variance above the nugget;
        must be non-negative
      - `:range` (double): the range at which the sill is reached exactly; must be positive

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[spherical]], [[pentaspherical]], [[gaussian]], [[circular]]."
  [{:keys [^double nugget ^double psill ^double range]}]
  (check-sill-range! :cubic nugget psill range)
  (fn ^double [^double x]
    (cond
      (m/zero? x) 0.0
      (m/< range x) (m/+ nugget psill)
      :else (let [xr (m// x range)
                  xr2 (m/sq xr)
                  xr3 (m/* xr2 xr)
                  xr5 (m/* xr3 xr2)
                  xr7 (m/* xr5 xr2)]
              (->> (m/+ (m/* 7.0 xr2)
                        (m/* -8.75 xr3)
                        (m/* 3.5 xr5)
                        (m/* -0.75 xr7))
                   (m/* psill)
                   (m/+ nugget))))))

(defn pentaspherical
  "Pentaspherical semivariogram model.

  The formula is:
  - `x = 0`: `0.0`
  - `x <= range`: `nugget + psill * (1.875*(x/range) - 1.25*(x/range)^3 + 0.375*(x/range)^5)`
  - `x > range`: `nugget + psill`

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the structured spatial variance above the nugget;
        must be non-negative
      - `:range` (double): the range at which the sill is reached exactly; must be positive

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[spherical]], [[cubic]], [[circular]], [[linear]]."
  [{:keys [^double nugget ^double psill ^double range]}]
  (check-sill-range! :pentaspherical nugget psill range)
  (fn ^double [^double x]
    (cond
      (m/zero? x) 0.0
      (m/< range x) (m/+ nugget psill)
      :else (let [xr (m// x range)
                  xr2 (m/* xr xr)]
              (->> (m/* xr2 0.375)
                   (m/+ -1.25)
                   (m/* xr2)
                   (m/+ 1.875)
                   (m/* psill xr)
                   (m/+ nugget))))))

(defn circular
  "Circular semivariogram model.

  The formula is:
  - `x = 0`: `0.0`
  - `x <= range`: `nugget + psill * (2/pi) * ((x/range) * sqrt(1 - (x/range)^2) + arcsin(x/range))`
  - `x > range`: `nugget + psill`

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the structured spatial variance above the nugget;
        must be non-negative
      - `:range` (double): the range at which the sill is reached exactly; must be positive

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[spherical]], [[pentaspherical]], [[cubic]], [[linear]]."
  [{:keys [^double nugget ^double psill ^double range]}]
  (check-sill-range! :circular nugget psill range)
  (fn ^double [^double x]
    (cond
      (m/zero? x) 0.0
      (m/< range x) (m/+ nugget psill)
      :else (let [xr (m// x range)]
              (->> (m/* m/TWO_INV_PI
                        (m/+ (m/* xr (m/sqrt (m/- 1.0 (m/* xr xr))))
                             (m/asin xr)))
                   (m/* psill)
                   (m/+ nugget))))))

(defn linear
  "Linear semivariogram model.

  The formula is:
  - `x = 0`: `0.0`
  - `x <= range`: `nugget + psill * (x / range)`
  - `x > range`: `nugget + psill`

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the structured spatial variance above the nugget;
        must be non-negative
      - `:range` (double): the range at which the sill is reached exactly; must be positive

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[spherical]], [[power]], [[tpower]], [[circular]]."
  [{:keys [^double nugget ^double psill ^double range]}]
  (check-sill-range! :linear nugget psill range)
  (fn ^double [^double x]
    (cond
      (m/zero? x) 0.0
      (m/< range x) (m/+ nugget psill)
      :else (->> (m// x range)
                 (m/* psill)
                 (m/+ nugget)))))

(defn tpower
  "Truncated power semivariogram model.

  The formula is:
  - `x = 0`: `0.0`
  - `x <= range`: `nugget + psill * (1 - (1 - x / range)^beta)`
  - `x > range`: `nugget + psill`

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the structured spatial variance above the nugget;
        must be non-negative
      - `:range` (double): the range at which the sill is reached exactly; must be positive
      - `:beta` (double): the power exponent controlling the shape of the rise within the
        range; must be positive; `1.0` gives a linear rise

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[power]], [[linear]], [[spherical]], [[expower]]."
  [{:keys [^double nugget ^double psill ^double range ^double beta]}]
  (check-sill-range! :tpower nugget psill range)
  (check-beta! :tpower beta ##Inf)
  (fn ^double [^double x]
    (cond
      (m/zero? x) 0.0
      (m/< range x) (m/+ nugget psill)
      :else (->> (m/pow (m/- 1.0 (m// x range)) beta)
                 (m/- 1.0)
                 (m/* psill)
                 (m/+ nugget)))))

(defn hole
  "Hole effect semivariogram model.

  The formula is: `nugget + psill * (1 - sinc(x / range))`

  where `sinc(t) = sin(pi * t) / (pi * t)` is the normalised sinc, with `sinc(0) = 1`. The model equals the gstat `Wav` model (the gstat `Hol` model uses the unnormalised sinc, so it corresponds to `range / pi` here). The semivariance oscillates around the sill and tends to `nugget + psill` for an infinite lag.

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the amplitude of the oscillation around the sill;
        must be non-negative
      - `:range` (double): the scale parameter controlling the period of oscillation (the
        semivariance reaches the sill for the first time at `range`); must be positive

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[exponential]], [[cauchy]], [[bessel]]."
  [{:keys [^double nugget ^double psill ^double range]}]
  (check-sill-range! :hole nugget psill range)
  (fn ^double [^double x]
    (cond
      (m/zero? x) 0.0
      (m/pos-inf? x) (m/+ nugget psill)
      :else (->> (m// x range)
                 (m/sinc)
                 (m/- 1.0)
                 (m/* psill)
                 (m/+ nugget)))))

(defn cauchy
  "Generalized Cauchy semivariogram model.

  The formula is: `nugget + psill * (1 - (1 + (x / range)^2)^(-beta))`

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the maximum contribution of the structured
        spatial variance above the nugget; must be non-negative
      - `:range` (double): the range parameter controlling the horizontal scale; must be positive
      - `:beta` (double): the decay exponent controlling the rate of convergence to the sill
        and the heaviness of the correlation tail; must be positive for a valid model

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[rational]], [[gaussian]], [[exponential]], [[matern]]."
  [{:keys [^double nugget ^double psill ^double range ^double beta]}]
  (check-sill-range! :cauchy nugget psill range)
  (check-beta! :cauchy beta ##Inf)
  (fn ^double [^double x]
    (if (m/zero? x) 0.0
        (->> (-> (m// x range)
                 (m/sq)
                 (m/inc)
                 (m/pow beta)
                 (m//))
             (m/- 1.0)
             (m/* psill)
             (m/+ nugget)))))

;; Correlation deficits (1 - correlation) of the Bessel and Matern families, computed so that they stay
;; finite and accurate for large orders and tiny lags, where the textbook product of a power and a Bessel
;; function is 0 * Inf (see the docstrings of `->bessel` and `->matern`).

(defn- bessel-series
  "1 - 0F1(; order + 1; -d^2 / 4) as the alternating series without the leading one."
  ^double [^double order ^double d]
  (let [q (m/* 0.25 d d)
        c (m/inc order)]
    (loop [k (long 1)
           mag 1.0
           sign 1.0
           s 0.0]
      (let [mag (m/* mag q (m// (m/* (double k) (m/+ c (m/dec k)))))
            s (m/+ s (m/* sign mag))]
        (if (or (m/< mag (m/* 1.0e-18 s)) (m/> k 200))
          s
          (recur (m/inc k) mag (m/- sign) s))))))

(defn- bessel-deficit
  "1 - Gamma(order + 1) (2 / d)^order J_order(d)."
  ^double [^double order ^double g ^double d]
  (cond
    (m/pos-inf? d) 1.0
    (m/< (m/* d d) (m/* 4.0 (m/inc order))) (bessel-series order d)
    :else (let [scale (if (m/< order 100.0)
                        (m/* g (m/pow (m// 2.0 d) order))
                        (m/exp (m/+ (special/log-gamma (m/inc order))
                                    (m/* order (m/log (m// 2.0 d))))))]
            (m/max 0.0 (m/- 1.0 (m/* scale (special/bessel-J order d)))))))

(def ^:const ^:private matern-large-order 20.0)

(defn- matern-series
  "1 - 2^(1-order) / Gamma(order) d^order K_order(d) for small d as minus the series of the regular
  part 0F1(; 1 - order; d^2 / 4) (the singular part is of the order d^(2 order)); order >= 20."
  ^double [^double order ^double d]
  (let [q (m/* 0.25 d d)]
    (loop [k (long 1)
           t 1.0
           s 0.0]
      (let [kd (double k)
            t (m/* t q (m// (m/* kd (m/- kd order))))
            s (m/- s t)]
        (if (or (m/< (m/abs t) (m/* 1.0e-18 (m/abs s))) (m/> k 9))
          s
          (recur (m/inc k) t s))))))

(defn- matern-debye
  "2^(1-order) / Gamma(order) d^order K_order(d) from the uniform (Debye) expansion of K in the log
  domain; no overflow for large orders."
  ^double [^double order ^double d]
  (let [z (m// d order)
        zs (m/hypot 1.0 z)
        s (if (m/< z 1.0) (m// (m/* z z) (m/inc zs)) (m/dec zs))
        psi (m/- (m/log1p (m/* 0.5 s)) s)
        t (m// zs)
        ^Vec2 res (spoly/split-poly (m/- (m// t order)) (spoly/uk-poly-10 (m/* t t)))
        w (m// (m/* order order))
        stirling (m/* (m// order)
                      (m/muladd w (m/muladd w (m/muladd w (m/muladd w (m// 1.0 1188.0) (m// -1.0 1680.0))
                                              (m// 1.0 1260.0))
                                      (m// -1.0 360.0))
                                (m// 1.0 12.0)))]
    (m/exp (m/+ (m/* order psi)
                (m/* -0.5 (m/log zs))
                (m/- stirling)
                (m/log (.y res))))))

(defn- matern-deficit
  "1 - 2^(1-order) / Gamma(order) d^order K_order(d); g is the constant 2^(1-order) / Gamma(order)."
  ^double [^double order ^double g ^double d]
  (cond
    (m/nan? d) ##NaN
    (m/pos-inf? d) 1.0
    (m/>= order matern-large-order) (if (m/< (m/* 0.25 d d) (m/* 0.01 order))
                                      (matern-series order d)
                                      (m/- 1.0 (matern-debye order d)))
    :else (let [v (m/* g (m/pow d order) (special/bessel-K order d))]
            (if (Double/isFinite v) (m/max 0.0 (m/- 1.0 v)) 0.0))))

(defn ->bessel
  "Creates a BesselJ semivariogram model constructor for a given Bessel function `order`.

  The formula is: `nugget + psill * (1 - Γ(order+1) * (2 * range / x)^order * J_order(x / range))`

  where `Γ` is the gamma function and `J_order` is the Bessel function of the first kind. The factor
  `Γ(order+1) * (2 / d)^order * J_order(d)`, `d = x / range`, tends to `1.0` at the origin and is evaluated
  as a series for small lags (`d^2 < 4 * (order + 1)`), so the semivariance is finite for any order and lag (no
  `0 * Infinity`): it tends to `nugget` for tiny lags, also for orders far above 100. The semivariance
  oscillates around `nugget + psill` and tends to this value for an infinite lag.

  Parameters:

  - `order` (double, optional): the order of the Bessel function of the first kind;
    must be non-negative; defaults to `1.0`, the same default as [[bessel]] (and the order used by
    [[fit-params]] for `:bessel` when `:parameter` is not given)

  Returns a semivariogram model constructor `(fn [params])` that accepts a map with:
  - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
    variation or measurement error; must be non-negative
  - `:psill` (double): the partial sill; the amplitude of the oscillation; must be non-negative
  - `:range` (double): the scale parameter controlling the period of oscillation; must be positive

  The constructor itself returns a function `(fn ^double [^double x])` that computes the
  semivariance for a given lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[bessel]], [[->matern]], [[hole]]."
  ([] (->bessel 1.0))
  ([^double order]
   (check-nonneg! :bessel :order order)
   (let [g (special/gamma (m/inc order))]
     (fn [{:keys [^double nugget ^double psill ^double range]}]
       (check-sill-range! :bessel nugget psill range)
       (fn ^double [^double x]
         (if (m/zero? x) 0.0
             (->> (bessel-deficit order g (m// x range))
                  (m/* psill)
                  (m/+ nugget))))))))

(defn bessel
  "BesselJ semivariogram model.

  A convenience wrapper around [[->bessel]] that reads the Bessel function order directly
  from the `params` map. When `:order` is not provided, defaults to `1.0`.

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the amplitude of the oscillation; must be non-negative
      - `:range` (double): the scale parameter controlling the period of oscillation; must be positive
      - `:order` (double, optional): the order of the Bessel function of the first kind;
        must be non-negative; defaults to `1.0`

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[->bessel]], [[matern]], [[hole]]."
  [params]
  ((->bessel (or (:order params) 1.0)) params))

(defn ->matern
  "Creates a Matérn semivariogram model constructor for a given smoothness `order`.

  The formula is: `nugget + psill * (1 - C(order) * (sqrt(2*order) * x / range)^order * K_order(sqrt(2*order) * x / range))`

  where `C(order) = 2^(1-order) / Γ(order)` is a normalisation constant ensuring the
  kernel equals `1.0` at the origin, `Γ` is the gamma function and `K_order` is the
  modified Bessel function of the second kind of the given order.

  The scale is `sqrt(2 * order) * x / range`, so `range` is the length scale of the textbook form (Rasmussen and Williams), which differs from the gstat `Mat` and `Ste` ranges. `order = 0.5` is [[exponential]] with the same `range`. The limit of large `order` is a Gaussian model with the length scale `range`, that is [[gaussian]] with the range `range * sqrt(2)` (`exp(-(x / range)^2 / 2)` is `exp(-(x / (range * sqrt 2))^2)`).

  The product of the power and the Bessel function is evaluated without overflow: for `order >= 20` the factor is computed from a series (small lags) or a log-domain uniform expansion, so the semivariance is finite for any order and lag and tends to `nugget` for tiny lags. For orders below 20 a product which is not finite (lags below about `1e-30`) is replaced by its limit, so the semivariance is `nugget`.

  Parameters:

  - `order` (double, optional): the smoothness order of the Matérn model; must be positive;
    defaults to `1.0`

  Returns a semivariogram model constructor `(fn [params])` that accepts a map with:
  - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
    variation or measurement error; must be non-negative
  - `:psill` (double): the partial sill; the maximum contribution of the structured
    spatial variance above the nugget; must be non-negative
  - `:range` (double): the range parameter controlling the horizontal scale; must be positive

  The constructor itself returns a function `(fn ^double [^double x])` that computes the
  semivariance for a given lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[matern]], [[->bessel]], [[gaussian]], [[exponential]]."
  ([] (->matern 1.0))
  ([^double order]
   (check-pos! :matern :order order)
   (let [g (m// (m/pow 2.0 (m/- 1.0 order)) (special/gamma order))
         f (m/sqrt (m/* 2.0 order))]
     (fn [{:keys [^double nugget ^double psill ^double range]}]
       (check-sill-range! :matern nugget psill range)
       (fn ^double [^double x]
         (if (m/zero? x) 0.0
             (->> (matern-deficit order g (m/* f (m// x range)))
                  (m/* psill)
                  (m/+ nugget))))))))

(defn matern
  "Matérn semivariogram model.

  A convenience wrapper around [[->matern]] that reads the smoothness order directly from
  the `params` map. When `:order` is not provided, defaults to `1.0`.

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the maximum contribution of the structured
        spatial variance above the nugget; must be non-negative
      - `:range` (double): the range parameter controlling the horizontal scale; must be positive
      - `:order` (double, optional): the smoothness order; must be positive; defaults to `1.0`;
        `0.5` is equivalent to the [[exponential]] model; large values approach the [[gaussian]] model
        with the range `range * sqrt(2)`, not with `range`

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[->matern]], [[bessel]], [[gaussian]], [[exponential]]."
  [params]
  ((->matern (or (:order params) 1.0)) params))

;; https://gmd.copernicus.org/preprints/gmd-2021-301/gmd-2021-301.pdf
(defn ->hyperspherical
  "Creates a hyperspherical semivariogram model constructor for a given sphere dimension `order`.

  The formula is:
  - `x = 0`: `0.0`
  - `x <= range`: `nugget + psill * (h * 2F1(0.5, d, 1.5, h^2) / 2F1(0.5, d, 1.5, 1.0))`
  - `x > range`: `nugget + psill`

  where `h = x / range`, `d = -0.5 * (order - 1)` and `2F1` is the Gauss hypergeometric function.
  The denominator `2F1(0.5, d, 1.5, 1.0)` is a normalisation constant ensuring the model
  reaches exactly `psill` at `range`.

  Accuracy limits (measured against 80-digit mpmath): the evaluation of the hypergeometric function `2F1` loses digits for large orders. For orders giving a negative integer parameter `d` (odd dimensions here, integer orders for [[->superspherical]]) the absolute error is below `1e-12` up to `d = -15` (dimension 31, order 15), about `4e-12` at `d = -20` and `3e-9` at `d = -30`. For the other orders the error stays below `3e-14` up to the dimension 60 (order 29.5 for [[->superspherical]]) and the result is unreliable from the dimension 75 (error `6e-7`).

  Parameters:

  - `order` (double, optional): the sphere dimension; must be a positive real number;
    `2.0` recovers [[circular]], `3.0` recovers [[spherical]]; defaults to `2.0`

  Returns a semivariogram model constructor `(fn [params])` that accepts a map with:
  - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
    variation or measurement error; must be non-negative
  - `:psill` (double): the partial sill; the structured spatial variance above the nugget;
    must be non-negative
  - `:range` (double): the range at which the sill is reached exactly; must be positive

  The constructor itself returns a function `(fn ^double [^double x])` that computes the
  semivariance for a given lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[hyperspherical]], [[->superspherical]], [[spherical]], [[circular]]."
  ([] (->hyperspherical 2.0))
  ([^double order]
   (check-pos! :hyperspherical :order order)
   (let [d (m/* (m/dec order) -0.5)
         denom (special/hypergeometric-2F1 0.5 d 1.5 1.0)]
     (fn [{:keys [^double nugget ^double psill ^double range]}]
       (check-sill-range! :hyperspherical nugget psill range)
       (fn ^double [^double x]
         (cond
           (m/zero? x) 0.0
           (m/< range x) (m/+ nugget psill)
           :else (let [h (m// x range)]
                   (->> (m/* h (m// (special/hypergeometric-2F1 0.5 d 1.5 (m/* h h))
                                    denom))
                        (m/* psill)
                        (m/+ nugget)))))))))

(defn hyperspherical
  "Hyperspherical semivariogram model.

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the structured spatial variance above the nugget;
        must be non-negative
      - `:range` (double): the range at which the sill is reached exactly; must be positive
      - `:order` (double, optional): the sphere dimension; must be a positive real number;
        defaults to `2.0`; `2.0` is equivalent to [[circular]], `3.0` to [[spherical]]

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[->hyperspherical]], [[superspherical]], [[spherical]], [[circular]]."
  [params]
  ((->hyperspherical (or (:order params) 2.0)) params))

(defn ->superspherical
  "Creates a superspherical semivariogram model constructor for a given `order`.

  The superspherical model is a bounded semivariogram from the same family as
  [[->hyperspherical]], using the Gauss hypergeometric function `2F1` with `d = -order`
  instead of `d = -0.5 * (order - 1)`. It is a compact reparametrisation of the
  hyperspherical family: superspherical `order` corresponds to hyperspherical dimension
  `2 * order + 1`. Notable special cases are `order = 0.0`, which recovers the [[linear]]
  model, `order = 0.5`, which recovers the [[circular]] model, and `order = 1.0` (default),
  which recovers the standard [[spherical]] model. Like all members of the spherical family,
  it reaches the sill (`nugget + psill`) exactly at `range` and remains constant beyond it.

  The formula is:
  - `x = 0`: `0.0`
  - `x <= range`: `nugget + psill * (h * 2F1(0.5, -order, 1.5, h^2) / 2F1(0.5, -order, 1.5, 1.0))`
  - `x > range`: `nugget + psill`

  where `h = x / range` and `2F1` is the Gauss hypergeometric function.
  The denominator `2F1(0.5, -order, 1.5, 1.0)` is a normalisation constant ensuring the
  model reaches exactly `psill` at `range`.

  Accuracy limits (measured against 80-digit mpmath): the evaluation of the hypergeometric function `2F1` loses digits for large orders. For orders giving a negative integer parameter `d` (odd dimensions here, integer orders for [[->superspherical]]) the absolute error is below `1e-12` up to `d = -15` (dimension 31, order 15), about `4e-12` at `d = -20` and `3e-9` at `d = -30`. For the other orders the error stays below `3e-14` up to the dimension 60 (order 29.5 for [[->superspherical]]) and the result is unreliable from the dimension 75 (error `6e-7`).

  Parameters:

  - `order` (double, optional): the shape order; must be non-negative; `0.0` gives [[linear]],
    `0.5` gives [[circular]], `1.0` gives [[spherical]]; defaults to `1.0`

  Returns a semivariogram model constructor `(fn [params])` that accepts a map with:
  - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
    variation or measurement error; must be non-negative
  - `:psill` (double): the partial sill; the structured spatial variance above the nugget;
    must be non-negative
  - `:range` (double): the range at which the sill is reached exactly; must be positive

  The constructor itself returns a function `(fn ^double [^double x])` that computes the
  semivariance for a given lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[superspherical]], [[->hyperspherical]], [[spherical]], [[circular]], [[linear]]."
  ([] (->superspherical 1.0))
  ([^double order]
   (check-nonneg! :superspherical :order order)
   (let [d (m/- order)
         denom (special/hypergeometric-2F1 0.5 d 1.5 1.0)]
     (fn [{:keys [^double nugget ^double psill ^double range]}]
       (check-sill-range! :superspherical nugget psill range)
       (fn ^double [^double x]
         (cond
           (m/zero? x) 0.0
           (m/< range x) (m/+ nugget psill)
           :else (let [h (m// x range)]
                   (->> (m/* h (m// (special/hypergeometric-2F1 0.5 d 1.5 (m/* h h))
                                    denom))
                        (m/* psill)
                        (m/+ nugget)))))))))

(defn superspherical
  "Superspherical semivariogram model.

  A convenience wrapper around [[->superspherical]] that reads the shape order directly
  from the `params` map. When `:order` is not provided, defaults to `1.0`, which is
  equivalent to the standard [[spherical]] model. Setting `:order` to `0.0` recovers
  [[linear]] and `0.5` recovers [[circular]].

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the structured spatial variance above the nugget;
        must be non-negative
      - `:range` (double): the range at which the sill is reached exactly; must be positive
      - `:order` (double, optional): the shape order; must be non-negative; defaults to `1.0`;
        `0.0` is equivalent to [[linear]], `0.5` to [[circular]], `1.0` to [[spherical]]

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[->superspherical]], [[hyperspherical]], [[spherical]], [[circular]], [[linear]]."
  [params]
  ((->superspherical (or (:order params) 1.0)) params))

(defn ->tplstable
  "Creates a Truncated-Power-Law (TPL) Stable semivariogram model constructor for a given Hurst exponent `H`.

  The TPL Stable model is a bounded (asymptotic) semivariogram that captures both
  short-range power-law growth and long-range correlation decay. The Hurst exponent `H`
  controls the fractal dimension and long-range dependence of the underlying spatial
  process; values closer to `1.0` indicate stronger long-range persistence, while values
  closer to `0.0` indicate stronger short-range variability. The `beta` parameter (provided
  at construction time via `params`) shapes the power-law growth rate similarly to [[power]]
  and [[expower]]. The semivariance asymptotically approaches the sill (`nugget + psill`)
  as lag distance increases.

  The formula is: `nugget + psill * (1 - Ha * E_{1+Ha}(h^beta))`

  where `h = x / range`, `Ha = 2H / beta` and `E_n` is the generalised exponential integral.
  The factor `Ha * E_{1+Ha}(0) = 1` ensures the semivariance is `0.0` at the origin.

  Parameters:

  - `H` (double, optional): the Hurst exponent controlling long-range dependence;
    must be in `(0.0, 1.0)`; defaults to `0.5`

  Returns a semivariogram model constructor `(fn [params])` that accepts a map with:
  - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
    variation or measurement error; must be non-negative
  - `:psill` (double): the partial sill; the maximum contribution of the structured
    spatial variance above the nugget; must be non-negative
  - `:range` (double): the scale parameter controlling the horizontal extent; must be positive
  - `:beta` (double): the power exponent controlling the short-range growth rate;
    must be in `(0.0, 2.0]` for a valid conditionally negative-definite model

  The constructor itself returns a function `(fn ^double [^double x])` that computes the
  semivariance for a given lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[tplstable]], [[power]], [[expower]], [[matern]]."
  ([] (->tplstable 0.5))
  ([^double H]
   (check-hurst! :tplstable H)
   (let [H2 (m/* 2.0 H)]
     (fn [{:keys [^double nugget ^double psill ^double range ^double beta]}]
       (check-sill-range! :tplstable nugget psill range)
       (check-beta! :tplstable beta 2.0)
       (fn ^double [^double x]
         (if (m/zero? x) 0.0
             (let [h (m// x range)
                   Ha (m// H2 beta)]
               (->> (m/pow h beta)
                    (special/En (m/inc Ha))
                    (m/* Ha)
                    (m/- 1.0)
                    (m/* psill)
                    (m/+ nugget)))))))))

(defn tplstable
  "Truncated-Power-Law (TPL) Stable semivariogram model.

  A convenience wrapper around [[->tplstable]] that reads the Hurst exponent directly
  from the `params` map via the `:H` key. When `:H` is not provided, defaults to `0.5`.

  Parameters:

  - `params` (map): a map of model parameters with the following keys:
      - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
        variation or measurement error; must be non-negative
      - `:psill` (double): the partial sill; the maximum contribution of the structured
        spatial variance above the nugget; must be non-negative
      - `:range` (double): the scale parameter controlling the horizontal extent; must be positive
      - `:beta` (double): the power exponent controlling the short-range growth rate;
        must be in `(0.0, 2.0]` for a valid conditionally negative-definite model
      - `:H` (double, optional): the Hurst exponent controlling long-range dependence;
        must be in `(0.0, 1.0)`; defaults to `0.5`

  Returns a function `(fn ^double [^double x])` that computes the semivariance for a given
  lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[->tplstable]], [[power]], [[expower]], [[matern]]."
  [params]
  ((->tplstable (or (:H params) 0.5)) params))

(defn rbf->variogram
  "Builds a semivariogram model constructor from a Radial Basis Function (RBF) kernel.

  Converts any RBF kernel from `fastmath.kernel.rbf` into a semivariogram model by
  treating the kernel as the correlation function: semivariance equals
  `nugget + psill * (1 - kernel(x / range))`. For this transformation to produce a valid
  semivariogram, the kernel must satisfy two conditions: it must equal `1.0` at the origin
  and decay to `0.0` at infinity. The kernel is evaluated at `x / range`, so the `range` of the model
  multiplies the `scale` of the RBF kernel: the kernel with `:scale` 2 and the `range` 3 gives the same
  model as the kernel with `:scale` 1 and the `range` 6 (for the gaussian kernel it is [[gaussian]] with
  the range 6).

  The formula is: `nugget + psill * (1 - kernel(x / range))`

  Suitable kernels from `fastmath.kernel.rbf` (value `1.0` at zero, decaying to `0.0`):

  - `:gaussian` - smooth Gaussian decay
  - `:matern` - Matérn family; differentiability controlled by order
  - `:truncated-power` - compactly supported power kernel
  - `:wendland` - compactly supported, C^(2k) smooth
  - `:gneiting` - compactly supported Gneiting family
  - `:wu` - compactly supported Wu family
  - `:whittaker` - Whittaker family, normalised to `1.0` at zero
  - `:poisson` - Bessel J family; oscillates around `0.0`, use with care; the kernel has to be normalised to `1.0` at zero (for dimensions above 2 it is only after the normalisation of the kernel in `fastmath.kernel.rbf`)

  Unsuitable kernels (do not satisfy the required conditions):

  - `:thin-plate` / `:thin-plate-splines` - grows without bound
  - `:linear` - zero at origin, grows without bound
  - `:radial-powers` - zero or infinite at origin depending on exponent
  - `:generalized-multiquadratic` - not normalised to `+1.0` at zero
  - `:shifted-surface-splines` - grows with distance for positive exponent

  Parameters:

  - `kernel` (function): an RBF kernel function `(fn ^double [^double x])` returning `1.0`
    at `x = 0` and decaying to `0.0` as `x` increases; typically constructed via
    `fastmath.kernel.rbf` functions. The value of the kernel at zero is not used (the semivariance
    at lag zero is `0.0`), so a kernel which is not `1.0` at zero gives a model which is not continuous
    as the lag goes to zero from `nugget`.

  Returns a semivariogram model constructor `(fn [params])` that accepts a map with:
  - `:nugget` (double): the nugget effect; the y-intercept representing micro-scale
    variation or measurement error; must be non-negative
  - `:psill` (double): the partial sill; the maximum contribution of the structured
    spatial variance above the nugget; must be non-negative
  - `:range` (double): the horizontal scale, multiplies the RBF `scale`; must be positive

  The constructor itself returns a function `(fn ^double [^double x])` that computes the
  semivariance for a given lag distance `x`. Returns `0.0` when `x` is zero.

  The parameters are checked once, when the model is built (not at every evaluation): a missing, non-numeric or non-finite parameter, or a value outside the stated range, throws `ex-info` with `:type :invalid-input` and the keys `:model`, `:property`, `:value` and `:expected` in the data. The lag `x` is expected to be non-negative; `NaN` gives `NaN`.

  See also [[exponential]], [[gaussian]], [[matern]], [[fit-params]]."
  [kernel]
  (fn [{:keys [^double nugget ^double psill ^double range]}]
    (check-sill-range! :rbf nugget psill range)
    (fn ^double [^double x]
      (if (m/zero? x) 0.0
          (->> (m// x range)
               (kernel)
               (double)
               (m/- 1.0)
               (m/* psill)
               (m/+ nugget))))))

(defn variogram->covariance
  "Converts a semivariogram function to the covariance function with the given total sill.

  The covariance is `C(x) = sill - gamma(x)`, the relation of a second order stationary process. The total sill is the variance of the process, `nugget + psill` of the models of this namespace: then `C(0)` equals the variance, `C(x)` tends to `psill` (the nugget is lost) for a tiny positive lag and `C(x) = 0` beyond the `range` of the bounded models.

  Parameters:

  - `semivariogram` (function): a semivariogram `(fn ^double [^double x])`, for example a model of this namespace built from a parameter map or the result of [[fit]]; it should return `0.0` at lag zero, which all the models do
  - `sill` (double): the total sill (the variance), it must be finite and non-negative; use `(+ nugget psill)` of the parameters of the model, or from the result of [[fit-params]] `(let [{:keys [nugget psill]} (:params result)] (+ nugget psill))`

  Returns a function `(fn ^double [^double x])` computing `sill - gamma(x)`; `NaN` gives `NaN`.

  Corner cases: for the models without a sill (see [[power]]) the covariance is negative from some lag and is not a covariance of a stationary process; an unbounded semivariogram makes only sense as the generalised covariance. Throws `ex-info` with `:type :invalid-input` and the key `:property` equal `:sill` for a `sill` which is not a finite non-negative number, and with the property `:semivariogram` when the first argument is not a function.

  See also [[fit-params]], [[exponential]], [[spherical]], [[gaussian]]."
  [semivariogram sill]
  (when-not (ifn? semivariogram)
    (throw (ex-info "The semivariogram has to be a function."
                    {:type :invalid-input :property :semivariogram :value semivariogram})))
  (when-not (and (number? sill) (Double/isFinite (double sill)) (m/not-neg? (double sill)))
    (throw (ex-info "The sill has to be a finite number >= 0."
                    {:type :invalid-input :property :sill :value sill})))
  (let [sill (double sill)]
    (fn ^double [^double x]
      (m/- sill (double (semivariogram x))))))

(def ^:private semivariograms {:linear linear
                               :pentaspherical pentaspherical
                               :spherical spherical
                               :gaussian gaussian
                               :exponential exponential
                               :power power
                               :tpower tpower
                               :expower expower
                               :hole hole
                               :circular circular
                               :cubic cubic
                               :cauchy cauchy
                               :rational rational
                               :bessel bessel
                               :matern matern
                               :superspherical superspherical
                               :hyperspherical hyperspherical
                               :tplstable tplstable})

(def ^:private power-semivariograms #{:power :expower :cauchy :rational :tpower :tplstable})

(def ^:private ->semivariograms {:bessel ->bessel
                               :matern ->matern
                               :superspherical ->superspherical
                               :hyperspherical ->hyperspherical
                               :tplstable ->tplstable})

;; Highly Robust Variogram Estimation
;; Marc G. Genton

;; Basic Steps in Geostatistics: The Variogram and Kriging
;; Margaret A. Oliver, Richard Webster

(defn matheron-estimator
  "Matheron (classical) empirical semivariogram estimator.

  Parameters:

  - `n` (long): number of pairs in the bin (not used)
  - `diffs` (sequence of doubles): differences of the values of the pairs in the bin

  Returns the mean of the squared differences, the estimate of the variogram (twice the semivariance; [[empirical]] halves it). `NaN` for an empty `diffs`.

  See also [[empirical]], [[cressie-estimator]], [[dowd-estimator]]."
  ^double [_ diffs]
  (stats/mean (map m/sq diffs)))

(defn cressie-estimator
  "Cressie and Hawkins (robust) empirical semivariogram estimator.

  The estimate of the variogram is `(mean(sqrt(abs(d))))^4 / (0.457 + 0.494 / N + 0.045 / N^2)`, where `d` are the differences of the pair values and `N` is the number of pairs. The constant is the one of Hawkins and Cressie (1984), which removes the bias of the fourth power for the Gaussian data; it is the one of scikit-gstat. The gstat and geoR packages omit the last term `0.045 / N^2`: their value is larger by 4.5% for `N` equal 1 and by 7.4e-4 for `N` equal 11 (the difference vanishes for large bins).

  Parameters:

  - `n` (long): number of pairs in the bin
  - `diffs` (sequence of doubles): differences of the values of the pairs in the bin

  Returns the estimate of the variogram (twice the semivariance; [[empirical]] halves it).

  See also [[empirical]], [[matheron-estimator]], [[highly-robust-estimator]], [[dowd-estimator]]."
  [^long n diffs]
  (-> (->> (map (fn [^double diff] (m/sqrt (m/abs diff))) diffs)
           (stats/mean))
      (m/fpow 4)
      (m// (m/+ 0.457 (m// 0.494 n) (m// 0.045 (m/* n n))))))

(defn highly-robust-estimator
  "Genton's highly robust empirical semivariogram estimator.

  The estimate of the variogram is `(2.2191 * Q)^2`, where `Q` is the quantile of the absolute differences `|d_i - d_j|` (`i < j`) of all pairs of the differences in the bin, at the level `k / C(N, 2)` with `k = C(floor(N / 2) + 1, 2)` and `N` the number of pairs. The quantile is the interpolated sample quantile of [[fastmath.stats/quantile]], not the exact k-th order statistic of the original estimator, so the result differs slightly from it (relative differences of a few percent for small bins, e.g. 4.8e-2 in the research data).

  Parameters:

  - `n` (long): number of pairs in the bin
  - `diffs` (sequence of doubles): differences of the values of the pairs in the bin

  Returns the estimate of the variogram (twice the semivariance; [[empirical]] halves it). The estimator is not defined for a bin with a single pair (`N` equal 1: no pair of pairs) and returns `NaN` there; [[empirical]] drops such bins. The cost is quadratic in the size of a bin.

  See also [[empirical]], [[cressie-estimator]], [[dowd-estimator]], [[robust-m-estimator]]."
  [^long n diffs]
  (let [idiffs (map-indexed vector diffs)
        ds (for [[^long i1 ^double d1] idiffs
                 [^long i2 ^double d2] idiffs
                 :when (m/< i1 i2)]
             (m/abs (m/- d1 d2)))
        q (m// (m/combinations (inc (int (m/* 0.5 n))) 2)
               (count ds))]
    (m/sq (m/* 2.2191 (stats/quantile ds q)))))

(defn dowd-estimator
  "Dowd's empirical semivariogram estimator.

  The estimate of the variogram is `2.198109... * (median(abs(d)))^2`, where `d` are the differences of the pair values (Dowd 1984); the constant is the reciprocal of the median of the chi-squared distribution with one degree of freedom, so the estimate is consistent for Gaussian data. For an even number of pairs the median of the absolute differences is averaged before squaring; this is not the same as the median of the squared differences, which is what the quantile estimator [[->quantile-estimator]] with `p` equal `0.5` computes.

  Parameters:

  - `n` (long): number of pairs in the bin (not used)
  - `diffs` (sequence of doubles): differences of the values of the pairs in the bin

  Returns the estimate of the variogram (twice the semivariance; [[empirical]] halves it). `NaN` for an empty `diffs`.

  See also [[empirical]], [[->quantile-estimator]], [[cressie-estimator]], [[matheron-estimator]]."
  [_ diffs]
  (m/* 2.198109338317728 (m/sq (stats/median (map (fn [^double d] (m/abs d)) diffs)))))

(defn ->quantile-estimator
  "Creates a quantile (Armstrong and Delfiner) empirical semivariogram estimator for given `p`.

  The estimate of the variogram is the `p`-quantile of the squared differences divided by the `p`-quantile of the chi-squared distribution with one degree of freedom.

  Parameters:

  - `p` (double): the quantile level, must be in `(0.0, 1.0)`

  Returns a function `(fn [n diffs])` (the arguments are the number of pairs in the bin, not used, and the differences of the pair values) which returns the estimate of the variogram (twice the semivariance; [[empirical]] halves it). Throws `ex-info` for `p` outside `(0.0, 1.0)`.

  See also [[empirical]], [[dowd-estimator]]."
  [^double p]
  (when-not (and (m/> p 0.0) (m/< p 1.0))
    (throw (ex-info "The quantile level has to be in (0, 1)."
                    {:type :invalid-input :property :quantile :value p :expected "a number in (0, 1)"})))
  (let [fact (m// (double (r/icdf (r/distribution :chi-squared) p)))]
    (fn [_ diffs]
      (m/* fact (stats/quantile (map m/sq diffs) p)))))

(defn- psi
  [^double t]
  (if (m/>= (m/abs t) 4.0)
    0.0
    (m/* t (m/sq (m/- 16.0 (m/sq t))))))

(defn- dpsi
  [^double t]
  (if (m/>= (m/abs t) 4.0)
    0.0
    (let [t2 (m/sq t)]
      (m/+ (m/* 5.0 (m/sq t2))
           (m/* -96.0 t2)
           256.0))))

;; ROBUST SEMIVARIOGRAM ESTIMATION 
;; IN THE PRESENCE OF INFLUENTIAL 
;; SPATIAL DATA VALUES
;; Richard F. Gunst and Molly I. Hartfield 

(defn robust-m-estimator
  "Robust M-estimator (Gunst and Hartfield) for empirical semivariogram.

  A biweight (Tukey) M-estimate of location of the square roots of the absolute differences, raised to the fourth power and divided by the Hawkins and Cressie constant `0.457 + 0.494 / N + 0.045 / N^2` (see [[cressie-estimator]]).

  Parameters:

  - `n` (long): number of pairs in the bin
  - `diffs` (sequence of doubles): differences of the values of the pairs in the bin

  Returns the estimate of the variogram (twice the semivariance; [[empirical]] halves it). The scale of the estimator is the median absolute deviation, which is `0.0` (the result is `NaN`) for a bin with a single pair or with identical absolute differences; [[empirical]] drops such bins. The implementation is checked only against an independent re-implementation, not against a reference package.

  See also [[empirical]], [[cressie-estimator]], [[highly-robust-estimator]]."
  [^long n diffs]
  (let [adiffs (map (fn [^double d] (m/sqrt (m/abs d))) diffs)
        m (stats/median adiffs)
        s (m// (stats/median-absolute-deviation adiffs m) 0.6745)
        y (map (fn [^double az] (m// (m/- az m) s)) adiffs)
        spsi (stats/sum (map psi y))
        sdpsi (stats/sum (map dpsi y))]
    (m// (m/fpow (m/+ m (m/* s (m// spsi sdpsi))) 4)
         (m/+ 0.457 (m// 0.494 n) (m// 0.045 (m/* n n))))))

(defn bounding-box-diagonal
  "Length of the diagonal of bounding box of spatial points.

  Parameters:

  - `dist` (function): a binary distance function `(fn [x1 x2])`
  - `xs` (sequence): the points, scalars or vectors; at least one point

  Returns the distance between the points made of the minimum and the maximum coordinates. It is `0.0` for a single point.

  See also [[empirical]]."
  ^double [dist xs]
  (let [mn (reduce v/emn xs)
        mx (reduce v/emx xs)]
    (dist mn mx)))

(defn remove-outliers-fence
  "Removes outliers using Tukey's fences (`k`=1.5).

  Parameters:

  - `combined` (sequence of `[id x y]` triples): the points, as built by [[empirical]]
  - `ys` (sequence of doubles): the values `y` of the triples

  Returns the triples with the value between `Q1 - 1.5 * IQR` and `Q3 + 1.5 * IQR` of `ys` (inclusive).

  See also [[maybe-remove-outliers]], [[remove-outliers-mad]]."
  [combined ys]
  (let [[^double q1 ^double q3] (stats/quantiles ys [0.25 0.75])
        iqr (m/* 1.5 (m/- q3 q1))
        lif-thr (m/- q1 iqr)
        uif-thr (m/+ q3 iqr)]
    (filter (fn [v] (m/<= lif-thr (double (v 2)) uif-thr)) combined)))

(defn remove-outliers-mad
  "Removes outliers using median absolute deviation.

  Parameters:

  - `combined` (sequence of `[id x y]` triples): the points, as built by [[empirical]]
  - `ys` (sequence of doubles): the values `y` of the triples

  Returns the triples for which the modified z-score `|y - median| / (MAD / 0.6745)` is at most `3.0`. When the median absolute deviation is `0.0` all values different from the median are removed.

  See also [[maybe-remove-outliers]], [[remove-outliers-fence]]."
  [combined ys]
  (let [m (stats/median ys)
        s (m// (stats/median-absolute-deviation ys m) 0.6745)]
    (filter (fn [v] (m/<= (m/abs (m// (m/- (double (v 2)) m) s)) 3.0)) combined)))

(defn maybe-remove-outliers
  "Optionally removes outliers from the points before the pairs are formed in [[empirical]].

  The outliers are detected in the observed values of the points, not in the differences of the values of the pairs.

  Parameters:

  - `combined` (sequence of `[id x y]` triples): the points
  - `ys` (sequence of doubles): the values `y` of the triples
  - `remove-outliers?`: `:mad` for [[remove-outliers-mad]], any other truthy value for [[remove-outliers-fence]], `false` or `nil` to keep all the points

  Returns the sequence of the kept triples.

  See also [[empirical]], [[remove-outliers-fence]], [[remove-outliers-mad]]."
  [combined ys remove-outliers?]
  (cond
    (= :mad remove-outliers?) (remove-outliers-mad combined ys)
    remove-outliers? (remove-outliers-fence combined ys)
    :else combined))

(def ^:private estimators {:classical matheron-estimator
                           :matheron matheron-estimator
                           :cressie cressie-estimator
                           :genton highly-robust-estimator
                           :highly-robust highly-robust-estimator
                           :dowd dowd-estimator
                           :m-robust robust-m-estimator})

(defn- invalid-input
  [message property value expected]
  (throw (ex-info message {:type :invalid-input :property property :value value :expected expected})))

(defn- resolve-estimator
  [estimator quantile]
  (cond
    (= :quantile estimator) (->quantile-estimator quantile)
    (keyword? estimator) (or (estimators estimator)
                             (invalid-input (str "Unknown empirical estimator " estimator)
                                            :estimator estimator
                                            (vec (sort (conj (keys estimators) :quantile)))))
    (ifn? estimator) estimator
    :else (invalid-input (str "Unknown empirical estimator " (pr-str estimator))
                         :estimator estimator
                         (vec (sort (conj (keys estimators) :quantile))))))

(defn- finite-seq?
  [diffs]
  (every? (fn [^double d] (Double/isFinite d)) diffs))

(defn empirical
  "Computes the empirical (experimental) semivariogram from spatial data.

  All pairs of spatial locations are formed, filtered by the cutoff distance, and grouped
  into lag bins. For each bin, a chosen estimator computes the gamma value, the classical one is
  half the average squared difference between the paired observations. The result is the primary
  input to [[fit]] and [[fit-params]] for fitting a parametric semivariogram model.

  The bins divide the interval `(0, dmax]` into `:bins` equal parts, where `dmax` is the largest distance among the pairs kept by the cutoff (not the cutoff). A bin includes its right edge and excludes the left one. This is the same as the explicit `boundaries` of the gstat `variogram` function with the same edges; the classical estimator then equals gstat (counts, mean lags and gamma), but the default bins of gstat (the cutoff divided by 15) are different.

  Parameters:

  - `xs` (sequence): spatial positions of observations; each element may be a scalar or
    a vector depending on the dimensionality of the data
  - `ys` (sequence): observed values at positions `xs`; must have the same length as `xs`
  - `params` (map, optional): configuration map with the following keys (a `nil` value is the same as a missing key):
      - `:bins` (long): number of lag bins; controls the resolution of the semivariogram;
        defaults to `15`; must be at least `1`
      - `:cutoff` (double): maximum lag distance to include; pairs beyond this distance are
        excluded; defaults to the bounding box diagonal divided by `:diagonal-den`
      - `:diagonal-den` (double): divisor applied to the bounding box diagonal to derive
        the default cutoff; defaults to `3.0`; ignored when `:cutoff` is provided explicitly
      - `:estimator` (keyword or function): the estimator used to compute gamma per bin;
        defaults to `:classical`; built-in options are:
          - `:classical` / `:matheron` - Matheron classical estimator ([[matheron-estimator]])
          - `:cressie` - Cressie robust estimator ([[cressie-estimator]])
          - `:genton` / `:highly-robust` - Genton highly robust estimator ([[highly-robust-estimator]])
          - `:dowd` - Dowd estimator ([[dowd-estimator]])
          - `:quantile` - Armstrong and Delfiner quantile estimator at level `:quantile` ([[->quantile-estimator]])
          - `:m-robust` - Gunst and Hartfield robust M-estimator ([[robust-m-estimator]])
          - a custom `(fn [^long n diffs])` where `n` is the number of pairs in the bin and `diffs` is a
            sequence of pairwise differences; it returns the estimate of the variogram (twice the semivariance)
      - `:quantile` (double): quantile level used by the `:quantile` estimator;
        must be in `(0.0, 1.0)`; defaults to `0.9`
      - `:remove-outliers?`: whether to remove outliers from the raw data before binning;
        defaults to `false`; accepted values are:
          - `false` - no outlier removal
          - `true` - Tukey's fences criterion (IQR-based, k=1.5) via [[remove-outliers-fence]]
          - `:mad` - median absolute deviation criterion via [[remove-outliers-mad]]
        the outliers are searched in the values `ys`, not in the pair differences, see [[maybe-remove-outliers]]
      - `:distance` (function): a binary distance function `(fn [x1 x2])` applied to
        elements of `xs`; defaults to `fastmath.distance/euclidean`

  Returns a vector of maps sorted by ascending lag `h`, one per non-empty bin, each containing:
  - `:n` (long): number of point pairs in the bin (not the number of points)
  - `:h` (double): mean lag distance of pairs in the bin
  - `:gamma` (double): estimated semivariance for the bin

  Pairs of points at distance zero (duplicated locations) are dropped, as pairs with `NaN` distance. A bin whose estimate is not finite, although all the differences of its pairs are finite, is dropped: the robust estimators [[highly-robust-estimator]] and [[robust-m-estimator]] are not defined for a bin with a single pair. A `NaN` in `ys` is not dropped: it gives a `NaN` gamma (which [[fit-params]] rejects).

  Throws `ex-info` (with `:type :invalid-input` and the key `:property`) when `xs` and `ys` have different lengths, there are fewer than two points, `:bins` is below `1`, the estimator is unknown, no pair is left within the cutoff (all the points at one location, the cutoff too small, too many outliers removed) or no bin has a finite estimate.

  See also [[fit]], [[fit-params]], [[bounding-box-diagonal]],
  [[matheron-estimator]], [[cressie-estimator]], [[highly-robust-estimator]],
  [[dowd-estimator]], [[->quantile-estimator]], [[robust-m-estimator]]."
  ([xs ys] (empirical xs ys nil))
  ([xs ys params]
   (let [{:keys [cutoff diagonal-den bins estimator quantile remove-outliers? distance]
          :or {bins 15 diagonal-den 3.0 estimator :classical quantile 0.9
               remove-outliers? false distance dist/euclidean}}
         (into {} (remove (fn [[_ v]] (nil? v))) params)
         nxs (count xs)
         _ (when-not (== nxs (count ys))
             (invalid-input (str "xs and ys have different lengths: " nxs " and " (count ys))
                            :ys (count ys) (str "a sequence of the length " nxs)))
         _ (when (< nxs 2)
             (invalid-input "At least two points are needed to form a pair." :xs nxs "at least 2 points"))
         bins (long bins)
         _ (when (< bins 1)
             (invalid-input "The number of bins has to be at least 1." :bins bins "a number >= 1"))
         estimator-fn (resolve-estimator estimator quantile)
         cutoff (double (or cutoff (/ (bounding-box-diagonal distance xs) (double diagonal-den))))
         combined (-> (map vector (range) xs ys)
                      (maybe-remove-outliers ys remove-outliers?))
         distances (->> (for [[^long id1 x1 ^double y1] combined
                              [^long id2 x2 ^double y2] combined
                              :when (< id1 id2)
                              :let [d (double (distance x1 x2))]
                              :when (and (pos? d) (<= d cutoff))]
                          (Vec2. d (m/- y1 y2)))
                        (sort-by first))
         _ (when-not (seq distances)
             (invalid-input "No pair of points with a positive distance within the cutoff."
                            :cutoff cutoff "a cutoff which keeps at least one pair of distinct points"))
         max-dist (.x ^Vec2 (last distances))
         ;; the last bin has no upper limit, so the largest distance never falls out
         uppers (concat (butlast (rest (m/slice-range 0 max-dist (inc bins)))) [##Inf])
         result (loop [buff []
                       distances distances
                       uppers uppers]
                  (if (or (empty? distances) (empty? uppers))
                    buff
                    (let [c (double (first uppers))
                          [found others] (split-with (fn [^Vec2 v] (m/<= (.x v) c)) distances)]
                      (if (empty? found)
                        (recur buff distances (rest uppers))
                        (let [n (count found)
                              diffs (map second found)
                              gamma (m/* 0.5 (double (estimator-fn n diffs)))]
                          (recur (if (or (Double/isFinite gamma) (not (finite-seq? diffs)))
                                   (conj buff {:n n
                                               :h (stats/mean (map first found))
                                               :gamma gamma})
                                   buff)
                                 others
                                 (rest uppers)))))))]
     (when-not (seq result)
       (invalid-input "No bin has a finite estimate of the semivariance (the estimator is not defined for the bins)."
                      :estimator estimator "an estimator which is defined for bins of the given size"))
     result)))

(defn- weights-n-v1  [ns vs] (map (fn [^double n ^double v] (m// n (m/max m/EPSILON (m/abs v)))) ns vs))
(defn- weights-n-v2  [ns vs] (map (fn [^double n ^double v] (m// n (m/max m/EPSILON (m/sq v)))) ns vs))
(defn- weights-n-v3  [ns vs] (map (fn [^double n ^double v] (m// n (m/max m/EPSILON (m/cb v)))) ns vs))
(defn- est-sq  ^double [^double a ^double b] (m/sq (m/- a b)))
(defn- est-abs ^double [^double a ^double b] (m/abs (m/- a b)))
(defn- est-soft ^double [^double a ^double b] (m/* 2.0 (m/dec (m/sqrt (m/inc (est-sq a b))))))
(defn- est-huber ^double [^double a ^double b] (let [d (est-sq a b)]
                                                 (if (m/<= d 1.0) d (m/dec (m/* 2.0 (est-abs a b))))))
(defn- est-cauchy ^double [^double a ^double b] (m/log (m/inc (est-sq a b))))
(defn- est-atan ^double [^double a ^double b] (m/atan (est-sq a b)))

(def ^:private estimation-fns {:sq est-sq
                               :abs est-abs
                               :soft est-soft
                               :huber est-huber
                               :cauchy est-cauchy
                               :atan est-atan})

(def ^:private weight-schemes [:n :nh :nh+ :nhh :nhh+ :nexp :nsqrt :ngg :ngg2])

(def ^:private parameter-names [:nugget :psill :range :beta])

;; valid intervals of beta used as default search bounds (a model checks them when it is built)
(def ^:private beta-bounds {:power [m/EPSILON 2.0]
                            :expower [m/EPSILON 2.0]
                            :tplstable [m/EPSILON 2.0]
                            :tpower [m/EPSILON 30.0]
                            :rational [m/EPSILON 30.0]
                            :cauchy [m/EPSILON 30.0]})

(defn- infer-target-args
  [defaults params]
  (->> params
       (reduce (fn [buff t]
                 (if-not (t defaults)
                   (conj buff t)
                   buff)) [])))

(defn- check-parameter-keys!
  [property m]
  (when-let [unknown (seq (remove (set parameter-names) (keys m)))]
    (invalid-input (str "Unknown keys in " property ": " (vec unknown))
                   property (vec unknown) (str "keys from " parameter-names))))

(defn- check-empirical!
  [empirical-semivariogram]
  (when-not (seq empirical-semivariogram)
    (invalid-input "The empirical semivariogram is empty." :empirical-semivariogram empirical-semivariogram
                   "a non-empty sequence of maps with :h, :gamma and :n"))
  (doseq [{:keys [h gamma n] :as bin} empirical-semivariogram
          [property value] [[:h h] [:gamma gamma] [:n n]]]
    (when-not (and (number? value) (Double/isFinite (double value)))
      (invalid-input (str "The empirical semivariogram has a non-finite or missing " property " in the bin " (pr-str bin))
                     property value "a finite number in every bin"))))

(defn- weights-sequence
  [weights ns hs hs+]
  (let [nbins (long (count ns))]
    (cond
      (not weights) (repeat nbins 1.0)
      (sequential? weights) (if (and (== (count weights) nbins)
                                     (every? (fn [w] (and (number? w) (Double/isFinite (double w)))) weights))
                              weights
                              (invalid-input (str "Explicit weights have to be " nbins " finite numbers, one for every bin")
                                             :weights weights (str "a sequence of " nbins " finite numbers")))
      (fn? weights) (weights-n-v1 ns (map weights hs))
      :else (case weights
              :n ns
              :nhh (weights-n-v2 ns hs)
              :nhh+ (weights-n-v2 ns hs+)
              :nh (weights-n-v1 ns hs)
              ;; gstat python
              :nh+ (map m// ns hs+)
              :nexp (weights-n-v1 ns (v/exp (v/reciprocal hs+)))
              :nsqrt (weights-n-v1 ns (v/sqrt hs+))
              (:ngg :ngg2) weights
              (invalid-input (str "Unknown weights " (pr-str weights)) :weights weights
                             (str "nil, a sequence, a function or one of " weight-schemes))))))

(defn fit-params
  "Fits a parametric semivariogram model to an empirical semivariogram and returns a detailed result map.

  Minimises a weighted loss between the model predictions and the empirical gamma values
  using L-BFGS-B numerical optimisation with jitter-based multi-start via
  `fastmath.optimization/scan-and-minimize`. Parameters fixed via `:defaults` are excluded
  from the optimisation. When all parameters are fixed, no optimisation is performed.
  Returns a richer result than [[fit]], including the fitted parameters, weights, and final
  loss value. Useful for model diagnostics and manual inspection of the fit quality.

  Parameters:

  - `empirical-semivariogram` (sequence of maps): empirical semivariogram data as returned
    by [[empirical]]; each map must contain `:h`, `:gamma` and `:n` keys, all finite
  - `semivariogram-model` (keyword or function): the parametric model to fit; accepted values are:
      - a keyword naming a built-in model: `:linear`, `:spherical`, `:circular`,
        `:cubic`, `:pentaspherical`, `:gaussian`, `:exponential`, `:hole`, `:power`,
        `:tpower`, `:expower`, `:cauchy`, `:rational`, `:bessel`, `:matern`,
        `:hyperspherical`, `:superspherical`, `:tplstable`
      - an RBF kernel function (any function is taken as a kernel), which is automatically wrapped via [[rbf->variogram]]; a model constructor such as [[gaussian]] is a function too, but it is not a kernel and is not accepted: use its keyword
  - `params` (map, optional): configuration map with the following keys (a `nil` value is the same as a missing key, except for `:weights` where `nil` means the uniform weights):
      - `:estimation` (keyword): loss function used during fitting; defaults to `:sq`:
          - `:sq` - squared difference (weighted least squares, default)
          - `:abs` - absolute difference
          - `:soft` - soft L1 (pseudo-Huber)
          - `:huber` - Huber loss
          - `:cauchy` - Cauchy/Lorentzian loss
          - `:atan` - arctangent loss
      - `:weights`: weighting scheme applied to the per-bin loss terms; defaults to `:n`:
          - `nil` - uniform weights (all `1.0`)
          - `:n` - number of pairs in each bin (default)
          - `:nh` - N divided by lag `h`
          - `:nh+` - N divided by scaled lag `h+` where `h+ = (h+1)/(hmax+1)`
          - `:nhh` - N divided by squared lag `h^2`
          - `:nhh+` - N divided by squared scaled lag `h+^2`
          - `:nexp` - N divided by `exp(1/h+)`
          - `:nsqrt` - N divided by `sqrt(h+)`
          - `:ngg` - Cressie weights: N divided by squared model gamma; recomputed from the
            fitted model after optimisation
          - `:ngg2` - N times empirical gamma divided by cubed model gamma; recomputed
            from the fitted model after optimisation. The loss is degenerate for the models of the power family (`:power`, `:tpower`, `:tplstable`): the weights grow without bound when the model goes to zero, and the optimum collapses (range at the lower bound, beta far from the data), do not use it for them
          - a sequence of explicit per-bin weights (finite numbers, one for every bin)
          - a function `(fn ^double [^double h])` returning a scale per bin; the final
            weight is N divided by the absolute value of the function result
      - `:defaults` (map): fixed parameter values excluded from optimisation; keys are a
        subset of `:nugget`, `:psill`, `:range`, `:beta`; fitting is performed only for
        the remaining parameters. When all parameters are provided, no optimisation runs.
      - `:parameter` (double): the order or shape parameter passed to model constructors
        that support it; used as the `order` argument for `:bessel`, `:matern`,
        `:hyperspherical`, `:superspherical`, and as `H` for `:tplstable`; when it is not given the
        default of the constructor is used: `1.0` for `:bessel` (the same as [[bessel]]), `1.0` for `:matern`,
        `2.0` for `:hyperspherical`, `1.0` for `:superspherical` and `0.5` for `:tplstable`
      - `:bounds` (map): overrides for the default parameter search bounds; keys are a
        subset of `:nugget`, `:psill`, `:range`, `:beta`, each mapping to a `[lo hi]`
        vector; default bounds are:
          - `:nugget` - `[0.0, 5 * min-gamma]`, where `min-gamma` is the smallest non-zero gamma
          - `:psill` - `[0.0, 5 * max-gamma]`
          - `:range` - `[EPSILON, max-h]`, where `max-h` is the largest lag of the empirical semivariogram
          - `:beta` - the valid interval of the model: `[EPSILON, 2.0]` for `:power`, `:expower` and `:tplstable`,
            `[EPSILON, 30.0]` for `:tpower`, `:rational` and `:cauchy`; a model checks its parameters when it is built, so
            a `:beta` bound outside the valid interval makes the evaluation throw
        The default bounds bind often: the range is capped by the largest lag and the nugget by the smallest gamma,
        so when the true range is beyond the data the fit stops at the bound (the fitted range equals `max-h`). Wider
        bounds (for example `{:range [1e-10 (* 10 max-h)]}`) improve the fit of most models, in the research data the loss dropped
        by up to 195% for single bins and the median change was negligible, but they degrade `:linear`, whose loss is
        flat in the range direction (about 2.5 times the optimal loss), and extremely wide bounds degrade every model because
        the same number of scanned points covers a much larger box (increase `:N` then)
      - `:rng`: random number generator (see [[fastmath.random/rng]]) drawing the scanned starting points of the
        optimisation, passed to `fastmath.optimization/scan-and-minimize`. Seed it to reproduce the result: two fits with equal
        seeded generators in the same state give bitwise equal parameters. The generator is consumed by the fit. When it is not given a new
        random generator is used for every call and the fitted parameters vary from run to run (below `1e-6` for most models, more
        for `:hole`, `:bessel` and `:linear`, which have flat or multimodal losses)
      - `:N` (long): number of scanned starting points of the optimisation, passed to
        `fastmath.optimization/scan-and-minimize` (its default is `100`); increase it for wide bounds or many parameters

  Returns a map containing:
  - `:semivariogram` (function): the fitted semivariogram `(fn ^double [^double x])`
  - `:semivariogram-creator` (function): the model constructor used to build the semivariogram
  - `:params` (map): all model parameters, both fitted and fixed via `:defaults`
  - `:weights` (sequence): final per-bin weights used to compute the loss
  - `:loss` (double): the final weighted loss value after fitting

  Throws `ex-info` (with `:type :invalid-input` and the key `:property`) for an empty empirical semivariogram, a non-finite `:h`, `:gamma` or `:n`
  (a `NaN` in the observations reaches the fit this way), an empirical semivariogram with only zero gamma, an unknown model keyword, `:estimation` or
  `:weights` keyword (the message lists the accepted values), explicit weights of a wrong length or unknown keys in `:defaults` and `:bounds`. A model which
  rejects its parameters throws `ex-info` as well (see [[exponential]]). Throws `ex-info` with `:type :fit-failed` when all the optimisation runs fail.

  See also [[fit]], [[empirical]], [[rbf->variogram]], [[variogram->covariance]]."
  ([empirical-semivariogram semivariogram-model] (fit-params empirical-semivariogram semivariogram-model nil))
  ([empirical-semivariogram semivariogram-model params]
   (let [{:keys [estimation weights defaults parameter bounds rng N]
          :or {estimation :sq weights :n}}
         (into {} (remove (fn [[k v]] (and (nil? v) (not= k :weights)))) params)
         _ (check-empirical! empirical-semivariogram)
         _ (check-parameter-keys! :defaults defaults)
         _ (check-parameter-keys! :bounds bounds)
         est-fn (or (estimation-fns estimation)
                    (invalid-input (str "Unknown estimation " (pr-str estimation)) :estimation estimation
                                   (vec (keys estimation-fns))))
         semivariogram-fn (cond
                            (fn? semivariogram-model) (rbf->variogram semivariogram-model)
                            (keyword? semivariogram-model)
                            (if-let [->sv (->semivariograms semivariogram-model)]
                              (if parameter (->sv parameter) (->sv))
                              (or (semivariograms semivariogram-model)
                                  (invalid-input (str "Unknown semivariogram model " semivariogram-model)
                                                 :model semivariogram-model (vec (sort (keys semivariograms))))))
                            :else (invalid-input (str "Unknown semivariogram model " (pr-str semivariogram-model))
                                                 :model semivariogram-model
                                                 (str "a keyword from " (vec (sort (keys semivariograms))) " or an RBF kernel function")))

         target-args (infer-target-args defaults (if (power-semivariograms semivariogram-model)
                                                   [:nugget :psill :range :beta]
                                                   [:nugget :psill :range]))

         hs (map :h empirical-semivariogram)
         gammas (map :gamma empirical-semivariogram)

         ns (map :n empirical-semivariogram)
         ng (map m/* ns gammas)

         max-h (double (last hs))
         nonzero-gammas (remove m/zero? gammas)
         _ (when-not (seq nonzero-gammas)
             (invalid-input "All the gamma of the empirical semivariogram are zero, there is nothing to fit."
                            :gamma gammas "at least one non-zero gamma"))
         [^double min-gamma ^double max-gamma] (stats/extent nonzero-gammas false)

         ;; hs+ -> (h+1)/(hmax+1)
         max-h+ (m/inc max-h)
         hs+ (map (fn [^double h] (m// (m/inc h) max-h+)) hs)

         weights-seq (weights-sequence weights ns hs hs+)

         params (if-not (seq target-args)
                  defaults
                  ;; fitting is needed
                  (let [bounds-map (merge {:nugget [0.0 (m/* 5.0 min-gamma)]
                                           :psill [0.0 (m/* 5.0 max-gamma)]
                                           :range [m/EPSILON max-h]
                                           :beta (get beta-bounds semivariogram-model [-30.0 30.0])}
                                          bounds)

                        bounds (map bounds-map target-args)
                        ;; a numerical gradient steps slightly out of the bounds (beta 9.9998e-11 for the lower bound 1e-10), still valid for a model
                        ->params (fn [r] (merge (zipmap target-args r) defaults))

                        target (case weights
                                 ;; cressie
                                 :ngg (fn [r] (let [sv (semivariogram-fn (->params r))
                                                   svh (map sv hs)
                                                   weights-seq (weights-n-v2 ns svh)]
                                               (-> (map (fn [^double gamma- ^double gamma]
                                                          (est-fn gamma- gamma)) svh gammas)
                                                   (v/dot weights-seq))))


                                 :ngg2 (fn [r] (let [sv (semivariogram-fn (->params r))
                                                    svh (map sv hs)
                                                    weights-seq (weights-n-v3 ng svh)]
                                                (-> (map (fn [^double gamma- ^double gamma]
                                                           (est-fn gamma- gamma)) svh gammas)
                                                    (v/dot weights-seq))))

                                 (fn [r] (let [sv (semivariogram-fn (->params r))]
                                          (-> (map (fn [^double h ^double gamma]
                                                     (est-fn (sv h) gamma)) hs gammas)
                                              (v/dot weights-seq)))))

                        m (optim/scan-and-minimize :lbfgsb target
                                                   (cond-> {:bounds bounds :jitter 0.05}
                                                     rng (assoc :rng rng)
                                                     N (assoc :N N)))]
                    (when-not m
                      (throw (ex-info "All the optimisation runs of the semivariogram fit failed."
                                      {:type :fit-failed :model semivariogram-model})))
                    (->> m
                         (first)
                         (zipmap target-args)
                         (merge defaults))))
         semivariogram (semivariogram-fn params)
         final-weights (case weights
                         :ngg (weights-n-v2 ns (map semivariogram hs))
                         :ngg2 (weights-n-v3 ng (map semivariogram hs))
                         weights-seq)]
     {:semivariogram semivariogram
      :semivariogram-creator semivariogram-fn
      :params params
      :weights final-weights
      :loss (v/dot final-weights (map (fn [^double h ^double gamma]
                                        (est-fn (semivariogram h) gamma)) hs gammas))})))

(defn fit
  "Fits a parametric semivariogram model to an empirical semivariogram and returns the fitted model function.

  A convenience wrapper around [[fit-params]] that returns only the fitted semivariogram
  function. Use [[fit-params]] directly when access to the fitted parameters, weights, or
  loss value is needed for diagnostics.

  Parameters:

  - `empirical-variogram` (sequence of maps): empirical semivariogram data as returned by
    [[empirical]]; each map must contain `:h`, `:gamma` and `:n` keys
  - `semivariogram-model` (keyword or function): the parametric model to fit; accepted values are:
      - a keyword naming a built-in model: `:linear`, `:spherical`, `:circular`,
        `:cubic`, `:pentaspherical`, `:gaussian`, `:exponential`, `:hole`, `:power`,
        `:tpower`, `:expower`, `:cauchy`, `:rational`, `:bessel`, `:matern`,
        `:hyperspherical`, `:superspherical`, `:tplstable`
      - an RBF kernel function, automatically wrapped via [[rbf->variogram]]; a function is always taken as a kernel, so a model constructor such as [[gaussian]] is not accepted (use `:gaussian`)
  - `parameters` (map, optional): configuration map; all keys are described in full in
    [[fit-params]]; the most commonly used keys are:
      - `:estimation` (keyword): loss function; `:sq` least squares (default), `:abs`,
        `:soft`, `:huber`, `:cauchy`, `:atan`
      - `:weights`: per-bin weighting scheme; defaults to `:n` (number of pairs); accepts
        `:n`, `:nh`, `:nh+`, `:nhh`, `:nhh+`, `:nexp`, `:nsqrt`, `:ngg`, `:ngg2`,
        a sequence of explicit weights, or a function `(fn ^double [^double h])`
      - `:defaults` (map): fixed parameter values excluded from optimisation; keys are a
        subset of `:nugget`, `:psill`, `:range`, `:beta`
      - `:parameter` (double): order for `:bessel` (default `1.0`), `:matern`, `:hyperspherical`,
        `:superspherical`; Hurst exponent `H` for `:tplstable`
      - `:bounds` (map): overrides for default parameter search bounds; keys are a subset
        of `:nugget`, `:psill`, `:range`, `:beta`, each mapping to a `[lo hi]` vector; the default
        bounds bind often (see the trade-off in [[fit-params]])
      - `:rng`, `:N`: random number generator and number of the scanned starting points; a seeded
        generator makes the fit reproducible

  Returns the fitted semivariogram function `(fn ^double [^double x])`.

  Throws `ex-info` for invalid input, see [[fit-params]].

  See also [[fit-params]], [[empirical]], [[rbf->variogram]]."
  ([empirical-variogram semivariogram-model] (fit empirical-variogram semivariogram-model nil))
  ([empirical-variogram semivariogram-model parameters]
   (:semivariogram (fit-params empirical-variogram semivariogram-model parameters))))
